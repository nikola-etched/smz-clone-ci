//! Time evolution of a burning plasma under its own power balance.

use serde::{Deserialize, Serialize};

use crate::confinement::TokamakGeometry;
use crate::constants::{DT_MEAN_MASS_AMU, KEV_TO_JOULE, MW_TO_WATT};
use crate::error::{FusionError, Result, require_non_negative, require_positive};
use crate::plasma::{PlasmaState, PowerBalance};

/// Floor applied to the temperature while integrating, in keV.
const MIN_TEMPERATURE_KEV: f64 = 0.05;

/// A zero-dimensional burn scenario at fixed density.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BurnSimulation {
    /// Machine the plasma sits in.
    pub geometry: TokamakGeometry,
    /// Fuel ion density held constant through the burn `[m^-3]`.
    pub ion_density_m3: f64,
    /// Fraction of the fuel ions that are deuterium.
    pub deuterium_fraction: f64,
    /// Effective charge of the plasma.
    pub z_eff: f64,
    /// Mean fuel ion mass `[amu]`, which enters the confinement scaling.
    pub ion_mass_amu: f64,
    /// Multiplier on the ITER98(y,2) confinement time.
    pub h_factor: f64,
    /// Auxiliary heating power delivered to the plasma `[MW]`.
    pub auxiliary_power_mw: f64,
    /// Temperature at `t = 0` `[keV]`.
    pub initial_temperature_kev: f64,
}

impl BurnSimulation {
    /// Builds an ITER-like scenario that can then be adjusted field by field.
    pub fn iter_baseline() -> Self {
        Self {
            geometry: TokamakGeometry::iter(),
            ion_density_m3: 1.0e20,
            deuterium_fraction: 0.5,
            z_eff: 1.5,
            ion_mass_amu: DT_MEAN_MASS_AMU,
            h_factor: 1.0,
            auxiliary_power_mw: 50.0,
            initial_temperature_kev: 5.0,
        }
    }

    /// Plasma state at a given temperature `[keV]`.
    pub fn plasma_at(&self, temperature_kev: f64) -> Result<PlasmaState> {
        PlasmaState::new(
            self.ion_density_m3,
            temperature_kev.max(MIN_TEMPERATURE_KEV),
            self.deuterium_fraction,
            self.z_eff,
        )
    }

    /// Auxiliary heating expressed as a power density `[W/m^3]`.
    fn auxiliary_density(&self) -> f64 {
        self.auxiliary_power_mw * MW_TO_WATT / self.geometry.plasma_volume_m3()
    }

    /// Energy confinement time `[s]` for the plasma's current heating power.
    fn confinement_time_s(&self, plasma: &PlasmaState) -> f64 {
        let volume = self.geometry.plasma_volume_m3();
        let heating_density = plasma.alpha_heating_density() + self.auxiliary_density()
            - plasma.bremsstrahlung_density();
        let heating_mw = heating_density * volume / MW_TO_WATT;
        self.geometry.iter98y2_confinement_time_s(
            heating_mw,
            plasma.electron_density_m3(),
            self.ion_mass_amu,
            self.h_factor,
        )
    }

    /// Instantaneous diagnostics of the scenario at a given temperature.
    pub fn sample_at(&self, time_s: f64, temperature_kev: f64) -> Result<BurnSample> {
        let plasma = self.plasma_at(temperature_kev)?;
        let confinement_time_s = self.confinement_time_s(&plasma);
        let balance = plasma.power_balance(confinement_time_s, self.auxiliary_density());
        Ok(BurnSample {
            time_s,
            temperature_kev: plasma.temperature_kev,
            confinement_time_s,
            triple_product: plasma.triple_product(confinement_time_s),
            power_mw: balance.scale(self.geometry.plasma_volume_m3() / MW_TO_WATT),
        })
    }

    /// Rate of change of the temperature `[keV/s]` from the net heating.
    fn temperature_derivative(&self, temperature_kev: f64) -> Result<f64> {
        let plasma = self.plasma_at(temperature_kev)?;
        let confinement_time_s = self.confinement_time_s(&plasma);
        let balance = plasma.power_balance(confinement_time_s, self.auxiliary_density());
        Ok(balance.net_heating() / (3.0 * plasma.electron_density_m3() * KEV_TO_JOULE))
    }

    /// Integrates the power balance with a fixed-step fourth-order Runge-Kutta scheme.
    pub fn run(&self, duration_s: f64, step_s: f64) -> Result<BurnHistory> {
        require_positive("ion density", self.ion_density_m3)?;
        require_non_negative("auxiliary power", self.auxiliary_power_mw)?;
        require_positive("h factor", self.h_factor)?;
        require_positive("ion mass", self.ion_mass_amu)?;
        if !duration_s.is_finite()
            || !step_s.is_finite()
            || duration_s <= 0.0
            || step_s <= 0.0
            || step_s > duration_s
        {
            return Err(FusionError::InvalidTimeStep { duration_s, step_s });
        }

        let steps = (duration_s / step_s).ceil() as usize;
        let mut temperature = self.initial_temperature_kev;
        let mut samples = Vec::with_capacity(steps + 1);
        samples.push(self.sample_at(0.0, temperature)?);

        for step in 0..steps {
            let k1 = self.temperature_derivative(temperature)?;
            let k2 = self.temperature_derivative(temperature + 0.5 * step_s * k1)?;
            let k3 = self.temperature_derivative(temperature + 0.5 * step_s * k2)?;
            let k4 = self.temperature_derivative(temperature + step_s * k3)?;
            temperature += step_s / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            temperature = temperature.max(MIN_TEMPERATURE_KEV);
            samples.push(self.sample_at((step + 1) as f64 * step_s, temperature)?);
        }

        Ok(BurnHistory { samples })
    }
}

/// Diagnostics of the plasma at one instant of a burn.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BurnSample {
    /// Time since the start of the burn `[s]`.
    pub time_s: f64,
    /// Plasma temperature `[keV]`.
    pub temperature_kev: f64,
    /// Energy confinement time `[s]`.
    pub confinement_time_s: f64,
    /// Lawson triple product `[keV s m^-3]`.
    pub triple_product: f64,
    /// Machine-total power flows `[MW]`.
    pub power_mw: PowerBalance,
}

impl BurnSample {
    /// Fusion energy gain `Q` at this instant.
    pub fn q_factor(&self) -> Option<f64> {
        self.power_mw.q_factor()
    }

    /// Whether alpha heating alone sustains the plasma at this instant.
    pub fn is_ignited(&self) -> bool {
        self.power_mw.is_ignited()
    }
}

/// The full time trace produced by [`BurnSimulation::run`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BurnHistory {
    /// Samples in increasing time order, starting at `t = 0`.
    pub samples: Vec<BurnSample>,
}

impl BurnHistory {
    /// State the plasma reached at the end of the run.
    pub fn final_sample(&self) -> &BurnSample {
        self.samples.last().expect("a run always records t = 0")
    }

    /// Highest fusion power reached during the run `[MW]`.
    pub fn peak_fusion_power_mw(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.power_mw.fusion)
            .fold(0.0, f64::max)
    }

    /// Highest temperature reached during the run `[keV]`.
    pub fn peak_temperature_kev(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.temperature_kev)
            .fold(0.0, f64::max)
    }

    /// Total fusion energy released over the run `[GJ]`, by trapezoidal integration.
    pub fn fusion_energy_gj(&self) -> f64 {
        self.samples
            .windows(2)
            .map(|pair| {
                let dt = pair[1].time_s - pair[0].time_s;
                0.5 * (pair[0].power_mw.fusion + pair[1].power_mw.fusion) * dt
            })
            .sum::<f64>()
            / 1.0e3
    }

    /// Time at which alpha heating first covers every loss channel `[s]`.
    pub fn ignition_time_s(&self) -> Option<f64> {
        self.samples
            .iter()
            .find(|sample| sample.is_ignited())
            .map(|sample| sample.time_s)
    }

    /// At most `max_samples` evenly spaced samples, always keeping the first and the last.
    pub fn downsample(&self, max_samples: usize) -> Vec<&BurnSample> {
        if max_samples == 0 {
            return Vec::new();
        }
        if self.samples.len() <= max_samples {
            return self.samples.iter().collect();
        }
        let stride = (self.samples.len() - 1) as f64 / (max_samples - 1) as f64;
        (0..max_samples)
            .map(|index| {
                &self.samples
                    [((index as f64 * stride).round() as usize).min(self.samples.len() - 1)]
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn baseline_scenario_heats_up_and_settles() {
        let history = BurnSimulation::iter_baseline().run(60.0, 0.05).unwrap();
        let final_sample = history.final_sample();
        assert!(
            final_sample.temperature_kev > BurnSimulation::iter_baseline().initial_temperature_kev,
            "temperature did not rise: {}",
            final_sample.temperature_kev
        );
        assert!(
            (5.0..60.0).contains(&final_sample.temperature_kev),
            "temperature ran away to {} keV",
            final_sample.temperature_kev
        );
    }

    #[test]
    fn baseline_scenario_reaches_a_useful_gain() {
        let history = BurnSimulation::iter_baseline().run(60.0, 0.05).unwrap();
        let q = history.final_sample().q_factor().unwrap();
        assert!(q > 1.0, "gain was only {q}");
    }

    #[test]
    fn every_sample_is_physically_finite() {
        let history = BurnSimulation::iter_baseline().run(20.0, 0.05).unwrap();
        for sample in &history.samples {
            assert!(sample.temperature_kev >= MIN_TEMPERATURE_KEV);
            assert!(sample.confinement_time_s.is_finite() && sample.confinement_time_s > 0.0);
            assert!(sample.power_mw.fusion.is_finite() && sample.power_mw.fusion >= 0.0);
            assert!(sample.triple_product.is_finite());
        }
    }

    #[test]
    fn samples_are_ordered_in_time() {
        let history = BurnSimulation::iter_baseline().run(10.0, 0.1).unwrap();
        assert_eq!(history.samples.len(), 101);
        assert_relative_eq!(history.samples[0].time_s, 0.0, epsilon = 1e-12);
        assert_relative_eq!(history.final_sample().time_s, 10.0, max_relative = 1e-9);
        for pair in history.samples.windows(2) {
            assert!(pair[1].time_s > pair[0].time_s);
        }
    }

    #[test]
    fn integration_is_insensitive_to_the_step_size() {
        let scenario = BurnSimulation::iter_baseline();
        let coarse = scenario.run(30.0, 0.2).unwrap();
        let fine = scenario.run(30.0, 0.02).unwrap();
        assert_relative_eq!(
            coarse.final_sample().temperature_kev,
            fine.final_sample().temperature_kev,
            max_relative = 1e-3
        );
    }

    #[test]
    fn more_auxiliary_heating_gives_a_hotter_plasma_and_a_lower_gain() {
        let weak = BurnSimulation {
            auxiliary_power_mw: 30.0,
            ..BurnSimulation::iter_baseline()
        }
        .run(60.0, 0.05)
        .unwrap();
        let strong = BurnSimulation {
            auxiliary_power_mw: 90.0,
            ..BurnSimulation::iter_baseline()
        }
        .run(60.0, 0.05)
        .unwrap();

        assert!(strong.final_sample().temperature_kev > weak.final_sample().temperature_kev);
        assert!(
            strong.final_sample().q_factor().unwrap() < weak.final_sample().q_factor().unwrap()
        );
    }

    #[test]
    fn poor_confinement_quenches_the_burn() {
        let history = BurnSimulation {
            h_factor: 0.3,
            auxiliary_power_mw: 5.0,
            ..BurnSimulation::iter_baseline()
        }
        .run(60.0, 0.05)
        .unwrap();
        assert!(history.final_sample().power_mw.fusion < history.peak_fusion_power_mw());
        assert!(history.ignition_time_s().is_none());
    }

    #[test]
    fn fusion_energy_accumulates_over_the_run() {
        let history = BurnSimulation::iter_baseline().run(60.0, 0.05).unwrap();
        let energy = history.fusion_energy_gj();
        assert!(energy > 0.0);
        assert!(energy <= history.peak_fusion_power_mw() * 60.0 / 1.0e3);
    }

    #[test]
    fn downsampling_keeps_the_endpoints() {
        let history = BurnSimulation::iter_baseline().run(10.0, 0.01).unwrap();
        let rows = history.downsample(11);
        assert_eq!(rows.len(), 11);
        assert_eq!(rows[0].time_s, history.samples[0].time_s);
        assert_eq!(rows[10].time_s, history.final_sample().time_s);
        assert!(history.downsample(0).is_empty());
        assert_eq!(history.downsample(usize::MAX).len(), history.samples.len());
    }

    #[test]
    fn invalid_time_steps_are_rejected() {
        let scenario = BurnSimulation::iter_baseline();
        assert!(matches!(
            scenario.run(10.0, 0.0),
            Err(FusionError::InvalidTimeStep { .. })
        ));
        assert!(matches!(
            scenario.run(0.0, 0.1),
            Err(FusionError::InvalidTimeStep { .. })
        ));
        assert!(matches!(
            scenario.run(1.0, 5.0),
            Err(FusionError::InvalidTimeStep { .. })
        ));
    }
}
