//! Zero-dimensional power balance of a deuterium-tritium plasma.

use serde::{Deserialize, Serialize};

use crate::constants::{
    BREMSSTRAHLUNG_COEFFICIENT, DT_ALPHA_FRACTION, DT_ENERGY_MEV, KEV_TO_JOULE, MEV_TO_JOULE,
};
use crate::error::{FusionError, Result, require_fraction, require_positive};
use crate::reactivity::{Reaction, reactivity};

/// A uniform plasma with equal electron and ion temperatures.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlasmaState {
    /// Total fuel ion density `[m^-3]`.
    pub ion_density_m3: f64,
    /// Ion and electron temperature `[keV]`.
    pub temperature_kev: f64,
    /// Fraction of the fuel ions that are deuterium.
    pub deuterium_fraction: f64,
    /// Effective charge, which sets the bremsstrahlung enhancement over a pure hydrogen plasma.
    pub z_eff: f64,
}

impl PlasmaState {
    /// Builds a plasma state, rejecting unphysical parameters.
    pub fn new(
        ion_density_m3: f64,
        temperature_kev: f64,
        deuterium_fraction: f64,
        z_eff: f64,
    ) -> Result<Self> {
        require_positive("ion density", ion_density_m3)?;
        require_positive("temperature", temperature_kev)?;
        require_fraction("deuterium fraction", deuterium_fraction)?;
        if !z_eff.is_finite() || z_eff < 1.0 {
            return Err(FusionError::InvalidEffectiveCharge(z_eff));
        }
        Ok(Self {
            ion_density_m3,
            temperature_kev,
            deuterium_fraction,
            z_eff,
        })
    }

    /// A 50:50 deuterium-tritium plasma with no impurities.
    pub fn balanced_dt(ion_density_m3: f64, temperature_kev: f64) -> Result<Self> {
        Self::new(ion_density_m3, temperature_kev, 0.5, 1.0)
    }

    /// Returns the same state at a different temperature `[keV]`.
    pub fn at_temperature(self, temperature_kev: f64) -> Result<Self> {
        Self::new(
            self.ion_density_m3,
            temperature_kev,
            self.deuterium_fraction,
            self.z_eff,
        )
    }

    /// Electron density `[m^-3]`, following quasi-neutrality of a hydrogenic fuel mixture.
    pub fn electron_density_m3(&self) -> f64 {
        self.ion_density_m3
    }

    /// Deuterium density `[m^-3]`.
    pub fn deuterium_density_m3(&self) -> f64 {
        self.ion_density_m3 * self.deuterium_fraction
    }

    /// Tritium density `[m^-3]`.
    pub fn tritium_density_m3(&self) -> f64 {
        self.ion_density_m3 * (1.0 - self.deuterium_fraction)
    }

    /// D-T reaction rate density `[reactions m^-3 s^-1]`.
    pub fn reaction_rate_density(&self) -> f64 {
        self.deuterium_density_m3()
            * self.tritium_density_m3()
            * reactivity(Reaction::DeuteriumTritium, self.temperature_kev)
    }

    /// Total fusion power density, neutrons included `[W/m^3]`.
    pub fn fusion_power_density(&self) -> f64 {
        self.reaction_rate_density() * DT_ENERGY_MEV * MEV_TO_JOULE
    }

    /// Power density deposited in the plasma by alpha particles `[W/m^3]`.
    pub fn alpha_heating_density(&self) -> f64 {
        self.fusion_power_density() * DT_ALPHA_FRACTION
    }

    /// Power density carried out of the plasma by neutrons `[W/m^3]`.
    pub fn neutron_power_density(&self) -> f64 {
        self.fusion_power_density() * (1.0 - DT_ALPHA_FRACTION)
    }

    /// Bremsstrahlung radiation loss `[W/m^3]`.
    pub fn bremsstrahlung_density(&self) -> f64 {
        BREMSSTRAHLUNG_COEFFICIENT
            * self.z_eff
            * self.electron_density_m3().powi(2)
            * self.temperature_kev.sqrt()
    }

    /// Thermal energy stored per unit volume `[J/m^3]`, summing the electron and ion pressure.
    pub fn stored_energy_density(&self) -> f64 {
        3.0 * self.electron_density_m3() * self.temperature_kev * KEV_TO_JOULE
    }

    /// Energy lost to transport per unit volume `[W/m^3]` for a given confinement time.
    pub fn transport_loss_density(&self, confinement_time_s: f64) -> f64 {
        if confinement_time_s > 0.0 {
            self.stored_energy_density() / confinement_time_s
        } else {
            f64::INFINITY
        }
    }

    /// Lawson triple product `n T tau_E` `[keV s m^-3]`.
    pub fn triple_product(&self, confinement_time_s: f64) -> f64 {
        self.electron_density_m3() * self.temperature_kev * confinement_time_s
    }

    /// Power balance for a given confinement time and auxiliary heating density `[W/m^3]`.
    pub fn power_balance(&self, confinement_time_s: f64, auxiliary_density: f64) -> PowerBalance {
        let fusion = self.fusion_power_density();
        let alpha_heating = fusion * DT_ALPHA_FRACTION;
        let neutron = fusion - alpha_heating;
        let bremsstrahlung = self.bremsstrahlung_density();
        let transport_loss = self.transport_loss_density(confinement_time_s);
        PowerBalance {
            fusion,
            alpha_heating,
            neutron,
            bremsstrahlung,
            transport_loss,
            auxiliary: auxiliary_density,
        }
    }

    /// Smallest confinement time at which alpha heating covers every loss channel `[s]`.
    ///
    /// Returns `None` when bremsstrahlung alone already exceeds the alpha heating, which no
    /// amount of confinement can fix.
    pub fn ignition_confinement_time_s(&self) -> Option<f64> {
        let surplus = self.alpha_heating_density() - self.bremsstrahlung_density();
        (surplus > 0.0).then(|| self.stored_energy_density() / surplus)
    }
}

/// Power flows of a plasma, either as densities `[W/m^3]` or as totals after [`PowerBalance::scale`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PowerBalance {
    /// Total fusion power.
    pub fusion: f64,
    /// Alpha power deposited in the plasma.
    pub alpha_heating: f64,
    /// Power carried out by neutrons.
    pub neutron: f64,
    /// Bremsstrahlung radiation loss.
    pub bremsstrahlung: f64,
    /// Conduction and convection loss.
    pub transport_loss: f64,
    /// Externally applied heating.
    pub auxiliary: f64,
}

impl PowerBalance {
    /// Multiplies every power flow by `factor`, converting densities to totals.
    pub fn scale(&self, factor: f64) -> Self {
        Self {
            fusion: self.fusion * factor,
            alpha_heating: self.alpha_heating * factor,
            neutron: self.neutron * factor,
            bremsstrahlung: self.bremsstrahlung * factor,
            transport_loss: self.transport_loss * factor,
            auxiliary: self.auxiliary * factor,
        }
    }

    /// Net heating available to raise the plasma temperature.
    pub fn net_heating(&self) -> f64 {
        self.alpha_heating + self.auxiliary - self.bremsstrahlung - self.transport_loss
    }

    /// Power that leaves the plasma by transport and radiation, driving the confinement scaling.
    pub fn loss_power(&self) -> f64 {
        self.bremsstrahlung + self.transport_loss
    }

    /// Fusion energy gain `Q`, or `None` when no auxiliary heating is applied.
    pub fn q_factor(&self) -> Option<f64> {
        (self.auxiliary > 0.0).then(|| self.fusion / self.auxiliary)
    }

    /// Ratio of alpha heating to total losses; a plasma ignites once this reaches one.
    pub fn ignition_margin(&self) -> f64 {
        let losses = self.loss_power();
        if losses > 0.0 {
            self.alpha_heating / losses
        } else {
            f64::INFINITY
        }
    }

    /// Whether alpha heating alone sustains the plasma against its losses.
    pub fn is_ignited(&self) -> bool {
        self.alpha_heating >= self.loss_power()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn reference_plasma() -> PlasmaState {
        PlasmaState::balanced_dt(1.0e20, 10.0).unwrap()
    }

    #[test]
    fn fusion_power_density_is_sub_megawatt_at_iter_conditions() {
        let power = reference_plasma().fusion_power_density();
        assert_relative_eq!(power, 7.7e5, max_relative = 0.1);
    }

    #[test]
    fn alpha_and_neutron_power_split_the_yield() {
        let plasma = reference_plasma();
        let total = plasma.fusion_power_density();
        assert_relative_eq!(
            plasma.alpha_heating_density() + plasma.neutron_power_density(),
            total,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            plasma.alpha_heating_density() / total,
            0.2,
            max_relative = 0.01
        );
    }

    #[test]
    fn bremsstrahlung_matches_the_standard_coefficient() {
        assert_relative_eq!(
            reference_plasma().bremsstrahlung_density(),
            1.69e4,
            max_relative = 0.02
        );
    }

    #[test]
    fn bremsstrahlung_grows_with_impurity_content() {
        let clean = reference_plasma();
        let dirty = PlasmaState::new(1.0e20, 10.0, 0.5, 2.0).unwrap();
        assert_relative_eq!(
            dirty.bremsstrahlung_density() / clean.bremsstrahlung_density(),
            2.0,
            max_relative = 1e-12
        );
    }

    #[test]
    fn stored_energy_follows_three_n_t() {
        let plasma = reference_plasma();
        assert_relative_eq!(
            plasma.stored_energy_density(),
            3.0 * 1.0e20 * 10.0 * KEV_TO_JOULE,
            max_relative = 1e-12
        );
    }

    #[test]
    fn unbalanced_fuel_mix_reduces_fusion_power() {
        let balanced = reference_plasma().fusion_power_density();
        let lean = PlasmaState::new(1.0e20, 10.0, 0.8, 1.0)
            .unwrap()
            .fusion_power_density();
        assert!(lean < balanced);
        assert_relative_eq!(lean / balanced, 0.8 * 0.2 / 0.25, max_relative = 1e-12);
    }

    #[test]
    fn pure_deuterium_produces_no_dt_power() {
        let plasma = PlasmaState::new(1.0e20, 10.0, 1.0, 1.0).unwrap();
        assert_eq!(plasma.fusion_power_density(), 0.0);
    }

    #[test]
    fn q_factor_is_none_without_auxiliary_heating() {
        let balance = reference_plasma().power_balance(3.0, 0.0);
        assert!(balance.q_factor().is_none());
    }

    #[test]
    fn q_factor_is_the_ratio_of_fusion_to_auxiliary_power() {
        let plasma = reference_plasma();
        let auxiliary = 1.0e4;
        let balance = plasma.power_balance(3.0, auxiliary);
        assert_relative_eq!(
            balance.q_factor().unwrap(),
            plasma.fusion_power_density() / auxiliary,
            max_relative = 1e-12
        );
    }

    #[test]
    fn scaling_preserves_the_gain_factor() {
        let balance = reference_plasma().power_balance(3.0, 1.0e4);
        let totals = balance.scale(830.0);
        assert_relative_eq!(
            totals.q_factor().unwrap(),
            balance.q_factor().unwrap(),
            max_relative = 1e-12
        );
        assert_relative_eq!(totals.fusion, balance.fusion * 830.0, max_relative = 1e-12);
    }

    #[test]
    fn longer_confinement_reduces_transport_loss() {
        let plasma = reference_plasma();
        assert!(plasma.transport_loss_density(6.0) < plasma.transport_loss_density(3.0));
    }

    #[test]
    fn good_confinement_reaches_ignition() {
        let plasma = reference_plasma();
        let required = plasma.ignition_confinement_time_s().unwrap();
        assert!(plasma.power_balance(required * 1.01, 0.0).is_ignited());
        assert!(!plasma.power_balance(required * 0.99, 0.0).is_ignited());
    }

    #[test]
    fn cold_plasma_cannot_ignite() {
        let plasma = PlasmaState::balanced_dt(1.0e20, 1.0).unwrap();
        assert!(plasma.ignition_confinement_time_s().is_none());
    }

    #[test]
    fn triple_product_scales_with_each_factor() {
        let plasma = reference_plasma();
        assert_relative_eq!(
            plasma.triple_product(3.0),
            1.0e20 * 10.0 * 3.0,
            max_relative = 1e-12
        );
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        assert!(matches!(
            PlasmaState::balanced_dt(-1.0, 10.0),
            Err(FusionError::NonPositive { .. })
        ));
        assert!(matches!(
            PlasmaState::balanced_dt(1.0e20, 0.0),
            Err(FusionError::NonPositive { .. })
        ));
        assert!(matches!(
            PlasmaState::new(1.0e20, 10.0, 1.5, 1.0),
            Err(FusionError::FractionOutOfRange { .. })
        ));
        assert!(matches!(
            PlasmaState::new(1.0e20, 10.0, 0.5, 0.5),
            Err(FusionError::InvalidEffectiveCharge(_))
        ));
    }
}
