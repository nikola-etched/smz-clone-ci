//! Maxwellian-averaged fusion reactivities from the Bosch-Hale parameterization.
//!
//! H.-S. Bosch and G. M. Hale, *Improved formulas for fusion cross-sections and thermal
//! reactivities*, Nuclear Fusion **32** (1992) 611.

use serde::{Deserialize, Serialize};

use crate::constants::CM3_TO_M3;
use crate::error::{FusionError, Result};

/// A thermonuclear reaction channel covered by the Bosch-Hale fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Reaction {
    /// `T(d,n)4He`
    DeuteriumTritium,
    /// `D(d,p)T`
    DeuteriumDeuteriumProton,
    /// `D(d,n)3He`
    DeuteriumDeuteriumNeutron,
    /// `3He(d,p)4He`
    DeuteriumHelium3,
}

/// Coefficients of the Bosch-Hale reactivity fit for one reaction channel.
struct Fit {
    /// Gamow constant `[keV^1/2]`.
    gamow_kev_sqrt: f64,
    /// Reduced mass energy `[keV]`.
    reduced_mass_kev: f64,
    /// Padé coefficients `C1..C7`.
    c: [f64; 7],
    /// Lower bound of the fitted temperature range `[keV]`.
    min_kev: f64,
    /// Upper bound of the fitted temperature range `[keV]`.
    max_kev: f64,
}

const DT_FIT: Fit = Fit {
    gamow_kev_sqrt: 34.3827,
    reduced_mass_kev: 1_124_656.0,
    c: [
        1.173_02e-9,
        1.513_61e-2,
        7.518_86e-2,
        4.606_43e-3,
        1.350_00e-2,
        -1.067_50e-4,
        1.366_00e-5,
    ],
    min_kev: 0.2,
    max_kev: 100.0,
};

const DDP_FIT: Fit = Fit {
    gamow_kev_sqrt: 31.3970,
    reduced_mass_kev: 937_814.0,
    c: [
        5.657_18e-12,
        3.412_67e-3,
        1.991_67e-3,
        0.0,
        1.050_60e-5,
        0.0,
        0.0,
    ],
    min_kev: 0.2,
    max_kev: 100.0,
};

const DDN_FIT: Fit = Fit {
    gamow_kev_sqrt: 31.3970,
    reduced_mass_kev: 937_814.0,
    c: [
        5.433_60e-12,
        5.857_78e-3,
        7.682_22e-3,
        0.0,
        -2.964_00e-6,
        0.0,
        0.0,
    ],
    min_kev: 0.2,
    max_kev: 100.0,
};

const DHE3_FIT: Fit = Fit {
    gamow_kev_sqrt: 68.7508,
    reduced_mass_kev: 1_124_572.0,
    c: [
        5.510_36e-10,
        6.419_18e-3,
        -2.028_96e-3,
        -1.910_80e-5,
        1.357_76e-4,
        0.0,
        0.0,
    ],
    min_kev: 0.5,
    max_kev: 190.0,
};

impl Reaction {
    fn fit(self) -> &'static Fit {
        match self {
            Self::DeuteriumTritium => &DT_FIT,
            Self::DeuteriumDeuteriumProton => &DDP_FIT,
            Self::DeuteriumDeuteriumNeutron => &DDN_FIT,
            Self::DeuteriumHelium3 => &DHE3_FIT,
        }
    }

    /// Temperature range over which the fit was published, in keV.
    pub fn validity_range_kev(self) -> (f64, f64) {
        let fit = self.fit();
        (fit.min_kev, fit.max_kev)
    }

    /// Short label used by the command line interface and by serialized output.
    pub fn label(self) -> &'static str {
        match self {
            Self::DeuteriumTritium => "D-T",
            Self::DeuteriumDeuteriumProton => "D-D (p)",
            Self::DeuteriumDeuteriumNeutron => "D-D (n)",
            Self::DeuteriumHelium3 => "D-3He",
        }
    }
}

/// Maxwellian-averaged reactivity `<sigma v>` `[m^3/s]` at an ion temperature in keV.
///
/// The fit is evaluated for any positive temperature; outside
/// [`Reaction::validity_range_kev`] the result is an extrapolation. Use
/// [`reactivity_checked`] to reject those temperatures instead.
pub fn reactivity(reaction: Reaction, temperature_kev: f64) -> f64 {
    if !temperature_kev.is_finite() || temperature_kev <= 0.0 {
        return 0.0;
    }

    let fit = reaction.fit();
    let t = temperature_kev;
    let [c1, c2, c3, c4, c5, c6, c7] = fit.c;

    let numerator = t * (c2 + t * (c4 + t * c6));
    let denominator = 1.0 + t * (c3 + t * (c5 + t * c7));
    let theta = t / (1.0 - numerator / denominator);

    let xi = (fit.gamow_kev_sqrt.powi(2) / (4.0 * theta)).cbrt();

    let sigma_v_cm3 =
        c1 * theta * (xi / (fit.reduced_mass_kev * t.powi(3))).sqrt() * (-3.0 * xi).exp();

    sigma_v_cm3 * CM3_TO_M3
}

/// Reactivity `<sigma v>` `[m^3/s]`, rejecting temperatures outside the fitted range.
pub fn reactivity_checked(reaction: Reaction, temperature_kev: f64) -> Result<f64> {
    let (min_kev, max_kev) = reaction.validity_range_kev();
    if !(min_kev..=max_kev).contains(&temperature_kev) {
        return Err(FusionError::TemperatureOutOfRange {
            temperature_kev,
            min_kev,
            max_kev,
        });
    }
    Ok(reactivity(reaction, temperature_kev))
}

/// Temperature of peak reactivity `[keV]`, located by golden-section search over the fit range.
pub fn peak_temperature_kev(reaction: Reaction) -> f64 {
    const GOLDEN_RATIO_INVERSE: f64 = 0.618_033_988_749_895;
    const ITERATIONS: usize = 200;

    let (mut low, mut high) = reaction.validity_range_kev();
    for _ in 0..ITERATIONS {
        let span = high - low;
        let probe_low = high - GOLDEN_RATIO_INVERSE * span;
        let probe_high = low + GOLDEN_RATIO_INVERSE * span;
        if reactivity(reaction, probe_low) < reactivity(reaction, probe_high) {
            low = probe_low;
        } else {
            high = probe_high;
        }
    }
    0.5 * (low + high)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn deuterium_tritium_matches_published_value_at_10_kev() {
        assert_relative_eq!(
            reactivity(Reaction::DeuteriumTritium, 10.0),
            1.1e-22,
            max_relative = 0.05
        );
    }

    #[test]
    fn deuterium_tritium_peaks_near_64_kev() {
        let peak = peak_temperature_kev(Reaction::DeuteriumTritium);
        assert!((60.0..=70.0).contains(&peak), "peak was {peak} keV");
    }

    #[test]
    fn deuterium_tritium_dominates_deuterium_deuterium() {
        let dt = reactivity(Reaction::DeuteriumTritium, 10.0);
        let ddp = reactivity(Reaction::DeuteriumDeuteriumProton, 10.0);
        let ddn = reactivity(Reaction::DeuteriumDeuteriumNeutron, 10.0);
        let ratio = dt / (ddp + ddn);
        assert!((50.0..=500.0).contains(&ratio), "ratio was {ratio}");
    }

    #[test]
    fn deuterium_deuterium_branches_are_comparable() {
        let ddp = reactivity(Reaction::DeuteriumDeuteriumProton, 20.0);
        let ddn = reactivity(Reaction::DeuteriumDeuteriumNeutron, 20.0);
        assert_relative_eq!(ddp, ddn, max_relative = 0.2);
    }

    #[test]
    fn reactivity_increases_monotonically_below_the_peak() {
        let mut previous = 0.0;
        let mut temperature = 1.0;
        while temperature < 50.0 {
            let current = reactivity(Reaction::DeuteriumTritium, temperature);
            assert!(current > previous, "not increasing at {temperature} keV");
            previous = current;
            temperature += 1.0;
        }
    }

    #[test]
    fn non_positive_temperature_gives_zero() {
        assert_eq!(reactivity(Reaction::DeuteriumTritium, 0.0), 0.0);
        assert_eq!(reactivity(Reaction::DeuteriumTritium, -5.0), 0.0);
        assert_eq!(reactivity(Reaction::DeuteriumTritium, f64::NAN), 0.0);
    }

    #[test]
    fn checked_reactivity_rejects_temperatures_outside_the_fit() {
        let error = reactivity_checked(Reaction::DeuteriumTritium, 500.0).unwrap_err();
        assert!(matches!(error, FusionError::TemperatureOutOfRange { .. }));
        assert!(reactivity_checked(Reaction::DeuteriumTritium, 10.0).is_ok());
    }
}
