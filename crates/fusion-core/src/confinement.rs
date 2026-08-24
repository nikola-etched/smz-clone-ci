//! Tokamak geometry and the ITER98(y,2) energy confinement scaling.

use serde::{Deserialize, Serialize};

use crate::error::{Result, require_positive};

/// Smallest loss power `[MW]` used by the confinement scaling, which diverges as power goes to zero.
const MIN_LOSS_POWER_MW: f64 = 0.1;

/// Leading coefficient of the ITER98(y,2) H-mode scaling.
const ITER98Y2_COEFFICIENT: f64 = 0.0562;

/// Shape and field parameters of an axisymmetric tokamak.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TokamakGeometry {
    /// Major radius `R` `[m]`.
    pub major_radius_m: f64,
    /// Minor radius `a` `[m]`.
    pub minor_radius_m: f64,
    /// Plasma elongation `kappa`.
    pub elongation: f64,
    /// Plasma current `I_p` `[MA]`.
    pub plasma_current_ma: f64,
    /// Toroidal magnetic field on axis `B_t` `[T]`.
    pub toroidal_field_t: f64,
}

impl TokamakGeometry {
    /// Builds a machine description, rejecting unphysical parameters.
    pub fn new(
        major_radius_m: f64,
        minor_radius_m: f64,
        elongation: f64,
        plasma_current_ma: f64,
        toroidal_field_t: f64,
    ) -> Result<Self> {
        require_positive("major radius", major_radius_m)?;
        require_positive("minor radius", minor_radius_m)?;
        require_positive("elongation", elongation)?;
        require_positive("plasma current", plasma_current_ma)?;
        require_positive("toroidal field", toroidal_field_t)?;
        Ok(Self {
            major_radius_m,
            minor_radius_m,
            elongation,
            plasma_current_ma,
            toroidal_field_t,
        })
    }

    /// The ITER baseline: R = 6.2 m, a = 2.0 m, kappa = 1.7, 15 MA, 5.3 T.
    pub fn iter() -> Self {
        Self {
            major_radius_m: 6.2,
            minor_radius_m: 2.0,
            elongation: 1.7,
            plasma_current_ma: 15.0,
            toroidal_field_t: 5.3,
        }
    }

    /// Aspect ratio `R / a`.
    pub fn aspect_ratio(&self) -> f64 {
        self.major_radius_m / self.minor_radius_m
    }

    /// Inverse aspect ratio `a / R`, the `epsilon` of the confinement scaling.
    pub fn inverse_aspect_ratio(&self) -> f64 {
        self.minor_radius_m / self.major_radius_m
    }

    /// Plasma volume `[m^3]` of an elliptical cross-section torus.
    pub fn plasma_volume_m3(&self) -> f64 {
        2.0 * std::f64::consts::PI.powi(2)
            * self.major_radius_m
            * self.minor_radius_m.powi(2)
            * self.elongation
    }

    /// Plasma surface area `[m^2]`, using the Ramanujan approximation for the elliptical perimeter.
    pub fn plasma_surface_m2(&self) -> f64 {
        let a = self.minor_radius_m;
        let b = self.minor_radius_m * self.elongation;
        let h = ((a - b) / (a + b)).powi(2);
        let perimeter =
            std::f64::consts::PI * (a + b) * (1.0 + 3.0 * h / (10.0 + (4.0 - 3.0 * h).sqrt()));
        2.0 * std::f64::consts::PI * self.major_radius_m * perimeter
    }

    /// Energy confinement time `[s]` from the ITER98(y,2) H-mode scaling law.
    ///
    /// `h_factor` multiplies the scaling, so `1.0` reproduces the reference H-mode.
    pub fn iter98y2_confinement_time_s(
        &self,
        loss_power_mw: f64,
        electron_density_m3: f64,
        ion_mass_amu: f64,
        h_factor: f64,
    ) -> f64 {
        let power_mw = loss_power_mw.max(MIN_LOSS_POWER_MW);
        let density_1e19 = electron_density_m3 / 1.0e19;

        ITER98Y2_COEFFICIENT
            * h_factor
            * self.plasma_current_ma.powf(0.93)
            * self.toroidal_field_t.powf(0.15)
            * power_mw.powf(-0.69)
            * density_1e19.powf(0.41)
            * ion_mass_amu.powf(0.19)
            * self.major_radius_m.powf(1.97)
            * self.inverse_aspect_ratio().powf(0.58)
            * self.elongation.powf(0.78)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn iter_volume_matches_the_published_value() {
        assert_relative_eq!(
            TokamakGeometry::iter().plasma_volume_m3(),
            830.0,
            max_relative = 0.02
        );
    }

    #[test]
    fn iter_confinement_time_matches_the_design_point() {
        let tau = TokamakGeometry::iter().iter98y2_confinement_time_s(87.0, 1.0e20, 2.5, 1.0);
        assert_relative_eq!(tau, 3.7, max_relative = 0.1);
    }

    #[test]
    fn confinement_degrades_with_heating_power() {
        let iter = TokamakGeometry::iter();
        let low = iter.iter98y2_confinement_time_s(50.0, 1.0e20, 2.5, 1.0);
        let high = iter.iter98y2_confinement_time_s(150.0, 1.0e20, 2.5, 1.0);
        assert!(high < low);
        assert_relative_eq!(
            low / high,
            (150.0f64 / 50.0).powf(0.69),
            max_relative = 1e-9
        );
    }

    #[test]
    fn confinement_scales_linearly_with_the_h_factor() {
        let iter = TokamakGeometry::iter();
        let base = iter.iter98y2_confinement_time_s(87.0, 1.0e20, 2.5, 1.0);
        let boosted = iter.iter98y2_confinement_time_s(87.0, 1.0e20, 2.5, 1.3);
        assert_relative_eq!(boosted / base, 1.3, max_relative = 1e-12);
    }

    #[test]
    fn vanishing_loss_power_stays_finite() {
        let tau = TokamakGeometry::iter().iter98y2_confinement_time_s(0.0, 1.0e20, 2.5, 1.0);
        assert!(tau.is_finite() && tau > 0.0);
    }

    #[test]
    fn bigger_machines_confine_better() {
        let small = TokamakGeometry::new(1.65, 0.5, 1.8, 2.0, 12.0).unwrap();
        let iter = TokamakGeometry::iter();
        assert!(
            iter.iter98y2_confinement_time_s(87.0, 1.0e20, 2.5, 1.0)
                > small.iter98y2_confinement_time_s(87.0, 1.0e20, 2.5, 1.0)
        );
    }

    #[test]
    fn aspect_ratios_are_reciprocal() {
        let iter = TokamakGeometry::iter();
        assert_relative_eq!(iter.aspect_ratio(), 3.1, max_relative = 1e-12);
        assert_relative_eq!(
            iter.aspect_ratio() * iter.inverse_aspect_ratio(),
            1.0,
            max_relative = 1e-12
        );
    }

    #[test]
    fn surface_area_is_positive_and_grows_with_elongation() {
        let round = TokamakGeometry::new(6.2, 2.0, 1.0, 15.0, 5.3).unwrap();
        let elongated = TokamakGeometry::iter();
        assert!(round.plasma_surface_m2() > 0.0);
        assert!(elongated.plasma_surface_m2() > round.plasma_surface_m2());
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(TokamakGeometry::new(0.0, 2.0, 1.7, 15.0, 5.3).is_err());
        assert!(TokamakGeometry::new(6.2, -2.0, 1.7, 15.0, 5.3).is_err());
        assert!(TokamakGeometry::new(6.2, 2.0, 1.7, 15.0, f64::NAN).is_err());
    }
}
