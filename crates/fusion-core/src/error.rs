//! Error type shared by the plasma, confinement and burn models.

use std::fmt;

/// Errors produced when plasma or machine parameters fall outside the model's domain.
#[derive(Debug, Clone, PartialEq)]
pub enum FusionError {
    /// A quantity that must be strictly positive was not.
    NonPositive {
        /// Name of the offending quantity.
        name: &'static str,
        /// Value that was supplied.
        value: f64,
    },
    /// The temperature lies outside the validity range of the reactivity fit.
    TemperatureOutOfRange {
        /// Temperature that was requested `[keV]`.
        temperature_kev: f64,
        /// Lower bound of the fit `[keV]`.
        min_kev: f64,
        /// Upper bound of the fit `[keV]`.
        max_kev: f64,
    },
    /// A value that must lie in `[0, 1]` was outside that interval.
    FractionOutOfRange {
        /// Name of the offending quantity.
        name: &'static str,
        /// Value that was supplied.
        value: f64,
    },
    /// A quantity that must be zero or greater was negative.
    Negative {
        /// Name of the offending quantity.
        name: &'static str,
        /// Value that was supplied.
        value: f64,
    },
    /// The effective charge was below one, which no physical plasma can reach.
    InvalidEffectiveCharge(f64),
    /// The requested integration step does not divide into the requested duration.
    InvalidTimeStep {
        /// Duration of the requested run `[s]`.
        duration_s: f64,
        /// Integration step that was requested `[s]`.
        step_s: f64,
    },
}

impl fmt::Display for FusionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositive { name, value } => {
                write!(f, "{name} must be positive, got {value}")
            }
            Self::TemperatureOutOfRange {
                temperature_kev,
                min_kev,
                max_kev,
            } => write!(
                f,
                "temperature {temperature_kev} keV is outside the fit range \
                 [{min_kev}, {max_kev}] keV"
            ),
            Self::FractionOutOfRange { name, value } => {
                write!(f, "{name} must lie in [0, 1], got {value}")
            }
            Self::Negative { name, value } => {
                write!(f, "{name} must be zero or greater, got {value}")
            }
            Self::InvalidEffectiveCharge(value) => {
                write!(f, "effective charge must be at least 1, got {value}")
            }
            Self::InvalidTimeStep { duration_s, step_s } => {
                write!(
                    f,
                    "time step {step_s} s is not usable for a duration of {duration_s} s"
                )
            }
        }
    }
}

impl std::error::Error for FusionError {}

/// Result alias for fallible plasma calculations.
pub type Result<T> = std::result::Result<T, FusionError>;

pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FusionError::NonPositive { name, value })
    }
}

pub(crate) fn require_non_negative(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FusionError::Negative { name, value })
    }
}

pub(crate) fn require_fraction(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(FusionError::FractionOutOfRange { name, value })
    }
}
