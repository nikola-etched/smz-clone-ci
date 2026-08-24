//! Zero-dimensional power balance model for magnetically confined fusion plasmas.
//!
//! The model combines the Bosch-Hale fusion reactivities, a uniform-plasma power balance and
//! the ITER98(y,2) energy confinement scaling into a burn simulation that can be integrated
//! in time.
//!
//! ```
//! use fusion_core::{BurnSimulation, PlasmaState};
//!
//! let plasma = PlasmaState::balanced_dt(1.0e20, 15.0)?;
//! assert!(plasma.fusion_power_density() > 0.0);
//!
//! let history = BurnSimulation::iter_baseline().run(30.0, 0.05)?;
//! assert!(history.final_sample().power_mw.fusion > 0.0);
//! # Ok::<(), fusion_core::FusionError>(())
//! ```

#![deny(missing_docs)]

pub mod burn;
pub mod confinement;
pub mod constants;
pub mod error;
pub mod plasma;
pub mod reactivity;

pub use burn::{BurnHistory, BurnSample, BurnSimulation};
pub use confinement::TokamakGeometry;
pub use error::{FusionError, Result};
pub use plasma::{PlasmaState, PowerBalance};
pub use reactivity::{Reaction, peak_temperature_kev, reactivity, reactivity_checked};
