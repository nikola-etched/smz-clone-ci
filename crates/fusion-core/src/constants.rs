//! Physical constants and unit conversions.

/// Elementary charge `[C]`, which is also the number of joules in one electronvolt.
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Joules per kiloelectronvolt.
pub const KEV_TO_JOULE: f64 = 1.0e3 * ELEMENTARY_CHARGE;

/// Joules per megaelectronvolt.
pub const MEV_TO_JOULE: f64 = 1.0e6 * ELEMENTARY_CHARGE;

/// Cubic metres per cubic centimetre.
pub const CM3_TO_M3: f64 = 1.0e-6;

/// Watts per megawatt.
pub const MW_TO_WATT: f64 = 1.0e6;

/// Total energy released by one D-T fusion reaction `[MeV]`.
pub const DT_ENERGY_MEV: f64 = 17.59;

/// Energy carried away by the alpha particle of a D-T reaction `[MeV]`.
pub const DT_ALPHA_ENERGY_MEV: f64 = 3.52;

/// Energy carried away by the neutron of a D-T reaction `[MeV]`.
pub const DT_NEUTRON_ENERGY_MEV: f64 = DT_ENERGY_MEV - DT_ALPHA_ENERGY_MEV;

/// Fraction of the D-T yield deposited in the plasma by the alpha particle.
pub const DT_ALPHA_FRACTION: f64 = DT_ALPHA_ENERGY_MEV / DT_ENERGY_MEV;

/// Bremsstrahlung emission coefficient `[W m^3 keV^-1/2]` for a hydrogenic plasma.
pub const BREMSSTRAHLUNG_COEFFICIENT: f64 = 5.34e-37;

/// Average ion mass of a 50:50 deuterium-tritium mixture `[amu]`.
pub const DT_MEAN_MASS_AMU: f64 = 2.5;
