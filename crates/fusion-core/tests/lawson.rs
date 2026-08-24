//! End-to-end checks against textbook results for the Lawson criterion and burning plasmas.

use approx::assert_relative_eq;
use fusion_core::{BurnSimulation, PlasmaState, Reaction, TokamakGeometry, reactivity};

/// The ignition triple product of a D-T plasma bottoms out in the 10-20 keV window and is of
/// order 1e21 keV s m^-3 there.
#[test]
fn ignition_triple_product_minimum_sits_in_the_textbook_window() {
    let mut best_temperature = 0.0;
    let mut best_triple_product = f64::INFINITY;

    let mut temperature = 2.0;
    while temperature <= 60.0 {
        let plasma = PlasmaState::balanced_dt(1.0e20, temperature).unwrap();
        if let Some(tau) = plasma.ignition_confinement_time_s() {
            let triple_product = plasma.triple_product(tau);
            if triple_product < best_triple_product {
                best_triple_product = triple_product;
                best_temperature = temperature;
            }
        }
        temperature += 0.5;
    }

    assert!(
        (10.0..=30.0).contains(&best_temperature),
        "optimum was {best_temperature} keV"
    );
    assert!(
        (1.0e21..=1.0e22).contains(&best_triple_product),
        "triple product was {best_triple_product} keV s m^-3"
    );
}

/// Ignition requires a longer confinement time as the plasma is diluted with impurities.
#[test]
fn impurities_raise_the_ignition_requirement() {
    let clean = PlasmaState::new(1.0e20, 15.0, 0.5, 1.0).unwrap();
    let dirty = PlasmaState::new(1.0e20, 15.0, 0.5, 2.5).unwrap();

    let clean_tau = clean.ignition_confinement_time_s().unwrap();
    let dirty_tau = dirty.ignition_confinement_time_s().unwrap();
    assert!(
        dirty_tau > clean_tau,
        "{dirty_tau} should exceed {clean_tau}"
    );
}

/// The ITER baseline scenario lands on a burning plasma with a gain of a few and a fusion
/// power in the hundreds of megawatts.
#[test]
fn iter_baseline_produces_a_burning_plasma() {
    let history = BurnSimulation::iter_baseline().run(90.0, 0.05).unwrap();
    let final_sample = history.final_sample();

    assert!(
        (50.0..=2000.0).contains(&final_sample.power_mw.fusion),
        "fusion power was {} MW",
        final_sample.power_mw.fusion
    );
    assert!(
        final_sample.q_factor().unwrap() > 1.0,
        "gain was {:?}",
        final_sample.q_factor()
    );
    assert!(
        (1.0..=10.0).contains(&final_sample.confinement_time_s),
        "confinement time was {} s",
        final_sample.confinement_time_s
    );
}

/// Alpha heating is what makes the burn self-amplifying: removing it collapses the gain.
#[test]
fn alpha_heating_dominates_the_final_power_balance() {
    let history = BurnSimulation::iter_baseline().run(90.0, 0.05).unwrap();
    let final_sample = history.final_sample();
    assert!(
        final_sample.power_mw.alpha_heating > final_sample.power_mw.auxiliary,
        "alpha {} MW did not exceed auxiliary {} MW",
        final_sample.power_mw.alpha_heating,
        final_sample.power_mw.auxiliary
    );
}

/// A machine at half the linear size of ITER cannot reach the same gain at fixed density.
#[test]
fn a_smaller_machine_reaches_a_lower_gain() {
    let small = BurnSimulation {
        geometry: TokamakGeometry::new(3.1, 1.0, 1.7, 7.5, 5.3).unwrap(),
        ..BurnSimulation::iter_baseline()
    };
    let small_gain = small
        .run(90.0, 0.05)
        .unwrap()
        .final_sample()
        .q_factor()
        .unwrap();
    let iter_gain = BurnSimulation::iter_baseline()
        .run(90.0, 0.05)
        .unwrap()
        .final_sample()
        .q_factor()
        .unwrap();

    assert!(
        small_gain < iter_gain,
        "{small_gain} should be below {iter_gain}"
    );
}

/// Power density follows the square of the density at fixed temperature.
#[test]
fn fusion_power_scales_with_the_square_of_density() {
    let low = PlasmaState::balanced_dt(1.0e20, 15.0).unwrap();
    let high = PlasmaState::balanced_dt(2.0e20, 15.0).unwrap();
    assert_relative_eq!(
        high.fusion_power_density() / low.fusion_power_density(),
        4.0,
        max_relative = 1e-12
    );
}

/// D-3He needs a far hotter plasma than D-T to reach a comparable reactivity.
#[test]
fn deuterium_helium3_peaks_at_a_higher_temperature_than_deuterium_tritium() {
    let dt_peak = fusion_core::peak_temperature_kev(Reaction::DeuteriumTritium);
    let dhe3_peak = fusion_core::peak_temperature_kev(Reaction::DeuteriumHelium3);
    assert!(
        dhe3_peak > dt_peak,
        "{dhe3_peak} keV should exceed {dt_peak} keV"
    );
    assert!(
        reactivity(Reaction::DeuteriumTritium, 20.0)
            > 10.0 * reactivity(Reaction::DeuteriumHelium3, 20.0)
    );
}
