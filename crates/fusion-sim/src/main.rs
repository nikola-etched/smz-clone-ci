//! Command line front end for the zero-dimensional fusion plasma model.

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use fusion_core::{
    BurnHistory, BurnSample, BurnSimulation, PlasmaState, Reaction, TokamakGeometry,
    peak_temperature_kev, reactivity,
};
use serde::Serialize;

/// Simulate the power balance of a magnetically confined deuterium-tritium plasma.
#[derive(Debug, Parser)]
#[command(name = "fusion-sim", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Integrate a burn scenario in time.
    Burn(BurnArgs),
    /// Tabulate the Maxwellian-averaged reactivity against temperature.
    Reactivity(ReactivityArgs),
    /// Scan the confinement time a plasma needs in order to ignite.
    Lawson(LawsonArgs),
}

#[derive(Debug, Args)]
struct MachineArgs {
    /// Major radius `[m]`.
    #[arg(long, default_value_t = 6.2)]
    major_radius: f64,
    /// Minor radius `[m]`.
    #[arg(long, default_value_t = 2.0)]
    minor_radius: f64,
    /// Plasma elongation.
    #[arg(long, default_value_t = 1.7)]
    elongation: f64,
    /// Plasma current `[MA]`.
    #[arg(long, default_value_t = 15.0)]
    plasma_current: f64,
    /// Toroidal field on axis `[T]`.
    #[arg(long, default_value_t = 5.3)]
    toroidal_field: f64,
}

impl MachineArgs {
    fn geometry(&self) -> fusion_core::Result<TokamakGeometry> {
        TokamakGeometry::new(
            self.major_radius,
            self.minor_radius,
            self.elongation,
            self.plasma_current,
            self.toroidal_field,
        )
    }
}

#[derive(Debug, Args)]
struct BurnArgs {
    #[command(flatten)]
    machine: MachineArgs,
    /// Fuel ion density `[m^-3]`.
    #[arg(long, default_value_t = 1.0e20)]
    density: f64,
    /// Fraction of the fuel ions that are deuterium.
    #[arg(long, default_value_t = 0.5)]
    deuterium_fraction: f64,
    /// Effective charge of the plasma.
    #[arg(long, default_value_t = 1.5)]
    z_eff: f64,
    /// Multiplier on the ITER98(y,2) confinement time.
    #[arg(long, default_value_t = 1.0)]
    h_factor: f64,
    /// Auxiliary heating power `[MW]`.
    #[arg(long, default_value_t = 50.0)]
    aux_power: f64,
    /// Temperature at the start of the burn `[keV]`.
    #[arg(long, default_value_t = 5.0)]
    initial_temperature: f64,
    /// Duration of the burn `[s]`.
    #[arg(long, default_value_t = 60.0)]
    duration: f64,
    /// Integration step `[s]`.
    #[arg(long, default_value_t = 0.05)]
    step: f64,
    /// Number of rows to report.
    #[arg(long, default_value_t = 12)]
    rows: usize,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Text)]
    format: Format,
}

#[derive(Debug, Args)]
struct ReactivityArgs {
    /// Reaction channel.
    #[arg(long, value_enum, default_value_t = ReactionArg::Dt)]
    reaction: ReactionArg,
    /// Lowest temperature of the scan `[keV]`.
    #[arg(long, default_value_t = 1.0)]
    min_temperature: f64,
    /// Highest temperature of the scan `[keV]`.
    #[arg(long, default_value_t = 100.0)]
    max_temperature: f64,
    /// Number of points in the scan.
    #[arg(long, default_value_t = 20)]
    points: usize,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Text)]
    format: Format,
}

#[derive(Debug, Args)]
struct LawsonArgs {
    /// Fuel ion density `[m^-3]`.
    #[arg(long, default_value_t = 1.0e20)]
    density: f64,
    /// Effective charge of the plasma.
    #[arg(long, default_value_t = 1.0)]
    z_eff: f64,
    /// Lowest temperature of the scan `[keV]`.
    #[arg(long, default_value_t = 4.0)]
    min_temperature: f64,
    /// Highest temperature of the scan `[keV]`.
    #[arg(long, default_value_t = 50.0)]
    max_temperature: f64,
    /// Number of points in the scan.
    #[arg(long, default_value_t = 12)]
    points: usize,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Text)]
    format: Format,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Format {
    /// Human readable table.
    Text,
    /// Machine readable JSON.
    Json,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum ReactionArg {
    /// T(d,n)4He
    Dt,
    /// D(d,p)T
    DdProton,
    /// D(d,n)3He
    DdNeutron,
    /// 3He(d,p)4He
    DHe3,
}

impl From<ReactionArg> for Reaction {
    fn from(value: ReactionArg) -> Self {
        match value {
            ReactionArg::Dt => Reaction::DeuteriumTritium,
            ReactionArg::DdProton => Reaction::DeuteriumDeuteriumProton,
            ReactionArg::DdNeutron => Reaction::DeuteriumDeuteriumNeutron,
            ReactionArg::DHe3 => Reaction::DeuteriumHelium3,
        }
    }
}

#[derive(Debug, Serialize)]
struct BurnReport<'a> {
    scenario: &'a BurnSimulation,
    plasma_volume_m3: f64,
    peak_fusion_power_mw: f64,
    fusion_energy_gj: f64,
    ignition_time_s: Option<f64>,
    final_q_factor: Option<f64>,
    samples: Vec<&'a BurnSample>,
}

#[derive(Debug, Serialize)]
struct ReactivityPoint {
    temperature_kev: f64,
    reactivity_m3_s: f64,
}

#[derive(Debug, Serialize)]
struct ReactivityReport {
    reaction: Reaction,
    peak_temperature_kev: f64,
    points: Vec<ReactivityPoint>,
}

#[derive(Debug, Serialize)]
struct LawsonPoint {
    temperature_kev: f64,
    ignition_confinement_time_s: Option<f64>,
    ignition_triple_product: Option<f64>,
    fusion_power_density_w_m3: f64,
    bremsstrahlung_density_w_m3: f64,
}

#[derive(Debug, Serialize)]
struct LawsonReport {
    density_m3: f64,
    z_eff: f64,
    points: Vec<LawsonPoint>,
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Burn(args) => run_burn(&args),
        Command::Reactivity(args) => run_reactivity(&args),
        Command::Lawson(args) => run_lawson(&args),
    }
}

fn scan(min: f64, max: f64, points: usize) -> Vec<f64> {
    if points <= 1 {
        return vec![min];
    }
    let step = (max - min) / (points - 1) as f64;
    (0..points).map(|index| min + step * index as f64).collect()
}

fn run_burn(args: &BurnArgs) -> Result<()> {
    let scenario = BurnSimulation {
        geometry: args
            .machine
            .geometry()
            .context("invalid machine geometry")?,
        ion_density_m3: args.density,
        deuterium_fraction: args.deuterium_fraction,
        z_eff: args.z_eff,
        ion_mass_amu: fusion_core::constants::DT_MEAN_MASS_AMU,
        h_factor: args.h_factor,
        auxiliary_power_mw: args.aux_power,
        initial_temperature_kev: args.initial_temperature,
    };

    let history = scenario
        .run(args.duration, args.step)
        .context("burn simulation failed")?;

    let report = BurnReport {
        scenario: &scenario,
        plasma_volume_m3: scenario.geometry.plasma_volume_m3(),
        peak_fusion_power_mw: history.peak_fusion_power_mw(),
        fusion_energy_gj: history.fusion_energy_gj(),
        ignition_time_s: history.ignition_time_s(),
        final_q_factor: history.final_sample().q_factor(),
        samples: history.downsample(args.rows),
    };

    match args.format {
        Format::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        Format::Text => print_burn_text(&scenario, &history, &report),
    }
    Ok(())
}

fn print_burn_text(scenario: &BurnSimulation, history: &BurnHistory, report: &BurnReport<'_>) {
    println!(
        "Machine   R = {:.2} m, a = {:.2} m, kappa = {:.2}, Ip = {:.1} MA, Bt = {:.1} T",
        scenario.geometry.major_radius_m,
        scenario.geometry.minor_radius_m,
        scenario.geometry.elongation,
        scenario.geometry.plasma_current_ma,
        scenario.geometry.toroidal_field_t,
    );
    println!(
        "Plasma    n = {:.2e} m^-3, Zeff = {:.2}, H = {:.2}, volume = {:.0} m^3",
        scenario.ion_density_m3, scenario.z_eff, scenario.h_factor, report.plasma_volume_m3,
    );
    println!(
        "Heating   {:.1} MW auxiliary, starting from {:.1} keV\n",
        scenario.auxiliary_power_mw, scenario.initial_temperature_kev,
    );

    println!(
        "{:>8}  {:>9}  {:>9}  {:>10}  {:>10}  {:>10}  {:>8}",
        "t [s]", "T [keV]", "tau_E [s]", "P_fus [MW]", "P_a [MW]", "P_loss[MW]", "Q"
    );
    for sample in &report.samples {
        println!(
            "{:>8.2}  {:>9.2}  {:>9.2}  {:>10.1}  {:>10.1}  {:>10.1}  {:>8}",
            sample.time_s,
            sample.temperature_kev,
            sample.confinement_time_s,
            sample.power_mw.fusion,
            sample.power_mw.alpha_heating,
            sample.power_mw.loss_power(),
            sample
                .q_factor()
                .map(|q| format!("{q:.2}"))
                .unwrap_or_else(|| "inf".to_string()),
        );
    }

    let final_sample = history.final_sample();
    println!(
        "\nPeak fusion power   {:.1} MW\nFusion energy       {:.2} GJ\nFinal triple product {:.2e} keV s m^-3",
        report.peak_fusion_power_mw, report.fusion_energy_gj, final_sample.triple_product,
    );
    match report.ignition_time_s {
        Some(time) => println!("Ignition reached at {time:.2} s"),
        None => println!("Ignition not reached"),
    }
}

fn run_reactivity(args: &ReactivityArgs) -> Result<()> {
    let reaction = Reaction::from(args.reaction);
    let points = scan(args.min_temperature, args.max_temperature, args.points)
        .into_iter()
        .map(|temperature_kev| ReactivityPoint {
            temperature_kev,
            reactivity_m3_s: reactivity(reaction, temperature_kev),
        })
        .collect();

    let report = ReactivityReport {
        reaction,
        peak_temperature_kev: peak_temperature_kev(reaction),
        points,
    };

    match args.format {
        Format::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        Format::Text => {
            println!(
                "Reaction {} (peak at {:.1} keV)\n",
                reaction.label(),
                report.peak_temperature_kev
            );
            println!("{:>10}  {:>14}", "T [keV]", "<sigma v> [m^3/s]");
            for point in &report.points {
                println!(
                    "{:>10.2}  {:>14.4e}",
                    point.temperature_kev, point.reactivity_m3_s
                );
            }
        }
    }
    Ok(())
}

fn run_lawson(args: &LawsonArgs) -> Result<()> {
    let mut points = Vec::with_capacity(args.points);
    for temperature_kev in scan(args.min_temperature, args.max_temperature, args.points) {
        let plasma = PlasmaState::new(args.density, temperature_kev, 0.5, args.z_eff)
            .context("invalid plasma parameters")?;
        let confinement = plasma.ignition_confinement_time_s();
        points.push(LawsonPoint {
            temperature_kev,
            ignition_confinement_time_s: confinement,
            ignition_triple_product: confinement.map(|tau| plasma.triple_product(tau)),
            fusion_power_density_w_m3: plasma.fusion_power_density(),
            bremsstrahlung_density_w_m3: plasma.bremsstrahlung_density(),
        });
    }

    let report = LawsonReport {
        density_m3: args.density,
        z_eff: args.z_eff,
        points,
    };

    match args.format {
        Format::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        Format::Text => {
            println!(
                "Ignition requirement at n = {:.2e} m^-3, Zeff = {:.2}\n",
                report.density_m3, report.z_eff
            );
            println!(
                "{:>10}  {:>12}  {:>18}  {:>14}",
                "T [keV]", "tau_E [s]", "n T tau [keV s/m^3]", "P_fus [W/m^3]"
            );
            for point in &report.points {
                let tau = point
                    .ignition_confinement_time_s
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "-".to_string());
                let triple = point
                    .ignition_triple_product
                    .map(|value| format!("{value:.3e}"))
                    .unwrap_or_else(|| "-".to_string());
                println!(
                    "{:>10.2}  {:>12}  {:>18}  {:>14.3e}",
                    point.temperature_kev, tau, triple, point.fusion_power_density_w_m3
                );
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn command_line_definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn burn_defaults_describe_the_iter_baseline() {
        let cli = Cli::try_parse_from(["fusion-sim", "burn"]).unwrap();
        let Command::Burn(args) = cli.command else {
            panic!("expected the burn subcommand");
        };
        assert_eq!(args.machine.major_radius, 6.2);
        assert_eq!(args.density, 1.0e20);
        assert_eq!(args.format, Format::Text);
    }

    #[test]
    fn burn_accepts_overrides() {
        let cli = Cli::try_parse_from([
            "fusion-sim",
            "burn",
            "--aux-power",
            "80",
            "--duration",
            "10",
            "--format",
            "json",
        ])
        .unwrap();
        let Command::Burn(args) = cli.command else {
            panic!("expected the burn subcommand");
        };
        assert_eq!(args.aux_power, 80.0);
        assert_eq!(args.duration, 10.0);
        assert_eq!(args.format, Format::Json);
    }

    #[test]
    fn reaction_arguments_map_onto_the_model() {
        assert_eq!(Reaction::from(ReactionArg::Dt), Reaction::DeuteriumTritium);
        assert_eq!(
            Reaction::from(ReactionArg::DHe3),
            Reaction::DeuteriumHelium3
        );
    }

    #[test]
    fn unknown_subcommands_are_rejected() {
        assert!(Cli::try_parse_from(["fusion-sim", "melt"]).is_err());
    }

    #[test]
    fn scan_covers_the_requested_interval() {
        let values = scan(1.0, 10.0, 4);
        assert_eq!(values.len(), 4);
        assert_eq!(values[0], 1.0);
        assert!((values[3] - 10.0).abs() < 1e-12);
        assert_eq!(scan(3.0, 9.0, 1), vec![3.0]);
        assert_eq!(scan(3.0, 9.0, 0), vec![3.0]);
    }

    #[test]
    fn a_short_burn_runs_end_to_end() {
        let args = BurnArgs {
            machine: MachineArgs {
                major_radius: 6.2,
                minor_radius: 2.0,
                elongation: 1.7,
                plasma_current: 15.0,
                toroidal_field: 5.3,
            },
            density: 1.0e20,
            deuterium_fraction: 0.5,
            z_eff: 1.5,
            h_factor: 1.0,
            aux_power: 50.0,
            initial_temperature: 5.0,
            duration: 2.0,
            step: 0.1,
            rows: 5,
            format: Format::Json,
        };
        run_burn(&args).unwrap();
    }

    #[test]
    fn an_invalid_geometry_is_reported() {
        let args = BurnArgs {
            machine: MachineArgs {
                major_radius: 0.0,
                minor_radius: 2.0,
                elongation: 1.7,
                plasma_current: 15.0,
                toroidal_field: 5.3,
            },
            density: 1.0e20,
            deuterium_fraction: 0.5,
            z_eff: 1.5,
            h_factor: 1.0,
            aux_power: 50.0,
            initial_temperature: 5.0,
            duration: 2.0,
            step: 0.1,
            rows: 5,
            format: Format::Text,
        };
        assert!(run_burn(&args).is_err());
    }
}
