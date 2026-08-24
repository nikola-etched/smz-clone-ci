# fusion

A zero-dimensional power balance model for magnetically confined deuterium-tritium
plasmas, built with Bazel.

The model combines three standard pieces of tokamak physics:

- **Fusion reactivities** from the Bosch-Hale parameterization (Nuclear Fusion **32**
  (1992) 611) for D-T, both D-D branches and D-3He.
- **Power balance** of a uniform plasma: fusion power, alpha heating, neutron losses,
  bremsstrahlung radiation, transport losses and the Lawson triple product.
- **Energy confinement** from the ITER98(y,2) H-mode scaling law, coupled to the plasma
  through the heating power.

These are integrated in time with a fixed-step fourth-order Runge-Kutta scheme to follow a
plasma as it heats up, burns and settles into a steady state.

## Layout

| Path | Contents |
| --- | --- |
| `crates/fusion-core` | The physics library: reactivities, power balance, confinement, burn integration |
| `crates/fusion-sim` | Command line front end |

## Building and testing

The build is hermetic: Bazel fetches its own Rust toolchain, so nothing beyond
[Bazelisk](https://github.com/bazelbuild/bazelisk) is required.

```console
$ bazel build //...
$ bazel test //...
```

Formatting and lints run as Bazel aspects over the same targets:

```console
$ bazel build --config=rustfmt //...   # rustfmt --check
$ bazel build --config=clippy //...    # clippy, warnings are errors
$ bazel run @rules_rust//tools/rustfmt # apply formatting
```

Third-party crates are declared in `Cargo.toml` and turned into Bazel targets by
`crate_universe`. After changing a dependency, re-pin the lockfile:

```console
$ CARGO_BAZEL_REPIN=1 bazel build //...
```

The Cargo manifests are kept working, so `cargo test --workspace` also runs the suite for
contributors who prefer it. Bazel is the source of truth used by CI.

## Running a simulation

Integrate the ITER baseline scenario:

```console
$ bazel run -- //crates/fusion-sim burn --duration 60 --rows 8
Machine   R = 6.20 m, a = 2.00 m, kappa = 1.70, Ip = 15.0 MA, Bt = 5.3 T
Plasma    n = 1.00e20 m^-3, Zeff = 1.50, H = 1.00, volume = 832 m^3
Heating   50.0 MW auxiliary, starting from 5.0 keV

   t [s]    T [keV]  tau_E [s]  P_fus [MW]    P_a [MW]  P_loss[MW]         Q
    0.00       5.00       5.18        80.1        16.0        53.5      1.60
    ...
   60.00       8.67       2.84       458.2        91.7       141.7      9.16

Peak fusion power   458.2 MW
Fusion energy       23.07 GJ
Final triple product 2.47e21 keV s m^-3
```

Other subcommands tabulate the reactivity against temperature and scan the confinement
time a plasma needs in order to ignite:

```console
$ bazel run -- //crates/fusion-sim reactivity --reaction dt --points 10
$ bazel run -- //crates/fusion-sim lawson --density 1e20
```

Every subcommand takes `--format json` for machine-readable output.

## Continuous integration

`.github/workflows/ci.yaml` runs on GitHub-hosted runners only:

- `bazel build //...` and `bazel test //...` on `ubuntu-latest` and `macos-latest`
- rustfmt and clippy aspects on `ubuntu-latest`

Bazelisk downloads, the repository cache and the disk cache are all cached between runs.

## License

Apache License 2.0. See [LICENSE](LICENSE).
