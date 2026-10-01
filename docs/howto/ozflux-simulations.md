# OzFlux Site Simulations

This page describes the input file structure for the OzFlux site-level
simulations included in the repository under `benchmarks/ozflux/`.

## Directory Structure

```
benchmarks/ozflux/
├── settings.ins          # Top-level simulation switches
├── outputs.ins           # Output file declarations
├── pft_dave.ins          # PFT parameters (trees, grasses)
├── soiln.ins             # Soil nitrogen parameters
├── forcing/              # Forcing data (climate, soil, fire, etc)
├── <SiteName>/
│   ├── <SiteName>.ins    # Site-level instruction file (entry point)
│   └── gridlist.txt      # Lon/lat for this site
└── ...
```

Each site has its own subdirectory containing a site-level `.ins` file and a
`gridlist.txt`. The site-level `.ins` file is short: it sets the site-specific
paths and imports the shared configuration files.

## Shared Configuration Files

All sites import the same set of common `.ins` files from the `benchmarks/ozflux/`
directory. These control parameters shared across all sites:

| File | Contents |
|---|---|
| `settings.ins` | Top-level switches: timestep, enabled modules, scenario |
| `outputs.ins` | Available output variables; most are disabled by default |
| `pft_dave.ins` | PFT parameters for trees and grasses |
| `soiln.ins` | Soil nitrogen parameters |

To enable an output variable, open `outputs.ins` and uncomment the relevant
line. The comment character is `!`.

## Obtaining Climate Data

The forcing data is stored on HIE-Storage and must be downloaded separately.
See [HIE-Storage](hie-storage.md) for access details. The data should be
placed at `benchmarks/ozflux/forcing/`.

## Running a Site

For a complete worked example, start with [Running and Debugging via
VSCode](debugging-vscode.md). That page shows how to run `CumberlandPlain`,
where outputs are written, and what to check after the run. To run from the
command line, see [Running Serially](run-serial.md).

Outputs are written to `benchmarks/ozflux/<SiteName>/out/` by default.
The model log is `benchmarks/ozflux/<SiteName>/guess.log`. Check it first if
the run fails or an output looks wrong. `verbosity` in the site `.ins` file
controls its detail; the shared `settings.ins` sets it to `3`. See
[Start with `guess.log`](common-failures.md#start-with-guesslog) for the
recommended diagnostic steps.
