# Climate Input

The DAVE input module is the standard input module for DAVE simulations. It
reads subdaily, CF-compliant NetCDF forcing data into LPJ-GUESS.

## Background

The DAVE climate input stack is loosely based on the `cf` input module used by
upstream LPJ-GUESS, but adds support for **subdaily meteorological inputs**.
Subdaily inputs are required by the [DAVE photosynthesis module](../science/photosynthesis.md),
which operates at sub-daily time resolution.

## Input Architecture

At a high level:

- `DAVEInput` is the top-level input module selected by `-input nc`.
- `NCInput` is the climate-specific submodule used by `DAVEInput`.
- `NCInput` discovers required variables from NetCDF metadata (`standard_name`
  + `units`) and reads the corresponding timeseries.

## Selecting the Module

Pass `-input nc` on the command line:

```bash
./build/guess -input nc myfile.ins
```

Or in the context of a parallel run:

```bash
mpirun -np 8 ./build/guess -parallel -input nc <insfile>
```

`nc` is the current CLI key for DAVE climate input. A future breaking release
may introduce `dave` as the module key.

## Instruction File Parameters

All parameters are declared in `modules/daveinput.cpp`. Parameters with no
default listed are required and must be set explicitly.

### Simulation resolution

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `nsteps_day` | int | 1–86400 | — | Number of subdaily timesteps per day (e.g. 48 for half-hourly) |

### Input file paths

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file_met_forcing` | string | — | Path or shell-glob pattern for CF-compliant NetCDF met forcing files |
| `file_met_spinup` | string | *(same as forcing)* | Path or shell-glob pattern for NetCDF met data used during spinup. If omitted, `file_met_forcing` is used for both spinup and simulation. |
| `file_gridlist` | string | — | Path to gridlist file (one `lon lat` pair per line) |
| `file_soildata` | string | — | Path to soil input file |
| `file_ndep` | string | — | Path to nitrogen deposition archive file (used when N deposition is not read from NetCDF variables) |
| `file_co2` | string | — | Path to annual CO₂ text file (format: `year ppm`, one row per year). Used when CO₂ is not read from NetCDF |
| `file_thome` | string | — | Path to mean maximum temperature ($T_{home}$) input file. Required when `iftempacclimation` is enabled. |

### Humidity method

If humidity is provided as relative humidity (rather than direct VPD), specify
the saturation vapour pressure method:

| Parameter | Type | Options | Description |
|-----------|------|---------|-------------|
| `esat_method` | string | `sonntag_1990`, `alduchov_1996`, `allen_1998` | Method for computing saturation vapour pressure from temperature |

### Atmospheric pressure

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `fixed_pressure` | float | 0–200000 | Fixed constant atmospheric pressure to use (Pa). Default: 0 (disabled). |

### CO₂ concentration

| Parameter | Type | Description |
|-----------|------|-------------|
| `file_co2` | string | Path to a text file of annual CO₂ values (see above). If omitted, CO₂ may be read from NetCDF where present. |

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_mean_annual_co2` | bool | `false` | If true, each timestep uses the annual mean CO₂; otherwise the subdaily values from the file are averaged per day |

### Nitrogen deposition

| Parameter | Type | Options | Description |
|-----------|------|---------|-------------|
| `ndep_scenario` | string | `HISTORIC`, `RCP26`, `RCP45`, `RCP60`, `RCP85` | Scenario for the nitrogen deposition archive (`file_ndep`) |

### Spinup

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `spinup_met_cycle_length` | int | 0–1000000 | 0 | Number of years of met data to draw from during spinup. 0 means all available years are used. |

### Perturbations

These add a constant offset to the corresponding variable at every timestep.
Useful for sensitivity experiments.

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `temperature_delta` | float | 0–1×10⁶ | 0 | Offset added to temperature at every timestep (°C) |
| `co2_delta` | float | 0–1×10⁶ | 0 | Offset added to CO₂ at every timestep (ppm) |

### Performance and diagnostics

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cache_strategy` | int (0–2) | — | Caching strategy: `0` = cache entire file (fastest, highest memory), `1` = cache one year at a time (recommended for large files), `2` = no caching (slowest) |
| `nc_check_inputs` | bool | — | If enabled, performs sanity checks on met values at each timestep |
| `progress_interval` | int | -1 | Interval in seconds between progress log lines. `0` = yearly; `-1` = daily |
| `progress_newline` | bool | `false` | Terminate progress messages with a newline (recommended for MPI runs) |
| `strict_math` | bool | — | Treat invalid floating-point operations (e.g. division by zero) as fatal errors |

---

## NetCDF Variable Discovery and Metadata Requirements

`NCInput` discovers variables from NetCDF metadata. For each variable it uses,
the NetCDF variable must include:

- `standard_name`
- `units`

Variables without both attributes are ignored by the resolver.

### Required coordinate metadata

Coordinate variables must be discoverable with these `standard_name` values:

- `longitude`
- `latitude`
- `time`

The `time` coordinate must also include CF-compatible `units` (for example
`days since 1850-01-01`) and `calendar`.

### Climate variables and accepted `standard_name` values

The table below lists the names currently recognised by DAVE climate input.

| Input quantity | `standard_name` values |
|---|---|
| Air temperature (required) | `air_temperature` |
| Shortwave radiation (required) | `surface_downwelling_shortwave_flux_in_air`, `surface_downwelling_shortwave_flux` |
| Precipitation (required) | `thickness_of_rainfall_amount`, `precipitation_amount`, `precipitation_flux` |
| Humidity (one required source) | VPD: `water_vapor_saturation_deficit_in_air` or `water_vapor_saturation_deficit`; OR relative humidity: `relative_humidity`; OR specific humidity: `specific_humidity` |
| Surface pressure (required unless `fixed_pressure` > 0) | `surface_air_pressure` |
| CO₂ (required if `file_co2` not provided) | `mole_fraction_of_carbon_dioxide_in_air` |
| Wind speed (required when BLAZE is enabled) | `wind_speed` |
| Daily Tmax/Tmin (required for daily forcing with BLAZE) | `maximum_air_temperature`, `minimum_air_temperature` |
| NetCDF N deposition (optional alternative to archive input) | `ndep_nhx`, `ndep_noy` |

### Units

Units are validated when climate variables are opened. In practice, datasets
should use these units:

- Temperature: `degC` / `°C` / `℃`
- Shortwave radiation: `W/m2` or `W m-2`
- Precipitation: `mm`, `kg m-2`, or `kg/m2`
- VPD: `kPa`
- Relative humidity: `1`
- Specific humidity: `1` or `kg/kg`
- Pressure: `Pa`
- CO₂: `ppm`
- Wind speed: `m/s` or `m s-1`
- NetCDF N deposition: `kgN/m2/day`

### Uniqueness rule

Within a resolved input set, each required `standard_name` must map to a
single variable. Ambiguous duplicates for the same `standard_name` are treated
as an error.

## Spinup Climate Cycling

During spinup (`date.year < nyear_spinup`), the module maps each spinup year
to a calendar year of met data by cycling deterministically through the
available years. The goal is for the met years used in spinup to arrive at
the start of the main simulation data as naturally as possible — the year
immediately before the simulation data begins is always the last full year
of the spinup dataset.

### How years are selected

The module cycles through the spinup dataset in order, wrapping back to the
start when it reaches the end. The cycling is deterministic — the same run
configuration always produces the same sequence.

The wrap-around point is chosen so that the sequence ends exactly at the
transition into the simulation data. Once the calendar year falls within the
range covered by the spinup dataset, actual data years are used directly
with no cycling.

**Example** — suppose the spinup dataset covers 1990–1999 (10 years), the
spinup runs for 30 years, and the main simulation begins in 2000:

| Spinup year | Calendar year | Met data year | |
|-------------|---------------|---------------|-|
| 0 | 1970 | 1990 | cycling begins |
| 1 | 1971 | 1991 | |
| ... | ... | ... | |
| 8 | 1978 | 1998 | |
| 9 | 1979 | 1999 | end of cycle 1 |
| 10 | 1980 | 1990 | cycle 2 begins |
| 11 | 1981 | 1991 | |
| ... | ... | ... | |
| 19 | 1989 | 1999 | end of cycle 2 |
| 20 | 1990 | 1990 | within data range — direct |
| ... | ... | ... | |
| 28 | 1998 | 1998 | | |
| 29 | 1999 | 1999 | last spinup year |
| — | 2000 | 2000 | simulation begins |

The data cycles twice, and the transition into the simulation data at year
2000 occurs naturally at the end of the final cycle year.

If `spinup_met_cycle_length` is set to a value greater than zero, only the
**first N years** of the spinup dataset are used for cycling. This is useful
when only a short pre-industrial baseline period should be repeated during
spinup, even if a longer dataset is available.

### CO₂ is not cycled

CO₂ concentration is always taken from the actual calendar year being
simulated. It is not cycled with the meteorology, even during spinup.

### Separate spinup file

A separate NetCDF file can be specified via `file_met_spinup` for the spinup
period. This is useful when a longer historical baseline record exists
(e.g. a reanalysis product covering a pre-industrial period) that differs
from the main simulation forcing. If `file_met_spinup` is not set, or is
identical to `file_met_forcing`, the same file is used for both periods.

When a separate spinup file is used, the module opens the spinup file first,
determines the available years and timestep, and then switches to the main
forcing file at the start of the transient simulation.

---

## Multi-file Input Paths

`file_met_forcing` and `file_met_spinup` accept either:

- a single NetCDF path, or
- a shell-glob pattern (`*` and `?`) that resolves to multiple NetCDF files.

Examples:

```
file_met_forcing "/data/barra2/forcing_*.nc"
file_met_spinup  "/data/barra2/spinup_19??.nc"
```

Matched files are sorted and then scanned for variables using metadata
resolution (`standard_name` + `units`).

---

## NetCDF Input Format

Input files must follow **CF Conventions 1.6**:

- A `time` coordinate variable with `units` (e.g. `"days since 1850-01-01"`)
  and `calendar` attributes.
- Coordinate variables with `standard_name` values `longitude`, `latitude`,
  and `time`.
- `scale_factor` / `add_offset` and `missing_value` attributes are not
  supported.

The input timestep (steps/day in the forcing files) must be greater than or
equal to `nsteps_day` and must be an integer multiple of `nsteps_day`. For
example, 48 input steps/day can be aggregated to 24 requested steps/day.
