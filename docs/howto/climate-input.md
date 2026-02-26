# Climate Input

The NC input module is the standard input module for all DAVE simulations. It
reads subdaily, CF-compliant NetCDF forcing data into LPJ-GUESS.

## Background

The NC module is loosely based on the `cf` input module used by upstream
LPJ-GUESS, but adds support for **subdaily meteorological inputs**. Subdaily
inputs are required by the [DAVE photosynthesis module](../science/photosynthesis.md),
which operates at sub-daily time resolution.

## Selecting the Module

Pass `-input nc` on the command line:

```bash
./build/guess -input nc myfile.ins
```

Or in the context of a parallel run:

```bash
mpirun -np 8 ./build/guess -parallel -input nc <insfile>
```

## Instruction File Parameters

All parameters are declared in `modules/ncinput.cpp`. Parameters with no
default listed are required and must be set explicitly.

### Simulation resolution

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `nsteps_day` | int | 1–86400 | — | Number of subdaily timesteps per day (e.g. 48 for half-hourly) |

### Input file paths

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file_met_forcing` | string | — | Path to CF-compliant NetCDF met forcing file |
| `file_met_spinup` | string | *(same as forcing)* | Path to NetCDF met file used during spinup. If omitted, `file_met_forcing` is used for both spinup and simulation. |
| `file_gridlist` | string | — | Path to gridlist file (one `lon lat` pair per line) |
| `file_soildata` | string | — | Path to soil input file |
| `file_ndep` | string | — | Path to nitrogen deposition archive file. Not required if `variable_ndep_nhx`/`variable_ndep_noy` are set. |
| `file_co2` | string | — | Path to annual CO₂ text file (format: `year ppm`, one row per year). Required if `variable_co2` is not set. |
| `file_thome` | string | — | Path to mean maximum temperature ($T_{home}$) input file. Required when `iftempacclimation` is enabled. |

### NetCDF variable names

The module does not assume fixed variable names — each must be declared
explicitly to match your NetCDF file.

**Coordinate variables (required):**

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_lon` | string | Name of the longitude variable |
| `variable_lat` | string | Name of the latitude variable |
| `variable_time` | string | Name of the time variable |

**Meteorological variables (required):**

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_temp` | string | Mean air temperature (°C) |
| `variable_insol` | string | Downwelling shortwave radiation (W m⁻²) |
| `variable_precip` | string | Precipitation (mm or kg m⁻² s⁻¹) |

**Meteorological variables (optional):**

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_wind` | string | Wind speed (m s⁻¹) |
| `variable_tmax` | string | Daily maximum air temperature (°C) |
| `variable_tmin` | string | Daily minimum air temperature (°C) |

### Humidity (one required)

Provide one of the following:

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_vpd` | string | Vapour pressure deficit (kPa) |
| `variable_rh` | string | Relative humidity (dimensionless, 0–1) |

If `variable_rh` is used, also specify the saturation vapour pressure method:

| Parameter | Type | Options | Description |
|-----------|------|---------|-------------|
| `esat_method` | string | `sonntag_1990`, `alduchov_1996`, `allen_1998` | Method for computing saturation vapour pressure from temperature |

### Atmospheric pressure (one required)

Provide one of the following:

| Parameter | Type | Range | Description |
|-----------|------|-------|-------------|
| `variable_ps` | string | — | Name of the surface pressure variable in the NetCDF file (Pa) |
| `fixed_pressure` | float | 0–200000 | Fixed constant atmospheric pressure to use (Pa). Default: 0 (disabled). |

### CO₂ concentration (one required)

Provide one of the following:

| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_co2` | string | Name of the CO₂ variable in the NetCDF file (ppm) |
| `file_co2` | string | Path to a text file of annual CO₂ values (see above) |

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_mean_annual_co2` | bool | `false` | If true, each timestep uses the annual mean CO₂; otherwise the subdaily values from the file are averaged per day |

### Nitrogen deposition

| Parameter | Type | Options | Description |
|-----------|------|---------|-------------|
| `ndep_scenario` | string | `HISTORIC`, `RCP26`, `RCP45`, `RCP60`, `RCP85` | Scenario for the nitrogen deposition archive (`file_ndep`) |
| `variable_ndep_nhx` | string | — | Name of NHx deposition variable in the NetCDF file (alternative to `file_ndep`) |
| `variable_ndep_noy` | string | — | Name of NOy deposition variable in the NetCDF file (alternative to `file_ndep`) |

### Spinup

| Parameter | Type | Range | Default | Description |
|-----------|------|-------|---------|-------------|
| `spinup_met_cycle_length` | int | 0–1000000 | 0 | Number of years of met data to draw from during spinup. 0 means all available years are used. |
| `met_sampling_seed` | int | 0–INT_MAX | — | Random seed for met-year sampling during spinup. Use the same seed to reproduce identical spinup sequences. |

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

## Spinup Climate Cycling

During spinup (`date.year < nyear_spinup`), the module randomly samples
calendar years of meteorological data rather than replaying a fixed sequence.
This avoids introducing spurious trends from repeatedly cycling the same
short period at the start of available data.

### How years are sampled

At the start of each spinup year, a random calendar year is drawn from the
available range in the met file. If `spinup_met_cycle_length` is set to a
value greater than zero, sampling is limited to the first
`spinup_met_cycle_length` years of available data, which is useful when a
short pre-industrial baseline period is available separately from the main
forcing record.

The sampling is deterministic given the same `met_sampling_seed`, so spinup
behaviour is reproducible across runs.

### CO₂ is not cycled

CO₂ concentration is always taken from the actual calendar year being
simulated. It is not randomised or cycled with the meteorology, even during
spinup.

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

## NetCDF Input Format

Input files must follow **CF Conventions 1.6**:

- A `time` coordinate variable with `units` (e.g. `"days since 1850-01-01"`)
  and `calendar` attributes.
- Longitude and latitude coordinates identifiable via `standard_name`,
  `axis`, or CF auxiliary coordinates.
- `scale_factor` / `add_offset` and `missing_value` attributes are not
  supported.

The input timestep must be greater than or equal to `nsteps_day` and must
divide evenly into it. For example, if the forcing file is 3-hourly (8
steps/day) and the simulation requests 48 steps/day, the module will
interpolate accordingly.
