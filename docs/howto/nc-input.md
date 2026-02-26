# NC Input Module

The NC input module reads subdaily, CF-compliant NetCDF forcing data into
LPJ-GUESS. It is the standard input module for all DAVE simulation workflows.

## Background

The NC module is loosely based on the `cf` input module used by upstream
LPJ-GUESS, but adds support for **subdaily meteorological inputs**. Subdaily
inputs are required by the DAVE photosynthesis module
<!-- TODO: link to photosynthesis science docs when available -->, which
operates at sub-daily time resolution.

## Selecting the Module

Pass `-input nc` on the command line:

```bash
./build/guess -input nc benchmarks/ozflux/Howard_Springs/Howard_Springs.ins
```

Or in the context of a parallel run:

```bash
mpirun -np 8 ./build/guess -parallel -input nc <insfile>
```

## Required Instruction File Parameters

The following parameters must be present in the `.ins` file. See
`modules/ncinput.h` for the full parameter reference.

### Simulation resolution

```
nsteps_day <N>     ! Subdaily timesteps per day (e.g. 48 for half-hourly)
```

### Input files

```
file_met_forcing  "path/to/forcing.nc"    ! CF-compliant NetCDF met file
file_gridlist     "path/to/gridlist.txt"  ! lon lat pairs, one per line
file_soildata     "path/to/soil.nc"
file_ndep         "path/to/ndep.nc"
ndep_scenario     HISTORIC                ! HISTORIC | RCP26 | RCP45 | RCP60 | RCP85
```

### NetCDF variable names

The module does not assume fixed variable names — each must be declared
explicitly to match your input file:

```
variable_lon   "longitude"
variable_lat   "latitude"
variable_time  "time"

variable_temp   "Tair"       ! Mean temperature (°C)
variable_insol  "SWdown"     ! Downwelling shortwave (W/m²)
variable_precip "Precip"     ! Precipitation (mm or kg/m²)
```

### Humidity (one required)

```
variable_vpd  "VPD"    ! Vapour pressure deficit (kPa)
! — or —
variable_rh   "RH"     ! Relative humidity (dimensionless, 0–1)
```

If `variable_rh` is used, also set the saturation vapour pressure method:

```
esat_method sonntag_1990   ! sonntag_1990 | alduchov_1996 | allen_1998
```

### Atmospheric pressure (one required)

```
variable_ps   "PSurf"    ! Atmospheric pressure (Pa)
! — or —
fixed_pressure 101325     ! Fixed constant (Pa)
```

### CO₂ (one required)

```
variable_co2  "CO2air"          ! CO₂ variable in the NetCDF file (ppm)
! — or —
file_co2      "path/to/co2.txt" ! Text file: year CO2_ppm, one row per year
```

### Performance

```
cache_strategy 1   ! 0 = whole file (fastest), 1 = one year at a time, 2 = no cache
```

## NetCDF Input Format

Input files must follow **CF Conventions 1.6**:

- A `time` coordinate variable with `units` (e.g. `"days since 1850-01-01"`)
  and `calendar` attributes.
- Longitude and latitude coordinates identifiable via `standard_name`,
  `axis`, or CF auxiliary coordinates.
- `scale_factor` / `add_offset` and `missing_value` attributes are supported.

The input timestep must be greater than or equal to `nsteps_day` and must
divide evenly into it. For example, if the forcing file is 3-hourly (8
steps/day) and the simulation requests 48 steps/day, the module will
interpolate accordingly.
