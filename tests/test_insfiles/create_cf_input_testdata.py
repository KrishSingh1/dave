# This code generates the netCDF input files for the cfinput integration test
import xarray as xr
import numpy as np
import pandas as pd

total_days = 11323 # from 1970-01-01 to 2000-12-31 (inclusive leap days), need to have more than 30 years for the spinup to work
lon = [9.5]
lat = [49.5]
time_dt = pd.date_range(start="1970-01-01", periods=total_days, freq="D")
time = np.arange(total_days, dtype=np.int32)
time_unit = 'days since 1970-01-01'

def create_dataset(variable_name, standard_name, unit, data):
    ds = xr.Dataset({variable_name: data})
    ds['time'].attrs['units'] = time_unit
    ds['time'].attrs['calendar'] = 'standard'
    ds['time'].attrs['axis'] = 'T'
    ds['lon'].attrs['axis'] = 'X'
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees_east'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees_north'
    ds['lat'].attrs['axis'] = 'Y'


    ds[variable_name].attrs['units'] = unit
    ds[variable_name].attrs['standard_name'] = standard_name

    return ds

# -------------------------------------------------------------
# Helper functions for the variability
# -------------------------------------------------------------
def seasonal_cycle(doy, amplitude, phase=0):
    """Annual sinusoid (period = 365 d)."""
    return amplitude * np.sin(2 * np.pi * (doy - 1 + phase) / 365.0)

def yearly_offset(year_idx, sigma, seed):
    """Random offset that stays constant for an entire calendar year."""
    rng = np.random.default_rng(seed)
    # one random number per distinct year
    offsets = rng.normal(0.0, sigma, size=year_idx.max() + 1)
    return offsets[year_idx]

# -----------------------------------------------------------------
# Common time‑related arrays (derived from the datetime version)
# -----------------------------------------------------------------
day_of_year = time_dt.dayofyear.values                 # 1 … 365
year_index   = (time_dt.year - time_dt.year[0]).values # 0,1,2,...

# -------------------------------------------------------------
# Temperature (K)
# -------------------------------------------------------------
base_temp   = 285.0          # mean annual temperature
amp_temp    = 10.0           # seasonal amplitude (±10 K)
sigma_ytemp = 2.0            # inter‑annual std‑dev (K)
sigma_dtemp = 0.5            # day‑to‑day noise (K)

temp_seas = seasonal_cycle(day_of_year, amp_temp)
temp_year = yearly_offset(year_index, sigma_ytemp, seed=12345)
temp_noise = np.random.default_rng(2025).normal(0.0, sigma_dtemp, size=total_days)
temp_vals = base_temp + temp_seas + temp_year + temp_noise
# expand `temp_vals` (1D: time) to a 3D array with shape (time, lat, lon)
temp = xr.DataArray(
    temp_vals[:, None, None],
    dims=('time', 'lat', 'lon'),
    coords={'time': time, 'lat': lat, 'lon': lon},
)
create_dataset('temp', 'air_temperature', 'K', temp).to_netcdf('dummy_temp.nc')

# -------------------------------------------------------------
# Precipitation (kg m-2) – precip amount
# -------------------------------------------------------------
# A simple wet‑season that peaks around day 180 (mid‑year)
base_pr    = 1.2e-5           # kg m-2 s-1 -> ≈ 1 kg m⁻² day⁻¹
amp_pr     = 0.9e-5
sigma_ypr  = 5e-6
sigma_dpr  = 2e-6

pr_seas = seasonal_cycle(day_of_year, amp_pr, phase=180)   # shift peak to summer
pr_year = yearly_offset(year_index, sigma_ypr, seed=12346)
pr_noise = np.random.default_rng(2026).normal(0.0, sigma_dpr, size=total_days)

pr_vals = (base_pr + pr_seas + pr_year + pr_noise) * 86400 # from kg m-2 s-1 to kg m-2 day-1
# Ensure precipitation is non-negative (clip small negative values from noise)
pr_vals = np.maximum(pr_vals, 0.0)

# Add random dry days
rng = np.random.default_rng(42)
rain_probability = 0.35  # ~35% chance of rain per day
rain_mask = rng.random(total_days) < rain_probability
pr_vals = pr_vals * rain_mask  # Set dry days to zero

#pr = np.full((len(time), len(lat), len(lon)), 1.0)
pr = xr.DataArray(pr_vals[:, None, None], dims=('time', 'lat', 'lon'), coords={'time': time, 'lat': lat, 'lon': lon})
create_dataset('pr', 'precipitation_amount', 'kg m-2', pr).to_netcdf('dummy_pr.nc')

# -------------------------------------------------------------
# Short‑wave radiation (W m⁻²)
# -------------------------------------------------------------
base_rad    = 200.0
amp_rad     = 80.0
sigma_yrad  = 5.0
sigma_drad  = 2.0

rad_seas = seasonal_cycle(day_of_year, amp_rad)
rad_year = yearly_offset(year_index, sigma_yrad, seed=12347)
rad_noise = np.random.default_rng(2027).normal(0.0, sigma_drad, size=total_days)

rad_vals = base_rad + rad_seas + rad_year + rad_noise
#rad = np.full((len(time), len(lat), len(lon)), 340.0)
rad = xr.DataArray(rad_vals[:, None, None], dims=('time', 'lat', 'lon'), coords={'time': time, 'lat': lat, 'lon': lon})
create_dataset('rad', 'surface_downwelling_shortwave_flux', 'W m-2', rad).to_netcdf('dummy_rad.nc')
