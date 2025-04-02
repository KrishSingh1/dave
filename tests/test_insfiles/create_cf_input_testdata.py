# This code generates the netCDF input files for the cfinput integration test
import xarray as xr
import numpy as np

total_days = 31*365 # need to have more than 30 years for the spinup to work
lon = [9.5]
lat = [49.5]
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


temp = np.full((len(time), len(lat), len(lon)), 285.0)
temp = xr.DataArray(temp, dims=('time', 'lat', 'lon'), coords={'time': time, 'lat': lat, 'lon': lon})
create_dataset('temp', 'air_temperature', 'K', temp).to_netcdf('dummy_temp.nc')

pr = np.full((len(time), len(lat), len(lon)), 1.0)
pr = xr.DataArray(pr, dims=('time', 'lat', 'lon'), coords={'time': time, 'lat': lat, 'lon': lon})
create_dataset('pr', 'precipitation_amount', 'kg m-2', pr).to_netcdf('dummy_pr.nc')

rad = np.full((len(time), len(lat), len(lon)), 340.0)
rad = xr.DataArray(rad, dims=('time', 'lat', 'lon'), coords={'time': time, 'lat': lat, 'lon': lon})
create_dataset('rad', 'surface_downwelling_shortwave_flux', 'W m-2', rad).to_netcdf('dummy_rad.nc')
