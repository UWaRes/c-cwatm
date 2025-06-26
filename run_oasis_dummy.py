import pyoasis
import xarray as xr
import numpy as np

import sys
import os

# Add the relative path to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './cwatm/management_modules/')))

from pyoasis_cpl import oasis_specify_partition, oasis_define_grid_simple

# -- read grid data of REMO input file---
filepath = '/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/2000/'
ds = xr.open_dataset(filepath+'e100001n_c140_200001.nc')
lon_2d = ds['lon'].values
lat_2d = ds['lat'].values
nlon_forcing = lon_2d.shape[0]
nlat_forcing = lon_2d.shape[1]

# get landmask from soil water data
landmask_input = ds['WS'][0,:,:].values
landmask_input[landmask_input>0] = 1 


# TODO: put in function?
# --- 1) Initialization ---
### OASIS_INIT_COMP ###
comp = pyoasis.Component("forcing_component")

### get localcomm and define partition
partition, w_unit = oasis_specify_partition(comp,nlon=nlon_forcing,nlat=nlat_forcing)

# grid definition
print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=w_unit)
print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=w_unit)
w_unit.flush()
# function for writing oasis grid information
oasis_define_grid_simple(nlon_forcing,nlat_forcing,lon_2d,lat_2d,landmask_input,partition,'forcing_grid')


# Define the variable to send
#forcing_var = oasis.def_var("FORCING_FIELD", oasis.OUT)
#  DECLARATION OF THE COUPLING FIELDS
################## OASIS_DEF_VAR #################################
#var_id = [None]
# TODO -> "FIELD_RECV_ATM" needs to be declared in namcouple file
#var_id[0] = pyoasis.Var("FIELD_SEND_forcing", self.partition, OASIS.OUT)
#print(f' var_id FRECVATM, {var_id[0]._id}', file=self.w_unit)
#self.w_unit.flush()

### TERMINATION OF DEFINITION PHASE ###
print(' End of initialisation phase', file=w_unit)
w_unit.flush()

### OASIS_ENDDEF ###
comp.enddef()

# -------------------------------------------------------------------

# Load forcing data (e.g., NetCDF)
# maybe read from settings file??
filepath = '/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/2000/'
ds = xr.open_dataset(filepath+'e100001n_c142_200001.nc')
data_array = ds['APRL']  # Replace with your variable name

#print(lat.shape)

# Time loop
#nt = data_array.sizes['time']
for t in range(5): # same number of time loops as C-CWatM
    # Extract data for current timestep
    data_t = data_array.isel(time=t).values

    # send (and get data), see atmos.py example
    # include dummy send and put also in ccwatm

# Finalize
del comp

