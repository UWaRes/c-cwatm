import pyoasis
import xarray as xr
#import numpy as np

import sys
import os

# Add the relative path to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './cwatm/management_modules/')))

from pyoasis_cpl import oasis_specify_partition


# TODO: put in function?
# --- 1) Initialization ---
### OASIS_INIT_COMP ###
comp = pyoasis.Component("dummy_component")

### get localcomm and define partition
partition, w_unit = oasis_specify_partition(comp,nlat=433,nlon=433)

# Define the variable to send
#forcing_var = oasis.def_var("FORCING_FIELD", oasis.OUT)

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

