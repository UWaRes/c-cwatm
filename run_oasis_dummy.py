import pyoasis
from pyoasis import OASIS
import xarray as xr
import numpy as np

import sys
import os
#import time
import datetime

# Add the relative path to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './cwatm/management_modules/')))

from pyoasis_cpl import pyoasis_cpl
from coupling import MeteoForc2Var
from configuration import parse_configuration


def parse_settings_file(filepath):
    """
    get values from settingsfile
    """
    settings = {}
    current_section = None
    with open(filepath, 'r') as f:
        for line in f:
            # Remove comments and whitespace
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Detect section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                settings[current_section] = {}
                continue

            # Parse key-value pairs
            if '=' in line and current_section:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                settings[current_section][key] = value

    return settings


# -- get data path from settings file ---

binding = parse_settings_file('settings_CCWatM_5min_example.ini')
meteoforc = MeteoForc2Var(binding['COUPLING']['PathForc'],binding['COUPLING']['fmodel_flag'])

# read monthly files
startdate = binding['TIME-RELATED_CONSTANTS']['StepStart']
enddate = binding['TIME-RELATED_CONSTANTS']['StepEnd']
starttime = datetime.datetime.strptime(startdate, '%d/%m/%Y')
endtime = datetime.datetime.strptime(enddate, '%d/%m/%Y')
simulated_days = (endtime-starttime+datetime.timedelta(days=1)) / datetime.timedelta(days=1)

ds = xr.open_dataset('/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/e100001n_c105.nc')
FCAP = np.squeeze(ds['FCAP'].values)
ds = xr.open_dataset('/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/e100001n_c229.nc')
WSMX = np.squeeze(ds['WSMX'].values)
gridfile = ds.copy()

# -- read grid data of REMO input file---
#filepath = '/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/2000/'
#ds = xr.open_dataset(filepath+'e100001n_c140_200001.nc')
lon_2d = gridfile['lon'].values.T
lat_2d = gridfile['lat'].values.T
nlon_forcing = lon_2d.shape[0]
nlat_forcing = lon_2d.shape[1]

# TODO: create coordinate and landmask file to be provided with C-CWatM
# get landmask from soil water data
landmask_input = gridfile['WSMX'][0,:,:].values
landmask_input[landmask_input>0] = 1 

# --- derive grid corners ---
#lat_rot,lon_rot = np.meshgrid(gridfile['rlat'].values,gridfile['rlon'].values)
#rlon_corners,rlat_corners = calc_rotgrid_corners(lon_rot,lat_rot)

#rot_pole_lon = gridfile.rotated_latitude_longitude.grid_north_pole_longitude
#rot_pole_lat = gridfile.rotated_latitude_longitude.grid_north_pole_latitude
#grid_clon = unrot_lon(rlat_corners, rlon_corners, rot_pole_lat, rot_pole_lon)
#grid_clat = unrot_lat(rlat_corners, rlon_corners, rot_pole_lat, rot_pole_lon)


# TODO: put in function?
# --- 1) Initialization ---
### OASIS_INIT_COMP ###
comp = pyoasis.Component("forcing_component")

### get localcomm and define partition
partition, w_unit = pyoasis_cpl.oasis_specify_partition(comp,nlon=nlon_forcing,nlat=nlat_forcing)

# grid definition
print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=w_unit)
print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=w_unit)
w_unit.flush()
# function for writing oasis grid information
pyoasis_cpl.oasis_define_grid(nlon_forcing,nlat_forcing,lon_2d,lat_2d,1-landmask_input.T,partition,'forcing_grid')

#  DECLARATION OF THE COUPLING FIELDS
################## OASIS_DEF_VAR #################################
numcouple = 4 # number of coupling fields
var_id = [None] * numcouple
# TODO -> "FIELD_RECV_ATM" needs to be declared in namcouple file
var_id[0] = pyoasis.Var("FIELD_SEND_runoff", partition, OASIS.OUT)
var_id[1] = pyoasis.Var("FIELD_SEND_gwRecharge", partition, OASIS.OUT)
var_id[2] = pyoasis.Var("FIELD_SEND_EWRef", partition, OASIS.OUT)
var_id[3] = pyoasis.Var("FIELD_SEND_rootzoneSM", partition, OASIS.OUT)
print(f' var_id FRECVATM, {var_id[0]._id}', file=w_unit)
w_unit.flush()

### TERMINATION OF DEFINITION PHASE ###
print(' End of initialisation phase', file=w_unit)
w_unit.flush()

### OASIS_ENDDEF ###
comp.enddef()

# ----- Time loop -----

for t in np.arange(simulated_days): 
    seconds_passed = int(t * 86400.)
    print('dummy',seconds_passed)

    currenttime = starttime + datetime.timedelta(seconds=seconds_passed)
    
    meteoforc.read_forcing(currenttime,'runoff',binding['COUPLING']['RunoffName'])
    meteoforc.read_forcing(currenttime,'sum_gwRecharge',binding['COUPLING']['GWName'])
    meteoforc.read_forcing(currenttime,'EWRef',binding['COUPLING']['OWEName'])
    meteoforc.read_forcing(currenttime,'rootzoneSM',binding['COUPLING']['SMName'])
    # note: this is in percent, 
    # needs to be multiplied with soilWaterStorageCap later
    meteoforc.rootzoneSM = meteoforc.rootzoneSM * FCAP / WSMX

    # -------------- get and put -----------------
    #var_id[0].put(seconds_passed, data_t[::-1,:])
    #var_id[0].put(seconds_passed, data_t.T[:,::-1])
    # C-CWatM needs input in [m]
    var_id[0].put(seconds_passed, meteoforc.runoff.values)
    var_id[1].put(seconds_passed, meteoforc.sum_gwRecharge.values)
    var_id[2].put(seconds_passed, meteoforc.EWRef.values)
    var_id[3].put(seconds_passed, meteoforc.rootzoneSM.values)

# Finalize
del comp

