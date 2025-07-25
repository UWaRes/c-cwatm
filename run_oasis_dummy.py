# -------------------------------------------------------------------------
# Name:        run_oasis_dummy
# Purpose:     Read REMO forcing data for every time step and send them
#              to C-CWatM using the OASIS3-MCT coupler.
#              Can be adjusted for different forcing datasets.
#
# Author:      Amelie Schmitt
#
# Created:     25/07/2025
# Copyright:   (c) Amelie Schmitt 2025 
# -------------------------------------------------------------------------

import pyoasis
from pyoasis import OASIS
import xarray as xr
import numpy as np
import sys
import os
import time
import datetime

# Add the relative path to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './cwatm/management_modules/')))

from pyoasis_cpl import pyoasis_cpl
from coupling import MeteoForc2Var, parse_settings_file, read_modelrun_info, init_oasis_forcing
from cwatm.management_modules.grid_tools import grid_tools

def read_remo_info(binding,starttime):
    """
    Read grid, landmask, and other auxiliary information for REMO.

    TODO: 
      - include optional grid corners and area
      - check if files exist

    Parameters
    ----------
    binding (dict) : Dictionary from C-CWatM settingsfile with keys structured as [headingname][varname].
    starttime (datetime.datetime) : Date for which a soil water file (_c140_) exists.
      
    Returns
    -------
    nlon_forcing (int) : Number of columns in grid (longitudes).
    nlat_forcing (int) : Number of rows in grid (latitudes).
    lon_2d (ndarray) : 2D array of grid cell center longitudes.
    lat_2d (ndarray) : 2D array of grid cell center latitudes.
    landmask (ndarray) : 2D array representing the land-sea mask; with 0 for land and 1 for ocean.
    soilwat_factor (ndarray) : Factor to convert soil water content in m to percent, 
                               corresponds to soil depth.

    """
    # --- read auxiliary information ---
    # factor for soil water content conversion 
    filepath = binding['COUPLING']['PathForc']
    file_prefix = binding['COUPLING']['RunoffName']
    ds = xr.open_dataset(filepath+file_prefix+'_c105.nc')
    FCAP = np.squeeze(ds['FCAP'].values)
    ds = xr.open_dataset(filepath+file_prefix+'_c229.nc')
    WSMX = np.squeeze(ds['WSMX'].values)
    soilwat_factor = FCAP/WSMX
   
    # --- read grid data of REMO input file ---
    rdate = starttime.strftime('%Y%m')
    ryear = starttime.strftime('%Y')  
    gridfile = xr.open_dataset(filepath+ryear+'/'+file_prefix+'_c140_'+rdate+'.nc')
    lon_2d = gridfile['lon'].values
    lat_2d = gridfile['lat'].values
    nlon_forcing = lon_2d.shape[1]
    nlat_forcing = lon_2d.shape[0]
    
    # --- get landmask from soil water data ---
    landmask_input = gridfile['WS'][0,:,:].values
    landmask_input[landmask_input>0] = 1 

    # --- derive grid corners and grid cell area ---
    # optional: only required for certain regridding methods
    if 0:
        rot_pole_lon = gridfile.rotated_latitude_longitude.grid_north_pole_longitude
        rot_pole_lat = gridfile.rotated_latitude_longitude.grid_north_pole_latitude
        # TODO: check orientation
        lat_rot,lon_rot = np.meshgrid(gridfile['rlat'].values,gridfile['rlon'].values)
        grid_clon_rot, grid_clat_rot = grid_tools.compute_grid_corners(lon_rot,lat_rot)
        grid_clon, grid_clat = grid_tools.unrot_coordinates(grid_clon_rot, grid_clat_rot, rot_pole_lon, rot_pole_lat)
        # calculation takes about 1s
        cell_areas = grid_tools.compute_grid_cell_areas(grid_clon, grid_clat)

    return nlon_forcing,nlat_forcing,lon_2d,lat_2d,1-landmask_input, soilwat_factor


############################
# ----- Initialization -----
############################

# --- get info from settings file ---
settingsfile = sys.argv[1]
binding = parse_settings_file(settingsfile)
starttime, simulated_days = read_modelrun_info(binding)

# --- initialize reading of forcing data ---
meteoforc = MeteoForc2Var(binding['COUPLING']['PathForc'],binding['COUPLING']['fmodel_flag'])

# --- read remo grid information and other info ---
nlon_forcing,nlat_forcing,lon_2d,lat_2d,landmask_input,soilwat_factor = read_remo_info(binding['COUPLING']['PathForc'],starttime)

# --- initialize OASIS coupler ---
# transpose REMO grid to match C-CWatM orientation
comp,w_unit,var_id = init_oasis_forcing(nlon_forcing,nlat_forcing,lon_2d.T,lat_2d.T,landmask_input.T)


#######################
# ----- Time loop -----
#######################

for t in np.arange(simulated_days): 
    seconds_passed = int(t * 86400.)
    currenttime = starttime + datetime.timedelta(seconds=seconds_passed)

    # --- read REMO forcing ---
    meteoforc.read_forcing(currenttime,'runoff',binding['COUPLING']['RunoffName'])
    meteoforc.read_forcing(currenttime,'sum_gwRecharge',binding['COUPLING']['GWName'])
    meteoforc.read_forcing(currenttime,'EWRef',binding['COUPLING']['OWEName'])
    meteoforc.read_forcing(currenttime,'rootzoneSM',binding['COUPLING']['SMName'])
    # transform soil water content from [m] to [%]
    meteoforc.rootzoneSM = meteoforc.rootzoneSM * soilwat_factor

    # -------------- get and put -----------------
    # C-CWatM needs input in [m]
    # transpose REMO variables to match C-CWatM orientation
    var_id[0].put(seconds_passed, meteoforc.runoff.values.T)
    var_id[1].put(seconds_passed, meteoforc.sum_gwRecharge.values.T)
    var_id[2].put(seconds_passed, meteoforc.EWRef.values.T)
    var_id[3].put(seconds_passed, meteoforc.rootzoneSM.values.T)

# --- finalize OASIS ---
del comp

