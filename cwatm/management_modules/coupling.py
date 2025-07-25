# -------------------------------------------------------------------------
# Name:        Coupling
# Purpose:     Classes and functions used for reading forcing data, 
#              used for offline-coupling with OASIS3-MCT
#
# Author:      Amelie Schmitt
#
# Created:     09/10/2024
# Adjusted for OASIS: 23/07/2025
# Copyright:   (c) Amelie Schmitt (2025)
# -------------------------------------------------------------------------

import xarray as xr
import math
import datetime
import os
import pyoasis
from pyoasis import OASIS
import numpy as np

from cwatm.management_modules.pyoasis_cpl import pyoasis_cpl
from cwatm.management_modules.messages import *

class MeteoForc2Var:
    """
    A class to handle meteorological forcing data from climate models.
    Called by the run_read_forcing_modelflag.py model coupled via oasis.

    Attributes
    ----------
        inpath (str): Path to the directory containing forcing data files.
        fmodel_flag (str): Identifier for the climate model used (e.g., 'remo').
        varsdict (dict): Mapping of variable flags to variable names in the dataset.
        Optional:
        varsid (dict): Mapping of variable flags to identifier strings used in filenames.
    """
    def __init__(self,inpath,fmodel_flag):
        """
        Parameters
        ----------
        inpath (str): Path to the directory containing forcing data files.
        fmodel_flag (str): Model identifier.
            List of valid model identifiers:
                - 'remo' : REMO model output; monthly files of daily values; yearly folders 
        """

        self.inpath = inpath   # path where forcing data are stored 
        self.fmodel_flag = fmodel_flag
        # get variable names (and filenames) for specified model
        self.varnames_dict()

    def read_forcing(self,ctime,varflag,infile):
        """
        Read and extract meteorological forcing data for a specific date and variable.

        Parameters
        ----------
        ctime (datetime): The date for which data is to be extracted.
        varflag (str): C-CWatM variable name.
        infile (str): Base name of the input file as specified in the settings.

        Raises
        ------
        FileNotFoundError
            - If the corresponding NetCDF file does not exist.

        as/copilot    
        """
        # check if dataset exists
        forc_filename = self.get_filename(ctime,varflag,infile)
        if not(os.path.exists(forc_filename)):
            msg = "File " +forc_filename + " does not exist. \n"
            raise FileNotFoundError(msg)

        try:
            # check if dataset is already loaded
            ds = getattr(self, varflag+'_infile')
        except (KeyError, AttributeError):
            # otherwise: load file
            ds = xr.open_dataset(forc_filename)  
            if self.fmodel_flag == 'remo':
                # convert time format for REMO forcing
                ds = parse_dates(ds)
            setattr(self, varflag+'_infile', ds)

        # select data for specified date 
        data_for_date = ds.sel(time=ctime)
        # get variable
        forc_varname = self.varsdict[varflag]
        datavar = data_for_date[forc_varname]
        setattr(self, varflag, datavar)

    
    # --- functions used by read_forcing ---
    def varnames_dict(self):
        """
        Define the mapping of variable flags to dataset variable names and filename identifiers.

        Sets
        ----
        self.varsdict (dict): Maps variable flags to variable names in the dataset.
        Optional:
        self.varsid (dict): Maps variable flags to identifier strings used in filenames.

        Raises
        ------
        InvalidModelflagError
            - If the model flag is not recognized.

        as/copilot
        """
        
        if self.fmodel_flag == 'remo':
            self.varsdict = {'runoff':'RUNOFF' ,
                             'sum_gwRecharge':'DRAIN' , 
                             'EWRef':'EVAPW' ,
                             'rootzoneSM':'WS'}
            self.varsid = {'runoff':'c160' ,
                           'sum_gwRecharge':'c053' , 
                           'EWRef':'c064' ,
                           'rootzoneSM':'c140'}
        else:
            raise InvalidModelflagError('fmodel_flag',self.fmodel_flag)

    def get_filename(self,ctime,varflag,infile):
        """
        Construct the full path to the NetCDF file for a given variable and date.

        Parameters
        ----------
        ctime (datetime): The date for which the filename is constructed.
        varflag (str): Variable name.
        infile (str): Base name of the input file.

        Returns
        -------
        ffname (str) : Full path to the NetCDF file.

        Raises
        ------
        InvalidModelflagError
            - If the model flag is not recognized.

        as/copilot    
        """
        if self.fmodel_flag == 'remo':
            # the date format in the REMO filename is YYYYMM
            rdate = ctime.strftime('%Y%m')
            ryear = ctime.strftime('%Y')
            varnumber = self.varsid[varflag]
            ffname = self.inpath+'/'+ryear+'/'+infile+'_'+varnumber+'_'+rdate+'.nc'
            return ffname
        else:
            raise InvalidModelflagError('fmodel_flag',self.fmodel_flag)

           
# ----- functions used for 'remo' -----
def parse_dates(ds):
    """
    Converts the time axis of a REMO xarray dataset from absolute time values to datetime.
    Functionality based on the pyremo package.
    """
    ds = ds.copy()
    ds["time"] = [num2date(date) for date in ds.time]
    return ds

def num2date(num):
    """
    Convert a numeric absolute date value to a datetime object.
    Functionality based on the pyremo package.
    """
    frac, whole = math.modf(num)
    date_str = str(int(whole))
    date = datetime.datetime.strptime(date_str[0:8], "%Y%m%d")
    datetime0 = date + datetime.timedelta(seconds=datetime.timedelta(days=1).total_seconds() * frac)
    return datetime0

# ----- function used to read forcing -----
def parse_settings_file(filepath):
    """
    Get values from C-CWatM settingsfile.

    Parameters
    ----------
    filepath (str) : path and name of C-CWatM settingsfile

    Returns
    -------
    settings (dict) : dictionary with keys structured as [headingname][varname]
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

def read_modelrun_info(binding):
    """
    Extract the model run start time and number of simulated days from the C-CWatM settingsfile dictionary.

    Parameters
    ----------
    binding (dict) : A nested dictionary from C-CWatM settingsfile with keys structured as [headingname][varname]; 
                     with keys 'StepStart' and 'StepEnd' in the format '%d/%m/%Y'.

    Returns
    -------
    starttime (datetime.datetime) : The start date of the model run.
    simulated_days (float) : The total number of simulated days, including both start and end dates.

    as/copilot
    """
    startdate = binding['TIME-RELATED_CONSTANTS']['StepStart']
    enddate = binding['TIME-RELATED_CONSTANTS']['StepEnd']
    starttime = datetime.datetime.strptime(startdate, '%d/%m/%Y')
    endtime = datetime.datetime.strptime(enddate, '%d/%m/%Y')
    simulated_days = (endtime-starttime+datetime.timedelta(days=1)) / datetime.timedelta(days=1)

    return starttime, simulated_days

    
def init_oasis_forcing(nlon_forcing,nlat_forcing,lon_2d,lat_2d,landmask_input):
    """
    Initializes the OASIS coupling interface for the forcing component.

    This function sets up the OASIS component, defines the grid and partitioning,
    and declares the coupling fields to be exchanged. 

    Parameters
    ----------
    nlon_forcing (int) : Number of longitudinal grid points in the forcing data.
    nlat_forcing (int) : Number of latitudinal grid points in the forcing data.
    lon_2d (ndarray) : 2D array of longitudes for the forcing grid.
    lat_2d (ndarray) : 2D array of latitudes for the forcing grid.
    landmask_input (ndarray) : 2D array representing the land-sea mask; with 0 for land and 1 for ocean.

    Returns
    -------
    comp (pyoasis.Component) : The initialized OASIS component for the forcing model.
    w_unit (file-like) : The file-like object used for logging OASIS grid information.
    var_id (list) : List of OASIS variable objects representing the declared coupling fields.

    as/copilot
    """

    # --- OASIS_INIT_COMP ---
    comp = pyoasis.Component("forcing_component")
    
    # get localcomm and define partition
    partition, w_unit = pyoasis_cpl.oasis_specify_partition(comp,nlon=nlon_forcing,nlat=nlat_forcing)
    
    # grid definition 
    print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=w_unit)
    print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=w_unit)
    w_unit.flush()
    # write OASIS grid information
    pyoasis_cpl.oasis_define_grid(nlon_forcing,nlat_forcing,lon_2d,lat_2d,landmask_input,partition,'forcing_grid')
    
    # --- DECLARATION OF THE COUPLING FIELDS ---
    numcouple = 4 # number of coupling fields
    var_id = [None] * numcouple
    var_id[0] = pyoasis.Var("FIELD_SEND_runoff", partition, OASIS.OUT)
    var_id[1] = pyoasis.Var("FIELD_SEND_gwRecharge", partition, OASIS.OUT)
    var_id[2] = pyoasis.Var("FIELD_SEND_EWRef", partition, OASIS.OUT)
    var_id[3] = pyoasis.Var("FIELD_SEND_rootzoneSM", partition, OASIS.OUT)
    print(f' var_id FRECVATM, {var_id[0]._id}', file=w_unit)
    w_unit.flush()
    
    # --- TERMINATION OF DEFINITION PHASE ---
    print(' End of initialisation phase', file=w_unit)
    w_unit.flush()
    
    # --- OASIS_ENDDEF ---
    comp.enddef()

    return comp,w_unit,var_id