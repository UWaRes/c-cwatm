# -------------------------------------------------------------------------
# Name:        Coupling
# Purpose:     Classes and functions used for (offline-)coupling with a climate model
#
# Author:      Amelie Schmitt
#
# Created:     09/10/2024
# Copyright:   (c) Amelie Schmitt (2024)
# -------------------------------------------------------------------------

import xarray as xr
import math
import datetime
import os
#from scipy.interpolate import griddata

from cwatm.management_modules.messages import *

class MeteoForc2Var:
    """
    A class to handle meteorological forcing data from climate models.
    Called by the run_read_forcing_modelflag.py model coupled with oasis.

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

        # read dataset
        forc_filename = self.get_filename(ctime,varflag,infile)

        if not(os.path.exists(forc_filename)):
            msg = "File " +forc_filename + " does not exist. \n"
            raise FileNotFoundError(msg)

        try:
            # check if data is already loaded
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
                           'rootzoneSM':'c140',
                           'FCAP':'c105',
                           'WSMX':'c229'}
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



