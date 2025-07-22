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
    def __init__(self,inpath,fmodel_flag):
        """
        fmodel_flag : 'remo' (add others later) model/data used for meteorological forcing
        """
        self.inpath = inpath   # path where forcing data are stored 
        self.fmodel_flag = fmodel_flag
        # get variable names (and filenames) for specified model
        self.varnames_dict()

    def read_forcing(self,ctime,varflag,infile):
        """
        Read climate model output file with forcing data
        varflag : 'runoff', 'sum_gwRecharge', 'rootzoneSM' or 'EWRef'
        infile : given in settings file

        raises cwatmerror if file does not exist
        
        """

        # read dataset
        forc_filename = self.get_filename(ctime,varflag,infile)

        if not(os.path.exists(forc_filename)):
            msg = "Error : File " +forc_filename + " does not exist. \n"
            raise CWATMError(msg)
            # TODO: both coupled models need to exit!!

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
        self.varsdict - variable names in forcing files
        self.varsid - optional, str in filename
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
            msg = "Error : " +self.fmodel_flag + " is not a valid model name. \n"
            raise CWATMError(msg)

    
    def get_filename(self,ctime,varflag,infile):
        """
        return a string of the path and filenames of the forcing files
        """
        if self.fmodel_flag == 'remo':
            # the date format in the REMO filename is YYYYMM
            rdate = ctime.strftime('%Y%m')
            ryear = ctime.strftime('%Y')
            varnumber = self.varsid[varflag]
            ffname = self.inpath+'/'+ryear+'/'+infile+'_'+varnumber+'_'+rdate+'.nc'
            return ffname
        else:
            msg = "Error : " + self.fmodel_flag + " is not a valid model name. \n"
            raise CWATMError(msg)

           

# ----- functions used for 'remo' -----
def parse_dates(ds):
    """
    Updates the time axis of a REMO dataset containing an absolute time axis.
    Based on pyremo
    """
    ds = ds.copy()
    ds["time"] = [num2date(date) for date in ds.time]
    return ds

def num2date(num):
    """
    Convert a numeric absolute date value to a datetime object.
    Based on pyremo
    """
    frac, whole = math.modf(num)
    date_str = str(int(whole))
    date = datetime.datetime.strptime(date_str[0:8], "%Y%m%d")
    datetime0 = date + datetime.timedelta(seconds=datetime.timedelta(days=1).total_seconds() * frac)
    return datetime0



