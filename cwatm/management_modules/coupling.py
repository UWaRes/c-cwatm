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
from scipy.interpolate import griddata

class MeteoForc2Var:
    def __init__(self,clat,clon,ctime,inpath,fmodel_flag):
        """
        fmodel_flag : 'remo' (add others later) model used for meteorological forcing
        """
        self.clat = clat   # 1D array of CWatM latitudes
        self.clon = clon   # 1D array of CWatM longitudes
        self.ctime = ctime   # time at current CWatM step
        self.inpath = inpath   # path where forcing data are stored 
        self.fmodel_flag = fmodel_flag
        # get variable names for specified model
        self.varnames_dict()

    def read_forcing(self,varflag,infile):
        """
        Read climate model output file with forcing data
        varflag : 'runoff', 'sum_gwRecharge', 'rootzoneSM' or 'EWRef'
        infile : given in settings file
        """

        # read dataset
        # TODO: check if exists first??
        ds = xr.open_mfdataset(self.get_filename(infile))  

        # select data for specified date 
        if fmodel_flag == 'remo':
            # convert time format for REMO forcing
            ds = parse_dates(ds)
        
        data_for_date = ds.sel(time=self.ctime)
        
        # get variable
        forc_varname = self.varsdict[varflag]
        datavar = data_for_date[forc_varname]
        setattr(self, varflag, datavar)

        if self.fmodel_flag == 'remo' and varflag=='rootzoneSM':
            # read additional variables required for soil moisture conversion
            setattr(self, 'FCAP', data_for_date['FCAP'])           
            setattr(self, 'WSMX', data_for_date['WSMX'])  

        # get coordinates
        self.forclat = data_for_date[self.varsdict['lat']].values
        self.forclon = data_for_date[self.varsdict['lon']].values
    
    def regridding(self,varflag):
        """
        Convert a 2D forcing variable to a C-CWatM 2D variable.
        """
        # convert units if necessary
        self.convert_units(varflag)

        # prepare coordinates
        if self.fmodel_flag == 'remo':
            # remo has an unstructured grid with 2d coordinates
            forclat_flat = self.forclat.flatten()
            forclon_flat = self.forclon.flatten()
        forcdata_flat = getattr(self,varflag).values.flatten()
        clon_mesh, clat_mesh = np.meshgrid(self.clon, self.clat)
        
        # interpolate
        # TODO: check methods, mask nans?
        # TODO: find more time efficient function
        interpolated_data = griddata((forclat_flat, forclon_flat), forcdata_flat, (clat_mesh, clon_mesh), method='linear')
        setattr(self, varflag, interpolated_data)

    
    # --- functions used by read_forcing ---
    def varnames_dict(self):
        """

        """
        if self.fmodel_flag == 'remo':
            self.varsdict = {'runoff':'RUNOFF' ,
                             'sum_gwRecharge':'DRAIN' , 
                             'EWRef':'EVAPW' ,
                             'rootzoneSM':'WSECH' ,
                             'lat':'lat' ,
                             'lon':'lon' }
        else:
            msg = "Error : " + fmodel_flag + " is not a valid model name. \n"
            raise CWATMError(msg)

    
    def get_filename(self,infile):
        """
        return a string of the path and filenames of the forcing files
        """
        if self.fmodel_flag == 'remo':
            # the date format in the REMO filename is YYYYMM
            rdate = self.ctime.strftime('%Y%m')
            ffname = self.inpath+infile+'_*_'+rdate+'.nc'
            return ffname
        else:
            msg = "Error : " + self.fmodel_flag + " is not a valid model name. \n"
            raise CWATMError(msg)

    # --- functions used by regridding ---
    def convert_units(self,varflag):
        """
        Convert forcing data to match the required units. These are [m/day] for all
        runoff, sum_gwRecharge, EWRef, rootzoneSM
        """
        if self.fmodel_flag == 'remo':
            # runoff: mm to m/day - 0.001
            # sum_gwRecharge: mm/day to m/day - 0.001
            # EWRef: mm/day to m/day - 0.001
            # rootzoneSM: WSECH*FCAP/WSMX; m/m/day to (m/day)
            if varflag in ['runoff', 'sum_gwRecharge', 'EWRef']:
                varconv = getattr(self,varflag) * 0.001
                setattr(self, varflag, varconv)
            elif varflag=='rootzoneSM':
                # note: this is in percent, 
                # needs to be multiplied with soilWaterStorageCap later
                self.rootzoneSM = self.rootzoneSM * self.FCAP / self.WSMX
            else:
                pass
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



