# -------------------------------------------------------------------------
# Name:        READ Landsurface input maps
# Purpose:
#
# Author:      Peter Greve, Peter Burek
#
# Created:     23/08/2023
# Copyright:   (c) Peter Burek 2016, (c) Peter Greve 2023
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
import scipy.ndimage
from scipy.interpolate import RegularGridInterpolator

class readmeteo(object):
    """
    READ METEOROLOGICAL DATA

    reads all meteorological data from netcdf4 files

    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    DtDay                                  seconds in a timestep (default=86400)                                   s    
    con_precipitation                      conversion factor for precipitation                                     --   
    con_e                                  conversion factor for evaporation                                       --   
    ETRef                                  potential evapotranspiration rate from reference crop                   m    
    Precipitation                          Precipitation (input for the model)                                     m    
    only_radiation                                                                                                  --
    TMin                                   minimum air temperature                                                 K    
    TMax                                   maximum air temperature                                                 K    
    Tavg                                   Input, average air Temperature                                          K    
    Rsds                                   short wave downward surface radiation fluxes                            W/m2 
    EAct                                                                                                           --   
    Psurf                                  Instantaneous surface pressure                                          Pa   
    Qair                                   specific humidity                                                       kg/kg
    Rsdl                                   long wave downward surface radiation fluxes                             W/m2 
    Wind                                   wind speed                                                              m/s  
    EWRef                                  potential evaporation rate from water surface                           m    
    meteomapsscale                         if meteo maps have the same extend as the other spatial static maps ->  --   
    meteodown                              if meteo maps should be downscaled                                      --   
    InterpolationMethod                                                                                            --   
    buffer                                                                                                         --   
    glaciermeltMaps                                                                                                --   
    glacierrainMaps                                                                                                --   
    wc2_tavg                               High resolution WorldClim map for average temperature                   K    
    wc4_tavg                               upscaled to low resolution WorldClim map for average temperature        K    
    wc2_tmin                               High resolution WorldClim map for min temperature                       K    
    wc4_tmin                               upscaled to low resolution WorldClim map for min temperature            K    
    wc2_tmax                               High resolution WorldClim map for max temperature                       K    
    wc4_tmax                               upscaled to low resolution WorldClim map for max temperature            K    
    wc2_prec                               High resolution WorldClim map for precipitation                         m    
    wc4_prec                               upscaled to low resolution WorldClim map for precipitation              m    
    xcoarse_prec                                                                                                   --   
    ycoarse_prec                                                                                                   --   
    xfine_prec                                                                                                     --   
    yfine_prec                                                                                                     --   
    meshlist_prec                                                                                                  --   
    xcoarse_tavg                                                                                                   --   
    ycoarse_tavg                                                                                                   --   
    xfine_tavg                                                                                                     --   
    yfine_tavg                                                                                                     --   
    meshlist_tavg                                                                                                  --   
    meteo                                                                                                          --   
    prec                                   precipitation in m                                                      m    
    temp                                   average temperature in Celsius deg                                      °C   
    WtoMJ                                  Conversion factor from [W] to [MJ] for radiation: 86400 * 1E-6          --   
    includeGlaciers                                                                                                --   
    includeOnlyGlaciersMelt                                                                                        --   
    GlacierMelt                                                                                                    --   
    GlacierRain                                                                                                    --   
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.model = model
        self.var = model.var

    def initial(self):
        """
        Initial part of meteo

        read multiple file of input
        """

        # fit land surface forcing data to size and resolution of mask map
        #-------------------------------------------------------------------

        name = cbinding('RunoffMaps')
        nameall = glob.glob(os.path.normpath(name))
        if not nameall:
            msg = "Error 215: In readmeteo, cannot find runoff maps " 
            raise CWATMFileError(name, msg, sname='RunoffMaps') 
        namemeteo = nameall[0]
        latmeteo, lonmeteo, cell, invcellmeteo, rows, cols = readCoordNetCDF(namemeteo)

        nameldd = cbinding('Ldd')
        latldd, lonldd, cell, invcellldd, rows, cols = readCoord(nameldd)
        maskmapAttr['reso_mask_meteo'] = round(invcellldd / invcellmeteo)

        # if meteo maps have the same extend as the other spatial static maps -> meteomapsscale = True
        self.var.meteomapsscale = True
        if invcellmeteo != invcellldd:
            if (not(Flags['quiet'])) and (not(Flags['veryquiet'])) and (not(Flags['check'])):
                msg = "Resolution of meteo forcing is " + str(maskmapAttr['reso_mask_meteo']) + " times higher than base maps."
                print(msg)
            self.var.meteomapsscale = False

        cutmap[0], cutmap[1], cutmap[2], cutmap[3] = mapattrNetCDF(nameldd)
        for i in range(4): cutmapFine[i] = cutmap[i]

        
        self.var.meteodown = False
        self.var.buffer = False

        # in case other mapsets are used e.g. Cordex RCM meteo data
        if (latldd != latmeteo) or (lonldd != lonmeteo):
            cutmapFine[0], cutmapFine[1], cutmapFine[2], cutmapFine[3], cutmapVfine[0], cutmapVfine[1], cutmapVfine[2], cutmapVfine[3] = mapattrNetCDFMeteo(namemeteo)

        # -------------------------------------------------------------------
        self.var.includeGlaciers = False
        if 'includeGlaciers' in option:
            self.var.includeGlaciers = checkOption('includeGlaciers')
            self.var.includeOnlyGlaciersMelt = False
            if 'includeOnlyGlaciersMelt' in binding:
                self.var.includeOnlyGlaciersMelt = returnBool('includeOnlyGlaciersMelt')

        if binding['coupl_flag']=='no_coupl':
            # read all forcing data at once
            self.var.QMaps = 'RunoffMaps' 
            self.var.GWMaps = 'GWMaps'
            self.var.OWEMaps = 'OWEMaps'
            self.var.SMMaps = 'SMMaps'
            
            # Read landsurface maps
            meteomaps = [self.var.QMaps, self.var.GWMaps, self.var.SMMaps, self.var.OWEMaps] #Peter Greve Test
            multinetdf(meteomaps)
        
        if self.var.includeGlaciers:
            self.var.glaciermeltMaps = 'MeltGlacierMaps'
            if not self.var.includeOnlyGlaciersMelt:
                self.var.glacierrainMaps = 'PrecGlacierMaps'

        # use radiation term in snow melt
        self.var.snowmelt_radiation = False
        self.var.only_radiation = False
        

    def dynamic(self):
        """
        Dynamic part of the readmeteo module

        Read landsurface input maps from netcdf files


        """
        if Flags['warm']:
            # if warmstart use stored meteo variables
            no = dateVar['curr']-1
            self.var.runoff = self.var.meteo[0,no]
            self.var.sum_gwRecharge = self.var.meteo[1,no]
            self.var.rootzoneSM = self.var.meteo[2,no]
            self.var.EWRef = self.var.meteo[3,no]
            j = 3
            if self.var.includeGlaciers:
                self.var.GlacierMelt = self.var.meteo[j+1, no]
                if not self.var.includeOnlyGlaciersMelt:
                    self.var.GlacierRain = self.var.meteo[j+2, no]
            return


        if binding['coupl_flag']=='no_coupl':
            # extract forcing data from Maps (read in initial)
            # read runoff
            self.var.runoff, MaskMapBoundary = readmeteodata(self.var.QMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.runoff = np.maximum(0., self.var.runoff)
    
            # read ground water recharge
            self.var.sum_gwRecharge, MaskMapBoundary = readmeteodata(self.var.GWMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.sum_gwRecharge = np.maximum(0., self.var.sum_gwRecharge)
            
            # read rootzone soil moisture
            self.var.rootzoneSM, MaskMapBoundary = readmeteodata(self.var.SMMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.rootzoneSM = np.maximum(0., self.var.rootzoneSM)
            
            # read open water evaporation
            self.var.EWRef, MaskMapBoundary = readmeteodata(self.var.OWEMaps, dateVar['currDate'], addZeros=True, mapsscale = True)
            self.var.EWRef = self.var.EWRef * self.var.DtDay * self.var.con_e
            
        elif binding['coupl_flag']=='offline_coupl':   
            # for each time step: read forcing file and convert to C-CWatM grid
            # TODO: replace clat and clon with actual C-CWatM variables
            meteoforc = MeteoForc2Var(clat,clon,dateVar['currDate'],binding['PathForc'],binding['fmodel_flag'])
            
            # read runoff
            meteoforc.read_forcing('runoff',binding['RunoffName'])
            meteoforc.regridding('runoff')
            self.var.runoff = np.maximum(0., meteoforc.runoff)
            
            # read ground water recharge
            meteoforc.read_forcing('sum_gwRecharge',binding['GWName'])
            meteoforc.regridding('sum_gwRecharge')
            self.var.sum_gwRecharge = np.maximum(0., meteoforc.sum_gwRecharge)

            # read open water evaporation
            meteoforc.read_forcing('EWRef',binding['OWEName'])
            meteoforc.regridding('EWRef')
            self.var.EWRef = meteoforc.EWRef * self.var.DtDay * self.var.con_e # TODO: check conversion

            # read rootzone soil moisture
            meteoforc.read_forcing('rootzoneSM',binding['SMName'])
            meteoforc.regridding('rootzoneSM')
            self.var.rootzoneSM = np.maximum(0., meteoforc.rootzoneSM)            

            #if returnBool('soilwater_as_fract'):
            #    # convert soil moisture in percent to soil water content in m 
            #    soilWaterStorageCap = self.var.ws1[3] + self.var.ws2[3]
            #    self.var.rootzoneSM = self.var.rootzoneSM * soilWaterStorageCap

                

        if Flags['calib']:
            # if first clibration run, store all meteo data in a variable
            if dateVar['curr'] == 1:
                number = 4
                if self.var.snowmelt_radiation:
                    number = number + 2
                if self.var.includeGlaciers:
                    number = number + 1
                    if not self.var.includeOnlyGlaciersMelt:
                        number = number + 1

                self.var.meteo = np.zeros([number, 1 + dateVar["intEnd"] - dateVar["intStart"], len(self.var.Precipitation)])

            no = dateVar['curr'] -1
            self.var.meteo[0,no] = self.var.runoff
            self.var.meteo[1,no] = self.var.sum_gwRecharge
            self.var.meteo[2,no] = self.var.rootzoneSM
            self.var.meteo[3,no] = self.var.EWRef
            j =3

            if self.var.includeGlaciers:
                self.var.meteo[j+1, no] = self.var.GlacierMelt
                if not self.var.includeOnlyGlaciersMelt:
                    self.var.meteo[j+2, no] = self.var.GlacierRain
            ii =1

