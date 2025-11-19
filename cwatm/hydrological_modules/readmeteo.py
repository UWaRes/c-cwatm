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
    conv_evap                              conversion factor for evaporation                                       --   
    conv_runoff                            conversion factor for runoff                                            --
    conv_groundw                           conversion factor for groundwater recharge                              --
    conv_soilw                             conversion factor for soilwater content                                 --   
    only_radiation                                                                                                 --  
    EWRef                                  potential evaporation rate from water surface                           m    
    meteomapsscale                         if meteo maps have the same extend as the other spatial static maps ->  --   
    meteodown                              if meteo maps should be downscaled                                      --   
    InterpolationMethod                                                                                            --   
    buffer                                                                                                         --   
    glaciermeltMaps                                                                                                --   
    glacierrainMaps                                                                                                --       
    meteo                                                                                                          --       
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

        if binding['coupl_flag']=='no_coupl':
            # read all forcing data at once
            self.var.QMaps = 'RunoffMaps' 
            self.var.SMMaps = 'SMMaps'
            if returnBool('use_GWrecharge'):
                self.var.GWMaps = 'GWMaps'
            if binding['OWE_meth']=='OWE':  
                self.var.OWEMaps = 'OWEMaps'
            if binding['OWE_meth']=='T':  
                self.var.OWEMaps = 'TMaps'
            if binding['OWE_meth']=='Rnet':  
                self.var.OWEnMaps = 'RnMaps'
            
            # Read landsurface maps
            if returnBool('use_GWrecharge'):
                meteomaps = [self.var.QMaps, self.var.GWMaps, self.var.SMMaps, self.var.OWEMaps]
            else:
                meteomaps = [self.var.QMaps, self.var.SMMaps, self.var.OWEMaps]
                
            multinetdf(meteomaps)
            
        elif binding['coupl_flag']=='oasis_coupl':  
            # check if oasis has been initialized
            # oasisvar_id should contain at least the 4 forcing variables
            if len(self.var.oasisvar_id) < 4:
                # TODO: raise error
                print('Not all 4 required forcing variables are exchanged via OASIS.')
            else:
                print('Forcing data will be received using OASIS.')

        # use radiation term in snow melt
        self.var.snowmelt_radiation = False
        self.var.only_radiation = False

        # Load quantile mapping weights if bias correction is enabled
        if checkOption('bias_correct_runoff'):            
            csv_path = binding.get('bias_correction_file', None)
            if csv_path and os.path.exists(csv_path):
                try:
                    data = np.genfromtxt(csv_path, delimiter=',', invalid_raise=False)
                    data = data[~np.isnan(data).any(axis=1)]  # remove rows with NaNs
                    self.var.quantile_bounds = data[:, 0]
                    self.var.correction_factors = data[:, 1]
                except Exception as e:
                    print("Error reading bias correction file:", e)
            else:
                print("Bias correction file not found or not specified.")


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
            return


        if binding['coupl_flag']=='no_coupl':
            # extract forcing data from Maps (read in initial)
            # read runoff
            self.var.runoff, MaskMapBoundary = readmeteodata(self.var.QMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.runoff = np.maximum(0., self.var.runoff * self.var.conv_runoff)
            
            # read rootzone soil moisture
            self.var.rootzoneSM, MaskMapBoundary = readmeteodata(self.var.SMMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.rootzoneSM = np.maximum(0., self.var.rootzoneSM * self.var.conv_soilw)

            if returnBool('use_GWrecharge'):
                # read ground water recharge
                self.var.sum_gwRecharge, MaskMapBoundary = readmeteodata(self.var.GWMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
                self.var.sum_gwRecharge = np.maximum(0., self.var.sum_gwRecharge * self.var.conv_groundw)
            else:
                self.var.GWscale = loadmap('GWscale_factor')
                self.var.sum_gwRecharge = self.var.runoff*self.var.GWscale

            if binding['OWE_meth']=='OWE':
                # read open water evaporation
                self.var.EWRef, MaskMapBoundary = readmeteodata(self.var.OWEMaps, dateVar['currDate'], addZeros=True, mapsscale = True)
                self.var.EWRef = self.var.EWRef * self.var.conv_evap

            if binding['OWE_meth']=='T':
                # read temperature
                self.var.Toffset = loadmap('Temp_offset')
                self.var.Tin, MaskMapBoundary = readmeteodata(self.var.OWEMaps, dateVar['currDate'], addZeros=True, mapsscale = True)
                self.var.EWRef = 13.97 * 0.0495 * np.exp(0.062 * (self.var.Tin+self.var.Toffset) )/1000 # in m/day 
                # based on Hamon (1963) and pyet Python package

            if binding['OWE_meth']=='Rnet':
                # read net radiation
                self.var.Rnetfact = loadmap('Rnet_conversion')
                self.var.Rnet, MaskMapBoundary = readmeteodata(self.var.OWEMaps, dateVar['currDate'], addZeros=True, mapsscale = True)
                self.var.EWRef = ((0.8 * self.var.Rnet*self.var.Rnetfact)/28.94)/1000 # based on Milly and Dunne, 2016

        

        # Apply bias correction if enabled
        if checkOption('bias_correct_runoff'):
            self.var.runoff = apply_bias_correction(
                self.var.runoff,
                self.var.quantile_bounds,
                self.var.correction_factors
            )                

        if Flags['calib']:
            # if first clibration run, store all meteo data in a variable
            if dateVar['curr'] == 1:
                number = 4
                if self.var.snowmelt_radiation:
                    number = number + 2

                self.var.meteo = np.zeros([number, 1 + dateVar["intEnd"] - dateVar["intStart"], len(self.var.runoff)])

            no = dateVar['curr'] -1
            self.var.meteo[0,no] = self.var.runoff
            self.var.meteo[1,no] = self.var.sum_gwRecharge
            self.var.meteo[2,no] = self.var.rootzoneSM
            self.var.meteo[3,no] = self.var.EWRef


