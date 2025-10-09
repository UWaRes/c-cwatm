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
        elif binding['coupl_flag']=='oasis_coupl':  
            # check if oasis has been initialized
            # oasisvar_id should contain at least the 4 forcing variables
            if len(self.var.oasisvar_id) < 4:
                # TODO: raise error
                print('Not all 4 required forcing variables are exchanged via OASIS.')
            else:
                print('Forcing data will be received using OASIS.')
        
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
            self.var.runoff = np.maximum(0., self.var.runoff * self.var.conv_runoff)
    
            # read ground water recharge
            self.var.sum_gwRecharge, MaskMapBoundary = readmeteodata(self.var.GWMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.sum_gwRecharge = np.maximum(0., self.var.sum_gwRecharge * self.var.conv_groundw)
            
            # read rootzone soil moisture
            self.var.rootzoneSM, MaskMapBoundary = readmeteodata(self.var.SMMaps, dateVar['currDate'], addZeros=True, mapsscale = self.var.meteomapsscale, buffering= self.var.buffer)
            self.var.rootzoneSM = np.maximum(0., self.var.rootzoneSM * self.var.conv_soilw)
            
            # read open water evaporation
            self.var.EWRef, MaskMapBoundary = readmeteodata(self.var.OWEMaps, dateVar['currDate'], addZeros=True, mapsscale = True)
            self.var.EWRef = self.var.EWRef * self.var.conv_evap    
             

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

                self.var.meteo = np.zeros([number, 1 + dateVar["intEnd"] - dateVar["intStart"], len(self.var.runoff)])

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

