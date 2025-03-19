# -------------------------------------------------------------------------
# Name:        Waterdemand modules
# Purpose:
#
# Author:      PB, YS, MS, JdB, PG
#
# Created:     15/07/2016
# Copyright:   (c) PB 2016, PG 2025
# -------------------------------------------------------------------------

from cwatm.management_modules import globals
from cwatm.management_modules.data_handling import *
import numpy as np

# from cwatm.management_modules.data_handling import *  # luca for testing
# import matplotlib.pyplot as plt


# def decompress(map, nanvalue=None):
#    """
#    Decompressing CWatM maps from 1D to 2D with missing values
#
#    :param map: compressed map
#    :return: decompressed 2D map
#    """
#
#    dmap = maskinfo['maskall'].copy()
#    dmap[~maskinfo['maskflat']] = map[:]
#    if nanvalue is not None:
#        dmap.data[np.isnan(dmap.data)] = nanvalue
#
#    return dmap.data

class waterdemand_irrigation:
    """
    WATERDEMAND

    calculating water demand - irrigation
    Agricultural water demand based on water need by plants

    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    load_initial                           Settings initLoad holds initial conditions for variables                input
    cropKC                                 crop coefficient for each of the 4 different land cover types (forest,  --   
    topwater                               quantity of water above the soil (flooding)                             m    
    efficiencyPaddy                        Input, irrPaddy_efficiency, paddy irrigation efficiency, the amount of  frac 
    efficiencyNonpaddy                     Input, irrNonPaddy_efficiency, non-paddy irrigation efficiency, the am  frac 
    returnfractionIrr                      Input, irrigation_returnfraction, the fraction of non-efficient water   frac 
    alphaDepletion                         Input, alphaDepletion, irrigation aims to alphaDepletion of field capa  frac 
    minimum_irrigation                     Cover-specific irrigation in metres is 0 if less than this, currently   1/m2 
    pot_irrConsumption                     Cover-specific potential irrigation consumption                         m/m  
    fraction_IncreaseIrrigation_Nonpaddy   Input, fraction_IncreaseIrrigation_Nonpaddy, scales pot_irrConsumption  frac 
    irrPaddyDemand                         Paddy irrigation demand                                                 m    
    availWaterInfiltration                 quantity of water reaching the soil after interception, more snowmelt   m    
    ws1                                    Maximum storage capacity in layer 1                                     m    
    ws2                                    Maximum storage capacity in layer 2                                     m    
    wfc1                                   Soil moisture at field capacity in layer 1                              --   
    wfc2                                   Soil moisture at field capacity in layer 2                              --   
    wwp1                                   Soil moisture at wilting point in layer 1                               --   
    wwp2                                   Soil moisture at wilting point in layer 2                               --   
    arnoBeta                                                                                                       --   
    maxtopwater                            maximum heigth of topwater                                              m    
    totAvlWater                            Field capacity minus wilting point in soil layers 1 and 2               m    
    InvCellArea                            Inverse of cell area of each simulated mesh                             1/m2 
    totalPotET                             Potential evaporation per land use class                                m    
    w1                                     Simulated water storage in the layer 1                                  m    
    w2                                     Simulated water storage in the layer 2                                  m    
    fracVegCover                           Fraction of specific land covers (0=forest, 1=grasslands, etc.)         %    
    unmetDemand                            Unmet groundwater demand to determine potential fossil groundwaterwate  m    
    unmetDemandPaddy                       Unmet paddy demand                                                      m    
    unmetDemandNonpaddy                    Unmet nonpaddy demand                                                   m    
    irrDemand                              Cover-specific Irrigation demand                                        m/m  
    irrNonpaddyDemand                                                                                              --   
    totalIrrDemand                         Irrigation demand                                                       m    
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

    def initial(self):
        """
        Initial part of the water demand module
        irrigation

        """

        # init unmetWaterDemand -> to calculate actual one the the unmet water demand from previous day is needed
        self.var.unmetDemandPaddy = self.var.load_initial('unmetDemandPaddy', default=globals.inZero.copy())
        self.var.unmetDemandNonpaddy = self.var.load_initial('unmetDemandNonpaddy', default=globals.inZero.copy())
        # in case fossil water abstraction is allowed this will be filled
        self.var.unmetDemand = globals.inZero.copy()
        self.var.unmetDemand_runningSum = globals.inZero.copy()

        # irrigation efficiency
        # at the moment a single map, but will be replaced by map stack for every year
        self.var.efficiencyPaddy = loadmap("irrPaddy_efficiency")
        self.var.efficiencyNonpaddy = loadmap("irrNonPaddy_efficiency")
        self.var.returnfractionIrr = loadmap("irrigation_returnfraction")
        self.var.minCropKC= loadmap('minCropKC')

        self.var.alphaDepletion = loadmap('alphaDepl')

        # ignore demand if less than self.var.minimum_irrigation #1 m3
        self.var.minimum_irrigation = self.var.InvCellArea

    def dynamic(self):
        """
        Dynamic part of the water demand module

        * calculate the fraction of water from surface water vs. groundwater
        * get non-Irrigation water demand and its return flow fraction
        """
        
        if dateVar['newStart'] or (dateVar['currDate'].day in [1,11,21]):
            coverType = 'irrPaddy'
            self.var.cropKC[2] = readnetcdf2(coverType + '_cropCoefficientNC', dateVar['10day'], "10day")
            coverType = 'irrNonPaddy'
            self.var.cropKC[3] = readnetcdf2(coverType + '_cropCoefficientNC', dateVar['10day'], "10day")
            for No in [2,3]:
                self.var.cropKC[No] = np.maximum(self.var.cropKC[No], self.var.minCropKC)
                self.var.cropKC_landCover[No] = self.var.cropKC[No].copy()

        # -----------------
        # irrNonPaddy No=3
        # irrPaddy will follow
        # -----------------
        No = 3


        # load soil parameters
        soilWaterStorageCap = loadmap("wsp")
        self.var.wfcirr = loadmap("wfcp")
        self.var.wwpirr = loadmap("wwpp")
        
        # Non-paddy irrigation
        # Total available water (field capacity-wilting point)
        self.var.totAvlWater = np.maximum(0., self.var.wfcirr - self.var.wwpirr)
        
        # Compute water deficit
        #  ========================================
        # readily available water in rootzone
        
        relSat = self.var.rootzoneSM
        readAvlWater = (soilWaterStorageCap*relSat) - self.var.wwpirr


        # Potential infiltration capacity
        #  ========================================
        # depends on relSat and arnoBeta
        
        satAreaFrac = np.maximum(1 - (1 - relSat),0) ** self.var.arnoBeta[No]
        satAreaFrac = np.maximum(np.minimum(satAreaFrac, 1.0), 0.0)
        store = soilWaterStorageCap / (self.var.arnoBeta[No] + 1)
        potBeta = (self.var.arnoBeta[No] + 1) / self.var.arnoBeta[No]
        potInf = store - store * (1 - (1 - satAreaFrac) ** potBeta)


        # Potential irrigation amount
        #  ========================================
        # Consider cropKC to determine irrigation timing
        
        self.var.pot_irrConsumption[No] = np.where(self.var.cropKC[No] > 0.20, np.maximum(0.0, self.var.alphaDepletion * self.var.totAvlWater -readAvlWater), 0.)

        # should not be bigger than infiltration capacity
        self.var.pot_irrConsumption[No] = np.minimum(self.var.pot_irrConsumption[No], potInf)

        # ignore demand if less than self.var.minimum_irrigation
        self.var.pot_irrConsumption[No] = np.where(self.var.pot_irrConsumption[No] > self.var.minimum_irrigation, self.var.pot_irrConsumption[No], 0)
        self.var.irrDemand[No] = self.var.pot_irrConsumption[No] / self.var.efficiencyNonpaddy
        
        # Crop coefficients
        self.var.cropkcpad = self.var.cropKC[2]
        self.var.cropkcnpad = self.var.cropKC[3]
        
        # Sum up irrigation water demand with area fraction
        if checkOption('paddy_irrig'): # set paddy to nonpaddy (implemented later)
            self.var.irrNonpaddyDemand = self.var.fracVegCover[3] * self.var.irrDemand[3]
            self.var.irrPaddyDemand = self.var.fracVegCover[2] * self.var.irrDemand[2]
        else:
            self.var.irrNonpaddyDemand = (self.var.fracVegCover[2] + self.var.fracVegCover[3]) * self.var.irrDemand[3]
            self.var.irrPaddyDemand = globals.inZero.copy()
            
        self.var.totalIrrDemand = self.var.irrPaddyDemand + self.var.irrNonpaddyDemand
