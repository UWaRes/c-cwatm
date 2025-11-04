# -------------------------------------------------------------------------
# Name:        Land Cover Type Initialisation
# Purpose:
#
# Author:      Peter Greve, Peter Burek
#
# Created:     15/07/2016
# Copyright:   (c) Peter Burek 2016, (c) Peter Greve 2023
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *

def decompress(map, nanvalue=None):
    """
    Decompressing CWatM maps from 1D to 2D with missing values

    :param map: compressed map
    :return: decompressed 2D map
    """

    dmap = maskinfo['maskall'].copy()
    dmap[~maskinfo['maskflat']] = map[:]
    if nanvalue is not None:
        dmap.data[np.isnan(dmap.data)] = nanvalue

    return dmap.data

class landcoverType(object):

    """
    LAND COVER TYPE

    runs the 6 land cover types through soil procedures

    This routine calls the soil routine for each land cover type


    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  ===== 
    snowEvap                               total evaporation from snow for a snow layers                           m    
    load_initial                           Settings initLoad holds initial conditions for variables                input
    topwater                               quantity of water above the soil (flooding)                             m    
    waterBodyID                            lakes/reservoirs map with a single ID for each lake/reservoir           --   
    compress_LR                            boolean map as mask map for compressing lake/reservoir                  --   
    decompress_LR                          boolean map as mask map for decompressing lake/reservoir                --   
    MtoM3C                                 conversion factor from m to m3 (compressed map)                         --   
    waterBodyTypTemp                                                                                               --   
    maxGWCapRise                           influence of capillary rise above groundwater level                     m    
    minCropKC                              minimum crop factor (default 0.2)                                       --   
    irrigatedArea_original                                                                                         --   
    frac_totalnonIrr                       Fraction sown with specific non-irrigated crops                         %    
    frac_totalIrr_max                      Fraction sown with specific irrigated crops, maximum throughout simula  %    
    frac_totalnonIrr_max                   Fraction sown with specific non-irrigated crops, maximum throughout si  %    
    GeneralCrop_Irr                        Fraction of irrigated land class sown with generally representative cr  %    
    fallowIrr                              Fraction of fallowed irrigated land                                     %    
    fallowIrr_max                          Fraction of fallowed irrigated land, maximum throughout simulation      %    
    GeneralCrop_nonIrr                     Fraction of grasslands sown with generally representative crop          %    
    fallownonIrr                           Fraction of fallowed non-irrigated land                                 %    
    fallownonIrr_max                       Fraction of fallowed non-irrigated land, maximum throughout simulation  %    
    availableArableLand                    Fraction of land not currently planted with specific crops              %    
    sum_gwRecharge                         groundwater recharge                                                    m    
    minInterceptCap                        Maximum interception read from file for forest and grassland land cove  m    
    interceptStor                          simulated vegetation interception storage                               m    
    availWaterInfiltration                 quantity of water reaching the soil after interception, more snowmelt   m    
    lakeStorage                                                                                                    --   
    resStorage                                                                                                     --   
    riverbedExchangeM                      Flow from channel into groundwater                                      m    
    leakageIntoGw                          Canal leakage leading to groundwater recharge                           m    
    leakageIntoRunoff                      Canal leakage leading to runoff                                         m    
    dynamicLandcover                                                                                               --   
    staticLandCoverMaps                    1=staticLandCoverMaps in settings file is True, 0=otherwise             --   
    landcoverSum                                                                                                   --   
    sum_interceptStor                      Total of simulated vegetation interception storage including all landc  m    
    minTopWaterLayer                                                                                               --   
    maxRootDepth                                                                                                   --   
    rootDepth                                                                                                      --   
    KSat1                                                                                                          --   
    KSat2                                                                                                          --   
    KSat3                                                                                                          --   
    alpha1                                                                                                         --   
    alpha2                                                                                                         --   
    alpha3                                                                                                         --   
    lambda1                                                                                                        --   
    lambda2                                                                                                        --   
    lambda3                                                                                                        --   
    thetas1                                                                                                        --   
    thetas2                                                                                                        --   
    thetas3                                                                                                        --   
    thetar1                                                                                                        --   
    thetar2                                                                                                        --   
    thetar3                                                                                                        --   
    genuM1                                                                                                         --   
    genuM2                                                                                                         --   
    genuM3                                                                                                         --   
    genuInvM1                                                                                                      --   
    genuInvM2                                                                                                      --   
    genuInvM3                                                                                                      --   
    ws1                                    Maximum storage capacity in layer 1                                     m    
    ws2                                    Maximum storage capacity in layer 2                                     m    
    ws3                                    Maximum storage capacity in layer 3                                     m    
    wres1                                  Residual storage capacity in layer 1                                    m    
    wres2                                  Residual storage capacity in layer 2                                    m    
    wres3                                  Residual storage capacity in layer 3                                    m    
    wrange1                                                                                                        --   
    wrange2                                                                                                        --   
    wrange3                                                                                                        --   
    wfc1                                   Soil moisture at field capacity in layer 1                              --   
    wfc2                                   Soil moisture at field capacity in layer 2                              --   
    wfc3                                   Soil moisture at field capacity in layer 3                              --   
    wwp1                                   Soil moisture at wilting point in layer 1                               --   
    wwp2                                   Soil moisture at wilting point in layer 2                               --   
    wwp3                                   Soil moisture at wilting point in layer 3                               --   
    kUnSat3FC                                                                                                      --   
    kunSatFC12                                                                                                     --   
    kunSatFC23                                                                                                     --   
    rootFraction1                                                                                                  --   
    cropCoefficientNC_filename                                                                                     --   
    interceptCapNC_filename                                                                                        --   
    coverFractionNC_filename                                                                                       --   
    sum_topwater                           quantity of water on the soil (flooding) (weighted sum for all landcov  m    
    sum_soil                                                                                                       --   
    sum_w1                                                                                                         --   
    sum_w2                                                                                                         --   
    sum_w3                                                                                                         --   
    totalSto                               Total soil,snow and vegetation storage for each cell including all lan  m    
    arnoBetaOro                            chosen ModFlow model timestep (1day, 7days, 30days, etc.)               --   
    arnoBeta                                                                                                       --   
    adjRoot                                                                                                        --   
    maxtopwater                            maximum heigth of topwater                                              m    
    totAvlWater                            Field capacity minus wilting point in soil layers 1 and 2               m    
    fracGlacierCover                                                                                               --   
    pretotalSto                            Previous totalSto                                                       m    
    prefFlow_GW                            Preferential flow to groundwater. sum_prefFlow goes either to groundwa  m    
    sum_prefFlow                           Preferential flow from soil to groundwater (summed up for all land cov  m    
    sum_perc3toGW                          Percolation from 3rd soil layer to groundwater (summed up for all land  m    
    perc3toGW_GW                           Percolation from 3rd soil layer to groundwater. sum_perc3toGW goes eit  m    
    riverbedExchangeM3                                                                                             --   
    lakebedExchangeM                       Flow of water from lakes and reservoirs into groundwater                m    
    sum_actBareSoilEvap                                                                                            --   
    sum_openWaterEvap                                                                                              --   
    sum_runoff                             Runoff above the soil, more interflow, including all landcover types    m    
    sum_directRunoff                                                                                               --   
    sum_interflow                                                                                                  --   
    GWVolumeVariation                                                                                              --   
    sum_availWaterInfiltration                                                                                     --   
    sum_capRiseFromGW                      Capillary rise from groundwater to 3rd soil layer (summed up for all l  m    
    sum_act_irrConsumption                                                                                         --   
    cellArea                               Area of cell                                                            m2   
    MtoM3                                  Coefficient to change units                                             --   
    InvCellArea                            Inverse of cell area of each simulated mesh                             1/m2 
    Precipitation                          Precipitation (input for the model)                                     m    
    coverTypes                             land cover types - forest - grassland - irrPaddy - irrNonPaddy - water  --   
    SnowMelt                               total snow melt from all layers                                         m    
    Rain                                   Precipitation less snow                                                 m    
    prevSnowCover                          snow cover of previous day (only for water balance)                     m    
    SnowCover                              snow cover (sum over all layers)                                        m    
    ElevationStD                                                                                                   --   
    frac_totalIrr                          Fraction sown with specific irrigated crops                             %    
    soilLayers                             Number of soil layers                                                   --   
    soildepth                              Thickness of the first soil layer                                       m    
    w1                                     Simulated water storage in the layer 1                                  m    
    w2                                     Simulated water storage in the layer 2                                  m    
    w3                                     Simulated water storage in the layer 3                                  m    
    baseflow                               simulated baseflow (= groundwater discharge to river)                   m    
    capriseindex                                                                                                   --   
    soildepth12                            Total thickness of layer 2 and 3                                        m    
    leakageriver_factor                                                                                            --   
    leakagelake_factor                                                                                             --    
    urbanleak                                                                                                      --   
    fracVegCover                           Fraction of specific land covers (0=forest, 1=grasslands, etc.)         %    
    lakeVolumeM3C                          compressed map of lake volume                                           m3   
    lakeStorageC                                                                                                   --   
    reservoirStorageM3C                                                                                            --   
    lakeResStorageC                                                                                                --   
    lakeResStorage                                                                                                 --   
    act_SurfaceWaterAbstract               Surface water abstractions                                              m    
    readAvlChannelStorageM                                                                                         --   
    leakageCanals_M                                                                                                --   
    addtoevapotrans                        Irrigation application loss to evaporation                              m    
    act_irrWithdrawal                      Irrigation withdrawals                                                  m    
    act_nonIrrConsumption                  Non-irrigation consumption                                              m    
    returnFlow                                                                                                     --   
    totalET                                Total evapotranspiration for each cell including all landcover types    m    
    sum_actTransTotal                                                                                              --   
    sum_interceptEvap                                                                                              --   
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

    # noinspection PyTypeChecker
    def initial(self):
        """
        Initial part of the land cover type module
        Initialise the six land cover types

        * Forest No.0
        * Grasland/non irrigated land No.1
        * Paddy irrigation No.2
        * non-Paddy irrigation No.3
        * Sealed area No.4
        * Water covered area No.5

        And initialize the soil variables
        """
        # Load topoStD
        self.var.ElevationStD = loadmap('ElevationStD')

        
        # make land cover change from year to year or fix it to 1 year
        if returnBool('dynamicLandcover'):
            self.var.dynamicLandcover = True
        else:
            self.var.dynamicLandcover = False

        self.var.staticLandCoverMaps = False
        if "staticLandCoverMaps" in option:
            self.var.staticLandCoverMaps = checkOption('staticLandCoverMaps')
        
        self.var.coverTypes= list(map(str.strip, cbinding("coverTypes").split(",")))
        landcoverAll = ['fracVegCover','interceptStor']
        for variable in landcoverAll:  vars(self.var)[variable] = np.tile(globals.inZero, (6, 1))

        landcoverPara = ['minInterceptCap']
        # arrays stored as list not as numpy, because it can contain strings, single parameters or arrays
        # list is filled with append afterwards
        for variable in landcoverPara: vars(self.var)[variable] = []

        # fraction (m2) of a certain irrigation type over (only) total irrigation area ; will be assigned by the landSurface module
        # output variable per land cover class
        landcoverVars = ['totAvlWater','cropKC', 'cropKC_landCover','pot_irrConsumption','irrDemand','act_irrConsumption']
     
        # for 6 landcover types
        for variable in landcoverVars:  vars(self.var)[variable] = np.tile(globals.inZero,(6,1))


        # arnoBeta just for irrigation land cover types
        self.var.arnoBeta = globals.inZero.copy()

        soilChar = ['ws1','ws2','wfc1','wfc2','wwp1','wwp2']
        # For 3 soil layers and 4 landcover types
        for variable in soilChar:  vars(self.var)[variable]= globals.inZero.copy()
        
        # set aggregated storages to zero
        self.var.landcoverSum = ['interceptStor']
        for variable in self.var.landcoverSum: vars(self.var)["sum_"+variable] = globals.inZero.copy()

         # ----------------------------------------------------------
        # Load initial values and calculate basic soil parameters which are not changed in time

        self.var.act_SurfaceWaterAbstract = globals.inZero.copy()
        self.var.irrigatedArea_original = globals.inZero.copy()
        
        self.dynamic_fracIrrigation(init=True, dynamic = True)
        i = 0
        for coverType in self.var.coverTypes:
            self.var.minInterceptCap.append(loadmap(coverType + "_minInterceptCap"))
            # init values
            if coverType in ['forest', 'grassland', 'irrPaddy', 'irrNonPaddy','sealed']:
                self.var.interceptStor[i] = self.var.load_initial(coverType + "_interceptStor")

            # summarize the following initial storages:
            self.var.sum_interceptStor += self.var.fracVegCover[i] * self.var.interceptStor[i]
            
            i += 1


        # use non-paddy irrigation arnoBeta for both irrigation types
        self.var.arnoBetaOro = (self.var.ElevationStD - 10.0) / (self.var.ElevationStD + 1500.0)

        self.var.arnoBetaOro = self.var.arnoBetaOro + loadmap('arnoBeta_add')
        self.var.arnoBeta = np.minimum(1.2, np.maximum(0.01, self.var.arnoBetaOro))
        
        
        self.var.minCropKC= loadmap('minCropKC')
        self.var.minTopWaterLayer = loadmap("minTopWaterLayer")
        self.var.maxGWCapRise = loadmap("maxGWCapRise")
    # --------------------------------------------------------------------------

    def dynamic_fracIrrigation(self, init = False, dynamic = True):
        """
        Dynamic part of the irrigation land cover type module

        Calculating fraction of land cover

        * loads the fraction of landcover for each year from netcdf maps
        * calculate the fraction of 6 land cover types based on the maps
        * if used add glacier maps

        :param init: (optional) True: set for the first time of a run
        :param dynamic: used in the dynmic run not in the initial phase
        :return: -

        """

        # updating fracVegCover of landCover (for historical irrigation areas, done at yearly basis)
        # if first day of the year or first day of run

        if init and dynamic:

            if self.var.staticLandCoverMaps:

                self.var.fracVegCover[0] = loadmap('forest_fracVegCover')
                self.var.fracVegCover[2] = loadmap('irrPaddy_fracVegCover')
                self.var.fracVegCover[3] = loadmap('irrNonPaddy_fracVegCover')
                self.var.fracVegCover[4] = loadmap('sealed_fracVegCover')
                self.var.fracVegCover[5] = loadmap('water_fracVegCover')

            else:
                if self.var.dynamicLandcover:
                    landcoverYear = dateVar['currDate']
                else:
                    landcoverYear = datetime.datetime(int(binding['fixLandcoverYear']), 1, 1)

                i = 0
                for coverType in self.var.coverTypes:

                    self.var.fracVegCover[i] = readnetcdf2('fractionLandcover', landcoverYear, useDaily="yearly",  value= 'frac'+coverType)
                    i += 1

                if 'static_irrigation_map' in option:
                    if checkOption('static_irrigation_map'):
                        self.var.fracVegCover[3] = loadmap('irrNonPaddy_fracVegCover')


            # correction of grassland if sum is not 1.0
            sum = np.sum(self.var.fracVegCover,axis=0)
            self.var.fracVegCover[1] = np.maximum(0.,self.var.fracVegCover[1] + 1.0 - sum)
            sum = np.sum(self.var.fracVegCover, axis=0)
            self.var.fracVegCover[0] = np.maximum(0., self.var.fracVegCover[0] + 1.0 - sum)
            sum = np.sum(self.var.fracVegCover,axis=0)

            ### Irrigation
            self.var.irrigatedArea_original = self.var.fracVegCover[3].copy()

            # if irrigation is off every fraction of paddy and non paddy irrigation is put to land dcover 'grassland'
            if not(checkOption('includeIrrigation')):
                self.var.fracVegCover[1] = self.var.fracVegCover[1] + self.var.fracVegCover[2] + self.var.fracVegCover[3]
                self.var.fracVegCover[2] = 0.0
                self.var.fracVegCover[3] = 0.0


