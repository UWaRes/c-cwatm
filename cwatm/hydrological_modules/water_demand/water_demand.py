# -------------------------------------------------------------------------
# Name:        Water demand module
#
# Author:      PB, MS, LG, JdeB, DF
#
# Created:     15/07/2016
# Copyright:   (c) PB 2016
# -------------------------------------------------------------------------

import numpy as np
from cwatm.management_modules import globals

from cwatm.management_modules.replace_pcr import npareatotal, npareamaximum
from cwatm.management_modules.data_handling import returnBool, binding, cbinding, loadmap, divideValues, checkOption, \
    npareaaverage, readnetcdf2
from cwatm.hydrological_modules.water_demand.domestic import waterdemand_domestic
from cwatm.hydrological_modules.water_demand.industry import waterdemand_industry
from cwatm.hydrological_modules.water_demand.livestock import waterdemand_livestock
from cwatm.hydrological_modules.water_demand.irrigation import waterdemand_irrigation
from cwatm.hydrological_modules.water_demand.environmental_need import waterdemand_environmental_need

# PB1507
from cwatm.management_modules.data_handling import *

# processes: 
# includeDesal - desalination, only allowed with sectorSourceAbstractionFractions
# reservoir_transfers - inter-basin transfers, provided by Excel sheet


class water_demand:
    """
    WATERDEMAND

    Calculating water demand and attributing sources to satisfy demands
    Industrial, domestic, and livestock are based on precalculated maps
    Agricultural water demand based on water need by plants
    
    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    readAvlStorGroundwater                 same as storGroundwater but equal to 0 when inferior to a treshold      m    
    includeDesal                                                                                                   --   
    unlimitedDesal                                                                                                 --   
    desalAnnualCap                                                                                                 --   
    reservoir_transfers                    [['Giving reservoir'][i], ['Receiving reservoir'][i], ['Fraction of li  array
    loadInit                               Flag: if true initial conditions are loaded                             --   
    efficiencyPaddy                        Input, irrPaddy_efficiency, paddy irrigation efficiency, the amount of  frac 
    efficiencyNonpaddy                     Input, irrNonPaddy_efficiency, non-paddy irrigation efficiency, the am  frac 
    returnfractionIrr                      Input, irrigation_returnfraction, the fraction of non-efficient water   frac 
    irrPaddyDemand                         Paddy irrigation demand                                                 m    
    compress_LR                            boolean map as mask map for compressing lake/reservoir                  --   
    decompress_LR                          boolean map as mask map for decompressing lake/reservoir                --   
    waterBodyID_C                                                                                                  --   
    resYearC                               Compressed map of resYear                                               --   
    waterBodyTyp_unchanged                                                                                         --   
    resVolumeC                             compressed map of reservoir volume                                      Milli
    resId_restricted                                                                                               --   
    waterBodyBuffer                                                                                                --    
    reservoir_transfers_net_M3             net reservoir transfers, after exports, transfers, and imports          m3   
    reservoir_transfers_in_M3              water received into reservoirs                                          m3   
    reservoir_transfers_out_M3             water given from reservoirs                                             m3   
    reservoir_transfers_from_outside_M3    water received into reservoirs from Outside                             m3   
    reservoir_transfers_to_outside_M3      water given from reservoirs to the Outside                              m3   
    MtoM3C                                 conversion factor from m to m3 (compressed map)                         --   
    waterBodyTypCTemp                                                                                              --   
    pot_livestockConsumption                                                                                       --   
    cellArea                               Area of cell                                                            m2   
    MtoM3                                  Coefficient to change units                                             --   
    InvDtSec                                                                                                       --   
    InvCellArea                            Inverse of cell area of each simulated mesh                             1/m2 
    M3toM                                  Coefficient to change units                                             --     
    fracVegCover                           Fraction of specific land covers (0=forest, 1=grasslands, etc.)         %    
    nonFossilGroundwaterAbs                Non-fossil groundwater abstraction.                                     m    
    reservoir_transfers_net_M3C                                                                                    --   
    reservoir_transfers_in_M3C                                                                                     --   
    reservoir_transfers_out_M3C                                                                                    --   
    reservoir_transfers_from_outside_M3C                                                                           --   
    reservoir_transfers_to_outside_M3C                                                                             --   
    lakeVolumeM3C                          compressed map of lake volume                                           m3   
    lakeStorageC                                                                                                   --   
    reservoirStorageM3C                                                                                            --   
    lakeResStorageC                                                                                                --   
    lakeResStorage                                                                                                 --   
    smalllakeVolumeM3                                                                                              --   
    smalllakeStorage                                                                                               --   
    act_SurfaceWaterAbstract               Surface water abstractions                                              m    
    readAvlChannelStorageM                                                                                         --   
    leakageCanals_M                                                                                                --   
    addtoevapotrans                        Irrigation application loss to evaporation                              m    
    act_irrWithdrawal                      Irrigation withdrawals                                                  m    
    act_nonIrrConsumption                  Non-irrigation consumption                                              m    
    returnFlow                                                                                                     --   
    act_irrConsumption                     actual irrigation water consumption                                     m    
    act_irrNonpaddyWithdrawal              non-paddy irrigation withdrawal                                         m    
    act_irrPaddyWithdrawal                 paddy irrigation withdrawal                                             m    
    unmetDemand                            Unmet groundwater demand to determine potential fossil groundwaterwate  m    
    act_nonIrrWithdrawal                   Non-irrigation withdrawals                                              m    
    returnflowIrr                                                                                                  --   
    nonIrrReturnFlowFraction                                                                                       --   
    unmet_lost                             Fossil water that disappears instead of becoming return flow            m    
    channelStorage                         Channel water storage                                                   m3   
    act_totalWaterWithdrawal               Total water withdrawals                                                 m    
    act_bigLakeResAbst                     Abstractions to satisfy demands from lakes and reservoirs               m    
    act_smallLakeResAbst                   Abstractions from small lakes at demand location                        m    
    waterdemandFixed                                                                                               --   
    domesticDemand                         Domestic demand                                                         m      
    demand_unit                                                                                                    --   
    pot_domesticConsumption                                                                                        --   
    sectorSourceAbstractionFractions                                                                               --   
    dom_efficiency                                                                                                 --   
    envFlow                                                                                                        --   
    industryDemand                                                                                                 --   
    pot_industryConsumption                                                                                        --   
    ind_efficiency                                                                                                 --   
    unmetDemandPaddy                       Unmet paddy demand                                                      m    
    unmetDemandNonpaddy                    Unmet nonpaddy demand                                                   m    
    irrDemand                              Cover-specific Irrigation demand                                        m/m  
    irrNonpaddyDemand                                                                                              --   
    totalIrrDemand                         Irrigation demand                                                       m    
    livestockDemand                                                                                                --   
    liv_efficiency                                                                                                 --    
    includeIndusDomesDemand                Input, True if includeIndusDomesDemand = True                           bool 
    irrWithdrawalSW_max                                                                                            --   
    irrWithdrawalGW_max                                                                                            --   
    relax_abstraction_fraction_initial                                                                             --   
    waterdemandFixedYear                                                                                           --   
    swAbstractionFraction_Channel_Livesto  Input, Fraction of Livestock demands to be satisfied from Channels      %    
    swAbstractionFraction_Channel_Industr  Input, Fraction of Industrial water demand to be satisfied by Channels  %    
    swAbstractionFraction_Channel_Irrigat  Input, Fraction of Irrigation demand to be satisfied from Channels      %    
    swAbstractionFraction_Lake_Livestock   Input, Fraction of Livestock water demands to be satisfied by Lakes     %    
    swAbstractionFraction_Lake_Industry    Input, Fraction of Industrial water demand to be satisfied by Lakes     %    
    swAbstractionFraction_Lake_Irrigation  Input, Fraction of Irrigation demand to be satisfied by Lakes           %    
    swAbstractionFraction_Res_Livestock    Input, Fraction of Livestock water demands to be satisfied by Reservoi  %    
    swAbstractionFraction_Res_Industry     Input, Fraction of Industrial water demand to be satisfied by Reservoi  %    
    swAbstractionFraction_Res_Irrigation   Input, Fraction of Irrigation demand to be satisfied by Reservoirs      %    
    othAbstractionFraction_Desal_Domestic                                                                          --   
    othAbstractionFraction_Desal_Livestoc                                                                          --   
    othAbstractionFraction_Desal_Industry                                                                          --   
    othAbstractionFraction_Desal_Irrigati                                                                          --   
    gwAbstractionFraction_Livestock        Fraction of livestock water demand to be satisfied by groundwater       %    
    gwAbstractionFraction_Industry         Fraction of industrial water demand to be satisfied by groundwater      %    
    gwAbstractionFraction_Irrigation       Fraction of irrigation water demand to be satisfied by groundwater      %    
    using_reservoir_command_areas          True if using_reservoir_command_areas = True, False otherwise           bool 
    load_command_areas                                                                                             --   
    reservoir_command_areas                                                                                        --   
    Water_conveyance_efficiency                                                                                    --   
    segmentArea                                                                                                    --   
    canals                                                                                                         --   
    canalsArea                                                                                                     --   
    canalsAreaC                                                                                                    --   
    swAbstractionFraction_Lift_Livestock   Input, Fraction of Livestock water demands to be satisfied from Lift a  %    
    swAbstractionFraction_Lift_Industry    Input, Fraction of Industrial water demand to be satisfied from Lift a  %    
    swAbstractionFraction_Lift_Irrigation  Input, Fraction of Irrigation demand to be satisfied from Lift areas    %    
    using_lift_areas                       True if using_lift_areas = True in Settings, False otherwise            bool 
    lift_command_areas                                                                                             --   
    allocSegments                                                                                                  --   
    swAbstractionFraction                  Input, Fraction of demands to be satisfied with surface water           %    
    allocation_zone                                                                                                --   
    leakage                                Canal leakage leading to either groundwater recharge or runoff          m3   
    pumping                                                                                                        --   
    allowedPumping                                                                                                 --   
    ratio_irrWithdrawalGW_month                                                                                    --   
    ratio_irrWithdrawalSW_month                                                                                    --   
    Desal_Domestic                                                                                                 --   
    Desal_Industry                                                                                                 --   
    Desal_Livestock                                                                                                --   
    Desal_Irrigation                                                                                               --   
    Channel_Domestic                       Channel water abstracted for domestic                                   m    
    Channel_Industry                       Channel water abstracted for industry                                   m    
    Channel_Livestock                      Channel water abstracted for livestock                                  m    
    Channel_Irrigation                     Channel water abstracted for irrigation                                 m    
    Lift_Domestic                                                                                                  --   
    Lift_Industry                                                                                                  --   
    Lift_Livestock                                                                                                 --   
    Lift_Irrigation                                                                                                --   
    Lake_Domestic                                                                                                  --   
    Lake_Industry                                                                                                  --   
    Lake_Livestock                                                                                                 --   
    Lake_Irrigation                                                                                                --   
    Res_Domestic                                                                                                   --   
    Res_Industry                                                                                                   --   
    Res_Livestock                                                                                                  --   
    Res_Irrigation                                                                                                 --   
    GW_Domestic                            Groundwater withdrawals to satisfy domestic water requests              m    
    GW_Industry                            Groundwater withdrawals to satisfy Industrial water requests            m    
    GW_Livestock                           Groundwater withdrawals to satisfy Livestock water requests             m    
    GW_Irrigation                          Groundwater withdrawals for Irrigation                                  m    
    abstractedLakeReservoirM3              Abstractions from lakes and reservoirs at the location of the waterbod  m3   
    act_DesalWaterAbstractM                                                                                        --   
    act_totalIrrConsumption                Total irrigation consumption                                            m    
    act_totalWaterConsumption              Total water consumption                                                 m    
    act_indConsumption                     Industrial consumption                                                  m    
    act_domConsumption                     Domestic consumption                                                    m    
    act_livConsumption                     Livestock consumptions                                                  m    
    returnflowNonIrr                                                                                               --   
    nonIrruse                                                                                                      --   
    act_indDemand                          Industrial demand                                                       m    
    act_domDemand                          Domestic demand                                                         m    
    act_livDemand                          Livestock demands                                                       m    
    nonIrrDemand                                                                                                   --   
    totalWaterDemand                       Irrigation and non-irrigation demand                                    m    
    act_indWithdrawal                      Industrial withdrawal                                                   m    
    act_domWithdrawal                      Domestic withdrawal                                                     m    
    act_livWithdrawal                      Livestock withdrawals                                                   m    
    pot_GroundwaterAbstract                Potential groundwater abstraction.                                      m    
    WB_elec                                Fractions of live storage to be exported from basin                     366-d
    act_nonpaddyConsumption                Non-paddy irrigation consumption                                        m    
    act_paddyConsumption                   Paddy consumption                                                       m    
    pot_nonIrrConsumption                                                                                          --   
    act_DesalWaterAbstractM3                                                                                       --   
    AvlDesalM3                                                                                                     --   
    act_channelAbst                        Abstractions to satisfy demands from channels                           m    
    metRemainSegment_lift                                                                                          --   
    act_channelAbstract_Lift               Abstractions from the channel in lift areas at the location of the cha  m    
    abstractedLakeReservoirM3C             Compressed abstractedLakeReservoirM3                                    m3   
    remainNeed                                                                                                     --   
    act_lakeAbst                           Abstractions from lakes at demand location                              m    
    inZero_C                                                                                                       --   
    swAbstractionFraction_nonIrr           Input, Fraction of non-irrigation demands to be satisfied with surface  %    
    act_ResAbst                            Abstractions from reservoirs at demand location                         m    
    leakageC_daily                                                                                                 --   
    leakageCanalsC_M                                                                                               --   
    act_irrPaddyDemand                     paddy irrigation demand                                                 m    
    act_irrNonpaddyDemand                  non-paddy irrigation demand                                             m    
    Channel_Domestic_fromZone                                                                                      --   
    Channel_Livestock_fromZone                                                                                     --   
    Channel_Industry_fromZone                                                                                      --   
    Channel_Irrigation_fromZone                                                                                    --   
    GW_Domestic_fromZone                                                                                           --   
    GW_Livestock_fromZone                                                                                          --   
    GW_Industry_fromZone                                                                                           --   
    GW_Irrigation_fromZone                                                                                         --   
    unmet_lostirr                          Fossil water for irrigation that disappears instead of becoming return  m    
    unmet_lostNonirr                       Fossil water for non-irrigation that disappears instead of becoming re  m    
    waterabstraction                                                                                               --   
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

        self.domestic = waterdemand_domestic(model)
        self.industry = waterdemand_industry(model)
        self.livestock = waterdemand_livestock(model)
        self.irrigation = waterdemand_irrigation(model)
        self.environmental_need = waterdemand_environmental_need(model)

    def initial(self):
        """
        Initial part of the water demand module
        """

        if checkOption('includeWaterDemand'):
            self.var.includeIndusDomesDemand = True
        else:
            self.var.includeIndusDomesDemand = False

        # True if all demands are taken into account,
        # False if not only irrigation is considered
        # This variable has no impact if includeWaterDemand is False
        if "includeIndusDomesDemand" in option:
            self.var.includeIndusDomesDemand = checkOption('includeIndusDomesDemand')

        if checkOption('includeWaterDemand'):

            self.irrigation.initial()
            self.environmental_need.initial()
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.domestic.initial()
                self.industry.initial()
                self.livestock.initial()

            # if waterdemand is fixed it means it does not change between years.
            self.var.waterdemandFixed = False
            if "waterdemandFixed" in binding:
                if returnBool('waterdemandFixed'):
                    self.var.waterdemandFixed = True
                    self.var.waterdemandFixedYear = loadmap('waterdemandFixedYear')

            self.var.sectorSourceAbstractionFractions = False
            # Sector-,source-abstraction fractions facilitate designating the specific source for the specific sector
            # Sources: River, Lake, Reservoir, Groundwater, Sectors: Domestic, Industry, Livestock, Irrigation
            # Otherwise, one can distinguish only between surface and groundwater, irrigation and non-irrigation

            # initializations for sectorSourceAbstractionFractions (if True)
            self.init_SourceAbstractionFractions()

            self.var.using_reservoir_command_areas = False
            self.var.load_command_areas = False

            # initialization of water bodies related variables (if True)
            self.init_WaterBodies()

            # initialization of lift areas related variables (if True)
            self.init_lift_areas()

            # initialization of surface water abstraction fraction
            self.init_swAbstractionFraction()

            # non-irrigation input maps have for each month or year the unit m/day (True) or million m3/month (False)
            self.var.demand_unit = True
            if "demand_unit" in binding:
                self.var.demand_unit = returnBool('demand_unit')

            # initialization of allocation area
            self.init_allocation_area()

            # initialization of water demand related variables
            self.init_vars_WaterDemand()

        else:  # no water demand
            # initialization of water demand related variables
            self.init_vars_noWaterDemand()


    # --------------------------------------------------------------------
    # ====================================================================
    # --------------------------------------------------------------------

    
    def dynamic(self):
        """
        Dynamic part of the water demand module

        * calculate the fraction of water from surface water vs. groundwater
        * get non-Irrigation water demand and its return flow fraction
        """

        if checkOption('includeWaterDemand'):

            # ----------------------------------------------------
            # WATER DEMAND
            # ----------------------------------------------------

            # Fix year of water demand on predefined year
            wd_date = globals.dateVar['currDate']
            if self.var.waterdemandFixed:
                wd_date = wd_date.replace(day=1)
                wd_date = wd_date.replace(year=self.var.waterdemandFixedYear)

            # calculate water demand for each sector
            self.irrigation.dynamic()
            self.environmental_need.dynamic()
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.domestic.dynamic(wd_date)
                self.industry.dynamic(wd_date)
                self.livestock.dynamic(wd_date)

            # calculate total water demand
            totalDemand, frac_industry, frac_domestic, frac_livestock =  calculate_total_water_demand(self) 


            # ----------------------------------------------------
            # WATER DEMAND vs. WATER AVAILABILITY
            # ----------------------------------------------------

            # conversion m3 -> m # minus environmental flow
            self.var.readAvlChannelStorageM = np.maximum(0.,
                                                         self.var.channelStorage * self.var.M3toM - self.var.envFlow)  # in [m]

            if dateVar['newStart'] or dateVar['newMonth']:
                self.var.act_irrWithdrawalSW_month = globals.inZero.copy()
                self.var.act_irrWithdrawalGW_month = globals.inZero.copy()

                if 'commandAreasRelaxGwAbstraction' in binding and self.var.sectorSourceAbstractionFractions:
                    self.var.gwAbstractionFraction_Irrigation = np.where(self.var.reservoir_command_areas > 0,
                                                                         0.01,
                                                                         self.var.gwAbstractionFraction_Irrigation)

            if self.var.sectorSourceAbstractionFractions and 'commandAreasRelaxGwAbstraction' in binding and \
                    self.var.using_reservoir_command_areas:

                if dateVar['currDate'].day > 15:
                    self.var.gwAbstractionFraction_Irrigation = np.where(self.var.reservoir_command_areas > 0,
                                                                         loadmap('commandAreasRelaxGwAbstraction'),
                                                                         self.var.gwAbstractionFraction_Irrigation)
            
            # Desalination
            self.calculate_desalination_abstraction(dateVar)

            # Surface water abstraction from channel storage
            self.calculate_channel_abstraction(totalDemand)   
            
            # UNDER CONSTRUCTION
            if self.var.using_lift_areas:
                # Calculate abstraction of surface water via lift systems within command areas
                self.calculate_lift_abstraction()

            if checkOption('includeWaterBodies'):

                # Abstraction from lakes and reservoirs 
                self.calculate_lake_reservoir_abstraction()

                # Transfer water between reservoirs
                # Send storage between reservoirs using the Excel sheet reservoir_transfers within cwatm_settings.xlsx
                # Using the waterBodyIDs defined in the settings, designate
                # the Giver, the Receiver, and the daily fraction of live storage the Giver sends to the Receiver.
                # If the Receiver is already at capacity, the Giver does not send any storage.
                # Reservoirs can only send to one reservoir. Reservoirs can receive from several reservoirs.
                if 'reservoir_transfers' in option:
                    if checkOption('reservoir_transfers'):
                        # Process water transfers between reservoirs and lakes
                        self.process_reservoir_transfers()


                # Calculate remaining demand after surface water abstractions
                remainNeed2 = self.substract_surface_abstraction()       

                self.var.act_ResAbst = globals.inZero.copy()
                if self.var.sectorSourceAbstractionFractions:
                    # Allocate reservoir abstractions to sectors
                    self.allocate_reservoir_abstraction_sectors()

                # If sector- and source-specific abstractions are activated, then domestic, industrial, and
                # livestock demands were attempted to be satisfied in the previous step. Otherwise, total demands
                # not satisfied by previous sources is attempted.
                self.allocate_reservoir_abstraction()

            # remaining demand is taken from groundwater if possible
            self.calculate_groundwater_abstraction(totalDemand)

            if checkOption('limitAbstraction'):
                # if limitAbstraction from groundwater is True
                # fossil gwAbstraction and water demand may be reduced
                # variable to reduce/limit groundwater abstraction (> 0 if limitAbstraction = True)
                self.allocate_groundwater_limitAbstraction(totalDemand)

            else:
                # Fossil groundwater abstractions are allowed (act = pot)
                self.allocate_groundwater_abstraction()

            # calculate actual withdrawals and consumption
            self.calculate_water_consumption()

            # calculate return flows and evaporation losses
            self.calculate_return_flow()




    # ---------------------------------------------------------------------------
    # ===========================================================================
    # ---------------------------------------------------------------------------

    # ===========================================================================
    #
    # =================== functions used in initial =============================
    #
    # ===========================================================================

    def init_SourceAbstractionFractions(self):
        """
        Initialization of sector- and source-specific abstraction fractions
        """
        if 'sectorSourceAbstractionFractions' in option:
            if checkOption('sectorSourceAbstractionFractions'):
                #print('Sector- and source-specific abstraction fractions are activated (water_demand.py)')
                self.var.sectorSourceAbstractionFractions = True

                self.var.swAbstractionFraction_Channel_Domestic = loadmap(
                    'swAbstractionFraction_Channel_Domestic')
                self.var.swAbstractionFraction_Channel_Livestock = loadmap(
                    'swAbstractionFraction_Channel_Livestock')
                self.var.swAbstractionFraction_Channel_Industry = loadmap(
                    'swAbstractionFraction_Channel_Industry')
                self.var.swAbstractionFraction_Channel_Irrigation = loadmap(
                    'swAbstractionFraction_Channel_Irrigation')

                self.var.swAbstractionFraction_Lake_Domestic = loadmap(
                    'swAbstractionFraction_Lake_Domestic')
                self.var.swAbstractionFraction_Lake_Livestock = loadmap(
                    'swAbstractionFraction_Lake_Livestock')
                self.var.swAbstractionFraction_Lake_Industry = loadmap(
                    'swAbstractionFraction_Lake_Industry')
                self.var.swAbstractionFraction_Lake_Irrigation = loadmap(
                    'swAbstractionFraction_Lake_Irrigation')

                self.var.swAbstractionFraction_Res_Domestic = loadmap(
                    'swAbstractionFraction_Res_Domestic')
                self.var.swAbstractionFraction_Res_Livestock = loadmap(
                    'swAbstractionFraction_Res_Livestock')
                self.var.swAbstractionFraction_Res_Industry = loadmap(
                    'swAbstractionFraction_Res_Industry')
                self.var.swAbstractionFraction_Res_Irrigation = loadmap(
                    'swAbstractionFraction_Res_Irrigation')
                    
                if self.var.includeDesal:
                    self.var.othAbstractionFraction_Desal_Domestic = loadmap(
                        'othAbstractionFraction_Desal_Domestic')
                    self.var.othAbstractionFraction_Desal_Livestock = loadmap(
                        'othAbstractionFraction_Desal_Livestock')
                    self.var.othAbstractionFraction_Desal_Industry = loadmap(
                        'othAbstractionFraction_Desal_Industry')
                    self.var.othAbstractionFraction_Desal_Irrigation = loadmap(
                        'othAbstractionFraction_Desal_Irrigation')

                if not checkOption('limitAbstraction'):
                    self.var.gwAbstractionFraction_Domestic = 1 + globals.inZero.copy()
                    self.var.gwAbstractionFraction_Livestock = 1 + globals.inZero.copy()
                    self.var.gwAbstractionFraction_Industry = 1 + globals.inZero.copy()
                    self.var.gwAbstractionFraction_Irrigation = 1 + globals.inZero.copy()
                else:
                    self.var.gwAbstractionFraction_Domestic = loadmap(
                        'gwAbstractionFraction_Domestic')
                    self.var.gwAbstractionFraction_Livestock = loadmap(
                        'gwAbstractionFraction_Livestock')
                    self.var.gwAbstractionFraction_Industry = loadmap(
                        'gwAbstractionFraction_Industry')
                    self.var.gwAbstractionFraction_Irrigation = loadmap(
                        'gwAbstractionFraction_Irrigation')


    def init_WaterBodies(self):
        """
        Initialization of water bodies related variables
        """
        if checkOption('includeWaterBodies'):
                
            # initiate reservoir_command_areas 
            self.var.reservoir_command_areas = globals.inZero.copy()

            if 'reservoir_command_areas' in binding:
                self.var.load_command_areas = True

            self.var.Water_conveyance_efficiency = 1.0 + globals.inZero

            # load command areas
            if self.var.load_command_areas:
                self.var.reservoir_command_areas = loadmap('reservoir_command_areas').astype(int)
                self.var.reservoir_command_areas = np.where(self.var.reservoir_command_areas<0,
                                                            0,
                                                            self.var.reservoir_command_areas)
            else:
                self.var.reservoir_command_areas = self.var.waterBodyBuffer

            # Lakes/restricted reservoirs within command areas are removed from the command area
            self.var.reservoir_command_areas = np.where(self.var.waterBodyTyp_unchanged == 1,
                                                    0, np.where(self.var.resId_restricted > 0, 0, self.var.reservoir_command_areas))
            self.var.segmentArea = np.where(self.var.reservoir_command_areas > 0,
                                            npareatotal(self.var.cellArea,
                                                        self.var.reservoir_command_areas), self.var.cellArea)

            # Water abstracted from reservoirs leaks along canals related to conveyance efficiency.
            # Canals are a map where canal cells have the number of the command area they are associated with
            # Command areas without canals experience leakage equally throughout the command area

            if 'canals' in binding:
                self.var.canals = loadmap('canals').astype(int)
            else:
                self.var.canals = globals.inZero.copy().astype(int)

            # canals for reservoir conveyance and loss
            self.var.canals = np.where(self.var.canals != self.var.reservoir_command_areas, 0, self.var.canals)

            # When there are no set canals, the entire command area experiences leakage
            self.var.canals = np.where(npareamaximum(self.var.canals, self.var.reservoir_command_areas) == 0,
                                    self.var.reservoir_command_areas, self.var.canals)
            self.var.canalsArea = np.where(self.var.canals > 0, npareatotal(self.var.cellArea, self.var.canals),
                                        0)
            self.var.canalsAreaC = np.compress(self.var.compress_LR, self.var.canalsArea)


    def init_lift_areas(self):
        """
        Initialization of lift areas related variables
        """
        self.var.swAbstractionFraction_Lift_Domestic = globals.inZero.copy()
        self.var.swAbstractionFraction_Lift_Livestock = globals.inZero.copy()
        self.var.swAbstractionFraction_Lift_Industry = globals.inZero.copy()
        self.var.swAbstractionFraction_Lift_Irrigation = globals.inZero.copy()

        self.var.using_lift_areas = False
        if 'using_lift_areas' in option:
            if checkOption('using_lift_areas'):

                self.var.using_lift_areas = True
                self.var.lift_command_areas = loadmap('lift_areas').astype(int)

                if self.var.sectorSourceAbstractionFractions:
                    self.var.swAbstractionFraction_Lift_Domestic = loadmap(
                        'swAbstractionFraction_Lift_Domestic')
                    self.var.swAbstractionFraction_Lift_Livestock = loadmap(
                        'swAbstractionFraction_Lift_Livestock')
                    self.var.swAbstractionFraction_Lift_Industry = loadmap(
                        'swAbstractionFraction_Lift_Industry')
                    self.var.swAbstractionFraction_Lift_Irrigation = loadmap(
                        'swAbstractionFraction_Lift_Irrigation')


    def init_swAbstractionFraction(self):
        """
        Initialization of surface water abstraction fraction
        """
        # -------------------------------------------
        # partitioningGroundSurfaceAbstraction
        # partitioning abstraction sources: groundwater and surface water
        # partitioning based on local average baseflow (m3/s) and upstream average discharge (m3/s)
        # estimates of fractions of groundwater and surface water abstractions
        swAbstractionFraction = loadmap('swAbstractionFrac')

        if swAbstractionFraction < 0:

            averageBaseflowInput = loadmap('averageBaseflow')
            averageDischargeInput = loadmap('averageDischarge')
            # convert baseflow from m to m3/s
            if returnBool('baseflowInM'):
                averageBaseflowInput = averageBaseflowInput * self.var.cellArea * self.var.InvDtSec

            if checkOption('usingAllocSegments'):
                averageBaseflowInput = np.where(self.var.allocSegments > 0,
                                                npareaaverage(averageBaseflowInput, self.var.allocSegments),
                                                averageBaseflowInput)

                # averageUpstreamInput = np.where(self.var.allocSegments > 0,
                #                                npareamaximum(averageDischargeInput, self.var.allocSegments),
                #                                averageDischargeInput)

            swAbstractionFraction = np.maximum(0.0, np.minimum(1.0, averageDischargeInput / np.maximum(1e-20,
                                                                                                        averageDischargeInput + averageBaseflowInput)))
            swAbstractionFraction = np.minimum(1.0, np.maximum(0.0, swAbstractionFraction))

        # weighting by land cover fractions
        self.var.swAbstractionFraction = globals.inZero.copy()
        for No in range(4):
            self.var.swAbstractionFraction += self.var.fracVegCover[No] * swAbstractionFraction
        for No in range(4, 6):
            # The motivation is to avoid groundwater on sealed and water land classes
            # TODO: Groundwater pumping should be allowed over sealed land
            self.var.swAbstractionFraction += self.var.fracVegCover[No]


    def init_allocation_areas(self):
        """
        Initialization of allocation areas related variables.
        
        allocation zone, regular grid inside the 2d array, inner grid size
        """
        # 
        inner = 1
        if "allocation_area" in binding:
            inner = int(loadmap('allocation_area'))

        latldd, lonldd, cell, invcellldd, rows, cols = readCoord(cbinding('Ldd'))
        filename = os.path.splitext(cbinding('Ldd'))[0] + '.nc'
        if os.path.isfile(filename):
            cut0, cut1, cut2, cut3 = mapattrNetCDF(filename, check=False)
        else:
            filename = os.path.splitext(cbinding('Ldd'))[0] + '.tif'

            if not(os.path.isfile(filename)):
                filename = os.path.splitext(cbinding('Ldd'))[0] + '.map'
                
            nf2 = gdal.Open(filename, gdalconst.GA_ReadOnly)
            cut0, cut1, cut2, cut3 = mapattrTiff(nf2)

        arr = np.kron(np.arange(rows // inner * cols // inner).reshape((rows // inner, cols // inner)),
                        np.ones((inner, inner)))
        arr = arr[cut2:cut3, cut0:cut1].astype(int)
        self.var.allocation_zone = compressArray(arr)


    def init_vars_WaterDemand(self):
        """
        Initialization of water demand related variables
        """
        self.var.leakage = globals.inZero.copy()
        self.var.pumping = globals.inZero.copy()
        self.var.Pumping_daily = globals.inZero.copy()
        self.var.allowedPumping = globals.inZero.copy()
        self.var.leakageCanals_M = globals.inZero.copy()

        self.var.act_nonIrrWithdrawal = globals.inZero.copy()
        self.var.act_irrWithdrawalSW_month = globals.inZero.copy()
        self.var.act_irrWithdrawalGW_month = globals.inZero.copy()

        self.var.Desal_Domestic = globals.inZero.copy()
        self.var.Desal_Industry = globals.inZero.copy()
        self.var.Desal_Livestock = globals.inZero.copy()
        self.var.Desal_Irrigation = globals.inZero.copy()
        
        self.var.Channel_Domestic = globals.inZero.copy()
        self.var.Channel_Industry = globals.inZero.copy()
        self.var.Channel_Livestock = globals.inZero.copy()
        self.var.Channel_Irrigation = globals.inZero.copy()

        self.var.Lift_Domestic = globals.inZero.copy()
        self.var.Lift_Industry = globals.inZero.copy()
        self.var.Lift_Livestock = globals.inZero.copy()
        self.var.Lift_Irrigation = globals.inZero.copy()

        self.var.Lake_Domestic = globals.inZero.copy()
        self.var.Lake_Industry = globals.inZero.copy()
        self.var.Lake_Livestock = globals.inZero.copy()
        self.var.Lake_Irrigation = globals.inZero.copy()

        self.var.Res_Domestic = globals.inZero.copy()
        self.var.Res_Industry = globals.inZero.copy()
        self.var.Res_Livestock = globals.inZero.copy()
        self.var.Res_Irrigation = globals.inZero.copy()

        self.var.GW_Domestic = globals.inZero.copy()
        self.var.GW_Industry = globals.inZero.copy()
        self.var.GW_Livestock = globals.inZero.copy()
        self.var.GW_Irrigation = globals.inZero.copy()
        self.var.abstractedLakeReservoirM3 = globals.inZero.copy()

        self.var.ind_efficiency = 1.
        self.var.dom_efficiency = 1.
        self.var.liv_efficiency = 1

        self.var.act_DesalWaterAbstractM = globals.inZero.copy()
        
        self.var.act_nonIrrConsumption = globals.inZero.copy()
        self.var.act_totalIrrConsumption = globals.inZero.copy()
        self.var.act_totalWaterConsumption = globals.inZero.copy()
        self.var.act_indConsumption = globals.inZero.copy()
        self.var.act_domConsumption = globals.inZero.copy()
        self.var.act_livConsumption = globals.inZero.copy()
        self.var.returnflowIrr = globals.inZero.copy()
        self.var.returnflowNonIrr = globals.inZero.copy()
        self.var.pitLatrinToGW = globals.inZero.copy()
        self.var.act_irrNonpaddyWithdrawal = globals.inZero.copy()
        self.var.act_irrPaddyWithdrawal = globals.inZero.copy()

        self.var.ratio_irrWithdrawalGW_month = globals.inZero.copy()
        self.var.ratio_irrWithdrawalSW_month = globals.inZero.copy()


    def init_vars_noWaterDemand(self):
        """
        Initialization of water demand related variables when water demand is not included
        """
        self.var.ratio_irrWithdrawalGW_month = globals.inZero.copy()
        self.var.ratio_irrWithdrawalSW_month = globals.inZero.copy()

        self.var.nonIrrReturnFlowFraction = globals.inZero.copy()
        self.var.nonFossilGroundwaterAbs = globals.inZero.copy()
        self.var.nonIrruse = globals.inZero.copy()

        self.var.act_indDemand = globals.inZero.copy()
        self.var.act_domDemand = globals.inZero.copy()
        self.var.act_livDemand = globals.inZero.copy()
        self.var.nonIrrDemand = globals.inZero.copy()
        self.var.totalIrrDemand = globals.inZero.copy()
        self.var.totalWaterDemand = globals.inZero.copy()
        self.var.act_irrWithdrawal = globals.inZero.copy()
        self.var.act_nonIrrWithdrawal = globals.inZero.copy()
        self.var.act_totalWaterWithdrawal = globals.inZero.copy()
        self.var.act_indConsumption = globals.inZero.copy()
        self.var.act_domConsumption = globals.inZero.copy()
        self.var.act_livConsumption = globals.inZero.copy()

        self.var.act_indWithdrawal = globals.inZero.copy()
        self.var.act_domWithdrawal = globals.inZero.copy()
        self.var.act_livWithdrawal = globals.inZero.copy()

        self.var.act_totalIrrConsumption = globals.inZero.copy()
        self.var.act_totalWaterConsumption = globals.inZero.copy()
        self.var.unmetDemand = globals.inZero.copy()
        self.var.unmetDemand_runningSum = globals.inZero.copy()
        self.var.addtoevapotrans = globals.inZero.copy()
        self.var.returnflowIrr = globals.inZero.copy()
        self.var.returnflowNonIrr = globals.inZero.copy()
        self.var.returnFlow = globals.inZero.copy()
        self.var.unmetDemandPaddy = globals.inZero.copy()
        self.var.unmetDemandNonpaddy = globals.inZero.copy()
        self.var.ind_efficiency = 1.
        self.var.dom_efficiency = 1.
        self.var.liv_efficiency = 1
        self.var.act_bigLakeResAbst = globals.inZero.copy()

        self.var.leakage = globals.inZero.copy()
        self.var.pumping = globals.inZero.copy()
        self.var.unmet_lost = globals.inZero.copy()
        self.var.pot_GroundwaterAbstract = globals.inZero.copy()
        self.var.leakageCanals_M = globals.inZero.copy()

        self.var.WB_elec = globals.inZero.copy()
        
        self.var.Desal_Domestic = globals.inZero.copy()
        self.var.Desal_Industry = globals.inZero.copy()
        self.var.Desal_Livestock = globals.inZero.copy()
        self.var.Desal_Irrigation = globals.inZero.copy()
        
        self.var.Channel_Domestic = globals.inZero.copy()
        self.var.Channel_Industry = globals.inZero.copy()
        self.var.Channel_Livestock = globals.inZero.copy()
        self.var.Channel_Irrigation = globals.inZero.copy()

        self.var.Lift_Domestic = globals.inZero.copy()
        self.var.Lift_Industry = globals.inZero.copy()
        self.var.Lift_Livestock = globals.inZero.copy()
        self.var.Lift_Irrigation = globals.inZero.copy()

        self.var.Lake_Domestic = globals.inZero.copy()
        self.var.Lake_Industry = globals.inZero.copy()
        self.var.Lake_Livestock = globals.inZero.copy()
        self.var.Lake_Irrigation = globals.inZero.copy()

        self.var.Res_Domestic = globals.inZero.copy()
        self.var.Res_Industry = globals.inZero.copy()
        self.var.Res_Livestock = globals.inZero.copy()
        self.var.Res_Irrigation = globals.inZero.copy()

        self.var.GW_Domestic = globals.inZero.copy()
        self.var.GW_Industry = globals.inZero.copy()
        self.var.GW_Livestock = globals.inZero.copy()
        self.var.GW_Irrigation = globals.inZero.copy()

        self.var.abstractedLakeReservoirM3 = globals.inZero.copy()

        self.var.act_nonpaddyConsumption = globals.inZero.copy()
        self.var.act_paddyConsumption = globals.inZero.copy()
        self.var.act_irrNonpaddyWithdrawal = globals.inZero.copy()
        self.var.act_irrPaddyWithdrawal = globals.inZero.copy()

        self.var.Pumping_daily = globals.inZero.copy()

        self.var.act_irrPaddyDemand = globals.inZero.copy()
        self.var.act_irrNonpaddyDemand = globals.inZero.copy()
        self.var.domesticDemand = globals.inZero.copy()
        self.var.industryDemand = globals.inZero.copy()
        self.var.livestockDemand = globals.inZero.copy()
        
        self.var.act_DesalWaterAbstractM = globals.inZero.copy()

        self.var.act_nonIrrConsumption = globals.inZero.copy()


    # ===========================================================================
    #
    # =================== functions used in dynamic =============================
    #
    # ===========================================================================

    def calculate_total_water_demand(self):
        """
        Calculates total water demand including irrigation and non-irrigation components.
        Updates internal variables for demand, consumption, and return flow fractions.
        """
        # Check if industrial and domestic water demand should be included
        if self.var.includeIndusDomesDemand:
            # Recalculate demand if it's a new time step or if reservoir transfers are involved
            if globals.dateVar['newStart'] or globals.dateVar['newMonth'] or 'reservoir_transfers' in option:
                # Total potential non-irrigation water demand
                self.var.nonIrrDemand = self.var.domesticDemand + self.var.industryDemand + self.var.livestockDemand

                # Potential non-irrigation consumption (limited by availability)
                self.var.pot_nonIrrConsumption = np.minimum( self.var.nonIrrDemand,
                    self.var.pot_domesticConsumption + self.var.pot_industryConsumption +
                    self.var.pot_livestockConsumption )

                # Fraction of return flow from domestic and industrial demand
                self.var.nonIrrReturnFlowFraction = divideValues(
                    (self.var.nonIrrDemand - self.var.pot_nonIrrConsumption),
                    self.var.nonIrrDemand )

            # Calculate fractional contributions of each non-irrigation sector
            frac_industry = divideValues(self.var.industryDemand, self.var.nonIrrDemand)
            frac_domestic = divideValues(self.var.domesticDemand, self.var.nonIrrDemand)
            frac_livestock = divideValues(self.var.livestockDemand, self.var.nonIrrDemand)

            # Total water demand includes both irrigation and non-irrigation (in [m])
            totalDemand = self.var.nonIrrDemand + self.var.totalIrrDemand

        else:
            # Only irrigation demand is considered
            totalDemand = np.copy(self.var.totalIrrDemand)

            # Reset non-irrigation demand variables to zero
            self.var.domesticDemand = globals.inZero.copy()
            self.var.industryDemand = globals.inZero.copy()
            self.var.livestockDemand = globals.inZero.copy()
            self.var.nonIrrDemand = globals.inZero.copy()
            self.var.pot_nonIrrConsumption = globals.inZero.copy()
            self.var.nonIrrReturnFlowFraction = globals.inZero.copy()
            frac_industry = globals.inZero.copy()
            frac_domestic = globals.inZero.copy()
            frac_livestock = globals.inZero.copy()

        return totalDemand, frac_industry, frac_domestic, frac_livestock


    def calculate_lift_abstraction(self):
        """
        Calculates the abstraction of surface water via lift systems within command areas.
        Water is proportionally taken from available channel storage in each cell, based on demand.
        """

        # Step 1: Calculate potential lift abstraction for each sector
        # Limited by abstraction fraction and remaining demand after desalination and channel supply
        pot_Lift_Domestic = np.minimum(
            self.var.swAbstractionFraction_Lift_Domestic * self.var.domesticDemand,
            self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic
        )
        pot_Lift_Livestock = np.minimum(
            self.var.swAbstractionFraction_Lift_Livestock * self.var.livestockDemand,
            self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock
        )
        pot_Lift_Industry = np.minimum(
            self.var.swAbstractionFraction_Lift_Industry * self.var.industryDemand,
            self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry
        )
        pot_Lift_Irrigation = np.minimum(
            self.var.swAbstractionFraction_Lift_Irrigation * self.var.totalIrrDemand,
            self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation
        )

        # Total potential lift abstraction across all sectors
        pot_liftAbst = pot_Lift_Domestic + pot_Lift_Livestock + pot_Lift_Industry + pot_Lift_Irrigation
        remainNeed_afterLocal = pot_liftAbst.copy()

        # Step 2: Aggregate demand and available water per command area
        # Each cell in a command area holds the total demand and available channel storage
        demand_Segment_lift = np.where(
            self.var.lift_command_areas > 0,
            npareatotal(remainNeed_afterLocal * self.var.cellArea, self.var.lift_command_areas),
            0
        )
        available_Segment_lift = np.where(
            self.var.lift_command_areas > 0,
            npareatotal(self.var.readAvlChannelStorageM * self.var.cellArea, self.var.lift_command_areas),
            0
        )

        # Step 3: Calculate fraction of available water used to meet demand
        frac_used_Segment_lift = np.where(
            available_Segment_lift > 0,
            np.minimum(demand_Segment_lift / available_Segment_lift, 1.0),
            0
        )

        # Step 4: Apply abstraction and update channel storage
        self.var.act_channelAbst += (frac_used_Segment_lift * self.var.readAvlChannelStorageM)

        # Step 5: Calculate how much of the segment demand was met
        metRemainSegment_lift = np.where(
            demand_Segment_lift > 0,
            divideValues(frac_used_Segment_lift * available_Segment_lift, demand_Segment_lift),
            0
        )
        self.var.metRemainSegment_lift = metRemainSegment_lift.copy()

        # Step 6: Distribute met demand back to cells
        lift_abstractions = metRemainSegment_lift * remainNeed_afterLocal
        self.var.act_SurfaceWaterAbstract += lift_abstractions

        # Step 7: Update available channel storage after abstraction
        self.var.readAvlChannelStorageM -= (frac_used_Segment_lift * self.var.readAvlChannelStorageM)
        self.var.readAvlChannelStorageM = np.where(
            self.var.readAvlChannelStorageM < 0.02,
            0,
            self.var.readAvlChannelStorageM
        )

        # Step 8: Store lift abstraction for use in other modules (e.g., riverbed exchange)
        self.var.act_channelAbstract_Lift = frac_used_Segment_lift * self.var.readAvlChannelStorageM

        # Step 9: Disaggregate lift abstractions by sector if sector fractions are used
        if self.var.sectorSourceAbstractionFractions:
            self.var.Lift_Domestic = np.minimum(lift_abstractions, pot_Lift_Domestic)
            self.var.Lift_Livestock = np.minimum(
                lift_abstractions - self.var.Lift_Domestic,
                pot_Lift_Livestock
            )
            self.var.Lift_Industry = np.minimum(
                lift_abstractions - self.var.Lift_Domestic - self.var.Lift_Livestock,
                pot_Lift_Industry
            )
            self.var.Lift_Irrigation = np.minimum(
                lift_abstractions - self.var.Lift_Domestic - self.var.Lift_Livestock - self.var.Lift_Industry,
                pot_Lift_Irrigation
            )


    def calculate_desalination_abstraction(self, dateVar):
        """
        Calculates desalinated water abstraction based on sector-specific fractions and demand.
        Applies annual capacity limits unless unlimited desalination is enabled.

        Parameters:
        - dateVar: dictionary containing current date information (e.g., dateVar['currDate'].year)
        """

        # Initialize desalinated water abstraction to zero
        self.var.act_DesalWaterAbstractM3 = globals.inZero.copy()

        # Desalination is only allowed if sector-specific abstraction fractions are enabled
        if self.var.sectorSourceAbstractionFractions and self.var.includeDesal:
            # Step 1: Calculate potential desalination abstraction for each sector
            pot_Desal_Domestic = self.var.othAbstractionFraction_Desal_Domestic * self.var.domesticDemand
            pot_Desal_Livestock = self.var.othAbstractionFraction_Desal_Livestock * self.var.livestockDemand
            pot_Desal_Industry = self.var.othAbstractionFraction_Desal_Industry * self.var.industryDemand
            pot_Desal_Irrigation = self.var.othAbstractionFraction_Desal_Irrigation * self.var.totalIrrDemand

            # Total potential desalinated water abstraction
            pot_DesalAbst = pot_Desal_Domestic + pot_Desal_Livestock + pot_Desal_Industry + pot_Desal_Irrigation

            # Step 2: Apply annual desalination capacity limit if not unlimited
            if not self.var.unlimitedDesal:
                self.var.AvlDesalM3 = self.var.desalAnnualCap[dateVar['currDate'].year] / 365
                total_pot_DesalAbst = np.nansum(pot_DesalAbst * self.var.cellArea)
                abstractLimitCoeff = np.minimum(total_pot_DesalAbst, self.var.AvlDesalM3) / total_pot_DesalAbst
                self.var.act_DesalWaterAbstractM = pot_DesalAbst * abstractLimitCoeff
            else:
                # If unlimited desalination is allowed, use full potential abstraction
                self.var.act_DesalWaterAbstractM = pot_DesalAbst

            # Step 3: Distribute actual desalinated water abstraction to sectors
            self.var.Desal_Domestic = np.minimum(
                self.var.act_DesalWaterAbstractM,
                self.var.othAbstractionFraction_Desal_Domestic * self.var.domesticDemand
            )
            self.var.Desal_Livestock = np.minimum(
                self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic,
                self.var.othAbstractionFraction_Desal_Livestock * self.var.livestockDemand
            )
            self.var.Desal_Industry = np.minimum(
                self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic - self.var.Desal_Livestock,
                self.var.othAbstractionFraction_Desal_Industry * self.var.industryDemand
            )
            self.var.Desal_Irrigation = np.minimum(
                self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic - self.var.Desal_Livestock - self.var.Desal_Industry,
                self.var.othAbstractionFraction_Desal_Irrigation * self.var.totalIrrDemand
            )


    def calculate_channel_abstraction(self, totalDemand):
        """
        Calculates surface water abstraction from channel storage to meet water demand.
        Handles both sector-specific abstraction fractions and general abstraction when sectors are not specified.

        Parameters:
        - totalDemand: total water demand (irrigation + non-irrigation) [m]
        """
        # Check if sector-specific abstraction fractions are used
        if self.var.sectorSourceAbstractionFractions:
            # Calculate potential abstraction for each sector, limited by demand minus desalination
            self.var.pot_Channel_Domestic = np.minimum(
                self.var.swAbstractionFraction_Channel_Domestic * self.var.domesticDemand,
                self.var.domesticDemand - self.var.Desal_Domestic
            )
            self.var.pot_Channel_Livestock = np.minimum(
                self.var.swAbstractionFraction_Channel_Livestock * self.var.livestockDemand,
                self.var.livestockDemand - self.var.Desal_Livestock
            )
            self.var.pot_Channel_Industry = np.minimum(
                self.var.swAbstractionFraction_Channel_Industry * self.var.industryDemand,
                self.var.industryDemand - self.var.Desal_Industry
            )
            self.var.pot_Channel_Irrigation = np.minimum(
                self.var.swAbstractionFraction_Channel_Irrigation * self.var.totalIrrDemand,
                self.var.totalIrrDemand - self.var.Desal_Irrigation
            )

            # Total potential abstraction from all sectors
            pot_channelAbst = ( self.var.pot_Channel_Domestic + self.var.pot_Channel_Livestock +
                                self.var.pot_Channel_Industry + self.var.pot_Channel_Irrigation )

            # Actual abstraction is limited by available channel storage
            self.var.act_SurfaceWaterAbstract = np.minimum(self.var.readAvlChannelStorageM, pot_channelAbst)
            self.var.act_channelAbst = self.var.act_SurfaceWaterAbstract.copy()

            # Distribute actual abstraction back to sectors proportionally
            self.var.Channel_Domestic = np.minimum(
                self.var.act_channelAbst,
                self.var.swAbstractionFraction_Channel_Domestic * self.var.domesticDemand
            )
            self.var.Channel_Livestock = np.minimum(
                self.var.act_channelAbst - self.var.Channel_Domestic,
                self.var.swAbstractionFraction_Channel_Livestock * self.var.livestockDemand
            )
            self.var.Channel_Industry = np.minimum(
                self.var.act_channelAbst - self.var.Channel_Domestic - self.var.Channel_Livestock,
                self.var.swAbstractionFraction_Channel_Industry * self.var.industryDemand
            )
            self.var.Channel_Irrigation = np.minimum(
                self.var.act_channelAbst - self.var.Channel_Domestic - self.var.Channel_Livestock - self.var.Channel_Industry,
                self.var.swAbstractionFraction_Channel_Irrigation * self.var.totalIrrDemand 
            )

        else:
            # If no sector-specific fractions are used, apply a general abstraction fraction
            self.var.pot_SurfaceAbstract = totalDemand * self.var.swAbstractionFraction

            # Actual abstraction is limited by available channel storage
            self.var.act_SurfaceWaterAbstract = np.minimum(self.var.readAvlChannelStorageM, self.var.pot_SurfaceAbstract)
            self.var.act_channelAbst = self.var.act_SurfaceWaterAbstract.copy()

        # Update available channel storage after abstraction
        self.var.readAvlChannelStorageM -= self.var.act_SurfaceWaterAbstract

        
    def calculate_lake_reservoir_abstraction(self):
        """
        Calculates abstraction from large and small lakes and reservoirs to meet remaining surface water demand.
        Updates storage volumes and abstraction variables accordingly.
        """

        # Initialize abstraction from lakes and reservoirs to zero (compressed format)
        self.var.abstractedLakeReservoirM3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        # Calculate remaining surface water demand after initial abstraction
        remainNeed0 = np.maximum(self.var.pot_SurfaceAbstract - self.var.act_SurfaceWaterAbstract, 0)

        # Reinitialize abstraction from lakes/reservoirs
        self.var.abstractedLakeReservoirM3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        # Create mask for cells near water bodies (buffer zones)
        mskWtrBody_unrestricted = self.var.waterBodyBuffer > 0

        # Aggregate remaining demand around large lakes/reservoirs using buffer zones
        remainNeedBig = npareatotal(remainNeed0, self.var.waterBodyBuffer)
        remainNeedBigC = np.compress(self.var.compress_LR, remainNeedBig)

        # Get current storage for lakes and reservoirs (converted to meters)
        lakeResStorageC = np.where(
            self.var.waterBodyTypCTemp == 0, 0.,
            np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeStorageC, self.var.reservoirStorageM3C)
        ) / self.var.MtoM3C

        # Apply a conservative limit: only 98% of storage is considered available
        minlake = np.maximum(0., 0.98 * lakeResStorageC)

        # Determine actual abstraction from large lakes/reservoirs
        act_bigLakeAbstC = np.minimum(minlake, remainNeedBigC)

        # Subtract abstracted water from lake and reservoir storage (converted to m³)
        self.var.lakeStorageC -= act_bigLakeAbstC * self.var.MtoM3C
        self.var.lakeVolumeM3C -= act_bigLakeAbstC * self.var.MtoM3C
        self.var.reservoirStorageM3C -= act_bigLakeAbstC * self.var.MtoM3C
        self.var.lakeResStorageC -= act_bigLakeAbstC * self.var.MtoM3C  # Combined storage for water balance

        # Store abstraction result
        self.var.abstractedLakeReservoirM3C = act_bigLakeAbstC.copy() * self.var.MtoM3C

        # Decompress lake/reservoir storage back to full grid
        self.var.lakeResStorage = globals.inZero.copy()
        np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)

        # Calculate fraction of demand met by large lakes/reservoirs
        bigLakesFactorC = divideValues(act_bigLakeAbstC, remainNeedBigC)

        # Expand fraction back to full grid
        bigLakesFactor = globals.inZero.copy()
        np.put(bigLakesFactor, self.var.decompress_LR, bigLakesFactorC)

        # Aggregate maximum abstraction factor around lake buffer zones
        bigLakesFactorAllaroundlake = npareamaximum(bigLakesFactor, self.var.waterBodyBuffer)

        # Distribute abstraction from large lakes to surrounding cells
        self.var.act_bigLakeResAbst = remainNeed0 * mskWtrBody_unrestricted * bigLakesFactorAllaroundlake

        # Calculate remaining demand not met by large lakes
        remainNeed1 = remainNeed0 * (1 - bigLakesFactorAllaroundlake)

        # Abstraction from small lakes if enabled
        if returnBool('useSmallLakes'):
            # Apply conservative limit to small lake storage
            minlake = np.maximum(0., 0.98 * self.var.smalllakeStorage) * self.var.M3toM
            # Determine actual abstraction from small lakes
            self.var.act_smallLakeResAbst = np.minimum(minlake, remainNeed1)
            # Update small lake storage and volume
            self.var.smalllakeVolumeM3 -= self.var.act_smallLakeResAbst * self.var.MtoM3
        else:
            self.var.act_smallLakeResAbst = 0

        # Update total surface water abstraction to include lakes and reservoirs
        self.var.act_SurfaceWaterAbstract += self.var.act_bigLakeResAbst + self.var.act_smallLakeResAbst

        # Store total lake abstraction separately
        self.var.act_lakeAbst = self.var.act_bigLakeResAbst + self.var.act_smallLakeResAbst


    def process_reservoir_transfers(self):
        """
        Processes water transfers between reservoirs and lakes, including abstraction to/from outside the basin.
        Updates storage volumes and abstraction fractions accordingly.
        """

        for transfer in self.var.reservoir_transfers:
            # Initialize temporary storage change array (compressed format)
            self.var.inZero_C = np.compress(self.var.compress_LR, globals.inZero.copy())

            # Determine the current year based on dynamic or fixed reservoir construction
            year = dateVar['currDate'].year if returnBool('dynamicLakesRes') else loadmap('fixLakesResYear')

            # Check if giver reservoir is constructed
            giver_already_constructed = True
            if transfer[0] > 0:
                index_giver = np.where(self.var.waterBodyID_C == transfer[0])[0][0]
                giver_already_constructed = self.var.resYearC[index_giver] <= year

            # Check if receiver reservoir is constructed
            receiver_already_constructed = True
            if transfer[1] > 0:
                index_receiver = np.where(self.var.waterBodyID_C == transfer[1])[0][0]
                receiver_already_constructed = self.var.resYearC[index_receiver] <= year

            # Proceed only if both reservoirs are constructed
            if receiver_already_constructed and giver_already_constructed:
                reservoir_unused = self.var.resVolumeC - self.var.reservoirStorageM3C

                # Determine unused capacity of receiver
                reservoir_unused_receiver = (
                    reservoir_unused[index_receiver] if transfer[1] > 0 else 10e12
                )

                # Determine available storage from giver
                reservoir_storage_giver = (
                    self.var.reservoirStorageM3C[index_giver] if transfer[0] > 0
                    else self.var.resVolumeC[index_receiver]
                )

                # Calculate actual transfer volume
                reservoir_transfer_actual = np.minimum(
                    reservoir_unused_receiver * 0.95,
                    reservoir_storage_giver * transfer[2] if transfer[2] <= 1 else transfer[2]
                )

                # Apply transfer: subtract from giver, add to receiver
                if transfer[0] > 0:
                    self.var.inZero_C[index_giver] = -reservoir_transfer_actual
                    self.var.reservoir_transfers_out_M3C[index_giver] += reservoir_transfer_actual
                else:
                    self.var.reservoir_transfers_from_outside_M3C[index_receiver] += reservoir_transfer_actual

                if transfer[1] > 0:
                    self.var.inZero_C[index_receiver] = reservoir_transfer_actual
                    self.var.reservoir_transfers_in_M3C[index_receiver] += reservoir_transfer_actual
                else:
                    self.var.reservoir_transfers_to_outside_M3C[index_giver] += reservoir_transfer_actual

                # Update storage volumes
                self.var.lakeStorageC += self.var.inZero_C
                self.var.lakeVolumeM3C += self.var.inZero_C
                self.var.lakeResStorageC += self.var.inZero_C
                self.var.reservoirStorageM3C += self.var.inZero_C
                self.var.reservoir_transfers_net_M3C += self.var.inZero_C

                # If water is transferred out of the basin, adjust demand and abstraction
                if transfer[1] == 0:
                    to_outside_basin = globals.inZero.copy()
                    np.put(to_outside_basin, self.var.decompress_LR, self.var.inZero_C)

                    pot_Lake_Industry -= to_outside_basin * self.var.M3toM
                    self.var.act_lakeAbst -= to_outside_basin * self.var.M3toM
                    self.var.act_SurfaceWaterAbstract -= to_outside_basin * self.var.M3toM
                    self.var.act_bigLakeResAbst -= to_outside_basin * self.var.M3toM
                    self.var.industryDemand -= to_outside_basin * self.var.M3toM
                    self.var.pot_industryConsumption -= to_outside_basin * self.var.M3toM
                    self.var.ind_efficiency = divideValues(
                        self.var.pot_industryConsumption,
                        self.var.industryDemand
                    )

        # Decompress updated transfer arrays back to full grid
        np.put(self.var.reservoir_transfers_net_M3, self.var.decompress_LR, self.var.reservoir_transfers_net_M3C)
        self.var.reservoir_transfers_net_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        np.put(self.var.reservoir_transfers_in_M3, self.var.decompress_LR, self.var.reservoir_transfers_in_M3C)
        self.var.reservoir_transfers_in_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        np.put(self.var.reservoir_transfers_out_M3, self.var.decompress_LR, self.var.reservoir_transfers_out_M3C)
        self.var.reservoir_transfers_out_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        np.put(self.var.reservoir_transfers_from_outside_M3, self.var.decompress_LR, self.var.reservoir_transfers_from_outside_M3C)
        self.var.reservoir_transfers_from_outside_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        np.put(self.var.reservoir_transfers_to_outside_M3, self.var.decompress_LR, self.var.reservoir_transfers_to_outside_M3C)
        self.var.reservoir_transfers_to_outside_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())

        # Adjust abstraction fractions if sector-specific abstraction is enabled
        if self.var.sectorSourceAbstractionFractions:
            self.var.swAbstractionFraction_Res_Industry = np.where(
                self.var.reservoir_transfers_to_outside_M3 > 0,
                0,
                self.var.swAbstractionFraction_Res_Industry
            )
            self.var.gwAbstractionFraction_Industry = np.where(
                self.var.reservoir_transfers_to_outside_M3 > 0,
                0,
                self.var.gwAbstractionFraction_Industry
            )
        else:
            self.var.pot_SurfaceAbstract -= to_outside_basin
            self.var.swAbstractionFraction = np.where(
                self.var.reservoir_transfers_to_outside_M3 != 0,
                1,
                self.var.swAbstractionFraction_nonIrr
            )


    def substract_surface_abstraction(self):
        """
        Calculates the remaining surface water abstraction need after accounting for 
        lake and other water sources.

        Returns:
            remainNeed2 (float or array): The remaining surface water abstraction need.
        """
        # Check if sector-specific abstraction fractions are defined
        if self.var.sectorSourceAbstractionFractions:
            # --- A: Allocate lake abstraction to sectors based on potential and availability ---
            self.var.Lake_Domestic = np.minimum(self.var.act_lakeAbst, pot_Lake_Domestic)
            self.var.Lake_Livestock = np.minimum(self.var.act_lakeAbst - self.var.Lake_Domestic,
                                                pot_Lake_Livestock )
            self.var.Lake_Industry = np.minimum(self.var.act_lakeAbst - self.var.Lake_Domestic - 
                                                self.var.Lake_Livestock, pot_Lake_Industry )       
            self.var.Lake_Irrigation = np.minimum(self.var.act_lakeAbst - self.var.Lake_Domestic - 
                                                self.var.Lake_Livestock - self.var.Lake_Industry,
                                                pot_Lake_Irrigation )

            # --- B: Calculate potential reservoir abstraction for each sector ---
            self.var.pot_Res_Domestic = np.minimum(self.var.swAbstractionFraction_Res_Domestic * self.var.domesticDemand,
                self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic -
                self.var.Lift_Domestic - self.var.Lake_Domestic )

            self.var.pot_Res_Livestock = np.minimum(self.var.swAbstractionFraction_Res_Livestock * self.var.livestockDemand,
                self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock -
                self.var.Lift_Livestock - self.var.Lake_Livestock )

            self.var.pot_Res_Industry = np.minimum(self.var.swAbstractionFraction_Res_Industry * self.var.industryDemand,
                self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry -
                self.var.Lift_Industry - self.var.Lake_Industry )

            self.var.pot_Res_Irrigation = np.minimum(self.var.swAbstractionFraction_Res_Irrigation * self.var.totalIrrDemand,
                self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation -
                self.var.Lift_Irrigation - self.var.Lake_Irrigation )

            # Currently only irrigation is considered in the remaining need
            remainNeed2 = self.var.pot_Res_Irrigation
        else:
            remainNeed2 = self.var.pot_SurfaceAbstract - self.var.act_SurfaceWaterAbstract

        return remainNeed2    

    def allocate_reservoir_abstraction_sectors(self):
        """
        Allocates surface water abstraction from reservoirs to meet domestic, livestock, and industrial demands 
        before irrigation. It calculates the available reservoir storage, determines the maximum allowable 
        abstraction based on release rules, and updates the abstraction volumes accordingly.

        Only used when sector-specific abstraction fractions are defined (sectorSourceAbstractionFractions = True).
        """
        # Calculate total demand for domestic, livestock, and industry (before irrigation)
        remainNeedPre = self.var.pot_Res_Domestic + self.var.pot_Res_Livestock + self.var.pot_Res_Industry

        # Aggregate demand per command area
        demand_Segment = np.where(self.var.reservoir_command_areas > 0,
            npareatotal(remainNeedPre * self.var.cellArea, self.var.reservoir_command_areas), 0)

        # --- Reservoir associated with the Command Area ---
        # If there is more than one reservoir in a command area, the storage of the reservoir with
        # maximum storage in this time-step is chosen. The map resStorageTotal_alloc holds this
        # maximum reservoir storage within a command area in all cells within that command area
        # Filter valid reservoirs (exclude restricted ones)
        ReservoirsThatAreCurrentlyReservoirs = np.where( (self.var.waterBodyTypCTemp == 2) | 
            (self.var.waterBodyTypCTemp == 4), self.var.reservoirStorageM3C, 0)
        ReservoirsThatAreCurrentlyReservoirs = np.where( 
            np.compress(self.var.compress_LR, self.var.resId_restricted) == 0,
            ReservoirsThatAreCurrentlyReservoirs, 0)

        # Map reservoir storage to full grid
        reservoirStorageM3 = globals.inZero.copy()
        np.put(reservoirStorageM3, self.var.decompress_LR, ReservoirsThatAreCurrentlyReservoirs)

        # Determine max storage per command area
        resStorageTotal_alloc = np.where( self.var.reservoir_command_areas > 0,
            npareamaximum(reservoirStorageM3, self.var.reservoir_command_areas), 0)

        # Identify max-storage reservoirs within each command area
        resStorageTotal_allocC = np.compress(self.var.compress_LR, resStorageTotal_alloc)
        resStorageTotal_allocC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C, resStorageTotal_allocC)

        day_of_year = globals.dateVar['currDate'].timetuple().tm_yday

        # Load reservoir release rules
        if 'Reservoir_releases' in binding:
            resStorage_maxFracForIrrigation = readnetcdf2(
                'Reservoir_releases', day_of_year, useDaily='DOY', value='Fraction of Volume'
            )
            resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR, resStorage_maxFracForIrrigation)
        elif self.var.reservoir_releases_excel_option:
            resStorage_maxFracForIrrigation = globals.inZero.copy()
            resStorage_maxFracForIrrigationC = np.where(
                self.var.lakeResStorage_release_ratioC > -1,
                self.var.reservoir_supply[globals.dateVar['doy'] - 1],
                0.03
            )
        else:
            resStorage_maxFracForIrrigation = 0.03 + globals.inZero.copy()
            resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR, resStorage_maxFracForIrrigation)

        # Apply release rules only to max-storage reservoirs
        resStorage_maxFracForIrrigationC = np.multiply(
            resStorageTotal_allocC == self.var.reservoirStorageM3C,
            resStorage_maxFracForIrrigationC )
        np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR, resStorage_maxFracForIrrigationC)

        # Aggregate max release fraction per command area
        resStorage_maxFracForIrrigation_CA = np.where( self.var.reservoir_command_areas > 0,
            npareamaximum(resStorage_maxFracForIrrigation, self.var.reservoir_command_areas), 0 )

        # Calculate actual abstraction from reservoirs
        act_bigLakeResAbst_alloc = np.minimum( resStorage_maxFracForIrrigation_CA * resStorageTotal_alloc,
            demand_Segment / self.var.Water_conveyance_efficiency )

        # Compute abstraction factor (fraction of available storage used)
        ResAbstractFactor = np.where( resStorageTotal_alloc > 0,
            divideValues(act_bigLakeResAbst_alloc, resStorageTotal_alloc), 0 )
        ResAbstractFactorC = np.compress(self.var.compress_LR, ResAbstractFactor)
        ResAbstractFactorC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C, ResAbstractFactorC)

        # Update reservoir and lake storage volumes
        self.var.lakeStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.lakeVolumeM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.lakeResStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.reservoirStorageM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC

        # Track abstracted volumes
        self.var.abstractedLakeReservoirM3C += self.var.reservoirStorageM3C * ResAbstractFactorC
        np.put(self.var.abstractedLakeReservoirM3, self.var.decompress_LR, self.var.abstractedLakeReservoirM3C)

        self.var.lakeResStorage = globals.inZero.copy()
        np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)

        # Determine how much of the demand was met
        metRemainSegment = np.where( demand_Segment > 0,
            divideValues(act_bigLakeResAbst_alloc * self.var.Water_conveyance_efficiency, demand_Segment), 0 )

        # Update actual abstraction values
        self.var.act_bigLakeResAbst += remainNeedPre * metRemainSegment
        self.var.act_SurfaceWaterAbstract += remainNeedPre * metRemainSegment
        self.var.act_ResAbst = remainNeedPre * metRemainSegment

        # Distribute abstraction to sectors
        self.var.Res_Domestic = np.minimum(self.var.act_ResAbst, self.var.pot_Res_Domestic)
        self.var.Res_Livestock = np.minimum(self.var.act_ResAbst - self.var.Res_Domestic, self.var.pot_Res_Livestock)
        self.var.Res_Industry = np.minimum(self.var.act_ResAbst - self.var.Res_Domestic - self.var.Res_Livestock,
            self.var.pot_Res_Industry )


    def allocate_reservoir_abstraction(self, remainNeed2, pot_Res_Irrigation):
        """
        Allocates reservoir water for remaining demand  after surface water abstraction.

        Parameters:
            remainNeed2 (array): Remaining water demand per cell 
            pot_Res_Irrigation (array): Potential reservoir abstraction for irrigation 

        """
        #  Aggregate remaining demand per command area
        demand_Segment = np.where( self.var.reservoir_command_areas > 0,
            npareatotal(remainNeed2 * self.var.cellArea, self.var.reservoir_command_areas), 0)

        # Identify reservoirs (type 2 only) and map their storage to the full grid
        ReservoirsThatAreCurrentlyReservoirs = np.where( self.var.waterBodyTypCTemp == 2,
            self.var.reservoirStorageM3C, 0 )
        reservoirStorageM3 = globals.inZero.copy()
        np.put(reservoirStorageM3, self.var.decompress_LR, ReservoirsThatAreCurrentlyReservoirs)

        # Determine max storage per command area
        resStorageTotal_alloc = np.where( self.var.reservoir_command_areas > 0,
            npareamaximum(reservoirStorageM3, self.var.reservoir_command_areas), 0 )

        # Identify max-storage reservoirs within each command area
        resStorageTotal_allocC = np.compress(self.var.compress_LR, resStorageTotal_alloc)
        resStorageTotal_allocC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C, resStorageTotal_allocC)

        # Determine day of year
        day_of_year = globals.dateVar['currDate'].timetuple().tm_yday

        # Load reservoir release rules
        if 'Reservoir_releases' in binding:
            resStorage_maxFracForIrrigation = readnetcdf2(
                'Reservoir_releases', day_of_year, useDaily='DOY', value='Fraction of Volume' )
        elif self.var.reservoir_releases_excel_option:
            resStorage_maxFracForIrrigation = globals.inZero.copy()
            resStorage_maxFracForIrrigationC = np.where( self.var.lakeResStorage_release_ratioC > -1,
                self.var.reservoir_supply[globals.dateVar['doy'] - 1], 0.03 )
            np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR, resStorage_maxFracForIrrigationC)
        else:
            resStorage_maxFracForIrrigation = 0.03 + globals.inZero.copy()

        # Apply release rules only to max-storage reservoirs
        resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR, resStorage_maxFracForIrrigation)
        resStorage_maxFracForIrrigationC = np.multiply( resStorageTotal_allocC == self.var.reservoirStorageM3C,
            resStorage_maxFracForIrrigationC )
        np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR, resStorage_maxFracForIrrigationC)

        # Aggregate max release fraction per command area
        resStorage_maxFracForIrrigation_CA = np.where( self.var.reservoir_command_areas > 0,
            npareamaximum(resStorage_maxFracForIrrigation, self.var.reservoir_command_areas), 0  )

        # Calculate actual abstraction from reservoirs
        act_bigLakeResAbst_alloc = np.minimum( resStorage_maxFracForIrrigation_CA * resStorageTotal_alloc,
            demand_Segment / self.var.Water_conveyance_efficiency )

        # Compute abstraction factor (fraction of available storage used)
        ResAbstractFactor = np.where( resStorageTotal_alloc > 0,
            divideValues(act_bigLakeResAbst_alloc, resStorageTotal_alloc), 0 )
        ResAbstractFactorC = np.compress(self.var.compress_LR, ResAbstractFactor)
        ResAbstractFactorC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C, ResAbstractFactorC)

        # Update reservoir and lake storage volumes
        self.var.lakeStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.lakeVolumeM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.lakeResStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
        self.var.reservoirStorageM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC

        # Track abstracted volumes
        self.var.abstractedLakeReservoirM3C += self.var.reservoirStorageM3C * ResAbstractFactorC
        np.put(self.var.abstractedLakeReservoirM3, self.var.decompress_LR, self.var.abstractedLakeReservoirM3C)

        self.var.lakeResStorage = globals.inZero.copy()
        np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)

        # Determine how much of the demand was met
        metRemainSegment = np.where( demand_Segment > 0,
            divideValues(act_bigLakeResAbst_alloc * self.var.Water_conveyance_efficiency, demand_Segment), 0 )

        # Calculate leakage losses
        self.var.leakageC_daily = resStorageTotal_allocC * ResAbstractFactorC * (
            1 - np.compress(self.var.compress_LR, self.var.Water_conveyance_efficiency) )
        self.var.leakage = globals.inZero.copy()
        np.put(self.var.leakage, self.var.decompress_LR, self.var.leakageC_daily + self.var.leakage_wwtC_daily)

        divleak_canal = divideValues( self.var.leakageC_daily + self.var.leakage_wwtC_daily,
            self.var.canalsAreaC )
        self.var.leakageCanalsC_M = np.where(self.var.canalsAreaC > 0, divleak_canal, 0)

        self.var.leakageCanals_M = globals.inZero.copy()
        np.put(self.var.leakageCanals_M, self.var.decompress_LR, self.var.leakageCanalsC_M)
        self.var.leakageCanals_M = npareamaximum(self.var.leakageCanals_M, self.var.canals)

        # Update actual abstraction values
        self.var.act_bigLakeResAbst += remainNeed2 * metRemainSegment
        self.var.act_SurfaceWaterAbstract += remainNeed2 * metRemainSegment
        self.var.act_ResAbst += remainNeed2 * metRemainSegment

        # Allocate irrigation abstraction if sector-specific abstraction is enabled
        if self.var.sectorSourceAbstractionFractions:
            self.var.Res_Irrigation = np.minimum(remainNeed2 * metRemainSegment, pot_Res_Irrigation)


    def calculate_groundwater_abstraction(self, totalDemand):
        """
        Calculates potential groundwater abstraction to meet remaining water demand 
        after other sources (desalination, channel, lift, lake, and reservoir) have been used.

        Parameters:
            totalDemand (array or float): Total water demand across all sectors 
        """

        if self.var.sectorSourceAbstractionFractions:
            # Domestic groundwater abstraction potential
            self.var.pot_GW_Domestic = np.minimum( self.var.gwAbstractionFraction_Domestic * self.var.domesticDemand,
                self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic -
                self.var.Lift_Domestic - self.var.Lake_Domestic - self.var.Res_Domestic )

            # Livestock groundwater abstraction potential
            self.var.pot_GW_Livestock = np.minimum( self.var.gwAbstractionFraction_Livestock * self.var.livestockDemand,
                self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock -
                self.var.Lift_Livestock - self.var.Lake_Livestock - self.var.Res_Livestock )

            # Industry groundwater abstraction potential
            self.var.pot_GW_Industry = np.minimum( self.var.gwAbstractionFraction_Industry * self.var.industryDemand,
                self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry -
                self.var.Lift_Industry - self.var.Lake_Industry - self.var.Res_Industry )

            # Irrigation groundwater abstraction potential
            self.var.pot_GW_Irrigation = np.minimum( self.var.gwAbstractionFraction_Irrigation * self.var.totalIrrDemand,
                self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation -
                self.var.Lift_Irrigation - self.var.Lake_Irrigation - self.var.Res_Irrigation )

            # Total potential groundwater abstraction
            self.var.pot_GroundwaterAbstract =  self.var.pot_GW_Domestic + self.var.pot_GW_Livestock + \
                                                self.var.pot_GW_Industry + self.var.pot_GW_Irrigation
        else:
            # If no sector-specific abstraction, use remaining demand after surface water
            self.var.pot_GroundwaterAbstract = totalDemand - self.var.act_SurfaceWaterAbstract

        # Actual abstraction limited by available groundwater storage
        self.var.nonFossilGroundwaterAbs = np.maximum( 0., 
            np.minimum(self.var.readAvlStorGroundwater, self.var.pot_GroundwaterAbstract) )


    def allocate_groundwater_limitAbstraction(self, totalDemand):
        """
        Allocates groundwater abstraction based on remaining demand and available groundwater,
        and finalizes actual water withdrawals for all sectors (domestic, livestock, industry, irrigation).

        Only used when groundwater abstraction limits are enabled (limitAbstraction = True).

        Parameters:
            totalDemand (array): Total water demand across all sectors [m³].

        """
        # Calculate actual surface water abstraction fraction
        act_swAbstractionFraction = divideValues(self.var.act_SurfaceWaterAbstract, totalDemand)

        # Sector-specific abstraction logic
        if self.var.sectorSourceAbstractionFractions:
            # Allocate groundwater to sectors in priority order
            self.var.GW_Domestic = np.minimum(self.var.nonFossilGroundwaterAbs, self.var.pot_GW_Domestic)
            self.var.GW_Livestock = np.minimum(self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic,
                self.var.pot_GW_Livestock )
            self.var.GW_Industry = np.minimum( self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic - 
                self.var.GW_Livestock, self.var.pot_GW_Industry )
            self.var.GW_Irrigation = np.minimum( self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic - 
                self.var.GW_Livestock - self.var.GW_Industry, self.var.pot_GW_Irrigation )
            
            # Compute actual withdrawals
            self.var.act_nonIrrWithdrawal = ( 
                self.var.Desal_Domestic + self.var.Desal_Livestock + self.var.Desal_Industry +
                self.var.Channel_Domestic + self.var.Channel_Livestock + self.var.Channel_Industry +
                self.var.Lift_Domestic + self.var.Lift_Livestock + self.var.Lift_Industry +
                self.var.Lake_Domestic + self.var.Lake_Livestock + self.var.Lake_Industry +
                self.var.Res_Domestic + self.var.Res_Livestock + self.var.Res_Industry +
                self.var.GW_Domestic + self.var.GW_Livestock + self.var.GW_Industry )

            self.var.act_irrWithdrawal = (
                self.var.Desal_Irrigation + self.var.Channel_Irrigation + self.var.Lift_Irrigation +
                self.var.Lake_Irrigation + self.var.Res_Irrigation + self.var.GW_Irrigation )

            # Separate irrigation into surface and groundwater components
            act_irrWithdrawalSW = ( self.var.Desal_Irrigation + self.var.Channel_Irrigation + 
                self.var.Lift_Irrigation + self.var.Lake_Irrigation + self.var.Res_Irrigation )
            act_irrWithdrawalGW = self.var.GW_Irrigation

            # Split irrigation into paddy and non-paddy
            self.var.act_irrNonpaddyWithdrawal = np.minimum( self.var.act_irrWithdrawal,
                self.var.fracVegCover[3] * self.var.irrDemand[3] )
            self.var.act_irrPaddyWithdrawal = self.var.act_irrWithdrawal - self.var.act_irrNonpaddyWithdrawal

        elif self.var.includeIndusDomesDemand:
            # Non-irrigated demand
            act_nonIrrWithdrawalGW = self.var.nonIrrDemand * (1 - act_swAbstractionFraction)
            act_nonIrrWithdrawalGW = np.minimum(self.var.nonFossilGroundwaterAbs, act_nonIrrWithdrawalGW)
            act_nonIrrWithdrawalSW = act_swAbstractionFraction * self.var.nonIrrDemand
            self.var.act_nonIrrWithdrawal = act_nonIrrWithdrawalSW + act_nonIrrWithdrawalGW

            # Irrigated demand
            act_irrWithdrawalGW = self.var.totalIrrDemand * (1 - act_swAbstractionFraction)
            act_irrWithdrawalGW = np.minimum(self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW,
                act_irrWithdrawalGW )
            act_irrWithdrawalSW = act_swAbstractionFraction * self.var.totalIrrDemand
            self.var.act_irrWithdrawal = act_irrWithdrawalSW + act_irrWithdrawalGW

            # Non-paddy
            act_irrnonpaddyGW = self.var.fracVegCover[3] * (1 - act_swAbstractionFraction) * self.var.irrDemand[3]
            act_irrnonpaddyGW = np.minimum( self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW,
                act_irrnonpaddyGW )
            act_irrnonpaddySW = self.var.fracVegCover[3] * act_swAbstractionFraction * self.var.irrDemand[3]
            self.var.act_irrNonpaddyWithdrawal = act_irrnonpaddySW + act_irrnonpaddyGW

            # Paddy
            act_irrpaddyGW = self.var.fracVegCover[2] * (1 - act_swAbstractionFraction) * self.var.irrDemand[2]
            act_irrpaddyGW = np.minimum( self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW - act_irrnonpaddyGW,
                act_irrpaddyGW )
            act_irrpaddySW = self.var.fracVegCover[2] * act_swAbstractionFraction * self.var.irrDemand[2]
            self.var.act_irrPaddyWithdrawal = act_irrpaddySW + act_irrpaddyGW

        else:
            # Only irrigation is considered
            self.var.act_nonIrrWithdrawal = globals.inZero.copy()

            act_irrWithdrawalGW = self.var.totalIrrDemand * (1 - act_swAbstractionFraction)
            act_irrWithdrawalGW = np.minimum(self.var.nonFossilGroundwaterAbs, act_irrWithdrawalGW)
            act_irrWithdrawalSW = act_swAbstractionFraction * self.var.totalIrrDemand
            self.var.act_irrWithdrawal = act_irrWithdrawalSW + act_irrWithdrawalGW

            # Non-paddy
            act.var.irrDemand[3]
            act_irrnonpaddyGW = np.minimum(self.var.nonFossilGroundwaterAbs, act_irrnonpaddyGW)
            act_irrnonpaddySW = self.var.fracVegCover[3] * act_swAbstractionFraction * self.var.irrDemand[3]
            self.var.act_irrNonpaddyWithdrawal = act_irrnonpaddySW + act_irrnonpaddyGW

            # Paddy
            act_irrpaddyGW = self.var.fracVegCover[2] * (1 - act_swAbstractionFraction) * self.var.irrDemand[2]
            act.nonFossilGroundwaterAbs - act_irrnonpaddyGW,
                act_irrpaddyGW
            )
            act_irrpaddySW = self.var.fracVegCover[2] * act_swAbstractionFraction * self.var.irrDemand[2]
            self.var.act_irrPaddyWithdrawal = act_irrpaddySW + act_irrpaddyGW

        # Correct irrigation demand to avoid double-counting unmet demand from previous days
        self.var.act_irrPaddyDemand = np.maximum(0, self.var.irrPaddyDemand - self.var.unmetDemandPaddy)
        self.var.act_irrNonpaddyDemand = np.maximum(0, self.var.irrNonpaddyDemand - self.var.unmetDemandNonpaddy)

        # Calculate unmet demand
        if self.var.includeIndusDomesDemand:
            self.var.unmetDemand = ( (self.var.totalIrrDemand - self.var.act_irrWithdrawal) +
                (self.var.nonIrrDemand - self.var.act_nonIrrWithdrawal) )
        else:
            self.var.unmetDemand = (
                self.var.totalIrrDemand - self.var.act_irrWithdrawal - self.var.act_nonIrrWithdrawal )

        self.var.unmetDemandPaddy = self.var.irrPaddyDemand - self.var.act_irrPaddyDemand
        self.var.unmetDemandNonpaddy = self.var.irrNonpaddyDemand - self.var.act_irrNonpaddyDemand


    def allocate_groundwater_abstraction(self):
        """
        If zonal abstraction is not enabled, it assumes full fossil groundwater abstraction is allowed
        and sets sectoral abstractions equal to their potential values.
        """
        if 'zonal_abstraction' in option and checkOption('zonal_abstraction'):
            # Surface water available in channels
            left_sf = self.var.readAvlChannelStorageM

            if self.var.sectorSourceAbstractionFractions:
                # Calculate unmet demand per sector
                unmetChannel_Domestic = self.var.pot_Channel_Domestic - self.var.Channel_Domestic
                unmetChannel_Livestock = self.var.pot_Channel_Livestock - self.var.Channel_Livestock
                unmetChannel_Industry = self.var.pot_Channel_Industry - self.var.Channel_Industry
                unmetChannel_Irrigation = self.var.pot_Channel_Irrigation - self.var.Channel_Irrigation

                unmet_Domestic = self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic - \
                                self.var.Lift_Domestic - self.var.Lake_Domestic - self.var.Res_Domestic - self.var.GW_Domestic
                unmet_Livestock = self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock - \
                                self.var.Lift_Livestock - self.var.Lake_Livestock - self.var.Res_Livestock - self.var.GW_Livestock
                unmet_Industry = self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry - \
                                self.var.Lift_Industry - self.var.Lake_Industry - self.var.Res_Industry - self.var.GW_Industry
                unmet_Irrigation = self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation - \
                                self.var.Lift_Industry - self.var.Lake_Irrigation - self.var.Res_Irrigation - self.var.GW_Irrigation

                # Potential channel abstraction per sector
                self.var.pot_Channel_Domestic = np.minimum(unmetChannel_Domestic, unmet_Domestic)
                self.var.pot_Channel_Livestock = np.minimum(unmetChannel_Livestock, unmet_Livestock)
                self.var.pot_Channel_Industry = np.minimum(unmetChannel_Industry, unmet_Industry)
                self.var.pot_Channel_Irrigation = np.minimum(unmetChannel_Irrigation, unmet_Irrigation)

                unmet_Channel = self.var.pot_Channel_Domestic + self.var.pot_Channel_Livestock + 
                    self.var.pot_Channel_Industry + self.var.pot_Channel_Irrigation
                zoneDemand = npareatotal(unmet_Channel * self.var.cellArea, self.var.allocation_zone)
            else:
                zoneDemand = npareatotal(self.var.unmetDemand * self.var.cellArea, self.var.allocation_zone)

            # Surface water availability per zone
            zone_sf_avail = npareatotal(left_sf, self.var.allocation_zone)
            zone_sf_abstraction = np.minimum(zoneDemand, zone_sf_avail)

            # Distribute surface water abstraction to cells
            cell_sf_abstraction = np.maximum(0., divideValues(left_sf, zone_sf_avail) * zone_sf_abstraction)
            cell_sf_allocation = np.maximum(0., divideValues(self.var.unmetDemand, zoneDemand) * zone_sf_abstraction)

            self.var.act_SurfaceWaterAbstract += cell_sf_abstraction
            self.var.act_channelAbst += cell_sf_abstraction

            if self.var.sectorSourceAbstractionFractions:
                self.var.Channel_Domestic_fromZone = np.minimum(cell_sf_abstraction, self.var.pot_Channel_Domestic)
                self.var.Channel_Domestic += self.var.Channel_Domestic_fromZone

                self.var.Channel_Livestock_fromZone = np.minimum(
                    cell_sf_abstraction - self.var.Channel_Domestic_fromZone,
                    self.var.pot_Channel_Livestock )
                self.var.Channel_Livestock += self.var.Channel_Livestock_fromZone

                self.var.Channel_Industry_fromZone = np.minimum(
                    cell_sf_abstraction - self.var.Channel_Domestic_fromZone - self.var.Channel_Livestock_fromZone,
                    self.var.pot_Channel_Industry )
                self.var.Channel_Industry += self.var.Channel_Industry_fromZone

                self.var.Channel_Irrigation_fromZone = np.minimum(
                    cell_sf_abstraction - self.var.Channel_Domestic_fromZone - self.var.Channel_Livestock_fromZone - self.var.Channel_Industry_fromZone,
                    self.var.pot_Channel_Irrigation )
                self.var.Channel_Irrigation += self.var.Channel_Irrigation_fromZone

            # Update potential groundwater abstraction
            self.var.pot_GroundwaterAbstract = np.maximum(0., self.var.pot_GroundwaterAbstract - cell_sf_allocation)

            # Groundwater allocation
            left_gw_demand = np.maximum(0., self.var.pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs)
            left_gw_avail = self.var.readAvlStorGroundwater - self.var.nonFossilGroundwaterAbs
            zone_gw_avail = npareatotal(left_gw_avail * self.var.cellArea, self.var.allocation_zone)
            zone_gw_demand = zoneDemand - zone_sf_abstraction
            zone_gw_abstraction = np.minimum(zone_gw_demand, zone_gw_avail)

            # Distribute groundwater abstraction to cells
            cell_gw_abstraction = np.maximum(0., divideValues(left_gw_avail, zone_gw_avail) * zone_gw_abstraction)
            cell_gw_allocation = np.maximum(0., divideValues(left_gw_demand, zone_gw_demand) * zone_gw_abstraction)

            self.var.unmetDemand = np.maximum(0., left_gw_demand - cell_gw_allocation)
            self.var.nonFossilGroundwaterAbs += cell_gw_abstraction

            if self.var.sectorSourceAbstractionFractions:
                self.var.GW_Domestic_fromZone = np.minimum(self.var.nonFossilGroundwaterAbs, self.var.pot_GW_Domestic)
                self.var.GW_Domestic += self.var.GW_Domestic_fromZone.copy()

                self.var.GW_Livestock_fromZone = np.minimum(
                    self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic_fromZone,
                    self.var.pot_GW_Livestock )
                self.var.GW_Livestock += self.var.GW_Livestock_fromZone.copy()

                self.var.GW_Industry_fromZone = np.minimum(
                    self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic_fromZone - self.var.GW_Livestock_fromZone,
                    self.var.pot_GW_Industry )
                self.var.GW_Industry += self.var.GW_Industry_fromZone.copy()

                self.var.GW_Irrigation_fromZone = np.minimum(
                    self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic_fromZone - self.var.GW_Livestock_fromZone - self.var.GW_Industry_fromZone,
                    self.var.pot_GW_Irrigation )
                self.var.GW_Irrigation += self.var.GW_Irrigation_fromZone.copy()

        else:
            # No zonal abstraction: assume full fossil groundwater abstraction is allowed
            self.var.unmetDemand = self.var.pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs
            if self.var.sectorSourceAbstractionFractions:
                self.var.GW_Domestic = self.var.pot_GW_Domestic
                self.var.GW_Industry = self.var.pot_GW_Industry
                self.var.GW_Livestock = self.var.pot_GW_Livestock
                self.var.GW_Irrigation = self.var.pot_GW_Irrigation

        # Finalize withdrawals
        if self.var.includeIndusDomesDemand:
            self.var.act_nonIrrWithdrawal = np.copy(self.var.nonIrrDemand)
        self.var.act_irrWithdrawal = np.copy(self.var.totalIrrDemand)

        self.var.act_irrNonpaddyWithdrawal = self.var.fracVegCover[3] * self.var.irrDemand[3]
        self.var.act_irrPaddyWithdrawal = self.var.fracVegCover[2] * self.var.irrDemand[2]

        self.var.act_irrNonpaddyDemand = self.var.act_irrNonpaddyWithdrawal.copy()
        self.var.act_irrPaddyDemand = self.var.act_irrPaddyWithdrawal.copy()

  

def calculate_water_consumption(self, frac_domestic, frac_industry, frac_livestock):
    """
    Calculates actual water consumption and withdrawals for irrigation and non-irrigation sectors,
    based on actual withdrawals, efficiencies, and sectoral abstraction settings.

    Parameters:
        frac_domestic (float): Fraction of non-irrigation demand attributed to domestic use.
        frac_industry (float): Fraction of non-irrigation demand attributed to industry.
        frac_livestock (float): Fraction of non-irrigation demand attributed to livestock.
    """
    # Irrigation consumption
    # paddy
    self.var.act_irrConsumption[2] = divideValues(
        self.var.act_irrPaddyWithdrawal, self.var.fracVegCover[2] ) * self.var.efficiencyPaddy
    # non-paddy
    self.var.act_irrConsumption[3] = divideValues(
        self.var.act_irrNonpaddyWithdrawal, self.var.fracVegCover[3] ) * self.var.efficiencyNonpaddy

    if self.var.sectorSourceAbstractionFractions:
        # Non-irrigation withdrawals and consumption by sector 
        self.var.act_domWithdrawal = self.var.Channel_Domestic + self.var.Lift_Domestic + \
                                        self.var.Desal_Domestic + self.var.Lake_Domestic + \
                                        self.var.Res_Domestic + self.var.GW_Domestic
        self.var.act_livWithdrawal = self.var.Channel_Livestock + self.var.Lift_Livestock + \
                                        self.var.Desal_Livestock + self.var.Lake_Livestock + \
                                        self.var.Res_Livestock + self.var.GW_Livestock
        self.var.act_indWithdrawal = self.var.Channel_Industry + self.var.Lift_Industry + \
                                        self.var.Desal_Industry + self.var.Lake_Industry + \
                                        self.var.Res_Industry + self.var.GW_Industry

        self.var.act_indConsumption = self.var.ind_efficiency * self.var.act_indWithdrawal
        self.var.act_domConsumption = self.var.dom_efficiency * self.var.act_domWithdrawal
        self.var.act_livConsumption = self.var.liv_efficiency * self.var.act_livWithdrawal
        self.var.act_nonIrrConsumption = self.var.act_domConsumption + self.var.act_indConsumption + \
                                            self.var.act_livConsumption

    elif self.var.includeIndusDomesDemand:
        # Split non-irrigation withdrawal by fractions
        self.var.act_indWithdrawal = frac_industry * self.var.act_nonIrrWithdrawal
        self.var.act_domWithdrawal = frac_domestic * self.var.act_nonIrrWithdrawal
        self.var.act_livWithdrawal = frac_livestock * self.var.act_nonIrrWithdrawal

        # Consumption based on efficiency
        self.var.act_indConsumption = self.var.ind_efficiency * self.var.act_indWithdrawal
        self.var.act_domConsumption = self.var.dom_efficiency * self.var.act_domWithdrawal
        self.var.act_livConsumption = self.var.liv_efficiency * self.var.act_livWithdrawal

        self.var.act_nonIrrConsumption = (
            self.var.act_domConsumption + self.var.act_indConsumption + self.var.act_livConsumptionConsumption = globals.inZero.copy()
            
    else:  # only irrigation is considered
        self.var.act_nonIrrConsumption = globals.inZero.copy()

    # --- Total irrigation consumption ---
    self.var.act_totalIrrConsumption = ( self.var.fracVegCover[2] * self.var.act_irrConsumption[2] +
        self.var.fracVegCover[3] * self.var.act_irrConsumption[3] )
    self.var.act_paddyConsumption = self.var.fracVegCover[2] * self.var.act_irrConsumption[2]
    self.var.act_nonpaddyConsumption = self.var.fracVegCover[3] * self.var.act_irrConsumption[3]

    # --- Total water demand, withdrawal, and consumption ---
    if self.var.includeIndusDomesDemand:
        self.var.totalWaterDemand = (
            self.var.fracVegCover[2] * self.var.irrDemand[2] +
            self.var.fracVegCover[3] * self.var.irrDemand[3] +
            self.var.nonIrrDemand
        )
        self.var.act_totalWaterWithdrawal = (
            self.var.act_nonIrrWithdrawal + self.var.act_irrWithdrawal
        )
        self.var.act_totalWaterConsumption = (
            self.var.act_nonIrrConsumption + self.var.act_totalIrrConsumption
        )
    else:
        self.var.totalWaterDemand = (
            self.var.fracVegCover[2] * self.var.irrDemand[2] +
            self.var.fracVegCover[3] * self.var.irrDemand[3]
        )
        self.var.act_totalWaterWithdrawal = np.copy(self.var.act_irrWithdrawal)
        self.var.act_totalWaterConsumption = np.copy(self.var.act_totalIrrConsumption)


    def calculate_return_flow(self):        
        """
        Calculates return flows from irrigation and non-irrigation sectors, as well as evaporation losses,
        based on actual withdrawals, consumption, return fractions, and unmet demand. 
        """
        # --- calculate return flow
        # Sum up loss - difference between withdrawn and consumed - split into return flow and evaporation
        sumIrrLoss = self.var.act_irrWithdrawal - self.var.act_totalIrrConsumption

        self.var.returnflowIrr = self.var.returnfractionIrr * sumIrrLoss
        self.var.addtoevapotrans = (1 - self.var.returnfractionIrr) * sumIrrLoss

        if self.var.sectorSourceAbstractionFractions:
            self.var.returnflowNonIrr = self.var.act_nonIrrWithdrawal - self.var.act_nonIrrConsumption
        elif self.var.includeIndusDomesDemand:  # all demands are taken into account
            self.var.returnflowNonIrr = self.var.nonIrrReturnFlowFraction * self.var.act_nonIrrWithdrawal

        # limit return flow to not put all fossil groundwater back into the system, because it can lead to higher
        # river discharge than without water demand, as water is taken from fossil groundwater (out of system)
        unmet_div_ww = 1. - np.minimum(1, divideValues(self.var.unmetDemand,
                                                        self.var.act_totalWaterWithdrawal + self.var.unmetDemand))

        # 'fossil_water_treated_normally' means that there is no lost fossil water
        if 'fossil_water_treated_normally' in option:
            if checkOption('fossil_water_treated_normally'):
                unmet_div_ww = 1

        if checkOption('limitAbstraction'):
            unmet_div_ww = 1

        if self.var.includeIndusDomesDemand:  # all demands are taken into account
            self.var.unmet_lost = (self.var.returnflowIrr + self.var.returnflowNonIrr + self.var.addtoevapotrans) \
                                    * (1 - unmet_div_ww)
        else:  # only irrigation is considered
            self.var.unmet_lost = (self.var.returnflowIrr + self.var.addtoevapotrans) * (1 - unmet_div_ww)

        # self.var.waterDemandLost = self.var.act_totalWaterConsumption + self.var.addtoevapotrans
        self.var.unmet_lostirr = (self.var.returnflowIrr + self.var.addtoevapotrans) * (1 - unmet_div_ww)
        if self.var.includeIndusDomesDemand:  # all demands are taken into account
            self.var.unmet_lostNonirr = self.var.returnflowNonIrr * (1 - unmet_div_ww)

        self.var.returnflowIrr = self.var.returnflowIrr * unmet_div_ww
        self.var.addtoevapotrans = self.var.addtoevapotrans * unmet_div_ww
        if self.var.includeIndusDomesDemand:  # all demands are taken into account
            self.var.returnflowNonIrr = self.var.returnflowNonIrr * unmet_div_ww            

        # returnflow to river and to evapotranspiration
        if self.var.includeIndusDomesDemand:  # all demands are taken into account
            self.var.returnFlow = self.var.returnflowIrr + self.var.returnflowNonIrr
        else:  # only irrigation is considered
            self.var.returnFlow = self.var.returnflowIrr

        self.var.waterabstraction = self.var.nonFossilGroundwaterAbs + self.var.unmetDemand + \
                                    self.var.act_SurfaceWaterAbstract