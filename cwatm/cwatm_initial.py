# -------------------------------------------------------------------------
# Name:       C-CWATM Initial
# Purpose:
#
# Author:      Peter Greve, Peter Burek
#
# Created:     30/11/2023
# Copyright:   (c) Peter Greve, Peter Burek 2023
# -------------------------------------------------------------------------

from cwatm.hydrological_modules.miscInitial import miscInitial
from cwatm.hydrological_modules.initcondition import initcondition

from cwatm.hydrological_modules.readmeteo import readmeteo
from cwatm.hydrological_modules.inflow import inflow
from cwatm.hydrological_modules.landcoverType import landcoverType
from cwatm.hydrological_modules.groundwater import groundwater
from cwatm.hydrological_modules.water_demand.water_demand import water_demand
from cwatm.hydrological_modules.water_demand.wastewater import waterdemand_wastewater as wastewater
from cwatm.hydrological_modules.lakes_res_small import lakes_res_small
from cwatm.hydrological_modules.environflow import environflow
from cwatm.hydrological_modules.routing_reservoirs.routing_kinematic import routing_kinematic
from cwatm.hydrological_modules.lakes_reservoirs import lakes_reservoirs

from cwatm.management_modules.output import *
from cwatm.management_modules.data_handling import *
from cwatm.management_modules.coupling import *
import os, glob


class Variables:
    def load_initial(self, name, default=0.0, number=None):
        """
        First it is checked if the initial value is given in the settings file

        * if it is <> None it is used directly
        * if None it is loaded from the init netcdf file

        :param name: Name of the init value
        :param default: default value -> default is 0.0
        :return: spatial map or value of initial condition
        """

        if number is not None:
            name = name + str(number)

        if self.loadInit:
            map = readnetcdfInitial(self.initLoadFile, name)
            if Flags['calib']:
                self.initmap[name] = map
            return map
        else:
            return default

class Config:
    pass


class CWATModel_ini(DynamicModel):

    """
    C-CWATN initial part
    this part is to initialize the variables.
    It will call the initial part of the hydrological modules
    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self):
        """
        Init part of the initial part
        defines the mask map and the outlet points
        initialization of the hydrological modules
        """

        DynamicModel.__init__(self)

        self.var = Variables()
        self.conf = Config()

        # ----------------------------------------
        # include output of tss and maps
        self.output_module = outputTssMap(self)

        # include all the hydrological modules
        self.misc_module = miscInitial(self)
        self.init_module = initcondition(self)
        self.readmeteo_module = readmeteo(self)
        self.environflow_module = environflow(self)
        self.inflow_module = inflow(self)
        self.landcoverType_module = landcoverType(self)
        self.groundwater_module = groundwater(self)
        self.waterdemand_module = water_demand(self)
        self.wastewater_module = wastewater(self)
        self.lakes_res_small_module = lakes_res_small(self)
        self.routing_kinematic_module = routing_kinematic(self)
        self.lakes_reservoirs_module = lakes_reservoirs(self)

        # as: OASIS3-MCT coupler
        #if binding['coupl_flag']=='full_coupl':
        if binding['coupl_flag']=='no_coupl': # for testing
            from cwatm.management_modules.pyoasis_cpl import pyoasis_cpl
            self.pyoasis_cpl_module = pyoasis_cpl(self)


        # ----------------------------------------
        # reading of the metainformation of variables to put into output netcdfs
        metaNetCDF()

        # no ModFlow
        self.var.modflow = False

        ## MakMap: the maskmap is flexible e.g. col,row,x1,y1  or x1,x2,y1,y2
        # set the maskmap
        self.MaskMap = loadsetclone(self, 'MaskMap')
        
        # run intial misc to get all global variables
        self.misc_module.initial()
        self.init_module.initial()

        # as: initialize oasis coupling        
        # initialize oasis first, then check in readmeteo
        if binding['coupl_flag']=='oasis_coupl': 
            self.pyoasis_cpl_module.initial()
        
        self.readmeteo_module.initial()
        self.inflow_module.initial()
        self.groundwater_module.initial()
        self.landcoverType_module.initial()
        self.lakes_res_small_module.initial()
        self.routing_kinematic_module.initial()
        if checkOption('includeWaterBodies'):
            self.lakes_reservoirs_module.initWaterbodies()
            self.lakes_reservoirs_module.initial_lakes()
            self.lakes_reservoirs_module.initial_reservoirs()

        self.waterdemand_module.initial()
        self.output_module.initial()
        self.environflow_module.initial()

        # TODO: delete after testing
        # as: initialize oasis coupling        
        #if binding['coupl_flag']=='full_coupl':
        if binding['coupl_flag']=='no_coupl': # for testing
            self.pyoasis_cpl_module.initial()    



