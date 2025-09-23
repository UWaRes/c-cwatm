# -------------------------------------------------------------------------
# Name:        Groundwater module
# Purpose:
#
# Author:      PB
#
# Created:     15/07/2016
# Copyright:   (c) PB 2016
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *


class groundwater(object):
    """
    GROUNDWATER


    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    load_initial                           Settings initLoad holds initial conditions for variables                input
    storGroundwater                        Groundwater storage (non-fossil). This is primarily used when not usin  m    
    specificYield                          groundwater reservoir parameters used to comp                           m    
    recessionCoeff                         groundwater storage times this coefficient gives baseflow               frac 
    readAvlStorGroundwater                 same as storGroundwater but equal to 0 when inferior to a treshold      m    
    prestorGroundwater                     storGroundwater at the beginning of each step                           m    
    sum_gwRecharge                         groundwater recharge                                                    m    
    baseflow                               simulated baseflow (= groundwater discharge to river)                   m     
    nonFossilGroundwaterAbs                Non-fossil groundwater abstraction.                                     m    
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model
        
    def initial(self):
        """
        Initial part of the groundwater module

        * load parameters from settings file
        * initial groundwater storage
        """

        self.var.recessionCoeff = loadmap('recessionCoeff')

        # for CALIBRATION
        self.var.recessionCoeff = 1 / self.var.recessionCoeff * loadmap('recessionCoeff_factor')
        self.var.recessionCoeff = 1 / self.var.recessionCoeff

        self.var.specificYield = loadmap('specificYield')

        # init calculation recession coefficient, speciefic yield, ksatAquifer
        self.var.recessionCoeff = np.maximum(5.e-4, self.var.recessionCoeff)
        self.var.recessionCoeff = np.minimum(1.000, self.var.recessionCoeff)
        self.var.specificYield = np.maximum(0.010, self.var.specificYield)
        self.var.specificYield = np.minimum(1.000, self.var.specificYield)

        # initial conditions
        self.var.storGroundwater = self.var.load_initial('storGroundwater')
        if 'storGroundwater' in binding and not self.var.loadInit:
            self.var.storGroundwater = loadmap('storGroundwater')
        self.var.storGroundwater = np.maximum(0.0, self.var.storGroundwater) + globals.inZero

        # for water demand to have some initial value
        tresholdStorGroundwater = 0.00001  # 0.01 mm
        self.var.readAvlStorGroundwater = np.where(self.var.storGroundwater > tresholdStorGroundwater,
                                                   self.var.storGroundwater - tresholdStorGroundwater, 0.0)


# --------------------------------------------------------------------------

    def dynamic(self):
        """
        Dynamic part of the groundwater module
        Calculate groundwater storage and baseflow
        """

        # WATER DEMAND
        # update storGoundwater after self.var.nonFossilGroundwaterAbs
        self.var.storGroundwater = np.maximum(0., self.var.storGroundwater - self.var.nonFossilGroundwaterAbs)
        # PS: We assume only local groundwater abstraction can happen (only to satisfy water demand within a cell).
        # unmetDemand (m), satisfied by fossil gwAbstractions (and/or desalinization or other sources)
        # (equal to zero if limitAbstraction = True)

        # get net recharge (percolation-capRise) and update storage:
        self.var.storGroundwater = np.maximum(0., self.var.storGroundwater + self.var.sum_gwRecharge)

        # Calculate groundwater baseflow with linear storage function
        self.var.baseflow = np.maximum(0., np.minimum(self.var.storGroundwater, self.var.recessionCoeff * self.var.storGroundwater))
        # update groundwater storage
        self.var.storGroundwater = np.maximum(0., self.var.storGroundwater - self.var.baseflow)

        # to avoid small values and to avoid excessive abstractions from dry groundwater
        tresholdStorGroundwater = 0.00001  # 0.01 mm
        self.var.readAvlStorGroundwater = np.where(self.var.storGroundwater > tresholdStorGroundwater, self.var.storGroundwater - tresholdStorGroundwater,0.0)





