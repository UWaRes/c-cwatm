# -------------------------------------------------------------------------
# Name:       C-CWATM Model Dynamic
# Purpose:
#
# Author:      Peter Greve, Peter Burek
#
# Created:     30/11/2023
# Copyright:   (c) Peter Burek 2016, (c) Peter Greve 2023
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
from cwatm.management_modules.messages import *

import time


class CWATModel_dyn(DynamicModel):

    # =========== DYNAMIC ====================================================

    def dynamic(self):
        """
        Dynamic part of C-CWATM
        calls the dynamic part of the hydrological modules
        Looping through time and space

        Note:
            if flags set the output on the screen can be changed e.g.

            * v: no output at all
            * l: time and first gauge discharge
            * t: timing of different processes at the end
        """

        timestep_dynamic(self)

        del timeMes[:]
        timemeasure("Start dynamic")

        # ************************************************************
        """ up to here it was fun, now the real stuff starts
        """

        if checkOption('calc_environflow') and (returnBool('calc_ef_afterRun') == False):
            # if only the dis is used for calculation of EF
            self.environflow_module.dynamic()
            self.output_module.dynamic(ef=True)
            header = "\n\n ======================== CWATM ONLY EF calculation===========\n"
            print(header + "done with Environmental Flow\n")
            sys.exit(400)

        # ***** Read forcing ****************
        self.readmeteo_module.dynamic()
        timemeasure("Read meteo")  # 1. timing after read input maps

        if Flags['calib']:
            self.output_module.dynamic()
            return


        """ Here it starts with hydrological modules:
        """

        # ***** INFLOW HYDROGRAPHS (OPTIONAL)****************
        self.inflow_module.dynamic()
        
        # ***** Lakes and Reservoirs ****************
        self.lakes_reservoirs_module.dynamic()

        # ***** READ land use fraction maps***************************
        self.landcoverType_module.dynamic_fracIrrigation(init=dateVar['newYear'], dynamic=self.var.dynamicLandcover)

        # ***** Water demand ****************
        self.waterdemand_module.dynamic()

        # ***** Groundwater ****************
        self.groundwater_module.dynamic()
        timemeasure("Groundwater")  # 2. timing

        # ***** Small Lakes and Reservoirs****************
        self.lakes_res_small_module.dynamic()
        timemeasure("Small lakes")  # 3. timing

        # ***** River Routing****************
        self.routing_kinematic_module.dynamic()
        timemeasure("Routing_Kin")  # 4. timing
        

        # ------------------------------------------------------
        # End of calculation -----------------------------------
        # ------------------------------------------------------

        self.environflow_module.dynamic()
        # in case environmental flow is calculated last

        self.output_module.dynamic()
        timemeasure("Output")  # 5. timing

        self.init_module.dynamic()

        for i in range(len(timeMes)):
            if self.currentStep == self.firstStep:
                timeMesSum.append(timeMes[i] - timeMes[0])
            else:
                timeMesSum[i] += timeMes[i] - timeMes[0]


