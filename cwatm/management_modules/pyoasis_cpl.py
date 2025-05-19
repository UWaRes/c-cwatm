# -------------------------------------------------------------------------
# Name:        pyoasis module
# Purpose:
#
# Author:      Amelie Schmitt
#
# Created:     25/02/2025
# Copyright:   (c) Amelie Schmitt 2025 
# -------------------------------------------------------------------------

# import stuff for pyoasis
import pyoasis
#from pyoasis import OASIS 


class pyoasis_cpl(object):

    """
    OASIS3-MCT_5.0 coupler

    functions related to the pyOASIS coupler

    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====

    dummydummy                                                                                                     --   
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

    def initial(self):
        self.comp = pyoasis.Component("hydro_component")
        # initialize pyoasis component, see toy model
        # 1) get localcomm
        # 2) partition definition
        # 3) grid definition
        # 4) declaration of coupling fields

        ### OASIS_ENDDEF ###
        self.comp.enddef()

    def dynamic(self):
        # 1) get
        # 2) put
        pass




