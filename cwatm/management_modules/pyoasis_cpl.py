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

def def_local_partition(nlon, nlat, npes, mype):
    il_extentx = nlon
    il_extenty = nlat//npes
    if (mype == npes-1):
        il_extenty = nlat - nlat//npes * mype
    il_size = il_extentx * il_extenty
    il_offsetx = 0
    il_offsety = nlat//npes * mype
    il_offset = nlon * il_offsety
    return il_extentx, il_extenty, il_size, il_offsetx, il_offsety, il_offset

def def_paral(decomp, il_offset, il_size, il_extentx, il_extenty, nlon):
    if decomp == 'DECOMP_APPLE':
        return [1, il_offset, il_size]
    elif decomp == 'DECOMP_BOX':
        return [2, il_offset, il_extentx, il_extenty, nlon]
    else:
        raise Exception(f'Invalid CPPKEYDECOMP in {__file__}: {decomp}')


def oasis_specify_partition(comp,rlat,rlon):
        ### OASIS_GET_LOCALCOMM ###
        local_comm = comp.localcomm
        # Get rank in local communicator
        mype = comp.localcomm.rank
        npes = comp.localcomm.size
        # Unit for output messages : one file for each process
        w_unit = 100 + mype
        comp_out_dummy = comp.name + '.out_'+str(w_unit)
        w_unit = open(comp_out_dummy, 'w')
        print(' -----------------------------------------------------------', file=w_unit)
        print(f' I am atmos process with rank : {mype}', file=w_unit)
        print(f' in my local communicator gathering {npes} processes', file=w_unit)
        print(' ----------------------------------------------------------', file=w_unit)
        w_unit.flush()

        ###  PARTITION DEFINITION ###
        il_extentx, il_extenty, il_size, il_offsetx, il_offsety, il_offset = def_local_partition(rlon, rlat, npes, mype)
        print(' Local partition definition', file=w_unit)
        print(' il_extentx, il_extenty, il_size, il_offsetx, il_offsety, il_offset = ', il_extentx, il_extenty, il_size, il_offsetx, il_offsety, il_offset, file=w_unit)
        w_unit.flush()

        ### OASIS_DEF_PARTITION ###
        CPPKEYDECOMP = 'DECOMP_APPLE'
        ig_paral = def_paral(CPPKEYDECOMP, il_offset, il_size, il_extentx, il_extenty, rlon)
        partition = pyoasis.ApplePartition(il_offset, il_size)

        return w_unit

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
        comp = pyoasis.Component("hydro_component")
        #self.comp = pyoasis.Component("hydro_component")
        # initialize pyoasis component, see toy model -ok
        # 1) get localcomm -ok
        # 2) partition definition
        # 3) grid definition
        # 4) declaration of coupling fields

        ### get localcomm and define partition
        w_unit = oasis_specify_partition(comp,rlat=437,rlon=625)

        ### TERMINATION OF DEFINITION PHASE ###
        print(' End of initialisation phase', file=w_unit)
        w_unit.flush()

        ### OASIS_ENDDEF ###
        comp.enddef()

    def dynamic(self):
        # 1) get
        # 2) put
        pass




