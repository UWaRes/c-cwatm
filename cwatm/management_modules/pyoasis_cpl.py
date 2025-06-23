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

#def def_local_partition(nlon, nlat, lcomm_size, lcomm_rank):
#    il_extentx = nlon
#    il_extenty = nlat//lcomm_size
#    if (lcomm_rank == lcomm_size-1):
#        il_extenty = nlat - nlat//lcomm_size * lcomm_rank
#    segm_local_size = il_extentx * il_extenty#   
#    return segm_local_size, segm_global_offset

#def def_paral(decomp, segm_global_offset, segm_local_size, il_extentx, il_extenty, nlon):
#    if decomp == 'DECOMP_APPLE':
#        return [1, segm_global_offset, segm_local_size]
#    elif decomp == 'DECOMP_BOX':
#        return [2, segm_global_offset, il_extentx, il_extenty, nlon]
#    else:
#        raise Exception(f'Invalid CPPKEYDECOMP in {__file__}: {decomp}')


def oasis_specify_partition(oasis_component,nlat,nlon):
        # --- Get local communicator ---
        local_communicator = oasis_component.localcomm
        # Get rank in local communicator
        lcomm_rank = local_communicator.rank
        lcomm_size = local_communicator.size
        # Unit for output messages : one file for each process
        w_unit = 100 + lcomm_rank
        componenet_outfile= oasis_component.name + '.out_'+str(w_unit)
        w_unit = open(componenet_outfile, 'w')
        print(' -----------------------------------------------------------', file=w_unit)
        print(f' I am atmos process with rank : {lcomm_rank}', file=w_unit)
        print(f' in my local communicator gathering {lcomm_size} processes', file=w_unit)
        print(' ----------------------------------------------------------', file=w_unit)
        w_unit.flush()

        # --- partition definition ---
        #segm_local_size, segm_global_offset = def_local_partition(rlon, rlat, lcomm_size, lcomm_rank
        if (lcomm_rank == lcomm_size-1):
            il_extenty = nlat - nlat//lcomm_size * lcomm_rank
        else:
            il_extenty = nlat//lcomm_size
        segm_local_size = nlon * il_extenty

        segm_global_offset = nlat//lcomm_size * lcomm_rank * nlon

        print(' Local partition definition', file=w_unit)
        #print(' il_extentx, il_extenty, segm_local_size, segm_global_offsetx, segm_global_offsety, segm_global_offset = ', il_extentx, il_extenty, segm_local_size, segm_global_offsetx, segm_global_offsety, segm_global_offset, file=w_unit)
        w_unit.flush()

        #CPPKEYDECOMP = 'DECOMP_APPLE'
        #ig_paral = def_paral(CPPKEYDECOMP, segm_global_offset, segm_local_size, il_extentx, il_extenty, rlon)
        ig_paral = [1, segm_global_offset, segm_local_size]
        print(' ig_paral = ', ' '.join(map(str, ig_paral)), file=w_unit)
        partition = pyoasis.ApplePartition(segm_global_offset, segm_local_size)

        return partition, w_unit

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
        # initialize component
        oasis_component = pyoasis.Component("hydro_component")
        # initialize pyoasis component, see toy model - ok
        # 1) get localcomm - ok
        # 2) partition definition - apple
        # 3) grid definition
        # 4) declaration of coupling fields

        ### get localcomm and define partition
        self.partition, self.w_unit = oasis_specify_partition(oasis_component,nlat=437,nlon=625)

        ### TERMINATION OF DEFINITION PHASE ###
        print(' End of initialisation phase', file=self.w_unit)
        self.w_unit.flush()

        ### OASIS_ENDDEF ###
        oasis_component.enddef()

    def dynamic(self):
        # 1) get
        # 2) put
        print(' Time loop', file=self.w_unit)
        self.w_unit.flush()
        pass




