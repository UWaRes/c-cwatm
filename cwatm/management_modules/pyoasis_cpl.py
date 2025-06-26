# -------------------------------------------------------------------------
# Name:        pyoasis module
# Purpose:
#
# Author:      Amelie Schmitt
#
# Created:     25/02/2025
# Copyright:   (c) Amelie Schmitt 2025 
# -------------------------------------------------------------------------

import pyoasis
import numpy as np
from cwatm.management_modules.globals import maskmapAttr, maskinfo

def oasis_specify_partition(oasis_component,nlon,nlat):
    # get local communicator
    local_communicator = oasis_component.localcomm
    # get rank in local communicator
    lcomm_rank = local_communicator.rank
    lcomm_size = local_communicator.size
    # unit for output messages : one file for each process
    w_unit = 100 + lcomm_rank
    component_outfile= oasis_component.name + '.out_'+str(w_unit)
    w_unit = open(component_outfile, 'w')
    print(' -----------------------------------------------------------', file=w_unit)
    print(f' I am {oasis_component.name} process with rank : {lcomm_rank}', file=w_unit)
    print(f' in my local communicator gathering {lcomm_size} processes', file=w_unit)
    print(' ----------------------------------------------------------', file=w_unit)
    w_unit.flush()

    print(' Local partition definition', file=w_unit)
    w_unit.flush()

    # use serial partitioning for C-CWatM (and box or apple for REMO)
    partition = pyoasis.SerialPartition(nlon*nlat)
        
    return partition, w_unit


def oasis_define_grid_simple(nlon,nlat,grid_lon,grid_lat,landmask,partition,grid_name):	
    # define grid name and center coordinates 
    grid = pyoasis.Grid(grid_name, nlon, nlat, grid_lon, grid_lat, partition)
    # define land-sea mask TODO: find mask variable in cwatm
    grid.set_mask(landmask)

    # -- different grid specifications are requied depending on the rempaping method --
    # for first test, use nearest-neighbor or bilinear
    # coners are only required for some conservative or ESMF-based remapping methods 
    # that use cell shapes rather than just centers
    #grid.set_corners(grid_clon, grid_clat)
	
    grid.write()

	


# self.MaskMap = loadsetclone(self, 'MaskMap')
# ?? MaskMap is never used again...


def create_2d_ccwatm_grid():
    """
    Create 2D arrays with grid coordinates from given C-CWatM grid specifications
    """
    # TODO: check ascending or descending
    lon_1d = maskmapAttr['x'] + maskmapAttr['cell']*np.arange(maskmapAttr['col'])
    lat_1d = maskmapAttr['y'] + maskmapAttr['cell']*np.arange(maskmapAttr['row'])

    lon_2d,lat_2d = np.meshgrid(lon_1d,lat_1d)

    return lon_2d,lat_2d


def derive_regular_grid_corners(lon,lat):
    """
    Derive the coordinates of the 4 corners of each grid box for given
    center longitude and latitude.
    """
    dx = np.abs(lon[1,1] - lon[0,0])
    dy = np.abs(lat[1,1] - lat[0,0])   

    clon = pyoasis.asarray(np.zeros((lon.shape[0], lon.shape[1], 4), dtype=np.float64))
    clon[:, :, 0] = lon[:, :] - dx/2.0
    clon[:, :, 1] = lon[:, :] + dx/2.0
    clon[:, :, 2] = clon[:, :, 1]
    clon[:, :, 3] = clon[:, :, 0]
    clan = pyoasis.asarray(np.zeros((lon.shape[0], lon.shape[1], 4), dtype=np.float64))
    clan[:, :, 0] = lat[:, :] - dy/2.0
    clan[:, :, 1] = clan[:, :, 0]
    clan[:, :, 2] = lat[:, :] + dy/2.0
    clan[:, :, 3] = clan[:, :, 2]

    return clon, clan

# ---------------------------------------------

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
        # 1) initialize component
        oasis_component = pyoasis.Component("hydro_component")

        # 2) get local communicator and define partition
        nlon_cwatm = maskmapAttr['col']
        nlat_cwatm = maskmapAttr['row']
        self.partition, self.w_unit = oasis_specify_partition(oasis_component,nlon=nlon_cwatm,nlat=nlat_cwatm)

        # 3) grid definition
        # get 2d coordinates of c-cwatm grid
        lon_2d,lat_2d = create_2d_ccwatm_grid()
        print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=self.w_unit)
        print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=self.w_unit)
        self.w_unit.flush()
        # get grid corner locations
        corner_lon,corner_lat = derive_regular_grid_corners(lon_2d,lat_2d)
        # function for writing oasis grid information
        oasis_define_grid_simple(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,maskinfo['mask'],self.partition,'ccwatm_grid')

        # load ldd for land-sea mask
        # TODO: needs to be cut to correct domain
        #print('sum mask:',np.sum(~maskinfo['mask']))
        #print('mask:',maskinfo['mask'].shape)
        #print('maskflat:',maskinfo['maskflat'].shape)
        #print('sum maskflat:', np.sum(~maskinfo['maskflat']))
        #print('gridsize:', maskmapAttr['col']*maskmapAttr['row'])

        # -> maskinfo['mask'] is the land-sea mask


        # 4) declaration of coupling fields TODO

        # 5) termination of definition phase
        print(' End of initialisation phase', file=self.w_unit)
        self.w_unit.flush()
        oasis_component.enddef()

    def dynamic(self):
        # 1) get
        # 2) put
        print(' Time loop', file=self.w_unit)
        self.w_unit.flush()
        pass




