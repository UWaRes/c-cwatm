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
from pyoasis import OASIS
import numpy as np
from cwatm.management_modules.globals import maskmapAttr, maskinfo, dateVar

# for testing
from netCDF4 import Dataset

# ------- OASIS3-MCT functions -------

def oasis_specify_partition(oasis_component,nlon,nlat):
    """
    TODO
    """
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
    """
    TODO
    """
    # define grid name and center coordinates 
    grid = pyoasis.Grid(grid_name, nlon, nlat, grid_lon, grid_lat, partition)
    # define land-sea mask
    grid.set_mask(landmask)

    # -- different grid specifications are requied depending on the remapping method --
    # for first test, use nearest-neighbor or bilinear
    # corners are only required for some conservative or ESMF-based remapping methods 
    # that use cell shapes rather than just centers
    #grid.set_corners(grid_clon, grid_clat)

def oasis_define_grid(nlon,nlat,grid_lon,grid_lat,grid_clon,grid_clat,landmask,partition,grid_name):	
    """
    TODO
    """
    # define grid name and center coordinates 
    grid = pyoasis.Grid(grid_name, nlon, nlat, grid_lon, grid_lat, partition)
    # define land-sea mask
    grid.set_mask(landmask)

    # -- different grid specifications are requied depending on the remapping method --
    # for first test, use nearest-neighbor or bilinear
    # corners are only required for some conservative or ESMF-based remapping methods 
    # that use cell shapes rather than just centers
    grid.set_corners(grid_clon, grid_clat)
	
    grid.write()


# ------- gridding functions -------

def create_2d_ccwatm_grid():
    """
    Create 2D arrays with grid coordinates from given C-CWatM grid specifications
    """
    # TODO: check ascending or descending
    lon_1d = maskmapAttr['x'] + maskmapAttr['cell']*np.arange(maskmapAttr['col'])
    # lat starts with the northernmost latitude
    #lat_1d = maskmapAttr['y'] - maskmapAttr['cell']*np.arange(maskmapAttr['row'])
    lat_1d = maskmapAttr['y'] - maskmapAttr['cell']*np.arange(maskmapAttr['row'])

    #lon_2d,lat_2d = np.meshgrid(lon_1d,lat_1d[::-1]) # test1
    lat_2d,lon_2d = np.meshgrid(lat_1d[::-1],lon_1d)
    #lon_2d,lat_2d = np.meshgrid(lon_1d,lat_1d)
    #lon_2d,lat_2d = np.meshgrid(lon_1d[::-1],lat_1d) # test2
    #lon_2d,lat_2d = np.meshgrid(lon_1d,lat_1d[::-1]) # test3
    #lon_2d,lat_2d = np.meshgrid(lon_1d[::-1],lat_1d[::-1]) # test4

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
        print(nlat_cwatm,nlon_cwatm)
        self.partition, self.w_unit = oasis_specify_partition(oasis_component,nlon=nlon_cwatm,nlat=nlat_cwatm)

        # 3) grid definition
        # get 2d coordinates of c-cwatm grid
        lon_2d,lat_2d = create_2d_ccwatm_grid()
        mapnp1 = np.ma.masked_array(lon_2d, maskinfo['mask'])
        self.var.londummy = np.ma.compressed(mapnp1)
        mapnp1 = np.ma.masked_array(lat_2d, maskinfo['mask'])
        self.var.latdummy = np.ma.compressed(mapnp1)
        #self.var.lontest = lon_2d.copy()
        #self.var.lattest =lat_2d.copy()
        #print('lon shape',self.var.lontest.shape)
        #mapnp1 = np.ma.masked_array(lon_2d, maskinfo['mask'])
        #mapC = np.ma.compressed(mapnp1)
        #print('lon compressed shape',mapC.shape)

        grid_clon, grid_clat = derive_regular_grid_corners(lon_2d,lat_2d)
        
        print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=self.w_unit)
        print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=self.w_unit)
        self.w_unit.flush()
        # write oasis grid information
        # TODO: check 0 and 1 in mask
        dummymask = np.ones(maskinfo['mask'].T.data.shape)
        #oasis_define_grid_simple(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,dummymask,self.partition,'ccwatm_grid')
        #oasis_define_grid_simple(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,1-maskinfo['mask'].T.data,self.partition,'ccwatm_grid')
        # dummy
        #oasis_define_grid_simple(nlat_cwatm,nlon_cwatm,lon_2d,lat_2d,maskinfo['mask'],self.partition,'ccwatm_grid')

        # with corners
        #oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,grid_clon,grid_clat,1-maskinfo['mask'].T.data,self.partition,'ccwatm_grid') 
        #oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,grid_clon,grid_clat,1-maskinfo['mask'].T.data[:,::-1],self.partition,'ccwatm_grid')
        oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,grid_clon,grid_clat,maskinfo['mask'].T.data[:,::-1],self.partition,'ccwatm_grid')

        # 4) declaration of coupling fields
        self.oasisvar_id = [None]
        self.oasisvar_id[0] = pyoasis.Var("FIELD_RECV_hydro", self.partition, OASIS.IN)
        print(f' self.oasisvar_id FIELD_RECV_hydro, {self.oasisvar_id[0]._id}', file=self.w_unit)
        self.w_unit.flush()

        # 5) termination of definition phase
        print(' End of initialisation phase', file=self.w_unit)
        self.w_unit.flush()
        oasis_component.enddef()

    def dynamic(self):
        # Calculate the time passed in seconds
        seconds_passed = int((dateVar['currDate'] - dateVar['dateStart']).total_seconds())
        print('C-CWatM time loop',seconds_passed, file=self.w_unit)
        self.w_unit.flush()
        # 1) get
        field_recv_atmos = pyoasis.asarray(np.full((maskmapAttr['col'], maskmapAttr['row']), -1.0))
        #field_recv_atmos = pyoasis.asarray(np.full((maskmapAttr['row'], maskmapAttr['col']), -1.0))
        self.oasisvar_id[0].get(seconds_passed, field_recv_atmos)

        #mapnp1 = np.ma.masked_array(field_recv_atmos[:,::-1], maskinfo['mask'])
        mapnp1 = np.ma.masked_array(field_recv_atmos.T[:,::-1], maskinfo['mask']) # this gives non-zero output
        #mapnp1 = np.ma.masked_array(field_recv_atmos.T[:,::-1], maskinfo['mask'])
        #mapnp1 = np.ma.masked_array(field_recv_atmos, maskinfo['mask'])
        print('recv shape:',field_recv_atmos.shape)
        print('mask shape:',maskinfo['mask'].shape)
        #mapnp1 = np.ma.masked_array(field_recv_atmos, maskinfo['mask'].T[:,::-1])
        self.var.oasisdummy = np.ma.compressed(mapnp1)
        #mapnp1 = np.ma.masked_array(field_recv_atmos.T, maskinfo['mask'])
        #self.var.oasisdummy = np.ma.compressed(field_recv_atmos)
        #print('lon compressed shape',mapC.shape)

       
        nf1 = Dataset('test_dummy.nc', 'w', format='NETCDF4')
        row,col = field_recv_atmos.shape
        lat = nf1.createDimension('lat',row)
        lon = nf1.createDimension('lon',col)
        #nf1.createDimension('time', 5)
        value = nf1.createVariable('dummyvar', 'f4', ('lat', 'lon'), zlib=True, fill_value=1e20,chunksizes=(row,col))
        nf1.variables['dummyvar'][:, :] = field_recv_atmos
        nf1.close()

        # 2) put





