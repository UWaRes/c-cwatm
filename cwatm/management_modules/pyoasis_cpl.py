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
from cwatm.management_modules.globals import maskmapAttr

def oasis_specify_partition(oasis_component,nlat,nlon):
        # --- Get local communicator ---
        local_communicator = oasis_component.localcomm
        # Get rank in local communicator
        lcomm_rank = local_communicator.rank
        lcomm_size = local_communicator.size
        # Unit for output messages : one file for each process
        w_unit = 100 + lcomm_rank
        component_outfile= oasis_component.name + '.out_'+str(w_unit)
        w_unit = open(component_outfile, 'w')
        print(' -----------------------------------------------------------', file=w_unit)
        print(f' I am {oasis_component.name} process with rank : {lcomm_rank}', file=w_unit)
        print(f' in my local communicator gathering {lcomm_size} processes', file=w_unit)
        print(' ----------------------------------------------------------', file=w_unit)
        w_unit.flush()

        if 0:
           # old code, only used for apple partitioning
            # --- partition definition ---
            if (lcomm_rank == lcomm_size-1):
                il_extenty = nlat - nlat//lcomm_size * lcomm_rank
            else:
                il_extenty = nlat//lcomm_size
            segm_local_size = nlon * il_extenty
            segm_global_offset = nlat//lcomm_size * lcomm_rank * nlon

            print(' Local partition definition', file=w_unit)
            w_unit.flush()

            # what is ig_paral used for? -> necessary only in fortran
            ig_paral = [1, segm_global_offset, segm_local_size]
            print(' ig_paral = ', ' '.join(map(str, ig_paral)), file=w_unit)
            w_unit.flush()

            partition = pyoasis.ApplePartition(segm_global_offset, segm_local_size)

        # TODO: use serial partitioning for C-CWatM (and box or apple for REMO)
        partition = pyoasis.SerialPartition(nlon*lat)

        
        return partition, w_unit

def oasis_define_grid():	
	#  GRID DEFINITION
	# Reading local grid arrays from input file atmos_mesh.nc
	grid_lon_atmos, grid_lat_atmos, grid_clo_atmos, grid_cla_atmos, grid_srf_atmos, grid_msk_atmos = read_grid(il_offsetx, il_offsety, il_extentx, il_extenty, nc_atmos, 'atmos_mesh.nc')

        # grid_lon_atmos - dda_lon  2d - lon coordinates
        # grid_lat_atmos - dda_lat  2d - lat coordinates
        # grid_clo_atmos - dda_clo  3d - 4 corners of each lon coordinate
        # grid_cla_atmos - dda_cla  3d - 4 conerns of each lat coordinate
        # grid_srf_atmos - dda_srf  2d - grid cell area (not necessary for regular grid)
        # grid_msk_atmos - ida_mask  2d - mask


	# OASIS_WRITE_GRID  
	grid = pyoasis.Grid('lmdz', nlon_atmos, nlat_atmos, grid_lon_atmos, grid_lat_atmos, partition)
	grid.set_corners(grid_clo_atmos, grid_cla_atmos)
	grid.set_mask(grid_msk_atmos)
	grid.write()

	print(f' grid_lat_atmos maximum and minimum', '%.5f' % np.max(grid_lat_atmos), '%.5f' % np.min(grid_lat_atmos), file=w_unit)
	w_unit.flush()

def read_grid(id_begi, id_begj, id_lon, id_lat, nc, data_filename):
    gf = netCDF4.Dataset(data_filename, 'r')
    # Keep in mind Fortran is column-major while Python is row-major
    # We are switching to column-major
    dda_lon = np.transpose(gf['lon'][id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    dda_lat = np.transpose(gf['lat'][id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    dda_clo = np.transpose(gf['clo'][:nc,id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    dda_cla = np.transpose(gf['cla'][:nc,id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    dda_srf = np.transpose(gf['srf'][id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    ida_mask = np.transpose(gf['imask'][id_begj:id_begj+id_lat,id_begi:id_begi+id_lon])
    gf.close()

    # OASIS3 mask convention (1=masked, 0=not masked) is opposite to usual one)
    ida_mask ^= 1

    return dda_lon, dda_lat, dda_clo, dda_cla, dda_srf, ida_mask

def oasis_define_local_partition(nlon, nlat, lcomm_size, lcomm_rank):
    # divided only in lat direction ?? does this also work for REMO??
    extentx = nlon
    extenty = nlat//lcomm_size
    if (lcomm_rank == lcomm_size-1):
        extenty = nlat - nlat//lcomm_size * lcomm_rank
    local_size = extentx * extenty
    offsetx = 0
    offsety = nlat//lcomm_size * lcomm_rank
    local_offset = nlon * offsety
    return extentx, extenty, local_size, offsetx, offsety, local_offset



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

    clo = pyoasis.asarray(np.zeros((lon.shape[0], lon.shape[1], 4), dtype=np.float64))
    clo[:, :, 0] = lon[:, :] - dx/2.0
    clo[:, :, 1] = lon[:, :] + dx/2.0
    clo[:, :, 2] = clo[:, :, 1]
    clo[:, :, 3] = clo[:, :, 0]
    cla = pyoasis.asarray(np.zeros((lon.shape[0], lon.shape[1], 4), dtype=np.float64))
    cla[:, :, 0] = lat[:, :] - dy/2.0
    cla[:, :, 1] = cla[:, :, 0]
    cla[:, :, 2] = lat[:, :] + dy/2.0
    cla[:, :, 3] = cla[:, :, 2]

    return clo, cla



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
        clon = maskmapAttr['col']
        clat = maskmapAttr['row']
        self.partition, self.w_unit = oasis_specify_partition(oasis_component,nlat=clat,nlon=clon)

        # 3) grid definition

        # get 2d coordinates of c-cwatm grid
        lon_2d,lat_2d = create_2d_ccwatm_grid()
        # get grid corner locations
        clo,cla = derive_regular_grid_corners(lon_2d,lat_2d)
        # TODO: function for 

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




