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

        # --- partition definition ---
        if (lcomm_rank == lcomm_size-1):
            il_extenty = nlat - nlat//lcomm_size * lcomm_rank
        else:
            il_extenty = nlat//lcomm_size
        segm_local_size = nlon * il_extenty
        segm_global_offset = nlat//lcomm_size * lcomm_rank * nlon

        print(' Local partition definition', file=w_unit)
        w_unit.flush()

        # what is ig_paral used for??
        ig_paral = [1, segm_global_offset, segm_local_size]
        print(' ig_paral = ', ' '.join(map(str, ig_paral)), file=w_unit)
        w_unit.flush()

        partition = pyoasis.ApplePartition(segm_global_offset, segm_local_size)
        
        return partition, w_unit

def oasis_define_grid():	
	#  GRID DEFINITION
	# Reading local grid arrays from input file atmos_mesh.nc
	grid_lon_atmos, grid_lat_atmos, grid_clo_atmos, grid_cla_atmos, grid_srf_atmos, grid_msk_atmos = read_grid(il_offsetx, il_offsety, il_extentx, il_extenty, nc_atmos, 'atmos_mesh.nc')

	# OASIS_WRITE_GRID  
	grid = pyoasis.Grid('lmdz', nlon_atmos, nlat_atmos, grid_lon_atmos, grid_lat_atmos, partition)
	grid.set_corners(grid_clo_atmos, grid_cla_atmos)
	grid.set_mask(grid_msk_atmos)
	grid.write()

	print(f' grid_lat_atmos maximum and minimum', '%.5f' % np.max(grid_lat_atmos), '%.5f' % np.min(grid_lat_atmos), file=w_unit)
	w_unit.flush()


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

# self.MaskMap = loadsetclone(self, 'MaskMap')
# ?? MaskMap is never used again...

#latmeteo, lonmeteo, cell, invcellmeteo, rows, cols = readCoordNetCDF(namemeteo)

#nameldd = cbinding('Ldd')
#latldd, lonldd, cell, invcellldd, rows, cols = readCoord(nameldd)
#maskmapAttr['reso_mask_meteo'] = round(invcellldd / invcellmeteo)

#y, x, cell, invcell, rows, cols = readCoordNetCDF(filename)
#setmaskmapAttr(x, y, cols, rows, cell)

#filename = cbinding(name)
#nf2 = gdal.Open(filename, gdalconst.GA_ReadOnly)
#geotransform = nf2.GetGeoTransform()
#geotrans.append(geotransform)
#setmaskmapAttr( geotransform[0], geotransform[3], nf2.RasterXSize, nf2.RasterYSize, geotransform[1])

# Default behaviour: grid size is derived from location attributes of
# base maps. Requirements:
# - Maps are in some equal-area projection
# - Length units meters
# - All grid cells have the same size

# Area of pixel [m2]
#self.var.cellArea=np.empty(maskinfo['mapC'])
#self.var.cellArea.fill(maskmapAttr['cell'] ** 2)

#maskmapAttr['coordx'] = 'lon'
#maskmapAttr['coordy'] = 'lat'

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

        # 3) grid definition TODO

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




