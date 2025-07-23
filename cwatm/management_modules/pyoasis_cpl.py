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

class pyoasis_cpl(object):

    """
    Interface for coupling C-CWatM with other components via OASIS3-MCT_5.0.

    This class contains all functionality related to initializing and managing data exchange
    through the pyOASIS coupler. It handles grid setup, partitioning, and the declaration and
    transfer of coupling fields between C-CWatM and external models.

    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    conv_evap                              conversion factor for evaporation                                       --   
    conv_groundw                           conversion factor for groundwater recharge                              --
    conv_runoff                            conversion factor for runoff                                            --
    conv_soilw                             conversion factor for soilwater content                                 --
    EWRef                                  (potential) evaporation rate from water surface                         m   
    oasisvar_id
    runoff
    sum_gwRecharge                         groundwater recharge                                                    m 
    rootzoneSM                                         
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

    def initial(self):
        """
        Initialize the OASIS coupling interface for the C-CWatM hydrological component.

        This method performs the following steps:
        1. Initializes the OASIS component.
        2. Retrieves the local communicator and defines the grid partitioning.
        3. Constructs the grid coordinates, defines and writes OASIS grid.
        4. Declares the coupling fields expected to be received from the atmosphere model.

        The grid name and the declared coupling fields must match the one specified in the 'namcouple' file.

        Returns
        -------
        None : Defines the following variables
            - self.partition : A serial partition object defining the local grid partitioning for the component.
            - self.w_unit : A writable file object for logging output specific to the current process.
            - self.var.oasisvar_id : A list containing the variables that are exchanged via OASIS.
            
        as/copilot
        """

        # --- 1) initialize component ---
        oasis_component = pyoasis.Component("hydro_component")

        # --- 2) get local communicator and define partition ---
        nlon_cwatm = maskmapAttr['col']
        nlat_cwatm = maskmapAttr['row']
        print(nlat_cwatm,nlon_cwatm)
        self.partition, self.w_unit = oasis_specify_partition(oasis_component,nlon=nlon_cwatm,nlat=nlat_cwatm)

        # --- 3) grid definition ---
        # get 2d coordinates of c-cwatm grid
        lon_2d,lat_2d = create_2d_ccwatm_grid()
        grid_clon, grid_clat = derive_regular_grid_corners(lon_2d,lat_2d)
        print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=self.w_unit)
        print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=self.w_unit)
        self.w_unit.flush()
        # write oasis grid information
        oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,maskinfo['mask'].T.data[:,::-1],self.partition,'ccwatm_grid')

        # --- 4) declaration of coupling fields ---
        # needs to match namcouple
        numcouple = 4 # number of coupling fields
        self.var.oasisvar_id = [None] * numcouple
        self.var.oasisvar_id[0] = pyoasis.Var("ccwatm_recv_runoff", self.partition, OASIS.IN)
        self.var.oasisvar_id[1] = pyoasis.Var("ccwatm_recv_gwRecharge", self.partition, OASIS.IN)
        self.var.oasisvar_id[2] = pyoasis.Var("ccwatm_recv_EWRef", self.partition, OASIS.IN)
        self.var.oasisvar_id[3] = pyoasis.Var("ccwatm_recv_rootzoneSM", self.partition, OASIS.IN)
        print(f' self.var.oasisvar_id FIELD_RECV_hydro, {self.var.oasisvar_id[0]._id}', file=self.w_unit)
        self.w_unit.flush()
  
        # --- 5) termination of definition phase ---
        print(' End of initialisation phase', file=self.w_unit)
        self.w_unit.flush()
        oasis_component.enddef()

    
    def dynamic(self):
        """
        Receive and send data via OASIS.

        Receive forcing data for runoff, groundwater recharge, evaporation over water, and soil water content.
        (To be implemented ?? :) Send irriation water amount. 
            
        Uses global variables
        ---------------------
        dateVar['currDate'] (datetime object): Current date.
        dateVar['dateStart'] (datetime object) : Starting date of model run.
        maskmapAttr['col'], maskmapAttr['row'] (int) : Grid dimensions.
        maskinfo['mask'] : The landmask used to filter valid grid cells.

        Returns
        -------
        None : The method updates the following global variables in-place:
            - runoff
            - sum_gwRecharge
            - EWRef
            - rootzoneSM    
        """
        # TODO: might need to be divided into 2 parts for model forcing and irriwater - to be called at different places
        # in cwatm_dynamic. -> check later
        
        # Calculate the passed time in seconds
        seconds_passed = int((dateVar['currDate'] - dateVar['dateStart']).total_seconds())
        print('C-CWatM time loop',seconds_passed, file=self.w_unit)
        self.w_unit.flush()
        
        # ----- 1) OASIS get -----
        # same order as declared in pyoasis_cpl.initial(); needs to match namcouple
        field_recv_runoff = pyoasis.asarray(np.full((maskmapAttr['col'], maskmapAttr['row']), -1.0))
        self.var.oasisvar_id[0].get(seconds_passed, field_recv_runoff)
        field_recv_gwRecharge = pyoasis.asarray(np.full((maskmapAttr['col'], maskmapAttr['row']), -1.0))
        self.var.oasisvar_id[1].get(seconds_passed, field_recv_gwRecharge)
        field_recv_EWRef = pyoasis.asarray(np.full((maskmapAttr['col'], maskmapAttr['row']), -1.0))
        self.var.oasisvar_id[2].get(seconds_passed, field_recv_EWRef)
        field_recv_rootzoneSM = pyoasis.asarray(np.full((maskmapAttr['col'], maskmapAttr['row']), -1.0))
        self.var.oasisvar_id[3].get(seconds_passed, field_recv_rootzoneSM)

        # --- assign to C-CWatM variables ---
        # apply landmask and conversion factors
        # runoff
        mapnp1 = np.ma.masked_array(field_recv_runoff.T[::-1,:], maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.runoff = np.maximum(0., mapnp1 * self.var.conv_runoff)
        # groundwater recharge
        mapnp1 = np.ma.masked_array(field_recv_gwRecharge.T[::-1,:], maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.sum_gwRecharge = np.maximum(0., mapnp1 * self.var.conv_groundw)
        # evaporation over water
        mapnp1 = np.ma.masked_array(field_recv_EWRef.T[::-1,:], maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.EWRef = mapnp1 * self.var.conv_evap 
        # soil water content
        mapnp1 = np.ma.masked_array(field_recv_rootzoneSM.T[::-1,:], maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.rootzoneSM = np.maximum(0., mapnp1 * self.var.conv_soilw)
        # TODO: the variables should not be flipped here, but rather already before sending -> easier adjustment
        # for other models/data 1) adjust grid here 2) adjust sending grid 3) adjust sending variables

        # ----- 2) OASIS put -----

        # -> TODO: irrigation water



# ------- OASIS3-MCT functions -------

def oasis_specify_partition(oasis_component,nlon,nlat):
    """
    Defines the local partitioning for an OASIS component and sets up process-specific logging.

    This function assumes a serial partitioning scheme (no parallel decomposition is applied).
    It creates a separate log file for each process, named `<component_name>.out_<rank>`.

    Parameters
    ----------
    oasis_component : An instance of an OASIS component, expected to have attributes `localcomm` (MPI communicator)
            and `name` (component name).
    nlon (int) : Number of longitudinal grid points.
    nlat (int) : Number of latitudinal grid points.

    Returns
    -------
    partition : A serial partition object defining the local grid partitioning for the component.
    w_unit (file object) : A writable file object for logging output specific to the current process.

    Called in
    ---------
     - pyoasis_cpl.initial -> initalize C-CWatM OASIS componenet
     - run_oasis_dummy.py -> initialize OASIS component for reading of forcing 

    as/copilot
    """
    # get local communicator
    local_communicator = oasis_component.localcomm
    # get rank in local communicator
    lcomm_rank = local_communicator.rank
    lcomm_size = local_communicator.size

    # --- open logging file ---
    # unit for output messages : one file for each process
    w_unit_number = 100 + lcomm_rank
    component_outfile = oasis_component.name + '.out_'+str(w_unit_number)
    w_unit = open(component_outfile, 'w')
    print(' -----------------------------------------------------------', file=w_unit)
    print(f' I am {oasis_component.name} process with rank : {lcomm_rank}', file=w_unit)
    print(f' in my local communicator gathering {lcomm_size} processes', file=w_unit)
    print(' ----------------------------------------------------------', file=w_unit)
    print(' Local partition definition', file=w_unit)
    w_unit.flush()

    # --- create partition ---
    # use serial partitioning for C-CWatM (and box for REMO)
    partition = pyoasis.SerialPartition(nlon*nlat)
        
    return partition, w_unit


def oasis_define_grid(nlon,nlat,grid_lon,grid_lat,landmask,partition,grid_name,
                      grid_clon=None,grid_clat=None,grid_area=None):	
    """
    Define and configure a grid for OASIS coupling.

    This function initializes a 'pyoasis.Grid' object with the provided grid dimensions,
    coordinates, land-sea mask, and optional corner and area data. It writes the grid
    definition to file for use in OASIS coupling.

    Parameters
    ----------
    nlon (int) : Number of grid points in the longitudinal direction.
    nlat (int) : Number of grid points in the latitudinal direction.
    grid_lon (ndarray) : 2D array of longitudes at the center of each grid cell.
    grid_lat (ndarray) : 2D array of latitudes at the center of each grid cell.
    landmask (ndarray) : 2D array indicating land (0) and sea (1) points.
    partition : A pyoasis partition object defining the local grid partitioning for the component.
    grid_name (str) : Name of the grid. Must match the name used in the 'namcouple' file.
    Optional:
    grid_clon (ndarray) : 3D array of longitudes at the corners of each grid cell.
    grid_clat (ndarray) : 3D array of latitudes at the corners of each grid cell.
    grid_area (ndarray) : 2D array of grid cell areas.

    as/copilot
    """
    # define grid name and center coordinates 
    grid = pyoasis.Grid(grid_name, nlon, nlat, grid_lon, grid_lat, partition)
    # define land-sea mask
    grid.set_mask(landmask)

    # different grid specifications are requied depending on the remapping method 
    # set grid cell corners if given
    if grid_clon is not None and grid_clat is not None:
        grid.set_corners(grid_clon, grid_clat)
    # set grid cell areas if given
    if grid_area is not None:
        grid.set_area(grid_area)
	
    grid.write()


# ------- gridding functions -------

def create_2d_ccwatm_grid():
    """
    Generate 2D arrays of longitude and latitude coordinates for a C-CWatM grid.

    This function computes the 2D grid coordinates based on the attributes defined in
    'maskmapAttr'. The latitude array is constructed starting from the northernmost point.

    Uses global variables
    ---------------------
    maskmapAttr['x'] : upper left corner of longitudes ??
    maskmapAttr['y'] : upper left corner of latitudes ??
    maskmapAttr['cell'] : Grid cell length (degrees).
    maskmapAttr['col'], maskmapAttr['row'] (int) : Grid dimensions.

    Returns
    -------
    lon_2d (ndarray) : 2D array of center longitudes for each grid cell.
    lat_2d (ndarray) : 2D array of latitudes for each grid cell.

    as/copilot
    """

    # TODO: check whether upper left corner or center coordinates
    lon_1d = maskmapAttr['x'] + maskmapAttr['cell']*np.arange(maskmapAttr['col'])
    # lat starts with the northernmost latitude
    lat_1d = maskmapAttr['y'] - maskmapAttr['cell']*np.arange(maskmapAttr['row'])

    lat_2d,lon_2d = np.meshgrid(lat_1d[::-1],lon_1d)

    return lon_2d,lat_2d


def derive_regular_grid_corners(lon,lat):
    """
    Compute the corner coordinates of each grid cell from center coordinates on a regular grid.
    The corners are ordered counter-clockwise, as required by OASIS.

    Parameters
    ----------
    lon (ndarray) : 2D array of longitudes at the center of each grid cell.
    lat (ndarray) : 2D array of latitudes at the center of each grid cell.

    Returns
    -------
    clon (ndarray) : 3D array of longitudes at the corners of each grid cell.
    clan (ndarray) : 3D array of latitudes at the corners of each grid cell.

    as/copilot
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

