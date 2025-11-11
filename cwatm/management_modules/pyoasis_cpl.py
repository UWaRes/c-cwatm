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
import time

from cwatm.management_modules.globals import maskmapAttr, maskinfo, dateVar
from cwatm.management_modules.grid_tools import grid_tools

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
        self.partition, self.w_unit = self.oasis_specify_partition(oasis_component,nlon=nlon_cwatm,nlat=nlat_cwatm)

        # --- 3) grid definition ---
        # get 2d coordinates of c-cwatm grid
        ccgrid = grid_tools(maskmapAttr)
        lon_2d,lat_2d = ccgrid.create_2d_ccwatm_grid()
        # optional: get grid cell corners and area (only required for certain regridding methods)
        grid_clon, grid_clat = grid_tools.compute_grid_corners(lon_2d,lat_2d)
        #cell_areas = grid_tools.compute_grid_cell_areas(grid_clon, grid_clat)
        print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=self.w_unit)
        print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=self.w_unit)
        self.w_unit.flush()
        # write oasis grid information
        #self.oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,np.fliplr(maskinfo['mask'].T),self.partition,'ccwatm_grid')
        self.oasis_define_grid(nlon_cwatm,nlat_cwatm,lon_2d,lat_2d,np.fliplr(maskinfo['mask'].T),self.partition,'ccwatm_grid',grid_clon,grid_clat))

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
        (To be implemented ?? :) Send irrigation water amount. 
            
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
        # received fields need to be flipped and transposed to match the C-CWatM grid orientation
        # runoff
        mapnp1 = np.ma.masked_array(np.flipud(field_recv_runoff.T), maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.runoff = np.maximum(0., mapnp1 * self.var.conv_runoff)
        # groundwater recharge
        mapnp1 = np.ma.masked_array(np.flipud(field_recv_gwRecharge.T), maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.sum_gwRecharge = np.maximum(0., mapnp1 * self.var.conv_groundw)
        # evaporation over water
        mapnp1 = np.ma.masked_array(np.flipud(field_recv_EWRef.T), maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.EWRef = mapnp1 * self.var.conv_evap 
        # soil water content
        mapnp1 = np.ma.masked_array(np.flipud(field_recv_rootzoneSM.T), maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.rootzoneSM = np.maximum(0., mapnp1 * self.var.conv_soilw)

        # ----- 2) OASIS put -----

        # -> TODO: irrigation water
        # new flag: coupl_irri

        
        # --- write to cwatm variables ---
        mapnp1 = np.ma.masked_array(field_recv_runoff.T[::-1,:], maskinfo['mask'])
        mapnp1 = np.ma.compressed(mapnp1)
        self.var.runoff = np.maximum(0., mapnp1)

    # ------- OASIS3-MCT functions -------
    @staticmethod 
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
         - pyoasis_cpl.initial -> initalize C-CWatM OASIS component
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

    @staticmethod
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

        Called in
        ---------
         - pyoasis_cpl.initial -> write C-CWatM grid information for OASIS
         - run_oasis_dummy.py -> write forcing grid information for OASIS 

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





