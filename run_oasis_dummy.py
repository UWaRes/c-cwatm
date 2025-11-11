import pyoasis
from pyoasis import OASIS
import xarray as xr
import numpy as np

import sys
import os

# Add the relative path to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './cwatm/management_modules/')))

from pyoasis_cpl import oasis_specify_partition, oasis_define_grid_simple, oasis_define_grid, derive_regular_grid_corners


# --- function for grid cell corners ---

def unrot_lon(rlat, rlon, pole_lat, pole_lon):
    """
    Transform rotated longitude to regular non-rotated longitude lon(2D).
    """
    nrlat = np.shape(rlat)
    nrlon = np.shape(rlon)
    rla = rlat
    rlo = rlon

    rla = np.deg2rad(rla)
    rlo = np.deg2rad(rlo)
    s1 = np.sin(np.deg2rad(pole_lat))
    c1 = np.cos(np.deg2rad(pole_lat))
    s2 = np.sin(np.deg2rad(pole_lon))
    c2 = np.cos(np.deg2rad(pole_lon))

    tmp1 = s2*(-s1*np.cos(rlo)*np.cos(rla)+c1*np.sin(rla))-c2*np.sin(rlo)*np.cos(rla)
    tmp2 = c2*(-s1*np.cos(rlo)*np.cos(rla)+c1*np.sin(rla))+s2*np.sin(rlo)*np.cos(rla)

    lon = np.rad2deg(np.arctan(tmp1/tmp2))
    return lon

def unrot_lat(rlat, rlon, pole_lat, pole_lon):
    """
    Transform rotated latitude to regular non-rotated latitude lat(2D)
    """
    nrlat = np.shape(rlat)
    nrlon = np.shape(rlon)

    rla = rlat
    rlo = rlon

    rla = np.deg2rad(rla)
    rlo = np.deg2rad(rlo)
    s1 = np.sin(np.deg2rad(pole_lat))
    c1 = np.cos(np.deg2rad(pole_lat))
    lat = s1*np.sin(rla)+c1*np.cos(rla)*np.cos(rlo)
    lat = np.rad2deg(np.arcsin(lat))

    return lat

def calc_rotgrid_corners(rlon,rlat):
    """
    Calculate the corner coordinates in rotated coordinates
    """
    deltalon = rlon[0,1]-rlon[0,0]
    deltalat = rlat[0,1]-rlat[0,0]

    rlon_corners = np.zeros((rlon.shape[0],rlon.shape[1],4))
    rlat_corners = np.zeros((rlat.shape[0],rlat.shape[1],4))

    rlon_corners[:,:,0] = rlon - deltalon/2.
    rlon_corners[:,:,1] = rlon + deltalon/2.
    rlon_corners[:,:,2] = rlon + deltalon/2.
    rlon_corners[:,:,3] = rlon - deltalon/2.
                    
    rlat_corners[:,:,0] = rlat - deltalat/2. 
    rlat_corners[:,:,1] = rlat - deltalat/2. 
    rlat_corners[:,:,2] = rlat + deltalat/2. 
    rlat_corners[:,:,3] = rlat + deltalat/2. 

    return rlon_corners, rlat_corners


#unrot_lon(rlat_corners[0,:], rlon_corners[0,:], ds.rotated_latitude_longitude.grid_north_pole_latitude, ds.rotated_latitude_longitude.grid_north_pole_longitude)

# --- functions for grid cell area ---

def earth_radius_at_latitude(latitude_deg):
    # Convert latitude to radians
    lat = np.radians(latitude_deg)
    
    # WGS-84 ellipsoid constants
    a = 6378.137  # Equatorial radius in km
    b = 6356.752  # Polar radius in km

    # Compute radius at latitude using ellipsoidal model
    numerator = (a**2 * np.cos(lat))**2 + (b**2 * np.sin(lat))**2
    denominator = (a * np.cos(lat))**2 + (b * np.sin(lat))**2
    radius = np.sqrt(numerator / denominator)
    
    return radius  # in kilometers

def to_unit_vector(lat, lon):
    """Convert lat/lon in degrees to 3D unit vector."""
    lat = np.radians(lat)
    lon = np.radians(lon)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.array([x, y, z])

def angle_between(a, b, c):
    """Compute the spherical angle at point b using great-circle arcs."""
    ab = np.cross(b, a)
    cb = np.cross(b, c)
    ab /= np.linalg.norm(ab)
    cb /= np.linalg.norm(cb)
    angle = np.arccos(np.clip(np.dot(ab, cb), -1.0, 1.0))
    return angle

def spherical_triangle_area(a, b, c):
    """Compute the area of a spherical triangle on a unit sphere using spherical excess."""
    # Compute angles at each vertex
    A = angle_between(b, a, c)
    B = angle_between(c, b, a)
    C = angle_between(a, c, b)
    # Spherical excess
    E = A + B + C - np.pi
    return E

def spherical_polygon_area(coords, radius=1.0):
    """Compute area of spherical polygon given list of (lat, lon) pairs."""
    if len(coords) != 4:
        raise ValueError("This function only works with 4 points.")
    
    # Convert lat/lon to 3D unit vectors
    pts = [to_unit_vector(lat, lon) for lat, lon in coords]
    
    # Split into two triangles: [0,1,2] and [0,2,3]
    area1 = spherical_triangle_area(pts[0], pts[1], pts[2])
    area2 = spherical_triangle_area(pts[0], pts[2], pts[3])
    
    # Total area on unit sphere
    total_area = area1 + area2
    
    # Scale by radius squared if needed (e.g. Earth: radius=6371e3 for meters)
    return np.abs(total_area * radius ** 2)

# get values from settingsfile
def parse_settings_file(filepath):
    settings = {}
    current_section = None
    with open(filepath, 'r') as f:
        for line in f:
            # Remove comments and whitespace
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Detect section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                settings[current_section] = {}
                continue

            # Parse key-value pairs
            if '=' in line and current_section:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                settings[current_section][key] = value

    return settings

# -- get data path from settings file ---
#starttime = time.time()
binding = parse_settings_file('settings_CCWatM_5min_example.ini')
#print(time.time()-starttime)
#print('Keys:',list(binding['COUPLING'].keys()))
meteoforc = MeteoForc2Var(binding['COUPLING']['PathForc'],binding['COUPLING']['fmodel_flag'])

# read monthly files
# TODO: reload every month
starttime = binding['TIME-RELATED_CONSTANTS']['StepStart']
ctime = datetime.datetime.strptime(starttime, '%d/%m/%Y')


# -- read grid data of REMO input file---
filepath = '/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/2000/'
ds = xr.open_dataset(filepath+'e100001n_c140_200001.nc')
lon_2d = ds['lon'].values.T
lat_2d = ds['lat'].values.T
nlon_forcing = lon_2d.shape[0]
nlat_forcing = lon_2d.shape[1]

# test with dummy regular grid:
latmin, latmax = 34.851, 71.185
lonmin, lonmax = -10.668, 41.29
deltalat = (latmax-latmin)/nlat_forcing
deltalon = (lonmax-lonmin)/nlon_forcing

# -> same as for ccwatm
#lon_1d = lonmin + deltalon*np.arange(nlon_forcing)
#lat_1d = latmax - deltalat*np.arange(nlat_forcing)
#lat_2d,lon_2d = np.meshgrid(lat_1d[::-1],lon_1d)


# --- derive grid corners ---

# TODO: this is a dummy for a regular grid
##grid_clon, grid_clat = derive_regular_grid_corners(lon_2d.T,lat_2d.T)
#grid_clon, grid_clat = derive_regular_grid_corners(lon_2d,lat_2d)


#lon_rot,lat_rot = np.meshgrid(ds['rlon'].values,ds['rlat'].values)
lat_rot,lon_rot = np.meshgrid(ds['rlat'].values,ds['rlon'].values)
rlon_corners,rlat_corners = calc_rotgrid_corners(lon_rot,lat_rot)

rot_pole_lon = ds.rotated_latitude_longitude.grid_north_pole_longitude
rot_pole_lat = ds.rotated_latitude_longitude.grid_north_pole_latitude
grid_clon = unrot_lon(rlat_corners, rlon_corners, rot_pole_lat, rot_pole_lon)
grid_clat = unrot_lat(rlat_corners, rlon_corners, rot_pole_lat, rot_pole_lon)


# get landmask from soil water data
landmask_input = ds['WS'][0,:,:].values
landmask_input[landmask_input>0] = 1 
# dummy
#landmask_input = np.ones(landmask_input.shape)
#landmask_input[:,:] = 0 

# TODO: put in function?
# --- 1) Initialization ---
### OASIS_INIT_COMP ###
comp = pyoasis.Component("forcing_component")

### get localcomm and define partition
partition, w_unit = oasis_specify_partition(comp,nlon=nlon_forcing,nlat=nlat_forcing)

# grid definition
print(f' grid_lon maximum and minimum', '%.5f' % np.max(lon_2d), '%.5f' % np.min(lon_2d), file=w_unit)
print(f' grid_lat maximum and minimum', '%.5f' % np.max(lat_2d), '%.5f' % np.min(lat_2d), file=w_unit)
w_unit.flush()
# function for writing oasis grid information
# test1:
#oasis_define_grid_simple(nlon_forcing,nlat_forcing,lon_2d.T,lat_2d.T,landmask_input.T[:,::-1],partition,'forcing_grid')

# with corners
#oasis_define_grid(nlon_forcing,nlat_forcing,lon_2d.T,lat_2d.T,grid_clon,grid_clat,landmask_input.T[:,::-1],partition,'forcing_grid')
#oasis_define_grid(nlon_forcing,nlat_forcing,lon_2d,lat_2d,grid_clon,grid_clat,landmask_input.T,partition,'forcing_grid')
oasis_define_grid(nlon_forcing,nlat_forcing,lon_2d,lat_2d,grid_clon,grid_clat,1-landmask_input.T,partition,'forcing_grid')

# Define the variable to send
#forcing_var = oasis.def_var("FORCING_FIELD", oasis.OUT)

#  DECLARATION OF THE COUPLING FIELDS
################## OASIS_DEF_VAR #################################
numcouple = 4 # number of coupling fields
var_id = [None] * numcouple
# TODO -> "FIELD_RECV_ATM" needs to be declared in namcouple file
var_id[0] = pyoasis.Var("FIELD_SEND_runoff", partition, OASIS.OUT)
var_id[1] = pyoasis.Var("FIELD_SEND_gwRecharge", partition, OASIS.OUT)
var_id[2] = pyoasis.Var("FIELD_SEND_EWRef", partition, OASIS.OUT)
var_id[3] = pyoasis.Var("FIELD_SEND_rootzoneSM", partition, OASIS.OUT)
print(f' var_id FRECVATM, {var_id[0]._id}', file=w_unit)
w_unit.flush()

### TERMINATION OF DEFINITION PHASE ###
print(' End of initialisation phase', file=w_unit)
w_unit.flush()

### OASIS_ENDDEF ###
comp.enddef()

# -------------- get and put -----------------

# Load forcing data (e.g., NetCDF)
# maybe read from settings file??
#filepath = '/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/2000/'
#ds = xr.open_dataset(filepath+'e100001n_c160_200001.nc')
#data_array = ds['RUNOFF']  # Replace with your variable name

ds = xr.open_dataset('/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/e100001n_c105.nc')
FCAP = np.squeeze(ds['FCAP'].values)
ds = xr.open_dataset('/work/ch0636/projects/uwares/CWatM_forcing/Remo_ERA5_27lev/daily_means/e100001n_c229.nc')
WSMX = np.squeeze(ds['WSMX'].values)
#print(lat.shape)

# Time loop
#nt = data_array.sizes['time']

for t in range(5): # same number of time loops as C-CWatM
    seconds_passed = int(t * 86400.)
    print('dummy',seconds_passed)

    currenttime = ctime + datetime.timedelta(seconds=seconds_passed)
    meteoforc.read_forcing(currenttime,'runoff',binding['COUPLING']['RunoffName'])
    meteoforc.read_forcing(currenttime,'sum_gwRecharge',binding['COUPLING']['GWName'])
    meteoforc.read_forcing(currenttime,'EWRef',binding['COUPLING']['OWEName'])
    meteoforc.read_forcing(currenttime,'rootzoneSM',binding['COUPLING']['SMName'])
    # note: this is in percent, 
    # needs to be multiplied with soilWaterStorageCap later
    meteoforc.rootzoneSM = meteoforc.rootzoneSM * FCAP / WSMX
    
    # Extract data for current timestep
    #data_t = data_array.isel(time=t).values #/ 100.
    #data_t = meteoforc.runoff.values

    # send (and get data), see atmos.py example
    # include dummy send and put also in ccwatm
    #var_id[0].put(seconds_passed, data_t[::-1,:])
    #var_id[0].put(seconds_passed, data_t.T[:,::-1])
    var_id[0].put(seconds_passed, meteoforc.runoff.values/100.)
    var_id[1].put(seconds_passed, meteoforc.sum_gwRecharge.values/1000.)
    var_id[2].put(seconds_passed, meteoforc.EWRef.values)
    var_id[3].put(seconds_passed, meteoforc.rootzoneSM.values)

# Finalize
del comp

