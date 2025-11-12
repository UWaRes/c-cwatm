# Name:        grid_tools module
# Purpose:
#
# Author:      Amelie Schmitt
#
# Created:     24/07/2025
# Copyright:   (c) Amelie Schmitt 2025 
# -------------------------------------------------------------------------

import numpy as np

class grid_tools:
    """
    Functions related to the model grid that are required for the OASIS grids.nc file. 

    **Example usage for C-CWatM**

    ccgrid = grid_tools()
    lon_2d,lat_2d = ccgrid.create_2d_ccwatm_grid()
    lon_corners, lat_corners = grid_tools.compute_grid_corners(lon_2d,lat_2d)
    areas = compute_grid_cell_areas(lon_corners, lat_corners)
    
    """
    def __init__(self, maskmapAttr=None):
        self.maskmapAttr = maskmapAttr

    def create_2d_ccwatm_grid(self):
        """
        Generate 2D arrays of longitude and latitude coordinates for a C-CWatM regular grid.
    
        Latitude values are arranged from north to south and longitude values from west to east.
    
        Returns
        -------
        lon_2d (ndarray) : 2D array of longitudes at the center of each grid cell.
        lat_2d (ndarray) : 2D array of latitudes at the center of each grid cell.

        Uses global variables
        ---------------------
        maskmapAttr['x'] : x-coordinate (longitude) of upper left corner cell.
        maskmapAttr['y'] : y-coordinate (latitude) of upper left corner cell.
        maskmapAttr['cell'] : Cell length in both x and y directions.
        maskmapAttr['col'] : Number of columns (longitude).
        maskmapAttr['row'] : Number of rows (latitude).

        as/copilot
        """
        lon_1d = self.maskmapAttr['x'] + self.maskmapAttr['cell'] * np.arange(self.maskmapAttr['col'])
        lat_1d = self.maskmapAttr['y'] - self.maskmapAttr['cell'] * np.arange(self.maskmapAttr['row'])
        lat_2d, lon_2d = np.meshgrid(lat_1d[::-1], lon_1d)
        
        return lon_2d, lat_2d

    @staticmethod
    def compute_grid_corners(lon, lat):
        """
        Compute the corner coordinates of grid cells from center coordinates.

        The corners are ordered counter-clockwise, as required by OASIS.
    
        Parameters
        ----------
        lon (ndarray) : 2D array of longitudes at the center of each grid cell.
        lat (ndarray) : 2D array of latitudes at the center of each grid cell.

        Returns
        -------
        lon_corners (ndarray) : 3D array of longitudes at the corners of each grid cell.
        lat_corners (ndarray) : 3D array of latitudes at the corners of each grid cell.

        as/copilot
        """
        # grid spacing
        dx = np.abs(lon[1, 1] - lon[0, 0])
        dy = np.abs(lat[1, 1] - lat[0, 0])
    
        lon_corners = np.zeros((*lon.shape, 4))
        lat_corners = np.zeros((*lat.shape, 4))
    
        lon_corners[:, :, 0] = lon - dx / 2.
        lon_corners[:, :, 1] = lon + dx / 2.
        lon_corners[:, :, 2] = lon + dx / 2.
        lon_corners[:, :, 3] = lon - dx / 2.
    
        lat_corners[:, :, 0] = lat - dy / 2.
        lat_corners[:, :, 1] = lat - dy / 2.
        lat_corners[:, :, 2] = lat + dy / 2.
        lat_corners[:, :, 3] = lat + dy / 2.
    
        return lon_corners, lat_corners

    @staticmethod
    def compute_grid_cell_areas(lon_corners, lat_corners):
        """
        Compute the area of each grid cell on a spherical Earth using corner coordinates.

        Formulas based on Girard's theorem.
    
        Parameters
        ----------
        lat_corners (ndarray) : 3D array of latitudes at the corners of each grid cell (shape: [rows, cols, 4]).
        lon_corners (ndarray) : 3D array of longitudes at the corners of each grid cell (shape: [rows, cols, 4]).
    
        Returns
        -------
        areas (ndarray) : 2D array of grid cell areas in square meters (shape: [rows, cols]).
    
        as/copilot
        """
        def to_unit_vector(lat, lon):
            """Convert lat/lon to 3D unit vectors on a sphere"""
            lat_rad = np.radians(lat)
            lon_rad = np.radians(lon)
            x = np.cos(lat_rad) * np.cos(lon_rad)
            y = np.cos(lat_rad) * np.sin(lon_rad)
            z = np.sin(lat_rad)
            return np.stack([x, y, z], axis=-1)
    
        def angle_between(a, b, c):
            """Compute the spherical angle at point b between vectors a and c."""
            ab = np.cross(b, a)
            cb = np.cross(b, c)
            ab /= np.linalg.norm(ab, axis=-1, keepdims=True)
            cb /= np.linalg.norm(cb, axis=-1, keepdims=True)
            dot_product = np.sum(ab * cb, axis=-1)
            return np.arccos(np.clip(dot_product, -1.0, 1.0))
    
        def spherical_triangle_area(a, b, c):
            """Calculate area of a spherical triangle using spherical excess."""
            A = angle_between(b, a, c)
            B = angle_between(c, b, a)
            C = angle_between(a, c, b)
            return A + B + C - np.pi
    
        def earth_radius_at_latitude(latitude_deg):
            """Compute Earth's radius in m at a given latitude using the WGS-84 ellipsoid."""
            lat_rad = np.radians(latitude_deg)
            a = 6378137.
            b = 6356752.
            numerator = (a**2 * np.cos(lat_rad))**2 + (b**2 * np.sin(lat_rad))**2
            denominator = (a * np.cos(lat_rad))**2 + (b * np.sin(lat_rad))**2
            return np.sqrt(numerator / denominator)
    
        # Convert corners to unit vectors
        pts = to_unit_vector(lat_corners, lon_corners)
    
        # Split each cell into two triangles: [0,1,2] and [0,2,3]
        area1 = spherical_triangle_area(pts[:, :, 0], pts[:, :, 1], pts[:, :, 2])
        area2 = spherical_triangle_area(pts[:, :, 0], pts[:, :, 2], pts[:, :, 3])
        total_area_unit_sphere = area1 + area2
    
        # Compute Earth radius at cell center (mean of corner latitudes)
        mean_lat = np.mean(lat_corners, axis=-1)
        radius = earth_radius_at_latitude(mean_lat)
    
        # Scale area by radius squared
        areas = np.abs(total_area_unit_sphere * radius**2)
    
        return areas

    @staticmethod
    def unrot_coordinates(rlon, rlat, pole_lon, pole_lat):
        """
        Convert rotated latitude and longitude coordinates to unrotated (geographic) coordinates.

        This function transforms coordinates from a rotated pole grid system to standard 
        geographic coordinates (latitude and longitude in degrees), based on the position 
        of the rotated pole. Used e.g. for REMO forcing data.

        Parameters
        ----------
        rlon (ndarray) : Rotated longitude(s) in degrees.
        rlat (ndarray) : Rotated latitude(s) in degrees.
        pole_lon (float) : Longitude of the rotated pole in degrees.
        pole_lat (float) : Latitude of the rotated pole in degrees.

        Returns
        -------
        unrot_lon (ndarray) : Unrotated (geographic) longitudes in degrees.
        unrot_lat (ndarray) : Unrotated (geographic) latitudes in degrees.
        
        as/copilot
        """
        # convert input to radians
        rlat_rad = np.deg2rad(rlat)
        rlon_rad = np.deg2rad(rlon)
        pole_lat_rad = np.deg2rad(pole_lat)
        pole_lon_rad = np.deg2rad(pole_lon)
        
        # Precompute sine and cosine of coordinates
        sin_pole_lat = np.sin(pole_lat_rad)
        cos_pole_lat = np.cos(pole_lat_rad)
        sin_pole_lon = np.sin(pole_lon_rad)
        cos_pole_lon = np.cos(pole_lon_rad)
        cos_rlon = np.cos(rlon_rad)
        sin_rlon = np.sin(rlon_rad)
        cos_rlat = np.cos(rlat_rad)
        sin_rlat = np.sin(rlat_rad)
        
        # Compute intermediate values
        x = -sin_pole_lat * cos_rlon * cos_rlat + cos_pole_lat * sin_rlat
        tmp1 = sin_pole_lon * x - cos_pole_lon * sin_rlon * cos_rlat
        tmp2 = cos_pole_lon * x + sin_pole_lon * sin_rlon * cos_rlat
        temp3 = sin_pole_lat*sin_rlat+cos_pole_lat*cos_rlat*cos_rlon
        unrot_lon = np.rad2deg(np.arctan2(tmp1, tmp2))
        unrot_lat = np.rad2deg(np.arcsin(temp3))
        
        # Return unrotated longitude and latitude in degrees
        return unrot_lon, unrot_lat

