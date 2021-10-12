#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 12:28:13 2021

@author: lucas
@contact: lgvivanco96@gmail.com

# =============================================================================
# Script for clipping a netCDF multidimensional raster with a polygon inside an
# ESRI shapefile. 
# Example use: clip_ncWshp.py polygon.shp netCDF.nc output.nc
# Only valid for EPSG:4326 projection files. netCDF coordinates names must be 
# lat,lon
# =============================================================================

"""

#%%
# =============================================================================
# Import libraries and install if not installed.
# =============================================================================
import os
import sys

try:
    import xarray as xr
except ModuleNotFoundError as exception:
    print(exception)
    print("Module xarray not installed...")
    user = input("Should Install? (y/n): ")
    if user=="y":
        print("Installing...")
        os.system("conda install -c conda-forge -y xarray")
    else:
        print("Done")
        exit()
try:
    import geopandas as gpd
except ModuleNotFoundError as exception:
    print(exception)
    print("Module geopandas not installed...")
    user = input("Should Install? (y/n): ")
    if user=="y":
        print("Installing...")
        os.system("conda install -c conda-forge -y geopandas")
    else:
        print("Done")
        exit()


try:
    import regionmask
except ModuleNotFoundError as exception:
    print(exception)
    print("Module regionmask not installed...")
    user = input("Should Install? (y/n): ")
    if user=="y":
        print("Installing...")
        os.system("conda install -c conda-forge -y regionmask")
    else:
        print("Done")
        exit()
#%%
# =============================================================================
# Load data
# =============================================================================

print("Loading files...")

path_shp = sys.argv[1]
path_nc  = sys.argv[2]


vector   = gpd.read_file(path_shp)
raster   = xr.open_dataset(path_nc)

try:
    raster.crs
except:
    raster = raster.rio.write_crs("EPSG:4326")
    
try:
    vector.crs
except:
    vector = vector.to_crs(epsg=4326)
    
    
#%%
# =============================================================================
# Clip and save file
# =============================================================================

print("Creating shapefile mask...")
lat,lon = raster.lat.values,raster.lon.values
mask = regionmask.mask_geopandas(vector,lon,lat)
print("Clipping...")
clipped = raster.where(mask==0)
clipped.to_netcdf(sys.argv[3])

print("Done")
