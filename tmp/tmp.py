#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:29:04 2022

@author: lucas

# =============================================================================
# This Script is going to create a tabular dataset with ROS events and their
# caracteristics (land cover, % area, related basin, runoff, precipitation, etc)
# =============================================================================
"""

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
import regionmask
import geopandas as gpd
from glob import glob
# %%
# Target Year

yr = "2015"

print('Computing year number: '+yr)
# =============================================================================
# ROS dataset
# =============================================================================
ROS = xr.open_dataset('datos/ROS/CORTES_CR2MET_ERA5/ROS_'+yr+'.nc').ROS
# ROS = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*')).ROS

# =============================================================================
# Precipitation dataset
# =============================================================================

PR = xr.open_dataset('datos/cr2met/CR2MET_pr_'+yr+'.nc').pr
# PR = xr.open_mfdataset(glob('datos/cr2met/CR2MET_pr_*')).pr
# =============================================================================
# SWE dataset
# =============================================================================

SWE = xr.open_dataset(
    'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+yr+'.nc').SWE

# SWE = xr.open_mfdataset(
#     glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_*')).SWE

# =============================================================================
# SWE loss dataset
# =============================================================================
dSWE = xr.open_dataset(
    'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_'+yr+'.nc').SWE
# dSWE = xr.open_mfdataset(
#     glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_*')).SWE
dSWE.name = 'dSWE'

# =============================================================================
# freezing level dataset
# =============================================================================

FL = xr.open_dataset('datos/era5/H0_ERA5_'+yr+'.nc').deg0l
# FL = xr.open_mfdataset(glob('datos/era5/H0_ERA5_*.nc')).deg0l
FL = FL.reindex({'time': SWE.time}, method='nearest')
FL.name = 'FL'

# =============================================================================
# Land Cover dataset
# =============================================================================

LC = xr.open_dataset('datos/landcover/LC_CHILE_2014_b_final_regionaldomain.nc')
LC = LC.Band1.to_dataset(name='Land_Cover').Land_Cover

# %%
basins = ["Rio Aconcagua En Chacabuquito",
          "Rio Maipo En El Manzano",
          "Rio Colorado En Junta Con Palos",
          "Rio Putaendo En Resguardo Los Patos",
          "Rio Teno Despues De Junta Con Claro",
          "Rio Tinguiririca Bajo Los Briones",
          "Rio ÑUble En San Fabian"]


# %%


def grab_basinpolygon(basin_name):
    try:
        basin_name = basin_name.replace(' ', '')
        basins = gpd.read_file('datos/vector/'+basin_name+'.shp')
    except:
        basins = gpd.read_file('datos/vector/cuencas_CAMELS.gpkg')
        mask = basins["gauge_name"] == basin_name
        basins = basins[mask]
    return basins


def grab_basinhypso(basin_name):
    basin_name = basin_name.replace(' ', '')
    try:
        hypso = pd.read_csv('datos/topography/basins/hypso/' +
                            basin_name+'_hypso.csv')
    except:
        raise ValueError('Basin hypsometry not created')
    return hypso


def clip_nc(raster, vector):
    extent = vector.bounds.values[0]
    lons = sorted((extent[0], extent[2]))
    lats = sorted((extent[1], extent[3]))
    clipped = raster.sortby("lat").sortby("lon")
    clipped = clipped.sel(lon=slice(*lons), lat=slice(*lats))

    lat, lon = clipped.lat.values, clipped.lon.values
    mask = regionmask.mask_geopandas(vector, lon, lat)

    clipped = clipped.where(~np.isnan(mask))
    return clipped


def join_basindatasets(basin_name, list_arrays):
    basin_data = []
    vector = grab_basinpolygon(basin_name)
    for var in list_arrays:
        clipped = clip_nc(var, vector)
        basin_data.append(clipped)
    raster = xr.merge(basin_data, join='inner')
    raster.attrs = {'basin': basin_name}
    return raster


def compute_basin_SCA(data, datatype='ROS'):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    if datatype == 'basinwide':
        SCA = np.where(data.SWE > 10, 1, 0)
    elif datatype == 'ROS':
        SCA = np.where((data.SWE > 10) & (data.ROS == True), 1, 0)
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    SCA = SCA.sum(axis=1).sum(axis=1)/npixels
    SCA = np.clip(SCA, 0, 1)
    return pd.Series(SCA, index=data.time)


def compute_basin_RCA(data, datatype='ROS'):
    if datatype in ['ROS', 'basinwide']:
        npixels = ~np.isnan(data.Land_Cover)
        npixels = np.count_nonzero(npixels)
        RCA = np.where(data.ROS == True, 1, 0)
        RCA = RCA.sum(axis=1).sum(axis=1)/npixels
        RCA = np.clip(RCA, 0, 1)
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    return pd.Series(RCA, index=data.time)


def compute_basin_dSWE(data, datatype='ROS'):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    if datatype == 'ROS':
        dSWE = data.dSWE.where(data.ROS == True)
    elif datatype == 'basinwide':
        dSWE = data.dSWE
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    dSWE = dSWE.sum(dim='lat').sum(dim='lon')
    dSWE = dSWE/npixels
    return dSWE.to_series()


def compute_basin_SWE(data, datatype='ROS'):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    if datatype == 'ROS':
        SWE = data.SWE.where(data.ROS == True)
    elif datatype == 'basinwide':
        SWE = data.SWE
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    SWE = SWE.sum(dim='lat').sum(dim='lon')
    SWE = SWE/npixels
    return SWE.to_series()


def compute_basin_meanpr(data, datatype='ROS'):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    if datatype == 'ROS':
        PR = data.pr.where(data.ROS == True)
    elif datatype == 'basinwide':
        PR = data.pr
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    PR = PR.sum(dim='lat').sum(dim='lon')
    PR = PR/npixels
    return PR.to_series()


def compute_basin_pluvialarea(data, datatype='ROS'):
    if datatype in ['ROS', 'basinwide']:
        npixels = ~np.isnan(data.Land_Cover)
        npixels = np.count_nonzero(npixels)
        FL = data.FL
        PR = data.pr
        hypso = grab_basinhypso(data.attrs['basin'])
        total_area = hypso.iloc[:, 2].max()
        Ap = np.where((FL > 300) & (PR > 0), 1, 0)
        Ap = Ap.sum(axis=1).sum(axis=1)/npixels
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    return pd.Series(Ap*total_area, data.time)


def compute_ROS_LC(data, datatype='ROS'):
    land_uses = {"Cropland": (100, 200),
                 "Native Forest": (200, 229),
                 "Forest Plantation": (229, 300),
                 "Grassland": (300, 400),
                 "Shrubland": (400, 500),
                 "Wetland": (500, 600),
                 "Water Bodies": (600, 800),
                 "Waterproofs": (800, 900),
                 "Barren": (900, 1000),
                 "Ice/Snow": (1000, 1200),
                 "Clouds": (1200, 1500)}
    return


def compute_dataset(data, datatype='basinwide'):
    if datatype == 'basinwide':
        functions = ['SCA', 'RCA', 'dSWE', 'SWE', 'meanpr', 'pluvialarea']
        functions = [eval('compute_basin_'+x) for x in functions]
        series = [f(data, datatype) for f in functions]
        series = pd.concat(series, axis=1)
        names = ['SCA', 'RCA', 'dSWE', 'SWE', 'PR', 'PR_Area']
        series.columns = names
    elif datatype == 'ROS':
        functions = ['RCA', 'SCA', 'dSWE', 'SWE', 'meanpr', 'pluvialarea']
        functions = [eval('compute_basin_'+x) for x in functions]
        series = [f(data, datatype) for f in functions]
        series = pd.concat(series, axis=1)
        names = ['ROS_SCA', 'RCA', 'ROS_dSWE', 'ROS_SWE', 'ROS_PR', 'PR_Area']
        series.columns = names
    else:
        error = datatype+' Wrong Input: It can only be "ROS" or "basinwide"'
        raise ValueError(error)
    return series


# %%
print('Creating raster data for each basin...')

raster_data = [join_basindatasets(b, [SWE, PR, ROS, LC, dSWE, FL])
               for b in basins]

# %%
print('Computing time series...')
data = [compute_dataset(r, 'basinwide') for r in raster_data]
data = pd.concat(data, keys=basins)

data_ROS = [compute_dataset(r, 'ROS') for r in raster_data]
data_ROS = pd.concat(data_ROS, keys=basins)

# %%
print('Saving...')
data.to_csv('datos/ROS/Basinwide_Timeseries_'+yr+'.csv')
data.to_csv('datos/ROS/ROS_Timeseries_'+yr+'.csv')
print('Done')
