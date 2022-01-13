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

yr = "2013"

# =============================================================================
# ROS dataset
# =============================================================================
# ROS = xr.open_dataset('datos/ROS/CORTES_CR2MET_ERA5/ROS_'+yr+'.nc').ROS
ROS = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*')).ROS

# =============================================================================
# Precipitation dataset
# =============================================================================

# PR = xr.open_dataset('datos/cr2met/CR2MET_pr_'+yr+'.nc').pr
PR = xr.open_mfdataset(glob('datos/cr2met/CR2MET_pr_*')).pr
# =============================================================================
# SWE dataset
# =============================================================================

# SWE = xr.open_dataset(
#     'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+yr+'.nc').SWE

SWE = xr.open_mfdataset(
    glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_*')).SWE

# =============================================================================
# SWE loss dataset
# =============================================================================
# dSWE = xr.open_dataset(
# 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_'+yr+'.nc').SWE
dSWE = xr.open_mfdataset(
    glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_*')).SWE
dSWE.name = 'dSWE'

# =============================================================================
# freezing level dataset
# =============================================================================

# FL = xr.open_dataset('datos/era5/H0_ERA5_'+yr+'.nc').deg0l
FL = xr.open_mfdataset(glob('datos/era5/H0_ERA5_*.nc')).deg0l
FL = FL.reindex({'time': SWE.time}, method='nearest')
FL.name = 'FL'

# =============================================================================
# Land Cover dataset
# =============================================================================

LC = xr.open_dataset('datos/landcover/LC_CHILE_2014_b_final_regionaldomain.nc')
LC = LC.Band1.to_dataset(name='Land_Cover').Land_Cover

# %%
basins = ["Rio Aconcagua En Chacabuquito",
          "Rio Choapa En Lamahuida",
          "Rio Elqui En Algarrobal",
          "Rio Illapel En Huintil",
          "Rio Grande En Puntilla San Juan",
          "Rio Hurtado En Angostura De Pangue",
          "Rio Putaendo En Resguardo Los Patos",
          "Rio Mapocho En Los Almendros",
          "Rio Maipo En El Manzano",
          "Rio Cachapoal En Pte Termas De Cauquenes",
          "Rio Tinguiririca Bajo Los Briones",
          "Rio Teno Despues De Junta Con Claro",
          "Rio Colorado En Junta Con Palos",
          "Rio Maule En Armerillo",
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
        hypso = pd.read_csv('datos/topography/basins/hipso/' +
                            basin_name+'_Hipso.csv')
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
    return xr.merge(basin_data, join='inner')


def compute_basin_SCA(data):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    SCA = np.where(data.SWE > 10, 1, 0)
    SCA = SCA.sum(axis=1).sum(axis=1)/npixels
    SCA = np.clip(SCA, 0, 1)
    return pd.Series(SCA, index=data.time)


def compute_basin_RCA(data):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    RCA = np.where(data.ROS == True, 1, 0)
    RCA = RCA.sum(axis=1).sum(axis=1)/npixels
    RCA = np.clip(RCA, 0, 1)
    return pd.Series(RCA, index=data.time)


def compute_basin_dSWE(data):
    dSWE = data.dSWE.mean(dim='lat').mean(dim='lon')
    return dSWE.to_series()


def compute_basin_SWE(data):
    SWE = data.SWE.mean(dim='lat').mean(dim='lon')
    return SWE.to_series()


def compute_basin_meanpr(data):
    PR = data.pr.mean(dim='lat').mean(dim='lon')
    return PR.to_series()


def compute_basin_pluvialarea(data):
    npixels = ~np.isnan(data.Land_Cover)
    npixels = np.count_nonzero(npixels)
    FL = data.FL
    PR = data.pr
    Ap = np.where((FL > 300) & (PR > 0), 1, 0)
    Ap = Ap.sum(axis=1).sum(axis=1)/npixels
    Ap = np.clip(Ap, 0, 1)
    return pd.Series(Ap, data.time)


def compute_dataset(data):
    names = ['SCA', 'RCA', 'dSWE', 'SWE', 'meanpr', 'pluvialarea']
    functions = [eval('compute_basin_'+x) for x in names]
    series = [f(data) for f in functions]
    series = pd.concat(series, axis=1)
    series.columns = names
    return series


# %%
raster_data = join_basindatasets(basins[-2], [SWE, PR, ROS, LC, dSWE, FL])

data = compute_dataset(raster_data)
