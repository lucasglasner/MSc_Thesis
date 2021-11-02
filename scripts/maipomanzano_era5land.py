#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:42:55 2021

@author: lucas

# =============================================================================
# Maipo en el Manzano Basin. ERA5 Snow analisis and comparison with observations.
# =============================================================================

"""


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import regionmask
import geopandas as gpd
import matplotlib as mpl

#%%

#ERA5LAND topography on MaipoEnElManzano basin.
topo  = xr.open_dataset("datos/topography/basins/RioMaipoEnElManzano_ERA5LAND.nc")
topo  = topo.z.squeeze().assign_coords({"lon": (((topo.lon + 180) % 360) - 180)})

#BASIN shapefile
basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")

#Basin hypsometry curve
dem = xr.open_dataset("datos/topography/basins/RioMaipoEnElManzano_regridmodis.nc",chunks="auto").Band1
#Elevation band width in meters
dz    = 50
tpix  = (~np.isnan(dem.values)).sum()
tpix2 = (~np.isnan(topo.values)).sum()
#elevation bands
elevation_bands = np.arange(dem.min()//10*10-2*dz,
                            dem.max()//10*10+2*dz,dz,
                            dtype=np.int64)

#bands number of pixels
bands_size      = np.empty(len(elevation_bands))*np.nan
#hypsometric curve
hypso_curve          = np.empty(len(elevation_bands))*np.nan
hypso_curve_ERA5     = np.empty(len(elevation_bands))*np.nan
#hypsometric density (band number of pixels / total pixels)
hypso_density        = np.empty(len(elevation_bands))*np.nan
hypso_density_ERA5   = np.empty(len(elevation_bands))*np.nan
for i,height in enumerate(elevation_bands):
    mask  = (~np.isnan(dem.where(dem>height).where(dem<height+dz)))
    mask2 = (~np.isnan(topo.where(topo>height).where(topo<height+dz)))
    bands_size[i]    = mask.sum()
    hypso_density[i] = bands_size[i]/tpix
    hypso_density_ERA5[i] = mask2.sum()/tpix2
hypso_curve = np.cumsum(hypso_density)
hypso_curve_ERA5 = np.cumsum(hypso_density_ERA5)
#%%

SCA_ianigla = pd.read_csv("datos/ianigla/RioMaipoEnElManzano_SCA_s_comp.filtro_MA.3días.csv",index_col=0)
SCA_ianigla = SCA_ianigla.iloc[:,0]
SCA_ianigla.index = pd.to_datetime(SCA_ianigla.index)

#%%

SCOV_ERA5    = xr.open_dataset("datos/era5land/maipo/snow_cover.nc",chunks="auto")
SCOV_ERA5    = SCOV_ERA5.snowc[23:,:,:][::12,:,:][::2,:,:]
SCA_ERA5     = pd.Series(np.empty(len(SCOV_ERA5.time))*np.nan,index=SCOV_ERA5.time.values)

for i,time in enumerate(SCOV_ERA5.time.values):
    sca = (SCOV_ERA5[i,:,:]>95).values.sum()/49
    SCA_ERA5[time] = sca
    
SCA_ERA5 = SCA_ERA5.resample("d").mean()

#%%
SCOV_ERA5    = xr.open_dataset("datos/era5land/maipo/snow_depth_water_equivalent.nc",chunks="auto")
# SCOV_ERA5    = SCOV_ERA5.sd[23:,:,:][::12,:,:][::2,:,:]
SCOV_ERA5    = SCOV_ERA5.sd
SCA_ERA5     = pd.Series(np.empty(len(SCOV_ERA5.time))*np.nan,index=SCOV_ERA5.time.values)

for i,time in enumerate(SCOV_ERA5.time.values):
    sca = (SCOV_ERA5[i,:,:]>10*1e-3).values.sum()/49
    SCA_ERA5[time] = sca
    
SCA_ERA5 = SCA_ERA5.resample("d").mean()


#%%
tot = pd.concat([SCA_ERA5.reindex(SCA_ianigla.index),SCA_ianigla/100],axis=1)


cmaplist = pd.read_csv("terraincolormap.txt").values
cmaplist = [list(np.array(i)/255) for i in cmaplist]
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, len(cmaplist))
fig,ax=plt.subplots()
basin.boundary.plot(ax=ax,zorder=2,color="k")
topo.plot(ax=ax,cmap=cmap)
ax.axis("off")