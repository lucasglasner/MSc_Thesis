#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 10:12:44 2022

@author: lucas

# =============================================================================
# This Script does the spatial analysis of ROS days in the Andes mountain range
# =============================================================================
"""

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.feature as cf
import matplotlib.ticker as mticker
import geopandas as gpd
from glob import glob
import matplotlib as mpl
# %%
# =============================================================================
# Load Basin shapefiles
# =============================================================================
paths = glob('datos/vector/Rio*.shp')
cuencas = []
for path in paths:
    cuencas.append(gpd.read_file(path))
cuencas = pd.concat(cuencas)

# =============================================================================
# Load CR2MET precipitation and compute the mean rainy days per year
# =============================================================================
PR_CR2MET = xr.open_mfdataset('datos/cr2met/CR2MET_pr*').pr
PR_CR2MET = xr.where(PR_CR2MET > 3, 1, 0)
mean_rainydays = PR_CR2MET.resample(
    {"time": "y"}).sum().mean(dim='time').load()
del PR_CR2MET
# %%

# =============================================================================
# Load ROS dataset and compute anual cycle, annual frequency maps and ROS
# maximum timing (month)
# =============================================================================
ROS = xr.open_mfdataset('datos/ROS/CORTES_CR2MET_ERA5/ROS*').ROS
lat, lon = ROS.lat.values, ROS.lon.values
LON, LAT = np.meshgrid(lon, lat)
freq_CORTESCR2MET = ROS.resample({'time': 'y'}).sum()
freq_CORTESCR2MET = freq_CORTESCR2MET.mean(dim='time').load()
cycle_CORTESCR2MET = ROS.groupby('time.month').mean().load()
timing = cycle_CORTESCR2MET.to_dataframe()
timing = timing.groupby(["lat", "lon"]).idxmax().ROS
timing = timing.map(lambda x: x[0] if type(x) == tuple else x).unstack()

# %%
# =============================================================================
# Build a mask with snow areas (SWE_year_max>20mm)
# =============================================================================
maxSWE = []
for yr in range(1984, 2016, 1):
    path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+str(yr)+'.nc'
    SWE = xr.open_dataset(path).SWE.max(dim='time')
    maxSWE.append(SWE)
mask = xr.concat(maxSWE, 'max_years').mean(dim='max_years') > 20
del maxSWE

# %%
# =============================================================================
# Create Data for maps
# =============================================================================

ROS_days = freq_CORTESCR2MET.where(mask)
mean_rainydays = mean_rainydays.where(mask).reindex({'lat': ROS_days.lat,
                                                     'lon': ROS_days.lon},
                                                    method='nearest')
timing = timing.where(mask.values)


# %%
fig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
                       subplot_kw={"projection": ccrs.PlateCarree()},
                       figsize=(18, 6))
ax2 = []
plt.rcParams.update({'font.size': 18})
ax = ax.ravel()
for axis in ax:
    axis.set_extent([-74, -68, -26, -38])
    # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
    axis.coastlines(rasterized=True)
    axis.add_feature(cf.BORDERS)
    gl = axis.gridlines(linestyle=":")
    gl.xlocator = mticker.FixedLocator([])
    gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = False
    cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
    axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
    axis.add_feature(cf.OCEAN, rasterized=True)


mask = (ROS_days/mean_rainydays != 0).values
mapa0 = ax[0].pcolormesh(LON, LAT, mean_rainydays,
                         cmap="cividis_r",
                         transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa1 = ax[1].pcolormesh(LON, LAT, ROS_days,
                         cmap="viridis",
                         transform=ccrs.PlateCarree(), vmin=0, vmax=30,
                         rasterized=True)
mapa2 = ax[2].pcolormesh(LON, LAT,
                         (ROS_days/mean_rainydays),
                         cmap="BuPu", transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa3 = ax[3].pcolormesh(LON, LAT, timing.where(mask),
                         cmap='nipy_spectral',
                         transform=ccrs.PlateCarree(), rasterized=True)


ax[0].set_title("Rain Frequency\n"+r"($\frac{N°Rainy Days}{year}$)")
ax[1].set_title("ROS Frequency\n"+r"($\frac{N°ROS Days}{year}$)")
ax[2].set_title("ROS/Rain\nFrequency\nRatio (%)")
ax[3].set_title("Maximum ROS\nTiming")

cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40)
cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40)
cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40)
cb2 = fig.colorbar(mapa3, ax=ax[3], aspect=40, ticks=np.arange(1, 13, 1))
cb2.set_ticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                    "JUL", "AUG", "SEP", "OCT", "NOV", "DIC"])

gl = ax[0].gridlines(draw_labels=True, linestyle=":")
gl.xlocator = mticker.FixedLocator([])
gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
gl.top_labels = False
gl.right_labels = False
gl.xlines = False
gl.ylines = False
plt.savefig('plots/ROS_CORTESCR2METERA5_73W-68W-26S-38S_new.pdf', dpi=150,
            bbox_inches="tight")

# %%

fig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
                       subplot_kw={"projection": ccrs.PlateCarree()},
                       figsize=(18, 6))
ax2 = []
plt.rcParams.update({'font.size': 18})
ax = ax.ravel()
for axis in ax:
    axis.set_extent([-72, -69, -32.4, -37.3])
    # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
    axis.coastlines(rasterized=True)
    axis.add_feature(cf.BORDERS)
    gl = axis.gridlines(linestyle=":")
    gl.xlocator = mticker.FixedLocator([])
    gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = False
    cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
    axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
    axis.add_feature(cf.OCEAN, rasterized=True)


mask = (ROS_days/mean_rainydays != 0).values
mapa0 = ax[0].pcolormesh(LON, LAT, mean_rainydays,
                         cmap="cividis_r",
                         transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa1 = ax[1].pcolormesh(LON, LAT, ROS_days,
                         cmap="viridis",
                         transform=ccrs.PlateCarree(), vmin=0, vmax=30,
                         rasterized=True)
mapa2 = ax[2].pcolormesh(LON, LAT,
                         (ROS_days/mean_rainydays),
                         cmap="BuPu", transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa3 = ax[3].pcolormesh(LON, LAT, timing.where(mask),
                         cmap='nipy_spectral',
                         transform=ccrs.PlateCarree(), rasterized=True)


ax[0].set_title("Rain Frequency\n"+r"($\frac{N°Rainy Days}{year}$)")
ax[1].set_title("ROS Frequency\n"+r"($\frac{N°ROS Days}{year}$)")
ax[2].set_title("ROS/Rain\nFrequency\nRatio (%)")
ax[3].set_title("Maximum ROS\nTiming")

cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40)
cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40)
cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40)
cb2 = fig.colorbar(mapa3, ax=ax[3], aspect=40, ticks=np.arange(1, 13, 1))
cb2.set_ticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                    "JUL", "AUG", "SEP", "OCT", "NOV", "DIC"])

gl = ax[0].gridlines(draw_labels=True, linestyle=":")
gl.xlocator = mticker.FixedLocator([])
gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
gl.top_labels = False
gl.right_labels = False
gl.xlines = False
gl.ylines = False
plt.savefig('plots/ROS_CORTESCR2METERA5_72W-69W-32S-37S_new.pdf', dpi=150,
            bbox_inches="tight")
