#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:12:51 2021

@author: lucas
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
import gc
# %%
# =============================================================================
# Load Basin shapefiles
# =============================================================================
paths = glob('datos/vector/Rio*.shp')
cuencas = []
for path in paths:
    cuencas.append(gpd.read_file(path))
cuencas = pd.concat(cuencas)

# %%
# =============================================================================
# Load ERA5LAND precipitation and compute the mean rainy days per year
# =============================================================================
PR_ERA5LAND = xr.open_mfdataset('datos/era5land/total_precipitation*')
PR_ERA5LAND = PR_ERA5LAND.resample({'time': 'd'}).sum()
mean_rainydays = xr.where(PR_ERA5LAND > 3/1e3, 1, 0)
mean_rainydays = mean_rainydays.resample(
    {"time": "y"}).sum().mean(dim='time').load()
# del PR_ERA5LAND
gc.collect()
# %%
# =============================================================================
# Load ROS dataset and compute anual cycle, annual frequency maps and ROS
# maximum timing (month)
# =============================================================================

ROS_ERA5LAND = xr.open_mfdataset('datos/ROS/ERA5LAND/ROS*').ROS
ROS_ERA5LAND = ROS_ERA5LAND.resample({'time': 'd'}).sum() > 3
lat, lon = ROS_ERA5LAND.lat.values, ROS_ERA5LAND.lon.values
LON, LAT = np.meshgrid(lon, lat)
freq_ERA5LAND = ROS_ERA5LAND.resample({'time': 'y'}).sum()
freq_ERA5LAND = freq_ERA5LAND.mean(dim='time').load()
cycle_ERA5LAND = ROS_ERA5LAND.groupby('time.month').mean().load()
timing = cycle_ERA5LAND.to_dataframe()
timing = timing.groupby(["lat", "lon"]).idxmax().ROS
timing = timing.map(lambda x: x[0] if type(x) == tuple else x).unstack()
gc.collect()
# %%
# =============================================================================
# Build a mask with snow areas (SWE_year_max>20mm)
# =============================================================================
maxSWE = []
for yr in range(1979, 2021, 1):
    path = 'datos/era5land/snow_depth_water_equivalent_'+str(yr)+'.nc'
    SWE = xr.open_dataset(path).sd
    SWE = SWE.resample({'time': 'd'}).sum().max(dim='time')/1e3
    maxSWE.append(SWE)
mask = xr.concat(maxSWE, 'max_years').mean(dim='max_years')
del maxSWE
# %%
# =============================================================================
# Load elevation raster and build same mask for its resolution
# =============================================================================
# DEM = xr.open_dataset('datos/topography/ANDES_static_CORTES.nc')
DEM = xr.open_dataset('datos/topography/Andes_topo_01x01grad.nc')
latC, lonC = DEM.lat.values.squeeze(), DEM.lon.values.squeeze()
LATC, LONC = np.meshgrid(latC, lonC)
DEM = DEM.elevation.values.squeeze().T

maskC = mask.reindex({'lat': latC,
                      'lon': lonC},
                     method='nearest')
mask = mask > 20
maskC = maskC.T > 20
# DEM = DEM.reindex({'lat': lat, 'lon': lon}, method="nearest")

# %%
fig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
                       subplot_kw={"projection": ccrs.PlateCarree()},
                       figsize=(16, 8))
plt.rcParams.update({'font.size': 12})
ax = ax.ravel()
for axis in ax:
    axis.set_extent([-73, -68, -26, -38])
    axis.coastlines()
    axis.add_feature(cf.BORDERS)
    gl = axis.gridlines(linestyle=":")
    gl.xlocator = mticker.FixedLocator([])
    gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = False
    cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
    axis.add_feature(cf.OCEAN)

cmaplist = pd.read_csv("terraincolormap.txt").values
cmaplist = [list(np.array(i)/255) for i in cmaplist]
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap',
                                                    cmaplist,
                                                    len(cmaplist))


mapa0 = ax[0].pcolormesh(LONC, LATC, np.where(maskC.values, DEM, np.nan),
                         cmap=cmap,
                         transform=ccrs.PlateCarree(), vmin=0, vmax=5e3,
                         rasterized=True)

mapa1 = ax[1].pcolormesh(LON, LAT, mean_rainydays.where(mask.values),
                         cmap="cividis_r",
                         transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa2 = ax[2].pcolormesh(LON, LAT, freq_ERA5LAND.where(mask.values),
                         cmap="viridis", transform=ccrs.PlateCarree(),
                         rasterized=True)
mapa3 = ax[3].pcolormesh(LON, LAT, timing.where(mask.values),
                         cmap='nipy_spectral',
                         transform=ccrs.PlateCarree(), rasterized=True)

cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40,
                   label="Orography (m.a.s.l)")
cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40,
                   label="Rain Frequency\n(N°Rainy Days/year)")
cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40,
                   label="ROS Frequency\n(N°ROS/year)")
cb2 = fig.colorbar(mapa3, ax=ax[3], aspect=40, ticks=np.arange(1, 13, 1),
                   label="Maximum ROS Timing")
cb2.set_ticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                    "JUL", "AGO", "SEP", "OCT", "NOV", "DIC"])

gl = ax[0].gridlines(draw_labels=True, linestyle=":")
gl.xlocator = mticker.FixedLocator([])
gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
gl.top_labels = False
gl.right_labels = False
gl.xlines = False
# gl.ylines = False
plt.savefig('plots/ROS_ERA5LAND_73W-68W-26S-38S.pdf', dpi=150,
            bbox_inches="tight")
