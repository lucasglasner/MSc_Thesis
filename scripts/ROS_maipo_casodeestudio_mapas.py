#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 14:47:34 2022

@author: lucas

# =============================================================================
# ROS in Maipo Basin: 2013/08/11 case of study and spatial analysis
# =============================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib as mpl
import datetime
import xarray as xr
import geopandas as gpd
from glob import glob
import cmocean
# %%
# =============================================================================
# LOAD STATIC DATA AND SET TIME FIX
# =============================================================================
date = "2013-08-11"
interval = datetime.datetime.strptime(date, '%Y-%m-%d')
interval = slice(interval-datetime.timedelta(days=10),
                 interval+datetime.timedelta(days=10))
basin = gpd.read_file('datos/vector/RioMaipoEnElManzano.shp')

# %%
# =============================================================================
# LOAD MODIS DATA
# =============================================================================

fSCA_t = xr.open_dataset('datos/modis/MOD10A1_2000-2021.nc', chunks='auto')
fSCA_t = fSCA_t.fSCA.sel(time=interval).load()
fSCA_a = xr.open_dataset('datos/modis/MYD10A1_2000-2021.nc', chunks='auto')
fSCA_a = fSCA_a.fSCA.sel(time=interval).load()

# %%
# =============================================================================
# LOAD ERA5 DATA
# =============================================================================

H0_ERA5 = xr.open_mfdataset(glob('datos/era5/H0_ERA5_*.nc'), chunks='auto')
H0_ERA5 = H0_ERA5.deg0l.reindex(
    {'time': fSCA_t.time}, method='nearest').load()

# %%
# =============================================================================
# LOAD PRECIPITATION DATA
# =============================================================================

PR = xr.open_mfdataset(glob('datos/cr2met/CR2MET_pr_*.nc'), chunks='auto')
PR = PR.pr.reindex({'time': fSCA_t.time}, method='nearest').load()

# %%
# =============================================================================
# LOAD SWE DATA
# =============================================================================

SWE = xr.open_mfdataset(
    glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_[!WY]*.nc'),
    chunks='auto')
maxSWE = xr.where(SWE.SWE.max(dim='time') > 20, 1, 0).load()
SWE = SWE.SWE.reindex({'time': fSCA_t.time}, method='nearest').load()
SWE = SWE.diff(dim='time')

# %%
# =============================================================================
# LOAD ROS DATA
# =============================================================================

ROS = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*.nc'),
                        chunks='auto')
ROS = ROS.ROS.reindex({'time': fSCA_t.time}, method='nearest').load()

# %%

days = ["2013-08-05", "2013-08-07", "2013-08-09",
        "2013-08-11", "2013-08-13"]
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
fig, ax = plt.subplots(3, 5, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(7, 10))
plt.rc('font', size=18)
lon2d, lat2d = np.meshgrid(ROS.lon, ROS.lat)
for axis in ax.ravel():
    axis.set_extent([lon2d.min(), lon2d.max(), lat2d.min(), lat2d.max()],
                    crs=ccrs.PlateCarree())
    axis.coastlines()
    axis.add_feature(cf.BORDERS, ls=":", rasterized=True)
    axis.add_feature(cf.OCEAN, rasterized=True)
    axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
    basin.boundary.plot(ax=axis, transform=ccrs.PlateCarree(),
                        color='k', lw=0.5)

for i in range(len(days)):
    pr_plot = ax[0, i].pcolormesh(lon2d, lat2d,
                                  PR.sel(time=days[i]),
                                  rasterized=True,
                                  cmap='Blues',
                                  norm=mpl.colors.Normalize(0, 60))
    h0_plot = ax[1, i].pcolormesh(lon2d, lat2d,
                                  H0_ERA5.sel(time=days[i])-300,
                                  rasterized=True,
                                  cmap='summer',
                                  norm=mpl.colors.Normalize(0, 3e3))
    SWE_plot = ax[2, i].pcolormesh(lon2d, lat2d,
                                   SWE.sel(time=days[i]),
                                   rasterized=True,
                                   cmap='RdBu',
                                   norm=mpl.colors.TwoSlopeNorm(vmin=-10.,
                                                                vcenter=0.,
                                                                vmax=40))
    ax[1, i].scatter(np.where(ROS.sel(time=days[i]) == 1,
                              lon2d,
                              np.nan)[:-1, :-1]+0.05/2,
                     np.where(ROS.sel(time=days[i]) == 1,
                              lat2d,
                              np.nan)[:-1, :-1]+0.05/2,
                     color='red',
                     s=0.01,
                     rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)
box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()
box3 = ax[2, -1].get_position()

cax1 = fig.add_axes([box1.xmax*1.05, box1.ymin, 0.025, box1.ymax-box1.ymin])
cax2 = fig.add_axes([box2.xmax*1.05, box2.ymin, 0.025, box2.ymax-box2.ymin])
cax3 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])

fig.colorbar(pr_plot, cax=cax1, label='Precipitation\n$(mm/day)$')
fig.colorbar(h0_plot, cax=cax2, label='Freezing Level\n$(m.a.s.l)$')
fig.colorbar(SWE_plot, cax=cax3,
             label='Snow Water Equivalent\nGain/Loss $(mm/day)$')

for axis in [ax[-1, 0], ax[-1, 1], ax[-1, 2], ax[-1, 3], ax[-1, 4]]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([-70.5])
    gl.ylocator = mpl.ticker.FixedLocator([])
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
for axis in [ax[0, 0], ax[1, 0], ax[2, 0]]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([])
    gl.ylocator = mpl.ticker.FixedLocator([-37, -35, -33, -31, -29, -27])
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
plt.savefig('plots/caseofstudy/pr_swe_fl_maps.pdf',
            dpi=150, bbox_inches='tight')
