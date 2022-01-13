#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 14:47:34 2022

@author: lucas

# =============================================================================
# ROS in Maipo Basin: 2013/08/11 case of study and spatial analysis
# =============================================================================
"""


from functions import add_labels
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
import sys
sys.path.append('functions.py')
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

SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv', index_col=0)
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA['IANIGLA'][interval]

# %%
# =============================================================================
# Make plot of maipo manzano fSCA
# =============================================================================
days = ["2013-08-05", "2013-08-07", "2013-08-09",
        "2013-08-11", "2013-08-13"]
lon2d, lat2d = np.meshgrid(fSCA_t.lon, fSCA_t.lat)
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
titles[0] = '2013\n'+titles[0]
titles = [p+'\n'+str(sca)+'%' for p, sca in zip(titles, SCA[days].values)]
fig, ax = plt.subplots(2, 5, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(7, 5))
plt.rc('font', size=18)

for axis in ax.ravel():
    axis.coastlines()
    axis.set_extent([-70.5, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
    axis.add_feature(cf.BORDERS, ls=":")
    axis.add_feature(cf.LAND, alpha=0.2, color='k', rasterized=True)

for i in range(len(days)):
    map1 = ax[0, i].pcolormesh(lon2d, lat2d,
                               fSCA_t.sel(time=days[i]),
                               rasterized=True,
                               cmap=cmocean.cm.ice,
                               norm=mpl.colors.Normalize(0, 100),
                               transform=ccrs.PlateCarree())
    ax[0, i].scatter(np.where(fSCA_t.where(fSCA_t > 100).sel(time=days[i]) > 100,
                              lon2d,
                              np.nan)[:-1, :-1],
                     np.where(fSCA_t.where(fSCA_t > 100).sel(time=days[i]) > 100,
                              lat2d,
                              np.nan)[:-1, :-1],
                     color='red',
                     s=0.0005,
                     rasterized=True)
    map1 = ax[1, i].pcolormesh(lon2d, lat2d,
                               fSCA_a.sel(time=days[i]),
                               rasterized=True,
                               cmap=cmocean.cm.ice,
                               norm=mpl.colors.Normalize(0, 100),
                               transform=ccrs.PlateCarree())
    ax[1, i].scatter(np.where(fSCA_a.where(fSCA_a > 100).sel(time=days[i]) > 100,
                              lon2d,
                              np.nan)[:-1, :-1],
                     np.where(fSCA_a.where(fSCA_a > 100).sel(time=days[i]) > 100,
                              lat2d,
                              np.nan)[:-1, :-1],
                     color='red',
                     s=0.0005,
                     rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)
    basin.boundary.plot(ax=ax[0, i], transform=ccrs.PlateCarree(),
                        colors='forestgreen')
    basin.boundary.plot(ax=ax[1, i], transform=ccrs.PlateCarree(),
                        colors='darkblue')

ax[0, 0].scatter([], [], s=100, marker='s', color='red', label='No Data')
ax[0, 0].scatter([], [], s=100, marker='s',
                 color='forestgreen', label='MODIS/TERRA')
ax[0, 0].scatter([], [], s=100, marker='s',
                 color='darkblue', label='MODIS/AQUA')
ax[0, 0].legend(frameon=False, loc=(-0.14, 1.5), ncol=3,
                fontsize=13)

add_labels(ax, xticks=[-70], yticks=[-33.2, -33.7, -34.2], linewidth=0)
box1, box2 = ax[0, -1].get_position(), ax[-1, -1].get_position()
cax = fig.add_axes([box1.xmax*1.05, box2.ymin, 0.025, box1.ymax-box2.ymin])
fig.text(0.09, 0.905, 'SCA: ', ha='center', va='center', fontsize=14)
fig.colorbar(map1, cax=cax, label='Fraction of Snow\nCover Area (%)')
plt.savefig('plots/caseofstudy_Aug2013/modis_maps.pdf',
            dpi=150, bbox_inches='tight')


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
SWEdiff = xr.open_mfdataset(
    glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_[!WY]*.nc'),
    chunks='auto')
SWEdiff = SWEdiff.SWE.reindex({'time': fSCA_t.time}, method='nearest').load()

# %%
# =============================================================================
# LOAD ROS DATA
# =============================================================================

ROS = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*.nc'),
                        chunks='auto')
ROS = ROS.ROS.reindex({'time': fSCA_t.time}, method='nearest').load()

# %%

days = ["2013-08-05", "2013-08-07", "2013-08-09",
        "2013-08-10", "2013-08-11", "2013-08-13"]
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
titles[0] = '2013\n'+titles[0]
fig, ax = plt.subplots(3, 6, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(10, 10))
plt.rc('font', size=18)
lon2d, lat2d = np.meshgrid(ROS.lon, ROS.lat)
for axis in ax.ravel():
    axis.set_extent([-74, -68, -26, -38],
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
                                   SWEdiff.sel(time=days[i]),
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
                     s=0.1,
                     rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)
box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()
box3 = ax[2, -1].get_position()

cax1 = fig.add_axes([box1.xmax*1.05, box1.ymin, 0.025, box1.ymax-box1.ymin])
cax2 = fig.add_axes([box2.xmax*1.05, box2.ymin, 0.025, box2.ymax-box2.ymin])
cax3 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])

fig.colorbar(pr_plot, cax=cax1, label='Precipitation\n$(mm/day)$')
fig.colorbar(h0_plot, cax=cax2, label='Freezing Level\n$(m.a.g.l)$')
fig.colorbar(SWE_plot, cax=cax3,
             label='Snow Water\nEquivalent Change\n$(mm/day)$')

for axis in ax[-1, :]:
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
plt.savefig('plots/caseofstudy_Aug2013/pr_swe_fl_maps.pdf',
            dpi=150, bbox_inches='tight')


# %%
days = ["2013-08-05", "2013-08-07", "2013-08-09",
        "2013-08-10", "2013-08-11", "2013-08-13"]
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
titles[0] = '2013\n'+titles[0]
fig, ax = plt.subplots(3, 6, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(10, 10))
plt.rc('font', size=18)
lon2d, lat2d = np.meshgrid(ROS.lon, ROS.lat)
for axis in ax.ravel():
    axis.set_extent([-72, -69, -32.4, -37.3])
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
                                   SWEdiff.sel(time=days[i]),
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
                     s=1,
                     rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)
box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()
box3 = ax[2, -1].get_position()

cax1 = fig.add_axes([box1.xmax*1.05, box1.ymin, 0.025, box1.ymax-box1.ymin])
cax2 = fig.add_axes([box2.xmax*1.05, box2.ymin, 0.025, box2.ymax-box2.ymin])
cax3 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])

fig.colorbar(pr_plot, cax=cax1, label='Precipitation\n$(mm/day)$')
fig.colorbar(h0_plot, cax=cax2, label='Freezing Level\n$(m.a.g.l)$')
fig.colorbar(SWE_plot, cax=cax3,
             label='Snow Water\nEquivalent Change\n$(mm/day)$')

for axis in ax[-1, :]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([-70.5])
    gl.ylocator = mpl.ticker.FixedLocator([])
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
for axis in [ax[0, 0], ax[1, 0], ax[2, 0]]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([])
    gl.ylocator = mpl.ticker.FixedLocator([-37, -36, -35, -34, -33])
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
plt.savefig('plots/caseofstudy_Aug2013/pr_swe_fl_maps_maipomanzano.pdf',
            dpi=150, bbox_inches='tight')
