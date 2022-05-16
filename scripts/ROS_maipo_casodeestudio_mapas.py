#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 14:47:34 2022

@author: lucas

# =============================================================================
# ROS in Maipo Basin: case of study and spatial analysis
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
# %%
date = "2008-06-04"
# date = "%YR%"
yr, month, day = [int(n) for n in date.split("-")]
interval = slice(datetime.datetime(yr, month, day)-datetime.timedelta(days=12),
                 datetime.datetime(yr, month, day)+datetime.timedelta(days=12))

basins = pd.concat([gpd.read_file(p)
                   for p in glob('datos/vector/basins/mains/*.shp')])
basins.index = basins.gauge_id
maipobasin = gpd.read_file('datos/vector/basins/RioMaipoEnElManzano.shp')

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
# SWEdiff = xr.open_mfdataset(
#     glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_[!WY]*.nc'),
#     chunks='auto')
# SWEdiff = SWEdiff.SWE.reindex({'time': fSCA_t.time}, method='nearest').load()
SWEdiff = (SWE.shift(time=-1)-SWE.shift(time=1))/2

SWE_mm = xr.open_dataset(
    'datos/ANDES_SWE_Cortes/maipomanzano/ANDES_SWE_WY'+str(interval.start.year+1)+'.nc')
SWE_mm = SWE_mm.SWE.sel(time=interval)

cortes_SCA = xr.where(SWE_mm > 100, 1, 0).sum(dim=['lat', 'lon'])
cortes_SCA = cortes_SCA.to_series()/467798

#%%

# fig,ax = plt.subplots(1,1)
# (SCA/100).plot(ax=ax)
# cortes_SCA.plot(ax=ax)
# 
# %%
# =============================================================================
# LOAD ROS DATA
# =============================================================================

# ROS = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*.nc'),
#                         chunks='auto')
# ROS = ROS.ROS.reindex({'time': fSCA_t.time}, method='nearest').load()
dem = xr.open_dataset('datos/topography/Andes_topo_005x005grad.nc').elevation




dem = dem.reindex({'lon': SWE.lon, 'lat': SWE.lat}, method='nearest')

from metpy.units import units
from metpy.interpolate import cross_section
start = (-50, -76.3)
end = (-20, -71.3)

H0_ERA5_coast = H0_ERA5.copy().to_dataset(name='h0').metpy.parse_cf()

H0_ERA5_coast = cross_section(H0_ERA5_coast, start, end).set_coords(('lat', 'lon')).h0
H0_ERA5_coast = H0_ERA5_coast.swap_dims({'index':'lat'})
H0_ERA5_coast = H0_ERA5_coast.reindex({'lat':SWE.lat}, method='nearest').drop('lon')
H0_ERA5_coast = [H0_ERA5_coast for l in SWE.lon]
H0_ERA5_coast = xr.concat(H0_ERA5_coast,dim=SWE.lon)

# H0_ERA5_coast = H0_ERA5_coast.to_series().unstack().T[::-1]
# H0_ERA5_coast.index = H0_ERA5_coast.index.get_level_values(0)

dem = [dem for t in SWE.time]
dem = xr.concat(dem,SWE.time)

# freeze = []
# for t in H0_ERA5_coast.index:
#     mask = xr.where(dem < H0_ERA5_coast.loc[t]-300, True, False)
#     freeze.append(mask)

freeze = xr.where(dem<H0_ERA5_coast-300,True,False)


ROS = xr.where((SWE>10) & (PR>10) & (H0_ERA5>300) & (SWEdiff/SWE<=0.05),
               True, False)
ROS = ROS.reindex(lat=SWE.lat,lon=SWE.lon, method='nearest')

# # %%
# # =============================================================================
# # Make plot of maipo manzano fSCA
# # =============================================================================
# days = ["2008-05-23","2008-05-29","2008-06-02","2008-06-04","2008-06-06"]
# titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
#           for d in days]
# titles[0] = '2008\n'+titles[0]
# titles = [p+'\n'+str(sca)+'%' for p, sca in zip(titles, SCA[days].values)]
# fig, ax = plt.subplots(2, 5, subplot_kw={'projection': ccrs.PlateCarree()},
#                         figsize=(8, 7))
# plt.rc('font', size=18)

# for axis in ax.ravel():
#     axis.set_extent([-70.5, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
#     axis.add_feature(cf.BORDERS, ls=":")
#     axis.add_feature(cf.LAND, alpha=0.2, color='k', rasterized=True)

# for i in range(len(days)):
#     lon2d, lat2d = np.meshgrid(fSCA_t.lon, fSCA_t.lat)
#     map1 = ax[0, i].pcolormesh(lon2d, lat2d,
#                                 fSCA_t.sel(time=days[i]),
#                                 rasterized=True,
#                                 cmap=cmocean.cm.ice,
#                                 norm=mpl.colors.Normalize(0, 100),
#                                 transform=ccrs.PlateCarree())
#     ax[0, i].scatter(np.where(fSCA_t.where(fSCA_t > 100).sel(time=days[i]) > 100,
#                               lon2d,
#                               np.nan)[:-1, :-1],
#                       np.where(fSCA_t.where(fSCA_t > 100).sel(time=days[i]) > 100,
#                               lat2d,
#                               np.nan)[:-1, :-1],
#                       color='red',
#                       s=0.0005,
#                       rasterized=True)
#     map1 = ax[1, i].pcolormesh(lon2d, lat2d,
#                                 fSCA_a.sel(time=days[i]),
#                                 rasterized=True,
#                                 cmap=cmocean.cm.ice,
#                                 norm=mpl.colors.Normalize(0, 100),
#                                 transform=ccrs.PlateCarree())
#     ax[1, i].scatter(np.where(fSCA_a.where(fSCA_a > 100).sel(time=days[i]) > 100,
#                               lon2d,
#                               np.nan)[:-1, :-1],
#                       np.where(fSCA_a.where(fSCA_a > 100).sel(time=days[i]) > 100,
#                               lat2d,
#                               np.nan)[:-1, :-1],
#                       color='red',
#                       s=0.0005,
#                       rasterized=True)

#     ax[0, i].set_title(titles[i], fontsize=14)
#     maipobasin.boundary.plot(ax=ax[0, i], transform=ccrs.PlateCarree(),
#                               colors='forestgreen')
#     maipobasin.boundary.plot(ax=ax[1, i], transform=ccrs.PlateCarree(),
#                               colors='darkblue')


# ax[0, 0].scatter([], [], s=100, marker='s', color='red', label='No Data')
# ax[0, 0].scatter([], [], s=100, marker='s',
#                   color='forestgreen', label='MODIS/TERRA')
# ax[0, 0].scatter([], [], s=100, marker='s',
#                   color='darkblue', label='MODIS/AQUA')
# ax[0, 0].legend(frameon=False, loc=(-1.1, 1.5), ncol=4,
#                 fontsize=13)

# add_labels(ax, xticks=[-70], yticks=[-33.2, -33.7, -34.2], linewidth=0)
# box1, box2 = ax[0, -1].get_position(), ax[1, -1].get_position()
# cax = fig.add_axes([box1.xmax*1.05, box2.ymin, 0.025, box1.ymax-box2.ymin])
# fig.text(0.09, 0.905, 'SCA: ', ha='center', va='center', fontsize=14)
# fig.colorbar(map1, cax=cax, label='Fraction of Snow\nCover Area (%)')

# # box3 = ax[-1, -1].get_position()
# # cax2 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])

# # plt.savefig('plots/caseofstudy_Aug2013/modis_maps.pdf',
#             # dpi=150, bbox_inches='tight')

# %%


days = ["2008-05-25","2008-05-26", "2008-05-27", "2008-05-31",
        "2008-06-03", "2008-06-04", "2008-06-05"]

# days = ["2013-08-05", "2013-08-06", "2013-08-07","2013-08-08",
#         "2013-08-10","2013-08-11","2013-08-12"]
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
titles[0] = date[:4]+'\n'+titles[0]
fig, ax = plt.subplots(3, 7, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(13, 10))
plt.rc('font', size=18)
lon2d, lat2d = np.meshgrid(ROS.lon, ROS.lat)
for axis in ax.ravel():
    axis.add_feature(cf.OCEAN, rasterized=True)
    axis.set_extent([-74, -68, -26, -38],
                    crs=ccrs.PlateCarree())
    axis.set_xticks([])
    axis.set_yticks([])
    axis.coastlines()
    axis.add_feature(cf.BORDERS, ls=":", rasterized=True)
    
    axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
    # basins.boundary.plot(ax=axis, color='k', lw=0.5)
    # maipobasin.boundary.plot(ax=axis,
                              # color='tab:red', lw=0.5)
for i in range(len(days)):
    pr_plot = ax[0, i].pcolormesh(lon2d, lat2d,
                                  PR.sel(time=days[i]),
                                  rasterized=True,
                                  cmap='Blues',
                                  norm=mpl.colors.Normalize(0, 80),
                                  transform=ccrs.PlateCarree())
    h0_plot = ax[1, i].contourf(lon2d, lat2d,
                                np.clip(H0_ERA5.sel(
                                    time=days[i])-300, 0, 4200),
                                rasterized=True,
                                cmap='summer',
                                norm=mpl.colors.Normalize(0, 3e3),
                                transform=ccrs.PlateCarree(),
                                levels=np.arange(0, 5500, 300))
    # SWE_plot = ax[2, i].pcolormesh(lon2d, lat2d,
    #                                SWE.sel(time=days[i]),
    #                                rasterized=True,
    #                                cmap=cmocean.cm.ice,
    #                                norm=mpl.colors.LogNorm(1, 1e3),
    #                                transform=ccrs.PlateCarree())
    dSWE_plot = ax[2, i].pcolormesh(lon2d, lat2d,
                                    SWEdiff.sel(time=days[i]),
                                    rasterized=True,
                                    cmap='RdBu',
                                    norm=mpl.colors.TwoSlopeNorm(vmin=-10.,
                                                                  vcenter=0.,
                                                                  vmax=40),
                                    transform=ccrs.PlateCarree())
    # ax[2, i].contour(lon2d, lat2d, dem,
    #                   colors='tab:red',
    #                   levels=[1500],
    # #                   linewidths=0.5)
    ax[1, i].scatter(np.where(ROS.sel(time=days[i]) == 1,
                              lon2d,
                              np.nan)[:-1, :-1]+0.05/2,
                      np.where(ROS.sel(time=days[i]) == 1,
                              lat2d,
                              np.nan)[:-1, :-1]+0.05/2,
                      color='purple',
                      s=0.1,
                      rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)
box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()
box3 = ax[2, -1].get_position()
# box4 = ax[3, -1].get_position()

cax1 = fig.add_axes([box1.xmax*1.05, box1.ymin, 0.025, box1.ymax-box1.ymin])
cax2 = fig.add_axes([box2.xmax*1.05, box2.ymin, 0.025, box2.ymax-box2.ymin])
cax3 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])
# cax4 = fig.add_axes([box4.xmax*1.05, box4.ymin, 0.025, box4.ymax-box4.ymin])


fig.colorbar(pr_plot, cax=cax1, label='Precipitation\n$(mm/day)$')
norm = mpl.colors.Normalize(0, 3500)
im = mpl.cm.ScalarMappable(norm=norm, cmap='summer')

fig.colorbar(im, cax=cax2, label='Freezing Level\n$(m.a.g.l)$',
             ticks=np.arange(0, 4200, 600), boundaries=np.linspace(0, 3600))
# fig.colorbar(SWE_plot, cax=cax3, ticks=[1, 1e1, 1e2, 1e3],
#               label='Snow Water\n Equivalent\n (mm)')
fig.colorbar(dSWE_plot, cax=cax3,
             label='Snow Water\nEquivalent Change\n$(mm/day)$')

for axis in ax[-1, :]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([-70.5])
    gl.ylocator = mpl.ticker.FixedLocator([])
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
for axis in ax[:, 0]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([])
    gl.ylocator = mpl.ticker.FixedLocator([-37, -35, -33, -31, -29, -27])
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
plt.savefig('plots/caseofstudy_Jun2008/pr_swe_fl_maps.pdf',
            dpi=150, bbox_inches='tight')


# %%
days = ["2008-05-25","2008-05-26", "2008-05-27", "2008-05-30",
        "2008-06-03", "2008-06-04", "2008-06-05"]

# days = ["2013-08-06", "2013-08-07", "2013-08-08","2013-08-09",
        # "2013-08-10","2013-08-11","2013-08-12"]
titles = [datetime.datetime.strptime(d, '%Y-%m-%d').strftime('%b-%d')
          for d in days]
titles[0] = date[:4]+'\n'+titles[0]
fig, ax = plt.subplots(3, 7, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(13, 10))
# fig.tight_layout(pad=1)
plt.rc('font', size=18)
lon2d, lat2d = np.meshgrid(ROS.lon, ROS.lat)
for axis in ax.ravel():
    axis.set_extent([-72.1, -69.5, -32.4, -37])
    axis.coastlines()
    axis.add_feature(cf.BORDERS, ls=":", rasterized=True)
    axis.add_feature(cf.OCEAN, rasterized=True)
    axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
    basins.boundary.plot(ax=axis, color='k', lw=0.5)
    maipobasin.boundary.plot(ax=axis,
                             color='tab:red', lw=0.5)

for i in range(len(days)):
    pr_plot = ax[0, i].pcolormesh(lon2d, lat2d,
                                  PR.sel(time=days[i]),
                                  rasterized=True,
                                  cmap='Blues',
                                  norm=mpl.colors.Normalize(0, 80),
                                  transform=ccrs.PlateCarree())
    h0_plot = ax[1, i].contourf(lon2d, lat2d,
                                np.clip(H0_ERA5.sel(
                                    time=days[i])-300, 0, 4200),
                                rasterized=True,
                                cmap='summer',
                                norm=mpl.colors.Normalize(0, 3e3),
                                transform=ccrs.PlateCarree(),
                                levels=np.arange(0, 5500, 300))
    dSWE_plot = ax[2, i].pcolormesh(lon2d, lat2d,
                                    SWEdiff.sel(time=days[i]),
                                    rasterized=True,
                                    cmap='RdBu',
                                    norm=mpl.colors.TwoSlopeNorm(vmin=-10.,
                                                                 vcenter=0.,
                                                                 vmax=40),
                                    transform=ccrs.PlateCarree())
    # SWE_plot = ax[3, i].pcolormesh(lon2d, lat2d,
    #                                 SWE.sel(time=days[i]),
    #                                 rasterized=True,
    #                                 cmap=cmocean.cm.ice,
    #                                 norm=mpl.colors.LogNorm(1, 1e3),
    #                                 transform=ccrs.PlateCarree())
    # c = ax[3, i].contour(lon2d, lat2d, dem,
    #                       colors='tab:red',
    #                       levels=[1500],
    #                       linewidths=0.8)
    ax[1, i].scatter(np.where(ROS.sel(time=days[i]) == 1,
                              lon2d,
                              np.nan)[:-1, :-1]+0.05/2,
                      np.where(ROS.sel(time=days[i]) == 1,
                              lat2d,
                              np.nan)[:-1, :-1]+0.05/2,
                      color='purple',
                      s=0.5,
                      rasterized=True)
    ax[0, i].set_title(titles[i], fontsize=14)


box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()
box3 = ax[2, -1].get_position()
# box4 = ax[3, -1].get_position()

norm = mpl.colors.Normalize(0, 3500)
im = mpl.cm.ScalarMappable(norm=norm, cmap='summer')


cax1 = fig.add_axes([box1.xmax*1.05, box1.ymin, 0.025, box1.ymax-box1.ymin])
cax2 = fig.add_axes([box2.xmax*1.05, box2.ymin, 0.025, box2.ymax-box2.ymin])
cax3 = fig.add_axes([box3.xmax*1.05, box3.ymin, 0.025, box3.ymax-box3.ymin])
# cax4 = fig.add_axes([box4.xmax*1.05, box4.ymin, 0.025, box4.ymax-box4.ymin])

fig.colorbar(pr_plot, cax=cax1, label='Precipitation\n$(mm/day)$')

fig.colorbar(im, cax=cax2, label='Freezing Level\n$(m.a.g.l)$',
             ticks=np.arange(0, 4200, 600), boundaries=np.linspace(0, 3600))
fig.colorbar(dSWE_plot, cax=cax3,
             label='Snow Water\nEquivalent Change\n$(mm/day)$')
# fig.colorbar(SWE_plot, cax=cax4, ticks=[1, 1e1, 1e2, 1e3],
#              label='Snow Water\n Equivalent\n (mm)')

for axis in ax[-1, :]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([-70.5])
    gl.ylocator = mpl.ticker.FixedLocator([])
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
for axis in ax[:, 0]:
    gl = axis.gridlines(linewidth=0, draw_labels=True)
    gl.xlocator = mpl.ticker.FixedLocator([])
    gl.ylocator = mpl.ticker.FixedLocator([-37, -36, -35, -34, -33])
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
# plt.savefig('plots/caseofstudy_Jun2008/pr_swe_fl_maps_maipomanzano.pdf',
            # dpi=150, bbox_inches='tight')
plt.savefig('plots/caseofstudy_Jun2008/pr_swe_fl_maps_maipomanzano.pdf',
            dpi=150, bbox_inches='tight')
# 