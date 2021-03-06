#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 10:12:44 2022

@author: lucas

# =============================================================================
# This Script does the spatial analysis of ROS days in the Andes mountain range
# =============================================================================
"""

from scipy.stats import linregress
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
import cmocean
# %%
# =============================================================================
# Load Basin shapefiles
# =============================================================================
paths = glob('datos/vector/basins/mains/Rio*.shp')
cuencas = []
for path in paths:
    cuencas.append(gpd.read_file(path))
cuencas = pd.concat(cuencas)

# =============================================================================
# Load CR2MET precipitation and compute the mean rainy days per year
# =============================================================================
PR_CR2MET = xr.open_mfdataset('datos/cr2met/CR2MET_pr*', chunks=None).pr
mean_rainydays = xr.where(PR_CR2MET > 10, 1, 0)
mean_rainydays = mean_rainydays.resample(
    {"time": "y"}).sum().mean(dim='time').load()

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
del cycle_CORTESCR2MET
# %%
# =============================================================================
# Build a mask with snow areas (SWE_year_max>20mm)
# =============================================================================
maxSWE = []
for yr in range(1984, 2016, 1):
    path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+str(yr)+'.nc'
    SWE = xr.open_dataset(path).SWE.max(dim='time')
    maxSWE.append(SWE)
mask = xr.where(xr.concat(maxSWE, 'max_years').mean(dim='max_years') > 20 ,1,0)
# mask = (mask) & (mean_rainydays>10)
del maxSWE, SWE
mask1 = mask
mask = xr.where((freq_CORTESCR2MET>1),1,0)
# %%
# =============================================================================
# Load delta SWE data
# =============================================================================
SWE = xr.open_mfdataset(
    glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_*')).SWE


dSWE = (SWE.shift(time=-1)-SWE.shift(time=1))/2


# compute mean of yearly maximums swe loss on ros days
dSWEmax = dSWE.where(ROS.reindex({'time': dSWE.time}) == True)
dSWEmax = dSWE.resample({'time': 'y'}).min().mean(dim='time')
dSWEmax = dSWEmax.compute()
# compute swe loss vs pr in ros days
dswe_vs_pr = -1*dSWE.where((dSWE<0)&(ROS==True)).mean(dim='time')
dswe_vs_pr = dswe_vs_pr/(PR_CR2MET.where(ROS==True).mean(dim='time'))
# dswe_vs_pr = dswe_vs_pr.mean(dim='time').compute()
dswe_vs_pr = dswe_vs_pr.where((dswe_vs_pr<1) & (dswe_vs_pr>0)).compute()
dswe_vs_pr = dswe_vs_pr.reindex({'lat':ROS.lat,'lon':ROS.lon})
# dSWE
# del dSWE
# %%
# =============================================================================
# COmpute ROS trend
# =============================================================================
trend = ROS.resample({'time': 'y'}).sum().load()
trend_new = np.empty((trend.shape[1], trend.shape[2]))
pvalues = np.empty((trend.shape[1], trend.shape[2]))

for i in range(trend.shape[1]):
    for j in range(trend.shape[2]):
        t = linregress(trend.time.to_series().index.year,
                       trend[:, i, j].to_series())
        trend_new[i, j] = t.slope
        pvalues[i, j] = t.pvalue
trend = xr.DataArray(trend_new, coords=[trend.lat, trend.lon],
                     dims=['lat', 'lon'])
pvalues = xr.DataArray(pvalues, coords=[trend.lat, trend.lon],
                       dims=['lat', 'lon'])


# %%
# =============================================================================
# Create Data for maps
# =============================================================================
MROS = xr.where(dSWE<0,ROS,0)
MROS1 = MROS.copy()
MROS = MROS.resample({'time': 'y'}).sum()
MROS = MROS.mean(dim='time').load()
ROS_days = freq_CORTESCR2MET.where(mask)
mean_rainydays = mean_rainydays.reindex({'lat': ROS_days.lat,
                                         'lon': ROS_days.lon},
                                         method='nearest')
timing = timing.where(mask.values.astype(bool))

# ROS_rainratio = PR_CR2MET.where(ROS == True).where(
#     PR_CR2MET > 3).mean(dim='time')
# ROS_rainratio = ROS_rainratio / \
#     PR_CR2MET.where(PR_CR2MET > 3).compute().quantile(0.99, dim='time')
# ROS_rainratio = ROS_rainratio.reindex({'lon': ROS.lon, 'lat': ROS.lat},
#                                       method='nearest')
# ROS_rainratio = ROS_rainratio.compute()
#%%

swe_pr = -dSWE.where(MROS1 == True).mean(dim='time')
swe_pr = swe_pr/PR_CR2MET.where(MROS1 == True).mean(dim='time')
# swe_pr = swe_pr.groupby('time.year').max()
swe_pr = swe_pr.load()*100
swe_pr = swe_pr.where(swe_pr<100)

dem = xr.open_dataset('datos/topography/Andes_topo_005x005grad.nc')
# # %%
# fig, ax = plt.subplots(1, 6, sharex=True, sharey=True,
#                         subplot_kw={"projection": ccrs.PlateCarree()},
#                         figsize=(18, 6))


# plt.rcParams.update({'font.size': 18})
# ax = ax.ravel()
# for axis in ax:
#     axis.set_extent([-74, -68, -26, -38])
#     # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
#     axis.coastlines(rasterized=True)
#     axis.add_feature(cf.BORDERS)
#     gl = axis.gridlines(linestyle=":")
#     gl.xlocator = mticker.FixedLocator([])
#     gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
#     gl.top_labels = False
#     gl.right_labels = False
#     gl.xlines = False
#     cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
#     axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
#     axis.add_feature(cf.OCEAN, rasterized=True)
#     # axis.set_extent([-72, -69, -32.4, -37.3])
#     # # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
#     # axis.coastlines(rasterized=True)
#     # axis.add_feature(cf.BORDERS)
#     # gl = axis.gridlines(linestyle=":")
#     # gl.xlocator = mticker.FixedLocator([])
#     # gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
#     # gl.top_labels = False
#     # gl.right_labels = False
#     # gl.xlines = False
#     # cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
#     # axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
#     # axis.add_feature(cf.OCEAN, rasterized=True)


# mapa0 = ax[0].pcolormesh(LON, LAT, dSWEmax.where(mask),
#                           cmap=cmocean.cm.deep,
#                           transform=ccrs.PlateCarree(),
#                           rasterized=True,
#                           norm=mpl.colors.Normalize(-25, 0))
# mapa1 = ax[1].pcolormesh(LON, LAT, dswe_vs_pr.where(mask),
#                           rasterized=True,
#                           transform=ccrs.PlateCarree(),
#                           cmap=cmocean.cm.ice_r,
#                           norm=mpl.colors.Normalize(0, 0.5))
# mapa2 = ax[2].pcolormesh(LON, LAT, ROS_rainratio.where(mask),
#                           rasterized=True,
#                           transform=ccrs.PlateCarree(),
#                           cmap=cmocean.cm.rain,
#                           norm=mpl.colors.Normalize(0, 1))

# mapa3 = ax[3].pcolormesh(LON, LAT,
#                           trend.where(mask).where(trend != 0),
#                           cmap="RdBu_r",
#                           transform=ccrs.PlateCarree(),
#                           rasterized=True,
#                           norm=mpl.colors.Normalize(-12, 12))
# # mapa2 = ax[2].pcolormesh(LON, LAT,
# #                          pvalues.where(mask).where(trend != 0),
# #                          cmap='Purples',
# #                          transform=ccrs.PlateCarree(),
# #                          rasterized=True,
# #                          norm=mpl.colors.Normalize(0, 1))
# # mapa3 = ax[3].pcolormesh(LON, LAT,
# #                          ROS_rainratio.where(mask),
# #                          cmap='PRGn',
# #                          transform=ccrs.PlateCarree(),
# #                          rasterized=True,
# #                          norm=mpl.colors.TwoSlopeNorm(1, 0, 2))

# # ax[1].pcolormesh(LON, LAT,
# #                  pvalues.where(mask).where(pvalues > 0.05).where(trend != 0),
# #                  cmap="bone", alpha=0.1, vmin=1, vmax=1,
# #                  rasterized=True, linewidth=0)


# ax[0].set_title("Yearly mean of\nmaximum SWE loss\non ROS days\n"+r"($mm$)")
# ax[1].set_title("SWE contribution\nto total aviable\nwater for runoff\n(-)")
# ax[2].set_title("ROS mean\nrainfall over\nextreme rainfall\n(-)")
# ax[3].set_title("Yearly ROS\ndays Trend\n"+r"($\frac{N??ROS Days}{year}$)")
# # ax[2].set_title("Yearly ROS\ndays Trend\npvalue"+"\n(-)")
# # ax[3].set_title("")

# cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40)
# cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40)
# cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40)
# cb3 = fig.colorbar(mapa3, ax=ax[3], aspect=40)

# gl = ax[0].gridlines(draw_labels=True, linestyle=":")
# gl.xlocator = mticker.FixedLocator([])
# gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
# gl.top_labels = False
# gl.right_labels = False
# gl.xlines = False
# gl.ylines = False

# # gl = ax[0].gridlines(draw_labels=True, linestyle=":")
# # gl.xlocator = mticker.FixedLocator([])
# # gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
# # gl.top_labels = False
# # gl.right_labels = False
# # gl.xlines = False
# # gl.ylines = False
# plt.savefig('plots/ROS_STATS_final.pdf', dpi=150, bbox_inches='tight')

# #%%

# fig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
#                        subplot_kw={"projection": ccrs.PlateCarree()},
#                        figsize=(18, 6))


# plt.rcParams.update({'font.size': 18})
# ax = ax.ravel()
# for axis in ax:
#     axis.set_extent([-72, -69, -32.4, -37.3])
#     # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
#     axis.coastlines(rasterized=True)
#     axis.add_feature(cf.BORDERS)
#     gl = axis.gridlines(linestyle=":")
#     gl.xlocator = mticker.FixedLocator([])
#     gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
#     gl.top_labels = False
#     gl.right_labels = False
#     gl.xlines = False
#     cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
#     axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
#     axis.add_feature(cf.OCEAN, rasterized=True)


# mapa0 = ax[0].pcolormesh(LON, LAT, dSWEmax.where(mask),
#                          cmap=cmocean.cm.deep,
#                          transform=ccrs.PlateCarree(),
#                          rasterized=True,
#                          norm=mpl.colors.Normalize(-25, 0))
# mapa1 = ax[1].pcolormesh(LON, LAT, dswe_vs_pr.interpolate_na(dim='lon').where(mask),
#                          rasterized=True,
#                          transform=ccrs.PlateCarree(),
#                          cmap=cmocean.cm.ice_r,
#                          norm=mpl.colors.Normalize(0, 0.5))
# mapa2 = ax[2].pcolormesh(LON, LAT, ROS_rainratio.where(mask),
#                          rasterized=True,
#                          transform=ccrs.PlateCarree(),
#                          cmap=cmocean.cm.haline,
#                          norm=mpl.colors.Normalize(0, .2))

# mapa3 = ax[3].pcolormesh(LON, LAT,
#                          trend.where(mask).where(trend != 0),
#                          cmap="RdBu_r",
#                          transform=ccrs.PlateCarree(),
#                          rasterized=True,
#                          norm=mpl.colors.Normalize(-12, 12))
# # mapa2 = ax[2].pcolormesh(LON, LAT,
# #                          pvalues.where(mask).where(trend != 0),
# #                          cmap='Purples',
# #                          transform=ccrs.PlateCarree(),
# #                          rasterized=True,
# #                          norm=mpl.colors.Normalize(0, 1))
# # mapa3 = ax[3].pcolormesh(LON, LAT,
# #                          ROS_rainratio.where(mask),
# #                          cmap='PRGn',
# #                          transform=ccrs.PlateCarree(),
# #                          rasterized=True,
# #                          norm=mpl.colors.TwoSlopeNorm(1, 0, 2))

# # ax[1].pcolormesh(LON, LAT,
# #                  pvalues.where(mask).where(pvalues > 0.05).where(trend != 0),
# #                  cmap="bone", alpha=0.1, vmin=1, vmax=1,
# #                  rasterized=True, linewidth=0)


# ax[0].set_title("Yearly mean of\nmaximum SWE loss\non ROS days\n"+r"($mm$)")
# ax[1].set_title("SWE contribution\nover total\nprecipitation\n(-)")
# ax[2].set_title("ROS mean\nrainfall over\nextreme rainfall\n(-)")
# ax[3].set_title("Yearly ROS\ndays Trend\n"+r"($\frac{N??ROS Days}{year}$)")
# # ax[2].set_title("Yearly ROS\ndays Trend\npvalue"+"\n(-)")
# # ax[3].set_title("")

# cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40)
# cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40)
# cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40)
# cb3 = fig.colorbar(mapa3, ax=ax[3], aspect=40)


# gl = ax[0].gridlines(draw_labels=True, linestyle=":")
# gl.xlocator = mticker.FixedLocator([])
# gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
# gl.top_labels = False
# gl.right_labels = False
# gl.xlines = False
# gl.ylines = False
# plt.savefig('plots/ROS_STATS_final_small.pdf', dpi=150, bbox_inches='tight')
# # # %%
# fig, ax = plt.subplots(1, 4, sharex=True, sharey=True,
#                        subplot_kw={"projection": ccrs.PlateCarree()},
#                        figsize=(18, 6))
# ax2 = []
# plt.rcParams.update({'font.size': 18})
# ax = ax.ravel()
# for axis in ax:
#     axis.set_extent([-74, -68, -26, -38])
#     # axis.set_extent([-74, -69.7, -33, -34.32], crs=ccrs.PlateCarree())
#     axis.coastlines(rasterized=True)
#     axis.add_feature(cf.BORDERS)
#     gl = axis.gridlines(linestyle=":")
#     gl.xlocator = mticker.FixedLocator([])
#     gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
#     gl.top_labels = False
#     gl.right_labels = False
#     gl.xlines = False
#     cuencas.boundary.plot(ax=axis, color="k", lw=0.5)
#     axis.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)
#     axis.add_feature(cf.OCEAN, rasterized=True)


# # mask = (ROS_days/mean_rainydays != 0).values
# mapa0 = ax[0].pcolormesh(LON, LAT, mean_rainydays,
#                          cmap="cividis_r",
#                          transform=ccrs.PlateCarree(),
#                          rasterized=True)
# mapa1 = ax[1].pcolormesh(LON, LAT, ROS_days,
#                          cmap="viridis",
#                          transform=ccrs.PlateCarree(), vmin=0, vmax=30,
#                          rasterized=True)
# mapa2 = ax[2].pcolormesh(LON, LAT,
#                          (ROS_days/mean_rainydays),
#                          cmap="BuPu", transform=ccrs.PlateCarree(),
#                          rasterized=True)
# mapa3 = ax[3].pcolormesh(LON, LAT, timing.where(mask.values.astype(bool)),
#                          cmap='nipy_spectral',
#                          transform=ccrs.PlateCarree(), rasterized=True)


# ax[0].set_title("Rain Frequency\n"+r"($\frac{N??Rainy Days}{year}$)")
# ax[1].set_title("ROS Frequency\n"+r"($\frac{N??ROS Days}{year}$)")
# ax[2].set_title("ROS/Rain\nFrequency\nRatio (%)")
# ax[3].set_title("Maximum ROS\nTiming")

# cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=40)
# cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=40)
# cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=40)
# cb2 = fig.colorbar(mapa3, ax=ax[3], aspect=40, ticks=np.arange(1, 13, 1))
# # cb2 = fig.colorbar(mapa3, ax=ax[3], aspect=40, ticks=np.linspace(1, 365, 12))
# cb2.set_ticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
#                     "JUL", "AUG", "SEP", "OCT", "NOV", "DIC"])


# gl = ax[0].gridlines(draw_labels=True, linestyle=":")
# gl.xlocator = mticker.FixedLocator([])
# gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])
# gl.top_labels = False
# gl.right_labels = False
# gl.xlines = False
# gl.ylines = False
# plt.savefig('plots/ROS_CORTESCR2METERA5_74W-68W-26S-38S_final.pdf', dpi=150,
#             bbox_inches="tight")

# %%

dz = 500
z = np.arange(0,6000+dz,dz)

x = [xr.where((dem>z[i])&(dem<z[i+1]),ROS_days,np.nan).mean().item()
     for i in range(len(z)-1)]

#%%
fig=plt.figure(figsize=(8,3))
ax = fig.add_subplot(111)
plt.rc('font',size=18)
ax.bar(np.arange(len(x)),x, ec='k')
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(["[0,500)", "[500,1000)","[1000,1500)",
                    "[1500,2000)","[2000,2500)","[2500,3000)",
                    "[3000,3500)","[3500,4000)","[4000,4500)",
                    "[4500,5000)","[5000,)",""])
ax.tick_params(rotation=90, axis="x", labelsize=14)

ax.set_ylabel('Mean Ros Days per Year')
ax.set_xlabel('Elevation band (m)')

plt.savefig('plots/elevation_band.pdf',dpi=150,bbox_inches='tight')
#%%

fig, ax = plt.subplots(1, 6, sharex=True, sharey=True,
                       subplot_kw={"projection": ccrs.PlateCarree()},
                       figsize=(20, 8))
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


# mask = (ROS_days/mean_rainydays != 0).values

mapa0 = ax[0].pcolormesh(LON, LAT, ROS_days.where(mask),
                         cmap="viridis",
                         transform=ccrs.PlateCarree(), vmin=0, vmax=30,
                         rasterized=True)
mapa1 = ax[1].pcolormesh(LON, LAT,
                          (ROS_days/mean_rainydays).where(mask)*100,
                          cmap="BuPu", transform=ccrs.PlateCarree(),
                          rasterized=True)
mapa2 = ax[2].pcolormesh(LON, LAT, 
                         (MROS/mean_rainydays).where(mask)*100,
                         cmap=cmocean.cm.tempo,vmin=0,
                         transform=ccrs.PlateCarree(), rasterized=True)

mapa3 = ax[3].pcolormesh(LON, LAT, timing.where(mask.values.astype(bool)),
                          cmap='nipy_spectral',
                          transform=ccrs.PlateCarree(), rasterized=True)
mapa4 = ax[4].pcolormesh(swe_pr.lon, swe_pr.lat, swe_pr.where(mask),
                         cmap='ocean',vmax=50,vmin=0,
                         transform=ccrs.PlateCarree(), rasterized=True)
mapa5 = ax[5].pcolormesh(LON,LAT, trend.where(mask),
                         cmap='coolwarm',norm=mpl.colors.TwoSlopeNorm(0),
                         transform=ccrs.PlateCarree(), rasterized=True)

ax[5].pcolor(LON,LAT,trend.where(mask).where(pvalues<0.05),
             hatch=".....", alpha=0, rasterized=True)

ax[0].set_title("(a)\nROS Frequency\n"+r"($\frac{N??ROS Days}{year}$)")
ax[1].set_title("(b)\nROS/Rain\nFrequency\nRatio (%)")
ax[2].set_title("(c)\nMROS/Rain\nFrequency\nRatio (%)")
ax[3].set_title("(d)\nMaximum ROS\nTiming (month)")

ax[4].set_title("(e)\nMelt snow\non ROS days\n(% of liquid\nprecipitation)")
ax[5].set_title("(f)\nROS days trend")

cb0 = fig.colorbar(mapa0, ax=ax[0], aspect=15, shrink=1, location='bottom')
cb0.set_ticks([0,10,20,30])
cb1 = fig.colorbar(mapa1, ax=ax[1], aspect=15, shrink=1, location='bottom')
cb1.set_ticks([0,15,30,45,60])
cb2 = fig.colorbar(mapa2, ax=ax[2], aspect=15, shrink=1, location='bottom')
cb2.set_ticks([0,5,10,15,20,25])
cb3 = fig.colorbar(mapa3, ax=ax[3], aspect=15, shrink=1, location='bottom',
                    ticks=np.arange(1, 13, 1))
cb3.set_ticklabels(["J", "F", "M", "A", "M", "J",
                    "J", "A", "S", "O", "N", "D"])
cb4 = fig.colorbar(mapa4, ax=ax[4], aspect=15, shrink=1, location='bottom',
                   ticks=[0,10,20,30,40,50])
cb5 = fig.colorbar(mapa5, ax=ax[5], aspect=15, shrink=1, location='bottom')
cb5.set_ticks([-0.6,-0.3,0,0.3])
# cb5.set_ticklabels([])

gl = ax[0].gridlines(draw_labels=True, linestyle=":")
gl.xlocator = mticker.FixedLocator([])
gl.ylocator = mticker.FixedLocator([-37, -36, -35, -34, -33])
gl.top_labels = False
gl.right_labels = False
gl.xlines = False
gl.ylines = False
plt.savefig('plots/ROS_CORTESCR2METERA5_FINAL.pdf', dpi=150,
            bbox_inches="tight")
