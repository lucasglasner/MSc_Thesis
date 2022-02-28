#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 10:00:56 2022

@author: lucas

# =============================================================================
# This Script makes the synoptic analysis of the 2013/08/13 case of study
# =============================================================================
"""

from functions import add_labels
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs

import sys

sys.path.append("functions.py")

# %%

surface_vars = xr.open_dataset('datos/era5/era5_Ago2013_surface.nc')
upper_vars = xr.open_dataset('datos/era5/era5_Ago2013_upper.nc')

surface_vars = surface_vars.isel(time=slice(4, 14))
upper_vars = upper_vars.isel(time=slice(4, 14))


times = surface_vars.time.to_series().index.strftime('%Y-%b-%d')

lon, lat = surface_vars.lon, surface_vars.lat
lon2d, lat2d = np.meshgrid(lon, lat)

ivt = np.sqrt(surface_vars["p71.162"]**2+surface_vars["p72.162"])
# %%
# =============================================================================
# dynamics
# =============================================================================

plt.rc('font', size=18)
fig, ax = plt.subplots(2, 5, sharex=True, sharey=True,
                       subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(18, 7))


add_labels(ax, [-25, -35, -45], [-90, -80, -70], linewidth=0)
for i, axis in enumerate(ax.ravel()):
    axis.coastlines(rasterized=True)
    axis.set_extent([-100, -60, -20, -50], crs=ccrs.PlateCarree())
    # axis.gridlines(linestyle=":")
    axis.set_title(times[i])


for i, axis in enumerate(ax.ravel()):
    CL = axis.contour(lon2d, lat2d, surface_vars.msl[i, :, :]/100,
                      levels=np.arange(1000, 1033, 3),
                      colors='k', alpha=0.25, rasterized=True)
    axis.clabel(CL, [1006, 1021], fmt='%i', fontsize=11)

    CL = axis.contour(lon2d, lat2d, upper_vars.z[i, 0, :, :]/9.8, colors='k',
                      levels=np.arange(5200, 5850, 80))
    axis.clabel(CL, [5360, 5760], fmt='%i', fontsize=11)

    temp = axis.pcolormesh(lon2d, lat2d, upper_vars.t[i, 1, :, :]-273.15,
                           rasterized=True,
                           cmap='coolwarm')

    # vort = axis.pcolormesh(lon2d,lat2d, upper_vars.vo[i,0,:,:]*1e4,
    #                        rasterized=True, cmap='RdBu_r', vmax=1, vmin=-2)

box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()

cax = fig.add_axes([box2.xmax*1.03, box2.ymin, 0.01, box1.ymax-box2.ymin])
fig.colorbar(temp, cax=cax, label='850hPa Temperature (°C)')

plt.savefig('plots/caseofstudy_Aug2013/synoptic_dynamicanalysis.pdf', dpi=150,
            bbox_inches='tight')


# %%
# =============================================================================
# water
# =============================================================================


plt.rc('font', size=18)
fig, ax = plt.subplots(2, 5, sharex=True, sharey=True,
                       subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(18, 7))


add_labels(ax, [-25, -35, -45], [-90, -80, -70], linewidth=0)
for i, axis in enumerate(ax.ravel()):
    axis.coastlines(rasterized=True)
    axis.set_extent([-100, -60, -20, -50], crs=ccrs.PlateCarree())
    # axis.gridlines(linestyle=":")
    axis.set_title(times[i])


for i, axis in enumerate(ax.ravel()):
    # CL = axis.contour(lon2d, lat2d, surface_vars.tcrw[i, :, :],
    #                   levels=np.arange(1000, 1033, 3),
    #                   colors='k', alpha=0.25, rasterized=True)
    # axis.clabel(CL, [1006, 1021], fmt='%i', fontsize=11)

    # CL = axis.contour(lon2d, lat2d, upper_vars.z[i, 0, :, :]/9.8, colors='k',
    #                   levels=np.arange(5200, 5850, 80))
    # axis.clabel(CL, [5360, 5760], fmt='%i', fontsize=11)

    pw = axis.contourf(lon2d, lat2d, ivt[i, :, :],
                       rasterized=True, cmap='twilight',
                       levels=np.arange(0, 1250+125, 125))

    axis.quiver(lon2d, lat2d, surface_vars["p71.162"][i, :, :].values,
                surface_vars["p72.162"][i, :, :].values, regrid_shape=20,
                transform=ccrs.PlateCarree(), color='k')

    # axis.quiver(lon2d,lat2d,upper_vars.u[i,0,:,:].values,
    #             upper_vars.v[i,0,:,:].values, regrid_shape=20,
    #             transform=ccrs.PlateCarree(), color='k')
    # vort = axis.pcolormesh(lon2d,lat2d, upper_vars.vo[i,0,:,:]*1e4,
    #                        rasterized=True, cmap='RdBu_r', vmax=1, vmin=-2)

box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()

cax = fig.add_axes([box2.xmax*1.03, box2.ymin, 0.01, box1.ymax-box2.ymin])
fig.colorbar(pw, cax=cax,
             label='Integrated Water Vapor\nTransport $(kgm^{-1}s^{-1})$')

plt.savefig('plots/caseofstudy_Aug2013/synoptic_ivtanalysis.pdf', dpi=150,
            bbox_inches='tight')
