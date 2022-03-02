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
import geopandas as gpd
import matplotlib as mpl

import sys

sys.path.append("functions.py")

# %%

surface_vars = xr.open_dataset(
    'datos/era5/caseofstudy_Ago2013/era5_Ago2013_surface.nc')
upper_vars = xr.open_dataset(
    'datos/era5/caseofstudy_Ago2013/era5_Ago2013_upper.nc')


days = pd.date_range('2013-08-04T00:00:00', '2013-08-14T00:00:00', freq='d')


times = days.strftime('%b-%d\n%H:%MZ')
times = list(times)
times[0] = '2013\n'+times[0]

lon, lat = surface_vars.lon, surface_vars.lat
lon2d, lat2d = np.meshgrid(lon, lat)

ivt = np.sqrt(surface_vars["p71.162"]**2+surface_vars["p72.162"])
# %%
# =============================================================================
# dynamics, temperatura 900hPa, Z500, SLP
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
    axis.set_title(times[i], loc='left')


for i, axis in enumerate(ax.ravel()):
    CL = axis.contour(lon2d, lat2d, surface_vars.msl.sel(time=days)[i, :, :]/100,
                      levels=np.arange(1000, 1033, 3),
                      colors='k', alpha=0.25, rasterized=True)
    axis.clabel(CL, [1006, 1021], fmt='%i', fontsize=11)

    CL = axis.contour(lon2d, lat2d, upper_vars.sel(level=500,
                                                   time=days).z[i, :, :]/9.8,
                      colors='k',
                      levels=np.arange(5200, 5850, 80))
    axis.clabel(CL, [5360, 5760], fmt='%i', fontsize=11)

    temp = axis.pcolormesh(lon2d, lat2d,
                           upper_vars.sel(
                               level=900, time=days).t[i, :, :]-273.15,
                           rasterized=True,
                           cmap='RdBu_r',
                           norm=mpl.colors.TwoSlopeNorm(vcenter=0,
                                                        vmin=-4,
                                                        vmax=20))

    # vort = axis.pcolormesh(lon2d,lat2d, upper_vars.vo[i,0,:,:]*1e4,
    #                        rasterized=True, cmap='RdBu_r', vmax=1, vmin=-2)

box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()

ax[-1, 0].plot([], [], color='k', label='500hPa Geopotential\nHeight (m)')
ax[-1, 0].plot([], [], color='k', alpha=0.25,
               label='Mean Sea Level Pressure (mb)')
ax[-1, 0].legend(frameon=False, loc=(0, -.5), fontsize=14, ncol=2)

cax = fig.add_axes([box2.xmax*1.03, box2.ymin, 0.01, box1.ymax-box2.ymin])
fig.colorbar(temp, cax=cax, label='900hPa Temperature (°C)',
             ticks=[-4, 0, 4, 8, 12, 16, 20])
plt.savefig('plots/caseofstudy_Aug2013/synoptic_generaldynamics.pdf', dpi=150,
            bbox_inches='tight')

# %%
# =============================================================================
# dynamics, jet stream
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
    axis.set_title(times[i], loc='left')


for i, axis in enumerate(ax.ravel()):

    CL = axis.contour(lon2d, lat2d, upper_vars.sel(level=300,
                                                   time=days).z[i, :, :]/9.8,
                      colors='k',
                      levels=np.arange(8500, 9600, 100), alpha=0.5)
    axis.clabel(CL, [8800, 9500], fmt='%i', fontsize=11)

    ws = (upper_vars.sel(level=300, time=days).u)**2
    ws = ws+(upper_vars.sel(level=300, time=days).v)**2
    ws = np.sqrt(ws)
    winds = axis.pcolormesh(lon2d, lat2d, ws[i, :, :]*3.6,
                            rasterized=True,
                            cmap='PuBu',
                            norm=mpl.colors.Normalize(0, 250))

    # vort = axis.pcolormesh(lon2d,lat2d, upper_vars.vo[i,0,:,:]*1e4,
    #                        rasterized=True, cmap='RdBu_r', vmax=1, vmin=-2)

box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()

ax[-1, 0].plot([], [], color='k', label='300hPa Geopotential\nHeight (m)')
ax[-1, 0].legend(frameon=False, loc=(0, -.5), fontsize=14, ncol=2)

cax = fig.add_axes([box2.xmax*1.03, box2.ymin, 0.01, box1.ymax-box2.ymin])
fig.colorbar(winds, cax=cax, label='300hPa Wind Speed (km/h)',
             ticks=np.arange(0, 250+25, 25))
plt.savefig('plots/caseofstudy_Aug2013/synoptic_jetstream.pdf', dpi=150,
            bbox_inches='tight')


# %%
# =============================================================================
# water advection1: IVT magnitude and IVT vectors
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
    axis.set_title(times[i], loc='left')


for i, axis in enumerate(ax.ravel()):
    # CL = axis.contour(lon2d, lat2d, surface_vars.tcrw[i, :, :],
    #                   levels=np.arange(1000, 1033, 3),
    #                   colors='k', alpha=0.25, rasterized=True)
    # axis.clabel(CL, [1006, 1021], fmt='%i', fontsize=11)

    # CL = axis.contour(lon2d, lat2d, upper_vars.z[i, 0, :, :]/9.8, colors='k',
    #                   levels=np.arange(5200, 5850, 80))
    # axis.clabel(CL, [5360, 5760], fmt='%i', fontsize=11)

    pw = axis.contourf(lon2d, lat2d, ivt.sel(time=days)[i, :, :],
                       rasterized=True, cmap='twilight',
                       levels=np.arange(0, 1250+125, 125))

    q = axis.quiver(lon2d, lat2d, surface_vars["p71.162"].sel(time=days)[i, :, :].values,
                    surface_vars["p72.162"].sel(
        time=days)[i, :, :].values, regrid_shape=20,
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

# %%
# =============================================================================
# water advection2: 800hpa winds and TPW
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
    axis.set_title(times[i], loc='left')


for i, axis in enumerate(ax.ravel()):

    pw = axis.contourf(lon2d, lat2d, surface_vars.tcw.sel(time=days)[i, :, :],
                       rasterized=True, cmap='twilight_r',
                       levels=np.arange(0, 40, 3))

    q = axis.quiver(lon2d, lat2d, upper_vars.u.sel(time=days, level=800)[i, :, :].values,
                    upper_vars.v.sel(time=days, level=800)[
        i, :, :].values, regrid_shape=20,
        transform=ccrs.PlateCarree(), color='k')


box1 = ax[0, -1].get_position()
box2 = ax[1, -1].get_position()

cax = fig.add_axes([box2.xmax*1.03, box2.ymin, 0.01, box1.ymax-box2.ymin])
fig.colorbar(pw, cax=cax,
             label='Precipitable Water $(mm/h)$')

plt.savefig('plots/caseofstudy_Aug2013/synoptic_tpwanalysis.pdf', dpi=150,
            bbox_inches='tight')

# %%
# =============================================================================
# ivt coastal profile time vs lat
# =============================================================================


chile = gpd.read_file('datos/vector/cl_continental_geo.shp')
fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(14, 4))
# fig.tight_layout(pad=0.5)

ivt_cut = ivt.sel(lon=-74, method='nearest').interpolate_na(dim='time')
pr_cut = surface_vars.tp.sel(lon=-71, method='nearest')*1000
lat, time = ivt_cut.lat, ivt_cut.time
time2d, lat2d = np.meshgrid(time, lat)

# ax[0].set_ylabel('Latitude (°S)')
ax[0].set_xlabel('2013-Ago', loc='left')
ax[1].set_xlabel('2013-Ago', loc='left')
ax[0].grid(True, ls=":", which='both', axis='x')
ax[1].grid(True, ls=":", which='both', axis='x')

map1 = ax[0].contourf(time2d, lat2d, ivt_cut.T, cmap='twilight',
                      levels=np.arange(0, 1250+125, 125), rasterized=True)
# ax[0].contour(time2d, lat2d, ivt_cut.T, colors='grey',levels=np.arange(0,1250+125,125))

map2 = ax[1].contourf(time2d, lat2d, pr_cut.T, cmap='Blues',
                      rasterized=True)


ax[0].xaxis.set_major_formatter(mpl.dates.DateFormatter('%d'))
ax[0].xaxis.set_major_locator(
    mpl.dates.DayLocator(bymonthday=np.arange(5, 15, 1)))


ax[0].xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
ax[0].xaxis.set_minor_locator(mpl.dates.HourLocator(byhour=[6, 12, 18]))

fig.colorbar(
    map1, ax=ax[0], label='Integrated Water Vapor\nTransport $(kgm^{-1}s^{-1})$')
fig.colorbar(map2, ax=ax[1], label='Precipitation $(mm/h)$')

box = ax[-1].get_position()
ax2 = fig.add_axes([box.xmax*1.1, box.ymin, 0.1, box.ymax-box.ymin])
chile.plot(ax=ax2, color='k', alpha=0.5)
chile.boundary.plot(ax=ax2, color='k', lw=0.5)
ax2.set_ylim(-50, -20)
# ax2.axis('off')
ax2.set_yticklabels([])
ax2.set_xticks([-74])
ax2.set_xticklabels([-74])
ax2.axvline(-74, color='tab:red', ls="--", lw=1.2)


plt.savefig('plots/caseofstudy_Aug2013/atmospheric_river_evolution.pdf', dpi=150,
            bbox_inches='tight')
