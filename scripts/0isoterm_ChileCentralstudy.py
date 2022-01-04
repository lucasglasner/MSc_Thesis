#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 11:35:21 2022

@author: lucas

# =============================================================================
#
# =============================================================================
"""

# %%
import xarray as xr
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import geopandas as gpd
from glob import glob
import cartopy.feature as cf
import datetime as dt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import cmocean as cm
# %%
# =============================================================================
# Load H0, height and precipitation data !!
# =============================================================================
H0 = xr.open_mfdataset(glob('datos/era5/H0_ERA5_*.nc'), chunks='auto').deg0l
dem = xr.open_dataset('datos/topography/Andes_topo_005x005grad.nc').elevation
dem = dem.reindex({'lat': H0.lat, 'lon': H0.lon}, method='nearest')

pr_cr2met = xr.open_mfdataset(
    glob('datos/cr2met/CR2MET_pr*'), chunks='auto').pr
# %%

# =============================================================================
# Load Chilean border data and build coastline and international border vectors
# =============================================================================
chile = gpd.read_file('datos/vector/cl_continental_geo.shp')
chile = gpd.GeoSeries(list(list(chile.geometry)[0]))
coords = pd.DataFrame(list(chile.exterior[0].coords), columns=['lon', 'lat'])
coords = coords[(coords.lat > H0.lat.min().item()) &
                (coords.lat < H0.lat.max().item())]
coords['type'] = np.empty(len(coords))*np.nan
for i in range(len(coords)):
    if coords.lon.iloc[i] > -70.61:
        coords['type'].iloc[i] = 'andes'
    elif coords.lat.iloc[i] < -34.5 and coords.lon.iloc[i] > -72:
        coords['type'].iloc[i] = 'andes'
    else:
        coords['type'].iloc[i] = 'coast'

coast = coords[coords['type'] == 'coast'].reset_index().drop(
    ['index', 'type'], axis=1)
andes = coords[coords['type'] == 'andes'].reset_index().drop(
    ['index', 'type'], axis=1)

del coords, chile

# %%
# =============================================================================
# Grab zero degree level values and elevation for the coastline and mountain
# =============================================================================
andes['elevation'] = [dem.sel(lat=lat, lon=lon, method='nearest').item()
                      for lat, lon in zip(andes.lat, andes.lon)]
andes = andes.drop_duplicates('elevation').reset_index(drop=True)
coast['elevation'] = [dem.sel(lat=lat, lon=lon, method='nearest').item()
                      for lat, lon in zip(coast.lat, coast.lon)]
coast = coast.drop_duplicates('elevation').reset_index(drop=True)
coast_shift = [dem.sel(lat=lat, lon=lon+1, method='nearest').item()
               for lat, lon in zip(coast.lat, coast.lon)]


H0_coast = [H0.sel(lat=lat, lon=lon, method='nearest').to_series()
            for lat, lon in zip(coast.lat, coast.lon)]
H0_coast = pd.concat(H0_coast, axis=1)
H0_coast.columns = list(zip(coast.lat, coast.lon))


pr_andes = [pr_cr2met.sel(lat=lat, lon=lon+1, method='nearest').to_series()
            for lat, lon in zip(coast.lat, coast.lon)]
pr_andes = pd.concat(pr_andes, axis=1)
pr_andes.columns = list(zip(coast.lat, coast.lon))


# %%
# =============================================================================
# grab dataset from some target latitudes along coastline and mountain
# =============================================================================
target_lats = np.arange(-27, -38, -2)
idx_coast = []
idx_andes = []
for lat in target_lats:
    idx_coast.append(np.argwhere(
        np.diff(np.sign(coast.lat-lat))).flatten().item())
    idx_andes.append(np.argwhere(
        np.diff(np.sign(andes.lat-lat))).flatten().item())

# %%
# =============================================================================
# Plot H0 coast profile
# =============================================================================


fig = plt.figure(figsize=(8, 4))
plt.rc('font', size=16)
ax = fig.add_subplot(111)
box = ax.get_position()
ax2 = fig.add_axes([box.xmax-0.1, box.ymin+0.35, 0.2, box.ymax-box.ymin],
                   projection=ccrs.Orthographic(-70, -30))
ax2.gridlines(linestyles=":")
ax2.coastlines(lw=0.2)
# ax2.set_ylim(H0.lat.min(),H0.lat.max())
ax2.set_global()
# ax2.plot(andes.lon,andes.lat,transform=ccrs.PlateCarree(), color='k', linewidth=2)
ax2.plot(coast.lon, coast.lat, transform=ccrs.PlateCarree(),
         color='tab:red', linewidth=2)
ax2.add_feature(cf.LAND)
ax2.add_feature(cf.OCEAN)


ax.plot(coast.lat, H0_coast.mean(axis=0), color='tab:red', label='All Days')
H0_prlat = [H0_coast.iloc[:, i][(pr_andes.iloc[:, i] > 3).values].mean()
            for i in range(len(H0_coast.columns))]
ax.plot(coast.lat, H0_prlat, color='purple', label='Precip. Days')


ax.fill_between(andes.lat, andes.elevation, y2=0,
                color="grey", zorder=0, alpha=0.8, lw=0)
ax.plot(andes.lat, andes.elevation, color='k', alpha=0.5, lw=1, zorder=1)
ax.plot([], [], color='grey', label='Andes Elevation')
ax.legend(fontsize=10, frameon=False, loc='upper left')

for i in range(len(target_lats)):
    v = H0_coast.iloc[:, idx_coast[i]]
    ax.boxplot(v, positions=[target_lats[i]-0.2], sym="", widths=0.25,
               patch_artist=True, medianprops={'color': 'k'},
               boxprops={'facecolor': 'tab:red'}, zorder=3)

    v = H0_coast.iloc[:, idx_coast[i]][(
        pr_andes.iloc[:, idx_coast[i]] > 3).values]
    ax.boxplot(v, positions=[target_lats[i]+0.2], sym="", widths=0.25,
               patch_artist=True, medianprops={'color': 'k'},
               boxprops={'facecolor': 'blueviolet'}, zorder=3)


ax.set_xlim(coast.lat.max(), coast.lat.min())
ax.set_ylim(0, 8e3)
ax.set_xticks(target_lats)
ax.set_xticklabels(target_lats)
ax.set_ylabel('Zero Degree Level\nabove the coastline (m)')
ax.set_xlabel('Latitude (°S)')

plt.savefig('plots/H0_coastprofile.pdf', dpi=150, bbox_inches='tight')

# %%
seasons = H0.time.to_series().index.month


def fmap(month):
    if month in [12, 1, 2]:
        return 'summer'
    elif month in [3, 4, 5]:
        return 'autumn'
    elif month in [6, 7, 8]:
        return 'winter'
    else:
        return 'spring'


seasons = seasons.map(fmap)
scycle = H0.assign_coords(season=seasons)
scycle = scycle.groupby('season').mean().compute()


scycle_rainy = H0.where(pr_cr2met.assign_coords(
    {'time': pd.date_range('1984-01-01T06:00:00', '2015-12-31T06:00:00')}) > 3)
scycle_rainy = scycle_rainy.assign_coords(season=seasons)
scycle_rainy = scycle_rainy.groupby('season').mean().compute()

# %%

fig, ax = plt.subplots(1, 4, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(6, 5))

lon2d, lat2d = np.meshgrid(H0.lon, H0.lat)
for axis in ax.ravel():
    axis.coastlines()
    axis.set_xlim(-72, -68)
    axis.set_ylim(H0.lat.min(), H0.lat.max())

pos = [0, 3, 1, 2]
for i in range(4):
    p1 = ax[i].pcolormesh(lon2d, lat2d, scycle[pos[i], :, :],
                          rasterized=True,
                          cmap='Spectral_r',
                          norm=colors.Normalize(0, 4.5e3))
    ax[i].set_title(scycle.season[pos[i]].item())
    gl = ax[i].gridlines(linestyle=":", draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    gl.left_labels = False
    gl.xlocator = mticker.FixedLocator([-70])

gl = ax[0].gridlines(linestyle=":", draw_labels=True)
gl.right_labels = False
gl.top_labels = False
gl.xlocator = mticker.FixedLocator([-70])
gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])

box1 = ax[-1].get_position()
cax1 = fig.add_axes([box1.xmax+0.05, box1.ymin, 0.02, box1.ymax-box1.ymin])

fig.colorbar(p1, cax=cax1, label='Zero Degree Level\n(Meters Above Surface)')

plt.savefig('plots/H0_climatology.pdf', dpi=150, bbox_inches='tight')
