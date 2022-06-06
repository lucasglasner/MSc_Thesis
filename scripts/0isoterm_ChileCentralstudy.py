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
from sklearn.decomposition import PCA
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
dem = xr.open_dataset(
    'datos/topography/ERA5_orography_005x005grad.nc').elevation
dem = dem.reindex({'lat': H0.lat, 'lon': H0.lon}, method='nearest').load()

# dem_era5 = xr.open_dataset(
#     'datos/topography/ERA5_orography_005x005grad.nc').elevation
# dem_era5 = dem_era5.reindex({'lat': H0.lat, 'lon': H0.lon},
#                             method='nearest').load()
H0 = (H0+dem).where(H0 >= 300).load()
# H0 = (H0+dem).load()

pr_cr2met = xr.open_mfdataset(
    glob('datos/cr2met/CR2MET_pr*'), chunks='auto').pr


H0_stodomingo = pd.read_csv('datos/isotermas0_maipomanzano.csv',
                            index_col=0)
H0_stodomingo = H0_stodomingo['STODOMINGO']
H0_stodomingo.index = pd.to_datetime(H0_stodomingo.index)
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


H0_coast = []
coords = []
for lat, lon in zip(coast.lat, coast.lon):
    h0 = H0.sel(lat=lat, lon=lon, method='nearest')
    lt, ln = h0.lat.item(), h0.lon.item()
    coords.append((lt, ln))
    H0_coast.append(h0.to_series())
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
# coast-montain h0 cuts for sto domingo
# =============================================================================

h0cuts = [H0[:, 87, :].sel(time=t).to_series()
          for t in H0.time]
h0cuts = pd.concat(h0cuts, axis=1)
h0cuts.columns = H0.time.values
h0cuts = h0cuts.T
# h0cuts = h0cuts.iloc[:,maskmont.values]

mask = pr_cr2met.sel(lat=-33.64, method='nearest').to_dataframe()
mask = mask.unstack()['pr'] > 3
h0cuts_pr = h0cuts[mask.values]


# %%
# =============================================================================
# Plot H0 coast profile and coast-montain profile
# =============================================================================

# figure
fig = plt.figure(figsize=(16, 4))
plt.rc('font', size=18)
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122, sharey=ax)

ax.text(0,1.04,'(a)',transform=ax.transAxes)
ax1.text(0.94,1.04,'(b)',transform=ax1.transAxes)
# coast profile
ax.plot(coast.lat, H0_coast.mean(axis=0).rolling(5, center=True).mean(),
        color='tab:red', label='All Days')
H0_prlat = [H0_coast.iloc[:, i][(pr_andes.iloc[:, i] > 3).values].mean()
            for i in range(len(H0_coast.columns))]
ax.plot(coast.lat, pd.Series(H0_prlat).rolling(5, center=True).mean(),
        color='purple', label='Precipitation Days')
ax.fill_between(andes.lat, andes.elevation.rolling(5, center=True).mean(),
                y2=0,
                color="grey", zorder=0, alpha=0.8, lw=0)
ax.plot(andes.lat, andes.elevation.rolling(5, center=True).mean(),
        color='k', alpha=0.5, lw=1, zorder=1)
ax.plot([], [], color='grey', label='Andes Elevation')

ax.scatter(-33.64, H0_stodomingo.mean(), color='tab:red',
           edgecolor='k', zorder=100, s=100)
ax.scatter(-33.64, H0_stodomingo.reindex(pr_andes.index,
                                         method='nearest')[pr_andes.iloc[:, 57] > 3].mean(),
           edgecolor='k', color='blueviolet', zorder=100, s=100)
ax.scatter([], [], color='white', edgecolor='k',
           label='Sto. Domingo\nRadiosonde')
ax.legend(fontsize=12, frameon=False, loc=(0, 1), ncol=2)
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
ax.set_ylim(0, 6e3)
ax.axvline(-33.64, ls=":", color='k')
ax.set_xticks(target_lats)
ax.set_xticklabels(target_lats)
ax.set_ylabel('Zero Degree Level (m)')
ax.set_xlabel('Latitude along the coastline (°S)')
ax.set_xlim(-26, -38)

# coast montain
target_lons = [-73.5, -72.4, -71, -68.5]
ax1.fill_between(dem.lon, dem.sel(lat=-33.64, method='nearest'),
                 color='grey', zorder=10, lw=0, alpha=0.8)
ax1.plot(dem.lon, dem.sel(lat=-33.64, method='nearest'),
         color='k', lw=1, alpha=0.5)
# h0cuts.T[::1000].plot(ax=ax1, color='tab:red', alpha=0.5, legend=False)
mask = (h0cuts.columns > -70.3) & (h0cuts.columns < -69.25)
# ax1.fill_between(h0cuts.columns,
#                  h0cuts.mean(axis=0).where(~mask)-h0cuts.std(axis=0).where(~mask),
#                  h0cuts.mean(axis=0).where(~mask)+h0cuts.std(axis=0).where(~mask),
#                  color='tab:red', alpha=0.25)
# ax1.fill_between(h0cuts.columns,
#                  h0cuts_pr.mean(axis=0).where(~mask)-h0cuts_pr.std(axis=0).where(~mask),
#                  h0cuts_pr.mean(axis=0).where(~mask)+h0cuts_pr.std(axis=0).where(~mask),
#                  color='blueviolet', alpha=0.25)
ax1.plot(h0cuts.columns,
         h0cuts.mean(axis=0).where(~mask),
         color='tab:red', lw=2)
ax1.plot(h0cuts.columns,
         h0cuts_pr.mean(axis=0).where(~mask),
         color='blueviolet', lw=2)
ax1.scatter(-71.62, H0_stodomingo.mean(), color='tab:red',
            edgecolor='k', zorder=100, s=100)
ax1.scatter(-71.62, H0_stodomingo.reindex(pr_andes.index,
                                          method='nearest')[pr_andes.iloc[:, 57] > 3].mean(),
            edgecolor='k', color='blueviolet', zorder=100, s=100)
for i, lon in enumerate(target_lons):
    pos = np.argwhere(np.diff(np.sign(h0cuts.columns.values-lon)))
    v = h0cuts.iloc[:, pos.item()].dropna()
    ax1.boxplot(v, positions=[target_lons[i]+0.1], widths=0.125,
                sym="", patch_artist=True, medianprops={'color': 'k'},
                boxprops={'facecolor': 'tab:red'}, zorder=11)

    pos = np.argwhere(np.diff(np.sign(h0cuts_pr.columns.values-lon)))
    v = h0cuts_pr.iloc[:, pos.item()].dropna()
    ax1.boxplot(v, positions=[target_lons[i]-0.1], widths=0.125,
                sym="", patch_artist=True, medianprops={'color': 'k'},
                boxprops={'facecolor': 'blueviolet'}, zorder=11)

ax1.set_ylim(0, 6e3)
ax1.set_xlim(-74, -68)
ax1.set_xlabel('Longitude along the 33.64°S parallel (°W)')
# ax1.grid(axis='y',ls=":")
ax11 = ax1.twinx()
ax11.set_yticks(ax1.get_yticks())
# ax1.set_yticklabels([])
ax1.set_xticks([-74, -73, -72, -71, -70, -69, -68])
ax1.set_xticklabels([-74, -73, -72, -71, -70, -69, -68])

# inset
box = ax.get_position()
ax2 = fig.add_axes([box.xmax-0.015, box.ymin+0.45, 0.1, box.ymax-box.ymin],
                   projection=ccrs.Orthographic(-70, -30))
ax2.gridlines(linestyles=":")
ax2.coastlines(lw=0.2)
ax2.set_global()
ax2.plot(coast.lon, coast.lat, transform=ccrs.PlateCarree(),
         color='tab:red', linewidth=2)
ax2.add_feature(cf.LAND)
ax2.add_feature(cf.OCEAN)


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
    {'time': pd.date_range('1984-01-01T12:00:00', '2015-12-31T12:00:00')}) > 3)
scycle_rainy = scycle_rainy.assign_coords(season=seasons)
scycle_rainy = scycle_rainy.groupby('season').mean().compute()
scycle_rainy = scycle_rainy.reindex({'lon': scycle.lon})

# %%

fig, ax = plt.subplots(1, 4, subplot_kw={'projection': ccrs.PlateCarree()},
                       figsize=(6, 5))

lon2d, lat2d = np.meshgrid(scycle.lon, scycle.lat)
for axis in ax.ravel():
    axis.coastlines()
    axis.add_feature(cf.BORDERS, ls=":")
    axis.set_xlim(-74, -68)
    axis.set_ylim(H0.lat.min(), H0.lat.max())

pos = [0, 3, 1, 2]
for i in range(4):
    p1 = ax[i].pcolormesh(lon2d, lat2d, scycle_rainy[pos[i], :, :],
                          rasterized=True,
                          shading='auto',
                          cmap='summer',
                          norm=colors.Normalize(2.5e3, 6e3),
                          linewidth=0,
                          antialiased=False)
    p1.set_edgecolor('face')
    ax[i].set_title(scycle.season[pos[i]].item())
    gl = ax[i].gridlines(draw_labels=True, linewidth=0)
    gl.right_labels = False
    gl.top_labels = False
    gl.left_labels = False
    gl.xlocator = mticker.FixedLocator([-70])

gl = ax[0].gridlines(draw_labels=True, linewidth=0)
gl.right_labels = False
gl.top_labels = False
gl.xlocator = mticker.FixedLocator([-70])
gl.ylocator = mticker.FixedLocator([-37, -35, -33, -31, -29, -27])

box1 = ax[-1].get_position()
cax1 = fig.add_axes([box1.xmax+0.05, box1.ymin, 0.02, box1.ymax-box1.ymin])

fig.colorbar(p1, cax=cax1, label='Zero Degree Level\n(Meters Above Sea Level)')

plt.savefig('plots/H0_climatology_precipdays.pdf',
            dpi=150, bbox_inches='tight')

# %% H0 MEAN CUTS

var1 = scycle.mean(dim='lat').to_dataframe(name='hi').unstack()
var2 = scycle_rainy.mean(dim='lat').to_dataframe(
    name='hi').interpolate().unstack()

var1_std = scycle.std(dim='lat').to_dataframe(name='hi').unstack()
var2_std = scycle_rainy.std(dim='lat').to_dataframe(
    name='hi').interpolate().unstack()

fig, ax = plt.subplots(1, 4, sharex=True, figsize=(15, 3), sharey=True)

titles = ['autumn', 'winter', 'spring', 'summer']
for i, s in enumerate(titles):
    ax[i].plot(H0.lon, var1.loc[s], color='tab:red', label='All Days')
    ax[i].plot(H0.lon, var2.loc[s], color='blueviolet',
               label='Precipitation\nDays')
    ax[i].fill_between(H0.lon, var1.loc[s]+var1_std.loc[s],
                       var1.loc[s]-var1_std.loc[s],
                       color='tab:red', lw=2, alpha=0.25)
    ax[i].fill_between(H0.lon, var2.loc[s]+var2_std.loc[s],
                       var2.loc[s]-var2_std.loc[s],
                       color='blueviolet', lw=2, alpha=0.25)
    ax[i].fill_between(H0.lon, dem.mean(dim='lat'), color='grey')
    ax[i].set_title(s)
    ax[i].set_ylim(0, 8e3)
    ax[i].set_xlim(-74, -68)

ax[0].legend(frameon=False, fontsize=16, loc=(0, 1.2))
ax[0].set_ylabel('Zero Degree Level (m)')
fig.text(0.5, -0.09, 'Longitude (°W)', ha='center', va='center')
plt.savefig('plots/H0_meanlat.pdf', dpi=150, bbox_inches='tight')
