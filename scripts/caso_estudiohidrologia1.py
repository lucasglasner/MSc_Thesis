#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:27:42 2022

@author: lucas

# =============================================================================
# This Script does some hydrological analysis in the case of study
# =============================================================================
"""

from functions import add_labels
from glob import glob
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import datetime
from scipy.interpolate import interp1d
import scipy.stats as st
from scipy.ndimage.filters import minimum_filter1d, generic_filter
from scipy.ndimage.measurements import label
from scipy.signal import argrelextrema
import xarray as xr
import cartopy.crs as ccrs
import geopandas as gpd
import cartopy.feature as cf


# %%
# =============================================================================
# big time interval and graph time interval
# =============================================================================

interval = slice(datetime.datetime(2005, 10, 15),
                 datetime.datetime(2005, 10, 25))
interval2 = slice(datetime.datetime(2015, 7, 28),
                  datetime.datetime(2015, 8, 10))


# %%
# =============================================================================
# RUNOFF DATA IN DIFFERENT BASINS
# =============================================================================

basin_attributes = pd.read_csv('datos/basins_attributes.csv', index_col=0)
basin_attributes = basin_attributes.sort_values(
    by='gauge_name', ascending=True)
runoff = pd.read_csv('datos/runoff_gauges_dataset.csv', index_col=0)
runoff.index = pd.to_datetime(runoff.index)
runoff2 = runoff[basin_attributes.gauge_name]
runoff = runoff[interval]
runoff = runoff.T.dropna(how='all').T


# %%
q95 = np.empty(runoff2.shape[1])
for i, c in enumerate(runoff2.columns):
    q95[i] = np.percentile(runoff2[c].dropna(), 90)

q95 = pd.Series(q95, index=runoff2.columns)
# %%
# =============================================================================
# PRECIPITATION, FREEZING LEVEL AND HYPSOMETRIC DATA
# =============================================================================

pr = pd.read_csv('datos/pr_cr2met_mainbasins.csv', index_col=0)
pr.index = pd.to_datetime(pr.index)
pr = pr[interval2]
pr.columns = [int(p) for p in pr.columns]
pr = pr[basin_attributes.index]

paths = glob('datos/topography/basins/hypso/*.csv')
hypso = [pd.read_csv(p, index_col=0) for p in paths]
hypso = pd.concat(hypso, keys=basin_attributes.index,
                  axis=1)


areas = [hypso[b].Area_km2.max() for b in pr.columns]
areas = pd.Series(areas, index=pr.columns)
int_func = [interp1d(hypso[b].index, hypso[b].fArea) for b in pr.columns]
# == == == == == == == == == == == == == == == == == == == == == == == == == == == ==

# %%
H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)-datetime.timedelta(hours=4)
H0_mm = H0_mm.resample('d').mean().reindex(pr.index)

pluv_area = [int_func[i](H0_mm-300)*areas.values[i] for i in range(len(areas))]
# nonpluv_area = [(1-int_func[i](H0_mm-300))*areas[i] for i in range(3)]
pluv_area = pd.DataFrame(pluv_area, columns=pr.index, index=areas.index).T
# nonpluv_area = pd.DataFrame(nonpluv_area, columns=pr.index, index=pr.columns).T

max_pluv_area = pluv_area.where(pr > 0.1).max()
# max_pluv_area.index = basin_attributes.gauge_name
# %%
# =============================================================================
# BASIN POLYGONS
# =============================================================================
paths = glob('datos/vector/basins/*.shp')
polygons = pd.concat([gpd.read_file(p) for p in paths])
polygons.index = polygons.gauge_id
polygons = polygons.loc[basin_attributes.index]
polygons.gauge_name = basin_attributes.gauge_name

polygons['max_pluv_area'] = (max_pluv_area/areas).values


# %%
plt.rc('font', size=18)
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-72.1, -69.5, -32.4, -37])
ax.coastlines()
ax.add_feature(cf.BORDERS, rasterized=True)
ax.add_feature(cf.OCEAN, rasterized=True)
ax.add_feature(cf.LAND, rasterized=True)

polygons.index = polygons.gauge_name

n = 0.25
polygons.plot(column='max_pluv_area',
              ax=ax,
              transform=ccrs.PlateCarree(),
              cmap='BuPu_r',
              legend=True,
              legend_kwds={'orientation': 'horizontal',
                           'fraction': 0.021,
                           'pad': 0.08,
                           'label': 'Maximum pluvial\narea during rain $(-)$',
                           'ticks': [0.3, 0.5, 0.7, 0.9]})


# polygons.plot(polygons.gauge_name, ax=ax,facecolor=None)
gauges = ['Rio Maipo En El Manzano',
          'Rio Teno Despues De Junta Con Claro',
          'Rio Colorado En Junta Con Palos',
          'Rio Claro En Camarico',
          'Rio Achibueno En La Recova',
          'Rio Uble En San Fabian N 2']

polygons.boundary.plot(ax=ax, lw=0.5, color='k', transform=ccrs.PlateCarree())


colors = plt.cm.tab10(np.linspace(0, 1, 10))

# polygons.loc[set2].boundary.plot(ax=ax,edgecolor='tab20')
# polygons.T[sorted(set2)].T.geometry
for i, g in enumerate(gauges):
    gpd.GeoSeries(polygons.loc[g].geometry).boundary.plot(color=colors[i, :],
                                                          ax=ax, zorder=11,
                                                          lw=2)
    # polygons.boundary.loc[g].plot(ax=ax, color=colors[i,:],zorder=10,
    # transform=ccrs.PlateCarree())


gl = ax.gridlines(linestyle=":", draw_labels=True)


box = ax.get_position()
box_h = box.ymax-box.ymin
box_w = box.xmax-box.xmin


ax0 = fig.add_axes([box.xmin-0.4, box.ymin-0.05, 0.3, 0.15])
ax1 = fig.add_axes([box.xmin-0.4, box.ymin+0.275, 0.3, 0.15])
ax2 = fig.add_axes([box.xmin-0.4, box.ymin+0.6, 0.3, 0.15])

ax3 = fig.add_axes([box.xmax+0.11, box.ymin-0.05, 0.3, 0.15])
ax4 = fig.add_axes([box.xmax+0.11, box.ymin+0.275, 0.3, 0.15])
ax5 = fig.add_axes([box.xmax+0.11, box.ymin+0.6, 0.3, 0.15])

axes = [ax2, ax1, ax0, ax5, ax4, ax3]
del ax0, ax1, ax2, ax3, ax4, ax5


# # colors = pd.DataFrame(colors,index=runoff.columns)
for i, g in enumerate(gauges):
    axis = axes[i]
    axis.plot(runoff[g], color=colors[i, :], lw=2)
    axis.set_xlim(interval.start, interval.stop)
    axis.set_xticks(runoff.index[:: 48])
    axis.set_xticklabels([])
    axis.tick_params(axis='x', which='major', rotation=45)
    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
    axis.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=6))
    axis.grid(axis="x", which='major', ls=":")
    title = g.split(" ")
    title = " ".join(title[:len(title)//2])+"\n" + \
        " ".join(title[len(title)//2:])
    # if i == 1:
    # title = 'Rio Cachapoal En\nPte Termas De\nCauquenes'
    axis.set_title(title, loc='left', fontsize=17)
    if i > 2:
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position("right")


axes[2].xaxis.set_major_formatter(mpl.dates.DateFormatter('%m-%d'))
axes[2].xaxis.set_major_locator(mpl.dates.DayLocator(interval=2))


axes[2].xaxis.set_major_formatter(mpl.dates.DateFormatter('%m-%d'))
axes[2].xaxis.set_major_locator(mpl.dates.DayLocator(interval=2))

axes[-1].xaxis.set_major_formatter(mpl.dates.DateFormatter('%m-%d'))
axes[-1].xaxis.set_major_locator(mpl.dates.DayLocator(interval=2))

axes[1].set_ylabel('Runoff $(m^3/s)$')
axes[-2].set_ylabel('Runoff $(m^3/s)$')
# ax1.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
# axes[2].xaxis.set_minor_locator(mpl.dates.HourLocator(byhour=np.arange(0, 24, 6)))
# axes.tick_params(axis='x', which='major', rotation=45)

# plt.savefig('plots/caseofstudy_Jun2010/runoff_basins.pdf',
#             dpi=150, bbox_inches='tight')
