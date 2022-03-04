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
interval = slice(datetime.datetime(2013, 8, 4),
                 datetime.datetime(2013, 8, 16))


# %%
# =============================================================================
# RUNOFF DATA IN DIFFERENT BASINS
# =============================================================================

basin_attributes = pd.read_csv('datos/basins_attributes.csv', index_col=0)
basin_attributes = basin_attributes.sort_values(
    by='gauge_name', ascending=False)
runoff = pd.read_csv('datos/runoff_gauges_dataset.csv', index_col=0)
runoff.index = pd.to_datetime(runoff.index)
runoff2 = runoff[basin_attributes.gauge_name]
runoff = runoff[interval]
runoff = runoff.T.dropna(how='all').T
# %%
# =============================================================================
# BASIN POLYGONS
# =============================================================================
paths = glob('datos/vector/basins/*.shp')
polygons = pd.concat([gpd.read_file(p) for p in paths])
polygons.index = polygons.gauge_id
polygons = polygons.loc[basin_attributes.index]
polygons.gauge_name = basin_attributes.gauge_name


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

n = 1250
polygons.plot(column='area_km2',
              ax=ax,
              transform=ccrs.PlateCarree(),
              cmap='BuPu',
              legend=True,
              legend_kwds={'orientation': 'horizontal',
                           'fraction': 0.021,
                           'pad': 0.08,
                           'label': 'Basin Total Area $(km^2)$',
                           'ticks': np.arange(500, 4500+n, n)})


# polygons.plot(polygons.gauge_name, ax=ax,facecolor=None)
gauges = ['Rio Maipo En El Manzano',
          'Rio Cachapoal En Pte Termas De Cauquenes',
          'Rio Teno Despues De Junta Con Claro',
          'Rio Colorado En Junta Con Palos',
          'Rio Achibueno En La Recova',
          'Rio Uble En San Fabian N 2']

polygons.boundary.plot(ax=ax, lw=0.5, color='k', transform=ccrs.PlateCarree())


colors = plt.cm.tab10(np.linspace(0, 1, 10))

# polygons.loc[set2].boundary.plot(ax=ax,edgecolor='tab20')
# polygons.T[sorted(set2)].T.geometry
for i, g in enumerate(gauges):
    gpd.GeoSeries(polygons.loc[g].geometry).boundary.plot(color=colors[i, :],
                                                          ax=ax, zorder=11)
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
    axis.set_xlim('2013-08-04', '2013-08-16')
    axis.set_xticks(runoff.index[:: 48])
    axis.set_xticklabels([])
    axis.tick_params(axis='x', which='major', rotation=45)
    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
    axis.xaxis.set_minor_locator(mpl.dates.HourLocator(interval=6))
    axis.grid(axis="x", which='major', ls=":")
    title = g.split(" ")
    title = " ".join(title[:len(title)//2])+"\n" + \
        " ".join(title[len(title)//2:])
    if i == 1:
        title = 'Rio Cachapoal En\nPte Termas De\nCauquenes'
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

plt.savefig('plots/caseofstudy_Aug2013/runoff_basins.pdf',
            dpi=150, bbox_inches='tight')
