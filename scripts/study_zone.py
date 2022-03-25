#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:26:51 2022

@author: lucas

# =============================================================================
# This Script does a quick analysis of the study zone and the different domains
# =============================================================================
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import cartopy.feature as cf
import matplotlib as mpl
import numpy as np
import pandas as pd
from glob import glob
import geopandas as gpd
# %%
# =============================================================================
# Define domains
# =============================================================================
# domain0 = [-72, -69, -32.4, -37.3]
domain0 = [-72.1, -69.5, -32.4, -37]
domain1 = [-74, -68, -26, -38]
domain2 = [-100, -60, -20, -50]


def draw_domain(ax, d, text=None, **kwargs):
    ax.plot([d[0], d[0], d[1], d[1], d[0]],
            [d[2], d[3], d[3], d[2], d[2]],
            **kwargs)
    ax.text(d[0]+0.05, d[2]+0.05, text, color='k')
    return None

# %%
# =============================================================================
# data
# =============================================================================


topo = xr.open_dataset('datos/topography/ncar.0.25-deg.nc').data
topo = topo.where(topo > 0).sel(lon=slice(360-120, 360-50), lat=slice(-5, -70))
hrestopo = xr.open_dataset('datos/topography/Andes_topo_005x005grad.nc')
hrestopo = hrestopo.where(hrestopo > 0).elevation

paths = glob('datos/vector/basins/Rio*.shp')
cuencas = []
for path in paths:
    cuencas.append(gpd.read_file(path))
cuencas = pd.concat(cuencas)
cuencas.index = cuencas.gauge_name.map(lambda x: x.replace(" ",""))
cuencas = cuencas.loc[['RioAconcaguaEnChacabuquito',
                       'RioMapochoEnLosAlmendros',
                       'RioMaipoEnElManzano',
                       'RioCachapoalEnPteTermasDeCauquenes',
                       'RioTinguiriricaBajoLosBriones',
                       'RioTenoDespuesDeJuntaConClaro',
                       'RioColoradoEnJuntaConPalos',
                       'RioMauleEnArmerillo',
                       'RioUbleEnSanFabianN2']]

ros = xr.open_dataset('datos/ROS/CORTES_CR2MET_ERA5/ROS_1984.nc')
landcover = xr.open_dataset('datos/landcover/LC_CHILE_2014_b_final.nc').Band1
landcover = landcover.sel(lon=slice(-74, -68), lat=slice(-39, -25))
landcover = landcover[::50, ::50]
# landcover = landcover.reindex({'lat':ros.lat,'lon':ros.lon},method='nearest')
# landcover = landcover.where(landcover!=0)
# %%
# =============================================================================
# Define land cover classes
# =============================================================================
# Land use values
land_uses = {"Cropland": (100, 200),
             "Native Forest": (200, 229),
             "Forest Plantation": (229, 300),
             "Grassland": (300, 400),
             "Shrubland": (400, 500),
             "Wetland": (500, 600),
             "Water Bodies": (600, 800),
             "Waterproofs": (800, 900),
             "Barren": (900, 1000),
             "Ice/Snow": (1000, 1200),
             "Clouds": (1200, 1500)}

# %%
# =============================================================================
# Create figure
# =============================================================================
fig = plt.figure(figsize=(5, 5))
plt.rc('font', size=18)


# =============================================================================
# map with synoptic domain
# =============================================================================
ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
ax.set_extent([-105, -55, -15, -55])
ax.coastlines('50m', rasterized=True)
draw_domain(ax, domain0, lw=2, color='k', text='d1',
            transform=ccrs.PlateCarree())
draw_domain(ax, domain1, lw=2, color='k', text='d2',
            transform=ccrs.PlateCarree())
draw_domain(ax, domain2, lw=2, color='k', text='d3',
            transform=ccrs.PlateCarree())
ax.add_feature(cf.BORDERS, rasterized=True, ls=":")
ax.add_feature(cf.OCEAN, rasterized=True)

lon2d, lat2d = np.meshgrid(topo.lon, topo.lat)
cmaplist = pd.read_csv("terraincolormap.txt", header=None).values
cmaplist = [list(np.array(i)/255) for i in cmaplist]
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, len(cmaplist))

map1 = ax.contourf(lon2d, lat2d, np.clip(topo.values.squeeze(),
                                         0, 5e3),
                   rasterized=True,
                   transform=ccrs.PlateCarree(), cmap=cmap,
                   norm=mpl.colors.Normalize(0, 5e3), levels=30)
box = ax.get_position()
cax = fig.add_axes([box.xmin, box.ymin-0.15, box.xmax-box.xmin, 0.05])
cb = fig.colorbar(map1, cax=cax, orientation='horizontal',
                  norm=mpl.colors.Normalize(0, 5e3),
                  ticks=np.arange(.75e3, 5.5e3, .75e3),
                  label='Orography (m.a.s.l)')


gl = ax.gridlines(linestyle=":", draw_labels=True, color='k')
gl.right_labels = False
gl.top_labels = False
gl.xlocator = mpl.ticker.FixedLocator([-100, -80, -60, -40])
# gl.ylocator = mpl.ticker.FixedLocator([-37, -36, -35, -34, -33])
# =============================================================================
# maps with nested domains
# =============================================================================
ax1 = fig.add_axes([1, 0, 1, 1], projection=ccrs.PlateCarree())
ax1.set_extent([-74.5, -67.5, -25.5, -38.5])
ax1.coastlines('10m')
cuencas.boundary.plot(ax=ax1, transform=ccrs.PlateCarree(), lw=0)
ax1.add_feature(cf.BORDERS, rasterized=True, ls=":")
ax1.add_feature(cf.OCEAN, rasterized=True)


lon2d, lat2d = np.meshgrid(hrestopo.lon, hrestopo.lat)
map1 = ax1.pcolormesh(lon2d, lat2d, np.clip(hrestopo.values.squeeze(),
                                            0, 5e3),
                      rasterized=True,
                      transform=ccrs.PlateCarree(), cmap=cmap,
                      norm=mpl.colors.Normalize(0, 5e3))

gl = ax1.gridlines(linewidth=0, draw_labels=True)
gl.right_labels = False
gl.top_labels = False
draw_domain(ax1, domain1, lw=2, color='k', text='d2')
draw_domain(ax1, domain0, lw=2, color='k', text='d1')


cmap = plt.cm.tab10
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap',
                                                    cmaplist, cmap.N)
# define the bins and normalize
bounds = [land_uses[name][0] for name in land_uses.keys()]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

ax11 = fig.add_axes([1.55, 0, 1, 1], projection=ccrs.PlateCarree())
# ax11.sharex(ax1)
# ax11.sharey(ax1)
ax11.set_extent([-72.1, -69.5, -32.4, -37], crs=ccrs.PlateCarree())
ax11.coastlines('10m')
ax11.add_feature(cf.BORDERS, rasterized=True, ls=":")
ax11.add_feature(cf.OCEAN, rasterized=True)
ax11.add_feature(cf.LAND, color='k', alpha=0.2, rasterized=True)

lon2d, lat2d = np.meshgrid(landcover.lon, landcover.lat)
map1 = ax11.pcolormesh(lon2d, lat2d, landcover.where(landcover != 0),
                       rasterized=True, cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       norm=norm)
box = ax11.get_position()
cax = fig.add_axes([box.xmax+0.05, box.ymin, 0.05, box.ymax-box.ymin])
ticks = [np.mean([land_uses[name][1], land_uses[name][0]])
         for name in land_uses.keys()]
cb = fig.colorbar(map1, cax=cax, ticks=ticks)
cb.set_ticklabels(list(land_uses.keys()))

cuencas.boundary.plot(ax=ax11, transform=ccrs.PlateCarree(), color='k')
gl = ax11.gridlines(linewidth=0, draw_labels=True)
gl.left_labels = False
gl.right_labels = False
gl.top_labels = False


# draw_domain(ax1,domain, lw=2, color='k')

# =============================================================================
# inset
# =============================================================================
box = ax.get_position()
ax2 = fig.add_axes([box.xmin-0.25, box.ymax-0.15, 0.3, 0.3],
                   projection=ccrs.Orthographic(-90, -20))

ax2.coastlines('110m')
ax2.set_global()
ax2.add_feature(cf.LAND, rasterized=True, color='k', alpha=0.2)
ax2.add_feature(cf.OCEAN, rasterized=True)
ax2.gridlines(linestyle=":", rasterized=True, color='k')
# =============================================================================
# Common stuff
# =============================================================================

ax.set_title('(a)', loc='right')
ax1.set_title('(b)', loc='right')
ax11.set_title('(c)', loc='right')

plt.savefig('plots/ZONA_DE_ESTUDIO.pdf', dpi=150, bbox_inches='tight')
