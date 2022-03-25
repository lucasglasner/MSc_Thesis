#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:30:28 2022

@author: lucas
"""

import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from glob import glob

#%%

data1 = xr.open_dataset('datos/modis/snapshot-2008-06-02T00_00_00Z.nc')
data2 = xr.open_dataset('datos/modis/snapshot-2008-06-06T00_00_00Z.nc') 

data1 = data1.sel(lon=slice(-70.7,-70),lat=slice(-34,-33))
data2 = data2.sel(lon=slice(-70.7,-70),lat=slice(-34,-33))

basins = [gpd.read_file(p) for p in glob('datos/vector/basins/mains/*.shp')]
basins = pd.concat(basins)
basin = gpd.read_file('datos/vector/basins/mains/RioMaipoEnElManzano.shp')
bounds = basin.geometry.total_bounds

lat,lon = data1.lat,data1.lon

data11 = np.moveaxis(np.stack([data1.Band1,data1.Band2,data1.Band3]),0,-1)[::-1]
data22 = np.moveaxis(np.stack([data2.Band1,data2.Band2,data2.Band3]),0,-1)[::-1]



dem = xr.open_dataset('datos/topography/Andes_topo_001x001grad.nc')
dem = dem.elevation

#%%
plt.rc('font',size=18)
fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,5),
                      subplot_kw={'projection':ccrs.PlateCarree()})
ax = ax.ravel()
extent = [lon.min(),lon.max(),lat.min(),lat.max()]
for a in ax:
    basins.boundary.plot(ax=a,color='k',lw=0.8)
    a.set_extent(extent)
    # a.axis('off')
    a.set_xticks([])
    a.set_yticks([])
    
# ax[0].set_title('2008-06-02',loc='left')
ax[0].imshow(data11, extent=extent)
ax[1].imshow(data22,extent=extent)
# ax[1].set_title('2008-06-06', loc='left')


# ax[0].plot([],[],label='1990 (m)',color='r')
# ax[0].legend(frameon=False,loc='lower left',
#              fontsize=12, labelcolor='w')

ax[0].contour(dem.lon,dem.lat,dem, levels=[1990], colors='r',
              linewidths=0.7)
ax[1].contour(dem.lon,dem.lat,dem, levels=[1990], colors='r',
              linewidths=0.7)

plt.savefig('plots/caseofstudy_Jun2008/comparacion2.pdf',dpi=150,
            bbox_inches='tight')