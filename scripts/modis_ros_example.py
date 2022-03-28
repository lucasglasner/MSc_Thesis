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

data = [xr.open_dataset(p) for p in glob('datos/modis/snapshot-2008-06-*.nc')]
# data4 = xr.open_dataset('datos/modis/snapshot-2008-06-06T00_00_00Z.nc') 
data = [d.sel(lon=slice(-70.8,-70),lat=slice(-33.9,-33)) for d in data]

data00 = [np.moveaxis(np.stack([d.Band1,d.Band2,d.Band3]),0,-1)[::-1] for d in data]
data1,data2,data3,data4 = data


basins = [gpd.read_file(p) for p in glob('datos/vector/basins/mains/*.shp')]
basins = pd.concat(basins)
basin = gpd.read_file('datos/vector/basins/mains/RioMaipoEnElManzano.shp')
bounds = basin.geometry.total_bounds

lat,lon = data1.lat,data1.lon

data11,data22,data33,data44 = data00

dem = xr.open_dataset('datos/topography/Andes_topo_001x001grad.nc')
dem = dem.elevation

#%%
plt.rc('font',size=18)
fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(8,3),
                      subplot_kw={'projection':ccrs.PlateCarree()})
ax = ax.ravel()
extent = [lon.min(),lon.max(),lat.min(),lat.max()]
for a in ax:
    basins.boundary.plot(ax=a,color='k',lw=0.8)
    a.set_extent(extent)
    # a.axis('off')
    a.set_xticks([])
    a.set_yticks([])
    a.contour(dem.lon,dem.lat,dem, levels=[1990], colors='r',
              linewidths=0.7)
    
# ax[0].set_title('2008-06-02',loc='left')
# ax[1].set_title('2008-06-06',loc='left')
# ax[2].set_title('2008-06-05',loc='left')
# ax[3].set_title('2008-06-06',loc='left')

ax[0].imshow(data11, extent=extent)
ax[1].imshow(data44, extent=extent)
# ax[2].imshow(data33, extent=extent)
# ax[3].imshow(data44,extent=extent)
# ax[1].set_title('2008-06-06', loc='left')


# ax[0].plot([],[],label='1990 (m)',color='r')
# ax[0].legend(frameon=False,loc='lower left',
              # fontsize=12, labelcolor='w')



plt.savefig('plots/caseofstudy_Jun2008/comparacion1.pdf',dpi=150,
            bbox_inches='tight')