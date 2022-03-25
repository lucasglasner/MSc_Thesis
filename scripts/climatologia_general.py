#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 10:43:11 2022

@author: lucas
"""

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
from glob import glob
import geopandas as gpd
#%%

domain1 = [-74, -68, -26, -38]
pr = xr.open_mfdataset(glob('datos/cr2met/complete_cr2met/CR2MET_pr_v2.0_*')).pr
pr = pr.sel(lat=slice(domain1[3],domain1[2]),
            lon=slice(domain1[0],domain1[1]))
pr = pr.groupby('time.year').sum().mean(dim='year').load()


temp = xr.open_mfdataset(glob('datos/cr2met/complete_cr2met/CR2MET_t2m*')).t2m
temp = temp.sel(lat=slice(domain1[3],domain1[2]),
                lon=slice(domain1[0],domain1[1]))
temp = temp.mean(dim='time').load()


swe = xr.open_mfdataset(glob('datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_*'))
swemax = swe.max(dim='time').swe
swe = swe.sel(lat=slice(domain1[3],domain1[2]),
              lon=slice(domain1[0],domain1[1]))
swe = swe.SWE.groupby('time.year').sum().mean(dim='year').load()

#%%
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
#%%

data = pd.read_csv('datos/estaciones/dgf/DATOSUTC_2004-2019.csv',index_col=0)
data = data[['9','5']]
data.index = pd.to_datetime(data.index)
pr_dgf = data['9']
pr_dgf = pr_dgf.groupby([pr_dgf.index.month,pr_dgf.index.year]).sum()
pr_dgf = pr_dgf.unstack().T.mean()

t_dgf = data['5']
t_dgf = t_dgf.resample('d').mean()
t_dgf = t_dgf[t_dgf>-40]
t_dgf = t_dgf.groupby(t_dgf.index.month).mean()
sca = pd.read_csv('datos/ianigla/RioMaipoEnElManzano_SCA_s_comp.filtro_MA.3días.csv')
sca.index = pd.to_datetime(sca.fecha)
sca = sca.iloc[:,1]
sca_ac = sca.groupby(sca.index.month).mean()

#%%
import matplotlib as mpl
import cartopy.feature as cf
import cmocean as cm
plt.rc('font',size=18)
fig = plt.figure(figsize=(12,6))
ax0 = fig.add_subplot(141,projection=ccrs.PlateCarree())
ax0.set_extent(domain1)
ax1 = fig.add_subplot(142,projection=ccrs.PlateCarree())
ax1.set_extent(domain1)

ax3 = fig.add_subplot(143,projection=ccrs.PlateCarree())
ax3.set_extent(domain1)
ax4 = fig.add_subplot(244)
ax5 = fig.add_subplot(248)
# ax4 = fig.add_subplot(236)

fig.tight_layout(pad=1.5)
p = ax0.pcolormesh(pr.lon,pr.lat,pr,cmap='twilight',rasterized=True)


listas = ['(a)','(b)','(c)','(d)','(e)']


for i,axis in enumerate([ax0,ax1,ax3,ax4,ax5]):
    axis.text(.86,1.02,listas[i],transform=axis.transAxes)
    if i<3:
        axis.coastlines()
        axis.add_feature(cf.BORDERS)    
        cuencas.boundary.plot(ax=axis,transform=ccrs.PlateCarree(), color='k',
                              lw=0.5)
        gpd.GeoSeries(cuencas.boundary.loc['RioMaipoEnElManzano']).plot(ax=axis,transform=ccrs.PlateCarree(),
                                                                        color='darkred',lw=0.8)
# ax0.add_feature(cf.OCEAN)+
cb = fig.colorbar(p,ax=ax0,orientation='horizontal', shrink=0.8, pad=0.02,
                  ticks=[0,800,1600,2400],
                  label='Mean anual\nprecipitation (mm)')
cb.ax.tick_params(rotation=45)
t = ax1.pcolormesh(temp.lon,temp.lat,temp,cmap='RdBu_r',
                   norm=mpl.colors.TwoSlopeNorm(0),rasterized=True)
cb=fig.colorbar(t,ax=ax1,orientation='horizontal', shrink=0.8,pad=0.02,
                label='Mean anual\ntemperature (°C)')
cb.ax.tick_params(rotation=45)
s = ax3.pcolormesh(swe.lon,swe.lat,swe.where(swemax>100),cmap=cm.cm.ice,
                   norm=mpl.colors.LogNorm(),rasterized=True
                   )
fig.colorbar(s,ax=ax3,orientation='horizontal', shrink=0.8,pad=0.02,
             label='Mean Anual\nSWE (mm)',ticks=[1e2,1e4])

import numpy as np
ax4.bar(np.arange(12),pr_dgf,edgecolor='k')
ax4.set_title('Santiago',loc='left')
ax4.set_ylabel('Precipitation (mm)')
ax44 = ax4.twinx()
ax44.plot(np.arange(12),t_dgf,marker='o',ls="--",color='tab:red')
ax44.set_ylabel('Temperature (°C)')
ax4.set_xticks(np.arange(12))
ax5.set_xticks(np.arange(12))
ax5.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'],
                    fontsize=14)
ax5.tick_params(axis='x',rotation=45)
ax4.set_xticklabels(['']*12)

ax5.plot(np.arange(12),sca_ac)
ax5.set_ylabel('Rio Maipo En\nEl Manzano SCA (%)')
plt.savefig('plots/clima_chilecentral.pdf',dpi=150,bbox_inches='tight')
