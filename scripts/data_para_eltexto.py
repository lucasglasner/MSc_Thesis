#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 09:58:47 2022

@author: lucas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import geopandas as gpd

#%%
# =============================================================================
# maipo
# =============================================================================

H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0)
H0_mm.index = pd.to_datetime(H0_mm.index)-pd.Timedelta(hours=4)

SL_mm = pd.read_csv('datos/snowlimits_maipomanzano.csv',index_col=0)
SL_mm.index = pd.to_datetime(SL_mm.index)

SCA_mm = pd.read_csv('datos/snowcovers_maipomanzano.csv',index_col=0)
SCA_mm.index = pd.to_datetime(SCA_mm.index)

datos_dgf = pd.read_csv('datos/estaciones/dgf/DATOSUTC_2004-2019.csv',index_col=0)
datos_dgf.index = pd.to_datetime(datos_dgf.index)-pd.Timedelta(hours=4)

pr_dgf = datos_dgf.iloc[:,9]

#%%
# =============================================================================
# basins
# =============================================================================
from glob import glob
names = glob('datos/era5/horarias/*.csv')
names = [n.split("/")[-1].split("_")[1] for n in names]
H0_era5 = [pd.read_csv(p,index_col=0) for p in glob('datos/era5/horarias/*.csv')]
H0_era5 = pd.concat(H0_era5,axis=1)
H0_era5.index = pd.to_datetime(H0_era5.index)-pd.Timedelta(hours=4)
H0_era5.columns = names

basins = names
#%%
names = glob('datos/estaciones/dga/qinst_*.csv')
names = [n.split("/")[-1].split("_")[1].split(".")[0] for n in names]
runoff = [pd.read_csv(p,index_col=0) for p in glob('datos/estaciones/dga/qinst_*.csv')]
runoff = pd.concat(runoff,axis=1)
runoff.index = pd.to_datetime(runoff.index)
runoff.columns = names
runoff['RioMauleEnArmerillo'] = runoff.RioMeladoEnElSalto
runoff = runoff[basins]

#%%
names = glob('datos/topography/basins/hypso/*.csv')
hypso = [pd.read_csv(p) for p in names]

names = [n.split("/")[-1].split("_")[0].split(".")[0] for n in names]
hypso = pd.concat(hypso,axis=1,keys=names).T
hypso = hypso.loc[basins]

#%%
from scipy.interpolate import interp1d
int_func = [interp1d(hypso.loc[b].T.height,hypso.loc[b].T.fArea)
            for b in basins]
int_func = pd.Series(int_func,index=basins)

pluv_area = [int_func.loc[b](H0_era5[b]-300) for b in basins]
pluv_area = pd.DataFrame(pluv_area,columns=H0_era5.index,index=basins).T
#%%
names = glob('datos/ianigla/*_SCA_*.csv')
SCA = [pd.read_csv(p,index_col=0).iloc[:,1] for p in names]

names = [n.split("/")[-1].split("_")[0].split(".")[0] for n in names]
SCA = pd.concat(SCA,axis=1,keys=names).T
SCA.columns = pd.to_datetime(SCA.columns)
SCA = SCA.loc[basins].T.dropna()
SCA = SCA.resample('h').interpolate('cubicspline')


#%%
names = glob('datos/estaciones/cr2met/pr*.csv')
pr_cr2met = [pd.read_csv(p,index_col=0) for p in names]
names = [n.split("/")[-1].split("_")[1].split(".")[0] for n in names]
pr_cr2met = pd.concat(pr_cr2met,axis=1)
pr_cr2met.columns = basins
pr_cr2met.index = pd.to_datetime(pr_cr2met.index)
pr_cr2met = pr_cr2met.resample('h').ffill()/24
# pr_cr2met = pr_cr2met.reindex(SCA.index)q
#%%
pluv_area = pluv_area.where(pr_cr2met>10/24)
ros_area = SCA/100-(1-pluv_area)
ros_area = ros_area.where(ros_area>0).where(pr_cr2met>10/24)

#%%
from scipy.ndimage import gaussian_filter

interval = slice("2008-06-02","2008-06-06")

cuenca = gpd.read_file('datos/vector/basins/mains/RioMaipoEnElManzano.shp')
dem = xr.open_dataset('datos/topography/basins/RioMaipoEnElManzano.nc')

fig,ax = plt.subplots(1,1,sharex=True,sharey=True)

# for axis in [ax]:
ax.axis('off')
cuenca.boundary.plot(ax=ax,color='k')



ax.contourf(dem.lon,dem.lat,dem.Band1,colors='royalblue',
               levels=[SL_mm.IANIGLA_HYPSO["2008-06-02"],
                       7000])

ax.contourf(dem.lon,dem.lat,dem.Band1,colors='tab:red',
               levels=[0,
                       H0_mm.squeeze().resample('d').mean()["2008-06-03":"2008-06-04"].mean()-300])

ax.contourf(dem.lon,dem.lat,dem.Band1,colors='cadetblue',
               levels=[SL_mm.IANIGLA_HYPSO["2008-06-02"],
                       H0_mm.squeeze().resample('d').mean()["2008-06-03":"2008-06-04"].mean()-300])


ax.scatter([],[],color='royalblue',marker='s',label='Snow area')
ax.scatter([],[],color='tab:red',marker='s',label='Pluvial area')
ax.scatter([],[],color='cadetblue',marker='s',label='ROS area')
ax.scatter([],[],color='gold',marker='s',label='Melted area')
ax.legend(frameon=False,loc=(1,0.5))


# ax[1].contourf(dem.lon,dem.lat,dem.Band1,colors='royalblue',
#                levels=[SL_mm.IANIGLA_HYPSO["2008-06-06"],
#                        7000])
ax.contourf(dem.lon,dem.lat,dem.Band1,colors='gold',
               levels=[SL_mm.IANIGLA_HYPSO["2008-06-02"],
                       SL_mm.IANIGLA_HYPSO["2008-06-06"]],
               linewidths=0.5)


# ax.contour(dem.lon,dem.lat,dem.Band1,colors='k',
#            levels=[SL_mm.IANIGLA_HYPSO["2008-06-02"],SL_mm.IANIGLA_HYPSO["2008-06-06"],
#                    H0_mm.squeeze().resample('d').mean()["2008-06-03":"2008-06-04"].mean()-300,
#                    ],linewidths=0.2)


# ax[1].contour(dem.lon,dem.lat,dem.Band1,colors='tab:blue',
#                levels=[]],
#                linewidths=0.5, linestyle="--")









