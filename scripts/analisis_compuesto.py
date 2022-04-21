#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:38:15 2022

@author: lucas
"""

import xarray as xr
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib as mpl

#%%
b = 'RioUbleEnSanFabianN2'
data_daily = pd.read_csv('datos/TABLAS_DIARIAS_CRECIDAS_CC.csv',index_col=[0,1])
# data_daily = data_daily.loc['RioMaipoEnElManzano']
data_daily = data_daily.loc[b]
data_daily.index = pd.to_datetime(data_daily.index)
data_daily['ros_day'] = (data_daily['ros_area']>0.1) & (data_daily['FL']>data_daily['SL'])


data_events = pd.read_csv('datos/TABLAS_EVENTOS_CRECIDAS.csv',index_col=[0,1])
data_events = data_events.loc[b]
data_events['ros_event'] = data_events['max_ros_area']>0.1


hypso = pd.read_csv('datos/topography/basins/hypso/'+b+'_hypso.csv',
                    index_col=0)

pr = pd.read_csv('datos/estaciones/vismet/pr_'+b+'.csv',
                 index_col=0)['Valor']
pr.index = pd.to_datetime(pr.index)

q = pd.read_csv('datos/estaciones/dga/qinst_'+b+'.csv',
                index_col=0).squeeze()
q.index = pd.to_datetime(q.index)

SCA = pd.read_csv('datos/ianigla/'+b+'_SCA_s_comp.filtro_MA.3días.csv',
                  index_col=0)['SCA(%)']
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA.resample('h').interpolate('cubicspline')
SCA = SCA.where(SCA>0).fillna(0)
SCA = SCA.where(SCA<100).fillna(100)

H0 = pd.read_csv('datos/era5/horarias/H0_'+b+'_2000-2020.csv',
                 index_col=0).deg0l
H0.index = pd.to_datetime(H0.index)

int_func = interp1d(hypso.index,hypso.fArea)
int_func2 = interp1d(hypso.fArea,hypso.index)
SL = pd.Series(int_func2(SCA.values/100),index=SCA.index)
ap = pd.Series(int_func(H0-300).squeeze(),index=H0.index)


basin_data = pd.concat([pr,q,SCA/100,H0,SL,ap],axis=1).dropna()
basin_data.columns = ['pr','q','sca','h0','sl','pluv_area']
basin_data['ros_area'] = basin_data['sca']-(1-basin_data['pluv_area'])
basin_data['ros_area'] = basin_data['ros_area'].where(basin_data['pr']>0)
basin_data['ros_area'] = basin_data['ros_area'].where(basin_data['ros_area']>0.1)
#%%

# times = data_daily[data_daily['ros_day']]
# t1 = times.sort_values(by='ros_area',ascending=False).head(20).ros_area
# t2 = times.SCA_trend.where(times.SCA_trend<0).dropna()
# # t3 = times.sort_values(by='Qmaxd',ascending=False).head(300)
# # t4 = times.sort_values(by='pr_cr2met',ascending=False).tail(150)
# times = data_events[data_events['max_ros_area']>0.1]
# times = times[times.delta_sca<0]
# t1 = times.delta_sca.sort_values().head(100)
# t2 = times.max_ros_area.sort_values(ascending=False).head(15)
# t3 = times.quickflow.sort_values(ascending=False).head(100)
# # times = t1.index
# times = set(t1.index) & set(t2.index) & set(t3.index)
# # times = t1.index
# times = np.array(list(times))

# times = data_events.loc[times].start
# times = pd.to_datetime(times.values)
# print(times.shape)

times = data_daily[data_daily.ros_area>0.1].dropna()
times = times[times.SCA_trend<0]
t1 = times.ros_area.sort_values().tail(15)
t2 = times.quickflow.sort_values().tail(300)
t3 = times.SCA_trend.sort_values().head(300)

times = set(t1.index) & set(t2.index) & set(t3.index)
times = list(times)
times = pd.to_datetime(times)
print(len(times))


#%%


# start = data_events[data_events['ros_event']].start
# start = pd.to_datetime(start)
# start = times
# ros_events = []
# for s in start:
#     b = basin_data.loc[s-pd.Timedelta(days=2):s+pd.Timedelta(days=2)]
#     b = b.reset_index().drop('index',axis=1)
#     if b.ros_area.max()>0.1:
#         ros_events.append(b)


# ros_events = pd.concat(ros_events,keys=range(len(ros_events)))
# ros_comp = ros_events.unstack().mean()
#%%
surface = xr.open_dataset('datos/era5/daily/era5_surface.nc',
                          chunks='auto')
surface = surface.sel(expver=1).reindex({'time':data_daily.index},
                                        method='nearest')

ivtx = surface['p71.162'].compute()
ivty = surface['p72.162'].compute()
surface = surface[['msl','tcw']]
surface_anomaly = surface.groupby('time.dayofyear')-surface.groupby('time.dayofyear').mean(dim='time')
surface_anomaly = surface_anomaly.compute()

msl_anomaly = surface_anomaly.msl
tpw_anomaly = surface_anomaly.tcw

del surface, surface_anomaly

#%%


upper1 = xr.open_dataset('datos/era5/daily/era5upper1.nc',chunks='auto')
upper2 = xr.open_dataset('datos/era5/daily/era5upper2.nc', chunks='auto')


upper1 = upper1.sel(expver=1)
upper2 = upper2.sel(expver=1)

upper = xr.concat([upper2,upper1],dim='time').drop_duplicates(dim='time')
upper = upper.reindex({'time':data_daily.index},method='nearest')
del upper1,upper2

z500 = xr.open_dataset('datos/era5/daily/geopotential.nc',chunks='auto')
z500 = z500.sel(expver=1).z.compute().reindex({'time':data_daily.index},
                                              method='nearest')
t900 = upper.t.sel(level=900)
u800,v800 = upper.u.sel(level=800),upper.v.sel(level=800)
u300,v300 = upper.u.sel(level=300),upper.v.sel(level=300)

del upper

t900_anomaly = t900.groupby('time.dayofyear')-t900.groupby('time.dayofyear').mean(dim='time')
t900_anomaly = t900_anomaly.compute()

z500_anomaly = z500.groupby('time.dayofyear')-z500.groupby('time.dayofyear').mean(dim='time')


u800,v800 = u800.compute(),v800.compute()
u300,v300 = u300.compute(),v300.compute()


#%%
fig,ax = plt.subplots(4,3, sharex=True,sharey=True, figsize=(13,9),
                      subplot_kw={'projection':ccrs.PlateCarree()})

from functions import add_labels
add_labels(ax,yticks=[-25,-35,-45],xticks=[-90,-70],
           linewidth=0)
for axis in ax.ravel():
    axis.set_extent([-100,-60,-20,-50], crs=ccrs.PlateCarree())
    axis.coastlines()
    
plt.rc('font',size=14)

lon,lat = u300.longitude,u300.latitude
lon2d,lat2d = np.meshgrid(lon,lat)
mask = (lat2d>-48) & (lat2d<-25)
for i in range(3):
    t = times+pd.Timedelta(days=i-1)
    ax[0,i].set_title(str(24*(i-1))+"hrs", loc='left')
    cs=ax[0,i].contour(lon2d,lat2d,
                       1e-3*msl_anomaly.sel(time=t).mean('time'),
                       transform=ccrs.PlateCarree(),
                       colors='k', levels=np.arange(-3,3.2,0.2),
                       alpha=0.6)
    # ax[0,i].clabel(cs,cs.levels[::2])
    mapat = ax[0,i].pcolormesh(lon2d,lat2d,
                               t900_anomaly.sel(time=t).mean('time'),
                               rasterized=True, cmap='RdBu_r',
                               norm=mpl.colors.TwoSlopeNorm(0,-4,4))
    # ax[0,i].streamplot(lon2d[mask],lat2d[mask],
    #                    u300.sel(time=t).mean('time').values[mask],
    #                    v300.sel(time=t).mean('time').values[mask],
    #                    transform=ccrs.PlateCarree(),
    #                    color='k', density=0.5,arrowsize=0.5,
    #                    linewidth=0.9)
    ax[1,i].contour(lon2d,lat2d,
                    z500_anomaly.sel(time=t).mean('time')/9.8,
                    colors='k',alpha=0.8,levels=np.arange(-500,530,30))
    
    
    uv = np.sqrt(u300.sel(time=t).mean('time')**2+v300.sel(time=t).mean('time')**2)
    mapawinds = ax[1,i].pcolormesh(lon2d,lat2d,
                                   uv,
                                   cmap='BuPu',
                                   rasterized=True)
    
    ax[1,i]
    
    mapaw = ax[2,i].pcolormesh(lon2d,lat2d,
                               tpw_anomaly.sel(time=t).mean('time'),
                               rasterized=True,
                               cmap='BrBG',
                               norm=mpl.colors.TwoSlopeNorm(0,vmin=-6,vmax=6))
    ax[2,i].quiver(lon2d,lat2d,
                   u800.sel(time=t).mean('time').values,
                   v800.sel(time=t).mean('time').values,
                   transform=ccrs.PlateCarree(),
                   regrid_shape=8)
    
    ivtxx = ivtx.sel(time=t).mean('time')
    ivtyy = ivty.sel(time=t).mean('time')
    ivt = (ivtxx**2+ivtyy)**0.5
    mapa_ivt = ax[3,i].pcolormesh(lon2d,lat2d,
                                  ivt,
                                  cmap='twilight',
                                  rasterized=True)
    ax[3,i].quiver(lon2d,lat2d,
                   ivtxx.values,ivtyy.values,
                   transform=ccrs.PlateCarree(),
                   regrid_shape=8)
            
box = ax[0,-1].get_position()
cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])                  
fig.colorbar(mapat,cax=cax,label='900hPa Temp.\nanomaly (°C)\n MSLP anomaly\n(contours)',
             ticks=[-4,-2,0,2,4])
# cax.set_ylim(-2,4)

box = ax[1,-1].get_position()
cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])                  
fig.colorbar(mapawinds,cax=cax,label='300hPa Winds\n (m/s)\n Z500 anomaly\n(contours)')

box = ax[2,-1].get_position()
cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])
fig.colorbar(mapaw,cax=cax, label='Precipitable\nwater anomaly\n(mm)\n800hPa winds\n(vectors)',
             ticks=[-6,-3,0,3,6])

box = ax[3,-1].get_position()
cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])
fig.colorbar(mapa_ivt,cax=cax, label='Integrated\nwater vapor\ntransport\n($kg m^{-1} s^{-1}$)',
             ticks=[0,100,200,300])

plt.savefig('plots/COMPUESTO_ROS.pdf',dpi=150,bbox_inches='tight')