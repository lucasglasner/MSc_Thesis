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
b = 'RioTenoDespuesDeJuntaConClaro'
data_daily = pd.read_csv('datos/TABLAS_DIARIAS_CRECIDAS_CC.csv',index_col=[0,1])
# data_daily = data_daily.loc['RioMaipoEnElManzano']
data_daily = data_daily.loc[b]
data_daily.index = pd.to_datetime(data_daily.index)
data_daily['ros_day'] = (data_daily['ros_area']>0.1) & (data_daily['FL']>data_daily['SL'])
#%%

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
SCA_trend = (SCA.shift(-2)-SCA.shift(2))/4
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
basin_data['SCA_trend'] = SCA_trend.resample('h').fillna('ffill')

data_daily['ros_area'] = basin_data['ros_area'].resample('d').mean()
#%%

# ros_times = data_daily[data_daily['ros_day']]
# t1 = ros_times.sort_values(by='ros_area',ascending=False).head(20).ros_area
# t2 = ros_times.SCA_trend.where(ros_times.SCA_trend<0).dropna()
# # t3 = ros_times.sort_values(by='Qmaxd',ascending=False).head(300)
# # t4 = ros_times.sort_values(by='pr_cr2met',ascending=False).tail(150)
# ros_times = data_events[data_events['max_ros_area']>0.1]
# ros_times = ros_times[ros_times.delta_sca<0]
# t1 = ros_times.delta_sca.sort_values().head(100)
# t2 = ros_times.max_ros_area.sort_values(ascending=False).head(15)
# t3 = ros_times.quickflow.sort_values(ascending=False).head(100)
# # ros_times = t1.index
# ros_times = set(t1.index) & set(t2.index) & set(t3.index)
# # ros_times = t1.index
# ros_times = np.array(list(ros_times))

# ros_times = data_events.loc[times].start
# times = pd.to_datetime(times.values)
# print(times.shape)



#%%


# start = data_events[data_events['ros_event']].start
# start = pd.to_datetime(start)
# start = ros_times
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


# %%

rain_times = data_daily[data_daily.pr_cr2met>50]
ros_times = rain_times[rain_times.ros_area>0.1]
ros_times = ros_times[ros_times.delta_ros>0]
ros_times = ros_times.sort_values(by='delta_ros')
ros_times = ros_times.tail(20)

rain_times = rain_times[(rain_times.ros_area<0.1)|(np.isnan(rain_times.ros_area))]
rain_times = rain_times.sort_values(by='delta_ros').head(20)

#%%

ros_times = ros_times.index
rain_times = rain_times.index
print(ros_times.shape)
print(rain_times.shape)
#%%
        
from scipy.ndimage import gaussian_filter
fig,ax = plt.subplots(2,4, sharex=True,sharey=True, figsize=(14,4),
                      subplot_kw={'projection':ccrs.PlateCarree()})

from functions import add_labels
add_labels(ax,yticks=[-25,-35,-45],xticks=[-90,-70],
            linewidth=0)
for axis in ax.ravel():
    axis.set_extent([-100,-60,-20,-50], crs=ccrs.PlateCarree())
    axis.coastlines()
    
plt.rc('font',size=18)

lon,lat = u300.longitude,u300.latitude
lon2d,lat2d = np.meshgrid(lon,lat)

vmin=[(-3,0,-5,0),(-3,0,-5,0)]
vmax=[(6,60,10,500),(6,60,10,500)]
ticks=[([-3,0,3,6],[0,20,40,60],[-5,0,5,10],[0,250,500]),
        ([-3,0,3,6],[0,20,40,60],[-5,0,5,10],[0,250,500])]

for i,t in enumerate([ros_times,rain_times]):

    cs=ax[i,0].contour(lon2d,lat2d,
                       1e-3*msl_anomaly.sel(time=t).mean('time'),
                       transform=ccrs.PlateCarree(),
                       colors='k', levels=np.arange(-3,3.2,0.2),
                       alpha=0.6)
    # ax[i,0].clabel(cs,cs.levels[::2])

    mapat = ax[i,0].pcolormesh(lon2d,lat2d,
                               t900_anomaly.sel(time=t).mean('time'),
                               rasterized=True, cmap='RdBu_r',
                               norm=mpl.colors.TwoSlopeNorm(0,vmin[i][0],vmax[i][0]),
                               shading='auto')
    fig.colorbar(mapat, ax=ax[i,0], ticks=ticks[i][0])
    
    
    
    ax[i,1].contour(lon2d,lat2d,
                    z500_anomaly.sel(time=t).mean('time')/9.8,
                    colors='k',alpha=0.8,levels=np.arange(-500,530,20))
    
    
    uv = np.sqrt(u300.sel(time=t).mean('time')**2+v300.sel(time=t).mean('time')**2)
    mapawinds = ax[i,1].pcolormesh(lon2d,lat2d,
                                   uv,vmin=vmin[i][1],vmax=vmax[i][1],
                                   cmap='BuPu',
                                   rasterized=True,
                                   shading='auto')
    fig.colorbar(mapawinds,ax=ax[i,1], ticks=ticks[i][1])

    mapaw = ax[i,2].pcolormesh(lon2d,lat2d,
                                tpw_anomaly.sel(time=t).mean('time'),
                                rasterized=True,
                                cmap='BrBG',
                                norm=mpl.colors.TwoSlopeNorm(0,vmin[i][2],vmax[i][2]),
                                shading='auto')
    fig.colorbar(mapaw, ax=ax[i,2], ticks=ticks[i][2])
    ax[i,2].quiver(lon2d,lat2d,
                    u800.sel(time=t).mean('time').values,
                    v800.sel(time=t).mean('time').values,
                    transform=ccrs.PlateCarree(),
                    regrid_shape=8)
    
    ivtxx = ivtx.sel(time=t).mean('time')
    ivtyy = ivty.sel(time=t).mean('time')
    ivt = (ivtxx**2+ivtyy)**0.5
    mapa_ivt = ax[i,3].pcolormesh(lon2d,lat2d,
                                  ivt,vmin=vmin[i][3],vmax=vmax[i][3],
                                  cmap='twilight',
                                  rasterized=True,
                                  shading='auto')
    fig.colorbar(mapa_ivt,ax=ax[i,3], ticks=ticks[i][3])
    ax[i,3].quiver(lon2d,lat2d,
                    ivtxx.values,ivtyy.values,
                    transform=ccrs.PlateCarree(),
                    regrid_shape=8)
            
# print(len(ro
    for j in range(4):
        ax[i,j].scatter(-70.8206, -34.9968,  marker='*',edgecolor='k',linewidth=0.5,
                        c='red', s=150, zorder=150)

fig.text(0.06,0.75,'ROS',ha='center',va='center',rotation=90)
fig.text(0.06,0.25,'NO-ROS',ha='center',va='center',rotation=90)
ax[0,0].set_title('900hPa Temp. \nanomaly (°C)\n MSLP anomaly\n(contours)',fontsize=14,loc='right')
ax[0,1].set_title('300hPa Winds (m/s)\n Z500 anomaly\n(contours)',fontsize=14,loc='right')
ax[0,2].set_title('Precipitable water\nanomaly (mm)\n800hPa winds\n(vectors)',fontsize=14,loc='right')
ax[0,3].set_title('Integrated water\nvapor transport\n($kg\cdot m^{-1} s^{-1}$)',fontsize=14,loc='right')
# box = ax[0,0].get_position()
# cax = fig.add_axes([box.xmin,box.ymin-0.065,box.xmax-box.xmin,0.025])                  
# fig.colorbar(mapat,cax=cax,#label='900hPa Temp.\nanomaly (°C)\n MSLP anomaly\n(contours)',
#               ticks=[-4,-2,0,2,4], orientation='horizontal')

# box = ax[0,1].get_position()
# cax = fig.add_axes([box.xmin,box.ymin-0.065,box.xmax-box.xmin,0.025])                  
# fig.colorbar(mapawinds,cax=cax, orientation='horizontal',
#              ticks=[0,15,30,45])

# box = ax[0,2].get_position()
# cax = fig.add_axes([box.xmin,box.ymin-0.065,box.xmax-box.xmin,0.025])                  
# fig.colorbar(mapaw,cax=cax, orientation='horizontal',
#              ticks=[-6,-3,0,3,6])

# box = ax[0,3].get_position()
# cax = fig.add_axes([box.xmin,box.ymin-0.065,box.xmax-box.xmin,0.025])                  
# fig.colorbar(mapa_ivt,cax=cax, orientation='horizontal')
# # cax.set_ylim(-2,4)

# box = ax[1,-1].get_position()
# cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])                  
# fig.colorbar(mapawinds,cax=cax,label='300hPa Winds\n (m/s)\n Z500 anomaly\n(contours)')

# box = ax[2,-1].get_position()
# cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])
# fig.colorbar(mapaw,cax=cax, label='Precipitable\nwater anomaly\n(mm)\n800hPa winds\n(vectors)',
#              ticks=[-6,-3,0,3,6])

# box = ax[3,-1].get_position()
# cax = fig.add_axes([box.xmax*1.01,box.ymin,0.01,box.ymax-box.ymin])
# fig.colorbar(mapa_ivt,cax=cax, label='Integrated\nwater vapor\ntransport\n($kg m^{-1} s^{-1}$)',
#              ticks=[0,100,200,300])

plt.savefig('plots/COMPUESTO_ROS.pdf',dpi=150,bbox_inches='tight')