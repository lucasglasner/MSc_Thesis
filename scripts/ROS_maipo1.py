#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:07:52 2022

@author: lucas


# =============================================================================
# ROS in Maipo Basin
# =============================================================================

"""

import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime as dt
import scipy.stats as st
import xarray as xr
from scipy.stats import linregress
from tqdm import trange
from glob import glob
import gc
import geopandas as gpd
import scipy.stats as st
import sys
# %%
sys.path.append('functions.py')
from functions import seasonal_decompose, local_minimum_filter
# %%
# =============================================================================
# load data from observations....
# =============================================================================

# hypsometry
hypso = pd.read_csv(
    'datos/topography/basins/hypso/RioMaipoEnElManzano_hypso.csv')

# pr cr2met
pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv", index_col=0)
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met.drop("date", axis=1, inplace=True)
pr_cr2met = pr_cr2met.squeeze()

pr_dgf = pd.read_csv(
    "datos/estaciones/dgf/DATOSUTC_2004-2019.csv", index_col=0)
pr_dgf.index = pd.to_datetime(
    pr_dgf.index.values)-datetime.timedelta(hours=4)
pr_dgf = pr_dgf.iloc[:, 9].resample('h').sum()

# pr la obra
pr_sanjose = pd.read_csv('datos/estaciones/pr_SanJoseMaipo.csv', dtype=str)
pr_sanjose.index = pd.to_datetime(
    pr_sanjose.iloc[:, 0]+"-"+pr_sanjose.iloc[:, 1]+"-"+pr_sanjose.iloc[:, 2])
pr_sanjose = pr_sanjose.iloc[:, 3]
pr_sanjose = pd.to_numeric(pr_sanjose).reindex(pr_cr2met.index)


# isoterma0
H0 = pd.read_csv(
    "datos/isotermas0_maipomanzano.csv", index_col=0).dropna(how="all")
H0.index = pd.to_datetime(H0.index)

H0 = H0['STODOMINGO'].dropna()
H0_h = H0.resample('h').interpolate('cubicspline')
H0_h = H0_h.where(H0_h>500).where(H0_h<6e3)
H0 = H0.resample("d").mean()
H0 = H0.reindex(pr_cr2met.index)


# area pluvial
int_func = interp1d(hypso.height, hypso.fArea)
pluv_area = int_func(H0-300)
pluv_area_h = int_func(H0_h-300)
pluv_area_h = pd.Series(pluv_area_h, index=H0_h.index)
pluv_area = pd.Series(pluv_area, index=H0.index)

# snow cover
SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0).dropna(how='all')
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA['IANIGLA'].reindex(pr_cr2met.index)/100


SCA_trend = SCA.shift(-3)-SCA
# snow limit

SL = pd.read_csv('datos/snowlimits_maipomanzano.csv',
                 index_col=0).dropna(how='all')
SL.index = pd.to_datetime(SL.index)
# SL = SL['MODIS_H50'].reindex(pr_cr2met.index)
for i in range(3):
    y = SL.iloc[:, i].dropna()
    x = SL["IANIGLA"].dropna().reindex(y.index)
    m = st.linregress(x, y)
    for j in range(len(SL)):
        if np.isnan(SL.iloc[j, i]):
            SL.iloc[j, i] = m.slope * \
                SL["IANIGLA"].values[j]+m.intercept

SL = SL['MODIS_H50'].reindex(pr_cr2met.index)
# SL_trend = pd.Series(np.empty(SL.shape[0]), index=SL.index)
SL_trend = SL.shift(-3)-SL

# area ros
ros_area = SCA-(1-pluv_area)
# ros_area = pd.Series(int_func(H0)-int_func(SL),SL.index)
ros_area = ros_area.where(pr_sanjose>3).where(ros_area>0)
# caudal
qinst_mm = pd.read_csv(
    "datos/estaciones/qinst_RioMaipoEnElManzano.csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)


quickflow = local_minimum_filter(qinst_mm, 25)[1]
baseflow  = qinst_mm-quickflow
qmax_d = qinst_mm.resample('d').max().reindex(pr_cr2met.index)
quickflow = quickflow.resample('d').mean().reindex(pr_cr2met.index)


# %%
# =============================================================================
# load reanalysis data
# =============================================================================
paths = "datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES_SWE*"
paths = glob(paths)
SWE = xr.open_mfdataset(paths, chunks='auto').SWE
# paths = 'datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES_dSWE_*.nc'
# dSWE = xr.open_mfdataset(paths).SWE
# dSWE = dSWE.reindex({'time': SWE.time.to_series().index})

dSWE = (SWE.shift({'time': -1})-SWE.shift({'time':1}))/2
paths = "datos/era5/H0_ERA5_*.nc"
H0_era5 = xr.open_mfdataset(paths, chunks='auto').deg0l
H0_era5 = H0_era5.reindex({'time': SWE.time.to_series().index},
                          method='nearest')


H0_era5_coast = H0_era5.sel(lon=-71.65, lat=-33.65, method='nearest')
H0_era5_coast = H0_era5_coast.to_series()

H0_era5 = H0_era5.reindex({"time": SWE.time.to_series().index,
                           'lat': SWE.lat,
                           'lon': SWE.lon},
                          method='nearest')

dem = 'datos/topography/basins/RioMaipoEnElManzano_CR2MET.nc'
dem = xr.open_dataset(dem).Band1

freeze = []
for t in H0_era5_coast.index:
    mask = xr.where(dem < H0_era5_coast.loc[t]-300, True, False)
    freeze.append(mask)

freeze = xr.concat(freeze, H0_era5_coast.index)

paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_pr_1979-2020.nc"
PR = xr.open_dataset(paths, chunks='auto').pr
PR = PR.reindex({"time": SWE.time.to_series().index},
                method='nearest')

# %%
# =============================================================================
# daily data in a single table
# =============================================================================
interval = slice(datetime.datetime(2000, 2, 25),
                 datetime.datetime(2020, 4, 30))
volume_h = (pluv_area_h.reindex(pr_dgf.index)*pr_dgf).resample(
    'd').sum()*4843.19/1e3
data_daily = pd.concat([pr_cr2met, pr_sanjose, SCA, pluv_area,
                        ros_area, H0, SL, SCA_trend, SL_trend, qmax_d,
                        quickflow, pr_dgf.resample('d').max(),
                        pr_dgf.resample('d').sum(), volume_h], axis=1)
data_daily.columns = ['pr_cr2met', 'pr_sanjose', 'sca', 'pluvarea',
                      'ros_area', 'fl', 'sl', 'sca_trend', 'sl_trend',
                      'qmax_d', 'quickflow','max_pr_int','pr_dgf',
                      'volume_h']
data_daily['date'] = pr_cr2met.index
data_daily['delta_ros'] = data_daily['fl']-data_daily['sl']
data_daily['new_events'] = data_daily.groupby(
    'pr_sanjose')['date'].apply(lambda s: s.diff().dt.days < 2)
data_daily['new_events'] = data_daily['new_events'].rolling(
    2, center=True).min().astype(bool)
data_daily['events'] = data_daily['new_events'].cumsum()
data_daily = data_daily.groupby([data_daily.index, data_daily.events]).mean()
# data_daily = data_daily[data_daily['pr_sanjose'] > 0]
# data_daily = data_daily[data_daily['ros_area']>0.1]
# data_daily = data_daily[interval]
# data_daily = data_daily["2004":]



# %%
# =============================================================================
# Group precipitation days by event
# =============================================================================
data_events = data_daily['pr_sanjose'].unstack().sum(axis=0)
data_events = pd.concat([data_events], axis=1)
data_events.columns = ['pr_sanjose']
data_events = data_events[data_events['pr_sanjose']>3]
data_events['max_pluv_area'] = data_daily.pluvarea.groupby(
    'events').apply(lambda x: x.max())
data_events['start'] = data_daily.sca.groupby(
    'events').apply(lambda x: x.index[0][0])
data_events['end'] = data_daily.sca.groupby(
    'events').apply(lambda x: x.index[-1][0])
data_events['start_sca'] = data_daily.sca.groupby(
    'events').apply(lambda x: x[0])
data_events['end_sca'] = data_daily.sca.groupby(
    'events').apply(lambda x: x[-1])
data_events['duration'] = data_events['end']-data_events['start']
# data_events['duration'] = data_events['duration'].map(lambda x: x.days+0)
data_events['delta_sca'] = data_events['end_sca']-data_events['start_sca']
# data_events['middate'] = data_events.start + \
#     data_events['duration'].map(lambda x: pd.Timedelta(days=x//2))

data_events['qmax'] = data_daily.qmax_d.groupby(
    'events').apply(lambda x: x.max())
data_events['return_period'] = data_events['qmax'].apply(
    lambda x: (1-st.percentileofscore(qinst_mm,x)/100))
data_events['max_ros_area'] = data_daily.ros_area.groupby(
    'events').apply(lambda x: x.max())
data_events['pr_cum'] = data_daily.pr_dgf.groupby(
    'events').apply(lambda x: x.sum())
data_events['pr_max'] = data_daily.max_pr_int.groupby(
    'events').apply(lambda x: x.max())
data_events['tot_vol_h'] = data_daily.volume_h.groupby(
    'events').apply(lambda x: x.sum())
data_events['tot_vol'] = data_daily.groupby(
    'events').apply(lambda x: 4843.19*np.sum(x['pr_sanjose']*x['pluvarea'])/1e3)
data_events['vol_dir'] = data_daily.groupby(
    'events').apply(lambda x: np.sum(x['quickflow']*3600*24)/1e6)
data_events['return_period'] = 1/data_events['return_period']
data_events.duration = data_events.duration.map(lambda x: x.days-2)
# data_events = data_events[data_events['return_period']>10]
data_events = data_events[data_events['tot_vol']>0]

#%%
plt.rc('font',size=18)
mask = (data_events['max_ros_area']>0.1) & (data_events['delta_sca']<0)
color = data_events[mask]['max_ros_area']
# mask = np.logical_and(mask,data_events['return_period'])
fig,ax = plt.subplots(1,1)
ax.set_xscale('log')
ax.set_yscale('log')
sc = ax.scatter(data_events[mask]['tot_vol'],data_events[mask]['vol_dir'],
                c = color, cmap='tab10',
                edgecolor='k',zorder=2,vmin=0,vmax=1)
ax.scatter(data_events[~mask]['tot_vol'], data_events[~mask]['vol_dir'],
            color='grey', alpha=0.3)
ax.plot([0,1e3],[0,1e3],color='k',ls=":")
ax.set_ylim(0.01,1000)
ax.set_xlim(0.01,1000)
ax.set_xlabel(r'$\int PR \cdot A_p dt$ $(hm^3)$')
ax.set_ylabel(r'$V_d$ $(hm^3)$')
fig.colorbar(sc, ticks=np.arange(0,1.1,0.1))

# ax.scatter(data_events.loc[6686]['tot_vol'],data_events.loc[6686]['vol_dir'],
#            color='purple',zorder=3)
# # plt.plot([0,.2],[0,.2],color='purple')
# plt.xlim(0,0.3)
# plt.ylim(0,0.3)
# plt.yscale('log')
# plt.xscale('log')
#%%
fig,ax = plt.subplots(1,2, figsize=(15,5))


dSWE = (SWE.shift({'time':-1})-SWE.shift({'time':1}))/2
import matplotlib as mpl
d = "2013-08-12"
dSWE.sel(time=d).plot(ax=ax[0], cmap='RdBu',
                      norm=mpl.colors.TwoSlopeNorm(0,-10,40))
PR.sel(time=d).plot(ax=ax[1])



#%%


#%%
# =============================================================================
# Ros on the reanalyisis mixture
# =============================================================================

th = 0.05
pr_th = 3
pr_th_cr2met = 10



ROS = xr.where((SWE > 100) & (H0_era5>300) & (PR > pr_th_cr2met)
                & ((dSWE/SWE < 0.05)),
               True, False).load()
ROS_CCE = np.empty(ROS.shape[0])
for i in range(ROS.shape[0]):
    ROS_CCE[i] = ROS[i, :, :].sum()

ROS_CCE = pd.Series(ROS_CCE, index=SWE.time.values)
ROS_CCE = ROS_CCE/191
# del SWE, H0_era5, PR, paths

ROS_CCE = ROS_CCE.reindex(ros_area.index)
ROS_CCE = ROS_CCE.where(ROS_CCE>0).where(pr_cr2met>pr_th_cr2met)
ROS_CCE = ROS_CCE["1985":"2014"]


#%%
#=============================================================================#
#
#=============================================================================
fig, ax = plt.subplots(2, 2, figsize=(12, 4), sharex='col')
# fig.tight_layout(pad=1.5)
plt.rc('font', size=18)
ax = ax.ravel()

common_period = set(ros_area.dropna().index)
common_period = list(common_period.intersection(ROS_CCE.dropna().index))
common_period = pd.Index(common_period)
common_period = slice(common_period.min(),common_period.max())
# common_period = slice('1900','2030')


ros_area1 = ros_area.where(ros_area>0).where(pr_sanjose>pr_th).dropna()
ros_area1 = (ros_area1.copy()>th).astype(float)

ROS_CCE1 = ROS_CCE.where(ROS_CCE>0).where(pr_cr2met>pr_th_cr2met).dropna()
ROS_CCE1 = (ROS_CCE1.copy()>th).astype(float)

ac1 = (ros_area1).groupby(ros_area1.index.month).sum()/(20)
ac1 = ac1.reindex(np.arange(1,13))

ac2 = (ROS_CCE1).groupby(ROS_CCE1.index.month).sum()/30
ac2 = ac2.reindex(np.arange(1,13))


# pr = PR.sel(lat=-33.64,lon=-70.35,method='nearest').to_series()
pr = pr_cr2met
sanjose_yrs = pr_sanjose.dropna().index[-1].year-pr_sanjose.dropna().index[0].year+1
rd1 = (pr_sanjose>pr_th).groupby(pr_sanjose.index.month).sum()/sanjose_yrs
rd2 = (pr>pr_th_cr2met).groupby(pr.index.month).sum()/(2015-1984+1)



# sc=ax[0].scatter(ros_area1,ROS_CCE1, edgecolor='k',
# #                  s=50, alpha=0.8)
# bias = np.sum(ros_area1-ROS_CCE1)/np.sum(ros_area1)
# pearson = pd.concat([ros_area1,ROS_CCE1],axis=1).dropna()
# pearson = st.pearsonr(pearson.iloc[:,0],pearson.iloc[:,1])
# ax[0].text(0.02,0.75,'pBias: '+'{:.2%}'.format(bias)+'\n'+r'$R^2: $'+'{:.2f}'.format(pearson[0]**2),
#             transform=ax[0].transAxes)


ax[3].bar(np.arange(1,13),rd1,width=0.25,align='edge',color='royalblue',
          edgecolor='k')

ax[3].bar(np.arange(1,13),rd2,width=-0.25,align='edge',color='purple',
          edgecolor='k')

ax[3].set_xticks(np.arange(1,13))
ax[3].set_xlim(0,13)
ax[3].set_xticklabels(['JAN','FEB','MAR','APR','MAY','JUN',
                        'JUL','AUG','SEP','OCT','NOV','DEC'])
ax[3].tick_params(axis="x",rotation=45)



ax[1].bar(np.arange(1,13),ac1.reindex(np.arange(1,13)),width=0.25,edgecolor='k',
          align='edge', color='royalblue')
ax[1].bar(np.arange(1,13),ac2,width=-0.25,edgecolor='k',
            align='edge', color='purple')
ax[1].set_xticks(np.arange(1,13))
ax[1].set_xlim(0,13)
ax[1].set_xticklabels(['JAN','FEB','MAR','APR','MAY','JUN',
                        'JUL','AUG','SEP','OCT','NOV','DEC'])
ax[1].tick_params(axis="x",rotation=45)

ycum1 = (ros_area.dropna()>th).groupby(ros_area.dropna().index.year).sum()
ycum2 = (ROS_CCE.dropna()>th).groupby(ROS_CCE.dropna().index.year).sum()


prycum1 = (pr_sanjose>pr_th).groupby(pr_sanjose.index.year).sum()
prycum2 = (pr>pr_th_cr2met).groupby(pr.index.year).sum()

ax[0].plot(ycum1, color='royalblue', label='Observations')
ax[0].plot(ycum2, color='purple', label='Reanalysis')
ax[0].set_xlim(1984,2019)
ax[0].set_xticks(np.arange(1984,2020,4))
ax[2].tick_params(axis="x", rotation=45)

ax[2].plot(prycum1, color='royalblue', label='Treshold: 3mm SanJose')
ax[2].plot(prycum2, color='purple', label='Treshold: 10mm CR2MET\nbasinwide mean.')
ax[2].legend(frameon=False,fontsize=12, loc=(-0.09,-1.1))


ax[1].text(0.02,0.85, r"$\langle N \rangle$: "+'{:.1f}'.format(ycum2.mean()),
           transform=ax[1].transAxes, color='purple')
ax[1].text(0.02,0.7, r"$\langle N \rangle$: "+'{:.1f}'.format(ycum1.mean()),
           transform=ax[1].transAxes, color='royalblue')
ax[0].legend(frameon=False, loc=(0.62,-2.2), fontsize=12, ncol=1)



ax[3].text(0.02,0.85, r"$\langle N \rangle$: "+'{:.1f}'.format(prycum2.mean()),
           transform=ax[3].transAxes, color='purple')
ax[3].text(0.02,0.7, r"$\langle N \rangle$: "+'{:.1f}'.format(prycum1.mean()),
           transform=ax[3].transAxes, color='royalblue')
# ax[3].yaxis.tick_right()
# ax[1].yaxis.tick_right()

ax[0].set_yticks([0,3,6,9,12,15,18])
ax[2].set_ylim(0,40)
ax[2].set_yticks([0,10,20,30,40,50])
ax[1].set_yticks([0,0.5,1,1.5])
ax[3].set_yticks([0,1,2,3,4,5,6])

ax[0].set_ylabel('#ROS Days')
ax[2].set_ylabel('#Rain Days')

ax[0].set_title('Yearly Trends', loc='left')
ax[1].set_title('Mean anual cycle', loc='left')

plt.savefig('plots/maipomanzano/ROS_maipo_3.pdf')