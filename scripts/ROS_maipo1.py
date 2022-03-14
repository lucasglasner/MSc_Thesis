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
from functions import seasonal_decompose, local_minimum_filter
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


# pr la obra
pr_laobra = pd.read_csv('datos/estaciones/pr_laobra.csv', dtype=str)
pr_laobra.index = pd.to_datetime(
    pr_laobra.iloc[:, 0]+"-"+pr_laobra.iloc[:, 1]+"-"+pr_laobra.iloc[:, 2])
pr_laobra = pr_laobra.iloc[:, 3]
pr_laobra = pd.to_numeric(pr_laobra).reindex(pr_cr2met.index)


# isoterma0
H0 = pd.read_csv(
    "datos/isotermas0_maipomanzano.csv", index_col=0).dropna(how="all")
H0.index = pd.to_datetime(H0.index)
H0 = H0.resample("d").mean()
H0 = H0.reindex(pr_cr2met.index)
H0 = H0['STODOMINGO']

# area pluvial
int_func = interp1d(hypso.height, hypso.fArea)
pluv_area = int_func(H0-300)
pluv_area = pd.Series(pluv_area, index=H0.index)
pluv_area = pluv_area.where(pr_laobra > 0)

# snow cover
SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0).dropna(how='all')
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA['IANIGLA'].reindex(pr_cr2met.index)/100

SCA_trend = pd.Series(np.empty(SCA.shape[0]), index=SCA.index)
for i in range(len(SCA_trend)):
    if i == len(SCA_trend)-1:
        SCA_trend[i] = SCA[i]-SCA[i-1]
    else:
        SCA_trend[i] = SCA[i+1]-SCA[i]

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
SL_trend = pd.Series(np.empty(SL.shape[0]), index=SL.index)
for i in range(len(SL_trend)):
    if i == len(SL_trend)-1:
        SL_trend[i] = SL[i]-SL[i-1]
    else:
        SL_trend[i] = SL[i+1]-SL[i]

# area ros
ros_area = SCA-(1-pluv_area)
ros_area = ros_area.where(ros_area > 0).fillna(0)

# caudal
qinst_mm = pd.read_csv(
    "datos/estaciones/qinst_RioMaipoEnElManzano.csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)
qinst_mm = qinst_mm.resample('d').max().reindex(pr_cr2met.index)

q_anomaly = local_minimum_filter(qinst_mm, 40)[1]

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

dSWE = SWE.shift({'time': -1})-SWE
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


ROS_CCE = xr.where((SWE > 10) & (freeze) & (PR > 3) & (dSWE <= 0),
                   True, False).load()
ROS = np.empty(ROS_CCE.shape[0])
for i in range(ROS_CCE.shape[0]):
    ROS[i] = ROS_CCE[i, :, :].sum()

ROS_CCE = pd.Series(ROS, index=SWE.time.values)
ROS_CCE = ROS_CCE/191
# del SWE, H0_era5, PR, paths

ROS_CCE = ROS_CCE.reindex(ros_area.index)

# %%
# =============================================================================
# ros daily data in a single table
# =============================================================================
interval = slice(datetime.datetime(2000, 2, 25),
                 datetime.datetime(2020, 4, 30))

data_daily = pd.concat([pr_cr2met, pr_laobra, SCA, pluv_area,
                        ros_area, H0, SL, SCA_trend, SL_trend, qinst_mm,
                        q_anomaly], axis=1)
data_daily.columns = ['pr_cr2met', 'pr_laobra', 'sca', 'pluvarea',
                      'ros_area', 'fl', 'sl', 'sca_trend', 'sl_trend',
                      'qmax_d', 'quickflow']
data_daily['date'] = pr_cr2met.index
data_daily['delta_ros'] = data_daily['fl']-data_daily['sl']
data_daily['new_events'] = data_daily.groupby(
    'pr_laobra')['date'].apply(lambda s: s.diff().dt.days < 2)
data_daily['new_events'] = data_daily['new_events'].rolling(
    2, center=True).min().astype(bool)
data_daily['events'] = data_daily['new_events'].cumsum()
data_daily = data_daily.groupby([data_daily.index, data_daily.events]).mean()
data_daily = data_daily[data_daily['pr_laobra'] > 3]
data_daily = data_daily[interval]


# %%
# =============================================================================
# Group precipitation days by event
# =============================================================================
data_events = data_daily['pr_laobra'].unstack().sum(axis=0)
data_events = pd.concat([data_events], axis=1)
data_events.columns = ['pr_laobra']
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
data_events['duration'] = data_events['duration'].map(lambda x: x.days+1)
data_events['delta_sca'] = data_events['end_sca']-data_events['start_sca']
data_events['middate'] = data_events.start + \
    data_events['duration'].map(lambda x: pd.Timedelta(days=x//2))


data_events['max_ros_area'] = data_daily.ros_area.groupby(
    'events').apply(lambda x: x.max())


# %%
# =============================================================================
# compute some ploteable stats
# =============================================================================

# ROS = 

# %%
# fig, ax = plt.subplots(2, 2, figsize=(14, 7))
# fig.tight_layout(pad=1.5)
# plt.rc('font', size=18)
# ax = ax.ravel()
# # month = ROS.resample("m").sum().applymap(lambda x: np.nan if x == 0 else x)
# for pair in pairs:
#     var = year[pair].dropna()
#     m = st.linregress(range(len(var)), var)
#     trend = "{:.2f}".format(m.slope)
#     pvalue = "{:.2f}".format(m.pvalue)
#     ax[1].plot(var, label="Trend: "+trend+" days/year",
#                alpha=0.8)
# ax[1].plot(year[pairs[0]].dropna(), color='tab:blue', lw=2)
# ax[1].legend(loc=(0, 1), frameon=False, fontsize=16)
# # ax[1].set_xticklabels()
# # ax[1].sharex(ax[3])
# ax[1].tick_params(axis='x', rotation=45)
# meanROS.plot.bar(ax=ax[0], alpha=0.8, edgecolor='k')


# for i in meanROS.columns:
#     ax[0].plot(np.arange(0, 12, 1), meanROS[i], ls=":")
# #     var = ROS.groupby(ROS.index.dayofyear).sum()[pair]
# #     var = seasonal_decompose(var,365,6,1)[0]
# #     ax[0].plot(np.linspace(1,12,len(var)),var)
# ax[0].legend(frameon=False, loc=(0, 1), fontsize=16)
# ax[0].set_ylim(0, 7)
# ax[0].set_xticks(np.arange(0, 12, 1))
# ax[0].set_xticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
#                        "AGO", "SEP", "OCT", "NOV", "DIC"])
# ax[0].tick_params(axis='x', rotation=45)
# ax[0].grid(ls=":", axis="x")
# ax[0].set_ylabel("Monthly mean\nROS Events)")
# ax[1].set_ylabel("N° ROS days")

# ax[2].scatter(ROS1.iloc[:, 0], ROS1.iloc[:, 1], color='tab:orange', ec='k',
#               alpha=0.5)
# ax[2].scatter(ROS1.iloc[:, 0], ROS1.iloc[:, 2], color='tab:green', ec='k',
#               alpha=0.5)
# ax[2].set_ylim(0, 1)
# ax[2].set_xlim(0, 1)
# ax[2].plot([0, 1], [0, 1], 'k')
# ax[2].set_ylabel('Model ROS (REANALYSIS)\nFractional Area (-)')
# ax[2].set_xlabel('Observed ROS (MODIS/STODOMINGO)\n Fractional Area (-)')


# durations.plot.bar(ax=ax[3], legend=False, ec='k')
# ax[3].set_xlim(0, 15)
# ax[3].set_xlabel('ROS Events Duration (days)')
# ax[3].set_xticks(np.arange(1, 15))
# ax[3].tick_params(axis='x', rotation=0)
# ax[3].set_ylabel('N°Events\nover total Events')
# ax[3].grid(axis='x')
# ax[3].set_ylim(0, 1)
