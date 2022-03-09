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
# load data
# =============================================================================

# hypsometry
hypso = pd.read_csv(
    'datos/topography/basins/hypso/RioMaipoEnElManzano_hypso.csv')

# pr cr2met
pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv", index_col=0)
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met.drop("date", axis=1, inplace=True)
pr_cr2met = pr_cr2met.squeeze()

# isoterma0
H0 = pd.read_csv(
    "datos/isotermas0_maipomanzano.csv", index_col=0).dropna(how="all")
H0.index = pd.to_datetime(H0.index)
H0 = H0.resample("d").fillna(method="ffill")
H0 = H0.reindex(pr_cr2met.index)
H0 = H0['STODOMINGO']

# area pluvial
int_func = interp1d(hypso.height, hypso.fArea)
pluv_area = int_func(H0-300)
pluv_area = pd.Series(pluv_area, index=H0.index)

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

# pr la obra
pr_laobra = pd.read_csv('datos/estaciones/pr_laobra.csv', dtype=str)
pr_laobra.index = pd.to_datetime(
    pr_laobra.iloc[:, 0]+"-"+pr_laobra.iloc[:, 1]+"-"+pr_laobra.iloc[:, 2])
pr_laobra = pr_laobra.iloc[:, 3]
pr_laobra = pd.to_numeric(pr_laobra).reindex(pr_cr2met.index)

# area ros
ros_area = SCA-(1-pluv_area)

# caudal
qinst_mm = pd.read_csv(
    "datos/estaciones/qinst_RioMaipoEnElManzano.csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)
qinst_mm = qinst_mm.resample('d').max().reindex(pr_cr2met.index)

q_anomaly = local_minimum_filter(qinst_mm, 40)[1]
q_anomaly.plot()
# %%
# =============================================================================
# ros daily data in a single table
# =============================================================================
data_daily = pd.concat([pr_cr2met, pr_laobra, SCA, pluv_area,
                        ros_area, H0, SL], axis=1)
data_daily.columns = ['pr_cr2met', 'pr_laobra', 'sca', 'pluvarea',
                      'ros_area', 'fl', 'sl']
data_daily['date'] = pr_cr2met.index
data_daily['delta_ros'] = data_daily['fl']-data_daily['sl']
data_daily['new_events'] = data_daily.groupby(
    'pr_laobra')['date'].apply(lambda s: s.diff().dt.days < 2)
data_daily['new_events'] = data_daily['new_events'].rolling(
    2, center=True).min().astype(bool)
data_daily['events'] = data_daily['new_events'].cumsum()
data_daily['sca_trend'] = SCA_trend
data_daily['sl_trend'] = SL_trend
data_daily['qmax_d'] = qinst_mm
data_daily['qanomaly'] = q_anomaly
data_daily = data_daily.groupby([data_daily.index, data_daily.events]).mean()
data_daily = data_daily[data_daily['pr_laobra'] > 3]
# data_daily = data_daily[data_daily['ros_area'] > 0.1]
# data_daily = data_daily[data_daily['sca_trend'] <= 0.1]
# data_daily = data_daily[data_daily['sl_trend'] >= 0]


# %%
data_events = data_daily['pr_laobra'].unstack().sum(axis=0)
data_events = pd.concat(
    [data_events], axis=1)
data_events.columns = ['pr_laobra']
data_events['max_pluv_area'] = data_daily.pluvarea.groupby(
    'events').apply(lambda x: x.mean())
data_events['max_ros_area'] = data_daily.ros_area.groupby(
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
# data_events = data_events.where(data_events['delta_sca'] < 0).dropna()
# data_events = data_events.where(data_events['max_ros_area'] > 0).dropna()
data_events = data_events.where(data_events['pr_laobra'] > 10).dropna()
# data_events = data_events.sort_values(by='sl_trend')

# data_events['sca_max'] = data_daily.sca.unstack().max(axis=0)
# data_events['start'] =
# data_events = data_events.where(data_events['pr_laobra'] > 10).dropna()


# %%
