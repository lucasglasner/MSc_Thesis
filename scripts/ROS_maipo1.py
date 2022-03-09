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
from functions import seasonal_decompose, sliding_interval_filter
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
isotermas0 = pd.read_csv(
    "datos/isotermas0_maipomanzano.csv", index_col=0).dropna(how="all")
isotermas0.index = pd.to_datetime(isotermas0.index)
isotermas0 = isotermas0.resample("d").fillna(method="ffill")
isotermas0 = isotermas0.reindex(pr_cr2met.index)

# area pluvial
int_func = interp1d(hypso.height, hypso.fArea)
pluv_area = int_func(isotermas0['STODOMINGO']-300)
pluv_area = pd.Series(pluv_area, index=isotermas0.index)

# snow cover
SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0).dropna(how='all')
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA['CORTES_10SWE'].reindex(pr_cr2met.index)/100

# pr la obra
pr_laobra = pd.read_csv('datos/estaciones/pr_laobra.csv', dtype=str)
pr_laobra.index = pd.to_datetime(
    pr_laobra.iloc[:, 0]+"-"+pr_laobra.iloc[:, 1]+"-"+pr_laobra.iloc[:, 2])
pr_laobra = pr_laobra.iloc[:, 3]
pr_laobra = pd.to_numeric(pr_laobra).reindex(pr_cr2met.index)

ros_area = SCA-(1-pluv_area)

# caudal
qinst_mm = pd.read_csv(
    "datos/estaciones/qinst_RioMaipoEnElManzano.csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)
qinst_mm = qinst_mm.resample('d').max().reindex(pr_cr2met.index)

# data days
data_daily = pd.concat([pr_cr2met, pr_laobra, SCA, pluv_area], axis=1)
data_daily.columns = ['pr_cr2met', 'pr_laobra', 'sca', 'pluvarea']
data_daily['date'] = data_daily.index
data_daily['ros_area'] = ros_area
data_daily['new_events'] = data_daily.groupby(
    'pr_laobra')['date'].apply(lambda s: s.diff().dt.days < 2)
data_daily['new_events'] = data_daily['new_events'].rolling(
    2, center=True).min().astype(bool)
data_daily['events'] = data_daily['new_events'].cumsum()
data_daily = data_daily.groupby([data_daily.index, data_daily.events]).mean()
data_daily = data_daily[data_daily['pr_laobra'] > 3]
data_daily = data_daily[data_daily['pluvarea'] > 0.3]

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
data_events = data_events.where(data_events['delta_sca'] < 0).dropna()
data_events = data_events.where(data_events['max_ros_area'] > 0.3).dropna()
data_events = data_events.where(data_events['pr_laobra'] > 10).dropna()
data_events = data_events.sort_values(by='delta_sca')

# data_events['sca_max'] = data_daily.sca.unstack().max(axis=0)
# data_events['start'] =
# data_events = data_events.where(data_events['pr_laobra'] > 10).dropna()


# %%
