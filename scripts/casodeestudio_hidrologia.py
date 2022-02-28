#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:27:42 2022

@author: lucas

# =============================================================================
# This Script does some hydrological analysis in the case of study
# =============================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import datetime
import scipy.stats as st

# %%

interval = slice(datetime.datetime(2013, 8, 3),
                 datetime.datetime(2013, 8, 17))


# =============================================================================
# SL_mm
# =============================================================================

SL_mm = pd.read_csv('datos/snowlimits_maipomanzano.csv',
                    index_col=0).squeeze()
SL_mm.index = pd.to_datetime(SL_mm.index)+datetime.timedelta(hours=12)
SL_mm = SL_mm.dropna(how='all')

times = pd.date_range("2000-01-01", "2021-12-31", freq="d")
for i in range(3):
    y = SL_mm.iloc[:, i].dropna()
    x = SL_mm["IANIGLA"].dropna().reindex(y.index)
    m = st.linregress(x, y)
    for j in range(len(SL_mm)):
        if np.isnan(SL_mm.iloc[j, i]):
            SL_mm.iloc[j, i] = m.slope * \
                SL_mm["IANIGLA"].values[j]+m.intercept

SL_mm = SL_mm["MODIS_H50"]

SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0)['IANIGLA']
SCA.index = pd.to_datetime(SCA.index)

# =============================================================================
# ISOTERMAS 0
# =============================================================================

H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)-datetime.timedelta(hours=4)
# H0_mm = H0_mm.dropna(how='all')
# =============================================================================
# Q y pr maipo manzano
# =============================================================================

qinst_mm = pd.read_csv("datos/estaciones/qinst_RioMaipoEnElManzano.csv",
                       index_col=0)
qinst_mm.index = pd.to_datetime(qinst_mm.index)
qinst_mm = qinst_mm.squeeze()

pr_mm = pd.read_csv('datos/estaciones/pr_RioMaipoEnElManzano_2013-08.csv')
pr_mm.index = pd.to_datetime(pr_mm['Fecha'])
pr_mm = pr_mm['Valor'].drop_duplicates()
pr_mm = pr_mm.reindex(pd.date_range(interval.start,
                                    interval.stop, freq='h')).fillna(0)
# =============================================================================
# Estacion dgf
# =============================================================================
datos_dgf = pd.read_csv(
    "datos/estaciones/dgf/DATOSUTC_2004-2019.csv", index_col=0)
datos_dgf.index = pd.to_datetime(
    datos_dgf.index.values)-datetime.timedelta(hours=4)

# %%
datos = []
for dat in [SL_mm, SCA, H0_mm,
            qinst_mm.resample('d').max(),
            datos_dgf.iloc[:, 5].resample('d').mean(),
            datos_dgf.iloc[:, 5].resample('d').max(),
            datos_dgf.iloc[:, 5].resample('d').min(),
            datos_dgf.iloc[:, 9].resample('d').sum(),
            pr_mm.resample('d').sum()]:
    datos.append(dat[interval])

datos = pd.concat(datos, axis=1)
datos = datos.resample('d').mean().iloc[:-1, :]
datos.columns = ["SL", "SCA", "H0", "Qmax",
                 "Tprom_DGF", "Tmax_DGF", "Tmin_DGF", "PR_DGF",
                 "PR_MM"]
datos = np.round(datos, 1)
