#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:40:32 2022

@author: lucas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
#%%

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