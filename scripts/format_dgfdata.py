#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 15:26:38 2021

@author: lucas
@contact: lgvivanco96@gmail.com

# =============================================================================
# Script for formatting DGF station raw data. Grab each year fix julian date
# and hour to real UTC date, then cat everything and save a single .csv file.
# =============================================================================

"""

import os
import pandas as pd
import datetime

#%%
#Open data
data = {key:None for key in list(map(lambda x: str(x),range(2004,2020)))}
for yr in range(2004,2020):
    data_path = "datos/estaciones/dgf/DATOS_"+str(yr)+".DAT"
    data[str(yr)] = pd.read_csv(data_path,header=None)
    print(str(yr)+":   "+str(data[str(yr)].shape))
#%%
#paste date
for yr in range(2004,2020):
    print(yr)
    year = datetime.datetime(yr,1,1)
    dates = []
    for i in range(len(data[str(yr)])):
        dataDGF=data[str(yr)]
        jd     = dataDGF.iloc[i,1].item()-1
        string = str(dataDGF.iloc[i,2].item())
        minute = string[-2:]
        hours  = string[:-2]
        if len(minute) == 0:
            minute = 0
        else:
            minute = int(minute)
        if len(hours) == 0:
            hours = 0
        else:
            hours = int(hours)
        dates.append(year+datetime.timedelta(days=jd,hours=hours+4,minutes=minute))   
    dataDGF.index = dates
    data[str(yr)] = dataDGF
#%%
#concatenate years fix index and save
data = pd.concat(data.values())
data = data.reindex(pd.date_range(str(data.index[0]),str(data.index[-1]),freq="15min"))
data.to_csv("datos/estaciones/dgf/DATOSUTC_2004-2019.csv") 