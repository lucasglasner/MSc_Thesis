#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:04:52 2021

@author: lucas
"""

from datetime import datetime
from siphon.simplewebservice.wyoming import WyomingUpperAir
import pandas as pd
import os
import numpy as np
from urllib.error import HTTPError
#%%
# =============================================================================
# Create datetime object for the sounding and string of the station identifier 
# Antofagasta  -> "SCFA"
# Sto. Domingo -> "SCSN"
# Pto. Montt   -> "SCTE"
# =============================================================================

dates   = pd.date_range("1980-01-01","1999-12-24",freq="12h")
station = "SCSN"

#%%
# =============================================================================
# Request data 
# =============================================================================

soundings = {date:np.nan for date in dates}
for i in range(len(dates)):
    date = dates[i]
    print(date)
    try: 
        name = "datos/stodomingo/radiosonde_"+datetime.strftime(date,"%Y-%m-%d_%H:%M:%S"+".csv")
        soundings[date] = WyomingUpperAir.request_data(date, station)
    except:
        print("pass")
        
clean_dict = {k: soundings[k] for k in soundings if not type(soundings[k]) == float}

np.save("datos/stodomingo/radiosondas.npy",clean_dict)
