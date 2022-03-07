#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: lucas

# =============================================================================
# download vismet data
# =============================================================================
"""
# %%
# https://vismet.cr2.cl/api/measure/by-station/132/by-measure-type/1/by-timestamp/1357056000/by-interval/72/csv
import datetime
import pandas as pd
import os
url = "https://vismet.cr2.cl/api/measure/by-station/"

dates = pd.date_range("2013-04-01", "2014-01-01",
                      freq="w")
station_code = 132
# timestamp = pd.to_datetime("2013-01-01T03:00:00").value
# timestamp = str(timestamp)[:10]
interval = 168

url = url + str(station_code)+"/"
url = url + "by-measure-type/1/by-timestamp/"

output = "station_"+str(station_code)+"_"
output = output+dates[0].strftime("%Y-%m-%d")+"_"
output = output+dates[-1].strftime("%Y-%m-%d")+".csv"
os.system('touch '+output)
# os.system('echo SID,Nombre,Código,Organización,Latitud,Longitud,Altura,Medida,Valor,Fecha,Confianza >> '+output)
for date in dates:
    timestamp = str(date.value)[:10]
    url = url + str(timestamp)+"/"
    url = url + "by-interval/"+str(interval)+"/csv"

    print(url)
    os.system("curl "+url+' | tail -n +2>> '+output)
