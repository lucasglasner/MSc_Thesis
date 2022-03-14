#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: lucas

# =============================================================================
# download vismet data
# =============================================================================
"""
# %%
import datetime
import pandas as pd
import os
# https://vismet.cr2.cl/api/measure/by-station/242/by-measure-type/1/by-timestamp/1646967600/by-interval/72/csv
# https://vismet.cr2.cl/api/measure/by-station/182/by-measure-type/1/by-timestamp/1646967600/by-interval/72/csv
#https://vismet.cr2.cl/api/measure/by-station/143/by-measure-type/1/by-timestamp/1646967600/by-interval/72/csv
base_url = "https://vismet.cr2.cl/api/measure/by-station/"
dates = pd.date_range("2000-01-21", "2021-01-01",
                      freq="7d")
station_code = 242

interval = 168

url = base_url + str(station_code)+"/"
url = url + "by-measure-type/1/by-timestamp/"

output = "station_"+str(station_code)+"_"
output = output+dates[0].strftime("%Y-%m-%d")+"_"
output = output+dates[-1].strftime("%Y-%m-%d")+".csv"
os.system('rm -rf '+output)
os.system('touch '+output)
os.system('echo SID,Nombre,Código,Organización,Latitud,Longitud,Altura,Medida,Valor,Fecha,Confianza >> '+output)
for date in dates:
    os.system("sleep 1")
    timestamp = str(int(str(date.value)[:10])+14400)
    url = url + str(timestamp)+"/"
    url = url + "by-interval/"+str(interval)+"/csv"

    print(url)
    print(date.strftime("%Y-%m-%d"))
    os.system("curl "+url+' | tail -n +2>> '+output)
    url = base_url + str(station_code)+"/"
    url = url + "by-measure-type/1/by-timestamp/"