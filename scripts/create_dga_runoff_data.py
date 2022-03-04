#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:23:48 2022

@author: lucas

# =============================================================================
#
# =============================================================================
"""


from functions import corregir_qintdga
import os
import pandas as pd
import sys
from glob import glob
sys.path.append("functions.py")

# %%
# basins = glob('datos/estaciones/dga/*[!.csv]')
# basins = [b.split("/")[-1] for b in basins]
basins = ["RioUpeoEnUpeo", "RioClaroEnCamarico", "RioLircayEnPuenteLasRastras"]


# corregir_qintdga('datos/estaciones/dga/'+basins[0]+"/"+"1984.xlsx","1984")
for basin in basins:
    data = []
    for yr in range(1984, 2016):
        try:
            sheet = corregir_qintdga('datos/estaciones/dga/' +
                                     basin+'/'+str(yr)+'.xlsx', str(yr))
        except:
            sheet = pd.Series([], dtype=float)
        data.append(sheet)

    data = pd.concat(data, axis=0)
    # data = data.reindex(pd.date_range('1984-01-01','2015-12-31', freq='h'))
    data.to_csv('datos/estaciones/dga/qinst_'+basin+'.csv')

# %%
