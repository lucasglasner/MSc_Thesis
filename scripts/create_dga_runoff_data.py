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
import pandas as pd
import sys
sys.path.append("functions.py")

# %%
basins = ['AconcaguaEnChacabuquito']


# corregir_qintdga('datos/estaciones/dga/'+basins[0]+"/"+"1984.xlsx","1984")
for basin in basins:
    data = []
    for yr in range(1984, 2016):
        sheet = corregir_qintdga('datos/estaciones/dga/' +
                                 basin+'/'+str(yr)+'.xlsx', str(yr))
        data.append(sheet)

data = pd.concat(data, axis=0)
