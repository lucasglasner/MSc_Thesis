#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 10:45:29 2022

@author: lucas

# =============================================================================

# =============================================================================
"""

import pandas as pd
from glob import glob
import os
# %%
basins = ["Rio Aconcagua En Chacabuquito",
          "Rio Maipo En El Manzano",
          "Rio Colorado En Junta Con Palos",
          "Rio Putaendo En Resguardo Los Patos",
          "Rio Teno Despues De Junta Con Claro",
          "Rio Tinguiririca Bajo Los Briones",
          "Rio ÑUble En San Fabian"]

paths = sorted(glob('datos/ROS/%%TYPE%%*'))
data = [pd.read_csv(p, index_col=[0, 1]) for p in paths]

# %%
data_basins = []
for b in basins:
    data_basin = [d.loc[b] for d in data]
    data_basin = pd.concat(data_basin, axis=0)
    data_basin.index = pd.to_datetime(data_basin.index)
    data_basins.append(data_basin)

data_basins = pd.concat(data_basins, keys=basins)

data_basins.to_csv('datos/ROS/%%TYPE%%_Timeseries.csv')
os.system('rm -rf datos/ROS/%%TYPE%%_Timeseries_*.csv')
