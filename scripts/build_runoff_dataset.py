#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:45:05 2022

@author: lucas

# =============================================================================
# This Script grabs runoff time series and build a single dataset for the main
# hydrological basins studied in this work.
# =============================================================================
"""

import pandas as pd
import geopandas as gpd
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

# %%
# =============================================================================
# LOAD CATCHMENT ATTRIBUTES
# =============================================================================

basin_attributes = pd.read_csv('datos/basins_attributes.csv',
                               index_col=0)

basin_attributes = basin_attributes.sort_values(by='gauge_lat')

blon, blat = basin_attributes.gauge_lon, basin_attributes.gauge_lat

# %%
# =============================================================================
# LOAD RUNOFF DATA
# =============================================================================
paths = glob('datos/estaciones/dga/qinst*')  # paths
basins = [p.split("_")[-1].split(".")[0] for p in paths]  # basin name
runoff_data = [pd.read_csv(p, index_col=0) for p in paths]  # load data
for d in runoff_data:
    d.index = pd.to_datetime(d.index)  # fix time index
# first timestamp
fd = np.min([d.index[0] for d in runoff_data])
# last timestamp
ld = np.max([d.index[-1] for d in runoff_data])
tot_interval = pd.date_range(fd, ld, freq='h')

reindex_data = [d.reindex(tot_interval) for d in runoff_data]
runoff_data = pd.concat(reindex_data, axis=1)
# drop data on dates where any gauge have data
runoff_data = runoff_data.dropna(how='all')
runoff_data.columns = basins

# %%
runoff_data.columns = ["Rio Achibueno En La Recova",
                       "Rio Aconcagua En Chacabuquito",
                       "Rio Ancoa En El Morro",
                       "Rio Cachapoal En Pte Termas De Cauquenes",
                       "Rio Claro En El Valle",
                       "Rio Claro En Hacienda Las Nieves",
                       "Rio Claro En Rauquen",
                       "Rio Colorado En Junta Con Palos",
                       "Rio Longavi En El Castillo",
                       "Rio Mapocho En Los Almendros",
                       "Rio Melado En El Salto",
                       "Rio Palos En Junta Con Colorado",
                       "Rio Perquilauquen En San Manuel",
                       "Rio Teno Despues De Junta Con Claro",
                       "Rio Uble En San Fabian N 2",
                       "Rio Tinguiririca Bajo Los Briones",
                       "Rio Maipo En El Manzano"]

runoff_data = runoff_data[basin_attributes.gauge_name]
# runoff_data.columns = basin_attributes.index
runoff_data.to_csv('datos/runoff_gauges_dataset.csv')

# #%%
# polygons = gpd.read_file('datos/vector/catchments_camels_cl_v1.3.shp',index_col=0)
# polygons.index = polygons.gauge_id
# polygons = polygons.loc[basin_attributes.index]

# #%%

# polygons.plot("gauge_name")
# polygons.boundary.plot()
