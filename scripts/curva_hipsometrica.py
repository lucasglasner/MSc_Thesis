#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:02:18 2021

@author: lucas

# =============================================================================
# Compute basin hypsometry from clipped DEM data.
# =============================================================================

"""

import xarray as xr
import numpy as np
import pandas as pd


basins = ["Rio Aconcagua En Chacabuquito", 
          "Rio Choapa En Lamahuida",
          "Rio Elqui En Algarrobal",
          "Rio Illapel En Huintil",
          "Rio Grande En Puntilla San Juan",
          "Rio Hurtado En Angostura De Pangue",
          "Rio Putaendo En Resguardo Los Patos",
          "Rio Mapocho En Los Almendros",
          "Rio Maipo En El Manzano",
          "Rio Cachapoal En Pte Termas De Cauquenes",
          "Rio Tinguiririca Bajo Los Briones",
          "Rio Teno Despues De Junta Con Claro",
          "Rio Colorado En Junta Con Palos",
          "Rio Maule En Armerillo",
          "Rio ÑUble En San Fabian"]

#%%
# =============================================================================
# Inputs
# =============================================================================
for name in basins:
    name = name.replace(" ","")
    print(name)
    path_basinDEM = "datos/topography/basins/"+name+".nc"
    basin         = xr.open_dataset(path_basinDEM, chunks="auto").Band1
    
    #Compute Pixel Area
    earth_radius  = 6371
    lon2d,lat2d   = np.meshgrid(basin.lon,basin.lat)
    dlon, dlat    = np.deg2rad(np.diff(basin.lon)[0]), np.deg2rad(np.diff(basin.lat)[0])
    areas         = earth_radius**2*np.cos(np.deg2rad(lat2d))*dlon*dlat


# =============================================================================
# Compute hypsometry
# =============================================================================
    
    heights    = np.linspace(0,7e3,500)
    basin_area = np.empty(heights.shape)
    for i in range(len(heights)):
        basin_area[i] = areas[basin.values < heights[i]].sum()
        
    fraction_area = basin_area/basin_area.max()

# =============================================================================
# Save hypsometry data
# =============================================================================

    data = pd.DataFrame([heights,fraction_area, basin_area],index=["height","fArea", "Area_km2"]).T
    np.round(data,3).to_csv(path_basinDEM[:-3].replace(name,"")+"hipso/"+name+
                            "_Hipso.csv", index=False)

