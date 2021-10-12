#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 18:26:09 2021

@author: lucas

# =============================================================================
# Clip DEM raster to basin polygons.
# =============================================================================

"""

import geopandas as gpd
import os

#%%
dem            = "datos/topography/Chile_DEM_0001x0001grad.nc"
basin_polygons = gpd.read_file("datos/vector/cuencas_CAMELS.gpkg")

#%%
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

polygons = [basin_polygons[basin_polygons["gauge_name"] == b] for b in basins]

#%%

for polygon in polygons:
    name = polygon.gauge_name.item().replace(" ","")
    print(name)
    path = "tmp/"+name+".shp"
    polygon.to_file(path)
    if ~os.path.isfile("datos/vector/"+name+".shp"):
        polygon.to_file("datos/vector/"+name+".shp")
    os.system("gdalwarp -cutline "+path+" -crop_to_cutline "+dem+" datos/topography/basins/"+name+".nc")
os.system("rm -rf tmp/*")
