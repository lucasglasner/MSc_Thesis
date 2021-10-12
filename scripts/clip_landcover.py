#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 21:54:21 2021

@author: lucas

# =============================================================================
# Clip LandCover raster to basin polygons.
# =============================================================================

"""

import geopandas as gpd
import os

#%%

landcover_path = "datos/landcover/LC_CHILE_2014_b.tif"
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
    os.system("gdalwarp -of netCDF -cutline "+path+" -crop_to_cutline "+landcover_path+
              " tmp/"+name+"_tmp.nc")
    os.system("gdalwarp -t_srs '+proj=longlat +datum=WGS84 +nodefs'"+
              " tmp/"+name+"_tmp.nc "+
              " datos/landcover/basins/"+name+".nc")
os.system("rm -rf tmp/*")
