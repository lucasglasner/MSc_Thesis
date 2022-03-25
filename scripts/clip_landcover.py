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

# %%

landcover_path = "datos/landcover/LC_CHILE_2014_b_final.nc"
basin_polygons = gpd.read_file("datos/vector/cuencas_CAMELS.gpkg")


# %%
basins = [5410002, 7355002, 6008005, 6027001, 6013001, 5722002, 7317005,
          7330001, 8106002, 6028001, 7354002, 5710001, 7350003, 7104002,
          7115001, 7112001, 7116001, 7374001, 7372001]
# basins = ['Rio Choapa En Salamanca']

basins=[7321002]
# %%

polygons = [basin_polygons[basin_polygons["gauge_id"] == b] for b in basins]

# %%

for polygon in polygons:
    name = polygon.gauge_name.item().replace(" ", "")
    print(name)
    path = "tmp/"+name+".shp"
    polygon.to_file(path)
    if ~os.path.isfile("datos/vector/"+name+".shp"):
        polygon.to_file("datos/vector/"+name+".shp")
    os.system("gdalwarp -of netCDF -cutline "+path+" -crop_to_cutline " +
              landcover_path+" tmp/"+name+"_tmp.nc")
    os.system("gdalwarp -t_srs '+proj=longlat +datum=WGS84 +nodefs'" +
              " tmp/"+name+"_tmp.nc " +
              " datos/landcover/basins/"+name+".nc")
os.system("rm -rf tmp/*")
