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

# %%
dem = "datos/topography/Chile_DEM_0001x0001grad.nc"
basin_polygons = gpd.read_file("datos/vector/cuencas_CAMELS.gpkg")

# %%

basins = [8106002, 7330001, 7350003, 7354002, 7355002, 7317005, 7379002,
          7112001, 7115001, 7104002, 6028001, 6027001, 6013001, 6008005,
          5710001, 5722002, 5410002]

# basins=[7303000]
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
    os.system("gdalwarp -cutline "+path+" -crop_to_cutline " +
              dem+" datos/topography/basins/"+name+".nc")
os.system("rm -rf tmp/*")
