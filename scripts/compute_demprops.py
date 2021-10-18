#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 11:16:12 2021

@author: lucas
python 
# =============================================================================
# Compute orentation pixel raster (aspect) from basins dem
# =============================================================================

"""

import os

#%%
inputfiles_path   = "datos/topography/basins/"
outputfiles_path  = "datos/topography/basins/"
basins            =  ["Rio Aconcagua En Chacabuquito", 
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
rasters=["aspect","accumulation","flowdirection"]
#%%
# os.system("rm -rf "+outputfiles_path+"*.nc")
for basin in basins:
    try:
        basin = basin.replace(" ","")
        if ~os.path.isfile(outputfiles_path+"aspect/"+basin+".nc"):
            program="gdaldem aspect "
            os.system(program+
                      inputfiles_path+basin+".nc "+
                      outputfiles_path+"aspect/"+basin+".nc")
        if ~os.path.isfile(outputfiles_path+"slope/"+basin+".nc"):
            program="gdaldem slope "
            os.system(program+
                      inputfiles_path+basin+".nc "+
                      outputfiles_path+"slope/"+basin+".nc -p")
    except:
        pass