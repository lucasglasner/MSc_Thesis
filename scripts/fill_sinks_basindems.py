#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 11:16:12 2021

@author: lucas

# =============================================================================
# Fix dem topography filling sinks so the drainage network makes sense.
# Uses saga Wang & Liu algorithm and return a filled dem
# =============================================================================

"""

import os
import gdal

#%%
inputfiles_path   = "datos/topography/basins/"
os.system("mkdir "+inputfiles_path+"filled_sinks")
outputfiles_path  = "datos/topography/basins/filled_sinks/"
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
#%%
for basin in basins:
    try:
        basin = basin.replace(" ","")
        program     = "saga_cmd ta_preprocessor 4 "
        #Create filled dems
        order = program+"-ELEV "+inputfiles_path+basin+".nc"+" -FILLED "+outputfiles_path+basin
        os.system(order)
        #Save output as netcdf
        dem = gdal.Open(outputfiles_path+basin+".sdat")
        gdal.Translate(outputfiles_path+basin+".nc",dem)
        os.system("rm -rf "+outputfiles_path+"*[!.nc]")
        os.system("rm -rf "+inputfiles_path+basin+".nc")
        os.system("mv "+outputfiles_path+basin+".nc "+inputfiles_path+basin+".nc")
    except:
        pass
        os.system("rm -rf "+outputfiles_path+"*[!.nc]")
        os.system("rm -rf "+inputfiles_path+basin+".nc")
        os.system("mv "+outputfiles_path+basin+".nc "+inputfiles_path+basin+".nc")
    
os.system("rm -rf "+outputfiles_path)
