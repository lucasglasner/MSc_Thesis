#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 10:08:01 2021

@author: lucas

# =============================================================================
# 
# =============================================================================


"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import geopandas as gpd
from glob import glob
from tqdm import trange
from scipy.interpolate import interp1d

# %%

# =============================================================================
# Basin hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso = pd.read_csv(
    "datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2", inplace=True)
basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")

dem = xr.open_dataset("datos/topography/basins/" +
                      cuenca+"_Cortes.nc").Band1

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv")
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met = pr_cr2met["5710001"]

# %%
# =============================================================================
# Cortes reanalisis SWE
# =============================================================================

try:
    SWE_bands = pd.read_csv(
        "datos/ANDES_SWE_Cortes/snowbands_maipomadsnzano.csv", index_col=0)
    SWE_bands.columns = pd.to_datetime(SWE_bands.columns)
except:
    for yr in range(1986, 2015+1):
        # yr = 1985
        SWE = xr.open_dataset(
            "datos/ANDES_SWE_Cortes/ANDES_SWE_WY"+str(yr)+".nc").SWE

        # =============================================================================
        # Binarize Cortes SWE based on threshold
        # =============================================================================
        SWE_treshold = 10
        SWE = SWE > SWE_treshold

        # =============================================================================
        # Create elevation bands and masks
        # =============================================================================
        dz = 50
        elevation_bands = np.arange(dem.min(), dem.max()+dz, dz)

        masks = []
        for j in range(len(elevation_bands)-1):
            z0 = elevation_bands[j]
            z1 = elevation_bands[j+1]
            mask = (dem.where((dem > z0) & (dem < z1)) > 0).values
            masks.append(mask)

        elevation_bands = elevation_bands[:-1]
        # =========================================================================
        # Apply band mask, and compute % of snow cover distribution by band
        # =========================================================================

        SWE_bands = np.empty((len(elevation_bands), len(SWE.time.values)))
        SWE_bands = pd.DataFrame(SWE_bands,
                                 index=elevation_bands,
                                 columns=SWE.time.values)
        for i in trange(len(SWE.time.values)):
            tile_date = SWE.time.values[i]
            for j in range(len(masks)):
                band = masks[j]
                tile_band = SWE.sel(time=tile_date).where(band)
                snow_band = tile_band.sum().item()
                SWE_bands.iloc[j, i] = snow_band/np.count_nonzero(band)
        SWE_bands.to_csv("datos/ANDES_SWE_Cortes/snowbands_maipomanzano_" +
                         str(yr)+".csv")

# mask = np.empty(SWE_bands.shape[1],dtype=bool)
# for i in range(SWE_bands.shape[1]):
#     if "." in SWE_bands.columns[i]:
#         mask[i]=False
#     else:
#         mask[i]=True

# SWE_bands = SWE_bands.iloc[:,mask]
# %%


var = SWE_bands  # .loc[SWE_bands.index<4.5e3,:]
H50 = np.ones((var.shape[1]))*np.nan
H20 = np.ones((var.shape[1]))*np.nan
H80 = np.ones((var.shape[1]))*np.nan


for i in trange(var.shape[1]):
    idx1 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.2*var.max(axis=0)[i]))).flatten()
    idx2 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.5*var.max(axis=0)[i]))).flatten()
    idx3 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.8*var.max(axis=0)[i]))).flatten()

    if var.iloc[:, i].mean() != 0.0:
        interp1 = interp1d([var.iloc[idx1[0], i], var.iloc[idx1[0]+1, i]],
                           [var.index[idx1[0]], var.index[idx1[0]+1]])
        interp2 = interp1d([var.iloc[idx2[0], i], var.iloc[idx2[0]+1, i]],
                           [var.index[idx2[0]], var.index[idx2[0]+1]])
        interp3 = interp1d([var.iloc[idx3[0], i], var.iloc[idx3[0]+1, i]],
                           [var.index[idx3[0]], var.index[idx3[0]+1]])
        try:
            H20[i] = interp1(0.2)
        except:
            H20[i] = np.nan
        try:
            H50[i] = interp2(0.5)
        except:
            H50[i] = np.nan
        try:
            H80[i] = interp3(0.8)
        except:
            H80[i] = np.nan

    # H20[i] = var.index[idx1[0]]
    # H50[i] = var.index[idx2[0]]
    # H80[i] = var.index[idx3[0]]
H20 = pd.Series(H20, index=SWE_bands.columns)
H50 = pd.Series(H50, index=SWE_bands.columns)
H80 = pd.Series(H80, index=SWE_bands.columns)
# H80  = SWE_bands.max(axis=0)*0.8
# dH = pd.DataFrame(np.gradient(SWE_bands)[0]).rolling(10,center=True,min_periods=2).mean().max(axis=0).values
dH = H80-H20


# %%
time = pd.date_range("1984-04-01", "2021-10-10", freq="d")
snowlimits = pd.read_csv("datos/snowlimits_maipomanzano.csv", index_col=0)
names = ["CORTES_H20", "CORTES_H50", "CORTES_H80"]
snowlimits = snowlimits.reindex(time)
snowlimits[names[0]] = H20.reindex(time)
snowlimits[names[1]] = H50.reindex(time)
snowlimits[names[2]] = H80.reindex(time)
