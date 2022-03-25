#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 09:50:34 2022

@author: lucas

# =============================================================================
# Analysis of CORTES SCA in different basins.
# =============================================================================
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from glob import glob
import re
import xarray as xr

#%%
# =============================================================================
# ianigla
# =============================================================================

ianigla_paths = glob('datos/ianigla/Rio*as.csv')
paths = [ianigla_paths[p] for p in [0,4,3,1,7,6,2,5,8]]
ianigla_paths = paths
basins = [p.split("/")[-1].split("_")[0] for p in ianigla_paths]

SCA_modis = [pd.read_csv(p) for p in ianigla_paths]
for file in SCA_modis:
    file.index = pd.to_datetime(file.fecha)
    
SCA_modis = pd.concat(SCA_modis,axis=1)['SCA(%)'].dropna()/100
# SCA_modis = SCA_modis.drop_duplicates()
SCA_modis.columns = basins

SCA_modis_trend = SCA_modis.diff()

#%%
# =============================================================================
# cortes
# =============================================================================

th = 100
SCA_cortes = []
for b in basins:
    path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/'+b+'/ANDES_SWE_1985-2015.nc'
    SWE = xr.open_dataset(path).SWE
    tpix = ~np.isnan(SWE.sel(time='2013/08/09'))
    tpix = np.count_nonzero(tpix)
    sca = xr.where(SWE>th, 1, 0).sum(dim=['lon','lat']).to_series()/tpix
    SCA_cortes.append(sca)

SCA_cortes = pd.concat(SCA_cortes,axis=1).dropna()
SCA_cortes.columns = basins

#%%
time_range = pd.date_range('2000-02-25','2015-03-31',freq='d')

diffs = (SCA_cortes.reindex(time_range)-SCA_modis.reindex(time_range)).dropna()
diffs = diffs.where(SCA_modis_trend.reindex(time_range)!=0).dropna(how='all')

seasons = np.array(diffs.index.dayofyear)

# "day of year" ranges for the northern hemisphere
spring = range(80, 172)
summer = range(172, 264)
fall = range(264, 355)
# winter = everything else


for i,doy in enumerate(seasons):
    if doy in spring:
      season = 1
    elif doy in summer:
      season = 2
    elif doy in fall:
      season = 3
    else:
      season = 4
    seasons[i] = season
     
seasons = list(map(lambda x: str(x),seasons))
seasons = ".".join(seasons).replace("1","spring")
seasons = seasons.replace("2","summer").replace("3","fall")
seasons = seasons.replace("4","winter").split(".")
seasons = np.array(seasons)

diffs = diffs.groupby([diffs.index,seasons]).mean()

#%%
fig,ax = plt.subplots(1,1, figsize=(12,3))

plt.rc('font',size=18)
pad = [-0.5,0,0.5,1]
colors = ['tab:blue','tab:green','tab:red','tab:orange']
for i,b in enumerate(basins):
    basin_data = diffs[b]
    for j,s in zip([0,1,2,3],
                   ['winter','spring','summer','fall']):
        ax.boxplot(basin_data.unstack()[s].dropna(),
                   positions=[i*4+pad[j]], sym="",
                   widths=0.5,
                   patch_artist=True,
                   boxprops={'facecolor':colors[j]},
                   medianprops={'color':'k'})
ax.set_xticks(np.arange(0,9,1)*4)
ax.set_xticklabels(['\n'.join(re.findall('[A-Z][^A-Z]*', b))
                    for b in basins])
ax.tick_params(axis="both", labelsize=14)
ax.grid(axis='y', ls=":", which='major')
ax.set_yticks(np.arange(-1,1+0.3,0.2))
ax.set_ylim(-0.85,0.85)
ax.set_ylabel('SCA Bias (-)')
ax.axhline(0, color='k', ls="--", lw=0.8)


for c,label in zip(colors,['winter','spring','summer','autumn']):
    ax.scatter([],[], marker="s", edgecolor='k',
               color=c, label=label)
ax.legend(frameon=False, loc=(0,1), ncol=4, fontsize=14)

plt.savefig('plots/CORTES_SCA_ANALYSIS_BASINS.pdf', dpi=150,
            bbox_inches='tight')










