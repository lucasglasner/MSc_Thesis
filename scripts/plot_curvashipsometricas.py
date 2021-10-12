#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 19:38:39 2021

@author: lucas

# =============================================================================
# Hypsometry quick plot
# =============================================================================

"""

import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
#%%
paths = glob("datos/topography/basins/hipso/*.csv")
data  = {paths[i].replace("datos/topography/basins/hipso/","").replace("_Hipso.csv","")
         :pd.read_csv(paths[i])
         for i in range(len(paths))}
#%%
names = list(data.keys())
areas = [data[name]["Area_km2"].values[-1] for name in names]
areas = np.array(areas)
scale_area = (areas-np.min(areas))/(np.max(areas)-np.min(areas))
scale_area = pd.Series(scale_area,index=names).sort_values()
#%%
fig,ax = plt.subplots(1,1,figsize=(10,5))
colors = plt.cm.brg_r(np.linspace(0.1,0.9,len(names)))
for i in range(len(names)):
    name = scale_area.index[i]
    ax.plot(data[name]["height"],data[name]["fArea"],
            label=name, color=colors[i],lw=scale_area[i]*2.5+1.5)
ax.legend(loc=(1.05,0), frameon=False)
ax.set_xlabel("Height above sea level (m)",fontsize=14)
ax.set_ylabel("%Area below Height",fontsize=14)
ax.set_title("Hypsometry",fontsize=15,loc="left")

plt.savefig("plots/terrain/curvas_hipsometricas.pdf",dpi=150,bbox_inches="tight")

#%%
