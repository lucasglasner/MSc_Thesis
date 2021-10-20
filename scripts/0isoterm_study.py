#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:44:52 2021

@author: lucas

# =============================================================================
# Compare 0°C Isotherm Height from Wyoming radiosonde data, amdar data base, and
# station data
# =============================================================================

"""


import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, splev, splrep, BSpline
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
#%%
# =============================================================================
# Load data
# =============================================================================

#load data from sto domingo radiosonde
H0_stodomingo = pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)
H0_stodomingo.index = pd.to_datetime(H0_stodomingo.index)
H0_stodomingo = H0_stodomingo.where(H0_stodomingo<7e3).dropna()["H0_StoDomingo"]

#load data from AMDAR database
H0_amdar = pd.read_csv("datos/amdar/isoterma0SCEL20172019.csv",index_col=0)
func = lambda x: datetime.fromordinal(int(x)) + timedelta(days=x%1) - timedelta(days = 366)
H0_amdar.index = H0_amdar["timeregh"].map(func)
H0_amdar = H0_amdar["zt0"]
H0_amdar.index = pd.date_range("2017-01-01T01:00:00","2020-01-01",freq="h").dropna()
# H0_amdar.index = H0_amdar.index.map(lambda x: x+timedelta(hours=4))

#load data from precipitation
pr_qn = pd.read_csv("datos/estaciones/pr_quintanormal.csv",dtype=str)
pr_qn.index = pr_qn.iloc[:,0]+"-"+pr_qn.iloc[:,1]+"-"+pr_qn.iloc[:,2]
pr_qn.index = pd.to_datetime(pr_qn.index)
pr_qn = pd.to_numeric(pr_qn.drop(pr_qn.columns[[0,1,2]],axis=1).iloc[:,0])

pr_lo = pd.read_csv("datos/estaciones/pr_laobra.csv",dtype=str)
pr_lo.index = pr_lo.iloc[:,0]+"-"+pr_lo.iloc[:,1]+"-"+pr_lo.iloc[:,2]
pr_lo.index = pd.to_datetime(pr_lo.index)
pr_lo = pd.to_numeric(pr_lo.drop(pr_lo.columns[[0,1,2]],axis=1).iloc[:,0])

data_dgf = pd.read_csv("datos/estaciones/dgf/DATOSUTC_2004-2019.csv",index_col=0)
data_dgf.index = pd.to_datetime(data_dgf.index.values)
#%%

fig = plt.figure(num=0,figsize=(8,6),dpi=150)
fig.tight_layout(pad=2)
ax  = fig.add_subplot(211)
ax1 = fig.add_subplot(223)
ax2 = fig.add_subplot(224)
fig.tight_layout(pad=3)
ax.plot(H0_stodomingo,color="powderblue")
ax.plot(H0_amdar,color="wheat")
ax.set_ylabel("Isotherm 0°C height (m)")



hamd = H0_amdar.resample("12h").interpolate()
hsto = H0_stodomingo.resample("12h").asfreq().reindex(hamd.index)
mask = data_dgf["9"].resample("12h").sum().reindex(hsto.index)>1

ax1.boxplot(hsto.dropna(),sym="",positions=[0],showmeans=True,meanline=True,
            patch_artist=True,
            boxprops={"facecolor":"powderblue"},
            meanprops={"linestyle":"--","color":"k"},
            medianprops={"color":"tab:red"})

ax1.boxplot(hsto[mask].dropna(),
            sym="",positions=[0.75],showmeans=True,meanline=True,
            patch_artist=True,
            boxprops={"facecolor":"powderblue"},
            meanprops={"linestyle":"--","color":"k"},
            medianprops={"color":"tab:red"})

ax1.boxplot(hamd.dropna(),sym="",positions=[0.25],showmeans=True,meanline=True,
            patch_artist=True,
            boxprops={"facecolor":"wheat"},
            meanprops={"linestyle":"--","color":"k"},
            medianprops={"color":"tab:red"})


ax1.boxplot(hamd[mask].dropna(),
            sym="",positions=[1],showmeans=True,meanline=True,
            patch_artist=True,
            boxprops={"facecolor":"wheat"},
            meanprops={"linestyle":"--","color":"k"},
            medianprops={"color":"tab:red"})
ax1.scatter([],[],color="powderblue",marker="s",label="Sto. Domingo\nRadiosonde",edgecolor="k")
ax1.scatter([],[],color="wheat",marker="s",label="AMDAR",edgecolor="k")
ax1.legend(frameon=False,ncol=2,loc=(0,1))
ax1.grid(True,axis="y",ls=":")
ax1.set_xticks([0.25/2,(0.75)+0.25/2])
ax1.set_xticklabels(["H0","H0\nPrecip.>1mm"])
ax1.set_ylabel("Isotherm 0°C height (m)")

ax2.plot([1e1,7e3],[1e1,7e3],color="tab:red",ls="--",zorder=1)
ax2.scatter(x=hamd[mask].values,
            y=hsto[mask].values,
            edgecolor="k",zorder=2,alpha=.9,color="gold",s=20,label="Precip.>1mm")
ax2.scatter(x=hamd.values,
            y=hsto.values,
            edgecolor="k",zorder=1,alpha=.5,color="tab:blue",s=20)
ax2.set_ylabel("Sto. Domingo")
ax2.set_xlabel("AMDAR")
ax2.legend(frameon=False)



ax3 = fig.add_subplot(3,9,27)
diff = hamd[mask]-hsto[mask]
diff = diff.dropna()
ax3.boxplot(diff,sym="",positions=[0],patch_artist=True,
            boxprops={"facecolor":"gold"},
            medianprops={"color":"k"})
ax3.set_title("$\Delta H0$",loc="center")
ax3.set_xticks([])
ax3.set_xticklabels([])
ax3.set_yticks(np.arange(-150,750,150))
# ax3.grid(True,axis="y",ls=":")
ax3.yaxis.tick_right()
plt.savefig("plots/maipomanzano/isotherm0_study.pdf",dpi=150,bbox_inches="tight")






