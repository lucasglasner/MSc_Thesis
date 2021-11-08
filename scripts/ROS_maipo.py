#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 18:48:01 2021

@author: lucas

# =============================================================================
# ROS in Maipo Basin
# =============================================================================

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime as dt
import scipy.stats as st
#%%
# =============================================================================
# Cargar datos
# =============================================================================
cuenca           = "RioMauleEnArmerillo"
curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)

snow_limit       =  pd.read_csv("datos/ianigla/"+cuenca+"_SCA_s_comp.filtro_MA.3días.csv",index_col=0)/100
snow_limit.index =  pd.to_datetime(snow_limit.index)

isoterma0        =  pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)-300
isoterma0.index  = pd.to_datetime(isoterma0.index)
isoterma0        = isoterma0.where(isoterma0<6000).dropna()["H0_StoDomingo"]
isoterma0 = pd.Series(isoterma0.groupby([isoterma0.index.year,
                                          isoterma0.index.month,
                                          isoterma0.index.day]).mean().ravel(),
                      index=np.unique(pd.to_datetime(isoterma0.index.date)),
                      name = "H0_StoDomingo")
isoterma0        = isoterma0.resample("d").interpolate("linear")

pr_quintanormal  =  pd.read_csv("datos/estaciones/pr_quintanormal.csv").applymap(lambda x: str(x))
pr_quintanormal["fecha"] = pr_quintanormal["agno"]+"-"+pr_quintanormal[" mes"]+"-"+pr_quintanormal[" dia"]
pr_quintanormal.index = pd.to_datetime(pr_quintanormal["fecha"])
pr_quintanormal.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
pr_quintanormal  = pd.to_numeric(pr_quintanormal[" valor"])
pr_quintanormal.name = "pr"

q_maipomanzano   = pd.read_csv("datos/estaciones/q_"+cuenca+".csv").applymap(lambda x: str(x))
q_maipomanzano["fecha"] = q_maipomanzano["agno"]+"-"+q_maipomanzano[" mes"]+"-"+q_maipomanzano[" dia"]
q_maipomanzano.index = pd.to_datetime(q_maipomanzano["fecha"])
q_maipomanzano.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
q_maipomanzano  = pd.to_numeric(q_maipomanzano[" valor"])
q_maipomanzano.name = "q"
#%%
# =============================================================================
# Calcular linea de nieves
# =============================================================================

interp = interp1d(1-curva_hipso["fArea"],curva_hipso["height"])
snow_limit["SL"] = list(map(lambda x: interp(x).item(),snow_limit["SCA(%)"]))

#%%
daterange = pd.date_range("1950-01-01","2021-08-22", freq="d")
data = pd.concat((pr_quintanormal.reindex(daterange),
                  isoterma0.reindex(daterange),
                  snow_limit.reindex(daterange)),axis=1)
data = data.where(data["pr"]>1).dropna() #Dias con lluvia
data = data.where(data["SL"]<data["H0_StoDomingo"]-300).dropna()
data["q"] = q_maipomanzano.reindex(data.index)

#%%
fig,ax = plt.subplots(2,2,figsize=(10,5))
fig.tight_layout(pad=2)
fig.text(0.75,1,"ROS 2000->2017", ha="center", va="center",fontsize=20)
ax = ax.ravel()

ax[0].hist(isoterma0,bins="auto",density=True,alpha=0.7,color="goldenrod",
           label="sin_pr")
ax[0].hist(isoterma0.reindex(pr_quintanormal.index)[pr_quintanormal>1],bins="auto",
           density=True,alpha=0.7,color="royalblue", label="con_pr")
ax[0].hist(data["H0_StoDomingo"], alpha=0.7, color="skyblue", density=True,
           label="ROS")
ax[0].legend(loc=(0,1.05), frameon=False, ncol=3)
ax[0].set_xlabel("isoterma0_maipomanzano (m)")

ax[1].hist(snow_limit["SL"],bins="auto", color="goldenrod",density=True , alpha=0.7)
ax[1].hist(snow_limit["SL"].reindex(pr_quintanormal.index)[pr_quintanormal>1],
           bins="auto", density=True, alpha=0.7, color="royalblue" )
ax[1].hist(data["SL"], bins="auto" , density=True, alpha=0.7, color="skyblue")
ax[1].set_xlabel("snow_limit_maipomanzano (m)")


ax[2].hist(pr_quintanormal[pr_quintanormal>1], bins="auto", density=True, alpha=0.7,
           color="royalblue" )
ax[2].hist(data["pr"], bins="auto" , density=True, alpha=0.7, color="skyblue")
ax[2].set_xlabel("pr_quintanormal (mm)")



ax[3].hist(q_maipomanzano, bins="auto", density=True, alpha=0.7,
           color="goldenrod" )
ax[3].hist(q_maipomanzano.reindex(pr_quintanormal.index)[pr_quintanormal>1], bins="auto", density=True, alpha=0.7,
           color="royalblue")
ax[3].hist(data["q"], bins="auto" , density=True, alpha=0.7, color="skyblue")
ax[3].set_xlabel("Q_maipomanzano (m3/s)")


ax[3].axvline(st.expon.ppf(1-1/25,*st.expon.fit(q_maipomanzano)), color="k", label="25, 50, 100, 500\nreturn periods", ls=":")
ax[3].axvline(st.expon.ppf(1-1/50,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
ax[3].axvline(st.expon.ppf(1-1/100,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
ax[3].axvline(st.expon.ppf(1-1/500,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
ax[3].legend()


plt.savefig("plots/maipomanzano/ROS_maipo.pdf",dpi=150,bbox_inches="tight")

#%%

fig,ax = plt.subplots(2,1,sharex=True,figsize=(10,5))
ax = ax.ravel()
date = "2009-09-06"
step = 40
sliced = slice(pd.to_datetime(date)-dt.timedelta(days=step),
               pd.to_datetime(date)+dt.timedelta(days=step))
# xticks = pd.date_range(sliced.start, sliced.stop, freq='d')
snow_limit["SL"].loc[sliced].plot(ax=ax[0],color="royalblue",label="snow_limit")
isoterma0.loc[sliced].plot(ax=ax[0],color="dimgrey",label="isoterma 0")
ax[1].bar(pr_quintanormal.loc[sliced].index,pr_quintanormal.loc[sliced], color="teal")
ax[1].set_ylabel("pr (mm)")
ax[0].set_ylabel("z (m.s.n.m)")
ax2 = ax[1].twinx()
ax2.set_ylabel("Q (m3/s)")
q_maipomanzano.loc[sliced].rolling(window=1).mean().plot(ax=ax2, color="blueviolet", lw=2)
ax[0].legend(loc=(0,1.1), ncol=2, frameon=False)
# ax[1].legend()
# ax[2].legend()
for axis in ax:
    axis.set_xlabel("")
# ax[0].set_xticklabels([x.strftime('%m-%d') for x in xticks]);

plt.savefig("plots/maipomanzano/ROS"+date+"_maipo.pdf",dpi=150,bbox_inches="tight")
