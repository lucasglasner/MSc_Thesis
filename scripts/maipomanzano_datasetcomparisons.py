#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 11:14:35 2021

@author: lucas

# =============================================================================
# Comparison of datasets in relationship with the best representation of 
# real snow limit and freezing level on RioMaipoEnElManzano basin.
# =============================================================================


"""


import sys
sys.path.append('scripts/')
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from taylorDiagram import TaylorDiagram, test1, test2
#%%

isotermas0 = pd.read_csv("datos/isotermas0_maipomanzano.csv",index_col=0)
isotermas0.index = pd.to_datetime(isotermas0.index)

snowlimits = pd.read_csv("datos/snowlimits_maipomanzano.csv",index_col=0)
snowlimits.index = pd.to_datetime(snowlimits.index)

snowcovers = pd.read_csv('datos/snowcovers_maipomanzano.csv',index_col=0)
snowcovers.index = pd.to_datetime(snowcovers.index)


#%%

from scipy.stats import linregress


fig,ax = plt.subplots(2,3,sharex=True,sharey=True,figsize=(10,5))
fig.tight_layout(pad=3)
fig.text(0,0.5,"Snow Limit height (m)",ha="center",va="center",rotation=90)

# box1 = []
box2 = []

for axis in ax.ravel()[:-1]:
    bbox = axis.get_position()
    # box1.append(fig.add_axes([bbox.xmin,bbox.ymax,bbox.xmax-bbox.xmin,0.1],sharex=axis))
    box2.append(fig.add_axes([bbox.xmax,bbox.ymin,0.05,bbox.ymax-bbox.ymin],sharey=axis))
    

var = [snowlimits.iloc[:,i] for i in range(6)]
names = ["MODIS_H20","MODIS_H50","MODIS_H80","DGF","HYPSOMETRY\nIANIGLA"]
colors = plt.cm.get_cmap("tab10",len(names))(np.linspace(0,1,len(names)))
ax = ax.ravel()
# ax[4].plot([],[],ls=":",color="k",label="Linear\nRegression")
ax[5].plot([],[],ls=":",color="r",label="$y\sim x$")
linmodels=[]
for i in range(len(ax)-1):
    x,y = var[5].dropna(),var[i].dropna()
    x   = x.reindex(y.index)
    y   = y.reindex(x.index)
    ax[i].scatter(x,
                  y,edgecolor="k",alpha=0.6,color=colors[i],zorder=3,
                  lw=0.4,label=names[i])
    m = linregress(x,y)
    ax[i].plot(np.arange(800,7.5e3),
                np.arange(800,7.5e3)*m.slope+m.intercept,ls=":",color="k")
    t=ax[i].text(x=0.5,y=0.02,
                  s="y~"+"{:.2f}".format(m.slope)+"x"+"{:.2f}".format(m.intercept)+"\n$R^2$: "+"{:.1%}".format(m.rvalue**2),
                  transform=ax[i].transAxes)
    ax[i].plot([1e2,7.5e3],[1e2,7.5e3],color="red",ls=":")
    ax[i].set_xlim(0,7.5e3)
    ax[i].set_ylim(0,7.5e3)
    # ax[i].set_title(names[i],loc="left")
    ax[i].grid(True,ls=":")
    # ax[5].scatter([],[],label=names[i],color=colors[i],edgecolor="k",lw=0.4)
    
    # box1[i].axis("off")
    box2[i].axis("off")
    # if i==len(ax)-2:
        # box1[i].boxplot(x,vert=False)
    box2[i].boxplot(y,vert=True)
    linmodels.append(m)
    ax[i].legend(frameon=False,loc=(0,1))
    # box1[i].boxplot(x,vert=True)

# for i in [1,2,4]:
#     ax[i].set_yticklabels("")
ax[5].axis("off")
# ax[4].legend(frameon=False,loc=(1.85,0.77))
ax[5].legend(frameon=False,loc=(0.1,0.8),ncol=1)
# ax[5].scatter(sl_ianigla.loc[tile_props.index],dH,edgecolor="k",alpha=0.6)
# ax[5].set_title("dH",loc="left")

ax[4].set_xlabel("Snow Limit IANIGLA (m)")

plt.savefig("plots/maipomanzano/datasetcomparison/scatterplots_snowlimitS.pdf",dpi=150,bbox_inches="tight")