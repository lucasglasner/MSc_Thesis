#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 16:23:54 2021

@author: lucas

# =============================================================================
# Maipo en el Manzano Basin. fSCA analisis and relation with snow limit for 
# 09/09/2009 case of study
# =============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from scipy.interpolate import interp1d
from scipy.cluster.vq import whiten
import matplotlib.colors as mplcolors
from cmcrameri import cm
import datetime
import matplotlib as mpl
import cmocean
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import geopandas as gpd
import matplotlib.patches as mpatches
#%%
# =============================================================================
# basin
# =============================================================================

basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")
# =============================================================================
# fSCA modis
# =============================================================================

fSCA_terra = xr.open_dataset("datos/modis/MOD10A1_2000-2021.nc",chunks="auto")
fSCA_aqua  = xr.open_dataset("datos/modis/MYD10A1_2000-2021.nc",chunks="auto")

# =============================================================================
# Basin hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)

# =============================================================================
# Snow limit
# =============================================================================
sl_ianigla       =  pd.read_csv("datos/ianigla/"+cuenca+"_SCA_s_comp.filtro_MA.3días.csv",index_col=0)/100
sl_ianigla.index =  pd.to_datetime(sl_ianigla.index)
interp = interp1d(1-curva_hipso["fArea"],curva_hipso["height"])
sl_ianigla["SL"] = list(map(lambda x: interp(x).item(),sl_ianigla["SCA(%)"]))
sl_ianigla = sl_ianigla["SL"]
sl_ianigla2 = pd.read_csv("datos/ianigla/RioMaipoEnElManzano_lim_nieve_ianigla_2000-2015.csv",index_col=1).iloc[:,1]
sl_ianigla2.index = pd.to_datetime(sl_ianigla2.index)


sl_dgf           =  pd.read_csv("datos/modis/MAIPO.txt",sep=" ",header=None)
sl_dgf.index     = pd.to_datetime(sl_dgf[0].values,format="%Y%m%d")
sl_dgf.drop([0,1,2,3],axis=1,inplace=True)
sl_dgf = sl_dgf[4]
# =============================================================================
# Basin orography
# =============================================================================
dem = xr.open_dataset("datos/topography/basins/"+cuenca+"_regridmodis.nc").Band1

# =============================================================================
# santo domingo
# =============================================================================

H0_sd        =  pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)-300
H0_sd.index  = pd.to_datetime(H0_sd.index)
H0_sd        = H0_sd.where(H0_sd<6000).dropna()["H0_StoDomingo"]


#%%
# =============================================================================
# Build elevation band masks
# =============================================================================

#Number of elevation bands
dz              = 50
elevation_bands = np.arange(dem.min(),dem.max()+dz,dz)
elevation_bands = np.hstack((elevation_bands))

masks = []
for i in range(len(elevation_bands)-1):
    mask = dem.where(dem>elevation_bands[i]).where(dem<elevation_bands[i+1])
    mask = ~np.isnan(mask)
    masks.append(mask.values)

#%%
# =============================================================================
# Load data to memory
# =============================================================================
date_sd = "2005-10-09"
date = "2005-10-09"
good_days1  = [datetime.datetime(int(date.split("-")[0]),
                                 int(date.split("-")[1]),
                                 int(date.split("-")[2]))]
good_days2  = [datetime.datetime(int(date.split("-")[0]),
                                 int(date.split("-")[1]),
                                 int(date.split("-")[2]))]

good_terra  = fSCA_terra.fSCA.sel(time=date).load()
good_aqua   = fSCA_aqua.fSCA.sel(time=date).load()


#%%
# =============================================================================
# Apply masks
# =============================================================================


fSCA_limits = [10,40,70]
# fSCA_limits = np.arange(30,80,10)
terra_singleimage = np.empty((len(masks),len(fSCA_limits)))
aqua_singleimage  = np.empty((len(masks),len(fSCA_limits)))
for i in range(len(masks)):
    print("Elevation Band "+str(i)+": "+"{:.1f}".format(elevation_bands[i])+"-"+"{:.1f}".format(elevation_bands[i+1]))
    mask1 = np.moveaxis(np.tile(np.expand_dims(masks[i],2),len(good_days1)),2,0)
    mask2 = np.moveaxis(np.tile(np.expand_dims(masks[i],2),len(good_days2)),2,0)
    for j in range(len(fSCA_limits)):
        tile_band1 = good_terra.values[mask1.squeeze()]
        tile_band2 = good_aqua.values[mask2.squeeze()]
        snow_band1 = np.where((tile_band1<100) & (tile_band1>fSCA_limits[j]),1,0)
        snow_band2 = np.where((tile_band2<100) & (tile_band2>fSCA_limits[j]),1,0)
        # fSCA_bands.iloc[j,i] = snow_band.sum()/len(snow_band)
        terra_singleimage[i,j] = snow_band1.sum()/len(snow_band1)
        aqua_singleimage[i,j]  = snow_band2.sum()/len(snow_band2)
        # terra_singleimage[i,j] = np.where(((np.where(mask1,good_terra,0)>fSCA_limits[j]) &
        #                              ((np.where(mask1,good_terra,0)>fSCA_limits[j]))),1,0).sum(axis=1).sum(axis=1)/mask1.sum()
        # aqua_singleimage[i,j]  = np.where(((np.where(mask2,good_aqua,0)>fSCA_limits[j]) &
        #                              ((np.where(mask2,good_aqua,0)>fSCA_limits[j]))),1,0).sum(axis=1).sum(axis=1)/mask2.sum()
        
#%%
# =============================================================================
# Plot8
# =============================================================================

fig,ax=plt.subplots(2,1,sharex=True,sharey=True,figsize=(5,8))
fig.text(0.03,0.5,"Probability of snow covered band",rotation=90,ha="center",va="center",fontsize=16)
fig.text(0.5,0.95,"%Snow cover band distribution",ha="center",va="center",fontsize=16)

colors = cmocean.cm.ice(np.linspace(0.15,.85,len(fSCA_limits)))
y = list(map(lambda x: int(x)-1,elevation_bands))
x=pd.cut(y,y).categories.values
x=list(map(lambda j: str(j),x))
x=elevation_bands[:-1]
for i in range((terra_singleimage.shape[1])):
    ax[0].step(x,terra_singleimage[:,i],color=colors[i],alpha=1,where="post")
    lb = "$fSCA_{xy}\geq$"+str(fSCA_limits[i])+"%"
    ax[1].step(x,aqua_singleimage[:,i],color=colors[i],alpha=1,where="post",label=lb)

for axis in ax:
    plt.xticks(rotation=45)
    axis.xaxis.set_major_locator(MultipleLocator(250))
    axis.xaxis.set_minor_locator(MultipleLocator(50))
    axis.grid(True,ls=":",which="major")
    axis.set_xlim(1e3,4.5e3)
    axis.set_yticks(np.arange(0,1.1,0.1))
lg=ax[1].legend(loc="lower right",title="Snow threshold",frameon=True)

ax[0].set_title("MODIS/TERRA",loc="left")
ax[1].set_title("MODIS/AQUA",loc="left")
ax[1].set_xlabel("Elevation Band (m)",fontsize=16)
axc = ax[0].get_position()
ax1 = fig.add_axes([axc.xmax*1.2,axc.ymin,(axc.xmax-axc.xmin)/3,axc.ymax-axc.ymin])
good_terra.plot(ax=ax1,cmap=cmocean.cm.ice,vmin=0,vmax=100,add_colorbar=False,rasterized=True)
good_terra.where(good_terra>100).plot.pcolormesh(ax=ax1,add_colorbar=False,vmin=100,vmax=100,cmap="bwr",rasterized=True)
basin.boundary.plot(ax=ax1,color="k",lw=0.5)
ax1.axis("off")
fig.text((ax1.get_position().xmax-ax1.get_position().xmin)/2+ax1.get_position().xmin,
         ax1.get_position().ymax*1.05,
         "TIME: "+date,
         ha="center",
         va="center")
patch=mpatches.Patch(color="red",label="NoData")
ax1.legend(frameon=False,loc=(0,1),handles=[patch])

ax1.annotate('', xy=(-0.7, 1),  xycoords='axes fraction',
             xytext=(0., 0.5), textcoords='axes fraction',
             arrowprops=dict(facecolor='black',arrowstyle="->"),
             horizontalalignment='right', verticalalignment='top',
             )
ax1.annotate('', xy=(-0.7, 0),  xycoords='axes fraction',
             xytext=(0., 0.5), textcoords='axes fraction',
             arrowprops=dict(facecolor='black',arrowstyle="->"),
             horizontalalignment='right', verticalalignment='top',
             )
ax1.set_title("")


axc = ax[1].get_position()
ax2 = fig.add_axes([axc.xmax*1.2,axc.ymin,(axc.xmax-axc.xmin)/3,axc.ymax-axc.ymin])
good_aqua.plot(ax=ax2,cmap=cmocean.cm.ice,vmin=0,vmax=100,add_colorbar=False,
               rasterized=True)
good_aqua.where(good_aqua>100).plot.pcolormesh(ax=ax2,add_colorbar=False,
                                               vmin=100,vmax=100,cmap="bwr",rasterized=True)
basin.boundary.plot(ax=ax2,color="k",lw=0.5)
ax2.axis("off")
ax2.set_title("")

ax2.annotate('', xy=(-0.7, 1),  xycoords='axes fraction',
             xytext=(0.0, 0.5), textcoords='axes fraction',
             arrowprops=dict(facecolor='black',arrowstyle="->"),
             horizontalalignment='right', verticalalignment='top',
             )
ax2.annotate('', xy=(-0.7, 0),  xycoords='axes fraction',
             xytext=(0., 0.5), textcoords='axes fraction',
             arrowprops=dict(facecolor='black',arrowstyle="->"),
             horizontalalignment='right', verticalalignment='top',
             )


cax = fig.add_axes([axc.xmax*1.6,axc.ymin,(axc.xmax-axc.xmin)/25,2*(axc.ymax-axc.ymin)*1.1])

cmap = cmocean.cm.ice
norm = mplcolors.Normalize(vmin=0,vmax=100)
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label("Fraction of Snow Cover Area (%)",fontsize=15)

plt.savefig("plots/maipomanzano/snowthreshold_fSCAstudy.pdf",dpi=150,bbox_inches="tight")
#%%
# =============================================================================
# Plot2
# =============================================================================


#%%
# =============================================================================
# Save output
# =============================================================================
# columns = []
# for i in range(len(elevation_bands)-1):
#     columns.append("("+str(elevation_bands[i])+"-"+str(elevation_bands[i+1])+")")
    
# # aqua_singleimage  = pd.DataFrame(aqua_singleimage.T,index=good_days2,columns=columns)
# # terra_singleimage = pd.DataFrame(terra_singleimage.T,index=good_days1,columns=columns)

# aqua_singleimage  = pd.DataFrame(aqua_singleimage.T,index=fSCA_limits,columns=columns).T
# terra_singleimage = pd.DataFrame(terra_singleimage.T,index=fSCA_limits,columns=columns).T





