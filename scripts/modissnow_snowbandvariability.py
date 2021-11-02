#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 15:50:06 2021

@author: lucas

# =============================================================================
# Maipo en el Manzano Basin. fSCA analisis and relation with snow band
# distribution
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
import geopandas as gpd
from tqdm import trange
import numba
import cmocean
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.patches as mpatches
#%%
# =============================================================================
# fSCA modis
# =============================================================================
sensor = "terra"
if sensor=="terra":
    fSCA = xr.open_dataset("datos/modis/MOD10A1_2000-2021.nc").fSCA
elif sensor=="aqua":
    fSCA  = xr.open_dataset("datos/modis/MYD10A1_2000-2021.nc").fSCA
else:
    raise RuntimeError("Sensor can only be terra or aqua")

# =============================================================================
# Basin hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)
basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")
# =============================================================================
# snow limit
# =============================================================================
# sl_ianigla2 = [3298,3079,2718]


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
# Basin aspect
# =============================================================================

aspect = xr.open_dataset("datos/topography/basins/aspect/"+cuenca+"_regridmodis.nc").Band1

# =============================================================================
# Precipitation in stgo
# =============================================================================
pr_lo = pd.read_csv("datos/estaciones/pr_laobra.csv").applymap(lambda x: str(x))
pr_lo.index = pr_lo.iloc[:,0]+"-"+pr_lo.iloc[:,1]+"-"+pr_lo.iloc[:,2]
pr_lo.index = pd.to_datetime(pr_lo.index)
pr_lo       = pd.to_numeric(pr_lo.iloc[:,3])

#%%
# =============================================================================
# Read file with "good images", or make it.
# Good images are defined as the ones with less than 20% cloud cover,
# 70% of basin area covered with snow and where last rainfall happen in
# less than 14 days
# =============================================================================

try:
    tnans   = 35811
    total   = 82501
    # =========================================================================
    # Read file with tile props and filter image collection
    # =========================================================================
    tile_props = pd.read_csv("datos/modis/modis_"+sensor+"_tileprops.csv",
                             index_col=0)
    tile_props.index = pd.to_datetime(tile_props.index)
    
except:
    # =========================================================================
    # compute cloudiness and snow cover for each tile
    # =========================================================================
    tile_props = np.empty(len(fSCA.time))
    tile_props = pd.DataFrame(tile_props,index=fSCA.time.values,columns=["%clouds"])
    tile_props["%SCA"] = np.empty(len(fSCA.time))
    
    tnans   = 35811
    total   = 82501
    for i in trange(len(tile_props.index)):
        index = tile_props.index[i]
        if not np.isnan(fSCA[i,100,100].item()):
            basin_mask = ~np.isnan(fSCA[i,:,:].values)
            n_clouds   = (fSCA[i,:,:].values[basin_mask]>100).sum()
            n_snow     = np.logical_and(fSCA[i,:,:].values[basin_mask]>30,
                                        fSCA[i,:,:].values[basin_mask]<100).sum()
            tile_props.loc[index,"%clouds"] = (n_clouds)/(total-tnans)
            tile_props.loc[index,"%SCA"]    = (n_snow)/(total-tnans)
            del basin_mask,n_clouds,n_snow
            
        else:
            tile_props.loc[index] = np.nan
    # =========================================================================
    # compute time to last rainfall in the valley
    # =========================================================================
    rainy_days = pr_lo[pr_lo>5].index
    rainy_days = np.array(list(map(lambda x: x.toordinal(),rainy_days)))
    
    snowy_days = tile_props.index
    snowy_days = np.array(list(map(lambda x: x.toordinal(),snowy_days)))
    pairwise_difference = np.array([[j-i for i in rainy_days] for j in snowy_days])
    
    good_days = np.where(pairwise_difference>0,pairwise_difference,9999).min(axis=1)
    
    tile_props["timetolastrain"] = good_days
    tile_props["rainfall_day"]   = [tile_props.index[i]-datetime.timedelta(days=int(good_days[i]))
                                    for i in range(len(good_days))]
    tile_props = tile_props.sort_values(by="timetolastrain")
    tile_props.to_csv("datos/modis/modis_"+sensor+"_tileprops.csv")
    
# =========================================================================
# Grab tiles with less than 20% cloud cover, more than 70% snow cover,
# and where last rainfall in the basin outlet happen in less than 30 days
# =========================================================================

tile_props = tile_props[tile_props["%clouds"]<0.05].dropna()
tile_props = tile_props[tile_props["%SCA"]>0.5].dropna()
tile_props = tile_props[tile_props["timetolastrain"]<30].dropna()
fSCA   = fSCA.sel(time=tile_props.index)
# fSCA=fSCA.sel(time=slice(*["2005-10","2005-10"]))


#%%
# =============================================================================
# Read file with snow cover by band or make it.
# =============================================================================
try: 
    fSCA_bands = pd.read_csv("datos/modis/modis_"+sensor+"_fSCA_bands.csv",
                             index_col=0)
    fSCA_bands.columns = pd.to_datetime(fSCA_bands.columns)
except:
    # =========================================================================
    # Build elevation band masks
    # =========================================================================
    dz              = 50
    elevation_bands = np.arange(dem.min(),dem.max()+dz,dz)
    
    masks = []
    for j in range(len(elevation_bands)-1):
            z0   = elevation_bands[j]
            z1   = elevation_bands[j+1]
            mask = (dem.where((dem>z0) & (dem<z1))>0).values
            masks.append(mask)
    
    elevation_bands = elevation_bands[:-1]
    # =========================================================================
    # Apply band mask, and compute % of snow cover distribution by band
    # =========================================================================
    
    fSCA_bands      = np.empty((len(elevation_bands),len(tile_props.index)))
    fSCA_bands      = pd.DataFrame(fSCA_bands,
                                   index = elevation_bands,
                                   columns=tile_props.index)
    for i in trange(len(tile_props.index)):
        tile_date = tile_props.index[i]
        for j in range(len(masks)):
            band = masks[j] 
            tile_band = fSCA.sel(time=tile_date).values[band]
            snow_band = np.where((tile_band<100) & (tile_band>30),1,0)
            fSCA_bands.iloc[j,i] = snow_band.sum()/len(snow_band)
    fSCA_bands.to_csv("datos/modis/modis_"+sensor+"_fSCA_bands.csv")
    
fSCA_bands = fSCA_bands.T.loc[tile_props.index].T
fSCA_bands = fSCA_bands.iloc[:,(fSCA_bands.min()<0.1).values]
elevation_bands = fSCA_bands.index.values
tile_props = tile_props.loc[fSCA_bands.columns]
#%%

#Dates of tiles to show, and respective snowband distribution
dates = ["2015-12-10","2006-11-10","2005-10-08"]
var   = fSCA_bands.loc[:,dates]

cmap = cmocean.cm.ice #color for maps
# colors=cmocean.cm.ice(np.linspace(0.2,0.8,var.shape[1]))
colors=["darkorange","limegreen","darkviolet"] #color for plots

#for x tick as intervals
y = list(map(lambda x: int(x)-1,elevation_bands))
x=pd.cut(y,y).categories.values
x=list(map(lambda j: str(j),x))
x = elevation_bands[:-1].astype(int)//10*10
#create figure
fig = plt.figure(figsize=(5,5))

#and axis
ax0 = fig.add_axes([0.0,0.23,0.7,.6])
cax = ax0.get_position()
cax = fig.add_axes([1,cax.ymin,0.03,cax.ymax-cax.ymin])

ax1 = fig.add_axes([0.6,0.00,0.5,0.33])
ax2 = fig.add_axes([0.6,0.33,0.5,0.33])
ax3 = fig.add_axes([0.6,0.66,0.5,0.33])

lon,lat = fSCA.lon,fSCA.lat
lon2d,lat2d = np.meshgrid(lon,lat)
for i in range(3):
    #plot maps
    ax = eval("ax"+str(i+1))
    tile = fSCA.sel(time=dates[i]).squeeze()
    ax.pcolormesh(lon2d,lat2d,tile,cmap=cmap,vmin=0,vmax=100,
                  rasterized=True,shading="auto")
    ax.pcolormesh(lon2d,lat2d,tile.where(tile>100),vmin=100,vmax=100,
                  cmap="bwr_r",shading="auto",rasterized=True)
    ax.axis("off")
    basin.boundary.plot(ax=ax,color=colors[i],lw=1)


    #plot curves
    label = dates[i]+": $\Delta t_{LAST\_RAIN}$: "
    label = label+str(tile_props["timetolastrain"].loc[dates[i]].item())+" days"
    ax0.step(x,var.iloc[:-1,i].values,where="post",color=colors[i],label=label)
    
    #plot hipsometry method SL ianigla
    sl_data   = sl_ianigla.loc[dates[i]]
    interp    = interp1d(elevation_bands,var.iloc[:,i])   
    ax0.scatter(sl_data,interp(sl_data),color=colors[i],
                edgecolor="k",marker="D",s=50,zorder=3,alpha=0.8)

    #plot SL ianigla
    sl_data   = sl_ianigla2.loc[dates[i]]
    interp    = interp1d(elevation_bands,var.iloc[:,i])   
    ax0.scatter(sl_data,interp(sl_data),color=colors[i],
                edgecolor="k",marker="s",s=50,zorder=3,alpha=0.8)
    
    #plot hypsometry method DGF
    sl_data   = sl_dgf.loc[dates[i]]
    interp    = interp1d(elevation_bands,var.iloc[:,i])   
    ax0.scatter(sl_data,interp(sl_data),color=colors[i],
                edgecolor="k",marker="o",s=50,zorder=3,alpha=0.8)
    
    #H50
    # sl_data   = H50.loc[dates[i]]
    # interp    = interp1d(elevation_bands,var.iloc[:,i])   
    # interp2   = interp1d(var.iloc[:-1,i],np.arange(len(x)))
    # ax0.scatter(interp2(interp(sl_data)),interp(sl_data),color=colors[i],
    #             edgecolor="k",marker="P",s=100,zorder=3,alpha=0.8)
    
    
ax1.scatter([],[],marker="o",label="DGF",color="w",edgecolor="k")
ax1.scatter([],[],marker="s",label="IANIGLA",color="w",edgecolor="k")
ax1.scatter([],[],marker="D",label="HYPSOMETRY\nIANIGLA",color="w",edgecolor="k")
# ax1.scatter([],[],marker="P",label="50%",color="w",edgecolor="k")

ax1.legend(frameon=False,loc=(-2.35,2.5),labelspacing=0.7)

# ax0.scatter([],[],marker="o",label="SL_Hypsometry_DGF")
  
# ax0.set_xticks(x)
# ax0.set_xticklabels(x)
ax0.xaxis.set_major_locator(MultipleLocator(250))
ax0.xaxis.set_minor_locator(MultipleLocator(50))
ax0.grid(True,ls=":",which="major")
ax0.set_xlim(1e3,4.5e3)
ax0.tick_params(axis="x",rotation=45)
ax0.legend(frameon=False,loc=(-0.18,1),ncol=1)
ax0.set_xlabel("Elevation Band (m)",fontsize=12)
ax0.set_ylabel("Probability of snow covered band",fontsize=12)
ax0.set_yticks(np.arange(0,1.1,0.1))
# ax0.axhline(0.5,ls=":",color="grey")
# ax0.axhline(0.8)
# ax0.axhline(0.2)

patch=mpatches.Patch(color="red",label="NoData")
ax3.legend(frameon=False,loc=(1.35,0.7),handles=[patch])

norm = mplcolors.Normalize(vmin=0,vmax=100)
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label("Fraction of Snow Cover Area (%)",fontsize=10)


# plt.savefig("plots/maipomanzano/timedifferences_fSCAstudy.pdf",dpi=150,bbox_inches="tight")
    
    

#%%

# =============================================================================
# Cut distribution curves to 5000 height maximum and minimum band with less
# less than 0.2 snow cover (to stay with tiles only where the snow limit
# is IN the basin)
# =============================================================================

# fSCA_bands = fSCA_bands.rolling(5,center=True,min_periods=2,axis=0).mean()


# fSCA_bands = fSCA_bands.iloc[fSCA_bands.index<=5.0e3,:]


var = fSCA_bands.iloc[:,:].iloc[fSCA_bands.index<=4.5e3,:]
H50 = np.ones((var.shape[1]))*np.nan
H20 = np.ones((var.shape[1]))*np.nan
H80 = np.ones((var.shape[1]))*np.nan

for i in range(var.shape[1]):
    interp = interp1d(var.iloc[:,i].values,
                      var.index,assume_sorted=False,fill_value="extrapolate")
    # if var.iloc[:,i].values.min()<0.2:
    #     if var.iloc[:,i].values.max()>0.8:
    H50[i] = interp(0.5*var.max(axis=0)[i])
    H20[i] = interp(0.2*var.max(axis=0)[i])
    H80[i] = interp(var.max(axis=0)[i]*0.8)
    # H80[i] = interp(var.iloc[:,i].sort_values().max()*0.8)

H50 = pd.Series(H50,index=tile_props.index)
H20 = pd.Series(H20,index=tile_props.index)
H80 = pd.Series(H80,index=tile_props.index)
# H80  = fSCA_bands.max(axis=0)*0.8
# dH = pd.DataFrame(np.gradient(fSCA_bands)[0]).rolling(10,center=True,min_periods=2).mean().max(axis=0).values
dH = H80-H20

#%%
# =============================================================================
# 
# =============================================================================
from scipy.stats import linregress

tile_props["dH"]  = dH
tile_props["H50"] = H50
tile_props["H80"] = H80
tile_props["H20"] = H20


fig,ax = plt.subplots(2,3,sharex=True,sharey=True,figsize=(8,4))
fig.tight_layout(pad=1)
var = [H20,H50,H80,sl_dgf.loc[tile_props.index],sl_ianigla.loc[tile_props.index]]
names = ["H20","H50","H80","DGF","HYPSOMETRY\nIANIGLA"]
colors = plt.cm.get_cmap("tab10",len(names))(np.linspace(0,1,len(names)))
ax = ax.ravel()
# ax[4].plot([],[],ls=":",color="k",label="Linear\nRegression")
ax[5].plot([],[],ls=":",color="r",label="$y\sim x$")

for i in range(len(ax)-1):
    ax[i].scatter(sl_ianigla2.loc[tile_props.index],
                  var[i],edgecolor="k",alpha=0.6,color=colors[i],zorder=3,
                  lw=0.4)
    # m = linregress(sl_ianigla2.loc[tile_props.index],var[i])
    # ax[i].plot(np.arange(800,4.5e3),
    #            np.arange(800,4.5e3)*m.slope+m.intercept,ls=":",color="k")
    # t=ax[i].text(x=0.02,y=0.75,
    #              s="y~"+"{:.2f}".format(m.slope)+"x"+"{:.2f}".format(m.intercept)+"\n$R^2$: "+"{:.1%}".format(m.rvalue**2),
    #              transform=ax[i].transAxes)
    ax[i].plot([1e3,7e3],[1e3,7e3],color="red",ls=":")
    ax[i].set_xlim(800,6.5e3)
    ax[i].set_ylim(800,6.5e3)
    # ax[i].set_title(names[i],loc="left")
    ax[i].grid(True,ls=":")
    ax[5].scatter([],[],label=names[i],color=colors[i],edgecolor="k",lw=0.4)


# for i in [1,2,4]:
#     ax[i].set_yticklabels("")
ax[5].axis("off")
# ax[4].legend(frameon=False,loc=(1.85,0.77))
ax[5].legend(frameon=False,loc=(-0.05,0.12),ncol=1)
# ax[5].scatter(sl_ianigla.loc[tile_props.index],dH,edgecolor="k",alpha=0.6)
# ax[5].set_title("dH",loc="left")

ax[4].set_xlabel("Snow Limit IANIGLA (m)")

plt.savefig("plots/maipomanzano/scatterplots_snowlimitS.pdf",dpi=150,bbox_inches="tight")
#%%























    
    