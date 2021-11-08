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
import xarray as xr
import geopandas as gpd
from tqdm import trange
from scipy.stats import linregress
#%%
# =============================================================================
# Load data
# =============================================================================

# =============================================================================
# Basin hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)
basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")

# =============================================================================
# #load data from sto domingo radiosonde
# =============================================================================
H0_stodomingo = pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)
H0_stodomingo.index = pd.to_datetime(H0_stodomingo.index)
H0_stodomingo = H0_stodomingo.where(H0_stodomingo<7e3).dropna()["H0_StoDomingo"]

# =============================================================================
# #load data from AMDAR database
# =============================================================================
H0_amdar = pd.read_csv("datos/amdar/isoterma0SCEL20172019.csv",index_col=0)
func = lambda x: datetime.fromordinal(int(x)) + timedelta(days=x%1) - timedelta(days = 366)
H0_amdar.index = H0_amdar["timeregh"].map(func)
H0_amdar = H0_amdar["zt0"]
H0_amdar.index = pd.date_range("2017-01-01T01:00:00","2020-01-01",freq="h").dropna()
# H0_amdar.index = H0_amdar.index.map(lambda x: x+timedelta(hours=4))

# =============================================================================
# #load precipitation data
# =============================================================================
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

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv")
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met = pr_cr2met["5710001"]

precip_days = pr_cr2met>5
precip_days = precip_days[precip_days].index
# =============================================================================
# 
# #load cr2met terrain and temperature data
# =============================================================================
dem        = xr.open_dataset("datos/topography/basins/RioMaipoEnElManzano_CR2MET.nc").Band1
t2m_cr2met = xr.open_dataset("datos/cr2met/CR2MET_t2m_v2.0_day_1979_2020_005deg_RioMaipoEnElManzano.nc",
                             chunks="auto").t2m


# =============================================================================
# load era5 zero degree level
# =============================================================================

H0_ERA5 = pd.read_csv("datos/era5/H0_stodomingo.csv")
H0_ERA5 = H0_ERA5[H0_ERA5.expver==1].dropna()
H0_ERA5.index = pd.to_datetime(H0_ERA5.time)
H0_ERA5.drop(["time","lon","lat","expver"],axis=1,inplace=True)
H0_ERA5 = H0_ERA5.deg0l

#%%

# =============================================================================
# Build freezing level height for maipo en el manzano basin based upon a simple
# linear regression model between elevation and temperature of each day
# =============================================================================
x = dem.values.ravel()
linmodels = []
for i,date in enumerate(t2m_cr2met.time):
    y=t2m_cr2met.sel(time=date).values.ravel()
    mask = ~np.isnan(x) & ~np.isnan(y)
    m=linregress(x[mask],y[mask])
    linmodels.append(m)
    
# pearsons       = [m.rvalue for m in linmodels]
H0_cr2met_OLR  = [-m.intercept/m.slope for m in linmodels] 
H0_cr2met_OLR  = pd.Series(H0_cr2met_OLR,index=t2m_cr2met.time.values)

#%%
# =============================================================================
# Define a function for calculating the number of pixels with T<0 below the
# freezing level and the number of pixels with T>0 above the freezing level
# =============================================================================

def level_error(timeseries,raster,dem,field_limit=0):
    timeseries = timeseries.reindex(raster.time.values).dropna()
    raster = raster.sel(time=timeseries.index)  
    tpix   = np.count_nonzero(~np.isnan(dem))
    metric = np.empty(len(timeseries))
    for i,date in enumerate(timeseries.index):
        height    = timeseries.loc[date]
        mask      = dem<height
        below     = raster[i,:,:].where(mask)<field_limit
        above     = raster[i,:,:].where(~mask)>field_limit
        metric[i] = (np.count_nonzero(below)+np.count_nonzero(above))/(tpix)
    return pd.Series(metric,index=timeseries.index)
        
#%%

# =============================================================================
# compute the metric for each method/dataset
# =============================================================================

metric_stodomingo = level_error(H0_stodomingo,t2m_cr2met,dem)
metric_amdar      = level_error(H0_amdar,t2m_cr2met,dem)
metric_cr2met_OLR = level_error(H0_cr2met_OLR,t2m_cr2met,dem)
metric_era5       = level_error(H0_ERA5,t2m_cr2met,dem)


#%%

fig,ax = plt.subplots(1,3,sharex=True,sharey=True,figsize=(12,3))
# plt.scatter(H0_amdar,H0_cr2met_OLR.reindex(H0_amdar.index))

mask = pd.Series(H0_stodomingo.index.date).map(lambda x: x in precip_days.date).values
for i,var in enumerate([H0_amdar,H0_cr2met_OLR,H0_ERA5]):
    ax[i].scatter(H0_stodomingo,var.reindex(H0_stodomingo.index),alpha=0.7)
    ax[i].scatter(H0_stodomingo[mask],var.reindex(H0_stodomingo.index)[mask],alpha=0.7)
    ax[i].plot([0,8e3],[0,8e3],"k--")
    ax[i].set_xlim(0,8e3)
    ax[i].set_ylim(0,8e3)
    ax[i].set_yticks(ax[i].get_xticks())
    ax[i].grid(True,ls=":")
    
    
ax[0].scatter([],[],color="tab:orange",label="Pr>5mm")
ax[0].legend(frameon=False)
ax[2].plot([],[],color="k",ls="--",label="y~x")
ax[2].legend(frameon=False)

ax[1].set_xlabel("H0_amdar (m)",fontsize=12)
ax[0].set_title("H0_StoDomingo (m)")
ax[1].set_title("H0_CR2MET_LR (m)")
ax[2].set_title("H0_ERA5 (m)")

plt.savefig("plots/maipomanzano/datasetcomparison/isotherm0scatters.pdf",dpi=150,bbox_inches="tight")

# # =============================================================================
# # Read file with freezing temperature by band or make it.
# # =============================================================================
# # precip_days=precip_days[:100]
# try: 
#     flevel_bands = pd.read_csv("datos/cr2met/freezinglevel_t2m_bdsands.csv",
#                              index_col=0)
#     flevel_bands.columns = pd.to_datetime(flevel_bands.columns)
# except:
#     # =========================================================================
#     # Build elevation band masks
#     # =========================================================================
#     dz              = 200
#     elevation_bands = np.arange(dem.min(),dem.max()+dz,dz)
    
#     masks = []
#     for j in range(len(elevation_bands)-1):
#             z0   = elevation_bands[j]
#             z1   = elevation_bands[j+1]
#             mask = (dem.where((dem>z0) & (dem<z1))>0).values
#             masks.append(mask)
    
#     elevation_bands = elevation_bands[:-1]
#     # =========================================================================
#     # Apply band mask, and compute % of pixels with temperature < 0°C by band
#     # =========================================================================
    
#     flevel_bands      = np.empty((len(elevation_bands),len(precip_days)))
#     flevel_bands      = pd.DataFrame(flevel_bands,
#                                      index = elevation_bands,
#                                      columns=precip_days)
#     for i in trange(len(precip_days)):
#         tile_date = precip_days[i]
#         for j in range(len(masks)):
#             band = masks[j] 
#             tile_band     = t2m_cr2met.sel(time=tile_date).where(band)
#             freezing_band = np.where((tile_band<0),1,0)
#             flevel_bands.iloc[j,i] = freezing_band.sum()/np.count_nonzero(band)
#     flevel_bands.to_csv("datos/cr2met/freezinglevel_t2m_bands.csv")




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
mask = pr_cr2met.reindex(hamd.index)>5

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
plt.savefig("plots/maipomanzano/datasetcomparison/isotherm0_study.pdf",dpi=150,bbox_inches="tight")


