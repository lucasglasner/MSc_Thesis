#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 18:48:01 2021

@author: lucas

# =============================================================================
# ROS in Maipo Basin
# =============================================================================

"""


from functions import seasonal_decompose, sliding_interval_filter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime as dt
import scipy.stats as st
import xarray as xr
from scipy.stats import linregress
from tqdm import trange
from glob import glob
import gc
import geopandas as gpd
import scipy.stats as st
import sys
# %%
sys.path.append('functions.py')
# %%
# =============================================================================
# load data
# =============================================================================

hypso = pd.read_csv(
    'datos/topography/basins/hypso/RioMaipoEnElManzano_hypso.csv')

isotermas0 = pd.read_csv(
    "datos/isotermas0_maipomanzano.csv", index_col=0).dropna(how="all")
isotermas0.index = pd.to_datetime(isotermas0.index)
isotermas0 = isotermas0.resample("d").fillna(method="ffill")

int_func = interp1d(hypso.height, hypso.fArea)
pluv_area = int_func(isotermas0['STODOMINGO']-300)
pluv_area = pd.Series(pluv_area, index=isotermas0.index)


snowlimits = pd.read_csv(
    "datos/snowlimits_maipomanzano.csv", index_col=0).dropna(how="all")
snowlimits.index = pd.to_datetime(snowlimits.index)


SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0).dropna(how='all')
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA['IANIGLA']/100


ROS_AREA = SCA-(1-pluv_area)

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv", index_col=0)
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met.drop("date", axis=1, inplace=True)

pr_mm = pd.read_csv('datos/estaciones/pr_laobra.csv', dtype=str)
pr_mm.index = pd.to_datetime(
    pr_mm.iloc[:, 0]+"-"+pr_mm.iloc[:, 1]+"-"+pr_mm.iloc[:, 2])
pr_mm = pr_mm.iloc[:, 3]
pr_mm = pd.to_numeric(pr_mm)

ROS_AREA = np.clip(ROS_AREA.reindex(pr_mm.index), 1, 0)
ROS_AREA = ROS_AREA.where(pr_mm > 3)


qinst_mm = pd.read_csv(
    "datos/estaciones/qinst_RioMaipoEnElManzano.csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)

mm_centroid = (-70.08756, -33.70859)
# %%
# =============================================================================
# Fill MODIS data with linear regression against ianigla
# =============================================================================
times = pd.date_range("2000-01-01", "2021-12-31", freq="d")
for i in trange(3):
    y = snowlimits.iloc[:, i].dropna()
    x = snowlimits["IANIGLA"].dropna().reindex(y.index)
    m = linregress(x, y)
    for j in range(len(snowlimits)):
        if np.isnan(snowlimits.iloc[j, i]):
            snowlimits.iloc[j, i] = m.slope * \
                snowlimits["IANIGLA"].values[j]+m.intercept
gc.collect()

# %% Rasters and ROS
# =============================================================================
# CORTES + PR_CR2MET + ERA5_H0
# =============================================================================
how = "whole"
if how == "whole":
    paths = "datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES*"
    paths = glob(paths)
    SWE = xr.open_mfdataset(paths, chunks='auto').SWE
    paths = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_*.nc'
    dSWE = xr.open_mfdataset(paths).SWE
    dSWE = dSWE.reindex({'time': SWE.time.to_series().index})
    paths = "datos/era5/H0_ERA5_*.nc"
    H0 = xr.open_mfdataset(paths, chunks='auto').deg0l
    H0 = H0.reindex({"time": SWE.time.to_series().index,
                     'lat': SWE.lat,
                     'lon': SWE.lon},
                    method='nearest')

    paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_pr_1979-2020.nc"
    PR = xr.open_dataset(paths, chunks='auto').pr
    PR = PR.reindex({"time": SWE.time.to_series().index},
                    method='nearest')

    ROS_CCE = np.where((SWE > 10) & (H0 > 300) & (PR > 3) & (dSWE < 0),
                       True, False)
    ROS = np.empty(ROS_CCE.shape[0])
    for i in range(ROS_CCE.shape[0]):
        ROS[i] = ROS_CCE[i, :, :].sum()

    ROS_CCE = pd.Series(ROS, index=SWE.time.values)
    ROS_CCE = ROS_CCE/191
    del SWE, H0, PR, paths
elif how == "centroid":
    paths = "datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES*"
    paths = glob(paths)
    SWE = xr.open_mfdataset(paths).SWE
    SWE = SWE.sel(lat=mm_centroid[1], lon=mm_centroid[0],
                  method="nearest").to_series()

    paths = "datos/era5/H0_ERA5_*.nc"
    H0 = xr.open_mfdataset(paths, chunks='auto').deg0l
    H0 = H0.sel(lat=mm_centroid[1], lon=mm_centroid[0],
                method="nearest").to_series()
    H0 = H0.reindex(SWE.index, method='nearest')

    paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_pr_1979-2020.nc"
    PR = xr.open_dataset(paths).pr
    PR = PR.sel(lat=mm_centroid[1], lon=mm_centroid[0],
                method="nearest")
    PR = PR.reindex({"time": SWE.index}).to_series()

    ROS_CCE = np.where((SWE > 10) & (H0 > 300) & (PR > 3),
                       True, False)
    ROS_CCE = pd.Series(ROS_CCE,
                        index=SWE.index).dropna()
    del H0, SWE, PR, paths
else:
    raise Exception("'how' method not valid")

gc.collect()


# # =============================================================================
# # ERA5LAND
# # =============================================================================
# how = "whole"
# if how == "whole":
#     paths = sorted(glob("datos/era5land/RioMaipoEnElManzano/*.nc"))
#     paths = paths[1:4:2]
#     ERA5LAND = xr.open_mfdataset(paths)

#     ROS_ERA5LAND = np.where((ERA5LAND.tp > 3/1e3) & (ERA5LAND.sd > 10/1e3),
#                             True, False)
#     ROS = np.empty(ROS_ERA5LAND.shape[0])
#     for i in range(ROS_ERA5LAND.shape[0]):
#         ROS[i] = ROS_ERA5LAND[i, :, :].sum()

#     ROS_ERA5LAND = pd.Series(ROS, index=ERA5LAND.time.values)
#     ROS_ERA5LAND = ROS_ERA5LAND/49
#     del ERA5LAND, paths
# elif how == "centroid":
#     paths = sorted(glob("datos/era5land/RioMaipoEnElManzano/*.nc"))
#     paths = paths[1:4:2]
#     ERA5LAND = xr.open_mfdataset(paths)
#     ERA5LAND = ERA5LAND.sel(lat=mm_centroid[1], lon=mm_centroid[0],
#                             method="nearest")
#     ERA5LAND = ERA5LAND.to_dataframe().drop(["lon", "lat"], axis=1).dropna()
#     ERA5LAND = ERA5LAND.resample("d").sum()
#     ROS_ERA5LAND = np.where((ERA5LAND.tp > 3*1e3) & (ERA5LAND.sd > 10*1e3),
#                             True, False)
#     ROS_ERA5LAND = pd.Series(ROS_ERA5LAND, index=ERA5LAND.index)
#     del ERA5LAND, paths
# else:
#     raise Exception("'how' method not valid")
# gc.collect()
# # %%
# # =============================================================================
# # CORTES+CR2MET
# # =============================================================================
# how = "whole"
# if how == "whole":
#     paths = "datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES*"
#     paths = glob(paths)
#     SWE = xr.open_mfdataset(paths).SWE

#     paths = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_*.nc'
#     dSWE = xr.open_mfdataset(paths).SWE
#     dSWE = dSWE.reindex({'time': SWE.time.to_series().index})
#     paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_t2m_1979-2020.nc"
#     T2M = xr.open_dataset(paths).t2m
#     T2M = T2M.reindex({"time": SWE.time.to_series().index})

#     paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_pr_1979-2020.nc"
#     PR = xr.open_dataset(paths).pr
#     PR = PR.reindex({"time": SWE.time.to_series().index})

#     ROS_CORTESCR2MET = np.where((SWE > 10) & (T2M > 0) & (PR > 3) & (dSWE < 0),
#                                 True, False)
#     ROS = np.empty(ROS_CORTESCR2MET.shape[0])
#     for i in range(ROS_CORTESCR2MET.shape[0]):
#         ROS[i] = ROS_CORTESCR2MET[i, :, :].sum()

#     ROS_CORTESCR2MET = pd.Series(ROS, index=SWE.time.values)
#     ROS_CORTESCR2MET = ROS_CORTESCR2MET/191
#     del SWE, T2M, PR, paths
# elif how == "centroid":
#     paths = "datos/ANDES_SWE_Cortes/regrid_cr2met/RioMaipoEnElManzano/ANDES*"
#     paths = glob(paths)
#     SWE = xr.open_mfdataset(paths).SWE
#     SWE = SWE.sel(lat=mm_centroid[1], lon=mm_centroid[0],
#                   method="nearest").to_series()

#     paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_t2m_1979-2020.nc"
#     T2M = xr.open_dataset(paths).t2m
#     T2M = T2M.sel(lat=mm_centroid[1], lon=mm_centroid[0],
#                   method="nearest")
#     T2M = T2M.reindex({"time": SWE.index}).to_series()

#     paths = "datos/cr2met/RioMaipoEnElManzano_CR2MET_pr_1979-2020.nc"
#     PR = xr.open_dataset(paths).pr
#     PR = PR.sel(lat=mm_centroid[1], lon=mm_centroid[0],
#                 method="nearest")
#     PR = PR.reindex({"time": SWE.index}).to_series()

#     ROS_CORTESCR2MET = np.where((SWE > 10) & (T2M > 0) & (PR > 3),
#                                 True, False)
#     ROS_CORTESCR2MET = pd.Series(ROS_CORTESCR2MET,
#                                  index=SWE.index).dropna()
#     del T2M, SWE, PR, paths
# else:
#     raise Exception("'how' method not valid")

# gc.collect()
# %%
# =============================================================================
# COMPUTE ROS BASED ON SNOWLIMIT AND FREEZING LEVEL METHOD.
# =============================================================================
times = pd.date_range("1985-01-01", "2014-12-31", freq="d")
SL_cond = snowlimits.reindex(pr_cr2met.index)[
    pr_cr2met.values > 3].dropna(how="all").reindex(times)
FL_cond = isotermas0.reindex(pr_cr2met.index)[
    pr_cr2met.values > 3].reindex(times)-300

ROS1 = {sl: {fl: None for fl in FL_cond.columns} for sl in SL_cond.columns}
for SL in SL_cond.columns:
    for FL in FL_cond.columns:
        interp = interp1d(hypso['height'], hypso['Area_km2'],
                          fill_value='extrapolate')
        snow_area = interp(SL_cond[SL])
        liquid_area = interp(FL_cond[FL])
        ROS_area = (liquid_area-snow_area)/hypso['Area_km2'].iloc[-1]
        ROS1[SL][FL] = pd.Series(ROS_area,
                                 index=SL_cond[SL].index)
        ROS1[SL][FL].name = SL+" - "+FL
ROS1 = pd.concat([pd.concat(list(ROS1.values())[i].values(), axis=1)
                  for i in range(len(SL_cond.columns))], axis=1)

pairs = ["MODIS_H50 - STODOMINGO",
         # "CORTES - CR2MET",
         # "CORTES_H50 - CR2MET_H50_MM",
         "CORTES - CR2MET - ERA5"]

ROS1["CORTES - CR2MET - ERA5"] = ROS_CCE.reindex(times)
# ROS1["CORTES - CR2MET"] = ROS_CORTESCR2MET.reindex(times)
ROS1 = ROS1[pairs]

ROS = ROS1 > 0
ROS["CORTES - CR2MET - ERA5"] = ROS_CCE.reindex(times) > 0.1
# ROS["CORTES - CR2MET"] = ROS_CORTESCR2MET.reindex(times) > 0.1
ROS = ROS[pairs]

# ROS["ERA5LAND"] = ROS_ERA5LAND.reindex(times) > 0.2

nanmask = np.isnan(ROS1['MODIS_H50 - STODOMINGO'])

ROS['MODIS_H50 - STODOMINGO'] = ROS['MODIS_H50 - STODOMINGO'].where(~nanmask)
meanROS = ROS.groupby([ROS.index.year, ROS.index.month]).sum()
meanROS = meanROS.unstack().mean(axis=0).unstack().T
meanROS['MODIS_H50 - STODOMINGO'] = ROS['MODIS_H50 - STODOMINGO'].dropna().groupby([ROS['MODIS_H50 - STODOMINGO'].dropna().index.year,
                                                                                    ROS['MODIS_H50 - STODOMINGO'].dropna().index.month]).sum().unstack().mean(axis=0)
meanROS = meanROS.iloc[:, [-1, 0, 1]]
for pair in pairs:
    ROS[pair] = ROS[pair].map(lambda x: False if np.isnan(x) else bool(x))

# meanROS = ROS.reindex(ROS1['MODIS_H50 - STODOMINGO'].dropna().index)
ROS1 = ROS1.where((~np.isnan(ROS1)) & (ROS1 >= 0), np.nan)

gc.collect()
# %%
events = []
for i in range(len(ROS.columns)):
    y = ROS.iloc[:, i]
    x = y.groupby([y, (y != y.shift()).cumsum()]).size()
    events.append(x.unstack().iloc[1, :].dropna().reset_index(drop=True))


durations = [(events[j].groupby(events[j]).sum()/events[j].shape).reindex(np.arange(20))
             for j in range(len(ROS.columns))]
durations = pd.concat(durations, axis=1)
durations.columns = ROS.columns
# %%
fig, ax = plt.subplots(2, 2, figsize=(14, 7))
fig.tight_layout(pad=1.5)
plt.rc('font', size=18)
ax = ax.ravel()
year = ROS.resample("y").sum().applymap(lambda x: np.nan if x == 0 else x)
# month = ROS.resample("m").sum().applymap(lambda x: np.nan if x == 0 else x)
for pair in pairs:
    var = year[pair].dropna()
    m = st.linregress(range(len(var)), var)
    trend = "{:.2f}".format(m.slope)
    pvalue = "{:.2f}".format(m.pvalue)
    ax[1].plot(var, label="Trend: "+trend+" days/year",
               alpha=0.8)
ax[1].plot(year[pairs[0]].dropna(), color='tab:blue', lw=2)
ax[1].legend(loc=(0, 1), frameon=False, fontsize=16)
# ax[1].set_xticklabels()
# ax[1].sharex(ax[3])
ax[1].tick_params(axis='x', rotation=45)
meanROS.plot.bar(ax=ax[0], alpha=0.8, edgecolor='k')


for i in meanROS.columns:
    ax[0].plot(np.arange(0, 12, 1), meanROS[i], ls=":")
#     var = ROS.groupby(ROS.index.dayofyear).sum()[pair]
#     var = seasonal_decompose(var,365,6,1)[0]
#     ax[0].plot(np.linspace(1,12,len(var)),var)
ax[0].legend(frameon=False, loc=(0, 1), fontsize=16)
ax[0].set_ylim(0, 7)
ax[0].set_xticks(np.arange(0, 12, 1))
ax[0].set_xticklabels(["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
                       "AGO", "SEP", "OCT", "NOV", "DIC"])
ax[0].tick_params(axis='x', rotation=45)
ax[0].grid(ls=":", axis="x")
ax[0].set_ylabel("Monthly mean\nROS Events)")
ax[1].set_ylabel("N° ROS days")

ax[2].scatter(ROS1.iloc[:, 0], ROS1.iloc[:, 1], color='tab:orange', ec='k',
              alpha=0.5)
ax[2].scatter(ROS1.iloc[:, 0], ROS1.iloc[:, 2], color='tab:green', ec='k',
              alpha=0.5)
ax[2].set_ylim(0, 1)
ax[2].set_xlim(0, 1)
ax[2].plot([0, 1], [0, 1], 'k')
ax[2].set_ylabel('Model ROS (REANALYSIS)\nFractional Area (-)')
ax[2].set_xlabel('Observed ROS (MODIS/STODOMINGO)\n Fractional Area (-)')


durations.plot.bar(ax=ax[3], legend=False, ec='k')
ax[3].set_xlim(0, 15)
ax[3].set_xlabel('ROS Events Duration (days)')
ax[3].set_xticks(np.arange(1, 15))
ax[3].tick_params(axis='x', rotation=0)
ax[3].set_ylabel('N°Events\nover total Events')
ax[3].grid(axis='x')
ax[3].set_ylim(0, 1)

# ax[2].plot(ROS)
# plt.savefig("plots/maipomanzano/ROS_maipo.pdf", dpi=150, bbox_inches="tight")
#
# %%
plt.rc('font', size=10)
# =============================================================================
# ros time series simple
# =============================================================================

ROS_mm = xr.open_mfdataset(glob('datos/ROS/CORTES_CR2MET_ERA5/ROS_*.nc'))


fig, ax = plt.subplots(1, 2, sharex=False, sharey=True, figsize=(10, 3))
ax = ax.ravel()

# ax[0].plot(ROS1.iloc[:, 0])
ax[0].plot(ROS1.iloc[:, 2], alpha=0.5)

yr = '2013'

ROS1[yr].iloc[:, 2].plot(ax=ax[1], alpha=0.5)
ax[1].axhline(0.15, c='k', ls=':')
# ax[1].plot(ROS1.iloc[:,0]["2013"])
# ax[1].plot(ROS1.iloc[:, 2][yr], alpha=0.5)
#
# %%
# =============================================================================
# ROS consequences on runoff
# =============================================================================


# %%
# minims = ROS11.min()
# dates  = []
# for i in range(minims.shape[0]):
#     v=np.where(ROS11.iloc[:,i]==minims.values[i])
#     dates.append(ROS11.index.values[v[0][0]])

# dates=pd.to_datetime(dates)

# #%%
# fig,ax=plt.subplots(2,1)
# ts=[]
# p=[]
# for date in dates[:1]:
#     d1 = date+dt.timedelta(days=5)
#     d2 = date-dt.timedelta(days=5)
#     t  = qinst_mm[d2.strftime("%Y-%m-%d"):d1.strftime("%Y-%m-%d")]
#     pr = pr_cr2met[d2.strftime("%Y-%m-%d"):d1.strftime("%Y-%m-%d")]
#     ind = t.index-date
#     ts.append(pd.Series(t.values,index=ind))
#     p.append(pd.Series(pr.values.squeeze(),index=pr.index-date))
# for t in ts:
#     ax[0].plot(t)
# for pr in p:
#     ax[1].bar(np.arange(len(pr)),pr,alpha=0.5)
# %%
# =============================================================================
# load era5land data
# =============================================================================
# pr_era5land  = xr.open_dataset("datos/era5land/RioMaipoEnElManzano/total_precipitation.nc",chunks="auto").tp*1e3
# t2m_era5land = xr.open_dataset("datos/era5land/RioMaipoEnElManzano/2m_temperature.nc",chunks="auto").t2m-273.15
# swe_era5land = xr.open_dataset("datos/era5land/RioMaipoEnElManzano/snow_depth_water_equivalent.nc",chunks="auto").sd*1e3


# cond1 = pr_era5land>1
# cond2 = t2m_era5land<0
# cond3 = swe_era5land>10
# ros_era5land = cond1 & cond2 & cond3
# ros_era5land = ros_era5land.resample({"time":"d"}).sum()//24

# ros_era5land = ros_era5land.groupby("time.year").sum()
# try:
#     pr_era5land = pd.read_csv("datos/era5land/RioMaipoEnElManzano/total_precipitation_basinmean.csv",index_col=0)
#     pr_era5land.index = pd.to_datetime(pr_era5land.index)
# except:
#
#     pr_era5land = pr_era5land.mean(dim="lat").mean(dim="lon").to_series()*1e3
#     pr_era5land.to_csv("datos/era5land/RioMaipoEnElManzano/total_precipitation_basinmean.csv")

# %%

# #%%
# # =============================================================================
# # Cargar datos
# # =============================================================================
# cuenca           = "RioMauleEnArmerillo"
# curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
# curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)

# snow_limit       =  pd.read_csv("datos/ianigla/"+cuenca+"_SCA_s_comp.filtro_MA.3días.csv",index_col=0)/100
# snow_limit.index =  pd.to_datetime(snow_limit.index)

# isoterma0        =  pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)-300
# isoterma0.index  = pd.to_datetime(isoterma0.index)
# isoterma0        = isoterma0.where(isoterma0<6000).dropna()["H0_StoDomingo"]
# isoterma0 = pd.Series(isoterma0.groupby([isoterma0.index.year,
#                                           isoterma0.index.month,
#                                           isoterma0.index.day]).mean().ravel(),
#                       index=np.unique(pd.to_datetime(isoterma0.index.date)),
#                       name = "H0_StoDomingo")
# isoterma0        = isoterma0.resample("d").interpolate("linear")

# pr_quintanormal  =  pd.read_csv("datos/estaciones/pr_quintanormal.csv").applymap(lambda x: str(x))
# pr_quintanormal["fecha"] = pr_quintanormal["agno"]+"-"+pr_quintanormal[" mes"]+"-"+pr_quintanormal[" dia"]
# pr_quintanormal.index = pd.to_datetime(pr_quintanormal["fecha"])
# pr_quintanormal.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
# pr_quintanormal  = pd.to_numeric(pr_quintanormal[" valor"])
# pr_quintanormal.name = "pr"

# q_maipomanzano   = pd.read_csv("datos/estaciones/q_"+cuenca+".csv").applymap(lambda x: str(x))
# q_maipomanzano["fecha"] = q_maipomanzano["agno"]+"-"+q_maipomanzano[" mes"]+"-"+q_maipomanzano[" dia"]
# q_maipomanzano.index = pd.to_datetime(q_maipomanzano["fecha"])
# q_maipomanzano.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
# q_maipomanzano  = pd.to_numeric(q_maipomanzano[" valor"])
# q_maipomanzano.name = "q"
# #%%
# # =============================================================================
# # Calcular linea de nieves
# # =============================================================================

# interp = interp1d(1-curva_hipso["fArea"],curva_hipso["height"])
# snow_limit["SL"] = list(map(lambda x: interp(x).item(),snow_limit["SCA(%)"]))

# #%%
# daterange = pd.date_range("1950-01-01","2021-08-22", freq="d")
# data = pd.concat((pr_quintanormal.reindex(daterange),
#                   isoterma0.reindex(daterange),
#                   snow_limit.reindex(daterange)),axis=1)
# data = data.where(data["pr"]>1).dropna() #Dias con lluvia
# data = data.where(data["SL"]<data["H0_StoDomingo"]-300).dropna()
# data["q"] = q_maipomanzano.reindex(data.index)

# #%%
# fig,ax = plt.subplots(2,2,figsize=(10,5))
# fig.tight_layout(pad=2)
# fig.text(0.75,1,"ROS 2000->2017", ha="center", va="center",fontsize=20)
# ax = ax.ravel()

# ax[0].hist(isoterma0,bins="auto",density=True,alpha=0.7,color="goldenrod",
#            label="sin_pr")
# ax[0].hist(isoterma0.reindex(pr_quintanormal.index)[pr_quintanormal>1],bins="auto",
#            density=True,alpha=0.7,color="royalblue", label="con_pr")
# ax[0].hist(data["H0_StoDomingo"], alpha=0.7, color="skyblue", density=True,
#            label="ROS")
# ax[0].legend(loc=(0,1.05), frameon=False, ncol=3)
# ax[0].set_xlabel("isoterma0_maipomanzano (m)")

# ax[1].hist(snow_limit["SL"],bins="auto", color="goldenrod",density=True , alpha=0.7)
# ax[1].hist(snow_limit["SL"].reindex(pr_quintanormal.index)[pr_quintanormal>1],
#            bins="auto", density=True, alpha=0.7, color="royalblue" )
# ax[1].hist(data["SL"], bins="auto" , density=True, alpha=0.7, color="skyblue")
# ax[1].set_xlabel("snow_limit_maipomanzano (m)")


# ax[2].hist(pr_quintanormal[pr_quintanormal>1], bins="auto", density=True, alpha=0.7,
#            color="royalblue" )
# ax[2].hist(data["pr"], bins="auto" , density=True, alpha=0.7, color="skyblue")
# ax[2].set_xlabel("pr_quintanormal (mm)")


# ax[3].hist(q_maipomanzano, bins="auto", density=True, alpha=0.7,
#            color="goldenrod" )
# ax[3].hist(q_maipomanzano.reindex(pr_quintanormal.index)[pr_quintanormal>1], bins="auto", density=True, alpha=0.7,
#            color="royalblue")
# ax[3].hist(data["q"], bins="auto" , density=True, alpha=0.7, color="skyblue")
# ax[3].set_xlabel("Q_maipomanzano (m3/s)")


# ax[3].axvline(st.expon.ppf(1-1/25,*st.expon.fit(q_maipomanzano)), color="k", label="25, 50, 100, 500\nreturn periods", ls=":")
# ax[3].axvline(st.expon.ppf(1-1/50,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
# ax[3].axvline(st.expon.ppf(1-1/100,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
# ax[3].axvline(st.expon.ppf(1-1/500,*st.expon.fit(q_maipomanzano)), color="k", ls=":")
# ax[3].legend()


# plt.savefig("plots/maipomanzano/ROS_maipo.pdf",dpi=150,bbox_inches="tight")

# #%%

# fig,ax = plt.subplots(2,1,sharex=True,figsize=(10,5))
# ax = ax.ravel()
# date = "2009-09-06"
# step = 40
# sliced = slice(pd.to_datetime(date)-dt.timedelta(days=step),
#                pd.to_datetime(date)+dt.timedelta(days=step))
# # xticks = pd.date_range(sliced.start, sliced.stop, freq='d')
# snow_limit["SL"].loc[sliced].plot(ax=ax[0],color="royalblue",label="snow_limit")
# isoterma0.loc[sliced].plot(ax=ax[0],color="dimgrey",label="isoterma 0")
# ax[1].bar(pr_quintanormal.loc[sliced].index,pr_quintanormal.loc[sliced], color="teal")
# ax[1].set_ylabel("pr (mm)")
# ax[0].set_ylabel("z (m.s.n.m)")
# ax2 = ax[1].twinx()
# ax2.set_ylabel("Q (m3/s)")
# q_maipomanzano.loc[sliced].rolling(window=1).mean().plot(ax=ax2, color="blueviolet", lw=2)
# ax[0].legend(loc=(0,1.1), ncol=2, frameon=False)
# # ax[1].legend()
# # ax[2].legend()
# for axis in ax:
#     axis.set_xlabel("")
# # ax[0].set_xticklabels([x.strftime('%m-%d') for x in xticks]);

# plt.savefig("plots/maipomanzano/ROS"+date+"_maipo.pdf",dpi=150,bbox_inches="tight")
