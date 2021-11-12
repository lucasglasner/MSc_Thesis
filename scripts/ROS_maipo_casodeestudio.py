#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 16:16:54 2021

@author: lucas

# =============================================================================
# ROS in Maipo Basin: case of study
# =============================================================================

"""

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy.interpolate import interp1d
import datetime
from scipy.ndimage.filters import minimum_filter1d, generic_filter
from scipy.ndimage.measurements import label
from scipy.signal import argrelextrema

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError

    if x.size < window_len:
        raise ValueError


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def minimum_filter(ts, **kwargs):
    """Return stationary base flow
    
    The base flow is set to the minimum observed flow.
    
    :param ts: 
    :return: 
    """
    minimum = min(ts)
    out_values = minimum * np.ones(len(ts))
    baseflow = pd.Series(data=out_values, index=ts.index)
    quickflow = ts - baseflow
    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'
    return baseflow, quickflow


def fixed_interval_filter(ts, size):
    """USGS HYSEP fixed interval method
    
    The USGS HYSEP fixed interval method as described in `Sloto & Crouse, 1996`_.
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.
    
    :param size: 
    :param ts: 
    :return: 
    """
    intervals = np.arange(len(ts)) // size
    baseflow = pd.Series(data=ts.groupby(intervals).transform('min'), index=ts.index)
    quickflow = ts - baseflow

    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


def sliding_interval_filter(ts, size):
    """USGS HYSEP sliding interval method
    
        The USGS HYSEP sliding interval method as described in `Sloto & Crouse, 1996`_.
        
        The flow series is filter with scipy.ndimage.genericfilter1D using np.nanmin function
        over a window of size `size`
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.
    
    :param size: 
    :param ts: 
    :return: 
    """
    # TODO ckeck the presence of nodata
    if (ts.isnull()).any():
        blocks, nfeatures = label(~ts.isnull())
        block_list = [ts[blocks == i] for i in range(1, nfeatures + 1)]
        na_df = ts[blocks == 0]
        block_bf = [pd.Series(data=minimum_filter1d(block, size, mode='reflect'), index=block.index) for block in
                    block_list]
        baseflow = pd.concat(block_bf + [na_df], axis=0)
        baseflow.sort_index(inplace=True)
    else:
        baseflow = pd.Series(data=minimum_filter1d(ts, size, mode='reflect'), index=ts.index)

    quickflow = ts - baseflow

    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


def _local_minimum(window):
    win_center_ix = len(window) / 2
    win_center_val = window[win_center_ix]
    win_minimum = np.min(window)
    if win_center_val == win_minimum:
        return win_center_val
    else:
        return np.nan

def local_minimum_filter(ts, size):
    
    
    """USGS HYSEP local minimum method
    
        The USGS HYSEP local minimum method as described in `Sloto & Crouse, 1996`_.
    
    .. _Slot & Crouse, 1996:
        Sloto, Ronald A., and Michele Y. Crouse. “HYSEP: A Computer Program for Streamflow Hydrograph Separation and 
        Analysis.” USGS Numbered Series. Water-Resources Investigations Report. Geological Survey (U.S.), 1996. 
        http://pubs.er.usgs.gov/publication/wri964040.
    
    :param size: 
    :param ts: 
    :return: 
    """

    origin = int(size) / 2
    baseflow_min = pd.Series(generic_filter(ts, _local_minimum, footprint=np.ones(size)), index=ts.index)
    baseflow = baseflow_min.interpolate(method='linear')
    # interpolation between values may lead to baseflow > streamflow
    errors = (baseflow > ts)
    while errors.any():
        print('hello world')
        error_labelled, n_features = label(errors)
        error_blocks = [ts[error_labelled == i] for i in range(1, n_features + 1)]
        error_local_min = [argrelextrema(e.values, np.less)[0] for e in error_blocks]
        print(error_local_min)
        break
    quickflow = ts - baseflow
    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow



#%%

# =============================================================================
# Leer hipsometrica
# =============================================================================
cuenca           = "RioMaipoEnElManzano"
curva_hipso      =  pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2",inplace=True)

# =============================================================================
# snow limit ianigla
# =============================================================================

sl_ianigla       =  pd.read_csv("datos/ianigla/"+cuenca+"_SCA_s_comp.filtro_MA.3días.csv",index_col=0)/100
sl_ianigla.index =  pd.to_datetime(sl_ianigla.index)
interp = interp1d(1-curva_hipso["fArea"],curva_hipso["height"])
sl_ianigla["SL"] = list(map(lambda x: interp(x).item(),sl_ianigla["SCA(%)"]))
sl_ianigla2 = pd.read_csv("datos/ianigla/RioMaipoEnElManzano_lim_nieve_ianigla_2000-2015.csv",
                          index_col=1).iloc[:,1]
sl_ianigla2.index = pd.to_datetime(sl_ianigla2.index.values)

# sl_ianigla.index = sl_ianigla.index+datetime.timedelta(hours=12)

# =============================================================================
# snow limit dgf
# =============================================================================
sl_dgf           =  pd.read_csv("datos/modis/MAIPO.txt",sep=" ",header=None)
sl_dgf.index     = pd.to_datetime(sl_dgf[0].values,format="%Y%m%d")
sl_dgf.drop([0,1,2,3],axis=1,inplace=True)
# sl_dgf.index     = sl_dgf.index+datetime.timedelta(hours=12)

# =============================================================================
# isoterma 0 radiosonda santo domingo
# =============================================================================
H0_sd        =  pd.read_csv("datos/stodomingo/isoterma0.csv",index_col=0)-300
H0_sd.index  = pd.to_datetime(H0_sd.index)
H0_sd        = H0_sd.where(H0_sd<6000).dropna()["H0_StoDomingo"]
# H0_sd = pd.Series(H0_sd.groupby([H0_sd.index.year,
#                                           H0_sd.index.month,
#                                           H0_sd.index.day]).mean().ravel(),
#                       index=np.unique(pd.to_datetime(H0_sd.index.date)),
#                       name = "H0_StoDomingo")
# H0_sd        = H0_sd.resample("d").interpolate("linear")

# =============================================================================
# pp quinta normal y Q maipo manzano
# =============================================================================
pr_qn  =  pd.read_csv("datos/estaciones/pr_quintanormal.csv").applymap(lambda x: str(x))
pr_qn["fecha"] = pr_qn["agno"]+"-"+pr_qn[" mes"]+"-"+pr_qn[" dia"]
pr_qn.index = pd.to_datetime(pr_qn["fecha"])
pr_qn.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
pr_qn  = pd.to_numeric(pr_qn[" valor"])
pr_qn.name = "pr"

qmd_mm   = pd.read_csv("datos/estaciones/qmdiario_"+cuenca+".csv").applymap(lambda x: str(x))
qmd_mm["fecha"] = qmd_mm["agno"]+"-"+qmd_mm[" mes"]+"-"+qmd_mm[" dia"]
qmd_mm.index = pd.to_datetime(qmd_mm["fecha"])
qmd_mm.drop(["fecha","agno"," mes"," dia"], inplace=True, axis=1)
qmd_mm  = pd.to_numeric(qmd_mm[" valor"])
qmd_mm.name = "q"

def corregir_qintdga(excel,yr):
# qinst_mm = pd.read_csv("datos/estaciones/qinst_"+cuenca+"_2009.csv",header=None)
    qinst_mm = pd.read_excel(excel,header=None,index_col=0,dtype=str)
    dummy = []
    month_start = np.where(qinst_mm.iloc[:,0].map(lambda x: x=="MES:"))[0]
    for pos in range(len(month_start)):
        if pos != len(month_start)-1:
            i = month_start[pos]
            j = month_start[pos+1]
            table = qinst_mm.iloc[i:j,[0,1,2,4,7,8,9,10,15,16,17,18]]
            table = table.iloc[2:,:]
            table = np.vstack((table.iloc[:,[0,1,2,3]],
                                table.iloc[:,[4,5,6,7]],
                                table.iloc[:,[8,9,10,11]]))
            table = np.hstack((np.expand_dims(np.tile(str(pos+1),len(table)),axis=1),table))
            table = np.hstack((np.expand_dims(np.tile(yr,len(table)),axis=1),table))
            dummy.append(pd.DataFrame(table))
        else:
            table = qinst_mm.iloc[j:,[0,1,2,4,7,8,9,10,15,16,17,18]]
            table = table.iloc[2:,:]
            table = np.vstack((table.iloc[:,[0,1,2,3]],
                                table.iloc[:,[4,5,6,7]],
                                table.iloc[:,[8,9,10,11]]))
            table = np.hstack((np.expand_dims(np.tile(str(pos+1),len(table)),axis=1),table))
            table = np.hstack((np.expand_dims(np.tile(yr,len(table)),axis=1),table))
            dummy.append(pd.DataFrame(table))
    qinst_mm = pd.concat(dummy,axis=0)
    qinst_mm["fecha"] = qinst_mm.iloc[:,0]+"-"+qinst_mm.iloc[:,1]+"-"+qinst_mm.iloc[:,2]+"T"+qinst_mm.iloc[:,3]
    index = []
    for i in range(len(qinst_mm.index)):
        try:
            pd.to_datetime(qinst_mm["fecha"].values[i])
            index.append(True)
        except:
            index.append(False)
    qinst_mm = qinst_mm[index]
    qinst_mm.index = pd.to_datetime(qinst_mm["fecha"])
    qinst_mm.drop(qinst_mm.columns[[0,1,2,3,4,6]],axis=1,inplace=True)
    qinst_mm = pd.to_numeric(qinst_mm[5]).rename("qinst_mm")
    qinst_mm = qinst_mm.resample("1h").max().dropna()
    return qinst_mm
qinst_mm = pd.read_csv("datos/estaciones/qinst_"+cuenca+".csv",index_col=0).qinst_mm
qinst_mm.index=pd.to_datetime(qinst_mm.index)
# =============================================================================
# Estacion dgf
# =============================================================================

datos_dgf = pd.read_csv("datos/estaciones/dgf/DATOSUTC_2004-2019.csv",index_col=0)
datos_dgf.index = pd.to_datetime(datos_dgf.index.values)
#%%
# =============================================================================
# fSCA modis
# =============================================================================

# fSCA_terra = xr.open_dataset("datos/modis/MOD10A1_2000-2021.nc",chunks="auto")
# fSCA_aqua  = xr.open_dataset("datos/modis/MYD10A1_2000-2021.nc",chunks="auto")

#%%
# =============================================================================
# era5land
# =============================================================================

# era5land_pr = xr.open_mfdataset("datos/era5land/total_precipitation*",
#                                 chunks="auto")
# era5land_pr = era5land_pr.sel(lat=-33.457236,lon=-70.661693,time="2009",
#                               method="nearest").to_dataframe()


# era5land_temp = xr.open_mfdataset("datos/era5land/2m*",
#                                   chunks="auto")
# era5land_temp = era5land_temp.sel(lat=-33.457236,lon=-70.661693,time="2009",
#                                   method="nearest").to_dataframe()
#%%
# =============================================================================
# era5land maipo
# =============================================================================

# era5land_maipo = xr.open_mfdataset("datos/era5land/maipo/*_maipo.nc",chunks="auto")

# era5land_sro   = era5land_maipo.sro.squeeze().to_series()
# era5land_ro    = era5land_maipo.ro.squeeze().to_series()

# # =============================================================================
# # 
# # =============================================================================
# era5land_sca = xr.open_mfdataset("datos/era5land/maipo/snow_cover.nc",chunks="auto")
# era5land_sca = era5land_sca.snowc.sel(time="2009").load()

#%%
# sl_era5land = np.where(era5land_sca>=70,1,0).sum(axis=1).sum(axis=1)/49
# sl_era5land = pd.Series(sl_era5land,index=era5land_sca.time.values)
# interp = interp1d(1-curva_hipso["fArea"],curva_hipso["height"])
# sl_era5land = list(map(lambda x: interp(x).item(),sl_era5land))
# sl_era5land = pd.Series(sl_era5land,index=era5land_sca.time.values)

#%%

#Case of study
cos   = pd.DataFrame([])
# dates = pd.date_range("2009-09-01","2009-09-13",freq="12h")
dates = pd.date_range("2005-06-20","2005-07-10",freq="12h")
# dates = pd.date_range("2005-08-20","2005-09-5",freq="12h")
cos["H0_sd"]      = H0_sd.reindex(dates).resample("12h").interpolate("spline",order=3)
cos["FL_sd"]      = cos["H0_sd"]-300
cos["pr_qn"]      = pr_qn.reindex(dates)
cos["qmd_mm"]     = qmd_mm.reindex(dates)
cos["temp_dgf"]   = datos_dgf["4"].resample("12h").mean().reindex(dates)
cos["pr_dgf"]     = datos_dgf["9"].resample("12h").sum().reindex(dates)
cos["sl_ianigla"] = sl_ianigla["SL"].reindex(dates)
cos["sl_dgf"]     = sl_dgf.reindex(dates) 

#%%

# =============================================================================
# Plot isotherm, snow limit heights
# =============================================================================
import matplotlib.dates as mdates
fig,ax = plt.subplots(3,1,sharex=True,figsize=(10,6))
fig.tight_layout(pad=2)
ax = ax.ravel()

cos["sl_ianigla"].rename("").dropna().plot(ax=ax[0],marker="o",lw=0,markersize=0)


# cos["sl_ianigla"].dropna().plot(ax=ax[0],marker="o",lw=0,markeredgecolor="k",
                                # color="royalblue",alpha=0.8,label="SnowLimit_IANIGLA")
ax[0].plot(cos["sl_ianigla"].dropna().index+datetime.timedelta(hours=12),
           cos["sl_ianigla"].dropna(),marker="o",
           markeredgecolor="k",color="powderblue",lw=0,alpha=0.8,label="SnowLimit_Hypsometry_IANIGLA")
# ax[0].plot(sl_era5land.index,
#            sl_era5land,
#            color="k",lw=1,alpha=0.8,label="SnowLimit_IANIGLA")
ax[0].plot(sl_ianigla2.index+datetime.timedelta(hours=12),
           sl_ianigla2,marker="o",
           markeredgecolor="k",color="royalblue",lw=0,alpha=0.8,label="SnowLimit_IANIGLA")

ax[0].plot(cos["sl_dgf"].dropna().index+datetime.timedelta(hours=12),
           cos["sl_dgf"].dropna(),marker="o",
           markeredgecolor="k",color="teal",lw=0,alpha=0.8,label="SnowLimit_Hypsometry_DGF")

ax[0].errorbar(cos["H0_sd"].dropna().index,cos["H0_sd"].dropna(),
                    marker="o",yerr=(np.ones(len(cos["H0_sd"].dropna()))*300,
                                      np.zeros(len(cos["H0_sd"].dropna()))),fmt="o",
                    capsize=3,alpha=0.8,label="Sto. Domingo Radiosonde",
                    markeredgecolor="k",color="tab:red",lw=0.6)
ax[0].set_ylim((0,5000))
ax[0].legend(loc=(0,1),ncol=2,frameon=False)
ax[0].set_ylabel("Height (m)")
ax[0].grid(True,axis="y",ls="--")

# =============================================================================
# precipitation and temperature
# =============================================================================

# ax[1].bar(cos["pr_qn"].index,cos["pr_qn"]/24,alpha=0.5,color="tab:red",width=1/4,
#           label="Precipitation_QN")
ax[1].bar(datos_dgf["9"].resample("1h").sum().reindex(pd.date_range(str(cos.index[0]),
                                                                  str(cos.index[-1]),
                                                                  freq="1h")).index,
        datos_dgf["9"].resample("1h").sum().reindex(pd.date_range(str(cos.index[0]),
                                                                  str(cos.index[-1]),
                                                                  freq="1h")),
        width=1/24,alpha=0.5,label="Precipitation_DGF")
ax1twin = ax[1].twinx()
ax1twin.plot(datos_dgf["4"].resample("1h").mean().reindex(pd.date_range(str(cos.index[0]),
                                                                        str(cos.index[-1]),
                                                                        freq="1h")).index,
              datos_dgf["4"].resample("1h").mean().reindex(pd.date_range(str(cos.index[0]),
                                                                        str(cos.index[-1]),
                                                                        freq="1h")),
              color="black",label="Temperature_DGF")
ax[1].legend(loc=(0,1),frameon=False,ncol=2)
# ax1twin.plot(era5land["t2m"]-273.15,color="tab:red",label="Temperature_ERA5Land")
ax1twin.legend(loc=(.22,1),frameon=False,ncol=2)
ax1twin.set_ylabel("Temperature (°C)")
ax[1].set_ylabel("Precipitation\nIntensity (mm/hr)")

# ax[1].plot(era5land["tp"])

# =============================================================================
# runoff
# =============================================================================



# ax[2].plot(cos["qmd_mm"].dropna(),marker="o",markeredgecolor="purple",lw=0,
#             color="purple")

q = qinst_mm
ax[2].plot(q,label="Direct Runoff", color="darkblue")
ax[2].plot(sliding_interval_filter(q,48)[0], label="Baseflow", color="forestgreen")
# ax[2].plot(q[argrelextrema(q.rolling(4).mean().values, np.less)[0]].index,
#             q[argrelextrema(q.rolling(4).mean().values, np.less)[0]],
#             color="forestgreen",label="Baseflow")

# ax[2].plot((81e5/3600*era5land_sro.reindex(pd.date_range(str(cos.index[0]),
#                                                          str(cos.index[-1]),
#                                                          freq="1h"))))
ax[2].set_ylabel("Runoff (m3/s)")
ax[2].set_ylim(0,500)
ax[2].legend()


# =============================================================================
# general stuff
# =============================================================================
loc = mdates.MonthLocator(interval=1)
fmt = mdates.DateFormatter('%b\n%Y')
for axis in ax:
    axis.grid(True,axis="x",which="both",ls="--")
    axis.xaxis.set_major_locator(loc)
    axis.xaxis.set_major_formatter(fmt)


plt.savefig("plots/maipomanzano/timeseries_case.pdf",dpi=150,bbox_inches="tight")

#%%
