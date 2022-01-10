#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 16:16:54 2021

@author: lucas

# =============================================================================
# ROS in Maipo Basin: 2013/08/11 case of study and time series analysis
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
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cf
import scipy.stats as st
import matplotlib.colors as colors
import matplotlib.ticker as mticker
from metpy.units import units
import matplotlib as mpl
from metpy.calc import mixing_ratio_from_relative_humidity


def smooth(x, window_len=11, window='hanning'):
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

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError

    s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
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
    baseflow = pd.Series(data=ts.groupby(
        intervals).transform('min'), index=ts.index)
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
        baseflow = pd.Series(data=minimum_filter1d(
            ts, size, mode='reflect'), index=ts.index)

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
    baseflow_min = pd.Series(generic_filter(
        ts, _local_minimum, footprint=np.ones(size)), index=ts.index)
    baseflow = baseflow_min.interpolate(method='linear')
    # interpolation between values may lead to baseflow > streamflow
    errors = (baseflow > ts)
    while errors.any():
        print('hello world')
        error_labelled, n_features = label(errors)
        error_blocks = [ts[error_labelled == i]
                        for i in range(1, n_features + 1)]
        error_local_min = [argrelextrema(e.values, np.less)[
            0] for e in error_blocks]
        print(error_local_min)
        break
    quickflow = ts - baseflow
    baseflow.name = 'baseflow'
    quickflow.name = 'quickflow'

    return baseflow, quickflow


# %%
date = "2013-08-11"
# date = "2002-06-01"
interval = datetime.datetime.strptime(date, '%Y-%m-%d')
interval = slice(interval-datetime.timedelta(days=8),
                 interval+datetime.timedelta(days=5))

# =============================================================================
# basin polygon and hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso = pd.read_csv("datos/topography/basins/hipso/"+cuenca+"_Hipso.csv")
curva_hipso.drop_duplicates(subset="Area_km2", inplace=True)
basin = gpd.read_file('datos/vector/'+cuenca+'.shp')


# =============================================================================
# SL_mm
# =============================================================================

SL_mm = pd.read_csv('datos/snowlimits_maipomanzano.csv',
                    index_col=0).squeeze()
SL_mm.index = pd.to_datetime(SL_mm.index)+datetime.timedelta(hours=12)
SL_mm = SL_mm.dropna(how='all')

times = pd.date_range("2000-01-01", "2021-12-31", freq="d")
for i in range(3):
    y = SL_mm.iloc[:, i].dropna()
    x = SL_mm["IANIGLA"].dropna().reindex(y.index)
    m = st.linregress(x, y)
    for j in range(len(SL_mm)):
        if np.isnan(SL_mm.iloc[j, i]):
            SL_mm.iloc[j, i] = m.slope * \
                SL_mm["IANIGLA"].values[j]+m.intercept


# =============================================================================
# ISOTERMAS 0
# =============================================================================

H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)+datetime.timedelta(hours=4)
H0_mm = H0_mm.dropna(how='all')
# =============================================================================
# Q y pr maipo manzano
# =============================================================================

qinst_mm = pd.read_csv("datos/estaciones/qinst_" +
                       cuenca+".csv", index_col=0).qinst_mm
qinst_mm.index = pd.to_datetime(qinst_mm.index)

pr_mm = pd.read_csv('datos/estaciones/pr_RioMaipoEnElManzano_2013-08.csv')
pr_mm.index = pd.to_datetime(pr_mm['Fecha'])
pr_mm = pr_mm['Valor'].drop_duplicates()
pr_mm = pr_mm.reindex(pd.date_range(interval.start,
                                    interval.stop, freq='h')).fillna(0)
# =============================================================================
# Estacion dgf
# =============================================================================
datos_dgf = pd.read_csv(
    "datos/estaciones/dgf/DATOSUTC_2004-2019.csv", index_col=0)
datos_dgf.index = pd.to_datetime(
    datos_dgf.index.values)-datetime.timedelta(hours=4)
# =============================================================================
# pr cr2met
# =============================================================================

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv", index_col=0)
pr_cr2met.index = pd.to_datetime(
    pr_cr2met["date"])+datetime.timedelta(hours=12)
pr_cr2met.drop("date", axis=1, inplace=True)
pr_cr2met = pr_cr2met.squeeze()


# %%
h0_mm = H0_mm.dropna().reindex(pr_cr2met.index, method='nearest')-300
sl_mm = SL_mm['MODIS_H20'].dropna().reindex(pr_cr2met.index)
ROS = h0_mm-sl_mm
ROS = ROS[pr_cr2met > 3].dropna()
ROS = pd.concat([ROS, pr_cr2met.reindex(ROS.index)], axis=1)
ROS.columns = ["dH", 'pr']

# %%
# date = "2010-06-23"


fig, ax = plt.subplots(4, 1, sharex=True, figsize=(8, 8))
fig.tight_layout(pad=0.8)
plt.rc('font', size=18)
# =============================================================================
# plot snowlimit and freezinglevel
# =============================================================================
ax[0].plot(SL_mm['MODIS_H50'][interval], color='lightblue', marker='o', mec='k',
           ms=7, ls=":", label='Snow Limit', zorder=10)
# ax[0].plot(h0_mm[interval], color="tab:red", ls=":",
#             marker="o", ms=7, mec="k", label='Zero Degree Level')
ax[0].errorbar(H0_mm[interval].index, H0_mm[interval],
               yerr=300, color='tab:red', ls=':',
               marker='o', ms=7, mec='k', label='Zero Degree Level',
               uplims=True, zorder=30)
ax[0].set_yticks(np.arange(1e3, 5e3, 1e3))
ax[0].set_ylim(0, 5e3)
ax[0].set_ylabel('Height\n$(m.a.s.l)$')
ax[0].legend(frameon=False, loc=(0, 1.01), ncol=2, fontsize=16)
ax00 = ax[0].twinx()
ax00.plot(datos_dgf[interval].iloc[:, 7], color='k', alpha=0.3, label='SWR')
ax00.set_ylim(0, 5e2)
ax00.set_yticks(np.arange(100, 500, 100))
ax00.set_ylabel('SWR\n$(W/m^2)$')

# =============================================================================
# plot precipitation and temperature
# =============================================================================
ax[1].bar(pr_mm[interval].index,
          datos_dgf[interval].iloc[:, 9].resample('h').sum(),
          width=1/24,
          zorder=1, color='royalblue', label='DGF-Roof')

# ax[1].bar(pr_cr2met[interval].index, pr_cr2met[interval], alpha=0.5,
#           zorder=0, color='cadetblue', label='CR2MET Basin Mean',
#           width=1, edgecolor='k')

ax[1].bar(pr_mm[interval].index, pr_mm[interval], alpha=0.5,
          zorder=0, color='cadetblue', label='MaipoEnElManzano',
          width=1/24)
ax[1].set_ylim(0, 8)
ax[1].set_yticks(np.arange(2, 8, 2))
ax[1].set_ylabel('Precipitation\n$(mm/hour)$')
ax[1].legend(loc=(0, 1.01), fontsize=16, frameon=False, ncol=2)
ax11 = ax[1].twinx()
ax11.plot(datos_dgf[interval].index, datos_dgf[interval].iloc[:, 5],
          color="k", lw=1, alpha=0.8)
ax11.set_ylim(0, 30)
ax11.set_yticks(np.arange(5, 30, 5))
ax11.set_ylabel('Temperature\n$(°C)$')

# =============================================================================
# plot wind and humidity
# =============================================================================

ax[2].plot(datos_dgf[interval].iloc[:, 8]+1e3, color='darkorchid', zorder=10)
ax[2].set_ylim(947, 965)
ax[2].set_yticks(np.arange(950, 965, 3))
ax[2].set_ylabel('Barometric\nPressure $(mb)$')
ax22 = ax[2].twinx()
ax22.plot(datos_dgf[interval].iloc[:, 12], color='tab:green')
ax22.set_ylim(0, 6)
ax22.set_yticks(np.arange(1, 6, 1))
ax22.set_ylabel('Wind Speed\n$(m/s)$')

# =============================================================================
# plot streamflow
# =============================================================================

ax[3].plot(qinst_mm[interval], color='darkblue', label='Surface Runoff')
ax[3].set_ylim(30, 90)
ax[3].set_yticks(np.arange(40, 90, 15))
ax[3].plot(sliding_interval_filter(qinst_mm[interval], 40)[0],
           color='chocolate', label='Base Flow')
ax[3].set_ylabel('Runoff\n$(m^3/s)$')
ax[3].legend(loc=(0, 1.01), fontsize=16, ncol=2, frameon=False)

for axis in ax:
    axis.axvspan("2013-08-06", "2013-08-09", color='k', alpha=0.15)
    axis.axvspan("2013-08-11", "2013-08-13", color='k', alpha=0.15)
    axis.set_xlim(interval.start, interval.stop)
    # axis.tick_params(axis='x',rotation=45)
    axis.grid(which='both', ls=":", axis="x")
    axis.xaxis.set_minor_locator(mpl.dates.DayLocator(interval=1))
    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%d'))
    axis.xaxis.set_major_locator(mpl.dates.DayLocator(interval=19))
    axis.xaxis.set_major_formatter(mpl.dates.DateFormatter('\n%d\n%b-%Y'))
    for maj in axis.xaxis.get_major_ticks():
        maj.label.set_fontsize(18)
    for m in axis.xaxis.get_minor_ticks():
        m.label.set_fontsize(16)

plt.savefig('plots/caseofstudy_Aug2013/ROS_maipo_caseofstudy.pdf', dpi=150,
            bbox_inches='tight')
