#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 10:32:14 2022

@author: lucas

# =============================================================================
# This Script contains functions with random purpose
# =============================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import sys
from scipy.ndimage.filters import minimum_filter1d, generic_filter
from scipy.ndimage.measurements import label
from scipy.signal import argrelextrema


def corregir_qintdga(excel, yr):
    """
    Parameters
    ----------
    excel [str]: Path to MS Excel Spreadsheet with DGA data.
    yr [str]: Year

    Returns
    -------
    qinst_mm [DataFrame]: Pandas DataFrame with the data in human style.

    """
    # qinst_mm = pd.read_csv("datos/estaciones/qinst_"+cuenca+"_2009.csv",header=None)
    qinst_mm = pd.read_excel(excel, header=None, index_col=0, dtype=str)
    dummy = []
    month_start = np.where(qinst_mm.iloc[:, 0].map(lambda x: x == "MES:"))[0]
    for pos in range(len(month_start)):
        if pos != len(month_start)-1:
            i = month_start[pos]
            j = month_start[pos+1]
            table = qinst_mm.iloc[i:j, [
                0, 1, 2, 4, 7, 8, 9, 10, 15, 16, 17, 18]]
            table = table.iloc[2:, :]
            table = np.vstack((table.iloc[:, [0, 1, 2, 3]],
                               table.iloc[:, [4, 5, 6, 7]],
                               table.iloc[:, [8, 9, 10, 11]]))
            table = np.hstack(
                (np.expand_dims(np.tile(str(pos+1), len(table)), axis=1), table))
            table = np.hstack(
                (np.expand_dims(np.tile(yr, len(table)), axis=1), table))
            dummy.append(pd.DataFrame(table))
        else:
            table = qinst_mm.iloc[j:, [
                0, 1, 2, 4, 7, 8, 9, 10, 15, 16, 17, 18]]
            table = table.iloc[2:, :]
            table = np.vstack((table.iloc[:, [0, 1, 2, 3]],
                               table.iloc[:, [4, 5, 6, 7]],
                               table.iloc[:, [8, 9, 10, 11]]))
            table = np.hstack(
                (np.expand_dims(np.tile(str(pos+1), len(table)), axis=1), table))
            table = np.hstack(
                (np.expand_dims(np.tile(yr, len(table)), axis=1), table))
            dummy.append(pd.DataFrame(table))
    qinst_mm = pd.concat(dummy, axis=0)
    qinst_mm["fecha"] = qinst_mm.iloc[:, 0]+"-"+qinst_mm.iloc[:, 1] + \
        "-"+qinst_mm.iloc[:, 2]+"T"+qinst_mm.iloc[:, 3]
    index = []
    for i in range(len(qinst_mm.index)):
        try:
            pd.to_datetime(qinst_mm["fecha"].values[i])
            index.append(True)
        except:
            index.append(False)
    qinst_mm = qinst_mm[index]
    qinst_mm.index = pd.to_datetime(qinst_mm["fecha"])
    qinst_mm.drop(qinst_mm.columns[[0, 1, 2, 3, 4, 6]], axis=1, inplace=True)
    qinst_mm = pd.to_numeric(qinst_mm[5]).rename("qinst_mm")
    qinst_mm = qinst_mm.resample("1h").max().dropna()
    return qinst_mm


def add_labels(geoaxes, yticks=[], xticks=[], **kwargs):
    if "mpl" in dir():
        pass
    else:
        import matplotlib as mpl

    if len(geoaxes.shape) == 2:
        # gridlines = []
        for axis in geoaxes.ravel():
            gl = axis.gridlines(**kwargs)
            # gridlines.append(gl)
        for axis in geoaxes[-1, :]:
            gl = axis.gridlines(draw_labels=True, **kwargs)
            gl.xlocator = mpl.ticker.FixedLocator(xticks)
            gl.ylocator = mpl.ticker.FixedLocator([])
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = False
        for axis in geoaxes[:, 0]:
            gl = axis.gridlines(draw_labels=True, **kwargs)
            gl.xlocator = mpl.ticker.FixedLocator([])
            gl.ylocator = mpl.ticker.FixedLocator(yticks)
            gl.top_labels = False
            gl.right_labels = False
            gl.bottom_labels = False
        # return np.array(gridlines).reshape(geoaxes.shape)
    if len(geoaxes.shape) == 1:
        gl = axis.gridlines(draw_labels=True,**kwargs)


def seasonal_decompose(ts, period, nharmonics=3, bandwidth=2):
    """
    Parameters
    ----------
    ts : Time series data in a pandas series format, with timestamps
         in the index.
    period : period of the season
    nharmonics : Number of harmonics to remove, default is 3.

    Returns
    -------
    season : Seasonal component of the time series.
    anomaly : The time series anomaly without the seasonal cycle.
    """
    n = len(ts)
    ft = np.fft.fft(ts)
    ft[0] = 0  # Remove mean#
    for i in range(nharmonics):  # Filter cycle#
        pos = n//(period//(i+1))
        ft[pos-bandwidth:pos+bandwidth] = 0
        ft[n-pos-bandwidth:n-pos+bandwidth] = 0
        # ft[pos]=0
        # ft[n-pos]=0
    anomaly = np.fft.ifft(ft).real
    anomaly = pd.Series(anomaly, index=ts.index)
    season = ts-anomaly
    return season, anomaly


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
    win_center_ix = len(window) // 2
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

    origin = int(size) // 2
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
