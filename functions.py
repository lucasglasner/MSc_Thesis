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
        gl = axis.gridlines(draw_labels=True)
