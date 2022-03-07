#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:27:42 2022

@author: lucas

# =============================================================================
# This Script does some hydrological analysis in the case of study
# =============================================================================
"""

from functions import add_labels
from glob import glob
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import datetime
from scipy.interpolate import interp1d
import scipy.stats as st
from scipy.ndimage.filters import minimum_filter1d, generic_filter
from scipy.ndimage.measurements import label
from scipy.signal import argrelextrema
import xarray as xr
import cartopy.crs as ccrs
import geopandas as gpd
import cartopy.feature as cf


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
# =============================================================================
# big time interval and graph time interval
# =============================================================================
interval = slice(datetime.datetime(2013, 8, 1),
                 datetime.datetime(2013, 8, 16))
interval2 = slice(datetime.datetime(2013, 8, 11),
                  datetime.datetime(2013, 8, 16))


# %%
# =============================================================================
# PRECIPITATION ON BASIN OUTLETS
# =============================================================================


pr_mm = pd.read_csv('datos/estaciones/pr_RioMaipoEnElManzano_2013-08.csv')
pr_mm.index = pd.to_datetime(pr_mm['Fecha'])
pr_mm = pr_mm['Valor'].drop_duplicates()
pr_mm = pr_mm.reindex(pd.date_range(interval.start,
                                    interval.stop, freq='h')).fillna(0)

pr_teno = pd.read_csv(
    'datos/estaciones/pr_RioTenoDespuesDeJuntaConClaro_2013-08.csv')
pr_teno.index = pd.to_datetime(pr_teno.Fecha)
pr_teno = pr_teno['Valor'].reindex(pr_mm.index)[interval]


pr_uble = pd.read_csv('datos/estaciones/pr_RioUbleEnSanFabian_2013-08.csv')
pr_uble.index = pd.to_datetime(pr_uble.Fecha)
pr_uble = pr_uble['Valor'].reindex(pr_mm.index)[interval]


pr = pd.concat([pr_mm, pr_teno, pr_uble], axis=1)
pr.columns = columns = ['Rio Maipo En El Manzano',
                        'Rio Teno Despues De Junta Con Claro',
                        'Rio Uble En San Fabian N 2']
del pr_mm, pr_teno, pr_uble

# %%
# =============================================================================
# hypsometric curves
# =============================================================================


hypso0 = pd.read_csv(
    'datos/topography/basins/hypso/RioMaipoEnElManzano_hypso.csv')
hypso1 = pd.read_csv(
    'datos/topography/basins/hypso/RioTenoDespuesDeJuntaConClaro_hypso.csv')
hypso2 = pd.read_csv(
    'datos/topography/basins/hypso/RioUbleEnSanFabianN2_hypso.csv')


hypso = pd.concat([hypso0, hypso1, hypso2], keys=pr.columns, axis=1)

areas = [hypso[b].Area_km2.max() for b in pr.columns]
areas = pd.Series(areas, index=pr.columns)
int_func = [interp1d(hypso[b].height, hypso[b].fArea) for b in pr.columns]

# %%
# =============================================================================
# ZERO DEGREE LEVEL AND PLUVIAL AREA
# =============================================================================

H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)-datetime.timedelta(hours=4)
H0_mm = H0_mm.reindex(pr.index)[interval].interpolate('cubicspline')

pluv_area = [int_func[i](H0_mm-300)*areas[i] for i in range(3)]
nonpluv_area = [(1-int_func[i](H0_mm-300))*areas[i] for i in range(3)]
pluv_area = pd.DataFrame(pluv_area, columns=pr.index, index=pr.columns).T
nonpluv_area = pd.DataFrame(nonpluv_area, columns=pr.index, index=pr.columns).T

# %%
# =============================================================================
# SNOW COVER
# =============================================================================


SCA0 = pd.read_csv(
    'datos/ianigla/RioMaipoEnElManzano_SCA_s_comp.filtro_MA.3días.csv',
    index_col=0)
SCA1 = pd.read_csv(
    'datos/ianigla/RioTenoDespuesDeJuntaConClaro_SCA_s_comp.filtro_MA.3días.csv',
    index_col=0)
SCA2 = pd.read_csv(
    'datos/ianigla/RioUbleEnSanFabianN2_SCA_s_comp.filtro_MA.3días.csv.csv',
    index_col=0)

SCA0.index = pd.to_datetime(SCA0.index)
SCA1.index = pd.to_datetime(SCA1.index)
SCA2.index = pd.to_datetime(SCA2.index)

SCA = pd.concat([s.iloc[:, 1]
                for s in [SCA0, SCA1, SCA2]], axis=1).reindex(pr.index)
SCA.columns = pr.columns
SCA = SCA[interval].interpolate('cubicspline')/100

snow_area = [SCA[b]*areas.loc[b] for b in pr.columns]
snow_area = pd.concat(snow_area, axis=1)

# %%
# =============================================================================
# ROS AREA
# =============================================================================

ros_area = np.clip(snow_area-nonpluv_area, 0, 7e3)

# %%
# =============================================================================
# RUNOFF DATA IN DIFFERENT BASINS
# =============================================================================


basin_attributes = pd.read_csv('datos/basins_attributes.csv', index_col=1)
basin_attributes = basin_attributes.sort_values(
    by='gauge_name', ascending=False)
basin_attributes = basin_attributes.loc[pr.columns]
runoff2 = pd.read_csv('datos/runoff_gauges_dataset.csv', index_col=0)
runoff2.index = pd.to_datetime(runoff2.index)
runoff2 = runoff2[pr.columns]
runoff = runoff2[interval]
runoff = runoff.T.dropna(how='all').T


baseflows = [local_minimum_filter(runoff[b], 30)[0] for b in pr.columns]
baseflows = pd.concat(baseflows, axis=1)
baseflows.columns = pr.columns

quickflows = [runoff[c]-baseflows[c] for c in pr.columns]
quickflows = pd.concat(quickflows, axis=1)
quickflows = quickflows.where(quickflows > 0).fillna(0)
# %%
# =============================================================================
# RASTER DATA, SWE, TOPOGRAPHY AND PLUVIAL AREA MASKS
# =============================================================================

SWE0 = xr.open_dataset(
    'datos/ANDES_SWE_Cortes/maipomanzano/ANDES_SWE_WY2014.nc')
SWE1 = xr.open_dataset('datos/ANDES_SWE_Cortes/ANDES_SWE_WY2014_Teno.nc')
SWE2 = xr.open_dataset('datos/ANDES_SWE_Cortes/ANDES_SWE_WY2014_Uble.nc')

SWE = [SWE0.SWE.sel(time=interval),
       SWE1.SWE.sel(time=interval),
       SWE2.SWE.sel(time=interval)]

del SWE0, SWE1, SWE2

dSWE = [swe.diff('time') for swe in SWE]

# dem0 = xr.open_dataset('datos/topography/basins/RioMaipoEnElManzano.nc')
# dem1 = xr.open_dataset(
#     'datos/topography/basins/RioTenoDespuesDeJuntaConClaro.nc')
# dem2 = xr.open_dataset('datos/topography/basins/RioUbleEnSanFabianN2.nc')

# dem = [dem0.Band1, dem1.Band1, dem2.Band1]
# del dem0, dem1, dem2

# new_dem = []
# for d, swe in zip(dem, SWE):
#     nd = d.reindex({'lat': swe.lat, 'lon': swe.lon}, method='nearest')
#     new_dem.append(nd)
# del dem
# dem = new_dem

melt = [dswe.where(dswe < 0).mean(dim=['lat', 'lon']).to_series()
        for dswe in dSWE]
melt = pd.concat(melt, axis=1)*-1
melt.columns = pr.columns
melt = melt.reindex(pr.index).interpolate(method='cubicspline')/24
# %%
# =============================================================================
# compute flood stats
# =============================================================================
qmax = runoff[interval2].max()
pr_max = pr[interval2].max()
pr_cum = pr[interval2].sum()

rain = pr[interval2].where(pr[interval2] > 0).dropna(how='all')
start_rain = [rain[c].dropna().index[0] for c in rain.columns]
end_rain = [rain[c].dropna().index[-1] for c in rain.columns]
start_rain = np.array(start_rain)
end_rain = np.array(end_rain)
rain_duration = end_rain-start_rain
# rain_duration = np.array([rd.total_seconds() for rd in rain_duration])/3600
rain_duration = pd.Series(rain_duration, index=qmax.index)
rain_duration = rain_duration.map(lambda x: int(x.total_seconds()/3600))
tau = (runoff)[interval2].idxmax()-start_rain
tau = tau.map(lambda x: int(x.total_seconds()/3600))

melt_cum = melt[interval2][pr[interval2] > 0].sum()


stats = pd.concat([qmax, pr_max, pr_cum, melt_cum, rain_duration, tau], axis=1)
stats.columns = ['Qmax', 'PRmax', 'PRcum', 'Snowmelt', 'RainDuration', 'tau']

stats = stats.round(2)
# %%
plt.rc("font", size=18)
fig, ax = plt.subplots(2, 3, figsize=(14, 4), sharex='col',
                       gridspec_kw={'height_ratios': [1, 2]},
                       sharey='row')
titles = ['Rio Maipo\nEn El Manzano',
          'Rio Teno Despues\nDe Junta Con Claro',
          'Rio Uble En\nSan Fabian N 2']

# these are matplotlib.patch.Patch properties
props = dict(facecolor='white', lw=1)
names = ['$Q_{max}$: ', '$PR_{max}$: ', '$PR_{cum}$: ', '$Melt_{cum}$: ',
         '$RainDuration$: ', r'$\tau_{peak}$: ']
units = ['$m^3/s$', '$mm/h$', '$mm$', '$mm$', '$hr$', '$hr$']

# # place a text box in upper left in axes coords
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
#         verticalalignment='top', bbox=props)

for i, axis in enumerate(ax[0, :]):
    axis.bar(pr.index-datetime.timedelta(hours=1),
             pr.iloc[:, i], color='cadetblue',
             width=0.05, align='edge', label='Precipitation',
             edgecolor='k', linewidth=0.0, zorder=1)
    axis.bar(melt.index-datetime.timedelta(hours=1),
             melt.iloc[:, i], edgecolor='k', linewidth=0.0, color='darkviolet',
             width=0.05, align='edge', label='Snowmelt',
             bottom=pr.iloc[:, i].fillna(0), zorder=2)
    axis.set_yticks(np.arange(0, 6+1.5, 1.5))
    axis.set_yticklabels([0, "", 3, "", 6])

    axis.grid(True, which='major', ls=":", zorder=0)
    axis.set_title(titles[i], loc='left', pad=25)

    box = stats.iloc[i, :]
    box = ['{:.1f}'.format((round(box[i], 1))) for i in range(len(box))]

    textstr = '\n'.join([n+b+u for n, b, u in zip(names, box, units)])
    axis.text(0.6, 1.2, textstr, transform=axis.transAxes, fontsize=10,
              verticalalignment='top', bbox=props)
    # axis.set_ylim(0,6)


for i, axis in enumerate(ax[1, :]):
    axis.plot(runoff.index,
              quickflows.iloc[:, i],
              color='darkblue')
    axis.set_ylim(0, quickflows.max().max()*1.05)
    axis1 = axis.twinx()
    axis1.set_ylim(0, 1)
    axis1.set_yticks([0, 0.25, 0.5, 0.75, 1])

    axis1.plot(pr.index, pluv_area.iloc[:, i]/areas.iloc[i],
               color='tab:red',
               ls="--", alpha=0.2)
    axis1.plot(pr.index, pluv_area.where(pr > 0).iloc[:, i]/areas.iloc[i],
               color='tab:red')
    ros = ros_area.iloc[:, i]
    ros = ros.where(ros > 0)
    ros = ros.where(pr.iloc[:, i] > 0)
    axis1.plot(pr.index,
               ros/areas.iloc[i],
               color='tab:green', label='ROS Area')

    if i < 2:
        axis1.set_yticklabels([""]*5)
    axis.grid(axis='x', which='major', ls=":")
    axis.set_yticks([0, 25, 50, 75, 100])
    axis.xaxis.set_major_formatter(mpl.dates.DateFormatter('\n\n%d'))
    axis.xaxis.set_major_locator(mpl.dates.DayLocator(interval=1))

    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
    axis.xaxis.set_minor_locator(
        mpl.dates.HourLocator(byhour=np.arange(0, 24, 6)))
    axis.tick_params(axis='x', which='minor', rotation=45)

    for maj in axis.xaxis.get_major_ticks():
        maj.label.set_fontsize(18)
    for m in axis.xaxis.get_minor_ticks():
        m.label.set_fontsize(10)
    axis.set_xlim([15927.7, 15933])
ax[-1, 0].set_ylabel('$(m^3/s)$')
ax[0, 0].set_ylabel('$(mm/h)$')
ax[1, 0].plot([], [], color='tab:green', label='ROS Area')
ax[1, 0].plot([], [], color='tab:red', label='Pluvial Area')
ax[1, 0].legend(frameon=False, fontsize=14, loc='upper left')

axis1.set_ylabel('Fraction of\ntotal Area (-)')


box = ax[1, 0].get_position()
fig.text(box.xmin*1.15, box.ymin*-1.15,
         '\nAug\n2013', ha='center', va='center')

plt.savefig('plots/caseofstudy_Aug2013/flood_study.pdf',
            dpi=150, bbox_inches='tight')
# ax[1, 1].plot(pr.index, (q1-baseflow1)[interval2], color='darkblue')
# ax[1, 2].plot(pr.index, (q2-baseflow2)[interval2], color='darkblue')
# # ax1 = ax[1].twinx()
# # ax1.set_ylim(0, 1)
# # ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
# # ax1.plot(pr.index, ap/area, color='tab:red', alpha=0.2, ls="--")
# # ax1.plot(pr.index, ap.where(pr > 0)/area, color='tab:red',
# #          label='Pluvial Area')
# # ax1.plot(pr.index, ros_area.where(ros_area > 0).where(pr > 0)/area,
# #          color='purple', label='ROS Area')
# # ax1.set_ylabel('Fraction of\ntotal Area (%)')
# # ax1.legend(frameon=False, fontsize=12, loc='upper right')
# # # ax1.set_yticklabels(ax1.get_yticks()*100)
# # # ax1.yaxis.set_major_formatter(FormatStrFormatter('%i'))
# # # ax[1].xaxis.set_major_formatter(
# # #       mpl.dates.ConciseDateFormatter(ax[1].xaxis.get_major_locator()))


# %%
# SWE = xr.open_dataset(
#     'datos/ANDES_SWE_Cortes/maipomanzano/ANDES_SWE_WY2014.nc')
# SWE = SWE.SWE.sel(time=interval)
# dSWE = SWE.diff('time')


# masks = []
# # Pluvial area masks
# pluv_area_daily = pluv_area.resample('d').mean()
# for i, pa in enumerate(pluv_area_daily):
#     if np.isnan(pa):
#         mask = dem*0
#         masks.append(mask)
#     else:
#         mask = xr.where((dem < pa), 1, 0)
#         masks.append(mask)

# masks = xr.concat(masks, pluv_area_daily.index).rename(
#     {'timestamp_UTC': 'time'})

# %%
# =============================================================================
# Build SWE loss/gain series
# =============================================================================
# parea = (6.4e6)**2*np.cos(np.deg2rad(SWE.lat.mean()))
# parea = parea*np.deg2rad(0.001)*np.deg2rad(0.001)
# parea = parea.item()

# melt = dSWE.where(dSWE < 0).mean(dim=['lat', 'lon']).to_series()
# gain = dSWE[:-1, :, :].where(masks[1:, :, :])
# gain = gain.where(gain > 0).mean(dim=['lat', 'lon']).to_series()

# %%
# # =============================================================================
# # Build table for document
# # =============================================================================
# datos = []
# for dat in [SL_mm, SCA, melt*-1, H0_mm,
#             datos_dgf.iloc[:, 9].resample('d').sum(),
#             pr_mm.resample('d').sum(),
#             qinst_mm.resample('d').max(),
#             datos_dgf.iloc[:, 5].resample('d').mean(),
#             datos_dgf.iloc[:, 5].resample('d').max(),
#             datos_dgf.iloc[:, 5].resample('d').min()]:
#     datos.append(dat[interval])

# datos = pd.concat(datos, axis=1)
# datos = datos.resample('d').mean().iloc[:-1, :]
# datos.columns = ["SL", "SCA", "MELT", "H0", "PR_DGF", "PR_MM", "Qmax", "T",
#                  "Tmax", "Tmin"]
# datos = np.round(datos, 1)
# datos = datos["2013-08-03":"2013-08-16"]

# # %%
# # =============================================================================
# # flood data
# # =============================================================================


# pr = pr_mm
# q = qinst_mm
# ap = pluv_area.reindex(q.index).interpolate(method='cubicspline')

# q1 = runoff['Rio Teno Despues De Junta Con Claro']
# q2 = runoff['Rio Uble En San Fabian N 2']

# snow_area = SCA.reindex(pluv_area.index, method='nearest')/100*area
# snow_area = snow_area.reindex(q.index).interpolate(method='cubicspline')
# ros_area = ap-(area-snow_area)
# # ros_area = ros_area/ap
# mlt = -1*melt.reindex(q.index).interpolate(method='cubicspline')/24
# gin = gain.reindex(q.index).interpolate(method='cubicspline')/24

# baseflow = sliding_interval_filter(q, 40)[0]
# baseflow1 = sliding_interval_filter(q1, 40)[0]
# baseflow2 = sliding_interval_filter(q2, 40)[0]

# mask = pr > 0.01


# pr = pr[interval2]
# q = q[interval2]
# ap = ap[interval2]
# snow_area = snow_area[interval2]
# ros_area = ros_area[interval2]
# mlt = mlt[interval2]
# gin = gin[interval2]
# baseflow = baseflow[interval2]
# mask = mask[interval2]

# rainpulse = slice('2013-08-11T14:00:00', '2013-08-12T15:00:00')
# floodpulse = slice('2013-08-11T18:00:00', '2013-08-13T09:00:00')


# # %%
# plt.rc("font", size=18)
# fig, ax = plt.subplots(2, 3, figsize=(14, 4), sharex=True,
#                        gridspec_kw={'height_ratios': [1, 2]})

# ax[0, 0].bar(pr.index-datetime.timedelta(hours=1),
#              pr, width=0.037, color='cadetblue',
#              align='edge', label='Precipitation')
# ax[0, 0].set_yticks([0, 1, 2, 3])
# ax[0, 0].set_ylim(0, 3)
# # # ax[0].bar(melt.index, melt*-1, edgecolor="k")
# ax[0, 0].bar(pr.index-datetime.timedelta(hours=1),
#              mlt,
#              width=0.037, color='darkviolet',
#              align='edge', label='Snowmelt', bottom=pr)
# # # ax[0].bar(pr.index-datetime.timedelta(hours=1),
# # #             gin,
# # #             width=0.037, edgecolor='k', color='orange')
# ax[0, 0].set_ylabel('(mm/h)')
# ax[0, 0].legend(frameon=False, loc=(0.01, 0.95), ncol=2, fontsize=12)
# # # ax[0].set_yscale("log")

# titles = ['Rio Maipo\nEn El Manzano',
#           'Rio Teno Despues\nDe Junta Con Claro',
#           'Rio Uble En\nSan Fabian N 2']
# for i, axis in enumerate(ax[0, :]):
#     axis.spines['top'].set_visible(False)
#     axis.spines['right'].set_visible(False)
#     axis.set_title(titles[i], loc='left', pad=30)
# ax[1, 0].plot(pr.index, q-baseflow, label='Direct\nRunoff', color='darkblue')
# # ax[1, 0].set_yticks([0, 7, 14, 21, 28])
# # ax[1, 0].set_ylim(0, 28)
# # ax[1].scatter((q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]].index,
# #               (q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]],
# #               ec="k",zorder=10)
# # ax[1, 0].axvspan("2013-08-11T18:00:00",
# # "2013-08-13T09:00:00", alpha=0.15, color='k')
# # ax[1].axvline("2013-08-11T14:00",color='k',alpha=0.5, ls=":")
# ax[1, 0].set_ylabel('$(m^3/s)$')
# ax[1, 0].legend(loc='upper right', frameon=False, fontsize=12)


# ax[1, 1].plot(pr.index, (q1-baseflow1)[interval2], color='darkblue')
# ax[1, 2].plot(pr.index, (q2-baseflow2)[interval2], color='darkblue')
# # ax1 = ax[1].twinx()
# # ax1.set_ylim(0, 1)
# # ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
# # ax1.plot(pr.index, ap/area, color='tab:red', alpha=0.2, ls="--")
# # ax1.plot(pr.index, ap.where(pr > 0)/area, color='tab:red',
# #          label='Pluvial Area')
# # ax1.plot(pr.index, ros_area.where(ros_area > 0).where(pr > 0)/area,
# #          color='purple', label='ROS Area')
# # ax1.set_ylabel('Fraction of\ntotal Area (%)')
# # ax1.legend(frameon=False, fontsize=12, loc='upper right')
# # # ax1.set_yticklabels(ax1.get_yticks()*100)
# # # ax1.yaxis.set_major_formatter(FormatStrFormatter('%i'))
# # # ax[1].xaxis.set_major_formatter(
# # #       mpl.dates.ConciseDateFormatter(ax[1].xaxis.get_major_locator()))


# # dtFmt =  # define the formatting
# for axis in ax[1, :]:
#     axis.xaxis.set_major_formatter(mpl.dates.DateFormatter('\n\n%d'))
#     axis.xaxis.set_major_locator(mpl.dates.DayLocator(interval=1))

#     axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
#     axis.xaxis.set_minor_locator(
#         mpl.dates.HourLocator(byhour=np.arange(0, 24, 6)))
#     axis.tick_params(axis='x', which='minor', rotation=45)

#     for maj in axis.xaxis.get_major_ticks():
#         maj.label.set_fontsize(18)
#     for m in axis.xaxis.get_minor_ticks():
#         m.label.set_fontsize(12)

# box = ax[1, 0].get_position()
# fig.text(box.xmin, box.ymin*-1.15, '2013-Aug', ha='center', va='center')
# # plt.show()

# # ax[1].set_xlim([15928.4, 15931])

# # tau_peak = int((q.idxmax()-pr[pr > 0].index[0]).seconds/60/60)

# # textstr = '\n'.join([r'$\tau_{peak}: 12h$',
# #                      r'$Rain Duration: 25h$',
# #                      r'$Q_{max}: 71.35 m^3\cdot s^{-1}$',
# #                      r'$PR_{max}: 2.55 mm\cdot h^{-1}$',
# #                      r'$PR_{cum}: 24mm$',
# #                      r'$Snowmelt^{ROS}_{cum}: 1.57mm$'])

# # textstr2 = '\n'.join([r'$FloodVolume: 1.67hm^3$',
# #                      r'$RainVolume: 21.6hm^3$',
# #                       r'$MeltedVolume: 1.32hm^3$'])

# # these are matplotlib.patch.Patch properties
# # props = dict(facecolor='teal', alpha=0.1, linewidth=0)

# # place a text box in upper left in axes coords
# # ax[0].text(0.6, 1.2, textstr, transform=ax[0].transAxes, fontsize=12,
# #            verticalalignment='top', bbox=props)
# # ax[0].text(0.886, 1.2, textstr2, transform=ax[0].transAxes, fontsize=12,
# #            verticalalignment='top', bbox=props)

# # plt.savefig('plots/caseofstudy_Aug2013/flood_study.pdf', dpi=150,
# #             bbox_inches='tight')
