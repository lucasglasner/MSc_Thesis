#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:27:42 2022

@author: lucas

# =============================================================================
# This Script does some hydrological analysis in the case of study
# =============================================================================
"""

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
# =============================================================================
# big time interval and graph time interval
# =============================================================================
interval = slice(datetime.datetime(2013, 8, 4),
                 datetime.datetime(2013, 8, 16))
interval2 = slice(datetime.datetime(2013, 8, 11),
                  datetime.datetime(2013, 8, 14))


# %%
# =============================================================================
# SNOW LIMIT AND SNOW COVER
# =============================================================================

SL_mm = pd.read_csv('datos/snowlimits_maipomanzano.csv',
                    index_col=0).squeeze()
SL_mm.index = pd.to_datetime(SL_mm.index)+datetime.timedelta(hours=12)
SL_mm = SL_mm.dropna(how='all')[interval]

for i in range(3):
    y = SL_mm.iloc[:, i].dropna()
    x = SL_mm["IANIGLA"].dropna().reindex(y.index)
    m = st.linregress(x, y)
    for j in range(len(SL_mm)):
        if np.isnan(SL_mm.iloc[j, i]):
            SL_mm.iloc[j, i] = m.slope * \
                SL_mm["IANIGLA"].values[j]+m.intercept

SL_mm = SL_mm["MODIS_H50"]

SCA = pd.read_csv('datos/snowcovers_maipomanzano.csv',
                  index_col=0)['IANIGLA']
SCA.index = pd.to_datetime(SCA.index)
SCA = SCA[interval]

# %%
# =============================================================================
# DGF ROOF STATIION DATA (SANTIAGO CITY)
# =============================================================================
datos_dgf = pd.read_csv(
    "datos/estaciones/dgf/DATOSUTC_2004-2019.csv", index_col=0)
datos_dgf.index = pd.to_datetime(
    datos_dgf.index.values)-datetime.timedelta(hours=4)
datos_dgf = datos_dgf[interval]

# %%
# =============================================================================
# CR2MET PRECIPITATION BASIN WIDE MEAN
# =============================================================================

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv", index_col=0)
pr_cr2met.index = pd.to_datetime(
    pr_cr2met["date"])+datetime.timedelta(hours=12)
pr_cr2met.drop("date", axis=1, inplace=True)
pr_cr2met = pr_cr2met.squeeze()[interval]


# %%
# =============================================================================
# hypsometric curves
# =============================================================================


hypso = pd.read_csv(
    'datos/topography/basins/hypso/RioMaipoEnElManzano_hypso.csv')
area = hypso.Area_km2.values[-1]

int_func = interp1d(hypso.iloc[:, 0], hypso.iloc[:, 1])

# %%
# =============================================================================
# ZERRO DEGREE LEVEL AND PLUVIAL AREA
# =============================================================================

H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)-datetime.timedelta(hours=4)
H0_mm = H0_mm[interval]

pluv_area = H0_mm.map(lambda x: int_func(x))*area


# %%
# =============================================================================
# PRECIPITATION IN MAIPO EN EL MANZANO BASIN OUTLET
# =============================================================================


pr_mm = pd.read_csv('datos/estaciones/pr_RioMaipoEnElManzano_2013-08.csv')
pr_mm.index = pd.to_datetime(pr_mm['Fecha'])
pr_mm = pr_mm['Valor'].drop_duplicates()
pr_mm = pr_mm.reindex(pd.date_range(interval.start,
                                    interval.stop, freq='h')).fillna(0)


# %%
# =============================================================================
# RUNOFF DATA IN DIFFERENT BASINS
# =============================================================================
qinst_mm = pd.read_csv("datos/estaciones/qinst_RioMaipoEnElManzano.csv",
                       index_col=0)
qinst_mm.index = pd.to_datetime(qinst_mm.index)
qinst_mm = qinst_mm.squeeze()[interval]

runoff = pd.read_csv('datos/runoff_gauges_dataset.csv', index_col=0)
runoff.index = pd.to_datetime(runoff.index)
runoff = runoff[interval]

# %%
# =============================================================================
# BASIN POLYGONS
# =============================================================================
paths = glob('datos/vector/basins/*.shp')
polygons = pd.concat([gpd.read_file(p) for p in paths])


# %%
# =============================================================================
# RASTER DATA, SWE AND TOPOGRAPHY
# =============================================================================
dem = xr.open_dataset(
    'datos/topography/basins/RioMaipoEnElManzano_Cortes.nc').Band1
SWE = xr.open_dataset(
    'datos/ANDES_SWE_Cortes/maipomanzano/ANDES_SWE_WY2014.nc')
SWE = SWE.SWE.sel(time=interval)
dSWE = SWE.diff('time')

masks = []
# Pluvial area masks
pluv_area_daily = pluv_area.resample('d').mean()
for i, pa in enumerate(pluv_area_daily):
    if np.isnan(pa):
        mask = dem*0
        masks.append(mask)
    else:
        mask = xr.where((dem < pa), 1, 0)
        masks.append(mask)

masks = xr.concat(masks, pluv_area_daily.index).rename(
    {'timestamp_UTC': 'time'})

# %%
# =============================================================================
# Build SWE loss/gain series
# =============================================================================
parea = (6.4e6)**2*np.cos(np.deg2rad(SWE.lat.mean()))
parea = parea*np.deg2rad(0.001)*np.deg2rad(0.001)
parea = parea.item()

melt = dSWE.where(dSWE < 0).mean(dim=['lat', 'lon']).to_series()
gain = dSWE[:-1, :, :].where(masks[1:, :, :])
gain = gain.where(gain > 0).mean(dim=['lat', 'lon']).to_series()

# %%
# =============================================================================
# Build table for document
# =============================================================================
datos = []
for dat in [SL_mm, SCA, melt*-1, H0_mm,
            datos_dgf.iloc[:, 9].resample('d').sum(),
            pr_mm.resample('d').sum(),
            qinst_mm.resample('d').max(),
            datos_dgf.iloc[:, 5].resample('d').mean(),
            datos_dgf.iloc[:, 5].resample('d').max(),
            datos_dgf.iloc[:, 5].resample('d').min()]:
    datos.append(dat[interval])

datos = pd.concat(datos, axis=1)
datos = datos.resample('d').mean().iloc[:-1, :]
datos.columns = ["SL", "SCA", "MELT", "H0", "PR_DGF", "PR_MM", "Qmax", "T",
                 "Tmax", "Tmin"]
datos = np.round(datos, 1)
datos = datos["2013-08-03":"2013-08-16"]

# %%
# =============================================================================
# flood data
# =============================================================================


pr = pr_mm[interval]
q = qinst_mm[interval]
ap = pluv_area.reindex(q.index).interpolate(method='cubicspline')


snow_area = SCA.reindex(pluv_area.index, method='nearest')/100*area
snow_area = snow_area.reindex(q.index).interpolate(method='cubicspline')
ros_area = ap-(area-snow_area)
# ros_area = ros_area/ap
mlt = -1*melt.reindex(q.index).interpolate(method='cubicspline')/24
gin = gain.reindex(q.index).interpolate(method='cubicspline')/24

baseflow = sliding_interval_filter(q, 40)[0]

mask = pr > 0.01


pr = pr[interval2]
q = q[interval2]
ap = ap[interval2]
snow_area = snow_area[interval2]
ros_area = ros_area[interval2]
mlt = mlt[interval2]
gin = gin[interval2]
baseflow = baseflow[interval2]
mask = mask[interval2]

rainpulse = slice('2013-08-11T14:00:00', '2013-08-12T15:00:00')
floodpulse = slice('2013-08-11T18:00:00', '2013-08-13T09:00:00')

# %%
plt.rc("font", size=18)
fig, ax = plt.subplots(2, 1, figsize=(9, 4), sharex=True)

ax[0].bar(pr.index-datetime.timedelta(hours=1),
          pr, width=0.037, edgecolor="k", color='cadetblue',
          align='edge', label='Precipitation')
ax[0].set_yticks([0, 1, 2, 3])
ax[0].set_ylim(0, 4)
# ax[0].bar(melt.index, melt*-1, edgecolor="k")
ax[0].bar(pr.index-datetime.timedelta(hours=1),
          mlt,
          width=0.037, edgecolor='k', color='violet',
          align='edge', label='Snowmelt', bottom=pr)
# ax[0].bar(pr.index-datetime.timedelta(hours=1),
#             gin,
#             width=0.037, edgecolor='k', color='orange')
ax[0].set_ylabel('(mm/h)')
ax[0].legend(frameon=False, loc=(0.01, 0.95), ncol=2, fontsize=12)
# ax[0].set_yscale("log")


ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)

ax[1].plot(pr.index, q-baseflow, label='Direct\nRunoff')
ax[1].set_yticks([0, 7, 14, 21, 28])
ax[1].set_ylim(0, 28)
# ax[1].scatter((q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]].index,
#               (q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]],
#               ec="k",zorder=10)
ax[1].axvspan("2013-08-11T18:00:00",
              "2013-08-13T09:00:00", alpha=0.15, color='k')
# ax[1].axvline("2013-08-11T14:00",color='k',alpha=0.5, ls=":")
ax[1].set_ylabel('$(m^3/s)$')
ax[1].legend(loc='upper left', frameon=False, fontsize=12)


ax1 = ax[1].twinx()
ax1.set_ylim(0, 1)
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax1.plot(pr.index, ap/area, color='tab:red', alpha=0.2, ls="--")
ax1.plot(pr.index, ap.where(pr > 0)/area, color='tab:red',
         label='Pluvial Area')
ax1.plot(pr.index, ros_area.where(ros_area > 0).where(pr > 0)/area,
         color='purple', label='ROS Area')
ax1.set_ylabel('Fraction of\ntotal Area (%)')
ax1.legend(frameon=False, fontsize=12, loc='upper right')
# ax1.set_yticklabels(ax1.get_yticks()*100)
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%i'))
# ax[1].xaxis.set_major_formatter(
#       mpl.dates.ConciseDateFormatter(ax[1].xaxis.get_major_locator()))


# dtFmt =  # define the formatting
ax[1].xaxis.set_major_formatter(mpl.dates.DateFormatter('\n\n%b-%d'))
ax[1].xaxis.set_major_locator(mpl.dates.DayLocator(interval=1))

ax[1].xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
ax[1].xaxis.set_minor_locator(
    mpl.dates.HourLocator(byhour=np.arange(0, 24, 3)))
ax[1].tick_params(axis='x', which='minor', rotation=45)

for maj in ax[1].xaxis.get_major_ticks():
    maj.label.set_fontsize(18)
for m in ax[1].xaxis.get_minor_ticks():
    m.label.set_fontsize(12)
# plt.show()

ax[1].set_xlim([15928.4, 15931])

tau_peak = int((q.idxmax()-pr[pr > 0].index[0]).seconds/60/60)

textstr = '\n'.join([r'$\tau_{peak}: 12h$',
                     r'$Rain Duration: 25h$',
                     r'$Q_{max}: 71.35 m^3\cdot s^{-1}$',
                     r'$PR_{max}: 2.55 mm\cdot h^{-1}$',
                     r'$PR_{cum}: 24mm$',
                     r'$Snowmelt^{ROS}_{cum}: 1.57mm$'])

textstr2 = '\n'.join([r'$FloodVolume: 1.67hm^3$',
                     r'$RainVolume: 21.6hm^3$',
                      r'$MeltedVolume: 1.32hm^3$'])

# these are matplotlib.patch.Patch properties
props = dict(facecolor='teal', alpha=0.1, linewidth=0)

# place a text box in upper left in axes coords
ax[0].text(0.6, 1.2, textstr, transform=ax[0].transAxes, fontsize=12,
           verticalalignment='top', bbox=props)
ax[0].text(0.886, 1.2, textstr2, transform=ax[0].transAxes, fontsize=12,
           verticalalignment='top', bbox=props)

plt.savefig('plots/caseofstudy_Aug2013/flood_study.pdf', dpi=150,
            bbox_inches='tight')

# %%

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-72.3, -69.5, -32.4, -37])
ax.coastlines()
ax.add_feature(cf.BORDERS, rasterized=True)
ax.add_feature(cf.OCEAN, rasterized=True)
ax.add_feature(cf.LAND, rasterized=True)
polygons.plot(polygons.gauge_name, ax=ax, transform=ccrs.PlateCarree(),
              cmap='tab20')
polygons.boundary.plot(ax=ax, lw=1, color='k', transform=ccrs.PlateCarree())
