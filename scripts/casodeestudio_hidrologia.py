#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:27:42 2022

@author: lucas

# =============================================================================
# This Script does some hydrological analysis in the case of study
# =============================================================================
"""

import cartopy.feature as cf
import geopandas as gpd
from scipy.signal import argrelextrema
from scipy.ndimage.measurements import label
from scipy.ndimage.filters import minimum_filter1d, generic_filter
import scipy.stats as st
from scipy.interpolate import interp1d
import datetime
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from glob import glob
from functions import add_labels, local_minimum_filter
import sys
sys.path.append('functions.py')

# %%


# %%
# =============================================================================
# big time interval and graph time interval
# =============================================================================
date = "2008-06-04"
# date = "%YR%"
yr, month, day = [int(n) for n in date.split("-")]
interval = slice(datetime.datetime(yr, month, day)-datetime.timedelta(days=9),
                 datetime.datetime(yr, month, day)+datetime.timedelta(days=3))

date_interval = pd.date_range(interval.start, interval.stop, freq='h')
basin_attributes = pd.read_csv('datos/basins_attributes.csv')
basin_attributes.index = basin_attributes.gauge_name
# %%
# =============================================================================
# PRECIPITATION ON BASIN OUTLETS
# =============================================================================

basins = ['Rio Maipo En El Manzano', 'Rio Teno Despues De Junta Con Claro',
          'Rio Uble En San Fabian N 2']

# pr_mm = pd.read_csv('datos/estaciones/pr_RioMaipoEnElManzano_2013-08.csv')
# pr_mm.index = pd.to_datetime(pr_mm['Fecha'])
# pr_mm = pr_mm['Valor'].drop_duplicates()
# pr_mm = pr_mm.reindex(pd.date_range(interval.start,
#                                     interval.stop, freq='h')).fillna(0)

# pr_teno = pd.read_csv(
#     'datos/estaciones/pr_RioTenoDespuesDeJuntaConClaro_2013-08.csv')
# pr_teno.index = pd.to_datetime(pr_teno.Fecha)
# pr_teno = pr_teno['Valor'].reindex(pr_mm.index)[interval]


# pr_uble = pd.read_csv('datos/estaciones/pr_RioUbleEnSanFabian_2013-08.csv')
# pr_uble.index = pd.to_datetime(pr_uble.Fecha)
# pr_uble = pr_uble['Valor'].reindex(pr_mm.index)[interval]


# pr = pd.concat([pr_mm, pr_teno, pr_uble], axis=1)
# pr.columns = columns = ['Rio Maipo En El Manzano',
#                         'Rio Teno Despues De Junta Con Claro',
#                         'Rio Uble En San Fabian N 2']
# del pr_mm, pr_teno, pr_uble

pr = []
for basin in basins:
    pp = pd.read_csv('datos/estaciones/station_' +
                     basin.replace(" ", "")+'_2000-01-21_2021-01-01.csv', index_col=0)
    pp.index = pd.to_datetime(pp['Fecha'])
    pp = pp['Valor'].drop_duplicates()
    pr.append(pp)

pr = pd.concat(pr, axis=1).reindex(date_interval)
pr.columns = basins
pr = pr[interval]

pr_laobra = pd.read_csv('datos/estaciones/pr_laobra.csv', dtype=str)
pr_laobra.index = pd.to_datetime(
    pr_laobra.iloc[:, 0]+"-"+pr_laobra.iloc[:, 1]+"-"+pr_laobra.iloc[:, 2])
pr_laobra = pd.to_numeric(pr_laobra.iloc[:, 3])
# %%
# =============================================================================
# hypsometric curves
# =============================================================================

hypso = [pd.read_csv('datos/topography/basins/hypso/' +
                     b.replace(" ", "")+"_hypso.csv") for b in basins]

hypso = pd.concat(hypso, keys=pr.columns, axis=1)

areas = [hypso[b].Area_km2.max() for b in pr.columns]
areas = pd.Series(areas, index=pr.columns)
int_func = [interp1d(hypso[b].height, hypso[b].fArea) for b in pr.columns]

# %%
# =============================================================================
# ZERO DEGREE LEVEL AND PLUVIAL AREA
# =============================================================================

# santo domingo
H0_mm = pd.read_csv('datos/stodomingo/isoterma0.csv',
                    index_col=0).squeeze()
H0_mm.index = pd.to_datetime(H0_mm.index)-datetime.timedelta(hours=4)
H0_mm = H0_mm[date[:-5]].resample("h").interpolate('cubicspline')
H0_mm = H0_mm[interval]

# era 5
H0 = xr.open_dataset('datos/era5/H0_ERA5_'+str(yr)+'.nc').deg0l
H0_m = []
for b in basins:
    lat, lon = basin_attributes.loc[b].gauge_lat, basin_attributes.loc[b].gauge_lon-1
    h0 = H0.sel(lat=lat, lon=lon, method='nearest').to_series()
    H0_m.append(h0[date[:-5]])

H0_m = pd.concat(H0_m, axis=1)
H0_m.columns = basins
H0_m = H0_m.resample('h').interpolate('cubicspline')[interval]

del H0

#%%
# area pluvial con santo domingo
pluv_area = [int_func[i](H0_mm-300)*areas[i]
             for i in range(pr.shape[1])]
nonpluv_area = [(1-int_func[i](H0_mm-300))*areas[i]
                for i in range(pr.shape[1])]
pluv_area = pd.DataFrame(pluv_area, columns=pr.index, index=pr.columns).T
nonpluv_area = pd.DataFrame(nonpluv_area, columns=pr.index, index=pr.columns).T

# area pluvial con era5
pluv_area_era5 = []
for i in range(H0_m.shape[1]):
    pluv_area_era5.append(int_func[i]((H0_m.iloc[:, i]-300))*areas[i])

pluv_area_era5 = pd.DataFrame(
    pluv_area_era5, index=pluv_area.columns, columns=pluv_area.index).T
nonpluv_area_era5 = areas-pluv_area_era5
# pluv_area_era5.columns = pluv_area.columns
# # %%
# fig, ax = plt.subplots(1, 1)
# (pluv_area).plot(ax=ax)
# (pluv_area_era5).plot(ax=ax, legend=False)


# %%
# =============================================================================
# SNOW COVER
# =============================================================================


SCA = []
for b in basins:
    p = 'datos/ianigla/'+b.replace(" ", "")+'_SCA_s_comp.filtro_MA.3días.csv'
    sca = pd.read_csv(p, index_col=0)
    sca.index = pd.to_datetime(sca.index)
    SCA.append(sca.iloc[:, 0]/100)

SCA = pd.concat(SCA, axis=1)
SCA.columns = pr.columns


SCA = SCA[str(yr)].resample('h').interpolate('cubicspline')[interval]

snow_area = [SCA[b]*areas.loc[b] for b in pr.columns]
snow_area = pd.concat(snow_area, axis=1)

# %%
# =============================================================================
# ROS AREA
# =============================================================================
# ROS = xr.open_dataset('datos/ROS/CORTES_CR2MET_ERA5/ROS_UBLE2013.nc')
# ROS = ROS.ROS.sel(time=interval)

ros_area = np.clip(snow_area-nonpluv_area, 0, 7e3)
# ros_area = ros_area.where(pr > 3/24).where(SCA.diff() <= 0)


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

SWE = [SWE0.SWE.sel(time="2013-08"),
       SWE1.SWE.sel(time="2013-08"),
       SWE2.SWE.sel(time="2013-08")]

del SWE0, SWE1, SWE2

dSWE = [swe.diff('time') for swe in SWE]

tot_pix = np.array([451481, 91890, 93014])
SCA_cortes = [xr.where(swe > 50, 1, 0) for swe in SWE]
SCA_cortes = [sca.sum(dim=['lat', 'lon']).to_series() for sca in SCA_cortes]
SCA_cortes = pd.concat(SCA_cortes, axis=1)
SCA_cortes.columns = SCA.columns


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
interval2 = slice(datetime.datetime(2013, 8, 10),
                  datetime.datetime(2013, 8, 14))

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


runoff3 = runoff2.resample('d').max()
q_ext = np.empty(runoff3.shape[1])
for i in range(len(q_ext)):
    q_ext[i] = np.percentile(runoff3.iloc[:, i].dropna(), 90)
q_ext = pd.Series(q_ext, index=runoff3.columns)

stats = pd.concat([qmax, q_ext, pr_max, pr_cum,
                  melt_cum, rain_duration, tau], axis=1)
stats.columns = ['Qmax', 'Q90', 'PRmax',
                 'PRcum', 'Snowmelt', 'RainDuration', 'tau']

stats = stats.round(2)
# %%
# =============================================================================
# MAIPO EN EL MANZANO
# =============================================================================
# interval3 = slice(datetime.datetime(2006, 8, 7, 6),
                  # datetime.datetime(2006, 8, 13, 3))
plt.rc("font", size=18)
fig, ax = plt.subplots(2, 1, figsize=(12, 4), sharex=True,
                       gridspec_kw={'height_ratios': [1, 2]})

basin = 'Rio Maipo En El Manzano'
props = dict(facecolor='white', lw=1)
names = ['$Q_{max}$: ', '$Q_{90}$: ', '$PR_{max}$: ', '$PR_{cum}$: ', '$Melt_{cum}$: ',
         '$RainDuration$: ', r'$\tau_{peak}$: ']
units = ['$m^3/s$', '$m^3/s$', '$mm/h$', '$mm$', '$mm$', '$hr$', '$hr$']

ax[0].bar(pr.index-datetime.timedelta(hours=1),
          pr[basin], width=pd.Timedelta(hours=1), color='cadetblue',
          align='edge', label='Precipitation')
ax[0].grid(True, which='both', ls=":")
ax[0].set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
ax[0].set_yticklabels([0, "", "", 3, "", "", 6, "", "", 9, "", "", 12])
ax[0].set_ylim(0, 7)
# ax[0].bar(melt.index, melt[basin], edgecolor="k")
# ax[0].bar(melt.index-datetime.timedelta(hours=1),
          # melt[basin], edgecolor='k', linewidth=0.0, color='darkviolet',
          # width=pd.Timedelta(hours=1), align='edge', label='Snowmelt',
          # bottom=pr[basin].fillna(0), zorder=2)
# # ax[0].bar(pr.index-datetime.timedelta(hours=1),
# #             gin,
# #             width=0.037, edgecolor='k', color='orange')
ax[0].set_ylabel('(mm/h)')
ax[0].legend(frameon=False, loc=(0.01, 0.95), ncol=2, fontsize=12)
ax[0].set_title(basin, loc='right')

# box = stats.loc[basin]
# box = ['{:.1f}'.format((round(box[i], 1))) for i in range(len(box))]

# textstr = '\n'.join([n+b+u for n, b, u in zip(names, box, units)])
# ax[0].text(0.4, 1.5, textstr, transform=ax[0].transAxes, fontsize=10,
           # verticalalignment='top', bbox=props)


ax[1].plot(pr.index, quickflows[basin],
           label='Direct Runoff', color='darkblue')
ax[1].grid(axis="x", which='both', ls=":")
# ax[1, 0].set_yticks([0, 7, 14, 21, 28])
# ax[1, 0].set_ylim(0, 28)
# ax[1].scatter((q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]].index,
#               (q-baseflow)[interval2][np.where((q-baseflow)[interval2]<1)[0][:2]],
#               ec="k",zorder=10)
# ax[1, 0].axvspan("2013-08-11T18:00:00",
# "2013-08-13T09:00:00", alpha=0.15, color='k')
# ax[1].axvline("2013-08-11T14:00",color='k',alpha=0.5, ls=":")
ax[1].set_ylabel('$(m^3/s)$')
ax[1].legend(loc=(0.8, 1), frameon=False, fontsize=12)


# ax[1].plot(pr.index, (quickflows)[interval2], color='darkblue')
# ax[1].plot(pr.index, (quickflows)[interval2], color='darkblue')
ax1 = ax[1].twinx()
ax1.set_ylim(0, 1)
ax1.set_yticks(np.arange(0, 1+0.1, 0.1))
ax1.set_yticklabels([0, "", .2, "", .4, "", .6, "", .8, "", 1])
ax1.plot(pr.index, pluv_area[basin]/areas[basin], color='tab:red',
         alpha=0.2, ls="--")
ax1.plot(pr.index, pluv_area[basin].where(pr[basin] > 0)/areas[basin],
         color='tab:red',
         label='Pluvial Area')
ax1.plot(pr.index, ros_area[basin].where(ros_area[basin] > 0).where(pr[basin] > 0)/areas[basin],
         color='tab:green', label='ROS Area')
ax1.plot(pr.index, SCA[basin].where(pr[basin] > 0), color='tab:blue',
         label='SCA')
ax1.plot(pr.index, SCA[basin], color='tab:blue', alpha=0.2,
         ls="--")

ax1.set_ylabel('Fraction of\ntotal Area (%)')
ax1.legend(frameon=False, fontsize=12, loc=(0, 1), ncol=3)


for axis in [ax[1]]:
    axis.set_xlim(interval.start, interval.stop)
    axis.xaxis.set_major_formatter(mpl.dates.DateFormatter('\n\n%d'))
    axis.xaxis.set_major_locator(mpl.dates.DayLocator(interval=1))

    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
    axis.xaxis.set_minor_locator(
        mpl.dates.HourLocator(byhour=np.arange(0, 24, 3)))
    axis.tick_params(axis='x', which='minor', rotation=45)

    for maj in axis.xaxis.get_major_ticks():
        maj.label.set_fontsize(18)
    for m in axis.xaxis.get_minor_ticks():
        m.label.set_fontsize(12)

box = ax[1].get_position()
fig.text(box.xmin, box.ymin*-1.15, '2013-Aug', ha='center', va='center')
# plt.savefig('plots/caseofstudy_Aug2013/flood_study_'+basin.replace(" ", "")+'.pdf',
#             dpi=150, bbox_inches='tight')

# %%
plt.rc("font", size=18)
fig, ax = plt.subplots(2, 3, figsize=(14, 4), sharex='col',
                       gridspec_kw={'height_ratios': [1, 2]},
                       sharey='row')
titles = ['Rio Maipo\nEn El\nManzano',
          'Rio Teno Despues\nDe Junta\nCon Claro',
          'Rio Ñuble En\nSan Fabian\nN 2']

# these are matplotlib.patch.Patch properties
props = dict(facecolor='white', lw=1)
names = ['$Q_{max}$: ', '$Q_{90}$: ', '$PR_{max}$: ', '$PR_{cum}$: ', '$Melt_{cum}$: ',
         '$RainDuration$: ', r'$\tau_{peak}$: ']
units = ['$m^3/s$', '$m^3/s$', '$mm/h$', '$mm$', '$mm$', '$hr$', '$hr$']

# # place a text box in upper left in axes coords
# ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
#         verticalalignment='top', bbox=props)

for i, axis in enumerate(ax[0, :]):
    axis.bar(pr.index-datetime.timedelta(hours=1),
             pr.iloc[:, i], color='cadetblue',
             width=pd.Timedelta(hours=1),
             align='edge', label='Precipitation',
             edgecolor='k', linewidth=0.0, zorder=1)
    # axis.bar(melt.index-datetime.timedelta(hours=1),
             # melt.iloc[:, i], edgecolor='k', linewidth=0.0, color='darkviolet',
             # width=pd.Timedelta(hours=1), align='edge', label='Snowmelt',
             # bottom=pr.iloc[:, i].fillna(0), zorder=2)
    axis.set_yticks(np.arange(0, 9+1.5, 1.5))
    axis.set_yticklabels([0, "", 3, "", 6, "", 9])

    axis.grid(True, which='major', ls=":", zorder=0)
    axis.set_title(titles[i], loc='left')
    if i == 0:
        axis.legend(frameon=False, ncol=1, loc=(0, 0.58), fontsize=9)

    # box = stats.iloc[i, :]
    # box = ['{:.1f}'.format((round(box[i], 1))) for i in range(len(box))]

    # textstr = '\n'.join([n+b+u for n, b, u in zip(names, box, units)])
    # axis.text(0.6, 1.5, textstr, transform=axis.transAxes, fontsize=10,
    #            verticalalignment='top', bbox=props)
    # axis.set_ylim(0,6)


for i, axis in enumerate(ax[1, :]):
    
    axis.set_xlim(interval.start, interval.stop)
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
    axis1.plot(pr.index, SCA.iloc[:, i].where(pr.iloc[:, i] > 0), color='tab:blue',
               label='SCA')
    axis1.plot(pr.index, SCA.iloc[:, i], color='tab:blue', alpha=0.2,
               ls="--")

    if i < 2:
        axis1.set_yticklabels([""]*5)
    axis.grid(axis='x', which='major', ls=":")
    axis.set_yticks([0, 25, 50, 75, 100])
    axis.xaxis.set_major_formatter(mpl.dates.DateFormatter('\n\n%d'))
    axis.xaxis.set_major_locator(mpl.dates.DayLocator(interval=1))

    axis.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%H:%M'))
    axis.xaxis.set_minor_locator(
        mpl.dates.HourLocator(byhour=np.arange(0, 24, 8)))
    axis.tick_params(axis='x', which='minor', rotation=45)

    for maj in axis.xaxis.get_major_ticks():
        maj.label.set_fontsize(18)
    for m in axis.xaxis.get_minor_ticks():
        m.label.set_fontsize(10)
    # axis.set_xlim([15927.7, 15933])
ax[-1, 0].set_ylabel('$(m^3/s)$')
ax[0, 0].set_ylabel('$(mm/h)$')
ax[1, 0].plot([], [], color='tab:green', label='ROS Area')
ax[1, 0].plot([], [], color='tab:red', label='Pluvial Area')
ax[1, 0].plot([], [], color='tab:blue', label='SCA')
ax[1, 0].legend(frameon=False, fontsize=9, loc=(0, 0.98), ncol=3)

axis1.set_ylabel('Fraction of\ntotal Area (-)')


box = ax[1, 0].get_position()
fig.text(box.xmin*1.15, box.ymin*-1.15,
         '\nAug\n2013', ha='center', va='center')

# plt.savefig('plots/caseofstudy_Aug2013/flood_study_stodomingo.pdf',
#             dpi=150, bbox_inches='tight')
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
