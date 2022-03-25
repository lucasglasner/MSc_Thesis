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

# =============================================================================
# Load data
# =============================================================================

# =============================================================================
# Basin hypsometry
# =============================================================================
cuenca = "RioMaipoEnElManzano"
curva_hipso = pd.read_csv("datos/topography/basins/hypso/"+cuenca+"_hypso.csv")
curva_hipso.drop_duplicates(subset="Area_km2", inplace=True)
basin = gpd.read_file("datos/vector/basins/RioMaipoEnElManzano.shp")

# =============================================================================
# #load data from sto domingo radiosonde
# =============================================================================
H0_stodomingo = pd.read_csv("datos/stodomingo/isoterma0.csv", index_col=0)
H0_stodomingo.index = pd.to_datetime(H0_stodomingo.index)
H0_stodomingo = H0_stodomingo.where(
    H0_stodomingo < 7e3).dropna()["H0_StoDomingo"]

# =============================================================================
# #load data from AMDAR database
# =============================================================================
H0_amdar = pd.read_csv("datos/amdar/isoterma0SCEL20172019.csv", index_col=0)


def func(x): return datetime.fromordinal(int(x)) + \
    timedelta(days=x % 1) - timedelta(days=366)


H0_amdar.index = H0_amdar["timeregh"].map(func)
H0_amdar = H0_amdar["zt0"]
H0_amdar.index = pd.date_range(
    "2017-01-01T01:00:00", "2020-01-01", freq="h").dropna()
# H0_amdar.index = H0_amdar.index.map(lambda x: x+timedelta(hours=4))

# =============================================================================
# #load precipitation data
# =============================================================================
pr_qn = pd.read_csv("datos/estaciones/pr_quintanormal.csv", dtype=str)
pr_qn.index = pr_qn.iloc[:, 0]+"-"+pr_qn.iloc[:, 1]+"-"+pr_qn.iloc[:, 2]
pr_qn.index = pd.to_datetime(pr_qn.index)
pr_qn = pd.to_numeric(pr_qn.drop(pr_qn.columns[[0, 1, 2]], axis=1).iloc[:, 0])

pr_lo = pd.read_csv("datos/estaciones/pr_laobra.csv", dtype=str)
pr_lo.index = pr_lo.iloc[:, 0]+"-"+pr_lo.iloc[:, 1]+"-"+pr_lo.iloc[:, 2]
pr_lo.index = pd.to_datetime(pr_lo.index)
pr_lo = pd.to_numeric(pr_lo.drop(pr_lo.columns[[0, 1, 2]], axis=1).iloc[:, 0])

data_dgf = pd.read_csv(
    "datos/estaciones/dgf/DATOSUTC_2004-2019.csv", index_col=0)
data_dgf.index = pd.to_datetime(data_dgf.index.values)

pr_cr2met = pd.read_csv("datos/cr2met/pr_RioMaipoEnElManzano.csv")
pr_cr2met.index = pd.to_datetime(pr_cr2met["date"])
pr_cr2met = pr_cr2met["5710001"]

precip_days = pr_cr2met > 5
precip_days = precip_days[precip_days].index


# =============================================================================
# load era5 zero degree level
# =============================================================================

H0_ERA5_stodom = pd.read_csv(
    "datos/era5/H0_ERA5_stodomingo.csv", index_col=0).squeeze()
H0_ERA5_stodom.index = pd.to_datetime(H0_ERA5_stodom.index)
H0_ERA5_stodom = H0_ERA5_stodom[~H0_ERA5_stodom.index.duplicated()].dropna()
# %%

# =============================================================================
# Build freezing level height for maipo en el manzano basin based upon a simple
# linear regression model between elevation and temperature of each day
# =============================================================================

try:
    H0_cr2met = pd.read_csv(
        "datos/cr2met/freezinglevel_t2m_linregress.csv", index_col=0)
    H0_cr2met.index = pd.to_datetime(H0_cr2met.index)
    H0_era5land = pd.read_csv(
        "datos/era5land/RioMaipoEnElManzano/freezinglevel_t2m_linregress.csv", index_col=0)
    H0_era5land.index = pd.to_datetime(H0_era5land.index)

    lm_cr2met = pd.read_csv("datos/cr2met/t2m_linregress.csv", index_col=0)
    lm_cr2met.index = pd.to_datetime(lm_cr2met.index)
    lm_era5land = pd.read_csv(
        "datos/era5land/RioMaipoEnElManzano/t2m_linregress.csv", index_col=0)
    lm_era5land.index = pd.to_datetime(lm_era5land.index)

    H0_cr2met = H0_cr2met[lm_cr2met.rsquared > 0.8]
    H0_era5land = H0_era5land[lm_era5land.rsquared > 0.8]
except:
    # =============================================================================
    #
    # #load cr2met terrain and temperature data
    # =============================================================================
    dem = xr.open_dataset(
        "datos/topography/basins/RioMaipoEnElManzano_CR2MET.nc").Band1
    t2m_cr2met = xr.open_dataset("datos/cr2met/CR2MET_t2m_v2.0_day_1979_2020_005deg_RioMaipoEnElManzano.nc",
                                 chunks="auto").t2m

    # =============================================================================
    # load era5 land temperature and terrain data
    # =============================================================================

    dem2 = xr.open_dataset(
        "datos/topography/basins/RioMaipoEnElManzano_ERA5LAND.nc").z
    t2m_era5land = xr.open_dataset(
        "datos/era5land/RioMaipoEnElManzano/2m_temperature.nc", chunks="auto").t2m

    x = dem.values.ravel()
    x1 = dem2.values.ravel()
    linmodels = []
    linmodels2 = []
    for i, date in enumerate(t2m_cr2met.time):
        y = t2m_cr2met.sel(time=date).values.ravel()
        mask = ~np.isnan(x) & ~np.isnan(y)
        m = linregress(x[mask], y[mask])
        linmodels.append(m)

    for i, date in enumerate(t2m_era5land.time):
        y1 = t2m_era5land.sel(time=date).values.ravel()-273.15
        mask = ~np.isnan(x1) & ~np.isnan(y1)
        m1 = linregress(x1[mask], y1[mask])
        linmodels2.append(m1)

        lm_cr2met = pd.DataFrame([(m.slope, m.intercept, m.rvalue**2) for m in linmodels],
                                 index=t2m_cr2met.time.values, columns=["slope", "intercept", "rsquared"])
        lm_era5land = pd.DataFrame([(m.slope, m.intercept, m.rvalue**2) for m in linmodels2],
                                   index=t2m_era5land.time.values, columns=["slope", "intercept", "rsquared"])

        H0_cr2met = [-m.intercept/m.slope for m in linmodels]
        H0_cr2met = pd.Series(H0_cr2met, index=t2m_cr2met.time.values)

        H0_era5land = [-m1.intercept/m1.slope for m1 in linmodels2]
        H0_era5land = pd.Series(H0_era5land, index=t2m_era5land.time.values)
        del mask, m, linmodels, linmodels2, x, x1, y, y1
# %%
# =============================================================================
# Define a function for calculating the number of pixels with T<0 below the
# freezing level and the number of pixels with T>0 above the freezing level
# =============================================================================


def level_error(timeseries, raster, dem, field_limit=0):
    timeseries = timeseries.reindex(raster.time.values).dropna()
    raster = raster.sel(time=timeseries.index)
    tpix = np.count_nonzero(~np.isnan(dem))
    metric = np.empty(len(timeseries))
    for i, date in enumerate(timeseries.index):
        height = timeseries.loc[date]
        mask = dem < height
        below = raster[i, :, :].where(mask) < field_limit
        above = raster[i, :, :].where(~mask) > field_limit
        metric[i] = (np.count_nonzero(below)+np.count_nonzero(above))/(tpix)
    return pd.Series(metric, index=timeseries.index)

# %%

# =============================================================================
# compute the metric for each method/dataset
# =============================================================================

# metric_stodomingo = level_error(H0_stodomingo,t2m_cr2met,dem)
# metric_amdar      = level_error(H0_amdar,t2m_cr2met,dem)
# metric_cr2met_OLR = level_error(H0_cr2met,t2m_cr2met,dem)
# metric_era5       = level_error(H0_ERA5_stodom,t2m_cr2met,dem)


# %%
# =============================================================================
# Read file with freezing temperature by band or make it.
# =============================================================================
# precip_days=precip_days[:100]
try:
    flevel_bands = pd.read_csv("datos/cr2met/freezinglevel_t2m_bands.csv",
                               index_col=0)
    flevel_bands.columns = pd.to_datetime(flevel_bands.columns)
except:
    # =========================================================================
    # Build elevation band masks
    # =========================================================================
    dz = 200
    elevation_bands = np.arange(dem.min(), dem.max()+dz, dz)

    masks = []
    for j in range(len(elevation_bands)-1):
        z0 = elevation_bands[j]
        z1 = elevation_bands[j+1]
        mask = (dem.where((dem > z0) & (dem < z1)) > 0).values
        masks.append(mask)

    elevation_bands = elevation_bands[:-1]
    # =========================================================================
    # Apply band mask, and compute % of pixels with temperature < 0°C by band
    # =========================================================================

    flevel_bands = np.empty((len(elevation_bands), len(t2m_cr2met.time)))
    flevel_bands = pd.DataFrame(flevel_bands,
                                index=elevation_bands,
                                columns=t2m_cr2met.time.values)
    for i in trange(len(t2m_cr2met.time.values)):
        # tile_date = precip_days[i]
        tile_date = t2m_cr2met.time.values[i]
        for j in range(len(masks)):
            band = masks[j]
            tile_band = t2m_cr2met.sel(time=tile_date).where(band)
            freezing_band = np.where((tile_band < 0), 1, 0)
            flevel_bands.iloc[j, i] = freezing_band.sum() / \
                np.count_nonzero(band)
    flevel_bands.to_csv("datos/cr2met/freezinglevel_t2m_bands.csv")

# %%

var = flevel_bands  # .loc[fSCA_bands.index<4.5e3,:]
H50 = np.ones((var.shape[1]))*np.nan
H20 = np.ones((var.shape[1]))*np.nan
H80 = np.ones((var.shape[1]))*np.nan

for i in trange(var.shape[1]):
    idx1 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.2*var.max(axis=0)[i]))).flatten()
    idx2 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.5*var.max(axis=0)[i]))).flatten()
    idx3 = np.argwhere(
        np.diff(np.sign(var.iloc[:, i]-0.8*var.max(axis=0)[i]))).flatten()
    if ((len(idx1) != 0) & (len(idx2) != 0) & (len(idx2) != 0)):
        interp1 = interp1d([var.iloc[idx1[0], i], var.iloc[idx1[0]+1, i]],
                           [var.index[idx1[0]], var.index[idx1[0]+1]])
        interp2 = interp1d([var.iloc[idx2[0], i], var.iloc[idx2[0]+1, i]],
                           [var.index[idx2[0]], var.index[idx2[0]+1]])
        interp3 = interp1d([var.iloc[idx3[0], i], var.iloc[idx3[0]+1, i]],
                           [var.index[idx3[0]], var.index[idx3[0]+1]])
        try:
            H20[i] = interp1(0.2)
        except:
            H20[i] = np.nan
        try:
            H50[i] = interp2(0.5)
        except:
            H50[i] = np.nan
        try:
            H80[i] = interp3(0.8)
        except:
            H80[i] = np.nan
H20 = pd.Series(H20, index=flevel_bands.columns)
H50 = pd.Series(H50, index=flevel_bands.columns)
H80 = pd.Series(H80, index=flevel_bands.columns)


# %%
timerange = pd.date_range("1979-01-01T00:00:00",
                          "2021-12-31T00:00:00", freq="h")
isotermas0 = [data.squeeze().reindex(timerange) for data in [
    H0_stodomingo, H0_amdar, H0_cr2met, H0_era5land, H0_ERA5_stodom, H20, H50, H80]]
isotermas0 = pd.concat(isotermas0, axis=1)

isotermas0.columns = ["stodomingo", "amdar", "cr2met", "era5land",
                      "era5_NNstodomingo", "H20_cr2met", "H50_cr2met", "H80_cr2met"]
isotermas0 = isotermas0[isotermas0 < 10e3]
isotermas0 = isotermas0[isotermas0 > 0]

# %%

rainy_mask = pd.Series(isotermas0.index.date).map(
    lambda x: x in precip_days.date).values
mask = pd.Series(H0_stodomingo.index.date).map(
    lambda x: x in precip_days.date).values
# %%


# %%
fig, ax = plt.subplots(2, 5, sharex="row", sharey="row", figsize=(12, 6))
# plt.scatter(H0_amdar,H0_cr2met.reindex(H0_amdar.index))

ax = ax.ravel()
var = isotermas0
for i in range(len(ax)-5):
    ax[i].scatter(var.iloc[:, 0], var.iloc[:, i], alpha=0.5)
    ax[i].scatter(var.iloc[:, 0][rainy_mask],
                  var.iloc[:, i][rainy_mask], alpha=0.5)
    ax[i].plot([0, 8e3], [0, 8e3], "k--")
    # ax[i].set_xlim(0,8.5e3)
    # ax[i].set_ylim(0,8.5e3)
    ax[i].set_aspect('equal')
    ax[i].set_yticks(np.arange(0, 8.5e3+1e3, 1e3))
    ax[i].set_xticks(np.arange(0, 8.5e3+1e3, 1e3)[::2])
    ax[i].grid(True, ls=":")
    ax[i].set_title(var.iloc[:, i].name)

    ax[i+5].boxplot(var.iloc[:, i].dropna(), positions=[0],
                    sym="", patch_artist=True, medianprops={"color": "k"})
    ax[i+5].boxplot(var.iloc[:, i][rainy_mask].dropna(), positions=[0.5], sym="", patch_artist=True,
                    boxprops={"facecolor": "tab:orange"}, medianprops={"color": "k"})
    ax[i+5].grid(True, ls=":")
    ax[i+5].set_xticks([])
    ax[i+5].set_xticklabels([])
    ax[i+5].set_yticks(np.arange(0, 8.5e3+1e3, 1e3))


ax[0].scatter([], [], color="tab:orange", label="Pr>5mm")
ax[0].legend(frameon=False)
ax[2].plot([], [], color="k", ls="--", label="y~x")
ax[2].legend(frameon=False)

ax[2].set_xlabel("H0_StoDomingo (m)", fontsize=12)
# ax[1].set_title("H0_Amdar (m)")
# ax[2].set_title("H0_CR2MET_LR (m)")
# ax[3].set_title("H0_ERA5_stodom-Land (m)")
# ax[4].set_title("H0_ERA5_stodom (m)")

# plt.savefig("plots/maipomanzano/datasetcomparison/isotherm0scatters.pdf",dpi=150,bbox_inches="tight")


# %%
fig = plt.figure(num=0, figsize=(8, 6), dpi=150)
fig.tight_layout(pad=2)
ax = fig.add_subplot(211)
ax1 = fig.add_subplot(223)
ax2 = fig.add_subplot(224)
fig.tight_layout(pad=3)
ax.plot(H0_stodomingo, color="powderblue")
ax.plot(H0_amdar, color="wheat")
ax.set_ylabel("Isotherm 0°C height (m)")


hamd = H0_amdar.resample("12h").interpolate()
hsto = H0_stodomingo.resample("12h").asfreq().reindex(hamd.index)
mask = pr_cr2met.reindex(hamd.index) > 5

ax1.boxplot(hsto.dropna(), sym="", positions=[0], showmeans=True, meanline=True,
            patch_artist=True,
            boxprops={"facecolor": "powderblue"},
            meanprops={"linestyle": "--", "color": "k"},
            medianprops={"color": "tab:red"})

ax1.boxplot(hsto[mask].dropna(),
            sym="", positions=[0.75], showmeans=True, meanline=True,
            patch_artist=True,
            boxprops={"facecolor": "powderblue"},
            meanprops={"linestyle": "--", "color": "k"},
            medianprops={"color": "tab:red"})

ax1.boxplot(hamd.dropna(), sym="", positions=[0.25], showmeans=True, meanline=True,
            patch_artist=True,
            boxprops={"facecolor": "wheat"},
            meanprops={"linestyle": "--", "color": "k"},
            medianprops={"color": "tab:red"})


ax1.boxplot(hamd[mask].dropna(),
            sym="", positions=[1], showmeans=True, meanline=True,
            patch_artist=True,
            boxprops={"facecolor": "wheat"},
            meanprops={"linestyle": "--", "color": "k"},
            medianprops={"color": "tab:red"})
ax1.scatter([], [], color="powderblue", marker="s",
            label="Sto. Domingo\nRadiosonde", edgecolor="k")
ax1.scatter([], [], color="wheat", marker="s", label="AMDAR", edgecolor="k")
ax1.legend(frameon=False, ncol=2, loc=(0, 1))
ax1.grid(True, axis="y", ls=":")
ax1.set_xticks([0.25/2, (0.75)+0.25/2])
ax1.set_xticklabels(["H0", "H0\nPrecip.>1mm"])
ax1.set_ylabel("Isotherm 0°C height (m)")

ax2.plot([1e1, 7e3], [1e1, 7e3], color="tab:red", ls="--", zorder=1)
ax2.scatter(x=hamd[mask].values,
            y=hsto[mask].values,
            edgecolor="k", zorder=2, alpha=.9, color="gold", s=20, label="Precip.>1mm")
ax2.scatter(x=hamd.values,
            y=hsto.values,
            edgecolor="k", zorder=1, alpha=.5, color="tab:blue", s=20)
ax2.set_ylabel("Sto. Domingo")
ax2.set_xlabel("AMDAR")
ax2.legend(frameon=False)


ax3 = fig.add_subplot(3, 9, 27)
diff = hamd[mask]-hsto[mask]
diff = diff.dropna()
ax3.boxplot(diff, sym="", positions=[0], patch_artist=True,
            boxprops={"facecolor": "gold"},
            medianprops={"color": "k"})
ax3.set_title("$\Delta H0$", loc="center")
ax3.set_xticks([])
ax3.set_xticklabels([])
ax3.set_yticks(np.arange(-150, 750, 150))
# ax3.grid(True,axis="y",ls=":")
ax3.yaxis.tick_right()
plt.savefig("plots/maipomanzano/datasetcomparison/isotherm0_study.pdf",
            dpi=150, bbox_inches="tight")
