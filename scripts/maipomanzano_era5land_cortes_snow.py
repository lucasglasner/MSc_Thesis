#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:42:55 2021

@author: lucas

# =============================================================================
# Maipo en el Manzano Basin. ERA5 Snow analisis and comparison with observations.
# =============================================================================

"""


from tqdm import trange
import datetime
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import regionmask
import geopandas as gpd
import matplotlib as mpl
import datetime
# %%

# ERA5LAND topography on MaipoEnElManzano basin.
topo = xr.open_dataset(
    "datos/topography/basins/RioMaipoEnElManzano_ERA5LAND.nc")
topo = topo.z.squeeze().assign_coords(
    {"lon": (((topo.lon + 180) % 360) - 180)})
topo_cortes = xr.open_dataset(
    "datos/topography/basins/RioMaipoEnElManzano_Cortes.nc").Band1
# BASIN shapefile
basin = gpd.read_file("datos/vector/RioMaipoEnElManzano.shp")

# Basin hypsometry curve
dem = xr.open_dataset(
    "datos/topography/basins/RioMaipoEnElManzano_regridmodis.nc", chunks="auto").Band1
# Elevation band width in meters
dz = 50
tpix = (~np.isnan(dem.values)).sum()
tpix2 = (~np.isnan(topo.values)).sum()
tpix3 = (~np.isnan(topo_cortes.values)).sum()
# elevation bands
elevation_bands = np.arange(dem.min()//10*10-2*dz,
                            dem.max()//10*10+2*dz, dz,
                            dtype=np.int64)

# bands number of pixels
bands_size = np.empty(len(elevation_bands))*np.nan
# hypsometric curve
hypso_curve = np.empty(len(elevation_bands))*np.nan
hypso_curve_ERA5 = np.empty(len(elevation_bands))*np.nan
hypso_curve_cortes = np.empty(len(elevation_bands))*np.nan
# hypsometric density (band number of pixels / total pixels)
hypso_density = np.empty(len(elevation_bands))*np.nan
hypso_density_ERA5 = np.empty(len(elevation_bands))*np.nan
hypso_density_cortes = np.empty(len(elevation_bands))*np.nan
for i, height in enumerate(elevation_bands):
    mask = (~np.isnan(dem.where(dem > height).where(dem < height+dz)))
    mask2 = (~np.isnan(topo.where(topo > height).where(topo < height+dz)))
    mask3 = (~np.isnan(topo_cortes.where(
        topo_cortes > height).where(topo_cortes < height+dz)))
    bands_size[i] = mask.sum()
    hypso_density[i] = bands_size[i]/tpix
    hypso_density_ERA5[i] = mask2.sum()/tpix2
    hypso_density_cortes[i] = mask3.sum()/tpix3
hypso_curve = np.cumsum(hypso_density)
hypso_curve_ERA5 = np.cumsum(hypso_density_ERA5)
hypso_curve_cortes = np.cumsum(hypso_density_cortes)
# %%

SCA_ianigla = pd.read_csv(
    "datos/ianigla/RioMaipoEnElManzano_SCA_s_comp.filtro_MA.3días.csv", index_col=0)
SCA_ianigla = SCA_ianigla.iloc[:, 0]
SCA_ianigla.index = pd.to_datetime(SCA_ianigla.index)

# %%


SCOV_ERA5 = xr.open_dataset(
    "datos/era5land/RioMaipoEnElManzano/snow_cover.nc", chunks="auto")
SCOV_ERA5 = SCOV_ERA5.snowc[23:, :, :][::12, :, :][::2, :, :]
SCA_ERA5 = pd.Series(np.empty(len(SCOV_ERA5.time)) *
                     np.nan, index=SCOV_ERA5.time.values)

for th in [20, 45, 70, 95]:
    t1 = datetime.datetime.now()
    for i, time in enumerate(SCOV_ERA5.time.values):
        sca = (SCOV_ERA5[i, :, :] > th).values.sum()/49
        SCA_ERA5[time] = sca

    SCA_ERA5 = SCA_ERA5.resample("d").mean()*100
    SCA_ERA5 = SCA_ERA5.reindex(SCA_ianigla.index)

    t2 = datetime.datetime.now()
    dt = t2-t1
    print(dt)
    SCA_ERA5.dropna().to_csv(
        "datos/era5land/RioMaipoEnElManzano/SCA_ERA5LAND_"+str(th)+".csv")

s = [pd.read_csv(i, index_col=0) for i in glob(
    "datos/era5land/RioMaipoEnElManzano/SCA_ERA5LAND*")]
s = pd.concat(s, axis=1)
s.columns = ["20", "45", "70", "95"]
s.to_csv("datos/era5land/RioMaipoEnElManzano/SCA_ERA5LAND_threshold.csv")

# %%


try:
    sca_cortes = pd.read_csv(
        "datos/ANDES_SWE_Cortes/maipomanzano_SCA.csv", index_col=0)
    sca_cortes.index = pd.to_datetime(sca_cortes.index)
except:
    SWE_cortes = [xr.open_dataset(path)
                  for path in glob("datos/ANDES_SWE_Cortes/*.nc")]
    SWE_cortes = [d.assign_coords({"time": np.arange(1, len(d.time)+1)})
                  for d, yr in zip(SWE_cortes, np.arange(1985, 2016))]
    # SWE_cortes    = xr.concat(SWE_cortes, "time")
    # SWE_cortes    = SWE_cortes.assign_coords({"time":np.arange(SWE_cortes.time.shape[0])})
    # SCA_cortes    = pd.Series(np.empty(len(SWE_cortes.time))*np.nan,index=SWE_cortes.time.values)

    models = []
    for swe_limit in [10, 25, 50, 75, 100, 150, 200, 250]:
        SCA_cortes = []
        for j, file in enumerate(SWE_cortes):
            for i in trange(file.time.shape[0]):
                sca = (file.SWE[i, :, :] > swe_limit).sum().load().item()
                SCA_cortes.append(sca)

        SCA_cortes = pd.Series(np.array(SCA_cortes))/467798
        SCA_cortes.index = pd.date_range(
            "1984-04-01", "2015-04-01", freq="1d")[:-1]
        SCA_cortes = SCA_cortes*100
        models.append(SCA_cortes)

    sca_cortes = pd.concat(models, axis=1)
    sca_cortes.columns = [10, 25, 50, 75, 100, 150, 200, 250]
    sca_cortes.to_csv("datos/ANDES_SWE_Cortes/maipomanzano_SCA.csv")
    # for i,time in enumerate(SCA_cortes.index):
    #     sca = (SWE_cortes.SWE[i,:,:]>10).sum().load().item()
    #     SCA_cortes[time] = sca

    # SCA_cortes = SCA_ERA5.resample("d").mean()*100
    # SCA_cortes = SCA_ERA5.reindex(SCA_ianigla.index)

# %%
# SCOV_ERA5    = xr.open_dataset("datos/era5land/maipo/snow_depth_water_equivalent.nc",chunks="auto")
# # SCOV_ERA5    = SCOV_ERA5.sd[23:,:,:][::12,:,:][::2,:,:]
# SCOV_ERA5    = SCOV_ERA5.sd
# SCA_ERA5     = pd.Series(np.empty(len(SCOV_ERA5.time))*np.nan,index=SCOV_ERA5.time.values)

# for i,time in enumerate(SCOV_ERA5.time.values):
#     sca = (SCOV_ERA5[i,:,:]>10*1e-3).values.sum()/49
#     SCA_ERA5[time] = sca

# SCA_ERA5 = SCA_ERA5.resample("d").mean()

# %%


fig = plt.figure(num=0, figsize=(8, 8), dpi=150)
# fig.tight_layout(pad=2)
ax = fig.add_subplot(311)
ax1 = fig.add_subplot(325)
ax2 = fig.add_subplot(326)
ax3 = fig.add_subplot(312)
# ax4 = fig.add_subplot(243)
fig.tight_layout(pad=3)

SCA_cortes = sca_cortes.iloc[:, 0]

SCA_ianigla.plot(ax=ax, alpha=0.7, label="IANIGLA")
SCA_cortes.plot(ax=ax, alpha=0.7, label="Cortes")

ax.legend(loc=(0, 1), ncol=2, frameon=False)
ax.set_xlabel("")
ax.set_ylabel("Fractional Snow\nCovered Area (%)")


mask1 = SCA_ianigla.index.month.map(lambda x: x in [12, 1, 2]).values
mask2 = SCA_ianigla.index.month.map(lambda x: x in [3, 4, 5]).values
mask3 = SCA_ianigla.index.month.map(lambda x: x in [6, 7, 8]).values
mask4 = SCA_ianigla.index.month.map(lambda x: x in [9, 10, 11]).values

ax1.scatter(SCA_ianigla[mask1], SCA_cortes.reindex(SCA_ianigla.index)[
            mask1], alpha=0.3, s=5, color="tab:red", label="summer")
ax1.scatter(SCA_ianigla[mask2], SCA_cortes.reindex(SCA_ianigla.index)[
            mask2], alpha=0.3, s=5, color="tab:brown", label="autumn")
ax1.scatter(SCA_ianigla[mask3], SCA_cortes.reindex(SCA_ianigla.index)[
            mask3], alpha=0.3, s=5, color="tab:cyan", label="winter")
ax1.scatter(SCA_ianigla[mask4], SCA_cortes.reindex(SCA_ianigla.index)[
            mask4], alpha=0.3, s=5, color="tab:green", label="spring")


ax3.boxplot([SCA_ianigla[mask] for mask in [mask1, mask2, mask3, mask4]], patch_artist=True, medianprops={"color": "k"},
            positions=np.arange(4), widths=0.1, sym="", boxprops={"facecolor": "tab:blue"})
ax3.boxplot([SCA_cortes.reindex(SCA_ianigla.index)[mask].dropna() for mask in [mask1, mask2, mask3, mask4]],
            positions=np.arange(4)+0.15, widths=0.1, sym="", patch_artist=True, medianprops={"color": "k"},
            boxprops={"facecolor": "tab:orange"})
ax3.set_xticks(np.arange(4)+0.05)
ax3.set_xticklabels(["summer", "autumn", "winter", "spring"])
ax3.set_ylabel("Fractional Snow\nCovered Area (%)")

# ax3.scatter(SCA_ianigla[mask1],SCA_cortes.reindex(SCA_ianigla.index)[mask1],alpha=0.3,s=5,color="tab:red",label="summer")
# ax3.scatter(SCA_ianigla[mask2],SCA_cortes.reindex(SCA_ianigla.index)[mask2],alpha=0.3,s=5,color="tab:brown",label="autumn")
# ax3.scatter(SCA_ianigla[mask3],SCA_cortes.reindex(SCA_ianigla.index)[mask3],alpha=0.3,s=5,color="tab:cyan",label="winter")
# ax3.scatter(SCA_ianigla[mask4],SCA_cortes.reindex(SCA_ianigla.index)[mask4],alpha=0.3,s=5,color="tab:green",label="spring")


ax1.plot([0, 100], [0, 100], "k--")
ax1.set_ylabel("fSCA_Cortes")
ax1.set_xlabel("fSCA_IANIGLA")
ax1.legend(frameon=False, loc=(0, 1), ncol=2)


data1 = SCA_ianigla.groupby(
    [SCA_ianigla.index.year, SCA_ianigla.index.dayofyear]).mean().unstack().T
data2 = SCA_cortes.reindex(SCA_ianigla.index).groupby(
    [SCA_ianigla.index.year, SCA_ianigla.index.dayofyear]).mean().unstack().T

data1.plot(legend=False, ax=ax2, color="tab:blue", alpha=0.1)
data2.plot(legend=False, ax=ax2, color="tab:orange", alpha=0.1)

data1.mean(axis=1).plot(ax=ax2)
data2.mean(axis=1).plot(ax=ax2)
ax2.set_xlabel("Day of year")
ax2.set_xticks(np.arange(0, 366)[::50])
ax2.set_ylabel("fSCA Annual Cycle (%)")

plt.savefig("plots/maipomanzano/cortes_SCA_study_" +
            str(SCA_cortes.name)+"_.pdf", dpi=150, bbox_inches="tight")

# %%

fig = plt.figure(num=0, figsize=(8, 8), dpi=150)
# fig.tight_layout(pad=2)
ax = fig.add_subplot(311)
ax1 = fig.add_subplot(325)
ax2 = fig.add_subplot(326)
ax3 = fig.add_subplot(312)
# ax4 = fig.add_subplot(243)
fig.tight_layout(pad=3)

SCA_ianigla.plot(ax=ax, alpha=0.7, label="IANIGLA")
SCA_ERA5.plot(ax=ax, alpha=0.7, label="ERA5-Land")

ax.legend(loc=(0, 1), ncol=2, frameon=False)
ax.set_xlabel("")
ax.set_ylabel("Fractional Snow\nCovered Area (%)")

mask1 = SCA_ianigla.index.month.map(lambda x: x in [12, 1, 2]).values
mask2 = SCA_ianigla.index.month.map(lambda x: x in [3, 4, 5]).values
mask3 = SCA_ianigla.index.month.map(lambda x: x in [6, 7, 8]).values
mask4 = SCA_ianigla.index.month.map(lambda x: x in [9, 10, 11]).values

ax1.scatter(SCA_ianigla[mask1], SCA_ERA5[mask1],
            alpha=0.3, s=5, color="tab:red", label="summer")
ax1.scatter(SCA_ianigla[mask2], SCA_ERA5[mask2],
            alpha=0.3, s=5, color="tab:brown", label="autumn")
ax1.scatter(SCA_ianigla[mask3], SCA_ERA5[mask3],
            alpha=0.3, s=5, color="tab:cyan", label="winter")
ax1.scatter(SCA_ianigla[mask4], SCA_ERA5[mask4],
            alpha=0.3, s=5, color="tab:green", label="spring")

ax1.plot([0, 100], [0, 100], "k--")
ax1.set_ylabel("fSCA_ERA5LAND")
ax1.set_xlabel("fSCA_IANIGLA")
ax1.legend(frameon=False)


ax3.boxplot([SCA_ianigla[mask] for mask in [mask1, mask2, mask3, mask4]], patch_artist=True, medianprops={"color": "k"},
            positions=np.arange(4), widths=0.1, sym="", boxprops={"facecolor": "tab:blue"})
ax3.boxplot([SCA_ERA5.reindex(SCA_ianigla.index)[mask].dropna() for mask in [mask1, mask2, mask3, mask4]],
            positions=np.arange(4)+0.15, widths=0.1, sym="", patch_artist=True, medianprops={"color": "k"},
            boxprops={"facecolor": "tab:orange"})
ax3.set_xticks(np.arange(4)+0.05)
ax3.set_xticklabels(["summer", "autumn", "winter", "spring"])
ax3.set_ylabel("Fractional Snow\nCovered Area (%)")

data1 = SCA_ianigla.groupby(
    [SCA_ianigla.index.year, SCA_ianigla.index.dayofyear]).mean().unstack().T
data2 = SCA_ERA5.groupby(
    [SCA_ianigla.index.year, SCA_ianigla.index.dayofyear]).mean().unstack().T

data1.plot(legend=False, ax=ax2, color="tab:blue", alpha=0.1)
data2.plot(legend=False, ax=ax2, color="tab:orange", alpha=0.1)

data1.mean(axis=1).plot(ax=ax2)
data2.mean(axis=1).plot(ax=ax2)
ax2.set_xlabel("Day of year")
ax2.set_xticks(np.arange(0, 366)[::50])
ax2.set_ylabel("fSCA Annual Cycle (%)")

plt.savefig("plots/maipomanzano/ERA5land_SCA_study.pdf",
            dpi=150, bbox_inches="tight")

# tot = pd.concat([SCA_ERA5.reindex(SCA_ianigla.index),SCA_ianigla/100],axis=1)


# cmaplist = pd.read_csv("terraincolormap.txt").values
# cmaplist = [list(np.array(i)/255) for i in cmaplist]
# cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, len(cmaplist))
# fig,ax=plt.subplots()
# basin.boundary.plot(ax=ax,zorder=2,color="k")
# topo.plot(ax=ax,cmap=cmap)
# ax.axis("off")
