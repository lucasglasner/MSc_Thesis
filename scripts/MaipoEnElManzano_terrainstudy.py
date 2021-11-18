#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 16:51:53 2021

@author: lucas


# =============================================================================
# Hypsometry study on Maipo en el Manzano River Basin
# =============================================================================

"""

import matplotlib as mpl
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs

import cartopy.feature as cfeature
#%%
# =============================================================================
# Open files
# =============================================================================

chile = gpd.read_file("datos/vector/cl_continental_geo.shp")

cuenca = "RioMaipoEnElManzano"
#Basin shapefile
basin = gpd.read_file("datos/vector/"+cuenca+".shp")

#Basin DEM
dem_full   = xr.open_dataset("datos/topography/basins/"+cuenca+".nc",chunks="auto").Band1
dem        = xr.open_dataset("datos/topography/basins/"+cuenca+"_regridmodis.nc",chunks="auto").Band1
dem_era5   = xr.open_dataset("datos/topography/basins/"+cuenca+"_ERA5LAND.nc",chunks="auto").z
dem_cr2met = xr.open_dataset("datos/topography/basins/"+cuenca+"_CR2MET.nc",chunks="auto").Band1
dem_cortes = xr.open_dataset("datos/topography/basins/"+cuenca+"_Cortes.nc",chunks="auto").Band1

#Basin River Network
rivers = gpd.read_file("datos/vector/red_hidrografica.shp")
rivers = gpd.clip(rivers,basin.to_crs(rivers.crs))
rivers = rivers.to_crs(basin.crs)

#LAND USE
landcover = xr.open_dataset("datos/landcover/basins/"+cuenca+".nc",
                            chunks="auto").Band1

#Aspect
aspect    = xr.open_dataset("datos/topography/basins/aspect/"+cuenca+".nc",chunks="auto").Band1

#slope

slope     = xr.open_dataset("datos/topography/basins/slope/"+cuenca+"_regridmodis.nc",chunks="auto").Band1  
slope     = slope.where(slope<1e7).where(slope>0)/1e5




#%%
# =============================================================================
# Define land cover classes
# =============================================================================
#Land use values
land_uses = {"Cropland":(100,200),
             "Native Forest":(200,229),
             "Forest Plantation":(229,300),
             "Grassland":(300,400),
             "Shrubland":(400,500),
             "Wetland":(500,600),
             "Water Bodies":(600,800),
             "Waterproofs":(800,900),
             "Barren":(900,1000),
             "Ice/Snow":(1000,1200),
             "Clouds":(1200,1500)}  
            

#%%
# =============================================================================
# Hypsometric curve and Hypsometric density
# =============================================================================

#Elevation band width in meters
dz   = 150
tpix = (~np.isnan(dem.values)).sum()
tpix_full   = (~np.isnan(dem_full.values)).sum()
tpix_era5   = (~np.isnan(dem_era5.values)).sum()
tpix_cr2met = (~np.isnan(dem_cr2met.values)).sum()
tpix_cortes = (~np.isnan(dem_cortes.values)).sum()
#elevation bands
elevation_bands = np.arange(dem.min()//10*10-2*dz,
                            dem.max()//10*10+2*dz,dz,
                            dtype=np.int)

#bands number of pixels
bands_size      = np.empty(len(elevation_bands))*np.nan
#band mean slope
bands_slope     = np.empty(len(elevation_bands))*np.nan
#hypsometric curve
hypso_curve      = np.empty(len(elevation_bands))*np.nan
hypso_curve1     = np.empty(len(elevation_bands))*np.nan
hypso_curve2     = np.empty(len(elevation_bands))*np.nan
hypso_curve3     = np.empty(len(elevation_bands))*np.nan
hypso_curve4     = np.empty(len(elevation_bands))*np.nan
#hypsometric density (band number of pixels / total pixels)
hypso_density0   = np.empty(len(elevation_bands))*np.nan
hypso_density1   = np.empty(len(elevation_bands))*np.nan
hypso_density2   = np.empty(len(elevation_bands))*np.nan
hypso_density3   = np.empty(len(elevation_bands))*np.nan
hypso_density4   = np.empty(len(elevation_bands))*np.nan
for i,height in enumerate(elevation_bands):
    mask1  = (~np.isnan(dem.where(dem>height).where(dem<height+dz)))
    mask = (~np.isnan(dem_full.where(dem_full>height).where(dem_full<height+dz)))
    mask4 = (~np.isnan(dem_era5.where(dem_era5>height).where(dem_era5<height+dz)))
    mask3 = (~np.isnan(dem_cr2met.where(dem_cr2met>height).where(dem_cr2met<height+dz)))
    mask2 = (~np.isnan(dem_cortes.where(dem_cortes>height).where(dem_cortes<height+dz)))
    bands_size[i]    = mask.sum()
    hypso_density0[i] = bands_size[i]/tpix_full
    hypso_density1[i]= mask1.sum()/tpix
    hypso_density2[i]= mask2.sum()/tpix_cortes
    hypso_density3[i]= mask3.sum()/tpix_cr2met
    hypso_density4[i]= mask4.sum()/tpix_era5
    bands_slope[i]   = np.nanmean(slope.where(mask).values)
hypso_curve0 = np.cumsum(hypso_density0)
hypso_curve1 = np.cumsum(hypso_density1)
hypso_curve2 = np.cumsum(hypso_density2)
hypso_curve3 = np.cumsum(hypso_density3)
hypso_curve4 = np.cumsum(hypso_density4)

hypso_curve = hypso_curve0
hypso_density = hypso_density0
#%% compare terrain models

fig,ax = plt.subplots(2,4,sharex=True,sharey="row",figsize=(10,4))

titles = ["MODIS 500m","Cortes 1km","CR2MET 5km",  "ERA5-Land 12km"]
colors = plt.cm.tab10(np.linspace(0,1,10))
for i in range(4):
    data1 = eval("hypso_curve"+str(i+1))
    data2 = eval("hypso_density"+str(i+1))
    ax[0,i].plot(elevation_bands,data1)
    ax[0,i].plot(elevation_bands,hypso_curve0,"r--",label="SRTM 30m")
    ax[1,i].plot(elevation_bands,data2)
    ax[1,i].plot(elevation_bands,hypso_density0,"r--")
    ax[0,i].set_title(titles[i])
    ax[1,i].set_xticks(np.arange(0,7e3,1e3))
    ax[1,i].tick_params(axis="x",rotation=45)
    ax[1,i].set_xlim(dem.min(),dem.max())
ax[0,0].legend(frameon=False)
ax[0,0].set_ylabel("Basin Hypsometry")
ax[1,0].set_ylabel("Basin Height Density")
fig.text(0.5,-0.05,"Elevation (m)",ha="center",va="center")

plt.savefig("plots/maipomanzano/datasetcomparison/terrainmodels.pdf",dpi=150,bbox_inches="tight")
#%%
# =============================================================================
# compute percentages of land cover clases and pixel orientations
# =============================================================================

landcover_percent   = np.empty(len(land_uses))
orientation_percent = np.empty(4) 

for i,landclass in enumerate(land_uses.keys()):
    ic,fc = land_uses[landclass][0],land_uses[landclass][1]
    landcover_percent[i]=np.where((landcover>ic) & (landcover<fc),1,0).sum()/6016322
    
degrees = {"N":(-45,45),
           "E":(45,135),
           "S":(135,225),
           "W":(225,315)}
for i,direct in enumerate(degrees.keys()):
    ic,fc = degrees[direct][0],degrees[direct][1]
    orientation_percent[i] = np.where((aspect-45>ic) & (aspect-45<fc),1,0).sum()/463940


#%%

flecha = mpl.image.imread("flechanorte.png")

fig = plt.figure(figsize=(6,5))
axf = fig.add_axes([0.72,1.6,0.2,0.2])
axf.imshow(flecha)
axf.axis("off")
plt.rcParams.update({"font.size":15})
axn  = fig.add_axes([0.6,0.0,1.0,1.0],projection=ccrs.Mercator())
ax0  = fig.add_axes([0.0,0.0,1.0,1.0],projection=ccrs.Mercator())
ax1  = fig.add_axes([0.0,1.0,1.0,1.0],projection=ccrs.Mercator())
ax2  = fig.add_axes([0.6,1.0,1.0,1.0],projection=ccrs.Mercator())

for ax,color in zip([ax0,ax1,ax2,axn],["green","red","blue","purple"]):
    # ax=eval("ax"+str(i))
    ax.axis("off")
    basin.boundary.plot(ax=ax,transform=ccrs.PlateCarree(),color=color,lw=2)
    
ax3  = fig.add_axes([ax2.get_position().xmax+0.3,1.45,0.5,0.5])
ax4  = fig.add_axes([ax2.get_position().xmax+0.3,0.7,0.5,0.5],sharex=ax3)
ax8  = fig.add_axes([ax2.get_position().xmax+0.3,-0.05,0.5,0.5],sharex=ax3)

for ax in [ax3,ax4,ax8]:
    ax.grid(True,ls=":")

ax3.plot(elevation_bands,hypso_density*100)
ax3.set_xticks(elevation_bands[::25])
ax3.set_xticklabels(elevation_bands[::25],rotation=45)
ax3.set_title("Basin area\nheight distribution",loc="left")
ax3.set_ylabel("%Total Area")

ax4.plot(elevation_bands,hypso_curve*100)
ax4.set_xticks(elevation_bands[::25])
ax4.set_xticklabels(elevation_bands[::25],rotation=45)
ax4.set_title("Hypsometric\nCurve",loc="left")
ax4.set_ylabel("%Area below Height")

ax8.plot(elevation_bands,pd.Series(bands_slope).rolling(5).mean())
ax8.set_title("Basin slope\nheight distribution",loc="left")
ax8.set_xticks(elevation_bands[::25])
ax8.set_xticklabels(elevation_bands[::25],rotation=45)
ax8.set_ylabel("Mean slope (%)")
ax8.set_xlabel("Elevation (m)")
# ax8.axvline(5e3,ls="--",color="k")

ax5  = fig.add_axes([ax0.get_position().xmax-0.9,1.3,0.5,0.5])
ax6  = fig.add_axes([ax0.get_position().xmax-0.9,0.2,0.5,0.5])

colors = plt.cm.get_cmap("nipy_spectral_r",4)(np.linspace(0,1,4))
for i in range(len(orientation_percent)):
    ax5.bar(i,orientation_percent[i]*100,color=colors[i])
    
ax5.set_ylim(0,80)
ax5.set_xticks(range(len(orientation_percent)))
ax5.set_xticklabels(list(degrees.keys()))
ax5.set_ylabel("%")
ax5.set_title("Basin inner\norientation",loc="left")

colors = plt.cm.get_cmap("tab10",11)(np.linspace(0,1,11))
for i in range(len(landcover_percent)):
    ax6.bar(i,landcover_percent[i]*100,color=colors[i])
ax6.set_xticks(range(len(landcover_percent)))
ax6.set_xticklabels(list(land_uses.keys()),rotation=90)
ax6.set_ylabel("%")
ax6.set_title("Basin Land Cover",loc="left")
ax6.set_ylim(0,80)


ax7 = fig.add_axes([0.55,0.7,0.6,0.6])
# ax7.axis("off")
ax7.spines["left"].set_visible(False)
ax7.spines["right"].set_visible(False)
ax7.spines["bottom"].set_visible(False)
ax7.spines["top"].set_visible(False)
ax7.grid(True,ls=":")
ax7.set_xticks([])
ax7.set_xticklabels("")
ax7.set_yticks([-20,-30,-40,-50])
ax7.set_yticklabels(["20°S","30°S","40°S","50°S"],fontsize=14)
# ax7.tick_params(color="green",labelcolor="green")
chile.plot(ax=ax7,lw=0.3,color="grey")
basin.plot(ax=ax7,color="red",lw=2)

for ax, color in zip([ax3, ax4, ax5, ax6,ax8], ['blue', 'blue', 'red', 'green','purple']):
    plt.setp(ax.spines.values(), color=color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)
#%
# ax7.set_ylim(-35,-30)
                   #projection=ccrs.Orthographic(-70,-25))
# # ax7.set_extent([-80,-60,-20,-40])
# ax7.set_global()
# # ax7.coastlines()
# ax7.gridlines(ls="--")
# ax7.scatter(-70,-33,transform=ccrs.PlateCarree(),color="red",marker="s")
# ax7.add_feature(cfeature.LAND,facecolor="grey")

# ax7 = fig.add_axes([0.75,
#                     1.2,0.3,0.3],projection=ccrs.Orthographic(-70,-25))
# # ax7.set_extent([-80,-60,-20,-40])
# ax7.set_global()
# # ax7.coastlines()
# ax7.gridlines(ls="--")
# ax7.scatter(-70,-33,transform=ccrs.PlateCarree(),color="red",marker="s")
# ax7.add_feature(cfeature.LAND,facecolor="grey")

# =============================================================================
# Plot basin orography and river network
# =============================================================================

cmaplist = pd.read_csv("terraincolormap.txt").values
cmaplist = [list(np.array(i)/255) for i in cmaplist]
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, len(cmaplist))
topo = ax2.pcolormesh(dem.lon,dem.lat,dem,transform=ccrs.PlateCarree(),
                      cmap=cmap,rasterized=True)
rivers.plot(ax=ax2,transform=ccrs.PlateCarree(),color="darkblue",lw=0.5)

cax = fig.add_axes([ax2.get_position().xmin+0.02,
                    ax2.get_position().ymax+0.04,
                    ax2.get_position().xmax-ax2.get_position().xmin,
                    0.02])
fig.colorbar(topo,cax=cax,orientation="horizontal")
cax.set_title("Elevation (m)",loc="left")

# # # =============================================================================
# # # plot aspect
# # # =============================================================================

# cmap = plt.cm.jet  # define the colormap
# # extract all colors from the .jet map
# cmaplist = [cmap(i) for i in range(cmap.N)]
# # force the first color entry to be grey
# # cmaplist[0] = (.5, .5, .5, 1.0)

# # create the new map
# cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
# # define the bins and normalize
# bounds = [0,90,180,270,360]
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

azimuth = ax1.pcolormesh(aspect.lon,aspect.lat,aspect-45,transform=ccrs.PlateCarree(),
                      cmap=plt.cm.get_cmap("nipy_spectral_r",4),rasterized=True)

cax = fig.add_axes([ax1.get_position().xmin+0.02,
                    ax1.get_position().ymax+0.07,
                    ax1.get_position().xmax-ax1.get_position().xmin,
                    0.02])
cb=fig.colorbar(azimuth,cax=cax,orientation="horizontal",ticks=[0,90,180,270])
cax.set_title("Orientation",loc="left")
cb.set_ticklabels(["N","E","S","W"])

# =============================================================================
# plot land cover
# =============================================================================
cmap = plt.cm.tab10  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
# cmaplist[0] = (.5, .5, .5, 1.0)

# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
bounds = [land_uses[name][0] for name in land_uses.keys()]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

land = ax0.pcolormesh(landcover.lon,landcover.lat,landcover,
                      transform=ccrs.PlateCarree(),cmap=cmap,
                      norm=norm,rasterized=True)

cax = fig.add_axes([ax0.get_position().xmin,
                    ax0.get_position().ymin-0.04,
                    ax0.get_position().xmax-ax0.get_position().xmin,
                    0.02])

ticks = [np.mean([land_uses[name][1],land_uses[name][0]]) for name in land_uses.keys()]
cb = fig.colorbar(land,cax=cax,orientation="horizontal",ticks=ticks)
cb.set_ticklabels(list(land_uses.keys()))
cb.ax.tick_params(rotation=90)

cax.set_title("Land Cover Class",loc="left")

# =============================================================================
# plot slope
# =============================================================================

slopes = axn.pcolormesh(slope.lon,slope.lat,slope,transform=ccrs.PlateCarree(),cmap="bone",rasterized=True)

cax = fig.add_axes([axn.get_position().xmin,
                    axn.get_position().ymin-0.04,
                    axn.get_position().xmax-axn.get_position().xmin,
                    0.02])


cb = fig.colorbar(slopes,cax=cax,orientation="horizontal",ticks=np.arange(0,120,20))
# cb.set_ticklabels(list(land_uses.keys()))
# cb.ax.tick_params(rotation=90)

cax.set_title("Slope (%)",loc="left")


plt.savefig("plots/maipomanzano/terrainstudy.pdf",dpi=150,bbox_inches="tight")