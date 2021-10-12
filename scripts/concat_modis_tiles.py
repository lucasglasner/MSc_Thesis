#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:16:08 2021

@author: lucas
"""

import xarray as xr
import sys
import pandas as pd
import datetime
import os
import glob
import numpy as np

#%%
# =============================================================================
# Path to files, output name, defitions and collect metadata
# =============================================================================
# Some names
basin        = "RioMaipoEnElManzano"
modis_sensor = "MYD10A1"
modis_path   = "datos/modis/"+modis_sensor+"/"
modis_files  = glob.glob(modis_path+"*["+basin+"].nc")
iy,fy        = 2000,2021
output_path  = "datos/modis/"

# Create metadata and grab coordinates
metadata     = xr.open_dataset(modis_files[0])
crs          = metadata.crs.copy()
coords       = metadata.coords
dims         = metadata.dims
crs.attrs.pop("long_name")
metadata.attrs.pop("history")
attrs        = {**{"author":"Lucas Glasner 2021",
                   "contact":"lgvivanco96@gmail.com",
                   "short_name":"fSCA"},
                **metadata.Band1.attrs,
                **metadata.attrs,**crs.attrs}
attrs["long_name"]  = "Fraction of Snow Cover Area"
del metadata,crs
#%%
# =============================================================================
# Load data, build dates and stack
# =============================================================================
time = []
data = []
for file in modis_files:
    date       = file.replace(modis_path,"").replace(modis_sensor,"")[:9][2:]
    year       = int(date[:4])
    if iy<=year<=fy:
        julian_day = int(date[-3:])-1
        date       = datetime.datetime(year,1,1)+datetime.timedelta(julian_day)
        print(date.strftime("%Y-%m-%d"))
        time.append(date)
        data.append(xr.open_dataset(file,chunks="auto").Band1.data)
data = np.stack(data,axis=0) #Stack in new dimension

#%%
# =============================================================================
# #Sort time and tiles in ascending order.
# =============================================================================
data = data[np.array(pd.Series(time).sort_values().index),:,:]
time = np.array(time)[np.array(pd.Series(time).sort_values().index)]

#%%
# =============================================================================
# Build xarray multidimensional raster and save to output
# =============================================================================

data = xr.DataArray(data=data,
                    dims=["time","lat", "lon"],
                    coords=dict(lon=(["lon"], coords.variables["lon"]),
                                lat=(["lat"], coords.variables["lat"]),
                                time=time,
                                ),
                    attrs=attrs,
                    name="fSCA")

#Fill with NaN days without image
data = data.reindex(time=pd.date_range(str(iy),str(fy+1),freq="1d")[:-1]) 
#%%
print("Saving output...")
for yr in list(map(lambda x: str(x),range(iy,fy+1))):
    print(yr)
    data.sel(time=yr).to_netcdf(output_path+modis_sensor+"_"+yr+".nc")
# print("Mergetime...")
# os.system("cdo mergetime "+
#           output_path+modis_sensor+"**.nc "+
#           output_path+modis_sensor+"_"+str(iy)+"-"+str(fy)+".nc")
# os.system("rm -rf "+output_path+modis_sensor+"*[!"+str(iy)+"-"+str(fy)+"]*.nc")
print("Done")