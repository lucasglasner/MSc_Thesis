#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 11:54:35 2021

@author: lucas

# =============================================================================
# Compute ROS events for different methods.
# =============================================================================
"""

import xarray as xr
import numpy as np
import os


compute = True
# %%
if compute:
    # method = 'CORTES - CR2MET'
    # method = "ERA5LAND"
    method = 'CORTES - CR2MET - ERA5'
    yr = "%YR%"
    # yr = "2013"
    print("Current year: "+yr)
    if method == 'ERA5LAND':
        outdir = 'datos/ROS/ERA5LAND/'
        swe_path = 'datos/era5land/snow_depth_water_equivalent_'+yr+'.nc'
        lpr_path = 'datos/era5land/liquid_precipitation_'+yr+'.nc'

        SWE = xr.open_dataset(swe_path, chunks="auto").sd
        LPR = xr.open_dataset(lpr_path, chunks="auto").tp

        print("Making boolean masks...")
        ROS = xr.where((SWE > 10/1e3) & (LPR > 3/1e3), 1.0, 0.0)
        ROS = ROS.reindex({'lon': SWE.lon, 'lat': SWE.lat}, method='nearest')
        ROS = xr.where(SWE > 10/1e3, ROS, -9999)
        ROS = ROS.compute()
        ROS.attrs = {"short_name": "ROS",
                     "long_name": "Rain-over-snow events",
                     "author": "lucasG",
                     "parent_dataset": "ERA5-Land"}
        print("Saving output...")
        ROS = ROS.to_dataset(name='ROS')
        ROS.to_netcdf(outdir+'Rain-Over-Snow_'+yr+'.nc')
        print('Setting missing values...')
        os.system('cdo -P 8 -z zip_9 setmissval,-9999 '+outdir +
                  'Rain-over-Snow_'+yr+'.nc '+outdir+'ROS_'+yr+'.nc')
        os.system('rm -rf '+outdir+'Rain-over-Snow_'+yr+'.nc')
        print("Done")
    if method == "CORTES - CR2MET":
        outdir = 'datos/ROS/CORTES_CR2MET/'
        swe_path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+yr+'.nc'
        pr_path = 'datos/cr2met/CR2MET_pr_'+yr+'.nc'
        t2m_path = 'datos/cr2met/CR2MET_t2m_'+yr+'.nc'

        SWE = xr.open_dataset(swe_path, chunks='auto').SWE
        PR = xr.open_dataset(pr_path, chunks='auto').pr.reindex(
            {"time": SWE.time})
        T2M = xr.open_dataset(t2m_path, chunks='auto').t2m.reindex(
            {"time": SWE.time})

        ROS = xr.where((SWE > 10) & (PR > 3) & (T2M > 0), 1, 0)
        ROS = ROS.reindex({'lon': SWE.lon, 'lat': SWE.lat}, method='nearest')
        ROS = xr.where(SWE > 10, ROS, -9999)
        ROS = ROS.compute()
        ROS.attrs = {"short_name": "ROS",
                     "long_name": "Rain-over-snow events",
                     "author": "lucasG",
                     "parent_dataset": "SWE: Cortes & Margulis et al (2017);\
PR & T2M from CR2MET Boisier et al"}
        print("Saving output...")
        ROS = ROS.to_dataset(name='ROS')
        ROS.to_netcdf(outdir+'Rain-Over-Snow_'+yr+'.nc')
        print('Setting missing values...')
        os.system('cdo -P 8 -z zip_9 setmissval,-9999 '+outdir +
                  'Rain-over-Snow_'+yr+'.nc '+outdir+'ROS_'+yr+'.nc')
        os.system('rm -rf '+outdir+'Rain-over-Snow_'+yr+'.nc')
        print("Done")
    if method == 'CORTES - CR2MET - ERA5':
        outdir = 'datos/ROS/CORTES_CR2MET_ERA5/'
        swe_path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_SWE_'+yr+'.nc'
        pr_path = 'datos/cr2met/CR2MET_pr_'+yr+'.nc'
        H0_path = 'datos/era5/H0_ERA5_'+yr+'.nc'
        dswe_path = 'datos/ANDES_SWE_Cortes/regrid_cr2met/ANDES_dSWE_'+yr+'.nc'

        SWE = xr.open_dataset(swe_path, chunks='auto').SWE
        # dSWE = xr.open_dataset(dswe_path, chunks='auto').SWE.reindex(
        #     {'time': SWE.time}, method='nearest')
        dSWE = (SWE.shift({'time':-1})-SWE.shift({'time':1}))/2
        PR = xr.open_dataset(pr_path, chunks='auto').pr.reindex(
            {'time': SWE.time}, method='nearest')
        FL = xr.open_dataset(H0_path, chunks='auto').deg0l.reindex(
            {'time': SWE.time}, method='nearest')

        ROS = xr.where((SWE > 10) & (PR > 10) & (FL > 300) &
                       (dSWE/SWE.where(SWE>100)<0.05),
                       1, 0)
        ROS = ROS.reindex({'lon': SWE.lon, 'lat': SWE.lat}, method='nearest')
        ROS = xr.where(SWE > 10, ROS, -9999)
        ROS = ROS.compute()
        ROS.attrs = {"short_name": "ROS",
                     "long_name": "Rain-over-snow events",
                     "author": "lucasG",
                     "parent_dataset": "SWE: Cortes & Margulis et al (2017);\
PR: CR2MET, Boisier et al; FL: Zero Degree Level ERA5-ECMWF Reanalysis."}
        print('Saving output...')
        ROS = ROS.to_dataset(name='ROS')
        ROS.to_netcdf(outdir+'/Rain-over-Snow_'+yr+'.nc')
        print('Setting missing values...')
        os.system('cdo -P 8 -z zip_9 setmissval,-9999 '+outdir +
                  'Rain-over-Snow_'+yr+'.nc '+outdir+'ROS_'+yr+'.nc')
        os.system('rm -rf '+outdir+'Rain-over-Snow_'+yr+'.nc')
        print('Done')
# %%
# os.system('for yr in {1984..2015}; do cat scripts/ROS_ChileCentral_computenetcdf.py | sed "s/%YR%/${yr}/g" > tmp/tmpscript; python3 tmp/tmpscript; done')
