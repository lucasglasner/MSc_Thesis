#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:03:31 2021

@author: lucas
"""

import cdsapi

c = cdsapi.Client()

# variables = ['runoff', 'snow_cover', 'snow_depth_water_equivalent',
#               'surface_runoff', 'total_precipitation','2m_temperature']
variables = ['2m_temperature']
for var in variables:
    for yr in list(map(lambda x: str(x),range(1985,2022))):
        print("Requesting: "+var)
        c.retrieve(
            'reanalysis-era5-land',
            {
                'variable': var,
                'year': yr,
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'day': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    '13', '14', '15',
                    '16', '17', '18',
                    '19', '20', '21',
                    '22', '23', '24',
                    '25', '26', '27',
                    '28', '29', '30',
                    '31',
                ],
                'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
                ],
                'area': [
                    -27, -75, -40,
                    -65,
                ],
                'format': 'netcdf',
            },
            'datos/era5land/hourly/'+var+"_"+yr+'.nc')


