import cdsapi

c = cdsapi.Client()


c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'mean_sea_level_pressure',
            'total_column_water',
            'total_precipitation',
            'vertical_integral_of_eastward_water_vapour_flux',
            'vertical_integral_of_northward_water_vapour_flux',
        ],
        'year': '2008',
        'month': '05',
        'day': [
            '20', '21', '22',
            '23', '24', '25',
            '26', '27', '28',
            '29', '30', '31',
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
            -20, -100, -50,
            -60,
        ],
    },
    'download_era5_surface.nc')
