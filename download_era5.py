import cdsapi

c = cdsapi.Client()


c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'geopotential', 'temperature', 'u_component_of_wind',
            'v_component_of_wind',
        ],
        'pressure_level': [
            '300', '500', '700',
            '800', '850', '900',
            '1000',
        ],
        'year': '2013',
        'month': '08',
        'day': [
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14',
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
    'download_era5_upper.nc')
