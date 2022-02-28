import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'geopotential', 'temperature', 'u_component_of_wind',
            'v_component_of_wind', 'vorticity',
        ],
        'pressure_level': [
            '500', '850', '1000',
        ],
        'year': '2013',
        'month': '08',
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
        'time': '00:00',
        'area': [
            -20, -100, -50,
            -60,
        ],
    },
    'download_era5_Ago2013_upper.nc')
