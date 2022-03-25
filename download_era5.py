import cdsapi


# p=-32.83,-71.50,"RioAconcanguaEnChacabuquito" #aconcagua
# p=-33.32,-71.72,"RioMapochoEnLosAlmendros" #mapocho
# p=-33.71,-71.68,"RioMaipoEnElManzano" #maipo
# p=-34.33,-72.01,"RioCachapoalEnPteTermasDeCauquenes" #cachapoal
# p=-34.77,-72.15,"RioTinguiriricaBajoLosBriones" #tinguiririca
# p=-35.09,-72.25,"RioTenoDespuesDeJuntaConClaro" #teno
# p=-35.36,-72.50,"RioColoradoEnJuntaConPalos" #colorado
# p=-35.94,-72.78,"RioMauleEnArmerillo" #maule
p=-36.68,-73.01,"RioUblewEnSanFabianN2" #ñuble

# print(p[2])
# print('==> TANDA2 <===')
c = cdsapi.Client()

for yr in range(2011,2021,1):
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': 'zero_degree_level',
            'year': str(yr),
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
                p[0]-1e-3, p[1]-1e-3, p[0]+1e-3,
                p[1]+1e-3,
            ],
            'format': 'netcdf',
        },
        'download_'+p[2]+'_'+str(yr)+'.nc')

