#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:02:18 2021

@author: lucas

# =============================================================================
# Compute basin hypsometry from clipped DEM data.
# =============================================================================

"""

import xarray as xr
import numpy as np
import pandas as pd


basins = ['Rio Aconcagua En Chacabuquito', 'Rio Ancoa En El Morro',
          'Rio Cachapoal En Pte Termas De Cauquenes',
          'Rio Claro En El Valle', 'Rio Claro En Hacienda Las Nieves',
          'Rio Mapocho En Los Almendros', 'Rio Melado En El Salto',
          'Rio Perquilauquen En San Manuel', 'Rio Uble En San Fabian N 2',
          'Rio Tinguiririca Bajo Los Briones', 'Rio Achibueno En La Recova',
          'Rio Maipo En El Manzano', 'Rio Longavi En El Castillo',
          'Rio Teno Despues De Junta Con Claro',
          'Rio Palos En Junta Con Colorado',
          'Rio Colorado En Junta Con Palos', 'Rio Upeo En Upeo',
          'Rio Lircay En Puente Las Rastras', 'Rio Claro En Camarico']

# basins = ['Rio Upeo En Upeo']
# %%
# =============================================================================
# Inputs
# =============================================================================
for name in basins:
    name = name.replace(" ", "")
    print(name)
    path_basinDEM = "datos/topography/basins/"+name+".nc"
    basin = xr.open_dataset(path_basinDEM, chunks="auto").Band1

    # Compute Pixel Area
    earth_radius = 6371
    lon2d, lat2d = np.meshgrid(basin.lon, basin.lat)
    dlon, dlat = np.deg2rad(np.diff(basin.lon)[0]), np.deg2rad(
        np.diff(basin.lat)[0])
    areas = earth_radius**2*np.cos(np.deg2rad(lat2d))*dlon*dlat


# =============================================================================
# Compute hypsometry
# =============================================================================

    heights = np.linspace(0, 7e3, 500)
    basin_area = np.empty(heights.shape)
    for i in range(len(heights)):
        basin_area[i] = areas[basin.values < heights[i]].sum()

    fraction_area = basin_area/basin_area.max()

# =============================================================================
# Save hypsometry data
# =============================================================================

    data = pd.DataFrame([heights, fraction_area, basin_area], index=[
                        "height", "fArea", "Area_km2"]).T
    np.round(data, 3).to_csv(path_basinDEM[:-3].replace(name, "")+"hypso/"+name +
                             "_hypso.csv", index=False)
