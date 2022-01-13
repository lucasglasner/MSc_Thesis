# !usr/bin/bash
#==============================================================================
# This script runs the other python scripts so that a basin static data is
# created. This means that a basin DEM, Land Cover, shapefile polygon and 
# hypsometric curve is created.
#==============================================================================

basins="RioChoapaEnSalamanca RioMaipoEnElManzano"

#Path to Topography Dataset
topodir=../datos/topography/
DEM=${topodir}/Chile_DEM_0001x0001grad.nc

#Path to LandCover Dataset
lcdir=../datos/landcover
LC=${lcdir}/LC_CHILE_2014_b_final.nc

#Path to vector data folder and basin polygons
vector=../datos/vector
polygons=${vector}/cuencas_CAMELS.gpkg


for basin in basins; do
    if test -f "${topodir}/basins/${basin}.nc"; then
        echo ${basin} DEM already exists.
    else
        echo Creating ${basin} DEM ...
        python3 clip_dembasins.py $DEM $polygons $basin
        echo Done
    fi
done

exit
