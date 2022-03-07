#!bin/bash
# =============================================================================
# Use other scripts for clip a dataset with all basins
# =============================================================================

basin_path=$(ls datos/vector/basins/*.shp)
data_folder=datos/ROS/CORTES_CR2MET_ERA5/
data_to_clip=$(ls $data_folder/ROS*)


for vector in $basin_path; do
    basin_name=$(echo $vector | cut -c21-100)
    echo ${basin_name%.shp}
    if [[ -d "${data_folder}/${basin_name}" ]]; then
        :
    else
        mkdir ${data_folder}/${basin_name}
    fi
    for file in $data_to_clip; do
        :
    done
done

exit
