# ==============================================================================
# Run build ros dataset for each year
# ==============================================================================


#Run for each year
iyr=1984
fyr=2015
for yr in $(seq $iyr $fyr);do
    file="datos/ROS/ROS_Timeseries_${yr}.csv"
    if test -f $file;then
        echo $file already exists...
    else
        sed "s/%YR%/$yr/g" scripts/build_ROS_dataset.py > tmp/tmp.py
        python3 tmp/tmp.py
        rm tmp/tmp.py
    fi
done

#Mix and clean
sed "s/%%TYPE%%/Basinwide/g" scripts/mixandclean_ROSdataset.py > tmp/tmp.py
python3 tmp/tmp.py
sed "s/%%TYPE%%/ROS/g" scripts/mixandclean_ROSdataset.py > tmp/tmp.py
python3 tmp/tmp.py
rm tmp/tmp.py
