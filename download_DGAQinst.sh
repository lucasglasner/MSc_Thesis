#Download data from DGA website

outdir=$1
iyr=1984
fyr=2015

dga=$(cat DGA_cURL)
str1=$(echo $dga | cut -d "&" -f14)
str2=$(echo $dga | cut -d "&" -f15)
str3=$(echo $dga | cut -d "&" -f16)
str4=$(echo $dga | cut -d "&" -f17)

nstr1=$(echo $str1 | sed 's/'$(echo $str1 | cut -d "=" -f2)'//g')
nstr2=$(echo $str2 | sed 's/'$(echo $str2 | cut -d "=" -f2)'//g')
nstr3=$(echo $str3 | sed 's/'$(echo $str3 | cut -d "=" -f2)'//g')
nstr4=$(echo $str4 | sed 's/'$(echo $str4 | cut -d "=" -f2)'//g')


for yr in $(seq $iyr $fyr);do
    echo Downloading Runoff Data from year $yr...
    n1=${nstr1}01%2F01%2F${yr}
    n2=${nstr2}01%2F${yr}
    n3=${nstr3}31%2F12%2F${yr}
    n4=${nstr4}01%2F${yr}

    cURL=$(echo $dga | sed "s/$str1/$n1/g")
    cURL=$(echo $cURL | sed "s/$str2/$n2/g")
    cURL=$(echo $cURL | sed "s/$str3/$n3/g")
    cURL=$(echo $cURL | sed "s/$str4/$n4/g")
    
    echo $cURL --output $yr.xlsx | bash
    mv $yr.xlsx $outdir/$yr.xlsx
done

echo Done


exit


