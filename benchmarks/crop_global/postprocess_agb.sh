#!/bin/bash 

# get above ground biomass data
function prepare_agb {
	model_input=$1
	output=$2
	dvar=$3
	c=$(head -1 $model_input | awk -v d=$dvar '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	#awk -v c="$c" ' BEGIN { FS = " " } ; FNR==NR{a[$1$2]=$c;next}BEGIN{OFS=" "};{if($1$2 in a){print $3,a[$1$2]*0.7}}' $1 >> $2
	awk -v c="$c" '{print $1,$2,$c}' $model_input > $output
}
# Liu AGB 

dpath="/lunarc/nobackup/users/x_larni/RUNDIR/DDON_crop"
lpath="/home/x_larni/DATA/Liu_AGB_2015"

describe_benchmark "LPJ-GUESS - Global Benchmarks for crops Above Ground Biomass"
source scatter_plot.sh

tslice ${dpath}/cpool.out -f 1993 -t 2012 -o cpool1993-2012.dat
prepare_agb cpool1993-2012.dat cpool1993-2012_agb.dat VegC
joyn cpool1993-2012_agb.dat ${lpath}/Global_mean_ABC_1993-2012_Liu2015_SI.dat -i Lon Lat -fast -o cpool1993-2012_joyned.dat

# delta plot Liu cpool VegC
echo 
awk '{if(FNR==1){print $1,$2, $3} else {print $1,$2, $4}}' cpool1993-2012_joyned.dat > cpool1993-2012_joyned_Liu.dat
awk '{print $1,$2, $3}' cpool1993-2012_joyned.dat > cpool1993-2012_joyned_VegC.dat
delta cpool1993-2012_joyned_VegC.dat cpool1993-2012_joyned_Liu.dat -i Lon Lat -o delta_cpool1993-2012_joyned.dat

# Scatterplot Liu cpool VegC
awk '(FNR>1){print $3, $4}' cpool1993-2012_joyned.dat > scat_cpool.dat
scatter_plot "AGB" "Liu et al. " "LPJ-GUESS" scat_cpool.dat agb.png
describe_image agb.png "Modelled compared to Liu et al. data. Units: kg m-2." embed
 
 

