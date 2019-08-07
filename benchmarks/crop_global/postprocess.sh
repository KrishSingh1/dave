#!/bin/bash

# Function for creating a scatter plot using gnuplot.
#
# Parameters:
# $1 model input 
# $2 observation input
# $3 output file
# $4 crop name
function prepareyielddata {
	model_input=$1
	obs_input=$2
	output=$3
	crop=$4
	c=$(head -1 $model_input | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	awk -v c="$c" ' BEGIN { FS = " " } ; FNR==NR{a[$1$2]=$c;next}BEGIN{OFS=" "};{if($1$2 in a){print $3,a[$1$2]*1.3}}' $1 $2 > $3
}

# get above ground biomass data
function prepare_agb {
	model_input=$1
	output=$2
	dvar=$3
	c=$(head -1 $model_input | awk -v d=$dvar '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	awk -v c="$c" '(FNR==1){print $1,$2,$c}' $model_input > $output
	awk -v c="$c" '(FNR>1){print $1,$2,$c*0.7}' $model_input >> $output
}

describe_benchmark "LPJ-GUESS - Global Benchmarks for crops"
source scatter_plot.sh

# link data-dirs for biomass and fire
if [ $ARCH == "aurora" ]
then
    datapath=/lunarc/nobackup/projects/lpjguess/data/
else
    datapath=/data/
fi
ln -sf ${datapath}/biomass/Global_mean_ABC_1993-2012_Liu2015_SI.dat
ln -sf ${datapath}/biomass/Pan
ln -sf ${datapath}/fire/gfed40_c-emissions_1997-2016.dat
ln -sf ${datapath}/fire/gfed_regions0.5.dat
ln -sf ${datapath}/benchmark_data/2015_12_14/landuse

common1961to1990.sh

tslice cflux.out -o cflux1990to2000.txt -f 1990 -t 2000 -lon 1 -lat 2 -y 3
aslice cflux1961to1990.txt -o cflux1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cflux1961to1990_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1961 to 1990. Units: Pg C/y"
aslice cflux1990to2000.txt -o cflux1990to2000_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cflux1990to2000_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1990 to 2000. Units: Pg C/y"

tslice nflux.out -o nflux1990to2000.txt -f 1990 -t 2000 -lon 1 -lat 2 -y 3
aslice nflux1961to1990.txt -o nflux1961to1990_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile nflux1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Fluxes, 1961 to 1990. Units: Tg N/y"
aslice nflux1990to2000.txt -o nflux1990to2000_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile nflux1990to2000_areaaverage.txt "Global Terrestrial Nitrogen Fluxes, 1990 to 2000. Units: Tg N/y"

aslice cpool1961to1990.txt -o cpool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg'  -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cpool1961to1990_areaaverage.txt "Global Terrestrial Carbon Pools, 1961 to 1990. Units: Pg C/y"

aslice npool1961to1990.txt -o npool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile npool1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Pools, 1961 to 1990. Units: Pg N/y"

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0   
describe_textfile tot_runoff1961to1990_areaaverage.txt "Global Runoff, 1961 to 1990. Units: km3 yr-1"

compute cpool.out -n -o cpool_total.out -i Lon Lat Year Total
compute cflux.out -n -o cflux_nee.out -i Lon Lat Year NEE
balance -pool cpool_total.out -flux cflux_nee.out -start 1901 -end 2006 -matter C
describe_textfile Cbalance_totalerror_GtC.txt "Global Terrestrial Carbon Uptake, 1901 to 2006. /
Determined using C pools (pool_GtC), Cumulative C fluxes (flux_GtC), and their absolute difference (absdiff_GtC)"

compute npool.out -n -o npool_total.out -i Lon Lat Year Total
compute nflux.out -n -o nflux_nee.out -i Lon Lat Year 'nee_m2=NEE/10000'
balance -pool npool_total.out -flux nflux_nee.out -start 1901 -end 2006 -matter N
describe_textfile Nbalance_totalerror_GtN.txt "Global Terrestrial Nitrogen Uptake, 1901 to 2006. /
Determined using N pools (pool_GtN), Cumulative N fluxes (flux_GtN), and their absolute difference (absdiff_GtN)"

tslice yield.out -f 1996 -t 2005 -o yield1996to2005.txt
prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_maize.dat temp_maize.dat TeCo
scatter_plot "Maize yields" "SPAM" "LPJ-GUESS" temp_maize.dat maize_yield.png
describe_image maize_yield.png "Modelled compared to SPAM data set. Units: kg m-2." embed

prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_wheat.dat temp_wheat.dat TeWW
scatter_plot "Wheat yields" "SPAM" "LPJ-GUESS" temp_wheat.dat wheat_yield.png
describe_image wheat_yield.png "Modelled compared to SPAM data set. Units: kg m-2." embed

#===============================================================================
# Above Ground Biomass    
# If benchmarks are run on Aurora or Simba link Liu-AGB 
# and gfed benchmarks into crop_global dir
if [ -f Global_mean_ABC_1993-2012_Liu2015_SI.dat ]
then
	tslice cpool.out -f 1993 -t 2012 -o cpool1993-2012.dat
	prepare_agb cpool1993-2012.dat cpool1993-2012_agb.dat VegC
	joyn Global_mean_ABC_1993-2012_Liu2015_SI.dat cpool1993-2012_agb.dat -i Lon Lat -fast -o cpool1993-2012_joyned.dat
    
	# delta plot Liu cpool VegC 
	awk '{if(FNR==1){print $1,$2, $4} else {print $1,$2, $3}}' cpool1993-2012_joyned.dat > cpool1993-2012_joyned_Liu.dat
	awk '{print $1,$2, $4}' cpool1993-2012_joyned.dat > cpool1993-2012_joyned_VegC.dat
	delta cpool1993-2012_joyned_VegC.dat cpool1993-2012_joyned_Liu.dat -i Lon Lat -o delta_cpool1993-2012_joyned.dat
	gmap delta_cpool1993-2012_joyned.dat -i VegC -lon 1 -lat 2 -landscape -s -20 2 20  -o delta_cpool1993-2012_joyned.png -t "VegC LPJ-GUESS - Liu kg(C)/m2" -c BLUE RED
	convert -geometry 25%x25% delta_cpool1993-2012_joyned.png tmp.png
	convert -rotate 90 tmp.png delta_cpool1993-2012_joyned.png
	describe_image  delta_cpool1993-2012_joyned.png "Modelled minus Liu et al. data. Units: kg m-2." embed
	
	. postprocess_AGB.sh
	# delta plot Liu cpool VegC against Jackson 
	joyn lu_cmass_agb_1993-2012_tot.dat cpool1993-2012_joyned.dat -i Lon Lat -o lu_cmass_agb_tot_1993-2012_joyned.dat
	awk '{if(FNR==1){print $1,$2, "VegC"} else {print $1,$2, $(NF-1)}}' lu_cmass_agb_tot_1993-2012_joyned.dat > lu_cmass_agb_1993-2012_tot.dat_Liu.dat
	awk '{print $1,$2, $NF}' lu_cmass_agb_tot_1993-2012_joyned.dat > cpool1993-2012_joyned_VegC.dat
	delta cpool1993-2012_joyned_VegC.dat lu_cmass_agb_1993-2012_tot.dat_Liu.dat -i Lon Lat -o delta_cpool1993-2012_joyned_jackson.dat
	gmap delta_cpool1993-2012_joyned_jackson.dat -i VegC -lon 1 -lat 2 -landscape -s -20 2 20  -o delta_cpool1993-2012_joyned_jackson.png -t "VegC LPJ-GUESS - Liu kg(C)/m2" -c BLUE RED
	convert -geometry 25%x25% delta_cpool1993-2012_joyned_jackson.png tmp.png
	convert -rotate 90 tmp.png delta_cpool1993-2012_joyned_jackson.png
	describe_image  delta_cpool1993-2012_joyned_jackson.png "Jackson AGB Modelled minus Liu et al. data. Units: kg m-2." embed
    
	# Scatterplot Liu cpool VegC (Jackson AGB)
	awk '(FNR>1){print $(NF-1),$(NF-2)}' lu_cmass_agb_tot_1993-2012_joyned.dat > scat_cpool2.dat
	scatter_plot "AGB" "Liu et al. " "LPJ-GUESS following Jackson" scat_cpool2.dat agb.png
	describe_image agb.png "Jackson AGB Modelled compared to Liu et al. data. Units: kg m-2." embed

	# Scatterplot Liu cpool VegC
	awk '(FNR>1){print $3, $4}' cpool1993-2012_joyned.dat > scat_cpool.dat
	scatter_plot "AGB" "Liu et al. " "LPJ-GUESS" scat_cpool.dat agb_0.7.png
	describe_image agb_0.7.png "Modelled compared to Liu et al. data. Units: kg m-2." embed
  
	# remove intermediate files
	rm -f  cpool1993-2012.dat cpool1993-2012_joyned.dat cpool1993-2012_joyned_Liu.dat delta_cpool1993-2012_joyned.dat cpool1993-2012_joyned_VegC.dat scat_cpool.dat
else
	# data files only available on Simba and Aurora
	echo "Dataset 'Liu Above-Ground-Biomass' not found. Skipping..."
fi

# pan biomass 

. pan_regional_biomass.sh

#===============================================================================
#Fire related benchmarks
if [ -f gfed40_c-emissions_1997-2016.dat ]; then
	tslice cflux.out -f 1997 -t 2016 -o cflux1997-2016.dat
	joyn cflux1997-2016.dat gfed40_c-emissions_1997-2016.dat -i Lon Lat -fast -o cflux1997-2016_joyned.dat

	# Plot fire emissions 
	gmap cflux1997-2016_joyned.dat -i Fire -lon 1 -lat 2 -landscape -o cflux1997-2016_blaze.png \
	-legend common/legend_fire_emis.txt -t "BLAZE mean annual C-emissions [kg(C)/m2a]"
	convert -geometry 25%x25% cflux1997-2016_blaze.png tmp.png
	convert -rotate 90 tmp.png cflux1997-2016_blaze.png
	describe_image  cflux1997-2016_blaze.png "BLAZE Mean annual C-emissions 1997-2016 [kg(C)/m2a]" embed
	
	# Plot gfed 4.0 emissions
	gmap cflux1997-2016_joyned.dat -i C-Emis -lon 1 -lat 2 -landscape -o cflux1997-2016_gfed4.png \
	-legend common/legend_fire_emis.txt -t "GFED 4.0 mean annual C-emissions [kg(C)/m2a]"
	convert -geometry 25%x25% cflux1997-2016_gfed4.png tmp.png
	convert -rotate 90 tmp.png cflux1997-2016_gfed4.png
	describe_image  cflux1997-2016_gfed4.png "GFED 4.0 C-emissions kg(C)/m2a." embed
	
	# delta plot gfed4 cflux
	awk '{print $1,$2, $6}' cflux1997-2016_joyned.dat > cflux1997-2016_joyned_Fire.dat
	awk '{if(FNR==1){print $1,$2, $6} else {print $1,$2, $13}}' cflux1997-2016_joyned.dat > cflux1997-2016_joyned_gfed.dat
	delta  cflux1997-2016_joyned_Fire.dat cflux1997-2016_joyned_gfed.dat -i Lon Lat -o delta_cflux1997-2016_joyned.dat
	gmap delta_cflux1997-2016_joyned.dat -i Fire -lon 1 -lat 2 -landscape \
	-legend common/legend_delta_fire_emis.txt -o delta_cflux1997-2016_joyned.png \
	-t "Fire C flux LPJ-GUESS - Gfed kg(C)/m2/a" -c BLUE RED -vert
	convert -geometry 25%x25% delta_cflux1997-2016_joyned.png tmp.png
	convert -rotate 90 tmp.png delta_cflux1997-2016_joyned.png
	describe_image  delta_cflux1997-2016_joyned.png "Modelled minus GFED 4.0 data. Units: kg(C)/m2a." embed
	
	# Scatterplot GFED C-emis 
	awk '(FNR>1){print $6, $13}' cflux1997-2016_joyned.dat > scat_fire_cflux.dat
	scatter_plot "Fire C-Flux" "GFED4.0" "LPJ-GUESS" scat_fire_cflux.dat scat_fire_cflux.png
	describe_image scat_fire_cflux.png "Modelled compared to GFED 4.0 C-Emissions Units: kg(C)/m2a" embed
	
	# A-slicing over regions 0.5 deg res
	GFEDreg=(BONA TENA CEAM NHSA SHSA EURO MIDE NHAF SHAF BOAS TEAS CEAS EQAS AUST)

	tot_lpjg=0.
	tot_gfed=0.
	for ((x=1; x<=14; x++)); do
		((xx=$x-1))
		creg=${GFEDreg[${xx}]} 
		awk -v reg=$x '(FNR==1 || $3==reg){print $0}' ~/DATA/gfed_regions0.5.dat > reg.dat
		joyn cflux1997-2016_joyned.dat reg.dat -i Lon Lat -fast -o cflux_reg_${x}_joyned.dat  
		aslice cflux_reg_${x}_joyned.dat -n -lon Lon -lat Lat  -sum "kg/m2->Pg" -o tot_cflux_reg_${x}.dat
		if [ $x -eq 1 ]; then
			echo "Region LPJ-GUESS GFED 4.0 "	> tot_cflux_reg.dat
		fi
		awk -v reg=$creg '(FNR==2){printf "%s	  %6.2f   %6.2f \n",reg,$4*1000,$11*1000}' \
		tot_cflux_reg_${x}.dat >> tot_cflux_reg.dat
		# remove intermediate files
		rm -f tot_cflux_reg_${x}.dat cflux_reg_${x}_joyned.dat reg.dat
	done
	tot_lpjg=$(awk '(FNR>1){sum+=$2} END {print sum}' tot_cflux_reg.dat)
	tot_gfed=$(awk '(FNR>1){sum+=$3} END {print sum}' tot_cflux_reg.dat)
	printf "Total	%6.2f  %6.2f\n" $tot_lpjg $tot_gfed >> tot_cflux_reg.dat
	describe_textfile tot_cflux_reg.dat "Fire C-emissions per GFED - region [Pg/a]"
else
	# data files only available on Simba and Aurora
	echo "Dataset for GFED fire benchmark not found. Skipping..."
fi
