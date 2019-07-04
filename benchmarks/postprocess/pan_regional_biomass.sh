#!/bin/bash -x

# set regions used in fig 1 in Pan et al. 2011, "A large and persistent Carbon sink in the World's forests"
# doi:10.1126/science.1201609
# Russia taken as one region, "other temperate countries" skipped
regnames="Russia Canada N_Europe USA Europe China Japan South_Korea Australia New_Zealand South_Asia Africa Americas"
panyears="2007"

root=$(dirname $0)
for y in $panyears; do 
    case $y in
	"1990")
	    cnt=1
	    ;;
	"2000")
	    cnt=2
	    ;;
	"2007")
	    cnt=3
	    ;;
	*)
	    echo"Wrong year in pan_regional_biomass.sh"
	    exit -1
	    ;;
    esac
    echo " C-pools compared to Pan et al. for $y" > pan_${y}.dat
    echo "  VegC LitterC    DWD  SoilC Total" >> pan_${y}.dat
    # read pan_regional data file and skip header 
    (( headlines = 4+cnt*5 ))
    head -n $headlines Pan/pan_regional_data.txt | tail -n 5 > pan_tmp_$y
    source pan_tmp_$y
    tslice cpool_natural.out -f $y -t $y -o cpool_natural_${y}.dat
    nreg=0
    for reg in $regnames; do 
	# joyn regional "grid_lists" with output files
	joyn Pan/gridlist_${reg}.txt cpool_natural_${y}.dat -i Lon Lat -o cpool_natural_${y}_${reg}.dat
	aslice cpool_natural_${y}_${reg}.dat -n -sum "kg/m2->Pg" -o cpool_natural_${y}_${reg}_tot.dat
	# here now weboutput
	echo "$reg" >> pan_${y}.dat
	awk -v t=${TLB[$nreg]} -v d=${DWD[$nreg]} -v l=${LIT[$nreg]} -v s=${SOI[$nreg]} '{OFS="\t"; if (FNR>1) {printf "%6.2f %6.2f %6.2f %6.2f %6.2f %10s\n%6.2f %6.2f %6.2f %6.2f %6.2f %10s", $1, $2, 0, $3, $1+$2+$3,"LPJ-GUESS",t,l,d,s,t+l+d+s," Pan et al.\n"}} ' cpool_natural_${y}_${reg}_tot.dat >> pan_${y}.dat
	# remove intermediate files
	rm -f cpool_natural_${y}_${reg}.dat cpool_natural_${y}_${reg}_tot.dat
	((nreg=nreg+1))
    done
    describe_textfile pan_${y}.dat "Forest Carbon Pools compared to Pan et al. regional dataset [Pg(C)]"
    # remove intermediate files
    rm -f pan_tmp_$y cpool_natural_${y}.dat 
    ((cnt=cnt+1))
done
