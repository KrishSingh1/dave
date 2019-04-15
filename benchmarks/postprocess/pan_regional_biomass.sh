#!/bin/bash -x

# set regions used in fig 1 in Pan et al. 2011, "A large and persistent Carbon sink in the World's forests"
# doi:10.1126/science.1201609
# Russia taken as one region, "other temperate countries" skipped
regnames="Russia Canada N_Europe USA Europe China Japan South_Korea Australia New_Zealand South_Asia Africa Americas"
regnames="Russia Canada"
panyears="1990 2000 2007"

trunkpath="/scratch/johan/Benchmarks/trunk_6296-/trunk_7068/output_trunk7068all_crgpp/crop_global/"



root=$(dirname $0)
cnt=1
for y in $panyears
  do 
  echo " C-pools compared to Pan et al. for $y" > pan_${y}.dat
  # read pan_regional data file and skip header 
  (( headlines = 4+cnt*5 ))
  head -n $headlines $root/../postprocess/pan_regional_data.txt | tail -n 5 > pan_tmp_$y
  source pan_tmp_$y
  tslice ${trunkpath}/cpool_natural.out -f $y -t $y -o cpool_natural_${y}.dat
  nreg=0
  for reg in $regnames
    do 
    # joyn regional "grid_lists" with output files
    joyn Pan/gridlist_${reg}.txt cpool_natural_${y}.dat -i Lon Lat -o cpool_natural_${y}_${reg}.dat
    aslice cpool_natural_${y}_${reg}.dat -n -sum "kg/m2->Pg" -o cpool_natural_${y}_${reg}_tot.dat
    # here now weboutput
    echo "$reg" >> pan_${y}.dat
    awk -v t=${TLB[$nreg]} -v d=${DWD[$nreg]} -v l=${LIT[$nreg]} -v s=${SOI[$nreg]} '{OFS="\t"; if (FNR==1) {print $1, $2, $3, $1+$2+$3} else {print $1, $2, $3, $1+$2+$3" LPJ-GUESS" ,"\n"t,d+l,s,t+l+d+s" Pan"}} ' cpool_natural_${y}_${reg}_tot.dat >> pan_${y}.dat

#    awk '(
#    echo ""
    nreg+=1
  done
  describe_textfile pan_${y}.dat "Forest Carbon Pools compared to Pan et al. regional dataset [Pg(C)]"
  cnt+=1
done