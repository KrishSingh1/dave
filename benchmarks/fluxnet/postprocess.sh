#!/bin/bash

DIR=$(grep flux_dir common/../paths.ins | cut -d\" -f4)

function compare {
  obs_idx=$(( $3 + 4 ))
  awk -v obs="$obs_idx" '{print $1,$2,$3,$4,$obs}' ${DIR}monthly.csv > tmp

  x=$1
  c=$2
  unit=$4
  awk -v c="$c" '{if(NR==1){print "lon","lat","year","month",c} else { for(i=4;i<=NF;i++){j=i-3;print $1,$2,$3,j,$i}}}' $1 > $c"_col.dat"

  awk 'BEGIN {FS = " "};FNR==NR{a[$1$2$3$4]=$5;next}BEGIN{OFS=" "};{if($1$2$3$4 in a){print $1,$2,$3,$4,$5,a[$1$2$3$4]}}' $c"_col.dat" tmp > output.txt

  more output.txt | grep -v "\-9999" > plot

  gplot plot -o $c"_scatter.png" -x 5 -y 6 -xt "Observed" -yt "Modelled" -scatter -eq -t "compared vs. modelled $c ($unit)"

  rm tmp plot output.txt
}
describe_benchmark "LPJ-GUESS - FLUXNET Benchmarks (global PFTs)"

# Convert from kgC m-2 month-1 -> gC m-2 day-1
compute mnee.out -o mnee.txt -n -i Lon Lat Year Janu="Jan/31*1000" Febru="Feb/28*1000" March="Mar/31*1000" April="Apr/30*1000" Ma="May/31*1000" June="Jun/30*1000" July="Jul/31*1000" Augu="Aug/31*1000" Sept="Sep/30*1000" Octo="Oct/31*1000" Nove="Nov/30*1000" Dece="Dec/31*1000"
compute mgpp.out -o mgpp.txt -n -i Lon Lat Year Janu="Jan/31*1000" Febru="Feb/28*1000" March="Mar/31*1000" April="Apr/30*1000" Ma="May/31*1000" June="Jun/30*1000" July="Jul/31*1000" Augu="Aug/31*1000" Sept="Sep/30*1000" Octo="Oct/31*1000" Nove="Nov/30*1000" Dece="Dec/31*1000"

compare mnee.txt NEE 2 "gC m-2 d-1"
compare mgpp.txt GPP 3 "gC m-2 d-1"
compare maet.out AET 4 "w m-2"
