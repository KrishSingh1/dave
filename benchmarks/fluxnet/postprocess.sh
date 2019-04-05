#!/bin/sh

# Function that flips rows of months to columns
#
# Parameters:
# $1 input filename
# $2 output filename
# $3 variable name
function flipMonthlyData {
  awk -v c=$3 '{if(NR==1){print "Lon","Lat","Year","Month", c} else { for(i=4;i<=NF;i++){j=i-3;print $1,$2,$3,j,$i}}}' $1 > $2
}

# Function that selects two columns from an input file
#
# Parameters:
# $1 input filename
# $2 output filename
# $3 obsvar
# $4 modvar
function selectData {
  obs=$(head -1 $1 | awk -v d=$3 '{for(i=1;i<=NF;i++){if($i==d){print i}}}')
   mod=$(head -1 $1 | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i}}}')

   awk -v mod=$mod -v obs=$obs '{if(FNR==1){next;}else { print $obs,$mod }}' $1 > $2
}

# Function that produce a scatterplot
#
# Parameters:
# $1 input data
# $2 output filename
# $3 variable
# $4 unit
function produceScatterplot {
  # remove missing data from the fluxnet dataset
  grep "\-9999" $1 -v  > plot

  gplot plot -o $2 -x 1 -y 2 -xt "Observed" -yt "Modelled" -scatter -eq -t "compared vs. modelled $3 ($4)"
  describe_image $2 "Observed fluxnet data vs modelled $3: Units: $4" embed
}

describe_benchmark "LPJ-GUESS - Fluxnet benchmark"

DIR=$(grep flux_dir common/../paths.ins | cut -d\" -f4)

# Convert from kgC m-2 month-1 -> gC m-2 day-1
compute mnee.out -o mnee.txt -n -i Lon Lat Year Janu="Jan/31*1000" Febru="Feb/28*1000" March="Mar/31*1000" April="Apr/30*1000" Ma="May/31*1000" June="Jun/30*1000" July="Jul/31*1000" Augu="Aug/31*1000" Sept="Sep/30*1000" Octo="Oct/31*1000" Nove="Nov/30*1000" Dece="Dec/31*1000"
compute mgpp.out -o mgpp.txt -n -i Lon Lat Year Janu="Jan/31*1000" Febru="Feb/28*1000" March="Mar/31*1000" April="Apr/30*1000" Ma="May/31*1000" June="Jun/30*1000" July="Jul/31*1000" Augu="Aug/31*1000" Sept="Sep/30*1000" Octo="Oct/31*1000" Nove="Nov/30*1000" Dece="Dec/31*1000"

# Flip monthly data from one column per month to Lon Lat Year Month var
flipMonthlyData mnee.txt nee_col.txt NEEmod
flipMonthlyData mgpp.txt gpp_col.txt GPPmod
flipMonthlyData maet.out aet_col.txt aet
flipMonthlyData mintercep.out intercep_col.txt intcpt

# Add AET and intercept to get corresponding LE value in Fluxnet
joyn aet_col.txt intercep_col.txt -o tmp -i 1 2 3 4
compute tmp -o le_col.txt -n -i Lon Lat Year Month LEmod="aet+intcpt"

# Join all dataframes into one big dataframe
joyn ${DIR}/monthly.txt nee_col.txt -o nee_joyn.txt -i 1 2 3 4
joyn nee_joyn.txt le_col.txt -o aet_joyn.txt -i 1 2 3 4
joyn aet_joyn.txt gpp_col.txt -o all_data.txt -i 1 2 3 4

# Iterate over the gridlist and plot each station separately
while read f; do
   lon=$(echo $f | awk '{ print $1 }')
   lat=$(echo $f | awk '{ print $2 }')
   site=$(echo $f | awk '{ print $3 }')

   # Extract and select data
   extract all_data.txt -o ${site}.csv -x "Lon==${lon} && Lat==${lat}"
   selectData ${site}.csv nee_scatter NEE NEEmod
   selectData ${site}.csv gpp_scatter GPP GPPmod
   selectData ${site}.csv le_scatter LE LEmod

   # Plot variables
   produceScatterplot nee_scatter "${site}_nee.png" "NEE" "gC m-2 day-1"
   produceScatterplot gpp_scatter "${site}_gpp.png" "GPP" "gC m-2 day-1"
   produceScatterplot le_scatter "${site}_le.png" "LE" "w m-2"

   # Cleanup
   rm nee_scatter gpp_scatter le_scatter plot ${site}.csv
done < gridlist.txt

# Cleanup
rm tmp
