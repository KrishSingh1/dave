#!/bin/bash
#
# Common post processing for the benchmarks where simulation years 3060 to 3089
# correspond to real years 1961 to 1990 (3000 year spin up for nitrogen version
# and CRU data from 1901).

# We will run tslice on these files                                                                              
files_to_tslice="cmass lai dens anpp cflux nflux cpool npool nleach nsources nuptake cton_leaf cton_veg firert tot_runoff \
vmaxnlim mnpp mlai mrh mgpp mra mnee maet mpet mevap mintercep mrunoff mwcont_upper mwcont_lower"

# Go through each file in the list and run tslice                                                                
for file in $files_to_tslice ; do
    tslice ${file}.out -o ${file}1961to1990.txt -f 3060 -t 3089 -lon 1 -lat 2 -y 3
done

dominance lai1961to1990.txt lai1961to1990max.txt

dominance cmass1961to1990.txt cmass1961to1990max.txt
