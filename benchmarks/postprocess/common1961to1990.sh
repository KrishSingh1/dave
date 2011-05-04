#!/bin/bash
#
# Common post processing for the benchmarks where simulation years 560 to 589
# correspond to real years 1961 to 1990 (500 year spin up and CRU data from 
# 1901).

tslice cmass.out -o cmass1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile cmass1961to1990.txt "Gridcell PFT/species biomass, 1961-90 averages. Units: kgC m-2"

tslice lai.out -o lai1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile lai1961to1990.txt "Gridcell PFT/species LAI, 1961-90 averages. Units: kgC m-2"

tslice dens.out -o dens1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile dens1961to1990.txt "Gridcell tree PFT/species density, 1961-90 averages. Units: trees m-2"

tslice anpp.out -o anpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile anpp1961to1990.txt "Gridcell PFT/species annual NPP, 1961-90 averages. Units: kgC m-2"

tslice cflux.out -o cflux1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile cflux1961to1990.txt "Gridcell annual carbon fluxes, 1961-90 averages. Units: kgC m-2"

tslice cpool.out -o cpool1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile cpool1961to1990.txt "Gridcell carbon pools, 1961-90 averages. Units: kgC m-2"

tslice firert.out -o firert1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile firert1961to1990.txt "Gridcell annual fire return times, 1961-90 averages. Units: years"

tslice tot_runoff.out -o tot_runoff1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile tot_runoff1961to1990.txt "Gridcell annual runoff, 1961-90 averages. Units: mm yr-1"

tslice mnpp.out -o mnpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mnpp1961to1990.txt "Gridcell monthly NPP, 1961-90 averages. Units: kgC m-2"

tslice mlai.out -o mlai1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mlai1961to1990.txt "Gridcell monthly LAI, 1961-90 averages. Units: m2 m-2"

tslice mrh.out -o mrh1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mrh1961to1990.txt "Gridcell monthly heterotrophic respiration, 1961-90 averages. Units: kgC m-2"

tslice mgpp.out -o mgpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mgpp1961to1990.txt "Gridcell monthly GPP, 1961-90 averages. Units: kgC m-2"

tslice mra.out -o mra1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mra1961to1990.txt "Gridcell monthly autotrophic respiration, 1961-90 averages. Units: kgC m-2"

tslice mnee.out -o mnee1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mnee1961to1990.txt "Gridcell monthly NEE, 1961-90 averages. Units: kgC m-2"

tslice maet.out -o maet1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile maet1961to1990.txt "Gridcell monthly AET, 1961-90 averages. Units: mm"

tslice mpet.out -o mpet1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mpet1961to1990.txt "Gridcell monthly PET, 1961-90 averages. Units: mm"

tslice mevap.out -o mevap1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mevap1961to1990.txt "Gridcell monthly evaporation, 1961-90 averages. Units: mm"

tslice mintercep.out -o mintercep1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mintercep1961to1990.txt "Gridcell monthly interception, 1961-90 averages. Units: mm"

tslice mrunoff.out -o mrunoff1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mrunoff1961to1990.txt "Gridcell monthly runoff, 1961-90 averages. Units: mm"

tslice mwcont_upper.out -o mwcontupper1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mwcontupper1961to1990.txt "Gridcell monthly water content in UPPER soil layer, 1961-90 averages (0-1)"

tslice mwcont_lower.out -o mwcontlower1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
describe_textfile mwcontlower1961to1990.txt "Gridcell monthly water content in LOWER soil layer, 1961-90 averages (0-1)"

dominance lai1961to1990.txt lai1961to1990max.txt
describe_textfile lai1961to1990max.txt "Dominant PFT/Species (with respect to 1961-90 average LAI). Units: m2 m-2"

dominance cmass1961to1990.txt cmass1961to1990max.txt
describe_textfile cmass1961to1990max.txt "Dominant PFT/Species (with respect to 1961-90 average cmass). Units: kgC m-2"