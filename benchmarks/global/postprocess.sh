#!/bin/bash
describe_benchmark "LPJ-GUESS - Global Benchmarks"

common1961to1990.sh
common1961to1990gmapall.sh -portrait

tslice cflux.out -o cflux1990to2000.txt -f 589 -t 599 -lon 1 -lat 2 -y 3
aslice cflux1961to1990.txt -o cflux1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
describe_textfile cflux1961to1990_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1961 to 1990. Units: kgC m-2"
aslice cflux1990to2000.txt -o cflux1990to2000_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
describe_textfile cflux1990to2000_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1990 to 2000. Units: kgC m-2"

aslice cpool1961to1990.txt -o cpool1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
describe_textfile cpool1961to1990_areaaverage.txt "Global Terrestrial Carbon Pools, 1961 to 1990. Units: kgC m-2"

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
describe_textfile tot_runoff1961to1990_areaaverage.txt "Global Runoff, 1961 to 1990. Units: mm yr-1"

gmap lai1961to1990max.txt -t 'Dominant PFT (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_global.txt -o maxLAI.jpg
describe_image maxLAI.jpg "PFT With the Highest LAI in Each Gridcell (1961-90 average)" embed

cbalance -spinup 500 -ncells 59191 -path ./ -start 500 -end 605
describe_textfile cbalance_totalerror_GtC.txt "Global Terrestrial Carbon Uptake, 1901 to 2006. /
Determined using C pools (cpool_GtC), Cumulative C fluxes (cflux_GtC), and their absolute difference (absdiff_GtC)"