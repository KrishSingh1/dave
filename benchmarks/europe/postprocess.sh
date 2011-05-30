#!/bin/bash
describe_benchmark "LPJ-GUESS - European Benchmarks"

common1961to1990.sh
common1961to1990gmapall.sh

gmap lai1961to1990max.txt -t 'Dominant Species (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxLAIeurope.jpg
describe_image maxLAIeurope.jpg "Species With the Highest LAI in Each Gridcell (1961-90 average)"

gmap cmass1961to1990max.txt -t 'Dominant Species (greatest cmass)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxCMASSeurope.jpg
describe_image maxCMASSeurope.jpg "Species With the Highest CMASS in Each Gridcell (1961-90 average)"

cbalance -spinup 500 -ncells 3464 -path ./ -start 500 -end 605
describe_textfile cbalance_totalerror_GtC.txt "European Terrestrial Carbon Uptake, 1901 to 2006. /
Determined using C pools (cpool_GtC), Cumulative C fluxes (cflux_GtC), and their absolute difference (absdiff_GtC)"