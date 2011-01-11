#!/bin/bash
common1961to1990.sh
common1961to1990gmapall.sh -portrait
tslice cflux.out -o cflux1990to2000.txt -f 589 -t 599 -lon 1 -lat 2 -y 3
aslice cflux1961to1990.txt -o cflux1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice cflux1990to2000.txt -o cflux1990to2000_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice cpool1961to1990.txt -o cpool1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
gmap lai1961to1990max.txt -t 'Dominant PFT (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_global.txt -o maxLAI.jpg
cbalance -spinup 500 -ncells 59191 -path ./ -start 500 -end 605
