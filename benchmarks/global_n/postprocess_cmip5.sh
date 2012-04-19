#!/bin/bash

mkdir results

common1961to1990.sh
tslice cflux.out -o results/cflux1990to2000.txt -f 589 -t 599 -lon 1 -lat 2 -y 3
aslice cflux1961to1990.txt -o results/cflux1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice cflux1990to2000.txt -o results/cflux1990to2000_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice cpool1961to1990.txt -o results/cpool1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
aslice tot_runoff1961to1990.txt -o results/tot_runoff1961to1990_areaaverage.txt -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.25 0.25   
gmap lai1961to1990max.txt -t 'Dominant PFT (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_global.txt -o results/maxLAI.jpg
cbalance -spinup 500 -ncells 59191 -path ./ -start 500 -end 605

#gunzip cpool.out.gz
#gunzip npool.out.gz
#gunzip nsources.out.gz
sed 's/nan/0/g' nsources.out > tmp_n.txt
#gunzip lai.out.gz

aslice cpool.out -o 'results/cpool_area.txt' -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice npool.out -o 'results/npool_area.txt' -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice tmp_n.txt -o 'results/nsources_area.txt' -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25

# N sources
tslice tmp_n.txt -o tmp1.txt -f 500 -t 505
tslice tmp_n.txt -o tmp2.txt -f 745 -t 750

awk '{$1=$1}1' tmp1.txt > 'nsources_1850.txt' 
mv 'nsources_1850.txt' results/
awk '{$1=$1}1' tmp2.txt > 'nsources_2100.txt'
mv 'nsources_2100.txt' results/

# N sources
tslice cpool.out -o tmp3.txt -f 500 -t 505
tslice cpool.out -o tmp4.txt -f 745 -t 750

awk '{$1=$1}1' tmp3.txt > 'cpool_1850.txt' 
mv 'cpool_1850.txt' results/
awk '{$1=$1}1' tmp4.txt > 'cpool_2100.txt'
mv 'cpool_2100.txt' results/
	
rm -rf tmp*

for ((j=500;j<=750;j=j+10))
do
	tslice lai.out -o 'lai_'$j'.txt' -f $j -t $j
	/home/david/post/biomes 'lai_'$j'.txt'
	mv 'biomes_lai_'$j'.txt' results/
	rm -rf 'lai_'$j'.txt'
done





	