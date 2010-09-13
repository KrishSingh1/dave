#!/bin/bash
tslice cmass.out -o cmass1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice lai.out -o lai1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice dens.out -o dens1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice anpp.out -o anpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice cflux.out -o cflux1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice cpool.out -o cpool1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice firert.out -o firert1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
/home/joe/work/benchmarks/gmapall/gmapall lai1961to1990.txt
/home/joe/work/benchmarks/dominance/dominance lai1961to1990.txt lai1961to1990max.txt
