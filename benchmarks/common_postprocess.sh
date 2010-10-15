#!/bin/bash
tslice cmass.out -o cmass1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice lai.out -o lai1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice dens.out -o dens1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice anpp.out -o anpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice cflux.out -o cflux1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice cpool.out -o cpool1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice firert.out -o firert1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice tot_runoff.out -o tot_runoff1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mnpp.out -o mnpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mlai.out -o mlai1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mrh.out -o mrh1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mgpp.out -o mgpp1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mra.out -o mra1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mnee.out -o mnee1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice maet.out -o maet1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mpet.out -o mpet1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mevap.out -o mevap1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mintercep.out -o mintercep1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mrunoff.out -o mrunoff1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mwcont_upper.out -o mwcontupper1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
tslice mwcont_lower.out -o mwcontlower1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3
gmapall lai1961to1990.txt
dominance lai1961to1990.txt lai1961to1990max.txt
gmapall cmass1961to1990.txt
dominance cmass1961to1990.txt cmass1961to1990max.txt
