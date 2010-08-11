#!/bin/bash
tslice lai.out -o lai500to509.txt -f 500 -t 509 -lon 1 -lat 2 -y 3
/home/joe/work/benchmarks/gmapall/gmapall lai500to509.txt
