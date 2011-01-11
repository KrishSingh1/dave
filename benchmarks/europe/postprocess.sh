#!/bin/bash
common1961to1990.sh
common1961to1990gmapall.sh
gmap lai1961to1990max.txt -t 'Dominant Species (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxLAIeurope.jpg
gmap cmass1961to1990max.txt -t 'Dominant Species (greatest cmass)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxCMASSeurope.jpg
cbalance -spinup 500 -ncells 3464 -path ./ -start 500 -end 605
