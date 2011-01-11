#!/bin/bash
#
# Runs gmapall (running gmap on all species/PFTs) generating maps
# for LAI and CMASS for 1961 to 1990 (see also common1961to1990.sh)
#
# All parameters to this script is passed along to gmap.

# Render the maps in scalar mode, giving each file an lai or
# cmass prefix.
gmapall lai1961to1990.txt -P lai_ -s $@
gmapall cmass1961to1990.txt -P cmass_ -s $@
