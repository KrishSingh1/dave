#!/bin/bash

#AGB percentages following Jackson et al. 1996
BF=0.75        # boreal forest
CROPS=0.9
DESERT=0.18
ScS=0.45       # Sclerophyllous Shrubs
TeCF=0.85      # Temp coniferous forest
TeDF=0.81      # Temp deciduous forest
TeG=0.21       # Temp grassland
TrDF=0.75      # Tropical deciduous forest
TrEF=0.84      # Tropical evergreen forest
TrGS=0.59      # Tropical grassland savanna
TU=0.13        # Tundra
# Hickler Coding -> Jackson Coding:
#  1 Boreal decid forest  -> BF
#  2 Boreal ever forest   -> BF
#  3 Temp/boreal mix fo.  -> TeDF
#  4 Temp conifer forest  -> TeCF
#  5 Temp decid forest    -> TeDF
#  6 Temp broad ever fo.  -> TeDF
#  7 Temp mixed forest    -> 1/2(TeDF+TeCF)
#  8 Trop season forest   -> 1/2(TrDF+TrEF)
#  9 Trop rain forest     -> TrEF
# 10 Trop decid forest    -> TrDF
# 11 Moist savannas       -> TrGS 
# 12 Dry savannas         -> TrGS 
# 13 Tall grassland       -> TeG
# 14 Dry grassland        -> TeG
# 15 Xeric wood/shrub     -> ScS
# 16 Arid shrub/steppe    -> ScS
# 17 Desert               -> DESERT
# 18 Arctic/alpine tundra -> TU 

# perform averaging over Liu time  
tslice cmass.out -f 1993 -t 2012 -o cmass1993-2012.dat
tslice landuse/landuse_hurtt_1901_2006_global.txt -f 1993 -t 2012 -o lu_1993-2012.dat
tslice lai.out   -f 1993 -t 2012 -o lai_1993-2012.dat

# remove crops and pasture from lai.out (i.e. use only first 16 cols)
awk '{ORS=" ";for (i=1;i<=16; i++) print $i; print "\n"}' lai_1993-2012.dat > lai_nat_1993-2012.dat

# get dominat PFT and compute biomes
biomes lai_nat_1993-2012.dat  

root=$(dirname $0)/..
#choose root2shoot depending on biome and Jackson Coding (see header)
awk -f $root/postprocess/agb.awk biomes_lai_nat_1993-2012.dat > agb.dat

# paste agb-fractions into cmass file
awk '{OFS="\t"; {print $23, $24, $25, $26}}' cmass1993-2012.dat | paste agb.dat - > cmass_agb.dat
awk '{OFS="\t"; if (FNR==1){print "pasture_agb_frac"} else {if($2<24. && $2>-24.){print "0.59"} else {print "0.21"}}}' lu_cmass_agb_1993-2012.dat | paste cmass_agb.dat - > cmass_agb1.dat

# loyn landuse, cmass and agb-frac
joyn lu_1993-2012.dat cmass_agb1.dat -o lu_cmass_agb_1993-2012.dat
# now compute gridcell wide contribution of each landuse type 
compute lu_cmass_agb_1993-2012.dat -i 'Crop_agb=Crop_sum*CROPLAND*0.9' 'Pasture_agb=Pasture_sum*PASTURE*pasture_agb_frac' 'Natural_agb=Natural_sum*NATURAL*agb_frac' 'Total_agb=Crop_sum*CROPLAND*0.9+Pasture_sum*PASTURE*pasture_agb_frac+Natural_sum*NATURAL*agb_frac' -o lu_cmass_agb_1993-2012_tot.dat
