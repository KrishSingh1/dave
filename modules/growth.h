///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Vegetation C allocation, litter production, tissue turnover
//                        leaf phenology, allometry and growth
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
// Header file name:      growth.h
// Source code file name: growth.cpp
// Written by:            Ben Smith
// Version dated:         2002-08-11
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GROWTH_H
#define LPJ_GUESS_GROWTH_H

#include "guess.h"

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv);
void leaf_phenology(Patch& patch,Climate& climate);
bool allometry(Individual& indiv); // void to bool - guess2008 - 080827
void allocation_init(double bminit,double ltor,Individual& indiv);
void growth(Stand& stand,Patch& patch);

#endif // LPJ_GUESS_GROWTH_H
