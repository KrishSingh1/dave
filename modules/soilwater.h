///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Soil hydrology and snow (version including evaporation from
//                        soil surface, based on work by Dieter Gerten, Sibyll
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
//                        Schaphoff and Wolfgang Lucht, Potsdam)
// Header file name:      soilwater.h
// Source code file name: soilwater.cpp
// Written by:            Ben Smith
// Version dated:         2002-08-11
// Updated:               2010-11-22


// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOILWATER_H
#define LPJ_GUESS_SOILWATER_H

#include "guess.h"

void soilwater(Climate& climate,Patch& patch);

#endif // LPJ_GUESS_SOILWATER_H
