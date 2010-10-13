///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Vegetation dynamics and disturbance
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
// Header file name:      vegdynam.h
// Source code file name: vegdynam.cpp
// Written by:            Ben Smith
// Version dated:         2002-08-11
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_VEGDYNAM_H
#define LPJ_GUESS_VEGDYNAM_H

#include "guess.h"

//void vegetation_dynamics(Stand& stand,Patch& patch,Pftlist& pftlist);

// guess2008 - euroflux
void vegetation_dynamics_euroflux(Stand& stand,Patch& patch,Pftlist& pftlist,int century_year);


#endif // LPJ_GUESS_VEGDYNAM_H
