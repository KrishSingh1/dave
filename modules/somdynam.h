///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Soil organic matter dynamics
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
// Header file name:      somdynam.h
// Source code file name: somdynam.cpp
// Written by:            Ben Smith
// Version dated:         2002-08-11
// Updated:               2010-11-22


// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOMDYNAM_H
#define LPJ_GUESS_SOMDYNAM_H

#include "guess.h"

void som_dynamics(Patch& patch,Pftlist& pftlist);

#endif // LPJ_GUESS_SOMDYNAM_H
