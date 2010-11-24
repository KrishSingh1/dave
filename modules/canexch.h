<<<<<<< .working
///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Vegetation-atmosphere exchange of H2O and CO2 via
//                        production, respiration and evapotranspiration
//                        *************************************************************
//                        "Fast" version, revised December 2002 by Ben Smith
//                        Modified according to code changes by Dieter Gerten 021216
//                        *************************************************************
// Header file name:      canexch.h
// Source code file name: canexch.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

void canopy_exchange(Patch& patch);
=======
///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Vegetation-atmosphere exchange of H2O and CO2 via
//                        production, respiration and evapotranspiration
//                        *************************************************************
//                        "Fast" version, revised December 2002 by Ben Smith
//                        Modified according to code changes by Dieter Gerten 021216
//                        *************************************************************
// Header file name:      canexch.h
// Source code file name: canexch.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
// Updated:               2010-11-22


// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_CANEXCH_H
#define LPJ_GUESS_CANEXCH_H

#include "guess.h"

void canopy_exchange(Patch& patch);

#endif // LPJ_GUESS_CANEXCH_H
>>>>>>> .merge-right.r118
