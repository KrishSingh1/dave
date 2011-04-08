///////////////////////////////////////////////////////////////////////////////////////
/// \file soilwater.h
/// \brief Soil hydrology and snow
///
/// Version including evaporation from soil surface, based on work by Dieter Gerten, 
/// Sibyll Schaphoff and Wolfgang Lucht, Potsdam
///
/// Version adapted for analyses for joint conceptual paper for
/// OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOILWATER_H
#define LPJ_GUESS_SOILWATER_H

#include "guess.h"

void soilwater(Climate& climate,Patch& patch);

#endif // LPJ_GUESS_SOILWATER_H
