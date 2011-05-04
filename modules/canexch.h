///////////////////////////////////////////////////////////////////////////////////////
/// \file canexch.h
/// \brief The canopy exchange module header file
///
/// Vegetation-atmosphere exchange of H2O and CO2 via
/// production, respiration and evapotranspiration.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_CANEXCH_H
#define LPJ_GUESS_CANEXCH_H

#include "guess.h"

void canopy_exchange(Patch& patch);

// Constants for photosynthesis calculations
const double CQ=4.6E-6;      // conversion factor for solar radiation at 550
                             // nm from J/m2 to mol_quanta/m2 (E=mol quanta); mol J-1
const double ALPHA_C3=0.08;  // intrinsic quantum efficiency of CO2 uptake for C3 plants
const double ALPHA_C4=0.053; // intrinsic quantum efficiency of CO2 uptake for C4 plants
const double PO2=2.09E4;     // O2 partial pressure (Pa)

#endif // LPJ_GUESS_CANEXCH_H
