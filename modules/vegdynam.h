///////////////////////////////////////////////////////////////////////////////////////
/// \file vegdynam.h
/// \brief Vegetation dynamics and disturbance
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_VEGDYNAM_H
#define LPJ_GUESS_VEGDYNAM_H

#include "guess.h"

void vegetation_dynamics(Stand& stand,Patch& patch);
void vegetation_dynamics_daily(Patch& patch, Climate& climate);
void planting(Patch& patch);
double calc_temp_dormancy(double temp_dormancy_base, double mtemp_min20);

#endif // LPJ_GUESS_VEGDYNAM_H
