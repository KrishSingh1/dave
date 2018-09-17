///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.h
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOMDYNAM_H
#define LPJ_GUESS_SOMDYNAM_H

#include "guess.h"

//WK Why do the following two need to be declared here explicitly?
//WK Explain what they are and what they are used for.

double metabolic_litter_fraction(double lton);

double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr);

void som_dynamics(Patch& patch);

#endif // LPJ_GUESS_SOMDYNAM_H
