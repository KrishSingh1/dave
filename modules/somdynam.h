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
//RLN because they are needed in blaze.cpp
//WK Explain what they are and what they are used for. 
// computes the fraction of leaf and root that goes to metabolic litter (used by BLAZE)
double metabolic_litter_fraction(double lton);

// determine lignin to N ratio for leaf and root litter (used by BLAZE)
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr);

void som_dynamics(Patch& patch);

#endif // LPJ_GUESS_SOMDYNAM_H
