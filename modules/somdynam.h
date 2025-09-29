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

// computes the decomposition temperature modifier (in range 0-1)
double temperature_modifier(double temp);

// computes the decomposition moisture modifier (in range 0-1)
double moisture_modifier(double wfps);

void som_dynamics(Patch& patch, Climate& climate);

// computes the fraction of leaf and root that goes to metabolic litter (used by BLAZE)
double metabolic_litter_fraction(double lton);

// determine lignin to N ratio for leaf and root litter (used by BLAZE)
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr);

void leaf_to_litter(Soil& soil, const Pft &pft, double litterme[NSOMPOOL], double fireresist[NSOMPOOL], double cmass_litter_leaf, double nmass_litter_leaf);
void root_to_litter(Soil& soil, const Pft &pft, double litterme[NSOMPOOL], double fireresist[NSOMPOOL], double cmass_litter_root, double nmass_litter_root);
void sap_to_litter(Soil& soil, const Pft &pft, double litterme[NSOMPOOL], double fireresist[NSOMPOOL], double cmass_litter_sap, double nmass_litter_sap);
void heart_to_litter(Soil& soil, const Pft &pft, double litterme[NSOMPOOL], double fireresist[NSOMPOOL], double cmass_litter_heart, double nmass_litter_heart);

void init_litterme(const Soil &soil, double litterme[NSOMPOOL], double fireresist[NSOMPOOL]);

/*
Update litterme and fireresist in the soil object iff GLOBFIRM is in use.
*/
void set_litterme(Soil &soil, double litterme[NSOMPOOL], double fireresist[NSOMPOOL]);

// Leaf, root and wood litter lignin fractions
// Leaf and root fractions: Comins & McMurtrie 1993; Friend et al 1997
const double LIGCFRAC_LEAF = 0.2;
const double LIGCFRAC_ROOT = 0.16;
const double LIGCFRAC_WOOD = 0.3;	// TODO Check wood fraction

#endif // LPJ_GUESS_SOMDYNAM_H
