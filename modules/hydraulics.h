///////////////////////////////////////////////////////////////////////////////////////
/// \file hydraulics.h
/// \brief The hydraulics module header file
///
/// Calculation of loss of hydraulic conductance from soil water potential,
/// gs and parameters (s1, s2, kmax).
/// See Sabot et al. (2022). One Stomatal Model to Rule Them All?
/// Under conditions of water stress (stomata closed), switch to consumption from
/// plant water store calculated as per Paschalis et al. (2024)
/// Also calculated predawn lwp for use in canexch.cpp stomatal conductance,
/// uses Zhou et al. 2014 version ifdavehydraulics
///
/// \author Clare Stephens
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_HYDRAULICS_H
#define LPJ_GUESS_HYDRAULICS_H

#include "guess.h"

void tree_hydraulics(const Day &day, Individual &indiv);

void grass_hydraulics(const Day &day, Individual &indiv);

void calc_soil_water_potenial_dave(Soil &soil);

double calc_weighted_psi(Soil &soil, Pft &pft);

void hydraulic_mortality(Patch &patch);

void adjust_kmax(Patch &patch);

#endif // LPJ_GUESS_HYDRAULICS_H

