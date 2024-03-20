///////////////////////////////////////////////////////////////////////////////////////
/// \file growth.h
/// \brief The growth module header file
///
/// Vegetation C allocation, litter production, tissue turnover
/// leaf phenology, allometry and growth.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GROWTH_H
#define LPJ_GUESS_GROWTH_H

#include "guess.h"

typedef enum {dormant, vegetative, reproductive} growthstage;

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv);
void leaf_phenology(Patch& patch,Climate& climate);
bool allometry(Individual& indiv); // guess2008 - now returns bool instead of void
void establishment_biomass(double cmass_init, double ltor, Individual& indiv, double nmass_init = 0.0);
void allocation_init(double bminit,double ltor,Individual& indiv);
void growth(Stand& stand,Patch& patch);
void growth_daily(Patch& patch);
void growth_daily_grass(Patch& patch);
double photo_growth(double reldaylength, double daylength_diff);
void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
	lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_sap,
	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
	double& nmass_heart, double& litter_leaf, double& litter_root, double& cmass_leaf_root_turnover,
	double& nmass_litter_leaf, double& nmass_litter_root,
	double& longterm_nstore, double &max_n_storage,
	bool alive);
// Minimum growth demand after reproductive phase
const double fdev_growth_min = 0.005;
void allometry_dave_tree(Individual& indiv);
double temp_growth(double temp, double Ha_growth, double T0_growth, double A_growth, double temp_dorm, double k);

#endif // LPJ_GUESS_GROWTH_H
