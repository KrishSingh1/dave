///////////////////////////////////////////////////////////////////////////////////////
/// \file growth.h
/// \brief The growth module header file
///
/// Vegetation C allocation, litter production, tissue turnover
/// leaf phenology, allometry and growth.
///
/// \author Ben Smith
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GROWTH_H
#define LPJ_GUESS_GROWTH_H

#include "guess.h"
#include "equationsolvers.h"

typedef enum {dormant, vegetative, reproductive} growthstage;

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv);
void leaf_phenology(Patch& patch,Climate& climate);
bool allometry(Individual& indiv); // guess2008 - now returns bool instead of void
void establishment_biomass(double cmass_init, double ltor, Individual& indiv, double nmass_init = 0.0);
void allocation_init(double bminit,double ltor,Individual& indiv,double patch_fpc = 0);
void growth(Stand& stand,Patch& patch);
void growth_daily(Patch& patch);
void growth_daily_grass_warlind(Patch& patch);
double photo_growth(double reldaylength, double daylength_diff);
double calc_dphen_activity_max(double dphen_activity_base, double dphen_activity_pot,
	   double k, size_t window_length, size_t window_length_min);
void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
	lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_sap,
	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
	double& nmass_heart, double& litter_leaf, double& litter_root, double& cmass_leaf_root_turnover,
	double& nmass_litter_leaf, double& nmass_litter_root,
	double& longterm_nstore, double& max_n_storage,
	bool alive, double& NH4_mass);
// Minimum growth demand after reproductive phase
const double fdev_growth_min = 0.002;
void allometry_dave_tree(Individual& indiv);
double cstorage_capacity(const Individual& indiv);
double temp_growth(double temp, double Ha_growth, double T0_growth, double A_growth, double T_dormancy);
double wscal_growth(double wscal, double wscal0, double b);
double get_cmass_init_min(Patchpft& ppft);

// ---------------------------------
// Grass parameters
// ---------------------------------
// Minimum lai below which the plant is either killed or cannot establish.
const double lai_min = 0.01;

// Maximum size of C storage pool as a fraction of biomass (pool-dependent!)
const double cstorage_max_temp = 0.3;   // leaves and roots
const double cstorage_max_perm = 0.75;  // crowns

// Water response parameters for growth and developmental processes
const double wscal0 = 0.7;     // inflection point
const double bgrowth = 0.07;   // steepness around wscal0

#endif // LPJ_GUESS_GROWTH_H
