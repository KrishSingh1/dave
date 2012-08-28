///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H

#include "guess.h"

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell);

void check_crop_temp_limits(Climate& climate, Gridcellpft& gridcellpft);
void calc_crop_dates_20y_mean(Climate& climate, Gridcellpft& gridcellpft);
void set_sdatecalc_temp(Climate& climate, Gridcellpft& gridcellpft);
void set_sdatecalc_prec(Climate& climate, Gridcellpft& gridcellpft);
void calc_sowing_windows(Gridcell& gridcell,Pftlist& pftlist);
void calc_m_climate_20y_mean(Climate& climate);
static double variation_coefficient(double data[],int n);
void calc_seasonality(Gridcell& gridcell);
void crop_sowing_gridcell(Gridcell& gridcell,Pftlist& pftlist);
void Crop_sowing_date_temp(Patch& patch, Pft& pft);
void Crop_sowing_date_prec(Patch& patch, Pft& pft);
void Crop_sowing_date_rice(Patch& patch, Pft& pft);
//void Crop_sowing_date_forced(Patch& patch, Pft& pft);
void Crop_sowing_date(Patch& patch, Pft& pft);
void Crop_sowing_date_new(Patch& patch, Pft& pft);
void crop_sowing_patch(Patch& patch, Pftlist& pftlist);

void leaf_phenology_crop(Pft& pft,Climate& climate,double wscal,double aphen, double& phen, Gridcellpft& gridcellpft, bool isirrigated, Patch& patch);
void fpar_crop(Patch& patch);
void growth_crop_daily(Patch& patch);void allocation_crop(double bminc,double cmass_leaf,double cmass_root,double cmass_ho,double ltor,double& cmass_plant_inc,double& cmass_leaf_inc,
	double& cmass_root_inc,double& cmass_ho_inc,double& cmass_agpool_inc, double& litter_leaf_inc,double& litter_root_inc, Individual& indiv);
void harvest_crop(double& cmass_plant,double& cmass_leaf,double& cmass_root,double& cmass_ho,double& cmass_agpool,double& litter_leaf,double& litter_root,
	double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
void harvest_natural(double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,double& cmass_debt,double& litter_leaf,double& litter_root,double& litter_wood,
	double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
void harvest_pasture(double& cmass_leaf,double& cmass_root,double& litter_leaf,double& litter_root,double& acflux_harvest,double& harvested_products_slow,Individual& indiv);

#endif // LPJ_GUESS_LANDCOVER_H