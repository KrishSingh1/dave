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
#include "growth.h"
#include "inputmodule.h"

///	Creates stands for landcovers present in the gridcell first year of the simulation
void landcover_init(Gridcell& gridcell, InputModule* input_module);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module);

/// Monitors climate history relevant for sowing date calculation. Calculates initial sowing dates/windows
void crop_sowing_gridcell(Gridcell& gridcell);

/// handles sowing date calculations for crop pft:s on patch level
void crop_sowing_patch(Patch& patch);

void phu_init(cropphen_struct& ppftcrop, Gridcellpft& gridcellpft, Patch& patch);
//void crop_phenology(Pft& pft, Patch& patch);
void crop_phenology(Patch& patch);
/// Updates crop phen from yesterday's lai_daily
void leaf_phenology_crop(Pft& pft, Patch& patch);
/// Transfers patchpft.cropphen lai and fpc-values to individuals
void update_indiv_lai_fpc(Patch& patch);
/// Updates patch.members fpc_total and fpc_rescale for crops (to be called after crop_phenology())
void update_patch_fpc(Patch& patch);
/// Calculates crop fpar for crops
void fpar_crop(Patch& patch);
/// Handles daily crop allocation and daily lai calculation
void crop_growth_daily(Patch& patch);
/// Transfer of this year's growth (ycmass_xxx) to cmass_xxx_inc
void growth_crop_year(double cmass_leaf,double cmass_root,double cmass_ho,double cmass_agpool,
	double& cmass_leaf_inc,double& cmass_root_inc,double& cmass_ho_inc,double& cmass_agpool_inc);
/// Yield function for true crops and intercrop grass
void yield_crop(Individual& indiv);
/// Yield function for pasture grass grown in cropland landcover
void yield_pasture(Individual& indiv, double cmass_leaf_inc);
/// Harvest function for cropland, including true crops, intercrop grass 
void harvest_crop(double& cmass_leaf,double& cmass_root,double& cmass_ho,double& cmass_agpool,
	double& nmass_leaf,double& nmass_root,double& nmass_ho,double& nmass_agpool,
	double& nmass_litter_leaf,double& nmass_litter_root,double& anflux_harvest,double& harvested_products_slow_nmass, double& retransn,
	double& litter_leaf,double& litter_root, double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
/// Harvest function used for clearing natural vegetation at land use change.
void harvest_natural(double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,double& cmass_debt,
	double& nmass_leaf,double& nmass_root,double& nmass_sap,double& nmass_heart,
	double& nmass_litter_leaf,double& nmass_litter_root, double& nmass_litter_sap, double& nmass_litter_heart, double& anflux_harvest,double& harvested_products_slow_nmass,
	double& litter_leaf,double& litter_root,double& litter_sap,double& litter_heart,double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
/// Harvest function for pasture, representing grazing (previous year).
void harvest_pasture(double& cmass_leaf,double& cmass_root,
	double& nmass_leaf,double& nmass_root,
	double& nmass_litter_leaf,double& nmass_litter_root,double& anflux_harvest,double& harvested_products_slow_nmass, double& retransn,
	double& litter_leaf,double& litter_root,double& acflux_harvest,double& harvested_products_slow,Individual& indiv);

/// struct storing carbon, nitrogen and water during landcover change
struct landcover_change_transfer {

	double *transfer_litter_leaf;
	double *transfer_litter_sap;
	double *transfer_litter_heart;
	double *transfer_litter_root;
	double *transfer_litter_repr;
	double *transfer_harvested_products_slow;
	double *transfer_nmass_litter_leaf;
	double *transfer_nmass_litter_sap;
	double *transfer_nmass_litter_heart;
	double *transfer_nmass_litter_root;
	double *transfer_harvested_products_slow_nmass;

	double transfer_acflux_harvest;
	double transfer_anflux_harvest;

	double transfer_cpool_fast;
	double transfer_cpool_slow;
	double transfer_wcont[NSOILLAYER];
	double transfer_wcont_evap;
	double transfer_decomp_litter_mean;
	double transfer_k_soilfast_mean;
	double transfer_k_soilslow_mean;
	Sompool transfer_sompool[NSOMPOOL];
	double transfer_nmass_avail;
	double transfer_snowpack;
	double transfer_snowpack_nmass;


	landcover_change_transfer();
	~landcover_change_transfer();
	void allocate();
};

#endif // LPJ_GUESS_LANDCOVER_H
