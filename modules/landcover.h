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

struct Harvest_CN;

///	Creates stands for landcovers present in the gridcell first year of the simulation
void landcover_init(Gridcell& gridcell, InputModule* input_module);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module);

/// Updates dynamic management options each year
void getmanagement(Gridcell& gridcell, InputModule* input_module);

/// Monitors climate history relevant for sowing date calculation. Calculates initial sowing dates/windows
void crop_sowing_gridcell(Gridcell& gridcell);

/// handles sowing date calculations for crop pft:s on patch level
void crop_sowing_patch(Patch& patch);

void phu_init(cropphen_struct& ppftcrop, Gridcellpft& gridcellpft, Patch& patch);
//void crop_phenology(Pft& pft, Patch& patch);
void crop_phenology(Patch& patch);
/// Updates crop phen from yesterday's lai_daily
void leaf_phenology_crop(Pft& pft, Patch& patch);
/// Updates patch.members fpc_total and fpc_rescale for crops (to be called after crop_phenology())
void update_patch_fpc(Patch& patch);
/// Handles daily crop allocation and daily lai calculation
void crop_growth_daily(Patch& patch);
/// Updates crop rotation status
void crop_rotation(Stand& stand, int firsthistyear);
/// Transfer of this year's growth (ycmass_xxx) to cmass_xxx_inc
void growth_crop_year(Individual& indiv, double& cmass_leaf_inc,double& cmass_root_inc,double& cmass_ho_inc,double& cmass_agpool_inc);
/// Yield function for true crops and intercrop grass
void yield_crop(Individual& indiv);
/// Yield function for pasture grass grown in cropland landcover
void yield_pasture(Individual& indiv, double cmass_leaf_inc);
/// Harvest function for cropland, including true crops, intercrop grass 
void harvest_crop(Harvest_CN& indiv_cp, Pft& pft, bool alive, bool isintercropgrass);
/// Harvest function for cropland, including true crops, intercrop grass 
void harvest_crop(Individual& indiv, Pft& pft, bool alive, bool isintercropgrass, bool harvest_grs);
/// Harvest function used for managed forest and for clearing natural vegetation at land use change.
void harvest_wood(Harvest_CN& indiv_cp,Pft& pft, bool alive, double frac_cut);
/// Harvest function used for managed forest and for clearing natural vegetation at land use change.
void harvest_wood(Individual& indiv,Pft& pft, bool alive, double frac_cut);
/// Harvest function for pasture, representing grazing.
void harvest_pasture(Harvest_CN& indiv_cp, Pft& pft, bool alive);
/// Harvest function for pasture, representing grazing.
void harvest_pasture(Individual& indiv, Pft& pft, bool alive);
void harvest_forest(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed);
/// Turnover function for continuous grass.
void turnover_grass(Individual& indiv);
/// Transfers all carbon and nitrogen from living tissue to litter.
void kill_remaining_vegetation(Harvest_CN& indiv_cp, Pft& pft, bool alive, bool istruecrop_or_intercropgrass, bool burn);
/// Step n days from a date.
int stepfromdate(int day, int step);

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

struct Harvest_CN {

	double cmass_leaf;
	double cmass_root;
	double cmass_sap;
	double cmass_heart;
	double cmass_debt;
	double cmass_ho;
	double cmass_agpool;
	double nmass_leaf;
	double nmass_root;
	double nmass_sap;
	double nmass_heart;
	double nmass_ho;
	double nmass_agpool;
	double nstore_longterm;
	double nstore_labile;
	double max_n_storage;

	double litter_leaf;
	double litter_root;
	double litter_sap;
	double litter_heart;
	double nmass_litter_leaf;
	double nmass_litter_root;
	double nmass_litter_sap;
	double nmass_litter_heart;
	double acflux_harvest;
	double anflux_harvest;
	double harvested_products_slow;
	double harvested_products_slow_nmass;

	Harvest_CN() {

		cmass_leaf = cmass_root = cmass_sap = cmass_heart = cmass_debt = cmass_ho = cmass_agpool = 0.0;
		nmass_leaf = nmass_root = nmass_sap = nmass_heart = nmass_ho = nmass_agpool = nstore_longterm = nstore_labile = max_n_storage = 0.0;
		litter_leaf = litter_root = litter_sap = litter_heart = 0.0;
		nmass_litter_leaf = nmass_litter_root = nmass_litter_sap = nmass_litter_heart = 0.0;
		acflux_harvest = anflux_harvest = 0.0;
		harvested_products_slow = harvested_products_slow_nmass = 0.0;
	}

	/// Copies C and N values from individual and patchpft tp struct. 
	void copy_from_indiv(Individual& indiv, bool copy_grsC = false, bool copy_dead_C = true) {

		Patch& patch = indiv.vegetation.patch;
		Patchpft& ppft = patch.pft[indiv.pft.id];

		if(copy_grsC) {

			if(indiv.cropindiv) {

				cmass_leaf = indiv.cropindiv->grs_cmass_leaf;
				cmass_root = indiv.cropindiv->grs_cmass_root;

				if(indiv.pft.landcover == CROPLAND) {
					cmass_ho = indiv.cropindiv->grs_cmass_ho;
					cmass_agpool = indiv.cropindiv->grs_cmass_agpool;
				}
			}
		}
		else {

		cmass_leaf = indiv.cmass_leaf;
		cmass_root = indiv.cmass_root;
		cmass_sap = indiv.cmass_sap;
		cmass_heart = indiv.cmass_heart;
		cmass_debt = indiv.cmass_debt;

			if(indiv.pft.landcover == CROPLAND) {
				cmass_ho = indiv.cropindiv->cmass_ho;
				cmass_agpool = indiv.cropindiv->cmass_agpool;
			}
		}


		nmass_leaf = indiv.nmass_leaf;
		nmass_root = indiv.nmass_root;
		nmass_sap = indiv.nmass_sap;
		nmass_heart = indiv.nmass_heart;
		nstore_longterm = indiv.nstore_longterm;
		nstore_labile = indiv.nstore_labile;
		max_n_storage = indiv.max_n_storage;

		if(indiv.pft.landcover == CROPLAND) {
			nmass_ho = indiv.cropindiv->nmass_ho;
			nmass_agpool = indiv.cropindiv->nmass_agpool;
		}

		if(copy_dead_C) {

		litter_leaf = ppft.litter_leaf;
		litter_root = ppft.litter_root;
		litter_sap = ppft.litter_sap;
		litter_heart = ppft.litter_heart;

		nmass_litter_leaf = ppft.nmass_litter_leaf;
		nmass_litter_root = ppft.nmass_litter_root;
		nmass_litter_sap = ppft.nmass_litter_sap;
		nmass_litter_heart = ppft.nmass_litter_heart;

		// acflux_harvest and anflux_harvest only for output
		harvested_products_slow = ppft.harvested_products_slow;
		harvested_products_slow_nmass = ppft.harvested_products_slow_nmass;
	}
	}

	/// Copies C and N values from struct to individual, patchpft and patch (fluxes).
	void copy_to_indiv(Individual& indiv, bool copy_grsC = false) {

		Patch& patch = indiv.vegetation.patch;
		Patchpft& ppft = patch.pft[indiv.pft.id];

		if(copy_grsC) {

			indiv.cropindiv->grs_cmass_leaf = cmass_leaf;
			indiv.cropindiv->grs_cmass_root = cmass_root;

			if(indiv.pft.landcover == CROPLAND) {
				indiv.cropindiv->grs_cmass_ho = cmass_ho;
				indiv.cropindiv->grs_cmass_agpool = cmass_agpool;
			}
		}
		else {

		indiv.cmass_leaf = cmass_leaf;
		indiv.cmass_root = cmass_root;
		indiv.cmass_sap = cmass_sap;
		indiv.cmass_heart = cmass_heart;
		indiv.cmass_debt = cmass_debt;

			if(indiv.pft.landcover == CROPLAND) {
				indiv.cropindiv->cmass_ho = cmass_ho;
				indiv.cropindiv->cmass_agpool = cmass_agpool;
			}
		}

		indiv.nmass_leaf = nmass_leaf;
		indiv.nmass_root = nmass_root;
		indiv.nmass_sap = nmass_sap;
		indiv.nmass_heart = nmass_heart;
		indiv.nstore_longterm = nstore_longterm;
		indiv.nstore_labile = nstore_labile;

		if(indiv.pft.landcover == CROPLAND) {
			indiv.cropindiv->nmass_ho = nmass_ho;
			indiv.cropindiv->nmass_agpool = nmass_agpool;
		}

		ppft.litter_leaf = litter_leaf;
		ppft.litter_root = litter_root;
		ppft.litter_sap = litter_sap;
		ppft.litter_heart = litter_heart;
		ppft.nmass_litter_leaf = nmass_litter_leaf;
		ppft.nmass_litter_root = nmass_litter_root;
		ppft.nmass_litter_sap = nmass_litter_sap;
		ppft.nmass_litter_heart = nmass_litter_heart;

		patch.fluxes.report_flux(Fluxes::HARVESTC, acflux_harvest);
		patch.fluxes.report_flux(Fluxes::HARVESTN, anflux_harvest);

		ppft.harvested_products_slow = harvested_products_slow;
		ppft.harvested_products_slow_nmass = harvested_products_slow_nmass;
	}
};

#endif // LPJ_GUESS_LANDCOVER_H
