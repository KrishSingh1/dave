///////////////////////////////////////////////////////////////////////////////////////
/// \file growth.cpp
/// \brief The growth module
///
/// Vegetation C allocation, litter production, tissue turnover
/// leaf phenology, allometry and growth
///
/// (includes updated FPC formulation as required for "fast"
/// cohort/individual mode - see canexch.cpp)
///
/// \author Ben Smith
/// $Date$
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "growth.h"
#include "canexch.h"
#include "landcover.h"
#include "ncompete.h"
#include <assert.h>
#include "somdynam.h"
#include "driver.h"


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Maximum number of equivalent days with full leaf cover per growing season
// for summergreen PFTs
const double APHEN_MAX = 210.0;

/// External seed input (kgC m-2 d-1)
const double ext_seed = 2.0e-6;

// ---------------------------------
// Grass parameters
// ---------------------------------
// Maximum size of C storage pool as a fraction of biomass (pool-dependent!)
const double cstorage_max_temp = 0.3;   // leaves and roots
const double cstorage_max_perm = 0.75;  // crowns

// Fraction of maximum C storage not available for growth.
const double frac_cstorage_res = 0.15;

// Maximum mobilisation rate of C in storage
const double max_mobil_cstorage = 0.2;

// Minimum physiological activity during dormancy (fraction)
const double fphen_activity_min = 0.1;

// Maximum physiological up/down-regulation in and out of dormancy (fraction d-1)
const double dphen_activity_max = 0.3;

// Error tolerance (expressed as number of significant figures) used in mass balance checks.
const double EPS = -12;

class DaveMassBalance {
public:
	void init(const Vegetation &vegetation);
	void check(const Vegetation &vegetation);
#ifndef NDEBUG
	~DaveMassBalance();
private:
	bool should_check(const Individual &);
	class IndivState {
public:
		IndivState(const Individual &indiv);
		void check(const Individual &indiv) const;
		int get_id() const { return id; };
private:
		double n_start;
		double c_start;
		int id;
		void check_cton(double cmass, double nmass, double cton_min, double cton_max) const;
	};
	std::vector<const IndivState*> indivs;
	void clear_indivs();
#endif // NDEBUG
};

#ifndef NDEBUG

DaveMassBalance::~DaveMassBalance() {
	clear_indivs();
}

void DaveMassBalance::clear_indivs() {
	for (const IndivState *indiv : indivs) {
		delete indiv;
	}
	indivs.clear();
}

bool DaveMassBalance::should_check(const Individual &indiv) {
	return indiv.pft.dave_allocation();
}

double get_nmass(const Individual &indiv) {
	if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
		return indiv.nmass_storage +
	       indiv.nmass_crown +
		   indiv.nmass_leafs[green] +
		   indiv.nmass_leafs[brown] +
		   indiv.nmass_root +
		   indiv.nmass_repr;
	}
	return indiv.nmass_storage +
		   indiv.nmass_leaf +
		   indiv.nmass_root +
		   indiv.nmass_sap +
		   indiv.nmass_repr +
		   indiv.nmass_heart;
}

double get_cmass(const Individual &indiv) {
	if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
		return indiv.cmass_storage +
			indiv.cmass_crown +
			indiv.cmass_leafs[green] +
			indiv.cmass_leafs[brown] +
			indiv.cmass_root +
			indiv.cmass_repr;
	}
	return indiv.cmass_storage +
		indiv.cmass_leaf +
		indiv.cmass_root +
		indiv.cmass_sap +
		indiv.cmass_repr +
		indiv.cmass_heart;
}

double get_cflux(const Individual &indiv) {
	return indiv.dnpp - indiv.ctolitter;
}

double get_nflux(const Individual &indiv) {
	return indiv.nmass_missing - indiv.ntolitter;
}

DaveMassBalance::IndivState::IndivState(const Individual &indiv) {
	n_start = get_nmass(indiv);
	c_start = get_cmass(indiv);
	id = indiv.id;
}

void DaveMassBalance::IndivState::check_cton(double cmass, double nmass, double cton_min, double cton_max) const {
	if (is_zero(cmass) || is_zero(nmass)) {
		return;
	}
	double cton = cmass / nmass;
	assert(greater_than_or_equal(cton, cton_min * 0.999));
	// assert(cton >= cton_min);
 	assert(less_than_or_equal(cton, cton_max * 1.001));
	// assert(cton <= cton_max);
}

void DaveMassBalance::IndivState::check(const Individual &indiv) const {
	double c_flux = get_cflux(indiv);
	double c_end = get_cmass(indiv) + indiv.cexcess;
	assert(equal(c_end, c_start + c_flux));

	double n_flux = get_nflux(indiv);
	double n_end = get_nmass(indiv);
	assert(equal(n_end, n_start + n_flux));

	check_cton(indiv.get_cmass_leaf(), indiv.get_nmass_leaf(), indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max);
	check_cton(indiv.cmass_root, indiv.nmass_root, indiv.pft.cton_root_min, indiv.pft.cton_root_max);
	check_cton(indiv.cmass_repr, indiv.nmass_repr, indiv.pft.cton_repr_min, indiv.pft.cton_repr_max);

	if (!indiv.pft.is_grass()) {
		check_cton(indiv.cmass_sap, indiv.nmass_sap, indiv.pft.cton_sap_min, indiv.pft.cton_sap_max);
	}
}

#endif // NDEBUG

void DaveMassBalance::init(const Vegetation &vegetation) {
#ifndef NDEBUG
	clear_indivs();
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &individual = vegetation[i];
		if (should_check(individual)) {
			IndivState *indiv = new IndivState(individual);
			indivs.push_back(indiv);
		}
	}
#endif // NDEBUG
}

void DaveMassBalance::check(const Vegetation &vegetation) {
#ifndef NDEBUG
	for (int i = 0; i < vegetation.nobj && i < indivs.size(); i++) {
		// If an individual is killed during the day, it will be removed from
		// the vegetation object. Ergo, we cannot rely on the order or size of
		// the vegetation array matching the order and size of our internal
		// vector. To allow for this we also check the indiv ID.
		Individual &indiv = vegetation[i];
		if (should_check(indiv) && indivs[i]->get_id() == indiv.id) {
			indivs[i]->check(indiv);
		}
	}
#endif
}


///////////////////////////////////////////////////////////////////////////////////////
// LEAF PHENOLOGY AND GROWTH
// Call function leaf_phenology each simulation day prior to calculation of FPAR, to
// calculate fractional leaf-out for each PFT and individual.
// Function leaf_phenology_pft is not intended to be called directly by the framework,


// Temperature scalar used for growth and dormancy
// Relationship and parameters from Parent & Tardieu 2012, New Phytologist
// representative species:
// 		C3s: wheat (similar to barley, A. thaliana, canola)
// 		C4s: sorghum (similar to millet, maize)
// Function includes a modified response at temperatures lower than dormancy temperature.
// from Fatichi et al. 2014 (New Phytologist): 
// "a 5–6°C threshold has often been identified to limit growth in cold adapted
// species, irrespective of photosynthetic activity which typically
// ceases only at freezing point (Koerner, 2008)"
double temp_growth(double temp, double Ha_growth, double T0_growth, double A_growth, double temp_dorm) {
	
	double alpha   = 3.5;                   // unitless
	double Ha      = Ha_growth * J_PER_KJ;  // J mol-1
	double T0      = T0_growth + K2degC;    // K
	double Tleaf   = temp + K2degC;         // K   
	double Tdorm   = temp_dorm + K2degC;    // K

    // Note: A_growth was rescaled so that temp_growth ranges from 0-1
	double ftemp;

	ftemp = A_growth * (Tleaf * exp(-Ha / (R_gas * Tleaf))) / 
			(1.0 + pow(exp(-Ha / (R_gas * Tleaf)),(alpha * (1.0 - Tleaf / T0))) );

	// low temperature adjustment to avoid simulating unrealistic growth demand
    // at very low or even negative temperatures (see also comment above).
	// When above tdorm, ftemp = ftemp
	// When below tdorm, ftemp *= 1 - tdef * 0.5
	// So, when 2 degrees below tdorm, ftemp = 0
	//     when 0 degrees below tdorm, ftemp = ftemp
    double Tdef = Tdorm - Tleaf;
    ftemp = Tdef < 0.0 ? ftemp : ftemp * min(1.0, max(0.0, (1.0 - Tdef * 0.5)));	// 0.5 for 0 growth 2degC below Tdorm.	

	return max(0.0,min(1.0,ftemp));
}

// Water stress scalar used for growth and dormancy
double wscal_growth(double wscal, double q) {
	double wscal_min = 0.4;  // no growth under severe water stress (ideally consistent with pft.wscal_min)
	return min(1.0,max(0.0,pow(((wscal - wscal_min) / (1.0 - wscal_min)),q)));
}

// Photoperiod affecting development and thus growth demand
// Effects likely limited to high latitude areas.
double photo_growth(double reldaylength, double daylength_diff) {
	return 1.0;
}


// Calculates dphen_activity, which is used for calculation of:
// 1) the rate at which grasses enter or release dormancy (dormancy = true).
// 2) the rate at which grasses resume growth after completion of reproductive cycle (dormancy = false).
double dphen_activity(Individual &indiv, Climate& climate, bool dormancy) {

	Pft& pft = indiv.pft;
	Patchpft &ppft = indiv.patchpft();
	
	// different water stress thresholds for dormancy and vegetative regrowth.
	double wscal_min = dormancy ? pft.wscal_min : 0.6;

	// Threshold temperature for dormancy (base temperature same as dormancy threshold, depends on photosyn. pathway)
	double grass_temp = indiv.vegetation.patch.soil.T_soil[7];   // Temperature of first soil layer
	double wscal_mean = ppft.wscal_10.size() > 0 ? ppft.wscal_10.mean() : 1.0;
	double tscal_crit = temp_growth(indiv.pft.temp_dormancy,indiv.pft.Ha_growth,indiv.pft.T0_growth,indiv.pft.A_growth,indiv.pft.temp_dormancy);

	//double wscal_actual = wscal_growth(wscal_mean);
	double wscal_actual = wscal_mean;
	double tscal_actual = temp_growth(ppft.soiltemp_10.mean(),indiv.pft.Ha_growth,indiv.pft.T0_growth,indiv.pft.A_growth,indiv.pft.temp_dormancy);

	// Changes in physiological activity due to changes in dormancy status
	double denom_wscal = wscal_actual >= wscal_min ? pow((1.0 - wscal_min), 1.5) : wscal_min;
	double num_wscal   = wscal_actual >= wscal_min ? pow((wscal_actual - wscal_min), 1.5) : (wscal_actual - wscal_min);

	double tscal = tscal_actual > tscal_crit ? tscal_actual : (tscal_actual - tscal_crit) / tscal_crit;
	tscal = tscal > 0.0 ? pow(tscal, 0.5) : tscal;
	double wscal = num_wscal / denom_wscal;

	assert(tscal >= -1.0 && tscal <= 1.0);
	assert(wscal >= -1.0 && wscal <= 1.0);

	// avoid greening just before winter in temperature-limited ecosystems.
	// TODO: interactions with tscal and pscal?
	double pscal = 1.0;

	if (climate.mtemp_min20 < (indiv.pft.temp_dormancy - 5.0) &&
	    ppft.daylength_10.periodictrend(1) < 0.0 && !dormancy) {
			
			double daylength_crit = 0.9;
			double daylength_min  = 0.6;

			pscal = 2.0 * min(1.0,max(0.0,(ppft.reldaylength - daylength_min) / 
				    				      (daylength_crit - daylength_min))) - 1.0;
	}

	return min(min(tscal,pscal),wscal) * dphen_activity_max;
}


// Development effects on relative growth rate.
// Reduction of growth demand in reproductive phase.
void dev_growth(Individual &indiv, Climate &climate, double frepr_cstruct) {
	
	Patchpft &ppft = indiv.patchpft();

	if (indiv.growth_stage == reproductive || indiv.growth_stage == dormant) {

		// TODO: better (and more interpretable) characterisation of a1 parameter.
		// TODO: check if a2 parameter is needed.
		// TODO: seeding too early when reprfrac is high.
		//  a1: higher values: later (at higher frepr_cstruct) reduction in growth demand.
		//   --> a1 should decrease with reprfrac. Same reprfrac --> earlier seeding with higher a1.
		// a1 parameter should depend on the size of the plant at flowering time.
		double a1m = indiv.pft.lifeform == GRASS_ANNUAL ? 1.75 : 0.03;  
		double a1 = max(0.5,indiv.lai_flowering) * a1m;
		double a2 = indiv.pft.lifeform == GRASS_ANNUAL ? 1.0 : 0.25;
		//double fac2 = min(1.0, max(0.0, 1.0 - (a1 * 0.4 * frepr_cstruct - 1.3 * indiv.pft.reprfrac)));
		//indiv.fdev_growth *= min(1.0, max(0.0, 1.0 - (a1 * frepr_cstruct - a2 * indiv.pft.reprfrac)));
		indiv.fdev_growth *= min(1.0, max(0.0, 1.0 - a1 * (frepr_cstruct / indiv.pft.reprfrac - a2)));


	} else if (indiv.growth_stage == vegetative && indiv.fdev_growth < 1.0) {

		// Phase after switching from reproductive period back to vegetative period.
        // check if conditions are suitable to resume vegetative growth.
		// If not, maintain reduced growth rate or enter dormancy.

		// TODO: deal with dphen_activity called twice.
		double dphen_growth = dphen_activity(indiv, climate, false);

		// trigger dormancy at higher soil moisture thresholds just after reproductive phase is finished.
		if (ppft.wscal_5.mean() < 0.6) { 
			indiv.fphen_activity = min(1.0,max(fphen_activity_min,indiv.fphen_activity * (1.0 + dphen_growth)));
		}
		indiv.fdev_growth = min(1.0,max(fdev_growth_min * indiv.fphen_activity,indiv.fdev_growth * (1.0 + dphen_growth)));

	}
}


void phenology_grass_dave(Vegetation& vegetation, Climate& climate) {

	// Calculates phenological cycle of grasses. 
	// Returns phenological activity factor fphen_activity [0-1].
	// Implemented for grasses only. 
    for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		Pft &pft = pftlist.getobj();
		Patchpft &ppft = vegetation.patch.pft[pft.id];
		ppft.soiltemp_10.add(vegetation.patch.soil.T_soil[7]);   // temp of upper soil layer
		ppft.daylength_10.add(climate.daylength);
		ppft.reldaylength = climate.daylength / climate.daylength_max;
	}


	// Loop through cohorts (=individuals)
	vegetation.firstobj();
	while (vegetation.isobj) {

		// Retrieve information about cohort
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft &ppft = vegetation.patch.pft[indiv.pft.id];

		if (pft.is_grass()) { 

			// 1) Dormancy
			if (indiv.pft.lifeform == GRASS_PERENNIAL) {

				double dphen_dormancy = dphen_activity(indiv, climate, true);
				indiv.fphen_activity = min(1.0,max(fphen_activity_min,indiv.fphen_activity * (1.0 + dphen_dormancy)));
			} else {
			    // no dormancy in annuals
				indiv.fphen_activity = 1.0;
			}


			// 2) Development factor for growth (returns indiv.fdev_growth, which is used in
			//    growth_demand() function returning growth demand).

			// Peak aboveground structural biomass and ratio of reproductive biomass to peak aboveground structural biomass
			//indiv.cmass_struct = indiv.get_cmass_leaf() > indiv.cmass_struct ? indiv.get_cmass_leaf() : indiv.cmass_struct;
			indiv.cmass_struct = max(1.0e-9,indiv.get_cmass_leaf());
			indiv.frepr_cstruct = min(1.0e3,indiv.cmass_repr / indiv.cmass_struct);

			dev_growth(indiv, climate, indiv.frepr_cstruct);
			

			// 3) Switches between growth stages (indiv.growth_stage):
			// 		0: dormant stage
			// 		1: vegetative stage
			// 		2: reproductive stage
			double repr_trend = 0.0;  // relative trend in repr. tissue
			double cmass_storage_level = indiv.cstorage_max > 0.0 ? indiv.cmass_storage / indiv.cstorage_max : 0.0; 

			// TODO: get rid of growth_stage_days when they aren't needed anymore.
			// TODO: maybe increase dormancy switch threshold to 0.5 or similar?
		   	if (indiv.growth_stage == dormant && indiv.fphen_activity > fphen_activity_min) {
				indiv.growth_stage = vegetative;
				indiv.fdev_growth = 1.0;
				indiv.growth_stage_days = 0;
				indiv.cmass_struct = 0.0;
			}
			else if (indiv.growth_stage == vegetative) {
				if (indiv.cmass_storage_20.size() > 5 && indiv.cmass_storage_20.periodictrend(15) > 0.0 &&
				    indiv.growth_stage_days > 15 && indiv.fdev_growth > 0.5 && indiv.lai > 0.0){
					indiv.growth_stage = reproductive;
					indiv.growth_stage_days = 0;
					indiv.lai_flowering = indiv.lai;
				}
				if (equal(indiv.fphen_activity, fphen_activity_min)) {		
					indiv.growth_stage = dormant;
					indiv.growth_stage_days = 0;
				}
			} 
			else if (indiv.growth_stage == reproductive && indiv.pft.lifeform == GRASS_PERENNIAL) {
				// stop reproductive phase either 
				// a) at the onset of dormancy
                // b) repr pool does not increase any longer
                repr_trend = indiv.cmass_repr_10.size() > 9 ? indiv.cmass_repr_10.periodictrend(10)/indiv.cmass_repr_10.mean() : 1.0;

                // dormancy
                if (equal(indiv.fphen_activity, fphen_activity_min)) {
                    indiv.growth_stage = dormant;
                    indiv.growth_stage_days = 0;
                }

                // no further increases in repr. 
                // Note: does not go in there during freenyears. Is that a problem?
                if (repr_trend <= 0.001) {
                   indiv.growth_stage = vegetative;
                   indiv.growth_stage_days = 0;
                   indiv.fdev_growth = fdev_growth_min;
				   indiv.cmass_struct = 0.0;
                }
 
                // flush reproductive tissue to seedbank.
                if (indiv.growth_stage != reproductive && indiv.cmass_repr > 0.0) {
                    indiv.flush_repr_litter_DAVE();
				}
			}
			indiv.growth_stage_days += 1;
		} 
		else {
			indiv.fphen_activity = 1.0;
		}

		vegetation.nextobj();
	}
}

void leaf_phenology_pft(Pft& pft, Climate& climate, double wscal, double aphen,
	double& phen) {

	// DESCRIPTION
	// Calculates leaf phenological status (fractional leaf-out) for a individuals of
	// a given PFT, given current heat sum and length of chilling period (summergreen
	// PFTs) and water stress coefficient (raingreen PFTs)

	// INPUT PARAMETER
	// wscal = water stress coefficient (0-1; 1=minimum stress)
	// aphen = sum of daily fractional leaf cover (equivalent number of days with
	//         full leaf cover) so far this growing season

	// OUTPUT PARAMETER
	// phen = fraction of full leaf cover for any individual of this PFT

	// Assaf: changed for as imple solution for phen for evergreens.
	
	double aphen_max = APHEN_MAX;

	bool raingreen = pft.phenology == RAINGREEN || pft.phenology == ANY;
	bool summergreen = pft.phenology == SUMMERGREEN || pft.phenology == ANY;

	phen = 1.0;
	
	// TODO: Cuurently copied as-is from the Grote model. Consider changing 
	// evergreen phenology altogether.
	if (pft.phenology == EVERGREEN) {
		phen = 1.0;
		
		// handle others here as well e.g. evergreen (=needleleaf?) if it has a phengdd5ramp
		// not tropical with year round leaf senescence (SEN_BACKGROUND)
		// CHECK: this can be enabled for annual allocation, but it messes up the C:N in the early offseason in the year
		//        on 1st Jan, N seems to be allocated based on cmass_leaf (max) in the year, not based on phen*cmass_leaf (might just not be written out right)
		if (pft.lifeform == TREE && pft.phengdd5ramp > 0.0) {
			
			if (pft.gdd0[climate.chilldays] < 0.0) {
				pft.gdd0[climate.chilldays] = pft.k_chilla +
					pft.k_chillb * exp(-pft.k_chillk * (double)climate.chilldays);
			}
			// either based on aphen<(210+30) or aphen_leafout<210
			if (climate.gdd5 > pft.gdd0[climate.chilldays] && (aphen < aphen_max)) {
				phen = min(1.0, (1.0-pft.turnover_leaf) + pft.turnover_leaf * 
				     (climate.gdd5 - pft.gdd0[climate.chilldays]) / pft.phengdd5ramp);
			} else {
				// evergreen start phen with (1-turnover_leaf), e.g. BNE 1-0.33=0.67
				phen = (1.0-pft.turnover_leaf);
				
			}
		}
	}

	if (summergreen) {

		// Summergreen PFT - phenology based on GDD5 sum

		if (pft.lifeform == TREE) {

			// Calculate GDD base value for this PFT (if not already known) given
			// current length of chilling period (Sykes et al 1996, Eqn 1)

			if (pft.gdd0[climate.chilldays] < 0.0)
				pft.gdd0[climate.chilldays] = pft.k_chilla +
					pft.k_chillb * exp(-pft.k_chillk * (double)climate.chilldays);

			if (climate.gdd5 > pft.gdd0[climate.chilldays] && aphen < APHEN_MAX)
				phen = min(1.0,
					(climate.gdd5 - pft.gdd0[climate.chilldays]) / pft.phengdd5ramp);
			else
				phen = 0.0;

		}
		else if (pft.is_grass() || pft.lifeform == MOSS) {

			// Summergreen grasses have no maximum number of leaf-on days per
			// growing season, and no chilling requirement

			phen = min(1.0, climate.gdd5 / pft.phengdd5ramp);
		}
	}

	if (raingreen && wscal < pft.wscal_min) {

		// Raingreen phenology based on water stress threshold
		phen = 0.0;
	}
}

// decides on phenological status of the tree individual for any given phen value. 
// Follows the graph in the grote document (except for evergreen)
// input: phen, winter_to_summer, leaf turnover
// output: indiv.leaf_out, indiv.leaf_senesced, indiv.leaf_full_switch_summer ,
// indiv.leaf_growth
void tree_phenology_grote(Individual& indiv) {
	// taken from Grote code.
	Climate climate = indiv.vegetation.patch.get_climate();

	// handle daily phen leaf state
	bool leaf_out = false;
	bool leaf_senesced = false;
	bool leaf_full_on = false;
	bool leaf_full_switch_summer = false;
	bool leaf_growth = false;
	switch (indiv.pft.phenology) {
		// they are all quite similar, but need to handle RAINGREEN better
		// FIXME: are tropical evergreens with a general background turnover handled in here correctly?
		case SUMMERGREEN:
			if (indiv.phen > 0.0 && indiv.phen < indiv.leaf_out_phen_maxlimit) {
				leaf_out = true;
			}
			leaf_senesced = negligible(indiv.phen) && !indiv.senescence_active;
			leaf_full_on = negligible(indiv.phen - 1.0);
			leaf_full_switch_summer = leaf_full_on && !climate.winter_to_summer;
			leaf_growth = !negligible(indiv.phen) && !leaf_full_on;
			break;
		case EVERGREEN:
			// true in case phen is between its minimum value (1 - turnover_leaf) and 0.9, else, false
			if ((indiv.phen  - (1.0 - indiv.pft.turnover_leaf)) > 0.0 && indiv.phen < indiv.leaf_out_phen_maxlimit) {
				leaf_out = true;
			}
			// true in case phen is at its minimum value and senescence is active, 
			// else false
			leaf_senesced = negligible(indiv.phen - (1.0 - indiv.pft.turnover_leaf)) && !indiv.senescence_active;
			// true in case phen is at its maximum value (=1), else false
			leaf_full_on = negligible(indiv.phen - 1.0);
			// true in case phen = 1 and it's after the summer solstice (days 
			// shortening), else false
			leaf_full_switch_summer = leaf_full_on && !climate.winter_to_summer;
			// true in case phen is not at its minimum value and not when it is 
			// at its maximum value, else false
			leaf_growth = !negligible(indiv.phen  - (1.0 - indiv.pft.turnover_leaf)) && !leaf_full_on;
			break;
		case RAINGREEN:
			// FIXME: raingreen is tricky, it can have multiple drought events and
			//        we should handle RAINGREEN differently: having the drop of leaves based on drought strength, not drop all leaves
			if ((indiv.phen  - (1.0-indiv.pft.turnover_leaf)) > 0.0 && indiv.phen < indiv.leaf_out_phen_maxlimit) {
				leaf_out = true;
			}
			leaf_senesced = false; // it is a whole year growth PFT // negligible(indiv.phen) && !indiv.senescence.active;
			leaf_full_on = negligible(indiv.phen - 1.0);
			leaf_full_switch_summer = leaf_full_on && !climate.winter_to_summer;
			leaf_growth = true; // it is a whole year growth PFT
			break;
		default:
			fail("%s: missing implementation of phen leaf state handling for phenology type %d\n", indiv.pft.phenology);
	}
	
	indiv.leaf_out = leaf_out;
	indiv.leaf_senesced = leaf_senesced;
	indiv.leaf_full_switch_summer = leaf_full_switch_summer;
	indiv.leaf_growth = leaf_growth;
#if VERBOSE > 1
	if (date.get_calendar_year() >= 1901)
		dprintf("%d leave: lo %d lg %d lfsw %d sd %d sen %d\n",
						date.day,
						(int) indiv.leaf_out,
						(int) indiv.leaf_growth,
						(int) indiv.leaf_full_switch_summer,
						(int) indiv.leaf_senesced,
						(int) indiv.senescence.active);
#endif

}

void leaf_phenology(Patch& patch, Climate& climate) {

	// DESCRIPTION
	// Updates leaf phenological status (fractional leaf-out) for Patch PFT objects and
	// all individuals in a particular patch.

	// Updated by Ben Smith 2002-07-24 for compatability with "fast" canopy exchange
	// code (phenology assigned to patchpft for all vegetation modes)

	// Obtain reference to Vegetation object
	Vegetation& vegetation = patch.vegetation;

	// INDIVIDUAL AND COHORT MODES
	// Calculate phenology for each PFT at this patch

	if(ifdavegrowth) {
		
		phenology_grass_dave(vegetation, climate);

	}

	// Loop through patch-PFTs
		
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& ppft = patch.pft.getobj();

		// For this PFT ...
		if(patch.stand.pft[ppft.id].active) {
			if(ppft.pft.landcover == CROPLAND && patch.stand.landcover == CROPLAND)
				leaf_phenology_crop(ppft.pft, patch);
			else	//natural, urban, pasture, forest and peatland stands/pft:s
				leaf_phenology_pft(ppft.pft, climate, ppft.wscal, ppft.aphen, ppft.phen);

			// Update monthly leaf-on sum
			ppft.mphen[date.month] += ppft.phen / date.ndaymonth[date.month];

			// Update annual leaf-on sum
			if ( (climate.lat >= 0.0 && date.day == COLDEST_DAY_NHEMISPHERE) ||
				 (climate.lat < 0.0 && date.day == COLDEST_DAY_SHEMISPHERE) ) {
				ppft.aphen = 0.0;
			}
			ppft.aphen += ppft.phen;

			// ... on to next PFT
		}
		patch.pft.nextobj();
	}


	// Copy PFT-specific phenological status to individuals of each PFT

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// For this individual ...
		if (!(indiv.pft.is_grass() && ifdavegrowth && ifdailygrass)) {
			indiv.phen = patch.pft[indiv.pft.id].phen;
		} else {
			// Dave daily grass - phenology accounted for elsewhere.
			indiv.phen = 1.0;
		}

		// Update annual leaf-day sum (raingreen PFTs)
		if (date.day == 0) indiv.aphen_raingreen = 0;
		indiv.aphen_raingreen += (indiv.phen != 0.0);

		if (indiv.pft.lifeform == TREE) {
			// decide on phenological state (from Grote)
			tree_phenology_grote(indiv);
		}
		// ... on to next individual
		vegetation.nextobj();
	}
}

/// Calculates nitrogen retranslocation fraction
/* Calculates actual nitrogen retranslocation fraction so maximum
 * nitrogen storage capacity is not exceeded
 */
double calc_nrelocfrac(lifeformtype lifeform, double turnover_leaf, double nmass_leaf,
	double turnover_root, double nmass_root, double turnover_sap, double nmass_sap, double max_n_storage, double longterm_nstore) {

	double turnover_nmass = turnover_leaf * nmass_leaf + turnover_root * nmass_root;

	if (lifeform == TREE)
		turnover_nmass += turnover_sap * nmass_sap;

	if (max_n_storage < longterm_nstore)
		return 0.0;
	else if (max_n_storage < longterm_nstore + turnover_nmass * nrelocfrac && !negligible(turnover_nmass))
		return (max_n_storage - longterm_nstore) / (turnover_nmass);
	else
		return nrelocfrac;
}

/// Calculates nitrogen retranslocation fraction for daily grass
/* Calculates actual nitrogen retranslocation fraction so maximum
* nitrogen storage capacity is not exceeded too much for daily grass
* At 70% of maximum storage capacity the fraction relocated is dropped
* At 100% it is zero
*/
double calc_nrelocfrac_daily(double cton_min, double cton_max, double cmass_storage, double nmass_storage) {

	// 20% decrease or increase in N retranslocation depending on if you are at 
	// the maximum storage capacity or minimum
	double nreloc_max = nrelocfrac * 1.2;
	double nreloc_min = nrelocfrac / 1.2;

	if (is_zero(cmass_storage)) {
		return nreloc_min;
	}

	// Current N storage level
	double ntoc_current = nmass_storage / cmass_storage;

	return max(nreloc_min, min(nreloc_max, nreloc_max - (nreloc_max - nreloc_min) *	
		(ntoc_current - 1.0 / cton_max) / (1.0 / cton_min - 1.0 / cton_max)));
}

///////////////////////////////////////////////////////////////////////////////////////
// TURNOVER
// Internal function (do not call directly from framework)

void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
	lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_sap,
	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
	double& nmass_heart, double& litter_leaf, double& litter_root, double& cmass_leaf_root_turnover,
	double& nmass_litter_leaf, double& nmass_litter_root,
	double& longterm_nstore, double &max_n_storage,
	bool alive) {

	// DESCRIPTION
	// Transfers carbon from leaves and roots to litter, and from sapwood to heartwood
	// Only turnover from 'alive' individuals is transferred to litter (Ben 2007-11-28)

	// INPUT PARAMETERS
	// turnover_leaf = leaf turnover per time period as a proportion of leaf C biomass
	// turnover_root = root turnover per time period as a proportion of root C biomass
	// turnover_sap  = sapwood turnover to heartwood per time period as a proportion of
	//                 sapwood C biomass
	// lifeform      = PFT life form class (TREE or GRASS)
	// alive         = signifies new Individual object if false (see vegdynam.cpp)

	// INPUT AND OUTPUT PARAMETERS
	// cmass_leaf    = leaf C biomass (kgC/m2)
	// cmass_root    = fine root C biomass (kgC/m2)
	// cmass_sap     = sapwood C biomass (kgC/m2)
	// nmass_leaf    = leaf nitrogen biomass (kgN/m2)
	// nmass_root    = fine root nitrogen biomass (kgN/m2)
	// nmass_sap     = sapwood nitrogen biomass (kgN/m2)

	// OUTPUT PARAMETERS
	// litter_leaf			= new leaf C litter (kgC/m2)
	// litter_root			= new root C litter (kgC/m2)
	// nmass_litter_leaf	= new leaf nitrogen litter (kgN/m2)
	// nmass_litter_root	= new root nitrogen litter (kgN/m2)
    // cmass_leaf_root_turnover = litter produced in turnover (kgC/m2)
	// cmass_heart			= heartwood C biomass (kgC/m2)
	// nmass_heart			= heartwood nitrogen biomass (kgC/m2)
	// longterm_nstore		= longterm nitrogen storage (kgN/m2)

	double turnover = 0.0;
	// Calculate actual nitrogen retranslocation so maximum nitrogen storage capacity is not exceeded
	double actual_nrelocfrac = calc_nrelocfrac(lifeform, turnover_leaf, nmass_leaf, turnover_root, nmass_root,
	                                           turnover_sap, nmass_sap, max_n_storage, longterm_nstore);

	// TREES AND GRASSES:

	// Leaf turnover
	turnover = turnover_leaf * cmass_leaf;
	cmass_leaf -= turnover;
	if (alive) litter_leaf += turnover;
	if (alive) cmass_leaf_root_turnover += turnover;

	turnover = turnover_leaf * nmass_leaf;
	nmass_leaf -= turnover;
	assert(nmass_leaf >= 0);
	nmass_litter_leaf += turnover * (1.0 - actual_nrelocfrac);
	longterm_nstore += turnover * actual_nrelocfrac;

	// Root turnover
	turnover = turnover_root * cmass_root;
	cmass_root -= turnover;
	if (alive) litter_root += turnover;
	if (alive) cmass_leaf_root_turnover += turnover;

	turnover = turnover_root * nmass_root;
	nmass_root -= turnover;
	nmass_litter_root += turnover * (1.0 - actual_nrelocfrac);
	longterm_nstore += turnover * actual_nrelocfrac;

	if (lifeform == TREE) {

		// TREES ONLY:

		// Sapwood turnover by conversion to heartwood
		turnover = turnover_sap * cmass_sap;
		cmass_sap -= turnover;
		cmass_heart += turnover;

		// NB: assumes nitrogen is translocated from sapwood prior to conversion to
		//     heartwood and that this is the same fraction that is conserved
		//     in conjunction with leaf and root shedding

		turnover = turnover_sap * nmass_sap;
		nmass_sap -= turnover;
		nmass_heart += turnover * (1.0 - actual_nrelocfrac);
		longterm_nstore += turnover * actual_nrelocfrac;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// REPRODUCTION
// Internal function (do not call directly from framework)

void reproduction(double reprfrac, double npp, double& bminc, double& cmass_repr) {

	// DESCRIPTION
	// Allocation of net primary production (NPP) to reproduction and calculation of
	// assimilated carbon available for production of new biomass

	// INPUT PARAMETERS
	// reprfrac = fraction of NPP for this time period allocated to reproduction
	// npp      = NPP (i.e. assimilation minus maintenance and growth respiration) for
	//            this time period (kgC/m2)

	// OUTPUT PARAMETER
	// bminc    = carbon biomass increment (component of NPP available for production
	//            of new biomass) for this time period (kgC/m2)

	if (npp >= 0.0) {
		cmass_repr = npp * reprfrac;
		bminc = npp - cmass_repr;
		return;
	}

	// Negative NPP - no reproduction cost

	cmass_repr = 0.0;
	bminc = npp;
}


///////////////////////////////////////////////////////////////////////////////////////
// ALLOCATION
// Function allocation is an internal function (do not call directly from framework);
// function allocation_init may be called to distribute initial biomass among tissues
// for a new individual.

// File scope global variables: used by function f below (see function allocation)

static double k1, k2, k3, b;
static double ltor_g;
static double cmass_heart_g;
static double cmass_leaf_g;

inline double f(double& cmass_leaf_inc) {

	// Returns value of f(cmass_leaf_inc), given by:
	//
	// f(cmass_leaf_inc) = 0 =
	//   k1 * (b - cmass_leaf_inc - cmass_leaf_inc/ltor + cmass_heart) -
	//   [ (b - cmass_leaf_inc - cmass_leaf_inc/ltor)
	//   / (cmass_leaf + cmass_leaf_inc )*k3 ] ** k2
	//
	// See function allocation (below), Eqn (13)

	return k1 * (b - cmass_leaf_inc - cmass_leaf_inc / ltor_g + cmass_heart_g) -
		pow((b - cmass_leaf_inc - cmass_leaf_inc / ltor_g) / (cmass_leaf_g + cmass_leaf_inc) * k3,
		k2);
}


void allocation(double bminc, Individual& indiv, double ltor, double& cmass_leaf_inc, double& cmass_root_inc,
		double& cmass_sap_inc, double& cmass_debt_inc, double& cmass_heart_inc,
		double& litter_leaf_inc, double& litter_root_inc, double& exceeds_cmass) {

	// leaf C biomass for last time period on individual basis (kgC)
	double cmass_leaf = indiv.get_cmass_leaf() / indiv.densindiv;
	// root C biomass for last time period on individual basis (kgC)
	double cmass_root = indiv.cmass_root / indiv.densindiv;
	// sapwood C biomass for last time period on individual basis (kgC)
	double cmass_sap = indiv.cmass_sap / indiv.densindiv;
	// C debt for last time period on individual basis (kgC)
	double cmass_debt = indiv.cmass_debt / indiv.densindiv;
	// heartwood C biomass for last time period on individual basis (kgC)
	double cmass_heart = indiv.cmass_heart / indiv.densindiv;
	// individual height (m)
	double height = indiv.height;
	// specific leaf area (PFT-specific constant) (m2/kgC)
	double sla = indiv.sla();
	// wood density (PFT-specific constant) (kgC/m3)
	double wooddens = indiv.pft.wooddens;
	// life form class (TREE or GRASS)
	lifeformtype lifeform = indiv.pft.lifeform;
	// ratio of leaf area to sapwood cross-sectional area (PFT-specific constant)          
	double k_latosa = indiv.pft.k_latosa;
	// constant in allometry equations
	double k_allom2 = indiv.pft.k_allom2;
	// constant in allometry equations
	double k_allom3 = indiv.pft.k_allom3;

	// DESCRIPTION
	// Calculates changes in C compartment sizes (leaves, roots, sapwood, heartwood)
	// and litter for a plant individual as a result of allocation of biomass increment.
	// Assumed allometric relationships are given in function allometry below.

	// INPUT PARAMETERS
	// bminc       = biomass increment this time period on individual basis (kgC)
	// ltor        = leaf to root mass ratio following allocation

	// OUTPUT PARAMETERS
	// cmass_leaf_inc  = increment (may be negative) in leaf C biomass following
	//                   allocation (kgC)
	// cmass_root_inc  = increment (may be negative) in root C biomass following
	//                   allocation (kgC)
	// cmass_sap_inc   = increment (may be negative) in sapwood C biomass following
	//                   allocation (kgC)
	// cmass_heart_inc = increment in heartwood C biomass following allocation (kgC)
	// litter_leaf_inc = increment in leaf litter following allocation, on individual
	//                   basis (kgC)
	// litter_root_inc = increment in root litter following allocation, on individual
	//                   basis (kgC)
	// exceeds_cmass      = negative increment that exceeds existing biomass (kgC)

	// MATHEMATICAL DERIVATION FOR TREE ALLOCATION
	// Allocation attempts to distribute biomass increment (bminc) among the living
	// tissue compartments, i.e.
	//   (1) bminc = cmass_leaf_inc + cmass_root_inc + cmass_sap_inc
	// while satisfying the allometric relationships (Shinozaki et al. 1964a,b; Waring
	// et al 1982, Huang et al 1992; see also function allometry, below) [** =
	// raised to the power of]:
	//   (2) (leaf area) = k_latosa * (sapwood xs area)
	//   (3) cmass_leaf = ltor * cmass_root
	//   (4) height = k_allom2 * (stem diameter) ** k_allom3
	// From (1) and (3),
	//   (5) cmass_sap_inc = bminc - cmass_leaf_inc -
	//         (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root
	// Let diam_new and height_new be stem diameter and height following allocation.
	// Then (see allometry),
	//   (6) diam_new = 2 * [ ( cmass_sap + cmass_sap_inc + cmass_heart )
	//         / wooddens / height_new / PI ]**(1/2)
	// From (4), (6) and (5),
	//   (7) height_new**(1+2/k_allom3) =
	//         k_allom2**(2/k_allom3) * 4 * [cmass_sap + bminc - cmass_leaf_inc
	//         - (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root + cmass_heart]
	//         / wooddens / PI
	// Now,
	//   (8) wooddens = cmass_sap / height / (sapwood xs area)
	// From (8) and (2),
	//   (9) wooddens = cmass_sap / height / sla / cmass_leaf * k_latosa
	// From (9) and (1),
	//  (10) wooddens = (cmass_sap + bminc - cmass_leaf_inc -
	//         (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root)
	//          / height_new / sla / (cmass_leaf + cmass_leaf_inc) * k_latosa
	// From (10),
	//  (11) height_new**(1+2/k_allom3) =
	//         [ (cmass_sap + bminc - cmass_leaf_inc - (cmass_leaf + cmass_leaf_inc)
	//           / ltor + cmass_root) / wooddens / sla
	//           / (cmass_leaf + cmass_leaf_inc ) * k_latosa ] ** (1+2/k_allom3)
	//
	// Combining (7) and (11) gives a function of the unknown cmass_leaf_inc:
	//
	//  (12) f(cmass_leaf_inc) = 0 =
	//         k_allom2**(2/k_allom3) * 4/PI * [cmass_sap + bminc - cmass_leaf_inc
	//         - (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root + cmass_heart]
	//         / wooddens -
	//         [ (cmass_sap + bminc - cmass_leaf_inc - (cmass_leaf + cmass_leaf_inc)
	//           / ltor + cmass_root) / (cmass_leaf + cmass_leaf_inc)
	//           / wooddens / sla * k_latosa] ** (1+2/k_allom3)
	//
	// Let k1 = k_allom2**(2/k_allom3) * 4/PI / wooddens
	//     k2 = 1+2/k_allom3
	//     k3 = k_latosa / wooddens / sla
	//     b  = cmass_sap + bminc - cmass_leaf/ltor + cmass_root
	//
	// Then,
	//  (13) f(cmass_leaf_inc) = 0 =
	//         k1 * (b - cmass_leaf_inc - cmass_leaf_inc/ltor + cmass_heart) -
	//         [ (b - cmass_leaf_inc - cmass_leaf_inc/ltor)
	//         / (cmass_leaf + cmass_leaf_inc )*k3 ] ** k2
	//
	// Numerical methods are used to solve Eqn (13) for cmass_leaf_inc

	const int NSEG=20; // number of segments (parameter in numerical methods)
	const int JMAX=40; // maximum number of iterations (in numerical methods)
	const double XACC=0.0001; // threshold x-axis precision of allocation solution
	const double YACC=1.0e-10; // threshold y-axis precision of allocation solution
	const double CDEBT_MAXLOAN_DEFICIT=0.8; // maximum loan as a fraction of deficit
	const double CDEBT_MAXLOAN_MASS=0.2; // maximum loan as a fraction of (sapwood-cdebt)

	double cmass_leaf_inc_min;
	double cmass_root_inc_min;
	double x1,x2,dx,xmid,fx1,fmid,rtbis,sign;
	int j;
	double cmass_deficit,cmass_loan;

	// initialise
	litter_leaf_inc = 0.0;
	litter_root_inc = 0.0;
	exceeds_cmass   = 0.0;
	cmass_leaf_inc  = 0.0;
	cmass_root_inc  = 0.0;
	cmass_sap_inc   = 0.0;
	cmass_heart_inc = 0.0;
	cmass_debt_inc = 0.0;

	if (!largerthanzero(ltor, -10)) {

		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_leaf_inc = 0.0;

		// Make sure we don't end up with negative cmass_root
		if (bminc < -cmass_root) {
			exceeds_cmass = -(cmass_root + bminc);
			cmass_root_inc = -cmass_root;
		}
		else {
			cmass_root_inc=bminc;
		}

		if (lifeform==TREE) {
			cmass_sap_inc=-cmass_sap;
			cmass_heart_inc=-cmass_sap_inc;
		}
	}
	else if (lifeform==TREE) {

		// TREE ALLOCATION

		cmass_heart_inc=0.0;

		// Calculate minimum leaf increment to maintain current sapwood biomass
		// Given Eqn (2)

		if (height>0.0)
			cmass_leaf_inc_min=k_latosa*cmass_sap/(wooddens*height*sla)-cmass_leaf;
		else
			cmass_leaf_inc_min=0.0;

		// Calculate minimum root increment to support minimum resulting leaf biomass
		// Eqn (3)

		if (height>0.0)
			cmass_root_inc_min=k_latosa*cmass_sap/(wooddens*height*sla*ltor)-
				cmass_root;
		else
			cmass_root_inc_min=0.0;

		if (cmass_root_inc_min<0.0) { // some roots would have to be killed

			cmass_leaf_inc_min=cmass_root*ltor-cmass_leaf;
			cmass_root_inc_min=0.0;
		}

		// BLARP! C debt stuff
		if (ifcdebt) {
			cmass_deficit=cmass_leaf_inc_min+cmass_root_inc_min-bminc;
			if (cmass_deficit>0.0) {
				cmass_loan=max(min(cmass_deficit*CDEBT_MAXLOAN_DEFICIT,
					(cmass_sap-cmass_debt)*CDEBT_MAXLOAN_MASS),0.0);
				bminc+=cmass_loan;
				cmass_debt_inc=cmass_loan;
			}
			else cmass_debt_inc=0.0;
		}
		else cmass_debt_inc=0.0;

		if ( (cmass_root_inc_min >= 0.0 && cmass_leaf_inc_min >= 0.0 &&
		      cmass_root_inc_min + cmass_leaf_inc_min <= bminc) || bminc<=0.0) {

			// Normal allocation (positive increment to all living C compartments)

			// Calculation of leaf mass increment (lminc_ind) satisfying Eqn (13)
			// using bisection method (Press et al 1986)

			// Set values for global variables for reuse by function f

			k1 = pow(k_allom2, 2.0 / k_allom3) * 4.0 / PI / wooddens;
			k2 = 1.0 + 2 / k_allom3;
			k3 = k_latosa / wooddens / sla;
			b = cmass_sap + bminc - cmass_leaf / ltor + cmass_root;
			ltor_g = ltor;
			cmass_leaf_g = cmass_leaf;
			cmass_heart_g = cmass_heart;

			x1 = 0.0;
			x2 = (bminc - (cmass_leaf / ltor - cmass_root)) / (1.0 + 1.0 / ltor);
			dx = (x2 - x1) / (double)NSEG;

			if (cmass_leaf < 1.0e-10) x1 += dx; // to avoid division by zero

			// Evaluate f(x1), i.e. Eqn (13) at cmass_leaf_inc = x1

			fx1 = f(x1);

			// Find approximate location of leftmost root on the interval
			// (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
			// change in sign of f(xmid) relative to f(x1).

			fmid = f(x1);

			xmid = x1;

			while (fmid * fx1 > 0.0 && xmid < x2) {

				xmid += dx;
				fmid = f(xmid);
			}

			x1 = xmid - dx;
			x2 = xmid;

			// Apply bisection to find root on new interval (x1,x2)

			if (f(x1) >= 0.0) sign = -1.0;
			else sign = 1.0;

			rtbis = x1;
			dx = x2 - x1;

			// Bisection loop
			// Search iterates on value of xmid until xmid lies within
			// xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

			fmid = 1.0; // dummy value to guarantee entry into loop
			j = 0; // number of iterations so far

			while (dx >= XACC && fabs(fmid) > YACC && j <= JMAX) {

				dx *= 0.5;
				xmid = rtbis + dx;

				fmid = f(xmid);

				if (fmid * sign <= 0.0) rtbis = xmid;
				j++;
			}

			// Now rtbis contains numerical solution for cmass_leaf_inc given Eqn (13)

			cmass_leaf_inc = rtbis;

			// Calculate increments in other compartments

			cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
			cmass_sap_inc = bminc - cmass_leaf_inc - cmass_root_inc; // Eqn (1)

			// guess2008 - extra check - abnormal allocation can still happen if ltor is very small
			if ((cmass_root_inc > 50 || cmass_root_inc < -50) && ltor < 0.0001) {
				cmass_leaf_inc = 0.0;
				cmass_root_inc = bminc;
				cmass_sap_inc = -cmass_sap;
				cmass_heart_inc = -cmass_sap_inc;
			}

			// Negative sapwood increment larger than existing sapwood or
			// if debt becomes larger than existing woody biomass
			if (cmass_sap < -cmass_sap_inc || cmass_sap + cmass_sap_inc + cmass_heart < cmass_debt + cmass_debt_inc) {

				// Abnormal allocation: reduction in some biomass compartment(s) to
				// satisfy allometry

				// Attempt to distribute this year's production among leaves and roots only
				// Eqn (3)

				cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
				cmass_root_inc = bminc - cmass_leaf_inc;

				// Make sure we don't end up with negative cmass_leaf
				cmass_leaf_inc = max(-cmass_leaf, cmass_leaf_inc);

				// Make sure we don't end up with negative cmass_root
				cmass_root_inc = max(-cmass_root, cmass_root_inc);

				// If biomass of roots and leaves can't meet biomass decrease then
				// sapwood also needs to decrease
				cmass_sap_inc = bminc - cmass_leaf_inc - cmass_root_inc;

				// No sapwood turned into heartwood
				cmass_heart_inc = 0.0;

				// Make sure we don't end up with negative cmass_sap
				if (cmass_sap_inc < -cmass_sap) {
					exceeds_cmass = -(cmass_sap + cmass_sap_inc);
					cmass_sap_inc = -cmass_sap;
				}

				// Comment: Can happen that biomass decrease is larger than biomass in all compartments.
				// Then bminc is more negative than there is biomass to respire
			}
		}
		else {

			// Abnormal allocation: reduction in some biomass compartment(s) to
			// satisfy allometry

			// Attempt to distribute this year's production among leaves and roots only
			// Eqn (3)

			cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);

			if (cmass_leaf_inc > 0.0) {

				// Positive allocation to leaves

				cmass_root_inc = bminc - cmass_leaf_inc; // Eqn (1)

				// Add killed roots (if any) to litter

				// guess2008 - back to LPJF method in this case
				// if (cmass_root_inc<0.0) litter_root_inc=-cmass_root_inc;
				if (cmass_root_inc < 0.0) {
					cmass_leaf_inc = bminc;
					cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
				}

			}
			else {

				// Negative or zero allocation to leaves
				// Eqns (1), (3)

				cmass_root_inc = bminc;
				cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf;
			}

			// Make sure we don't end up with negative cmass_leaf
			if (cmass_leaf_inc < -cmass_leaf) {
				exceeds_cmass += -(cmass_leaf + cmass_leaf_inc);
				cmass_leaf_inc = -cmass_leaf;
			}

			// Make sure we don't end up with negative cmass_root
			if (cmass_root_inc < -cmass_root) {
				exceeds_cmass += -(cmass_root + cmass_root_inc);
				cmass_root_inc = -cmass_root;
			}

			// Add killed leaves to litter
			litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

			// Add killed roots to litter
			litter_root_inc = max(-cmass_root_inc, 0.0);

			// Calculate increase in sapwood mass (which must be negative)
			// Eqn (2)
			cmass_sap_inc = (cmass_leaf_inc + cmass_leaf) * wooddens * height * sla / k_latosa -
				cmass_sap;

			// Make sure we don't end up with negative cmass_sap
			if (cmass_sap_inc < -cmass_sap) {
				exceeds_cmass += -(cmass_sap + cmass_sap_inc);
				cmass_sap_inc = -cmass_sap;
			}

			// Convert killed sapwood to heartwood
			cmass_heart_inc = -cmass_sap_inc;
		}
	}
	else {

		// GRASS and MOSS ALLOCATION
		// Allocation attempts to distribute biomass increment (bminc) among leaf
		// and root compartments, i.e.
		//   (14) bminc = cmass_leaf_inc + cmass_root_inc
		// while satisfying Eqn(3)

		cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
		cmass_root_inc = bminc - cmass_leaf_inc;

		if (bminc >= 0.0) {

			// Positive biomass increment

			if (cmass_leaf_inc < 0.0) {

				// Positive bminc, but ltor causes negative allocation to leaves,
				// put all of bminc into roots

				cmass_root_inc = bminc;
				cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf; // Eqn (3)
			}
			else if (cmass_root_inc < 0.0) {

				// Positive bminc, but ltor causes negative allocation to roots,
				// put all of bminc into leaves

				cmass_leaf_inc = bminc;
				cmass_root_inc = (cmass_leaf + bminc) / ltor - cmass_root;
			}

			// Make sure we don't end up with negative cmass_leaf
			if (cmass_leaf_inc < -cmass_leaf) {
				exceeds_cmass += -(cmass_leaf + cmass_leaf_inc);
				cmass_leaf_inc = -cmass_leaf;
			}

			// Make sure we don't end up with negative cmass_root
			if (cmass_root_inc < -cmass_root) {
				exceeds_cmass += -(cmass_root + cmass_root_inc);
				cmass_root_inc = -cmass_root;
			}

			// Add killed leaves to litter
			litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

			// Add killed roots to litter
			litter_root_inc = max(-cmass_root_inc, 0.0);
		}
		else {

			// Abnormal allocation: negative biomass increment

			// Negative increment is respiration (neg bminc) or/and increment in other
			// compartments leading to no litter production

			if (bminc < -(cmass_leaf + cmass_root)) {

				// Biomass decrease is larger than existing biomass

				exceeds_cmass = -(bminc + cmass_leaf + cmass_root);

				cmass_leaf_inc = -cmass_leaf;
				cmass_root_inc = -cmass_root;
			}
			else if (cmass_root_inc < 0.0) {

				// Negative allocation to root
				// Make sure we don't end up with negative cmass_root

				if (cmass_root < -cmass_root_inc) {
					cmass_leaf_inc = bminc + cmass_root;
					cmass_root_inc = -cmass_root;
				}
			}
			else if (cmass_leaf_inc < 0.0) {

				// Negative allocation to leaf
				// Make sure we don't end up with negative cmass_leaf

				if (cmass_leaf < -cmass_leaf_inc) {
					cmass_root_inc = bminc + cmass_leaf;
					cmass_leaf_inc = -cmass_leaf;
				}
			}
		}
	}

	// Check C budget after allocation

	// maximum carbon mismatch
	double EPS = 1.0e-12;

	assert(fabs(bminc + exceeds_cmass - (cmass_leaf_inc + cmass_root_inc + cmass_sap_inc + cmass_heart_inc + litter_leaf_inc + litter_root_inc)) < EPS);
}

// Capacity of C storage (kg C m-2)
double cstorage_capacity(Individual& indiv) {

	if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
		return cstorage_max_temp * (indiv.get_cmass_leaf() + indiv.cmass_root) + cstorage_max_perm * indiv.cmass_crown;
	} else {
		// TODO: what about seedlings? cstorage_max cn be much higher than the total cmass of all the pools
		const int nyear_turnover = 3;
		double turnover_sap = indiv.cmass_sap * indiv.pft.turnover_sap;
		double turnover_leaf = indiv.cmass_leaf * indiv.pft.turnover_leaf;
		double turnover_root = indiv.cmass_root * indiv.pft.turnover_root;
		double cstorage_opt = nyear_turnover * (turnover_sap + turnover_leaf + turnover_root);
		return max(cstorage_opt, indiv.cstorage_max);
		// return cstorage_max_temp * (indiv.get_cmass_leaf() + indiv.cmass_root) + cstorage_max_perm * indiv.cmass_sap;
	}

}

// Capacity of N storage (kg N m-2)
double nstorage_capacity(Individual& indiv) {

	return indiv.cstorage_max / indiv.pft.cton_root_min;  // TODO: Rethink upper limit for storage.
														  //		--> prevent C:N ratio in storage from becoming too crazy through N storage limit. 

}

/*
Given a potential C and N growth/increment amounts for a particular pool,
constrain them such that their ratio doesn't exceed the bounds defined by
cton_min and cton_max. Any excess C or N which can't be allocated will be added
into the creturn and nreturn variables, and cgrow or ngrow will be reduced by
the same amount.

Notice that cgrow, ngrow, creturn, and nreturn are all passed by reference (ie
they are in/out parameters).

@param cgrow: The potential C growth for this pool.
@param ngrow: The potential N growth for this pool.
@param cton_min: Minimum C:N ratio of this pool.
@param cton_max: Maximum C:N ratio of this pool.
@param creturn: Any C in cgrow which can't be allocated to this pool (ie any C
which, if allocated, would cause the C:N ratio to exceed cton_max) will be added
into this variable. This must be initialised before calling this function.
@param nreturn: Any N in ngrow which can't be allocated to this pool (ie any N
which, if allocated, would cause the C:N ratio to fall below cton_min) will be
added into this variable. This must be initialised before calling this function.
*/
void constrain_growth(double &cgrow, double &ngrow, double cton_min, double cton_max, double &creturn, double &nreturn) {
	// When the upper/lower C:N limits are reached, only allow new growth to
	// move the C:N ratio to within this fraction of the limits. If we don't
	// have this limits, growth will move the C:N ratio exactly onto the limits,
	// and rounding errors can then push the C:N ratio beyond the limits.
	const double FGROWTH_MAX = 0.999;

	if ( ngrow < (cgrow / cton_max)) {
		// C:N ratio in new growth exceeds max C:N ratio.Limit C in new growth
		// to the max C:N ratio and put the rest into creturn.
		double cgrow_max = ngrow * cton_max * FGROWTH_MAX;
		creturn += cgrow - cgrow_max;
		cgrow = cgrow_max;
	} else if ( ngrow > (cgrow / cton_min) ) {
		// N:C ratio in new growth exceeds max N:C ratio. Limit N in new growth
		// to the max N:C ratio and put the rest into nreturn.
		double ngrow_max = FGROWTH_MAX * cgrow / cton_min;
		nreturn += ngrow - ngrow_max;
		ngrow = ngrow_max;
	}
}

void establishment_biomass(double cmass_init, double ltor, Individual& indiv, double nmass_init) {

	double dval = 0;
	double cmass_leaf_ind = 0;
	double cmass_root_ind = 0;
	double cmass_sap_ind = 0;

	double cmass_before = indiv.ccont(1.0, false, true);
	double nmass_before = indiv.ncont(1.0, false, true);
	
	if (indiv.pft.is_grass() && ifdailygrass) {

		if (ifdavephotosyn){
			// allocate initial cmass and nmass to plant pools

			double init_cton = cmass_init / nmass_init;

			// following equations solved for cshoot:
			// cmass_init = cshoot + croot + cstorage
			// croot = cshoot / ltor
			// cstorage = frac_cstorage_max * cmass_init
			double cmass_leaf = (cmass_init - cstorage_max_temp * cmass_init) * (ltor / (ltor + 1.0));
			double cmass_root = cmass_leaf / ltor;
			double cmass_storage = cmass_init - cmass_leaf - cmass_root;

			double nmass_leaf = cmass_leaf / init_cton;
			double nmass_root = cmass_root / init_cton;
			double nmass_storage = cmass_storage / init_cton;

			constrain_growth(
				cmass_leaf, nmass_leaf,
				indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
				cmass_storage, nmass_storage);

			constrain_growth(
				cmass_root, nmass_root,
				indiv.pft.cton_root_min, indiv.pft.cton_root_max,
				cmass_storage, nmass_storage);

			indiv.cmass_leafs[green] += cmass_leaf;
			indiv.cmass_root    += cmass_root;
			indiv.cmass_storage += cmass_storage;
			indiv.cmass_crown   += 0.0;

			assert(equal(cmass_init, cmass_leaf + cmass_root + cmass_storage));
			
			// Initial N mass
			indiv.nmass_leafs[green] += nmass_leaf;
			indiv.nmass_root += nmass_root;
			indiv.nmass_storage += nmass_storage;
			indiv.nmass_crown += 0.0;

			assert(equal(nmass_init, nmass_leaf + nmass_root + nmass_storage));

			// Capacity of C and N storage
			indiv.cstorage_max = cstorage_capacity(indiv);
			indiv.nstorage_max = nstorage_capacity(indiv);

			// Growth stage is vegetative after establishment.
			// note: don't set growth_stage here as it could change the status of an already established plant.
			// indiv.growth_stage = vegetative;

			// Need nonzero LAI for fapar calculation.
			indiv.alive = true;
			allometry(indiv);
		} else {
			indiv.cmass_storage += cmass_init;
		 	indiv.nmass_storage += nmass_init;
		}

		indiv.alive = true;
		indiv.initialn = indiv.get_nmass_leaf() + indiv.nmass_root + indiv.nmass_storage + indiv.nmass_crown;

		double cmass_after = indiv.ccont();
		double nmass_after = indiv.ncont();
		assert(equal(cmass_before + cmass_init, cmass_after));
		assert(equal(nmass_before + nmass_init, nmass_after));
	}
	else {
		// trees (and mosses?)
		if (ifdavegrowth) {
			indiv.alive = true;
			// proportion of the initial biomass that goes to storage
			const double fract_stor = 0.5;
			indiv.cmass_storage += cmass_init * fract_stor;
			double cmass_total = cmass_init * (1 - fract_stor);
			Pft &pft = indiv.pft;

			// Initial height of saplings in m.
			const double INIT_HEIGHT = 0.5;
			double height = max(INIT_HEIGHT, indiv.height);

			// constant
			double x = (1 / indiv.pft.ltor_max) + 
			((pft.sla * height * pft.wooddens) / pft.k_latosa) + 1;
			cmass_leaf_ind = cmass_total / x;
			cmass_root_ind = cmass_leaf_ind / pft.ltor_max;
			cmass_sap_ind = cmass_total - (cmass_leaf_ind + cmass_root_ind);

			if (is_zero(indiv.height)) {
				indiv.height = height;
				// height = k_allom2 * diam ** k_allom3
				// diam = (height / k2) ^ (1 / k3)
				indiv.diam = pow(height / pft.k_allom2, 1.0 / pft.k_allom3);
				indiv.crownarea = min(indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp)
			    , indiv.pft.crownarea_max);
			}
		} else {
			// TODO: make sure tree establishment uses the C:N ratio of the seeds
			allocation(cmass_init, indiv, ltor, cmass_leaf_ind, cmass_root_ind, 
			cmass_sap_ind, dval, dval, dval, dval, dval);
			cmass_leaf_ind *= indiv.densindiv;
			cmass_root_ind *= indiv.densindiv;
			cmass_sap_ind *= indiv.densindiv;
		}
		
		indiv.cmass_leaf += cmass_leaf_ind;
		indiv.cmass_root += cmass_root_ind;
		indiv.cmass_sap += cmass_sap_ind;
		indiv.lai = indiv.cmass_leaf * indiv.pft.sla;
		// Capacity of C and N storage on establishment
		if (ifdavegrowth && is_zero(indiv.cstorage_max)) {
			indiv.cstorage_max = cstorage_capacity(indiv);
			indiv.nstorage_max = nstorage_capacity(indiv);
		}
		// initialise limits
		if (is_zero(indiv.cmass_sap_limit)) {
			indiv.cmass_sap_limit = indiv.cmass_sap;
			indiv.cmass_leaf_limit = indiv.cmass_leaf;
			indiv.cmass_root_limit = indiv.cmass_root;
		}

		double nmass_remaining = nmass_init;

		if (nmass_init > 0) {
			if (ifdavephotosyn) {
				// we need initial values for nmass_leaf when calling ifdavephotosyn
				double leafndemand = indiv.cmass_leaf / indiv.pft.cton_leaf_min;
				double rootndemand = indiv.cmass_root / indiv.pft.cton_root_min;
				double sapndemand = indiv.cmass_sap / indiv.pft.cton_sap_min;
				if (leafndemand + rootndemand + sapndemand > nmass_init) { 
					leafndemand = indiv.cmass_leaf / indiv.pft.cton_leaf_avr;
					rootndemand = indiv.cmass_root / indiv.pft.cton_root_avr;
					sapndemand = indiv.cmass_sap / indiv.pft.cton_sap_avr;
				}
				if (leafndemand + rootndemand + sapndemand > nmass_init) {
					leafndemand = indiv.cmass_leaf / indiv.pft.cton_leaf_max;
					rootndemand = indiv.cmass_root / indiv.pft.cton_root_max;
					sapndemand = indiv.cmass_sap / indiv.pft.cton_sap_max;
				}

				// Remember, the destination pools may be nonzero already if
				// this "establishment" event is contributing biomass towards an
				// existing individual.

				double nalloc_leaf = min(nmass_remaining, leafndemand);
				indiv.nmass_leaf += nalloc_leaf;
				nmass_remaining -= nalloc_leaf;
				assert(indiv.nmass_leaf >= 0);

				double nalloc_root = min(nmass_remaining, rootndemand);
				indiv.nmass_root += nalloc_root;
				nmass_remaining -= nalloc_root;
				assert(indiv.nmass_root >= 0);

				if (indiv.pft.lifeform == TREE && nmass_remaining > 0) {
					double nalloc_sap = min(nmass_remaining, sapndemand);
					indiv.nmass_sap += nalloc_sap;
					nmass_remaining -= nalloc_sap;
					assert(indiv.nmass_sap >= 0);
				}
			}
			indiv.nmass_storage += nmass_remaining;
			assert(indiv.nmass_storage >= 0);
		}

		// We need to report the establishment C and N quantities as fluxes for
		// trees, in much the same way we did for grasses. However, trees with
		// annual allocation are not considered "alive" until 1 year after
		// establishment, so we can't report the establishment flux until then.
		// That is handled inside growth(), which reports the live C content as
		// the establishment C, but it can't use the live N content, because
		// that's based on a year's worth of uptake, which is already reported
		// separately. Therefore, we need to record the establishment N content
		// now, in order to report it as an establishment flux later.
		indiv.initialn = indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_sap + indiv.nmass_heart;

		if (indiv.pft.dave_allocation()) {
			indiv.initialn += indiv.nmass_storage;
		}
	}
}

void establishment_sapling_dave(Individual& indiv, double ltor) {

	indiv.alive = true;
	// initial height (m)
	const double INIT_HEIGHT = 1;
	// initial storage fraction of capacity
	const double INIT_STORAGE_FRACTION = 1;

	indiv.height = INIT_HEIGHT;

	indiv.diam = pow(indiv.height / indiv.pft.k_allom2, 1 / indiv.pft.k_allom3);
	
	// C for new sapling
	indiv.cmass_sap = indiv.height * indiv.diam * indiv.diam / 4 * PI * indiv.pft.wooddens;

	double fns = indiv.pft.k_latosa / (indiv.pft.sla * indiv.pft.wooddens * indiv.height);

	indiv.cmass_leaf = indiv.cmass_sap * fns;

	indiv.cmass_root = indiv.cmass_leaf / ltor;

	indiv.cstorage_max = cstorage_capacity(indiv);

	indiv.cmass_storage = INIT_STORAGE_FRACTION * indiv.cstorage_max;

	// N of new sapling
	
	indiv.nmass_sap = indiv.cmass_sap / indiv.pft.cton_sap_min;

	indiv.nmass_leaf = indiv.cmass_leaf / indiv.pft.cton_leaf_min;

	indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_min;

	indiv.nstorage_max = nstorage_capacity(indiv);

	indiv.nmass_storage = INIT_STORAGE_FRACTION * indiv.nstorage_max;

	indiv.crownarea = min(indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp), 
			indiv.pft.crownarea_max);
	// calculate allometry ( current state of the individual should be 0, including
	// cgrow_sap)
	allometry_dave_tree(indiv);

	// for nitrogen balance mass balance checks
	indiv.initialn = indiv.nmass_sap + indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_storage;

}

void allocation_init(double bminit, double ltor, Individual& indiv) {

	// DESCRIPTION
	// Allocates initial biomass among tissues for a new individual (tree or grass),
	// assuming standard LPJ allometry (see functions allocation, allometry).

	// INPUT PARAMETERS
	// bminit = initial total biomass (kgC)
	// ltor   = initial leaf:root biomass ratio
	//
	// Note: indiv.densindiv (density of individuals across patch or modelled area)
	//       should be set to a meaningful value before this function is called

	if (ifdavegrowth && ifdailytree && indiv.pft.lifeform == TREE) {
		establishment_sapling_dave(indiv, ltor);
	} else {
		establishment_biomass(bminit, ltor, indiv);
	}


	if ( (indiv.pft.is_grass() && ifdailygrass) || indiv.pft.dave_allocation()) {

		indiv.report_flux(Fluxes::ESTC, -(indiv.get_cmass_leaf() + indiv.cmass_root + indiv.cmass_storage + indiv.cmass_crown + indiv.cmass_sap));
		indiv.report_flux(Fluxes::ESTN, -indiv.initialn);
	}
}



///////////////////////////////////////////////////////////////////////////////////////
// ALLOMETRY
// Should be called to update allometry, FPC and FPC increment whenever biomass values
// for a vegetation individual change.

bool allometry(Individual& indiv) {

	// DESCRIPTION
	// Calculates tree allometry (height and crown area) and fractional projective
	// given carbon biomass in various compartments for an individual.

	// Returns true if the allometry is normal, otherwise false - guess2008

	// TREE ALLOMETRY
	// Trees aboveground allometry is modelled by a cylindrical stem comprising an
	// inner cylinder of heartwood surrounded by a zone of sapwood of constant radius,
	// and a crown (i.e. foliage) cylinder of known diameter. Sapwood and heartwood are
	// assumed to have the same, constant, density (wooddens). Tree height is related
	// to sapwood cross-sectional area by the relation:
	//   (1) height = cmass_sap / wooddens / (sapwood xs area)
	// Sapwood cross-sectional area is also assumed to be a constant proportion of
	// total leaf area (following the "pipe model"; Shinozaki et al. 1964a,b; Waring
	// et al 1982), i.e.
	//   (2) (leaf area) = k_latosa * (sapwood xs area)
	// Leaf area is related to leaf biomass by specific leaf area:
	//   (3) (leaf area) = sla * cmass_leaf
	// From (1), (2), (3),
	//   (4) height = cmass_sap / wooddens / sla / cmass_leaf * k_latosa
	// Tree height is related to stem diameter by the relation (Huang et al 1992)
	// [** = raised to the power of]:
	//   (5) height = k_allom2 * diam ** k_allom3
	// Crown area may be derived from stem diameter by the relation (Zeide 2001):
	//   (6) crownarea = min ( k_allom1 * diam ** k_rp , crownarea_max )
	// Bole height (individual/cohort mode only; currently set to 0):
	//   (7) boleht = 0

	// FOLIAR PROJECTIVE COVER (FPC)
	// The same formulation for FPC (Eqn 8 below) is now applied in all vegetation
	// modes (Ben Smith 2002-07-23). FPC is equivalent to fractional patch/grid cell
	// coverage for the purposes of canopy exchange calculations and, in population
	// mode, vegetation dynamics calculations.
	//
	//   FPC on the modelled area (stand, patch, "grid-cell") basis is related to mean
	//   individual leaf area index (LAI) by the Lambert-Beer law (Monsi & Saeki 1953,
	//   Prentice et al 1993) based on the assumption that success of a PFT population
	//   in competition for space will be proportional to competitive ability for light
	//   in the vertical profile of the forest canopy:
	//     (8) fpc = crownarea * densindiv * ( 1.0 - exp ( -0.5 * lai_ind ) )
	//   where
	//     (9) lai_ind = cmass_leaf/densindiv * sla / crownarea
	//
	//   For grasses,
	//    (10) fpc = ( 1.0 - exp ( -0.5 * lai_ind ) )
	//    (11) lai_ind = cmass_leaf * sla

	double fpc_new; // updated FPC

	// guess2008 - max tree height allowed (metre).
	const double HEIGHT_MAX = 150.0;


	if (indiv.pft.lifeform == TREE) {

		// TREES

		// Height (Eqn 4)

		// guess2008 - new allometry check
		if (!negligible(indiv.cmass_leaf)) {

			indiv.height = indiv.cmass_sap / indiv.cmass_leaf / indiv.sla() * indiv.pft.k_latosa / indiv.pft.wooddens;

			// Stem diameter (Eqn 5)
			indiv.diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);

			// Stem volume
			double vol = indiv.height * PI * indiv.diam * indiv.diam * 0.25;

			if (indiv.age && (indiv.cmass_heart + indiv.cmass_sap) / indiv.densindiv / vol < indiv.pft.wooddens * 0.9) {
				Patch& patch = indiv.vegetation.patch;
				Patchpft& patchpft = patch.pft[indiv.pft.id];

				patchpft.cmass_mort += indiv.ccont();
				return false;
			}
		}
		else {
			indiv.height = 0.0;
			indiv.diam = 0.0;
			return false;
		}


		// guess2008 - extra height check
		if (indiv.height > HEIGHT_MAX) {
			indiv.height = 0.0;
			indiv.diam = 0.0;
			return false;
		}


		// Crown area (Eqn 6)
		indiv.crownarea = min(indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp),
			indiv.pft.crownarea_max);

		if (!negligible(indiv.crownarea)) {

			// Individual LAI (Eqn 9)
			indiv.lai_indiv = indiv.cmass_leaf / indiv.densindiv *
				indiv.sla() / indiv.crownarea;

			// FPC (Eqn 8)

			fpc_new = indiv.crownarea * indiv.densindiv *
				(1.0 - lambertbeer(indiv.lai_indiv));

			// Increment deltafpc
			indiv.deltafpc += fpc_new - indiv.fpc;
			indiv.fpc = fpc_new;
			indiv.rpc = indiv.fpc;  // TODO Assaf: decide if rpc is needed for trees. RPC is used in Nuptake function.
		}
		else {
			indiv.lai_indiv = 0.0;
			indiv.fpc = 0.0;
		}

		// Bole height (Eqn 7)
		indiv.boleht = 0.5;

		// Stand-level LAI
		indiv.lai = indiv.cmass_leaf * indiv.sla();
	}
	else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {

		// GRASSES AND MOSSES

		if (indiv.pft.landcover != CROPLAND) {

			// guess2008 - bugfix - added if
			if (!negligible(indiv.get_cmass_leaf())) {

				// Grass "individual" LAI (Eqn 11). Green LAI only!
				indiv.lai_indiv = indiv.get_cmass_leaf() * indiv.sla();
				if (ifdavegrowth) {
					indiv.lai_indiv *= (1.0 - indiv.pft.culmfrac);
				}
				
				// FPC (Eqn 10)
				indiv.fpc = 1.0 - lambertbeer(indiv.lai_indiv);
				// Root projective cover (RPC) introduced to capture feedback between 
				// belowground C investment and available N. eff_rpc represents the efficiency of 
				// roots to cover volume, related to specfic root length (SRL in m g-1). SRL is ca. 4 times 
				// higher in annuals than perennials (Roumet et al. 2006, New Phytologist).
				double eff_rpc = indiv.pft.lifeform == GRASS_PERENNIAL ? 4.0 : 16.0; 
				indiv.rpc = 1.0 - lambertbeer(indiv.lai_indiv * min(10.0,(eff_rpc * indiv.cmass_root / indiv.cmass_leafs[green])));
				// TODO: Replace with better solution. Does this work with competition??

				// Stand-level LAI
				indiv.lai = indiv.lai_indiv;
			}
			else {
				indiv.lai_indiv = 0.0;
				indiv.fpc = 0.0;
				indiv.lai = 0.0;
				return false;
			}
		}
		else {

			// True crops use cmass_leaf_max, cover-crop grass uses lai of stands with whole-year grass growth
			allometry_crop(indiv);
		}
	}

	// guess2008 - new return value (was void)
	return true;
}



///////////////////////////////////////////////////////////////////////////////////////
// RELATIVE CHANGE IN BIOMASS
// Call this function to calculate the change in biomass on a grid cell area basis
// associated with a specified change in FPC

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv) {

	// DESCRIPTION
	// Calculates and returns new biomass as a fraction of old biomass given an FPC
	// reduction from fpc_high to fpc_low, assuming LPJ allometry (see function
	// allometry)

	// guess2008 - check
	if (fpc_high < fpc_low)
		fail("fracmass_lpj: fpc_high < fpc_low");

	if (indiv.pft.lifeform==TREE) {

		if (negligible(fpc_high)) return 1.0;

		// else
		return fpc_low/fpc_high;
	}
	else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {		// grass and moss

		if (fpc_high>=1.0 || fpc_low>=1.0 || negligible(indiv.cmass_leaf)) return 1.0;

		// else
		return 1.0+2.0/indiv.cmass_leaf/indiv.sla()*
			(log(1.0-fpc_high)-log(1.0-fpc_low));
	}
	else {
		fail("fracmass_lpj: unknown lifeform");
		return 0;
	}

	// This point will never be reached in practice, but to satisfy more pedantic
	// compilers ...

	return 1.0;
}

void flush_litter_repr(Patch& patch) {

	// Returns nitrogen-free reproduction "litter" to atmosphere

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& ppft = patch.pft.getobj();

		if (!ppft.pft.dave_allocation()) {
			// Updated soil fluxes
			patch.fluxes.report_flux(Fluxes::REPRC, ppft.cmass_litter_repr);
			ppft.cmass_litter_repr = 0.0;
		}

		patch.pft.nextobj();
	}
}

// Set standpft-level reproductive cmass to 0 if this is the first patch.
// This should be called once per year. Called from growth() iff dailygrass is
// disabled. Called from growth_daily_grass() on 1st day of year iff davegrowth
// is disabled. Called from growth_daily_grass_dave() on 1st day of year
// otherwise. Standpft-level reproductive cmass is used for tree establishment
// when the spatial mass effect (ifsme) is enabled, so failure to reset this to
// 0 results in absurdly large initial tree biomass at establishment. Resetting
// to 0 on 1st day of year kind of works because establishment happens on the
// last day of year (ie the day before).
void flush_cmass_repr(Patch &patch) {
	if (!patch.id) {
		for (int p=0; p < npft; p++) {
			patch.stand.pft[p].cmass_repr = 0.0;
		}
	}
}

/// GROWTH
/** Tissue turnover and allocation of fixed carbon to reproduction and new biomass
 *	Accumulated NPP (assimilation minus maintenance and growth respiration) on
 * 	patch or modelled area basis assumed to be given by 'anpp' member variable for
 *	each individual.
 *  Should be called by framework at the end of each simulation year for modelling
 *  of turnover, allocation and growth, prior to vegetation dynamics and disturbance
 */
void growth(Stand& stand, Patch& patch) {

	// minimum carbon mass allowed (kgC/m2)
	const double MINCMASS = 1.0e-8;
	// maximum carbon mass allowed (kgC/m2)
	const double MAXCMASS = 1.0e8;

	const double CDEBT_PAYBACK_RATE = 0.2;

	// carbon biomass increment (component of NPP available for production of
	// new biomass) for this time period on modelled area basis (kgC/m2)
	double bminc = 0.0;
	// C allocated to reproduction this time period on modelled area basis (kgC/m2)
	double cmass_repr = 0.0;
	// increment in leaf C biomass following allocation, on individual basis (kgC)
	double cmass_leaf_inc = 0.0;
	// increment in root C biomass following allocation, on individual basis (kgC)
	double cmass_root_inc = 0.0;
	// increment in sapwood C biomass following allocation, on individual basis (kgC)
	double cmass_sap_inc = 0.0;
	// increment in heartwood C biomass following allocation, on individual basis (kgC)
	double cmass_heart_inc = 0.0;
	// increment in heartwood C biomass following allocation, on individual basis (kgC)
	double cmass_debt_inc = 0.0;
	// increment in crop harvestable organ C biomass following allocation, on individual basis (kgC)
	double cmass_ho_inc = 0.0;
	// increment in crop above-ground pool C biomass following allocation, on individual basis (kgC)
	double cmass_agpool_inc = 0.0;
	// increment in crop stem C biomass following allocation, on individual basis (kgC)
	double cmass_stem_inc = 0.0;
	// increment in leaf litter following allocation, on individual basis (kgC)
	double litter_leaf_inc = 0.0;
	// increment in root litter following allocation, on individual basis (kgC)
	double litter_root_inc = 0.0;
	// negative increment that exceeds existing biomass following allocation, on individual basis (kgC)
	double exceeds_cmass;
	// C biomass of leaves in "excess" of set allocated last year to raingreen PFT last year (kgC/m2)
	double cmass_excess = 0.0;
	// Raingreen nitrogen demand for leaves dropped during the year
	double raingreen_ndemand;
	// Nitrogen stress scalar for leaf to root allocation
	double nscal;
	// Leaf C:N ratios before growth
	double cton_leaf_bg;
	// Root C:N ratios before growth
	double cton_root_bg;
	// Sap C:N ratios before growth
	double cton_sap_bg;

	double dval = 0.0;
	int p;

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation = patch.vegetation;
	Gridcell& gridcell = vegetation.patch.stand.get_gridcell();

	// On first call to function growth this year (patch #0), initialise stand-PFT
	// record of summed allocation to reproduction

	// Daily carbon allocation does this in growth_pasture_daily
	flush_cmass_repr(patch);

	double total_ndemand = 0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patchpft& patchpft = patch.pft[indiv.pft.id];

		if (indiv.pft.dave_allocation()) {
			// This individual is using dave allocation.
			// TODO: need to skip if it's a grass, and davegrowth is false but
			// ifdailygrass is true (ie the DW version).
			vegetation.nextobj();
			continue;
		}

		bool killed = false;

		// Standard yearly growth functionalities

		if (!(ifdailygrass && indiv.pft.is_grass()) || indiv.istruecrop_or_intercropgrass()) {
	
			// For this individual

			// Calculate vegetation carbon and nitrogen mass before growth to determine vegetation C:N ratios
			indiv.cmass_veg = indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_wood();
			indiv.nmass_veg = indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_wood();

			// Save real compartment C:N ratios before growth
			// phen is switched off for leaf and root
			// Leaf
			cton_leaf_bg = indiv.cton_leaf(false);

			// Root
			cton_root_bg = indiv.cton_root(false);

			// Sap
			cton_sap_bg = indiv.cton_sap();

			// Save leaf annual average C:N ratio
			if (!negligible(indiv.nday_leafon))
				indiv.cton_leaf_aavr /= indiv.nday_leafon;
			else
				indiv.cton_leaf_aavr = indiv.pft.cton_leaf_max;

			// Nitrogen stress scalar for leaf to root allocation (adopted from Zaehle and Friend 2010 SM eq 19)
			double cton_leaf_aopt = max(ifdavephotosyn ? indiv.cton_leaf() : indiv.cton_leaf_aopt, indiv.pft.cton_leaf_avr);

			if (ifnlim)
				nscal = min(1.0, cton_leaf_aopt / indiv.cton_leaf_aavr);
			else
				nscal = 1.0;

			// Set leaf:root mass ratio based on water stress parameter
			// or nitrogen stress scalar
			indiv.ltor = min(indiv.wscal_mean(), nscal) * indiv.pft.ltor_max;

			// Move leftover compartment nitrogen storage to longterm storage
			if(!indiv.has_daily_turnover())	{
				indiv.nstore_longterm += indiv.nstore_labile;
				indiv.nstore_labile = 0.0;
			}

			indiv.deltafpc = 0.0;

			if (negligible(indiv.densindiv))
				fail("growth: negligible densindiv for %s",(char*)indiv.pft.name);
			else {

				// Allocation to reproduction
				if(!indiv.istruecrop_or_intercropgrass())
					reproduction(indiv.pft.reprfrac,indiv.anpp,bminc,cmass_repr);

				if(indiv.alive)
					patchpft.cmass_repr += cmass_repr;

				raingreen_ndemand = 0.0;

				// added bminc check. Otherwise we get -ve litter_leaf for grasses when indiv.anpp < 0.
				//
				// cmass_excess-code inactivated for grass pft:s, due to frequent oscillations between high bminc and zero bminc in
				// certain grasslands using updated leaflong-value for grass-pft:s (0.5)
				if (bminc >= 0 && indiv.pft.phenology == RAINGREEN) {

					// Raingreen PFTs: reduce biomass increment to account for NPP
					// allocated to extra leaves during the past year.
					// Excess allocation to leaves given by:
					//   aphen_raingreen / ( leaf_longevity * year_length) * cmass_leaf -
					//   cmass_leaf
					
					// BLARP! excess allocation to roots now also included (assumes leaf longevity = root longevity)

					cmass_excess = max((double)indiv.aphen_raingreen /
							(indiv.pft.leaflong * (double)date.year_length()) * (indiv.cmass_leaf + indiv.cmass_root) -
							indiv.cmass_leaf - indiv.cmass_root, 0.0);

					if (cmass_excess > bminc) cmass_excess = bminc;

					// Transfer excess leaves to litter
					// only for 'alive' individuals
					if (indiv.alive) {
						patch.pft[indiv.pft.id].cmass_litter_leaf += cmass_excess;
						if (!negligible(cton_leaf_bg))
							raingreen_ndemand = min(indiv.nmass_leaf, cmass_excess / cton_leaf_bg);
						else
							raingreen_ndemand = 0.0;

						patch.pft[indiv.pft.id].nmass_litter_leaf += raingreen_ndemand * (1.0 - nrelocfrac);
						indiv.nstore_longterm += raingreen_ndemand * nrelocfrac;
						indiv.nmass_leaf -= raingreen_ndemand;
						assert(indiv.nmass_leaf >= 0);
					}

					// Deduct from this year's C biomass increment
					// added alive check
					if (indiv.alive) bminc -= cmass_excess;
					if (indiv.alive) patchpft.cmass_mort += cmass_excess;
				}

				// All yearly harvest events
				killed = harvest_year(indiv);

				if (!killed) {

					if(!indiv.has_daily_turnover()) {
						// Tissue turnover and associated litter production
						turnover(indiv.pft.turnover_leaf, indiv.pft.turnover_root,
								indiv.pft.turnover_sap, indiv.pft.lifeform, indiv.pft.landcover,
								indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap, indiv.cmass_heart,
								indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap, indiv.nmass_heart,
								patch.pft[indiv.pft.id].cmass_litter_leaf,
								patch.pft[indiv.pft.id].cmass_litter_root, patch.pft[indiv.pft.id].cmass_leaf_root_turnover,
								patch.pft[indiv.pft.id].nmass_litter_leaf,
								patch.pft[indiv.pft.id].nmass_litter_root,
								indiv.nstore_longterm,indiv.max_n_storage,
								indiv.alive);
					}
					// Update stand record of reproduction by this PFT
					stand.pft[indiv.pft.id].cmass_repr += cmass_repr / (double)stand.npatch();

					// Transfer reproduction straight to litter
					// only for 'alive' individuals
					if (indiv.alive) {
						patch.pft[indiv.pft.id].cmass_litter_repr += cmass_repr;
					}

					if (indiv.pft.lifeform == TREE) {

						// TREE GROWTH

						// BLARP! Try and pay back part of cdebt

						if (ifcdebt && bminc > 0.0) {
							double cmass_payback = min(indiv.cmass_debt * CDEBT_PAYBACK_RATE, bminc);
							bminc -= cmass_payback;
							indiv.cmass_debt -= cmass_payback;
						}

						// Allocation: note conversion of mass values from grid cell area
						// to individual basis

						allocation(bminc / indiv.densindiv, indiv, indiv.ltor,
								cmass_leaf_inc, cmass_root_inc, cmass_sap_inc, cmass_debt_inc,
								cmass_heart_inc, litter_leaf_inc, litter_root_inc, exceeds_cmass);

						// Update carbon pools and litter (on area basis)
						// (litter not accrued for not 'alive' individuals - Ben 2007-11-28)

						// Leaves
						indiv.cmass_leaf += cmass_leaf_inc * indiv.densindiv;

						// Roots
						indiv.cmass_root += cmass_root_inc * indiv.densindiv;

						// Sapwood
						indiv.cmass_sap += cmass_sap_inc * indiv.densindiv;

						// Heartwood
						indiv.cmass_heart += cmass_heart_inc * indiv.densindiv;

						indiv.cmass_wood_inc_5.add((cmass_sap_inc + cmass_heart_inc - cmass_debt_inc) * indiv.densindiv);

						// If negative sap growth, then nrelocfrac of nitrogen will go to heart wood and
						// (1.0 - nreloctrac) will go to storage
						double nmass_sap_inc = cmass_sap_inc * indiv.densindiv / cton_sap_bg;

						if (cmass_sap_inc < 0.0 && indiv.nmass_sap >= -nmass_sap_inc){
							indiv.nmass_sap += nmass_sap_inc;
							indiv.nmass_heart -= nmass_sap_inc * nrelocfrac;
							indiv.nstore_longterm -= nmass_sap_inc * (1.0 - nrelocfrac);
							assert(indiv.nmass_sap >= 0.0);
						}

						// C debt
						indiv.cmass_debt += cmass_debt_inc * indiv.densindiv;

						// alive check before ensuring C balance
						if (indiv.alive) {
							patch.pft[indiv.pft.id].cmass_litter_leaf += litter_leaf_inc * indiv.densindiv;
							patch.pft[indiv.pft.id].cmass_litter_root += litter_root_inc * indiv.densindiv;
							patchpft.cmass_mort += (litter_leaf_inc + litter_root_inc) * indiv.densindiv;

							// C litter exceeding existing biomass
							indiv.report_flux(Fluxes::NPP, exceeds_cmass * indiv.densindiv);
							indiv.report_flux(Fluxes::RA, -exceeds_cmass * indiv.densindiv);
						}

						double leaf_inc = min(indiv.nmass_leaf, litter_leaf_inc * indiv.densindiv / cton_leaf_bg);
						double root_inc = min(indiv.nmass_root, litter_root_inc * indiv.densindiv / cton_root_bg);

						// Nitrogen litter always return to soil litter and storage
						// Leaf
						patch.pft[indiv.pft.id].nmass_litter_leaf += leaf_inc * (1.0 - nrelocfrac);
						indiv.nstore_longterm += leaf_inc * nrelocfrac;

						// Root
						patch.pft[indiv.pft.id].nmass_litter_root += root_inc * (1.0 - nrelocfrac);
						indiv.nstore_longterm += root_inc * nrelocfrac;

						// Subtracting litter nitrogen from individuals
						indiv.nmass_leaf -= leaf_inc;
						indiv.nmass_root -= root_inc;

						assert(indiv.nmass_leaf >= 0);
						assert(indiv.nmass_root >= 0);

						// Update individual age

						indiv.age++;

						// Kill individual and transfer biomass to litter if any biomass
						// compartment negative

						if (indiv.cmass_leaf < MINCMASS || indiv.cmass_root < MINCMASS ||
								indiv.cmass_sap < MINCMASS) {

							if(indiv.alive)
								patchpft.cmass_mort += indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_wood();

							if (date.year >= nyear_spinup && _ifdavediagnostics) {
								printf("%04d-%02d-%02d (year %d) %s killed by carbon starvation...\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, date.year, (const char *)indiv.pft.name);
							}
							indiv.kill();

							vegetation.killobj();
							killed = true;
						}
					}
					else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {

					// GRASS AND MOSS GROWTH

						//True crops do not use bminc.or cmass_leaf etc.
						if(indiv.istruecrop_or_intercropgrass()) {
							// transfer crop cmass increase values to common variables
							growth_crop_year(indiv, cmass_leaf_inc, cmass_root_inc, cmass_ho_inc, cmass_agpool_inc, cmass_stem_inc);

							exceeds_cmass = 0.0;	//exceeds_cmass not used for true crops
						}
						else {
							allocation(bminc, indiv, indiv.ltor, cmass_leaf_inc, cmass_root_inc, dval, dval, dval,
									litter_leaf_inc, litter_root_inc, exceeds_cmass);
						}

						if(indiv.pft.landcover == CROPLAND) {
							if(indiv.istruecrop_or_intercropgrass())
								yield_crop(indiv);
							else
								yield_pasture(indiv, cmass_leaf_inc);
						}

						// Update carbon pools and litter (on area basis)
						// only litter in the case of 'alive' individuals

						// Leaves
						indiv.cmass_leaf += cmass_leaf_inc;

						// Roots
						indiv.cmass_root += cmass_root_inc;

						if(indiv.pft.landcover == CROPLAND)	{
							indiv.cropindiv->cmass_ho += cmass_ho_inc;
							indiv.cropindiv->cmass_agpool += cmass_agpool_inc;
							indiv.cropindiv->cmass_stem += cmass_stem_inc;
						}

						if(indiv.pft.phenology != CROPGREEN && !(indiv.has_daily_turnover() && indiv.continous_grass())) {

							// alive check before ensuring C balance
							if (indiv.alive && !indiv.istruecrop_or_intercropgrass()) {

								patch.pft[indiv.pft.id].cmass_litter_leaf += litter_leaf_inc;
								patch.pft[indiv.pft.id].cmass_litter_root += litter_root_inc;

								patchpft.cmass_mort += litter_leaf_inc + litter_root_inc;

								// C litter exceeding existing biomass
								indiv.report_flux(Fluxes::NPP, exceeds_cmass * indiv.densindiv);
								indiv.report_flux(Fluxes::RA, -exceeds_cmass * indiv.densindiv);
							}

							// Nitrogen always return to soil litter and storage
							// Leaf
							if (indiv.nmass_leaf > 0.0){
								double nmass = min(indiv.nmass_leaf, litter_leaf_inc * indiv.densindiv / cton_leaf_bg);

								patch.pft[indiv.pft.id].nmass_litter_leaf += nmass * (1.0 - nrelocfrac);
								indiv.nstore_longterm += nmass * nrelocfrac;

								// Subtracting litter nitrogen from individuals
								indiv.nmass_leaf -= nmass;
							}
							assert(indiv.nmass_leaf >= 0);

							// Root
							if (indiv.nmass_root > 0.0){
								double nmass = min(indiv.nmass_root, litter_root_inc * indiv.densindiv / cton_root_bg);

								patch.pft[indiv.pft.id].nmass_litter_root += nmass * (1.0 - nrelocfrac);
								indiv.nstore_longterm += nmass * nrelocfrac;
								// Subtracting litter nitrogen from individuals
								indiv.nmass_root -= nmass;
							}
						}
						// Kill individual and transfer biomass to litter if either biomass
						// compartment negative

						// Note - what happens if cmass_root = 0?
						if ((indiv.cmass_leaf < MINCMASS || indiv.cmass_root < MINCMASS) && !indiv.istruecrop_or_intercropgrass()) {

							if(indiv.alive)
								patchpft.cmass_mort += indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_wood();

							if (date.year >= nyear_spinup && _ifdavediagnostics) {
								printf("%04d-%02d-%02d (year %d) %s killed by carbon starvation...\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, date.year, (const char *)indiv.pft.name);
							}
							indiv.kill();

							vegetation.killobj();
							killed = true;
						}

					} // grass or moss

					if (!killed && indiv.pft.phenology != CROPGREEN && !(indiv.has_daily_turnover() && indiv.continous_grass())) {
						// Update nitrogen longtime storage

						// Nitrogen approx retranslocated next year
						double retransn_nextyear = indiv.cmass_leaf * indiv.pft.turnover_leaf / cton_leaf_bg * nrelocfrac +
								indiv.cmass_root * indiv.pft.turnover_root / cton_root_bg * nrelocfrac;

						if (indiv.pft.lifeform == TREE)
							retransn_nextyear += indiv.cmass_sap * indiv.pft.turnover_sap / cton_sap_bg * nrelocfrac;

						// Assume that raingreen will lose same amount of N through extra leaves next year
						if (indiv.alive && bminc >= 0 && (indiv.pft.phenology == RAINGREEN || indiv.pft.phenology == ANY))
							retransn_nextyear -= min(raingreen_ndemand, retransn_nextyear);

						// Max longterm nitrogen storage
						if (indiv.pft.lifeform == TREE)
							indiv.max_n_storage = max(0.0, max(indiv.cmass_sap * indiv.pft.fnstorage / cton_leaf_bg, retransn_nextyear));
						else // GRASS
							indiv.max_n_storage = max(0.0, max(indiv.cmass_root * indiv.pft.fnstorage / cton_leaf_bg, retransn_nextyear));

						// Scale this year productivity to max storage
						if (indiv.anpp > 0.0) {
							indiv.scale_n_storage = max(indiv.max_n_storage * 0.1, indiv.max_n_storage - retransn_nextyear) * cton_leaf_bg / indiv.anpp;
						} // else use last years scaling factor
					}
				} else if (ifdavephotosyn) {
					// This tree was killed and newphotosyn is enabled. Let's reset
					// ndemand variables just to be on the safe side.
					indiv.ndemand = 0;
					indiv.andemand = 0;
					indiv.andemand_leaf = 0;
					indiv.andemand_root = 0;
					indiv.andemand_sap = 0;
				}
			}


			if (!killed) {

				if (!allometry(indiv)) {

					if (date.year >= nyear_spinup && _ifdavediagnostics) {
						printf("%04d-%02d-%02d (year %d) %s killed by carbon starvation...\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, date.year, (const char *)indiv.pft.name);
					}
					indiv.kill();

					vegetation.killobj();
					killed = true;
				}

				if (!killed) {
					if (!indiv.alive) {
						// The individual has survived its first year...
						indiv.alive = true;

						// ...now we can start counting its fluxes,
						// debit current biomass as establishment flux
						if (!indiv.istruecrop_or_intercropgrass()) {
							indiv.report_flux(Fluxes::ESTC,
							              - (indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap +
							              indiv.cmass_heart - indiv.cmass_debt));
							indiv.report_flux(Fluxes::ESTN, -indiv.initialn);

							patchpft.cmass_est -= (indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap +
							              indiv.cmass_heart - indiv.cmass_debt);
						}
					}

					// Move long-term nitrogen storage pool to labile storage pool for usage next year
					if(!indiv.has_daily_turnover()) {
						indiv.nstore_labile = indiv.nstore_longterm;
						indiv.nstore_longterm = 0.0;
					}

					// ... on to next individual
					vegetation.nextobj();
				}
			}
		} 
		else if (ifdailygrass && indiv.pft.is_grass()) {
			
			// Check daily grass in the end of the year after all processes has been completed

			if (negligible(indiv.densindiv)) {
				fail("growth: negligible densindiv for %s", (char*)indiv.pft.name);
			}
			else {

				// Kill individual and transfer biomass to litter if either biomass
				// compartment negative or if yearly growth this year has been low

				if (indiv.cmass_leaf_ymax < MINCMASS || indiv.cmass_root_ymax < MINCMASS || indiv.cmass_leaf_ygrowth < MINCMASS) {

					if (ifdavegrowth) {
						indiv.flush_repr_litter_DAVE();
					}

					if (date.year >= nyear_spinup && _ifdavediagnostics) {
						printf("%04d-%02d-%02d (year %d) %s killed by carbon starvation...\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, date.year, (const char *)indiv.pft.name);
					}
					indiv.kill();

					vegetation.killobj();
					killed = true;
				}
				else {
					// ... on to next individual
					vegetation.nextobj();
				}
			}
		}
		else {  // not tree and dcarb turned on

			// ... on to next individual
			if (!killed) {
				vegetation.nextobj();
			}
		}
	}

	// Flush nitrogen free litter from reproduction straight to atmosphere
	flush_litter_repr(patch);
}

/// GROWTH DAILY GRASS
/** Daily grass tissue turnover and allocation of fixed carbon to reproduction and new biomass
 *  Should be called by framework at the end of each simulation day for modelling
 *  of turnover, allocation and growth.
 *  Daily grass growth follow the model of Johnson & Thornley (1983)
 */
//
void growth_daily_grass(Patch& patch) {

	Stand& stand = patch.stand;

	// transfer constant between first and second compartment constant, to make that transfer faster following Johnson & Thornley 1983
	const double LAMBDA = 2.0;
	// maximum daily growth rate
	const double MAX_GROWTH_RATE = 0.2;
	// maximum gpp trend growth modifier
	const double GPP_FAC_MAX = 1.0;
	// storage to biomass ratio growth modifier exponent
	const int SB_FAC_POW = 2;
	// storage to biomass ratio growth modifier maximum
	const double SB_FAC_MAX = 5.0;
	// water scalar minimum for growth
	const double WSCAL_MIN = 0.1;
	// fraction of N limitation that C growth will experience
	const double FRAC_C_OF_N_LIM = 0.8;
	// C retranslocation fraction for C3 grass senescense
	const double C3_RETRANS = 0.1;
	// C retranslocation fraction for C4 grass senescense
	const double C4_RETRANS = 0.05;

	// new biomass for this time period on modelled area basis (kgC/m2)
	double bminc = 0.0;
	// C allocated to reproduction this time period on modelled area basis (kgC/m2)
	double cmass_repr = 0.0;
	// leaf C growth this time period
	double cgrow_leaf = 0.0;
	// leaf N growth this time period
	double ngrow_leaf = 0.0;
	// leaf N update this time period
	double nup_leaf = 0.0;
	// root C growth this time period
	double cgrow_root = 0.0;
	// root N growth this time period
	double ngrow_root = 0.0;
	// Retranslocation fraction of C
	double crelocfrac_daily;
	// Retranslocation fraction of N
	double nrelocfrac_daily;

	// Movment factors for growth
	double g_fac;
	// Movment factors for senescense and movement to litter
	double s_fac;

	// Temporary variables to move carbon and nitrogen between compartments
	double c1, c2, c3, c4;
	double n1, n2, n3, n4;

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation = patch.vegetation;
	Gridcell& gridcell = vegetation.patch.stand.get_gridcell();

	// Air temperature or bottom snow layer if snow present
	double grass_temp = patch.soil.T_soil[6];

	if (date.day == 0) {
		flush_cmass_repr(patch);
	}

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		Patchpft& ppft = patch.pft[indiv.pft.id];

		bool killed=false;

		if (indiv.pft.is_grass() && !indiv.istruecrop_or_intercropgrass()) {
			
			cgrow_leaf = cgrow_root = ngrow_leaf = nup_leaf = ngrow_root = 0.0;
			c1 = c2 = c3 = c4 = 0.0;
			n1 = n2 = n3 = n4 = 0.0;

			// Initiate yearly maximum
			if (date.day == 0) {
				indiv.cmass_leaf_ygrowth = 0.0;
				indiv.cmass_root_ygrowth = 0.0;
				indiv.lai_ymax = 0.0;
				indiv.fpc_ymax = 0.0;
				indiv.cmass_leaf_ymax = 0.0;
				indiv.nmass_leaf_ymax = 0.0;
				indiv.cmass_root_ymax = 0.0;
				indiv.nmass_root_ymax = 0.0;

				// scale indiv for landcover change
				scale_indiv(indiv, false);
			}

			// Reproduction cost
			reproduction(indiv.pft.reprfrac, indiv.dnpp, bminc, cmass_repr);

			// Update pft reproduction and flush it directly to litter
			if (cmass_repr > 0.0) {
				stand.pft[indiv.pft.id].cmass_repr += cmass_repr / (double)stand.npatch();
				ppft.cmass_litter_repr += cmass_repr;
				flush_litter_repr(patch);
			}

			// Check C budget
			double initc = bminc + indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_storage +
				indiv.cmass_root + ppft.cmass_litter_leaf + ppft.cmass_litter_root;

			// Update state of biomass and storage pools

			// Storage 
			// Positive bminc is added to storage pool for retranslocation out to leaves/roots further down
			if (bminc >= 0.0) {
				indiv.cmass_storage += bminc;
			}
			else {

				// Negative bminc

				// Respiration initially from labile (storage) pool
				// Labile pool can meet demand
				if (indiv.cmass_storage + bminc > 0.0) {
					indiv.cmass_storage += bminc;
					bminc = 0.0;
				}
				else {
					bminc += indiv.cmass_storage;
					indiv.cmass_storage = 0.0;
				}
				
				// If bminc still negative respiration from existing leaf and root
				// biomass and retranslocate corresponding N to keep leaf and root 
				// C:N ratio

				if (bminc < 0.0) {
					
					// Enough C to compensate for negative C storage
					if (indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_root + bminc > 0.0) {

						// Root fraction to be reduced
						double froot_dec = 0.0;
						if (indiv.cmass_root > 0.0) {
							double current_ltor = indiv.cmass_root > 0 ? (indiv.get_cmass_leaf() + indiv.cmass_leafs[3]) / indiv.cmass_root : 0.0;
							froot_dec = 1.0 / (current_ltor + 1.0);
						}

						// Reduce root
						double ntoc = indiv.nmass_root / indiv.cmass_root;
						double cdec = -bminc * froot_dec;
						double ndec = cdec * ntoc;
						indiv.cmass_root -= cdec;
						indiv.nmass_root -= ndec;
						bminc += cdec;
						indiv.nmass_storage += ndec;

						// Reduce leaf
						int comp = 3;
						while (bminc < 0.0 && comp > -1) {
							if (indiv.cmass_leafs[comp] > 0.0) {
								double ntoc = indiv.nmass_leafs[comp] / indiv.cmass_leafs[comp];
								double cdec = min(indiv.cmass_leafs[comp], -bminc);
								double ndec = cdec * ntoc;
								indiv.cmass_leafs[comp] -= cdec;
								indiv.nmass_leafs[comp] -= ndec;
								bminc += cdec;
								indiv.nmass_storage += ndec;
							}
							comp--;
						}
					}
					else {
						// Not enough C to compensate (should never happen)
						double diff = indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_root + bminc;
						for (int comp = 0; comp < indiv.cmass_leafs.size(); comp++)
							indiv.cmass_leafs[comp] = 0.0;
						indiv.cmass_root = 0.0;
						bminc = 0.0;
						fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
					}
				}
			}

			// If no N limitation or free N years add N to match minimum C:N ratio of the storage pool
			// Added N is reported as fertilisation to keep N mass balance
			if ((!ifnlim || date.year < freenyears) && indiv.nmass_storage < indiv.cmass_storage / indiv.pft.cton_leaf_min) {
				double nmass_missing = indiv.cmass_storage / indiv.pft.cton_leaf_min - indiv.nmass_storage;
				indiv.nmass_storage += nmass_missing;
				patch.fluxes.report_flux(Fluxes::NFERT, nmass_missing);
			}

			// Calculate leaf to root ratio based on N and water limitation and current leaf to root status 
			// or only water limitation and current leaf to root status if N limitation is turned off

			// Current leaf to root status
			double ltor_status = (!is_zero(indiv.get_cmass_leaf()) && !is_zero(indiv.cmass_root)) ? max(0.5, min(2.0, indiv.pft.ltor_max / (indiv.get_cmass_leaf() / indiv.cmass_root))) : 1.0;
			double ltor_min = 0.5 * indiv.pft.ltor_max;
			if (ifnlim && date.year > freenyears) {
				if (!negligible(indiv.get_cmass_leaf())) {
					indiv.nscal_5.add(min(1.0, indiv.cton_leaf_opt / indiv.cton_leaf()));
				}
				
				// leaf to root ratio of new growth (Johnson 2008 eqn 4.13b)
				//indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * min(indiv.wscal_running.mean(), indiv.nscal_running.mean()) * ltor_status;
				indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * min(ppft.wscal_5.mean(), indiv.nscal_5.mean()) * ltor_status;
			}
			else {
				indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * ppft.wscal_5.mean() * ltor_status;
			}
			
			// Growth rate modifiers

			// GPP trend growth modifier 
			// Negative gpp trend less growth more to storage
			double gpp_trend = indiv.gpp_15.trend();
			double gpp_fac = gpp_trend > 0.0 ? GPP_FAC_MAX : GPP_FAC_MAX * 0.05;
			indiv.gpptrend_15.add(gpp_fac);
			indiv.gpp_fac = indiv.gpptrend_15.mean();

			// Storage to biomass growth modifier (adopted from eqn 13 Johnson & Thornley 1983)
			// Lower storage fraction less growth
			double sb_ratio = !negligible(indiv.get_cmass_leaf() + indiv.cmass_root) ? indiv.cmass_storage / (indiv.get_cmass_leaf() + indiv.cmass_root) : SB_FAC_MAX;
			indiv.sb_fac = min(SB_FAC_MAX, pow(sb_ratio, SB_FAC_POW));

			// Fraction of leaves that die during one day due to frost (eqn 2 and 3 H�glind 2016, error in paper eqn)
			// Temperature where half the leaves would be killed. Should be flexible due to hardening, but here set to a fixed value
			// also should use snow temperature if snowpack covers grass
			double lethal_temp_50 = -15.0;
			indiv.frost_sen = grass_temp < 0.0 ? 1.0 - max(0.5, 1.0 / (1.0 + exp(lethal_temp_50 - grass_temp))) : 0.0;
			
			// Water stress growth modifier 
			indiv.wscal_fac = max(0.0, (ppft.wscal_5.mean() - WSCAL_MIN) / (1.0 - WSCAL_MIN));
			
			// Temperature growth modifier
			double low_temp = (indiv.pft.pstemp_min + indiv.pft.pstemp_low) / 2.0;
			double min_temp = max(0.0, indiv.pft.pstemp_min);
			indiv.tscal_fac = grass_temp > indiv.pft.pstemp_min ? (1.0 - 0.01 * exp(4.6 * (grass_temp - indiv.pft.pstemp_high) / (indiv.pft.pstemp_max - indiv.pft.pstemp_high))) /
				(1.0 + exp(4.6 * (low_temp - grass_temp) / (low_temp - indiv.pft.pstemp_min))) : 0.0;
			
			// Senescence modifiers

			// Water limitation increases senescence
			indiv.wscal_sen = 1.0 + 1.0 * (1.0 - indiv.wscal_fac);

			// Temperature senescence modifier
			indiv.tscal_sen = 1.0 + 0.5 * (1.0 - indiv.tscal_fac);

			// Growth and senescence rates (adopted from eqn 24 Johnson & Thornley 1983, H�glind 2016)
			double const g_max = 1.0; // CALIBRATE
			indiv.g_fac = min(MAX_GROWTH_RATE, g_max * indiv.gpp_fac * indiv.sb_fac * indiv.wscal_fac * indiv.tscal_fac);
			indiv.s_fac = max(indiv.frost_sen, indiv.pft.transferconst * indiv.wscal_sen * indiv.tscal_sen);

			// Storage used for new leaf and root growth and updating old leaves and roots
			if (!negligible(indiv.cmass_storage) && !negligible(indiv.nmass_storage)) {

				// Todays growth
				double cgrow = indiv.g_fac * indiv.cmass_storage;

				// Leaf and root C growth
				cgrow_root = cgrow / (1.0 + indiv.ltor);
				cgrow_leaf += cgrow - cgrow_root;

				// N demand for growth
				double ngrow_leaf_demand = cgrow_leaf / indiv.cton_leaf_opt;
				double ngrow_root_demand = cgrow_root / (indiv.cton_leaf_opt * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);

				// N demand for update current tissue to optimal value (10% day-1)
				double nup_leaf_demand = 0.0;
				double nup_root_demand = 0.0;
				if (!negligible(indiv.nmass_storage) && indiv.cmass_storage / indiv.nmass_storage < indiv.cton_leaf_opt) {
					nup_leaf_demand = max(0.0, indiv.get_cmass_leaf() / indiv.cton_leaf_opt - indiv.get_nmass_leaf()) * 0.1;
					nup_root_demand = max(0.0, indiv.cmass_root / (indiv.cton_leaf_opt * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr) - indiv.nmass_root) * 0.1;
				}

				// Total N demand
				double ndemand_tot = ngrow_leaf_demand + ngrow_root_demand + nup_leaf_demand + nup_root_demand;

				// N limitation on growth
				double nlim = 1.0;
				if (ndemand_tot > 1.0e-12) {
					nlim = min(1.0, MAX_GROWTH_RATE * indiv.nmass_storage / ndemand_tot);
				}

				// Too small N growth of too strong N limitation-> no growth
				if (ndemand_tot <= 1.0e-12 || nlim < 1.0e-6) {
					nlim = 1.0;
					cgrow = cgrow_root = cgrow_leaf = 0.0;
					ngrow_leaf = ngrow_root = nup_leaf = 0.0;
					ngrow_leaf_demand = ngrow_root_demand = nup_leaf_demand = nup_root_demand = 0.0;
				}
				// No new growth but N demand from current tissue
				else if (is_zero(cgrow) && nup_leaf_demand + nup_root_demand > 1.0e-12) {
					cgrow = cgrow_root = cgrow_leaf = 0.0;
					ngrow_leaf_demand = ngrow_root_demand = 0.0;

					// N limitation on N update on current biomass
					double nup_lim = min(1.0, (MAX_GROWTH_RATE * indiv.nmass_storage) / (nup_leaf_demand + nup_root_demand));
					ngrow_leaf = 0.0; 
					ngrow_root = nup_root_demand * nup_lim;
					nup_leaf = nup_leaf_demand * nup_lim;
				}
				// N limitation on growth
				else if (nlim < 1.0 && ifnlim && date.year > freenyears) {

					// Optimal C:N ratio of growth
					double cton_grow_opt = cgrow / ndemand_tot;

					// New growth C:N ratio without any effect on C growth
					double cton_grow_fix = cgrow / (nlim * ndemand_tot);

					// New growth C:N ratio with effect on C growth
					double cton_grow_new = 1.0 / (FRAC_C_OF_N_LIM / cton_grow_fix + (1.0 - FRAC_C_OF_N_LIM) / cton_grow_opt);

					// New fraction of storage N used for growth
					double fNmax = MAX_GROWTH_RATE * min(1.0, (indiv.nmass_storage / indiv.cmass_storage + 1.0 / indiv.pft.cton_leaf_min) / (2.0 / indiv.pft.cton_leaf_min));
					
					// N growth
					double ngrow = fNmax * indiv.nmass_storage;

					// C growth limitation
					double cgrow_lim = cton_grow_new * min(ndemand_tot, ngrow) / cgrow;

					// New leaf and root C growth
					cgrow *= cgrow_lim;
					cgrow_root *= cgrow_lim;
					cgrow_leaf *= cgrow_lim;

					// N growth limitation
					double ngrow_lim = ngrow / ndemand_tot;
					
					// New leaf and root N growth and update
					ngrow_leaf = ngrow_leaf_demand * ngrow_lim;
					ngrow_root = (ngrow_root_demand + nup_root_demand) * ngrow_lim;
					nup_leaf = nup_leaf_demand * ngrow_lim;
				}
				// Full growth with no N limitation
				else {
					ngrow_leaf = ngrow_leaf_demand;
					ngrow_root = (ngrow_root_demand + nup_root_demand);
					nup_leaf = nup_leaf_demand;
				}

				// Leaf N mass update according to compartment size
				if (!negligible(indiv.get_cmass_leaf()) && !negligible(nup_leaf)) {
					for (int lf = 0; lf < indiv.cmass_leafs.size() - 1; lf++) {
						double fnup = indiv.cmass_leafs[lf] / indiv.get_cmass_leaf();
						indiv.nmass_leafs[lf] += nup_leaf * fnup;
					}
				}
			}

			// Save yearly leaf and root growth
			indiv.cmass_leaf_ygrowth += cgrow_leaf;
			indiv.cmass_root_ygrowth += cgrow_root;

			// Root growth
			indiv.cmass_root += cgrow_root;
			indiv.nmass_root += ngrow_root;
			
			// Storage reduction
			indiv.cmass_storage -= cgrow_leaf + cgrow_root;
			indiv.nmass_storage -= ngrow_leaf + nup_leaf + ngrow_root;

			// Calculate daily C and N retranslocation fraction
			crelocfrac_daily = indiv.pft.pathway == C3 ? C3_RETRANS : C4_RETRANS;
			nrelocfrac_daily = calc_nrelocfrac_daily(indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max, indiv.cmass_storage, indiv.nmass_storage);

			// Movement variables for C and N (eqn 26a-d Johnson & Thornley 1983)
			c1 = min(1.0, LAMBDA * indiv.s_fac) * indiv.cmass_leafs[0];
			n1 = min(1.0, LAMBDA * indiv.s_fac) * indiv.nmass_leafs[0];
			c2 = indiv.s_fac * indiv.cmass_leafs[1];
			n2 = indiv.s_fac * indiv.nmass_leafs[1];
			c3 = indiv.s_fac * indiv.cmass_leafs[2];
			n3 = indiv.s_fac * indiv.nmass_leafs[2];
			c4 = indiv.s_fac * indiv.cmass_leafs[3];
			n4 = indiv.s_fac * indiv.nmass_leafs[3];

			// Growing leaves (eqn 26a Johnson & Thornley 1983)
			indiv.cmass_leafs[0] += cgrow_leaf - c1;
			indiv.nmass_leafs[0] += ngrow_leaf - n1;

			// First fully expanded leaves (eqn 26b Johnson & Thornley 1983)
			indiv.cmass_leafs[1] += c1 - c2;
			indiv.nmass_leafs[1] += n1 - n2;

			// Second fully expanded leaves (eqn 26c Johnson & Thornley 1983)
			indiv.cmass_leafs[2] += c2 - c3;
			indiv.nmass_leafs[2] += n2 - n3;

			// Senescing leaves (eqn 26d Johnson & Thornley 1983)
			indiv.cmass_leafs[3] += (1.0 - crelocfrac_daily) * c3 - c4;
			indiv.cmass_storage += crelocfrac_daily * c3;
			indiv.nmass_leafs[3] += (1.0 - nrelocfrac_daily) * n3 - n4;
			indiv.nmass_storage += nrelocfrac_daily * n3;

			// Grazing
			grazing_daily(indiv);

			// Litter - if we have leaf litter then all compartments have litter to be dropped.
			if (c4 > 0.0) {

				patch.is_litter_day = true;

				// Maximum turnover rate (to prevent biomass becoming negligible)
				double turnover_max = 0.1;

				// Leaf turnover
				ppft.cmass_litter_leaf += c4;
				ppft.nmass_litter_leaf += n4;

				// Daily leaf and root turnover
				double total_cmass_leaf = indiv.get_cmass_leaf() + c3;
				double leaves_turnover_frac = min(turnover_max, !negligible(total_cmass_leaf) ? c3 / total_cmass_leaf : 0.0);

				// Root C turnover (use same fraction as for leaves)
				double turnov = leaves_turnover_frac * indiv.cmass_root;
				indiv.cmass_root -= turnov;
				ppft.cmass_litter_root += turnov;

				// Root N turnover
				turnov = leaves_turnover_frac * indiv.nmass_root;
				indiv.nmass_root -= turnov;
				ppft.nmass_litter_root += turnov * (1.0 - nrelocfrac_daily);
				indiv.nmass_storage += turnov * nrelocfrac_daily;
			}

			// C budget check
			double finalc = indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_storage +
				indiv.cmass_root + ppft.cmass_litter_leaf + ppft.cmass_litter_root;

			if (fabs(initc - finalc) > 1.0e-12)
				dprintf("Year %d Day %d Diff %g Init %g Final %g\n",date.year,date.day, initc - finalc, initc, finalc);

			// Update allometry - updates LAI and FPC based on cmass
			allometry(indiv);

			// Determine the yearly maximum values which is reported
			indiv.lai_ymax = max(indiv.lai_ymax, indiv.lai);
			indiv.fpc_ymax = max(indiv.fpc_ymax, indiv.fpc);

			if (indiv.get_cmass_leaf() > indiv.cmass_leaf_ymax) {
				indiv.cmass_leaf_ymax = indiv.get_cmass_leaf();
				indiv.nmass_leaf_ymax = indiv.get_nmass_leaf();
			}

			if (indiv.cmass_root > indiv.cmass_root_ymax) {
				indiv.cmass_root_ymax = indiv.cmass_root;
				indiv.nmass_root_ymax = indiv.nmass_root;
			}
			
			// Report yearly data
			if (date.islastday && date.islastmonth) {

				harvest_year(indiv);

				// Calculate vegetation carbon and nitrogen mass and average C:N ratios for yearly outputs
				indiv.cmass_veg = indiv.cmass_leaf_ymax + indiv.cmass_root_ymax;
				indiv.nmass_veg = indiv.nmass_leaf_ymax + indiv.nmass_root_ymax;

				// Save leaf annual average C:N ratio
				if (!negligible(indiv.nday_leafon)) {
					indiv.cton_leaf_aavr /= indiv.nday_leafon;
				}
				else {
					indiv.cton_leaf_aavr = indiv.pft.cton_leaf_max;
				}

				indiv.deltafpc = 0.0;
			} // Last day of year
		} // Grass

		if (!killed) {
			vegetation.nextobj();
		}

	} // while vegetation is object

} // growth_daily() 





// --------------------------------------
// New DAVE functions
// --------------------------------------

// Determine growth demand (aka 'sink strength') in kgC m-2 d-1
double growth_demand(Individual& indiv) {

	Vegetation& veg = indiv.vegetation;
	Pft& pft = indiv.pft;
    Patchpft &ppft = veg.patch.pft[indiv.pft.id]; 

	// Reduction of RGR in reproductive phase
	const double fRGR_repr = 0.75;
	// Sensitivity of water stress to growth
	// growth is more sensitive to water stress than photosynthesis (see Fatichi et al. 2014)
	const double qgrowth = 1.5;

	double fwater_growth = 0.0;
	double ftemp_growth  = 0.0;
	double fphoto_growth = 0.0;
	double growth_demand = 0.0;

    // ------------------------------
	// Calculate ancillary metrics
	// ------------------------------
	// Air temperature or bottom snow layer if snow present
	double grass_temp = veg.patch.soil.T_soil[6];

	// Total biomass used for estimation of potential growth demand.
	//double cmass_total  = indiv.cmass_root + 0.2 * indiv.cmass_crown + indiv.get_cmass_leaf();
	double cmass_total  = indiv.cmass_root + indiv.cmass_storage + indiv.get_cmass_leaf();

	// Fill status of C and N storage
    double cstorage_fill = max(0.0, (indiv.cmass_storage - (frac_cstorage_res * indiv.cstorage_max)) / (indiv.cstorage_max - frac_cstorage_res * indiv.cstorage_max));
    double nstorage_fill = max(0.0, (indiv.nmass_storage / indiv.nstorage_max));


    // ------------------------------
	// Calculate growth demand
	// ------------------------------
	// Factors determining growth demand
	fwater_growth = wscal_growth(ppft.wscal, qgrowth);
	ftemp_growth  = temp_growth(grass_temp, pft.Ha_growth, pft.T0_growth, pft.A_growth, pft.temp_dormancy);
	fphoto_growth = photo_growth(ppft.reldaylength, ppft.daylength_10.periodictrend(1));  // set to 1.0 
	//dev_growth(indiv, indiv.frepr_cstruct); // indiv.fdev_growth now calculated in phenology_grass_dave beforehand!
	double frepr = pft.is_grass() ? indiv.frepr_cstruct : fRGR_repr; // TODO Assaf: is this needed here?

	// NOTE: alternatively, RGR_max could be found as a function of sla!

	growth_demand = indiv.pft.RGR_pot * cmass_total * indiv.fdev_growth * ftemp_growth * fwater_growth * fphoto_growth;

	// feedbacks with C storage and N storage
	// i.e. reduce growth demand when C or N are scarce.
	//growth_demand *= (1.0 - exp(-20.0 * cstorage_fill));
	double qfeedback = indiv.growth_stage == vegetative ? 0.1 : 1.5 - indiv.fdev_growth;
	growth_demand *= pow(cstorage_fill,qfeedback);
	growth_demand *= pow(nstorage_fill,qfeedback);

	// reduced growth under dormancy (total cessation in full dormancy)
	if (indiv.fphen_activity < 1.0){
		growth_demand *= (indiv.fphen_activity - fphen_activity_min) / (1.0 - fphen_activity_min);
	}

	return growth_demand;

}

void storage_growth_partitioning_grass(Individual &indiv) {
	// Local:
	double cstore = 0.0;

	// C storage not available for growth.
	double cstorage_res = frac_cstorage_res * indiv.cstorage_max;

	// C source: Potential growth rate (kgC m-2 d-1)
	double Gpot = max(0.0, max_mobil_cstorage * max(0.0,(indiv.cmass_storage - cstorage_res)) + indiv.dnpp); 

	// C sink: growth demand (kgC m-2 d-1).
	indiv.growth_demand = growth_demand(indiv);

	if (indiv.dnpp >= 0.0) {
		
		// Positive NPP

		// TODO: make clear this only works for daily time steps (NPP in kgC m-2 d-1)	
		indiv.cgrow = min(indiv.growth_demand,Gpot);
		cstore = indiv.dnpp - indiv.cgrow;  // = Gpot - growth - (Gpot - NPP) 
											// cstore = re-allocation if negative

		// limit allocation to C storage if it is almost at capacity
		double crefill_max = max(0.0, (indiv.cstorage_max - indiv.cmass_storage) * max_mobil_cstorage);
		if (cstore > crefill_max) {
			indiv.cexcess += cstore - crefill_max;
			cstore = crefill_max; 
		}

		// If C storage exceed its capacity take track of the excess C and reduce cstore
		//if (indiv.cmass_storage + cstore > indiv.cstorage_max) {
	    //	indiv.cexcess += indiv.cmass_storage + cstore - indiv.cstorage_max;
		//	cstore = indiv.cstorage_max - indiv.cmass_storage;
		//}

	} else {

		// Negative NPP

		// Respiration initially from storage (labile) pool
		// Storage pool can meet demand. Still check if growth should occur.
		if (indiv.cmass_storage + indiv.dnpp >= 0.0) {
			indiv.cgrow = min(indiv.growth_demand,Gpot);
			cstore = indiv.dnpp - indiv.cgrow;
		}
		else {
			
			// Not enough C in storage to fulfill respiration demands --> kill individual
			double cmass_deficit = -(indiv.cmass_storage + indiv.dnpp);
			indiv.cmass_storage = 0.0;

			// The respiration and NPP fluxes reported in canexch are
			// now incorrect, so we need to fix them. Note that the RA
			// flux is not used for mass balance checks (it's only used
			// for outputs). The NPP flux is used in mass balance checks
			// however, so it's important to correct it here.
			indiv.report_flux(Fluxes::RA, -cmass_deficit);
				// (npp has opposite sign to respiration)
			indiv.report_flux(Fluxes::NPP, cmass_deficit);

			indiv.kill();  // flushes all pools to litter
			indiv.vegetation.killobj();  // removes individual object
								
			//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
			// Don't think we need to call nextobj() here, as we've just
			// removed an element from the list; when this happens, the
			// "current" pointer internal to the list object points to
			// the element that was previously the "next" element.
			return;
		}
	}

	assert(positive(indiv.cmass_storage));
	assert(equal(indiv.dnpp, indiv.cgrow + cstore + indiv.cexcess));

	// update cmass_storage
	indiv.cmass_storage += cstore;
}

// finds the carbon pool optiniums 
void adjust_cpool_optimum(Individual &indiv) {
	// calculate leaf limit (optimum)
	const double MAX_HEIGHT = 150;
	const double MIN_HEIGHT = 0.05;
	// limit height to below max and above min (to avoid FNS > 1)(from David's 
	// code). 
	// TODO: what about seedlings????
	double height = min(max(MIN_HEIGHT, indiv.height), MAX_HEIGHT);
	double fns = indiv.pft.k_latosa / (indiv.pft.sla * indiv.pft.wooddens * height); 
	
	// proportion of cgrow that goes to reproduction
	const double repr_prop = 0.1;
	// minimum age for reproduction in years.
	// TODO: fix to be a function of tree size and not age, and a PFT parameter.
	// also, check that age is counted in years and not days.
	const double min_age_for_repr = 10;
	// initialise reproduction limit (optimum)
	double crepr_opt = 0.0;

	double csap = indiv.cmass_sap;
	// initialise leaf and root growth limits (optimum)
	double cleaf_opt = csap * fns;
	double croot_opt = indiv.ltor>0.0?cleaf_opt / indiv.ltor:0.0;
	// initialise reserve limit (optimum)
	double creserve_opt = 0.0;

	// initialise sapwood optimum. Here start with last calculated csap limit
	double csap_opt = indiv.cmass_sap_limit;

	if (indiv.leaf_growth) {
		// at the start of the leaf growth period, estimate a new cmass_sap_limit 
		// to allow some sapwood growth
		if (negligible(indiv.cmass_sap_limit)) { 
			// ensures it is only set ones at the start of leaf growth
			// at least want to replace what will be turned over in the year
			double ratio = 1.0 + indiv.pft.turnover_sap;
			csap_opt = ratio * indiv.cmass_sap;
		}
	} else {
		// this ensures that in the next year's leaf_growth period the limit  
		// will be set again
		csap_opt = 0.0;
	}
	
	if (indiv.leaf_full_switch_summer || indiv.senescence_active) {
		// we keep the latest calculated leaf and root limits in the leaf full 
		// on after the summer solstice and also during the active senescence phase
		cleaf_opt = indiv.cmass_leaf_limit;
		croot_opt = indiv.cmass_root_limit;
	}
	
	if (indiv.leaf_senesced || indiv.senescence_active) {
		// prevent allocation to leaf and root in leaves off season
		// FIXME: need to handle the disabling of allocation differently,
		// e.g. indiv.nsc.leaf_alloc_enabled = [true|false], root_ sap_ ...
		cleaf_opt = 0.0;
		croot_opt = 0.0;
	}

	// only allow reserve filling when we have full leaf out after summer to 
	// winter transition and in senesced period
	if (indiv.leaf_senesced || indiv.leaf_full_switch_summer || indiv.senescence_active) {
		// optimal carbon storage. either to replace the amount of leaf 
		// and root turned over within a year or the maximum capacity.
		// TODO: check if  cstorage_max has enough to turn over x3 year-worth 
		// of leaves and roots.
		creserve_opt = indiv.pft.turnover_leaf * cleaf_opt + indiv.pft.turnover_root * croot_opt;
		creserve_opt = min(creserve_opt, indiv.cstorage_max);
		// this will disable reserve allocation in all other periods (e.g. during leave flushout)
	} 
	
	if (indiv.leaf_full_switch_summer || indiv.leaf_out) {
		// assuming that the reproduction occurs only in spring and summer
		// TODO: check if assumption is true.
		if (indiv.age >= min_age_for_repr) {
			// assuming that repr_prop% of the daily cgrow will potentialyly go  
			// to reproduction in case 
			crepr_opt = indiv.cmass_repr + (indiv.cgrow * repr_prop);
			// else: crepr_opt stays 0
		}
		// else: crepr_opt stays 0
	}

	// outputs:
	indiv.cmass_leaf_limit = cleaf_opt;
	indiv.cmass_root_limit = croot_opt;
	indiv.cmass_sap_limit = csap_opt;
	indiv.cmass_repr_limit = crepr_opt;
	indiv.cmass_storage_limit = creserve_opt;
}

// calculate the alpha factors, which determines the proportion of NPP that will 
// be allocated to the different C pools, depening on their phenology and the
// sink strength. 
void decide_on_allocation(Individual &indiv) {
	// Inputs: cmass_<pool>, cmass_<pool>_limit, phenology  
	// Outputs: alpha_<pool>

	// unless in the senescence period (autumn), calculate delta for leaf (Grote)
	double d_leaf = 0;
	if (indiv.cmass_leaf_limit > 0 && !indiv.senescence_active) {
		d_leaf = max(0.0, (indiv.cmass_leaf_limit - indiv.cmass_leaf) / indiv.cmass_leaf_limit);
	}
	// unless in the senescence period (autumn), calculate delta for root (Grote)
	double d_root = 0;
	if (indiv.cmass_root_limit > 0 && !indiv.senescence_active) {
		d_root = max(0.0, (indiv.cmass_root_limit - indiv.cmass_root) / indiv.cmass_root_limit);
	}
	// calculate delta for sapwood (Grote)
	double d_sap = 0.0;
	if (indiv.cmass_sap_limit > 0) {
		d_sap = divide(indiv.cmass_sap_limit - indiv.cmass_sap, indiv.cmass_sap_limit);
	}
	// calculate delta for storage pool (Grote)	
	double d_reserve = 0;
	if (indiv.cmass_storage_limit > 0) {
	 	d_reserve = divide(indiv.cmass_storage_limit - indiv.cmass_storage , indiv.cmass_storage_limit);
	}
	// calculate delta for reproduction
	double d_repr = 0;
	if (indiv.cmass_repr_limit > 0) {
		d_repr = divide(indiv.cmass_repr_limit - indiv.cmass_repr ,indiv.cmass_repr_limit);
	}
	// TODO: decide how to include reproduction in this scheme.
	d_repr = 0.0; 

	// sum all the pool sink strngths.
	double d_summed = d_leaf + d_root + d_sap + d_reserve + d_repr;

	// if nobody wants anything (ie, in case C pools are close to their limit,   
	// and/or during senescence phase):
	if (is_zero(d_summed)) { 	
		// since cstorage seems to increase to hugh values, lets dump things into stem
		// which will trigger higher optimums in other c-pools
		// this allows exess carbon to be sent to sapwood (and stem growth)
		d_sap = 1.0;
		d_summed = d_leaf + d_root + d_sap + d_reserve + d_repr;
	}
	// outputs:
	indiv.cgrow_leaf = indiv.cgrow * (d_leaf / d_summed);
	indiv.cgrow_root = indiv.cgrow * (d_root / d_summed);
	indiv.cgrow_sap = indiv.cgrow * (d_sap / d_summed);
	// not necessary when using Juergen's update pools
	// // indiv.cgrow_reserve = indiv.cgrow * (d_reserve / d_summed); 
	indiv.cgrow_repr = indiv.cgrow * (d_repr / d_summed);
	
}

// returns how much available storage there is (kgC/m2)
double available_storage(Individual &indiv) {

	// Fraction of storage capacity unavailable for growth.
	const double cstorage_reserve = 0.1;

	// Maximum daily mobilisation rate of storage as a fraction of storage
	// capacity.
	const double max_mobil_rate = 0.05;

	// maximum amount of storage that can be extracted today
	double max_mobil = max_mobil_rate * indiv.cstorage_max;
	// the amount of storage that has to be held in the reserve
	double reserve = cstorage_reserve * indiv.cstorage_max;
	// in case there is enough storage to spare today, take the maximum that you can
	if ((indiv.cmass_storage - max_mobil) > reserve) {
		return max_mobil;
	}
	// in case there isn't enough storage to spare today, take the maximum possible 
	// without overdawing.
	return max(0.0, indiv.cmass_storage - reserve);

}

// calculate allowance (cgrow) depending on cmass_storage
double calculate_cgrow_tree(Individual &indiv) {

	Climate climate = indiv.vegetation.patch.get_climate();


	// parameters of "saver" responce
	// TODO: calibrate to express differences between "saver" and "spender" responces.
	// On top of this (otr instead), add temperature and soil moisture limitation to NSC vs growth
	// (Fatichi et al, 2014), which limits growth in certain conditions (whilst saving NSCs),
	// such as in low temperature and low soil moisture (low soil water potential)
	double const shape_pie = 16.0;
	double const inflection_pie_saver = 0.7; // Saver
	double const inflection_pie_int = 0.5; // intermediate
	double const inflection_pie_spender = 0.2; // Spender
	// maximum daily mobilisation rate of c in the storage. 
	const double max_mobil_cstorage_tree = 0.1;
	// fraction of cmass_storage that is not available to growth
	const double frac_cstorage_res_tree = 0.15;
	// actual amount that is not available to growth (kgC/m^2)
	double cstorage_res = frac_cstorage_res_tree * indiv.cstorage_max;
	// amount of c available for growth today (kgC/m^2)
	double cmass_avail = max(0.0, (indiv.cmass_storage - cstorage_res));
	
	double store_state = 0;
	// sigmoid function that will change depending on the saver-spender spectrum of the pft.
	double piesize_coeff = 0;
	if (indiv.cstorage_max > 0) {
		
		store_state = (indiv.cmass_storage / indiv.cstorage_max);
		
		piesize_coeff = 1 - ((1 / (1 + exp(shape_pie * (store_state - inflection_pie_int)))));
	}
	
	double cgrow = piesize_coeff * max_mobil_cstorage_tree * cmass_avail;
	
	// temperature scalar for decoupling potosynthesis from growth (see Figure 9 in  Korner's 2-20 book)
	double ftemp = temp_growth(indiv.patchpft().airtemp_10.mean(), 
			indiv.pft.Ha_growth, indiv.pft.T0_growth, indiv.pft.A_growth, 
			indiv.pft.temp_dormancy);
	assert(ftemp <= 1);
	assert(ftemp >= 0);
	// water scalar for decoupling potosynthesis from growth (Fatichi et al 2014)
	double fwater = 1;
	// daylength scalar (currently off)
	double fphoto = 1;

	// rescale cgrow aacording to environmental conditions (temperature, soil water and light)
	cgrow *= ftemp * fwater * fphoto;
	// actual allowance of c for growth today (kgC/m^2)
	return piesize_coeff * max_mobil_cstorage_tree * cmass_avail;

}

// Determine how much of the available C (NPP + C fraction from storage) goes to
// growth vs. storage. Equivalent to calc_fc_growth() in the grote code
void storage_growth_partitioning_tree(Individual &indiv) {
	// Inputs: indiv.dnpp
	// Outputs: indiv.cmass_storage, indiv.cgrow (potential growth in kgC/m2)
	// Optional outputs which the grass code use which we may also set:
	// indiv.growth_demand, indiv.cexcess, indiv.nexcess, indiv.nmass_missing

	// Set indiv.cgrow, indiv.cexcess, and indiv.nexcess to 0.
	indiv.cgrow = 0.0;
	indiv.cexcess = 0.0;
	indiv.nexcess = 0.0;
	indiv.nmass_missing = 0.0;
	indiv.ctolitter = 0;
	indiv.ntolitter = 0;

	// first add NPP to storage.
	indiv.cmass_storage += indiv.dnpp;
	
	// if this causes storage to be negitive - kill tree and return.
	if (indiv.cmass_storage < 0) {
		// Not enough C in storage to fulfill respiration demands --> kill individual
		double cmass_deficit = -indiv.cmass_storage;
		indiv.cmass_storage = 0.0;

		// The respiration and NPP fluxes reported in canexch are
		// now incorrect, so we need to fix them. Note that the RA
		// flux is not used for mass balance checks (it's only used
		// for outputs). The NPP flux is used in mass balance checks
		// however, so it's important to correct it here.
		indiv.report_flux(Fluxes::RA, -cmass_deficit);
			// (npp has opposite sign to respiration)
		indiv.report_flux(Fluxes::NPP, cmass_deficit);

		indiv.kill();  // flushes all pools to litter
		indiv.vegetation.killobj();  // removes individual object
							
		//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
		// Don't think we need to call nextobj() here, as we've just
		// removed an element from the list; when this happens, the
		// "current" pointer internal to the list object points to
		// the element that was previously the "next" element.
		return;
	}

	// calculate allowance (cgrow) depending on the enviro conditions and a  
	// maximumvalue (eg, 10% of storage).
	// Note: cgrow can be > 0 even if storage has gon wodn (ie dnpp < 0).
	indiv.cgrow = calculate_cgrow_tree(indiv);

	assert(indiv.cgrow < indiv.cmass_storage);
	assert(positive(indiv.cgrow));
	
	// remove cgrow from storage
	indiv.cmass_storage -= indiv.cgrow;
}

// Determine how much of the available C (NPP + C fraction from storage) goes to
// growth vs. storage. Equivalent to calc_fc_growth() in the grote code
void storage_growth_partitioning_tree_grote(Individual &indiv) {
	// Inputs: indiv.dnpp
	// Outputs: indiv.cmass_storage, indiv.cgrow (potential growth in kgC/m2)
	// Optional outputs which the grass code use which we may also set:
	// indiv.growth_demand, indiv.cexcess, indiv.nexcess, indiv.nmass_missing
	
	if (indiv.dnpp < 0) {

		if ((indiv.cmass_storage + indiv.dnpp) > 0) {
			// deficit can be satisfied by the storage
			indiv.cmass_storage += indiv.dnpp; // it's a + becuase npp is negative
			// so that growth will be entirely from storage (see below).
			indiv.dnpp = 0;
		} else {
			// Not enough C in storage to fulfill respiration demands --> kill individual
			double cmass_deficit = -(indiv.cmass_storage + indiv.dnpp);
			indiv.cmass_storage = 0.0;

			// The respiration and NPP fluxes reported in canexch are
			// now incorrect, so we need to fix them. Note that the RA
			// flux is not used for mass balance checks (it's only used
			// for outputs). The NPP flux is used in mass balance checks
			// however, so it's important to correct it here.
			indiv.report_flux(Fluxes::RA, -cmass_deficit);
				// (npp has opposite sign to respiration)
			indiv.report_flux(Fluxes::NPP, cmass_deficit);

			indiv.kill();  // flushes all pools to litter
			indiv.vegetation.killobj();  // removes individual object
								
			//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
			// Don't think we need to call nextobj() here, as we've just
			// removed an element from the list; when this happens, the
			// "current" pointer internal to the list object points to
			// the element that was previously the "next" element.
			return;
		}

	
	}
	// bminc is set to today's npp.
	double bminc = indiv.dnpp;

	// Length (in days) of leaf flushout period.
	// TODO: this is a PFT parameter in Grote. Perhaps it should be in dave too.
	const int flushout_ndays = 45;

	// Are we in the leaf flushout period?
	// TODO: decide on phenology
	bool leaf_flushout = indiv.phen > 0;

	if (leaf_flushout) {
		assert(flushout_ndays > 0);
		// TODO: create function to set cmass optimum for pools.
		double leaf_turnover = indiv.cmass_leaf_limit * indiv.pft.turnover_leaf;
		double root_turnover = indiv.cmass_root_limit * indiv.pft.turnover_root;
		// TODO: currently the maximum supply is the daily turnover of leaves 
		// and roots. Possibly add some more from the storage (if available) with 
		// a maximum daily fraction that the tree can get from storage. In Juergen's 
		// grass it is 5% of storage per day.
		double fc_growth_suggested = (leaf_turnover + root_turnover) / flushout_ndays;
		fc_growth_suggested = max(0.0, fc_growth_suggested);

		if (fc_growth_suggested > indiv.dnpp) {
			// We want to grow more than can be sustained by today's NPP.
			// Therefore we draw from stroage.
	
			double from_storage = available_storage(indiv);
			if ((indiv.dnpp + from_storage) > fc_growth_suggested) {
				// demand limited: we can meet growth demand by drawing from 
				// storage
				from_storage = fc_growth_suggested - indiv.dnpp;
			} else {
				// supply limited: we can not meet growth demand by drawing 
				// from storage.
				fc_growth_suggested = from_storage + indiv.dnpp;
			}
			// take from storage
			indiv.cmass_storage -= from_storage;
		} else {
			// Today's NPP is greater than our potential growth. Therefore we
			// send some to storage.
			// TODO: according to the lpj-grote code, allocation to storage 
			// depends on phenology (ie after summer solstice).  
			double cstore = indiv.dnpp - fc_growth_suggested;
			indiv.cmass_storage += cstore;
			bminc = fc_growth_suggested;
		}
	}

	indiv.cgrow = bminc;
}

// Determine how much of the available C (NPP + C fraction from storage) goes to growth vs. storage.
// Function calculates cgrow (daily growth), cstore (change in cmass_storage) 
// and updates indiv.cmass_storage. 
// Note: cgrow and cstore are fluxes (kgC/m2/d), cmass_storage is a pool (kg/m2).
void storage_growth_partitioning_dave(Vegetation& vegetation) {

	// Loop through individuals
	// Note - we sometimes kill individuals inside this loop, but vegetation.killobj()
	// requires the internal pointer to be set appropriately, which effectively
	// prevents us from using a range-based loop, and instead forces us to use
	// the clunky mess you see below. todo: refactor the container classes!
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			vegetation.nextobj();
			continue;
		}
		size_t nobj = vegetation.nobj;
		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			storage_growth_partitioning_grass(indiv);
		} else {
			storage_growth_partitioning_tree(indiv);
		}

		// If this individual was killed, vegetation.nobj will go down by one.
		// In this case, getobj() will return the individual which was
		// previously after the individual which was killed and removed from the
		// vegetation, and so we don't want to call nextobj() if that happened
		// (if we do, we'd be skipping the next individual).
		//
		// This is really janky, so in the long term we should think about
		// restructuring this, especially if we later decide that other
		// functions should be able to kill an individual.
		if (vegetation.nobj == nobj) {
			// indiv was not killed
			vegetation.nextobj();
		}
	}
}

void carbon_allocation_grass(Individual &indiv) {
	double calloc_leaf  = 0.0;
	double calloc_root  = 0.0;
	double calloc_crown = 0.0;
	double calloc_repr  = 0.0;

	Patchpft &ppft = indiv.patchpft();

	if (!negligible(indiv.cgrow)) {
		
		// 1) Crowns
		if (indiv.age < 1.0 && indiv.pft.lifeform != GRASS_ANNUAL) {
			calloc_crown = 2.0 * indiv.pft.crownfrac;
		} else if (indiv.pft.lifeform == GRASS_ANNUAL && indiv.growth_stage == reproductive) {
			calloc_crown = 0.0;
		} else {
			calloc_crown = indiv.pft.crownfrac;
		}
		
		// 2) Reproduction
		if (indiv.growth_stage == reproductive) {
			if (indiv.age < 1.0 && indiv.pft.lifeform != GRASS_ANNUAL) {
				calloc_repr = min(0.2 * indiv.pft.reprfrac,(1.0 - calloc_crown));
			} else {
				calloc_repr = min(indiv.pft.reprfrac,(1.0 - calloc_crown));
			}
		} else {
			calloc_repr = 0.0;
		}

		// 3) Shoots (Leaves + culms)
		// express ltor_max as maximum fractional allocation to shoots (leaves + culms)
		double calloc_leaf_max = 1.0 - 1.0 / (1.0 + indiv.pft.ltor_max);

		// calculate fraction of C allocated to shoots (leaves + culms) following Johnson 2016 Eq. 3.80
		calloc_leaf = calloc_leaf_max * pow(indiv.nscal_5.mean() * ppft.wscal_5.mean(), 0.5);
		calloc_leaf *= (1.0 - calloc_repr - calloc_crown);
		
		// 4) Roots
		calloc_root = 1.0 - calloc_leaf - calloc_crown - calloc_repr;

	}
	else {
		calloc_leaf  = 0.0;
		calloc_root  = 0.0;
		calloc_crown = 1.0;
		calloc_repr  = 0.0;
	}

	assert(equal(1.0, calloc_leaf + calloc_root + calloc_crown + calloc_repr));

	indiv.cgrow_leaf  = calloc_leaf * indiv.cgrow;
	indiv.cgrow_root  = calloc_root * indiv.cgrow;
	indiv.cgrow_crown = calloc_crown * indiv.cgrow;
	indiv.cgrow_repr  = calloc_repr * indiv.cgrow;
}

// Determine the sink strength of each carbon pool in the individual.
// This probably encapsulates adjust_cpool_optimum() and decide_on_allocation()
// in the grote code.
void carbon_allocation_tree(Individual &indiv) {
	// Inputs: cmass_<pools>, growing conditions maybe?
	// Outputs: indiv.cgrow_leaf, indiv.cgrow_root, indiv.cgrow_sap, indiv.cgrow_repr

	// Steps: 
	// 1) calculate the target growth (a-la Grote) for each pool.
	// 2) calculate the alpha factor (decide_on_allocation) using sink strength 
	//    (ie the difference between the actual pool size and optimal growth)

	Pft &pft = indiv.pft;
	// get the patchpft object for this individual's pft
	Patchpft &ppft = indiv.patchpft();	
	// 1) calculate the target growth (a-la Grote) for each pool.
	// calculate leaf limit (optimum)
	const double MAX_HEIGHT = 150;
	const double MIN_HEIGHT = 0.05;
	// limit height to below max and above min (to avoid FNS > 1)(from David's 
	// code). 
	// TODO: what about seedlings????
	double height = min(max(MIN_HEIGHT, indiv.height), MAX_HEIGHT);
	
	// proportion of cgrow that goes to reproduction. 
	// TODO: find literature values
	const double repr_prop = 0.1;
	// minimum age for reproduction in years.
	// TODO: fix to be a function of tree size and not age, and a PFT parameter.
	// also, check that age is counted in years and not days.
	const double min_age_for_repr = 10;
	// maximum proportion of growth that will be allocated to reproduction
	// TODO: check this out. guestimated.
	const double max_alpha_repr_factor = 0.2;

	// cmass sapwood target (to grow at least the turnover of today) (kgC/m^2)
	indiv.cmass_sap_limit = indiv.cmass_sap * (1 + (pft.turnover_sap));
	// cmass leaf target (to grow to supply leaves to the current sapwood ) (kgC/m^2)
	indiv.cmass_leaf_limit = max(indiv.cmass_sap * (pft.k_latosa / (indiv.sla() * pft.wooddens * height)) , indiv.cmass_leaf);
	//  calculate the ltor according to the leaf optimum and belogrownd stress
	// TODO: think what to do when nscal or wscal is zero. Now set a minimum ltor
	// ...to .5 of ltor_max

	indiv.ltor = max(min(ppft.wscal_5.mean(), indiv.nscal_5.mean()) * pft.ltor_max, 0.5 * pft.ltor_max);
	// cmass root target according to the leaf optimum and ltor (kgC/m^2)
	indiv.cmass_root_limit = 0.0;
	if (indiv.ltor > 0.0) {
		indiv.cmass_root_limit = max(indiv.cmass_leaf_limit / indiv.ltor, indiv.cmass_root);
	}
	// cmass reproduction target
	indiv.cmass_repr_limit = 0.0;
	if (indiv.age > min_age_for_repr) {
		indiv.cmass_repr_limit = repr_prop * indiv.cgrow + indiv.cmass_repr;
	}

	// Calculate the alpha factor (what proportion to allocate for each pool) according
	// to the distance of actual C to it's target.
	double relative_alpha_leaf = 0;
	if (!is_zero(indiv.cmass_leaf_limit)) {
		relative_alpha_leaf = divide(indiv.cmass_leaf_limit - indiv.cmass_leaf, indiv.cmass_leaf_limit);
		// relative_alpha_leaf = max(0.05, indiv.fpar);
	}

	double relative_alpha_root = 0;
	if (!is_zero(indiv.cmass_root_limit)) {
		relative_alpha_root = divide(indiv.cmass_root_limit - indiv.cmass_root, indiv.cmass_root_limit);
		// relative_alpha_root = max(0.05, 1 - min(ppft.wscal_5.mean(), indiv.nscal_5.mean()));
	}
	
	double relative_alpha_sap =  0.0;
	if (!is_zero(indiv.cmass_sap_limit)) {
		relative_alpha_sap = divide(indiv.cmass_sap_limit - indiv.cmass_sap, indiv.cmass_sap_limit);
		// relative_alpha_sap = max(0.05, 1 - indiv.fpar);
	}

	// To prevent allocation to reproduction if the tree is stressed (lower priority.),
	// 1. calculate the mean of the alphas or leaf, root and sapwood, that will 
	// be used ased as a threshold for allocating to reproduction. 
	// double mean_relative_alphas = (relative_alpha_leaf + relative_alpha_root + relative_alpha_sap + relative_alpha_storage) / 4.0;
	double mean_relative_alphas = (relative_alpha_leaf + relative_alpha_root + relative_alpha_sap) / 3.0;
	// logistic function of the relative strength of  allocation to reproduction
	// depending on how far the other pools are from their target.
	// TODO: Unsupported assumption: logistic function assumed and it's parameters 
	// guestimated.  Change that.

	double alpha_repr_factor = (max_alpha_repr_factor / 
                       (1 + exp(20 * (mean_relative_alphas - .25))));
	// calculate alpha reproduction "tunned down" by the factor that depends on
	// the conditions of the other pools.
	double relative_alpha_repr = 0.0;
	if (indiv.cmass_repr_limit > 0) {
		relative_alpha_repr = alpha_repr_factor * 
			divide(indiv.cmass_repr_limit - indiv.cmass_repr, indiv.cmass_repr_limit);
	}
	// caps the relative_alpha_repr to 0.1 to avoid over reprentation of cgrow_repr
	// over the rest of the pools
	// TODO: Fix this (0.1???)
	relative_alpha_repr = min(relative_alpha_repr, 0.1);
	//  sum of all the relative alphas so it will be used to calculate the alphas
	double sum_relative_alphas = relative_alpha_leaf + relative_alpha_root + 
						relative_alpha_sap + relative_alpha_repr;								

	// calculate the alpha for each pool relative to the alpha of the rest of the pools.
	if (is_zero(sum_relative_alphas)) {
		indiv.alpha_leaf = 0;
		indiv.alpha_root = 0;
		indiv.alpha_sap  = 0;
		indiv.alpha_repr = 0;
	} else {
		// No need for protected divde here. If we refactor out the above check
		// for a zero value in the denominator, we should reinstate the use of a
		// protected division here.
		indiv.alpha_leaf = relative_alpha_leaf / sum_relative_alphas;
		indiv.alpha_root = relative_alpha_root / sum_relative_alphas;
		indiv.alpha_sap = relative_alpha_sap / sum_relative_alphas;
		indiv.alpha_repr = relative_alpha_repr / sum_relative_alphas;
	}

	// outputs:
	indiv.cgrow_leaf = indiv.cgrow * indiv.alpha_leaf;
	indiv.cgrow_root = indiv.cgrow * indiv.alpha_root;
	indiv.cgrow_sap = indiv.cgrow * indiv.alpha_sap;
	// not necessary when using Juergen's update pools
	// // indiv.cgrow_reserve = indiv.cgrow * (d_reserve / d_summed); 
	indiv.cgrow_repr = indiv.cgrow * indiv.alpha_repr;

	assert(positive(indiv.alpha_leaf));
	assert(positive(indiv.alpha_root));
	assert(positive(indiv.alpha_sap));
	assert(positive(indiv.alpha_repr));
	assert(equal(1.0, indiv.alpha_leaf + indiv.alpha_root + indiv.alpha_sap + indiv.alpha_repr));
}


// Determine the sink strength of each carbon pool in the individual.
// This probably encapsulates adjust_cpool_optimum() and decide_on_allocation()
// in the grote code.
void carbon_allocation_tree_grote(Individual &indiv) {
	// Inputs: Phenology, growing conditions maybe?
	// Outputs: indiv.cgrow_leaf, indiv.cgrow_root, indiv.cgrow_sap, indiv.cgrow_repr

	// Steps: 
	// 1) calculate the optimal growth (a-la Grote) for each pool.
	// 2) calculate the alpha factor (decide_on_allocation) using sink strength 
	//    (ie the difference between the actual pool size and optimal growth)
	
	// 1) calculate the optimal growth (a-la Grote) for each pool.
	adjust_cpool_optimum(indiv);
	// 2) calculate the alpha factor (decide_on_allocation) using sink strength 
	//    (ie the difference between the actual pool size and optimal growth)
	decide_on_allocation(indiv);

}

// Determine C allocation fractions to shoots, crowns, roots, and reproductive tissue in grasses.
// Note that leaf pool represents the shoot and includes leaves and culms.
void carbon_allocation_dave(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}

		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			carbon_allocation_grass(indiv);
		} else {
			carbon_allocation_tree(indiv);
		}
	}
}

/*
Get the C:N ratio for new N growth which would maintain the current C:N ratio.
If the current C:N ratio is already outside the allowed C:N bounds, calculate
a C:N ratio for new growth which would put the C:N ratio back within the legal
upper/lower bounds.

@param cmass: Current mass of the C pool (kgC/m2).
@param nmass: Current mass of the N pool (kgC/m2).
@param cton_min: Minimum allowed C:N ratio.
@param cton_max: Maximum allowed C:N ratio.
@param cgrow: Today's C growth (kgC/m2).

Note: technically the units don't really matter, as long as they're consistent.
*/
double constrain_cton(double cmass, double nmass, double cton_min, double cton_max, double cgrow) {
	// Can happen if initialised to 0
	if (nmass == 0 || cgrow == 0) {
		return (cton_min + cton_max) / 2.0;
	}

	// Get C:N ratio of existing biomass. If it's within the allowed bounds then
	// there's nothing more to do.

	double cton = cmass / nmass;

	if (greater_than_or_equal(cton, cton_max) || less_than_or_equal(cton, cton_min)) {
		// Current C:N ratio is outside the allowed C:N ratio. We need to
		// up- or down-regulate the N demand in order to move towards the
		// nearest C:N ratio bounds.
		double target_cton = cton > cton_max ? cton_max : cton_min;

		// Get the N content which would put us at the target C:N ratio after
		// today's growth.
		double c_new = cmass + cgrow;
		double n_new = c_new / target_cton;

		// Get the delta between current and target N content.
		double ngrow = max(0.0, n_new - nmass);

		// Calculate the C:N ratio of new growth which will put us at the target
		// C:N ratio.
		cton = cgrow / ngrow;
	} else if (cton < cton_min) {
		// Current C:N ratio is below the minimum C:N ratio. In this case, we
		// simply request N at the minimum C:N ratio. This will move us toward
		// the minimum C:N ratio, though it won't get us all the way there.
		cton = cton_min;
	}

	assert(positive(cton));
	return cton;
}

void ndemand_grass(Individual &indiv) {
	const double frac_nchange = 0.01;   // TODO: better to express as C:N? Zaehle as reference?
	const double max_update_growth = 0.1;

	// update growth depending on N stress of plant
	// TODO: double check cases where nmass_storage > nstorage_max!!
	// TODO: constrain cton ratio ;) 
	double nstorage_fill = min(0.999999, max(0.000001, (indiv.nmass_storage / indiv.nstorage_max)));

	// Local	
	double cton_demand_leaf = 0.0;
	double cton_demand_root = 0.0;
	double cton_demand_crown = 0.0;
	double cton_demand_repr = 0.0;

	double ngrow_demand_leaf  = 0.0;
	double ngrow_demand_root  = 0.0;
	double ngrow_demand_repr  = 0.0;
	double ngrow_demand_crown = 0.0;

	// N demand of new growth
	if (!negligible(indiv.cgrow)) {
		// Current C:N ratios
		double cton_leaf = constrain_cton(
			indiv.get_cmass_leaf(), indiv.get_nmass_leaf(),
			indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
			indiv.cgrow_leaf);
		double cton_root = constrain_cton(
			indiv.cmass_root, indiv.nmass_root,
			indiv.pft.cton_root_min, indiv.pft.cton_root_max,
			indiv.cgrow_root);
		double cton_crown = constrain_cton(
			indiv.cmass_crown, indiv.nmass_crown,
			indiv.pft.cton_root_min, indiv.pft.cton_root_max,
			indiv.cgrow_crown);
		double cton_repr = constrain_cton(
			indiv.cmass_repr, indiv.nmass_repr,
			indiv.pft.cton_repr_min, indiv.pft.cton_repr_max,
			indiv.cgrow_repr);
		
		cton_demand_leaf = (1.0 - nstorage_fill) * (indiv.pft.cton_leaf_max - indiv.pft.cton_leaf_min) + indiv.pft.cton_leaf_min;
		cton_demand_root = (1.0 - nstorage_fill) * (indiv.pft.cton_root_max - indiv.pft.cton_root_min) + indiv.pft.cton_root_min;
		cton_demand_crown = (1.0 - nstorage_fill) * (indiv.pft.cton_root_max - indiv.pft.cton_root_min) + indiv.pft.cton_root_min;
		cton_demand_repr = (1.0 - nstorage_fill) * (indiv.pft.cton_repr_max - indiv.pft.cton_repr_min) + indiv.pft.cton_repr_min;

		ngrow_demand_leaf  = indiv.cgrow_leaf / cton_demand_leaf;
		ngrow_demand_root  = indiv.cgrow_root / cton_demand_root;
		ngrow_demand_crown = indiv.cgrow_crown / cton_demand_crown;
		ngrow_demand_repr  = indiv.cgrow_repr / cton_demand_repr;
	} else {
		ngrow_demand_leaf  = 0.0;
		ngrow_demand_root  = 0.0;
		ngrow_demand_crown = 0.0;
		ngrow_demand_repr  = 0.0;
	}

	// N demand of existing biomass (weight by current activity (growth) of the tissue?)
	double nupdate_demand_leaf  = max(0.0, indiv.get_cmass_leaf() / indiv.pft.cton_leaf_min - indiv.get_nmass_leaf()) * frac_nchange;
	double nupdate_demand_root  = max(0.0, indiv.cmass_root / indiv.pft.cton_root_min - indiv.nmass_root) * frac_nchange;
	double nupdate_demand_crown = max(0.0, indiv.cmass_crown / indiv.pft.cton_root_min - indiv.nmass_crown) * frac_nchange;
	double nupdate_demand_repr  = max(0.0, indiv.cmass_repr / indiv.pft.cton_repr_min - indiv.nmass_repr) * frac_nchange;

	// Make sure that nupdate_demands do not exceed N demands for growth (i.e. prioritise growth over existing biomass)
	// Otherwise plants could not grow if C:N of the existing biomass is high.
	// Alternatively, one could ask for the min C:N for new growth and not update existing biomass.
	double ngrow_demand_growth_total = ngrow_demand_leaf + ngrow_demand_root + ngrow_demand_crown + ngrow_demand_repr;
	double ngrow_demand_update_total = max(0.0, nupdate_demand_leaf + nupdate_demand_root + nupdate_demand_crown + nupdate_demand_repr);
	double growth_update_ratio = min(max_update_growth, ngrow_demand_growth_total / ngrow_demand_update_total);

	// Take N demand of storage into account
	indiv.storendemand = max(0.0, indiv.nstorage_max - indiv.nmass_storage);

	// Total biomass N demand of individual
	indiv.leafndemand  = ngrow_demand_leaf + nupdate_demand_leaf * growth_update_ratio;
	indiv.rootndemand  = ngrow_demand_root + nupdate_demand_root * growth_update_ratio;
	indiv.crownndemand = ngrow_demand_crown + nupdate_demand_crown * growth_update_ratio;
	indiv.reprndemand  = ngrow_demand_repr + nupdate_demand_repr * growth_update_ratio;
	indiv.ndemand      = indiv.leafndemand + indiv.rootndemand + indiv.crownndemand + indiv.reprndemand + indiv.storendemand;

}

// Calculate Nitrogen demand for a tree.
void ndemand_tree(Individual &indiv) {
	// Inputs: Whatever is required. Grasses use cgrow and attempt to maintain a
	// certain c:n ratio
	// Ouptuts: indiv.leafndemand, indiv.rootndemand, indiv.sapndemand,
	// indiv.reprndemand, all of which are in kgN/m2
	const double frac_nchange = 0.01;   // TODO: better to express as C:N? Zaehle as reference?
	const double max_update_growth = 0.1;

	// Local	
	double ngrow_demand_leaf  = 0.0;
	double ngrow_demand_root  = 0.0;
	double ngrow_demand_repr  = 0.0;
	double ngrow_demand_sap   = 0.0;
	

	// N demand of new growth
	if (!negligible(indiv.cgrow)) {
		// Current C:N ratios
		double cton_leaf = constrain_cton(
			indiv.get_cmass_leaf(), indiv.get_nmass_leaf(),
			indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
			indiv.cgrow_leaf);
		double cton_root = constrain_cton(
			indiv.cmass_root, indiv.nmass_root,
			indiv.pft.cton_root_min, indiv.pft.cton_root_max,
			indiv.cgrow_root);
		double cton_sap = constrain_cton(
			indiv.cmass_sap, indiv.nmass_sap,
			indiv.pft.cton_sap_min, indiv.pft.cton_sap_max,
			indiv.cgrow_sap);
		double cton_repr = constrain_cton(
			indiv.cmass_repr, indiv.nmass_repr,
			indiv.pft.cton_repr_min, indiv.pft.cton_repr_max,
			indiv.cgrow_repr);
		ngrow_demand_leaf  = indiv.cgrow_leaf / cton_leaf;
		ngrow_demand_root  = indiv.cgrow_root / cton_root;
		ngrow_demand_sap   = indiv.cgrow_sap / cton_sap;
		ngrow_demand_repr  = indiv.cgrow_repr / cton_repr;
	} else {
		ngrow_demand_leaf  = 0.0;
		ngrow_demand_root  = 0.0;
		ngrow_demand_sap   = 0.0;
		ngrow_demand_repr  = 0.0;
	}

	// N demand of existing biomass (weight by current activity (growth) of the tissue?)
	// TODO: this implies that growth is limited if current C:N is high, which shouldn't be the case. Make growth a priority?
	double nupdate_demand_leaf  = max(0.0, indiv.get_cmass_leaf() / indiv.pft.cton_leaf_min - indiv.get_nmass_leaf()) * frac_nchange;
	double nupdate_demand_root  = max(0.0, indiv.cmass_root / indiv.pft.cton_root_min - indiv.nmass_root) * frac_nchange;
	double nupdate_demand_sap   = max(0.0, indiv.cmass_sap / indiv.pft.cton_sap_min - indiv.nmass_sap) * frac_nchange;
	double nupdate_demand_repr  = max(0.0, indiv.cmass_repr / indiv.pft.cton_repr_min - indiv.nmass_repr) * frac_nchange;

	// Make sure that nupdate_demands do not exceed N demands for growth (i.e. prioritise growth over existing biomass)
	// Alternatively, one could ask for the min C:N for new growth and not update existing biomass.
	double ngrow_demand_growth_total = ngrow_demand_leaf + ngrow_demand_root + ngrow_demand_sap + ngrow_demand_repr;
	double ngrow_demand_update_total = max(0.0, nupdate_demand_leaf + nupdate_demand_root + nupdate_demand_sap + nupdate_demand_repr);
	double growth_update_ratio = min(max_update_growth, ngrow_demand_growth_total / ngrow_demand_update_total);

	// Total biomass N demand of individual
	indiv.leafndemand  = ngrow_demand_leaf + nupdate_demand_leaf * growth_update_ratio;
	indiv.rootndemand  = ngrow_demand_root + nupdate_demand_root * growth_update_ratio;
	indiv.sapndemand = ngrow_demand_sap + nupdate_demand_sap * growth_update_ratio;
	indiv.reprndemand  = ngrow_demand_repr + nupdate_demand_repr * growth_update_ratio;
	indiv.ndemand      = indiv.leafndemand + indiv.rootndemand + indiv.sapndemand + indiv.reprndemand;

}

// Calculate N demand based on growth and N status of the plant
void ndemand_dave(Vegetation& vegetation) {
	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}
		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			ndemand_grass(indiv);
		} else {
			ndemand_tree(indiv);
		}

		assert(positive(indiv.ndemand));
	}
}

// Calculate N uptake which is limited by competition between individuals and
// availability of N in the soil.
void nuptake_dave(Vegetation& vegetation) {
	Patch& patch = vegetation.patch;
	Soil& soil = patch.soil;
	Gridcell& gridcell = patch.stand.get_gridcell();

	// Create vector of CompetingIndividual objects, which will be passed into
	// ncompete(). These are populated from the Individuals in the patch.
	std::vector<NCompetingIndividual> competing_individuals;
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		// Calculate whole-plant N:C status 
		// TODO: revisit this calculation and check where cton_status is used.
		double ntoc = (indiv.get_nmass_leaf() + indiv.nmass_root + indiv.nmass_crown + indiv.nmass_sap) / (indiv.get_cmass_leaf() + indiv.cmass_root + indiv.cmass_crown + indiv.cmass_sap);
		double ntoc_max = 1.0 / indiv.pft.cton_leaf_min;
		double ntoc_avr = 1.0 / indiv.pft.cton_leaf_avr;
		indiv.cton_status = max(0.0, (ntoc - ntoc_max) / (ntoc_avr - ntoc_max));  // C:N relative to average. Input to nitrogen_uptake_strength()
		NCompetingIndividual competing_individual;
		competing_individual.ndemand = indiv.ndemand;
		competing_individual.fnuptake = 1;
		competing_individual.strength = nitrogen_uptake_strength(indiv);
		competing_individuals.push_back(competing_individual);
	}

	// Calculate N competition.
	double nmass_avail = soil.nmass_avail() * min(1.0, patch.rpc_total);

	ncompete(competing_individuals, nmass_avail);  // TODO: check if OK/efficient to call ncompete even if ndemand < nsupply. 
												   // --> if (nmass_avail < patch_total_ndemand) {

	// Copy results back into Individual objects in the patch, and calculate
	// total patch N demand.
	patch.ndemand = 0;

	// Calculate N limitation uptake due to limited N availability in soil
	// soil available mineral nitrogen (kgN/m2)
	const double nmin_avail = soil.nmass_avail();

	/// Rate of nitrogen uptake not associated with Michaelis-Menten Kinetics (Zaehle and Friend 2010)
	double kNmin = 0.05;

	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		// Nitrogen availablilty scalar due to saturating Michealis-Menten kinetics
		double nmin_scale = kNmin + nmin_avail / (nmin_avail + gridcell.pft[indiv.pft.id].Km);
        
		// Nitrogen availability based on root projective cover
		double max_indiv_avail = min(1.0, indiv.rpc);

		// Nitrogen available for uptake
		double favail_nuptake = min(nmin_scale, max_indiv_avail);

		indiv.fnuptake = competing_individuals[i].fnuptake * favail_nuptake;
		indiv.ndemand = competing_individuals[i].ndemand;

		patch.ndemand += indiv.ndemand;

		// Taken out of soil in somdynam::vegetation_n_uptake():
		// indiv.nmass_storage += indiv.ndemand * indiv.fnuptake;
		// This happens right after growth in simulate_day.

	}

	// If daily tree or grass allocation is disabled, we should probably call
	// nstore_usage. This will have no effect on individuals for whom the dave
	// daily allocation is enabled.
	//
	// This function sets indiv.nstress, and moves N from the nstore_labile into
	// the other pools if today's N demand exceeds the supply. This is not
	// required if using dave allocation.
	//
	// In the trunk version of the model, this function is called from canopy
	// exchange, but when dave growth is enabled we call it from growth, because
	// that's where our N demand and uptake occurs.
	if (!ifdailytree || !ifdailygrass) {
		nstore_usage(vegetation);
	}
}



// Update C and N pools of vegetation  
void update_pools_dave(Vegetation& vegetation) { 
	
	Patch& patch = vegetation.patch;
	
	// new parameters. TODO: move to pft parameters(?)
	// maximum mobilisation rate of N in storage
	//const double max_mobil_nstorage = 0.9;
	
	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (!pft.dave_allocation()) {
			// Dave allocation is not enabled for this individual. This would
			// occur, for example, if it's a grass and ifdailygrass is false.
			continue;
		}

		double creturn = 0.0; // kgC m-2 d-1
		double cexcess = 0.0;
		double cgrow_unlim = indiv.cgrow;

		// Always fill up N storage during freenyears
		if ((!ifnlim || date.year < freenyears) && indiv.nmass_storage < indiv.nstorage_max ) {
			// We need to reduce nmass_missing by the uptake amount, otherwise
			// the storage pool will exceed its capacity when uptake occurs.
			// However, if the storage pool is already close to capacity, then
			// subtracting this will result in negative nmass_missing.
			indiv.nmass_missing = max(0.0, indiv.nstorage_max - indiv.nmass_storage - (indiv.fnuptake * indiv.ndemand));
			indiv.nmass_storage += indiv.nmass_missing;
			patch.fluxes.report_flux(Fluxes::NFERT, indiv.nmass_missing);
		} else {
			indiv.nmass_missing = 0;
		}

		// Make sure N storage does not exceed its capacity. This mostly happens
		// when nstorage_max decreases as the plant 'shrinks' (extreme stress)
		// Flush to root litter for now, but better alternatives could exist
		if (indiv.nmass_storage > indiv.nstorage_max) {
			double nexcess = indiv.nmass_storage - indiv.nstorage_max;
			indiv.ntolitter += nexcess;
			ppft.nmass_litter_root += nexcess;
			indiv.nmass_storage = indiv.nstorage_max;
		}

		// TODO: Check if there is enough N in storage. At the moment Nstorage needs to be > 0.
		// Alternatively, one could implement a N:C storage > N:C leaf min constraint.
		
		double fnavail = 0.0;
		double nstress = 0.0;
		if (!negligible(indiv.ndemand)) {
			//fnavail = max(0.0,min(1.0,(indiv.nmass_storage / indiv.ndemand)));
			// TODO: define max_mobil_nstorage based on Yang et al. 2023, Journal of Ecology.
			fnavail = max(0.0,min(1.0,((max_mobil_cstorage * indiv.nmass_storage) / (indiv.ndemand - indiv.storendemand))));
			nstress = max(0.0,min(1.0,(indiv.nmass_storage / (indiv.ndemand - indiv.storendemand))));
		}
		
		// calculate N stress factor 
		if (!ifnlim || date.year < freenyears) {
			indiv.nscal_5.add(1.0);
		} else {
			indiv.nscal_5.add(nstress);
		}

		double navail_leaf  = indiv.leafndemand * fnavail;
		double navail_root  = indiv.rootndemand * fnavail;
		double navail_crown = indiv.crownndemand * fnavail;
		double navail_sap  = indiv.sapndemand * fnavail; 
		double navail_repr  = indiv.reprndemand * fnavail; 


		// Check with fnavail, not fnuptake. There are cases where fnuptake can
		// be 1 but fnavail < 1. This can happen when n storage is low, because
		// the uptake doesn't occur until after growth. If we were to check
		// fnuptake in this case, we would incorrectly assume that we are able
		// to meet the N demands today, and potentially allocate too much C,
		// causing the C:N ratio to exceed its max.
		//if (negligible(fnavail - 1.0)) {
			
			// N supply >= N demand

		//	indiv.cmass_leafs[green] += indiv.cgrow_leaf;
		//	indiv.cmass_root += indiv.cgrow_root;
		//	indiv.cmass_crown += indiv.cgrow_crown;
		//	indiv.cmass_repr += indiv.cgrow_repr;
			
		//}
		//else {

		// N supply < N demand
		double nreturn = 0;
			
		// Leaves
		constrain_growth(
			indiv.cgrow_leaf, navail_leaf,
			indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
			creturn, nreturn);

		// Roots
		constrain_growth(
			indiv.cgrow_root, navail_root,
			indiv.pft.cton_root_min, indiv.pft.cton_root_max,
			creturn, nreturn);

		// Crowns (for grasses and mosses only).
		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			constrain_growth(
				indiv.cgrow_crown, navail_crown,
				indiv.pft.cton_root_min, indiv.pft.cton_root_max,
				creturn, nreturn);
		}

		// sapwood
		if (indiv.pft.lifeform == TREE) {
			constrain_growth(
				indiv.cgrow_sap, navail_sap,
				indiv.pft.cton_sap_min, indiv.pft.cton_sap_max,
				creturn, nreturn);
		}

		// Reproductive tissue
		constrain_growth(
			indiv.cgrow_repr, navail_repr,
			indiv.pft.cton_repr_min, indiv.pft.cton_repr_max,
			creturn, nreturn);

		// Make sure storage does not overflow.
		// Maximum C storage pool size
		if (indiv.cmass_storage + creturn > indiv.cstorage_max) {
			cexcess = indiv.cmass_storage + creturn - indiv.cstorage_max;
			indiv.cexcess += cexcess;
			creturn = indiv.cstorage_max - indiv.cmass_storage;
		}
		indiv.cmass_storage += creturn;

		// nreturn doesn't need to go back into the storage pool, because we
		// haven't actually taken anything out of the storage pool yet, and
		// the amount which is taken out is navail_{leaf,root,repr}, which
		// have already been reduced by nreturn. The only reason we have an
		// nreturn variable is for the mass balance checks.

		if (pft.is_grass() || pft.lifeform == MOSS) {
			indiv.cmass_leafs[green] += indiv.cgrow_leaf;
		} else {
			indiv.cmass_leaf += indiv.cgrow_leaf;
		}
		
		indiv.cmass_root  += indiv.cgrow_root;
		indiv.cmass_crown += indiv.cgrow_crown;
		indiv.cmass_sap += indiv.cgrow_sap;
		indiv.cmass_repr  += indiv.cgrow_repr;

		assert(equal((indiv.ndemand - indiv.storendemand) * fnavail, navail_leaf + navail_root + navail_crown + navail_sap + navail_repr + nreturn));

		// Always take up available N
		if (pft.is_grass() || pft.lifeform == MOSS) {
			indiv.nmass_leafs[green] += navail_leaf;
		} else {
			indiv.nmass_leaf += navail_leaf;
		}
		indiv.nmass_root  += navail_root;
		indiv.nmass_crown += navail_crown;
		indiv.nmass_sap += navail_sap;
		indiv.nmass_repr  += navail_repr;
		indiv.nmass_storage = max(0.0, indiv.nmass_storage - (navail_leaf + navail_root + navail_crown + navail_sap + navail_repr));

		// Save yearly leaf and root growth. 
		// TODO: check if crown and repr need to be added
		indiv.cmass_leaf_ygrowth += indiv.cgrow_leaf;
		indiv.cmass_root_ygrowth += indiv.cgrow_root;

		// update indiv.cgrow (actual growth)
		indiv.cgrow = indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_crown + indiv.cgrow_sap + indiv.cgrow_repr;

		assert(positive(indiv.cgrow_leaf));
		assert(positive(indiv.cgrow_root));
		assert(positive(indiv.cgrow_crown));
		assert(positive(indiv.cgrow_sap));
		assert(positive(indiv.cgrow_repr));
		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));
		assert(less_than_or_equal(indiv.nmass_storage, indiv.nstorage_max));
		assert(equal(cgrow_unlim, indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_crown + indiv.cgrow_sap + indiv.cgrow_repr + creturn + cexcess));

		// Checking C:N ratios

		// Leaf
		assert(less_than_or_equal(indiv.cton_leaf(), indiv.pft.cton_leaf_max * 1.001));
		assert(greater_than_or_equal(indiv.cton_leaf(), indiv.pft.cton_leaf_min * 0.999));

		// Root
		assert(less_than_or_equal(indiv.cton_root(), indiv.pft.cton_root_max * 1.001));
		assert(greater_than_or_equal(indiv.cton_root(), indiv.pft.cton_root_min * 0.999));

		// Sapwood (trees only!).
		if (indiv.pft.lifeform == TREE) {
			assert(indiv.cmass_sap >= 0.0);
			assert(less_than_or_equal(indiv.cton_sap(), indiv.pft.cton_sap_max * 1.001));
			assert(greater_than_or_equal(indiv.cton_sap(), indiv.pft.cton_sap_min * 0.999));
			assert(positive(indiv.cmass_sap));
		}

		if (indiv.pft.is_grass()) {
			assert(less_than_or_equal(indiv.cton_crown(), indiv.pft.cton_root_max * 1.001));
			assert(greater_than_or_equal(indiv.cton_crown(), indiv.pft.cton_root_min * 0.999));
			assert(positive(indiv.cmass_crown));
		}
	}
}

void senescence_grass(Individual &indiv) {

	// Parameters
	const double crelocfrac = 0.0;      // TODO: move to constants?
	const double nrelocfrac = 0.5;
	const double senesc_base = 0.00022;         // Baseline senescence rate (d-1)
	const double senesc_scal_dorm = 300.0;     // senescence scalar for dormancy     --> x * 100% increase in senescence rate under max. temp. stress
	// senescence scalar for water stress --> x * 100% increase in senescence rate under max. water stress
	const double senesc_scal_water = indiv.pft.lifeform == GRASS_PERENNIAL ? 60.0 : 60.0;
	// senescence scalar for reproductive phase
	const double senesc_scal_repr = indiv.pft.lifeform == GRASS_PERENNIAL ? 40.0 : 40.0;
	// higher value --> less (more) senescence at beginning (end) of reproductive phase
	const double qdev_senesc   = indiv.pft.lifeform == GRASS_PERENNIAL ? 5.0 : 20.0;    
    // higher value --> lower sensitivity to water stress at high wscal
	const double qwater_senesc = indiv.pft.lifeform == GRASS_PERENNIAL ? 1.0 : 2.0;     
	// increase in senescence rate per 1degC (Johnson & Thornley 1983)
	const double temp_senesc1 = 1.0;    // for temperatures below 0 degC
	const double temp_senesc2 = 0.04;   // for temperatures above 0 degC  

	// TODO: is there an active retranslocation of storage C and N in senescing annual plants?

	Vegetation &vegetation = indiv.vegetation;
	Patchpft ppft = indiv.patchpft();
	Pft &pft = indiv.pft;

	// Air temperature or bottom snow layer if snow present
	double grass_temp = vegetation.patch.soil.T_soil[6];

	// Baseline (minimum) senescence rate from leaf longevity under unstressed conditions.
	// double senesc_base = 1.0 / (indiv.pft.leaflong * date.year_length());
	// Note: leaflong not used anymore for senescence as it is conceptually difficult

	double senesc = senesc_base;

	if (indiv.growth_stage == reproductive || indiv.growth_stage == dormant || indiv.fdev_growth < 1.0) {

		// stronger increase of senescence in annuals
		//double senesc_rate_repr = indiv.pft.lifeform == GRASS_PERENNIAL ? 1.015 : 1.08;
		//indiv.fdev_senesc = max(1.0, indiv.fdev_senesc * senesc_rate_repr);
		//indiv.fdev_senesc  = min(25.0, pow(indiv.fdev_growth, -0.1));
		//indiv.fdev_senesc = senesc_scal_repr * (1.0 - indiv.fdev_growth);  // senescing too early
		indiv.fdev_senesc = senesc_scal_repr * pow((1.0 - indiv.fdev_growth),qdev_senesc);

	} else {

		indiv.fdev_senesc = 0.0;

	}

	// Stress factors increase senescence
	// Soil moisture
	double senesc_fwater = senesc_scal_water * (1.0 - pow(ppft.wscal,qwater_senesc));
	
	// High and low temperatures
	// linear increase as in Johnson & Thornley 1983.
	double senesc_ftemp = 0.0;
	if (grass_temp < 0.0){
		senesc_ftemp = -temp_senesc1 * grass_temp;
	} else {
		if (grass_temp < indiv.pft.temp_dormancy) {
			senesc_ftemp = temp_senesc2 * (indiv.pft.temp_dormancy - grass_temp);
		} else if (grass_temp >= indiv.pft.temp_dormancy) {
			senesc_ftemp = temp_senesc2 * (grass_temp - 25.0);
		}
	}

	// Dormancy-induced
	double senesc_fdorm = 0.0;
	if (indiv.fphen_activity < 1.0) {

		//senesc_fdorm = senesc_scal_dorm * (1.0 - indiv.fphen_activity);
        senesc_fdorm = senesc_scal_dorm * (1.0 - ((indiv.fphen_activity - fphen_activity_min) / (1.0 - fphen_activity_min)));
	}

	// Total senescence rate
	//senesc *= ((1.0 + max(indiv.fdev_senesc,senesc_fdorm)) * (1.0 + senesc_fwater) * (1.0 + senesc_ftemp));
	senesc *= (1.0 + max(indiv.fdev_senesc,senesc_fdorm) + (1.0 + senesc_fwater)) * (1.0 + senesc_ftemp);
	senesc = min(1.0, senesc);
	indiv.senesc_root = senesc_base;

	// Retranslocation and pool updates
	// TODO: maybe call cstorage_capacity once at the beginning to be consistent across routines!
	// TODO: calculate retranslocation coefficients based on current C:N ratios?

	// Carbon
	double csenesc  = indiv.cmass_leafs[green] * senesc;
	double cstore   = csenesc * crelocfrac;
	double ctobrown = csenesc * (1.0 - crelocfrac);

	// If C storage pool overflows, re-adjust cstore and ctobrown
	if (indiv.cmass_storage + cstore > indiv.cstorage_max) {
		cstore = indiv.cstorage_max - indiv.cmass_storage;
		ctobrown = csenesc - cstore;
	}
	
	indiv.cmass_leafs[green] -= csenesc;   // C from green leaf pool
	indiv.cmass_leafs[brown] += ctobrown;  // C to brown leaf pool
	indiv.cmass_storage      += cstore;    // C to storage
	

	// Nitrogen
	double nsenesc  = indiv.nmass_leafs[green] * senesc;
	double nstore   = nsenesc * nrelocfrac;
	double ntobrown = nsenesc * (1.0 - nrelocfrac);

	// If N storage pool overflows, re-adjust nstore and ntobrown
	if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
		nstore = indiv.nstorage_max - indiv.nmass_storage;
		ntobrown = nsenesc - nstore;
	}

	indiv.nmass_leafs[green] -= nsenesc;   // N from green leaf pool
	indiv.nmass_leafs[brown] += ntobrown;  // N to brown leaf pool
	indiv.nmass_storage      += nstore;    // N to storage
		
	// TODO: make sure N:C ratio in storage does not fall outside its ranges
	//       understand what the N:C in storage looks like (values). How flexible is N:C in storage?
	//		 sufficient to only account for Nstorage_max as done here?

	assert(equal(csenesc, cstore + ctobrown));
	assert(equal(nsenesc, ntobrown + nstore));
}

// Transfer of green leaf pool to brown leaf pool.
void senescence_tree(Individual &indiv) {
	// Inputs: indiv.cmass_leafs, indiv.nmass_leafs
	// Outputs: indiv.cmass_leafs, indiv.nmass_leafs
	// Do we even need this? If using a single leaf pool, perhaps we just use
	// turnover instead.

}

// Transfer of C and N from green (live) leaf pool to brown (dead) leaf pool
// Accounts for retranslocation of C and N to storage prior to senescence. 
// Note: Leaves and culms are assumed to senesce at the same rate.
void senescence_dave(Vegetation& vegetation) {
	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (!pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}
		if (pft.is_grass() || pft.lifeform == MOSS) {
			senescence_grass(indiv);
		} else {
			senescence_tree(indiv);
		}

		assert(positive(indiv.cmass_leafs[green]));
		assert(positive(indiv.nmass_leafs[green]));
		assert(positive(indiv.cmass_leafs[brown]));
		assert(positive(indiv.nmass_leafs[brown]));

		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));
	}
}


// External seed input 
void seed_dispersal_dave(Vegetation& vegetation) {

	Patch& patch = vegetation.patch;
	for (int i = 0; i < patch.pft.nobj; i++) {
		Patchpft &ppft = patch.pft[i];

		if (!ppft.pft.dave_allocation()) {
			// This pft is not using dave allocation, and is therefore not using
			// the daily establishment. This means the seedbank and seed
			// dispersal will have no effect on this pft.
			continue;
		}

		double ext_seed_c = ext_seed;

		// test for initial 
		if (date.year == 0 && ppft.pft.lifeform == TREE) {
			ext_seed_c *= 10;
		}

		double ext_seed_n = ext_seed_c / ((ppft.pft.cton_leaf_avr + ppft.pft.cton_root_avr) / 2.0);

		ppft.cmass_litter_repr += ext_seed_c;
		ppft.nmass_litter_repr += ext_seed_n;

        patch.fluxes.report_flux(Fluxes::ESTC, -ext_seed_c);
		patch.fluxes.report_flux(Fluxes::ESTN, -ext_seed_n);
	}
}


// Seed turnover (decrease in viability over time).
void turnover_seed(Vegetation& vegetation) {
	
	Patch& patch = vegetation.patch;
	for (int i = 0; i < patch.pft.nobj; i++) {
		Patchpft &ppft = patch.pft[i];

		if (!ppft.pft.dave_allocation()) {
			// This PFT is not using the dave allocation, and therefore the
			// seedbank is handled separately.
			continue;
		}

		// Seed turnover
		
		// Air temperature or bottom snow layer if snow present and soil moisture of first soil layer.
		double surf_temp = vegetation.patch.soil.T_soil[6];
		double surf_water = patch.soil.get_layer_soil_water(0);
		double HA_seed = 50.0;  // activation energy. TODO: obtain from data
		
		// Base turnover rate modified by temperature and water availability
		double base_turnover_seed = ppft.pft.turnover_repr / date.year_length();
		double frac_turnover_seed = Arrhenius_temp_response(base_turnover_seed,surf_temp,HA_seed,0.0,0.0) * surf_water;

		// Carbon
		double cturnover_repr = frac_turnover_seed * ppft.cmass_litter_repr;
		ppft.cmass_litter_repr -= cturnover_repr;
		ppft.cmass_litter_root += cturnover_repr;

		// Nitrogen
		double nturnover_repr = frac_turnover_seed * ppft.nmass_litter_repr;
		ppft.nmass_litter_repr -= nturnover_repr;
		ppft.nmass_litter_root += nturnover_repr;
	}
}


//double turnover_dave(double turnover_rate, double& source, double& sink, double max_sink_size){

//}

// Note: cropallocation.cpp contains a function called turnover_grass(). I am
// therefore avoiding reusing this name here.
void turnover_dave_grass(Individual &indiv) {
	Vegetation &vegetation = indiv.vegetation;
	Patch& patch = vegetation.patch;
	Pft& pft = indiv.pft;
	Patchpft &ppft = patch.pft[indiv.pft.id];

	// Determine turnover rate (same for all individuals in a PFT).
	// Assumed daily rates here (fraction d-1)!
	// Air temperature or bottom snow layer if snow present and soil moisture of first soil layer.
	double surf_temp = vegetation.patch.soil.T_soil[6];
	double surf_water = patch.soil.get_layer_soil_water(0);

	// crown turnover relative to root turnover.
	double scal_crown_turnover = 0.35;

	// Turnover of living tissues (root and crown)
	// senesc_root equals the baseline senescence rate of leaves.
	double frac_turnover_root = 20.0 * indiv.senesc_root;	
	double frac_turnover_crown = scal_crown_turnover * frac_turnover_root;

	// Turnover of brown shoots
	double ftemp = 0.0, fwater = 0.0;
	// Temperature response
	double kleaf = 0.01;  // base respiration rate of brown shoots at 10degC
	respiration_temperature_response(surf_temp, ftemp);
	// Soil moisture response (Parton et al. 2001). 
	// Note that Parton et al. 2001 used PPT/PET instead of wcont.
	//fwater = pow(surf_water,1.0);
	fwater = 1.0 / (1.0 + 30.0 * exp(-8.5 * surf_water));
	double frac_turnover_leaf = kleaf * ftemp * fwater;


	// Dead shoot turnover (brown shoot --> leaf litter)
	if (indiv.cmass_leafs[brown] > 0.0) {

		// Carbon
		double cturnover_leaf = indiv.cmass_leafs[brown] * frac_turnover_leaf;

		indiv.cmass_leafs[brown]  -= cturnover_leaf;
		ppft.cmass_litter_leaf    += cturnover_leaf;
		indiv.ctolitter           += cturnover_leaf;
		indiv.ctolitter_leaf      =  cturnover_leaf; // Diagnostic output



		// Nitrogen
		double nturnover_leaf = indiv.nmass_leafs[brown] * frac_turnover_leaf;
		
		indiv.nmass_leafs[brown]   -= nturnover_leaf;
		ppft.nmass_litter_leaf     += nturnover_leaf;
		indiv.ntolitter            += nturnover_leaf;
		indiv.ntolitter_leaf       =  nturnover_leaf; // Diagnostic output
	}

	// Crown turnover
	// TODO: put into function (same for roots and crowns)!	
	if (indiv.cmass_crown > 0.0) {

		// baseline turnover rate
		double cturnover_crown = indiv.cmass_crown * frac_turnover_crown;

		indiv.cmass_crown        -= cturnover_crown;
		ppft.cmass_litter_root   += cturnover_crown;
		indiv.ctolitter          += cturnover_crown;
		indiv.ctolitter_crown    =  cturnover_crown; // Diagnostic output

		// Nitrogen (including retranslocation)
		double nturnover_crown = indiv.nmass_crown * frac_turnover_crown;
		double nstore    = nturnover_crown * nrelocfrac;
		double ntolitter = nturnover_crown * (1.0 - nrelocfrac);

		if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
			nstore    = indiv.nstorage_max - indiv.nmass_storage;	
			ntolitter = nturnover_crown - nstore;
		}

		indiv.nmass_crown        -= nturnover_crown;
		ppft.nmass_litter_root   += ntolitter;
		indiv.nmass_storage      += nstore;
		indiv.ntolitter          += ntolitter;
		indiv.ntolitter_crown    = ntolitter; // Diagnostic output

		assert(equal(nturnover_crown, nstore + ntolitter));
	}


	// Root turnover (roots --> root litter)
	// TODO: should this be in senescence as well? Dead root pool?
	if (indiv.cmass_root > 0.0) { 

		// TODO: include C resorption from roots?
		// Carbon
		double cturnover_root = indiv.cmass_root * frac_turnover_root;
		
		indiv.cmass_root        -= cturnover_root;
		ppft.cmass_litter_root  += cturnover_root;
		indiv.ctolitter         += cturnover_root;
		indiv.ctolitter_root    =  cturnover_root; // Diagnostic output

		// Nitrogen (including retranslocation)
		double nturnover_root = indiv.nmass_root * frac_turnover_root;
		double nstore    = nturnover_root * nrelocfrac;
		double ntolitter = nturnover_root * (1.0 - nrelocfrac);

		if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
			nstore    = indiv.nstorage_max - indiv.nmass_storage;	
			ntolitter = nturnover_root - nstore;
		}

		indiv.nmass_root        -= nturnover_root;
		ppft.nmass_litter_root  += ntolitter;
		indiv.ntolitter         += ntolitter;
		indiv.ntolitter_root    =  ntolitter; // Diagnostic output

		indiv.nmass_storage     += nstore;

		assert(equal(nturnover_root, nstore + ntolitter));
	}
}

// Turnover of green pools to litter. 
// TODO: If only using a single leaf pool, we could remove the senescence function altogether.
void turnover_tree(Individual &indiv) {
	// Inputs:  indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap
	//          indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap,
	//          indiv.pft.leaflong
	//
	// Outputs: indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap
	//          indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap
	//          ppft.cmass_litter_leaf, ppft.cmass_litter_root, indiv.cmass_heart

	// fraction of senesced N that goes to storage (-). 
	// TODO: also in turnover_dave() function. Further, check for literature 
	// evidence.
	const double nrelocfrac = 0.5;
	// maximum relative change in leaf and root senescence as a result of  
	// water stress (made up)
	const double senesc_scal_water = 4; 
	// slope of curve of the change in leaf and root senescence as a result of 
	// water stress (made up)
	const double qwater_senesc = 0.7; 
	// how much turnover should be reduced when tree age is younger than leaf 
	// longevity OR young trees under water strees
	const double ftreeage = 2;
	// shape of function of fage (age-based modifier to scenescence) as a function of age
	// shp_fage = 1 gives a linear increase from 0 to age = leaflong
	const double shp_fage = 1.0;
	// turnover of leaves initialisation
	double turnover_leaf = 0;
	// turnover of root initialisation
	double turnover_root = 0;
	// turnover of leaves initialisation
	double turnover_sap = 0;

	Patchpft& ppft = indiv.patchpft();
	
	// fraction of change in leaf and root turnover as a result of water stress.
	// leaf (root) turnover will increase (decrease) by the same fraction 
	// in high soil water stress conditions.
	// TODO: function seems too extereme (too much turnover for low water stress).
	double fwater = max( 1.0, senesc_scal_water * ( 1 - pow(ppft.wscal_30.mean(), qwater_senesc) ) );
	// double fwater = 1 + (senesc_scal_water * pow(( 1- ppft.wscal_30.mean()), qwater_senesc)); // OLD

	//
	double cturnover_leaf = 0.0;
	double cturnover_root = 0.0;
	double cturnover_sap = 0.0;
	// scenescence of leaves (also depending on water stress (fwater)). 
	// TODO: 1. no regards to what this would do to actual k_latosa tomorrow; 
	// 2. Should the base turnover (pft.turnover_leaf) be more related to leaf longevity?; 
	// 3. leaf turnover shoudl be turned off in very young trees (until the leaves 
	// reach leaf longevity). Alternatively, reduce base leaf turnover by a fraction
	// in the first few years of the tree. For example, in the first years, turnover_leaf
	// will be 1/2 of it's value when the tree is older than leaf longevity:
	// (e.g: turnover_leaf = indiv.pft.turnover_leaf / 2 * fwater / Date::MAX_YEAR_LENGTH)
	if (indiv.cmass_leaf > 0) {		
		// TODO: fix this. Currently, turnover will bee too small when there is water 
		// stress (fwater > 1), and age is below leaf longevity. Possible to add 
		// another variable to increases the weight of the effect of water  
		// stress at young cohort age. for exmplae: 
		// if (indiv.age < indiv.pft.leaflong), fwater will be pow(fwater, 3)
		// modifier of turnover a a function of cohort age.
		double fage = 1;
		// fage (unlimited) (ie can be higher than 1)
		fage = pow(((indiv.age) / indiv.pft.leaflong), shp_fage);
		// limit fage
		fage = min(fage, 1.0);
		fage = 1; // temporary disabled age effect
		// fraction of leaf C that is turned over today
		turnover_leaf = indiv.pft.turnover_leaf * fwater * fage / Date::MAX_YEAR_LENGTH;

		// Leaf C amount to turnover (kgC/m^2)
		cturnover_leaf = turnover_leaf * indiv.cmass_leaf; 
		// Leaf N amount to turnover (kgN/m^2)
		double nturnover_leaf = turnover_leaf * indiv.nmass_leaf;
		// Leaf N that goes to litter (kgN/m^2)
		double n2litter_leaf = nturnover_leaf * (1.0 - nrelocfrac);
		// Leaf N amount to storage (kgC/m^2)
		double n2storage_leaf = nturnover_leaf * nrelocfrac;
		// prevent overflow of N storage:
		if ((indiv.nmass_storage + n2storage_leaf) > indiv.nstorage_max) {
			// only fill up storage to capacity
			n2storage_leaf = indiv.nstorage_max - indiv.nmass_storage;
			// update N that goes to litter (when there is N overflow above capacity)
			n2litter_leaf = nturnover_leaf - n2storage_leaf;
		}

		// updete pools
		ppft.cmass_litter_leaf += cturnover_leaf;
		indiv.ctolitter += cturnover_leaf;
		indiv.cmass_leaf -= cturnover_leaf;

		ppft.nmass_litter_leaf += n2litter_leaf;
		indiv.ntolitter += n2litter_leaf;
		indiv.nmass_leaf -= nturnover_leaf;
		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_leaf;

		assert(equal(nturnover_leaf, n2storage_leaf + n2litter_leaf));
	}
	// save daily leaf turnover flux (kgC/m^2/day)
	indiv.cturnover_leaf = turnover_leaf;

	if (indiv.cmass_root > 0) {
		// fraction of root biomass that is turned over today. 
		turnover_root = indiv.pft.turnover_root / fwater / Date::MAX_YEAR_LENGTH;

		// root C amount to turnover (kgC/m^2)
		cturnover_root = turnover_root * indiv.cmass_root; 
		// root N amount to turnover (kgN/m^2)
		double nturnover_root = turnover_root * indiv.nmass_root;
		// Lean N that goes to litter (kgN/m^2)
		double n2litter_root = nturnover_root * (1.0 - nrelocfrac);
		// root N amount to storage (kgC/m^2)
		double n2storage_root = nturnover_root * nrelocfrac;
		// prevent overflow of N storage:
		if ((indiv.nmass_storage + n2storage_root) > indiv.nstorage_max) {
			// only fill up storage to capacity
			n2storage_root = indiv.nstorage_max - indiv.nmass_storage;
			// update N that goes to litter (when there is N overflow above capacity)
			n2litter_root = nturnover_root - n2storage_root;
		}

		// updete pools
		ppft.cmass_litter_root += cturnover_root;
		indiv.ctolitter += cturnover_root;
		indiv.cmass_root -= cturnover_root;

		ppft.nmass_litter_root += n2litter_root;
		indiv.ntolitter += n2litter_root;
		indiv.nmass_root -= nturnover_root;
		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_root;

		assert(equal(nturnover_root, n2litter_root + n2storage_root));
	}
	// save daily root turnover flux (kgC/m^2/day)
	indiv.cturnover_root = cturnover_root;

	if (indiv.cmass_sap > 0) {
		// if there isn't any leaf turnover today, then there shouldn't be sapwood 
		// turnover becasue current sapwood (cmass_sap) is needed to maintain the 
		// current leaves.
		// TODO: check if plausible (?). indiv.pft.k_latosa > 0 is only checked 
		// to prevent dividing by 0 (which would happen if k_latosa is a variable)
		
		if ((turnover_leaf > 0) && (indiv.pft.k_latosa > 0)) {
			// turnover of sapwood will be calculated from the turnover of leaf
			// in order to maintain k_latosa (derivation in notebook 22/11/2023)
			
			double lambda = (indiv.height * indiv.sla() * indiv.pft.wooddens)/indiv.pft.k_latosa;
			lambda = min(lambda, 1.0);
			lambda = max(lambda, 0.0);

			// sapwood turnover 
			turnover_sap =  lambda * turnover_leaf;
			// double lambda = (indiv.height * indiv.sla() * indiv.pft.wooddens) / indiv.pft.k_latosa;
			// // sapwood turnover 
			// turnover_sap =  lambda * turnover_leaf; // oldish
			// turnover_sap =  indiv.pft.turnover_sap / 365; // trunk
			// turnover_sap =  min( lambda * turnover_leaf, lambda * turnover_leaf );

		} else {
			// no turnover sap today
			turnover_sap = 0;
		}

		// sapwood C amount to turnover (kgC/m^2)
		cturnover_sap = turnover_sap * indiv.cmass_sap; 
		// sapwood N amount to turnover (kgN/m^2)
		double nturnover_sap = turnover_sap * indiv.nmass_sap;
		// Lean N that goes to heartwood (kgN/m^2)
		double n2heart_sap = nturnover_sap * (1.0 - nrelocfrac);
		// sapwood N amount to storage (kgC/m^2)
		double n2storage_sap = nturnover_sap * nrelocfrac;
		// prevent overflow of N storage:
		if ((indiv.nmass_storage + n2storage_sap) > indiv.nstorage_max) {
			// only fill up storage to capacity
			n2storage_sap = indiv.nstorage_max - indiv.nmass_storage;
			// update N that goes to heartwood (when there is N overflow above capacity)
			n2heart_sap = nturnover_sap - n2storage_sap;
		}

		// updete pools
		indiv.cmass_heart += cturnover_sap;
		indiv.cmass_sap -= cturnover_sap;

		indiv.nmass_heart += n2heart_sap;
		indiv.nmass_sap -= nturnover_sap;
		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_sap;
		assert(equal(nturnover_sap, n2storage_sap + n2heart_sap));
	}
	// save daily sapwood turnover flux (kgC/m^2/day)
	indiv.cturnover_sap = cturnover_sap;

	// // drop reproductive pool when trends goes down
	// // TODO: we need a better solution here to account for flushing of the 
	// // reproductiove pool into the litter. possible do a trend rather than all 
	// // at once.
	// if (indiv.cgrow_repr_20.periodictrend(20) < 0) {
	// 	indiv.flush_repr_litter_DAVE();
	// }

	assert(indiv.cmass_leaf >= 0);
	assert(indiv.cmass_root >= 0);
	assert(indiv.cmass_sap >= 0);
	assert(indiv.cmass_heart >= 0);
	assert(indiv.cmass_storage >= 0);

}

// Transfer of crown pool, brown shoot pool, and fine root pool to leaf and root litter.
// Transfer reproductive tissue to litter if end of growth cycle is reached.
void turnover_dave(Vegetation& vegetation) {

    // TODO: move to constants
	const double nrelocfrac = 0.5;

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}
		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			turnover_dave_grass(indiv);
		} else {
			turnover_tree(indiv);
		}
	}

	turnover_seed(vegetation);
}


// Allocate excess C and N to target pool.
void allocate_excess(double &cexcess, double &nexcess, double &cpool_dest, double &npool_dest,
					 double &cmass_storage, double &nmass_storage, double cton_min, double cton_max,
					 double calloc_max) {

	if ((nmass_storage + nexcess) > 0.0) {

		double cdest = 0.0;
		double ndest = 0.0;

		// if there is excess N, try to allocate all of it to the destination pool
		if (nexcess > 0) {
			cdest = nexcess * cton_min;
			ndest = nexcess;
		}

		// allocate as much of cexcess as we can (given N excess and C:N ratios)
		if (cexcess > 0) {
			cdest = max(cdest, cexcess);
			ndest = max(ndest, cdest / cton_max);
		}

		// case where cdest exceeds available C storage and C excess (could happen under large N excess)
		// Note: half of max_mobil_cstorage is taken as cdest is allocated to two pools.
		if (cdest > (cmass_storage * 0.5 * max_mobil_cstorage + cexcess)) {
			cdest = (cmass_storage * 0.5 * max_mobil_cstorage + cexcess);
			ndest = cdest / cton_min;
		}

		// case where ndest exceeds available N storage and N excess (could happen under large C excess)
		if (ndest > (nmass_storage * 0.5 * max_mobil_cstorage + nexcess)) {
			ndest = (nmass_storage * 0.5 * max_mobil_cstorage + nexcess);
			cdest = ndest * cton_min;
		}

		if (cdest > calloc_max) {
			cdest = calloc_max;
			ndest = calloc_max / cton_max;
		}

		// add allocated cexcess to indiv.cgrow for use in respiration routine
		// TODO: is this needed?
		// indiv.cgrow += cdest; 

		// TODO Drew: verify floating point accuracy. Perhaps create an else-if branch.
		if (cexcess <= cdest) {
			// cannot be fulfillled by cexcess (could happen under high N excess)
			cmass_storage -= (cdest - cexcess);
			cexcess = 0.0;
		} else {
			// can be fulfilled by cexcess
			cexcess -= cdest;
		}

		if (nexcess <= ndest) {
			// cannot be fulfillled by nexcess (could happen under high C excess)
			nmass_storage -= (ndest - nexcess);
			nexcess = 0.0;
		} else {
			// can be fulfilled by nexcess
			nexcess -= ndest;
		}

		cpool_dest += cdest;
		npool_dest += ndest;

	}

}


// Determine what to do with excess C when storage exceeds its capacity.
// This routine is also used to actively redistribute C storage to 
// reproduction in some cases. Strictly, this should be done somewhere else,
// but it's more practical to deal with this in one place.
void storage_overflow(Individual& indiv) {

	Pft& pft = indiv.pft;
	Patchpft& ppft = indiv.patchpft();

	// for annual grasses in late reproductive phase, redistribute C storage to reproduction.
	// If LAI is very low, move all storage to reproduction to kill the annual.
	// If this is not done, there is a risk that the plant will continue living with a very
	// low storage for a very long time. 
	// TODO: have a look at respiration of reproductive pool. This could solve the problem as well.
	if (indiv.pft.lifeform == GRASS_ANNUAL && indiv.growth_stage == reproductive) {
		// TODO: limit allocation to repr for annuals to avoid unrealistic growth rates? 
		// TODO: refine frac to be more consistent with the rest of the code!
		double frac = max_mobil_cstorage * indiv.cmass_repr / (indiv.cmass_repr + indiv.get_cmass_leaf()); 
		if (indiv.lai < 0.001) {
			frac = 1.0;
		}

		indiv.cexcess += frac * indiv.cmass_storage; 
		indiv.nexcess += frac * indiv.nmass_storage;

		indiv.cmass_storage -= frac * indiv.cmass_storage;
		indiv.nmass_storage -= frac * indiv.nmass_storage;
	}

	if (indiv.cexcess > 0.0 || indiv.nexcess > 0.0) {

		// decide where to put excess C (or excess N).
		// If the plant is in the reproductive phase, put excess into repr
		if (indiv.growth_stage == reproductive) {
			
			double cexcess_repr = min(indiv.cexcess, indiv.cgrow_repr);

			// determine amount of C excess allocated to reproductive pool
			// Note: allocate as much as possible of the C and N excess, but it needs to be 
			//       constrained by stoichiometry. Crepr is the amount of C that was actually
			//		 allocated to the reproductive pool. Function will draw down storage if 
			//		 there is a strong imbalance between cexcess and nexcess.
			allocate_excess(indiv.cexcess, indiv.nexcess, indiv.cmass_repr, indiv.nmass_repr, indiv.cmass_storage, 
							indiv.nmass_storage, pft.cton_repr_min, pft.cton_repr_max, cexcess_repr);	
		}

		// Put as much of possible of the rest of C excess into the crown (perennials only).
		double clim_crown = indiv.pft.lifeform == GRASS_PERENNIAL ? 0.5 : 0.0;
		double cexcess_crown = min(indiv.cexcess, clim_crown * indiv.cgrow);
		allocate_excess(indiv.cexcess, indiv.nexcess, indiv.cmass_crown, indiv.nmass_crown, indiv.cmass_storage, 
						indiv.nmass_storage, pft.cton_root_min, pft.cton_root_max, cexcess_crown);

	}

	// Now deal with excess C and N that cannot be allocated because of stoichiometry or
	// because they would lead to excessive growth.
	if (indiv.cexcess > 0.0) {
		// TODO: decide what to do with excess C if you are not in the reproductive phase
		// Option 1: respire
		indiv.report_flux(Fluxes::RA, indiv.cexcess);
		indiv.report_flux(Fluxes::NPP, -indiv.cexcess);

		// Option 2; root exudates
		// TODO: yet to be implemented. Alternatively just allocate to roots.

		// Option 3?

		// Reset cexcess to 0.0
		// Now done at the beginning of the day in storage_growth_partitioning_dave() so that mass balance checks work.
		// indiv.cexcess = 0.0;
	}

	if (indiv.nexcess > 0.0) {
		// flush to litter
		ppft.nmass_litter_root += indiv.nexcess;
		indiv.ntolitter += indiv.nexcess;
	}
}


// Photosynthetic sink limitation.
// Calculates downregulation factors applied to photosynthetic parameters.
void sink_source_feedback(Individual& indiv) {
	
	double dmax_pslim = 0.015;
	// Do not let ps_sink_lim vary by more than dmax_pslim d-1
	double const dpslim = dmax_pslim * (1.0 - indiv.fdev_growth);

	// TODO: revisit cexcess / dgpp
	if (indiv.cexcess > 0.0 && indiv.dgpp > 0.0 && is_zero(indiv.fphen_activity - 1.0)) {
		indiv.ps_sink_lim = max((1.0 - dpslim) * indiv.ps_sink_lim, (1.0 - min(1.0, indiv.cexcess / indiv.dgpp)));
	} 
	else {
		// assuming faster recovery than decline
		indiv.ps_sink_lim = min(1.0, indiv.ps_sink_lim * (1.0 + dpslim * 2.0)); 
	}
} 


// Housekeeping of C storage. Handles overflow and feedback with source dynamics.
void storage_dynamics_dave(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		
		Individual& indiv = vegetation[i];
		Patchpft& ppft = indiv.patchpft();

		if (!indiv.pft.dave_allocation()) {
			// This PFT is not using the dave daily allocation scheme. Therefore
			// we should skip it; the storage dynamics are handled in growth().
			continue;
		}

		// TODO: make sla age-dependent or make culmfrac dynamic

		// Update capacity of C and N storage
		indiv.cstorage_max = cstorage_capacity(indiv);
		indiv.nstorage_max = nstorage_capacity(indiv);

		if (indiv.cstorage_max < indiv.cmass_storage) {
			indiv.cexcess += indiv.cmass_storage - indiv.cstorage_max;  // TODO: check if cexcess here comes only from decreasing plant size.
			indiv.cmass_storage = indiv.cstorage_max;
		}
		if (indiv.nstorage_max < indiv.nmass_storage) {
			indiv.nexcess += indiv.nmass_storage - indiv.nstorage_max;
			indiv.nmass_storage = indiv.nstorage_max;
		}

		// Determine in how far C storage induces a negative feedback on photosynthetic capacity.
		sink_source_feedback(indiv);

		// Handle excess C and N.
		storage_overflow(indiv);

		// keep track of changes in C storage and reproduction.
		if ((indiv.age < 0.5 && indiv.frac_germ > 0.01) || (indiv.age < 1.0 / 365.0)) {
			// do not 'count' days with germination in annuals and newly established perennials.
			// Also do not include the very first day.
			// Otherwise repr switch will happen too early.
			// Note: right now it checks germination of the previous day, which is OK in this context.
			indiv.cmass_storage_20.add(1.0);
		} else {
			indiv.cmass_storage_20.add(indiv.cmass_storage);
		}
		indiv.cmass_repr_10.add(indiv.cmass_repr);
		indiv.frac_germ = 0.0;  // reset in vegdynam.cpp if > 0.0

		indiv.cgrow_repr_20.add(indiv.cgrow_repr);

	}
}


void calc_tree_height_radius_increment(Individual& indiv, double& dh, double& dr) {
	/* we now incrementally grow trees based on the cmass stem (sap) increase
	 Still assuming sap+heart wood is a cylinder with height h and radious R
	 if I have a dWcyl (weight increase of cylinder), how much is dr and dh? We know the current R, and h and Vcyl
	 two unknowns: dr and dh?
	 need two equations:
	 (1) Wcyl= pi * wooddens * R^2 * h
	 
	 partial derivative Wcyl:
	 dWcyl = @Wcyl/@r  * dr + @Wcyl/@h * dh
	 @Wcyl/@r = pi*wooddens * 2 * R * h
	 @Wcyl/@h = pi*wooddens * R^2
	 
	 (2) dWcyl = Fr(R,h) * dr + Fh(R)  * dh
	 (2a) Fr(R,h)   = pi*wooddens * 2 * R * h
	 (2b) Fh(R)     = pi*wooddens * R^2
	 we can calculate the Fr(R,h) and Fh(R) for a current tree.
	 
	 Need a second equation that relates R and h, also have:
	 (3) h = k_allom2 * (2*R)^k_allom3    see allometry() Eq. 4 with stem diameter = 2*R
	 = k_allom2 * 2^k_allom3 * R^k_allom3
	 derivative of (3) with respect to R:	 Note: most PFTs k_allom3 = 0.67, k_allom2 = 60
	 (4) dh(R) = k_allom2 * 2^k_allom3 * k_allom3 * R^(k_allom3-1) * dr = Fallom(R) * dr
	 (4a) Fallom(R) = k_allom2 * 2^k_allom3 * k_allom3 * R^(k_allom3-1)
	 Combine (2) and (4):
	 dWcyl = (Fr(R,h) + Fh(R) * Fallom(R)) * dr
	 Solve for dr and then dh:
	 (5) dr = dWcyl/(Fr(R,h) + Fh(R) * Fallom(R))
	 (6) dh = Fallom(R) * dr
	 So have dr and dh based on dWcyl.
	 */
	
	// need to estimate initial height and diam in allocation_init()
	// might work already since allometry() is called after allocation_init() with indiv.age=0
	double h = indiv.height;
	double R = indiv.diam / 2.0;
	// this can not work for a non existing tree
	assert(!negligible(h) && !negligible(R));
	
	double wooddens = indiv.pft.wooddens;
	double k_allom2 = indiv.pft.k_allom2;
	double k_allom3 = indiv.pft.k_allom3;
	
	double Fr = PI * wooddens * 2.0 * R * h;
	double Fh = PI * wooddens * R * R;
	double Fallom = k_allom2 * pow(2.0, k_allom3) * k_allom3 * pow(R, k_allom3 - 1.0);
	// should calc once: k_allom2 * pow(2.0,k_allom3) * k_allom3
	// we need to calc change for individual trees in cohort
	if (indiv.cgrow_sap < 0) {
		printf("asd");
	}
	double dWcyl = indiv.cgrow_sap / indiv.densindiv;
	dr = dWcyl / (Fr + Fh * Fallom);
	dh = Fallom * dr;

	assert(dr >= 0);
	assert(dh >= 0);

}

void allometry_dave_tree(Individual& indiv) {
		
	Vegetation& vegetation = indiv.vegetation;

	// maximum tree height 
	// TODO: make pft parameter.
	double const MAX_HEIGHT = 100.0; 
	
	// calculate change in height and radius
	if (!negligible(indiv.cgrow_sap)) {
		
		double dh;
		double dr;
		calc_tree_height_radius_increment(indiv, dh, dr);

		indiv.d_height = dh;
		indiv.d_diam = 2.0 * dr;

		indiv.height += indiv.d_height;
		indiv.diam += indiv.d_diam;
		
		// Stem volume
		// double vol = indiv.height * 3.1415 * indiv.diam * indiv.diam * 0.25;
		// Crown area (Eqn 6)

		if (!is_zero(indiv.diam)) {
			indiv.crownarea = min(indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp)
			    , indiv.pft.crownarea_max);
		} 
	}

	if (!is_zero(indiv.cmass_leaf)) {

		if (!is_zero(indiv.crownarea)) {
			
			// Individual LAI (Eqn 9)
			indiv.lai_indiv = indiv.cmass_leaf / indiv.densindiv * indiv.pft.sla / indiv.crownarea;
			
			// FPC (Eqn 8)
			double fpc_new = indiv.crownarea * indiv.densindiv *	(1.0 - lambertbeer(indiv.lai_indiv));
			
			// Increment deltafpc
			indiv.deltafpc += fpc_new - indiv.fpc;
			indiv.fpc = fpc_new;
			
		}
		else {
			indiv.lai_indiv = 0.0;
			indiv.fpc = 0.0;
		}
		// Patch-level LAI
		indiv.lai = indiv.cmass_leaf * indiv.pft.sla;

		// calculate leaf to sapwood area (unitless vatrible)
		if (!is_zero(indiv.height)) {
			double const LATOSA_MAX = 8000.0;
			// sapwood area (m^2)
			indiv.sa = (indiv.cmass_sap) / (indiv.height * indiv.pft.wooddens);
			if (!is_zero(indiv.sa)) {
				// leaf to sapwood area (unitelss)
				indiv.latosa = min(LATOSA_MAX, indiv.lai / indiv.sa);
			} else {
				indiv.latosa = LATOSA_MAX;
			}
			
		} else {
			indiv.sa = 0.0;
			indiv.latosa = 0.0;
		}
	} else {
		indiv.lai_indiv = 0.0;
		indiv.lai = 0.0;
		indiv.fpc = 0.0;
		indiv.sa = 0.0;
		indiv.latosa = 0.0;
	}
	// TODO: to change when we decide that the allocation prority to rots that affects 
	// rpc to be pft-specific.
	indiv.rpc = indiv.fpc;

	assert(indiv.cmass_sap >= 0.0);

}

// Updates LAI and FPC for grasses and treesw.
// Equivalent to default allometry() function at the moment,
// but might develop. (e.g. more than 1 indiv for grasses, separate green and brown lai, 
// specific extinction coefficient etc.)
void allometry_dave(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}

		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			// Call trunk allometry.
			allometry(indiv);
		} else {
			// TODO: Assaf
			allometry_dave_tree(indiv);
		}
	}
}

// Aggregate fluxes for annual output.
void annual_accounting_dave(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			// This pft is not using dave allocation. This could happen, for
			// example, if it's a tree and ifdailytree is false.
			continue;
		}

		// Increase age by one day
		indiv.age += 1.0 / Date::MAX_YEAR_LENGTH;

		// Initiate yearly maximum
		if (date.day == 0) {
			indiv.cmass_leaf_ygrowth = 0.0;
			indiv.cmass_root_ygrowth = 0.0;
			indiv.lai_ymax = 0.0;
			indiv.fpc_ymax = 0.0;
			indiv.cmass_leaf_ymax = 0.0;
			indiv.nmass_leaf_ymax = 0.0;
			indiv.cmass_root_ymax = 0.0;
			indiv.nmass_root_ymax = 0.0;

			// scale indiv for landcover change
			// scale_indiv(indiv, false);  // TODO: Is this needed?
		}

		// Determine the yearly maximum values which is reported
		indiv.lai_ymax = max(indiv.lai_ymax, indiv.lai);
		indiv.fpc_ymax = max(indiv.fpc_ymax, indiv.fpc);

		if (indiv.get_cmass_leaf() > indiv.cmass_leaf_ymax) {
			indiv.cmass_leaf_ymax = indiv.get_cmass_leaf();
			indiv.nmass_leaf_ymax = indiv.get_nmass_leaf();
		}

		if (indiv.cmass_root > indiv.cmass_root_ymax) {
			indiv.cmass_root_ymax = indiv.cmass_root;
			indiv.nmass_root_ymax = indiv.nmass_root;
		}

		// Report yearly data
		if (date.islastday && date.islastmonth) {

			// harvest_year(indiv); // called once a year --> to be replaced in DAVE

			// Calculate vegetation carbon and nitrogen mass and average C:N ratios for yearly outputs
			// TODO: why report maximum and not average?
			indiv.cmass_veg = indiv.cmass_leaf_ymax + indiv.cmass_root_ymax;
			indiv.nmass_veg = indiv.nmass_leaf_ymax + indiv.nmass_root_ymax;

			// Save leaf annual average C:N ratio
			if (!negligible(indiv.nday_leafon)) {
				indiv.cton_leaf_aavr /= indiv.nday_leafon;
			}
			else {
				indiv.cton_leaf_aavr = indiv.pft.cton_leaf_max;
			}

			indiv.deltafpc = 0.0;
		} // Last day of year
	}
}

void growth_daily_accounting_dave(Vegetation &vegetation) {
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		// Set indiv.cgrow, indiv.cexcess, and indiv.nexcess to 0.
		indiv.cgrow           = 0.0;
		indiv.cexcess         = 0.0;
		indiv.nexcess         = 0.0;
		indiv.nmass_missing   = 0.0;
		indiv.ctolitter       = 0.0;
		indiv.ntolitter       = 0.0;
		indiv.ctolitter_leaf  = 0.0;
		indiv.ntolitter_leaf  = 0.0;
		indiv.ctolitter_root  = 0.0;
		indiv.ntolitter_root  = 0.0;
		indiv.ctolitter_repr  = 0.0;
		indiv.ntolitter_repr  = 0.0;
		indiv.ctolitter_crown = 0.0;
		indiv.ntolitter_crown = 0.0;
	}
}

void growth_daily_grass_dave(Patch& patch) {


	// Obtain reference to Vegetation and Stand object for this patch
	Vegetation& vegetation = patch.vegetation;
	Stand& stand = patch.stand;

	// Record mass balance at start of day.
	DaveMassBalance mass_balance;
	mass_balance.init(vegetation);

	growth_daily_accounting_dave(vegetation);

		if (date.day == 0) {
			flush_cmass_repr(patch);
		}

	// External seed input
	seed_dispersal_dave(vegetation);

	// Calculate total daily growth and update C storage 
	storage_growth_partitioning_dave(vegetation);

	// Allocate available C for growth to different plant pools
	carbon_allocation_dave(vegetation);
	
	// Calculate N demand based on growth and stoichiometry.
	ndemand_dave(vegetation);

	// Calculate actual N uptake after taking competition among individuals into account
	nuptake_dave(vegetation);

	// update pools (C + N) to account for growth
	update_pools_dave(vegetation);
	mass_balance.check(vegetation);

	// Leaf senescence and retranslocation.
	senescence_dave(vegetation);
	mass_balance.check(vegetation);

	// Dead leaf and root to litter
	turnover_dave(vegetation);	
	mass_balance.check(vegetation);

	// C storage housekeeping and feedback to photosynthesis
	storage_dynamics_dave(vegetation);
	mass_balance.check(vegetation);

	// Allometry (Updates LAI and FPC)
	allometry_dave(vegetation);
	mass_balance.check(vegetation);

	// Management
	// Grazing and mowing
	grazing_daily_dave(vegetation);

	// yearly reporting. 
	// TODO: Shouldn't be here. Check if needed at all.
	annual_accounting_dave(vegetation);

	// Ensure that daily mass balance is closed.
	mass_balance.check(vegetation);

} // growth_daily_grass_dave() 


/// Handles daily crop allocation and daily lai calculation
/** Simple allocation based on heat unit accumulation.
 *  LAI is set directly after allocation from leaf carbon mass.
 */
void growth_daily(Patch& patch) {

	if (patch.stand.landcover == CROPLAND) {

		// allocate daily npp to leaf, roots and harvestable organs
		growth_daily_crop(patch);

		// update patchpft.lai_daily and fpc_daily
		lai_crop(patch);
	}

	// Daily C allocation.
	if (ifdavegrowth) {
		// Dave daily allocation.
		growth_daily_grass_dave(patch);
	} else if (ifdailygrass) {
		// David Warlind daily allocation.
		growth_daily_grass(patch);
	}

	// necessary updates after daily change to fpc in daily grass
	update_patch_fpc(patch);
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2001
// Huang, S, Titus, SJ & Wiens, DP (1992) Comparison of nonlinear height-diameter
//   functions for major Alberta tree species. Canadian Journal of Forest Research 22:
//   1297-1304
// Johnson, I. R., & Thornley, J. H. M. (1983). Vegetative crop growth model incorporating
//   leaf area expansion and senescence, and applied to grass. Plant Cell and Environment, 6(9)
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Prentice, IC, Sykes, MT & Cramer W (1993) A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Press, WH, Teukolsky, SA, Vetterling, WT & Flannery, BT. (1986) Numerical
//   Recipes in FORTRAN, 2nd ed. Cambridge University Press, Cambridge
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Shinozaki, K, Yoda, K, Hozumi, K & Kira, T (1964) A quantitative analysis of
//   plant form - the pipe model theory. I. basic analyses. Japanese Journal of
//   Ecology 14: 97-105
// Shinozaki, K, Yoda, K, Hozumi, K & Kira, T (1964) A quantitative analysis of
//   plant form - the pipe model theory. II. further evidence of the theory and
//   its application in forest ecology. Japanese Journal of Ecology 14: 133-139
// Sykes, MT, Prentice IC & Cramer W 1996 A bioclimatic model for the potential
//   distributions of north European tree species under present and future climates.
//   Journal of Biogeography 23: 209-233.
// Waring, RH Schroeder, PE & Oren, R (1982) Application of the pipe model theory
//   to predict canopy leaf area. Canadian Journal of Forest Research 12:
//   556-560
// Zaehle, S. & Friend, A. D. 2010. Carbon and nitrogen cycle dynamics in the O-CN
//   land surface model: 1. Model description, site-scale evaluation, and sensitivity
//   to parameter estimates. Global Biogeochemical Cycles, 24.
// Zeide, B. 2001. Natural thinning and environmental change: an ecological process model. 
// Forest Ecology and Management 154: 165-177.
