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
#include "guess.h"
#include "guessmath.h"
#include "landcover.h"
#include "ncompete.h"
#include <assert.h>
#include "somdynam.h"
#include "driver.h"
#include "logging.h"


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Maximum number of equivalent days with full leaf cover per growing season
// for summergreen PFTs
const double APHEN_MAX = 210.0;

/// External seed input (kgC m-2 d-1)
const double ext_seed = 10.0e-6;

// ---------------------------------
// Grass parameters
// ---------------------------------
// Maximum size of C storage pool as a fraction of biomass (pool-dependent!)
//const double cstorage_max_temp = 0.3;   // leaves and roots
//const double cstorage_max_perm = 0.75;  // crowns

// Fraction of maximum C storage not available for growth.
const double frac_cstorage_res = 0.15;

// Minimum physiological activity during dormancy (fraction)
const double fphen_activity_min = 0.1;

// SLA of culms relative to leaves in grasses
const double frac_sla_culm = 0.2;

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
	double nmass_uptake = indiv.pft.store_uptake() ? 0.0 : indiv.nmass_uptake;
	return indiv.nmass_storage +
		   indiv.nmass_leaf +
		   indiv.nmass_root +
		   indiv.nmass_sap +
		   indiv.nmass_repr +
		   indiv.nmass_heart +
		   nmass_uptake;
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
	return indiv.nmass_missing - indiv.ntolitter - indiv.ntoavail;
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


double wscal_growth_trees(double wscal) {
	double b = 0.1;
	double wscal_o = 0.5;
	const double THRESHOLD = 1e-4; 

  // return 1 / (1 + exp( (wscal_o - wscal) / b) ); // OLD
  double wscal_growth = ((1 + exp( (wscal_o - 1) / b) ) / (1 + exp( (wscal_o - wscal) / b) ));
  
  // force wscal_growth to be 0 when very small
  if (fabs(wscal_growth) < THRESHOLD) { wscal_growth = 0; }
  
//   wscal_growth = min(1.0, (1.0 / 0.2) * wscal);
//   wscal_growth = 1.0;

  return wscal_growth;
}

// Fill fraction of C storage (0-1)
double Cstorage_fillfrac(Individual& indiv) {
	return max(0.0, (indiv.cmass_storage - (frac_cstorage_res * indiv.cstorage_max)) / (indiv.cstorage_max - frac_cstorage_res * indiv.cstorage_max));
}

// Fill fraction of N storage (0-1)
double Nstorage_fillfrac(const Individual& indiv) {
	return max(0.0, (indiv.nmass_storage / indiv.nstorage_max));
}


// Temperature scalar used for growth and dormancy
// Relationship and parameters from Parent & Tardieu 2012, New Phytologist
// representative species:
// 		C3s: wheat (similar to barley, A. thaliana, canola)
// 		C4s: sorghum (similar to millet, maize)
// The function includes a modification at low temperatures.
// Purpose is to shut down physiological functioning below the dormancy threshold
// We need this function because common temperature response curves are not measured 
// at these temperatures and applying the function would likely overestimate activity.
// from Fatichi et al. 2014 (New Phytologist): 
// "a 5–6°C threshold has often been identified to limit growth in cold adapted
// species, irrespective of photosynthetic activity which typically
// ceases only at freezing point (Koerner, 2008)"
double temp_growth(double temp, double Ha_growth, double T0_growth, double A_growth, double T_dormancy) {
	
	double alpha   = 3.5;                   // unitless
	double Ha      = Ha_growth * J_PER_KJ;  // J mol-1
	double T0      = T0_growth + K2degC;    // K
	double Tleaf   = temp + K2degC;         // K   

	// ASSAF : TEMPORATY PARAMETERS FOR TESTING. AIM: TO GET TREES TO STORE IN WINTER BY
	// GETTING ftemp TO BE 0.
	// T0 = 30+ K2degC;

    // Note: A_growth was rescaled so that temp_growth ranges from 0-1
	double ftemp;
	double k = 1.0;
	double offset = 2.0;

	ftemp = A_growth * (Tleaf * exp(-Ha / (R_gas * Tleaf))) / 
			(1.0 + pow(exp(-Ha / (R_gas * Tleaf)),(alpha * (1.0 - Tleaf / T0))) );

	// attenuation at cold temperatures
	ftemp /= (1.0 + exp(-(temp - (T_dormancy - offset)) / k));	

	return max(0.0,min(1.0,ftemp));
}

// Water stress scalar used for growth and dormancy
double wscal_growth(double wscal, double wscal0, double b) {
	return (1 + exp((wscal0 - 1) / b)) / (1 + exp((wscal0 - wscal) / b));
	//double wscal_min = 0.4;  // no growth under severe water stress (ideally consistent with pft.wscal_min)
	//return min(1.0,max(0.0,pow(((wscal - wscal_min) / (1.0 - wscal_min)),q)));
}

// Photoperiod affecting development and thus growth demand
// Effects likely limited to high latitude areas.
double photo_growth(double reldaylength, double daylength_diff) {
	return 1.0;
}

// Functions calculating traits based on the environment (i.e. making them adaptive)
double calc_dphen_activity_max(double dphen_activity_base, double dphen_activity_pot,
	                           double k, size_t window_length, size_t window_length_min){
	return dphen_activity_base + (dphen_activity_pot - dphen_activity_base) * exp(-k * (window_length - window_length_min));
}

double calc_max_mobil_cstorage(double max_mobil_base, double max_mobil_pot,
	double k, size_t window_length, size_t window_length_min){
	return max_mobil_base + (max_mobil_pot - max_mobil_base) * exp(-k * (window_length - window_length_min));
}

double calc_temp_dormancy(double temp_dormancy_base, double mtemp_min20){
	return temp_dormancy_base + min(0.0, max(-3.0,mtemp_min20 * 0.3));
}


// Calculates dphen_activity, the change of metabolic activity related to dormancy
double dphen_activity(Individual &indiv, Climate& climate) {

	Pft& pft = indiv.pft;
	Patchpft &ppft = indiv.patchpft();
	
	// higher water stress thresholds for dormancy and vegetative regrowth
	// compared to the vegetative phase.
	double wscal_min = pft.wscal_min_dormant;
	if (indiv.growth_stage == vegetative && is_zero(1.0 - indiv.fphen_activity) && indiv.fdev_growth > 0.1) {
		wscal_min = pft.wscal_min;
	}

	const Historic<double,10> &wscal_10 = indiv.get_wscal_10();

	// Threshold temperature for dormancy (base temperature same as dormancy threshold, depends on photosyn. pathway)
	//double grass_temp = indiv.vegetation.patch.soil.T_soil[7];   // Temperature of first soil layer
	double grass_temp = climate.temp;
	double wscal_mean = indiv.get_wscal_moving_average();
	double tscal_crit = temp_growth(indiv.temp_dormancy,indiv.pft.Ha_growth,indiv.pft.T0_growth,indiv.pft.A_growth,indiv.temp_dormancy);

	//double wscal_actual = wscal_growth(wscal_mean);
	double wscal_actual = wscal_mean;
	double tscal_actual = temp_growth(indiv.soiltemp_dynamic.mean(), indiv.pft.Ha_growth, indiv.pft.T0_growth, indiv.pft.A_growth,indiv.temp_dormancy);

	// Changes in physiological activity due to changes in dormancy status
	double denom_wscal = wscal_actual >= wscal_min ? pow((1.0 - wscal_min), 1.5) : wscal_min;
	double num_wscal   = wscal_actual >= wscal_min ? pow((wscal_actual - wscal_min), 1.5) : (wscal_actual - wscal_min);

	double tscal = tscal_actual > tscal_crit ? tscal_actual : (tscal_actual - tscal_crit) / tscal_crit;
	tscal = tscal > 0.0 ? pow(tscal, 0.25) : tscal;
	double wscal = num_wscal / denom_wscal;

	assert(tscal >= -1.0 && tscal <= 1.0);
	assert(wscal >= -1.0 && wscal <= 1.0);

	// avoid greening just before winter in temperature-limited ecosystems.
	// TODO: interactions with tscal and pscal?
	double pscal = 1.0;

	// (climate.mtemp_min20 < (indiv.pft.temp_dormancy - 5.0) &&
	//  ppft.daylength_10.periodictrend(1) < 0.0 && indiv.growth_stage != dormant) {
			
	//double daylength_crit = 0.9;
	//double daylength_min  = 0.6;

	//pscal = 2.0 * min(1.0,max(0.0,(ppft.reldaylength - daylength_min) / 
	//	    				      (daylength_crit - daylength_min))) - 1.0;
	//

	return min(min(tscal,pscal),wscal) * indiv.dphen_activity_max;

}


// Development effects on growth demand.
// Reduction of growth demand in reproductive phase.
void dev_growth(Individual &indiv, Climate &climate, double frepr_cstorage) {
	
	Patchpft &ppft = indiv.patchpft();

	if (indiv.growth_stage == reproductive || indiv.growth_stage == dormant) {

		double k_repr = 12.0;
		double seed_target = k_repr * indiv.pft.reprfrac * indiv.cmass_crown_flowering;

		indiv.fdev_growth = max(0.0, min(1.0, (1.0 - indiv.cmass_repr/seed_target)));

	} else if (indiv.growth_stage == vegetative && indiv.fdev_growth < 1.0) {

		// Phase after switching from reproductive period back to vegetative period.
        // check if conditions are suitable to resume vegetative growth.
		// If not, maintain reduced growth rate or enter dormancy.

		indiv.fdev_growth = min(1.0,max(fdev_growth_min * indiv.fphen_activity,indiv.fdev_growth * (1.0 + indiv.dphen_activity)));
	}
}

double move_toward(double current, double target, double max_change) {
	// Change in response unit, in absolute number of days.
	double change = fabs(target - current);

	// Limit change to change_rate.
	if (change > max_change) {
		change = max_change;
	}

	double sign = current >= target ? -1 : 1;
	return current + sign * change;
}

// Adjust window length for moving averages of soil temperature and
// water scalar.
void adjust_response_window(Individual &indiv, Climate &climate) {
	Patchpft &ppft = indiv.patchpft();

	// Only grasses respond dynamically to the growing conditions (for
	// now!). TODO: Assaf
	if (!indiv.pft.is_grass()) {
		return;
	}

	if (date.islastday && date.islastmonth) {
		double dormancy_depth = ppft.dormancy_depth_historical.mean();
		double wscal_avg = ppft.wscal_historical.mean();

		// Only allow the response window to change in size by 1 day per year.
		// TODO: parameterise this?
		const size_t MAX_CHANGE_RATE = 1;

		size_t response_window_new = get_response_window(dormancy_depth);
		size_t response_window_old = indiv.get_wscal_dynamic_capacity();

		// Move toward the new response window length by at most change_rate.
		double response_window = move_toward(response_window_old, response_window_new, MAX_CHANGE_RATE);
		assert(response_window > 0);

		indiv.set_wscal_dynamic_capacity(response_window);
		indiv.soiltemp_dynamic.set_capacity(response_window);

		// Adjust root beta factor (and therefore root distribution) based on average wscal.
		double root_beta_new = calc_root_beta(wscal_avg);
		const double MAX_ROOT_BETA_CHANGE = 0.0005;
		indiv.root_beta = move_toward(indiv.root_beta, root_beta_new, MAX_ROOT_BETA_CHANGE);
		calc_root_distribution(indiv.root_beta, indiv.rootdist);

		// Calculate dphen_activity_max (Maximum physiological up/down-regulation in and out of dormancy)
		// based on window length (quicker responses with shorter window lengths) as well as
		// max. C mobilisation rate (also faster with shorter window lengths).
		double k = 0.4;  // decline rate
        indiv.dphen_activity_max = calc_dphen_activity_max(indiv.pft.dphen_activity_max, 0.4,
			                                               k, response_window, 2);
		indiv.max_mobil_cstorage = calc_max_mobil_cstorage(indiv.pft.max_mobil_cstorage, 0.4,
														   k, response_window, 2);
											
		// lower dormancy threshold temperature for C3 grasses in cold environments
		if (indiv.pft.pathway == C3){  // only adjust for C3s for now. Could implement something similar for C4s.
			indiv.temp_dormancy = calc_temp_dormancy(indiv.pft.temp_dormancy,climate.mtemp_min20);
		} else {
			indiv.temp_dormancy = indiv.pft.temp_dormancy;
		}

	}
}

void state_accounting_dave(Vegetation &vegetation, Climate &climate) {
	// Keeping track of long-term stress factors and the response window are only required if the
	// adaptive root distribution for grasses is enabled.
	// Note: this wouldn't work with multiple cohorts at the moment!
	if (!ifadaptive_rootdist) {
		return;
	}

	for (size_t i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		Patchpft &ppft = indiv.patchpft();

		// Accumulate dormancy depth this year for this cohort as well as water stress.
		ppft.dormancy_depth += (1.0 - indiv.fphen_activity);
		ppft.wscal_cum += indiv.get_wscal();
		ppft.ndays_alive += 1;
	}

	// Add new value to moving average of dormancy depth and wscal, for each PFT.
	// Note: this is done for all PFTs, not just those currently established.
	if (date.islastday && date.islastmonth) {
		for (size_t i = 0; i < pftlist.nobj; i++) {
			Patchpft &ppft = vegetation.patch.pft[i];
			//double dormancy_depth_add = ppft.dormancy_depth + (Date::MAX_YEAR_LENGTH - ppft.ndays_alive);
			//double wscal_add = ppft.wscal_cum + (Date::MAX_YEAR_LENGTH - ppft.ndays_alive) * ppft.pft.wscal_min;
			if (ppft.ndays_alive < 90){
				ppft.dormancy_depth_historical.add(ppft.dormancy_depth_historical.mean());
				ppft.wscal_historical.add(ppft.wscal_historical.mean());
			} else {
				ppft.dormancy_depth_historical.add(ppft.dormancy_depth / ppft.ndays_alive);
				ppft.wscal_historical.add(ppft.wscal_cum / ppft.ndays_alive);
			}
		}
	}

	// Now adjust dynamic parameter accordingly, based on the long-term averages.
	for (size_t i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		Patchpft &ppft = indiv.patchpft();

		// Adjust window length for moving averages of soil temperature and
		// water scalar.
		adjust_response_window(indiv, climate);
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

		// Driving temperature for grass phenology (air vs soil temp)
		//indiv.soiltemp_dynamic.add(vegetation.patch.soil.T_soil[7]);  // temperature of first soil layer
		indiv.soiltemp_dynamic.add(climate.temp);

		if (pft.is_grass()) { 

			// 1) Dormancy
			if (indiv.pft.lifeform == GRASS_PERENNIAL) {

				indiv.dphen_activity = dphen_activity(indiv, climate);
				indiv.fphen_activity = min(1.0,max(fphen_activity_min,indiv.fphen_activity * (1.0 + indiv.dphen_activity)));

			} else {
			    // no dormancy in annuals
				indiv.fphen_activity = 1.0;
			}


			// 2) Development factor for growth (returns indiv.fdev_growth, which is used in
			//    growth_demand() function returning growth demand).

			// Ratio of reproductive biomass to storage capacity, which is fed into dev_growth
			indiv.frepr_cstorage = indiv.cmass_repr / max(indiv.cstorage_max,1e-9);
			dev_growth(indiv, climate, indiv.frepr_cstorage);
			

			// 3) Switches between growth stages (indiv.growth_stage):
			// 		0: dormant stage
			// 		1: vegetative stage
			// 		2: reproductive stage
			double repr_trend = 0.0;  // relative trend in repr. tissue

			// TODO: get rid of growth_stage_days when they aren't needed anymore.
		   	if (indiv.growth_stage == dormant && indiv.fphen_activity > fphen_activity_min) {
				indiv.growth_stage = vegetative;
				indiv.fdev_growth = 1.0;
				indiv.growth_stage_days = 0;
				indiv.lai_trend_max = 0.0;
			}
			else if (indiv.growth_stage == vegetative) {
				if (indiv.cmass_storage_20.size() > 5 && indiv.cmass_storage_20.periodictrend(20) > 0.0 &&
				    indiv.cmass_storage_20.periodictrend(5) > 0.0 &&
				    is_zero(1.0 - indiv.fdev_growth) && Cstorage_fillfrac(indiv) < 0.9 &&
				    indiv.growth_stage_days > 20 && indiv.lai_20.periodictrend(5) > 0.05 * indiv.lai_trend_max){  
					// last two conditions above are there to prevent too early flowering.
					if (indiv.lai > lai_min) {
						indiv.growth_stage = reproductive;
						indiv.growth_stage_days = 0;
						indiv.cmass_crown_flowering = indiv.cmass_crown;
					} else {
						indiv.kill();
						vegetation.killobj();
						continue;
					}
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
                }
 
                // flush reproductive tissue to seedbank.
                if (indiv.growth_stage != reproductive && indiv.cmass_repr > 0.0) {
                    indiv.flush_repr_litter_DAVE(1.0);
				}
			}
			indiv.growth_stage_days += 1;
		} 
		else {
			indiv.fphen_activity = 1.0;
		}

		vegetation.nextobj();
	}

	state_accounting_dave(vegetation, climate);
}

void leaf_phenology_pft(Pft& pft, Climate& climate, double wscal, double aphen,
	double& phen) {

	// if (pft.dave_allocation()) {
	// 	return; // TODO: rationalise
	// }

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
	const Climate &climate = indiv.vegetation.patch.get_climate();

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
				leaf_phenology_pft(ppft.pft, climate, ppft.wscal, ppft.aphen, ppft.phen); // note wscal is different in hydraulics version when calculated at PPFT versus indiv level

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

		if (indiv.pft.lifeform == TREE && ifdailytree) {
			// decide on phenological state (from Grote)
			tree_phenology_grote(indiv);
		}

		if (!indiv.pft.dave_allocation()) {
			// Update mlai. Better here than in npp(). No effect on dave indivs
			// because they will update this after calling allometry each day.
			if (indiv.lai_today() > indiv.mlai_max[date.month]) {
				indiv.mlai_max[date.month] = indiv.lai_today();
			}

			// This is the same as indiv.lai for grasses, but it includes phen for trees/mosses.
			indiv.mlai[date.month] += indiv.lai_today() / (double)date.ndaymonth[date.month];
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
	bool alive, double& cmass_sap_turnover) {

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
		cmass_sap_turnover = turnover;

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

	//assert(fabs(bminc + exceeds_cmass - (cmass_leaf_inc + cmass_root_inc + cmass_sap_inc + cmass_heart_inc + litter_leaf_inc + litter_root_inc)) < EPS);
}

double calc_creserve_opt_grote(const Individual &indiv) {

	double nsc_leaf = indiv.get_cmass_leaf() * indiv.pft.max_leaf_nsc_conc;
	double nsc_wood = indiv.cmass_sap * indiv.pft.max_wood_nsc_conc;

	return nsc_leaf + nsc_wood;
}

// Capacity of C storage (kg C m-2)
double cstorage_capacity(Individual& indiv) {

	if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
		return cstorage_max_temp * (indiv.get_cmass_leaf() + indiv.cmass_root) + cstorage_max_perm * indiv.cmass_crown;
	} else {

		// TODO
		if (ifgrote) {
			return calc_creserve_opt_grote(indiv);
		}

		// TODO: what about seedlings? cstorage_max cn be much higher than the total cmass of all the pools
		// const int nyear_turnover = 3;
		// double turnover_sap = indiv.cmass_sap * indiv.pft.turnover_sap;
		// double turnover_leaf = indiv.cmass_leaf * indiv.pft.turnover_leaf;
		// double turnover_root = indiv.cmass_root * indiv.pft.turnover_root;
		// double cstorage_opt = nyear_turnover * (turnover_sap + turnover_leaf + turnover_root);

		// return max(cstorage_opt, indiv.cstorage_max);

		// define the maximum storage capacity of the plant live carbon pool. 
		// TO DO: set as PFT parameters !!!!
		// Literature values vary significantly , depending on the type of measurement. 
		// For example, Smith et al. (2018) measured 4.4% of NSC from dry mass of
		//  E. Obliqua (14yo). See Quentin et al (2015) how the differences
		// in NSCs depends on the laboratory and measurement method (!)
		// Smith et al. (2018)  https://doi.org/10.1093/treephys/tpx141
		// Quentin et al. (2015) https://doi.org/10.1093/treephys/tpv073

		// const double cstorage_fraction_sap = 0.5; // old simulations 
		// const double cstorage_fraction_leaf = 0.2; // old simulations
		// const double cstorage_fraction_root = 0.2; // old simulations

		const double cstorage_fraction_sap = 0.27; // good simulation
		const double cstorage_fraction_leaf = 0.1; // good simulation
		const double cstorage_fraction_root = 0.1; // good simulation

		// maximum storage capacity
		return (indiv.cmass_sap * cstorage_fraction_sap) + 
			(indiv.cmass_leaf * cstorage_fraction_leaf) + 
			(indiv.cmass_root * cstorage_fraction_root);

		// return cstorage_max_temp * (indiv.get_cmass_leaf() + indiv.cmass_root) + cstorage_max_perm * indiv.cmass_sap;
	}

}

// Capacity of N storage (kg N m-2)
double nstorage_capacity(Individual& indiv) {

	return indiv.cstorage_max / indiv.pft.cton_root_min;

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
@param frac_resp: Fraction of growth to be respired.
@param resp_growth: (Output) Growth respiration.
*/
void constrain_growth(
		double &cgrow,
		double &ngrow,
		double cton_min,
		double cton_max,
		double &creturn,
		double &nreturn,
		double frac_resp,
		double &resp_growth) {
	assert(frac_resp >= 0.0);
	assert(frac_resp < 1.0);

	// When the upper/lower C:N limits are reached, only allow new growth to
	// move the C:N ratio to within this fraction of the limits. If we don't
	// have this limits, growth will move the C:N ratio exactly onto the limits,
	// and rounding errors can then push the C:N ratio beyond the limits.
	const double FGROWTH_MAX = 0.999;

	double cgrow_pot = cgrow;
	cgrow = cgrow_pot / (1.0 + frac_resp);
	resp_growth = cgrow_pot - cgrow;

	if ( ngrow < (cgrow / cton_max)) {
		// C:N ratio in new growth exceeds max C:N ratio.Limit C in new growth
		// to the max C:N ratio and put the rest into creturn.
		double cgrow_max = ngrow * cton_max * FGROWTH_MAX;
		double overflow = cgrow_pot - cgrow_max;

		// Growth respiration can be paid by the C we can't grow.
		resp_growth = cgrow_max * frac_resp;

		// We were already able to pay for the growth respiration with
		// cgrow_pot - cgrow, and now that difference has increased (because
		// cgrow has decreased).
		assert(overflow > resp_growth);

		overflow -= resp_growth;
		creturn += overflow;
		cgrow = cgrow_max;
	} else if ( ngrow > (cgrow / cton_min) ) {
		// N:C ratio in new growth exceeds max N:C ratio. Limit N in new growth
		// to the max N:C ratio and put the rest into nreturn.
		double ngrow_max = FGROWTH_MAX * cgrow / cton_min;
		nreturn += ngrow - ngrow_max;
		ngrow = ngrow_max;

		// cgrow is unchanged, so no need to adjust resp_growth.
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

			// TODO: do we want growth respiration during establishment? The
			// carbon is already accounted for via a magical establishment flux.
			// Should we then respire some of that?
			double resp;
			constrain_growth(
				cmass_leaf, nmass_leaf,
				indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
				cmass_storage, nmass_storage, 0, resp);

			constrain_growth(
				cmass_root, nmass_root,
				indiv.pft.cton_root_min, indiv.pft.cton_root_max,
				cmass_storage, nmass_storage, 0, resp);

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
			assert(indiv.cstorage_max > 0);
			assert(indiv.nstorage_max > 0);

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
		if (!indiv.pft.dave_allocation()) {
			indiv.initialn = indiv.get_nmass_leaf() + indiv.nmass_root + indiv.nmass_storage + indiv.nmass_crown;
		}

		double cmass_after = indiv.ccont();
		double nmass_after = indiv.ncont();
		assert(equal(cmass_before + cmass_init, cmass_after));
		assert(equal(nmass_before + nmass_init, nmass_after));
	}
	else {
		// trees (and mosses?)
		if (ifdavegrowth && ifdailytree) {
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

		if (ifforcedlai) {
			indiv.lai = fixed_lai;
		} else {
			indiv.lai = indiv.cmass_leaf * indiv.pft.sla;
		}

		if (indiv.pft.dave_allocation()) {
			// Capacity of C and N storage on establishment
			if (is_zero(indiv.cstorage_max)) {
				indiv.cstorage_max = cstorage_capacity(indiv);
				indiv.nstorage_max = nstorage_capacity(indiv);
			}
			assert(indiv.cstorage_max > 0);
			assert(indiv.nstorage_max > 0);

			// initialise limits
			if (is_zero(indiv.cmass_sap_limit)) {
				indiv.cmass_sap_limit = indiv.cmass_sap;
				indiv.cmass_leaf_limit = indiv.cmass_leaf;
				indiv.cmass_root_limit = indiv.cmass_root;
			}
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

typedef enum { FIXEDHEIGHT, FIXEDLAI, BMINIT } est_alloc_strategy_t;

// Initial height in m (for FIXEDHEIGHT strategy).
const double INIT_HEIGHT = 2;

// Initial LAI in m2/m2 (for FIXEDLAI strategy).
const double INIT_LAI = 1;

// Initial storage fraction of capacity.
const double INIT_STORAGE_FRACTION = 1;

double solve_csap(double bminit, const Pft &pft, double ltor, double init_storfrac, double EPS = 1e-6, size_t NITER_MAX = 10000) {

	// TODO: consider how to handle ltor=0.
	if (ltor == 0) {
		fail("solve_csap(): ltor cannot be zero.");
	}

	if (is_zero(bminit)) {
		return 0;
	}

	// Parameters.
	double k   = pft.k_latosa;
	double k2  = pft.k_allom2;
	double k3  = pft.k_allom3;
	double rho = pft.wooddens;
	double sla = pft.sla;
	double ks = init_storfrac;

	double lowbound = 0;
	double highbound = bminit;

	for (int i = 0; i < NITER_MAX; i++) {
		double csap = (lowbound + highbound) / 2.0;
		double height = pow(2.0 * csap * pow(k2, 2.0 / k3) / (PI * rho), 1.0 / (1.0 + 2.0 / k3));
		double diam   = pow(height / k2, 1.0 / k3);
		// v = h * pi * r^2 = h * pi * d^2 / 4
		double cleaf  = k * csap / (height * rho * sla);

		double lai = cleaf * sla;
		double sa  = csap / (height * rho);
		double latosa = lai / sa;

		double croot  = cleaf / ltor;
		double cstor  = ks * (cleaf / 5.0 + csap / 2.0 + croot / 2.0);

		double bmtot  = cleaf + croot + csap + cstor;
		double cerr   = bminit - bmtot;

		// If cerr is positive, we have allocated too little biomass.
		if (cerr > 0) {
			lowbound = csap;
		} else {
			// If cerr is negative, we have allocated too much biomass.
			highbound = csap;
		}

		if (fabs(cerr) < EPS) {
			return csap;
		}
	}
	xtring msg;
	msg.printf("solve_csap(): Failed to find a solution for csap;\nbminit=%f, k_latosa=%f, k_allom2=%f, k_allom3=%f, wooddens=%f, sla=%f, init_storfrac=%f, ltor=%f",
	           bminit, k, k2, k3, rho, sla, ks, ltor);
	fail(msg);
	return -1; // Can never reach here: fail() calls exit().
}

double get_la2sa(const Individual &indiv) {
	// double area = indiv.diam * indiv.diam * PI / 4.0;
	// double area_heart = indiv.cmass_heart / indiv.pft.wooddens;
	// double area_sap = area - area_heart;
	double sa = indiv.cmass_sap / (indiv.height * indiv.pft.wooddens);
	return indiv.lai / sa;
}

void establishment_sapling_dave_dynamic(double bminit, Individual &indiv, double ltor) {
	// Constants.
	const Pft &pft = indiv.pft;
	double k   = pft.k_latosa;
	double k2  = pft.k_allom2;
	double k3  = pft.k_allom3;
	double rho = pft.wooddens;
	double sla = pft.sla;
	double ks  = INIT_STORAGE_FRACTION;

	// Maximum number of iterations to be used for solving for csap.
	const size_t IT_MAX = 8192;
	const double TOL = 1e-6;
	const double EPS = TOL;

	bminit *= indiv.densindiv; // convert to kgC/m^2

	double csap = solve_csap(bminit, indiv.pft, ltor, ks);

	// Eqn 12 in derivation.
	double height = pow(4.0 * (csap / indiv.densindiv) * pow(k2, 2.0 / k3) / (PI * rho), 1.0 / (1.0 + 2.0 / k3));
	double cleaf = k * csap / (height * rho * sla);
	double croot = cleaf / ltor;
	double cstor = ks * (cleaf / 5.0 + csap / 2.0 + croot / 2.0);

	indiv.cmass_leaf    = cleaf;
	indiv.cmass_root    = croot;
	indiv.cmass_sap     = csap;
	indiv.cmass_storage = cstor;
	indiv.cstorage_max  = cstor;
	indiv.meanleafage   = pft.turnover_leaf; // assuming sapling is 1 year old
	indiv.meansapage   = pft.turnover_sap; // assuming sapling is 1 year old
	indiv.meanrootage   = pft.turnover_root; // assuming sapling is 1 year old

	indiv.nmass_leaf    = cleaf / pft.cton_leaf_min;
	indiv.nmass_root    = croot / pft.cton_root_min;
	indiv.nmass_sap     = csap  / pft.cton_sap_min;
	indiv.nstorage_max  = nstorage_capacity(indiv);
	indiv.nmass_storage = ks * indiv.nstorage_max;

	// This is important for accounting purposes.
	indiv.initialn = indiv.nmass_sap + indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_storage;

	indiv.lai = indiv.cmass_leaf * pft.sla;
	indiv.height = height;
	indiv.diam = pow(indiv.height / k2, 1.0 / k3);
	double latosa = get_la2sa(indiv);

	assert(equal(bminit, cleaf + croot + csap + cstor, 0.1));
	assert(equal(latosa, pft.k_latosa));
	assert(equal(PI * indiv.diam * indiv.diam * 0.25 * indiv.height * indiv.pft.wooddens * indiv.densindiv, indiv.cmass_sap));
}

void establishment_sapling_dave_fixedlai(double bminit, Individual &indiv, double ltor) {
	// lai = cmass_leaf * sla
	indiv.cmass_leaf = INIT_LAI / indiv.pft.sla;
	indiv.cmass_root = indiv.cmass_leaf / ltor;
	// leaf = sap * pft.k_latosa / (indiv.pft.sla * pft.wooddens * indiv.height)
	// h = k2 * (4 * v / pi) ^ (1 / (1 + 2 / k3))
	// v = sap / rho

	// leaf * sla = sap / (rho * h) * klatosa
	// h = sap * klatosa / (leaf * sla * rho)
	/*
    leaf = sap * klatosa / (sla * rho * k2 * (4 * sap / (rho * pi)) ^ (1 / (1 + 2 / k3)))
	*/
	indiv.height = 1.0;
	// indiv.cmass_sap = solve_sap(indiv.cmass_leaf, 1e-6, 1e-6, 1e6, indiv.cmass_leaf, indiv.pft.sla, indiv.pft.wooddens, indiv.pft.k_allom2, indiv.pft.k_allom3, indiv.pft.k_latosa);
	indiv.cmass_sap = indiv.cmass_leaf * indiv.pft.sla * indiv.pft.wooddens * indiv.height / indiv.pft.k_latosa;
	// double volume = indiv.cmass_sap / indiv.pft.wooddens;
	indiv.diam = pow(indiv.height / indiv.pft.k_allom2, 1 / indiv.pft.k_allom3);
	indiv.lai = indiv.cmass_leaf * indiv.pft.sla;


	// B = cs(1 + k / (hp.SLA) * (ks/5 + ks/2lr + 1 + 1/lr))

	// Update capacity of C and N storage
	indiv.cstorage_max = cstorage_capacity(indiv);
	indiv.nstorage_max = nstorage_capacity(indiv);
	indiv.cmass_storage = indiv.cstorage_max;

	assert(indiv.cstorage_max > 0);
	assert(indiv.nstorage_max > 0);

	// Nitrogen.
	indiv.nmass_leaf = indiv.cmass_leaf / indiv.pft.cton_leaf_min;
	indiv.nmass_sap = indiv.cmass_sap / indiv.pft.cton_sap_min;
	indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_min;
	double nstorage_cap = nstorage_capacity(indiv);
	indiv.nmass_storage = nstorage_cap * INIT_STORAGE_FRACTION;
	indiv.crownarea = indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp);
	indiv.crownarea = min(indiv.crownarea, indiv.pft.crownarea_max);
	indiv.initialn = indiv.nmass_sap + indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_storage;
}

void establishment_sapling_dave(double bminit, Individual& indiv, double ltor) {

	const est_alloc_strategy_t strategy = BMINIT;
	// const est_alloc_strategy_t strategy = FIXEDHEIGHT;

	indiv.alive = true;
	assert(bminit > 0);

	// ##### hypothetical allocation

	if (strategy == BMINIT) {
		establishment_sapling_dave_dynamic(bminit, indiv, ltor);
		return;
	}

	if (strategy == FIXEDLAI) {
		establishment_sapling_dave_fixedlai(bminit, indiv, ltor);
		return;
	}

	// ##### end hypothetical allocation

	// initial height (m)

	indiv.height = INIT_HEIGHT;

	indiv.diam = pow(indiv.height / indiv.pft.k_allom2, 1 / indiv.pft.k_allom3);
	
	// C for new sapling
	indiv.cmass_sap = indiv.height * indiv.diam * indiv.diam / 4 * PI * indiv.pft.wooddens * indiv.densindiv;

	double fns = indiv.pft.k_latosa / (indiv.pft.sla * indiv.pft.wooddens * indiv.height);

	indiv.cmass_leaf = indiv.cmass_sap * fns;

	indiv.cmass_root = indiv.cmass_leaf / ltor;

	indiv.cstorage_max = cstorage_capacity(indiv);
	assert(indiv.cstorage_max > 0);

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
		establishment_sapling_dave(bminit, indiv, ltor);
	} else {
		establishment_biomass(bminit, ltor, indiv);
	}


	if ( (indiv.pft.is_grass() && ifdailygrass) || indiv.pft.dave_allocation()) {

		indiv.report_flux(Fluxes::ESTC, -(indiv.get_cmass_leaf() + indiv.cmass_root + indiv.cmass_storage + indiv.cmass_crown + indiv.cmass_sap));
		indiv.report_flux(Fluxes::ESTN, -indiv.initialn);
	}
}

double calc_bole_height(const Individual &indiv) {
	// boleht_frac should be 0 for trunk runs.
	return indiv.pft.boleht_frac * indiv.height;
}

double calc_rpc_tree(Individual &indiv) {
	return 4.0 * indiv.fpc;
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
	double lai_leaf=0.0;        // grass only
	double lai_culm=0.0;        // grass only
	double lai_leaf_brown=0.0;  // grass only
	double lai_culm_brown=0.0;  // grass only

	// guess2008 - max tree height allowed (metre).
	const double HEIGHT_MAX = 150.0;


	if (indiv.pft.lifeform == TREE) {

		if (indiv.pft.dave_allocation()) {
			
			double cgrow_sap = indiv.cgrow_sap;
			
			indiv.cgrow_sap = 0;

			// redirect to dave tree allomtery
			allometry_dave_tree(indiv);
			
			indiv.cgrow_sap = cgrow_sap;

			return true;
		}

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
			indiv.rpc = calc_rpc_tree(indiv);  // TODO Assaf: decide if rpc is needed for trees. RPC is used in Nuptake function.
		}
		else {
			indiv.lai_indiv = 0.0;
			indiv.fpc = 0.0;
		}

		// Bole height (Eqn 7)
		indiv.boleht = calc_bole_height(indiv);

		// Stand-level LAI
		if (ifforcedlai) {
			indiv.lai = fixed_lai;
		} else {
			indiv.lai = indiv.cmass_leaf * indiv.sla();
		}
	}
	else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {

		// GRASSES AND MOSSES

		if (indiv.pft.landcover != CROPLAND) {

			// guess2008 - bugfix - added if
			if (!negligible(indiv.get_cmass_leaf())) {

				// Grass "individual" LAI (Eqn 11).
				indiv.lai_indiv = indiv.get_cmass_leaf() * indiv.sla();
				if (ifdavegrowth) {
					lai_leaf = indiv.get_cmass_leaf() * (1.0 - indiv.pft.culmfrac) * indiv.sla();
					lai_culm = indiv.get_cmass_leaf() * indiv.pft.culmfrac * (indiv.sla() * frac_sla_culm);
					lai_leaf_brown = indiv.get_cmass_leaf_brown() * (1.0 - indiv.pft.culmfrac) * indiv.sla();
					lai_culm_brown = indiv.get_cmass_leaf_brown() * indiv.pft.culmfrac * (indiv.sla() * frac_sla_culm);
					indiv.lai_indiv = lai_leaf + lai_culm;
					indiv.lai_brown_indiv = lai_leaf_brown + lai_culm_brown;
					indiv.lai_20.add(indiv.lai_indiv);
					if (indiv.lai_20.periodictrend(5) > indiv.lai_trend_max){
						indiv.lai_trend_max = indiv.lai_20.periodictrend(5);
					}
					// indiv.lai_indiv *= (1.0 - indiv.pft.culmfrac);
				}
				
				// FPC (Eqn 10)
				indiv.fpc = 1.0 - lambertbeer(indiv.lai_indiv);
				indiv.fpc_greenbrown = 1.0 - lambertbeer(indiv.lai_indiv + indiv.lai_brown_indiv);
				// Root projective cover (RPC) introduced to capture feedback between 
				// belowground C investment and available N. eff_rpc represents the efficiency of 
				// roots to cover volume, related to specfic root length (SRL in m g-1). SRL is ca. 4 times 
				// higher in annuals than perennials (Roumet et al. 2006, New Phytologist).
				double eff_rpc = indiv.pft.lifeform == GRASS_PERENNIAL ? 6.0 : 20.0; 
				indiv.rpc = 1.0 - lambertbeer(indiv.lai_indiv * min(10.0,(eff_rpc * indiv.cmass_root / indiv.cmass_leafs[green])));
				// TODO: Replace with better solution. Does this work with competition??

				// Stand-level LAI
				if (ifforcedlai) {
					indiv.lai = fixed_lai;
				} else {
					indiv.lai = indiv.lai_indiv;
				}
			}
			else {
				indiv.lai_indiv = 0.0;
				indiv.fpc = 0.0;
				if (ifforcedlai) {
					indiv.lai = fixed_lai;
				} else {
					indiv.lai = 0.0;
				}
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
// is disabled. Called from growth_daily_dave() on 1st day of year
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

		indiv.cmass_sap_yesterday = indiv.cmass_sap;

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
			double cton_leaf_aopt = max(indiv.cton_leaf_aopt, indiv.pft.cton_leaf_avr);

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

				indiv.acalloc_repr = cmass_repr;
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
							(indiv.pft.leaflong() * (double)date.year_length()) * (indiv.cmass_leaf + indiv.cmass_root) -
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
								indiv.alive, indiv.cmass_sap_turnover);
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
						indiv.acalloc_leaf = cmass_leaf_inc * indiv.densindiv;

						// Roots
						indiv.cmass_root += cmass_root_inc * indiv.densindiv;
						indiv.acalloc_root = cmass_root_inc * indiv.densindiv;

						// Sapwood
						indiv.cmass_sap += cmass_sap_inc * indiv.densindiv;
						indiv.acalloc_sap = cmass_sap_inc * indiv.densindiv;

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

							log_info("%s killed by carbon starvation...", (const char *)indiv.pft.name);
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

							log_info("%s killed by carbon starvation...", (const char *)indiv.pft.name);
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

					log_info("%s killed by carbon starvation...", (const char *)indiv.pft.name);
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
						indiv.flush_repr_litter_DAVE(1.0);
					}

					log_info("%s killed by carbon starvation...", (const char *)indiv.pft.name);
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
void growth_daily_grass_warlind(Patch& patch) {

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
					indiv.nscal_30.add(min(1.0, indiv.cton_leaf_opt / indiv.cton_leaf()));
				}
				
				// leaf to root ratio of new growth (Johnson 2008 eqn 4.13b)
				//indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * min(indiv.wscal_running.mean(), indiv.nscal_running.mean()) * ltor_status;
				indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * min(indiv.get_wscal_5().mean(), indiv.nscal_5.mean()) * ltor_status;
			}
			else {
				indiv.ltor = ltor_min + (indiv.pft.ltor_max - ltor_min) * indiv.get_wscal_5().mean() * ltor_status;
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
			indiv.wscal_fac = max(0.0, (indiv.get_wscal_5().mean() - WSCAL_MIN) / (1.0 - WSCAL_MIN));
			
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

	// Sensitivity of water stress to growth
	// growth is more sensitive to water stress than photosynthesis (see Fatichi et al. 2014)

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
	//double cmass_total  = indiv.cmass_root + indiv.cmass_crown + indiv.get_cmass_leaf() + indiv.cmass_storage;
	double cmass_total  = indiv.cmass_root + indiv.cmass_storage + indiv.get_cmass_leaf();

	// Fill status of C and N storage
    double cstorage_fill = Cstorage_fillfrac(indiv);
    double nstorage_fill = Nstorage_fillfrac(indiv);


    // ------------------------------
	// Calculate growth demand
	// ------------------------------
	// Factors determining growth demand
	fwater_growth = wscal_growth(indiv.get_wscal(), wscal0, bgrowth);
	ftemp_growth  = temp_growth(grass_temp, pft.Ha_growth, pft.T0_growth, pft.A_growth,indiv.temp_dormancy);
	fphoto_growth = photo_growth(ppft.reldaylength, ppft.daylength_10.periodictrend(1));  // set to 1.0 
	//dev_growth(indiv, indiv.frepr_cstruct); // indiv.fdev_growth now calculated in phenology_grass_dave beforehand!

	growth_demand = indiv.fdev_growth * ftemp_growth * fwater_growth * fphoto_growth;

	// feedbacks with C storage (reduce growth demand when NSCs are low)
	double qf1 = indiv.growth_stage == vegetative ? 0.1 : 1.1 - indiv.fdev_growth;
	growth_demand *= pow(cstorage_fill,qf1);

	// reduced growth under dormancy (total cessation in full dormancy)
	growth_demand *= (indiv.fphen_activity - fphen_activity_min) / (1.0 - fphen_activity_min);

	return growth_demand;

}



void storage_growth_partitioning_grass(Individual &indiv) {
	// Local:
	double cstore = 0.0;

	// C storage not available for growth.
	double cstorage_res = frac_cstorage_res * indiv.cstorage_max;

	// C source: Potential growth rate (kgC m-2 d-1)
	double Gpot = max(0.0, indiv.max_mobil_cstorage * max(0.0,(indiv.cmass_storage - cstorage_res)) + indiv.dnpp); 

	// C sink: growth demand (0-1).
	indiv.growth_demand = growth_demand(indiv);

	if (indiv.dnpp >= 0.0) {
		
		// Positive NPP
		indiv.cgrow = Gpot * indiv.growth_demand;

		cstore = indiv.dnpp - indiv.cgrow;  // = Gpot - growth - (Gpot - NPP) 
											// cstore = re-allocation if negative

		// limit allocation to C storage if it is almost at capacity
		double crefill_max = max(0.0, (indiv.cstorage_max - indiv.cmass_storage) * indiv.max_mobil_cstorage);
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

			indiv.cgrow = Gpot * indiv.growth_demand;
			cstore = indiv.dnpp - indiv.cgrow;
		}
		else {
			log_info("%s killed by carbon starvation: negative NPP cannot be fulfilled by storage. NPP: %.9g kgC/m2, cmass_storage: %.9g kgC/m2",
				    (const char *)indiv.pft.name, indiv.dnpp, indiv.cmass_storage);
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

// A new function usedin Drake et al (2019) from June, Evans and Farquhar, 2004.
// The eqn is:
// k(T) = kTopt * exp -( (T - Topt) / omega )^2
// where Topt is the optimal temperature for growth (˚c), kTopt is the maximum k
// at Topt, and Omega, which is the temperature where k = e(-1)*Topt.
// The only change to the equation is the removal of kTopt, which makes this equation
// scale between 0 - 1. 
double temp_growth_JEF(double Temp) {
	// TODO: make PFT parameters and possibly allow to change (acllimation).
	// const double Topt = 28.3; // CUP measured Optimal temperature (degC)
	// const double Omega = 8.3; // CUP measured temperature where k = e(-1)*Topt
	const double Topt = 17.5; // WOMBAT (for E. obliqua; Hinko-Najera et al. 2019)
	// const double Omega = 6.0; // WOMBAT (Guestimated from Figure 3)
	const double Omega = 10.0; // WARRA (Guestimated)
	// const double Topt = 28.3; // CUP measured (for E. obliqua; Hinko-Najera et al. 2019)
	// const double Omega = 8.0; // CUP measured (Guestimated from Figure 3)

	return exp ( - pow( ((Temp - Topt) / Omega ) , 2));
	// return 1.0;
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
	double const inflection_pie_saver = 0.67; // Saver Hoch and Korner 2003 (https://doi.org/10.1046/j.0016-8025.2003.01032.x)
	double const inflection_pie_int = 0.5; // intermediate
	double const inflection_pie_spender = 0.35; // Spender
	// maximum daily mobilisation rate of c in the storage. Was 0.1
	double max_mobil_cstorage_tree = 0.2;
	// fraction of cmass_storage that is not available to growth
	const double frac_cstorage_res_tree = 0.15; // ASSAF BLARP
	// actual amount that is not available to growth (kgC/m^2)
	double cstorage_res = frac_cstorage_res_tree * indiv.cstorage_max;
	// amount of c available for growth today (kgC/m^2)
	double cmass_avail = max(0.0, (indiv.cmass_storage - cstorage_res));
	
	double store_state = 0;
	// sigmoid function that will change depending on the saver-spender spectrum of the pft.
	double fstorage = 0;
	if (indiv.cstorage_max > 0) {
		
		store_state = min(indiv.cmass_storage / indiv.cstorage_max, 1.0);
		
		fstorage = 1 - ((1 / (1 + exp(shape_pie * (store_state - inflection_pie_int)))));
	}
	// fstorage = 1.0;
	// temperature scalar for decoupling potosynthesis from growth (see Figure 9 in  Korner's 2-20 book)
	// double ftemp = temp_growth(indiv.patchpft().airtemp_10.mean(), 
	// 		indiv.pft.Ha_growth, indiv.pft.T0_growth, indiv.pft.A_growth, 
	// 		indiv.pft.temp_dormancy);
	double ftemp = temp_growth_JEF(indiv.patchpft().airtemp_10.mean());
	assert(ftemp <= 1);
	assert(ftemp >= 0);
	// water scalar for decoupling potosynthesis from growth (Fatichi et al 2014)
	double fwater = wscal_growth_trees(indiv.get_wscal_10().mean());
	// fwater = max(1.0, 0.5);
	// double fwater = wscal_10.mean();
	// daylength scalar (currently off)
	// double fphoto = climate.par / climate.par_max;
	double fphoto = 1.0;
	// in the first year, it comes up as negative
	if (date.year > 0) {
		// fphoto = (climate.daylength - climate.daylength_min) / (climate.daylength_max - climate.daylength_min);
		fphoto = climate.daylength  / climate.daylength_max ;
	}
	
	// TO DO: trying to implement storage strategy, where more storage 
	// (over growth) is encountered in some species (ie "savers").
	// NOT ENCOURAGE FILLING OF STORAGE, BUT RATHER EXCERBRATE IT!)
	// double xthreshold = 0.7;
	// double slope = 1/(1-xthreshold);
	// double yintercept = -slope * xthreshold; // solving for the intercept
	// fstorage = max(slope * store_state + yintercept, 0.0);
	// fstorage = store_state; // TEST THIS - LINEAR FSTORAGE
	// fstorage = 1.0; // TEST THIS - NO EFFECT OF STORAGE (REASON: LOW STORAGE WILL

	// This one is trying to rescale the maximum daily mobilisation rate of c 
	// in the storage using the fstorage, where more is allowed to be taken out 
	// of the storage when storage is more full.
	max_mobil_cstorage_tree *= fstorage; // max will be 0.1
	double cgrow = max_mobil_cstorage_tree * cmass_avail;
	// double cgrow = fstorage * max_mobil_cstorage_tree * cmass_avail;
	
	// cgrow potential which depends on environmental conditions(unitless). the higher 
	// the value, higher proportion of cgrow will go to growth. 0 = no growth.
	// fstorage = min(0.5, fstorage);
	double cgrow_pot = 1.0;
	if (iftreephen) {
		cgrow_pot = min(min(ftemp, fwater), fphoto);
	} 
	
	// Scale the carbon that is potentially headed for tree growth by the minimum of 
	// temperature and a water factor, assuming that either one can be the bottelneck 
	// for growth.
	cgrow *= cgrow_pot;	
	// if (store_state < xthreshold){
	// 	cgrow = 0;
	// } 
	// save the coeficients (for debugging)
	indiv.ftemp = ftemp;
	indiv.fwater = fwater;
	indiv.fstorage = fstorage;
	indiv.cgrow_pot = cgrow_pot;
	indiv.fphoto = fphoto;
	// cgrow *= min(ftemp, fwater);
	// actual allowance of c for growth today (kgC/m^2)
	return cgrow;
	
}

void reduce_tree_biomass(Individual &indiv, double frac) {
	// TODO DH: reconcile with Individual::reduce_biomass(). We don't want to
	// remove cmass_storage here - that's the main difference. We just want to
	// take leaves, roots, and sapwood.
	Patchpft &ppft = indiv.patchpft();

	double cmass_leaf_litter = frac * indiv.get_cmass_leaf();
	double cmass_root_litter = frac * indiv.cmass_root;
	double cmass_sap_litter = frac * indiv.cmass_sap;
	double cmass_repr_litter = frac * indiv.cmass_repr;

	double nmass_leaf_litter = frac * indiv.get_nmass_leaf();
	double nmass_root_litter = frac * indiv.nmass_root;
	double nmass_sap_litter = frac * indiv.nmass_sap;
	double nmass_repr_litter = frac * indiv.nmass_repr;

	Soil &soil = indiv.vegetation.patch.soil;
	Pft &pft = indiv.pft;

	leaf_to_litter(soil, pft, cmass_leaf_litter, nmass_leaf_litter);
	root_to_litter(soil, pft, cmass_root_litter, nmass_root_litter);

	// Sapwood turnover to heartwood.
	indiv.cmass_heart += cmass_sap_litter;
	indiv.nmass_heart += nmass_sap_litter;

	// See comment in turnover_dave_tree(). How should repr be handled? For now,
	// treat it the same as root litter.
	root_to_litter(soil, pft, cmass_repr_litter, nmass_repr_litter);

	// Used by DAVE mass balance checks; not an actual litter pool.
	indiv.ctolitter += cmass_leaf_litter + cmass_root_litter + cmass_repr_litter;
	indiv.ntolitter += nmass_leaf_litter + nmass_root_litter + nmass_repr_litter;

	double frac_remain = 1.0 - frac;

	indiv.cmass_leaf *= frac_remain;
	indiv.cmass_root *= frac_remain;
	indiv.cmass_sap *= frac_remain;
	indiv.cmass_repr *= frac_remain;

	indiv.nmass_leaf *= frac_remain;
	indiv.nmass_root *= frac_remain;
	indiv.nmass_sap *= frac_remain;
	indiv.nmass_repr *= frac_remain;

	indiv.densindiv *= frac_remain;
}

// Reduce density of indoviduals *daily* when carbon storage reduces below a threshold.
// To consider doing:
// 1. Do we want mortality to be storchastic or deterministic (depending on ifstochmort)?
//    currently it ignores ifstochmort completely.
bool carbon_starvation_mortality(Individual &indiv) {

	const Pft &pft = indiv.pft;
	Patchpft &ppft = indiv.patchpft();

	// Maximum fraction of cohort that can die per day.
	double max_mort_frac = pft.cstarv_max / Date::MAX_YEAR_LENGTH;

	// Negative NPP has already been subtracted from storage at this point.
	double storage_state = 0.0;
	if (indiv.cmass_storage > 0) {
		if (indiv.cstorage_max > 0) {
			storage_state = min(indiv.cmass_storage / indiv.cstorage_max, 1.0);
		}
	} else {
		log_info("%s killed by carbon starvation: negative NPP cannot be fulfilled by storage. NPP: %.9g kgC/m2, cmass_storage: %.9g kgC/m2",
				(const char *)indiv.pft.name, indiv.dnpp, indiv.cmass_storage);
		// next check: if adding the NPP to the storage caused storage to be 
		// negitive - kill tree and return.
		// Not enough C in storage to fulfill respiration demands --> kill individual
		double cmass_deficit = -indiv.cmass_storage;
		indiv.cmass_storage = 0.0;
		
		// ASSAF: this code was coppied from storage_growth_partitioning_tree() function... 
		// The respiration and NPP fluxes reported in canexch are
		// now incorrect, so we need to fix them. Note that the RA
		// flux is not used for mass balance checks (it's only used
		// for outputs). The NPP flux is used in mass balance checks
		// however, so it's important to correct it here.
		indiv.report_flux(Fluxes::RA, -cmass_deficit);
			// (npp has opposite sign to respiration)
		indiv.report_flux(Fluxes::NPP, cmass_deficit);
		// report carbon starvation mortality 
		indiv.mort_cstarve += 1.0; // annual
		indiv.mort_dcstarve = 1.0; // 100% mortality of this cohort (daily)
		ppft.amort_cstarve += 1.0;
		ppft.amort_cstarve_ntree += indiv.densindiv * patcharea;
		ppft.dmort_cstarve_ntree += indiv.densindiv * patcharea;
		indiv.kill();  // flushes all pools to litter
		indiv.vegetation.killobj();  // removes individual object
							
		//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
		// Don't think we need to call nextobj() here, as we've just
		// removed an element from the list; when this happens, the
		// "current" pointer internal to the list object points to
		// the element that was previously the "next" element.
		return false;
	}
		
	// calculate probability of mortality due to carbon starvation (unitless)
	double pmort_cstarve = max_mort_frac / (1.0 + pow(storage_state / pft.cstarv_inflect, 5.0));
	
    if (ifdailycstarvmort) {
		// Start thinning.
        if (indiv.densindiv > pft.cstarv_killdens) {
			// Accounting. Need to do this - particularly ntree-related
			// accounting - before reducing biomass.
			indiv.mort_dcstarve = pmort_cstarve;
			indiv.mort_cstarve += pmort_cstarve;
			ppft.amort_cstarve += pmort_cstarve;
			ppft.amort_cstarve_ntree += pmort_cstarve * indiv.densindiv * patcharea;
			ppft.dmort_cstarve_ntree += pmort_cstarve * indiv.densindiv * patcharea;

			reduce_tree_biomass(indiv, pmort_cstarve);
        } else {
			// kill individual
			log_info("%s killed by carbon starvation: Cohort density (%g/m2) has fallen below kill threshold of %g/m2",
				    (const char *)indiv.pft.name, indiv.densindiv, pft.cstarv_killdens);
			ppft.dmort_cstarve_ntree += indiv.densindiv * patcharea;
			ppft.amort_cstarve_ntree += indiv.densindiv * patcharea;
			ppft.amort_cstarve += 1.0;
			// No need to update indiv-level properties because it has died.
            indiv.kill();
            indiv.vegetation.killobj();
			return false;
        }
    }
	// We don't need to update allometry for the individual becasue 
	// it hasn't occured yet today
	
	return true;
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
	
	// save cmass_storage from time t-1
	double cmass_storage_before = indiv.cmass_storage;
	// first add NPP to storage.
	indiv.cmass_storage += indiv.dnpp;

	// implement *daily* carbon starvation mortality/thinning here (ASSAF).
	// note: we should keep the killing of the whole cohort when cmass_storage <0 
	// (below), but hopefully the thinning will alleviate the C starvation stress.
	if (!carbon_starvation_mortality(indiv)) {
		// in case the individual was killed there is no use in the next steps.
		return;
	}

	// In case the individual wasn't killed (only lost a couple of trees), continue...
	
	// calculate allowance (cgrow) depending on the enviro conditions and a  
	// maximum value (eg, 10% of storage).
	// Note: cgrow can be > 0 even if storage has gon wodn (ie dnpp < 0).
	indiv.cgrow = calculate_cgrow_tree(indiv);

	assert(indiv.cgrow < indiv.cmass_storage);
	assert(positive(indiv.cgrow));
	
	// remove cgrow from storage
	indiv.cmass_storage -= indiv.cgrow;

	// save cmass_storage from time t
	double cmass_storage_after = indiv.cmass_storage;

	// calculate the cmass storage dynamcs (> 0 if saving C)
	indiv.cmass_storage_dynam = cmass_storage_after - cmass_storage_before;

}

// Determine how much of the available C (NPP + C fraction from storage) goes to
// growth vs. storage. Equivalent to calc_fc_growth() in the grote code
// This is never used.
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

void storage_growth_partitioning_grote(Individual &indiv) {
	// This function calculates the potential growth.
	// Outputs: indiv.cgrow
	//
	// Anything we plan to spend must come out of storage here.
	// However, NPP has not yet gone into storage either.
	if (!ifgrote || !ifdailytree || !ifdavegrowth || indiv.pft.lifeform != TREE) {
		return;
	}

	double bminc = indiv.dnpp;
	if (bminc < 0) {
		// Negative NPP.
		indiv.cmass_storage += bminc;
		bminc = 0;
	}

	// implement *daily* carbon starvation mortality/thinning here (ASSAF).
	// note: we should keep the killing of the whole cohort when cmass_storage <0 
	// (below), but hopefully the thinning will alleviate the C starvation stress.
	if (!carbon_starvation_mortality(indiv)) {
		// in case the individual was killed there is no use in the next steps.
		return;
	}
	// In case the individual wasn't killed (only lost a couple of trees), continue...

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
			if (ifgrote) {
				storage_growth_partitioning_grote(indiv); 
				// TRY TO RUN THE NEW CGROW CALCULATION WITH NO PHENOLOGY (ONLY FSTORAGE)
				// storage_growth_partitioning_tree(indiv);
			} else {
				storage_growth_partitioning_tree(indiv);
			}
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
		double estab_scal_crown = indiv.pft.lifeform == GRASS_ANNUAL ? 0.0 : exp(-indiv.age_active / 1.0);
		if (indiv.pft.lifeform == GRASS_ANNUAL && indiv.growth_stage == reproductive) {
			calloc_crown = 0.0;
		} else {
			calloc_crown = min(indiv.pft.crownfrac * (1.0 + estab_scal_crown),1.0);
		}
		
		// 2) Reproduction
		if (indiv.growth_stage == reproductive) {
			double extra_crown_alloc = indiv.pft.crownfrac * estab_scal_crown;
			calloc_repr = min(max(0.0,indiv.pft.reprfrac - extra_crown_alloc),(1.0 - calloc_crown));
		} else {
			calloc_repr = 0.0;
		}

		// 3) Shoots (Leaves + culms)
		// express ltor_max as maximum fractional allocation to shoots (leaves + culms)
		double calloc_leaf_max = 1.0 - 1.0 / (1.0 + indiv.pft.ltor_max);

		// calculate fraction of C allocated to shoots (leaves + culms) following Johnson 2016 Eq. 3.80
		calloc_leaf = calloc_leaf_max * pow(indiv.nscal_5.mean() * indiv.get_wscal_5().mean(), 0.5);
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


// Calculate the leaf to root ratio of an individual for carbon allocation.
// deviates from the original formulation where there was a linear increase in 
// ltor (from 0.5-1) from min(wscal, nscal) of 0.5 times the ltor_max.
// original:
// (indiv.ltor = max(min(ppft.wscal_5.mean(), indiv.nscal_5.mean()) * pft.ltor_max, 0.5 * pft.ltor_max);)
// Assaf: The purpose of trying that is to reduce the senitivity of ltor at high 
// min(wscal, nscal) values (which significantly affects carbon allocation).
// Comment: conside including expon, mid_scal and ltor_min as PFT parameter values.
double calculate_ltor_tree(Individual &indiv){
	
	Pft &pft = indiv.pft;
	// get the patchpft object for this individual's pft
	Patchpft &ppft = indiv.patchpft();	
	// wscal_5 or nscal_5 which would represent the inflection point of the sigmoid fuction
	const double mid_scal = 0.5;
	// exponent of the ltor function
	const double expon = -10;
	const double ltor_min = 0.5 * pft.ltor_max;
	// state of soil stress
	double stress_state = min(indiv.get_wscal_5().mean(), indiv.nscal_5.mean());
	// 
	double ltor = ltor_min + ( (pft.ltor_max - ltor_min) / (1 + pow( (stress_state / mid_scal), expon) ) );
	
	// notrmalise to have a maximum of ltor_max (and not above/below)
	ltor *= ltor_min / ( (pft.ltor_max - ltor_min) / (1 + pow( (1 / mid_scal), expon) ) ) * pft.ltor_max;

	return ltor;
	
}

// Calculates the mean leaf age relative to the daily leaf increment (Holzworth et al. in press). 
double calculate_leaf_age(Individual &indiv) {
	// yesterday's leaf biomass (kgC/m2)
	double cmass_leaf_prev = indiv.cmass_leaf - indiv.cgrow_leaf;
	// day increment
	double deltat = 1.0 / Date::MAX_YEAR_LENGTH;

	double rel_daily_leaf_growth = cmass_leaf_prev / indiv.cmass_leaf;
	// return new leaf age
	return (rel_daily_leaf_growth * (indiv.meanleafage + deltat) ) ;
}

void calculate_pool_age(Individual &indiv) {
	
	// yesterday's pool biomass (kgC/m2)
	double cmass_leaf_prev = indiv.cmass_leaf - indiv.cgrow_leaf;
	double cmass_sap_prev = indiv.cmass_sap - indiv.cgrow_sap;
	double cmass_root_prev = indiv.cmass_root - indiv.cgrow_root;
	
	// day increment
	double deltat = 1.0 / Date::MAX_YEAR_LENGTH;

	double rel_daily_leaf_growth = cmass_leaf_prev / indiv.cmass_leaf;
	double rel_daily_sap_growth = cmass_sap_prev / indiv.cmass_sap;
	double rel_daily_root_growth = cmass_root_prev / indiv.cmass_root;

	// return new leaf age
	indiv.meanleafage =  rel_daily_leaf_growth * (indiv.meanleafage + deltat);
	indiv.meansapage =  rel_daily_sap_growth * (indiv.meansapage + deltat);
	indiv.meanrootage =  rel_daily_root_growth * (indiv.meanrootage + deltat);
}

void rescale_alphas(Individual &indiv) {

	//  sum of all the relative alphas so it will be used to calculate the alphas
	double sum_relative_alphas = indiv.absolute_alpha_leaf + indiv.absolute_alpha_root + 
						indiv.absolute_alpha_sap;
	if (ifgrote) {
		sum_relative_alphas += indiv.absolute_alpha_storage;
	}

	// calculate the alpha for each pool relative to the alpha of the rest of the pools.
	if (is_zero(sum_relative_alphas)) {
		indiv.alpha_leaf = 0;
		indiv.alpha_root = 0;
		indiv.alpha_sap  = 0;
		indiv.alpha_storage = 0;
	} else {
		// No need for protected divde here. If we refactor out the above check
		// for a zero value in the denominator, we should reinstate the use of a
		// protected division here.
		indiv.alpha_leaf = indiv.absolute_alpha_leaf / sum_relative_alphas;
		indiv.alpha_root = indiv.absolute_alpha_root / sum_relative_alphas;
		indiv.alpha_sap = indiv.absolute_alpha_sap / sum_relative_alphas;
		if (ifgrote) {
			indiv.alpha_storage = indiv.absolute_alpha_storage / sum_relative_alphas;
		}
		assert(equal(1.0, indiv.alpha_leaf + indiv.alpha_root + indiv.alpha_sap + (ifgrote ? indiv.alpha_storage : 0)));
	}
}

void calc_cgrow_pools_tree(Individual &indiv){
	// outputs:
	indiv.cgrow_leaf = indiv.cgrow * indiv.alpha_leaf;
	indiv.cgrow_root = indiv.cgrow * indiv.alpha_root;
	indiv.cgrow_sap = indiv.cgrow * indiv.alpha_sap;
	if (ifgrote) {
		indiv.cgrow_storage = indiv.cgrow * indiv.alpha_storage;
	}
}

// This variant of carbon_allocation_tree is using the Grote approach (of calulating
// sink strength using the difference between the limit and theactual caron for each pool),
// but driven by the leaf and not sapwood turnover. Allocation ot leaves is driven by 
// the ratio between par and par_max, which peaks at 1 in summer. Then, the sapwood limit 
// is calculated to maintain k_latosa.
void carbon_allocation_tree_assaf(Individual &indiv) {
	Pft &pft = indiv.pft;
	// get the patchpft object for this individual's pft
	Patchpft &ppft = indiv.patchpft();	
	// get today's climate 
	const Climate &climate = indiv.vegetation.patch.get_climate();

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
	const double repr_prop = 0.01;
	// minimum age for reproduction in years.
	// TODO: fix to be a function of tree size and not age, and a PFT parameter.
	// also, check that age is counted in years and not days.
	const double min_age_for_repr = 10;
	// maximum proportion of growth that will be allocated to reproduction
	// TODO: check this out. guestimated.
	const double max_alpha_repr_factor = 0.2;

	// calculate the potential potential proportion of cgrow that will go to leaves (cgrow_leaf)
	// maximum fraction of cgrow that can go to leaves.
	double alpha_leaf_pot_max = .3;
	// exponent 
	double alpha_leaf_pot_expon = -4.0; // was 0.4
	// Photosynthetic active readiation as a proportion of the maximum at this gridcell.
	// Represents the *light conditions*.
	double rel_par = max(0.05, climate.rel_par);
	
	// potential potential proportion of cgrow that will go to leaves
	double alpha_leaf_pot = alpha_leaf_pot_max * (1 - exp(alpha_leaf_pot_expon * rel_par));
	
	// // !!!! TRYING A NEW APPROACH. LEAF ALLLOCATION IS RELATED TO LEAF WATER POTENTIAL (HIGHER - MORE ALLOCATION)
	// alpha_leaf_pot = alpha_leaf_pot_max * (1 - exp(-2.0 * indiv.get_wscal_30().mean()));
	
	// calculate the growth demand of leaves, depending on the *light conditions*
	indiv.cmass_leaf_limit = indiv.cmass_leaf + alpha_leaf_pot * indiv.cgrow;

	double klatosa = pft.k_latosa;
	// double klatosa = indiv.latosa;
	double sla = indiv.sla();
	double rho = pft.wooddens;

	// calculate growth demand for the sapwood by maintaining k_latosa (responding to leaves)
	indiv.cmass_sap_limit = max( (sla * indiv.cmass_leaf * rho * height) / klatosa , indiv.cmass_sap);
	
	// an alternative ltor equation (for testing). Reason: ltor too sensitive to wscal,
	// which can increase the proportion of C going to root unnecesseraly.
	double ltor = calculate_ltor_tree(indiv);
	indiv.ltor = max(min(1.0 , ltor), 0.5 * pft.ltor_max);
	// indiv.ltor = 1.0; // TEST

	assert(indiv.ltor <= 1 && indiv.ltor >= 0);

	// cmass root target according to the leaf optimum and ltor (kgC/m^2)
	indiv.cmass_root_limit = 0.0;
	if (indiv.ltor > 0.0) {
		indiv.cmass_root_limit = max(indiv.cmass_leaf / indiv.ltor, indiv.cmass_root);
		// indiv.cmass_root_limit = indiv.cmass_root + indiv.cmass_root * (pft.turnover_root / Date::MAX_YEAR_LENGTH / fwater);
	}
	
	// cmass reproduction target
	// TODO: cmass_repr_limit has been deleted. Will need to refactor this if we
	// ever want to revive it.
	// indiv.cmass_repr_limit = 0.0;
	// if (indiv.age > min_age_for_repr) {
	// 	indiv.cmass_repr_limit = repr_prop * indiv.cgrow + indiv.cmass_repr;
	// }

	// Calculate the alpha factor (what proportion to allocate for each pool) according
	// to the distance of actual C to it's target.
	if (!is_zero(indiv.cmass_leaf_limit)) {
		indiv.absolute_alpha_leaf = divide(indiv.cmass_leaf_limit - indiv.cmass_leaf, indiv.cmass_leaf_limit);
	}

	if (!is_zero(indiv.cmass_root_limit)) {
		indiv.absolute_alpha_root = divide(indiv.cmass_root_limit - indiv.cmass_root, indiv.cmass_root_limit);
	}
	
	double relative_alpha_sap_temp = 0;
	if (!is_zero(indiv.cmass_sap_limit)) {
		indiv.absolute_alpha_sap = divide(indiv.cmass_sap_limit - indiv.cmass_sap, indiv.cmass_sap_limit);	
	}

	// To prevent allocation to reproduction if the tree is stressed (lower priority.),
	// 1. calculate the mean of the alphas or leaf, root and sapwood, that will 
	// be used ased as a threshold for allocating to reproduction. 
	// double mean_relative_alphas = (relative_alpha_leaf + relative_alpha_root + relative_alpha_sap + relative_alpha_storage) / 4.0;
	double mean_relative_alphas = (indiv.absolute_alpha_leaf + indiv.absolute_alpha_root + indiv.absolute_alpha_sap) / 3.0;
	// logistic function of the relative strength of  allocation to reproduction
	// depending on how far the other pools are from their target.
	// TODO: Unsupported assumption: logistic function assumed and it's parameters 
	// guestimated.  Change that.

	double alpha_repr_factor = (max_alpha_repr_factor / 
                       (1 + exp(20 * (mean_relative_alphas - .25))));
	// calculate alpha reproduction "tunned down" by the factor that depends on
	// the conditions of the other pools.
	// TODO: cmass_repr_limit has been deleted. Will need to refactor this if we
	// ever want to revive it.
	// if (indiv.cmass_repr_limit > indiv.cmass_repr) {
	// 	// indiv.relative_alpha_repr = alpha_repr_factor * 
	// 	// 	(indiv.cmass_repr_limit - indiv.cmass_repr) / indiv.cmass_repr_limit;
	// 		indiv.relative_alpha_repr = 0.1 * (indiv.relative_alpha_leaf + 
	// 			indiv.relative_alpha_root + indiv.relative_alpha_sap);
	// }
	// caps the relative_alpha_repr to 0.1 to avoid over reprentation of cgrow_repr
	// over the rest of the pools
	// TODO: Fix this (0.1???)
	// TODO: relative_alpha_repr has been deleted. Will need to refactor this if we
	// ever want to revive it.
	// indiv.relative_alpha_repr = min(indiv.relative_alpha_repr, 0.1);

	// rescale the alphas according to the demand (output: indiv.alpha_*)
	rescale_alphas(indiv);
	// // ASSAF: try fixed fractions
	// indiv.alpha_leaf = 0.3;
	// indiv.alpha_root = 0.3;
	// indiv.alpha_sap = 0.4;
	// indiv.alpha_repr = 0.0;

	// calculate cgrow for each of the pools (output: indiv.cgrow_*)
	calc_cgrow_pools_tree(indiv);

	
	// // rescale alpha and then cgrow according to k_latosa
	// double cmass_leaf_pot = indiv.cmass_leaf + indiv.cgrow_leaf;
	// double cmass_sap_pot = indiv.cmass_sap + indiv.cgrow_sap;
	// double k2 = pft.k_allom2;
	// double k3 = pft.k_allom3;
	
	
	// const double k_latosa_max = 8000;
	// double height_pot = pow(4.0 * (cmass_sap_pot / indiv.densindiv) * pow(k2, 2.0 / k3) / (PI * rho) , 1.0 / (1.0 + 2.0 / k3));
	// height_pot = indiv.height;
	// double k_latosa_pot = (cmass_leaf_pot * pft.sla * height_pot * rho) / cmass_sap_pot;
	// const double F = 10;
	// double smax = (k_latosa_max - pft.k_latosa) / pft.k_latosa;
	// const double f_scalar = F / smax;
	// if (indiv.latosa > pft.k_latosa) {
	// 	double scalar = 1.0 + f_scalar * (indiv.latosa - pft.k_latosa) / pft.k_latosa;
	// 	indiv.relative_alpha_sap *= scalar;
	// 	rescale_alphas(indiv);
	// 	calc_cgrow_pools_tree(indiv);

	// } else if (indiv.latosa < pft.k_latosa) {
	// 	double scalar = 1.0 + f_scalar * (pft.k_latosa - indiv.latosa) / pft.k_latosa;
	// 	indiv.relative_alpha_leaf *= scalar;
	// 	rescale_alphas(indiv);
	// 	calc_cgrow_pools_tree(indiv);
	// }

	assert(positive(indiv.alpha_leaf));
	assert(positive(indiv.alpha_root));
	assert(positive(indiv.alpha_sap));
	// assert(positive(indiv.alpha_repr));
	assert(positive(indiv.cgrow_leaf));
	assert(positive(indiv.cgrow_root));
	assert(positive(indiv.cgrow_sap));
	// assert(positive(indiv.cgrow_repr));
}

double calc_alpha(double limit, double actual, double turnover_rate) {
	// TODO: could make this a function of the turnover rate. If we do this, we
	// need to remember that turnover_rate could be 0 (e.g. for storage).
	if (is_zero(limit)) {
		return 0;
	}
	double alpha = divide(limit - actual, limit); // relative alpha
	// double alpha = limit - actual; // absolute alpha

	// if (!is_zero(turnover_rate)) {
	// 	alpha *= turnover_rate;
	// }

	// Constrain alpha to [0, 1].
	alpha = min(1.0, alpha);
	alpha = max(0.0, alpha);
	return alpha;
}

// Determine the sink strength of each carbon pool in the individual.
// This probably encapsulates adjust_cpool_optimum() and decide_on_allocation()
// in the grote code.
void carbon_allocation_tree_grote(Individual &indiv) {
	// Inputs: cmass_<pools>, growing conditions maybe?
	// Outputs: indiv.cgrow_leaf, indiv.cgrow_root, indiv.cgrow_sap, indiv.cgrow_repr

	// Steps: 
	// 1) calculate the target growth (a-la Grote) for each pool.
	// 2) calculate the alpha factor (decide_on_allocation) using sink strength 
	//    (ie the difference between the actual pool size and optimal growth)

	const Pft &pft = indiv.pft;
	// get the patchpft object for this individual's pft
	const Patchpft &ppft = indiv.patchpft();	
	// get today's climate 
	const Climate &climate = indiv.vegetation.patch.get_climate();

	// 1) calculate the target growth (a-la Grote) for each pool.
	// calculate leaf limit (optimum)
	const double MAX_HEIGHT = 150;
	const double MIN_HEIGHT = 0.05;

	// minimum age for reproduction in years.
	// TODO: fix to be a function of tree size and not age, and a PFT parameter.
	// also, check that age is counted in years and not days.
	const double MIN_AGE_FOR_REPR = 10;

	// minimum storage fraction for reproduction
	const double MIN_STORFRAC_FOR_REPR = 0.5;

	// Number of days to look ahead for the targets.
	const int NDAY_OUTLOOK = 365;

	// limit height to below max and above min (to avoid FNS > 1)(from David's 
	// code). 
	// TODO: what about seedlings????
	double height = min(max(MIN_HEIGHT, indiv.height), MAX_HEIGHT);
	
	// Tying the form factor with the wood density (just for convinience). This
	// is appropriate here becuase in the FNS calculation, sapwood biomass is
	// calculated using the assumption of a cylindrical stem.
	double wooddens_FF = pft.wooddens * indiv.pft.formfactor;

	// Calculate the ltor according to the leaf optimum and below-ground stress.
	// TODO: think what to do when nscal or wscal is zero. Now set a minimum
	// ltor to .5 of ltor_max.
	//
	// An alternative ltor equation to test. Reason: ltor too sensitive to
	// wscal, which can increase the proportion of C going to root
	// unnecesseraly. indiv.ltor = calculate_ltor_tree(indiv);
	indiv.ltor = max(min(indiv.get_wscal_30().mean(), indiv.nscal_30.mean()) * pft.ltor_max, 0.5 * pft.ltor_max);
	assert(indiv.ltor <= pft.ltor_max && indiv.ltor >= 0);

	// Update storage capacity.
	// TODO: is this necessary? Pretty sure this is handled elsewhere.
	indiv.cstorage_max = calc_creserve_opt_grote(indiv);
	assert(indiv.cstorage_max > 0);

	// Allocate to reproduction when:
	// - age > MIN_AGE_FOR_REPR
	// - storage_fraction > MIN_STORFRAC_FOR_REPR
	double storage_fraction = indiv.cmass_storage / indiv.cstorage_max;
	if (indiv.age > MIN_AGE_FOR_REPR && storage_fraction > MIN_STORFRAC_FOR_REPR) {
		indiv.cgrow_repr = indiv.pft.reprfrac * indiv.cgrow;
	}

	// Remaining carbon available for growth (kgC/m2).
	double cgrow = indiv.cgrow - indiv.cgrow_repr;

	// Sapwood target (to grow at least the turnover of N days) (kgC/m^2). This
	// drives the growth of sapwood AND the leaves (in the next day). TODO: to
	// get the observed bi-modal sapwood growth and unimodal leaf growth, we
	// need to change that.
	double growth_rate = 1.0 + pft.turnover_sap / Date::MAX_YEAR_LENGTH;
	indiv.cmass_sap_limit = indiv.cmass_sap * pow(growth_rate, 1.0 * NDAY_OUTLOOK);

	// Leaf target (to grow to supply leaves to the current sapwood ) (kgC/m^2).
	indiv.cmass_leaf_limit = max(indiv.cmass_sap_limit * (pft.k_latosa / (indiv.sla() * wooddens_FF * height)) , indiv.cmass_leaf);

	// Root target according to the leaf optimum and ltor (kgC/m^2).
	if (indiv.ltor > 0.0) {
		indiv.cmass_root_limit = max(indiv.cmass_leaf_limit / indiv.ltor, indiv.cmass_root);
	} else {
		indiv.cmass_root_limit = indiv.cmass_root;
	}

	// Storage target (kgC/m^2). This is 80% of capacity, so we don't try to
	// constantly fill up, which results in frequent overflows (e.g. when
	// capacity decreases slightly in a low-growth period).
	indiv.cmass_storage_limit = 0.8 * indiv.cstorage_max;

	// Calculate the alpha factor (what proportion to allocate for each pool)
	// according to the distance of actual C to it's target. Alphas are zeroed
	// at the start of each day.
	indiv.absolute_alpha_leaf = calc_alpha(indiv.cmass_leaf_limit, indiv.cmass_leaf, pft.turnover_leaf);
	indiv.absolute_alpha_root = calc_alpha(indiv.cmass_root_limit, indiv.cmass_root, pft.turnover_root);
	indiv.absolute_alpha_sap = calc_alpha(indiv.cmass_sap_limit, indiv.cmass_sap, pft.turnover_sap);
	indiv.absolute_alpha_storage = calc_alpha(indiv.cmass_storage_limit, indiv.cmass_storage, 1 - indiv.cmass_storage / indiv.cstorage_max);

	// Rescale the alphas according to the demand.
	rescale_alphas(indiv);

	indiv.alpha_leaf_trend.add(indiv.alpha_leaf);
	indiv.alpha_root_trend.add(indiv.alpha_root);
	indiv.alpha_sap_trend.add(indiv.alpha_sap);
	indiv.alpha_storage_trend.add(indiv.alpha_storage);

	indiv.alpha_leaf = indiv.alpha_leaf_trend.mean();
	indiv.alpha_root = indiv.alpha_root_trend.mean();
	indiv.alpha_sap = indiv.alpha_sap_trend.mean();
	indiv.alpha_storage = indiv.alpha_storage_trend.mean();

	// Calculate cgrow for each of the pools. We need to use cgrow (not
	// indiv.cgrow), because the alphas are fractions of the remaining cgrow,
	// after allocation to reproduction.
	indiv.cgrow_leaf = cgrow * indiv.alpha_leaf;
	indiv.cgrow_root = cgrow * indiv.alpha_root;
	indiv.cgrow_sap = cgrow * indiv.alpha_sap;
	indiv.cgrow_storage = cgrow * indiv.alpha_storage;

	assert(positive(indiv.alpha_leaf));
	assert(positive(indiv.alpha_root));
	assert(positive(indiv.alpha_sap));
	assert(positive(indiv.alpha_storage));

	assert(indiv.cgrow_leaf <= indiv.cgrow);
	assert(indiv.cgrow_root <= indiv.cgrow);
	assert(indiv.cgrow_sap <= indiv.cgrow);
	assert(indiv.cgrow_storage <= indiv.cgrow);
	assert(indiv.cgrow_repr <= indiv.cgrow);

	assert(equal(indiv.cgrow, indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_sap + indiv.cgrow_storage + indiv.cgrow_repr));
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
			// Choose between 3 approaches (temporary)
			if (ifgrote) {
				// 1. Grote - driven by sapwood turnover (grote-Peter Anthony)
				carbon_allocation_tree_grote(indiv);
			} else {
				// 2. Driven by leaf turnover
				carbon_allocation_tree_assaf(indiv);
			}
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
	
	const Pft &pft = indiv.pft;

	// Parameters
	const double frac_nchange = 0.05;

	// Local	
	double cton_demand_leaf = 0.0;
	double cton_demand_root = 0.0;
	double cton_demand_crown = 0.0;
	double cton_demand_repr = 0.0;

	double ngrow_demand_leaf  = 0.0;
	double ngrow_demand_root  = 0.0;
	double ngrow_demand_repr  = 0.0;
	double ngrow_demand_crown = 0.0;

	// update growth depending on N stress of plant
	// TODO JK: double check cases where nmass_storage > nstorage_max!!
	// TODO JK: constrain cton ratio ;) 
	double nstorage_fill = min(0.999999, max(0.000001, (indiv.nmass_storage / indiv.nstorage_max)));

	//double fN_storage = (1.0 - sqrt(nstorage_fill));
	double fN_storage = (1.0 - nstorage_fill);
	cton_demand_leaf = fN_storage * (indiv.pft.cton_leaf_max - indiv.pft.cton_leaf_min) + indiv.pft.cton_leaf_min;
	cton_demand_root = fN_storage * (indiv.pft.cton_root_max - indiv.pft.cton_root_min) + indiv.pft.cton_root_min;
	cton_demand_crown = fN_storage * (indiv.pft.cton_root_max - indiv.pft.cton_root_min) + indiv.pft.cton_root_min;
	cton_demand_repr = fN_storage * (indiv.pft.cton_repr_max - indiv.pft.cton_repr_min) + indiv.pft.cton_repr_min;

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
		
		// N demand for new growth
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

	// N demand of existing biomass
	// We allow this value to be negative to allow for active retranslocation of N within the plant
	// This often applies to cases where new growth is limited by N but existing tissues have plenty.
	double nupdate_demand_leaf  = (indiv.get_cmass_leaf() / cton_demand_leaf - indiv.get_nmass_leaf()) * frac_nchange;
	double nupdate_demand_root  = (indiv.cmass_root / cton_demand_root - indiv.nmass_root) * frac_nchange;
	double nupdate_demand_crown = (indiv.cmass_crown / cton_demand_crown - indiv.nmass_crown) * frac_nchange;
	double nupdate_demand_repr  = (indiv.cmass_repr / cton_demand_repr - indiv.nmass_repr) * frac_nchange * 0.1;

	// If demand is negative, give N back to storage and set demand to 0.
	if (nupdate_demand_leaf < 0.0) {
		indiv.nmass_storage += -nupdate_demand_leaf;
		if (pft.is_grass() || pft.lifeform == MOSS) {
			indiv.nmass_leafs[green] += nupdate_demand_leaf;
		} else {
			indiv.nmass_leaf += nupdate_demand_leaf;
		}
		nupdate_demand_leaf = 0.0;
	}
	if (nupdate_demand_root < 0.0) {
		indiv.nmass_storage += -nupdate_demand_root;
		indiv.nmass_root += nupdate_demand_root;
		nupdate_demand_root = 0.0;
	}
	if (nupdate_demand_crown < 0.0) {
		indiv.nmass_storage += -nupdate_demand_crown;
		indiv.nmass_crown += nupdate_demand_crown;
		nupdate_demand_crown = 0.0;
	}
	if (nupdate_demand_repr < 0.0) {
		indiv.nmass_storage += -nupdate_demand_repr;
		indiv.nmass_repr += nupdate_demand_repr;
		nupdate_demand_repr = 0.0;
	}
	// TODO: do we need to check for stoichiometry here or can we rely on the one in update_pools?

	// Calculate total N demand from growth and redistribution
	double ntotal_demand_leaf = ngrow_demand_leaf + nupdate_demand_leaf;
	double ntotal_demand_root = ngrow_demand_root + nupdate_demand_root;
	double ntotal_demand_crown = ngrow_demand_crown + nupdate_demand_crown;
	double ntotal_demand_repr = ngrow_demand_repr + nupdate_demand_repr;

	// Take N demand of storage into account
	indiv.storendemand = max(0.0, indiv.nstorage_max - indiv.nmass_storage);

	// Total biomass N demand of individual
	indiv.leafndemand  = ntotal_demand_leaf;
	indiv.rootndemand  = ntotal_demand_root;
	indiv.crownndemand = ntotal_demand_crown;
	indiv.reprndemand  = ntotal_demand_repr;
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
	// TODO JK: this implies that growth is limited if current C:N is high, which shouldn't be the case. Make growth a priority?
	double nupdate_demand_leaf  = max(0.0, indiv.get_cmass_leaf() / indiv.pft.cton_leaf_min - indiv.get_nmass_leaf()) * frac_nchange;
	double nupdate_demand_root  = max(0.0, indiv.cmass_root / indiv.pft.cton_root_min - indiv.nmass_root) * frac_nchange;
	double nupdate_demand_sap   = max(0.0, indiv.cmass_sap / indiv.pft.cton_sap_min - indiv.nmass_sap) * frac_nchange;
	double nupdate_demand_repr  = max(0.0, indiv.cmass_repr / indiv.pft.cton_repr_min - indiv.nmass_repr) * frac_nchange;

	// Make sure that nupdate_demands do not exceed N demands for growth (i.e. prioritise growth over existing biomass)
	// Alternatively, one could ask for the min C:N for new growth and not update existing biomass.
	double ngrow_demand_total = ngrow_demand_leaf + ngrow_demand_root + ngrow_demand_sap + ngrow_demand_repr;
	double nupdate_demand_total = max(0.0, nupdate_demand_leaf + nupdate_demand_root + nupdate_demand_sap + nupdate_demand_repr);
	double growth_update_ratio = min(max_update_growth, ngrow_demand_total / nupdate_demand_total);

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
		// TODO JK: revisit this calculation and check where cton_status is used.
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

        // Nitrogen availability based on root projective cover
        double max_avail = min(1.0, indiv.rpc) * nmin_avail;
 
        double fnuptake = competing_individuals[i].fnuptake;
        double ndemand = competing_individuals[i].ndemand;
        double nuptake = ndemand * fnuptake;
 
        // Check if N supply is limited by RPC. If so, modify fnuptake such that
        // the total uptake equals max_avail. This can only reduce total uptake.
        if (max_avail < nuptake) {
            competing_individuals[i].fnuptake = max_avail / ndemand;
        }

		// TODO: verify that this belongs here.
		indiv.fnuptake = competing_individuals[i].fnuptake;
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

/*
Calculate the fraction of today's N demand which can be met by supply.
*/
double fnavail_dave(const Individual &indiv) {
	//fnavail = max(0.0,min(1.0,(indiv.nmass_storage / indiv.ndemand)));
	// TODO: define max_mobil_nstorage based on Yang et al. 2023, Journal of Ecology.
	double supply = indiv.max_mobil_cstorage * indiv.nmass_storage;

	// If N uptake is going into its own pool, this is part of the supply.
	if (!indiv.pft.store_uptake()) {
		supply += indiv.nmass_uptake;
	}

	double demand = indiv.ndemand - indiv.storendemand;

	if (is_zero(demand)) {
		return 1;
	}

	return max(0.0, min(1.0, supply / demand));
}

/*
Calculate the N stress for grass. 0 = max stress, 1 = no stress.
*/
double nstress_dave_grass(const Individual &indiv) {
	// N stress is calculated based on how many days of N reserves we have, given the current N demand (reserve_days).
	double reserve_days = indiv.nmass_storage / (max(1e-15,indiv.ndemand - indiv.storendemand));
	double reserve_days_50 = 2.0; // number of reserve days that give stress of 0.5
	double nstress = 1.0 / (1.0 + exp(-(reserve_days - reserve_days_50)));
	return max(0.0, min(1.0, nstress));
}

/*
Calculate the N stress for trees. 0 = max stress, 1 = no stress.
*/
double nstress_dave_tree(const Individual &indiv) {
	// N stress for trees based on how far the leaf N is below the max N:C.
	// N stress is 1 when at max N:C, and 0 when at min N:C.
	double ntoc = indiv.get_nmass_leaf() / indiv.get_cmass_leaf();
	double ntoc_max = 1.0 / indiv.pft.cton_leaf_min;
	double ntoc_min = 1.0 / indiv.pft.cton_leaf_max;
	double nstress = (ntoc_max - ntoc) / (ntoc_max - ntoc_min);
	return max(0.0, min(1.0, nstress));
}

/*
Calculate the N stress and N availability for an individual.

@param fnavail: (output) The fraction of today's N demand which can be met by supply.
@param nstress: (output) The N stress for the individual.
*/
void nstress_dave(const Individual &indiv, double &fnavail, double &nstress) {
	// If no demand, reuse yesterday's N stress.
	if (negligible(indiv.ndemand)) {
		nstress = indiv.nscal_5.lastadd();
		fnavail = 0;
		return;
	}

	fnavail = fnavail_dave(indiv);

	// No N stress in free N years.
	if (date.year < freenyears) {
		nstress = 1;
	} else if (indiv.pft.is_grass()) {
		nstress = nstress_dave_grass(indiv);
	} else {
		// nstress = nstress_dave_tree(indiv);
		nstress = nstress_dave_grass(indiv);
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
			if (indiv.initialn > 0) {
				double nvoid = min(indiv.initialn, nexcess);
				indiv.initialn -= nvoid;
				indiv.ntolitter += nvoid;
				nexcess -= nvoid;
				vegetation.patch.fluxes.report_flux(Fluxes::ESTN, nvoid);
			}
			indiv.ntolitter += nexcess;
			ppft.nmass_litter_root += nexcess;
			indiv.nmass_storage = indiv.nstorage_max;
		}

		// TODO: Check if there is enough N in storage. At the moment Nstorage needs to be > 0.
		// Alternatively, one could implement a N:C storage > N:C leaf min constraint.

		// Calculate N stress factor.
		double fnavail = 0.0;
		double nstress = 0.0;
		nstress_dave(indiv, fnavail, nstress);

		assert(fnavail >= 0.0);
		assert(fnavail <= 1.0);
		assert(nstress >= 0.0);
		assert(nstress <= 1.0);

		// nstress will be 1 in freenyears.
		indiv.nscal_5.add(nstress);
		indiv.nscal_30.add(nstress);

		// Get available N for each pool.
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

		// Growth respiration as fraction of new growth.
		const double RESP_GROWTH_FRAC = 0.25;

		double resp_growth_leaf = 0;
		double resp_growth_root = 0;
		double resp_growth_crown = 0;
		double resp_growth_sap = 0;
		double resp_growth_repr = 0;

		// N supply < N demand
		double nreturn = 0;

		// Leaves
		constrain_growth(
			indiv.cgrow_leaf, navail_leaf,
			indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
			creturn, nreturn, RESP_GROWTH_FRAC, resp_growth_leaf);

		// Roots
		constrain_growth(
			indiv.cgrow_root, navail_root,
			indiv.pft.cton_root_min, indiv.pft.cton_root_max,
			creturn, nreturn, RESP_GROWTH_FRAC, resp_growth_root);

		// Crowns (for grasses and mosses only).
		if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			constrain_growth(
				indiv.cgrow_crown, navail_crown,
				indiv.pft.cton_root_min, indiv.pft.cton_root_max,
				creturn, nreturn, RESP_GROWTH_FRAC, resp_growth_crown);
		}

		// sapwood
		if (indiv.pft.lifeform == TREE) {
			constrain_growth(
				indiv.cgrow_sap, navail_sap,
				indiv.pft.cton_sap_min, indiv.pft.cton_sap_max,
				creturn, nreturn, RESP_GROWTH_FRAC, resp_growth_sap);
		}

		// Reproductive tissue
		constrain_growth(
			indiv.cgrow_repr, navail_repr,
			indiv.pft.cton_repr_min, indiv.pft.cton_repr_max,
			creturn, nreturn, RESP_GROWTH_FRAC, resp_growth_repr);

		assert(positive(creturn));
		indiv.resp_growth = resp_growth_leaf + resp_growth_root + resp_growth_crown + resp_growth_sap + resp_growth_repr;
		indiv.dnpp -= indiv.resp_growth;
		indiv.report_flux(Fluxes::NPP, -indiv.resp_growth);
		indiv.report_flux(Fluxes::RA, indiv.resp_growth);

		if (ifgrote) {
			// cgrow_storage goes back into storage.
			creturn += indiv.cgrow_storage;
		}

		// Make sure storage does not overflow.
		// Maximum C storage pool size
		if (indiv.cmass_storage + creturn > indiv.cstorage_max) {
			cexcess = indiv.cmass_storage + creturn - indiv.cstorage_max;
			indiv.cexcess += cexcess;
			creturn = indiv.cstorage_max - indiv.cmass_storage;
		}
		indiv.cmass_storage += creturn;
		indiv.creturn = creturn;

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

		indiv.acalloc_leaf += indiv.cgrow_leaf;
		indiv.acalloc_root += indiv.cgrow_root;
		indiv.acalloc_sap += indiv.cgrow_sap;
		indiv.acalloc_repr += indiv.cgrow_repr;
		indiv.acalloc_storage += indiv.cgrow_storage;
		indiv.acalloc_crown += indiv.cgrow_crown;

		// Verify that growth respiration is correctly related to growth.
		if (indiv.cgrow_leaf > 0) {
			assert(equal(RESP_GROWTH_FRAC, resp_growth_leaf / indiv.cgrow_leaf));
		}
		if (indiv.cgrow_root > 0) {
			assert(equal(RESP_GROWTH_FRAC, resp_growth_root / indiv.cgrow_root));
		}
		if (indiv.cgrow_sap > 0) {
			assert(equal(RESP_GROWTH_FRAC, resp_growth_sap / indiv.cgrow_sap));
		}
		if (indiv.cgrow_repr > 0) {
			assert(equal(RESP_GROWTH_FRAC, resp_growth_repr / indiv.cgrow_repr));
		}
		if (indiv.cgrow_crown > 0) {
			assert(equal(RESP_GROWTH_FRAC, resp_growth_crown / indiv.cgrow_crown));
		}

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

		double ngrow = navail_leaf + navail_root + navail_crown + navail_sap + navail_repr;
		indiv.nstorage_drawdown = 0;
		if (!pft.store_uptake()) {
			// First take ngrow out of uptake pool.
			indiv.nmass_uptake -= ngrow;
			if (indiv.nmass_uptake < 0) {
				// Supply remaining N from storage.
				double drawdown = -indiv.nmass_uptake; // positive sign now
				indiv.nmass_storage = max(0.0, indiv.nmass_storage - drawdown);
				indiv.nmass_uptake = 0;
				indiv.nstorage_drawdown = drawdown;
			} else if (indiv.nmass_uptake > 0) {
				// N supply was met by uptake. If any uptake remains, put it
				// back into the available soil N, as if it had never been
				// taken up.
				//
				// We can't subtract from dnuptake because it's set later in the
				// day in vegetation_n_uptake(). Therefore, we account for this
				// using a separate variable: ntoavail. The same is true for the
				// annual variables (anuptake, antoavail). This isn't really too
				// important, except to ensure that we've accounted for all
				// fluxes so that the nitrogen budget outputs make sense.
				indiv.ntoavail += indiv.nmass_uptake;
				indiv.antoavail += indiv.nmass_uptake;

				vegetation.patch.soil.nmass_inc(indiv.nmass_uptake);
				indiv.nmass_uptake = 0;
			}
		} else {
			// Draw down on storage.
			indiv.nstorage_drawdown = ngrow;
			indiv.nmass_storage = max(0.0, indiv.nmass_storage - ngrow);
		}

		if (date.day == 0) {
			indiv.anstorage_drawdown = 0;
		}
		indiv.anstorage_drawdown += indiv.nstorage_drawdown;

		// Save yearly leaf and root growth. 
		// TODO JK: check if crown and repr need to be added
		indiv.cmass_leaf_ygrowth += indiv.cgrow_leaf;
		indiv.cmass_root_ygrowth += indiv.cgrow_root;

		// Update stand-level cmass_repr used for annual establishment.
		// This variable is zeroed in Standpft::annual_accounting().
		patch.stand.pft[pft.id].cmass_repr += indiv.cgrow_repr;

		// update indiv.cgrow (actual growth)
		indiv.cgrow = indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_crown + indiv.cgrow_sap + indiv.cgrow_repr;
		if (ifgrote) {
			indiv.cgrow += indiv.cgrow_storage;
		}

		// update mean pool age (years)
		// indiv.meanleafage = calculate_leaf_age(indiv);
		calculate_pool_age(indiv);

		assert(positive(indiv.cgrow_leaf));
		assert(positive(indiv.cgrow_root));
		assert(positive(indiv.cgrow_crown));
		assert(positive(indiv.cgrow_sap));
		assert(positive(indiv.cgrow_repr));
		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));
		assert(less_than_or_equal(indiv.nmass_storage, indiv.nstorage_max));
		assert(equal(cgrow_unlim, indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_crown + indiv.cgrow_sap + indiv.cgrow_repr + creturn + cexcess + indiv.resp_growth));

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

		if (!negligible(indiv.nmass_repr)) {
			assert(less_than_or_equal(indiv.cmass_repr / indiv.nmass_repr, indiv.pft.cton_repr_max * 1.001));
		}
	}
}

void senescence_dave_grass(Individual &indiv) {

	const Pft &pft = indiv.pft;
	Vegetation &vegetation = indiv.vegetation;

	// TODO: is there an active retranslocation of storage C and N in senescing
	// annual plants?

	// Air temperature or bottom snow layer if snow present
	// Note that grass_temp is daily, not half-hourly. 
	// At the moment, we account for that (imperfectly) by introducing Tfreeze below.
	// TODO: replace with actual minimimum subdaily T_soil[6].
	double grass_temp = vegetation.patch.soil.T_soil[6];

	// Baseline (minimum) senescence rate from leaf longevity under unstressed conditions.
	// double senesc_base = 1.0 / (indiv.pft.leaflong * date.year_length());
	// Note: leaflong not used anymore for senescence as it is conceptually difficult

	double senesc = pft.senesc_base;

	// Developmental effects on senescence (stronger senescence in reproductive phase)
	if (indiv.growth_stage == reproductive || indiv.fdev_growth < 1.0) {
		// stronger increase of senescence in annuals
		//double senesc_rate_repr = indiv.pft.lifeform == GRASS_PERENNIAL ? 1.015 : 1.08;
		//indiv.fdev_senesc = max(1.0, indiv.fdev_senesc * senesc_rate_repr);
		//indiv.fdev_senesc  = min(25.0, pow(indiv.fdev_growth, -0.1));
		//indiv.fdev_senesc = senesc_scal_repr * (1.0 - indiv.fdev_growth);  // senescing too early
		// Note: multiplied by fphen_activity as fdev_senesc transitions to senesc_fdorm during dormancy
		indiv.fdev_senesc = pft.senesc_scal_repr * pow((1.0 - indiv.fdev_growth),
													   pft.qdev_senesc) * indiv.fphen_activity;
	} else {
		indiv.fdev_senesc = 0.0;
	}

	// Stress factors increase senescence
	// Soil moisture
	// TODO: multiply this by the fraction of soil that is not frozen!
	double senesc_fwater = pft.senesc_scal_water * pow(1.0 - indiv.get_wscal(),pft.qwater_senesc);
	
	// Senescence due to freezing and hot temperatures 
	double senesc_ftemp = 0.0;
	double Tfreeze = 2.0;  // mean daily temperature below which freezing is assumed to occur (see comment above)
    double Tcold = max(0.0, Tfreeze - grass_temp);               // freezing stress
	double Theat = max(0.0, grass_temp - indiv.pft.T0_growth);   // heat stress
	double senesc_freeze_scal = indiv.pft.pathway == C3 ? max(0.0,min(1.0,1.0 - 0.25 * (pft.temp_dormancy - indiv.temp_dormancy))) : 1.0;   // less freezing stress for C3s
	senesc_ftemp = pft.senesc_scal_temp * (senesc_freeze_scal * (1.0 - exp(-Tcold / 3.0)) + (1.0 - exp(-pow(Theat / 12.0, 2.5))));

	// Dormancy-induced
	double senesc_fdorm = 0.0;
	if (indiv.fphen_activity < 1.0 && indiv.dphen_activity < 0.0) {
		// could include dphen_activity here (i.e. stronger senescence with more stress)
        // senesc_fdorm = pft.senesc_scal_dorm * (1.0 - ((indiv.fphen_activity - fphen_activity_min) / (1.0 - fphen_activity_min)));
		senesc_fdorm = pft.senesc_scal_dorm * (1.0 - indiv.fphen_activity);
	}

	// Self-shading
	double kshade = 6.5;
	double senesc_fshade = kshade * indiv.fpc;
	double senesc_fshade_norm = senesc_fshade / (1.0 + senesc_fshade);

	// Total senescence rate
	//senesc *= ((1.0 + max(indiv.fdev_senesc,senesc_fdorm)) * (1.0 + senesc_fwater) * (1.0 + senesc_ftemp));
	senesc *= (1.0 + indiv.fdev_senesc * (1.0 - 0.9 * senesc_fshade_norm) + senesc_fdorm + senesc_fwater) 
	           * (1.0 + senesc_ftemp) 
			   * (1.0 + senesc_fshade);
	senesc = min(1.0, senesc);
	indiv.senesc_root = pft.senesc_base;

	// Retranslocation and pool updates
	// TODO: maybe call cstorage_capacity once at the beginning to be consistent across routines!
	// TODO: calculate retranslocation coefficients based on current C:N ratios?

	// Carbon
	double csenesc  = indiv.cmass_leafs[green] * senesc;
	double cstore   = csenesc * pft.crelocfrac;
	double ctobrown = csenesc * (1.0 - pft.crelocfrac);
	indiv.csenesc            =  csenesc;

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
void senescence_dave_tree(Individual &indiv) {
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
			senescence_dave_grass(indiv);

			assert(positive(indiv.cmass_leafs[green]));
			assert(positive(indiv.nmass_leafs[green]));
			assert(positive(indiv.cmass_leafs[brown]));
			assert(positive(indiv.nmass_leafs[brown]));
		} else {
			senescence_dave_tree(indiv);
		}

		assert(positive(indiv.get_cmass_leaf()));
		assert(positive(indiv.get_nmass_leaf()));
		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));
	}
}

/*
Check if this PFT is established (and out of probation) in this patch.
*/
bool is_established(const Pft &pft, const Vegetation &vegetation) {
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		if (indiv.pft.id == pft.id) {
			return true;
		}
	}
	return false;
}

// calculate the minimum cmass that leads to an establishment of lai >= lai_min
double get_cmass_init_min(Patchpft& ppft) {

	return lai_min / (ppft.pft.sla * (1.0 - cstorage_max_temp) * ppft.pft.ltor_max / (1.0 + ppft.pft.ltor_max));

}



// External seed input 
// This should only be called if the PFT is not yet established at that location.
// Otherwise PFTs need to rely on their local seedbank.
void seed_dispersal_dave(Vegetation& vegetation) {

	Patch& patch = vegetation.patch;
	for (int i = 0; i < patch.pft.nobj; i++) {
		Patchpft &ppft = patch.pft[i];

		double cmass_init_min = get_cmass_init_min(ppft);

		if (!ppft.pft.dave_allocation() || is_established(ppft.pft, vegetation) ||
			ppft.cmass_litter_repr > cmass_init_min) {
			// This pft is not using dave allocation, and is therefore not using
			// the daily establishment. This means the seedbank and seed
			// dispersal will have no effect on this pft.

			// This routine is also not used (i.e. no external seeds added)
			// if either an individual of the PFT is alive or it still 
			// has enough seeds in the seedbank to establish a plant > lai_min.
			continue;
		}

		// divide the external seed influx by the number of PFTs simulated
		// to not inflate N influx with higher number of PFTs.
		double ext_seed_c = ext_seed / patch.pft.nobj;

		// test for initial 
		if (date.year == 0 && ppft.pft.lifeform == TREE) {
			ext_seed_c *= 10;
		}

		double ext_seed_n = ext_seed_c / ((ppft.pft.cton_leaf_avr + ppft.pft.cton_root_avr) / 2.0);

		ppft.cmass_ext_seed += ext_seed_c;
		ppft.nmass_ext_seed += ext_seed_n;
		//ppft.cmass_litter_repr += ext_seed_c;
		//ppft.nmass_litter_repr += ext_seed_n;

        //patch.fluxes.report_flux(Fluxes::ESTC, -ext_seed_c);
		//patch.fluxes.report_flux(Fluxes::ESTN, -ext_seed_n);
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

		// Turnover rate of seeds

		// Air temperature or bottom snow layer if snow present and soil moisture of first soil layer.
		double surf_temp = vegetation.patch.soil.T_soil[6];
		double surf_water = patch.soil.get_layer_soil_water(0);
		
		// Base turnover rate modified by temperature and water availability
		double ftemp = 0.0, fwater = 0.0;
		double base_turnover_seed = ppft.pft.turnover_repr_litter / date.year_length();
		respiration_temperature_response(surf_temp, ftemp); 
		fwater = 1.0 / (1.0 + 30.0 * exp(-8.5 * surf_water));
		double frac_turnover_seed = base_turnover_seed * ftemp * fwater;

		// Turnover of the actual seedbank

		// Carbon
		double cturnover_repr = frac_turnover_seed * ppft.cmass_litter_repr;
		ppft.cmass_litter_repr -= cturnover_repr;
		ppft.cmass_litter_root += cturnover_repr;

		// Nitrogen
		double nturnover_repr = frac_turnover_seed * ppft.nmass_litter_repr;
		ppft.nmass_litter_repr -= nturnover_repr;
		ppft.nmass_litter_root += nturnover_repr;

		
		// Turnover of the hypothetical external seedbank that is used for establishment
		// of new PFTs. Note that this pool does not turn over into litter but 
		// simply ceases to exist just like it appeared out of nowhere.

		// Carbon
		double cturnover_ext_seed = frac_turnover_seed * ppft.cmass_ext_seed;
		ppft.cmass_ext_seed -= cturnover_ext_seed;
		//ppft.cmass_litter_root += cturnover_ext_seed;

		// Nitrogen
		double nturnover_ext_seed = frac_turnover_seed * ppft.nmass_ext_seed;
		ppft.nmass_ext_seed -= nturnover_ext_seed;
		//ppft.nmass_litter_root += nturnover_repr;

	}
}


/*
@param ctolitter: Cumulative C sent to litter
@param ntolitter: Cumulative N sent to litter
*/
void turnover_pool(double turnover_frac, double crelocfrac, double nrelocfrac,
	double &cpool, double &npool, double &cdest, double &ndest,
	double &ctolitter, double &ntolitter, double &initialc, double &initialn,
	double &cmass_storage, double &nmass_storage, double cstorage_max, double nstorage_max,
	Patch &patch) {

	if (cpool > 0.0) {

		// Carbon

		double cturnover = cpool * turnover_frac;

		// Send some C to storage and some to litter.
		double cstore = cturnover * crelocfrac;
		double clitter = cturnover * (1.0 - crelocfrac);

		if ( (cmass_storage + cstore) > cstorage_max) {
			// Storage allocation would cause storage to exceed capacity.
			// Reduce cstore to fill c storage and send the remainder to litter.
			cstore = cstorage_max - cmass_storage;
			clitter = cturnover - cstore;
		}

		// When this individual was established, some biomass appeared out of
		// thin air (out of the void). This amount is initialc/initialn. We
		// don't want this going into the soil, because this would result in
		// long-term fertilisation of the system. Therefore, we reduce the
		// biomass sent to litter until this imaginary biomass is returned to
		// the void.
		double cvoid = 0;
		if (initialc > 0) {
			cvoid = min(initialc, clitter);
			clitter -= cvoid;

			// This C entered the system as a negative establishment flux.
			patch.fluxes.report_flux(Fluxes::ESTC, cvoid);
			initialc -= cvoid;
		}

		cpool      -= cturnover; // Split into clitter and cstore
		cdest      += clitter;
		ctolitter  += clitter; // Cumulative across pools
		cmass_storage += cstore;

		assert(equal(cturnover, cstore + clitter + cvoid));
		assert(positive(cpool));
		ctolitter += cvoid; // temporary hack for mass balance checks


		// Nitrogen

		double nturnover = npool * turnover_frac;

		// Send some N to storage and some to litter.
		double nstore  = nturnover * nrelocfrac;
		double nlitter = nturnover * (1.0 - nrelocfrac);

		if (nmass_storage + nstore > nstorage_max) {
			// Storage allocation would cause storage to exceed capacity.
			// Reduce nstore to fill nstorage and send the remainder to litter.
			nstore  = nstorage_max - nmass_storage;
			nlitter = nturnover - nstore;
		}

		// See comment on cvoid for explanation.
		double nvoid = 0;
		if (initialn > 0) {
			nvoid = min(initialn, nlitter);
			nlitter -= nvoid;

			// This N entered the system as a negative establishment flux.
			patch.fluxes.report_flux(Fluxes::ESTN, nvoid);
			initialn -= nvoid;
		}

		// Move nitrogen between pools.
		npool           -= nturnover; // Split into nlitter and nstore
		ndest           += nlitter;
		ntolitter 		+= nlitter; // Cumulative between pools
		nmass_storage   += nstore; // Move N to storage

		assert(equal(nturnover, nstore + nlitter + nvoid));
		ntolitter += nvoid; // temporary hack
	}
}

// Note: cropallocation.cpp contains a function called turnover_grass(). I am
// therefore avoiding reusing this name here.
void turnover_dave_grass(Individual &indiv) {
	Vegetation &vegetation = indiv.vegetation;
	Patch& patch = vegetation.patch;
	Pft& pft = indiv.pft;
	Patchpft &ppft = patch.pft[indiv.pft.id];
	Soil &soil = patch.soil;

	// Determine turnover rate (same for all individuals in a PFT).
	// Assumed daily rates here (fraction d-1)!
	// Air temperature or bottom snow layer if snow present and soil moisture of first soil layer.
	double surf_temp = vegetation.patch.soil.T_soil[6];
	double surf_water = patch.soil.get_layer_soil_water(0);

	// crown turnover relative to root turnover.
	double scal_crown_turnover = 0.35;

	// Turnover of living tissues (root and crown)
	// senesc_root equals the baseline senescence rate of leaves.
	double frac_turnover_root = 30.0 * indiv.senesc_root * indiv.fphen_activity;	
	double frac_turnover_crown = scal_crown_turnover * frac_turnover_root;

	// Turnover of brown shoots (same temperature and moisture responses as for seeds)
	double ftemp = 0.0, fwater = 0.0;
	double kleaf = 0.01;  // base respiration rate of brown shoots at 10degC
	// Temperature response
	respiration_temperature_response(surf_temp, ftemp);  // Lloyd & Taylor 1994
	// Soil moisture response (Parton et al. 2001). 
	// Note that Parton et al. 2001 used PPT/PET instead of wcont.
	//fwater = pow(surf_water,1.0);
	fwater = 1.0 / (1.0 + 30.0 * exp(-8.5 * surf_water));
	double frac_turnover_leaf = kleaf * ftemp * fwater;

	// Dead shoot turnover (brown shoot --> leaf litter)
	double cturnover_leaf = 0;
	double nturnover_leaf = 0;

	// Work out how much goes to litter.
	turnover_pool(frac_turnover_leaf, 0, 0, indiv.cmass_leafs[brown],
                  indiv.nmass_leafs[brown], cturnover_leaf,
                  nturnover_leaf, indiv.ctolitter, indiv.ntolitter,
                  indiv.initialc, indiv.initialn, indiv.cmass_storage,
                  indiv.nmass_storage, indiv.cstorage_max, indiv.nstorage_max, patch);

	// Send leaf turnover to litter.
	// This will go directly into the SURFSTRUCT and SURFMETA SOM pools.
	leaf_to_litter(soil, pft, cturnover_leaf, nturnover_leaf);
	indiv.cturnover_leaf = cturnover_leaf;
	indiv.nturnover_leaf = nturnover_leaf;

	// Crown turnover (goes into root litter).
	// TODO: put into function (same for roots and crowns)!
	double cturnover_root = 0;
	double nturnover_root = 0;
	turnover_pool(frac_turnover_crown, 0, nrelocfrac, indiv.cmass_crown,
                  indiv.nmass_crown, cturnover_root,
                  nturnover_root, indiv.ctolitter, indiv.ntolitter,
                  indiv.initialc, indiv.initialn, indiv.cmass_storage,
                  indiv.nmass_storage, indiv.cstorage_max, indiv.nstorage_max, patch);
    indiv.cturnover_crown = cturnover_root;
	indiv.nturnover_crown = nturnover_root;

	// Root turnover (roots --> root litter)
	turnover_pool(frac_turnover_root, 0, nrelocfrac, indiv.cmass_root,
                  indiv.nmass_root, cturnover_root,
                  nturnover_root, indiv.ctolitter, indiv.ntolitter,
                  indiv.initialc, indiv.initialn, indiv.cmass_storage,
                  indiv.nmass_storage, indiv.cstorage_max, indiv.nstorage_max, patch);
	indiv.cturnover_root = cturnover_root - indiv.cturnover_crown;
	indiv.nturnover_root = nturnover_root - indiv.nturnover_crown;

	// Send crown and root turnover to litter.
	// This will go directly into the SOILSTRUCT and SOILMETA SOM pools.
	root_to_litter(soil, pft, cturnover_root, nturnover_root);

	// Reproductive turnover
	//double repr_turnover_time = 75.0;    // d
	//if (indiv.growth_stage == reproductive && indiv.cmass_repr > indiv.cmass_crown){
	//	indiv.flush_repr_litter_DAVE(0.01);
	//}

}

// calculate the coefficient of water stress for the turnover of leaves and roots in trees
double calculate_fwater_turnover_tree(const Individual& indiv) {
	// fraction of change in leaf and root turnover as a result of water stress.
	// leaf (root) turnover will increase (decrease) by the same fraction 
	// in high soil water stress conditions.
	// TODO: function seems too extereme (too much turnover for low water stress).
	// // Option 1 (older version)
	// double fwater = 1 + (senesc_scal_water * pow(( 1- ppft.wscal_30.mean()), qwater_senesc)); // OLD
	// // Option 2 (newer version of the same function)
	// double fwater = max( 1.0, senesc_scal_water * ( 1 - pow(ppft.wscal_30.mean(), qwater_senesc) ) ); // newer
	// // Option 3
	const double maxy = 4.0;
	const double miny = 1.0;
	const double fwater_expo = 3.0;
	const double fwater_mid = 0.1;
	double fwater_min = (miny +((maxy - miny) / (1 + pow((1.0 / fwater_mid), fwater_expo))));
	double fwater = miny + ((maxy - miny) / (1 + pow((indiv.get_wscal_30().mean() / fwater_mid), fwater_expo)));
	fwater /= fwater_min;
	fwater = 1.0;
	return fwater;
}

// linear function between minimum and maximum turnover as a function of growth potnetial.
// The idea being that turonover is reduced to a minimum when growth conditions are 
// not good. This should potentially create high growth demand int the growing season
// by increasing the turnover.
double calculate_fcgrow_tree(double min_turnover_fraction, double cgrow_pot) {

	double fturnover_max = 1.0;
	// shape of the exponential curve between the turnover coefficient to cgrow_pot (ie growth potential)
	const double shapefactor = 1;

	double xmax = 0.2;

	// fcgrow_tree
	return  min( 1.0, min_turnover_fraction + pow(((fturnover_max - min_turnover_fraction) / xmax * cgrow_pot), shapefactor) ) ;
	// return  fturnover_min + pow(((fturnover_max - fturnover_min) / xmax * cgrow_pot), shapefactor);

}


// Turnover of green pools to litter. 
// TODO: If only using a single leaf pool, we could remove the senescence function altogether.
void turnover_dave_tree(Individual &indiv) {
	// Inputs:  indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap
	//          indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap,
	//          indiv.pft.leaflong
	//
	// Outputs: indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap
	//          indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap
	//          ppft.cmass_litter_leaf, ppft.cmass_litter_root, indiv.cmass_heart

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
	

	// NEW: scalar to scale the turnover across the pools with the phenological state,
	// which is represented via cgrow_pot (function of temeprature and soil moisture).
	// The scalar has a minimum, which is 0.1 * turnover_pool, assuming a 10% of the
	// turnover occuring when temperature (in the temperrate zone) and soil mositure 
	// (in areas with a distinct dry/wet season).
	// proportion of daily turnover that sets the minimum turnover for a specific pool
	// TODO: turn into pft parameter (if it works)
	double min_tunrover_proportion = 0.1;
	// Calculate the proportion of the daily turnover as a function of growth condition 
	// this variable is actually alows to isolate the effect of growth conditions on
	// turnover (and not on growth) .
	double growth_lim = min(indiv.fphoto, min(indiv.fwater, indiv.ftemp));

	double fcgrow = calculate_fcgrow_tree(min_tunrover_proportion, growth_lim);
	// fcgrow = 1.0;

	// if (growth_lim == 0.0) {
	// 	fcgrow = 0.0;
	// }
	// bypass phenology
	if (!iftreephen) {
		fcgrow = 1.0;
	}
	// fcgrow = pow(growth_lim, 4.0);
	Soil &soil = indiv.vegetation.patch.soil;
	Patchpft& ppft = indiv.patchpft();
	const Pft &pft = indiv.pft;
	
	// fraction of change in leaf and root turnover as a result of water stress.
	// leaf (root) turnover will increase (decrease) by the same fraction 
	// in high soil water stress conditions.
	double fwater = calculate_fwater_turnover_tree(indiv);

	double cturnover_leaf = 0.0, nturnover_leaf = 0.0;
	double cturnover_root = 0.0, nturnover_root = 0.0;
	double cturnover_sap = 0.0, nturnover_sap = 0.0;
	double cturnover_repr = 0.0, nturnover_repr = 0.0;

	double nretranslocated = 0;

	// scenescence of leaves (also depending on water stress (fwater)). 
	// TODO: 1. no regards to what this would do to actual k_latosa tomorrow; 
	// 2. Should the base turnover (pft.turnover_leaf) be more related to leaf longevity?; 
	// 3. leaf turnover shoudl be turned off in very young trees (until the leaves 
	// reach leaf longevity). Alternatively, reduce base leaf turnover by a fraction
	// in the first few years of the tree. For example, in the first years, turnover_leaf
	// will be 1/2 of it's value when the tree is older than leaf longevity:
	// (e.g: turnover_leaf = indiv.pft.turnover_leaf / 2 * fwater / Date::MAX_YEAR_LENGTH)
	if (indiv.cmass_leaf > 0) {		
		// scales leaf turnover according to mean leaf age
		double fleafage = 1.0;

		// fleafage = 1.0 + indiv.meanleafage / indiv.pft.leaflong;
		double leaflong = 1.0 / pft.turnover_leaf;
		if (iftreephen) {
			fleafage = indiv.meanleafage / leaflong;
		}
		// fleafage = 1.0 + indiv.meanleafage / indiv.pft.leaflong;
		// fleafage = 0.5 + 0.002 * pow(indiv.meanleafage, 5.0);
		// fleafage = 1.0;

		// the proportion of the turnover that is implemented when cgrow is 0
		// TEST (ASSAF): turn over a minimum amount if there is no growth today
		// maximum daily leaf turnover
		double maximum_daily_turnover_leaf = (indiv.pft.turnover_leaf / Date::MAX_YEAR_LENGTH);

		// // // Option 1 - set turnover to a minimum when cgrow == 0
		// // // fraction of leaf C that is turned over today. Tuned down to allow lower
		// // // when cgrow is 0, or when the tree prioritises storge over growth.
		// turnover_leaf = maximum_daily_turnover_leaf * min_tunrover_leaf_proportion;
		// if (indiv.cgrow > 0) {
		// 	turnover_leaf = fleafage * fwater * indiv.pft.turnover_leaf / Date::MAX_YEAR_LENGTH;
		// }

		// // Option 2 - sacaling of leaf turnover according to cgrow
		// // Scale all leaf turnover according to cgrow potential scalar (indiv.cgrow_pot; 
		// // which is between 0-1 and is affected by temperature, soil moisture and the storage state). 
		// // This formulation will scale all leaf turnover between a minimum and leaf turnover.
		// // Reasoning: turnover is very low when the trees are not growing leaves and storing carbon.
		// // For example, in drought, the trees want to keep their leaves to maximise photosynthesis, 
		// // but up to a point where keeping the leaves is a liability in terms of water.
		// double turnover_leaf_max = indiv.pft.turnover_leaf / Date::MAX_YEAR_LENGTH;
		// double turnover_leaf_min = turnover_leaf_max * min_tunrover_leaf_proportion;	
		// turnover_leaf = fleafage * fwater * (turnover_leaf_min + ((turnover_leaf_max - turnover_leaf_min) * indiv.cgrow_pot));

		turnover_leaf = /*fwater * fleafage * fcgrow * */ maximum_daily_turnover_leaf;
		// indiv.rel_turnover_leaf = turnover_leaf / turnover_leaf_max;
		// Leaf C amount to turnover (kgC/m^2)
		cturnover_leaf = turnover_leaf * indiv.cmass_leaf;
		// // test (ASSAF): increasing the turnover significantly when mean leaf age is
		// // higher thean the leaf longevity. Purpose: to increase the alpha_leaf
		// // in the allocation 
		// if (indiv.meanleafage > indiv.pft.leaflong) {
		// 	// option 1: reduce the carbon above leaf longevity (approx)
		// 	cturnover_leaf = indiv.cmass_leaf * (indiv.meanleafage - indiv.pft.leaflong);
		// 	turnover_leaf = cturnover_leaf / indiv.cmass_leaf;
		// }

		// save turnover (TODO: remove the local turnover_leaf variable and use this instead)
		indiv.pturnover_leaf = turnover_leaf;
		// Leaf N amount to turnover (kgN/m^2)
		nturnover_leaf = turnover_leaf * indiv.nmass_leaf;
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

		// Update pools
		indiv.ctolitter += cturnover_leaf;
		indiv.cmass_leaf -= cturnover_leaf;

		indiv.ntolitter += n2litter_leaf;
		indiv.nmass_leaf -= nturnover_leaf;

		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_leaf;
		nretranslocated += n2storage_leaf;

		// ppft.cmass_litter_leaf += cturnover_leaf;
		// ppft.nmass_litter_leaf += n2litter_leaf;
		// Let somdynam decide what to do with leaf litter.
		leaf_to_litter(soil, pft, cturnover_leaf, n2litter_leaf);

		assert(equal(nturnover_leaf, n2storage_leaf + n2litter_leaf));
	}
	// save daily leaf turnover flux (kgC/m^2/day)
	indiv.cturnover_leaf = cturnover_leaf;
	indiv.nturnover_leaf = nturnover_leaf;

	if (indiv.cmass_root > 0) {
		// fraction of root biomass that is turned over today. 
		turnover_root = /*fcgrow * 1.0 / fwater */ indiv.pft.turnover_root / Date::MAX_YEAR_LENGTH;
		// root C amount to turnover (kgC/m^2)
		cturnover_root = turnover_root * indiv.cmass_root; 
		// root N amount to turnover (kgN/m^2)
		nturnover_root = turnover_root * indiv.nmass_root;
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

		// Update pools
		indiv.ctolitter += cturnover_root;
		indiv.cmass_root -= cturnover_root;

		indiv.ntolitter += n2litter_root;
		indiv.nmass_root -= nturnover_root;

		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_root;
		nretranslocated += n2storage_root;

		// ppft.cmass_litter_root += cturnover_root;
		// ppft.nmass_litter_root += n2litter_root;
		// Let somdynam decide what to do with repr litter (treating it as leaf
		// litter).
		root_to_litter(soil, pft, cturnover_root, n2litter_root);

		assert(equal(nturnover_root, n2litter_root + n2storage_root));
	}
	// save daily root turnover flux (kgC/m^2/day)
	indiv.cturnover_root = cturnover_root;
	indiv.nturnover_root = nturnover_root;

	if (indiv.cmass_sap > 0) {

		// TODO: CHECK!
		// if there isn't any leaf turnover today, then there shouldn't be sapwood 
		// turnover becasue current sapwood (cmass_sap) is needed to maintain the 
		// current leaves.
		double fsapage = 1.0;
		if (iftreephen) {
			fsapage = indiv.meansapage / (1.0 / indiv.pft.turnover_sap);
		}
		if ((turnover_leaf > 0) && (indiv.pft.k_latosa > 0)) {
			turnover_sap =  /* fsapage * fcgrow * */ indiv.pft.turnover_sap / Date::MAX_YEAR_LENGTH;
		} else {
			// no turnover sap today
			turnover_sap = 0;
		}

		// sapwood C amount to turnover (kgC/m^2)
		cturnover_sap = turnover_sap * indiv.cmass_sap; 
		// sapwood N amount to turnover (kgN/m^2)
		nturnover_sap = turnover_sap * indiv.nmass_sap;
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

		// Update pools
		indiv.cmass_heart += cturnover_sap;
		indiv.cmass_sap -= cturnover_sap;

		indiv.nmass_heart += n2heart_sap;
		indiv.nmass_sap -= nturnover_sap;

		// Update N storage (kgN/m^2)
		indiv.nmass_storage += n2storage_sap;
		nretranslocated += n2storage_sap;

		assert(equal(nturnover_sap, n2storage_sap + n2heart_sap));
	}

	if (indiv.cmass_repr > 0) {
		double turnover_repr = indiv.pft.turnover_repr / Date::MAX_YEAR_LENGTH;

		cturnover_repr = turnover_repr * indiv.cmass_repr;
		nturnover_repr = turnover_repr * indiv.nmass_repr;

		// Update pools
		indiv.cmass_repr -= cturnover_repr;
		indiv.nmass_repr -= nturnover_repr;

		indiv.ctolitter += cturnover_repr;
		indiv.ntolitter += nturnover_repr;

		// Let somdynam decide what to do with reproductive litter.
		// TODO: should repr be treated as leaf or root litter?
		// For now, use root (this reproduces previous behaviour).
		root_to_litter(soil, pft, cturnover_repr, nturnover_repr);
	}

	indiv.cturnover_repr = cturnover_repr;
	indiv.nturnover_repr = nturnover_repr;

	// save daily sapwood turnover flux (kgC/m^2/day)
	indiv.cturnover_sap = cturnover_sap;
	indiv.nturnover_sap = nturnover_sap;
	indiv.nretranslocated = nretranslocated;
	indiv.anretranslocated += nretranslocated;
	// // drop reproductive pool when trends goes down
	// // TODO: we need a better solution here to account for flushing of the 
	// // reproductiove pool into the litter. possible do a trend rather than all 
	// // at once.
	// if (indiv.cgrow_repr_20.periodictrend(20) < 0) {
	// 	indiv.flush_repr_litter_DAVE(1.0);
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

	Soil &soil = vegetation.patch.soil;

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
			turnover_dave_tree(indiv);
		}
	}

	turnover_seed(vegetation);
}

// Allocate excess C and N to target pool.
void allocate_excess(double &cexcess, double &nexcess, double &cpool_dest, double &npool_dest,
					 double &cmass_storage, double &nmass_storage, double max_mobil_cstorage,
					 double cton_min, double cton_max, double calloc_max) {

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
	if (indiv.pft.lifeform == GRASS_ANNUAL && indiv.growth_stage == reproductive) {
		// TODO JK: limit allocation to repr for annuals to avoid unrealistic growth rates? 
		// TODO JK: refine frac to be more consistent with the rest of the code!
		double frac = indiv.max_mobil_cstorage * indiv.cmass_repr / (indiv.cmass_repr + indiv.get_cmass_leaf()); 
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
							indiv.nmass_storage, indiv.max_mobil_cstorage, pft.cton_repr_min, 
							pft.cton_repr_max, cexcess_repr);	
		}

		// Put as much of possible of the rest of C excess into the crown (perennials only).
		double clim_crown = indiv.pft.lifeform == GRASS_PERENNIAL ? 0.5 : 0.0;
		double cexcess_crown = min(indiv.cexcess, clim_crown * indiv.cgrow);
		allocate_excess(indiv.cexcess, indiv.nexcess, indiv.cmass_crown, indiv.nmass_crown, indiv.cmass_storage, 
						indiv.nmass_storage, indiv.max_mobil_cstorage, pft.cton_root_min, 
						pft.cton_root_max, cexcess_crown);

	}

	// Now deal with excess C and N that cannot be allocated because of stoichiometry or
	// because they would lead to excessive growth.
	if (indiv.cexcess > 0.0) {
		
		// part of it is respired and the other part goes to soil litter as root exudates
		double fexcess_respire = 0.3;
		double cexcess_to_resp = fexcess_respire * indiv.cexcess;
		double cexcess_to_soil = indiv.cexcess - cexcess_to_resp;

		// 1) autotrophic respiration
		indiv.dnpp -= cexcess_to_resp;
		indiv.aresp += cexcess_to_resp;

		indiv.report_flux(Fluxes::RA, cexcess_to_resp);
		indiv.report_flux(Fluxes::NPP, -cexcess_to_resp);

		// 2) root exudates / rhizodeposition (approximated by input to root litter)
		ppft.cmass_litter_leaf += cexcess_to_soil;
		indiv.ctolitter += cexcess_to_soil;

		// Reset cexcess to 0.0
		indiv.cexcess = 0.0;
	}

	if (indiv.nexcess > 0.0) {
		// flush to litter
		// TODO: return this to the available N pool (e.g. via root exudates).

		double nvoid = min(indiv.nexcess, indiv.initialn);
		double nexcess = indiv.nexcess - nvoid;

		ppft.nmass_litter_root += nexcess;
		indiv.ntolitter += indiv.nexcess;

		// This N entered the system as a negative establishment flux.
		indiv.vegetation.patch.fluxes.report_flux(Fluxes::ESTN, nvoid);
		indiv.initialn -= nvoid;
	}
}


// Photosynthetic sink limitation.
// Calculates downregulation factors applied to photosynthetic parameters.
void sink_source_feedback(Individual& indiv) {
	
	//double dpslim = 0.00;
	//double qf2 = 100.0;
	//double cstorage_fill = Cstorage_fillfrac(indiv);
	double slim_max = 0.0;  // if 0.0, feedback is turned off.

	//indiv.ps_sink_lim = min(1.0, (1.0 + dpslim - 2.0 * dpslim * pow(indiv.fphen_activity*cstorage_fill,qf2)) * indiv.ps_sink_lim);
	double excess_strength = 0;
	if (indiv.cexcess_npp_5.size() > 1) {
		indiv.ps_sink_lim = max(0.0, 1.0 - min(indiv.cexcess_npp_5.mean(),slim_max));
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

		// TODO JK: make sla age-dependent or make culmfrac dynamic

		// Update capacity of C and N storage
		indiv.cstorage_max = cstorage_capacity(indiv);
		indiv.nstorage_max = nstorage_capacity(indiv);
		assert(indiv.cstorage_max > 0);
		assert(indiv.nstorage_max > 0);

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
			indiv.cstorage_dynam_14.add(0.0);
		} else {
			indiv.cmass_storage_20.add(indiv.cmass_storage);
			indiv.cstorage_dynam_14.add(indiv.cmass_storage_dynam);

		}
		indiv.cmass_repr_10.add(indiv.cmass_repr);
		indiv.frac_germ = 0.0;  // reset in vegdynam.cpp if > 0.0

		indiv.cexcess_npp_5.add(indiv.cexcess / max(1.0e-12,indiv.dnpp));

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
	// added the form factor here to PI. It belongs im the next equations (Fr and Fh)
	// since they are derived from the derivative of the equation of the mass of a cylinder (ie:
	// Wcyl = PI * (R^2) * H * wooddens, once for change in radius with constant 
	// height, and once for change in height for constnat  radius )

	// Added the form factor here to PI. It belongs im the next equations (Fr
	// and Fh) since they are derived from the derivative of the equation of the
	// mass of a cylinder (ie: Wcyl = PI * (R^2) * H * wooddens, once for change
	// in radius with constant height, and once for change in height for
	// constnat  radius).
	double PI_FF = PI * indiv.pft.formfactor;

	double Fr = PI_FF * wooddens * 2.0 * R * h;
	double Fh = PI_FF * wooddens * R * R;
	double Fallom = k_allom2 * pow(2.0, k_allom3) * k_allom3 * pow(R, k_allom3 - 1.0);
	// should calc once: k_allom2 * pow(2.0,k_allom3) * k_allom3
	// we need to calc change for individual trees in cohort
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

		// save increment of dimentsion (m)
		indiv.d_height = dh;
		indiv.d_diam = 2.0 * dr;

		double basalarea_t1 = PI * pow((indiv.diam / 2.0), 2.0);

		// add increments to dimentions
		indiv.height += indiv.d_height;
		indiv.diam += indiv.d_diam;
		
		// calculate annual (cumulative) diameter incremnt (reset at Jan 1st of every year)
		if (date.day == 0){
			indiv.c_d_diam = 0.0;
		}
		indiv.c_d_diam += indiv.d_diam;
		
		// new basal area (m2)
		indiv.basalarea = PI * pow((indiv.diam / 2.0), 2.0);
		
		// basal area increment (m2/day)
		indiv.basalarea_inc = indiv.basalarea - basalarea_t1;
	
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
		if (ifforcedlai) {
			indiv.lai = fixed_lai;
		} else {
			indiv.lai = indiv.cmass_leaf * indiv.pft.sla;
		}

		// calculate leaf to sapwood area (unitless vatrible)
		if (!is_zero(indiv.height)) {
			double const LATOSA_MAX = 8000.0;
			// sapwood area (m^2)
			indiv.sa = (indiv.cmass_sap) / (indiv.height * indiv.pft.formfactor * indiv.pft.wooddens);
			if (!is_zero(indiv.sa)) {
				// leaf to sapwood area (unitelss)
				indiv.latosa = indiv.lai / indiv.sa;
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

	// calculate new root projected cover, but aassumes that rpc does not change 
	// if lai decreases due to fire/drought). 
	// TODO: in the future
	// 1. set this to be a function of root biomass/volume to increase competativeness
	// 2.to change when we decide that the allocation prority to roots that affects 
	// rpc to be pft-specific.
	indiv.rpc = max(indiv.rpc, calc_rpc_tree(indiv));

	// Bole height (Eqn 7)
	indiv.boleht = calc_bole_height(indiv);

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

		if (indiv.lai_today() > indiv.mlai_max[date.month]) {
			indiv.mlai_max[date.month] = indiv.lai_today();
		}

		// This is the same as indiv.lai for grasses, but it includes phen for trees/mosses.
		indiv.mlai[date.month] += indiv.lai_today() / (double)date.ndaymonth[date.month];
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
		indiv.age_active += 1.0 / Date::MAX_YEAR_LENGTH * indiv.fphen_activity;

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

			if (indiv.pft.lifeform == TREE) {
				indiv.flush_repr_litter_DAVE(1.0);
			}
		} // Last day of year
	}

	// Reset nday_active to 0 at start of year for all PFTs, regardless of
	// whether they're currently established in this patch, to avoid nday_active
	// from carrying over from one year to the next.
	for (size_t i = 0; i < pftlist.nobj; i++) {
		Patchpft &ppft = vegetation.patch.pft[i];
		if (date.day == 0) {
			ppft.dormancy_depth = 0.0;
			ppft.wscal_cum = 0.0;
			ppft.ndays_alive = 0;
		}
	}
}

void growth_daily_accounting_dave(Vegetation &vegetation) {
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		if (!indiv.pft.dave_allocation()) {
			continue;
		}

		// Set indiv.cgrow, indiv.cexcess, and indiv.nexcess to 0.
		indiv.cgrow                  = 0.0;
		indiv.cexcess                = 0.0;
		indiv.nexcess                = 0.0;
		indiv.nmass_missing          = 0.0;
		indiv.ctolitter              = 0.0;
		indiv.ntolitter              = 0.0;
		indiv.ntoavail               = 0.0;
		indiv.cturnover_leaf         = 0.0;
		indiv.nturnover_leaf         = 0.0;
		indiv.cturnover_root         = 0.0;
		indiv.nturnover_root         = 0.0;
		indiv.cturnover_repr         = 0.0;
		indiv.nturnover_repr         = 0.0;
		indiv.cturnover_crown        = 0.0;
		indiv.nturnover_crown        = 0.0;
		indiv.cturnover_sap          = 0.0;
		indiv.nturnover_sap          = 0.0;
		indiv.absolute_alpha_leaf    = 0.0;
		indiv.absolute_alpha_root    = 0.0;
		indiv.absolute_alpha_sap     = 0.0;
		indiv.absolute_alpha_storage = 0.0;
		indiv.cgrow_leaf		     = 0.0;
		indiv.cgrow_root		     = 0.0;
		indiv.cgrow_sap		         = 0.0;
		indiv.cgrow_repr		     = 0.0;
		indiv.cgrow_storage		     = 0.0;
		indiv.cgrow_crown		     = 0.0;
		indiv.csenesc                = 0.0;

		if (indiv.pft.is_grass()) {
			assert(positive(indiv.cmass_leafs[brown]));
		}

		if (date.day == 0) {
			indiv.acalloc_leaf = 0;
			indiv.acalloc_root = 0;
			indiv.acalloc_sap = 0;
			indiv.acalloc_repr = 0;
			indiv.acalloc_storage = 0;
			indiv.acalloc_crown = 0;
			indiv.anretranslocated = 0;
		}
	}
}

void growth_daily_dave(Patch& patch) {

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

} // growth_daily_dave()

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
		growth_daily_dave(patch);
	} else if (ifdailygrass) {
		// David Warlind daily allocation.
		growth_daily_grass_warlind(patch);
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
