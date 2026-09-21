///////////////////////////////////////////////////////////////////////////////////////
/// \file vegdynam.cpp
/// \brief Vegetation dynamics and disturbance
///
/// \author Ben Smith
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
#include "vegdynam.h"
#include "growth.h"
#include "driver.h"
#include "logging.h"
#include "landcover.h"
#include "management.h"
#include "blaze.h"

// Alternatives for forestry:

/// Whether to use smaller sapsize in managed forests after cutting (defined by PLANTSIZE) when not using planting()
const bool SMALL_SAPSIZE_POST_CUT = false;

/// Fixed initial individual plantsize (kgC) in tree planting after clearcut if SMALL_SAPSIZE_POST_CUT == true
const double PLANTSIZE = 0.25;

/// Upper LAI limit for wetland species. No limit: 0 //TODO remove this after daily allocation.
const double wetlandlailimit = 2.0;

/// Internal help function for splitting up nitrogen fire fluxes into components
void report_fire_nfluxes(Fluxes& fluxes, double nflux_fire) {
	fluxes.report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire);
	fluxes.report_flux(Fluxes::NOx_FIRE, Fluxes::NOx_FIRERATIO * nflux_fire);
	fluxes.report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);
	fluxes.report_flux(Fluxes::N2_FIRE,  Fluxes::N2_FIRERATIO  * nflux_fire);
}

/*
Get the number of established cohorts of this pft in this patch.
*/
int get_nindiv(const Patch &patch, const Pft &pft) {
	int n = 0;
	const Vegetation &vegetation = patch.vegetation;
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		if (indiv.pft.id == pft.id) {
			n++;
		}
	}
	return n;
}

///////////////////////////////////////////////////////////////////////////////////////
// RANDPOISSON
// Internal functions for generating random numbers

int randpoisson(double expectation, long& seed) {

	// DESCRIPTION
	// Returns a random integer drawn from the Poisson distribution with specified
	// expectation (for computational reasons, the Gaussian normal distribution is
	// used as an approximation of the Poisson distribution for expected values >100)

	double p, q, r;
	int n;

	if (expectation <= 100) {

		// For expected values up to 100, calculate a true Poisson number

		p = exp(-expectation);
		q = p;
		r=randfrac(seed);

		n = 0;
		while (q < r) {
			n++;
			p *= expectation / (double)n;
			q += p;
		}
		return n;
	}

	// For higher expected values than 100, approximate the Poisson distribution
	// by the Gaussian normal distribution with mean equal to the expected value,
	// and standard deviation the square root of this value

	do {
		r=randfrac(seed)*8.0-4.0;
		p = exp(-r * r / 2.0);
	} while (randfrac(seed)>p);

	return max(0, (int)(r * sqrt(expectation) + expectation + 0.5));
}


///////////////////////////////////////////////////////////////////////////////////////
// BIOCLIMATIC LIMITS ON ESTABLISHMENT AND SURVIVAL
// Internal functions (do not call directly from framework)

bool establish(Patch& patch, const Climate& climate, Pft& pft) {

	// DESCRIPTION
	// Determines whether specified PFT is within its bioclimatic limits for
	// establishment in a given patch and climate. Returns true if PFT can establish
	// under specified conditions, false otherwise

	// The following limits are implemented:
	//   tcmin_est   = minimum coldest month mean temperature for the last 20 years
	//   tcmax_est   = maximum coldest month mean temperature for the last 20 years
	//   twmin_est   = minimum warmest month mean temperature
	//   gdd5min_est = minimum growing degree day sum on 5 deg C base
	// Special rules for forestry after clearcut
	ManagementType& mt = patch.stand.get_current_management();
	Standpft& spft = patch.stand.pft[pft.id];
	bool pft_selection = spft.plant;

	// Bypass temperature limits of establishment for selected species
	bool relaxed_establishment = patch.plant_this_year && mt.relaxed_establishment && pft_selection;

	if (!relaxed_establishment && (climate.mtemp_min20 < pft.tcmin_est ||
		climate.mtemp_min20 > pft.tcmax_est ||
		climate.mtemp_max < pft.twmin_est ||
		climate.agdd5 < pft.gdd5min_est)) return false;

	if (!relaxed_establishment) {
		// Wolf et al. (2008) bioclimatic limits related to snow depth and GDD0
		if (patch.soil.snow.snowpack_avrmax3mth() < pft.min_snow ||
			climate.agdd0_20.mean() <= pft.gdd0_min || climate.agdd0_20.mean() >= pft.gdd0_max) {
			return false;
		}
	}

	/* Bypass all environmental limits except temperature for selected tree species after clearcut
	 * (par_grass_mean not updated yet)
	 */
	if(patch.plant_this_year && (pft_selection || mt.planting_system == ""))
		return true;

	if(patch.stand.landcover != CROPLAND) {
		if (vegmode != POPULATION && patch.par_grass_mean < pft.parff_min) return false;
	}

	// guess2008 - drought limited establishment
    if (ifdroughtlimitedestab && patch.stand.landcover != PEATLAND) {
		// Compare this PFT's/species' drought_tolerance with the average wcont over the
		// growing season, in this patch. Higher drought_tolerance values (set in the .ins file)
		// lead to greater restrictions on establishment.
		double awcontupper = patch.soil.var_depth(0.0, SOILDEPTH_UPPER, &Soillayer::awcont, &Soillayer::awc);
        if (pft.drought_tolerance > awcontupper) {
           return false;
        }
    }
	// else

	return true;
}


bool survive(const Climate& climate, Pft& pft) {

	// DESCRIPTION
	// Determines whether specified PFT is within its bioclimatic limits for survival
	// in a given climate. Returns true if PFT can survive in the specified climate,
	// false otherwise

	// The following limit is implemented:
	//   tcmin_surv = minimum coldest month temperature for the last 20 years (deg C)

	if (climate.mtemp_min20 < pft.tcmin_surv ||
		climate.mtemp_max20 - climate.mtemp_min20 < pft.twminusc) return false;

	// else

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////
// ESTABLISHMENT
// Internal functions (do not call directly from framework)

void establishment_lpj(Stand& stand,Patch& patch) {

	// DESCRIPTION
	// Establishment in population (standard LPJ) mode.
	// Simulates population increase through establishment each simulation year for
	// trees and grasses and introduces new PFTs if climate conditions become suitable.
	// This function assumes each Individual object represents the average individual
	// of a PFT population, and that there is (at most) one Individual object per PFT
	// per modelled area (stand and patch).
	// This is the only function in which new average individuals are added to the
	// vegetation in population mode.

	const double K_EST=0.12; // maximum overall sapling establishment rate (indiv/m2)

	double fpc_tree;  // summed fractional cover for tree PFTs
	double fpc_grass; // summed fractional cover for grass PFTs
	double est_tree;  // overall establishment rate for trees on modelled area basis (indiv/m2)

	// Establishment rate for a particular PFT on modelled area basis (for trees,
	// indiv/m2; for grasses, fraction of modelled area colonised)
	double est_pft;

	int ntree_est;	  // number of establishing tree PFTs
	int ngrass_est;	  // number of establishing grass PFTs
	bool present;

	// Obtain reference to Vegetation object

	Vegetation& vegetation=patch.vegetation;

	// Loop through all possible PFTs to determine whether there are any not currently
	// present but within their bioclimatic limits for establishment

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft=pftlist.getobj();

		if (stand.pft[pft.id].active) {	// standpft.active is set in landcover_init according to rules for each stand

		// Is this PFT already represented?

			present=false;
			vegetation.firstobj();
			while (vegetation.isobj && !present) {
				Individual& indiv=vegetation.getobj();
				if (indiv.pft.id==pft.id) present=true;
				vegetation.nextobj();
			}

			if (!present) {
				if (establish(patch,stand.get_climate(),pft)) {

					// Not present but can establish, so introduce as new average
					// individual

					Individual& indiv=vegetation.createobj(pft,vegetation);
					if (pft.is_grass() || pft.lifeform == MOSS) {
						indiv.height=0.0;
						indiv.crownarea=1.0; // (value not used)
						indiv.densindiv=1.0;
						indiv.fpc=0.0;
					}
				}
			}
		}
		pftlist.nextobj(); // ... on to next PFT
	}

	// Calculate total FPC and number of PFT's (i.e. average individuals) establishing
	// this year for trees and grasses respectively

	fpc_tree=0.0;
	fpc_grass=0.0;
	ntree_est=0;
	ngrass_est=0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		if (indiv.pft.lifeform==TREE) {
			if (establish(patch, stand.get_climate(), indiv.pft)) ntree_est++;
			fpc_tree+=indiv.fpc;
		}
		else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {
			fpc_grass+=indiv.fpc;
			if (establish(patch, stand.get_climate(), indiv.pft)) ngrass_est++;
		}

		vegetation.nextobj(); // ... on to next individual
	}

	// Calculate overall establishment rate for trees
	// Trees can establish in unoccupied regions of the stand, and in regions
	// occupied by grasses. Establishment reduced in response to canopy closure as
	// tree cover approaches 100%
	// Smith et al 2001, Eqn (17)

	est_tree=K_EST*(1.0-exp(-5.0*(1.0-fpc_tree)))*(1.0-fpc_tree);

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {

		// For this individual ...

		Individual& indiv=vegetation.getobj();

		if (indiv.pft.lifeform==TREE && establish(patch, stand.get_climate(), indiv.pft)) {

			// ESTABLISHMENT OF NEW TREE SAPLINGS
			// Partition overall establishment equally among establishing PFTs

			est_pft=est_tree/(double)ntree_est;

			if (est_pft<0.0)
				fail("establishment_area: negative establishment rate");

			// Adjust density of individuals across modelled area to account for
			// addition of new saplings

			indiv.densindiv+=est_pft;

			// Account for flux from the atmosphere to new saplings
			// (flux is downward and therefore negative)
			if(!indiv.istruecrop_or_intercropgrass()) {
				indiv.report_flux(Fluxes::ESTC,
								  -(indiv.pft.regen.cmass_leaf+
									indiv.pft.regen.cmass_root+indiv.pft.regen.cmass_sap+
									indiv.pft.regen.cmass_heart)*est_pft);
			}

			// Adjust average individual C biomass based on average biomass and density
			// of the new saplings

			indiv.cmass_leaf+=indiv.pft.regen.cmass_leaf*est_pft;
			indiv.cmass_root+=indiv.pft.regen.cmass_root*est_pft;
			indiv.cmass_sap+=indiv.pft.regen.cmass_sap*est_pft;
			indiv.cmass_heart+=indiv.pft.regen.cmass_heart*est_pft;

			// Calculate storage pool size
			indiv.max_n_storage = max(indiv.cmass_sap * indiv.pft.fnstorage / indiv.pft.cton_leaf_avr,
				indiv.cmass_leaf * indiv.pft.turnover_leaf / indiv.pft.cton_leaf_avr + indiv.cmass_root * indiv.pft.turnover_root / indiv.pft.cton_root_avr);
			if (!indiv.alive) {
				indiv.scale_n_storage = indiv.max_n_storage /
					((indiv.pft.regen.cmass_leaf + indiv.pft.regen.cmass_root +
						indiv.pft.regen.cmass_sap + indiv.pft.regen.cmass_heart) * est_pft);
			}
		}
		else if (((indiv.pft.is_grass() || indiv.pft.lifeform==MOSS) && indiv.pft.phenology!=CROPGREEN) &&
			establish(patch, stand.get_climate(), indiv.pft)) {

			// ESTABLISHMENT OF GRASSES AND MOSSES
			// Grasses establish throughout unoccupied regions of the grid cell
			// Overall establishment partitioned equally among establishing PFTs

			est_pft=(1.0-fpc_tree-fpc_grass)/(double)ngrass_est;

			// Account for flux from atmosphere to grass regeneration
			if (!indiv.istruecrop_or_intercropgrass()) {
				indiv.report_flux(Fluxes::ESTC,
			                  -(indiv.pft.regen.cmass_leaf+
			                    indiv.pft.regen.cmass_root)*est_pft);
			}

			// Add regeneration biomass to overall biomass

			indiv.cmass_leaf+=est_pft*indiv.pft.regen.cmass_leaf;
			indiv.cmass_root+=est_pft*indiv.pft.regen.cmass_root;

			// Calculate storage pool size
			indiv.max_n_storage = max(0.0, max(indiv.cmass_root * indiv.pft.fnstorage / indiv.pft.cton_root_avr,
				indiv.cmass_leaf * indiv.pft.turnover_leaf / indiv.pft.cton_leaf_avr + indiv.cmass_root * indiv.pft.turnover_root / indiv.pft.cton_root_avr));
			if (!indiv.alive) {
				indiv.scale_n_storage = indiv.max_n_storage / ((indiv.pft.regen.cmass_leaf + indiv.pft.regen.cmass_root) * est_pft);
			}
		}

		// Update allometry to account for changes in average individual biomass

		allometry(indiv);

		vegetation.nextobj(); // ... on to next individual
	}
}


Individual* get_beneficiary(Patch& patch, Pft& pft) {
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		Individual &indiv = patch.vegetation[i];
		// Check if this individual is of the correct pft.
		if (indiv.pft.id == pft.id) {
			if (pft.lifeform == TREE && indiv.age > 0) {
				// This indiv is the correct PFT, but it's a tree with age > 1 year.
				// This tree will therefore not be the beneficiary of the new biomass.
				continue;
			}
			// Return a pointer to this indiv.
			return &indiv;
		}
	}
	// No indiv found of this pft.
	return NULL;
}


double germination_fraction_grass(Patchpft& ppft, Patch& patch, Climate& climate) {

	const Soil &soil = patch.soil;

	// Water content and temperature of the first soil layer (0-10cm).
	// wcont as fraction of available water holding capacity.
	//double tdorm = ppft.pft.temp_dormancy; 
	double tdorm = calc_temp_dormancy(ppft.pft.temp_dormancy,climate.mtemp_min20);
	double fgerm_max25 = 0.6; // maximum daily germination rate at temperature optimum (fraction) 

	// get average conditions over the last 3-5 (currently 4) days.
	// wcont at 30cm depth.
	double wcont = soil.soillayers[0].wcont();
	ppft.wcont_germ_30.add(wcont);  // longer term wcont used for seed dormancy
	ppft.wcont_germ_4.add(soil.soillayers[0].wcont()); 
	ppft.tsoil_germ_4.add(soil.soillayers[0].temp);
	wcont = ppft.wcont_germ_4.mean();
	double tsoil = ppft.tsoil_germ_4.mean();
	
	// TODO: parameterise this relationship from data
	double fwater = wscal_growth(wcont, wscal0, bgrowth);
	
	double ftemp = temp_growth(tsoil,ppft.pft.Ha_growth,ppft.pft.T0_growth,ppft.pft.A_growth,tdorm);
	double frac_germ = fgerm_max25 * ftemp * fwater; 
	frac_germ = frac_germ > 0.01 ? frac_germ : 0.0;  // avoid very low germination rates

	// Seed dormancy
    // Do not germinate in autumn or winter under temperature-limited conditions
    // In particular, if germination occurs because of increased water availability, do not 
    // check for temperature (or daylength) as in most cases this won't matter.
    // If there is always enough water, determine if it's the right season to germinate
    // by checking daylength and temperature.

    // add up temperature factor after winter solstice
    ppft.ftemp_cum += ftemp;
    if (ppft.daylength_10.periodictrend(1) >= 0.0 && ppft.daylength_10.periodictrend(2) <= 0.0) {
        ppft.ftemp_cum = 0.0;
    }

    // do not germinate if seeds are dormant (see criteria above)
    if (ppft.wcont_germ_30.size() > 10 && ppft.wcont_germ_30.mean() > 0.7) {
        if ((ppft.daylength_10.periodictrend(1) < 0.0 && ppft.reldaylength < 0.9) ||
            (ppft.daylength_10.periodictrend(1) > 0.0 && ppft.ftemp_cum < 10.0)) {
                frac_germ = 0.0;
        }
    }

	return frac_germ;
}

// Calculates the fraction of the seedbank that is germinating today.
double germination_fraction(Patchpft& ppft, Patch& patch, Climate& climate) {

	if (ppft.pft.is_grass()) {
		return germination_fraction_grass(ppft, patch, climate);
	}
	return 0;
}

double nseedlings(double cmass_seed, Pft &pft) {
	// seedling mass (kgC/m^2)
	// TODO: set as PFT parameter
	double cmass_seedling = 0.005;
	
	// grass and moss pft has only 1 individual
	if (pft.is_grass() || pft.lifeform == MOSS) {
		return 1;
	}

	return max(1.0, floor(cmass_seed / cmass_seedling));

}

// call every day
void establishment_daily(Patch& patch, Climate& climate) {
	// check conditions if germination is likely to occur 
	// if yes, invest some of the seedbank into initial biomass
	// if there is no seedbank, maybe don't try every day.
	Vegetation &vegetation = patch.vegetation;
	// initial height of trees (m)
	const double init_height = 0.05; 

	if (ifdailygrass && ifdavegrowth) {
		for (int i = 0; i < pftlist.nobj; i++) {
			Pft &pft = pftlist[i];

			if (!pft.daily_establishment()) {
				// This pft is not using the dave allocation. This could happen,
				// for example, if it's a tree and ifdailytree is false.
				continue;
			}

			Patchpft &ppft = vegetation.patch.pft[pft.id];

			// Check conditions for germination (germination fraction)
			// TODO: Add germination fraction calculation for TREES.
			double frac_germ = germination_fraction(ppft, patch, climate);

			// calculate the C that would germinate
			// Purpose of this is to avoid very low absolute germination rates 
			// (leading to lai < lai_min at establishment)
			// that would 'clutter' the simulations with very small individuals.
			// cmass_leaf = cmass_init * (1 - cstorage_max_temp) * ltor / (ltor + 1)
			// cmass_init = cmass_leaf / ((1 - cstorage_max_temp) * ltor / (ltor + 1))
			double cmass_init_min = get_cmass_init_min(ppft);
			double cmass_init = frac_germ * (ppft.cmass_litter_repr + ppft.cmass_ext_seed);

			// Do conditions allow for germination?
			if (frac_germ > 0.0) {
				Individual *indiv = get_beneficiary(patch, pft);
				if (indiv == NULL) {
					// only establish a new individual if it results in a 
					// sufficiently high lai.
					if (cmass_init < cmass_init_min) {
						continue;
					}
					// number of seedling (ntrees/m^2)
					int nseedlings = 0;
					// in case of TREES, calcualt the number of seedlings
					// by dividing the total cmass_litter_repr by the average
					// mass of a seedling 
					if (pft.lifeform == TREE) {
						// start calculating avrage mass of seedling:
						double height = init_height;
						// calcualte tree diameter (m) using the allometric relationships
						double diam = pow((height / pft.k_allom2), (1 / pft.k_allom3));
						// calculate cmass_sap (kgC/m^2; assuming seedling stems are all sapwood)
						// by multiplying the setm volume and wood density.
						double seedling_cmass_sap = height * pow(diam / 2, 2) * PI * pft.wooddens;
						// use sapwood mass (kgC/m^2) and k_latosa to calculate 
						// cmass_leaf of the new seedling
						double seedling_cmass_leaf = (pft.k_latosa * seedling_cmass_sap) / 
								(pft.sla * height * pft.wooddens);
						// use ltor_max to calculate the cmass_root (kgC/m^2)
						// of the new seedling
						double seedling_cmass_root = seedling_cmass_leaf / pft.ltor_max;
						// total seedling mass (kgC/m^2)
						double seedling_cmass = seedling_cmass_sap + seedling_cmass_root 
								+ seedling_cmass_leaf;
						// number of seeds = int(total seedling mass / mass of one seedling).
						// TODO: make new function to calculate frac_germ that 
						// will fit TREES.
						nseedlings = (frac_germ * ppft.cmass_litter_repr) / seedling_cmass;
						// check if there is not enough mass for establishing even one tree
						if (nseedlings < 1) {
							// go to next pft
							continue;
						}
					}

					// Establish new individual.
					log_info("Establishing %s...", (const char *)pft.name);
					Individual& individual = vegetation.createobj(pft, vegetation);
					indiv = &individual; // & gets address of the variable (* is the opposite)
					//indiv->height=0.0;  // -> access properties of a pointer (pointer is the address of the object)
					indiv->crownarea=1.0; // (value not used). this is the same as (*indiv).crownarea
					indiv->densindiv=1.0;
					if (pft.lifeform == TREE) {
						// in case of trees, set the new individual's height, 
						// diameter, canopy area and seeding density (calculated
						// above)
						indiv->height = init_height;
						indiv->diam = pow((indiv->height / pft.k_allom2), (1 / pft.k_allom3));
						indiv->crownarea = pft.k_allom1 * pow(indiv->diam, pft.k_rp); // (value not used). this is the same as (*indiv).crownarea
						indiv->densindiv = nseedlings;
						indiv->age = 1.0 / 365.0;
					}

					indiv->fpc=1.0;  // TODO: is that right?
					indiv->age=0.0;
					indiv->age_active=0.0;
					indiv->nscal_5.add(1.0);
					indiv->growth_stage = vegetative;
					indiv->frac_germ = frac_germ;
				} else {
					// If there is already an individual, then:
					// 1) reduce germination rate by fpc in the vegetative phase. I.e. do not grow new plants in space
					//    that is already occupied.
					// 2) avoid germination when main individual is in the reproductive phase,
					// 	  assuming that seeds are dormant in that time period, avoiding unfavourable growth conditions.
					//    E.g. no germination in a wet and warm autumn.
					//    TODO: implement seed dormancy to account for this when there is no existing individual?
					if (indiv->growth_stage == vegetative) {
						// reduced by fpc_ymax * 2.0 to account for green and brown leaves and rpc.
						// TODO: account for fpc_ymax of all individuals (i.e. all pfts).
						//       and create new variable.
						indiv->frac_germ = frac_germ * (1.0 - indiv->fpc_ymax * 2.0);
					} 
					if (indiv->growth_stage == reproductive || is_zero(indiv->fdev_growth - fdev_growth_min) || indiv->growth_stage == dormant) {
						indiv->frac_germ = 0.0;
					}
				}

				if (indiv->frac_germ > 0.0) {
					double cmass_ext_seed = ppft.cmass_ext_seed * indiv->frac_germ;
					double nmass_ext_seed = ppft.nmass_ext_seed * indiv->frac_germ;

					double cmass_seed = ppft.cmass_litter_repr * indiv->frac_germ;
					double nmass_seed = ppft.nmass_litter_repr * indiv->frac_germ;

					double cmass_total_seed = cmass_seed + cmass_ext_seed;
					double nmass_total_seed = nmass_seed + nmass_ext_seed;

					if (is_zero(cmass_total_seed) || is_zero(nmass_total_seed)) {
						// no seeds to germinate
						continue;
					}

					double cmass_before = ppft.cmass_litter_repr + indiv->ccont() + ppft.cmass_ext_seed;
					double nmass_before = ppft.nmass_litter_repr + indiv->ncont() + ppft.nmass_ext_seed;
			
					// Take germinated biomass out of the seedbank
					ppft.cmass_litter_repr -= cmass_seed;
					ppft.nmass_litter_repr -= nmass_seed;

					ppft.cmass_ext_seed -= cmass_ext_seed;
					ppft.nmass_ext_seed -= nmass_ext_seed;

					// With a germination event, the hypothetical seed bank comes into existence
					patch.fluxes.report_flux(Fluxes::ESTC, -cmass_ext_seed);
					patch.fluxes.report_flux(Fluxes::ESTN, -nmass_ext_seed);
				
					// Add some new biomass into this individual.
					establishment_biomass(cmass_total_seed, pft.ltor_max, *indiv, nmass_total_seed); // * is a dereference (gets the value at the address)

					indiv->cmass_germ = cmass_ext_seed;
					indiv->acmass_germ += cmass_ext_seed;

					double cmass_after = ppft.cmass_litter_repr + indiv->ccont(1.0, false, true) + ppft.cmass_ext_seed;
					double nmass_after = ppft.nmass_litter_repr + indiv->ncont(1.0, false, true) + ppft.nmass_ext_seed;
					assert(equal(cmass_before, cmass_after));
					assert(equal(nmass_before, nmass_after));
				}
			}
		}
	}
}

double cmass_repr_tot(const Stand &stand, const Pft &pft) {
	double total = 0;
	for (size_t i = 0; i < stand.nobj; i++) {
		const Vegetation &vegetation = stand[i].vegetation;
		for (size_t j = 0; j < vegetation.nobj; j++) {
			const Individual &indiv = vegetation[j];
			if (indiv.pft.id == pft.id) {
				total += indiv.cmass_repr;
			}
		}
	}
	return total / (double)stand.npatch();
}

void establishment_guess(Stand& stand,Patch& patch) {

	// DESCRIPTION
	// Establishment in cohort or individual mode.
	// Establishes new tree saplings and simulates grass population increase each
	// establishment year (every year in individual mode). New grass PFTs are
	// introduced (any year) if climate conditions become suitable.
	// Saplings are given an amount of carbon proportional to the hypothetical forest
	// floor NPP for this year and are "moulded" by calls to the standard allocation
	// function (growth module). The coefficient SAPSIZE controls the relative amount
	// of carbon given each sapling and should be set to a value resulting in saplings
	// around breast height in the first year (before development of a shading canopy).

	// For tree PFTs within their bioclimatic limits for establishment, the expected
	// number of saplings (est) established in each patch is given by:
	//
	// for model initialisation:
	//   (1) est = est_max*patcharea
	// if spatial mass effect enabled (stand-level "propagule pool" influences
	// establishment):
	//   (2) est = c*(kest_repr*cmass_repr+kest_bg)
	// if no spatial mass effect and propagule pool exists (cmass_repr non-negligible):
	//   (3) est = c*(kest_pres+kest_bg)
	// if no spatial mass effect and no propagule pool:
	//   (4) est = c*kest_bg
	// where
	//   (5) c = mu(anetps_ff/anetps_ff_max)*est_max*patcharea
	//   (6) mu(F) = exp(alphar*(1-1/F)) (Fulton 1991, Eqn 4)
	// where
	//   kest_repr     = PFT-specific constant (kest_repr*cmass_repr should
	//                   approximately equal 1 at the highest plausible value of
	//                   cmass_repr for the PFT)
	//   kest_bg       = PFT-specific constant in range 0-1 (usually << 1)
	//                   (set to zero if background establishment disabled)
	//   kest_pres     = PFT-specific constant in range 0-1 (usually ~1)
	//   cmass_repr    = net C allocated to reproduction for this PFT in all patches of
	//                   this stand this year (kgC/m2)
	//   est_max       = (nominal) maximum establishment rate for this PFT
	//                  (saplings/m2/year)
	//   patcharea     = patch area (m2)
	//   anetps_ff     = potential annual net assimilation at forest floor for this
	//                   PFT (kgC/m2/year)
	//   anetps_ff_max = maximum value of anetps_ff for this PFT in this stand so far
	//                   in the simulation (kgC/m2/year)
	//   alphar        = PFT-specific constant (parameter capturing non-linearity in
	//                   recruitment rate relative to understorey growing conditions)
	//
	// In individual mode, and in cohort mode if stochastic establishment is enabled,
	// the actual number of new saplings established is a number drawn from the Poisson
	// distribution with expectation 'est'. In cohort mode, with stochastic
	// establishment disabled, a cohort representing exactly 'est' individuals (may be
	// not-integral) is established.

	// Coefficient in calculation of initial sapling size and initial
	// grass biomass (see comment above)
	const double SAPSIZE = 0.1;

	// Minimum carbon mass allowed (kgC/m2)
	const double MINCMASS = 1.0e-8;
	bool present = false; // whether PFT already present in this patch
	double c = 0.0; // constant in equation for number of new saplings (Eqn 5)
	double est = 0.0; // expected number of new saplings for PFT in this patch

	// Actual number of new saplings for PFT in this patch (may include a
	// fractional part in cohort mode with stochastic establishment disabled)
	double nsapling = 0.0;
	double bminit = 0.0; // initial sapling biomass (kgC) or new grass biomass (kgC/m2)
	double ltor = 0.0; // leaf to fine root mass ratio for new saplings or grass
	int newindiv = 0; // number of new Individual objects to add to vegetation for this PFT
	double kest_bg = 0.0;
	int i = 0;

	// Obtain reference to Vegetation object

	Vegetation& vegetation=patch.vegetation;

	const bool establish_active_pfts_before_management = true;
	const bool cloned_or_changed_man = stand.cloned || stand.current_rot;

	ManagementType& mt = patch.stand.get_current_management();
	if(patch.plant_this_year && mt.set_planting_density && mt.planting_system != "") {
		planting(patch);	// NB planting() gives small fixed sapsize
		return;
	}

	// guess2008 - determine the number of woody PFTs that can establish
	// Thomas Hickler
	int nwoodypfts_estab=0;
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft=pftlist.getobj();
		Standpft& standpft=stand.pft[pft.id];

		// Force planting if planting() not used
		bool force_planting = 
			patch.plant_this_year && standpft.plant && !(mt.set_planting_density && mt.planting_system != "");
		bool est_this_year;
		if(establish_active_pfts_before_management) {
			est_this_year = !patch.management_started && !cloned_or_changed_man || !patch.plant_this_year && standpft.reestab;
		}
		else {
			est_this_year = !run_landcover || !patch.plant_this_year && standpft.reestab;
		}

		if (establish(patch, stand.get_climate(), pft) && pft.lifeform == TREE && standpft.active 
				&& (est_this_year || force_planting)) {
			nwoodypfts_estab++;
		}
		pftlist.nextobj();
	}


	// Loop through PFTs
	
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft=pftlist.getobj();

		// define sapling size as teh constant default.
		double sapsize = SAPSIZE;

		if (ifdavegrowth){
			// change sapling size to what is defined in the PFT (only with dave growth).
			sapsize = pft.sapsize;
		}
		// from this place forward, use onlt sapsize instaed of SAPSIZE

		if (pft.dave_allocation() && pft.lifeform != TREE) {
			// This pft is using dave daily allocation and establishment.
			pftlist.nextobj();
			continue;
		}

		Standpft& standpft=stand.pft[pft.id];
		Patchpft& ppft = patch.pft[pft.id];

		// For this PFT ...

		// Stands cloned this year to be treated here as first year
		bool init_clone = date.year == stand.clone_year && pft.landcover == stand.landcover;

		// Force planting of pfts in selection in managed forest stands after clear-cut if planting() not used
		// No grass establishment will occur during planting year (force_planting is always false for grass and 
		// est_this_year is false when patch.plant_this_year is true)
		bool force_planting = 
			patch.plant_this_year && standpft.plant && !(mt.set_planting_density && mt.planting_system != "");
		bool est_this_year;
		if(establish_active_pfts_before_management) {
			est_this_year = !patch.management_started && !cloned_or_changed_man || !patch.plant_this_year && standpft.reestab;
		}
		else {
			est_this_year = !run_landcover || !patch.plant_this_year && standpft.reestab;
		}

		if (stand.pft[pft.id].active) {
			if (patch.age==0 || init_clone) {
				ppft.anetps_ff_est=ppft.anetps_ff;
				ppft.wscal_mean_est=ppft.wscal_mean;

				// BLARP
				if (date.year==0 || date.year==stand.first_year || init_clone || force_planting) // also after clearcut
					ppft.anetps_ff_est_initial=ppft.anetps_ff;

			}
			else {
				ppft.anetps_ff_est+=ppft.anetps_ff;
				ppft.wscal_mean_est+=ppft.wscal_mean;
			}

			if (establish(patch, stand.get_climate(), pft) && (est_this_year || force_planting)) {

				if (pft.is_grass() || pft.lifeform==MOSS) {

					// ESTABLISHMENT OF GRASSES AND MOSSES

					// Each grass PFT represented by just one individual in each patch
					// Check whether this grass PFT already represented ...

					present=false;
					vegetation.firstobj();
					while (vegetation.isobj && !present) {
						if (vegetation.getobj().pft.id==pft.id) present=true;
						vegetation.nextobj();
					}

					if (!present) {

						// ... if not, add it
						log_info("Establishing %s...", (const char *)pft.name);
						Individual& indiv=vegetation.createobj(pft,vegetation);
						indiv.height=0.0;
						indiv.crownarea=1.0; // (value not used)
						indiv.densindiv=1.0;
						indiv.fpc=1.0;
						indiv.nscal_5.add(1.0);

						// Initial grass biomass proportional to potential forest floor
						// net assimilation this year on patch area basis

						if(pft.phenology == CROPGREEN) {
							bminit = sapsize * 0.01;
						}
						else if(patch.disturbed) {	// Removing has_disturbances() to allow this for all stands disturbed in disturbe_force()
							bminit = sapsize * patch.pft[pft.id].anetps_ff_est_initial;
						} else
							bminit = sapsize * patch.pft[pft.id].anetps_ff;

						if (ifdailygrass && pft.phenology != CROPGREEN)
							bminit *= 0.1;

						// Initial leaf to fine root biomass ratio based on
						// hypothetical value of water stress parameter

						ltor=ppft.wscal_mean*pft.ltor_max;

						// Allocate initial biomass
						if(!indiv.istruecrop_or_intercropgrass())
							allocation_init(bminit,ltor,indiv);

						// Calculate initial allometry

						allometry(indiv);

						// Calculate storage pool size
						if(indiv.pft.phenology == CROPGREEN) {
							indiv.max_n_storage = CMASS_SEED / indiv.pft.cton_leaf_avr;
							indiv.scale_n_storage = indiv.max_n_storage / CMASS_SEED;
						}
						else {
							indiv.max_n_storage = max(0.0, max(indiv.cmass_root * indiv.pft.fnstorage * indiv.pft.cton_leaf_avr,
								indiv.cmass_leaf * indiv.pft.turnover_leaf / indiv.pft.cton_leaf_avr + indiv.cmass_root * indiv.pft.turnover_root / indiv.pft.cton_root_avr));
							indiv.scale_n_storage = indiv.max_n_storage / bminit;
						}

						// Establishment flux is not debited for 'new' Individual
						// objects - their carbon is debited in function growth()
						// if they survive the first year

						// If established individuals initial biomass is too small 
						// then don't establish (e.g. anetps_ff extremely low or even zero)
						if (bminit < MINCMASS) {
							indiv.kill(MORT_OTHER);

							vegetation.killobj();
						}
					}
				}
				else if (pft.lifeform==TREE) {

					// ESTABLISHMENT OF NEW TREE SAPLINGS

					if (patch.age == 0 || init_clone){

						// First simulation year - initialising patch
						// Eqn 1

						est=pft.est_max*patcharea;
					}
					else {

						// Every year except year 1
						// Eqns 5, 6

						if (ppft.anetps_ff>0.0 &&
							!negligible(ppft.anetps_ff)) {

							c=exp(pft.alphar-pft.alphar/ppft.anetps_ff*
								stand.pft[pft.id].anetps_ff_max)*pft.est_max*patcharea;
						}
						else {
							c=0.0;
						}

						// Background establishment enabled?

						if (ifbgestab)
							kest_bg=pft.kest_bg;
						else
							kest_bg=0.0;

						// Spatial mass effect enabled?
						// Eqns 2, 3, 4
						double cmass_repr;
						if (est_live_repr)
							cmass_repr = cmass_repr_tot(stand, pft);
						else
							cmass_repr = stand.pft[pft.id].cmass_repr;

						if (ifsme)
							est=c*(pft.kest_repr*cmass_repr+kest_bg);
						else if (!negligible(stand.pft[pft.id].cmass_repr))
							est=c*(pft.kest_pres+kest_bg);
						else
							est=c*kest_bg;
					}


					// guess2008 - scale est by the number of woody PFTs/species that can establish
					// Otherwise, simply adding more PFTs or species would increase est
					est*=min(1.0,3.0/double(nwoodypfts_estab));
					//est*=3.0/double(nwoodypfts_estab);


					// Have a value for expected number of new saplings (est)
					// Actual number of new saplings drawn from the Poisson distribution
					// (except cohort mode with stochastic establishment disabled)

					if (ifstochestab && mt.stochestab && !force_planting || vegmode==INDIVIDUAL) {
						nsapling=randpoisson(est, stand.seed);
					}
					else { 
						nsapling=est;
					}

					if (vegmode==COHORT) {

						patch.pft[pft.id].nsapling_year = nsapling;

						// BLARP added for OECD experiment (is this sensible?)
						if (patch.disturbed || force_planting) {	// Removing has_disturbances() to allow this for all stands disturbed in disturbe_force()

							ppft.anetps_ff_est=
								ppft.anetps_ff_est_initial;
							ppft.wscal_mean_est=ppft.wscal_mean;
							newindiv=!negligible(nsapling / patcharea);
						}
						else if (patch.age%estinterval && !init_clone || patch.plant_this_year) {

							// Not an establishment year - save sapling count for
							// establishment the next establishment year
							// Establishment each year in managed patches

							ppft.nsapling+=nsapling;
							newindiv=0;
						}
						else {
							if (patch.age && !init_clone) {// all except first year after disturbance
								nsapling+=ppft.nsapling;
								ppft.anetps_ff_est/=(double)estinterval;
								ppft.wscal_mean_est/=(double)estinterval;
							}
							newindiv=!negligible(nsapling / patcharea);
								// round down to 0 if nsapling very small
							if (pft.dave_allocation() && nsapling < 1) {
								// Only establish with at least 1 sapling.
								newindiv = 0;
							}
						}
					}
					else if (vegmode == INDIVIDUAL) {
						newindiv=(int)(nsapling+0.5); // round up to be on the safe side
					}

					if (pft.lifeform == TREE && nindiv_max > 0 && date.year > freenyears) {
						int nindiv = get_nindiv(patch, pft);
						if (nindiv >= nindiv_max) {
							if (newindiv > 0) {
								// Reset running totals.
								patch.pft[pft.id].anetps_ff_est = 0;
								patch.pft[pft.id].wscal_mean_est = 0;
								patch.pft[pft.id].nsapling = 0;
							}
							pftlist.nextobj();
							continue;
						}
					}

					// Now create 'newindiv' new Individual objects

					if (newindiv > 0) {
						log_info("Establishing new %s cohort with %.0f tree%s...", (const char *)pft.name, nsapling, nsapling > 1 ? "s" : "");
					}

					// Initial biomass proportional to potential forest floor net
					// assimilation for this PFT in this patch

					if (patch.has_been_cut && (SMALL_SAPSIZE_POST_CUT && force_planting && mt.planting_system != "")) {
						bminit = PLANTSIZE; // Fixed sap size after first cutting
					}
					else {
						bminit = sapsize*patch.pft[pft.id].anetps_ff_est;
					}

					// Initial leaf to fine root biomass ratio based on hypothetical
					// value of water stress parameter
					ltor=patch.pft[pft.id].wscal_mean_est*pft.ltor_max;

					// If l:r ratio is zero, solve_csap() will fail to find a
					// solution for initial biomass allocation. In this case,
					// we just skip establishment. If wscal and ltor are really
					// zero, any established saplings would quickly die anyway.
					if (bminit == 0 || fabs(ltor) < 1e-4) {
						log_diagnostic("%s: Establishment skipped; ltor=%f, wscal_mean_est=%f, bminit=%f",
									   (const char *)pft.name,
									   ltor,
									   patch.pft[pft.id].wscal_mean_est,
									   bminit);
						patch.pft[pft.id].anetps_ff_est = 0;
						patch.pft[pft.id].wscal_mean_est = 0;
						patch.pft[pft.id].nsapling = 0;
						pftlist.nextobj();
						continue;
					}

					for (i=0;i<newindiv;i++) {

						// Create average individual for a new cohort (cohort mode)
						// or an actual individual (individual mode)

						Individual& indiv=vegetation.createobj(pft,vegetation);

						if (vegmode==COHORT)
							indiv.densindiv=nsapling/patcharea;
						else if (vegmode==INDIVIDUAL)
							indiv.densindiv=1.0/patcharea;

						indiv.age=0.0;

						// Initial biomass proportional to potential forest floor net
						// assimilation for this PFT in this patch

						// Height will be set in allometry very soon, but it's
						// accessed in allocation_init(), so let's initialize it
						// to 0 here rather than rely on uninitialized memory.
						indiv.height = 0;

						// set the initial ltor to the indiv object
						indiv.ltor = ltor;

						if (ifnewalloc) {
							// A different made-up number. Not sure why.
							bminit = 1.0 / indiv.pft.sla;
						}

						// Allocate initial biomass
						allocation_init(bminit,ltor,indiv);

						// Calculate initial allometry
						if (pft.dave_allocation()) {
							allometry_dave_tree(indiv);
						} else {
							allometry(indiv);
						}

						// Calculate storage pool size
						// These are not used by dave cohorts, and some of the
						// parameters used to calculate them are unavailable.
						// Therefore it's best to not even attempt to calculate
						// them for dave cohorts.
						if (!pft.dave_allocation()) {
							indiv.max_n_storage = max(indiv.cmass_sap * indiv.pft.fnstorage / indiv.pft.cton_leaf_avr,
								indiv.cmass_leaf * indiv.pft.turnover_leaf / indiv.pft.cton_leaf_avr + indiv.cmass_root * indiv.pft.turnover_root / indiv.pft.cton_root_avr);
							indiv.scale_n_storage = indiv.max_n_storage / bminit;
						}

						// Establishment flux is not debited for 'new' Individual
						// objects - their carbon is debited in function growth()
						// if they survive the first year

						// If established individuals initial biomass is too small 
						// then don't establish (e.g. anetps_ff extremely low or even zero)
						if (bminit < MINCMASS) {
							indiv.kill(MORT_OTHER);

							vegetation.killobj();
						}
					}
				}
			}

			// Reset running sums for next year (establishment years only in cohort mode)

			if (vegmode!=COHORT || !(patch.age%estinterval)) {
				ppft.nsapling=0.0;
				ppft.wscal_mean_est=0.0;
				ppft.anetps_ff_est=0.0;
			}

		}

		// ... on to next PFT
		pftlist.nextobj();
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// MORTALITY
// Internal functions (do not call directly from framework)

void mortality_lpj(Stand& stand, Patch& patch, const Climate& climate, double fireprob) {

	// DESCRIPTION
	// Mortality in population (standard LPJ) mode.
	// Simulates population decrease through mortality each simulation year for trees
	// and grasses; "kills" PFT populations if climate conditions become unsuitable for
	// survival.
	// This function assumes each Individual object represents the average individual
	// of a PFT population, and that there is (at most) one Individual object per PFT
	// per modelled area (stand and patch).

	// For tree PFTs, the fraction of the population killed is given by:
	//
	//   (1) mort = min ( mort_greff + mort_shade + mort_fire, 1)
	// where mort_greff = mortality preempted by low growth efficiency, given by:
	//   (2) mort_greff = K_MORT1 / (1 + K_MORT2*greff)
	//   (3) greff = annual_npp / leaf_area
	//   (4) leaf_area = cmass_leaf * sla
	// mort_shade = mortality due to shading ("self thinning") as total tree cover
	//   approaches 1 (see code)
	// mort_fire = mortality due to fire when using fire model GLOBFIRM. 
	//   Fire model BLAZE handles fire mortality separately. See blaze.cpp.
	//   the fraction of the modelled area affected by
	//   fire (fireprob) is calculated in function fire; actual mortality is influenced
	//   by PFT-specific resistance to burning (see code).
	//
	// Grasses are subject to shading and fire mortality only

	// References: Sitch et al 2001, Smith et al 2001, Thonicke et al 2001

	// INPUT PARAMETER
	// fireprob = fraction of modelled area affected by fire

	// Constant in mortality equation [c.f. mort_max; LPJF]
	const double K_MORT1=0.01; 
	// Constant in mortality equation [c.f. k_mort; LPJF]; the value here differs
	// from LPJF's k_mort in that growth efficiency here based on annual NPP, c.f.
	// net growth increment; LPJF	
	const double K_MORT2=35.0;
	const double FPC_TREE_MAX=0.95; // maximum summed tree fractional projective cover

	double fpc_dec; // required reduction in FPC
	double fpc_tree; // summed FPC for all tree PFTs
	double fpc_grass; // summed FPC for all grasses
	double deltafpc_tree_total; // total tree increase in FPC this year
	
	// Total mortality for PFT (fraction of current FPC)
	double mort;
	// Background mortality plus mortality preempted by low growth efficiency
	// (fraction of current FPC)
	double mort_greffic;
	// Mortality associated with light competition (fraction of current FPC)
	double mort_shade;
	double mort_fire;
	bool killed;

	// Obtain reference Vegetation object

	Vegetation& vegetation=patch.vegetation;

	// Calculate total tree and grass FPC and total increase in FPC this year for trees

	fpc_tree=0.0;
	fpc_grass=0.0;
	deltafpc_tree_total=0.0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		if (indiv.pft.lifeform==TREE) {

			fpc_tree+=indiv.fpc;
			if (indiv.deltafpc<0.0) indiv.deltafpc=0.0;
			deltafpc_tree_total+=indiv.deltafpc;

		}
		else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) fpc_grass+=indiv.fpc;

		vegetation.nextobj(); // ... on to next individual
	}

	// KILL PFTs BEYOND BIOCLIMATIC LIMITS FOR SURVIVAL

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		killed=false;

		if (!survive(climate,indiv.pft)) {

			// Remove completely if PFT beyond its bioclimatic limits for survival

			Patchpft& patchpft = patch.pft[indiv.pft.id];
			patchpft.cmass_mort += indiv.ccont();
			// patchpft.cmass_mort_wood += indiv.cmass_wood();
			// if (indiv.diam >= 0.1)
			// 	patchpft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

			indiv.kill(MORT_OTHER);

			vegetation.killobj();
			killed=true;
		}

		if (!killed) vegetation.nextobj(); // ... on to next individual
	}

	// CALCULATE AND IMPLEMENT MORTALITY

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		killed=false;

		if (indiv.pft.lifeform==TREE) {

			// TREE MORTALITY

			// Mortality associated with low growth efficiency (includes
			// parameterisation of population background mortality)
			// NOTE: growth efficiency quantified as (annual NPP)/(leaf area)
			//       c.f. LPJF, which uses net growth increment instead of NPP

			mort_greffic=K_MORT1/(1.0+K_MORT2*max(indiv.anpp,0.0)/indiv.cmass_leaf/
				indiv.sla());

			// Mortality associated with light competition
			// Self thinning imposed when total tree cover above FPC_TREE_MAX,
			// partitioned among tree PFTs in proportion to this year's FPC increment

			if (fpc_tree>FPC_TREE_MAX) {
				if (!negligible(deltafpc_tree_total))
					fpc_dec=(fpc_tree-FPC_TREE_MAX)*indiv.deltafpc/deltafpc_tree_total;
				else
					fpc_dec=0.0;
				mort_shade=1.0-fracmass_lpj(indiv.fpc-fpc_dec,indiv.fpc,indiv);
			}
			else
				mort_shade=0.0;

			// Mortality due to fire when using fire model GLOBFIRM (BLAZE handles it separately)
			if (patch.has_fires() && firemodel == GLOBFIRM) {
				mort_fire=fireprob*(1.0-indiv.pft.fireresist);
			}
			else { 
				mort_fire=0.0; 
			}

			// Sum mortality components to give total mortality (maximum 1)

			mort=min(mort_greffic+mort_shade+mort_fire,1.0);

			// Reduce population density and C biomass on modelled area basis
			// to account for loss of killed individuals
			double closs, nloss;
			indiv.reduce_biomass(mort, mort_fire, mort_greffic, mort_shade, 0.0, 0.0, false, closs, nloss);
			if (mort > 0) {
				indiv.closs_globfirm += closs * mort_fire / mort;
				indiv.nloss_globfirm += nloss * mort_fire / mort;
				indiv.closs_mort += closs * (1.0 - mort_fire / mort);
				indiv.nloss_mort += nloss * (1.0 - mort_fire / mort);
			}
		}
		else if (indiv.pft.is_grass() || indiv.pft.lifeform == MOSS) {

			// GRASS AND MOSS MORTALITY

			if (indiv.pft.landcover==CROPLAND)
				fpc_grass=indiv.fpc;

			// Shading mortality: grasses can persist only on regions not occupied
			// by trees

			if (fpc_grass>1.0-min(fpc_tree,FPC_TREE_MAX)) {
				fpc_dec=(fpc_grass-1.0+min(fpc_tree,FPC_TREE_MAX))*indiv.fpc/fpc_grass;
				mort_shade=1.0-fracmass_lpj(indiv.fpc-fpc_dec,indiv.fpc,indiv);
			}
			else
				mort_shade=0.0;

			// Mortality due to fire

			if (patch.has_fires() && firemodel == GLOBFIRM) {
				mort_fire=fireprob*(1.0-indiv.pft.fireresist);
			}
			else {
				mort_fire=0.0;
			}

			// Sum mortality components to give total mortality (maximum 1)

			mort=min(mort_shade+mort_fire,1.0);

			// Reduce C biomass on modelled area basis to account for biomass lost
			// through mortality
			double closs, nloss;
			indiv.reduce_biomass(mort, mort_fire, 0.0, mort_shade, 0.0, 0.0, false, closs, nloss);
			if (mort > 0) {
				indiv.closs_globfirm += closs * mort_fire / mort;
				indiv.nloss_globfirm += nloss * mort_fire / mort;
				indiv.closs_mort += closs * (1.0 - mort_fire / mort);
				indiv.nloss_mort += nloss * (1.0 - mort_fire / mort);
			}
		}

		// Remove this PFT population completely if all individuals killed

		if (negligible(indiv.densindiv)) {
			vegetation.killobj();
			killed=true;
		}
		else {
			// Update allometry
			if (!allometry(indiv)) {
				// Handle allometry related mortality
				Patchpft& patchpft = patch.pft[indiv.pft.id];
				patchpft.cmass_mort += indiv.ccont();
				// patchpft.cmass_mort_wood += indiv.cmass_wood();
				// if (indiv.diam >= 0.1)
				// 	patchpft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

				indiv.kill(MORT_OTHER);

				vegetation.killobj();
				killed = true;
			}
		}

		if (!killed) vegetation.nextobj(); // ... on to next individual
	}
}

void mortality_guess(Stand& stand, Patch& patch, const Climate& climate, double fireprob) {

	// DESCRIPTION
	// Mortality in cohort and individual modes.
	// Simulates death of individuals or reduction in individual density through
	// mortality (including fire mortality) each simulation year for trees, and
	// reduction in biomass due to fire for grasses; kills individuals/cohorts if
	// climate conditions become unsuitable for survival.
	//
	// For trees, death can result from the following factors:
	//
	//   (1) a background mortality rate related to PFT mean longevity;
	//   (2) when mean growth efficiency over an integration period (five years) falls
	//       below a PFT-specific threshold;
	//   (3) stress associated with high temperatures (a "soft" bioclimatic limit,
	//       intended mainly for boreal tree PFTs; Sitch et al 2001, Eqn 55);
	//   (4) fire (probability calculated in function fire)
	//   (5) when climatic conditions are such that all individuals of the PFT are
	//       killed.
	//
	// For grasses, (4) and (5) above are modelled only.
	//
	// For trees and in cohort mode, mortality may be stochastic or deterministic
	// (according to the value of global variable ifstochmort). In stochastic mode,
	// expected mortality rates are imposed as stochastic probabilities of death; in
	// deterministic mode, cohort density is reduced by the fraction represented by
	// the mortality rate.

	// INPUT PARAMETER
	// fireprob = probability of fire in this patch

	// Overall mortality (excluding fire mortality): fraction of cohort killed, or:
	// probability of individual being killed
	double mort;
	// Background component of overall mortality (see 'mort')
	double mort_min;
	// Component of overall mortality associated with low growth efficiency
	double mort_greff;
	// Expected fraction of cohort killed (or: probability of individual being
	// killed) due to fire
	double mort_fire;
	// Fraction of cohort (or individual) surviving
	double frac_survive;
	// Growth efficiency for individual/cohort this year (kgC/m2 leaf/year)
	double greff;
	// Five-year-mean growth efficiency (kgC/m2 leaf/year)
	double greff_mean;
	// Number of individuals (remaining) in cohort
	int nindiv; 
	// Number of individuals in cohort prior to mortality
	int nindiv_prev; 
	int i;
	bool killed;

	// Value of mort_greff when growth efficiency below PFT-specific threshold
	const double KMORTGREFF=0.3;
	// Coefficient in calculation of background mortality (negated natural log of
	// fraction of population expected to survive to age 'longevity'; see Eqn 14
	// below)
	const double KMORTBG_LNF=-log(0.001);
	// Exponent in calculation of background mortality (shape parameter for
	// relationship between mortality and age (0=constant mortality; 1=linear
	// increase; >1->exponential; steepness increases for increasing positive
	// values)
	const double KMORTBG_Q=2.0;

	// Obtain reference to Vegetation object for this patch

	Vegetation& vegetation=patch.vegetation;

	ManagementType& mt = stand.get_current_management();

	// FPC on peatlands
	double fpc_grass = 0.0;
	double fpc_moss = 0.0;
	double fpc_other = 0.0;

	if (patch.stand.landcover == PEATLAND && wetlandlailimit > 0) {

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();

			// For this individual ...
			if (indiv.pft.lifeform == MOSS)
				fpc_moss += indiv.fpc;
			else if (indiv.pft.is_grass())
				fpc_grass += indiv.fpc;
			else
				fpc_other += indiv.fpc;

			vegetation.nextobj(); // ... on to next individual
		}
	}


	// FIRE MORTALITY GLOBFIRM
	// For BLAZE there is a separate call in simulate_day (framework.cpp)
	if (patch.has_fires() && firemodel == GLOBFIRM) {

		// Impose fire in this patch with probability 'fireprob'

		if (randfrac(stand.seed)<fireprob) {

			// Loop through individuals

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();

				// For this individual ...

				killed=false;

				// Fraction of biomass destroyed by fire

				mort_fire=1.0-indiv.pft.fireresist;

				if (indiv.pft.is_grass() || indiv.pft.lifeform==MOSS) {

					// GRASS AND MOSS PFTs

					// Reduce individual biomass
					indiv.reduce_biomass(mort_fire, mort_fire, 0.0, 0.0, 0.0,
										 0.0, false,
						                 indiv.closs_globfirm,
										 indiv.nloss_globfirm);//MORT

					// Update allometry

					allometry(indiv);
				}
				else {

					// TREE PFT

					if (ifstochmort && mt.stochmort) {

						// Impose stochastic mortality
						// Each individual in cohort dies with probability 'mort_fire'

						// Number of individuals represented by 'indiv'
						// (round up to be on the safe side)

						nindiv=(int)(indiv.densindiv*patcharea+0.5);
						nindiv_prev=nindiv;

						for (i=0;i<nindiv_prev;i++)
							if (randfrac(stand.seed)>indiv.pft.fireresist) nindiv--;

						if (nindiv_prev)
							frac_survive=(double)nindiv/(double)nindiv_prev;
						else
							frac_survive=0.0;
					}

					// Deterministic mortality (cohort mode only)

					else frac_survive=1.0-mort_fire;

					// Reduce individual biomass on patch area basis
					// to account for loss of killed individuals
					indiv.reduce_biomass(1.0 - frac_survive, 1.0 - frac_survive,
										 0.0, 0.0, 0.0, 0.0, false,
										 indiv.closs_globfirm,
										 indiv.nloss_globfirm);//MORT

					// Remove this cohort completely if all individuals killed
					// (in individual mode: removes individual if killed)

					if (negligible(indiv.densindiv)) {
						log_info("%s killed by fire...", (const char *)indiv.pft.name);
						vegetation.killobj();
						killed=true;
					}
				}

				if (!killed) vegetation.nextobj(); // ... on to next individual
			}
		}
	}

	// MORTALITY NOT ASSOCIATED WITH FIRE

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
		Patchpft& ppft = patch.pft[indiv.pft.id];

		// For this individual ...

		killed=false;

		// KILL INDIVIDUALS/COHORTS BEYOND BIOCLIMATIC LIMITS FOR SURVIVAL

		if (!survive(climate,indiv.pft)) {

			// Kill cohort/individual, transfer biomass to litter

			ppft.cmass_mort += indiv.ccont();

			log_info("%s killed by exceedance of bioclimatic limits...", (const char *)indiv.pft.name);
			indiv.kill(MORT_OTHER); //MORT

			vegetation.killobj();
			killed=true;
		}
		else {

			if (indiv.pft.lifeform==TREE) {

				// Calculate this year's growth efficiency
				// Eqn 31, Smith et al 2001

				if (!negligible(indiv.cmass_leaf))
					greff=max(indiv.anpp,0.0)/indiv.cmass_leaf/indiv.sla();
				else
					greff=0.0;

				// Calculate 5 year mean growth efficiency

				indiv.greff_5.add(greff);
				greff_mean = indiv.greff_5.mean();

				// BACKGROUND MORTALITY
				//
				// This is now modelled as an increasing probability of death with
				// age. The expected likelihood of death by "background" factors at a
				// given age is given by the general relationship [** = raised to the
				// power of]:
				//  (1) mort(age) = min( Z * ( age / longevity ) ** Q , 1)
				// where Z         is a constant (value derived below);
				//       Q         is a positive constant
				//                 (or zero for constant mortality with age)
				//       longevity is the age at which the fraction of the initial
				//                 cohort expected to survive is a known value, F
                // It is possible to derive the value of Z given values for longevity,
				// Q and F, by integration:
				//
				// The rate of change of cohort size (P) at any age is given by:
				//  (2) dP/dage = -mort(age) * P
				// Combining (1) and (2),
				//  (3) dP/dage = -Z.(age/longevity)**Q.P
				// From (3),
				//  (4) (1/P) dP = -Z.(age/longevity)**Q.dage
				// Integrating the LHS of eqn 4 from P_0 (initial cohort size) to
				// P_longevity (cohort size at age=longevity); and the RHS of eqn 4
				// from 0 to longevity:
				//   (5) LHS = ln(P_longevity) - ln(P_0)
				//   (6)     = ln(P_longevity/P_0)
				//   (7)     = ln(P_0*F/P_0)
				//   (8)     = ln(F)
				//   (9) RHS = integral[0,longevity] -Z.(age/longevity)**Q.dage
				//  (10)     = integral[0,longevity] -Z.(1/longevity)**Q.age**Q.dage
				//  (11)     = -Z.(1/longevity)**Q.integral[0,longevity] age**Q.dage
				//  (12)     = -Z.(1/longevity)**Q.longevity**(Q+1)/(Q+1)
				//             (Zwillinger 1996, p 360)
				// Combining (8) and (12),
				//  (13) Z = - ln(F) * (Q+1) / longevity
				// Combining (1) and (13),
				//  (14) mort(age) =
				//         min ( -ln(F) * (Q+1) / longevity * (age/longevity)**P, 1)

				// Eqn 14
				mort_min=min(1.0,KMORTBG_LNF*(KMORTBG_Q+1)/indiv.pft.longevity*
					pow(indiv.age/indiv.pft.longevity,KMORTBG_Q));


				// Growth suppression mortality
				// Smith et al 2001; c.f. Pacala et al 1993, Eqn 5

				// guess2008 - introduce a smoothly-varying mort_greff - 5 is the exponent in the global validation
				if (ifdailycstarvmort) {
					mort_greff = 0.0;
				} else if (ifsmoothgreffmort)
					mort_greff=KMORTGREFF/(1.0+pow((greff_mean/(indiv.pft.greff_min)),5.0));
				else {
					// Standard case, as in guess030124
					if (greff_mean<indiv.pft.greff_min)
						mort_greff=KMORTGREFF;
					else
						mort_greff=0.0;
				}

				//MORT
				double mort_greff_before_self_thin = mort_greff;

				// Increase growth efficiency mortality if summed crown area within
				// cohort exceeds 1 (to ensure self-thinning for shade-tolerant PFTs)

				if (vegmode==COHORT) {

					if (indiv.crownarea*indiv.densindiv>1.0)
						mort_greff=max((indiv.crownarea*indiv.densindiv-1.0)/
							(indiv.crownarea*indiv.densindiv),mort_greff);
				}

				indiv.mort_greff = mort_greff;
				indiv.mort_age = mort_min;

				// Overall mortality: c.f. Eqn 29, Smith et al 2001

				mort=mort_min+mort_greff-mort_min*mort_greff;


				// guess2008 - added safety check
				if (mort > 1.0 || mort < 0.0)
					fail("error in mortality_guess: bad mort value");

				if (ifstochmort && mt.stochmort) {

					// Impose stochastic mortality
					// Each individual in cohort dies with probability 'mort'

					nindiv=(int)(indiv.densindiv*patcharea+0.5);
					nindiv_prev=nindiv;

					for (i=0;i<nindiv_prev;i++)
						if (randfrac(stand.seed)<mort) nindiv--;

					if (nindiv_prev)
						frac_survive=(double)nindiv/(double)nindiv_prev;
					else
						frac_survive=0.0;
				}

				// Deterministic mortality (cohort mode only)

				else frac_survive=1.0-mort;

				ppft.cmass_mort += (1.0 - frac_survive) * indiv.ccont();

				// Reduce individual density and biomass on patch area basis
				// to account for loss of killed individuals

				//MORT
				double mort_greff_frac = mort ? mort_greff_before_self_thin / mort : 0.0;
				double mort_selfthin_frac = mort ? (mort_greff - mort_greff_before_self_thin) / mort : 0.0;
				double mort_age_frac = mort ? mort_min / mort : 0.0;
				double mort_other_frac = 1.0 - mort_greff_frac - mort_selfthin_frac - mort_age_frac;
				indiv.reduce_biomass(1.0 - frac_survive,
                                     0.0,
                                     mort_greff_frac * (1.0 - frac_survive),
                                     mort_selfthin_frac * (1.0 - frac_survive),
                                     (1.0-mort_greff_frac-mort_selfthin_frac) * (1.0 - frac_survive),
									 0.0, false,
				                     indiv.closs_mort, indiv.nloss_mort); // Calculate the age fraction as a residual to avoid rounding errors creeping in

				// Remove this cohort completely if all individuals killed
				// (in individual mode: removes individual if killed)

				if (negligible(indiv.densindiv)) {
					log_info("%s killed by background mortality. mort_greff=%f, mort_age=%f, total_mort_frac=%f...",
						    (const char *)indiv.pft.name, mort_greff, mort_min,
							mort);
					vegetation.killobj();
					killed=true;
				}

				// Update allometry

				else {
					if (!allometry(indiv)) {
						// Handle allometry related mortality
						ppft.cmass_mort += indiv.ccont();
						// ppft.cmass_mort_wood += indiv.cmass_wood();
						// if (indiv.diam >= 0.1)
						// 	ppft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

						indiv.kill(MORT_OTHER); //MORT

						vegetation.killobj();
						killed = true;
					}
				}
			} 
			else { 
				// MORTALITY in PEATLAND stands
				// i.e. shading mortality when LAI gets too high (cf mortality_lpj above), whereupon
				// biomass is reduced (ends up in litter) to be consistent with the given LAI limit (maxlai_peatland)
				if (patch.stand.landcover == PEATLAND && wetlandlailimit > 0.0) {

					double mort_shade = 0.0;

					double fpc_peat = fpc_moss + fpc_grass + fpc_other;

					if (!negligible(indiv.cmass_leaf)) {

						// Max values
						double maxlai_peatland = wetlandlailimit; // set from .ins file
						
						if (patch.stand.is_true_wetland_stand())
							maxlai_peatland *= 2; // Allow higher LAI values for wetlands south of PEATLAND_WETLAND_LATITUDE_LIMIT N as these are more productive  

						double maxfpc_peatland = 1.0 - exp(-0.5 * maxlai_peatland);

						// This individual

						if (indiv.fpc > maxfpc_peatland) {
							// fracmass_lpj calculates and returns new biomass as a fraction of old biomass given an FPC
							// reduction from fpc_high to fpc_low, assuming LPJ allometry
							double fpc_dec = (indiv.fpc - maxfpc_peatland) * indiv.fpc / fpc_peat;
							mort_shade = 1.0 - fracmass_lpj(indiv.fpc - fpc_dec, indiv.fpc, indiv);
						}

						// Reduce C biomass to account for biomass lost through shading mortality
						if (mort_shade> 0.0) {
							indiv.reduce_biomass(0.0, 0.0, 0.0, mort_shade, 0.0,
												 0.0, false,
											     indiv.closs_shade,
												 indiv.nloss_shade);//MORT			: first parameter is zero !! -> NO MORTALITY
							if (!allometry(indiv)) {
								// Handle allometry related mortality
								ppft.cmass_mort += indiv.ccont();
								// ppft.cmass_mort_wood += indiv.cmass_wood();
								// if (indiv.diam >= 0.1)
								// 	ppft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

								indiv.kill(MORT_OTHER); //MORT

								vegetation.killobj();
								killed = true;
							}
						}

					} // correct stand and grass type

				} // limit lai only
			}
		}

		// ... on to next individual

		if (!killed) vegetation.nextobj();
	}
}

void mortality_daily(Patch& patch){
	// TODO Assaf :)
}

void insect_mortality_guess(Stand& stand, Patch& patch) {
	if (!ifinsect) {
		return;
	}

	// Fraction of cohort (or individual) surviving
	double frac_survive;
	bool killed;
	Vegetation& vegetation=patch.vegetation;

	// INSECT MORTALITY
	//With reference to tables, figures and equations from Lagergren et al. 2025
	// https://doi.org/10.5194/gmd-18-8071-2025

	double mortality_insect;

	// Parameter setting in the "default expert judgement" simulation (Table 2)
	//const double insect_factor = 0.0030; //Callibrated factor to scale damage, set to 0.003 with SSC and 0.002 without
	//const double insect_diam_lim = 0.15; // 15 cm is concidered as a low limit for Ips Typographus
	//const double kgc_min = -2.0; // Min of f(population index gridcell)
	//const double kgc_max = 1.0; // Max of f(population index gridcell)
	//const double k1 = 0.05; // Parameter related to intercept of f(population index gridcell)
	//const double k2 = 0.25; // Parameter related to slope of f(population index gridcell)
	//const double kbase_bm = 0.0001; // Base level of insect stem C
	//const double kp_min = -3.8 - kgc_min; // Min of f(population index patch / brood material)
	//const double kp_max = 1.0; // Max of f(population index patch / brood material)
	//const double k3 = 0.03; // Parameter related to intercept of f(population index patch / brood material)
	//const double k4 = 0.4; // Parameter related to slope of f(population index patch / brood material)
	//const double k_pyw_wscal = 3.0; // Weight of previous year's wscal in relation to current year
	//const double kdrought_max = 2.0; // Max of f(water stress)
	//const double k5 = 8.0; // Slope of f(water stress)
	//// Constant in f(phenology)const double k6 = 3.33 in guess.cpp;
	//// Constant in f(phenology)const double k7 = 75.0 in guess.cpp;
	// Salvage cutting parameters (right column of page 8076), set a first_salvcutyear value outside simulated years to disable
	//const int first_salvcutyear = 1920; // First year to reduce storm litter assuming salvage harvest of storm felled trees
	//const double min_salvcut = 5.0 * 0.0001 * 400 * 0.5; // Min C mass for salvage cutting (m3 ha-1 * ha m-2 * kg m-3 kg C kg-1)
	//const double salvcutcapacity = 0.5; // Fraction of 10-yr annual average harvest capacity that can be used for salvage cutting
	//const double salvcutfrac = 0.9; // Fraction of wood removed in salvage cutting
	//const double sankillfrac = 0.25; // Fraction of bark beetles killed in sanitary counter-measures
	//const double sanmin = 0.01; // Min fraction of killed spruce biomass for sanitary cutting to take place

	// Calibrated parameters using SSC (Table 3)
	const double insect_factor = 0.005124; //Callibrated factor to scale damage, set to 0.003 with SSC and 0.002 without
	const double insect_diam_lim = 0.15; // 15 cm is concidered as a low limit for Ips Typographus
	const double kgc_min = -1.75; // Min of f(population index gridcell)
	const double kgc_max = 1.0; // Max of f(population index gridcell)
	const double k1 = 0.005; // Parameter related to intercept of f(population index gridcell)
	const double k2 = 0.5; // Parameter related to slope of f(population index gridcell)
	const double kbase_bm = 0.0001; // Base level of insect stem C
	const double kp_min = -3.8 - kgc_min; // Min of f(population index patch / brood material)
	const double kp_max = 1.0; // Max of f(population index patch / brood material)
	const double k3 = 0.15; // Parameter related to intercept of f(population index patch / brood material)
	const double k4 = 0.8; // Parameter related to slope of f(population index patch / brood material)
	const double k_pyw_wscal = 0.25; // Weight of previous year's wscal in relation to current year
	const double kdrought_max = 2.0; // Max of f(water stress)
	const double k5 = 5.0; // Slope of f(water stress)
	// Constant in f(phenology)const double k6 = 3.33 in guess.cpp;
	// Constant in f(phenology)const double k7 = 75.0 in guess.cpp;
	// Salvage cutting parameters (right column of page 8076), set a first_salvcutyear value outside simulated years to disable
	const int first_salvcutyear = 1920; // First year to reduce storm litter assuming salvage harvest of storm felled trees, was 1920 in 230113 run
	const double min_salvcut = 5.0 * 0.0001 * 400.0 * 0.5; // Min C mass for salvage cutting (m3 ha-1 * ha m-2 * kg m-3 kg C kg-1)
	const double salvcutcapacity = 0.5; // Fraction of 10-yr annual average harvest capacity that can be used for salvage cutting
	const double salvcutfrac = 0.9; // Fraction of wood removed in salvage cutting
	const double sankillfrac = 0.25; // Fraction of bark beetles killed in sanitary counter-measures
	const double sanmin = 0.01; // Min fraction of killed spruce biomass for sanitary cutting to take place

	patch.insect_dam_wood = 0.0;

	// Loop over ift
	for (int iftid = 0;iftid < nift;iftid++) {
		Ift& ift = iftlist[iftid];
		Gridcellift& gcift = stand.get_gridcell().ift[iftid];
		Patchift& pift = patch.ift[iftid];

		double insect_litter = 0.0;
		double drought_ind = 0.0;
		double wscal_patchinsect = 0.0;

		// First loop over pft to get sum of insect stem litter
		for (int p = 0; p < npft; p++) {
			// Retrieve next patch PFT
			Patchpft& ppft = patch.pft[p];

			if (!patch.stand.pft[p].active)
				continue;

			if (ift.pftinhosts(ppft.pft.name)) {
				insect_litter += ppft.insect_stem_litter_prev;
			}
		}

		// Ruduce insect litter if salvage cutting is applied (right column of page 8076)
		double insect_litter_before = insect_litter; //Insect litter before SC, needed for weighting below
		if (date.get_calendar_year() >= first_salvcutyear && insect_litter > min_salvcut) {
			double maxsalvcut = stand.get_gridcell().flux_charvest10y.mean() * salvcutcapacity;
			double dam_wood_available = stand.get_gridcell().storm_dam_wood_prev * salvcutfrac;
			insect_litter *= 1.0 - min(salvcutfrac, salvcutfrac * maxsalvcut / dam_wood_available);
		}

		// Second loop over pfts to get drought index (Fig 1D) weighted by insect stem litter
		for (int p = 0; p < npft; p++) {
			// Retrieve next patch PFT
			Patchpft& ppft = patch.pft[p];

			if (!patch.stand.pft[p].active)
				continue;

			if (ift.pftinhosts(ppft.pft.name)) {
				// Drought index (0-2) calculated from average wscal May-July from previous and current year,
				// preivios year is given a weight of three following Marini et al 2017 (Table 2, 1st ranked model)
				// Linear treshold function between 0.6 and 1.0
				double pft_weight; //If possible weight by litter before salvage cutting
				if (insect_litter_before > 0.0)
					pft_weight = ppft.insect_stem_litter_prev / insect_litter_before;
				else
					pft_weight = 1.0 / (double)gcift.npfthosts;
				double wscal_insect = ((ppft.wscal_mean_month_prev[4] + ppft.wscal_mean_month_prev[5] + ppft.wscal_mean_month_prev[6]) * k_pyw_wscal +
					(ppft.wscal_mean_month[4] + ppft.wscal_mean_month[5] + ppft.wscal_mean_month[6])) /
					(3.0 * k_pyw_wscal + 3.0); //Eq. 11 without first making average over the current and previous year
				drought_ind += min(kdrought_max, (1.0 - wscal_insect) * k5) * pft_weight; //Eq. 10
				wscal_patchinsect += wscal_insect * pft_weight;

			}
			ppft.insect_stem_litter_prev = ppft.insect_stem_litter;
			ppft.insect_stem_litter = 0.0;
		}
		drought_ind = max(0.0, min(2.0, drought_ind)); //Make sure that it is within the 0-2 range

		// Get stem C mass available to insect for host trees with diameter > insect_diam_lim
		double insect_stem = 0.0;
		// Loop through individuals
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();
			Patchpft& patchpft = patch.pft[indiv.pft.id];
			// For this individual ...
			if (ift.pftinhosts(indiv.pft.name)) {

				if (indiv.diam > insect_diam_lim)
					insect_stem += (indiv.cmass_sap + indiv.cmass_heart - indiv.cmass_debt);
			}
			vegetation.nextobj();
		}

		// Update insect patch population index

		// Substitute for wind damage cmass Can be used for simulation with random insect litter without connection to real mortality
		const double max_stormfrac = 0.50; // Maximum fraction of available stem C that can be used after storm, higher if patch level storm seed?
		//double storm_sub = max(0.0, pow((randfrac(stand.seed) - 0.75) * 4.0, 3.0)) * insect_stem * max_stormfrac; // Unique seed for each pft-ipf
		double storm_sub = max(0.0, pow((gcift.storm_rand - 0.75) * 4.0, 3.0)) * insect_stem * max_stormfrac; // Common seed for ift at gridcell
		storm_sub = 0.0; // Outcomment to use. Should be a real switch!!

		// Index of general conditions for population growth related to population size and availability of dead trees (Fig 1C)
		double base_insect_ind = min(kp_max, max(kp_min, -log(pift.pop_ind / (insect_litter + storm_sub + kbase_bm) * k3) * k4)); // Eq. 8 & 9

		// Index related to average gridcell level average population index (Fig 1B)
		double gridcell_ind = min(kgc_max, max(kgc_min, -log(gcift.pop_ind * k1) * k2)); // Eq. 7

		// Growt rate based on the components
		double R = drought_ind + base_insect_ind + gridcell_ind + gcift.phen_ind; // Eq. 4
		pift.pop_ind *= exp(R); // New population size, Eq. 5

		// Reduce population if sanitary cutting is applied (right column of page 8076)
		if (date.get_calendar_year() >= first_salvcutyear &&
			(insect_factor * pift.pop_ind / insect_stem) > sanmin &&
			insect_stem > 0.0) {
			pift.pop_ind *= (1.0 - sankillfrac);
		}

		//Contribution to mean gridcell level pop_ind used next year
		gcift.pop_ind_new += pift.pop_ind / (double)stand.npatch() * stand.get_gridcell_fraction();


		// Calculate mortality insect for this ift and pft if there is available C mass
		mortality_insect = 0.0;

		if (insect_stem > 0.0) {
			const double max_insect_mort = 0.75;

			mortality_insect = min(max_insect_mort * insect_stem, insect_factor * pift.pop_ind) / insect_stem; // Eq. 6
			frac_survive = 1.0 - mortality_insect;

			if (insect_mortality && mortality_insect > 0.0
				&& (date.get_calendar_year() >= insect_mortality_firstyear) && (date.get_calendar_year() <= insect_mortality_lastyear)
				&& !(readdist_barkbeetle && date.get_calendar_year() >= disturbance_input_firstyear && date.get_calendar_year() <= disturbance_input_lastyear)) {

				// Apply mortality insect if d>limit && indiv.pft == ift.host_pft
				// Loop through individuals
				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv = vegetation.getobj();
					Patchpft& patchpft = patch.pft[indiv.pft.id];

					// For this individual ...
					killed = false;
					if (ift.pftinhosts(indiv.pft.name)) {

						if (indiv.diam > insect_diam_lim) {

							patchpft.cmass_barkbeetle_killed += mortality_insect * indiv.ccont();
							patchpft.cmass_mort +=  mortality_insect * indiv.ccont();
							// patchpft.cmass_mort_wood += mortality_insect * indiv.cmass_wood();
							if (indiv.diam >= 0.1) {
								patchpft.cmass_barkbeetle_killed_plus10cm += mortality_insect * indiv.ccont();
								// patchpft.cmass_mort_wood_plus10cm += mortality_insect * indiv.cmass_wood();
							}

							gcift.insect_damage += indiv.cmass_wood() * (1.0 - frac_survive) / (double)stand.npatch() * stand.get_gridcell_fraction();
							patch.insect_dam_wood += indiv.cmass_wood() * (1.0 - frac_survive);

							// Reduce individual density and biomass on patch area basis
							// to account for loss of killed individuals
							double closs, nloss;
							indiv.reduce_biomass(0.0, 0.0, 0.0, 0.0, 0.0,
												 1.0 - frac_survive, false, closs, nloss);

							// Remove this cohort completely if all individuals killed
							// (in individual mode: removes individual if killed)

							if (negligible(indiv.densindiv)) {
								vegetation.killobj();
								killed = true;
							}

							// Update allometry

							else {
								// Update allometry
								if (!allometry(indiv)) {
									// Handle allometry related mortality
									patchpft.cmass_mort += indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_wood();
									// patchpft.cmass_mort_wood += indiv.cmass_wood();
									// if (indiv.diam >= 0.1)
									// 	patchpft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

									indiv.kill(MORT_OTHER);

									vegetation.killobj();
									killed = true;
								}
							}
						}

					}

					if (!killed) vegetation.nextobj();
				}
			}
		}
	}
}

void storm_mortality_guess(Stand& stand, Patch& patch) {
	if (!ifstorm) {
		return;
	}

	// Fraction of cohort (or individual) surviving
	double frac_survive;
	bool killed;
	Vegetation& vegetation=patch.vegetation;

	// STORM MORTALITY

	patch.storm_dam_wood = 0.0;
	patch.storm_dam_wood_pot = 0.0;
	patch.storm_cmass_wood_save = 0.0;
	patch.storm_sens_ind = 0.0;

	Gridcell& gridcell = stand.get_gridcell();

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patchpft& patchpft = patch.pft[indiv.pft.id];

		// For this individual ...

		killed = false;

		double mortality_storm = max(0.0, indiv.dam_frac - indiv.man_strength); // Also reduce fraction for other mortality?
		double mortality_storm_pot = max(0.0, indiv.dam_frac_pot - indiv.man_strength);	// deduct man_strength also here ?

		// Updatate potential damaged wood at different levels, used for output
		patch.storm_dam_wood_pot += mortality_storm_pot * indiv.cmass_wood();
		patch.storm_cmass_wood_save += indiv.cmass_wood();
		patch.storm_sens_ind += indiv.sens_ind_tot * indiv.cmass_wood();				// cmass_wood weight

		frac_survive = 1.0 - mortality_storm;

		if (storm_mortality && mortality_storm > 0.0 
			&& (date.get_calendar_year() >= storm_mortality_firstyear) && (date.get_calendar_year() <= storm_mortality_lastyear)
			&& !(readdist_wind && date.get_calendar_year() >= disturbance_input_firstyear && date.get_calendar_year() <= disturbance_input_lastyear)) {

			double mort_cmass = (1.0 - frac_survive) * indiv.ccont();
			double mort_cmass_wood = (1.0 - frac_survive) * indiv.cmass_wood();
				
			patchpft.cmass_wind_killed += mort_cmass;
			patchpft.cmass_mort += mort_cmass;
			// patchpft.cmass_mort_wood += mort_cmass_wood;
			if (indiv.diam >= 0.1) {
				patchpft.cmass_wind_killed_plus10cm += mort_cmass;
				// patchpft.cmass_mort_wood_plus10cm += mort_cmass_wood;
			}

			// Reduce individual density and biomass on patch area basis
			// to account for loss of killed individuals
			double closs, nloss;
			indiv.reduce_biomass(1.0 - frac_survive, 0.0, 0.0, 0.0, 0.0, 0.0,
							     false, closs, nloss);

			// Updatate damaged wood at different levels, used for output
			patch.storm_dam_wood += mort_cmass_wood;
			gridcell.storm_dam_wood += mort_cmass_wood / (double)stand.npatch() * stand.get_gridcell_fraction();	// zeroed in storm(), needed for saving previous year's value

			//dprintf("Stormres %5.2f %5.2f %i %i %i %i %i %3.1f %6.5f\n", 
			//	stand.get_gridcell().get_lon(), stand.get_gridcell().get_lat(), date.get_calendar_year(), stand.id, patch.id, indiv.id,
			//	indiv.age, indiv.height, mortality_storm);

			// Remove this cohort completely if all individuals killed
			// (in individual mode: removes individual if killed)

			if (negligible(indiv.densindiv)) {
				vegetation.killobj();
				killed = true;
			}

			// Update allometry

			else {
				// Update allometry
				if (!allometry(indiv)) {
					// Handle allometry related mortality
					patchpft.cmass_mort += indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_wood();
					// patchpft.cmass_mort_wood += indiv.cmass_wood();
					// if (indiv.diam >= 0.1)
					// 	patchpft.cmass_mort_wood_plus10cm += indiv.cmass_wood();

					indiv.kill(MORT_OTHER);

					vegetation.killobj();
					killed = true;
				}
			}

		}

		// ... on to next individual

		if (!killed) vegetation.nextobj();
	}

	if (patch.cmass_wood() > 0.0) {
		patch.storm_sens_ind /= patch.cmass_wood();	// cmass_wood weight
	}
	else {
		patch.storm_sens_ind = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// Simplified update of insect dynamics after disturbance, done in mortality_guess if no disturbance
// Internal function (do not call directly from framework)

void insectupdate(Stand& stand, Patch& patch) {
	if (!ifinsect) {
		return;
	}

	for (int iftid = 0;iftid < nift;iftid++) {
		Ift& ift = iftlist[iftid];
		Gridcellift& gcift = stand.get_gridcell().ift[iftid];
		Patchift& pift = patch.ift[iftid];

		double insect_litter = 0.0;
		double drought_ind = 0.0;

		// Loop over pft to get insect stem litter
		for (int p = 0; p < npft; p++) {
			// Retrieve next patch PFT
			Patchpft& ppft = patch.pft[p];

			if (!patch.stand.pft[p].active)
				continue;

			ppft.insect_stem_litter_prev = ppft.insect_stem_litter;
			ppft.insect_stem_litter = 0.0;
		}

		const double kgc_min = -1.75; // Min of f(population index gridcell)
		const double kgc_max = 1.0; // Max of f(population index gridcell)
		const double k1 = 0.005; // Parameter related to intercept of f(population index gridcell)
		const double k2 = 0.5; // Parameter related to slope of f(population index gridcell)
		const double kbase_bm = 0.0001; // Base level of insect stem C
		const double kp_min = -3.8 - kgc_min; // Min of f(population index patch / brood material)
		const double kp_max = 1.0; // Max of f(population index patch / brood material)
		const double k3 = 0.15; // Parameter related to intercept of f(population index patch / brood material)
		const double k4 = 0.8; // Parameter related to slope of f(population index patch / brood material)

		const double storm_sub = 0.0; // Outcomment to use. Should be a real switch!!


		// Index of general conditions for population growth related to population size and availability of dead trees (Fig 1C)
		double base_insect_ind = min(kp_max, max(kp_min, -log(pift.pop_ind / (insect_litter + storm_sub + kbase_bm) * k3) * k4)); // Eq. 8 & 9

		// Index related to average gridcell level average population index (Fig 1B)
		double gridcell_ind = min(kgc_max, max(kgc_min, -log(gcift.pop_ind * k1) * k2)); // Eq. 7

		// Growt rate based on the components
		double R = drought_ind + base_insect_ind + gridcell_ind + gcift.phen_ind; // Eq. 4
		pift.pop_ind *= exp(R); // New population size, Eq. 5

		//Contribution to mean gridcell level pop_ind used next year
		gcift.pop_ind_new += pift.pop_ind / (double)stand.npatch() * stand.get_gridcell_fraction();
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// FIRE DISTURBANCE
// Internal function (do not call directly from framework)

void fire(Patch& patch, double& fireprob) {

	// DESCRIPTION
	// Calculates probability/incidence of disturbance by fire; implements
	// volatilisation of above-ground litter due to fire.
	// Mortality due to fire implemented in mortality functions, not here.
	//
	// Reference: Thonicke et al 2001, with updated Eqn 9

	// OUTPUT PARAMETER
	// fireprob = probability of fire in this patch this year
	//            (in population mode: fraction of modelled area affected by fire)

	// Minimum total aboveground litter required for fire (kgC/m2)
	const double MINFUEL = 0.2;
	// Total aboveground litter (kgC/m2) [litter_ag_total/1000, fuel/1000; LPJF]
	double litter_ag;
	// Mean litter flammability moisture threshold [=moisture of extinction, me;
	// Thonicke et al 2001 Eqn 2; moistfactor; LPJF]
	double me_mean;
	// Probability of fire on a given day [p(m); Thonicke et al 2001, Eqn 2;
	// fire_prob; LPJF]
	double pm;
	// Summed length of fire season (days) [N; Thonicke et al 2001, Eqn 4;
	// fire_length; LPJF]
	double n;
	// Annual mean daily probability of fire [s; Thonicke et al 2001, Eqn 8;
	// fire_index; LPJF]
	double s;

	double sm; // s-1
	double mort_fire; // fire mortality as fraction of current FPC
	int p;

	// Calculate fuel load (total aboveground litter)

	litter_ag = 0.0;

	// Soil
	Soil& soil = patch.soil;

	// Top soil layer
	Soillayer& surflayer = soil.soillayers[0];

	// Loop through PFTs

	for (p = 0; p < npft; p++) {
		Patchpft& ppft = patch.pft[p];
		litter_ag += ppft.cmass_litter_leaf + ppft.cmass_litter_sap + ppft.cmass_litter_heart +
			ppft.cmass_litter_repr;
	}

	// Surface Litter
	litter_ag += surflayer.sompool[SURFSTRUCT].cmass + surflayer.sompool[SURFMETA].cmass +
		surflayer.sompool[SURFFWD].cmass + surflayer.sompool[SURFCWD].cmass;

	// Prevent fire if fuel load below prescribed threshold

	if (litter_ag < MINFUEL) {
		fireprob = 0.0;
		return;
	}

	// Calculate mean litter flammability moisture threshold

	me_mean = 0.0;

	// Loop through PFTs

	for (p = 0; p < npft; p++) {
		Patchpft& ppft = patch.pft[p];
		me_mean += (ppft.cmass_litter_leaf + ppft.cmass_litter_sap + ppft.cmass_litter_heart +
			ppft.cmass_litter_repr) * ppft.pft.litterme / litter_ag;
	}

	// Surface litter
	me_mean += (surflayer.sompool[SURFSTRUCT].cmass * surflayer.sompool[SURFSTRUCT].litterme +
		surflayer.sompool[SURFMETA].cmass * surflayer.sompool[SURFMETA].litterme +
		surflayer.sompool[SURFFWD].cmass * surflayer.sompool[SURFFWD].litterme +
		surflayer.sompool[SURFCWD].cmass * surflayer.sompool[SURFCWD].litterme) / litter_ag;

	// Calculate length of fire season in days
	// Eqn 2, 4, Thonicke et al 2001
	// NOTE: error in Eqn 2, Thonicke et al - multiplier should be PI/4, not PI

	n = 0.0;
	for (int mth = 0; mth < 12; mth++) {

		// No fires unless there is some soil melt.
		if (hydrology == GERTEN || (soil.mthaw[mth] > 0.0 && date.year > FIRST_FREEZE_YEAR + 10)) {
			pm = exp(-PI * surflayer.mwcont[mth] / me_mean * surflayer.mwcont[mth] / me_mean) * date.ndaymonth[mth];		// Eqn 2
		}
		else {
			pm = 0.0;
		}


		// Eqn 4
		n += pm;
	}

	// Calculate fraction of grid cell burnt
	// Thonicke et al 2001, Eqn 9

	s = n / (double)date.year_length();
	sm = s - 1;

	fireprob = s * exp(sm / (0.45 * sm * sm * sm + 2.83 * sm * sm + 2.96 * sm + 1.04));

	if (fireprob > 1.0)
		fail("fire: probability of fire >1.0");
	else if (fireprob < 0.001) fireprob = 0.001; // c.f. LPJF

	// Calculate expected flux from litter due to fire
	// (fluxes from vegetation calculated in mortality functions)

	Fluxes& fluxes = patch.fluxes;

	for (p = 0; p < npft; p++) {

		Patchpft& ppft = patch.pft[p];

		// Assume litter is as fire resistant as the vegetation it originates from

		mort_fire = fireprob*(1.0 - ppft.pft.fireresist);

		// Calculate flux from burnt litter

		fluxes.report_flux(Fluxes::FIREC,
			mort_fire*(ppft.cmass_litter_leaf + ppft.cmass_litter_sap +
				ppft.cmass_litter_heart + ppft.cmass_litter_repr));

		report_fire_nfluxes(fluxes, mort_fire * (ppft.nmass_litter_leaf +
			ppft.nmass_litter_sap + ppft.nmass_litter_heart));

		// Account for burnt above ground litter

		ppft.cmass_litter_leaf  *= (1.0 - mort_fire);
		ppft.cmass_litter_sap   *= (1.0 - mort_fire);
		ppft.cmass_litter_heart *= (1.0 - mort_fire);
		ppft.cmass_litter_repr  *= (1.0 - mort_fire);
		ppft.nmass_litter_leaf  *= (1.0 - mort_fire);
		ppft.nmass_litter_sap   *= (1.0 - mort_fire);
		ppft.nmass_litter_heart *= (1.0 - mort_fire);
	}

	// Soil litter
	double mort_fire_struct = fireprob * (1.0 - surflayer.sompool[SURFSTRUCT].fireresist);
	double mort_fire_meta = fireprob * (1.0 - surflayer.sompool[SURFMETA].fireresist);
	double mort_fire_fwd = fireprob * (1.0 - surflayer.sompool[SURFFWD].fireresist);
	double mort_fire_cwd = fireprob * (1.0 - surflayer.sompool[SURFCWD].fireresist);

	// Calculate flux from burnt soil litter
	fluxes.report_flux(Fluxes::FIREC, surflayer.sompool[SURFSTRUCT].cmass * mort_fire_struct +
		surflayer.sompool[SURFMETA].cmass * mort_fire_meta +
		surflayer.sompool[SURFFWD].cmass * mort_fire_fwd +
		surflayer.sompool[SURFCWD].cmass * mort_fire_cwd);

	// Account for burnt above ground litter
	surflayer.sompool[SURFSTRUCT].cmass *= (1.0 - mort_fire_struct);
	surflayer.sompool[SURFMETA].cmass *= (1.0 - mort_fire_meta);
	surflayer.sompool[SURFFWD].cmass *= (1.0 - mort_fire_fwd);
	surflayer.sompool[SURFCWD].cmass *= (1.0 - mort_fire_cwd);

	// Calculate nitrogen flux from burnt soil litter
	double nflux_fire = surflayer.sompool[SURFSTRUCT].nmass * mort_fire_struct +
		surflayer.sompool[SURFMETA].nmass * mort_fire_meta +
		surflayer.sompool[SURFFWD].nmass * mort_fire_fwd +
		surflayer.sompool[SURFCWD].nmass * mort_fire_cwd;

	report_fire_nfluxes(fluxes, nflux_fire);

	surflayer.sompool[SURFSTRUCT].nmass *= (1.0 - mort_fire_struct);
	surflayer.sompool[SURFMETA].nmass *= (1.0 - mort_fire_meta);
	surflayer.sompool[SURFFWD].nmass *= (1.0 - mort_fire_fwd);
	surflayer.sompool[SURFCWD].nmass *= (1.0 - mort_fire_cwd);
}


///////////////////////////////////////////////////////////////////////////////////////
// DISTURBANCE
// Generic patch-destroying disturbance with a prescribed probability
// Internal function - do not call from framework

void disturbance(Patch& patch, double disturb_prob) {

	// DESCRIPTION
	// Destroys all biomass in a patch with a certain stochastic probability.
	// Biomass enters the litter, which is not affected by the disturbance.
	// NB: cohort and individual mode only
	// This is the generic patch-destroying disturbance pathway used by
	// stochastic disturbance and by Dave deterministic disturbance. Do not
	// confuse this with prescribed_disturbance(), which allocates prescribed
	// disturbed area across stands and usually routes non-fire disturbance
	// through harvest_forest().

	// INPUT PARAMETER
	// disturb_prob = the probability of a disturbance this year

	if (randfrac(patch.stand.seed)<disturb_prob) {

		Vegetation& vegetation = patch.vegetation;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();
			Patchpft& ppft = patch.pft[indiv.pft.id];

			if(indiv.alive)
				ppft.cmass_dist += indiv.ccont();

			indiv.mortflux(1.0, MORT_DIST, patch.stand.get_gridcell_fraction() / (double)patch.stand.nobj, true);
            		
			if(ifback_disturb_as_harv){
                		harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0.1);
            		}
            		indiv.kill(MORT_DIST); //MORT mortmech = Disturbance (4)

			vegetation.killobj();
		}

		patch.disturbed = true;
		patch.age = 0;
		log_info("Disturbance...");
	}
	else {
		patch.disturbed = false;
	}
}

/// PLANTING
// Planting of trees in managed forest

void planting(Patch& patch) {

	double ltor;
	Stand& stand = patch.stand;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();

		Standpft& spft = patch.stand.pft[pft.id];
		bool pft_selection = spft.plant;

		if(pft_selection) {

			Standpft& spft = stand.pft[pft.id];
			// plants per ha., could be modified by stand productivity
			double plantdensity = 1000;	// temporary default value if not specified for pft
			if(spft.plantdensity >= 0.0) {
				plantdensity = spft.plantdensity;
			}
			else if(pft.plantdensity) {
				plantdensity = pft.plantdensity;
			}

			if(plantdensity) {
				Individual& indiv = patch.vegetation.createobj(pft, patch.vegetation);
				indiv.densindiv = plantdensity / 10000.0;	// per ha to per m2
				ltor = pft.ltor_max;
				allocation_init(PLANTSIZE, ltor, indiv);
				allometry(indiv);

				// Calculate storage pool size
				indiv.max_n_storage = indiv.cmass_sap * indiv.pft.fnstorage / indiv.pft.cton_leaf_avr;
				indiv.scale_n_storage = indiv.max_n_storage * indiv.pft.cton_leaf_avr / PLANTSIZE;
			}
		}

		patch.pft[pft.id].nsapling=0.0;
		patch.pft[pft.id].wscal_mean_est=0.0;
		patch.pft[pft.id].anetps_ff_est=0.0;

		pftlist.nextobj();
	}
}

void establish_input(Patch& patch, bool pool_soil = false) {

	Stand& stand = patch.stand;
	Gridcellst& gcst = stand.get_gridcell().st[stand.stid];

	if (!gcst.prescribed_establishment || gcst.establishment_list.establishment_year != date.get_calendar_year())
		return;

	bool print_init_info = false;
	double cmass_veg_tree_tot_old = check_harvest_cmass(patch);

	if(!pool_soil) {

		Vegetation& vegetation = patch.vegetation;
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();
			Patchpft& patchpft = patch.pft[indiv.pft.id];

			// Output of killed tree cmass in connection with tree initialisation turned off for now, probably not needed.
//			if(indiv.alive)
//				patchpft.cmass_dist += indiv.ccont();

			harvest_wood(indiv, 1.0, 1.0, 1.0, 1.0, false, true);
			indiv.kill(MORT_OTHER);	// Only grass biomass left at this stage.

			vegetation.killobj();
		}
	}

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();
		Gridcellpft& gcpft = stand.get_gridcell().pft[pft.id];

		for (int i = 0;i < MAX_NESTABLISH_INDIV;i++) {

			xtring pftname = gcst.establishment_list.est_array[i].pft_name;

			if (pftname == "")
				break;
			
			if (pft.name == pftname && pft.is_grass()) {

				if(print_init_info)
					dprintf("\n%s individual no. %d created in st %s\n", (char*)pftname, i, (char*)stlist[stand.stid].name);

				const bool use_latosa_20y_mean = true;
				const bool use_qsap_20y_mean = true;
				const bool use_age_20y_mean = false;
				const bool use_cton_20y_mean = true;

				Individual& indiv = patch.vegetation.createobj(pft, patch.vegetation);

				const double M_TO_CM = 100.0;
				const double C_TO_TOTAL_BIOMASS = 2.0;
				const double TOTAL_TO_AG_BIOMASS = 0.8;
				const double TON_TO_KG = 1000.0;

				double diam = stand.get_gridcell().st[stand.stid].establishment_list.est_array[i].diam / M_TO_CM;	// cm to m
				double vol = (PI / 4.0) * pft.k_allom2 * pow(diam, 2.0 + pft.k_allom3);
				double cmass_wood_ind = vol * indiv.pft.wooddens;

				if(forest_init_from_biomass) {
					cmass_wood_ind = stand.get_gridcell().st[stand.stid].establishment_list.est_array[i].biomass / C_TO_TOTAL_BIOMASS / TOTAL_TO_AG_BIOMASS * TON_TO_KG;	// Total cmass_wood from AG biomass)
					diam = pow((4.0 * cmass_wood_ind /(PI * pft.k_allom2 * indiv.pft.wooddens)), (1.0 / (2.0 + pft.k_allom3)));
					vol = cmass_wood_ind / indiv.pft.wooddens;
				}

				double height = pft.k_allom2 * pow(diam, pft.k_allom3);
				double treedens = gcst.establishment_list.est_array[i].dens;
				double f_sap = 0.6625 - 0.4372 * diam;
				double q_sap = ( 2.0 * f_sap - (f_sap * f_sap)) / 4.0;

				if(print_init_info)
					dprintf("initial q_sap calculated (diam=%.2f): %f\n", diam, q_sap);

				if (use_qsap_20y_mean) {
					int diam_bin = -1;
					int diam_bin_reserve = -1;
					for (int d=(N_DIAM_CLASSES_SAMP - 1); d>=0; d--) {
						if(gcpft.qsap_avr_20[d].size())
							diam_bin_reserve = d;
						if(diam * 100 > d * DIAM_CLASS_SIZE_SAMP && gcpft.qsap_avr_20[d].size()) {
							diam_bin = d;
							break;
						}
					}
					if(diam_bin > -1) {
						q_sap = gcpft.qsap_avr_20[diam_bin].mean();
						double diam_bin_mean = ((double)diam_bin + 0.5) * DIAM_CLASS_SIZE_SAMP;
						if((diam * 100) < diam_bin_mean) {
							if(gcpft.qsap_avr_20[diam_bin - 1].size() && (diam_bin - 1) >= 0) {
								double diam_bin_mean_prev = ((double)(diam_bin - 1) + 0.5) * DIAM_CLASS_SIZE_SAMP;
								double q_sap_prev = gcpft.qsap_avr_20[diam_bin - 1].mean();
								q_sap -= (diam_bin_mean - (diam * 100)) / (diam_bin_mean - diam_bin_mean_prev) * (q_sap - q_sap_prev);
							}
						}
						else {
							if(gcpft.qsap_avr_20[diam_bin + 1].size() && (diam_bin + 1) < N_DIAM_CLASSES_SAMP) {
								double diam_bin_mean_next = ((double)(diam_bin + 1) + 0.5) * DIAM_CLASS_SIZE_SAMP;
								double q_sap_next = gcpft.qsap_avr_20[diam_bin + 1].mean();
								q_sap -= (diam_bin_mean - (diam * 100)) / (diam_bin_mean - diam_bin_mean_next) * (q_sap - q_sap_next);
							}
						}
						if(print_init_info)
							dprintf("q_sap used from diam bin %d: %f\n", diam_bin, q_sap);
					}
					else if(diam_bin_reserve > -1) {
						q_sap = gcpft.qsap_avr_20[diam_bin_reserve].mean();
						if(print_init_info)
							dprintf("q_sap used from diam bin %d: %f\n", diam_bin_reserve, q_sap);
					}
				}	// else use the value calculated above

				int age = gcst.establishment_list.est_array[i].age;

				if (!age) {		// age value set to 0 in getforest_init() if not in input file

					age = (int)(100.0 * diam);

					if (use_age_20y_mean) {
						int diam_bin = -1;
						int diam_bin_reserve = -1;

						/// Using running means of N_DIAM_CLASSES_SAMP diameter classes
						for (int d=(N_DIAM_CLASSES_SAMP - 1); d>=0; d--) {
							if (gcpft.treeage_avr_20[d].size())
								diam_bin_reserve = d;
							if (diam * 100 > d * DIAM_CLASS_SIZE_SAMP && gcpft.treeage_avr_20[d].size()) {
								diam_bin = d;
								break;
							}
						}
						if (diam_bin > -1) {
							age = (int)gcpft.treeage_avr_20[diam_bin].mean();

						double diam_bin_mean = ((double)diam_bin + 0.5) * DIAM_CLASS_SIZE_SAMP;
						if((diam * 100) < diam_bin_mean) {
							if(gcpft.treeage_avr_20[diam_bin - 1].size() && (diam_bin - 1) >= 0) {
								double diam_bin_mean_prev = ((double)(diam_bin - 1) + 0.5) * DIAM_CLASS_SIZE_SAMP;
								double age_prev = gcpft.treeage_avr_20[diam_bin - 1].mean();
								age -= (int)((diam_bin_mean - (diam * 100)) / (diam_bin_mean - diam_bin_mean_prev) * (age - age_prev));
							}
						}
						else {
							if(gcpft.treeage_avr_20[diam_bin + 1].size() && (diam_bin + 1) < N_DIAM_CLASSES_SAMP) {
								double diam_bin_mean_next = ((double)(diam_bin + 1) + 0.5) * DIAM_CLASS_SIZE_SAMP;
								double age_next = gcpft.treeage_avr_20[diam_bin + 1].mean();
								age -= (int)((diam_bin_mean - (diam * 100)) / (diam_bin_mean - diam_bin_mean_next) * (age - age_next));
							}
						}
							if(print_init_info)
								dprintf("age used from diam bin %d: %d\n", diam_bin, age);
						}
						else if (diam_bin_reserve > -1) {
							age = (int)gcpft.treeage_avr_20[diam_bin_reserve].mean();
							if(print_init_info)
								dprintf("age used from diam bin %d: %d\n", diam_bin_reserve, age);
						}
					}
				}

				double cmass_sap_ind = q_sap * cmass_wood_ind;
				double cmass_heart_ind = cmass_wood_ind - cmass_sap_ind;
				double latosa = indiv.pft.k_latosa;
				if(use_latosa_20y_mean && gcpft.latosa_avr_20.sum() != 0)	// else use the value calculated above
					latosa = gcpft.latosa_avr_20.mean();

				if(print_init_info) {
					dprintf("latosa_avr_20.size for %s: %d \n", (char*)pft.name, gcpft.latosa_avr_20.size());
					dprintf("latosa_avr_20[] for %s: \n", (char*)pft.name);
					for(unsigned int x=0;x<gcpft.cton_sap_avr_20.size();x++) {
						dprintf("%.3f ", gcpft.latosa_avr_20[x]);
					}
				}

				double cmass_leaf_ind = cmass_sap_ind * latosa / (height * indiv.pft.sla * indiv.pft.wooddens);
				double ltor = pft.ltor_max;
				double cmass_root_ind = cmass_leaf_ind / ltor;

				indiv.diam = diam;
				indiv.height = height;
				indiv.densindiv = treedens;
				indiv.cmass_sap = cmass_sap_ind * indiv.densindiv;
				indiv.cmass_heart = cmass_heart_ind * indiv.densindiv;
				indiv.cmass_leaf = cmass_leaf_ind * indiv.densindiv;
				indiv.cmass_root = cmass_root_ind * indiv.densindiv;
				indiv.crownarea = indiv.pft.k_allom1 * pow(indiv.diam, indiv.pft.k_rp);
				indiv.lai = indiv.cmass_leaf * indiv.pft.sla;
				indiv.lai_indiv = indiv.cmass_leaf / indiv.densindiv * indiv.pft.sla / indiv.crownarea;
				indiv.boleht = indiv.height * indiv.pft.boleht_frac;
				indiv.fpc = indiv.crownarea * indiv.densindiv *	(1.0 - exp(-0.5 * indiv.lai_indiv));
				indiv.age = (double)age;

				// NB. If setting to alive, we need to report fluxes here. Put into seed fluxes
				indiv.alive = true;
				patch.fluxes.report_flux(Fluxes::SEEDC, -(indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap + indiv.cmass_heart));

				double cton_leaf = indiv.pft.cton_leaf_avr;
				double cton_root = indiv.pft.cton_root_avr;
				double cton_sap = indiv.pft.cton_sap_avr;
				double cton_heart = indiv.pft.cton_sap_max;

				if(print_init_info) {	
					dprintf("cton_sap_avr_20.size for %s: %d \n", (char*)pft.name, gcpft.cton_sap_avr_20.size());
					dprintf("cton_sap_avr_20[] for %s: \n", (char*)pft.name);
					for(unsigned int x=0;x<gcpft.cton_sap_avr_20.size();x++) {
						dprintf("%.3f ", gcpft.cton_sap_avr_20[x]);
					}
				}

				if (use_cton_20y_mean && gcpft.cton_sap_avr_20.sum() != 0) {
					if(print_init_info) {
						dprintf("\nUsing cton_leaf_avr_20 for %s: \n", (char*)pft.name);
						for(unsigned int x=0;x<gcpft.cton_leaf_avr_20.size();x++) {
							dprintf("%.3f ", gcpft.cton_leaf_avr_20[x]);
						}
						dprintf("\n");
						dprintf("Using cton_sap_avr_20 for %s: \n", (char*)pft.name);
						for(unsigned int x=0;x<gcpft.cton_sap_avr_20.size();x++) {
							dprintf("%.3f ", gcpft.cton_sap_avr_20[x]);
						}
						dprintf("\n");
					}
					cton_leaf = gcpft.cton_leaf_avr_20.mean();
					cton_root = cton_leaf;
					cton_sap = gcpft.cton_sap_avr_20.mean();
					cton_heart = cton_sap;
				}
				if(print_init_info) {
					dprintf("Values used for for %s: \n", (char*)pft.name);
					dprintf("cton_leaf=%.3f\n", cton_leaf);
					dprintf("cton_root=%.3f\n", cton_root);
					dprintf("cton_sap=%.3f\n", cton_sap);
					dprintf("cton_heart=%.3f\n", cton_heart);
					dprintf("\n");
				}

				if(cton_leaf)
					indiv.nmass_leaf = indiv.cmass_leaf / cton_leaf;
				if(cton_root)
					indiv.nmass_root = indiv.cmass_root / cton_root;
				if(cton_sap)
					indiv.nmass_sap = indiv.cmass_sap / cton_sap;
				if(cton_heart)
					indiv.nmass_heart = indiv.cmass_heart / cton_heart;

				// Bypass alive conditional in indiv:report:flux (?)
				patch.fluxes.report_flux(Fluxes::SEEDN, -(indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_sap + indiv.nmass_heart));

				if(print_init_info) {
					dprintf("Creating %s individual id %d\n", (char*)pftname, indiv.id);
					dprintf("latosa=%f cmass=%f f_sap=%f q_sap=%f cton_leaf=%f cton_sap=%f\n", latosa, cmass_sap_ind + cmass_heart_ind + cmass_leaf_ind + cmass_root_ind, f_sap, q_sap, cton_leaf, cton_sap);
					dprintf("cmass_leaf=%f cmass_root=%f cmass_sap=%f cmass_heart=%f\n", cmass_leaf_ind, cmass_root_ind, cmass_sap_ind, cmass_heart_ind);
					dprintf("nmass_leaf=%f nmass_root=%f nmass_sap=%f nmass_heart=%f\n", indiv.nmass_leaf / indiv.densindiv, indiv.nmass_root / indiv.densindiv, indiv.nmass_sap / indiv.densindiv, indiv.nmass_heart / indiv.densindiv);
					dprintf("age=%d diam=%f height=%f densindiv=%f crownarea=%f lai_indiv=%f lai=%f\n\n", age, diam, height, indiv.densindiv, indiv.crownarea, indiv.lai_indiv, indiv.lai);
				}

				// Calculate storage pool size
				indiv.max_n_storage = indiv.cmass_sap * indiv.pft.fnstorage / cton_leaf;
				if(indiv.cmass_sap + indiv.cmass_heart + indiv.cmass_leaf + indiv.cmass_root)
					indiv.scale_n_storage = indiv.max_n_storage * cton_leaf / (indiv.cmass_sap + indiv.cmass_heart + indiv.cmass_leaf + indiv.cmass_root);
			}
		}
		patch.pft[pft.id].nsapling = 0.0;
		patch.pft[pft.id].wscal_mean_est = 0.0;
		patch.pft[pft.id].anetps_ff_est = 0.0;

		pftlist.nextobj();
	}

	double cmass_veg_tree_tot_new = check_harvest_cmass(patch);

	int scale_soil = 0;		// 0 = don't scale soil; 1 = scale all pools; 2 = scale slow pools

	if(scale_soil) {

		double updated_cmass_veg_ratio = 0.0;
		if(cmass_veg_tree_tot_old)
			updated_cmass_veg_ratio = cmass_veg_tree_tot_new / cmass_veg_tree_tot_old;
		for (size_t i = 0; i < patch.soil.nsoillayers; i++) {
			Soillayer &layer = patch.soil.soillayers[i];
			for (int r = 0; r < NSOMPOOL; r++) {

				if(scale_soil == 1) {
					layer.sompool[r].cmass *= updated_cmass_veg_ratio;
					layer.sompool[r].nmass *= updated_cmass_veg_ratio;
				}
				else if(scale_soil == 2) {	// NB: these pools are not included in the SoilC column in cpool.out.
					if (r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
						layer.sompool[r].cmass *= updated_cmass_veg_ratio;	
						layer.sompool[r].nmass *= updated_cmass_veg_ratio;
					}
					else if (r == SURFFWD || r == SURFCWD) {
						layer.sompool[r].cmass *= updated_cmass_veg_ratio;
						layer.sompool[r].nmass *= updated_cmass_veg_ratio;
					}
				}
			}
		}
		pftlist.firstobj();
		while (pftlist.isobj) {

			Pft& pft=pftlist.getobj();
			Patchpft& patchpft = patch.pft[pft.id];
			patchpft.cmass_litter_leaf *= updated_cmass_veg_ratio;
			patchpft.cmass_litter_root *= updated_cmass_veg_ratio;
			patchpft.cmass_litter_sap *= updated_cmass_veg_ratio;
			patchpft.cmass_litter_heart *= updated_cmass_veg_ratio;
			patchpft.cmass_litter_repr *= updated_cmass_veg_ratio;
			patchpft.nmass_litter_leaf *= updated_cmass_veg_ratio;
			patchpft.nmass_litter_root *= updated_cmass_veg_ratio;
			patchpft.nmass_litter_sap *= updated_cmass_veg_ratio;
			patchpft.nmass_litter_heart *= updated_cmass_veg_ratio;

			pftlist.nextobj();
		}
	}
}

/// Function distributing soil C & N equally among all initialised stand types after killing vegetation. LC input must have equal area for all 
/*  these stand types the following year.
*/
void establish_input_setup_stands(Gridcell& gridcell) {

	landcover_change_transfer transfer;
	double receiving_fraction = 0.0;

	// Set required stand variables
	for(unsigned int s=0;s<gridcell.nbr_stands();s++) {
		Stand& stand = gridcell[s];
		Gridcellst& gcst = stand.get_gridcell().st[stand.stid];

		if(gcst.prescribed_establishment) {
			stand.gross_frac_decrease = stand.get_gridcell_fraction();
			stand.transfer_area_st[stand.stid] = stand.gross_frac_decrease;
		}
		else {
			stand.gross_frac_decrease = 0.0;	// To avoid problems from possible LCC for other stands this year
		}
	}

	// Scale by area fractions for this year and next year.
	double init_frac = 0.0;
	double init_frac_nextyear = 0.0;

	for(unsigned int s=0; s<stlist.nobj; s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[st.id];
		if(gcst.prescribed_establishment) {
			init_frac += gcst.frac;
			init_frac_nextyear += gcst.frac_nextyear;
		}
	}

	// Set required st variables
	for(unsigned int s=0; s<stlist.nobj; s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[st.id];

		if (gcst.prescribed_establishment) {
			double st_new_fraction = gcst.frac_nextyear;

			if(init_frac_nextyear)
				st_new_fraction *= (init_frac / init_frac_nextyear);
			else
				fail("In establish_input_setup_stands(): st input file not consistent with forest_init_poolsoil=true. st sum of initiaised stands %d can not have zero area year %d.\n\n", 
					(int)st.landcover, date.get_calendar_year() + 1);
			gcst.frac_old = gcst.frac;
			gcst.frac = st_new_fraction;
			receiving_fraction += st_new_fraction;
		}
	}

	lc_change_harvest_params harvest_params;

	harvest_params.harv_eff = 1.0;
	harvest_params.res_outtake_twig = 1.0;
	harvest_params.res_outtake_coarse_root = 1.0;

	// Kill vegetation in all stands, pool soil, ((change forest area to next year's value)).
	donor_stand_change(gridcell, receiving_fraction, transfer, -1, -1, -1, -1, -1, &harvest_params, true);

	// Kill all stands that have prescribed establishment, create one new stand:
	for(unsigned int s=0; s<stlist.nobj; s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[st.id];

		if (!gcst.prescribed_establishment)
			continue;

		Gridcell::iterator gc_itr = gridcell.begin();
		while (gc_itr != gridcell.end()) {
			Stand& stand_kill = *gc_itr;

			if(stand_kill.stid == st.id) {
				dprintf("Killed st %s stand %d\n", (char*)st.name, stand_kill.id);	
				gc_itr = gridcell.delete_stand(gc_itr);
			}
			else {
				++gc_itr;
			}
		}
		if(gcst.frac) {
			Stand& new_stand = gridcell.create_stand_lu(st, gcst.frac, npatch_treeinit, false);
			receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand.id);
			dprintf("%s stand no.%d created\n", (char*)st.name, new_stand.id);
		}
	}


}

void establish_input_gridcell(Gridcell& gridcell) {

	if (!readforest_init || !(date.islastday && date.islastmonth))
		return;

	int establishment_year = -1;
	bool establish_thisyear = false;

	// Check that establishment year is the same for all stand types with option forest_init_poolsoil
	for(unsigned int s=0; s<stlist.nobj; s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[st.id];

		int establishment_year_temp = gcst.establishment_list.establishment_year;
		if(gcst.prescribed_establishment && establishment_year_temp != -1) {
			if(forest_init_poolsoil && establishment_year != -1 && establishment_year != establishment_year_temp)
				fail("With option forest_init_poolsoil, all stand types must have the same establishement year !\n\n");
			establishment_year = establishment_year_temp;
			if(establishment_year == date.get_calendar_year())
				establish_thisyear = true;
		}
	}

	if(establish_thisyear)
		gridcell.year_treeinit = date.get_calendar_year();
	else
		return;			

	if(forest_init_poolsoil)
		establish_input_setup_stands(gridcell);

	for(unsigned int s=0;s<gridcell.nbr_stands();s++) {
		Stand& stand = gridcell[s];
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			establish_input(patch, forest_init_poolsoil);
			stand.nextobj();
		}
	}

	bool use_rotation = true;	// Keeping the option without rotation for now, but remove later if confident in functionality
	// Turned on again 240708. The solution is that stand types must have separate secondary management types.

	for(unsigned int s=0; s<stlist.nobj; s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[st.id];

		if(gcst.prescribed_establishment) {

			ManagementType& mt = st.get_management(1);
//			ManagementType& mt = st.get_management();		// if(!use_rotation) : set manually

			char selection_cp[200] = {0};

			for(int p=0; p<npft; p++) {
				Pft& pft = pftlist[p];
				for (int i = 0;i < MAX_NESTABLISH_INDIV;i++) {
					xtring pftname = gcst.establishment_list.est_array[i].pft_name;
					if(!strcmp(pft.name, (char*)pftname)) {
						strcat(selection_cp, pftname);
						strcat(selection_cp, " ");
						break;
					}
				}
			}
			mt.selection = selection_cp;
			if(!use_rotation) {
				mt.relaxed_establishment = true;
				mt.cutfirstyear = 0;
				mt.harvest_system = "";
			}
		}
	}

	// Redo management setup for stands after updating parameters
	for(unsigned int s=0;s<gridcell.nbr_stands();s++) {
		Stand& stand = gridcell[s];
		Gridcellst& gcst = stand.get_gridcell().st[stand.stid];
		if(gcst.prescribed_establishment) {
			if(!use_rotation)
				stand.set_management();
			else
				// Use second rotation (set rotation year to far into the future), also calls set_management()
				stand.rotate(1);
		}
	}
}

void force_prescribed_disturbance(Gridcell& gridcell, Rank_stands& stand_class, disturbancetype disttype) {

	double disturbance_fraction = 0.0;

	lc_change_harvest_params harvest_params;

	// Default forestry harvest parameters: pft values
	harvest_params.harv_eff = 0.9;
	harvest_params.res_outtake_twig = 0.4;
	harvest_params.res_outtake_coarse_root = 0.1;

	if(disttype == CLEARCUT_PRESCRIBED) {
		disturbance_fraction = gridcell.dist_clearcut_gc;
		stand_class.sort_diameter_rms();
	}
	else if(disttype == FIRE) {
		disturbance_fraction = gridcell.dist_fire_gc;
		stand_class.sort_fire_sensitivity();
	}
	else if(disttype == WIND) {

		disturbance_fraction = gridcell.dist_wind_gc;
		stand_class.sort_storm_sensitivity();

        // Currently unchanged from above - i.e. assume a salvage logging with the same intensity as typically used for wood harvest
		harvest_params.harv_eff = 0.9;
		harvest_params.res_outtake_twig = 0.4;
		harvest_params.res_outtake_coarse_root = 0.1;
	}
	else if(disttype == BARKBEETLE) {
		disturbance_fraction = gridcell.dist_barkbeetle_gc;
		stand_class.sort_barkbeetle_sensitivity();

        // Currently unchanged from above - i.e. assume a salvage logging with the same intensity as typically used for wood harvest
		harvest_params.harv_eff = 0.9;
		harvest_params.res_outtake_twig = 0.4;
		harvest_params.res_outtake_coarse_root = 0.1;
	}

	if(!disturbance_fraction)
		return;
//	else
//		dprintf("\nYear %d: Disturbance type %d gridcell fraction = %f\n\n", date.get_calendar_year(), (int)disttype, disturbance_fraction);

	// With the assumption that NATURAL is forest, might be a problem with a natural grassland
	double disturbance_fraction_remain = disturbance_fraction * (gridcell.landcover.frac[FOREST] + gridcell.landcover.frac[NATURAL]);
	
	for(unsigned int i=0;  i<stand_class.nstands_initial && disturbance_fraction_remain > 0.0; i++) {

		int stand_index = stand_class.get_index(stand_class.nstands_initial -1 - i);
		Stand& stand = gridcell[stand_index];
		int stand_age = (date.year - max(stand.first_year, stand.clone_year));

		// Avoid disturbing the the same stand twice
		if(!stand.disturbed && (stand.landcover == FOREST || stand.landcover == NATURAL)) {

//			dprintf("Area fraction of stand %d = %f, vegetation C = %f, age = %d\n\n", stand.id, stand.get_gridcell_fraction(), (stand.ccont() - stand.ccont(0.0)), stand_age);
			ManagementType& mt = stand.get_current_management();

			// Kill vegetation for the whole stand
			if(disturbance_fraction_remain >= stand.get_gridcell_fraction()) {

//				dprintf("Year %d: Disturbance type %d: st %s, Stand %d (area fraction %f) harvested/disturbed\n", date.get_calendar_year(), (int)disttype, (char*)stlist[stand.stid].name, stand.id, stand.get_gridcell_fraction());

				stand.firstobj();

				while (stand.isobj) {

					Patch& patch = stand.getobj();

					if(!patch.age) {
						stand.nextobj();
						continue;
					}
					patch.man_strength = 1.0;

					// blaze function at patch level
					if(disttype == FIRE) {
						if(firemodel != BLAZE)
							fail("fire disturbance input only implemented for BLAZE fore model\n");
						blaze(patch, gridcell.climate, true);
					}

					// Non-fire prescribed disturbance is implemented as a
					// full-strength harvest/disturbance treatment, not by
					// calling the generic patch-destroying disturbance().
					// Trees are removed by harvest_forest(); grasses are
					// only removed when the management type requests this.
					Vegetation& vegetation = patch.vegetation;
					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv = vegetation.getobj();

						if(disttype != FIRE) {
							bool killed = false;
							harvest_forest(indiv, indiv.pft, indiv.alive, 0.0, killed, &harvest_params, disttype);
							if(killed) {
//								dprintf("Year %d: ST %d, Stand %d, Patch %d; Individual %d killed in disturbance no. %d\n", date.get_calendar_year(), stand.stid, stand.id, patch, indiv.id, disttype);
							}
							else {
								vegetation.nextobj();	// Grass not killed if mt.killgrass_at_cc is false!
							}
						}
						else {
							vegetation.nextobj();
						}
					}
					patch.disturbed = true;
					patch.age = 1;				// Needs to be set to 1 to behave like clearcut in harvest_forest day 0, since this is after vegetation_dynamics()
					// Try these setting, perhaps they are needed:
					patch.plant_this_year = true;
					patch.prescribed_disturbed = true;	// Will delay resetting of plant_this_year
					stand.nextobj();
				}
				disturbance_fraction_remain -= stand.get_gridcell_fraction();
				stand.disturbed = true;
			}
			else {

				// Use landcover_dynamics functionality here to kill part of the stand's vegetation (area)

				double stand_frac = stand.get_gridcell_fraction();

				// From reduce_stands():
				stand.frac_change -= disturbance_fraction_remain;
				stand.gross_frac_decrease += disturbance_fraction_remain;
				stand.transfer_area_st[stand.stid] += disturbance_fraction_remain;
				stand.set_gridcell_fraction(stand_frac - disturbance_fraction_remain);

				// new stands may have several patches; if mt.cutfirstear == true, they may be collapsed to a single patch-stands below.
				Stand& new_stand = stand.clone(stlist[stand.stid], disturbance_fraction_remain, false, true);

				new_stand.cloned = true;
				new_stand.cloned_fraction = disturbance_fraction_remain;
				new_stand.lc_origin = stand.landcover;
				new_stand.st_origin = stand.stid;
				new_stand.frac_temp = disturbance_fraction_remain;	// Needed in stand patch reduction below

				stand.frac_temp -= disturbance_fraction_remain;

//				dprintf("Year %d: Disturbance type %d: st %s, stand %d, Fraction of %f put into new stand no.%d\n", date.get_calendar_year(), (int)disttype, (char*)stlist[stand.stid].name, stand.id, disturbance_fraction_remain, new_stand.id);

				new_stand.firstobj();

				while (new_stand.isobj) {

					Patch& patch = new_stand.getobj();

					patch.man_strength = 1.0;
					patch.plant_this_year = true;
					patch.prescribed_disturbed = true;	// Will delay resetting of plant_this_year

					if(disttype == FIRE) {		// blaze function at patch level
						if(firemodel != BLAZE)
							fail("fire disturbance input only implemented for BLAZE fore model\n");
						blaze(patch, gridcell.climate, true);	// What happens to grass ?
					}
					else {
						// This cloned stand represents the disturbed
						// fraction of the original stand. As above,
						// non-fire prescribed disturbance is a
						// full-strength harvest/disturbance treatment
						// rather than generic patch destruction.
						Vegetation& vegetation = patch.vegetation;
						vegetation.firstobj();
						while (vegetation.isobj) {
							Individual& indiv = vegetation.getobj();
							Patchpft& ppft = patch.pft[indiv.pft.id];

							if(disttype != FIRE) {
								// Alt.1. use for clearcut:
								bool killed = false;
								harvest_forest(indiv, indiv.pft, indiv.alive, 0.0, killed, &harvest_params, disttype);	// use anpp ?
								if(killed) {
//									dprintf("Year %d: ST %d, Stand %d, Individual %d killed in disturbance no. %d\n", date.get_calendar_year(), stand.stid, stand.id, indiv.id, disttype);
								}
								else {
									vegetation.nextobj();	// Grass not killed if mt.killgrass_at_cc is false!
								}
							}
							else {
								dprintf("Add code for new disturbance type in disturb_force()\n");
								vegetation.nextobj();
							}
						}
						stand.nextobj();
					}
					patch.disturbed = true;
					patch.age = 1;				// Needs to be set to 1 to behave like clearcut in harvest_forest day 0, since this is after vegetation_dynamics()

					new_stand.nextobj();
				}

				// We can only reduce the number of patches in stands where all vegetation is killed in normal forestry wood harvest (
				if(mt.killgrass_at_cc && stand.nobj > 1) {

					landcover_change_transfer transfer;
					dprintf("Converting stand to single-patch stand\n");

					// 1. Copy soil biomass and water from a stand with several patches to a new stand with a single patch after killing remaining grass:
					new_stand.transfer_area_st[new_stand.stid] += disturbance_fraction_remain;
					donor_stand_change(gridcell, disturbance_fraction_remain, transfer, -1,-1, new_stand.stid, -1, new_stand.id, &harvest_params);
					Stand& new_stand2 = gridcell.create_stand_lu(stlist[new_stand.stid], disturbance_fraction_remain, npatch_disturbedstand, false, false, stand.current_rot);
					receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand2.id);
					dprintf("Single-patch stand created\n");

					while(new_stand2.isobj) {

						Patch& patch = new_stand2.getobj();
						patch.disturbed = true;
						patch.plant_this_year = true;
						patch.prescribed_disturbed = true;	// Will delay resetting of plant_this_year
						new_stand2.nextobj();
					}

					// 2. Remove new stand with many patches
					int killed_stand_id = new_stand.id;

					Gridcell::iterator gc_itr = gridcell.begin();
					while (gc_itr != gridcell.end()) {
						Stand& stand_kill = *gc_itr;
						if(stand_kill.id == new_stand.id) {
							gc_itr = gridcell.delete_stand(gc_itr);
						}
						else {
							++gc_itr;
						}
					}
					dprintf("Killed temporary stand %d\n\n", killed_stand_id);
				}
				stand.gross_frac_decrease -= disturbance_fraction_remain;
				disturbance_fraction_remain = 0.0;
			}
		}
	}
}

void prescribed_disturbance(Gridcell& gridcell) {

	// Prescribed disturbance input allocates externally supplied disturbed
	// area fractions to ranked stands, then applies disturbance-specific
	// harvest or fire handling. This is not the same as the generic
	// patch-destroying disturbance() pathway: non-fire prescribed
	// disturbances remove trees via harvest_forest(), may leave grass, and
	// set patch.age to 1.
	if (!(date.islastday && date.islastmonth)
		|| !(readdist_clearcut || readdist_fire || readdist_wind || readdist_barkbeetle) 
		|| (date.get_calendar_year() < disturbance_input_firstyear || date.get_calendar_year() > disturbance_input_lastyear))
		return;
	
	// Sort stands
	Rank_stands stand_class(gridcell);

	// Loop through the different disturbance types in the order specified in the disturbancetype enumeration
	for(int i=0; i<NDISTURBANCETYPES;i++) {
		force_prescribed_disturbance(gridcell, stand_class, (disturbancetype)i);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// VEGETATION DYNAMICS
// Should be called by framework at the end of each simulation year, after vegetation,
// climate and soil attributes have been updated

void vegetation_dynamics(Stand& stand,Patch& patch) {

	// DESCRIPTION
	// Implementation of fire disturbance and population dynamics (establishment and
	// mortality) at end of simulation year. Bioclimatic constraints to survival and
	// establishment are imposed within mortality and establishment functions
	// respectively.

	double fireprob = 0.0;
		// probability of fire in this patch this year
		// (in population mode: fraction of modelled area affected by fire this year)

	// Calculate fire probability and volatilise litter
	if (patch.has_fires() && firemodel == GLOBFIRM) {
		fire(patch,fireprob);
	}
	patch.fireprob = fireprob;

	if (vegmode == POPULATION) {

		// POPULATION MODE

		// Mortality
		mortality_lpj(stand, patch, stand.get_climate(), fireprob);

		// Establishment
		establishment_lpj(stand,patch);

	}
	else {

		// INDIVIDUAL AND COHORT MODES

		// Patch-destroying disturbance

		// Disturbance when N limitation is switched on to get right pft composition under N limitation faster
		if (ifcentury && ifnlim && (date.year == patch.soil.solvesomcent_beginyr - 2 || date.year == patch.soil.solvesomcent_endyr - 2)){
			disturbance(patch, 1.0);
			if (patch.disturbed) {
				insectupdate(stand, patch); // Simplified updatde of insect dynamics normally done in mortality_guess
				return; // no mortality or establishment this year
			}
		}

		// Generic patch-destroying disturbance. This includes the Dave
		// deterministic mode, which is a fixed-year wrapper around
		// disturbance(patch, 1.0).
		//
		// During active prescribed-disturbance input years, this pathway is
		// suppressed because prescribed_disturbance() has already allocated
		// the externally supplied disturbed area for the gridcell. Keeping
		// the two pathways mutually exclusive avoids disturbing the same
		// patch through two different disturbance concepts in one year.
		bool disturbance_model = !(readdist_clearcut || readdist_fire || readdist_wind || readdist_barkbeetle)
			|| (date.get_calendar_year() < disturbance_input_firstyear || date.get_calendar_year() > disturbance_input_lastyear);

		if (patch.has_disturbances() && disturbance_model) {

			// We don't allow disturbance while documenting for calculation of Century equilibrium
			bool during_century_solvesom = ifcentury &&
										   date.year >= patch.soil.solvesomcent_beginyr &&
										   date.year <= patch.soil.solvesomcent_endyr;
			if (!during_century_solvesom && patch.age) {
				if (disttype == DETERMINISTIC) {
					// Dave deterministic disturbance is a debugging and
					// reproducibility mode: predefined years force the
					// generic patch-destroying disturbance with
					// probability 1.0.
					const Gridcell &gridcell = stand.get_gridcell();
					const std::set<int> &years = gridcell.static_disturbances;
					if (years.find(date.get_calendar_year()) != years.end()) {
						// Deterministic disturbance occurs in this year.
						disturbance(patch, 1.0);
					}
					else {
						patch.disturbed = false;
					}
				} else {
					// Stochastic disturbances.
					double distinterval_ = stand.get_distinterval();	
					disturbance(patch,1.0 / distinterval_);
				}
			}
		}

		if (patch.disturbed) {
			insectupdate(stand, patch); // Simplified updatde of insect dynamics normally done in mortality_guess
			return; // no mortality or establishment this year
		}

		// Mortality
		mortality_guess(stand, patch, stand.get_climate(), fireprob);
		insect_mortality_guess(stand, patch);
		storm_mortality_guess(stand, patch);

		bool establishment_input = stand.get_gridcell().st[stand.stid].establishment_list.establishment_year == date.get_calendar_year();

		// Establishment
		if (!establishment_input) {
			establishment_guess(stand, patch);
		}
	}

	patch.age++;

	for (unsigned int i=0; i<patch.vegetation.nobj; i++) {

		Individual& indiv = patch.vegetation[i];
		if(indiv.pft.lifeform == TREE) {
			indiv.cmass_wood_inc_5.add(indiv.cmass_wood() - indiv.cmass_wood_old);
			indiv.cmass_wood_old = indiv.cmass_wood();
		}
	}
}

void vegetation_dynamics_daily(Patch& patch, Climate& climate) {
	if (ifdavegrowth){
		establishment_daily(patch, climate);
		mortality_daily(patch); // TODO Assaf
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2001
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Pacala SW, Canham, CD & Silander JA Jr 1993 Forest models defined by field
//   measurements: I. The design of a northeastern forest simulator. Canadian Journal
//   of Forest Research 23: 1980-1988.
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65:
//   51-70.
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Smith, B, Prentice, IC & Sykes, M (2001) Representation of vegetation dynamics in
//   the modelling of terrestrial ecosystems: comparing two contrasting approaches
//   within European climate space. Global Ecology and Biogeography 10: 621-637
// Thonicke, K, Venevsky, S, Sitch, S & Cramer, W (2001) The role of fire disturbance
//   for global vegetation dynamics: coupling fire into a Dynamic Global Vegetation
//   Model. Global Ecology and Biogeography 10: 661-677.
// Wolf, A., Callaghan T.V., & Larson K. (2008) Future changes in vegetation and ecosystem 
//   function of the Barents Region. Climatic Change, 87:51-73 DOI 10.1007/s10584-007-9342-
// Zwillinger, D 1996 CRC Standard Mathematical Tables and Formulae, 30th ed. CRC
//   Press, Boca Raton, Florida.
