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


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Maximum number of equivalent days with full leaf cover per growing season
// for summergreen PFTs
const double APHEN_MAX = 210.0;

// ---------------------------------
// Grass parameters
// ---------------------------------
// Maximum size of C storage pool as a fraction of total biomass (leaf + root).
const double frac_cstorage_max = 0.3;

// Fraction of C storage not available for growth expressed on a total biomass (leaf + root) basis.
const double frac_cstorage_res = 0.05;

// Maximum mobilisation rate of C in storage
const double max_mobil_cstorage = 0.2;

// Dormancy temperature threshold C3 grasses (degC)
const double temp_dormancy_C3 = 5.0;

// Dormancy temperature threshold C4 grasses (degC)
const double temp_dormancy_C4 = 15.0;

// Minimum physiological activity during dormancy (fraction)
const double fphen_activity_min = 0.1;

// Maximum physiological up/down-regulation in and out of dormancy (fraction)
const double dphen_activity_max = 0.3;

// Fraction of reproductive pool available for reproduction
const double fracseed = 0.5;

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
	return indiv.pft.lifeform == GRASS;
}

double get_nmass(const Individual &indiv) {
	return indiv.nmass_storage +
		   indiv.nmass_leafs[green] +
		   indiv.nmass_leafs[brown] +
		   indiv.nmass_root +
		   indiv.nmass_repr +
		   indiv.nmass_litter_leaf +
		   indiv.nmass_litter_root;
}

double get_cmass(const Individual &indiv) {
	return indiv.cmass_storage +
		   indiv.cmass_leafs[green] +
		   indiv.cmass_leafs[brown] +
		   indiv.cmass_root +
		   indiv.cmass_repr +
		   indiv.cmass_litter_leaf +
		   indiv.cmass_litter_repr +
		   indiv.cmass_litter_root;
}

DaveMassBalance::IndivState::IndivState(const Individual &indiv) {
	n_start = get_nmass(indiv);
	c_start = get_cmass(indiv);
	id = indiv.id;
}

void DaveMassBalance::IndivState::check(const Individual &indiv) const {
	double c_flux = indiv.dnpp;
	double c_end = get_cmass(indiv) + indiv.cexcess;
	assert(equal(c_end, c_start + c_flux));

	double n_flux = indiv.nmass_missing;
	double n_end = get_nmass(indiv);
	assert(equal(n_end, n_start + n_flux));
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

void flush_repr_litter_DAVE(Individual &indiv);

///////////////////////////////////////////////////////////////////////////////////////
// LEAF PHENOLOGY
// Call function leaf_phenology each simulation day prior to calculation of FPAR, to
// calculate fractional leaf-out for each PFT and individual.
// Function leaf_phenology_pft is not intended to be called directly by the framework,

// Water stress scalar used for growth and dormancy
double wscal_growth(double wscal) {
	double qgrowth = 1.1;
	return min(1.0,max(0.0,pow(wscal,qgrowth)));
}

// Temperature stress scalar used for growth and dormancy
double temp_growth(double temp) {
	return Arrhenius_temp_response(1.78,temp,76.8,285.0,933)/2.23;  // Note: parameterise properly (now: Parent et al. 2010)
}

		

void phenology_grass_dave(Vegetation& vegetation, Climate& climate) {

	// Calculates phenological cycle of grasses. 
	// Returns phenological state of grasses as factor fphen [0-1].
	// Implemented for grasses only. 
    for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		Pft &pft = pftlist.getobj();
		Patchpft &ppft = vegetation.patch.pft[pft.id];
		ppft.airtemp_10.add(climate.temp);
	}
	
	
	// Loop through cohorts (=individuals)
	vegetation.firstobj();
	while (vegetation.isobj) {

		// Retrieve information about cohort
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft &ppft = vegetation.patch.pft[indiv.pft.id];

		if (pft.lifeform == GRASS) { 

			// Threshold temperature for dormancy (base temperature same as dormancy threshold, depends on photosyn pathway)
			double temp_dormancy = indiv.pft.pathway == C3 ? temp_dormancy_C3 : temp_dormancy_C4;
			double grass_temp = indiv.vegetation.patch.soil.T_soil[6]; // Air temperature or bottom snow layer if snow present
			double tscal_crit = temp_growth(temp_dormancy);
			double wscal_mean = ppft.wscal_10.size() > 0 ? ppft.wscal_10.mean() : 1.0;

			double wscal_actual = max(pft.wscal_min, wscal_growth(wscal_mean));
			double tscal_actual = max(tscal_crit, temp_growth(ppft.airtemp_10.mean()));

			// Changes in physiological activity due to changes in dormancy status
			double tscal = (tscal_actual - tscal_crit) / (1.0 - tscal_crit);
			double wscal = (wscal_actual - pft.wscal_min) / (1.0 - pft.wscal_min);
			assert(tscal >= 0);
			assert(wscal >= 0);

			double dphen_activity = min(tscal,wscal)*dphen_activity_max;

			indiv.fphen_activity = min(1.0,max(fphen_activity_min,indiv.fphen_activity * (1.0 + dphen_activity)));

			// Determine if grasses enter dormancy
			//if ( (ppft.wscal_10.size() > 0 && ppft.wscal_10.mean() < pft.wscal_min) || 
			//	 (ppft.airtemp_10.size() > 0 && ppft.airtemp_10.mean() < temp_dormancy) ) { 
				
			//	// physiological downregulation factor due to dormancy (fphen_activity_min = full dormancy, 1 = no dormancy)
			//	indiv.fphen_activity = max(fphen_activity_min,indiv.fphen_activity * (1.0 + dphen_activity));
			//}
			//else {

				// gradually release dormancy if conditions allow
			//	indiv.fphen_activity = min(1.0,indiv.fphen_activity * (1.0 + dphen_activity));
			//}	

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

	bool raingreen = pft.phenology == RAINGREEN || pft.phenology == ANY;
	bool summergreen = pft.phenology == SUMMERGREEN || pft.phenology == ANY;

	phen = 1.0;

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
		else if (pft.lifeform == GRASS || pft.lifeform == MOSS) {

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
		indiv.phen = patch.pft[indiv.pft.id].phen;

		// Update annual leaf-day sum (raingreen PFTs)
		if (date.day == 0) indiv.aphen_raingreen = 0;
		indiv.aphen_raingreen += (indiv.phen != 0.0);

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

	double dval;
	double cmass_leaf_ind;
	double cmass_root_ind;
	double cmass_sap_ind;

	if (indiv.pft.lifeform == GRASS && ifdailygrass) {
		 indiv.cmass_storage = bminit * indiv.densindiv;
		 indiv.nmass_storage = indiv.cmass_storage / indiv.pft.cton_leaf_avr;

		if (ifdavephotosyn){
			// allocate cmass and nmass in storage to leaves and roots
			// Initial C mass: set C storage to capacity, rest to leaves and roots.
			double cmass_init = indiv.cmass_storage;

			// following 3 equations solved for cleaf:
			// cmass_init = cleaf + croot + cstorage
			// croot = cleaf / ltor
			// cstorage = (cleaf + croot) * frac_cstorage_max
			double cmass_leaf = cmass_init / (1.0 + 1.0 / ltor + frac_cstorage_max + frac_cstorage_max / ltor );
			if (ifdailygrass) {
				indiv.cmass_leafs[green] = cmass_leaf;
			} else {
				indiv.cmass_leaf = cmass_leaf;
			}
			indiv.cmass_root = cmass_leaf / ltor;
			indiv.cmass_storage = cmass_init - cmass_leaf - indiv.cmass_root;

			// Initial N mass
			double nmass_leaf = cmass_leaf / indiv.pft.cton_leaf_avr;
			if (ifdailygrass) {
				indiv.nmass_leafs[green] = nmass_leaf;
			} else {
				indiv.nmass_leaf = nmass_leaf;
				assert(indiv.nmass_leaf >= 0);
			}
			indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_avr;
			indiv.nmass_storage = indiv.cmass_storage / indiv.pft.cton_leaf_avr;
		}

		indiv.alive = true;
		indiv.initialn = indiv.nmass_storage + indiv.get_nmass_leaf() + indiv.nmass_root;
		indiv.report_flux(Fluxes::ESTC, -(indiv.cmass_storage + indiv.get_cmass_leaf() + indiv.cmass_root) );
		indiv.report_flux(Fluxes::ESTN, -indiv.initialn);
	}
	else {

		allocation(bminit, indiv, ltor, cmass_leaf_ind, cmass_root_ind, cmass_sap_ind, dval, dval, dval, dval, dval);

		indiv.cmass_leaf = cmass_leaf_ind * indiv.densindiv;
		indiv.cmass_root = cmass_root_ind * indiv.densindiv;

		if (ifdavephotosyn){
			// we need initial values for nmass_leaf when calling ifdavephotosyn
			indiv.nmass_leaf = indiv.cmass_leaf / indiv.pft.cton_leaf_avr;
			assert(indiv.nmass_leaf >= 0);
			indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_avr;
		}
		else {
			indiv.nmass_leaf = 0.0;
			indiv.nmass_root = 0.0;
		}

		if (indiv.pft.lifeform == TREE) {
			indiv.cmass_sap = cmass_sap_ind * indiv.densindiv;
			indiv.nmass_sap = indiv.cmass_sap / indiv.pft.cton_sap_avr;
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
		}
		else {
			indiv.lai_indiv = 0.0;
			indiv.fpc = 0.0;
		}

		// Bole height (Eqn 7)
		indiv.boleht = 0.0;

		// Stand-level LAI
		indiv.lai = indiv.cmass_leaf * indiv.sla();
	}
	else if (indiv.pft.lifeform == GRASS || indiv.pft.lifeform == MOSS) {

		// GRASSES AND MOSSES

		if (indiv.pft.landcover != CROPLAND) {

			// guess2008 - bugfix - added if
			if (!negligible(indiv.get_cmass_leaf())) {

				// Grass "individual" LAI (Eqn 11). Green LAI only!
				indiv.lai_indiv = indiv.get_cmass_leaf() * indiv.sla();

				// FPC (Eqn 10)
				indiv.fpc = 1.0 - lambertbeer(indiv.lai_indiv);

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
	else if (indiv.pft.lifeform == GRASS || indiv.pft.lifeform == MOSS) {		// grass and moss

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
		Patchpft& pft = patch.pft.getobj();

		// Updated soil fluxes
		patch.fluxes.report_flux(Fluxes::REPRC, pft.cmass_litter_repr);
		pft.cmass_litter_repr = 0.0;

		patch.pft.nextobj();
	}
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		patch.vegetation[i].cmass_litter_repr = 0;
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

		bool killed = false;

		// Standard yearly growth functionalities

		if (!(ifdailygrass && indiv.pft.lifeform == GRASS) || indiv.istruecrop_or_intercropgrass()) {
	
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

						if (ifdavephotosyn) {
							// Set annual N demand when running in newphotosyn mode.
							// fixme: this is just a placeholder, we need a way to run trees on
							// an annual timestep. For now, just request enough N to meet the
							// expected turnover this year.
							indiv.andemand_leaf = indiv.get_nmass_leaf() * indiv.pft.turnover_leaf;
							indiv.andemand_root = indiv.nmass_root * indiv.pft.turnover_root;
							indiv.andemand_sap = indiv.nmass_sap * indiv.pft.turnover_sap;
							indiv.andemand = indiv.andemand_leaf + indiv.andemand_root + indiv.andemand_sap;
							indiv.ndemand = indiv.andemand / 365; // fixme!! (Assaf)
							total_ndemand += indiv.ndemand;

							indiv.leaffndemand = divide(indiv.andemand_leaf, indiv.andemand);
							indiv.rootfndemand = divide(indiv.andemand_root, indiv.andemand);
							indiv.sapfndemand = divide(indiv.andemand_sap, indiv.andemand);

							assert(indiv.leaffndemand >= 0);
							assert(indiv.rootfndemand >= 0);
							assert(indiv.sapfndemand >= 0);

							indiv.leafndemand = indiv.ndemand * indiv.leafndemand;
							indiv.rootndemand = indiv.ndemand * indiv.rootfndemand;
							indiv.sapndemand = indiv.ndemand * indiv.sapfndemand;
							indiv.storendemand = 0;
						}

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

							indiv.kill();

							vegetation.killobj();
							killed = true;
						}
					}
					else if (indiv.pft.lifeform == GRASS || indiv.pft.lifeform == MOSS) {

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
		else if (ifdailygrass && indiv.pft.lifeform == GRASS) {
			
			// Check daily grass in the end of the year after all processes has been completed

			if (negligible(indiv.densindiv)) {
				fail("growth: negligible densindiv for %s", (char*)indiv.pft.name);
			}
			else {

				// Kill individual and transfer biomass to litter if either biomass
				// compartment negative or if yearly growth this year has been low

				if (indiv.cmass_leaf_ymax < MINCMASS || indiv.cmass_root_ymax < MINCMASS || indiv.cmass_leaf_ygrowth < MINCMASS) {

					if (ifdavegrowth) {
						flush_repr_litter_DAVE(indiv);
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

	if (ifdavephotosyn) {
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual &indiv = vegetation.getobj();
			if (indiv.pft.lifeform == TREE) {
				indiv.fnuptake = indiv.ndemand / total_ndemand;
			}
			vegetation.nextobj();
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

		if (indiv.pft.lifeform == GRASS && !indiv.istruecrop_or_intercropgrass()) {
			
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
// C storage unavailable for growth (kg C m-2)
double cstorage_reserve(Individual& indiv){

	return frac_cstorage_res * (indiv.get_cmass_leaf() + indiv.cmass_root);

} 
// Determine growth demand (aka 'sink strength') 
// Returns fraction of available C used for growth vs. storage.
// Put fraction into growth and 1-fraction into storage.
// --> 1.0: put all available C of the day into growth.
// --> 0.0: put all available C of the day into storage.
double growth_demand(Individual& indiv) {

	Vegetation& veg = indiv.vegetation;
	Pft& pft = indiv.pft;
    Patchpft &ppft = veg.patch.pft[indiv.pft.id]; 

	double growth_rate = 0.0;
	double growth_rate_pot = 0.2;  // TODO: parameterise from data. Different between annuals and perennials
	double growth_demand = 0.0;
	double fwater_growth = 0.0;
	double ftemp_growth = 0.0;
	double fdev_growth = 0.0;
	double qgrowth = 1.1;
	double ag_cmass = 0.0;
	bool simple_formulation = false;

	// Test: do we need to account for growth_rate_pot to get different behaviour of annuals vs perennials
	//       or can we get that with SLA alone?

	if (!simple_formulation) {
		// Air temperature or bottom snow layer if snow present
		double grass_temp = veg.patch.soil.T_soil[6];

		// Factors determining growth demand
		fwater_growth = wscal_growth(ppft.wscal);
		ftemp_growth = temp_growth(grass_temp);
		
		ag_cmass = (indiv.get_cmass_leaf() + indiv.cmass_repr);
		ag_cmass = ag_cmass > 0.0 ? ag_cmass : 0.0001; 
		//fdev_growth = ag_cmass/pow(ag_cmass,1.5);
		fdev_growth = divide(pow(ag_cmass,0.75),ag_cmass);

        // growth demand (=sink activity) as a fraction of existing biomass
		// if greater than one, plant is sink limited and growth demand is set to 1 so that growth = NPP + max_mob_rate * storage
		growth_demand = min(1.0,growth_rate_pot * fdev_growth * ftemp_growth * fwater_growth);

		// TODO: make reproduction flux dependent on fdev_growth

	} else {

		// simpler alternative:
		// assumes that growth follows the same environmental responses than NPP
		//if (negligible(indiv.cgrow_leaf)) {
			//growth_rate = min(0.1,growth_rate_pot * (indiv.fphen_activity - fphen_activity_max) * (1.0/(1.0 - fphen_activity_max)));	
		// } else {
       		//growth_rate = max(0.0,indiv.dnpp) / pow((indiv.get_cmass_leaf() + indiv.cmass_repr),1.5);
			//growth_demand = min(1.0, growth_rate / growth_rate_pot);
		//}
		// TODO: does storage itself also respires?
	}

	if (indiv.fphen_activity < 1.0){
		growth_demand *= (indiv.fphen_activity - fphen_activity_min) / (1.0 - fphen_activity_min);
	}

	return growth_demand;

}


// Capacity of C storage (kg C m-2)
double cstorage_capacity(Individual& indiv) {

	double cstorage_new = frac_cstorage_max * (indiv.get_cmass_leaf() + indiv.cmass_root);
    
	if (cstorage_new > indiv.cstorage_max) {
		return cstorage_new;
	} else {
		return indiv.cstorage_max;
	}
	
}



// Capacity of N storage (kg N m-2)
double nstorage_capacity(Individual& indiv) {

	return indiv.cstorage_max / indiv.pft.cton_leaf_min;  // TODO: OK to use leaf C:N here? Different limit for storage?

}

// Flush reproductive carbon and nitrogen to ppft-level reproductive litter (=seedbank).
// Happens every year at the onset of dormancy
// TODO: deal with case of no dormancy
void flush_repr_litter_DAVE(Individual& indiv) {

	indiv.patchpft().cmass_litter_repr += indiv.cmass_repr * fracseed; // TODO: turnover!!
	indiv.patchpft().nmass_litter_repr += indiv.nmass_repr * fracseed; // TODO: turnover!!
	indiv.patchpft().cmass_litter_leaf += indiv.cmass_repr * (1.0 - fracseed);
	indiv.patchpft().nmass_litter_leaf += indiv.nmass_repr * (1.0 - fracseed);
	indiv.cmass_repr = 0.0;
	indiv.nmass_repr = 0.0;
}

void kill_grass(Individual& indiv) {

	if (ifdavegrowth) {
		flush_repr_litter_DAVE(indiv);
	}

	indiv.cmass_storage = 0.0;
	indiv.kill();  // flushes all pools to litter (except for reproductive C, which is handled above)
	indiv.vegetation.killobj();  // removes individual object
}

// Determine how much of the available C (NPP + C fraction from storage) goes to growth vs. storage.
// Function calculates cgrow (daily growth), cstore (change in cmass_storage) 
// and updates indiv.cmass_storage. 
// Note: cgrow and cstore are fluxes (kgC/m2/d), cmass_storage is a pool (kg/m2).
void storage_growth_partitioning(Vegetation& vegetation) {

	// Loop through individuals
	// Note - we sometimes kill individuals inside this loop, but vegetation.killobj()
	// requires the internal pointer to be set appropriately, which effectively
	// prevents us from using a range-based loop, and instead forces us to use
	// the clunky mess you see below. todo: refactor the container classes!
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (indiv.pft.lifeform != GRASS) {
			vegetation.nextobj();
			continue;
		}

		// Local:
		double cstore = 0.0;
		double cexcess = 0.0;
		double cstorage_res = 0.0;
		double Gpot = 0.0;

		double ceaten = 0.0;

		// Maximum C and N storage pool size
		indiv.cstorage_max = cstorage_capacity(indiv);
		indiv.nstorage_max = nstorage_capacity(indiv);

		// C storage not available for growth
		cstorage_res = cstorage_reserve(indiv);

		// Potential growth rate
		Gpot = max(0.0, max_mobil_cstorage * max(0.0,(indiv.cmass_storage - cstorage_res)) + indiv.dnpp); 

		// Determine growth demand based on growth stage
		indiv.growth_demand = growth_demand(indiv);

		// Set indiv.cgrow and indiv.cexcess to 0
		indiv.cgrow = 0.0;
		indiv.cexcess = 0.0;
		indiv.nmass_missing = 0;

		if (indiv.dnpp >= 0.0) {
			
			// Positive NPP

			// TODO: make clear this only works for daily time steps (NPP in gC m-2 d-1)	
			indiv.cgrow = indiv.growth_demand * Gpot;
			cstore = indiv.dnpp - indiv.cgrow;  // = Gpot - growth - (Gpot - NPP) 
												// cstore = re-allocation if negative

			// would it help to have an absolute growth demand instead of a relative growth_demand?
			
			// If C storage exceed its capacity take track of the excess C and reduce cstore
			if (indiv.cmass_storage + cstore > indiv.cstorage_max) {
				cexcess += indiv.cmass_storage + cstore - indiv.cstorage_max;
				indiv.cexcess += cexcess;
				cstore = indiv.cstorage_max - indiv.cmass_storage;
			}

		} else {

			// Negative NPP

			// Respiration initially from storage (labile) pool
			// Storage pool can meet demand. Still check if growth should occur.
			if (indiv.cmass_storage + indiv.dnpp >= 0.0) {
				indiv.cgrow = indiv.growth_demand * Gpot;
				cstore = indiv.dnpp - indiv.cgrow;
			}
			else {
				// If C balance still negative respiration from existing leaf and root
				// biomass and retranslocate corresponding N to keep leaf and root 
				// C:N ratio
				double cmass_deficit = -(indiv.cmass_storage + indiv.dnpp);
				cstore = -indiv.cmass_storage;
		
				// Enough C to compensate for negative C storage
				// Only allow 10% of biomass loss and call mortality if more than that.
				// 10% Consistent with senescence retranslocation but needs a reference...
				// TODO: replace 0.1 with constant.
				if ( 0.1 * (indiv.get_cmass_leaf() + indiv.cmass_root) > cmass_deficit && 
					 0.1 * (indiv.get_cmass_leaf() + indiv.cmass_root) > (indiv.cdeficit_5.sum() + cmass_deficit) ) {

					indiv.cdeficit_5.add(cmass_deficit); 

					// Reduce leaf to root ratio
					double froot_dec = 0.0;
					if (indiv.cmass_root > 0.0) {
						double current_ltor = indiv.cmass_root > 0 ? indiv.get_cmass_leaf() / indiv.cmass_root : 0.0;
						froot_dec = 1.0 / (current_ltor + 1.0);
					}

					// Reduce root biomass
					double ntoc = indiv.nmass_root / indiv.cmass_root;
					double cdec = cmass_deficit * froot_dec;  // take half from roots if ltor is 1
					double ndec = cdec * ntoc;
					indiv.cmass_root -= cdec;
					indiv.nmass_root -= ndec;
					cmass_deficit -= cdec;
					indiv.nmass_storage += ndec;
					ceaten += cdec;

					assert(positive(indiv.cmass_root));
					assert(positive(indiv.nmass_root));

					// Reduce leaf biomass (green pool only)
					if (indiv.cmass_leafs[green] > 0.0) {
						double ntoc = indiv.nmass_leafs[green] / indiv.cmass_leafs[green];
						double cdec = cmass_deficit;
						double ndec = cdec * ntoc;
						indiv.cmass_leafs[green] -= cdec;
						indiv.nmass_leafs[green] -= ndec;
						cmass_deficit -= cdec;
						indiv.nmass_storage += ndec;
						ceaten += cdec;
					}

					assert(positive(indiv.cmass_leafs[green]));
					assert(positive(indiv.nmass_leafs[green]));
					assert(equal(cmass_deficit, 0));

				}
				else {
					// Not enough C to compensate --> kill individual

					// The respiration and NPP fluxes reported in canexch are
					// now incorrect, so we need to fix them. Note that the RA
					// flux is not used for mass balance checks (it's only used
					// for outputs). The NPP flux is used in mass balance checks
					// however, so it's important to correct it here.
					indiv.report_flux(Fluxes::RA, -cmass_deficit);
					 // (npp has opposite sign to respiration)
					indiv.report_flux(Fluxes::NPP, cmass_deficit);

					kill_grass(indiv);
										
					//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
					// Don't think we need to call nextobj() here, as we've just
					// removed an element from the list; when this happens, the
					// "current" pointer internal to the list object points to
					// the element that was previously the "next" element.
					continue;
				}
			}
		}

		assert(positive(indiv.cmass_storage));
		assert(equal(indiv.dnpp, indiv.cgrow + cstore + cexcess - ceaten));

		// update cmass_storage
		indiv.cmass_storage += cstore;

		vegetation.nextobj();
	}

}


// Determine C allocation fractions to roots, leaves, and reproductive tissue in grasses.
// Note: 'leaf' includes shoot and leaves and 'root' includes coarse and fine roots.
void carbon_allocation_grass(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (pft.lifeform != GRASS) {
			continue;
		}

		// Leaf nitrogen status
		if (!negligible(indiv.get_cmass_leaf())) {
			if (!ifnlim || date.year < freenyears) {
				indiv.nscal_5.add(1.0);
			} else {
				indiv.nscal_5.add(min(1.0, indiv.pft.cton_leaf_avr / indiv.cton_leaf()));
			}
		}
		

		double calloc_repr = 0;
		double calloc_leaf = 0;
		double calloc_root = 0;
		if (!negligible(indiv.cgrow)) {
			// determine if some of the new carbon goes to reproduction
			if (indiv.growth_demand < 0.5 && indiv.lai > 2.0){  // TODO: define threshold
				calloc_repr = indiv.pft.reprfrac;    // TODO: move to constants
			} else {
				calloc_repr = 0.0;
			}

			// express ltor_max as maximum fractional allocation to leaves
			double calloc_leaf_max = 1.0 - 1.0/(1.0 + indiv.pft.ltor_max);

			// calculate fraction of C allocated to leaves following Johnson 2016 Eq. 3.80
			calloc_leaf = calloc_leaf_max * sqrt(indiv.nscal_5.mean() * ppft.wscal_5.mean());

			// allocation to leaves and roots
			calloc_leaf = calloc_leaf * (1.0 - calloc_repr);
			calloc_root = 1.0 - calloc_leaf - calloc_repr;
		}
		else {
			calloc_leaf = 0.0;
			calloc_root = 1.0;
			calloc_repr = 0.0;
		}

		assert(equal(1.0, calloc_leaf + calloc_root + calloc_repr));

		indiv.cgrow_leaf = calloc_leaf * indiv.cgrow;
		indiv.cgrow_root = calloc_root * indiv.cgrow;
		indiv.cgrow_repr = calloc_repr * indiv.cgrow;
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
	// Can this happen on the first day?
	if (nmass == 0 || cgrow == 0) {
		return 1;
	}

	// Get C:N ratio of existing biomass. If it's within the allowed bounds then
	// there's nothing more to do.
	double cton = cmass / nmass;

	if (cton > cton_max || cton < cton_min) {
		// Current C:N ratio is outside the allowed C:N ratio. We need to
		// up- or down-regulate the N demand in order to move towards the
		// nearest C:N ratio bounds.
		double target_cton = cton > cton_max ? cton_max : cton_min;

		// Get the N content which would put us at the target C:N ratio after
		// today's growth.
		double c_new = cmass + cgrow;
		double n_new = c_new / target_cton;

		// Get the delta between current and target N content.
		double ngrow = n_new - nmass;

		// Calculate the C:N ratio of new growth which will put us at the target
		// C:N ratio.
		cton = cgrow / ngrow;
	}
	assert(positive(cton));
	return cton;
}


// Calculate N demand based on growth and N status of the plant
	
void ndemand(Vegetation& vegetation) {
	const double frac_nchange = 0.01;   // TODO: better to express as C:N? Zaehle as reference?

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		// Local	
		double ngrow_demand_repr = 0;
		double ngrow_demand_leaf = 0;
		double ngrow_demand_root = 0;

		// double cton_root_min = 2.0 * indiv.pft.cton_root_avr - indiv.pft.cton_root_max;
		// cton_root = max(cton_root, cton_root_min);

		// N demand of new growth
		if (!negligible(indiv.cgrow)) {
			// Current C:N ratios
			double cton_leaf = constrain_cton(
				indiv.get_cmass_leaf(), indiv.get_nmass_leaf(),
				indiv.pft.cton_leaf_min, indiv.pft.cton_leaf_max,
				indiv.cgrow_leaf);
			double cton_root = constrain_cton(
				indiv.cmass_root, indiv.nmass_root,
				0, indiv.pft.cton_root_max,
				indiv.cgrow_root);
			ngrow_demand_leaf = indiv.cgrow_leaf / cton_leaf;
			ngrow_demand_root = indiv.cgrow_root / cton_root;
			ngrow_demand_repr = indiv.cgrow_repr / cton_leaf;  // C:N repr = C:N leaf for now.
		} else {
			ngrow_demand_leaf = 0.0;
			ngrow_demand_root = 0.0;
			ngrow_demand_repr = 0.0;
		}

		// N demand of existing biomass (will also be applied to new biomass)
		double nupdate_demand_leaf = max(0.0, indiv.get_cmass_leaf() / indiv.pft.cton_leaf_min - indiv.get_nmass_leaf()) * frac_nchange;
		double nupdate_demand_root = max(0.0, indiv.cmass_root / (indiv.pft.cton_leaf_min * indiv.pft.cton_root_avr /indiv.pft.cton_leaf_avr)
									    - indiv.nmass_root) * frac_nchange;

		// Total biomass N demand of individual
		indiv.leafndemand = ngrow_demand_leaf + nupdate_demand_leaf;
		indiv.rootndemand = ngrow_demand_root + nupdate_demand_root;
		indiv.reprndemand = ngrow_demand_repr;
		indiv.ndemand = indiv.leafndemand + indiv.rootndemand + indiv.reprndemand;

		assert(positive(indiv.ndemand));
	}
}



// Calculate N uptake which is limited by competition between individuals and
// availability of N in the soil.
void nuptake(Vegetation& vegetation) {
	Patch& patch = vegetation.patch;
	Soil& soil = patch.soil;

	// Create vector of CompetingIndividual objects, which will be passed into
	// ncompete(). These are populated from the Individuals in the patch.
	std::vector<NCompetingIndividual> competing_individuals;
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		// Calculate whole-plant N:C status 
		double ntoc = (indiv.get_nmass_leaf() + indiv.nmass_root) / (indiv.get_cmass_leaf() + indiv.cmass_root);
		double ntoc_min = 1.0 / indiv.pft.cton_leaf_min;
		double ntoc_avr = 1.0 / indiv.pft.cton_leaf_avr;
		indiv.cton_status = max(0.0, (ntoc - ntoc_min) / (ntoc_avr - ntoc_min));

		NCompetingIndividual competing_individual;
		competing_individual.ndemand = indiv.ndemand;
		competing_individual.fnuptake = 1;
		competing_individual.strength = nitrogen_uptake_strength(indiv);
		competing_individuals.push_back(competing_individual);
	}

	// Calculate N competition.
	double nmass_avail = soil.nmass_avail() * min(1.0, patch.fpc_total);
	ncompete(competing_individuals, nmass_avail);  // TODO: check if OK/efficient to call ncompete even if ndemand < nsupply. 
												   // --> if (nmass_avail < patch_total_ndemand) {

	// Copy results back into Individual objects in the patch, and calculate
	// total patch N demand.
	patch.ndemand = 0;
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		indiv.fnuptake = competing_individuals[i].fnuptake;
		indiv.ndemand = competing_individuals[i].ndemand;

		patch.ndemand += indiv.ndemand;

		// Taken out of soil in somdynam::vegetation_n_uptake():
		// indiv.nmass_storage += indiv.ndemand * indiv.fnuptake;
		// This happens right after growth in simulate_day.
	}
}



// Update C and N pools of vegetation  
void update_pools(Vegetation& vegetation) { 
	
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

		if (pft.lifeform != GRASS) {
			continue;
		}

		double creturn = 0.0; // kgC m-2 d-1
		double cexcess = 0.0, nexcess = 0.0;
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
			nexcess = indiv.nmass_storage - indiv.nstorage_max;
			indiv.nmass_litter_root += nexcess;
			ppft.nmass_litter_root += nexcess;
			indiv.nmass_storage = indiv.nstorage_max;
		}

		// TODO: Check if there is enough N in storage. At the moment Nstorage needs to be > 0.
		// Alternatively, one could implement a N:C storage > N:C leaf min constraint.
		
		double fnavail = 0.0;
		if (!negligible(indiv.ndemand)) {
			// fnavail = max(0.0,min(1.0,(indiv.nmass_storage - indiv.cmass_storage / indiv.pft.cton_leaf_max) / indiv.ndemand));
			fnavail = max(0.0,min(1.0,(indiv.nmass_storage / indiv.ndemand)));
		}

		double navail_leaf = indiv.leafndemand * fnavail;
		double navail_root = indiv.rootndemand * fnavail;	
		double navail_repr = indiv.reprndemand * fnavail; 

		// Check with fnavail, not fnuptake. There are cases where fnuptake can
		// be 1 but fnavail < 1. This can happen when n storage is low, because
		// the uptake doesn't occur until after growth. If we were to check
		// fnuptake in this case, we would incorrectly assume that we are able
		// to meet the N demands today, and potentially allocate too much C,
		// causing the C:N ratio to exceed its max.
		if (negligible(fnavail - 1.0)) {
			
			// N supply >= N demand

			indiv.cmass_leafs[green] += indiv.cgrow_leaf;
			indiv.cmass_root += indiv.cgrow_root;
			indiv.cmass_repr += indiv.cgrow_repr;
			
		}
		else {

			// N supply < N demand
			double nreturn = 0;
			if ( navail_leaf < (indiv.cgrow_leaf / indiv.pft.cton_leaf_max)) {
				
				// Maximum possible leaf growth given N constraints. Put rest of C back into storage (Creturn)
				double growth_max_leaf = navail_leaf * indiv.pft.cton_leaf_max;
				creturn += indiv.cgrow_leaf - growth_max_leaf;
				indiv.cgrow_leaf = growth_max_leaf;
			} else if ( navail_leaf > (indiv.cgrow_leaf / indiv.pft.cton_leaf_min) ) {
				// N:C ratio in new leaf exceeds max N:C ratio. Limit N in new
				// leaf to the max N:C ratio and put the rest into nreturn.
				double max_ngrow_leaf = indiv.cgrow_leaf / indiv.pft.cton_leaf_min;
				nreturn += navail_leaf - max_ngrow_leaf;
				navail_leaf = max_ngrow_leaf;
			}

			double cton_root_min = indiv.pft.cton_root_avr - (indiv.pft.cton_root_max - indiv.pft.cton_root_avr);
			if ( navail_root < (indiv.cgrow_root / indiv.pft.cton_root_max)) {
				
				// Maximum possible root growth given N constraints. Put rest of C back into storage (Creturn)
				double growth_max_root = navail_root * indiv.pft.cton_root_max;
				creturn += indiv.cgrow_root - growth_max_root;
				indiv.cgrow_root = growth_max_root;
			} else if (navail_root > (indiv.cgrow_root / cton_root_min)) {
				// N:C ratio in new root exceeds max N:C ratio. Limit N in new
				// root to the max N:C ratio and put the rest into nreturn.
				double max_ngrow_root = indiv.cgrow_root / cton_root_min;
				nreturn += navail_root - max_ngrow_root;
				navail_root = max_ngrow_root;
			}

			if (navail_repr < (indiv.cgrow_repr / indiv.pft.cton_leaf_max)) {

				// Maximum possible reproductive growth given N constraints. Put rest of C back into storage (Creturn)
				double growth_max_repr = navail_repr * indiv.pft.cton_leaf_max;
				creturn += indiv.cgrow_repr - growth_max_repr;
				indiv.cgrow_repr = growth_max_repr;
			} else if (navail_repr > (indiv.cgrow_repr / indiv.pft.cton_leaf_min)) {
				// N:C ratio in new reproductive biomass exceeds max N:C ratio.
				// Limit N in new reproductive biomass to the max N:C ratio and
				// put the rest into nreturn.
				double max_ngrow_repr = indiv.cgrow_repr / indiv.pft.cton_leaf_min;
				nreturn += navail_repr - max_ngrow_repr;
				navail_repr = max_ngrow_repr;
			}

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

			indiv.cmass_leafs[green] += indiv.cgrow_leaf;
			indiv.cmass_root += indiv.cgrow_root;
			indiv.cmass_repr += indiv.cgrow_repr;

			assert(equal(indiv.ndemand * fnavail, navail_leaf + navail_root + navail_repr + nreturn));
		}

		// Always take up available N
		indiv.nmass_leafs[green] += navail_leaf;
		indiv.nmass_root += navail_root;
		indiv.nmass_repr += navail_repr;
		indiv.nmass_storage = max(0.0, indiv.nmass_storage - (navail_leaf + navail_root + navail_repr));

		// Save yearly leaf and root growth
		indiv.cmass_leaf_ygrowth += indiv.cgrow_leaf;
		indiv.cmass_root_ygrowth += indiv.cgrow_root;

		// update indiv.cgrow (actual growth)
		indiv.cgrow = indiv.cgrow_leaf + indiv.cgrow_root + indiv.cgrow_repr;

        assert(positive(indiv.cgrow_repr));
		assert(positive(indiv.cgrow_leaf));
		assert(positive(indiv.cgrow_root));
		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));
		assert(less_than_or_equal(indiv.nmass_storage, indiv.nstorage_max));
		assert(equal(cgrow_unlim, indiv.cgrow_repr + indiv.cgrow_leaf + indiv.cgrow_root + creturn + cexcess));
		assert(less_than_or_equal(indiv.cton_leaf(), indiv.pft.cton_leaf_max));
		assert(greater_than_or_equal(indiv.cton_leaf(), indiv.pft.cton_leaf_min));
		assert(less_than_or_equal(indiv.cton_root(), indiv.pft.cton_root_max));
	}
}

// Transfer of C and N from green (live) leaf pool to brown (dead) leaf pool
// Accounts for retranslocation of C and N to storage prior to senescence. 
void senescence(Vegetation& vegetation) {
	
	// Parameters
	const double crelocfrac_c3 = 0.1;  // TODO: move to constants?
	const double crelocfrac_c4 = 0.05; 
	const double nrelocfrac = 0.5;
	//const double k = 10.0;            // Shape parameter in senescence function (steepness)
	const double qwater_senesc = 0.5;

	// Air temperature or bottom snow layer if snow present
	double grass_temp = vegetation.patch.soil.T_soil[6];

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (pft.lifeform != GRASS) {
			continue;
		}

		// Senescence

		// Baseline (minimum) senescence rate from leaf longevity under unstressed conditions.
		double senesc_base = 1.0 / (indiv.pft.leaflong * date.year_length());
		
		// Stress factors increase senescence 
		// Soil moisture
		double senesc_fwater = senesc_base * pow((1 - ppft.wscal),qwater_senesc) * 3.0;
		
		// High and low temperatures
		double senesc_ftemp  = 0.0;  // TODO: replace Temp response
		if (grass_temp < 10.0) {
			senesc_ftemp = senesc_base * (10.0 - grass_temp) / 10.0;
		}
		else if (grass_temp > 35.0) {
			senesc_ftemp = senesc_base * (35.0 + grass_temp) / 10.0;
		}

		// Dormancy-induced (drought dormancy only at the moment)
		double senesc_fdorm = 0;
		if (indiv.fphen_activity < 1.0 && (ppft.wscal_10.size() > 0 && ppft.wscal_10.mean() < pft.wscal_min)) {

			senesc_fdorm = senesc_base * (1.0 - indiv.fphen_activity) * 5.0;

		}

		// Total senescence rate
		double senesc = senesc_base + senesc_fwater + senesc_ftemp + senesc_fdorm;

		// Total senescence rate
		//senesc = senesc_max / (1 + exp( -k * (min(1.0,indiv.growth_stage - growth_stage0)) ));
		
		// Retranslocation and pool updates
		// TODO: maybe call cstorage_capacity once at the beginning to be consistent across routines!
		// TODO: calculate retranslocation coefficients based on current C:N ratios?


		// TODO: move fractions to parameters??
		double crelocfrac = indiv.pft.pathway == C3 ? crelocfrac_c3 : crelocfrac_c4;

		// Carbon
		double csenesc  = indiv.cmass_leafs[green] * senesc;
		double cstore   = csenesc * crelocfrac;
		double ctobrown = csenesc * (1.0 - crelocfrac);

		// If C storage pool overflows, re-adjust cstore and ctobrown
		if (indiv.cmass_storage + cstore > indiv.cstorage_max) {
			cstore = indiv.cstorage_max - indiv.cmass_storage;
			ctobrown = csenesc - cstore;
		}

		indiv.cmass_leafs[green] -= csenesc;   // green leaf pool
		indiv.cmass_leafs[brown] += ctobrown;  // brown leaf pool
		indiv.cmass_storage  += cstore;    // C storage

		// Nitrogen
		double nsenesc  = indiv.nmass_leafs[green] * senesc;
		double nstore   = nsenesc * nrelocfrac;
		double ntobrown = nsenesc * (1.0 - nrelocfrac);

		// If N storage pool overflows, re-adjust nstore and ntobrown
		if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
			nstore = indiv.nstorage_max - indiv.nmass_storage;
			ntobrown = nsenesc - nstore;
		}

		// TODO: make sure N:C ratio in storage does not fall outside its ranges
		//       understand what the N:C in storage looks like (values). How flexible is N:C in storage?
		//		 sufficient to only account for Nstorage_max as done here?

		indiv.nmass_leafs[green] -= nsenesc;   // green leaf pool
		indiv.nmass_leafs[brown] += ntobrown;  // brown leaf pool
		indiv.nmass_storage  += nstore;        // N storage

		assert(positive(indiv.cmass_leafs[green]));
		assert(positive(indiv.nmass_leafs[green]));
		assert(positive(indiv.cmass_leafs[brown]));
		assert(positive(indiv.nmass_leafs[brown]));

		assert(positive(indiv.cmass_storage));
		assert(positive(indiv.nmass_storage));

		assert(equal(csenesc, cstore + ctobrown));
		assert(equal(nsenesc, ntobrown + nstore));
	}
}


// Transfer of brown leaf pool and fine root pool to leaf and root litter.
// Transfer reproductive tissue to litter if end of growth cycle is reached.
void turnover_grass(Vegetation& vegetation) {

    // TODO: move to constants
	const double nrelocfrac = 0.5;

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (pft.lifeform != GRASS) {
			continue;
		}

		// Determine turnover rate (same for all individuals in a PFT).
		// Assumed daily rates here (fraction d-1)!
		// TODO: turnover rates should be moisture and temperature dependent
		double frac_turnover_leaf = indiv.pft.turnover_leaf / date.year_length();
		//double frac_turnover_root = (pow((1 - ppft.wscal),0.4) * 4.0) * indiv.pft.turnover_root / date.year_length();
		double frac_turnover_root = 3.0 * indiv.pft.turnover_root / date.year_length();

		// Dead leaf turnover (brown leaves --> leaf litter)
		if (indiv.cmass_leafs[brown] > 0.0) {

			// Carbon
			double cturnover = indiv.cmass_leafs[brown] * frac_turnover_leaf;

			indiv.cmass_leafs[brown]   -= cturnover;
			indiv.cmass_litter_leaf    += cturnover;
			ppft.cmass_litter_leaf     += cturnover;


			// Nitrogen
			double nturnover = indiv.nmass_leafs[brown] * frac_turnover_leaf;
			
			indiv.nmass_leafs[brown]   -= nturnover;
			indiv.nmass_litter_leaf    += nturnover;
            ppft.nmass_litter_leaf     += nturnover;

		}


		// Root turnover (fine roots --> root litter)
		// TODO: should this be in senescence as well? Dead root pool?
		if (indiv.cmass_root > 0.0) { 

			// TODO: include C resorption from roots?
			// Carbon
			double cturnover = indiv.cmass_root * frac_turnover_root;
			
			indiv.cmass_root        -= cturnover;
			indiv.cmass_litter_root += cturnover;
			ppft.cmass_litter_root  += cturnover;


			// Nitrogen (including retranslocation)
			double nturnover = indiv.nmass_root * frac_turnover_root;
			double nstore    = nturnover * nrelocfrac;
			double ntolitter = nturnover * (1.0 - nrelocfrac);

		    if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
				nstore    = indiv.nstorage_max - indiv.nmass_storage;	
				ntolitter = nturnover - nstore;
			}

			indiv.nmass_root        -= nturnover;
			indiv.nmass_litter_root += ntolitter;
            ppft.nmass_litter_root  += ntolitter;

			indiv.nmass_storage    += nstore;

			assert(equal(nturnover, nstore + ntolitter));
		}

		// Move reproductive C and N to litter at the beginning of dormancy
		if (indiv.fphen_activity >= 1.0 && indiv.cmass_repr > 0.0) {
			ppft.cmass_litter_repr += indiv.cmass_repr;
			indiv.cmass_litter_repr += indiv.cmass_repr;
			indiv.nmass_litter_leaf += indiv.nmass_repr;
			ppft.nmass_litter_leaf += indiv.nmass_repr; // TODO: decide where to put N from repr. 
														//       Added to leaf litter for now.
														//       Call report flux function??
			// TODO: for annuals, flush everything to litter at the end of growing season.
			//		 Simulate seeds as cmass_storage and nmass_storage.
			// TODO: call establishment function every day for grasses?
			indiv.cmass_repr = 0;
			indiv.nmass_repr = 0;
		}
	}
}


// Determine what to do with excess C when storage exceeds its capacity.
void cstorage_overflow(Individual& indiv) {

	// TODO: decide what to do with excess C.
	// Option 1: respire
	indiv.report_flux(Fluxes::RA, indiv.cexcess);

	// Option 2; root exudates
	// TODO: yet to be implemented

	// Option 3?

	// Reset cexcess to 0.0
	// Now done at the beginning of the day in storage_growth_partitioning() so that mass balance checks work.
	// indiv.cexcess = 0.0;

}


// Photosynthetic sink limitation.
// Calculates downregulation factors applied to photosynthetic parameters.
void sink_source_feedback(Individual& indiv) {
	
	// Do not let ps_sink_lim vary by more than 10% d-1
	// TODO: revisit cexcess / dgpp
	if (indiv.cexcess > 0.0 && indiv.dgpp > 0.0) {
		indiv.ps_sink_lim = max(0.9 * indiv.ps_sink_lim, (1.0 - min(1.0, indiv.cexcess / indiv.dgpp)));
	} 
	else {
		indiv.ps_sink_lim = min(1.0, indiv.ps_sink_lim * 1.1); 
	}

} 


// Housekeeping of C storage. Handles overflow and feedback with source dynamics.
void storage_dynamics(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (indiv.pft.lifeform != GRASS) {
			continue;
		}
	
		// Handle excessive C 
		if (indiv.cexcess > 0.0){
			cstorage_overflow(indiv);
		}

		// Determine in how far C storage induces a negative feedback on photosynthetic capacity.
		sink_source_feedback(indiv);
	}
}


// Updates LAI and FPC for grasses.
// Equivalent to default allometry() function at the moment,
// but might develop. (e.g. more than 1 indiv for grasses, separate green and brown lai, 
// specific extinction coefficient etc.)
void grass_allometry(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		if (indiv.pft.lifeform != GRASS) {
			continue;
		}
		allometry(indiv);
	}
}


// Aggregate fluxes for annual output.
void annual_accounting(Vegetation& vegetation) {

	// Loop through individuals
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];

		if (indiv.pft.lifeform != GRASS) {
			continue;
		}

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

void growth_daily_grass_dave(Patch& patch) {

	// Obtain reference to Vegetation and Stand object for this patch
	Vegetation& vegetation = patch.vegetation;
	Stand& stand = patch.stand;

	// Record mass balance at start of day.
	DaveMassBalance mass_balance;
	mass_balance.init(vegetation);

	if (date.day == 0) {
		flush_cmass_repr(patch);
	}

	// Calculate total daily growth and update C storage 
	storage_growth_partitioning(vegetation);

	// Allocate available C for growth to different plant pools
	carbon_allocation_grass(vegetation);
	
	// Calculate N demand based on growth and stoichiometry.
	ndemand(vegetation);

	// Calculate actual N uptake after taking competition among individuals into account
	nuptake(vegetation);

	// update pools (C + N) to account for growth
	update_pools(vegetation);
	mass_balance.check(vegetation);

	// Leaf senescence and retranslocation.
	senescence(vegetation);
	mass_balance.check(vegetation);

	// Dead leaf and root to litter
	turnover_grass(vegetation);	
	mass_balance.check(vegetation);

	// C storage housekeeping and feedback to photosynthesis
	storage_dynamics(vegetation);
	mass_balance.check(vegetation);

	// Management
	// grazing()
	//grazing_daily(indiv);

	// Allometry (Updates LAI and FPC)
	grass_allometry(vegetation);
	mass_balance.check(vegetation);

	// C balance checks

	// yearly reporting
	annual_accounting(vegetation);

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

	if (ifdailygrass) {
		
		// daily C allocation (grass only)
		if (ifdavegrowth) {
			growth_daily_grass_dave(patch);
		} else {
			growth_daily_grass(patch);
		}
		// necessary updates after daily change to fpc in daily grass
		update_patch_fpc(patch);
	}
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
