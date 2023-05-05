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

// Maximum size of C storage pool as a fraction of total biomass (leaf + root). Grasses only.
const double frac_cstorage_max = 0.2;

// maximum mobilisation rate of C in storage
const double max_mobil_cstorage = 0.2;

class DaveMassBalance {
public:
	~DaveMassBalance();
	void init(Vegetation &vegetation);
	void check(Vegetation &vegetation);
private:
	class IndivState {
public:
		IndivState(const Individual &indiv);
		void check(const Individual &indiv) const;
private:
		double cmass_storage;
		double cmass_leaf_green;
		double cmass_leaf_brown;
		double cmass_root;
		double cmass_repr;
		double cmass_litter_leaf;
		double cmass_litter_root;
		double cmass_litter_repr;
		double nmass_storage;
		double nmass_leaf_green;
		double nmass_leaf_brown;
		double nmass_root;
		double nmass_repr;
		double nmass_litter_leaf;
		double nmass_litter_root;
	};
	std::vector<const IndivState*> indivs;
	void clear_indivs();
};

DaveMassBalance::~DaveMassBalance() {
	clear_indivs();
}

void DaveMassBalance::clear_indivs() {
	for (const IndivState *indiv : indivs) {
		delete indiv;
	}
	indivs.clear();
}

void DaveMassBalance::init(Vegetation &vegetation) {
	clear_indivs();
	for (vegetation.firstobj(); vegetation.isobj; vegetation.nextobj()) {
		IndivState *indiv = new IndivState(vegetation.getobj());
		indivs.push_back(indiv);
	}
}

void DaveMassBalance::check(Vegetation &vegetation) {
	vegetation.firstobj();
	int i = 0;
	while (vegetation.isobj) {
		indivs[i]->check(vegetation.getobj());
		vegetation.nextobj();
		i++;
	}
}

DaveMassBalance::IndivState::IndivState(const Individual &indiv) {
	cmass_storage = indiv.cmass_storage;
	cmass_leaf_green = indiv.cmass_leafs[green];
	cmass_leaf_brown = indiv.cmass_leafs[brown];
	cmass_root = indiv.cmass_root;
	cmass_repr = indiv.cmass_repr;
	nmass_storage = indiv.nmass_storage;
	nmass_leaf_green = indiv.nmass_leafs[green];
	nmass_leaf_brown = indiv.nmass_leafs[brown];
	nmass_root = indiv.nmass_root;
	nmass_repr = indiv.nmass_repr;
	cmass_litter_leaf = indiv.cmass_litter_leaf;
	cmass_litter_root = indiv.cmass_litter_root;
	cmass_litter_repr = indiv.cmass_litter_repr;
	nmass_litter_leaf = indiv.nmass_litter_leaf;
	nmass_litter_root = indiv.nmass_litter_root;
}

void DaveMassBalance::IndivState::check(const Individual &indiv) const {
	double c_start = cmass_storage + cmass_leaf_green + cmass_leaf_brown + cmass_root + cmass_repr + cmass_litter_leaf + cmass_litter_repr + cmass_litter_root;
	double c_flux = indiv.dnpp;
	double c_end = indiv.cmass_storage + indiv.cmass_leafs[green] + indiv.cmass_leafs[brown] + indiv.cmass_root + indiv.cmass_repr + indiv.cmass_litter_leaf + indiv.cmass_litter_repr + indiv.cmass_litter_root + indiv.cexcess;
	assert(negligible(c_end - (c_start + c_flux), -14));

	double n_start = nmass_storage + nmass_leaf_green + nmass_leaf_brown + nmass_root + nmass_repr + nmass_litter_leaf + nmass_litter_root;
	double n_end = indiv.nmass_storage + indiv.nmass_leafs[green] + indiv.nmass_leafs[brown] + indiv.nmass_root + indiv.nmass_repr + indiv.nmass_litter_leaf + indiv.nmass_litter_root;
	assert(negligible(n_end - n_start, -14));
}

///////////////////////////////////////////////////////////////////////////////////////
// LEAF PHENOLOGY
// Call function leaf_phenology each simulation day prior to calculation of FPAR, to
// calculate fractional leaf-out for each PFT and individual.
// Function leaf_phenology_pft is not intended to be called directly by the framework,


void phenology_grass_dave(Vegetation& vegetation, Climate& climate) {

	// Calculates phenological cycle of grasses. 
	// Returns phenological state of grasses as factor fphen [0-1].
	// Implemented for grasses only.

	// Loop through cohorts (=individuals)
	vegetation.firstobj();
	while (vegetation.isobj) {

		// Retrieve information about cohort
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft &ppft = vegetation.patch.pft[indiv.pft.id];

		if (pft.lifeform == GRASS) { 

			// Sum up growing degree days (5degC base temperature)
			double gdd = climate.temp - 5.0;

			// Reset growth stage if maturity has been reached the day before.
			if (indiv.growth_stage >= 1.0) {
				indiv.growth_stage = 0.0;
			} 
			// Determine if grass is dormant based on climate conditions of the last 10 days.
			// If not, determine phenological state (0=germination/start growing season to 1 = full senescence].
			else if ((ppft.wscal_5.size() > 0 && ppft.wscal_5.mean() < pft.wscal_min) || climate.temp < 0.0) { 
				indiv.growth_stage = 0.0;
			} else {
				indiv.growth_stage = min(1.0, indiv.growth_stage + gdd / pft.phengdd5mat);
			}
		} else {
			indiv.growth_stage = 1.0;
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

	if (negligible(cmass_storage, -12)) {
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
			}
			indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_avr;
			indiv.nmass_storage = indiv.cmass_storage / indiv.pft.cton_leaf_avr;
		}

		indiv.alive = true;
		indiv.report_flux(Fluxes::ESTC, -(indiv.cmass_storage + indiv.get_cmass_leaf() + indiv.cmass_root) );
		indiv.report_flux(Fluxes::ESTN, -(indiv.nmass_storage + indiv.get_nmass_leaf() + indiv.nmass_root) );
	}
	else {

		allocation(bminit, indiv, ltor, cmass_leaf_ind, cmass_root_ind, cmass_sap_ind, dval, dval, dval, dval, dval);

		indiv.cmass_leaf = cmass_leaf_ind * indiv.densindiv;
		indiv.cmass_root = cmass_root_ind * indiv.densindiv;

		if (ifdavephotosyn){
			// we need initial values for nmass_leaf when calling ifdavephotosyn
			indiv.nmass_leaf = indiv.cmass_leaf / indiv.pft.cton_leaf_avr;
			indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_avr;
		}
		else {
			indiv.nmass_leaf = 0.0;
			indiv.nmass_root = 0.0;
		}
	}

	if (indiv.pft.lifeform == TREE) {
		indiv.cmass_sap = cmass_sap_ind * indiv.densindiv;
		indiv.nmass_sap = 0.0;
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

				// Grass "individual" LAI (Eqn 11)
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
	if (!patch.id && !ifdailygrass){ 
	
		for (p=0; p<npft; p++){
			stand.pft[p].cmass_repr = 0.0;
		}
	}

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
							indiv.andemand_leaf = cmass_leaf_inc * indiv.cton_leaf();
							indiv.andemand_root = cmass_root_inc * indiv.cton_root();
							indiv.andemand_sap = cmass_sap_inc * indiv.cton_sap();
							indiv.andemand = indiv.andemand_leaf + indiv.andemand_root;
							// Note: sap_inc can be negative.
							if (indiv.andemand_sap > 0) {
								indiv.andemand += indiv.andemand_sap;
							}
							indiv.ndemand = indiv.andemand / 365; // fixme!! (Assaf)
							total_ndemand += indiv.ndemand;

							indiv.leaffndemand = indiv.andemand_leaf / indiv.andemand;
							indiv.rootfndemand = indiv.andemand_root / indiv.andemand;
							indiv.sapfndemand = indiv.andemand_sap > 0 ? indiv.andemand_sap / indiv.andemand : 0;

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
		if (!patch.id)
			for (int p=0; p < npft; p++)
				stand.pft[p].cmass_repr = 0.0;
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
			double ltor_status = (!negligible(indiv.get_cmass_leaf(), -6) && !negligible(indiv.cmass_root, -6)) ? max(0.5, min(2.0, indiv.pft.ltor_max / (indiv.get_cmass_leaf() / indiv.cmass_root))) : 1.0;
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
				else if (negligible(cgrow, -10) && nup_leaf_demand + nup_root_demand > 1.0e-12) {
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

// Determine growth demand (aka 'sink strength') 
// Returns fraction of available C used for growth vs. storage.
// Put fraction into growth and 1-fraction into storage.
// --> 1.0: put all available C of the day into growth.
// --> 0.0: put all available C of the day into storage.
double growth_demand(Individual& indiv) {

    // first approximation: growth demand inversely related to growth stage
	if (negligible(indiv.growth_stage)) { // dormant
		return 0.0; 
	} 
	else { 
		return 1.0 - indiv.growth_stage;
	}
	// TODO: do NOT call growth functions if grasses are dormant.
}

// Capacity of C storage (kg C m-2)
double cstorage_capacity(Individual& indiv) {

	return frac_cstorage_max * (indiv.get_cmass_leaf() + indiv.cmass_root);

}

// Capacity of N storage (kg N m-2)
double nstorage_capacity(Individual& indiv) {

	return indiv.cmass_storage / indiv.pft.cton_leaf_min;  // TODO: OK to use leaf C:N here? Different limit for storage?

}

// Determine how much of the available C (NPP + C fraction from storage) goes to growth vs. storage.
// Function calculates cgrow (daily growth), cstore (change in cmass_storage) 
// and updates indiv.cmass_storage. 
// Note: cgrow and cstore are fluxes (kgC/m2/d), cmass_storage is a pool (kg/m2).
void storage_growth_partitioning(Vegetation& vegetation) {
		
	// Local:
	double cstore = 0.0;
	double cexcess = 0.0;
	double Gpot = 0.0;

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		double ceaten = 0.0;

		// Determine growth demand based on growth stage
		indiv.growth_demand = growth_demand(indiv);

		// Potential growth rate
		Gpot = max(0.0, max_mobil_cstorage * indiv.cmass_storage + indiv.dnpp); 

		// Maximum C and N storage pool size
		indiv.cstorage_max = cstorage_capacity(indiv);
		indiv.nstorage_max = nstorage_capacity(indiv);

		// Set indiv.cgrow and indiv.cexcess to 0
		indiv.cgrow = 0.0;
		indiv.cexcess = 0.0;

		if (date.year == 1 && date.day == 182){
			printf("");
		}

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
				// Only allow 50% of biomass loss and call mortality if more than that
				if ( 0.5 * (indiv.get_cmass_leaf() + indiv.cmass_root) > cmass_deficit) {

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

					assert(indiv.cmass_root > 0);
					assert(indiv.nmass_root > 0);

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

					assert(indiv.cmass_leafs[green] > 0.0);
					assert(indiv.nmass_leafs[green] > 0.0);
					assert(negligible(cmass_deficit, -16));

				}
				else {
					// Not enough C to compensate --> kill individual

					indiv.patchpft().cmass_litter_repr += indiv.cmass_repr;
					indiv.patchpft().nmass_litter_leaf += indiv.nmass_repr;
					indiv.cmass_repr = 0.0;
					indiv.nmass_repr = 0.0;

					indiv.cmass_storage = 0.0;

					indiv.report_flux(Fluxes::RA, -cmass_deficit);
					indiv.kill();  // flushes all pools to litter (except for reproductive C, which is handled above)
					vegetation.killobj();  // removes individual object
					// TODO: think about whether we want to kill grasses or just keep them dormant.
					//       annuals might be a different case than perennials.
					//fail("Year %d Day %d Daily Grass Negative %g !!!!!!!!!\n", date.year, date.day, diff);
					vegetation.nextobj();
					continue;
				}
			}
		}

		assert(indiv.cmass_storage >= 0.0);
		assert(negligible(indiv.dnpp - (indiv.cgrow + cstore + cexcess - ceaten), -16));

		// update cmass_storage
		indiv.cmass_storage += cstore;

		vegetation.nextobj();
	}

}


// Determine C allocation fractions to roots, leaves, and reproductive tissue in grasses.
// Note: 'leaf' includes shoot and leaves and 'root' includes coarse and fine roots.
void carbon_allocation_grass(Vegetation& vegetation) {

	double calloc_repr, calloc_leaf, calloc_root;
	double calloc_leaf_max;
	
	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		// Leaf nitrogen status
		indiv.nscal_5.add(min(1.0, indiv.pft.cton_leaf_avr / indiv.cton_leaf()));

		if (!negligible(indiv.cgrow)) {
			// determine if some of the new carbon goes to reproduction
			if (indiv.growth_stage > 0.5 && indiv.growth_stage < 0.9){
				calloc_repr = indiv.pft.reprfrac;    // TODO: move to constants
			} else {
				calloc_repr = 0.0;
			}

			// express ltor_max as maximum fractional allocation to leaves
			calloc_leaf_max = 1.0 - 1.0/(1.0 + indiv.pft.ltor_max);

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

		assert(negligible(1.0 - (calloc_leaf + calloc_root + calloc_repr), -15));

		indiv.cgrow_leaf = calloc_leaf * indiv.cgrow;
		indiv.cgrow_root = calloc_root * indiv.cgrow;
		indiv.cgrow_repr = calloc_repr * indiv.cgrow;

		vegetation.nextobj();
	}
}



// Calculate N demand based on growth and N status of the plant
void ndemand(Vegetation& vegetation) {
	
	const double frac_nchange = 0.01;   // TODO: better to express as C:N? Zaehle as reference?

	// Local	
	double ngrow_demand_repr, ngrow_demand_leaf, ngrow_demand_root;
	double nupdate_demand_leaf, nupdate_demand_root;
	double cton_leaf, cton_root;

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// Current C:N ratios
		cton_leaf = !negligible(indiv.get_nmass_leaf()) ? indiv.get_cmass_leaf() / indiv.get_nmass_leaf() : indiv.pft.cton_leaf_avr;
		cton_root = !negligible(indiv.nmass_root) ? indiv.cmass_root / indiv.nmass_root : indiv.pft.cton_root_avr;

		// N demand of new growth
		if (!negligible(indiv.cgrow)) {
			ngrow_demand_leaf = indiv.cgrow_leaf / cton_leaf;
			ngrow_demand_root = indiv.cgrow_root / cton_root;
			ngrow_demand_repr = indiv.cgrow_repr / cton_leaf;  // C:N repr = C:N leaf for now.
		} else {
			ngrow_demand_leaf = 0.0;
			ngrow_demand_root = 0.0;
			ngrow_demand_repr = 0.0;
		}

		// N demand of existing biomass (will also be applied to new biomass)
		nupdate_demand_leaf = max(0.0, indiv.get_cmass_leaf() / indiv.pft.cton_leaf_min - indiv.get_nmass_leaf()) * frac_nchange;
		nupdate_demand_root = max(0.0, indiv.cmass_root / (indiv.pft.cton_leaf_min * indiv.pft.cton_root_avr /indiv.pft.cton_leaf_avr) -
									    - indiv.nmass_root) * frac_nchange;

		// Total biomass N demand of individual
		indiv.leafndemand = ngrow_demand_leaf + nupdate_demand_leaf;
		indiv.rootndemand = ngrow_demand_root + nupdate_demand_root;
		indiv.reprndemand = ngrow_demand_repr;
		indiv.ndemand = indiv.leafndemand + indiv.rootndemand + indiv.reprndemand;

		assert(indiv.ndemand >= 0);

		vegetation.nextobj();
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
	for (vegetation.firstobj(); vegetation.isobj; vegetation.nextobj()) {
		Individual& indiv = vegetation.getobj();

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
	int i = 0;
	for (vegetation.firstobj(); vegetation.isobj; vegetation.nextobj()) {
		Individual& indiv = vegetation.getobj();

		indiv.fnuptake = competing_individuals[i].fnuptake;
		indiv.ndemand = competing_individuals[i].ndemand;

		patch.ndemand += indiv.ndemand;

		// Taken out of soil in somdynam::vegetation_n_uptake():
		// indiv.nmass_storage += indiv.ndemand * indiv.fnuptake;
		// This happens right after growth in simulate_day.
		
		i++;
	}
}



// Update C and N pools of vegetation  
void update_pools(Vegetation& vegetation) { 

	// new parameters. TODO: move to pft parameters(?)
	// maximum mobilisation rate of N in storage
	//const double max_mobil_nstorage = 0.9;

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

        // ensure Nstorage does not fall below its minimum N:C (-> fnavail)
		double fnavail = 0.0;
		if (!negligible(indiv.ndemand)) {
			fnavail = max(0.0,min(1.0,(indiv.nmass_storage - indiv.cmass_storage / indiv.pft.cton_leaf_max) / indiv.ndemand));
		}
		double navail_leaf = indiv.leafndemand * fnavail * indiv.fnuptake;
		double navail_root = indiv.rootndemand * fnavail * indiv.fnuptake;	
		double navail_repr = indiv.reprndemand * fnavail * indiv.fnuptake; 
		double creturn = 0.0; // kgC m-2 d-1
		double growth_max_leaf = 0.0, growth_max_root = 0.0, growth_max_repr = 0.0;
		double cexcess = 0.0;


		if (negligible(indiv.fnuptake - 1.0)) {
			
			// N supply >= N demand

			indiv.cmass_leafs[green] += indiv.cgrow_leaf;
			indiv.cmass_root += indiv.cgrow_root;
			indiv.cmass_repr += indiv.cgrow_repr;
			
		}
		else {

			// N supply < N demand

			if ( navail_leaf < (indiv.cgrow_leaf / indiv.pft.cton_leaf_max)) {
				// TODO: add this later (no check at beginning of spinup) || (!ifnlim || date.year < freenyears)
				// Enough N for new growth
				
				// Maximum possible leaf growth given N constraints. Put rest of C back into storage (Creturn)
				growth_max_leaf = navail_leaf * indiv.pft.cton_leaf_max;
				creturn += indiv.cgrow_leaf - growth_max_leaf;
				indiv.cgrow_leaf = growth_max_leaf;
			}

			if ( navail_root < (indiv.cgrow_root / indiv.pft.cton_root_max)) {
				
				// Maximum possible root growth given N constraints. Put rest of C back into storage (Creturn)
				growth_max_root = navail_root * indiv.pft.cton_root_max;
				creturn += indiv.cgrow_root - growth_max_root;
				indiv.cgrow_root = growth_max_root;
			}

			if (navail_repr < (indiv.cgrow_repr / indiv.pft.cton_leaf_max)) {

				// Maximum possible reproductive growth given N constraints. Put rest of C back into storage (Creturn)
				growth_max_repr = navail_repr * indiv.pft.cton_leaf_max;
				creturn += indiv.cgrow_repr - growth_max_repr;
				indiv.cgrow_repr = growth_max_repr;
			}

			// Make sure storage does not overflow.
			// Maximum C storage pool size 
			if (indiv.cmass_storage + creturn > indiv.cstorage_max) {
				cexcess = indiv.cmass_storage + creturn - indiv.cstorage_max;
				indiv.cexcess += cexcess;
				creturn = indiv.cstorage_max - indiv.cmass_storage;
			} 
			indiv.cmass_storage += creturn;

			indiv.cmass_leafs[green] += indiv.cgrow_leaf;
			indiv.cmass_root += indiv.cgrow_root;
			indiv.cmass_repr += indiv.cgrow_repr;
		
		}
		
		// Always take up available N
		indiv.nmass_leafs[green] += navail_leaf;
		indiv.nmass_root += navail_root;
		indiv.nmass_repr += navail_repr;
		indiv.nmass_storage = max(0.0, indiv.nmass_storage - (navail_leaf + navail_root + navail_repr));

        assert(indiv.cgrow_repr >= 0.0);
		assert(indiv.cgrow_leaf >= 0.0);
		assert(indiv.cgrow_root >= 0.0);
		assert(indiv.cmass_storage >= 0.0);
		//assert(indiv.nmass_storage >= indiv.cmass_storage / indiv.pft.cton_leaf_max);
		assert(indiv.nmass_storage >= 0.0);
		assert(negligible(indiv.cgrow - (indiv.cgrow_repr + indiv.cgrow_leaf + indiv.cgrow_root + creturn + cexcess), -16));
		assert(indiv.cton_leaf() <= indiv.pft.cton_leaf_max);
		assert(indiv.cton_leaf() >= indiv.pft.cton_leaf_min);
		assert(indiv.cton_root() <= indiv.pft.cton_root_max);

		// Save yearly leaf and root growth
		indiv.cmass_leaf_ygrowth += indiv.cgrow_leaf;
		indiv.cmass_root_ygrowth += indiv.cgrow_root;
		
		vegetation.nextobj();
	}
}


// Transfer of C and N from green (live) leaf pool to brown (dead) leaf pool
// Accounts for retranslocation of C and N to storage prior to senescence. 
void senescence(Vegetation& vegetation) {
	
	// Parameters
	const double crelocfrac_c3 = 0.1;  // TODO: move to constants?
	const double crelocfrac_c4 = 0.05; 
	const double nrelocfrac = 0.5;
	const double k = 10.0;            // Shape parameter in senescence function (steepness)
	const double growth_stage0 = 0.5; // Shape parameter (mid point)
	const double qwater_senesc = 1.0;

	double senesc = 0.0, senesc_base = 0.0, senesc_fwater = 0.0, senesc_ftemp = 0.0, senesc_max = 0.0;
	double csenesc = 0.0, cstore = 0.0, ctobrown = 0.0;
	double nsenesc = 0.0, nstore = 0.0, ntobrown = 0.0;
	double nstorage_max = 0.0;

	// Air temperature or bottom snow layer if snow present
	double grass_temp = vegetation.patch.soil.T_soil[6];

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		// Senescence

		// Baseline (minimum) senescence rate from leaf longevity under unstressed conditions.
		senesc_base = 1.0 / (indiv.pft.leaflong * date.year_length());
		
		// Stress factors increase senescence 
		senesc_fwater = pow((1 - ppft.wscal),qwater_senesc);
		
		senesc_ftemp  = 0.0;  // TODO: replace Temp response
		if (grass_temp < 10.0) {
			senesc_ftemp = (10.0 - grass_temp) / 20.0;
		}
		else if (grass_temp > 35.0) {
			senesc_ftemp = (35.0 + grass_temp) / 20.0;
		}

		// Maximum senescence rate
		senesc_max = senesc_base + (senesc_fwater + senesc_ftemp) * senesc_base;

		// Total senescence rate
		senesc = senesc_max / (1 + exp( -k * (min(1.0,indiv.growth_stage - growth_stage0)) ));
		
		// Retranslocation and pool updates
		// TODO: maybe call cstorage_capacity once at the beginning to be consistent across routines!
		// TODO: calculate retranslocation coefficients based on current C:N ratios?

		// TODO: move fractions to parameters??
		double crelocfrac = indiv.pft.pathway == C3 ? crelocfrac_c3 : crelocfrac_c4;

		// Carbon
		csenesc  = indiv.cmass_leafs[green] * senesc;
		cstore   = csenesc * crelocfrac;
		ctobrown = csenesc * (1.0 - crelocfrac);

		// If C storage pool overflows, re-adjust cstore and ctobrown
		if (indiv.cmass_storage + cstore > indiv.cstorage_max) {
			cstore = indiv.cstorage_max - indiv.cmass_storage;
			ctobrown = csenesc - cstore;
		}

		indiv.cmass_leafs[green] -= csenesc;   // green leaf pool
		indiv.cmass_leafs[brown] += ctobrown;  // brown leaf pool
		indiv.cmass_storage  += cstore;    // C storage


		// Nitrogen
		nsenesc  = indiv.nmass_leafs[green] * senesc;
		nstore   = nsenesc * nrelocfrac;
		ntobrown = nsenesc * (1.0 - nrelocfrac);

		// If N storage pool overflows, re-adjust nstore and ntobrown
		if (indiv.nmass_storage + nstore > indiv.nstorage_max) {
			nstore = indiv.nstorage_max - indiv.nmass_storage;
			ntobrown = nsenesc - nstore;
		}

		// TODO: make sure N:C ratio in storage does not fall outside its ranges
		//       understand what the N:C in storage looks like (values)

		indiv.nmass_leafs[green] -= nsenesc;   // green leaf pool
		indiv.nmass_leafs[brown] += ntobrown;  // brown leaf pool
		indiv.nmass_storage  += nstore;        // N storage

		assert(indiv.cmass_leafs[green] >= 0.0);
		assert(indiv.nmass_leafs[green] >= 0.0);
		assert(indiv.cmass_leafs[brown] >= 0.0);
		assert(indiv.nmass_leafs[brown] >= 0.0);

		assert(indiv.cmass_storage >= 0.0);
		assert(indiv.nmass_storage >= 0.0);

		assert(negligible(csenesc - (cstore + ctobrown), -16));
		assert(negligible(nsenesc - (ntobrown + nstore), -16));

		vegetation.nextobj();
	}
}


// Transfer of brown leaf pool and fine root pool to leaf and root litter.
// Transfer reproductive tissue to litter if end of growth cycle is reached.
void turnover_grass(Vegetation& vegetation) {

    // TODO: move to constants
	const double nrelocfrac = 0.5;

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patch& patch = vegetation.patch;
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		// Determine turnover rate (same for all individuals in a PFT).
		// Assumed daily rates here (fraction d-1)!
		// TODO: turnover rates should be moisture and temperature dependent
		double frac_turnover_leaf = indiv.pft.turnover_leaf / date.year_length();
		double frac_turnover_root = indiv.pft.turnover_root / date.year_length();

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

			assert(negligible(nturnover - (nstore + ntolitter), -16));
		}

		// Move reproductive C and N to litter at the end of the growing season.
		if (indiv.growth_stage >= 1.0 && indiv.cmass_repr > 0.0) {
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
		
		vegetation.nextobj();
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
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
	
		// Handle excessive C 
		if (indiv.cexcess > 0.0){
			cstorage_overflow(indiv);
		}

		// Determine in how far C storage induces a negative feedback on photosynthetic capacity.
		sink_source_feedback(indiv);

		vegetation.nextobj();
	}
}


// Updates LAI and FPC for grasses.
// Equivalent to default allometry() function at the moment,
// but might develop. (e.g. more than 1 indiv for grasses, separate green and brown lai, 
// specific extinction coefficient etc.)
void grass_allometry(Vegetation& vegetation) {

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		allometry(indiv);

		vegetation.nextobj();
	}
}


// Aggregate fluxes for annual output.
void annual_accounting(Vegetation& vegetation) {

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

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

		vegetation.nextobj();
	}
}


void growth_daily_grass_dave(Patch& patch) {

	// Obtain reference to Vegetation and Stand object for this patch
	Vegetation& vegetation = patch.vegetation;
	Stand& stand = patch.stand;

	// Record mass balance at start of day.
	DaveMassBalance mass_balance;
	mass_balance.init(vegetation);

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
