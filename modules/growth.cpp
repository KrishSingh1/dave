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
#include <assert.h>


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

const double APHEN_MAX = 210.0;
// Maximum number of equivalent days with full leaf cover per growing season
// for summergreen PFTs

///////////////////////////////////////////////////////////////////////////////////////
// LEAF PHENOLOGY
// Call function leaf_phenology each simulation day prior to calculation of FPAR, to
// calculate fractional leaf-out for each PFT and individual.
// Function leaf_phenology_pft is not intended to be called directly by the framework,


void leaf_phenology_pft(Pft& pft, Climate& climate, double wscal, double aphen,
		double& phen, Patch& patch) {

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
		else if (pft.lifeform == GRASS) {

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

	// Loop through patch-PFTs

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& ppft = patch.pft.getobj();

		// For this PFT ...
		if(patch.stand.pft[ppft.id].active) {
			if(ppft.pft.landcover == CROPLAND && patch.stand.landcover == CROPLAND)
				leaf_phenology_crop(ppft.pft, patch);
			else	//natural, urban, pasture, forest and peatland stands/pft:s
				leaf_phenology_pft(ppft.pft, climate, ppft.wscal, ppft.aphen, ppft.phen, patch);

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
double calc_nrelocfrac_daily(double cton_min, double cmass_storage, double nmass_storage) {
	
	// Current N storage level
	double cton_current = cmass_storage / nmass_storage;

	// If nearing maximum N storage level
	if (1.0 / cton_current > 1.0 / cton_min * 0.7) 
		return max(0.0, nrelocfrac * (1.0 - (1.0 / cton_current - 1.0 / cton_min * 0.7) / (1.0 / cton_min * (1.0 - 0.7))));
	else
		return nrelocfrac;
}

///////////////////////////////////////////////////////////////////////////////////////
// TURNOVER
// Internal function (do not call directly from framework)

void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
		lifeformtype lifeform, landcovertype landcover, double& cmass_leaf, double& cmass_root, double& cmass_sap,
		double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
		double& nmass_heart, double& litter_leaf, double& litter_root,
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

	turnover = turnover_leaf * nmass_leaf;
	nmass_leaf -= turnover;
	nmass_litter_leaf += turnover * (1.0 - actual_nrelocfrac);
	longterm_nstore += turnover * actual_nrelocfrac;

	// Root turnover
	turnover = turnover_root * cmass_root;
	cmass_root -= turnover;
	if (alive) litter_root += turnover;

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
	cmass_debt_inc	= 0.0;

	if (!largerthanzero(ltor, -10)) {

		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_root_inc = bminc;

		// Make sure we don't end up with negative cmass_root or cmass_root_sg

		if(cmass_root_inc < -cmass_root){
			exceeds_cmass += -(cmass_root_inc + cmass_root);
			cmass_root_inc = -cmass_root;
		}

		if (lifeform==TREE) {
			cmass_sap_inc=-cmass_sap;
			cmass_heart_inc=-cmass_sap_inc;
		}
	}
	else if (lifeform==TREE) {

		// TREE ALLOCATION

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

				// If biomass of roots and leafs can't meet biomass decrease then
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
	else if (lifeform == GRASS) {

		// GRASS ALLOCATION
		// Allocation attempts to distribute biomass increment (bminc) among leaf
		// and root compartments, i.e.
		//   (14) bminc = cmass_leaf_inc + cmass_root_inc
		// while satisfying Eqn(3)

		cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
		cmass_root_inc = bminc - cmass_leaf_inc;

		if (bminc >= 0.0) {

			// Positive biomass increment (and bminc == 0.0 for dailygrass)

			if (cmass_leaf_inc < 0.0) {

				// Positive bminc, but ltor causes negative allocation to leaves,
				// put all of bminc into roots and storage

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

		 indiv.alive = true;
		 indiv.report_flux(Fluxes::ESTC, -indiv.cmass_storage);
		 indiv.report_flux(Fluxes::ESTN, -indiv.nmass_storage);
	}
	else {

		allocation(bminit, indiv, ltor, cmass_leaf_ind, cmass_root_ind, cmass_sap_ind, dval, dval, dval, dval, dval);

		indiv.cmass_leaf = cmass_leaf_ind * indiv.densindiv;
		indiv.cmass_root = cmass_root_ind * indiv.densindiv;

		indiv.nmass_leaf = 0.0;
		indiv.nmass_root = 0.0;
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
	//   (1) height = cmass_sap / (sapwood xs area)
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
	// Crown area may be derived from stem diameter by the relation (Zeide 1993):
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

	double diam; // stem diameter (m)
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
			diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);

			// Stem volume
			double vol = indiv.height * PI * diam * diam * 0.25;

			if (indiv.age && (indiv.cmass_heart + indiv.cmass_sap) / indiv.densindiv / vol < indiv.pft.wooddens * 0.9) {
				return false;
			}
		}
		else {
			indiv.height = 0.0;
			diam = 0.0;
			return false;
		}


		// guess2008 - extra height check
		if (indiv.height > HEIGHT_MAX) {
			indiv.height = 0.0;
			diam = 0.0;
			return false;
		}


		// Crown area (Eqn 6)
		indiv.crownarea = min(indiv.pft.k_allom1 * pow(diam, indiv.pft.k_rp),
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
	else if (indiv.pft.lifeform == GRASS) {

		// GRASSES

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
	else if (indiv.pft.lifeform==GRASS) { // grass

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
		patch.fluxes.report_flux(Fluxes::REPRC, pft.litter_repr);
		pft.litter_repr = 0.0;

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

	// Set forest management intensity for this year
	patch.man_strength = cut_fraction(patch);

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		bool killed = false;

		// Standard yearly growth functionalities

		if (!ifdailygrass || indiv.pft.lifeform != GRASS || indiv.istruecrop_or_intercropgrass()) {

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
							(indiv.pft.leaflong * date.year_length()) * (indiv.cmass_leaf + indiv.cmass_root) -
							indiv.cmass_leaf - indiv.cmass_root, 0.0);

					if (cmass_excess > bminc) cmass_excess = bminc;

					// Transfer excess leaves to litter
					// only for 'alive' individuals
					if (indiv.alive) {
						patch.pft[indiv.pft.id].litter_leaf += cmass_excess;
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
								patch.pft[indiv.pft.id].litter_leaf,
								patch.pft[indiv.pft.id].litter_root,
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
						patch.pft[indiv.pft.id].litter_repr += cmass_repr;
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
							patch.pft[indiv.pft.id].litter_leaf += litter_leaf_inc * indiv.densindiv;
							patch.pft[indiv.pft.id].litter_root += litter_root_inc * indiv.densindiv;

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

							indiv.kill();

							vegetation.killobj();
							killed = true;
						}
					}
					else if (indiv.pft.lifeform == GRASS) {

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

								patch.pft[indiv.pft.id].litter_leaf += litter_leaf_inc;
								patch.pft[indiv.pft.id].litter_root += litter_root_inc;

								// C litter exceeding existing biomass
								indiv.report_flux(Fluxes::NPP, exceeds_cmass * indiv.densindiv);
								indiv.report_flux(Fluxes::RA, -exceeds_cmass * indiv.densindiv);
							}

							// Nitrogen always return to soil litter and storage
							// Leaf
							if (indiv.nmass_leaf > 0.0){
								patch.pft[indiv.pft.id].nmass_litter_leaf += litter_leaf_inc * indiv.densindiv /
										cton_leaf_bg * (1.0 - nrelocfrac);
								indiv.nstore_longterm += litter_leaf_inc * indiv.densindiv / cton_leaf_bg * nrelocfrac;
								// Subtracting litter nitrogen from individuals
								indiv.nmass_leaf -= min(indiv.nmass_leaf, litter_leaf_inc * indiv.densindiv / cton_leaf_bg);
							}

							// Root
							if (indiv.nmass_root > 0.0){
								patch.pft[indiv.pft.id].nmass_litter_root += litter_root_inc * indiv.densindiv /
										cton_root_bg * (1.0 - nrelocfrac);
								indiv.nstore_longterm += litter_root_inc / cton_root_bg * nrelocfrac;

								// Subtracting litter nitrogen from individuals
								indiv.nmass_root -= min(indiv.nmass_root, litter_root_inc * indiv.densindiv / cton_root_bg);
							}
						}
						// Kill individual and transfer biomass to litter if either biomass
						// compartment negative

						if (((indiv.cmass_leaf < MINCMASS || indiv.cmass_root < MINCMASS)) && !indiv.istruecrop_or_intercropgrass() ) {

							indiv.kill();

							vegetation.killobj();
							killed = true;
						}
					}

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
							indiv.max_n_storage = max(0.0, min(indiv.cmass_sap * indiv.pft.fnstorage / cton_leaf_bg, retransn_nextyear));
						else // GRASS
							indiv.max_n_storage = max(0.0, min(indiv.cmass_root * indiv.pft.fnstorage / cton_leaf_bg, retransn_nextyear));

						// Scale this year productivity to max storage
						if (indiv.anpp > 0.0) {
							indiv.scale_n_storage = max(indiv.max_n_storage * 0.1, indiv.max_n_storage - retransn_nextyear) * cton_leaf_bg / indiv.anpp;
						} // else use last years scaling factor
					}
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
	const double MAX_GROWTH_RATE = 0.05;
	// maximum gpp trend growth modifier
	const double GPP_FAC_MAX = 0.1;
	// storage to biomass ratio growth modifier max
	const double SB_FAC_MAX = 1.0;
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

	// Air temperature. Should be snowpack temperature if grass is covered in snow.
	double air_temp = patch.get_climate().temp;

	if (date.day == 0) {
		if (!patch.id)
			for (int p=0; p < npft; p++)
				stand.pft[p].cmass_repr = 0.0;
	}

	// Loop through individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

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
				patch.pft[indiv.pft.id].litter_repr += cmass_repr;
				flush_litter_repr(patch);
			}

			// Check C budget
			double initc = bminc + indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_storage +
				indiv.cmass_root + patch.pft[indiv.pft.id].litter_leaf + patch.pft[indiv.pft.id].litter_root;

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

			// Calculate leaf to root ratio based on N and water limitation and current leaf to root status 
			// or only water limitation and current leaf to root status if N limitation is turned off

			// Current leaf to root status
			double ltor_status = (!negligible(indiv.get_cmass_leaf()) && !negligible(indiv.cmass_root)) ? indiv.pft.ltor_max / (indiv.get_cmass_leaf() / indiv.cmass_root) : 1.0;

			if (ifnlim && date.year > freenyears) {
				if (!negligible(indiv.get_cmass_leaf())) {
					indiv.nscal_running.add(min(1.0, indiv.cton_leaf_opt / indiv.cton_leaf()));
				}
				
				// leaf to root ratio of new growth (Johnson 2008 eqn 4.13b)
				indiv.ltor = min(indiv.wscal_running.mean(), indiv.nscal_running.mean()) * indiv.pft.ltor_max * ltor_status;
			}
			else {
				indiv.ltor = indiv.wscal_running.mean() * indiv.pft.ltor_max * ltor_status;
			}
			
			// Growth rate modifiers

			// GPP trend growth modifier 
			// Negative gpp trend less growth more to storage
			double gpp_trend = indiv.gpptolai_15.trend();
			double gpp_fac = gpp_trend > 0.0 ? GPP_FAC_MAX : GPP_FAC_MAX * 0.5;

			// Storage to biomass growth modifier 
			// Lower storage fraction less growth
			double sb_ratio = !negligible(indiv.get_cmass_leaf() + indiv.cmass_root) ? min(5.0, indiv.cmass_storage / (indiv.get_cmass_leaf() + indiv.cmass_root)) : SB_FAC_MAX;
			double sb_fac = sb_ratio > 1.0 && gpp_trend > 0.0 ? sb_ratio : min(1.0, 1.0 - (SB_FAC_MAX - min(SB_FAC_MAX, sb_ratio)) * (SB_FAC_MAX - min(SB_FAC_MAX, sb_ratio)) / (SB_FAC_MAX * SB_FAC_MAX));

			// Fraction of leaves that die during one day due to frost (eqn 2 and 3 Höglind 2016, error in paper eqn)
			// Temperature where half the leaves would be killed. Should be flexible due to hardening, but here set to a fixed value
			// also should use snow temperature if snowpack covers grass.
			double lethal_temp_50 = -15.0;
			double frost_sen = air_temp < 0.0 ? 1.0 - max(0.5, 1.0 / (1.0 + exp(lethal_temp_50 - air_temp))) : 0.0;

			// Water stress growth modifier 
			double wscal_fac = (patch.soil.temp > 0.0 && patch.pft[indiv.pft.id].wscal_5.mean() > 0.2) ? 
				(patch.pft[indiv.pft.id].wscal_5.mean() - 0.2) / 0.8 : 0.0;
			
			// Water limitation increases senescence
			double wscal_sen = 1.0 + (1.0 - wscal_fac);

			// Temperature growth modifier
			double opt_temp = (indiv.pft.pstemp_high + indiv.pft.pstemp_low) / 2.0;
			double min_temp = max(0.0, indiv.pft.pstemp_min);
			double tscal_fac = max(0.0, (air_temp - min_temp) / (opt_temp - min_temp));
			
			// Temperature senescence modifier
			double tscal_sen = max(0.1, tscal_fac);

			// Growth and senescence rates (adopted from eqn 24 Johnson & Thornley 1983, Höglind 2016)
			g_fac = min(MAX_GROWTH_RATE, gpp_fac * sb_fac * wscal_fac * tscal_fac);
			s_fac = max(frost_sen, indiv.pft.transferconst * wscal_sen * tscal_sen);
			
			// Storage used for new leaf and root growth and updating old leaves and roots
			if (!negligible(indiv.cmass_storage) && !negligible(indiv.nmass_storage)) {

				// Todays growth
				double cgrow = g_fac * indiv.cmass_storage;

				// Leaf and root C growth
				cgrow_root = cgrow / (1.0 + indiv.ltor);
				cgrow_leaf += cgrow - cgrow_root;

				// N demand for growth
				double ngrow_leaf_demand = cgrow_leaf / indiv.cton_leaf_opt;
				double ngrow_root_demand = cgrow_root / (indiv.cton_leaf_opt * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);

				// N demand for update current tissue to optimal value (10% day-1)
				double nup_leaf_demand = max(0.0, indiv.get_cmass_leaf() / indiv.cton_leaf_opt - indiv.get_nmass_leaf()) * 0.1;
				double nup_root_demand = max(0.0, indiv.cmass_root / (indiv.cton_leaf_opt * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr) - indiv.nmass_root) * 0.1;

				// Total N demand
				double ndemand_tot = ngrow_leaf_demand + ngrow_root_demand + nup_leaf_demand + nup_root_demand;
				
				// Growth and Update N fraction of storage
				double guN_fac = ndemand_tot / indiv.nmass_storage;

				// Fration of growth and update N demand meet
				double fndemand = 0.0;
				double fupndemand = 0.0;

				// Only N growth or update if there is C growth
				if (g_fac > 0.0) {

					// Fraction of N storage used is less that C storage used -> no N limitation
					if (guN_fac <= g_fac) {
						fndemand = 1.0;
						fupndemand = 1.0;
					}
					// N demand higher than total storage pool -> strong N limitation
					else if (guN_fac > 1.0) {
						fndemand = 0.05 * indiv.nmass_storage / ndemand_tot;
						fupndemand = fndemand;
					}
					else {
						// N limitation determined from difference between C storage fraction to N storage fraction used
						double fnlim = 0.9;
						fndemand = 0.05 + 0.95 * pow(fnlim, guN_fac / g_fac - 1.0);
						fupndemand = fndemand;
					}
				}

				// Leaf and root N growth and update
				ngrow_leaf = ngrow_leaf_demand * fndemand;
				ngrow_root = ngrow_root_demand * fndemand + nup_root_demand * fupndemand;
				nup_leaf = nup_leaf_demand * fupndemand;

				// Leaf N mass update according to compartment size
				if (!negligible(indiv.get_cmass_leaf())) {
					for (int lf = 0; lf < indiv.cmass_leafs.size() - 1; lf++) {
						double fnup = indiv.cmass_leafs[lf] / indiv.get_cmass_leaf();
						indiv.nmass_leafs[lf] += nup_leaf * fnup;
					}
				}

				// Reduce C depending on fndemand to not get too high C:N ratios
				if (ifnlim && date.year > freenyears && !negligible(fndemand - 1.0)) {
					cgrow_leaf *= min(1.0, 1.0 - (1.0 - fndemand) * FRAC_C_OF_N_LIM);
					cgrow_root *= min(1.0, 1.0 - (1.0 - fndemand) * FRAC_C_OF_N_LIM);

					// New tissue can't have a C:N ratio higher than allowed PFT specific maximum
					if (cgrow_leaf / ngrow_leaf > indiv.pft.cton_leaf_max) {
						cgrow_leaf = ngrow_leaf * indiv.pft.cton_leaf_max;
						cgrow_root = ngrow_root * indiv.pft.cton_leaf_max * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr;
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
			nrelocfrac_daily = calc_nrelocfrac_daily(indiv.cton_leaf_aopt, indiv.cmass_storage, indiv.nmass_storage);

			// Movement variables for C and N (eqn 26a-d Johnson & Thornley 1983)
			c1 = min(1.0, LAMBDA * s_fac) * indiv.cmass_leafs[0];
			n1 = min(1.0, LAMBDA * s_fac) * indiv.nmass_leafs[0];
			c2 = s_fac * indiv.cmass_leafs[1];
			n2 = s_fac * indiv.nmass_leafs[1];
			c3 = s_fac * indiv.cmass_leafs[2];
			n3 = s_fac * indiv.nmass_leafs[2];
			c4 = s_fac * indiv.cmass_leafs[3];
			n4 = s_fac * indiv.nmass_leafs[3];

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
				patch.pft[indiv.pft.id].litter_leaf += c4;
				patch.pft[indiv.pft.id].nmass_litter_leaf += n4;

				// Daily leaf and root turnover
				double total_cmass_leaf = indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + c4;
				double leaves_turnover_frac = min(turnover_max, !negligible(total_cmass_leaf) ? c4 / total_cmass_leaf : 0.0);

				// Root C turnover (use same fraction as for leaves)
				double turnov = leaves_turnover_frac * indiv.cmass_root;
				indiv.cmass_root -= turnov;
				patch.pft[indiv.pft.id].litter_root += turnov;

				// Root N turnover
				turnov = leaves_turnover_frac * indiv.nmass_root;
				indiv.nmass_root -= turnov;
				patch.pft[indiv.pft.id].nmass_litter_root += turnov * (1.0 - nrelocfrac_daily);
				indiv.nmass_storage += turnov * nrelocfrac_daily;
			}

			// C budget check
			double finalc = indiv.get_cmass_leaf() + indiv.cmass_leafs[3] + indiv.cmass_storage +
				indiv.cmass_root + patch.pft[indiv.pft.id].litter_leaf + patch.pft[indiv.pft.id].litter_root;

			if (fabs(initc - finalc) > 1.0e-12)
				dprintf("Year %d Day %d Diff %g Init %g Final %g\n",date.year,date.day, initc - finalc, initc, finalc);

			// Update allometry - updates LAI and FPC based on cmass
			allometry(indiv);

			// Determine the yearly maximum values which is reported
			indiv.lai_ymax = max(indiv.lai_ymax, indiv.lai);

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
		growth_daily_grass(patch);

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
// Zeide, B (1993) Primary unit of the tree crown. Ecology 74: 1598-1602.
