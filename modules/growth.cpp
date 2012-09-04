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
	double& phen) {

	// DESCRIPTION
	// Calculates leaf phenological status (fractional leaf-out) for a individuals of
	// a given PFT, given current heat sum and length of chilling period (summergreen
	// PFTs) and water stress coefficient (raingreen PFTs)

	// INPUT PARAMETER
	// wscal = water stress coefficient (0-1; 1=maximum stress)
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

	// guess2008
	bool leafout = true; // CHILLDAYS

	// Obtain reference to Vegetation object
	Vegetation& vegetation = patch.vegetation;

	// INDIVIDUAL AND COHORT MODES
	// Calculate phenology for each PFT at this patch

	// Loop through patch-PFTs

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft = patch.pft.getobj();

		// For this PFT ...
		if(patch.stand.pft[pft.id].active)
		{
			leaf_phenology_pft(pft.pft, climate, pft.wscal, pft.aphen, pft.phen);

			if (pft.pft.lifeform == TREE && (pft.pft.phenology == SUMMERGREEN || pft.pft.phenology == ANY))
				if (pft.phen < 1.0) leafout = false; // CHILLDAYS
		}
		// Update annual leaf-on sum
		if (climate.lat >= 0.0 && date.day == COLDEST_DAY_NHEMISPHERE ||
			climate.lat < 0.0 && date.day == COLDEST_DAY_SHEMISPHERE) pft.aphen = 0.0;
		pft.aphen += pft.phen;

		// ... on to next PFT
		patch.pft.nextobj();
	}


	if (leafout) {
		climate.ifsensechill = true; // CHILLDAYS
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


///////////////////////////////////////////////////////////////////////////////////////
// TURNOVER
// Internal function (do not call directly from framework)

void turnover(double turnover_leaf, double turnover_root, double turnover_sap,
	lifeformtype lifeform, double& cmass_leaf, double& cmass_root, double& cmass_sap,
	double& cmass_heart, double& nmass_leaf, double& nmass_root, double& nmass_sap,
	double& nmass_heart, double& litter_leaf, double& litter_root,
	double& nmass_litter_leaf, double& nmass_litter_root,
	double& nstore, Fluxes& fluxes, bool alive, double& nmass_avail,
	landcovertype landcover, Gridcell& gridcell) {

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

	// OUTPUT PARAMETERS
	// litter_leaf   = new leaf litter (kgC/m2)
	// litter_root   = new root litter (kgC/m2)
	// cmass_heart   = heartwood C biomass (kgC/m2)

	double turnover = 0.0;
	double scale = 1.0;

	if(run_landcover && gridcell.LC_updated) {
		//scale harvest products of stands with increased area by (old area/new area) if landcover change has occurred:
		scale = gridcell.landcoverfrac_old[landcover] / gridcell.landcoverfrac[landcover];

		if(scale >= 1.0)
			scale = 1.0;
	}

	// TREES AND GRASSES:

	// Leaf turnover
	turnover = turnover_leaf * cmass_leaf * scale;
	cmass_leaf -= turnover;
	if (alive) litter_leaf += turnover;
	
	turnover = turnover_leaf * nmass_leaf * scale;
	nmass_leaf -= turnover;
	if (alive) {
		nmass_litter_leaf += turnover * (1.0 - nrelocfrac);
		nstore += turnover * nrelocfrac;
	}
	else {
		nmass_avail += turnover * (1.0 - nrelocfrac);	
		nstore += turnover * nrelocfrac;
	}

	// Root turnover
	turnover = turnover_root * cmass_root * scale;
	cmass_root -= turnover;
	if (alive) litter_root += turnover;

	turnover = turnover_root * nmass_root * scale;
	nmass_root -= turnover;
	if (alive) {
		nmass_litter_root += turnover * (1.0 - nrelocfrac);
		nstore += turnover * nrelocfrac;
	}
	else {
		nmass_avail += turnover * (1.0 - nrelocfrac);
		nstore += turnover * nrelocfrac;
	}

	if (lifeform == TREE) {

		// TREES ONLY:

		// Sapwood turnover by conversion to heartwood
		turnover = turnover_sap * cmass_sap * scale;
		cmass_sap -= turnover;
		cmass_heart += turnover;

		// NB: assumes N is translocated from sapwood prior to conversion to
		//     heartwood and that this is the same fraction that is conserved
		//     in conjunction with leaf and root shedding
		
		turnover = turnover_sap * nmass_sap * scale;
		nmass_sap -= turnover;
		nmass_heart += turnover * (1.0 - nrelocfrac);
		nstore += turnover * nrelocfrac;
	}	
	
}


///////////////////////////////////////////////////////////////////////////////////////
// REPRODUCTION
// Internal function (do not call directly from framework)

void reproduction(double reprfrac, double reprCN, double npp, double& bminc, double& cmass_repr, double& nmass_repr, double& nstore) {

	// DESCRIPTION
	// Allocation of net primary production (NPP) to reproduction and calculation of
	// assimilated carbon available for production of new biomass

	// INPUT PARAMETERS
	// reprfrac = fraction of NPP for this time period allocated to reproduction
	// reprCN   = reproduction C:N ratio
	// npp      = NPP (i.e. assimilation minus maintenance and growth respiration) for
	//            this time period (kgC/m2)

	// OUTPUT PARAMETER
	// bminc    = carbon biomass increment (component of NPP available for production
	//            of new biomass) for this time period (kgC/m2)

	if (npp >= 0.0) {
		cmass_repr = npp * reprfrac;
		bminc = npp - cmass_repr;

		/*nmass_repr=cmass_repr/reprCN;
		if(nmass_repr>nstore)
			nmass_repr=nstore;

		nstore-=nmass_repr;*/
		return;
	}

	// Negative NPP - no reproduction cost

	cmass_repr = 0.0;
	nmass_repr = 0.0;
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


void allocation(double bminc, double cmass_leaf, double cmass_root, double cmass_sap,
	double cmass_debt, double cmass_heart, double ltor, double height, double sla,
	double wooddens, lifeformtype lifeform, double k_latosa, double k_allom2,
	double k_allom3, double& cmass_leaf_inc, double& cmass_root_inc,
	double& cmass_sap_inc, double& cmass_debt_inc, double& cmass_heart_inc,
	double& litter_leaf_inc, double& litter_root_inc) {

	// DESCRIPTION
	// Calculates changes in C compartment sizes (leaves, roots, sapwood, heartwood)
	// and litter for a plant individual as a result of allocation of biomass increment.
	// Assumed allometric relationships are given in function allometry below.

	// INPUT PARAMETERS
	// bminc       = biomass increment this time period on individual basis (kgC)
	// cmass_leaf  = leaf C biomass for last time period on individual basis (kgC)
	// cmass_root  = root C biomass for last time period on individual basis (kgC)
	// cmass_sap   = sapwood C biomass for last time period on individual basis (kgC)
	// cmass_heart = heartwood C biomass for last time period on individual basis (kgC)
	// ltor        = leaf to root mass ratio following allocation
	// height      = individual height (m)
	// sla         = specific leaf area (PFT-specific constant) (m2/kgC)
	// wooddens    = wood density (PFT-specific constant) (kgC/m3)
	// lifeform    = life form class (TREE or GRASS)
	// k_latosa    = ratio of leaf area to sapwood cross-sectional area (PFT-specific
	//               constant)
	// k_allom2    = constant in allometry equations
	// k_allom3    = constant in allometry equations

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

	const int NSEG = 20; // number of segments (parameter in numerical methods)
	const int JMAX = 40; // maximum number of iterations (in numerical methods)
	const double XACC = 0.0001; // threshold x-axis precision of allocation solution
	const double YACC = 1.0e-10; // threshold y-axis precision of allocation solution
	const double CDEBT_MAXLOAN_DEFICIT = 0.8; // maximum loan as a fraction of deficit
	const double CDEBT_MAXLOAN_MASS = 0.2; // maximum loan as a fraction of (sapwood-cdebt)

	double cmass_leaf_inc_min;
	double cmass_root_inc_min;
	double x1, x2, dx, xmid, fx1, fmid, rtbis, sign;
	int j;
	double cmass_deficit, cmass_loan;

	litter_leaf_inc = 0.0;
	litter_root_inc = 0.0;
	cmass_root_inc = 0.0; // guess2008 - initialise 

	if (ltor < 1.0e-10) {
		
		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_leaf_inc = 0.0;
		cmass_root_inc = bminc;

		if (lifeform == TREE) {
			cmass_sap_inc =- cmass_sap;
			cmass_heart_inc =- cmass_sap_inc;
		}

		dprintf("Year %d 111 ltor %g No leaf production possible\n", date.year, ltor);

		return;
	}

	if (lifeform == TREE) {

		// TREE ALLOCATION

		cmass_heart_inc = 0.0;

		// Calculate minimum leaf increment to maintain current sapwood biomass
		// Given Eqn (2)

		if (height > 0.0)
			cmass_leaf_inc_min = k_latosa * cmass_sap / (wooddens * height * sla) - cmass_leaf;
		else
			cmass_leaf_inc_min = 0.0;

		// Calculate minimum root increment to support minimum resulting leaf biomass
		// Eqn (3)

		if (height > 0.0)
			cmass_root_inc_min = k_latosa * cmass_sap / (wooddens * height * sla * ltor) -
				cmass_root;
		else
			cmass_root_inc_min = 0.0;

		if (cmass_root_inc_min < 0.0) { // some roots would have to be killed

			cmass_leaf_inc_min = cmass_root * ltor - cmass_leaf;
			cmass_root_inc_min = 0.0;
		}

		// BLARP! C debt stuff
		if (ifcdebt) {
			cmass_deficit = cmass_leaf_inc_min + cmass_root_inc_min - bminc;
			if (cmass_deficit > 0.0) {
				cmass_loan = max(min(cmass_deficit * CDEBT_MAXLOAN_DEFICIT,
					(cmass_sap - cmass_debt) * CDEBT_MAXLOAN_MASS), 0.0);
				bminc += cmass_loan;
				cmass_debt_inc = cmass_loan;
			}
			else cmass_debt_inc = 0.0;
		}
		else cmass_debt_inc = 0.0;

		if (cmass_root_inc_min >= 0.0 && cmass_leaf_inc_min >= 0.0 &&
			cmass_root_inc_min + cmass_leaf_inc_min <= bminc) {

			// Normal allocation (positive increment to all living C compartments)
			// NOTE: includes allocation of zero or negative NPP, c.f. LPJF

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

			if (f(x1) >= 0.0) sign =- 1.0;
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
		}
		else {

			// Abnormal allocation: negativ bminc

			if (bminc < 0) {
				
				cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
				cmass_root_inc = bminc - cmass_leaf_inc;

				if (cmass_leaf_inc > 0.0) {
					cmass_leaf_inc = 0.0;
					cmass_root_inc = bminc;
				}
				else if (cmass_root_inc > 0.0) {
					cmass_root_inc = 0.0;
					cmass_leaf_inc = bminc;
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
			}

			// Make sure we don't end up with negative cmass_leaf
			cmass_leaf_inc = max(-cmass_leaf, cmass_leaf_inc);

			// Make sure we don't end up with negative cmass_root
			cmass_root_inc = max(-cmass_root, cmass_root_inc);

			// Add killed leaves to litter
			litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

			// Add killed roots to litter
			litter_root_inc = max(-cmass_root_inc, 0.0);

			// Calculate increase in sapwood mass (which must be negative)
			// Eqn (2)

			cmass_sap_inc = (cmass_leaf_inc + cmass_leaf) * wooddens * height * sla / k_latosa -
				cmass_sap;

			// Convert killed sapwood to heartwood
			if (cmass_sap_inc < 0.0)  
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

		if (cmass_leaf_inc < 0.0 && bminc > 0.0) {

			// Positive bminc, but ltor causes negative allocation to leaves,
			// put all of bminc into roots

			cmass_root_inc = bminc;
			cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf; // Eqn (3)
		}
		else if (cmass_root_inc < 0.0 && bminc > 0.0) {

			// Positive bminc, but ltor causes negative allocation to roots,
			// put all of bminc into leaves

			cmass_leaf_inc = bminc;
			cmass_root_inc = (cmass_leaf + bminc) / ltor - cmass_root;
		}
		else if (bminc < 0) {
			if (cmass_leaf_inc > 0.0) {
				cmass_leaf_inc = 0.0;
				cmass_root_inc = bminc;
			}

			if (cmass_root_inc > 0.0) {
				cmass_root_inc = 0.0;
				cmass_leaf_inc = bminc;
			}
		}

		// Make sure we don't end up with negative cmass_leaf
		cmass_leaf_inc = max(-cmass_leaf, cmass_leaf_inc);

		// Make sure we don't end up with negative cmass_root
		cmass_root_inc = max(-cmass_root, cmass_root_inc);

		// Add killed leaves to litter
		litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

		// Add killed roots to litter
		litter_root_inc = max(-cmass_root_inc, 0.0);
	}
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

	allocation(bminit, 0.0, 0.0, 0.0, 0.0, 0.0, ltor, 0.0, indiv.pft.sla, indiv.pft.wooddens,
		indiv.pft.lifeform, indiv.pft.k_latosa, indiv.pft.k_allom2, indiv.pft.k_allom3,
		cmass_leaf_ind, cmass_root_ind, cmass_sap_ind, dval, dval, dval, dval);

	indiv.cmass_leaf = cmass_leaf_ind * indiv.densindiv;
	indiv.cmass_root = cmass_root_ind * indiv.densindiv;

	indiv.nmass_leaf = indiv.cmass_leaf / indiv.pft.cton_leaf_avr;
	indiv.nmass_root = indiv.cmass_root / indiv.pft.cton_root_avr;
		
	if (indiv.pft.lifeform == TREE) {
		indiv.cmass_sap = cmass_sap_ind * indiv.densindiv;
		indiv.nmass_sap = indiv.cmass_sap / indiv.pft.cton_sap_avr; 

		indiv.cton_growth = (indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap) /
			(indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_sap);
	}
	else {
		indiv.cton_growth = (indiv.cmass_leaf + indiv.cmass_root) /
			(indiv.nmass_leaf + indiv.nmass_root);
	}
}

// As raingreen has a portion of its ANPP which is certain, this portion has priority to nitrogen
void raingreen_ndemand(Vegetation& vegetation, double& ndemand_patch, double& nsupply_patch, double& fnuptake_patch) {

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (indiv.pft.phenology == RAINGREEN) {
			nsupply_patch -= indiv.raingreen_ndemand;
			ndemand_patch -= indiv.raingreen_ndemand;
			indiv.ndemand -= indiv.raingreen_ndemand;
			indiv.nstore += indiv.raingreen_ndemand;
		}
		vegetation.nextobj();
	} 

	fnuptake_patch = nsupply_patch / ndemand_patch;

}

/// Calculates individual nitrogen demand
/** Calculates this year's N demand by doing a fake growth with tissue turnover
 *  and allocation of fixed carbon to reproduction and new biomass, and by doing 
 *  so determining the N demand. 
 *  Accumulated NPP (assimilation minus maintenance and growth respiration) on
 *  patch or modelled area basis assumed to be given by 'anpp' member variable for
 *  each individual.
 */
void indiv_ndemand(Gridcell& gridcell, landcovertype landcover, Fluxes& fluxes, Individual& indiv) {

	double CDEBT_PAYBACK_RATE = 0.2;
	double cmass_payback;
	double cmass_excess;

	double cmass_leaf_tmp = indiv.cmass_leaf;
	double cmass_root_tmp = indiv.cmass_root;
	double cmass_sap_tmp = indiv.cmass_sap;
	double cmass_heart_tmp = indiv.cmass_heart;
	double cmass_debt_tmp = indiv.cmass_debt;

	double nmass_leaf_tmp = indiv.nmass_leaf;
	double nmass_root_tmp = indiv.nmass_root;
	double nmass_sap_tmp = indiv.nmass_sap;
	double nmass_heart_tmp = indiv.nmass_heart;

	double densindiv_tmp = indiv.densindiv;
 
	double nstore_turnover = 0.0;	// N retranslocated in turnover
	indiv.raingreen_ndemand = 0.0;
	//double repr_ndemand;

	// Thomas this year's ndemand: declaration
	double bminc;
	double dval;
	Fluxes fluxes_tmp = fluxes;

	// C:N ratios
	double cton_root_new, cton_sap_new;

	// Increases
	double cmass_leaf_inc, cmass_root_inc, cmass_sap_inc, cmass_heart_inc;
	double cmass_debt_inc;

	double nscal;
	double ltor;

	// reproduction
	if (indiv.anpp >= 0.0) {
		bminc = indiv.anpp * (1.0 - indiv.pft.reprfrac);
		// repr_ndemand = indiv.anpp * indiv.pft.reprfrac / indiv.pft.reprCN;
	}
	else {
		bminc = indiv.anpp;
		// repr_ndemand = 0.0;
	}	

	// Determine N demand of extra raingreen leaf
	if (bminc >= 0 && indiv.pft.phenology == RAINGREEN) {

		// Raingreen PFTs: reduce biomass increment to account for NPP
		// allocated to extra leaves during the past year.
		// Excess allocation to leaves given by:
		//   aphen_raingreen / ( leaf_longevity * 365) * cmass_leaf -
		//   cmass_leaf

		// BLARP! excess allocation to roots now also included (assumes leaf longevity = root longevity)

		cmass_excess = max((double)indiv.aphen_raingreen /
			(indiv.pft.leaflong * 365.0) * (cmass_leaf_tmp + cmass_root_tmp) -
			cmass_leaf_tmp - cmass_root_tmp, 0.0);

		if (cmass_excess > bminc) cmass_excess = bminc;

		// Deduct from this year's C biomass increment
		if (indiv.alive) {
			bminc -= cmass_excess;
			indiv.raingreen_ndemand = cmass_excess / (cmass_leaf_tmp / nmass_leaf_tmp);
		}
	}

	// N stress scalar for leaf to root allocation (based on Zaehle 2010 SM eq 19)
	nscal = min(1.0, indiv.pft.cton_leaf_avr / (cmass_leaf_tmp / nmass_leaf_tmp));

	// Set leaf:root mass ratio based on water stress parameter 
	// or N stress scalar 
	ltor = min(indiv.wscal_mean, nscal) * indiv.pft.ltor_max;

	// C:N ratio for new and current biomass		
	if (ifnlim && date.year > freenyears) {

		if (!negligible(indiv.leafn_mean) && !negligible(cmass_leaf_tmp))
			indiv.cton_leaf_new = cmass_leaf_tmp / indiv.leafn_mean; 
				// actual mean leafN based on Vmax from photosynthesis
		else
			indiv.cton_leaf_new = indiv.pft.cton_leaf_avr;

		if (date.year <= freenyears)
			indiv.cton_leaf_new = indiv.pft.cton_leaf_avr;

		indiv.cton_leaf_opt = indiv.cton_leaf_new;

		// C:N ratio can't be outside of pft min max range
		indiv.cton_leaf_new = min(indiv.pft.cton_leaf_max, max(indiv.pft.cton_leaf_min, indiv.cton_leaf_new));

		cton_root_new=
			indiv.cton_leaf_new * (indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);
		cton_sap_new=
			indiv.cton_leaf_new * (indiv.pft.cton_sap_avr / indiv.pft.cton_leaf_avr);
	}
	else {
		indiv.cton_leaf_new = indiv.pft.cton_leaf_avr;
		cton_root_new = indiv.pft.cton_root_avr;
		cton_sap_new = indiv.pft.cton_sap_avr;
	}

	turnover(indiv.pft.turnover_leaf, indiv.pft.turnover_root, indiv.pft.turnover_sap, indiv.pft.lifeform,
			cmass_leaf_tmp, cmass_root_tmp, cmass_sap_tmp, cmass_heart_tmp,
			nmass_leaf_tmp, nmass_root_tmp, nmass_sap_tmp, nmass_heart_tmp,
			dval, dval, dval, dval, nstore_turnover, fluxes_tmp, indiv.alive,
			dval, landcover, gridcell);

	// N demand not associated with growth
	indiv.ndemand = indiv.raingreen_ndemand - nstore_turnover;

	if (indiv.pft.lifeform == TREE) { 

		// pay back part of cdebt
		if (ifcdebt && bminc > 0.0) {
			cmass_payback = min(cmass_debt_tmp * CDEBT_PAYBACK_RATE, bminc);
			bminc -= cmass_payback;
			cmass_debt_tmp -= cmass_payback;
		}

		allocation(bminc / densindiv_tmp, cmass_leaf_tmp / densindiv_tmp,
				cmass_root_tmp / densindiv_tmp, cmass_sap_tmp / densindiv_tmp,
				cmass_debt_tmp / densindiv_tmp, cmass_heart_tmp / densindiv_tmp,
					ltor, indiv.height, indiv.pft.sla, indiv.pft.wooddens, TREE,
					indiv.pft.k_latosa, indiv.pft.k_allom2, indiv.pft.k_allom3,
					cmass_leaf_inc, cmass_root_inc, cmass_sap_inc, cmass_debt_inc,
					cmass_heart_inc, dval, dval);

		// Calculate N needed for this new biomass
		indiv.ndemand += 
			max(0.0, cmass_leaf_inc) * densindiv_tmp / indiv.cton_leaf_new+
			max(0.0, cmass_root_inc) * densindiv_tmp / cton_root_new+
			max(0.0, cmass_sap_inc) * densindiv_tmp / cton_sap_new;
		
		// Calculate fractions of new biomass that goes to leaf and root
		if (cmass_leaf_inc > 0.0 && cmass_root_inc > 0.0) {
			indiv.bminc_leaf_frac = cmass_leaf_inc / (cmass_leaf_inc + cmass_root_inc + max(0.0, cmass_sap_inc));
			indiv.bminc_root_frac = cmass_root_inc / (cmass_leaf_inc + cmass_root_inc + max(0.0, cmass_sap_inc));
		}

		// No negative N demand
		indiv.ndemand = max(indiv.ndemand, 0.0);
	}
	else {

		allocation(bminc, cmass_leaf_tmp, cmass_root_tmp,
					0.0, 0.0, 0.0, ltor, 0.0, 0.0, 0.0, GRASS, 0.0,
					0.0, 0.0, cmass_leaf_inc, cmass_root_inc, dval, dval, dval,
					dval, dval);

		// Calculate N needed for this new biomass
		indiv.ndemand += 
			max(0.0, cmass_leaf_inc) * densindiv_tmp / indiv.cton_leaf_new +
			max(0.0, cmass_root_inc) * densindiv_tmp / cton_root_new;

		// Calculate fractions of new biomass that goes to leaf and root
		if (!negligible(max(0.0, cmass_leaf_inc)) && !negligible(max(0.0, cmass_root_inc))) {
			indiv.bminc_leaf_frac = max(0.0, max(0.0, cmass_leaf_inc) / (max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc)));
			indiv.bminc_root_frac = max(0.0, max(0.0, cmass_root_inc) / (max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc)));
		}

		// No negative N demand
		indiv.ndemand = max(indiv.ndemand, 0.0); 
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// ALLOCATION_NLIM
// Nitrogen-limited allocation

double f_nlim(double& cmass_leaf_inc,
	double nstore, double cmass_leaf, double cmass_root, double cmass_sap, double cmass_heart,
	double cton_leaf, double cton_root, double cton_sap, double ltor,
	double wooddens, double sla, double k_allom2, double k_allom3, double k_latosa,
	Pft& pft, int place, double x2, int count) {

	// Returns value of f(cmass_leaf_inc), given by:

	// f(cmass_leaf_inc) = 0 =
	//          ( 4 * (cmass_sap+(nstore-cmass_leaf_inc/cton_leaf-((cmass_leaf+cmass_leaf_inc)/ltor-cmass_root)/cton_root)*
	//				cton_sap+cmass_heart) / wooddens / PI / k_allom2 ) ^ 1/(k_allom3+2) 
	//			-
	//			((cmass_sap+(nstore-cmass_leaf_inc/cton_leaf-((cmass_leaf+cmass_leaf_inc)/ltor-cmass_root)/cton_root)*cton_sap) /
	//				(cmass_leaf+cmass_leaf_inc) / k_allom2 / sla / wooddens  * k_latosa ) ^ (1/k_allom3)
	//
	// See function allocation_nlim (below), Eqn (15)

	const double PI = 3.1415926536;
	double op0, op00, op1, op2, result;

	// Validate arguments
	if (negligible(cton_leaf))
		fail("Year %d f_nlim for %s at %d: cton_leaf=%g count=%d", date.year, (char*)pft.name, place, cton_leaf, count);
	
	if (negligible(ltor)) 
		fail("Year %d f_nlim for %s at %d: ltor=%g count=%d", date.year, (char*)pft.name, place, ltor, count);

	if (negligible(cmass_leaf + cmass_leaf_inc)) 
		fail("Year %d f_nlim for %s at %d: (cmass_leaf+cmass_leaf_inc)=%g cmass_leaf=%g cmass_leaf_inc=%g count=%d",
			date.year, (char*)pft.name, place,
			cmass_leaf + cmass_leaf_inc, cmass_leaf, cmass_leaf_inc, count);

	// op0 is difference between available N and proposed N allocation to leaves+roots
	op00 = cmass_leaf_inc / cton_leaf + ((cmass_leaf + cmass_leaf_inc) / ltor - cmass_root) / cton_root;
	op0 = nstore - op00;

	if (op0 < -1.0e-12) 
		fail("Year %d f_nlim for %s at %d: op0=%g count=%d\n", date.year, (char*)pft.name, place, op0, count);

	op1 = (cmass_sap + op0 * cton_sap + cmass_heart) / wooddens / PI / k_allom2;	

	op2 = (cmass_sap + op0 * cton_sap) / (cmass_leaf + cmass_leaf_inc) / k_allom2 / sla / wooddens * k_latosa;

	// Validate partial expressions
	if (op1 < 0.0) 
		fail("f_nlim for %s at %d: op1=%g C:N sap %g cmass_heart %g count=%d\n", (char*)pft.name, place, op1, cton_sap, cmass_heart, count);

	if (op2 < 0.0) 
		fail("f_nlim for %s at %d: op2=%g count=%d\n", (char*)pft.name, place, op2, count);

	result = pow(4.0 * op1, 1.0 / (k_allom3 + 2.0)) - pow(op2, 1.0 / k_allom3);

	return result;
}

/// Allocation with N constraint
/** Calculates changes in C compartment sizes (leaves, roots, sapwood, heartwood), biomass
 *  increment and litter for a plant individual, given the prescribed C:N ratios for the
 *  living biomass compartments (leaves, roots, sapwood) and the total amount of N available
 *  in the living compartments plus 'storage'.
 *  This version of allocation must only be called when standard allocation fails due to
 *  insufficient N availability, as the C biomass increment is predicted, not prescribed.
 */
void allocation_nlim(Pft& pft, double nstore, double cton_leaf, double cton_root, double cton_sap,
	double bminc, double cmass_leaf, double cmass_root, double cmass_sap, double cmass_heart,
	double ltor, double height,
	double& cmass_leaf_inc, double& cmass_root_inc, double& cmass_sap_inc,
	double& cmass_heart_inc, double& litter_leaf_inc, double& litter_root_inc,
	double densindiv) { 

	// Assumed allometric relationships are given in function allometry below.

	// INPUT PARAMETERS
	// (all C and N masses are on individual, not area, basis)
	// bminc         = MAXIMUM biomass increment this time period (kgC)
	// cmass_leaf    = current leaf C biomass (kgC)
	// cmass_root    = current root C biomass (kgC)
	// cmass_sap     = current sapwood C biomass (kgC)
	// cmass_heart   = current heartwood C biomass (kgC)
	// cton_leaf     = target C:N for new leaves following allocation
	// cton_root     = target C:N for new roots following allocation
	// cton_sap      = target C:N for new sapwood following allocation
	// ltor          = leaf to root mass ratio following allocation
	// height        = current individual height (m)
	// pft = Pft object including the members:
	//   wooddens    = wood density (PFT-specific constant) (kgC/m3)
	//   k_allom2    = constant in allometry equations
	//   k_allom3    = constant in allometry equations
	//   k_latosa    = ratio of leaf area to sapwood cross-sectional area (PFT-specific
	//                 constant)

	// OUTPUT PARAMETERS
	// bminc           = biomass increment this time period on individual basis (kgC)
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

	// MATHEMATICAL DERIVATION FOR TREE ALLOCATION

	// Available nitrogen provides overall constraint on structure following allocation:
	//   (1) nstore = cmass_leaf_inc/cton_leaf + 
	//               cmass_root_inc/cton_root +
	//               cmass_sap_inc/cton_sap
	// Prescribed leaf:root mass ratio:
	//   (2) cmass_leaf+cmass_leaf_inc = ltor * (cmass_root+cmass_root_inc)
	// Combining (1) and (2):
	//   (3) nstore = cmass_leaf_inc/cton_leaf + 
	//               ((cmass_leaf+cmass_leaf_inc)/ltor - cmass_root)/cton_root +
	//               cmass_sap_inc/cton_sap
	// Rearranging (3):
	//   (4) cmass_sap_inc = ( nstore - cmass_leaf_inc/cton_leaf -
	//						((cmass_leaf+cmass_leaf_inc)/ltor - cmass_root)/cton_root ) * cton_sap
	//
	// Let V = stem volume following allocation
	//     H = stem height following allocation
	//     D = stem diameter following allocation
	//     PI = ratio of the circumference to the radius of a circle
	// Stems are cylindrical:
	//   (5) V = H * PI * (D/2)^2
	//   (6) V = (cmass_sap+cmass_sap_inc+cmass_heart) / wooddens
	// Combining (5) and (6):
	//   (7) H = 4 * (cmass_sap+cmass_sap_inc+cmass_heart) / wooddens / PI / D^2
	// Allometry equation (5):
	//   (8) H = k_allom2 * D^k_allom3
	// Combining (7) and (8):
	//   (9) D = ( 4 * (cmass_sap+cmass_sap_inc+cmass_heart) / wooddens / PI / k_allom2 ) ^ 1/(k_allom3+2)
	// Allometry equation (2):
	//  (10) (leaf area) = k_latosa * (sapwood xs area)
	//  (11) wooddens = (cmass_sap+cmass_sap_inc) / H / (sapwood xs area) 
	// From (10) and (11):
	//  (12) wooddens = (cmass_sap+cmass_sap_inc) / H / sla / (cmass_leaf+cmass_leaf_inc) * k_latosa
	// Combining (8) and (12):
	//  (13) D = ( (cmass_sap+cmass_sap_inc) / (cmass_leaf+cmass_leaf_inc) /
	//            k_allom2 / sla / wooddens  * k_latosa ) ^ (1/k_allom3)
	// Combining (9) and (13):
	//  (14) ( 4 * (cmass_sap+cmass_sap_inc+cmass_heart) / wooddens / PI / k_allom2 ) ^ 1/(k_allom3+2) -
	//       ( (cmass_sap+cmass_sap_inc) / (cmass_leaf+cmass_leaf_inc) /
	//            k_allom2 / sla / wooddens  * k_latosa ) ^ (1/k_allom3) = 0
	// Substituting (cmass_sap_inc) from (4) into (14):
	//  (15) f(cmass_leaf_inc) = 0 =
	//          ( 4 * (cmass_sap+(nstore-cmass_leaf_inc/cton_leaf-((cmass_leaf+cmass_leaf_inc)/ltor-cmass_root)/cton_root)*
	//				cton_sap+cmass_heart) / wooddens / PI / k_allom2 ) ^ 1/(k_allom3+2) -
	//					((cmass_sap+(nstore-cmass_leaf_inc/cton_leaf-((cmass_leaf+cmass_leaf_inc)/ltor-cmass_root)/cton_root)*
	//						cton_sap) /	(cmass_leaf+cmass_leaf_inc) / k_allom2 / sla / wooddens  * k_latosa ) ^ (1/k_allom3)			
	//
	// Numerical methods are used to solve Eqn (15) for cmass_leaf_inc

	const int NSEG = 20; // number of segments (parameter in numerical methods)
	const int JMAX = 40; // maximum number of iterations (in numerical methods)
	const double XACC = 0.0001; // threshold x-axis precision of allocation solution
	const double YACC = 1.0e-10; // threshold y-axis precision of allocation solution
	const double PI = 3.14159265;

	double cmass_leaf_inc_min;
	double cmass_root_inc_min;
	double x1, x2, dx, xmid, fx1, fmid, rtbis, sign, x2_store;
	int j;

	/// if abnormal allocation because of C limitation (bminc)
	bool ifabnormal_alloc_Climit = false;
	/// if abnormal allocation because of N limitation (bminc_n)
	bool ifabnormal_alloc_Nlimit = false;
	/// dummies
	double op0, op00;
	/// bminc limited by nstore	
	double bminc_n;
	/// bminc limited by carbon
	double bminc_c;

	litter_leaf_inc = 0.0;
	litter_root_inc = 0.0;
	cmass_root_inc = 0.0; // guess2008 - initialise

	if (ltor < 1.0e-10) {
		
		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_leaf_inc = 0.0;
		cmass_root_inc = bminc;

		if (pft.lifeform == TREE) {
			cmass_sap_inc = -cmass_sap;
			cmass_heart_inc = -cmass_sap_inc;
		}

		dprintf("Year %d 222 ltor %g No leaf production possible\n", date.year, ltor);

		return;
	}

	if (pft.lifeform == TREE) {

		// TREE ALLOCATION

		cmass_heart_inc = 0.0;

		// Calculate minimum leaf increment to maintain current sapwood biomass. Given Eqn (2)
		// Calculate minimum root increment to support minimum resulting leaf biomass. Eqn (3)

		if (height > 0.0) {
			cmass_leaf_inc_min = pft.k_latosa * cmass_sap / (pft.wooddens * height * pft.sla) - cmass_leaf;
			cmass_root_inc_min = pft.k_latosa * cmass_sap / (pft.wooddens * height * pft.sla * ltor) - cmass_root;
		}
		else {
			cmass_leaf_inc_min = 0.0;
			cmass_root_inc_min = 0.0;
		}			

		if (cmass_root_inc_min < 0.0) { // some roots would have to be killed
			cmass_leaf_inc_min = cmass_root * ltor - cmass_leaf;
			cmass_root_inc_min = 0.0;
		}

		// Check if bminc and nstore sufficient for normal allocation,

		// First bminc, i.e. carbon limitation
		if (cmass_root_inc_min >= 0.0 && cmass_leaf_inc_min >= 0.0 &&
			cmass_root_inc_min + cmass_leaf_inc_min <= bminc) // Note that this includes negative bminc!
			ifabnormal_alloc_Climit = false;
		else 
			ifabnormal_alloc_Climit = true;

		// Now nstore, i.e. N limitation
		// op0 is difference between available N and nesseceary N for minimum leaf and root increment
		op00 = cmass_leaf_inc_min / cton_leaf + cmass_root_inc_min / cton_root;
		op0 = nstore - op00;	
		
		if (op0 < 0.0) {
			ifabnormal_alloc_Nlimit = true;
			// Calculate bminc as limited by nmass
			bminc_n = nstore * cton_leaf;		
		}

		if (!ifabnormal_alloc_Climit && !ifabnormal_alloc_Nlimit) {

			if (cmass_root_inc_min >= 0.0 && cmass_leaf_inc_min >= 0.0) {
				
				// Normal allocation (positive increment to all living C compartments)

				// allocation fix make sure while ends before exceeding max
				int count = 0;	// number of iterations so far

				// Calculation of leaf mass increment (lminc_ind) satisfying Eqn (13)
				// using bisection method (Press et al 1986)

				x1 = 0.0;

				// Maximum bound for cmass_leaf_inc is given by allocation of all available
				// N to leaves and roots (the actual bound would be lower than this because
				// of sapwood allocation, but this is the limit for numerical stability in f_nlim)
			
				// nstore = nmass_leaf_inc + nmass_root_inc
				//       = cmass_leaf_inc/cton_leaf + cmass_leaf_inc/ltor/cton_root	
				// cmass_leaf + cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor
				// Rearranging:
				//	cmass_leaf_inc = ( (cmass_root + nstore * cton_leaf)*ltor - cmass_leaf) / (ltor + 1)

				x2 = ((cmass_root + nstore * cton_leaf) * ltor - cmass_leaf) / (ltor + 1);

				x2_store = x2;

				dx = (x2 - x1) / (double)NSEG;

				if (cmass_leaf < 1.0e-10) {
					x1 += dx; // to avoid division by zero
					count++;
				}

				// Evaluate f(x1), i.e. Eqn (15) at cmass_leaf_inc = x1

				fx1 = f_nlim(x1, nstore, cmass_leaf, cmass_root, cmass_sap, cmass_heart, cton_leaf, cton_root, cton_sap,
					ltor, pft.wooddens, pft.sla, pft.k_allom2, pft.k_allom3, pft.k_latosa,
					pft, 1, x2_store, 0);

				// Find approximate location of leftmost root on the interval
				// (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
				// change in sign of f_nlim(xmid) relative to f_nlim(x1).

				fmid = f_nlim(x1, nstore, cmass_leaf, cmass_root, cmass_sap, cmass_heart, cton_leaf, cton_root, cton_sap,
					ltor, pft.wooddens, pft.sla, pft.k_allom2, pft.k_allom3, pft.k_latosa,
					pft, 2, x2_store, 0);

				xmid = x1;

				while (fmid * fx1 > 0.0 && xmid < x2 && count < NSEG) {
					count++;
					xmid += dx;
					fmid = f_nlim(xmid, nstore, cmass_leaf, cmass_root, cmass_sap, cmass_heart, cton_leaf, cton_root, cton_sap,
						ltor, pft.wooddens, pft.sla, pft.k_allom2, pft.k_allom3, pft.k_latosa,
						pft, 3, x2_store, count);
				}

				x1 = xmid - dx;
				x2 = xmid;

				// Apply bisection to find root on new interval (x1,x2)

				if (f_nlim(x1, nstore, cmass_leaf, cmass_root, cmass_sap, cmass_heart, cton_leaf, cton_root, cton_sap,
					ltor, pft.wooddens, pft.sla, pft.k_allom2, pft.k_allom3, pft.k_latosa, pft, 4, x2_store, 0) >= 0.0)
						sign = -1.0;
				else sign = 1.0;

				rtbis = x1;
				dx = x2 - x1;

				// Bisection loop
				// Search iterates on value of xmid until xmid lies within
				// xacc of the root, i.e. until |xmid-x|<xacc where f_nlim(x)=0

				fmid = 1.0; // dummy value to guarantee entry into loop
				j = 0; // number of iterations so far

				while (dx >= XACC && fabs(fmid) > YACC && j <= JMAX) {

					dx *= 0.5;
					xmid = rtbis + dx;

					fmid = f_nlim(xmid, nstore, cmass_leaf, cmass_root, cmass_sap, cmass_heart, cton_leaf, cton_root, cton_sap,
						ltor, pft.wooddens, pft.sla, pft.k_allom2, pft.k_allom3, pft.k_latosa, pft, 5, x2_store, j);

					if (fmid * sign <= 0.0) rtbis = xmid;
					j++;
				}

				// Now rtbis contains numerical solution for cmass_leaf_inc given Eqn (15)

				cmass_leaf_inc = rtbis;

				// Calculate increments in other compartments

				cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
				
				// Rearranging Eqn (1):
				// cmass_sap_inc = ( nstore -
				//                   cmass_leaf_inc/cton_leaf -
				//                   cmass_root_inc/cton_root ) * cton_sap

				// But also
				// bminc >= cmass_leaf_inc+cmass_root_inc+cmass_sap_inc

				cmass_sap_inc = (nstore - cmass_leaf_inc / cton_leaf - cmass_root_inc / cton_root ) * cton_sap;

				if (cmass_leaf_inc + cmass_root_inc + cmass_sap_inc > bminc) {

					// C increment higher than bminc
					//
					// diff = cmass_leaf_inc+cmass_root_inc+cmass_sap_inc-bminc
					//
					// A = leaf_inc_new + root_inc_new
					//
					// A + sap_inc_new + diff = cmass_leaf_inc+cmass_root_inc+cmass_sap_inc
					//
					// A/cton_leaf_new + B/cton_sap = (cmass_leaf_inc+cmass_root_inc)/cton_leaf + cmass_sap_inc/cton_sap
					//
					// sap_inc_new = ( (cmass_leaf_inc+cmass_root_inc)/cton_leaf + cmass_sap_inc/cton_sap -
					//		(cmass_leaf_inc+cmass_root_inc+cmass_sap_inc-diff)/cton_leaf ) * (1/cton_sap-1/cton_leaf)
					//
					// A = cmass_leaf_inc+cmass_root_inc+cmass_sap_inc - (sap_inc_new + diff)
					//
					// cmass_leaf+leaf_inc_new = ltor * (cmass_root+A-leaf_inc_new)
					//
					// leaf_inc_new = ( ltor*(cmass_root+A) - cmass_leaf )/(1+ltor)

					double diff = cmass_leaf_inc + cmass_root_inc + cmass_sap_inc-bminc;

					cmass_sap_inc = ((cmass_leaf_inc + cmass_root_inc) / cton_leaf + cmass_sap_inc / cton_sap -
							(cmass_leaf_inc + cmass_root_inc + cmass_sap_inc - diff) / cton_leaf) / (1.0 / cton_sap - 1.0 / cton_leaf);

					double A = cmass_leaf_inc + cmass_root_inc + cmass_sap_inc - (cmass_sap_inc + diff);

					cmass_leaf_inc = (ltor * (cmass_root + A) - cmass_leaf ) / (1 + ltor);

					cmass_root_inc = A - cmass_leaf_inc;
				}
				
				if (cmass_sap_inc < 0.0)
					cmass_heart_inc = -cmass_sap_inc;
			}
			else fail("allocation_nlim: cmass_leaf_inc_min=%g cmass_root_inc_min=%g",
				cmass_leaf_inc_min, cmass_root_inc_min);

		} 
		else {
			// Now carry out abnormal allocation under N limitation,
			if (ifabnormal_alloc_Nlimit && !ifabnormal_alloc_Climit) {

				// Abnormal allocation: reduction in some biomass compartment(s) to
				// satisfy allometry

				// Attempt to distribute this year's production among leaves and roots only

				if (bminc_n < 0.0) {

					// Try and divid up the negative bminc_n between leaf and root
					//
					// bminc_n = cmass_leaf_inc + cmass_root_inc
					//
					// cmass_leaf + cmass_leaf_inc = ltor (cmass_root + cmass_root_inc)
					//
					//	->
					//
					// cmass_root_inc = (bminc_n + cmass_leaf - ltor*cmass_root)/(ltor + 1.0)
					//
					// cmass_leaf_inc = bminc_n - cmass_root_inc

					cmass_root_inc = (bminc_n + cmass_leaf - ltor * cmass_root) / (ltor + 1.0);

					cmass_leaf_inc = bminc_n - cmass_root_inc;

					// compartment increment can't be positive
					if (cmass_leaf_inc > 0.0) {
						cmass_leaf_inc = 0.0;
						cmass_root_inc = bminc_n;
					}

					if (cmass_root_inc > 0.0) {
						cmass_root_inc = 0.0;
						cmass_leaf_inc = bminc_n;
					}
				}
				else {

					cmass_leaf_inc = (bminc_n - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);

					if (cmass_leaf_inc > 0.0) {

						// Positive allocation to leaves

						cmass_root_inc = bminc_n - cmass_leaf_inc; // Eqn (1)

						if (cmass_root_inc < 0.0) {
							cmass_leaf_inc = bminc_n;
							cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
						}
					}
					else {

						// Negative or zero allocation to leaves
						// Eqns (1), (3)

						cmass_root_inc = bminc_n;
						cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf;
					}
				}
			}
			// Now carry out abnormal allocation under C limitation, 
			else if (!ifabnormal_alloc_Nlimit && ifabnormal_alloc_Climit) {

				// Attempt to distribute this year's production with C:N ratio as leaves and roots to
				// prevent over allocation

				if (bminc > 0.0 && bminc / cton_leaf > nstore)
					bminc_c = nstore * cton_leaf;
				else 
					bminc_c = bminc;

				if (bminc_c >= 0.0) {

					// Abnormal allocation: reduction in some biomass compartment(s) to
					// satisfy allometry

					// Attempt to distribute this year's production among leaves and roots only

					cmass_leaf_inc = (bminc_c - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
	
					if (cmass_leaf_inc > 0.0) {

						// Positive allocation to leaves

						cmass_root_inc = bminc_c - cmass_leaf_inc; // Eqn (1)

						// Add killed roots (if any) to litter

						if (cmass_root_inc < 0.0) {
							cmass_leaf_inc = bminc_c;
							cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
						}
					}
					else {

						// Negative or zero allocation to leaves
						// Eqns (1), (3)

						cmass_root_inc = bminc_c;
						cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf;
					}
				}
				else {	// Negative bminc_c

					cmass_leaf_inc = (bminc_c - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
	
					cmass_root_inc = bminc_c - cmass_leaf_inc;

					// compartment increment can't be positive
					if (cmass_leaf_inc > 0.0) {
						cmass_leaf_inc = 0.0;
						cmass_root_inc = bminc_c;
					}

					if (cmass_root_inc > 0.0) {
						cmass_root_inc = 0.0;
						cmass_leaf_inc = bminc_c;
					}
				}
			} 
			// Now carry out abnormal allocation under both C and N limitation, 
			// for now done in same way as under C limitation

			else if (ifabnormal_alloc_Nlimit && ifabnormal_alloc_Climit) {

				if (bminc > 0.0 && bminc / cton_leaf > nstore)
					bminc_c = nstore * cton_leaf;
				else 
					bminc_c = bminc;

				if (bminc_c >= 0.0) {
					
					// Abnormal allocation: reduction in some biomass compartment(s) to
					// satisfy allometry

					// Attempt to distribute this year's production among leaves and roots only

					cmass_leaf_inc = (bminc_c - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
	
					if (cmass_leaf_inc > 0.0) {

						// Positive allocation to leaves

						cmass_root_inc = bminc_c - cmass_leaf_inc; // Eqn (1)

						if (cmass_root_inc < 0.0) {
							cmass_leaf_inc = bminc_c;
							cmass_root_inc = (cmass_leaf_inc + cmass_leaf) / ltor - cmass_root; // Eqn (3)
						}
					}
					else {

						// Negative or zero allocation to leaves
						// Eqns (1), (3)

						cmass_root_inc = bminc_c;
						cmass_leaf_inc = (cmass_root + cmass_root_inc) * ltor - cmass_leaf;
					}
				}
				else {
				
					cmass_leaf_inc = (bminc_c - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
	
					cmass_root_inc = bminc_c - cmass_leaf_inc; 

					// compartment increment can't be positive
					if (cmass_leaf_inc > 0.0) {
						cmass_leaf_inc = 0.0;
						cmass_root_inc = bminc_c;	
					}

					if (cmass_root_inc > 0.0) {
						cmass_root_inc = 0.0;
						cmass_leaf_inc = bminc_c;
					}
				}
			}

			// Make sure we don't end up with negative cmass_leaf
			cmass_leaf_inc = max(-cmass_leaf, cmass_leaf_inc);

			// Make sure we don't end up with negative cmass_root
			cmass_root_inc = max(-cmass_root, cmass_root_inc);

			// Add killed leaves to litter
			litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

			// Add killed roots to litter
			litter_root_inc = max(-cmass_root_inc, 0.0);

			// Calculate increase in sapwood mass
			// Eqn (12)
			cmass_sap_inc = (max(cmass_leaf_inc + cmass_leaf, 0.0)) * pft.wooddens * height * pft.sla / pft.k_latosa -
				cmass_sap;

			// Convert killed sapwood to heartwood
			if (cmass_sap_inc < 0.0)  {
				if (-cmass_sap_inc > cmass_sap) 
					cmass_sap_inc = -cmass_sap;

				cmass_heart_inc = -cmass_sap_inc;
			}

		}
		//end abnormal allocation

		// Check that total increment does not exceed available carbon (should never do so but just to be sure ...)

		if (cmass_leaf_inc + cmass_root_inc + cmass_sap_inc - 1.0e-14 > bminc)
			dprintf("Year %d TREE %s allocation_nlim: total increment (%g) exceeds available carbon (%g) diff (%g) available bm (%g) sap_inc %g\n",
				date.year, (char*)pft.name, cmass_leaf_inc + cmass_root_inc + cmass_sap_inc, bminc,
				cmass_leaf_inc + cmass_root_inc + cmass_sap_inc - bminc, cmass_root + cmass_leaf + cmass_sap, cmass_sap_inc);
	}
	else if (pft.lifeform == GRASS) {

		// GRASS ALLOCATION
		// Allocation attempts to distribute available nitrogen (nmass) among leaf
		// and root compartments, 
		//   (16)  nstore = cmass_leaf_inc/cton_leaf + cmass_root_inc/cton_root
		// while satisfying:
		//   (17)  cmass_leaf+cmass_leaf_inc = ltor * (cmass_root+cmass_root_inc)
		// Combining (16) and (17):
		//   (19)  cmass_leaf_inc = ltor * (cmass_root + nstore*cton_leaf) / (1+ltor)         

		if (bminc >= 0.0) {

			cmass_leaf_inc = ltor * (cmass_root + nstore * cton_leaf) / (1 + ltor) ;
			
			// Eqn (17)			

			cmass_root_inc = (cmass_leaf + cmass_leaf_inc) / ltor - cmass_root;

			if (bminc < cmass_leaf_inc + cmass_root_inc) {

				//   (14) bminc = cmass_leaf_inc + cmass_root_inc
				// while satisfying Eqn(3)

				cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
				cmass_root_inc = bminc - cmass_leaf_inc;		
			}
		}
		else {	
				 
			cmass_leaf_inc = (bminc - cmass_leaf / ltor + cmass_root) / (1.0 + 1.0 / ltor);
			cmass_root_inc = bminc - cmass_leaf_inc;

			// compartment increment can't be positive
			if (cmass_leaf_inc > 0.0) {
				cmass_leaf_inc = 0.0;
				cmass_root_inc = bminc;
			}

			if (cmass_root_inc > 0.0) {
				cmass_root_inc = 0.0;
				cmass_leaf_inc = bminc;
			}
		}

		// Make sure we don't end up with negative cmass_leaf
		cmass_leaf_inc = max(-cmass_leaf, cmass_leaf_inc);

		// Make sure we don't end up with negative cmass_root
		cmass_root_inc = max(-cmass_root, cmass_root_inc);

		// Add killed leaves to litter
		litter_leaf_inc = max(-cmass_leaf_inc, 0.0);

		// Add killed roots to litter
		litter_root_inc = max(-cmass_root_inc, 0.0);

		// Check that total increment does not exceed available carbon (should never do so but just to be sure ...)

		if (cmass_leaf_inc + cmass_root_inc > bminc)
			dprintf("Year %d GRASS %s allocation_nlim: total increment (%g) exceeds available carbon bminc (%g) available bm (%g)\n",
				date.year, (char*)pft.name, cmass_leaf_inc + cmass_root_inc, bminc, cmass_root + cmass_leaf);
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
	double vol; // stem volume (m^2)
	double PI = 3.1415927;
	double fpc_new; // updated FPC

	// guess2008 - max tree height allowed (metre).
	const double HEIGHT_MAX = 150.0; 

	if (indiv.pft.lifeform == TREE) {

		// TREES

		// Height (Eqn 4)

		// guess2008 - new allometry check 
		if (!negligible(indiv.cmass_leaf)) {

			indiv.height = indiv.cmass_sap / indiv.cmass_leaf / indiv.pft.sla * indiv.pft.k_latosa / indiv.pft.wooddens;

			// Stem diameter (Eqn 5)
			diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);

			// Stem volume
			vol = indiv.height * PI * diam * diam * 0.25;

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
				indiv.pft.sla / indiv.crownarea;
			
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
		indiv.lai = indiv.cmass_leaf * indiv.pft.sla;
	}
	else if (indiv.pft.lifeform == GRASS) {
		
		// GRASSES

		// guess2008 - bugfix - added if 
		if (!negligible(indiv.cmass_leaf)) {

			// Grass "individual" LAI (Eqn 11)
			indiv.lai_indiv = indiv.cmass_leaf * indiv.pft.sla;

			// FPC (Eqn 10)
			indiv.fpc = 1.0 - lambertbeer(indiv.lai_indiv);

			// Stand-level LAI
			indiv.lai = indiv.lai_indiv;
		} 
		else {
			return false;
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
		return 1.0+2.0/indiv.cmass_leaf/indiv.pft.sla*
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

	// Returns N-free reproduction "litter" to atmosphere

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft = patch.pft.getobj();
		
		patch.fluxes.acflux_soil += pft.litter_repr;
		patch.fluxes.mcflux_soil[date.month] += pft.litter_repr;
		pft.litter_repr = 0.0;

		for (int d=0; d<365; d++)
			patch.fluxes.dcflux_soil[d] += pft.litter_repr / 365.0;

		//patch.fluxes.aNrepr+=pft.nlitter_repr;

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

	/// minimum carbon mass allowed (kgC/m2)
	const double MINCMASS = 1.0e-8;
	/// maximum carbon mass allowed (kgC/m2)
	const double MAXCMASS = 1.0e8;

	const double CDEBT_PAYBACK_RATE = 0.2;

	const double EPS = 1.0e-15;

	/// carbon biomass increment (component of NPP available for production of
	/// new biomass) for this time period on modelled area basis (kgC/m2)
	double bminc;
	/// C allocated to reproduction this time period on modelled area basis (kgC/m2)	
	double cmass_repr;
	/// N allocated to reproduction this time period on modelled area basis (kgN/m2)	
	double nmass_repr;
	/// increment in leaf C biomass following allocation, on individual basis (kgC)	
	double cmass_leaf_inc;
	/// increment in root C biomass following allocation, on individual basis (kgC)	
	double cmass_root_inc;
	/// increment in sapwood C biomass following allocation, on individual basis (kgC)	
	double cmass_sap_inc;
	/// increment in heartwood C biomass following allocation, on individual basis (kgC)
	double cmass_heart_inc;
	/// increment in heartwood C biomass following allocation, on individual basis (kgC)
	double cmass_debt_inc = 0.0; 
	/// increment in leaf litter following allocation, on individual basis (kgC)	
	double litter_leaf_inc = 0.0;
	/// increment in root litter following allocation, on individual basis (kgC)	
	double litter_root_inc = 0.0;
	/// C biomass of leaves in "excess" of set allocated last year to raingreen PFT last year (kgC/m2)
	double cmass_excess;
	/// RAINGREEN N demand for leaves dropped during the year
	double raingreen_ndemand;
	/// N stress scalar for leaf to root allocation
	double nscal;
	/// C:N ratio of leafs that will result in no N limitation
	double cton_leaf_full_growth;

	double dval;
	int p;
	bool killed;

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation = patch.vegetation;
	Gridcell& gridcell = vegetation.patch.stand.gridcell;

	// On first call to function growth this year (patch #0), initialise stand-PFT
	// record of summed allocation to reproduction

	if (!patch.id)
		for (p=0; p<npft; p++)
			stand.pft[p].cmass_repr = 0.0;

	// Loop through individuals	

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		// For this individual 

		// N stress scalar for leaf to root allocation (adopted from Zaehle 2010 SM eq 19) 		
		if (ifnlim && date.year > freenyears)
			nscal = min(1.0, indiv.pft.cton_leaf_avr / (indiv.cmass_leaf / indiv.nmass_leaf));
		else
			nscal = 1.0;

		// Set leaf:root mass ratio based on water stress parameter 
		// or N stress scalar 
		indiv.ltor = min(indiv.wscal_mean, nscal) * indiv.pft.ltor_max;

		indiv.deltafpc = 0.0;

		killed = false;

		if (negligible(indiv.densindiv))
			fail("growth: negligible densindiv for %s", (char*)indiv.pft.name);
		else {

			// Calculate compartment C:N ratios
			if (!negligible(indiv.nmass_leaf))
				indiv.cton_leaf_old = indiv.cmass_leaf / indiv.nmass_leaf;
			else
				indiv.cton_leaf_old = indiv.pft.cton_leaf_avr;
		
			if (!negligible(indiv.nmass_root))
				indiv.cton_root_old = indiv.cmass_root / indiv.nmass_root;
			else
				indiv.cton_root_old = indiv.pft.cton_root_avr;

			if (!negligible(indiv.nmass_sap))
				indiv.cton_sap_old = indiv.cmass_sap / indiv.nmass_sap;
			else
				indiv.cton_sap_old = indiv.pft.cton_sap_avr;

			// Allocation to reproduction
			reproduction(indiv.pft.reprfrac, indiv.pft.reprCN / indiv.fnuptake, indiv.anpp, bminc, cmass_repr, nmass_repr, indiv.nstore);

			// added bminc check. Otherwise we get -ve litter_leaf for grasses when indiv.anpp < 0.
			if (bminc >= 0 && (indiv.pft.phenology == RAINGREEN || indiv.pft.phenology == ANY)) {

				// Raingreen PFTs: reduce biomass increment to account for NPP
				// allocated to extra leaves during the past year.
				// Excess allocation to leaves given by:
				//   aphen_raingreen / ( leaf_longevity * 365) * cmass_leaf -
				//   cmass_leaf

				// BLARP! excess allocation to roots now also included (assumes leaf longevity = root longevity)

				cmass_excess = max((double)indiv.aphen_raingreen /
					(indiv.pft.leaflong * 365.0) * (indiv.cmass_leaf + indiv.cmass_root) -
					indiv.cmass_leaf - indiv.cmass_root, 0.0);

				if (cmass_excess > bminc) cmass_excess = bminc;

				// Transfer excess leaves to litter
				// only for 'alive' individuals
				if (indiv.alive) {
					patch.pft[indiv.pft.id].litter_leaf += cmass_excess;
					
					raingreen_ndemand = min(indiv.nstore, cmass_excess / indiv.cton_leaf_old);
					patch.pft[indiv.pft.id].nmass_litter_leaf += raingreen_ndemand;
					indiv.nstore -= raingreen_ndemand;
				}

				// Deduct from this year's C biomass increment
				// added alive check
				if (indiv.alive) bminc -= cmass_excess;
			}

			// Tissue turnover and associated litter production

			turnover(indiv.pft.turnover_leaf, indiv.pft.turnover_root,
				indiv.pft.turnover_sap, indiv.pft.lifeform,
				indiv.cmass_leaf, indiv.cmass_root, indiv.cmass_sap, indiv.cmass_heart,
				indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap, indiv.nmass_heart,
				patch.pft[indiv.pft.id].litter_leaf,
				patch.pft[indiv.pft.id].litter_root,
				patch.pft[indiv.pft.id].nmass_litter_leaf,
				patch.pft[indiv.pft.id].nmass_litter_root,
				indiv.nstore, patch.fluxes, indiv.alive, patch.soil.nmass_avail,
				indiv.pft.landcover, gridcell);

			// Update stand record of reproduction by this PFT
			stand.pft[indiv.pft.id].cmass_repr += cmass_repr / (double)stand.npatch();
			stand.pft[indiv.pft.id].nmass_repr += nmass_repr / (double)stand.nobj;

			// Transfer reproduction straight to litter
			// only for 'alive' individuals
			if (indiv.alive) 
				patch.pft[indiv.pft.id].litter_repr += cmass_repr;

			//patch.pft[indiv.pft.id].nlitter_repr+=nmass_repr;

			// C:N ratio for new biomass
			if (ifnlim && date.year > freenyears) {	

				// A simple allocation with fractions of biomass going to leafs, roots and sap 
				// determined from the ndemand allocation without any N limitation (indiv_ndemand())
				// to determine C:N ratio of leafs that will result in no N limitation

				if (bminc > 0.0 && indiv.ndemand > 0.0) {
					double A;
					if (indiv.pft.lifeform == TREE) {
						A = bminc * (indiv.bminc_leaf_frac + indiv.bminc_root_frac * indiv.pft.cton_leaf_avr / indiv.pft.cton_root_avr +
							(1.0 - indiv.bminc_leaf_frac - indiv.bminc_root_frac) * indiv.pft.cton_leaf_avr / indiv.pft.cton_sap_avr);
					}
					else {
						A = bminc * (indiv.bminc_leaf_frac + indiv.bminc_root_frac * indiv.pft.cton_leaf_avr / indiv.pft.cton_root_avr);
					}
					cton_leaf_full_growth = A / indiv.nstore;
				}
				else
					cton_leaf_full_growth = indiv.cton_leaf_new;	

				// Determine C:N of new tissue
				indiv.cton_leaf_new = 1.0 / (1.0 / indiv.cton_leaf_opt - (1.0 / indiv.cton_leaf_opt - 1.0 / cton_leaf_full_growth) * indiv.avmaxnlim);

				// C:N ratio can't be outside of pft min max range and not lower than the optimal value
				indiv.cton_leaf_new = min(indiv.pft.cton_leaf_max, max(indiv.pft.cton_leaf_min, max(indiv.cton_leaf_opt, indiv.cton_leaf_new)));
				
				indiv.cton_root_new =
					indiv.cton_leaf_new * (indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);
		
				indiv.cton_sap_new =
					indiv.cton_leaf_new * (indiv.pft.cton_sap_avr / indiv.pft.cton_leaf_avr);
			}
			else {
				indiv.cton_leaf_new = indiv.pft.cton_leaf_avr;
				indiv.cton_root_new = indiv.pft.cton_root_avr;
				indiv.cton_sap_new = indiv.pft.cton_sap_avr;
				indiv.avmaxnlim = 1.0;
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

				allocation(bminc / indiv.densindiv, indiv.cmass_leaf / indiv.densindiv,
					indiv.cmass_root / indiv.densindiv, indiv.cmass_sap / indiv.densindiv,
					indiv.cmass_debt / indiv.densindiv,
					indiv.cmass_heart / indiv.densindiv, indiv.ltor,
					indiv.height, indiv.pft.sla, indiv.pft.wooddens, TREE,
					indiv.pft.k_latosa, indiv.pft.k_allom2, indiv.pft.k_allom3,
					cmass_leaf_inc, cmass_root_inc, cmass_sap_inc, cmass_debt_inc,
					cmass_heart_inc,
					litter_leaf_inc, litter_root_inc);

				indiv.ndemand =
					max(0.0, cmass_leaf_inc) * indiv.densindiv / indiv.cton_leaf_new +
					max(0.0, cmass_root_inc) * indiv.densindiv / indiv.cton_root_new +
					max(0.0, cmass_sap_inc) * indiv.densindiv / indiv.cton_sap_new;

				indiv.ndemand_no_nlim = indiv.ndemand;
				
				// Compute limitation factor based on balance between individual N demand and supply
				// (NB: this overwrites the alternative factor calculated in canexch.cpp, but this
				// one is better!)

				if (ifnlim && (indiv.nstore + EPS < indiv.ndemand || indiv.nstore < 0.0) && indiv.ndemand > 0.0) { 
					if (indiv.nstore < 0.0) {
						indiv.limnfact = 0.0;
					}
					else {
						// Use N from long-term storage to reduce N limitation
						double diff = indiv.ndemand - indiv.nstore;
						indiv.nstore += diff * min(indiv.nmass_reserve / diff, 2.0) / 2.0;
						indiv.nmass_reserve -= diff * min(indiv.nmass_reserve / diff, 2.0) / 2.0;

						indiv.limnfact = indiv.nstore / indiv.ndemand;
					}
				}
				else 
					indiv.limnfact = 1.0;

				// Nitrogen limitation of production
				if (ifnlim && indiv.limnfact < 1.0 && date.year > freenyears) {

					double bminc_nlim, bminc_dec, gpp_dec, agpp;

					// C debt
					cmass_debt_inc = 0.0;

					// Redo allocation, this time with N constraint	
					allocation_nlim(indiv.pft, indiv.nstore / indiv.densindiv, indiv.cton_leaf_new, indiv.cton_root_new,	
						indiv.cton_sap_new, bminc / indiv.densindiv, indiv.cmass_leaf / indiv.densindiv,
						indiv.cmass_root / indiv.densindiv, indiv.cmass_sap / indiv.densindiv,
						indiv.cmass_heart / indiv.densindiv, indiv.ltor, indiv.height,
						cmass_leaf_inc, cmass_root_inc, cmass_sap_inc, cmass_heart_inc,
						litter_leaf_inc, litter_root_inc, indiv.densindiv); 

					// Update N demand
					indiv.ndemand =
						max(0.0, cmass_leaf_inc) * indiv.densindiv / indiv.cton_leaf_new +
						max(0.0, cmass_root_inc) * indiv.densindiv / indiv.cton_root_new +
						max(0.0, cmass_sap_inc) * indiv.densindiv / indiv.cton_sap_new;
					
					// Calculate new biomass increment
					bminc_nlim = max(0.0, cmass_leaf_inc) * indiv.densindiv +
							   max(0.0, cmass_root_inc) * indiv.densindiv +
							   max(0.0, cmass_sap_inc) * indiv.densindiv; 

					// Update accumulate C fluxes according to N limitation
					bminc_dec = bminc - bminc_nlim;
					agpp = 0.0;

					for (int mon=0; mon<12; mon++) 
						agpp += indiv.mgpp[mon];

					for (int month=0; month<12; month++) {
						gpp_dec = indiv.mgpp[month] * bminc_dec / agpp;
						indiv.mgpp[month] -= gpp_dec;
						indiv.mnpp[month] -= gpp_dec;
						if (indiv.alive) {
							patch.fluxes.mcflux_gpp[month] -= gpp_dec;
							patch.fluxes.acflux_veg += gpp_dec;
						}
					}

					// Update annual npp
					// NB: this is important because it affects growth efficiency and
					//     therefore mortality and litter fluxes (in vegdynam.cpp).
					bminc -= bminc_dec;
					indiv.anpp -= bminc_dec;
					indiv.frac_agpp = (agpp - bminc_dec) / agpp;
				}

				// Subtract used nitrogen for biomass increment from storage pool
				if (ifnlim && date.year > freenyears)
					indiv.nstore -= indiv.ndemand;
				else 
					indiv.nstore = 0.0;

				// Optimal growth C:N
				if (!negligible(max(0.0,cmass_leaf_inc)+max(0.0,cmass_root_inc)+max(0.0,cmass_sap_inc)))
					indiv.cton_growth = (max(0.0, cmass_leaf_inc) * max(indiv.pft.cton_leaf_min, indiv.cton_leaf_opt) +
						max(0.0, cmass_root_inc) * max(indiv.pft.cton_leaf_min, indiv.cton_leaf_opt) *
						(indiv.cton_root_new / indiv.cton_leaf_new) +
						max(0.0, cmass_sap_inc) * max(indiv.pft.cton_leaf_min,indiv.cton_leaf_opt) * 
						(indiv.cton_sap_new / indiv.cton_leaf_new)) /
						(max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc) + max(0.0, cmass_sap_inc));

				// Update carbon pools and litter (on area basis)
				// (litter not accrued for not 'alive' individuals - Ben 2007-11-28)

				// Leaves
				indiv.cmass_leaf += cmass_leaf_inc * indiv.densindiv;
				indiv.nmass_leaf += cmass_leaf_inc * indiv.densindiv / ((cmass_leaf_inc > 0) ? indiv.cton_leaf_new : indiv.cton_leaf_old);

				// Roots
				indiv.cmass_root += cmass_root_inc * indiv.densindiv;
				indiv.nmass_root += cmass_root_inc * indiv.densindiv / ((cmass_root_inc > 0) ? indiv.cton_root_new : indiv.cton_root_old);

				// Sapwood
				indiv.cmass_sap += cmass_sap_inc * indiv.densindiv;
				double nmass_sap_inc = cmass_sap_inc * indiv.densindiv / ((cmass_sap_inc > 0) ? indiv.cton_sap_new : indiv.cton_sap_old);
				indiv.nmass_sap += nmass_sap_inc;

				// Heartwood
				indiv.cmass_heart += cmass_heart_inc * indiv.densindiv;
				indiv.nmass_heart -= min(0.0, nmass_sap_inc) * nrelocfrac;				

				// C debt
				indiv.cmass_debt += cmass_debt_inc * indiv.densindiv;

				// N longtime reserves
				indiv.max_n_reserve_old = indiv.max_n_reserve;
				indiv.max_n_reserve = indiv.pft.n_reserve * indiv.cmass_sap / indiv.cton_leaf_new;

				// Allocation fractions
				if (!negligible(max(0.0, cmass_leaf_inc)) && !negligible(max(0.0, cmass_root_inc)) && !negligible(max(0.0, cmass_sap_inc))) {
					indiv.bminc_leaf_frac = max(0.0, cmass_leaf_inc) / (max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc) + max(0.0, cmass_sap_inc));
					indiv.bminc_root_frac = max(0.0, cmass_root_inc) / (max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc) + max(0.0, cmass_sap_inc));
				}

				if (indiv.alive) {
					patch.pft[indiv.pft.id].litter_leaf += litter_leaf_inc * indiv.densindiv;
					patch.pft[indiv.pft.id].litter_root += litter_root_inc * indiv.densindiv;

					patch.pft[indiv.pft.id].nmass_litter_leaf += litter_leaf_inc * indiv.densindiv /
						indiv.cton_leaf_old * (1.0 - nrelocfrac);
					indiv.nstore += litter_leaf_inc * indiv.densindiv / indiv.cton_leaf_old * nrelocfrac;

					patch.pft[indiv.pft.id].nmass_litter_root += litter_root_inc * indiv.densindiv /
						indiv.cton_root_old * (1.0 - nrelocfrac);
					indiv.nstore += litter_root_inc * indiv.densindiv / indiv.cton_root_old * nrelocfrac;
											
					// if sapwood gets killed transfer 50% of N into storage,
					// the other 50% going into heartwood
					if (cmass_sap_inc < 0.0)	
						indiv.nstore -= nmass_sap_inc*(1.0 - nrelocfrac);
				}
				else {	// return N to soil so N budget is preserved
					patch.soil.nmass_avail += litter_leaf_inc * indiv.densindiv / indiv.cton_leaf_old +
						litter_root_inc * indiv.densindiv / indiv.cton_root_old -
						min(0.0, nmass_sap_inc) * (1.0 - nrelocfrac); // nrelocfrac gone to heartwood above
				}

				// Update individual age

				indiv.age++;

				// Kill individual and transfer biomass to litter if any biomass
				// compartment negative

				if (indiv.cmass_leaf < MINCMASS || indiv.cmass_root < MINCMASS || 
					indiv.cmass_sap < MINCMASS) {

					// alive check
					if (indiv.alive) {

						// catches small, negative values too
						patch.pft[indiv.pft.id].litter_leaf += indiv.cmass_leaf;
						patch.pft[indiv.pft.id].litter_root += indiv.cmass_root;

						patch.pft[indiv.pft.id].litter_wood += indiv.cmass_sap;
						patch.pft[indiv.pft.id].litter_wood += indiv.cmass_heart - indiv.cmass_debt;
					
						patch.pft[indiv.pft.id].nmass_litter_leaf += max(indiv.nmass_leaf, 0.0);
						patch.pft[indiv.pft.id].nmass_litter_root += max(indiv.nmass_root, 0.0);

						patch.pft[indiv.pft.id].nmass_litter_wood += max(indiv.nmass_sap, 0.0) +
							max(indiv.nmass_heart, 0.0);
						
						// Transfer N storage to wood N litter for now
						patch.pft[indiv.pft.id].nmass_litter_wood += max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
					} 
					else {	// return N to soil so N budget is preserved
						patch.soil.nmass_avail += max(indiv.nmass_leaf, 0.0) + max(indiv.nmass_root, 0.0) + max(indiv.nmass_sap, 0.0) +
							max(indiv.nmass_heart, 0.0) + max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
					}

					vegetation.killobj();
					killed = true;
				}
			}
			else if (indiv.pft.lifeform == GRASS) {

				// GRASS GROWTH

				// initial grass cmass
				double indiv_mass_before = indiv.cmass_leaf + indiv.cmass_root;
	
				allocation(bminc, indiv.cmass_leaf, indiv.cmass_root,
					0.0, 0.0, 0.0, indiv.ltor, 0.0, 0.0, 0.0, GRASS, 0.0,
					0.0, 0.0, cmass_leaf_inc, cmass_root_inc, dval, dval, dval,
					litter_leaf_inc, litter_root_inc);

				// Calculate N needed for this new biomass
				indiv.ndemand =
						max(0.0, cmass_leaf_inc) * indiv.densindiv / indiv.cton_leaf_new +
						max(0.0, cmass_root_inc) * indiv.densindiv / indiv.cton_root_new;

				indiv.ndemand_no_nlim = indiv.ndemand;

				// Compute limitation factor based on balance between individual N demand and supply
				// (NB: this overwrites the alternative factor calculated in canexch.cpp, but this
				// one is better!)

				if (ifnlim && (indiv.nstore + EPS < indiv.ndemand || indiv.nstore < 0.0) && indiv.ndemand > 0.0) {
					if (indiv.nstore < 0.0)	
						indiv.limnfact = 0.0;
					else {
						// Use N from long-term storage to reduce N limitation
						double diff = indiv.ndemand - indiv.nstore;
						indiv.nstore += diff * min(indiv.nmass_reserve / diff, 2.0) / 2.0;
						indiv.nmass_reserve -= diff * min(indiv.nmass_reserve / diff, 2.0) / 2.0;

						indiv.limnfact = indiv.nstore / indiv.ndemand;
					}
				}
				else 
					indiv.limnfact = 1.0;

				// Nitrogen limitation of production
				if (ifnlim && date.year > freenyears && indiv.limnfact < 1.0) {

					double bminc_nlim, bminc_dec, gpp_dec, agpp;
					
					// Redo allocation, this time with N constraint	
					allocation_nlim(indiv.pft, indiv.nstore / indiv.densindiv, indiv.cton_leaf_new, indiv.cton_root_new,	
						0.0, bminc / indiv.densindiv, indiv.cmass_leaf / indiv.densindiv,
						indiv.cmass_root / indiv.densindiv, 0.0,
						0.0, indiv.ltor, 0.0,
						cmass_leaf_inc, cmass_root_inc, dval, dval,
						litter_leaf_inc, litter_root_inc, indiv.densindiv); 

					// Update N demand
					indiv.ndemand =
						max(0.0, cmass_leaf_inc) * indiv.densindiv / indiv.cton_leaf_new +
						max(0.0, cmass_root_inc) * indiv.densindiv / indiv.cton_root_new;
					
					// Calculate new biomass increment
					bminc_nlim = cmass_leaf_inc * indiv.densindiv +
							cmass_root_inc * indiv.densindiv; 

					// update accumulate annual C flux according to N limitation

					bminc_dec = bminc - bminc_nlim;
					agpp = 0.0;

					for (int mon=0; mon<12; mon++) 
						agpp += indiv.mgpp[mon];

					for (int month=0; month<12; month++) {
						gpp_dec = indiv.mgpp[month] * bminc_dec / agpp;
						indiv.mgpp[month] -= gpp_dec;
						indiv.mnpp[month] -= gpp_dec;
						if (indiv.alive) {
							patch.fluxes.mcflux_gpp[month] -= gpp_dec;
							patch.fluxes.acflux_veg += gpp_dec;
						}
					}	

					// Temporary: save biomass increment as new npp
					// NB: this is important because it affects growth efficiency and
					//     therefore mortality and litter fluxes (in vegdynam.cpp).

					bminc -= bminc_dec;
					indiv.anpp -= bminc_dec; 
				}

				// Subtract used nitrogen for biomass increment from storage pool
				if (ifnlim && date.year > freenyears)
					indiv.nstore -= indiv.ndemand;
				else 
					indiv.nstore = 0.0;

				// Optimal growth C:N
				if (!negligible(max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc)))
					indiv.cton_growth = (max(0.0, cmass_leaf_inc) * indiv.cton_leaf_opt +
						max(0.0, cmass_root_inc) * indiv.cton_leaf_opt * (indiv.cton_root_new / indiv.cton_leaf_new)) /
						(max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc));

				// Update carbon pools and litter (on area basis)
				// only litter in the case of 'alive' individuals

				// Leaves
				indiv.cmass_leaf += cmass_leaf_inc;
				indiv.nmass_leaf += cmass_leaf_inc / ((cmass_leaf_inc > 0) ? indiv.cton_leaf_new : indiv.cton_leaf_old);

				// Roots
				indiv.cmass_root += cmass_root_inc;
				indiv.nmass_root += cmass_root_inc / ((cmass_root_inc > 0) ? indiv.cton_root_new : indiv.cton_root_old);

				// N longtime reserves
				indiv.max_n_reserve_old = indiv.max_n_reserve;
				indiv.max_n_reserve = indiv.pft.n_reserve * indiv.cmass_root / indiv.cton_root_new;

				// Allocation fractions
				if (!negligible(max(0.0, cmass_leaf_inc)) && !negligible(max(0.0, cmass_root_inc))) {
					indiv.bminc_leaf_frac = max(0.0, max(0.0, cmass_leaf_inc) / (max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc)));
					indiv.bminc_root_frac = 1.0 - indiv.bminc_leaf_frac;
				}

				// Determine the (small) mass imbalance (kgC) for this individual. 
				// This can arise in the event of numerical errors in the allocation routine.
				double indiv_mass_after = indiv.cmass_leaf + indiv.cmass_root + litter_leaf_inc + litter_root_inc;
				double indiv_cmass_diff = (indiv_mass_before + bminc - indiv_mass_after);	

				// alive check before ensuring C balance
				if (indiv.alive) {

					patch.pft[indiv.pft.id].litter_leaf += litter_leaf_inc + indiv_cmass_diff / 2.0;
					patch.pft[indiv.pft.id].litter_root += litter_root_inc + indiv_cmass_diff / 2.0;

					patch.pft[indiv.pft.id].nmass_litter_leaf += litter_leaf_inc * indiv.densindiv /
						indiv.cton_leaf_old * (1.0 - nrelocfrac);
					indiv.nstore += litter_leaf_inc * indiv.densindiv / indiv.cton_leaf_old * nrelocfrac;

					patch.pft[indiv.pft.id].nmass_litter_root += litter_root_inc * indiv.densindiv /
						indiv.cton_root_old * (1.0 - nrelocfrac);
					indiv.nstore += litter_root_inc / indiv.cton_root_old * nrelocfrac;
				}
				else {	// return N to soil so N budget is preserved
					patch.soil.nmass_avail += litter_leaf_inc / indiv.cton_leaf_old + litter_root_inc / indiv.cton_root_old;
				}

				// Kill individual and transfer biomass to litter if either biomass
				// compartment negative

				if (indiv.cmass_leaf < MINCMASS || indiv.cmass_root < MINCMASS) {

					// alive check
					if (indiv.alive) {

						patch.pft[indiv.pft.id].litter_leaf += indiv.cmass_leaf;
						patch.pft[indiv.pft.id].litter_root += indiv.cmass_root;

						patch.pft[indiv.pft.id].nmass_litter_leaf += max(indiv.nmass_leaf, 0.0);
						patch.pft[indiv.pft.id].nmass_litter_root += max(indiv.nmass_root, 0.0); 
						
						// Transfer N storage to root N litter for now
						patch.pft[indiv.pft.id].nmass_litter_root += max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
					} 
					else {	// return N to soil so N budget is preserved
						patch.soil.nmass_avail += max(indiv.nmass_leaf, 0.0) + max(indiv.nmass_root, 0.0) +
							max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
					}

					vegetation.killobj();
					killed = true;
				}
			}
		}

		if (!killed) {

			if (!allometry(indiv)) {

				// alive check
				if (indiv.alive) {
					patch.pft[indiv.pft.id].litter_leaf += max(0.0, indiv.cmass_leaf);
					patch.pft[indiv.pft.id].litter_root += max(0.0, indiv.cmass_root);

					patch.pft[indiv.pft.id].litter_wood += max(0.0, indiv.cmass_sap) +
						indiv.cmass_heart - indiv.cmass_debt;

					patch.pft[indiv.pft.id].nmass_litter_leaf += max(indiv.nmass_leaf, 0.0);
					patch.pft[indiv.pft.id].nmass_litter_root += max(indiv.nmass_root, 0.0);

					patch.pft[indiv.pft.id].nmass_litter_wood += max(indiv.nmass_sap, 0.0) +
						max(indiv.nmass_heart, 0.0);
					
					// Transfer N storage to root N litter
					patch.pft[indiv.pft.id].nmass_litter_root += max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
				}
				else {	// return N to soil so N budget is preserved
					patch.soil.nmass_avail += max(indiv.nmass_leaf, 0.0) + max(indiv.nmass_root, 0.0) + max(indiv.nmass_sap, 0.0) +
						max(indiv.nmass_heart, 0.0) + max(indiv.nstore, 0.0) + max(indiv.nmass_reserve, 0.0);
				}

				vegetation.killobj();
				killed = true;
			}

			if (!killed) {
				if (!indiv.alive) {
					patch.fluxes.acflux_est -= indiv.cmass_leaf + indiv.cmass_root +
						indiv.cmass_sap + indiv.cmass_heart - indiv.cmass_debt;
					indiv.alive = true;
				}

				// N left, put in longterm storage
				if (indiv.nstore > 0.0 || indiv.nmass_reserve > indiv.max_n_reserve) {

					// individual longterm N storage
					if (date.year > freenyears) {

						// check if max nmass storage pool will be exceeded with the addition of nstore
						if (indiv.max_n_reserve < indiv.nmass_reserve) {
							indiv.nstore += indiv.nmass_reserve - indiv.max_n_reserve;
							indiv.nmass_reserve = indiv.max_n_reserve;
						}
						else if (indiv.nstore > indiv.max_n_reserve - indiv.nmass_reserve){
							indiv.nstore -= indiv.max_n_reserve - indiv.nmass_reserve;
							indiv.nmass_reserve = indiv.max_n_reserve;
						}
						else {
							indiv.nmass_reserve += indiv.nstore;
							indiv.nstore = 0.0;
						}
					}
					// Return the rest back to the soil
					patch.soil.nmass_avail += indiv.nstore;
					indiv.nstore = 0.0;	
				}
				else if (indiv.nstore < 0.0 && indiv.nmass_reserve > 0.0) {

					if (indiv.nmass_reserve > -indiv.nstore) {
						indiv.nmass_reserve += indiv.nstore;
						indiv.nstore = 0.0;
					}
					else {
						indiv.nstore += indiv.nmass_reserve;
						indiv.nmass_reserve = 0.0;
					}
				}

				// ... on to next individual
				vegetation.nextobj();
			}
		}
	}

	// Flush litter from reproduction straight to atmosphere
	flush_litter_repr(patch);
}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2001
// Huang, S, Titus, SJ & Wiens, DP (1992) Comparison of nonlinear height-diameter
//   functions for major Alberta tree species. Canadian Journal of Forest Research 22:
//   1297-1304
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
// Zaehle, S. & Friend, A. D. 2010. Carbon and nitrogen cycle dynamics in the O-CN land surface 
//   model: 1. Model description, site-scale evaluation, and sensitivity to parameter estimates. 
//   Global Biogeochemical Cycles, 24.
// Zeide, B (1993) Primary unit of the tree crown. Ecology 74: 1598-1602.
