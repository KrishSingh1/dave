///////////////////////////////////////////////////////////////////////////////////////
/// \file canexch.cpp
/// \brief The canopy exchange module
///
/// Vegetation-atmosphere exchange of H2O and CO2 via
/// production, respiration and evapotranspiration.
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
#include "canexch.h"
#include "driver.h"
#include "q10.h"
#include "bvoc.h"
#include "ncompete.h"
#include "somdynam.h"
#include <assert.h>


// Anonymous namespace for variables with file scope
namespace {

/// leaf nitrogen (kgN/kgC) not associated with photosynthesis
/** (value given by Haxeltine & Prentice 1996b) */
const double N0 = 7.15 * G_PER_MG;

// Lookup tables for parameters with Q10 temperature responses

/// lookup table for Q10 temperature response of CO2/O2 specificity ratio
LookupQ10 lookup_tau(0.57, 2600.0);

/// lookup table for Q10 temperature response of Michaelis constant for O2
LookupQ10 lookup_ko(1.2, 3.0e4);

/// lookup table for Q10 temperature response of Michaelis constant for CO2
LookupQ10 lookup_kc(2.1, 30.0);

}

///////////////////////////////////////////////////////////////////////////////////////
// INTERCEPTION

/// Daily loss of water and energy through evaporation of rain or snow intercepted by the vegetation canopy
/** Gerten et al. (2004) Eq 2-4.

  */
void interception(Patch& patch,Climate& climate) {

	// Calculates daily loss of water and energy through evaporation of rainfall
	// intercepted by the vegetation canopy

	double scap; // canopy storage capacity (mm)
	double fwet; // fraction of day that canopy is wet (Kergoat 1996)
	double pet; // potential evapotranspiration (mm)

	pet=climate.eet*PRIESTLEY_TAYLOR;

	// Retrieve Vegetation object
	Vegetation& vegetation=patch.vegetation;

	patch.intercep=0.0;

	// Loop through individuals ...

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		if (!negligible(pet)) {

			if (indiv.alive) {
				// Storage capacity for precipitation by canopy (point scale)
				scap=climate.prec*min(indiv.lai_indiv_today()*indiv.pft.intc,0.999);

				// Fraction of day that canopy remains wet
				fwet=min(scap/pet,patch.fpc_rescale);

				// Calculate interception by this individual, and increment patch total

				indiv.intercep=fwet*pet*indiv.fpc;
				patch.intercep+=indiv.intercep;
			}
			else {
				indiv.intercep=0.0;
			}
		}
		else {

			indiv.intercep=0.0;
			patch.intercep=0.0;
		}

		// ... on to next individual
		vegetation.nextobj();
	}

	// Calculate net EET for vegetated parts of patch (deducting loss to interception)

	patch.eet_net_veg=max(climate.eet-patch.intercep,0.0);

	// Interception accounting for patch
	patch.aintercep+=patch.intercep;
	patch.mintercep[date.month]+=patch.intercep;

}

double get_veg_height(Vegetation &vegetation) {
	double height_veg = 0;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		if (indiv.growingseason()) {
			if (indiv.height>height_veg) {
				height_veg=indiv.height;
			}
		}

		vegetation.nextobj(); // ... on to next individual
	}
	return height_veg;
}

///////////////////////////////////////////////////////////////////////////////////////
// FPAR
// Internal function - not intended to be called by framework

void fpar(Patch& patch) {

	// DESCRIPTION
	// Calculates daily fraction of incoming PAR (FPAR) taken up by individuals in a
	// particular patch over their projective areas, given current leaf phenological
	// status. Calculates PAR and FPAR at top of grass canopy (individual and cohort
	// modes). Calculates fpar assuming leaf-on (phen=1) for all vegetation.
	//
	// Note: In order to compensate for the additional computational cost of
	//       calculating fpar_leafon in cohort/individual mode, the grain of the
	//       integration of FPAR through the canopy has been increased from 1 to 2 m
	//
	// NEW ASSUMPTIONS CONCERNING FPC AND FPAR (Ben Smith 2002-02-20)
	// FPAR = average individual fraction of PAR absorbed on patch basis today,
	//        including effect of current leaf phenology (this differs from previous
	//        versions of LPJ-GUESS in which FPAR was on an FPC basis)
	// FPC =  PFT population (population mode), cohort (cohort mode) or individual
	//        (individual mode) fractional projective cover as a fraction of patch area
	//        (in population mode, corresponds to LPJF variable fpc_grid). Updated
	//        annually based on leaf-out LAI (see function allometry in growth module).
	//        (FPC was previously equal to summed crown area as a fraction of patch
	//        area in cohort/individual mode)
	//
	// Population mode: FPAR on patch (grid cell) area basis assumed to be equal to fpc
	// under full leaf cover; i.e.
	//     (1) fpar        = fpc*phen
	//     (2) fpar_leafon = fpc
	//
	// Individual and cohort modes: FPAR calculated assuming trees shade themselves
	//   and all individuals below them according to the Lambert-Beer law (Prentice
	//   et al 1993, Eqn 27; Monsi & Saeki 1953):
	//     (3) fpar = integral [0-tree height] exp ( -k * plai(z) )
	//   where
	//       k       = extinction coefficient;
	//       plai(z) = summed leaf-area index for leaves of all individuals, above
	//                 canopy depth z, taking account of current phenological status
 
	const double PHEN_GROWINGSEASON=0.5;
		// minimum expected vegetation leaf-on fraction for growing season

	double plai; // cumulative leaf-area index (LAI) for patch (m2 leaf/m2 ground)
	double plai_leafon;
		// cumulative LAI for patch assuming full leaf cover for all individuals
	double plai_layer; // summed LAI by layer for patch
	double plai_leafon_layer;
		// summed LAI by layer for patch assuming full leaf cover for all individuals
	double plai_grass; // summed LAI for grasses
	double plai_leafon_grass;
		// summed LAI for grasses assuming full leaf cover for all individuals
	double flai; // fraction of total grass LAI represented by a particular grass
	double fpar_layer_top; // FPAR by layer
	double fpar_leafon_layer_top;
		// FPAR by layer assuming full leaf cover for all individuals
	double fpar_layer_bottom;
	double fpar_leafon_layer_bottom;
	double fpar_grass; // FPAR at top of grass canopy
	double fpar_leafon_grass;
		// FPAR at top of grass canopy assuming full leaf cover for all individuals
	double fpar_ff; // FPAR at forest floor (beneath grass canopy)
	double fpar_leafon_ff;
		// FPAR at forest floor assuming full leaf cover for all individuals
	double frac;
		// vertical fraction of layer occupied by crown cylinder(s) of a particular
		// individual or cohort
	double atoh; // term in calculating LAI sum for a given layer
	double height_veg; // maximum vegetation height (m)
	int toplayer; // number of vertical layers of width vstep in vegetation (minus 1)
	int layer; // layer number (0=lowest)
	double lowbound; // lower bound of current layer (m)
	double highbound; // upper bound of current layer (m)
	double fpar_min; // minimum FPAR required for grass growth
	double par_grass; // PAR reaching top of grass canopy (J/m2/day)
	double phen_veg; // LAI-weighted mean fractional leaf-out for vegetation
	//variables needed for "S�kes" FPAR scheme
	double fpar_uptake_layer;
	double fpar_uptake_leafon_layer;

	// Obtain reference to Vegetation object
	Vegetation& vegetation=patch.vegetation;

	// And to Climate object
	const Climate& climate = patch.get_climate();

	if (vegmode==POPULATION) {

		// POPULATION MODE

		// Loop through individuals

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			// For this individual ...

			indiv.fpar=indiv.fpc_today(); // Eqn 1
			indiv.fpar_leafon=indiv.fpc * indiv.growingseason(); // Eqn 2

			vegetation.nextobj(); // ... on to next individual
		}
	}

	else {

		// INDIVIDUAL OR COHORT MODE

		// Initialise individual FPAR, find maximum height of vegetation, calculate
		// individual LAI given current phenology, calculate summed LAI for grasses

		plai=0.0;
		plai_leafon=0.0;
		plai_grass=0.0;
		plai_leafon_grass=0.0;
		phen_veg=0.0;
		height_veg=0.0;

		// Loop through individuals

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			// For this individual ...

			if (indiv.growingseason()) {

				indiv.fpar=0.0;
				indiv.fpar_leafon=0.0;
				// Empty the daily layerwise-fpar vector.
				indiv.fpar_layerwise.clear();
				indiv.lai_layerwise.clear();
				if (indiv.height>height_veg) height_veg=indiv.height;
				plai_leafon+=indiv.lai;

				if (indiv.pft.lifeform==GRASS || indiv.pft.lifeform == MOSS) {
					plai_leafon_grass+=indiv.lai;
					plai_grass += indiv.lai_today();
				}

				// Accumulate LAI-weighted sum of individual leaf-out fractions
				phen_veg+=indiv.lai_today();
			}

			vegetation.nextobj(); // ... on to next individual
		}

		// Calculate LAI-weighted mean leaf-out fraction for vegetation
		// guess2008 - bugfix - was: if (!negligible(plai))
		if (!negligible(plai_leafon))
			phen_veg/=plai_leafon;
		else
			phen_veg=1.0;

		// width of vertical layers for canopy-area integration (m)
		double vstep = nlayer_canopy ? height_veg / nlayer_canopy : layer_thickness_canopy;

		// Calculate number of layers (minus 1) from ground surface to top of canopy
		toplayer = (int)(height_veg/vstep-0.0001);

		// Calculate FPAR by integration from the top of the canopy (Eqn 2)
		plai=0.0;
		plai_leafon=0.0;

		// Set FPAR for bottom of layer above (initially 1 at top of canopy)

		fpar_layer_bottom=1.0;
		fpar_leafon_layer_bottom=1.0;

		for (layer=toplayer;layer>=0;layer--) {

			lowbound=(double)layer*vstep;
			highbound=lowbound+vstep;

			// FPAR at top of this layer = FPAR at bottom of layer above

			fpar_layer_top=fpar_layer_bottom;
			fpar_leafon_layer_top=fpar_leafon_layer_bottom;

			plai_layer=0.0;
			plai_leafon_layer=0.0;

			// Loop through individuals

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();

				// For this individual ...

				if (indiv.pft.lifeform==TREE) {
					if (indiv.height>lowbound && indiv.boleht<highbound &&
						!negligible(indiv.height-indiv.boleht)) {

						// Calculate vertical fraction of current layer occupied by
						// crown cylinders of this cohort

						frac=1.0;
						if (indiv.height<highbound)
							frac-=(highbound-indiv.height)/vstep;
						if (indiv.boleht>lowbound)
							frac-=(indiv.boleht-lowbound)/vstep;

						// Calculate summed LAI of this cohort in this layer

						atoh=indiv.lai/(indiv.height-indiv.boleht);
						indiv.lai_leafon_layer=atoh*frac*vstep;
						plai_layer+=indiv.lai_leafon_layer*indiv.phen;
						plai_leafon_layer+=indiv.lai_leafon_layer;
					}
					else {
						indiv.lai_layer=0.0;
						indiv.lai_leafon_layer=0.0;
					}
				}

				// ... on to next individual
				vegetation.nextobj();
			}

			// Update cumulative LAI for this layer and above
			plai+=plai_layer;
			plai_leafon+=plai_leafon_layer;

			// Calculate FPAR at bottom of this layer
			// Eqn 27, Prentice et al 1993

			fpar_layer_bottom = lambertbeer(plai);
			fpar_leafon_layer_bottom = lambertbeer(plai_leafon);

			// Total PAR uptake in this layer

			fpar_uptake_layer=fpar_layer_top-fpar_layer_bottom;
			fpar_uptake_leafon_layer=fpar_leafon_layer_top-fpar_leafon_layer_bottom;

			// Partition PAR for this layer among trees,

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();

				// For this individual ...

				if (indiv.pft.lifeform==TREE) {
					if (!negligible(plai_leafon_layer))

						// FPAR partitioned according to the relative amount
						// of leaf area in this layer for this individual

						indiv.fpar_leafon+=fpar_uptake_leafon_layer*
							indiv.lai_leafon_layer/plai_leafon_layer;

					if (!negligible(plai_layer)) {
						double fpar_layer = fpar_uptake_layer*
							(indiv.lai_leafon_layer*indiv.phen)/plai_layer;
						indiv.fpar+=fpar_layer;
						// Store layer-wise fpar in the individual.
						indiv.fpar_layerwise.push_back(fpar_layer);
						indiv.lai_layerwise.push_back(indiv.lai_leafon_layer);
					}
				}

				// ... on to next individual
				vegetation.nextobj();
			}

		}

		// FPAR reaching grass canopy
		fpar_grass = lambertbeer(plai);
		fpar_leafon_grass = lambertbeer(plai_leafon);

		// Add grass LAI to calculate PAR reaching forest floor
		// BLARP: Order changed Ben 050301 to overcome optimisation bug in pgCC

		//plai+=plai_grass;
		fpar_ff = lambertbeer(plai+plai_grass);
		plai+=plai_grass;

		// Save this
		patch.fpar_ff=fpar_ff;

		plai_leafon+=plai_leafon_grass;
		fpar_leafon_ff = lambertbeer(plai_leafon);

		// FPAR for grass PFTs is difference between relative PAR at top of grass canopy
		// canopy and at forest floor, or lower if FPAR at forest floor below threshold
		// for grass growth. PAR reaching the grass canopy is partitioned among grasses
		// in proportion to their LAI (a somewhat simplified assumption)

		// Loop through individuals

		double fpar_tree_total=0.0;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			// For this individual ...

			if (indiv.pft.lifeform==GRASS || indiv.pft.lifeform==MOSS) {

				// Calculate minimum FPAR for growth of this grass

				// Fraction of total grass LAI represented by this grass

				if (!negligible(plai_grass))
					flai=indiv.lai_today()/plai_grass;
				else
					flai=1.0;

				// If using dave photosynthesis, just ignore parff_min and
				// fpar_min. We check in photosynthesis_gs() whether there's
				// enough PAR for photosynthesis to occur.
				double t2 = fpar_ff * flai;
				if (!ifdavephotosyn) {
					if (!negligible(climate.par))
						fpar_min=min(indiv.pft.parff_min/climate.par,1.0);
					else
						fpar_min=1.0;

					t2 = max(t2, fpar_min);
				}
				indiv.fpar=max(0.0,fpar_grass*flai-t2);
				indiv.fpar_layerwise.push_back(indiv.fpar);
				indiv.lai_layerwise.push_back(indiv.lai);
				// Repeat assuming full leaf cover for all individuals

				if (!negligible(plai_leafon_grass))
					flai=indiv.lai/plai_leafon_grass;
				else
					flai=1.0;

				indiv.fpar_leafon=max(0.0,fpar_leafon_grass*flai-
					max(fpar_leafon_ff*flai,fpar_min));
			}

			if (indiv.pft.lifeform==TREE) fpar_tree_total+=indiv.fpar;

			vegetation.nextobj();
		}

		// Save grass canopy FPAR and update mean growing season grass canopy PAR
		// Growing season defined here as days when mean vegetation leaf-on fraction
		// exceeds 50% and we're in the light half of the year (daylength >= 11).
		//
		// The daylength condition was added because sites with evergreens can have
		// a mean vegetation leaf-on fraction over 50% even during polar night.
		// 11 hours was chosen because some sites never reach exactly 12 hours, the
		// exact limit shouldn't matter much.

		patch.fpar_grass=fpar_grass;
		par_grass=fpar_grass*climate.par;

		if (date.day==0) {
			patch.par_grass_mean=0.0;
			patch.nday_growingseason=0;
		}

		if (phen_veg>PHEN_GROWINGSEASON && patch.get_climate().daylength >= 11.0) {
			patch.par_grass_mean+=par_grass;
			patch.nday_growingseason++;
		}

		// Convert from sum to mean on last day of year

		if (date.islastday && date.islastmonth && patch.nday_growingseason) {
			patch.par_grass_mean/=(double)patch.nday_growingseason;
		}
	}
}

double alphaa(const Pft& pft) {

	if (pft.phenology == CROPGREEN)
		return ifnlim ? ALPHAA_CROP_NLIM : ALPHAA_CROP;
	else
		return ifnlim ? ALPHAA_NLIM : ALPHAA;
}

/// Non-water stressed rubisco capacity, with or without nitrogen limitation
void vmax(double b, double c1, double c2, double apar, double tscal,
		  double daylength, double temp, double nactive, bool ifnlimvmax, double& vm, double& vmaxnlim, double& nactive_opt) {

	// Calculation of non-water-stressed rubisco capacity assuming leaf nitrogen not
	// limiting (Eqn 11, Haxeltine & Prentice 1996a)
	// Calculation of sigma is based on Eqn 12 (same source)

	double s =  24.0 / daylength * b;
	double sigma = sqrt(max(0., 1. - (c2 - s) / (c2 - THETA * s)));
	vm = 1 / b * CMASS * CQ * c1 / c2 * tscal * apar *
						(2. * THETA * s * (1. - sigma) - s + c2 * sigma);

	// Calculate nitrogen-limited Vmax for current leaf nitrogen
	// Haxeltine & Prentice 1996b Eqn 28

	const double M = 25.0; // corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996b

	// Conversion factor in calculation of leaf nitrogen: includes conversion of:
	//		- Vm from gC/m2/day to umolC/m2/sec
	//      - nitrogen from mg/m2 to kg/m2

	double CN = 1.0 / (3600 * daylength * CMASS);

	double tfac = exp(-0.0693 * (temp - 25.0));
	double vm_max = nactive / (M * CN * tfac);

	// Calculate optimal leaf nitrogen based on [potential] Vmax (Eqn 28 Haxeltine & Prentice 1996b)
	nactive_opt = M * vm * CN * tfac;

	if (vm > vm_max && ifnlimvmax) {
		vmaxnlim = vm_max / vm;	// Save vmax nitrogen limitation
		vm = vm_max;
	}
	else {
		vmaxnlim = 1.0;
	}
}

/// Total daily gross photosynthesis
/** Calculation of total daily gross photosynthesis and leaf-level net daytime
 *  photosynthesis given degree of stomatal closure (as parameter lambda).
 *  Includes implicit scaling from leaf to plant projective area basis.
 *  Adapted from Farquhar & von Caemmerer (1982) photosynthesis model, as simplified
 *  by Collatz et al (1991), Collatz et al (1992), Haxeltine & Prentice (1996a,b)
 *  and Sitch et al. (2000).
 *
 *  To calculate vmax call w/ daily averages of temperature and par.
 *  Vmax is to be calculated daily and only with lambda == lambda_max.
 *  lambda values greater than lambda_max are forbidden.
 *  In sub-daily mode daylength should be 24 h, to obtain values in daily units.
 *
 *  INPUT PARAMETERS
 * 
 *  \param PhotosynthesisEnvironment struct containing the following public members:
 *   - co2        atmospheric ambient CO2 concentration (ppmv)
 *   - temp       mean air temperature today (deg C)
 *   - par        total daily photosynthetically-active radiation today (J/m2/day)
 *   - daylength  day length, must equal 24 in diurnal mode (h)
 *   - fpar       fraction of PAR absorbed by foliage
 *
 *  \param PhotosynthesisStresses struct containing the following members:
 *   - ifnlimvmax			- whether nitrogen should limit Vmax
 *   - moss_ps_limit		- limit to moss photosynthesis. [0,1], where 1 means no limit
 *   - graminoid_ps_limit	- limit to graminoid photosynthesis. [0,1], where 1 means no limit
 *   - inund_stress			- limit to photosynthesis due to inundation, where 1 means no limit
 *
 *  \param pft        Pft object containing the following public members:
 *   - pathway         biochemical pathway for photosynthesis (C3 or C4)
 *   - pstemp_min      approximate low temperature limit for photosynthesis (deg C)
 *   - pstemp_low      approximate lower range of temperature optimum for
 *                     photosynthesis (deg C)
 *   - pstemp_high     approximate upper range of temperature optimum for photosynthesis
 *                     (deg C)
 *   - pstemp_max      maximum temperature limit for photosynthesis (deg C)
 *   - lambda_max      non-water-stressed ratio of intercellular to ambient CO2 pp
 *
 *  \param lambda     ratio of intercellular to ambient partial pressure of CO2
 *
 *  \param nactive    nitrogen available for photosynthesis
 *
 *  \param vm         pre-calculated value of Vmax for this stand for this day if
 *                    available, otherwise calculated
 *
 * OUTPUT PARAMETERS
 *
 * \param ps_result      see documentation of PhotosynthesisResult struct
 * 
 *
 * IMPORTANT for users adding new call parameters to the list above:
 *
 * Never place new call parameters in the proper photosynthesis() function header, instead
 * place new parameters insided the structs PhotosynthesisEnvironment and PhotosynthesisStresses,
 * or in PhotosynthesisResult if it is a result.
 */
void photosynthesis(const PhotosynthesisEnvironment& ps_env, 
					const PhotosynthesisStresses& ps_stresses,
					const Pft& pft,
					double lambda, 
					double nactive, 
					double vm, 
					PhotosynthesisResult& ps_result) {

	// NOTE: This function is identical to LPJF subroutine "photosynthesis" except for
	// the formulation of low-temperature inhibition coefficient tscal (tstress; LPJF).
	// The function adopted here draws down metabolic activity in approximately the
	// temperature range pstemp_min-pstemp_low but does not affect photosynthesis
	// at high temperatures.

	// HISTORY
	// Ben Smith 18/1/2001: Tested in comparison to LPJF subroutine "photosynthesis":
	// function showed identical behaviour except at temperatures >= c. 35 deg C where
	// LPJF temperature inhibition function results in lower photosynthesis.

	// Make sure that only two alternative modes are possible:
	//  * daily non-water stressed (forces Vmax calculation)
	//  * with pre-calculated Vmax (sub-daily and water-stressed)
	assert(vm >= 0 || lambda == pft.lambda_max);
	assert(lambda <= pft.lambda_max);

	// Get the environmental variables
	double temp = ps_env.get_temp();
	double co2 = ps_env.get_co2();
	double fpar = ps_env.get_fpar();
	double par = ps_env.get_par();
	double daylength = ps_env.get_daylength();

	// Get the stresses
	bool ifnlimvmax = ps_stresses.get_ifnlimvmax();

	// No photosynthesis during polar night, outside of temperature range or no RuBisCO activity
	if (negligible(daylength) || negligible(fpar) || temp > pft.pstemp_max || temp < pft.pstemp_min || !vm) {
		ps_result.clear();
		return;
	}


	// Scale fractional PAR absorption at plant projective area level (FPAR) to
	// fractional absorption at leaf level (APAR)
	// Eqn 4, Haxeltine & Prentice 1996a
	double apar = par * fpar * alphaa(pft);
	double b, c1, c2;


	// Calculate temperature-inhibition coefficient
	// This function (tscal) is mathematically identical to function tstress in LPJF.
	// In contrast to earlier versions of modular LPJ and LPJ-GUESS, it includes both
	// high- and low-temperature inhibition.
	double k1 = (pft.pstemp_min+pft.pstemp_low) / 2.0;
	double tscal = (1. - .01*exp(4.6/(pft.pstemp_max-pft.pstemp_high)*(temp-pft.pstemp_high)))/
							(1.0+exp((k1-temp)/(k1-pft.pstemp_min)*4.6));

	if (pft.pathway == C3) 	{			// C3 photosynthesis

		// Calculate CO2 compensation point (partial pressure)
		// Eqn 8, Haxeltine & Prentice 1996a
		double gammastar = PO2 / 2.0 / lookup_tau[temp];

		// Intercellular partial pressure of CO2 given stomatal opening (Pa)
		// Eqn 7, Haxeltine & Prentice 1996a
		double pi_co2 = lambda * co2 * PATMOS * CO2_CONV;

		// Calculation of C1_C3, Eqn 4, Haxeltine & Prentice 1996a
		// High-temperature inhibition modelled by suppression of LUE by decreased
		// relative affinity of rubisco for CO2 with increasing temperature (Table 3.7,
		// Larcher 1983)
		// Notes: - there is an error in Eqn 4, Haxeltine & Prentice 1996a (missing
		//          2.0* in denominator) which is fixed here (see Eqn A2, Collatz
		//          et al 1991)
		//        - the explicit low temperature inhibition function has been removed
		//          and replaced by a temperature-dependent upper limit on V_m, see
		//          below
		//        - the reduction in maximum photosynthesis due to leaf age (phi_c)
		//          has been removed
		//        - alpha_a, accounting for reduction in PAR utilisation efficiency
		//          from the leaf to ecosystem level, appears in the calculation of
		//          apar (above) instead of here
		//        - C_mass, the atomic weight of carbon, appears in the calculation
		//          of V_m instead of here
		c1 = (pi_co2 - gammastar) / (pi_co2 + 2.0 * gammastar) * ALPHA_C3;

		// Calculation of C2_C3, Eqn 6, Haxeltine & Prentice 1996a
		c2 = (pi_co2 - gammastar) / (pi_co2 + lookup_kc[temp] * (1.0 + PO2/lookup_ko[temp]));

		// see Wania et al. 2009b
		if (pft.lifeform == MOSS) {
			b = BC_moss;
		}
		else {
			b = BC3;
		}
	}
	else {							// C4 photosynthesis
		// Calculation of C1_C4 given actual pi (lambda)
		// C1_C4 incorporates term accounting for effect of intercellular CO2
		// concentration on photosynthesis (Eqn 14, 16, Haxeltine & Prentice 1996a)
		c1 = min(lambda/LAMBDA_SC4, 1.0) * ALPHA_C4;
		c2 = 1;
		b = BC4;
	}
	if (vm < 0) {

		// Calculation of non-water-stressed rubisco capacity (Eqn 11, Haxeltine & Prentice 1996a)
		vmax(b, c1, c2, apar, tscal, daylength, temp, nactive, ifnlimvmax, ps_result.vm, ps_result.vmaxnlim, ps_result.nactive_opt);
	}
	else {
		ps_result.vm = vm;			// reuse existing Vmax
	}
	// Calculation of daily leaf respiration
	// Eqn 10, Haxeltine & Prentice 1996a
	ps_result.rd_g = ps_result.vm * b;

	// PAR-limited photosynthesis rate (gC/m2/h)
	// Eqn 3, Haxeltine & Prentice 1996a
	ps_result.je = c1 * tscal * apar * CMASS * CQ / daylength;

	// Rubisco-activity limited photosynthesis rate (gC/m2/h)
	// Eqn 5, Haxeltine & Prentice 1996a
	double jc = c2 * ps_result.vm / 24.0;

	// Calculation of daily gross photosynthesis
	// Eqn 2, Haxeltine & Prentice 1996a
	// Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing
	// 			theta in 4*theta*je*jc term) which is fixed here
	ps_result.agd_g = (ps_result.je + jc - sqrt((ps_result.je + jc) * (ps_result.je + jc) - 4.0 * THETA * ps_result.je * jc)) /
														(2.0 * THETA) * daylength;


	if (!iftwolayersoil) {

		// LIMITS TO PHOTOSYNTHESIS
		// On wetlands, both agd_g and rd_g are scaled in the event of inundation (all PFTS),
		// or, for mosses and graminoids only, in the event of the water table dropping below
		// the PFT's optimal water table depth (dessication). See Wania et al. (2009b)

		// Get the stresses
		double moss_ps_limit = ps_stresses.get_moss_ps_limit();
		double graminoid_ps_limit = ps_stresses.get_graminoid_ps_limit();
		double inund_stress = ps_stresses.get_inund_stress();

		// 1) Inundation stress
		// Reduce GPP if there is inundation stress
		// (possibility of) inund stress (i.e. values < 1) only on PEATLAND stands and when ifinundationstress is true
		ps_result.agd_g *= inund_stress; 
		ps_result.rd_g *= inund_stress;

		// 2a) Moss dessication
		if (pft.lifeform == MOSS) {
			// Reduce agd_g using moss_wtp_limit (lies between [0.3, 1.0])
			ps_result.agd_g *= moss_ps_limit;
			ps_result.rd_g *= moss_ps_limit;
		}
	
		// 2b) Graminoid dessication (NB! all other PFTs have graminoid_ps_limit == 1)
		ps_result.agd_g *= graminoid_ps_limit;
		ps_result.rd_g *= graminoid_ps_limit;
	
	}

    // Leaf-level net daytime photosynthesis (gC/m2/day)
	// Based on Eqn 19, Haxeltine & Prentice 1996a
	double adt = ps_result.agd_g - daylength / 24.0 * ps_result.rd_g;

	// Convert to CO2 diffusion units (mm/m2/day) using ideal gas law
	ps_result.adtmm = adt / CMASS * 8.314 * (temp + K2degC) / PATMOS * 1e3;
}

/// Calculate value for canopy conductance component associated with photosynthesis (mm/s)
/** Eqn 21, Haxeltine & Prentice 1996
 *  includes conversion of daylight from hours to seconds
 */
inline double gpterm(double adtmm, double co2, double lambda, double daylength) {
	if (adtmm <= 0) {
		return 0;
	}
	return 1.6 / CO2_CONV / 3600 * adtmm / co2 / (1 - lambda) / daylength;
}


/// Determine CO2 in peatland water
/**
 * Updated daily
 */
double get_co2(Patch& p, Climate& climate, Pft& pft) {

	double pftco2 = climate.co2;

	if (p.stand.is_highlatitude_peatland_stand() && pft.ismoss())
		pftco2 = p.soil.acro_co2; // override for peat mosses

	return pftco2;
}

/// Determine inundation stress
/**
 * Updated daily
 */
double get_inund_stress(Patch& p, Patchpft& ppft) {

	double inund_stress = 1.0; // No stress by default

	if (ifinundationstress && p.stand.is_highlatitude_peatland_stand())
		inund_stress = ppft.inund_stress; // Determine inundation stress for this Patchpft

	return inund_stress;
}

/// Determine limits on graminoid photosynthesis due to WTP
/**
 * Updated daily
 */
double get_graminoid_wtp_limit(Patch& p, Pft& pft) {

	double graminoid_wtp_limit = 1.0; // No limit by default

	// Update limit if this is a graminoid
	if (!pft.ismoss() && pft.iswetlandspecies() && p.stand.is_highlatitude_peatland_stand())
		graminoid_wtp_limit = p.soil.dgraminoid_wtp_limit;

	return graminoid_wtp_limit;
}

/// Determine limits on moss photosynthesis due to WTP
/**
 * Updated daily
 */
double get_moss_wtp_limit(Patch& p, Pft& pft) {

	double moss_wtp_limit = 1.0; // No limit by default

	// Update limit if this is a moss
	if (pft.ismoss() && p.stand.is_highlatitude_peatland_stand())
		moss_wtp_limit = p.soil.dmoss_wtp_limit;

	return moss_wtp_limit;
}


/// Pre-calculate Vmax and no-stress assimilation and canopy conductance
/**
 * Vmax is calculated on a daily scale (w/ daily averages of temperature and par)
 * Subdaily values calculated if needed
 */
void photosynthesis_nostress(Patch& patch, Climate& climate) {

	PhotosynthesisEnvironment ps_env;					
	PhotosynthesisStresses ps_stress;
	ps_stress.no_stress();

	// If this is the first patch, calculate no-stress assimilation for
	// each Standpft, assuming FPAR=1. This is then later used in
	// forest_floor_conditions.

	double pftco2 = climate.co2; // will override for peat mosses

	if (!patch.id) {

		for (int p=0; p<npft; p++) {
			Standpft& spft = patch.stand.pft[p];
			Patchpft& ppft = patch.pft[p];

			if (spft.active) {

				pftco2 = get_co2(patch, climate, spft.pft);
				ps_env.set(pftco2, climate.temp, climate.par, 1.0, climate.daylength);

				ps_stress.set(false, get_moss_wtp_limit(patch, spft.pft), get_graminoid_wtp_limit(patch, spft.pft), get_inund_stress(patch, ppft));

				// Call photosynthesis assuming stomates fully open (lambda = lambda_max)
				photosynthesis(ps_env, ps_stress, spft.pft, 
							   spft.pft.lambda_max, 1.0, -1, 
							   spft.photosynthesis);
			}
		}
	}

	// Pre-calculation of no-stress assimilation for each individual
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();

	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[indiv.pft.id];

		pftco2 = get_co2(patch, climate, pft);
		ps_env.set(pftco2, climate.temp, climate.par, indiv.fpar, climate.daylength);

		ps_stress.set(false, get_moss_wtp_limit(patch, pft), get_graminoid_wtp_limit(patch, pft), get_inund_stress(patch, ppft));

		// Individual photosynthesis with no nitrogen limitation
		photosynthesis(ps_env, ps_stress, pft,
		               pft.lambda_max, 1.0, -1,
		               indiv.photosynthesis);

		indiv.gpterm = gpterm(indiv.photosynthesis.adtmm, pftco2, pft.lambda_max, climate.daylength);

		if (date.diurnal()) {

			indiv.gpterms.assign(date.subdaily, 0);
			PhotosynthesisResult res;
			indiv.phots.assign(date.subdaily, res);

			for (int i=0; i<date.subdaily; i++) {

				PhotosynthesisResult& ps_result = indiv.phots[i];
				// Update temperature and PAR
				ps_env.set(pftco2, climate.temps[i], climate.pars[i], indiv.fpar, 24);

				photosynthesis(ps_env, ps_stress, pft,
				               pft.lambda_max, 1.0, indiv.photosynthesis.vm,
				               ps_result);

				indiv.gpterms[i] = gpterm(ps_result.adtmm, climate.co2, pft.lambda_max, 24);
			}
		}
		vegetation.nextobj();
	}
}

/// Calculates individual fnuptake based on surface of fine root
/** Calculates individual fraction nitrogen uptake based on surface of fine root
 *  Roots are cone formed with height == radius.
 *  V = PI * r^3 / 3
 *  A = (2^1/2 + 1) * PI * r^2
 *  -> A = const * cmass_root^2/3
 */
double nitrogen_uptake_strength(const Individual& indiv) {
	return pow(max(0.0, indiv.cmass_root_today()) * indiv.pft.nupscoeff * indiv.cton_status / indiv.densindiv, 2.0 / 3.0) * indiv.densindiv;
}

/// Individual nitrogen uptake fraction
/** Determining individual nitrogen uptake as a fraction of its nitrogen demand.
 *
 *  \see ncompete
 *
 *  Function nitrogen_uptake_strength() determines how good individuals are at
 *  acquiring nitrogen.
 */
void fnuptake(Vegetation& vegetation, double nmass_avail) {

	// Create vector describing the individuals to ncompete()
	std::vector<NCompetingIndividual> individuals(vegetation.nobj);

	for (unsigned int i = 0; i < vegetation.nobj; i++) {
		individuals[i].ndemand = vegetation[i].ndemand;
		individuals[i].strength = nitrogen_uptake_strength(vegetation[i]);
	}

	// Let ncompete() do the actual distribution
	ncompete(individuals, nmass_avail);

	// Get the results, nitrogen uptake fraction for each individual
	for (unsigned int i = 0; i < vegetation.nobj; i++) {
		vegetation[i].fnuptake = individuals[i].fnuptake;
	}
}


/// Use nitrogen storage to limit stress
/** Retranslocated nitrogen from last year is used to
 *  limit nitrogen stress in leaves, roots, and sap wood
 */
void nstore_usage(Vegetation& vegetation) {

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		if (indiv.pft.lifeform == GRASS && ifdailygrass) {

			indiv.nstress = ifnlim && date.year > freenyears;

			// Check if current tissue N mass matches optimal level otherwise N limitation on photosynthesis
			if (indiv.get_nmass_leaf() < indiv.get_cmass_leaf() / indiv.cton_leaf_opt) {

				// Nitrogen stressed photosynthesis is allowed only if optimal leaf nitrogen is above allowed level 
				indiv.nstress = indiv.n_opt_isabovelim;
			}
		}
		else {

			// individual excess nitrogen demand after uptake
			double excess_ndemand = (indiv.leafndemand + indiv.rootndemand) * (1.0 - indiv.fnuptake)
				+ indiv.leafndemand_store + indiv.rootndemand_store;

			// if individual is in need of using its labile nitrogen storage
			if (!negligible(excess_ndemand)) {

				// if labile nitrogen storage is larger than excess nitrogen demand
				if (excess_ndemand <= indiv.nstore_labile) {

					// leaf nitrogen demand
					double leaf_ndemand = indiv.leafndemand * (1.0 - indiv.fnuptake) + indiv.leafndemand_store;
					indiv.nmass_leaf += leaf_ndemand;
					assert(indiv.nmass_leaf >= 0);
					indiv.nstore_labile -= leaf_ndemand;

					// root nitrogen demand
					double root_ndemand = indiv.rootndemand * (1.0 - indiv.fnuptake) + indiv.rootndemand_store;
					indiv.nmass_root += root_ndemand;
					indiv.nstore_labile -= root_ndemand;

					// nitrogen stressed photosynthesis is allowed only if optimal leaf nitrogen is above allowed level 
					indiv.nstress = indiv.n_opt_isabovelim;
				}
				else {

					if (!negligible(indiv.nstore_labile)) {

						// calculate total nitrogen mass
						double tot_nmass = indiv.nmass_leaf + indiv.nmass_root + indiv.fnuptake * (indiv.leafndemand + indiv.rootndemand) + indiv.nstore_labile;

						// new leaf C:N ratio
						double cton_leaf = (indiv.cmass_leaf_today() + indiv.cmass_root_today() * (indiv.pft.cton_leaf_avr / indiv.pft.cton_root_avr)) / tot_nmass;

						// nitrogen added to leaf from storage
						double labile_nto_leaf = indiv.cmass_leaf_today() / cton_leaf - (indiv.nmass_leaf + indiv.fnuptake * indiv.leafndemand);

						// new leaf nitrogen
						indiv.nmass_leaf += labile_nto_leaf;
						assert(indiv.nmass_leaf >= 0);

						// new root nitrogen
						indiv.nmass_root += indiv.nstore_labile - labile_nto_leaf;

						indiv.nstore_labile = 0.0;
					}

					// nitrogen stressed photosynthesis is allowed only when nitrogen limitation is turned on
					indiv.nstress = ifnlim && date.year > freenyears;

				}
			}
			else {
				// photosynthesis will not be nitrogen stresses unless optimal leaf N is above maximum limit
				indiv.nstress = indiv.n_opt_isabovelim;
			}
		}

		vegetation.nextobj();
	}
}

/// Nitrogen demand
/** Determines nitrogen demand based on vmax for leaves.
 *  Roots and sap wood nitrogen concentration follows leaf
 *  nitrogen concentration.
 *  Also determines individual nitrogen uptake capability
 */
void ndemand(Patch& patch, Vegetation& vegetation) {

	Gridcell& gridcell = patch.stand.get_gridcell();
	Soil& soil = patch.soil;

	/// daily nitrogen demand for patch (kgN/m2)
	patch.ndemand = 0.0;

	// soil available mineral nitrogen (kgN/m2)
	const double nmin_avail = soil.nmass_avail(NO);
	// Scalar to soil temperature (Eqn A9, Comins & McMurtrie 1993) for nitrogen uptake
	double soilT = patch.soil.get_soil_temp_25();
	double temp_scale = temperature_modifier(soilT); 
	
	/// Rate of nitrogen uptake not associated with Michaelis-Menten Kinetics (Zaehle and Friend 2010)
	double kNmin = 0.05;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// Rescaler of nitrogen uptake
		indiv.fnuptake = 1.0;

		// Assume that optimal leaf nitrogen isn't above allowed limit
		indiv.n_opt_isabovelim = false;

		// Starts with no nitrogen stress
		indiv.nstress = false;

		// Optimal leaf nitrogen content
		double leafoptn;

		// If grass with daily allocation
		bool daily_grass = indiv.pft.lifeform == GRASS && ifdailygrass;

		// Calculate optimal leaf nitrogen content and demand
		if ((!negligible(indiv.phen) && !daily_grass) || (!negligible(indiv.cmass_leaf_today()) && daily_grass)) {

			indiv.nday_leafon++;

			// Added a scalar depending on individual lai to slow down light optimization of newly shaded leaves
			// Peltoniemi et al. 2012
			indiv.nextin = exp(0.12 * min(10.0*indiv.phen, indiv.lai_indiv_today()));

			// Calculate optimal leaf nitrogen associated with photosynthesis and none photosynthetic
			// active nitrogen (Haxeltine et al. 1996 eqn 27/28)

			leafoptn = indiv.photosynthesis.nactive_opt * indiv.nextin + N0 * indiv.cmass_leaf_today();

			// Can not have higher nitrogen concentration than minimum leaf C:N ratio
			if (indiv.cmass_leaf_today() / leafoptn < indiv.pft.cton_leaf_min) {
				leafoptn = indiv.cmass_leaf_today() / indiv.pft.cton_leaf_min;

				// Optimal leaf N above limit -> always N limitation on Vmax
				indiv.n_opt_isabovelim = ifnlim && date.year > freenyears;
			}
			// Can not have lower nitrogen concentration than maximum leaf C:N ratio
			else if (indiv.cmass_leaf_today() / leafoptn > indiv.pft.cton_leaf_max) {
				leafoptn = indiv.cmass_leaf_today() / indiv.pft.cton_leaf_max;
			}

			// Updating annual optimal leaf C:N ratio
			indiv.cton_leaf_aopt = min(indiv.cmass_leaf_today() / leafoptn, indiv.cton_leaf_aopt);

			// Leaf nitrogen demand
			indiv.leafndemand = !daily_grass ? max(leafoptn - indiv.get_nmass_leaf(), 0.0) : 0.0;

			// Setting daily optimal leaf C:N ratio
			if (indiv.leafndemand || daily_grass) {
				indiv.cton_leaf_opt = indiv.cmass_leaf_today() / leafoptn;
			}
			else {
				indiv.cton_leaf_opt = max(indiv.pft.cton_leaf_min, indiv.cton_leaf());
			}
		}
		else {
			indiv.leafndemand = 0.0;
			indiv.cton_leaf_opt = indiv.cton_leaf();
		}

		// Nitrogen demand

		// Total nitrogen demand of individual
		double ndemand_tot;

		// Annual growth
		if (!daily_grass) {

			// Root nitrogen demand
			indiv.rootndemand = max(0.0, indiv.cmass_root_today() / (indiv.cton_leaf_opt * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr) - indiv.nmass_root);

			// Sap wood nitrogen demand. Demand is ramped up throughout the year.
			if (indiv.pft.lifeform == TREE) {
				indiv.sapndemand = max(0.0, indiv.cmass_sap / (indiv.cton_leaf_opt * indiv.pft.cton_sap_avr / indiv.pft.cton_leaf_avr) - indiv.nmass_sap) * ((1.0 + (double)date.day) / date.year_length());
			}

			// Labile nitrogen storage demand
			indiv.storendemand = indiv.ndemand_storage(indiv.cton_leaf_opt);

			//TODO HO demand
			indiv.hondemand = 0.0;

			// Total nitrogen demand
			ndemand_tot = max(0.0, indiv.leafndemand + indiv.rootndemand + indiv.sapndemand + indiv.storendemand + indiv.hondemand);
		}
		// Daily growth
		else {
			ndemand_tot = indiv.ndemand_storage(indiv.cton_leaf_opt);
		}

		// Calculate scalars to possible nitrogen uptake

		// Current plant mobile nitrogen concentration
		double ntoc = !negligible(indiv.cmass_leaf_today() + indiv.cmass_root_today()) ? (indiv.get_nmass_leaf() + indiv.nmass_root) / (indiv.cmass_leaf_today() + indiv.cmass_root_today()) : 0.0;

		// Scale to maximum nitrogen concentrations
		indiv.cton_status = max(0.0, (ntoc - 1.0 / indiv.pft.cton_leaf_min) / (1.0 / indiv.pft.cton_leaf_avr - 1.0 / indiv.pft.cton_leaf_min));

		// Nitrogen availablilty scalar due to saturating Michealis-Menten kinetics
		double nmin_scale = kNmin + nmin_avail / (nmin_avail + gridcell.pft[indiv.pft.id].Km);

		// Maximum available soil mineral nitrogen for this individual is base on its root area.
		// This is considered to be related to FPC which is proportional to crown area which is approx
		// 4 times smaller than the root area
		double max_indiv_avail = min(1.0, indiv.fpc * 4.0) * nmin_avail;

		// Maximum nitrogen uptake due to all scalars (times 2 because considering both NO3- and NH4+ uptake)
		// and soil available nitrogen within individual projectived coverage
		double maxnup = min(2.0 * indiv.pft.nuptoroot * nmin_scale * temp_scale * indiv.cton_status * indiv.cmass_root_today(), max_indiv_avail);

		// Nitrogen demand limitation due to maximum nitrogen uptake capacity
		double fractomax = ndemand_tot > 0.0 ? min(maxnup/ndemand_tot,1.0) : 0.0;

		// Annual growth
		if (!daily_grass) {

			// Root and leaf demand from storage pools
			indiv.leafndemand_store = indiv.leafndemand * (1.0 - fractomax);
			indiv.rootndemand_store = indiv.rootndemand * (1.0 - fractomax);

			// Nitrogen demand after adjustment to maximum uptake capacity
			indiv.leafndemand *= fractomax;
			indiv.rootndemand *= fractomax;
			indiv.sapndemand *= fractomax;
			indiv.storendemand *= fractomax;

			// Sum total nitrogen demand individual is capable of taking up
			indiv.ndemand = indiv.leafndemand + indiv.rootndemand + indiv.sapndemand + indiv.storendemand;

			// Negative nitrogen demand not allowed
			if (indiv.ndemand <= 0.0) {
				indiv.ndemand = 0.0;

				// Compartments fraction of total nitrogen demand
				indiv.leaffndemand = 0.0;
				indiv.rootfndemand = 0.0;
				indiv.sapfndemand = 0.0;
				indiv.storefndemand = 0.0;
			}
			else {

				// Compartments fraction of total nitrogen demand
				if (indiv.storendemand >= 0.0) {
					indiv.leaffndemand = indiv.leafndemand / indiv.ndemand;
					indiv.rootfndemand = indiv.rootndemand / indiv.ndemand;
					indiv.sapfndemand = indiv.sapndemand / indiv.ndemand;
					indiv.storefndemand = max(0.0, 1.0 - (indiv.leaffndemand + indiv.rootfndemand + indiv.sapfndemand));
				}
				else {
					indiv.leaffndemand = indiv.leafndemand / (indiv.leafndemand + indiv.rootndemand);
					indiv.rootfndemand = indiv.rootndemand / (indiv.leafndemand + indiv.rootndemand);
					indiv.sapfndemand = 0.0;
					indiv.storefndemand = 0.0;
				}
			}
		}
		// Daily growth
		else {
			indiv.ndemand = max(0.0, min(maxnup, ndemand_tot));
		}

		// Sum total patch nitrogen demand
		patch.ndemand += indiv.ndemand;

		vegetation.nextobj();
	}
}

/// Recalculation of photosynthesis under vmax nitrogen stress
/** If nitrogen supply is not able to meet demand it will lead
 *  to down-regulation of vmax resulting in lower photosynthesis
 */
void vmax_nitrogen_stress(Patch& patch, Climate& climate, Vegetation& vegetation) {

	// Supply function for nitrogen and determination of nitrogen stress leading
	// to down-regulation of vmax.

	// Nitrogen within projective cover of all individuals
	double tot_nmass_avail = patch.soil.nmass_avail(NO) * min(1.0, patch.fpc_total);

	// Calculate individual uptake fraction of nitrogen demand
	if (patch.ndemand > tot_nmass_avail) {

		// Determine individual nitrogen uptake fractions
		fnuptake(vegetation, tot_nmass_avail);
	}

	// Resolve nitrogen stress with longterm stored nitrogen
	nstore_usage(vegetation);

	// Calculate leaf nitrogen associated with photosynthesis, nitrogen limited photosynthesis,
	// and annual otimal leaf nitrogen content and nitrogen limitation on vmax

	vegetation.firstobj();

	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[pft.id];

		// Calculate leaf nitrogen associated with photosynthesis (Haxeltine et al. 1996 eqn 27/28)
		// Added difference between needleleaved and broadleaved mentioned in Friend et al. 1997
		// Should be done on FPC basis, but is not as it does not matter mathematically
		// Needs to be calculated for each individual due to possible water stress

		if (indiv.pft.lifeform == GRASS && ifdailygrass) {
			indiv.nactive = max(0.0, indiv.get_nmass_leaf() - N0 * indiv.get_cmass_leaf());
		}
		else {
			// Todays leaf nitrogen after uptake
			double nmass_leaf = indiv.get_nmass_leaf() + indiv.leafndemand * indiv.fnuptake;

			if (indiv.phen > 0.0) {
				indiv.nactive = max(0.0, nmass_leaf - N0 * indiv.cmass_leaf_today());
			}
			else {
				indiv.nactive = 0.0;
			}
		}

		// Individuals photosynthesis is nitrogen stressed
		if (indiv.nstress) {

			double pftco2 = get_co2(patch, climate, pft);
			PhotosynthesisEnvironment ps_env;
			ps_env.set(pftco2, climate.temp, climate.par, indiv.fpar, climate.daylength);

			// Set stresses
			PhotosynthesisStresses ps_stress;
			ps_stress.set(true, get_moss_wtp_limit(patch, pft), get_graminoid_wtp_limit(patch, pft), get_inund_stress(patch, ppft));

			// Individual photosynthesis
			photosynthesis(ps_env, ps_stress, pft,
						   pft.lambda_max, indiv.nactive / indiv.nextin, -1, 
						   indiv.photosynthesis);

			indiv.gpterm = gpterm(indiv.photosynthesis.adtmm, pftco2, pft.lambda_max, climate.daylength);

			if (date.diurnal()) {
				for (int i=0; i<date.subdaily; i++) {
					PhotosynthesisResult& ps_result = indiv.phots[i];

					ps_env.set(pftco2, climate.temps[i], climate.pars[i], indiv.fpar, 24);

					photosynthesis(ps_env, ps_stress, pft,
								   pft.lambda_max, indiv.nactive / indiv.nextin, indiv.photosynthesis.vm, 
								   ps_result);

					indiv.gpterms[i] = gpterm(ps_result.adtmm, climate.co2, pft.lambda_max, 24);
				}
			}
		}

		// Sum annual average nitrogen limitation on vmax
		if (indiv.phen)
			indiv.avmaxnlim += indiv.photosynthesis.vmaxnlim;

		// On last day of year determined average nitrogen limitation
		// based on days with leaf on
		if (date.islastday && date.islastmonth) {

			if (!negligible(indiv.nday_leafon)) {
				indiv.avmaxnlim /= (double)indiv.nday_leafon;
			}
			else {
				indiv.avmaxnlim = 0.0;
			}
		}
		vegetation.nextobj();
	}
}

/// Transpirative demand
/** Two alternative parameterisations of aet_monteith are available:
 *      AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL
 *  \see canexch.h
 */
void wdemand(Patch& patch, Climate& climate, Vegetation& vegetation, const Day& day) {

	// Determination of transpirative demand based on a Monteith parameterisation of
	// boundary layer dynamics, i.e. demand = f(EET, conductance)

	double gp_patch = 0.0;
		// non-water-stressed canopy conductance for patch, patch vegetated area
		// basis (mm/s)
	double gp_leafon_patch = 0.0;
		// non-water-stressed canopy conductance assuming full leaf cover, patch
		// vegetated area basis (mm/s)

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		Pft& pft = indiv.pft;

		// Calculate non-water-stressed canopy conductance assuming full leaf cover
		//        - include canopy-conductance component not linked to
		//          photosynthesis (diffusion through leaf cuticle etc); this is
		//          assumed to be proportional to leaf-on fraction

		// Call photosynthesis for individual assuming stomates fully open
		// (lambda = lambda_max)

		if (indiv.growingseason()) {

			PhotosynthesisResult leafon_photosynthesis;

			// Call photosynthesis first with fpar_leafon to get gp_leafon below.
			// Should hopefully not be needed in future, demand_leafon only used
			// by raingreen phenology.

			double temp = date.diurnal() ? climate.temps[day.period] : climate.temp;
			double par = date.diurnal() ? climate.pars[day.period] : climate.par;
			double daylength = date.diurnal() ? 24 : climate.daylength;
			double pftco2 = get_co2(patch, climate, pft);

			PhotosynthesisEnvironment ps_env;
			ps_env.set(pftco2, temp, par, indiv.fpar_leafon, daylength);

			PhotosynthesisStresses ps_stress;
			Patchpft& ppft = patch.pft[indiv.pft.id];
			ps_stress.set(false, get_moss_wtp_limit(patch, pft), get_graminoid_wtp_limit(patch, pft), get_inund_stress(patch, ppft));  

			// No nitrogen limitation when calculating gp_leafon
			photosynthesis(ps_env, ps_stress, pft,
			               pft.lambda_max, 1.0, -1, 
			               leafon_photosynthesis);

			double gp_leafon = gpterm(leafon_photosynthesis.adtmm, pftco2, pft.lambda_max, daylength) + pft.gmin * indiv.fpc;

			// Increment patch sums of non-water-stressed gp by individual value
			gp_patch +=  (date.diurnal() ? indiv.gpterms[day.period] : indiv.gpterm) + pft.gmin * indiv.fpc_today();
			gp_leafon_patch += gp_leafon;
		}

		vegetation.nextobj();
	}

	// Calculate transpirational demand on patch vegetated area basis
	// Eqn 23, Haxeltine & Prentice 1996
	if (!negligible(gp_patch) && !negligible(patch.fpc_total)) {
		gp_patch /= patch.fpc_total;
		patch.wdemand = aet_monteith(patch.eet_net_veg, gp_patch);
	}
	else
		patch.wdemand = 0.0;

	patch.wdemand_day += patch.wdemand;
	if (day.isend) {
		patch.wdemand_day /= date.subdaily;
	}

	if (!negligible(gp_leafon_patch) && !negligible(patch.fpc_total)) {
		gp_leafon_patch /= patch.fpc_total;
		patch.wdemand_leafon = aet_monteith(patch.eet_net_veg, gp_leafon_patch);
	}
	else
		patch.wdemand_leafon = 0.0;
}

/// Plant water uptake
/**
 * Returns plant water uptake (point scale, or mean for patch) as a fraction of
 * maximum possible (daily basis).
 *
 * Supports alternative parameterisations of plant water uptake:
 *
 * WCONT           = uptake rate coupled to water content and vertical
 *                   root distribution (as in earlier versions of LPJ-GUESS and LPJF)
 * ROOTDIST        = uptake rate independent of water content (to wilting point)
 *                   but with fractional uptake from different layers according
 *                   to prescribed root distribution
 * SMART           = uptake rate independent of water content (to wilting point),
 *                   fractional uptake from different layers according to layer
 *                   water content for trees, according to prescribed root
 *                   distribution for grasses
 * SPECIESSPECIFIC = uptake rate is species specific, with more drought
 *                   tolerance species (lower species_drought_tolerance values)
 *                   having greater relative uptake rates.
 */
double water_uptake(double wcont[NSOILLAYER], double awc[NSOILLAYER],
	double rootdist[NSOILLAYER], double emax, double fpc_rescale,
	double fwuptake[NSOILLAYER], bool ifsmart, double species_drought_tolerance) {
	// INPUT PARAMETERS:
	//   wcont       = water content of soil layers as fraction between wilting point
	//                 (0) and available water holding capacity (1)
	//   awc         = available water holding capacity of each soil layer (mm)
	//   rootdist    = plant root distribution (fraction in each soil layer)
	//   emax        = maximum evapotranspiration rate (mm/day)
	//   fpc_rescale = scaling factor for foliar projective cover (complement of patch
	//                 summed FPC overlap)
	//   ifsmart     = whether plants can freely adapt root profile to distribution of
	//                 available water among layers (required for "smart" mode)
	//   species_drought_tolerance = used only if the SPECIESSPECIFIC option is specified.


	// OUTPUT PARAMETER:
	//   fwuptake     = fraction of total uptake originating from each layer

	double wr=0.0;
	int s;

	switch (wateruptake) {
	case WR_WCONT:

		// LPJ "standard" formulation with linear scaling of uptake to water content
		// and weighting by plant root profiles

		for (s=0; s<NSOILLAYER; s++) {
			fwuptake[s] = rootdist[s] * wcont[s] * fpc_rescale;
			wr += fwuptake[s];
		}
		break;

	// guess2008 - drought/water uptake changes - new option
	case WR_SPECIESSPECIFIC:

		// Uptake rate is species specific, with more drought tolerance species (lower species_drought_tolerance
		// values) having greater relative uptake rates.
		// Reduces to WCONT if species_drought_tolerance = 0.5

		for (s=0; s<NSOILLAYER; s++) {
			double max_rel_uptake = pow(wcont[s], 2.0 * 0.1); // Upper limit. Limits C3 grass uptake
			fwuptake[s] = rootdist[s] * min(pow(wcont[s], 2.0 * species_drought_tolerance), max_rel_uptake) * fpc_rescale;
			wr += fwuptake[s];
		}
		break;
	case WR_ROOTDIST:

		// Uptake rate independent of water content (to wilting point) but with fractional
		// uptake from different layers according to prescribed root distribution

		for (s=0; s<NSOILLAYER; s++) {
			fwuptake[s] = min(wcont[s] * awc[s] * fpc_rescale, emax * rootdist[s]) / emax;
			wr += fwuptake[s];
		}
		break;
	case WR_SMART:
		{

			// Uptake rate independent of water content (to wilting point), fractional uptake
			// from different layers according to layer water content for trees, and according
			// to prescribed root distribution for grasses

			double wcsum = 0.0;
			double wcfrac;

			for (s=0; s<NSOILLAYER; s++) wcsum += wcont[s];

			if (negligible(wcsum))
				for (s=0; s<NSOILLAYER; s++) fwuptake[s] = 0.0;
			else {
				for (s=0; s<NSOILLAYER; s++) {
					wcfrac = wcont[s] / wcsum;
					if (ifsmart)
						fwuptake[s] = min(wcont[s] * awc[s] * wcfrac * fpc_rescale, emax * wcfrac) / emax;
					else
						fwuptake[s] = min(wcont[s] * awc[s] * fpc_rescale, emax * rootdist[s]) / emax;
					wr += fwuptake[s];
				}
			}
		}
		break;
	default:
		// Should never happen
		fail("Unsupported wateruptake type");
	}

	if (!negligible(wr))
		for (s=0; s<NSOILLAYER; s++)
			fwuptake[s] /= wr;

	return wr;
}


/// Plant water uptake (two layer)
/**
 * Returns plant water uptake (point scale, or mean for patch) as a fraction of
 * maximum possible (daily basis).
 *
 * Supports alternative parameterisations of plant water uptake:
 *
 * WCONT           = uptake rate coupled to water content and vertical
 *                   root distribution (as in earlier versions of LPJ-GUESS and LPJF)
 * ROOTDIST        = uptake rate independent of water content (to wilting point)
 *                   but with fractional uptake from different layers according
 *                   to prescribed root distribution
 * SMART           = uptake rate independent of water content (to wilting point),
 *                   fractional uptake from different layers according to layer
 *                   water content for trees, according to prescribed root
 *                   distribution for grasses
 * SPECIESSPECIFIC = uptake rate is species specific, with more drought
 *                   tolerance species (lower species_drought_tolerance values)
 *                   having greater relative uptake rates.
 */
double water_uptake_twolayer(double wcont[NSOILLAYER], double awc[NSOILLAYER],
	double rootdist[NSOILLAYER], double emax, double fpc_rescale,
	double fwuptake[NSOILLAYER], bool ifsmart, double species_drought_tolerance) {
	// INPUT PARAMETERS:
	//   wcont       = water content of soil layers as fraction between wilting point
	//                 (0) and available water holding capacity (1)
	//   awc         = available water holding capacity of each soil layer (mm)
	//   rootdist    = plant root distribution (fraction in each soil layer)
	//   emax        = maximum evapotranspiration rate (mm/day)
	//   fpc_rescale = scaling factor for foliar projective cover (complement of patch
	//                 summed FPC overlap)
	//   ifsmart     = whether plants can freely adapt root profile to distribution of
	//                 available water among layers (required for "smart" mode)
	//   species_drought_tolerance = used only if the SPECIESSPECIFIC option is specified.


	// OUTPUT PARAMETER:
	//   fwuptake     = fraction of total uptake originating from each layer

	double wr=0.0;
	int s;

	double gawc[2] = { 0.0, 0.0 };				// AWC for each Gerten-type layer 
	double grootdist[2] = { 0.0, 0.0 };			// rootdist for each Gerten-type layer
	double gwcont[2] = { 0.0, 0.0 };			// wcont for each Gerten-type soil layer
	double gfwuptake[2] = { 0.0, 0.0 };			// fwuptake for each Gerten-type soil layer

	double cap_upper = 0.0;
	double cap_lower = 0.0;

	for (int i=0;i<NSOILLAYER;i++) {
		if (i < NSOILLAYER_UPPER) {
			gawc[0] += awc[i];
			grootdist[0] += rootdist[i];
			gwcont[0] += wcont[i] * awc[i];
			cap_upper += awc[i];
		}
		else {
			gawc[1] += awc[i];
			grootdist[1] += rootdist[i];
			gwcont[1] += wcont[i] * awc[i];
			cap_lower += awc[i];
		}
	}

	gwcont[0] /= cap_upper;
	gwcont[1] /= cap_lower;


	switch (wateruptake) {
	case WR_WCONT:

		// LPJ "standard" formulation with linear scaling of uptake to water content
		// and weighting by plant root profiles

		for (s=0; s<NSOILLAYER_SIMPLE; s++) {
			gfwuptake[s] = grootdist[s] * gwcont[s] * fpc_rescale;
			wr += gfwuptake[s];
		}
		break;

	// guess2008 - drought/water uptake changes - new option
	case WR_SPECIESSPECIFIC:

		// Uptake rate is species specific, with more drought tolerance species (lower species_drought_tolerance
		// values) having greater relative uptake rates.
		// Reduces to WCONT if species_drought_tolerance = 0.5

		for (s=0; s<NSOILLAYER_SIMPLE; s++) {
			double max_rel_uptake = pow(gwcont[s], 2.0 * 0.1); // Upper limit. Limits C3 grass uptake
			gfwuptake[s] = grootdist[s] * min(pow(gwcont[s], 2.0 * species_drought_tolerance), max_rel_uptake) * fpc_rescale;
			wr += gfwuptake[s];
		}
		break;
	case WR_ROOTDIST:

		// Uptake rate independent of water content (to wilting point) but with fractional
		// uptake from different layers according to prescribed root distribution

		for (s=0; s<NSOILLAYER_SIMPLE; s++) {
			gfwuptake[s] = min(gwcont[s] * gawc[s] * fpc_rescale, emax * grootdist[s]) / emax;
			wr += gfwuptake[s];
		}
		break;
	case WR_SMART:
		{

			// Uptake rate independent of water content (to wilting point), fractional uptake
			// from different layers according to layer water content for trees, and according
			// to prescribed root distribution for grasses

			double wcsum = 0.0;
			double wcfrac;

			for (s=0; s<NSOILLAYER_SIMPLE; s++) wcsum += gwcont[s];

			if (negligible(wcsum))
				for (s=0; s<NSOILLAYER_SIMPLE; s++) gfwuptake[s] = 0.0;
			else {
				for (s=0; s<NSOILLAYER_SIMPLE; s++) {
					wcfrac = gwcont[s] / wcsum;
					if (ifsmart)
						gfwuptake[s] = min(gwcont[s] * gawc[s] * wcfrac * fpc_rescale, emax * wcfrac) / emax;
					else
						gfwuptake[s] = min(gwcont[s] * gawc[s] * fpc_rescale, emax * grootdist[s]) / emax;
					wr += gfwuptake[s];
				}
			}
		}
		break;
	default:
		// Should never happen
		fail("Unsupported wateruptake type");
	}

	if (!negligible(wr))
		for (s=0; s<NSOILLAYER_SIMPLE; s++)
			gfwuptake[s] /= wr;


	// Ensure equal water uptake from each Gerten-type layer is set in the finer layers
	for (int i=0;i<NSOILLAYER;i++) {
		if (i < NSOILLAYER_UPPER)
			fwuptake[i] = gfwuptake[0] / NSOILLAYER_UPPER;
		else 
			fwuptake[i] = gfwuptake[1] / NSOILLAYER_LOWER;
	}

	return wr;
}



/// Plant water uptake for irrigated crops
/**
 * Returns plant water uptake (point scale, or mean for patch) as a fraction of
 * maximum possible (daily basis), after adding required water to obtain maximum
 * water uptake.
 * Irrigation water is added to the soil in hydrology_lpjf
 *
 * Only ROOTDIST currently supported plant water uptake parameterisation:
 *
 * ROOTDIST        = uptake rate independent of water content (to wilting point)
 *                   but with fractional uptake from different layers according
 *                   to prescribed root distribution
 */
double irrigated_water_uptake(Patch& patch, Pft& pft, const Day& day) {
	Patchpft& ppft = patch.pft[pft.id];
	
	// Prepare for water balance tests
	double initial_water_in_column = 0.0;

	for (int ly = 0; ly < NSOILLAYER; ly++) {

		double Faw_layer = patch.soil.get_layer_soil_water(ly) * patch.soil.soiltype.awc[ly]; // mm
		double ice_layer = patch.soil.Frac_ice[ly + patch.soil.IDX] * patch.soil.Dz[ly + patch.soil.IDX]; // mm
		double layerwater = Faw_layer + ice_layer; // mm
		// for water balance checks
		initial_water_in_column += layerwater;

	} // for loop (ly)



	double awc0 = patch.soil.soiltype.gawc[0]; // available water holding capacity; see Gerten et al., 2004  
	double awc1 = patch.soil.soiltype.gawc[1]; // see Gerten et al., 2004
	double grootdist[2] = { 0.0, 0.0 };

	double wcont_cp[NSOILLAYER];
	for (int i=0;i<NSOILLAYER;i++) {
		wcont_cp[i] = patch.soil.get_layer_soil_water(i);
		if (i < NSOILLAYER_UPPER)
			grootdist[0] += pft.rootdist[i];
		else
			grootdist[1] += pft.rootdist[i];
	}

	if (day.isstart) {
		ppft.water_deficit_d = 0.0;
		if (date.day == 0) {
			ppft.water_deficit_y = 0.0;
		}
	}

	double wcont_opt[NSOILLAYER_UPPER];
	bool add_water[NSOILLAYER_UPPER];
	for (int i = 0; i<NSOILLAYER_UPPER; i++) add_water[i] = true;

	if (patch.soil.get_soil_water_upper()<0.9 && ppft.phen > 0.0)	{
		double wcont_0_opt = 0.0;
		double wr_opt = min(1.0, patch.wdemand / ppft.phen / pft.emax);

		if (wateruptake == WR_ROOTDIST) {

			wcont_0_opt = (wr_opt * pft.emax - min(patch.soil.get_soil_water_lower() * awc1 * patch.fpc_rescale, pft.emax * grootdist[1])) / awc0 / patch.fpc_rescale;
			if (wcont_0_opt * awc0 * patch.fpc_rescale > pft.emax * grootdist[0]) {
				wcont_0_opt = pft.emax * grootdist[0] / awc0 / patch.fpc_rescale;
			}

			if (!iftwolayersoil) {

				// Array of optimal wcont_values instead.
				for (int i = 0; i<NSOILLAYER_UPPER; i++) {
					wcont_opt[i] = (wr_opt * pft.emax - min(patch.soil.get_soil_water_lower() * awc1 * patch.fpc_rescale, pft.emax * grootdist[1])/ NSOILLAYER_UPPER) / patch.soil.soiltype.awc[i] / patch.fpc_rescale;
			
					if (wcont_opt[i] * patch.soil.soiltype.awc[i] * patch.fpc_rescale > pft.emax * pft.rootdist[i]) {
						wcont_opt[i] = pft.emax * pft.rootdist[i] / patch.soil.soiltype.awc[i] / patch.fpc_rescale;
					}

					if (wcont_cp[i] > wcont_opt[i]) {
						// Do not change wcont for layers that already have enough available water, or add water to them
						wcont_opt[i] = wcont_cp[i];
						add_water[i] = false;
					}
				}

			}

		}
		else {
			fail("Irrigation soil water only balanced for WR_ROOTDIST currently !\n");
		}

		bool irrigate_soil = false;

		if (iftwolayersoil) {
			if (wcont_0_opt > patch.soil.get_soil_water_upper()) 
				irrigate_soil = true;
		}
		else {
			for (int i = 0; i < NSOILLAYER_UPPER; i++) {
				if (add_water[i])
					irrigate_soil = true; // Allows irrigation if even one layer needs more water 
			}
		}

		if (irrigate_soil && patch.soil.dsnowdepth <= 0.001) { // No irrigation when there is snow on the ground

			// No irrigation when there is ice left in the top 50cm
			if (!patch.soil.ice_in_top_layer()) {

				// Irrigation water for this PFT:

				double water_to_add = 0.0;

				double water_to_add_ly[NSOILLAYER_UPPER];
				for (int i = 0; i < NSOILLAYER_UPPER; i++) water_to_add_ly[i] = 0.0;

				if (iftwolayersoil)
					water_to_add = (wcont_0_opt - patch.soil.get_soil_water_upper()) * awc0;
				else {
					for (int i = 0; i<NSOILLAYER_UPPER; i++) {
						if (add_water[i]) {
							water_to_add_ly[i] = max((wcont_opt[i] - wcont_cp[i]) * patch.soil.soiltype.awc[i], 0.0); // ensures that we add water
							water_to_add += water_to_add_ly[i];
						}
					}
				}

				ppft.water_deficit_d += water_to_add;

				// Local storage

				// available water for each soil layer (mm)
				double Faw_layer[NSOILLAYER_UPPER];
				// water that can still be added to each soil layer (mm)
				double potential_layer[NSOILLAYER_UPPER];

				// initialise to 0 mm
				for (int sl = 0; sl<NSOILLAYER_UPPER; sl++) {
					Faw_layer[sl] = 0.0;
					potential_layer[sl] = 0.0;
				}

				Soil& soil = patch.soil;
				double potential_water = 0.0;

				for (int ly = 0; ly < NSOILLAYER_UPPER; ly++) {
					if (add_water[ly]) {
						// Use: wcont[ly - IDX] = Faw_layer[ly - IDX] / soiltype.awc[ly - IDX];
						Faw_layer[ly] = soil.get_layer_soil_water(ly) * soil.soiltype.awc[ly]; // mm

						// Water in this layer
						potential_layer[ly] = soil.aw_max[ly] - Faw_layer[ly];
						potential_water += potential_layer[ly];
					}
				} // for loop (ly)

				// to test the new wcont
				double total_available_water = 0.0;
				double total_capacity = 0.0;
	
				for (int s=0; s<NSOILLAYER_UPPER; s++) {

					if (add_water[s]) {

						double water_input_ly = 0.0;
						
						if (iftwolayersoil)
							water_input_ly = water_to_add * (potential_layer[s] / potential_water);
						else
							water_input_ly = water_to_add_ly[s];

						Faw_layer[s] += water_input_ly;
						potential_layer[s] -= water_input_ly;
						wcont_cp[s] = Faw_layer[s] / soil.soiltype.awc[s];
				
						// Update wcont (and wcont_evap, Frac_water etc.) too
						soil.set_layer_soil_water(s, Faw_layer[s] / soil.soiltype.awc[s]);

					}
					
					total_available_water += wcont_cp[s] * soil.soiltype.awc[s]; // mm
					total_capacity += soil.soiltype.awc[s]; // mm

					oob_check_wcont(wcont_cp[s]);

					// Error check
					if (!iftwolayersoil) {
						if (fabs(wcont_cp[s] - wcont_opt[s]) > 0.0001) {
							dprintf("irrigated_water_uptake - error in the water balance!\n");
							return -9999.0;
						}
					}

				}

				// Error check
				if (iftwolayersoil) {
					double new_wcont = total_available_water/ total_capacity;
			
					if (fabs(new_wcont - wcont_0_opt) > 0.0001) {
						dprintf("irrigated_water_uptake - error in the water balance!\n");
						return -9999.0;
					}
				}
			}
		}

		if (day.isend) {
			ppft.water_deficit_d /= date.subdaily;
			ppft.water_deficit_y += ppft.water_deficit_d;
		}
	}


	// Prepare for water balance test
	double final_water_in_column = 0.0;

	for (int ly = 0; ly < NSOILLAYER; ly++) {

		double Faw_layer = patch.soil.get_layer_soil_water(ly) * patch.soil.soiltype.awc[ly]; // mm
		double ice_layer = patch.soil.Frac_ice[ly + patch.soil.IDX] * patch.soil.Dz[ly + patch.soil.IDX]; // mm
		double layerwater = Faw_layer + ice_layer; // mm
		// for water balance checks
		final_water_in_column += layerwater;

	} // for loop (ly)

	// is water in + initial storage = water out + final storage?
	double water_in_storage_in = initial_water_in_column + ppft.water_deficit_d;
	double water_out_storage_out = final_water_in_column;

	const double maxerr = 0.0001; // mm - max error allowed

	if (fabs(water_in_storage_in - water_out_storage_out) > maxerr) {
		fail("irrigated_water_uptake - error in the water balance!\n");
	}


	if (iftwolayersoil) 
		return water_uptake_twolayer(wcont_cp, patch.soil.soiltype.awc, pft.rootdist, pft.emax, patch.fpc_rescale,
				ppft.fwuptake, pft.lifeform == TREE, pft.drought_tolerance);
	else
		return water_uptake(wcont_cp, patch.soil.soiltype.awc, pft.rootdist, pft.emax, patch.fpc_rescale,
				ppft.fwuptake, pft.lifeform == TREE, pft.drought_tolerance); // i.e. use upland awc instead of awc_peat - never irrigate on peatlands

};


/// Actual evapotranspiration and water stress
/** Soil water supply at the roots available to meet the transpirational demand
 *  Fundamentally, water stress = supply < demand
 */
void aet_water_stress(Patch& patch, Vegetation& vegetation, const Day& day) {

	// Supply function for evapotranspiration and determination of water stress leading
	// to down-regulation of stomatal conductance. Actual evapotranspiration determined
	// as smaller of supply and transpirative demand (see function demand).
	// Base value for actual canopy conductance calculated here for water-stressed
	// individuals and used to derive actual photosynthesis in function npp (below)

	// Calculate common point supply for each PFT in this patch
	for (int p=0; p<npft; p++) {

		Standpft& spft = patch.stand.pft[p];
		if (!spft.active)
			continue;

		// Retrieve next patch PFT
		Patchpft& ppft = patch.pft[p];

		// Retrieve PFT
		Pft& pft = ppft.pft;

		if (day.isstart || patch.stand.isirrigated) {

			// Calculate effective water supply from plant roots
			// Rescale available water by patch FPC if exceeds 1
			// (this then represents the average amount of water available over an
			// individual's FPC, assuming individuals are equal in competition for water)
			double wr;
 
			if (patch.stand.isirrigated) {
				wr = irrigated_water_uptake(patch, pft, day);
			} 
			else {
				double wcont_local[NSOILLAYER];
				for (int i = 0; i<NSOILLAYER; i++) {
					wcont_local[i] = patch.soil.get_layer_soil_water(i);
				}

			if (iftwolayersoil) 
				wr = water_uptake_twolayer(wcont_local, patch.soil.soiltype.awc,
							pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
							pft.lifeform == TREE, pft.drought_tolerance);
			else {

				if (patch.stand.is_highlatitude_peatland_stand()) // Use awc_peat
					wr = water_uptake(wcont_local, patch.soil.soiltype.awc_peat,
							pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
							pft.lifeform == TREE, pft.drought_tolerance);
				else
					wr = water_uptake(wcont_local, patch.soil.soiltype.awc,
							pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
							pft.lifeform == TREE, pft.drought_tolerance);

				}
			}

			// Calculate supply (Eqn 24, Haxeltine & Prentice 1996)
			if (patch.stand.landcover != CROPLAND && !(ifdailygrass && ppft.pft.lifeform == GRASS)){
				ppft.wsupply_leafon = pft.emax * wr;
				ppft.wsupply = ppft.wsupply_leafon * ppft.phen;
			}
			else if (ifdailygrass && ppft.pft.lifeform == GRASS) {
				ppft.wsupply_leafon = pft.emax * wr;
				ppft.wsupply = ppft.wsupply_leafon;
			}
			else {  // crop specific water supply
				// Temporarily removed the scaling with fpc due to 
				// too high water stress for young crops. (r10394)
				// TODO: Make the water supply dependant on root allocation
				if (ppft.cropphen->growingseason) {
					ppft.wsupply_leafon = pft.emax * wr;
					ppft.wsupply = ppft.wsupply_leafon;
				}
				else {
					ppft.wsupply_leafon = 0.0;
					ppft.wsupply = 0.0;
				}
			}
		}

		if (date.get_calendar_year() == 1853 && date.day == 200)
			int sch = 0;

		ppft.wstress = ppft.wsupply < patch.wdemand && !negligible(ppft.phen) && !(pft.phenology==CROPGREEN && !largerthanzero(patch.wdemand-ppft.wsupply, -10));

		if (ifdailygrass && ppft.pft.lifeform == GRASS) {
			ppft.wstress = ppft.wsupply < patch.wdemand;
		}

		// Calculate water-stressed canopy conductance on FPC basis assuming
		// FPAR=1 and deducting canopy conductance component not associated
		// with CO2 uptake; valid for all individuals of this PFT in this patch
		// today.
		// Eqn 25, Haxeltine & Prentice 1996

		// Fix, valid for monocultures, for faulty equation, manifesting itself in problems with crops in high scenario CO2-levels.
		// No fix for natural vegetation yet.
		double gmin = pft.phenology==CROPGREEN ? ppft.phen * pft.gmin : pft.gmin;

		ppft.gcbase = ppft.wstress ? max(gc_monteith(ppft.wsupply, patch.eet_net_veg)-
					gmin * ppft.wsupply / patch.wdemand, 0.0) : 0;

		if (!date.diurnal()) {
			ppft.wstress_day = ppft.wstress;
			ppft.gcbase_day = ppft.gcbase;
		}
		else if (day.isend) {

			ppft.wstress_day = ppft.wsupply < patch.wdemand_day && !negligible(ppft.phen) && !(pft.phenology==CROPGREEN && !largerthanzero(patch.wdemand-ppft.wsupply, -10));

			ppft.gcbase_day = ppft.wstress_day ? max(gc_monteith(ppft.wsupply,
					patch.eet_net_veg) - gmin * ppft.wsupply / patch.wdemand_day, 0.0) : 0;
		}
	}

	// Calculate / transfer supply to individuals
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		Patchpft& ppft = patch.pft[indiv.pft.id];
		if (day.isstart) {
			indiv.aet = 0;

			if (date.day == 0)
				indiv.aaet = 0.0;
		}

		indiv.wstress = ppft.wstress;

		if (indiv.alive) {
			if (indiv.wstress) {
				indiv.aet += ppft.wsupply;
			}
			else {
				indiv.aet += negligible(indiv.phen) ? 0.0 : patch.wdemand;
			}
		}
		if (day.isend) {
			indiv.aet *= indiv.fpc / date.subdaily;
		}

		if (day.isend) {
			indiv.aaet += indiv.aet;
		}

		vegetation.nextobj();
	}
}

/// Water scalar
void water_scalar(Patch& patch, Vegetation& vegetation, const Day& day) {

	// Derivation of daily and annual versions of water scalar (wscal, omega)
	// Daily version is used to determine leaf onset and abscission for raingreen PFTs.
	// Annual version determines relative allocation to roots versus leaves for
	// subsequent year

	for (int p=0; p<npft; p++) {

		// Retrieve next patch PFT
		Patchpft& ppft = patch.pft[p];

		if (!patch.stand.pft[p].active)
			continue;

		if (day.isstart) {
			ppft.wscal = 0;
			if (date.day == 0) {
				ppft.wscal_mean = 0;
				if (ppft.pft.phenology==CROPGREEN || ppft.pft.isintercropgrass)
					ppft.cropphen->growingdays_y=0;
			}
		}

		// Calculate patch PFT water scalar value

		if (!negligible(patch.wdemand_leafon)) {
			ppft.wscal += min(1.0, ppft.wsupply_leafon/patch.wdemand_leafon);
		}
		else {
			ppft.wscal += 1.0;
		}

		if (day.isend) {
			ppft.wscal /= (double)date.subdaily;

			if (patch.stand.landcover!=CROPLAND										//natural, urban, pasture, forest and peatland stands
					|| ppft.pft.phenology==ANY && ppft.pft.id==patch.stand.pftid) {	//normal grass growth
				ppft.wscal_mean += ppft.wscal;

				// Convert from sum to mean on last day of year
				if (date.islastday && date.islastmonth) {
					ppft.wscal_mean /= date.year_length();
				}
			}
			else if (ppft.cropphen->growingseason									// true crops and cover-crop grass
					|| ppft.pft.phenology == CROPGREEN && date.day == ppft.cropphen->hdate
					|| ppft.pft.isintercropgrass && date.day == patch.pft[patch.stand.pftid].cropphen->eicdate) {

				ppft.cropphen->growingdays_y++;
				ppft.wscal_mean = ppft.wscal_mean + (ppft.wscal - ppft.wscal_mean) / ppft.cropphen->growingdays_y;
			}
		}

		ppft.wscal_5.add(ppft.wscal);
		ppft.wscal_10.add(ppft.wscal);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// ASSIMILATION_WSTRESS
// Internal function (do not call directly from framework)

void assimilation_wstress(const Pft& pft, double co2, double temp, double par,
			double daylength, double fpar, double fpc, double gcbase,
			double vmax, PhotosynthesisResult& phot_result, double& lambda,
			double nactive, bool ifnlimvmax, double moss_wtp_limit, double graminoid_wtp_limit, double inund_stress) {

	// DESCRIPTION
	// Calculation of net C-assimilation under water-stressed conditions
	// (demand>supply; see function canopy_exchange). Utilises a numerical
	// iteration procedure to find the level of stomatal aperture (characterised by
	// lambda, the ratio of leaf intercellular to ambient CO2 concentration) which
	// satisfies simulataneously a canopy-conductance based and light-based
	// formulation of photosynthesis (Eqns 2, 18 and 19, Haxeltine & Prentice (1996)).

	// Numerical method is a tailored implementation of the bisection method,
	// assuming root (f(lambda)=0) bracketed by f(0.02)<0 and
	// f(lambda_max)>0 (Press et al 1986)

	// The bisection method terminates when we're close enough to a root
	// (absolute value of f(lambda) < EPS), or after a maximum number of
	// iterations.

	// Note that the function sometimes doesn't search for a lambda,
	// and returns zero assimilation (for instance if there is no
	// root within the valid interval, or if daylength is zero).
	// So if zero assimilation is returned, the returned lambda should
	// not be used!

	// OUTPUT PARAMETER
	// phot_result = result of photosynthesis for the found lambda
	// lambda      = the lambda found by the bisection method (see above)


	// Set lambda to something for cases where we don't actually search for
	// a proper lambda. This value shouldn't be used (see documentation
	// above), but we'll set it to something anyway so we don't return
	// random garbage.
	lambda = -1;

	if (negligible(fpc) || negligible(fpar) || negligible(gcbase * daylength * 3600)) {
		// Return zero assimilation
		phot_result.clear();
		return;
	}

	// Canopy conductance component associated with photosynthesis on a
	// daily basis (mm / m2 / day)
	double gcphot = gcbase * daylength * 3600 / 1.6 * co2 * CO2_CONV;

	// At this point the function f(x) = g(x) - h(x) can be calculated as:
	//
	// g(x) = phot_result.adtmm / fpc (after a call to photosynthesis with lambda x)
	// h(x) = gcphot * (1 - x)

	// Evaluate f(lambda_max) to see if there's a root
	// in the interval we're searching

	PhotosynthesisEnvironment ps_env;
	ps_env.set(co2, temp, par, fpar, daylength);

	PhotosynthesisStresses ps_stress;
	ps_stress.set(ifnlimvmax, moss_wtp_limit, graminoid_wtp_limit, inund_stress);  

	photosynthesis(ps_env, ps_stress, pft, pft.lambda_max, nactive, vmax, phot_result);
	double f_lambda_max = phot_result.adtmm / fpc - gcphot * (1 - pft.lambda_max);

	if (f_lambda_max <= 0) {
		// Return zero assimilation
		phot_result.clear();
		return;
	}

	const double EPS = 0.1; // minimum precision of solution in bisection method

	double xmid = 0.0;

	// Implement numerical solution

	double x1 = 0.02;                      // minimum bracket of root
	double x2 = pft.lambda_max;            // maximum bracket of root
	double rtbis = x1;                     // root of the bisection
	double dx = x2 - x1;

	const int MAXTRIES = 6; // maximum number of iterations towards a solution
	int b = 0;              // number of tries so far towards solution

	double fmid = EPS + 1.0;

	ps_env.set(co2, temp, par, fpar, daylength);

	while (fabs(fmid) > EPS && b <= MAXTRIES) {

		b++;
		dx *= 0.5;
		xmid = rtbis + dx;				// current guess for lambda

		// Call function photosynthesis to calculate alternative value
		// for total daytime photosynthesis according to Eqns 2 & 19,
		// Haxeltine & Prentice (1996), and current guess for lambda

		photosynthesis(ps_env, ps_stress, pft, 
					   xmid, nactive, vmax, 
					   phot_result);

		// Evaluate fmid at the point lambda=xmid
		// fmid will be an increasing function of xmid, with a solution
		// (fmid=0) between x1 and x2

		// Second term is total daytime photosynthesis (mm/m2/day) implied by
		// canopy conductance and current guess for lambda (xmid)
		// Eqn 18, Haxeltine & Prentice 1996
		fmid = phot_result.adtmm / fpc - gcphot * (1 - xmid);

		if (fmid < 0) {
			rtbis = xmid;
		}
	}

	// bvoc
	lambda=xmid;
}

///////////////////////////////////////////////////////////////////////////////////////
// AUTOTROPHIC RESPIRATION
// Internal function (do not call directly from framework)

void respiration(double gtemp_air, double gtemp_soil, lifeformtype lifeform,
	double respcoeff, double cton_sap, double cton_root,
	double cmass_sap, double cmass_root_today, double growth_demand, double fphen_activity,
	double assim, double& resp) {

	// DESCRIPTION
	// Calculation of daily maintenance and growth respiration for individual with
	// specified life form, phenological state, tissue C:N ratios and daily net
	// assimilation, given current air and soil temperatures.
	// Sitch et al. (2000), Lloyd & Taylor (1994), Sprugel et al (1996).

	// NOTE: leaf respiration is not calculated here, but included in the calculation
	// of net assimilation (function production above) as a proportion of rubisco
	// capacity (Vmax).

	// INPUT PARAMETERS
	// gtemp_air  = respiration temperature response incorporating damping of Q10
	//              response due to temperature acclimation (Eqn 11, Lloyd & Taylor
	//              1994); Eqn B2 below
	// gtemp_soil = as gtemp_air given soil temperature
	// lifeform   = PFT life form class (TREE or GRASS)
	// respcoeff  = PFT respiration coefficient
	// cton_sap   = PFT sapwood C:N ratio
	// cton_root  = PFT root C:N ratio
	// phen       = vegetation phenological state (fraction of potential leaf cover)
	// cmass_sap  = sapwood C biomass on grid cell area basis (kgC/m2)
	// cmass_root = fine root C biomass on grid cell area basis (kgC/m2)
	// assim      = net assimilation on grid cell area basis (kgC/m2/day)

	// OUTPUT PARAMETER
	// resp       = sum of maintenance and growth respiration on grid cell area basis
	//              (kgC/m2/day)

	// guess2008 - following a comment by Annett Wolf, the following parameter value was changed:
	// const double K=0.0548; // OLD value
	const double K=0.095218;  // NEW parameter value in respiration equations
	// See the comment after Eqn (4) below.

	double resp_sap = 0;    // sapwood respiration (kg/m2/day)
	double resp_root = 0;   // root respiration (kg/m2/day)
	double resp_growth = 0; // growth respiration (kg/m2/day)

	// Calculation of maintenance respiration components for each living tissue:
	//
	// Based on the relations
	//
	// (A) Tissue respiration response to temperature
	//     (Sprugel et al. 1996, Eqn 7)
	//
	//     (A1) Rm = 7.4e-7 * N * f(T)
	//     (A2) f(T) = EXP (beta * T)
	//
	//       where Rm   = tissue maintenance respiration rate in mol C/sec
	//             N    = tissue nitrogen in mol N
	//             f(T) = temperature response function
	//             beta = ln Q10 / 10
	//             Q10  = change in respiration rate with a 10 K change
	//                    in temperature
	//             T    = tissue absolute temperature in K
	//
	// (B) Temperature response of soil respiration across ecosystems
	//     incorporating damping of Q10 response due to temperature acclimation
	//     (Lloyd & Taylor 1994, Eqn 11)
	//
	//     (B1) R = R10 * g(T)
	//     (B2) g(T) = EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))]
	//
	//       where R    = respiration rate
	//             R10  = respiration rate at 10 K
	//             g(T) = temperature response function at 10 deg C
	//             T    = soil absolute temperature in K
	//
	// Mathematical derivation:
	//
	// For a tissue with C:N mass ratio cton, and C mass, c_mass, N concentration
	// in mol given by
	//  (1) N = c_mass / cton / atomic_mass_N
	// Tissue respiration in gC/day given by
	//  (2) R = Rm * atomic_mass_C * seconds_per_day
	// From (A1), (1) and (2),
	//  (3) R = 7.4e-7 * c_mass / cton / atomic_mass_N * atomic_mass_C
	//          * seconds_per_day * f(T)
	// Let
	//  (4) k = 7.4e-7 * atomic_mass_C / atomic_mass_N * seconds_per_day
	//        = 0.0548

	// guess2008 - there is an ERROR here, spotted by Annett Wolf
	// If we calculate the respiration at 20 degC using g(T) and compare it to
	// Sprugel's eqn 3, for 1 mole tissue N, say, we do NOT get the same result with this
	// k value. This is because g(T) = 1 at 10 degC, not 20 degC. Changing k from 0.0548
	// to 0.095218 gives exactly the same results as Sprugel at 20 degC. The scaling factor
	// 7.4e-7 used here is taken from Sprugel's eqn. (7), but they used f(T), not g(T), and
	// these are defined on different bases.

	// from (3), (4)
	//  (5) R = k * c_mass / cton * f(T)
	// substituting ecosystem temperature response function g(T) for f(T) (Eqn B2),
	//  (6) R = k * c_mass / cton * g(T)
	// incorporate PFT-specific respiration coefficient to model acclimation
	// of respiration rates to average (temperature) conditions for PFT (Ryan 1991)
	//  (7) R_pft = respcoeff_pft * k * c_mass / cton * g(T)

	if (lifeform == TREE) {

		// Sapwood respiration (Eqn 7)

		resp_sap = respcoeff * K * cmass_sap / cton_sap * gtemp_air;

		// Root respiration (Eqn 7)
		// Assumed that root phenology follows leaf phenology

		resp_root = respcoeff * K * cmass_root_today / cton_root * gtemp_soil;

		// Growth respiration = 0.25 ( GPP - maintenance respiration)
		resp_growth = (assim - resp_sap - resp_root) * 0.25;

		// guess2008 - disallow negative growth respiration
		// (following a comment (060823) from Annett Wolf)
		if(resp_growth < 0.0) resp_growth = 0.0;

		// Total respiration is sum of maintenance and growth respiration

		resp = resp_sap + resp_root + resp_growth;
	}
	else if (lifeform == GRASS || lifeform == MOSS) {

		if (ifdavegrowth) {

			// Handle respiration differently in DAVE: 
			// -> scale down resp_root during dormancy
			// -> account for growth_demand in Rgrowth

			// Root respiration
			// Reduction factor of resp_root during dormancy
			resp_root = (respcoeff * K * cmass_root_today / cton_root * gtemp_soil) * fphen_activity;

			// Growth respiration
			const double Rg_max = 0.4;  // maximum Rg as fraction of GPP - Rm. TODO: adjust value
			resp_growth = (assim - resp_sap - resp_root) * (growth_demand * Rg_max);

		} else {
			
			// Root respiration
			resp_root = respcoeff * K * cmass_root_today / cton_root * gtemp_soil;
			
			// Growth respiration (see above)
			resp_growth = (assim - resp_root) * 0.25;
		}

		// guess2008 - disallow negative growth respiration
		// (following a comment (060823) from Annett Wolf)
		if(resp_growth < 0.0) resp_growth = 0.0;

		// Total respiration (see above)

		resp = resp_root + resp_growth;
	}
	else fail ("Canopy exchange function respiration: unknown life form");
}

/// Net Primary Productivity
/** Includes BVOC calculations \see bvoc.cpp
 */
void npp(Patch& patch, Climate& climate, Vegetation& vegetation, const Day& day) {

	// Determination of daily NPP. Leaf level net assimilation calculated for non-
	// water-stressed individuals (i.e. with fully-open stomata) using base value
	// from function demand (above); for water-stressed individuals using base value
	// for canopy conductance by a simultaneous solution of light-based and canopy
	// conductance-based equations for net daily photosynthesis (see function
	// assimilation wstress above). The latter uses the PFT-specific base value for
	// conductance from function aet_water_stress (above).
	// Plant respiration obtained by a call to function respiration (above).

	double par, temp, assim, resp, rad, gtemp;
	double hours = 24;			// diurnal "daylength" to convert to daily units

	if (date.diurnal()) {
		par   = climate.pars[day.period];
		temp  = climate.temps[day.period];
		rad   = climate.rads[day.period];
		gtemp = climate.gtemps[day.period];
	}
	else {
		par   = climate.par;
		temp  = climate.temp;
		hours = climate.daylength;
		rad   = climate.rad;
		gtemp = climate.gtemp;
	}

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// For this individual ...

		// Retrieve PFT and patch PFT

		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[pft.id];

		// Don't do calculations for crops outside their growingseason
		if (!indiv.growingseason()) {
			indiv.dnpp = 0.0;
			indiv.lambda = 0.0;
			vegetation.nextobj();
			continue;
		}

		double pftco2 = get_co2(patch, climate, pft);
		double inund_stress = get_inund_stress(patch, ppft);
		double graminoid_wtp_limit = get_graminoid_wtp_limit(patch, pft);
		double moss_wtp_limit = get_moss_wtp_limit(patch, pft);

		PhotosynthesisResult phot = date.diurnal() ? indiv.phots[day.period] : indiv.photosynthesis;

		if (indiv.wstress) {

			// Water stress - derive assimilation by simultaneous solution
			// of light- and conductance-based equations of photosynthesis

			assimilation_wstress(pft, pftco2, temp, par, hours, indiv.fpar, indiv.fpc,
				ppft.gcbase, phot.vm, phot, indiv.lambda,
				indiv.nactive / indiv.nextin, ifnlim, moss_wtp_limit, graminoid_wtp_limit, inund_stress);

			// Make sure that lambda isn't negative (-1) if something went wrong in the photosynthesis call
			indiv.lambda = max(0.0, indiv.lambda);
		}
		else {
			// No water limitation
			indiv.lambda = pft.lambda_max;
		}

		assim = phot.net_assimilation();

		if (ifbvoc) {
			PhotosynthesisResult phot_nostress = date.diurnal() ? indiv.phots[day.period] : indiv.photosynthesis;
			bvoc(temp, hours, rad, climate, patch, indiv, pft, phot_nostress, phot.adtmm, day);
		}

		// Calculate autotrophic respiration
		double cmass_root;
		if (indiv.cropindiv && indiv.cropindiv->isintercropgrass && indiv.phen == 0.0)
			cmass_root = 0.0;
		else
			cmass_root = indiv.cmass_root_today();

		// Static root and sap wood C:N ratio if no N limitation
		// to not let N affect respiration for C only version of model
		// Also not allow lower C:N ratios than minimum which
		// would affect respiration severely (happens for deciduous trees
		// in the beginning of the season)

		double cton_sap, cton_root;
		if (ifnlim) {
			cton_sap = max(indiv.cton_sap(), indiv.pft.cton_leaf_min * indiv.pft.cton_sap_avr / indiv.pft.cton_leaf_avr);
			cton_root = max(indiv.cton_root(), indiv.pft.cton_leaf_min * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);
		}
		else {
			cton_sap = pft.cton_sap_avr;
			cton_root = pft.cton_root_avr;
		}

		respiration(gtemp, patch.soil.gtemp, indiv.pft.lifeform,
			indiv.pft.respcoeff, cton_sap, cton_root,
			indiv.cmass_sap, cmass_root, 0.0, 0.0, assim, resp);

		// Convert to averages for this period for accounting purposes
		assim /= date.subdaily;
		resp /= date.subdaily;

		// Update accumulated annual NPP and daily vegetation-atmosphere flux
		indiv.dgpp = assim;
		indiv.dnpp = assim - resp;
		indiv.anpp += indiv.dnpp;

		// Save historical GPP for daily growth of grass
		indiv.gpp_15.add(assim);

		indiv.report_flux(Fluxes::NPP, indiv.dnpp);
		indiv.report_flux(Fluxes::GPP, assim);
		indiv.report_flux(Fluxes::RA, resp);

		if (indiv.lai_today() > indiv.mlai_max[date.month])
			indiv.mlai_max[date.month] = indiv.lai_today();

		if (day.isend) {
			indiv.mlai[date.month] += indiv.lai_today() / (double)date.ndaymonth[date.month];
		}

		vegetation.nextobj();
	}
}

/// Leaf senescence for crops Eqs. 8,9,13 and 14 in Olin 2015
void leaf_senescence(Vegetation& vegetation) {

	if (!(vegetation.patch.stand.is_true_crop_stand() && ifnlim)) {
		return;
	}

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// Age dependent N retranslocation, Sec. 2.1.3 Olin 2015
		if (indiv.patchpft().cropphen->dev_stage > 1.0) {
			const double senNr = 0.1;
			double senN = senNr * (indiv.get_nmass_leaf()-indiv.cmass_leaf_today() / (indiv.pft.cton_leaf_max));

			// Senescence is not done during the period without no N-limitation
			if (date.year > freenyears && senN > 0) {
				indiv.nmass_leaf -= senN;
				assert(indiv.nmass_leaf >= 0);
				indiv.cropindiv->nmass_agpool += senN;
			}
		}

		double r = 0.0;

		// N dependent C mass loss, with an inertia of 1/10, Eq. 13 Olin 2015
		if (indiv.cmass_leaf_today() > 0.0) {
			double Ln = indiv.lai_nitrogen_today();
			double Lnld = indiv.lai_today();
			r = (Lnld - min(Lnld, Ln))/indiv.sla()/10.0;
		}

		// No senescence during the initial growing period
		if (indiv.patchpft().cropphen->fphu < 0.05) {
			indiv.daily_cmass_leafloss = 0.0;
		} 
		else {
			indiv.daily_cmass_leafloss = max(0.0, r);
		}
		indiv.daily_nmass_leafloss = 0.0;

		vegetation.nextobj();
	}
}

/// Forest-floor conditions
/** Called in cohort/individual mode (not population mode) to quantify growth
 *  conditions at the forest floor for each PFT
 *  Calculates net assimilation at top of grass canopy (or at soil surface if
 *  there is none).
 */
void forest_floor_conditions(Patch& patch) {

	const Climate& climate = patch.get_climate();
	double lambda;			// not used here
	PhotosynthesisResult phot;

	for (int p=0; p<npft; p++) {

		Patchpft& ppft = patch.pft[p];
		Standpft& spft = patch.stand.pft[p];
		if (!spft.active) {
			continue;
		}
		Pft& pft = spft.pft;

		// peatland limits on photosynthesis 
		double pftco2 = get_co2(patch, patch.stand.get_gridcell().climate, pft);
		double inund_stress = get_inund_stress(patch, ppft);
		double graminoid_wtp_limit = get_graminoid_wtp_limit(patch, pft);
		double moss_wtp_limit = get_moss_wtp_limit(patch, pft);

		// Initialise net photosynthesis sum on first day of year
		if (date.day == 0) {
			ppft.anetps_ff = 0.0;
		}
		if (patch.stand.landcover != CROPLAND || pft.phenology != CROPGREEN && ppft.cropphen->growingseason) {
			double assim = 0;

			if (ppft.wstress_day) {

				//PhotosynthesisStresses ps_stress;
				//ps_stress.set(false, get_moss_wtp_limit(patch, pft), get_graminoid_wtp_limit(patch, pft), get_inund_stress(patch, ppft));

				assimilation_wstress(pft, pftco2, climate.temp, climate.par,
					climate.daylength, patch.fpar_grass * ppft.phen, 1., ppft.gcbase_day,
					spft.photosynthesis.vm, phot, lambda, 1.0, false, 
					moss_wtp_limit, graminoid_wtp_limit, inund_stress);
				assim = phot.net_assimilation();
			}
			else {
				assim = spft.photosynthesis.net_assimilation() * ppft.phen * patch.fpar_grass;
			}

			// Accumulate annual value
			ppft.anetps_ff += assim;
		}

		if (date.islastmonth && date.islastday) {

			// Avoid negative ppft.anetps_ff
			ppft.anetps_ff = max(0.0, ppft.anetps_ff);

			if (ppft.anetps_ff > spft.anetps_ff_max) {
				spft.anetps_ff_max = ppft.anetps_ff;
			}
		}
	}
}

/// Initiate required variables for the module
void init_canexch(Patch& patch, Climate& climate, Vegetation& vegetation) {

	if (date.day == 0) {
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();

			indiv.anpp           = 0.0;

			indiv.leafndemand    = 0.0;
			indiv.rootndemand    = 0.0;
			indiv.sapndemand     = 0.0;
			indiv.storendemand   = 0.0;
			indiv.hondemand      = 0.0;

			indiv.nday_leafon    = 0;
			indiv.avmaxnlim      = 1.0;
			indiv.cton_leaf_aavr = 0.0;

			if (!negligible(indiv.get_cmass_leaf()) && !negligible(indiv.get_nmass_leaf()))
				indiv.cton_leaf_aopt = indiv.get_cmass_leaf() / indiv.get_nmass_leaf();
			else
				indiv.cton_leaf_aopt = indiv.pft.cton_leaf_max;

 			for (int m=0; m<12; m++) {
				indiv.mlai[m] = 0.0;
				indiv.mlai_max[m] = 0.0;
			}

			vegetation.nextobj();
		}
	}

	patch.wdemand_day = 0;
}


// Calculate Vcmax and Jmax from leaf N following Ellsworth et al. 2022
void vcmax_from_N(Vegetation& vegetation) {

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

        double nmass_leaf;
		double leaf_NtoP = 22.9; // from EucFACE -> TODO: move to constants
		
        // leaf N of green leaves per leaf area (kgN m-2 leaf)
		// take green leaves only and calculate leaf N per unit leaf area
		// (not total plant LAI) to not double-account for canopy gradients in leaf N.
		if (indiv.pft.lifeform == GRASS) {
			nmass_leaf = (!negligible(indiv.nmass_leafs[0])) ? 
						indiv.nmass_leafs[0] / (indiv.cmass_leafs[0] * indiv.sla()) : 0.0;
		} else {
			nmass_leaf = indiv.nmass_leaf;
		}

		assert(nmass_leaf >= 0);

		// leaf P in first fully expanded leaves per unit leaf area (kgP m-2 leaf) 
		double pmass_leaf = nmass_leaf / leaf_NtoP;         

		if (!negligible(pmass_leaf, -6)) {
			indiv.ps_subdaily.vcmax25 = exp(4.308 + 0.298 * log(pmass_leaf * G_PER_KG) + 0.197 * 
											log(nmass_leaf * G_PER_KG)) / UMOL_PER_MOL;
			indiv.ps_subdaily.jmax25 = exp(5.139 + 0.325 * log(pmass_leaf * G_PER_KG) + 0.112 * 
											log(nmass_leaf * G_PER_KG)) / UMOL_PER_MOL;
		}

		vegetation.nextobj();
	}

}


// Coupled model of photosynthesis and stomatal conductance for C3 and C4 vegetation.
void photosynthesis_gs(Patch& patch, Climate& climate, Vegetation& vegetation) {

    double par = 0.0, temp = 0.0, rad = 0.0, vpd = 0.0, ca = 0.0;
	double lai_layer = 0.0, apar_layer = 0.0, lai_cum_bottom = 0.0, lai_cum = 0.0, Nscal = 0.0, pslim = 0.0;
	double g0 = 0.0, g1 = 0.0, X = 0.0, beta = 0.0, gamma = 0.0, gamstar = 0.0, Ko = 0.0, Kc = 0.0;
	double Rd = 0.0, vcmax = 0.0, jmax = 0.0;
	double Rm =0.0, Rs = 0.0, Je = 0.0, Vpmax = 0.0, Vp = 0.0, Kp = 0.0, Km = 0.0;
	double Anc = 0.0, Anj = 0.0, Anp = 0.0, Ag = 0.0, An = 0.0, gsc = 0.0, ci = 0.0;
	double ci_diff = 0.0, ci_diff_max = 0.0, ci_new = 0.0;
	double a_c = 0.0, b_c = 0.0, c_c = 0.0, a_j = 0.0, b_j = 0.0, c_j = 0.0; 
	int    nit = 0.0, nit_max = 0.0;
	
	double height_veg = get_veg_height(vegetation);

	// Loop through cohorts (=individuals)
	vegetation.firstobj();
	while (vegetation.isobj) {

		// Retrieve information about cohort
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		// get PFT-specfic parameters
		g0 = pft.g0;  // mol m-2 s-1
		g1 = pft.g1;  // kPa^0.5

		// Atmospheric CO2 concentration (mol mol-1)
		// TODO: check if annual global means would make more sense 
 		ca = get_co2(patch, climate, pft) / UMOL_PER_MOL;  // mol mol-1

		indiv.photosynthesis.agd_g = 0.0;
		indiv.photosynthesis.rd_g = 0.0;
		indiv.photosynthesis.je = 0.0;

		// Loop over time steps per day
        for (Day day; day.period != date.subdaily; day.next()) {

			// initialise timestep-level variables
			double Rd_canopy = 0.0;
			double An_canopy = 0.0;
			double Ag_canopy = 0.0;
			double Anc_canopy = 0.0;
			double Anj_canopy = 0.0;
			double gsw_canopy = 0.0;
			double ci_canopy = 0.0;
			double vcmax_canopy = 0.0;
			double jmax_canopy = 0.0;

			// get subdiurnal climate variables
			par   = climate.pars[day.period];    // W m-2; 0.5 * climate.insols
			temp  = climate.temps[day.period];   // degC
			vpd   = climate.vpds[day.period];    // kPa	

			// Stomatal term (Medlyn et al. 2011)
			// X needs to represent gs for CO2, hence no 1.6 factor
    		X = divide(divide(1.0 + ppft.wscal * g1, sqrt(max(vpd,0.1))), ca); 

			// initialise cumulative LAI
			lai_cum_bottom = 0.0;

			// Loop over canopy layers (from top layer to bottom layer)
		    for (int i = 0; i < indiv.fpar_layerwise.size(); i++) {

				// determine fpar and apar of the layer
		        double fpar_layer = indiv.fpar_layerwise[i]; // fPAR in this layer
				double lai_layer = indiv.lai_layerwise[i];   // LAI of the layer and individual

				if (negligible(fpar_layer, -12)) {
					if (negligible(lai_layer, -12)) {
						// LAI and fpar both 0 - no photosynthesis calculation 
						// tree may be shorter than this layer.
						continue;
					} else if (!negligible(par)) {
						// fpar is 0, lai and par both non zero. This means there
						// is par and leaf in this layer, but no fpar.
						// For grasses, this can happen if PAR is below the threshold
						// required for grass growth.
						// if (pft.lifeform != GRASS && pft.lifeform != MOSS) {
						// 	fail("Lai (%.2f) and par (%.2f) are both nonzero but fpar is 0 in layer %d day %d year %d patch %d",
						// 	lai_layer, par, i, date.day, date.year, patch.id);
						// }
					}
				}

				// determine LAI and cumulative LAI of the layer
				lai_cum_bottom += lai_layer;                 	// cumulative LAI at the bottom of current layer
				lai_cum = lai_cum_bottom - 0.5 * lai_layer;  	// cumulative LAI at the middle of current layer
				apar_layer = lambertbeer(lai_cum) * par;	    // available PAR in the middle of the layer (W m-2), TODO: check for double counting of absorption

				// determine N scalar of the layer
				//Nscal = (i == 0) ? 1.0 : max(1.0e-12, exp(-kn * lai_cum) / exp(-kn * 0.5 * indiv.lai_layerwise[0]));
				Nscal = max(1.0e-12, exp(-kn * lai_cum));
				// TODO discretise grass layers

          		if (pft.pathway == C3) 	{	// C3 photosynthesis (largely Farquhar et al. 1980)
					
					// Daytime (mitochondrial) leaf respiration (mol m-2 s-1)
					Rd = indiv.fphen_activity * Nscal * ppft.wscal * BC3 * indiv.ps_subdaily.vcmax25 * Arrhenius_temp_response(1.0,temp,HA_J,0.0,0.0);  // TODO: check temperature response

					if (negligible(apar_layer, -9)) {
						// Night
						An = -Rd;
						gsc = g0;
					} else {

						// Michaelis-Menten parameters and photorespiratory CO2 compensation point (mol mol-1)
						Ko      = Arrhenius_temp_response(KO25,temp,HA_KO,0.0,0.0);
						Kc      = Arrhenius_temp_response(KC25,temp,HA_KC,0.0,0.0); 
						gamstar = Arrhenius_temp_response(GAMSTAR25,temp,HA_GS,0.0,0.0);

						// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax
						pslim = indiv.fphen_activity * indiv.ps_sink_lim * Nscal * ppft.wscal;
						vcmax = pslim * Arrhenius_temp_response(indiv.ps_subdaily.vcmax25,temp,HA_V,HD,DS_V);
						jmax  = pslim * Arrhenius_temp_response(indiv.ps_subdaily.jmax25,temp,HA_J,HD,DS_J);

						// 1) Rubisco-limited gross photosynthesis
						gamma = vcmax;
						beta  = Kc * (1.0 + divide(O2atm,Ko));
						Anc   = photosyn_solver(g0,g1,X,gamma,beta,gamstar,Rd,ca);

						// 2) RuBP regeneration-limited gross photosynthesis
						Je    = electron_transport_rate(apar_layer,jmax);   
						gamma = 0.25 * Je;    // assuming demand of 4 electrons per CO2 molecule
						beta  = 2.0 * gamstar;
						Anj   = photosyn_solver(g0,g1,X,gamma,beta,gamstar,Rd,ca);

						// 3) TPU-limited gross photosynthesis
						// Ap = 3TPU; TPU = 0.167Vcmax => Ap = 0.5 Vcmax
						Anp = 0.5 * vcmax - Rd;

						// Net photosynthesis as minimum of 1), 2), and 3) (mol m-2 s-1)
						An = min(Anc,min(Anj,Anp));

						// Stomatal conductance to CO2 (mol m-2 s-1)
						gsc = g0 + An * X;
						if (An < 0) {
							gsc = g0;
						}
					}

					// Re-calculate intercellular CO2 concentration (mol mol-1)
					ci = ca - divide(An, gsc);

				}
				else {  // C4 photosynthesis (von Caemmerer 2000)

					//  Daytime (mitochondrial) leaf respiration (mol m-2 s-1)
					// TODO: check BC4 (half the value in von Caemmerer 2000) and Vcmax
					Rd = indiv.fphen_activity * Nscal * ppft.wscal * BC4 * Vcmax25_C4 * Arrhenius_temp_response(1.0,temp,HA_RD_C4,0.0,0.0);
					
					Rm = 0.5 * Rd; // Rd in the mesophyll cells
					Rs = Rd - Rm;  // Rd in the bundle-sheath cells

					if (negligible(par, -9)) {
						// No photosynthesis at night.
						An = -Rd;
						gsc = g0;
						ci = ca - divide(An, gsc);
					} else {
						// C4/day time

						// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax
						// TODO: different temp dependency in C3 and C4 plants??
						//vcmax = Nscal * ppft.wscal * Arrhenius_temp_response(Vcmax25_C4,temp,HA_V4,HD_V4,DS_V4);
						pslim = indiv.fphen_activity * indiv.ps_sink_lim * Nscal * ppft.wscal;
						vcmax = pslim * Arrhenius_temp_response(Vcmax25_C4,temp,HA_V4,0.0,0.0);
						//jmax  = Nscal * ppft.wscal * Arrhenius_temp_response(Jmax25_C4,temp,HA_J4,HD_J4,DS_J4);
						jmax  = pslim * (Jmax25_C4 * exp(-(pow((temp - 43.0)/26.0,2))));  // TODO: replace with equivalent Arrhenius equation.

						// Temperature dependency of Vpmax and Michaelis-Menten constants
						Vpmax = pslim * Arrhenius_temp_response(Vpmax25,temp,HA_VP,HD_VP,DS_VP);
						Kp    = Arrhenius_temp_response(KP25,temp,HA_KP,0.0,0.0);
						Kc    = Arrhenius_temp_response(KC25_C4,temp,HA_KC_C4,0.0,0.0);
						Ko    = Arrhenius_temp_response(KO25_C4,temp,HA_KO_C4,0.0,0.0);
						Km    = Kc * (1.0 + divide(O2atm,Ko));

						// Iteration to find An, gs, and ci
						ci = 0.5 * ca;
						ci_diff_max = 0.5e-6; // mol mol-1
						ci_diff = 2.0 * ci_diff_max;
						nit = 0;
						nit_max = 20;
						while (abs(ci_diff) > ci_diff_max && nit < nit_max) {
							if (nit == nit_max - 1) {
								// No solution found. Use initial guess.
								ci = 0.5 * ca;
							}
							// PEP carboxylation rate
							Vp = min(divide(ci * Vpmax, ci + Kp),Vpr);

							// Enzyme-limited net photosynthesis rate (mol m-2 s-1)
							if (ifapproximate_C4){
								Anc = min((Vp + gbs * ci - Rm),(vcmax - Rd));
							} 
							else {
								a_c = 1.0 - alpha_PSII/a0 * Kc/Ko;
								b_c = -((Vp - Rm + gbs * ci) + (vcmax - Rd) + gbs*Km + alpha_PSII/a0 
								 	  * (lgamstar * vcmax + Rd * Kc/Ko));
								c_c = (vcmax - Rd) * (Vp - Rm + gbs * ci) - 
								  	  (vcmax * gbs * lgamstar * O2atm + Rd * gbs * Km);

								Anc = (-b_c - sqrt(pow(b_c,2) - 4.0 * a_c * c_c)) / (2 * a_c);
							}
							// TODO: check value of alpha_PSII.
							// TODO: check differences to von Caemmerer 2021.

							// Light-limited net photosynthesis rate (mol m-2 s-1)
							// electron transport rate (same as C3) (mol (electrons) m-2 s-1)
							Je = electron_transport_rate(apar_layer,jmax);
							if (ifapproximate_C4){
								Anj = min(((xj * Je)/2.0 - Rm + gbs * ci),(((1.0 - xj) * Je)/3.0 - Rd));
							} 
							else {
								a_j = 1.0 - (7.0 * lgamstar * alpha_PSII) / (3.0 * a0);
								b_j = -(((xj * Je)/2.0 - Rm + gbs * ci) + (((1.0 - xj) * Je)/3.0 - Rd) 
								  	  + gbs * ((7.0 * lgamstar * O2atm)/3.0) + 
								  	  (alpha_PSII * lgamstar)/a0 * (((1.0 - xj) * Je)/3.0 + (7.0 * Rd)/3.0));
								c_j = (((xj * Je)/2.0 - Rm + gbs * ci) * (((1.0 - xj) * Je)/3.0 - Rd)) - 
								  	  gbs * lgamstar * O2atm * (((1.0 - xj) * Je)/3.0 + (7.0 * Rd)/3.0);

								Anj = (-b_j - sqrt(pow(b_j,2) - 4.0 * a_j * c_j)) / (2 * a_j);
							}
							// TODO: check xj vs. z in von Caemmerer 2021

							// Net photosynthesis rate (mol m-2 s-1)
							An = min(Anc,Anj);

							// Stomatal conductance to CO2 (mol m-2 s-1)
							gsc = g0 + An * X;

							if (An < 0) {
								gsc = g0;
							}

							// Intercellular CO2 concentration (mol mol-1)
							ci_new = ca - divide(An, gsc);

							ci_diff = ci - ci_new;
							if (nit != nit_max - 1) {
								// If we can't find a solution, use initial guess
								// on last iteration.
								ci -= ci_diff / 2.0;
							}
							nit += 1;
						} // end solver
						if (nit == nit_max) {
							// If we didn't find a solution, use initial guess.
							ci = 0.5 * ca;
							gsc = max(gsc, 0.0);
						}
					} // end if c4+daytime

				}  // end photosynthetic pathway switch

				// sum variables over canopy layers
				Rd_canopy  += Rd * lai_layer;
				An_canopy  += An * lai_layer;
				Ag_canopy  += (An + Rd) * lai_layer;
				Anc_canopy += Anc * lai_layer;
				Anj_canopy += Anj * lai_layer;
				gsw_canopy += (gsc * DWDC) * lai_layer;
				ci_canopy  += ci * lai_layer;
				vcmax_canopy += vcmax * lai_layer;
				jmax_canopy  += jmax * lai_layer;
				
			} // end canopy layer loop

			// calculate mean for some variables
			ci_canopy /= indiv.lai;
			vcmax_canopy /= indiv.lai;
			jmax_canopy /= indiv.lai;

			// save key variables to PhotosynthesisSubdaily Structure
			indiv.ps_subdaily.Rd[day.period]    = Rd_canopy;
			indiv.ps_subdaily.An[day.period]    = An_canopy;
			indiv.ps_subdaily.Anc[day.period]   = Anc_canopy;
			indiv.ps_subdaily.Anj[day.period]   = Anj_canopy;
			indiv.ps_subdaily.gsw[day.period]   = gsw_canopy;
			indiv.ps_subdaily.ci[day.period]    = ci_canopy;
			indiv.ps_subdaily.vcmax[day.period] = vcmax_canopy;
			indiv.ps_subdaily.jmax[day.period]  = jmax_canopy;

			// sum up over the day
			indiv.photosynthesis.agd_g += Ag_canopy;
			indiv.photosynthesis.rd_g  += Rd_canopy;
			indiv.photosynthesis.je    += indiv.ps_subdaily.Anj[day.period];

		} // end subdaily loop

		// Calculate daily means for the PhotosynthesisResult Structure that goes to output
		indiv.photosynthesis.agd_g /= date.subdaily;
		indiv.photosynthesis.rd_g  /= date.subdaily;
		indiv.photosynthesis.je    /= date.subdaily;

		// Convert units from molar to mass units and from seconds to day
		indiv.photosynthesis.agd_g *= MOL2G_C * SECS_PER_DAY;  // gC m-2 d-1
		indiv.photosynthesis.rd_g *= MOL2G_C * SECS_PER_DAY;   // gC m-2 d-1
		// Give values all variables even they aren't used (yet)
		indiv.photosynthesis.vm          = indiv.ps_subdaily.vcmax25 * UMOL_PER_MOL;  // umol m-2 s-1; not used
		indiv.photosynthesis.adtmm       = 0.0;    // not used
		indiv.photosynthesis.nactive_opt = 0.0;    // not used
		indiv.photosynthesis.vmaxnlim    = 0.0;    // not used 

		vegetation.nextobj();

	} // end cohort loop
}


// System of equations to solve coupled photosynthesis-stomatal conductance scheme.
// Wang & Leuning 1998
double photosyn_solver(double g0, double g1, double X, double gamma, double beta, 
					   double gamstar, double Rd, double Ca) {
	
	double b2, b1, b0, Ci, An;

	b2 = g0 + X * (gamma - Rd);
	b1 = (1.0 - Ca * X) * (gamma - Rd) + g0 * (beta - Ca) - X * (gamma * gamstar + beta * Rd);
	b0 = -(1.0 - Ca * X) * (gamma * gamstar + beta * Rd) - g0 * beta * Ca;

	assert(b2 != 0);
	Ci = (-b1 + sqrt(pow(b1,2) - 4.0 * b0 * b2)) / (2.0 * b2); 
	
	An = gamma * divide(Ci - gamstar, Ci + beta) - Rd;

	return(An);  // net photosynthesis (mol m-2 s-1)

}

// Temperature response function following a modified Arrhenius equation (Johnson et al. 1942) 
double Arrhenius_temp_response(double k25, double Tl, double Ha, double Hd, double dS) {

	/// Inputs
	// k25 = process rate at 25 degC
	// Tl  = leaf temperature (deg C)
	// Ha  = activation energy (kJ mol-1)
	// Hd  = deactivation energy (kJ mol-1)
	// dS  = entropy term (J mol-1 K-1)
	
	// Define reference temperature and convert units
	double Tref = 25.0 + K2degC; // Reference Temperature (25 degC in K)
    Tl += K2degC;
    Ha *= J_PER_KJ; 
	Hd *= J_PER_KJ; 

	double kTl = k25 * exp(Ha * (Tl - Tref)/(Tref * R_gas * Tl)) * 
				      (1.0 + exp((Tref * dS - Hd)/(Tref * R_gas))) / 
				      (1.0 + exp((Tl * dS - Hd)/(Tl * R_gas)));

	return(kTl);  // process rate at leaf Temperature Tl
}


// Calculate electron transport rate from PAR
// e.g. von Caemmerer 2013
double electron_transport_rate(double PAR,double Jmax) {

	/// Inputs
	// PAR  = incident photosynthetically active radiation (W m-2) 
    // Jmax = max. electron transport rate (mol m-2 s-1)

	double PPFD = PAR * JTOMOL;  // photosynthetic photon flux density (mol (photons) m-2 s-1) 
	double alpha = 0.3; // quantum yield of electron transport (mol (electrons) mol-1 (photon))
						// value from Medlyn et al. 2002
						// von Caemmerer 2013: ~0.36
	double J = (alpha * PPFD + Jmax - sqrt(pow(alpha * PPFD + Jmax,2) - 
				4.0 * THETA * alpha * PPFD * Jmax)) / (2.0 * THETA);  

	return(J);   // electron transport rate in mol (electrons) m-2 s-1

}

/*
DAVE water scalar routine (see water_scalar()).

@param patch: The patch.
*/
void water_stress(Patch &patch) {
	// Must set:
	// - ppft.wscal_mean (scalar, total annual water stress, convert to mean at end of year)
	// - ppft.cropphen->growingdays_y -> increment once per day for true crops and cover-crop grass
	// - ppft.wscal -> mean daily wscal value in range [0, day_len]
	// - indiv.wscal_running (call .add() to add today's value in once at end of day)
	// - ppft.wscal_5 and ppft.wscal_10 -> 5 and 10 day running mean of wscal

	// Iterate over individuals in patch.
	patch.vegetation.firstobj();
	while (patch.vegetation.isobj) {
		Individual &indiv = patch.vegetation.getobj();
		Patchpft &ppft = patch.pft[indiv.pft.id];
		Pft &pft = ppft.pft; 
		
		double wstress_layer[NSOILLAYER];
		double wstress = 0.0;
		double rootdist_sum = 0.0;

		for (int i = 0; i<NSOILLAYER; i++) {
			// get actual water content and available water content per layer
			double wcont = patch.soil.get_layer_soil_water(i);

			// calculate water stress factor per layer and weight with root fraction
			wstress_layer[i] = max(min(wcont / wcont_crit,1.0),0.0) * pft.rootdist[i];
			wstress += wstress_layer[i];
			rootdist_sum += pft.rootdist[i]; // rootdist should add up to 1 but summed here just in case...
		}
		// water stress factor [0-1] in the root zone
		wstress /= rootdist_sum;

		// Mean daily water stress value.
		ppft.wscal = wstress;

		if (date.day == 0 && (ppft.pft.phenology==CROPGREEN || ppft.pft.isintercropgrass) ) {
			ppft.cropphen->growingdays_y = 0;
		}

		if (patch.stand.landcover!=CROPLAND										//natural, urban, pasture, forest and peatland stands
			|| ppft.pft.phenology==ANY &&
			ppft.pft.id==patch.stand.pftid) {
			// Normal grass growth.
			ppft.wscal_mean += ppft.wscal;

			// Convert from sum to mean on last day of year.
			if (date.islastday && date.islastmonth) {
				ppft.wscal_mean /= date.year_length();
			}
		}
		else if (ppft.cropphen->growingseason
				|| ppft.pft.phenology == CROPGREEN && date.day == ppft.cropphen->hdate
				|| ppft.pft.isintercropgrass && date.day == patch.pft[patch.stand.pftid].cropphen->eicdate) {
			// True crops and cover-crop grass.
			ppft.cropphen->growingdays_y++;
			ppft.wscal_mean = ppft.wscal_mean + (ppft.wscal - ppft.wscal_mean) / ppft.cropphen->growingdays_y;
		} else {
			// No wscal_mean for you...
		}

		ppft.wscal_5.add(ppft.wscal);
		ppft.wscal_10.add(ppft.wscal);

		patch.vegetation.nextobj();
	}
}


// Canopy transpiration
void transpiration(Patch &patch, Climate& climate, Vegetation& vegetation) {

	// this formulation assumes infinite Ga
	// TODO: include Ga --> Penman-Monteith (?)

	double wsupply = 0;
	double wcont_local[NSOILLAYER];
	for (int i = 0; i < NSOILLAYER; i++) {
		wcont_local[i] = patch.soil.get_layer_soil_water(i);
		wsupply += wcont_local[i];
	}

	// Loop through cohorts (=individuals)
	vegetation.firstobj();
	while (vegetation.isobj) {

		// Retrieve information about cohort
		Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (date.day == 0)
			indiv.aaet = 0.0;

		double transp = 0;

		// Loop over time steps per day
		for (Day day; day.period != date.subdaily; day.next()) {	

			// T = Gc / (pressure/VPD) 
			// Gc in mol m-2 s-1
			// VPD and pressure in Pa
			double patmos = climate.pressure[day.period];
			double vpd = climate.vpds[day.period];
			double gsw = indiv.ps_subdaily.gsw[day.period];
			if (vpd > 0 && patmos > 0) {
				transp += gsw / (patmos / (vpd * 1000.0));
			}

		} // end sub-daily loop
		
		// daily transpiration rate in mol m-2 s-1 
		transp /= date.subdaily;

		// daily transpiration rate in mm d-1
		indiv.aet = min(wsupply, transp * MOL2G_H2O * KG_PER_G * SECS_PER_DAY); // assuming 1kg m-2 = 1l m-2 == 1mm
        
		// add to annual transpiration
		indiv.aaet += indiv.aet;

		// Note that in standard LPJ-GUESS transpiration is described as evapotranspiration (AET),
		// which is incorrect as soil evaporation is calculated separately.

		// Limit water uptake.
		double wr;

		Standpft& spft = patch.stand.pft[indiv.pft.id];

		if (spft.irrigated && pft.id == patch.stand.pftid) {
			Day day;
			wr = irrigated_water_uptake(patch, pft, day);
		}
		else {
			double wcont_local[NSOILLAYER];
			for (int i = 0; i<NSOILLAYER; i++) {
				wcont_local[i] = patch.soil.get_layer_soil_water(i);
			}

			if (iftwolayersoil) {
				wr = water_uptake_twolayer(wcont_local, patch.soil.soiltype.awc,
							pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
							pft.lifeform == TREE, pft.drought_tolerance);
			} else {
				wr = water_uptake(wcont_local, patch.soil.soiltype.awc,
							pft.rootdist, pft.emax, patch.fpc_rescale, ppft.fwuptake,
							pft.lifeform == TREE, pft.drought_tolerance);
			}
		}

		vegetation.nextobj();
	} // end cohort loop

	// Ensure water uptake doesn't exceed availability in any layer.
	double uptake[NSOILLAYER];
	double total_uptake = 0;
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		const Patchpft &ppft = indiv.patchpft();
		for (int j = 0; j < NSOILLAYER; j++) {
			uptake[j] += indiv.aet * ppft.fwuptake[j];
			total_uptake += uptake[j];
		}
	}
	for (int i = 0; i < NSOILLAYER; i++) {
		if (uptake[i] > wcont_local[i]) {
			// Demand in this layer exceeds supply.
			double fwavail = wcont_local[i] / uptake[i];
			for (int j = 0; j < vegetation.nobj; j++) {
				vegetation[j].patchpft().fwuptake[i] *= fwavail;
			}
		}
	}
}



// This is a modification of the original npp routine of standard LPJ-GUESS.
// Modifications include:
// - called once a day (no sub-daily option, day structure not an input)
// - no more calculation of photosynthesis
// - bvoc function not implemented (could be added later)
void npp_dave(Patch& patch, Climate& climate, Vegetation& vegetation) {

	double assim, resp, lambda, gtemp;

	gtemp = climate.gtemp;  // TODO: find out what that is.

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// For this individual ...

		// Retrieve PFT and patch PFT

		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[pft.id];

		//Don't do calculations for crops outside their growingseason
		if (!indiv.growingseason()) {
			indiv.dnpp=0.0;
			indiv.lambda = 0;
			vegetation.nextobj();
			continue;
		}

		// get daily net assimilation
		assim = indiv.photosynthesis.net_assimilation();
		
		// Calculate autotrophic respiration
		double cmass_root;
		if (indiv.cropindiv && indiv.cropindiv->isintercropgrass && indiv.phen == 0.0)
			cmass_root = 0.0;
		else
			cmass_root = indiv.cmass_root_today();

		// Static root and sap wood C:N ratio if no N limitation
		// to not let N affect respiration for C only version of model
		// Also not allow lower C:N ratios than minimum which
		// would affect respiration severely (happens for deciduous trees
		// in the beginning of the season)

		double cton_sap, cton_root;
		if (ifnlim) {
			cton_sap = max(indiv.cton_sap(), indiv.pft.cton_leaf_min * indiv.pft.cton_sap_avr / indiv.pft.cton_leaf_avr);
			cton_root = max(indiv.cton_root(), indiv.pft.cton_leaf_min * indiv.pft.cton_root_avr / indiv.pft.cton_leaf_avr);
		}
		else {
			cton_sap = pft.cton_sap_avr;
			cton_root = pft.cton_root_avr;
		}

		respiration(gtemp, patch.soil.gtemp, indiv.pft.lifeform,
			indiv.pft.respcoeff, cton_sap, cton_root, indiv.cmass_sap, 
			cmass_root, indiv.growth_demand, indiv.fphen_activity, assim, resp);


		// Update accumulated annual NPP and daily vegetation-atmosphere flux
		// Note that resp does not include leaf day respiration.
		indiv.dgpp = assim;
		indiv.dnpp = assim - resp;
		indiv.anpp += indiv.dnpp;

		// Save historical GPP for daily growth of grass
		indiv.gpp_15.add(assim);

		indiv.report_flux(Fluxes::NPP, indiv.dnpp);
		indiv.report_flux(Fluxes::GPP, assim);
		indiv.report_flux(Fluxes::RA, resp);

		if (indiv.lai_today() > indiv.mlai_max[date.month])
			indiv.mlai_max[date.month] = indiv.lai_today();

		//if (day.isend) {
		indiv.mlai[date.month] += indiv.lai_today() / (double)date.ndaymonth[date.month];
		//}

		vegetation.nextobj();
	}
}

/// Canopy exchange
/** Vegetation-atmosphere exchange of CO2 and water including calculations
 *  of actual evapotranspiration (AET), canopy conductance, carbon assimilation
 *  and autotrophic respiration.
 *  Should be called each simulation day for each modelled area or patch,
 *  following update of leaf phenology and soil temperature and prior to update
 *  of soil water.
 */
void canopy_exchange(Patch& patch, Climate& climate) {

	// NEW ASSUMPTIONS CONCERNING FPC AND FPAR (Ben Smith 2002-02-20)
	// FPAR = average individual fraction of PAR absorbed on patch basis today,
	//        including effect of current leaf phenology (this differs from previous
	//        versions of LPJ-GUESS in which FPAR was on an FPC basis)
	// FPC =  PFT population (population mode), cohort (cohort mode) or individual
	//        (individual mode) fractional projective cover as a fraction of patch area
	//        (in population mode, corresponds to LPJF variable fpc_grid). Updated
	//        annually based on leaf-out LAI (see function allometry in growth module).
	//        (FPC was previously equal to summed crown area as a fraction of patch
	//        area in cohort/individual mode)

	// Retrieve Vegetation and Climate objects for this patch
	Vegetation& vegetation = patch.vegetation;

	// Initial no-stress canopy exchange processes
	init_canexch(patch, climate, vegetation);

	// Canopy exchange processes
	fpar(patch);


    if (ifdavephotosyn) {

		// Supply function for nitrogen and determination of nitrogen stress leading
		// to down-regulation of vmax.

		/*
		// Nitrogen within projective cover of all individuals
		double tot_nmass_avail = patch.soil.nmass_avail * min(1.0, patch.fpc_total);

		if (patch.stand.landcover == CROPLAND && ifnlim) { // Also for other landcovers ??
			// Take soil wcont into account
			tot_nmass_avail *= patch.soil.get_soil_water_upper() * 0.9 + patch.soil.get_soil_water_lower() * 0.1;
		}

		// Nitrogen demand
		ndemand(patch, vegetation);

		// Calculate individual uptake fraction of nitrogen demand
		if (patch.ndemand > tot_nmass_avail) {

			// Determine individual nitrogen uptake fractions
			fnuptake(vegetation, tot_nmass_avail);
		}

		// Resolve nitrogen stress with longterm stored nitrogen
		nstore_usage(vegetation);

		*/

		// Calculate photosynthetic capacity from leaf N
		vcmax_from_N(vegetation);

		// Calculate water stress factor
		water_stress(patch);

		// Calculate subdaily photosynthesis and stomatal conductance 
		photosynthesis_gs(patch, climate, vegetation);

		// Transpiration from VPD and canopy conductance
		transpiration(patch, climate, vegetation);

		// Calculate respiration and NPP on a daily basis
		npp_dave(patch, climate, vegetation);

		leaf_senescence(vegetation);

	}
	else {
		// Calculates no-stress daily values of photosynthesis and gpterm
		photosynthesis_nostress(patch, climate);

		// Nitrogen demand
		ndemand(patch, vegetation);

		// Nitrogen stress
		vmax_nitrogen_stress(patch, climate, vegetation);

		// Only these processes are affected in diurnal mode
		for (Day day; day.period != date.subdaily; day.next()) {

			wdemand(patch, climate, vegetation, day);
			aet_water_stress(patch, vegetation, day);
			water_scalar(patch, vegetation, day);
			npp(patch, climate, vegetation, day);
		}

		leaf_senescence(vegetation);

		// Forest-floor conditions
		forest_floor_conditions(patch);
	}

	// Total potential evapotranspiration for patch (mm, patch basis)
	// is a sum of: (1) potential transpirative demand of the vegetation;
	// (2) evaporation of canopy-intercepted precipitation; and
	// (3) evaporation from the soil surface of non-vegetated parts of patch -
	//     currently with boleht at 0, a patch has only two surfaces - vegetated
	//     and non-vegetated.
	// This value is only diagnostic, it is not to be used in further calculations.
	// Correct value should use daily value of patch.demand_leafon.
	double pet_patch = patch.wdemand_day * patch.fpc_total + patch.intercep +
				climate.eet * PRIESTLEY_TAYLOR * max(1.0-patch.fpc_total, 0.0);

	patch.apet += pet_patch;
	patch.mpet[date.month] += pet_patch;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2001
// Collatz, GJ, Ball, JT, Grivet C & Berry, JA 1991 Physiological and
//   environmental regulation of stomatal conductance, photosynthesis and
//   transpiration: a model that includes a laminar boundary layer. Agricultural
//   and Forest Meteorology 54: 107-136
// Collatz, GJ, Ribas-Carbo, M & Berry, JA 1992 Coupled photosynthesis-stomatal
//   conductance models for leaves of C4 plants. Australian Journal of Plant
//   Physiology 19: 519-538
// Comins, H. N. & McMurtrie, R. E. 1993. Long-Term Response of Nutrient-Limited
//   Forests to CO2 Enrichment - Equilibrium Behavior of Plant-Soil Models.
//   Ecological Applications, 3, 666-681.
// Farquhar GD & von Caemmerer 1982 Modelling of photosynthetic response to
//   environmental conditions. In: Lange, OL, Nobel PS, Osmond CB, Ziegler H
//   (eds) Physiological Plant Ecology II: Water Relations and Carbon
//   Assimilation, Vol 12B. Springer, Berlin, pp 549-587.
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A
//   process-based, terrestrial biosphere model of ecosystem dynamics
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Gerten, D., Schaphoff, S., Haberlandt, U., Lucht, W. & Sitch, S. 2004.
//   Terrestrial vegetation and water balance - hydrological evaluation of a
//   dynamic global vegetation model. Journal of Hydrology 286: 249-270.
// Haxeltine A & Prentice IC 1996a BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Haxeltine A & Prentice IC 1996b A general model for the light-use efficiency
//   of primary production. Functional Ecology 10: 551-561
// Huntingford, C & Monteith, JL 1998. The behaviour of a mixed-layer model of the
//   convective boundary layer coupled to a big leaf model of surface energy
//   partitioning. Boundary Layer Meteorology 88: 87-101
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Monteith, JL, 1995. Accomodation between transpiring vegetation and the convective
//   boundary layer. Journal of Hydrology 166: 251-263.
// S. Olin, G. Schurgers, M. Lindeskog, D. Wårlind, B. Smith, P. Bodin, J. Holmér, and A. Arneth. 2015
//   Biogeosciences Discuss., 12, 1047-1111. The impact of atmospheric CO2 and N management on yields
//   and tissue C:N in the main wheat regions of Western Europe
// Peltoniemi, MS, Duursma, RA & Medlyn, BE. 2012. Co-optimal distribution of leaf
//   nitrogen and hydraulic conductance in plant canopies. Tree Physiology, 32, 510-519.
// Prentice, IC, Sykes, MT & Cramer W (1993) A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Press, WH, Teukolsky, SA, Vetterling, WT & Flannery, BT. 1986. Numerical
//   Recipes in FORTRAN, 2nd ed. Cambridge University Press, Cambridge
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Sprugel, DG, Ryan MG, Renee Brooks, J, Vogt, KA & Martin, TA (1996) Respiration
//   from the organ level to the stand. In: Smith, WK & Hinckley, TM (eds),
//   Physiological Ecology of Coniferous Forests.
// Wania, R., Ross, I., & Prentice, I.C. (2009a) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: I. Evaluation and sensitivity of physical 
//   land surface processes. Global Biogeochemical Cycles, 23, GB3014, doi:10.1029/2008GB003412
// Wania, R., Ross, I., & Prentice, I.C. (2009b) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: II. Evaluation and sensitivity of vegetation 
//   and carbon cycle processes. Global Biogeochemical Cycles, 23, GB015, doi:10.1029/2008GB003413
// Zaehle, S. & Friend, A. D. 2010. Carbon and nitrogen cycle dynamics in the O-CN
//   land surface model: 1. Model description, site-scale evaluation, and sensitivity
//   to parameter estimates. Global Biogeochemical Cycles, 24.

