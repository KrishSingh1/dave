///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Vegetation-atmosphere exchange of H2O and CO2 via
//                        production, respiration and evapotranspiration
//                        *************************************************************
//                        "Fast" version, revised December 2002 by Ben Smith
//                        Modified according to code changes by Dieter Gerten 021216
//                        Includes updated FPAR formulation in cohort/individual mode
//                        based on changes suggested by Soenke Zaehle
//                        *************************************************************
//                        Corrected error in forest_floor_conditions(): wstress sums
//                        in Patchpft were not converted to means, 2005-01-25
//                        * 2005-03-01: Corrected problem caused by optimisation in pgCC
//                        that resulted in zero FPC for grasses (see function fpar)
// Header file name:      canexch.h
// Source code file name: canexch.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16/2005-01-25
// Updated:               2010-11-22


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


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// GUESSN:
const double N0=7.15*0.001; //7.15*2.0*0.001; // leaf N (kgN/kgC) not associated with photosynthesis
	// (value given by Haxeltine & Prentice 1996, multiplied by 2 to convert from
	// leaf dry weight to leaf carbon basis)
	// BLARP! unclear from Haxeltine & prentice whether the the value No=7.15 is
	// mgN/g (as stated in the paper) or mgN/gC. The value is the intercept of a figure
	// in Field & Mooney 1986, but this intercept appear to be c. 0.25 mmolN/g which
	// equates to 3.5 mgN/g or 7 mgN/gC. If the latter interpretation is correct,
	// the 2.0 should be removed from the above term.
// end GUESSN

///////////////////////////////////////////////////////////////////////////////////////
// PROCESS SWITCHES
// This module contains alternative formulations for several processes. Each of a set
// of alternative formulations is represented by a keyword declared by a #define
// directive (below). All but one of the #define's relating to a particular process
// should be commented out prior to a build (this causes the preprocesser to exclude
// the unwanted code before sending the remaining code on to the compiler)

// Alternative formulations for transpirational demand
//   DEMAND_PATCH = a single demand calculated for the entire patch / grid cell
//   DEMAND_INDIV = a separate demand calculated for each individual (as in LPJF)
// Comment out one of the following two lines:

#define DEMAND_PATCH
//#define DEMAND_INDIV

// Check:
#if defined(DEMAND_PATCH) && defined(DEMAND_INDIV)
#error Only one of DEMAND_PATCH and DEMAND_INDIV should be #defined
#elif !defined(DEMAND_PATCH) && !defined(DEMAND_INDIV)
#error One of DEMAND_PATCH and DEMAND_INDIV must be #defined
#endif

// Alternative parameterisations of the convective boundary layer
//   AET_MONTEITH_HYPERBOLIC = hyperbolic parameterisation (Huntington & Monteith 1998)
//   AET_MONTEITH_EXPONENTIAL = exponential parameterisation (Monteith 1995)
// Comment out one of the following two lines:

#define AET_MONTEITH_HYPERBOLIC
//#define AET_MONTEITH_EXPONENTIAL

// Check:
#if defined(AET_MONTEITH_HYPERBOLIC) && defined(AET_MONTEITH_EXPONENTIAL)
#error Only one of AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL should be #defined
#elif !defined(AET_MONTEITH_HYPERBOLIC) && !defined(AET_MONTEITH_EXPONENTIAL)
#error One of AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL must be #defined
#endif

// Alternative parameterisations of plant water uptake

// guess2008 - drought/water uptake changes - added WR_SPECIESSPECIFIC option

//   WR_WCONT = uptake rate coupled to water content and vertical root distribution
//              (as in earlier versions of LPJ-GUESS and LPJF)
//   WR_ROOTDIST = uptake rate independent of water content (to wilting point) but
//                 with fractional uptake from different layers according to prescribed
//                 root distribution
//   WR_SMART = uptake rate independent of water content (to wilting point), fractional
//              uptake from different layers according to layer water content for
//              trees, according to prescribed root distribution for grasses
//	 WR_SPECIESSPECIFIC = uptake rate is species specific, with more drought tolerance species 
//            = (lower species_drought_tolerance values) having greater relative uptake rates. 

// Comment out all but one of the following three lines:


// guess2008 - drought/water uptake changes - added WR_SPECIESSPECIFIC option
//#define WR_WCONT
#define WR_ROOTDIST
//#define WR_SMART
//#define WR_SPECIESSPECIFIC


#if defined(WR_WCONT) && defined(WR_ROOTDIST)
#error Only one of WR_WCONT, WR_ROOTDIST and WR_SMART should be #defined
#elif defined(WR_WCONT) && defined(WR_SMART)
#error Only one of WR_WCONT, WR_ROOTDIST and WR_SMART should be #defined
#elif defined(WR_ROOTDIST) && defined(WR_SMART)
#error Only one of WR_WCONT, WR_ROOTDIST and WR_SMART should be #defined
#elif defined(WR_SPECIESSPECIFIC) && (defined(WR_SMART) || defined(WR_ROOTDIST) || defined(WR_WCONT)) // guess2008
#error Only one of WR_SPECIESSPECIFIC, WR_WCONT, WR_ROOTDIST and WR_SMART should be #defined
#elif !defined(WR_WCONT) && !defined(WR_ROOTDIST) && !defined(WR_SMART) && !defined(WR_SPECIESSPECIFIC)
#error One of WR_WCONT, WR_SPECIESSPECIFIC, WR_ROOTDIST and WR_SMART should be #defined
#endif


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES



///////////////////////////////////////////////////////////////////////////////////////
// INTERCEPTION

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

			// Storage capacity for precipitation by canopy (point scale)
			scap=climate.prec*min(indiv.lai_indiv*indiv.phen*indiv.pft.intc,0.999);

			// Fraction of day that canopy remains wet
			fwet=min(scap/pet,patch.fpc_rescale);

			// Calculate interception by this individual, and increment patch total

			indiv.intercep=fwet*pet*indiv.fpc;
			patch.intercep+=indiv.intercep;
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

	const double VSTEP=2.0; // width of vertical layers for canopy-area integration (m)
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
	int toplayer; // number of vertical layers of width VSTEP in vegetation (minus 1)
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
	Climate& climate=patch.stand.climate;

	if (vegmode==POPULATION) {
		
		// POPULATION MODE

		// Loop through individuals

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();
		
			// For this individual ...

			indiv.fpar=indiv.fpc*indiv.phen; // Eqn 1
			indiv.fpar_leafon=indiv.fpc; // Eqn 2

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

			indiv.fpar=0.0;
			indiv.fpar_leafon=0.0;
			if (indiv.height>height_veg) height_veg=indiv.height;
			plai_leafon+=indiv.lai;
			
			if (indiv.pft.lifeform==GRASS || indiv.pft.lifeform==CROP) {
				plai_leafon_grass+=indiv.lai;
				plai_grass+=indiv.lai*indiv.phen;
			}

			// Accumulate LAI-weighted sum of individual leaf-out fractions
			phen_veg+=indiv.phen*indiv.lai;

			vegetation.nextobj(); // ... on to next individual
		}

		// Calculate LAI-weighted mean leaf-out fraction for vegetation
		// guess2008 - bugfix - was: if (!negligible(plai))
		if (!negligible(plai_leafon))
			phen_veg/=plai_leafon;
		else
			phen_veg=1.0;

		// Calculate number of layers (minus 1) from ground surface to top of canopy
		toplayer=(int)(height_veg/VSTEP-0.0001);

		// Calculate FPAR by integration from the top of the canopy (Eqn 2)
		plai=0.0;
		plai_leafon=0.0;

		// Set FPAR for bottom of layer above (initially 1 at top of canopy)

		fpar_layer_bottom=1.0;
		fpar_leafon_layer_bottom=1.0;
		
		for (layer=toplayer;layer>=0;layer--) {

			lowbound=(double)layer*VSTEP;
			highbound=lowbound+VSTEP;

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
							frac-=(highbound-indiv.height)/VSTEP;
						if (indiv.boleht>lowbound)
							frac-=(indiv.boleht-lowbound)/VSTEP;

						// Calculate summed LAI of this cohort in this layer

						atoh=indiv.lai/(indiv.height-indiv.boleht);
						indiv.lai_leafon_layer=atoh*frac*VSTEP;
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

			fpar_layer_bottom=exp(-LAMBERTBEER_K*plai);
			fpar_leafon_layer_bottom=exp(-LAMBERTBEER_K*plai_leafon);

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

					else 
						indiv.fpar_leafon=0.0;

					if (!negligible(plai_layer))
						indiv.fpar+=fpar_uptake_layer*
							(indiv.lai_leafon_layer*indiv.phen)/plai_layer;
					else
						indiv.fpar=0.0;

				}

				// ... on to next individual
				vegetation.nextobj();
			}

		}

		// FPAR reaching grass canopy
		fpar_grass=exp(-LAMBERTBEER_K*plai);
		fpar_leafon_grass=exp(-LAMBERTBEER_K*plai_leafon);

		// Add grass LAI to calculate PAR reaching forest floor
		// BLARP: Order changed Ben 050301 to overcome optimisation bug in pgCC

		//plai+=plai_grass;
		fpar_ff=exp(-LAMBERTBEER_K*(plai+plai_grass));
		plai+=plai_grass;

		// Save this
		patch.fpar_ff=fpar_ff;

		plai_leafon+=plai_leafon_grass;
		fpar_leafon_ff=exp(-LAMBERTBEER_K*plai_leafon);

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

			if (indiv.pft.lifeform==GRASS || indiv.pft.lifeform==CROP) {

				// Calculate minimum FPAR for growth of this grass

				// Fraction of total grass LAI represented by this grass

				if (!negligible(plai_grass))
					flai=indiv.lai*indiv.phen/plai_grass;
				else
					flai=1.0;

				if (!negligible(climate.par))
					fpar_min=min(indiv.pft.parff_min/climate.par,1.0);
				else
					fpar_min=1.0;

				indiv.fpar=max(0.0,fpar_grass*flai-max(fpar_ff*flai,fpar_min));

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
		// exceeds 50%

		patch.fpar_grass=fpar_grass;
		par_grass=fpar_grass*climate.par;

		if (date.day==0) {
			patch.par_grass_mean=0.0;
			patch.nday_growingseason=0;
		}

		if (phen_veg>PHEN_GROWINGSEASON) {
			patch.par_grass_mean+=par_grass;
			patch.nday_growingseason++;
		}

		// Convert from sum to mean on last day of year

		if (date.islastday && date.islastmonth && patch.nday_growingseason) {
			patch.par_grass_mean/=(double)patch.nday_growingseason;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// PHOTOSYNTHESIS
// Internal function (do not call directly from framework)


// Constants required for Q10 lookup tables used by photosynthesis

const double LOOKUPQ10_MINTEMP=-70;
	// minimum temperature ever (deg C)
const double LOOKUPQ10_MAXTEMP=70;
	// maximum temperature ever (deg C)
const double LOOKUPQ10_PRECISION=0.01;
	// rounding precision for temperature in Q10 lookup tables
const int LOOKUPQ10_NDATA=(LOOKUPQ10_MAXTEMP-LOOKUPQ10_MINTEMP+1.0)/
	LOOKUPQ10_PRECISION+0.5;
	// maximum number of values to store in each lookup table
	

// Definition of Q10 lookup table class

class LookupQ10 {

private:
	double* data;

public:
	inline int element(double& temp) {

		// Returns element number corresponding to a particular temperature

		if (temp<LOOKUPQ10_MINTEMP) temp=LOOKUPQ10_MINTEMP;
		else if (temp>LOOKUPQ10_MAXTEMP) temp=LOOKUPQ10_MAXTEMP;

		return (temp-LOOKUPQ10_MINTEMP)/LOOKUPQ10_PRECISION+0.5;
	}

	LookupQ10(double q10,double base25) {
		
		// Constructor (initialises lookup table)
		
		double temp;

		data=new double[LOOKUPQ10_NDATA];
		if (!data) fail("LookupQ10: out of memory creating array");

		for (temp=LOOKUPQ10_MINTEMP;temp<=LOOKUPQ10_MAXTEMP;
			temp+=LOOKUPQ10_PRECISION) {

			data[element(temp)]=base25*pow(q10,(temp-25.0)/10.0);
		}
	}

	double& operator[](double& temp) {
		
		// "Array element" operator (returns temperature-adjusted value
		// based on Q10 and 25-degree base value)
		
		return data[element(temp)];
	}

	// guess2008 - new destructor added
	~LookupQ10() {
		
		delete[] data;
	}
};


// Constants for parameters with Q10 temperature responses used in photosynthesis
// calculations

const double Q10KO=1.2;
	// Q10 for temperature dependency of Michaelis constant for O2 (ko)
const double Q10KC=2.1;
	// Q10 for temperature dependency of Michaelis constant for CO2 (kc)
const double Q10TAU=0.57;
	// Q10 for temperature dependency of CO2/O2 specificity ratio (tau)
const double KO25=3.0E4; // value of ko at 25 deg C (Pa)
const double KC25=30.0; // value of kc at 25 deg C (Pa)
const double TAU25=2600.0; // value of tau at 25 deg C


// Lookup tables for parameters with Q10 temperature responses

LookupQ10 lookup_ko(Q10KO,KO25);
	// lookup table for Q10 temperature response of Michaelis constant for O2
LookupQ10 lookup_kc(Q10KC,KC25);
	// lookup table for Q10 temperature response of Michaelis constant for CO2
LookupQ10 lookup_tau(Q10TAU,TAU25);
	// lookup table for Q10 temperature response of CO2/O2 specificity ratio


void photosynthesis(double co2,double temp,double par,double daylength,
	double fpar,double lambda,pathwaytype pathway,double pstemp_min,
	double pstemp_low,double pstemp_high,double pstemp_max,double lambda_max,
	double& agd,double& adtmm,double& rd,double& na,bool ifnlimvmax) {

	// DESCRIPTION
	// Calculation of total daily gross photosynthesis and leaf-level net daytime
	// photosynthesis given degree of stomatal closure (as parameter lambda).
	// Includes implicit scaling from leaf to plant projective area basis.
	// Adapted from Farquhar & von Caemmerer (1982) photosynthesis model, as simplified
	// by Collatz et al (1991), Collatz et al (1992), Haxeltine & Prentice (1996a,b)
	// and Sitch et al. (2000).

	// NOTE: This function is identical to LPJF subroutine "photosynthesis" except for
	// the formulation of low-temperature inhibition coefficient tscal (tstress; LPJF).
	// The function adopted here draws down metabolic activity in approximately the
	// temperature range pstemp_min-pstemp_low but does not affect photosynthesis
	// at high temperatures.

	// HISTORY
	// Ben Smith 18/1/2001: Tested in comparison to LPJF subroutine "photosynthesis":
	// function showed identical behaviour except at temperatures >= c. 35 deg C where
	// LPJF temperature inhibition function results in lower photosynthesis.

	// INPUT PARAMETERS
	// co2          = atmospheric ambient CO2 concentration (ppmv)
	// temp         = mean air temperature today (deg C)
	// par          = total daily photosynthetically-active radiation today (J/m2/day)
	// daylength    = day length (h)
	// fpar         = fraction of PAR absorbed by foliage
	// lambda       = ratio of intercellular to ambient partial pressure of CO2
	// pathway      = biochemical pathway for photosynthesis (C3 or C4)
	// pstemp_min   = approximate low temperature limit for photosynthesis (deg C)
	// pstemp_low   = approximate lower range of temperature optimum for
	//                photosynthesis (deg C)
	// pstemp_high  = approximate upper range of temperature optimum for photosynthesis
	//                (deg C)
	// pstemp_max   = maximum temperature limit for photosynthesis (deg C)
	// lambda_max   = non-water-stressed ratio of intercellular to ambient CO2 pp

	// OUTPUT PARAMETERS
	// agd    = gross daily photosynthesis (kgC/m2/day)
	// adtmm  = leaf-level net daytime photosynthesis expressed in CO2 diffusion units
	//          (mm/m2/day)
	// rd     = leaf respiration (kgC/m2/day)
	// na     = leaf N associated with photosynthesis (kgN/m2) (GUESSN)

	// guess2008 - ALPHAA value chosen to give global carbon pool and flux values that 
	// agree with published estimates.
	const double ALPHAA=0.5;
		// scaling factor for PAR absorption from leaf to plant projective area level
		// alias "twigloss"
		// Should normally be in the range 0-1

	const double CO2_CONV=1.0E-6;
		// conversion factor for CO2 from ppmv to mole fraction
	const double PO2=2.09E4; // O2 partial pressure (Pa)
	const double PATMOS=1.0E5; // atmospheric pressure (Pa)
	const double LAMBDA_SC4=0.4;
		// 'saturation' ratio of intercellular to ambient CO2 partial pressure for C4
		// plants
	const double ALPHA_C3=0.08;
		// intrinsic quantum efficiency of CO2 uptake for C3 plants
	const double ALPHA_C4=0.053;
		// intrinsic quantum efficiency of CO2 uptake for C4 plants
	const double TMC3=45.0; // maximum temperature for C3 photosynthesis (deg C)
	const double TMC4=55.0; // maximum temperature for C4 photosynthesis (deg C)
	const double BC3=0.015;
		// leaf respiration as fraction of maximum rubisco capacity for C3 plants
	const double BC4=0.02;
		// leaf respiration as fraction of maximum rubisco capacity for C4 plants
	const double THETA=0.7; // colimitation (shape) parameter
	const double CMASS=12.0; // atomic mass of carbon
	const double CQ=4.6E-6;
		// conversion factor for solar radiation at 550 nm from J/m2 to E/m2
		// (E=mol quanta)

	// GUESSN
	const double M=25.0; // corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996a
	// end GUESSN

	double tscal; // temperature scaling coefficient
	double tk; // temperature in kelvin units
	double ko; // Michaelis constant of rubisco for O2 (Pa)
	double kc; // Michaelis constant of rubisco for CO2 (Pa)
	double tau; // CO2/O2 specificity ratio
	double gammastar; // CO2 compensation point in partial pressure units (Pa)
	double pa_co2; // ambient partial pressure of CO2 (Pa)
	double pi_co2_opt; // non-water-stressed intercellular partial pressure of CO2 (Pa)
	double c1_c3_opt; // term in photosynthesis equations
	double c2_c3_opt; // term in photosynthesis equations
	double sc3; // term in photosynthesis equations
	double sc4; // term in photosynthesis equations
	double sigma_c3; // term in photosynthesis equations
	double sigma_c4; // term in photosynthesis equations
	double vm; // rubisco capacity (gC/m2/day)
	double apar; // fraction of PAR absorbed at leaf level
	double pi_co2; // intercellular partial pressure of CO2 (Pa)
	double phi_pi;
		// factor accounting for effect of intercellular CO2 concentration on C4
		// photosynthesis
	double je; // PAR-limited photosynthetic rate (molC/m2/h)
	double jc; // rubisco-activity limited photosynthetic rate (molC/m2/h)
	double agd_g; // gross photosynthesis (gC/m2/day)
	double rd_g; // leaf respiration (gC/m2/day)
	double adt; // leaf-level net daytime photosynthesis (gC/m2/day)
	double k1; // parameter in calculation of temperature inhibition function
	double c1,c2;

	// GUESSN
	double tfac;
	double tempfac;
	double vm_max;

	// Conversion factor in calculation of leaf N: includes conversion of:
	//		- Vm from gC/m2/day to umolC/m2/sec
	//      - N from mg/m2 to kg/m2

	double CN=1.0/(3600*daylength*CMASS);

	// end GUESSN

	// No photosynthesis during polar night

	if (negligible(daylength) || negligible(fpar)) {
		agd=0.0;
		adtmm=0.0;
		rd=0.0;
		na=0.0; // GUESSN
		return;
	}

	// Convert temperature to Kelvin

	tk=temp+273.0;

	// Scale fractional PAR absorption at plant projective area level (FPAR) to
	// fractional absorption at leaf level (APAR)
	// Eqn 4, Haxeltine & Prentice 1996a

	apar=par*fpar*ALPHAA;

	// Calculate temperature-inhibition coefficient
	// This function (tscal) is mathematically identical to function tstress in LPJF.
	// In contrast to earlier versions of modular LPJ and LPJ-GUESS, it includes both
	// high- and low-temperature inhibition.

	if (temp<pstemp_max) {

		// GUESSN
		if (ifnlimvmax)	
			tscal=1.0; 
		else {
			k1=(pstemp_min+pstemp_low)/2.0;
			tscal=(1.0-0.01*exp(4.6/(pstemp_max-pstemp_high)*(temp-pstemp_high)))/
				(1.0+exp((k1-temp)/(k1-pstemp_min)*4.6));
			if (tscal<1.0e-2) tscal=0.0;
		}
		// end GUESSN
	}
	else tscal=0.0;

	if (pathway==C3) { // C3 photosynthesis

		// Calculate temperature-adjusted values of Michaelis constants of rubisco
		// for O2 and CO2, and CO2/O2 specificity ratio
		// Eqn 22, Haxeltine & Prentice 1996b

		ko=lookup_ko[temp];
		kc=lookup_kc[temp];
		tau=lookup_tau[temp];

		// Calculate CO2 compensation point (partial pressure)
		// Eqn 8, Haxeltine & Prentice 1996a

		gammastar=PO2/2.0/tau;

		// Convert ambient CO2 from ppmv to Pa

		pa_co2=co2*CO2_CONV*PATMOS;
		
		// Calculate non-water-stressed intercellular CO2 partial pressure
		// Eqn 7, Haxeltine & Prentice 1996a

		pi_co2_opt=lambda_max*pa_co2;

		// Calculation of non-water-stressed C1_C3, Eqn 4, Haxeltine & Prentice 1996a
		// High-temperature inhibition modelled by suppression of LUE by decreased
		// relative affinity of rubisco for CO2 with increasing temperature, plus
		// a step function to prohibit photosynthesis above 45 deg C (Table 3.7,
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

		if (temp<=TMC3)
			c1_c3_opt=tscal*ALPHA_C3*(pi_co2_opt-gammastar)/(pi_co2_opt+2.0*gammastar);
		else
			c1_c3_opt=0.0;

		// Calculation of non-water-stressed C2_C3, Eqn 6, Haxeltine & Prentice 1996a

		c2_c3_opt=(pi_co2_opt-gammastar)/(pi_co2_opt+kc*(1.0+PO2/ko));

		// Calculation of s, Eqn 13, Haxeltine & Prentice 1996a

		sc3=(24.0/daylength)*BC3;

		// Calculation of sigma, Eqn 12, Haxeltine & Prentice 1996a

		sigma_c3=sqrt(max(0.0,1.0-(c2_c3_opt-sc3)/(c2_c3_opt-THETA*sc3)));

		// Calculation of non-water-stressed rubisco capacity assuming leaf N not
		// limiting (Eqn 11, Haxeltine & Prentice 1996a)

		vm=1.0/BC3*c1_c3_opt/c2_c3_opt*((2.0*THETA-1.0)*sc3-
			(2.0*THETA*sc3-c2_c3_opt)*sigma_c3)*apar*CMASS*CQ;

		// GUESSN: calculate N-limited Vmax for current leaf nitrogen
		// Haxeltine & Prentice 1996b Eqn 28

		tfac=1.0/exp(-0.0693*(temp-25.0));
		vm_max=na/M/CN*tfac;

		// Calculate leaf N based on [potential] Vmax (Eqn 28 Haxeltine & Prentice 1996b)
		na=M*vm*CN/tfac;

		if (vm>vm_max && ifnlimvmax)
		  	vm=vm_max;
		// end GUESSN

		// Calculation of intercellular partial pressure of CO2 given stomatal opening
		// Eqn 7, Haxeltine & Prentice 1996a

		pi_co2=lambda*pa_co2;

		// Calculation of C1_C3, C2_C3 given actual pi

		if (temp<=TMC3)
			c1=tscal*ALPHA_C3*(pi_co2-gammastar)/(pi_co2+2.0*gammastar);
		else
			c1=0.0;

		c2=(pi_co2-gammastar)/(pi_co2+kc*(1.0+PO2/ko));

	}
	else { // C4 photosynthesis

		// Calculation of s, Eqn 13, Haxeltine & Prentice 1996a

		sc4=(24.0/daylength)*BC4;

		// Calculation of sigma, Eqn 12, Haxeltine & Prentice 1996a

		sigma_c4=sqrt(max(0.0,1.0-(1.0-sc4)/(1.0-THETA*sc4)));

		if (temp<=TMC4) {
			vm=1.0/BC4*tscal*ALPHA_C4*((2.0*THETA-1.0)*sc4-
				(2.0*THETA*sc4-1.0)*sigma_c4)*apar*CMASS*CQ;

			// GUESSN
			// Calculate leaf N based on Vmax (Eqn 28 Haxeltine & Prentice 1996b)

			na=M*vm*CN*exp(-0.0693*(temp-25.0));

			// end GUESSN
		}
		else {
			vm=0.0;
			na=0.0; // GUESSN
		}
		
		// Calculation of C1_C4, C2_C4 given actual pi
		// C1_C4 incorporates term accounting for effect of intercellular CO2
		// concentration on photosynthesis (Eqn 14, 16, Haxeltine & Prentice 1996a)

		if (temp<=TMC4) {
			phi_pi=min(lambda/LAMBDA_SC4,1.0);
			c1=tscal*phi_pi*ALPHA_C4;
		}
		else c1=0.0;

		c2=1.0;
	}

	// Calculation of PAR-limited photosynthesis rate
	// Eqn 3, Haxeltine & Prentice 1996a

	je=c1*apar*CMASS*CQ/daylength;

	// Calculation of rubisco-activity-limited photosynthesis rate
	// Eqn 5, Haxeltine & Prentice 1996a

	jc=c2*vm/24.0;

	// Calculation of daily gross photosynthesis
	// Eqn 2, Haxeltine & Prentice 1996a
	// Notes: - there is an error in Eqn 2, Haxeltine & Prentice 1996a (missing theta
	//          in 4*theta*je*jc term) which is fixed here

	agd_g=(je+jc-sqrt((je+jc)*(je+jc)-4.0*THETA*je*jc))/(2.0*THETA)*daylength;

	// Calculation of daily leaf respiration
	// Eqn 10, Haxeltine & Prentice 1996a

	if (pathway==C3)
		rd_g=BC3*vm;
	else
		rd_g=BC4*vm;

	// Calculation of leaf-level net daytime photosynthesis
	// Based on Eqn 19, Haxeltine & Prentice 1996a

	adt=agd_g-daylength/24.0*rd_g;

	// Convert to CO2 diffusion units (mm/m2/day) using ideal gas law

	adtmm=adt/CMASS*8.314*tk/PATMOS*1000.0;

	// Convert gross photosynthesis from gC to kgC units

	agd=agd_g/1000.0;
	rd=rd_g/1000.0;
}

// GUESSN - if Vmax is to be limited by actual leaf nitrogen

void photosynthesis_limvmax(double co2,double temp,double par,double daylength,
	double fpar,double fpc,double lambda,pathwaytype pathway,double pstemp_min,
	double pstemp_low,double pstemp_high,double pstemp_max,double lambda_max,
	double nmass_leaf,double cmass_leaf,
	double& agd,double& adtmm,double& rd,double& na) {

	// Call this version of photosynthesis if Vmax is to be limited by actual
	// leaf nitrogen. Performs calculations on an individual FPC basis but returns
	// output on patch basis.
	// fpar = FPAR on grid-cell basis

	double fpar_fpc,na_send;

	if (negligible(fpc)) {
		agd=0.0;
		adtmm=0.0;
		rd=0.0;
		na=0.0;
		return;
	}

	// Calculate leaf N associated with photosynthesis on FPC basis

	na_send=nmass_leaf/fpc-N0*cmass_leaf/fpc;

	if (na_send<0.0) na_send=0.0;

	// FPAR on FPC basis

	fpar_fpc=fpar/fpc;

	photosynthesis(co2,temp,par*fpar_fpc,daylength,
		1.0,lambda,pathway,pstemp_min,pstemp_low,pstemp_high,pstemp_max,
		lambda_max,agd,adtmm,rd,na_send,true);

	// Convert results back to patch basis

	agd*=fpc;
	adtmm*=fpc;
	rd*=fpc;
	na=na_send*fpc;
}

// end GUESSN


///////////////////////////////////////////////////////////////////////////////////////
// CONVECTIVE BOUNDARY LAYER
// Generic Monteith (1995) or Huntingford & Monteith (1998) CBL parameterisation

#if defined(AET_MONTEITH_EXPONENTIAL)

// Empirical parameters (exponential parameterisation)
const double ALPHAM=1.4;
const double GM=5.0;

#elif defined(AET_MONTEITH_HYPERBOLIC)

// Empirical parameters (hyperbolic parameterisation)
const double ALPHAM=1.391;
const double GM=3.26;

#endif

inline double aet_monteith(double& eet,double& gc) {

	// Returns AET given equilibrium evapotranspiration and
	// canopy conductance

#if defined(AET_MONTEITH_EXPONENTIAL)

	// Exponential version of function

	if (negligible(gc)) return 0.0;
	else return eet*ALPHAM*(1.0-exp(-gc/GM));

#elif defined(AET_MONTEITH_HYPERBOLIC)

	// Hyperbolic version of function

	return eet*ALPHAM*gc/(gc+GM);

#endif

}


inline double gc_monteith(double& aet,double& eet) {

	// Returns canopy conductance given AET and equilibrium evapotranspiration

#if defined(AET_MONTEITH_EXPONENTIAL)

	// Exponential version of function

	double t;

	if (negligible(eet)) return 0.0;
	t=aet/eet/ALPHAM;
	if (t>=1.0) fail("gc_monteith: invalid value for aet/eet/ALPHAM");

	return -GM*log(1.0-aet/eet/ALPHAM);

#elif defined(AET_MONTEITH_HYPERBOLIC)

	// Hyperbolic version of function

	return (aet*GM)/(eet*ALPHAM-aet);

#endif

}


///////////////////////////////////////////////////////////////////////////////////////
// TRANSPIRATIVE DEMAND AND NON-WATER-STRESSED PHOTOSYNTHESIS

void demand(Patch& patch) {

	// Determination of transpirative demand based on a Monteith parameterisation of
	// boundary layer dynamics, i.e. demand = f(conductance, EET) (see alternative
	// parameterisations in function aet_monteith).
	// A base value for non-water-stressed photosynthesis is calculated here (as a 
	// biproduct of the calculation of canopy conductance) and stored for reuse later.

	double adtmm;
		// leaf-level net daytime photosynthesis expressed as CO2 diffusion (mm/m2/d)
	double agd;
		// gross daily photosynthesis (kgC/m2/day)
	double rd;
		// leaf respiration (kgC/m2/day)
	double gp_patch;
		// non-water-stressed canopy conductance for patch, patch vegetated area
		// basis (mm/s)
	double gp_leafon_patch;
		// non-water-stressed canopy conductance assuming full leaf cover, patch
		// vegetated area basis (mm/s)
	double gp_indiv;
		// non-water-stressed canopy conductance for individual/cohort/population,
		// FPC basis

	// Retrieve Stand, Climate and Vegetation objects for this patch

	Stand& stand=patch.stand;
	Climate& climate=stand.climate;
	Vegetation& vegetation=patch.vegetation;

	gp_patch=0.0;
	gp_leafon_patch=0.0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		// Retrieve PFT
		Pft& pft=indiv.pft;

		if (!negligible(climate.daylength)) {

			// GUESSN
			if (iflimvmax) {

				// Call photosynthesis with actual FPAR assuming stomates fully open
				// (lambda = lambda_max)

				photosynthesis_limvmax(climate.co2,climate.temp,climate.par,climate.daylength,
					indiv.fpar_leafon,indiv.fpc,pft.lambda_max,pft.pathway,pft.pstemp_min,
					pft.pstemp_low,pft.pstemp_high,pft.pstemp_max,pft.lambda_max,
					indiv.nmass_leaf,indiv.cmass_leaf,
					agd,adtmm,rd,indiv.na_fpar);

				indiv.gp_leafon=444.4*adtmm/climate.co2/(1.0-pft.lambda_max)/
						climate.daylength+pft.gmin*indiv.fpc;

				// Save assimilation in case this turns out to be a non-water-stress day

				indiv.assim_nowstress=(agd-rd)*indiv.phen;

			}
			else {
				if (!stand.pft[pft.id].have_phot) {

					// Call photosynthesis with FPAR=1 and assuming stomates fully open
					// (lambda = lambda_max)

					photosynthesis(climate.co2,climate.temp,climate.par,climate.daylength,
						1.0,pft.lambda_max,pft.pathway,pft.pstemp_min,pft.pstemp_low,
						pft.pstemp_high,pft.pstemp_max,pft.lambda_max,agd,adtmm,rd,
						stand.pft[pft.id].na,false);

					// Eqn 21, Haxeltine & Prentice 1996
					// NB: includes conversion of daylight from hours to seconds (*3600),
					//     and CO2 from ppmv to mole fraction (*1.0e-6);
					//     scalar multiplier = 1.6 / 1.0e-6 / 3600 = 444.4
	
					stand.pft[pft.id].gpterm=444.4*adtmm/climate.co2/(1.0-pft.lambda_max)/
						climate.daylength;

					// Store net C-assimilation (gross photosynthesis minus leaf
					// respiration); valid for all individuals of this PFT given today's
					// climate and FPAR=1 assuming no water stress

					stand.pft[pft.id].assim_term=agd-rd;

					stand.pft[pft.id].have_phot=true;
				}

				// Calculate non-water-stressed canopy conductance assuming full leaf cover
				//        - include canopy-conductance component not linked to
				//          photosynthesis (diffusion through leaf cuticle etc); this is
				//          assumed to be proportional to leaf-on fraction

				indiv.gp_leafon=stand.pft[pft.id].gpterm*indiv.fpar_leafon+
					pft.gmin*indiv.fpc;
			}
			// end GUESSN


#if defined(DEMAND_INDIV)

			if (!negligible(indiv.fpc*indiv.phen)) {

				// Individual conductance and demand assuming full leaf cover,
				// FPC basis

				gp_indiv=indiv.gp_leafon*indiv.phen/indiv.fpc;
				indiv.demand=aet_monteith(patch.eet_net_veg,gp_indiv);

				// Actual conductance and demand, FPC basis

				gp_indiv=indiv.gp_leafon/indiv.fpc;
				indiv.demand_leafon=aet_monteith(patch.eet_net_veg,gp_indiv);
			}
			else {
				indiv.demand=0.0;
				indiv.demand_leafon=0.0;
			}
#endif
		}
		else {

			// special case if daylength=0
			
			indiv.gp_leafon=0.0;
			stand.pft[pft.id].assim_term=0.0;
			stand.pft[pft.id].na=0.0; // GUESSN
		}

		// Increment patch sums of non-water-stressed gp by individual value

		gp_patch+=indiv.gp_leafon*indiv.phen;
		gp_leafon_patch+=indiv.gp_leafon;

		// ... on to next individual
		vegetation.nextobj();
	}

	// Calculate transpirational demand on patch vegetated area basis
	// Eqn 23, Haxeltine & Prentice 1996

	// guess2008 - added fpc_total check
	if (!negligible(gp_patch) && !negligible(patch.fpc_total)) {
		gp_patch/=patch.fpc_total;
		patch.demand=aet_monteith(patch.eet_net_veg,gp_patch);
	}
	else
		patch.demand=0.0;

	// guess2008 - added fpc_total check
	if (!negligible(gp_leafon_patch) && !negligible(patch.fpc_total)) {
		gp_leafon_patch/=patch.fpc_total;
		patch.demand_leafon=aet_monteith(patch.eet_net_veg,gp_leafon_patch);
	}
	else
		patch.demand_leafon=0.0;
}



///////////////////////////////////////////////////////////////////////////////////////
// PLANT WATER UPTAKE

inline double water_uptake(double wcont[NSOILLAYER],double awc[NSOILLAYER],
	double rootdist[NSOILLAYER],double& emax,double& fpc_rescale,
	double fuptake[NSOILLAYER],bool ifsmart, double species_drought_tolerance) {

	// Returns plant water uptake (point scale, or mean for patch) as a fraction of
	// maximum possible (daily basis)

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

	// guess2008
	// species_drought_tolerance = used only if the WR_SPECIESSPECIFIC option is specified.
	

	// OUTPUT PARAMETER:
	//   fuptake     = fraction of total uptake originating from each layer

	double wr;
	int s;

#if defined(WR_WCONT)

	// LPJ "standard" formulation with linear scaling of uptake to water content
	// and weighting by plant root profiles

	wr=0.0;
	for (s=0;s<NSOILLAYER;s++) {
		fuptake[s]=rootdist[s]*wcont[s]*fpc_rescale;
		wr+=fuptake[s];
	}

// guess2008 - drought/water uptake changes - new option
#elif defined(WR_SPECIESSPECIFIC)

	// Uptake rate is species specific, with more drought tolerance species (lower species_drought_tolerance
	// values) having greater relative uptake rates. 
	// Reduces to WR_WCONT if species_drought_tolerance = 0.5
	
	wr=0.0;
	for (s=0;s<NSOILLAYER;s++) {
		double max_rel_uptake = pow(wcont[s],2.0*0.1); // Upper limit. Limits C3 grass uptake
		fuptake[s]=rootdist[s]*min(pow(wcont[s],2.0*species_drought_tolerance),max_rel_uptake)*fpc_rescale;
		wr+=fuptake[s];
	}

#elif defined(WR_ROOTDIST)


	// Uptake rate independent of water content (to wilting point) but with fractional
	// uptake from different layers according to prescribed root distribution

	wr=0.0;
	for (s=0;s<NSOILLAYER;s++) {
		fuptake[s]=min(wcont[s]*awc[s]*fpc_rescale,emax*rootdist[s])/emax;
		wr+=fuptake[s];
	}

#elif defined(WR_SMART)

	// Uptake rate independent of water content (to wilting point), fractional uptake
	// from different layers according to layer water content for trees, and according
	// to prescribed root distribution for grasses

	double wcsum=0.0;
	double wcfrac;

	for (s=0;s<NSOILLAYER;s++) wcsum+=wcont[s];

	wr=0.0;
	if (negligible(wcsum))
		for (s=0;s<NSOILLAYER;s++) fuptake[s]=0.0;
	else {
		for (s=0;s<NSOILLAYER;s++) {
			wcfrac=wcont[s]/wcsum;
			if (ifsmart)
				fuptake[s]=min(wcont[s]*awc[s]*wcfrac*fpc_rescale,emax*wcfrac)/emax;
			else
				fuptake[s]=min(wcont[s]*awc[s]*fpc_rescale,emax*rootdist[s])/emax;
			wr+=fuptake[s];
		}
	}

#endif

	if (!negligible(wr))
		for (s=0;s<NSOILLAYER;s++)
			fuptake[s]/=wr;

	return wr;
}




///////////////////////////////////////////////////////////////////////////////////////
// ACTUAL EVAPOTRANSPIRATION AND WATER STRESS

void aet_water_stress(Patch& patch) {

	// Supply function for evapotranspiration and determination of water stress leading
	// to down-regulation of stomatal conductance. Actual evapotranspiration determined
	// as smaller of supply and transpirative demand (see function demand).
	// Base value for actual canopy conductance calculated here for water-stressed
	// individuals and used to derive actual photosynthesis in function npp (below)

	int p;
	double wr;
	double gcbase;

	// Retrieve Climate object for patch
	Climate& climate=patch.stand.climate;

	// Calculate common point supply for each PFT in this patch

	for (p=0;p<npft;p++) {

		// Retrieve next patch PFT
		Patchpft& ppft=patch.pft[p];

		// Retrieve PFT
		Pft& pft=ppft.pft;

		if (ifdailynpp) ppft.gcbase=0.0;
		else if (date.dayofmonth==0) {

			// On first day of month, monthly mode, initialise cumulative
			// environmental drivers and counter for water-stress days

			ppft.temp_wstress=0.0;
			ppft.par_wstress=0.0;
			ppft.daylength_wstress=0.0;
			ppft.co2_wstress=0.0;
			ppft.nday_wstress=0;
			ppft.fpar_grass_wstress=0.0;
			ppft.gcbase=0.0;
		}

		// Calculate effective water supply from plant roots
		// Rescale available water by patch FPC if exceeds 1
		// (this then represents the average amount of water available over an
		// individual's FPC, assuming individuals are equal in competition for water)

		// ----------------------------------------
		// guess2008 - specieds specific drought/water uptake changes
		double species_drought_tolerance = 0.5; 
		// default, ensures that WR_SPECIESSPECIFIC gives identical results to WR_WCONT 
		
		// override with species value (always <= 0.5) iff ifspeciesspecificwateruptake == 1
		if (ifspeciesspecificwateruptake) 
			species_drought_tolerance = pft.drought_tolerance;

		wr=water_uptake(patch.soil.wcont,patch.soil.soiltype.awc,pft.rootdist,pft.emax,
			patch.fpc_rescale,ppft.fuptake,pft.lifeform==TREE,species_drought_tolerance);
		// ----------------------------------------

		// Calculate supply (Eqn 24, Haxeltine & Prentice 1996)

		ppft.supply_leafon=pft.emax*wr;
		ppft.supply=ppft.supply_leafon*ppft.phen;

		if (ppft.supply<patch.demand && !negligible(ppft.phen)) {

			// DAY WITH WATER STRESS FOR THIS PFT

			ppft.ifwstress=true;

			if (!ifdailynpp) {

				// Monthly mode - increment cumulative means for environmental
				// drivers and counter for water-stressed days

				ppft.temp_wstress+=climate.temp;
				ppft.par_wstress+=climate.par;
				ppft.daylength_wstress+=climate.daylength;
				ppft.co2_wstress+=climate.co2;
				ppft.nday_wstress++;
				ppft.fpar_grass_wstress+=patch.fpar_grass*ppft.phen;
			}

			// Calculate water-stressed canopy conductance on FPC basis assuming
			// FPAR=1 and deducting canopy conductance component not associated
			// with CO2 uptake; valid for all individuals of this PFT in this patch
			// today.
			// Eqn 25, Haxeltine & Prentice 1996

			gcbase=max(gc_monteith(ppft.supply,patch.eet_net_veg)-pft.gmin*
				ppft.supply/patch.demand,0.0);

			if (ifdailynpp)
				ppft.gcbase=gcbase;
			else
				ppft.gcbase+=gcbase;

		}
		else {

			// NO WATER STRESS TODAY FOR THIS PFT

			ppft.ifwstress=false;
		}

		// On last day of month ...

		if (date.islastday && !ifdailynpp && ppft.nday_wstress)
			ppft.gcbase/=(double)ppft.nday_wstress;
	}

	// Calculate / transfer supply to individuals

	Vegetation& vegetation=patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		Pft& pft=indiv.pft;
		Patchpft& ppft=patch.pft[pft.id];

		// Initialise on first day of month

		if (!ifdailynpp && date.dayofmonth==0) {

			indiv.fpar_wstress=0.0;
			indiv.temp_wstress=0.0;
			indiv.par_wstress=0.0;
			indiv.daylength_wstress=0.0;
			indiv.co2_wstress=0.0;
			indiv.nday_wstress=0;
		}

		indiv.supply=ppft.supply;
		indiv.supply_leafon=ppft.supply_leafon;

#if defined(DEMAND_PATCH)

		if (ppft.ifwstress) {

#elif defined(DEMAND_INDIV)
		
		if (indiv.supply<indiv.demand) {

#endif

			// WATER STRESS DAY FOR INDIVIDUAL

			indiv.ifwstress=true;

			indiv.aet=indiv.supply*indiv.fpc;

			// Record FPAR for this individual today

			if (ifdailynpp)
				indiv.fpar_wstress=indiv.fpar;
			else {
				indiv.fpar_wstress+=indiv.fpar;
				indiv.temp_wstress+=climate.temp;
				indiv.par_wstress+=climate.par;
				indiv.daylength_wstress+=climate.daylength;
				indiv.co2_wstress+=climate.co2;
				indiv.nday_wstress++;
			}
		}
		else {

			// NON-WATER STRESS DAY FOR INDIVIDUAL

			indiv.ifwstress=false;

#if defined(DEMAND_PATCH)

                        if (negligible(indiv.phen))
                                indiv.aet=0.0;
                        else
                                indiv.aet=patch.demand*indiv.fpc;

#elif defined(DEMAND_INDIV)

			if (negligible(indiv.phen))
				indiv.aet=0.0;
			else
				indiv.aet=indiv.demand*indiv.fpc;
#endif
		}

		vegetation.nextobj();
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// WATER SCALAR

void water_scalar(Patch& patch) {

	// Derivation of daily and annual versions of water scalar (wscal, omega)
	// Daily version is used to determine leaf onset and abscission for raingreen PFTs.
	// Annual version determines relative allocation to roots versus leaves for
	// subsequent year

	int p;

	Vegetation& vegetation=patch.vegetation;

	for (p=0;p<npft;p++) {

		// Retrieve next patch PFT
		Patchpft& ppft=patch.pft[p];

		// Calculate patch PFT water scalar value

		if (!negligible(patch.demand_leafon))
			ppft.wscal=min(1.0,ppft.supply_leafon/patch.demand_leafon);
		else
			ppft.wscal=1.0;

		// Update annual mean water scalar

		if (date.day==0)
			ppft.wscal_mean=ppft.wscal;
		else
			ppft.wscal_mean+=ppft.wscal;

		// Convert from sum to mean on last day of year
		if (date.islastday && date.islastmonth) ppft.wscal_mean/=365.0;

	}

	// Calculate individual water scalars

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

#if defined(DEMAND_PATCH)

		indiv.wscal=patch.pft[indiv.pft.id].wscal;

#elif defined(DEMAND_INDIV)

		if (!negligible(indiv.demand_leafon))
			indiv.wscal=min(1.0,indiv.supply_leafon/indiv.demand_leafon);
		else
			indiv.wscal=1.0;
#endif

		if (date.day==0)
			indiv.wscal_mean=indiv.wscal;
		else
			indiv.wscal_mean+=indiv.wscal;
		
		if (date.islastday && date.islastmonth)
			indiv.wscal_mean/=365.0;

		vegetation.nextobj();
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// ASSIMILATION_WSTRESS
// Internal function (do not call directly from framework)

void assimilation_wstress(Pft& pft,Patchpft& ppft,double co2,double temp,double par,
	double daylength,double fpar,double fpc,double& assim,double& na_fpar) {

	// GUESSN: this is the standard version, without N-limitation of Vmax

	// DESCRIPTION
	// Calculation of net C-assimilation under water-stressed conditions
	// (demand>supply; see function canopy_exchange_fast). Utilises a numerical
	// iteration procedure to find the level of stomatal aperture (characterised by
	// lambda, the ratio of leaf intercellular to ambient CO2 concentration) which
	// satisfies simulataneously a canopy-conductance based and light-based
	// formulation of photosynthesis (Eqns 2, 18 and 19, Haxeltine & Prentice (1996)).

	// Numerical method is a tailored implementation of the bisection method, with a
	// fixed 10 bisections, assuming root (f(lambda)=0) bracketed by f(0.02)<0 and
	// f(lambda_max+0.05)>0 (Press et al 1986)

	// To increase the efficiency with which the iteration is performed in cohort
	// and individual modes, dynamic lookup tables of class Lookup_lambda are used
	// to store output from successive calls to function photosynthesis. A separate
	// lookup table is maintained by each Patchpft object. Stored values remain valid
	// for one simulation day. This implies that function photosynthesis need be called
	// just once for a particular patch, PFT, day and value of lambda. This technique
	// will not influence computational speed in population mode, where each PFT is
	// represented by just one Individual object in each grid cell.

	// OUTPUT PARAMETER
	// assim = net assimilation on patch area basis (kgC/m2/day)

	const double EPS=0.1; // minimum precision of solution in bisection method
	const int MAXTRIES=6;
		// maximum number of iterations towards a solution in bisection method

	double gcphot;
		// canopy conductance component associated with photosynthesis on FPC basis
		// (mm/day)
	double fpar_fpc; // fraction of PAR absorbed on FPC basis
	double na; // GUESSN
	double agd; // gross daily photosynthesis (kgC/m2/day)
	double rd; // daily leaf respiration (kgC/m2/day)
	double ca; // ambient CO2 concentration in molar units
	double adt1,adt2;
	double x1,x2,xmid,rtbis,dx,fmid;
	int b;

	// Retrieve lookup table for this PFT in this patch
	Lookup_lambda& lookup_lambda=ppft.lookup_lambda;

	// Get canopy conductance component associated with photosynthesis and
	// convert from second to daily basis

	gcphot=ppft.gcbase*daylength*3600.0;

	if (negligible(fpc) || negligible(fpar) || negligible(gcphot)) {
		
		// Return zero assimilation

		assim=0.0;
		na_fpar=0.0; // GUESSN no assimilation no N associated with photosynthesis 
		return;
	}

	// Convert fpar from patch to fpc basis
	fpar_fpc=fpar/fpc;

	// convert CO2 from ppmv to mole fraction
	ca=co2*1.0e-6; 

	// Implement numerical solution

	x1=0.02; // minimum bracket of root
	x2=pft.lambda_max+0.05; // maximum bracket of root
	rtbis=x1; // root of the bisection
	dx=x2-x1;

	b=0; // number of tries so far towards solution

	fmid=EPS+1.0;

	lookup_lambda.newsearch();

	while (fabs(fmid)>EPS && b<=MAXTRIES) {

		b++;
		dx*=0.5;
		xmid=rtbis+dx; // current guess for lambda

		// Calculate total daytime photosynthesis (mm/m2/day) implied by
		// canopy conductance and current guess for lambda (xmid)
		// Eqn 18, Haxeltine & Prentice 1996

		adt1=gcphot/1.6*ca*(1.0-xmid);

		// Call function photosynthesis to calculate alternative value
		// for total daytime photosynthesis according to Eqns 2 & 19,
		// Haxeltine & Prentice (1996), and current guess for lambda
		// Value from lookup table used if available

		if (!lookup_lambda.getdata(date.year,date.day,adt2,agd,rd,na)) {	// GUESSN

			photosynthesis(co2,temp,par,daylength,1.0,xmid,pft.pathway,pft.pstemp_min,
				pft.pstemp_low,pft.pstemp_high,pft.pstemp_max,pft.lambda_max,agd,adt2,
				rd,na,false);
			
			lookup_lambda.setdata(date.year,date.day,adt2,agd,rd,na);	// GUESSN
		}

		// Evaluate fmid at the point lambda=xmid
		// fmid will be an increasing function of xmid, with a solution
		// (fmid=0) between x1 and x2

		fmid=adt2*fpar_fpc-adt1;
		
		if (fmid<0.0) {
			rtbis=xmid;
			lookup_lambda.increase();
		}
		else lookup_lambda.decrease();
	}

	// Calculate net assimilation, i.e. gross primary production minus leaf
	// respiration, including conversion from FPC to grid cell basis

	assim=(agd-rd)*fpar;

	// GUESSN: Calculate correct value of nmass_term

	na_fpar=na*fpar;

	if (na_fpar < 0.0 || na_fpar > 0.05)
		dprintf("Year %d Day %d %d na_fpar %g\n",
			date.year,date.day,1,na_fpar);
}

// GUESSN
///////////////////////////////////////////////////////////////////////////////////////
// ASSIMILATION_WSTRESS
// Internal function (do not call directly from framework)

void assimilation_wstress_limvmax(Pft& pft,Patchpft& ppft,double co2,double temp,double par,
	double daylength,double fpar,double fpc,double nmass_leaf,double cmass_leaf,
	double& assim,double& na_fpar) {

	// Version to call if Vmax reduced under leaf-level N limitation

	const double EPS=0.1; // minimum precision of solution in bisection method
	const int MAXTRIES=6;
		// maximum number of iterations towards a solution in bisection method

	double gcphot;
		// canopy conductance component associated with photosynthesis on FPC basis
		// (mm/day)
	double fpar_fpc; // fraction of PAR absorbed on FPC basis
	double agd; // gross daily photosynthesis (kgC/m2/day)
	double rd; // daily leaf respiration (kgC/m2/day)
	double ca; // ambient CO2 concentration in molar units
	double adt1,adt2;
	double x1,x2,xmid,rtbis,dx,fmid;
	int b;

	// Retrieve lookup table for this PFT in this patch
	Lookup_lambda& lookup_lambda=ppft.lookup_lambda;

	// Get canopy conductance component associated with photosynthesis and
	// convert from second to daily basis

	gcphot=ppft.gcbase*daylength*3600.0;

	if (negligible(fpc) || negligible(fpar) || negligible(gcphot)) {
		
		// Return zero assimilation

		assim=0.0;
		na_fpar=0.0; // GUESSN no assimilation no N associated with photosynthesis 
		return;
	}

	// Convert fpar from patch to fpc basis
	fpar_fpc=fpar/fpc;

	// convert CO2 from ppmv to mole fraction
	ca=co2*1.0e-6; 

	// Implement numerical solution

	x1=0.02; // minimum bracket of root
	x2=pft.lambda_max+0.05; // maximum bracket of root
	rtbis=x1; // root of the bisection
	dx=x2-x1;

	b=0; // number of tries so far towards solution

	fmid=EPS+1.0;

	lookup_lambda.newsearch();

	while (fabs(fmid)>EPS && b<=MAXTRIES) {

		b++;
		dx*=0.5;
		xmid=rtbis+dx; // current guess for lambda

		// Calculate total daytime photosynthesis (mm/m2/day) implied by
		// canopy conductance and current guess for lambda (xmid)
		// Eqn 18, Haxeltine & Prentice 1996

		adt1=gcphot/1.6*ca*(1.0-xmid);

		// Call function photosynthesis to calculate alternative value
		// for total daytime photosynthesis according to Eqns 2 & 19,
		// Haxeltine & Prentice (1996), and current guess for lambda

		photosynthesis_limvmax(co2,temp,par,daylength,fpar,fpc,xmid,pft.pathway,pft.pstemp_min,
				pft.pstemp_low,pft.pstemp_high,pft.pstemp_max,pft.lambda_max,nmass_leaf,cmass_leaf,
				agd,adt2,rd,na_fpar);
			
		// Evaluate fmid at the point lambda=xmid
		// fmid will be an increasing function of xmid, with a solution
		// (fmid=0) between x1 and x2

		fmid=adt2/fpc-adt1;
		
		if (fmid<0.0) rtbis=xmid;
	}

	// Calculate net assimilation, i.e. gross primary production minus leaf
	// respiration, including conversion from FPC to grid cell basis

	assim=agd-rd;
}

// end GUESSN


///////////////////////////////////////////////////////////////////////////////////////
// AUTOTROPHIC RESPIRATION
// Internal function (do not call directly from framework)

void respiration(double gtemp_air,double gtemp_soil,lifeformtype lifeform,
	double respcoeff,double cton_sap,double cton_root,
	double phen,double cmass_sap,double cmass_root,double assim,double& resp) {

// GUESSNFIX - Should it differ between if there is N limitation activated or not?
//void respiration(double gtemp_air,double gtemp_soil,lifeformtype lifeform,
//	double respcoeff,double cton_sap,double cton_root,
//	double phen,double nmass_sap,double nmass_root,double assim,double& resp) {

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
	// cton_sap   = PFT sapwood C:N ratio	// GUESSNFIX
	// cton_root  = PFT root C:N ratio
	// phen       = vegetation phenological state (fraction of potential leaf cover)
	// cmass_sap  = sapwood C biomass on grid cell area basis (kgC/m2) // GUESSFIX
	// cmass_root = fine root C biomass on grid cell area basis (kgC/m2)
	// nmass_sap  = sapwood N biomass on grid cell area basis (kgN/m2)  // GUESSN
	// nmass_root = fine root N biomass on grid cell area basis (kgN/m2) // GUESSN
	// assim      = net assimilation on grid cell area basis (kgC/m2/day)

	// OUTPUT PARAMETER
	// resp       = sum of maintenance and growth respiration on grid cell area basis
	//              (kgC/m2/day)

	// guess2008 - following a comment by Annett Wolf, the following parameter value was changed: 
	// const double K=0.0548; // OLD value
	const double K=0.095218;  // NEW parameter value in respiration equations 
	// See the comment after Eqn (4) below.

	double resp_sap;    // sapwood respiration (kg/m2/day)
	double resp_root;   // root respiration (kg/m2/day)
	double resp_growth; // growth respiration (kg/m2/day)

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

	if (lifeform==TREE) {

		// Sapwood respiration (Eqn 7)

		resp_sap=respcoeff*K*cmass_sap/cton_sap*gtemp_air;
		//resp_sap=respcoeff*K*nmass_sap*gtemp_air; // GUESSNFIX

		// Root respiration (Eqn 7)
		// Assumed that root phenology follows leaf phenology

		resp_root=respcoeff*K*cmass_root/cton_root*gtemp_soil*phen;
		//resp_root=respcoeff*K*nmass_root*gtemp_soil*phen; // GUESSN

		// Growth respiration = 0.25 ( GPP - maintenance respiration)

		resp_growth=(assim-resp_sap-resp_root)*0.25;
		
		// guess2008 - disallow negative growth respiration 
		// (following a comment (060823) from Annett Wolf)
		if(resp_growth<0.0) resp_growth = 0.0;

		// Total respiration is sum of maintenance and growth respiration

		resp=resp_sap+resp_root+resp_growth;
	}
	else if (lifeform==GRASS || lifeform==CROP) {

		// Root respiration

		resp_root=respcoeff*K*cmass_root/cton_root*gtemp_soil*phen;
		//resp_root=respcoeff*K*nmass_root*gtemp_soil*phen; // GUESSN

		// Growth respiration (see above)

		resp_growth=(assim-resp_root)*0.25;

		// guess2008 - disallow negative growth respiration 
		// (following a comment (060823) from Annett Wolf)
		if(resp_growth<0.0) resp_growth = 0.0;

		// Total respiration (see above)

		resp=resp_root+resp_growth;
	}
	else fail ("Canopy exchange function respiration: unknown life form");
}

///////////////////////////////////////////////////////////////////////////////////////
// NET PRIMARY PRODUCTIVITY

void npp(Patch& patch) {

	// Determination of daily NPP. Leaf level net assimilation calculated for non-
	// water-stressed individuals (i.e. with fully-open stomata) using base value
	// from function demand (above); for water-stressed individuals using base value
	// for canopy conductance by a simultaneous solution of light-based and canopy
	// conductance-based equations for net daily photosynthesis (see function
	// assimilation wstress above). The latter uses the PFT-specific base value for
	// conductance from function aet_water_stress (above).
	// Plant respiration obtained by a call to function respiration (above).
	
	double assim; // leaf-level net assimilation today

	// GUESSN
	double leafn; // actual leaf N today (kgN/m2)
	double na_fpar;
	// end GUESSN

	// Retrieve Vegetation, Stand and Climate objects for patch

	Vegetation& vegetation=patch.vegetation;
	Stand& stand=patch.stand;
	Climate& climate=stand.climate;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		// For this individual ...

		// Retrieve PFT and patch PFT

		Pft& pft=indiv.pft;
		Patchpft& ppft=patch.pft[pft.id];

		// GUESSN: Phenology
		if (!negligible(indiv.phen)) 
			indiv.nday_leafon++;
		// end GUESSN

		if (ifdailynpp) {

			// DAILY NPP MODE

			// GUESSN
			if (indiv.ifwstress) {

				// Water-stress day - derive assimilation by simultaneous solution
				// of light- and conductance-based equations of photosynthesis

				// GUESSN
				if (iflimvmax) {
					
					assimilation_wstress_limvmax(pft,ppft,climate.co2,climate.temp,
						climate.par,climate.daylength,indiv.fpar,indiv.fpc,
						indiv.nmass_leaf,indiv.cmass_leaf,indiv.assim,na_fpar);
				}
				else {

					assimilation_wstress(pft,ppft,climate.co2,climate.temp,climate.par,
						climate.daylength,indiv.fpar,indiv.fpc,indiv.assim,na_fpar);
				}
			}
			else {

				// Non-water-stress day 
				
				// GUESSN
				if (iflimvmax) {
				
					// Assimilation as calculated in function demand
					// (assuming no water stress)

					indiv.assim=indiv.assim_nowstress;
					na_fpar=indiv.na_fpar;

					if (na_fpar <= 0.0 || na_fpar > 0.05)
						dprintf("Year %d Day %d %d na_fpar %g\n",
							date.year,date.day,2,na_fpar);

				}
				else {

					// Use base value for non-water-stressed
					// assimilation, scaling to patch by FPAR
				
					indiv.assim=stand.pft[pft.id].assim_term*indiv.fpar;

					// GUESSN
					// Calculate leaf N today

					na_fpar=stand.pft[pft.id].na*indiv.fpar;
				}
			}

			// Calculate leaf nitrogen today on patch area basis

			if (!negligible(indiv.phen)) {
				leafn=na_fpar+N0*indiv.cmass_leaf; // (N0 is defined near top of canexch.cpp)
				indiv.leafn+=leafn;
			}
			// end GUESSN

			// Calculate respiration response to air and soil temperature
			// (if not already known for this day)

			if (climate.last_gtemp!=date.day) {
				respiration_temperature_response(climate.temp,climate.gtemp);
				climate.last_gtemp=date.day;
			}

			if (patch.soil.last_gtemp!=date.day) {
				respiration_temperature_response(patch.soil.temp,patch.soil.gtemp);
				patch.soil.last_gtemp=date.day;
			}

			// Calculate autotrophic respiration

			respiration(climate.gtemp,patch.soil.gtemp,indiv.pft.lifeform,
				indiv.pft.respcoeff,indiv.pft.cton_sap,indiv.pft.cton_root,
				indiv.phen,indiv.cmass_sap,indiv.cmass_root,indiv.assim,indiv.resp);

		// GUESSNFIX
		//	respiration(climate.gtemp,patch.soil.gtemp,indiv.pft.lifeform,
		//		indiv.pft.respcoeff,indiv.pft.cton_sap,indiv.pft.cton_root,
		//		indiv.phen,indiv.nmass_sap,indiv.nmass_root,indiv.assim,indiv.resp);

			// Update accumulated annual NPP and daily vegetation-atmosphere flux

			indiv.dassim[date.day]=indiv.assim;	// GUESSN
			indiv.anpp+=indiv.assim-indiv.resp;

			// guess2008
			if (indiv.alive)
				patch.fluxes.dcflux_veg+=indiv.resp-indiv.assim;

			// Monthly NPP and LAI

			indiv.mnpp[date.month]+=indiv.assim-indiv.resp;
			// guess2008 - changed indiv.phen_mean to indiv.phen here. mlai is always 0 otherwise 
			indiv.mlai[date.month]+=indiv.lai*indiv.phen;

			// guess2008 - update monthly arrays
			indiv.mgpp[date.month]+=indiv.assim;
			indiv.mra[date.month]+=indiv.resp;
			patch.fluxes.mcflux_gpp[date.month]+=indiv.assim;
			patch.fluxes.mcflux_ra[date.month]+=indiv.resp;

			// On last day of month - convert monthly LAI from sum to mean

			if (date.islastday) {
				indiv.mlai[date.month]/=(double)date.ndaymonth[date.month];

				// GUESSN: Monthly leaf N

				if (indiv.nday_leafon!=0.0) 
					indiv.leafn/=(double)indiv.nday_leafon;
				else
					indiv.leafn=0.0;

				if (indiv.leafn>indiv.leafn_max)
					indiv.leafn_max=indiv.leafn;

				indiv.leafn_mean+=indiv.leafn/12.0;

				// Reset leaf-on-day counter for following month
				indiv.nday_leafon=0;

				// end GUESSN
			}
		}
		else {

			// MONTHLY NPP MODE

			// Accumulate fractional leaf cover for month

			if (date.dayofmonth==0)
				indiv.phen_mean=indiv.phen;
			else
				indiv.phen_mean+=indiv.phen;

			// Non-water-stressed photosynthesis - use daily value and scale to patch
			// by FPAR

			if (!indiv.ifwstress) {

				// GUESSN
				if (iflimvmax) {

					// Assimilation as calculated in function demand
					// (assuming no water stress)

					indiv.assim+=indiv.assim_nowstress;
					na_fpar=indiv.na_fpar;

					if (na_fpar <= 0.0 || na_fpar > 0.05)
						dprintf("Year %d Day %d %d na_fpar %g\n",
							date.year,date.day,3,na_fpar);
				}
				else {

					indiv.assim+=stand.pft[pft.id].assim_term*indiv.fpar;
					na_fpar=stand.pft[pft.id].na*indiv.fpar;

				}

				// Calculate leaf nitrogen today on patch area basis

				if (!negligible(indiv.phen)) {

					leafn=na_fpar+N0*indiv.cmass_leaf; // (N0 is defined near top of canexch.cpp
					
					indiv.leafn+=leafn;
				}
				// end GUESSN
			}

			if (date.islastday) {

				// On last day of month

				// Convert fractional leaf cover to mean for this month
				indiv.phen_mean/=(double)date.ndaymonth[date.month];
				
				if (indiv.nday_wstress) {

					// Water-stressed photosynthesis

					// Obtain means for drivers of photosynthesis over water-stress
					// days this month

					double frac_wstress=1.0/indiv.nday_wstress;

					indiv.fpar_wstress*=frac_wstress;
					indiv.temp_wstress*=frac_wstress;
					indiv.par_wstress*=frac_wstress;
					indiv.daylength_wstress*=frac_wstress;
					indiv.co2_wstress*=frac_wstress;

					// Calculate mean water-stressed photosynthesis for water-stress
					// days this month by simulataneous solution of light- and
					// conductance-based equations for photosynthesis

					// GUESSN
					 if (iflimvmax) {
					
						 assimilation_wstress_limvmax(pft,ppft,indiv.co2_wstress,
							indiv.temp_wstress,indiv.par_wstress,indiv.daylength_wstress,
							indiv.fpar_wstress,indiv.fpc,indiv.nmass_leaf,indiv.cmass_leaf,
							assim,na_fpar);
					 }
					 else {

						assimilation_wstress(indiv.pft,ppft,indiv.co2_wstress,
							indiv.temp_wstress,indiv.par_wstress,indiv.daylength_wstress,
							indiv.fpar_wstress,indiv.fpc,assim,na_fpar);
					 }
					 // end GUESSN

					// Convert from mean to sum over water-stress-days
					indiv.assim+=assim*(double)indiv.nday_wstress;

					// GUESSN: Calculate leaf N today
					leafn=na_fpar+N0*indiv.cmass_leaf;

					indiv.leafn+=leafn*(double)indiv.nday_wstress;
					// end GUESSN
				}

				// Calculate respiration response to mean monthly air and soil temperature
				// (if not already known for this month)

				if (climate.last_mgtemp!=date.month) {
					respiration_temperature_response(climate.mtemp,climate.mgtemp);
					climate.last_mgtemp=date.month;
				}

				if (patch.soil.last_mgtemp!=date.month) {
					respiration_temperature_response(patch.soil.mtemp,patch.soil.mgtemp);
					patch.soil.last_mgtemp=date.month;
				}

				// Calculate autotrophic respiration

				assim=indiv.assim/(double)date.ndaymonth[date.month];
					// average daily assimilation for this month

				respiration(climate.mgtemp,patch.soil.mgtemp,indiv.pft.lifeform,
					indiv.pft.respcoeff,indiv.pft.cton_sap,indiv.pft.cton_root,
					indiv.phen_mean,indiv.cmass_sap,indiv.cmass_root,assim,indiv.resp);


				// GUESSNFIX
			//	respiration(climate.mgtemp,patch.soil.mgtemp,indiv.pft.lifeform,
			//		indiv.pft.respcoeff,indiv.pft.cton_sap,indiv.pft.cton_root,
			//		indiv.phen_mean,indiv.nmass_sap,indiv.nmass_root,assim,indiv.resp);


				indiv.resp*=(double)date.ndaymonth[date.month];

				// Update accumulated annual NPP and daily vegetation-atmosphere flux

				// GUESSN
				// Distribute monthly assimilation through all days of this month
				// (required by SOM dynamics to distribute annual N uptake at end of year)
				for (int d=0;d<date.ndaymonth[date.month];d++)
					indiv.dassim[date.day]=assim;
				// end GUESSN

				indiv.anpp+=indiv.assim-indiv.resp;

				// guess2008
				if (indiv.alive) // Ben 2007-11-28	
					patch.fluxes.dcflux_veg+=indiv.resp-indiv.assim;

				// Monthly NPP and LAI

				indiv.mnpp[date.month]=indiv.assim-indiv.resp;
				indiv.mlai[date.month]=indiv.lai*indiv.phen_mean;

				// guess2008 - update monthly arrays
				indiv.mgpp[date.month]+=indiv.assim;
				indiv.mra[date.month]+=indiv.resp;
				patch.fluxes.mcflux_gpp[date.month]+=indiv.assim; // ANDERS A TRENDY
				patch.fluxes.mcflux_ra[date.month]+=indiv.resp; // ANDERS A TRENDY

				// GUESSN: Monthly leaf N
				if (indiv.nday_leafon)
					indiv.leafn/=(double)indiv.nday_leafon;
				if (indiv.leafn>indiv.leafn_max) {
					indiv.leafn_max=indiv.leafn;
				}

				indiv.leafn_mean+=indiv.leafn/12.0;
				// end GUESSN

				// Reinitialise for next month
				indiv.assim=0.0;
				indiv.fpar_wstress=0.0;
				indiv.leafn=0.0; // GUESSN
				indiv.nday_leafon=0; // GUESSN
			}
		}

		vegetation.nextobj();
	}

	// Update annual and monthly vegetation-atmosphere flux

	patch.fluxes.acflux_veg+=patch.fluxes.dcflux_veg;
	patch.fluxes.mcflux_veg[date.month]+=patch.fluxes.dcflux_veg;
}


///////////////////////////////////////////////////////////////////////////////////////
// FOREST-FLOOR CONDITIONS

void forest_floor_conditions(Patch& patch) {

	// DESCRIPTION
	// Called in cohort/individual mode (not population mode) to quantify growth
	// conditions at the forest floor for each PFT

	double netps_ff;
		// net daily photosynthesis (kgC/m2/day)
	double agd;
		// gross daily photosynthesis (kgC/m2/day)
	double rd;
		// daily leaf respiration (kgC/m2/day)
	double adtmm;
		// leaf-level net daytime photosynthesis expressed as CO2 diffusion (mm/m2/d)
	int p;
	double na; // GUESSN

	// Retrieve Stand and Climate objects for patch

	Stand& stand=patch.stand;
	Climate& climate=stand.climate;

	// Loop through PFTs

	for (p=0;p<npft;p++) {

		// Retrieve patch PFT
		Patchpft& ppft=patch.pft[p];

		// Initialise net photosynthesis sum on first day of year
		if (date.day==0) ppft.anetps_ff=0.0;

		// WATER-STRESSED ASSIMILATION

		if (ifdailynpp && ppft.ifwstress) {

			// Daily mode

			// GUESSN
			if (iflimvmax) {
			
				assimilation_wstress_limvmax(ppft.pft,ppft,climate.co2,climate.temp,climate.par,
					climate.daylength,patch.fpar_grass*ppft.phen,1.0,1.0,ppft.pft.cton_leaf,
					netps_ff,na);
			}
			else {
				 
				assimilation_wstress(ppft.pft,ppft,climate.co2,climate.temp,
					climate.par,climate.daylength,patch.fpar_grass*ppft.phen,
					1.0,netps_ff,na);
			}
			// end GUESSN
	
			ppft.anetps_ff+=netps_ff;
		}
		else if (date.islastday && ppft.nday_wstress) {

			// Monthly mode
			
			ppft.temp_wstress/=(double)ppft.nday_wstress;
			ppft.par_wstress/=(double)ppft.nday_wstress;
			ppft.daylength_wstress/=(double)ppft.nday_wstress;
			ppft.fpar_grass_wstress/=(double)ppft.nday_wstress;
			ppft.co2_wstress/=(double)ppft.nday_wstress;

			// GUESSN
			if (iflimvmax) {
			
				assimilation_wstress_limvmax(ppft.pft,ppft,climate.co2,climate.temp,climate.par,
					climate.daylength,patch.fpar_grass*ppft.phen,1.0,1.0,ppft.pft.cton_leaf,
					netps_ff,na);
			}
			else {
				assimilation_wstress(ppft.pft,ppft,ppft.co2_wstress,
					ppft.temp_wstress,ppft.par_wstress,ppft.daylength_wstress,
					ppft.fpar_grass_wstress,1.0,netps_ff,na);
			}
			// end GUESSN
			
			ppft.anetps_ff+=netps_ff*(double)ppft.nday_wstress;
		}

		// NON-WATER-STRESSED ASSIMILATION

		if (!ppft.ifwstress) {

			Pft& pft=patch.pft[p].pft;

			// GUESSN
			if (iflimvmax) {				
			
				photosynthesis_limvmax(climate.co2,climate.temp,climate.par,climate.daylength,
					patch.fpar_grass*ppft.phen,1.0,pft.lambda_max,pft.pathway,pft.pstemp_min,
					pft.pstemp_low,pft.pstemp_high,pft.pstemp_max,pft.lambda_max,
					1.0,pft.cton_leaf,agd,adtmm,rd,na);

				ppft.anetps_ff+=agd-rd;
			}
			else {
				if (!stand.pft[p].have_phot) {

					// Call photosynthesis with FPAR=1 and assuming stomates fully open
					// (lambda = lambda_max)

					photosynthesis(climate.co2,climate.temp,climate.par,climate.daylength,
						1.0,pft.lambda_max,pft.pathway,pft.pstemp_min,pft.pstemp_low,
						pft.pstemp_high,pft.pstemp_max,pft.lambda_max,agd,adtmm,rd,na,false);

					// Store net C-assimilation (gross photosynthesis minus leaf
					// respiration); valid for all individuals of this PFT given today's
					// climate and FPAR=1 assuming no water stress

					stand.pft[pft.id].assim_term=agd-rd;
				}

				// Calculate net assimilation at top of grass canopy (or at soil surface
				// if there is none)

				ppft.anetps_ff+=stand.pft[p].assim_term*patch.fpar_grass*
					ppft.phen;
			}

			// end GUESSN
		}

		// On last day of year ...

		if (date.islastday && date.islastmonth) {

			// guess2008 - avoid negative ppft.anetps_ff
			if(ppft.anetps_ff < 0.0) ppft.anetps_ff = 0.0;

			if (ppft.anetps_ff>patch.stand.pft[p].anetps_ff_max)
				patch.stand.pft[p].anetps_ff_max=ppft.anetps_ff;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// CANOPY EXCHANGE
// Should be called each simulation day for each modelled area or patch, following
// update of leaf phenology and soil temperature and prior to update of soil water.


void canopy_exchange(Patch& patch) {

	// DESCRIPTION
	// Vegetation-atmosphere exchange of CO2 and water including calculations
	// of actual evapotranspiration (AET), canopy conductance, carbon assimilation
	// and autotrophic respiration.

	// This "fast" function attempts to maximise computational efficiency, mainly by
	// restricting the number of calls to the computation-intensive function
	// photosynthesis. Unfortunately this is done somewhat at the expense of code
	// readability.

	// NB: The global variable 'ifdailynpp' determines whether respiration and
	// some photosynthesis calculations are performed every day, or on the last day of
	// each month, based on average conditions for the month (the latter mode is much
	// faster).

	// In "monthly" mode, following LPJF, AET and canopy conductance are calculated
	// daily, while carbon assimilation is calculated daily under non-water-stress
	// conditions (when transpirational demand for water is met by plant-regulated
	// supply) and on the last day of the month for water-stress days (demand>supply).
	// This results in the fastest overall simulation speed, since non-water-stressed
	// photosynthesis is calculated as a biproduct of daily AET calculations, while
	// water-stressed photosynthesis is calculated using a computationally-intensive
	// numerical iteration procedure. Note that, in monthly mode, the accumulated
	// annual NPP and vegetation C flux values are valid only for the last day of each
	// month, following the call to this function. Daily flux values are never valid in
	// "monthly" mode.

	// If you require daily output, use "daily" mode

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

	Vegetation& vegetation=patch.vegetation;
	Climate& climate=patch.stand.climate;

	double pet_s;
		// potential evapotranspiration over non-vegetated parts of patch (mm,
		// patch basis)
	double pet_patch;
		// total potential evapotranspiration for patch
	int m;

	if (date.day==0) {
		
		// On first day of year ...

		// Calculate total FPC and initialise sums for each individual

		patch.fpc_total=0.0;
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			patch.fpc_total+=indiv.fpc;
			indiv.anpp=0.0;

			// GUESSN
			indiv.leafn=0.0;
			indiv.leafn_max=0.0;
			indiv.leafn_mean=0.0;
			indiv.nday_leafon=0;
			// end GUESSN
 
			for (m=0;m<12;m++) {
				indiv.mnpp[m]=0.0;
				indiv.mlai[m]=0.0;

				// guess2008 - initialise
				indiv.mgpp[m]=0.0;
				indiv.mra[m]=0.0;
			}

			vegetation.nextobj();
		}

		// Calculate rescaling factor to account for overlap between populations/
		// cohorts/individuals (i.e. total FPC > 1)

		if (patch.fpc_total>1.0)
			patch.fpc_rescale=1.0/patch.fpc_total;
		else
			patch.fpc_rescale=1.0;
	}

	// Canopy exchange processes

	interception(patch,climate);
	fpar(patch);
	demand(patch);
	aet_water_stress(patch);
	water_scalar(patch);
	npp(patch);
	forest_floor_conditions(patch);

	// Interception for patch

	patch.aintercep+=patch.intercep;
	patch.mintercep[date.month]+=patch.intercep;

	// Potential evapotranspiration for patch

	pet_s=climate.eet*PRIESTLEY_TAYLOR*max(1.0-patch.fpc_total,0.0);
	pet_patch=pet_s+patch.demand*patch.fpc_total+patch.intercep;
	patch.apet+=pet_patch;
	patch.mpet[date.month]+=pet_patch;
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
// Farquhar GD & von Caemmerer 1982 Modelling of photosynthetic response to
//   environmental conditions. In: Lange, OL, Nobel PS, Osmond CB, Ziegler H
//   (eds) Physiological Plant Ecology II: Water Relations and Carbon
//   Assimilation, Vol 12B. Springer, Berlin, pp 549-587.
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
