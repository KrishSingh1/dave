///////////////////////////////////////////////////////////////////////////////////////
/// \file soilwater.cpp
/// \brief Soil hydrology and snow
///
/// Version including evaporation from soil surface, based on work by Dieter Gerten, 
/// Sibyll Schaphoff and Wolfgang Lucht, Potsdam
///
/// Also includes modified code compatible with "fast" cohort/
/// individual mode - see canexch.cpp
/// Includes weather generator and Dieter G:s latest updates
/// Includes baseflow runoff
/// Note: error in runoff calculation corrected 030120
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
#include "soilwater.h"


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANT

const bool ifevap=true; // whether to model evaporation from soil surface


///////////////////////////////////////////////////////////////////////////////////////
// SNOW
// Internal function (do not call directly from framework)

void snow(double prec,double temp,double& snowpack,double& rain,double& melt) {

	// DESCRIPTION
	// Daily calculation of snowfall and rainfall from precipitation and snow melt from
	// snow pack; update of snow pack with new snow and snow melt and snow melt

	// INPUT PARAMETERS
	// prec = precipitation today (mm)
	// temp = air temperature today (deg C)
	
	// INPUT AND OUTPUT PARAMETER
	// snowpack = stored snow (rainfall mm equivalents)

	// OUTPUT PARAMETERS
	// rain = rainfall today (mm)
	// melt = snow melt today (mm)

	const double TSNOW=0.0;
		// maximum temperature for precipitation as snow (deg C)
		// previously 2 deg C; new value suggested by Dieter Gerten 2002-12
	const double KM=3.0;
		// coefficient in snowmelt function
		// changed from 0.7, S. Sitch, pers. comm, 2001-11-28
	const double SNOWPACK_MAX=10000.0;
		// maximum size of snowpack (mm) (S. Sitch, pers. comm. 2001-11-28)

	double inc;

	if (temp<TSNOW) { // snowing today
		inc=min(prec,SNOWPACK_MAX-snowpack);
		snowpack+=inc;
		rain=prec-inc;
		melt=0.0;
	}
	else { // melting today
		rain=prec;

		// New snow melt formulation
		// Dieter Gerten 021121
		// Ref: Choudhury et al 1998

		melt=min((1.5+0.007*prec)*(temp-TSNOW),snowpack);

		snowpack-=melt;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// HYDROLOGY
// Internal function (do not call directly from framework)

void hydrology_lpjf(Patch& patch,double pet,double rain,double melt,
	double perc_base,double perc_exp,double awc[NSOILLAYER],double fevap,
	double awcont[NSOILLAYER],double wcont[NSOILLAYER],double& wcont_evap,double& runoff,double snowpack) {


	// DESCRIPTION
	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET) and percolation between layers;
	// calculation of runoff

	// guess2008 - DLE - new function signature. We now take awcont[] as an argument

	// INPUT PARAMETERS
	// pet        = potential evapotranspiration today (mm)
	// rain       = rainfall today (mm)
	// melt       = snowmelt today (mm)
	// perc_base  = coefficient in percolation calculation (K in Eqn 31, Haxeltine
	//              & Prentice 1996)
	// perc_exp   = exponent in percolation calculation (=4 in Eqn 31, Haxeltine &
	//              Prentice 1996)
	// awc        = array containing available water holding capacity of soil
	//              layers (mm rainfall) [0=upper layer]
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// INPUT AND OUTPUT PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont     = wcont averaged over the growing season - guess2008


	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	const double SOILDEPTH_EVAP=200.0;
		// depth of sublayer at top of upper soil layer, from which evaporation is
		// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double BASEFLOW_FRAC=0.5;
		// Fraction of standard percolation amount from lower soil layer that is
		// diverted to baseflow runoff
	const double K_DEPTH=0.4;
	const double K_AET=0.52;
		// Fraction of total (vegetation) AET from upper soil layer that is derived
		// from the top K_DEPTH (fraction) of the upper soil layer
		// (parameters for calculating K_AET_DEPTH below)
	const double K_AET_DEPTH=(SOILDEPTH_UPPER/SOILDEPTH_EVAP-1.0)*
		(K_AET/K_DEPTH-1.0)/(1.0/K_DEPTH-1.0)+1.0;
		// Weighting coefficient for AET flux from evaporation layer, assuming active
		//   root density decreases with soil depth
		// Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
		//   K_DEPTH=0.4, K_AET=0.52

	int s;
	double aet; // AET for a particular layer and individual (mm)
	double aet_layer[NSOILLAYER]; // total AET for each soil layer (mm)
	double aet_total; // total AET (mm)
	double runoff_surf; // runoff from upper soil layer (mm)
	double runoff_drain; // runoff (drainage) from lower soil layers (mm)
	double runoff_baseflow; // base flow (mm)
	double evap; // evaporation from soil surface (mm)
	double perc,perc_frac,perc_baseflow;
	double influx; // inward water flux to soil today (rain + snowmelt) (mm)

	// Retrieve Vegetation for this patch
	Vegetation& vegetation=patch.vegetation;

	for (s=0;s<NSOILLAYER;s++) aet_layer[s]=0.0;
	aet_total=0.0;

	// Sum AET for across all vegetation individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		for (s=0;s<NSOILLAYER;s++) {
			aet=patch.pft[indiv.pft.id].fuptake[s]*indiv.aet;
			aet_layer[s]+=aet;
			aet_total+=aet;
		}
		vegetation.nextobj();
	}

	// Evaporation from soil surface
	
	// guess2008 - changed to wcont_evap**2, as in LPJ-mL 
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	if (snowpack < 10.0) // i.e. evap only if snow depth < 10mm
		evap=pet*PRIESTLEY_TAYLOR*wcont_evap*wcont_evap*fevap;
	else
		evap = 0.0;

	// Implement in- and outgoing fluxes to upper soil layer
	// BLARP: water content can become negative, though apparently only very slightly
	//    - quick fix implemented here, should be done better later

	wcont[0]+=(rain+melt-aet_layer[0]-evap)/awc[0];
	if (wcont[0]!=0.0 && wcont[0] < 0.0001) // guess2008 - bugfix
		wcont[0]=0.0;


	// Surface runoff

	if (wcont[0]>1.0) {
		runoff_surf=(wcont[0]-1.0)*awc[0];
		wcont[0]=1.0;
	}
	else runoff_surf=0.0;

	// Update water content in evaporation layer for tomorrow

	influx=rain+melt;

	wcont_evap+=(influx-aet_layer[0]*SOILDEPTH_EVAP*K_AET_DEPTH/SOILDEPTH_UPPER-evap)
		/awc[0];

	if (wcont_evap>wcont[0]) wcont_evap=wcont[0];

	// Percolation from evaporation layer

	if (influx>0.1)
		perc=min(SOILDEPTH_EVAP/SOILDEPTH_UPPER*perc_base*pow(wcont_evap,perc_exp),
			influx);
	else
		perc=0.0;

	wcont_evap-=perc/awc[0];

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	runoff_drain=0.0;

	for (s=1;s<NSOILLAYER;s++) {

		// Percolation
		// Allow only on days with rain or snowmelt (Dieter Gerten, 021216)

		if (influx>=0.1)
			perc=min(perc_base*pow(wcont[s-1],perc_exp),influx);
		else
			perc=0.0;

		perc_frac=min(perc/awc[s-1],wcont[s-1]);
		
		wcont[s-1]-=perc_frac;
		wcont[s]+=perc_frac*awc[s-1]/awc[s];
		if (wcont[s]>1.0) {
			runoff_drain+=(wcont[s]-1.0)*awc[s];
			wcont[s]=1.0;
		}

		// Deduct AET from this soil layer
		// BLARP! Quick fix here to prevent negative soil water

		wcont[s]-=aet_layer[s]/awc[s];
		if (wcont[s]<0.0) wcont[s]=0.0;
	}

	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)

	if (influx>=0.1) {
		perc_baseflow=BASEFLOW_FRAC*perc_base*pow(wcont[NSOILLAYER-1],perc_exp);
		// guess2008 - Added "&& influx >= runoff_surf" to guarantee nonnegative baseflow.
		if (perc_baseflow>influx-runoff_surf && influx >= runoff_surf) perc_baseflow=influx-runoff_surf;

		// Deduct from water content of bottom soil layer

		perc_frac=min(perc_baseflow/awc[NSOILLAYER-1],wcont[NSOILLAYER-1]);
		wcont[NSOILLAYER-1]-=perc_frac;
		runoff_baseflow=perc_frac*awc[NSOILLAYER-1];
	}
	else runoff_baseflow=0.0;

	runoff=runoff_surf+runoff_drain+runoff_baseflow;


	patch.arunoff+=runoff;
	patch.aaet+=aet_total;
	patch.aevap+=evap;

	patch.maet[date.month]+=aet_total;
	patch.mevap[date.month]+=evap;
	patch.mrunoff[date.month]+=runoff;


	// guess2008 - DLE - update awcont
	// Original algorithm by Thomas Hickler
    for (s=0;s<NSOILLAYER;s++) {
        
		// Reset the awcont array on the first day of every year
		if (date.day==0) {
            awcont[s]=0.0;
            if (s==0) patch.growingseasondays=0;
        }

		// If it's warm enough for growth, update awcont with this day's wcont
        if (patch.stand.gridcell.climate.temp>5.0) {
            awcont[s]+=wcont[s];
            if (s==0) patch.growingseasondays++;
        }

		// Do the averaging on the last day of every year
        if (date.islastday && date.islastmonth)
            awcont[s]/=(double)patch.growingseasondays;
		
		// In case it's never warm enough:
		if (patch.growingseasondays<1)
            awcont[s]=1.0;
    }

}

///////////////////////////////////////////////////////////////////////////////////////
// SOIL WATER
// Call this function each simulation day for each modelled area or patch, following
// calculation of vegetation production and evapotranspiration and before soil organic
// matter and litter dynamics

void soilwater(Climate& climate,Patch& patch) {

	// DESCRIPTION
	// Performs daily accounting of soil water and snow pack

	double rain; // rainfall today (mm)
	double melt; // snow melt today (mm)
	double fpc_phen_total; // phenology-weighted FPC sum for patch

	// Retrieve Vegetation and Soil objects for patch

	Vegetation& vegetation=patch.vegetation;
	Soil& soil=patch.soil;

	// Update snowpack and derive actual water input to soil, taking into account
	// interception and snowmelt

	snow(climate.prec-patch.intercep,climate.temp,soil.snowpack,rain,melt);

	// Sum vegetation phenology-weighted FPC

	fpc_phen_total=0.0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		fpc_phen_total+=indiv.fpc*indiv.phen;

		vegetation.nextobj(); // ... on to next individual
	}

	// Hydrology
	// Fraction of grid cell subject to evaporation from soil surface is
	// complement of summed vegetation projective cover (FPC)

	// guess2008 - DLE - added soil.awcont & soil.snowpack to the function call
	hydrology_lpjf(patch,climate.eet,rain,melt,soil.soiltype.perc_base,
		soil.soiltype.perc_exp,soil.soiltype.awc,max(1.0-fpc_phen_total,0.0),
		soil.awcont,soil.wcont,soil.wcont_evap,soil.runoff,soil.snowpack);

}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709

// guess2008 - new references:

// BONDEAU, A., SMITH, P. C., ZAEHLE, S., SCHAPHOFF, S., LUCHT, W., CRAMER, W., GERTEN, D., 
//   LOTZE-CAMPEN, H., MÜLLER, C., REICHSTEIN, M. and SMITH, B. (2007), 
//   Modelling the role of agriculture for the 20th century global terrestrial carbon balance. 
//   Global Change Biology, 13: 679–706. doi: 10.1111/j.1365-2486.2006.01305.x

// Rost, S., D. Gerten, A. Bondeau, W. Luncht, J. Rohwer, and S. Schaphoff (2008), 
//   Agricultural green and blue water consumption and its influence on the global 
//   water system, Water Resour. Res., 44, W09405, doi:10.1029/2007WR006331

