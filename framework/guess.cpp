///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.cpp
/// \brief LPJ-GUESS Combined Modular Framework
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
vegmodetype vegmode; // vegetation mode (population, cohort or individual)
int npatch; // number of patches in each stand (should always be 1 in population mode); cropland stands always have 1 patch
double patcharea; // patch area (m2) (individual and cohort mode only)
bool ifdailynpp; // whether NPP calculations performed daily (alt: monthly)
bool ifdailydecomp;
	// whether soil decomposition calculations performed daily (alt: monthly)
bool ifbgestab; // whether background establishment enabled (individual, cohort mode)
bool ifsme;
	// whether spatial mass effect enabled for establishment (individual, cohort mode)
bool ifstochestab; // whether establishment stochastic (individual, cohort mode)
bool ifstochmort; // whether mortality stochastic (individual, cohort mode)
bool iffire; // whether fire enabled
bool ifdisturb;
	// whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
bool ifcalcsla; // whether SLA calculated from leaf longevity (alt: prescribed)
int estinterval; // establishment interval in cohort mode (years)
double distinterval;
	// generic patch-destroying disturbance interval (individual, cohort mode)
int npft; // number of possible PFTs
bool iffast;
bool ifcdebt;

/// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
bool ifcentury;
/// whether plant growth limited by available N	
bool ifnlim;
/// number of years to allow spinup without N limitation	
int freenyears;
/// fraction of N relocated by plants from roots and leaves
double nrelocfrac;
/// whether to allow N leaching	
bool ifleachn;
/// whether to allow individual fractional N uptake
bool ifindiv_fnuptake;
/// first term in N fixation eqn
double nfix_a;
/// second term in N fixation eqn
double nfix_b;
/// whether N deposition data available from a file	
bool ifndepdata;

// CMIP5
bool ifcmip5;
bool iflandusesimple;
bool iflandusechange;

// N budget check
double somfluxnerror;

// guess2008 - new inputs from the .ins file
bool ifsmoothgreffmort;				// smooth growth efficiency mortality
bool ifdroughtlimitedestab;			// whether establishment affected by growing season drought
bool ifrainonwetdaysonly;			// rain on wet days only (1, true), or a little every day (0, false); 
// bvoc
bool ifbvoc; // BVOC calculations included

wateruptaketype wateruptake;

bool run_landcover;
bool run[NLANDCOVERTYPES];
bool lcfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;				// If a slow harvested product pool is included in patchpft.
int nyear_spinup;		


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX,Pftlist& pftlist):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL;

	for(p=0; p<pftlist. nobj;p++) {
		pft.createobj(pftlist[p]);
	}


	if(landcover == CROPLAND || landcover == PASTURE || landcover == URBAN || landcover == PEATLAND) {
		npatchL=1;
	}
	else if(landcover == NATURAL || landcover == FOREST) {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	for (p=0; p<npatchL; p++) {
		createobj(*this,pftlist,gc.soiltype);
	}

	first_year = date.year;
}

double Stand::get_gridcell_fraction() const {
	return frac * gridcell.landcoverfrac[landcover];
}

double Stand::get_landcover_fraction() const {
	return frac;
}

void Stand::set_landcover_fraction(double fraction) {
	frac = fraction;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp = 0.0;
	fpc = 0.0;
	densindiv = 0.0;
	cmass_leaf = 0.0;
	cmass_root = 0.0;
	cmass_sap = 0.0;
	cmass_heart = 0.0;
	cmass_debt = 0.0;
	phen = 0.0;
	aphen = 0.0;
	deltafpc = 0.0;
	fpar_wstress = 0.0;
	assim = 0.0;
	resp = 0.0;
	assim_nowstress = 0.0;

	nmass_leaf = 0.0;
	nmass_root = 0.0;
	nmass_sap = 0.0;
	nmass_heart = 0.0;
	nmass_reserve = 0.0;

	nstore = 0.0;
	nuptake = 0.0;
	ndemand = 0.0;
	ndemand_no_nlim = 0.0;
	fnuptake = 1.0;
	n_reserve_uptake = 0.0;
	max_n_reserve = 0.0;
	raingreen_ndemand = 0.0;

	frac_agpp = 1.0;

	// additional initialisation
	age = 0.0;
	fpar = 0.0;
	aphen_raingreen = 0;
	demand = 0.0;
	supply = 0.0;
	intercep = 0.0;
	phen_mean = 0.0;
	temp_wstress = 0.0;
	par_wstress = 0.0;
	daylength_wstress = 0.0;
	co2_wstress = 0.0; 
	nday_wstress = 0; 
	wstress = false;
	lai = 0.0;
	lai_layer = 0.0;
	lai_indiv = 0.0;
	alive = false;

	int m;
	for (m=0; m<12; m++) {
		mnpp[m] = mlai[m] = mgpp[m] = mra[m] = 0.0;
	}

	// bvoc
	monstor = 0.;
	iso = 0.;
	mon = 0.;
	aiso = 0.;
	amon = 0.;
	fvocseas = 1.;
	dtr_wstress = 0.;
	eet_wstress = 0.;
	agdd5_wstress = 0.;
	rad_wstress = 0.;		

	int d;
	for (d=0; d<365; d++) {
		dassim[d] = vmax_lim[d] = 0.0;
	}
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcell member functions
////////////////////////////////////////////////////////////////////////////////

double Gridcell::get_lon() const {
	return lon;
}

double Gridcell::get_lat() const {
	return lat;
}

void Gridcell::set_coordinates(double longitude, double latitude) {
	lon = longitude;
	lat = latitude;
}
