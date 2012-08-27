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
bool cftfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;				// If a slow harvested product pool is included in patchpft.
bool ifintercropgrass;
int ncft=0; // number of CFTs in Pftlist, set in plib_callback()
int nyear_spinup;		

bool forcesowingdates;
bool forceharvestdates;

////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions

////////////////////////////////////////////////////////////////////////////////

//const cropphen_struct* Patchpft::get_cropphen() 
cropphen_struct* Patchpft::get_cropphen() 
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropphen;
}

cropphen_struct* Patchpft::set_cropphen()
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropphen;
}


//Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX,Pftlist& pftlist):id(i),pftid(-1),gridcell(gc),isirrigated(false),hasgrassintercrop(false),landcover(landcoverX),frac(1.0) {
Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX,Pftlist& pftlist):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}

#ifdef MATS_TEST
	dprintf("Stand N:o %d, (landcover:%d) created year %d.\n", id,landcover, ::date.year);
#endif

	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND) {
		npatchL=1;
	}
	else if(landcover==NATURAL || landcover==FOREST) {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	for (p=0;p<npatchL;p++) {
		createobj(*this,pftlist,gc.soiltype);
	}

	first_year=date.year;
	natural_frac_change=0.0;
	seed=12345678;

	pftid=-1;
	cftid=-1;
	isirrigated=false;
	hasgrassintercrop=false;
	gdd0_intercrop=0.0;
}

double Stand::get_gridcell_fraction() const {
	return frac;
}

double Stand::get_landcover_fraction() const {
	if(gridcell.landcoverfrac[landcover])
		return frac/gridcell.landcoverfrac[landcover];
	else
		return 0.0;
}

void Stand::set_gridcell_fraction(double fraction) {
	frac = fraction;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp=0.0;
	fpc=0.0;
	densindiv=0.0;
	cmass_leaf=0.0;
	cmass_root=0.0;
	cmass_sap=0.0;
	cmass_heart=0.0;
	cmass_debt=0.0;
	phen=0.0;
	aphen=0.0;
	deltafpc=0.0;
	fpar_wstress=0.0;
	assim=0.0;

	// guess2008 - additional initialisation
	age=0.0;
	fpar=0.0;
	aphen_raingreen=0;
	demand=0.0;
	supply=0.0;
	intercep=0.0;
	phen_mean=0.0;
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
	for (m=0;m<12;m++) {
		mnpp[m]=mlai[m]=mlai_max[m]=mgpp[m]=mra[m]=0.0;
	}

	// bvoc
	monstor=0.;
	iso=0.;
	mon=0.;
	aiso=0.;
	amon=0.;
	fvocseas=1.;
	dtr_wstress=0.;
	eet_wstress=0.;
	agdd5_wstress=0.;
	rad_wstress=0.;		

	dnpp=0.0;
	cropindiv=NULL;

	if(pft.landcover==CROPLAND)
	{
		cropindiv=new cropindiv_struct;
//		vegetation.patch.pft[pft.id].cropphen->est_year=cropindiv->est_year;
	}

#ifdef MATS_TEST
//	dprintf("Year %d: Individual in stand %d created:id=%d, pft=%s\n", ::date.year-nyear_spinup+1901,vegetation.patch.stand.id,id,(char*)pft.name);
#endif
}

Individual::~Individual()
{
#ifdef MATS_TEST
//	dprintf("Year %d: Individual  in stand %d destroyed:id=%d, pft=%s\n",::date.year-nyear_spinup+1901,vegetation.patch.stand.id,id,(char*)pft.name);
#endif
	if(cropindiv)
		delete cropindiv;
}

//const cropindiv_struct* Individual::get_cropindiv() 
cropindiv_struct* Individual::get_cropindiv() 
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropindiv;
}

cropindiv_struct* Individual::set_cropindiv()
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropindiv;
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
