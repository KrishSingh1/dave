///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK SOURCE CODE FILE
//
// Framework:             LPJ-GUESS Combined Modular Framework
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
// Header file name:      guess.h
// Source code file name: guess.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
// Updated:               2010-11-22


#include "config.h"
#include "guess.h"

#include "guessio.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"
#include "landcover.h"
#include "bvoc.h"

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
bool ifspeciesspecificwateruptake;	// water uptake is species specific 
// bvoc
bool ifbvoc; // BVOC calculations included

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
		npatchL=npatch;
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

Individual::Individual(int i,Pft& p,Vegetation& v):id(i),pft(p),vegetation(v) {

	anpp=0.0;
	fpc=0.0;
	densindiv=0.0;
	cmass_leaf=0.0;
	cmass_root=0.0;
	cmass_sap=0.0;
	cmass_heart=0.0;
	cmass_debt=0.0;
	wscal=1.0;
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
	ifwstress = false;
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

///////////////////////////////////////////////////////////////////////////////////////
// THE FRAMEWORK
// The 'mission control' of the model, responsible for maintaining the primary model
// data structures and containing all explicit loops through space (grid cells/stands)
// and time (days and years).

int framework(int argc,char* argv[]) {

	bool dogridcell;

#if defined DYNAMIC_LANDCOVER_INPUT
remove("LUdata.old");
rename("LUdata.out", "LUdata.old");
remove("CFTdata.old");
rename("CFTdata.out", "CFTdata.old");
#endif

	// The one and only linked list of Pft objects	
	Pftlist pftlist;

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(argc,argv,pftlist);

#ifdef MATS_TEST
	dprintf("\n");

	pftlist.firstobj();
	while(pftlist.isobj)
	{
		dprintf("pft n:o %d: %s",pftlist.getobj().id,(char*)pftlist.getobj().name);
		if(pftlist.getobj().landcover==CROPLAND)
			dprintf(", cftid %d",pftlist.getobj().cftid);			
		dprintf("\n");
		pftlist.nextobj();
	}

	if(run_landcover)
	{
		dprintf("\nLandcover version.\n");
		if(run[URBAN])
			dprintf("Urban stand simulated (landcover type %d)\n", URBAN);
		if(run[CROPLAND])
			dprintf("Crop stands simulated (landcover type %d)\n", CROPLAND);
		if(run[PASTURE])
			dprintf("Pasture stand simulated (landcover type %d)\n", PASTURE);
		if(run[FOREST])
			dprintf("Forest stand simulated (landcover type %d)\n", FOREST);
		if(run[NATURAL])
			dprintf("Natural stand simulated (landcover type %d)\n", NATURAL);
		if(run[PEATLAND])
			dprintf("Peatland stand simulated (landcover type %d)\n", PEATLAND);
	}
	dprintf("npft=%d\n",npft);
	if(run_landcover)
	{
		dprintf("ncft=%d\n",ncft);
//		if(run[CROPLAND] && minimizecftlist && !cftfrac_fixed)
//			dprintf("Pftlist is minimized\n");
	}

	dprintf("\n");
#endif


	// bvoc
	if(ifbvoc){
	  initbvoc(pftlist);
	}


	// Assume there is at least one grid cell to simulate
	dogridcell=true;

	while (dogridcell) {

		// START OF LOOP THROUGH GRID CELLS

		// Initialise global variable date
		// (argument nyear not used in this implementation)
		date.init(1);

		// Create and initialise a new Gridcell object for each locality
		Gridcell gridcell(pftlist);	

		// Call input/output to obtain latitude and soil driver data for this grid cell.
		// Function getgridcell returns false if no further grid cells remain to be simulated

		if (getgridcell(gridcell)) {

			// Initialise certain climate and soil drivers
			gridcell.climate.initdrivers(gridcell.climate.lat);

			if(run_landcover) {
				//Read static landcover and cft fraction data from ins-file and/or from data files for the spinup peroid and create stands.
				landcover_init(gridcell,pftlist);
			}
			
			// Call input/output to obtain climate, insolation and CO2 for this
			// day of the simulation. Function getclimate returns false if last year
			// has already been simulated for this grid cell

			while (getclimate(gridcell)) {

				// START OF LOOP THROUGH SIMULATION DAYS

				// Update daily climate drivers etc
				dailyaccounting_gridcell(gridcell,pftlist);

				if (run_landcover && run[CROPLAND])
					crop_sowing_gridcell(gridcell,pftlist);

				// Calculate daylength, insolation and potential evapotranspiration
				daylengthinsoleet(gridcell.climate);

				if(run_landcover && date.day==0) {
					// Update dynamic landcover and crop fraction data during historical period and create/kill stands.
					if(date.year>=nyear_spinup)
						landcover_dynamics(gridcell,pftlist);
/*
					if(run[CROPLAND] && forcesowingdates)
						getsowingdates(gridcell,pftlist);
					if(run[CROPLAND] && forceharvestdates)
						getharvestdates(gridcell,pftlist);
*/
				}

				gridcell.firstobj();
				while (gridcell.isobj) {

					// START OF LOOP THROUGH STANDS

					Stand& stand=gridcell.getobj();

					dailyaccounting_stand(stand,pftlist);

					stand.firstobj();
					while (stand.isobj) {
						// START OF LOOP THROUGH PATCHES

						// Get reference to this patch
						Patch& patch=stand.getobj();
						// Update daily soil drivers including soil temperature
						dailyaccounting_patch(patch,pftlist);

						if(stand.landcover==CROPLAND)
							crop_sowing_patch(patch, pftlist);

						// Leaf phenology for PFTs and individuals
						leaf_phenology(patch,gridcell.climate);
						// Photosynthesis, respiration, evapotranspiration
						canopy_exchange(patch);
						// Soil water accounting, snow pack accounting
						soilwater(gridcell.climate,patch);
						// Soil organic matter and litter dynamics
						som_dynamics(patch);

						if (stand.landcover==CROPLAND)
							growth_crop_daily(patch);

						if (date.islastday && date.islastmonth) {

							// LAST DAY OF YEAR
							// Tissue turnover, allocation to new biomass and reproduction,
							// updated allometry
							growth(stand,patch);
						}
						stand.nextobj();
					}// End of loop through patches

					if (date.islastday && date.islastmonth)
					{
						// LAST DAY OF YEAR
						stand.firstobj();
						while (stand.isobj) {

							// For each patch ...
							Patch& patch=stand.getobj();
							// Establishment, mortality and disturbance by fire
							vegetation_dynamics(stand,patch,pftlist);
							stand.nextobj();
						}
					}

					gridcell.nextobj();			
				}	// End of loop through stands

				if (date.islastday && date.islastmonth) {
					// LAST DAY OF YEAR
					// Call input/output module to output results for end of year
					// or end of simulation for this grid cell
					outannual(gridcell,pftlist);

					// Check whether to abort
					if (abort_request_received()) {
						termio();
						return 99;
					}
				}

				// Advance timer to next simulation day
				date.next();

				// End of loop through simulation days
			}//while (getclimate())
		}//if getgridcell()
		else dogridcell=false; // no more grid cells to simulate

		int test = 0;

		// End of loop through grid cells
	}

	// Call to input/output module to perform any necessary clean up
	termio();

	// END OF SIMULATION

	return 0;
}
