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
#include "bvoc.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
vegmodetype vegmode; // vegetation mode (population, cohort or individual)
int npatch; // number of patches in each stand (should always be 1 in population mode)
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


///////////////////////////////////////////////////////////////////////////////////////
// THE FRAMEWORK
// The 'mission control' of the model, responsible for maintaining the primary model
// data structures and containing all explicit loops through space (grid cells/stands)
// and time (days and years).

int framework(int argc,char* argv[]) {

	bool dostand;
	int p;

	// The one and only linked list of Pft objects	
	Pftlist pftlist;

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(argc,argv,pftlist);

	// bvoc
	if(ifbvoc){
	  initbvoc(pftlist);
	}

	// Assume there is at least one stand to simulate
	dostand=true;

	while (dostand) {

		// START OF LOOP THROUGH STANDS

		// Create and initialise a new Stand object for each grid cell / locality
		Stand stand(pftlist);

		// Call input/output to obtain latitude and soil driver data for this stand.
		// Function getstand returns false if no further stands remain to be simulated

		if (getstand(stand)) {

			// Initialise certain climate and soil drivers
			stand.climate.initdrivers(stand.climate.lat);
			initsoildrivers(stand);

			// Initialise global variable date
			// (argument nyear not used in this implementation)
			date.init(1);

			// Call input/output to obtain climate, insolation and CO2 for this
			// day of the simulation. Function getclimate returns false if last year
			// has already been simulated for this stand


			while (getclimate(stand)) {

				// START OF LOOP THROUGH SIMULATION DAYS

				// Update daily climate drivers etc
				dailyaccounting_stand(stand,pftlist);

				// Calculate daylength, insolation and potential evapotranspiration
				daylengthinsoleet(stand.climate);

				for (p=0;p<npatch;p++) {

					// START OF LOOP THROUGH PATCHES

					// Get reference to this patch
					Patch& patch=stand[p];

					// Update daily soil drivers including soil temperature
					dailyaccounting_patch(patch);

					// Leaf phenology for PFTs and individuals
					leaf_phenology(patch,stand.climate);

					// Photosynthesis, respiration, evapotranspiration
					canopy_exchange(patch);

					// Soil water accounting, snow pack accounting
					soilwater(stand.climate,patch);

					// Soil organic matter and litter dynamics
					som_dynamics(patch);

					if (date.islastday && date.islastmonth) {

						// LAST DAY OF YEAR

						// Tissue turnover, allocation to new biomass and reproduction,
						// updated allometry
						growth(stand,patch);
					}

					// End of loop through patches
				}


				if (date.islastday && date.islastmonth) {

					// LAST DAY OF YEAR

					for (p=0;p<npatch;p++) {

						// For each patch ...
						Patch& patch=stand[p];

						// Establishment, mortality and disturbance by fire
						vegetation_dynamics(stand,patch,pftlist);
					}

					// Call input/output module to output results for end of year
					// or end of simulation for this stand
					outannual(stand,pftlist);

					// Check whether to abort
					if (abort_request_received()) {
						termio();
						return 99;
					}
				}

				// Advance timer to next simulation day
				date.next();

				// End of loop through simulation days
			}
		}
		else dostand=false; // no more stands to simulate

		int test = 0;

		// End of loop through stands
	}

	// Call to input/output module to perform any necessary clean up
	termio();

	// END OF SIMULATION

	return 0;
}
