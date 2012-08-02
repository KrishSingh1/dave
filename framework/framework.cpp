///////////////////////////////////////////////////////////////////////////////////////
/// \file framework.cpp
/// \brief Implementation of the framework() function
///
/// \author Ben Smith
/// $Date: 2012-01-24 11:33:51 +0100 (Tue, 24 Jan 2012) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "framework.h"

#include "guessio.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"
#include "landcover.h"
#include "bvoc.h"


int framework(int argc,char* argv[]) {

	bool dogridcell;

	// The one and only linked list of Pft objects	
	Pftlist pftlist;

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(argc,argv,pftlist);

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

				// Calculate daylength, insolation and potential evapotranspiration
				daylengthinsoleet(gridcell.climate);

				if(run_landcover && date.day==0) {
					// Update dynamic landcover and crop fraction data during historical period and create/kill stands.
					if(date.year>=nyear_spinup)
						landcover_dynamics(gridcell,pftlist);
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
						// Leaf phenology for PFTs and individuals
						leaf_phenology(patch,gridcell.climate);
						// Interception
						interception(patch, gridcell.climate);
						initial_infiltration(patch, gridcell.climate);
						// Photosynthesis, respiration, evapotranspiration
						canopy_exchange(patch, gridcell.climate);
						// Soil water accounting, snow pack accounting
						soilwater(patch, gridcell.climate);
						// Soil organic matter and litter dynamics						
						som_dynamics(patch,pftlist);

						if (date.islastday && date.islastmonth) {

							// LAST DAY OF YEAR
							// Tissue turnover, allocation to new biomass and reproduction,
							// updated allometry
							growth(stand,patch);
						}
						stand.nextobj();
					}// End of loop through patches

					if (date.islastday && date.islastmonth) {
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
		
	}		// End of loop through grid cells

	// Call to input/output module to perform any necessary clean up
	termio();

	// END OF SIMULATION

	return 0;
}
