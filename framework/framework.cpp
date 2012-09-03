///////////////////////////////////////////////////////////////////////////////////////
/// \file framework.cpp
/// \brief Implementation of the framework() function
///
/// \author Ben Smith
/// $Date$
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
#include "guessserializer.h"


int framework(int argc, char* argv[]) {

	// The 'mission control' of the model, responsible for maintaining the 
	// primary model data structures and containing all explicit loops through 
	// space (grid cells/stands) and time (days and years).

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(argc, argv);

	// bvoc
	if (ifbvoc) {
	  initbvoc();
	}

	GuessSerializer* serializer = 0;
    GuessDeserializer* deserializer = 0;

	if (save && restart)
		fail("Serialization: Can't save and restart at the same time");

    if (save) {
        serializer = new GuessSerializer(state_path, 0);
    }
    if (restart) {
        deserializer = new GuessDeserializer(state_path);
    } 

	while (true) {

		// START OF LOOP THROUGH GRID CELLS

		// Initialise global variable date
		// (argument nyear not used in this implementation)
		date.init(1);

		// Create and initialise a new Gridcell object for each locality
		Gridcell gridcell;	

		// Call input/output to obtain latitude and soil driver data for this grid cell.
		// Function getgridcell returns false if no further grid cells remain to be simulated

		if (!getgridcell(gridcell)) {
			break;
		}

		// Initialise certain climate and soil drivers
		gridcell.climate.initdrivers(gridcell.get_lat());

		if(run_landcover) {
			//Read static landcover and cft fraction data from ins-file and/or from data files for the spinup peroid and create stands.
			landcover_init(gridcell);
		}

		if (restart) {
			deserializer->deserialize_gridcell(gridcell);
			date.year = start_year;
		}

		// Call input/output to obtain climate, insolation and CO2 for this
		// day of the simulation. Function getclimate returns false if last year
		// has already been simulated for this grid cell

		while (getclimate(gridcell)) {

			// START OF LOOP THROUGH SIMULATION DAYS

			// Update daily climate drivers etc
			dailyaccounting_gridcell(gridcell);

			// Calculate daylength, insolation and potential evapotranspiration
			daylengthinsoleet(gridcell.climate);

			if(run_landcover && date.day == 0 && date.year >= nyear_spinup) {
				// Update dynamic landcover and crop fraction data during historical
				// period and create/kill stands.
				landcover_dynamics(gridcell);
			}

			gridcell.firstobj();
			while (gridcell.isobj) {

				// START OF LOOP THROUGH STANDS

				Stand& stand = gridcell.getobj();

				dailyaccounting_stand(stand);

				stand.firstobj();
				while (stand.isobj) {
					// START OF LOOP THROUGH PATCHES

					// Get reference to this patch
					Patch& patch = stand.getobj();
					// Update daily soil drivers including soil temperature
					dailyaccounting_patch(patch);
					// Leaf phenology for PFTs and individuals
					leaf_phenology(patch, gridcell.climate);
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
						growth(stand, patch);
					}
					stand.nextobj();
				}// End of loop through patches

				if (date.islastday && date.islastmonth) {

					// FACE plantation
					if (has_FACE_clim) {
						if (!ifduke) {
							stand.plantyear=nyear_spinup+106-28; // 33 2073 should be 2087=1988
							stand.distyear2=stand.plantyear-45;
						}
						else {
							stand.plantyear=nyear_spinup+106-25; // 2074 should be 2082 (i.e. 1983), but forest needs more time to grow 2074
							stand.distyear2=stand.plantyear-283;
						}
					}
					
					// LAST DAY OF YEAR
					stand.firstobj();
					while (stand.isobj) {

						// For each patch ...
						Patch& patch = stand.getobj();
						// Establishment, mortality and disturbance by fire
						vegetation_dynamics(stand, patch);
						stand.nextobj();
					}
				}

				gridcell.nextobj();			
			}	// End of loop through stands

			if (date.islastday && date.islastmonth) {
				// LAST DAY OF YEAR
				// Call input/output module to output results for end of year
				// or end of simulation for this grid cell
				outannual(gridcell);

				if (date.year == start_year-1 && save) {
					serializer->serialize_gridcell(gridcell);
				}

				// Check whether to abort
				if (abort_request_received()) {
					termio();
					return 99;
				}
			}

			// Advance timer to next simulation day
			date.next();

			// End of loop through simulation days
		}	//while (getclimate())
	}		// End of loop through grid cells

	delete serializer;
    delete deserializer; 

	// Call to input/output module to perform any necessary clean up
	termio();

	// END OF SIMULATION

	return 0;
}
