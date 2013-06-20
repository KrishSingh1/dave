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
#include "commandlinearguments.h"
#include "guessserializer.h"
#include "parallel.h"

#include "cruinputmodule.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"
#include "landcover.h"
#include "bvoc.h"
#include "commonoutput.h"

#include <memory>

int framework(const CommandLineArguments& args) {

	// The 'mission control' of the model, responsible for maintaining the 
	// primary model data structures and containing all explicit loops through 
	// space (grid cells/stands) and time (days and years).

	using std::auto_ptr;

	const char* input_module_name;

#ifdef USE_CRU_IO
	input_module_name = "cru";
#else
	input_module_name = "demo";
#endif

	auto_ptr<InputModule> input_module(InputModuleRegistry::get_instance().create_input_module(input_module_name));

	GuessOutput::OutputModuleContainer output_modules;
	output_modules.add(new GuessOutput::CommonOutput());

	// Read the instruction file to obtain PFT static parameters and
	// simulation settings
	read_instruction_file(args.get_instruction_file());

	// Initialise input/output
	input_module->init();
	output_modules.init();

	// Nitrogen limitation
	if (ifnlim && !ifcentury) {
		fail("\n\nIf nitrogen limitation is switched on then century soil module also needs to be switched on!");
	}

	// bvoc
	if (ifbvoc) {
	  initbvoc();
	}

	// Create objects for (de)serializing grid cells
	auto_ptr<GuessSerializer> serializer;
	auto_ptr<GuessDeserializer> deserializer;

	if (save_state) {
		serializer = auto_ptr<GuessSerializer>(new GuessSerializer(state_path, GuessParallel::get_rank()));
	}

	if (restart) {
		deserializer = auto_ptr<GuessDeserializer>(new GuessDeserializer(state_path));
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

		if (!input_module->getgridcell(gridcell)) {
			break;
		}

		// Initialise certain climate and soil drivers
		gridcell.climate.initdrivers(gridcell.get_lat());

		if(run_landcover) {
			//Read static landcover and cft fraction data from ins-file and/or from data files for the spinup peroid and create stands.
			landcover_init(gridcell, input_module.get());
		}

		if (restart) {
			// Get the whole grid cell from file...
			deserializer->deserialize_gridcell(gridcell);
			// ...and jump to the restart year
			date.year = state_year;
		}
		
		// Call input/output to obtain climate, insolation and CO2 for this
		// day of the simulation. Function getclimate returns false if last year
		// has already been simulated for this grid cell

		while (input_module->getclimate(gridcell)) {

			// START OF LOOP THROUGH SIMULATION DAYS
			
			// Update daily climate drivers etc
			dailyaccounting_gridcell(gridcell);

			// Calculate daylength, insolation and potential evapotranspiration
			daylengthinsoleet(gridcell.climate);

			if(run_landcover && date.day == 0 && date.year >= nyear_spinup) {
				// Update dynamic landcover and crop fraction data during historical
				// period and create/kill stands.
				landcover_dynamics(gridcell, input_module.get());
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
					som_dynamics(patch);

					if (date.islastday && date.islastmonth) {

						// LAST DAY OF YEAR
						// Tissue turnover, allocation to new biomass and reproduction,
						// updated allometry
						growth(stand, patch);
					}
					stand.nextobj();
				}// End of loop through patches

				if (date.islastday && date.islastmonth) {
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

			output_modules.outdaily(gridcell);

			if (date.islastday && date.islastmonth) {
				// LAST DAY OF YEAR
				// Call output module to output results for end of year
				// or end of simulation for this grid cell
				output_modules.outannual(gridcell);
				
				// Time to save state?
				if (date.year == state_year-1 && save_state) {
					serializer->serialize_gridcell(gridcell);
				}

				// Check whether to abort
				if (abort_request_received()) {
					return 99;
				}
			}

			// Advance timer to next simulation day
			date.next();

			// End of loop through simulation days
		}	//while (getclimate())
	}		// End of loop through grid cells

	// END OF SIMULATION

	return 0;
}
