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

#include "inputmodule.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"
#include "blaze.h"
#include "simfire.h"
#include "landcover.h"
#include "bvoc.h"
#include "commonoutput.h"
#include "soilmethane.h"

#include <memory>

#include <chrono>
using namespace std;
using namespace std::chrono;

/// Prints the date and time together with the name of this simulation
void print_logfile_heading() {
	xtring datetime;
	unixtime(datetime);

	xtring header = xtring("[LPJ-GUESS  ") + datetime + "]\n\n";
	dprintf((char*)header);

	// Print the title of this run
	std::string dashed_line(50, '-');
	dprintf("\n\n%s\n%s\n%s\n",
	        dashed_line.c_str(), (char*)title, dashed_line.c_str());
}

/// Simulate one day for a given Gridcell
/**
 * The climate object in the gridcell needs to be set up with
 * the day's forcing data before calling this function.
 *
 * \param gridcell            The gridcell to simulate
 * \param input_module        Used to get land cover fractions
 */
void simulate_day(Gridcell& gridcell, InputModule* input_module) {

	// Update daily climate drivers etc
	dailyaccounting_gridcell(gridcell);

	// Calculate daylength, insolation and potential evapotranspiration
	daylengthinsoleet(gridcell.climate);

	// Update crop sowing date calculation framework
	crop_sowing_gridcell(gridcell);

	// Dynamic landcover and crop fraction data during historical
	// period and create/kill stands.
	landcover_dynamics(gridcell, input_module);

	// Update dynamic management options
	input_module->getmanagement(gridcell);

	// Set forest management for all stands this year
	manage_forests(gridcell);

	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		// START OF LOOP THROUGH STANDS
		Stand& stand = *gc_itr;

		// Get starting timepoint
		auto time00 = high_resolution_clock::now();

		dailyaccounting_stand(stand);

		// Get starting timepoint
		auto time01 = high_resolution_clock::now();

		stand.firstobj();
		while (stand.isobj) {

			// START OF LOOP THROUGH PATCHES

			// Get reference to this patch
			Patch& patch = stand.getobj();

			// Get starting timepoint
			auto time0 = high_resolution_clock::now();

			// Update daily soil drivers including soil temperature
			dailyaccounting_patch(patch);

			// Get starting timepoint
			auto time1 = high_resolution_clock::now();

			// Determine nitrogen fertilisation amount 
			if(run_landcover)
				nfert(patch);

			// Calculate crop sowing dates
			crop_sowing_patch(patch);

			// Crop phenology
			crop_phenology(patch);

			// Get starting timepoint
			auto time2 = high_resolution_clock::now();

			// Leaf phenology for PFTs and individuals
			leaf_phenology(patch, gridcell.climate);

			// Get starting timepoint
			auto time3 = high_resolution_clock::now();

			// Interception
			interception(patch, gridcell.climate);

			// Get starting timepoint
			auto time4 = high_resolution_clock::now();

			// Initial infiltration
			initial_infiltration(patch, gridcell.climate);

			// Get starting timepoint
			auto time5 = high_resolution_clock::now();

			// Photosynthesis, respiration, evapotranspiration
			canopy_exchange(patch, gridcell.climate);

			// Get starting timepoint
			auto time6 = high_resolution_clock::now();

			// Sum total required irrigation
			irrigation(patch);

			// Get starting timepoint
			auto time7 = high_resolution_clock::now();
			
			// Soil water and snow pack accounting
			soilwater(patch, gridcell.climate);

			// Get starting timepoint
			auto time8 = high_resolution_clock::now();

			// Daily C allocation (cropland)
			growth_daily(patch);

			// Get starting timepoint
			auto time9 = high_resolution_clock::now();

			// Soil organic matter and litter dynamics
			som_dynamics(patch, gridcell.climate);

			// Get starting timepoint
			auto time10 = high_resolution_clock::now();

			// Methane production/consumption on wetlands and peatlands (no methane dynamics for other stand types at present) 
			methane_dynamics(patch);

			// Get starting timepoint
			auto time11 = high_resolution_clock::now();

			// BLAZE fire model
			blaze_driver(patch,gridcell.climate);

			// Get starting timepoint
			auto time12 = high_resolution_clock::now();

			// Get duration. Substart timepoints to 
			// get durarion. To cast it to proper unit
			// use duration cast method
			auto duration0 = duration_cast<microseconds>(time01 - time00);
			auto duration1 = duration_cast<microseconds>(time1 - time0);
			auto duration2 = duration_cast<microseconds>(time2 - time1);
			auto duration3 = duration_cast<microseconds>(time3 - time2);
			auto duration4 = duration_cast<microseconds>(time4 - time3);
			auto duration5 = duration_cast<microseconds>(time5 - time4);
			auto duration6 = duration_cast<microseconds>(time6 - time5);
			auto duration7 = duration_cast<microseconds>(time7 - time6);
			auto duration8 = duration_cast<microseconds>(time8 - time7);
			auto duration9 = duration_cast<microseconds>(time9 - time8);
			auto duration10 = duration_cast<microseconds>(time10 - time9);
			auto duration11 = duration_cast<microseconds>(time11 - time10);
			auto duration12 = duration_cast<microseconds>(time12 - time11);

			dprintf("dailyac_s %d dailyac_p %d crop %d phen %d inter %d infil %d canex %d irri %d soilw %d dgrowth %d somdyn %d methane %d blaze %d\n",
				duration0.count(), duration1.count(), duration2.count(), duration3.count(), duration4.count(), duration5.count(), duration6.count(),
				duration7.count(), duration8.count(), duration9.count(), duration10.count(), duration11.count(), duration12.count());

			/*dprintf("Year %d Day %d dailyac_s %d dailyac_p %d crop %d phen %d inter %d infil %d canex %d irri %d soilw %d dgrowth %d somdyn %d methane %d blaze %d\n",
				date.year, date.day, duration0.count(), duration1.count(), duration2.count(), duration3.count(), duration4.count(), duration5.count(), duration6.count(),
				duration7.count(), duration8.count(), duration9.count(), duration10.count(), duration11.count(), duration12.count());*/

			if (date.islastday && date.islastmonth) {

				// LAST DAY OF YEAR
				// Tissue turnover, allocation to new biomass and reproduction,
				// updated allometry
				growth(stand, patch);
			}
			stand.nextobj();
		}// End of loop through patches

		// Update crop rotation status
		crop_rotation(stand);

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

		++gc_itr;
	}	// End of loop through stands
}


int framework(const CommandLineArguments& args) {

	// The 'mission control' of the model, responsible for maintaining the
	// primary model data structures and containing all explicit loops through
	// space (grid cells/stands) and time (days and years).

	using std::auto_ptr;

	const char* input_module_name = args.get_input_module();

	auto_ptr<InputModule> input_module(InputModuleRegistry::get_instance().create_input_module(input_module_name));

	GuessOutput::OutputModuleContainer output_modules;
	GuessOutput::OutputModuleRegistry::get_instance().create_all_modules(output_modules);

	// Read the instruction file to obtain PFT static parameters and
	// simulation settings
	read_instruction_file(args.get_instruction_file());

	// Initialise input/output

	input_module->init();
	output_modules.init();

	print_logfile_heading();

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
		serializer = auto_ptr<GuessSerializer>(new GuessSerializer(state_path, GuessParallel::get_rank(), GuessParallel::get_num_processes()));
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

		// Call input module to obtain latitude and driver data for this grid cell.
		if (!input_module->getgridcell(gridcell)) {
			break;
		}

		// Initialise certain climate and soil drivers
		gridcell.climate.initdrivers(gridcell.get_lat());

		// Read landcover and cft fraction data from 
		// data files for the spinup period and create stands
		landcover_init(gridcell, input_module.get());

		// Set soil depth dependent properties
		gridcell.set_soil_depth();

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

			simulate_day(gridcell, input_module.get());

			output_modules.outdaily(gridcell);

			if (date.islastday && date.islastmonth) {
				// LAST DAY OF YEAR
				if(printseparatestands)
					output_modules.openlocalfiles(gridcell);
				// Call output module to output results for end of year
				// or end of simulation for this grid cell
				output_modules.outannual(gridcell);

				gridcell.balance.check_year(gridcell);

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

		if(printseparatestands)
			output_modules.closelocalfiles(gridcell);

		gridcell.balance.check_period(gridcell);

	}		// End of loop through grid cells

	// END OF SIMULATION

	return 0;
}
