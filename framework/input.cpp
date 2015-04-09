////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file input.cpp
/// \brief Master class for all environmental input		
/// \author Mats Lindeskog
/// $Date: $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "input.h"
#include "driver.h"
#include "inputdefinitions.h"

// Switch to keep CO2 level at first historical year
const bool fixedco2_hist = 0;

////////////////////////////////////////////////////////////////////////////////
// Implementation of Input member functions
////////////////////////////////////////////////////////////////////////////////

Input::Input(const char* climate_input_module_name, const char* landcover_input_module_name) 
	:
	  ndep_input(*this),
	  soil_input(*this) 
	  {

	firsthistyear = -1;
	lasthistyear = -1;
	nyear_hist = 0;
	co2_fixed = 0;
	ngridcell = 0;
	gridlist_spatial_resolution = 0.5;

	climate_input_module = auto_ptr<InputModule>(InputModuleRegistry::get_instance().create_input_module(climate_input_module_name, *this));
	landcover_input_module = auto_ptr<LandcoverInputModule>(new LandcoverInputModule(landcover_input_module_name, *this));
	management_input_module = auto_ptr<ManagementInputModule>(new ManagementInputModule(*this));

	declare_parameter("firsthistyear", &firsthistyear, 1, 10000, "First historic year after spinup");
	declare_parameter("lasthistyear", &lasthistyear, 1, 10000, "Last historic year of simulation");
	declare_parameter("nyear_hist", &nyear_hist, 1, 10000, "Number of simulation years to run after spinup");
	declare_parameter("co2_fixed", &co2_fixed, 0.0, 2000.0,"Fixed CO2 value");

}

Input::~Input() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
}

void Input::read_gridlist() {

		// Reads list of grid cells and (optional) description text from grid list file
		// This file should consist of any number of one-line records in the format:
		//   <longitude> <latitude> [<description>]

		double dlon, dlat;
		bool eof = false;
		xtring descrip;

		// Read list of grid coordinates and store in Coord object 'gridlist'

		// Retrieve name of grid list file as read from ins file
		xtring file_gridlist=param["file_gridlist"].str;

		FILE* in_grid = fopen(file_gridlist,"r");
		if (!in_grid) fail("initio: could not open %s for input", (char*)file_gridlist);

		ngridcell = 0;
		while (!eof) {
			
			// Read next record in file
			eof =! readfor(in_grid, "f,f,a#", &dlon, &dlat, &descrip);

			if (!eof && !(dlon == 0.0 && dlat == 0.0)) { // ignore blank lines at end (if any)
				Coord& c = gridlist.createobj(); // add new coordinate to grid list

				c.lon = dlon;
				c.lat = dlat;
				c.descrip = descrip;
				ngridcell++;
			}
		}
		fclose(in_grid);
		gridlist.firstobj();

		// Parse spatial resolution
		double precision = 100;
		double dif_lon;
		double dif_lat;
		const int maxnsample = 200;
		int nsample = min((unsigned)maxnsample, gridlist.nobj);

		for (int i=0; i<nsample; i++) {
			for (int j=0; j<nsample; j++) {

				dif_lon = fabs(gridlist[i].lon - gridlist[j].lon);
				dif_lat = fabs(gridlist[i].lat - gridlist[j].lat);
				if(dif_lon > 1.0e-12)
					precision = min(precision, dif_lon);
				if(dif_lat > 1.0e-12)
					precision = min(precision, dif_lat);
			}
		}
		gridlist_spatial_resolution = precision;
}

void Input::init() {

	climate_input_module->init(); // Allow climate module to read gridlist

	if(!ngridcell)
		read_gridlist();

	if(co2_fixed == 0.0)
		// Read CO2 data from file if fixed CO2 not defined in instruction file
		co2.load_file(param["file_co2"].str);

	landcover_input_module->init();
	management_input_module->init();
#ifndef SOIL_INPUT_IN_CLIMATE_MODULE
	soil_input.init();
#endif
	// Setting of firsthistyear and nyear_hist; firsthistyear and lasthistyear are initiated to -1, nyear_hist to 0.

	// First look for firsthistyear and lasthistyear in instruction file
	if(firsthistyear > -1 && lasthistyear > -1) {
		nyear_hist = lasthistyear - firsthistyear + 1;
	}
	else if(firsthistyear < 0 || nyear_hist == 0) {

		// If firsthistyear or nyear_hist not defined in instruction file, use values in climate input module
		if(firsthistyear < 0)
			firsthistyear = climate_input_module->getfirsthistyear();
		if(nyear_hist == 0) {
			if(climate_input_module->getnyear_hist() > 0)
				nyear_hist = (climate_input_module->getfirsthistyear() + climate_input_module->getnyear_hist() - 1) - firsthistyear + 1;
			else
				nyear_hist = 0;
		}

		// If not already found, use values from landcover input module
		if(firsthistyear < 0)
			firsthistyear = landcover_input_module->getfirsthistyear();
		if(nyear_hist == 0) {
			if(landcover_input_module->getnyear_hist() > 0)
				nyear_hist = (landcover_input_module->getfirsthistyear() + landcover_input_module->getnyear_hist() - 1) - firsthistyear + 1;
			else
				nyear_hist = 0;
		}

		lasthistyear = firsthistyear + nyear_hist - 1;	// Not used further

		// firsthistyear and nyear_hist must be defined by now
		if(firsthistyear < 0 || nyear_hist == 0)
			fail("firsthistyear or nyear_hist not defined\n");

	}
//	date.set_first_calendar_year(firsthistyear - nyear_spinup);	// Will print historical years in output files

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

bool Input::getgridcell(Gridcell& gridcell) {

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	static bool first_call = true;

	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		bool gridfound = false;
		bool LUerror = false;

		while(!gridfound) {

			if(gridlist.isobj) {

				gridfound = climate_input_module->getgridcell(gridcell);
#ifndef SOIL_INPUT_IN_CLIMATE_MODULE
				if(gridfound)
					gridfound = soil_input.getgridcell(gridcell);
#endif				
				if (run_landcover && gridfound) {
					LUerror = landcover_input_module->getgridcell(gridcell);
					if(!LUerror)
						LUerror = management_input_module->getgridcell(gridcell);
				}

				if(!gridfound || LUerror) {
					if(!gridfound)
						dprintf("\nError: could not find stand at (%g,%g) in climate data files\n", gridlist.getobj().lon,gridlist.getobj().lat);
					else if(LUerror)
						dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", gridlist.getobj().lon,gridlist.getobj().lat);
					gridfound = false;
					gridlist.nextobj();
				}
			}
			else return false;
		}
		double lon = gridlist.getobj().lon;
		double lat = gridlist.getobj().lat;

		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Load nitrogen deposition data
#ifndef NDEP_INPUT_IN_CLIMATE_MODULE
		ndep_input.getgridcell(gridcell);
#endif
#ifndef SOIL_INPUT_IN_CLIMATE_MODULE
		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype, soil_input.getsoilcode());
#endif
		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

bool Input::getclimate(Gridcell& gridcell) {

	if (date.day == 0) {

		// Return false if last year was the last for the simulation
		if (date.year == nyear_spinup + nyear_hist)
			return false;

		// Progress report to user and update timer

		double progress;
		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id * (nyear_spinup + nyear_hist)
				+ date.year) / (double)(ngridcell * (nyear_spinup + nyear_hist));
			tprogress.setprogress(progress);
//			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress * 100.0),
			printf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress * 100.0),
				tprogress.elapsed.str, tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	// Send environmental values for today to framework
	climate_input_module->getclimate(gridcell);
	getco2(gridcell);
#ifndef NDEP_INPUT_IN_CLIMATE_MODULE
	getndep(gridcell);
#endif

	return true;
}

bool Input::getsoil(Gridcell& gridcell, const int soilmap_index) {

	bool error = soil_input.getgridcell(gridcell);

	return error;
}

double Input::getco2(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;

	if(co2_fixed > 0.0) {
		climate.co2 = co2_fixed;
	}
	else {
		
		if(fixedco2_hist)
			climate.co2 = co2[firsthistyear];
		else
			climate.co2 = co2[firsthistyear + date.year - nyear_spinup];
	}

	return climate.co2;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of NdepInput member functions
////////////////////////////////////////////////////////////////////////////////

NdepInput::NdepInput(Input& in)
	: input(in), 
	  gridlist(in.gridlist) {

	  ndep_fixed = 0;
	  declare_parameter("ndep_fixed", &ndep_fixed, 0.0, 2000.0,"Fixed ndep value");
}

bool NdepInput::getgridcell(Gridcell& gridcell) {

	double offset = 0.0;
	offset = search_for_centre_of_gridcell * (input.gridlist_spatial_resolution - Lamarque::SPATIAL_RESOLUTION) / 2.0;
	double lon = gridlist.getobj().lon + offset;
	double lat = gridlist.getobj().lat + offset;

	ndep.getndep(param["file_ndep"].str, lon, lat);

	return true;	// No info about success from NDepData
}

double NdepInput::getndep(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;

	if(ndep_fixed) {
		climate.dndep = ndep_fixed / (365.0 * 10000.0);	// Compatibility with the old demo input module. No distribution by prec.
	}
	else {

		if(date.day == 0) {

			int calender_year = date.year - nyear_spinup + input.getfirsthistyear();

			// Extract N deposition to use for this year,
			// monthly means to be distributed into daily values further down
			double mndrydep[12], mnwetdep[12];
			ndep.get_one_calendar_year(calender_year, mndrydep, mnwetdep);

			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, input.get_climate_module()->getdprec(), dndep);
		}

		climate.dndep = dndep[date.day];
	}

	return climate.dndep;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of SoilInput member functions
////////////////////////////////////////////////////////////////////////////////

SoilInput::SoilInput(Input& in)
	: input(in), 
	  gridlist(in.gridlist) {
}

void SoilInput::init() {

	file_soilcode = param["file_soilcode"].str;
	if(!soilcode.Open(file_soilcode, gridlist))
		fail("initio: could not open %s for input", (char*)file_soilcode);
}

bool SoilInput::loadsoilcode(Gridcell& gridcell, Coord c) {

	bool gridfound = false;
	double offset = search_for_centre_of_gridcell * input.gridlist_spatial_resolution / 2.0;

	if(!soilcode.Load(c, offset)) {
		dprintf("Problems with soil code input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n", c.lon, c.lat);
		gridfound = false;	// skip this stand
	}
	else
		gridfound = true;

	return gridfound;
}

bool SoilInput::getgridcell(Gridcell& gridcell) {
	
	Coord& c = gridlist.getobj();
	bool gridfound = loadsoilcode(gridcell, c);
	// New soil data input here

	if(!gridfound)
		dprintf("\nError: could not find stand at (%g,%g) in soil input data files\n", gridlist.getobj().lon,gridlist.getobj().lat);

	return gridfound;
}

int SoilInput::getsoilcode() {

	return (int)soilcode.Get(0, "soilcode");
}

