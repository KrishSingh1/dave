///////////////////////////////////////////////////////////////////////////////////////
/// \file ncinput.h
/// \brief Input module for NetCDF files with subdaily data.
///
/// \author Drew Holzworth
/// $Date: 2022-11-15 11:14:08 +1100 (Tue, 15 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_NCINPUT_H
#define LPJ_GUESS_NCINPUT_H

#ifdef HAVE_NETCDF

#include <memory>
#include <limits>
#include <stdint.h>
#include <string>
#include <vector>

#include "externalinput.h"
#include "guessnc.h"
#include "inputmodule.h"
#include "lamarquendep.h"
#include "ncutil.h"
#include "ncvariable.h"
#include "soilinput.h"

/*
NC Input: An input module which reads gridded met forcing data from a single
NetCDF file.

Requirements:

- A single .nc file containing met forcing data
- Requires the following variables in the NetCDF file:
  - Mean Temperature (°C)
  - Solar Insolation (W.m^-2)
  - Precipitation (mm)
  - CO2 (ppm)
  - Maximum temperature (°C)
  - Minimum temperature (°C)

Input data can be daily or subdaily. The following options must be specified in
the instruction file (run guess -help for descriptions, otherwise consult the
description in the code below):
  - nsteps_day
  - file_met_forcing
  - file_gridlist
  - file_soildata
  - file_ndep
  - ndep_scenario
  - has_leap_years
  - cache_strategy
  - variable_lon
  - variable_lat
  - variable_time
  - variable_temp
  - variable_insol
  - variable_precip
  - variable_co2
  - variable_tmax
  - variable_tmin
*/
class NCInput : public InputModule {
public:
	/*
	NCInput Constructor. This is called to create an NCInput instance.
	*/
	NCInput();

	/*
	Destructor.
	*/
	~NCInput();

	/*
	Initialise this input module by checking settings, opening input files, etc.
	This is called once at the start of the simulation.
	*/
	void init(void);

	/*
	Called once per gridcell. Obtains coordinates and soil static parameters for
	the gridcell before it is run. See base class documentation for more info.
	*/
	bool getgridcell(Gridcell&);

	/*
	Called once per day to load met forcing data for this day and store it in
	the gridcell's climate object. See base class documentation for more info.

	@return: True if there is another day to simulate. False otherwise.
	*/
	bool getclimate(Gridcell&);

	/*
	Obtain land transitions for this gridcell. See base class documentation for
	more info.
	*/
	void getlandcover(Gridcell&);

	/*
	Obtain land management data for one day. Called once per day. See base class
	documentation for more info.
	*/
	void getmanagement(Gridcell&);

private:

	/*
	Structure used to store coordinate information.
	*/
	struct Coord {
		/*
		Longitude of this coordinate.
		*/
		double lon;

		/*
		Latitude of this coordinate.
		*/
		double lat;

		/*
		Optional description of this coordinate, if one is provided.
		*/
		xtring description;
	};

	/*
	Create a new NCVariable object on the heap for the specified variable name.

	@param variable_name: Name of the variable.
	@param limits: Upper/lower bounds of the variable.
	@param units: Expected units of the variable.
	@param require_valid_units: Iff true, an error will be thrown upon
								encountering invalid units.
	*/
	std::unique_ptr<NCVariable> create_variable(xtring variable_name
		, bounds limits, const std::vector<xtring>& units
		, bool require_valid_units) const;

	/*
	Parse grid list from the specified plain text file.

	@param filename: Path to the gridlist file.
	*/
	std::vector<Coord> parse_gridlist(xtring filename) const;

	/*
	Calculate daily values of ndep for the current year and store them in dndep.

	@param lon: Longitude of the current gridcell.
	@param lat: Latitude of the current gridcell.
	@param year: Year of data to be used.
	*/
	void update_ndep(double lon, double lat, int year);

	/*
	Get the number of complete years of data which exist in the input file.
	*/
	int get_nyear_data() const;

	/*
	Get the progress of the current simulation in range [0, 1].

	@param since_spinup: If true, only consider progress since the restart.
	*/
	double get_progress(bool since_restart);

	/*
	Land cover input module.
	*/
	LandcoverInput landcover_input;

	/*
	Management input module.
	*/
	ManagementInput management_input;

	/*
	Soil input module.
	*/
	SoilInput soil_input;

	/*
	Nitrogen deposition input module.
	*/
	Lamarque::NDepData ndep_input;

	/*
	Grid cells to be simulated.
	*/
	std::vector<Coord> gridlist;

	/*
	The current grid cell.
	*/
	std::vector<Coord>::iterator current_gridcell;

	/*
	Number of time steps per day. This is a user-configurable setting in the
	instruction file.
	*/
	int32_t nsteps_day;

	/*
	Path to the NetCDF file. This is a user-configurable setting in the
	instruction file.
	*/
	xtring file_met_forcing;

	/*
	Path to gridlist file. This is a user-configurable setting in the
	instruction file.
	*/
	xtring file_gridlist;

	/*
	Path to the soil data file. This is a user-configurable setting in the
	instruction file.
	*/
	xtring file_soildata;

	/*
	Path to the nitrogen deposition fastarchive file. This is a
	user-configurable setting in the instruction file.
	*/
	xtring file_ndep;

	/*
	Scenario to be used for nitrogen deposition data. This is a
	user-configurable setting in the instruction file. Valid values are
	HISTORIC, RCP26, RCP45, RCP60, RCP85, and RCP85.
	*/
	xtring ndep_scenario;

	/*
	True iff the input data contains leap days. This is a user-configurable
	setting in the instruction file.
	*/
	bool has_leap_years;

	/*
	Caching strategy to be used when reading NetCDF data. This is a
	user-configurable setting in the instruction file.
	*/
	int cache_strategy;

	/*
	If true, met inputs will be checked for sensibility. This is a
	user-configurable setting in the instruction file.
	*/
	bool nc_check_inputs;

	/*
	Seed used for random sampling of met data during spinup period.
	*/
	int met_sampling_seed;

	/*
	Name of the longitude variable in the NetCDF file.
	*/
	xtring variable_lon;

	/*
	Name of the latitude variable in the NetCDF file.
	*/
	xtring variable_lat;

	/*
	Name of the time variable in the NetCDF file.
	*/
	xtring variable_time;

	/*
	Name of the mean temperature variable in the NetCDF file.
	*/
	xtring variable_temp;

	/*
	Name of the insolation variable in the NetCDF file.
	*/
	xtring variable_insol;

	/*
	Name of the precipitation variable in the NetCDF file.
	*/
	xtring variable_precip;

	/*
	Name of the co2 variable in the NetCDF file.
	*/
	xtring variable_co2;

	/*
	Name of the max temperature variable in the NetCDF file.
	*/
	xtring variable_tmax;

	/*
	Name of the min temperature variable in the NetCDF file.
	*/
	xtring variable_tmin;

	/*
	Name of the VPD variable in the NetCDF file.
	*/
	xtring variable_vpd;

	/*
	Name of the atmospheric pressure variable in the NetCDF file.
	*/
	xtring variable_ps;

	/*
	If true, each day's CO2 value will be the annual mean co2. Otherwise, a mean
	of that day's subdaily data will be used.
	*/
	bool use_mean_annual_co2;

	/*
	Temp variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> temp;

	/*
	Insolation variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> insol;

	/*
	Precipitation variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> precip;

	/*
	CO2 variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> co2;

	/*
	Tmax variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> tmax;

	/*
	Tmin variable in the NetCDF file.
	*/
	std::unique_ptr<NCVariable> tmin;

	/*
	VPD variable in the NetCDF file (kPa).
	*/
	std::unique_ptr<NCVariable> vpd;

	/*
	Atmospheric pressure variable in the NetCDF file (Pa).
	*/
	std::unique_ptr<NCVariable> ps;

	/*
	Number of seconds per day. Leap seconds are a myth.
	*/
	const int SECONDS_PER_DAY = 60 * 60 * 24;

	/*
	The NetCDF file ID.
	*/
	int nc_file;

	/*
	The first timestep in the input data.
	*/
	timepoint first_timestep;

	/*
	The last timestep in the input data. This is used to terminate the
	simulation when we reach the end date.
	*/
	timepoint final_timestep;

	/*
	Daily N deposition data for the current year.
	*/
	double dndep[Date::MAX_YEAR_LENGTH];

	/*
	Upper and lower bounds for prec.
	*/
	bounds bounds_prec;

	/*
	Upper and lower bounds for temp.
	*/	
	bounds bounds_temp;

	/*
	Upper and lower bounds for insol.
	*/	
	bounds bounds_insol;

	/*
	Upper and lower bounds for co2.
	*/	
	bounds bounds_co2;

	/*
	Upper and lower bounds for VPD.
	*/
	bounds bounds_vpd;

	/*
	Upper and lower bounds for atmospheric pressure.
	*/
	bounds bounds_ps;

	/*
	The year of data currently being sampled in the spinup. This will have a
	value of -1 when not in the spinup period.
	*/
	int spinup_year;

	/*
	Timer used to track simulation progress/runtime.
	*/
	Timer progress_timer;

	/*
	Timer used to display a progress message every MUTESEC seconds. (In other 
	words, it tracks the progress toward showing the next progress message - it
	tracks the progress progress.)
	*/
	Timer progress_progress_timer;

	/*
	Minimum number of sec to wait between progress messages
	*/
	static const int PROGRESS_INTERVAL=20;

	/*
	Total walltime spent in ncinput (cumulatively) (ms).
	*/
	double ncinput_walltime_use;

	/*
	Annual mean CO2 concentration (ppm) for the current simulation year.
	*/
	double mean_annual_co2;
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCINPUT_H
