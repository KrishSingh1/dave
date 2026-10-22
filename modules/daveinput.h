///////////////////////////////////////////////////////////////////////////////////////
/// \file ncinput.h
/// \brief Input module for NetCDF files with subdaily data.
///
/// \author Drew Holzworth
/// $Date: 2022-11-15 11:14:08 +1100 (Tue, 15 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_DAVEINPUT_H
#define LPJ_GUESS_DAVEINPUT_H

#include "driver.h"
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
#include "globalco2file.h"
#include "thome.h"
#include "ncinput.h"
#include "stopwatch.h"

/*
NC Input: An input module which reads gridded met forcing data from a single
NetCDF file.

Requirements:

- A single .nc file containing met forcing data
- Requires the following variables in the NetCDF file:
  - Mean Temperature (°C)
  - Solar Insolation (W.m^-2)
  - Precipitation (mm)
  - CO2 (ppm) (optional if file_co2 is provided instead)
  - Maximum temperature (°C) (optional if temperature data is subdaily or if not using blaze)
  - Minimum temperature (°C) (optional if temperature data is subdaily or if not using blaze)
  - Nitrogen Deposition (MgN.m^-2)

Input data can be daily or subdaily. The following options must be specified in
the instruction file (run guess -help for descriptions, otherwise consult the
description in the code below):
  - nsteps_day
  - file_met_forcing
  - file_gridlist
  - file_soildata
  - file_ndep
  - ndep_scenario
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
class DaveInput : public InputModule {
public:
	/*
	DaveInput Constructor. This is called to create an DaveInput instance.
	*/
	DaveInput();

	/*
	Destructor.
	*/
	~DaveInput();

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
	bool getlandcover(Gridcell&, double*);

	/*
	Obtain land management data for one day. Called once per day. See base class
	documentation for more info.
	*/
	void getmanagement(Gridcell&);

	/*
	Obtains additional environmental data that are not dynamic (e.g. elevation).
	*/
	void getmiscinput_static(Gridcell& gridcell);

	/*
	Obtains additional environmental data (e.g. disturbance) for one year.
	*/
	void getmiscinput_yearly(Gridcell& gridcell);

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
	Parse grid list from the specified plain text file.

	@param filename: Path to the gridlist file.
	*/
	std::vector<Coord> parse_gridlist(xtring filename) const;

	/*
	Calculate daily values of ndep for the current year and store them in dndep.

	@param year Year of data to be used.
	*/
	void update_ndep(int year);

	/*
	Get the progress of the current simulation in range [0, 1].

	@param since_spinup: If true, only consider progress since the restart.
	*/
	double get_progress(bool since_restart);

	/*
	Get an atmospheric co2 concentration for the specified location on the
	specified date.
	*/
	double get_co2();

	/*
	Get the year of data to be used for the current year of the spinup.
	*/
	int get_spinup_year() const;

	/*
	Resize the climate arrays in the gridcell to match the specified value.
	*/
	void resize_climate_arrays(Gridcell &gridcell, size_t size) const;

	/*
	Parse the specified esat_method name into an esat_method_t.
	*/
	esat_method_t parse_esat_method(xtring name) const;

	/*
	Land cover input module.
	*/
	LandcoverInput landcover_input;

	/*
	Management input module.
	*/
	ManagementInput management_input;

	/*
	Additional text data input module.
	*/
	MiscInput misc_input;

	/*
	Soil input module.
	*/
	SoilInput soil_input;

	/*
	THome input module (use for temperature acclimation of photosynthesis).
	*/
	THomeInput thome_input;

	/*
	Nitrogen deposition input module.
	*/
	Lamarque::NDepData ndep_input;

	/*
	True iff the Lamarque N deposition data source should be used.
	*/
	bool use_lamarque;

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
	Path to the NetCDF file used for the spinup. This is a user-configurable
	setting in the instruction file.
	*/
	xtring file_met_spinup;

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
	Optional path to text file containing co2 concentration by year. If not set,
	variable_co2 must be set, allowing co2 data to be read from the netcdf file.
	*/
	xtring file_co2;

	/*
	Path to the nitrogen deposition fastarchive file. This is a
	user-configurable setting in the instruction file.
	*/
	xtring file_ndep;

	/*
	Path to the mean maximum daily temperature (T_home) input file, used for
	DAVE temperature acclimation calculations.
	*/
	xtring file_thome;

	/*
	Scenario to be used for nitrogen deposition data. This is a
	user-configurable setting in the instruction file. Valid values are
	HISTORIC, RCP26, RCP45, RCP60, RCP85, and RCP85.
	*/
	xtring ndep_scenario;

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
	User input from the instruction file. The number of years of data to recycle
	during the spinup. If <= 0, all data will be recycled during freenyears.
	*/
	int spinup_met_cycle_length;

	/*
	Constant to be added to daily temperature (°C). Set in instruction file.
	*/
	double temperature_delta;

	/*
	Constant to be added to daily co2 (ppm). Set in instruction file.
	*/
	double co2_delta;

	/*
	Scalar to multiply daily precipitation by. Set in instruction file.
	*/
	double precip_scalar;

	/*
	Constant wind speed value to be used if it's unavailable in the input data
	(m s-1).
	*/
	double wind_speed;

	/*
	If true, each day's CO2 value will be the annual mean co2. Otherwise, a mean
	of that day's subdaily data will be used.
	*/
	bool use_mean_annual_co2;

	/*
	Optional co2 text file reader.
	*/
	std::unique_ptr<GlobalCO2File> co2_file;

	/*
	Climate input submodule.
	*/
	NCInput ncinput;

	/*
	Number of seconds per day. Leap seconds are a myth.
	*/
	const int SECONDS_PER_DAY = 60 * 60 * 24;

	/*
	Daily NH4 deposition data for the current year.
	*/
	double nh4_dep[Date::MAX_YEAR_LENGTH];

	/*
	Daily NO3 deposition data for the current year.
	*/
	double no3_dep[Date::MAX_YEAR_LENGTH];

	/*
	Upper and lower bounds for prec.
	*/
	bounds bounds_prec;

	/*
	Upper and lower bounds for wind speed (m/s).
	*/
	bounds bounds_wind;

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
	Upper and lower bounds for relative humidity.
	*/
	bounds bounds_rh;

	/*
	Upper and lower bounds for atmospheric pressure.
	*/
	bounds bounds_ps;

	/*
	Upper and lower bounds for nitrogen deposition data.
	*/
	bounds bounds_ndep;

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
	Minimum amount of time (in seconds) to wait between progress messages.
	*/
	int progress_interval;

	/*
	Force the use of newline characters to terminate progress messages. Useful
	for MPI runs. If false, \r would be used if stdout is a terminal, and \n
	would be used otherwise.
	*/
	bool progress_newline;

	/*
	If enabled, invalid mathematical operations (e.g. division by 0) will result
	in a fatal error.
	*/
	bool strict_math;

	/*
	Stopwatch used to track the total time spent in this input module.
	*/
	Stopwatch stopwatch;

	/*
	Fixed/constant air pressure to use, in Pa.
	*/
	double fixed_pressure;

	/*
	Name of the method to use for estimating saturation vapor pressure.
	*/
	xtring esat_method_name;

	/*
	Method to use for estimating saturation vapor pressure.
	*/
	esat_method_t esat_method;
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_DAVEINPUT_H
