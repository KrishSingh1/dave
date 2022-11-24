#ifdef HAVE_NETCDF

#include <exception>
#include <limits.h>
#include <netcdf>
#include <stdexcept>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "driver.h"
#include "guessnc.h"
#include "guessstring.h"
#include "ncinput.h"
#include "ncutil.h"
#include "parameters.h"

REGISTER_INPUT_MODULE("nc", NCInput)

NCInput::NCInput() {
	double ctor_start = get_wall_time();

	// Register instruction file parameters.

	// Setting max timesteps per day to the number of seconds in a day.
	// Technically we could support more; this limitation is quite arbitrary.
	declare_parameter("nsteps_day", &nsteps_day, 1, SECONDS_PER_DAY
	, "Number of subdaily timesteps");

	// Again, this maximum file name length is pretty arbitrary and not really a
	// hard limit.
	declare_parameter("file_met_forcing", &file_met_forcing, PATH_MAX
	, "Path to NetCDF file containing met forcing data");
	declare_parameter("file_gridlist", &file_gridlist, PATH_MAX
	, "Path to Gridlist file containing lines of the form \"lon lat\"");
	declare_parameter("file_soildata", &file_soildata, PATH_MAX
	, "Path to soil input file");
	declare_parameter("file_ndep", &file_ndep, PATH_MAX
	, "Path to nitrogen deposition fastarchive input file");

	declare_parameter("cache_strategy", &cache_strategy, 0, 2
	, "Cache strategy to be used when reading NetCDF data (0 = cache entire file (fastest), 1 = cache one year at a time (slower, but uses less memory), 2 = no caching (slowest, not recommended)");

	declare_parameter("ndep_scenario", &ndep_scenario, 256
	, "Scenario to be used for nitrogen deposition data. E.g. HISTORIC, RCP26, RCP45, RCP60, RCP85, RCP85");

	declare_parameter("has_leap_years", &has_leap_years
	, "1 iff input data contains leap years.");

	declare_parameter("met_sampling_seed", &met_sampling_seed, 0, __INT32_MAX__, "Seed used for random sampling of met data during spinup period");

	declare_parameter("nc_check_inputs", &nc_check_inputs, "1 to enable sensibility checking of met inputs, 0 otherwise");

	// Names of variables in the NetCDF file.
	declare_parameter("variable_lon", &variable_lon, NC_MAX_NAME
	, "Name of the longitude variable in the NetCDF file");
	declare_parameter("variable_lat", &variable_lat, NC_MAX_NAME
	, "Name of the latitude variable in the NetCDF file");
	declare_parameter("variable_time", &variable_time, NC_MAX_NAME
	, "Name of the time variable in the NetCDF file");
	declare_parameter("variable_temp", &variable_temp, NC_MAX_NAME
	, "Name of the mean temperature variable in the NetCDF file");
	declare_parameter("variable_insol", &variable_insol, NC_MAX_NAME
	, "Name of the insolation variable in the NetCDF file");
	declare_parameter("variable_precip", &variable_precip, NC_MAX_NAME
	, "Name of the precipitation variable in the NetCDF file");
	declare_parameter("variable_co2", &variable_co2, NC_MAX_NAME
	, "Name of the co2 variable in the NetCDF file");
	declare_parameter("variable_tmax", &variable_tmax, NC_MAX_NAME
	, "Name of the max temperature variable in the NetCDF file");
	declare_parameter("variable_tmin", &variable_tmin, NC_MAX_NAME
	, "Name of the min temperature variable in the NetCDF file");

	bounds_prec.max = 1000;
	bounds_prec.min = 0;

	bounds_insol.max = 1e6;
	bounds_insol.min = 0;

	bounds_co2.max = 900;
	bounds_co2.min = 0;

	bounds_temp.max = 100;
	bounds_temp.min = -100;

	double ctor_end = get_wall_time();
	double ctor_duration = (ctor_end - ctor_start) * 1000;
	ncinput_walltime_use = ctor_duration;
}

NCInput::~NCInput() {
	GuessNC::close_ncdf(nc_file);

	// The NCVariable instances are dynamically allocated, but they're wrapped
	// in smart pointers, so we don't need to manually free them.
}

void NCInput::init(void) {
	try {
		double init_start = get_wall_time();

		// Open the met forcing file (in readonly mode).
		nc_file = GuessNC::open_ncdf((const char *)file_met_forcing);

		// Parse the gridlist.
		gridlist = parse_gridlist(file_gridlist);

		// Open soil data file.
		soil_input.init(file_soildata);

		std::vector<xtring> temp_units = {
			"degC",
		};
		std::vector<xtring> insol_units = {
			"W/m2"
			// Should provide units for other insol types as well...
		};
		std::vector<xtring> precip_units = {
			"mm"
		};
		std::vector<xtring> co2_units = {
			"ppm"
		};

		// Create NCVariable objects for each variable of interest.
		// This will also do some basic sanity checking.
		temp = create_variable(variable_temp, bounds_temp, temp_units, true);
		insol = create_variable(variable_insol, bounds_insol, insol_units, true);
		precip = create_variable(variable_precip, bounds_prec, precip_units, true);
		co2 = create_variable(variable_co2, bounds_co2, co2_units, true);
		if (variable_tmax != "") {
			tmax = create_variable(variable_tmax, bounds_temp, temp_units, true);
		}
		if (variable_tmin != "") {
			tmin = create_variable(variable_tmin, bounds_temp, temp_units, true);
		}

		// Move to first gridcell in the gridlist.
		current_gridcell = gridlist.begin();

		// Record length of input data.
		first_timestep = temp->get_first_date();
		final_timestep = temp->get_last_date();

		// Initialise PRNG.
		srand(met_sampling_seed);

		int nyear_data = get_nyear_data();
		if (nyear_data < 1) {
			throw std::runtime_error("Insufficient met forcing data (<1 year)");
		}

	// Set timers.
	progress_timer.init();
	progress_progress_timer.init();

	progress_timer.settimer();
	progress_progress_timer.settimer(PROGRESS_INTERVAL);

	double init_end = get_wall_time();
	double init_duration = (init_end - init_start) * 1000;
	ncinput_walltime_use += init_duration;

	} catch (const std::exception& error) {
		// todo: do we need to manually unwrap nested exceptions?
		fail(get_nested_error_msg(error));
	}
}

bool NCInput::getgridcell(Gridcell& gridcell) {
	try {
		double getgridcell_start = get_wall_time();

		if (current_gridcell == gridlist.end()) {
			return false;
		}

		double lon = current_gridcell->lon;
		double lat = current_gridcell->lat;
		gridcell.set_coordinates(lon, lat);

		// We can get instype from insol's standard_name attribute.
		xtring insol_name = insol->get_standard_name();
		gridcell.climate.instype = parse_insoltype((const char*)insol_name);

		// Get soil data for this grid cell.
		soil_input.get_soil(lon, lat, gridcell);

		// Read N deposition data for this grid cell.
		ndep_input.getndep(file_ndep, lon, lat
		, Lamarque::parse_timeseries((const char *)ndep_scenario));

		dprintf("\nCommencing simulation for gridcell at (%g,%g)", lon, lat);
		if (current_gridcell->description != "") {
			dprintf(" (%s)", (const char *)current_gridcell->description);
		}
		dprintf("\n");

		// Increment gridlist iterator (it now points to the next grid cell).
		current_gridcell++;

		double getgridcell_end = get_wall_time();
		double getgridcell_duration = (getgridcell_end - getgridcell_start) * 1000;
		ncinput_walltime_use += getgridcell_duration;

		return true;
	} catch (const std::exception& e) {
		fail(get_nested_error_msg(e));
		return false; // To eliminate a compiler warning.
	}
}

bool NCInput::getclimate(Gridcell& gridcell) {
	try {
	double getclimate_start = get_wall_time();

	// Shorthand access to some variables.
	Climate& climate = gridcell.climate;
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();

	// First, check if we have enough data to simulate another day.
	if (year == final_timestep.year) {
		if (final_timestep.day < Date::MAX_YEAR_LENGTH - 1) {
			// We are in the final year, which is missing some days of data.
			// Because the N deposition disaggregation requires an entire years'
			// worth of precipitation data, we can't simulate part of a year.
			return false;
		}

		if (final_timestep.timestep < nsteps_day - 1) {
			// We are in the final year, and the final day is missing some
			// timesteps.
			return false;
		}
	}

	if (date.year < nyear_spinup) {
		// Shuffle years of met data during spinup. Whenever we move to a new
		// year, select another year at random.
		if (date.day == 0) {
			int first_year = first_timestep.year;
			if (!is_start_of_year(first_timestep, nsteps_day)) {
				first_year++;
			}
			spinup_year = rand() % get_nyear_data() + first_year;
		}
		year = spinup_year;
	} else {
		// We are now out of the spinup. Reset spinup_year to avoid confusion.
		spinup_year = -1;
	}

	// Generate daily N deposition values on first day of year.
	if (date.day == 0) {
		update_ndep(lon, lat, year);
	}

	// Read data from variables (not necessarily from the file - these variable
	// objects will handle caching internally).
	std::vector<double> temps = temp->read(lon, lat, year, date.day, 1);
	std::vector<double> precs = precip->read(lon, lat, year, date.day, 1);
	std::vector<double> insols = insol->read(lon, lat, year, date.day, 1);
	std::vector<double> co2s = co2->read(lon, lat, year, date.day, 1);

	// Aggregate subdaily values to daily where necessary, and store them in the
	// climate object. Note that sum inputs are rates and sum are quantities;
	// We can use the mean of the rates and the sum of the quantities. Check the
	// units for more info.
	climate.temp = mean(temps);
	climate.prec = sum(precs);
	climate.insol = mean(insols);
	climate.co2 = mean(co2s);
	climate.temps = temps;
	climate.insols = insols;
	climate.dndep = dndep[date.day];

	if (tmax && tmin) {
		std::vector<double> tmaxs = tmax->read(lon, lat, year, date.day, 1);
		std::vector<double> tmins = tmin->read(lon, lat, year, date.day, 1);
		climate.dtr = mean(tmaxs) - mean(tmins);
	} else if (temps.size() > 1) {
		climate.dtr = max(temps) - min(temps);
	} else {
		throw std::runtime_error("Missing tmax/tmin data, and day length is 1");
	}

	// rads and pars are calculated automatically in daylengthinsoleet().
	// gtemps are calculated in respiration_temperature_response().

	double getclimate_duration = (get_wall_time() - getclimate_start) * 1000;
	ncinput_walltime_use += getclimate_duration;

	// Progress reporting.
	if (date.day == 0 && progress_progress_timer.getprogress() >= 1.0) {

		// Number of gridcells which have been run.
		int cells_done = distance(gridlist.begin(), current_gridcell) - 1;

		// Number of years to simulate per gridcell.
		int years_to_simulate = nyear_spinup + get_nyear_data();

		// Total number of years to be simulated across all gridcells.
		int total_years = years_to_simulate * gridlist.size();

		// Number of years which have been simulated (date.year starts at 0).
		int years_simulated = years_to_simulate * cells_done + date.year;

		double progress = 1.0 * years_simulated / total_years;
		progress_timer.setprogress(progress);
		dprintf("%3d%% complete, %s elapsed, %s remaining\n"
			, (int)(progress * 100), progress_timer.elapsed.str
			, progress_timer.remaining.str);
		progress_progress_timer.settimer(PROGRESS_INTERVAL);
	}

	// We have enough data to simulate another day.
	return true;
	} catch (const std::exception& e) {
		fail(get_nested_error_msg(e));
		return false; // To eliminate a compiler warning.
	}
}

void NCInput::getlandcover(Gridcell& gridcell) {
	// Proxy this request through to the landcover input module.
	landcover_input.getlandcover(gridcell);
}

void NCInput::getmanagement(Gridcell& gridcell) {
	// Proxy this request through to the management input module.
	management_input.getmanagement(gridcell);
}

std::unique_ptr<NCVariable> NCInput::create_variable(xtring variable_name
	, bounds limits, const std::vector<xtring>& units
	, bool require_valid_units) const {
	return std::unique_ptr<NCVariable>(new NCVariable(
		nc_file,
		variable_name,
		variable_lon,
		variable_lat,
		variable_time,
		nsteps_day,
		has_leap_years,
		(CacheStrategy)cache_strategy,
		nc_check_inputs,
		limits,
		units,
		require_valid_units
	));
}

std::vector<NCInput::Coord> NCInput::parse_gridlist(xtring filename) const {
	// Open the file for reading.
	std::ifstream ifs(filename, std::ifstream::in);

	std::vector<Coord> coordinates;

	// Iterate over the file line by line.
	std::string line;
	while (getline(ifs, line)) {
		std::istringstream iss(line);
		Coord coord;
		std::string description;

		// Parse latitude/longitude from text.
		if (iss >> coord.lon >> coord.lat) {
			// Store the remainder of the line into the description field.
			getline(iss, description);
			coord.description = xtring(trim(description).c_str());
		}
		coordinates.push_back(coord);
	}

	return coordinates;
}

void NCInput::update_ndep(double lon, double lat, int year) {
	// Read monthly N deposition values from input file.
	double ndep_dry[12];
	double ndep_wet[12];
	ndep_input.get_one_calendar_year(year, ndep_dry, ndep_wet);

	// Converting N deposition data to daily values requires daily
	// precipitation data.

	// Read precipitation data for entire year.
	std::vector<double> subdaily_prec = precip->read(lon, lat, year, 0
	, Date::MAX_YEAR_LENGTH);

	// Aggregate to daily values.
	double daily_prec[Date::MAX_YEAR_LENGTH];

	// Iterate over each day of year and calculate a mean of the subdaily
	// values.
	for (int i = 0; i < Date::MAX_YEAR_LENGTH; i++) {
		double *start = &(subdaily_prec[i]);
		daily_prec[i] = mean(start, nsteps_day);
	}

	// Now aggregate to generate daily N deposition values.
	distribute_ndep(ndep_dry, ndep_wet, daily_prec, dndep);
}

int NCInput::get_nyear_data() const {
	int nyear = final_timestep.year - first_timestep.year;

	// If final timestep is not the last timestep in the year, then we don't
	// have a full year of data for the final year. Therefore decrement nyear.
	if (!is_end_of_year(final_timestep, nsteps_day, has_leap_years)) {
		nyear--;
	}

	// If first timestep is not the first timestep in the year, then we don't
	// have a full year of data for the first year. Therefore decrement nyear.
	if (!is_start_of_year(first_timestep, nsteps_day)) {
		nyear--;
	}

	// Clamp to lower bound of 0.
	return max(nyear, 0);
}

#endif // HAVE_NETCDF
