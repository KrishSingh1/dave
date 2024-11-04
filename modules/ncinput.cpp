#ifdef HAVE_NETCDF

#include <exception>
#include <limits>
#include <limits.h>
#include <netcdf.h>
#include <stdexcept>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <fenv.h>

#include "driver.h"
#include "guessnc.h"
#include "guessstring.h"
#include "ncinput.h"
#include "ncutil.h"
#include "parameters.h"
#include "guessmath.h"
#include "util.h" // Required by MSVC

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

REGISTER_INPUT_MODULE("nc", NCInput)

NCInput::NCInput() {
	ctor_start = get_wall_time();

	// Register instruction file parameters.

	// Setting max timesteps per day to the number of seconds in a day.
	// Technically we could support more; this limitation is quite arbitrary.
	declare_parameter("nsteps_day", &nsteps_day, 1, SECONDS_PER_DAY
	, "Number of subdaily timesteps");

	// Again, this maximum file name length is pretty arbitrary and not really a
	// hard limit.
	declare_parameter("file_met_forcing", &file_met_forcing, PATH_MAX
	, "Path to NetCDF file containing met forcing data");
	declare_parameter("file_met_spinup", &file_met_spinup, PATH_MAX
	, "(Optional) Path to NetCDF file containing met forcing data to be used for spinup");
	declare_parameter("file_gridlist", &file_gridlist, PATH_MAX
	, "Path to Gridlist file containing lines of the form \"lon lat\"");
	declare_parameter("file_soildata", &file_soildata, PATH_MAX
	, "Path to soil input file");
	declare_parameter("file_ndep", &file_ndep, PATH_MAX
	, "Path to nitrogen deposition fastarchive input file");
	declare_parameter("file_co2", &file_co2, PATH_MAX
	, "Path to global co2 concentration input file. Either this or variable_co2 must be provided.");

	declare_parameter("cache_strategy", &cache_strategy, 0, 2
	, "Cache strategy to be used when reading NetCDF data (0 = cache entire file (fastest), 1 = cache one year at a time (slower, but uses less memory), 2 = no caching (slowest, not recommended)");

	declare_parameter("ndep_scenario", &ndep_scenario, 256
	, "Scenario to be used for nitrogen deposition data. E.g. HISTORIC, RCP26, RCP45, RCP60, RCP85, RCP85");

	declare_parameter("spinup_met_cycle_length", &spinup_met_cycle_length, 0, 1000000, "User input from the instruction file. The number of years of data to recycle during the freenyears. If 0, all data will be recycled during freenyears.");

	declare_parameter("met_sampling_seed", &met_sampling_seed, 0, std::numeric_limits<int>::max(), "Seed used for random sampling of met data during spinup period");

	declare_parameter("nc_check_inputs", &nc_check_inputs, "1 to enable sensibility checking of met inputs, 0 otherwise");

	declare_parameter("use_mean_annual_co2", &use_mean_annual_co2, "If true, each day's CO2 value will be the annual mean co2. Otherwise, a mean of that day's subdaily data will be used.");

	declare_parameter("progress_interval", &progress_interval, -1, std::numeric_limits<int>::max(), "Interval in seconds between progress reports; 0 for yearly or -1 for daily");
	declare_parameter("progress_newline", &progress_newline, "Force the use of newline characters to terminate progress messages. Useful for MPI runs.");

	declare_parameter("strict_math", &strict_math, "Iff enabled, invalid mathematical operations (such as division by 0) will result in a fatal error.");

	declare_parameter("temperature_delta", &temperature_delta, 0, 1e6, "Optional constant to be added to daily temperature (°C).");
	declare_parameter("co2_delta", &co2_delta, 0, 1e6, "Optional constant to be added to daily co2 (ppm).");

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
	declare_parameter("variable_wind", &variable_wind, NC_MAX_NAME
	, "Optional name of the wind speed variable in the NetCDF file");
	declare_parameter("variable_co2", &variable_co2, NC_MAX_NAME
	, "Name of the co2 variable in the NetCDF file. Optional if file_co2 is provided instead.");
	declare_parameter("variable_vpd", &variable_vpd, NC_MAX_NAME
	, "Name of the VPD variable in the NetCDF file");
	declare_parameter("variable_ps", &variable_ps, NC_MAX_NAME
	, "Name of the atmospheric pressure variable in the NetCDF file");
	declare_parameter("variable_tmax", &variable_tmax, NC_MAX_NAME
	, "Name of the max temperature variable in the NetCDF file");
	declare_parameter("variable_tmin", &variable_tmin, NC_MAX_NAME
	, "Name of the min temperature variable in the NetCDF file");
	declare_parameter("variable_ndep_nhx", &variable_ndep_nhx, NC_MAX_NAME
	, "Optional name of the nhx nitrogen deposition variable in the NetCDF file. If not set, file_ndep must be set instead.");
	declare_parameter("variable_ndep_noy", &variable_ndep_noy, NC_MAX_NAME
	, "Optional name of the noy nitrogen deposition variable in the NetCDF file. If not set, file_ndep must be set instead.");

	bounds_prec.max = 1000;
	bounds_prec.min = 0;

	bounds_wind.max = 100;
	bounds_wind.min = 0;

	bounds_insol.max = 1e6;
	bounds_insol.min = 0;

	bounds_co2.max = 900;
	bounds_co2.min = 0;

	bounds_temp.max = 100;
	bounds_temp.min = -100;

	bounds_vpd.max = 100;
	bounds_vpd.min = 0;

	temperature_delta = 0;
	co2_delta = 0;

	// Are these too generous?
	bounds_ps.max = 150000;
	bounds_ps.min = 50000;

	bounds_ndep.max = 1000;
	bounds_ndep.min = 0;

	use_lamarque = true;
	use_mean_annual_co2 = false;

	double ctor_end = get_wall_time();
	double ctor_duration = (ctor_end - ctor_start) * 1000;
	ncinput_walltime_use = ctor_duration;

	nc_file = -1;
	nc_temp = -1;
	nc_insol = -1;
	nc_precip = -1;
	nc_vpd = -1;
	nc_huss = -1;
	nc_ps = -1;
	nc_wind = -1;
	nc_co2 = -1;
	nc_tmax = -1;
	nc_tmin = -1;
	nc_nhx = -1;
	nc_noy = -1;

	progress_newline = false;
}

NCInput::~NCInput() {
	close_all();

	// The NCVariable instances are dynamically allocated, but they're wrapped
	// in smart pointers, so we don't need to manually free them.
}

void NCInput::init(void) {
	try {
		double init_start = get_wall_time();

		// Parse the gridlist.
		gridlist = parse_gridlist(file_gridlist);

		// Open soil data file.
		soil_input.init(file_soildata);

		// Open misc input file.
		misc_input.init(file_gridlist);

		// Open management input files.
		management_input.init(file_gridlist);

		// Open landcover input file.
		landcover_input.init(file_gridlist);

		bool separate_spinup = file_met_spinup != "" && file_met_spinup != file_met_forcing;
		open_netcdf(separate_spinup ? file_met_spinup : file_met_forcing);

		// Ensure that input timestep is an integer multiple of simulation
		// timestep. Note that if using a separate .nc file for spinup and
		// simulation, the timestep may be different between the two files. In
		// this case, we use the input timestep for the spinup, and then use
		// the nsteps_day timestep for the simulation.
		int input_timestep = temp->get_day_length();
		date.subdaily = separate_spinup ? input_timestep : nsteps_day;

		spinup_start = first_timestep = temp->get_first_date();
		spinup_end = final_timestep = temp->get_last_date();

		// Ensure that there is at least one year of data for the spinup.
		verify_data_length(first_timestep, final_timestep, input_timestep, temp->get_calendar());

		// Record length of input data.
		if (separate_spinup) {
			// The temperature variable currently loaded comes from the NetCDF
			// file used for the spinup, and its start/end dates don't reflect
			// the actual simulation start/end dates, which depend upon the
			// NetCDF file used for the main simulation.
			xtring fn;
			fn.printf(file_met_forcing, (const char *)variable_temp);
			int nc = GuessNC::open_ncdf(fn);

			// Read metadata from the temperature variable.
			// todo: verify that this is consistent across variables.
			bounds limits;
			limits.max = 100;
			limits.min = -100;
			std::vector<xtring> units(1);
			units[0] = "degC";
			NCVariable tmp(nc, variable_temp, variable_lon, variable_lat
				, variable_time, (CacheStrategy)cache_strategy
				, nc_check_inputs, limits, units, false);

			// Update simulation timesteps.
			first_timestep = tmp.get_first_date();
			final_timestep = tmp.get_last_date();

			// Set input_timestep to the value from file_met_forcing.
			input_timestep = tmp.get_day_length();

			// These three parameters all now come from file_met_forcing; last
			// time we called this function, they came from file_met_spinup.
			// Therefore we should check the main file. Better to fail early.
			verify_data_length(first_timestep, final_timestep, input_timestep, tmp.get_calendar());

			// Close file_met_forcing for now.
			GuessNC::close_ncdf(nc);
		}

		// Move to first gridcell in the gridlist.
		current_gridcell = gridlist.begin();

		// input_timestep is the timestep from file_met_forcing. Ensure that it
		// may be aggregated to nsteps_day.
		verify_timestep(input_timestep);

		// Now we must determine the first calendar year of the simulation (not
		// the spinup).
		int first_year = first_timestep.year;

		// Set the calendar year at start of spinup.
		date.set_first_calendar_year(first_year - nyear_spinup);

		// Initialise PRNG.
		srand(met_sampling_seed);

	if (strict_math) {
#ifndef _WIN32
		// feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
	}

	// Set timers.
	progress_timer.init();
	progress_progress_timer.init();

	progress_timer.settimer();
	if (progress_interval > 0) {
		progress_progress_timer.settimer(progress_interval);
	}

	double init_end = get_wall_time();
	double init_duration = (init_end - init_start) * 1000;
	ncinput_walltime_use += init_duration;

	} catch (const std::exception& error) {
		fail(get_nested_error_msg(error));
	}
}

bool NCInput::getgridcell(Gridcell& gridcell) {
	try {
		double getgridcell_start = get_wall_time();

		if (file_met_spinup != "" && current_gridcell != gridlist.begin()) {
			open_netcdf(file_met_spinup);

			// Ensure that input timestep is an integer multiple of simulation
			// timestep. Note that if using a separate .nc file for spinup and
			// simulation, the timestep may be different between the two files. In
			// this case, we use the input timestep for the spinup, and then use
			// the nsteps_day timestep for the simulation.
			int input_timestep = temp->get_day_length();
			date.subdaily = input_timestep;
		}

		if (current_gridcell == gridlist.end()) {
			double total_walltime = (getgridcell_start - ctor_start) * 1000;
			double prop = ncinput_walltime_use / total_walltime;
			dprintf("Simulations completed. Total time: %.2fs.                     \nTotal time spent in dave  input: %.2fs\n", total_walltime / 1000, ncinput_walltime_use / 1000);
			return false;
		}

		double lon = current_gridcell->lon;
		double lat = current_gridcell->lat;
		gridcell.set_coordinates(lon, lat);

		set_insoltype(gridcell);

		// Get soil data for this grid cell.
		soil_input.get_soil(lon, lat, gridcell);

		// Read N deposition data for this grid cell.
		if (use_lamarque) {
			ndep_input.getndep(file_ndep, lon, lat
			, Lamarque::parse_timeseries((const char *)ndep_scenario));
		}

		if(readdisturbance || readdisturbance_st || readelevation_st) {
			// Not all gridcells have to be included in input file.
			misc_input.loaddisturbance(lon, lat);
			misc_input.loadelevation(lon, lat);
		}

		// Initialise subdaily photosynthesis structures with the correct number
		// of elements corresponding to number of subdaily timesteps.
		resize_climate_arrays(gridcell, date.subdaily);

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
		if (current_gridcell > gridlist.begin() && current_gridcell < gridlist.end()) {
			fprintf(stderr, "Failed to read data for gridcell (%f, %f)", current_gridcell->lon, current_gridcell->lat);
			if (current_gridcell->description != "") {
				fprintf(stderr, " (%s)", (const char *)current_gridcell->description);
			}
			fprintf(stderr, "\n");
		}
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
	if (year > final_timestep.year) {
		return false;
	}
	if (date.day == 0 && year == final_timestep.year) {
		if (final_timestep.day < Date::MAX_YEAR_LENGTH - 1) {
			// We are in the final year, which is missing some days of data.
			// Because the N deposition disaggregation requires an entire years'
			// worth of precipitation data, we can't simulate part of a year.
			return false;
		}

		if (final_timestep.timestep < temp->get_day_length() - 1) {
			// We are in the final year, and the final day is missing some
			// timesteps.
			return false;
		}
	}

	if (date.year < nyear_spinup) {
		// Shuffle years of met data during spinup. Whenever we move to a new
		// year, select another year at random.
		if (date.day == 0) {
			spinup_year = get_spinup_year();
		}
		year = spinup_year;
	} else {
		// We are now out of the spinup. Reset spinup_year to avoid confusion.
		spinup_year = -1;

		if (date.day == 0 && date.year == nyear_spinup && file_met_spinup != "") {
			// Open the NetCDF file used for the main simulation.
			open_netcdf(file_met_forcing);
			set_insoltype(gridcell);
			// Input timestep for this .nc file was validated in init().
			date.subdaily = nsteps_day;
			resize_climate_arrays(gridcell, nsteps_day);
		}
	}

	climate.met_year = year;

	int co2_year = date.get_calendar_year() < spinup_start.year ? spinup_start.year : date.get_calendar_year();

	// Generate daily N deposition values on first day of year.
	if (date.day == 0) {
		update_ndep(lon, lat, year);

		if (use_mean_annual_co2 && co2) {
			int year_len = get_year_len(year, co2->get_calendar());
			std::vector<double> co2_annual = co2->read(lon, lat, co2_year, 0, year_len);
			mean_annual_co2 = mean(co2_annual);
		}
	}

	// Read data from variables (not necessarily from the file - these variable
	// objects will handle caching internally).
	// Note that these arrays will still be on the timestep of the input data,
	// which is not necessarily the same as that of the nsteps_day parameter.
	std::vector<double> temps = temp->read(lon, lat, year, date.day, 1);
	std::vector<double> precs = precip->read(lon, lat, year, date.day, 1);
	std::vector<double> insols = insol->read(lon, lat, year, date.day, 1);
	std::vector<double> vpds = vpd->read(lon, lat, year, date.day, 1);
	std::vector<double> pressure = ps->read(lon, lat, year, date.day, 1);

	// Aggregate subdaily values to daily where necessary, and store them in the
	// climate object. Note that sum inputs are rates and sum are quantities;
	// We can use the mean of the rates and the sum of the quantities. Check the
	// units for more info.
	climate.temp = mean(temps);
	climate.prec = sum(precs);
	climate.insol = mean(insols);

	// Aggregate from one subdaily timescale to another. Note that this will
	// only work if the input timestep is larger than the output timestep, a
	// condition which we have already verified is true in init().
	//
	// nsteps_day is the desired simulation timestep.
	// temp->get_day_length is the input data timestep.

	double maxt = max(temps);
	double mint = min(temps);

	if (date.subdaily != temp->get_day_length()) {
		temps = aggregate(temps, date.subdaily, mean);
		insols = aggregate(insols, date.subdaily, mean);
		vpds = aggregate(vpds, date.subdaily, mean);
		pressure = aggregate(pressure, date.subdaily, mean);
		precs = aggregate(precs, date.subdaily, sum);
		// If required, we can aggregate using other functions as well. E.g.
		// precs = aggregate(precs, date.subdaily, sum);
	}

	if (wind) {
		std::vector<double> winds = wind->read(lon, lat, year, date.day, 1);
		climate.u10 = mean(winds);
	}

	climate.temps = temps;
	climate.insols = insols;
	climate.vpds = vpds;
	climate.precs = precs;
	gridcell.dNH4dep = nh4_dep[date.day];
	gridcell.dNO3dep = no3_dep[date.day];
	climate.vpds = vpds;
	climate.pressure = pressure;

	climate.co2 = get_co2(lon, lat, co2_year, date.day);

	if (temp->get_day_length() > 1) {
		climate.dtr = maxt - mint;
		climate.tmax = maxt;
		climate.tmin = mint;
	} else if (tmax && tmin) {
		// It only really makes sense to provide separate variables for tmin and
		// tmax if input data is daily. Therefore these arrays will usually be
		// of length 1.
		std::vector<double> tmaxs = tmax->read(lon, lat, year, date.day, 1);
		std::vector<double> tmins = tmin->read(lon, lat, year, date.day, 1);

		climate.tmax = mean(tmaxs);
		climate.tmin = mean(tmins);
		climate.dtr = climate.tmax - climate.tmin;
	} else {
		throw std::runtime_error("Missing tmax/tmin data, and day length is 1");
	}

	if (temperature_delta > 0) {
		for (int i = 0; i < date.subdaily; i++) {
			climate.temps[i] += temperature_delta;
		}
		climate.temp += temperature_delta;
	}

	if (co2_delta > 0) {
		climate.co2 += co2_delta;
	}

	// rads and pars are calculated automatically in daylengthinsoleet().
	// gtemps are calculated in respiration_temperature_response().

	double getclimate_duration = (get_wall_time() - getclimate_start) * 1000;
	ncinput_walltime_use += getclimate_duration;

	// Progress reporting.
	if (progress_interval == -1 ||
		(progress_interval == 0 && date.day == 0) ||
		(date.day == 0 && progress_progress_timer.getprogress() >= 1.0) ) {
		double overall_progress = get_progress(false);
		double progress = overall_progress;
		if (restart) {
			progress = get_progress(true);
		}
		progress_timer.setprogress(progress);
		const char *eol = isatty(fileno(stdout)) ? "\r" : "\n";
		if (progress_newline) {
			eol = "\n";
		}
		dprintf("%3d%% complete, %s elapsed, %s remaining%s"
			, (int)(overall_progress * 100), progress_timer.elapsed.str
			, progress_timer.remaining.str, eol);
		if (progress_interval > 0) {
			progress_progress_timer.settimer(progress_interval);
		}
	}

	// We have enough data to simulate another day.
	return true;
} catch (const std::exception& e) {
	if (current_gridcell > gridlist.begin() && current_gridcell < gridlist.end()) {
		fprintf(stderr, "Failed to read data for gridcell (%f, %f)", current_gridcell->lon, current_gridcell->lat);
		if (current_gridcell->description != "") {
			fprintf(stderr, " (%s)", (const char *)current_gridcell->description);
		}
		fprintf(stderr, "\n");
	}
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
	management_input.getmanagement(gridcell, landcover_input);
}

void NCInput::getmiscinput_static(Gridcell& gridcell) {
	misc_input.getmiscinput_static(gridcell);
}

void NCInput::getmiscinput_yearly(Gridcell& gridcell) {
	misc_input.getmiscinput_yearly(gridcell);
}

std::unique_ptr<NCVariable> NCInput::create_variable(xtring netcdf_file
	, xtring variable_name, int &ncid, bounds limits, const std::vector<xtring>& units
	, bool require_valid_units) const {

	// nc is the ID of the file from which we will read.
	int nc;

	if (ncid != -1) {
		// We have already opened a variable-specific file. This should probably
		// never happen.
		nc = ncid;
	} else {
		// Attempt to get a variable-specific file by treating the file name as
		// a printf-style format string.
		xtring filename;
		filename.printf((const char *)netcdf_file, (const char *)variable_name);

		// nc_open() can be called multiple times on the same file, and will
		// return the same file ID each time.
		nc = GuessNC::open_ncdf((const char *)filename);

		// Store the opened file ID in ncid, but only if it's different to the
		// 'primary' open file. If it's the same file, just store -1 in ncid.
		ncid = nc == nc_file ? -1 : nc;
	}

	return std::unique_ptr<NCVariable>(new NCVariable(
		nc,
		variable_name,
		variable_lon,
		variable_lat,
		variable_time,
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
	if (use_lamarque) {
		// Read monthly N deposition values from input file.
		double nhx_dry[12];
		double noy_dry[12];
		double nhx_wet[12];
		double noy_wet[12];
		ndep_input.get_one_calendar_year(year, nhx_dry, noy_dry, nhx_wet, noy_wet);

		// Converting N deposition data to daily values requires daily
		// precipitation data.

		// Read precipitation data for entire year.
		int n = Date::MAX_YEAR_LENGTH;
		std::vector<double> subdaily_prec = precip->read(lon, lat, year, 0, n);

		// Aggregate to daily values. Note that precip is currently only allowed to
		// be in mm, so we can sum the subdaily values to get a daily one.
		std::vector<double> daily_prec = aggregate(subdaily_prec, n, sum);

		// Distribute N deposition
		distribute_ndep(nhx_dry, noy_dry, nhx_wet, noy_wet,
						&(daily_prec.front()), nh4_dep, no3_dep);
	} else {
		// todo: leap years.
		int year_len = Date::MAX_YEAR_LENGTH;

		// Read n deposition data from NetCDF input.
		std::vector<double> nhx_subdaily = ndep_nhx->read(lon, lat, year, 0, year_len);
		std::vector<double> noy_subdaily = ndep_noy->read(lon, lat, year, 0, year_len);

		// Aggregate to daily timestep. Input units are mgN/m2/day, so we take
		// the mean rate over the course of a day, and that will give us the
		// amount of N/m2 deposited in a day.
		std::vector<double> nhx = aggregate(nhx_subdaily, year_len, mean);
		std::vector<double> noy = aggregate(noy_subdaily, year_len, mean);

		// Copy data into the nh4_dep and no3_dep arrays.
		// nhx = nh4, noy = no3a
		memcpy(nh4_dep, &nhx.front(), year_len * sizeof(double));
		memcpy(no3_dep, &noy.front(), year_len * sizeof(double));
	}
}

int NCInput::get_nyear_data(timepoint start, timepoint end, int day_length, GuessNC::CF::CalendarType calendar) const {
	int nyear = end.year - start.year;

	// If final timestep is the last timestep in a year, then we can add that
	// simulate that entire year (so increment nyear).
	if (is_end_of_year(end, day_length, calendar)) {
		nyear++;
	}

	// If first timestep is not the first timestep in the year, then we don't
	// have a full year of data for the first year. Therefore decrement nyear.
	if (!is_start_of_year(start, day_length)) {
		nyear--;
	}

	// Clamp to lower bound of 0.
	return max(nyear, 0);
}

double NCInput::get_progress(bool since_restart) {
	// Number of gridcells which have been run.
		int cells_done = distance(gridlist.begin(), current_gridcell) - 1;

		// Number of years to simulate per gridcell.
		int years_to_simulate = nyear_spinup + get_nyear_data(first_timestep, final_timestep, nsteps_day, temp->get_calendar());

		int skip = since_restart && restart ? state_year : 0;
		years_to_simulate -= skip;

		// Total number of years to be simulated across all gridcells.
		int total_years = years_to_simulate * gridlist.size();

		// Number of years which have been simulated (date.year starts at 0).
		int years_simulated = years_to_simulate * cells_done + date.year - skip;

		double progress = 1.0 * years_simulated / total_years;
		return progress;
}

double NCInput::get_co2(double lon, double lat, int year, int day) {
	if (co2) {
		return use_mean_annual_co2
			? mean_annual_co2
			: mean(co2->read(lon, lat, year, date.day, 1));
	} else {
		// Note: using calendar year, rather than met year to avoid cycling co2
		// concentration during the spinup period.
		return (*co2_file)[date.get_calendar_year()];
	}
}

void NCInput::close_if_open(int &nc) {
	if (nc != -1) {
		GuessNC::close_ncdf(nc);
		nc = -1;
	}
}

void NCInput::close_all() {
	close_if_open(nc_file);
	close_if_open(nc_temp);
	close_if_open(nc_insol);
	close_if_open(nc_precip);
	close_if_open(nc_vpd);
	close_if_open(nc_huss);
	close_if_open(nc_ps);
	close_if_open(nc_wind);
	close_if_open(nc_co2);
	close_if_open(nc_tmax);
	close_if_open(nc_tmin);
	close_if_open(nc_nhx);
	close_if_open(nc_noy);
}

void NCInput::open_netcdf(xtring netcdf_file) {

	// Close the previously open file(s).
	close_all();

	// Open the met forcing file (in readonly mode).
	xtring tmp;
	tmp.printf(netcdf_file, (const char *)variable_temp);
	if (tmp == netcdf_file) {
		nc_file = GuessNC::open_ncdf((const char *)netcdf_file);
	}

	std::vector<xtring> temp_units(1);
	temp_units[0] = "degC";
	std::vector<xtring> insol_units(3);
	insol_units[0] = "W/m2";
	insol_units[1] = "W m-2";
	insol_units[2] = "w/m^2";
	// Should provide units for other insol types as well...
	std::vector<xtring> precip_units(3);
	precip_units[0] = "mm";
	precip_units[1] = "kg m-2";
	precip_units[2] = "kg/m2";

	std::vector<xtring> wind_units(3);
	wind_units[0] = "m/s";
	wind_units[1] = "m s-1";
	wind_units[2] = "m s^-1";

	std::vector<xtring> co2_units(1);
	co2_units[0] = "ppm";

	std::vector<xtring> vpd_units(1);
	vpd_units[0] = "kPa";

	std::vector<xtring> ndep_units(1);
	ndep_units[0] = "kgN/m2/day";

	if (nsteps_day >= SIZE_MAX) {
		throw_xtring("nsteps_day is too large: %d", nsteps_day);
	}

	std::vector<xtring> ps_units(1);
	ps_units[0] = "Pa";

	// Create NCVariable objects for each variable of interest.
	// This will also do some basic sanity checking.
	temp = create_variable(netcdf_file, variable_temp, nc_temp, bounds_temp, temp_units, true);
	insol = create_variable(netcdf_file, variable_insol, nc_insol, bounds_insol, insol_units, true);
	precip = create_variable(netcdf_file, variable_precip, nc_precip, bounds_prec, precip_units, true);
	vpd = create_variable(netcdf_file, variable_vpd, nc_vpd, bounds_vpd, vpd_units, true);
	ps = create_variable(netcdf_file, variable_ps, nc_ps, bounds_ps, ps_units, true);


	if (variable_wind != "" && variable_exists(netcdf_file, variable_wind)) {
		wind = create_variable(netcdf_file, variable_wind, nc_wind, bounds_wind, wind_units, true);
	} else {
		wind = NULL;
		if (firemodel == BLAZE) {
			fail("When using BLAZE, wind speed is not optional. Please specify variable_wind in the instruction file.");
		}
	}

	int day_length = temp->get_day_length();

	if (variable_co2 == "") {
		if (file_co2 == "") {
			fail("Must set either file_co2 or variable_co2 in ins file");
		}
		co2_file = std::unique_ptr<GlobalCO2File>(new GlobalCO2File());
		co2_file->load_file((const char *)file_co2);
		co2 = NULL;
	} else {
		co2 = create_variable(netcdf_file, variable_co2, nc_co2, bounds_co2, co2_units, true);
	}
	if (variable_tmax != "" && variable_exists(netcdf_file, variable_tmax)) {
		tmax = create_variable(netcdf_file, variable_tmax, nc_tmax, bounds_temp, temp_units, true);
	} else if (day_length < 2) {
		fail("Variable '%s' does not exist in %s and timestep is daily", (const char *)variable_tmax, (const char *)netcdf_file);
	} else {
		tmax.release();
		tmax = NULL;
	}
	if (variable_tmin != "" && variable_exists(netcdf_file, variable_tmin)) {
		tmin = create_variable(netcdf_file, variable_tmin, nc_tmin, bounds_temp, temp_units, true);
	} else if (day_length < 2) {
		fail("Variable '%s' does not exist in %s and timestep is daily", (const char *)variable_tmax, (const char *)netcdf_file);
	} else {
		tmin.release();
		tmin = NULL;
	}
	if (variable_ndep_nhx == "") {
		use_lamarque = true;
	} else {
		use_lamarque = false;
		ndep_nhx = create_variable(netcdf_file, variable_ndep_nhx, nc_nhx, bounds_ndep, ndep_units, true);
	}
	if (variable_ndep_noy == "") {
		use_lamarque = true;
	} else {
		use_lamarque = false;
		ndep_noy = create_variable(netcdf_file, variable_ndep_noy, nc_noy, bounds_ndep, ndep_units, true);
	}
}

void NCInput::set_insoltype(Gridcell &gridcell) const {
		// We can get instype from insol's standard_name attribute.
		xtring insol_name = insol->get_standard_name();
		gridcell.climate.instype = parse_insoltype((const char*)insol_name);
}

void NCInput::verify_timestep(int timestep) const {
	if (timestep < nsteps_day) {
		throw_xtring("Input timestep (%d/day) must be >= simulation timestep (%d/day)"
		, timestep, nsteps_day);
	}
	if (timestep % nsteps_day != 0) {
		throw_xtring("Input timestep (%d/day) must be an integer multiple of simulation timestep (%d/day)"
		, timestep, nsteps_day);
	}
}

void NCInput::resize_climate_arrays(Gridcell &gridcell, size_t size) const {
	for (Stand& stand : gridcell) {
		for (size_t i = 0; i < stand.nobj; i++) {
			Patch &patch = stand[i];
			Vegetation &vegetation = patch.vegetation;
			for (size_t j = 0; j < vegetation.nobj; j++) {
				Individual &indiv = vegetation[j];
				indiv.ps_subdaily.resize(size);
			}
		}
	}
}

void NCInput::verify_data_length(timepoint start, timepoint end, int timestep, GuessNC::CF::CalendarType calendar) const {
	// Ensure the simulation file contains at least 1 year of data.
	int nyear_data = get_nyear_data(start, end, timestep, calendar);
	if (nyear_data < 1) {
		throw std::runtime_error("Insufficient met forcing data (<1 year)");
	}
}

bool NCInput::variable_exists(xtring netcdf_file, xtring name) const {
	// Note: the NetCDF API does accept NULL for the pointer argument, but in
	// that case it will return 0 even if the variable doesn't exist.
	xtring var_file_name;
	var_file_name.printf(netcdf_file, (const char *)name);
	int ncid = nc_file;
	if (var_file_name != netcdf_file) {
		// The variables are stored in separate files.
		if (!exists(var_file_name)) {
			// ...however no file exists for this variable.
			return false;
		}
		// ... and a file exists for this variable.
		ncid = GuessNC::open_ncdf(var_file_name);
	}

	int id;
	int status = nc_inq_varid(ncid, (const char *)name, &id);

	if (var_file_name != netcdf_file) {
		GuessNC::close_ncdf(ncid);
	}

	if (status == NC_NOERR) {
		return true;
	}

	if (status == NC_ENOTVAR) {
		return false;
	}

	// Maybe we should return false here. I'm not sure what else could cause an
	// error here anyway.
	fail("NetCDF error while checking if variable %s exists: %s", (const char *)name, nc_strerror(status));

	// Execution will never reach here.
	return false;
}

int NCInput::get_spinup_year() const {

	int calendar_year = date.get_calendar_year();
	if (calendar_year >= first_timestep.year ||
		(calendar_year >= spinup_start.year && calendar_year <= spinup_end.year) ) {
		// We are past the start of the data. Therefore use data for the current
		// calendar year.
		return calendar_year;
	} else if (calendar_year > spinup_end.year && calendar_year < first_timestep.year) {
		return spinup_end.year;
	}

	// Skip the first year if the first timestep lies partway through the year.
	int first_year = spinup_start.year;
	if (!is_start_of_year(spinup_start, temp->get_day_length())) {
		first_year++;
	}

	// Get the number of years of data we have for the spinup.
	int nyear_data = get_nyear_data(spinup_start, spinup_end, temp->get_day_length(), temp->get_calendar());

	// User wants to limit the spinup to using only the first N years of data.
	if (spinup_met_cycle_length > 0) {
		nyear_data = min(nyear_data, spinup_met_cycle_length);
	}

	// The offset from the start of the met data to be used in the first year.
	int cycle_length = nyear_spinup - (first_timestep.year - spinup_start.year);
	int offset = nyear_data - modulo(cycle_length, nyear_data);

	// Add offset to current year, ensure it lies within range (0..nyear_data],
	// and add this to the first year of data. This gives us the (calendar) year
	// of data to be used this year.
	return (date.year + offset) % nyear_data + first_year;
}

#endif // HAVE_NETCDF
