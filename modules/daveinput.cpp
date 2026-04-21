#ifdef HAVE_NETCDF

#include <array>
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
#include "daveinput.h"
#include "ncutil.h"
#include "parameters.h"
#include "guessmath.h"
#include "util.h" // Required by MSVC

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

REGISTER_INPUT_MODULE("nc", DaveInput)

DaveInput::DaveInput() {
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
	declare_parameter("file_thome", &file_thome, PATH_MAX
	, "Path to the mean maximum daily temperature (T_home) input file, used for DAVE temperature acclimation calculations. This is optional if iftempacclimation is disabled.");

	declare_parameter("cache_strategy", &cache_strategy, 0, 2
	, "Cache strategy to be used when reading NetCDF data (0 = cache entire file (fastest), 1 = cache one year at a time (slower, but uses less memory), 2 = no caching (slowest, not recommended)");

	declare_parameter("ndep_scenario", &ndep_scenario, 256
	, "Scenario to be used for nitrogen deposition data. E.g. HISTORIC, RCP26, RCP45, RCP60, RCP85, RCP85");

	declare_parameter("spinup_met_cycle_length", &spinup_met_cycle_length, 0, 1000000, "User input from the instruction file. The number of years of data to recycle during the freenyears. If 0, all data will be recycled during freenyears.");

	declare_parameter("nc_check_inputs", &nc_check_inputs, "1 to enable sensibility checking of met inputs, 0 otherwise");

	declare_parameter("use_mean_annual_co2", &use_mean_annual_co2, "If true, each day's CO2 value will be the annual mean co2. Otherwise, a mean of that day's subdaily data will be used.");

	declare_parameter("progress_interval", &progress_interval, -1, std::numeric_limits<int>::max(), "Interval in seconds between progress reports; 0 for yearly or -1 for daily");
	declare_parameter("progress_newline", &progress_newline, "Force the use of newline characters to terminate progress messages. Useful for MPI runs.");

	declare_parameter("strict_math", &strict_math, "Iff enabled, invalid mathematical operations (such as division by 0) will result in a fatal error.");

	declare_parameter("temperature_delta", &temperature_delta, 0, 1e6, "Optional constant to be added to daily temperature (°C).");
	declare_parameter("co2_delta", &co2_delta, 0, 1e6, "Optional constant to be added to daily co2 (ppm).");
	declare_parameter("precip_scalar", &precip_scalar, -1e6, 1e6, "Optional scalar to multiply daily precipitation by.");

	declare_parameter("fixed_pressure", &fixed_pressure, 0, 200000
	, "Optional fixed/constant air pressure to use in Pa");
	declare_parameter("esat_method", &esat_method_name, 256
	, "Method to use for estimating saturation vapor pressure");

	fixed_pressure = 0;

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
	bounds_vpd.min = 0.01; // CS updated from 0

	bounds_rh.max = 1;
	bounds_rh.min = 0;

	temperature_delta = 0;
	co2_delta = 0;
	precip_scalar = 1.0;

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

	progress_newline = false;
}

DaveInput::~DaveInput() {
	
}

void DaveInput::init(void) {
	try {
		double init_start = get_wall_time();

		// Parse the gridlist.
		if (file_gridlist == "") {
			fail("No gridlist file specified. Please set file_gridlist in the instruction file.");
		}

		gridlist = parse_gridlist(file_gridlist);

		// Open soil data file.
		soil_input.init(file_soildata);

		// Open misc input file.
		misc_input.init(file_gridlist);

		// Open management input files.
		management_input.init(file_gridlist);

		// Open landcover input file.
		landcover_input.init(file_gridlist);

		if (file_thome != "") {
			if (iftempacclimation) {
				thome_input.init(file_thome);
			} else {
				printf("WARNING: Temperature acclimation file specified but temperatue acclimation is disabled, so the file will not be used.\n");
			}
		} else {
			if (iftempacclimation) {
				fail("Temperature acclimation is enabled, but t_home input file is not specified. Please set file_thome in the instruction file.");
			}
			// else not specified and not needed - no problem.
		}

		esat_method = parse_esat_method(esat_method_name);

		// Move to first gridcell in the gridlist.
		current_gridcell = gridlist.begin();

		NCInputOptions options;
		options.file_met_forcing = file_met_forcing;
		options.file_met_spinup = file_met_spinup;
		options.nsteps_day = nsteps_day;
		options.cache_strategy = (CacheStrategy)cache_strategy;
		options.nc_check_inputs = nc_check_inputs;
		options.blaze_enabled = firemodel == BLAZE;
		options.allow_fixed_pressure = fixed_pressure > 0;
		options.allow_netcdf_co2 = file_co2 == "";
		options.allow_netcdf_ndep = !use_lamarque;
		options.constant_ps = fixed_pressure;
		options.use_mean_annual_co2 = use_mean_annual_co2;
		options.esat_method = esat_method;
		ncinput.init(options);

		if (!options.allow_netcdf_co2) {
			co2_file = std::unique_ptr<GlobalCO2File>(new GlobalCO2File());
			co2_file->load_file((const char *)file_co2);
		}

		date.subdaily = nsteps_day;

		// Set the calendar year at start of spinup.
		int first_year = ncinput.get_simulation_start_year();
		date.set_first_calendar_year(first_year - nyear_spinup);

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

bool DaveInput::getgridcell(Gridcell& gridcell) {
	try {
		double getgridcell_start = get_wall_time();

		if (current_gridcell == gridlist.end()) {
			double total_walltime = (getgridcell_start - ctor_start) * 1000;
			double prop = ncinput_walltime_use / total_walltime;
			printf("Simulations completed. Total time: %.2fs.                     \nTotal time spent in dave  input: %.2fs\n", total_walltime / 1000, ncinput_walltime_use / 1000);
			return false;
		}

		double lon = current_gridcell->lon;
		double lat = current_gridcell->lat;
		gridcell.set_coordinates(lon, lat);

		ncinput.begin_gridcell(lon, lat);

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

		if (iftempacclimation) {
			// Temperature acclimation is enabled - read thome value for this
			// gridcell.
			gridcell.climate.thome = thome_input.read_gridcell(lat, lon);
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

bool DaveInput::getclimate(Gridcell& gridcell) {
try {
	double getclimate_start = get_wall_time();

	// Shorthand access to some variables.
	Climate& climate = gridcell.climate;
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();

	if (date.year < nyear_spinup) {
		if (date.day == 0) {
			spinup_year = get_spinup_year();
		}
		year = spinup_year;
	} else {
		// We are now out of the spinup. Reset spinup_year to avoid confusion.
		spinup_year = -1;
	}

	// When reading data from the input files, we specify which day of the year
	// to read. The 29th of February is the 60th day of the year in leap years,
	// which would correspond to date.day == 59. This model does not simulate
	// leap years, so date.day = 59 in this model ALWAYS means the first of
	// March. If we pass 59 into NCVariable::read(), we will get data for the
	// 29th of Feb in a leap year (if the data contains leap years). This is
	// actually convenient, as it makes the outputs more directly interpretable.
	// Whatever we do, we need to skip a day of data in that year, and the
	// current implementation means the data for 31st December will not be used.
	// This actually corresponds better with the day of year value that is
	// written to the output files, and makes them easier to compare with other
	// data. Therefore, I will leave this message here as a warning.
	int day = date.day;

	climate.met_year = year;

	// Read main climate drivers from NCInput.
	if (!ncinput.read_day(year, day, climate)) {
		return false;
	}

	// Conditionally use external N Deposition data.
	update_ndep(year);
	if (use_lamarque) {
		climate.gridcell.dNH4dep = nh4_dep[day];
		climate.gridcell.dNO3dep = no3_dep[day];
	}

	// Conditionally use external CO2 data.
	if (!ncinput.has_co2()) {
		climate.co2 = get_co2();
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

	if (precip_scalar != 1.0) {
		for (int i = 0; i < date.subdaily; i++) {
			climate.precs[i] *= precip_scalar;
		}
		climate.prec *= precip_scalar;
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
		printf("%3d%% complete, %s elapsed, %s remaining%s"
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
		fprintf(stderr, "Failed to read data for gridcell (%f, %f)",
				gridcell.get_lon(), gridcell.get_lat());
		auto prev_gridcell = current_gridcell - 1;
		if (prev_gridcell >= gridlist.begin() && prev_gridcell->description != "") {
			fprintf(stderr, " (%s)", (const char *)prev_gridcell->description);
		}
		fprintf(stderr, "\n");
	}
	fail(get_nested_error_msg(e));
	return false; // To eliminate a compiler warning.
}
}

void DaveInput::getlandcover(Gridcell& gridcell) {
	// Proxy this request through to the landcover input module.
	landcover_input.getlandcover(gridcell);
}

void DaveInput::getmanagement(Gridcell& gridcell) {
	// Proxy this request through to the management input module.
	management_input.getmanagement(gridcell, landcover_input);
}

void DaveInput::getmiscinput_static(Gridcell& gridcell) {
	misc_input.getmiscinput_static(gridcell);
}

void DaveInput::getmiscinput_yearly(Gridcell& gridcell) {
	misc_input.getmiscinput_yearly(gridcell);
}

std::vector<DaveInput::Coord> DaveInput::parse_gridlist(xtring filename) const {
	// Open the file for reading.
	std::ifstream ifs(filename, std::ifstream::in);

	// Validate that the file was opened successfully.
	if (!ifs.is_open()) {
		throw std::runtime_error("Failed to open gridlist file: " + std::string(filename));
	}

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
			coordinates.push_back(coord);
		}
	}

	return coordinates;
}

void DaveInput::update_ndep(int year) {
	// Non-lamarque = NetCDF input, which is handled by NCInput.
	if (!use_lamarque || date.day != 0) {
		return;
	}

	// Read monthly N deposition values from input file.
	double nhx_dry[12];
	double noy_dry[12];
	double nhx_wet[12];
	double noy_wet[12];
	ndep_input.get_one_calendar_year(year, nhx_dry, noy_dry, nhx_wet, noy_wet);

	// Converting N deposition data to daily values requires daily
	// precipitation data.

	// Read precipitation data for entire year.
	std::vector<double> subdaily_prec = ncinput.read_annual_precip(year);

	// Aggregate to daily values. Note that precip is currently only allowed to
	// be in mm, so we can sum the subdaily values to get a daily one.
	int n = Date::MAX_YEAR_LENGTH;
	std::vector<double> daily_prec = aggregate(subdaily_prec, n, sum);

	assert(daily_prec.size() == Date::MAX_YEAR_LENGTH);

	// Distribute N deposition
	distribute_ndep(nhx_dry, noy_dry, nhx_wet, noy_wet,
					&(daily_prec.front()), nh4_dep, no3_dep);
}

double DaveInput::get_progress(bool since_restart) {
	// Number of gridcells which have been run.
		int cells_done = distance(gridlist.begin(), current_gridcell) - 1;

		// Number of years to simulate per gridcell.
		int years_to_simulate = nyear_spinup + ncinput.get_simulation_length();

		int skip = since_restart && restart ? state_year : 0;
		years_to_simulate -= skip;

		// Total number of years to be simulated across all gridcells.
		int total_years = years_to_simulate * gridlist.size();

		// Number of years which have been simulated (date.year starts at 0).
		int years_simulated = years_to_simulate * cells_done + date.year - skip;

		double progress = 1.0 * years_simulated / total_years;
		return progress;
}

double DaveInput::get_co2() {
	// Note: using calendar year, rather than met year to avoid cycling co2
	// concentration during the spinup period.
	return (*co2_file)[date.get_calendar_year()];
}

void DaveInput::resize_climate_arrays(Gridcell &gridcell, size_t size) const {
	for (Stand& stand : gridcell) {
		for (size_t i = 0; i < stand.nobj; i++) {
			Patch &patch = stand[i];
			Vegetation &vegetation = patch.vegetation;
			for (size_t j = 0; j < vegetation.nobj; j++) {
				Individual &indiv = vegetation[j];
				indiv.ps_subdaily.resize(size);
				indiv.lwp.resize(size);
			}
		}
	}
}

int DaveInput::get_spinup_year() const {

	int calendar_year = date.get_calendar_year();

	timepoint first_timestep = ncinput.get_first_timestep();
	timepoint spinup_start = ncinput.get_spinup_start();
	timepoint spinup_end = ncinput.get_spinup_end();

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
	if (!is_start_of_year(spinup_start)) {
		first_year++;
	}

	// Get the number of years of data we have for the spinup.
	int nyear_data = ncinput.get_spinup_length();

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

esat_method_t DaveInput::parse_esat_method(xtring name) const {
	name = name.lower();
	if (name == "sonntag_1990") {
		return SONNTAG_1990;
	} else if (name == "alduchov_1996") {
		return ALDUCHOV_1996;
	} else if (name == "allen_1998") {
		return ALLEN_1998;
	} else {
		fail("Invalid esat_method: %s", (const char *)name);
	}

	// Silence compiler warning (execution can never reach here).
	return (esat_method_t)-1;
}

#endif // HAVE_NETCDF
