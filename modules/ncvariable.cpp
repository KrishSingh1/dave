#ifdef HAVE_NETCDF

#include <algorithm>
#include <netcdf.h>
#include <vector>

#include "cftime.h"
#include "guessmath.h"
#include "guessnc.h"
#include "gutil.h"
#include "ncutil.h"
#include "ncvariable.h"

int NCVariable::get_dimension_id(xtring dim_name) const {
	int dim_id;
	int status = nc_inq_dimid(id_file, (const char *)dim_name, &dim_id);

	// Handle an error (if an error occurred).
	handle_nc_error(status, "Failed to get ID of dimension '%s'", (const char *)dim_name);

	return dim_id;
}

int NCVariable::get_variable_id(xtring var_name) const {
	int var_id;
	int status = nc_inq_varid(id_file, (const char *)var_name, &var_id);

	// Handle an error (if an error occurred).
	handle_nc_error(status, "Failed to get variable '%s'", (const char *)var_name);

	return var_id;
}

int NCVariable::get_num_dimensions(int var_id) const {
	int ndims;
	int status = nc_inq_varndims(id_file, var_id, &ndims);

	handle_nc_error(status, "Failed to read number of dimensions of variable %s"
	, (const char *)var_name);

	return ndims;
}

void NCVariable::check_num_dims(int var_id, int ndim_expected) const {
	// Get the number of dimensions for this variable.
	int ndims = get_num_dimensions(var_id);

	// Throw error if not correct.
	if (ndims != ndim_expected) {
		char name[NC_MAX_NAME];
		nc_inq_varname(id_file, var_id, name);
		const char *fmt = "Invalid number of dimensions in variable %s; expected %d but was %d";
		throw_xtring(fmt, (const char *)name, ndim_expected, ndims);
	}
}

xtring NCVariable::get_units() const {
	std::string result;
	if (GuessNC::get_attribute(id_file, id_var, ATTR_UNITS, result)) {
		return xtring(result.c_str());
	}
	else {
		// The units attribute isn't required by the CF standard
		throw_xtring("Variable %s doesn't have a units attribute", (const char *)var_name);
	}
}

void NCVariable::check_units(const std::vector<xtring>& expected, bool fatal) const {
	xtring actual = get_units().lower();
	for (const xtring& expect : expected) {
		xtring cpy = expect;
		cpy = cpy.lower();
		if (cpy == actual) {
			return;
		}
	}
	xtring msg;
	msg.printf("Variable %s has incorrect units. Expected %s but was %s"
	, (const char *)var_name, (const char *)expected[0]
	, (const char *)actual);

	// todo: automatic unit conversion (where feasible).
	dprintf("TBI: Automatic unit conversion\n");
	fatal = true; // fixme

	if (fatal) {
		throw std::runtime_error((const char *)msg);
	} else {
		dprintf("WARNING: %s. Values will be converted to %s, but consider fixing the input data for faster simulations."
		, (const char *)msg, (const char *)actual);
	}
}

NCVariable::NCVariable(int file_id, xtring name, xtring lon_name
, xtring lat_name, xtring time_name, bool has_leap
, CacheStrategy chunk_mode, bool sensibility_checks, bounds b
, std::vector<xtring> units, bool force_legal_units) {
	try {
		id_file = file_id;
		var_name = name;
		cache_strategy = chunk_mode;
		id_var = get_variable_id(name);
		has_leap_years = has_leap;
		check_bounds = sensibility_checks;
		limits = b;

		ndim = get_num_dimensions(id_var);

		// Get dimensions of this variable.
		std::vector<int> dimensions(ndim);
		int status = nc_inq_vardimid(id_file, id_var, &(dimensions.front()));
		GuessNC::handle_error(status, "Failed to get dimensions");

		// Get longitude variable ID.
		id_lon = get_variable_id(lon_name);
		check_num_dims(id_lon, 1);
		int id_lon_dim = get_dimension_id(lon_name);

		// Get latitude variable ID.
		id_lat = get_variable_id(lat_name);
		check_num_dims(id_lat, 1);
		int id_lat_dim = get_dimension_id(lat_name);

		// Get time variable ID.
		id_time = get_variable_id(time_name);
		int id_time_dim = get_dimension_id(time_name);
		check_num_dims(id_time, 1);

		// Check total number of timesteps. This will be useful later.
		ntime = get_dimension_length(id_time_dim);

		// Determine and validate the calendar type used by this NetCDF file.
		calendar_type = get_calendar_type();

		// Determine and validate time units.
		time_spec = get_time_spec();

		// Calculate and store the number of timesteps per day.
		day_length = calculate_day_length();
		if (day_length == -1) {
			throw std::runtime_error("File contains only a single day of data");
		}

		// Get the first time step.
		data_start = get_first_date();

		// Get pointer to end of the dimensions array.
		int *dims_end = &(dimensions.front()) + ndim * sizeof(int);

		// Get the ordering of the dimensions for this variable.
		for (int i = 0; i < ndim; i++) {
			if (dimensions[i] == id_lon_dim) {
				dim_lon = i;
			} else if (dimensions[i] == id_lat_dim) {
				dim_lat = i;
			} else if (dimensions[i] == id_time_dim) {
				dim_time = i;
			}
			// Technically, there could be other dimensions as well. A common
			// one would be height. We're going to ignore these for now.
		}

		// Unit checking.
		check_units(units, force_legal_units);

		// Initialise cache.
		cache_prev.year = -1;
		cache_prev.day = -1;
		cache_prev.timestep = -1;
	} catch (...) {
		// Rethrow the exception with some extra context added to the message.
		throw_xtring_nested("Failed to initialise variable %s", (const char *)name);
	}
}

std::vector<double> NCVariable::read(double lon, double lat, int year, int day
, int nday) {
	try {
		int result_len = day_length * nday;
		std::vector<double> result(result_len);
		int read = 0;

		// Reading in chunks. Remember, cache size may be less than the result
		// length.
		while (read < result_len) {
			// If the cache is not valid for this combination of lon/lat/date,
			// update the cache.
			if (!is_cache_valid(lon, lat, year, day)) {
				// Load data from file according to the chosen caching strategy.
				cache = load(lon, lat, year, day);

				// Update starting date of cached data.
				cache_start.year = cache_strategy == All ? data_start.year : year;
				cache_start.day = cache_strategy == Year ? 0 : day;
				cache_start.timestep = 0;
				longitude = lon;
				latitude = lat;
			}

			timepoint date;
			date.year = year;
			date.day = day;
			date.timestep = 0;

			int start_index = get_cache_offset(date);

			// This is an exclusive index.
			int end_index = min(start_index + result_len, cache.size());

			cache_prev.year = year;
			cache_prev.day = day;
			cache_prev.timestep = 0;
			cache_prev_index = start_index;

			if (start_index < 0 || start_index + day_length > cache.size()) {
				throw std::runtime_error("Insufficient data in file");
			}

			// Copy N elements (up to cache size) into the result vector.
			std::vector<double>::iterator start = cache.begin() + start_index;
			std::vector<double>::iterator end = cache.begin() + end_index;
			std::vector<double>::iterator result_start = result.begin() + read;
			std::copy(start, end, result_start);

			read += end_index - start_index;
		}

		if (check_bounds) {
			try {
				bounds_check(limits, result);
			} catch (...) {
				throw_xtring_nested("Bounds check failure for variable %s"
					, (const char *)var_name);
			}
		}

		return result;
	} catch (...) {
		throw_xtring_nested(
			"Unable to read %s data for gridpoint (%.2f,%.2f) year %d day %d",
			(const char *)var_name, lon, lat, year, day);
	}
}

std::vector<double> NCVariable::load(double lon, double lat, int year, int day)
const {
	if (cache_strategy == None) {
		// No caching. Read the data from the file.
		return load_day(lon, lat, year, day);
	} else if (cache_strategy == Year) {
		return load_year(lon, lat, year);
	} else if (cache_strategy == All) {
		return load_all(lon, lat);
	}
	throw std::runtime_error("Unknown caching strategy: " + std::to_string((int)cache_strategy));
}

std::vector<double> NCVariable::load_day(double lon, double lat, int year, int day) const {
	timepoint t;
	t.year = year;
	t.day = day;
	t.timestep = 0;
	int start = get_timesteps_between(data_start, t);
	return load_data(lon, lat, start, day_length);
}

std::vector<double> NCVariable::load_year(double lon, double lat, int year) const {
	timepoint t;
	t.year = year;
	t.day = t.timestep = 0;
	int start = get_timesteps_between(data_start, t);
	int len = get_year_len(year, has_leap_years) * day_length;
	return load_data(lon, lat, start, len);
}

std::vector<double> NCVariable::load_all(double lon, double lat) const {
	return load_data(lon, lat, 0, ntime);
}

std::vector<double> NCVariable::load_data(double lon, double lat, int time_start
, int ntimesteps) const {

	int index_lat = index_of(lat, id_lat);
	int index_lon = index_of(lon, id_lon);

	if (index_lat == -1 || index_lon == -1) {
		throw_xtring("Grid point (%.2f, %.2f) not found", lon, lat);
	}

	std::vector<size_t> start(ndim);
	std::vector<size_t> count(ndim);
	std::vector<ptrdiff_t> imap(ndim);
	for (int i = 0; i < ndim; i++) {
		start[i] = 0;
		count[i] = 1;
		imap[i] = 1;
	}
	start[dim_lon] = index_lon;
	start[dim_lat] = index_lat;
	start[dim_time] = time_start;

	count[dim_time] = ntimesteps;

	std::vector<double> data(ntimesteps);
	int status = nc_get_varm_double(id_file, id_var, &(start.front())
		, &(count.front()), NULL, &(imap.front()), &(data.front()));
	GuessNC::handle_error(status, "Failed to load data from file");
	return data;
}

int NCVariable::get_cache_offset(timepoint date) const {
	if (cache_strategy == None) {
		return date.timestep;
	}
	// Shortcut - use previous index + day_length.
	if (date.timestep == 0
	&&	cache_prev.timestep == 0
	// Same year, next day.
	&&	( (cache_prev.year == date.year
	&&	cache_prev.day + 1 == date.day)
	// Or cache_prev is last day of year and date is 1st day of next year.
	||	(cache_prev.year == date.year - 1
	&&	cache_prev.day == get_year_len(cache_prev.year, has_leap_years) - 1
	&&	date.day == 0
	&&  cache_strategy == All) )) {
		return cache_prev_index + day_length;
	}
	// Year offset is complicated by the possibility of leap years.
	return get_timesteps_between(cache_start, date);
}

int NCVariable::get_timesteps_between(timepoint date0, timepoint date1) const {

	if (date1.year < date0.year
	|| (date1.year == date0.year && date1.day < date0.day)
	|| (date1.year == date0.year && date1.day == date0.day && date1.timestep < date0.timestep)) {
		throw_xtring("Date %04d/%d/%d is earlier than %04d/%d/%d", date0.year
		, date0.day, date0.timestep, date1.year, date1.day, date1.timestep);
	}

	int offset = 0;
	for (int i = date0.year; i < date1.year; i++) {
		offset += get_year_len(i, has_leap_years) * day_length;
	}

	offset += (date1.day - date0.day) * day_length + date1.timestep - date0.timestep;

	return offset;
}

xtring NCVariable::get_standard_name() const {
	std::string attribute;
	if (!GuessNC::get_attribute(id_file, id_var, "standard_name", attribute)) {
		throw std::runtime_error("Unable to read standard_name attribute");
	}
	return xtring(attribute.c_str());
}

bool NCVariable::is_cache_valid(double lon, double lat, int year, int day) const {
	// If cache is empty, or a different location is cached, it's invalid.
	if (cache.size() == 0
	|| !floats_equal(lon, longitude)
	|| !floats_equal(lat, latitude)) {
		return false;
	}

	if (cache_strategy == All) {
		return true;
	}

	if (cache_strategy == Year) {
		return cache_start.year == year;
	}

	return cache_start.day == day && cache_start.year == year \
		&& cache_start.timestep == 0;
}

xtring NCVariable::get_dimension_name(int id_dim) const {
	char name[NC_MAX_NAME + 1];
	int status = nc_inq_dimname(id_file, id_dim, name);
	return xtring(name);
}

size_t NCVariable::get_dimension_length(int id_dim) const {
	size_t result;
	int status = nc_inq_dimlen(id_file, id_dim, &result);

	handle_nc_error(status
	, "Unable to read length of dimension %s"
	, (const char *)get_dimension_name(id_dim));

	return result;
}

int NCVariable::index_of(double x, int var_id) const {
	// We're assuming that the variable is 1-dimensional. First, let's verify
	// that assumption.
	check_num_dims(var_id, 1);

	// Now get the ID of the one and only dimension.
	int id_dim;
	int status = nc_inq_vardimid(id_file, var_id, &id_dim);
	handle_nc_error(status, "Failed to get dimensions of variable %s"
	, GuessNC::get_variable_name(id_file, var_id).c_str());

	// Get the length of the dimension.
	size_t length = get_dimension_length(id_dim);

	nc_type var_type;
	status = nc_inq_vartype(id_file, var_id, &var_type);
	handle_nc_error(status, "Failed to get type of variable %s"
	, GuessNC::get_variable_name(id_file, var_id).c_str());

	if (var_type == NC_DOUBLE) {
		// Read all data from this variable along its only dimension.
		std::vector<double> data(length);
		status = nc_get_var_double(id_file, var_id, &(data.front()));
		handle_nc_error(status, "Failed to read data for %s"
		, GuessNC::get_variable_name(id_file, var_id).c_str());

		// Finally we can search the data for the specified input value (x).
		return index_of_double(x, &(data.front()), data.size());
	}

	if (var_type == NC_FLOAT) {
		// Read all data from this variable along its only dimension.
		std::vector<float> data(length);
		status = nc_get_var_float(id_file, var_id, &(data.front()));
		handle_nc_error(status, "Failed to read data for %s"
		, GuessNC::get_variable_name(id_file, var_id).c_str());

		// Finally we can search the data for the specified input value (x).
		return index_of_float(x, &(data.front()), data.size());
	}

	// Not found.
	return -1;
}

timepoint NCVariable::to_timepoint(const GuessNC::CF::DateTime& date) const {
	timepoint result;
	result.year = date.get_year();

	const double days_in_month[12] =
		{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	// Calculate day of year. Note that get_day() returns a 1-31 number, but we
	// want a 0-based index, so we need to subtract 1.
	int day = date.get_day() - 1;
	for (int i = 0; i < date.get_month() - 1; i++) {
		int month_length = days_in_month[i];
		if (i == 1 && is_leap(result.year)) {
			// February.
			month_length = 29;
		}
		day += days_in_month[i];
	}
	result.day = day;

	// Get the time as fractions of an hour (e.g. 08:30:45 = 8.7625).
	double hour = date.get_hour() + date.get_minute() / MINUTES_PER_HOUR
				+ date.get_second() / SECONDS_PER_HOUR;
	double timestep_width = 1.0 * HOURS_PER_DAY / day_length;
	result.timestep = round_to_nearest_multiple(hour, timestep_width);
	if (result.timestep >= day_length) {
		result.timestep = 0;
		result.day++;
		if (result.day >= get_year_len(result.year, has_leap_years)) {
			result.day = 0;
			result.year++;
		}
	}

	return result;
}

GuessNC::CF::DateTime NCVariable::get_nth_date(size_t n) const {
	double value;
	int status = nc_get_var1_double(id_file, id_time, &n, &value);
	GuessNC::handle_error(status, "Failed to read from time variable");

	// Parse this time value according to the specified calendar type.
	return time_spec.get_date_time(value, calendar_type);
}

timepoint NCVariable::get_first_date() const {
	return to_timepoint(get_nth_date(0));
}

timepoint NCVariable::get_last_date() const {
	return to_timepoint(get_nth_date(ntime - 1));
}

int NCVariable::calculate_day_length() const {
	int first_day = get_nth_date(0).get_day();

	// Keep moving forward in time until we reach a different day.
	for (size_t i = 1; i < ntime; i++) {
		GuessNC::CF::DateTime timestep = get_nth_date(i);
		if (timestep.get_day() != first_day) {
			return i;
		}
	}
	return -1;
}

int NCVariable::get_day_length() const {
	return day_length;
}

GuessNC::CF::CalendarType NCVariable::get_calendar_type() const {
	// Verify datatype - we require one of the NetCDF classic numeric types,
	// even though that's not a CF requirement. The 64 bit integral types in
	// NetCDF4 would require special care since we can't represent that with
	// a double.
	nc_type type;
	int status = nc_inq_vartype(id_file, id_time, &type);
	GuessNC::handle_error(status, "Failed to get type of time variable");

	if (!(type == NC_BYTE
	 || type == NC_SHORT
	 || type == NC_INT
	 || type == NC_FLOAT
	 || type == NC_DOUBLE)) {
		auto& msg = "Time coordinates must be a NetCDF classic numeric type";
		throw std::runtime_error(msg);
	}

	// Get information about calendar and time unit
	std::string calendar_attribute;
	if (!GuessNC::get_attribute(id_file, id_time, CALENDAR, calendar_attribute)) {
		auto& msg = "Time coordinate variable didn't have a calendar attribute";
		throw std::runtime_error(msg);
	}

	return GuessNC::CF::parse_calendar(calendar_attribute.c_str());
}

GuessNC::CF::TimeUnitSpecification NCVariable::get_time_spec() const {
	std::string units;
	// Get the units on the time variable.std::string attribute;
	if (!GuessNC::get_attribute(id_file, id_time, ATTR_UNITS, units)) {
		throw std::runtime_error("Unable to read units attribute on time variable");
	}

	// Parse units (this will throw if invalid).
	return GuessNC::CF::TimeUnitSpecification(units.c_str());
}

#endif
