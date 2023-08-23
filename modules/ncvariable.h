///////////////////////////////////////////////////////////////////////////////////////
/// \file ncvariable.h
/// \brief Class which encapsulates a variable in a NetCDF file.
///
/// \author Drew Holzworth
/// $Date: 2022-11-16 10:08:14 +1100 (Tue, 16 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_NCVARIABLE_H
#define LPJ_GUESS_NCVARIABLE_H

#ifdef HAVE_NETCDF

#include <stdint.h>
#include <vector>

#include "cftime.h"
#include "gutil.h"
#include "ncutil.h"

/*
Represents a single time point for a single gridcell.
*/
struct timepoint {
	/*
	Gregorian Year (not necessarily an index).
	*/
	int year;

	/*
	Day of year (in range [0, 365-ish]).
	*/
	int day;

	/*
	Subdaily timestep (in range [0, day_length]).
	*/
	int timestep;
};

/*
Represents upper/lower bounds for a variable.
*/
struct bounds {
	/*
	Smallest allowed value of the variable.
	*/
	double min;

	/*
	Largest allowed value of the variable.
	*/
	double max;
};

/*
This enum defines the cache strategy used by this variable.
*/
enum CacheStrategy {
	// Cache the entire timeseries. Fastest but most memory intensive.
	All,
	// Cache one year of data at a time. Reduces memory usage.
	Year,
	// No caching. Read from file every time. Not recommended.
	None
};

/*
This class encapsulates a variable in a NetCDF file. It provides an easy API
for reading data and all caching is handled internally.
*/
class NCVariable {
public:
	/*
	Initialise a new NCVariable instance, perform basic sanity checking.

	@param file_id: NetCDF file ID.
	@param name: Name of the variable.
	@param dim_lon: Name of the longitude dimension.
	@param dim_lat: Name of the latitude dimension.
	@param has_leap: True iff the input data contains leap days.
	@param dim_time: Name of the time dimension.
	@param sensibility_checks: True if bounds checking should occur.
	@param limits: Upper/lower bounds of valid data.
	@param units: Expected units.
	@param check_units: Iff true, an error will be generated for invalid units.
	*/
	NCVariable(int file_id, xtring name, xtring dim_lon, xtring dim_lat
	, xtring dim_time, bool has_leap, CacheStrategy chunk_mode
	, bool sensibility_checks, bounds limits, std::vector<xtring> units
	, bool check_units);

	/*
	Read one day's data for a particular gridpoint. Note that this method may
	not necessarily read data from the file (internal cache may be used).

	@param lon: Longitude of the data to be read.
	@param lat: Latitude of the data to be read.
	@param year: Start year of data to be read.
	@param day: Start day of year of data to be read.
	@param nday: Number of days' data to be read.
	*/
	std::vector<double> read(double lon, double lat, int year, int day
	, int nday);

	/*
	Get the standard name of the variable, as stored in the "standard_name"
	attribute.
	*/
	xtring get_standard_name() const;

	/*
	Get the date/time of the first timestep in the input data.
	*/
	timepoint get_first_date() const;

	/*
	Get the date/time of the final timestep in the input data.
	*/
	timepoint get_last_date() const;

	/*
	Get the number of timesteps per day in this NetCDF file. Note that this does
	not calculate the value.
	*/
	int get_day_length() const;

private:
	/*
	Get the length of the specified dimension.

	@param dim_id: Dimension ID.
	*/
	size_t get_dimension_length(int dim_id) const;

	/*
	Get the name of the specified dimension.
	*/
	xtring get_dimension_name(int dim_id) const;

	/*
	Check if the currently cached data is valid for a particular day of the
	specified grid point.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param year: Year of data to be read.
	@param day: Day of year of data to be read.
	*/
	bool is_cache_valid(double lon, double lat, int year, int day) const;

	/*
	Read data for the specified day from the cache.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param year: Year of data to be read.
	@param day: Day of year of data to be read.
	*/
	std::vector<double> read_from_cache(double lon, double lat, int year
	, int day) const;

	/*
	Load data from the file using the specified caching strategy for a
	particular day of the specified gridpoint.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param year: Year of data to be read.
	@param day: Day of year of data to be read.
	*/
	std::vector<double> load(double lon, double lat, int year, int day) const;

	/*
	Load one day of data from the input file for a particular grid point.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param year: Year of data to be read.
	@param day: Day of year of data to be read.
	*/
	std::vector<double> load_day(double lon, double lat, int year, int day) const;

	/*
	Load one year of data from the input file for a particular grid point.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param year: Year of data to be read.
	*/
	std::vector<double> load_year(double lon, double lat, int year) const;

	/*
	Load all data from the input file for a particular grid point.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	*/
	std::vector<double> load_all(double lon, double lat) const;

	/*
	Load data from the file for a particular gridpoint.

	@param lon: Longitude of data to be read.
	@param lat: Latitude of data to be read.
	@param time_start: Start index along time dimension of data to be read.
	@param ntimesteps: Number of timesteps to be read.
	*/
	std::vector<double> load_data(double lon, double lat, int time_start
	, int ntimesteps) const;

	/*
	Get the number of timesteps between the specified year/days.
	This is variable depending on day_length and has_leap_years.

	The first date must be less (ie earlier) than the second date.

	@param date0: The first date (must be earlier than the second).
	@param date1: The second date (must be later than the first).
	*/
	int get_timesteps_between(timepoint date0, timepoint date1) const;

	/*
	Calculate the offset of the specified date in the cache.

	@param date: A date.
	*/
	int get_cache_offset(timepoint date) const;

	/*
	Get the ID of the variable with the specified name. Throw if it doesn't
	exist.

	@param var_name: Name of the variable.
	*/
	int get_variable_id(xtring var_name) const;

	/*
	Get the number of dimensions for the specified variable.

	@param var_id: ID of the variable.
	*/
	int get_num_dimensions(int var_id) const;

	/*
	Check the number of dimensions for the specified variable and throw an
	exception if they are not equal to the specified number.

	@param vad_id: ID of the variable.
	@param ndim_expected: Expected number of dimensions.
	*/
	void check_num_dims(int var_id, int ndim_expected) const;

	/*
	Get the units of this variable.
	*/
	xtring get_units() const;

	/*
	Ensure that the units of the variable match the expected units. A number of
	"synonyms" can be passed in via the expected units argument, to allow for a
	certain degree of flexibility in naming (e.g. m/s vs ms^-1). If the units do
	not match the expected, an error or warning will be emitted, depending on
	the value of the fatal argument.

	@param expected: Expected units (and synonyms thereof).
	@param fatal: If true, a fatal error will be triggered by invalid units. If
	              false, a warning will be written instead.
	*/
	void check_units(const std::vector<xtring>& expected, bool fatal) const;

	/*
	Return the index of the given value in the specified dimension

	@param x: The variable value.
	@param var_id: ID of the variable. Will throw if this is invalid.
	*/
	int index_of(double x, int dim_id) const;

	/*
	Convert a GuessNC::CF::DateTime object to a timepoint.

	@param date: The date object to be converted.
	*/
	timepoint to_timepoint(const GuessNC::CF::DateTime& date) const;

	/*
	Calculate the number of timesteps per day based on the contents of the
	time variable. Returns -1 if file only contains a single day (or less).
	*/
	int calculate_day_length() const;

	/*
	Get the calendar type used by this NetCDF.
	*/
	GuessNC::CF::CalendarType get_calendar_type() const;

	/*
	Get the units of the time variable used in this file.
	*/
	GuessNC::CF::TimeUnitSpecification get_time_spec() const;

	/*
	Get the nth date/time step in the NetCDF file.
	*/
	GuessNC::CF::DateTime get_nth_date(size_t n) const;

	/*
	NetCDF file ID.
	*/
	int id_file;

	/*
	Variable ID.
	*/
	int id_var;

	/*
	Number of dimensions of this variable.
	*/
	int ndim;

	/*
	ID of the longitude variable.
	*/
	int id_lon;

	/*
	Index of the longitude dimension in this variable's dimensions.
	*/
	int dim_lon;

	/*
	ID of the latitude variable.
	*/
	int id_lat;

	/*
	Index of the latitude dimension in this variable's dimensions.
	*/
	int dim_lat;

	/*
	ID of the time variable.
	*/
	int id_time;

	/*
	Total number of timesteps in the file.
	*/
	size_t ntime;

	/*
	Index of the time dimension in this variable's dimensions.
	*/
	int dim_time;

	/*
	Variable name.
	*/
	xtring var_name;

	/*
	Number of timesteps per day.
	*/
	int day_length;

	/*
	True if input data contains leap years. False otherwise. We assume that the
	input data never contains other anomalies such as leap seconds.
	*/
	bool has_leap_years;

	/*
	Defines how much data should be cached.
	*/
	CacheStrategy cache_strategy;

	/*
	Cached data.
	*/
	std::vector<double> cache;

	/*
	Longitude of cached data.
	*/
	double longitude;

	/*
	Latitude of cached data.
	*/
	double latitude;

	/*
	The starting date of the cached data.
	*/
	timepoint cache_start;

	/*
	The previously-read timepoint of the cached data. This is used as an
	optimisation, to avoid needless seek operations.
	*/
	timepoint cache_prev;

	/*
	The index in the cached data of cache_prev.
	*/
	int cache_prev_index;

	/*
	The starting date of the data in the input file.
	*/
	timepoint data_start;

	/*
	The calendar type used by this variable.
	*/
	GuessNC::CF::CalendarType calendar_type;

	/*
	The way that time data is represented in the file (ie time units).
	*/
	GuessNC::CF::TimeUnitSpecification time_spec;

	/*
	Upper and lower bounds for data.
	*/
	bounds limits;

	/*
	True if bounds checking should occur.
	*/
	bool check_bounds;

	/*
	Name of the units attribute on the time variable. This is part of the CF
	spec.
	*/
	const std::string ATTR_UNITS = "units";

	/*
	Name of the calendar attribute required on the time variable. This should be
	of the form "X since yyyy-MM-dd", where X may be days, seconds, years, etc.
	*/
	const std::string CALENDAR = "calendar";
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCVARIABLE_H