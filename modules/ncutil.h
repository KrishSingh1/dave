#ifndef LPJ_GUESS_NCUTIL_H
#define LPJ_GUESS_NCUTIL_H

#ifdef HAVE_NETCDF

#include <limits.h>

#include "guess.h"
#include "ncvariable.h"

/*
Throw an exception using the specified format string and arguments.

@param fmt: printf-style format string.
@param ...: Arguments to the format string.
*/
#define throw_xtring(fmt, ...) \
	xtring tmp_msg; \
	tmp_msg.printf((fmt), __VA_ARGS__); \
	throw std::runtime_error((const char *)tmp_msg);

/*
Throw a nested exception using the specified format string and arguments.

@param fmt: printf-style format string.
@param ...: Arguments to the format string.
*/
#define throw_xtring_nested(fmt, ...) \
	xtring tmp_msg; \
	tmp_msg.printf((fmt), __VA_ARGS__); \
	std::throw_with_nested(std::runtime_error((const char *)tmp_msg));

/*
This macro will check the status code returned by functions in the NetCDF API,
and throw an exception if the status code does not indicate success. The thrown
exception will contain the specified error message.

This macro exists so that the string concatenation does not occur unless the
status code actually indicates an error.

@param s: The status code (an integer).
@param fmt: printf-style format string.
@param ...: Arguments to the format string.
*/
#define handle_nc_error(s, fmt, ...) if ((s) != NC_NOERR) { \
	xtring tmp_msg; \
	tmp_msg.printf((fmt), __VA_ARGS__); \
	GuessNC::handle_error((s), (const char *)tmp_msg); \
}

/*
Number of seconds in a minute.
*/
const int SECONDS_PER_MINUTE = 60;

/*
Number of minutes in an hour.
*/
const int MINUTES_PER_HOUR = 60;

/*
Number of seconds in an hour.
*/
const int SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR;

/*
Number of hours in a day.
*/
const int HOURS_PER_DAY = 24;

/*
Converts a string to one of our known insolation types, or calls fail() if the
name is invalid.

@param standard_name: The input string.
*/
insoltype parse_insoltype(const std::string& standard_name);

/*
Calculate the sum of a vector of doubles.

@param data: The input data.
*/
double sum(const std::vector<double>& data);

/*
Calculate the mean of a vector of doubles.

@param data: The input data.
*/
double mean(const std::vector<double>& data);

/*
Return the smaller of the two numbers.

@param x: Any real number.
@param y: Any real number.
*/
inline double min(double x, double y) {
	return x < y ? x : y;
}

/*
Return the larger of the two numbers.

@param x: Any real number.
@param y: Any real number.
*/
inline double max(double x, double y) {
	return x > y ? x : y;
}

/*
Return the smallest element in a vector.

@param values: The vector.
*/
inline double min(std::vector<double> values) {
	double min = __DBL_MAX__;
	for (const double& value : values) {
		if (value < min) {
			min = value;
		}
	}
	return min;
}

/*
Return the largest element in a vector.

@param values: The vector.
*/
inline double max(std::vector<double> values) {
	double max = __DBL_MIN__;
	for (const double& value : values) {
		if (value > max) {
			max = value;
		}
	}
	return max;
}

/*
Check if two floating point numbers are equal (or very close to) by performing
a float-friendly equality check.
*/
bool floats_equal(double x, double y);

/*
Round x to the nearest multiple of y.

@param x: The value to be rounded.
@param y: The roundee.
*/
double round_to_nearest_multiple(double x, double y);

/*
Check if the given timepoint represents the first time step in the year.

@param date: The date.
@param day_length: Number of timesteps per day.
*/
bool is_start_of_year(timepoint date, int day_length);

/*
Check if the given timepoint represents the first timestep in the year.

@param date: The date.
@param day_length: Number of timesteps per day.
@param has_leaps: True iff leap years should be considered.
*/
bool is_end_of_year(timepoint date, int day_length, bool has_leaps);

/*
Check if the specified year is a leap year.

@param year: The year.
*/
bool is_leap(int year);

/*
Get the number of days in the specified year.

@param year: The year.
@param consider_leaps: True to consider leap years. False otherwise.
*/
int get_year_len(int year, bool consider_leaps);

/*
Bounds check all values in the vector.

@param limits: Upper/lower bounds for the data.
@param data: The data to be checked.
*/
void bounds_check(const bounds& limits, const std::vector<double>& data);

/*
Get the nested error message for the given exception.

@param e: The exception.
*/
xtring get_nested_error_msg(const std::exception& e);

double get_wall_time();

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCUTIL_H
