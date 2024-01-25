#ifndef LPJ_GUESS_NCUTIL_H
#define LPJ_GUESS_NCUTIL_H

#ifdef HAVE_NETCDF

#include <limits>

#include "guess.h"
#include "ncvariable.h"
#include "util.h"

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
Calculate the sum of a c-style array of data.

@param arr: The array.
@param length: Length of the array.
*/
double sum(const double *arr, size_t length);

/*
Calculate the sum of a vector of doubles.

@param data: The input data.
*/
double sum(const std::vector<double>& data);

/*
Calculate the sum of a c-style array of data.

@param arr: The array.
@param length: Length of the array.
*/
double mean(const double *arr, size_t length);

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
	double min = std::numeric_limits<double>::max();
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
	double max = std::numeric_limits<double>::min();
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

/*
This function aggregates data over a vector in a given chunk size using the
specified aggregator function.

@param in: The input vector.
@param step_size: Desired length of the output vector. An exception will be
				  thrown if length of the input vector is not an integer
				  multiple of this value.
@param aggregator: Function which takes a C-style array and length as input, and
                   aggregates the contents into a single output value.
*/
std::vector<double> aggregate(const std::vector<double>& in, size_t out_length
	, double (*aggregator)(const double *, size_t));

/*
Find the index of the specified element in an array. Return -1 if not found.

@param needle: The element to search for.
@param haystack: The array to be searched.
@param n: Size of the array.
*/
size_t index_of_double(double needle, double *haystack, size_t n);

/*
Find the index of the specified element in an array. Return -1 if not found.

@param needle: The element to search for.
@param haystack: The array to be searched.
@param n: Size of the array.
*/
size_t index_of_float(float needle, float *haystack, size_t n);

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCUTIL_H
