#ifdef HAVE_NETCDF

#include <stdlib.h>

#include "ncutil.h"
#include "ncvariable.h"

const int GREGORIAN_START = 1582;

insoltype parse_insoltype(const std::string& standard_name) {
	if (standard_name == "surface_downwelling_shortwave_flux_in_air" ||
	    standard_name == "surface_downwelling_shortwave_flux") {
		return SWRAD_TS;
	}
	else if (standard_name == "surface_net_downward_shortwave_flux") {
		return NETSWRAD_TS;
	}
	else if (standard_name == "cloud_area_fraction") {
		return SUNSHINE;
	}
	else {
		fail("Unknown insolation type: %s", standard_name.c_str());
		return SUNSHINE; // To avoid compiler warning
	}
}

double sum(const double *arr, size_t length) {
	double sum = 0;
	for (size_t i = 0; i < length; i++) {
		sum += arr[i];
	}
	return sum;
}

double mean(const double *arr, size_t length) {
	return sum(arr, length) / length;
}

double sum(const std::vector<double>& data) {
	return sum(&(data.front()), data.size());
}

double mean(const std::vector<double>& data) {
	return mean(&(data.front()), data.size());
}

bool floats_equal(double x, double y) {
	return negligible(abs(x - y));
}

size_t index_of_double(double needle, double *haystack, size_t n) {
	for (size_t i = 0; i < n; i++) {
		if (equal(needle, haystack[i])) {
			return i;
		}
	}
	return -1;
}

size_t index_of_float(float needle, float *haystack, size_t n) {
	for (size_t i = 0; i < n; i++) {
		if (fequal(needle, haystack[i])) {
			return i;
		}
	}
	return -1;
}

double round_to_nearest_multiple(double x, double y) {
	// 9, 24 => 0
	// 9, 3  => 3
	// 10, 3 => 3
	// 11, 3 => 4
	if (floats_equal(x, 0) || floats_equal(y, 0)) {
		return 0;
	}

	double z = x / y;
	return (int)(z < 0 ? z - 0.5 : z + 0.5);
}

bool contains(double x, const std::vector<double> &vec) {
	for (const double &xi : vec) {
		if (x == xi) {
			return true;
		}
	}
	return false;
}

bool is_start_of_year(timepoint date) {
	return date.day == 0 && date.timestep == 0;
}

/*
Check if this year would be a leap year according to the specified CF calendar.
*/
bool is_leap_year(int year, GuessNC::CF::CalendarType calendar) {
	if (calendar == GuessNC::CF::CalendarType::NO_LEAP) {
		// No leap years in this calendar.
		return false;
	}

	if (calendar == GuessNC::CF::CalendarType::STANDARD && year < GREGORIAN_START) {
		// Gregorian calendar, but this year is before the Gregorian calendar
		// was introduced, so this year has 365 days.
		return false;
	}

	// This is a gregorian calendar year, so we need to check if it would be
	// considered a leap year.
	return is_leap(year);
}

int get_year_len(int year, GuessNC::CF::CalendarType calendar) {
	return is_leap_year(year, calendar) ? 366 : 365;
}

int get_month_len(int year, int month, GuessNC::CF::CalendarType calendar) {

	const double days_in_month[12] =
		{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	if (month == 1 && is_leap_year(year, calendar)) {
		// February.
		return 29;
	}

	return days_in_month[month];
}

int get_month(timepoint timestamp, GuessNC::CF::CalendarType calendar) {
	int day = 0;
	for (int mon = 0; mon < 12; mon++) {
		day += get_month_len(timestamp.year, mon, calendar);
		if (day > timestamp.day) {
			return mon;
		}
	}
	fail("get_month(): Invalid timestamp: %04d-%02d", timestamp.year, timestamp.day);
	return -1; // Cannot happen; including this to silence some compilers.
}

int get_day_of_month(timepoint timestamp, GuessNC::CF::CalendarType calendar) {
	int day = 0;
	for (int mon = 0; mon < 12; mon++) {
		int len = get_month_len(timestamp.year, mon, calendar);
		if ( (day + len) > timestamp.day) {
			return timestamp.day - day;
		}
		day += len;
	}
	fail("get_month(): Invalid timestamp: %04d-%02d", timestamp.year, timestamp.day);
	return -1; // Cannot happen; including this to silence some compilers.
}

bool is_end_of_year(timepoint date, int day_length, GuessNC::CF::CalendarType calendar) {
	int days_in_year = get_year_len(date.year, calendar);
	// Remember, day is [0, 365) and timestep is [0, N)
	return date.day == days_in_year - 1 && date.timestep == day_length - 1;
}

bool is_leap(int year) {
	if (year % 4 != 0) {
		return false;
	}
	if (year % 100 != 0) {
		return true;
	}
	if (year % 400 != 0) {
		return false;
	}
	return true;
}

void bounds_check(const bounds& limits, const std::vector<double>& data) {
	for (const double& value : data) {
		if (value < limits.min) {
			throw_xtring("Value %.2f exceeds lower bound of %.2f", value
				, limits.min);
		}
		if (value > limits.max) {
			throw_xtring("Value %.2f exceeds upper bound of %.2f", value
				, limits.max);
		}
	}
}

xtring rep(const char *c, int n) {
	xtring res(n);
	res += c;
	return res;
}

std::vector<double> aggregate(const std::vector<double>& in, size_t out_length
	, double (*aggregator)(const double *, size_t)) {

	if (in.size() % out_length != 0) {
		throw_xtring("Unable to aggregate data; input size (%d) must be an integer multiple of out size (%d)"
		, in.size(), out_length);
	}
	int step_size = in.size() / out_length;

	std::vector<double> out(out_length);
	for (size_t i = 0; i < out_length; i++) {
		out[i] = aggregator(&(in[i * step_size]), step_size);
	}

	return out;
}

bool exists(const char *filename) {
	if (FILE *file = fopen(filename, "r")) {
		fclose(file);
		return true;
	}
	return false;
}

#endif // HAVE_NETCDF
