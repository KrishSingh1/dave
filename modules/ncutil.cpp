#ifdef HAVE_NETCDF

#include <stdlib.h>

#include "ncutil.h"
#include "ncvariable.h"

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

bool is_start_of_year(timepoint date, int day_length) {
	return date.day == 0 && date.timestep == 0;
}

int get_year_len(int year, bool consider_leaps) {
	return consider_leaps && is_leap(year) ? 366 : 365;
}

bool is_end_of_year(timepoint date, int day_length, bool has_leaps) {
	int days_in_year = get_year_len(date.year, has_leaps);
	// Remember, day is [0, 365) and timestep is [0, N)
	return date.day == days_in_year - 1 == date.timestep == day_length - 1;
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

xtring get_nested_error_msg(const std::exception& e, int level = 0) {
	xtring msg;
	std::string spaces(level, ' ');
	msg.printf("%s%s\n", spaces.c_str(), e.what());
	try {
		std::rethrow_if_nested(e);
	} catch (const std::exception& nested) {
		msg += get_nested_error_msg(nested, level + 1);
	} catch(...) {
		// Ignore any errors which may be thrown while generating nested msg.
	}
	return msg;
}

xtring get_nested_error_msg(const std::exception& e) {
	return get_nested_error_msg(e, 0);
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
