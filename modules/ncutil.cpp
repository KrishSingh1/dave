#ifdef HAVE_NETCDF

#include <time.h>
#include <sys/time.h>

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

double sum(const std::vector<double>& data) {
	double sum = 0;
	for (int i = 0; i < data.size(); i++) {
		sum += data[i];
	}
	return sum;
}

double mean(const std::vector<double>& data) {
	return sum(data) / data.size();
}

bool floats_equal(double x, double y) {
	return negligible(abs(x - y));
}

double round_to_nearest_multiple(double x, double y) {
	if (floats_equal(x, 0) || floats_equal(y, 0)) {
		return 0;
	}

	double ratio = y / x;
	double upper = y * ceil(ratio);
	double lower = y * floor(ratio);
	double midpoint = (upper + lower) / 2;

	if (x >= midpoint) {
		return upper;
	}
	return lower;
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

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

#endif // HAVE_NETCDF
