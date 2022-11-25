#ifdef HAVE_NETCDF

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

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

#ifdef _WIN32
#define mkdir(path, mode) _mkdir((path))
#endif

/* Make a directory; already existing dir okay */
void maybe_mkdir(const char* path, mode_t mode)
{
    errno = 0;

    /* Try to make the directory */
    if (mkdir(path, mode) == 0)
        return;

    /* If it fails for any reason but EEXIST, fail */
    if (errno != EEXIST) {
		char *reason = strerror(errno);
		errno = 0;
		throw_xtring("Unable to create directory %s: %s", path, reason);
	}

    /* Check if the existing path is a directory */
    struct stat st;
    if (stat(path, &st) != 0) {
		char *reason = strerror(errno);
		errno = 0;
		throw_xtring("Unable to stat %s: %s", path, reason);
	}

    /* If not, fail with ENOTDIR */
    if (!S_ISDIR(st.st_mode)) {
		std::string reason = "Path already exists but is not a directory";
        throw_xtring("Unable to create directory %s: %s", path, reason.c_str());
    }
}

void mkdir_p(const char *path)
{
    /*
	Adapted from
	https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
	*/
    mode_t mode = 0777;

    /* Copy string so it's mutable */
    char *_path = strdup(path);
	if (_path == NULL) {
		throw std::runtime_error("Out of memory");
	}

	try {
		// Iterate over the path.
		for (char *p = _path + 1; *p; p++) {
			if (*p == '/') {
				// Temporarily truncate the string.
				*p = '\0';

				// Attempt to create this directory (this can throw).
				maybe_mkdir(_path, mode);

				// Un-truncate the path.
				*p = '/';
			}
		}   

		// This is probably redundant.
		maybe_mkdir(_path, mode);

		free(_path);

	} catch (...) {
		free(_path);
		throw;
	}
}

#endif // HAVE_NETCDF
