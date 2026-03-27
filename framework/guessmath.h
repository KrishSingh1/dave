////////////////////////////////////////////////////////////////////////////////
/// \file guessmath.h
/// \brief Mathematical utilites header file.
///
/// This header file contains:
///  (1) Definitions of constants and common functions used throughout LPJ-GUESS
///
/// \author Michael Mischurow
/// $Date$
///
////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESSMATH_H
#define LPJ_GUESS_GUESSMATH_H

#include <assert.h>
#include <float.h>
#include "archive.h"
#include "shell.h"

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
double const PI = 4 * atan(1.0);
#else
double const PI = M_PI;
#undef M_PI
#endif
const double DEGTORAD = PI / 180.;

const double SECS_PER_DAY	= 24.0*60.0*60.0;
const double M3_PER_MM3		= 1E-9;
const double MM3_PER_M3		= 1E9;
const double HA_PER_M2		= 1E-4;
const double M2_PER_HA		= 1E4;
const double CM2_PER_M2		= 1E4;
const double MM2_PER_M2		= 1E6;
const double SQ_M			= 1.0;		// 1 m2
const double M_PER_MM		= 1E-3;
const double CM_PER_M		= 1E2;
const double CM_PER_MM		= 1E-1;
const double MM_PER_M		= 1E3;
const double MM_PER_CM		= 1E1;
const double KG_PER_G		= 1E-3;
const double G_PER_KG		= 1E3;
const double G_PER_MG		= 1E-3;
const double MG_PER_G		= 1E3;		// Milligrams
const double KG_PER_MT		= 1E3;		// Metric tonnes
const double MMOL_PER_MOL	= 1E3;
const double UMOL_PER_MOL   = 1E6;
const double J_PER_KJ		= 1E3;
const double KMH_PER_MS         = 3.6;
const double FRACT_TO_PERCENT   = 100.;
const double PERCENT_TO_FRACT   = 0.01;
const double R_EARTH            = 6371.2213; // mean earth-radius[km]

inline bool negligible(double dval, int limit = 0) {
	// If limit is nonzero, returns true if |dval| < 10^limit.
	// Otherwise, returns true if |dval| < EPSILON, otherwise false
	return limit ? fabs(dval) < pow(10.0, limit) : fabs(dval) < 1.0e-30;
}

inline double roundoff(double dval, int limit) {
	// Returns a double where decimals smaller than 10^(-limit) are rounded off.
	// i.e. returns a double where the number of decimals is 'limit'.
	bool neg = dval < 0;

	int ival = int((fabs(dval) * pow(10.0, limit)) + 0.5);
	dval = (double(ival) / pow(10.0, limit));

	if (neg) return -dval; else return dval;
}

inline bool largerthanzero(double dval, int limit = 0) {
	// Returns true if |dval| < EPSILON, otherwise false
	return limit ? dval > pow(10.0, limit) : dval > 1.0e-30;
}

/*
Return true iff x is sufficiently close to 0.

This is an optimised form of equal() which runs faster with a y-value of 0. It's
faster because it skips the check against a relative epsilon, which is fairly
meaningless when comparing against 0.

Do not use this to check if the difference between 2 numbers is 0. This function
will not work perfectly for that purpose, especially if one or both of the
numbers is distant from 0. Use equal() for that purpose.

@param x: The number to check.
@param eps: The tolerance to be used for the comparison. If this really must be
changed, make it a multiple of DBL_EPSILON.
*/
inline bool is_zero(double x, double eps = DBL_EPSILON) {
	return fabs(x) < eps;
}

/*
Returns true if x is sufficiently close to y.

The tolerance increases with the magnitude of the operands, to allow for the
loss in precision in floating-point numbers distant from 0.

@param x: Any real number.
@param y: Any real number.
@param eps: Maximum difference between x and y at which they may be
            considered equal. If the absolute or relative difference is less
			than this number, this function will return true.
*/
inline bool equal(double x, double y, double eps = DBL_EPSILON) {
	double diff = fabs(x - y);

	// This absolute check is important when the numbers are close to 0.
	if (diff <= eps) {
		return true;
	}

	// If the numbers are distant from zero, a relative difference is useful.
	x = fabs(x);
	y = fabs(y);
	double magnitude = (x > y) ? x : y;

	// BLARP! It's best for the tolerance to be a multiple of DBL_EPSILON, but
	// the question is, what multiple? Increasing this scalar will increase the
	// tolerance of the checks. If you're reading this and are tempted to
	// increase SCALAR, it may be better to restructure your calculations in
	// such a way as to reduce your roundoff error.
	const double SCALAR = 8;

	if (diff <= magnitude * eps * SCALAR) {
		return true;
	}

	// The numbers appear to be distant from each other.
	return false;
}

/*
Returns true if x is sufficiently close to y.

The tolerance increases with the magnitude of the operands, to allow for the
loss in precision in floating-point numbers distant from 0.

@param x: Any real number.
@param y: Any real number.
@param eps: Maximum difference between x and y at which they may be
            considered equal. If the absolute or relative difference is less
			than this number, this function will return true.
*/
inline bool fequal(float x, float y, float eps = FLT_EPSILON) {
	float diff = fabs(x - y);

	// This absolute check is important when the numbers are close to 0.
	if (diff <= eps) {
		return true;
	}

	// If the numbers are distant from zero, a relative difference is useful.
	x = abs(x);
	y = abs(y);
	float magnitude = (x > y) ? x : y;

	// BLARP! It's best for the tolerance to be a multiple of DBL_EPSILON, but
	// the question is, what multiple? Increasing this scalar will increase the
	// tolerance of the checks. If you're reading this and are tempted to
	// increase SCALAR, it may be better to restructure your calculations in
	// such a way as to reduce your roundoff error.
	const float SCALAR = 8;

	if (diff <= magnitude * eps * SCALAR) {
		return true;
	}

	// The numbers appear to be distant from each other.
	return false;
}

/*
Return true iff x is greater than or equal to y, within certain tolerance.

@param x: Any real number.
@param y: Any real number.
@param eps: Maximum absolute difference between x and y at which they may be
            considered equal.
*/
inline bool greater_than_or_equal(double x, double y, double eps = DBL_EPSILON) {
	return x > y || equal(x, y, eps);
}

/*
Return true iff x is greater than or equal to 0, within certain tolerance.

@param x: Any real number.
@param y: Any real number.
@param eps: Maximum absolute difference between x and y at which they may be
            considered equal.
*/
inline bool positive(double x, double eps = DBL_EPSILON) {
	return x > 0 || is_zero(x, eps);
}

/*
Return true iff x is less than or equal to y, within certain tolerance.

@param x: Any real number.
@param y: Any real number.
@param eps: Maximum absolute difference between x and y at which they may be
            considered equal.
*/
inline bool less_than_or_equal(double x, double y, double eps = DBL_EPSILON) {
	return x < y || equal(x, y, eps);
}

inline double mean(double* array, int nitem) {

	// Returns arithmetic mean of 'nitem' values in 'array'

	double sum=0.0;
	for (int i=0; i<nitem; sum += array[i++]);
	return sum / (double)nitem;
}

/// Gives the mean of just two values
inline double mean(double x, double y) {
	return (x+y)/2.0;
}

// Divide x by y. Raise an error if denominator is 0.
#define divide(x, y) safe_divide((x), (y), __FILE__, __LINE__, __func__)

inline double safe_divide(double x, double y,
                          const char *file, int line, const char *func) {
	if (is_zero(y, 0)) {
		fail("Divide by 0: %f / %f at %s:%d (%s)", x, y, file, line, func);
	}
	return x / y;
}

/// Calculates modulus
inline size_t modulo(int a, int b) {
	return ((a % b) + b) % b;
}

/// Calculates variation coefficient of values in an array
inline double variation_coefficient(double* data, int n) {
	// 0 and 1 will give division with zero.
	if (n <= 1) {
		return -1;
	}

	double avg = mean(data, n);
	double dev = 0;
	for (int i=0; i<n; i++) {
		dev += (data[i] - avg) * (data[i] - avg);
	}
	double std = sqrt(dev / (n-1));

	if (std > 0 && avg > 0) {	// check that data appear in the array
		return std / avg;
	}
	return 0;
}

/// A short version of Richards curve where:
/** a is the lower asymptote,
 *  b is the upper asymptote. If a=0 then b is called the carrying capacity,
 *  c the growth rate,
 *  d is the time of maximum growth
 *  Source: https://en.wikipedia.org/wiki/Generalised_logistic_function, 2013-11-11
 */
inline double richards_curve(double a, double b, double c, double d, double x) {
	return a + (b - a) / (1 + exp(-c * (x - d)));
}

inline void regress(double* x, double* y, int n, double& a, double& b) {

	// Performs a linear regression of array y on array x (n values)
	// returning parameters a and b in the fitted model: y=a+bx
	// Source: Press et al 1986, Sect 14.2

	double sx,sy,sxx,sxy,delta;
	sx = sy = sxy = sxx = 0.0;

	for (int i=0; i<n; i++) {
		sx += x[i];
		sy += y[i];
		sxx+= x[i]*x[i];
		sxy+= x[i]*y[i];
	}
	delta = (double)n*sxx - sx*sx;
	a = (sxx*sy - sx*sxy)/delta;
	b = ((double)n*sxy-sx*sy)/delta;
}

// Forward declarations needed for the friendship between Historic and operator&

template<typename T, size_t capacity>
class Historic;

template<typename T, size_t capacity>
ArchiveStream& operator&(ArchiveStream& stream,
                         Historic<T, capacity>& data);


/// Keeps track of historic values of some variable
/** Useful for calculating running means etc.
 *
 *  The class behaves like a queue with a fixed size,
 *  when a new value is added, and the queue is full,
 *  the oldest value is overwritten.
 */
template<typename T, size_t capacity>
class Historic {
public:

	/// The maximum number of elements stored, given as template parameter
	static const size_t CAPACITY = capacity;

	Historic()
		: current_index(0), full(false) {
		
		for (size_t i = 0; i < capacity; ++i) {
			values[i] = 0.0;
		}
	}

	/// Adds a value, overwriting the oldest if full
	void add(double value) {
		values[current_index] = value;

		current_index = (current_index+1) % capacity;

		if (current_index == 0) {
			full = true;
		}
	}

	/// Returns the number of values stored (0-CAPACITY)
	size_t size() const {
		return full ? capacity : current_index;
	}

	/// Calculates arithmetic mean of the stored values
	T mean() const {
		const size_t nvalues = size();

		assert(nvalues != 0);

		return sum()/nvalues;
	}

	/// Sum of stored values
	T sum() const {
		T result = 0.0;

		const size_t nvalues = size();
		for (size_t i = 0; i < nvalues; ++i) {
			result += values[i];
		}
		return result;
	}

	/// Returnes the latest of the stored values
	T lastadd() const {	//latest
		return (*this)[size()-1];
	}

	/// Returns the maximum of the stored values
	T max() const {
		T result = -999999.0;

		const size_t nvalues = size();
		for (size_t i = 0; i < nvalues; ++i) {
			if (values[i]>result) {
				result = values[i];
			}
		}
		return result;
	}

	/// Returns the minimum of the stored values
	T min() const {
		T result = 999999.0;

		const size_t nvalues = size();
		for (size_t i = 0; i < nvalues; ++i) {
			if (values[i]<result) {
				result = values[i];
			}
		}
		return result;
	}

	/// Calculates arithmetic mean of the stored values for a period from current position and backwards nsteps
	T periodicmean(const size_t nsteps) const {

		assert(nsteps != 0);
		if (nsteps > size()){
			return mean();
		} else {
			return periodicsum(nsteps)/nsteps;
		}
	}

	/// Sum of stored values for a period from current position and backwards nsteps
	T periodicsum(const size_t nsteps) const {
		T result = 0.0;

		if (nsteps >= size()){
			return sum();
		} else {
			for (size_t i = size()-1; i >=size()-nsteps; --i) {
				result += (*this)[i];
			}
		}

		return result;
	}

	/// Periodic sum of differences 
    T periodictrend(size_t nsteps) const {
        size_t current_size = size();
        if (current_size == 0) {
            return 0;
        }
 
        nsteps = std::min(nsteps, current_size - 1);
 
        T trend = 0.0;
        for (size_t n = 0; n < nsteps; n++) {
            size_t index = modulo(current_index - (n + 1), current_size);
            size_t index_prior = modulo(current_index - (n + 2), current_size);
            // todo: use modulus here (much faster but need to test it first)
            T current = values[index];
            T previous = values[index_prior];
            trend += (current - previous);
        }
        return trend;
    }

	/// Trend/Slope of a linear regression of values
	T trend() const {
		T result = 0.0;

		if (size() == 1)
			return result;
		else {

			const size_t nvalues = size();
			std::vector<double> x(nvalues);
			std::vector<double> y(nvalues);
			double x_avr = 0.0;
			double y_avr = 0.0;

			for (size_t i = 0; i < nvalues; ++i) {

				int p = (current_index - nvalues + i) % capacity;
				x[i] = i;
				y[i] = values[p];
				x_avr += i / (double)nvalues;
				y_avr += values[p] / (double)nvalues;
			}

			double numerator = 0.0;
			double denominator = 0.0;

			for (size_t i = 0; i < nvalues; ++i) {
				numerator += (x[i] - x_avr) * (y[i] - y_avr);
				denominator += (x[i] - x_avr) * (x[i] - x_avr);
			}

			result = numerator / denominator;

			return result;
		}
	}

	/// Returns a single value
	/** The values are ordered by age, the oldest value has index 0.
	 *
	 *  \param pos  Index of value to retrieve (must be less than size())
	 */
	T operator[](size_t pos) const {
		assert(pos < size());

		if (full) {
			return values[(current_index+pos)%capacity];
		}
		else {
			return values[pos];
		}
	}

	/// Writes all values to a plain buffer
	/** The values will be ordered by age, the oldest value has index 0.
	 *
	 *  \param buffer   Array to write to, must have room for at least size() values
	 */
	void to_array(T* buffer) const {
		const size_t first_position = full ? current_index : 0;
		const size_t nvalues = size();

		for (size_t i = 0; i < nvalues; ++i) {
			buffer[i] = values[(first_position+i)%capacity];
		}
	}

	friend ArchiveStream& operator&<T, capacity>(ArchiveStream& stream,
	                                             Historic<T, capacity>& data);

private:
	/// The stored values
	T values[capacity];

	/// The next position (in the values array) to write to
	size_t current_index;

	/// Whether we've stored CAPACITY values yet
	bool full;
};


/// New Historic class with the possibility to define 
// length as input.
template<typename T>
class HistoricDynamic {
public:
	HistoricDynamic() {
		capacity = 0;
		values = NULL;
		full = false;
		current_index = 0;
	}

	// Set the capacity to the specified value.
	void set_capacity(size_t capacity) {
		if (capacity == this->capacity) {
			return;
		}
		// TODO Drew: double check this logic.
		T *new_list = (T *)malloc(sizeof(T) * capacity);
		size_t n = size();
		if (capacity < this->capacity) {
			size_t start = modulo(current_index - std::min(n,capacity), this->capacity);
			for (size_t i = 0; i < capacity; i++) {
				new_list[i] = values[(start + i) % this->capacity];
			}
		} else {
			for (size_t i = 0; i < capacity; i++) {
				new_list[i] = i < n ? values[i] : 0;
			}
		}
		if (values) {
			free(values);
		}
		values = new_list;
		this->capacity = capacity;
		// Current index should wrap around to the start of the list if
		// necessary.
		current_index = current_index % capacity;
	}

	size_t get_capacity() const {
		return capacity;
	}

	/// Adds a value, overwriting the oldest if full
	void add(double value) {
		assert(capacity > 0);

		values[current_index] = value;

		current_index = (current_index + 1) % capacity;

		if (current_index == 0) {
			full = true;
		}
	}

	/// Get the mean value stored in the vector.
	double mean() const {
		// Get the number of values currently stored, in case the vector is not
		// yet full.
		size_t n = size();

		double sum = 0;
		for (size_t i = 0; i < n; i++) {
			sum += values[i];
		}
		return sum / n;
	}

	/// Get the mean value stored in the vector.
	double weighted_mean() const {
		// Get the number of values currently stored, in case the vector is not
		// yet full.
		size_t n = size();

		double sum = 0;
		double total_weight = 0;
		for (size_t i = 0; i < n; i++) {
			double weight = 1.0 - i / n;
			sum += values[i] * weight;
			total_weight += weight;
		}
		return sum / (n * total_weight);
	}

	/// Get the number of values currently stored.
	size_t size() const {
		return full ? capacity : current_index;
	}

private:
	/// The number of values this object can currently hold.
	size_t capacity;

	/// The stored values.
	T *values;

	/// The next position (in the values array) to write to.
	size_t current_index;

	/// Whether we've reached the maximum number of values.
	bool full;
};


/// Serialization support for Historic
/** We have a friend serialization operator instead of
 *  implementing Serializable, to avoid overhead of a
 *  vtable in Historic (since serialize() is virtual).
 *
 *  This could be implemented without friend, by only
 *  using size(), operator[] and add(). The order of
 *  the elements in the buffer will however not be
 *  preserved then, and for instance the sum() function
 *  will then not give _exactly_ the same results
 *  (due to limited floating point precision).
 */
template<typename T, size_t capacity>
ArchiveStream& operator&(ArchiveStream& stream,
                         Historic<T, capacity>& data) {
	stream & data.values
		& data.current_index
		& data.full;

	return stream;
}


// used in soil temperature calculation
inline bool verysmall(double dval) {

	// Returns true if dval < EPSILON, otherwise false

	const double EPSILON = 1.0e-12;
	if (dval>EPSILON) return false;
	if (dval<0.0 && dval<-EPSILON) return false;
	return true;
}

#endif // LPJ_GUESS_GUESSMATH_H
