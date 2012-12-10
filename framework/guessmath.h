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
#include <archive.h>

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
double const PI = 4 * atan(1.0);
#else
double const PI = M_PI;
#undef M_PI
#endif
const double DEGTORAD = PI / 180.;

inline bool negligible(double dval) {
	// Returns true if |dval| < EPSILON, otherwise false
	return fabs(dval) < 1.0e-30;
}

inline bool equal(double dval1, double dval2) {
	// Returns true if |dval1-dval2| < EPSILON, otherwise false
	return negligible(dval1 - dval2);
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


/// Keeps track of historic values of some variable
/** Useful for calculating running means etc.
 *
 *  The class behaves like a queue with a fixed size,
 *  when a new value is added, and the queue is full,
 *  the oldest value is overwritten.
 */
template<typename T, int capacity>
class Historic {
public:

	/// The maximum number of elements stored, given as template parameter
	static const int CAPACITY = capacity;

	Historic() 
		: current_index(0), full(false) {
	}

	/// Adds a value, overwriting the oldest if full
	void add(double value) {
		values[current_index] = value;

		current_index = (current_index+1) % CAPACITY;

		if (current_index == 0) {
			full = true;
		}
	}

	/// Returns the number of values stored (0-CAPACITY)
	int size() const {
		return full ? CAPACITY : current_index;
	}

	/// Calculates arithmetic mean of the stored values
	T mean() const {
		const int nvalues = size();

		assert(nvalues != 0);

		return sum()/nvalues;
	}

	/// Sum of stored values
	T sum() const {
		T result = 0.0;

		const int nvalues = size();
		for (int i = 0; i < nvalues; ++i) {
			result += values[i];
		}

		return result;
	}

	/// Returns a single value
	/** The values are ordered by age, the oldest value has index 0.
	 *
	 *  \param pos  Index of value to retrieve (must be less than size())
	 */
	T operator[](size_t pos) const {
		assert(pos < size());

		if (full) {
			return values[(current_index+pos)%CAPACITY];
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
		const int first_position = full ? current_index : 0;
		const int nvalues = size();

		for (int i = 0; i < nvalues; ++i) {
			buffer[i] = values[(first_position+i)%CAPACITY];
		}
	}

private:
	/// The stored values
	T values[CAPACITY];

	/// The next position (in the values array) to write to
	int current_index;

	/// Whether we've stored CAPACITY values yet
	bool full;
};

/// Serialization support for Historic
template<typename T, int capacity>
ArchiveStream& operator&(ArchiveStream& stream,
                         Historic<T, capacity>& data) {
	if (stream.save()) {
		size_t size = data.size();
		stream & size;
		
		for (int i = 0; i < data.size(); ++i) {
			double value = data[i];
			stream & value;
		}
	}
	else {
		size_t size;
		stream & size;

		for (int i = 0; i < size; ++i) {
			double value;
			stream & value;
			data.add(value);
		}
	}
	return stream;
}

#endif // LPJ_GUESS_GUESSMATH_H
