///////////////////////////////////////////////////////////////////////////////////////
/// \file natexptemp.h
/// \brief exp calculations 
///
/// Calculations of exp values 
///
/// \author David Wårlind (based on LPJ-GUESS 2.1 / Ben Smith)
/// $Date: 2012-07-11 18:10:03 +0200 (Wed, 11 Jul 2012) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions
// defined in the module that are to be accessible to the calling framework or
// to other modules.

#ifndef LPJ_GUESS_NATEXPTEMP_H
#define LPJ_GUESS_NATEXPTEMP_H

#include "guessmath.h"
#include <vector>

// Constants required for NatExp lookup tables used by photosynthesis
const double NatExpTemp_MINTEMP = -70;	// minimum temperature ever (deg C)
const double NatExpTemp_MAXTEMP = 70;	// maximum temperature ever (deg C)
const double NatExpTemp_PRECISION = 0.01;	// rounding precision for temperature
const int NatExpTemp_NDATA = static_cast<int>((NatExpTemp_MAXTEMP- NatExpTemp_MINTEMP)/ NatExpTemp_PRECISION + 1.5);
	// maximum number of values to store in each lookup table

/// NatExpTemp lookup table class
/** Stores pre-calculated temperature-adjusted values based on a, c, and d values.
 */
class LookupNatExpTemp {

private:
	/// The temperature-adjusted values
	std::vector<double> data;

public:
	/// Creates a lookup table
	/** \param a  The base value
	 *  \param c  
	 *  \param d  
	 */
	LookupNatExpTemp(double a, double c, double d) : data(NatExpTemp_NDATA) {

		for (int i=0; i< NatExpTemp_NDATA; i++) {
			data[i] = a * exp(c * (NatExpTemp_MINTEMP + i * NatExpTemp_PRECISION - d));
		}
	}

	/// "Array element" operator
	/** \param temp  Temperature (deg C)
	 *  \returns     Temperature-adjusted value based on natural exponetial  
	 */
	double& operator[](double& temp) {
		// Element number corresponding to a particular temperature
		if (temp < NatExpTemp_MINTEMP) {
			temp = NatExpTemp_MINTEMP;
		} else if (temp > NatExpTemp_MAXTEMP) {
			temp = NatExpTemp_MAXTEMP;
		}

		int i = static_cast<int>((temp- NatExpTemp_MINTEMP)/ NatExpTemp_PRECISION+0.5);

		return data[i];
	}

};

#endif // LPJ_GUESS_NATEXPTEMP_H
