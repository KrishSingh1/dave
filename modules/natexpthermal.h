///////////////////////////////////////////////////////////////////////////////////////
/// \file natexpthermal.h
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

#ifndef LPJ_GUESS_NATEXPTHERMAL_H
#define LPJ_GUESS_NATEXPTHERMAL_H

#include "guessmath.h"
#include "soil.h"
#include <vector>

// Constants required for NatExpThermal lookup tables used by soillayer
const double NatExpThermal_MINLogK = Kpeat_log;	// minimum log(K)
const double NatExpThermal_MAXLogK = Kice_log;	// maximum log(K)
const double NatExpThermal_PRECISION = 0.001;	// rounding precision for LogK
const int NatExpThermal_NDATA = static_cast<int>((NatExpThermal_MAXLogK- NatExpThermal_MINLogK)/ NatExpThermal_PRECISION + 1.5);
	// maximum number of values to store in each lookup table

/// NatExpThermal lookup table class
/** Stores pre-calculated LogK-adjusted values.
 */
class LookupNatExpThermal {

private:
	/// The Thermal-adjusted values
	std::vector<double> data;

public:
	/// Creates a lookup table
	LookupNatExpThermal(double a, double b) : data(NatExpThermal_NDATA) {

		for (int i = 0; i < NatExpThermal_NDATA; i++) {
			data[i] = a * exp(b * (NatExpThermal_MINLogK + i * NatExpThermal_PRECISION));
		}
	}

	/// "Array element" operator
	/** \param logK  Log(K)
	 *  \returns     Log(K)-adjusted value based on natural exponetial  
	 */
	double& operator[](double& logK) {
		// Element number corresponding to a particular Log(K)
		if (logK < NatExpThermal_MINLogK) {
			logK = NatExpThermal_MINLogK;
		}
		else if (logK > NatExpThermal_MAXLogK) {
			logK = NatExpThermal_MAXLogK;
		}

		int i = static_cast<int>((logK - NatExpThermal_MINLogK) / NatExpThermal_PRECISION + 0.5);

		return data[i];
	}

};

#endif // LPJ_GUESS_NATEXPTHERMAL_H
