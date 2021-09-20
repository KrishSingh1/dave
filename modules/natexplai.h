///////////////////////////////////////////////////////////////////////////////////////
/// \file natexplai.h
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

#ifndef LPJ_GUESS_NATEXPLAI_H
#define LPJ_GUESS_NATEXPLAI_H

#include "guessmath.h"
#include <vector>

// Constants required for NatExpLAI lookup tables used by photosynthesis
const double NatExpLAI_MINLAI = 0;	// minimum LAI
const double NatExpLAI_MAXLAI = 10;	// maximum LAI
const double NatExpLAI_PRECISION = 0.001;	// rounding precision for LAI
const int NatExpLAI_NDATA = static_cast<int>((NatExpLAI_MAXLAI- NatExpLAI_MINLAI)/ NatExpLAI_PRECISION + 1.5);
	// maximum number of values to store in each lookup table

/// NatExpLAI lookup table class
/** Stores pre-calculated LAI-adjusted values based on a, c, and d values.
 */
class LookupNatExpLAI {

private:
	/// The LAI-adjusted values
	std::vector<double> data;

public:
	/// Creates a lookup table
	/** \param a  The base value
	 *  \param c  
	 *  \param d  
	 */
	LookupNatExpLAI(double a, double c) : data(NatExpLAI_NDATA) {

		for (int i=0; i< NatExpLAI_NDATA; i++) {
			data[i] = a * exp(c * (NatExpLAI_MINLAI + i * NatExpLAI_PRECISION));
		}
	}

	/// "Array element" operator
	/** \param lai   LAI
	 *  \returns     LAI-adjusted value based on natural exponetial  
	 */
	double& operator[](double& lai) {
		// Element number corresponding to a particular LAI
		if (lai < NatExpLAI_MINLAI) {
			lai = NatExpLAI_MINLAI;
		} else if (lai > NatExpLAI_MAXLAI) {
			lai = NatExpLAI_MAXLAI;
		}

		int i = static_cast<int>((lai- NatExpLAI_MINLAI)/ NatExpLAI_PRECISION+0.5);

		return data[i];
	}

};

#endif // LPJ_GUESS_NATEXPLAI_H
