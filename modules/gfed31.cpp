///////////////////////////////////////////////////////////////////////////////////////
/// \file fged31.cpp
/// \brief provide GFED 3.1 data for use in BLAZE
///
/// \author Lars Nieradzik
/// $Date: 2015-05-06 13:19:07 +0200 (Wed, 06 May 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "guess.h"
#include "gfed31.h"

const static int gfed_start_year = 1997;
const static int gfed_end_year = 2011;


/// SIMFIRE biome mapping
double gfed31_ba(Gridcell& gridcell) {
	
	int cyear = date.get_calendar_year(); // current year

	// starting in 01/1997 as gfed month 0 here we get the months
	// Month's index
	int midx  = (cyear - gfed_start_year) * 12 + date.month;

	double ba = 0.0;
	if ( midx >= 0 ) {
		ba = gridcell.monthly_GFED31_ba[midx];
	} 
	else {
		return 0.0;
	}

	if ( blaze_tstep == DAILY ) {
		fail("DAILY TIMESTEP NOT AVAILABLE in gfed31 at the moment!");
		//CLN
		// ba *= daily_GFED31_frac[didx];
	} 
	else {
		ba *= 1./(double)date.ndaymonth[date.month];
	}
	return ba;
}
	
	
