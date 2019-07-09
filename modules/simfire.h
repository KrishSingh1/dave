///////////////////////////////////////////////////////////////////////////////////////
/// \file simfire.cpp
/// \brief SIMFIRE - SIMple FIRE module to compute burnt area  
///
/// \author Lars Nieradzik
/// $Date: 2015-08-25 09:19:28 +0200 (Tue, 25 Aug 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SIMFIRE_H
#define LPJ_GUESS_SIMFIRE_H

#include "guess.h"

// read input data for SIMFIRE
void getsimfiredata(Gridcell& gridcell);
// daily accounting for SIMFIRE
void simfire_accounting_gridcell(Gridcell& gridcell);
// SIMFIRE
double simfire_ba(Climate& climate, Gridcell& gridcell); 
#endif // LPJ_GUESS_SIMFIRE_H
