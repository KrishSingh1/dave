///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.h
/// \brief BLAZE fire simulation and combustion
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 16:02:51 +0100 (Tue, 24 Jan 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_BLAZE_H
#define LPJ_GUESS_BLAZE_H

#include "guess.h"

void blaze_accounting_gridcell(Climate& climate);	

void blaze_driver(Patch& patch, Climate& climate);

#endif 


