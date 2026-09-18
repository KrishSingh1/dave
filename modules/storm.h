////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file storm.h
/// \Calculates storm sensitivity and damaged fraction of individuals, which is also scaled to stand and gridcell
/// \Equations are taken from Lagergren et al 2012, Ecological Modelling 247:71-82
/// \author Fredrik Lagergren
/// $Date: 2018-05-24 08:03:00 +0200 (Thu, 24 May 2018) $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_STORM_H
#define LPJ_GUESS_STORM_H

#include "guess.h"

void storm_init(Gridcell& gridcell);
void storm(Gridcell& gridcell);

#endif // LPJ_GUESS_STORM_H

