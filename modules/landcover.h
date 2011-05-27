///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H

#include "guess.h"

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell,Pftlist& pftlist);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell,Pftlist& pftlist);

#endif // LPJ_GUESS_LANDCOVER_H
