///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.h
/// \brief BLAZE fire simulation and combustion
///
/// \author Lars Nieradzik
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_BLAZE_H
#define LPJ_GUESS_BLAZE_H

#include "guess.h"

/// Do daily accounting of BLAZE relevant parameters
/** Accounting of long-term averages needed for BLAZE as well as
 *  the computation of daily burned area and fire-specific parameteers like
 *  the Keetch-Byram Drought-index and Forest Fire Danger Index (FFDI)
 */
void blaze_accounting_gridcell(Climate& climate);	

/// The driver routine for BLAZE
/**This is the driver routine for BLAZE. It retrieves potential Fire-Line-Intensity
 * and calls the blaze main routine patch-wise
 */
void blaze_driver(Patch& patch, Climate& climate);

bool blaze(Patch& patch, Climate& climate, bool force_prescribed_fire_today = false);

double available_fuel (Patch& patch,int fli_index, double k_tun_litter);

extern const double TURNOVERFRACT[13][5];

/* Compute area of Gridcell
 * Taken from aslice.cpp of gutils. 
 * Returns area in square km of a pixel of a given size at a given point
 * on the world. The formula applied is the surface area of a segment of
 * a hemisphere of radius r from the equator to a parallel (circular)
 * plane h vertical units towards the pole: S=2*pi*r*h .
 * latpos    latitude position (see postype)
 * longsize  longitude range in degrees
 * latsize   latitude range in degrees
 * postype   declares which part of the pixel latpos
 *           refer to:
 *           0 = centre
 *           1 = NW corner
 *           2 = NE corner
 *           3 = SW corner
 *           4 = SE corner
 */
double pixelsize(double latpos,double longsize,double latsize,int postype);

#endif 


