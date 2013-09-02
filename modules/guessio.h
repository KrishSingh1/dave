///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio.h
/// \brief LPJ-GUESS input/output module with input from instruction script
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GUESSIO_H
#define LPJ_GUESS_GUESSIO_H

#include "guess.h"

//Moved here temporarily from guessio_cru pending merge with trunk r2758
struct Coord {

	// Type for storing grid cell longitude, latitude and description text

	int id;
	double lon;
	double lat;
	xtring descrip;

};
//

void initio(const xtring& insfilename);
bool getgridcell(Gridcell& gridcell);
bool getclimate(Gridcell& gridcell);
void getlandcover(Gridcell& gridcell);
void getsowingdates(Gridcell& gridcell,Pftlist& pftlist);
void getharvestdates(Gridcell& gridcell,Pftlist& pftlist);
void outannual(Gridcell& gridcell);
void termio();
void printhelp();

#endif // LPJ_GUESS_GUESSIO_H
