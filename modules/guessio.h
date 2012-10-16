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

void initio(int argc,char* argv[],Pftlist& pftlist);
bool getgridcell(Gridcell& gridcell);
bool getclimate(Gridcell& gridcell);
void getlandcover(Gridcell& gridcell,Pftlist& pftlist);
void outannual(Gridcell& gridcell,Pftlist& pftlist);
void termio();

#endif // LPJ_GUESS_GUESSIO_H
