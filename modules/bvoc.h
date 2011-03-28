///////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Calculation of VOC production and emission by 
//                        vegetation 
//                        *****************************************************
// Header file name:      voccalc.h
// Source code file name: voccalc.cpp
// Written by:            Guy Schurgers (using Almut's previous attempts) 
// Version dated:         August 2006
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions 
// defined in the module that are to be accessible to the calling framework or 
// to other modules.

#ifndef LPJ_GUESS_BVOC_H
#define LPJ_GUESS_BVOC_H

#include "guess.h"

void bvoc(double daylength, double temp, double tempamp,
	  const PhotosynthesisResult& photosynthesis, 
	  double co2, double lambda, double eet, double agdd5, int nday,
	  double rs_day, double lai, const Pft& pft,
	  double& iso, double& mon, double& dmonstor, double& dleaftemp,
	  double& fvocseas);
void initbvoc(Pftlist& pftlist);

#endif // LPJ_GUESS_BVOC_H
