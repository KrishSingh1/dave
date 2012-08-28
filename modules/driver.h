///////////////////////////////////////////////////////////////////////////////////////
/// \file driver.h
/// \brief Environmental driver calculation/transformation
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_DRIVER_H
#define LPJ_GUESS_DRIVER_H

#include "guess.h"

void setseed(long init);
double randfrac(long& seed);
void soilparameters(Soiltype& soiltype,int soilcode);
void interp_monthly_means(double mvals[12], double dvals[365]);
void interp_monthly_totals(double mvals[12], double dvals[365]);
void prdaily(double mval_prec[12],double dval_prec[365],double mval_wet[12]);
void dailyaccounting_gridcell(Gridcell& gridcell);
void dailyaccounting_stand(Stand& stand);
void dailyaccounting_patch(Patch& patch);
void respiration_temperature_response(double temp,double& gtemp);
void daylengthinsoleet(Climate& climate);
void soiltemp(Climate& climate,Soil& soil);

#endif // LPJ_GUESS_DRIVER_H
