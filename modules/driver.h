///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                Environmental driver calculation/transformation
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
//                        Includes weather generator and Dieter G:s latest updates
// Header file name:      driver.h
// Source code file name: driver.cpp
// Written by:            Ben Smith
// Version dated:         2002-11-22
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

void setseed(long init);
double randfrac();
void soilparameters(Soiltype& soiltype,int soilcode, double pawc); // guess2008 - emdi - new declaration
void initsoildrivers(Stand& stand);
void interp_climate(double mtemp[12],double mprec[12],double msun[12],
	double dtemp[365],double dprec[365],double dsun[365]);
void prdaily(double mval_prec[12],double dval_prec[365],double mval_wet[12]);
void dailyaccounting_stand(Stand& stand,Pftlist& pftlist);
void dailyaccounting_patch(Patch& patch);
void respiration_temperature_response(double temp,double& gtemp);
void daylengthinsoleet(Climate& climate);
void soiltemp(Climate& climate,Soil& soil);
