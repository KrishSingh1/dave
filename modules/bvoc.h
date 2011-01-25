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

void bvoc(double,double,double,double,double,double,double,double,
	     double,double,double,double,double,double,double,Pft,
	     double&,double&,double&,double&,double&);
double dayT(double,double,double);
double leafT(double,double,double,double,double,double,double,double,double,
	     double);
void isoprmonot1(double,double,double,double,double,
		 double,double,double,double,
		 Pft,double,double,
		 double&,double&,double&);
void initbvoc(Pftlist&);
double vocseas(double&,double,double,double,Pft);
