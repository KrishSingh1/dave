///////////////////////////////////////////////////////////////////////////////////////
/// \file gwgen.h
/// \brief Global Weather GENerator 
///
/// \author Lars Nieradzik
/// $Date: 2017-11-24 15:04:09 +0200 (Fri, 24 Nov 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

//#ifndef LPJ_GUESS_CRUINPUT_H
//#define LPJ_GUESS_CRUINPUT_H
//
#include "guess.h"
//#include "inputmodule.h"
//#include <vector>
//#include "gutil.h"
#ifndef GWGEN_H
#define GWGEN_H

#include <limits>

/// A weathergenerator for the use with e.g. BLAZE when wind etc is needed
/** This module generates daily weather data from monthly
 *  CRU-NCEP (1901-2015). 

 ENTER MORE HERE !!!!!
 */

//CLN// number of days in each month
const int ndaymonth[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };

// Classes (see below)

class GWGen /*: public Serializable */{

public:
	// MEMBER VARIABLES

	// new ones
	int month;

	// Derived datatype for the monthly weather generator input

	double mprec; // monthly total precipitation amount (mm)
	double mwetd; // number of days in month with precipitation
	double mwetf; // fraction of days in month with precipitation

	double mtmin; // minumum temperture (C)
	double mtmax; // maximum temperture (C)
	double mcldf; // cloud fraction (0=clear sky, 1=overcast) (fraction)
	double mwind; // wind speed (m/s)

	// CLN use them both for in and output
	bool pday[2]; //precipitation status: true if the day was a rain day
	//type(randomstate)      :: rndst   !state of the random number generator
	double resid[4];   //previous day's weather residuals

	//end type metvars_in

	//type metvars_out
	// Derived datatype for the daily weather generator output

	double dprec; // 24 hour total precipitation (mm)
	double dtmin; // 24 hour mean minimum temperature (degC)
	double dtmax; // 24 hour mean maximum temperature (degC)
	double dcldf; // 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
	double dwind; // wind speed (m s-1)
	double drhum; // relative humidity (CLN units!!!)
	//logical, dimension(2)  :: pday    !precipitation state
	//type(randomstate)      :: rndst   !state of the random number generator, 15 elements
	//real(sp), dimension(4) :: resid   !previous day's weather residuals
	double unorm[4];

	//end type metvars_out

	double tmn;
	double tmx;
	double wnd;
	double cld;

	//    type daymetvars
	// Derived datatype for monthly climate variables

	double dmtmax_mn; // maximum temperature monthly mean (degC)
	double dmtmin_mn; // minimum temperature mothly mean (degC)
	double dmcldf_mn; // mean cloud fraction (fraction)
	double dmwind_mn; // wind speed

	double dmtmax_sd; // standard deviation of corresponding variable above
	double dmtmin_sd; // ------- " ------
	double dmcldf_sd; // ------- " ------
	double dmwind_sd; // ------- " ------

	//end type daymetvars

	// the following parameters are computed by the cloud_params subroutine
	double cldf_w1, cldf_w2, cldf_w3, cldf_w4, cldf_d1, cldf_d2, cldf_d3, cldf_d4;
	double cldf_sd_w, cldf_sd_d;

	//void serialize(ArchiveStream& arch);  JN

	/// Constructor function: initialise cell member
};

void ran_seed(unsigned int sval, RnDst& state);
double ran_normal(RnDst& state);
int ranu(RnDst& state) ;
double ranur(RnDst& state);
int refill(RnDst& state) ;
double ran_gamma(RnDst& state,bool first, double shape, double scale);
double ran_gamma_gp(RnDst& state,bool first,double shape,double scale,double thresh,double shape_gp,double scale_gp) ;
double ran_gp(RnDst& state,double shape,double scale,double loc) ;
void calc_cloud_params(GWGen& gwgen);
void temp_sd(GWGen& gwgen);
void meansd(GWGen& gwgen) ;
double qchisq_appr(double p, double nu, double g, double tol);
double gamma_cdf_inv(double p, double alpha, double scale) ;
double gamma_cdf ( double x, double a, double b, double c);
double gamma_pdf ( double x, double a, double b, double c);
double gamma_inc ( double p, double x );
double gamma_log( double x );
double r8_gamma ( double x );
void normal_cdf_inv ( double cdf, double a, double b, double x );
void normal_01_cdf_inv (double p,double x);
double r8poly_value_horner ( int m, double *c, double x );
void normal_01_cdf ( double x, double cdf );
void matmul(double AA[4][4], double B[4], double CC[4]);
double roundto(double val,int precision);
void rmsmooth(int lm, double *m,int *dmonth,double bcond[2],int lr, double *r);
void init_weathergen(GWGen& gwgen, RnDst& rndst);
double cldf2rad(double input, double lat, int doy, bool cldf2rad);
void gwgen_get_daily_met(GWGen& gwgen, RnDst& rndst);
void gwgen_get_met(Gridcell& gridcell, double* in_mtemp, double* in_mprec, 
		   double* in_mwetm, double* in_msol, double* in_mdtr, 
		   double* in_mwind, double* in_rhum, double* out_temp,
		   double* out_dprec,double* out_dsol,double* out_ddtr,
		   double* out_dwind,double* out_rhum);


//class GWGen {
//	// MEMBER VARIABLES
//
//public:
//
//	// new ones
//	int month;
//
//        // Derived datatype for the monthly weather generator input
//
//        double mprec ; // monthly total precipitation amount (mm)
//        double mwetd ; // number of days in month with precipitation
//        double mwetf ; // fraction of days in month with precipitation
//	
//        double mtmin ; // minumum temperture (C)
//        double mtmax ; // maximum temperture (C)
//        double mcldf ; // cloud fraction (0=clear sky, 1=overcast) (fraction)
//        double mwind ; // wind speed (m/s)
//
//	// CLN use them both for in and output
//        bool pday[2]; //precipitation status: true if the day was a rain day
//        //type(randomstate)      :: rndst   !state of the random number generator
//        double resid[4];   //previous day's weather residuals
//
//	//end type metvars_in
//
//	//type metvars_out
//        // Derived datatype for the daily weather generator output
//
//        double dprec ; // 24 hour total precipitation (mm)
//        double dtmin ; // 24 hour mean minimum temperature (degC)
//        double dtmax ; // 24 hour mean maximum temperature (degC)
//        double dcldf ; // 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
//        double dwind ; // wind speed (m s-1)
//
//        //logical, dimension(2)  :: pday    !precipitation state
//        //type(randomstate)      :: rndst   !state of the random number generator, 15 elements
//        //real(sp), dimension(4) :: resid   !previous day's weather residuals
//        double unorm[4];
//
//	//end type metvars_out
//
//	double tmn;
//	double tmx;
//	double wnd;
//	double cld;
//	
//	//    type daymetvars
//        // Derived datatype for monthly climate variables
//
//        double dmtmax_mn ; // maximum temperature monthly mean (degC)
//        double dmtmin_mn ; // minimum temperature mothly mean (degC)
//        double dmcldf_mn ; // mean cloud fraction (fraction)
//        double dmwind_mn ; // wind speed
//
//        double dmtmax_sd ; // standard deviation of corresponding variable above
//        double dmtmin_sd ; // ------- " ------
//        double dmcldf_sd ; // ------- " ------
//        double dmwind_sd ; // ------- " ------
//
//	//end type daymetvars
//	
//	// the following parameters are computed by the cloud_params subroutine
//	double cldf_w1, cldf_w2, cldf_w3, cldf_w4, cldf_d1, cldf_d2, cldf_d3, cldf_d4;
//	double cldf_sd_w,cldf_sd_d;
//	// seed for rondom generation
//	int seed;
//};
//
//class RnDst {
//	// MEMBER VARIABLES
//
//public:
//	
//	int q[10];
//        int carry ;
//        int xcng  ;
//        unsigned int xs    ; //!default seed
//        int indx  ;
//	bool have  ;
//	double gamma_vals[2];
//
//};
	
#endif // GWGEN_H
