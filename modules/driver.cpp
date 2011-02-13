///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Environmental driver calculation/transformation
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
//                        cohort/individual mode - see canexch.cpp)
//                        Includes weather generator and Dieter G:s latest updates
// Header file name:      driver.h
// Source code file name: driver.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
// Updated:               2010-11-22


// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "driver.h"


static long seed=12345678; // seed for random number generator (see randfrac)


// guess2008
extern int nyear_spinup; 
	// allows access to the value declared guessio_cru.cpp
 

///////////////////////////////////////////////////////////////////////////////////////
// RANDFRAC
// Internal function for generating random numbers

void setseed(long init) {

	seed=init;
}

double randfrac() {

	// DESCRIPTION
	// Returns a random floating-point number in the range 0-1.
	// Uses and updates the global variable 'seed' which may be initialised to any
	// positive integral value (the same initial value will result in the same sequence
	// of returned values on subsequent calls to randfloat every time the program is
	// run)

	// Reference: Park & Miller 1988 CACM 31: 1192

	const long modulus=2147483647;
	const double fmodulus=modulus;
	const long multiplier=16807;
	const long q=127773;
	const long r=2836;

	seed=multiplier*(seed%q)-r*seed/q;
	if (!seed) seed++; // increment seed to 1 in unlikely event of 0 value
	else if (seed<0) seed+=modulus;
	return (double)seed/fmodulus;
}

///////////////////////////////////////////////////////////////////////////////////////
// SOILPARAMETERS
// May be called from input/output module to initialise stand Soiltype objects when
// soil data supplied as LPJ soil code rather than soil physical parameter values


void soilparameters(Soiltype& soiltype,int soilcode) {

	// DESCRIPTION
	// Derivation of soil physical parameters given LPJ soil code

	// INPUT AND OUTPUT PARAMETER
	// soil = patch soil

	const double PERC_EXP=2.0;
		// exponent in percolation equation [k2; LPJF]
		// (Eqn 31, Haxeltine & Prentice 1996)
		// Changed from 4 to 2 (Sitch, Thonicke, pers comm 26/11/01)

	double data[9][5]= {

		//    0  empirical parameter in percolation equation (k1) (mm/day)
		//    1  volumetric water holding capacity at field capacity minus vol water
		//       holding capacity at wilting point (Hmax), as fraction of soil layer
		//       depth
		//    2  thermal diffusivity (mm2/s) at wilting point (0% WHC)
		//    3  thermal diffusivity (mm2/s) at 15% WHC
		//    4  thermal diffusivity at field capacity (100% WHC)
		//       Thermal diffusivities follow van Duin (1963),
		//       Jury et al (1991), Fig 5.11.

		//    0      1      2      3      4   soilcode
		//  ------------------------------------------

		{   5.0, 0.110,   0.2, 0.800,   0.4 },   // 1
		{   4.0, 0.150,   0.2, 0.650,   0.4 },   // 2
		{   3.0, 0.120,   0.2, 0.500,   0.4 },   // 3
		{   4.5, 0.130,   0.2, 0.725,   0.4 },   // 4
		{   4.0, 0.115,   0.2, 0.650,   0.4 },   // 5
		{   3.5, 0.135,   0.2, 0.575,   0.4 },   // 6
		{   4.0, 0.127,   0.2, 0.650,   0.4 },   // 7
		{   9.0, 0.300,   0.1, 0.100,   0.1 },   // 8
		{   0.2, 0.100,   0.2, 0.500,   0.4 }    // 9
	};

	if (soilcode<1 || soilcode>9)
		fail("soilparameters: invalid LPJ soil code (%d)",soilcode);

	
	soiltype.perc_base=data[soilcode-1][0];
	soiltype.perc_exp=PERC_EXP;
	soiltype.awc[0]=SOILDEPTH_UPPER*data[soilcode-1][1];
	soiltype.awc[1]=SOILDEPTH_LOWER*data[soilcode-1][1];
	soiltype.thermdiff_0=data[soilcode-1][2];
	soiltype.thermdiff_15=data[soilcode-1][3];
	soiltype.thermdiff_100=data[soilcode-1][4];

	// guess2008 - override the default SOM years with 70-80% of the spin-up period
	soiltype.updateSolveSOMvalues(nyear_spinup);
}


///////////////////////////////////////////////////////////////////////////////////////
// CLIMATE INTERPOLATION FROM MONTHLY TO QUASI-DAILY VALUES
// May be called from input/output module to generate daily climate values when raw
// data are on monthly basis

void interp_climate(double mtemp[12],double mprec[12],double msun[12],
	double dtemp[365],double dprec[365],double dsun[365]) {

	// DESCRIPTION
	// Interpolates monthly climate data in arrays mtemp, mprec and msun to daily
	// values in arrays dtemp, dprec and dsun

	// INPUT PARAMETERS
	// mtemp = mean monthly temperatures for this year (deg C)
	// mprec = monthly precipitation totals for this year (mm)
	// msun  = mean monthly percentage sunshine values for this year

	// OUTPUT PARAMETERS
	// dtemp = mean daily temperatures for this year (deg C)
	// dprec = mean daily precipitation values for this year (mm)
	// dsun  = mean daily percentage sunshine values for this year

	Date date; // Date object used for interpolation (local to this function)
	double nday,dayct;
	int thismonth,lastmonth,m;
	double mprec_daily[12];

	// Convert monthly precipitation from monthly totals to mean daily values

	for (m=0;m<12;m++)
		mprec_daily[m]=mprec[m]/(double)date.ndaymonth[m];

	date.init(1);

	nday=(double)(date.middaymonth[0]-(date.middaymonth[11]-365));
	thismonth=0;
	lastmonth=11;
	dayct=(double)(366-date.middaymonth[11]);

	// Perform interpolations

	while (date.year==0) {
		if (date.day==date.middaymonth[date.month]) {
			if (date.month==11) // December
				nday=(double)(date.middaymonth[0]+365-date.middaymonth[11]);
			else
				nday=(double)(date.middaymonth[date.nextmonth()]-
					date.middaymonth[date.month]);
			thismonth=date.nextmonth();
			lastmonth=date.month;
			dayct=0.0;
		}
		dtemp[date.day]=(mtemp[thismonth]-mtemp[lastmonth])/nday*dayct+
			mtemp[lastmonth];
		dsun[date.day]=(msun[thismonth]-msun[lastmonth])/nday*dayct+msun[lastmonth];
		dprec[date.day]=(mprec_daily[thismonth]-mprec_daily[lastmonth])/nday*dayct+
			mprec_daily[lastmonth];
		date.next();
		dayct++;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
//  PRDAILY
//  Distribution of monthly precipitation totals to quasi-daily values
//  (From Dieter Gerten 021121)

void prdaily(double mval_prec[12],double dval_prec[365],double mval_wet[12]) {

	// mval_prec = total rainfall (mm) for month
	// dval_prec = actual rainfall (mm) for each day of year
	// mval_wet  = expected number of rain days for month

	const double c1=1.0; // normalising coefficient for exponential distribution
	const double c2=1.2; // power for exponential distribution

	int m,d,dy,dyy,dy_hold;
	int daysum;
	double prob_rain; // daily probability of rain for this month
	double mprec; // average rainfall per rain day for this month
	double mprec_sum; // cumulative sum of rainfall for this month
		// (= mprecip in Dieter's code)
	double prob;

	dy=0;
	daysum=0;

	for (m=0;m<12;m++) {

		if (negligible(mval_prec[m])) {

			// Special case if no rainfall expected for month

			for (d=0;d<date.ndaymonth[m];d++) {
				dval_prec[dy]=0.0;
				dy++;
			}
		}
		else {

			mprec_sum=0.0;

			if (negligible(mval_wet[m])) mval_wet[m]=1.0;
				// force at least one rain day per month

			prob_rain=mval_wet[m]/(double)date.ndaymonth[m];

			mprec=mval_prec[m]/mval_wet[m];
			
			dy_hold=dy;

			while (negligible(mprec_sum)) {

				dy=dy_hold;

				for (d=0;d<date.ndaymonth[m];d++) {
				
					// Transitional probabilities (Geng et al 1986)

					if (dy==0) { // first day of year only
						prob=0.75*prob_rain;
					}
					else {
						if (dval_prec[dy-1]<0.1)
							prob=0.75*prob_rain;
						else
							prob=0.25+(0.75*prob_rain);
					}

					// Determine wet days randomly and use Krysanova/Cramer estimates of
					// parameter values (c1,c2) for an exponential distribution

					if (randfrac()>prob)
						dval_prec[dy]=0.0;
					else {
						double x=randfrac();
						dval_prec[dy]=pow(-log(x),c2)*mprec*c1;
						if (dval_prec[dy]<0.1) dval_prec[dy]=0.0;
					}

					mprec_sum+=dval_prec[dy];
					dy++;
				}

				// Normalise generated precipitation by prescribed monthly totals

				if (!negligible(mprec_sum)) {
					for (d=0;d<date.ndaymonth[m];d++) {
						dyy=daysum+d;
						dval_prec[dyy]*=mval_prec[m]/mprec_sum;
						if (dval_prec[dyy]<0.1) dval_prec[dyy]=0.0;
					}
				}
			}
		}

		daysum+=date.ndaymonth[m];
	}
}

///////////////////////////////////////////////////////////////////////////////////////
//  SOIL TEMPERATURES
//  Call each simulation day following update of daily air temperature prior to canopy
//  exchange and SOM dynamics

inline void regress(double* x,double* y,int n,double& a,double& b) {

	// Performs a linear regression of array y on array x (n values)
	// returning parameters a and b in the fitted model: y=a+bx
	// (Used by function soiltemp)
	// Source: Press et al 1986, Sect 14.2

	int i;
	double sx,sy,sxx,sxy,delta;

	sx=0.0;
	sy=0.0;
	sxx=0.0;
	sxy=0.0;
	for (i=0;i<n;i++) {
		sx+=x[i];
		sy+=y[i];
		sxx+=x[i]*x[i];
		sxy+=x[i]*y[i];
	}
	delta=(double)n*sxx-sx*sx;
	a=(sxx*sy-sx*sxy)/delta;
	b=((double)n*sxy-sx*sy)/delta;
}


void soiltemp(Climate& climate,Soil& soil) {

	// DESCRIPTION
	// Calculation of soil temperature at 0.25 m depth (middle of upper soil layer).
	// Soil temperatures are assumed to follow surface temperatures according to an
	// annual sinusoidal cycle with damped oscillation about a common mean, and a
	// temporal lag.
	
	// For a sinusoidal cycle, soil temperature at depth z and time t from beginning
	// of cycle given by (Carslaw & Jaeger 1959; Eqn 52; Jury et al 1991):
	//
	//   (1) t(z,t) = t_av + a*exp(-z/d)*sin(omega*t - z/d)
	//
	//   where
	//     t_av      = average (base) air/soil temperature
	//     a         = amplitude of air temp fluctuation
	//     exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z,
	//                 relative to surface temp fluctuation
	//     z/d       = oscillation lag in angular units at soil depth z
	//     z         = soil depth
	//     d         = sqrt(2*k/omega), damping depth
	//     k         = soil thermal diffusivity
	//     omega     = 2*PI/tau, angular frequency of oscillation (radians)
	//     tau       = oscillation period (365 days)
	//
	// Here we assume a sinusoidal cycle, but estimate soil temperatures based on
	// a lag (z/d, converted from angular units to days) relative to air temperature
	// and damping ( exp(-z/d) ) of soil temperature amplitude relative to air
	// temperature amplitude. A linear model for change in air temperature with time
	// during the last 31 days is used to estimate 'lagged' air temperature. Soil
	// temperature today is thus given by:
	//
	//   (2) temp_soil = atemp_mean + exp( -alag ) * ( temp_lag - temp_mean )
	//
	//   where
	//     atemp_mean = mean of monthly mean temperatures for the last year (deg C)
	//     alag       = oscillation lag in angular units at depth 0.25 m
	//                  (corresponds to z/d in Eqn 1)
	//     temp_lag   = air temperature 'lag' days ago (estimated from linear model)
	//                  where 'lag' = 'alag' converted from angular units to days
	//
	// Soil thermal diffusivity (k) is sensitive to soil water content and is estimated
	// monthly based on mean daily soil water content for the past month, interpolating
	// between estimates for 0, 15% and 100% AWHC (Van Duin 1963; Jury et al 1991,
	// Fig 5.11).

	const double DIFFUS_CONV=0.0864;
		// conversion factor for soil thermal diffusivity from mm2/s to m2/day
	const double HALF_OMEGA=8.607E-3; // corresponds to omega/2 = pi/365 (Eqn 1)
	const double DEPTH=SOILDEPTH_UPPER*0.0005;
		// soil depth at which to estimate temperature (m)
	const double LAG_CONV=58.09;
		// conversion factor for oscillation lag from angular units to days (=365/(2*PI))

	double a,b; // regression parameters
	double k; // soil thermal diffusivity (m2/day)
	double temp_lag; // air temperature 'lag' days ago (see above; deg C)
	double day[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
		16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};

	if ((date.year==0 || date.year==soil.patch.stand.first_year) && date.month==0 && !date.islastday) {

		// First month of simulation, use air temperature for soil temperature

		soil.temp=climate.temp;
	}
	else {

		if (date.islastday) {

			// Linearly interpolate soil thermal diffusivity given mean
			// soil water content

			if (soil.mwcontupper<0.15)
				k=((soil.soiltype.thermdiff_15-soil.soiltype.thermdiff_0)/0.15*
					soil.mwcontupper+soil.soiltype.thermdiff_0)*DIFFUS_CONV;
			else
				k=((soil.soiltype.thermdiff_100-soil.soiltype.thermdiff_15)/0.85*
					(soil.mwcontupper-0.15)+soil.soiltype.thermdiff_15)*DIFFUS_CONV;

			// Calculate parameters alag and exp(-alag) from Eqn 2

			soil.alag=DEPTH/sqrt(k/HALF_OMEGA); // from Eqn 1
			soil.exp_alag=exp(-soil.alag);

		}

		// Every day, calculate linear model for trend in daily air
		// temperatures for the last 31 days: temp_day = a + b * day

		regress(day,climate.dtemp_31,31,a,b);

		// Calculate soil temperature

		temp_lag=a+b*(30.0-soil.alag*LAG_CONV);
		soil.temp=climate.atemp_mean+soil.exp_alag*(temp_lag-climate.atemp_mean);
			// Eqn 2
	}
}


inline double mean(double* array,int nitem) {

	// Returns arithmetic mean of 'nitem' values in 'array'
	// (Used by function soiltemp)

	double sum=0.0;
	int i;

	for (i=0;i<nitem;sum+=array[i++]);
	return sum/(double)nitem;
}

/// Called each simulation day before any other driver or process functions
void dailyaccounting_gridcell(Gridcell& gridcell,Pftlist& pftlist) {

	// DESCRIPTION
	// Updates daily climate parameters including growing degree day sums and
	// exponential temperature response term (gtemp, see below). Maintains monthly
	// and longer term records of variation in climate variables. PFT-specific
	// degree-day sums in excess of damaging temperatures are also calculated here.

	const double W11DIV12=11.0/12.0;
	const double W1DIV12=1.0/12.0;
	int d,y,startyear;

	// guess2008 - changed this from an int to a double
	double mtemp_last;

	Climate& climate=gridcell.climate;

	// On first day of year ...

	if (date.day==0) {
		// ... reset annual GDD5 counter
		climate.agdd5=0.0;

		if (date.year==0) {
			// First day of simulation - initialise running annual mean temperature and daily temperatures for the last month
			for (d=0;d<31;d++)
				climate.dtemp_31[d]=climate.temp;
			climate.atemp_mean=climate.temp;
		}
	}
	else if (climate.lat>=0.0 && date.day==COLDEST_DAY_NHEMISPHERE ||
		climate.lat<0.0 && date.day==COLDEST_DAY_SHEMISPHERE) {
		// In midwinter, reset GDD counter for summergreen phenology
		climate.gdd5=0.0;
		climate.ifsensechill=false; // guess2008 - CHILLDAYS
	}

	// Update GDD counters and chill day count
	climate.gdd5+=max(0.0,climate.temp-5.0);
	climate.agdd5+=max(0.0,climate.temp-5.0);
	if (climate.temp<5.0 && climate.chilldays<=365)
		climate.chilldays++;

///	if (run_landuse && run_crop)
///		dailyaccounting_gridcell_crop(gridcell,pftlist);

	// Save yesterday's mean temperature for the last month
	mtemp_last=climate.mtemp;

	// Update daily temperatures, and mean overall temperature, for last 31 days
	climate.mtemp=climate.temp;
	for (d=0;d<30;d++) {
		climate.dtemp_31[d]=climate.dtemp_31[d+1];
		climate.mtemp+=climate.dtemp_31[d];
	}
	climate.dtemp_31[30]=climate.temp;
	climate.mtemp/=31.0;

	// Reset GDD and chill day counter if mean monthly temperature falls below base
	// temperature
	if (mtemp_last>=5.0 && climate.mtemp<5.0 && climate.ifsensechill) { // guess2008 - CHILLDAYS
		climate.gdd5=0.0;
		climate.chilldays=0;
	}

	// On last day of month ...

	if (date.islastday) {
		// Update mean temperature for the last 12 months
		// atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)
		climate.atemp_mean=climate.atemp_mean*W11DIV12+climate.mtemp*W1DIV12;
		
		// Record minimum and maximum monthly temperatures
		if (date.month==0) {
			climate.mtemp_min=climate.mtemp;
			climate.mtemp_max=climate.mtemp;
		}
		else {
			if (climate.mtemp<climate.mtemp_min)
				climate.mtemp_min=climate.mtemp;
			if (climate.mtemp>climate.mtemp_max)
				climate.mtemp_max=climate.mtemp;
		}

		// On 31 December update records of minimum monthly temperatures for the last
		// 20 years and find mean of minimum monthly temperatures for the last 20 years
		if (date.islastmonth) {
			startyear=20-(int)min(19,date.year);
			climate.mtemp_min20=climate.mtemp_min;
			climate.mtemp_max20=climate.mtemp_max;

			for (y=startyear;y<20;y++) {
				climate.mtemp_min_20[y-1]=climate.mtemp_min_20[y];
				climate.mtemp_min20+=climate.mtemp_min_20[y];
				climate.mtemp_max_20[y-1]=climate.mtemp_max_20[y];
				climate.mtemp_max20+=climate.mtemp_max_20[y];
			}

			climate.mtemp_min20/=(double)(21-startyear);
			climate.mtemp_max20/=(double)(21-startyear);
			climate.mtemp_min_20[19]=climate.mtemp_min;
			climate.mtemp_max_20[19]=climate.mtemp_max;
		}
	}
}

void dailyaccounting_stand(Stand& stand,Pftlist& pftlist)
{		
	// Loop through PFTs
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft=pftlist.getobj();
		// For this PFT ...

		// [BEGIN CEFAST0207]
		// Flag used by evapotranspiration_fast in canopy exchange module ...
		stand.pft[pft.id].have_phot=false;

		// [END CEFAST0207]

		// ... on to next PFT
		pftlist.nextobj();
	}
}

void dailyaccounting_patch_lc(Patch& patch, Pftlist& pftlist) {
	if(date.day==0) {
		Fluxes& fluxes=patch.fluxes;

		if(!patch.stand.gridcell.LC_updated) {	// NB. landcover_dynamics() is called before this function !
			fluxes.acflux_harvest=0.0;
		}

		if(ifslowharvestpool) {
			pftlist.firstobj();
			while(pftlist.isobj) {
				Pft& pft=pftlist.getobj();
				Patchpft& patchpft=patch.pft[pft.id];

				fluxes.acflux_harvest+=patchpft.harvested_products_slow*pft.turnover_harv_prod;
				patchpft.harvested_products_slow=patchpft.harvested_products_slow*(1-pft.turnover_harv_prod);

				pftlist.nextobj();
			}
		}
	}
}

void dailyaccounting_patch(Patch& patch, Pftlist& pftlist) {
	// DESCRIPTION
	// Updates daily soil parameters including exponential temperature response terms
	// (gtemp, see below). Maintains monthly and longer term records of variation in
	// soil variables. Initialises flux sums at start of simulation year.

	// INPUT AND OUTPUT PARAMETER
	// soil   = patch soil
	// fluxes = current and accumulated C fluxes for patch

	//int p;
	Soil& soil=patch.soil;
	Fluxes& fluxes=patch.fluxes;

	if (date.day==0) {

		// Reset fluxes
		fluxes.acflux_soil=0.0;
		fluxes.acflux_veg=0.0;
		fluxes.acflux_est=0.0;
		fluxes.acflux_fire=0.0;

		patch.aaet=0.0;
		patch.aevap=0.0;
		patch.arunoff=0.0;
		patch.aintercep=0.0;
		patch.apet=0.0;
	}

	if (date.dayofmonth==0) {

		fluxes.mcflux_veg[date.month]=0.0;

		patch.maet[date.month]=0.0;
		patch.mevap[date.month]=0.0;
		patch.mrunoff[date.month]=0.0;
		patch.mintercep[date.month]=0.0;
		patch.mpet[date.month]=0.0;

		// guess2008 - reset month C budget arrays each month
		fluxes.mcflux_gpp[date.month] = 0.0;
		fluxes.mcflux_ra[date.month] = 0.0;

	}

	fluxes.dcflux_veg=0.0;

	if(run_landcover)
		dailyaccounting_patch_lc(patch, pftlist);
	
	// Store daily soil water in upper layer
	soil.dwcontupper[date.day]=soil.wcont[0];

	// Store daily soil water in lower layer - guess2008
	soil.dwcontlower[date.day]=soil.wcont[1];

	// On last day of month, calculate mean content of upper soil layer

	if (date.islastday) {

		soil.mwcontupper=mean(soil.dwcontupper+date.day-date.ndaymonth[date.month]+1,
			date.ndaymonth[date.month]);
		
		// guess2008 - record water in lower layer too, and then update mwcont  
		soil.mwcontlower=mean(soil.dwcontlower+date.day-date.ndaymonth[date.month]+1,
			date.ndaymonth[date.month]);
		
		soil.mwcont[date.month][0] = soil.mwcontupper;
		soil.mwcont[date.month][1] = soil.mwcontlower;

	}

	// Calculate soil temperatures
	soiltemp(patch.stand.gridcell.climate,soil);

	// On last day of month, calculate mean soil temperature for last month

	soil.dtemp[date.dayofmonth]=soil.temp;

	if (date.islastday)
		soil.mtemp=mean(soil.dtemp,date.ndaymonth[date.month]);
}


///////////////////////////////////////////////////////////////////////////////////////
// RESPIRATION TEMPERATURE RESPONSE
// Called by canopy exchange and soil organic matter dynamics module to calculate
// response of respiration to temperature

void respiration_temperature_response(double temp,double& gtemp) {

	// DESCRIPTION
	// Calculates g(T), response of respiration rate to temperature (T), based on
	// empirical relationship for temperature response of soil temperature across
	// ecosystems, incorporating damping of Q10 response due to temperature
	// acclimation (Eqn 11, Lloyd & Taylor 1994)
	//
	//   r    = r10 * g(t)
	//   g(T) = EXP [308.56 * (1 / 56.02 - 1 / (T - 227.13))] (T in Kelvin)

	// INPUT PARAMETER
	// temp = air or soil temperature (deg C)

	// OUTPUT PARAMETER
	// gtemp = respiration temperature response

	if (temp>=-40.0)
		gtemp=exp(308.56*(1.0/56.02-1.0/(temp+46.02))); // NB: temperature in deg C
	else
		gtemp=0.0;
}


///////////////////////////////////////////////////////////////////////////////////////
// DAYLENGTH, INSOLATION AND POTENTIAL EVAPOTRANSPIRATION
// Called by framework each simulation day following update of daily air temperature 
// and before canopy exchange processes

void daylengthinsoleet(Climate& climate) {

	// DESCRIPTION
	// Calculation of daylength, insolation and equilibrium evapotranspiration
	// for each day, given mean daily temperature, insolation (as percentage
	// of full sunshine or mean daily instantaneous downward shortwave
	// radiation flux, W/m2), latitude and day of year

	// INPUT AND OUTPUT PARAMETER
	// climate = stand climate

	const double QOO=1360.0;
	const double PI=3.1415927;
	const double BETA=0.17;
	const double A=107.0;
	const double B=0.2;
	const double C=0.25;
	const double D=0.5;
	const double K=13750.98708;
	const double DEGTORAD=0.01745329;
	const double FRADPAR=0.5;
		// fraction of net incident shortwave radiation that is photosynthetically
		// active (PAR)

	double delta; // solar declination angle (radians)
	double rs_day; // daily net downward shortwave radiation sum (J/m2/day)
	double rl; // instantaneous net upward longwave radiation flux (W/m2)
	double eet_day; // equilibrium evapotranspiration sum (mm/day)
	double w,gamma,lambda,s,uu,vv,hn;

	//	CALCULATION OF NET DOWNWARD SHORT-WAVE RADIATION FLUX
	//	Refs: Prentice et al 1993, Monteith & Unsworth 1990,
	//	      Henderson-Sellers & Robinson 1986

	//	 (1) rs = (c + d*ni) * (1 - beta) * Qo * cos Z * k
	//	       (Eqn 7, Prentice et al 1993)
	//	 (2) Qo = Qoo * ( 1 + 2*0.01675 * cos ( 2*pi*(i+0.5)/365) )
	//	       (Eqn 8, Prentice et al 1993; angle in radians)
	//	 (3) cos Z = sin(lat) * sin(delta) + cos(lat) * cos(delta) * cos h
	//	       (Eqn 9, Prentice et al 1993)
	//	 (4) delta = -23.4 * pi / 180 * cos ( 2*pi*(i+10.5)/365 )
	//	       (Eqn 10, Prentice et al 1993, angle in radians)
	//	 (5) h = 2 * pi * t / 24 = pi * t / 12

	//	     where rs    = instantaneous net downward shortwave radiation
	//	                   flux, including correction for terrestrial shortwave albedo
	//	                   (W/m2 = J/m2/s)
	//	           c, d  = empirical constants (c+d = clear sky
	//	                   transmissivity)
	//	           ni    = proportion of bright sunshine
	//	           beta  = average 'global' value for shortwave albedo
	//	                   (not associated with any particular vegetation)
	//	           i     = julian day, (0-364, 0=1 Jan)
	//	           Qoo   = solar constant, 1360 W/m2
	//	           Z     = solar zenith angle (angular distance between the
	//	                   sun's rays and the local vertical)
	//	           k     = conversion factor from solar angular units to
	//	                   seconds, 12 / pi * 3600
	//	           lat   = latitude (+=N, -=S, in radians)
	//	           delta = solar declination (angle between the orbital
	//	                   plane and the Earth's equatorial plane) varying
	//	                   between +23.4 degrees in northern hemisphere
	//	                   midsummer and -23.4 degrees in N hemisphere
	//	                   midwinter
	//	           h     = hour angle, the fraction of 2*pi (radians) which
	//	                   the earth has turned since the local solar noon
	//	           t     = local time in hours from solar noon

	//	From (1) and (3), shortwave radiation flux at any hour during the
	//	day, any day of the year and any latitude given by
	//	 (6) rs = (c + d*ni) * (1 - beta) * Qo * ( sin(lat) * sin(delta) +
	//	          cos(lat) * cos(delta) * cos h ) * k
	//	Solar zenith angle equal to -pi/2 (radians) at sunrise and pi/2 at
	//	sunset.  For Z=pi/2 or Z=-pi/2,
	//	 (7) cos Z = 0
	//	From (3) and (7),
	//	 (8)  cos hh = - sin(lat) * sin(delta) / ( cos(lat) * cos(delta) )
	//	      where hh = half-day length in angular units
	//	Define
	//	 (9) u = sin(lat) * sin(delta)
	//	(10) v = cos(lat) * cos(delta)
	//	Thus 
	//	(11) hh = acos (-u/v)
	//	To obtain the daily net downward short-wave radiation sum, integrate
	//	equation (6) from -hh to hh with respect to h,
	//	(12) rs_day = 2 * (c + d*ni) * (1 - beta) * Qo *
	//	              ( u*hh + v*sin(hh) )
	//	Define
	//	(13) w = (c + d*ni) * (1 - beta) * Qo
	//	From (12) & (13), and converting from angular units to seconds
	//	(14) rs_day = 2 * w * ( u*hh + v*sin(hh) ) * k

	if (!climate.doneday[date.day]) {

		// Calculate values of saved parameters for this day

		climate.qo[date.day]=QOO*(1.0+2.0*0.01675
			*cos(2.0*PI*((double)date.day+0.5)/365.0)); // Eqn 2
		delta=-23.4*DEGTORAD*cos(2.0*PI*((double)date.day+10.5)/365.0); // Eqn 4
		climate.u[date.day]=climate.sinelat*sin(delta); // Eqn 9
		climate.v[date.day]=climate.cosinelat*cos(delta); // Eqn 10

		if (climate.u[date.day]>=climate.v[date.day])
			climate.hh[date.day]=PI; // polar day
		else if (climate.u[date.day]<=-climate.v[date.day])
			climate.hh[date.day]=0.0; // polar night
		else climate.hh[date.day]=
			acos(-climate.u[date.day]/climate.v[date.day]); // Eqn 11

		climate.sinehh[date.day]=sin(climate.hh[date.day]);

		// Calculate daylength in hours from hh

		climate.daylength_save[date.day]=24.0*climate.hh[date.day]/PI;

		climate.doneday[date.day]=true;
	}

	if (climate.instype==SUNSHINE) { // insolation provided as percentage sunshine
		
		w=(C+D*climate.insol/100.0)*(1.0-BETA)*climate.qo[date.day]; // Eqn 13
		rs_day=2.0*w*(climate.u[date.day]*climate.hh[date.day]
			+climate.v[date.day]*climate.sinehh[date.day])*K; // Eqn 14

	}
	else { // insolation provided as instantaneous downward shortwave radiation flux 
		   // Replace climate.daylength_save[date.day] with 24 when using radiation data 
		   // representing the average time period radiation.  

		if (climate.instype==NETSWRAD) // net radiation known
			rs_day=climate.insol*climate.daylength_save[date.day]*3600.0;
		else // include correction for albedo
			rs_day=climate.insol*(1.0-BETA)*climate.daylength_save[date.day]*3600.0;

		// guess2008 - special case for polar night
		if (climate.sinehh[date.day]<0.001) {
			w=0.0 ; // polar night
		}
		else {
			w=rs_day/2.0/(climate.u[date.day]*climate.hh[date.day]
				+climate.v[date.day]*climate.sinehh[date.day])/K; // from Eqn 14
		}
	}

	//	CALCULATION OF DAILY EQUILIBRIUM EVAPOTRANSPIRATION
	//	(EET, or evaporative demand)
	//	Refs: Jarvis & McNaughton 1986, Prentice et al 1993

	//	(15) eet = ( s / (s + gamma) ) * rn / lambda
	//	       (Eqn 5, Prentice et al 1993)
	//	(16) s = 2.503E+6 * exp ( 17.269 * temp / (237.3 + temp) ) /
	//	         (237.3 + temp)**2
	//	       (Eqn 6, Prentice et al 1993)
	//	(17) rn = rs - rl
	//	(18) rl = ( b + (1-b) * ni ) * ( a - temp )
	//	       (Eqn 11, Prentice et al 1993)

	//	     where eet    = instantaneous evaporative demand (mm/s)
	//	           gamma  = psychrometer constant, c. 65 Pa/K
	//	           lambda = latent heat of vapourisation of water,
	//	                    c. 2.5E+6 J/kg
	//	           temp   = temperature (deg C)
	//	           rl     = net upward longwave radiation flux
	//	                    ('terrestrial radiation') (W/m2)
	//	           rn     = net downward radiation flux (W/m2)
	//	           a, b   = empirical constants

	//	Note: gamma and lambda are weakly dependent on temperature. Simple
	//	      linear functions are used to obtain approximate values for a
	//	      given temperature

	//	From (13) & (18),
	//	(19) rl = ( b + (1-b) * ( w / Qo / (1 - beta) - c ) / d ) * ( a - temp )

	//	Define
	//	(20) uu = w * u - rl
	//	(21) vv = w * v

	//	Daily EET sum is instantaneous EET integrated over the period
	//	during which rn >= 0.  Limits for the integration (half-period
	//	hn) are obtained by solving for

	//	(22) rn=0
	//	From (17) & (22),
	//	(23) rs - rl = 0
	//	From (6), (20), (21) and (23),
	//	(24) uu + vv * cos hn = 0
	//	From (24),
	//	(25) hn = acos ( -uu/vv )

	//	Integration of (15) w.r.t. h in the range -hn to hn leads to the
	//	following formula for total daily EET (mm)

	//	(26) eet_day = 2 * ( s / (s + gamma) / lambda ) *
	//	               ( uu*hn + vv*sin(hn) ) * k

	rl=(B+(1.0-B)*(w/climate.qo[date.day]/(1.0-BETA)-C)/D)*(A-climate.temp); // Eqn 19

	//	Calculate gamma and lambda

	gamma=65.05+climate.temp*0.064;
	lambda=2.495E6-climate.temp*2380.0;

	s=2.503E6*exp(17.269*climate.temp/(237.3+climate.temp))/
		(237.3+climate.temp)/(237.3+climate.temp); // Eqn 16

	uu=w*climate.u[date.day]-rl; // Eqn 20
	vv=w*climate.v[date.day]; // Eqn 21

	//	Calculate half-period with positive net radiation, hn
	//	In Eqn (25), hn defined for uu in range -vv to vv
	//	For uu >= vv, hn = pi (12 hours, i.e. polar day)
	//	For uu <= -vv, hn = 0 (i.e. polar night)

	if (uu>=vv) hn=PI; // polar day
	else if (uu<=-vv) hn=0.0; // polar night
	else hn=acos(-uu/vv); // Eqn 25

	//	Calculate total EET for this day

	eet_day=2.0*(s/(s+gamma)/lambda)*(uu*hn+vv*sin(hn))*K; // Eqn 26

	// Transfer output to member variables of Climate object

	climate.rad=rs_day;
	climate.eet=eet_day;
	climate.daylength=climate.daylength_save[date.day];

	// Calculate PAR from radiation
	// Eqn A1, Haxeltine & Prentice 1996

	climate.par=rs_day*FRADPAR;
}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Carslaw, HS & Jaeger JC 1959 Conduction of Heat in Solids, Oxford University
//   Press, London
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Henderson-Sellers, A & Robinson, PJ 1986 Contemporary Climatology. Longman,
//   Essex.
// Jarvis, PG & McNaughton KG 1986 Stomatal control of transpiration: scaling up
//   from leaf to region. Advances in Ecological Research 15: 1-49
// Jury WA, Gardner WR & Gardner WH 1991 Soil Physics 5th ed, John Wiley, NY
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65:
//   51-70.
// Press, WH, Teukolsky, SA, Vetterling, WT & Flannery, BT. (1986) Numerical
//   Recipes in FORTRAN, 2nd ed. Cambridge University Press, Cambridge
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Monteith, JL & Unsworth, MH 1990 Principles of Environmental Physics, 2nd ed,
//   Arnold, London
// van Duin, RHA 1963 The influence of soil management on the temperature
//   wave near the surface. Tech Bull 29 Inst for Land and Water Management
//   Research, Wageningen, Netherlands
