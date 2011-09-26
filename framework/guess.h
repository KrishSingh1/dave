///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK HEADER FILE
//
// Framework:             LPJ-GUESS Combined Modular Framework
//                        *************************************************************
//                        "Fast" version, revised December 2002 by Ben Smith
//                        Modified according to code changes by Dieter Gerten 021216
//                        *************************************************************
// Header file name:      guess.h
// Source code file name: guess.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
// Updated:               2010-11-22


// WHAT SHOULD THIS FILE CONTAIN?
// Framework header files should contain:
//   (1) definitions of all classes used by the framework and modules. Modules may
//       require classes to contain certain member variables and functions (see module
//       source files for details).
//   (2) other type, constant and function definitions to be accessible throughout the
//       model code.
//   (3) a forward declaration of the framework function if this is not the main
//       function.


///////////////////////////////////////////////////////////////////////////////////////
// #INCLUDES FOR LIBRARY HEADER FILES
// C/C++ libraries required for member functions of classes defined in this file.
// These libraries will also be available globally (so omit these #includes from source
// files). In addition to various standard C/C++ runtime libraries, the framework
// requires the following libraries (individual modules may use additional libraries)
//
// GUTIL
//   Includes class xtring, providing functionality for pointer-free dynamic handling
//   of character strings; wherever possible in LPJ-GUESS, strings are represented as
//   objects of type xtring rather than simple arrays of type char. GUTIL also provides
//   templates for dynamic collection classes (list arrays of various types), argument
//   processing for printf-style functions, timing functions and other utilities.

#ifndef LPJ_GUESS_GUESS_H
#define LPJ_GUESS_GUESS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gutil.h>


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

typedef enum {NOLIFEFORM,TREE,GRASS,CROP} lifeformtype;
	// Life form class for PFTs (trees, grasses)

typedef enum {NOPHENOLOGY,EVERGREEN,RAINGREEN,SUMMERGREEN,ANY} phenologytype;
	// Phenology class for PFTs

typedef enum {NOPATHWAY,C3,C4} pathwaytype;
	// Biochemical pathway for photosynthesis (C3 or C4)

typedef enum {NOINSOL,SUNSHINE,NETSWRAD,SWRAD} insoltype;
	// Units for insolation driving data (percentage sunshine, net instantaneous
	// downward shortwave radiation flux [W/m2], total [i.e. with no correction for
	// surface albedo] instantaneous downward shortwave radiation flux [W/m2])

typedef enum {NOVEGMODE,INDIVIDUAL,COHORT,POPULATION} vegmodetype;
	// Vegetation 'mode', i.e. what each Individual (see below) object represents;
	// either: (1) the average characteristics of all individuals comprising a PFT
	// population over the modelled area (standard LPJ mode); (2) a cohort of
	// individuals of a PFT that are roughly the same age; (3) an individual plant.

// GUESSN
typedef enum {SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFCWD,
	SOILMETA,SLOWSOM,PASSIVESOM,LEACHED,NSOMPOOL} pooltype;	
	// CENTURY pool names, NSOMPOOL number of SOM pools 
// end GUESSN

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

const int NSOILLAYER=2;
	// number of soil layers modelled

// FACE DAVID ORNL
//const double SOILDEPTH_UPPER=500.0; // soil upper layer depth (mm)
//const double SOILDEPTH_LOWER=1500.0;// soil lower layer depth (mm)

// FACE DAVID Duke
//const double SOILDEPTH_UPPER=400.0; // soil upper layer depth (mm)
//const double SOILDEPTH_LOWER=100.0; // soil lower layer depth (mm)

// NORMAL VALUES
const double SOILDEPTH_UPPER=500.0; // soil upper layer depth (mm)
const double SOILDEPTH_LOWER=1000.0; // soil lower layer depth (mm)


	// guess2008 - new default SOM values
const int SOLVESOM_END=400;
	// year at which to calculate equilibrium soil carbon
const int SOLVESOM_BEGIN=350;
	// year at which to begin documenting means for calculation of equilibrium
	// soil carbon
const double LAMBERTBEER_K=0.50;
	// Lambert-Beer extinction coefficient (Prentice et al 1993; Monsi & Saeki 1953)
const int NYEARGREFF=5; 
	// number of years to average growth efficiency over in function mortality
const int COLDEST_DAY_NHEMISPHERE=14;
	// day at which to start counting GDD's and leaf-on days for summergreen phenology
	// in N hemisphere (January 15)
const int COLDEST_DAY_SHEMISPHERE=195;
	// day at which to start counting GDD's and leaf-on days for summergreen phenology
	// in S hemisphere (July 15)
const int OUTPUT_MAXAGECLASS=40;
	// maximum number of age classes in age structure plots produced by function
	// outannual

// FACE DAVID climate David
const int NYEAR_NDEP = 258;			// Number of years with ndep data
const double KgTOg = 1000;			// Convert kg to g
const bool ifplantation=true; // not used yet
const bool ifdisturb_init=true; // disturbance during spin up to help trees establishing in
	// competition with grasses
const int distyear=50;

	// guess2008 - this is now a global, constant variable Previously, we had duplicate definitions in 
	// both canexch.cpp and soilwater.cpp
const double PRIESTLEY_TAYLOR=1.32;
	// Priestley-Taylor coefficient (conversion factor from equilibrium
	// evapotranspiration to PET)


///////////////////////////////////////////////////////////////////////////////////////
// FORWARD DECLARATIONS OF CLASSES DEFINED IN THIS FILE
// Forward declarations of classes used as types (e.g. for reference variables in some
// classes) before they are actually defined

class Date;
class Stand;
class Patch;
class Vegetation;


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are defined in the framework source code file, and are accessible
// throughout the code

extern Date date; // object describing timing stage of simulation
extern vegmodetype vegmode;
	// vegetation mode (population, cohort or individual)
extern int npatch;
	// number of patches in each stand (should always be 1 in population mode)
extern double patcharea;
	// patch area (m2) (individual and cohort mode only)
extern bool ifdailynpp;
	// whether photosynthesis calculations performed daily (alt: monthly)
extern bool ifdailydecomp;
	// whether soil decomposition calculations performed daily (alt: monthly)
extern bool ifbgestab;
	// whether background establishment enabled (individual, cohort mode)
extern bool ifsme;
	// whether spatial mass effect enabled for establishment (individual, cohort mode)
extern bool ifstochestab;
	// whether establishment stochastic (individual, cohort mode)
extern bool ifstochmort;
	// whether mortality stochastic (individual, cohort mode)
extern bool iffire; // whether fire enabled
extern bool ifdisturb;
	// whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
extern double distinterval;
	// generic patch-destroying disturbance interval (individual, cohort mode)
extern bool ifcalcsla; // whether SLA calculated from leaf longevity (alt: prescribed)
extern int estinterval; // establishment interval in cohort mode (years)
extern int npft; // number of possible PFTs
extern bool iffast; // whether to run in "fast" mode
extern bool ifcdebt; // whether C debt (storage between years) permitted

// GUESSN
extern bool ifcentury;
	// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
extern bool ifnlim;
	// whether plant growth limited by available N
extern int freenyears;
	// number of years to allow spinup without N limitation
extern bool iflimvmax;
	// whether Vmax limited by leaf N content
extern double nrelocfrac;
	// fraction of N relocated by plants from roots and leaves
extern bool ifvarycn;
	// whether leaf and tissue C:N ratios are adjusted according to photosynthetic
	// demand (i.e. Vmax)
extern bool ifnlimvarycn;
	// whether leaf and tissue C:N ratios are adjusted according to N limitation
extern double cwdtransfer;
	// fraction of woody debris transferred to SOM each year
extern double nmass_avail_max;
	// max N:C ratio in the soil (should be 0.002 (Parton et al 1993, Fig. 4))
extern bool ifleachn;
	// whether to allow N leaching
extern bool ifindiv_fnuptake;
	// whether to allow individual fractional N uptake
extern int ifnfix;
	// whether to include an estimate for N fixation
extern bool ifndepdata;
	// whether N deposition data available from a file
extern double andep;
	// annual N deposition (used only if ifndepdata=false)
extern double minndep;
	// minimum annual N deposition
extern bool ifdailysetntoc;
	// if to use daily version of setntoc (set N:C ratio of som pools)
extern bool ifnstorage;
	// if to use a N storage for each individual
extern bool ifndemand_new_est;
	// if to use N limitation on new establishment 
// end GUESSN

// FACE David
extern int FYEAR_SCENARIO_FACE;
extern int FACE_ring;				// Which FACE ring examined
extern int ifduke;					// Duke or Oak Ridge
extern int has_FACE_clim;			// If we are using FACE clim data

// CANIF David
extern int has_CANIF_clim;			// If we are using CANIF clim data
extern int WSITE;

///////////////////////////////////////////////////////////////////////////////////////
// guess2008 - new input variables, from the .ins file
extern bool ifsmoothgreffmort;				
	// whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard 
	// step-function (0)
extern bool ifdroughtlimitedestab;			
	// whether establishment is limited by growing season drought 
extern bool ifrainonwetdaysonly;			
	// rain on wet days only (1, true), or a little every day (0, false); 
extern bool ifspeciesspecificwateruptake;	
	// whether water uptake is species specific 



///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL FUNCTION DECLARATIONS
// These functions are defined in the framework source file or main module (if
// a separate main module is implemented), and are accessible throughout the code

void dprintf(xtring format,...);
	// a printf-style function for text output to the screen and/or a log file. To
	// maintain portability of the modular code, please use this function for general
	// output instead of the standard C++ printf function
void fail(xtring format,...);
	// a printf-style function that sends output to the screen and/or a log file, then
    // terminates execution.
void plot(xtring window_name,xtring series_name,double x,double y);
	// adds data point (x,y) to series 'series_name' of line graph 'window_name'. If
	// the series and/or line graph do not yet exist, they are created. Functional only
	// when the framework is built as a DLL and linked to the LPJ-GUESS Windows Shell
	// (the function may still be called in other implementations, but will have no
	// effect).
void resetwindow(xtring window_name);
	// 'forgets' series and data for line graph 'window_name' created using function
	// plot (above). Functional only when the framework is built as a DLL and
	// linked to the LPJ-GUESS Windows Shell.
void clear_all_graphs();
	// 'forgets' series and data for all currently-defined line graphs created using
	// function plot (above). Functional only when the framework is built as a DLL and
	// linked to the LPJ-GUESS Windows Shell.
bool abort_request_received();
	// May be called by framework to respond to abort request from Windows shell
	// (returns true if shell has sent an abort request, otherwise false)


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL FUNCTION DEFINITIONS
// Small inline ("macro") functions, accessible throughout the code

inline bool negligible(double dval) {

	// Returns true if dval < EPSILON, otherwise false

	const double EPSILON=1.0e-30;
	if (dval>EPSILON) return false;
	if (dval<0.0 && dval<-EPSILON) return false;
	return true;
}

inline bool equal(double dval1,double dval2) {

	// Returns true if |dval1-dval2| < EPSILON, otherwise false

	const double EPSILON=1.0e-30;
	if (dval1==dval2) return true;
	if (dval1>dval2) {
		if (dval1-dval2<EPSILON) return true;
	}
	else {
		if (dval2-dval1<EPSILON) return true;
	}
	return false;
}

inline double max(double dval1,double dval2) {

	// Returns the larger of dval1 and dval2

	if (dval1>dval2) return dval1;
	return dval2;
}

inline double min(double dval1,double dval2) {

	// Returns the smaller of dval1 and dval2

	if (dval1<dval2) return dval1;
	return dval2;
}


///////////////////////////////////////////////////////////////////////////////////////
// DATE
// General purpose object for handling simulation timing. In general, frameworks should
// use a single Date object for all simulation timing.
//
// Date provides the following functionality:
//
// Date()
//   Constructor function called automatically when Date object is created (do not call
//   explicitly). Initialises some member variables.
//
// void init(int nyearsim)
//   Call to initialise date to day 0 of year 0 and set intended number of simulation
//   years to nyearsim (used only to set islastyear flag - actual simulation may be
//   longer or shorter than nyearsim)
//
// void next()
//   Call at end of every simulation day to update member variables.
//
// int prevmonth()
//   Returns index (0-11) of previous month (11 if currently month 0).
//
// int nextmonth()
//   Returns index of next month (0 if currently month 11)
//
// Member variables of the class (see below) provide various kinds of calender and
// timing information, assuming init has been called to initialise the object, and
// next() has been called at the end of each simulation day.

class Date {

	// MEMBER VARIABLES

public:

	int ndaymonth[12];
		// number of days in each month (0=January - 11=December)
	int day;
		// julian day of year (0-364; 0=Jan 1)
	int dayofmonth;
		// day of current month (0=first day)
	int month;
		// month number (0=January - 11=December)
	int year;
		// year since start of simulation (0=first simulation year)
	int middaymonth[12];
		// julian day for middle day of each month
	bool islastyear;
		// true if last year of simulation, false otherwise
	bool islastmonth;
		// true if last month of year, false otherwise
	bool islastday;
		// true if last day of month, false otherwise
	bool ismidday;
		// true if middle day of month, false otherwise

private:

	int nyear;

	// MEMBER FUNCTIONS

public:
	
	Date() {
		const int data[]={31,28,31,30,31,30,31,31,30,31,30,31};
		int month;
		int dayct=0;
		for (month=0;month<12;month++) {
			ndaymonth[month]=data[month];
			middaymonth[month]=dayct+data[month]/2;
			dayct+=data[month];
		}
	}
	
	void init(int nyearsim)	{
		nyear=nyearsim;
		day=month=year=dayofmonth=0;
		islastmonth=islastday=ismidday=false;
		if (nyear==1) islastyear=true;
		else islastyear=false;
	}

	void next() {
		if (islastday) {
			if (islastmonth) {
				dayofmonth=0;
				day=0;
				month=0;
				year++;
				if (year==nyear-1) islastyear=true;
				islastmonth=false;
			}
			else {
				day++;
				dayofmonth=0;
				month++;
				if (month==11) islastmonth=true;
			}
			islastday=false;
		}
		else {
			day++;
			dayofmonth++;
			if (dayofmonth==ndaymonth[month]/2) ismidday=true;
			else {
				ismidday=false;
				if (dayofmonth==ndaymonth[month]-1) islastday=true;
			}
		}
	}

	int prevmonth() {
		if (month>0) return month-1;
		return 11;
	}

	int nextmonth() {
		if (month<11) return month+1;
		return 0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// CLIMATE
// Stores all static and variable data relating to climate parameters, as well as 
// latitude, atmospheric CO2 concentration and daylength for a stand (corresponding to
// a modelled locality or grid cell). Includes a reference to the parent Stand object
// (defined below). Initialised by a call to initdrivers.

class Climate {

	// MEMBER VARIABLES

public:
	Stand& stand;
		// reference to parent Stand object
	double temp;
		// mean air temperature today (deg C)
	double rad;
		// total daily net downward shortwave solar radiation today (J/m2/day)
	double par;
		// total daily photosynthetically-active radiation today (J/m2/day)
	double prec;
		// precipitation today (mm)
	double daylength;
		// day length today (h)
	double co2;
		// atmospheric ambient CO2 concentration today (ppmv)
	double lat;
		// latitude (degrees; +=north, -=south)
	double insol;
		// insolation today
	insoltype instype;
		// units in which insol expressed:
		// SUNSHINE = percentage of full sunshine
		// NETSWRAD = net downward shortwave radiation flux (albedo corrected) (W/m2)
		// SWRAD    = total downward shortwave radiation flux (W/m2)
	double eet;
		// equilibrium evapotranspiration today (mm/day)
	double mtemp;
		// mean temperature for the last 31 days (deg C)
	double mtemp_min20;
		// lowest mean monthly temperature for the last 20 years (deg C)
	double mtemp_max20;
	double mtemp_max;
		// highest mean monthly temperature for the last 12 months (deg C)
	double gdd5;
		// accumulated growing degree day sum on 5 degree base (reset when temperatures
		// fall below 5 deg C)
	double agdd5; // total gdd5 (accumulated) for this year (reset 1 January)
	int chilldays;
		// number of days with temperatures <5 deg C (reset when temperatures fall
		// below 5 deg C; maximum value 365)
	bool ifsensechill;
		// guess2008 - CHILLDAYS - true if chill day count may be reset by temperature 
		// fall below 5 deg C
	double gtemp;
		// respiration response to today's air temperature incorporating damping of Q10
		// due to temperature acclimation (Lloyd & Taylor 1994)
	int last_gtemp;
		// the last day (0-364) for which gtemp was calculated
	double mgtemp;
		// gtemp (see above) calculated for this month's average temperature
	int last_mgtemp;
		// the last month (0-11) for which mgtemp was calculated
	double dtemp_31[31];
		// daily temperatures for the last 31 days (deg C)
	double mtemp_min_20[20];
		// minimum monthly temperatures for the last 20 years (deg C)
	double mtemp_max_20[20];
	double mtemp_min;
		// minimum monthly temperature for the last 12 months (deg C)
	double atemp_mean;
		// mean of monthly temperatures for the last 12 months (deg C)

	// GUESSN
	double andep_1860;
		// annual nitrogen deposition (kgN/m2/year)
	double andep_1993;
		// annual nitrogen deposition (kgN/m2/year)
	double andep_2050;
		// annual nitrogen deposition (kgN/m2/year)
	double andep;
		// annual nitrogen deposition (kgN/m2/year)
	// end GUESSN

	// Monthly sums (converted to means) used by canopy exchange module

	double temp_mean;
		// accumulated mean temperature for this month (deg C)
	double par_mean;
		// accumulated mean daily net PAR sum (J/m2/day) for this month
	double co2_mean;
		// accumulated mean CO2 for this month (ppmv)
	double daylength_mean;
		// accumulated mean daylength for this month (h)

	// Saved parameters used by function daylengthinsolpet

	double sinelat;
	double cosinelat;
	double qo[365],u[365],v[365],hh[365],sinehh[365];
	double daylength_save[365];
	bool doneday[365];
		// indicates whether saved values exist for this day

	// MEMBER FUNCTIONS

public:

	Climate(Stand& s):stand(s) {};
		// constructor function: initialises stand member

	void initdrivers(double latitude) {

		// Initialises certain member variables
		// Should be called before Climate object is applied to a new grid cell

		const double DEGTORAD=0.01745329;
		int day,year;

		for (year=0;year<20;year++) {
			mtemp_min_20[year]=0.0;
			mtemp_max_20[year]=0.0;
		}
		mtemp=0.0;
		gdd5=0.0;
		chilldays=0;
		ifsensechill=true; //  guess2008 - CHILLDAYS
		atemp_mean=0.0;
		last_gtemp=-1;
		last_mgtemp=-1;

		lat=latitude;
		for (day=0;day<365;day++) doneday[day]=false;
		sinelat=sin(lat*DEGTORAD);
		cosinelat=cos(lat*DEGTORAD);
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// FLUXES
// Stores daily and accumulated annual fluxes (currently only C fluxes). Upward fluxes
// (from vegetation to atmosphere or soil to atmosphere) are positive values, downward
// fluxes (from atmosphere to vegetation) are negative values. Accumulated fluxes
// should be initialised where appropriate in the model code - initialisation is not
// provided as part of the class functionality. One Fluxes object is defined for each
// patch (see below).

class Fluxes {

	// MEMBER VARIABLES
	// (all fluxes on stand area basis, kgC/m2)

public:

	Patch& patch;
		// reference to patch to which this Fluxes object belongs
	double acflux_veg;
		// annual flux to vegetation (=total vegetation annual NPP)
	double acflux_fire;
		// annual carbon flux to atmosphere from burnt vegetation and litter
	double acflux_soil;
		// annual carbon flux to atmosphere from soil respiration
	double acflux_est;
		// annual flux from atmosphere to vegetation associated with establishment
	double dcflux_soil;
		// daily carbon flux to atmosphere from soil respiration
		// NB: not implemented by som_dynamics_monthly
	double mcflux_soil[12];
		// monthly C flux to atmosphere from soil respiration
	double mcflux_veg[12];
		// monthly C flux to vegetation from atmosphere
	double dcflux_veg;
		// daily net carbon flux to vegetation (respiration-assimilation)
		// NB: not implemented by canopy_exchange_monthly

	// guess2008 - new carbon budget arrays
	double mcflux_gpp[12];
		// monthly GPP
	double mcflux_ra[12];
		// monthly autotrophic respiration
	double dcflux_gpp[365];
		// daily GPP

	// MEMBER FUNCTIONS

public:

	Fluxes(Patch& p):patch(p) {}
		// constructor: initialises patch member

	double anee() {
		
		// If called following update of annual accumulated fluxes on last day of
		// simulation year, returns annual net ecosystem exchange (NEE)

		return acflux_veg+acflux_fire+acflux_soil+acflux_est;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// PFT
// Holds static functional parameters for a plant functional type (PFT). There should
// be one Pft object for each potentially occurring PFT. The same Pft object may be
// referenced (via the pft member of the Individual object; see below) by different
// average individuals. Member functions are included for initialising SLA given leaf
// longevity, and for initialising sapling/regen characteristics (required for
// population mode).

class Pft {

	// MEMBER VARIABLES

public:
	int id;
		// id code (should be zero based and sequential, 0...npft-1)
	xtring name;
		// name of PFT
	lifeformtype lifeform;
		// life form (tree or grass)
	phenologytype phenology;
		// leaf phenology (raingreen, summergreen, evergreen, rain+summergreen)
	double phengdd5ramp;
		// growing degree sum on 5 degree base required for full leaf cover
	double wscal_min;
		// water stress threshold for leaf abscission (range 0-1; raingreen PFTs)
	pathwaytype pathway;
		// biochemical pathway for photosynthesis (C3 or C4)
	double pstemp_min;
		// approximate low temperature limit for photosynthesis (deg C)
	double pstemp_low;
		// approximate lower range of temperature optimum for photosynthesis (deg C)
	double pstemp_high;
		// approximate upper range of temperature optimum for photosynthesis (deg C)
	double pstemp_max;
		// maximum temperature limit for photosynthesis (deg C)
	double lambda_max;
		// non-water-stressed ratio of intercellular to ambient CO2 partial pressure
	double rootdist[NSOILLAYER];
		// vegetation root profile (array containing fraction of roots in each soil
		// layer, [0=upper layer])
	double gmin;
		// canopy conductance component not associated with photosynthesis (mm/s)
	double emax;
		// maximum evapotranspiration rate (mm/day)
	double respcoeff;
		// maintenance respiration coefficient (0-1)
	// GUESSN
	double cton_leaf_min;
		// minimum leaf C:N mass ratio
	double cton_leaf_max;
		// maximum leaf C:N mass ratio
	double cton_leaf_avr;
		// average leaf C:N mass ratio
	double cton_root_avr;
		// average fine root C:N mass ratio
	double cton_sap_avr;
		// average sapwood C:N mass ratio
	double n_reserve;
		// N storage organ in relation to sapwood carbon
	// end GUESSN
	double reprfrac;
		// fraction of NPP allocated to reproduction
	double turnover_leaf;
		// annual leaf turnover as a proportion of leaf C biomass
	double turnover_root;
		// annual fine root turnover as a proportion of fine root C biomass
	double turnover_sap;
		// annual sapwood turnover as a proportion of sapwood C biomass
	double wooddens;
		// sapwood and heartwood density (kgC/m3)
	double crownarea_max;
		// maximum tree crown area (m2)
	double k_allom1;
		// constant in allometry equations
	double k_allom2;
		// constant in allometry equations
	double k_allom3;
		// constant in allometry equations
	double k_rp;
		// constant in allometry equations
	double k_latosa;
		// tree leaf to sapwood area ratio
	double sla;
		// specific leaf area (m2/kgC)
	double leaflong;
		// leaf longevity (years)
	double ltor_max;
		// leaf to root mass ratio under non-water-stressed conditions
	double litterme;
		// litter moisture flammability threshold (fraction of AWC)
	double fireresist;
		// fire resistance (0-1)
	double parff_min;
		// minimum forest-floor PAR level for growth (grasses) or establishment (trees)
		// (J/m2/day) (individual and cohort modes)
	double alphar;
		// parameter capturing non-linearity in recruitment rate relative to
		// understorey growing conditions for trees (Fulton 1991) (individual and
		// cohort modes)
	double est_max;
		// maximum sapling establishment rate (saplings/m2/year) (individual and cohort
		// modes)
	double kest_repr;
		// constant used in calculation of sapling establishment rate when spatial
		// mass effect enabled (individual and cohort modes)
	double kest_bg;
		// constant affecting amount of background establishment (when enabled)
		// (individual and cohort modes)
	double kest_pres;
		// constant used in calculation of sapling establishment rate when spatial
		// mass effect disabled (individual and cohort modes)
	double longevity;
		// expected longevity under non-stressed conditions (individual and cohort
		// modes)
	double greff_min;
		// threshold growth efficiency for imposition of growth suppression mortality
		// (kgC/m2 leaf/year) (individual and cohort modes)

	// Bioclimatic limits (all temperatures deg C)

	double tcmin_surv;
		// minimum 20-year coldest month mean temperature for survival
	double tcmax_est;
		// maximum 20-year coldest month mean temperature for establishment
	double gdd5min_est;
		// minimum degree day sum on 5 deg C base for establishment
	double tcmin_est;
		// minimum 20-year coldest month mean temperature for establishment
	double twmin_est;
		// minimum warmest month mean temperature for establishment
	double twminusc;
		// continentality parameter for boreal summergreen trees
	double k_chilla;
		// constant in equation for budburst chilling time requirement (Sykes et al 1996)
	double k_chillb;
		// coefficient in equation for budburst chilling time requirement
	double k_chillk;
		// exponent in equation for budburst chilling time requirement
	double gdd0[366];
		// array containing values for GDD0(c) given c=number of chill days (0-365)
		// (Sykes et al 1996, Eqn 1)
	double intc;
		// interception coefficient (unitless)

	// guess2008 - drought-limited establishment (DLE)
	double drought_tolerance;
		// Drought tolerance level (0 = very -> 1 = not at all) (unitless)

	// Sapling/regeneration characteristics (used only in population mode):
	// for trees, on sapling individual basis (kgC); for grasses, on stand area basis,
	// kgC/m2

	struct {
		double cmass_leaf;
			// leaf C biomass
		double cmass_root;
			// fine root C biomass
		double cmass_sap;
			// sapwood C biomass
		double cmass_heart;
			// heartwood C biomass
	} regen;

	// Variables used by new hydrology (Dieter Gerten 2002-07)
	
	// MEMBER FUNCTIONS

public:

	Pft() {
		
		// Constructor (initialises array gdd0)
		
		int y;
		for (y=0;y<366;y++)
			gdd0[y]=-1.0; // value<0 signifies "unknown"; see function phenology()

		// guess2008 - DLE
		drought_tolerance=0.0; // Default, means that the PFT will never be limited by drought.
	}

	void initsla() {

		// Calculates SLA given leaf longevity
		// Reich et al 1997, Fig 1f (includes conversion x2.0 from m2/kg_dry_weight to
		// m2/kgC)

		sla=0.2*exp(6.15-0.46*log(leaflong*12.0));
	}

	void initregen() {

		// Initialises sapling/regen characteristics in population mode
		// following LPJF formulation; see function allometry in growth module.
		// Note: primary PFT parameters, including SLA, must be set before this
		//       function is called
	
		const double PI=3.14159265;
		const double REGENLAI_TREE=1.5;
		const double REGENLAI_GRASS=0.001;
		const double SAPLINGHW=0.2;

		if (lifeform==TREE) {

			// Tree sapling characteristics

			regen.cmass_leaf=pow(REGENLAI_TREE*k_allom1*pow(1.0+SAPLINGHW,k_rp)*
				pow(4.0*sla/PI/k_latosa,k_rp*0.5)/sla,2.0/(2.0-k_rp));

			regen.cmass_sap=wooddens*k_allom2*pow((1.0+SAPLINGHW)*
				sqrt(4.0*regen.cmass_leaf*sla/PI/k_latosa),k_allom3)*
				regen.cmass_leaf*sla/k_latosa;

			regen.cmass_heart=SAPLINGHW*regen.cmass_sap;
		}
		else if (lifeform==GRASS || lifeform==CROP) {

			// Grass regeneration characteristics

			regen.cmass_leaf=REGENLAI_GRASS/sla;
		}

		regen.cmass_root=1.0/ltor_max*regen.cmass_leaf;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// PFTLIST
// Functionality for building, maintaining, referencing and destroying a list array of
// Pft objects. In general, frameworks should define a single Pftlist object containing
// a single list of PFTs. Pft objects within the list are then referenced by the pft
// member of each Individual object.
//
// Functionality is inherited from the ListArray_id template type in the GUTIL
// Library. Sequential Pft objects can be referenced as array elements by id:
//
//   Pftlist pftlist;
//   ...
//   for (i=0; i<npft; i++) {
//     Pft& thispft=pftlist[i];
//     /* query or modify object thispft here */
//   }
//
// or by iteration through the linked list:
//
//   pftlist.firstobj();
//   while (pftlist.isobj) {
//     Pft& thispft=pftlist.getobj();
//     /* query or modify object thispft here */
//     pftlist.nextobj();
//   }

class Pftlist : public ListArray_id<Pft> {};



///////////////////////////////////////////////////////////////////////////////////////
// INDIVIDUAL
// State variables for a vegetation individual. In population mode this is the average
// individual of a PFT population; in cohort mode: the average individual of a cohort;
// in individual mode: an individual plant. Each grass PFT is represented as a single
// individual in all modes. Individual objects are collected within list arrays of
// class Vegetation (defined below), of which there is one for each patch, and include
// a reference to their 'parent' Vegetation object. Use the createobj member function
// of class Vegetation to add new individuals.

class Individual {

public:
	Pft& pft;
		// reference to Pft object containing static parameters for this individual
	Vegetation& vegetation;
		// reference to Vegetation object to which this Individual belongs
	int id;
		// id code (0-based, sequential)
	double cmass_leaf;
		// leaf C biomass on modelled area basis (kgC/m2)
	double cmass_root;
		// fine root C biomass on modelled area basis (kgC/m2)
	double cmass_sap; 
		// sapwood C biomass on modelled area basis (kgC/m2)
	double cmass_heart;
		// heartwood C biomass on modelled area basis (kgC/m2)
	double cmass_debt;
		// C "debt" (retrospective storage) (kgC/m2)

	// GUESSN
	double nmass_leaf;
		// N content of leaves on patch area basis (kgN/m2)
	double nmass_root;
		// N content of roots on patch area basis (kgN/m2)
	double nmass_sap;
		// N content of sapwood on patch area basis (kgN/m2)
	double nmass_heart;
		// N content of heartwood on patch area basis (kgN/m2)
	double nmass_reserve;
		// N content of storage on patch area basis (kgN/m2)
	// end GUESSN

	double fpc;
		// foliar projective cover (FPC) under full leaf cover as fraction of modelled
		// area
	double fpar;
		// fraction of PAR absorbed by foliage over projective area today, taking
		// account of leaf phenological state
	double densindiv;
		// average density of individuals over patch (indiv/m2)
	double phen;
		// vegetation phenological state (fraction of potential leaf cover)
	double aphen;
		// annual sum of daily fractional leaf cover (equivalent number of days with
		// full leaf cover) (population mode only; reset on expected coldest day of
		// year)
	int aphen_raingreen;
		// annual number of days with full leaf cover) (raingreen PFTs only; reset on
		// 1 January)
	double assim;
		// daily net assimilation (GPP-leaf respiration) on modelled area basis
		// (kgC/m2/day)
	double resp;
		// daily maintenance respiration (not including leaf respiration) and growth
		// respiration on modelled area basis (kgC/m2/day)
	double anpp;
		// accumulated NPP over modelled area (kgC/m2/year); = annual NPP following
		// call to growth module on last day of simulation year
	double aet;
		// actual evapotranspiration over projected area (mm/day)
	double ltor;
		// leaf to root mass ratio
	double height;
		// plant height (m)
	double crownarea;
		// plant crown area (m2)
	double deltafpc;
		// increment in fpc since last simulation year
	double wscal;
		// water stress parameter (0-1 range; 1=minimum stress) (updated daily)
	double wscal_mean;
		// running sum (converted to annual mean) for wscal
	double boleht;
		// bole height, i.e. height above ground of bottom of crown cylinder (m)
		// (individual and cohort modes only)
	double lai;
		// patch-level lai for this individual or cohort (function fpar)
	double lai_layer;
		// patch-level lai for cohort in current vertical layer (function fpar)
	double lai_indiv;
		// individual leaf area index (individual and cohort modes only)
	double greff_5[NYEARGREFF];
		// growth efficiency (NPP/leaf area) for each of the last five simulation years
		// (kgC/m2/yr)
	double age;
		// individual/cohort age (years)
	double mnpp[12];
		// monthly NPP (kgC/m2/month)
	double mlai[12];
		// monthly LAI (including phenology component)

	double mgpp[12];
		// monthly GPP-leafresp (kgC/m2/month)
	double mra[12];
		// monthly respiration

	// Variables used by "fast" canopy exchange code (Ben Smith 2002-07)

	double fpar_wstress;
		// FPAR for days with water stress (see canopy exchange module)
	double fpar_leafon;
		// FPAR assuming full leaf cover for all vegetation
	double lai_leafon_layer;
		// LAI for current layer in canopy (cohort/individual mode; see function fpar)
	double gp_leafon;
		// non-water-stressed canopy conductance on FPC basis (mm/s)
	double demand;
		// transpirative demand on FPC basis (mm/day)
	double demand_leafon;
		// transpirative demand assuming full leaf cover on FPC basis (mm/day)
	double supply;
		// supply function of AET, FPC basis (mm/day)
	double supply_leafon;
		// supply function of AET assuming full leaf cover, FPC basis (mm/day)
	double intercep;
		// interception associated with this individual today (patch basis)


	// Monthly sums (converted to means) maintained by function canopy_exchange

	double phen_mean;
		// accumulated mean fraction of potential leaf cover

	// Means for driving parameters of photosynthesis required for "individual" demand
	// mode (see canexch.cpp)

	double temp_wstress; // temperature (deg C)
	double par_wstress; // PAR (J/m2/day)
	double daylength_wstress; // daylength (h)
	double co2_wstress; // CO2 (ppmv)
	int nday_wstress; // number of water-stress days for month
	bool ifwstress; // whether individual subject to water stress today

	// GUESSN

	double nstore;
		// relocated N from leaves and roots and accumulated uptake from soil mineral N pool
	double leafn;
		// cumulative mean (calculated at end of each month) of daily leaf N (kgN/m2)
		// (leaf N demand calculated from Vmax)
	double leafn_max;
		// largest monthly value of leafn for year
	double leafn_mean;
		// mean monthly value of leafn for year
	double ndemand;
		// annual N demand (used in growth)
	double ndemand_uptake;
		// annual N demand (used in vegetation_n_uptake)
	double fnuptake;
		// fractional N uptake of indiv demand
	double n_reserve_uptake;
		// fraction extra N uptake to reserve pool
	double max_n_reserve;
		// maximum size of N reserve
	double max_n_reserve_old;
	double limnfact;
		// actual fractional N available to indiv N demand in Growth()
	double na_fpar;
		// leaf N associated with photosynthesis 
	double assim_nowstress;
		// saved assimilation in case it turns out to be a non-water-stress day
	int nday_leafon;	
		// Number of days with non-negligible phenology this month
	double dassim[365];
		// daily net assimilation (kgC/m2/yr) - used by SOM dynamics to distribute
		// plant N uptake through the year
	double aassim;
		// annual sum of positive dassim (above) - used by SOM dynamics
	double vmax_lim[365];
		// daily N limitation to vmax
	double nopt;
		// N limitation on vmax
	double cton_leaf_new;
		// C:N ratio for new biomass (leaf)
	double cton_root_new;
		// C:N ratio for new biomass (root)
	double cton_sap_new;
		// C:N ratio for new biomass (sap)
	double cton_leaf_old;
		// C:N ratio of old (current) biomass (leaf)
	double cton_root_old;
		// C:N ratio of old (current) biomass (root)
	double cton_sap_old;
		// C:N ratio of old (current) biomass (sap)
	double cton_leaf_opt;
		// optimal (photosynthesis) C:N ratio for new biomass (leaf) 

	double dnupnpp;	// Daily N uptake variables	guessnfix 
	double nstore_daily;
	double bminc_leaf_frac;	
	double bminc_root_frac;
	
	// end GUESSN

	bool alive; 
		// guess2008 - whether this individual is truly alive. Set to false for first year 
		// after the Individual object is created, then true.


	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v):id(i),pft(p),vegetation(v) {

		anpp=0.0;
		fpc=0.0;
		densindiv=0.0;
		cmass_leaf=0.0;
		cmass_root=0.0;
		cmass_sap=0.0;
		cmass_heart=0.0;
		cmass_debt=0.0;
		wscal=1.0;
		phen=0.0;
		aphen=0.0;
		deltafpc=0.0;
		fpar_wstress=0.0;
		assim=0.0;
		assim_nowstress=0.0;

		// GUESSN
		nmass_leaf=0.0;
		nmass_root=0.0;
		nmass_sap=0.0;
		nmass_heart=0.0;
		nmass_reserve=0.0;

		nstore=0.0;
		fnuptake=1.0;
		n_reserve_uptake=0.0;
		max_n_reserve=0.0;

		// end GUESSN
	
		// guess2008 - additional initialisation
		age=0.0;
		fpar=0.0;
		aphen_raingreen=0;
		demand=0.0;
		supply=0.0;
		intercep=0.0;
		phen_mean=0.0;
		temp_wstress = 0.0;
		par_wstress = 0.0;
		daylength_wstress = 0.0;
		co2_wstress = 0.0; 
		nday_wstress = 0; 
		ifwstress = false;
		lai = 0.0;
		lai_layer = 0.0;
		lai_indiv = 0.0;
		alive = false;

		int m;
		for (m=0;m<12;m++) {
			mnpp[m]=mlai[m]=mgpp[m]=mra[m]=0.0;
		}
		int d;
		for (d=0;d<365;d++) {
			dassim[d]=vmax_lim[d]=0.0;
		}
	};
};


///////////////////////////////////////////////////////////////////////////////////////
// VEGETATION
// Functionality for building, maintaining, referencing and destroying a list array of
// Individual objects. A single Vegetation object is defined for each patch. A
// reference to the parent Patch object (defined below) is included as a member
// variable.
//
// Functionality is inherited from the ListArray_idin1 template type in the GUTIL
// Library. Sequential Individual objects can be referenced as array elements by id,
// or by iteration through the linked list:
//
//   Vegetation vegetation
//   ...
//   vegetation.firstobj();
//   while (vegetation.isobj) {
//     Individual& thisindiv=vegetation.getobj();
//     /* query or modify object thisindiv here */
//     vegetation.nextobj();
//   }

class Vegetation : public ListArray_idin2<Individual,Pft,Vegetation> {

public:
	// MEMBER VARIABLES

	Patch& patch; // reference to parent Patch object

	// MEMBER FUNCTIONS

	Vegetation(Patch& p):patch(p) {};
		// constructor (initialises member variable patch)
};


///////////////////////////////////////////////////////////////////////////////////////
// SOILTYPE
// Stores static parameters for soils and the snow pack. One Soiltype object is defined
// for each stand. State variables for soils are held by objects of class Soil, of
// which there is one for each patch (see below).

class Soiltype {

	// MEMBER VARIABLES

public:

	double awc_frac;
		// available water holding capacity as fraction of soil volume
	double awc[2];
		// available water holding capacity of soil layers [0=upper layer] (mm)
	double perc_base;
		// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice
		// 1996)
	double perc_exp;
		// exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice
		// 1996)
	double thermdiff_0;
		// thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_15;
		// thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_100;
		// thermal diffusivity at 100% WHC (mm2/s)
	double wp[2];
		// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double f_FC[2];
		// ratio between saturation capacity and field capacity. Cosby et al 1984
	int solvesom_end;
		// year at which to calculate equilibrium soil carbon
	int solvesom_begin;
		// year at which to begin documenting means for calculation of equilibrium
		// soil carbon

	// GUESSN: For CENTURY ...
	
	double sand_frac;
		// fraction of soil that is sand
	double clay_frac;
		// fraction of soil that is clay
	double silt_frac;
		// fraction of soil that is silt plus clay
	// end GUESSN

	// MEMBER FUNCTIONS

public:

	Soiltype() {

		// Constructor: initialises certain member variables

		solvesom_end=SOLVESOM_END;
		solvesom_begin=SOLVESOM_BEGIN;

		// GUESSN
		sand_frac=0.3;
		clay_frac=0.3;
		silt_frac=0.2;
		// end GUESSN
	}

	// guess2008 - override the default SOM years with 70-80% of the spin-up period length
	void updateSolveSOMvalues(const int& nyrspinup) {
		
		solvesom_end=0.8*nyrspinup;
		solvesom_begin=0.7*nyrspinup;

	}
};

///////////////////////////////////////////////////////////////////////////////////////
// CENTURY SOIL POOL (GUESSN)

class Sompool {

public:
	double cmass;
		// C mass in pool kgC/m2
	double nmass;
		// N mass in pool kgN/m2
	double cdec; // (potential) decrease in C following decomposition today (kgC/m2)
	double ndec; // (potential) decrease in N following decomposition today (kgN/m2)
	double delta_cmass,delta_nmass;
	double ligcfrac;
	double frc;
	double ntoc;

	void init() {
		
		// Initialise pool
		
		cmass=0.0;
		nmass=0.0;
		ligcfrac=0.0;
	};
};

// end GUESSN

///////////////////////////////////////////////////////////////////////////////////////
// SOIL
// Stores state variables for soils and the snow pack. Initialised by a call to
// initdrivers. One Soil object is defined for each patch. A reference to the parent
// Patch object (defined below) is included as a member variable. Soil static
// parameters are stored as objects of class Soiltype, of which there is one for each
// stand. A reference to the Soiltype object holding the static parameters for this
// soil is included as a member variable.

class Soil {

	// MEMBER VARIABLES

public:

	Patch& patch;
		// reference to parent Patch object
	Soiltype& soiltype;
		// reference to Soiltype object holding static parameters for this soil
	double wcont[NSOILLAYER];
		// water content of soil layers [0=upper layer] as fraction of available water
		// holding capacity;
	double awcont[NSOILLAYER];
		// guess2008 - DLE - the average wcont over the growing season, for each soil layer
	double wcont_evap;
		// water content of sublayer of upper soil layer for which evaporation from
		// the bare soil surface is possible (fraction of available water holding
		// capacity)
	double dwcontupper[365];
		// daily water content in upper soil layer for each day of year
	double mwcontupper;
		// mean water content in upper soil layer for last month
		// (valid only on last day of month following call to daily_accounting_patch)
	double snowpack;
		// stored snow as average over modelled area (mm rainfall equivalents)
	double runoff;
		// total runoff today (mm/day)
	double temp;
		// soil temperature today at 0.25 m depth (deg C)
	double dtemp[31];
		// daily temperatures for the last month (deg C)
		// (valid only on last day of month following call to daily_accounting_patch)
	double mtemp;
		// mean soil temperature for the last month (deg C)
		// (valid only on last day of month following call to daily_accounting_patch)
	double gtemp;
		// respiration response to today's soil temperature at 0.25 m depth
		// incorporating damping of Q10 due to temperature acclimation (Lloyd & Taylor
		// 1994)
	int last_gtemp;
		// the last day (0-364) for which gtemp was calculated
	double mgtemp;
		// gtemp (see above) calculated for this month's average temperature
	int last_mgtemp;
		// the last month (0-11) for which mgtemp was calculated
	double cpool_slow;
		// soil organic matter (SOM) pool with c. 1000 yr turnover (kgC/m2)
	double cpool_fast;
		// soil organic matter (SOM) pool with c. 33 yr turnover (kgC/m2)

	// Running sums (converted to long term means) maintained by SOM dynamics module

	double decomp_litter_mean;
		// mean annual litter decomposition (kgC/m2/yr)
	double k_soilfast_mean;
		// mean value of decay constant for fast SOM fraction
	double k_soilslow_mean;
		// mean value of decay constant for slow SOM fraction

	// Parameters used by function soiltemp and updated monthly

	double alag,exp_alag;


	// guess2008 - 3 new soil water variables
	double mwcont[12][NSOILLAYER];
		// water content of soil layers [0=upper layer] as fraction of available water
		// holding capacity;
	double dwcontlower[365];
		// daily water content in lower soil layer for each day of year
	double mwcontlower;
		// mean water content in lower soil layer for last month
		// (valid only on last day of month following call to daily_accounting_patch)

//////////////////////////////////////////////////////////////////////////////////
// GUESSN: CENTURY SOM pools and other variables

	Sompool sompool[NSOMPOOL];

	double dperc;				// daily percolation (mm)
	double dbaseflow;			// daily baseflow (mm)
	double wcontmm_yesterday;	// ...

	double nmin_daily[365];		// daily N mineralisation (kgN/m2)
	double nimmob_daily[365];	// daily N immobilisation (kgN/m2)
	double leachfrac_daily[365]; // fraction of excess mineral N leached each day;
	double nmass_avail;			// soil mineral N pool (kgN/m2)

	double nmin_annual;			// annual sum of N mineralisation
	double nimmob_annual;		// annual sum of N immobilisation
	double nleach_annual;		// annual leaching from available N pool
	double ndep_annual;			// annual N deposition

	double setntoc_nmass_avail;	// soil mineral N pool (kgN/m2) (used in daily setntoc)
	double daily_minimmndep;	// sum of mineralization, immobilization and N deposition (used in daily setntoc)

	double N_fix;				// total annual N fixation

	double nmass_avail_daily;	// soil mineral N pool (kgN/m2) (used when trying to do daily N uptake)
	double daily_leaching[365];	// daily N uptake leaching 

// end GUESSN


	// MEMBER FUNCTIONS

public:

	Soil(Patch& p,Soiltype& s):patch(p),soiltype(s) {}
		// constructor (initialises member variable patch)

	void initdrivers() {

		// Initialises certain member variables

		alag=0.0;
		exp_alag=1.0;
		cpool_slow=0.0;
		cpool_fast=0.0;
		decomp_litter_mean=0.0;
		k_soilfast_mean=0.0;
		k_soilslow_mean=0.0;
		wcont[0]=0.0;
		wcont[1]=0.0;
		wcont_evap=0.0;
		snowpack=0.0;
		last_gtemp=-1;
		last_mgtemp=-1;


		// guess2008 - extra initialisation
		mwcontupper = 0.0;
		mwcontlower = 0.0;
		for (int mth = 0; mth < 12; mth++) {
			mwcont[mth][0] = 0.0;
			mwcont[mth][1] = 0.0;
		}

		for (int d=0; d<365; d++) {
			dwcontupper[d] = 0.0;
			dwcontlower[d] = 0.0;

			// GUESSN
			nmin_daily[d]=0.0;	
			nimmob_daily[d]=0.0;	
			leachfrac_daily[d]=0.0;
			// end GUESSN
		}

		/////////////////////////////////////////////////////
		// GUESSN: Initialise CENTURY pools

		for (int p=0;p<NSOMPOOL;p++)
			sompool[p].init();

		// Set initial CENTURY pool N:C ratios 
		// Parton et al 1993, Fig 4

		sompool[SOILMICRO].ntoc=1.0/15.0;
		sompool[SLOWSOM].ntoc=1.0/20.0;
		sompool[PASSIVESOM].ntoc=1.0/10.0;
		sompool[SURFMICRO].ntoc=1.0/20.0;

		nmass_avail=0.0;

		nmin_annual=0.0;			
		nimmob_annual=0.0;		
		nleach_annual=0.0;		
		ndep_annual=0.0;
		N_fix=0.0;

		dperc=0.0;
		dbaseflow=0.0;
		wcontmm_yesterday=0.0;

		setntoc_nmass_avail=0.0;
		daily_minimmndep=0.0;	

		nmass_avail_daily=0.0;

		// end GUESSN

	}

};

///////////////////////////////////////////////////////////////////////////////////////
// CLASS LOOKUP_LAMBDA
// Lookup table for photosynthesis parameters (required for "fast" version of canopy
// exchange code; see canexch.cpp)

const int LOOKUP_LAMBDA_MAXITEM=130;

struct Lookup_lambda_item {
	double adtmm;
	double agd;
	double rd;
	double nmass_term; // GUESSN
	int year;
	int day;

	Lookup_lambda_item()
			: adtmm(0.0), agd(0.0), rd(0.0), nmass_term(0.0), year(-1), day(0) {
	}
};


class Lookup_lambda {

private:
	Lookup_lambda_item data[LOOKUP_LAMBDA_MAXITEM];
	int position;

public:

	void newsearch() {
		position=0;
	}

	bool getdata(int year,int day,double& adtmm,double& agd,double& rd,double& nmass_term) {
		if (position>=LOOKUP_LAMBDA_MAXITEM)
			fail("class Lookup_lambda: exceeded dimension of lookup table");
		Lookup_lambda_item& thisitem=data[position];
		if (thisitem.year==year && thisitem.day==day) {
			adtmm=thisitem.adtmm;
			agd=thisitem.agd;
			rd=thisitem.rd;
			nmass_term=thisitem.nmass_term;	// GUESSN
			return true;
		}
		// else
		return false;
	}

	void setdata(int year,int day,double adtmm,double agd,double rd, double nmass_term) {
		if (position>=LOOKUP_LAMBDA_MAXITEM)
			fail("class Lookup_lambda: exceeded dimension of lookup table");
		Lookup_lambda_item& thisitem=data[position];
		thisitem.year=year;
		thisitem.day=day;
		thisitem.adtmm=adtmm;
		thisitem.agd=agd;
		thisitem.rd=rd;
		thisitem.nmass_term=nmass_term;	// GUESSN
	}

	bool increase() {
		position+=position+1;
		if (position>=LOOKUP_LAMBDA_MAXITEM) return false;
		return true;
	}

	bool decrease() {
		position+=position+2;
		if (position>=LOOKUP_LAMBDA_MAXITEM) return false;
		return true;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// PATCHPFT
// State variables common to all individuals of a particular PFT in a particular patch
// Used in individual and cohort modes only.

class Patchpft {

	// MEMBER VARIABLES:

public:

	int id;
		// id code (equal to value of member variable id in corresponding Pft object)
	Pft& pft;
		// reference to corresponding Pft object in PFT list
	double anetps_ff;
		// potential annual net assimilation (leaf-level net photosynthesis) at forest
		// floor (kgC/m2/year)
	double wscal;
		// water stress parameter (0-1 range; 1=minimum stress)
	double wscal_mean;
		// running sum (converted to annual mean) for wscal
	double anetps_ff_est;
		// potential annual net assimilation at forest floor averaged over
		// establishment interval (kgC/m2/year)
	double anetps_ff_est_initial;
		// first-year value of anetps_ff_est
	double wscal_mean_est;
		// annual mean wscal averaged over establishment interval
	double phen;
		// vegetation phenological state (fraction of potential leaf cover)
		// updated daily
	double aphen;
		// annual sum of daily fractional leaf cover (equivalent number of days with
		// full leaf cover) (reset on expected coldest day of year)
	bool establish;
		// whether PFT can establish in this patch under current conditions
	double nsapling;
		// running total for number of saplings of this PFT to establish (cohort mode)
	double litter_leaf;
		// leaf-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_root;
		// fine root-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_wood;
		// heartwood and sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_repr;
		// litter derived from allocation to reproduction for PFT on modelled area
		// basis (kgC/m2)

	// Variables used by "fast" canopy exchange code (Ben Smith 2002-07)

	double gcbase;
		// non-FPC-weighted canopy conductance value for PFT under water-stress
		// conditions (mm/s)
	double temp_wstress;
		// cumulative mean temperature for water stress days this month (deg C)
	double par_wstress;
		// cumulative mean PAR for water stress days this month (J/m2/day)
	double daylength_wstress;
		// cumulative mean day length for water stress days this month (h)
	double co2_wstress;
		// cumulative mean atmospheric CO2 concentration for water stress days this
		// month (ppmv)
	int nday_wstress;
		// cumulative number of water stress days this month
	double fpar_grass_wstress;
		// mean FPAR at top of grass canopy for days with water stress for this PFT
		// in this patch
	double supply;
		// evapotranspirational "supply" function for this PFT today (mm/day)
	double supply_leafon;
	double fuptake[NSOILLAYER];
		// fractional uptake of water from each soil layer today
	bool ifwstress;
		// whether water-stress conditions for this PFT today
	Lookup_lambda lookup_lambda;
		// lookup table for values of lambda (parameter in photosynthesis calculations)
		// today (see canexch.cpp)

	double nsapling_yearly;
	int no_cohorts; // no cohorts used for averaging greffmort for monitoring
	double greff_mort_fraction;

	// GUESSN
	double nmass_litter_leaf;
	double nmass_litter_root;
	double nmass_litter_wood;

	double nstore_est;
		// N store for establishment
	double nsapling_nuptake;
		// number of saplings of this PFT established in vegetation_n_uptake() (cohort mode)
	// end GUESSN

	// MEMBER FUNCTIONS:

	Patchpft(int i,Pft& p):id(i),pft(p) {

		// Constructor: initialises id, pft and data members

		litter_leaf=0.0;
		litter_root=0.0;
		litter_wood=0.0;
		litter_repr=0.0;
		nday_wstress=0;
		wscal=1.0;
		wscal_mean=0.0;
		anetps_ff=0.0;
		aphen=0.0;

		// GUESSN
		nmass_litter_leaf=0.0;
		nmass_litter_root=0.0;
		nmass_litter_wood=0.0;

		nsapling_nuptake=0.0;

		nstore_est=0.0;
		
		// end GUESSN
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// PATCH
// Stores data for a patch. In cohort and individual modes, replicate patches are
// required in each stand to accomodate stochastic variation; in population mode there
// should be just one Patch object, representing average conditions for the entire
// stand. A reference to the parent Stand object (defined below) is included as a
// member variable.

class Patch {

public:

	// MEMBER VARIABLES

	int id;
		// id code in range 0-npatch for patch
	Stand& stand;
		// reference to parent Stand object
	ListArray_idin1<Patchpft,Pft> pft;
		// list array [0...npft-1] of Patchpft objects (initialised in constructor)
	Vegetation vegetation;
		// vegetation for this patch
	Soil soil;
		// soil for this patch
	Fluxes fluxes;
		// fluxes for this patch
	double fpar_grass;
		// FPAR at top of grass canopy today
	double fpar_ff;
		// FPAR at soil surface today
	double par_grass_mean;
		// mean growing season PAR at top of grass canopy (J/m2/day)
	int nday_growingseason;
		// number of days in growing season, estimated from mean vegetation leaf-on
		// fraction (see function fpar in canopy exchange module)
	double fpc_total;
		// total patch FPC
	bool disturbed;
		// whether patch was disturbed last year
	int age;
		// patch age (years since last disturbance)
	double fireprob;
		// probability of fire this year

	int growingseasondays;
		// guess2008 - DLE - the number of days over which wcont is averaged for this 
		// patch, i.e. those days for which daily temp > 5.0 degC


	// Variables used by new hydrology (Dieter Gerten 2002-07)

	double intercep;
		// interception by vegetation today on patch basis (mm)
	double aaet;
		// annual sum of AET (mm/year)
	double aevap;
		// annual sum of soil evaporation (mm/year)
	double aintercep;
		// annual sum of interception (mm/year)
	double arunoff;
		// annual sum of runoff (mm/year)
	double apet;
		// annual sum of potential evapotranspiration (mm/year)

	double eet_net_veg;
		// equilibrium evapotranspiration today, deducting interception (mm)
	double demand;
		// transpirative demand for patch today, mm/day, patch vegetative area basis
	double demand_leafon;
		// transpirative demand for patch assuming full leaf cover today, mm/day,
		// patch vegetative area basis
	double fpc_rescale;
		// rescaling factor to account for spatial overlap between individuals/cohorts
		// populations

	double maet[12];
		// monthly AET (mm/month)
	double mevap[12];
		// monthly soil evaporation (mm/month)
	double mintercep[12];
		// monthly interception (mm/month)
	double mrunoff[12];
		// monthly runoff (mm/month)
	double mpet[12];
		// monthly PET (mm/month)

	// GUESSN
	double fnuptake;
		// fractional N uptake of patch demand
	double ndemand;
		// yearly N demand
	double nsupply;
		// yearly N supply
	double new_est_ndemand;
		// last years N demand for new establishments
	double est_ndemand;
		// cumulative N demand for new establishments
	// end GUESSN

	// MEMBER FUNCTIONS

	Patch(int i,Stand& s,Pftlist& pftlist,Soiltype& st):
		id(i),stand(s),fluxes(*this),vegetation(*this),soil(*this,st) {
		
		// Constructor: initialises various members and builds list array
		// of Patchpft objects.

		pftlist.firstobj();
		while (pftlist.isobj) {
			pft.createobj(pftlist.getobj());
			pftlist.nextobj();
		}

		age=0;
		disturbed=false;
		
		// guess2008 - initialise
		growingseasondays=0;

		fireprob=0.0;
		est_ndemand=0.0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// STANDPFT
// State variables common to all individuals of a particular PFT in a stand. Used in
// individual and cohort modes only.

class Standpft {

public:

	// MEMBER VARIABLES

	int id;
	Pft& pft;
	double cmass_repr;
		// net C allocated to reproduction for this PFT in all patches of this stand
		// this year (kgC/m2)
	double anetps_ff_max;
		// maximum value of anetpsff (potential annual net assimilation at forest
		// floor) for this PFT in this stand so far in the simulation (kgC/m2/year)
	double addtw;
		// annual degree day sum above threshold damaging temperature (used in
		// calculation of heat stess mortality; Sitch et al 2000, Eqn 55)
	double gterm;
		// term in calculation of potential canopy conductance (mm/s)
	bool have_gterm;
		// true if value of gterm available for this PFT today, otherwise false

	// Variables used only by input/output module

	double cmass_total;
		// sum/mean across patches for carbon biomass (kgC/m2)
	double anpp_total;
		// sum/mean across patches for annual NPP (kgC/m2/year)
	double lai_total;
		// sum/mean across patches for 'grid-cell' LAI
	double densindiv_total;
		// sum/mean across patches for density of (true) individuals (indiv/m2)
		// (meaningful in cohort/individual mode only)
	double densindiv_ageclass[OUTPUT_MAXAGECLASS];
		// stem density by age class (cohort/individual mode only; used by function
		// outannual)
	double greff_mort_total;
		// sum/mean across patches for saplings per PFT
	double nsapling_total;
		// sum/mean across patches for saplings per PFT

	// Variables used by "fast" canopy exchange code (Ben Smith 2002-07)

	double gpterm;
		// non-FPAR-weighted value for canopy conductance component associated with
		// photosynthesis for PFT under non-water-stress conditions (mm/s)
	double assim_term;
		// non-FPAR-weighted leaf-level net photosynthesis value for PFT under non-
		// water-stress conditions (kgC/m2/day)
	bool have_phot;
		// whether gpterm and assim_term values are valid today
	double fpc_total;
		// FPC sum for this PFT as average for stand (used by some versions of
		// guessio.cpp)

	// GUESSN
	double na;
		// Leaf nitrogen associated with photosynthesis today, patch basis kgN/m2
	double cmass_repr_nuptake;
		// net C allocated to reproduction for this PFT in all patches of this stand
		// this year (kgC/m2)
	// end GUESSN

	// MEMBER FUNCTIONS

	Standpft(int i,Pft& p):id(i),pft(p) {
		
		// Constructor: initialises various data members
		
		anetps_ff_max=0.0;
		addtw=0.0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// STAND
// Stores data for a stand, corresponding to a modelled locality or grid cell.
// Member variables include an object of type Climate (holding climate, insolation and
// CO2 data), a object of type Soiltype (holding soil static parameters) and a list
// array of Patch objects. Soil objects (holding soil state variables) are associated
// with patches, not stands. A separate Stand object must be declared for each modelled
// locality or grid cell.

class Stand : public ListArray_idin3<Patch,Stand,Pftlist,Soiltype> {

public:

	// MEMBER VARIABLES

	ListArray_idin1<Standpft,Pft> pft;
		// list array [0...npft-1] of Standpft (initialised in constructor)
	Climate climate;
		// climate, insolation and CO2 for this stand
	Soiltype soiltype;
		// soil static parameters for this stand

	// MEMBER FUNCTIONS

	// FACE DAVID Thomas plantation general
	int plantyear; // year of plantation of FACE forest
	int distyear2; // cutting natural forest to be replaced by grassland: 1700 for Duke, Oak Ridge unclear, was also field
	// Thomas2 this year's ndemand: checking, see grwoth.cpp

	Stand(Pftlist& pftlist):climate(*this) {
		
		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
		int p;

		pftlist.firstobj();
		while (pftlist.isobj) {
			pft.createobj(pftlist.getobj());
			pftlist.nextobj();
		}

		for (p=0;p<npatch;p++) createobj(*this,pftlist,soiltype);
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK FUNCTION DECLARATION

int framework(int argc,char* argv[]);

#endif // LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Prentice, IC, Sykes, MT & Cramer W (1993) A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Reich, PB, Walters MB & Ellsworth DS 1997 From tropics to tundra: global
//   convergence in plant functioning. Proceedings of the National Academy of Sciences
//   USA 94: 13730-13734.
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Sykes, MT, Prentice IC & Cramer W 1996 A bioclimatic model for the potential
//   distributions of north European tree species under present and future climates.
//   Journal of Biogeography 23: 209-233.
