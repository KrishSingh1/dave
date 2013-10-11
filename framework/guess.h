///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.h
/// \brief Framework header file, LPJ-GUESS Combined Modular Framework
///
/// This header file contains:
///  (1) definitions of all main classes used by the framework and modules. Modules may
///      require classes to contain certain member variables and functions (see module
///      source files for details).
///  (2) other type, constant and function definitions to be accessible throughout the
///      model code.
///  (3) a forward declaration of the framework function if this is not the main
///      function.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESS_H
#define LPJ_GUESS_GUESS_H

// Defines for landcover version:
#define DYNAMIC_LANDCOVER_INPUT		// Reads landcover data from text files, using the TimeDataD class.
#define LUTOMEMORY					// Write land use fraction data to memory; enables efficient usage of randomized gridlists for parallell runs on Simba.
#define NEWSOWINGDATE				// Use sowing date method based on climate seasonality (modified version of Waha et al. 2012), as opposed to old method used in Bondeau et al. 2007.
#define IRRIGATION					// Crop irrigation on
#define NOPASTURESTOCH				// Undefine for fire and disturbance for pasture grass. Number of patches will be the same as for natural stands.
//#define GRASSFORCROP				// Transfer cropland to pasture landcover for simplified crop definition (harvested competing c3/c4 grass).

const bool SUPPRESSLARGEOUTPUT=true;

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gutil.h>
#include <vector>
#include "shell.h"
#include "guessmath.h"
#include "archive.h"
#include "parameters.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

/// Life form class for PFTs (trees, grasses)
typedef enum {NOLIFEFORM, TREE, GRASS} lifeformtype;

/// Phenology class for PFTs
typedef enum {NOPHENOLOGY,EVERGREEN,RAINGREEN,SUMMERGREEN,CROPGREEN,ANY} phenologytype;

/// Biochemical pathway for photosynthesis (C3 or C4)
typedef enum {NOPATHWAY, C3, C4} pathwaytype;

/// Leaf physiognomy types for PFTs
typedef enum {NOLEAFTYPE, NEEDLELEAF, BROADLEAF} leafphysiognomytype;

/// Units for insolation driving data
/** Insolation can be expressed as:
 *
 *  - Percentage sunshine
 *  - Net instantaneous downward shortwave radiation flux (W/m2)
 *  - Total (i.e. with no correction for surface albedo) instantaneous downward
 *    shortwave radiation flux (W/m2)
 *
 *  Radiation flux can be interpreted as W/m2 during daylight hours, or averaged
 *  over the whole time step which it represents (24 hours in daily mode). For
 *  this reason there are two enumerators for these insolation types (e.g. SWRAD
 *  and SWRAD_TS).
 */
typedef enum {
	/// No insolation type chosen
	NOINSOL,
	/// Percentage sunshine
	SUNSHINE,
	/// Net shortwave radiation flux during daylight hours (W/m2)
	NETSWRAD,
	/// Total shortwave radiation flux during daylight hours (W/m2)
	SWRAD,
	/// Net shortwave radiation flux during whole time step (W/m2)
	NETSWRAD_TS,
	/// Total shortwave radiation flux during whole time step (W/m2)
	SWRAD_TS
} insoltype;

/// CENTURY pool names, NSOMPOOL number of SOM pools
typedef enum {SURFSTRUCT, SOILSTRUCT, SOILMICRO, SURFHUMUS, SURFMICRO, SURFMETA, SURFFWD, SURFCWD,
	SOILMETA, SLOWSOM, PASSIVESOM, LEACHED, NSOMPOOL} pooltype;	

/// Irrigation type for PFTs
typedef enum {RAINFED, IRRIGATED} hydrologytype;
/// Intercrop type for PFTs
typedef enum {NOINTERCROP, NATURALGRASS} intercroptype;

/// Seasonality type of gridcell
/** 0:SEASONALITY_NO			No seasonality
 *  1:SEASONALITY_PREC			Precipitation seasonality only
 *  2:SEASONALITY_PRECTEMP		Both temperature and precipitation seasonality, but "weak" temperature seasonality (coldest month > 10degC)
 *  3:SEASONALITY_TEMP			Temperature seasonality only
 *  4:SEASONALITY_TEMPPREC		Both temperature and precipitation seasonality, but temperature most important (coldest month < 10degC)
 *  5:SEASONALITY_TEMPWARM		Temperature seasonality, always above 10 degrees (currently not used)
 */
typedef enum {SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC} seasonality_type;

/// Precipitation seasonality type of gridcell
/** 0:DRY						(minprec_pet20<=0.5 && maxprec_pet20<=0.5)
 *  1:DRY_INTERMEDIATE			(minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
 *  2:DRY_WET					(minprec_pet20<=0.5 && maxprec_pet20>1.0)
 *  3:INTERMEDIATE				(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
 *  4:INTERMEDIATE_WET			(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)
 *  5:WET						(minprec_pet20>1.0 && maxprec_pet20>1.0)
 */
typedef enum {DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET} prec_seasonality_type;


/// Temperature seasonality type of gridcell
/** 0:COLD						(mtemp_max20<=10)
 *  1:COLD_WARM					(mtemp_min20<=10 && mtemp_max20>10 && mtemp_max20<=30)
 *  2:COLD_HOT					(mtemp_min20<=10 && mtemp_max20>30)
 *  3:WARM						(mtemp_min20>10 && mtemp_max20<=30)
 *  4:WARM_HOT					(mtemp_min20>10 && mtemp_max20>30)
 *  5:HOT						(mtemp_min20>30)
 */
typedef enum {COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT} temp_seasonality_type;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

/// number  of soil layers modelled
const int NSOILLAYER = 2;

// SOIL DEPTH VALUES

/// soil upper layer depth (mm)
const double SOILDEPTH_UPPER = 500.0;
/// soil lower layer depth (mm)
const double SOILDEPTH_LOWER = 1000.0;

/// Year at which to calculate equilibrium soil carbon
const int SOLVESOM_END=400;

/// Year at which to begin documenting means for calculation of equilibrium soil carbon
const int SOLVESOM_BEGIN = 350;

/// Number of years to average growth efficiency over in function mortality
const int NYEARGREFF = 5;

/// Coldest day in N hemisphere (January 15)
/** Used to decide when to start counting GDD's and leaf-on days 
 *  for summergreen phenology.
 */
const int COLDEST_DAY_NHEMISPHERE = 14;

/// Coldest day in S hemisphere (July 15)
/** Used to decide when to start counting GDD's and leaf-on days 
 *  for summergreen phenology.
 */
const int COLDEST_DAY_SHEMISPHERE = 195;

/// number of years to average aaet over in function soilnadd
const int NYEARAAET = 5;

/// Maximum number of age classes in age structure plots produced by function outannual
const int OUTPUT_MAXAGECLASS = 40;

/// Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
const double PRIESTLEY_TAYLOR = 1.32;

// Solving Century SOM pools 

/// fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
const double SOLVESOMCENT_SPINBEGIN  = 0.1;
/// fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
const double SOLVESOMCENT_SPINEND    = 0.3;

/// Kelvin to deg c conversion
const double K2degC = 273.15;
const int NCROPSTANDS_MAX=26;

/// Conversion factor for CO2 from ppmv to mole fraction
const double CO2_CONV = 1.0e-6;


///////////////////////////////////////////////////////////////////////////////////////
// FORWARD DECLARATIONS OF CLASSES DEFINED IN THIS FILE
// Forward declarations of classes used as types (e.g. for reference variables in some
// classes) before they are actually defined

class Date;
class Stand;
class Patch;
class Vegetation;
class Gridcell;
class Patchpft;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are defined in the framework source code file, and are accessible
// throughout the code

/// Object describing timing stage of simulation
extern Date date;

/// Number of possible PFTs
extern int npft;
extern int ncft;


/// General purpose object for handling simulation timing.
/** In general, frameworks should use a single Date object for all simulation
 *  timing.
 *
 *  Member variables of the class (see below) provide various kinds of calender
 *  and timing information, assuming init has been called to initialise the
 *  object, and next() has been called at the end of each simulation day.
 */
class Date {

	// MEMBER VARIABLES

public:

	/// number of days in each month (0=January - 11=December)
	int ndaymonth[12];

	/// julian day of year (0-364; 0=Jan 1)
	int day;

	/// day of current month (0=first day)
	int dayofmonth;

	/// month number (0=January - 11=December)
	int month;

	/// year since start of simulation (0=first simulation year)
	int year;

	/// number of subdaily periods in a day (to be set in IO module)
	int subdaily;

	/// julian day for middle day of each month
	int middaymonth[12];

	/// true if last year of simulation, false otherwise
	bool islastyear;

	/// true if last month of year, false otherwise
	bool islastmonth;

	/// true if last day of month, false otherwise
	bool islastday;

	/// true if middle day of month, false otherwise
	bool ismidday;


private:

	int nyear;

	// MEMBER FUNCTIONS

public:

	/// Constructor function called automatically when Date object is created
	/** Do not call explicitly. Initialises some member variables. */
	Date() {
		const int data[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
		int month;
		int dayct = 0;
		for (month=0; month<12; month++) {
			ndaymonth[month] = data[month];
			middaymonth[month] = dayct + data[month] / 2;
			dayct += data[month];
		}
		subdaily = 1;
	}

	/// Initialises date to day 0 of year 0 and sets intended number of simulation years
	/** Intended number of simulation years is only used to set islastyear flag,
	 *  actual simulation may be longer or shorter.
	 *
	 *  \param nyearsim  Intended number of simulation years
	 */
	void init(int nyearsim)	{
		nyear = nyearsim;
		day = month=year = dayofmonth = 0;
		islastmonth = islastday = ismidday = false;
		if (nyear == 1) islastyear = true;
		else islastyear = false;
	}

	/// Call at end of every simulation day to update member variables.
	void next() {
		if (islastday) {
			if (islastmonth) {
				dayofmonth = 0;
				day = 0;
				month = 0;
				year++;
				if (year == nyear - 1) islastyear = true;
				islastmonth = false;
			}
			else {
				day++;
				dayofmonth = 0;
				month++;
				if (month == 11) islastmonth = true;
			}
			islastday = false;
		}
		else {
			day++;
			dayofmonth++;
			if (dayofmonth == ndaymonth[month] / 2) ismidday = true;
			else {
				ismidday = false;
				if (dayofmonth == ndaymonth[month] - 1) islastday = true;
			}
		}
	}

	// \returns index (0-11) of previous month (11 if currently month 0).
	int prevmonth() {
		if (month > 0) return month - 1;
		return 11;
	}

	/// \returns index of next month (0 if currently month 11)
	int nextmonth() {
		if (month < 11) return month+1;
		return 0;
	}

	/// Check if the year is leap
	/** \param year        Calendar year
	*   The algorith is as follows: only year that are divisible by 4 could
	*   potentially be leap (e.g., 1904), however, not if they're divisble by
	*   100 (e.g., 1900 is not leap), unless they're divisble by 400 (e.g., 2000
	*   is still leap).
	*/
	static bool is_leap(int year) {
		return (!(year % 4) && (year % 100 | !(year % 400)));
	}
	
	/// Whether the current mode is diurnal
	bool diurnal() const { return subdaily > 1; }
};

/// Object describing sub-daily periods
class Day {
public:
	/// Whether sub-daily period first/last within day (both true in daily mode)
	bool isstart, isend;
	
	/// Ordinal number of the sub-daily period [0, date.subdaily)
	int period;

	/// Constructs beginning of the day period (the only one in daily mode)
	Day() {
		isstart = true;
		isend = !date.diurnal();
		period = 0;
	}

	/// Advances to the next sub-daily period
	void next() {
		period++;
		isstart = false;
		isend = period == date.subdaily - 1;
	}
};

/// This struct contains the result of a photosynthesis calculation.
/** \see photosynthesis */  
struct PhotosynthesisResult : public Serializable {
	/// Constructs an empty result
	PhotosynthesisResult() {
		clear();
	}

	/// Clears all members
	/** This is returned by the photosynthesis function when no photosynthesis
	 *  takes place.
	 */
	void clear() {
		agd_g       = 0;
		adtmm       = 0;
		rd_g        = 0;
		vm          = 0;
		je          = 0;
		nactive_opt = 0.0;
		vmaxnlim    = 1.0;
	}

	/// RuBisCO capacity (gC/m2/day)
	double vm;

	/// gross daily photosynthesis (gC/m2/day)
	double agd_g;

	/// leaf-level net daytime photosynthesis
	/** expressed in CO2 diffusion units (mm/m2/day) */
    double adtmm;

	/// leaf respiration (gC/m2/day)
	double rd_g;

	/// PAR-limited photosynthesis rate (gC/m2/h)
    double je;

	/// optimal leaf nitrogen associated with photosynthesis (kgN/m2)
	double nactive_opt;

	/// nitrogen limitation on vm
	double vmaxnlim;

	/// net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
    double net_assimilation() const {
		return (agd_g - rd_g) * 1e-3;
    }

	void serialize(ArchiveStream& arch);
};


/// The Climate for a grid cell
/** Stores all static and variable data relating to climate parameters, as well as
 *  latitude, atmospheric CO2 concentration and daylength for a grid cell. Includes
 *  a reference to the parent Gridcell object (defined below). Initialised by a
 *  call to initdrivers.
 */
class Climate : public Serializable {

	// MEMBER VARIABLES

public:
	Gridcell& gridcell;
		// reference to parent Gridcell object
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
	double lon; 
		// longitude 

	/// Insolation today, see also instype
	double insol;

	/// Type of insolation
	/** This decides how to interpret the variable insol,
	 *  see also documentation for the insoltype enum.
	 */
	insoltype instype;

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
	double gdd5_pasture;	// gdd5 for pasture grass (not sensitive to ifsensechill)
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
	Historic<double, 31> dtemp_31;
		// daily temperatures for the last 31 days (deg C)
	double mtemp_min_20[20];
		// minimum monthly temperatures for the last 20 years (deg C)
	double mtemp_max_20[20];
	double mtemp_min;
		// minimum monthly temperature for the last 12 months (deg C)
	double atemp_mean;
		// mean of monthly temperatures for the last 12 months (deg C)

	/// annual nitrogen deposition (kgN/m2/year)
	double andep;
	/// daily nitrogen deposition (kgN/m2)
	double dndep;

	/// annual nitrogen fertilization (kgN/m2/year)
	double anfert;
	/// daily nitrogen fertilization (kgN/m2/year)
	double dnfert;

	// Monthly sums (converted to means) used by canopy exchange module

	double temp_mean;
		// accumulated mean temperature for this month (deg C)
	double par_mean;
		// accumulated mean daily net PAR sum (J/m2/day) for this month
	double co2_mean;
		// accumulated mean CO2 for this month (ppmv)
	double daylength_mean;
		// accumulated mean daylength for this month (h)

	// Saved parameters used by function daylengthinsoleet

	double sinelat;
	double cosinelat;
	double qo[365], u[365], v[365], hh[365], sinehh[365];
	double daylength_save[365];
	bool doneday[365];
		// indicates whether saved values exist for this day

	double dtr; // diurnal temperature range, used in daily/monthly BVOC (deg C)

	std::vector<double> temps, insols, pars, rads, gtemps;
		// containers for sub-daily values of temperature, short-wave downward
		// radiation, par, rad and gtemp (equivalent to temp, insol, par, rad and gtemp)
		// NB: units of these variable are the same as their daily counterparts,
		// i.e. representing daily averages (e.g. pars [J/m2/day])

	/// Variables used for crop sowing date or seasonality calculation

	/// daily precipitations for the last 10 days (mm)
	double dprec_10[10];
	/// daily 10 day-sums of precipitations for today and yesterday (mm)
	double sprec_2[2];
	/// max temperature during the last test period
	double maxtemp;
	/// summer day when we test last year's crossing of sowing temperature limits; NH:June 30(day 180), SH:Dec.31(day 364), set in getgridcell()
	int testday_temp;
	/// last day of dry month when we test last year's crossing of sowing precipitation limits; NH:Dec.31(day 364), SH:June 30(day 180), set in getgridcell()
	int testday_prec;
	/// date used for sowing if no frost or spring occured during the year between the testmonths; NH:14, SH:195, set in getgridcell()
	int coldestday;
	/// used to adapt equations to hemisphere, set in getgridcell()
	int adjustlat;

	/// copy of monthly historical or spinup temperature values for one year (from local guessio.cpp-variables hist_mtemp or spinup_mtemp).
	double mtemp_year[12];
	/// copy of monthly historical or spinup precipitation values for one year (from local guessio.cpp-variables hist_mtemp or spinup_mtemp).
	double mprec_year[12];
	/// accumulated monthly pet values for this year
	double mpet_year[12];

	/// past 20 years monthly temperature values
	double mtemp_20[20][12];
	/// past 20 years monthly precipitation values
	double mprec_20[20][12];
	/// past 20 years monthly PET values
	double mpet_20[20][12];
	/// past 20 years monthly precipitation to PET ratios
	double mprec_pet_20[20][12];
	/// past 20 years minimum of monthly precipitation to PET ratios
	double mprec_petmin_20[20];
	/// past 20 years maximum of monthly precipitation to PET ratios
	double mprec_petmax_20[20];

	/// 20-year running average monthly temperature values 
	double mtemp20[12];
	/// 20-year running average monthly precipitation values 
	double mprec20[12];
	/// 20-year running average monthly PET values 
	double mpet20[12];
	/// 20-year running average monthly precipitation to PET ratios 
	double mprec_pet20[12];

	/// 20-year running average of minimum monthly precipitation to PET ratios
	double mprec_petmin20;
	/// 20-year running average of maximum monthly precipitation to PET ratios
	double mprec_petmax20;

	/// seasonality type (SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC)
	seasonality_type seasonality;

	/// precipitation seasonality type (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the extremes of the 20-year monthly means
	 */
	prec_seasonality_type prec_seasonality;

	/// precipitation range (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the average of the 20-year monthly extremes
	 */
	prec_seasonality_type prec_range;

	/// temperature seasonality (COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT)
	temp_seasonality_type temp_seasonality;

	/// whether several months with precipitation maxima exists (remains to be implemented)
	bool biseasonal;

	/// variation coefficient of 20-year mean monthly temperatures
	double var_prec;
	/// variation coefficient of 20-year mean monthly precipitation to PET ratios
	double var_temp;

	/// annual precipitation sum
	double aprec;

	// MEMBER FUNCTIONS

public:
	Climate(Gridcell& gc):gridcell(gc) {
		memset(mtemp_20, 0, sizeof(double)*20*12);
		memset(mtemp20, 0, sizeof(double)*12);
		memset(mtemp_year, 0, sizeof(double)*12);

		memset(mprec_20, 0, sizeof(double)*20*12);
		memset(mprec20, 0, sizeof(double)*12);
		memset(mprec_year, 0, sizeof(double)*12);

		memset(mpet_20, 0, sizeof(double)*20*12);
		memset(mpet20, 0, sizeof(double)*12);
		memset(mpet_year, 0, sizeof(double)*12);

		memset(mprec_pet_20, 0, sizeof(double)*20*12);
		memset(mprec_pet20, 0, sizeof(double)*12);

		memset(mprec_petmin_20, 0, sizeof(double)*20);		//111005
		mprec_petmin20=0.0;		//111005
		memset(mprec_petmax_20, 0, sizeof(double)*20);		//111115
		mprec_petmax20=0.0;		//111115

		seasonality=SEASONALITY_NO;
		prec_seasonality=DRY;
		prec_range=DRY;
		temp_seasonality=COLD;
		biseasonal=false;

		eet=0.0;
	};
		// constructor function: initialises gridcell member

	void initdrivers(double latitude) {

		// Initialises certain member variables
		// Should be called before Climate object is applied to a new grid cell

		int day, year;

		for (year=0; year<20; year++) {
			mtemp_min_20[year] = 0.0;
			mtemp_max_20[year] = 0.0;
		}
		mtemp_min20 = 0.0;
		mtemp_max20 = 0.0;
		mtemp = 0.0;
		maxtemp = 0.0;
		gdd5 = 0.0;
		gdd5_pasture = 0.0;
		chilldays = 0;
		ifsensechill = true; //  guess2008 - CHILLDAYS
		atemp_mean = 0.0;

		lat = latitude;
		for (day=0; day<365; day++) doneday[day] = false;
		sinelat = sin(lat * DEGTORAD);
		cosinelat = cos(lat * DEGTORAD);
	}

	void serialize(ArchiveStream& arch);
};


/// Stores accumulated monthly and annual fluxes.
/** This class handles the storage and accounting of fluxes for a single patch.
 *  Different fluxes can be stored in different ways, depending on what kind of
 *  flux it is and what kind of output we want. The details of whether fluxes
 *  are stored per PFT or just as a patch total, or per day, month or only a
 *  yearly sum, is hidden from the 'scientific' code, which merely reports the
 *  fluxes generated.
 */
class Fluxes : public Serializable {

public:

	/// Fluxes stored as totals for the whole patch
	enum PerPatchFluxType { 
		/// Carbon flux to atmosphere from burnt vegetation and litter (kgC/m2)
		FIREC,
		/// Carbon flux to atmosphere from soil respiration (kgC/m2)
		SOILC,
		/// Flux from atmosphere to vegetation associated with establishment (kgC/m2)
		ESTC,
		/// Flux to atmosphere from consumed harvested products (kgC/m2)
		HARVESTC,
		/// Flux from atmosphere to vegetation associated with sowing (kgC/m2)
		SEEDC,
		/// Flux to atmosphere from consumed harvested products (kgN/m2)
		HARVESTN,
		/// NH3 flux to atmosphere from fire
		NH3_FIRE,
		/// NO flux to atmosphere from fire	
		NO_FIRE,
		/// NO2 flux to atmosphere from fire
		NO2_FIRE,
		/// N2O flux to atmosphere from fire	
		N2O_FIRE,
		/// N flux from soil
		N_SOIL,
		/// Number of types, must be last
		NPERPATCHFLUXTYPES
	};

	/// Fluxes stored per pft
	enum PerPFTFluxType { 
		/// NPP (kgC/m2)
		NPP,
		/// GPP (kgC/m2)
		GPP,
		/// Autotrophic respiration (kgC/m2)
		RA,
		/// Isoprene (mgC/m2)
		ISO,
		/// Monoterpene (mgC/m2)
		MON,
		/// Number of types, must be last
		NPERPFTFLUXTYPES
	};

	// emission ratios from fire (NH3, NO, NO2, N2O) Delmas et al. 1995
	// values in .cpp file

	static const double NH3_FIRERATIO;
	static const double NO_FIRERATIO;
	static const double NO2_FIRERATIO;
	static const double N2O_FIRERATIO;	


	/// Reference to patch to which this Fluxes object belongs
	Patch& patch;

	// MEMBER FUNCTIONS

public:
	/// constructor: initialises members
	Fluxes(Patch& p);

	/// Sets all fluxes to zero (call at the beginning of each year)
	void reset();

	void serialize(ArchiveStream& arch);

	/// Report flux for a certain flux type
	void report_flux(PerPFTFluxType flux_type, int pft_id, double value);

	/// Report flux for a certain flux type
	void report_flux(PerPatchFluxType flux_type, double value);

	/// \returns flux for a given month and flux type (for all PFTs)
	double get_monthly_flux(PerPFTFluxType flux_type, int month) const;

	/// \returns flux for a given month and flux type
	double get_monthly_flux(PerPatchFluxType flux_type, int month) const;

	/// \returns annual flux for a given PFT and flux type
	double get_annual_flux(PerPFTFluxType flux_type, int pft_id) const;

	/// \returns annual flux for a given flux type (for all PFTs)
	double get_annual_flux(PerPFTFluxType flux_type) const;

	/// \returns annual flux for a given flux type
	double get_annual_flux(PerPatchFluxType flux_type) const;

private:

	/// Stores one flux value per PFT and flux type
	std::vector<std::vector<double> > annual_fluxes_per_pft;

	/// Stores one flux value per month and flux type
	/** For the fluxes only stored as totals for the whole patch */
	double monthly_fluxes_patch[12][NPERPATCHFLUXTYPES];

	/// Stores one flux value per month and flux type
	/** For the fluxes stored per pft for annual values */
	double monthly_fluxes_pft[12][NPERPFTFLUXTYPES];
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
		// leaf phenology (NOPHENOLOGY,EVERGREEN,RAINGREEN,SUMMERGREEN,CROPGREEN, ANY)
	leafphysiognomytype leafphysiognomy;
		// leaf physiognomy (needleleaf, broadleaf)
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

	/// minimum leaf C:N mass ratio allowed when nitrogen demand is determined
	double cton_leaf_min;
	/// maximum leaf C:N mass ratio	allowed when nitrogen demand is determined
	double cton_leaf_max;
	/// average leaf C:N mass ratio (between min and max)
	double cton_leaf_avr;
	/// average fine root C:N mass ratio (connected cton_leaf_avr)
	double cton_root_avr;
	/// maximum fine root C:N mass ratio (used when mass is negligible) 	
	double cton_root_max;
	/// average sapwood C:N mass ratio (connected cton_leaf_avr)
	double cton_sap_avr;
	/// maximum sapwood C:N mass ratio (used when mass is negligible) 
	double cton_sap_max;
	/// reference fine root C:N mass ratio	
	double cton_root;
	/// reference sapwood C:N mass ratio	
	double cton_sap;
	/// Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
	double nuptoroot;
	/// coefficient to compensate for vertical distribution of fine root on nitrogen uptake
	double nupscoeff;
	/// fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
	double fnstorage;

	/// Michaelis-Menten kinetic parameters 
	/** Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
	double km_volume;

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

	// bvoc
	double ga;
	        // aerodynamic conductance (m s-1)
	double eps_iso;
 	        // isoprene emission capacity (ug C g-1 h-1)
	bool seas_iso;
	        // whether (1) or not (1) isoprene emissions show a seasonality
	double eps_mon;
	        // monoterpene emission capacity (ug C g-1 h-1)
	double storfrac_mon;
	        // fraction of monoterpene production that goes into storage pool (-)



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

	/// specifies type of landcover
	/** \see landcovertype */
	landcovertype landcover;	// specifies type of landcover (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND); initialized in constructor
	/// hydrology (RAINFED,IRRIGATED) 
	hydrologytype hydrology;
	/// intercrop (NOINTERCROP,NATURALGRASS)
	intercroptype intercrop;

	/// fraction of residue outtake at harvest
	double res_outtake;
	/// harvest efficiency
	double harv_eff;
	/// harvest efficiency of intercrop grass
	double harv_eff_ic;
	/// fraction of harvested products that goes into patchpft.harvested_products_slow
	double harvest_slow_frac;
	/// yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
	double turnover_harv_prod;
	/// crop pft id defined in the insfile (numbered 0-25, corresponding to columns 0-25 in crop cover input file)
	int cftid;
	/// whether pft may grow as cover crop 
	bool isintercropgrass;
	/// whether sowing date is calculated
	bool ifsdcalc;
	/// whether temperature dependent sowing date is calculated
	bool ifsdtemp;
	/// whether autumn temperature dependent sowing date is calculated
	bool ifsdautumn;
	/// whether spring temperature dependent sowing date is calculated
	bool ifsdspring;
	/// whether precipitation dependent sowing date is calculated
	bool ifsdprec;
	/// upper temperature limit for autumn sowing
	double tempautumn;
	/// lower temperature limit for spring sowing
	double tempspring;
	/// default sowing date in the northern hemisphere (julian day)
	int sdatenh;
	/// default sowing date in the southern hemisphere
	int sdatesh;
	/// latest date for harvesting in the northern hemisphere
	int hlimitdatenh;
	/// latest date for harvesting in the southern hemisphere
	int hlimitdatesh;
	/// default base temperature (°C) for heat unit (hu) calculation
	double tb;
	/// early sowing limit for precipitation-limited crops (northern hemisphere); only used in Crop_sowing_date_prec()
	int firstsowdatenh_prec;
	/// early sowing limit for precipitation-limited crops (southern hemisphere); only used in Crop_sowing_date_prec()
	int firstsowdatesh_prec;
	/// temperature under which vernalisation is possible (°C)
	double trg;
	/// default number of vernalising days required
	int pvd;
    /// sensitivity to the photoperiod effect [0-1]
	double psens;
	/// basal photoperiod (h) (pb<ps for longer days plants)
	double pb;
	/// saturating photoperiod (h) (ps<pb for shorter days plants)
	double ps;
	/// default potential heat units required for crop maturity (degree-days)
	double phu;
	/// fraction of growing season (phu) at which senescence starts [0-1]
	double fphusen;
	/// type of senescence curve (see Bondeau et al. 2007)
	bool shapesenescencenorm;
	/// fraction of maximal LAI still present at harvest [0-1]
	double flaimaxharvest;
	/// default maximum LAI (only used for intercrop grass in the case where no pasture is present in any stand)
	double laimax; 
	/// whether harvestable organs are above ground
	bool aboveground_ho;
	/// optimum harvest index
	double hiopt;
	/// minimum harvest index
	double himin;
	/// initial fraction of growing season's npp allocated to roots
	double frootstart;
	/// final fraction of growing season's npp allocated to roots
	double frootend;
	/// whether sowing dates are read from input file
	bool forcesowingdate;
	/// whether harvest dates are read from input file
	bool forceharvestdate;

	// MEMBER FUNCTIONS

public:

	Pft() {

		// Constructor (initialises array gdd0)

		int y;
		for (y=0; y<366; y++)
			gdd0[y] = -1.0; // value<0 signifies "unknown"; see function phenology()

		// guess2008 - DLE
		drought_tolerance = 0.0; // Default, means that the PFT will never be limited by drought.

		res_outtake = 0.0;
		harv_eff = 0.0;
		harv_eff_ic = 0.0;
		turnover_harv_prod = 1.0;	// default 1 year turnover time

		cftid=-1;
		isintercropgrass=false;
		ifsdcalc=false;
		ifsdtemp=false;
		ifsdautumn=false;
		ifsdspring=false;
		ifsdprec=false;
		sdatenh=-1;
		sdatesh=-1;
		hlimitdatenh=-1;
		hlimitdatesh=-1;
		tb=-999.9;
		trg=-999.9;
		pvd=-1;
		psens=-1.0;
		pb=-1.0;
		ps=-1.0;
		phu=-1.0;
		fphusen=-1.0;
		shapesenescencenorm=0;
		flaimaxharvest=-1.0;
		laimax=0.0;
		aboveground_ho=true;
		frootstart=0.0;
		frootend=0.0;
		forcesowingdate=false;
		forceharvestdate=false;
	}

	/// Calculates SLA given leaf longevity
	void initsla() {

		// Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
		// m2/kgC)

		if (leafphysiognomy == BROADLEAF) {
			sla = 0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * leaflong));
		}
		else if (leafphysiognomy == NEEDLELEAF) {
			sla = 0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * leaflong));
		}
	}

	void init_cton_min() {

		// Calculates minimum leaf C:N ratio given leaf longevity
		// Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
		// kgN/kgC)

		if (leafphysiognomy == BROADLEAF)
			cton_leaf_min = 500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * leaflong));
		else if (leafphysiognomy == NEEDLELEAF)
			cton_leaf_min = 500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * leaflong));
	}

	void init_cton_limits() {

		// Fraction between min and max C:N ratio White et al. 2000
		double frac_mintomax = 2.78;

		// Fraction between leaf and root C:N ratio
		double frac_leaftoroot = 1.16; // Friend et al. 1997
		
		// Fraction between leaf and sap wood C:N ratio
		double frac_leaftosap = 6.9;   // Friend et al. 1997

		// Max leaf C:N ratio
		cton_leaf_max = cton_leaf_min * frac_mintomax;
		
		// Average leaf C:N ratio
		cton_leaf_avr = 1.0 / ((1.0 / cton_leaf_min + 1.0 / cton_leaf_max) / 2.0);

		// Average fine root C:N ratio
		cton_root_avr = cton_leaf_avr * frac_leaftoroot;

		// Maximum fine root C:N ratio
		cton_root_max = cton_leaf_min * frac_leaftoroot * frac_mintomax;

		// Average sap C:N ratio
		cton_sap_avr  = cton_leaf_avr * frac_leaftosap;

		// Maximum sap C:N ratio
		cton_sap_max  = cton_leaf_min * frac_leaftosap * frac_mintomax;

		if (lifeform == GRASS)
			respcoeff /= 2.0 * cton_root / (cton_root_avr + cton_leaf_min * frac_leaftoroot);
		else
			respcoeff /= cton_root / (cton_root_avr + cton_leaf_min * frac_leaftoroot) +
			             cton_sap  / (cton_sap_avr  + cton_leaf_min * frac_leaftosap);
	}

	void init_nupscoeff() {

		// Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
		
		// Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
		// Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving 
		// an approximate advantage of 2 of having more roots in the upper soil layer
		const double upper_adv = 2.0;

		nupscoeff = rootdist[0] * upper_adv + rootdist[1];

	}

	void initregen() {

		// Initialises sapling/regen characteristics in population mode
		// following LPJF formulation; see function allometry in growth module.
		// Note: primary PFT parameters, including SLA, must be set before this
		//       function is called

		const double PI = 3.14159265;
		const double REGENLAI_TREE = 1.5;
		const double REGENLAI_GRASS = 0.001;
		const double SAPLINGHW = 0.2;

		if (lifeform == TREE) {

			// Tree sapling characteristics

			regen.cmass_leaf = pow(REGENLAI_TREE * k_allom1 * pow(1.0 + SAPLINGHW, k_rp) *
				pow(4.0 * sla / PI / k_latosa, k_rp * 0.5) / sla, 2.0 / (2.0 - k_rp));

			regen.cmass_sap = wooddens * k_allom2 * pow((1.0 + SAPLINGHW) *
				sqrt(4.0 * regen.cmass_leaf * sla / PI / k_latosa), k_allom3) *
				regen.cmass_leaf * sla / k_latosa;

			regen.cmass_heart = SAPLINGHW * regen.cmass_sap;
		}
		else if (lifeform == GRASS) {

			// Grass regeneration characteristics

			regen.cmass_leaf = REGENLAI_GRASS / sla;
		}

		regen.cmass_root = 1.0 / ltor_max * regen.cmass_leaf;
	}
};


/// A list of PFTs
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Pft objects. In general, there should be a single Pftlist object containing
 *  a single list of PFTs. Pft objects within the list are then referenced by the pft
 *  member of each Individual object.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential Pft objects can be referenced as array elements by id:
 *
 *   Pftlist pftlist;
 *   ...
 *   for (i=0; i<npft; i++) {
 *     Pft& thispft=pftlist[i];
 *     // query or modify object thispft here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   pftlist.firstobj();
 *   while (pftlist.isobj) {
 *     Pft& thispft=pftlist.getobj();
 *     // query or modify object thispft here
 *     pftlist.nextobj();
 *   }
 */
class Pftlist : public ListArray_id<Pft> {};

/// The one and only linked list of Pft objects	
extern Pftlist pftlist;

///////////////////////////////////////////////////////////////////////////////////////
//cropindiv_struct
//container for crop-specific data at the individual level

class cropindiv_struct : public Serializable {

public:

	//Plant carbon biomass variables are all on patch area basis (kgC/m2)

	/// year's harvestable organ C biomass (= ycmass_plant)
	double cmass_ho;
	/// above-ground pool C biomass (when calculating daily cmass_leaf from lai_crop) (= ycmass_agpool)
	double cmass_agpool;		
	/// nitrogen content of harvestable organs
	double nmass_ho;
	/// nitrogen content of above-ground pool
	double nmass_agpool;


	/// year's maximum value of leaf C biomass
	double cmass_leaf_max;
	/// grs_cmass_leaf value saved at day before senescence (for LAI-calculation in allometry)
	double cmass_leaf_sen;

	/// today's increase of whole plant C biomass
	double dcmass_plant;
	/// today's increase of leaf C biomass
	double dcmass_leaf;
	/// today's increase of root C biomass
	double dcmass_root;
	/// today's increase of harvestable organ C biomass
	double dcmass_ho;
	/// today's increase of above-ground pool C biomass
	double dcmass_agpool;

	/// daily updated whole plant C biomass, reset at harvest day
	double grs_cmass_plant;
	/// daily updated leaf C biomass, reset at harvest day
	double grs_cmass_leaf;
	/// daily updated root C biomass, reset at harvest day
	double grs_cmass_root;
	/// daily updated harvestable organ C biomass, reset at harvest day
	double grs_cmass_ho;
	/// daily updated above-ground pool C biomass, reset at harvest day
	double grs_cmass_agpool;

	/// daily updated whole plant C biomass, reset at day 0
	double ycmass_plant;
	/// daily updated leaf C biomass, reset at day 0
	double ycmass_leaf;
	/// daily updated root C biomass, reset at day 0
	double ycmass_root;
	/// daily updated harvestable organ C biomass, reset at day 0
	double ycmass_ho;
	/// daily updated above-ground pool C biomass, reset at day 0
	double ycmass_agpool;

	/// year's whole plant C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_plant;
	/// year's leaf C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_leaf;
	/// year's root C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_root;
	/// year's harvestable organ C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_ho;
	/// year's above-ground pool C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_agpool;

	/// dry weight crop yield harvested this year (cumulative if several harvest events), based on harv_cmass_xx
	double harv_yield;

	/// harvestable organ C biomass at the last two harvest events this year
	double cmass_ho_harvest[2];
	/// dry weight crop yield at the last two harvest events this year
	double yield_harvest[2];	

	/// dry weight crop yield grown this year (cumulative if several harvest events), based on ycmass_xx
	double yield;				

	/// whether this pft is the main crop in the stand (pft.id==stand.pftid)
	bool isprimarycrop;
	/// whether this pft is allowed to compete with the main crop during the same growing period (for future use)
	bool isprimarycovegetation;
	/// whether this pft is grown during a second growing period, different from the primary (main) crop (for future use)
//	bool issecondarycrop;

	/// set to true if pft.isintercropgrass is true and the stand's main crop pft.intercrop is "naturalgrass"
	bool isintercropgrass;

	cropindiv_struct()
	{
		cmass_ho=0.0;
		cmass_agpool=0.0;
		nmass_ho=0.0;
		nmass_agpool=0.0;
		cmass_leaf_max=0.0;
		cmass_leaf_sen=0.0;
		yield=0.0;
		yield_harvest[0]=0.0;
		yield_harvest[1]=0.0;
		dcmass_leaf=0.0;
		dcmass_root=0.0;
		dcmass_plant=0.0;
		dcmass_ho=0.0;
		dcmass_agpool=0.0;
		grs_cmass_leaf=0.0;
		grs_cmass_root=0.0;
		grs_cmass_plant=0.0;
		grs_cmass_ho=0.0;
		grs_cmass_agpool=0.0;
		ycmass_leaf=0.0;
		ycmass_root=0.0;
		ycmass_plant=0.0;
		ycmass_ho=0.0;
		ycmass_agpool=0.0;
		harv_cmass_leaf=0.0;
		harv_cmass_root=0.0;
		harv_cmass_root=0.0;
		harv_cmass_ho=0.0;
		harv_yield=0.0;
		harv_cmass_agpool=0.0;
		cmass_ho_harvest[0]=0.0;
		cmass_ho_harvest[1]=0.0;

		isprimarycrop=false;
		isprimarycovegetation=false;
//		issecondarycrop=false;
		isintercropgrass=false;
	}

	void serialize(ArchiveStream& arch);
};


///////////////////////////////////////////////////////////////////////////////////////
// INDIVIDUAL
// State variables for a vegetation individual. In population mode this is the average
// individual of a PFT population; in cohort mode: the average individual of a cohort;
// in individual mode: an individual plant. Each grass PFT is represented as a single
// individual in all modes. Individual objects are collected within list arrays of
// class Vegetation (defined below), of which there is one for each patch, and include
// a reference to their 'parent' Vegetation object. Use the createobj member function
// of class Vegetation to add new individuals.

class Individual : public Serializable {

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

	/// nitrogen content of leaves on patch area basis (kgN/m2)
	double nmass_leaf;
	/// nitrogen content of roots on patch area basis (kgN/m2)	
	double nmass_root;
	/// nitrogen content of sapwood on patch area basis (kgN/m2)	
	double nmass_sap;
	/// nitrogen content of heartwood on patch area basis (kgN/m2)
	double nmass_heart;	

	double fpc;
		// foliar projective cover (FPC) under full leaf cover as fraction of modelled
		// area
	double fpc_daily;
		// foliar projective cover (FPC) this day as fraction of modelled
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

	/// Photosynthesis values for this individual under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;

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
	double lai_daily;
		// daily individual leaf area index (individual and cohort modes only)
	double lai_indiv_daily;
		// daily individual leaf area index (individual and cohort modes only)
	Historic<double, NYEARGREFF> greff_5;
		// growth efficiency (NPP/leaf area) for each of the last five simulation years
		// (kgC/m2/yr)
	double age;
		// individual/cohort age (years)
	double mlai[12];
		// monthly LAI (including phenology component)
	double mlai_max[12];
		// monthly maximum LAI (including phenology component)
	double fpar_leafon;
		// FPAR assuming full leaf cover for all vegetation
	double lai_leafon_layer;
		// LAI for current layer in canopy (cohort/individual mode; see function fpar)
	double gpterm;
		// non-water-stressed canopy conductance on FPC basis (mm/s)
	std::vector<double> gpterms;		// sub-daily version of the above variable (mm/s)
	double intercep;
		// interception associated with this individual today (patch basis)


	// Monthly sums (converted to means) maintained by function canopy_exchange

	double phen_mean;
		// accumulated mean fraction of potential leaf cover

	bool wstress; // whether individual subject to water stress

	/// leaf nitrogen that is photosyntetic active
	double nactive;
	/// Nitrogen extinction scalar
	/** Scalar to account for leaf nitrogen not following the optimal light 
	  * extinction, but is shallower.
	  */
	double nextin;
	/// long-term storage of labile nitrogen
	double nstore_longterm;
	/// storage of labile nitrogen
	double nstore_labile;
	/// daily total nitrogen demand
	double ndemand;
	/// fraction of individual nitrogen demand available for uptake
	double fnuptake;
	/// annual nitrogen uptake
	double anuptake;
	/// maximum size of nitrogen storage
	double max_n_storage;
	/// scales annual npp to maximum nitrogen storage
	double scale_n_storage;
	/// annual nitrogen limitation on vmax
	double avmaxnlim;
	/// annual optimal leaf C:N ratio
	double cton_leaf_aopt;
	/// annual average leaf C:N ratio
	double cton_leaf_aavr;
	/// plant mobile nitrogen status
	double cton_status;
	/// total carbon in compartments before growth
	double cmass_veg;
	/// total nitrogen in compartments before growth
	double nmass_veg;
	/// whether individual subject to nitrogen stress
	bool nstress;
	/// daily leaf nitrogen demand calculated from Vmax (kgN/m2)
	double leafndemand;
	/// daily root nitrogen demand based on leafndemand
	double rootndemand;
	/// daily sap wood nitrogen demand based on leafndemand
	double sapndemand;
	/// daily labile nitrogen demand based on npp
	double storendemand;
	/// leaf fraction of total nitrogen demand
	double leaffndemand;
	/// root fraction of total nitrogen demand
	double rootfndemand;
	/// sap fraction of total nitrogen demand
	double sapfndemand;
	/// store fraction of total nitrogen demand
	double storefndemand;
	/// daily leaf nitrogen demand over possible uptake (storage demand)
	double leafndemand_store;
	/// daily root nitrogen demand over possible uptake (storage demand)
	double rootndemand_store;
		
	int nday_leafon;	
		// Number of days with non-negligible phenology this month
	bool alive;
		// guess2008 - whether this individual is truly alive. Set to false for first year
		// after the Individual object is created, then true.
	double dnpp;

	// bvoc
	double iso; // isoprene production (mg C m-2 d-1)
	double mon; // monoterpene production (mg C m-2 d-1)
	double monstor; // monoterpene storage pool (mg C m-2)
	double fvocseas; // isoprene seasonality factor (-)

//private:
	cropindiv_struct *cropindiv;

	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v);
	~Individual();

//	const cropindiv_struct* get_cropindiv();
	cropindiv_struct* get_cropindiv();
	cropindiv_struct* set_cropindiv();

	void serialize(ArchiveStream& arch);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPFTFluxType flux_type, double value);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPatchFluxType flux_type, double value);

	inline bool istruecrop_or_intercropgrass() const {
		return (pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || cropindiv->isintercropgrass));
	}

	/// Reduce current biomass due to mortality and/or fire
	/** The removed biomass is put into litter pools and/or goes to fire fluxes.
	 *
	 *  \param mortality      fraction of Individual's biomass killed due to
	 *                        mortality (including fire)
	 *  \param mortality_fire fraction of Individual's biomass killed due to
	 *                        fire only
	 */
	void reduce_biomass(double mortality, double mortality_fire);

	/// Total storage of nitrogen
	double nstore() const {
		return nstore_longterm + nstore_labile;
	}

	/// Total carbon wood biomass
	double cmass_wood() const {
		return cmass_sap + cmass_heart - cmass_debt;
	}

	/// Total nitrogen wood biomass
	double nmass_wood() const {
		return nmass_sap + nmass_heart;
	}

	/// Current leaf C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_leaf(bool use_phen = true) const;

	/// Current fine root C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_root(bool use_phen = true) const;

	/// Current sap C:N ratio
	double cton_sap() const;

	/// Gets the individual's Patchpft
	Patchpft& patchpft();

	/// Transfers the individual's biomass (C and N) to litter and harvest pools/fluxes
	/** 
	 *  \param harvest Set to true if some of the biomass should be harvested,
	 *                 harvest will be done according to the PFT's harvest efficiency
	 *                 and residue outtake.
	 */
	void kill(bool harvest = false);

	/// Gets the individual's daily cmass_leaf value
	double cmass_leaf_phen() const {
		double cmass_leaf_phen = pft.phenology==CROPGREEN ? cropindiv->grs_cmass_leaf : cmass_leaf * phen;
		return cmass_leaf_phen;
	}
	/// Gets the individual's daily cmass_leaf value
	double cmass_root_phen() const {
		double cmass_root_phen = pft.phenology==CROPGREEN ? cropindiv->grs_cmass_root : cmass_root * phen;
		return cmass_root_phen;
	}
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

class Vegetation : public ListArray_idin2<Individual,Pft,Vegetation>, public Serializable {

public:
	// MEMBER VARIABLES

	Patch& patch; // reference to parent Patch object

	// MEMBER FUNCTIONS

	Vegetation(Patch& p):patch(p) {};
		// constructor (initialises member variable patch)

	void serialize(ArchiveStream& arch);
};


/// Soiltype stores static parameters for soils and the snow pack.
/** One Soiltype object is defined for each Gridcell. State variables for soils
 *  are held by objects of class Soil, of which there is one for each patch
 *  (see below).
 */
class Soiltype {

	// MEMBER VARIABLES

public:

	double awc_frac;
		// available water holding capacity as fraction of soil volume
	double awc[NSOILLAYER];
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
	double wp[NSOILLAYER];
		// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double wsats[NSOILLAYER];
		// saturation point. Cosby et al 1984
	int solvesom_end;
		// year at which to calculate equilibrium soil carbon
	int solvesom_begin;
		// year at which to begin documenting means for calculation of equilibrium
		// soil carbon

	/// water holding capacity plus wilting point for whole soil volume
	double wtot; 

	// For CENTURY ...
	/// fraction of soil that is sand
	double sand_frac;
	/// fraction of soil that is clay
	double clay_frac;
	/// fraction of soil that is silt plus clay	
	double silt_frac;

	// MEMBER FUNCTIONS

public:

	Soiltype() {

		// Constructor: initialises certain member variables

		solvesom_end = SOLVESOM_END;
		solvesom_begin = SOLVESOM_BEGIN;

		sand_frac = 0.4;
		clay_frac = 0.4;
		silt_frac = 0.2;
	}

	// guess2008 - override the default SOM years with 70-80% of the spin-up period length
	void updateSolveSOMvalues(const int& nyrspinup) {

		solvesom_end = static_cast<int>(0.8 * nyrspinup);
		solvesom_begin = static_cast<int>(0.7 * nyrspinup);

	}
};

/// CENTURY SOIL POOL
class Sompool : public Serializable {

public:

	/// Constructor
	Sompool() {
		
		// Initialise pool
		
		cmass = 0.0;
		nmass = 0.0;
		ligcfrac = 0.0;
		delta_cmass = 0.0;
		delta_nmass = 0.0;
		fracremain = 0.0;
		litterme = 0.0;
		fireresist = 0.0;

		for (int m = 0; m < 12; m++) {
			mfracremain_mean[m] = 0.0;
		}
	}

	/// C mass in pool kgC/m2
	double cmass;
	/// Nitrogen mass in pool kgN/m2
	double nmass;
	/// (potential) decrease in C following decomposition today (kgC/m2)
	double cdec; 
	/// (potential) decrease in nitrogen following decomposition today (kgN/m2)
	double ndec; 
	/// daily change in carbon and nitrogen
	double delta_cmass,delta_nmass;
	/// lignin fractions
	double ligcfrac;
	/// fraction of pool remaining after decomposition
	double fracremain;
	/// nitrogen to carbon ratio
	double ntoc;

	// Fire
	/// soil litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// soil litter fire resistance (0-1)
	double fireresist;

	// Fast SOM spinup variables

	/// monthly mean fraction of carbon pool remaining after decomposition
	double mfracremain_mean[12];

	void serialize(ArchiveStream& arch);
};

/// This struct contains litter for solving Century SOM pools.
/** \see equilsom() */  
struct LitterSolveSOM : public Serializable {
	/// Constructs an empty result
	LitterSolveSOM() {
		clear();
	}

	/// Clears all members
	void clear() {
		for (int p = 0; p < NSOMPOOL; p++) {
			clitter[p] = 0.0;
			nlitter[p] = 0.0;
		}
	}

	/// Add litter
    void add_litter(double cvalue, double nvalue, int pool) {
		clitter[pool] += cvalue;
		nlitter[pool] += nvalue;
    }

	double get_clitter(int pool) {
		return clitter[pool];
	}
	double get_nlitter(int pool) {
		return nlitter[pool];
	}

	void serialize(ArchiveStream& arch);

private:
	// Carbon litter
	double clitter[NSOMPOOL];
	
	// Nitrogen litter
	double nlitter[NSOMPOOL];
};

/// Soil stores state variables for soils and the snow pack.
/** Initialised by a call to initdrivers. One Soil object is defined for each patch.
 *  A reference to the parent Patch object (defined below) is included as a member
 *  variable. Soil static parameters are stored as objects of class Soiltype, of which
 *  there is one for each grid cell. A reference to the Soiltype object holding the
 *  static parameters for this soil is included as a member variable.
 */
class Soil : public Serializable {

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

	double alag, exp_alag;


	// guess2008 - 3 new soil water variables
	double mwcont[12][NSOILLAYER];
		// water content of soil layers [0=upper layer] as fraction of available water
		// holding capacity;
	double dwcontlower[365];
		// daily water content in lower soil layer for each day of year
	double mwcontlower;
		// mean water content in lower soil layer for last month
		// (valid only on last day of month following call to daily_accounting_patch)

	double rain_melt;						// rainfall and snowmelt today (mm)
	double max_rain_melt;					// upper limit for percolation (mm)
	bool percolate;							// whether to percolate today

//////////////////////////////////////////////////////////////////////////////////
// CENTURY SOM pools and other variables

	Sompool sompool[NSOMPOOL];

	/// daily percolation (mm)
	double dperc;
	/// fraction of decayed organic nitrogen leached each day;
	double orgleachfrac;
	/// soil mineral nitrogen pool (kgN/m2)
	double nmass_avail;		
	/// soil nitrogen input (kgN/m2)
	double ninput;
	/// annual sum of nitrogen mineralisation
	double anmin;			
	/// annual sum of nitrogen immobilisation
	double animmob;			
	/// annual leaching from available nitrogen pool
	double aminleach;		
	/// annual leaching of organics from active nitrogen pool
	double aorgleach;		
	/// total annual nitrogen fixation 
	double anfix;
	/// calculated annual mean nitrogen fixation
	double anfix_calc;
	
	// Variables for fast spinup of SOM pools

	/// monthly fraction of available mineral nitrogen taken up
	double fnuptake_mean[12];
	/// monthly fraction of organic carbon/nitrogen leached
	double morgleach_mean[12];
	/// monthly fraction of available mineral nitrogen leached
	double mminleach_mean[12];
	/// annual nitrogen fixation
	double anfix_mean;

	// Solving Century SOM pools 

	/// years at which to begin documenting for calculation of Century equilibrium
	int solvesomcent_beginyr;
	/// years at which to end documentation and start calculation of Century equilibrium
	int solvesomcent_endyr;

	std::vector<LitterSolveSOM> solvesom;

	/// stored nitrogen deposition in snowpack
	double snowpack_nmass;

	// MEMBER FUNCTIONS

public:
	/// constructor (initialises member variable patch)
	Soil(Patch& p,Soiltype& s):patch(p),soiltype(s) {
			initdrivers();
	}

	void initdrivers() {

		// Initialises certain member variables

		alag = 0.0;
		exp_alag = 1.0;
		cpool_slow = 0.0;
		cpool_fast = 0.0;
		decomp_litter_mean = 0.0;
		k_soilfast_mean = 0.0;
		k_soilslow_mean = 0.0;
		wcont[0] = 0.0;
		wcont[1] = 0.0;
		wcont_evap = 0.0;
		snowpack = 0.0;
		orgleachfrac = 0.0;


		// guess2008 - extra initialisation
		mwcontupper = 0.0;
		mwcontlower = 0.0;
		for (int mth=0; mth<12; mth++) {
			mwcont[mth][0] = 0.0;
			mwcont[mth][1] = 0.0;
			fnuptake_mean[mth] = 0.0;
			morgleach_mean[mth] = 0.0;
			mminleach_mean[mth] = 0.0;
		}

		for (int d=0; d<365; d++) {
			dwcontupper[d] = 0.0;
			dwcontlower[d] = 0.0;
		}

		/////////////////////////////////////////////////////
		// Initialise CENTURY pools

		// Set initial CENTURY pool N:C ratios 
		// Parton et al 1993, Fig 4

		sompool[SOILMICRO].ntoc = 1.0 / 15.0;
		sompool[SURFHUMUS].ntoc = 1.0 / 15.0;
		sompool[SLOWSOM].ntoc = 1.0 / 20.0;
		sompool[SURFMICRO].ntoc = 1.0 / 20.0;

		// passive has a fixed value
		sompool[PASSIVESOM].ntoc = 1.0 / 9.0;

		nmass_avail = 0.0;
		ninput = 0.0;
		anmin = 0.0;			
		animmob = 0.0;		
		aminleach = 0.0;
		aorgleach = 0.0;
		anfix = 0.0;
		anfix_calc = 0.0;
		anfix_mean = 0.0;
		snowpack_nmass = 0.0;
		dperc = 0.0;

		solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
		solvesomcent_endyr   = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);
	}

	void serialize(ArchiveStream& arch);
};


class cropphen_struct : public Serializable {

public:

	/// latest sowing date
	int sdate;

	/// sowing date of growing period ending in latest harvest this year
	int sdate_harv;

	/// sowing dates of growing periods ending in the two latest harvests this year
	int sdate_harvest[2];

	/// sowing dates of growing periods starting this year
	int sdate_thisyear[2];

	/// number of sowings this year
	int nsow;

	/// latest harvest date
	int hdate;

	/// two latest harvest dates this year
	int hdate_harvest[2];

	/// last date for harvest
	int hlimitdate;

	/// last day of heat unit sampling period, set in Crop_sowing_date_new()
	int hucountend;

	/// number of harvests this year
	int nharv;

	/// whether sdate_harvest[0] happened last year
	bool sownlastyear;

	/// latest senescence start date this year
	int sendate;

	/// latest beginning of intercropseason (2 weeks after the harvest date)
	int bicdate;

	/// latest end of intercropseason (2 weeks before the sowing date)
	int eicdate;

	/// number of growing days this growing period
	int growingdays;

	/// number of growing days this year (used for wscal_mean calculation)
	int growingdays_y;

	/// length of growingseason ending in last harvest
	int lgp;

	/// base temp for heat unit calculation (°C)
	double tb;

	/// number of vernalising days required
	int pvd;

	/// number of accumulated vernalizing days
	int vdsum;

	/// heat unit reduction factor due to vernalization [0-1]
	double vrf;

	/// heat unit reduction factor due to photoperiodism [0-1]
	double prf;

	/// potential heat units required for crop maturity (°Cd)
	double phu;

	/// potential heat units that would have been used without dynamic phu calculation
	double phu_old;

	/// heat unit sum aquired during last growing period (°Cd)
	double husum;

	/// heat unit sum aquired durin sampling period, starting with sdate
	double husum_sampled;

	/// this year's heat unit sum aquired from sdate to hucountend
	double husum_max;

	/// running mean of recent past's husum_max
	double husum_max_10;

	/// number of heat units sampling years
	int nyears_hu_sample;

	/// fraction of growing season [0-1] (husum/phu)
	double fphu;

	/// fraction of growing season at latest harvest
	double fphu_harv;

	/// fraction of growing season at the two latest harvests this year
//	double fphu_harvest[2];	

	bool hu_samplingperiod;

	int hu_samplingdays;

	/// harvest index today [0-1, >1 if below-ground ho], harvestable organ/above-ground C for above-ground harvestable organs, dependent on fphu, reduced by water stress 
	double hi;

	/// harvest index yesterday
	double hi_ystd;

	/// fraction of harvest index today
	double fhi;

	/// phenology (fphu) contribution of fraction of harvest index today
	double fhi_phen;	//Phenology (fPHU) compoment of fhi

	/// water stress contribution of fraction of harvest index today
	double fhi_water;

	/// fraction of harvest index at latest harvest
	double fhi_harv;

	/// acheived fraction of harvest index at the two latest harvests this year
//	double fhi_harvest[2];

	/// sum of crop patch demand (patch.wdemand) during crop growing period, reset on harvest day
	double demandsum_crop;

	/// sum of crop supply (patchpft.wsupply) during crop growing period, reset on harvest day
	double supplysum_crop;

	/// copy of indiv.lai, could be used in crop_phenology() to set phen
	double lai;

	/// daily lai value, set in lai_crop()
	double lai_daily;

	/// copy of indiv.fpc, used in crop_phenology() to set phen 
	double fpc;

	/// daily fpc value, set in crop_phenology()
	double fpc_daily;

	/// whether inside crop/intercrop grass growing period
	bool growingseason;

	/// whether yesterday was inside crop/intercrop grass growing period
	bool growingseason_ystd;

	/// whether inside crop senescence
	bool senescence;

	/// whether yesterday was inside crop senescence
	bool senescence_ystd;

	/// whether inside intercrop crass growing period (main crop pft variable)
	bool intercropseason;

	/// used to distinguish the two growing seasons for rice in some geographical regions
	bool maincrop;

	cropphen_struct()
	{
		lai=0.0;
		lai_daily=0.0;
		fpc=0.0;
		fpc_daily=0.0;
		sdate=-1;
		sdate_harv=-1;
		nsow=0;
		sownlastyear=false;	
		sendate=-1;		
		hdate=-1;				
		hlimitdate=-1;
		hucountend=-1;
		nharv=0;
		tb=0.0;
		pvd=0;
		vdsum=0;
		vrf=1.0;
		phu=0.0;
		phu_old=0.0;
		husum_max=0.0;
		husum_sampled=0.0;
		husum_max_10=0.0;
		nyears_hu_sample = 0;
		prf=1.0;
		husum=0.0;
		fphu=0.0;
		fphu_harv=0.0;
		hu_samplingdays=0;
		hu_samplingperiod=false;

		hi=0.0;
		hi_ystd=0.0;
		fhi=0.0;
		fhi_phen=0.0;
		fhi_water=1.0;
		fhi_harv=0.0;
		demandsum_crop=0.0;
		supplysum_crop=0.0;

		growingseason=false;	//Initialized to true for normal grass growth (CC3G & CC4G) in establishment
		growingseason_ystd=false;
		senescence=false;
		senescence_ystd=false;
		intercropseason=false;
		bicdate=-1;
		eicdate=-1;
		maincrop=true;
		growingdays=0;
		growingdays_y=0;
		lgp=0;

		for(int j=0;j<2;j++)			
		{
			sdate_harvest[j]=-1;	
			hdate_harvest[j]=-1;
			sdate_thisyear[j]=-1;		
//			fhi_harvest[j]=-1.0;		
//			fphu_harvest[j]=-1.0;		
		}
	}

	void serialize(ArchiveStream& arch);
};


///////////////////////////////////////////////////////////////////////////////////////
// PATCHPFT
// State variables common to all individuals of a particular PFT in a particular patch
// Used in individual and cohort modes only.

class Patchpft : public Serializable {

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
	double litter_sap;
		// sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_heart;
		// heartwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_repr;
		// litter derived from allocation to reproduction for PFT on modelled area
		// basis (kgC/m2)

	/// leaf-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_leaf;
	/// root-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_root;
	/// sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_sap;
	/// heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_heart;

	double gcbase;
		// non-FPC-weighted canopy conductance value for PFT under water-stress
		// conditions (mm/s)
	double gcbase_day;				// daily value of the above variable (mm/s)

	double wsupply;
		// evapotranspirational "supply" function for this PFT today (mm/day)
	double wsupply_leafon;
	double fwuptake[NSOILLAYER];
		// fractional uptake of water from each soil layer today
	bool wstress;				// whether water-stress conditions for this PFT
	bool wstress_day;			// daily version of the above variable

	/// carbon depository for long-lived products like wood
	double harvested_products_slow;	
	/// nitrogen depository for long-lived products like wood
	double harvested_products_slow_nmass; 
	/// first and last day of crop sowing window, calculated in crop_sowing_patch() or Crop_sowing_date_new()
	int swindow[2];
	/// daily value of water deficit, calculated in irrigated_water_uptake()
	double water_deficit_d;
	/// yearly sum of water deficit
	double water_deficit_y;

//private:

	/// Struct for crop-specific variables
	cropphen_struct *cropphen;

	// MEMBER FUNCTIONS:

	Patchpft(int i,Pft& p):id(i),pft(p) {

		// Constructor: initialises id, pft and data members

		litter_leaf = 0.0;
		litter_root = 0.0;
		litter_sap   = 0.0;
		litter_heart = 0.0;
		litter_repr = 0.0;

		nmass_litter_leaf  = 0.0;
		nmass_litter_root  = 0.0;
		nmass_litter_sap   = 0.0;
		nmass_litter_heart = 0.0;

		wscal = 1.0;
		wscal_mean = 1.0;
		anetps_ff = 0.0;
		aphen = 0.0;
		phen = 0.0;
		wsupply = 0.0;
		wsupply_leafon = 0.0;

		for(int i=0;i<NSOILLAYER;i++)
			fwuptake[i]=0.0;

		cropphen = NULL;
		harvested_products_slow = 0.0;
		harvested_products_slow_nmass = 0.0;

		swindow[0]=-1;
		swindow[1]=-1;

		if(pft.landcover==CROPLAND)
		{
			cropphen=new cropphen_struct;
		}
	}

	~Patchpft()
	{
		if(cropphen)
			delete cropphen;
	}

	/// safe method to obtain cropphen_struct pointer
	cropphen_struct* get_cropphen();
	/// safe method to set cropphen_struct variables
	cropphen_struct* set_cropphen();

	void serialize(ArchiveStream& arch);
};


///////////////////////////////////////////////////////////////////////////////////////
// PATCH
// Stores data for a patch. In cohort and individual modes, replicate patches are
// required in each stand to accomodate stochastic variation; in population mode there
// should be just one Patch object, representing average conditions for the entire
// stand. A reference to the parent Stand object (defined below) is included as a
// member variable.

class Patch : public Serializable {

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
	Historic<double, NYEARAAET> aaet_5;
		// annual sum of AET (mm/year) for each of the last five simulation years
	double aevap;
		// annual sum of soil evaporation (mm/year)
	double aintercep;
		// annual sum of interception (mm/year)
	double asurfrunoff;
		// annual sum of runoff (mm/year)
	double adrainrunoff;
		// annual sum of runoff (mm/year)
	double abaserunoff;
		// annual sum of runoff (mm/year)
	double arunoff;
		// annual sum of runoff (mm/year)
	double apet;
		// annual sum of potential evapotranspiration (mm/year)

	double eet_net_veg;
		// equilibrium evapotranspiration today, deducting interception (mm)

	double wdemand;
		// transpirative demand for patch, patch vegetative area basis (mm/day)
	double wdemand_day;			
		// daily average of the above variable (mm/day)
	double wdemand_leafon;
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

	/// daily nitrogen demand
	double ndemand;

	/// daily value of irrigation water, set in irrigation(), derived from water_deficit_d
	double irrigation_d;
	/// yearly sum of irrigation water
	double irrigation_y;

	// MEMBER FUNCTIONS

	Patch(int i,Stand& s,Soiltype& st):
		id(i),stand(s),vegetation(*this),soil(*this,st),fluxes(*this) {

		// Constructor: initialises various members and builds list array
		// of Patchpft objects.

		for(unsigned int p = 0; p < pftlist.nobj; p++) {
			pft.createobj(pftlist[p]);
		}

		age = 0;
		disturbed = false;
		wdemand = 0.0;
		wdemand_leafon = 0.0;
		
		// guess2008 - initialise
		growingseasondays = 0;

		fireprob = 0.0;
		ndemand = 0.0;
	}

	void serialize(ArchiveStream& arch);
};

/// Container for variables common to individuals of a particular PFT in a stand.
/** Used in individual and cohort modes only
 */
class Standpft : public Serializable {

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
	double gpterm;
		// non-FPAR-weighted value for canopy conductance component associated with
		// photosynthesis for PFT under non-water-stress conditions (mm/s)
	std::vector<double> gpterms;		// sub-daily version of the above variable (mm/s)

	double fpc_total;
		// FPC sum for this PFT as average for stand (used by some versions of
		// guessio.cpp)

	/// Photosynthesis values for this PFT under non-water-stress conditions
	PhotosynthesisResult photosynthesis;
	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;

	/// Is this PFT allowed to grow in this stand?
	bool active;

	// MEMBER FUNCTIONS

	Standpft(int i,Pft& p):id(i),pft(p) {

		// Constructor: initialises various data members
		anetps_ff_max = 0.0;
		active = !run_landcover;
	}

	void serialize(ArchiveStream& arch);
};


/// The stand class corresponds to a modelled area of a specific landcover type in a grid cell.
/** There may be several stands of the same landcover type (but with different settings).
 */
class Stand : public ListArray_idin2<Patch,Stand,Soiltype>, public Serializable {

public:

	// MEMBER VARIABLES

	/// list array [0...npft-1] of Standpft (initialised in constructor)
	ListArray_idin1<Standpft,Pft> pft;

	/// A number identifying this Stand within the grid cell
	int id;

	/// pft id of main crop
	int pftid;

	/// crop id of main crop pft (crop pft:s numbered 0-25, corresponding to columns 0-25 in crop cover input file)
	int cftid;

	// true if main crop pft.hydrology==irrigated
	bool isirrigated;

	/// true if the stand's main crop pft intercrop==naturalgrass and a pft with isintercrop==true is in the pftlist.
	bool hasgrassintercrop;

	/// gdd5_pasture-value at first intercrop grass growth
	double gdd0_intercrop;
	
	/// fraction removed from natural stand when converted to other landcover type
	double natural_frac_change;

	/// Seed for generating random numbers within this Stand
	/** The reason why Stand has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when using
	 *  different land cover types.
	 *
	 *  Randomness not associated with a specific stand, but rather a whole
	 *  grid cell should instead use the seed in the Gridcell class.
	 *
	 *  \see randfrac()
	 */
	long seed;

	/// reference to parent object
	Gridcell& gridcell;

	/// type of landcover
	/** \see landcovertype
	 *  initialised in constructor
	 */
	landcovertype landcover;

	/// The year when this stand was created.
	/** Will typically be year zero unless running with dynamic
	 *  land cover.
	 *
	 *  Needed to set patchpft.anetps_ff_est_initial
	 */
	int first_year;

	/// scaling factor for stands that have grown in area this year (old fraction/new fraction)
	double scale_LC_change;

	// MEMBER FUNCTIONS

	/// Constructs a Stand
	/** \param i         The id for the stand within the grid cell
	 *  \param gc        The parent grid cell
	 *  \param landcover The type of landcover to use for this stand
	 */
	Stand(int i, Gridcell& gc,landcovertype landcover); 

	/// Gives the fraction of this Stand relative to the whole grid cell
	double get_gridcell_fraction() const;

	/// Gives the fraction of this Stand relative to its land cover type; NB: unsafe to use within landcover_dynamics() !
	double get_landcover_fraction() const;

	/// Set the fraction of this Stand relative to the gridcell
	void set_gridcell_fraction(double fraction);

	/// Returns the number of patches in this Stand
	unsigned int npatch() const { return nobj; }

	/// Returns true if stand is tru crop stand, as opposed to pasture grass grown on cropland
	inline bool is_true_crop_stand() {
		return landcover==CROPLAND && pft[pftid].pft.phenology==CROPGREEN;
	}

	bool ifnlim_stand() const {
		return ifnlim && ifnlim_lc[landcover];
	}

	void serialize(ArchiveStream& arch);

private:

	/// Fraction of this stand relative to the gridcell
	/** used by crop stands; initialized in constructor to 1,
	 *  set in landcover_init()
	 */
	double frac;
};



/// State variables common to all individuals of a particular PFT in a GRIDCELL.
class Gridcellpft : public Serializable {

public:

	// MEMBER VARIABLES

	/// A number identifying this object within its list array
	int id;

	/// A reference to the Pft object for this Gridcellpft
	Pft& pft;

	/// annual degree day sum above threshold damaging temperature
	/** used in calculation of heat stess mortality; Sitch et al 2000, Eqn 55
	 */
	double addtw;

	/// Michaelis-Menten kinetic parameters 
	/** Half saturation concentration for N uptake (Rothstein 2000, Macduff 2002) 
	 */
	double Km;

	///Crop-specific variables:
	/// whether the daily temperature has fallen below the autumn temperature limit (tempautumn) this year
	bool autumnoccurred;
	/// whether the daily temperature has risen above the spring temperature limit (tempspring) this year
	bool springoccurred;
	/// whether the daily temperature has fallen below the vernalization limit (trg) this year
	bool vernstartoccurred;
	/// whether the daily temperature rises over the vernalization limit (trg) this year
	bool vernendoccurred;
	/// whether the precipitation sum rises over the limit (defined in set_sdatecalc_prec), not used if NEWSOWINGDATE is defined
	bool precoccurred;
	/// first day when temperature fell below the autumn temperature limit (tempautumn) this year
	int first_autumndate;
	/// 20-year mean
	int first_autumndate20;
	/// memory of the last 20 years' values
	int first_autumndate_20[20];
	/// last day when temperature rose above the spring temperature limit (tempspring) this year
	int last_springdate;
	/// 20-year mean
	int last_springdate20;
	/// memory of the last 20 years' values
	int last_springdate_20[20];
	/// last day when temperature has fallen below the vernilisation temperature limit (trg) this year (if vernstartoccurred==true)
	int last_verndate;
	/// 20-year mean
	int last_verndate20; 
	/// memory of the last 20 years' values
	int last_verndate_20[20];
	/// first day when 10-day precipitation sum rises over the limit (defined in set_sdatecalc_prec), not used if NEWSOWINGDATE is defined
	int first_precdate;
	/// default sowing date (pft.sdatenh/sdatesh)
	int sdate_default;
	/// calculated sowing date from temperature limits
	int sdatecalc_temp;
	/// calculated sowing date from precipitation limits
	int sdatecalc_prec;
	/// sowing date from input file
	int sdate_force;
	/// harvest date from input file
	int hdate_force;
	/// default harvest date (pft.hlimitdatenh/hlimitdatesh)
	int hlimitdate_default;
	/// whether autumn sowing is either calculated or prescribed
	bool wintertype;
	/// whether only one sowing season per year is allowed (currently used only for rice, based on geograhical limits defined in getgridcell() )
	bool singlecrop;
	/// first and last day of crop sowing window, calculated in calc_sowing_windows()
	int swindow[2];

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellpft object
	/** \param i   The id for this object
	 *  \param p   A reference to the Pft for this Gridcellpft
	 */
	Gridcellpft(int i,Pft& p):id(i),pft(p) {
		addtw = 0.0;
		Km = 0.0;

		autumnoccurred=false;
		springoccurred=false;
		vernstartoccurred=false;
		vernendoccurred=false;
		precoccurred=false;
		first_autumndate=-1;
		first_autumndate20=-1;
		last_springdate=-1;
		last_springdate20=-1;
		last_verndate=-1;
		last_verndate20=-1;
		first_precdate=-1;
		for (int year=0;year<20;year++) 
		{
			first_autumndate_20[year]=-1;
			last_springdate_20[year]=-1;
			last_verndate_20[year]=-1;
		}
		sdate_default=-1;
		sdate_force=-1;
		hdate_force=-1;
		sdatecalc_temp=-1;
		sdatecalc_prec=-1;
		hlimitdate_default=-1;
		wintertype=false;
		singlecrop=true;
		swindow[0]=-1;
		swindow[1]=-1;
	}

	void serialize(ArchiveStream& arch);
};


/// The Gridcell class corresponds to a modelled locality or grid cell.
/** Member variables include an object of type Climate (holding climate, insolation and
 *  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
 *  array of Stand objects. Soil objects (holding soil state variables) are associated
 *  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
 *  locality or grid cell.
 */
class Gridcell : public ListArray_idin2<Stand,Gridcell,landcovertype>, public Serializable {

public:

	// MEMBER VARIABLES

	/// climate, insolation and CO2 for this grid cell
	Climate climate;

    /// soil static parameters for this grid cell
	Soiltype soiltype;

	/// The fractions of the different land cover types.
	/** landcoverfrac is read in from land cover input file or from
	 *  instruction file in getlandcover().
	 */
	double landcoverfrac[NLANDCOVERTYPES];
	double cftfrac[NCROPSTANDS_MAX];

	/// The land cover fractions from the previous year
	/** Used to keep track of the changes when running with dynamic
	 *  land cover.
	 */
	double landcoverfrac_old[NLANDCOVERTYPES];
	double cftfrac_old[NCROPSTANDS_MAX];

	/// Whether the land cover fractions changed for this grid cell this year
	/** \see landcover_dynamics
	 */
	bool LC_updated;

	/// Gridcell-level flux from slow harvested products
	double acflux_harvest_slow;
	/// Gridcell-level flux from harvest associated with landcover change
	double acflux_landuse_change;
	/// Landcover-level flux from slow harvested products (donating landcover)
	double acflux_harvest_slow_lc[NLANDCOVERTYPES];	
	/// Landcover-level flux from harvest associated with landcover change (donating landcover)
	double acflux_landuse_change_lc[NLANDCOVERTYPES];

	/// list array [0...npft-1] of Gridcellpft (initialised in constructor)
	ListArray_idin1<Gridcellpft,Pft> pft;

	/// Seed for generating random numbers within this Gridcell
	/** The reason why Gridcell has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when for
	 *  instance changing the order in which the simulation proceeds. It also
	 *  gets serialized together with the rest of the Gridcell state to make it
	 *  possible to get exactly identical results after a restart.
	 *
	 *  \see randfrac()
	 */
	long seed;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcell object
	Gridcell():climate(*this) {
		landcovertype landcover;
		LC_updated = false;

		for(unsigned int p=0; p<pftlist.nobj; p++) {
			pft.createobj(pftlist[p]);
		}

		memset(landcoverfrac, 0, sizeof(double) * NLANDCOVERTYPES);
		memset(landcoverfrac_old, 0, sizeof(double) * NLANDCOVERTYPES);
		memset(cftfrac, 0, sizeof(double) * NCROPSTANDS_MAX);
		memset(cftfrac_old, 0, sizeof(double) * NCROPSTANDS_MAX);
		acflux_harvest_slow=0.0;
		acflux_landuse_change=0.0;
		memset(acflux_harvest_slow_lc, 0, sizeof(double)*NLANDCOVERTYPES);
		memset(acflux_landuse_change_lc, 0, sizeof(double)*NLANDCOVERTYPES);

		if(!run_landcover) {
			landcover = NATURAL;
			createobj(*this,landcover);
			landcoverfrac[NATURAL] = 1.0;
		}

		seed = 12345678;
	}

	/// Creates new stand
	void create_stand_lu(landcovertype landcover, double fraction, int cftid=-2);

	/// Longitude for this grid cell
	double get_lon() const;

	/// Latitude for this grid cell
	double get_lat() const;

	/// Set longitude and latitude for this grid cell
	void set_coordinates(double longitude, double latitude);

	void serialize(ArchiveStream& arch);

private:

	/// Longitude for this grid cell
	double lon;

	/// Latitude for this grid cell
	double lat;

};


#endif // LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995. 
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991 
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical 
//   exploration of the relationships of soil moisture characteristic to the 
//   physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Franzlubbers, AJ & Stuedemann, JA 2009 Soil-profile organic carbon and total 
//   nitrogen during 12 years of pasture management in the Southern Piedmont USA. 
//   Agriculture Ecosystems & Environment, 129, 28-36.
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A 
//   process-based, terrestrial biosphere model of ecosystem dynamics 
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Macduff, JH, Humphreys, MO & Thomas, H 2002. Effects of a stay-green mutation on
//   plant nitrogen relations in Lolium perenne during N starvation and after 
//   defoliation. Annals of Botany, 89, 11-21.
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Reich, PB, Walters MB & Ellsworth DS 1992 Leaf Life-Span in Relation to Leaf,
//   Plant, and Stand Characteristics among Diverse Ecosystems. 
//   Ecological Monographs 62: 365-392.
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Sykes, MT, Prentice IC & Cramer W 1996 A bioclimatic model for the potential
//   distributions of north European tree species under present and future climates.
//   Journal of Biogeography 23: 209-233.
// White, M A, Thornton, P E, Running, S. & Nemani, R 2000 Parameterization and 
//   Sensitivity Analysis of the BIOME-BGC Terrestrial Ecosystem Model: Net Primary 
//   Production Controls. Earth Interactions, 4, 1-55.
