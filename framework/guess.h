///////////////////////////////////////////////////////////////////////////////////////
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
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESS_H
#define LPJ_GUESS_GUESS_H

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

#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "gutil.h"
#include <vector>
#include <algorithm>
#include <set>
#include "shell.h"
#include "guessmath.h"
#include "archive.h"
#include "parameters.h"
#include "guesscontainer.h"
#include "soil.h"
#include "guessstring.h"
#include "interpolator.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

/// Life form class for PFTs (trees, grasses, mosses)
typedef enum {NOLIFEFORM, TREE, GRASS_PERENNIAL, GRASS_ANNUAL, MOSS} lifeformtype;

/// Phenology class for PFTs
typedef enum {NOPHENOLOGY, EVERGREEN, RAINGREEN, SUMMERGREEN, CROPGREEN, ANY} phenologytype;

/// Biochemical pathway for photosynthesis (C3 or C4)
typedef enum {NOPATHWAY, C3, C4} pathwaytype;

/// Leaf physiognomy types for PFTs
typedef enum {NOLEAFTYPE, NEEDLELEAF, BROADLEAF} leafphysiognomytype;

/// The level of verbosity of LPJ-GUESS. Decides the amount of information that is written to the log-file.
typedef enum {ERROR, WARNING, INFO, DEBUG_WARNING} verbositylevel;

/// Names of leaf pools in DAVE version
typedef enum {green, brown} leafpool;

/// Array size of mowdates
const int num_mowdates = 5;

/// Soil input type
typedef enum { LPJ, MINERAL, VARYING } soilinputtype;

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
typedef enum { SURFMETA, SURFSTRUCT, SURFFWD, SURFCWD, SURFMICRO, SURFHUMUS, 
SOILMETA, SOILSTRUCT, SOILFWD, SOILCWD, SOILMICRO, SLOWSOM, PASSIVESOM, NSOMPOOL } pooltype;

/// CENTURY SOM pool group names
typedef enum { GALL, GSURFL, GHUMUS, GSOM, GHUMUSSOM } grouppooltype;

/// Ways of calculating mineral fractions
typedef enum { ABS, ABSEFF, MINEFF } mincalctype;

/// Irrigation type for PFTs
typedef enum {RAINFED, IRRIGATED} waterusetype;
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

/// Gas type (used in methane code)
/** 
  */
typedef enum {O2gas, CO2gas, CH4gas} gastype;

/// Nitrogen preferance
typedef enum {NO, NH4, NO3} n_pref_type;

/// Gross landcover change input type
typedef enum {NO_GROSS_INPUT, LC_GROSS_INPUT, ST_GROSS_INPUT} gross_LCC_input_type;

/// thinselectpft types (integers in instruction file)
#define SELECT_PFT_NOPREF 0
#define SELECT_PFT_UNSEL 1
#define SELECT_PFT_SEL 2
#define SELECT_PFT_UNSEL_SEL_SEPARATE 3
#define SELECT_PFT_SHADEINTOL 4

/// thinselectage types (integers in instruction file)
#define SELECT_AGE_NOPREF 0
#define SELECT_AGE_YOUNG 1
#define SELECT_AGE_OLD 2

/// thinselectdiam types (integers in instruction file)
#define SELECT_DIAM_NOPREF 0
#define SELECT_DIAM_SMALL 1
#define SELECT_DIAM_LARGE 2
#define SELECT_DIAM_LIMIT 3

/// Number of forest structure (age- and diameter) classes
const int NFOREST_STRUCTURAL_CLASSES = 31;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

/// maximum number of thinnings in a forest management scheme
const int NTHINNINGS = 5;
/// maximum number of thinning loops in a forest management scheme (e.g. regrowth and continuous loops)
const int NTHINNINGLOOPS = 2;

/// bvoc: number of monoterpene species used
const int NMTCOMPOUNDS=NMTCOMPOUNDTYPES;

// SOIL DEPTH VALUES

/// Maximum soil depth - 50m (mm)
const double SOILDEPTH_MAX = 50000.0;

/// Minimum soil depth - 0.2m (mm)
const double SOILDEPTH_MIN = 200.0;

/// soil upper layer depth (mm)
const double SOILDEPTH_UPPER = 500.0;
/// soil lower layer depth (mm)
const double SOILDEPTH_LOWER = 1000.0;

/// Depth of sublayer at top of upper soil layer, from which evaporation is
//  possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
//  Must be a multiple of Dz_soil
const double SOILDEPTH_EVAP_TRUNK = 200.0;
extern double SOILDEPTH_EVAP;

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

/// Warmest day in N hemisphere (same as COLDEST_DAY_SHEMISPHERE)
const int WARMEST_DAY_NHEMISPHERE = COLDEST_DAY_SHEMISPHERE;

/// Warmest day in S hemisphere (same as COLDEST_DAY_NHEMISPHERE)
const int WARMEST_DAY_SHEMISPHERE = COLDEST_DAY_NHEMISPHERE;

/// number of years to average aaet over in function soilnadd
const int NYEARAAET = 5;

/// number of years to average max snow depth over in function soilnadd
const int NYEARMAXSNOW = 20;

/// Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
const double PRIESTLEY_TAYLOR = 1.32;

// Solving Century SOM pools

/// Fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
const double SOLVESOMCENT_SPINBEGIN = 0.1;
/// Fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
const double SOLVESOMCENT_SPINEND = 0.3;

/// Kelvin to deg C conversion
const double K2degC = 273.15;

// The gas constant [J⋅mol−1 K-1]
const double R = 8.3145;

/// Maximum number of crop rotation items
const int NROTATIONPERIODS_MAX = 3;

/// Conversion factor for CO2 from ppmv to mole fraction
const double CO2_CONV = 1.0e-6;

/// Initial carbon allocated to crop organs at sowing, kg m-2
const double CMASS_SEED = 0.01;

/// Averaging interval for average maximum annual fapar (SIMFIRE)
const int AVG_INTERVAL_FAPAR = 3;

/// Averaging interval for biome averaging (SIMFIRE)
const int N_YEAR_BIOMEAVG = 3;

/// Default year in the future never expected to be covered during the simulation
const int FAR_FUTURE_YEAR = 100000;
/// Default year in the past never expected to be covered during the simulation
const int FAR_PREHISTORIC_YEAR = -100000;

// Max error allowed in water balance (mm). A value of 1e-9 is too strict, and
// results in mass balance errors, likely due to the many additions and
// subtractions of small floating point numbers when we aggregate across the
// layers of the soil. A value of 1e-4 was used in the trunk version, so if this
// value turns out to be too strict, we should relax it back to 1e-4.
const double EPS_WATER = 1e-6;

///////////////////////////////////////////////////////////////////////////////////////
// FORWARD DECLARATIONS OF CLASSES DEFINED IN THIS FILE
// Forward declarations of classes used as types (e.g. for reference variables in some
// classes) before they are actually defined

class Date;
class Stand;
class Patch;
class Vegetation;
class Individual;
class Gridcell;
class Patchpft;
class Soil;
class Soillayer;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are defined in the framework source code file, and are accessible
// throughout the code

/// Object describing timing stage of simulation
extern Date date;

/// Number of possible PFTs
extern int npft;
/// Number of stand types in stlist
extern int nst;
/// Number of stand types per land cover
extern int nst_lc[NLANDCOVERTYPES];
/// Number of management types in stlist
extern int nmt;

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

	/// Maximum number of days in an LPJ-GUESS simulation year
	/** The standard version doesn't yet support leap years. */
	static const int MAX_YEAR_LENGTH = 365;

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

	/// The calendar year corresponding to simulation year 0
	int first_calendar_year;

private:

	int nyear;
	bool initialised;

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
		first_calendar_year = 0;
		nyear = 0;
		ismidday = false;
		initialised = false;
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
		initialised = true;
	}

	/// Call at end of every simulation day to update member variables.
	void next() {
		if (islastday) {
			if (islastmonth) {
				day = 0;
				month = 0;
				year++;
				if (year == nyear - 1) islastyear = true;
				islastmonth = false;
			}
			else {
				day++;
				month++;
				if (month == 11) islastmonth = true;
			}
			dayofmonth = 0;
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

	/// Sets calendar year for simulation year 0
	/** Astronomical year numbering is used, so year 1 BC is represented by 0,
	 *  2 BC = -1 etc. See ISO 8601.
	 */
	void set_first_calendar_year(int calendar_year) {
		first_calendar_year = calendar_year;
	}

	/// Returns the calendar year corresponding to the current simulation year
	/** Astronomical year numbering is used, so year 1 BC is represented by 0,
	 *  2 BC = -1 etc. See ISO 8601.
	 */
	int get_calendar_year() const {
		return year + first_calendar_year;
	}

	/// \returns the number of days in the current simulation year
	/** For this function to work properly in simulations with varying number
	 *  of days per year, the set_first_calendar_year must have been called first.
	 *
	 *  Currently there is no support for leap years, so this function
	 *  always returns 365. */
	int year_length() const {
		return MAX_YEAR_LENGTH;
	}

	/// Step n days from a date.
	/** Current implementation does not consider leap days, and the same 
	  * apply for the current use of the function through-out the model.
      */
	static int stepfromdate(int day, int step) {

		if(day < 0)			// a negative value should not be a valid day
			return -1;
		else if(day + step > 0)
			return (day + step) % MAX_YEAR_LENGTH;
		else if(day + step < 0)
			return day + step + MAX_YEAR_LENGTH;
		else
			return 0;
	}

	bool is_initialised() const {
		return initialised;
	}
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

/// Object updating gridcell mass balance; currently used in framework()
class MassBalance : public Serializable  {

	int start_year;

	double ccont;
	double ccont_zero;
	double ccont_zero_scaled;
	double ccont_yesterday;

	double cflux;
	double cflux_zero;
	double cflux_yesterday;

	double ncont;
	double ncont_zero;
	double ncont_zero_scaled;
	double ncont_yesterday;

	double nflux;
	double nflux_zero;
	double nflux_yesterday;

	double water_cont;
	double water_cont_zero;
	double water_cont_zero_scaled;
	double water_cont_yesterday;

	double water_flux;
	double water_flux_zero;
	double water_flux_yesterday;

public:
	MassBalance() {

		start_year = nyear_spinup;
		ccont = 0.0;
		ccont_zero = 0.0;
		ccont_zero_scaled = 0.0;
		cflux = 0.0;
		cflux_zero = 0.0;
		ncont = 0.0;
		ncont_zero = 0.0;
		ncont_zero_scaled = 0.0;
		nflux = 0.0;
		nflux_zero = 0.0;
		water_cont = 0.0;
		water_cont_zero = 0.0;
		water_cont_zero_scaled = 0.0;
		water_flux = 0.0;
		water_flux_zero = 0.0;
		ccont_yesterday = 0.0;
		cflux_yesterday = 0.0;
		ncont_yesterday = 0.0;
		nflux_yesterday = 0.0;
		water_cont_yesterday = 0.0;
		water_flux_yesterday = 0.0;
	}

	MassBalance(int start_yearX) {

		start_year = start_yearX;
		ccont = 0.0;
		ccont_zero = 0.0;
		ccont_zero_scaled = 0.0;
		cflux = 0.0;
		cflux_zero = 0.0;
		ncont = 0.0;
		ncont_zero = 0.0;
		ncont_zero_scaled = 0.0;
		nflux = 0.0;
		nflux_zero = 0.0;
		water_cont = 0.0;
		water_cont_zero = 0.0;
		water_cont_zero_scaled = 0.0;
		water_flux = 0.0;
		water_flux_zero = 0.0;
		ccont_yesterday = 0.0;
		cflux_yesterday = 0.0;
		ncont_yesterday = 0.0;
		nflux_yesterday = 0.0;
		water_cont_yesterday = 0.0;
		water_flux_yesterday = 0.0;
	}

	void init(Gridcell& gridcell);
	void check(Gridcell& gridcell);
	// indiv and patch-level functions are for use with true crop stands only
	void init_indiv(Individual& indiv);
	bool check_indiv(Individual& indiv, bool check_harvest = false);
	bool check_indiv_C(Individual& indiv, bool check_harvest = false);
	bool check_indiv_N(Individual& indiv, bool check_harvest = false);
	
	void init_patch(Patch& patch);
	bool check_patch(Patch& patch, bool check_harvest = false);
	bool check_patch_C(Patch& patch, bool check_harvest = false);
	bool check_patch_N(Patch& patch, bool check_harvest = false);
	bool check_patch_water(Patch& patch);

	void check_year(Gridcell& gridcell); // calls both check_year_C and check_year_N
	void check_year_C(Gridcell& gridcell);
	void check_year_N(Gridcell& gridcell);
	void check_year_water(Gridcell& gridcell);
	void check_period(Gridcell& gridcell);

	void check_day(Gridcell&);

	void serialize(ArchiveStream& arch);

private:
	void check_day_c(Gridcell&);
	void check_day_n(Gridcell&);
	void check_day_water(Gridcell&);
};

/// Struct containing wood harvest information from from forest stands created by LUC functionality
struct Wood_harvest_struct {

	double prim_frac;
	double prim_cmass;
	double sec_mature_frac;
	double sec_mature_cmass;
	double sec_young_frac;
	double sec_young_cmass;

	/// Constructs an empty Wood_harvest_struct struct
	Wood_harvest_struct() {
		zero();
	}

	void zero() {
		prim_frac = 0.0;
		prim_cmass = 0.0;
		sec_mature_frac = 0.0;
		sec_mature_cmass = 0.0;
		sec_young_frac = 0.0;
		sec_young_cmass = 0.0;
	}
};

struct CanopyLayer {
	double lai;
	double fipar;
};

class Canopy {
public:
	/*
	Called once per day to invalidate/clear yesterday's values without altering
	layer structure.
	*/
	void clear();

	/*
	Set LAI and fipar in the given layer.

	@param layer zero-based index of layer to set (0=top layer).
	@param lai leaf area index of this layer (m2 m-2 ground).
	@param fipar fraction of PAR incident to this layer (0-1).
	*/
	void set(size_t layer, double lai, double fipar);

	/*
	Access the i-th canopy layer.

	@param i zero-based index of layer to access (0=top layer).
	*/
	CanopyLayer & operator[] (size_t i);

	/*
	Read-only access to the i-th canopy layer.

	@param i zero-based index of layer to access (0=top layer).
	*/
	const CanopyLayer & operator[] (size_t i) const;

	/*
	Get the number of layers in the canopy.
	*/
	size_t size() const;
private:
	/*
	Ensure that layer i exists, and if not, add layers up to and including i.

	@param layer index (0=top).
	*/
	void verify_size(size_t i);

	/*
	The vector of canopy layers. First layer is top, last is bottom.
	*/
	std::vector<CanopyLayer> layers;
};

/// Inputs required to calculate photosynthesis for a single individual in a
/// single canopy layer at a single subdaily timestep (DAVE photosynthesis).
class PhotosynthesisInputs {
public:
	// todo: do we want to protect these members?
	// It seems like a sensible thing to do, but it makes things a bit uglier
	// from the caller's perspective (which seems to be the only thing we care
	// about).

	// Residual stomatal conductance term
	double gdark;

	// LAI of this individual in this layer.
	double lai;

	// Cumulative LAI above this layer.
	double lai_cum;

	// Physiological downregulation due to dormancy.
	double fphen_activity;

	// Photosynthetic limitation due to sink activity [0-1].
	double ps_sink_lim;

	// Maximum carboxylation rate at 25 degC (μmol m-2 s-1).
	double vcmax25;

	// Maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1).
	double vpmax25;

	// Maximum electron transport rate at 25 degC (μmol m-2 s-1).
	double jmax25;

	// Nitrogen scalar to account for the nitrogen gradient within the canopy.
	double nscal;

	// Water stress parameter (0-1 range; 1=minimum stress).
	double wscal;

	// Leaf temperature in °C.
	double tleaf;

	// activation energy for vcmax (kJ mol-1)
	double ha_v;

	// activation energy for jmax (kJ mol-1)
	double ha_j;

	// Entropy terms (deltaS) for vcmax (J mol-1 K-1)
	double ds_v;

	// Entropy terms (deltaS) for jmax (J mol-1 K-1)
	double ds_j;

	// Quantum yield of electron transport (mol e- mol-1 photon).
	double alpha;

	// Above-canopy radiation this timestep (W m-2).
	double par;

	// Fraction of PAR incident to this layer (0-1).
	double fipar;

	// Fraction of PAR absorbed by leaf area (0-1).
	double absorptance_par;

	// Fraction of NIR absorbed by leaf area (0-1).
	double absorptance_nir;

	// Leaf width (m).
	double leaf_width;

	// PAR incident to this layer (W m-2).
	double ipar() const {
		return par * fipar;
	}

	// Absorbed PAR in the middle of the layer (W m-2).
	double apar() const {
		return ipar() * absorptance_par;
	}

	// VPD in kPa.
	double vpd;

	// Atmospheric CO2 concentration in mol mol-1.
	double ca;

	// Medlyn stomatal term.
	double X;

	// Whether stomata are open (at all). If FALSE, stomata are fully closed.
	bool stomata_open;

	// Predawn leaf water potential (MPa)
	double predawn_lwp;

	// Sensitivity parameter indicating steepness of decline in Vcmax with water stress (Zhou 2014)
	double s_f;

	// Reference value at which Vcmax or Jmax is half it's non-water-stresed value (Zhou 2014)
	double psi_f;

	/// Tgrowth, mean air temperature of preceding 30 d.
	double tgrowth;

    /// The acclimation sensitivity coefficient for respiration (1/degC).
	double racclim;
};

/// Outputs from the photosynthesis function for a single individual in a single
/// canopy layer at a single subdaily timestep (DAVE photosynthesis).
class PhotosynthesisOutputs {
public:
	PhotosynthesisOutputs() {
		clear();
	}
	void clear() {
		rd = 0;
		an = 0;
		anc = 0;
		anj = 0;
		anp = 0;
		gsc = 0;
		ci = 0;
		vcmax = 0;
		jmax = 0;
		tleaf = 0;
	}

	// Respiration (mol m-2 s-1).
	double rd;

	// Net photosynthesis (mol m-2 s-1).
	double an;

	// Rubisco-limited net photosynthesis (mol m-2 s-1).
	double anc;

	// RuBP regeneration-limited photosynthesis (mol m-2 s-1).
	double anj;

	// PEP carboxylation-limited photosynthesis (mol m-2 s-1) [C3 only].
	double anp;

	// Stomatal conductance to CO2 (mol m-2 s-1).
	double gsc;

	// Intercellular CO2 concentration (mol mol-1).
	double ci;

	// Leaf temperature (degC).
	double tleaf;

	// Maximum carboxylation rate (mol m-2 s-1).
	double vcmax;

	// Maximum electron transport rate (mol m-2 s-1).
	double jmax;
};

class PhotosynthesisDiagnostics {
public:
	double g1;
	double g0;
	double gcutic_adj;
};

class EnergyBalance {
public:
	EnergyBalance(double tleaf, size_t niter, double delta_tleaf,
						bool converged, double transp, double rnet) :
		tleaf(tleaf), niter(niter), delta_tleaf(delta_tleaf),
		converged(converged), transp(transp), rnet(rnet) {
	}
	double get_tleaf() const { return tleaf; }
	size_t get_niter() const { return niter; }
	double get_delta_tleaf() const { return delta_tleaf; }
	bool get_converged() const { return converged; }
	double get_transp() const { return transp; }
	double get_rnet() const { return rnet; }
private:
	// Leaf temperature (degC).
	double tleaf;
	// Number of iterations taken to converge.
	size_t niter;
	// Difference in leaf temperature between last two iterations (degC).
	double delta_tleaf;
	// Did the solver converge within the maximum number of iterations?
	bool converged;
	// Transpiration rate (mol H2O m-2 s-1).
	double transp;
	// Net radiation (W m-2).
	double rnet;
};

/*
Inputs to the C4 photosynthesis solver function for a single layer/timestep.
*/
struct C4Biochemistry {
	// Scaled maximum Rubisco carboxylation rate (mol m-2 s-1).
	double vcmax;
	// Scaled maximum electron transport rate (mol m-2 s-1).
	double jmax;
	// Scaled maximum PEP carboxylation rate (mol m-2 s-1).
	double vpmax;
	// PEP regeneration rate (mol m-2 s-1).
	double vpr;
	// Daytime (mitochondrial) leaf respiration (mol m-2 s-1).
	double rd;
	// Rd in the mesophyll cells (mol m-2 s-1).
	double rm;
	// Rd in the bundle-sheath cells (mol m-2 s-1).
	double rs;
	// Michaelis-Menten constant of PEPC for CO2 (mol mol-1).
	double kp;
	// Michaelis-Menten constant of Rubisco for CO2 (mol mol-1).
	double kc;
	// Michaelis-Menten constant of Rubisco for O2 (mol mol-1).
	double ko;
	// Effective Michaelis-Menten constant of Rubisco for CO2 in the presence of O2 (mol mol-1).
	double km;
	// Ratio of solubility and diffusivity of O2 to CO2.
	double a0;
	// Half the reciprocal of Rubisco specificity.
	double lgamstar;
};

struct C4IterationState {
	// Intercellular CO2 concentration (mol mol-1).
	double ci;
	// PEP carboxylation rate (mol m-2 s-1).
	double vp;
};

/*
Outputs from the C4 photosynthesis solver function for a single layer/timestep.
*/
struct C4SolverResult {
	// Rubisco-limited net photosynthesis (mol m-2 s-1).
	double ac;
	// RuBP regeneration-limited photosynthesis (mol m-2 s-1).
	double aj;
	// Net photosynthesis rate (mol m-2 s-1).
	double an;
	// Stomatal conductance to CO2 (mol m-2 s-1).
	double gsc;
	// Intercellular CO2 concentration (mol mol-1).
	double ci;
	// Number of iterations taken to find a solution.
	int iterations;
	// Whether the solution converged within the maximum number of iterations.
	bool converged;
	// Difference between ci in the last two iterations (mol mol-1).
	double ci_diff;
};

/// Structure containing variables for sub-daily photosynthesis calculation.
struct PhotosynthesisSubdaily : public Serializable {
public:
	// Constructs an empty structure
	PhotosynthesisSubdaily() {
		// Ensure capacity is initialised before it's accessed for the 1st time.
		capacity = 0;
		resize(date.subdaily);
	}

	/*
	Resize the arrays to the specified size.
	@param size: Desired size of the arrays.
	*/
	void resize(int size) {
		if (capacity == size) {
			return;
		}

		capacity = size;

		clear();
	}

	// clears all variables
	void clear(){
		vcmax25   = 0.0;
		jmax25    = 0.0;
		vpmax25   = 0.0;

		vcmax.clear();
		jmax.clear();
		gamstar.clear();
		Rd.clear();
		An.clear();
		Anc.clear();
		Anj.clear();
		gsw.clear();
		ci.clear();
		for (int i = 0; i < capacity; i++) {
			vcmax.push_back(0.0); 
			jmax.push_back(0.0);
			gamstar.push_back(0.0);
			Rd.push_back(0.0);
			An.push_back(0.0);
			Anc.push_back(0.0);
			Anj.push_back(0.0);
			gsw.push_back(0.0);   
			ci.push_back(0.0);
		}
	}

	// Maximum carboxylation rate at 25 degC (μmol m-2 s-1)
	double vcmax25;

	// Maximum electron transport rate at 25 degC (μmol m-2 s-1)
	double jmax25;  

	// Maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1)
	double vpmax25;

	// Maximum carboxylation rate (mol m-2 s-1)
	std::vector<double> vcmax;

	// Maximum electron transport rate (mol m-2 s-1)
	std::vector<double> jmax;

	// Photorespiratory CO2 compensation point (mol mol-1)
	std::vector<double> gamstar;

	// Leaf day respiration (mol m-2 s-1)
    std::vector<double> Rd; 

	// Net photosynthesis (mol m-2 s-1)
	std::vector<double> An;

	// Rubisco-limited net photosynthesis (mol m-2 s-1)
	std::vector<double> Anc;

	// RuBP regeneration-limited photosynthesis (mol m-2 s-1)
	std::vector<double> Anj;

	// Stomatal conductance to water vapour (mol m-2 s-1)
	std::vector<double> gsw;

	// Intercellular CO2 concentration (mol mol-1)
	std::vector<double> ci;


    void serialize(ArchiveStream& arch);

private:
	int capacity;
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
		nmass_leaf_opt = 0.0;
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

	/// optimal leaf nitrogen (kg N m-2)
	double nmass_leaf_opt;

	/// nitrogen limitation on vm
	double vmaxnlim;

	/// Daily Vmax for the last 5 days [nmol g-1 d.w. s-1]
	Historic<double, 5> vmax_5;

	/// net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
    double net_assimilation() const {
		return (agd_g - rd_g) * KG_PER_G;
    }

	/// Gross daily C assimilation (kgC/m2/day)
	double gross_assimilation() const {
		return agd_g * KG_PER_G;
	}

	void serialize(ArchiveStream& arch);
};

/// Class containing serializable variables for Weathergenerator GWGen
/** 
Variables for the build-in random-generator and to keep track of whether past days 
were rain days
*/
class WeatherGenState : public Serializable {

public:
	/// Random state variable q
	int q[10];
	/// Random state variable carry
	int carry;
	/// Random state variable xcng
	int xcng;
	/// Random state variable xs
	unsigned int xs; 
	/// Random state variable indx
	int indx;
	/// Random state variable have
	bool have;
	/// Random state gamma	
	double gamma_vals[2];
	/// Indicator for whether the recent two days were rein-days
	bool pday[2];
	/// Random state's residuals
	double resid[4];

	WeatherGenState() {
		for (int i = 0; i<10; i++) q[i] = 0;
		carry = 0;
		xcng = 0;
		xs = 0;
		indx = 0;
		have = false;
		for (int i = 0; i<2; i++) gamma_vals[i] = 0.0;
		for (int i = 0; i<2; i++) pday[i] = false;
		for (int i = 0; i<4; i++) resid[i] = 0.0;
	}

	void serialize(ArchiveStream& arch);
};

/// This struct contains the environmental input to a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisEnvironment {

	/// Constructs an empty result
	PhotosynthesisEnvironment() {
		clear();
	}

	/// Clears all members
	/** Nonsense values to cause a crash if used
	 */
	void clear() {
		co2 = 0;
		temp = 0;
		par = 0;
		fpar = 0;
		daylength = 0;
	}

private:
	/// atmospheric ambient CO2 concentration (ppmv)
	double co2;

	/// mean air temperature today (deg C)
	double temp;

	/// total daily photosynthetically-active radiation today (J / m2 / day) (ALPHAA not yet accounted for)
	double par;

	/// fraction of PAR absorbed by foliage
	double fpar;

	/// day length, must equal 24 in diurnal mode(h)
	double daylength;

public:
	/// set 
	void set(double co2_env, double temp_env, double par_env, double fpar_env, double daylength_env) {
		co2 = co2_env;
		temp = temp_env;
		par = par_env;
		fpar = fpar_env;
		daylength = daylength_env;
	}

	double get_apar() const {
		return par * fpar;
	}

	double get_par() const {
		return par;
	}

	double get_fpar() const {
		return fpar;
	}

	double get_temp() const {
		return temp;
	}

	double get_co2() const {
		return co2;
	}

	double get_daylength() const {
		return daylength;
	}

};

/// This struct contains the stresses used in a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisStresses {

	/// Constructs an empty result
	PhotosynthesisStresses() {
		no_stress();
	}

	/// All members set to no stress values
	/** Default values indicating no stress
	 */
	void no_stress() {

		ifopt = true;
		iflimvmax = false;
		moss_ps_limit = 1.0;
		graminoid_ps_limit = 1.0;
		inund_stress = 1.0;
	}

private:
	/// whether to calculate optimal leaf nitrogen
	bool ifopt;

	/// whether nutrients should limit Vmax
	bool iflimvmax;

	///  limit to moss photosynthesis. [0,1], where 1 means no limit
	double moss_ps_limit;

	/// limit to graminoid photosynthesis. [0,1], where 1 means no limit
	double graminoid_ps_limit;

	/// limit to photosynthesis due to inundation, where 1 means no limit
	double inund_stress;

public:
	/// Set the stresses
    void set(bool thisopt, bool thislimvmax, double this_moss_ps_limit, double this_graminoid_ps_limit, double this_inund_stress) {

		ifopt = thisopt;
		iflimvmax = thislimvmax;
		moss_ps_limit = this_moss_ps_limit;
		graminoid_ps_limit = this_graminoid_ps_limit;
		inund_stress = this_inund_stress;
    }

	bool get_ifopt() const {
		return ifopt;
	}

	bool get_iflimvmax() const {
		return iflimvmax;
	}

	double get_moss_ps_limit() const {
		return moss_ps_limit;	
	}

	double get_graminoid_ps_limit() const {
		return graminoid_ps_limit;	
	}

	double get_inund_stress() const {
		return inund_stress;	
	}
};

/// This struct contains the leaf properties input to a photosynthesis calculation.
/** \see photosynthesis */
struct PhotosynthesisLeafstatus {

	/// Constructs an empty result
	PhotosynthesisLeafstatus() {
		noinfluence();
	}

	/// For PFTs use non sensible values as they don't have any influence
	/** For PFTs use non sensible values as they don't have any influence
	 */
	void noinfluence() {
		cmass_leaf = 0.0;
		nmass_leaf = 0.0;
		cton_min = 0.0;
		cton_max = 0.0;
		nextin = 1.0;
	}

private:
	/// Todays leaf C mass (kgC/m2)
	double cmass_leaf;

	/// Todays leaf N mass (kgN/m2)
	double nmass_leaf;

	/// Minimum leaf C:N ratio
	double cton_min;

	/// Maximum leaf C:N ratio
	double cton_max;

	/// Scalar depending on individual lai to slow down light optimization of newly shaded leaves
	double nextin;

public:
	/// set 
	void set(double cmass, double nmass, double ctonmin, double ctonmax, double nin) {
		cmass_leaf = cmass;
		nmass_leaf = nmass;
		cton_min = ctonmin;
		cton_max = ctonmax;
		nextin = nin;
	}

	double get_cmass_leaf() const {
		return cmass_leaf;
	}

	double get_nmass_leaf() const {
		return nmass_leaf;
	}

	double get_cton_min() const {
		return cton_min;
	}

	double get_cton_max() const {
		return cton_max;
	}

	double get_nextin() const {
		return nextin;
	}
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
	/// reference to parent Gridcell object
	Gridcell& gridcell;

	/// values for randomisation in Weathergenerator GWGEN
	WeatherGenState weathergenstate;

	/// mean air temperature today (deg C)
	double temp;

	/// total daily net downward shortwave solar radiation today (J/m2/day)
	double rad;

	/// total daily photosynthetically-active radiation today (J/m2/day)
	double par;

	/// precipitation today (mm)
	double prec;

	/// 10 m wind (m/s)
	double u10;

	/// rel. humidity (fract.)
	double relhum;

	/// min and max daily temperature (deg C)
	double tmin, tmax; 

	/// day length today (h)
	double daylength;

	/// maximum daylength in a year
	double daylength_max;

	/// minimum daylength in a year
	double daylength_min;

	/// atmospheric ambient CO2 concentration today (ppmv)
	double co2;

	/// latitude (degrees; +=north, -=south)
	double lat;

	/// gridcell mean elevation
	double mean_elevation;

	/// Insolation today, see also instype
	double insol;

	/// Type of insolation
	/** This decides how to interpret the variable insol,
	 *  see also documentation for the insoltype enum.
	 */
	insoltype instype;

	/// annual insolation counter
	double ainsol;

	/// equilibrium evapotranspiration today (mm/day)
	double eet;

	/// mean temperature for the last 31 days (deg C)
	double mtemp;

	/// mean of lowest mean monthly temperature for the last 20 years (deg C)
	double mtemp_min20;

	/// mean of highest mean monthly temperature for the last 20 years (deg C)
	double mtemp_max20;

	/// highest mean monthly temperature for the last 12 months (deg C)
	double mtemp_max;

	/// maximum par at a gridcell
	double par_max;

	/// relative par today (par/par_max)
	double rel_par;
	/// accumulated growing degree day sum on 5 degree base
	/** reset when temperatures fall below 5 deg C */
	double gdd5;

	/// total gdd5 (accumulated) for this year (reset 1 January)
	double agdd5;

	/// accumulated growing degree day sum on 0 degree base (Wolf et al. 2008)
	double gdd0;

	/// total gdd0 (accumulated) for this year (reset 1 January)
	double agdd0;

	/// total gdd0 (accumulated) over each of the last 20 years
	/// climate.agdd0_20.mean() gives the average total gdd0 (accumulated) over the last 20 years
	Historic<double, 20> agdd0_20;

	/// number of days with temperatures <5 deg C
	/** reset when temperatures fall below 5 deg C;
	 *  maximum value is number of days in the year */
	int chilldays;

	/// true if chill day count may be reset by temperature fall below 5 deg C
	bool ifsensechill;

	/** Respiration response to today's air temperature incorporating damping of Q10
	 *  due to temperature acclimation (Lloyd & Taylor 1994)
	 */
	double gtemp;

	/// daily temperatures for the last 31 days (deg C)
	Historic<double, 31> dtemp_31;

	/// daily precipitation for the last 31 days (deg C)
	Historic<double, 31> dprec_31;

	/// daily eet for the last 31 days (deg C)
	Historic<double, 31> deet_31;

	/// from Grote. indicates we are in the winter to summer transition, determined from daylength increases between days
	bool winter_to_summer;

	/// minimum monthly temperatures for the last 20 years (deg C)
	double mtemp_min_20[20];

	/// maximum monthly temperatures for the last 20 years (deg C)
	double mtemp_max_20[20];

	/// minimum monthly temperature for the last 12 months (deg C)
	double mtemp_min;

	/// mean of monthly temperatures for the last 12 months (deg C)
	double atemp_mean;

	/// first spinup years annual average air temperature (deg C)
	double temp_iaavr;

	// BLAZE
	/// average annual rainfall (mm/a)
	double rainfall_annual_avg;
	/// current sum of annual Rainfall (mm)
	double rainfall_cur;
	/// Accumulated last rainfall (mm)
	double last_rainfall;
	/// Days since last rainfall 
	double days_since_last_rainfall;
	/// Keetch-Byram-Drought-Index
	double kbdi;
	/// McArthur forest fire index (FFDI)
	double mcarthur_forest_fire_index;	
	/// To keep track of running months daily FFDI 
	double ffdi_monthly[30];	

	// Saved parameters used by function daylengthinsoleet

	double sinelat;
	double cosinelat;
	double qo[Date::MAX_YEAR_LENGTH];
	double u[Date::MAX_YEAR_LENGTH];
	double v[Date::MAX_YEAR_LENGTH];
	double hh[Date::MAX_YEAR_LENGTH];
	double sinehh[Date::MAX_YEAR_LENGTH];
	double daylength_save[Date::MAX_YEAR_LENGTH];
	/// indicates whether saved values exist for this day
	bool doneday[Date::MAX_YEAR_LENGTH];

	/// diurnal temperature range, used in daily/monthly BVOC (deg C)
	double dtr;

	// containers for sub-daily values of temperature, short-wave downward
	// radiation, par, rad and gtemp (equivalent to temp, insol, par, rad and gtemp)
	// NB: units of these variable are the same as their daily counterparts,
	// i.e. representing daily averages (e.g. pars [J/m2/day])

	/// Sub-daily temperature (deg C) (\see temp)
	std::vector<double> temps;

	/// Sub-daily insolation (\see insol)
	std::vector<double> insols;

	/// Sub-daily vapour pressure deficit (kPa).
	std::vector<double> vpds;

	/// Sub-daily PAR (J/m2/day) (\see par)
	std::vector<double> pars;

	/// Sub-daily net downward shortwave solar radiation (J/m2/day) (\see rad)
	std::vector<double> rads;

	/// Sub-daily respiration response (\see gtemp)
	std::vector<double> gtemps;

	/// Sub-daily precipitation (mm) (\see prec)
	std::vector<double> precs;

	/// Sub-daily atmospheric pressure values (Pa).
	std::vector<double> pressure;

	/// Sub-daily 10 m wind speed (m s-1) (\see u10)
	std::vector<double> winds;

	/// Variables used for crop sowing date or seasonality calculation

	/// max temperature during the last crop test period (between two summer seasons)
	double maxtemp;
	/// summer day when we test last year's crossing of sowing temperature limits; NH:June 30(day 180), SH:Dec.31(day 364), set in getgridcell()
	int testday_temp;
	/// last day of dry month when we test last year's crossing of sowing precipitation limits; NH:Dec.31(day 364), SH:June 30(day 180), set in getgridcell()
	int testday_prec;
	/// date used for sowing if no frost or spring occured during the year between the testmonths; NH:14, SH:195, set in getgridcell()
	int coldestday;
	/// used to adapt equations to hemisphere, set in getgridcell()
	int adjustlat;
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

	Historic<double, 20> hmtemp_20[12];
	Historic<double, 20> hmprec_20[12];
	Historic<double, 20> hmeet_20[12];

	/// seasonality type (SEASONALITY_NO, SEASONALITY_PREC, SEASONALITY_PRECTEMP, SEASONALITY_TEMP, SEASONALITY_TEMPPREC)
	seasonality_type seasonality;
	seasonality_type seasonality_lastyear;

	/// precipitation seasonality type (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the extremes of the 20-year monthly means
	 */
	prec_seasonality_type prec_seasonality;
	prec_seasonality_type prec_seasonality_lastyear;

	/// precipitation range (DRY, DRY_INTERMEDIATE, DRY_WET, INTERMEDIATE, INTERMEDIATE_WET, WET)
	/** based on the average of the 20-year monthly extremes
	 */
	prec_seasonality_type prec_range;
	prec_seasonality_type prec_range_lastyear;

	/// temperature seasonality (COLD, COLD_WARM, COLD_HOT, WARM, WARM_HOT, HOT)
	temp_seasonality_type temp_seasonality;
	temp_seasonality_type temp_seasonality_lastyear;

	/// whether several months with precipitation maxima exists (remains to be implemented)
	bool biseasonal;

	/// variation coefficient of 20-year mean monthly temperatures
	double var_prec;
	/// variation coefficient of 20-year mean monthly precipitation to PET ratios
	double var_temp;

	/// annual precipitation sum
	double aprec;
	/// annual average precipitation (last year) (mm)
	double aprec_lastyear;	

	/// annual PAR sum
	double apar;

	/// Current (calendar) year of met data from input file.
	int met_year;

	/// Long-term mean maximum temperature of the warmest month (°C).
	double thome;

	/*
	Get the growth temperature for the current day, calculated as the mean of the
	previous 30 days' mean daily air temperatures.
	*/
	double get_tgrowth() const;

	/*
	Get mean daily respiration response to today's soil temperature at 0.25 m
	depth.
	*/
	double get_gtemp_mean() const;
public:
	/// constructor function: initialises gridcell member
	Climate(Gridcell& gc):gridcell(gc) {

		var_prec = 0.0;
		var_temp = 0.0;
		apar = 0.0;
		aprec = 0.0;
		aprec_lastyear = 0.0;
		met_year = 0;

		winter_to_summer = true;

		par_max	= 0.0;
		rel_par  = 0.0;

		for(int m=0;m<12;m++) {

			mtemp20[m] = 0.0;
			mprec20[m] = 0.0;
			mpet20[m] = 0.0;
			mpet_year[m] = 0.0;
			mprec_pet20[m] = 0.0;

			for(int y=0;y<20;y++) {
				mtemp_20[y][m] = 0.0;
				mprec_20[y][m] = 0.0;
				mpet_20[y][m] = 0.0;
				mprec_pet_20[y][m] = 0.0;
			}
		}

		for(int y=0;y<20;y++) {

			mprec_petmin_20[y] = 0.0;
			mprec_petmax_20[y] = 0.0;
		}

		mprec_petmin20=0.0;
		mprec_petmax20=0.0;

		seasonality=SEASONALITY_NO;
		seasonality_lastyear=SEASONALITY_NO;
		prec_seasonality=DRY;
		prec_seasonality_lastyear=DRY;
		prec_range=DRY;
		prec_range_lastyear=DRY;
		temp_seasonality=COLD;
		temp_seasonality_lastyear=COLD;
		biseasonal=false;
		mean_elevation = 0;
		u10 = 0;

		eet=0.0;
		par = 0.0;

		gdd0 = 0.0;
		agdd0 = 0.0;

		chilldays = 0;
		daylength = 0;
		daylength_max = -9999;
		daylength_min = 9999;
	};

	/// Initialises certain member variables
	/** Should be called before Climate object is applied to a new grid cell */
	void initdrivers(double latitude) {

		std::fill_n(mtemp_min_20, 20, 0.0);
		std::fill_n(mtemp_max_20, 20, 0.0);

		mtemp_min20 = 0.0;
		mtemp_max20 = 0.0;
		mtemp = 0.0;
		maxtemp = 0.0;
		gdd5 = 0.0;
		chilldays = 0;
		ifsensechill = true;
		atemp_mean = 0.0;
		temp_iaavr = 0.0;

		lat = latitude;
		std::fill_n(doneday, Date::MAX_YEAR_LENGTH, false);
		sinelat = sin(lat * DEGTORAD);
		cosinelat = cos(lat * DEGTORAD);

		winter_to_summer = (lat > 0);

		// Set crop-specific members
		if (latitude >= 0) {
			testday_temp = 180;		//June 30(day 180)
			testday_prec = 364;		//Dec.31(day 364)
			coldestday = COLDEST_DAY_NHEMISPHERE;
			adjustlat = 0;
		}
		else {
			testday_temp = 364;		//Dec.31(day 364)
			testday_prec = 180;		//June 30(day 180)
			coldestday = COLDEST_DAY_SHEMISPHERE;
			adjustlat = 181;
		}

		// BLAZE related variables
		rainfall_annual_avg      = 0.0; // average annual rainfall [mm]
		days_since_last_rainfall = 0;   // #Days-since-last-rainfall >3mm
		last_rainfall            = 0.0; // rainfall of last day of previous year [mm]
		kbdi                     = 0.0; // Keetch-Byram-Drought-index []

		for (int x=0; x<30; x++) {
			ffdi_monthly[x] = 0.;
		}

	}

	void serialize(ArchiveStream& arch);
};

/*
Defines the period of a flux.
*/
enum class Period {
	DAY,
	MONTH,
	YEAR
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
		/// Flux to runoff associated with organic leaching (kgC/m2)
		LEACHORGC,
		/// Flux from atmosphere to vegetation associated with manure addition (kgC/m2)
		MANUREC,
		/// Flux to soil associated with biological nitrogen fixation (kgN/m2)
		NFIX,
		/// Flux to vegetation associated with manure addition (kgN/m2)
		MANUREN,
		/// Flux to vegetation associated with N addition (kgN/m2) 
		NFERT,
		/// Nitrogen flux from atmosphere to vegetation associated with establishment (kgN/m2)
		ESTN,
		/// Nitrogen flux to atmosphere from consumed harvested products (kgN/m2)
		HARVESTN,
		/// Nitrogen flux from atmosphere to vegetation associated with sowing (kgN/m2)
		SEEDN,
		/// Nitrogen flux to runoff associated with organic leaching (kgN/m2)
		LEACHORGN,
		/// Nitrogen flux to runoff associated with mineral leaching (kgN/m2)
		LEACHMINN,
		/// NH3 flux to atmosphere from fire
		NH3_FIRE,
		/// NOx flux to atmosphere from fire
		NOx_FIRE,
		/// N2O flux to atmosphere from fire
		N2O_FIRE,
		/// N2 flux to atmosphere from fire
		N2_FIRE,
		/// NH3 flux from soil (ntransform)
		NH3_SOIL,
		/// NO flux from soil (ntransform)
		NO_SOIL,
		/// N2O flux in soil (ntransform)
		N2O_SOIL,
		/// N2 flux from soil (ntransform)
		N2_SOIL,
		/// DOC flux from soil (ntransform)
		DOC_FLUX,
		/// Net nitrification (ntransform)
		NET_NITRIF,
		/// Net denitrification (ntransform)
		NET_DENITRIF,
		/// Gross nitrification (ntransform)
		GROSS_NITRIF,
		/// Gross denitrification (ntransform)
		GROSS_DENITRIF,
		/// Reproduction costs
		REPRC,
		/// Total (i.e. CH4C_DIFF + CH4C_PLAN + CH4C_EBUL) CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C,
		/// Diffused CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C_DIFF,
		/// Plant-mediated CH4 flux to atmosphere from peatland soils (gC/m2).
		CH4C_PLAN,
		/// CH4 flux to atmosphere from peatland soils due to ebullition (gC/m2).
		CH4C_EBUL,
		/// Surface runoff
		WRUNOFF_SURF,
		/// Surface standing water runoff
		WRUNOFF_SSW,
		/// Drainage runoff
		WRUNOFF_DRAIN,
		/// Baseflow runoff
		WRUNOFF_BASE,
		/// Soil evaporation
		WEVAP,
		/// Potential evapotranspiration
		WPET,
		/// Transpiration
		WTRAN,
		/// Interception
		WINTERCEP,
		/// Irrigation water
		WIRRI,
		/// Water added to wetlands
		WADDED,
		/// Number of types, must be last
		NPERPATCHFLUXTYPES
	};

	/// Fluxes stored per pft
	enum PerPFTFluxType {
		/// NPP (kgC/m2)
		NPP,
		/// GPP (kgC/m2) with nutrient and water stress
		GPP,
		/// GPP (kgC/m2) with no stress
		GPPno,
		/// GPP (kgC/m2) with nutrient stress
		GPPns,
		/// Autotrophic respiration (kgC/m2)
		RA,
		/// Isoprene (mgC/m2)
		ISO,
		/// Monoterpene (mgC/m2)
		MT_APIN,
		MT_BPIN,
		MT_LIMO,
		MT_MYRC,
		MT_SABI,
		MT_CAMP,
		MT_TRIC,
		MT_TBOC,
		MT_OTHR,
		/// Number of types, must be last
		NPERPFTFLUXTYPES
	};

	// emission ratios from fire (NH3, NOx, N2O, N2) Delmas et al. 1995
	// values in .cpp file

	static const double NH3_FIRERATIO;
	static const double NOx_FIRERATIO;
	static const double N2O_FIRERATIO;
	static const double N2_FIRERATIO;

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

	/// \returns daily flux for a given flux type and PFT
	double get_daily_flux(PerPFTFluxType flux_type, int day, int pft_id) const;

	/// \returns daily flux for a given flux type (for all PFTs)
	double get_daily_flux(PerPFTFluxType flux_type, int day) const;

	/// \returns daily flux for a given flux type
	double get_daily_flux(PerPatchFluxType flux_type, int day) const;

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

	/// \returns flux for a given period and flux type (for all PFTs)
	double get_flux(PerPatchFluxType flux_type, Period period) const;

	/// \returns flux for a given period and flux type
	double get_flux(PerPFTFluxType flux_type, Period period) const;
private:

	size_t daily_flux_index(int day, int pft_id, PerPFTFluxType flux_type) const;

	/// Stores one flux value per PFT and flux type
	std::vector<std::vector<double> > annual_fluxes_per_pft;

	/// Stores one flux value per month and flux type
	/** For the fluxes only stored as totals for the whole patch */
	double monthly_fluxes_patch[12][NPERPATCHFLUXTYPES];

	/// Stores one flux value per month and flux type
	/** For the fluxes stored per pft for annual values */
	double monthly_fluxes_pft[12][NPERPFTFLUXTYPES];

	/// Stores one flux value per day, flux type and PFT
	//  Flattened order: [day][pft][type]
	std::vector<double> daily_fluxes_per_pft;

	/// Stores one flux value per day and flux type
	double daily_fluxes_patch[365][NPERPATCHFLUXTYPES];

	/// Stores one flux value per day and flux type
	double daily_fluxes_pft[365][NPERPFTFLUXTYPES];
};

/// Storage class of crop management information for one rotation period for a stand type, read from the instruction file.
class ManagementType {

public:
	/// id code (should be zero based and sequential, 0...nst-1)
	int id;
	/// name of management type
	xtring name;
	/// First management year (calendar year)
	/** Sets time when common features for managed stands begin, e.g. relaxed establishment rules and absence of disturbance before harvest begins.
	 *  This currently only applies for stands with wood havest */
	int firstmanageyear;
	/// First year with wood harvest (calendar year)
	int firstcutyear;
	// Whether the reference age for thinning timing is time since firstcutyear rather than the pach age
	bool firstcutyear_is_referenceyear;
	/// First year with clearcut (calendar year)
	int firstclearcutyear;
	/// Number of years to distribute clearcut of patches that were due to be cut before firstclearcutyear (using ifclearcut_by_density)
	int delayduecutting;
	/// When to start cutting to reach pft target fractions (calendar year)
	int firsttargetyear;
	/// When to stop cutting to reach pft target fractions (calendar year)
	int lasttargetyear;
	/// type of planting system ("", "MONOCULTURE", "SELECTION", etc.)
	xtring planting_system;
	/// type of harvest system ("", "CLEARCUT", "CONTINUOUS")
	xtring harvest_system;
	/// name of crop pft 
	xtring pftname;
	/// pft planting density (seedlings/ha)
	double plantdensity_pft;
	/// string of pft selection
	xtring selection;
	/// string of pft planting densities (seedlings/ha)
	xtring plantdensity;
	/// string of pft target cmass fractions
	xtring targetfrac;
	/// pft target cmass fraction input file name
	xtring file_targetfrac_pft_mt;
	/// How to use values in file_targetfrac_pft_mt
	/*  (0): use values for target cutting and keep the selection defined in the instruction file, 
	 *  (1): use values for target cutting and restrict selection to pft:s in input file with non-zero fractions (overriding selection defined in the instruction file),
	 *  (2): only use target fractions to define selection (can be 0/1 values) as in (1), without using target cutting.
	 */
	int targetfrac_input_mode;
	/// Patch age when target cutting starts
	int targetstartage;
	/// Interval of target cuttings
	int targetcutinterval;
	/// Whether patch (1) or stand (2,3) deviations from pft cmass fraction targets used (see documentation for details)
	int targetcutmode;
	/// Whether young (1) or old (2) individuals are preferentially cut, or no preference (0) in target cuttings; overridden by thinselectdiam[] settings
	int	targetthinselectage;
	/// Whether small (1) or large (2) diameter individuals are preferentially cut, or no preference (0) in target cuttings.
	int	targetthinselectdiam;
	/// Whether to stop cutting to reach pft fraction targets when second (continuous) period starts
	bool suppress_second_target;
	/// Rotation time in years
	int cutinterval;
	/// Whether to use Reineke's rule-based automatic thinning
	bool ifthin_reineke;
	/// Self-thinning parameter for ifthin_reineke
	double alpha_st;
	/// Thinning "intensity" (low value more intense) when using ifthin_reineke
	double rdi_target;
	/// Whether to use tree density as a trigger for clearcut
	bool ifclearcut_by_density;
	/// Tree density target (trees/ha) below which a clearcut will occur for clearcut_by_density
	int dens_target_cc;
	/// Whether to use optimum rotation age as a trigger for clearcut
	bool ifclearcut_optimal_age;
	/// Whether to distribute patch ages in a new managed forest stand
	bool distribute_patch_ages;
	/// Minimum tree diameter limit (cm) for cutting (thinstrength*100)% of trees
	double diam_cut_low;
	/// Minimum tree diameter limit (cm) for cutting 100% of trees
	double diam_cut_high;
	/// Whether to adapt diam_cut_low to forest stands with small trees
	bool adapt_diam_limit;
	/// Timing of thinning events, relative to rotation period
	double thintime[NTHINNINGLOOPS][NTHINNINGS];
	/// Strength (fraction cut) of thinning events
	double thinstrength[NTHINNINGLOOPS][NTHINNINGS];
	/// Strength (fraction cut) of thinning events for unselected pft:s
	double thinstrength_unsel[NTHINNINGLOOPS][NTHINNINGS];
	/// Whether non-selected (1) or selected (2) pft:s are preferentially cut, unselected and selected cutting strengths specified separately (3), ...
	/** shrubs and shade-intolerant species preferentially cut (4) or no preference (0) */
	int thinselectpft[NTHINNINGLOOPS][NTHINNINGS];
	/// Whether young (1) or old (2) individuals are preferentially cut, or no preference (0); overridden by thinselectdiam[] settings
	int thinselectage[NTHINNINGLOOPS][NTHINNINGS];
	/// Whether small (1) or large (2) diameter individuals are preferentially cut, thinstrength of trees with diam between diam_cut_low and diam_cut_high...
	/** and 100% of trees with diam > diam_cut_high (3), or no preference (0) */
	int thinselectdiam[NTHINNINGLOOPS][NTHINNINGS];
	/// When to start contiuous cutting period (years after start of regeneration period)
	int secondintervalstart;
	/// Wood cutting interval in years in the contiuous cutting period
	int secondcutinterval;
	/// Whether to distribute cuttings evenly in time among patches in a managed forest stand
	bool distribute_cuttings_among_patches;
	/// Harvested forest area fraction, unless woodharv_frac read from file
	double woodharv_frac;
	/// Wood harvest C mass for gridcell, unless woodharv_cmass read from file
	double woodharv_cmass;
	/// Whether to ignore climate establishment limits
	bool relaxed_establishment;
	/// Whether to suppress fires
	bool suppress_fire;
	/// Whether to suppress disturbance
	bool suppress_disturbance;
	/// Whether to use tree pft planting densities after clearcut (plus bypass all envirinmental establishment limits)
	bool set_planting_density;
	/// Whether to clearcut firstmanageyear or rotation year (or first stand year); 0 = don't cut(clone), 1 = cut(don't clone), 2 = cut(clone)
	int cutfirstyear;
	/// Whether to cut pft:s outside of selection clone year or first year of new management in a rotation (if reestab "restricted" or "none")
	bool cutfirstyear_unsel;
	/// Whether to kill grass at clearcut
	bool killgrass_at_cc;
	/// Whether to use stochastic mortality
	bool stochmort;
	/// Whether to use stochastic establishment
	bool stochestab;
	/// Harvest efficiency during thinning for all tree pft:s
	double harv_eff_thin;
	/// Residue removal fraction during thinning for twigs and branches for all tree pft:s
	double res_outtake_twig_thin;
	/// Residue removal fraction during thinning for coarse roots and stumps for all tree pft:s
	double res_outtake_coarse_root_thin;
	/// Harvest efficiency during clearcut for all tree pft:s
	double harv_eff_cc;
	/// Residue removal fraction during clearcut for twigs and branches for all tree pft:s
	double res_outtake_twig_cc;
	/// Residue removal fraction during clearcut for coarse roots and stumps for all tree pft:s
	double res_outtake_coarse_root_cc;

	/// Water use (RAINFED,IRRIGATED) 
	waterusetype wateruse;
	/// irrigation efficiency
//	double firr;
	/// forced sowing date, unless sdate_force read from file
	int sdate;
	/// forced harvest date, unless hdate_force read from file
	int hdate;
	/// Nitrogen fertilisation amount, unless Nfert_read read from file
	double nfert;
	/// Tillage factor (1.0 - 4.5)
	double tillage_fact;
	/// Whether grass is grown in fallow
	bool fallow;
	/// Grazing intensity (daily grazing fraction)
	double grazeintens;
	/// Mowing date (julian day)
	int mowdates[num_mowdates];
	/// Mowing intensity (fraction of biomass removed on mowing date)
	double mowintens;

	ManagementType() {

		id = -1;
		firstmanageyear = FAR_FUTURE_YEAR;
		firstcutyear = FAR_FUTURE_YEAR;
		firstcutyear_is_referenceyear = false;
		firstclearcutyear = FAR_PREHISTORIC_YEAR;
		delayduecutting = 0;
		firsttargetyear = FAR_FUTURE_YEAR;
		lasttargetyear = FAR_FUTURE_YEAR;
		planting_system = "";
		harvest_system = "";
		pftname = "";
		plantdensity_pft = -1;
		selection = "";
		plantdensity = "";
		targetfrac = "";
		file_targetfrac_pft_mt = "";
		targetfrac_input_mode = 0;
		targetstartage = 10;
		targetcutinterval = 5;
		targetcutmode = 1;
		targetthinselectage = SELECT_AGE_NOPREF;
		targetthinselectdiam = SELECT_DIAM_LARGE;
		suppress_second_target = false;
		cutinterval = 0;
		ifthin_reineke = false;
		alpha_st = 0;
		rdi_target = 0;
		ifclearcut_by_density = false;
		dens_target_cc = 0;
		ifclearcut_optimal_age = false;
		distribute_patch_ages = false;
		for(int n=0;n<NTHINNINGLOOPS;n++) {
			for(int t=0;t<NTHINNINGS;t++) {
				thintime[n][t] = 0.0;
				thinstrength[n][t] = 0.0;
				thinstrength_unsel[n][t] = 0.0;
				thinselectpft[n][t] = SELECT_PFT_NOPREF;
				thinselectage[n][t] = SELECT_AGE_NOPREF;
				thinselectdiam[n][t] = SELECT_DIAM_NOPREF;
			}
		}
		diam_cut_low = 0.0;
		diam_cut_high = 1000.0;
		adapt_diam_limit = false;
		secondintervalstart = -1;
		secondcutinterval = 0;
		distribute_cuttings_among_patches = false;
		woodharv_frac = -1.0;
		woodharv_cmass = -1.0;
		relaxed_establishment = false;
		suppress_fire = false;
		suppress_disturbance = false;
		set_planting_density = false;
		cutfirstyear = 1;
		cutfirstyear_unsel = false;
		killgrass_at_cc = false;
		stochmort = true;
		stochestab = true;
		harv_eff_thin = -1.0;
		res_outtake_twig_thin = -1.0;
		res_outtake_coarse_root_thin = -1.0;
		harv_eff_cc = -1.0;
		res_outtake_twig_cc = -1.0;
		res_outtake_coarse_root_cc = -1.0;

		wateruse = RAINFED;
//		firr = 0.0;
		sdate = -1;
		hdate = -1;
		nfert = -1.0;
		tillage_fact = 1.0;
		fallow = false;
		grazeintens = 0.0;
		for (int n=0;n<num_mowdates;n++) {
			mowdates[n] = -1;
		}
		mowintens = 0.0;
	}

	bool is_managed() {

		// Add new management parameters here
		if(pftname != "" || planting_system != "" || selection != ""||  harvest_system != "" ||  wateruse == IRRIGATED 
			|| fallow || relaxed_establishment || suppress_fire || suppress_disturbance || nfert > -1.0 || firstmanageyear < FAR_FUTURE_YEAR)
			return true;
		else
			return false;
	}

	/// Returns true if pft is in pftselection.
	bool pftinselection(const char* name);
};

/// A list of management types
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  management types objects.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential management type objects can be referenced as array elements by id:
 *
 *   ManagementTypelist mtlist;
 *   ...
 *   for (i=0; i<nst; i++) {
 *     ManagementType& thismt=stlist[i];
 *     // query or modify object thismt here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   mtlist.firstobj();
 *   while (mtlist.isobj) {
 *    ManagementType& thismt=mtlist.getobj();
 *     // query or modify object thismt here
 *     mtlist.nextobj();
 *   }
 */
class ManagementTypelist : public ListArray_id<ManagementType> {

public:
	int getmtid(xtring mtname) {

		int id = -1;

		for(unsigned int i=0; i< this->nobj; i++) {

			ManagementType& mt = (*this)[i];
			if(mt.name == mtname) {
				id = mt.id;
				break;
			}
		}

		return id;
	}
};

/// The one and only linked list of ManagementType objects	
extern ManagementTypelist mtlist;

/// Storage class of management rotation information for a stand type, read from the instruction file.
struct Rotation {

	/// Number of managements in rotation
	int nmanagements;
	/// First crop rotation year
	int firstrotyear;
	/// Double cropping of one crop (e.g. rice)
	bool multicrop;

	Rotation() {
		nmanagements = 0;
		firstrotyear = 0;
		multicrop = false;
	}
};

/// Stand type class for storing both static parameters, read from the instruction file,
/*	and dynamic variables, updated in landcover_change()
 *  Active stand types are stored in the stlist analogous to the pftlist.
 */
class StandType {

public:
	/// id code (should be zero based and sequential, 0...nst-1)
	int id;
	/// name of stand type
	xtring name;

	/// specifies type of landcover
	/** \see landcovertype */
	landcovertype landcover;	// specifies type of landcover (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
	/// Rotation information, read from the instruction file
	Rotation rotation;
	/// Management struct (static); temporary storage of mt when defined directly in the instruction file stand type; copied to mt in mtlist
	ManagementType management;
	/// Management types in a rotation cycle
	xtring mtnames[NROTATIONPERIODS_MAX];
	/// Start of the managements in a rotation cycle (calendar year)
	int mtstartyear[NROTATIONPERIODS_MAX];
	/// Whether to wait for clearcut before moving to next mt in a forestry rotation
	bool rotation_wait_for_cc;
	/// Disturbance interval (years)
	double distinterval;

	/// intercrop (NOINTERCROP,NATURALGRASS)
	intercroptype intercrop;
	/// whether natural pft:s are allowed to grow in stand type
	xtring naturalveg; // "", "GRASSONLY", "ALL"
	/// whether planted pft:s or all active pft:s are allowed to established after planting in a forest stand ("", "RESTRICTED", "ALL")
	xtring reestab;

	/// used for output from separate stand types
	double anpp;
	double cmass;
	double cmass_tree;
	double cmass_tree_mort;
	double cmass_wood;
	double cmass_wood_potharv;
	double cmass_wood_potharv_products;
	double cmass_potfuel;
	double cmass_mort;
	double cmass_fire;
	double cmass_dist;
	double cmass_leaf_root_turnover;
	double cmass_repr;
	double cmass_est;
	double cmass_wood_harv;
	double cmass_wood_harv_toprod;
	double cmass_wood_clearcut;
	double cmass_harv_tolitter;
	double cmass_harv_killed;
	double densindiv;
	double diam_g;
	double csoil;
	double clitter;
	double csink;
	double lai;
	double lai_tree;
	double height;

	StandType() {

		intercrop = NOINTERCROP;
		naturalveg = "";
		reestab = "ALL";
		rotation_wait_for_cc = false;
		distinterval = 1.0e10;
		for(int m=0;m<NROTATIONPERIODS_MAX;m++)
			mtstartyear[m] = -1;
	}

	ManagementType& get_management(int rot = 0) {

		int mtid = mtlist.getmtid(mtnames[rot]);

		if(mtid < 0) {
			fail("get_management(): mt name not in mtlist; check instruction file.\n");
		}
		return mtlist[mtid];
	}

	/// Returns position of management in rotation list if present. Returns -1 if not.
	int mtinrotation(xtring name) {

		int mtno = -1;
		for(int i=0; i<rotation.nmanagements; i++) {
			if(name == mtnames[i])
				mtno = i;
		}

		return mtno;
	}

	/// Returns position of crop in rotation list if present. Returns -1 if not.
	int pftinrotation(xtring name) {

		int cropno = -1;
		for(int i=0; i<rotation.nmanagements; i++) {
			if(name == get_management(i).pftname)
				cropno = i;
		}

		return cropno;
	}
};

/// A list of stand types
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  stand types objects.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential stand type objects can be referenced as array elements by id:
 *
 *   StandTypelist stlist;
 *   ...
 *   for (i=0; i<nst; i++) {
 *     StandType& thisst=stlist[i];
 *     // query or modify object thisst here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   stlist.firstobj();
 *   while (stlist.isobj) {
 *     StandType& thisst=stlist.getobj();
 *     // query or modify object thisst here
 *     stlist.nextobj();
 *   }
 */
class StandTypelist : public ListArray_id<StandType> {

};

/// The one and only linked list of StandType objects
extern StandTypelist stlist;

/// TODO : Assaf copied and pasted this section as is from Grote branch (20/11/2023).
/// senescence helper to be able to have a within year senescene of leaves and roots
// is senescence pft specific oder individual specific?
// the senescence period should be maybe pft specicif, but the onset and end of senescecne
// would be indiv specific, e.g. TrBr indiv might start senescence based on their own wstress at some point

/// Bit coded senescence type:
// !bit: 0(1)=background allyear,1(2)=aphen<aphen_max,2(4)=daylen<39300s,
//       4(8)=TrBr drought (phen drop) w/ switch to background,5(16)=day==(coldestday-senlen)
//       6(32)=sense chill event; 7(64)=cmass excess handling for (TrBR, will modify background handling)
//       8(128=0x80)=SEN_TRBR_WSCAL_B
//       9(256=0x100)=SEN_APHEN_LO_MAX
//       10(512=0x200)=SEN_APHEN_MAX_GILL (currently not implemented)
//       11(1024=0x400)=SEN_WSCALMIN
//			 12(????=0x800)=SEN_DELPIERRE (not implemented)
enum sen_type_ { SEN_BACKGROUND, SEN_APHEN_MAX, SEN_DAYLENGTH, SEN_TRBR_WSCAL, SEN_COLDESTDAY, SEN_CHILLSENSE, SEN_CMASSEXCESS, SEN_TRBR_WSCAL_B, SEN_APHEN_LO_MAX, SEN_APHEN_MAX_GILL, SEN_WSCALMIN, NSEN_TYPE};
typedef std::bitset<NSEN_TYPE> sen_type;
//inline enum sen_type operator|(enum sen_type a, enum sen_type b)
//{
//    return static_cast<enum sen_type>(static_cast<int>(a) | static_cast<int>(b));
//}
//inline enum sen_type operator&(enum sen_type a, enum sen_type b)
//{
//    return static_cast<enum sen_type>(static_cast<int>(a) & static_cast<int>(b));
//}

/*
Calculate the root distribution given a root beta value.
*/
void calc_root_distribution(const Soil &soil, double root_beta, double rootdist[NSOILLAYERS_MAX]);

/// Holds static functional parameters for a plant functional type (PFT).
/** There should be one Pft object for each potentially occurring PFT. The same Pft object
 *  may be referenced (via the pft member of the Individual object; see below) by different
 *  average individuals. Member functions are included for initialising SLA given leaf
 *  longevity, and for initialising sapling/regen characteristics (required for
 * population mode).
 */
class Pft {

	// MEMBER VARIABLES

public:
	/// id code (should be zero based and sequential, 0...npft-1)
	int id;
	/// name of PFT
	xtring name;
	/// life form (tree or grass)
	lifeformtype lifeform;
	/// leaf phenology (raingreen, summergreen, evergreen, rain+summergreen, cropgreen)
	phenologytype phenology;
	/// leaf physiognomy (needleleaf, broadleaf)
	leafphysiognomytype leafphysiognomy;
	/// growing degree sum on 5 degree base required for full leaf cover
	double phengdd5ramp;
	/// growing degree sum on 5 degree base required for full maturity/senescence (grasses only)
	double phengdd5mat;
	/// water stress threshold for leaf abscission (range 0-1; raingreen PFTs)
	double wscal_min;
	/// water stress threshold when plants are dormant (grasses only)
	double wscal_min_dormant;
	/// Maximum physiological up/down-regulation in and out of dormancy (fraction d-1)
    double dphen_activity_max;
	/// biochemical pathway for photosynthesis (C3 or C4)
	pathwaytype pathway;
	/// approximate low temperature limit for photosynthesis (deg C)
	double pstemp_min;
	/// approximate lower range of temperature optimum for photosynthesis (deg C)
	double pstemp_low;
	/// approximate upper range of temperature optimum for photosynthesis (deg C)
	double pstemp_high;
	/// maximum temperature limit for photosynthesis (deg C)
	double pstemp_max;
	/// non-water-stressed ratio of intercellular to ambient CO2 partial pressure
	double lambda_max;
	/// vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
	/// Should be private, but parameters.cpp wants a pointer for parameter binding..
	/// This should generally never be used. Always use indiv.rootdist if possible.
	double input_rootdist[NSOILLAYERS_MAX];

	/// vegetation root profile input in an array containing fraction of roots in each soil layer, [0=upper layer]
	double root_dist[NSOILLAYERS_SIMPLE];
	/// shape parameter for initialisation of root distribtion
	double root_beta;
	/// canopy conductance component not associated with photosynthesis (mm/s)
	double gmin;
	/// Scalar for Vpmax25 (C4 photosynthesis only)
	double vpmax_scal;
    /// Stomatal intercept parameter (Medlyn et al. 2011) (mol m-2 s-1)
    double g0;
	/// Cuticular conductance (mol m-2 s-1)
	double gcutic;
    /// Stomatal slope parameter (Medlyn et al. 2011) (kPa^0.5)
    double g1;
	/// Stomatal sensitivity to water stress (0-1)
	double qs;
	/// Stomatal intercept parameter at leaf water potential = -0.3 Mpa (Zhou et al. 2014), (kPa^0.5)
    double g1x;
	/// maximum evapotranspiration rate (mm/day)
	double emax;
	/// maintenance respiration coefficient (0-1)
	double respcoeff;
	/// water potential at 12% loss of hydraulic conductance
	double stom_clos;
	/// water potential at 50% loss of hydraulic conductance
	double p50;
	/// water potential at 88% loss of hydraulic conductance
	double p88;
	/// shape of the vulnerability curve (unitless)
	double s2;
	/// sensitivity of the vulnerability curve (MPa)
	double s1;
	/// maximum hydraulic conductance per unit xylem area (mol/m2/s/MPa)
	double kmax;
	/// Nitrogen extinction coefficient within the canopy (Hikosaka et al. 2016)
	double extkn;
    /// The acclimation sensitivity coefficient for respiration (1/degC).
    /// Positive values reduce the respiration scalar when tgrowth is above 15
    /// degC and increase it when tgrowth is below 15 degC, with the strength of
    /// the effect determined by racclim. A value of zero means no acclimation
    /// effect.
	double racclim;

	// Temperature sensitivity coefficient (°C). Controls how rapidly damage
	// changes with temperature; smaller values produce a steeper damage
	// response curve.
	double z;

	// Reference temperature (°C) at which 50% leaf damage occurs after 1 minute
	// of exposure.
	double t50;

	/// Turgor loss point (MPa) - used for grasses only, where this often comes long before stomatal closure
	double tlp;
	// Sensitivity parameter indicating steepness of decline in Vcmax with water stress (Zhou 2014)
	double s_f;
	// Reference value at which Vcmax or Jmax is half it's non-water-stresed value (Zhou 2014)
	double psi_f;
	/// Leaf lignin fraction
	double ligcfrac_leaf;
	/// Root lignin fraction
	double ligcfrac_root;
	/// Wood lignin fraction
	double ligcfrac_wood;

	/// minimum leaf C:N mass ratio allowed when nitrogen demand is determined
	double cton_leaf_min;
	/// average leaf C:N mass ratio (between min and max)
	double cton_leaf_avr;
	/// maximum leaf C:N mass ratio	allowed when nitrogen demand is determined
	double cton_leaf_max;
	/// minimum fine root C:N mass ratio (connected to cton_root_avr and cton_root_max)
	double cton_root_min;
	/// average fine root C:N mass ratio (connected cton_leaf_avr)
	double cton_root_avr;
	/// maximum fine root C:N mass ratio (used when mass is negligible)
	double cton_root_max;
	/// minimum reproductive C:N mass ratio (connected to cton_root_avr and cton_root_max)
	double cton_repr_min;
	/// maximum reproductive C:N mass ratio (used when mass is negligible)
	double cton_repr_max;
	/// minimum sapwood C:N mass ratio (used when mass is negligible)
	double cton_sap_min;
	/// average sapwood C:N mass ratio (connected cton_leaf_avr)
	double cton_sap_avr;
	/// maximum sapwood C:N mass ratio (used when mass is negligible)
	double cton_sap_max;
	/// average heartwood C:N mass ratio (connected cton_leaf_avr)
	double cton_heart_avr;
	/// maximum heartwood C:N mass ratio (used when mass is negligible)
	double cton_heart_max;
	/// minimum heartwood C:N mass ratio (used when mass is negligible)
	double cton_heart_min;
	/// reference fine root C:N mass ratio
	double cton_root;
	/// reference sapwood C:N mass ratio
	double cton_sap;
	/// Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
	double nuptoroot;
	/// fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
	double fnstorage;

	/// Michaelis-Menten kinetic parameters
	/** Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
	double km_volume;

	/// Development and Growth parameters
    /// temperature below which plants fall dormant (degC)
	double temp_dormancy;
	/// Maximum NSC mobilisation rate from C storage (d-1)
	double max_mobil_cstorage;
    /// activation energy for growth temperature relationship (kJ mol-1) (Parent & Tardieu 2012)
	double Ha_growth;
	/// Parameter 'T0' in growth temperature relationship (K) (Parent & Tardieu 2012)
	double T0_growth;
	/// Parameter 'A' in growth temperature relationship (K-1) (Parent & Tardieu 2012)
	double A_growth;
	/// fraction of NPP allocated to reproduction
	double reprfrac;
	/// fraction of NPP allocated to crowns
	double crownfrac;
	// fraction of shoot allocated to culms (=culm / (culm + leaves))
	double culmfrac;
	/// Fraction of senesced carbon retranslocated back to storage.
	double crelocfrac;
	/// Baseline senescence rate for dave grasses (d-1).
	double senesc_base;
	/// Senescence coefficient for freezing temperatures in dave greasses
	double senesc_scal_temp;
	/// Senescence coefficient for dormancy in dave grasses. Senescence rate is
	/// increased by this factor when phenology is changing at the maximum rate.
	double senesc_scal_dorm;
	/// Senescence coefficient for water stress in dave grasses. Senescence rate
	/// is increased by this factor under max water stress.
	double senesc_scal_water;
	/// Senescence coefficient for reproductive phase in dave grasses.
	/// Senescence rate is increased by this factor during reproductive phase.
	double senesc_scal_repr;
	/// Exponent shaping development-driven senescence for dave grasses. Higher
	/// values delay senescence (less early, more late in the reproductive
	/// phase). Unitless.
	double qdev_senesc;
	/// Exponent shaping water stress-driven senescence for dave grasses. q = 1
	/// gives a linear response. q > 1 increases senescence at moderate stress;
	/// 0 < q < 1 reduces it. No effect at the extremes (wscal = 0 or 1).
	/// Unitless.
	double qwater_senesc;

	// Maximum NSC concentration in leaves.
	double max_leaf_nsc_conc;

	// Maximum NSC concentration in wood.
	double max_wood_nsc_conc;

	// Maximum depth of stems to which NSCs are stored (cm).
	double nsc_storage_depth;

	// Fraction of dead woody biomass that is bark.
	double bark_frac;

	/// annual leaf turnover as a proportion of leaf C biomass
	double turnover_leaf;
	/// annual fine root turnover as a proportion of fine root C biomass
	double turnover_root;
	/// annual turnover of reproductive litter to root litter as a proportion of reproductive litter C biomass (patch-basis)
	double turnover_repr_litter;
	/// annual turnover of "live" reproductive tissue to litter as a proportion of reproductive C biomass
	double turnover_repr;
	/// Scalar for initial tree biomass.
	double sapsize;
	/// annual sapwood turnover as a proportion of sapwood C biomass
	double turnover_sap;
	/// sapwood and heartwood density (kgC/m3)
	double wooddens;
	/// maximum tree crown area (m2)
	double crownarea_max;
	/// constant in allometry equations
	double k_allom1;
	/// constant in allometry equations
	double k_allom2;
	/// constant in allometry equations
	double k_allom3;
	/// constant in allometry equations
	double k_rp;
	/// tree leaf to sapwood area ratio
	double k_latosa;
	/// specific leaf area (m2/kgC)
	double sla;
	/// fraction of reproductive biomass consisting of seeds (the rest being supporting structures)
	double frac_seed;
	/// leaf to root mass ratio under non-water-stressed conditions
	double ltor_max;
	/// litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// fire resistance (0-1)
	double fireresist;
	/// minimum forest-floor PAR level for growth (grasses) or establishment (trees)
	/** J/m2/day, individual and cohort modes */
	double parff_min;
	/** parameter capturing non-linearity in recruitment rate relative to
	 *  understorey growing conditions for trees (Fulton 1991) (individual and
	 *  cohort modes)
	 */
	double alphar;

	/// Form factor to use in tree height calculations that takes into account
	/// the relationship between the difference in tree mass between the
	/// cylindrical stem hypothesis and measured mass of tree log with its
	/// branches.
	double formfactor;

	/// Whether vcmax is calculated from leaf nutrient content (alt: prescribed)
	bool ifcalcvcmax;

	/*
	Note: these 3 parameters are only used by dave photosynthesis, and only when
	ifcalcvcmax is false (this is a parameter in the .ins file). When false,
	these values must be provided in the .ins file for each PFT (only vpmax_25
	required for C4 grasses, only vcmax_25 and jmax_25 required for C3). When
	ifcalcvcmax is true, values for these will be estimated daily based on leaf
	N (and P) content, and stored in the indiv.ps_subdaily object.
	*/
	// (Constant) Maximum carboxylation rate at 25 degC (μmol m-2 s-1)
	double vcmax_25;
	// (Constant) Maximum electron transport rate at 25 degC (μmol m-2 s-1).
	double jmax_25;
	// (Constant) Maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1)
	double vpmax_25;

	// Slope of relationship between vcmax25 and leaf N content in log scale.
	double vcmax25n_slope;

	// Slope of relationship between vcmax25 and leaf P content in log scale.
	double vcmax25p_slope;

	// Intercept of relationship between vcmax25 and leaf nutrient content.
	double vcmax25_intercept;

	// Slope of relationship between jmax25 and leaf N content in log scale.
	double jmax25n_slope;

	// Slope of relationship between jmax25 and leaf P content in log scale.
	double jmax25p_slope;

	// Intercept of relationship between jmax25 and leaf nutrient content.
	double jmax25_intercept;

	// Slope of relationship between vpmax25 and leaf N content in log scale.
	double vpmax25n_slope;

	// Intercept of relationship between vpmax25 and leaf nutrient content.
	double vpmax25_intercept;

	// Activation energy for vcmax (kJ mol-1)
	double ha_v;
	// Activation energy for jmax (kJ mol-1)
	double ha_j;
	// Entropy terms (deltaS) for vcmax (J mol-1 K-1)
	double ds_v;
	// Entropy terms (deltaS) for jmax (J mol-1 K-1)
	double ds_j;

	// Bole height as a fraction of tree height (0-1).
	double boleht_frac;

	// Fraction of PAR absorbed by the leaf area (0-1).
	double absorptance_par;

	// Fraction of NIR absorbed by the leaf area (0-1).
	double absorptance_nir;

	// Quantum yield of electron transport (mol (electrons) mol-1 (photon)).
	double alpha;

	// Leaf width (m).
	double leaf_width;

	// Inflection point of sigmoidal carbon starvation as a function of storage
	// fraction full. Higher values mean carbon starvation starts at higher
	// storage levels.
	double cstarv_inflect;

	// Maximum carbon starvation rate (fraction of cohort killed per year) when
	// storage is low.
	double cstarv_max;

	// Theshold for carbon starvation mortality killing the entire cohort,
	// expressed in density terms (number of trees per m2).
	double cstarv_killdens;

	/// maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
	double est_max;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect enabled (individual and cohort modes)
	 */
	double kest_repr;
	/// constant affecting amount of background establishment
	/** \see ifbgestab */
	double kest_bg;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect disabled (individual and cohort modes)
	 */
	double kest_pres;
	/// expected longevity under non-stressed conditions (individual and cohort modes)
	double longevity;
	/// Above ground woody biomass fraction
	double ag_woody_frac;
	/// threshold growth efficiency for imposition of growth suppression mortality
	/** kgC/m2 leaf/year, individual and cohort modes */
	double greff_min;

	// Bioclimatic limits (all temperatures deg C)

	/// minimum 20-year coldest month mean temperature for survival
	double tcmin_surv;
	/// maximum 20-year coldest month mean temperature for establishment
	double tcmax_est;
	/// minimum degree day sum on 5 deg C base for establishment
	double gdd5min_est;
	/// minimum 20-year coldest month mean temperature for establishment
	double tcmin_est;
	/// minimum warmest month mean temperature for establishment
	double twmin_est;
	/// continentality parameter for boreal summergreen trees
	double twminusc;
	/// constant in equation for budburst chilling time requirement (Sykes et al 1996)
	double k_chilla;
	/// coefficient in equation for budburst chilling time requirement
	double k_chillb;
	/// exponent in equation for budburst chilling time requirement
	double k_chillk;
	/// array containing values for GDD0(c) given c=number of chill days
	/** Sykes et al 1996, Eqn 1
	 *  gdd0 has one element for each possible value for number of chill days
	 */
	double gdd0[Date::MAX_YEAR_LENGTH + 1];
	/// interception coefficient (unitless)
	double intc;

	/// nitrogen fertilisation amount (kgN/m2/year)
	double N_appfert;
	/// 0 - 1 how much of the fertiliser is applied at the two fertdates, default 1.
	double fertrate[2];
	/// fertilisation dates relative to sowing date
	int fertdates[2];
	/// development stage at the two fertilisation events after sowing
	double fert_stages[2];

	double T_vn_min;
	double T_vn_opt;
	double T_vn_max;

	double T_veg_min;
	double T_veg_opt;
	double T_veg_max;

	double T_rep_min;
	double T_rep_opt;
	double T_rep_max;

	double photo[3];

	double dev_rate_veg;
	double dev_rate_rep;

	double a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
	double cton_stem_avr;
	double cton_stem_max;

	/// Drought tolerance level (0 = very -> 1 = not at all) (unitless)
	/** Used to implement drought-limited establishment */
	double drought_tolerance;

	// bvoc

	/// aerodynamic conductance (m s-1)
	double ga;
	/// isoprene emission capacity (ug C g-1 h-1)
	double eps_iso;
	/// whether (1) or not (1) isoprene emissions show a seasonality
	bool seas_iso;
	/// monoterpene emission capacity (ug C g-1 h-1) per monoterpene species
	double eps_mon[NMTCOMPOUNDS];
	/// fraction of monoterpene production that goes into storage pool (-) per monoterpene species
	double storfrac_mon[NMTCOMPOUNDS];

	/// Bioclimatic limits parameters from Wolf et al. 2008

	/// snow max [mm]
	double max_snow;
	/// snow min [mm]
	double min_snow;
	/// GDD0 min
	double gdd0_min;
	/// GDD0 max
	double gdd0_max;

	/// New parameters from parameters from Wania et al. (2009a, 2009b, 2010)

	/// Days per month for which inundation is tolerated
	int inund_duration;
	/// Inundation stress is felt when the water table (mm) is above wtp_max
	double wtp_max;
	/// Whether this PFT has aerenchyma through which O2 and CH4 can be transported (Wania et al. 2010 - Sec 2.6)
	bool has_aerenchyma;

	/// Sapling/regeneration characteristics (used only in population mode)
	/** For trees, on sapling individual basis (kgC); for grasses, on stand area basis,
	 *  kgC/m2 */
	struct {
		/// leaf C biomass
		double cmass_leaf;
		/// fine root C biomass
		double cmass_root;
		/// sapwood C biomass
		double cmass_sap;
		/// heartwood C biomass
		double cmass_heart;
	} regen;

	/// specifies type of landcover pft is allowed to grow in (0 = URBAN, 1 = CROP, 2 = PASTURE, 3 = FOREST, 4 = NATURAL, 5 = PEATLAND)
	landcovertype landcover;
	/// pft selection
	xtring selection;
	/// fraction of residue outtake at harvest
	double res_outtake;
	/// plant density after clearcut (seedlings/ha)
	double plantdensity;
	/// harvest efficiency
	double harv_eff;
	/// harvest efficiency of intercrop grass
	double harv_eff_ic;
	/// fraction of harvested products that goes into patchpft.harvested_products_slow
	double harvest_slow_frac;
	/// yearly turnover fraction of patchpft.harvested_products_slow (goes to gridcell.acflux_harvest_slow)
	double turnover_harv_prod;
	/// Fraction of wood cmass that belongs to stems
	double stem_frac;
	/// Fraction of wood cmass that belongs to twigs
	double twig_frac;
	/// whether pft may grow as cover crop
	bool isintercropgrass;
	/// whether autumn temperature dependent sowing date is calculated
	bool ifsdautumn;
	/// upper temperature limit for autumn sowing
	double tempautumn;
	/// lower temperature limit for spring sowing
	double tempspring;
	/// default length of growing period
	int lgp_def;
	/// upper minimum temperature limit for crop sowing
	double maxtemp_sowing;
	/// default sowing date in the northern hemisphere (julian day)
	int sdatenh;
	/// default sowing date in the southern hemisphere
	int sdatesh;
	/// whether sowing date adjusting equation is used
	bool sd_adjust;
	/// parameter 1 in sowing date adjusting equation
	double sd_adjust_par1;
	/// parameter 2 in sowing date adjusting equation
	double sd_adjust_par2;
	/// parameter 3 in sowing date adjusting equation
	double sd_adjust_par3;
	/// latest date for harvesting in the northern hemisphere
	int hlimitdatenh;
	/// latest date for harvesting in the southern hemisphere
	int hlimitdatesh;
	/// default base temperature (°C) for heat unit (hu) calculation
	double tb;
	/// temperature under which vernalisation is possible (°C)
	double trg;
	/// default number of vernalising days required
	int pvd;
	/// sensitivity to the photoperiod effect [0-1]
	double psens;
	/// basal photoperiod (h) (pb<ps for longer days plants)
	double pb;
	/// lag in days after sowing before vernalization starts
	int vern_lag;
	/// saturating photoperiod (h) (ps<pb for shorter days plants)
	double ps;
	/// default potential heat units required for crop maturity (degree-days)
	double phu;
	/// whether quadratic equation used for calculating potential heat units (Bondeau method)
	bool phu_calc_quad;
	/// whether linear equation used for calculating potential heat units (Bondeau method)
	bool phu_calc_lin;
	/// minimum potential heat units required for crop maturity (Bondeau method) (degree-days)
	double phu_min;
	/// maximum potential heat units required for crop maturity (Bondeau method) (degree-days)
	double phu_max;
	/// reduction factor of potential heat units in spring crops (Bondeau method) (degree-days)
	double phu_red_spring_sow;
	/// number of days of phu decrease in the linear phu equation (Bondeau method)
	double ndays_ramp_phu;
	/// intercept for the linear phu equation (Bondeau method)
	double phu_interc;
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
	/// autumn/spring sowing of pft:s with tempautumn = 1
	int forceautumnsowing;	//0 = NOFORCING,  1 = AUTUMNSOWING, 2 = SPRINGSOWING

	double avg_cton(const double& min, const double& max) {
		return 2.0 / (1. / min + 1. / max);
	}

    bool is_grass() const {
		return lifeform == GRASS_PERENNIAL || lifeform == GRASS_ANNUAL;
	}

	bool store_uptake() const {
		return dave_allocation() && ifstoreuptake || is_grass();
	}

	/// Check whether DAVE allocation is enabled for this PFT.
	bool dave_allocation() const;

	/// Return true if this PFT is using daily establishment, or false otherwise
	/// (ie if it's using the annual (trunk) establishment).
	bool daily_establishment() const;

	/// transfer constant of leaf material between compartments (used for daily grass growth)
	double transferconst;

	// MEMBER FUNCTIONS

public:

	/// Constructor (initialises array gdd0)
	Pft() {

		std::fill_n(gdd0, Date::MAX_YEAR_LENGTH + 1, -1.0); // value<0 signifies "unknown"; see function phenology()

		root_beta = 0.0;

		drought_tolerance = 0.0; // Default, means that the PFT will never be limited by drought.
		res_outtake = 0.0;
		plantdensity = 0.0;
		harv_eff = 0.0;
		harv_eff_ic = 0.0;
		turnover_harv_prod = 1.0;	// default 1 year turnover time

		vcmax_25 = 0;
		vpmax_25 = 0;
		jmax_25  = 0;
		vcmax25n_slope = 0;
		vcmax25p_slope = 0;
		vcmax25_intercept = 0;
		jmax25n_slope = 0;
		jmax25p_slope = 0;
		jmax25_intercept = 0;
		vpmax25n_slope = 0;
		vpmax25_intercept = 0;
		vpmax_scal = 0;
		ha_v = 0;
		ha_j = 0;
		ds_v = 0;
		ds_j = 0;
		formfactor = 1.0;
		sapsize = 0.1;
		boleht_frac = 0;
		absorptance_par = 0.0;
		absorptance_nir = 0.0;
		alpha = 0.0;
		leaf_width = 0.0;
		cstarv_inflect = 0.0;
		cstarv_max = 0.0;
		cstarv_killdens = 0.0;
		max_mobil_cstorage = 0;

		isintercropgrass = false;
		ifsdautumn = false;
		maxtemp_sowing = 60;
		sdatenh = -1;
		sdatesh = -1;
		lgp_def = 190;
		hlimitdatenh = -1;
		hlimitdatesh = -1;
		tb = -999.9;
		trg = -999.9;
		pvd = -1;
		psens = -1.0;
		pb = -1.0;
		vern_lag=0;
		ps = -1.0;
		pstemp_high = 0.0;
		pstemp_low = 0.0;	
		pstemp_max = 0.0;
		pstemp_min = 0.0;
		phu = -1.0;
		phu_red_spring_sow = 1.0;
		fphusen = -1.0;
		shapesenescencenorm = 0;
		flaimaxharvest = -1.0;
		laimax = 0.0;
		aboveground_ho = true;
		frootstart = 0.0;
		frootend = 0.0;
		forceautumnsowing = 0;
		fertrate[0] = 0.0;
		fertrate[1] = 1.0;
		fertdates[0] = 0;
		fertdates[1] = 30;

		fert_stages[0] = 0.5;
		fert_stages[1] = 0.9;

		N_appfert = 0.0;

		T_vn_min=0.0;
		T_vn_opt=0.0;
		T_vn_max=0.0;
		T_veg_min=0.0;
		T_veg_opt=0.0;
		T_veg_max=0.0;
		T_rep_min=0.0;
		T_rep_opt=0.0;
		T_rep_max=0.0;

		a1=0.0;
		b1=0.0;
		c1=0.0;
		d1=0.0;
		a2=0.0;
		b2=0.0;
		c2=0.0;
		d2=0.0;
		a3=0.0;
		b3=0.0;
		c3=0.0;
		d3=0.0;

		g0=0.0;
		g1=0.0;
		gcutic=0.0;
		qs=0.0;

		stom_clos=0.0;
		p50=0.0;
		p88=0.0;
		kmax=0.0;
		g1x=0.0;
		tlp=0.0;
		ifcalcvcmax = false;

		vcmax_25=0.0;
		turnover_repr = 0.0;
		crelocfrac = 0.0;
		senesc_base = 0;
		senesc_scal_temp = 0;
		senesc_scal_dorm = 0;
		senesc_scal_water = 0;
		senesc_scal_repr = 0;
		qdev_senesc = 0;
		qwater_senesc = 0;
		Ha_growth = 0.0;
		T0_growth = 0.0;
		temp_dormancy = 0.0;
		extkn = 0;
		racclim = 0;
		for (int i=0; i<3; i++)
			photo[i]=0.0;

		for (int i = 0; i<NSOILLAYERS_SIMPLE; i++)
			root_dist[i] = 0.0;
	}

	/// Calculates SLA given leaf longevity
	void initsla() {

		// SLA has to be supplied in the insfile for crops with N limitation
		if (!(phenology == CROPGREEN && ifnlim)) {

			// Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
			// m2/kgC)

			if (leafphysiognomy == BROADLEAF) {
				sla = 0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * leaflong()));
			}
			else if (leafphysiognomy == NEEDLELEAF) {
				sla = 0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * leaflong()));
			}
		}
	}

	////
	void init_hydraulic_vulnerability_curves() {
		s2 = log(log(1.0-50.0/100.0)/log(1.0-88.0/100.0))/(log(abs(p50))-log(abs(p88)));
		s1 = abs(p50)/(pow(-log(1-50.0/100),(1/s2)));
	}

	// Set spline for proportion of conductivity.
	void setup_transpiration(double resolution);

	/// Calculates minimum leaf C:N ratio given leaf longevity
	void init_cton_min() {

		// cton_leaf_min has to be supplied in the insfile for crops with N limitation
		if (!(phenology == CROPGREEN && ifnlim)) {
			// Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
			// kgN/kgC)

			if (leafphysiognomy == BROADLEAF)
				cton_leaf_min = 500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * leaflong()));
			else if (leafphysiognomy == NEEDLELEAF)
				cton_leaf_min = 500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * leaflong()));
		}
	}

	// Calculate leaf longevity in years.
	double leaflong() {
		return 1.0 / turnover_leaf;
	}

	void init_cton_limits() {

		// Fraction between min and max C:N ratio White et al. 2000
		double frac_mintomax = (phenology == CROPGREEN && ifnlim) ? 5.0 : 2.78;	// Use value also without nlim ?

		// Fraction between leaf and root C:N ratio
		double frac_leaftoroot = 1.16; // Friend et al. 1997

		// Fraction between leaf and sap wood C:N ratio
		double frac_leaftosap = 6.9;   // Friend et al. 1997

		// Max leaf C:N ratio
		cton_leaf_max = cton_leaf_min * frac_mintomax;

		if (dave_allocation() && is_grass()) {
			// TODO: revisit C:N avg!
			// Bounds for fine root C:N ratio
			cton_root_min = cton_leaf_min * frac_leaftoroot;
			cton_root_max = cton_leaf_max * frac_leaftoroot;

			// Bounds for sap C:N ratio
			cton_sap_min = cton_leaf_min * frac_leaftosap;
			cton_sap_max = cton_leaf_max * frac_leaftosap;
		} else {
			// Tighter C:N ratio range for roots and sapwood: picked out thin air
			double frac_maxtomin = .9;

			// Maximum fine root C:N ratio
			cton_root_max = cton_leaf_max * frac_leaftoroot;
			cton_root_min = cton_root_max * frac_maxtomin;

			// Maximum sap C:N ratio
			cton_sap_max  = cton_leaf_max * frac_leaftosap;
			cton_sap_min = cton_sap_max * frac_maxtomin;
		}

		// Bounds for repr C:N ratio
		cton_repr_min = std::min(cton_root_min, cton_leaf_min);
		cton_repr_max = std::max(cton_root_max, cton_leaf_max);

		// Average leaf C:N ratio
		cton_leaf_avr = avg_cton(cton_leaf_min, cton_leaf_max);

		// Tighter C:N ratio range for roots and sapwood: picked out thin air
		double frac_fromavr = 1.05;

		// Average fine root C:N ratio
		cton_root_avr = cton_leaf_avr * frac_leaftoroot;

		// Minimum fine root C:N ratio
		cton_root_min = cton_root_avr / frac_fromavr;

		// Maximum fine root C:N ratio
		cton_root_max = cton_root_avr * frac_fromavr;

		// Average sap C:N ratio
		cton_sap_avr = cton_leaf_avr * frac_leaftosap;

		// Minimum sap C:N ratio
		cton_sap_min = cton_sap_avr / frac_fromavr;

		// Maximum sap C:N ratio
		cton_sap_max = cton_sap_avr * frac_fromavr;

		// Average heart C:N ratio
		cton_heart_avr = cton_sap_avr * 2.0;

		// Minimum heart C:N ratio
		cton_heart_min = cton_heart_avr / frac_fromavr;

		// Maximum heart C:N ratio
		cton_heart_max = cton_heart_avr * frac_fromavr;

		if (!dave_allocation()) {
			if (is_grass() || lifeform == MOSS) {
				respcoeff /= 2.0 * cton_root / (cton_root_avr + cton_root_min);
			}
			else {
				respcoeff /= cton_root / (cton_root_avr + cton_root_min) +
			             cton_sap  / (cton_sap_avr  + cton_sap_min);
			}
		}

		cton_stem_max = 1.0/(2.0*0.0034); //Maize params
		cton_stem_avr = 1.0/(2.0*0.0068);
	/*
		printf("%s: cton_leaf_min = %.2f\n", (const char *)name, cton_leaf_min);
		printf("%s: cton_leaf_avr = %.2f\n", (const char *)name, cton_leaf_avr);
		printf("%s: cton_leaf_max = %.2f\n", (const char *)name, cton_leaf_max);

		printf("%s: cton_root_min = %.2f\n", (const char *)name, cton_root_min);
		printf("%s: cton_root_avr = %.2f\n", (const char *)name, cton_root_avr);
		printf("%s: cton_root_max = %.2f\n", (const char *)name, cton_root_max);

		printf("%s: cton_sap_min = %.2f\n", (const char *)name, cton_sap_min);
		printf("%s: cton_sap_avr = %.2f\n", (const char *)name, cton_sap_avr);
		printf("%s: cton_sap_max = %.2f\n", (const char *)name, cton_sap_max);
	*/
	}

	/// Calculates leaf lignin fraction from leaf SLA
	void init_ligcfrac_leaf() {

		// The lignin content is assumed constant for all but the leaf tissues. For the latter, an empirical 
		// dependency between lignin content and SLA is used (White et al., 2000).
		const double ligcfrac_leaf_max = 0.344;
		const double kleaf2sla = -0.0052;

		ligcfrac_leaf = ligcfrac_leaf_max + kleaf2sla * sla;

		if (ligcfrac_leaf < 0.0) {
			fail("ERROR: Negative ligcfrac_leaf %g SLA too high %g for %s\n", ligcfrac_leaf, sla, (char*)name);
		}
	}

	/// Initialises sapling/regen characteristics in population mode following LPJF formulation
	void initregen() {

		// see function allometry in growth module.

		// Note: primary PFT parameters, including SLA, must be set before this
		//       function is called

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
		else if (is_grass() || lifeform==MOSS) {

			// Grass regeneration characteristics

			regen.cmass_leaf = REGENLAI_GRASS / sla;
		}

		regen.cmass_root = 1.0 / ltor_max * regen.cmass_leaf;
	}

	/// Initialises transfer constant of leaf material between compartments for daily growth of grass
	void init_transferconst() {

		if (lifeform == TREE) {

			// Not used for trees yet

			transferconst = 0.0;
		}
		else if (is_grass()) {

			// Used for daily grass growth. Results in correct leaf longevity under 0.75 * perfect conditions

			//transferconst = 0.01 / leaflong; // DW_COMMENT to be calibrated!!! 0.0091 / leaflong;
			if (pathway == C3) {
				transferconst = 0.02 / leaflong(); // DW_COMMENT to be calibrated!!! 0.0091 / leaflong;
				//transferconst = 0.005 / leaflong; // Sturt Plains: 0.005, Yanco: 0.02, Otway: 0.025
			} 
			else {
				transferconst = 0.022 / leaflong(); 
				//transferconst = 0.01 / leaflong; // Otway: 0.02, Yanco: 0.022
			}
		}
	}

	// Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
	void init_rootdist(double soil_lyr_depth[NSOILLAYERS_MAX], double soil_depth, int nsoillayers) {

		if (hydrology == GERTEN) {

			// Old two layer scheme
			input_rootdist[0] = root_dist[0];
			input_rootdist[1] = root_dist[1];
		}
		else {

			// Multi layer soil scheme uses Jackson root distribution
			// root_beta is calibrated for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.

			// Standard depth
			const double depth_standard = SOILDEPTH_UPPER + SOILDEPTH_LOWER;

			// Standard root depth
			double depth_root = depth_standard;

			// Fraction of fine roots to be homogenously distributed within root depth
			// is zero for MILLER hydrology
			double fhom_root = 0.0;

			if (hydrology == RICHARDS) {

				// Fraction of deep soil, beyond standard depth, that roots have access to
				double fdeep_root = rootfrac_of_deep_soil;

				// Fraction of fine roots to be homogenously distributed within root depth
				// is set in the ins file
				fhom_root = homrootdistfrac;

				// Grass has half the amount of homogenously distributed roots and depth of deep roots
				if (isgrass()) {
					fdeep_root *= 0.5;
					fhom_root *= 0.5;
				}

				// Moss don't have homogenously distributed roots nor deep roots
				if (ismoss()) {
					fdeep_root = 0.0;
					fhom_root = 0.0;
				}

				// If soil depth is deeper than standard depth, then update root depth
				if (soil_depth > depth_standard) {
					depth_root = depth_standard + fdeep_root * (soil_depth - depth_standard);
				}
			}

			double rootdist_new[NSOILLAYERS_MAX] = { 0.0 };

			calc_rootdist(rootdist_new, soil_lyr_depth, nsoillayers, depth_root, root_beta, fhom_root);

			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				input_rootdist[lyr] = rootdist_new[lyr];
			}
		}
	}

	// Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
	void calc_rootdist(double rootdist_new[NSOILLAYERS_MAX], double soil_lyr_depth[NSOILLAYERS_MAX], int nsoillayers,
		double depth_root, double root_beta, double fhom_root) const{

		// Multi layer soil scheme uses Jackson root distribution
		// root_beta is calibrated for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.

		// Root depth in cm
		double depth_root_cm = depth_root * CM_PER_MM;

		// Total root dist
		double tot = 0.0;

		// Total soil layer depth
		double tot_depth = soil_lyr_depth[0] * CM_PER_MM;

		// Total previous soil layer depth
		double tot_depth_prev = 0.0;

		// init first layer
		rootdist_new[0] = 1.0 - pow(root_beta, tot_depth); // - (1.0 - pow(root_beta, tot_depth_prev));
		tot += rootdist_new[0];

		for (int lyr = 1; lyr < nsoillayers; lyr++) {

			// Update depths
			tot_depth_prev = tot_depth;
			tot_depth += soil_lyr_depth[lyr] * CM_PER_MM;

			// Check so depth doesn't exceed rooting depth
			if (tot_depth >= depth_root_cm) {
				tot_depth = depth_root_cm;
			}

			double root_frac = 1.0 - pow(root_beta, tot_depth) - (1.0 - pow(root_beta, tot_depth_prev));

			// Don't allow too small root fractions
			rootdist_new[lyr] = root_frac > 0.001 ? root_frac : 0.0;
			tot += rootdist_new[lyr];
		}

		// Adjust so rootdist always add up to 1
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			rootdist_new[lyr] /= tot;
		}

		// If a fraction of root distribution should be homogenously distributed
		if (!negligible(fhom_root)) {

			// homogenous root fraction per depth [mm-1]
			double frac_per_depth = fhom_root / depth_root;

			// Reset total depth [mm]
			tot = 0.0;

			// Accomolated depth [mm]
			double acc_depth = 0.0;

			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				double depth = soil_lyr_depth[lyr];
				if (acc_depth >= depth_root) {
					depth = 0.0;
				}
				else if (acc_depth + depth > depth_root) {
					depth = depth_root - acc_depth;
				}

				// New root dist
				rootdist_new[lyr] = rootdist_new[lyr] * (1.0 - fhom_root) + frac_per_depth * depth;

				// Accomolated total root dist
				tot += rootdist_new[lyr];

				// Update accomolated depth
				acc_depth += depth;
			}

			// Adjust so rootdist always add up to 1
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				rootdist_new[lyr] /= tot;
			}
		}
	}

	bool ismoss() const {

		if (lifeform == MOSS)
			return true;
		else
			return false;
	}

	bool isgrass() const {

		if (is_grass())
			return true;
		else
			return false;
	}

	bool istree() const {

		if (lifeform == TREE)
			return true;
		else
			return false;
	}

	bool iswetlandspecies() const {

		if (ismoss() || has_aerenchyma)
			return true;
		else
			return false;
	}

	bool is_shrub() const {

		// Current limit for shrubs (10 m); may be changed.
		return crownarea_max < 10.0 && lifeform == TREE;
	}

	bool is_shade_intolerant_tree() const {

		// Current limit for shade intolerant trees (10 in europe.ins); may be changed.
		return alphar >= 10.0 && lifeform == TREE;
	}

	bool is_boreal() const {

		// Current value of pstemp_low for boreal trees and cold (C3) grass in global.ins and europe.ins is 10; may be changed.
		return pstemp_low <= 12;
	}

	bool is_temperate() const {

		// Current value of pstemp_low for temperate trees in global.ins and europe.ins is 15; may be changed.
		return pstemp_low > 12 && pstemp_low <= 20;
	}

	bool is_tropical() const {

		// Current value of pstemp_low for tropical trees and warm (C4) grass in global.ins and europe.ins is 20; may be changed.
		return pstemp_low > 20;
	}

	bool get_psi_from_transpiration(double transpiration, double &psi) const;
	bool get_transpiration_lower_bound(double psi, double &transpiration) const;
private:
	Interpolator psi_from_transpiration;
	Interpolator transpiration_from_psi;
	double proportion_of_conductivity(double psi) const;
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
class Pftlist : public ListArray_id<Pft> {

public:
	int getpftid(xtring pftname) {

		int id = -1;

		for(unsigned int i=0; i< this->nobj; i++) {

			Pft& pft = (*this)[i];
			if(pft.name == pftname) {
				id = pft.id;
				break;
			}
		}

		return id;
	}
};

/// The one and only linked list of Pft objects
extern Pftlist pftlist;

/// Container for crop-specific data at the individual level
struct cropindiv_struct : public Serializable {

	// Plant carbon biomass variables are all on patch area basis (kgC/m2)

	/// year's harvestable organ C biomass (= ycmass_plant)
	double cmass_ho;
	/// above-ground pool C biomass (when calculating daily cmass_leaf from lai_crop) (= ycmass_agpool)
	double cmass_agpool;
	double cmass_stem;

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
	double dcmass_stem;

	/// today's increase of leaf N biomass
	double dnmass_leaf;
	/// today's increase of root N biomass
	double dnmass_root;
	/// today's increase of harvestable organ N biomass
	double dnmass_ho;
	/// today's increase of above-ground pool N biomass
	double dnmass_agpool;

	/// CARBON
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
	/// daily updated dead leaf C biomass, reset at harvest day
	double grs_cmass_dead_leaf;
	/// daily updated stem pool C biomass, reset at harvest day
	double grs_cmass_stem;

	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_leaf_luc;
	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_root_luc;
	/// carbon content of harvestable organs saved on first day of land use change year
	double grs_cmass_ho_luc;
	/// carbon content of above-ground pool saved on first day of land use change year
	double grs_cmass_agpool_luc;
	/// carbon content of dead leaves saved on first day of land use change year
	double grs_cmass_dead_leaf_luc;
	/// carbon content of stem saved on first day of land use change year
	double grs_cmass_stem_luc;

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
	/// daily updated dead leaf C biomass, reset at day 0
	double ycmass_dead_leaf;
	/// daily updated stem C biomass, reset at day 0
	double ycmass_stem;

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
	/// year's stem C biomass at time of harvest (cumulative if several harvest events)
	double harv_cmass_stem;

	/// Nitrogen
	/// nitrogen content of harvestable organs
	double nmass_ho;
	/// nitrogen content of above-ground pool
	double nmass_agpool;
	/// nitrogen content of dead leaves
	double nmass_dead_leaf;

	/// nitrogen content of harvestable organs saved on first day of land use change year
	double nmass_ho_luc;
	/// nitrogen content of above-ground pool saved on first day of land use change year
	double nmass_agpool_luc;
	/// nitrogen content of dead leaves saved on first day of land use change year
	double nmass_dead_leaf_luc;

	/// daily updated leaf N biomass, reset at day 0
	double ynmass_leaf;
	/// daily updated root N biomass, reset at day 0
	double ynmass_root;
	/// daily updated harvestable organ N biomass, reset at day 0
	double ynmass_ho;
	/// daily updated above-ground pool N biomass, reset at day 0
	double ynmass_agpool;
	/// daily updated dead leaf N biomass, reset at day 0
	double ynmass_dead_leaf;

	/// year's leaf N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_leaf;
	/// year's root N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_root;
	/// year's harvestable organ N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_ho;
	/// year's above-ground pool N biomass at time of harvest (cumulative if several harvest events)
	double harv_nmass_agpool;

	/// dry weight crop yield harvested this year (cumulative if several harvest events), based on harv_cmass_xx
	double harv_yield;

	/// harvestable organ C biomass at the last two harvest events this year
	double cmass_ho_harvest[2];
	/// harvestable organ N biomass at the last two harvest events this year
	double nmass_ho_harvest[2];
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

	cropindiv_struct() {
		cmass_ho=0.0;
		cmass_agpool=0.0;
		cmass_stem = 0.0;
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
		dcmass_stem = 0.0;
		grs_cmass_leaf=0.0;
		grs_cmass_root=0.0;
		grs_cmass_plant=0.0;
		grs_cmass_ho=0.0;
		grs_cmass_agpool=0.0;
		grs_cmass_stem = 0.0;
		grs_cmass_dead_leaf = 0.0;
		grs_cmass_leaf_luc=0.0;
		grs_cmass_root_luc=0.0;
		grs_cmass_ho_luc=0.0;
		grs_cmass_agpool_luc=0.0;
		grs_cmass_dead_leaf_luc = 0.0;
		grs_cmass_stem_luc = 0.0;
		nmass_ho=0.0;
		nmass_agpool=0.0;
		nmass_dead_leaf = 0.0;
		ycmass_leaf=0.0;
		ycmass_root=0.0;
		ycmass_plant=0.0;
		ycmass_ho=0.0;
		ycmass_agpool=0.0;
		ycmass_stem = 0.0;
		ycmass_dead_leaf = 0.0;
		harv_cmass_leaf=0.0;
		harv_cmass_root=0.0;
		harv_cmass_root=0.0;
		harv_cmass_ho=0.0;
		harv_yield=0.0;
		harv_cmass_agpool=0.0;
		harv_cmass_stem = 0.0;
		cmass_ho_harvest[0]=0.0;
		cmass_ho_harvest[1]=0.0;

		//Nitrogen
		dnmass_leaf=0.0;
		dnmass_root=0.0;
		dnmass_ho=0.0;
		dnmass_agpool=0.0;
		ynmass_leaf=0.0;
		ynmass_root=0.0;
		ynmass_ho=0.0;
		ynmass_agpool=0.0;
		ynmass_dead_leaf = 0.0;
		harv_nmass_leaf=0.0;
		harv_nmass_root=0.0;
		harv_nmass_root=0.0;
		harv_nmass_ho=0.0;
		harv_nmass_agpool=0.0;
		nmass_dead_leaf_luc = 0.0;
		nmass_ho_harvest[0]=0.0;
		nmass_ho_harvest[1]=0.0;

		isprimarycrop=false;
		isprimarycovegetation=false;
//		issecondarycrop=false;
		isintercropgrass=false;
	}

	void serialize(ArchiveStream& arch);
};


/// A vegetation individual.
/** In population mode this is the average individual of a PFT population;
 *  in cohort mode: the average individual of a cohort;
 *  in individual mode: an individual plant. Each grass PFT is represented as a single
 *  individual in all modes. Individual objects are collected within list arrays of
 *  class Vegetation (defined below), of which there is one for each patch, and include
 *  a reference to their 'parent' Vegetation object. Use the createobj member function
 *  of class Vegetation to add new individuals.
 */
class Individual : public Serializable {

public:
	/// reference to Pft object containing static parameters for this individual
	Pft& pft;
	/// reference to Vegetation object to which this Individual belongs
	Vegetation& vegetation;
	/// id code (0-based, sequential)
	int id;
	/// leaf C biomass on modelled area basis (kgC/m2)
	double cmass_leaf;
	/// fine root C biomass on modelled area basis (kgC/m2)
	double cmass_root;
	/// crown C biomass (kgC/m2)
	double cmass_crown;
	/// reproductive tissue C biomass (kgC/m2)
	double cmass_repr;
	/// sapwood C biomass on modelled area basis (kgC/m2)
	double cmass_sap;
	/// Total annual C biomass allocated to leaf on modelled area basis (kgC/m2)
	double acalloc_leaf;
	/// Total annual C biomass allocated to fine root on modelled area basis (kgC/m2)
	double acalloc_root;
	/// Total annual C biomass allocated to sapwood on modelled area basis (kgC/m2)
	double acalloc_sap;
	/// Total annaul C biomass allocated to reproductive tissue on modelled area basis (kgC/m2)
	double acalloc_repr;
	/// Total annaul C biomass allocated to storage tissue on modelled area basis (kgC/m2)
	double acalloc_storage;
	/// Total annaul C biomass allocated to crown on modelled area basis (kgC/m2)
	double acalloc_crown;
	/// heartwood C biomass on modelled area basis (kgC/m2)
	double cmass_heart;
	/// C "debt" (retrospective storage) (kgC/m2)
	double cmass_debt;
	/// Total C mass at land use change (kgC/m2)
	double cmass_tot_luc;
	/// leaf C mass after tunrnover
	double cmass_leaf_post_turnover;
	/// root C mass after turnover
	double cmass_root_post_turnover;
	/// total C in sapwood yesterday (kgC/m2)
	double cmass_sap_yesterday;
	/// amount of sapwood lost to turnover today (kgC/m2)
	double cmass_sap_turnover;
	/// drought affected dynamic hydraulic conductance of sapwood (mmol/m2/s/MPa)
	double kmax_drought_affected;
	/// daily loss of maximum hydraulic conductivity due to xylem embolism
	double max_k_loss;
	/// Daily fraction of xylem lost to hydraulic stress.
	double xylem_loss_frac;
	/// Daily carbon lost to xylem embolism (kgC/m2).
	double closs_embolism;
	/// Daily nitrogen lost to xylem embolism (kgN/m2).
	double nloss_embolism;
	/// Daily carbon lost to hydraulic failure (kgC/m2).
	double closs_hydraulic;
	/// Daily nitrogen lost to hydraulic failure (kgN/m2).
	double nloss_hydraulic;
	/// Daily carbon lost to blaze fire (kgC/m2).
	double closs_blaze;
	/// Daily nitrogen lost to blaze fire (kgN/m2).
	double nloss_blaze;
	/// Daily carbon lost to globfirm fire (kgC/m2).
	double closs_globfirm;
	/// Daily nitrogen lost to globfirm fire (kgN/m2).
	double nloss_globfirm;
	/// Daily/annual carbon lost to mortality (kgC/m2).
	double closs_mort;
	/// Daily/annual nitrogen lost to mortality (kgN/m2).
	double nloss_mort;
	/// Daily/annual carbon lost to shading mortality (kgC/m2).
	double closs_shade;
	/// Daily/annual nitrogen lost to shading mortality (kgN/m2).
	double nloss_shade;
	/// Daily drawdown on storage used to satisfy negative NPP.
	double closs_drawdown;
	/// Carbon lost to starvation today (kgC/m2).
	double closs_starve;
	/// Nitrogen lost to starvation today (kgN/m2).
	double nloss_starve;
	/// target leaf C mass (kgC/m2) 
	double cmass_leaf_limit;
	/// target root C mass (kgC/m2) 
	double cmass_root_limit;
	/// target sapwood C mass (kgC/m2)
	double cmass_sap_limit;
	/// target storage C mass (kgC/m2)
	double cmass_storage_limit;
	/// realtive leaf sink strengthy (unitless)
	double absolute_alpha_leaf;
	/// realtive root sink strengthy (unitless)
	double absolute_alpha_root;
	/// realtive sapwwod sink strengthy (unitless)
	double absolute_alpha_sap;
	/// relative storage sink strengh (unitless)
	double absolute_alpha_storage;
	/// Latest tunover day for this individual
	int last_turnover_day;
	/// daily total growth rate (kgC/m2/d)
	double cgrow;
    /// daily leaf growth rate (kgC/m2/d)
	double cgrow_leaf;
	/// daily root growth rate (kgC/m2/d)
	double cgrow_root;
	/// daily crown growth rate (kgC/m2/d)
	double cgrow_crown;
	/// daily sapwood growth rate (kgC/m2/d)
	double cgrow_sap;
	/// daily reproductive growth rate (kgC/m2/d)
	double cgrow_repr;
	/// daily storage growth (kgC/m2/day)
	double cgrow_storage;
	/// daily storage growth plus excess from N limitation (kgC/m2/day)
	double creturn;
	/// temporary object keeping track of farea of growth demand
	double farea_growth;
	/// daily growth demand / sink strength (kgC/m2/d)
	double growth_demand;
	/// development factor used for calculating growth demand.
	double fdev_growth;
	/// cmass_crown at flowering (vegetative - reproductive switch).
	double cmass_crown_flowering;
	double lai_trend_max;
	/// fraction of reproductive biomass to storage capacity.
	double frepr_cstorage;
	/// development factor used for calculating senescence.
	double fdev_senesc;
	/// senescence rate of roots (d-1)
    double senesc_root;
	/// total daily carbon senesced by dave grasses (gC/m2/day)
	double csenesc;
	/// growth stage (dormant = 0, vegetative=1, reproductive=2)
	int growth_stage;
	/// number of days in current growth stage
	int growth_stage_days;
	/// fraction of seedbank that is germinating each day
	double frac_germ;
	/// C mass of germinated seedbank (kgC/m2/d).
	double cmass_germ;
	/// Total annual C mass of germinated seedbank (kgC/m2/year).
	double acmass_germ;
	/// C exceeding storage capacity (kgC/m2/d)
	double cexcess;
	/// N exceeding storage capacity (kgN/m2/d)
	double nexcess;
	/// photosynthetic limitation due to sink activity [0-1]
	double ps_sink_lim;
	/// Daily N mass added to the indiv during the free n years (kgN/m2/d).
	double nmass_missing;
	/// leaf N mass on modelled area basis (kgN/m2)
	double nmass_leaf;
	/// root N biomass on modelled area basis (kgN/m2)
	double nmass_root;
	/// crown N mass (kgN/m2)
	double nmass_crown;
	/// reproductive tissue biomass (kgN/m2)
	double nmass_repr;
	/// sap N biomass on modelled area basis (kgN/m2)
	double nmass_sap;
	/// heart N biomass on modelled area basis (kgN/m2)
	double nmass_heart;
	/// Last year's cmass_wood value for calculation of period annual increment post mortality
	double cmass_wood_old;

	/// Daily C sent to litter (kgC/m2). This is used for accounting purposes.
	double ctolitter;
	/// Daily N sent to litter (kgN/m2). This is used for accounting purposes.
	double ntolitter;
	/// Daily C lost to high temperature (kgC/m2). This is used for accounting purposes.
	double closs_temp;
	/// Daily N lost to high temperature (kgN/m2). This is used for accounting purposes.
	double nloss_temp;
	/// Daily fraction of leaves lost to high temperature (unitless). This is used for accounting purposes.
	double loss_temp;
	/// Cumulative annual C sent to litter (kgC/m2).
	double acloss_temp;
	/// Cumulative annual N sent to litter (kgN/m2).
	double anloss_temp;
	/// Cumulative annual fraction of leaves lost to high temperature (unitless).
	double aloss_temp;
	/// Daily N taken up but then sent back to available N pool (kgN/m2).
	double ntoavail;
	/// Cumulative annual N taken up but then sent back to available N pool (kgN/m2).
	double antoavail;

	/// leaf N biomass on modelled area basis saved on first day of land use change year
	double nmass_leaf_luc;
	/// root N biomass on modelled area basis on first day of land use change year
	double nmass_root_luc;
	/// sap N biomass on modelled area basis on first day of land use change year
	double nmass_sap_luc;
	/// heart N biomass on modelled area basis on first day of land use change year
	double nmass_heart_luc;
	/// total N biomass on modelled area basis on first day of land use change year
	double nmass_tot_luc;

	/// foliar projective cover (FPC) under full leaf cover as fraction of modelled area
	double fpc;
	/// foliar projective cover (FPC) under full leaf cover as fraction of modelled area accounting for green and brown pools
	double fpc_greenbrown;
	/// foliar projective cover (FPC) this day as fraction of modelled area
	double fpc_daily;
	/// root projective cover (RPC)
	double rpc;
	/// fraction of PAR absorbed by foliage over projective area today, taking account of leaf phenological state
	double fpar;

	// Canopy layers.
	Canopy canopy;

	/// average density of individuals over patch (indiv/m2)
	double densindiv;
	/// vegetation phenological state (fraction of potential leaf cover)
	double phen;
	/// annual sum of daily fractional leaf cover
	/** Equivalent number of days with full leaf cover
	 *  (population mode only; reset on expected coldest day of year)
	 */
	double aphen;
	/// annual number of days with full leaf cover) (raingreen PFTs only; reset on 1 January)
	int aphen_raingreen;
	/// metabolic downregulation due to dormancy 
	double fphen_activity;
	/// change of fphen_activity per day
	double dphen_activity;
	// grote_dave - summer. after leaf out
	bool leaf_full_switch_summer;
	// grote_dave - autumn. scenescence after summer
	bool senescence_active;
	// grote_dave - winter. leaf senesced period
	bool leaf_senesced;
	// grote_dave - spring. leaf flushout period
	bool leaf_growth;
	// grote_dave - spring to summer. leaf out period 
	bool leaf_out;
	// grote_dave - phen value in which leaf_flushout ends and leaf_out begins
	double leaf_out_phen_maxlimit;

    /// Photosynthesis variables for subdaily calculations
	PhotosynthesisSubdaily ps_subdaily;

	/// Per-subdaily-timestep, per-canopy-layer photosynthesis inputs.
	/// Indexed as [timestep][layer]. Reinitialised each day in photosynthesis_gs().
	std::vector<std::vector<PhotosynthesisInputs>> ps_inputs_layerwise;

	/// Per-subdaily-timestep, per-canopy-layer photosynthesis outputs.
	/// Indexed as [timestep][layer]. Reinitialised each day in photosynthesis_gs().
	std::vector<std::vector<PhotosynthesisOutputs>> ps_outputs_layerwise;

	/// Per-subdaily-timestep, per-canopy-layer energy balance outputs.
	/// Indexed as [timestep][layer]. Reinitialised each day in photosynthesis_gs().
	std::vector<std::vector<EnergyBalance>> energy_balance_layerwise;

	std::vector<std::vector<PhotosynthesisDiagnostics>> ps_diagnostics_layerwise;

	/// Photosynthesis values for this individual under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;

	// Stomatal slope parameter (Medlyn et al. 2011) (kPa^0.5)
	// dynamic g1 calculated based on predawn lwp (Zhou et al. 2014)
	// if zhou is disabled, this is the same as pft.g1.
	double g1;

	// Fraction of total photosynthesis limited by Anc and Anj, respectively.
	double anc_frac;
	double anj_frac;
	double anp_frac;

	/// Adaptive parameters (i.e. changing with the environment)
	/// Adaptive N-day moving avegage of soil temperature.
	HistoricDynamic<double> soiltemp_dynamic;
	/// Dynamic root beta value used in dave version.
	double root_beta;
	// Maximum physiological up/down-regulation in and out of dormancy (fraction d-1)
	double dphen_activity_max;
	// Maximum C storage remobilisation rate (d-1)
	double max_mobil_cstorage;
	// Dormancy temperature (deg C)
	double temp_dormancy;

	/// Dynamic root distribution used in dave version. This is the fraction of
	/// roots in each soil layer.
	/// If this individual is not using dave allocation, this will be the same
	/// pft.input_rootdist (ie the old pft.rootdist, as provided in the ins file
	/// or calculted from the root beta).
	double rootdist[NSOILLAYERS_MAX];

	/// accumulated NPP over modelled area (kgC/m2/year);
	/** annual NPP following call to growth module on last day of simulation year */
	double anpp;
	double agpp;
	double aresp;
	/// yearly lai for daily allocation
	double lai_ymax;
	/// yearly fpc for daily allocation
	double fpc_ymax;
	/// yearly maximum cmass_leaf for daily allocation
	double cmass_leaf_ymax;
	/// yearly maximum nmass_leaf for daily allocation
	double nmass_leaf_ymax;
	/// yearly maximum cmass_root for daily allocation
	double cmass_root_ymax;
	/// yearly maximum nmass_root for daily allocation
	double nmass_root_ymax;
	/// Carbon storage for leaf daily allocation (kgC/m2)
	double cmass_storage;
	/// Storage flux (>0 saving; <0 storing) (kgC/m2)
	double cmass_storage_dynam;
	/// coefficient that scales growth demand using air temperature (unitless)
	double ftemp;
	/// coefficient that scales growth demand using wscal (unitless)
	double fwater;
	/// coefficient that scales growth according to actualt storage state (unitless)
	double fstorage;
	/// coefficient that scales growth according to relative day length (unitless)
	double fphoto;
	/// a coefficient used to scale cgrow by taking into account ftem, fwater and fstrage (unitless)
	double cgrow_pot;
	/// mean leaf age of the individual (y)
	double meanleafage;
	/// mean sapwood age of the individual (y)
	double meansapage;
	/// mean root age of the individual (y)
	double meanrootage;

	/// 20-day record of LAI (0-1)
	Historic<double, 20> lai_20;
	/// 20-day record of NSC pool (kgC/m2/d)
	Historic<double, 20> cmass_storage_20;
	/// 10-day record of daily reproduction pool (kgC/m2/d)
	Historic<double, 10> cmass_repr_10;
	/// 5-day record of excess C relative to NPP (cexcess/npp)
	Historic<double, 5> cexcess_npp_5;
	/// 14-day record of storage dynamics (kgC/m2/d)
	Historic<double, 14> cstorage_dynam_14;
	/// Maximum capacity of C storage pool (kgC/m2)
	double cstorage_max;
	/// Maximum capacity of N storage pool (kgN/m2)
	double nstorage_max;
	/// Amount of N withdrawn from storage today for growth (kgN/m2)
	double nstorage_drawdown;
	/// Total annual amount of N withdrawn from storage for growth (kgN/m2)
	double anstorage_drawdown;
	/// leaf C biomass on modelled area basis (kgC/m2)
	std::vector<double> cmass_leafs;
	/// leaf N biomass on modelled area basis (kgN/m2)
	std::vector<double> nmass_leafs;
	/// sum of growth of leaves
	double cmass_leaf_ygrowth;
	/// sum of growth of leaves
	double cmass_root_ygrowth;
	/// Nitrogen storage for leaf daily allocation
	double nmass_storage;
	/// Nitrogen uptake today (0 when ifstoreuptake=true) (kgN/m2).
	double nmass_uptake_no3;
	double nmass_uptake_nh4;

	/// Nitrogen retranslocated today (kgN/m2).
	double nretranslocated;
	/// Cumulative total nitrogen retranslocated this year (kgN/m2).
	double anretranslocated;

	// green leaf c mass grazed today (kgC/m2).
	double charvest_green;

	// brown leaf c mass grazed today (kgC/m2).
	double charvest_brown;

	// repr c mass grazed today (kgC/m2).
	double charvest_repr;

	// green leaf c mass grazed today (kgC/m2).
	double nharvest_green;

	// brown leaf c mass grazed today (kgC/m2).
	double nharvest_brown;

	// repr c mass grazed today (kgC/m2).
	double nharvest_repr;

	/// running mean of nscal for last 5 days for each individual daily allocation
	Historic<double, 5> nscal_5;
	/// running mean of nscal for last 5 days for each individual daily allocation
	Historic<double, 30> nscal_30;
	/// actual transpiration over projected area (mm/day)
	double aet;
	/// annual actual transpiration over projected area (mm/year)
	double aaet;
	/// fractional water uptake from each soil layer
	double fwuptake[NSOILLAYERS_MAX];
	/// leaf to root mass ratio
	double ltor;
	/// plant height (m)
	double height;
	// stem diameter (m)
	double diam;
	/// tree height increment (m)
	double d_height;
	/// tree diameter increment (m)
	double d_diam;
	/// tree annual diamter increment (m)
	double c_d_diam; 
	/// basal area increment (m2/day)
	double basalarea_inc;
	/// tree baseal area  (m2)
	double basalarea;
	/// plant crown area (m2)
	double crownarea;
	/// sapwood area (m2 sapwood/m2 ground)
	double sa;
	/// Sapwood depth (m).
	double sapwood_depth;
	/// dynamic leaf to sapwood area (unitless)
	double latosa;
	/// increment in fpc since last simulation year
	double deltafpc;
	/// bole height, i.e. height above ground of bottom of crown cylinder (m)
	/** (individual and cohort modes only) */
	double boleht;
	/// patch-level lai for this individual or cohort (function fpar)
	double lai;
	/// patch-level lai for cohort in current vertical layer (function fpar)
	double lai_layer;
	/// individual leaf area index (individual and cohort modes only)
	double lai_indiv;
	/// individual leaf area index of brown leaves (individual and cohort modes only; grasses only)
	double lai_brown_indiv;
	/// patch-level individual leaf area index (individual and cohort modes only)
	double lai_daily;
	/// daily individual leaf area index (individual and cohort modes only)
	double lai_indiv_daily;
	/// growth efficiency (NPP/leaf area) for each of the last five simulation years (kgC/m2/yr)
	Historic<double, NYEARGREFF> greff_5;
	/// increment of wood C for each of the last five simulation years (kgC/m2/yr)
	Historic<double, 10> cmass_wood_inc_5;
	/// individual/cohort age (years)
	double age;
	/// idndividual/cohort age accounting for active/dormant phases (years)
	double age_active;
	/// monthly LAI (including phenology component)
	double mlai[12];
	/// monthly maximum LAI (including phenology component)
	double mlai_max[12];
	/// FPAR assuming full leaf cover for all vegetation
	double fpar_leafon;
	/// LAI for current layer in canopy (cohort/individual mode; see function fpar)
	double lai_leafon_layer;
	/// non-water-stressed canopy conductance on FPC basis (mm/s)
	double gpterm;
	/// sub-daily version of the above variable (mm/s)
	std::vector<double> gpterms;
	/// interception associated with this individual today (patch basis)
	double intercep;

	/// accumulated mean fraction of potential leaf cover
	double phen_mean;

	/// whether individual subject to water stress
	bool wstress;
	/// non-water-stressed ratio of intercellular to ambient CO2
	double lambda;

	/// Nitrogen extinction scalar - Scalar to account for leaf nitrogen not following the optimal light extinction, but is shallower.
	double nextin;
	/// long-term storage of labile nitrogen
	double nstore_longterm;
	/// storage of labile nitrogen
	double nstore_labile;
	/// long-term storage of labile nitrogen saved on first day of land use change year
	double nstore_longterm_luc;
	/// storage of labile nitrogen saved on first day of land use change year
	double nstore_labile_luc;
	/// daily grass leaf nitrogen saved on first day of land use change year
	double nmass_leafs_luc[4];
	/// storage of root storage nitrogen saved on first day of land use change year
	double nmass_storage_luc;
	/// daily total nitrogen demand
	double ndemand;
	/// Annual total leaf nitrogen demand. Only used for trees when ifdavephotosyn is 1
	double andemand_leaf;
	/// Annual total root nitrogen demand. Only used for trees when ifdavephotosyn is 1
	double andemand_root;
	/// Annual total sapwood nitrogen demand. Only used for trees when ifdavephotosyn is 1
	double andemand_sap;
	/// Total annual nitrogen demand. Only used for trees when ifdavephotosyn is 1.
	double andemand;
	/// fraction of individual nitrogen demand available for uptake
	double fnuptake;
	/// annual nitrogen uptake (kgN/m2)
	double anuptake;
	/// daily nitrogen uptake (kgN/m2)
	double dnuptake;
	/// maximum size of nitrogen storage
	double max_n_storage;
	/// scales annual npp to maximum nitrogen storage
	double scale_n_storage;
	/// annual nitrogen limitation on vmax
	double avmaxnlim;
	/// monthly nitrogen limitation on vmax
	double mvmaxnlim[12];
	/// annual lambda limitation on photosythesis
	double alambda;
	/// monthly lambda limitation on photosythesis
	double mlambda[12];
	/// daily optimal leaf C:N ratio
	double cton_leaf_dopt;
	/// annual optimal leaf C:N ratio
	double cton_leaf_aopt;
	/// annual average leaf C:N ratio
	double cton_leaf_aavr;
	/// daily optimal leaf C:N ratio
	double cton_leaf_opt;
	/// plant mobile nitrogen status
	double cton_status;
	/// total carbon in compartments before growth
	double cmass_veg;
	/// total nitrogen in compartments before growth
	double nmass_veg;
	/// whether individual optimal leaf nitrogen is above allowed limit
	bool n_opt_isabovelim;
	/// whether individual subject to nitrogen stress
	bool nstress;
	/// daily leaf nitrogen demand. Calculated from Vmax in orig. LPJ-GUESS (kgN/m2)
	double leafndemand;
	/// daily root nitrogen demand based on leafndemand
	double rootndemand;
	/// daily crown nitrogen demand in grasses (kgN/m2)
	double crownndemand;
	/// daily reproductive tissue nitrogen demand
	double reprndemand;
	/// daily sap wood nitrogen demand based on leafndemand
	double sapndemand;
	/// daily labile nitrogen demand based on npp
	double storendemand;
	/// daily harvestable organ nitrogen demand
	double hondemand;
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

	/// The daily C lossed from leaves due to senescense, only crops.
	double daily_cmass_leafloss;
	/// The daily N lossed from leaves due to senescense, only crops.
	double daily_nmass_leafloss;
	/// The daily C lossed from roots due to senescense, only crops.
	double daily_cmass_rootloss;
	/// The daily N lossed from roots due to senescense, only crops.
	double daily_nmass_rootloss;

	/// daily leaf turnover flux (kgC/m^2/day) (DAVE)
	double cturnover_leaf;
	/// daily root turnover flux (kgC/m^2/day) (DAVE)
	double cturnover_root;
	/// daily sapwood turnover flux (kgC/m^2/day) (DAVE)
	double cturnover_sap;
	/// daily repr turnover flux (kgC/m^2/day) (DAVE)
	double cturnover_repr;
	/// daily crown turnover flux (kgC/m^2/day) (DAVE)
	double cturnover_crown;
	/// daily leaf turnover flux (kgN/m^2/day) (DAVE)
	double nturnover_leaf;
	/// daily root turnover flux (kgN/m^2/day) (DAVE)
	double nturnover_root;
	/// daily sapwood turnover flux (kgN/m^2/day) (DAVE)
	double nturnover_sap;
	/// daily repr turnover flux (kgN/m^2/day) (DAVE)
	double nturnover_repr;
	/// daily crown turnover flux (kgN/m^2/day) (DAVE)
	double nturnover_crown;
	/// daily leaf turnover as a proportion of leaf C biomass  (DAVE)
	double pturnover_leaf;

	/// Number of days with non-negligible phenology this year
	int nday_leafon;
	// Whether this individual is truly alive.
	/** Set to false for first year after the Individual object is created, then true. */
	bool alive;

	/// GPP today (kgC/m2/day)
	/// Note: in the trunk version, this is gross photosynthesis minus leaf respiration.
	/// In the dave version, this is gross photosynthesis only.
	double dgpp;

	/// NPP this day (kgC/m2/day)
	double dnpp;
	/// GPP this day without any stresses
	double dgpp_no;
	/// GPP this day with only nutrient stresses
	double dgpp_ns;
	/// Montly GPP without any stresses
	double mgpp_no[12];
	/// Annual GPP without any stresses
	double agpp_no;
	/// Annual GPP with only N stress
	double agpp_ns;

	Historic<double, 15> gpp_15;
	Historic<double, 15> gpptrend_15;

	// Maintenance respiration today (kgC/m2)
	double resp_maintenance;
	// Growth respiration today (kgC/m2)
	double resp_growth;

	// Leaf respiration today (kgC/m2)
	double resp_leaf;

	// Culm respiration today (kgC/m2)
	double resp_culm;

	// Sapwood respiration today (kgC/m2)
	double resp_sap;

	// Root respiration today (kgC/m2)
	double resp_root;

	// Reproductive tissue respiration today (kgC/m2)
	double resp_repr;

	// Crown respiration today (kgC/m2)
	double resp_crown;

	// bvoc

	/// isoprene production (mg C m-2 d-1)
	double iso;
	/// monoterpene production (mg C m-2 d-1) per monoteprene species
	double mon[NMTCOMPOUNDS];
	/// monoterpene storage pool (mg C m-2) per monoterpene species
	double monstor[NMTCOMPOUNDS];
	/// isoprene seasonality factor (-)
	double fvocseas;

	/// Pointer to struct with crop-specific data
	cropindiv_struct *cropindiv;

	/// cutting intensity (fraction cut)
	double man_strength;

	/// Daily allocation
	double g_fac;
	double gpp_fac;
	double sb_fac;
	double wscal_fac;
	double tscal_fac;
	double s_fac;
	double frost_sen;
	double wscal_sen;
	double tscal_sen;

	bool stomata_open; // flag on whether stomata are open (at all) - if false, stomata are fully closed
	
	double weighted_psi; // soil water potential weighted by root fraction in each layer

	/// hydraulic resistance between soil and roots (mmol/m/s/MPa)
	double Rs_r[NSOILLAYERS_MAX];

	// predawn leaf water potential calculated using hydraulics module
	double predawn_lwp;
	// leaf water potential at subdaily timescale (MPa)
	std::vector<double> lwp;
	// leaf water potential (MPa) from previous timestep (for plant water storage)
	double lwp_previous;
	// loss of hydraulic conductance at subdaily timescale (mol/m2/s/MPa)
	std::vector<double> k_loss;

	// indicator of whether hydraulic failure has occurred
	bool hydraulic_failure;

	// individual-level transpiration, mm (includes cuticular conductance when roots are disconnected)
	std::vector<double> transpiration;

	private:
	/// Adaptive N-day moving avegage of wscal.
	HistoricDynamic<double> wscal_dynamic;
	// wscal calculated on an individual basis using root conductance
	double wscal;
	Historic<double, 5> wscal_5;
	Historic<double, 10> wscal_10;
	Historic<double, 30> wscal_30;
	double wscal_mean_year;

	public:

	double alpha_leaf;
	double alpha_root;
	double alpha_sap;
	double alpha_storage;

	Historic<double, 30> alpha_leaf_trend;
	Historic<double, 30> alpha_root_trend;
	Historic<double, 30> alpha_sap_trend;
	Historic<double, 30> alpha_storage_trend;

	double cdeficit_leaf;
	double cdeficit_root;
	double cdeficit_sap;
	double cdeficit_storage;

	// Fraction of individuals within this cohort killed this year due to age-
	// related mortality.
	double mort_age;

	// Fraction of individuals within this cohort killed this year due to
	// growth efficiency-related mortality.
	double mort_greff;

	// Fraction of individuals within this cohort killed this year due to
	// carbon starvation.
	double mort_cstarve;

	// Fraction of individuals within this cohort killed today due to
	// carbon starvation.
	double mort_dcstarve;
	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v);
	~Individual();

	/// Access functions for cropindiv:
	cropindiv_struct* get_cropindiv() const;
	cropindiv_struct* set_cropindiv();

	void serialize(ArchiveStream& arch);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPFTFluxType flux_type, double value);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPatchFluxType flux_type, double value);

	/// Whether an individual is either a true crop or a cover crop grass
	inline bool istruecrop_or_intercropgrass() const {
		return (pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || cropindiv->isintercropgrass));
	}

	// Flush reproductive pool to litter
	void flush_repr_litter_DAVE(double flush_fraction);

	/// Whether harvest and turnover is done on actual C and N on harvest or turnover day, which can occur any day of the year.
	bool has_daily_turnover() const;

	/// Distribute stored N between plant compartments to prevent unrealistic C:N ratios in litter when killed
	void partition_nstore(double& NH4_mass);

	/// Whether resetting of grs_cmass and turnover (if has_daily_turnover() returns true) of continuous grass is to be done this day.
	/** This should occur at the very end of the growing period */
	bool is_turnover_day() const;

	/// Get brown leaf biomass for this individual. (kgC/m2)
	/// This is 0 for trees, shrubs, and mosses.
	double get_cmass_leaf_brown() const;

	double get_cmass_leaf() const;

	double get_nmass_leaf() const;

	/// @brief Set the water stress value in this individual. This will update all moving-averages of water stress.
	/// @param wscal The water stress value (1=no stress, 0=max stress).
	void set_wscal(double wscal);
	// Return appropriate wscal (based on root resistance if using hydraulics, otherwise based on wcont)
	double get_wscal() const;
	const Historic<double, 5> &get_wscal_5() const;
	const Historic<double, 10> &get_wscal_10() const;
	const Historic<double, 30> &get_wscal_30() const;

	/// Reduce current biomass due to mortality and/or fire
	/** The removed biomass is put into litter pools and/or goes to fire fluxes.
	 *
	 *  \param mortality      fraction of Individual's biomass killed due to
	 *                        mortality (including fire)
	 *  \param mortality_fire fraction of Individual's biomass killed due to
	 *                        fire only
	 */
	void reduce_biomass(double mortality, double mortality_fire, double &closs, double &nloss);

	/// A version of the above reduce_biomass for the use with BLAZE
	void blaze_reduce_biomass(Patch& patch, double frac_survive, double &closs, double &nloss);

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

	// Total above-ground biomass (kg/m2).
	double aboveground_biomass() const;

	/// Total carbon biomass
	double ccont(double scale_indiv = 1.0, bool luc = false, bool ignore_alive = false) const;
	/// Total nitrogen biomass
	double ncont(double scale_indiv = 1.0, bool luc = false, bool ignore_alive = false) const;

	/// Whether grass growth is uninterrupted by crop growth.
	bool continous_grass() const;

	/// Checks whether any grs_cmass is negative, in which case it is zeroed and fluxes are corrected (only cropland).
	double check_C_mass();
	/// Checks whether any nmass is negative, in which case it is zeroed and fluxes are corrected (only cropland).
	double check_N_mass();

	/// Save cmass-values on first day of the year of land cover change in expanding stands
	void save_cmass_luc();
	/// Save nmass-values on first day of the year of land cover change in expanding stands
	void save_nmass_luc();

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

	/// Current crown C:N ratio (DAVE grasses only)
	double cton_crown() const;

	/// Current sap C:N ratio
	double cton_sap() const;

	/// Current heartwood C:N ratio.
	double cton_heart() const;

	/// Gets the individual's Patchpft
	Patchpft& patchpft() const;

	/// returns SLA value of indiv
	double sla() const;

	/// Transfers the individual's biomass (C and N) to litter and harvest pools/fluxes
	/**
	 *  \param harvest Set to true if some of the biomass should be harvested,
	 *                 harvest will be done according to the PFT's harvest efficiency
	 *                 and residue outtake.
	 */
	void kill(bool harvest = false);

	/// Annual mean wscal - water stress parameter (0-1 range; 1 = minimum stress)
	/** Value only valid at end of year, after call to canopy_exchange().
	 *
	 *  Currently, all Individuals belonging to a Patchpft share the same water stress.
	 */
	double wscal_mean() const;

	/// Gets the individual's daily cmass_leaf value
	double cmass_leaf_today() const;
	/// Gets the individual's daily cmass_root value
	double cmass_root_today() const;

	/// Return crown cmass for DAVE grasses, 0 for all other PFTs (kgC/m2).
	double get_cmass_crown() const;

	/// Return live reproductive cmass for DAVE PFTs, 0 for all other PFTs (kgC/m2).
	double get_cmass_repr() const;

	/// Gets the individual's daily LAI value (patch-level)
	/** Based on total leaf area for whatever the individual represents
	 *  (individual, cohort, population), over the whole patch.
	 */
	double lai_today() const;

	/// Gets the individual's daily LAI value (individual-level)
	/** Based on the leaf area for the average individual and
	 *  the average individual's crown area.
	 */
	double lai_indiv_today() const;

	/// Gets the Nitrogen limited LAI, Eq. 8 Olin 2015
	double lai_nitrogen_today() const;

	/// Gets the individual's daily fpc value
	double fpc_today() const;

	/// Gets the individual's daily fpc value accounting for green and brown leaves
	double fpc_greenbrown_today() const;

	/// Gets the growingseason status for crop individual. Non-crop individuals always return true.
	bool growingseason() const;

	/// The N demand of the storage, only used for PNV.
	double ndemand_storage(double cton_leaf_opt);

	/// Reset daily variables.
	void daily_accounting();

	/// Get the capacity of the dynamic-length moving average of water stress.
	size_t get_wscal_dynamic_capacity() const;

	/// Set the capacity of the dynamic-length moving average of water stress.
	/// @param capacity The new window length of the moving average.
	void set_wscal_dynamic_capacity(size_t capacity);

	/// Get the moving average of water stress.
	double get_wscal_moving_average() const;

	/// Ensure that canopy layers are internally consistent.
	void validate_canopy() const;

	/// Get the fraction of water uptake from a given soil layer.
	double get_fwuptake(size_t layer) const;
private:
	// Reset annual variables. Called from daily_accounting() on the first day
	// of the year.
	void annual_accounting();
	bool use_indiv_wscal() const;
};


/// The vegetation in a patch - a list of individuals
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Individual objects. A single Vegetation object is defined for each patch. A
 *  reference to the parent Patch object (defined below) is included as a member
 *  variable.
 *
 *  Functionality is inherited from the ListArray_idin1 template type in the GUTIL
 *  Library. Sequential Individual objects can be referenced as array elements by id,
 *  or by iteration through the linked list:
 *
 *    Vegetation vegetation
 *    ...
 *    vegetation.firstobj();
 *    while (vegetation.isobj) {
 *      Individual& thisindiv=vegetation.getobj();
 *      // query or modify object thisindiv here
 *      vegetation.nextobj();
 *    }
 */
class Vegetation : public ListArray_idin2<Individual,Pft,Vegetation>, public Serializable {

public:
	// MEMBER VARIABLES

	/// reference to parent Patch object
	Patch& patch;

	// MEMBER FUNCTIONS

	/// constructor (initialises member variable patch)
	Vegetation(Patch& p):patch(p) {};

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

	soilinputtype soil_input;

	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)
	double perc_base;
	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996) for the evaporation soil layer
	double perc_base_evap;
	/// exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice 1996)
	double perc_exp;

	/// thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0;
	/// thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15;
	/// thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100;

	/// Available water holding capacity of soil layers (mm)
	double awc_min;
	/// Wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double wp_min;
	/// Saturation point. Cosby et al 1984
	double wsats_min;

	/// Organic soil properties
	
	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)
	double perc_base_org;

	/// thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0_org;
	/// thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15_org;
	/// thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100_org;

	// Sand, silt and clay fractions, should always add up to 1.
	/// fraction of soil that is sand
	std::vector<double> sand_frac;
	/// fraction of soil that is clay
	std::vector<double> clay_frac;
	/// fraction of soil that is coarse fragmnets
	std::vector<double> coarse_frac;
	/// pH
	std::vector<double> pH;
	/// carbon content (kg C m-2)
	std::vector<double> soilC;

	/// Soilcode, 0 to 8
	int soilcode;

	/// Volumetric fraction of organic material (m3 m-3) (Hillel, 1998)
	std::vector<double> organic_frac;

	/// Bottom depth of texture data
	std::vector<double> text_depths;

	// MEMBER FUNCTIONS

public:

	/// Constructor: initialises certain member variables
	Soiltype() {
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
		ntoc = 0.0;
		cdec = 0.0; 
		ndec = 0.0; 
		cinc = 0.0;
		ninc = 0.0;
		respsum = 0.0;
		respsum_surf = 0.0;
		ligcfrac = 0.0;
		delta_cmass = 0.0;
		delta_nmass = 0.0;
		fracremain = 0.0;
		fireresist = 0.0;
		cmass_init = 0.0;
		multi = 1.0;
		sign_trend = 1;
		same_trend = 0;
		changed_trend = 0;
		litterme = 0.0;
		fireresist = 0.0;
		ntoc = 0.0;

		for (int m = 0; m < 12; m++) {
			mfracremain_mean[m] = 0.0;
		}
	}

	/// Destructor
	~Sompool() {
		// Performs memory deallocation, closing of files or other "cleanup" functions.
	}

	/// C mass in pool kgC/m2
	double cmass;
	/// Nitrogen mass in pool kgN/m2
	double nmass;
	/// (potential) decrease in C following decomposition today (kgC/m2)
	double cdec;
	/// (potential) decrease in nitrogen following decomposition today (kgN/m2)
	double ndec;
	/// (potential) increase in C following decomposition today (kgC/m2)
	double cinc;
	/// (potential) increase in nitrogen following decomposition today (kgN/m2)
	double ninc;
	/// (potential) microbial respiration following decomposition today (kgC/m2)
	double respsum;
	/// (potential) surface microbial respiration following decomposition today (kgC/m2)
	double respsum_surf;
	/// daily change in carbon and nitrogen
	double delta_cmass, delta_nmass;
	/// lignin fractions
	double ligcfrac;
	/// fraction of pool remaining after decomposition
	double fracremain;
	/// nitrogen to carbon ratio
	double ntoc;
	/// Initial C mass in pool (kgC/m2) used for fast spinup
	double cmass_init;
	/// Multiplier for passive SOM C mass trend
	double multi;
	/// Sign of passive SOM C mass trend
	int sign_trend;
	/// How many times the passive SOM C mass trend has have the same sign in a raw
	int same_trend;
	/// How many times the passive SOM C mass trend has have changed sign in a raw
	int changed_trend;

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

	/// Destructor
	~LitterSolveSOM() {
		// Performs memory deallocation, closing of files or other "cleanup" functions.
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
	/// Carbon litter
	double clitter[NSOMPOOL];

	/// Nitrogen litter
	double nlitter[NSOMPOOL];
};

/// Soillayer stores state variables for.
/**
*
*/
class Soillayer : public Serializable {

	// MEMBER VARIABLES

public:

	/// reference to parent soil object
	Soil* soil;

	/// Depth
	/// Soil layer thickness (mm)
	double thickness;
	/// Soil layer middle depth (mm)
	double depth_mid;
	/// Soil layer bottom depth (mm)
	double depth_bot;

	// Sand, silt and clay fractions, should always add up to 1.
	/// fraction of mineral soil that is sand
	double sand_frac;
	/// fraction of mineral soil that is clay
	double clay_frac;
	/// fraction of mineral soil that is silt
	double silt_frac;
	/// fraction of soil that is coarse fragments
	double coarse_frac;
	/// pH
	double pH;

	/// Water
	/// Water mass in soil layer (kg m-2)
	double wmass;
	/// Monthly water mass in soil layer (kg m-2)
	double mwmass[12];
	/// Suppercooled water mass in soil layer (kg m-2)
	double wmass_sc;
	/// Monthly water content as fraction of available water holding capacity
	double mwcont[12];
	/// Mean water content above current soil layer for last month
	double mwcontabove;
	/// Average wcont over the growing season
	double awcont;
	/// Saturation point of soil layer (mm) Cosby et al 1984
	double wsats;
	/// Available water holding capacity of soil layer (mm)
	double awc;
	/// Wilting point of soil layer (mm) Cosby et al 1984
	double wp;

	/// Ice
	/// Ice mass in soil layer (kg m-2)
	double imass;
	/// Monthly ice mass in soil layer (kg m-2)
	double mimass[12];

	/// Matrix potential (mm)
	double Psi;
	/// Derivative in matrix potential w.r.t. volumetric liquid water content (mm)
	double dPsidw;
	/// Soil matric potential at saturation (mm)
	double Psi_max;

	/// Soil moisture conductivity (mm s-1)
	double Kw;
	/// Derivative in hydraulic conductivity w.r.t relative saturation (mm s-1)
	double dKwdw;
	/// Soil moisture conductivity at saturation (mm s-1)
	double Kw_max;
	/// Exponential in Campbell (1974) parameterization for soil water retention curves
	double b;
	/// exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice 1996)
	double perc_exp;
	/// Percolation mass from this layer
	double perc_mass;
	/// Coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)
	double perc_base;
	/// Coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996) for mineral soils
	double perc_base_min;
	/// Coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996) for organic soils
	double perc_base_org;
	/// Lateral drainage runoff
	double runoff_drain;
	/// Monthly lateral drainage runoff
	double mdrain[12];
	/// Baseflow runoff
	double runoff_baseflow;

	/// Mineral available water holding fraction of soil layer (-)
	double awc_min;
	/// Mineral wilting point of soil layer (-)
	double wp_min;
	/// Mineral saturation point (-)
	double wsats_min;
	/// Mineral exponential in Campbell (1974) parameterization for soil water retention curves 
	double b_min;
	/// Mineral hydraulic conductivity at saturation (mm s-1)
	double Kw_max_min;
	/// Saturated mineral matter water suction (matrix potential) (mm)
	double Psi_max_min;

	/// Coarse fragment (rock) available water holding fraction of soil layer (-)
	double awc_rock;
	/// Coarse fragment (rock) wilting point of soil layer (-)
	double wp_rock;
	/// Coarse fragment (rock) saturation point (-)
	double wsats_rock;
	/// Coarse fragment (rock) exponential in Campbell (1974) parameterization for soil water retention curves 
	double b_rock;
	/// Coarse fragment (rock) hydraulic conductivity at saturation (mm s-1)
	double Kw_max_rock;
	/// Coarse fragment (rock) mineral matter water suction (matrix potential) (mm)
	double Psi_max_rock;
	
	/// Volume
	/// Coarse root volume (mm)
	double root_volume;
	/// Mineral volume (mm)
	double mineral_volume;
	/// Coarse fragment volume (mm)
	double coarse_volume;
	/// Organic matter volume (mm)
	double organic_volume;
	/// Effective volume when accounting for space occupied by roots (mm)
	double effective_volume;
	// Soil specific surface area (m2 g−1)
	double ssa;

	/// Bulk density (g cm-3)
	double bulk_density;
	
	/// Organic influence on water
	// Organic carbon density (kg C m-3)
	double rho_org;
	/// Organic available water holding fraction of soil layer (-)
	double awc_org;
	/// Organic wilting point of soil layer (-)
	double wp_org;
	/// Organic saturation point (-)
	double wsats_org;
	/// Organic exponential in Campbell (1974) parameterization for soil water retention curves 
	double b_org;
	/// Organic soil matric potential at saturation 
	double Psi_max_org;
	/// Hydraulic conductivity at saturation for organic matter (m s-1)
	double Kw_max_org;
	/// Volumetric fraction of peat material (m3 m-3)
	double peat_frac;
	/// if soil layer is in the acrotelm
	bool acrotelm;
	/// if soil layer is in the catotelm
	bool catotelm;
	/// Tiller area
	double tiller_area;

	/// Temperature
	/// Thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0;
	/// Thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15;
	/// Thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100;
	/// Thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0_min;
	/// Thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15_min;
	/// Thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100_min;
	/// Thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0_org;
	/// Thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15_org;
	/// Thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100_org;
	/// Thermal conductivities (J m-1 K-1 s-1)
	double Kh;
	/// Heat capacities (J m-3 K-1)
	double Ch;
	/// Heat capacity of mineral soil [J m-3 K-1]
	double Cp_min;
	/// Soil layer temperature (deg C)
	double temp;
	/// Soil layer temperature from previous phase change (deg C)
	double temp_old;
	/// Record the monthly average soil temp at SOILTEMPOUT layers [deg C]
	double mtemp[12];
	/** respiration response to today's soil temperature at 0.25 m depth
	*  incorporating damping of Q10 due to temperature acclimation (Lloyd & Taylor 1994)
	*/
	double gtemp;
	// Parameters used by function soiltemp and updated monthly
	double alag, exp_alag;

	/// Gases
	/// Dissolved CO2 concentration in each layer [g CO2-C layer-1 d-1]
	double CO2;
	/// Dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4;
	/// CH4 which bubbles out [g CH4-C layer-1]
	double CH4_ebull_ind;
	/// Gaseous CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4_gas;
	/// Dissolved CH4 concentration in each layer [g CH4-C layer-1 d-1]
	double CH4_diss;
	/// Dissolved O2 concentration in each layer [mol O2 layer-1 d-1]
	double O2;
	/// Todays heterotrophic respiration [kg C m-2]
	double respsum;
	/// Annual total respiration [kg C m-2]
	double arespsum;

	/// Annual fine litter C input
	double acmassfinput;
	/// Annual coarse litter C input
	double acmasscinput;

	//////////////////////////////////////////////////////////////////////////////////
	// Old Soil Organic Matter (SOM) scheme

	/// SOM pool with c. 1000 yr turnover (kgC/m2)
	double cpool_slow;
	/// SOM pool with c. 33 yr turnover (kgC/m2)
	double cpool_fast;
	/// mean annual litter decomposition (kgC/m2/yr)
	double decomp_litter_mean;
	/// mean value of decay constant for fast SOM fraction
	double k_soilfast_mean;
	/// mean value of decay constant for slow SOM fraction
	double k_soilslow_mean;


	//////////////////////////////////////////////////////////////////////////////////
	// CENTURY SOM pools and other variables

	Sompool sompool[NSOMPOOL];

	/// Decay rate modifier depending on depth profile (Koven et al 2013)
	double depth_mod;
	/// Soil organic C microorganism collocation factor
	double a;
	/// Annual average temperature rate modifier
	double atemp_mod;
	/// Annual average moisture rate modifier
	double amoist_mod;
	/// Annual average total rate modifier
	double atot_mod;
	
	/// Fraction of decayed organic matter leached each day;
	double orgleachfrac;
	/// Annual leaching from available nitrogen pool
	double aminleach;
	/// Annual leaching of organics from active nitrogen pool
	double aorgNleach;
	/// Annual leaching of organics nitrogen from carbon pool
	double aorgCleach;
	// Soil labile carbon availability daily (kgC/m2/day)
	double labile_carbon;
	// Age of organic matter
	double socage;

	/// Annual sum of nitrogen mineralisation
	double anmin;
	/// Annual sum of nitrogen immobilisation
	double animmob;
	/// Monthly sum of nitrogen net mineralisation
	double mnetmin[12];
	/// Annual C decomposition without and with nutrient limitation
	double acdec[2];
	/// Annual C mineralisation without and with nutrient limitation
	double acinc[2];

	/// Soil NH4 mass in pool (kgN/m2)
	double NH4_mass;
	/// Soil NO3 mass in pool (kgN/m2)
	double NO3_mass; 
	/// Soil NO2 mass in pool (kgN/m2)
	double NO2_mass;
	/// Soil NO mass in pool (kgN/m2)
	double NO_mass;
	/// Soil NO mass in pool (kgN/m2)
	double N2O_mass;
	/// Soil N2 mass in pool (kgN/m2)
	double N2_mass;

	// Variables for fast spinup of SOM pools

	/// monthly fraction of organic carbon/nitrogen leached
	double morgleach_mean[12];
	/// monthly fraction of available mineral nitrogen leached
	double mminleach_mean[12];

	// Solving Century SOM pools

	/// Cumulative litter pools for one year.
	LitterSolveSOM litterSolveSOM;

	std::vector<LitterSolveSOM> solvesom;

	double k_actual;

	// Output variables.
	double nvol;
	double anvol;
	double nloss_gas;
	double anloss_gas;
	// MEMBER FUNCTIONS

public:

	/// Constructs a Soillayer
	/**  \param d         The soil layer depth
	*/
	Soillayer() : soil(NULL) {};

	Soillayer(Soil* soilcol, double d, double dp, bool acrotelm, bool catotelm);

	// Destructor
	~Soillayer() {

		// Performs memory deallocation, closing of files or other "cleanup" functions.
	}

	// Mineral exponential in Campbell (1974) parameterization for soil 
	// water retention curves from Table 4, Cosby (1984)
	double mineral_exp(double clay_frac, double sand_frac);
	/// Set soil layer properties
	void set_properties();
	/// Initiate soil layer properties
	void initiate_properties();
	/// Update soil layer properties
	void update_properties(double& wmass_move, double& imass_move, bool must_set = false);
	/// C content of soil layer
	double ccont(grouppooltype group = GALL) const;
	/// N content of soil layer
	double ncont(grouppooltype group = GALL) const;
	/// Calculate organic matter density
	double calc_organic_density(grouppooltype pools = GALL);
	/// Calculate clay fraction of this layer
	double fclay(mincalctype type = MINEFF);
	/// Calculate sand fraction of this layer
	double fsand(mincalctype type = MINEFF);
	/// Calculate silt fraction of this layer
	double fsilt(mincalctype type = MINEFF);
	/// Freeze full layer
	void freeze_lyr();
	/// Phase change
	void phase_change();
	/// Water content of soil layer
	double wcont() const;
	/// Water content of soil layer
	double wcont(double wmass_added) const;
	/// Field capacity of soil layer
	double fc();
	/// Fraction air in soil layer
	double fracair();
	/// available N content of soil layer
	double nmass_avail(n_pref_type pref = NO, bool water_limited = false);
	void nmass_subtract(double nmass, n_pref_type pref = NO);
	void nmass_inc(double nmass, n_pref_type pref = NO);
	/// Return the amount of soil available nitrogen where SOM C:N ratio 
	/// reach their minimum (nitrogen saturation) [kg N m-3]
	double nmass_saturation(double NMASS_SAT);
	/// Water filled pore space ([0-1])
	double wfps();

	/// Maximum water capacity of soil layer
	double wmass_max();	
	/// Set target nutrient to carbon ratios for SOM pools
	void setntoc(double fac, double& ntoc, double cton_max, double cton_min, double fmin, double fmax);
	/// Calculates CENTURY instantaneous decay rates
	void decayrates_century(double tillage_fact, bool savesomcent, double temp_surf);
	/// Water modifier for decomposition (Friend et al 1997)
	double moisture_modifier();
	/// Water modifier for decomposition (Yan et al 20189
	double moisture_modifier_yan();

	void update_temp(double add_temp, double add_mass, double add_Ch);

	void update_matric_potential();

	void update_thermal_properties();

	void update_socage(double cmass, double age);

	double org_frac();

	/// Get soil water potential in this layer, in MPa.
	double psi_mpa() const;

	void serialize(ArchiveStream& arch);

	void daily_accounting();
private:
	void annual_accounting();
};

/// Container for soil padlayer parameters
struct Padlayer : public Serializable {

	Padlayer() {
		temp = 0.0;
		thickness = 0.0;
		// Bedrock parameters from Chadburn et al. (2015)
		// Thermal conductivity [W m-1 K-1 s-1]
		Kh = 8.6;
		// Heat capacity [J m-3 K-1]
		Ch = 2100000.0;

		for (int mth = 0; mth < 12; mth++) {
			mtemp[mth] = 0.0;
		}
	}

	/// Soil padding layer temperatures [deg C]
	double temp;
	/// Record the monthly average padding temperature [deg C]
	double mtemp[12];
	/// Soil padding layer thickness
	double thickness;
	/// Soil padding layer thermal conductivities (J m-1 K-1 s-1)
	double Kh;
	/// Soil padding layer heat capacities (J m-3 K-1)
	double Ch;
	/// Update thermal properties for padding layers
	void update_thermal_properties(double Kh_last_lyr = 0.0, double Ch_last_lyr = 0.0);
	/// Serialize padding layer variables
	void serialize(ArchiveStream& arch);
};

class Padding : public Serializable {

	// MEMBER VARIABLES

public:

	std::vector<Padlayer> padlayers;

	// MEMBER FUNCTIONS

public:

	/// Constructs an empty result
	Padding() {
		padlayers.resize(NPADLAYERS);
	}

	/// Set padding layers depth
	void set_depths(double depth_last_lyr);
	/// Update thermal properties for padding layers
	void update_thermal_properties(double Kh_last_lyr, double Ch_last_lyr);
	/// Serialize padding layer variables
	void serialize(ArchiveStream& arch);
};

struct Snowlayer : public Serializable {

	Snowlayer() {
		clear();
	}

	Snowlayer(double init_ice, double init_water, double init_dens, double init_temp) {
		clear();

		imass = init_ice;
		wmass = init_water;
		dens = init_dens;
		temp = init_temp;
	}

	/// water mass in snow layer (kg m-2)
	double wmass;

	/// ice mass in snow layer (kg m-2)
	double imass;

	/// Snow density per layer (kg m-3)
	double dens;

	/// Snow water holding capacity per layer
	double whc;

	/// Snow layer thickness [mm]
	double thickness;

	/// Snow temperature [C]
	double temp;

	/// Thermal conductivity of snow [J m-1 K-1 s-1]
	double Kh;

	/// Heat capacity of snow [J m-3 K-1]
	double Ch;

	void clear();

	void update_temp(double add_temp, double add_mass, double add_Ch);

	/// Update snow thermal properties
	void update_thermal_properties();

	/// Update snow water holding capacity
	void update_whc();

	/// Phase change
	void phase_change(bool toplayer = false);

	void serialize(ArchiveStream& arch);
};

/// This struct contains snow properties.
class Snow : public Serializable {

	// MEMBER VARIABLES

public:

	/// vector with snow layers
	std::vector<Snowlayer> snowlayers;
	/// number of snow layers
	int nsnowlayers;
	/// NH4 mass (kgN/m2)
	double NH4_mass;
	/// NO3 mass (kgN/m2)
	double NO3_mass;
	/// Monthly snow depth (average) [mm]
	double mdepth[12];
	/// Monthly snow cover (average) [-]
	double mcover[12];
	/// Number of snow days during the year
	int snowdays;
	/// Monthly snow layer temperature [C]
	double mtemp[NSNOWLAYERS_MAX][12];
	/// Excess snow depth fractions
	double snow_excess_frac[NSNOWLAYERS_MAX];

	// MEMBER FUNCTIONS

public:

	/// Constructs an empty result
	Snow() {
		init_state();
	}

	void init_state();

	/// Redistributing snow layers after melt or compaction
	void potential_depth(double snowpack_depth, double pot_ice[NSNOWLAYERS_MAX], int& pot_nsnowlayers);

	/// Update snow water holding capacity
	void update_whc();

	/// Reduce snowpack to its maximum depth, excess goes to runoff
	void reduce_snowpack(Patch& patch, double depth_reduction);

	/// Redistribute snow when a non-overfull layer increases/decreases due to melt or compaction
	void redist_layers(double snowpack_depth, double snowfall_depth = 0.0, 
		double snowfall_imass = 0.0, double snowfall_temp = 0.0, double snowfall_Ch = 0.0);

	/// Rain melts snowpack
	void rainwater_melt(double airT, double& prec, double& rain_melt, double& rain_melt_temp);

	/// Update snow fractions
	void phase_change();

	/// Soak snowpack with water
	void soak_snowpack(double& wmass_ssw, double& imass_ssw, double temp_ssw);

	/// Snow water over water holding capacity
	double wmass_excess();

	/// Snow compaction
	void compaction(Patch& patch, double& rain_melt, double& rain_melt_temp);

	/// Update snow thermal properties for each snow layer
	void update_thermal_properties();

	/// Balance ice and water in snow layers after melt
	void balance_snow(Patch& patch, double& rain_melt, double& rain_melt_temp);

	/// Mean snow depth for the three months with largest snowpack
	double snowpack_avrmax3mth();

	/// Returns total snow depth [mm]
	double depth();

	/// Returns nmass available [kgN/m2]
	double nmass_avail();

	double ccont() const;

	double ncont() const;

	/// Returns total snowpack (water and ice) [mm]
	double water_and_ice_content() const;

	/// Returns total snowpack water [mm]
	double wmass();

	/// Returns total snowpack ice [mm]
	double imass();

	void serialize(ArchiveStream& arch);
};

/// Container for surface standing water and ice
struct SSW : public Serializable {

	/// Surface standing water (kg m-2)
	double wmass;
	/// Surface standing ice (kg m-2)
	double imass;
	/// Surface standing water/ice temperature (deg C)
	double temp;
	/// Surface standing water/ice thermal conductivities (J m-1 K-1 s-1)
	double Kh;
	/// Surface standing water/ice heat capacities (J m-3 K-1)
	double Ch;
	/// Monthly surface standing water and ice [mm]
	double mwimass[12];
	/// Fraction of surface covered by surface standing water and ice
	double wmass_frac;
	// Fraction of the area that is inundated threshold below which no 
	// single connected inundated area spans the grid cell (-)
	double wmass_frac_thresh;
	// Surface water threshold below which no single connected  
	// inundated area spans the grid cell (mm)
	double wmass_thresh;
	/// Microtopography distribution (mm) (Ekcic et al. 2019; eqn 2)
	double slope_sigma;
	double sigma_sqrt2;
	double sigma_sqrt2pi;
	double sigma2;

	SSW() {
		wmass = 0.0;
		imass = 0.0;
		temp = 0.0;
		Kh = 0.0;
		Ch = 0.0;
		wmass_frac = 0.0;
		wmass_frac_thresh = 0.0;
		wmass_thresh = 0.0;
		slope_sigma = 0.0;
		sigma_sqrt2 = 0.0;
		sigma_sqrt2pi = 0.0;
		sigma2 = 0.0;

		for (int mth = 0; mth < 12; mth++) {
			mwimass[mth] = 0.0;
		}
	}

	/// Water and ice content of surface standing water
	double water_and_ice_content() const;
	/// Update surface standing water temperatures due to water fluxes
	void update_temp(double temp_add, double wmass_add, double imass_add = 0.0);
	/// Update thermal properties for surface standing water
	void update_thermal_properties();
	/// Calculate surface standing water and ice
	void phase_change();
	/// Serialize surface standing water variables
	void serialize(ArchiveStream& arch);
};

/// Soil stores state variables for soillayers and the snow pack.
/** Initialised by a call to initdrivers. One Soil object is defined for each patch.
*  A reference to the parent Patch object (defined below) is included as a member
*  variable. Soil static parameters are stored as objects of class Soiltype, of which
*  there is one for each grid cell. A reference to the Soiltype object holding the
*  static parameters for this soil is included as a member variable. Soil holds the
*  object array of soil layers
*/
class Soil : public Serializable {

	// MEMBER VARIABLES

public:

	/// reference to parent Patch object
	Patch& patch;
	/// reference to Soiltype object holding static parameters for this soil
	Soiltype& soiltype;
	/// Soil depth
	double depth;
	/// vector with soil layers
	std::vector<Soillayer> soillayers;
	/// array with soil layer thicknesses
	double soillayerthicks[NSOILLAYERS_MAX];
	/// number of soil layers
	int nsoillayers;
	/// Annual surface fine litter input [kg C m-2 yr-1]
	double asurfflitterinput;
	/// Annual surface coarse litter input [kg C m-2 yr-1]
	double asurfclitterinput;
	/// Annual surface litter respiration [kg C m-2 yr-1]
	double asurflitterresp;
	/// vector with SOM pool litter input
	double csom_litter_input[NSOMPOOL];
	/// number of acrotelm soil layers
	int nacrotelm;
	/// number of catotelm soil layers
	int ncatotelm;
	/// the average wcont over the growing season, for each of the upper soil layers. Used in drought limited establishment. 
	double awcont_upper;
	/// water content of sublayer of upper soil layer for which evaporation from the bare soil surface is possible
	/// fraction of available water holding capacity 
	double wcont_evap;
	/// water fluxes between soil layers (mm)
	double wfluxes[NSOILLAYERS_MAX + 1];
	/// rainfall and snowmelt today (mm)
	double rain_melt;
	/// rainfall and snowmelt temperature today (C)
	double rain_melt_temp;
	/// upper limit for percolation (mm)
	double max_rain_melt;
	/// water added to wetland today (mm)
	double wmass_wetland_added;
	/// whether to percolate today
	bool percolate;
	/// Soil temperature padding
	Padding pad;
	/// Snow
	Snow snow;
	/// Surface standing water
	SSW ssw;

	// Daily storage:

	/// Depth in mm of the maximum depth of thaw this year
	double maxthaw;
	/// daily thaw depth. The depth to the first soil layer with a temperature greater than 0 degrees C [mm]
	double thaw;
	/// monthly thawing depth full, where ALL the ice has melted [mm]
	double mthaw[12];
	//
	double mald_mean[12];
	//
	Historic<double, 10> maxthaw_10;

	/// acrotelm CO2 level [mimil L-1]
	double acro_co2;

	// Daily photosynthetic limits:

	/// Daily limit on moss photosynthetic activity due to dessication [0,1], but really [0.3,1]
	double dmoss_wtp_limit;
	/// Daily limit on graminoid photosynthetic activity as WTP drops.
	double dgraminoid_wtp_limit;

	// Peatland hydrology variables:

	/// daily water table position [mm]
	double wtp;
	/// monthly average water table position [mm]
	double mwtp[12];
	/// annual average water table position [mm]
	double awtp;
	/// Water table depth from yesterday and updated today = -wtp (range -100, +300) [mm]
	double wtd;

	// Gas diffusion variables:

	/// gas transport velocity of O2 [m d-1]
	double k_O2;
	/// gas transport velocity of CO2 [m d-1]
	double k_CO2;
	/// gas transport velocity of CH4 [m d-1]
	double k_CH4;
	/// Equilibrium concentration of O2 [mol L-1]
	double Ceq_O2;
	/// Equilibrium concentration of CO2 [mol L-1]
	double Ceq_CO2;
	/// Equilibrium concentration of CH4 [mol L-1]
	double Ceq_CH4;

	/// records the first time T is calculated
	bool firstTempCalc;

	// Variables for fast spinup of SOM pools

	/// daily percolation (mm)
	double dperc;
	/// fraction of decayed organic nitrogen leached each day;
	double orgleachfrac;
	/// soil NH4 mass in pool (kgN/m2)
	double NH4_mass;
	/// soil NO3 mass in pool (kgN/m2)
	double NO3_mass;
	/// soil NH4 mass input (kgN/m2)
	double NH4_input;
	/// soil NO3 mass input (kgN/m2)
	double NO3_input;
	double nmin_gross;
	double nimmob;
	double nimmob_pool[NSOMPOOL];
	double nmin_gross_pool[NSOMPOOL];
	double decay_reduction[NSOMPOOL];
	/// annual sum of gross nitrogen mineralisation
	double anmin;
	/// annual sum of nitrogen immobilisation
	double animmob;
	/// annual leaching from available nitrogen pool
	double aminleach;
	/// annual leaching of organics from active nitrogen pool
	double aorgNleach;
	/// total annual nitrogen fixation
	double anfix;
	/// mean annual nitrogen fixation used in equilsom()
	double anfix_mean;
	/// calculated annual mean nitrogen fixation
	double anfix_calc;
	/// annual leaching of organics from carbon pool
	double aorgCleach;	
	/// monthly fraction of available NH4 taken up
	double fNH4uptake_mean[12];
	/// monthly fraction of available NO3 taken up
	double fNO3uptake_mean[12];

	// Daily net mineralisation (kgN/m2/day)
	double net_mineralisation;

	// Daily mineral N leaching (kgN/m2/day)
	double dminleach;

	// Daily leaching of organics from active nitrogen pool (kgN/m2/day)
	double dorgNleach;

	/// Cumulative annual nitrogen flux to atmosphere (kgN/m2/year).
	double anloss_gas;
	// Daily nitrogen flux to atmosphere (kgN/m2/day).
	double nloss_gas;

	// Variables for fast spinup of SOM pools
	// Restart variables

	/// year at which to calculate equilibrium soil carbon
	int solvesom_end;
	/// year at which to begin documenting means for calculation of equilibrium soil carbon
	int solvesom_begin;

	/// years at which to begin documenting for calculation of Century equilibrium
	int solvesomcent_beginyr;
	/// years at which to end documentation and start calculation of Century equilibrium
	int solvesomcent_endyr;
	/// document for calculation of Century equilibrium
	bool savesomcent;
	/// Slope gradient [m m-1]
	double slope_grade_tan;
	double slope_grade_cos;
	/// soil water potential averaged over layers (MPa)
	double soil_psi;

private:

	// Private member variables.
	static const bool ifallowphasechanges = true;
	static const bool snowdensityconstant = false;


	// MEMBER FUNCTIONS

public:

	/// constructor (initialises member variable patch)
	Soil(Patch& p, Soiltype& s) : patch(p), soiltype(s) {
		initdrivers();
	}
	
	// Destructor
	~Soil() {

		// Performs memory deallocation, closing of files or other "cleanup" functions.
	}

	/// Initialize variables
	void initdrivers();
	/// Create soillayers
	void create_soillayers(double soil_layer_depths[NSOILLAYERS_MAX]);
	/// Add a soil layer
	void add_soillayer(double soil_layer_thick, double soil_layer_depth_mid, bool acrotelm, bool catotelm);
	/// Returns the Climate for this Patch
	/** This function returns a const reference to prevent code which operates
	 *  on a patch basis to modify the climate and thereby affect other
	 *  patches/stands.
	 */
	const Climate& get_climate() const;
	/// Set slope properties
	void set_slope_properties();
	/// Initiate soil water content
	void init_water();
	/// Returns surface litter C content
	double get_surflittercmass();
	/// Returns surface litter N content
	double get_surflitternmass();
	/// C content of soil 
	double ccont(grouppooltype group = GALL) const;
	/// N content of soil 
	double ncont(grouppooltype group = GALL) const;
	/// Water and ice content of the soil
	double water_and_ice_content() const;
	/// Get soil property totals
	double get_sum(double Soillayer::* var_member) const;
	/// Get soil property totals of a certain month
	double get_sum(double(Soillayer::* var_arr)[12], int m) const;
	/// Get soil layers available nmass fraction
	void get_nfrac(n_pref_type pref, double nfrac[NSOILLAYERS_MAX], double nmin, bool water_limited = false);
	/// Get soil layers root fractions
	void calc_rootfrac(double rootfrac[NSOILLAYERS_MAX]);
	/// Get soil layers root C mass
	void calc_rootmass(double cmass_root[NSOILLAYERS_MAX]);
	/// Get total soil nmass_avail (if water_limited true then water up to awc limits nmass availability) 
	double nmass_avail(n_pref_type pref = NO, bool water_limited = false);
	/// Nitrogen uptake from soil nmass_avail
	void nmass_avail_uptake(double uptake, n_pref_type pref);
	/// Nitrogen increment to soil nmass_avail
	void nmass_avail_inc(double inc, n_pref_type pref);
	/// Calculate root volume
	void calc_root_volume();
	/// Returns a fine vertical resolution C soil fraction
	std::vector<double> interp_soil_conserve(double fineDz);
	/// Vertical movement of mineral and organic matter
	void vertical_movement();
	/// Calculate this months average soil layer wcont
	void calc_mwcont();
	/// Calculate this months average soil layer temperature
	void calc_mtemp();
	/// Function to get temperature at a specific depth
	double temp_at_depth(double depth) const;
	/// gtemp weighted against root distribution
	double gtemp_rootdist(double rootdist[NSOILLAYERS_MAX]);
	/// Returns the soil temperature for a certain depth interval
	double temp_depth(double depth_start, double depth_end);
	/// Returns the month average soil temperature for a certain depth interval
	double temp_depth(double depth_start, double depth_end, int mth);
	/// Returns the value of some vector within a certain depth interval
	double var_depth(double depth_start, double depth_end, double Soillayer::* member) const;
	double var_depth(double depth_start, double depth_end, double (Soillayer::*getter)() const) const;
	double var_depth(double depth_start, double depth_end, double(Soillayer::* var_arr)[12], int m) const;
	/// Returns the value of some vector within a certain depth interval with a scalar
	double var_depth(double depth_start, double depth_end, double Soillayer::* var_member, double Soillayer::* scal_member) const;
	double var_depth(double depth_start, double depth_end, double(Soillayer::* var_arr)[12], double Soillayer::* scal_member, int mth) const;
	/// Returns the value of soil C within a certain depth interval
	double var_depth_ccont(double depth_start, double depth_end, grouppooltype group = GALL);
	/// Returns the value of soil N within a certain depth interval
	double var_depth_ncont(double depth_start, double depth_end, grouppooltype group = GALL);
	/// Determine snow or rainfall, and melt if simple snow scheme is used
	void snow_or_rainfall(double airT, double& prec, double u10, double& rain_melt, double& rain_melt_temp);
	/// Add snowfall
	void snowfall(double airT, double& prec, double u10);
	/// Calculate surface evaporation, infiltration and runoff
	void compute_evap_infilt_and_runoff(double dtsub, double qflx_rain_melt, double depth_infilt,
		double imass_infilt, double wsats_infilt, double imp_ice, double Kw_max, double fevap, double& qflx_evap, 
		double& qflx_surf_runoff, double& qflx_ssw_runoff, double& qflx_in_surf, double& qflx_infilt, double& qflx_ssw_drain); // DW_COMMENT
	/// Update matric potential
	void update_matric_potential();
	/// Update hydraulic conductivity and diffusivity
	void update_hydraulic_properties(int nsoillayers_sim);
	/// Update soil layer temperatures due to water fluxes
	void wflow_temp(double temp_infilt);
	/// Update heat capacity and conductivity for soil layers
	void update_thermal_properties();
	/// Soil water freezing allowed?
	bool can_freeze() const;
	/// Updates the CO2 concentration in the acrotelm water
	void update_acrotelm_co2(double atmo_co2);
	/// return snow depth [mm] 
	double snowdepth();

	/// Peatland hydrology routine. Implements the peatland hydrology scheme of Wania et al. (2008)
	void hydrology_peat(const Climate& climate, double fevap);
	/// If it is wet season for wetland hydrology
	bool is_wet_season();
	/// Mineral wetland hydrology (simplistic placeholder)
	void hydrology_wetland(const Climate& climate, double fevap);
	/// Two layer hydrology routine from LPJ-GUESS v4.0
	void hydrology_gerten(const Climate& climate, double fevap);
	/// 15 layer hydrology routine from LPJ-GUESS v4.1
	void hydrology_miller(const Climate& climate, double fevap);
	/// Flexible layer hydrology routine
	void hydrology_richards(const Climate& climate, double fevap);
	void compute_moisture_fluxes_and_derivs(double qflx_infl, double Kw[NSOILLAYERS_MAX], double Psi[NSOILLAYERS_MAX], 
		double dKwdw[NSOILLAYERS_MAX], double dPsidw[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX], double qout[NSOILLAYERS_MAX],
		double dqidw0[NSOILLAYERS_MAX], double dqidw1[NSOILLAYERS_MAX], double dqodw1[NSOILLAYERS_MAX], double dqodw2[NSOILLAYERS_MAX],
		double vws_liq[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], int nsoillayers_sim);
	void compute_RHS_moisture_form(double qflx_root[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX], 
		double qout[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX], double rmx[NSOILLAYERS_MAX], 
		int nsoillayers_sim);
	void compute_LHS_moisture_form(double dt_dz[NSOILLAYERS_MAX], double dqidw0[NSOILLAYERS_MAX], double dqidw1[NSOILLAYERS_MAX], 
		double dqodw1[NSOILLAYERS_MAX], double dqodw2[NSOILLAYERS_MAX], double amx[NSOILLAYERS_MAX], double bmx[NSOILLAYERS_MAX], 
		double cmx[NSOILLAYERS_MAX], int nsoillayers_sim);
	void TridiagonalCol(double a[NSOILLAYERS_MAX], double b[NSOILLAYERS_MAX], double c[NSOILLAYERS_MAX], 
		double r[NSOILLAYERS_MAX], double u[NSOILLAYERS_MAX], int nsoillayers_sim);
	void check_error(double dwat[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX], 
		double qout[NSOILLAYERS_MAX], double qflx_root[NSOILLAYERS_MAX], int nsoillayers_sim, double dtmin, double xTolerUpper, double& dtsub,
		double& errorMax, bool& calcsub);
	void update_wtp();
	/// initial percolation
	void initial_percolation();

	// CH4 methods
	/// Calculate and return diffusivities
	void calculate_gas_diffusivities(double CH4[NSOILLAYERS_MAX], double CO2[NSOILLAYERS_MAX], double O2[NSOILLAYERS_MAX]);
	/// Calculate CH4 and CO2 production 
	void calculate_ch4_and_co2_production();
	/// Update gas transport velocities[m d - 1] and equilibrium gas concentrations[mmol m - 3] for CH4, CO2 and O2
	void update_daily_gas_parameters();
	/// Generic gas diffusion method that works with CO2, CH4 and O2
	void diffuse_gas(double D[NSOILLAYERS_MAX], gastype thisgastype, double Ceq,
		double kgas, double &dailyDiff);
	/// Plant transport of O2 or CH4
	void plant_gas_transport(double Ceq, double kgas, gastype thisgastyp,
		double tillers, double& plantCH4TransportToday);
	/// Calculate the area of tillers
	void calculate_tiller_areas(double& tiller_density);
	/// CH4 oxidation
	void oxidation_ch4();
	/// Gas ebullition
	void calculate_gas_ebullition(double& ebull_today);
	/// Calculate methane dynamics today
	void methane(bool generatemethane);

	/// Calculate soil temperature
	void temp();
	/// Calculate multilayer soil temperature
	void calc_temp();
	/// Calculate analytic soil temperatures
	void calc_temp_analytic();
	/// Calculate water and ice for each layer
	void phase_change();
	/// Calculate thaw depth
	void thaw_depth();
	/// Calculate monthly snow depth
	void calc_msnow();
	void serialize(ArchiveStream & arch);

	void daily_accounting();
private:
	void annual_accounting();
};

/// Container for crop-specific data at patchpft level
struct cropphen_struct : public Serializable {

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
	/// whether in period of heat unit sampling
	bool hu_samplingperiod;
	/// number of heat unit sampling days
	int hu_samplingdays;
	/// harvest index today [0-1, >1 if below-ground ho], harvestable organ/above-ground C for above-ground harvestable organs, dependent on fphu, reduced by water stress
	double hi;
	/// fraction of harvest index today
	double fhi;
	/// phenology (fphu) contribution of fraction of harvest index today
	double fhi_phen;	//Phenology (fPHU) compoment of fhi
	/// water stress contribution of fraction of harvest index today
	double fhi_water;
	/// fraction of harvest index at latest harvest
	double fhi_harv;
	/// sum of crop patch demand (patch.wdemand) during crop growing period, reset on harvest day
	double demandsum_crop;
	/// sum of crop supply (patchpft.wsupply) during crop growing period, reset on harvest day
	double supplysum_crop;

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

	double vdsum_alloc;
	double vd;

	/// The fraction of the daily assimilates allocated to roots.
	double f_alloc_root;
	/// The fraction of the daily assimilates allocated to leaves.
	double f_alloc_leaf;
	/// The fraction of the daily assimilates allocated to harvestable organs, seeds.
	double f_alloc_horg;
	/// The fraction of the daily assimilates allocated to stem.
	double f_alloc_stem;
	/// Development stage (0-2) from Wang & Engel 1998
	double dev_stage;
	// A variable holding the memory of whether this field was fertilised or not.
	bool fertilised[3];

	cropphen_struct() {
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
		growingdays=0;
		growingdays_y=0;
		lgp=0;

		for(int j=0;j<2;j++) {
			sdate_harvest[j]=-1;
			hdate_harvest[j]=-1;
			sdate_thisyear[j]=-1;
		}

		vdsum_alloc=0.0;
		vd = 0.0;
		f_alloc_root=0.0;
		f_alloc_leaf=0.0;
		f_alloc_horg=0.0;
		f_alloc_stem=0.0;
		dev_stage = 0.0;

		fertilised[0] = false;
		fertilised[1] = false;
		fertilised[2] = false;
	}

	void serialize(ArchiveStream& arch);
};


/// State variables common to all individuals of a particular PFT in a particular patch
/** Used in individual and cohort modes only. */
class Patchpft : public Serializable {

	// MEMBER VARIABLES:

public:

	/// id code (equal to value of member variable id in corresponding Pft object)
	int id;
	/// reference to corresponding Pft object in PFT list
	Pft& pft;
	/// potential annual net assimilation (leaf-level net photosynthesis) at forest floor (kgC/m2/year)
	double anetps_ff;
	/// water stress parameter (0-1 range; 1=minimum stress)
	double wscal;
	/// five day record of wscal. Used in daily grass growth
	Historic<double, 5> wscal_5;
	/// ten day record of wscal. Used in DAVE daily grass growth for determination of dormancy.
	Historic<double, 10> wscal_10;
	/// ten day record of temperature of upper soil layer. Used in DAVE daily grass growth for determination of dormancy.
	Historic<double, 10> soiltemp_10;
	/// thirty day mean of wscal. Used in DAVE daily grass growth for determination of dormancy.
	Historic<double, 30> wscal_30;
	/// Populate all wscal fields with today's wscal value.
	void set_wscal(double wstress);
	/// Flag to track whether wscal has been set today, to avoid multiple
	/// updates per day.
	bool wscal_set_today;

	// Assaf: ???? Probably reason for senescence
	sen_type senescence_reason;
	/// ten day record of air temperature. Used in DAVE daily grass growth for determination of dormancy.
	Historic<double, 10> airtemp_10;
	/// ten day record of daylength. Used in DAVE daily grass growth for determination of growth stages.
	Historic<double, 10> daylength_10;
	// Daylength relative to maximum daylength in a year.
	double reldaylength;
	/// cumulative temperature factor for growth.
	double ftemp_cum;
	// Water content of upper soil layer as fraction of available water holding capacity over the last 4 days
	Historic<double, 4> wcont_germ_4;
	// Water content of upper soil layer as fraction of available water holding capacity over the last 30 days
	Historic<double, 30> wcont_germ_30;
	// Water content of upper soil layer as fraction of available water holding capacity over the last 3 days
	Historic<double, 14> wcont_germ_14;
	// Temperature of upper soil layer over the last 3 days
	Historic<double, 4> tsoil_germ_4;
	// Temperature of upper soil layer over the last 14 days
	Historic<double, 14> tsoil_germ_14;
	/// running sum (converted to annual mean) for wscal
	double wscal_mean;
	/// potential annual net assimilation at forest floor averaged over establishment interval (kgC/m2/year)
	double anetps_ff_est;
	/// first-year value of anetps_ff_est
	double anetps_ff_est_initial;
	/// annual mean wscal averaged over establishment interval
	double wscal_mean_est;
	/// vegetation phenological state (fraction of potential leaf cover), updated daily
	double phen;
	/// monthly sum of daily fractional leaf cover
	double mphen[12];
	/// annual sum of daily fractional leaf cover
	/** equivalent number of days with full leaf cover
	 *  (reset on expected coldest day of year)
	 */
	double aphen;
	/// Month of lowest phen for RAINGREEN
	int driest_mth;
	/// Monthly average leaf C mass for this patch PFT. Used for output
	double mcmass_leaf[12];
	/// Monthly average leaf N mass for this patch PFT. Used for output
	double mnmass_leaf[12];
	/// whether PFT can establish in this patch under current conditions
	bool establish;
	/// running total for number of saplings of this PFT to establish (cohort mode)
	double nsapling;
	/// The number of saplings of this PFT established this year.
	double nsapling_year;
	/// leaf-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_leaf;
	/// fine root-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_root;
	/// remaining sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_sap;
	/// remaining heartwood-derived litter for PFT on modelled area basis (kgC/m2)
	double cmass_litter_heart;
	/// C litter derived from allocation to reproduction for PFT on modelled area basis (kgC/m2)
	double cmass_litter_repr;
	/// C in hypothetical external seed bank (kgC/m2)
	double cmass_ext_seed;
	/// N in hypothetical external seed bank (kgN/m2)
	double nmass_ext_seed;
	/// carbon lost during mortality
	double cmass_mort;
	/// carbon removed from live vegetation by Individual::kill() today
	double dcmass_mort_kill;
	/// annual carbon removed from live vegetation by Individual::kill()
	double acmass_mort_kill;
	/// carbon lost during fire
	double cmass_fire;
	/// carbon lost during ´disturbance
	double cmass_dist;
	/// carbon lost during turnover
	double cmass_leaf_root_turnover;
	/// carbon lost to reproduction
	double cmass_repr;
	/// carbon added at establishment
	double cmass_est;
	/// carbon in biomass killed in harvest this year (kgC/m2)
	double cmass_harv_killed;
	/// harvested wood C before removing part to the product pool (kgC/m2)
	double cmass_wood_harv;
	/// harvested wood C removed to the product pool (kgC/m2)
	double cmass_wood_harv_toprod;
	/// trunk part of carbon harvested in clearcut
	double cmass_wood_clearcut;
	/// harvested tree C left as litter (kgC/m2)
	double cmass_harv_tolitter;

	/// leaf-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_leaf;
	/// root-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_root;
	/// remaining sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_sap;
	/// remaining heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_heart;
	/// N content of reproductive litter pool (seedbank) (kgN/m2)
	double nmass_litter_repr;
	/// remaining NH4-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_NH4;
	/// Half saturation concentration for N uptake (Rothstein 2000, Macduff 2002) Michaelis-Menten kinetic parameters
	double Km;

	/// non-FPC-weighted canopy conductance value for PFT under water-stress conditions (mm/s)
	double gcbase;
	/// daily value of the above variable (mm/s)
	double gcbase_day;

	/// Vegetation root profile in an array containing fraction of roots in each soil layer, [0=upper layer]
	double rootdist[NSOILLAYERS_MAX];

	/// evapotranspirational "supply" function for this PFT today (mm/day)
	double wsupply;
	double wsupply_leafon;
	/// fractional uptake of water from each soil layer today
	double fwuptake[NSOILLAYERS_MAX];

	/// whether water-stress conditions for this PFT
	bool wstress;
	/// daily version of the above variable
	bool wstress_day;

	/// carbon depository for long-lived products like wood
	double cmass_harvested_products_slow;
	/// nitrogen depository for long-lived products like wood
	double nmass_harvested_products_slow;
	/// first and last day of crop sowing window, calculated in crop_sowing_patch() or Crop_sowing_date_new()
	int swindow[2];
	/// daily value of water deficit, calculated in irrigated_water_uptake()
	double wmass_deficit_d;
	/// yearly sum of water deficit
	double wmass_deficit_y;

	/// Struct for crop-specific variables
	cropphen_struct *cropphen;

	// INUNDATION STRESS TERMS:
    
	/// Number of days a month that the water table is above this PFT's wtp_max (updated on the first day of the month)
	int inund_count;
	/// [0,1] - a measure of the inundation stress. Daily photosynthesis is reduced by this factor.
	double inund_stress;

	// Cumulative dormancy depth (1-fphen_activity) this year.
	double dormancy_depth;
	// Cumulative wscal this year.
	double wscal_cum;
	// Number of days this Pft has been established/alive.
	int ndays_alive;
	// Moving average over N years of the annual dormancy depth of this cohort.
	HistoricDynamic<double> dormancy_depth_historical;
    // Historical average over N years of the wscal of this cohort.
	HistoricDynamic<double> wscal_historical;

	// Maximum total LAI recorded on any day this year (m2/m2).
	double lai_ymax;

	// Maximum above-ground biomass recorded on any day this year (kg/m2).
	double agb_ymax;

	// Mean LAI over the year (m2/m2). Note: this is a rolling sum, divided by
	// year length only on last day of year.
	double lai_ymean;

	// Mean above-ground biomass over the year (kg/m2). Note: this is a rolling
	// sum, divided by year length only on last day of year.
	double agb_ymean;

	// Total fraction of cohort that was killed this year by carbon starvation
	// mortality, summed over all days of the year and cohorts of this PFT in
	// this patch.
	double amort_cstarve;

	// Total number of trees that were killed this year by carbon starvation
	// mortality, summed over all days of the year and cohorts of this PFT in
	// this patch.
	double amort_cstarve_ntree;

	// Daily number of trees that were killed by carbon starvation mortality,
	// summed over all cohorts of this PFT in this patch.
	double dmort_cstarve_ntree;

	// MEMBER FUNCTIONS:

	/// Constructor: initialises id, pft and data members
	Patchpft(int i,Pft& p):id(i),pft(p) {

		cmass_litter_leaf = 0.0;
		cmass_litter_root = 0.0;
		cmass_litter_sap  = 0.0;
		cmass_litter_heart = 0.0;
		cmass_litter_repr = 0.0;
		cmass_ext_seed = 0.0;
		dcmass_mort_kill = 0.0;
		acmass_mort_kill = 0.0;

		nmass_litter_leaf = 0.0;
		nmass_litter_root = 0.0;
		nmass_litter_sap = 0.0;
		nmass_litter_heart = 0.0;
		nmass_litter_repr  = 0.0;
		nmass_ext_seed     = 0.0;
		nmass_litter_NH4 = 0.0;

		wscal = 1.0;
		wscal_mean = 1.0;
		anetps_ff = 0.0;
		aphen = 0.0;
		phen = 0.0;
		driest_mth = 0;
		wsupply = 0.0;
		wsupply_leafon = 0.0;
		anetps_ff_est = 0.0;
		anetps_ff_est_initial = 0.0;
		wscal_mean_est = 0.0;
		nsapling = 0;
		nsapling_year = 0;

		for (int mth = 0; mth < 12; mth++) {
			mphen[mth] = 0.0;
			mcmass_leaf[mth] = 0.0;
			mnmass_leaf[mth] = 0.0;
		}

		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			rootdist[lyr] = 0.0;
			fwuptake[lyr] = 0.0;
		}
		rootdist[0] = 1.0;

		cropphen = NULL;
		cmass_harvested_products_slow = 0.0;
		nmass_harvested_products_slow = 0.0;

		swindow[0]=-1;
		swindow[1]=-1;

		if(pft.landcover==CROPLAND)
		{
			cropphen=new cropphen_struct;
		}

		wmass_deficit_d = 0.0;

		inund_count=0;
		inund_stress=1.0; // No stress by default

		cmass_wood_harv = 0.0;
		cmass_wood_harv_toprod = 0.0;
		cmass_wood_clearcut = 0.0;
		cmass_harv_tolitter = 0.0;
		cmass_harv_killed = 0.0;

		dormancy_depth = 0;

		// 30 years by default.
		dormancy_depth_historical.set_capacity(30);
		wscal_historical.set_capacity(30);

		// Initial values for long-term mean dormancy_depth and wscal
		dormancy_depth_historical.add(0.5);
		wscal_historical.add(0.7);

		senescence_reason.reset();
		lai_ymax = 0;
		agb_ymax = 0;
		lai_ymean = 0;
		agb_ymean = 0;
		wscal_set_today = false;
		amort_cstarve = 0.0;
		amort_cstarve_ntree = 0.0;
		dmort_cstarve_ntree = 0.0;
	}

	~Patchpft() {
		if(cropphen) {
			delete cropphen;
		}
	}

	/// safe method to obtain cropphen_struct pointer
	cropphen_struct* get_cropphen();
	/// safe method to set cropphen_struct variables
	cropphen_struct* set_cropphen();

	bool growingseason() const;

	void update_rootdist(double maxthaw, double soildepth, double soil_lyr_depth[NSOILLAYERS_MAX], int nsoillayers);

	void serialize(ArchiveStream& arch);

    double total_litter() const;

	/// Reset daily variables.
	void daily_accounting(const Patch &patch);
private:
	/// Reset annual variables. Called from daily_accounting() on the first day
	/// of the year.
	void annual_accounting(const Patch &patch);
};


/// Stores data for a patch.
/** In cohort and individual modes, replicate patches are
 *  required in each stand to accomodate stochastic variation; in population mode there
 *  should be just one Patch object, representing average conditions for the entire
 *  stand. A reference to the parent Stand object (defined below) is included as a
 *  member variable.
 */
class Patch : public Serializable {

public:

	// MEMBER VARIABLES

	/// id code in range 0-npatch for patch
	int id;
	/// reference to parent Stand object
	Stand& stand;
	/// list array [0...npft-1] of Patchpft objects (initialised in constructor)
	ListArray_idin1<Patchpft,Pft> pft;
	/// vegetation for this patch
	Vegetation vegetation;
	/// soil  for this patch
	Soil soil;
	/// fluxes for this patch
	Fluxes fluxes;
	/// FPAR at top of grass canopy today
	double fpar_grass;
	/// FPAR at soil surface today
	double fpar_ff;
	/// mean growing season PAR at top of grass canopy (J/m2/day)
	double par_grass_mean;
	/// number of days in growing season, estimated from mean vegetation leaf-on fraction
	/** \see function fpar in canopy exchange module */
	int nday_growingseason;
	/// total patch FPC
	double fpc_total;
	/// total patch RPC (root projective cover)
	double rpc_total;
	/// whether patch was disturbed last year
	bool disturbed;
	/// patch age (years since last disturbance)
	int age;
	/// probability of fire this year (GlobFIRM)
	double fireprob;

	/// BLAZE if patch has burnt this year
	bool burned;
	/// BLAZE Fire line intensity (kW/m)
	double fire_line_intensity;
	/// BLAZE available fuel (g/m2)
	double avail_fuel;

	// BLAZE fire related carbon fluxes
	/// BLAZE-fire carbon flux: live wood to atmosphere (kgC/m2)
	double wood_to_atm;
	/// BLAZE-fire carbon flux: leaves to atmosphere (kgC/m2)
	double leaf_to_atm;
	/// BLAZE-fire carbon flux: leaves to litter (kgC/m2)
	double leaf_to_lit;
	/// BLAZE-fire carbon flux: live wood to structural litter (kgC/m2)
	double wood_to_str;
	/// BLAZE-fire carbon flux: live wood to fine woody debris (kgC/m2)
	double wood_to_fwd;
	/// BLAZE-fire carbon flux: live wood to coarse woody debris (kgC/m2)
	double wood_to_cwd;
	/// BLAZE-fire carbon flux: fine litter to atmosphere (kgC/m2)
	double litf_to_atm;
	/// BLAZE-fire carbon flux: fine woody debris to atmosphere (kgC/m2)
	double lfwd_to_atm;
	/// BLAZE-fire carbon flux: coarse woody debris to atmosphere (kgC/m2)
	double lcwd_to_atm;

	// Storage for averaging of different Fapars for biome mapping in SIMFIRE
	/// SIMFIRE fapar: Grasses
	double fapar_grass_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Needle-leaf tree
	double fapar_ndlt_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Broad-leaf tree
	double fapar_brlt_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: TrBR (Savanna)
	double fapar_trbr_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Shrubs
	double fapar_shrub_avg[N_YEAR_BIOMEAVG];
	/// SIMFIRE fapar: Total fapar
	double fapar_total_avg[N_YEAR_BIOMEAVG];

	/// whether management has started on this patch (usually, no disturbance or fire)
	bool managed;
	/// whether cutting started on this patch
	bool has_been_cut;
	/// cutting intensity (initial fraction of trees cut, further selection at individual level is done in a separate function)
	double man_strength;

	/// Non-commercial thinning (harvest to litter)
	bool harvest_to_litter;
	/// Whether patch has been managed this year
	bool managed_this_year;
	/// Whether patch has been clearcut this year
	bool clearcut_this_year;
	/// Latest cutting interval (patch age at year of clearcut)
	int cutinterval_actual;
	/// Cutting interval if patch has been cut this year (patch age at year of clearcut)
	int cutinterval_actual_thisyear;
	/// Whether patch should be planted this year
	bool plant_this_year;
	/// Whether man_strength has been partitioned on individuals this year
	bool distributed_cutting;
	/// Whether tree density is below limit after firstclearcutyear (when using tree density as a trigger for clearcut)
	bool cut_due;
	/// Initial tree density (when using Reineke's rule-based automatic thinning)
	double dens_start;

	/// DLE - the number of days over which wcont is averaged for this patch
	/** i.e. those days for which daily temp > 5.0 degC */
	int growingseasondays;


	// Variables used by new hydrology (Dieter Gerten 2002-07)

	/// interception by vegetation today on patch basis (mm)
	double intercep;
	/// annual sum of AET (mm/year) for each of the last five simulation years
	Historic<double, NYEARAAET> aaet_5;
	/// daily transpiration (mm/day)
	double transp;
	/// water added to wetlands today (mm)
	double wetland_water_added_today;
	/// annual sum of water added to wetlands (mm/year)
	double awetland_water_added;
	/// annual sum of potential evapotranspiration (mm/year)
	double apet;

	/// equilibrium evapotranspiration today, deducting interception (mm)
	double eet_net_veg;

	/// transpirative demand for patch, patch vegetative area basis (mm/day)
	double wdemand;
	/// daily average of the above variable (mm/day)
	double wdemand_day;
	/// transpirative demand for patch assuming full leaf cover today
	/** mm/day, patch vegetative area basis	*/
	double wdemand_leafon;
	/// rescaling factor to account for spatial overlap between individuals/cohorts populations
	double fpc_rescale;

	/// daily nitrogen demand
	double ndemand;

	/// annual nitrogen fertilization (kgN/m2/year)
	double anfert;
	/// daily nitrogen fertilization (kgN/m2/day)
	double dnfert;

	/// whether litter is to be sent to the soil today
	bool is_litter_day;
	/// number of harvests and/or cover-crop killing or turnover events
	int nharv;
	/// whether today is a harvest day and/or cover-crop killing or turnover day
	bool isharvestday;

	// MEMBER FUNCTIONS

	/// Constructor: initialises various members and builds list array of Patchpft objects.
	Patch(int i,Stand& s,Soiltype& st);

	void serialize(ArchiveStream& arch);

	/// Returns the Climate for this Patch
	/** This function returns a const reference to prevent code which operates
	 *  on a patch basis to modify the climate and thereby affect other
	 *  patches/stands.
	 */
	const Climate& get_climate() const;

	/// Returns whether we should model fire in this patch
	bool has_fires() const;

	/// Returns whether we should model disturbances in this patch
	bool has_disturbances() const;
	/// Total patch carbon biomass and litter
	double ccont(double scale_indiv = 1.0, bool luc = false);
	/// Total patch nitrogen biomass and litter
	double ncont(double scale_indiv = 1.0, bool luc = false) const;
	/// Total patch water and ice mass
	double water_content() const;
	/// Total patch carbon fluxes so far this year
	double cflux(Period period) const;
	/// Total patch nitrogen fluxes so far this year
	double nflux(Period period) const;
	/// Total patch water fluxes so far this year
	double water_flux(Period period) const;
	
	/// Get 5-year mean of tree wood C mass increase (periodic annual increment)
	double get_tree_cmass_wood_inc_5() {
		double cmass_wood_inc_5_mean = 0.0;
		for (unsigned int i=0; i<vegetation.nobj; i++) {

			// Disregard shrubs (crownarea_max = 10)
			Individual& indiv = vegetation[i];
			if(indiv.pft.lifeform == TREE && !indiv.pft.is_shrub()) {
				if(indiv.cmass_wood_inc_5.size())
					cmass_wood_inc_5_mean += indiv.cmass_wood_inc_5.mean();
			}
		}
		return cmass_wood_inc_5_mean;
	}
	
	/// Get cmass_wood of all individuals in patch (including shrubs unless function parameter exclude_shrubs is true)
	double cmass_wood(bool exclude_shrubs = false) {
		double cmass_wood = 0.0;
		for (unsigned int i=0; i<vegetation.nobj; i++) {

			Individual& indiv = vegetation[i];
			if(!exclude_shrubs || !indiv.pft.is_shrub()) {
				cmass_wood += indiv.cmass_wood();
			}
		}
		return cmass_wood;
	}
};

/// Container for variables common to individuals of a particular PFT in a stand.
/** Used in individual and cohort modes only
 */
class Standpft : public Serializable {

public:

	// MEMBER VARIABLES

	int id;
	Pft& pft;
	/// net C allocated to reproduction for this PFT in all patches of this stand this year (kgC/m2)
	double cmass_repr;
	/// maximum value of Patchpft::anetps_ff for this PFT in this stand so far in the simulation (kgC/m2/year)
	double anetps_ff_max;

	/// FPC sum for this PFT as average for stand
	double fpc_total;

	/// Photosynthesis values for this PFT under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// Whether this PFT is allowed to grow in this stand
	bool active;
	/// Order of PFT in selection string
	int selection;
	/// Whether this PFT is planted in this stand
	bool plant;
	/// Whether this PFT is allowed to establish (after planting) in this stand
	bool reestab;
	/// Planting density (seedlings/ha)
	double plantdensity;
	/// target cmass fraciton
	double targetfrac;

	/// Whether this PFT is irrigated in this stand
	bool irrigated;
	/// sowing date specified in stand type or read from input file
	int sdate_force;
	/// harvest date specified in stand type or read from input file
	int hdate_force;

	// MEMBER FUNCTIONS

	/// Constructor: initialises various data members
	Standpft(int i,Pft& p):id(i),pft(p) {

		anetps_ff_max = 0.0;
		active = !run_landcover;
		selection = -1;
		plant = false;
		reestab = false;
		plantdensity = -1;
		targetfrac = 0;
		irrigated = false;
		sdate_force = -1;
		hdate_force = -1;
		cmass_repr = 0.0;
	}

	void serialize(ArchiveStream& arch);

	/// Reset daily variables.
	void daily_accounting();
private:
	/// Reset annual variables. Called from daily_accounting() on the first day
	/// of the year.
	void annual_accounting();
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

	/// stand type id
	int stid;

	/// pft id of main crop, updated during rotation
	int pftid;
	/// Number of PFTs in selection
	int npft_selection;
	/// current management rotation item
	int current_rot;
	/// number of years passed in current rotation item
	int nyears_in_rotation;
	/// number of days passed in current rotation item
	int ndays_in_rotation;
	/// Returns true if stand is in fallow (with cover crop grass)
	bool infallow;
	/// Returns true if crop rotation item is to be updated today
	bool isrotationday;
	/// Returns true if current crop management hydrology == irrigated, updated during rotation
	bool isirrigated;
	/// Returns true if the stand's main crop pft intercrop==naturalgrass and a pft with isintercrop==true is in the pftlist.
	bool hasgrassintercrop;
	/// gdd5-value at first intercrop grass growth
	double gdd5_intercrop;

	/// old fraction of this stand relative to the gridcell before update
	double frac_old;

	/// used during land cover change involving several calls to reveiving_stand_change()
	/** Set to frac_old in reduce_stands(), then modified in donor_stand_change() and receiving_stand_change().
	 */
	double frac_temp;
	/// fraction unavailable for transfer to other stand types
	double protected_frac;
	/// net stand fraction change
	double frac_change;
	/// gross fraction increase
	double gross_frac_increase;
	/// gross fraction decrease
	double gross_frac_decrease;
	/// fraction that has been cloned from another stand
	double cloned_fraction;
	/// Returns true if this stand is cloned from another stand
	bool cloned;
	/// pointer to array of fractions transferred from this stand to other stand types
	double *transfer_area_st;
	/// land cover origin of this stand
	/** Set to landcover for primary stands, for secondary stands: origin lc if created by transfer_to_new_stand..(),
	 *  otherwise to NLANDCOVERTYPES
	 */
	landcovertype lc_origin;
	/// stand type origin of this stand
	int st_origin; 
	//  Variables used for output from separate stands
	/// NPP
	double anpp;
	/// LAI
	double lai;
	/// C mass
	double cmass;
	/// Wood C
	double cmass_wood;
	/// Harvested wood C
	double cmass_wood_harv;
	/// C lost in mortality
	double cmass_mort;

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

	/// type of landcover
	/** \see landcovertype
	 *  initialised in constructor
	 */
	landcovertype landcover;

	/// The simulation year when this stand was created.
	/** Will typically be year zero unless running with dynamic
	 *  land cover.
	 *
	 *  Needed to set patchpft.anetps_ff_est_initial
	 */
	int first_year;
	// The simulation year this stand was cloned from another stand
	int clone_year;
	/// scaling factor for stands that have grown in area this year (old fraction/new fraction)
	double scale_LC_change;
	/// animal units (1 cow = 1.0) per ha
	double animal_units;

	/// The patches which are currently being grazed.
	std::vector<int> patches_grazed;

	/// Number of days that the currently-grazed patches have been grazed.
	int days_grazed;

	// MEMBER FUNCTIONS

	/// Constructs a Stand
	/** \param i         The id for the stand within the grid cell
	 *  \param gc        The parent grid cell
	 *  \param st        The soil type to be used within this Stand
	 *  \param landcover The type of landcover to use for this stand
	 */
	Stand(int i, Gridcell* gc, Soiltype& st, landcovertype landcover, int no_patch = 0);

	~Stand();

	/// Gives the fraction of this Stand relative to the whole grid cell
	double get_gridcell_fraction() const;

	/// Gives the fraction of this Stand relative to its land cover type; NB: unsafe to use within landcover_dynamics() !
	double get_landcover_fraction() const;

	double get_distinterval() const;

	/// Set the fraction of this Stand relative to the gridcell
	void set_gridcell_fraction(double fraction);

	/// Returns the number of patches in this Stand
	unsigned int npatch() const { return nobj; }

	/// Returns true if stand is true crop stand, as opposed to pasture grass grown on cropland or other land cover
	inline bool is_true_crop_stand() {
		return landcover==CROPLAND && pft[pftid].pft.phenology==CROPGREEN;	// OK also for fallow (pftid always cropgreen)
	}
	/// Moves crop rotation forward
	void rotate(int rot = -1);
	/// Sets management parameters
	void set_management();
	/// Sets parameters for the pft:s in selection
	void set_selection_params();
	/// Returns area transferred to other land cover during land cover change
	double transfer_area_lc(landcovertype to);
	/// Initiates new stand land cover settings
	void init_stand_lu(StandType& st, double fraction, bool suppress_disturbance = false);
	/// Total stand carbon biomass and litter
	double ccont(double scale_indiv = 1.0);
	/// Total stand nitrogen biomass and litter
	double ncont(double scale_indiv = 1.0);
	/// Total stand water and ice
	double water_content() const;
	/// Total stand carbon fluxes so far this year
	double cflux(Period period) const;
	/// Total stand nitrogen fluxes so far this year
	double nflux(Period period) const;
	/// Total stand water fluxes so far this year
	double water_flux(Period period) const;
	/// Returns true if stand is true high-latitude peatland stand, as opposed to a wetland < PEATLAND_WETLAND_LATITUDE_LIMIT N
	bool is_highlatitude_peatland_stand() const;
	/// Returns true if stand is wetland stand, as opposed to a peatland >= PEATLAND_WETLAND_LATITUDE_LIMIT N
	bool is_true_wetland_stand() const;

	inline ManagementType& get_current_management() {
		return stlist[stid].get_management(current_rot);
	}

    /// Creates a duplicate stand with a new landcovertype
    /** The new stand is added to this stand's gridcell.
     *
     *  \returns reference to the new stand
     */
    Stand& clone(StandType& st, double fraction, bool suppress_disturbance = false);

	void serialize(ArchiveStream& arch);

	/// Returns the Climate for this Stand
	/** This function returns a const reference to prevent code which operates
	 *  on a stand basis to modify the climate and thereby affect other
	 *  stands.
	 */
	const Climate& get_climate() const;

	/// Returns the Gridcell containing this Stand
	Gridcell& get_gridcell() const;

private:

	/// Pointer to parent object, could be a null pointer
	/** Prefer to access the gridcell through get_gridcell(), even internally
	 *  within the Stand class.
	 */
	Gridcell* gridcell;

	/// Soil type to be used in this Stand
	Soiltype& soiltype;

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

	///Crop-specific variables:
	/// whether the daily temperature has fallen below the autumn temperature limit (tempautumn) this year
	bool autumnoccurred;
	/// whether the daily temperature has risen above the spring temperature limit (tempspring) this year
	bool springoccurred;
	/// whether the daily temperature has fallen below the vernalization limit (trg) this year
	bool vernstartoccurred;
	/// whether the daily temperature rises over the vernalization limit (trg) this year
	bool vernendoccurred;
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
	/// N fertilization (total, including manure) from input file (kgN/m2/year)
	double Nfert_read;
	/// Manure N fertilization from input file (kgN/m2/year)
	double Nfert_man_read;
	/// default harvest date (pft.hlimitdatenh/hlimitdatesh)
	int hlimitdate_default;
	/// whether autumn sowing is either calculated or prescribed
	bool wintertype;
	/// first and last day of crop sowing window, calculated in calc_sowing_windows()
	int swindow[2];
	/// first and last day of crop sowing window for irrigated crops, calculated in calc_sowing_windows()
	int swindow_irr[2];
	/// flag to preclude crop sowing during fallow
	bool sowing_restriction;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellpft object
	/** \param i   The id for this object
	 *  \param p   A reference to the Pft for this Gridcellpft
	 */
	Gridcellpft(int i,Pft& p):id(i),pft(p) {
		addtw = 0.0;

		autumnoccurred=false;
		springoccurred=false;
		vernstartoccurred=false;
		vernendoccurred=false;
		first_autumndate=-1;
		first_autumndate20=-1;
		last_springdate=-1;
		last_springdate20=-1;
		last_verndate=-1;
		last_verndate20=-1;
		for (int year=0;year<20;year++) {
			first_autumndate_20[year]=-1;
			last_springdate_20[year]=-1;
			last_verndate_20[year]=-1;
		}
		sdate_default=-1;
		sdate_force=-1;
		hdate_force=-1;
		Nfert_read=-1;
		Nfert_man_read=-1;
		sdatecalc_temp=-1;
		sdatecalc_prec=-1;
		hlimitdate_default=-1;
		wintertype=false;
		swindow[0]=-1;
		swindow[1]=-1;
		swindow_irr[0] = -1;
		swindow_irr[1] = -1;
		sowing_restriction = false;
	}

	void serialize(ArchiveStream& arch);
};

/// State variables common to all individuals of a particular STANDTYPE in a GRIDCELL.
class Gridcellst : public Serializable {

public:

	// MEMBER VARIABLES

	/// A number identifying this object within its list array
	int id;

	/// A reference to the StandType object for this Gridcellst
	StandType& st;

	/// fraction of this stand type relative to the gridcell
	double frac;
	/// old fraction of this stand type relative to the gridcell before update
	double frac_old;
	/// original input value of old fraction of this stand type before rescaling
	double frac_old_orig;
	/// fraction unavailable for transfer to other stand types
	double protected_frac;

	/// net fraction change
	double frac_change;
	/// gross fraction increase
	double gross_frac_increase;
	/// gross fraction decrease
	double gross_frac_decrease;
	/// disturbance interval
	double distinterval_st;
	/// switch to reset cutinterval_st this year
	bool reset_cutinterval_st;
	/// cutting interval
	double cutinterval_st;
	/// stand type elevation
	double elevation_st;
	/// first management year
	int firstmanageyear_st;

	// current number of stands of this stand type
	int nstands;

	/// Nitrogen fertilisation amount (kgN/m2/year)
	double nfert;

	/// Harvested forest area fraction
	double woodharv_frac;
	/// Wood harvest C mass / m2 for gridcell
	double woodharv_cmass;
	/// Lower tree diameter limit (cm) for cutting (dynamic variable)
	double diam_cut_low;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellst object
	/** \param i   The id for this object
	 *  \param s   A reference to the StandType for this Gridcellst
	 */
	Gridcellst(int i,StandType& s):id(i),st(s) {
		frac = 1.0;
		frac_old = 0.0;
		frac_old_orig = 0.0;
		protected_frac = 0.0;
		frac_change = 0.0;
		gross_frac_increase = 0.0;
		gross_frac_decrease = 0.0;
		nstands = 0;
		nfert = -1.0;
		woodharv_frac = -1.0;
		woodharv_cmass = -1.0;
		diam_cut_low = 0.0;
		reset_cutinterval_st = false;
		cutinterval_st = 0.0;
		firstmanageyear_st = FAR_FUTURE_YEAR;
		distinterval_st = 1.0e10;
		elevation_st = 0;
	}

	void serialize(ArchiveStream& arch);
};

struct forest_lc_frac_transfer {

	double primary[NLANDCOVERTYPES][NLANDCOVERTYPES];
	double secondary_young[NLANDCOVERTYPES][NLANDCOVERTYPES];

	forest_lc_frac_transfer() {
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			for(int j=0;j<NLANDCOVERTYPES;j++) {
				primary[i][j] = 0.0;
				secondary_young[i][j] = 0.0;
			}
		}
	};
};

struct forest_st_frac_transfer {

	double* primary;
	double* secondary_young;

	forest_st_frac_transfer(int nst) {
		primary = NULL;
		secondary_young = NULL;
		primary = new double[nst * nst];
		secondary_young = new double[nst * nst];

		for(int i=0;i<nst*nst;i++) {
			primary[i] = 0.0;
			secondary_young[i] = 0.0;
		}

	};

	~forest_st_frac_transfer(){
		if(primary)
			delete[] primary;
		if(secondary_young)
			delete[] secondary_young;
	}
};

/// Storage of land cover fraction data and some land cover change-related pools and fluxes
struct Landcover : public Serializable {

	Landcover();

	/// The fractions of the different land cover types.
	/** landcoverfrac is read in from land cover input file or from
	 *  instruction file in getlandcover().
	 */
	double frac[NLANDCOVERTYPES];

	/// The land cover fractions from the previous year
	/** Used to keep track of the changes when running with dynamic
	 *  land cover.
	 */
	double frac_old[NLANDCOVERTYPES];

	double frac_change[NLANDCOVERTYPES];

	/// Wood harvest information read from input file.
	Wood_harvest_struct wood_harvest;

	/// Transfer matrices
	double frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES];
	forest_lc_frac_transfer  forest_lc_subset_transfer;

	// Peat duration
	int wet_startday;
	int wet_days;

	/// Whether the land cover fractions changed for this grid cell this year
	/** \see landcover_dynamics
	 */
	bool updated;

	/// Gridcell-level C flux from slow harvested products
	double acflux_harvest_slow;

	/// Gridcell-level C flux from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture
	double acflux_landuse_change;
	/// Gridcell-level C removal from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture
	double acflux_landuse_change_orig;
	/// Gridcell-level C flux from wood harvest at LUC among NATURAL and FOREST stands
	double acflux_wood_harvest;
	/// Gridcell-level C removal at wood harvest at LUC among NATURAL and FOREST stands
	double acflux_wood_harvest_orig;
	/// Gridcell-level C flux from clearing NATURAL and FOREST stands for cropland and pasture
	double acflux_clearing;
	/// Gridcell-level C removal at clearing NATURAL and FOREST stands for cropland and pasture
	double acflux_clearing_orig;
	/// Gridcell-level C removal of harvested tree stem
	double cmass_stem_harvest;
	/// Gridcell-level C harvested tree stem going to the product pool
	double cmass_stem_toprod;
	/// Gridcell-level C of killed trees in wood harvest
	double cmass_harv_killed;
	/// Gridcell-level C of killed trees in wood harvest left on site
	double cmass_harv_tolitter;

	/// Gridcell-level N flux from slow harvested products
	double anflux_harvest_slow;

	/// Gridcell-level N flux from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture
	double anflux_landuse_change;
	/// Gridcell-level N removal from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture
	double anflux_landuse_change_orig;
	/// Gridcell-level N flux from wood harvest at LUC among NATURAL and FOREST stands
	double anflux_wood_harvest;
	/// Gridcell-level N removal at wood harvest at LUC among NATURAL and FOREST stands
	double anflux_wood_harvest_orig;
	/// Gridcell-level N flux from clearing NATURAL and FOREST stands for cropland and pasture
	double anflux_clearing;
	/// Gridcell-level N removal at clearing NATURAL and FOREST stands for cropland and pasture
	double anflux_clearing_orig;

	/// Landcover-level C flux from slow harvested products (donating landcover)
	double acflux_harvest_slow_lc[NLANDCOVERTYPES];

	/// Landcover-level C flux from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture (donating landcover)
	double acflux_landuse_change_lc[NLANDCOVERTYPES];
	/// Landcover-level C flux from wood harvest at LUC among NATURAL and FOREST stands (donating landcover)
	double acflux_wood_harvest_lc[NLANDCOVERTYPES];
	/// Landcover-level C flux from clearing NATURAL and FOREST stands for cropland and pasture (donating landcover)
	double acflux_clearing_lc[NLANDCOVERTYPES];
	/// Landcover-level C flux transferred during cloning of stands
	double acflux_cloned_lc[NLANDCOVERTYPES];

	/// Landcover-level N flux from slow harvested products (donating landcover)
	double anflux_harvest_slow_lc[NLANDCOVERTYPES];

	/// Landcover-level N flux from harvest associated with landcover change other than wood harvest and clearing for cropland and pasture (donating landcover)
	double anflux_landuse_change_lc[NLANDCOVERTYPES];
	/// Landcover-level N flux from wood harvest at LUC among NATURAL and FOREST stands (donating landcover)
	double anflux_wood_harvest_lc[NLANDCOVERTYPES];
	/// Landcover-level N flux from clearing NATURAL and FOREST stands for cropland and pasture (donating landcover)
	double anflux_clearing_lc[NLANDCOVERTYPES];

	/// Which landcover types create new stands when area increases.
	bool expand_to_new_stand[NLANDCOVERTYPES];

	/// Whether to pool all transferred land from a donor landcover (overrides different landcover targets of different stand types and stands in a landcover)
	bool pool_to_all_landcovers[NLANDCOVERTYPES];

	/// Whether to pool transferred land to a receptor landcover (crop and pasture stands to new natural stand: pool!)
	bool pool_from_all_landcovers[NLANDCOVERTYPES];

	void serialize(ArchiveStream& arch);
};

/// The Gridcell class corresponds to a modelled locality or grid cell.
/** Member variables include an object of type Climate (holding climate, insolation and
 *  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
 *  array of Stand objects. Soil objects (holding soil state variables) are associated
 *  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
 *  locality or grid cell.
 */
class Gridcell : public GuessContainer<Stand>, public Serializable {

public:

	// MEMBER VARIABLES

	/// climate, insolation and CO2 for this grid cell
	Climate climate;

	/// soil static parameters for this grid cell
	Soiltype soiltype;

	/// landcover fractions and landcover-specific variables
	Landcover landcover;

	/// list array [0...npft-1] of Gridcellpft (initialised in constructor)
	ListArray_idin1<Gridcellpft,Pft> pft;

	/// list array [0...nst-1] of Gridcellst (initialised in constructor)
	ListArray_idin1<Gridcellst,StandType> st;

	/// object for keeping track of carbon, nitrogen and water balance
	MassBalance balance;

	/// Soil depth for gridcell (mm)
	double soil_depth;

	// SIMFIRE
	/// the region index to chosose from set of optimisations
	int simfire_region;
	/// timeseries of population density from the Hyde 3.1 dataset (inhabitants/km2)
	double hyde31_pop_density[57];
	/// current year's population density (inhabitants/km2)
	double pop_density;
	/// tuning factor for available litter
	double k_tun_litter;
	/// maximum annual Nesterov Index
	double nesterov_max;
	/// current Nexterov index
	double nesterov_cur;
	/// Monthly max Nexterov index (to keep track of running year)
	double nesterov_monthly_max[12];
	/// biome classification used in SIMFIRE
	int simfire_biome;
	/// Average maximum annual fAPAR (over avg_interv_fpar years)
	double fapar_annual_max;
	/// Average maximum annual fAPAR of recent years
	double fapar_recent_max[AVG_INTERVAL_FAPAR];
	/// maximum fapar of running year so far
	double fapar_cur_max;
	/// monthly fire risk (factor describing local monthly fire climatology)
	double monthly_fire_risk[12];
	/// current burned area from SIMFIRE (fract.)
	double burned_area;
	/// effectively burned area (after randFrac in BLAZE)
	double effective_burned_area;
	/// accumulated burned area from SIMFIRE for tstep < 1a (fract.)
	double simfire_annual_burned_area;
	/// Simple tracker to check whether at least one patch has enough fuel to burn
	int can_burn;
	/// annual burned area
	double annual_burned_area;
	/// monthly burned area
	double monthly_burned_area[12];

	// Nitrogen deposition
	/// annual NH4 deposition (kgN/m2/year)
	double aNH4dep;
	/// annual NO3 deposition (kgN/m2/year)
	double aNO3dep;
	/// daily NH4 deposition (kgN/m2)
	double dNH4dep;
	/// daily NO3 deposition (kgN/m2)
	double dNO3dep;

	double distinterval_gc;

	// Vector of (calendar) years in which a disturbance should occur.
	// This only contains data and is only used when disttype = DETERMINISTIC.
	std::set<int> static_disturbances;

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
	Gridcell();

	/// Longitude for this grid cell
	double get_lon() const;

	/// Latitude for this grid cell
	double get_lat() const;

	/// Set longitude and latitude for this grid cell
	void set_coordinates(double longitude, double latitude);

	/// Set initial soil layer properties
	void set_soillayer_properties();
	
	/// Calculates soil layer depths and number of soil layers
	void soillayer_thickness(double lyr_thickness[NSOILLAYERS_MAX], int& nsoillayers);

	void serialize(ArchiveStream& arch);

	/// Creates a new Stand in this grid cell
	Stand& create_stand(landcovertype lc, int no_patch = 0);

	/// Creates new stand and initiates land cover settings when run_landcover==true
	Stand& create_stand_lu(StandType& st, double fraction, int no_patch = 0, bool suppress_disturbance = false);

	/// Total gridcell carbon biomass and litter
	double ccont() const;
	/// Total gridcell nitrogen biomass and litter
	double ncont() const;
	/// Total gridcell water and ice
	double water_content() const;
	/// Total gridcell carbon fluxes so far this year
	double cflux(Period period) const;
	/// Total gridcell nitrogen fluxes so far this year
	double nflux(Period period) const;
	/// Total gridcell water fluxes so far this year
	double water_flux(Period period) const;

	/// Deletes the stand which the iterator is pointing at
	/** Returns an iterator pointing to the object following the erased object.
	 */
	iterator delete_stand(iterator itr);

	/// Returns number of stands
	unsigned int nbr_stands() const;

private:

	/// Longitude for this grid cell
	double lon;

	/// Latitude for this grid cell
	double lat;

};

/*
Calculate a value for A_growth such that the temperature response equation
described in <...> has a maximum value of 1.
*/
double calc_a_growth(const Pft *pft);

/*
Calculate the response window (in days) for phenological responses to growing
conditions.

@param nday_dormant: The number of dormant days in a year.
*/
size_t get_response_window(double dormancy_depth);

/*
Calculate the root beta value corresponding to the specified number of dormant
days in a year.
*/
double calc_root_beta(double wscal_avg);

int get_solvesomcent_endyr();

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
// Gerten, D., Schaphoff, S., Haberlandt, W., Lucht, W. & Sitch, S. 2004. 
//   Terrestrial vegetation and water balance—hydrological evaluation of a dynamic 
//   global vegetation model. Journal of Hydrology 286: 249-270.
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Jackson, R.B., Canadell, J., Ehleringer, J.R., Mooney, H.A., Sala O.E. & Schulze, E.D. 1996
//   A global analysis of root distributions for terrestrial biomes. 
//   Oecologia, Volume 108: 389–411
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Macduff, JH, Humphreys, MO & Thomas, H 2002. Effects of a stay-green mutation on
//   plant nitrogen relations in Lolium perenne during N starvation and after
//   defoliation. Annals of Botany, 89, 11-21.
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Olin S., G. Schurgers, M. Lindeskog, D. Wårlind, B. Smith, P. Bodin, J.
//    Holmér, and A. Arneth. 2015 Biogeosciences Discuss., 12, 1047-1111. The
//    impact of atmospheric CO2 and N management on yields and tissue C:N in
//    the main wheat regions of Western Europe
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W.
//   & Kelly, R. 2010. ForCent model development and testing using the Enriched
//   Background Isotope Study experiment. Journal of Geophysical
//   Research-Biogeosciences, 115.
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
// Wang, E, Engel, T, 1998 Simulation of phenological development of wheat crops
//   Agricultural Systems 58:1-24
// White, M A, Thornton, P E, Running, S. & Nemani, R 2000 Parameterization and
//   Sensitivity Analysis of the BIOME-BGC Terrestrial Ecosystem Model: Net Primary
//   Production Controls. Earth Interactions, 4, 1-55.
