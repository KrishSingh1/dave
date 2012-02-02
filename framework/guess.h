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
#include "shell.h"
#include "guessmath.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

typedef enum {NOLIFEFORM,TREE,GRASS} lifeformtype;
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

/// Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
typedef enum {URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, NLANDCOVERTYPES} landcovertype;

/// Water uptake parameterisations
/** \see water_uptake in canexch.cpp
  */
typedef enum {WR_WCONT, WR_ROOTDIST, WR_SMART, WR_SPECIESSPECIFIC} wateruptaketype;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

const int NSOILLAYER=2;
	// number of soil layers modelled
const double SOILDEPTH_UPPER=500.0; // soil upper layer depth (mm)
const double SOILDEPTH_LOWER=1000.0; // soil lower layer depth (mm)


	// guess2008 - new default SOM values
const int SOLVESOM_END=400;
	// year at which to calculate equilibrium soil carbon
const int SOLVESOM_BEGIN=350;
	// year at which to begin documenting means for calculation of equilibrium
	// soil carbon
const int NYEARGREFF=5;
	// number of years to average growth efficiency over in function mortality
const int COLDEST_DAY_NHEMISPHERE=14;
	// day at which to start counting GDD's and leaf-on days for summergreen phenology
	// in N hemisphere (January 15)
const int COLDEST_DAY_SHEMISPHERE=195;
	// day at which to start counting GDD's and leaf-on days for summergreen phenology
	// in S hemisphere (July 15)
const int OUTPUT_MAXAGECLASS=2000;
	// maximum number of age classes in age structure plots produced by function
	// outannual

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
class Gridcell;

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

/// Water uptake parameterisation
extern wateruptaketype wateruptake;

/// Whether other landcovers than natural vegetation are simulated.
extern bool run_landcover;

/// Whether a specific landcover type is simulated (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND).
extern bool run[NLANDCOVERTYPES];

/// Whether landcover fractions are read from ins-file.
extern bool lcfrac_fixed;

/// Set to false by initio( ) if fraction input files have yearly data.
extern bool all_fracs_const;

extern bool ifslowharvestpool; 	// If a slow harvested product pool is included in patchpft.
extern int nyear_spinup; // number of spinup years (ML)	Moved to guess.cpp to be accessed globally.

///////////////////////////////////////////////////////////////////////////////////////
// guess2008 - new input variables, from the .ins file
extern bool ifsmoothgreffmort;				
	// whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard 
	// step-function (0)
extern bool ifdroughtlimitedestab;			
	// whether establishment is limited by growing season drought 
extern bool ifrainonwetdaysonly;			
	// rain on wet days only (1, true), or a little every day (0, false); 
// bvoc
extern bool ifbvoc; 
        // whether BVOC calculations are included



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
		const int data[]={31,28,31,30,31,30,31,31,30,31,30,31};
		int month;
		int dayct=0;
		for (month=0;month<12;month++) {
			ndaymonth[month]=data[month];
			middaymonth[month]=dayct+data[month]/2;
			dayct+=data[month];
		}
	}

	/// Initialises date to day 0 of year 0 and sets intended number of simulation years
	/** Intended number of simulation years is only used to set islastyear flag,
	 *  actual simulation may be longer or shorter.
	 *
	 *  \param nyearsim  Intended number of simulation years
	 */
	void init(int nyearsim)	{
		nyear=nyearsim;
		day=month=year=dayofmonth=0;
		islastmonth=islastday=ismidday=false;
		if (nyear==1) islastyear=true;
		else islastyear=false;
	}

	/// Call at end of every simulation day to update member variables.
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

	// \returns index (0-11) of previous month (11 if currently month 0).
	int prevmonth() {
		if (month>0) return month-1;
		return 11;
	}

	/// \returns index of next month (0 if currently month 11)
	int nextmonth() {
		if (month<11) return month+1;
		return 0;
	}
};

/// This struct contains the result of a photosynthesis calculation.
/** \see photosynthesis */  
struct PhotosynthesisResult {
	/// Constructs an empty result
	PhotosynthesisResult() {
		clear();
	}

	/// Clears all members
	/** This is returned by the photosynthesis function when no photosynthesis
	 *  takes place.
	 */
	void clear() {
		agd_g      = 0.0;
		adtmm      = 0.0;
		rd_g       = 0.0;
		pi_co2_opt = 0.0;
		gammastar  = 0.0;
		apar       = 0.0;
		phi_pi     = 0.0;
	}

	/// gross daily photosynthesis (gC/m2/day)
	double agd_g;

	/// leaf-level net daytime photosynthesis 
	/** expressed in CO2 diffusion units (mm/m2/day) */
    double adtmm;

	/// leaf respiration (gC/m2/day)
	double rd_g;

	/// non-water-stressed intercellular partial pressure of CO2 (Pa)
	double pi_co2_opt;

	/// CO2 compensation point in partial pressure units (Pa)
    double gammastar;

	/// amount of PAR absorbed at leaf level (J m-2 d-1)
    double apar;

	/// factor accounting for effect of intercellular CO2 concentration on C4 photosynthesis
    double phi_pi;

	/// gross daily photosynthesis (kgC/m2/day)
    double agd() const {
        return agd_g/1000.0;
    }

	/// leaf-level net daytime photosynthesis (kgC/m2/day)
    double rd() const {
        return rd_g/1000.0;
    }

	/// net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
    double net_assimilation() const {
        return agd()-rd();
    }
};


/// The Climate for a grid cell
/** Stores all static and variable data relating to climate parameters, as well as 
 *  latitude, atmospheric CO2 concentration and daylength for a grid cell. Includes 
 *  a reference to the parent Gridcell object (defined below). Initialised by a 
 *  call to initdrivers.
 */
class Climate {

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
	double insol;
		// insolation today, see also instype
		// When instype is NETSWRAD or SWRAD insol is assumed to be W/m2 during
		// daylight hours. If input data is averaged over a 24 hour period, code
		// dealing with this variable needs to be changed 
		// (see function daylengthinsoleet).
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
	double qo[365],u[365],v[365],hh[365],sinehh[365];
	double daylength_save[365];
	bool doneday[365];
		// indicates whether saved values exist for this day

	// bvoc
	double dtr; // diurnal temperature range (oC)

	// MEMBER FUNCTIONS

public:

	Climate(Gridcell& gc):gridcell(gc) {};
		// constructor function: initialises gridcell member

	void initdrivers(double latitude) {

		// Initialises certain member variables
		// Should be called before Climate object is applied to a new grid cell

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
	// (all CO2 fluxes on stand area basis, kgC/m2 ;
        // BVOC fluxes (isoprene and monoterpenes) in gC/m2)

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
	double acflux_harvest; 
		// annual flux to atmosphere from consumed harvested products
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
	// bvoc
	double miso[12];
                // monthly isoprene flux (g C/m2/month)
	double mmon[12];
	        // monthly monoterpene flux (g C/m2/month)



	// MEMBER FUNCTIONS

public:
	/// constructor: initialises members
	Fluxes(Patch& p):patch(p) {
		acflux_veg=0.0;
		acflux_fire=0.0;
		acflux_soil=0.0;
		acflux_est=0.0;
		acflux_harvest=0.0;	
	}
		

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
	double cton_leaf;
		// leaf C:N mass ratio
	double cton_root;
		// fine root C:N mass ratio
	double cton_sap;
		// sapwood C:N mass ratio
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
	double Y_eps_iso;
	        // fraction of electron transport to isoprene production under standard conditions (-)
	bool seas_iso; 
	        // whether (1) or not (1) isoprene emissions show a seasonality
	double eps_mon;
	        // monoterpene emission capacity (ug C g-1 h-1)
	double Y_eps_mon;
	        // fraction of electron transport to monoterpene production under standard conditions (-)
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
	landcovertype landcover;

	double res_outtake;				// Fraction of residue outtake at harvest.
	double harv_eff;				// Harvest efficiency.
	double harvest_slow_frac;		// Fraction of harvested products that goes into patchpft.harvested_products_slow
	double turnover_harv_prod;		// Yearly turnover fraction of patchpft.harvested_products_slow (goes to fluxes.acflux_harvest).
	
	// MEMBER FUNCTIONS

public:

	Pft() {
		
		// Constructor (initialises array gdd0)
		
		int y;
		for (y=0;y<366;y++)
			gdd0[y]=-1.0; // value<0 signifies "unknown"; see function phenology()

		// guess2008 - DLE
		drought_tolerance=0.0; // Default, means that the PFT will never be limited by drought.

		res_outtake=0.0;
		harv_eff=0.0;
		turnover_harv_prod=1.0;	// default 1 year turnover time
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
		else if (lifeform==GRASS) {

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

	bool alive; 
		// guess2008 - whether this individual is truly alive. Set to false for first year 
		// after the Individual object is created, then true.


	// bvoc
	double iso; // isoprene production (mg C m-2 d-1)
	double mon; // monoterpene production (mg C m-2 d-1)
	double aiso; // annual isoprene emission (mg C m-2 y-1)
	double amon; // annual monoterpene emission (mg C m-2 y-1)
	double monstor; // monoterpene storage pool (mg C m-2)
	double fvocseas; // isoprene seasonality factor (-)
	double dtr_wstress; // diurnal temperature range (oC)
	double eet_wstress; // equilibrium evapotranspiration today (mm/day)
	double agdd5_wstress; // total gdd5 (accumulated) for this year (reset 1 January)
	double rad_wstress; // total daily net downward shortwave solar radiation today (J/m2/day)

	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v);
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
	int solvesom_end;
		// year at which to calculate equilibrium soil carbon
	int solvesom_begin;
		// year at which to begin documenting means for calculation of equilibrium
		// soil carbon

	// MEMBER FUNCTIONS

public:

	Soiltype() {

		// Constructor: initialises certain member variables

		solvesom_end=SOLVESOM_END;
		solvesom_begin=SOLVESOM_BEGIN;
	}

	// guess2008 - override the default SOM years with 70-80% of the spin-up period length
	void updateSolveSOMvalues(const int& nyrspinup) {
		
		solvesom_end=static_cast<int>(0.8*nyrspinup);
		solvesom_begin=static_cast<int>(0.7*nyrspinup);

	}
};


/// Soil stores state variables for soils and the snow pack. 
/** Initialised by a call to initdrivers. One Soil object is defined for each patch. 
 *  A reference to the parent Patch object (defined below) is included as a member 
 *  variable. Soil static parameters are stored as objects of class Soiltype, of which 
 *  there is one for each grid cell. A reference to the Soiltype object holding the 
 *  static parameters for this soil is included as a member variable.
 */
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

	double rain_melt;						// rainfall and snowmelt today (mm)
	double max_rain_melt;					// upper limit for percolation (mm)
	bool percolate;							// whether to percolate today

	// MEMBER FUNCTIONS

public:
	/// constructor (initialises member variable patch)
	Soil(Patch& p,Soiltype& s):patch(p),soiltype(s) {
			initdrivers();
	}

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
		}

	}

};

///////////////////////////////////////////////////////////////////////////////////////
// CLASS LOOKUP_LAMBDA
// Lookup table for photosynthesis parameters (required for "fast" version of canopy
// exchange code; see canexch.cpp)

const int LOOKUP_LAMBDA_MAXITEM=130;

struct Lookup_lambda_item {
	PhotosynthesisResult photosynthesis;
	int year;
	int day;

	Lookup_lambda_item()
			: photosynthesis(), year(-1), day(0) {
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

	bool getdata(int year,int day,PhotosynthesisResult& photosynthesis) {
		if (position>=LOOKUP_LAMBDA_MAXITEM)
			fail("class Lookup_lambda: exceeded dimension of lookup table");
		Lookup_lambda_item& thisitem=data[position];
		if (thisitem.year==year && thisitem.day==day) {
			photosynthesis = thisitem.photosynthesis;
			return true;
		}
		// else
		return false;
	}

	void setdata(int year,int day, const PhotosynthesisResult& photosynthesis) {
		if (position>=LOOKUP_LAMBDA_MAXITEM)
			fail("class Lookup_lambda: exceeded dimension of lookup table");
		Lookup_lambda_item& thisitem=data[position];
		thisitem.year=year;
		thisitem.day=day;
		thisitem.photosynthesis = photosynthesis;
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
	double harvested_products_slow;	//carbon depository for long-lived products like wood


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
		harvested_products_slow=0.0;
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
	double gterm;
		// term in calculation of potential canopy conductance (mm/s)
	bool have_gterm;
		// true if value of gterm available for this PFT today, otherwise false

	// Variables used by "fast" canopy exchange code (Ben Smith 2002-07)

	double gpterm;
		// non-FPAR-weighted value for canopy conductance component associated with
		// photosynthesis for PFT under non-water-stress conditions (mm/s)
	double assim_term;
		// non-FPAR-weighted leaf-level net photosynthesis value for PFT under non-
		// water-stress conditions (kgC/m2/day)
	double fpc_total;
		// FPC sum for this PFT as average for stand (used by some versions of
		// guessio.cpp)

	/// Photosynthesis values for this PFT under non-water-stress conditions
	PhotosynthesisResult photosynthesis;
	
	

	/// Is this PFT allowed to grow in this stand ?
	bool active;

	// MEMBER FUNCTIONS

	Standpft(int i,Pft& p):id(i),pft(p) {
		
		// Constructor: initialises various data members
		
		anetps_ff_max=0.0;

		if (run_landcover)
			active=false;
		else
			active=true;
	}
};


/// The stand class corresponds to a modelled area of a specific landcover type in a grid cell.
/** There may be several stands of the same landcover type (but with different settings).
 */
class Stand : public ListArray_idin3<Patch,Stand,Pftlist,Soiltype> {

public:

	// MEMBER VARIABLES

	/// list array [0...npft-1] of Standpft (initialised in constructor)
	ListArray_idin1<Standpft,Pft> pft;

	/// A number identifying this Stand within the grid cell
	int id;

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


	// MEMBER FUNCTIONS

	/// Constructs a Stand
	/** \param i         The id for the stand within the grid cell
	 *  \param gc        The parent grid cell
	 *  \param landcover The type of landcover to use for this stand
	 *  \param pftlist   The list of PFTs
	 */
	Stand(int i, Gridcell& gc,landcovertype landcover,Pftlist& pftlist); 

	/// Gives the fraction of this Stand relative to the whole grid cell
	double get_gridcell_fraction() const;

	/// Gives the fraction of this Stand relative to its land cover type
	double get_landcover_fraction() const;

	/// Set the fraction of this Stand relative to its land cover type
	void set_landcover_fraction(double fraction);

private:

	/// Fraction of this stand relative to its landcover
	/** used by crop stands; initialized in constructor to 1, 
	 *  set in landcover_init() 
	 */
	double frac;
};



/// State variables common to all individuals of a particular PFT in a GRIDCELL.
class Gridcellpft {

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

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellpft object
	/** \param i   The id for this object
	 *  \param p   A reference to the Pft for this Gridcellpft
	 */
	Gridcellpft(int i,Pft& p):id(i),pft(p) {
		addtw=0.0;
	}
};


/// The Gridcell class corresponds to a modelled locality or grid cell.
/** Member variables include an object of type Climate (holding climate, insolation and
 *  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
 *  array of Stand objects. Soil objects (holding soil state variables) are associated
 *  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
 *  locality or grid cell.
 */
class Gridcell : public ListArray_idin3<Stand,Gridcell,landcovertype, Pftlist> {

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

	/// The land cover fractions from the previous year
	/** Used to keep track of the changes when running with dynamic
	 *  land cover.
	 */
	double landcoverfrac_old[NLANDCOVERTYPES];

	/// Whether the land cover fractions changed for this grid cell this year
	/** \see landcover_dynamics
	 */
	bool LC_updated;

	/// list array [0...npft-1] of Gridcellpft (initialised in constructor)
	ListArray_idin1<Gridcellpft,Pft> pft;
    

	// MEMBER FUNCTIONS

	/// Constructs a Gridcell object
	/** \param pftlist    The list of plant functional types
	 */
	Gridcell(Pftlist& pftlist):climate(*this) {
		landcovertype landcover;
		LC_updated=false;

		for(unsigned int p=0;p<pftlist.nobj;p++) {
			Gridcellpft& gcpft=pft.createobj(pftlist[p]);
		}		

		memset(landcoverfrac, 0, sizeof(double)*NLANDCOVERTYPES);
		memset(landcoverfrac_old, 0, sizeof(double)*NLANDCOVERTYPES);

		if(!run_landcover) {
			landcover=NATURAL;
			createobj(*this,landcover,pftlist);
			landcoverfrac[NATURAL]=1.0;
		}
	}
};


#endif // LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
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
