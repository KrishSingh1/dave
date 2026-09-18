///////////////////////////////////////////////////////////////////////////////////////
/// \file parameters.h
/// \brief The parameters module is responsible for reading in the instruction file
///
/// This module defines and makes available a lot of the instruction file parameters
/// used by the model, but also lets other modules define their own parameters or
/// access "custom" parameters without defining them.
///
/// A new parameter can be added by creating a new global variable here (or a new
/// Pft member variable if it's a PFT parameter), and then declaring it in
/// plib_declarations in parameters.cpp. See the many existing examples, and
/// documentation in the PLIB library for further documentation about this.
///
/// Sometimes, adding a new parameter shouldn't (or can't) be done here however.
/// A parameter specific for a certain input module, should only be declared if
/// that input module is used. In this case the input module should declare its
/// own parameters when it is created. This can also be a good idea simply to
/// make modules more independent. For parameters like this, we can either use
/// the "custom" parameters (\see Paramlist) which don't need to be declared at
/// all, or the parameters can be declared with the declare_parameter family of
/// functions.
///
/// \author Joe Siltberg
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_PARAMETERS_H
#define LPJ_GUESS_PARAMETERS_H

#include "gutil.h"
#include <string>


///////////////////////////////////////////////////////////////////////////////////////
// Enums needed by some of the global instruction file parameters defined below


/// Vegetation 'mode', i.e. what each Individual object represents
/** Can be one of:
 *  1. The average characteristics of all individuals comprising a PFT
 *     population over the modelled area (standard LPJ mode)
 *  2. A cohort of individuals of a PFT that are roughly the same age
 *  3. An individual plant
 */
typedef enum {NOVEGMODE, INDIVIDUAL, COHORT, POPULATION} vegmodetype;

/// Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
/*  NB. set_lc_change_array() must be modified when adding new land cover types
 */
typedef enum {URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, BARREN, NLANDCOVERTYPES} landcovertype;

// TODO: reconcile with our disturbancetype
typedef enum {CLEARCUT, CLEARCUT_PRESCRIBED, FIRE, WIND, BARKBEETLE, NDISTURBANCETYPES} disturbancetype;

/// Water uptake parameterisations
/** \see water_uptake in canexch.cpp
  */
typedef enum {WR_WCONT, WR_ROOTDIST, WR_SMART, WR_SPECIESSPECIFIC, WR_ROOT_RESISTANCE} wateruptaketype;

///bvoc: define monoterpene species used
typedef enum {APIN, BPIN, LIMO, MYRC, SABI, CAMP, TRIC, TBOC, OTHR, NMTCOMPOUNDTYPES} monoterpenecompoundtype;

/// Fire model setting. Either use 
/**	One of
 *	BLAZE 		Use the BLAZE model to generate fire fluxes 
 *                      (must be accompanied by ignitionmode; DEFAULT)
 *	GLOBFIRM	fire parameterization following Thonicke et al. 2001
 *	NOFIRE		no fire model	
 */
typedef enum {BLAZE, GLOBFIRM, NOFIRE} firemodeltype;

/// SIMFIRE version setting. Either use 
/**	One of
 *	EUROPE 		Use European parameterization from 2024 
 *	GLOBAL		Use global version (as in trunk)
 */
typedef enum {SIMFIRE_GLOBAL, SIMFIRE_EUROPE} simfireversiontype;

/// Type of weathergenerator used 
/**     One of:
 *      GWGEN           Global Weather GENerator (needed by BLAZE, due to 
 *                      additional rel. humidity and wind; DEFAULT)
 *      INTERP          use standard interpolation scheme
 *      NONE            Should be set if daily input is used (e.g. in cfinput) 
 */
typedef enum {GWGEN, INTERP, NONE} weathergeneratortype;

/// Type of hydrology used 
/**     One of:
 *      GERTEN          Standard two-layer hydrology (0.5 and 1m layers) based on Gerten et al. (2004)
 *      MILLER          Extention of GERTEN where each layer is divided into 10 cm layers. What distribution
 *						within the two basic layers is based on Miller et al. (in prep.)
 *      RICHARDS        Flexible hydrology scheme based on Richards (1931)
 */
typedef enum { GERTEN, MILLER, RICHARDS } hydrologytype;

/// Type of energy calculation used 
/**     One of:
 *      CNSTEP_FULL
 *      CNSTEP
 *      BELDA
 */
typedef enum { CNSTEP_FULL, CNSTEP, BELDA} energytype;

/// Type of water modifier for decomposition of soil organic matter used 
/**     One of:
 *      FRIEND          Standard soil moisture decomposition function based on Friend et al (1997)
 *      YAN             A moisture function of soil heterotrophic respiration that incorporates 
 *                      microscale processes based on Yan et al (2018)
 */
typedef enum { FRIEND, YAN } moistdecaytype;

/// The type of disturbance regime.
/**		One of:
 * 		NONE			No disturbances enabled (DEFAULT)
 * 		STOCHASTIC		Random disturbances, varying by patch with mean interval
 *  					defined by distinterval (ie trunk disturbance behaviour)
 * 		DETERMINISTIC	Disturbances occurring in predetermined years specified
 * 						by file_disturb
 */
typedef enum {DIST_NONE, STOCHASTIC, DETERMINISTIC} davedisturbancetype;

///////////////////////////////////////////////////////////////////////////////////////
// Global instruction file parameters

/// Title for this run
extern xtring title;

/// Whether to produce daily output
extern bool ifdailyoutput; 

/*
True iff outputs should be written during spinup.
*/
extern bool spinup_output;

/// Vegetation mode (population, cohort or individual)
extern vegmodetype vegmode;

/// Default number of patches in each stand
/** Should always be 1 in population mode,
 *  cropland stands always have 1 patch.
 *  Actual patch number for stand objects may differ and 
 *  should always be queried by stand.npatch()
 */
extern int npatch;

/// Number of patches in each stand for secondary stands
extern int npatch_secondarystand;

/// Number of patches in stands created after prescribed disturbance (if not set, npatch_secondarystand is used)
extern int npatch_disturbedstand;

/// Number of patches simulated in stands created at whole tree initialisation with option forest_init_poolsoil
extern int npatch_treeinit;

/// Whether to reduce equal percentage of all stands of a stand type at land cover change
extern bool reduce_all_stands;

/// Minimum age of stands to reduce at land cover change
extern int age_limit_reduce;

/// Patch area (m2) (individual and cohort mode only)
extern double patcharea;

/// Whether background establishment enabled (individual, cohort mode)
extern bool ifbgestab;

/// Whether spatial mass effect enabled for establishment (individual, cohort mode)
extern bool ifsme;

/// Whether establishment stochastic (individual, cohort mode)
extern bool ifstochestab;

/// Whether mortality stochastic (individual, cohort mode)
extern bool ifstochmort;

/// Fire-model switch
extern firemodeltype firemodel;

/// SIMFIRE version switch
extern simfireversiontype burnedareamodel;

/// Startyear for prescribed burned area
extern int presc_burned_area_start;

/// Endyear for prescribed burned area
extern int presc_burned_area_end;

/// Weather Generator switch
extern weathergeneratortype weathergenerator;

/// Hydrology switch
extern hydrologytype hydrology;

/// Energy calculation switch
extern energytype energy;

/// Mater modifier for decomposition of soil organic matter switch
extern moistdecaytype moistdecay;

/// Whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
extern bool ifdisturb;

/// Generic patch-destroying disturbance interval (individual, cohort mode)
extern double distinterval;

/// The disturbance algorithm being used.
extern davedisturbancetype disttype;

/// Path to a file containing years in which disturbances occur, one per line.
extern xtring file_disturb;

/// Whether SLA calculated from leaf longevity (alt: prescribed)
extern bool ifcalcsla;

/// Whether vcmax acclimation to temperature is enabled (alt: prescribed)
extern bool ifvcmaxtempacclimation;

/// Whether maintenance respiration acclimation to temperature is enabled.
extern bool ifresptempacclimation;

/// Whether leaf damage due to high temperature is enabled.
extern bool ifleafdamage;

/// Whether leaf C:N ratio minimum calculated from leaf longevity (alt: prescribed)
extern bool ifcalccton;

/// Establishment interval in cohort mode (years)
extern int estinterval;

/// Whether C debt (storage between years) permitted
extern bool ifcdebt;

/// Water uptake parameterisation
extern wateruptaketype wateruptake;

/// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
extern bool ifcentury;

/// whether plant growth limited by available N
extern bool ifnlim;

/// Whether daily leaf N retranslocation is allowed
extern bool ifdailynretrans;

/// number of years to allow spinup without nitrogen limitation
extern int freenyears;

/// fraction of nitrogen relocated by plants from roots and leaves
extern double nrelocfrac;

/// first term in nitrogen fixation eqn (Cleveland et al 1999)
extern double nfix_a;

/// second term in nitrogen fixation eqn (Cleveland et al 1999)
extern double nfix_b;

/// whether to use nitrification/denitrification in CENTURY SOM dynamics
extern bool ifntransform;
/// Fraction of microbial respiration assumed to produce DOC, 0.0,0.3
extern double frac_labile_carbon;

/// Soil pH (used for calculating N-transformation), 3.5,8.5
extern double pH_soil;
/// Maximum nitrification rate, 0.03,0.15
extern double f_nitri_max;
/// Denitrification Michaelis–Menten constant for N oxides [kg N kg H2O-1], 0.000001,0.0001
extern double k_N;
/// Denitrification Michaelis–Menten constant for labile carbon [kg C kg H2O-1], 0.000005,0.00005
extern double k_C;
/// Maximum gaseus losses in nitrification
extern double f_nitri_gas_max;
/// Maximum fraction of NO3 converted to NO2
extern double f_denitri_max;
/// Maximum fraction of NO2 converted to gaseus N
extern double f_denitri_gas_max;

// Mapping of text input file data in index file.
extern bool map_text_file;

/// Whether to use the new allocation calculations (1), or the original (0)
extern bool ifnewalloc;
/// Whether to use the new tree volume to height calculations (1), or the original (0)
extern bool ifvolhght;
/// Whether to print which mortality function was triggered (1) or not (0)
extern bool ifprintallommort;

/// Number of vertical layers used for canopy-area integration.
/// Mutually exclusive with layer_thickness_canopy.
extern int nlayer_canopy;

/// Width of vertical layers for canopy-area integration (m) (default: 2).
/// Mutually exclusive with nlayer_canopy.
extern double layer_thickness_canopy;

///////////////////////////////////////////////////////////////////////////////////////
// Landuse and crop settings

/// Whether other landcovers than natural vegetation are simulated.
extern bool run_landcover;

/// Whether a specific landcover type is simulated (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, BARREN).
extern bool run[NLANDCOVERTYPES];

/// Whether landcover fractions are not read from input file.
extern bool lcfrac_fixed;

/// Whether peat durations are read from input file.
extern bool peatduration;

/// Whether fractions of stand types of a specific land cover are not read from input file.
extern bool frac_fixed[NLANDCOVERTYPES];

/// Whether BARREN landcover excluded from area fraction correction in cases of non-unity sum
extern bool no_barren_frac_corr;

/// Whether PEATLAND landcover excluded from area fraction correction in cases of non-unity sum
extern bool no_peatland_frac_corr;

/// Set to false by initio( ) if fraction input files have yearly data.
extern bool all_fracs_const;

/// Whether a fraction of harvested wood is put into a product pool
extern bool ifslowharvestpool;

/// Whether grass is allowed to grow between crop growingseasons
extern bool ifintercropgrass;

/// Whether to calculate dynamic potential heat units
extern bool ifcalcdynamic_phu;

/// Whether to use gross land transfer: read landcover transfer matrix input file (1), read stand type transfer matrix input file (2), or not (0)
extern int gross_land_transfer;

/// Whether gross land cover transfer input read for this gridcell
extern bool gross_LC_input_present;

/// Whether gross stand type transfer input read for this gridcell
extern bool gross_ST_input_present;

/// Whether to use primary/secondary land transition info in landcover transfer input file (1). or not (0)
extern bool ifprimary_lc_transfer;

/// Distinguish between primary and secondary natural stands at area reduction
extern bool use_primary_lc_transfer;

/// Whether to use primary-to-secondary land transition info (within land cover type) in landcover transfer input file (1). or not (0)
extern bool ifprimary_to_secondary_transfer;

/// Pooling level of land cover transitions; 0: one big pool; 1: land cover-level; 2: stand type-level
extern int transfer_level;

/// Whether to create new stands in transfer_to_new_stand() according to the rules in copy_stand_type()
extern bool iftransfer_to_new_stand;

/// Whether to suppress disturbance and fire in forestry stands created as NATURAL stands in transfer_to_new_stand_from_stand() or transfer_to_new_stand_from_st_lc() (eg. LUH2 input)
extern bool suppress_disturbance_in_forestry_stands;

/// Whether to harvest (remove) wood at natural-to-forest transitions
extern bool harvest_natural_to_forest;

/// Whether to limit dynamic phu calculation to a period specified by nyear_dyn_phu
extern bool ifdyn_phu_limit;

// Whether to use closest available n deposition data, should a particular
// gridcell be missing (ie when not running a 0.25 degree grid).
extern bool useClosestNdep;

// Number of years to calculate dynamic phu if dynamic_phu_limit is true
extern int nyear_dyn_phu;

/// number of spinup years
extern int nyear_spinup;

/// Whether to use sowingdates from input file
extern bool readsowingdates;

/// Whether to use harvestdates from input file
extern bool readharvestdates;

/// Whether to read N fertilization from input file
extern bool readNfert;

/// Whether to read grazing dates from input file.
extern bool read_grazing;

/// Whether to read manure N fertilization from input file
extern bool readNman;

/// Whether to read N fertilization (stand type level) from input file
extern bool readNfert_st;

/// Whether to use forest harvested fraction from input file (using LUC functionality)
extern bool readwoodharvest_frac;

/// Whether to use wood harvest C mass from input file (using LUC functionality)
extern bool readwoodharvest_cmass;

/// Whether to use wood harvest probability for broadleaves from input file
extern bool readwoodharv_prob_BL;
/// Whether to use wood harvest probability for needleleaves from input file
extern bool readwoodharv_prob_NL;
/// Whether to use clearcut probability for broadleaves from input file
extern bool readCC_prob_BL;
/// Whether to use clearcut probability for needleleaves from input file
extern bool readCC_prob_NL;
/// Whether to use mean wood harvest intensity for broadleaves from input file
extern bool readwoodharv_intens_mean_BL;
/// Whether to use mean wood harvest intensityfor needleleaves from input file
extern bool readwoodharv_intens_mean_NL;
/// Whether to use mean wood harvest intensity standard deviation for broadleaves from input file
extern bool readwoodharv_intens_sd_BL;
/// Whether to use mean wood harvest intensity standard deviation for needleleaves from input file
extern bool readwoodharv_intens_sd_NL;
/// Whether to use thinning size preference for broadleaves from input file
extern bool readwoodharv_size_preference_BL;
/// Whether to use thinning size preference for needleleaves from input file
extern bool readwoodharv_size_preference_NL;
/// Whether to read probabilistic wood harvest intensity scalar from input file
extern bool readwoodharv_intensity_scalar;
/// Whether to draw thinning probabilities from a normal distribution
extern bool ifthin_prob_distribution;
/// Whether to turn off clearcuts while using clearcut probabilities input to bypass thinning
extern bool clearcutprob_off;
/// First year to use clearcut probabilities input
extern int clearcutprob_input_firstyear;
/// Last year to use clearcut probabilities input
extern int clearcutprob_input_lastyear;
/// Number of basal area bins in harvest probablilty input
extern int nBA_bins;
/// Size of basal area bins in harvest probablilty input
extern int BA_bin_size;
/// Value of first (smallest) basal area bin in harvest probablilty input
extern int BA_bin_min;
/// Number of diameter bins in harvest probablilty input
extern int ndiam_bins;
/// Size of diameter bins in harvest probablilty input
extern int diam_bin_size;
/// Value of first (smallest) diameter bin in harvest probablilty input
extern int diam_bin_min;
/// Whether last (largest) bin includes all larger values
extern bool iflast_bin_inclusive;
/// Whether forest initialisation uses biomass rather than diameter from input
extern bool forest_init_from_biomass;

/// Whether generic patch-destroying disturbance should be enacted as wood harvest (0,1)
extern bool ifback_disturb_as_harv;

/// Whether to create new stands at clearcut of secondary stands when using wood harvest input (LUC functionality)
extern bool harvest_secondary_to_new_stand;

/// Whether to read disturbance intervals from input file
extern bool readdisturbance;

/// Whether to read disturbance intervals for stand types from input file
extern bool readdisturbance_st;

/// Whether to read windload from input file
extern bool readwindload;

/// Whether to read cutinterval for stand types from input file
extern bool readcutinterval_st;

/// Whether to read gridcell mean elevation from input file
extern bool readelevation_gc_mean;

/// Whether to read local (site) elevation from input file
extern bool readelevation_local;

/// Whether to read stand type elevation from input file
extern bool readelevation_st;

/// Whether to read clearcut disturbance biomass fraction from input file
extern bool readdist_clearcut;

/// Whether to read fire disturbance biomass fraction from input file
extern bool readdist_fire;

/// Whether to read wind disturbance biomass fraction from input file
extern bool readdist_wind;

/// Whether to read barkbeetle disturbance biomass fraction from input file
extern bool readdist_barkbeetle;

extern int disturbance_input_firstyear;
extern int disturbance_input_lastyear;

/// Whether to read wind disturbance probability from input file
extern bool readstorm_prob;

/// Whether to read burned area from input file
extern bool readpresc_burned_area;

extern bool readforest_init;

extern bool forest_init_poolsoil;

/// Whether to read firstmanageyear for stand types from input file
extern bool readfirstmanageyear_st;

// Whether to read target-cutting distribution for selection in mt from input file
extern bool readtargetcutting;

/// Whether to burn thin trees during tree harvest (ignoring pft.harvest_slow_frac)
extern bool harvest_burn_thin_trees;

/// Whether to print multiple stands within a stand type (except cropland) separately
extern bool printseparatestands;

/// Whether to simulate tillage by increasing soil respiration
extern bool iftillage;

/// Use silt/sand fractions per soiltype
extern bool textured_soil;

/// Whether pastures are affected by disturbance and fire (affects pastures' npatch)
extern bool disturb_pasture;

/// Whether to simulate cropland as pasture
extern bool grassforcrop;

/// source for the fire population density: 1=simfire.bin, 2=netCDF within cfxinput
extern int fire_popdens_method;

/// Keep popdens level at fixed_popdens_year (only applies within cfxinput)
// 1=popdens start in year fixed_popdens_year, 2=allways use popdens of fixed_popdens_year
extern int fixed_popdens_hist;

/// Year of popdens to use if fixed_popdens_hist>0 (only applies within cfxinput)
extern int fixed_popdens_year;

/// Soil depth
extern double depth_soil;

///////////////////////////////////////////////////////////////////////////////////////
// Settings controlling the saving and loading from state files

/// Location of state files
extern xtring state_path;

/// Whether to restart from state files
extern bool restart;

/// Whether to save state files
extern bool save_state;

/// Save/restart year
extern int state_year;

/// The level of verbosity
extern int verbosity;

/// whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard step-function (0)
extern bool ifsmoothgreffmort;

/// whether establishment is limited by growing season drought
extern bool ifdroughtlimitedestab;

/// rain on wet days only (1, true), or a little every day (0, false);
extern bool ifrainonwetdaysonly;

/// whether BVOC calculations are included
extern bool ifbvoc;

/// Storm damage settings and parameters
extern bool storm_mortality; // Apply calculated storm damage to mortality
extern int storm_mortality_firstyear; // Apply calculated storm damage to mortality from specified year when storm_mortality is true
extern int storm_mortality_lastyear; // Apply calculated storm damage to mortality to specified year when storm_mortality is true
extern bool insect_mortality; // Apply calculated insect damage to mortality
extern int insect_mortality_firstyear; // Apply calculated insect damage to mortality from specified year when insect_mortality is true
extern int insect_mortality_lastyear; // Apply calculated insect damage to mortality to specified year when insect_mortality is true
extern double m_stormfactor; // Callibration factor for storm damage from sensitivity and WL, Eq. 15 Lagergren et al
extern bool print_storm_stand;

/// whether daily carbon allocation for grasses is turned on
extern bool ifdailygrass;

/// whether daily carbon allocation for trees is turned on.
extern bool ifdailytree;

/// whether DAVE (empirical) photosynthesis scheme is turned on
extern bool ifdavephotosyn;

/// Turn on/off tree phenology (Dave)
extern bool iftreephen;

/// Turn on/off daily carbon starvation mortality in trees (Dave)
extern bool ifdailycstarvmort;

/// Turn on/off grote tree allocation.
extern bool ifgrote;

/// whether DAVE growth routine for grasses is turned on
extern bool ifdavegrowth;

/// whether DAVE hydrology changes are used (e.g. allow water content to drop below water potential = -1.5 MPa, allow percolation between layers on days without rain)
extern bool ifdavehydrology;

/// Whether to use penman-monteith equation to calculate transpiration, or the
/// relationship with VPD.
extern bool penman_monteith;

/// Whether establishment should use live reproductive biomss (true) or total
/// allocation to reproduction (false) for spatial mass effect.
extern bool est_live_repr;

/// Depth of sublayer at top of upper soil layer, from which evaporation is
/// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
/// Must be a multiple of Dz_soil
extern int soildepth_evap;

/// whether DAVE hydraulics scheme is turned on
extern bool ifdavehydraulics;
/// whether Zhou response of photosynthesis to water stress is turned on (requires ifdavehydrology).
extern bool ifzhou;
/// Whether LAI is distributed evenly across canopy layers for grasses.
extern bool ifevenlaidist;

/// Whether root distribution and response window are adaptive to the environment.
/// (Number of dormant+dead days)
extern bool ifadaptive_rootdist;

/// Whether tree N uptake should go into long-term storage pool.
extern bool ifstoreuptake;

/// Resolution for splines used to calculate psi from transpiration.
/// This is the number of control points used for the integration. Default: 100.
extern double epsiresolution;

/// whether approximated solution for C4 photosynthesis is used
extern bool ifapproximate_C4;

/// whether to compress output files
extern bool gzip_output;

/// Maximum number of cohorts (individuals) allowed per pft (0 = unlimited).
extern int nindiv_max;

/// Maximum number of iterations to use in DAVE energy balance.
extern int energy_balance_maxiter;

/// Convergence criterion for DAVE energy balance (K). Iteration for tleaf will
/// cease when the change in tleaf between iterations is less than this value.
extern double energy_balance_eps;

/// Number of days that animals are allowed to graze a particular set of
/// patches (default: 1).
extern int grazing_duration;
/// True if we are using forced LAI.
extern bool ifforcedlai;
/// The fixed LAI value (m2/m2) applied to all cohorts when ifforcedlai is 1.
extern double fixed_lai;

///////////////////////////////////////////////////////////////////////////////////////
// Soil, arctic and wetland inputs

/// Whether to use fixed slope value (>=0 degree) or input data (<0)
extern double slope;

/// Variable infiltration capacity parameter (Global value 0.2, higher value gives lower infiltration and more surface runoff)
extern double b_infilt;

/// Calibration parameter for lateral drainage (-)
extern double C_dw;

/// Whether to consider bedrock as lower boundary condition for multilayer soil hydrology (1), or not (0)
extern bool ifbedrock;

/// Whether to consider supercooled soil water. 0 - not, 1 - Niuet al.2006 (CLM), 2 - Jin et al. 2020
extern int ifsupercooled;

/// Fraction of fine roots to be homogenously distributed in soil column
extern double homrootdistfrac;

/// Fraction of deep soil, beyond 1500mm, that roots have access to
extern double rootfrac_of_deep_soil;

/// Fraction of organic matter to move from receiver layer
extern double receiver_frac;

/// Use multilayer snow scheme, or the original LPJ-GUESS v4 scheme
extern bool ifmultilayersnow;

/// Whether to soak snowpack with melted water (1), or not (0)
extern bool ifsoaksnow;

/// Additional thickness ontop of snowpack to account for radiation energy
extern double top_thickness;

/// whether to reduce GPP if there's inundation (1), or not (0)
extern bool ifinundationstress;

/// Whether to limit soilC decomposition below 0 degC in upland soils (1), or not (0)
extern bool ifcarbonfreeze;

/// Extra daily water input or output, in mm, to wetlands. Positive values are run ON, negative run OFF.
extern double wetland_runon;

/// Whether methane calculations are included
extern bool ifmethane;

/// Age of peat (0 peatland will be spun up to equilibrium)
extern int peat_age;

/// Whether to update physical soil properties using soil carbon content (1), or not (0)
extern bool iforganicsoilproperties;

/// Whether to take water from runoff to saturate low latitide wetlands
extern bool ifsaturatewetlands;

/// Whether soil layer volume is affected by coarse root biomass taking up space (1), or not (0)
extern bool ifeffvolume;

/// Whether to use depth soil organic matter decay rate modifier
extern bool ifdepthdecaymod;

/// Whether thaw depth affects root distribution
extern bool ifthawaffectrootdist;

/// Cali k_max
extern double surf_li_scalar;
extern double soil_om_scalar;
extern double moist_mod_anaerobic;

///////////////////////////////////////////////////////////////////////////////////////
// The Paramlist class (and Paramtype)
//

/// Represents one custom "param" item
/** \see Paramlist */
struct Paramtype {
	xtring name;
	xtring str;
	double num;
};

/// List for the "custom" parameters
/** Functionality for storing and retrieving custom "param" items from the instruction
 *  script. "Custom" parameters can be accessed by other modules without the need to
 *  define them beforehand. This of course also means there is no help text associated
 *  with these parameters, so the user can't get any documentation about them from
 *  the command line.
 *
 * Custom keywords may be included in the instruction script using syntax similar to
 * the following examples:
 *
 * \code
 *     param "co2" (num 340)
 *     param "file_gridlist" (str "gridlist.txt")
 * \endcode
 *
 * To retrieve the values associated with the "param" strings in the above examples,
 * use the following function calls (may appear anywhere in this file; instruction
 * script must have been read in first):
 *
 * \code
 *     param["co2"].num
 *     param["file_gridlist"].str
 * \endcode
 *
 * Each "param" item can store EITHER a number (int or double) OR a string, but not
 * both types of data. Function fail is called to terminate output if a "param" item
 * with the specified identifier was not read in.
 */
class Paramlist : public ListArray<Paramtype> {

public:
	/// Adds a parameter with a numeric value, overwriting if it already existed
	void addparam(xtring name,xtring value);

	/// Adds a parameter with a string value, overwriting if it already existed
	void addparam(xtring name,double value);

	/// Fetches a parameter from the list, aborts the program if it didn't exist
	Paramtype& operator[](xtring name);

	/// Tests if param exists
	bool isparam(xtring name);

	/// Tries to find the parameter in the list
	/** \returns 0 if it wasn't there. */
	Paramtype* find(xtring name);
};

/// The global Paramlist object
/** Contains all the custom parameters after reading in the instruction file */
extern Paramlist param;

/// Reads in the instruction file
/** Uses PLIB library functions to read instructions from file specified by
 * 'insfilename'.
 */
void read_instruction_file(const char* insfilename);

/// Displays documentation about the instruction file parameters to the user
void printhelp();


///////////////////////////////////////////////////////////////////////////////////////
// Interface for declaring parameters from other modules

/// Declares an xtring parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param maxlen   Maximum allowed length of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, xtring* param, int maxlen, const char* help = "");

/// Declares a std:string parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param maxlen   Maximum allowed length of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, std::string* param, int maxlen, const char* help = "");

/// Declares an int parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param min      Minimum allowed value of the parameter in the ins file
 *  \param max      Maximum allowed value of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, int* param, int min, int max, const char* help = "");

/// Declares a double parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param min      Minimum allowed value of the parameter in the ins file
 *  \param max      Maximum allowed value of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, double* param, double min, double max, const char* help = "");

/// Declares a bool parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, bool* param, const char* help = "");

/// Clears parameters declared by input/output modules.
/** This is needed before constructing modules for a new in-process framework()
 *  run, because module instances own some of the variables referenced by these
 *  declarations.
 */
void clear_parameter_declarations();

#endif // LPJ_GUESS_PARAMETERS_H
