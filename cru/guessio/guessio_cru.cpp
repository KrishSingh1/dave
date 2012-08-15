///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio_cru.cpp
/// \brief LPJ-GUESS input/output module with input from instruction script
///
/// This I/O module reads in CRU climate data in a customised binary format.
/// The binary files contain CRU half-degree global historical climate data
/// for 1901-2006.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

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

#ifdef USE_CRU_IO

#include "guessio.h"

#include "driver.h"
#include "outputchannel.h"
#include <plib.h>
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>
#include "globalco2file.h"


// guess2008 - header file for the CRU TS 3.0 data archives
#include "cru_1901_2006.h"
#include "cru_1901_2006misc.h"

// header file for reading binary data archive of global nitrogen deposition
#include "GlobalNitrogenDeposition.h"

//AA CMIP5
//#define CM5rcp85
//#if defined CM5rcp85
//	#include "ipsl_cm5a_lr_historical_1850_2005_r1i1p1.h"
//	#include "ipsl_cm5a_lr_rcp85_2006_2100_r1i1p1.h"
//#endif
#include "cmip5_hist.h"
#include "cmip5_scen.h"
#include "GlobalNitrogenDepositionRCP26.h"
#include "GlobalNitrogenDepositionRCP45.h"
#include "GlobalNitrogenDepositionRCP60.h"
#include "GlobalNitrogenDepositionRCP85.h"

///////////////////////////////////////////////////////////////////////////////////////
//
//                      SECTION: INPUT FROM INSTRUCTION SCRIPT
//
//  - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY -
//
// The first section of this module is concerned with reading simulation settings and
// PFT parameters from an instruction script using functionality from the PLIB library.
// In general model users should not modify this section of the input/output module.
// New instructions (PLIB keywords) may be added (this would require addition of a
// declareitem call in function plib_declarations, and possibly some additional code in
// function plib_callback). However, it is probably preferable to use the "param"
// keyword feature, as this does not require any changes to this section of the module.
//
// Custom keywords may be included in the instruction script using syntax similar to
// the following examples:
//
//   param "co2" (num 340)
//   param "file_gridlist" (str "gridlist.txt")
//
// To retrieve the values associated with the "param" strings in the above examples,
// use the following function calls (may appear anywhere in this file; instruction
// script must have been read in first):
//
//   param["co2"].num
//   param["file_gridlist"].str
//
// Each "param" item can store EITHER a number (int or double) OR a string, but not
// both types of data. Function fail is called to terminate output if a "param" item
// with the specified identifier was not read in.
//
///////////////////////////////////////////////////////////////////////////////////////

/// Represents one custom "param" item
struct Paramtype {
	xtring name;
	xtring str;
	double num;
};

/// List for the custom parameters
/** Functionality for storing and retrieving custom "param" items from the instruction
 *  script.
 */
class Paramlist : public ListArray<Paramtype> {

public:
	/// Adds a parameter with a numeric value, overwriting if it already existed
	void addparam(xtring name,xtring value) {
		Paramtype* p = find(name);
		if (p == 0) {
			p = &createobj();
		}
		p->name=name.lower();
		p->str=value;
	}

	/// Adds a parameter with a string value, overwriting if it already existed
	void addparam(xtring name,double value) {
		Paramtype* p = find(name);
		if (p == 0) {
			p = &createobj();
		}
		p->name=name.lower();
		p->num=value;
	}

	/// Fetches a parameter from the list, aborts the program if it didn't exist
	Paramtype& operator[](xtring name) {
		Paramtype* param = find(name);
		
		if (param == 0) {
			fail("Paramlist::operator[]: parameter \"%s\" not found",(char*)name);
		}
		
		return *param;
	}

private:
	/// Tries to find the parameter in the list
	/** \returns 0 if it wasn't there. */
	Paramtype* find(xtring name) {
		name = name.lower();
		firstobj();
		while (isobj) {
			Paramtype& p=getobj();
			if (p.name==name) return &p;
			nextobj();
		}
		// nothing found
		return 0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// ENUM DECLARATIONS OF INTEGER CONSTANTS FOR PLIB INTERFACE

enum {BLOCK_GLOBAL,BLOCK_PFT,BLOCK_PARAM};
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_LANDCOVER,CB_PHENOLOGY,CB_LEAF,CB_PATHWAY,	
	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM,CB_WATERUPTAKE};


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH FILE SCOPE

Paramlist param;

xtring title; // Title for this run
// guess2008 - new optional parameter
int searchradius; // search radius to use when finding CRU data

/// Landcover fractions read from ins-file (% area).
int lc_fixed_frac[NLANDCOVERTYPES]={0};

/// Whether gridcell is divided into equal active landcover fractions.
bool equal_landcover_area;

Pftlist* ppftlist; // pointer to PFT list
Pft* ppft; // pointer to Pft object currently being assigned to

xtring paramname;
xtring strparam;
double numparam;
bool ifhelp=false;
bool includepft;


// guess2008 - Now declare the output file xtrings here
// Output file names ...
xtring outputdirectory;
xtring file_cmass,file_anpp,file_dens,file_lai,file_cflux,file_cpool,file_runoff;
xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
xtring file_firert,file_speciesheights;
// bvoc
xtring file_aiso,file_miso,file_amon,file_mmon;

xtring file_dgpp,file_dlai,file_dleafN,file_dresp;

// GUESSN
xtring file_cton,file_nmass,file_nsources,file_npool,file_nleach,file_nuptake,file_anppn,file_vmaxnlim,file_nlim,file_nflux;
// end GUESSN

// GUESSN allometry
xtring file_allometry,file_canopyh,file_allometry_ind;
// end GUESSN

void initsettings() {

	// Initialises global settings
	// Parameters not initialised here must be set in instruction script

	iffire=true;
	ifcalcsla=true;
	ifdisturb=false;
	ifcalcsla=false;
	ifcdebt=false;
	distinterval=1.0e10;
	npatch=1;
	vegmode=COHORT;
	searchradius = 0;
	run_landcover = false;

	// guess2008 - initialise filenames here
	outputdirectory = "";
	file_cmass=file_anpp=file_lai=file_cflux=file_dens=file_runoff="";
	file_mnpp=file_mlai=file_maet=file_mpet=file_mevap=file_mrunoff=file_mintercep=file_mrh="";
	file_mgpp=file_mra=file_mnee=file_mwcont_upper=file_mwcont_lower="";
	file_cpool=file_firert=file_speciesheights="";
	// bvoc
	file_aiso=file_miso=file_amon=file_mmon="";

	file_dgpp="";

	// GUESSN
	file_cton=file_nmass=file_nsources=file_npool=file_nleach=file_nuptake=file_anppn=file_vmaxnlim=file_nlim=file_nflux="";
	// end GUESSN

	// GUESSN allometry
	file_allometry=file_canopyh=file_allometry_ind="";
	// end GUESSN
}

void initpft(Pft& pft,xtring& setname) {

	// Initialises a PFT object
	// Parameters not initialised here must be set in instruction script

	pft.name=setname;
	pft.lifeform=NOLIFEFORM;
	pft.phenology=NOPHENOLOGY;

	// Set bioclimatic limits so that PFT can establish and survive under all
	// conditions (may be overridden by settings in instruction script)

	pft.tcmin_surv=-1000.0;
	pft.tcmin_est=-1000.0;
	pft.tcmax_est=1000.0;
	pft.twmin_est=-1000.0;
	pft.gdd5min_est=1000.0;
	pft.twminusc=0.0;

	// Set chilling parameters so that no chilling period required for budburst

	pft.k_chilla=0.0;
	pft.k_chillb=0.0;
	pft.k_chillk=0.0;
}


///////////////////////////////////////////////////////////////////////////////////////
// INPUT FROM INSTRUCTION SCRIPT FILE
// The following code uses functionality from the PLIB library to process an
// instruction script (ins) file containing simulation settings and PFT parameters.
// Function readins() is called by the framework to initiate parsing of the script.
// Function printhelp() is called if GUESS is run with '-help' instead of an ins file
// name as a command line argument. Functions plib_declarations, plib_callback and
// plib_receivemessage comprise part of the interface to PLIB.

void plib_declarations(int id,xtring setname) {

	switch (id) {

	case BLOCK_GLOBAL:


		declareitem("title",&title,80,CB_NONE,"Title for run");
		declareitem("nyear_spinup",&nyear_spinup,1,10000,1,CB_NONE,"Number of simulation years to spinup for");
		declareitem("vegmode",&strparam,16,CB_VEGMODE,
			"Vegetation mode (\"INDIVIDUAL\", \"COHORT\", \"POPULATION\")");
		declareitem("ifdailynpp",&ifdailynpp,1,CB_NONE,
			"Whether photosynthesis calculated daily (alt monthly)");
		declareitem("ifdailydecomp",&ifdailydecomp,1,CB_NONE,
			"Whether soil decomposition calculated daily (alt monthly)");
		declareitem("ifbgestab",&ifbgestab,1,CB_NONE,
			"Whether background establishment enabled (0,1)");
		declareitem("ifsme",&ifsme,1,CB_NONE,
			"Whether spatial mass effect enabled for establishment (0,1)");
		declareitem("ifstochmort",&ifstochmort,1,CB_NONE,
			"Whether mortality stochastic (0,1)");
		declareitem("ifstochestab",&ifstochestab,1,CB_NONE,
			"Whether establishment stochastic (0,1)");
		declareitem("estinterval",&estinterval,1,10,1,CB_NONE,
			"Interval for establishment of new cohorts (years)");
		declareitem("distinterval",&distinterval,1.0,1.0e10,1,CB_NONE,
			"Generic patch-destroying disturbance interval (years)");
		declareitem("iffire",&iffire,1,CB_NONE,
			"Whether fire enabled (0,1)");
		declareitem("ifdisturb",&ifdisturb,1,CB_NONE,
			"Whether generic patch-destroying disturbance enabled (0,1)");
		declareitem("ifcalcsla",&ifcalcsla,1,CB_NONE,
			"Whether SLA calculated from leaf longevity");
		declareitem("ifcdebt",&ifcdebt,1,CB_NONE,
			"Whether to allow C storage");
		declareitem("npatch",&npatch,1,1000,1,CB_NONE,
			"Number of patches simulated");
		declareitem("patcharea",&patcharea,1.0,1.0e4,1,CB_NONE,
			"Patch area (m2)");
		declareitem("wateruptake", &strparam, 20, CB_WATERUPTAKE, 
			"Water uptake mode (\"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");

		// GUESSN
		declareitem("nrelocfrac",&nrelocfrac,0.0,1.0,1,CB_NONE,
			"Fractional N relocation from shed leaves & roots");
		declareitem("nfix_a",&nfix_a,0.0,0.4,1,CB_NONE,
			"first term in N fixation eqn");
		declareitem("nfix_b",&nfix_b,-10.0,10.,1,CB_NONE,
			"second term in N fixation eqn");

		declareitem("ifcentury",&ifcentury,1,CB_NONE,
			"Whether to use CENTURY SOM dynamics (default standard LPJ)");
		declareitem("ifnlim",&ifnlim,1,CB_NONE,
			"Whether plant growth limited by available N");
		declareitem("freenyears",&freenyears,0,1000,1,CB_NONE,
			"Number of years to spinup without N limitation");
		declareitem("ifleachn",&ifleachn,1,CB_NONE,
			"Whether to allow N leaching");
		declareitem("ifindiv_fnuptake",&ifindiv_fnuptake,1,CB_NONE,
			"Whether to allow individual fractional N uptake");
		// end GUESSN

		// guess2008
		// Annual output variables
		declareitem("outputdirectory",&outputdirectory,300,CB_NONE,"Directory for the output files");
		declareitem("file_cmass",&file_cmass,300,CB_NONE,"C biomass output file");
		declareitem("file_anpp",&file_anpp,300,CB_NONE,"Annual NPP output file");
		declareitem("file_lai",&file_lai,300,CB_NONE,"LAI output file");
		declareitem("file_cflux",&file_cflux,300,CB_NONE,"C fluxes output file");
		declareitem("file_dens",&file_dens,300,CB_NONE,"Tree density output file");
		declareitem("file_cpool",&file_cpool,300,CB_NONE,"Soil C output file");
		declareitem("file_runoff",&file_runoff,300,CB_NONE,"Runoff output file");
		declareitem("file_firert",&file_firert,300,CB_NONE,"Fire retrun time output file");
		
		// GUESSN
		declareitem("file_cton",&file_cton,300,CB_NONE,"Mean leaf C:N output file");
		declareitem("file_nmass",&file_nmass,300,CB_NONE,"N biomass output file");
		declareitem("file_nsources",&file_nsources,300,CB_NONE,"annual N sources output file");
		declareitem("file_npool",&file_npool,300,CB_NONE,"Soil N output file");
		declareitem("file_nleach",&file_nleach,300,CB_NONE,"Leached mineral N output file");
		declareitem("file_nuptake",&file_nuptake,300,CB_NONE,"annual N uptake output file");
		declareitem("file_anppn",&file_anppn,300,CB_NONE,"annual N usage output file");
		declareitem("file_vmaxnlim",&file_vmaxnlim,300,CB_NONE,"annual N limitation on vm output file");
		declareitem("file_nlim",&file_nlim,300,CB_NONE,"annual N limitation on growth output file");
		declareitem("file_nflux",&file_nflux,300,CB_NONE,"annual N fluxes output file");
		// end GUESSN

		// GUESSN allometry
		declareitem("file_allometry",&file_allometry,300,CB_NONE,"Allometry output file");
		declareitem("file_canopyh",&file_canopyh,300,CB_NONE,"Canopy height output file");
		declareitem("file_allometry_ind",&file_allometry_ind,300,CB_NONE,"Individual Allometry output file");
		// end GUESSN

		
		declareitem("file_speciesheights",&file_speciesheights,300,CB_NONE,"Mean species heights");

		// Monthly output variables
		declareitem("file_mnpp",&file_mnpp,300,CB_NONE,"Monthly NPP output file");
		declareitem("file_mlai",&file_mlai,300,CB_NONE,"Monthly LAI output file");
		declareitem("file_mgpp",&file_mgpp,300,CB_NONE,"Monthly GPP-LeafResp output file");
		declareitem("file_mra",&file_mra,300,CB_NONE,"Monthly autotrophic respiration output file");
		declareitem("file_maet",&file_maet,300,CB_NONE,"Monthly AET output file");
		declareitem("file_mpet",&file_mpet,300,CB_NONE,"Monthly PET output file");
		declareitem("file_mevap",&file_mevap,300,CB_NONE,"Monthly Evap output file");
		declareitem("file_mrunoff",&file_mrunoff,300,CB_NONE,"Monthly runoff output file");
		declareitem("file_mintercep",&file_mintercep,300,CB_NONE,"Monthly intercep output file");
		declareitem("file_mrh",&file_mrh,300,CB_NONE,"Monthly heterotrphic respiration output file");
		declareitem("file_mnee",&file_mnee,300,CB_NONE,"Monthly NEE output file");
		declareitem("file_mwcont_upper",&file_mwcont_upper,300,CB_NONE,"Monthly wcont_upper output file");
		declareitem("file_mwcont_lower",&file_mwcont_lower,300,CB_NONE,"Monthly wcont_lower output file");
		// bvoc
		declareitem("file_aiso",&file_aiso,300,CB_NONE,"annual isoprene flux output file");
		declareitem("file_miso",&file_miso,300,CB_NONE,"monthly isoprene flux output file");
		declareitem("file_amon",&file_amon,300,CB_NONE,"annual monoterpene flux output file");
		declareitem("file_mmon",&file_mmon,300,CB_NONE,"monthly monoterpene flux output file");

		declareitem("file_dgpp",&file_dgpp,300,CB_NONE,"Daily GPP output file");
		declareitem("file_dlai",&file_dlai,300,CB_NONE,"Daily LAI output file");
		declareitem("file_dleafN",&file_dleafN,300,CB_NONE,"Daily Leaf N output file");
		declareitem("file_dresp",&file_dresp,300,CB_NONE,"Daily Resp output file");

		// guess2008 - new options
		declareitem("ifsmoothgreffmort",&ifsmoothgreffmort,1,CB_NONE,
			"Whether to vary mort_greff smoothly with growth efficiency (0,1)");
		declareitem("ifdroughtlimitedestab",&ifdroughtlimitedestab,1,CB_NONE,
			"Whether establishment drought limited (0,1)");
		declareitem("ifrainonwetdaysonly",&ifrainonwetdaysonly,1,CB_NONE,
			"Whether it rains on wet days only (1), or a little every day (0);");
		declareitem("searchradius", &searchradius, 0, 100, 1, CB_NONE,
			"If specified, CRU data will be searched for in a circle");

		// bvoc 
		declareitem("ifbvoc",&ifbvoc,1,CB_NONE,
			"Whether or not BVOC calculations are performed (0,1)");
		declareitem("run_landcover",&run_landcover,1,CB_NONE,"Landcover version");
		declareitem("run_urban",&run[URBAN],1,CB_NONE,"Whether urban land is to be simulated");
		declareitem("run_crop",&run[CROPLAND],1,CB_NONE,"Whether crop-land is to be simulated");
		declareitem("run_pasture",&run[PASTURE],1,CB_NONE,"Whether pasture is to be simulated");
		declareitem("run_forest",&run[FOREST],1,CB_NONE,"Whether managed forest is to be simulated");
		declareitem("run_natural",&run[NATURAL],1,CB_NONE,"Whether natural vegetation is to be simulated");
		declareitem("run_peatland",&run[PEATLAND],1,CB_NONE,"Whether peatland is to be simulated");
		declareitem("ifslowharvestpool",&ifslowharvestpool,1,CB_NONE,"If a slow harvested product pool is included in patchpft.");
		declareitem("lcfrac_fixed",&lcfrac_fixed,1,CB_NONE,"Whether static landcover fractions are set in the ins-file (0,1)");
		declareitem("equal_landcover_area",&equal_landcover_area,1,CB_NONE,"Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
		declareitem("lc_fixed_urban",&lc_fixed_frac[URBAN],0,100,1,CB_NONE,"% lc_fixed_urban");
		declareitem("lc_fixed_cropland",&lc_fixed_frac[CROPLAND],0,100,1,CB_NONE,"% lc_fixed_cropland");
		declareitem("lc_fixed_pasture",&lc_fixed_frac[PASTURE],0,100,1,CB_NONE,"% lc_fixed_pasture");
		declareitem("lc_fixed_forest",&lc_fixed_frac[FOREST],0,100,1,CB_NONE,"% lc_fixed_forest");
		declareitem("lc_fixed_natural",&lc_fixed_frac[NATURAL],0,100,1,CB_NONE,"% lc_fixed_natural");
		declareitem("lc_fixed_peatland",&lc_fixed_frac[PEATLAND],0,100,1,CB_NONE,"% lc_fixed_peatland");

		// CMIP5 - land use input
		declareitem("ifcmip5",&ifcmip5,1,CB_NONE,
			    "Whether CMIP5 climate should be used");
		declareitem("iflandusesimple",&iflandusesimple,1,CB_NONE,
			    "Whether to apply a simple land use representation (0,1)");
		declareitem("iflandusechange",&iflandusechange,1,CB_NONE,
			    "Whether land use is static (0) or dynamic (1)");

		declareitem("pft",BLOCK_PFT,CB_NONE,"Header for block defining PFT");
		declareitem("param",BLOCK_PARAM,CB_NONE,"Header for custom parameter block");
		callwhendone(CB_CHECKGLOBAL);

		break;
	
	case BLOCK_PFT:

		if (!ifhelp) {

			// Create and initialise a new Pft object and obtain a reference to it
			
			ppft=&ppftlist->createobj();
			initpft(*ppft,setname);
			includepft=true;
		}

		declareitem("include",&includepft,1,CB_NONE,"Include PFT in analysis");
		declareitem("lifeform",&strparam,16,CB_LIFEFORM,
			"Lifeform (\"TREE\" or \"GRASS\")");
		declareitem("landcover",&strparam,16,CB_LANDCOVER,
			"Landcovertype (\"URBAN\", \"CROP\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
		declareitem("phenology",&strparam,16,CB_PHENOLOGY,
			"Phenology (\"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
		declareitem("phengdd5ramp",&ppft->phengdd5ramp,0.0,1000.0,1,CB_NONE,
			"GDD on 5 deg C base to attain full leaf cover");
		declareitem("wscal_min",&ppft->wscal_min,0.0,1.0,1,CB_NONE,
			"Water stress threshold for leaf abscission (raingreen PFTs)");
		declareitem("pathway",&strparam,16,CB_PATHWAY,
			"Biochemical pathway (\"C3\" or \"C4\")");
		declareitem("pstemp_min",&ppft->pstemp_min,-50.0,50.0,1,CB_NONE,
			"Approximate low temp limit for photosynthesis (deg C)");
		declareitem("pstemp_low",&ppft->pstemp_low,-50.0,50.0,1,CB_NONE,
			"Approx lower range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_high",&ppft->pstemp_high,0.0,60.0,1,CB_NONE,
			"Approx higher range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_max",&ppft->pstemp_max,0.0,60.0,1,CB_NONE,
			"Maximum temperature limit for photosynthesis (deg C)");
		declareitem("lambda_max",&ppft->lambda_max,0.1,0.99,1,CB_NONE,
			"Non-water-stressed ratio of intercellular to ambient CO2 pp");
		declareitem("rootdist",ppft->rootdist,0.0,1.0,NSOILLAYER,CB_ROOTDIST,
			"Fraction of roots in each soil layer (first value=upper layer)");
		declareitem("gmin",&ppft->gmin,0.0,1.0,1,CB_NONE,
			"Canopy conductance not assoc with photosynthesis (mm/s)");
		declareitem("emax",&ppft->emax,0.0,50.0,1,CB_NONE,
			"Maximum evapotranspiration rate (mm/day)");
		// guess2008 - increased the upper limit to possible respcoeff values (was 1.2)
		declareitem("respcoeff",&ppft->respcoeff,0.0,3,1,CB_NONE,
			"Respiration coefficient (0-1)");
		declareitem("leaftype",&strparam,16,CB_LEAF,
			"Leaf type (\"BROAD\" or \"NEEDLE\")");
		

		// GUESSN
		declareitem("cton_leaf_min",&ppft->cton_leaf_min,1.0,1.0e4,1,CB_NONE,
			"Min Leaf C:N mass ratio");
		declareitem("cton_leaf_max",&ppft->cton_leaf_max,1.0,1.0e4,1,CB_NONE,
			"Max Leaf C:N mass ratio");
		declareitem("cton_leaf_avr",&ppft->cton_leaf_avr,1.0,1.0e4,1,CB_NONE,
			"Average Leaf C:N mass ratio");
		declareitem("cton_root_avr",&ppft->cton_root_avr,1.0,1.0e4,1,CB_NONE,
			"Average Fine root C:N mass ratio");
		declareitem("cton_sap_avr",&ppft->cton_sap_avr,1.0,1.0e4,1,CB_NONE,
			"Average Sapwood C:N mass ratio");
		declareitem("n_reserve",&ppft->n_reserve,0.0,1.0,1,CB_NONE,
			"N storage organ in relation to sapwood carbon");
		declareitem("a0",&ppft->a0,0.0,2.0,1,CB_NONE,
			"Intercept parameter in the relation between leaf N not associated with photosynthesis and total leaf N");
		// end GUESSN

		declareitem("reprfrac",&ppft->reprfrac,0.0,1.0,1,CB_NONE,
			"Fraction of NPP allocated to reproduction");
		declareitem("reprCN",&ppft->reprCN,1.0,1.0e4,1,CB_NONE,
			"Reproduction CN ratio");
		declareitem("turnover_leaf",&ppft->turnover_leaf,0.0,1.0,1,CB_NONE,
			"Leaf turnover (fraction/year)");
		declareitem("turnover_root",&ppft->turnover_root,0.0,1.0,1,CB_NONE,
			"Fine root turnover (fraction/year)");
		declareitem("turnover_sap",&ppft->turnover_sap,0.0,1.0,1,CB_NONE,
			"Sapwood turnover (fraction/year)");
		declareitem("wooddens",&ppft->wooddens,10.0,1000.0,1,CB_NONE,
			"Sapwood and heartwood density (kgC/m3)");
		declareitem("crownarea_max",&ppft->crownarea_max,1.0,1000.0,1,CB_NONE,
			"Maximum tree crown area (m2)");
		declareitem("k_allom1",&ppft->k_allom1,10.0,1000.0,1,CB_NONE,
			"Constant in allometry equations");
		// guess2008 - changed lower limit for k_allom2 to 1 from 10. This is needed
		// for the shrub allometries.
		declareitem("k_allom2",&ppft->k_allom2,1.0,1.0e4,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_allom3",&ppft->k_allom3,0.1,1.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_rp",&ppft->k_rp,1.0,2.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_latosa",&ppft->k_latosa,100.0,1.0e5,1,CB_NONE,
			"Tree leaf to sapwood xs area ratio");
		declareitem("sla",&ppft->sla,1.0,1000.0,1,CB_NONE,
			"Specific leaf area (m2/kgC)");
		declareitem("ltor_max",&ppft->ltor_max,0.1,10.0,1,CB_NONE,
			"Non-water-stressed leaf:fine root mass ratio");
		declareitem("litterme",&ppft->litterme,0.0,1.0,1,CB_NONE,
			"Litter moisture flammability threshold (fraction of AWC)");
		declareitem("fireresist",&ppft->fireresist,0.0,1.0,1,CB_NONE,
			"Fire resistance (0-1)");
		declareitem("tcmin_surv",&ppft->tcmin_surv,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for survival (deg C)");
		declareitem("tcmin_est",&ppft->tcmin_est,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for establishment (deg C)");
		declareitem("tcmax_est",&ppft->tcmax_est,-50.0,1000.0,1,CB_NONE,
			"Max 20-year coldest month mean temp for establishment (deg C)");
		declareitem("twmin_est",&ppft->twmin_est,-1000.0,50.0,1,CB_NONE,
			"Min warmest month mean temp for establishment (deg C)");
		declareitem("twminusc",&ppft->twminusc,0,100,1,CB_NONE,
			"Stupid larch parameter");
		declareitem("gdd5min_est",&ppft->gdd5min_est,0.0,5000.0,1,CB_NONE,
			"Min GDD on 5 deg C base for establishment");
		declareitem("k_chilla",&ppft->k_chilla,0.0,5000.0,1,CB_NONE,
			"Constant in equation for budburst chilling time requirement");
		declareitem("k_chillb",&ppft->k_chillb,0.0,5000.0,1,CB_NONE,
			"Coefficient in equation for budburst chilling time requirement");
		declareitem("k_chillk",&ppft->k_chillk,0.0,1.0,1,CB_NONE,
			"Exponent in equation for budburst chilling time requirement");
		declareitem("parff_min",&ppft->parff_min,0.0,1.0e7,1,CB_NONE,
			"Min forest floor PAR for grass growth/tree estab (J/m2/day)");
		declareitem("alphar",&ppft->alphar,0.01,100.0,1,CB_NONE,
			"Shape parameter for recruitment-juv growth rate relationship");
		declareitem("est_max",&ppft->est_max,1.0e-4,1.0,1,CB_NONE,
			"Max sapling establishment rate (indiv/m2/year)");
		declareitem("kest_repr",&ppft->kest_repr,1.0,1000.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_bg",&ppft->kest_bg,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_pres",&ppft->kest_pres,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("longevity",&ppft->longevity,0.0,3000.0,1,CB_NONE,
			"Expected longevity under lifetime non-stressed conditions (yr)");
		declareitem("greff_min",&ppft->greff_min,0.0,1.0,1,CB_NONE,
			"Threshold for growth suppression mortality (kgC/m2 leaf/yr)");
		declareitem("leaflong",&ppft->leaflong,0.1,100.0,1,CB_NONE,
			"Leaf longevity (years)");
		declareitem("intc",&ppft->intc,0.0,1.0,1,CB_NONE,"Interception coefficient");
		
		// guess2008 - DLE
		declareitem("drought_tolerance",&ppft->drought_tolerance,0.0,1.0,1,CB_NONE,
			"Drought tolerance level (0 = very -> 1 = not at all) (unitless)");

		// bvoc
		declareitem("ga",&ppft->ga,0.0,1.0,1,CB_NONE,
			"aerodynamic conductance (m/s)");
		declareitem("eps_iso",&ppft->eps_iso,0.,100.,1,CB_NONE,
			"isoprene emission capacity (ug C g-1 h-1)");
		declareitem("seas_iso",&ppft->seas_iso,1,CB_NONE,
			"whether (1) or not (0) isoprene emissions show seasonality");
		declareitem("eps_mon",&ppft->eps_mon,0.,100.,1,CB_NONE,
			"monoterpene emission capacity (ug C g-1 h-1)");
		declareitem("storfrac_mon",&ppft->storfrac_mon,0.,1.,1,CB_NONE,
			"fraction of monoterpene production that goes into storage pool (-)");
		
		declareitem("harv_eff",&ppft->harv_eff,0.0,1.0,1,CB_NONE,"Harvest efficiency");
		declareitem("harvest_slow_frac",&ppft->harvest_slow_frac,0.0,1.0,1,CB_NONE,
			"Fraction of harvested products that goes into carbon depository for long-lived products like wood");
		declareitem("turnover_harv_prod",&ppft->turnover_harv_prod,0.0,1.0,1,CB_NONE,"Harvested products turnover (fraction/year)");
		declareitem("res_outtake",&ppft->res_outtake,0.0,1.0,1,CB_NONE,"Fraction of residue outtake at harvest");

		callwhendone(CB_CHECKPFT);
		
		break;

	case BLOCK_PARAM:

		paramname=setname;
		declareitem("str",&strparam,300,CB_STRPARAM,
			"String value for custom parameter");
		declareitem("num",&numparam,-1.0e38,1.0e38,1,CB_NUMPARAM,
			"Numerical value for custom parameter");
		
		break;
	}
}

void badins(xtring missing) {

	xtring message=(xtring)"Missing mandatory setting: "+missing;
	sendmessage("Error",message);
	plibabort();
}

void plib_callback(int callback) {

	xtring message;
	int i;
	double numval;

	switch (callback) {

	case CB_VEGMODE:
		if (strparam.upper()=="INDIVIDUAL") vegmode=INDIVIDUAL;
		else if (strparam.upper()=="COHORT") vegmode=COHORT;
		else if (strparam.upper()=="POPULATION") vegmode=POPULATION;
		else {
			sendmessage("Error",
				"Unknown vegetation mode (valid types: \"INDIVIDUAL\",\"COHORT\", \"POPULATION\")");
			plibabort();
		}
		break;
	case CB_WATERUPTAKE:
		if (strparam.upper() == "WCONT") wateruptake = WR_WCONT;
		else if (strparam.upper() == "ROOTDIST") wateruptake = WR_ROOTDIST;
		else if (strparam.upper() == "SMART") wateruptake = WR_SMART;
		else if (strparam.upper() == "SPECIESSPECIFIC") wateruptake = WR_SPECIESSPECIFIC;
		else {
			sendmessage("Error",
				"Unknown water uptake mode (valid types: \"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");
		}
		break;
	case CB_LIFEFORM:
		if (strparam.upper()=="TREE") ppft->lifeform=TREE;
		else if (strparam.upper()=="GRASS") ppft->lifeform=GRASS;
		else {
			sendmessage("Error",
				"Unknown lifeform type (valid types: \"TREE\", \"GRASS\")");
			plibabort();
		}
		break;
	case CB_LANDCOVER:
		if (strparam.upper()=="NATURAL") ppft->landcover=NATURAL;
		else if (strparam.upper()=="URBAN") ppft->landcover=URBAN;
		else if (strparam.upper()=="CROPLAND") ppft->landcover=CROPLAND;
		else if (strparam.upper()=="PASTURE") ppft->landcover=PASTURE;
		else if (strparam.upper()=="FOREST") ppft->landcover=FOREST;			
		else if (strparam.upper()=="PEATLAND") ppft->landcover=PEATLAND;
		else {
			sendmessage("Error",
				"Unknown landcover type (valid types: \"URBAN\", \"CROPLAND\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
			plibabort();
		}
		break;
	case CB_PHENOLOGY:
		if (strparam.upper()=="SUMMERGREEN") ppft->phenology=SUMMERGREEN;
		else if (strparam.upper()=="RAINGREEN") ppft->phenology=RAINGREEN;
		else if (strparam.upper()=="EVERGREEN") ppft->phenology=EVERGREEN;
		else if (strparam.upper()=="ANY") ppft->phenology=ANY;
		else {
			sendmessage("Error",
				"Unknown phenology type\n  (valid types: \"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
			plibabort();
		}
		break;
	case CB_LEAF:
		if (strparam.upper()=="BROAD") ppft->leaf=BROAD;
		else if (strparam.upper()=="NEEDLE") ppft->leaf=NEEDLE;
		else {
			sendmessage("Error",
				"Unknown leaf type\n  (valid types: \"BROAD\" or \"NEEDLE\")");
			plibabort();
		}
		break;
	case CB_PATHWAY:
		if (strparam.upper()=="C3") ppft->pathway=C3;
		else if (strparam.upper()=="C4") ppft->pathway=C4;
		else {
			sendmessage("Error",
				"Unknown pathway type\n  (valid types: \"C3\" or \"C4\")");
			plibabort();
		}
		break;
	case CB_ROOTDIST:
		numval=0.0;
		for (i=0;i<NSOILLAYER;i++) numval+=ppft->rootdist[i];
		if (numval<0.99 || numval>1.01) {
			sendmessage("Error","Specified root fractions do not sum to 1.0");
			plibabort();
		}
		ppft->rootdist[NSOILLAYER-1]+=1.0-numval;
		break;
	case CB_STRPARAM:
		param.addparam(paramname,strparam);
		break;
	case CB_NUMPARAM:
		param.addparam(paramname,numparam);
		break;
	case CB_CHECKGLOBAL:
		if (!itemparsed("title")) badins("title");
		if (!itemparsed("nyear_spinup")) badins("nyear_spinup");
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("ifdailynpp")) badins("ifdailynpp");
		if (!itemparsed("ifdailydecomp")) badins("ifdailydecomp");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");
		if (!itemparsed("wateruptake")) badins("wateruptake");

		// GUESSN
		if (!itemparsed("nrelocfrac")) badins("nrelocfrac");
		if (!itemparsed("nfix_a")) badins("nfix_a");
		if (!itemparsed("nfix_b")) badins("nfix_b");

		if (!itemparsed("ifcentury")) badins("ifcentury");
		if (!itemparsed("ifnlim")) badins("ifnlim");
		if (!itemparsed("freenyears")) badins("freenyears");
		if (!itemparsed("ifleachn")) badins("ifleachn");
		if (!itemparsed("ifindiv_fnuptake")) badins("ifindiv_fnuptake");
		// end GUESSN

		// guess2008
		if (!itemparsed("outputdirectory")) badins("outputdirectory");
		if (!itemparsed("ifsmoothgreffmort")) badins("ifsmoothgreffmort");
		if (!itemparsed("ifdroughtlimitedestab")) badins("ifdroughtlimitedestab");
		if (!itemparsed("ifrainonwetdaysonly")) badins("ifrainonwetdaysonly");
		// bvoc
		if (!itemparsed("ifbvoc")) badins("ifbvoc");

		if (!itemparsed("run_landcover")) badins("run_landcover");
		if (run_landcover) {
			if (!itemparsed("lcfrac_fixed")) badins("lcfrac_fixed");
			if (!itemparsed("equal_landcover_area")) badins("equal_landcover_area");
			if (!itemparsed("lc_fixed_urban")) badins("lc_fixed_urban");
			if (!itemparsed("lc_fixed_cropland")) badins("lc_fixed_cropland");
			if (!itemparsed("lc_fixed_pasture")) badins("lc_fixed_pasture");
			if (!itemparsed("lc_fixed_forest")) badins("lc_fixed_forest");
			if (!itemparsed("lc_fixed_natural")) badins("lc_fixed_natural");
			if (!itemparsed("lc_fixed_peatland")) badins("lc_fixed_peatland");
			if (!itemparsed("run_natural")) badins("run_natural");
			if (!itemparsed("run_crop")) badins("run_crop");
			if (!itemparsed("run_forest")) badins("run_forest");
			if (!itemparsed("run_urban")) badins("run_urban");
			if (!itemparsed("run_pasture")) badins("run_pasture");
			if (!itemparsed("ifslowharvestpool")) badins("ifslowharvestpool");
		}

		// CMIP5 - land use input
		if (!itemparsed("ifcmip5")) badins("ifcmip5");
		if (!itemparsed("iflandusesimple")) badins("iflandusesimple");
		if (!itemparsed("iflandusechange")) badins("iflandusechange");

		if (!itemparsed("pft")) badins("pft");
		if (vegmode==COHORT || vegmode==INDIVIDUAL) {
			if (!itemparsed("ifbgestab")) badins("ifbgestab");
			if (!itemparsed("ifsme")) badins("ifsme");
			if (!itemparsed("ifstochmort")) badins("ifstochmort");
			if (!itemparsed("ifstochestab")) badins("ifstochestab");
			if (itemparsed("ifdisturb") && !itemparsed("distinterval"))
				badins("distinterval");
			if (!itemparsed("npatch")) badins("npatch");
			if (!itemparsed("patcharea")) badins("patcharea");
			if (!itemparsed("estinterval")) badins("estinterval");
		}
		else if (vegmode==POPULATION && npatch!=1) {
			sendmessage("Information",
				"Value specified for npatch ignored in population mode");
			npatch=1;
		}
		break;
	case CB_CHECKPFT:
		if (!itemparsed("lifeform")) badins("lifeform");
		if (!itemparsed("phenology")) badins("phenology");
		if (ppft->phenology==SUMMERGREEN || ppft->phenology==ANY)
			if (!itemparsed("phengdd5ramp")) badins("phengdd5ramp");
		if (ppft->phenology==RAINGREEN || ppft->phenology==ANY)
			if (!itemparsed("wscal_min")) badins("wscal_min");
		if (!itemparsed("pathway")) badins("pathway");
		if (!itemparsed("pstemp_min")) badins("pstemp_min");
		if (!itemparsed("pstemp_low")) badins("pstemp_low");
		if (!itemparsed("pstemp_high")) badins("pstemp_high");
		if (!itemparsed("pstemp_max")) badins("pstemp_max");
		if (!itemparsed("lambda_max")) badins("lambda_max");
		if (!itemparsed("rootdist")) badins("rootdist");
		if (!itemparsed("gmin")) badins("gmin");
		if (!itemparsed("emax")) badins("emax");
		if (!itemparsed("respcoeff")) badins("respcoeff");
		if (!itemparsed("sla") && !ifcalcsla) badins("sla");

		// GUESSN
		if (!itemparsed("cton_leaf_min")) badins("cton_leaf_min");
		if (!itemparsed("cton_leaf_max")) badins("cton_leaf_max");
		if (!itemparsed("cton_leaf_avr")) badins("cton_leaf_avr");
		if (!itemparsed("cton_root_avr")) badins("cton_root_avr");
		if (!itemparsed("n_reserve")) badins("n_reserve");
		if (!itemparsed("a0")) badins("a0");
		// end GUESSN

		if (!itemparsed("reprfrac")) badins("reprfrac");
		if (!itemparsed("reprCN")) badins("reprCN");
		if (!itemparsed("turnover_leaf")) badins("turnover_leaf");
		if (!itemparsed("turnover_root")) badins("turnover_root");
		if (!itemparsed("ltor_max")) badins("ltor_max");
		if (!itemparsed("intc")) badins("intc");

		if (run_landcover)
		{
			if (!itemparsed("landcover")) badins("landcover");
			if (!itemparsed("turnover_harv_prod")) badins("turnover_harv_prod");
			if (!itemparsed("harvest_slow_frac")) badins("harvest_slow_frac");
			if (!itemparsed("harv_eff")) badins("harv_eff");
			if (!itemparsed("res_outtake")) badins("res_outtake");
		}

		// guess2008 - DLE
		if (!itemparsed("drought_tolerance")) badins("drought_tolerance");

		// bvoc
		if(ifbvoc){
		  if (!itemparsed("ga")) badins("ga");
		  if (!itemparsed("eps_iso")) badins("eps_iso");
		  if (!itemparsed("seas_iso")) badins("seas_iso");
		  if (!itemparsed("eps_mon")) badins("eps_mon");
		  if (!itemparsed("storfrac_mon")) badins("storfrac_mon");
		}

		if (ppft->lifeform==TREE) {
			if (!itemparsed("cton_sap_avr")) badins("cton_sap_avr");
			if (!itemparsed("turnover_sap")) badins("turnover_sap");
			if (!itemparsed("wooddens")) badins("wooddens");
			if (!itemparsed("crownarea_max")) badins("crownarea_max");
			if (!itemparsed("k_allom1")) badins("k_allom1");
			if (!itemparsed("k_allom2")) badins("k_allom2");
			if (!itemparsed("k_allom3")) badins("k_allom3");
			if (!itemparsed("k_rp")) badins("k_rp");
			if (!itemparsed("k_latosa")) badins("k_latosa");
			if (vegmode==COHORT || vegmode==INDIVIDUAL) {
				if (!itemparsed("kest_repr")) badins("kest_repr");
				if (!itemparsed("kest_bg")) badins("kest_bg");
				if (!itemparsed("kest_pres")) badins("kest_pres");
				if (!itemparsed("longevity")) badins("longevity");
				if (!itemparsed("greff_min")) badins("greff_min");		
				if (!itemparsed("alphar")) badins("alphar");
				if (!itemparsed("est_max")) badins("est_max");
			}
		}
		if (iffire) {
			if (!itemparsed("litterme")) badins("litterme");
			if (!itemparsed("fireresist")) badins("fireresist");
		}
		if (ifcalcsla) {
			if (!itemparsed("leaflong")) {
				sendmessage("Error",
					"Value required for leaflong when ifcalcsla enabled");
				plibabort();
			}
			if (itemparsed("sla"))
				sendmessage("Warning",
				"Specified sla value not used when ifcalcsla enabled");

			// Calculate SLA
			ppft->initsla();
		}
		if (vegmode==COHORT || vegmode==INDIVIDUAL) {
			if (!itemparsed("parff_min")) badins("parff_min");	
		}

		// Calculate regeneration characteristics for population mode
		ppft->initregen();

		ppft->id=npft++;
			// VERY IMPORTANT (cannot rely on internal id counter of collection class)

		//	delete unused pft:s from pftlist

		if (ppft->landcover!=NATURAL) {
			if (!run_landcover || !run[ppft->landcover])
				includepft=0;
		}
		else if (run_landcover && !run[NATURAL]) {
			if (ppft->landcover==NATURAL)
				includepft=0;
		}

		// If "include 0", remove this PFT from list, and set id to correct value

		if (!includepft) {
			ppftlist->killobj();
			npft--;
		}

		break;
	}
}

void plib_receivemessage(xtring text) {

	// Output of messages to user sent by PLIB

	dprintf((char*)text);
}

bool readins(xtring filename,Pftlist& pftlist) {

	// DESCRIPTION
	// Uses PLIB library functions to read instructions from file specified by
	// 'filename', returning true if file could be successfully opened and read, and
	// no errors were encountered.

	// OUTPUT PARAMETERS
	// pftlist  = initialised list array of PFT parameters

	// Store global pointer to pftlist
	ppftlist=&pftlist;

	// Initialise PFT count
	npft=0;

	// Initialise certain parameters
	initsettings();
	param.killall();

	// Call PLIB
	return plib(filename);
}

void printhelp() {

	// Calls PLIB to output help text

	ifhelp=true;
	plibhelp();
	ifhelp=false;
}


///////////////////////////////////////////////////////////////////////////////////////
//
//             SECTION: INPUT OF ENVIRONMENTAL DRIVING DATA FOR SIMULATION
//                            OUTPUT OF SIMULATION RESULTS
//
// In general it is the responsibility of the user of the model to provide code for
// this section of the input/output module. The following functions are called by the
// framework at various stages of the simulation and should contain appropriate code:
//
// void initio(int argc,char* argv[],Pftlist& pftlist)
//   Initialises input/output (e.g. opening files), sets values for the global
//   simulation parameter variables (currently vegmode, npatch, patcharea, ifdailynpp,
//   ifdailydecomp, ifbgestab, ifsme, ifstochestab, ifstochmort, iffire, estinterval,
//   npft), initialises pftlist (the one and only list of PFTs and their static
//   parameters for this run of the model). Normally all of the above parameters,
//   and possibly others, are read from the ins file (see above). Function readins
//   should be called to input settings from the ins file. The syntax for this call
//   should be similar to the following (note that readins returns false in the event
//   of an error in the ins file; normally this should result in program termination):
//
//   xtring insfilename=argv[1];
//   if (!readins(insfilename,pftlist))
//       fail("\nUsage: %s <instruction-script-filename> | -help",argv[0]);
//
//   Arguments argc and argv normally correspond to the command-line arguments
//   imported from the main function (main module, usually main.cpp). The first
//   command line argument (argv[0]) is the name of the binary executable (e.g.
//   guess, guess.exe); the second (argv[1]) should normally be the ins file name.
//   This demonstration version of initio also implements "-help" as an alternative
//   command-line argument, resulting in output of a brief description of the
//   keywords recognised in the ins file, instead of a model run.
//
// bool getgridcell(Gridcell& gridcell)
//   Obtains coordinates and soil static parameters for the next grid cell to
//   simulate. The function should return false if no grid cells remain to be simulated,
//   otherwise true. Currently the following member variables of gridcell should be
//   initialised: longitude, latitude and climate.instype; the following members of
//   member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
//   thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
//   code (Sitch et al 2000) by a call to function soilparameters in the driver
//   module (driver.cpp):
//
//   soilparameters(gridcell.soiltype,soilcode);
//
//   If the model is to be driven by quasi-daily values of the climate variables
//   derived from monthly means, this function may be the appropriate place to
//   perform the required interpolations. The utility functions interp_monthly_means
//   and interp_monthly_totals in driver.cpp may be called for this purpose.
//
// bool getclimate(Gridcell& gridcell)
//   Obtains climate data (including atmospheric CO2 and insolation) for this day.
//   The function should return false if the simulation is complete for this grid cell,
//   otherwise true. This will normally require querying the year and day member
//   variables of the global class object date:
//
//   if (date.day==0 && date.year==nyear_spinup) return false;
//   // else
//   return true;
//
//   Currently the following member variables of the climate member of gridcell must be
//   initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
//   values of the climate variables derived from monthly means, this day's values
//   will presumably be extracted from arrays containing the interpolated daily
//   values (see function getgridcell):
//
//   gridcell.climate.temp=dtemp[date.day];
//   gridcell.climate.prec=dprec[date.day];
//   gridcell.climate.insol=dsun[date.day];
//
//   Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
//   BVOC:
//   gridcell.climate.dtr=ddtr[date.day]; 
//
//   If model is run in diurnal mode, which requires appropriate climate forcing data, 
//   additional members of the climate must be initialised: temps, insols. Both of the
//   variables must be of type std::vector. The length of these vectors should be equal
//   to value of date.subdaily which also needs to be set either in getclimate or 
//   getgridcell functions. date.subdaily is a number of sub-daily period in a single 
//   day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
//   diurnal mode.
//
// void outannual(Gridcell& gridcell,Pftlist& pftlist)
//   Called at the end of the last day of each simulation year to permit output of
//   model results.
//
// termio()
//   Called after simulation is complete for all gridcells to allow memory deallocation,
//   closing of files or other cleanup functions.
//
///////////////////////////////////////////////////////////////////////////////////////

struct Coord {

	// Type for storing grid cell longitude, latitude and description text

	int id;
	double lon;
	double lat;
	xtring descrip;

};


ListArray_id<Coord> gridlist;
	// Will maintain a list of Coord objectsc ontaining coordinates
	// of the grid cells to simulate

int ngridcell; // the number of grid cells to simulate

class Spinup_data {

	// Class for management of climate data for spinup
	// (derived from first few years of historical climate data)

private:
	int nyear;
	int thisyear;
	double* data;
	bool havedata;

	// guess2008 - this array holds the climatology for the spinup period
	double dataclim[12];


public:
	Spinup_data(int nyear_loc) {
		nyear=nyear_loc;
		havedata=false;
		data=new double[nyear*12];
		if (!data) fail("Spinup_data::Spinup_data: out of memory");
		thisyear=0;
		havedata=true;
		reset_clim(); // guess2008
	}

	~Spinup_data() {
		if (havedata) delete[] data;
	}

	double& operator[](int month) {

		return data[thisyear*12+month];
	}

	void nextyear() {
		if (thisyear==nyear-1) thisyear=0;
		else thisyear++;
	}

	void firstyear() {
		thisyear=0;
	}

	void get_data_from(double source[][12]) {
		
		int y,m;
		thisyear=0; // guess2008 - ML bugfix
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m]=source[y][m];
			}
		}
	}

	// guess2008 - NEW METHODS 

	void reset_clim() {
		for (int ii = 0; ii < 12; ii++) dataclim[ii] = 0.0;
	}


	void make_clim() {
		
		reset_clim(); // Always reset before calculating

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				dataclim[m] += data[y*12+m] / (double)nyear;
			}
		}
	}


	bool extract_data(double source[][12], const int& startyear, const int& endyear) {
		
		// Populate data with data from the middle of source. 
		// Condition: endyear - startyear + 1 == nyear
		// if startyear == 1 and endyear == 30 then this function is identical to get_data_from above.

		if (endyear < startyear) return false;
		if (endyear - startyear + 1 == nyear) {

			int y,m;
			for (y=startyear-1;y<endyear;y++) {
				for (m=0;m<12;m++) {
					data[(y-(startyear-1))*12+m]=source[y][m];
				}
			}

		} else return false;

		return true;
	}


	void adjust_data(double anom[12], bool additive) {
		
		// Adjust the spinup data to the conditions prevailing at a particular time, as given by 
		// the (additive or multiplicative) anomalies in anom 
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (additive)	
					data[y*12+m] += anom[m];
				else
					data[y*12+m] *= anom[m];
			}
		}

	}
	
	
	// Replace interannual data with the period's climatology.
	void use_clim_data() {
	
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m];
			}
		}
	}


	// Alter variability about the mean climatology
	void adjust_data_variability(const double& factor) {
	
		// factor == 0 gives us the climatology (i.e. generalises use_clim_data above)
		// factor == 1 leaves everything unchanged
		// Remember to check the for negative precip or cloudiness values etc. 
		// after calling this method.

		if (factor == 1.0) return;

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m] + (data[y*12+m] - dataclim[m]) * factor;
			}
		}
	}


	void limit_data(double minval, double maxval) {

		// Limit data to a range
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < minval) data[y*12+m] = minval;
				if (data[y*12+m] > maxval) data[y*12+m] = maxval;
			}
		}

	}
	
	
	void set_min_val(const double& oldval, const double& newval) {

		// Change values < oldval to newval
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < oldval) data[y*12+m] = newval;
			}
		}

	}

	// guess2008 - END OF NEW METHODS


	void detrend_data() {

		int y,m;
		double a,b,anomaly;
		double* annual_mean=new double[nyear];
		double* year_number=new double[nyear];

		if (!annual_mean || !year_number)
			fail("Spinup_driver::detrend_data: out of memory");

		for (y=0;y<nyear;y++) {
			annual_mean[y]=0.0;
			for (m=0;m<12;m++) annual_mean[y]+=data[y*12+m];
			annual_mean[y]/=12.0;
			year_number[y]=y;
		}

		regress(year_number,annual_mean,nyear,a,b);

		for (y=0;y<nyear;y++) {
			anomaly=b*(double)y;
			for (m=0;m<12;m++)
				data[y*12+m]-=anomaly;
		}
		
		// guess2008 - added [] - Clean up
		delete[] annual_mean;
		delete[] year_number;
	}
};

// Constants associated with historical climate data set

// AA CMIP5
const int NYEAR_CMIP5_HIST=156; //1850-01 - 2005-12
const int NYEAR_CMIP5=251; //1850-01 - 2100-12
const int FIRSTHISTYEAR_CMIP5=1850;

// CRU
const int NYEAR_CRU=106;
const int FIRSTHISTYEAR_CRU=1901;

// guess2008
const int NYEAR_HIST=NYEAR_CRU; // guess2008 - CRU TS 3.0 has 106 years of data (1901-2006)
	// number of years of historical climate in CRU and CO2 files (see below)
const int FIRSTHISTYEAR=FIRSTHISTYEAR_CRU;
	// calender year corresponding to first year in CRU climate data set
const int NYEAR_SPINUP_DATA=30;
	// number of years to use for temperature-detrended spinup data set
	// (not to be confused with the number of years to spinup model for, which
	// is read from the ins file)

// Stream pointer to binary CRU historical climate data file (read from ins file)
FILE *in_cru;

// Full pathname of ASCII file containing annual CO2 values (read from ins file)
xtring file_co2;

using namespace GuessOutput;

/// The output channel through which all output is sent
OutputChannel* output_channel;

// GUESSN
// Full pathname of ASCII file containing annual N deposition values (read from ins file)
xtring file_ndep;
// end GUESSN

// Output tables
Table out_cmass, out_anpp, out_dens, out_lai, out_cflux, out_cpool, out_firert, out_runoff, out_speciesheights;

Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
Table out_mrh, out_mnee, out_mwcont_upper, out_mwcont_lower;

// bvoc
Table out_aiso, out_miso, out_amon, out_mmon;

// GUESSN
Table out_cton, out_nmass, out_nsources, out_npool, out_nleach, out_nuptake, out_anppn, out_vmaxnlim, out_nlim, out_nflux;

// GUESSN allometry
Table out_canopyh, out_allometry_ind;
// end GUESSN


// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

/// Yearly CO2 data read from file
/**
 * This object is indexed with calendar years, so to get co2 value for
 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
 * more information.
 */
GlobalCO2File co2;

// Monthly temperature, precipitation and sunshine data for current grid cell
// and historical period
double hist_mtemp[NYEAR_HIST][12];
double hist_mprec[NYEAR_HIST][12];
double hist_msun[NYEAR_HIST][12];

// guess2008
// Monthly frost days, precipitation days and DTR data for current grid cell
// and historical period
double hist_mfrs[NYEAR_HIST][12];
double hist_mwet[NYEAR_HIST][12];
double hist_mdtr[NYEAR_HIST][12];

// GUESSN
// Monthly data on daily dry NHx deposition (kgN/m2/day)
double NHxDryDep[NYEAR_HIST][12];
// Monthly data on daily wet NHx deposition (kgN/m2/day)
double NHxWetDep[NYEAR_HIST][12];
// Monthly data on daily dry NOy deposition (kgN/m2/day)
double NOyDryDep[NYEAR_HIST][12];
// Monthly data on daily wet NOy deposition (kgN/m2/day)
double NOyWetDep[NYEAR_HIST][12];

// CMIP5 - land use input
double hist_frluse[NYEAR_HIST];

//AA CMIP5 Monthly temperature, precipitation and shortwave radiation
double mtemp_cmip5[NYEAR_CMIP5][12];
double mprec_cmip5[NYEAR_CMIP5][12];
double mswrad_cmip5[NYEAR_CMIP5][12];

double mtemp_cru[NYEAR_CRU][12];
double mprec_cru[NYEAR_CRU][12];
double msun_cru[NYEAR_CRU][12];
double mwet_cru[NYEAR_CRU][12];

double mswrad_cru[NYEAR_CRU][12];
//AA CMIP5 Montly CRU climatology 
double clim_mtemp_cru[12];
double clim_mprec_cru[12];
double clim_msun_cru[12];
double clim_swrad_cru[12];
double clim_mwet_cru_1901_1930[12];
double clim_mwet_cru_1961_1990[12];

//AA CMIP5 Montly CMIP5 climatology 
double clim_mtemp_cmip5[12];
double clim_mprec_cmip5[12];
double clim_msun_cmip5[12];
double clim_swrad_cmip5[12];

// AA CMIP5 global strings 
xtring correctionmethod;
xtring gcm;
xtring rcp;
xtring path_cmip5_co2;

// Spinup data sets for current grid cell
Spinup_data spinup_mtemp(NYEAR_SPINUP_DATA);
Spinup_data spinup_mprec(NYEAR_SPINUP_DATA);
Spinup_data spinup_msun(NYEAR_SPINUP_DATA);

// guess2008
// Spinup data sets for monthly frost days, precipitation days and DTR data for 
// current grid cell
Spinup_data spinup_mfrs(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwet(NYEAR_SPINUP_DATA);
Spinup_data spinup_mdtr(NYEAR_SPINUP_DATA);


// Daily temperature, precipitation and sunshine for one year
double dtemp[365],dprec[365],dsun[365];
// bvoc
// Daily diurnal temperature range for one year
double ddtr[365];

// guess2008 - make file_cru and file_cru_misc global variables
xtring file_cru;
xtring file_cru_misc;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double mtemp[12], double mprec[12], double msun[12], double mdtr[12],
					double dtemp[365], double dprec[365], double dsun[365], double ddtr[365]) {
	interp_monthly_means(mtemp, dtemp);
	interp_monthly_totals(mprec, dprec);
	interp_monthly_means(msun, dsun);
	interp_monthly_means(mdtr, ddtr);
}

//AA CMIP5 - climate input
xtring file_cmip5hist;
xtring file_cmip5scen;

//Landuse:

//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
//TimeDataD input code may be put here
TimeDataD LUdata(LOCAL_YEARLY);
TimeDataD Peatdata;
#endif
xtring file_lu, file_peat;
const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHLANDUSE
// Determine land use fraction
// CMIP5 - land use input

bool searchlanduse(double dlon,double dlat,
		   double frluse[NYEAR_HIST]){ 
  
  double rlon;           // longitude in land use input file
  double rlat;           // latitude in land use input file
  double rfhist[NYEAR_CMIP5_HIST]; // land use fraction in historical land use input file
  double rfscen[NYEAR_CMIP5-NYEAR_CMIP5_HIST]; // land use fraction in scenario land use input file
  double rf_const;       // land use fraction in land use input file
  double ds=0.0001;      // allowed difference between gridlist file and land use file 
                         // (rounding/binary differences)
  int iy;                // year counter  
  int ic;                // grid cell counter

  // Open historical land use input file
  xtring file_landusehist=param["file_landusehist"].str;
  FILE* inhist=fopen(file_landusehist,"rt");
  if(!inhist)fail("searchlanduse: could not open land use file %s for input",
		  (char*)file_landusehist);
  // Open scenario land use input file
  xtring file_landusescen=param["file_landusescen"].str;
  FILE* inscen=fopen(file_landusescen,"rt");
  
  // Find (dlon,dlat) and read land use 
  bool grf=false;
  bool grfhist=false;
  bool grfscen=false;
  while(!grfhist&&!grfscen){
    for(ic=0;ic<59191;ic++){
      if(iflandusechange){
	// Check for presence scenario file
	if(!inscen)fail("searchlanduse: could not open land use file %s for input",
			(char*)file_landusescen);
  	// Read historical land use input file
	xtring readline="f,f,";
	char chyear[3];
	sprintf(chyear,"%d",NYEAR_CMIP5_HIST);
	readline+=(xtring)chyear;
	readline+="f";
	readfor(inhist,readline,&rlon,&rlat,rfhist);
	if(rlon>dlon-ds&&rlon<dlon+ds&&rlat>dlat-ds&&rlat<dlat+ds){
	  for(iy=0;iy<NYEAR_CMIP5_HIST;iy++){
	    frluse[iy]=rfhist[iy];
	  }
	  grfhist=true;
	}
	// Read scenario land use input file
	readline="f,f,";
	sprintf(chyear,"%d",NYEAR_CMIP5-NYEAR_CMIP5_HIST);
	readline+=(xtring)chyear;
	readline+="f";
	readfor(inscen,readline,&rlon,&rlat,rfscen);
       	if(rlon>dlon-ds&&rlon<dlon+ds&&rlat>dlat-ds&&rlat<dlat+ds){
	  for(iy=NYEAR_CMIP5_HIST;iy<NYEAR_CMIP5;iy++){
	    frluse[iy]=rfscen[iy-NYEAR_CMIP5_HIST];
	  }
  	  grfscen=true;
	}
      }
      else{
	xtring readline="f,f,f";
	readfor(inhist,readline,&rlon,&rlat,&rf_const);
	if(rlon>dlon-ds&&rlon<dlon+ds&&rlat>dlat-ds&&rlat<dlat+ds){
	  for(iy=0;iy<NYEAR_HIST;iy++){
	    frluse[iy]=rf_const;
	  }
	  grfhist=true;
	  grfscen=true;
	}
      }
    }
  }
  fclose(inhist);
  fclose(inscen);

  if(grfhist&&grfscen)grf=true;
  
  return grf;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 
// CMIP5IO 
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Regress_data
// //AA cmip5
void regress_data(double* x,double* y,int n,double& a,double& b) {

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
	delta=0.0;
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



///////////////////////////////////////////////////////////////////////////////////////
// calculate_swrad
// Calculates swrad from sunshine//AA cmip5
void calculate_swrad(double msun[NYEAR_CRU][12], double lon, double lat, double swrad[NYEAR_CRU][12])
{
	//{31,28,31,30,31,30,31,31,30,31,30,31};
	const double days_per_month[12]={31,28,31,30,31,30,31,31,30,31,30,31};
	const double day_start_month[13]={0,31,59,90,120,151,181,212,243,273,304,334,365};
	
	int y, m;
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
	

	double dummy1[1][12], dummy2[1][12];//, ddummy1[365], ddummy2[365], dsun[365];
	for (y=0; y<NYEAR_CRU; y++)
	{
		for (m=0; m<12; m++)
		{
			swrad[y][m]=0.0;
			dummy1[0][m]=m;
			dummy2[0][m]=m;
		}
	}

	double delta; // solar declination angle (radians)
	double rs_day;
	double w;	
			//variabler tidigare i climate:

	double qo[365];
	double u[365];
	double v[365];
	double hh[365];
	double sinehh[365];
	//double daylength_save;

	for (int year = 0; year < NYEAR_CRU; year++) 
	{

	/*	interp_climate(dummy1[0],
				dummy2[0],msun[year],
				ddummy1,ddummy2,dsun);
	*/	
	/*	if (year==0) //debug
		{
			double msuntot=0.0, dsuntot=0.0;
			for (m=0;m<12;m++)
				msuntot+=msun[year][m];
			
			for (int d=0;d<365;d++)
				dsuntot+=dsun[d];

			printf("msuntot: %7.3f dsuntot: %7.3f\n",msuntot/12, dsuntot/365);
		}
		*/

		for (m = 0; m < 12; m++) 
		{
			for (int nday=(int)day_start_month[m]; nday<(int)day_start_month[m+1]; nday++) //(int nday=0; nday<=365; nday++)
			{
				if (year==0)
				{
					qo[nday]=QOO*(1.0+2.0*0.01675*cos(2.0*PI*((double)nday+0.5)/365.0)); // Eqn 2
						
					delta=-23.4*DEGTORAD*cos(2.0*PI*((double)nday+10.5)/365.0); // Eqn 4
						
					u[nday]=sin(lat*DEGTORAD)*sin(delta); // Eqn 9
					v[nday]=cos(lat*DEGTORAD)*cos(delta); // Eqn 10

					if (u[nday]>=v[nday])
						hh[nday]=PI; // polar day
					else if (u[nday]<=-v[nday])
						hh[nday]=0.0; // polar night
					else 
						hh[nday]=acos(-u[nday]/v[nday]); // Eqn 11

					sinehh[nday]=sin(hh[nday]);
				}
				
				w=(C+D*msun[year][m]/100.0)*qo[nday]; // Eqn 13 dsun[nday]

				rs_day=2.0*w*(u[nday]*hh[nday]+v[nday]*sinehh[nday])*K; // Eqn 14
				
				//if (daylength_save>0)
				//{
					swrad[year][m]+=((rs_day)/(24*3600))/days_per_month[m];//daylength_save*3600.0); //cswrad[m]+=((rs_day/(days_per_month[m]*30))/((1.0-BETA)*daylength_save*3600.0));
				//}	
			}			
		}
	}

	for (y=0; y<NYEAR_CRU; y++)
	{
		for (m=0; m<12; m++)
		{
			if (swrad[y][m]<0.0) swrad[y][m]=0.0; //Anders A debug
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// CREATECLIMATOLOGY_CRU
// Creates a 30-year climatology from CRU and creates swrad climatology from CRU sunshine//AA cmip5

void createclimatology_cru(double mtemp[NYEAR_CRU][12],double mprec[NYEAR_CRU][12],
					   double msun[NYEAR_CRU][12], double mwet[NYEAR_CRU][12],double cruswrad[NYEAR_CRU][12],
					   double ctemp[12],double cprec[12], double csun[12], double cswrad[12], 
					   double cwet_1901_1930[12], double cwet_1961_1990[12]) 
{

	int y ,m;

	for (m = 0; m < 12; m++) 
	{
		ctemp[m] = 0.0; 
		cprec[m] = 0.0;  
		csun[m] = 0.0;
		cswrad[m] = 0.0; 
		cwet_1901_1930[m]= 0.0;
	    cwet_1961_1990[m]=0.0;
	}

	for (y = 60; y < 90; y++) 
	{
		for (m = 0; m < 12; m++) 
		{
			ctemp[m] += mtemp[y][m] / 30.0; 
			cprec[m]+= mprec[y][m] / 30.0;
			csun[m] += msun[y][m] / 30.0;
			cswrad[m] += cruswrad[y][m] / 30.0;
			cwet_1961_1990[m]+=mwet[y][m] / 30.0;
		}
	}

	for (y = 0; y < 30; y++) 
	{
		for (m = 0; m < 12; m++) 
		{
			cwet_1901_1930[m]+=mwet[y][m] / 30.0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// CREATECLIMATOLOGY_CMIP5
// Creates a 30-year climatology from CRU and creates swrad climatology from CRU sunshine//AA cmip5

void createclimatology_cmip5(double mtemp[NYEAR_HIST][12],double mprec[NYEAR_HIST][12],
					   double mswrad[NYEAR_HIST][12], double ctemp[12],double cprec[12],
					   double cswrad[12]) 
{
	int m, y;
	
	for (m = 0; m < 12; m++) 
	{
		ctemp[m] = 0.0; 
		cprec[m] = 0.0;  
		cswrad[m] = 0.0; 		
	}

	for (y=111;y<141;y++) //1961-1990
	{
		for (m = 0; m < 12; m++) 
		{
			ctemp[m] += mtemp[y][m] / 30.0; 
			cprec[m]+= mprec[y][m] / 30.0;
			cswrad[m] += mswrad[y][m] / 30.0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCMIP5HIST
// Determine temp, precip and shortwave radiation //AA cmip5
 
bool searchcmip5hist(char* cmip5histark,double dlon,double dlat,
	double mtemp[NYEAR_CMIP5][12],double mprec[NYEAR_CMIP5][12],
	double mswrad[NYEAR_CMIP5][12]) {

	// Archive object. Definition in cmip5 header file, 
	//#if defined CM5rcp85
	//	Ipsl_cm5a_lr_historical_1850_2005_r1i1p1Archive ark;
	//#endif
		Cmip5_histArchive ark;

	int target_ilon=dlon*10.0;
	int target_ilat=dlat*10.0;

	int y,m;

	// Try block to catch any unexpected errors
	try {
		//#if defined CM5rcp85
		//	Ipsl_cm5a_lr_historical_1850_2005_r1i1p1 data;
		//#endif
		
		Cmip5_hist data;
	//	Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cmip5histark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldnt rewind
				return false;
			}
		}
		else
			return false;

		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		//soilcode=(int)data.soilcode[0];

		for (y=0;y<NYEAR_CMIP5_HIST;y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[y*12+m];//*0.1; // now degC
				mprec[y][m] = data.mprec[y*12+m];//*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				mswrad[y][m]  = data.mswrad[y*12+m];//*0.1;   // % sun 
			}
		}

		// Close the archive
		ark.close();

		return true;
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCMIP5SCEN
// Determine temp, precip and shortwave radiation //AA cmip5
 
bool searchcmip5scen(char* cmip5scenark,double dlon,double dlat,
	double mtemp[NYEAR_CMIP5][12],double mprec[NYEAR_CMIP5][12],
	double mswrad[NYEAR_CMIP5][12]) {

	// Archive object. Definition in cmip5 header file, 
	//#if defined CM5rcp85
	//	Ipsl_cm5a_lr_rcp85_2006_2100_r1i1p1Archive ark;
	//#endif
	Cmip5_scenArchive ark;

	int target_ilon=dlon*10.0;
	int target_ilat=dlat*10.0;

	int y,m;

	// Try block to catch any unexpected errors
	try {		
		//#if defined CM5rcp85
		//	Ipsl_cm5a_lr_rcp85_2006_2100_r1i1p1 data;
		//#endif
		
		Cmip5_scen data;
	//	Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cmip5scenark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldnt rewind
				return false;
			}
		}
		else
			return false;

		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		//soilcode=(int)data.soilcode[0];

		int year=0;
		for (y=NYEAR_CMIP5_HIST;y<NYEAR_CMIP5;y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[year*12+m];//*0.1; // now degC
				mprec[y][m] = data.mprec[year*12+m];//*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				mswrad[y][m]  = data.mswrad[year*12+m];//*0.1;   // % sun 
			}
			year++;
		}

	//				for ( m = 0; m < 12; m++) 
	//	{
	//	printf("m: %d make: %4.2f\n",m, mswrad[170][m]);
	//	}

		// Close the archive
		ark.close();

		return true;
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// makeCMIP5data
// Create scenario and historical temp, precip and shortwave radiation  //AA cmip5

void makeCMIP5data(double cmip5temp[NYEAR_CMIP5][12],double cmip5prec[NYEAR_CMIP5][12], double cmip5swrad[NYEAR_CMIP5][12],
				   double NHxW[NYEAR_CMIP5][12],double NHxD[NYEAR_CMIP5][12],double NOyW[NYEAR_CMIP5][12],double NOyD[NYEAR_CMIP5][12],
				   double ctemp_cmip5[12], double cprec_cmip5[12],double cswrad_cmip5[12],
				   double ctemp_cru[12], double cprec_cru[12],double csun_cru[12],double cswrad_cru[12], double cwet_cru_1901_1930[12], double cwet_cru_1961_1990[12],
				   double crutemp[NYEAR_CRU][12],double cruprec[NYEAR_CRU][12],double crusun[NYEAR_CRU][12], double cruwet[NYEAR_CRU][12], double cruswrad[NYEAR_CRU][12],
				   double temp[NYEAR_CMIP5][12],double prec[NYEAR_CMIP5][12],double sun[NYEAR_CMIP5][12], double wet[NYEAR_CMIP5][12] ) 
{
	bool cmip5data_no_correction=false;
	bool correct_yearly=false;
	bool correct_monthly=false;
	bool merge_at_2006=false;
	bool cmip5_trend_cru_var=false;
	
	/*if (correctionmethod=="c1")
		cmip5data_no_correction=true;
	else if (correctionmethod=="c2")
		correct_yearly=true;
	else if (correctionmethod=="c3")
		correct_monthly=true;
	else if (correctionmethod=="c4")
		merge_at_2006=true;
	else if (correctionmethod=="c5")
		cmip5_trend_cru_var=true;
	else
		fail("\nNo valid correctionmethod choice\n");
*/
	int y;
	int m;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Cmip5 data, no correction
	if (correctionmethod=="c1")
	{
		//printf("in c1\n");
		for (y = 0; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m];
				prec[y][m]=cmip5prec[y][m];
				sun[y][m]=cmip5swrad[y][m];
			}
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// correctyearly
	// Creates CMIP5 historical and scenario data.
	// The correction is based on the annual difference 
	// between the 1961-1990 climatologies between CRU and CMIP5 historical. 
	else if (correctionmethod=="c2")
	{	
		double cruatemp=0.0;
		double cruaprec=0.0;
		double cruaswrad=0.0;
		double cmip5atemp=0.0;
		double cmip5aprec=0.0;
		double cmip5aswrad=0.0;

		//find yearly climatologies
		for (m = 0; m < 12; m++) 
		{
			cruatemp+=ctemp_cru[m]/12;
			cruaprec+=cprec_cru[m]/12;
			cruaswrad+=cswrad_cru[m]/12;
			
			cmip5atemp+=ctemp_cmip5[m]/12;
			cmip5aprec+=cprec_cmip5[m]/12;
			cmip5aswrad+=cswrad_cmip5[m]/12;
		}
	
		// Correct the CMIP5 data 
		for (y = 0; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m]-cmip5atemp+cruatemp;

				if (cmip5aprec <= 5 || cruaprec ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
					prec[y][m]=cmip5prec[y][m]-cmip5aprec+cruaprec;
				else
					prec[y][m]=(cmip5prec[y][m]/cmip5aprec)*cruaprec;

				if (cmip5aswrad == 0 || cruaswrad ==0)
					sun[y][m]=cmip5swrad[y][m]-cmip5aswrad+cruaswrad;
				else
					sun[y][m]=(cmip5swrad[y][m]/cmip5aswrad)*cruaswrad;

				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
	}// end correctyearly

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// correctmonthly
	// Creates CMIP5 historical and scenario data. 
	// The correction is based on the monthly difference / ratio
	// between the 1961-1990 climatologies between CRU and CMIP5 historical.
	else if (correctionmethod=="c3")
	{
		//printf("in correctmonthly\n");
		for (y = 0; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m];
				
				if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
					prec[y][m]=cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m];
				else
					prec[y][m]=(cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m];

				if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
					sun[y][m]=cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m];
				else	
					sun[y][m]=(cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m];

				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}

	}// end correctmonthly

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Merge at 2006
	// Creates CMIP5 historical and scenario data. 
	// CMIP5 Anomalies are superimposed on cru climatology 
	else if (correctionmethod=="c4")
	{
		//printf("in merge at 2006\n");
		// Calculate swrad from cloudiness
		const int n=30;;
		double x[n], tempy[n]; //, swrady[n], precy[n]; 
		double a_temp, b_temp;//, a_prec, b_prec, a_swrad, b_swrad;
		double anom_temp;//, anom_prec, anom_swrad;
		
		//debug
		/*for(y=0;y<NYEAR_CRU;y++)
		{
			for(m=0;m<12;m++)
			{
				cruswrad[y][m]=crusun[y][m];
			}
		}
		*/
		// Create cru data for 1850-1900 use detrended 1901-1930 climatology twice
		
		// create vectors for regression
		for (y=0;y<30;y++)
		{
			x[y]=(double)y;
			
			tempy[y]=0.0;
			//precy[y]=0.0;
			//swrady[y]=0.0;

			for (m = 0; m < 12; m++) 
			{
				tempy[y]+=crutemp[y][m]/12;
				//precy[y]+=cruprec[y][m]/12;
				//swrady[y]+=cruswrad[y][m]/12;
			}	
		}
		
		// regress data to remove trend
		regress_data(x,tempy,n,a_temp,b_temp);
		//regress_data(x,precy,n,a_prec,b_prec);
		//regress_data(x,swrady,n,a_swrad,b_swrad);

		// remove trend and fill 1850-1879 with detrended cru data
		for (y=0;y<30;y++)
		{
			anom_temp=(double)y*b_temp;
			//anom_prec=(double)y*b_prec;
			//anom_swrad=(double)y*b_swrad;
			
			for (m = 0; m < 12; m++) 
			{
				//1850-1879
				temp[y][m]=crutemp[y][m]-anom_temp;
				prec[y][m]=cruprec[y][m];// no detrend -anom_prec;
				sun[y][m]=cruswrad[y][m];//no detrend -anom_swrad;
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		// remove trend and fill 1880-1900 with detrended cru data
		for (y=30;y<51;y++)
		{
			anom_temp=(double)(y-30)*b_temp;
			//anom_prec=(double)(y-30)*b_prec;
			//anom_swrad=(double)(y-30)*b_swrad;
			
			for (m = 0; m < 12; m++) 
			{
				//1880-1900
				temp[y][m]=crutemp[y-30][m]-anom_temp;
				prec[y][m]=cruprec[y-30][m];//no detrend -anom_prec;
				sun[y][m]=cruswrad[y-30][m];//no detrend -anom_swrad;
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		
		// Fill 1901-2005 with raw CRU data
		for (y = 51; y < 156; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=crutemp[y-51][m];
				prec[y][m]=cruprec[y-51][m];
				sun[y][m]=cruswrad[y-51][m];
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
			}
		}

		// Fill with bias corrected cmip5 data 2006-2100
		for (y = 156; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m];
				
				if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
					prec[y][m]=cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m];
				else
					prec[y][m]=(cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m];

				if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
					sun[y][m]=cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m];
				else	
					sun[y][m]=(cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m];


				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
					
			}
		}
	}// end merge at 2006
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Cmip5 trend cru variability
	// Creates CMIP5 historical and scenario data. 
	// 1850-2005 from on CRU as above, 2006-2085 based on offset corrected CMIP5 data (as above) running average.
	// On top of the running average CRU 1961-1990 anomalies are added and recycled. Simulation stops early because of the running average calculation.
	else if (correctionmethod=="c5" || correctionmethod=="c6")
	{	
		//1 repeat almost all steps in merge 2006
		///////////////////////////////////////////////////
		//from merge at 2006
		
		// Calculate swrad from cloudiness
		const int n=30;
		double x[n], tempy[n], precy[n], swrady[n];
		double a_temp, b_temp, a_prec, b_prec, a_swrad, b_swrad;
		double anom_temp, anom_prec, anom_swrad;
				
		// Create cru data for 1850-1900 use detrended 1901-1930 climatology twice
		
		// create vectors for regression
		for (y=0;y<30;y++)
		{
			x[y]=(double)y;
			
			tempy[y]=0.0;
			//precy[y]=0.0;
			//swrady[y]=0.0;

			for (m = 0; m < 12; m++) 
			{
				tempy[y]+=crutemp[y][m]/12;
				//precy[y]+=cruprec[y][m]/12;
				//swrady[y]+=cruswrad[y][m]/12;
			}	
		}		
		
		// regress data to remove trend
		regress_data(x,tempy,n,a_temp,b_temp);
		//regress_data(x,precy,n,a_prec,b_prec);
		//regress_data(x,swrady,n,a_swrad,b_swrad);

		// remove trend and fill 1850-1879 with detrended cru data
		for (y=0;y<30;y++)
		{
			anom_temp=(double)y*b_temp;
			//anom_prec=(double)y*b_prec;
			//anom_swrad=(double)y*b_swrad;
			
			for (m = 0; m < 12; m++) 
			{
				//1850-1879
				temp[y][m]=crutemp[y][m]-anom_temp;
				prec[y][m]=cruprec[y][m];//-anom_prec; no detrend
				sun[y][m]=cruswrad[y][m];//-anom_swrad; no detrend
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		// remove trend and fill 1880-1900 with detrended cru data
		for (y=30;y<51;y++)
		{
			anom_temp=(double)(y-30)*b_temp;
			//anom_prec=(double)(y-30)*b_prec;
			//anom_swrad=(double)(y-30)*b_swrad;
			
			for (m = 0; m < 12; m++) 
			{
				//1880-1900
				temp[y][m]=crutemp[y-30][m]-anom_temp;
				prec[y][m]=cruprec[y-30][m];//-anom_prec; //no detrend
				sun[y][m]=cruswrad[y-30][m];//-anom_swrad; // no detrend
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		
		// Fill 1901-2005 with raw CRU data
		for (y = 51; y < 156; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=crutemp[y-51][m];
				prec[y][m]=cruprec[y-51][m];
				sun[y][m]=cruswrad[y-51][m];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////
		//end repeat steps of merge at 2006
				
		// Create offset annual corrected cmip5 data 1850-2100
		double no_temp[NYEAR_CMIP5], no_prec[NYEAR_CMIP5],no_swrad[NYEAR_CMIP5]; //no offset cmip5 climate
		double av_temp[NYEAR_CMIP5], av_prec[NYEAR_CMIP5],av_swrad[NYEAR_CMIP5]; // yearly moving average cmip5 climate
		double cruanom_temp[30][12], cruanom_prec[30][12],cruanom_swrad[30][12]; // cru 61-90 anomalies
		double trend_temp;//, trend_prec, trend_swrad;
		double average_prec=0.0, average_swrad=0.0;		
		
		if (correctionmethod=="c5")
		{
			for (y = 0; y < NYEAR_CMIP5; y++) 
			{  
				no_temp[y]=0.0;
				no_prec[y]=0.0;
				no_swrad[y]=0.0;
				
				for (m = 0; m < 12; m++) 
				{
					no_temp[y]+=( cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m] )/12;

					if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
						no_prec[y]+=( cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m] )/12;
					else
						no_prec[y]+=( (cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m] )/12;

					if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
						no_swrad[y]+=( cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m] )/12;
					else	
						no_swrad[y]+=( (cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m] )/12;
	
				}
			}
		
			// calculate moving average	
			for (y=15; y< NYEAR_CMIP5-15; y++)
			{
				av_temp[y]=0.0;
				av_prec[y]=0.0;
				av_swrad[y]=0.0;

				for (int i=y-15; i<=y+15;i++)
				{
					av_temp[y]+=no_temp[i]/31;
					av_prec[y]+=no_prec[i]/31;
					av_swrad[y]+=no_swrad[i]/31;
				}
			}

			// calculate 1961-1990 cru anomalies
			// create vectors for regression
			for (y=60;y<90;y++)
			{
				x[y-60]=(double)(y-60);
				
				tempy[y-60]=0.0;
				//precy[y-60]=0.0;
				//swrady[y-60]=0.0;

				for (m = 0; m < 12; m++) 
				{
					tempy[y-60]+=crutemp[y][m]/12;
					//precy[y-60]+=cruprec[y][m]/12;
					//swrady[y-60]+=cruswrad[y][m]/12;
					average_prec+=cruprec[y][m]/360;
					average_swrad+=cruswrad[y][m]/360;
				}	
			}
			
			// regress data
			regress_data(x,tempy,n,a_temp,b_temp);
			//regress_data(x,precy,n,a_prec,b_prec);
			//regress_data(x,swrady,n,a_swrad,b_swrad);

			// find anomalies by removing the trend
			for (y=60;y<90;y++)
			{
				trend_temp=(double)(y-60)*b_temp+a_temp;
				//trend_prec=(double)(y-60)*b_prec+a_prec;
				//trend_swrad=(double)(y-60)*b_swrad+a_swrad;
				
				for (m = 0; m < 12; m++) 
				{
					cruanom_temp[y-60][m]=crutemp[y][m]-trend_temp;
					cruanom_prec[y-60][m]=cruprec[y][m]-average_prec;// no detrend trend_prec;
					cruanom_swrad[y-60][m]=cruswrad[y][m]-average_swrad;//no detrend trend_swrad;
				}
			}
			
			//for (y=0;y<30;y++) printf("cruanom_prec: %5.3f\n",cruanom_prec[y][10]);

			// add cru 61-90 anomalies to 2006-2085 CMIP5 "trend"
			int anomyear=0;
			for (y=156; y<NYEAR_CMIP5-15;y++)
			{

				for (m=0;m<12;m++)
				{
					temp[y][m]=av_temp[y]+cruanom_temp[anomyear][m];
					prec[y][m]=av_prec[y]+cruanom_prec[anomyear][m];
					sun[y][m]=av_swrad[y]+cruanom_swrad[anomyear][m];

					// Limit very low precip amounts because negligible precipitation causes problems 
					// in the prdaily function (infinite loops). 
					if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
					if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
					
				}
				anomyear++;
				if (anomyear>=30) anomyear=0;
			}
		}


		///////////////////////////////////////////////
		// difference from c5 is that moving averages are calculated on monthly basis. Also, CRU climatology is now de seasonalized before addition to cmip5 trend
		if (correctionmethod=="c6")	{

			double mtrend_temp[12];
			double a_mtemp[12],b_mtemp[12];
			double average_mprec[12], average_mswrad[12];
			double mav_temp[NYEAR_CMIP5][12], mav_prec[NYEAR_CMIP5][12],mav_swrad[NYEAR_CMIP5][12]; //monthly moving average cmip5 climate
			double mno_temp[NYEAR_CMIP5][12], mno_prec[NYEAR_CMIP5][12],mno_swrad[NYEAR_CMIP5][12]; // no offset monthly cmip5 climate
			
			//for (y=0;y<30;y++)
			//{
				for (m=0;m<12;m++)
				{
					//mtrend_temp[m]=0.0;
					//mtempy[y][m]=0.0;
					average_mprec[m]=0.0;
					average_mswrad[m]=0.0;
				}
			//}
			
			for (y = 0; y < NYEAR_CMIP5; y++) 
			{  
				
				for (m = 0; m < 12; m++) 
				{
					mno_temp[y][m]=( cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m] );

					if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
						mno_prec[y][m]=( cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m] );
					else
						mno_prec[y][m]=( (cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m] );

					if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
						mno_swrad[y][m]=( cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m] );
					else	
						mno_swrad[y][m]=( (cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m] );
				}
			}
		
			// calculate moving average	
			for (y=15; y< NYEAR_CMIP5-15; y++)
			{
				for (m=0; m<12; m++)
				{
					mav_temp[y][m]=0.0;
					mav_prec[y][m]=0.0;
					mav_swrad[y][m]=0.0;

					for (int i=y-15; i<=y+15;i++)
					{
						mav_temp[y][m]+=mno_temp[i][m]/31;
						mav_prec[y][m]+=mno_prec[i][m]/31;
						mav_swrad[y][m]+=mno_swrad[i][m]/31;
					}
				}
			}
		

			// calculate 1961-1990 cru anomalies
			// create vectors for regression
			for (y=60;y<90;y++)
			{
				x[y-60]=(double)(y-60);

				for (m = 0; m < 12; m++) 
				{
					average_mprec[m]+=cruprec[y][m]/30;
					average_mswrad[m]+=cruswrad[y][m]/30;
				}	
			}
				
			// regress data, now 12 times
			for (m=0;m<12;m++)
			{
				for (y=60;y<90;y++)
					tempy[y-60]=crutemp[y][m];

				regress_data(x,tempy,n,a_mtemp[m],b_mtemp[m]);
			}
			//regress_data(x,precy,n,a_prec,b_prec);
			//regress_data(x,swrady,n,a_swrad,b_swrad);

			// find anomalies by removing the trend
			for (y=60;y<90;y++)
			{	
				for (m = 0; m < 12; m++) 
				{
					trend_temp=(double)(y-60)*b_mtemp[m]+a_mtemp[m];
					cruanom_temp[y-60][m]=crutemp[y][m]-trend_temp;
					cruanom_prec[y-60][m]=cruprec[y][m]-average_mprec[m];// no detrend trend_prec;
					cruanom_swrad[y-60][m]=cruswrad[y][m]-average_mswrad[m];//no detrend trend_swrad;
				}
			}
			
			//for (y=0;y<30;y++) printf("cruanom_prec: %5.3f\n",cruanom_prec[y][10]);

			// add cru 61-90 anomalies to 2006-2085 CMIP5 "trend"

			int anomyear=0;
			for (y=156; y<NYEAR_CMIP5-15;y++)
			{
				for (m=0;m<12;m++)
				{
					temp[y][m]=mav_temp[y][m]+cruanom_temp[anomyear][m];
					prec[y][m]=mav_prec[y][m]+cruanom_prec[anomyear][m];
					sun[y][m]=mav_swrad[y][m]+cruanom_swrad[anomyear][m];

					// Limit very low precip amounts because negligible precipitation causes problems 
					// in the prdaily function (infinite loops). 
					if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
					if (sun[y][m] <= 0.0) sun[y][m] = 0.0;				
				}
				anomyear++;
				if (anomyear>=30) anomyear=0;
			}	
		}

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// end simulation 15 years short!!!

	}// end cmip5 trend cru var

	// Standard CRU simulation using cloudiness/sunshine instead of shortwave radiation as above.
	// The daily interpolation does not preserve the sums of the variables, therefore using CRU shortwave radiation
	// is not identical to using cloudiness/sunshine
	else if (correctionmethod=="c7")
	{
				// Calculate swrad from cloudiness
		const int n=30;
		double x[n], tempy[n];
		double a_temp, b_temp, a_prec, b_prec, a_swrad, b_swrad;
		double anom_temp, anom_prec, anom_swrad;
			
		// Create cru data for 1850-1900 use detrended 1901-1930 climatology twice
		
		// create vectors for regression
		for (y=0;y<30;y++)
		{
			x[y]=(double)y;
			
			tempy[y]=0.0;
			//precy[y]=0.0;
			//swrady[y]=0.0;

			for (m = 0; m < 12; m++) 
			{
				tempy[y]+=crutemp[y][m]/12;
				//precy[y]+=cruprec[y][m]/12;
				//swrady[y]+=cruswrad[y][m]/12;
			}	
		}	
		
		// regress data to remove trend
		regress_data(x,tempy,n,a_temp,b_temp);

		// remove trend and fill 1850-1879 with detrended cru data
		for (y=0;y<30;y++)
		{
			anom_temp=(double)y*b_temp;
			
			for (m = 0; m < 12; m++) 
			{
				//1850-1879
				temp[y][m]=crutemp[y][m]-anom_temp;
				prec[y][m]=cruprec[y][m];//-anom_prec; no detrend
				sun[y][m]=crusun[y][m];//-anom_swrad; no detrend
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		// remove trend and fill 1880-1900 with detrended cru data
		for (y=30;y<51;y++)
		{
			anom_temp=(double)(y-30)*b_temp;
			
			for (m = 0; m < 12; m++) 
			{
				//1880-1900
				temp[y][m]=crutemp[y-30][m]-anom_temp;
				prec[y][m]=cruprec[y-30][m];//-anom_prec; //no detrend
				sun[y][m]=crusun[y-30][m];//-anom_swrad; // no detrend
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
		
		// Fill 1901-2005 with raw CRU data
		for (y = 51; y < 156; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=crutemp[y-51][m];
				prec[y][m]=cruprec[y-51][m];
				sun[y][m]=crusun[y-51][m];
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// correctmonthly
	// Creates CMIP5 historical and scenario data. 
	// The correction is based on the monthly difference / ratio
	// between the 1961-1990 climatologies between CRU and CMIP5 historical.
	// c8: correct annual and seasonal offset (c3) and keeping temperature and CO2 at pre industrial level (1861-1870).
	// c9: correct annual and seasonal offset (c3) and keeping temperature and N deposition at pre industrial level (1861-1870).
	// c10: correct annual and seasonal offset (c3) and keeping CO2 and N deposition at pre industrial level (1861-1870). 
	else if (correctionmethod=="c8" || correctionmethod=="c9" || correctionmethod=="c10") {
		
		//printf("in correctmonthly\n");
		for (y = 0; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m];
				
				if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
					prec[y][m]=cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m];
				else
					prec[y][m]=(cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m];

				if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
					sun[y][m]=cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m];
				else	
					sun[y][m]=(cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m];

				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}

		if (correctionmethod=="c9" || correctionmethod=="c10" ) {

			// Create pre industrial N deposition data for 1871-2100 using detrended 1861-1870 N deposition

			const int n=10;
			double x[n], ndepNHxWy[n], ndepNHxDy[n], ndepNOyWy[n], ndepNOyDy[n];
			double a_ndepNHxW, b_ndepNHxW, a_ndepNHxD, b_ndepNHxD, a_ndepNOyW, b_ndepNOyW, a_ndepNOyD, b_ndepNOyD;
			double anom_ndepNHxW, anom_ndepNHxD, anom_ndepNOyW, anom_ndepNOyD;
		
			// create vectors for regression
			for (y=0;y<n;y++) {

				x[y]=(double)y+11.0;
			
				ndepNHxWy[y]=0.0;
				ndepNHxDy[y]=0.0;
				ndepNOyWy[y]=0.0;
				ndepNOyDy[y]=0.0;

				for (m = 0; m < 12; m++) {

					// 1861 == year 11
					
					ndepNHxWy[y]+=NHxW[11+y][m]/12.0;
					ndepNHxDy[y]+=NHxD[11+y][m]/12.0;
					ndepNOyWy[y]+=NOyW[11+y][m]/12.0;
					ndepNOyDy[y]+=NOyD[11+y][m]/12.0;
				}	
			}
		
			// regress data to remove trend
			regress_data(x,ndepNHxWy,n,a_ndepNHxW,b_ndepNHxW);
			regress_data(x,ndepNHxDy,n,a_ndepNHxD,b_ndepNHxD);
			regress_data(x,ndepNOyWy,n,a_ndepNOyW,b_ndepNOyW);
			regress_data(x,ndepNOyDy,n,a_ndepNOyD,b_ndepNOyD);

			// remove trend and fill 1871-2100 with detrended pre industrial data
			for (y=21;y<NYEAR_CMIP5;y++) {

				anom_ndepNHxW=(double)(y%10)*b_ndepNHxW;
				anom_ndepNHxD=(double)(y%10)*b_ndepNHxD;
				anom_ndepNOyW=(double)(y%10)*b_ndepNOyW;
				anom_ndepNOyD=(double)(y%10)*b_ndepNOyD;
			
				for (m = 0; m < 12; m++) {

					NHxW[y][m]=NHxW[11+y%10][m]-anom_ndepNHxW;
					NHxD[y][m]=NHxD[11+y%10][m]-anom_ndepNHxD;
					NOyW[y][m]=NOyW[11+y%10][m]-anom_ndepNOyW;
					NOyD[y][m]=NOyD[11+y%10][m]-anom_ndepNOyD;
				}
			}
		}
		if (correctionmethod=="c8" || correctionmethod=="c10") {
	
			// Create pre industrial CO2 data for 1871-2100 using detrended 1861-1870 CO2

			const int n=10;
			double x[n], co2y[n];
			double a_co2, b_co2;
			double anom_co2;
		
			// create vectors for regression
			for (y=0;y<n;y++) {

				x[y]=(double)y+11.0;
				
				co2y[y]=co2[1861+y];
			}
		
			// regress data to remove trend
			regress_data(x,co2y,n,a_co2,b_co2);

			// remove trend and fill 1871-2100 with detrended pre industrial data
			for (y=21;y<NYEAR_CMIP5;y++) {

				anom_co2=(double)(y%10)*b_co2;

				co2[1850+y]=co2[1861+y%10]-anom_co2;
			}
		}	
		if (correctionmethod=="c8" || correctionmethod=="c9") {
	
			// Create pre industrial temperature data for 1871-2100 using detrended 1861-1870 temperature

			const int n=10;
			double x[n], tempy[n];
			double a_temp, b_temp;
			double anom_temp;
		
			// create vectors for regression
			for (y=0;y<n;y++) {

				x[y]=(double)y+11.0;
			
				tempy[y]=0.0;

				for (m = 0; m < 12; m++) {

					// 1861 == year 11
					tempy[y]+=temp[11+y][m]/12;
				}	
			}
		
			// regress data to remove trend
			regress_data(x,tempy,n,a_temp,b_temp);

			// remove trend and fill 1871-2100 with detrended pre industrial data
			for (y=21;y<NYEAR_CMIP5;y++) {

				anom_temp=(double)(y%10)*b_temp;
			
				for (m = 0; m < 12; m++) {

					temp[y][m]=temp[11+y%10][m]-anom_temp;
				}
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// correctmonthly
	// Creates CMIP5 historical and scenario data. 
	// The correction is based on the monthly difference / ratio
	// between the 1961-1990 climatologies between CRU and CMIP5 historical.
	// c8: correct annual and seasonal offset (c3) and keeping temperature and CO2 at pre industrial level (1861-1870).
	// c9: correct annual and seasonal offset (c3) and keeping temperature and N deposition at pre industrial level (1861-1870).
	// c10: correct annual and seasonal offset (c3) and keeping CO2 and N deposition at pre industrial level (1861-1870). 
	else if (correctionmethod=="c11") {
		
		//printf("in correctmonthly\n");
		for (y = 0; y < NYEAR_CMIP5; y++) 
		{  
			for (m = 0; m < 12; m++) 
			{
				temp[y][m]=cmip5temp[y][m]-ctemp_cmip5[m]+ctemp_cru[m];
				
				if (cprec_cmip5[m] <= 5 || cprec_cru[m] ==0) // divide by zero fix, also solves problem with low climatology precip/swrad
					prec[y][m]=cmip5prec[y][m]-cprec_cmip5[m]+cprec_cru[m];
				else
					prec[y][m]=(cmip5prec[y][m]/cprec_cmip5[m])*cprec_cru[m];

				if (cswrad_cmip5[m] == 0 || cswrad_cru[m] ==0)
					sun[y][m]=cmip5swrad[y][m]-cswrad_cmip5[m]+cswrad_cru[m];
				else	
					sun[y][m]=(cmip5swrad[y][m]/cswrad_cmip5[m])*cswrad_cru[m];

				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (prec[y][m] <= 1.0) prec[y][m] = 0.0;
				if (sun[y][m] <= 0.0) sun[y][m] = 0.0;
			}
		}
	
		// Create pre industrial temperature data for 1871-2100 using detrended 1861-1870 temperature

		const int n=10;
		double x[n], tempy[n];
		double a_temp, b_temp;
		double anom_temp;
		
		// create vectors for regression
		for (y=0;y<n;y++) {

			x[y]=(double)y+11.0;
			
			tempy[y]=0.0;

			for (m = 0; m < 12; m++) {

				// 1861 == year 11
				tempy[y]+=temp[11+y][m]/12;
			}	
		}
		
		// regress data to remove trend
		regress_data(x,tempy,n,a_temp,b_temp);

		// remove trend and fill 1871-2100 with detrended pre industrial data
		for (y=21;y<NYEAR_CMIP5;y++) {

			anom_temp=(double)(y%10)*b_temp;
			
			for (m = 0; m < 12; m++) {

				temp[y][m]=temp[11+y%10][m]-anom_temp;
			}
		}
	}
	else fail("\nNot valid correctionmethod choice\n");

	///////////////////////////////////////////////////////////////////////////////////////////
	// Create mwet ranging CMIP5 period, all years are set to 1961-1990 climatology
	// Not integers.
	for (y = 0; y < NYEAR_HIST; y++) 
	{  
		for (m = 0; m < 12; m++) 
		{
				wet[y][m]=cwet_cru_1961_1990[m];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 
// guess2008 - new functions for reading CRU TS 3.0 binary files.
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU
// Determine temp, precip, sunshine & soilcode
 
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
	double mtemp[NYEAR_CRU][12],double mprec[NYEAR_CRU][12],
	double msun[NYEAR_CRU][12]) {

	// !!!! NEW VERSION OF THIS FUNCTION - guess2008 - NEW VERSION OF THIS FUNCTION !!!!
	// Please note the new function signature. 

	// Archive object. Definition in new header file, cru.h
	Cru_1901_2006Archive ark;

	int target_ilon=(int)(dlon*10.0);
	int target_ilat=(int)(dlat*10.0);

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		soilcode=(int)data.soilcode[0];


		for (y=0;y<NYEAR_CRU;y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[y*12+m]*0.1; // now degC
				mprec[y][m] = data.mprec[y*12+m]*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				msun[y][m]  = data.msun[y*12+m]*0.1;   // % sun 

			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU_MISC
// Determine elevation, frs frq, wet frq & DTR

bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
	double mfrs[NYEAR_HIST][12],double mwet[NYEAR_CRU][12],
	double mdtr[NYEAR_HIST][12]) {
	
	// Please note the new function signature. 

	// Archive object
	Cru_1901_2006miscArchive ark; 
	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006misc data;

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;

		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays.
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000

		for (y=0;y<NYEAR_CRU;y++) { 
			for (m=0;m<12;m++) {

				// guess2008 - catch rounding errors 
				mfrs[y][m] = data.mfrs[y*12+m]*0.01; // days
				if (mfrs[y][m] < 0.1) 
					mfrs[y][m] = 0.0; // Catches rounding errors

				mwet[y][m] = data.mwet[y*12+m]*0.01; // days
				if (mwet[y][m] <= 0.1) 
					mwet[y][m] = 0.0; // Catches rounding errors

				mdtr[y][m] = data.mdtr[y*12+m]*0.1;  // degC

				/*
				If vapour pressure is needed:
				mvap[y][m] = data.mvap[y*12+m]*0.01;
				*/
			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


// guess2008
// Utility function that returns the CRU data from the nearest cell to (lon,lat) within
// a given search radius
bool findnearestCRUdata(int searchradius, char* cruark, double& lon, double& lat, 
                        int& scode, double hist_mtemp1[NYEAR_CRU][12], 
                        double hist_mprec1[NYEAR_CRU][12], 
                        double hist_msun1[NYEAR_CRU][12]) {

	// First try the exact coordinate
	if (searchcru(cruark, lon, lat, scode, hist_mtemp1, hist_mprec1, hist_msun1)) {
		return true;
	}
	
	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon*2)/2;
	double center_lat = floor(lat*2)/2;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.5;

	for (double y = center_lon-searchradius; y <= center_lon+searchradius; y += STEP) {
		for (double x = center_lat-searchradius; x <= center_lat+searchradius; x += STEP) {
			double xdist = x-center_lat;
			double ydist = y-center_lon;
			double dist = sqrt(xdist*xdist + ydist*ydist);
			
			if (dist <= searchradius) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found in CRU
	for (unsigned int i = 0; i < search_points.size(); i++) {
		point search_point = search_points[i].second;
		double search_lon = search_point.first;
		double search_lat = search_point.second;

		if (searchcru(cruark, search_lon, search_lat, scode, 
		              hist_mtemp1, hist_mprec1, hist_msun1)) {
			lon = search_lon;
			lat = search_lat;
			return true;
		}
	}

	return false;
}

/// Help function to define_output_tables, creates one output table
void create_output_table(Table& table, const char* file, const ColumnDescriptors& columns) {
	 table = output_channel->create_table(TableDescriptor(file, columns));
}

/// Defines all output tables
/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void define_output_tables(Pftlist& pftlist) {
	// create a vector with the pft names
	std::vector<std::string> pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		 Pft& pft=pftlist.getobj();

		 pfts.push_back((char*)pft.name);

		 pftlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if (run_landcover) {
		 const char* landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum", "Forest_sum", "Natural_sum", "Peatland_sum"};
		 for (int i=0; i<NLANDCOVERTYPES; i++) {
			  if(run[i]) {
					landcovers.push_back(landcover_string[i]);
			  }
		 }
	}

	// Create the month columns
	ColumnDescriptors month_columns;
	ColumnDescriptors month_columns_wide;
	xtring months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
	for (int i = 0; i < 12; i++) {
		 month_columns      += ColumnDescriptor(months[i], 8,  3);
		 month_columns_wide += ColumnDescriptor(months[i], 10, 3);
	}

	// Create the columns for each output file

	// CMASS
	ColumnDescriptors cmass_columns;
	cmass_columns += ColumnDescriptors(pfts, 8, 3);
	cmass_columns += ColumnDescriptor("Total", 8, 3);
	cmass_columns += ColumnDescriptors(landcovers, 13, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;

	// DENS
	ColumnDescriptors dens_columns;
	dens_columns += ColumnDescriptors(pfts, 8, 4);
	dens_columns += ColumnDescriptor("Total", 8, 4);
	dens_columns += ColumnDescriptors(landcovers, 13, 4);

	// LAI
	ColumnDescriptors lai_columns = dens_columns;

	// CFLUX
	ColumnDescriptors cflux_columns;
	cflux_columns += ColumnDescriptor("Veg",     8, 3);
	cflux_columns += ColumnDescriptor("Soil",    8, 3);
	cflux_columns += ColumnDescriptor("Fire",    8, 3);
	cflux_columns += ColumnDescriptor("Est",     8, 3);
	if (run_landcover) {
		 cflux_columns += ColumnDescriptor("Harvest", 9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",    10, 5);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",   8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",  8, 3);
		cpool_columns += ColumnDescriptor("SoilfC", 8, 3);
		cpool_columns += ColumnDescriptor("SoilsC", 8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LittVC", 8, 3);
		cpool_columns += ColumnDescriptor("LittSC", 8, 3);
		cpool_columns += ColumnDescriptor("CwdC",   8, 3);
		cpool_columns += ColumnDescriptor("MicroC", 8, 3);
		cpool_columns += ColumnDescriptor("HumusC", 8, 3);
		cpool_columns += ColumnDescriptor("SoilC",  8, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC", 10, 3);
	}
	cpool_columns += ColumnDescriptor("Total", 10, 3);

	// FIRERT
	ColumnDescriptors firert_columns;
	firert_columns += ColumnDescriptor("FireRT", 8, 1);

	// RUNOFF
	ColumnDescriptors runoff_columns;
	runoff_columns += ColumnDescriptor("Surf", 8, 1);
	runoff_columns += ColumnDescriptor("Drain", 8, 1);
	runoff_columns += ColumnDescriptor("Base", 8, 1);
	runoff_columns += ColumnDescriptor("Total", 8, 1);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts, 8, 2);

	// AISO
	ColumnDescriptors aiso_columns;
	aiso_columns += ColumnDescriptors(pfts, 10, 3);
	aiso_columns += ColumnDescriptor("Total", 10, 3);
	aiso_columns += ColumnDescriptors(landcovers, 13, 3);

	// AMON
	ColumnDescriptors amon_columns = aiso_columns;



	// GUESSN

	//TODO Fix these for landcover

	// CTON
	ColumnDescriptors cton_columns;
	cton_columns += ColumnDescriptors(pfts, 8, 3);
	cton_columns += ColumnDescriptor("Total", 8, 4);
	cton_columns += ColumnDescriptors(landcovers, 13, 3);

	// NMASS
	ColumnDescriptors nmass_columns;
	nmass_columns += ColumnDescriptors(pfts, 8, 3);
	nmass_columns += ColumnDescriptor("Total", 8, 3);
	nmass_columns += ColumnDescriptors(landcovers, 13, 3);

	// NSOURCES
	ColumnDescriptors nsources_columns;
	nsources_columns += ColumnDescriptor("dep",     8, 3);
	nsources_columns += ColumnDescriptor("fix",     8, 3);
	nsources_columns += ColumnDescriptor("input",   8, 3);
	nsources_columns += ColumnDescriptor("min",     8, 3);
	nsources_columns += ColumnDescriptor("imm",     8, 3);
	nsources_columns += ColumnDescriptor("min-imm", 8, 3);
	nsources_columns += ColumnDescriptor("Total",   8, 3);
	nsources_columns += ColumnDescriptor("Ndemand", 8, 3);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",   8, 4);
	npool_columns += ColumnDescriptor("LittVN", 8, 4);
	npool_columns += ColumnDescriptor("LittSN", 8, 4);
	npool_columns += ColumnDescriptor("CwdN",   8, 4);
	npool_columns += ColumnDescriptor("MicroN", 8, 4);
	npool_columns += ColumnDescriptor("HumusN", 8, 4);
	npool_columns += ColumnDescriptor("SoilN",  8, 4);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN", 8, 4);
	}

	npool_columns += ColumnDescriptor("Total", 10, 4);

	// NLEACH
	ColumnDescriptors nleach_columns;
	nleach_columns += ColumnDescriptor("Min",   8, 3);
	nleach_columns += ColumnDescriptor("Org",   8, 3);
	nleach_columns += ColumnDescriptor("Total", 8, 3);
	nleach_columns += ColumnDescriptor("Temp",  8, 3);
	nleach_columns += ColumnDescriptor("Prec",  8, 3);
	nleach_columns += ColumnDescriptor("Insol", 8, 3);

	// NUPTAKE
	ColumnDescriptors nuptake_columns;
	nuptake_columns += ColumnDescriptors(pfts, 9, 5);
	nuptake_columns += ColumnDescriptor("Total", 9, 5);
	nuptake_columns += ColumnDescriptors(landcovers, 13, 4);

	// ANPPN
	ColumnDescriptors anppn_columns;
	anppn_columns += ColumnDescriptors(pfts, 9, 5);
	anppn_columns += ColumnDescriptor("Total", 8, 3);
	anppn_columns += ColumnDescriptors(landcovers, 13, 3);

	// VMAXNLIM
	ColumnDescriptors vmaxnlim_columns;
	vmaxnlim_columns += ColumnDescriptors(pfts, 8, 3);
	vmaxnlim_columns += ColumnDescriptor("Total", 8, 3);
	vmaxnlim_columns += ColumnDescriptors(landcovers, 13, 4);

	// NLIM
	ColumnDescriptors nlim_columns;
	nlim_columns += ColumnDescriptors(pfts, 8, 3);
	nlim_columns += ColumnDescriptor("Total", 8, 3);
	nlim_columns += ColumnDescriptors(landcovers, 13, 4);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("NH3",   10, 6);
	nflux_columns += ColumnDescriptor("NO",    10, 6);
	nflux_columns += ColumnDescriptor("NO2",   10, 6);
	nflux_columns += ColumnDescriptor("N2O",   10, 6);
	nflux_columns += ColumnDescriptor("Total", 10, 6);
	nflux_columns += ColumnDescriptor("Nconc", 10, 4);

	// CANOPYH
	ColumnDescriptors canopyh_columns;
	canopyh_columns += ColumnDescriptor("Height", 8, 1);

	// ALLOMETRY_IND
	ColumnDescriptors allometry_ind_columns;
	allometry_ind_columns += ColumnDescriptor("PFT",     9, 0);
	allometry_ind_columns += ColumnDescriptor("Age",     9, 0);
	allometry_ind_columns += ColumnDescriptor("N",       9, 0);
	allometry_ind_columns += ColumnDescriptor("Mfol",    9, 3);
	allometry_ind_columns += ColumnDescriptor("Mstem",   9, 2);
	allometry_ind_columns += ColumnDescriptor("Mcroot",  9, 2);
	allometry_ind_columns += ColumnDescriptor("Mfroot",  9, 3);
	allometry_ind_columns += ColumnDescriptor("Mroot",   9, 2);
	allometry_ind_columns += ColumnDescriptor("Mwood",   9, 2);
	allometry_ind_columns += ColumnDescriptor("Mactive", 9, 3);
	allometry_ind_columns += ColumnDescriptor("M",       9, 2);
	allometry_ind_columns += ColumnDescriptor("H",       9, 3);
	allometry_ind_columns += ColumnDescriptor("D",       9, 3);

	// end GUESSN



	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_cmass,          file_cmass,          cmass_columns);
	create_output_table(out_anpp,           file_anpp,           anpp_columns);
	create_output_table(out_dens,           file_dens,           dens_columns);
	create_output_table(out_lai,            file_lai,            lai_columns);
	create_output_table(out_cflux,          file_cflux,          cflux_columns);
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
	create_output_table(out_firert,         file_firert,         firert_columns);
	create_output_table(out_runoff,         file_runoff,         runoff_columns);
	create_output_table(out_speciesheights, file_speciesheights, speciesheights_columns);
	create_output_table(out_aiso,           file_aiso,           aiso_columns);
	create_output_table(out_amon,           file_amon,           amon_columns);

	// GUESSN
	create_output_table(out_cton,           file_cton,           cton_columns);
	create_output_table(out_nmass,          file_nmass,          nmass_columns);
	create_output_table(out_nsources,       file_nsources,       nsources_columns);
	create_output_table(out_npool,          file_npool,          npool_columns);
	create_output_table(out_nleach,         file_nleach,         nleach_columns);
	create_output_table(out_nuptake,        file_nuptake,        nuptake_columns);
	create_output_table(out_anppn,          file_anppn,          anppn_columns);
	create_output_table(out_vmaxnlim,       file_vmaxnlim,       vmaxnlim_columns);
	create_output_table(out_nlim,           file_nlim,           nlim_columns);
	create_output_table(out_nflux,          file_nflux,          nflux_columns);
	create_output_table(out_canopyh,        file_canopyh,        canopyh_columns);
	create_output_table(out_allometry_ind,  file_allometry_ind,  allometry_ind_columns);

	// end GUESSN

	// *** MONTHLY OUTPUT VARIABLES ***

	create_output_table(out_mnpp,         file_mnpp,         month_columns);
	create_output_table(out_mlai,         file_mlai,         month_columns);
	create_output_table(out_mgpp,         file_mgpp,         month_columns);
	create_output_table(out_mra,          file_mra,          month_columns);
	create_output_table(out_maet,         file_maet,         month_columns);
	create_output_table(out_mpet,         file_mpet,         month_columns);
	create_output_table(out_mevap,        file_mevap,        month_columns);
	create_output_table(out_mrunoff,      file_mrunoff,      month_columns_wide);
	create_output_table(out_mintercep,    file_mintercep,    month_columns);
	create_output_table(out_mrh,          file_mrh,          month_columns);
	create_output_table(out_mnee,         file_mnee,         month_columns);
	create_output_table(out_mwcont_upper, file_mwcont_upper, month_columns);
	create_output_table(out_mwcont_lower, file_mwcont_lower, month_columns);
	create_output_table(out_miso,         file_miso,         month_columns_wide);
	create_output_table(out_mmon,         file_mmon,         month_columns_wide);
}

///////////////////////////////////////////////////////////////////////////////////////
// INITIO
// Called by the framework at the start of the model run

void initio(int argc,char* argv[],Pftlist& pftlist) {

	// DESCRIPTION
	// Initialises input/output (e.g. opening files), sets values for the global
	// simulation parameter variables (currently vegmode, npatch, patcharea,
	// ifdailynpp, ifdailydecomp, ifbgestab, ifsme, ifstochestab, ifstochmort, iffire,
	// estinterval, npft), initialises pftlist (the one and only list of PFTs and their
	// static parameters for this run of the model). Normally all of the above
	// parameters, and possibly others, are read from the ins file (see above).
	// Function readins should be called to input settings from the ins file. The
	// syntax for this call should be similar to the following (note that readins
	// returns false in the event of an error in the ins file; normally this should
	// result in program termination):
	//
	// xtring insfilename=argv[1];
	// if (!readins(insfilename,pftlist))
	//     fail("\nUsage: %s <instruction-script-filename> | -help",argv[0]);
	//
	// Arguments argc and argv normally correspond to the command-line arguments
	// imported from the main function (main module, usually main.cpp). The first
	// command line argument (argv[0]) is the name of the binary executable (e.g.
	// guess, guess.exe); the second (argv[1]) should normally be the ins file name.
	// This demonstration version of initio also implements "-help" as an alternative
	// command-line argument, resulting in output of a brief description of the
	// keywords recognised in the ins file, instead of a model run.

	///////////////////////////////////////////////////////////////////////////////////
	// GENERIC SECTION - DO NOT MODIFY

	bool abort;
	xtring insfilename;
	xtring header;


	unixtime(header);
	header=(xtring)"[LPJ-GUESS  "+header+"]\n\n";
	dprintf((char*)header);

	abort=false;
	if (argc>1) {
		insfilename=argv[1];
		if (insfilename[0]=='-') {
			if (insfilename.lower()=="-help") {
				printhelp();
				abort=true;
			}
			else {
				dprintf("Unknown option \"%s\"\n",(char*)insfilename);
				abort=true;
			}
		}
		else if (!fileexists(insfilename)) {
			dprintf("Error: could not open %s for input\n",(char*)insfilename);
			abort=true;
		}

		// Initialise simulation settings and PFT parameters from instruction script
		// Call to readins() returns false if file could not be opened for reading
		// or contained errors (including missing parameters)

		else if (!readins(insfilename,pftlist)){
			abort=true;
		}
	}
	else {
		abort=true;
	}

	if (abort) fail("\nUsage: %s <instruction-script-filename> | -help",argv[0]);

	// Print the title of this run
	dprintf("\n\n-----------------------------------------------\n%s\n-----------------------------------------------\n",(char*)title);

	///////////////////////////////////////////////////////////////////////////////////
	// USER-SPECIFIC SECTION (Modify as necessary or supply own code)
	//
	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon,dlat;
	bool eof=false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;
	
	ngridcell=0;

	while (!eof) {
		
	// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
			ngridcell++;
		}
	}

	fclose(in_grid);

	// Read CO2 data from file
	if (!ifcmip5) {
		co2.load_file(param["file_co2"].str);
	}
	else {
		// Retrieve name of CO2 file from ins file
		path_cmip5_co2=param["path_cmip5_co2"].str;
		xtring filename=path_cmip5_co2;

		if (rcp=="26")
			filename+="co2_1850_2100_hist_rcp26.txt";
		else if (rcp=="45")
			filename+="co2_1850_2100_hist_rcp45.txt"; 
		else if (rcp=="60")
			filename+="co2_1850_2100_hist_rcp60.txt";
		else if (rcp=="85")
			filename+="co2_1850_2100_hist_rcp85.txt";
		else fail("CO2 file not valid");

		co2.load_file(filename);
	}

	//AA CMIP5
	correctionmethod=param["correctionmethod"].str;
	gcm=param["gcm"].str;
	rcp=param["rcp"].str;

	file_cmip5hist=param["path_cmip5hist"].str;
	file_cmip5scen=param["path_cmip5scen"].str;

	if (rcp=="26")
		file_cmip5scen+="mpi_esm_lr_rcp26_r1i1p1/";
	else if (rcp=="45")
		file_cmip5scen+="mpi_esm_lr_rcp45_r1i1p1/";
	else if (rcp=="60")
		file_cmip5scen+="mpi_esm_lr_rcp60_r1i1p1/";
	else if (rcp=="85")
		file_cmip5scen+="mpi_esm_lr_rcp85_r1i1p1/";
	else fail("CMIP5 scenario file not valid");

	file_cmip5hist+="cmip5_hist.bin";
	file_cmip5scen+="cmip5_scen.bin";

	// GUESSN
	file_ndep=param["file_ndep"].str;
	if (file_ndep=="")
		ifndepdata=false;
	else {
		xtring file_ndep_hist=file_ndep+".bin";
		FILE* in_ndep=fopen(file_ndep_hist,"rt");
		if (!in_ndep)
			fail("initio: could not open %s for input",(char*)file_ndep_hist);

		fclose(in_ndep);
		ifndepdata=true;
	}
	// end GUESSN

	// CMIP5
	if (ifcmip5) {
		xtring file_ndep_cmip5;
		if (rcp=="26")
			file_ndep_cmip5=file_ndep+"RCP26.bin";
		else if (rcp=="45")
			file_ndep_cmip5=file_ndep+"RCP45.bin";
		else if (rcp=="60")
			file_ndep_cmip5=file_ndep+"RCP60.bin";
		else if (rcp=="85")
			file_ndep_cmip5=file_ndep+"RCP85.bin";
		else fail("N dep file not valid");

		FILE* in_ndep=fopen(file_ndep_cmip5,"rt");
		if (!in_ndep)
			fail("initio: could not open %s for input",(char*)file_ndep_cmip5);

		fclose(in_ndep);
	}

	if (run_landcover) {
		all_fracs_const=true;	//If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year


		//Retrieve file names for landcover files and open them if static values from ins-file are not used !
		if (!lcfrac_fixed) {	//This version does not support dynamic landcover fraction data

			if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
				file_lu=param["file_lu"].str;
#if defined DYNAMIC_LANDCOVER_INPUT
				if(!LUdata.Open(file_lu))				//Open Bondeau area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_lu);
				else if(LUdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}

			if (run[PEATLAND]) {	//special case for peatland: separate fraction file
				file_peat=param["file_peat"].str;
#if defined DYNAMIC_LANDCOVER_INPUT				
				if(!Peatdata.Open(file_peat))			//Open peatland area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_peat);
				else if(Peatdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}
		}
	}

	// We MUST have an output directory
	if (outputdirectory=="") {
		fail("No output directory given in the .ins file!");
	}

	// Create the output channel
	const int COORDINATES_PRECISION = 1; // decimal places for coords in output
	output_channel = new FileOutputChannel((char*)outputdirectory,
														COORDINATES_PRECISION);

	// Define all output tables and their formats
	define_output_tables(pftlist);

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

///	Loads landcover area fraction data from file(s) for a gridcell.
/** Called from getgridcell() if run_landcover is true. 
  */
bool loadlandcover(Gridcell& gridcell, Coord c)	{
	bool LUerror=false;

	if (!lcfrac_fixed) {
		// Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
		// transferred to gridcell.landcoverfrac each year in getlandcover()

		if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
#if defined DYNAMIC_LANDCOVER_INPUT					
			if (!LUdata.Load(c))		//Load area fraction data from Bondeau input file to data object
			{
				dprintf("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;		// skip this stand
			}
#endif
		}

		if (run[PEATLAND] && !LUerror) {
#if defined DYNAMIC_LANDCOVER_INPUT
			if(!Peatdata.Load(c))	//special case for peatland: separate fraction file
			{
				dprintf("Problems with natural fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand						
			}
#endif
		}
	}

	return LUerror;
}

// GUESSN
/// Retrieves nitrogen deposition for a particular gridcell
/** The values are either taken from the andep parameter in the instruction
 *  file, or from a binary archive file.
 *
 *  The binary archive has nitrogen deposition in gN/m2 on a monthly timestep
 *  for 16 years (Galloway et. al., 2004).
 *
 *  Returned values will not be smaller than minndep.
 *
 *  \param  filename    The file name of the binary archive
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 */
bool getndep(xtring filename,double lon,double lat,Climate& climate) {

	int y,m;
	double dval;
	double dailyndep=2000.0/(4.0*365.0);	// pre-industrial N depostion [gN ha-1] (2 kgN/ha/year)
	double convert=0.0000001;				// converting from gN ha-1 to kgN m-2
	double NHxWetDep_10[26][12]={0.0};
	double NHxDryDep_10[26][12]={0.0};
	double NOyWetDep_10[26][12]={0.0};
	double NOyDryDep_10[26][12]={0.0};

	if (!ifndepdata) {
		for (y=0;y<16;y++) {
			for (m=0;m<12;m++) {
				NHxDryDep_10[y][m]=dailyndep;	
				NHxWetDep_10[y][m]=dailyndep;	
				NOyDryDep_10[y][m]=dailyndep;	
				NOyWetDep_10[y][m]=dailyndep;	
			}
		}
	}

	xtring historic_filename=filename+".bin";

	GlobalNitrogenDepositionArchive ark;
	if (!ark.open(historic_filename)) {
		 fail("Could not open %s for input",(char*)historic_filename);
		 return false;
	}

	GlobalNitrogenDeposition rec;
	rec.longitude = lon;
	rec.latitude = lat;
	
	if (!ark.getindex(rec)) {
		 // The coordinate wasn't found in the archive
		 ark.close();
		 return false;
	}
	else {
		 // Found the record, get the values
		for (y=0;y<16;y++) {
			for (m=0;m<12;m++) {
				NHxDryDep_10[y][m]=rec.NHxDry[y*12+m];
				NHxWetDep_10[y][m]=rec.NHxWet[y*12+m];	
				NOyDryDep_10[y][m]=rec.NOyDry[y*12+m];	
				NOyWetDep_10[y][m]=rec.NOyWet[y*12+m];
			}
		}

		ark.close();
	}

	if (ifcmip5) {

		if (!ifndepdata) {
			for (y=16;y<26;y++) {
				for (m=0;m<12;m++) {
					NHxDryDep_10[y][m]=dailyndep;	// use pre-industrial N depostion 2 kgN/ha/year
					NHxWetDep_10[y][m]=dailyndep;	// use pre-industrial N depostion 2 kgN/ha/year
					NOyDryDep_10[y][m]=dailyndep;	// use pre-industrial N depostion 2 kgN/ha/year
					NOyWetDep_10[y][m]=dailyndep;	// use pre-industrial N depostion 2 kgN/ha/year
				}
			}
		}
		else {

			if (rcp=="26") {

				xtring scenario_filename=filename+"RCP26.bin";

				GlobalNitrogenDepositionRCP26Archive ark_sce;
				if (!ark_sce.open(scenario_filename)) {
					 fail("Could not open %s for input",(char*)scenario_filename);
					 return false;
				}

				GlobalNitrogenDepositionRCP26 rec_sce;

				rec_sce.longitude = lon;
				rec_sce.latitude = lat;

				if (!ark_sce.getindex(rec_sce)) {
					 // The coordinate wasn't found in the archive
					 ark_sce.close();
					 return false;
				}
				else {
					// Found the record, get the values
					for (y=15;y<26;y++) {
						for (m=0;m<12;m++) {
							if (y==15) { // Scenario and hist data has the same year -> avr
								dval=(NHxDryDep_10[y][m]+rec_sce.NHxDry[(y-15)*12+m])/2.0;
								NHxDryDep_10[y][m]=dval;
								
								dval=(NHxWetDep_10[y][m]+rec_sce.NHxWet[(y-15)*12+m])/2.0;
								NHxWetDep_10[y][m]=dval;
								
								dval=(NOyDryDep_10[y][m]+rec_sce.NOyDry[(y-15)*12+m])/2.0;
								NOyDryDep_10[y][m]=dval;
								
								dval=(NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m])/2.0;
								NOyWetDep_10[y][m]=dval;
							}
							else {
								NHxDryDep_10[y][m]=rec_sce.NHxDry[(y-15)*12+m];
								NHxWetDep_10[y][m]=rec_sce.NHxWet[(y-15)*12+m];	
								NOyDryDep_10[y][m]=rec_sce.NOyDry[(y-15)*12+m];	
								NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m];
							}
						}
					}
					ark_sce.close();
				}
			}
			else if (rcp=="45") {

				xtring scenario_filename=filename+"RCP45.bin";

				GlobalNitrogenDepositionRCP45Archive ark_sce;
				if (!ark_sce.open(scenario_filename)) {
					 fail("Could not open %s for input",(char*)scenario_filename);
					 return false;
				}

				GlobalNitrogenDepositionRCP45 rec_sce;

				rec_sce.longitude = lon;
				rec_sce.latitude = lat;

				if (!ark_sce.getindex(rec_sce)) {
					 // The coordinate wasn't found in the archive
					 ark_sce.close();
					 return false;
				}
				else {
					// Found the record, get the values
					for (y=15;y<26;y++) {
						for (m=0;m<12;m++) {
							if (y==15) { // Scenario and hist data has the same year -> avr
								dval=(NHxDryDep_10[y][m]+rec_sce.NHxDry[(y-15)*12+m])/2.0;
								NHxDryDep_10[y][m]=dval;
								
								dval=(NHxWetDep_10[y][m]+rec_sce.NHxWet[(y-15)*12+m])/2.0;
								NHxWetDep_10[y][m]=dval;
								
								dval=(NOyDryDep_10[y][m]+rec_sce.NOyDry[(y-15)*12+m])/2.0;
								NOyDryDep_10[y][m]=dval;
								
								dval=(NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m])/2.0;
								NOyWetDep_10[y][m]=dval;
							}
							else {
								NHxDryDep_10[y][m]=rec_sce.NHxDry[(y-15)*12+m];
								NHxWetDep_10[y][m]=rec_sce.NHxWet[(y-15)*12+m];	
								NOyDryDep_10[y][m]=rec_sce.NOyDry[(y-15)*12+m];	
								NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m];
							}
						}
					}
					ark_sce.close();
				}
			}
			else if (rcp=="60") {

				xtring scenario_filename=filename+"RCP60.bin";

				GlobalNitrogenDepositionRCP60Archive ark_sce;
				if (!ark_sce.open(scenario_filename)) {
					 fail("Could not open %s for input",(char*)scenario_filename);
					 return false;
				}

				GlobalNitrogenDepositionRCP60 rec_sce;

				rec_sce.longitude = lon;
				rec_sce.latitude = lat;

				if (!ark_sce.getindex(rec_sce)) {
					 // The coordinate wasn't found in the archive
					 ark_sce.close();
					 return false;
				}
				else {
					// Found the record, get the values
					for (y=15;y<26;y++) {
						for (m=0;m<12;m++) {
							if (y==15) { // Scenario and hist data has the same year -> avr
								dval=(NHxDryDep_10[y][m]+rec_sce.NHxDry[(y-15)*12+m])/2.0;
								NHxDryDep_10[y][m]=dval;
								
								dval=(NHxWetDep_10[y][m]+rec_sce.NHxWet[(y-15)*12+m])/2.0;
								NHxWetDep_10[y][m]=dval;
								
								dval=(NOyDryDep_10[y][m]+rec_sce.NOyDry[(y-15)*12+m])/2.0;
								NOyDryDep_10[y][m]=dval;
								
								dval=(NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m])/2.0;
								NOyWetDep_10[y][m]=dval;
							}
							else {
								NHxDryDep_10[y][m]=rec_sce.NHxDry[(y-15)*12+m];
								NHxWetDep_10[y][m]=rec_sce.NHxWet[(y-15)*12+m];	
								NOyDryDep_10[y][m]=rec_sce.NOyDry[(y-15)*12+m];	
								NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m];
							}
						}
					}
					ark_sce.close();
				}
			}
			else if (rcp=="85") {

				xtring scenario_filename=filename+"RCP85.bin";

				GlobalNitrogenDepositionRCP85Archive ark_sce;
				if (!ark_sce.open(scenario_filename)) {
					 fail("Could not open %s for input",(char*)scenario_filename);
					 return false;
				}

				GlobalNitrogenDepositionRCP85 rec_sce;

				rec_sce.longitude = lon;
				rec_sce.latitude = lat;

				if (!ark_sce.getindex(rec_sce)) {
					 // The coordinate wasn't found in the archive
					 ark_sce.close();
					 return false;
				}
				else {
					// Found the record, get the values
					for (y=15;y<26;y++) {
						for (m=0;m<12;m++) {
							if (y==15) { // Scenario and hist data has the same year -> avr
								dval=(NHxDryDep_10[y][m]+rec_sce.NHxDry[(y-15)*12+m])/2.0;
								NHxDryDep_10[y][m]=dval;
								
								dval=(NHxWetDep_10[y][m]+rec_sce.NHxWet[(y-15)*12+m])/2.0;
								NHxWetDep_10[y][m]=dval;
								
								dval=(NOyDryDep_10[y][m]+rec_sce.NOyDry[(y-15)*12+m])/2.0;
								NOyDryDep_10[y][m]=dval;
								
								dval=(NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m])/2.0;
								NOyWetDep_10[y][m]=dval;
							}
							else {
								NHxDryDep_10[y][m]=rec_sce.NHxDry[(y-15)*12+m];
								NHxWetDep_10[y][m]=rec_sce.NHxWet[(y-15)*12+m];	
								NOyDryDep_10[y][m]=rec_sce.NOyDry[(y-15)*12+m];	
								NOyWetDep_10[y][m]=rec_sce.NOyWet[(y-15)*12+m];
							}
						}
					}
					ark_sce.close();
				}
			}	
		}
	}

	// interpolate to all hist and scenario years

	int years[]={5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255};
	int interyear[2]={0.0};
	int yy=0;

	for (y=0;y<NYEAR_HIST;y++) {

		bool found=false;
		while (!found){
			if (y<=years[0]){
				interyear[0]=0;
				interyear[1]=0;
				found=true;
			}
			else if (y<=years[yy]){
				interyear[0]=yy-1;
				interyear[1]=yy;
				found=true;
			}
			else
				yy++;
		}

		for (m=0;m<12;m++){

			NHxWetDep[y][m]=(NHxWetDep_10[interyear[0]][m]+((double)(y-years[interyear[0]]))/10.0*
									(NHxWetDep_10[interyear[1]][m]-NHxWetDep_10[interyear[0]][m]))*convert;
			NHxDryDep[y][m]=(NHxDryDep_10[interyear[0]][m]+((double)(y-years[interyear[0]]))/10.0*
									(NHxDryDep_10[interyear[1]][m]-NHxDryDep_10[interyear[0]][m]))*convert;
			NOyWetDep[y][m]=(NOyWetDep_10[interyear[0]][m]+((double)(y-years[interyear[0]]))/10.0*
									(NOyWetDep_10[interyear[1]][m]-NOyWetDep_10[interyear[0]][m]))*convert;
			NOyDryDep[y][m]=(NOyDryDep_10[interyear[0]][m]+((double)(y-years[interyear[0]]))/10.0*
									(NOyDryDep_10[interyear[1]][m]-NOyDryDep_10[interyear[0]][m]))*convert;
		}
	}
	return true;
}
// end GUESSN

///////////////////////////////////////////////////////////////////////////////////////
// GETGRIDCELL
// Called by the framework at the start of the simulation for a particular grid cell
bool getgridcell(Gridcell& gridcell) {

	// DESCRIPTION
	// Obtains coordinates and soil static parameters for the next grid cell to
	// simulate. The function should return false if no grid cells remain to be simulated,
	// otherwise true. Currently the following member variables of Gridcell should be
	// initialised: longitude, latitude and climate.instype; the following members of
	// member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	// thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	// code (Sitch et al 2000) by a call to function soilparameters in the driver
	// module (driver.cpp):
	//
	// soilparameters(gridcell.soiltype,soilcode);
	//
	// If the model is to be driven by quasi-daily values of the climate variables
	// derived from monthly means, this function may be the appropriate place to
	// perform the required interpolations. The utility functions interp_monthly_means
	// and interp_monthly_totals in driver.cpp may be called for this purpose.

	// Select coordinates for next grid cell in linked list
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	bool gridfound;
	bool LUerror=false;

	// to ensure an identical random number sequence for each gridcell.
	setseed(12345678);

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	static bool first_call = true;

	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		
		// guess2008 - New searchcru functions takee the CRU filenames as their first 
		// argument, i.e. cru_1901_2002.bin and cru_1901_2002_misc.bin

		// New code:

		double lon = gridlist.getobj().lon;
		double lat = gridlist.getobj().lat;

		// GUESSN
		if (!getndep(file_ndep,lon,lat,gridcell.climate)) {
			dprintf("Year %d 13\n",date.year);
			fail("Grid cell Lat %g Long %g not found in %s",lat,lon,(char*)file_ndep);
		}
		// end GUESSN

		if (!ifcmip5) {
			gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
			                               hist_mtemp, hist_mprec, hist_msun);

			if (gridfound) // Get more historical CRU data for this grid cell
				gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);

			if (run_landcover) {
				Coord& c=gridlist.getobj();
				LUerror=loadlandcover(gridcell, c);
			}
			if (LUerror)
				gridfound=false;

			while (!gridfound) {

				if (run_landcover && LUerror)
					dprintf("\nError: could not find stand at (%g,%g) in landcover data file\n", gridlist.getobj().lon,gridlist.getobj().lat);
				else
					dprintf("\nError: could not find stand at (%g,%g) in CRU data file\n", gridlist.getobj().lon,gridlist.getobj().lat);

				gridlist.nextobj();
				if (gridlist.isobj) {
					double lon = gridlist.getobj().lon;
					double lat = gridlist.getobj().lat;
					gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
					                               hist_mtemp, hist_mprec, hist_msun);
			  
					if (gridfound) // Get more historical CRU data for this grid cell
						gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation,
						                           hist_mfrs, hist_mwet, hist_mdtr);

					if (run_landcover) {
						Coord& c=gridlist.getobj();
						LUerror=loadlandcover(gridcell, c);
					}
					if (LUerror)
						gridfound=false;
				}
				else return false;
				}
		}
		else {	// CMIP5

			gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
			                               mtemp_cru, mprec_cru, msun_cru);

			if (gridfound) // Get more historical CRU data for this grid cell
				gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
				                           hist_mfrs, mwet_cru, hist_mdtr);

			while (!gridfound) {

				dprintf("\nError: could not find stand at (%g,%g) in CRU data file\n",
					gridlist.getobj().lon,gridlist.getobj().lat);

				gridlist.nextobj();
				if (gridlist.isobj) {
					double lon = gridlist.getobj().lon;
					double lat = gridlist.getobj().lat;
					gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
					                               mtemp_cru, mprec_cru, msun_cru);
			  
					if (gridfound) // Get more historical CRU data for this grid cell
						gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation,
						                           hist_mfrs, mwet_cru, hist_mdtr);

				}
				else return false;
			}

			// CMIP5 - land use input
			if(iflandusesimple)
			  gridfound=searchlanduse(lon,lat,
						  hist_frluse);

			// CMIP5 AA - cmip 5 historical input
			gridfound=searchcmip5hist(file_cmip5hist, lon, lat, mtemp_cmip5, mprec_cmip5, mswrad_cmip5);
		
			// CMIP5 AA - cmip5 scenario input
			gridfound=searchcmip5scen(file_cmip5scen,lon,lat,mtemp_cmip5,mprec_cmip5,mswrad_cmip5);

			// CMIP5 AA - calculate swrad from CRU cloudiness
			calculate_swrad(msun_cru, gridlist.getobj().lon, gridlist.getobj().lat, mswrad_cru);

			createclimatology_cru(mtemp_cru,mprec_cru,msun_cru, mwet_cru,mswrad_cru, 
				clim_mtemp_cru, clim_mprec_cru,clim_msun_cru, clim_swrad_cru, clim_mwet_cru_1901_1930, clim_mwet_cru_1961_1990); 

			createclimatology_cmip5(mtemp_cmip5, mprec_cmip5, mswrad_cmip5, clim_mtemp_cmip5, clim_mprec_cmip5, clim_swrad_cmip5);

			makeCMIP5data(mtemp_cmip5,mprec_cmip5,mswrad_cmip5, 
				NHxWetDep,NHxDryDep,NOyWetDep,NOyDryDep,
				clim_mtemp_cmip5, clim_mprec_cmip5, clim_swrad_cmip5,
				clim_mtemp_cru,clim_mprec_cru,clim_msun_cru,clim_swrad_cru, clim_mwet_cru_1901_1930, clim_mwet_cru_1961_1990,
				mtemp_cru,mprec_cru,msun_cru, mwet_cru, mswrad_cru,
				hist_mtemp,hist_mprec,hist_msun, hist_mwet);
			// hist_* is now the new data created of CRU and CMIP5
		}		   

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);
		spinup_mdtr.detrend_data();

		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);
		
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		if (!ifcmip5)
			gridcell.climate.instype=SUNSHINE;
		else {
			gridcell.climate.instype=SWRAD_TS;
			if (correctionmethod=="c7") gridcell.climate.instype=SUNSHINE; //AA CMIP5
		}

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

///////////////////////////////////////////////////////////////////////////////////////
// GETLANDCOVER
// Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
void getlandcover(Gridcell& gridcell,Pftlist& pftlist) {
	int i, year;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;

	if(date.year<nyear_spinup)					//Use values for first historic year during spinup period !
		year=0;
	else if(date.year>=nyear_spinup+NYEAR_LU)	//AR4 adaptation
	{
//		dprintf("setting LU data for scenario period\n");
		year=NYEAR_LU-1;
	}
	else
		year=date.year-nyear_spinup;

	if(lcfrac_fixed)	// If area fractions are set in the ins-file.
	{
		if(date.year==0) // called by landcover_init
		{
			int nactive_landcovertypes=0;

			if(equal_landcover_area)
			{
				for(i=0;i<NLANDCOVERTYPES;i++)
				{
					if(run[i])
						nactive_landcovertypes++;
				}
			}

			for(i=0;i<NLANDCOVERTYPES;i++)
			{
				if(equal_landcover_area)
				{
					sum_active+=gridcell.landcoverfrac[i]=1.0*run[i]/(double)nactive_landcovertypes;	// only set fractions that are active !
					sum_tot=sum_active;
				}
				else
				{
					sum_tot+=gridcell.landcoverfrac[i]=(double)lc_fixed_frac[i]/100.0;					//count sum of all fractions (should be 1.0)

					if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)					//discard unreasonable values
					{
						if(date.year==0)
							dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
						sum_tot-=gridcell.landcoverfrac[i];
						gridcell.landcoverfrac[i]=0.0;
					}

					sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];				//only set fractions that are active !
				}
			}
			
			if(sum_tot<0.99 || sum_tot>1.01)	// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions
				if(date.year==0)
					dprintf("WARNING ! landcover fixed fraction sum is %4.2f, rescaling landcover fractions !\n", sum_tot);

				for(i=0;i<NLANDCOVERTYPES;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}

			//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
			if(sum_active<0.99)	//if landcover types are turned off in the ini-file, always <=1.0 here
			{
				if(date.year==0)
					dprintf("WARNING ! landcover active fraction sum is %4.2f.\n", sum_active);

				if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
				{
					if(date.year==0)
						dprintf("Inactive fractions (%4.2f) transferred to NATURAL fraction.\n", 1.0-sum_active);

					gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference 1.0-(sum of active landcover fractions) are added to the natural fraction
				}
				else
				{
/*					if(date.year==0)
						dprintf("Rescaling landcover fractions !\n");
					for(i=0;i<NLANDCOVERTYPES;i++)
						gridcell.landcoverfrac[i]/=sum_active;			// if NATURAL not simulated, rescale active fractions to 1.0
*/					if(date.year==0)
						dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
				}
																	
			}
		}
	}
	else	//area fractions are read from input file(s);
	{
		if(run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST])
		{	

			for(i=0;i<PEATLAND;i++)		//peatland fraction data is not in this file, otherwise i<NLANDCOVERTYPES.
			{	
#if defined DYNAMIC_LANDCOVER_INPUT
				sum_tot+=gridcell.landcoverfrac[i]=LUdata.Get(year,i);					//count sum of all fractions (should be 1.0)
#endif
				if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)			//discard unreasonable values
				{		
					if(date.year==0)
						dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
					sum_tot-=gridcell.landcoverfrac[i];
					gridcell.landcoverfrac[i]=0.0;
				}

				sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];
			}

			if(sum_tot!=1.0)		// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions

				if(sum_tot<0.99 || sum_tot>1.01)
				{
					if(date.year==0)
					{
						dprintf("WARNING ! landcover fraction sum is %4.2f for year %d\n", sum_tot, year+FIRSTHISTYEAR);
						dprintf("Rescaling landcover fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);
					}
				}
				else				//added scaling to sum=1.0 (sum often !=1.0)
					dprintf("Rescaling landcover fractions year %d ! (sum is within 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);

				for(i=0;i<PEATLAND;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}
		}
		else
			gridcell.landcoverfrac[NATURAL]=0.0;

		if(run[PEATLAND])
		{
#if defined DYNAMIC_LANDCOVER_INPUT
			sum_active+=gridcell.landcoverfrac[PEATLAND]=Peatdata.Get(year,"PEATLAND");			//peatland fraction data is currently in a separate file !
#endif
		}

		//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
		if(sum_active!=1.0)		//if landcover types are turned off in the ini-file, or if more landcover types are added in other input files, can be either less or more than 1.0
		{
			if(date.year==0)
				dprintf("Landcover fraction sum not 1.0 !\n");

			if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
			{
				if(date.year==0)
				{
					if(sum_active<1.0)
						dprintf("Inactive fractions (%4.3f) transferred to NATURAL fraction.\n", 1.0-sum_active);
					else
						dprintf("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).\n", sum_active-1.0, gridcell.landcoverfrac[NATURAL]);
				}

				gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					dprintf("New NATURAL fraction is %4.3f.\n", gridcell.landcoverfrac[NATURAL]);

				sum_active=1.0;		//sum_active should now be 1.0

				if(gridcell.landcoverfrac[NATURAL]<0.0)	//If new landcover type fraction is bigger than the natural fraction (something wrong in the distribution of input file area fractions)
				{										
					if(date.year==0)
						dprintf("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.\n");

					sum_active-=gridcell.landcoverfrac[NATURAL];	//fraction not possible to transfer moved back to sum_active, which will now be >1.0 again
					gridcell.landcoverfrac[NATURAL]=0.0;

					for(i=0;i<NLANDCOVERTYPES;i++)
					{
						gridcell.landcoverfrac[i]/=sum_active;		//fraction rescaled to unity sum
						if(run[i])
							if(date.year==0)
								dprintf("Landcover type %d fraction is %4.3f\n", i, gridcell.landcoverfrac[i]);
					}
				}
			}
			else
			{
//				if(date.year==0)
//					dprintf("Rescaling landcover fractions !\n");
//				for(i=0;i<NLANDCOVERTYPES;i++)
//					gridcell.landcoverfrac[i]/=sum_active;						// if NATURAL not simulated, rescale active fractions to 1.0
				if(date.year==0)
					dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// GETCLIMATE
// Called by the framework each simulation day before any process modelling is performed for this day
// Obtains climate data (including atmospheric CO2 and insolation) for this day. 
bool getclimate(Gridcell& gridcell) {

	// DESCRIPTION
	// The function should returns false if the simulation is complete for this grid cell,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of gridcell must be
	// initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	// values of the climate variables derived from monthly means, this day's values
	// will presumably be extracted from arrays containing the interpolated daily
	// values (see function getgridcell):
	//
	// gridcell.climate.temp=dtemp[date.day];
	// gridcell.climate.prec=dprec[date.day];
	// gridcell.climate.insol=dsun[date.day];
	// 
	// Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
	// BVOC:
	// gridcell.climate.dtr=ddtr[date.day];
	//
	// If model is run in diurnal mode, which requires appropriate climate forcing data, 
	// additional members of the climate must be initialised: temps, insols. Both of the
	// variables must be of type std::vector. The length of these vectors should be equal
	// to value of date.subdaily which also needs to be set either in getclimate or 
	// getgridcell functions. date.subdaily is a number of sub-daily period in a single 
	// day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
	// diurnal mode.

	double progress;

	// guess2008 - changed name from mwet to mwet_all
	double mwet_all[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // number of rain days per month
	int dd;
	Climate& climate=gridcell.climate;

	if (date.day==0) {

		// First day of year ..
		
		if (date.year<nyear_spinup) {

			// During spinup period

			int m;
			double mtemp[12],mprec[12],msun[12];
			double mfrs[12],mwet[12],mdtr[12];

			for (m=0;m<12;m++) {
				mtemp[m]=spinup_mtemp[m];
				mprec[m]=spinup_mprec[m];
				msun[m]=spinup_msun[m];

				// guess2008
				mfrs[m]=spinup_mfrs[m];
				mwet[m]=spinup_mwet[m];
				mdtr[m]=spinup_mdtr[m];
			}

			// Interpolate monthly spinup data to quasi-daily values
			interp_climate(mtemp,mprec,msun,mdtr,dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(mprec,dprec,mwet);
			}

			spinup_mtemp.nextyear();
			spinup_mprec.nextyear();
			spinup_msun.nextyear();

			// guess2008
			spinup_mfrs.nextyear();
			spinup_mwet.nextyear();
			spinup_mdtr.nextyear();

		}
		else if (date.year<nyear_spinup+NYEAR_HIST) {

			// Historical period

			// Interpolate this year's monthly data to quasi-daily values
			interp_climate(hist_mtemp[date.year-nyear_spinup],
				hist_mprec[date.year-nyear_spinup],hist_msun[date.year-nyear_spinup],
					   hist_mdtr[date.year-nyear_spinup],
				       dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year-nyear_spinup],dprec,hist_mwet[date.year-nyear_spinup]);
			}
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}

		if (ifcmip5) {
			// CMIP5 - land use input
			if (date.year<nyear_spinup){ 
				if(iflandusesimple)
				  climate.frluse=hist_frluse[0];
			}
			else if (date.year<nyear_spinup+NYEAR_HIST){
				if(iflandusesimple)
				  climate.frluse=hist_frluse[date.year-nyear_spinup];
			}
		}

		climate.andep=0.0;
		int m;
		if (date.year<nyear_spinup){
			dd=0;
			for (m=0;m<12;m++) {
				climate.andep+=(NHxDryDep[0][m]+NOyDryDep[0][m]+
					NHxWetDep[0][m]+NOyWetDep[0][m])*date.ndaymonth[m];

				for (int dm=0;dm<date.ndaymonth[m];dm++) {
					climate.dndep[dd]=(NHxDryDep[0][m]+
						NOyDryDep[0][m]+NHxWetDep[0][m]+
						NOyWetDep[0][m]);
					dd++;
				}
			}
		}
		else {
			dd=0;
			for (m=0;m<12;m++) {
				climate.andep+=(NHxDryDep[date.year-nyear_spinup][m]+
					NOyDryDep[date.year-nyear_spinup][m]+
					NHxWetDep[date.year-nyear_spinup][m]+
					NOyWetDep[date.year-nyear_spinup][m])*date.ndaymonth[m];

				for (int dm=0;dm<date.ndaymonth[m];dm++) {
					climate.dndep[dd]=(NHxDryDep[date.year-nyear_spinup][m]+
						NOyDryDep[date.year-nyear_spinup][m]+
						NHxWetDep[date.year-nyear_spinup][m]+
						NOyWetDep[date.year-nyear_spinup][m]);
					dd++;
				}
			}
		}
	}

	// Send environmental values for today to framework

	climate.co2 = co2[FIRSTHISTYEAR + date.year - nyear_spinup];

	// FACE
	//if (date.year > nyear_spinup+NYEAR_HIST-10)
	//	climate.co2=550.0;

	climate.temp=dtemp[date.day];
	climate.prec=dprec[date.day];
	climate.insol=dsun[date.day];

	// bvoc
	if(ifbvoc){
	  climate.dtr=ddtr[date.day];
	}

	// First day of year only ...

	if (date.day==0) {

		if (ifcmip5) {
			if ((correctionmethod=="c5" || correctionmethod=="c6") && date.year==nyear_spinup+NYEAR_HIST-15) return false; //AA CMIP5 break at 2085

			if (correctionmethod=="c7" && date.year==nyear_spinup+NYEAR_CMIP5_HIST) return false; //AA CMIP5 
		}

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+NYEAR_HIST)
				+date.year)/(double)(ngridcell*(nyear_spinup+NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}
	
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////
// Canopy Height
//
double canopy_height(Stand& stand) {

	// DESCRIPTION
	// Determines canopy height from the top 100
	// tallest trees in each patch (Nakai 2010). 

	double accumulated_canopy_height=0.0;

	for (int p=0;p<npatch;p++) {

		double tree_height[10];
		double dens[10];

		for (int k=0;k<10;k++) {
			tree_height[k]=0.0;	//indiv.height
			dens[k]=0.0;		//indiv.densindiv
		}

		Patch& patch=stand[p];
		Vegetation& vegetation=patch.vegetation;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			if (indiv.height > tree_height[9]) {
				for (int i=0;i<10;i++){
					if (indiv.height > tree_height[i]) {
						for (int j=8;j>=i;j--) {
							tree_height[j+1]=tree_height[j];
							dens[j+1]=dens[j];
						}
						tree_height[i]=indiv.height;
						dens[i]=indiv.densindiv;

						i=10;
					}
				}
			}

			vegetation.nextobj();
		}

		double tree_height_patch=0.0;
		double density_patch=0.0;
		int l=0;

		while (density_patch<0.1 && l<10){

			if (density_patch+dens[l]>0.1)
				dens[l]=0.1-density_patch;

			tree_height_patch+=tree_height[l]*dens[l];
			density_patch+=dens[l];
			l++;
		}

		if (!negligible(density_patch))
			accumulated_canopy_height+=tree_height_patch/density_patch;
	}
	return accumulated_canopy_height/(double)npatch;
}


///////////////////////////////////////////////////////////////////////////////////////
// OUTANNUAL
// Called by the framework at the end of the last day of each simulation year
void outannual(Gridcell& gridcell,Pftlist& pftlist) {

	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int c, m, nclass;
	double flux_veg, flux_soil, flux_fire, flux_est, flux_harvest;
	double c_litter, c_fast, c_slow, c_harv_slow; 

	// GUESSN
	double surfsoillitterc,surfsoillittern,cwdc,cwdn,microc,micron,humusc,humusn,centuryc,centuryn,n_harv_slow;
	double flux_nh3,flux_no,flux_no2,flux_n2o,flux_ntot,flux_nconc;
	// end GUESSN

	// guess2008 - hold the monthly average across patches
	double mnpp[12];
	double mgpp[12];
	double mlai[12];
	double maet[12];
	double mpet[12];
	double mevap[12];
	double mintercep[12];
	double mrunoff[12];
	double mrh[12];
	double mra[12];
	double mnee[12];
	double mwcont_upper[12];
	double mwcont_lower[12];
	// bvoc
	double miso[12];
	double mmon[12];

	if (vegmode==COHORT)
		nclass=min(date.year/estinterval+1,OUTPUT_MAXAGECLASS);

	// guess2008 - yearly output after spinup

	if (date.year>=nyear_spinup) {

		double lon = gridcell.get_lon();
		double lat = gridcell.get_lat();

		// The OutputRows object manages the next row of output for each
		// output table
		OutputRows out(output_channel, lon, lat, date.year);

		// guess2008 - reset monthly average across patches each year
		for (m=0;m<12;m++)
			mnpp[m]=mlai[m]=mgpp[m]=mra[m]=maet[m]=mpet[m]=mevap[m]=mintercep[m]=mrunoff[m]=mrh[m]=mnee[m]=mwcont_upper[m]=mwcont_lower[m]=miso[m]=mmon[m]=0.0;

		double landcover_cmass[NLANDCOVERTYPES]={0.0};
		double landcover_anpp[NLANDCOVERTYPES]={0.0};
		double landcover_lai[NLANDCOVERTYPES]={0.0};
		double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};
		double landcover_aiso[NLANDCOVERTYPES]={0.0};
		double landcover_amon[NLANDCOVERTYPES]={0.0};

		double gcpft_cmass=0.0;
		double gcpft_anpp=0.0;
		double gcpft_lai=0.0;
		double gcpft_densindiv_total=0.0;
		double gcpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double gcpft_aiso=0.0;
		double gcpft_amon=0.0;

		double cmass_gridcell=0.0;
		double anpp_gridcell=0.0;
		double lai_gridcell=0.0;
		double surfrunoff_gridcell=0.0;
		double drainrunoff_gridcell=0.0;
		double baserunoff_gridcell=0.0;
		double runoff_gridcell=0.0;
		double dens_gridcell=0.0;
		double firert_gridcell=0.0;
		double aiso_gridcell=0.0;
		double amon_gridcell=0.0;

		double canopyheight_gridcell=0.0;

		double standpft_cmass=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double standpft_aiso=0.0;
		double standpft_amon=0.0;

		// GUESSN

		double landcover_cmass_leaf[NLANDCOVERTYPES]={0.0};
		double landcover_nmass_leaf[NLANDCOVERTYPES]={0.0};
		double landcover_nmass[NLANDCOVERTYPES]={0.0};
		double landcover_nuptake[NLANDCOVERTYPES]={0.0};
		double landcover_anppn[NLANDCOVERTYPES]={0.0};
		double landcover_vmaxnlim[NLANDCOVERTYPES]={0.0};
		double landcover_nlim[NLANDCOVERTYPES]={0.0};

		double gcpft_cmass_leaf=0.0;	// nitrogen mass of leafs
		double gcpft_nmass_leaf=0.0;	// carbon mass of leafs
		double gcpft_nmass=0.0;			// sum/mean across patches for nitrogen biomass (kgN/m2)
		double gcpft_nuptake=0.0;		// sum across patches for nitrogen uptake (kgN/m2)
		double gcpft_anppn=0.0;			// sum across patches for nitrogen ANPP usage (kgN/m2)
		double gcpft_vmaxnlim=0.0;		// N limitation on vm
		double gcpft_nlim=0.0;			// N limitation on growth

		double cmass_leaf_gridcell=0.0;
		double nmass_leaf_gridcell=0.0;
		double nmass_gridcell=0.0;
		double nuptake_gridcell=0.0;
		double anppn_gridcell=0.0;
		double vmaxnlim_gridcell=0.0;
		double nlim_gridcell=0.0;
		double anpp_no_nlim_gridcell=0.0;

		double standpft_cmass_leaf=0.0;
		double standpft_nmass_leaf=0.0;
		double standpft_nmass=0.0;
		double standpft_nuptake=0.0;
		double standpft_anppn=0.0;
		double standpft_vmaxnlim=0.0;
		double standpft_nlim=0.0;
		double standpft_anpp_no_nlim=0.0;

		surfsoillitterc,surfsoillittern,cwdc,cwdn,microc,micron,humusc,humusn,centuryc,centuryn=0.0;
		double n_litter=0.0;
		double andep_gridcell=0.0;
		double anmin_gridcell=0.0;
		double animm_gridcell=0.0;
		double anfix_gridcell=0.0;
		double nsupply_gridcell=0.0;
		double ndemand_gridcell=0.0;
		double n_min_leach_gridcell=0.0;
		double n_org_leach_gridcell=0.0;
		// end GUESSN

		// *** Loop through PFTs ***

		double nmass_indiv=0.0;			// GUESSN N budget
		double nmass_sompools=0.0;		// GUESSN N budget
		double nmass_litterpools=0.0;	// GUESSN N budget

		pftlist.firstobj();
		while (pftlist.isobj) {

			Pft& pft=pftlist.getobj();
			Gridcellpft& gridcellpft=gridcell.pft[pft.id];

			// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs		
			gcpft_cmass=0.0;
			gcpft_anpp=0.0;
			gcpft_lai=0.0;
			gcpft_densindiv_total=0.0;
			gcpft_aiso=0.0;
			gcpft_amon=0.0;

			// GUESSN
			gcpft_cmass_leaf=0.0;
			gcpft_nmass_leaf=0.0;
			gcpft_nmass=0.0;
			gcpft_nuptake=0.0;
			gcpft_anppn=0.0;
			gcpft_vmaxnlim=0.0;
			gcpft_nlim=0.0;

			// end GUESSN

			double heightindiv_total = 0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
				standpft_cmass=0.0;
				standpft_anpp=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total=0.0;
				standpft_aiso=0.0;
				standpft_amon=0.0;
				//GUESSN
				standpft_cmass_leaf=0.0;
				standpft_nmass_leaf=0.0;
				standpft_nmass=0.0;
				standpft_nuptake=0.0;
				standpft_anppn=0.0;
				standpft_vmaxnlim=0.0;
				standpft_nlim=0.0;
				standpft_anpp_no_nlim=0.0;
				//GUESSN END

				// Initialise age structure array

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++){
						standpft_densindiv_ageclass[c]=0.0;
					}
				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch=stand.getobj();
					Vegetation& vegetation=patch.vegetation;

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv=vegetation.getobj();

						// guess2008 - alive check added
						if (indiv.id!=-1 && indiv.alive) {

							if (indiv.pft.id==pft.id) {
								standpft_cmass+=indiv.cmass_leaf+
									indiv.cmass_root+indiv.cmass_sap+indiv.cmass_heart-indiv.cmass_debt;

								standpft_anpp+=indiv.anpp;
								standpft_lai+=indiv.lai;
								standpft_aiso+=indiv.aiso;
								standpft_amon+=indiv.amon;

								// GUESSN N budget
								nmass_indiv+=(indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap+
									indiv.nmass_heart+indiv.nmass_reserve+indiv.nstore)/(double)npatch;
								// GUESS N
								standpft_cmass_leaf+=indiv.cmass_leaf*indiv.densindiv;
								standpft_nmass_leaf+=indiv.nmass_leaf*indiv.densindiv;
								standpft_vmaxnlim+=indiv.avmaxnlim*indiv.cmass_leaf*indiv.densindiv;
								standpft_anpp_no_nlim+=indiv.anpp/indiv.limnfact;
								standpft_nuptake+=indiv.nuptake;
								standpft_anppn+=indiv.ndemand;
								standpft_nmass+=indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap+
									indiv.nmass_heart+indiv.nmass_reserve;
								// end GUESSN

								// WOLF
								if (date.year == 605 && pft.lifeform==TREE && indiv.age>0) {

									int angio;
									if (pft.id < 4)
										angio = 1;
									else
										angio = 2;

									output_channel->add_value(out_allometry_ind, angio);
									output_channel->add_value(out_allometry_ind, indiv.age);														  // Age
									output_channel->add_value(out_allometry_ind, indiv.densindiv*10000.0/(double)stand.npatch());                         // N
									output_channel->add_value(out_allometry_ind, indiv.cmass_leaf/indiv.densindiv);                                   // Mfol
									output_channel->add_value(out_allometry_ind, (indiv.cmass_sap+indiv.cmass_heart)*0.75/indiv.densindiv);           // Mstem
									output_channel->add_value(out_allometry_ind, (indiv.cmass_sap+indiv.cmass_heart)*0.25/indiv.densindiv);	          // Mcroot
									output_channel->add_value(out_allometry_ind, indiv.cmass_root/indiv.densindiv);                                   // Mfroot
									output_channel->add_value(out_allometry_ind, 
									                          ((indiv.cmass_sap+indiv.cmass_heart)*0.25+indiv.cmass_root)/indiv.densindiv);           // Mroot
									output_channel->add_value(out_allometry_ind, (indiv.cmass_sap+indiv.cmass_heart)/indiv.densindiv);                // Mwood
									output_channel->add_value(out_allometry_ind, (indiv.cmass_leaf+indiv.cmass_root)/indiv.densindiv);		          // Mactive
									output_channel->add_value(out_allometry_ind, 
									                          (indiv.cmass_sap+indiv.cmass_heart+indiv.cmass_leaf+indiv.cmass_root)/indiv.densindiv); // M
									output_channel->add_value(out_allometry_ind, indiv.height);                                                       // Height
									output_channel->add_value(out_allometry_ind, pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3));        // Diameter

									output_channel->finish_row(out_allometry_ind, lon, lat, date.year);
								}

								if (vegmode==COHORT || vegmode==INDIVIDUAL) {

									// Age structure

									c=(int)(indiv.age/estinterval); // guess2008
									if (c<OUTPUT_MAXAGECLASS)
										standpft_densindiv_ageclass[c]+=indiv.densindiv;

									// guess2008 - only count trees with a trunk above a certain diameter
									if (pft.lifeform==TREE && indiv.age>0) {
										double diam=pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3);
										if (diam>0.03) {
											standpft_densindiv_total+=indiv.densindiv; // indiv/m2

											heightindiv_total+=indiv.height * indiv.densindiv;
										}
									}
								}
							}

						} // alive?
						vegetation.nextobj();
					}

					stand.nextobj();
				} // end of patch loop

				standpft_cmass/=(double)stand.npatch();
				standpft_anpp/=(double)stand.npatch();
				standpft_lai/=(double)stand.npatch();
				standpft_densindiv_total/=(double)stand.npatch();
				standpft_aiso/=(double)stand.npatch();
				standpft_amon/=(double)stand.npatch();

				// GUESSN
				standpft_nuptake/=(double)stand.npatch();
				standpft_anppn/=(double)stand.npatch();
				standpft_nmass/=(double)stand.npatch();
				standpft_anpp_no_nlim/=(double)stand.npatch();

				if (standpft_anpp_no_nlim>0.0 && standpft_anpp>0.0)
					standpft_nlim=standpft_anpp/standpft_anpp_no_nlim;
				else
					standpft_nlim=1.0;
				
				gcpft_cmass_leaf+=standpft_cmass_leaf;
				gcpft_nmass_leaf+=standpft_nmass_leaf;
				gcpft_vmaxnlim+=standpft_vmaxnlim;
				gcpft_nlim+=standpft_nlim;

				cmass_leaf_gridcell+=standpft_cmass_leaf;
				nmass_leaf_gridcell+=standpft_nmass_leaf;
				vmaxnlim_gridcell+=standpft_vmaxnlim;
				nlim_gridcell+=standpft_nlim;

				if (!negligible(standpft_cmass_leaf))
					standpft_vmaxnlim/=standpft_cmass_leaf;
				
				// end GUESSN

				heightindiv_total/=(double)stand.npatch();

				//Update landcover totals
				landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
				landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
				landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();
				landcover_aiso[stand.landcover]+=standpft_aiso*stand.get_landcover_fraction();
				landcover_amon[stand.landcover]+=standpft_amon*stand.get_landcover_fraction();

				landcover_cmass_leaf[stand.landcover]+=standpft_cmass_leaf*stand.get_landcover_fraction();
				landcover_nmass_leaf[stand.landcover]+=standpft_nmass_leaf*stand.get_landcover_fraction();
				landcover_nmass[stand.landcover]+=standpft_nmass*stand.get_landcover_fraction();
				landcover_nuptake[stand.landcover]+=standpft_nuptake*stand.get_landcover_fraction();
				landcover_anppn[stand.landcover]+=standpft_anppn*stand.get_landcover_fraction();
				landcover_vmaxnlim[stand.landcover]+=standpft_vmaxnlim*stand.get_landcover_fraction();
				landcover_nlim[stand.landcover]+=standpft_nlim*stand.get_landcover_fraction();

				//Update pft totals
				gcpft_cmass+=standpft_cmass;
				gcpft_anpp+=standpft_anpp;
				gcpft_lai+=standpft_lai;
				gcpft_densindiv_total+=standpft_densindiv_total;
				gcpft_aiso+=standpft_aiso;
				gcpft_amon+=standpft_amon;

				gcpft_nmass+=standpft_nmass;
				gcpft_nuptake+=standpft_nuptake;
				gcpft_anppn+=standpft_anppn;

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++)
						gcpft_densindiv_ageclass[c]+=standpft_densindiv_ageclass[c];

				// Update gridcell totals
				double fraction_of_gridcell = stand.get_gridcell_fraction();
				cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
				anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
				lai_gridcell+=standpft_lai*fraction_of_gridcell;
				dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
				aiso_gridcell+=standpft_aiso*fraction_of_gridcell;
				amon_gridcell+=standpft_amon*fraction_of_gridcell;

				nmass_gridcell+=standpft_nmass*fraction_of_gridcell;
				nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
				anppn_gridcell+=standpft_anppn*fraction_of_gridcell;
				anpp_no_nlim_gridcell+=standpft_anpp_no_nlim*fraction_of_gridcell;

				if (!out_canopyh.invalid()) {
					canopyheight_gridcell=canopy_height(stand);
				}

				// Graphical output every 10 years
				// (Windows shell only - "plot" statements have no effect otherwise)
				if (true) {
					plot("cmass",pft.name,date.year,gcpft_cmass);
					plot("anpp",pft.name,date.year,gcpft_anpp);
					plot("lai",pft.name,date.year,gcpft_lai);
					if (gcpft_cmass_leaf>0.0 && ifnlim) {
						plot("vmax N lim",pft.name,date.year,gcpft_vmaxnlim/gcpft_cmass_leaf);
						plot("N lim on growth",pft.name,date.year,gcpft_nlim);
						plot("leaf C:N ratio",pft.name,date.year,gcpft_cmass_leaf/gcpft_nmass_leaf);
					}
				}
				gridcell.nextobj();
			}//End of loop through stands

			// Print PFT sums to files

			// GUESSN
			double gcpft_cton_leaf=0.0;
			if (gcpft_cmass_leaf>0.0) {
				gcpft_cton_leaf=gcpft_cmass_leaf/gcpft_nmass_leaf;
				gcpft_vmaxnlim/=gcpft_cmass_leaf;
			}
			
			out.add_value(out_cton,     gcpft_cton_leaf);
			out.add_value(out_vmaxnlim, gcpft_vmaxnlim);
			out.add_value(out_nlim,     gcpft_nlim);
			out.add_value(out_nmass,    gcpft_nmass);
			out.add_value(out_nuptake,  gcpft_nuptake);
			out.add_value(out_anppn,    gcpft_anppn);
			// end GUESSN

			out.add_value(out_cmass, gcpft_cmass);
			out.add_value(out_anpp,  gcpft_anpp);
			out.add_value(out_dens,  gcpft_densindiv_total);
			out.add_value(out_lai,   gcpft_lai);

			// print species heights
			double height = 0.0;
			if (gcpft_densindiv_total > 0.0)
				height = heightindiv_total/gcpft_densindiv_total;
			
			out.add_value(out_speciesheights, height);

			out.add_value(out_aiso, gcpft_aiso);
			out.add_value(out_amon, gcpft_amon);

			pftlist.nextobj();

		} // *** End of PFT loop ***

		flux_veg=flux_soil=flux_fire=flux_est=flux_harvest=0.0;

		// guess2008 - carbon pools
		c_litter=c_fast=c_slow=c_harv_slow=0.0;

		// GUESSN
		surfsoillitterc=surfsoillittern=cwdc=cwdn=microc=micron=humusc=humusn=centuryc=centuryn=n_litter=n_harv_slow=0.0;
		andep_gridcell=anmin_gridcell=animm_gridcell=anfix_gridcell=nsupply_gridcell=ndemand_gridcell=0.0;
		n_org_leach_gridcell=n_min_leach_gridcell=0.0;
		flux_nh3=flux_no=flux_no2=flux_n2o=flux_ntot=0.0;

		// end GUESSN

		// Sum C fluxes, dead C pools and runoff across patches

		gridcell.firstobj();

		// Loop through Stands
		while (gridcell.isobj) {
			Stand& stand=gridcell.getobj();
			stand.firstobj();

			//Loop through Patches
			while (stand.isobj) {
				Patch& patch=stand.getobj();

				double to_gridcell_average = stand.get_gridcell_fraction()/(double)stand.npatch();

				flux_veg+=patch.fluxes.acflux_veg*to_gridcell_average;
				flux_soil+=patch.fluxes.acflux_soil*to_gridcell_average;
				flux_fire+=patch.fluxes.acflux_fire*to_gridcell_average;
				flux_est+=patch.fluxes.acflux_est*to_gridcell_average;
				flux_harvest+=patch.fluxes.acflux_harvest*to_gridcell_average;

				// GUESSN
				flux_nh3+=patch.fluxes.aNH3_fire*to_gridcell_average;
				flux_no+=patch.fluxes.aNO_fire*to_gridcell_average;
				flux_no2+=patch.fluxes.aNO2_fire*to_gridcell_average;
				flux_n2o+=patch.fluxes.aN2O_fire*to_gridcell_average;	
				flux_ntot+=(patch.fluxes.aNH3_fire+patch.fluxes.aNO_fire+patch.fluxes.aNO2_fire+
					patch.fluxes.aN2O_fire)*to_gridcell_average;
				
				c_fast+=patch.soil.cpool_fast*to_gridcell_average;
				c_slow+=patch.soil.cpool_slow*to_gridcell_average;

				// Sum all litter
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft=patch.pft[q];
					c_litter+=(patchpft.litter_leaf+patchpft.litter_root+patchpft.litter_wood+patchpft.litter_repr)*to_gridcell_average;
					n_litter+=(patchpft.nmass_litter_leaf+patchpft.nmass_litter_root+patchpft.nmass_litter_wood)*to_gridcell_average;
				}

				//Sum slow pools of harvested products
				if(run_landcover && ifslowharvestpool)
				{
					for (int q=0;q<npft;q++)
					{
						Patchpft& patchpft=patch.pft[q];
						c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
						n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
					}
				}

				surfrunoff_gridcell+=patch.asurfrunoff*to_gridcell_average;
				drainrunoff_gridcell+=patch.adrainrunoff*to_gridcell_average;
				baserunoff_gridcell+=patch.abaserunoff*to_gridcell_average;
				runoff_gridcell+=patch.arunoff*to_gridcell_average;

				// Fire return time
				if (!iffire || patch.fireprob < 0.001)
					firert_gridcell+=1000.0/(double)stand.npatch(); // Set a limit of 1000 years
				else
					firert_gridcell+=(1.0/patch.fireprob)/(double)stand.npatch();

				// GUESSN
				andep_gridcell+=patch.soil.andep/(double)stand.npatch()*10000.0;	// convert from m2 to ha
				anmin_gridcell+=patch.soil.anmin/(double)stand.npatch()*10000.0;	// convert from m2 to ha
				animm_gridcell+=patch.soil.animmob/(double)stand.npatch()*10000.0; // convert from m2 to ha
				anfix_gridcell+=patch.soil.anfix/(double)stand.npatch()*10000.0;		// convert from m2 to ha
				n_min_leach_gridcell+=patch.soil.aminleach/(double)stand.npatch()*10000.0;	// convert from m2 to ha
				n_org_leach_gridcell+=patch.soil.aorgleach/(double)stand.npatch()*10000.0;	// convert from m2 to ha
				nsupply_gridcell+=patch.nsupply/(double)stand.npatch()*10000.0;			// convert from m2 to ha
				ndemand_gridcell+=patch.ndemand/(double)stand.npatch()*10000.0;			// convert from m2 to ha

				for (int r=0;r<NSOMPOOL;r++) {
					if (patch.soil.sompool[r].nmass > 0.0 && r<NSOMPOOL-1) {
						if(r==SURFMETA||r==SURFSTRUCT||r==SOILMETA||r==SOILSTRUCT){
							surfsoillitterc+=patch.soil.sompool[r].cmass/(double)stand.npatch();
							surfsoillittern+=patch.soil.sompool[r].nmass/(double)stand.npatch();
						}
						else if (r==SURFCWD) {
							cwdc+=patch.soil.sompool[r].cmass/(double)stand.npatch();
							cwdn+=patch.soil.sompool[r].nmass/(double)stand.npatch();
						}
						else {
							if (r==SURFMICRO||r==SOILMICRO) {
								microc+=patch.soil.sompool[r].cmass/(double)stand.npatch();
								micron+=patch.soil.sompool[r].nmass/(double)stand.npatch();
							} 

							if (r==SURFHUMUS){
								humusc+=patch.soil.sompool[r].cmass/(double)stand.npatch();
								humusn+=patch.soil.sompool[r].nmass/(double)stand.npatch();
							}
							
							centuryc+=patch.soil.sompool[r].cmass/(double)stand.npatch();
							centuryn+=patch.soil.sompool[r].nmass/(double)stand.npatch();
						}
					}
				}
				// end GUESSN

				// Monthly output variables

				for (m=0;m<12;m++) {
					maet[m] += patch.maet[m]*to_gridcell_average;
					mpet[m] += patch.mpet[m]*to_gridcell_average;
					mevap[m] += patch.mevap[m]*to_gridcell_average;
					mintercep[m] += patch.mintercep[m]*to_gridcell_average;
					mrunoff[m] += patch.mrunoff[m]*to_gridcell_average;
					mrh[m] += patch.fluxes.mcflux_soil[m]*to_gridcell_average;
					mwcont_upper[m] += patch.soil.mwcont[m][0]*to_gridcell_average;
					mwcont_lower[m] += patch.soil.mwcont[m][1]*to_gridcell_average;

					// guess2008 - average across stands to get mgpp and mra here.
					mgpp[m] += patch.fluxes.mcflux_gpp[m]*to_gridcell_average;
					mra[m] += patch.fluxes.mcflux_ra[m]*to_gridcell_average;
					// bvoc
					miso[m]+=patch.fluxes.miso[m]*to_gridcell_average;
					mmon[m]+=patch.fluxes.mmon[m]*to_gridcell_average;
				}


				// Calculate monthly NPP and LAI

				Vegetation& vegetation=patch.vegetation;

				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv=vegetation.getobj();

					// guess2008 - alive check added
					if (indiv.id!=-1 && indiv.alive) {

						for (m=0;m<12;m++) {
							mlai[m] += indiv.mlai[m]*to_gridcell_average;
						}

					} // alive?

					vegetation.nextobj();

				} // while/vegetation loop
				stand.nextobj();
			} // patch loop
			gridcell.nextobj();
		} // stand loop

		// In contrast to annual NEE, monthly NEE does not include fire
		// or establishment fluxes
		for (m=0;m<12;m++) {
			mnpp[m] = mgpp[m]-mra[m];
			mnee[m] = mnpp[m]-mrh[m];
		}

		// Print gridcell totals to files

		double cton_leaf_gridcell=0.0;
		if (cmass_leaf_gridcell>0.0) {
			cton_leaf_gridcell=cmass_leaf_gridcell/nmass_leaf_gridcell;
			vmaxnlim_gridcell/=cmass_leaf_gridcell;
		}
		flux_nconc=(flux_fire>0.0) ? flux_ntot/flux_fire : 0.0;

		out.add_value(out_cmass,  cmass_gridcell);
		out.add_value(out_anpp,   anpp_gridcell);
		out.add_value(out_dens,   dens_gridcell);
		out.add_value(out_lai,    lai_gridcell);
		out.add_value(out_firert, firert_gridcell);
		out.add_value(out_runoff, surfrunoff_gridcell);
		out.add_value(out_runoff, drainrunoff_gridcell);
		out.add_value(out_runoff, baserunoff_gridcell);
		out.add_value(out_runoff, runoff_gridcell);
		out.add_value(out_aiso,   aiso_gridcell);
		out.add_value(out_amon,   amon_gridcell);

		//GUESSN
		double nlim_grid=0.0;
		if (!negligible(anpp_no_nlim_gridcell))
			nlim_grid=anpp_gridcell/anpp_no_nlim_gridcell;

		out.add_value(out_nmass,     nmass_gridcell);
		out.add_value(out_anppn,     anppn_gridcell);
		out.add_value(out_cton,      cton_leaf_gridcell);
		out.add_value(out_vmaxnlim,  vmaxnlim_gridcell);
		out.add_value(out_nlim,      nlim_grid);
		out.add_value(out_canopyh,   canopyheight_gridcell);
		out.add_value(out_nuptake,   nuptake_gridcell);

		out.add_value(out_nsources,	andep_gridcell);
		out.add_value(out_nsources, anfix_gridcell);
		out.add_value(out_nsources, andep_gridcell+anfix_gridcell);
		out.add_value(out_nsources, anmin_gridcell);
		out.add_value(out_nsources, animm_gridcell);
		out.add_value(out_nsources, anmin_gridcell-animm_gridcell);
		out.add_value(out_nsources, andep_gridcell+anmin_gridcell-animm_gridcell+anfix_gridcell);
		out.add_value(out_nsources, ndemand_gridcell);

		out.add_value(out_nleach, n_min_leach_gridcell);
		out.add_value(out_nleach, n_org_leach_gridcell);
		out.add_value(out_nleach, n_min_leach_gridcell+n_org_leach_gridcell);
		out.add_value(out_nleach, gridcell.climate.atemp_mean);
		out.add_value(out_nleach, gridcell.climate.aprec);
		out.add_value(out_nleach, gridcell.climate.asun);

		// end GUESSN

		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {
					out.add_value(out_cmass, landcover_cmass[i]);
					out.add_value(out_anpp,  landcover_anpp[i]);
					out.add_value(out_dens,  landcover_densindiv_total[i]);
					out.add_value(out_lai,   landcover_lai[i]);
					out.add_value(out_aiso,  landcover_aiso[i]);
					out.add_value(out_amon,  landcover_amon[i]);

					//GUESSN
					out.add_value(out_nmass, landcover_nmass[i]);
					out.add_value(out_anppn, landcover_anppn[i]);

					double landcover_cton_leaf=0.0;
					if (landcover_cmass_leaf[i]>0.0) {
						landcover_cton_leaf=landcover_cmass_leaf[i]/landcover_nmass_leaf[i];
						landcover_vmaxnlim[i]/=landcover_cmass_leaf[i];
					}

					out.add_value(out_cton, landcover_cton_leaf);
					out.add_value(out_vmaxnlim, landcover_vmaxnlim[i]);
					out.add_value(out_nlim, landcover_nlim[i]);
					out.add_value(out_nuptake, landcover_nuptake[i]);
					// end GUESSN
				}
			}
		}

		// Print monthly output variables
		for (m=0;m<12;m++) {
			 out.add_value(out_mnpp,         mnpp[m]);
			 out.add_value(out_mlai,         mlai[m]);
			 out.add_value(out_mgpp,         mgpp[m]);
			 out.add_value(out_mra,          mra[m]);
			 out.add_value(out_maet,         maet[m]);
			 out.add_value(out_mpet,         mpet[m]);
			 out.add_value(out_mevap,        mevap[m]);
			 out.add_value(out_mrunoff,      mrunoff[m]);
			 out.add_value(out_mintercep,    mintercep[m]);
			 out.add_value(out_mrh,          mrh[m]);
			 out.add_value(out_mnee,         mnee[m]);
			 out.add_value(out_mwcont_upper, mwcont_upper[m]);
			 out.add_value(out_mwcont_lower, mwcont_lower[m]);
			 out.add_value(out_miso,         miso[m]);
			 out.add_value(out_mmon,         mmon[m]);
		}

		// Graphical output every 10 years
		// (Windows shell only - no effect otherwise)

		if (true) {
			gridcell.firstobj();
			if(gridcell.isobj)	//Fixed bug here if no stands were present.
			{
				Stand& stand=gridcell.getobj();
				plot("fluxes","flux_veg",date.year,flux_veg);
				plot("fluxes","flux_soil",date.year,flux_soil);
				plot("fluxes","flux_fire",date.year,flux_fire);
				plot("fluxes","flux_est",date.year,flux_est);
				plot("fluxes","NEE",date.year,flux_veg+flux_soil+flux_fire+flux_est);
				// GUESSN
				if (!ifcentury) {
					plot("soilc","slow",date.year,stand[0].soil.cpool_slow);
					plot("soilc","fast",date.year,stand[0].soil.cpool_fast);
				}
				else {
					plot("N addition (kgN/ha/yr)","Soil N fix",date.year,anfix_gridcell);
					plot("N addition (kgN/ha/yr)","N deposition",date.year,andep_gridcell);
					plot("N min-immob (kgN/ha/yr)","N",date.year,anmin_gridcell-animm_gridcell);

					plot("N demand/supply (kgN/ha/yr)","N supply",date.year,nsupply_gridcell);
					plot("N demand/supply (kgN/ha/yr)","N demand",date.year,ndemand_gridcell);
				}
				//end GUESSN
			}
		}

		// Write fluxes to file

		out.add_value(out_cflux, flux_veg);
		out.add_value(out_cflux, flux_soil);
		out.add_value(out_cflux, flux_fire);
		out.add_value(out_cflux, flux_est);
		if (run_landcover) {
			 out.add_value(out_cflux, flux_harvest);
		}
		out.add_value(out_cflux, flux_veg+flux_soil+flux_fire+flux_est+flux_harvest);

		out.add_value(out_cpool, cmass_gridcell);
		out.add_value(out_cpool, c_litter);
		if (!ifcentury) {
			out.add_value(out_cpool, c_fast);
			out.add_value(out_cpool, c_slow);
		}
		else {
			out.add_value(out_cpool, surfsoillitterc);
			out.add_value(out_cpool, cwdc);
			out.add_value(out_cpool, microc);
			out.add_value(out_cpool, humusc);
			out.add_value(out_cpool, centuryc);
		}
		
		if (run_landcover && ifslowharvestpool) {
			out.add_value(out_cpool, c_harv_slow);
		}

		// Calculate total cpool, starting with cmass and litter...
		double cpool_total = cmass_gridcell + c_litter;

		// Add SOM pools
		if (!ifcentury) {
			cpool_total += c_fast + c_slow;
		}
		else {
			cpool_total += centuryc+surfsoillitterc+cwdc;
		}

		// Add slow harvest pool if needed
		if (run_landcover && ifslowharvestpool) {
			cpool_total += c_harv_slow;
		}

		out.add_value(out_cpool, cpool_total);

		// GUESSN
		if (ifcentury) {
			out.add_value(out_npool, nmass_gridcell);
			out.add_value(out_npool, n_litter);
			out.add_value(out_npool, surfsoillittern);
			out.add_value(out_npool, cwdn);
			out.add_value(out_npool, micron);
			out.add_value(out_npool, humusn);
			out.add_value(out_npool, centuryn);

			if(run_landcover && ifslowharvestpool) {
				out.add_value(out_npool, n_harv_slow);
				out.add_value(out_npool, nmass_gridcell+n_litter+surfsoillittern+cwdn+centuryn+n_harv_slow);
			}
			else {
				out.add_value(out_npool, nmass_gridcell+n_litter+centuryn);
			}
		}

		out.add_value(out_nflux, flux_nh3);
		out.add_value(out_nflux, flux_no);
		out.add_value(out_nflux, flux_no2);
		out.add_value(out_nflux, flux_n2o);
		out.add_value(out_nflux, flux_ntot);
		out.add_value(out_nflux, flux_nconc);
		// end GUESSN

		// Output of age structure (Windows shell only - no effect otherwise)

		if (vegmode==COHORT || vegmode==INDIVIDUAL) {

			if (!(date.year%20) && date.year<2000) {

				resetwindow("age_structure");
				gridcell.firstobj();

				// Loop through Stands
				while (gridcell.isobj) {
					Stand& stand=gridcell.getobj();

					if(stand.landcover==NATURAL){

						while (pftlist.isobj) {
							Pft& pft=pftlist.getobj();
							if(pft.landcover==NATURAL){
								if (pft.lifeform==TREE) {

									Gridcellpft& gridcellpft=gridcell.pft[pft.id];
									//TREES ONLY
									for (c=0;c<nclass;c++){
										plot("age_structure",pft.name,c*estinterval+estinterval/2,gcpft_densindiv_ageclass[c]/(double)stand.npatch());
									}
								}
							}
							pftlist.nextobj();

						}
					}
					gridcell.nextobj();
				}
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// TERMIO
// Called at end of model run (i.e. following simulation of all stands)

void termio() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.
	delete output_channel;

	// Clean up
	gridlist.killall();
}

#endif // USE_CRU_IO


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Galloway, J. N., F. J. Dentener, D. G. Capone, E. W. Boyer, R. W. Howarth, S. P. Seitzinger,
//   G. P. Asner, C. Cleveland, P. Green, E. Holland, D. M. Karl, A. F. Michaels, J. H. Porter, 
//   A. Townsend, and C. V?r?smarty. 2004.
//   Nitrogen Cycles: Past, Present and Future. Biogeochemistry 70: 153-226.
// Nakai, T., Sumida, A., Kodama, Y., Hara, T., Ohta, T. (2010). A comparison between
//   various definitions of forest stand height and aerodynamic canopy height.
//   Agricultural and Forest Meteorology, 150(9), 1225-1233
