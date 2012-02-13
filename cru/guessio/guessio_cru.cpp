///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                LPJ-GUESS input/output module with input from instruction
//                        script
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
//                        Includes Dieter G:s latest updates 021121
//                        Version compatible with LPJ-GUESS version 2.1
//                        (excludes PFT paramter twmax)
//                        Updated 20050125: last line in output files ends in newline
// Header file name:      guessio.h
// Source code file name: guessio.cpp
// Written by:            Ben Smith
// Version dated:         2003-07-22/2005-01-25
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

#ifdef USE_CRU

#include "guessio.h"

#include "driver.h"
#include <plib.h>
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>


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
#include "GlobalNitrogenDepositionRCP45.h"
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


///////////////////////////////////////////////////////////////////////////////////////
// CLASS PARAMLIST
// Functionality for storing and retrieving custom "param" items from the instruction
// script

struct Paramtype {
	xtring name;
	xtring str;
	double num;
};

class Paramlist : public ListArray<Paramtype> {

public:
	void addparam(xtring& name,xtring& value) {
		Paramtype& p=createobj();
		p.name=name.lower();
		p.str=value;
	}

	void addparam(xtring& name,double value) {
		Paramtype& p=createobj();
		p.name=name.lower();
		p.num=value;
	}

	Paramtype& operator[](xtring name) {
		name.lower();
		firstobj();
		while (isobj) {
			Paramtype& p=getobj();
			if (p.name==name) return p;
			nextobj();
		}
		fail("Paramlist::operator[]: parameter \"%s\" not found",(char*)name);

		// This point cannot be reached in practice, but to satisfy more pedantic
		// compilers ...

		return getobj();
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// ENUM DECLARATIONS OF INTEGER CONSTANTS FOR PLIB INTERFACE

enum {BLOCK_GLOBAL,BLOCK_PFT,BLOCK_PARAM};
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_LANDCOVER,CB_PHENOLOGY,CB_PATHWAY,	
	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM};


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
xtring file_firert;

xtring file_dgpp;

// GUESSN
xtring file_cton,file_nmass,file_nsources,file_npool,file_nleach,file_nuptake,file_anppn,file_vmaxnlim,file_nlim;
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
	file_cpool=file_firert="";
	file_dgpp="";

	// GUESSN
	file_cton=file_nmass=file_nsources=file_npool=file_nleach=file_nuptake=file_anppn=file_vmaxnlim=file_nlim="";
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

		// GUESSN
		declareitem("nrelocfrac",&nrelocfrac,0.0,1.0,1,CB_NONE,
			"Fractional N relocation from shed leaves & roots");
		declareitem("ifnfix",&ifnfix,0,3,1,CB_NONE,
			"Whether to include an estimate for N fixation");

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

		// SENS
		declareitem("sens_cton_needle",&sens_cton_needle,0.5,2.0,1,CB_NONE,
			"Needleleaved C:N min change");
		declareitem("sens_cton_broad",&sens_cton_broad,0.5,2.0,1,CB_NONE,
			"Broadleaved C:N min change");
		declareitem("sens_cton_vmax",&sens_cton_vmax,0.5,2.0,1,CB_NONE,
			"vmax N limitation effect on leaf C:N");
		declareitem("sens_decayrate",&sens_decayrate,0.5,2.0,1,CB_NONE,
			"Change decay rates constant of som pools");

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
		// end GUESSN

		// GUESSN allometry
		declareitem("file_allometry",&file_allometry,300,CB_NONE,"Allometry output file");
		declareitem("file_canopyh",&file_canopyh,300,CB_NONE,"Canopy height output file");
		declareitem("file_allometry_ind",&file_allometry_ind,300,CB_NONE,"Individual Allometry output file");
		// end GUESSN

		
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

		declareitem("file_dgpp",&file_dgpp,300,CB_NONE,"Daily GPP output file");

		// guess2008 - new options
		declareitem("ifsmoothgreffmort",&ifsmoothgreffmort,1,CB_NONE,
			"Whether to vary mort_greff smoothly with growth efficiency (0,1)");
		declareitem("ifdroughtlimitedestab",&ifdroughtlimitedestab,1,CB_NONE,
			"Whether establishment drought limited (0,1)");
		declareitem("ifrainonwetdaysonly",&ifrainonwetdaysonly,1,CB_NONE,
			"Whether it rains on wet days only (1), or a little every day (0);");
		declareitem("ifspeciesspecificwateruptake",&ifspeciesspecificwateruptake,1,CB_NONE,
			"Whether or not there is species specific soil water uptake (0,1)");
		declareitem("searchradius", &searchradius, 0, 100, 1, CB_NONE,
			"If specified, CRU data will be searched for in a circle");

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
		// end GUESSN

		declareitem("reprfrac",&ppft->reprfrac,0.0,1.0,1,CB_NONE,
			"Fraction of NPP allocated to reproduction");
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

		declareitem("harv_eff",&ppft->harv_eff,0.0,1.0,1,CB_NONE,"Harvest efficiency");
		declareitem("harvest_slow_frac",&ppft->harvest_slow_frac,0.0,1.0,1,CB_NONE,
			"Fraction of harvested products that goes into carbon depository for long-lived products like wood");
		declareitem("turnover_harv_prod",&ppft->turnover_harv_prod,0.0,1.0,1,CB_NONE,"Harvested products turnover (fraction/year)");
		declareitem("res_outtake",&ppft->res_outtake,0.0,1.0,1,CB_NONE,"ï¿½Fraction of residue outtake at harvest");

		callwhendone(CB_CHECKPFT);
		
		break;

	case BLOCK_PARAM:

		paramname=setname;
		declareitem("str",&strparam,80,CB_STRPARAM,
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

		// GUESSN
		if (!itemparsed("nrelocfrac")) badins("nrelocfrac");
		if (!itemparsed("ifnfix")) badins("ifnfix");

		if (!itemparsed("ifcentury")) badins("ifcentury");
		if (!itemparsed("ifnlim")) badins("ifnlim");
		if (!itemparsed("freenyears")) badins("freenyears");
		if (!itemparsed("ifleachn")) badins("ifleachn");
		if (!itemparsed("ifindiv_fnuptake")) badins("ifindiv_fnuptake");
		// end GUESSN
	
		// SENS
		if (!itemparsed("sens_cton_needle")) badins("sens_cton_needle");
		if (!itemparsed("sens_cton_broad")) badins("sens_cton_broad");
		if (!itemparsed("sens_decayrate")) badins("sens_decayrate");
		if (!itemparsed("sens_cton_vmax")) badins("sens_cton_vmax");

		// guess2008
		if (!itemparsed("outputdirectory")) badins("outputdirectory");
		if (!itemparsed("ifsmoothgreffmort")) badins("ifsmoothgreffmort");
		if (!itemparsed("ifdroughtlimitedestab")) badins("ifdroughtlimitedestab");
		if (!itemparsed("ifrainonwetdaysonly")) badins("ifrainonwetdaysonly");
		if (!itemparsed("ifspeciesspecificwateruptake")) badins("ifspeciesspecificwateruptake");

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
		// end GUESSN

		if (!itemparsed("reprfrac")) badins("reprfrac");
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
// bool getstand(Stand& stand)
//   Obtains latitude and soil static parameters for the next stand (grid cell) to
//   simulate. The function should returns false if no stands remain to be simulated,
//   otherwise true. Currently the following member variables of stand should be
//   initialised: members lat and instype of member climate; the following members of
//   member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
//   thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
//   code (Sitch et al 2000) by a call to function soilparameters in the driver
//   module (driver.cpp):
//
//   soilparameters(stand.soiltype,soilcode);
//
//   If the model is to be driven by quasi-daily values of the climate variables
//   derived from monthly means, this function may be the appropriate place to
//   perform the required interpolations. The utility function interp_climate in
//   driver.cpp may be called for this purpose:
//
//   interp_climate(mtemp,mprec,msun,dtemp,dprec,dsun);
//
//   This assumes the following arrays are declared, presumably at file scope:
//
//   double mtemp[12]   monthly average temperature (deg C)
//   double mprec[12]   monthly precipitation sum (mm)
//   double msun[12]    monthly average sunshine (%)
//   double dtemp[365]  daily interpolated temperature (deg C)
//   double dprec[365]  daily interpolated rainfall (mm)
//   double dsun[365]   daily interpolated sunshine (%)
//
// bool getclimate(Stand& stand)
//   Obtains climate data (including atmospheric CO2 and insolation) for this day.
//   The function should returns false if the simulation is complete for this stand,
//   otherwise true. This will normally require querying the year and day member
//   variables of the global class object date:
//
//   if (date.day==0 && date.year==nyear_spinup) return false; // guess2008
//   // else
//   return true;
//
//   Currently the following member variables of the climate member of stand must be
//   initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
//   values of the climate variables derived from monthly means, this day's values
//   will presumably be extracted from arrays containing the interpolated daily
//   values (see function getstand):
//
//   gridcell.climate.temp=dtemp[date.day];
//   gridcell.climate.prec=dprec[date.day];
//   gridcell.climate.insol=dsun[date.day];
//
// void outannual(Stand& stand,Pftlist& pftlist)
//   Called at the end of the last day of each simulation year to permit output of
//   model results.
//
// termio()
//   Called after simulation is complete for all stands to allow memory deallocation,
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
bool firstgrid; // whether simulating first grid cell in linked list

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


	void regress(double* x,double* y,int n,double& a,double& b) {

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
const int NYEAR_HIST=NYEAR_CMIP5; // guess2008 - CRU TS 3.0 has 106 years of data (1901-2006)
	// number of years of historical climate in CRU and CO2 files (see below)
const int FIRSTHISTYEAR=FIRSTHISTYEAR_CMIP5;
	// calender year corresponding to first year in CRU climate data set
const int NYEAR_SPINUP_DATA=30;
	// number of years to use for temperature-detrended spinup data set
	// (not to be confused with the number of years to spinup model for, which
	// is read from the ins file)

// Stream pointer to binary CRU historical climate data file (read from ins file)
FILE *in_cru;

// Full pathname of ASCII file containing annual CO2 values (read from ins file)
xtring file_co2;

// Output streams
FILE *out_cmass,*out_anpp,*out_lai,*out_cflux,*out_cpool,*out_runoff,*out_dens;
FILE *out_mnpp,*out_mlai,*out_mgpp,*out_mra,*out_maet,*out_mpet,*out_mevap,*out_mrunoff,*out_mintercep,*out_mrh;
FILE *out_mnee,*out_mwcont_upper,*out_mwcont_lower; 
FILE *out_firert; 

FILE *out_dgpp;

// GUESSN
// Full pathname of ASCII file containing annual N deposition values (read from ins file)
xtring file_ndep;

FILE *out_cton,*out_nmass, *out_nsources, *out_npool, *out_nleach, *out_nuptake, *out_anppn, *out_vmaxnlim, *out_nlim;
// end GUESSN

// GUESSN allometry
FILE *out_allometry, *out_canopyh, *out_allometry_ind;
// end GUESSN

// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

// CO2 data for each year of historical data set
double co2[NYEAR_HIST];

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

bool annual_output;
	// whether output should occur each simulation year (true) or at end of simulation
	// for each grid cell only (false)

// guess2008 - make file_cru and file_cru_misc global variables
xtring file_cru;
xtring file_cru_misc;

//AA CMIP5 - climate input
xtring file_cmip5hist;
xtring file_cmip5scen;

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
	

	double dummy1[1][12], dummy2[1][12], ddummy1[365], ddummy2[365], dsun[365];
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
	double daylength_save;


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


			for (int nday=day_start_month[m]; nday<day_start_month[m+1]; nday++) //(int nday=0; nday<=365; nday++)
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

	int i, y ,m;


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

	int i, y;
	
	for (int m = 0; m < 12; m++) 
	{
		ctemp[m] = 0.0; 
		cprec[m] = 0.0;  
		cswrad[m] = 0.0; 		
	}

	for (y=111;y<141;y++) //1961-1990
	{
		for (int m = 0; m < 12; m++) 
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
				ark.close(); // I.e. we opened it but we couldn´t rewind
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
				ark.close(); // I.e. we opened it but we couldn´t rewind
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
		double x[n], tempy[n], precy[n], swrady[n];
		double a_temp, b_temp, a_prec, b_prec, a_swrad, b_swrad;
		double anom_temp, anom_prec, anom_swrad;
		

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
		double trend_temp, trend_prec, trend_swrad;
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
		if (correctionmethod=="c6")
				
		{
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
				
				co2y[y]=co2[y+11];
			}
		
			// regress data to remove trend
			regress_data(x,co2y,n,a_co2,b_co2);

			// remove trend and fill 1871-2100 with detrended pre industrial data
			for (y=21;y<NYEAR_CMIP5;y++) {

				anom_co2=(double)(y%10)*b_co2;

				co2[y]=co2[11+y%10]-anom_co2;
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

void readco2_cmip5() {

	// Reads in atmospheric CO2 concentrations for historical period
	// from ascii text file with records in format: <year> <co2-value>

	int year,calender_year;

	// Retrieve name of CO2 file from ins file
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

	FILE* in=fopen(filename,"rt");
	if (!in) fail("readco2: 111 could not open CO2 file %s for input",
		(char*)filename);

	for (year=0;year<NYEAR_HIST;year++) {
		readfor(in,"i,f",&calender_year,&co2[year]);
		if (calender_year!=FIRSTHISTYEAR+year)
			fail("readco2: 222 %s, line %d - incorrect year specified",
				(char*)filename,year+1);
	}

	fclose(in);
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

	int target_ilon=dlon*10.0;
	int target_ilat=dlat*10.0;

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldnï¿½t rewind
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
				ark.close(); // I.e. we opened it but we couldnï¿½t rewind
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
	for (int i = 0; i < search_points.size(); i++) {
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




void readco2() {

	// Reads in atmospheric CO2 concentrations for historical period
	// from ascii text file with records in format: <year> <co2-value>

	int year,calender_year;

	// Retrieve name of CO2 file from ins file
	xtring filename=param["file_co2"].str;

	FILE* in=fopen(filename,"rt");
	if (!in) fail("readco2: 333 could not open CO2 file %s for input",
		(char*)filename);

	for (year=0;year<NYEAR_CRU;year++) {
		readfor(in,"i,f",&calender_year,&co2[year]);
		if (calender_year!=FIRSTHISTYEAR+year)
			fail("readco2: 444 %s, line %d - incorrect year specified",
				(char*)filename,year+1);
	}

	fclose(in);
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

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;

	//AA CMIP5
	correctionmethod=param["correctionmethod"].str;
	gcm=param["gcm"].str;
	rcp=param["rcp"].str;

	file_cmip5hist=param["path_cmip5hist"].str;
	file_cmip5scen=param["path_cmip5scen"].str;

	if (rcp=="45")
		file_cmip5scen+="mpi_esm_lr_rcp45_r1i1p1\\";
	else if (rcp=="85")
		file_cmip5scen+="mpi_esm_lr_rcp85_r1i1p1\\";
	else fail("CMIP5 scenario file not valid");

	file_cmip5hist+="cmip5_hist.bin";
	file_cmip5scen+="cmip5_scen.bin";

	path_cmip5_co2=param["path_cmip5_co2"].str;

	if (!ifcmip5) {
		// Read CO2 data from file
		readco2();
	}
	else {
		// AA CMIP5
		readco2_cmip5();
	}

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
		if (rcp=="45")
			file_ndep_cmip5=file_ndep+"RCP45.bin";
		else if (rcp=="85")
			file_ndep_cmip5=file_ndep+"RCP85.bin";
		else fail("N dep file not valid");

		FILE* in_ndep=fopen(file_ndep_cmip5,"rt");
		if (!in_ndep)
			fail("initio: could not open %s for input",(char*)file_ndep_cmip5);

		fclose(in_ndep);
	}

	// Remember whether to produce output each year or not
	annual_output=param["annual_output"].num;


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

	// guess2008
	// Retrieve output file names as read from ins file

	// We MUST have an output directory
	if (outputdirectory=="") {
		fail("No output directory given in the .ins file!");
	}


	// *** ANNUAL OUTPUT VARIABLES ***

	if (file_cmass!="") {
		file_cmass = outputdirectory + file_cmass;
		out_cmass=fopen(file_cmass,"w");
		if (!out_cmass) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_cmass);
	}
	else out_cmass=NULL;

	if (file_anpp!="") {
		file_anpp = outputdirectory + file_anpp;
		out_anpp=fopen(file_anpp,"w");
		if (!out_anpp) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_anpp);
	}
	else out_anpp=NULL;

	if (file_dens!="") {
		file_dens = outputdirectory + file_dens;
		out_dens=fopen(file_dens,"w");
		if (!out_dens) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_dens);
	}
	else out_dens=NULL;

	if (file_lai!="") {
		file_lai = outputdirectory + file_lai;
		out_lai=fopen(file_lai,"w");
		if (!out_lai) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_lai);
	}
	else out_lai=NULL;

	if (file_cflux!="") {
		file_cflux = outputdirectory + file_cflux;
		out_cflux=fopen(file_cflux,"w");
		if (!out_cflux) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_cflux);
	}
	else out_cflux=NULL;

	if (file_cpool!="") {
		file_cpool = outputdirectory + file_cpool;
		out_cpool=fopen(file_cpool,"w");
		if (!out_cpool) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_cpool);
	}
	else out_cpool=NULL;

	if (file_firert!="") {
		file_firert = outputdirectory + file_firert;
		out_firert=fopen(file_firert,"w");
		if (!out_firert) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_firert);
	}
	else out_firert=NULL;
	

	if (file_runoff!="") {
		file_runoff = outputdirectory + file_runoff;
		out_runoff=fopen(file_runoff,"w");
		if (!out_runoff) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_runoff);
	}
	else out_runoff=NULL;

	// GUESSN
	if (file_cton!="") {
		file_cton = outputdirectory + file_cton;
		out_cton=fopen(file_cton,"w");
		if (!out_cton) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_cton);
	}
	else out_cton=NULL;

	if (file_nmass!="") {
		file_nmass = outputdirectory + file_nmass;
		out_nmass=fopen(file_nmass,"w");
		if (!out_nmass) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_nmass);
	}
	else out_nmass=NULL;

	if (file_nsources!="") {
		file_nsources = outputdirectory + file_nsources;
		out_nsources=fopen(file_nsources,"w");
		if (!out_nsources) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_nsources);
	}
	else out_nsources=NULL;

	if (file_npool!="" && ifcentury) {
		file_npool = outputdirectory + file_npool;
		out_npool=fopen(file_npool,"w");
		if (!out_npool) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_npool);
	}
	else out_npool=NULL;
	
	if (file_nleach!="" && ifcentury) {
		file_nleach = outputdirectory + file_nleach;
		out_nleach=fopen(file_nleach,"w");
		if (!out_nleach) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_nleach);
	}
	else out_nleach=NULL;

	if (file_nuptake!="") {
		file_nuptake = outputdirectory + file_nuptake;
		out_nuptake=fopen(file_nuptake,"w");
		if (!out_nuptake) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_nuptake);
	}
	else out_nuptake=NULL;

	if (file_anppn!="") {
		file_anppn = outputdirectory + file_anppn;
		out_anppn=fopen(file_anppn,"w");
		if (!out_anppn) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_anppn);
	}
	else out_anppn=NULL;

	if (file_vmaxnlim!="") {
		file_vmaxnlim = outputdirectory + file_vmaxnlim;
		out_vmaxnlim=fopen(file_vmaxnlim,"w");
		if (!out_vmaxnlim) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_vmaxnlim);
	}
	else out_vmaxnlim=NULL;

	if (file_nlim!="") {
		file_nlim = outputdirectory + file_nlim;
		out_nlim=fopen(file_nlim,"w");
		if (!out_nlim) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_nlim);
	}
	else out_nlim=NULL;
	// end GUESSN

	// GUESSN allometry
	if (file_allometry!="") {
		file_allometry = outputdirectory + file_allometry;
		out_allometry=fopen(file_allometry,"w");
		if (!out_allometry) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_allometry);
	}
	else out_allometry=NULL;

	if (file_allometry_ind!="") {
		file_allometry_ind = outputdirectory + file_allometry_ind;
		out_allometry_ind=fopen(file_allometry_ind,"w");
		if (!out_allometry_ind) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_allometry_ind);
	}
	else out_allometry_ind=NULL;
	

	if (file_canopyh!="") {
		file_canopyh = outputdirectory + file_canopyh;
		out_canopyh=fopen(file_canopyh,"w");
		if (!out_canopyh) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_canopyh);
	}
	else out_canopyh=NULL;
	// end GUESSN


	// *** MONTHLY OUTPUT VARIABLES ***

	if (file_mnpp!="") {
		file_mnpp = outputdirectory + file_mnpp;
		out_mnpp=fopen(file_mnpp,"w");
		if (!out_mnpp) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mnpp);
	}
	else out_mnpp=NULL;

	if (file_mlai!="") {
		file_mlai = outputdirectory + file_mlai;
		out_mlai=fopen(file_mlai,"w");
		if (!out_mlai) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mlai);
	}
	else out_mlai=NULL;

	if (file_mgpp!="") {
		file_mgpp = outputdirectory + file_mgpp;
		out_mgpp=fopen(file_mgpp,"w");
		if (!out_mgpp) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mgpp);
	}
	else out_mgpp=NULL;

	if (file_mra!="") {
		file_mra = outputdirectory + file_mra;
		out_mra=fopen(file_mra,"w");
		if (!out_mra) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mra);
	}
	else out_mra=NULL;

	if (file_maet!="") {
		file_maet = outputdirectory + file_maet;
		out_maet=fopen(file_maet,"w");
		if (!out_maet) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_maet);
	}
	else out_maet=NULL;

	if (file_mpet!="") {
		file_mpet = outputdirectory + file_mpet;
		out_mpet=fopen(file_mpet,"w");
		if (!out_mpet) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mpet);
	}
	else out_mpet=NULL;

	if (file_mevap!="") {
		file_mevap = outputdirectory + file_mevap;
		out_mevap=fopen(file_mevap,"w");
		if (!out_mevap) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mevap);
	}
	else out_mevap=NULL;

	if (file_mrunoff!="") {
		file_mrunoff = outputdirectory + file_mrunoff;
		out_mrunoff=fopen(file_mrunoff,"w");
		if (!out_mrunoff) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mrunoff);
	}
	else out_mrunoff=NULL;

	if (file_mintercep!="") {
		file_mintercep = outputdirectory + file_mintercep;
		out_mintercep=fopen(file_mintercep,"w");
		if (!out_mintercep) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mintercep);
	}
	else out_mintercep=NULL;

	if (file_mrh!="") {
		file_mrh = outputdirectory + file_mrh;
		out_mrh=fopen(file_mrh,"w");
		if (!out_mrh) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mrh);
	}
	else out_mrh=NULL;

	if (file_mnee!="") {
		file_mnee = outputdirectory + file_mnee;
		out_mnee=fopen(file_mnee,"w");
		if (!out_mnee) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mnee);
	}
	else out_mnee=NULL;

	if (file_mwcont_upper!="") {
		file_mwcont_upper = outputdirectory + file_mwcont_upper;
		out_mwcont_upper=fopen(file_mwcont_upper,"w");
		if (!out_mwcont_upper) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mwcont_upper);
	}
	else out_mwcont_upper=NULL;

	if (file_mwcont_lower!="") {
		file_mwcont_lower = outputdirectory + file_mwcont_lower;
		out_mwcont_lower=fopen(file_mwcont_lower,"w");
		if (!out_mwcont_lower) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_mwcont_lower);
	}
	else out_mwcont_lower=NULL;

	// Daily output
	if (file_dgpp!="") {
		file_dgpp = outputdirectory + file_dgpp;
		out_dgpp=fopen(file_dgpp,"w");
		if (!out_dgpp) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_dgpp);
	}
	else out_dgpp=NULL;

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);

	// Start at first object in linked list of grid cell coordinates ...
	firstgrid=true;
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

			if (rcp=="45") {

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
	// Obtains latitude and soil static parameters for the next grid cell to
	// simulate. The function should return false if no grid cells remain to be simulated,
	// otherwise true. Currently the following member variables of Gridcell should be
	// initialised: members lat and instype of member climate; the following members of
	// member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	// thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	// code (Sitch et al 2000) by a call to function soilparameters in the driver
	// module (driver.cpp):
	//
	// soilparameters(stand.soiltype,soilcode);
	//
	// If the model is to be driven by quasi-daily values of the climate variables
	// derived from monthly means, this function may be the appropriate place to
	// perform the required interpolations. The utility function interp_climate in
	// driver.cpp may be called for this purpose:
	//
	// interp_climate(mtemp,mprec,msun,dtemp,dprec,dsun);
	//
	// This assumes the following arrays are declared, presumably at file scope:
	//
	// double mtemp[12]   monthly average temperature (deg C)
	// double mprec[12]   monthly precipitation sum (mm)
	// double msun[12]    monthly average sunshine (%)
	// double dtemp[365]  daily interpolated temperature (deg C)
	// double dprec[365]  daily interpolated rainfall (mm)
	// double dsun[365]   daily interpolated sunshine (%)

	// Select coordinates for next grid cell in linked list
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	bool gridfound;
	bool LUerror=false;

	// guess2008 - run with the same randon number sequence each time
	setseed(12345678);

	if (firstgrid) {
		gridlist.firstobj();
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

			fail("Grid cell not found in %s",(char*)file_ndep);
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
		
		// Tell framework the latitude of this grid cell
		gridcell.climate.lat=gridlist.getobj().lat;
		
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		if (!ifcmip5)
			gridcell.climate.instype=SUNSHINE;
		else {
			gridcell.climate.instype=SWRAD;
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
void getlandcover(Gridcell& gridcell,Pftlist& pftlist)
{
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
				for(int i=0;i<NLANDCOVERTYPES;i++)
				{
					if(run[i])
						nactive_landcovertypes++;
				}
			}

			for(int i=0;i<NLANDCOVERTYPES;i++)
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
					for(int i=0;i<NLANDCOVERTYPES;i++)
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

					for(int i=0;i<NLANDCOVERTYPES;i++)
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
//				for(int i=0;i<NLANDCOVERTYPES;i++)
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
	// The function should returns false if the simulation is complete for this stand,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of stand must be
	// initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	// values of the climate variables derived from monthly means, this day's values
	// will presumably be extracted from arrays containing the interpolated daily
	// values (see function getstand):
	//
	// gridcell.climate.temp=dtemp[date.day];
	// gridcell.climate.prec=dprec[date.day];
	// gridcell.climate.insol=dsun[date.day];

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
			interp_climate(mtemp,mprec,msun,dtemp,dprec,dsun);

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
				dtemp,dprec,dsun);

			// guess2008 - only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year-nyear_spinup],dprec,hist_mwet[date.year-nyear_spinup]);
			}
		}

		if (!ifcmip5) {
			if (date.year<nyear_spinup)
				climate.co2=co2[0];
			else if (date.year<nyear_spinup+NYEAR_HIST)
				climate.co2=co2[date.year-nyear_spinup];
		}
		else {
			// CMIP5 - land use input
			if (date.year<nyear_spinup){ 
				climate.co2=co2[0];
				if(iflandusesimple)
				  climate.frluse=hist_frluse[0];
			}
			else if (date.year<nyear_spinup+NYEAR_HIST){
				climate.co2=co2[date.year-nyear_spinup];
				if(iflandusesimple)
				  climate.frluse=hist_frluse[date.year-nyear_spinup];
			}
		}

		climate.andep=0.0;
		int m;
		if (date.year<nyear_spinup){
			for (m=0;m<12;m++) {
				climate.andep+=(NHxDryDep[0][m]+NOyDryDep[0][m]+
					NHxWetDep[0][m]+NOyWetDep[0][m])*date.ndaymonth[m];
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

	climate.temp=dtemp[date.day];
	climate.prec=dprec[date.day];
	climate.insol=dsun[date.day];

	// First day of year only ...

	if (date.day==0) {

		// Return false if last year was the last for the simulation
		if (date.year==nyear_spinup+NYEAR_HIST) return false;

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

	int p,c,m,nclass;
	double flux_veg,flux_soil,flux_fire,flux_est,flux_harvest;
	double c_litter,c_fast,c_slow,c_harv_slow;

	// GUESSN
	double surfsoillitterc,surfsoillittern,cwdc,cwdn,microc,micron,humusc,humusn,centuryc,centuryn,n_harv_slow;
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
	
	// DGPP
	double dgpp[365][20];

	double lon,lat;

	if (vegmode==COHORT)
		nclass=min(date.year/estinterval+1,OUTPUT_MAXAGECLASS);

	if (date.year==0 && firstgrid) {

		// Very first time only

		// Print column labels
		// guess2008 - added runoff & dens

		const char* lonlatyearstr = "%8s%8s%8s"; // easier to change now.
		const char* lonlatyearstr_extended = "%8s%8s%8s%8s%8s%8s%8s%10s\n";

		if (out_cmass) fprintf(out_cmass,lonlatyearstr,"Lon","Lat","Year");
		if (out_anpp) fprintf(out_anpp,lonlatyearstr,"Lon","Lat","Year");
		if (out_lai) fprintf(out_lai,lonlatyearstr,"Lon","Lat","Year");
		if (out_runoff) fprintf(out_runoff,lonlatyearstr,"Lon","Lat","Year");
		if (out_dens) fprintf(out_dens,lonlatyearstr,"Lon","Lat","Year");
		if (out_cflux) {
			if(run_landcover)
				fprintf(out_cflux,"%8s%8s%8s%8s%8s%8s%8s%9s%10s\n","Lon","Lat","Year","Veg","Soil",
					"Fire","Est","Harvest","NEE");
			else
				fprintf(out_cflux,lonlatyearstr_extended,"Lon","Lat","Year","Veg","Soil",
					"Fire","Est","NEE");
		}
		if (out_cpool) {
			if (!ifcentury){
				if(run_landcover && ifslowharvestpool)
					fprintf(out_cpool,"%8s%8s%8s%8s%8s%8s%8s%10s%10s\n","Lon","Lat","Year","VegC","LittC",
						"SoilfC","SoilsC", "HarvSlowC","Total");
				else
					fprintf(out_cpool,lonlatyearstr_extended,"Lon","Lat","Year","VegC","LittC",
						"SoilfC","SoilsC","Total");
			}else{
				if(run_landcover && ifslowharvestpool)
					fprintf(out_cpool,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%10s\n","Lon","Lat","Year","VegC","LittVC",
						"LittSC","CwdC","MicroC","HumusC","SoilC","HarvSlowC","Total");
				else // GUESSN
					fprintf(out_cpool,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%10s\n","Lon","Lat","Year","VegC","LittVC",
						"LittSC","CwdC","MicroC","HumusC","SoilC","Total");
			}
		}

		if (out_firert) fprintf(out_firert,"%8s%8s%8s%8s\n","Lon","Lat","Year","FireRT");

		if (out_mnpp) fprintf(out_mnpp,lonlatyearstr,"Lon","Lat","Year");
		if (out_mlai) fprintf(out_mlai,lonlatyearstr,"Lon","Lat","Year");
		if (out_mgpp) fprintf(out_mgpp,lonlatyearstr,"Lon","Lat","Year");
		if (out_mra) fprintf(out_mra,lonlatyearstr,"Lon","Lat","Year");
		if (out_maet) fprintf(out_maet,lonlatyearstr,"Lon","Lat","Year");
		if (out_mpet) fprintf(out_mpet,lonlatyearstr,"Lon","Lat","Year");
		if (out_mevap) fprintf(out_mevap,lonlatyearstr,"Lon","Lat","Year");
		if (out_mintercep) fprintf(out_mintercep,lonlatyearstr,"Lon","Lat","Year");
		if (out_mrunoff) fprintf(out_mrunoff,lonlatyearstr,"Lon","Lat","Year");
		if (out_mrh) fprintf(out_mrh,lonlatyearstr,"Lon","Lat","Year");
		if (out_mnee) fprintf(out_mnee,lonlatyearstr,"Lon","Lat","Year");
		if (out_mwcont_upper) fprintf(out_mwcont_upper,lonlatyearstr,"Lon","Lat","Year");
		if (out_mwcont_lower) fprintf(out_mwcont_lower,lonlatyearstr,"Lon","Lat","Year");

		// DGPP
		if (out_dgpp) fprintf(out_dgpp,"%8s%8s%8s%8s","Lon","Lat","Year","Day");

		// GUESSN
		if (out_cton) fprintf(out_cton,lonlatyearstr,"Lon","Lat","Year");
		if (out_nmass) fprintf(out_nmass,lonlatyearstr,"Lon","Lat","Year");
		if (out_nsources) fprintf(out_nsources,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","Lon","Lat","Year","dep","fix","input","min",
				"imm","min-imm","Total","Ndemand");
		if (out_npool && ifcentury)
			fprintf(out_npool,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%10s\n","Lon","Lat","Year","VegN","LittVN",
			"LittSN","CwdN","MicroN","HumusN","SoilN","Total");
		if (out_nleach) fprintf(out_nleach,lonlatyearstr,"Lon","Lat","Year");
		if (out_nuptake) fprintf(out_nuptake,lonlatyearstr,"Lon","Lat","Year");
		if (out_anppn) fprintf(out_anppn,lonlatyearstr,"Lon","Lat","Year");
		if (out_vmaxnlim) fprintf(out_vmaxnlim,lonlatyearstr,"Lon","Lat","Year");
		if (out_nlim) fprintf(out_nlim,lonlatyearstr,"Lon","Lat","Year");
		if (out_canopyh) fprintf(out_canopyh,lonlatyearstr,"Lon","Lat","Year");
		// end GUESSN

		// Loop through PFT's and print PFT names as column labels

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pft=pftlist.getobj();
			if (out_cmass) fprintf(out_cmass,"%8s",(char*)pft.name);
			if (out_anpp) fprintf(out_anpp,"%8s",(char*)pft.name);
			if (out_lai) fprintf(out_lai,"%8s",(char*)pft.name);
			if (out_dens) fprintf(out_dens,"%8s",(char*)pft.name);
			// GUESSN
			if (out_cton) fprintf(out_cton,"%8s",(char*)pft.name);
			if (out_nmass) fprintf(out_nmass,"%8s",(char*)pft.name);
			if (out_nuptake) fprintf(out_nuptake,"%9s",(char*)pft.name);
			if (out_anppn) fprintf(out_anppn,"%9s",(char*)pft.name);
			if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8s",(char*)pft.name);
			if (out_nlim) fprintf(out_nlim,"%8s",(char*)pft.name);
			// end GUESSN
			
			// DGPP
			if (out_dgpp) fprintf(out_dgpp,"%9s",(char*)pft.name);
			pftlist.nextobj();
		}

		// Print labels for "Total" columns

		if (out_cmass) fprintf(out_cmass,"%8s","Total");
		if (out_anpp) fprintf(out_anpp,"%8s","Total");
		if (out_lai) fprintf(out_lai,"%8s","Total");
		if (out_runoff) fprintf(out_runoff,"%8s\n","Total");
		if (out_dens) fprintf(out_dens,"%8s\n","Total");

		// DGPP
		if (out_dgpp) fprintf(out_dgpp,"%9s%9s%9s\n","Total","Temp","N_dep");

		//TODO Fix these for landcover
		// GUESSN
		if (out_cton) fprintf(out_cton,"%8s\n","Total");
		if (out_nmass) fprintf(out_nmass,"%8s\n","Total");
		if (out_nleach) fprintf(out_nleach,"%8s%8s%8s\n","Min","Org","Total");	//(kgN/ha/yr)
		if (out_nuptake) fprintf(out_nuptake,"%9s\n","Total");
		if (out_anppn) fprintf(out_anppn,"%9s\n","Total");
		if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8s\n","Total");
		if (out_nlim) fprintf(out_nlim,"%8s\n","Total");
		if (out_allometry) fprintf(out_allometry,"\n");
		if (out_canopyh) fprintf(out_canopyh,"%8s\n","Height");
		// end GUESSN


		if (run_landcover) {
			xtring landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum", "Forest_sum", "Natural_sum", "Peatland_sum"};
			for (int i=0; i<NLANDCOVERTYPES; i++) {
				if(run[i]) {
					if (out_cmass) fprintf(out_cmass,"%13s",(char*)landcover_string[i]);
					if (out_anpp) fprintf(out_anpp,"%13s",(char*)landcover_string[i]);
					if (out_lai) fprintf(out_lai,"%13s",(char*)landcover_string[i]);
				}
			}
		}

		if (out_cmass) fprintf(out_cmass,"\n");
		if (out_anpp) fprintf(out_anpp,"\n");
		if (out_lai) fprintf(out_lai,"\n");



		// guess2008
		const char* monthstr = "%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n";
		const char* monthstr_long = "%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n";
		if (out_mnpp) fprintf(out_mnpp,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mlai) fprintf(out_mlai,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mgpp) fprintf(out_mgpp,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mra) fprintf(out_mra,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_maet) fprintf(out_maet,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mpet) fprintf(out_mpet,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mevap) fprintf(out_mevap,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mintercep) fprintf(out_mintercep,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mrunoff) fprintf(out_mrunoff,monthstr_long,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mrh) fprintf(out_mrh,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mnee) fprintf(out_mnee,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mwcont_upper) fprintf(out_mwcont_upper,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		if (out_mwcont_lower) fprintf(out_mwcont_lower,monthstr,"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");

		firstgrid=false;
	}

	// guess2008 - yearly output after spinup

	if (date.year>=nyear_spinup) {

		lon=gridlist.getobj().lon;
		lat=gridlist.getobj().lat;

		// Print longitude, latitude, year

		// guess2008
		const char* lonlatyeardatastr = "%8.1f%8.1f%8d"; // std CRU
		if (out_cmass) fprintf(out_cmass,lonlatyeardatastr,lon,lat,date.year);
		if (out_anpp) fprintf(out_anpp,lonlatyeardatastr,lon,lat,date.year);
		if (out_lai) fprintf(out_lai,lonlatyeardatastr,lon,lat,date.year);
		if (out_cflux) fprintf(out_cflux,lonlatyeardatastr,lon,lat,date.year);
		if (out_runoff) fprintf(out_runoff,lonlatyeardatastr,lon,lat,date.year);
		if (out_dens) fprintf(out_dens,lonlatyeardatastr,lon,lat,date.year);
		if (out_cpool) fprintf(out_cpool,lonlatyeardatastr,lon,lat,date.year);
		if (out_firert) fprintf(out_firert,lonlatyeardatastr,lon,lat,date.year);

		// GUESSN
		if (out_cton) fprintf(out_cton,lonlatyeardatastr,lon,lat,date.year);
		if (out_nmass) fprintf(out_nmass,lonlatyeardatastr,lon,lat,date.year);
		if (out_nsources) fprintf(out_nsources,lonlatyeardatastr,lon,lat,date.year);
		if (out_npool && ifcentury) fprintf(out_npool,lonlatyeardatastr,lon,lat,date.year);
		if (out_nleach) fprintf(out_nleach,lonlatyeardatastr,lon,lat,date.year);
		if (out_nuptake) fprintf(out_nuptake,lonlatyeardatastr,lon,lat,date.year);
		if (out_anppn) fprintf(out_anppn,lonlatyeardatastr,lon,lat,date.year);
		if (out_vmaxnlim) fprintf(out_vmaxnlim,lonlatyeardatastr,lon,lat,date.year);
		if (out_nlim) fprintf(out_nlim,lonlatyeardatastr,lon,lat,date.year);
		if (out_canopyh) fprintf(out_canopyh,lonlatyeardatastr,lon,lat,date.year);
		// end GUESSN

		if (out_mnpp) fprintf(out_mnpp,lonlatyeardatastr,lon,lat,date.year);
		if (out_mlai) fprintf(out_mlai,lonlatyeardatastr,lon,lat,date.year);
		if (out_mgpp) fprintf(out_mgpp,lonlatyeardatastr,lon,lat,date.year);
		if (out_mra) fprintf(out_mra,lonlatyeardatastr,lon,lat,date.year);
		if (out_maet) fprintf(out_maet,lonlatyeardatastr,lon,lat,date.year);
		if (out_mpet) fprintf(out_mpet,lonlatyeardatastr,lon,lat,date.year);
		if (out_mevap) fprintf(out_mevap,lonlatyeardatastr,lon,lat,date.year);
		if (out_mintercep) fprintf(out_mintercep,lonlatyeardatastr,lon,lat,date.year);
		if (out_mrunoff) fprintf(out_mrunoff,lonlatyeardatastr,lon,lat,date.year);
		if (out_mrh) fprintf(out_mrh,lonlatyeardatastr,lon,lat,date.year);
		if (out_mnee) fprintf(out_mnee,lonlatyeardatastr,lon,lat,date.year);
		if (out_mwcont_upper) fprintf(out_mwcont_upper,lonlatyeardatastr,lon,lat,date.year);
		if (out_mwcont_lower) fprintf(out_mwcont_lower,lonlatyeardatastr,lon,lat,date.year);

		// guess2008 - reset monthly average across patches each year
		for (m=0;m<12;m++)
			mnpp[m]=mlai[m]=mgpp[m]=mra[m]=maet[m]=mpet[m]=mevap[m]=mintercep[m]=mrunoff[m]=mrh[m]=mnee[m]=mwcont_upper[m]=mwcont_lower[m]=0.0;

		for (int day=0;day<365;day++)
			for (int ppfftt=0;ppfftt<20;ppfftt++)
				dgpp[day][ppfftt]=0.0;

		double landcover_cmass[NLANDCOVERTYPES]={0.0};
		double landcover_anpp[NLANDCOVERTYPES]={0.0};
		double landcover_lai[NLANDCOVERTYPES]={0.0};
		double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};

		double gcpft_cmass=0.0;
		double gcpft_anpp=0.0;
		double gcpft_lai=0.0;
		double gcpft_densindiv_total=0.0;
		double gcpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};

		double cmass_gridcell=0.0;
		double anpp_gridcell=0.0;
		double lai_gridcell=0.0;
		double runoff_gridcell=0.0;
		double dens_gridcell=0.0;
		double firert_gridcell=0.0;
		double canopyheight_gridcell=0.0;

		double standpft_cmass=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};

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
		double gcpft_nmass=0.0;		// sum/mean across patches for nitrogen biomass (kgN/m2)
		double gcpft_nuptake=0.0;	// sum across patches for nitrogen uptake (kgN/m2)
		double gcpft_anppn=0.0;		// sum across patches for nitrogen ANPP usage (kgN/m2)
		double gcpft_vmaxnlim=0.0;	// N limitation on vm
		double gcpft_nlim=0.0;	// N limitation on growth

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

		double surfsoillitterc=0.0;
		double surfsoillittern=0.0;
		double cwdc=0.0;
		double cwdn=0.0;
		double microc=0.0;
		double micron=0.0;
		double humusc=0.0;
		double humusn=0.0;
		double centuryc=0.0;
		double centuryn=0.0;
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

			// Sum C biomass, NPP and LAI across patches and PFTs
			gcpft_cmass=0.0;
			gcpft_anpp=0.0;
			gcpft_lai=0.0;
			gcpft_densindiv_total=0.0;

			// GUESSN
			gcpft_cmass_leaf=0.0;
			gcpft_nmass_leaf=0.0;
			gcpft_nmass=0.0;
			gcpft_nuptake=0.0;
			gcpft_anppn=0.0;
			gcpft_vmaxnlim=0.0;
			gcpft_nlim=0.0;

			// end GUESSN

			for (int day=0;day<365;day++)
				dgpp[day][pft.id]=0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				// Sum C biomass, NPP and LAI across patches and PFTs
				standpft_cmass=0.0;
				standpft_anpp=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total=0.0;
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

								// DGPP
								if (date.year>595)
									for (int day=0;day<365;day++)
										dgpp[day][pft.id]+=indiv.dassim[day]*indiv.limnfact;

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
										}
									}
								}

							}

						} // alive?
						vegetation.nextobj();
					}
					stand.nextobj();
				} // end of patch loop

				standpft_cmass/=(double)stand.nobj;
				standpft_anpp/=(double)stand.nobj;
				standpft_lai/=(double)stand.nobj;
				standpft_densindiv_total/=(double)stand.nobj;

				// GUESSN
				standpft_nuptake/=(double)stand.nobj;
				standpft_anppn/=(double)stand.nobj;
				standpft_nmass/=(double)stand.nobj;
				standpft_anpp_no_nlim/=(double)stand.nobj;

				if (!negligible(standpft_anpp_no_nlim))
					standpft_nlim=standpft_anpp/standpft_anpp_no_nlim;
				
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

				// DGPP
				if (date.year>595)
					for (int day=0;day<365;day++)
						dgpp[day][pft.id]/=(double)stand.nobj;

				//Update landcover totals
				landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
				landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
				landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();

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

				nmass_gridcell+=standpft_nmass*fraction_of_gridcell;
				nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
				anppn_gridcell+=standpft_anppn*fraction_of_gridcell;
				anpp_no_nlim_gridcell+=standpft_anpp_no_nlim*fraction_of_gridcell;

				if (out_canopyh)
					canopyheight_gridcell=canopy_height(stand);
				

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
			if (out_lai)	fprintf(out_lai,"%8.4f",gcpft_lai);
			if (out_dens)	fprintf(out_dens,"%8.4f",gcpft_densindiv_total);
			if (out_cmass)	fprintf(out_cmass,"%8.3f",gcpft_cmass);
			if (out_anpp)	fprintf(out_anpp,"%8.3f",gcpft_anpp);

			// GUESSN
			double gcpft_cton_leaf=0.0;
			if (gcpft_cmass_leaf>0.0) {
				gcpft_cton_leaf=gcpft_cmass_leaf/gcpft_nmass_leaf;
				gcpft_vmaxnlim/=gcpft_cmass_leaf;
			}
			
			if (out_cton)		fprintf(out_cton,"%8.3f",gcpft_cton_leaf);
			if (out_vmaxnlim)	fprintf(out_vmaxnlim,"%8.3f",gcpft_vmaxnlim);
			if (out_nlim)	fprintf(out_nlim,"%8.3f",gcpft_nlim);
			if (out_nmass)		fprintf(out_nmass,"%8.3f",gcpft_nmass);
			if (out_nuptake)	fprintf(out_nuptake,"%9.5f",gcpft_nuptake);
			if (out_anppn)		fprintf(out_anppn,"%9.5f",gcpft_anppn);
			
			// end GUESSN

			pftlist.nextobj();

		} // *** End of PFT loop ***

		if (date.year > 595) {
			for (int day=0;day<365;day++) {

				if (out_dgpp) fprintf(out_dgpp,lonlatyeardatastr,lon,lat,date.year+FIRSTHISTYEAR-nyear_spinup);
				int dateday=day+1;
				if (out_dgpp) fprintf(out_dgpp,"%9.0d",dateday);

				double total_dgpp=0.0;

				pftlist.firstobj();
				while (pftlist.isobj) {

					Pft& pft=pftlist.getobj();
	
					if (out_dgpp)	fprintf(out_dgpp,"%9.5f",dgpp[day][pft.id]);

					total_dgpp+=dgpp[day][pft.id];
	
					pftlist.nextobj();

				} // *** End of PFT loop ***

				if (out_dgpp)	fprintf(out_dgpp,"%9.5f%9.3f%9.5f\n",total_dgpp,dtemp[day],gridcell.climate.dndep[day]*10000.0);
			}
		}

		flux_veg=flux_soil=flux_fire=flux_est=flux_harvest=0.0;

		// guess2008 - carbon pools
		c_litter=c_fast=c_slow=c_harv_slow=0.0;

		// GUESSN
		surfsoillitterc=surfsoillittern=cwdc=cwdn=microc=micron=humusc=humusn=centuryc=centuryn=n_litter=n_harv_slow=0.0;
		andep_gridcell=anmin_gridcell=animm_gridcell=anfix_gridcell=nsupply_gridcell=ndemand_gridcell=0.0;
		n_org_leach_gridcell=n_min_leach_gridcell=0.0;

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

				double to_gridcell_average = stand.get_gridcell_fraction()/(double)stand.nobj;

				flux_veg+=patch.fluxes.acflux_veg*to_gridcell_average;
				flux_soil+=patch.fluxes.acflux_soil*to_gridcell_average;
				flux_fire+=patch.fluxes.acflux_fire*to_gridcell_average;
				flux_est+=patch.fluxes.acflux_est*to_gridcell_average;
				flux_harvest+=patch.fluxes.acflux_harvest*to_gridcell_average;

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

				runoff_gridcell+=patch.arunoff*to_gridcell_average;

				// Fire return time
				if (!iffire || patch.fireprob < 0.001)
					firert_gridcell+=1000.0/(double)stand.nobj; // Set a limit of 1000 years
				else
					firert_gridcell+=(1.0/patch.fireprob)/(double)stand.nobj;

				// GUESSN
				andep_gridcell+=patch.soil.ndep_annual/(double)stand.nobj*10000.0;	// convert from m2 to ha
				anmin_gridcell+=patch.soil.nmin_annual/(double)stand.nobj*10000.0;	// convert from m2 to ha
				animm_gridcell+=patch.soil.nimmob_annual/(double)stand.nobj*10000.0; // convert from m2 to ha
				anfix_gridcell+=patch.soil.N_fix/(double)stand.nobj*10000.0;		// convert from m2 to ha
				n_min_leach_gridcell+=patch.soil.n_min_leach_annual/(double)stand.nobj*10000.0;	// convert from m2 to ha
				n_org_leach_gridcell+=patch.soil.n_org_leach_annual/(double)stand.nobj*10000.0;	// convert from m2 to ha
				nsupply_gridcell+=patch.nsupply/(double)stand.nobj*10000.0;			// convert from m2 to ha
				ndemand_gridcell+=patch.ndemand/(double)stand.nobj*10000.0;			// convert from m2 to ha

				for (int r=0;r<NSOMPOOL;r++) {
					if (patch.soil.sompool[r].nmass > 0.0) {
						if(r==SURFMETA||r==SURFSTRUCT||r==SOILMETA||r==SOILSTRUCT){
							surfsoillitterc+=patch.soil.sompool[r].cmass/(double)stand.nobj;
							surfsoillittern+=patch.soil.sompool[r].nmass/(double)stand.nobj;
						}
						else if (r==SURFCWD) {
							cwdc+=patch.soil.sompool[r].cmass/(double)stand.nobj;
							cwdn+=patch.soil.sompool[r].nmass/(double)stand.nobj;
						}
						else if (r==SURFMICRO||r==SOILMICRO) {
							microc+=patch.soil.sompool[r].cmass/(double)stand.nobj;
							micron+=patch.soil.sompool[r].nmass/(double)stand.nobj;
						}
						else if (r==SURFHUMUS){
							humusc+=patch.soil.sompool[r].cmass/(double)stand.nobj;
							humusn+=patch.soil.sompool[r].nmass/(double)stand.nobj;
						}
						
						centuryc+=patch.soil.sompool[r].cmass/(double)stand.nobj;
						centuryn+=patch.soil.sompool[r].nmass/(double)stand.nobj;
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
		double testmnpp = 0.0;
		double testmlai = 0.0;

		for (m=0;m<12;m++) {
			mnpp[m] = mgpp[m]-mra[m];
			mnee[m] = mnpp[m]-mrh[m];
			testmnpp += mnpp[m];
			testmlai += mlai[m]/12.0;
		}

		// Print gridcell totals to files

		double cton_leaf_gridcell=0.0;
		if (cmass_leaf_gridcell>0.0) {
			cton_leaf_gridcell=cmass_leaf_gridcell/nmass_leaf_gridcell;
			vmaxnlim_gridcell/=cmass_leaf_gridcell;
		}

		if (out_cmass) fprintf(out_cmass,"%8.3f",cmass_gridcell);
		if (out_anpp) fprintf(out_anpp,"%8.3f",anpp_gridcell);
		if (out_lai) fprintf(out_lai,"%8.4f",lai_gridcell);
		if (out_runoff) fprintf(out_runoff,"%8.1f",runoff_gridcell);
		if (out_dens) fprintf(out_dens,"%8.4f",dens_gridcell);
		if (out_firert) fprintf(out_firert,"%8.1f",firert_gridcell);

		//GUESSN
		double nlim_grid=0.0;
		if (!negligible(anpp_no_nlim_gridcell))
			nlim_grid=anpp_gridcell/anpp_no_nlim_gridcell;

		if (out_nmass) fprintf(out_nmass,"%8.3f",nmass_gridcell);
		if (out_anppn) fprintf(out_anppn,"%8.3f",anppn_gridcell);
		if (out_cton) fprintf(out_cton,"%8.4f",cton_leaf_gridcell);
		if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8.3f",vmaxnlim_gridcell);
		if (out_nlim) fprintf(out_nlim,"%8.3f",nlim_grid);
		if (out_canopyh) fprintf(out_canopyh,"%8.1f",canopyheight_gridcell);
		if (out_nuptake) fprintf(out_nuptake,"%9.5f",nuptake_gridcell);
		if (out_nsources) fprintf(out_nsources,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
			andep_gridcell,anfix_gridcell,andep_gridcell+anfix_gridcell,anmin_gridcell,animm_gridcell,
			anmin_gridcell-animm_gridcell,andep_gridcell+anmin_gridcell-animm_gridcell+anfix_gridcell,ndemand_gridcell);
		if (out_nleach)	fprintf(out_nleach,"%8.3f%8.3f%8.3f\n",n_min_leach_gridcell,n_org_leach_gridcell,
			n_min_leach_gridcell+n_org_leach_gridcell);

		// end GUESSN

		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {
					if (out_cmass)
						fprintf(out_cmass,"%13.3f", landcover_cmass[i]);
					if (out_anpp)
						fprintf(out_anpp,"%13.3f", landcover_anpp[i]);
					if (out_lai)
						fprintf(out_lai,"%13.3f", landcover_lai[i]);
					if (out_dens)
						fprintf(out_dens,"%13.4f", landcover_densindiv_total[i]);
					//GUESSN
					if (out_nmass)
						fprintf(out_nmass,"%13.3f", landcover_nmass[i]);
					if (out_anppn)
						fprintf(out_anppn,"%13.3f", landcover_anppn[i]);

					double landcover_cton_leaf=0.0;
					if (landcover_cmass_leaf[i]>0.0) {
						landcover_cton_leaf=landcover_cmass_leaf[i]/landcover_nmass_leaf[i];
						landcover_vmaxnlim[i]/=landcover_cmass_leaf[i];
					}
					if (out_cton)
						fprintf(out_cton,"%13.3f", landcover_cton_leaf);
					if (out_vmaxnlim)
						fprintf(out_vmaxnlim,"%13.4f", landcover_vmaxnlim[i]);
					if (out_nlim)
						fprintf(out_nlim,"%13.4f", landcover_nlim[i]);
					if (out_nuptake)
						fprintf(out_nuptake,"%13.4f", landcover_nuptake[i]);
					// end GUESSN
				}
			}
		}

		if (out_cmass) fprintf(out_cmass,"\n");
		if (out_anpp) fprintf(out_anpp,"\n");
		if (out_lai) fprintf(out_lai,"\n");
		if (out_runoff) fprintf(out_runoff,"\n");
		if (out_dens) fprintf(out_dens, "\n");
		if (out_firert) fprintf(out_firert, "\n");
		//GUESSN
		if (out_nmass) fprintf(out_nmass,"\n");
		if (out_anppn) fprintf(out_anppn,"\n");
		if (out_cton) fprintf(out_cton,"\n");
		if (out_vmaxnlim) fprintf(out_vmaxnlim,"\n");
		if (out_nlim) fprintf(out_nlim,"\n");
		if (out_canopyh) fprintf(out_canopyh,"\n");
		if (out_nuptake) fprintf(out_nuptake, "\n");
		//end GUESSN
		// Print monthly output variables
		for (m=0;m<12;m++) {

			if (out_mnpp) fprintf(out_mnpp,"%8.3f",mnpp[m]);
			if (out_mlai) fprintf(out_mlai,"%8.3f",mlai[m]);
			if (out_mgpp) fprintf(out_mgpp,"%8.3f",mgpp[m]);
			if (out_mra) fprintf(out_mra,"%8.3f",mra[m]);
			if (out_maet) fprintf(out_maet,"%8.3f",maet[m]);
			if (out_mpet) fprintf(out_mpet,"%8.3f",mpet[m]);
			if (out_mevap) fprintf(out_mevap,"%8.3f",mevap[m]);
			if (out_mintercep) fprintf(out_mintercep,"%8.3f",mintercep[m]);
			if (out_mrunoff) fprintf(out_mrunoff,"%10.3f",mrunoff[m]);
			if (out_mrh) fprintf(out_mrh,"%8.3f",mrh[m]);
			if (out_mnee) fprintf(out_mnee,"%8.3f",mnee[m]);
			if (out_mwcont_upper) fprintf(out_mwcont_upper,"%8.3f",mwcont_upper[m]);
			if (out_mwcont_lower) fprintf(out_mwcont_lower,"%8.3f",mwcont_lower[m]);

			if (m==11) {
				if (out_mnpp) fprintf(out_mnpp,"\n");
				if (out_mlai) fprintf(out_mlai,"\n");
				if (out_mgpp) fprintf(out_mgpp,"\n");
				if (out_mra) fprintf(out_mra,"\n");
				if (out_maet) fprintf(out_maet,"\n");
				if (out_mpet) fprintf(out_mpet,"\n");
				if (out_mevap) fprintf(out_mevap,"\n");
				if (out_mintercep) fprintf(out_mintercep,"\n");
				if (out_mrunoff) fprintf(out_mrunoff,"\n");
				if (out_mrh) fprintf(out_mrh,"\n");
				if (out_mnee) fprintf(out_mnee,"\n");
				if (out_mwcont_upper) fprintf(out_mwcont_upper,"\n");
				if (out_mwcont_lower) fprintf(out_mwcont_lower,"\n");
			}

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

		if (out_cflux) {
			if(run_landcover)
				fprintf(out_cflux,"%8.3f%8.3f%8.3f%8.3f%9.3f%10.5f\n",flux_veg,flux_soil,flux_fire,
					flux_est,flux_harvest,flux_veg+flux_soil+flux_fire+flux_est+flux_harvest);
			else
				fprintf(out_cflux,"%8.3f%8.3f%8.3f%8.3f%10.5f\n",flux_veg,flux_soil,flux_fire,
					flux_est,flux_veg+flux_soil+flux_fire+flux_est);
		}



		// guess2008 - output carbon pools


		if (out_cpool) {
			if (!ifcentury){
				if(run_landcover && ifslowharvestpool)
					fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%10.3f%10.4f\n",cmass_gridcell,c_litter,c_fast,
						c_slow,c_harv_slow,cmass_gridcell+c_litter+c_fast+c_slow+c_harv_slow);
				else
					fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%10.4f\n",cmass_gridcell,c_litter,c_fast,
						c_slow,cmass_gridcell+c_litter+c_fast+c_slow);
			}
			else{
				if(run_landcover && ifslowharvestpool)
					fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10.3f\n",cmass_gridcell,c_litter,
						surfsoillitterc,cwdc,microc,humusc,c_litter+centuryc,c_harv_slow,
						cmass_gridcell+c_litter+centuryc+c_harv_slow);
				else // GUESSN
					fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10.3f\n",cmass_gridcell,c_litter,
						surfsoillitterc,cwdc,microc,humusc,c_litter+centuryc,
						cmass_gridcell+c_litter+centuryc);
			}
		}

		// GUESSN
		if (out_npool && ifcentury) {
			if(run_landcover && ifslowharvestpool) {
					fprintf(out_npool,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10.3f\n",nmass_gridcell,n_litter,
					surfsoillittern,cwdn,micron,humusn,n_litter+centuryn,n_harv_slow,
					nmass_gridcell+n_litter+centuryn+n_harv_slow);
			}
			else {
				fprintf(out_npool,"%8.3f%8.3f%8.4f%8.4f%8.4f%8.3f%8.3f%10.3f\n",nmass_gridcell,n_litter,
					surfsoillittern,cwdn,micron,humusn,n_litter+centuryn,
					nmass_gridcell+n_litter+centuryn);
			}
		}
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
										plot("age_structure",pft.name,c*estinterval+estinterval/2,gcpft_densindiv_ageclass[c]/(double)stand.nobj);
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

	// DESCRIPTION
	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Close output files if open

	if (annual_output) {
		if (out_cmass) fclose(out_cmass);
		if (out_anpp) fclose(out_anpp);
		if (out_lai) fclose(out_lai);
		if (out_cflux) fclose(out_cflux);
		if (out_runoff) fclose(out_runoff);
		if (out_dens) fclose(out_dens);
		if (out_cpool) fclose(out_cpool);
		if (out_firert) fclose(out_firert);

		if (out_mnpp) fclose(out_mnpp);
		if (out_mlai) fclose(out_mlai);
		if (out_mgpp) fclose(out_mgpp);
		if (out_mra) fclose(out_mra);
		if (out_maet) fclose(out_maet);
		if (out_mpet) fclose(out_mpet);
		if (out_mevap) fclose(out_mevap);
		if (out_mrunoff) fclose(out_mrunoff);
		if (out_mintercep) fclose(out_mintercep);
		if (out_mrh) fclose(out_mrh);
		if (out_mnee) fclose(out_mnee);
		if (out_mwcont_upper) fclose(out_mwcont_upper);
		if (out_mwcont_lower) fclose(out_mwcont_lower);
		// GUESSN
		if (out_cton) fclose(out_cton);
		if (out_nmass) fclose(out_nmass);
		if (out_nsources) fclose(out_nsources);
		if (out_npool && ifcentury) fclose(out_npool);
		if (out_nleach) fclose(out_nleach);
		if (out_nuptake) fclose(out_nuptake);
		if (out_anppn) fclose(out_anppn);
		if (out_vmaxnlim) fclose(out_vmaxnlim);
		if (out_nlim) fclose(out_nlim);
		if (out_canopyh) fclose(out_canopyh);
		// end GUESSN

		// DGPP
		if (out_dgpp) fclose(out_dgpp);
	}

	// Clean up

	gridlist.killall();
}



#endif // USE_CRU


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Galloway, J. N., F. J. Dentener, D. G. Capone, E. W. Boyer, R. W. Howarth, S. P. Seitzinger,
//   G. P. Asner, C. Cleveland, P. Green, E. Holland, D. M. Karl, A. F. Michaels, J. H. Porter, 
//   A. Townsend, and C. Vï¿½rï¿½smarty. 2004.
//   Nitrogen Cycles: Past, Present and Future. Biogeochemistry 70: 153-226.
// Nakai, T., Sumida, A., Kodama, Y., Hara, T., Ohta, T. (2010). A comparison between
//   various definitions of forest stand height and aerodynamic canopy height.
//   Agricultural and Forest Meteorology, 150(9), 1225-1233
