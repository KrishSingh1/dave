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

// ABOUT THIS I/O MODULE:
// This is a demonstration I/O module. It is compatible with the input data files
// distributed with LPJ-GUESS (in the data directory).

#include "config.h"

#ifdef USE_DEMO_IO

#include "guessio.h"

#include "driver.h"
#include <plib.h>
#include <stdio.h>


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
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_LANDCOVER,CB_PHENOLOGY,CB_PATHWAY,	
	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM};


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH FILE SCOPE

Paramlist param;

xtring title; // Title for this run
int nyear; // number of simulation years to run after spinup

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
	run_landcover = false;

	// guess2008 - initialise filenames here
	outputdirectory = "";
	file_cmass=file_anpp=file_lai=file_cflux=file_dens=file_runoff="";
	file_mnpp=file_mlai=file_maet=file_mpet=file_mevap=file_mrunoff=file_mintercep=file_mrh="";
	file_mgpp=file_mra=file_mnee=file_mwcont_upper=file_mwcont_lower="";
	file_cpool=file_firert=file_speciesheights="";

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
		declareitem("nyear",&nyear,1,10000,1,CB_NONE,"Number of simulation years to run after spinup");
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

		// guess2008 - new options
		declareitem("ifsmoothgreffmort",&ifsmoothgreffmort,1,CB_NONE,
			"Whether to vary mort_greff smoothly with growth efficiency (0,1)");
		declareitem("ifdroughtlimitedestab",&ifdroughtlimitedestab,1,CB_NONE,
			"Whether establishment drought limited (0,1)");
		declareitem("ifrainonwetdaysonly",&ifrainonwetdaysonly,1,CB_NONE,
			"Whether it rains on wet days only (1), or a little every day (0);");
		declareitem("ifspeciesspecificwateruptake",&ifspeciesspecificwateruptake,1,CB_NONE,
			"Whether or not there is species specific soil water uptake (0,1)");

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
		declareitem("cton_leaf",&ppft->cton_leaf,1.0,1.0e4,1,CB_NONE,
			"Leaf C:N mass ratio");
		declareitem("cton_root",&ppft->cton_root,1.0,1.0e4,1,CB_NONE,
			"Fine root C:N mass ratio");
		declareitem("cton_sap",&ppft->cton_sap,1.0,1.0e4,1,CB_NONE,
			"Sapwood C:N mass ratio");
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
		declareitem("res_outtake",&ppft->res_outtake,0.0,1.0,1,CB_NONE,"´Fraction of residue outtake at harvest");

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
		if (!itemparsed("nyear")) badins("nyear");
		if (!itemparsed("nyear_spinup")) badins("nyear_spinup");
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("ifdailynpp")) badins("ifdailynpp");
		if (!itemparsed("ifdailydecomp")) badins("ifdailydecomp");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");


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
		if (!itemparsed("cton_leaf")) badins("cton_leaf");
		if (!itemparsed("cton_root")) badins("cton_root");
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
			if (!itemparsed("cton_sap")) badins("cton_sap");
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
//   if (date.day==0 && date.year==nyear) return false; // guess2008
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

// File names for temperature, precipitation, sunshine and soil code driver files
xtring file_temp,file_prec,file_sun,file_soil;

// Open streams - guess2008
FILE *out_cmass,*out_anpp,*out_lai,*out_cflux,*out_cpool,*out_runoff,*out_dens;
FILE *out_mnpp,*out_mlai,*out_mgpp,*out_mra,*out_maet,*out_mpet,*out_mevap,*out_mrunoff,*out_mintercep,*out_mrh;
FILE *out_mnee,*out_mwcont_upper,*out_mwcont_lower; 
FILE *out_firert,*out_speciesheights; 



// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

double co2; // atmospheric CO2 concentration (ppmv) (read from ins file)

// Daily temperature, precipitation and sunshine for one year
double dtemp[365],dprec[365],dsun[365];

double cpool_sum;
	// sum of all C pools (including fireC) this/last year (for current gridcell)

// LPJ soil code
int soilcode;

bool annual_output;
	// whether output should occur each simulation year (true) or at end of simulation
	// for each grid cell only (false)


//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
//TimeDataD input code may be put here
TimeDataD LUdata(LOCAL_YEARLY);
TimeDataD Peatdata;
#endif
xtring file_lu, file_peat;
const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise
//

void readenv(Coord coord) {

	// Searches for environmental data in driver temperature, precipitation,
	// sunshine and soil code files for the grid cell whose coordinates are given by
	// 'coord'. Data are written to arrays mtemp, mprec, msun and the variable
	// soilcode, which are defined as global variables in this file

	// The temperature, precipitation and sunshine files (Cramer & Leemans,
	// unpublished) should be in ASCII text format and contain one-line records for the
	// climate (mean monthly temperature, mean monthly percentage sunshine or total
	// monthly precipitation) of a particular 0.5 x 0.5 degree grid cell. Elevation is
	// also included in each grid cell record.

	// The following sample record from the temperature file:
	//   " -4400 8300 293-298-311-316-239-105  -7  26  -3 -91-184-239-277"
	// corresponds to the following data:
	//   longitude 44 deg (-=W)
	//   latitude 83 deg (+=N)
	//   elevation 293 m
	//   mean monthly temperatures (deg C) -29.8 (Jan), -31.1 (Feb), ..., -27.7 (Dec)

	// The following sample record from the precipitation file:
	//   " 12750 -200 223 190 165 168 239 415 465 486 339 218 162 149 180"
	// corresponds to the following data:
	//   longitude 127.5 deg (+=E)
	//   latitude 20 deg (-=S)
	//   elevation 223 m
	//   monthly precipitation sum (mm) 190 (Jan), 165 (Feb), ..., 180 (Dec)

	// The following sample record from the sunshine file:
	//   "  2600 7000 293  0 20 38 37 31 28 28 25 21 17  9  7"
	// corresponds to the following data:
	//   longitude 26 deg (+=E)
	//   latitude 70 deg (+=N)
	//   elevation 293 m
	//   monthly mean %age of full sunshine 0 (Jan), 20 (Feb), ..., 7 (Dec)

	// The LPJ soil code file is in ASCII text format and contains one-line records for
	// each grid cell in the form:
	//   <lon> <lat> <soilcode>
	// where <lon>      = longitude as a floating point number (-=W, +=E)
	//       <lat>      = latitude as a floating point number (-=S, +=N)
	//       <soilcode> = integer in the range 0 (no soil) to 9 (see function
	//                    soilparameters in driver module)
	// The fields in each record are separated by spaces

	FILE* in;
	bool foundgrid;
	double dlon,dlat;
	int elev;
	double mtemp[12]; // monthly mean temperature (deg C)
	double mprec[12]; // monthly precipitation sum (mm)
	double msun[12]; // monthly mean percentage sunshine values

	double mwet[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // number of rain days per month


	// Search for record for this grid cell in temperature file

	in=fopen(file_temp,"r");
	if (!in) fail("readenv: could not open %s for input",(char*)file_temp);

	foundgrid=false;
	while (!feof(in) && !foundgrid) {

		// Read next record in file
		readfor(in,"f6.2,f5.2,i4,12f4.1",&dlon,&dlat,&elev,mtemp);
		
		if (equal(coord.lon,dlon) && equal(coord.lat,dlat)) foundgrid=true;
	}

	if (!foundgrid) fail("readenv: could not find record for (%g,%g) in %s",
		coord.lon,coord.lat,(char*)file_temp);

	fclose(in);

	// Search for record for this grid cell in precipitation file

	in=fopen(file_prec,"r");
	if (!in) fail("readenv: could not open %s for input",(char*)file_prec);

	foundgrid=false;
	while (!feof(in) && !foundgrid) {

		// Read next record in file
		readfor(in,"f6.2,f5.2,i4,12f4",&dlon,&dlat,&elev,mprec);
		
		if (equal(coord.lon,dlon) && equal(coord.lat,dlat)) foundgrid=true;
	}

	if (!foundgrid) fail("readenv: could not find record for (%g,%g) in %s",
		coord.lon,coord.lat,(char*)file_prec);

	fclose(in);

	// Search for record for this grid cell in sunshine file

	in=fopen(file_sun,"r");
	if (!in) fail("readenv: could not open %s for input",(char*)file_sun);

	foundgrid=false;
	while (!feof(in) && !foundgrid) {

		// Read next record in file
		readfor(in,"f6.2,f5.2,i4,12f3",&dlon,&dlat,&elev,msun);
		
		if (equal(coord.lon,dlon) && equal(coord.lat,dlat)) foundgrid=true;
	}

	if (!foundgrid) fail("readenv: could not find record for (%g,%g) in %s",
		coord.lon,coord.lat,(char*)file_sun);

	fclose(in);

	// Search for record for this grid cell in soil code file

	in=fopen(file_soil,"r");
	if (!in) fail("readenv: could not open %s for input",(char*)file_soil);

	foundgrid=false;
	while (!feof(in) && !foundgrid) {

		// Read next record in file
		readfor(in,"f,f,i",&dlon,&dlat,&soilcode);
		
		if (equal(coord.lon,dlon) && equal(coord.lat,dlat)) foundgrid=true;
	}

	if (!foundgrid) fail("readenv: could not find record for (%g,%g) in %s",
		coord.lon,coord.lat,(char*)file_soil);

	fclose(in);

	// Interpolate monthly values for environmental drivers to daily values
	// (relevant daily values will be sent to the framework each simulation
	// day in function getclimate, below)

	interp_climate(mtemp,mprec,msun,dtemp,dprec,dsun);

	// Recalculate precipitation values using weather generator
	// (from Dieter Gerten 021121)

	prdaily(mprec,dprec,mwet);
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

		else if (!readins(insfilename,pftlist))
			abort=true;
	}
	else abort=true;

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

	// Retrieve specified CO2 value as read from ins file
	co2=param["co2"].num;

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


	// Retrieve input file names as read from ins file

	file_temp=param["file_temp"].str;
	file_prec=param["file_prec"].str;
	file_sun=param["file_sun"].str;
	file_soil=param["file_soil"].str;

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
	
	if (file_speciesheights!="") {
		file_speciesheights = outputdirectory + file_speciesheights;
		out_speciesheights=fopen(file_speciesheights,"w");
		if (!file_speciesheights) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_speciesheights);
	}
	else out_speciesheights=NULL;

	if (file_runoff!="") {
		file_runoff = outputdirectory + file_runoff;
		out_runoff=fopen(file_runoff,"w");
		if (!out_runoff) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_runoff);
	}
	else out_runoff=NULL;


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
			if(!LUdata.Load(c))		//Load area fraction data from Bondeau input file to data object
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

/// Called by the framework at the start of the simulation for a particular grid cell
bool getgridcell(Gridcell& gridcell) 
{
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
	
	// guess2008 - elevation
	int elevation;

	bool gridfound=false;
	bool LUerror=false;

	// guess2008 - to ensure an identical random number sequence for each stand.
	setseed(12345678);

	if (firstgrid) {
		gridlist.firstobj();
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		while(!gridfound)
		{

			// Retrieve coordinate of next grid cell from linked list
			Coord& c=gridlist.getobj();

			// Load environmental data for this grid cell from files
			// (these will be the same for every year of the simulation, but must be sent
			// anew to the framework each year in function getclimate, below)

			if(run_landcover)
				LUerror=loadlandcover(gridcell, c);

			if (!LUerror) {
				readenv(c);
				gridfound=true;
			}
			else
				gridlist.nextobj();
		}

		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n");
		
		// Tell framework the latitude of this grid cell
		gridcell.climate.lat=gridlist.getobj().lat;
		
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		gridcell.climate.instype=SUNSHINE;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

///	Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
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
						dprintf("WARNING ! landcover fraction sum is %4.2f for year %d\n", sum_tot, year);
						dprintf("Rescaling landcover fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year);
					}
				}
				else				//added scaling to sum=1.0 (sum often !=1.0)
					dprintf("Rescaling landcover fractions year %d ! (sum is within 0.99-1.01)\n", date.year);

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


/// Called by the framework each simulation day before any process modelling is performed for this day
/** Obtains climate data (including atmospheric CO2 and insolation) for this day. */
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

	Climate& climate=gridcell.climate;

	// Send environmental values for today to framework

	climate.co2=co2;
	climate.temp=dtemp[date.day];
	climate.prec=dprec[date.day];
	climate.insol=dsun[date.day];

	// First day of year only ...

	if (date.day==0) {

		// Return false if last year was the last for the simulation
		if (date.year==nyear_spinup+nyear) return false;

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+nyear)
				+date.year)/(double)(ngridcell*(nyear_spinup+nyear));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


/// Called by the framework at the end of the last day of each simulation year
void outannual(Gridcell& gridcell,Pftlist& pftlist) {

	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int p,c,m,nclass;
	double flux_veg,flux_soil,flux_fire,flux_est,flux_harvest;
	double c_litter,c_fast,c_slow,c_harv_slow; 

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
			if(run_landcover && ifslowharvestpool)
				fprintf(out_cpool,"%8s%8s%8s%8s%8s%8s%8s%10s%10s\n","Lon","Lat","Year","VegC","LittC",
					"SoilfC","SoilsC", "HarvSlowC","Total");
			else
				fprintf(out_cpool,lonlatyearstr_extended,"Lon","Lat","Year","VegC","LittC",
					"SoilfC","SoilsC","Total");
		}

		if (out_firert) fprintf(out_firert,"%8s%8s%8s%8s\n","Lon","Lat","Year","FireRT");
		if (out_speciesheights) fprintf(out_speciesheights,lonlatyearstr,"Lon","Lat","Year");

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
		


		// Loop through PFT's and print PFT names as column labels

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pft=pftlist.getobj();
			if (out_cmass) fprintf(out_cmass,"%8s",(char*)pft.name);
			if (out_anpp) fprintf(out_anpp,"%8s",(char*)pft.name);
			if (out_lai) fprintf(out_lai,"%8s",(char*)pft.name);
			if (out_dens) fprintf(out_dens,"%8s",(char*)pft.name);
			if (out_speciesheights) fprintf(out_speciesheights,"%8s",(char*)pft.name);
			pftlist.nextobj();
		}

		// Print labels for "Total" columns

		if (out_cmass) fprintf(out_cmass,"%8s","Total");
		if (out_anpp) fprintf(out_anpp,"%8s","Total");
		if (out_lai) fprintf(out_lai,"%8s","Total");
		if (out_runoff) fprintf(out_runoff,"%8s\n","Total");
		if (out_dens) fprintf(out_dens,"%8s\n","Total");

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
		if (out_speciesheights) fprintf(out_speciesheights, "\n");


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
		if (out_speciesheights) fprintf(out_speciesheights,lonlatyeardatastr,lon,lat,date.year);

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

		double standpft_cmass=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};


		// *** Loop through PFTs ***

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Gridcellpft& gridcellpft=gridcell.pft[pft.id];

			// Sum C biomass, NPP and LAI across patches and PFTs		
			gcpft_cmass=0.0;
			gcpft_anpp=0.0;
			gcpft_lai=0.0;
			gcpft_densindiv_total=0.0;		

			double heightindiv_total = 0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				// Sum C biomass, NPP and LAI across patches and PFTs
				standpft_cmass=0.0;
				standpft_anpp=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total = 0.0;

				// Initialise age structure array

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++)
						standpft_densindiv_ageclass[c]=0.0;
		
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

				standpft_cmass/=(double)stand.nobj;
				standpft_anpp/=(double)stand.nobj;
				standpft_lai/=(double)stand.nobj;
				standpft_densindiv_total/=(double)stand.nobj;

				heightindiv_total/=(double)stand.nobj;

				//Update landcover totals
				landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
				landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
				landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();

				//Update pft totals
				gcpft_cmass+=standpft_cmass;
				gcpft_anpp+=standpft_anpp;
				gcpft_lai+=standpft_lai;
				gcpft_densindiv_total+=standpft_densindiv_total;

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++)
						gcpft_densindiv_ageclass[c]+=standpft_densindiv_ageclass[c];

				// Update gridcell totals
				double fraction_of_gridcell = stand.get_gridcell_fraction();
				cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
				anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
				lai_gridcell+=standpft_lai*fraction_of_gridcell;
				dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
			
				// Graphical output every 10 years
				// (Windows shell only - "plot" statements have no effect otherwise)
				if (!(date.year%10)) {
					plot("cmass",pft.name,date.year,gcpft_cmass);
					plot("anpp",pft.name,date.year,gcpft_anpp);
					plot("lai",pft.name,date.year,gcpft_lai);
				}
				gridcell.nextobj();
			}//End of loop through stands

			// Print PFT sums to files
			if (out_lai)
				fprintf(out_lai,"%8.4f",gcpft_lai);

			if (out_dens) fprintf(out_dens,"%8.4f",gcpft_densindiv_total);

			if (out_cmass)
				fprintf(out_cmass,"%8.3f",gcpft_cmass);

			if (out_anpp)
				fprintf(out_anpp,"%8.3f",gcpft_anpp);

			// print species heights
			double height = 0.0;
			if (gcpft_densindiv_total > 0.0)
				height = heightindiv_total/gcpft_densindiv_total;
			
			if (out_speciesheights)
				fprintf(out_speciesheights,"%8.2f",height);

			pftlist.nextobj();
		
		} // *** End of PFT loop ***


		flux_veg=flux_soil=flux_fire=flux_est=flux_harvest=0.0;

		// guess2008 - carbon pools
		c_litter=c_fast=c_slow=c_harv_slow=0.0;

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
				}

				//Sum slow pools of harvested products
				if(run_landcover && ifslowharvestpool)
				{
					for (int q=0;q<npft;q++) 
					{
						Patchpft& patchpft=patch.pft[q];
						c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
					}
				}

				runoff_gridcell+=patch.arunoff*to_gridcell_average;
	
				// Fire return time
				if (!iffire || patch.fireprob < 0.001)
					firert_gridcell+=1000.0/(double)stand.nobj; // Set a limit of 1000 years
				else	
					firert_gridcell+=(1.0/patch.fireprob)/(double)stand.nobj;


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

		if (out_cmass) fprintf(out_cmass,"%8.3f",cmass_gridcell);
		if (out_anpp) fprintf(out_anpp,"%8.3f",anpp_gridcell);
		if (out_lai) fprintf(out_lai,"%8.4f",lai_gridcell);
		if (out_runoff) fprintf(out_runoff,"%8.1f",runoff_gridcell);
		if (out_dens) fprintf(out_dens,"%8.4f",dens_gridcell);
		if (out_firert) fprintf(out_firert,"%8.1f",firert_gridcell);

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
				}
			}
		}

		if (out_cmass) fprintf(out_cmass,"\n");
		if (out_anpp) fprintf(out_anpp,"\n");
		if (out_lai) fprintf(out_lai,"\n");
		if (out_runoff) fprintf(out_runoff,"\n");
		if (out_dens) fprintf(out_dens, "\n");
		if (out_firert) fprintf(out_firert, "\n");
		if (out_speciesheights) fprintf(out_speciesheights, "\n");

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

		if (!(date.year%10)) {
			gridcell.firstobj();
			if(gridcell.isobj)	//Fixed bug here if no stands were present.
			{
				Stand& stand=gridcell.getobj();
				plot("fluxes","flux_veg",date.year,flux_veg);
				plot("fluxes","flux_soil",date.year,flux_soil);
				plot("fluxes","flux_fire",date.year,flux_fire);
				plot("fluxes","flux_est",date.year,flux_est);
				plot("fluxes","NEE",date.year,flux_veg+flux_soil+flux_fire+flux_est);
				plot("soilc","slow",date.year,stand[0].soil.cpool_slow);
				plot("soilc","fast",date.year,stand[0].soil.cpool_fast);
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
			if(run_landcover && ifslowharvestpool)
				fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%10.3f%10.4f\n",cmass_gridcell,c_litter,c_fast,
					c_slow,c_harv_slow,cmass_gridcell+c_litter+c_fast+c_slow+c_harv_slow);
			else
				fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%10.4f\n",cmass_gridcell,c_litter,c_fast,
					c_slow,cmass_gridcell+c_litter+c_fast+c_slow);
		}


		// Output of age structure (Windows shell only - no effect otherwise)

		if (vegmode==COHORT || vegmode==INDIVIDUAL) {

			if (!(date.year%20) && date.year<2000) {
			
				resetwindow("age_structure");

				pftlist.firstobj();
				while (pftlist.isobj) {
					Pft& pft=pftlist.getobj();

					if (pft.lifeform==TREE) {

						Gridcellpft& gridcellpft=gridcell.pft[pft.id];

						for (c=0;c<nclass;c++)
							plot("age_structure",pft.name,
								c*estinterval+estinterval/2,
								gcpft_densindiv_ageclass[c]/(double)npatch);
					}
					
					pftlist.nextobj();
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
		if (out_speciesheights) fclose(out_speciesheights);

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
	}

	// Clean up

	gridlist.killall();
}

#endif // USE_DEMO_IO
