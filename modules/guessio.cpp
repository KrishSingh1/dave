///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                LPJ-GUESS input/output module with input from instruction
//                        script
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
//                        Includes Dieter G:s latest updates 021121
//                        Updated 20050125: last line in output files ends in newline
// Header file name:      guessio.h
// Source code file name: guessio.cpp
// Written by:            Ben Smith
// Version dated:         2002-11-22/2005-01-25
//
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
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_PHENOLOGY,CB_PATHWAY,
	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM};


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH FILE SCOPE

Paramlist param;

xtring title; // Title for this run
int nyear; // number of simulation years

Pftlist* ppftlist; // pointer to PFT list
Pft* ppft; // pointer to Pft object currently being assigned to

xtring paramname;
xtring strparam;
double numparam;
bool ifhelp=false;
bool includepft;

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
		declareitem("nyear",&nyear,1,10000,1,CB_NONE,"Number of simulation years");
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
		declareitem("respcoeff",&ppft->respcoeff,0.0,1.2,1,CB_NONE,
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
		declareitem("k_allom2",&ppft->k_allom2,10.0,1.0e4,1,CB_NONE,
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
		else if (strparam.upper()=="CROP") ppft->lifeform=CROP;
		else {
			sendmessage("Error",
				"Unknown lifeform type (valid types: \"TREE\", \"GRASS\")");
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
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("ifdailynpp")) badins("ifdailynpp");
		if (!itemparsed("ifdailydecomp")) badins("ifdailydecomp");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");
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
//   if (date.day==0 && date.year==nyear) return false;
//   // else
//   return true;
//
//   Currently the following member variables of the climate member of stand must be
//   initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
//   values of the climate variables derived from monthly means, this day's values
//   will presumably be extracted from arrays containing the interpolated daily
//   values (see function getstand):
//
//   stand.climate.temp=dtemp[date.day];
//   stand.climate.prec=dprec[date.day];
//   stand.climate.insol=dsun[date.day];
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
	// Will maintain a list of Coord objects containing coordinates
	// of the grid cells to simulate

int ngridcell; // the number of grid cells to simulate
bool firstgrid; // whether simulating first grid cell in linked list

// File names for temperature, precipitation, sunshine and soil code driver files
xtring file_temp,file_prec,file_sun,file_soil;

// Output file names ...
xtring file_cmass,file_anpp,file_lai,file_flux;

// ... and streams
FILE *out_cmass,*out_anpp,*out_lai,*out_flux;

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
				dprintf("Unknown option \"%s\"\n",insfilename);
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

	gridlist.killall();

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	ngridcell=0;
	while (!eof) {
		
		// Read next record in file
		eof=!readfor(in_grid,"f,f,a",&dlon,&dlat,&descrip);

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

	// Retrieve input file names as read from ins file

	file_temp=param["file_temp"].str;
	file_prec=param["file_prec"].str;
	file_sun=param["file_sun"].str;
	file_soil=param["file_soil"].str;
	
	// Retrieve output file names as read from ins file

	file_cmass=param["file_cmass"].str;
	file_anpp=param["file_anpp"].str;
	file_lai=param["file_lai"].str;
	file_flux=param["file_flux"].str;

	// If producing annual output, open files now
	// (This means that files will be open throughout the simulation - if this is a
	// problem, switch annual output off, or store annual results in an array or
	// dynamic list and write them out at the end of the simulation for each grid
	// cell - see function outannual below)

	if (annual_output) {
		out_cmass=fopen(file_cmass,"w");
		out_anpp=fopen(file_anpp,"w");
		out_lai=fopen(file_lai,"w");
		out_flux=fopen(file_flux,"w");

		if (!out_cmass || !out_anpp || !out_lai || !out_flux)
			fail("initio: could not open file(s) for output");
	}

	// Set timers

	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);

	// Start at first object in linked list of grid cell coordinates ...
	firstgrid=true;
}


///////////////////////////////////////////////////////////////////////////////////////
// GETSTAND
// Called by the framework at the start of the simulation for a particular stand

bool getstand(Stand& stand) {

	// DESCRIPTION
	// Obtains latitude and soil static parameters for the next stand (grid cell) to
	// simulate. The function should returns false if no stands remain to be simulated,
	// otherwise true. Currently the following member variables of stand should be
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
	
	if (firstgrid) {
		gridlist.firstobj();
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		// Reset random number seed for stochastic rainday generator
		setseed(1234567);

		// Retrieve coordinate of next grid cell from linked list
		Coord& c=gridlist.getobj();

		// Load environmental data for this grid cell from files
		// (these will be the same for every year of the simulation, but must be sent
		// anew to the framework each year in function getclimate, below)

		readenv(c);

		dprintf("\nCommencing simulation for stand at (%g,%g)",c.lon,c.lat);
		if (c.descrip!="") dprintf(" (%s)\n",(char*)c.descrip);
		else dprintf("\n");
		
		// Tell framework the latitude of this grid cell
		stand.climate.lat=c.lat;
		
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		stand.climate.instype=SUNSHINE;

		// Tell framework the soil type of this grid cell
		soilparameters(stand.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}


///////////////////////////////////////////////////////////////////////////////////////
// GETCLIMATE
// Called by the framework each simulation day before any process modelling is
// performed for this day

bool getclimate(Stand& stand) {

	// DESCRIPTION
	// Obtains climate data (including atmospheric CO2 and insolation) for this day.
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
	// stand.climate.temp=dtemp[date.day];
	// stand.climate.prec=dprec[date.day];
	// stand.climate.insol=dsun[date.day];

	double progress;

	// Send environmental values for today to framework

	stand.climate.co2=co2;
	stand.climate.temp=dtemp[date.day];
	stand.climate.prec=dprec[date.day];
	stand.climate.insol=dsun[date.day];

	// First day of year only ...

	if (date.day==0) {

		// Return false if last year was the last for the simulation
		if (date.year==nyear) return false;

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*nyear+date.year)/
				(double)(ngridcell*nyear);
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////
// OUTANNUAL
// Called by the framework at the end of the last day of each simulation year

void outannual(Stand& stand,Pftlist& pftlist) {

	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int p,c,nclass;
	double cmass_stand,anpp_stand,lai_stand;
	double flux_veg,flux_soil,flux_fire,flux_est;
	double lon,lat;

	if (vegmode==COHORT)
		nclass=min(date.year/estinterval+1,OUTPUT_MAXAGECLASS);

	if (date.year==0 && firstgrid) {

		// Very first time only

		// Open output files if annual output suppressed
		// (otherwise they were already opened in function initio)

		if (!annual_output) {
			out_cmass=fopen(file_cmass,"w");
			out_anpp=fopen(file_anpp,"w");
			out_lai=fopen(file_lai,"w");
			out_flux=fopen(file_flux,"w");

			if (!out_cmass || !out_anpp || !out_lai || !out_flux)
				fail("outannual: could not open file(s) for output");
		}

		// Print column labels

		fprintf(out_cmass,"%6s%6s%6s","Lon","Lat","Year");
		fprintf(out_anpp,"%6s%6s%6s","Lon","Lat","Year");
		fprintf(out_lai,"%6s%6s%6s","Lon","Lat","Year");
		fprintf(out_flux,"%6s%6s%6s%8s%8s%8s%8s%8s\n","Lon","Lat","Year","Veg","Soil",
			"Fire","Est","NEE");

		// Loop through PFT's and print PFT names as column labels

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pft=pftlist.getobj();
			fprintf(out_cmass,"%8s",(char*)pft.name);
			fprintf(out_anpp,"%8s",(char*)pft.name);
			fprintf(out_lai,"%8s",(char*)pft.name);
			pftlist.nextobj();
		}

		// Print labels for "Total" columns

		fprintf(out_cmass,"%8s\n","Total");
		fprintf(out_anpp,"%8s\n","Total");
		fprintf(out_lai,"%8s\n","Total");

		// Close files if annual output suppressed

		if (!annual_output) {
			fclose(out_cmass);
			fclose(out_anpp);
			fclose(out_lai);
			fclose(out_flux);
		}

		firstgrid=false;
	}

	if (annual_output || date.year==nyear-1) {

		// Each year (if annual output), or last year of simulation

		// Open output files in append mode if annual output suppressed

		if (!annual_output) {
			out_cmass=fopen(file_cmass,"a");
			out_anpp=fopen(file_anpp,"a");
			out_lai=fopen(file_lai,"a");
			out_flux=fopen(file_flux,"a");

			if (!out_cmass || !out_anpp || !out_lai || !out_flux)
				fail("outannual: could not open file(s) for output");			
		}

		lon=gridlist.getobj().lon;
		lat=gridlist.getobj().lat;

		cmass_stand=0.0;
		anpp_stand=0.0;
		lai_stand=0.0;

		// Print longitude, latitude, year

		fprintf(out_cmass,"%6.1f%6.1f%6d",lon,lat,date.year);
		fprintf(out_anpp,"%6.1f%6.1f%6d",lon,lat,date.year);
		fprintf(out_lai,"%6.1f%6.1f%6d",lon,lat,date.year);
		fprintf(out_flux,"%6.1f%6.1f%6d",lon,lat,date.year);

		// Loop through PFTs

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Standpft& standpft=stand.pft[pft.id];

			// Sum C biomass, NPP and LAI across patches and PFTs
			
			standpft.cmass_total=0.0;
			standpft.anpp_total=0.0;
			standpft.lai_total=0.0;

			// Initialise age structure array

			if (vegmode==COHORT || vegmode==INDIVIDUAL)
				for (c=0;c<nclass;c++)
					standpft.densindiv_ageclass[c]=0.0;
	
			for (p=0;p<npatch;p++) {
				Patch& patch=stand[p];
				Vegetation& vegetation=patch.vegetation;
				Patchpft& patchpft=patch.pft[pft.id];

				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv=vegetation.getobj();
					
					if (indiv.pft.id==pft.id) {
						standpft.cmass_total+=indiv.cmass_leaf+
							indiv.cmass_root+indiv.cmass_sap+indiv.cmass_heart-indiv.cmass_debt;
						standpft.anpp_total+=indiv.anpp;
						standpft.lai_total+=indiv.lai;

						if (vegmode==COHORT || vegmode==INDIVIDUAL) {
						
							// Age structure
							
							c=indiv.age/estinterval;
							if (c<OUTPUT_MAXAGECLASS)
								standpft.densindiv_ageclass[c]+=indiv.densindiv;
						}
					}

					vegetation.nextobj();
				}
			}

			standpft.cmass_total/=(double)npatch;
			standpft.anpp_total/=(double)npatch;
			standpft.lai_total/=(double)npatch;

			// Update stand totals

			cmass_stand+=standpft.cmass_total;
			anpp_stand+=standpft.anpp_total;
			lai_stand+=standpft.lai_total;

			// Print PFT sums to files

			fprintf(out_cmass,"%8.3f",standpft.cmass_total);
			fprintf(out_anpp,"%8.3f",standpft.anpp_total);
			fprintf(out_lai,"%8.4f",standpft.lai_total);

			// Graphical output every 10 years
			// (Windows shell only - "plot" statements have no effect otherwise)

			if (!(date.year%10)) {
				plot("cmass",pft.name,date.year,stand.pft[pft.id].cmass_total);
				plot("anpp",pft.name,date.year,stand.pft[pft.id].anpp_total);
				plot("lai",pft.name,date.year,stand.pft[pft.id].lai_total);
			}

			pftlist.nextobj();
		}

		flux_veg=flux_soil=flux_fire=flux_est=0.0;

		// Sum C fluxes across patches

		for (p=0;p<npatch;p++) {
			flux_veg+=stand[p].fluxes.acflux_veg/(double)npatch;
			flux_soil+=stand[p].fluxes.acflux_soil/(double)npatch;
			flux_fire+=stand[p].fluxes.acflux_fire/(double)npatch;
			flux_est+=stand[p].fluxes.acflux_est/(double)npatch;
		}

		// Print stand totals to files

		fprintf(out_cmass,"%8.3f\n",cmass_stand);
		fprintf(out_anpp,"%8.3f\n",anpp_stand);
		fprintf(out_lai,"%8.4f\n",lai_stand);
		
		// Graphical output every 10 years
		// (Windows shell only - no effect otherwise)

		if (!(date.year%10)) {
			plot("fluxes","flux_veg",date.year,flux_veg);
			plot("fluxes","flux_soil",date.year,flux_soil);
			plot("fluxes","flux_fire",date.year,flux_fire);
			plot("fluxes","flux_est",date.year,flux_est);
			plot("fluxes","NEE",date.year,flux_veg+flux_soil+flux_fire+flux_est);
			plot("soilc","slow",date.year,stand[0].soil.cpool_slow);
			plot("soilc","fast",date.year,stand[0].soil.cpool_fast);
		}

		// Write fluxes to file

		fprintf(out_flux,"%8.3f%8.3f%8.4f%8.4f%8.3f\n",flux_veg,flux_soil,flux_fire,
			flux_est,flux_veg+flux_soil+flux_fire+flux_est);

		// Output of age structure (Windows shell only - no effect otherwise)

		if (vegmode==COHORT || vegmode==INDIVIDUAL) {

			if (!(date.year%20) && date.year<2000) {
			
				resetwindow("age_structure");

				pftlist.firstobj();
				while (pftlist.isobj) {
					Pft& pft=pftlist.getobj();

					if (pft.lifeform==TREE) {

						Standpft& standpft=stand.pft[pft.id];

						for (c=0;c<nclass;c++)
							plot("age_structure",pft.name,
								c*estinterval+estinterval/2,
								standpft.densindiv_ageclass[c]/(double)npatch);
					}
					
					pftlist.nextobj();
				}
			}
		}

		// Close files if annual output suppressed

		if (!annual_output) {
			fclose(out_cmass);
			fclose(out_anpp);
			fclose(out_lai);
			fclose(out_flux);
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

		fclose(out_cmass);
		fclose(out_anpp);
		fclose(out_lai);
		fclose(out_flux);
	}

	// Clean up

	gridlist.killall();

	dprintf("\nModel run terminated\n");
}

#endif // USE_DEMO_IO
