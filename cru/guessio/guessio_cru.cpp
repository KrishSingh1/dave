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
#include "GlobalNdep.h"

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
// guess2008 - changed from nyear to nyear_spinup
int nyear_spinup; // number of simulation years during spinup
// guess2008 - new optional parameter
int searchradius; // search radius to use when finding CRU data

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
xtring file_cton,file_nmass,file_nsources,file_npool,file_nleach,file_age,file_nuptake,file_anppn,file_vmaxnlim;
// end GUESSN

// GUESSN allometry
xtring file_allometry;
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

	// guess2008 - initialise filenames here
	outputdirectory = "";
	file_cmass=file_anpp=file_lai=file_cflux=file_dens=file_runoff="";
	file_mnpp=file_mlai=file_maet=file_mpet=file_mevap=file_mrunoff=file_mintercep=file_mrh="";
	file_mgpp=file_mra=file_mnee=file_mwcont_upper=file_mwcont_lower="";
	file_cpool=file_firert="";

	file_dgpp="";

	// GUESSN
	file_cton=file_nmass=file_nsources=file_npool=file_nleach=file_age=file_nuptake=file_anppn=file_vmaxnlim="";
	// end GUESSN

	// GUESSN allometry
	file_allometry="";
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

		// FACE David
		declareitem("FACE_ring",&FACE_ring,1,2,1,CB_NONE,"Which FACE ring to simulated");
		declareitem("FACE_DUKE_or_OAK",&ifduke,0,1,1,CB_NONE,"Which FACE site to simulated");
		declareitem("Has_FACE_clim",&has_FACE_clim,0,1,1,CB_NONE,"Using FACE clim data");
		
		// CANIF David
		declareitem("Has_CANIF_clim",&has_CANIF_clim,0,1,1,CB_NONE,"Using CANIF clim data");

		declareitem("title",&title,80,CB_NONE,"Title for run");
		// guess2008 - changed this input parameter name from nyear to nyear_spinup, which is more descriptive 
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

		declareitem("iflimvmax",&iflimvmax,1,CB_NONE,
			"Whether Vmax limited by actual leaf N");
		declareitem("nrelocfrac",&nrelocfrac,0.0,1.0,1,CB_NONE,
			"Fractional N relocation from shed leaves & roots");
		declareitem("ifvarycn",&ifvarycn,1,CB_NONE,
			"Whether leaf and tissue C:N adjusted depending on Vmax");
		declareitem("ifnlimvarycn",&ifnlimvarycn,1,CB_NONE,
			"Whether leaf and tissue C:N adjusted depending on N limitation");
		declareitem("ifnfix",&ifnfix,0,3,1,CB_NONE,
			"Whether to include an estimate for N fixation");
		declareitem("andep",&andep,0.0,100.0,1,CB_NONE,
			"Annual N deposition kgN/m2 (if not read from file)");
		declareitem("minndep",&minndep,0.0,0.001,1,CB_NONE,"Minimum N dep");

		declareitem("ifcentury",&ifcentury,1,CB_NONE,
			"Whether to use CENTURY SOM dynamics (default standard LPJ)");
		declareitem("ifnlim",&ifnlim,1,CB_NONE,
			"Whether plant growth limited by available N");
		declareitem("freenyears",&freenyears,0,1000,1,CB_NONE,
			"Number of years to spinup without N limitation");
		declareitem("ifdailysetntoc",&ifdailysetntoc,1,CB_NONE,
			"If to use daily version of setntoc (set N:C ratio of som pools)");
		declareitem("ifleachn",&ifleachn,1,CB_NONE,
			"Whether to allow N leaching");
		declareitem("ifindiv_fnuptake",&ifindiv_fnuptake,1,CB_NONE,
			"Whether to allow individual fractional N uptake");
		declareitem("cwdtransfer",&cwdtransfer,0.0,1.0,1,CB_NONE,
			"Fraction of woody debris transferred to SOM each year");
		declareitem("nmass_avail_max",&nmass_avail_max,0.0,0.1,1,CB_NONE,
			"Max N:C ratio in the soil (should be 0.002 (Parton et al 1993, Fig. 4))");
		declareitem("ifnstorage",&ifnstorage,1,CB_NONE,
			"If to use a N storage for each individual");
		declareitem("ifndemand_new_est",&ifndemand_new_est,1,CB_NONE,
			"N limitation of new establishment");
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
		declareitem("file_age",&file_age,300,CB_NONE,"Age structure output file");
		declareitem("file_nuptake",&file_nuptake,300,CB_NONE,"annual N uptake output file");
		declareitem("file_anppn",&file_anppn,300,CB_NONE,"annual N usage output file");
		declareitem("file_vmaxnlim",&file_vmaxnlim,300,CB_NONE,"annual N limitation on vm output file");
		// end GUESSN

		// GUESSN allometry
		declareitem("file_allometry",&file_allometry,300,CB_NONE,"Allometry output file");
		// GUESSN
		
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
		if (!itemparsed("nyear_spinup")) badins("nyear_spinup"); // guess2008
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("ifdailynpp")) badins("ifdailynpp");
		if (!itemparsed("ifdailydecomp")) badins("ifdailydecomp");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");

		// FACE David
		if (!itemparsed("FACE_ring")) badins("FACE_ring");
		if (!itemparsed("FACE_DUKE_or_OAK")) badins("FACE_DUKE_or_OAK");
		if (!itemparsed("Has_FACE_clim")) badins("Has_FACE_clim");

		// CANIF David
		if (!itemparsed("Has_CANIF_clim")) badins("Has_CANIF_clim");

		// GUESSN

		if (!itemparsed("iflimvmax")) badins("iflimvmax");
		if (!itemparsed("nrelocfrac")) badins("nrelocfrac");
		if (!itemparsed("ifvarycn")) badins("ifvarycn");
		if (!itemparsed("ifnlimvarycn")) badins("ifnlimvarycn");
		if (!itemparsed("ifnfix")) badins("ifnfix");
		if (!itemparsed("andep")) badins("andep");
		if (!itemparsed("minndep")) badins("minndep");

		if (!itemparsed("ifcentury")) badins("ifcentury");
		if (!itemparsed("ifnlim")) badins("ifnlim");
		if (!itemparsed("freenyears")) badins("freenyears");
		if (!itemparsed("ifdailysetntoc")) badins("ifdailysetntoc");
		if (!itemparsed("ifleachn")) badins("ifleachn");
		if (!itemparsed("ifindiv_fnuptake")) badins("ifindiv_fnuptake");
		if (!itemparsed("cwdtransfer")) badins("cwdtransfer");
		if (!itemparsed("nmass_avail_max")) badins("nmass_avail_max");
		if (!itemparsed("ifnstorage")) badins("ifnstorage");
		if (!itemparsed("ifndemand_new_est")) badins("ifndemand_new_est");

		// end GUESSN


		// guess2008
		if (!itemparsed("outputdirectory")) badins("outputdirectory");
		if (!itemparsed("ifsmoothgreffmort")) badins("ifsmoothgreffmort");
		if (!itemparsed("ifdroughtlimitedestab")) badins("ifdroughtlimitedestab");
		if (!itemparsed("ifrainonwetdaysonly")) badins("ifrainonwetdaysonly");
		if (!itemparsed("ifspeciesspecificwateruptake")) badins("ifspeciesspecificwateruptake");



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

// guess2008
const int NYEAR_HIST=106; // guess2008 - CRU TS 3.0 has 106 years of data (1901-2006)
	// number of years of historical climate in CRU and CO2 files (see below)
const int FIRSTHISTYEAR=1901;
	// calender year corresponding to first year in CRU climate data set
const int NYEAR_SPINUP_DATA=30;
	// number of years to use for temperature-detrended spinup data set
	// (not to be confused with the number of years to spinup model for, which
	// is read from the ins file)

// FACE DAVID climate dataset NYEAR_SCENARIO_FACE
//const int NYEAR_SCENARIO_FACE=11;	// ORNL
const int NYEAR_SCENARIO_FACE=12;	// Duke

// CANIF DAVID max years with climate data
const int MAXNYEAR_SCENARIO_CANIF=8;
const int NSITES=8;

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

FILE *out_cton,*out_nmass, *out_nsources, *out_npool, *out_nleach, *out_age, *out_nuptake, *out_anppn, *out_vmaxnlim;
// end GUESSN

// GUESSN allometry
FILE *out_allometry;
// end GUESSN

// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

// CO2 data for each year of historical data set
double co2[NYEAR_HIST];

// FACE DAVID climate CO2 data for each day of the FACE scenario
double dco2_FACE[NYEAR_SCENARIO_FACE][365];

// Monthly temperature, precipitation and sunshine data for current grid cell
// and historical period
double hist_mtemp[NYEAR_HIST][12];
double hist_mprec[NYEAR_HIST][12];
double hist_msun[NYEAR_HIST][12];

// FACE DAVID climate Daily temperature, precipitation and sunshine data for FACE scenario 
double dtemp_FACE[NYEAR_SCENARIO_FACE][365];
double dprecip_FACE[NYEAR_SCENARIO_FACE][365];
double dsun_FACE[NYEAR_SCENARIO_FACE][365];

// FACE DAVID climate Yearly ndep data for FACE scenario 
double yndep_FACE[NYEAR_NDEP];	// inputdata ranges from 1750-2007

// CANIF DAVID climate Daily temperature, precipitation, sunshine and N dep data for CANIF scenarios
double dtemp_CANIF[NSITES][MAXNYEAR_SCENARIO_CANIF][365];
double dprecip_CANIF[NSITES][MAXNYEAR_SCENARIO_CANIF][365];
double dsun_CANIF[NSITES][MAXNYEAR_SCENARIO_CANIF][365];
double lonlatyearsndep[NSITES][6];

// guess2008
// Monthly frost days, precipitation days and DTR data for current grid cell
// and historical period
double hist_mfrs[NYEAR_HIST][12];
double hist_mwet[NYEAR_HIST][12];
double hist_mdtr[NYEAR_HIST][12];


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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// FACE DAVID climate Reads in climate data for sceanrio period from temp, sun, precip and CO2 scenario files
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_FACE_clim(double FACE_dtemp[NYEAR_SCENARIO_FACE][365],double FACE_dprec[NYEAR_SCENARIO_FACE][365],
	double FACE_dsun[NYEAR_SCENARIO_FACE][365],double FACE_dco2[NYEAR_SCENARIO_FACE][365], 
	double FACE_yndep[NYEAR_NDEP],int NYEAR_SCENARIO_FACE,int NYEAR_NDEP)
{
	xtring filename_FACE;
	FILE* file_FACE_met;
	FILE* file_FACE_ndep;

	if (ifduke)
	{
		filename_FACE=param["met_face_duke"].str;
		file_FACE_met=fopen(filename_FACE,"rt");

		if (!file_FACE_met) 
			fail("FACE: could not open file %s for input",(char*)filename_FACE);

		filename_FACE=param["ndep_face_duke"].str;
		file_FACE_ndep=fopen(filename_FACE,"rt");

		if (!file_FACE_ndep) 
			fail("FACE: could not open file %s for input",(char*)filename_FACE);

		double met_duke[15];
		double ndep_duke[3];

		for (int year = 0; year < NYEAR_SCENARIO_FACE; year++)  // go through the years of data
		{
			for (int day = 0; day < 365; day++)
			{
				readfor(file_FACE_met, "15f", met_duke);				// read data

				FACE_dsun[year][day] = met_duke[1] / 4.56 * 1000000;	// val_duke is in umol/m2/day, dsun should be in J/m2/day 1.8umol/m2/day == 1J/m2/day
																		// http://www.hydro.co.nz/1_information/1_light_info/info_light.html
																		// Thomas got 4.56 from meeting in Santa Barbara!!! Better value
				FACE_dtemp[year][day] = met_duke[2];
				FACE_dprec[year][day] = met_duke[5];

				if (FACE_ring == 2)		// David co2 elevated concentration 
					FACE_dco2[year][day] = met_duke[8];
				else					// David co2 ambient concentration
					FACE_dco2[year][day] = met_duke[7];
			}
		}

		for (year = 0; year < NYEAR_NDEP; year++)
		{
			readfor(file_FACE_ndep, "3f", ndep_duke);

			FACE_yndep[year] = ndep_duke[1];	// Two different columns in the input file with ndep data (1,2)
		}
	}
	else // ORNL
	{
		filename_FACE=param["met_face_oak"].str;
		file_FACE_met=fopen(filename_FACE,"rt");

		if (!file_FACE_met) 
			fail("FACE: could not open file %s for input",(char*)filename_FACE);

		filename_FACE=param["ndep_face_oak"].str;
		file_FACE_ndep=fopen(filename_FACE,"rt");

		if (!file_FACE_ndep) 
			fail("FACE: could not open file %s for input",(char*)filename_FACE);

		double met_oak[11];
		double ndep_oak[3];

		for (int year = 0; year < NYEAR_SCENARIO_FACE; year++)  // go through the years of data
		{
			for (int day = 0; day < 365; day++)
			{
				readfor(file_FACE_met, "11f", met_oak);				// read data

				FACE_dsun[year][day] = met_oak[2] / 4.56 * 1000000;	// val_duke is in mol/m2/day, dsun should be in J/m2/day 1.8umol/m2/day == 1J/m2/day
																	// http://www.hydro.co.nz/1_information/1_light_info/info_light.html
																	// Thomas got 4.56 from meeting in Santa Barbara!!! Better value
				FACE_dtemp[year][day] = met_oak[3];
				FACE_dprec[year][day] = met_oak[10];
				if (FACE_ring == 2)	// David co2 elevated concentration 
					FACE_dco2[year][day] = (met_oak[6] + met_oak[7])/2.0;	// average of ring 1 and 2
				else				// David co2 ambient concentration 
					FACE_dco2[year][day] = (met_oak[8] + met_oak[9])/2.0;	// average of ring 4 and 5
			}
		}

		for (year = 0; year < NYEAR_NDEP; year++)
		{
			readfor(file_FACE_ndep, "3f", ndep_oak);

			FACE_yndep[year] = ndep_oak[1];	// Two different columns in the input file with ndep data (1,2)
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// FACE DAVID climate Reads in climate data for sceanrio period from temp, sun, precip and CO2 scenario files
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_CANIF_clim(double CANIF_dtemp[NSITES][MAXNYEAR_SCENARIO_CANIF][365],
	double CANIF_dprec[NSITES][MAXNYEAR_SCENARIO_CANIF][365],
	double CANIF_dsun[NSITES][MAXNYEAR_SCENARIO_CANIF][365],
	double lonlatyearsndep[NSITES][6],int NSITES,
	double lon,double lat)
{
	xtring filename_CANIF;
	FILE* file_CANIF_met;
	
	filename_CANIF=param["met_canif"].str;
	file_CANIF_met=fopen(filename_CANIF,"rt");

	if (!file_CANIF_met) 
		fail("CANIF: could not open file %s for input",(char*)filename_CANIF);

	double met_canif[5];

	for (int site_lonlat = 0; site_lonlat < NSITES; site_lonlat++) {
		readfor(file_CANIF_met, "6f", lonlatyearsndep[site_lonlat]);
		if (lonlatyearsndep[site_lonlat][0]==lon && lonlatyearsndep[site_lonlat][1]==lat)
			WSITE = site_lonlat;
	}

	for (int site = 0; site < NSITES; site++)  // go through the site data
	{
		for (int year = 0; year < lonlatyearsndep[site][2];year++)
		{
			for (int day = 0; day < 365; day++)
			{
				readfor(file_CANIF_met, "5f", met_canif);				// read data

				CANIF_dsun[site][year][day] = met_canif[2];	
				CANIF_dtemp[site][year][day] = met_canif[3];
				CANIF_dprec[site][year][day] = met_canif[4];
			}
		}
	}
	
}
//----------------------------------------------------------------------------------------------


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
	double mtemp[NYEAR_HIST][12],double mprec[NYEAR_HIST][12],
	double msun[NYEAR_HIST][12]) {

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
		soilcode=(int)data.soilcode[0];


		for (y=0;y<NYEAR_HIST;y++) {
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
	double mfrs[NYEAR_HIST][12],double mwet[NYEAR_HIST][12],
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
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000

		for (y=0;y<NYEAR_HIST;y++) { 
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
                        int& scode, double hist_mtemp1[NYEAR_HIST][12], 
                        double hist_mprec1[NYEAR_HIST][12], 
                        double hist_msun1[NYEAR_HIST][12]) {

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
	if (!in) fail("readco2: could not open CO2 file %s for input",
		(char*)filename);

	for (year=0;year<NYEAR_HIST;year++) {
		readfor(in,"i,f",&calender_year,&co2[year]);
		if (calender_year!=FIRSTHISTYEAR+year)
			fail("readco2: %s, line %d - incorrect year specified",
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

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;

	// GUESSN
	file_ndep=param["file_ndep"].str;
	if (file_ndep=="")
		ifndepdata=false;
	else {
		FILE* in_ndep=fopen(file_ndep,"rt");
		if (!in_ndep)
			fail("initio: could not open %s for input",(char*)file_ndep);

		fclose(in_ndep);
		ifndepdata=true;
	}
	// end GUESSN
	
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

	// Read CO2 data from file
	readco2();


	// Remember whether to produce output each year or not
	annual_output=param["annual_output"].num;


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

	if (file_age!="") {
		file_age = outputdirectory + file_age;
		out_age=fopen(file_age,"w");
		if (!out_age) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_age);
	}
	else out_age=NULL;

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
	// end GUESSN

	// GUESSN allometry
	if (file_allometry!="") {
		file_allometry = outputdirectory + file_allometry;
		out_allometry=fopen(file_allometry,"w");
		if (!out_allometry) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_allometry);
	}
	else out_allometry=NULL;
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


// GUESSN
/// Retrieves nitrogen deposition for a particular grid cell
/** The values are either taken from the andep parameter in the instruction
 *  file, or from a binary archive file.
 *
 *  The binary archive has nitrogen deposition in gN/m2/year for the years
 *  1860, 1993 and 2050 (Galloway et. al., 2004).
 *
 *  Returned values will not be smaller than minndep (ins file parameter).
 *
 *  \param  filename    The file name of the binary archive
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 *  \param  xandep1860  Nitrogen deposition for 1860 (kgN/m2/year)
 *  \param  xandep1993  Nitrogen deposition for 1993 (kgN/m2/year)
 *  \param  xandep2050  Nitrogen deposition for 2050 (kgN/m2/year)
 */
bool getndep(xtring filename,double lon,double lat,double &xandep1860,double &xandep1993,double &xandep2050) {

	if (!ifndepdata) {
		xandep1860=andep;
		xandep1993=andep;
		xandep2050=andep;
		return true;
	}

	GlobalNdepArchive ark;
	if (!ark.open(filename)) {
		 fail("Could not open %s for input",(char*)filename);
		 return false;
	}

	GlobalNdep rec;
	rec.longitude = lon;
	rec.latitude = lat;
	
	if (!ark.getindex(rec)) {
		 // The coordinate wasn't found in the archive
		 ark.close();
		 return false;
	}
	else {
		 // Found the record, get the values
		 // Convert from gN to kgN, don't allow values smaller than minndep
		 xandep1860 = max(minndep, rec.ndep1860[0]*0.001);
		 xandep1993 = max(minndep, rec.ndep1993[0]*0.001);
		 xandep2050 = max(minndep, rec.ndep2050[0]*0.001);

		 ark.close();
		 return true;
	}
}
// end GUESSN


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
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	bool gridfound;

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
		gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
		                               hist_mtemp, hist_mprec, hist_msun);

		if (gridfound) // Get more historical CRU data for this grid cell
			gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);


		// FACE DAVID climate Reading met data
		if (has_FACE_clim)		
			read_FACE_clim(dtemp_FACE,dprecip_FACE,dsun_FACE,dco2_FACE,yndep_FACE,NYEAR_SCENARIO_FACE,NYEAR_NDEP);

		// CANIF DAVID climate Reading met data
		if (has_CANIF_clim)	{
			read_CANIF_clim(dtemp_CANIF,dprecip_CANIF,dsun_CANIF,lonlatyearsndep,NSITES,lon,lat);
			stand.plantyear=(nyear_spinup+NYEAR_HIST-lonlatyearsndep[WSITE][4]);
			dprintf("Plant year %d \n",stand.plantyear);
		}

		while (!gridfound) {

			dprintf("\nError: could not find stand at (%g,%g) in CRU data file\n",
				gridlist.getobj().lon,gridlist.getobj().lat);

			gridlist.nextobj();
			if (gridlist.isobj) {
				double lon = gridlist.getobj().lon;
				double lat = gridlist.getobj().lat;
				gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
				                               hist_mtemp, hist_mprec, hist_msun);
			  
				if (gridfound) // Get more historical CRU data for this grid cell
					gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation,
					                           hist_mfrs, hist_mwet, hist_mdtr);

			}
			else return false;
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

		// GUESSN
		if (!getndep(file_ndep,lon,lat,stand.climate.andep_1860,
			stand.climate.andep_1993,stand.climate.andep_2050)) {

			fail("Grid cell not found in %s",(char*)file_ndep);
		}
		else {
		//	dprintf("N deposition for longitude (%g) and latitude (%g) is\n1860: (%g) 1993: (%g) and 2050: (%g) kgN/ha/yr\n",
		//		gridlist.getobj().lon,gridlist.getobj().lat,stand.climate.andep_1860*10000.0,
		//			stand.climate.andep_1993*10000.0,stand.climate.andep_2050*10000.0);
		}
		// end GUESSN

		// CANIF David
		if (has_CANIF_clim)
			stand.climate.andep_1993=lonlatyearsndep[WSITE][5];

		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the latitude of this grid cell
		stand.climate.lat=gridlist.getobj().lat;
		
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		stand.climate.instype=SUNSHINE;

		if (has_FACE_clim) 
			soilcode = 4; // FACE climate soilcode

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
// THISYEARSNDEP
// Called by getclimate(). Calculates this years N deposition from three known values.
// needs to be redone

double thisyearsndep(double ndep_1860, double ndep_1993, double ndep_2050, int year, int nyears_spin, int firsthistyear) {

	if (ndep_1860 == ndep_1993)
		return ndep_1860;

	int hist_year = firsthistyear - (nyears_spin - year);

	// CANIF David
	if(has_CANIF_clim) {
		double years_shift=(FIRSTHISTYEAR+NYEAR_HIST)-(lonlatyearsndep[WSITE][3]+lonlatyearsndep[WSITE][2]);
		hist_year -= (years_shift + lonlatyearsndep[WSITE][3] - 1993);
	}
		

	if (hist_year <= 1860)
		return ndep_1860;
	else if (hist_year <= 1993)
		return ndep_1860+(hist_year-1860)*((ndep_1993-ndep_1860)/133.0);
	else
		return ndep_1993+(hist_year-1993)*((ndep_2050-ndep_1993)/57.0);
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

	// guess2008 - changed name from mwet to mwet_all
	double mwet_all[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // number of rain days per month
	int dd;

	if (date.day==0) {

		// First day of year ...
		
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
			
			/*if (!(date.year%100) && date.day==0) {

				dprintf("\nClimate for year %d\n",date.year+1);
				dprintf("\n        Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec\n");
				dprintf("Temp %6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f\n",
					mtemp[0],mtemp[1],mtemp[2],mtemp[3],mtemp[4],mtemp[5],
					mtemp[6],mtemp[7],mtemp[8],mtemp[9],mtemp[10],mtemp[11]);
				dprintf("Prec %6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n",
					mprec[0],mprec[1],mprec[2],mprec[3],mprec[4],mprec[5],
					mprec[6],mprec[7],mprec[8],mprec[9],mprec[10],mprec[11]);
				dprintf("Sun  %6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n",
					msun[0],msun[1],msun[2],msun[3],msun[4],msun[5],
					msun[6],msun[7],msun[8],msun[9],msun[10],msun[11]);
			}*/
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

			/*if (!(date.year%20) && date.day==0) {

				dprintf("\nClimate for year %d\n",date.year+1);
				dprintf("\n        Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec\n");
				dprintf("Temp %6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f\n",
					hist_mtemp[date.year-nyear_spinup][0],hist_mtemp[date.year-nyear_spinup][1],
					hist_mtemp[date.year-nyear_spinup][2],hist_mtemp[date.year-nyear_spinup][3],
					hist_mtemp[date.year-nyear_spinup][4],hist_mtemp[date.year-nyear_spinup][5],
					hist_mtemp[date.year-nyear_spinup][6],hist_mtemp[date.year-nyear_spinup][7],
					hist_mtemp[date.year-nyear_spinup][8],hist_mtemp[date.year-nyear_spinup][9],
					hist_mtemp[date.year-nyear_spinup][10],hist_mtemp[date.year-nyear_spinup][11]);
				dprintf("Prec %6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n",
					hist_mprec[date.year-nyear_spinup][0],hist_mprec[date.year-nyear_spinup][1],
					hist_mprec[date.year-nyear_spinup][2],hist_mprec[date.year-nyear_spinup][3],
					hist_mprec[date.year-nyear_spinup][4],hist_mprec[date.year-nyear_spinup][5],
					hist_mprec[date.year-nyear_spinup][6],hist_mprec[date.year-nyear_spinup][7],
					hist_mprec[date.year-nyear_spinup][8],hist_mprec[date.year-nyear_spinup][9],
					hist_mprec[date.year-nyear_spinup][10],hist_mprec[date.year-nyear_spinup][11]);
				dprintf("Sun  %6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n",
					hist_msun[date.year-nyear_spinup][0],hist_msun[date.year-nyear_spinup][1],
					hist_msun[date.year-nyear_spinup][2],hist_msun[date.year-nyear_spinup][3],
					hist_msun[date.year-nyear_spinup][4],hist_msun[date.year-nyear_spinup][5],
					hist_msun[date.year-nyear_spinup][6],hist_msun[date.year-nyear_spinup][7],
					hist_msun[date.year-nyear_spinup][8],hist_msun[date.year-nyear_spinup][9],
					hist_msun[date.year-nyear_spinup][10],hist_msun[date.year-nyear_spinup][11]);
			}*/

		}
	}

	// Send environmental values for today to framework

	if (date.year<nyear_spinup)
		stand.climate.co2=co2[0];
	else if (date.year<nyear_spinup+NYEAR_HIST)
		stand.climate.co2=co2[date.year-nyear_spinup];

//	if (date.year > nyear_spinup+NYEAR_HIST-10)
//		stand.climate.co2=550;

	// GUESSN calculate annual N deposition value
	stand.climate.andep=thisyearsndep(stand.climate.andep_1860,stand.climate.andep_1993,
		stand.climate.andep_2050,date.year,nyear_spinup,FIRSTHISTYEAR);
	// end GUESSN

	//if (!(date.year%20) && date.day==0 && date.year>=nyear_spinup ||!(date.year%100) && date.day==0 && date.year<nyear_spinup) 
	//	dprintf("CO2  %6.0f         Ndep %7.3f (kgN/ha/yr)\n",stand.climate.co2,stand.climate.andep*10000.0);

		///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	//
	// FACE DAVID climate - Using Duke or ORNL data set instead. This code 
	//				just writes over old climate data
	//
	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////



	if (has_FACE_clim && (SOILDEPTH_LOWER != 100.0 || SOILDEPTH_UPPER != 400.0 || NYEAR_SCENARIO_FACE != 12) && ifduke)
		fail("DUKE input data is wrong!!!\n");

	if (has_FACE_clim && (SOILDEPTH_LOWER != 1500.0 || SOILDEPTH_UPPER != 500.0 || NYEAR_SCENARIO_FACE != 11) && !ifduke)
		fail("ORNL input data is wrong!!!\n");

	if (!has_FACE_clim && (SOILDEPTH_LOWER != 1000.0 || SOILDEPTH_UPPER != 500.0))
		fail("wrong soil depths!!!\n");

	FYEAR_SCENARIO_FACE = nyear_spinup+NYEAR_HIST-NYEAR_SCENARIO_FACE; 
	// Calculates which year we should start with the data so we end the simulation with the right order

	if (has_FACE_clim && date.year >= FYEAR_SCENARIO_FACE) //Last years the data is in the right order. 
	{
		stand.climate.instype=SUNSHINE;
		//stand.climate.instype=NETSWRAD;

		// Which year in the data set this year represent
		int scenario_year = date.year-FYEAR_SCENARIO_FACE;
	
		stand.climate.temp=dtemp_FACE[scenario_year][date.day]; // Temperature
		stand.climate.prec=dprecip_FACE[scenario_year][date.day];// Precipitation
		stand.climate.par=dsun_FACE[scenario_year][date.day];	// Photosynthetically-active radiation
		stand.climate.co2=dco2_FACE[scenario_year][date.day];	// CO2


		// Nitrogen deposition
		if (date.day == 0)
			stand.climate.andep=yndep_FACE[scenario_year+NYEAR_NDEP-NYEAR_SCENARIO_FACE] / 10000.0;  // from ha to m2

		//stand.climate.temp=dtemp[date.day];	// used when comparing clim data to cru data
		//stand.climate.prec=dprec[date.day];	// used when comparing clim data to cru data
		stand.climate.insol=dsun[date.day];

	}	
	else if (has_FACE_clim) // Just uses the clim data over and over again.
	{
		stand.climate.instype=SUNSHINE;
		//stand.climate.instype=NETSWRAD;

		// Which year in the data set this year represent
		int scenario_year;

		if (!ifduke){

			scenario_year = (date.year+6)%(NYEAR_SCENARIO_FACE);

			if (scenario_year == 7 || scenario_year == 0){
				scenario_year=2;
			}
		}

		if (ifduke){

			scenario_year = (date.year+6)%(NYEAR_SCENARIO_FACE);

			if (scenario_year == 4 || scenario_year == 3 || scenario_year == 9 || scenario_year == 7){
				scenario_year=2;
			}
		}

		stand.climate.temp=dtemp_FACE[scenario_year][date.day];
		stand.climate.prec=dprecip_FACE[scenario_year][date.day];	
		stand.climate.par=dsun_FACE[scenario_year][date.day];

		if(date.year < FYEAR_SCENARIO_FACE - 97)			
			stand.climate.co2=co2[0];
		else
			stand.climate.co2=co2[97+date.year-FYEAR_SCENARIO_FACE];

		// Nitrogen deposition
		if (date.day == 0)
		{
			if (date.year < nyear_spinup+NYEAR_HIST-NYEAR_NDEP)
				stand.climate.andep=yndep_FACE[0] / 10000;	// from ha to m2
			
			else 
				stand.climate.andep=yndep_FACE[date.year-(nyear_spinup+NYEAR_HIST-NYEAR_NDEP)] / 10000.0;	// from ha to m2
		}

		//stand.climate.temp=dtemp[date.day];	// used when comparing clim data to cru data
		//stand.climate.prec=dprec[date.day];	// used when comparing clim data to cru data
		stand.climate.insol=dsun[date.day];
		
	}

	///////////////////////////////////////////////////////////////////////////////
	/////////////////////////// End FACE climate //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////

	else {
		stand.climate.temp=dtemp[date.day];
		stand.climate.prec=dprec[date.day];
		stand.climate.insol=dsun[date.day];
	}

	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////// CANIF climate ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////


	FYEAR_SCENARIO_FACE = nyear_spinup+NYEAR_HIST-lonlatyearsndep[WSITE][2];

	double years_shift=(FIRSTHISTYEAR+NYEAR_HIST)-(lonlatyearsndep[WSITE][3]+lonlatyearsndep[WSITE][2]);

	if (has_CANIF_clim && date.year >= FYEAR_SCENARIO_FACE) //Last years the data is in the right order. 
	{
		stand.climate.instype=SUNSHINE;
		//stand.climate.instype=NETSWRAD;

		// Which year in the data set this year represent
		int scenario_year = date.year-FYEAR_SCENARIO_FACE;
	
		stand.climate.temp=dtemp_CANIF[WSITE][scenario_year][date.day]; // Temperature
		stand.climate.prec=dprecip_CANIF[WSITE][scenario_year][date.day];// Precipitation
		stand.climate.par=dsun_CANIF[WSITE][scenario_year][date.day];	// Photosynthetically-active radiation

		if (date.year<nyear_spinup)
			stand.climate.co2=co2[0];
		else if (date.year<nyear_spinup+NYEAR_HIST)
			stand.climate.co2=co2[date.year-nyear_spinup];

		stand.climate.insol=dsun[date.day];

	}	
	else if (has_CANIF_clim) // Just uses the clim data over and over again.
	{
		stand.climate.instype=SUNSHINE;
		//stand.climate.instype=NETSWRAD;

		int NSYR = (int)lonlatyearsndep[WSITE][2];

		// Which year in the data set this year represent
		int scenario_year = (date.year)%(NSYR);

		stand.climate.temp=dtemp_CANIF[WSITE][scenario_year][date.day]; // Temperature
		stand.climate.prec=dprecip_CANIF[WSITE][scenario_year][date.day];// Precipitation
		stand.climate.par=dsun_CANIF[WSITE][scenario_year][date.day];	// Photosynthetically-active radiation

		if (date.year<(nyear_spinup+years_shift))
			stand.climate.co2=co2[0];
		else if (date.year<(nyear_spinup+NYEAR_HIST+years_shift))
			stand.climate.co2=co2[date.year-nyear_spinup];

		stand.climate.insol=dsun[date.day];
		
	}

	// First day of year only ...

	if (date.day==0) {

		// Return false if last year was the last for the simulation
		if (date.year==nyear_spinup+NYEAR_HIST) return false;

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
// OUTANNUAL
// Called by the framework at the end of the last day of each simulation year

// guess2008 - many changes to this output routine. 
void outannual(Stand& stand,Pftlist& pftlist) {

	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int p,c,m,nclass;
	double cmass_stand,anpp_stand,lai_stand,runoff_stand,dens_stand;
	double flux_veg,flux_soil,flux_fire,flux_est;
	double c_litter,c_fast,c_slow; 
	double firert_stand; 

	// GUESSN
	double nmass_stand,n_litter,densindiv_ageclass_stand,nleach_stand,nuptake_stand,anppn_stand;
	double andep_stand,anmin_stand,animm_stand,anfix_stand,nsupply_stand,ndemand_stand,vmaxnlim_stand,total_dens;
	double surfsoillitterc,surfsoillittern,cwdc,cwdn,microc,micron,humusc,humusn,centuryc,centuryn;

	double stand_ageclass[OUTPUT_MAXAGECLASS];
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

	double dgpp[365];


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
		if (out_cflux) fprintf(out_cflux,lonlatyearstr_extended,"Lon","Lat","Year","Veg","Soil",
			"Fire","Est","NEE");
		if (out_cpool) 
			if (!ifcentury)
				fprintf(out_cpool,lonlatyearstr_extended,"Lon","Lat","Year","VegC","LittC",
					"SoilfC","SoilsC","Total");
			else // GUESSN
				fprintf(out_cpool,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%10s\n","Lon","Lat","Year","VegC","LittVC",
			"LittSC","CwdC","MicroC","HumusC","SoilC","Total");
			// end GUESSN

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

		if (out_dgpp) fprintf(out_dgpp,"%8s%8s%8s%8s%9s\n","Lon","Lat","Year","Day","GPP");
		
		// GUESSN
		if (out_cton) fprintf(out_cton,lonlatyearstr,"Lon","Lat","Year");
		if (out_nmass) fprintf(out_nmass,lonlatyearstr,"Lon","Lat","Year");
		if (out_nsources) fprintf(out_nsources,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","Lon","Lat","Year","dep","min",
				"imm","min-imm","fix","Total","Ndemand");
		if (out_npool && ifcentury) 
			fprintf(out_npool,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%10s\n","Lon","Lat","Year","VegN","LittVN",
			"LittSN","CwdN","MicroN","HumusN","SoilN","Total");
		if (out_nleach) fprintf(out_nleach,lonlatyearstr,"Lon","Lat","Year");
		if (out_age) fprintf(out_age,lonlatyearstr,"Lon","Lat","Year");
		if (out_nuptake) fprintf(out_nuptake,lonlatyearstr,"Lon","Lat","Year");
		if (out_anppn) fprintf(out_anppn,lonlatyearstr,"Lon","Lat","Year");
		if (out_vmaxnlim) fprintf(out_vmaxnlim,lonlatyearstr,"Lon","Lat","Year");
		// end GUESSN

		// GUESSN allometry
		if (out_allometry) fprintf(out_allometry,lonlatyearstr,"Lon","Lat","Year");
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
			// end GUESSN

			// GUESSN allometry
			if (out_allometry && pft.lifeform==TREE) fprintf(out_allometry,"%8s%8s%8s%8s%10s%8s%10s%8s%8s",(char*)pft.name,"N","Mfol","Mfroot","Mwood","Mactive","M","H","D");
			// end GUESSN

			pftlist.nextobj();
		}

		for(c=0;c<OUTPUT_MAXAGECLASS;c++)
			if (out_age) fprintf(out_age,"%7.1f ",(double)(c*estinterval+estinterval/2.0));

		// Print labels for "Total" columns

		if (out_cmass) fprintf(out_cmass,"%8s\n","Total");
		if (out_anpp) fprintf(out_anpp,"%8s\n","Total");
		if (out_lai) fprintf(out_lai,"%8s\n","Total");
		if (out_runoff) fprintf(out_runoff,"%8s\n","Total");
		if (out_dens) fprintf(out_dens,"%8s\n","Total");

		// GUESSN
		if (out_cton) fprintf(out_cton,"\n");
		if (out_nmass) fprintf(out_nmass,"%8s\n","Total");
		if (out_nleach) fprintf(out_nleach,"%8s\n","Total");	//(kgN/ha/yr)
		if (out_age) fprintf(out_age,"%8s\n","Total");	//(indiv/ha)
		if (out_nuptake) fprintf(out_nuptake,"%9s\n","Total");
		if (out_anppn) fprintf(out_anppn,"%9s\n","Total");
		if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8s\n","Total");
		if (out_allometry) fprintf(out_allometry,"\n");	//(indiv/ha)out_vmaxnlim
		// end GUESSN

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
		
	// If only yearly output between, say 1961 and 1990 is requred, use: 
	//	if (date.year>=nyear_spinup+60 && date.year<nyear_spinup+90) {

		if (date.year>=nyear_spinup) {

		lon=gridlist.getobj().lon;
		lat=gridlist.getobj().lat;

		cmass_stand=0.0;
		anpp_stand=0.0;
		lai_stand=0.0;
		runoff_stand=0.0;
		dens_stand=0.0;
		firert_stand=0.0;

		// GUESSN
		nmass_stand=0.0;
		vmaxnlim_stand=0.0;
		total_dens=0.0;
		nuptake_stand=0.0;
		anppn_stand=0.0;
		nleach_stand=0.0;
		densindiv_ageclass_stand=0.0;
		for (c=0;c<nclass;c++)
			stand_ageclass[c]=0.0;
		
		// end GUESSN

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
		if (out_age) fprintf(out_age,lonlatyeardatastr,lon,lat,date.year);
		if (out_nuptake) fprintf(out_nuptake,lonlatyeardatastr,lon,lat,date.year);
		if (out_anppn) fprintf(out_anppn,lonlatyeardatastr,lon,lat,date.year);
		if (out_vmaxnlim) fprintf(out_vmaxnlim,lonlatyeardatastr,lon,lat,date.year);
		// end GUESSN

		// GUESSN allometry
		if (out_allometry) fprintf(out_allometry,lonlatyeardatastr,lon,lat,date.year);
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


		for (int d=0;d<365;d++)
			dgpp[d]=0.0;

		// *** Loop through PFTs ***

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Standpft& standpft=stand.pft[pft.id];

			// Sum C biomass, NPP and LAI across patches and PFTs
			
			standpft.cmass_total=0.0;
			standpft.anpp_total=0.0;
			standpft.lai_total=0.0;
			standpft.densindiv_total = 0.0;
			standpft.greff_mort_total=0.0;
			standpft.nsapling_total=0.0;

			// GUESSN
			standpft.cton_leaf_avr=0.0;
			standpft.nmass_total=0.0;
			standpft.nuptake_total=0.0;
			standpft.anppn_total=0.0;
			standpft.vmaxnlim_avr=0.0;

			double nr_pft_indiv=0.0;
				// number of individuals of this pft. Used for avr calculation
			// end GUESSN

			// GUESSN allometry
			double allometry[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

			// Initialise age structure array

			if (vegmode==COHORT || vegmode==INDIVIDUAL)
				for (c=0;c<nclass;c++)
					standpft.densindiv_ageclass[c]=0.0;
	
			for (p=0;p<npatch;p++) {

				Patch& patch=stand[p];
				Vegetation& vegetation=patch.vegetation;
				Patchpft& patchpft=patch.pft[pft.id];

				standpft.nsapling_total+=patch.pft[pft.id].nsapling_yearly;
				
				if (patch.pft[pft.id].no_cohorts) {

					standpft.greff_mort_total+=patch.pft[pft.id].greff_mort_fraction/(patch.pft[pft.id].no_cohorts);
				}

				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv=vegetation.getobj();
					
					// guess2008 - alive check added
					if (indiv.id!=-1 && indiv.alive) { 
					
						if (indiv.pft.id==pft.id) {

							standpft.cmass_total+=indiv.cmass_leaf+
								indiv.cmass_root+indiv.cmass_sap+indiv.cmass_heart-indiv.cmass_debt;
							standpft.anpp_total+=indiv.anpp;
							standpft.lai_total+=indiv.lai;

							// GUESSN
							standpft.cton_leaf_avr+=indiv.cmass_leaf/indiv.nmass_leaf*indiv.densindiv;
							standpft.vmaxnlim_avr+=indiv.nopt*indiv.densindiv;
							vmaxnlim_stand+=indiv.nopt*indiv.densindiv;
							nr_pft_indiv+=indiv.densindiv;
							total_dens+=indiv.densindiv;
							standpft.nuptake_total+=indiv.fnuptake*indiv.ndemand_uptake;
							standpft.anppn_total+=indiv.ndemand;
							standpft.nmass_total+=indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap+
								indiv.nmass_heart+indiv.nmass_reserve;
							// end GUESSN

							if (vegmode==COHORT || vegmode==INDIVIDUAL) {

								// GUESSN allometry - only count trees with a trunk above a certain diameter 
								double diam=pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3);

								if (diam>0.03 && pft.lifeform==TREE) {
									// Number of individuals
									allometry[0]+=indiv.densindiv*10000.0;	
								
									// Total leaf C mass for these individuals
									allometry[1]+=(indiv.cmass_leaf/indiv.densindiv)*	// Leaf mass for one indiv
										(indiv.densindiv*10000.0);						// Number of individuals
								
									// Total fine root C mass
									allometry[2]+=(indiv.cmass_root/indiv.densindiv)*	
										(indiv.densindiv*10000.0);
								
									// Total wood C mass
									allometry[3]+=((indiv.cmass_sap+indiv.cmass_heart)/indiv.densindiv)*	
										(indiv.densindiv*10000.0);
								
									// Total active C mass
									allometry[4]+=((indiv.cmass_leaf+indiv.cmass_root)/indiv.densindiv)*	
										(indiv.densindiv*10000.0);
								
									// Total C mass
									allometry[5]+=((indiv.cmass_sap+indiv.cmass_heart+indiv.cmass_leaf+indiv.cmass_root)/indiv.densindiv)*
										(indiv.densindiv*10000.0);
								
									// Accumulated height of these individuals
									allometry[6]+=indiv.height*indiv.densindiv*10000.0;	
							
									// Accumulated diameter of these individuals
									allometry[7]+=(pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3))*indiv.densindiv*10000.0;
								}// end GUESSN
							
								// Age structure
								
								c=(int)(indiv.age/estinterval); // guess2008
								if (c<OUTPUT_MAXAGECLASS && indiv.pft.lifeform == TREE)
									standpft.densindiv_ageclass[c]+=indiv.densindiv*10000.0;

								// guess2008 - only count trees with a trunk above a certain diameter  
								if (diam>0.03 && pft.lifeform==TREE)
									standpft.densindiv_total+=indiv.densindiv; // indiv/m2
							}
						
						}

					} // alive


					vegetation.nextobj();
				}
			} // end of patch loop

			standpft.cmass_total/=(double)npatch;
			standpft.anpp_total/=(double)npatch;
			standpft.lai_total/=(double)npatch;
			standpft.densindiv_total/=(double)npatch;
			standpft.greff_mort_total/=(double)npatch;
			standpft.nsapling_total/=(double)npatch;
			standpft.nuptake_total/=(double)npatch;
			standpft.anppn_total/=(double)npatch;

			if (!negligible(nr_pft_indiv)) {
				standpft.cton_leaf_avr/=nr_pft_indiv;
				standpft.vmaxnlim_avr/=nr_pft_indiv;
			}
			else {
				standpft.cton_leaf_avr=0.0;
				standpft.vmaxnlim_avr=0.0;
			}
			standpft.nmass_total/=(double)npatch;
			for (c=0;c<nclass;c++) {
				standpft.densindiv_ageclass[c]/=(double)npatch;
				stand_ageclass[c]+=standpft.densindiv_ageclass[c];
				densindiv_ageclass_stand+=standpft.densindiv_ageclass[c];
			}

			for (int all=0;all<8;all++)
				allometry[all]/=(double)npatch;
			// end GUESSN

			// Update stand totals

			cmass_stand+=standpft.cmass_total;
			anpp_stand+=standpft.anpp_total;
			lai_stand+=standpft.lai_total;
			dens_stand+=standpft.densindiv_total;

			// GUESSN
			nmass_stand+=standpft.nmass_total;
			nuptake_stand+=standpft.nuptake_total;
			anppn_stand+=standpft.anppn_total;
			// end GUESSN
		
			// Print PFT sums to files

			if (out_cmass) fprintf(out_cmass,"%8.3f",standpft.cmass_total);
			if (out_anpp) fprintf(out_anpp,"%8.3f",standpft.anpp_total);
			if (out_lai) fprintf(out_lai,"%8.4f",standpft.lai_total);
			if (out_dens) fprintf(out_dens,"%8.4f",standpft.densindiv_total);

			// GUESSN
			if (out_cton) fprintf(out_cton,"%8.3f",standpft.cton_leaf_avr);
			if (out_nmass) fprintf(out_nmass,"%8.3f",standpft.nmass_total);
			if (out_nuptake) fprintf(out_nuptake,"%9.5f",standpft.nuptake_total);
			if (out_anppn) fprintf(out_anppn,"%9.5f",standpft.anppn_total);
			if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8.3f",standpft.vmaxnlim_avr);
			// end GUESSN
			
			// GUESSN allometry
			if (out_allometry && pft.lifeform==TREE)
				if (allometry[0]>0.0)
					fprintf(out_allometry,"%8d%8.0f%8.4f%8.4f%10.2f%8.4f%10.2f%8.2f%8.3f",
					standpft.pft.id,allometry[0],allometry[1]/allometry[0],
					allometry[2]/allometry[0],allometry[3]/allometry[0],allometry[4]/allometry[0],
					allometry[5]/allometry[0],allometry[6]/allometry[0],allometry[7]/allometry[0]);
				else
					fprintf(out_allometry,"%8d%8.0f%8.1f%8.1f%10.1f%8.1f%10.1f%8.1f%8.1f",
					(char*)standpft.pft.id,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0);

			if (has_CANIF_clim && date.year-100 >= nyear_spinup && (((WSITE==2 || WSITE==3 || WSITE==4) && pft.name == "TeBS") || (WSITE==0 && pft.name=="BNE") || ((WSITE==1 || WSITE==5 || WSITE==6 || WSITE==7) && pft.name=="TeNE")))
				dprintf("Year %d height %g dens %g anpp %g C:N %g\n",date.year,allometry[6]/allometry[0],allometry[0],standpft.anpp_total,standpft.cton_leaf_avr);

			// end GUESSN

			// Graphical output every 10 years
			// (Windows shell only - "plot" statements have no effect otherwise)
			

			//if (!(date.year%10)) {
				plot("cmass",pft.name,date.year,stand.pft[pft.id].cmass_total);
				plot("anpp",pft.name,date.year,stand.pft[pft.id].anpp_total);
				plot("lai",pft.name,date.year,stand.pft[pft.id].lai_total);
				plot("pft Vmax N limitation",pft.name,date.year,standpft.vmaxnlim_avr);
				if (standpft.cton_leaf_avr>10.0)
					plot("C:N ratios pft",pft.name,date.year,standpft.cton_leaf_avr);
			//}

			
			// if (!(date.year%10)) {
				if (pft.lifeform==TREE && vegmode==COHORT || pft.lifeform==TREE && vegmode==INDIVIDUAL) {
					plot("growth_efficiency_mort;averaged all cohorts",pft.name,date.year,stand.pft[pft.id].greff_mort_total);
					plot("tree_density/ha2",pft.name,date.year,standpft.densindiv_total*10000.0);
				    plot("no_saplings/patch",pft.name,date.year,stand.pft[pft.id].nsapling_total);
				}
			//}

				// FACE DAVID 
				double ORNL_amb[10] = {800.0,850.0,1000.0,1025.0,900.0,975.0,800.0,775.0,625.0,650.0};
				double ORNL_ele[10] = {1025.0,1000.0,1200.0,1300.0,1200.0,1125.0,1080.0,875.0,700.0,700.0};
				double Duke_amb[10] = {1020.0,1248.0,1330.0,1121.0,718.0,1046.0,1142.0,996.0,0.0,0.0};
				double Duke_ele[10] = {1315.0,1615.0,1711.0,1435.0,1017.0,1424.0,1537.0,1367.0,0.0,0.0};

				if (date.year >= 2096 && has_FACE_clim) { 
					plot("anpp output",pft.name,date.year,stand.pft[pft.id].anpp_total);
					if (!ifduke) {
						if (FACE_ring == 1)
							plot("anpp output","Real Amb",date.year,ORNL_amb[date.year-2096]/1000.0);
						else 
							plot("anpp output","Real Ele",date.year,ORNL_ele[date.year-2096]/1000.0);
					}
					else {
						if (FACE_ring == 1)
							plot("anpp output","Real Amb",date.year,Duke_amb[date.year-2096]/1000.0);
						else 
							plot("anpp output","Real Ele",date.year,Duke_ele[date.year-2096]/1000.0);
					}
				}

			pftlist.nextobj();
		
		} // *** End of PFT loop ***


		flux_veg=flux_soil=flux_fire=flux_est=0.0;

		// guess2008 - carbon pools
		c_litter=c_fast=c_slow=0.0;

		// GUESSN
		surfsoillitterc=surfsoillittern=cwdc=cwdn=microc=micron=humusc=humusn=centuryc=centuryn=n_litter=0.0;
		andep_stand=anmin_stand=animm_stand=anfix_stand=nsupply_stand=ndemand_stand=0.0;
		// end GUESSN

		// Sum C fluxes, dead C pools and runoff across patches

		for (p=0;p<npatch;p++) {
			flux_veg+=stand[p].fluxes.acflux_veg/(double)npatch;
			flux_soil+=stand[p].fluxes.acflux_soil/(double)npatch;
			flux_fire+=stand[p].fluxes.acflux_fire/(double)npatch;
			flux_est+=stand[p].fluxes.acflux_est/(double)npatch;

			c_fast+=stand[p].soil.cpool_fast/(double)npatch;
			c_slow+=stand[p].soil.cpool_slow/(double)npatch;

			// GUESSN
			andep_stand+=stand[p].soil.ndep_annual/(double)npatch*10000.0; // convert from m2 to ha
			anmin_stand+=stand[p].soil.nmin_annual/(double)npatch*10000.0; // convert from m2 to ha
			animm_stand+=stand[p].soil.nimmob_annual/(double)npatch*10000.0; // convert from m2 to ha
			anfix_stand+=stand[p].soil.N_fix/(double)npatch*10000.0; // convert from m2 to ha
			nleach_stand+=stand[p].soil.nleach_annual/(double)npatch*10000.0;	// convert from m2 to ha
			nsupply_stand+=stand[p].nsupply/(double)npatch*10000.0;			// convert from m2 to ha
			ndemand_stand+=stand[p].ndemand/(double)npatch*10000.0;			// convert from m2 to ha
			
			for (int r=0;r<NSOMPOOL;r++) {
				if (stand[p].soil.sompool[r].nmass > 0.0) {
					if(r==SURFMETA||r==SURFSTRUCT||r==SOILMETA||r==SOILSTRUCT){
						surfsoillitterc+=stand[p].soil.sompool[r].cmass/(double)npatch;
						surfsoillittern+=stand[p].soil.sompool[r].nmass/(double)npatch;
					}
					else if (r==SURFCWD) {
						cwdc+=stand[p].soil.sompool[r].cmass/(double)npatch;
						cwdn+=stand[p].soil.sompool[r].nmass/(double)npatch;
					}
					else if (r==SURFMICRO||r==SOILMICRO) {
						microc+=stand[p].soil.sompool[r].cmass/(double)npatch;
						micron+=stand[p].soil.sompool[r].nmass/(double)npatch;
					}
					else if (r==SURFHUMUS){
						humusc+=stand[p].soil.sompool[r].cmass/(double)npatch;
						humusn+=stand[p].soil.sompool[r].nmass/(double)npatch;
					}
					else {
						centuryc+=stand[p].soil.sompool[r].cmass/(double)npatch;
						centuryn+=stand[p].soil.sompool[r].nmass/(double)npatch;
					}
				}
			}
			// end GUESSN



			// Sum all litter
			for (int q=0;q<npft;q++) {
				Patchpft& pft=stand[p].pft[q];

				c_litter+=(pft.litter_leaf+pft.litter_root+pft.litter_wood+pft.litter_repr)/(double)npatch;
				n_litter+=(pft.nmass_litter_leaf+pft.nmass_litter_root+pft.nmass_litter_wood)/(double)npatch;
			}

			runoff_stand+=stand[p].arunoff/(double)npatch;
			
			// Fire return time
			if (!iffire || stand[p].fireprob < 0.001)
				firert_stand+=1000.0/(double)npatch; // Set a limit of 1000 years
			else	
				firert_stand+=(1.0/stand[p].fireprob)/(double)npatch;


			// Monthly output variables
			
			for (m=0;m<12;m++) {
				maet[m] += stand[p].maet[m]/(double)npatch;
				mpet[m] += stand[p].mpet[m]/(double)npatch;
				mevap[m] += stand[p].mevap[m]/(double)npatch;
				mintercep[m] += stand[p].mintercep[m]/(double)npatch;
				mrunoff[m] += stand[p].mrunoff[m]/(double)npatch;
				mrh[m] += stand[p].fluxes.mcflux_soil[m]/(double)npatch;
				
				mwcont_upper[m] += stand[p].soil.mwcont[m][0]/(double)npatch;
				mwcont_lower[m] += stand[p].soil.mwcont[m][1]/(double)npatch;

				// guess2008 - average across stands to get mgpp and mra here. 
				mgpp[m] += stand[p].fluxes.mcflux_gpp[m]/(double)npatch;
				mra[m] += stand[p].fluxes.mcflux_ra[m]/(double)npatch;

			}

			for (int ddddd=0;ddddd<365;ddddd++){
				dgpp[ddddd] += stand[p].fluxes.dcflux_gpp[ddddd]/(double)npatch;
			}




			// Calculate monthly NPP and LAI

			Vegetation& vegetation=stand[p].vegetation;

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();
				
				// guess2008 - alive check added
				if (indiv.id!=-1 && indiv.alive) { 

					for (m=0;m<12;m++) {
						mlai[m] += indiv.mlai[m]/(double)npatch;
					}

				} // alive?

				vegetation.nextobj();

			} // while/vegetation loop


		} // patch loop



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

		// GUESSN 
		if (!negligible(total_dens))
			vmaxnlim_stand/=total_dens;
		else
			vmaxnlim_stand=0.0;

		// Print stand totals to files

		if (out_cmass) fprintf(out_cmass,"%8.3f\n",cmass_stand);
		if (out_anpp) fprintf(out_anpp,"%8.3f\n",anpp_stand);
		if (out_lai) fprintf(out_lai,"%8.4f\n",lai_stand);
		if (out_runoff) fprintf(out_runoff,"%8.1f\n",runoff_stand);
		if (out_dens) fprintf(out_dens,"%8.4f\n",dens_stand);
		if (out_firert) fprintf(out_firert,"%8.1f\n",firert_stand);

		// GUESSN
		if (out_cton) fprintf(out_cton,"\n");
		if (out_nmass) fprintf(out_nmass,"%8.3f\n",nmass_stand);
		if (out_nleach) fprintf(out_nleach,"%8.3f\n",nleach_stand);
		if (out_nsources) fprintf(out_nsources,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
			andep_stand,anmin_stand,animm_stand,anmin_stand-animm_stand,
			anfix_stand,andep_stand+anmin_stand-animm_stand+anfix_stand,ndemand_stand);
		if (out_nuptake) fprintf(out_nuptake,"%9.5f\n",nuptake_stand);
		if (out_anppn) fprintf(out_anppn,"%9.5f\n",anppn_stand);
		if (out_vmaxnlim) fprintf(out_vmaxnlim,"%8.3f\n",vmaxnlim_stand);

		for (c=0;c<nclass;c++)
			if (out_age) fprintf(out_age,"%8.1f",stand_ageclass[c]);
		if (out_age) fprintf(out_age,"%9.1f\n",densindiv_ageclass_stand);
		// end GUESSN

		// GUESSN allometry
		if (out_allometry) fprintf(out_allometry,"\n");

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

		for (int dddd=0;dddd<365;dddd++) {

			if (out_dgpp) fprintf(out_dgpp,"%8.1f%8.1f%8d%8d%9.5f\n",lon,lat,date.year,dddd,dgpp[dddd]);//lon,lat,date.year,dddd);

		}



		// Graphical output every 10 years
		// (Windows shell only - no effect otherwise)

		if (!(date.year%10)) {
		//	plot("fluxes","flux_veg",date.year,flux_veg);
		//	plot("fluxes","flux_soil",date.year,flux_soil);
		//	plot("fluxes","flux_fire",date.year,flux_fire);
		//	plot("fluxes","flux_est",date.year,flux_est);
		//	plot("fluxes","NEE",date.year,flux_veg+flux_soil+flux_fire+flux_est);
			if (!ifcentury) { // GUESSN
				plot("soilc","slow",date.year,stand[0].soil.cpool_slow);
				plot("soilc","fast",date.year,stand[0].soil.cpool_fast);
			}
		}

		plot("N fixation (kgN/ha/yr)","Soil N fix",date.year,anfix_stand);
		plot("N min-immob (kgN/ha/yr)","N",date.year,anmin_stand-animm_stand);
		
		plot("N demand/supply (kgN/ha/yr)","N supply",date.year,nsupply_stand);
		plot("N demand/supply (kgN/ha/yr)","N demand",date.year,ndemand_stand);


		// Write fluxes to file

		if (out_cflux) fprintf(out_cflux,"%8.3f%8.3f%8.3f%8.3f%10.5f\n",flux_veg,flux_soil,flux_fire,
			flux_est,flux_veg+flux_soil+flux_fire+flux_est);

		// guess2008 - output carbon pools
		if (out_cpool)
			if (!ifcentury)
				fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%10.4f\n",cmass_stand,c_litter,c_fast,
					c_slow,cmass_stand+c_litter+c_fast+c_slow);
			else
				fprintf(out_cpool,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%10.3f\n",cmass_stand,c_litter,
					surfsoillitterc,cwdc,microc,humusc,centuryc,
					cmass_stand+c_litter+surfsoillitterc+cwdc+microc+humusc+centuryc);

		// GUESSN
		if (out_npool && ifcentury)
				fprintf(out_npool,"%8.3f%8.3f%8.4f%8.4f%8.4f%8.3f%8.3f%10.3f\n",nmass_stand,n_litter,
					surfsoillittern,cwdn,micron,humusn,centuryn,
					nmass_stand+n_litter+surfsoillittern+cwdn+micron+humusn+centuryn);
		// end GUESSN

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
								standpft.densindiv_ageclass[c]);
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

		// GUESSN
		if (out_cton) fclose(out_cton);
		if (out_nmass) fclose(out_nmass);
		if (out_nsources) fclose(out_nsources);
		if (out_npool && ifcentury) fclose(out_npool);
		if (out_nleach) fclose(out_nleach);
		if (out_age) fclose(out_age);
		if (out_nuptake) fclose(out_nuptake);
		if (out_anppn) fclose(out_anppn);
		if (out_vmaxnlim) fclose(out_vmaxnlim);
		// end GUESSN

		// GUESSN allometry
		if (out_allometry) fclose(out_allometry);
		// end GUESSN

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

		if (out_dgpp) fclose(out_dgpp);
	}

	// Clean up

	gridlist.killall();
}


#endif // USE_CRU


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Galloway, J. N., F. J. Dentener, D. G. Capone, E. W. Boyer, R. W. Howarth, S. P. Seitzinger,
// G. P. Asner, C. Cleveland, P. Green, E. Holland, D. M. Karl, A. F. Michaels, J. H. Porter, 
// A. Townsend, and C. Vörösmarty. 2004.
// Nitrogen Cycles: Past, Present and Future. Biogeochemistry 70: 153-226.
