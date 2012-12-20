///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.cpp
/// \brief LPJ-GUESS Combined Modular Framework
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "driver.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
vegmodetype vegmode; // vegetation mode (population, cohort or individual)
int npatch; // number of patches in each stand (should always be 1 in population mode); cropland stands always have 1 patch
double patcharea; // patch area (m2) (individual and cohort mode only)
bool ifbgestab; // whether background establishment enabled (individual, cohort mode)
bool ifsme;
	// whether spatial mass effect enabled for establishment (individual, cohort mode)
bool ifstochestab; // whether establishment stochastic (individual, cohort mode)
bool ifstochmort; // whether mortality stochastic (individual, cohort mode)
bool iffire; // whether fire enabled
bool ifdisturb;
	// whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
bool ifcalcsla; // whether SLA calculated from leaf longevity (alt: prescribed)
int estinterval; // establishment interval in cohort mode (years)
double distinterval;
	// generic patch-destroying disturbance interval (individual, cohort mode)
int npft; // number of possible PFTs
bool ifcdebt;

/// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
bool ifcentury;
/// whether plant growth limited by available nitrogen	
bool ifnlim;
/// number of years to allow spinup without nitrogen limitation	
int freenyears;
/// fraction of nitrogen relocated by plants from roots and leaves
double nrelocfrac;
/// whether to allow nitrogen leaching	
bool ifleachn;
/// first term in nitrogen fixation eqn
double nfix_a;
/// second term in nitrogen fixation eqn
double nfix_b;
/// whether nitrogen deposition data available from a file	
bool ifndepdata;


// guess2008 - new inputs from the .ins file
bool ifsmoothgreffmort;				// smooth growth efficiency mortality
bool ifdroughtlimitedestab;			// whether establishment affected by growing season drought
bool ifrainonwetdaysonly;			// rain on wet days only (1, true), or a little every day (0, false); 
// bvoc
bool ifbvoc; // BVOC calculations included

wateruptaketype wateruptake;

bool run_landcover;
bool run[NLANDCOVERTYPES];
bool lcfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;				// If a slow harvested product pool is included in patchpft.
int nyear_spinup;	

/// Solving Century SOM pools 
/// years at which to begin documenting for calculation of Century equilibrium
int solvesomcent_beginyr;
/// years at which to end documentation and start calculation of Century equilibrium
int solvesomcent_endyr;

xtring state_path;
bool restart;
bool save_state;
int state_year;

Pftlist pftlist;

////////////////////////////////////////////////////////////////////////////////
// Implementation of PhotosynthesisResult member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisResult::serialize(ArchiveStream& arch) {
	arch & agd_g
		& adtmm
		& rd_g
		& vm
		& je
		& nmass_term
		& vmaxnlim;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Climate member functions
////////////////////////////////////////////////////////////////////////////////


void Climate::serialize(ArchiveStream& arch) {
	arch & temp
		& rad
		& par
		& prec
		& daylength
		& co2
		& lat
		& insol
		& instype
		& eet
		& mtemp
		& mtemp_min20
		& mtemp_max20
		& mtemp_max
		& gdd5
		& agdd5 
		& chilldays
		& ifsensechill
		& gtemp
		& dtemp_31
		& mtemp_min_20
		& mtemp_max_20
		& mtemp_min
		& atemp_mean
		& temp_mean
		& par_mean
		& co2_mean
		& daylength_mean
		& sinelat
		& cosinelat
		& qo & u & v & hh & sinehh
		& daylength_save
		& doneday
		& andep
		& dndep
		& anfert
		& dnfert;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Fluxes member functions
////////////////////////////////////////////////////////////////////////////////

Fluxes::Fluxes(Patch& p) 		
  : patch(p), 

    annual_fluxes_per_pft(npft, std::vector<double>(NPERPFTFLUXTYPES)) {

	reset();
}

void Fluxes::reset() {
	for (size_t i = 0; i < annual_fluxes_per_pft.size(); ++i) {
		std::fill_n(annual_fluxes_per_pft[i].begin(), int(NPERPFTFLUXTYPES), 0);
	}

	for (int m = 0; m < 12; ++m) {
		std::fill_n(monthly_fluxes_pft[m], int(NPERPFTFLUXTYPES), 0);

		std::fill_n(monthly_fluxes_patch[m], int(NPERPATCHFLUXTYPES), 0);
	}
}

void Fluxes::serialize(ArchiveStream& arch) {
	arch & annual_fluxes_per_pft 
		& monthly_fluxes_patch
		& monthly_fluxes_pft;
}

void Fluxes::report_flux(PerPFTFluxType flux_type, int pft_id, double value) {
	annual_fluxes_per_pft[pft_id][flux_type] += value;
	monthly_fluxes_pft[date.month][flux_type] += value;
}

void Fluxes::report_flux(PerPatchFluxType flux_type, double value) {
	monthly_fluxes_patch[date.month][flux_type] += value;
}

double Fluxes::get_monthly_flux(PerPFTFluxType flux_type, int month) const {
	return monthly_fluxes_pft[month][flux_type];
}

double Fluxes::get_monthly_flux(PerPatchFluxType flux_type, int month) const {
	return monthly_fluxes_patch[month][flux_type];
}

double Fluxes::get_annual_flux(PerPFTFluxType flux_type, int pft_id) const {
	return annual_fluxes_per_pft[pft_id][flux_type];
}

double Fluxes::get_annual_flux(PerPFTFluxType flux_type) const {
	double sum = 0;
	for (size_t i = 0; i < annual_fluxes_per_pft.size(); ++i) {
		sum += annual_fluxes_per_pft[i][flux_type];
	}
	return sum;
}

double Fluxes::get_annual_flux(PerPatchFluxType flux_type) const {
	double sum = 0;
	for (int m = 0; m < 12; ++m) {
		sum += monthly_fluxes_patch[m][flux_type];
	}
	return sum;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Vegetation member functions
////////////////////////////////////////////////////////////////////////////////


void Vegetation::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		arch & nobj;

		for (unsigned int i = 0; i < nobj; i++) {
			Individual& indiv = (*this)[i];
			arch & indiv.pft.id
				& indiv;
		}
	}
	else {
		killall();
		unsigned int number_of_individuals;
		arch & number_of_individuals;

		for (unsigned int i = 0; i < number_of_individuals; i++) {
			int pft_id;
			arch & pft_id;
			Individual& indiv = createobj(pftlist[pft_id], *this);
			arch & indiv;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of LitterSolveSOM member functions
////////////////////////////////////////////////////////////////////////////////


void LitterSolveSOM::serialize(ArchiveStream& arch) {
	for (int p = 0; p<NSOMPOOL; p++) {
		arch & clitter[p]
		     & nlitter[p];
	} 
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soil member functions
////////////////////////////////////////////////////////////////////////////////


void Soil::serialize(ArchiveStream& arch) {
	arch & wcont
		& awcont
		& wcont_evap
		& dwcontupper
		& mwcontupper
		& snowpack
		& runoff
		& temp
		& dtemp
		& mtemp
		& gtemp
		& cpool_slow
		& cpool_fast
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean
		& alag
		& exp_alag
		& mwcont
		& dwcontlower
		& mwcontlower
		// probably shouldn't need to serialize these
		& rain_melt
		& max_rain_melt
		& percolate;

		for (int i = 0; i<NSOMPOOL; i++) {
			arch & sompool[i];
		} 

	arch & dperc		
		& orgleachfrac
		& nmass_avail		
		& anmin			
		& animmob			
		& aminleach		
		& aorgleach					
		& anfix
		& anfix_calc
		& anfix_mean
		& solvesomcent_beginyr
		& solvesomcent_endyr
		& solvesom
		& fnuptake_mean
		& morgleach_mean
		& mminleach_mean; 
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patchpft member functions
////////////////////////////////////////////////////////////////////////////////


void Patchpft::serialize(ArchiveStream& arch) {
	arch & anetps_ff
		& wscal
		& wscal_mean
		& anetps_ff_est
		& anetps_ff_est_initial
		& wscal_mean_est
		& phen
		& aphen
		& establish
		& nsapling
		& litter_leaf
		& litter_root
		& litter_wood
		& litter_repr
		& gcbase
		& gcbase_day
		& wsupply
		& wsupply_leafon
		& fwuptake
		& wstress
		& wstress_day
		& harvested_products_slow
		& nmass_litter_leaf
		& nmass_litter_root
		& nmass_litter_wood
		& harvested_products_slow_nmass;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patch member functions
////////////////////////////////////////////////////////////////////////////////


void Patch::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}
	}
	else {
		pft.killall();
				
		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			pft.createobj(pftlist[i]);
			arch & pft[i];
		}
	}

	arch & vegetation
		& soil
		& fluxes
		& fpar_grass
		& fpar_ff
		& par_grass_mean
		& nday_growingseason
		& fpc_total
		& disturbed
		& age
		& fireprob
		& growingseasondays
		& intercep
		& aaet
		& aaet_5
		& aevap
		& aintercep
		& arunoff
		& apet
		& eet_net_veg
		& wdemand
		& wdemand_day
		& wdemand_leafon
		& fpc_rescale
		& maet
		& mevap
		& mintercep
		& mrunoff
		& mpet
		& fnuptake
		& ndemand;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Standpft member functions
////////////////////////////////////////////////////////////////////////////////


void Standpft::serialize(ArchiveStream& arch) {
	arch & cmass_repr
		& anetps_ff_max
		& gpterm
		& assim_term
		& fpc_total
		& active
		& cmass_repr_nuptake;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL;

	for(p=0; p<pftlist. nobj;p++) {
		pft.createobj(pftlist[p]);
	}


	if(landcover == CROPLAND || landcover == PASTURE || landcover == URBAN || landcover == PEATLAND) {
		npatchL=1;
	}
	else if(landcover == NATURAL || landcover == FOREST) {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	for (p=0; p<npatchL; p++) {
		createobj(*this,gc.soiltype);
	}

	first_year = date.year;
	seed = 12345678;
}

double Stand::get_gridcell_fraction() const {
	return frac * gridcell.landcoverfrac[landcover];
}

double Stand::get_landcover_fraction() const {
	return frac;
}

void Stand::set_landcover_fraction(double fraction) {
	frac = fraction;
}

void Stand::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}

		arch & nobj;
		for (unsigned int k = 0; k < nobj; k++) {
			arch & (*this)[k];
		}
	}
	else {
		pft.killall();
		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			Standpft& standpft = pft.createobj(pftlist[i]);
			arch & standpft;
		}

		killall();
		unsigned int npatch;
		arch & npatch;
		for (unsigned int k = 0; k < npatch; k++) {
			Patch& patch = createobj(*this, gridcell.soiltype);
			arch & patch;
		}
	}

	arch & first_year
		& frac
		& seed;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp = 0.0;
	fpc = 0.0;
	densindiv = 0.0;
	cmass_leaf = 0.0;
	cmass_root = 0.0;
	cmass_sap = 0.0;
	cmass_heart = 0.0;
	cmass_debt = 0.0;
	phen = 0.0;
	aphen = 0.0;
	deltafpc = 0.0;
	assim = 0.0;
	resp = 0.0;
	assim_term = 0.0;

	nmass_leaf = 0.0;
	nmass_root = 0.0;
	nmass_sap = 0.0;
	nmass_heart = 0.0;
	cton_leaf_aopt = 0.0;
	cton_leaf_aavr = 0.0;
	cmass_veg = 0.0;
	nmass_veg = 0.0;

	nactive = 0.0;
	nstore_leaf = 0.0;
	nstore_root = 0.0;
	nstore_labile = 0.0;
	ndemand = 0.0;
	fnuptake = 1.0;
	anuptake = 0.0;
	max_n_storage = 0.0;
	scale_n_storage = 0.0;

	nstress = false;

	leafndemand = 0.0;
	rootndemand = 0.0;
	sapndemand = 0.0;
	storendemand = 0.0;
	for (int c=0; c<3; c++) {
		fndemand[c] = 0.0;
	}
	leafndemand_store = 0.0;
	rootndemand_store = 0.0;

	// additional initialisation
	age = 0.0;
	fpar = 0.0;
	aphen_raingreen = 0;
	wdemand = 0.0;
	wsupply = 0.0;
	intercep = 0.0;
	phen_mean = 0.0;
	wstress = false;
	lai = 0.0;
	lai_layer = 0.0;
	lai_indiv = 0.0;
	alive = false;

	int m;
	for (m=0; m<12; m++) {
		mlai[m]=0.0;
	}

	// bvoc
	monstor = 0.;
	iso = 0.;
	mon = 0.;
	fvocseas = 1.;
}

void Individual::serialize(ArchiveStream& arch) {
	arch & cmass_leaf
		& cmass_root
		& cmass_sap 
		& cmass_heart
		& cmass_debt
		& fpc
		& fpar
		& densindiv
		& phen
		& aphen
		& aphen_raingreen
		& assim
		& resp
		& anpp
		& aet
		& ltor
		& height
		& crownarea
		& deltafpc
		& wscal_mean
		& boleht
		& lai
		& lai_layer
		& lai_indiv
		& greff_5
		& age
		& mlai
		& fpar_leafon
		& lai_leafon_layer
		& wdemand
		& wdemand_leafon
		& wsupply
		& wsupply_leafon
		& intercep
		& phen_mean
		& wstress 
		& alive 
		& iso 
		& mon 
		& monstor 
		& fvocseas 
		& nmass_leaf
		& nmass_root
		& nmass_sap
		& nmass_heart
		& nactive
		& nstore_leaf
		& nstore_root
		& nstore_labile
		& ndemand
		& fnuptake
		& anuptake
		& max_n_storage
		& scale_n_storage
		& avmaxnlim
		& cton_leaf_aopt
		& cton_leaf_aavr
		& cmass_veg
		& nmass_veg

		& nstress
		& leafndemand
		& rootndemand
		& sapndemand
		& storendemand
		& fndemand
		& leafndemand_store
		& rootndemand_store
		
		& assim_term
		& nday_leafon;
}

void Individual::report_flux(Fluxes::PerPFTFluxType flux_type, double value) {
	if (alive) {
		vegetation.patch.fluxes.report_flux(flux_type, pft.id, value);
	}
}

void Individual::report_flux(Fluxes::PerPatchFluxType flux_type, double value) {
	if (alive) {
		vegetation.patch.fluxes.report_flux(flux_type, value);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellpft member functions
////////////////////////////////////////////////////////////////////////////////


void Gridcellpft::serialize(ArchiveStream& arch) {
	arch & addtw
		 & Km;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcell member functions
////////////////////////////////////////////////////////////////////////////////

double Gridcell::get_lon() const {
	return lon;
}

double Gridcell::get_lat() const {
	return lat;
}

void Gridcell::set_coordinates(double longitude, double latitude) {
	lon = longitude;
	lat = latitude;
}

void Gridcell::serialize(ArchiveStream& arch) {

	arch & climate
		& landcoverfrac
		& landcoverfrac_old
		& LC_updated
		& seed;

	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}

		arch & nobj;
		for (unsigned int s = 0; s < nobj; s++) {
			arch & (*this)[s].landcover
				& (*this)[s];
		}
	}
	else {
		pft.killall();

		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			pft.createobj(pftlist[i]);
			arch & pft[i];
		}

		killall();
		unsigned int number_of_stands;
		arch & number_of_stands;
				
		for (unsigned int s = 0; s < number_of_stands; s++) {
			landcovertype landcover;
			arch & landcover;
			createobj(*this, landcover);
			arch & (*this)[s];
		}
	}
}

void Sompool::serialize(ArchiveStream& arch) {
	arch & cmass
		& nmass
		& cdec 
		& ndec 
		& delta_cmass
		& delta_nmass
		& ligcfrac
		& fracremain
		& ntoc
		& litterme
		& fireresist
		& mfracremain_mean;
}
