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


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
vegmodetype vegmode; // vegetation mode (population, cohort or individual)
int npatch; // number of patches in each stand (should always be 1 in population mode); cropland stands always have 1 patch
double patcharea; // patch area (m2) (individual and cohort mode only)
bool ifdailynpp; // whether NPP calculations performed daily (alt: monthly)
bool ifdailydecomp;
	// whether soil decomposition calculations performed daily (alt: monthly)
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
bool iffast;
bool ifcdebt;

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

Pftlist pftlist;

////////////////////////////////////////////////////////////////////////////////
// Implementation of PhotosynthesisResult member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisResult::serialize(ArchiveStream& arch) {
	arch & agd_g
		& adtmm
		& rd_g
		& vm
		& je;
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
		& mgtemp
		& last_mgtemp
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
		& doneday;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Fluxes member functions
////////////////////////////////////////////////////////////////////////////////


void Fluxes::serialize(ArchiveStream& arch) {
	arch & acflux_veg
		& acflux_fire
		& acflux_soil
		& acflux_est
		& acflux_harvest 
		& dcflux_soil
		& mcflux_soil
		& mcflux_veg
		& dcflux_veg
		& mcflux_gpp
		& mcflux_ra
		& miso
		& mmon;
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
		& mgtemp
		& last_mgtemp
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
		& gcbase_wstress
		& temp_wstress
		& par_wstress
		& daylength_wstress
		& co2_wstress
		& nday_wstress
		& fpar_grass_wstress
		& gpterm_wstress
		& supply
		& supply_leafon
		& fuptake
		& wstress
		& wstress_day
		& harvested_products_slow
		& phot_wstress;
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
		& aevap
		& aintercep
		& arunoff
		& apet
		& eet_net_veg
		& demand
		& demand_day
		& demand_leafon
		& fpc_rescale
		& maet
		& mevap
		& mintercep
		& mrunoff
		& mpet;
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
		& active;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}


	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND) {
		npatchL=1;
	}
	else if(landcover==NATURAL || landcover==FOREST) {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	for (p=0;p<npatchL;p++) {
		createobj(*this,gc.soiltype);
	}

	first_year=date.year;
}

double Stand::get_gridcell_fraction() const {
	return frac*gridcell.landcoverfrac[landcover];
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
		& frac;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp=0.0;
	fpc=0.0;
	densindiv=0.0;
	cmass_leaf=0.0;
	cmass_root=0.0;
	cmass_sap=0.0;
	cmass_heart=0.0;
	cmass_debt=0.0;
	phen=0.0;
	aphen=0.0;
	deltafpc=0.0;
	fpar_wstress=0.0;
	assim=0.0;

	// guess2008 - additional initialisation
	age=0.0;
	fpar=0.0;
	aphen_raingreen=0;
	intercep=0.0;
	phen_mean=0.0;
	temp_wstress = 0.0;
	par_wstress = 0.0;
	daylength_wstress = 0.0;
	co2_wstress = 0.0; 
	nday_wstress = 0; 
	wstress = false;
	lai = 0.0;
	lai_layer = 0.0;
	lai_indiv = 0.0;
	alive = false;

	int m;
	for (m=0;m<12;m++) {
		mnpp[m]=mlai[m]=mgpp[m]=mra[m]=0.0;
	}

	// bvoc
	monstor=0.;
	iso=0.;
	mon=0.;
	aiso=0.;
	amon=0.;
	fvocseas=1.;
	dtr_wstress=0.;
	eet_wstress=0.;
	agdd5_wstress=0.;
	rad_wstress=0.;		
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
		& mnpp
		& mlai
		& mgpp
		& mra
		& fpar_wstress
		& fpar_leafon
		& lai_leafon_layer
		& demand
		& demand_leafon
		& supply
		& supply_leafon
		& intercep
		& phen_mean
		& temp_wstress 
		& par_wstress 
		& daylength_wstress 
		& co2_wstress 
		& nday_wstress 
		& wstress 
		& alive 
		& iso 
		& mon 
		& aiso 
		& amon 
		& monstor 
		& fvocseas 
		& dtr_wstress 
		& eet_wstress 
		& agdd5_wstress 
		& rad_wstress; 
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellpft member functions
////////////////////////////////////////////////////////////////////////////////


void Gridcellpft::serialize(ArchiveStream& arch) {
	arch & addtw;
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
		& LC_updated;

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
