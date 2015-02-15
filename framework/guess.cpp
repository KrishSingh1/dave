///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.cpp
/// \brief LPJ-GUESS Combined Modular Framework
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "config.h"
#include "guess.h"
#include "landcover.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
int npft; // number of possible PFTs
int nst;
int nst_lc[NLANDCOVERTYPES];

StandTypelist stlist;
Pftlist pftlist;

// emission ratios from fire (NH3, NO, NO2, N2O, N2) Levine et al. 1996

const double Fluxes::NH3_FIRERATIO = 0.236;
const double Fluxes::NO_FIRERATIO  = 0.303;
const double Fluxes::NO2_FIRERATIO = 0.076;
const double Fluxes::N2O_FIRERATIO = 0.035;
const double Fluxes::N2_FIRERATIO  = 0.350;


////////////////////////////////////////////////////////////////////////////////
// Implementation of MassBalance member functions
////////////////////////////////////////////////////////////////////////////////


void MassBalance::serialize(ArchiveStream& arch) {
	arch & start_year
		& ccont_zero
		& cflux_zero
		& ncont_zero
		& nflux_zero
		& cflux
		& nflux;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of PhotosynthesisResult member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisResult::serialize(ArchiveStream& arch) {
	arch & agd_g
		& adtmm
		& rd_g
		& vm
		& je
		& nactive_opt
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
		& lon
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
		& sinelat
		& cosinelat
		& qo & u & v & hh & sinehh
		& daylength_save
		& doneday
		& andep
		& dndep
		& dprec_10
		& sprec_2
		& maxtemp
		& testday_temp
		& testday_prec
		& coldestday
		& adjustlat
		& mtemp_20
		& mprec_20
		& mpet_20
		& mprec_pet_20
		& mprec_petmin_20
		& mprec_petmax_20
		& mtemp20
		& mprec20
		& mpet20
		& mprec_pet20
		& mprec_petmin20
		& mprec_petmax20
		& seasonality
		& prec_seasonality
		& var_prec
		& var_temp
		& aprec;
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

	for (int d = 0; d < 365; ++d) {
		std::fill_n(daily_fluxes_pft[d], int(NPERPFTFLUXTYPES), 0);
		std::fill_n(daily_fluxes_patch[d], int(NPERPATCHFLUXTYPES), 0);
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
	daily_fluxes_pft[date.day][flux_type] += value;	//Var = value ???
}

void Fluxes::report_flux(PerPatchFluxType flux_type, double value) {
	monthly_fluxes_patch[date.month][flux_type] += value;
	daily_fluxes_patch[date.day][flux_type] += value;
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
	arch & clitter
		& nlitter;
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
		& ninput
		& anmin			
		& animmob			
		& aminleach		
		& aorgleach					
		& anfix
		& anfix_calc
		& anfix_mean
		& snowpack_nmass
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
		& litter_sap
		& litter_heart
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
		& nmass_litter_sap
		& nmass_litter_heart
		& harvested_products_slow_nmass
		& swindow
		& water_deficit_y;
	if(pft.landcover==CROPLAND)
		arch & *cropphen;
		
}

cropphen_struct* Patchpft::get_cropphen() {
	if (pft.landcover != CROPLAND) {
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	}
	return cropphen;
}

cropphen_struct* Patchpft::set_cropphen() {
	if (pft.landcover != CROPLAND) {
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	}
	return cropphen;
}

void cropphen_struct::serialize(ArchiveStream& arch) {
	arch & sdate
		& sdate_harv
		& sdate_harvest
		& sdate_thisyear 
		& hdate
		& hdate_harvest
		& hlimitdate
		& hucountend
		& nharv
		& sendate
		& bicdate 
		& eicdate
		& tb
		& pvd
		& vdsum
		& vrf
		& prf
		& phu 
		& phu_old
		& husum_max
		& husum_sampled
		& husum_max_10
		& nyears_hu_sample
		& husum
		& fphu 
		& fphu_harv
		& demandsum_crop
		& supplysum_crop
		& growingseason 
		& growingseason_ystd
		& senescence
		& senescence_ystd
		& intercropseason;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patch member functions
////////////////////////////////////////////////////////////////////////////////


Patch::Patch(int i,Stand& s,Soiltype& st):
	id(i),stand(s),vegetation(*this),soil(*this,st),fluxes(*this) {

	for(unsigned int p = 0; p < pftlist.nobj; p++) {
		pft.createobj(pftlist[p]);
	}

	age = 0;
	disturbed = false;
	managed = false;
	wdemand = 0.0;
	wdemand_leafon = 0.0;
		
	growingseasondays = 0;

	fireprob = 0.0;
	ndemand = 0.0;
	dnfert = 0.0;
	anfert = 0.0;
	nharv = 0;
}

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
		& managed
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
		& ndemand
		& irrigation_y;
}

const Climate& Patch::get_climate() const {
	// All patches within a stand share the same climate
	return stand.get_climate();
}

bool Patch::has_fires() const {
#ifdef NOPASTURESTOCH
	return iffire && stand.landcover != CROPLAND && stand.landcover != PASTURE && !managed;
#else
	return iffire && stand.landcover != CROPLAND && !managed;
#endif
}

bool Patch::has_disturbances() const {
#ifdef NOPASTURESTOCH
	return ifdisturb && stand.landcover != CROPLAND && stand.landcover != PASTURE && !managed;
#else
	return ifdisturb && stand.landcover != CROPLAND && !managed;
#endif
}

/// C content of patch
/** 
 *  INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 */
double Patch::ccont(double scale_indiv, bool luc) {

	double ccont = 0.0;

	ccont += soil.cpool_fast;
	ccont += soil.cpool_slow;

	for(int i=0; i<NSOMPOOL-1; i++) {
		ccont += soil.sompool[i].cmass;
	}

	for(int i=0; i<npft; i++) {
		Patchpft& ppft = pft[i];
		ccont += ppft.litter_leaf;
		ccont += ppft.litter_root;
		ccont += ppft.litter_sap;
		ccont += ppft.litter_heart;
		ccont += ppft.harvested_products_slow;
	}

	for (unsigned int i=0; i<vegetation.nobj; i++) {

		Individual& indiv = vegetation[i];
		ccont += indiv.ccont(scale_indiv, luc);
	}

	return ccont;
}

/// N content of patch
/** 
 *  INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living N
 *  \param luc 					down-scales living N (used in N balance tests)
 */
double Patch::ncont(double scale_indiv, bool luc) {

	double ncont = 0.0;

	ncont += soil.nmass_avail;
	ncont += soil.snowpack_nmass;

	for(int i=0; i<NSOMPOOL-1; i++)
		ncont += soil.sompool[i].nmass;

	for(int i=0; i<npft; i++) {
		Patchpft& ppft = pft[i];
		ncont += ppft.nmass_litter_leaf;
		ncont += ppft.nmass_litter_root;
		ncont += ppft.nmass_litter_sap;
		ncont += ppft.nmass_litter_heart;
		ncont += ppft.harvested_products_slow_nmass;
	}

	for (unsigned int i=0; i<vegetation.nobj; i++) {

		Individual& indiv = vegetation[i];
		ncont += indiv.ncont(scale_indiv, luc);
	}

	return ncont;
}

/// C flux of patch
double Patch::cflux() {

	double cflux = 0.0;

	cflux += -fluxes.get_annual_flux(Fluxes::NPP);
	cflux += fluxes.get_annual_flux(Fluxes::REPRC);
	cflux += fluxes.get_annual_flux(Fluxes::SOILC);
	cflux += fluxes.get_annual_flux(Fluxes::FIREC);
	cflux += fluxes.get_annual_flux(Fluxes::ESTC);
	cflux += fluxes.get_annual_flux(Fluxes::SEEDC);
	cflux += fluxes.get_annual_flux(Fluxes::HARVESTC);

	return cflux;
}

/// N flux of patch
double Patch::nflux() {

	double nflux = 0.0;

	nflux += -stand.get_climate().andep;
	nflux += -anfert;
	nflux += -soil.anfix;
	nflux += soil.aminleach;
	nflux += soil.aorgleach;
	nflux += fluxes.get_annual_flux(Fluxes::HARVESTN);
	nflux += fluxes.get_annual_flux(Fluxes::SEEDN);
	nflux += fluxes.get_annual_flux(Fluxes::NH3_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::NO_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::NO2_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N2O_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N2_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N_SOIL);

	return nflux;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Standpft member functions
////////////////////////////////////////////////////////////////////////////////


void Standpft::serialize(ArchiveStream& arch) {
	arch & cmass_repr
		& anetps_ff_max
		& fpc_total
		& active;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell* gc, Soiltype& st, landcovertype landcoverX, int no_patch)
 : id(i),
   gridcell(gc),
   soiltype(st),
   landcover(landcoverX),
   origin(landcoverX),
   frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL = 1;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}

#if defined NOPASTURESTOCH
	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND || landcover== BARREN) {
#else
	if(landcover==CROPLAND || landcover==URBAN || landcover==PEATLAND || landcover== BARREN) {
#endif
		npatchL=1;
	}
	else {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	if(!(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND || landcover==NATURAL || landcover==FOREST || landcover== BARREN)) {
		// Someone has added a new landcover type, the code above needs to be updated and
		// npatchL needs to be set properly.
		fail("Unrecognized landcover type");
	}

	if(no_patch > 0)
		npatchL = no_patch;

	for (p=0;p<npatchL;p++) {
		createobj(*this, soiltype);
	}

	first_year = date.year;
	clone_year = -1;
	transfer_area_st = new double[nst];
	memset(transfer_area_st, 0, nst * sizeof(double));
	seed = 12345678;

	stid = -1;
	pftid = -1;
	current_rot = 0;
	ndays_inrotation = 0;
	infallow = false;
	isrotationday = false;
	isirrigated = false;
	hasgrassintercrop = false;
	gdd0_intercrop = 0.0;
	frac = 1.0;
	frac_old = 0.0;
	frac_temp = 0.0;
	protected_frac = 0.0;
	frac_change = 0.0;
	gross_frac_increase = 0.0;
	gross_frac_decrease = 0.0;
	cloned_fraction = 0.0;
	cloned = false;
	anpp = 0.0;
	cmass = 0.0;
	scale_LC_change = 1.0;
}

Stand::~Stand() {

	if(transfer_area_st)
		delete[] transfer_area_st;
}

double Stand::get_gridcell_fraction() const {
	return frac;
}

void Stand::rotate() {

	if(pftid >= 0 && stid >= 0) {

		ndays_inrotation = 0;

		current_rot = (current_rot + 1) % stlist[stid].rotation.ncrops;
		pftid = pftlist.getpftid(stlist[stid].management[current_rot].pftname);

		Standpft& standpft = pft[pftid];

		if(stlist[stid].management[current_rot].hydrology == IRRIGATED) {
			isirrigated = true;					
			standpft.irrigated = true;
		}
		else {
			isirrigated = false;					
			standpft.irrigated = false;
		}

		if(!readsowingdates)
			standpft.sdate_force = stlist[stid].management[current_rot].sdate;
		if(!readharvestdates)
			standpft.hdate_force = stlist[stid].management[current_rot].hdate;
	}
}

double Stand::transfer_area_lc(int to) {

	double area = 0.0;

	if(transfer_area_st) {

		for(int j=0; j<nst; j++) {

			if(stlist[j].landcover == to)		
				area += transfer_area_st[j];
		}
	}
	return area;
}

double Stand::ccont(double scale_indiv) {

	double ccont = 0.0;

	for(unsigned int p = 0; p < nobj; p++)
		ccont += (*this)[p].ccont(scale_indiv) / nobj;

	return ccont;
}

double Stand::ncont(double scale_indiv) {

	double ncont = 0.0;

	for(unsigned int p = 0; p < nobj; p++)
		ncont += (*this)[p].ncont(scale_indiv) / nobj;

	return ncont;
}

double Stand::cflux() {

	double cflux = 0.0;

	for(unsigned int p = 0; p < nobj; p++)
		cflux += (*this)[p].cflux() / nobj;

	return cflux;
}

double Stand::nflux() {

	double nflux = 0.0;

	for(unsigned int p = 0; p < nobj; p++)
		nflux += (*this)[p].nflux() / nobj;

	return nflux;
}

Stand& Stand::clone(StandType& st, double fraction) {

	// Serialize this stand to an in-memory stream
	std::stringstream ss;
	ArchiveOutStream aos(ss);
	serialize(aos);

	// Create a new stand in the gridcell...
	// NB: the patch number is always that of the old stand, even if the new stand is a pasture or crop stand
	Stand& new_stand = gridcell->create_stand(st.landcover, nobj);
	int new_seed = new_stand.seed;

	// ...and deserialize to that stand
	ArchiveInStream ais(ss);
	new_stand.serialize(ais);

	new_stand.clone_year = date.year;
//	new_stand.seed = new_seed;	// ?

	// Set land use settings for new stand
	new_stand.init_stand_lu(st, fraction);

//	for(unsigned int p = 0; p < nobj; p++)	// probably not what we want
//		new_stand[p].age = 0;

	return new_stand;
}

double Stand::get_landcover_fraction() const {
	if(get_gridcell().landcoverfrac[landcover])
		return frac / get_gridcell().landcoverfrac[landcover];
	else
		return 0.0;
}

void Stand::set_gridcell_fraction(double fraction) {
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
	else {	//NB. needs to be modified for use with land cover change (?)
		pft.killall();
		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			Standpft& standpft = pft.createobj(pftlist[i]);
			arch & standpft;
		}

		killall();
		unsigned int npatch;
		arch & npatch;
		for (unsigned int k = 0; k < npatch; k++) {
			Patch& patch = createobj(*this, soiltype);
			arch & patch;
		}
	}

	arch & first_year
		& frac
		& stid
		& pftid
		& current_rot
		& isirrigated
		& hasgrassintercrop
		& gdd0_intercrop
		& landcover
		& seed;
}

const Climate& Stand::get_climate() const {

	// In this implementation all stands within a grid cell
	// share the same climate. Note that this might not be
	// true in all versions of LPJ-GUESS, some may have
	// different climate per landcover type for instance.

	return get_gridcell().climate;
}

Gridcell& Stand::get_gridcell() const {
	assert(gridcell);
	return *gridcell;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp              = 0.0;
	fpc               = 0.0;
	fpc_daily		  = 0.0;
	densindiv         = 0.0;
	cmass_leaf        = 0.0;
	cmass_root        = 0.0;
	cmass_sap         = 0.0;
	cmass_heart       = 0.0;
	cmass_debt        = 0.0;
	cmass_leaf_post_turnover      = 0.0;
	cmass_root_post_turnover      = 0.0;
	cmass_tot_luc     = 0.0;
	phen              = 0.0;
	aphen             = 0.0;
	deltafpc          = 0.0;

	nmass_leaf        = 0.0;
	nmass_root        = 0.0;
	nmass_sap         = 0.0;
	nmass_heart       = 0.0;
	cton_leaf_aopt    = 0.0;
	cton_leaf_aavr    = 0.0;
	cton_status       = 0.0;
	cmass_veg         = 0.0;
	nmass_veg         = 0.0;
	nmass_tot_luc     = 0.0;

	nactive           = 0.0;
	nextin            = 1.0;
	nstore_longterm   = 0.0;
	nstore_labile     = 0.0;
	ndemand           = 0.0;
	fnuptake          = 1.0;
	anuptake          = 0.0;
	max_n_storage     = 0.0;
	scale_n_storage   = 0.0;

	leafndemand       = 0.0;
	rootndemand       = 0.0;
	sapndemand        = 0.0;
	storendemand      = 0.0;
	leaffndemand      = 0.0;
	rootfndemand      = 0.0;
	sapfndemand       = 0.0;
	storefndemand     = 0.0;
	leafndemand_store = 0.0;
	rootndemand_store = 0.0;

	nstress           = false;

	// additional initialisation
	age               = 0.0;
	fpar              = 0.0;
	aphen_raingreen   = 0;
	intercep          = 0.0;
	phen_mean         = 0.0;
	wstress           = false;
	lai               = 0.0;
	lai_layer         = 0.0;
	lai_indiv         = 0.0;
	lai_daily       = 0.0;
	lai_indiv_daily = 0.0;
	alive             = false;
	wscal_mean        = 1.0;

	int m;
	for (m=0; m<12; m++) {
		mlai[m] = 0.0;
		mlai_max[m] = 0.0;
	}

	// bvoc
	monstor           = 0.;
	iso               = 0.;
	mon               = 0.;
	fvocseas          = 1.;

	dnpp              = 0.0;
	cropindiv         = NULL;
	last_turnover_day = -1;

	Stand& stand = vegetation.patch.stand;

	if(pft.landcover==CROPLAND) {
		cropindiv=new cropindiv_struct;

		if (stand.pftid == pft.id) {
			cropindiv->isprimarycrop = true;
		}
		else if (stand.hasgrassintercrop && pft.isintercropgrass) {	// grass intercrop growth
			cropindiv->isintercropgrass = true;
		}
	}
//	dprintf("Year %d: Individual in stand %d created:id=%d, pft=%s\n", ::date.year-nyear_spinup+1901,vegetation.patch.stand.id,id,(char*)pft.name);
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
		& anpp
		& aet
		& aaet
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
		& nextin
		& nstore_longterm
		& nstore_labile
		& ndemand
		& fnuptake
		& anuptake
		& max_n_storage
		& scale_n_storage
		& avmaxnlim
		& cton_leaf_aopt
		& cton_leaf_aavr
		& cton_status
		& cmass_veg
		& nmass_veg

		& nstress
		& leafndemand
		& rootndemand
		& sapndemand
		& storendemand
		& leaffndemand
		& rootfndemand
		& sapfndemand
		& storefndemand
		& leafndemand_store
		& rootndemand_store
		& cmass_leaf_post_turnover
		& cmass_root_post_turnover
		& last_turnover_day
		& nday_leafon;

	if(pft.landcover==CROPLAND)
		arch & *cropindiv;
}

Individual::~Individual() {
	if(cropindiv)
		delete cropindiv;
//	dprintf("Year %d: Individual  in stand %d destroyed:id=%d, pft=%s, age=%.0f\n",::date.year-nyear_spinup+1901,vegetation.patch.stand.id,id,(char*)pft.name, age);
}

cropindiv_struct* Individual::get_cropindiv() {
	if (pft.landcover != CROPLAND) {
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	}
	return cropindiv;
}

cropindiv_struct* Individual::set_cropindiv() {
	if (pft.landcover != CROPLAND) {
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	}
	return cropindiv;
}

void cropindiv_struct::serialize(ArchiveStream& arch) {
	arch & grs_cmass_plant
		& grs_cmass_leaf
		& grs_cmass_root
		& grs_cmass_ho
		& grs_cmass_agpool 
		& isprimarycrop
		& isprimarycovegetation
		& isintercropgrass;
}



void Individual::report_flux(Fluxes::PerPFTFluxType flux_type, double value) {
	if (alive || istruecrop_or_intercropgrass()) {
		vegetation.patch.fluxes.report_flux(flux_type, pft.id, value);
	}
}

void Individual::report_flux(Fluxes::PerPatchFluxType flux_type, double value) {
	if (alive || istruecrop_or_intercropgrass()) {
		vegetation.patch.fluxes.report_flux(flux_type, value);
	}
}


/// Help function for reduce_biomass(), partitions nstore into leafs and roots
/** 
 *  As leaf and roots can have a very low N concentration after growth and allocation, 
 *  N in nstore() is split between them to saticfy relationship between their average C:N ratios
 */
void nstore_adjust(double& cmass_leaf,double& cmass_root, double& nmass_leaf, double& nmass_root, 
				   double nstore, double cton_leaf, double cton_root) {

	// (1) cmass_leaf / ((nmass_leaf + leaf_ndemand) * cton_leaf) = cmass_root / ((nmass_root + root_ndemand) * cton_root)
	// (2) leaf_ndemand + root_ndemand = nstore

	// (1) + (2) leaf_ndemand = (cmass_leaf * ratio (nmass_root + nstore) - cmass_root * nmass_leaf) / (cmass_root + cmass_leaf * ratio)
	//
	// where ratio = cton_root / cton_leaf

	double ratio = cton_root / cton_leaf;

	double leaf_ndemand = (cmass_leaf * ratio * (nmass_root + nstore) - cmass_root * nmass_leaf) / (cmass_root + cmass_leaf * ratio);
	double root_ndemand = nstore - leaf_ndemand;

	nmass_leaf += leaf_ndemand;
	nmass_root += root_ndemand;
}

void Individual::reduce_biomass(double mortality, double mortality_fire) {

	// This function needs to be modified if a new lifeform is added,
	// specifically to deal with nstore().
	assert(pft.lifeform == TREE || pft.lifeform == GRASS);

	if (!negligible(mortality)) {

		const double mortality_non_fire = mortality - mortality_fire;

		// Transfer killed biomass to litter
		// (above-ground biomass killed by fire enters atmosphere, not litter)

		Patchpft& ppft = patchpft();

		double cmass_leaf_litter = mortality * cmass_leaf;
		double cmass_root_litter = mortality * cmass_root;

		if(pft.landcover==CROPLAND) {
			if(pft.aboveground_ho)
				cmass_leaf_litter += mortality * cropindiv->cmass_ho;
			else
				cmass_root_litter += mortality * cropindiv->cmass_ho;

			cmass_leaf_litter += mortality * cropindiv->cmass_agpool;
		}

		ppft.litter_leaf += cmass_leaf_litter * mortality_non_fire / mortality;
		ppft.litter_root += cmass_root_litter;

		if (cmass_debt <= cmass_heart + cmass_sap) {
			if (cmass_debt <= cmass_heart) {
				ppft.litter_sap   += mortality_non_fire * cmass_sap;
				ppft.litter_heart += mortality_non_fire * (cmass_heart - cmass_debt);
			}
			else {
				ppft.litter_sap   += mortality_non_fire * (cmass_sap + cmass_heart - cmass_debt);
			}
		}
		else {
			double debt_excess = mortality_non_fire * (cmass_debt - (cmass_sap + cmass_heart));
			report_flux(Fluxes::NPP, debt_excess);
			report_flux(Fluxes::RA, -debt_excess);
		}

		double nmass_leaf_litter = mortality * nmass_leaf;
		double nmass_root_litter = mortality * nmass_root;

		if(pft.landcover==CROPLAND) {
			if(pft.aboveground_ho)
				nmass_leaf_litter += mortality * cropindiv->nmass_ho;
			else
				nmass_root_litter += mortality * cropindiv->nmass_ho;

			nmass_leaf_litter += mortality * cropindiv->nmass_agpool;
		}

		// stored N is partioned out to leaf and root biomass as new tissue after growth might have extremely low 
		// N content (to get closer to relationship between compartment averages (cton_leaf, cton_root, cton_sap))
		nstore_adjust(cmass_leaf_litter, cmass_root_litter, nmass_leaf_litter, nmass_root_litter,
			mortality * nstore(), pft.cton_leaf_avr,pft.cton_root_avr);

		ppft.nmass_litter_leaf  += nmass_leaf_litter * mortality_non_fire / mortality;
		ppft.nmass_litter_root  += nmass_root_litter;
		ppft.nmass_litter_sap   += mortality_non_fire * nmass_sap;
		ppft.nmass_litter_heart += mortality_non_fire * nmass_heart;

		// Flux to atmosphere from burnt above-ground biomass

		double cflux_fire = mortality_fire * (cmass_leaf_litter / mortality + cmass_wood());
		double nflux_fire = mortality_fire * (nmass_leaf_litter / mortality + nmass_wood());

		report_flux(Fluxes::FIREC,    cflux_fire);

		report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire); 
		report_flux(Fluxes::NO_FIRE,  Fluxes::NO_FIRERATIO  * nflux_fire);
		report_flux(Fluxes::NO2_FIRE, Fluxes::NO2_FIRERATIO * nflux_fire);
		report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);
		report_flux(Fluxes::N2_FIRE,  Fluxes::N2_FIRERATIO  * nflux_fire);

		// Reduce this Individual's biomass values

		const double remaining = 1.0 - mortality;

		if (pft.lifeform != GRASS) {
			densindiv *= remaining;
		}

		cmass_leaf      *= remaining;
		cmass_root      *= remaining;
		cmass_sap       *= remaining;
		cmass_heart     *= remaining;
		cmass_debt      *= remaining;
		if(pft.landcover==CROPLAND) {
			cropindiv->cmass_ho *= remaining;
			cropindiv->cmass_agpool *= remaining;
		}
		nmass_leaf      *= remaining;
		nmass_root      *= remaining;
		nmass_sap       *= remaining;
		nmass_heart     *= remaining;
		nstore_longterm *= remaining;
		nstore_labile   *= remaining;
		if(pft.landcover==CROPLAND) {
			cropindiv->nmass_ho *= remaining;
			cropindiv->nmass_agpool *= remaining;
		}
	}
}

double Individual::cton_leaf(bool use_phen /* = true*/) const {

	Stand& stand = vegetation.patch.stand;

	if (stand.ifnlim_stand()) {

		if(stand.is_true_crop_stand() && !negligible(cmass_leaf_today()) && !negligible(nmass_leaf)) {	//Detta kan möjligen tas bort
			return cmass_leaf_today() / nmass_leaf;
		}
		else if (!stand.is_true_crop_stand() && !negligible(cmass_leaf) && !negligible(nmass_leaf)) {
			if (use_phen) {
				if (!negligible(phen)) {
					return cmass_leaf_today() / nmass_leaf;
				}
				else {
					return pft.cton_leaf_avr;
				}
			}
			else {
				return cmass_leaf / nmass_leaf;
			}
		}
		else {
			return pft.cton_leaf_max;
		}
	}
	else {
		return pft.cton_leaf_avr;
	}
}

double Individual::cton_root(bool use_phen /* = true*/) const {

	Stand& stand = vegetation.patch.stand;

	if (stand.ifnlim_stand()) {
		if (!negligible(cmass_root) && !negligible(nmass_root)) { 
			if (use_phen) {
				if (!negligible(cmass_root_today())) {
					return cmass_root_today() / nmass_root;
				}
				else {
					return pft.cton_root_avr;
				}
			}
			else {
				return cmass_root / nmass_root;
			}
		}
		else {
			return pft.cton_root_max;
		}
	}
	else {
		return pft.cton_root_avr;
	}
}

double Individual::cton_sap() const {

	Stand& stand = vegetation.patch.stand;

	if (pft.lifeform == TREE) {
		if (stand.ifnlim_stand()) {
			if (!negligible(cmass_sap) && !negligible(nmass_sap))
				return cmass_sap / nmass_sap;
			else
				return pft.cton_sap_max;
		}
		else {
			return pft.cton_sap_avr;
		}
	}
	else {
		return 1.0;
	}
}

/// C content of individual
/** 
 *  INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 */
double Individual::ccont(double scale_indiv, bool luc) const {

	double ccont = 0.0;

	if (alive || istruecrop_or_intercropgrass()) {

		if(has_daily_turnover()) {	// Not taking into account future daily wood allocation/turnover

			if(cropindiv) {

				if(luc) {
					ccont += cropindiv->grs_cmass_leaf - cropindiv->grs_cmass_leaf_luc * (1.0 - scale_indiv);
					ccont += cropindiv->grs_cmass_root - cropindiv->grs_cmass_root_luc * (1.0 - scale_indiv);
				}
				else {
					ccont += cropindiv->grs_cmass_leaf * scale_indiv;
					ccont += cropindiv->grs_cmass_root * scale_indiv;
				}

				if(pft.phenology == CROPGREEN) {

					if(luc) {
						ccont += cropindiv->grs_cmass_ho - cropindiv->grs_cmass_ho_luc * (1.0 - scale_indiv);
						ccont += cropindiv->grs_cmass_agpool - cropindiv->grs_cmass_agpool_luc * (1.0 - scale_indiv);
						ccont += cropindiv->grs_cmass_dead_leaf - cropindiv->grs_cmass_dead_leaf_luc * (1.0 - scale_indiv);
						ccont += cropindiv->grs_cmass_stem - cropindiv->grs_cmass_stem_luc * (1.0 - scale_indiv);
					}
					else {
						ccont += cropindiv->grs_cmass_ho * scale_indiv;
						ccont += cropindiv->grs_cmass_agpool * scale_indiv;
						ccont += cropindiv->grs_cmass_dead_leaf * scale_indiv;
						ccont += cropindiv->grs_cmass_stem * scale_indiv;
					}
				}
			}
		}
		else {

			ccont += cmass_leaf * scale_indiv;
			ccont += cmass_root * scale_indiv;
			ccont += cmass_sap * scale_indiv;
			ccont += cmass_heart * scale_indiv;
			ccont -= cmass_debt * scale_indiv;

			if(pft.landcover == CROPLAND) {
				ccont += cropindiv->cmass_ho * scale_indiv;
				ccont += cropindiv->cmass_agpool * scale_indiv;
				// Yearly allocation not defined for crops with nlim
			}
		}
	}

	return ccont;
}

/// N content of individual
/** 
 *  INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living N
 *  \param luc 					down-scales living N (used in C balance tests)
 */
double Individual::ncont(double scale_indiv, bool luc) const {

	double ncont = 0.0;

	if(luc) {

		ncont += nmass_leaf - nmass_leaf_luc * (1.0 - scale_indiv);
		ncont += nmass_root - nmass_root_luc * (1.0 - scale_indiv);
		ncont += nmass_sap - nmass_sap_luc * (1.0 - scale_indiv);
		ncont += nmass_heart - nmass_heart_luc * (1.0 - scale_indiv);
		ncont += nstore_longterm - nstore_longterm_luc * (1.0 - scale_indiv);
		ncont += nstore_labile - nstore_labile_luc * (1.0 - scale_indiv);
	}
	else {
		ncont += nmass_leaf * scale_indiv;
		ncont += nmass_root * scale_indiv;
		ncont += nmass_sap * scale_indiv;
		ncont += nmass_heart * scale_indiv;
		ncont += nstore_longterm * scale_indiv;
		ncont += nstore_labile * scale_indiv;
	}

	if(pft.landcover == CROPLAND) {

		if(luc) {
			ncont += cropindiv->nmass_ho - cropindiv->nmass_ho_luc * (1.0 - scale_indiv);
			ncont += cropindiv->nmass_agpool - cropindiv->nmass_agpool_luc * (1.0 - scale_indiv);
			ncont += cropindiv->nmass_dead_leaf - cropindiv->nmass_dead_leaf_luc * (1.0 - scale_indiv);
		}
		else {
			ncont += cropindiv->nmass_ho * scale_indiv;
			ncont += cropindiv->nmass_agpool * scale_indiv;
			ncont += cropindiv->nmass_dead_leaf * scale_indiv;
		}
	}

	return ncont;
}

bool Individual::continous_grass() const {

	Stand& stand = vegetation.patch.stand;

	if(pft.landcover == CROPLAND) {

		StandType& st = stlist[stand.stid];
		bool sowing_restriction = true;

		for(int i=0; i<st.rotation.ncrops; i++) {
			int pftid = pftlist.getpftid(st.management[i].pftname);
			if(!stand.get_gridcell().pft[pftid].sowing_restriction)
				sowing_restriction = false;
		}

		if(cropindiv->isintercropgrass && sowing_restriction)
			return true;
		else
			return false;
	}
	else
		return false;
}

double Individual::ndemand_storage(double cton_leaf_opt) {

	if (vegetation.patch.stand.is_true_crop_stand() && ifnlim_lc[CROPLAND])	// only CROPGREEN, only ifnlim ?
		// analogous with root demand
		storendemand = max(0.0, cropindiv->grs_cmass_stem / (cton_leaf_opt * pft.cton_stem_avr / pft.cton_leaf_avr) - cropindiv->nmass_agpool);
	else
		storendemand = max(0.0, min(anpp * scale_n_storage / cton_leaf(), max_n_storage) - nstore());

	return storendemand;
}

/// Checks C mass and zeroes any negative value, balancing by adding to npp and reducing respiration
double Individual::check_C_mass() {

	if(pft.landcover != CROPLAND)
		return 0;

	double negative_cmass = 0.0;

	if(cropindiv->grs_cmass_leaf < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_leaf;
		cropindiv->ycmass_leaf -= cropindiv->grs_cmass_leaf;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_leaf;
		cropindiv->grs_cmass_leaf = 0.0;
	}
	if(cropindiv->grs_cmass_root < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_root;
		cropindiv->ycmass_root -= cropindiv->grs_cmass_root;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_root;
		cropindiv->grs_cmass_root = 0.0;
	}
	if(cropindiv->grs_cmass_ho < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_ho;
		cropindiv->ycmass_ho -= cropindiv->grs_cmass_ho;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_ho;
		cropindiv->grs_cmass_ho = 0.0;
	}
	if(cropindiv->grs_cmass_agpool < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_agpool;
		cropindiv->ycmass_agpool -= cropindiv->grs_cmass_agpool;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_agpool;
		cropindiv->grs_cmass_agpool = 0.0;
	}
	if(cropindiv->grs_cmass_dead_leaf < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->ycmass_dead_leaf -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->grs_cmass_dead_leaf = 0.0;
	}
	if(cropindiv->grs_cmass_stem < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_stem;
		cropindiv->ycmass_stem -= cropindiv->grs_cmass_stem;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_stem;
		cropindiv->grs_cmass_stem = 0.0;
	}

	if(negative_cmass > 1.0e-14) {
		anpp += negative_cmass;
		report_flux(Fluxes::NPP, negative_cmass);
		report_flux(Fluxes::RA, -negative_cmass);
//		dprintf("Year %d day %d Stand %d indiv %d: Negative C mass: %.15f\n", date.year, date.day, vegetation.patch.stand.id, id, negative_cmass);
	}

	return negative_cmass;
}

/// Checks C mass and zeroes any negative value, balancing by reducing C mass of other organs and (if needed) reducing anflux_landuse_change
double Individual::check_N_mass() {

	if(pft.landcover != CROPLAND && pft.landcover != PASTURE)
		return 0;

	double negative_nmass = 0.0;

	if(nmass_leaf < 0.0) {
		negative_nmass -= nmass_leaf;
		if(cropindiv)
			cropindiv->ynmass_leaf -= nmass_leaf;
		nmass_leaf = 0.0;
	}
	if(nmass_root < 0.0) {
		negative_nmass -= nmass_root;
		if(cropindiv)
			cropindiv->ynmass_root -= nmass_root;
		nmass_root = 0.0;
	}
	if(cropindiv) {
		if(cropindiv->nmass_ho < 0.0) {
			negative_nmass -= cropindiv->nmass_ho;
			cropindiv->ynmass_ho -= cropindiv->nmass_ho;
			cropindiv->nmass_ho = 0.0;
		}
		if(cropindiv->nmass_agpool < 0.0) {
			negative_nmass -= cropindiv->nmass_agpool;
			cropindiv->ynmass_agpool -= cropindiv->nmass_agpool;
			cropindiv->nmass_agpool = 0.0;
		}
		if(cropindiv->nmass_dead_leaf < 0.0) {
			negative_nmass -= cropindiv->nmass_dead_leaf;
			cropindiv->ynmass_dead_leaf -= cropindiv->nmass_dead_leaf;
			cropindiv->nmass_dead_leaf = 0.0;
		}
	}
	if(nstore_labile < 0.0) {
		negative_nmass -= nstore_labile;
		nstore_labile = 0,0;
	}	
	if(nstore_longterm < 0.0) {
		negative_nmass -= nstore_longterm;
		nstore_longterm = 0,0;
	}	

	if(negative_nmass > 1.0e-14) {
		double pos_nmass = ncont();
		if(pos_nmass > negative_nmass) {
			nmass_leaf -= negative_nmass * nmass_leaf / pos_nmass;
			nmass_root -= negative_nmass * nmass_root / pos_nmass;
			if(cropindiv) {
				cropindiv->nmass_ho -= negative_nmass * cropindiv->nmass_ho / pos_nmass;
				cropindiv->nmass_agpool -= negative_nmass * cropindiv->nmass_agpool / pos_nmass;
				cropindiv->nmass_dead_leaf -= negative_nmass * cropindiv->nmass_dead_leaf / pos_nmass;
			}
		}
		else {
			vegetation.patch.stand.get_gridcell().anflux_landuse_change -= (negative_nmass - pos_nmass) * vegetation.patch.stand.get_gridcell_fraction();
			nmass_leaf = 0,0;
			nmass_leaf = 0,0;
			if(cropindiv) {
				cropindiv->nmass_ho = 0,0;
				cropindiv->nmass_agpool = 0,0;
				cropindiv->nmass_dead_leaf = 0,0;
			}
		}
//		dprintf("Year %d day %d Stand %d indiv %d: Negative N mass: %.15f\n", date.year, date.day, vegetation.patch.stand.id, id, negative_nmass);
	}

	return negative_nmass;
}

bool Individual::is_turnover_day() const {

	if(patchpft().cropphen && patchpft().cropphen->growingseason) {

		const Climate& climate = vegetation.patch.get_climate();

		if(date.day == climate.testday_prec)
			return true;
		else
			return false;
	}
	else 
		return false;
}

Patchpft& Individual::patchpft() const {
	return vegetation.patch.pft[pft.id];
}

/// Save cmass-values on first day of the year of land cover change in expanding stands
void Individual::save_cmass_luc() {
	Stand& stand = vegetation.patch.stand;

	cmass_tot_luc = 0.0;

	if(cropindiv) {
		cropindiv->grs_cmass_leaf_luc = cropindiv->grs_cmass_leaf;
		cropindiv->grs_cmass_root_luc = cropindiv->grs_cmass_root;
		cropindiv->grs_cmass_ho_luc = cropindiv->grs_cmass_ho;
		cropindiv->grs_cmass_agpool_luc = cropindiv->grs_cmass_agpool;
		cropindiv->grs_cmass_dead_leaf_luc = cropindiv->grs_cmass_dead_leaf;
		cropindiv->grs_cmass_stem_luc = cropindiv->grs_cmass_stem;
	}
	cmass_tot_luc = ccont();
}

/// Save nmass-values on first day of the year of land cover change in expanding stands
void Individual::save_nmass_luc() {
	Stand& stand = vegetation.patch.stand;

	nmass_leaf_luc = nmass_leaf;
	nmass_root_luc = nmass_root;
	nmass_sap_luc = nmass_sap;
	nmass_heart_luc = nmass_heart;
	nstore_longterm_luc = nstore_longterm;
	nstore_labile_luc = nstore_labile;

	if(cropindiv) {
		cropindiv->nmass_ho_luc = cropindiv->nmass_ho;
		cropindiv->nmass_agpool_luc = cropindiv->nmass_agpool;
		cropindiv->nmass_dead_leaf_luc = cropindiv->nmass_dead_leaf;
	}
	nmass_tot_luc = ncont();
}

/// Gets the individual's daily cmass_leaf value
double Individual::cmass_leaf_today() const {

	if(istruecrop_or_intercropgrass()) {

		if(patchpft().cropphen->growingseason)
			return cropindiv->grs_cmass_leaf;
		else
			return 0.0;
	}
	else
		return cmass_leaf * phen;
}

/// Gets the individual's daily cmass_root value
double Individual::cmass_root_today() const {

	if(istruecrop_or_intercropgrass()) {

		if(patchpft().cropphen->growingseason)
			return cropindiv->grs_cmass_root;
		else
			return 0.0;
	}
	else
		return cmass_root * phen;
}

double Individual::fpc_today() const {

	if(pft.phenology == CROPGREEN) {

		if(patchpft().cropphen->growingseason)
			return fpc_daily;
		else
			return 0.0;
	}
	else
		return fpc * phen;
}

double Individual::lai_today() const {

	if(pft.phenology == CROPGREEN) {

		if(patchpft().cropphen->growingseason)
			return lai_daily;
		else
			return 0.0;
	}
	else
		return lai * phen;
}

double Individual::lai_indiv_today() const {

	if(pft.phenology == CROPGREEN) {

		if(patchpft().cropphen->growingseason)
			return lai_indiv_daily;
		else
			return 0.0;
	}
	else
		return lai_indiv * phen;
}

/// Get the Nitrigen limited LAI
double Individual::lai_nitrogen_today() const{
	if(pft.phenology==CROPGREEN) {

		double Ln = 0.0;
		if(patchpft().cropphen->growingseason){

			double k = 0.5;
			double ktn = 0.52*k+0.01; //Yin et al 2003
			double nb = 1/(pft.cton_leaf_max*pft.sla);
			if(cmass_leaf_today()>0.0){
				Ln = (1/ktn)*log(1+ktn*nmass_leaf/nb);
			}
		}
		return Ln;
	}
	else {
		return 1.0;
	}
}

bool Individual::growingseason() const {
	if(patchpft().cropphen)
		return patchpft().cropphen->growingseason;
	else
		return true;
}

bool Individual::has_daily_turnover() const {

#ifdef HARVEST_GRSC
		return istruecrop_or_intercropgrass();
#else
		return false;
#endif
}

/// Help function for kill(), partitions wood biomass into litter and harvest
/** 
 *  Wood biomass (either C or N) is partitioned into litter pools and
 *  harvest, according to PFT specific harvest fractions.
 *
 *  Biomass is sent in as sap and heart, any debt should already have been
 *  subtracted from these before calling this function.
 *
 *  \param mass_sap          Sapwood
 *  \param mass_heart        Heartwood
 *  \param harv_eff          Harvest efficiency (fraction of biomass harvested)
 *  \param harvest_slow_frac Fraction of harvested products that goes into slow depository
 *  \param res_outtake       Fraction of residue outtake at harvest
 *  \param litter_sap        Biomass going to sapwood litter pool
 *  \param litter_heart      Biomass going to heartwood litter pool
 *  \param fast_harvest      Biomass going to harvest flux
 *  \param slow_harvest      Biomass going to slow depository
 */
void partition_wood_biomass(double mass_sap, double mass_heart,
                            double harv_eff, double harvest_slow_frac, double res_outtake,
                            double& litter_sap, double& litter_heart,
                            double& fast_harvest, double& slow_harvest) {

	double sap_left = mass_sap;
	double heart_left = mass_heart;

	// Remove harvest
	double total_wood_harvest = harv_eff * (sap_left + heart_left);

	sap_left   *= 1 - harv_eff;
	heart_left *= 1 - harv_eff;

	// Partition wood harvest into slow and fast
	slow_harvest = total_wood_harvest * harvest_slow_frac;
	fast_harvest = total_wood_harvest * (1 - harvest_slow_frac);

	// Remove residue outtake
	fast_harvest += res_outtake * (sap_left + heart_left);
				
	sap_left   *= 1 - res_outtake;
	heart_left *= 1 - res_outtake;

	// The rest goes to litter
	litter_sap   = sap_left;
	litter_heart = heart_left;
}


void Individual::kill(bool harvest /* = false */) {
	Patchpft& ppft = patchpft();

	double charvest_flux = 0.0;
	double charvested_products_slow = 0.0;

	double nharvest_flux = 0.0;
	double nharvested_products_slow = 0.0;

	double harv_eff = 0.0;
	double harvest_slow_frac = 0.0;
	double res_outtake = 0.0;

	// The function always deals with harvest, but the harvest
	// fractions are zero when there is no harvest.
	if (harvest) {
		harv_eff = pft.harv_eff;

		if (ifslowharvestpool) {
			harvest_slow_frac = pft.harvest_slow_frac;
		}
		
		res_outtake = pft.res_outtake;
	}

	// C doesn't return to litter/harvest if the Individual isn't alive
	if(alive || istruecrop_or_intercropgrass()) {

		// For leaf and root, catches small, negative values too

		// Leaf: remove residue outtake and send the rest to litter
		if(has_daily_turnover() && cropindiv) {

			if(pft.lifeform == GRASS && pft.phenology != CROPGREEN) {
				charvest_flux += cropindiv->grs_cmass_leaf * harv_eff;
				cropindiv->grs_cmass_leaf *= (1 - harv_eff);
			}

			ppft.litter_leaf += cropindiv->grs_cmass_leaf * (1 - res_outtake);
			charvest_flux    += cropindiv->grs_cmass_leaf * res_outtake;
		}
		else {

			if(pft.lifeform == GRASS && pft.phenology != CROPGREEN) {
				charvest_flux += cmass_leaf * harv_eff;
				cmass_leaf *= (1 - harv_eff);
			}
			ppft.litter_leaf += cmass_leaf * (1 - res_outtake);
			charvest_flux    += cmass_leaf * res_outtake;
		}
		// Root: all goes to litter
		if(has_daily_turnover() && cropindiv)
			ppft.litter_root += cropindiv->grs_cmass_root;
		else
			ppft.litter_root += cmass_root;

		if(pft.landcover==CROPLAND) {

			if(has_daily_turnover()) {

				charvest_flux += cropindiv->grs_cmass_ho * harv_eff;
				cropindiv->grs_cmass_ho *= (1 - harv_eff);

				if(pft.aboveground_ho) {
					ppft.litter_leaf+=cropindiv->grs_cmass_ho * (1 - res_outtake);
					charvest_flux += cropindiv->grs_cmass_ho * res_outtake;
				}
				else {
					ppft.litter_root+=cropindiv->grs_cmass_ho;
				}
				ppft.litter_leaf+=cropindiv->grs_cmass_agpool * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_agpool * res_outtake;

				ppft.litter_leaf+=cropindiv->grs_cmass_dead_leaf * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_dead_leaf * res_outtake;

				ppft.litter_leaf+=cropindiv->grs_cmass_stem * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_stem * res_outtake;
			}
			else {

				charvest_flux += cropindiv->cmass_ho * harv_eff;
				cropindiv->cmass_ho *= (1 - harv_eff);

				if(pft.aboveground_ho) {
					ppft.litter_leaf+=cropindiv->cmass_ho * (1 - res_outtake);
					charvest_flux += cropindiv->cmass_ho * res_outtake;
				}
				else {
					ppft.litter_root+=cropindiv->cmass_ho;
				}
				ppft.litter_leaf+=cropindiv->cmass_agpool * (1 - res_outtake);
				charvest_flux += cropindiv->cmass_agpool * res_outtake;
			}
		}

		// Deal with the wood biomass and carbon debt for trees
		if (pft.lifeform == TREE) {

			// debt smaller than existing wood biomass
			if (cmass_debt <= cmass_sap + cmass_heart) {

				// before partitioning the biomass into litter and harvest,
				// first get rid of the debt so we're left with only
				// sap and heart
				double to_partition_sap   = 0.0;
				double to_partition_heart = 0.0;

				if (cmass_heart >= cmass_debt) {
					to_partition_sap   = cmass_sap;
					to_partition_heart = cmass_heart - cmass_debt;
				}
				else {
					to_partition_sap   = cmass_sap + cmass_heart - cmass_debt;
				}

				double clitter_sap, clitter_heart, cwood_harvest;

				partition_wood_biomass(to_partition_sap, to_partition_heart,
				                       harv_eff, harvest_slow_frac, res_outtake,
				                       clitter_sap, clitter_heart,
				                       cwood_harvest, charvested_products_slow);
				
				ppft.litter_sap   += clitter_sap;
				ppft.litter_heart += clitter_heart;

				charvest_flux += cwood_harvest;
			}
			// debt larger than existing wood biomass
			else {
				double debt_excess = cmass_debt - (cmass_sap + cmass_heart);
				report_flux(Fluxes::NPP, debt_excess);
				report_flux(Fluxes::RA, -debt_excess);
			}
		}
	}

	// Nitrogen always return to soil litter
	if (pft.lifeform == TREE) {

		double nlitter_sap, nlitter_heart, nwood_harvest;

		// Transfer nitrogen storage to sapwood nitrogen litter/harvest
		partition_wood_biomass(nmass_sap + nstore(), nmass_heart,
		                       harv_eff, harvest_slow_frac, res_outtake,
		                       nlitter_sap, nlitter_heart,
		                       nwood_harvest, nharvested_products_slow);

		ppft.nmass_litter_sap   += nlitter_sap;
		ppft.nmass_litter_heart += nlitter_heart;

		nharvest_flux += nwood_harvest;
	}
	else {
		// Transfer nitrogen storage to root nitrogen litter
		ppft.nmass_litter_root += nstore();
	}

	// Leaf: remove residue outtake and send the rest to litter
	ppft.nmass_litter_leaf += nmass_leaf * (1 - res_outtake);
	nharvest_flux          += nmass_leaf * res_outtake;

	// Root: all goes to litter
	ppft.nmass_litter_root += nmass_root;

	if(pft.landcover==CROPLAND) {
		if(pft.aboveground_ho) {
			ppft.nmass_litter_leaf+=cropindiv->nmass_ho * (1 - res_outtake);
			nharvest_flux += cropindiv->nmass_ho * res_outtake;
		}
		else
			ppft.litter_root+=cropindiv->nmass_ho;
		ppft.nmass_litter_leaf+=cropindiv->nmass_agpool * (1 - res_outtake);
		nharvest_flux += cropindiv->nmass_agpool * res_outtake;
		ppft.nmass_litter_leaf += cropindiv->nmass_dead_leaf * (1 - res_outtake);
		nharvest_flux          += cropindiv->nmass_dead_leaf * res_outtake;
	}

	// Report harvest fluxes
	report_flux(Fluxes::HARVESTC, charvest_flux);
	report_flux(Fluxes::HARVESTN, nharvest_flux);

	// Add to biomass depositories for long-lived products
	ppft.harvested_products_slow += charvested_products_slow;
	ppft.harvested_products_slow_nmass += nharvested_products_slow;
}

/// Should be used together with check_patch() e.g. in framework()
void MassBalance::init_indiv(Individual& indiv) {

	Patch& patch = indiv.vegetation.patch;
	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return;
	Gridcell& gridcell = stand.get_gridcell();

	double scale = 1.0;
	if(patch.stand.get_gridcell().LC_updated && (patch.nharv == 0 || date.day == 0))
		scale = stand.scale_LC_change;

	ccont_zero = indiv.ccont();
	ccont_zero_scaled = indiv.ccont(scale, true);
	// Add soil C
	ccont_zero += patch.ccont(0.0);
	ccont_zero_scaled += patch.ccont(0.0);
	cflux_zero = patch.cflux();

	ncont_zero = indiv.ncont();
	ncont_zero_scaled = indiv.ncont(scale, true);
	// Add soil N
	ncont_zero += patch.ncont(0.0);
	ncont_zero_scaled += patch.ncont(0.0);
	nflux_zero = patch.nflux();
}

bool MassBalance::check_indiv_C(Individual& indiv, bool check_harvest) {

	bool balance = true;
	Patch& patch = indiv.vegetation.patch;
	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return balance;
	Gridcell& gridcell = stand.get_gridcell();
	double ccont = indiv.ccont();
	ccont += patch.ccont(0.0);
	double cflux = patch.cflux();

	if(check_harvest && patch.isharvestday)
		ccont_zero = ccont_zero_scaled;
if(date.year >= nyear_spinup)
	if(fabs(ccont - ccont_zero + cflux - cflux_zero) > 1.0e-10) {
		dprintf("\nStand %d Patch %d Indiv %d C balance year %d day %d: %.10f\n", patch.stand.id, patch.id, indiv.id, date.year, date.day, ccont - ccont_zero + cflux - cflux_zero);
		dprintf("C pool change: %.10f\n", ccont - ccont_zero);
		dprintf("C flux: %.10f\n\n",  cflux - cflux_zero);
		balance = false;
	}

	return balance;
}

bool MassBalance::check_indiv_N(Individual& indiv, bool check_harvest) {

	bool balance = true;
	Patch& patch = indiv.vegetation.patch;
	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return balance;
	Gridcell& gridcell = stand.get_gridcell();
	double ncont = indiv.ncont();
	ncont += patch.ncont(0.0);
	double nflux = patch.nflux();

	if(check_harvest && patch.isharvestday)
		ncont_zero = ncont_zero_scaled;
if(date.year >= nyear_spinup)
	if(fabs(ncont - ncont_zero + nflux - nflux_zero) > 1.0e-14) {
		dprintf("\nStand %d Patch %d Indiv %d N balance year %d day %d: %.10f\n", patch.stand.id, patch.id, indiv.id, date.year, date.day, ncont - ncont_zero + nflux - nflux_zero);
		dprintf("N pool change: %.14f\n", ncont - ncont_zero);
		dprintf("N flux: %.14f\n\n",  nflux - nflux_zero);
		balance = false;
	}

	return balance;
}

/// Should be preceded by init_patch() e.g. i framework()
/** check_harvest must be true if growth_daily() is tested
 *  canopy_exchange() and growth_daily() and functions in between cannot be tested separately
 */
bool MassBalance::check_indiv(Individual& indiv, bool check_harvest) {

	bool balance = true;

	balance = check_indiv_C(indiv, check_harvest);
	balance = balance && check_indiv_N(indiv, check_harvest);

	return balance;
}

/// Should be used together with check_patch() e.g. in framework()
void MassBalance::init_patch(Patch& patch) {

	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return;
	Gridcell& gridcell = stand.get_gridcell();

	double scale = 1.0;
	if(patch.stand.get_gridcell().LC_updated && (patch.nharv == 0 || date.day == 0))
		scale = stand.scale_LC_change;

	ccont_zero = patch.ccont();
	ccont_zero_scaled = patch.ccont(scale, true);
	cflux_zero = patch.cflux();

	if(stand.get_gridcell_fraction())
		cflux_zero += gridcell.acflux_harvest_slow / stand.get_gridcell_fraction();

	ncont_zero = patch.ncont();
	ncont_zero_scaled = patch.ncont(scale, true);
	nflux_zero = patch.nflux();

	if(stand.get_gridcell_fraction())
		nflux_zero += gridcell.anflux_harvest_slow / stand.get_gridcell_fraction();
}

bool MassBalance::check_patch_C(Patch& patch, bool check_harvest) {

	bool balance = true;
	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return balance;
	Gridcell& gridcell = stand.get_gridcell();
	double ccont = patch.ccont();
	double cflux = patch.cflux();

	if(stand.get_gridcell_fraction())
		cflux += gridcell.acflux_harvest_slow / stand.get_gridcell_fraction();

	if(check_harvest && patch.isharvestday)
		ccont_zero = ccont_zero_scaled;
if(date.year >= nyear_spinup)
	if(fabs(ccont - ccont_zero + cflux - cflux_zero) > 1.0e-10) {
		dprintf("\nStand %d Patch %d C balance year %d day %d: %.10f\n", patch.stand.id, patch.id, date.year, date.day, ccont - ccont_zero + cflux - cflux_zero);
		dprintf("C pool change: %.10f\n", ccont - ccont_zero);
		dprintf("C flux: %.10f\n\n",  cflux - cflux_zero);
		balance = false;
	}

	return balance;
}

bool MassBalance::check_patch_N(Patch& patch, bool check_harvest) {

	bool balance = true;
	Stand& stand = patch.stand;
	if(!stand.is_true_crop_stand())
		return balance;
	Gridcell& gridcell = stand.get_gridcell();
	double ncont = patch.ncont();
	double nflux = patch.nflux();

	if(stand.get_gridcell_fraction())
		nflux += gridcell.anflux_harvest_slow / stand.get_gridcell_fraction();

	if(check_harvest && patch.isharvestday)
		ncont_zero = ncont_zero_scaled;
if(date.year >= nyear_spinup)
	if(fabs(ncont - ncont_zero + nflux - nflux_zero) > 1.0e-14) {
		dprintf("\nStand %d Patch %d N balance year %d day %d: %.14f\n", patch.stand.id, patch.id, date.year, date.day, ncont - ncont_zero + nflux - nflux_zero);
		dprintf("N pool change: %.14f\n", ncont - ncont_zero);
		dprintf("N flux: %.14f\n\n",  nflux - nflux_zero);
		balance = false;
	}

	return balance;
}

/// Should be preceded by init_patch() e.g. i framework()
/** check_harvest must be true if growth_daily() is tested
 *  canopy_exchange() and growth_daily() and functions in between cannot be tested separately
 */
bool MassBalance::check_patch(Patch& patch, bool check_harvest) {

	bool balance = true;

	balance = check_patch_C(patch, check_harvest);
	balance = balance && check_patch_N(patch, check_harvest);

	return balance;
}

void MassBalance::check_year(Gridcell& gridcell) {

	if(date.year >= start_year) {

		double ccont_year = gridcell.ccont();
		double cflux_year = gridcell.cflux();

		double ncont_year = gridcell.ncont();
		double nflux_year = gridcell.nflux();

		if(date.year >= start_year) {

			if(date.year == start_year) {
				ccont_zero = ccont_year;
				ncont_zero = ncont_year;
			}
			else if(date.year > start_year) {

				cflux += cflux_year;
				nflux += nflux_year;

				// C balance check:
				if(fabs(ccont_year - ccont + cflux_year) > 1.0e-9) {
					dprintf("\nC balance year %d: %.10f\n", date.year, ccont_year - ccont + cflux_year);
					dprintf("C pool change: %.5f\n", ccont_year - ccont);
					dprintf("C flux: %.5f\n",  cflux_year);
				}
if(ifnlim_lc[CROPLAND]) {
				// N balance check:
				if(fabs(ncont_year - ncont + nflux_year) > 1.0e-3) {
					dprintf("\nN balance year %d: %.4f\n", date.year, ncont_year - ncont + nflux_year);
					dprintf("N pool change: %.4f\n", ncont_year - ncont);
					dprintf("N flux: %.4f\n",  nflux_year);
				}
}
			}
			ccont = ccont_year;
			ncont = ncont_year;
		}
	}
}

void MassBalance::check_period() {

	// C balance check:
	if(fabs(ccont - ccont_zero + cflux) > 1.0e-9) {
		dprintf("\nWARNING: Period C balance: %.10f\n", ccont - ccont_zero + cflux);
		dprintf("C pool change: %.5f\n", ccont - ccont_zero);
		dprintf("C fluxes: %.5f\n",  cflux);
	}
if(ifnlim_lc[CROPLAND]) {
	// N balance check:
	if(fabs(ncont - ncont_zero + nflux) > 1.0e-3) {
		dprintf("\nWARNING: Period N balance: %.4f\n", ncont - ncont_zero + nflux);
		dprintf("N pool change: %.4f\n", ncont - ncont_zero);
		dprintf("N fluxes: %.4f\n",  nflux);
	}
}
}

void MassBalance::init(Gridcell& gridcell) {

//	start_year = date.year;
	ccont_zero = gridcell.ccont();
	cflux_zero = gridcell.cflux();
}

void MassBalance::check(Gridcell& gridcell) {

	double ccont = gridcell.ccont();
	double cflux = gridcell.cflux();

	if(fabs(ccont - ccont_zero + cflux) > 1.0e-5) {
		dprintf("\nC balance year %d: %.5f\n", date.year, ccont - ccont_zero + cflux);
		dprintf("C pool change: %.5f\n", ccont - ccont_zero);
		dprintf("C flux: %.5f\n\n",  cflux);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellpft member functions
////////////////////////////////////////////////////////////////////////////////


void Gridcellpft::serialize(ArchiveStream& arch) {
	arch & addtw
		& Km
		& autumnoccurred
		& springoccurred
		& vernstartoccurred
		& vernendoccurred 
		& precoccurred
		& first_autumndate
		& first_autumndate20
		& first_autumndate_20
		& last_springdate
		& last_springdate20
		& last_springdate_20 
		& last_verndate
		& last_verndate20
		& last_verndate_20
		& first_precdate
		& sdate_default
		& sdatecalc_temp
		& sdatecalc_prec 
		& sdate_force
		& hdate_force
		& hlimitdate_default
		& wintertype
		& multicrop
		& swindow	// ?
		& sowing_restriction;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcell member functions
////////////////////////////////////////////////////////////////////////////////

Gridcell::Gridcell():climate(*this) {
	landcovertype landcover;
	LC_updated = false;

	for(unsigned int p=0; p<pftlist.nobj; p++) {
		pft.createobj(pftlist[p]);
	}

	memset(landcoverfrac, 0, sizeof(double) * NLANDCOVERTYPES);
	memset(landcoverfrac_old, 0, sizeof(double) * NLANDCOVERTYPES);
	acflux_harvest_slow=0.0;
	acflux_landuse_change=0.0;
	anflux_harvest_slow=0.0;
	anflux_landuse_change=0.0;
	memset(acflux_harvest_slow_lc, 0, sizeof(double)*NLANDCOVERTYPES);
	memset(acflux_landuse_change_lc, 0, sizeof(double)*NLANDCOVERTYPES);

	for(int i=0; i<NLANDCOVERTYPES; i++) {		
		if(i == NATURAL || i == FOREST)
			expand_to_new_stand[i] = true;
		else
			expand_to_new_stand[i] = false;

		pool_to_all_landcovers[i] = false;		// from a donor landcover; alt.c
		pool_from_all_landcovers[i] = false;		// to a receptor landcover; alt.a

/*		if(i == CROPLAND) {
			pool_to_all_landcovers[i] = true;
			pool_to_all_standtypes[i] = true;
		}
		else {
			pool_to_all_landcovers[i] = false;
			pool_to_all_standtypes[i] = false;
		}
*/
	}

	if(!run_landcover) {
		landcover = NATURAL;
		create_stand(landcover);
		landcoverfrac[NATURAL] = 1.0;
	}

	seed = 12345678;
}

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

Stand& Gridcell::create_stand_lu(StandType& st, double fraction, int no_patch) {

	landcovertype lc = st.landcover;

	Stand& stand = create_stand(lc, no_patch);
	stand.init_stand_lu(st, fraction);

	return stand;
}

void Stand::init_stand_lu(StandType& st, double fraction) {

	int error = 0;
	landcovertype lc = st.landcover;
	landcover = lc;

	stid = st.id;
	set_gridcell_fraction(fraction);
	frac_old = 0.0;
	frac_change = fraction;
	gross_frac_increase = fraction;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pftx = pftlist.getobj();

		if(!st.restrictpfts && pftx.landcover == lc
			|| st.naturalveg && pftx.landcover == NATURAL
			|| st.naturalgrass && pftx.landcover == NATURAL && pftx.lifeform == GRASS) {

			pft[pftx.id].active = true;
		}
		else {
			pft[pftx.id].active = false;
		}
		pftlist.nextobj();
	}

	if(lc == CROPLAND) {

		pftid = pftlist.getpftid(st.management[0].pftname);	// First main crop, will change during crop rotation
		current_rot = 0;

#ifdef IRRIGATION
		if(st.management[0].hydrology == IRRIGATED) {
			isirrigated = true;								// First main crop, may change during crop rotation
			if(pftid >= 0)
				pft[pftid].irrigated = true;
		}
#endif
		if(st.intercrop==NATURALGRASS && ifintercropgrass) {
			hasgrassintercrop = true;

			for(unsigned int i=0; i<pftlist.nobj; i++) {
				if(pftlist[i].isintercropgrass)
					pft[pftlist[i].id].active = true;
			}
		}

		// Set standpft- and patchpft-variables for all active crops in all rotations
		for(int rot=0; rot<st.rotation.ncrops; rot++) {

			int id = pftlist.getpftid(st.management[rot].pftname);

			if(id >=0) {
				pft[id].active = true;

				if(rot == 0) {
					// Set crop cycle dates to default values only for first crop in a rotation.
					for(unsigned int p = 0; p < nobj; p++) {

						Gridcellpft& gcpft = get_gridcell().pft[id]; 
						Patchpft& ppft = (*this)[p].pft[id];

						ppft.set_cropphen()->sdate = gcpft.sdate_default;
						ppft.set_cropphen()->hlimitdate = gcpft.hlimitdate_default;
				
						if(pftlist[id].phenology == ANY)
							ppft.set_cropphen()->growingseason = true;
						else if(pftlist[id].phenology == CROPGREEN) {
							ppft.set_cropphen()->eicdate = stepfromdate(ppft.get_cropphen()->sdate, -15);
						}
					}
				}
			}
			else {
				dprintf("Warning: stand type %d pft %s not in pftlist !\n", stid, (char*)st.management[rot].pftname);;
				break;
			}
		}
	}
}

double Gridcell::ccont() {

	double ccont = 0.0;

	for(unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ccont += stand.ccont() * stand.get_gridcell_fraction();
	}

	return ccont;
}

double Gridcell::ncont() {

	double ncont = 0.0;

	for(unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ncont += stand.ncont() * stand.get_gridcell_fraction();
	}

	return ncont;
}

double Gridcell::cflux() {

	double cflux = 0.0;

	for(unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		cflux += stand.cflux() * stand.get_gridcell_fraction();
	}

	cflux += acflux_landuse_change;
	cflux += acflux_harvest_slow;

	return cflux;
}

double Gridcell::nflux() {

	double nflux = 0.0;

	for(unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		nflux += stand.nflux() * stand.get_gridcell_fraction();
	}

	nflux += anflux_landuse_change;
	nflux += anflux_harvest_slow;

	return nflux;
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

		unsigned int nstands = nbr_stands();
		arch & nstands;
		for (unsigned int s = 0; s < nstands; s++) {
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

		clear();
		unsigned int number_of_stands;
		arch & number_of_stands;
				
		for (unsigned int s = 0; s < number_of_stands; s++) {
			landcovertype landcover;
			arch & landcover;
			create_stand(landcover);
			arch & (*this)[s];
		}
	}
}

Stand& Gridcell::create_stand(landcovertype landcover, int no_patch) {
	Stand* stand = new Stand(get_next_id(), this, soiltype, landcover, no_patch);

	push_back(stand);

	return *stand;
}

Gridcell::iterator Gridcell::delete_stand(iterator itr) {
	return erase(itr);
}

unsigned int Gridcell::nbr_stands() const {
	return size();
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

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling 
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine, 
//   XXXV-XLIII, MIT Press, Mass.