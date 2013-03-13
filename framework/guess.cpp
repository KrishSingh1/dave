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
bool cftfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;				// If a slow harvested product pool is included in patchpft.
bool ifintercropgrass;
int ncft=0; // number of CFTs in Pftlist, set in plib_callback()
int nyear_spinup;		

xtring state_path;
bool restart;
bool save_state;
int state_year;

bool forcesowingdates;
bool forceharvestdates;
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
		& lon
		& insol
		& instype
		& eet
		& mtemp
		& mtemp_min20
		& mtemp_max20
		& mtemp_max
		& gdd5
		& gdd5_pasture
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
		& doneday
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
		& aprec
		& SOAsia;
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
		& phot_wstress
		& swindow
		& water_deficit_y;
	if(pft.landcover==CROPLAND)
		arch & *cropphen;
		
}

void cropphen_struct::serialize(ArchiveStream& arch) {
	arch & lai_crop_actual
		& sdate
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
		& husum_max_postharv
		& husum_max_hlim
		& husum_max_10
		& husum_h
		& husum
		& fphu 
		& fphu_harv
		& demandsum_crop
		& supplysum_crop
		& lai
		& fpc
		& growingseason 
		& growingseason_ystd
		& senescence
		& senescence_ystd
		& intercropseason
		& maincrop;
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
		& mpet
		& irrigation_y;
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

//const cropphen_struct* Patchpft::get_cropphen() 
cropphen_struct* Patchpft::get_cropphen() 
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropphen;
}

cropphen_struct* Patchpft::set_cropphen()
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropphen;
}


Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}

#if defined NOPASTURESTOCH
	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND) {
#else
	if(landcover==CROPLAND || landcover==URBAN || landcover==PEATLAND) {
#endif
		npatchL=1;
	}
	else {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}

	for (p=0;p<npatchL;p++) {
		createobj(*this,gc.soiltype);
	}

	first_year=date.year;
	natural_frac_change=0.0;
	seed=12345678;

	pftid=-1;
	cftid=-1;
	isirrigated=false;
	hasgrassintercrop=false;
	gdd0_intercrop=0.0;
}

double Stand::get_gridcell_fraction() const {
	return frac;
}

double Stand::get_landcover_fraction() const {
	if(gridcell.landcoverfrac[landcover])
		return frac/gridcell.landcoverfrac[landcover];
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
		& pftid
		& cftid
		& isirrigated
		& hasgrassintercrop
		& gdd0_intercrop
		& seed;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp=0.0;
	fpc=0.0;
	fpc_thisday=0.0;
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
	demand=0.0;
	supply=0.0;
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
	wscal_mean=1.0;

	int m;
	for (m=0;m<12;m++) {
		mlai[m]=mlai[m]=0.0;
	}

	// bvoc
	monstor=0.;
	iso=0.;
	mon=0.;
	fvocseas=1.;
	dtr_wstress=0.;
	eet_wstress=0.;
	agdd5_wstress=0.;
	rad_wstress=0.;		

	dnpp=0.0;
	cropindiv=NULL;

	if(pft.landcover==CROPLAND)
	{
		cropindiv=new cropindiv_struct;
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
		& monstor 
		& fvocseas 
		& dtr_wstress 
		& eet_wstress 
		& agdd5_wstress 
		& rad_wstress;
	if(pft.landcover==CROPLAND)
		arch & *cropindiv;
}

Individual::~Individual()
{
	if(cropindiv)
		delete cropindiv;

//	dprintf("Year %d: Individual  in stand %d destroyed:id=%d, pft=%s\n",::date.year-nyear_spinup+1901,vegetation.patch.stand.id,id,(char*)pft.name);
}

//const cropindiv_struct* Individual::get_cropindiv() 
cropindiv_struct* Individual::get_cropindiv() 
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
		return cropindiv;
}

cropindiv_struct* Individual::set_cropindiv()
{
	if(pft.landcover!=CROPLAND)
		fail("Only crop individuals have cropindiv struct. Re-write code !\n");
	else
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
	if (alive || pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || cropindiv->isintercropgrass)) {
		vegetation.patch.fluxes.report_flux(flux_type, pft.id, value);
	}
}

void Individual::report_flux(Fluxes::PerPatchFluxType flux_type, double value) {
	if (alive || pft.landcover==CROPLAND && (pft.phenology==CROPGREEN || cropindiv->isintercropgrass)) {
		vegetation.patch.fluxes.report_flux(flux_type, value);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellpft member functions
////////////////////////////////////////////////////////////////////////////////


void Gridcellpft::serialize(ArchiveStream& arch) {
	arch & addtw
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
		& singlecrop
		& swindow;
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
		& cftfrac
		& cftfrac_old
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
