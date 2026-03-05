///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.cpp
/// \brief LPJ-GUESS Combined Modular Framework
///
/// \author Ben Smith
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
int nst;  // number of possible stand types
int nst_lc[NLANDCOVERTYPES];  // number of possible stand types in each land cover type
int nmt;  // number of possible management types

ManagementTypelist mtlist;
StandTypelist stlist;
Pftlist pftlist;

// emission ratios from fire (NH3, NOx, N2O, N2) Delmas et al. 1995

const double Fluxes::NH3_FIRERATIO = 0.005;
const double Fluxes::NOx_FIRERATIO = 0.237;
const double Fluxes::N2O_FIRERATIO = 0.036;
const double Fluxes::N2_FIRERATIO  = 0.722;


////////////////////////////////////////////////////////////////////////////////
// Implementation of PhotosynthesisResult member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisResult::serialize(ArchiveStream& arch) {
	arch & agd_g
		& adtmm
		& rd_g
		& vm
		& je
		& vmaxnlim
		& vmax_5;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Climate member functions
////////////////////////////////////////////////////////////////////////////////

void Climate::serialize(ArchiveStream& arch) {
	arch & temp
		& rad
		& par
		& prec
		& aprec
		& aprec_lastyear
		& daylength
		& co2
		& lat
		& insol
		& instype	
		& ainsol
		& eet
		& mtemp
		& mtemp_min20
		& mtemp_max20
		& mtemp_max
		& gdd5
		& gdd0
		& agdd5
		& agdd0
		& agdd0_20
		& chilldays
		& ifsensechill
		& gtemp
		& dtemp_31
		& dprec_31
		& deet_31
		& mtemp_min_20
		& mtemp_max_20
		& mtemp_min
		& atemp_mean
		& sinelat
		& cosinelat
		& qo & u & v & hh & sinehh
		& daylength_save
		& doneday
		& maxtemp
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
		& hmtemp_20
		& hmprec_20
		& hmeet_20
		& seasonality
		& seasonality_lastyear
		& prec_seasonality
		& prec_seasonality_lastyear
		& prec_range
		& prec_range_lastyear
		& temp_seasonality
		& temp_seasonality_lastyear
		& var_prec
		& var_temp
		& rainfall_annual_avg
		& last_rainfall
		& days_since_last_rainfall
		& kbdi
		& ffdi_monthly
		& weathergenstate;
}

void WeatherGenState::serialize(ArchiveStream& arch) {

	arch & carry
		& xcng
		& xs
		& indx
		& have
		& gamma_vals
		& pday
		& resid
		& q;

}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Fluxes member functions
////////////////////////////////////////////////////////////////////////////////

Fluxes::Fluxes(Patch& p)
  : patch(p),
    annual_fluxes_per_pft(npft, std::vector<double>(NPERPFTFLUXTYPES)) {

	// Resize daily_fluxes_per_pft to required sizes
	// Order  daily_fluxes_per_pft[day][pft][type]
	daily_fluxes_per_pft.resize(date.year_length());
	for (int day = 0; day < date.year_length(); ++day) {
		daily_fluxes_per_pft[day].resize(npft);
		for (int pft = 0; pft < npft; ++pft) {
			daily_fluxes_per_pft[day][pft].resize(NPERPFTFLUXTYPES, 0.0);
		}
	}

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

	for (int d = 0; d < date.year_length(); ++d) {
		std::fill_n(daily_fluxes_pft[d], int(NPERPFTFLUXTYPES), 0);
		std::fill_n(daily_fluxes_patch[d], int(NPERPATCHFLUXTYPES), 0);
		for (int pft = 0; pft < npft; ++pft) {
			std::fill_n(daily_fluxes_per_pft[d][pft].begin(), int(NPERPFTFLUXTYPES), 0);
		}
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
	daily_fluxes_per_pft[date.day][pft_id][flux_type] += value;
}

void Fluxes::report_flux(PerPatchFluxType flux_type, double value) {
	monthly_fluxes_patch[date.month][flux_type] += value;
	daily_fluxes_patch[date.day][flux_type] += value;
}

double Fluxes::get_daily_flux(PerPFTFluxType flux_type, int day, int pft_id) const {
	return daily_fluxes_per_pft[day][pft_id][flux_type];
}

double Fluxes::get_daily_flux(PerPFTFluxType flux_type, int day) const {
	return daily_fluxes_pft[day][flux_type];
}

double Fluxes::get_daily_flux(PerPatchFluxType flux_type, int day) const {
	return daily_fluxes_patch[day][flux_type];
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
// Implementation of Pools member functions
////////////////////////////////////////////////////////////////////////////////

/*Pools::Pools(Soillayer& sl)
	: soillayer(sl) {

	reset();
}

void Pools::reset() {
	for (size_t i = 0; i < NPERPATCHPOOLTYPES; ++i) {

		annual_pools[i] = 0.0;

		for (int m = 0; m < 12; ++m) {
			monthly_pools[m][i] = 0.0;
		}

		for (int d = 0; d < date.year_length(); ++d) {
			daily_pools[d][i] = 0.0;
		}
	}
}

void Pools::serialize(ArchiveStream& arch) {
	arch & annual_pools
		& monthly_pools
		& daily_pools;
}

void Pools::report_flux(PerPatchPoolType pool_type, double value) {
	annual_pools[pool_type] += value;
	monthly_pools[date.month][pool_type] += value;
	daily_pools[date.day][pool_type] += value;
}

double Pools::get_daily_pool(PerPatchPoolType pool_type, int month) const {
	return monthly_pools[month][pool_type];
}

double Pools::get_monthly_pool(PerPatchPoolType pool_type, int month) const {
	return monthly_pools[month][pool_type];
}

double Pools::get_daily_pool(PerPatchPoolType pool_type, int day) const {
	return daily_pools[day][pool_type];
}*/

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
// Implementation of cropphen_struct member functions
////////////////////////////////////////////////////////////////////////////////

void cropphen_struct::serialize(ArchiveStream& arch) {
	arch & sdate
		& sdate_harv
		& sdate_harvest
		& hdate
		& hlimitdate
		& hucountend
		& sendate
		& bicdate
		& eicdate
		& growingdays
		& growingdays_y
		& lgp
		& tb
		& pvd
		& vdsum
		& vrf
		& prf
		& phu
		& phu_old
		& husum
		& husum_max
		& husum_sampled
		& husum_max_10
		& nyears_hu_sample
		& hu_samplingperiod
		& hu_samplingdays
		& fphu
		& fphu_harv
		& hi
		& fhi_harv
		& demandsum_crop
		& supplysum_crop
		& growingseason
		& growingseason_ystd
		& senescence
		& senescence_ystd
		& intercropseason
		& fertilised
		& vdsum_alloc
		& vd
		& dev_stage;
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
		& mphen
		& aphen
		& establish
		& nsapling
		& cmass_litter_leaf
		& cmass_litter_root
		& cmass_litter_sap
		& cmass_litter_heart
		& cmass_litter_repr
		& gcbase
		& gcbase_day
		& wsupply
		& wsupply_leafon
		& fwuptake
		& wstress
		& wstress_day
		& cmass_harvested_products_slow
		& nmass_litter_leaf
		& nmass_litter_root
		& nmass_litter_sap
		& nmass_litter_heart
		& nmass_harvested_products_slow
		& nmass_litter_NH4
		& swindow
		& wmass_deficit_y
		& inund_count
		& inund_stress;
	if (pft.landcover==CROPLAND)
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

/// Gets the growingseason status for crop patch pft. Non-crop patch pft:s always return true.
bool Patchpft::growingseason() const {
	if(cropphen)
		return cropphen->growingseason;
	else
		return true;
}

double Patchpft::total_litter() const {
	return cmass_litter_leaf + cmass_litter_root + cmass_litter_heart + cmass_litter_sap + cmass_litter_repr;
}

/// Update patch pfts root distribution if frozen soil layers
/** Update patch pfts root distribution if frozen soil layers
 * Root distribution follow Jackson et al 1996 and is limited to 
 * maximum thaw depth during the year
 */
void Patchpft::update_rootdist(double maxthaw, double soildepth, double soillyrdepth[NSOILLAYERS_MAX], int nsoillayers) {

	if (hydrology != GERTEN) {

		// Max thaw doesn't exceed first soil layer
		// All roots in top soil layer
		if (maxthaw <= soillyrdepth[0]) {
			rootdist[0] = 1.0;

			for (int lyr = 1; lyr < nsoillayers; lyr++) {
				rootdist[lyr] = 0.0;
			}
		}
		else {

			// Multi layer soil scheme uses Jackson root distribution
			// root_beta is calibrated for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.

			// Standard depth 1500 mm
			const double depth_standard = 1500.0;

			// Standard root depth
			double depth_root = depth_standard;

			// Fraction of fine roots to be homogenously distributed within root depth
			// is zero for MILLER hydrology
			double fhom_root = 0.0;

			if (hydrology == RICHARDS) {

				// Fraction of deep soil, beyond 1500mm, that roots have access to
				double fdeep_root = rootfrac_of_deep_soil;

				// Fraction of fine roots to be homogenously distributed within root depth
				// is set in the ins file
				fhom_root = homrootdistfrac;

				// Grass has half the amount of homogenously distributed roots and depth of deep roots
				if (pft.isgrass()) {
					fdeep_root *= 0.5;
					fhom_root *= 0.5;
				}

				// Moss don't have homogenously distributed roots nor deep roots
				if (pft.ismoss()) {
					fdeep_root = 0.0;
					fhom_root = 0.0;
				}

				// If soil depth is deeper than standard depth, then update root depth
				if (soildepth > depth_standard) {
					depth_root = depth_standard + fdeep_root * (soildepth - depth_standard);
				}
			}

			depth_root = min(maxthaw, depth_root);

			double rootdist_new[NSOILLAYERS_MAX] = { 0.0 };

			pft.calc_rootdist(rootdist_new, soillyrdepth, nsoillayers, depth_root, pft.root_beta, fhom_root);

			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				rootdist[lyr] = rootdist_new[lyr];
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patch member functions
////////////////////////////////////////////////////////////////////////////////


Patch::Patch(int i,Stand& s,Soiltype& st):
	id(i),stand(s),vegetation(*this),soilcolumn(*this, st),fluxes(*this) {

	for (unsigned int p = 0; p < pftlist.nobj; p++) {
		pft.createobj(pftlist[p]);
	}

	age = 0;
	disturbed = false;
	managed = false;
	has_been_cut = false;
	man_strength = 0.0;
	harvest_to_litter = false;
	clearcut_this_year = false;
	cutinterval_actual = 0;
	cutinterval_actual_thisyear = 0;
	managed_this_year = false;
	plant_this_year = false;
	distributed_cutting = false;
	cut_due = false;
	dens_start = 0.0;
	wdemand = 0.0;
	wdemand_leafon = 0.0;
	fpar_grass = 0.0;
	growingseasondays = 0;

	burned = false;
	fire_line_intensity = 0.0;
	fireprob = 0.0;
	ndemand = 0.0;
	dnfert = 0.0;
	anfert = 0.0;
	nharv = 0;
	for (int i = 0; i < NYEARAAET; i++) {
		aaet_5.add(0.0);
	}

	for (int i = 0; i < N_YEAR_BIOMEAVG; i++) {
		fapar_brlt_avg[i]  = 0.0;
		fapar_trbr_avg[i]  = 0.0;
		fapar_grass_avg[i] = 0.0;
		fapar_ndlt_avg[i]  = 0.0;
		fapar_shrub_avg[i] = 0.0;
		fapar_total_avg[i] = 0.0;
	}

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
		& soilcolumn
		& fluxes
		& fpar_grass
		& fpar_ff
		& par_grass_mean
		& nday_growingseason
		& fpc_total
		& disturbed
		& managed
		& has_been_cut
		& cut_due
		& dens_start
		& cutinterval_actual
		& age
		& fireprob
		& growingseasondays
		& intercep
		& aaet_5
		& eet_net_veg
		& wdemand
		& wdemand_day
		& wdemand_leafon
		& fpc_rescale
		& ndemand
		& fire_line_intensity
		& wood_to_atm
		& leaf_to_atm
		& leaf_to_lit
		& wood_to_str
		& wood_to_fwd
		& wood_to_cwd
		& litf_to_atm
		& lfwd_to_atm
		& lcwd_to_atm;
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_grass_avg[i];
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_ndlt_avg[i];
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_brlt_avg[i];
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_trbr_avg[i];
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_shrub_avg[i];
		for (unsigned int i=0; i < N_YEAR_BIOMEAVG; i++)
			arch & fapar_total_avg[i];
}

const Climate& Patch::get_climate() const {
	// All patches within a stand share the same climate
	return stand.get_climate();
}

bool Patch::has_fires() const {
	// Since the standard fire parameterization was not developed for wetland vegetation and wetland/peatland soils, including 
	// fires in tropical peatlands, we disallow this for now.
	return firemodel != NOFIRE && stand.landcover != CROPLAND && stand.landcover != PEATLAND
		&& !(managed && (stand.get_current_management().suppress_fire || suppress_disturbance_in_forestry_stands))
		&& (stand.landcover != PASTURE || disturb_pasture) && stand.landcover != BARREN && stand.landcover != URBAN;
}

bool Patch::has_disturbances() const {
	return ifdisturb && stand.landcover != CROPLAND && !(managed && (stand.get_current_management().suppress_disturbance
		|| suppress_disturbance_in_forestry_stands)) && (stand.landcover != PASTURE || disturb_pasture)
		&& stand.landcover != BARREN && stand.landcover != URBAN;
}

/// C content of patch
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 */
double Patch::ccont(double scale_indiv, bool luc) {

	double ccont = soilcolumn.ccont();

	ccont += soilcolumn.get_sum(&Soillayer::CH4) * KG_PER_G;

	ccont += soilcolumn.snow.ccont();

	for (int i=0; i<npft; i++) {
		Patchpft& ppft = pft[i];
		ccont += ppft.cmass_litter_leaf;
		ccont += ppft.cmass_litter_root;
		ccont += ppft.cmass_litter_sap;
		ccont += ppft.cmass_litter_heart;
		ccont += ppft.cmass_litter_repr;
		ccont += ppft.cmass_harvested_products_slow;
	}

	for (unsigned int i=0; i<vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		ccont += indiv.ccont(scale_indiv, luc);
	}

	return ccont;
}

/// N content of patch
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living N
 *  \param luc 					down-scales living N (used in N balance tests)
 */
double Patch::ncont(double scale_indiv, bool luc) {

	double ncont = soilcolumn.ncont();

	ncont += soilcolumn.snow.ncont();

	for (int i = 0; i<npft; i++) {
		Patchpft& ppft = pft[i];
		ncont += ppft.nmass_litter_leaf;
		ncont += ppft.nmass_litter_root;
		ncont += ppft.nmass_litter_sap;
		ncont += ppft.nmass_litter_heart;
		ncont += ppft.nmass_harvested_products_slow;
		ncont += ppft.nmass_litter_NH4;
	}

	for (unsigned int i = 0; i<vegetation.nobj; i++) {

		Individual& indiv = vegetation[i];
		ncont += indiv.ncont(scale_indiv, luc);
	}

	return ncont;
}

/// Water content of patch
double Patch::water_content() {

	double wcontent = soilcolumn.get_sum(&Soillayer::wmass) + soilcolumn.get_sum(&Soillayer::imass);

	wcontent += soilcolumn.snow.wcontent();

	wcontent += soilcolumn.ssw.wcontent();

	return wcontent;
}

/// C flux of patch
double Patch::cflux() {

	double cflux = 0.0;

	cflux += fluxes.get_annual_flux(Fluxes::LEACHORGC);
	cflux += -fluxes.get_annual_flux(Fluxes::NPP);
	cflux += fluxes.get_annual_flux(Fluxes::REPRC);
	cflux += fluxes.get_annual_flux(Fluxes::SOILC);
	cflux += fluxes.get_annual_flux(Fluxes::FIREC);
	cflux += fluxes.get_annual_flux(Fluxes::ESTC);
	cflux += fluxes.get_annual_flux(Fluxes::SEEDC);
	cflux += fluxes.get_annual_flux(Fluxes::MANUREC);
	cflux += fluxes.get_annual_flux(Fluxes::HARVESTC);
	cflux += fluxes.get_annual_flux(Fluxes::CH4C) * KG_PER_G; // Convert to kg CH4-C m-2 from g CH4-C m-2

	return cflux;
}

/// N flux of patch
double Patch::nflux() {

	double nflux = 0.0;

	nflux += -(stand.get_gridcell().aNH4dep + stand.get_gridcell().aNO3dep);
	nflux += -anfert;
	nflux += -fluxes.get_annual_flux(Fluxes::NFIX);
	nflux += fluxes.get_annual_flux(Fluxes::LEACHMINN);
	nflux += fluxes.get_annual_flux(Fluxes::LEACHORGN);
	nflux += fluxes.get_annual_flux(Fluxes::HARVESTN);
	nflux += fluxes.get_annual_flux(Fluxes::SEEDN);
	nflux += fluxes.get_annual_flux(Fluxes::NH3_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::NOx_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N2O_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N2_FIRE);
	nflux += fluxes.get_annual_flux(Fluxes::N2O_SOIL);
	nflux += fluxes.get_annual_flux(Fluxes::N2_SOIL);
	nflux += fluxes.get_annual_flux(Fluxes::NO_SOIL);
	nflux += fluxes.get_annual_flux(Fluxes::NH3_SOIL);

	return nflux;
}

/// water flux of patch
double Patch::water_flux() {

	double wflux = 0.0;

	// Influx
	wflux += -get_climate().aprec;
	wflux += -fluxes.get_annual_flux(Fluxes::WIRRI);
	wflux += -fluxes.get_annual_flux(Fluxes::WADDED);

	// Outflux
	wflux += fluxes.get_annual_flux(Fluxes::WINTERCEP);
	wflux += fluxes.get_annual_flux(Fluxes::WTRAN);
	wflux += fluxes.get_annual_flux(Fluxes::WEVAP);
	wflux += fluxes.get_annual_flux(Fluxes::WRUNOFF_SURF);
	wflux += fluxes.get_annual_flux(Fluxes::WRUNOFF_SSW);
	wflux += fluxes.get_annual_flux(Fluxes::WRUNOFF_DRAIN);
	wflux += fluxes.get_annual_flux(Fluxes::WRUNOFF_BASE);

	return wflux;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Standpft member functions
////////////////////////////////////////////////////////////////////////////////


void Standpft::serialize(ArchiveStream& arch) {
	arch & cmass_repr
		& anetps_ff_max
		& fpc_total
		& active
		& selection
		& plant
		& reestab
		& plantdensity
		& targetfrac
		& irrigated;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell* gc, Soiltype& st, landcovertype landcoverX, int npatch)
 : id(i),
   gridcell(gc),
   soiltype(st),
   landcover(landcoverX),
   lc_origin(landcoverX),
   st_origin(0),
   frac(1.0) {

	// Constructor: initialises reference member of climate and
	// builds list array of Standpft objects

	if (landcover >= NLANDCOVERTYPES) {
		fail("Unrecognized landcover type\n");
	}

	for (unsigned int p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}

	unsigned int num_patches = 1;
	if (landcover == FOREST || landcover == NATURAL || (disturb_pasture && landcover == PASTURE) || landcover == PEATLAND) {
		// stands with stochastic events
		if (npatch > 0) {
			num_patches = npatch;	// Use patch number provided by calling function
		}
		else {
			num_patches = ::npatch; // Use the global variable npatch
		}
	}

	for (unsigned int p=0;p<num_patches;p++) {
		createobj(*this, soiltype);
	}

	first_year = date.year;
	clone_year = -1;
	transfer_area_st = NULL;
	transfer_area_st = new double[nst];
	for(int i=0;i<nst;i++)
		transfer_area_st[i] = 0.0;
	seed = 12345678;

	stid = 0;
	pftid = -1;
	npft_selection = 0;
	current_rot = 0;
	nyears_in_rotation = 0;
	ndays_in_rotation = 0;
	infallow = false;
	isrotationday = false;
	isirrigated = false;
	hasgrassintercrop = false;
	gdd5_intercrop = 0.0;
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
	lai = 0.0;
	cmass = 0.0;
	cmass_wood = 0.0;
	cmass_mort = 0.0;
	scale_LC_change = 1.0;
}

Stand::~Stand() {

	if (transfer_area_st) {
		delete[] transfer_area_st;
	}
}

double Stand::get_gridcell_fraction() const {
	return frac;
}

/// Initiation of stand variables when run_landcover==true
void Stand::init_stand_lu(StandType& st, double fraction, bool suppress_disturbance) {

	landcovertype lc = st.landcover;
	landcover = lc;

	stid = st.id;
	set_gridcell_fraction(fraction);
	frac_old = 0.0;
	frac_change = fraction;
	gross_frac_increase = fraction;
	st_origin = st.id;

	ManagementType& mt0 = st.get_management(0);

	if(suppress_disturbance || date.get_calendar_year() >= mt0.firstmanageyear) {
		for(unsigned int i=0;i<npatch();i++)
			(*this)[i].managed = true;
	}

	pftid = pftlist.getpftid(mt0.pftname);	// First main crop, will change during crop rotation
	if(pftid < 0) {
			// In case rotation starts with fallow
			if(mt0.fallow && st.rotation.nmanagements > 1) {
				ManagementType& mt_prev = st.get_management(st.rotation.nmanagements - 1);
				pftid = pftlist.getpftid(mt_prev.pftname);
		}
	}
	current_rot = 0;
	set_management();

	if (st.intercrop==NATURALGRASS && ifintercropgrass) {
		hasgrassintercrop = true;

		for (unsigned int i=0; i<pftlist.nobj; i++) {
			if (pftlist[i].isintercropgrass)
				pft[pftlist[i].id].active = true;
		}
	}

	// Set standpft- and patchpft-variables for all crops in a rotation
	if(lc == CROPLAND) {
		for(int rot=0; rot<st.rotation.nmanagements; rot++) {

			ManagementType& mt = st.get_management(rot);
			if(mt.planting_system == "MONOCULTURE") {
				
				int id = pftlist.getpftid(mt.pftname);

				if(id >=0) {
					// Set active pft:s for all management types in rotation
					pft[id].active = true;
				}
				else if(!mt.fallow) {
					fail("Stand type %d pft %s not in pftlist; set to include 1 in instruction file !\n", stid, (char*)mt.pftname);
					break;
				}
			}
		}
	}
}

/// Setting of management parameters for PFT selections at stand creation and forest rotation
void Stand::set_selection_params() {

	ManagementType& mt = get_current_management();
	char selection_cp[200] = {0}, plantdensity_cp[200] = {0}, targetfrac_cp[200] = {0};
	char *p_sel = selection_cp, *p_dens = plantdensity_cp, *p_frac = targetfrac_cp;
	int count_sel = 0, count_dens = 0, count_frac = 0;

	strcpy(selection_cp, mt.selection);
	strcpy(plantdensity_cp, mt.plantdensity);
	strcpy(targetfrac_cp, mt.targetfrac);

	count_sel = split_string(selection_cp);
	count_dens = split_string(plantdensity_cp);
	count_frac = split_string(targetfrac_cp);
	npft_selection = count_sel;

	if(count_dens && count_dens != count_sel || count_frac && count_frac != count_sel)
		fail("Selection parameter number must correspond to number in selection\n");

	for(int i=0;i<count_sel;i++) {
		if(pftlist.getpftid(p_sel) != -1) {
			Standpft& spft = pft[pftlist.getpftid(p_sel)];
			spft.selection = i;
			if(count_dens)
				spft.plantdensity = strtod(p_dens, NULL);
			if(count_frac)
				spft.targetfrac = strtod(p_frac, NULL);
		}
		p_sel += strlen(p_sel) + 1;
		if(count_dens)
			p_dens += strlen(p_dens) + 1;
		if(count_frac)
			p_frac += strlen(p_frac) + 1;
	}
}

/// Setting of management parameters at stand creation and forest rotation
/** Rules for which PFT:s are allowed to establish are set in the instruction file by the parameters landcover
  * (allows all active PFT:s with the same landcovertype), naturalveg (allows none, natural grass or all natural pft:s)
  * and intercrop ("naturalgrass" allows dedicated covercrop grass pft:s).
  * If restrictpfts is true, further restriction of pft:s are specified in the management settings.
  * Rules for reestablishment (after sowing or planting) are set by the parameter reestab, "none", "restricted"(only planted pft:s)
  */
void Stand::set_management() {

	StandType& st = stlist[stid];
	ManagementType& mt = get_current_management();

	// Move variables to stand ?
	if(!readNfert_st)
		gridcell->st[stid].nfert = mt.nfert;
	gridcell->st[stid].diam_cut_low = mt.diam_cut_low;

	if(mt.wateruse == IRRIGATED) {
		isirrigated = true;					
	}
	else {
		isirrigated = false;					
	}

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pftx = pftlist.getobj();
		Standpft& standpft = pft[pftx.id];

		if(mt.wateruse == IRRIGATED) {
			standpft.irrigated = true;
		}
		else {	
			standpft.irrigated = false;
		}
		pftlist.nextobj();			
	}

	bool restrictpfts = (mt.planting_system != "");
	bool naturalveg = st.naturalveg == "ALL";
	bool naturalgrass = st.naturalveg == "ALL" || st.naturalveg == "GRASSONLY";

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pftx = pftlist.getobj();

		if(!restrictpfts && pftx.landcover == landcover
			|| !restrictpfts && naturalveg && pftx.landcover == NATURAL // Allow all natural pft:s to grow in e.g. forests
			|| naturalgrass && pftx.landcover == NATURAL && pftx.lifeform == GRASS // Allow natural grass pft:s to grow in e.g. forests
			|| pftx.landcover == landcover && landcover == FOREST && pftx.lifeform == GRASS) { // Always allow forest grass pft:s to grow in forests

			pft[pftx.id].active = true;
			pft[pftx.id].reestab = true;
			// If restrictpfts = false, plant all tree pft:s after clearcut
			if(pftx.lifeform == TREE)
				pft[pftx.id].plant = true;
			// Test case when PNV planted but not re-established
			if(st.reestab == "NONE" || st.reestab == "") {
				if(pftx.lifeform == TREE)
					pft[pftx.id].reestab = false;
			}
		}
		else {
			pft[pftx.id].active = false;
			pft[pftx.id].reestab = false;
			pft[pftx.id].plant = false;
		}
		pftlist.nextobj();
	}

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pftx = pftlist.getobj();

		for(unsigned int p = 0; p < nobj; p++) {
			Patch& patch = (*this)[p];
			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();
				Patchpft& ppft = patch.pft[indiv.pft.id];
				if(indiv.pft.id == pftx.id) {
					if(clone_year == date.year) {
						// vegetation C transferred during cloning (may be harvested below)
						double stand_frac = get_gridcell_fraction();
						get_gridcell().landcover.acflux_cloned_lc[lc_origin] +=  indiv.ccont() * stand_frac / (double)nobj;
						get_gridcell().landcover.acflux_cloned_lc[landcover] -=  indiv.ccont() * stand_frac / (double)nobj;
					}
				}
				vegetation.nextobj();
			}
		}
		pftlist.nextobj();
	}

	if(mt.cutfirstyear) {
		for(unsigned int p = 0; p < nobj; p++) {
			Patch& patch = (*this)[p];
			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			while (vegetation.isobj) {
				// Program only enters here during LUC when cutfirstyear == 2 (copy_type = CLONESTAND_KILLTREES)
				// and during rotation when cutfirstyear != 0.
				Individual& indiv = vegetation.getobj();
				Patchpft& ppft = patch.pft[indiv.pft.id];
				if(indiv.pft.lifeform == TREE) {
					ppft.cmass_wood_clearcut += check_harvest_cmass(indiv, true);
					ppft.cmass_harv_killed += indiv.ccont();
					harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0.1, clone_year == date.year);
					// frac_cut=1, harv_eff=0.9, res_outtake_twig=0.4, res_outtake_coarse_root=0.1
					indiv.vegetation.killobj();
				}
				else {
					if(mt.killgrass_at_cc) {
						kill_remaining_vegetation(indiv);
						indiv.vegetation.killobj();
					}
					else {
						vegetation.nextobj();
					}
				}
			}
			patch.age = 0;
			patch.plant_this_year = true;
			patch.clearcut_this_year = true;
		}
	}

	// Allow unselected individuals to stay alive if cutfirstyear_unsel is false after rotation or cloning
	for(unsigned int p = 0; p < nobj; p++) {
		Patch& patch = (*this)[p];
		Vegetation& vegetation = patch.vegetation;
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();
			pft[indiv.pft.id].active = true;
			vegetation.nextobj();
		}
	}

	if(!restrictpfts)
		return;

	if(mt.planting_system == "MONOCULTURE") {
		
		int id = pftlist.getpftid(mt.pftname);

		if(id >=0) {

			if(landcover == CROPLAND) {
				pft[id].active = true;

				// Only crop or first crop in rotation
				if(first_year == date.year) {
					// Set crop cycle dates to default values only for first crop in a rotation for a new stand.
					for(unsigned int p = 0; p < nobj; p++) {

						Gridcellpft& gcpft = get_gridcell().pft[id]; 
						Patchpft& ppft = (*this)[p].pft[id];

						ppft.set_cropphen()->sdate = gcpft.sdate_default;
						ppft.set_cropphen()->hlimitdate = gcpft.hlimitdate_default;
				
						if(pftlist[id].phenology == ANY)
							ppft.set_cropphen()->growingseason = true;
						else if(pftlist[id].phenology == CROPGREEN) {
							ppft.set_cropphen()->eicdate = Date::stepfromdate(ppft.get_cropphen()->sdate, -15);
						}
					}
				}
				// At crop rotation
				else {
					pftid = id;
				}

				Standpft& standpft = pft[pftid];
				Gridcellpft& gridcellpft = gridcell->pft[pftid];

				if(mt.wateruse == IRRIGATED) {
					standpft.irrigated = true;
				}
				else {	
					standpft.irrigated = false;
				}

				if(!readNfert)
					gridcellpft.Nfert_read = mt.nfert;
				if(!readsowingdates)
					standpft.sdate_force = mt.sdate;
				if(!readharvestdates)
					standpft.hdate_force = mt.hdate;
			}
			else {
				// Reset and set new values for planting and reestablishment for all pft:s
				pftlist.firstobj();
				while (pftlist.isobj) {
					Pft& pftx = pftlist.getobj();
					Standpft& spft = pft[pftx.id];

					if(pftx.lifeform == TREE) {
						pft[pftx.id].plant = false;
						pft[pftx.id].reestab = false;

						if(st.reestab == "ALL") {
							// Options here are only relevant when using both FOREST and NATURAL tree pfts and they need to be
							// distiguished in their abilities to re-establish when outside the monoculture when naturalveg is "ALL" 
							// (e.g. in the case of planted exotic FOREST pfts).
							// 1. reestablishment by both forest and natural pfts
							// if(pftx.landcover == landcover || st.naturalveg == "ALL" && pftx.landcover == NATURAL) {
							// 2. Reestablishment outside the monoculture only by natural pfts (when naturalveg "ALL") to avoid 
							// reestablishment of exotic trees everywhere.
							if(pftx.landcover == landcover && st.naturalveg != "ALL" || st.naturalveg == "ALL" && pftx.landcover == NATURAL) {
								pft[pftx.id].active = true;
								pft[pftx.id].reestab = true;
							}
						}
						if(mt.cutfirstyear_unsel) {
							for(unsigned int p = 0; p < nobj; p++) {
								Patch& patch = (*this)[p];
								Vegetation& vegetation = patch.vegetation;
								vegetation.firstobj();
								while (vegetation.isobj) {
									Individual& indiv = vegetation.getobj();
									Patchpft& ppft = patch.pft[indiv.pft.id];
									if(indiv.pft.id == pftx.id && pftx.id != id) {
										// cut at cloning (LUC) or at rotation
										ppft.cmass_harv_killed += indiv.ccont();
										harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0, clone_year == date.year);
										// frac_cut=1, harv_eff=0.9, res_outtake_twig=0.4, res_outtake_coarse_root=0
										indiv.vegetation.killobj();
									}
									else {
										vegetation.nextobj();
									}
								}
							}
						}
					}
					pftlist.nextobj();			
				}

				// Set parameters for pft in monoculture
				pft[id].active = true;

				pft[id].plant = true;
				if(st.reestab == "RESTRICTED") {
					pft[id].reestab = true;
				}

				pft[id].plantdensity = mt.plantdensity_pft;
			}
		}
		else if(!mt.fallow) {	// Absence of pftname only in cropland fallow management type.
			fail("Stand type %d pft %s not in pftlist; set to include 1 in instruction file !\n", stid, (char*)mt.pftname);
		}
	}
	else if(mt.planting_system == "SELECTION") {

		if(landcover == CROPLAND) {
			fail("planting system SELECTION not available for CROPLAND stands\n");
		}

		if(mt.selection != "") {

			set_selection_params();

			bool included_pft = false;
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pftx = pftlist.getobj();

				if(mt.pftinselection((const char*)pftx.name)) {
					pft[pftx.id].active = true;

					if(st.reestab == "NONE" || st.reestab == "") {
						pft[pftx.id].reestab = false;
					}
					else {
						pft[pftx.id].reestab = true;
					}
					if(pftx.lifeform == TREE)
						pft[pftx.id].plant = true;

					included_pft = true;
				}
				else if(pftx.lifeform == TREE) {	// Whether grass is allowed is specified in the generic code above

					pft[pftx.id].plant = false;
					pft[pftx.id].reestab = false;

					if(st.reestab == "ALL") {
						// Options here are only relevant when using both FOREST and NATURAL tree pfts and they need to be
						// distiguished in their abilities to re-establish when outside the selection when naturalveg is "ALL" 
						// (e.g. in the case of planted exotic FOREST pfts).
						// 1. Reestablishment by both forest and natural pfts
						// if(pftx.landcover == landcover || st.naturalveg == "ALL" && pftx.landcover == NATURAL) {
						// 2. Reestablishment outside the selection only by natural pfts (when naturalveg "ALL") to avoid 
						// reestablishment of exotic trees everywhere.
						if(pftx.landcover == landcover && st.naturalveg != "ALL" || st.naturalveg == "ALL" && pftx.landcover == NATURAL) {
							pft[pftx.id].active = true;
							pft[pftx.id].reestab = true;
						}
					}
					if(mt.cutfirstyear_unsel) {
						for(unsigned int p = 0; p < nobj; p++) {
							Patch& patch = (*this)[p];
							Vegetation& vegetation = patch.vegetation;
							vegetation.firstobj();
							while (vegetation.isobj) {
								Individual& indiv = vegetation.getobj();
								Patchpft& ppft = patch.pft[indiv.pft.id];
								if(indiv.pft.id == pftx.id) {
									// Cut at cloning (LUC) or at rotation
									ppft.cmass_harv_killed += indiv.ccont();
									harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0, clone_year == date.year);
									// frac_cut=1, harv_eff=0.9, res_outtake_twig=0.4, res_outtake_coarse_root=0
									indiv.vegetation.killobj();
								}
								else {
									vegetation.nextobj();
								}
							}
						}
					}
				}
				pftlist.nextobj();
			}

			if(!included_pft)
				fail("Set all PFTs in selection to include 1 in instruction file !\n");
		}
		else {
			dprintf("Warning: stand type %d planting selection not defined !\n", stid);
		}
	}
	else if(mt.planting_system != "") {

		// Planting systems (pft selections) defined here


		// Functional tree pft classes

		// Special case when both FOREST and NATURAL pfts used, but only FOREST pfts planted
		const bool natural_reestab_not_plant = false;

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pftx = pftlist.getobj();

			if(pftx.lifeform == TREE && (pftx.landcover == landcover || st.naturalveg == "ALL" && pftx.landcover == NATURAL)) {

				if(mt.planting_system == "NEEDLELEAF_EVERGREEN" &&
					pftx.leafphysiognomy == NEEDLELEAF && pftx.phenology == EVERGREEN ||
					mt.planting_system == "NEEDLELEAF_DECIDUOUS" &&
					pftx.leafphysiognomy == NEEDLELEAF && pftx.phenology == SUMMERGREEN ||
					mt.planting_system == "BROADLEAF_EVERGREEN" &&
					pftx.leafphysiognomy == BROADLEAF && pftx.phenology == EVERGREEN ||
					mt.planting_system == "BROADLEAF_DECIDUOUS" &&
					pftx.leafphysiognomy == BROADLEAF && (pftx.phenology == SUMMERGREEN || pftx.phenology == RAINGREEN)
					 && pftx.crownarea_max > 10) {

					pft[pftx.id].active = true;
					if(pftx.landcover == landcover || !natural_reestab_not_plant)					
						pft[pftx.id].plant = true;

					if(!(st.reestab == "NONE" || st.reestab == "")) {
						// Options here are only relevant when using both FOREST and NATURAL tree pfts and they need to be
						// distiguished in their abilities to re-establish within the forest class when naturalveg is "ALL" 
						// (e.g. in the case of planted exotic FOREST pfts).
						// 1. Reestablishment within the forest class by both forest and natural pfts (no conditional)
						// 2. Reestablishment within the forest class only by natural pfts (when naturalveg "ALL") to avoid 
						// reestablishment of exotic FOREST trees everywhere.
						if(pftx.landcover == landcover && st.naturalveg != "ALL" || st.naturalveg == "ALL" && pftx.landcover == NATURAL)
							pft[pftx.id].reestab = true;
					}
				}
				else {
					
					if(st.reestab == "ALL") {
						// 1. Reestablishment outside the forest class by both forest and natural pfts (no conditional)
						// 2. Reestablishment outside the forest class only by natural pfts (when naturalveg "ALL") to  
						// avoid reestablishment of exotic FOREST trees everywhere.
						if(pftx.landcover == landcover && st.naturalveg != "ALL" || st.naturalveg == "ALL" && pftx.landcover == NATURAL) {
							pft[pftx.id].active = true;
							pft[pftx.id].reestab = true;
						}
					}

					if(mt.cutfirstyear_unsel && pftx.lifeform == TREE) {
						for(unsigned int p = 0; p < nobj; p++) {
							Patch& patch = (*this)[p];
							Vegetation& vegetation = patch.vegetation;
							vegetation.firstobj();
							while (vegetation.isobj) {
								Individual& indiv = vegetation.getobj();
								Patchpft& ppft = patch.pft[indiv.pft.id];
								if(indiv.pft.id == pftx.id) {
									// cut at cloning (LUC) or at rotation
									ppft.cmass_harv_killed += indiv.ccont();
									harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0, clone_year == date.year);
									// frac_cut=1, harv_eff=0.9, res_outtake_twig=0.4, res_outtake_coarse_root=0
									indiv.vegetation.killobj();
								}
								else {
									vegetation.nextobj();
								}
							}
						}
					}
				}
			}
			pftlist.nextobj();
		}
	}
}

void Stand::rotate(int rot) {

	StandType& st = stlist[stid];

	if(st.rotation.nmanagements < 2)
		return;

	if(rot > -1)
		current_rot = rot;
	else
		current_rot = (current_rot + 1) % st.rotation.nmanagements;
	ManagementType& mt = get_current_management();

	set_management();

	nyears_in_rotation = 0;
	ndays_in_rotation = 0;
}

double Stand::transfer_area_lc(landcovertype to) {

	double area = 0.0;

	if (transfer_area_st) {

		for (int j=0; j<nst; j++) {

			if (stlist[j].landcover == to)
				area += transfer_area_st[j];
		}
	}
	return area;
}

double Stand::ccont(double scale_indiv) {

	double ccont = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		ccont += (*this)[p].ccont(scale_indiv) / nobj;

	return ccont;
}

double Stand::ncont(double scale_indiv) {

	double ncont = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		ncont += (*this)[p].ncont(scale_indiv) / nobj;

	return ncont;
}

double Stand::water_content() {

	double water_content = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		water_content += (*this)[p].water_content() / nobj;

	return water_content;
}

double Stand::cflux() {

	double cflux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		cflux += (*this)[p].cflux() / nobj;

	return cflux;
}

double Stand::nflux() {

	double nflux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		nflux += (*this)[p].nflux() / nobj;

	return nflux;
}

double Stand::water_flux() {

	double water_flux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		water_flux += (*this)[p].water_flux() / nobj;

	return water_flux;
}

/// Returns true if stand is true high-latitude peatland stand, as opposed to a wetland < PEATLAND_WETLAND_LATITUDE_LIMIT N
bool Stand::is_highlatitude_peatland_stand() const {

	double lat = gridcell->get_lat();

	return landcover==PEATLAND && lat >= PEATLAND_WETLAND_LATITUDE_LIMIT;
}

/// Returns true if stand is wetland stand, as opposed to a peatland >= PEATLAND_WETLAND_LATITUDE_LIMIT N
bool Stand::is_true_wetland_stand() const {

	double lat = gridcell->get_lat();

	return landcover==PEATLAND && lat < PEATLAND_WETLAND_LATITUDE_LIMIT;
}

Stand& Stand::clone(StandType& st, double fraction, bool suppress_disturbance) {

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

	// reset managed before setting in init_stand_lu()
	for (unsigned int p = 0; p < nobj; p++) {
//		new_stand[p].age = 0;				// probably not what we want
		if(st.landcover == NATURAL)
			new_stand[p].managed = false;
	}

	// Set land use settings for new stand
	new_stand.init_stand_lu(st, fraction, suppress_disturbance);

	return new_stand;
}

double Stand::get_landcover_fraction() const {
	if (get_gridcell().landcover.frac[landcover])
		return frac / get_gridcell().landcover.frac[landcover];
	else
		return 0.0;
}

double Stand::get_distinterval () const { 
	return get_gridcell().st[stid].distinterval_st;
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
			Patch& patch = createobj(*this, soiltype);
			arch & patch;
		}
	}

	arch & first_year
		& clone_year
		& frac
		& stid
		& pftid
		& npft_selection
		& current_rot
		& nyears_in_rotation
		& ndays_in_rotation
		& infallow
		& isirrigated
		& hasgrassintercrop
		& gdd5_intercrop
		& cloned
		& lc_origin
		& st_origin
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
// Implementation of cropindiv_struct member functions
////////////////////////////////////////////////////////////////////////////////

void cropindiv_struct::serialize(ArchiveStream& arch) {
	arch & grs_cmass_plant
		& grs_cmass_leaf
		& grs_cmass_root
		& grs_cmass_ho
		& grs_cmass_agpool
		& grs_cmass_dead_leaf
		& grs_cmass_stem
		& cmass_leaf_sen
		& nmass_ho
		& nmass_agpool
		& nmass_dead_leaf
		& isintercropgrass;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp              = 0.0;
	agpp              = 0.0;
	aet               = 0.0;
	fpc               = 0.0;
	fpc_daily		  = 0.0;
	densindiv         = 0.0;
	cmass_leaf        = 0.0;
	cmass_root        = 0.0;
	cmass_sap         = 0.0;
	cmass_heart       = 0.0;
	cmass_debt        = 0.0;
	cmass_wood_old	  = 0.0;
	cmass_veg         = 0.0;
	cmass_leaf_post_turnover      = 0.0;
	cmass_root_post_turnover      = 0.0;
	cmass_tot_luc     = 0.0;
	phen              = 0.0;
	aphen             = 0.0;
	deltafpc          = 0.0;

	// Nitrogen
	nmass_leaf        = 0.0;
	nmass_leaf_luc    = 0.0;
	nmass_root        = 0.0;
	nmass_sap         = 0.0;
	nmass_heart       = 0.0;
	cton_leaf_aopt    = 0.0;
	cton_leaf_aavr    = 0.0;
	cton_status       = 0.0;
	nmass_veg         = 0.0;
	nmass_tot_luc     = 0.0;

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

	n_opt_isabovelim  = false;
	nstress           = false;

	// additional initialisation
	age               = 0.0;
	height            = 0.0;
	diam              = 0.0;
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

	int m;
	for (m=0; m<12; m++) {
		mlai[m] = 0.0;
		mlai_max[m] = 0.0;
		mgpp_no[m] = 0.0;
		mvmaxnlim[m] = 0.0;
		mlambda[m] = 0.0;
	}

	// bvoc
	iso               = 0.;
	fvocseas          = 1.;
	for (int im=0; im<NMTCOMPOUNDS; im++){
		mon[im]		= 0.;
		monstor[im]	= 0.;
	}

	dnpp              = 0.0;
	dgpp_no           = 0.0;
	dgpp_ns           = 0.0;
	agpp_no           = 0.0;
	agpp_ns           = 0.0;
	cropindiv         = NULL;
	last_turnover_day = -1;

	Stand& stand = vegetation.patch.stand;

	if (pft.landcover==CROPLAND) {
		cropindiv=new cropindiv_struct;

		if (stand.pftid == pft.id) {
			cropindiv->isprimarycrop = true;
		}
		else if (stand.hasgrassintercrop && pft.isintercropgrass) {	// grass cover crop growth
			cropindiv->isintercropgrass = true;
		}
	}

	man_strength = 0.0;
}

void Individual::serialize(ArchiveStream& arch) {
	arch & cmass_leaf
		& cmass_root
		& cmass_sap
		& cmass_heart
		& cmass_debt
		& cmass_veg
		& cmass_leaf_post_turnover
		& cmass_root_post_turnover
		& last_turnover_day
		& fpc
		& fpc_daily
		& fpar
		& densindiv
		& phen
		& aphen
		& aphen_raingreen
		& anpp
		& agpp
		& aet
		& aaet
		& ltor
		& height
		& diam
		& crownarea
		& deltafpc
		& boleht
		& lai
		& lai_layer
		& lai_indiv
		& lai_daily
		& lai_indiv_daily
		& greff_5
		& cmass_wood_inc_5
		& age
		& mlai
		& fpar_leafon
		& lai_leafon_layer
		& intercep
		& phen_mean
		& wstress
		& alive
		& monstor
		& fvocseas
		& photosynthesis

		// Nitrogen
		& nmass_leaf
		& nmass_root
		& nmass_sap
		& nmass_heart
		& cmass_wood_old
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
		& nday_leafon;

	if (pft.landcover==CROPLAND)
		arch & *cropindiv;
}

Individual::~Individual() {
	if (cropindiv)
		delete cropindiv;
}

/// Access functions for cropindiv
cropindiv_struct* Individual::get_cropindiv() const {
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

/// Stored N is partioned out to make sure tissue C:N ratios are within limits
/**
 *  Help function for reduce_biomass(). Stored N is partioned out to make sure 
 *  tissue C:N ratios are within limits
 */ // DW_COMMENT don't I need & before variables?
void nstore_adjust(double cmass_leaf_litter, double cmass_root_litter, double cmass_sap_litter, double cmass_heart_litter,
	double& nmass_leaf_litter, double& nmass_root_litter, double& nmass_sap_litter, double& nmass_heart_litter,
	double cton_leaf_min, double cton_root_min, double cton_sap_min, double cton_heart_min,
	double cton_leaf_max, double cton_root_max, double cton_sap_max, double cton_heart_max,
	double nstore, double& NH4_mass) {

	double nbefore = nmass_leaf_litter + nmass_root_litter + nmass_sap_litter + nmass_heart_litter + nstore + NH4_mass;

	// If compartment are above maximum allowed N concentration, move excess N mass to storage
	if (nmass_leaf_litter > cmass_leaf_litter / cton_leaf_min) {
		nstore += nmass_leaf_litter - cmass_leaf_litter / cton_leaf_min;
		nmass_leaf_litter = cmass_leaf_litter / cton_leaf_min;
	}
	if (nmass_root_litter > cmass_root_litter / cton_root_min) {
		nstore += nmass_root_litter - cmass_root_litter / cton_root_min;
		nmass_root_litter = cmass_root_litter / cton_root_min;
	}
	if (nmass_sap_litter > cmass_sap_litter / cton_sap_min) {
		nstore += nmass_sap_litter - cmass_sap_litter / cton_sap_min;
		nmass_sap_litter = cmass_sap_litter / cton_sap_min;
	}
	if (nmass_heart_litter > cmass_heart_litter / cton_heart_min) {
		nstore += nmass_heart_litter - cmass_heart_litter / cton_heart_min;
		nmass_heart_litter = cmass_heart_litter / cton_heart_min;
	}

	// N to partition
	double ndist = nstore;

	// N demand to maximum N concentration in different compartments
	double ndemand_leaf = max(0.0, cmass_leaf_litter / cton_leaf_min - nmass_leaf_litter);
	double ndemand_root = max(0.0, cmass_root_litter / cton_root_min - nmass_root_litter);
	double ndemand_sap = max(0.0, cmass_sap_litter / cton_sap_min - nmass_sap_litter);
	double ndemand_heart = max(0.0, cmass_heart_litter / cton_heart_min - nmass_heart_litter);

	// Total N demand
	double ndemand = ndemand_leaf + ndemand_root + ndemand_sap + ndemand_heart;

	double nfrac = 0.0;

	if (ndemand > 0.0) {

		// N demand fraction that is meet
		nfrac = min(1.0, ndist / ndemand);

		// Subtract N demand from N to partition
		ndist -= nfrac * ndemand;

		// Partition stored N
		nmass_leaf_litter += ndemand_leaf * nfrac;
		nmass_root_litter += ndemand_root * nfrac;
		nmass_sap_litter += ndemand_sap * nfrac;
		nmass_heart_litter += ndemand_heart * nfrac;
	}

	// Add leftover stored N to NH4 pool
	NH4_mass += ndist;

	double nafter = nmass_leaf_litter + nmass_root_litter + nmass_sap_litter + nmass_heart_litter + NH4_mass;

	if (fabs(nbefore - nafter) > 1.0e-16)
		dprintf("Year %d ADJNLIM Diff %g B %g A %g\n", date.year, nbefore - nafter, nbefore, nafter);
}

void Individual::reduce_biomass(double mortality, double mortality_fire) {

	// This function needs to be modified if a new lifeform is added,
	// specifically to deal with nstore().
	assert(pft.lifeform == TREE || pft.lifeform == GRASS || pft.lifeform == MOSS);

	if (!negligible(mortality)) {

		const double mortality_non_fire = mortality - mortality_fire;
		const double mortality_fire_wood = mortality_fire * pft.ag_woody_frac;
		const double mortality_non_fire_wood = mortality - mortality_fire_wood;

		// Transfer killed biomass to litter
		// (above-ground biomass killed by fire enters atmosphere, not litter)

		Patchpft& ppft = patchpft();

		double cmass_leaf_litter = mortality * cmass_leaf;
		double cmass_root_litter = mortality * cmass_root;

		if (pft.landcover==CROPLAND) {
			if (pft.aboveground_ho)
				cmass_leaf_litter += mortality * cropindiv->cmass_ho;
			else
				cmass_root_litter += mortality * cropindiv->cmass_ho;

			cmass_leaf_litter += mortality * cropindiv->cmass_agpool;
		}

		ppft.cmass_litter_leaf += cmass_leaf_litter * mortality_non_fire / mortality;
		ppft.cmass_litter_root += cmass_root_litter;

		double cmass_sap_litter = 0.0;
		double cmass_heart_litter = 0.0;

		if (cmass_debt <= cmass_heart + cmass_sap) {
			if (cmass_debt <= cmass_heart) {
				cmass_sap_litter = mortality_non_fire_wood * cmass_sap;
				cmass_heart_litter = mortality_non_fire_wood * (cmass_heart - cmass_debt);
			}
			else {
				cmass_sap_litter = mortality_non_fire_wood * (cmass_sap + cmass_heart - cmass_debt);
			}
		}
		else {
			double debt_excess = mortality_non_fire_wood * (cmass_debt - (cmass_sap + cmass_heart));
			report_flux(Fluxes::NPP, debt_excess);
			report_flux(Fluxes::RA, -debt_excess);
		}

		ppft.cmass_litter_sap += cmass_sap_litter;
		ppft.cmass_litter_heart += cmass_heart_litter;

		// Nitrogen
		double nmass_leaf_litter = mortality * nmass_leaf;
		double nmass_root_litter = mortality * nmass_root;
		double nmass_sap_litter = mortality * nmass_sap;
		double nmass_heart_litter = mortality * nmass_heart;

		if (pft.landcover==CROPLAND) {
			if (pft.aboveground_ho)
				nmass_leaf_litter += mortality * cropindiv->nmass_ho;
			else
				nmass_root_litter += mortality * cropindiv->nmass_ho;

			nmass_leaf_litter += mortality * cropindiv->nmass_agpool;
		}

		// Stored N is partioned out to make sure tissue C:N ratios are within limits
		nstore_adjust(cmass_leaf_litter, cmass_root_litter, cmass_sap_litter, cmass_heart_litter,
			nmass_leaf_litter, nmass_root_litter, nmass_sap_litter, nmass_heart_litter,
			pft.cton_leaf_min, pft.cton_root_min, pft.cton_sap_min, pft.cton_heart_min,
			pft.cton_leaf_max, pft.cton_root_max, pft.cton_sap_max, pft.cton_heart_max,
			mortality * nstore(), ppft.nmass_litter_NH4);

		ppft.nmass_litter_leaf  += nmass_leaf_litter * mortality_non_fire / mortality;
		ppft.nmass_litter_root  += nmass_root_litter;
		ppft.nmass_litter_sap   += nmass_sap_litter * mortality_non_fire_wood / mortality;
		ppft.nmass_litter_heart += nmass_heart_litter * mortality_non_fire_wood / mortality;

		// Flux to atmosphere from burned above-ground biomass

		double cflux_fire = (mortality_fire * cmass_leaf_litter + mortality_fire_wood * (cmass_sap_litter + cmass_heart_litter)) / mortality;
		double nflux_fire = (mortality_fire * nmass_leaf_litter + mortality_fire_wood * (nmass_sap_litter + nmass_heart_litter)) / mortality;

		// Fire C fluxes
		report_flux(Fluxes::FIREC,    cflux_fire);

		// Fire N fluxes
		report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire);
		report_flux(Fluxes::NOx_FIRE, Fluxes::NOx_FIRERATIO * nflux_fire);
		report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);
		report_flux(Fluxes::N2_FIRE,  Fluxes::N2_FIRERATIO  * nflux_fire);

		// Reduce this Individual's biomass values

		ppft.cmass_fire += mortality_fire * ccont();

		const double remaining = 1.0 - mortality;

		if (pft.lifeform != GRASS && pft.lifeform != MOSS) {
			densindiv *= remaining;
		}

		cmass_leaf      *= remaining;
		cmass_root      *= remaining;
		cmass_sap       *= remaining;
		cmass_heart     *= remaining;
		cmass_debt      *= remaining;
		if (pft.landcover==CROPLAND) {
			cropindiv->cmass_ho *= remaining;
			cropindiv->cmass_agpool *= remaining;
		}
		nmass_leaf      *= remaining;
		nmass_root      *= remaining;
		nmass_sap       *= remaining;
		nmass_heart     *= remaining;
		nstore_longterm *= remaining;
		nstore_labile   *= remaining;
		if (pft.landcover==CROPLAND) {
			cropindiv->nmass_ho *= remaining;
			cropindiv->nmass_agpool *= remaining;
		}
	}
}

double Individual::cton_leaf(bool use_phen /* = true*/) const {

	Stand& stand = vegetation.patch.stand;

	if (stand.is_true_crop_stand() && !negligible(cmass_leaf_today()) && !negligible(nmass_leaf)) {
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

double Individual::cton_root(bool use_phen /* = true*/) const {

	if (!negligible(cmass_root) && !negligible(nmass_root)) {
		if (use_phen) {
			if (!negligible(cmass_root_today())) {
				return max(pft.cton_root_min, cmass_root_today() / nmass_root);
			}
			else {
				return pft.cton_root_avr;
			}
		}
		else {
			return max(pft.cton_root_min, cmass_root / nmass_root);
		}
	}
	else {
		return pft.cton_root_max;
	}
}

double Individual::cton_sap() const {

	if (pft.lifeform == TREE) {
		if (!negligible(cmass_sap) && !negligible(nmass_sap))
			return max(pft.cton_sap_min, cmass_sap / nmass_sap);
		else
			return pft.cton_sap_max;
	}
	else {
		return 1.0;
	}
}

double Individual::cton_heart() const {

	if (pft.lifeform == TREE) {
		if (!negligible(cmass_heart) && !negligible(nmass_heart))
			return max(pft.cton_heart_min, cmass_heart / nmass_heart);
		else
			return pft.cton_heart_max;
	}
	else {
		return 1.0;
	}
}

/// C content of individual
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 */
double Individual::ccont(double scale_indiv, bool luc) const {

	double ccont = 0.0;

	if (alive || istruecrop_or_intercropgrass()) {

		if (has_daily_turnover()) {	// Not taking into account future daily wood allocation/turnover

			if (cropindiv) {

				if (luc) {
					ccont += cropindiv->grs_cmass_leaf - cropindiv->grs_cmass_leaf_luc * (1.0 - scale_indiv);
					ccont += cropindiv->grs_cmass_root - cropindiv->grs_cmass_root_luc * (1.0 - scale_indiv);
				}
				else {
					ccont += cropindiv->grs_cmass_leaf * scale_indiv;
					ccont += cropindiv->grs_cmass_root * scale_indiv;
				}

				if (pft.phenology == CROPGREEN) {

					if (luc) {
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

			ccont += (cmass_leaf + cmass_root + cmass_sap + cmass_heart - cmass_debt) * scale_indiv;

			if (pft.landcover == CROPLAND) {
				ccont += (cropindiv->cmass_ho + cropindiv->cmass_agpool) * scale_indiv;
				// Yearly allocation not defined for crops with nlim
			}
		}
	}

	return ccont;
}

/// N content of individual
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living N
 *  \param luc 					down-scales living N (used in C balance tests)
 */
double Individual::ncont(double scale_indiv, bool luc) const {

	double ncont = 0.0;

	if (luc) {

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

	if (pft.landcover == CROPLAND) {

		if (luc) {
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

/// Whether grass growth is uninterrupted by crop growth.
bool Individual::continous_grass() const {

	if (pft.landcover != CROPLAND) {
		return false;
	}

	Stand& stand = vegetation.patch.stand;
	StandType& st = stlist[stand.stid];
	bool sowing_restriction = true;

	for (int i=0; i<st.rotation.nmanagements; i++) {
		int pftid = pftlist.getpftid(st.get_management(i).pftname);
		if (pftid > -1 && !stand.get_gridcell().pft[pftid].sowing_restriction) {
			sowing_restriction = false;
		}
	}

	return cropindiv->isintercropgrass && sowing_restriction;
}

double Individual::ndemand_storage(double cton_leaf_opt) {

	if (vegetation.patch.stand.is_true_crop_stand() && ifnlim)	// only CROPGREEN, only ifnlim ?
		// analogous with root demand
		storendemand = max(0.0, cropindiv->grs_cmass_stem / (cton_leaf_opt * pft.cton_stem_avr / pft.cton_leaf_avr) - cropindiv->nmass_agpool);
	else
		storendemand = max(0.0, min(max(0.0, agpp_no) * scale_n_storage, max_n_storage) - nstore());

	return storendemand;
}

/// Checks C mass and zeroes any negative value, balancing by adding to npp and reducing respiration
double Individual::check_C_mass() {

	if (pft.landcover != CROPLAND)
		return 0;

	double negative_cmass = 0.0;

	if (cropindiv->grs_cmass_leaf < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_leaf;
		cropindiv->ycmass_leaf -= cropindiv->grs_cmass_leaf;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_leaf;
		cropindiv->grs_cmass_leaf = 0.0;
	}
	if (cropindiv->grs_cmass_root < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_root;
		cropindiv->ycmass_root -= cropindiv->grs_cmass_root;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_root;
		cropindiv->grs_cmass_root = 0.0;
	}
	if (cropindiv->grs_cmass_ho < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_ho;
		cropindiv->ycmass_ho -= cropindiv->grs_cmass_ho;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_ho;
		cropindiv->grs_cmass_ho = 0.0;
	}
	if (cropindiv->grs_cmass_agpool < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_agpool;
		cropindiv->ycmass_agpool -= cropindiv->grs_cmass_agpool;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_agpool;
		cropindiv->grs_cmass_agpool = 0.0;
	}
	if (cropindiv->grs_cmass_dead_leaf < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->ycmass_dead_leaf -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_dead_leaf;
		cropindiv->grs_cmass_dead_leaf = 0.0;
	}
	if (cropindiv->grs_cmass_stem < 0.0) {
		negative_cmass -= cropindiv->grs_cmass_stem;
		cropindiv->ycmass_stem -= cropindiv->grs_cmass_stem;
		cropindiv->grs_cmass_plant -= cropindiv->grs_cmass_stem;
		cropindiv->grs_cmass_stem = 0.0;
	}

	if (largerthanzero(negative_cmass, -14)) {
		anpp += negative_cmass;
		report_flux(Fluxes::NPP, negative_cmass);
		report_flux(Fluxes::RA, -negative_cmass);
	}

	return negative_cmass;
}

/// Checks N mass and zeroes any negative value, balancing by reducing N mass of other organs and (if needed) reducing anflux_landuse_change
double Individual::check_N_mass() {

	if (pft.landcover != CROPLAND && pft.landcover != PASTURE)
		return 0;

	double negative_nmass = 0.0;

	if (nmass_leaf < 0.0) {
		negative_nmass -= nmass_leaf;
		if (cropindiv)
			cropindiv->ynmass_leaf -= nmass_leaf;
		nmass_leaf = 0.0;
	}
	if (nmass_root < 0.0) {
		negative_nmass -= nmass_root;
		if (cropindiv)
			cropindiv->ynmass_root -= nmass_root;
		nmass_root = 0.0;
	}
	if (cropindiv) {
		if (cropindiv->nmass_ho < 0.0) {
			negative_nmass -= cropindiv->nmass_ho;
			cropindiv->ynmass_ho -= cropindiv->nmass_ho;
			cropindiv->nmass_ho = 0.0;
		}
		if (cropindiv->nmass_agpool < 0.0) {
			negative_nmass -= cropindiv->nmass_agpool;
			cropindiv->ynmass_agpool -= cropindiv->nmass_agpool;
			cropindiv->nmass_agpool = 0.0;
		}
		if (cropindiv->nmass_dead_leaf < 0.0) {
			negative_nmass -= cropindiv->nmass_dead_leaf;
			cropindiv->ynmass_dead_leaf -= cropindiv->nmass_dead_leaf;
			cropindiv->nmass_dead_leaf = 0.0;
		}
	}
	if (nstore_labile < 0.0) {
		negative_nmass -= nstore_labile;
		nstore_labile = 0.0;
	}
	if (nstore_longterm < 0.0) {
		negative_nmass -= nstore_longterm;
		nstore_longterm = 0.0;
	}

	if (largerthanzero(negative_nmass, -14)) {
		double pos_nmass = ncont();
		if (pos_nmass > negative_nmass) {
			nmass_leaf -= negative_nmass * nmass_leaf / pos_nmass;
			nmass_root -= negative_nmass * nmass_root / pos_nmass;
			if (cropindiv) {
				cropindiv->nmass_ho -= negative_nmass * cropindiv->nmass_ho / pos_nmass;
				cropindiv->nmass_agpool -= negative_nmass * cropindiv->nmass_agpool / pos_nmass;
				cropindiv->nmass_dead_leaf -= negative_nmass * cropindiv->nmass_dead_leaf / pos_nmass;
			}
		}
		else {
			vegetation.patch.stand.get_gridcell().landcover.anflux_landuse_change -= (negative_nmass - pos_nmass) * vegetation.patch.stand.get_gridcell_fraction();
			nmass_leaf = 0.0;
			nmass_root = 0.0;
			if (cropindiv) {
				cropindiv->nmass_ho = 0.0;
				cropindiv->nmass_agpool = 0.0;
				cropindiv->nmass_dead_leaf = 0.0;
			}
		}
	}

	return negative_nmass;
}

/// Whether resetting of grs_cmass and turnover (if has_daily_turnover() returns true) of continuous grass is to be done this day.
bool Individual::is_turnover_day() const {

	if (patchpft().cropphen && patchpft().cropphen->growingseason) {

		const Climate& climate = vegetation.patch.get_climate();

		return date.day == climate.testday_prec;
	}
	else {
		return false;
	}
}

Patchpft& Individual::patchpft() const {
	return vegetation.patch.pft[pft.id];
}

/// Save cmass-values on first day of the year of land cover change in expanding stands
void Individual::save_cmass_luc() {

	cmass_tot_luc = 0.0;

	if (cropindiv) {
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

	nmass_leaf_luc = nmass_leaf;
	nmass_root_luc = nmass_root;
	nmass_sap_luc = nmass_sap;
	nmass_heart_luc = nmass_heart;
	nstore_longterm_luc = nstore_longterm;
	nstore_labile_luc = nstore_labile;

	if (cropindiv) {
		cropindiv->nmass_ho_luc = cropindiv->nmass_ho;
		cropindiv->nmass_agpool_luc = cropindiv->nmass_agpool;
		cropindiv->nmass_dead_leaf_luc = cropindiv->nmass_dead_leaf;
	}
	nmass_tot_luc = ncont();
}

/// Gets the individual's daily cmass_leaf value
double Individual::cmass_leaf_today() const {

	if (istruecrop_or_intercropgrass())
		return patchpft().cropphen->growingseason ? cropindiv->grs_cmass_leaf : 0;
	else
		return cmass_leaf * phen;
}

/// Gets the individual's daily cmass_root value
double Individual::cmass_root_today() const {

	if (istruecrop_or_intercropgrass())
		return patchpft().cropphen->growingseason ? cropindiv->grs_cmass_root : 0;
	else
		return cmass_root * phen;
}

/// Gets the individual's daily fpc value
double Individual::fpc_today() const {

	if (pft.phenology == CROPGREEN) {
		return patchpft().cropphen->growingseason ? fpc_daily : 0.0;
	}
	else if (phen == 1.0) {
		return fpc;
	}
	else if (phen == 0.0) {
		return 0.0;
	}
	else {
		if (pft.lifeform == TREE) {
			return crownarea * densindiv * (1.0 - lambertbeer(Kb, lai_indiv_today()));
		}
		else {
			// Grass and mosses
			return 1.0 - lambertbeer(Kb, lai_indiv_today());
		}
	}
}

/// Gets the individual's daily lai value
double Individual::lai_today() const {

	if (pft.phenology == CROPGREEN)
		return patchpft().cropphen->growingseason ? lai_daily : 0;
	else
		return lai * phen;
}

/// Gets the individual's daily lai_indiv value
double Individual::lai_indiv_today() const {

	if (pft.phenology == CROPGREEN)
		return patchpft().cropphen->growingseason ? lai_indiv_daily : 0;
	else
		return lai_indiv * phen;
}

/// Gets the Nitrigen limited LAI
double Individual::lai_nitrogen_today() const{
	if (pft.phenology==CROPGREEN) {

		double Ln = 0.0;
		if (patchpft().cropphen->growingseason && cmass_leaf_today() > 0.0) {
			const double k = 0.5;
			const double ktn = 0.52*k + 0.01; // Yin et al 2003
			double nb = 1/(pft.cton_leaf_max*pft.sla);
			Ln = (1/ktn) * log(1+ktn*nmass_leaf/nb);
		}
		return Ln;
	}
	else {
		return 1.0;
	}
}

/// Gets the growingseason status for crop individual. Non-crop individuals always return true.
bool Individual::growingseason() const {
	return patchpft().cropphen ? patchpft().cropphen->growingseason : true;
}

/// Whether harvest and turnover is done on actual C and N on harvest or turnover day, which can occur any day of the year.
bool Individual::has_daily_turnover() const {
	return istruecrop_or_intercropgrass();
}

/// Partition stored N between plant compartments to prevent unrealistic C:N ratios in litter when killed
void Individual::partition_nstore(double& NH4_mass) {

	double nbefore = ncont() + NH4_mass;

	// If compartment are above maximum allowed N concentration, move excess N mass to storage
	if (nmass_leaf > cmass_leaf / pft.cton_leaf_min) {
		nstore_longterm += nmass_leaf - cmass_leaf / pft.cton_leaf_min;
		nmass_leaf = cmass_leaf / pft.cton_leaf_min;
	}
	if (nmass_root > cmass_root / pft.cton_root_min) {
		nstore_longterm += nmass_root - cmass_root / pft.cton_root_min;
		nmass_root = cmass_root / pft.cton_root_min;
	}
	if (nmass_sap > cmass_sap / pft.cton_sap_min) {
		nstore_longterm += nmass_sap - cmass_sap / pft.cton_sap_min;
		nmass_sap = cmass_sap / pft.cton_sap_min;
	}
	if (nmass_heart > cmass_heart / pft.cton_heart_min) {
		nstore_longterm += nmass_heart - cmass_heart / pft.cton_heart_min;
		nmass_heart = cmass_heart / pft.cton_heart_min;
	}

	// N to partition
	double ndist = nstore();

	// N demand to maximum N concentration is different compartments
	double ndemand_leaf = max(0.0, cmass_leaf / pft.cton_leaf_min - nmass_leaf);
	double ndemand_root = max(0.0, cmass_root / pft.cton_root_min - nmass_root);
	double ndemand_sap = max(0.0, cmass_sap / pft.cton_sap_min - nmass_sap);
	double ndemand_heart = max(0.0, cmass_heart / (2.0 * pft.cton_sap_min) - nmass_heart);

	// Total N demand
	double ndemand = ndemand_leaf + ndemand_root + ndemand_sap + ndemand_heart;

	if (ndemand > 0.0) {

		// N demand fraction that is meet
		double nfrac = min(1.0, ndist / ndemand);

		// Subtract N demand from N to partition
		ndist -= nfrac * ndemand;

		// Partition stored N
		nmass_leaf += ndemand_leaf * nfrac;
		nmass_root += ndemand_root * nfrac;
		nmass_sap += ndemand_sap * nfrac;
		nmass_heart += ndemand_heart * nfrac;
	}

	// Reset storage
	nstore_longterm = 0.0;
	nstore_labile = 0.0;

	// Add leftover stored N to NH4 pool
	NH4_mass += ndist;

	double nafter = ncont() + NH4_mass;

	if (fabs(nbefore - nafter) > 1.0e-16)
		dprintf("Year %d %s IID %d Diff %g B %g A %g\n",date.year,(char*)pft.name, id, nbefore - nafter, nbefore, nafter);
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

	double clitter_sap = 0.0;
	double clitter_heart = 0.0;

	// Distribute stored N between plant compartments to prevent
	// unrealistic C:N ratios in litter when killed
	partition_nstore(ppft.nmass_litter_NH4);

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
	if (alive || istruecrop_or_intercropgrass()) {

		// For leaf and root, catches small, negative values too

		// Leaf: remove residue outtake and send the rest to litter
		if (has_daily_turnover() && cropindiv) {

			if (pft.lifeform == GRASS && pft.phenology != CROPGREEN) {
				charvest_flux += cropindiv->grs_cmass_leaf * harv_eff;
				cropindiv->grs_cmass_leaf *= (1 - harv_eff);
			}

			ppft.cmass_litter_leaf += cropindiv->grs_cmass_leaf * (1 - res_outtake);
			charvest_flux    += cropindiv->grs_cmass_leaf * res_outtake;
		}
		else {

			if (pft.lifeform == GRASS && pft.phenology != CROPGREEN) {
				charvest_flux += cmass_leaf * harv_eff;
				cmass_leaf *= (1 - harv_eff);
			}
			ppft.cmass_litter_leaf += cmass_leaf * (1 - res_outtake);
			charvest_flux    += cmass_leaf * res_outtake;
		}
		// Root: all goes to litter
		if (has_daily_turnover() && cropindiv)
			ppft.cmass_litter_root += cropindiv->grs_cmass_root;
		else
			ppft.cmass_litter_root += cmass_root;

		if (pft.landcover == CROPLAND) {

			if (has_daily_turnover()) {

				charvest_flux += cropindiv->grs_cmass_ho * harv_eff;
				cropindiv->grs_cmass_ho *= (1 - harv_eff);

				if (pft.aboveground_ho) {
					ppft.cmass_litter_leaf+=cropindiv->grs_cmass_ho * (1 - res_outtake);
					charvest_flux += cropindiv->grs_cmass_ho * res_outtake;
				}
				else {
					ppft.cmass_litter_root+=cropindiv->grs_cmass_ho;
				}
				ppft.cmass_litter_leaf+=cropindiv->grs_cmass_agpool * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_agpool * res_outtake;

				ppft.cmass_litter_leaf+=cropindiv->grs_cmass_dead_leaf * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_dead_leaf * res_outtake;

				ppft.cmass_litter_leaf+=cropindiv->grs_cmass_stem * (1 - res_outtake);
				charvest_flux += cropindiv->grs_cmass_stem * res_outtake;
			}
			else {

				charvest_flux += cropindiv->cmass_ho * harv_eff;
				cropindiv->cmass_ho *= (1 - harv_eff);

				if (pft.aboveground_ho) {
					ppft.cmass_litter_leaf+=cropindiv->cmass_ho * (1 - res_outtake);
					charvest_flux += cropindiv->cmass_ho * res_outtake;
				}
				else {
					ppft.cmass_litter_root+=cropindiv->cmass_ho;
				}
				ppft.cmass_litter_leaf+=cropindiv->cmass_agpool * (1 - res_outtake);
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

				double cwood_harvest;

				partition_wood_biomass(to_partition_sap, to_partition_heart,
				                       harv_eff, harvest_slow_frac, res_outtake,
				                       clitter_sap, clitter_heart,
				                       cwood_harvest, charvested_products_slow);

				ppft.cmass_litter_sap   += clitter_sap;
				ppft.cmass_litter_heart += clitter_heart;

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

	// Nitrogen always needs to be returned.
	// Alive and crop to litter/harvest
	if (alive || istruecrop_or_intercropgrass()) { // DW_COMMENT this if statement is not needed after r11208
		if (pft.lifeform == TREE) {

			double nlitter_sap, nlitter_heart, nwood_harvest;

			// Transfer nitrogen storage to sapwood nitrogen litter/harvest
			partition_wood_biomass(nmass_sap, nmass_heart,
				harv_eff, harvest_slow_frac, res_outtake,
				nlitter_sap, nlitter_heart,
				nwood_harvest, nharvested_products_slow);

			// Don't allow litter C:N below minimum level (needed due to cdebt)
			double nmass_NH4_sap = max(0.0, nlitter_sap - clitter_sap / pft.cton_sap_min);
			ppft.nmass_litter_NH4 += nmass_NH4_sap;
			ppft.nmass_litter_sap += nlitter_sap - nmass_NH4_sap;

			double nmass_NH4_heart = max(0.0, nlitter_heart - clitter_heart / pft.cton_heart_min);
			ppft.nmass_litter_NH4 += nmass_NH4_heart;
			ppft.nmass_litter_heart += nlitter_heart - nmass_NH4_heart;

			nharvest_flux += nwood_harvest;
		}

		// Leaf: remove residue outtake and send the rest to litter
		if (cmass_leaf > 0.0) {
			ppft.nmass_litter_leaf += nmass_leaf * (1.0 - res_outtake);
		}
		else {
			ppft.nmass_litter_NH4 += nmass_leaf * (1.0 - res_outtake);
		}
		nharvest_flux += nmass_leaf * res_outtake;

		// Root: all goes to litter
		if (cmass_root > 0.0) {
			ppft.nmass_litter_root += nmass_root;
		}
		else {
			ppft.nmass_litter_NH4 += nmass_root;
		}

		if (pft.landcover == CROPLAND) {
			if (pft.aboveground_ho) {
				ppft.nmass_litter_leaf += cropindiv->nmass_ho * (1.0 - res_outtake);
				nharvest_flux += cropindiv->nmass_ho * res_outtake;
			}
			else {
				ppft.nmass_litter_root += cropindiv->nmass_ho;
			}

			ppft.nmass_litter_leaf += cropindiv->nmass_agpool * (1.0 - res_outtake);
			nharvest_flux += cropindiv->nmass_agpool * res_outtake;
			ppft.nmass_litter_leaf += cropindiv->nmass_dead_leaf * (1.0 - res_outtake);
			nharvest_flux += cropindiv->nmass_dead_leaf * res_outtake;
		}
	}
	else {
		// Else to NH4 pool
		if (pft.lifeform == TREE) {
			ppft.nmass_litter_NH4 += nmass_sap + nmass_heart;
		}
		
		ppft.nmass_litter_NH4 += nmass_leaf + nmass_root + nstore();

		if (pft.landcover == CROPLAND) {
			ppft.nmass_litter_NH4 += cropindiv->nmass_ho + cropindiv->nmass_agpool + cropindiv->nmass_dead_leaf;
		}
	}

	// Report harvest fluxes
	report_flux(Fluxes::HARVESTC, charvest_flux);
	report_flux(Fluxes::HARVESTN, nharvest_flux);

	// Add to biomass depositories for long-lived products
	ppft.cmass_harvested_products_slow += charvested_products_slow;
	ppft.nmass_harvested_products_slow += nharvested_products_slow;
}

double Individual::wscal_mean() const {
	return patchpft().wscal_mean;
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
		& first_autumndate
		& first_autumndate20
		& first_autumndate_20
		& last_springdate
		& last_springdate20
		& last_springdate_20
		& last_verndate
		& last_verndate20
		& last_verndate_20
		& sdate_default
		& sdatecalc_temp
		& sdatecalc_prec
		& sdate_force
		& hdate_force
		& Nfert_read
		& hlimitdate_default
		& wintertype
		& swindow
		& swindow_irr
		& sowing_restriction;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellst member functions
////////////////////////////////////////////////////////////////////////////////

void Gridcellst::serialize(ArchiveStream& arch) {
	arch & frac
		& frac_old_orig
		& nstands
		& distinterval_st
		& diam_cut_low
		& nfert;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Landcover member functions
////////////////////////////////////////////////////////////////////////////////

Landcover::Landcover() {

	wet_startday = 0;
	wet_days = 0;

	updated = false;

	acflux_harvest_slow = 0.0;
	acflux_landuse_change = 0.0;
	acflux_landuse_change_orig = 0.0;
	acflux_wood_harvest = 0.0;
	acflux_wood_harvest_orig = 0.0;
	acflux_clearing = 0.0;
	acflux_clearing_orig = 0.0;
	cmass_stem_harvest = 0.0;
	cmass_stem_toprod = 0.0;
	cmass_harv_killed = 0.0;
	cmass_harv_tolitter = 0.0;
	anflux_harvest_slow = 0.0;
	anflux_landuse_change = 0.0;
	anflux_landuse_change_orig = 0.0;
	anflux_wood_harvest = 0.0;
	anflux_wood_harvest_orig = 0.0;
	anflux_clearing = 0.0;
	anflux_clearing_orig = 0.0;

	for (int i=0; i<NLANDCOVERTYPES; i++) {

		frac[i] = 0.0;
		frac_old[i] = 0.0;
		frac_change[i] = 0.0;
		acflux_harvest_slow_lc[i] = 0.0;
		acflux_wood_harvest_lc[i] = 0.0;
		acflux_clearing_lc[i] = 0.0;
		acflux_cloned_lc[i] = 0.0;
		acflux_landuse_change_lc[i] = 0.0;
		anflux_harvest_slow_lc[i] = 0.0;
		anflux_landuse_change_lc[i] = 0.0;
		anflux_wood_harvest_lc[i] = 0.0;
		anflux_clearing_lc[i] = 0.0;

		for(int j=0;j<NLANDCOVERTYPES;j++) {
			frac_transfer[i][j] = 0.0;
		}

		expand_to_new_stand[i] = (i == NATURAL || i == FOREST);

		pool_to_all_landcovers[i] = false;		// from a donor landcover; alt.c
		pool_from_all_landcovers[i] = false;	// to a receptor landcover; alt.a
	}
}

void Landcover::serialize(ArchiveStream& arch) {
	arch & frac;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcell member functions
////////////////////////////////////////////////////////////////////////////////

Gridcell::Gridcell():climate(*this) {

	for (unsigned int p=0; p<pftlist.nobj; p++) {
		pft.createobj(pftlist[p]);
	}

	for (unsigned int s=0; s<stlist.nobj; s++) {
		st.createobj(stlist[s]);
	}

	if (!run_landcover) {
		create_stand(NATURAL);
		landcover.frac[NATURAL] = 1.0;
	}
	
	// Initialise SIMFIRE variables
	for (int i = 0; i<AVG_INTERVAL_FAPAR; i++) {
		fapar_recent_max[i] = 0.5;
	}
	fapar_annual_max = 0.5;

	// Initialize Max annual Nesterov Index on first day of simulation
	for (int i = 0; i<12; i++) {
		nesterov_monthly_max[i] = 0.;
	}
	nesterov_cur = 0.;

	// Initialise BLAZE variables
	seed = 12345678;

	distinterval_gc = 1.0e10;
	for (int m=0;m<12;m++) {
		monthly_burned_area[m] = 0.0;
		monthly_fire_risk[m] = 0.0;

	}
	burned_area = 0.0;
	simfire_region = 0;
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

void Gridcell::set_soil_depth() {

	double soil_lyr_depth[NSOILLAYERS_MAX];
	int nsoillayers = 0;
	
	soil_layer_depths(soil_lyr_depth, nsoillayers);

	// Set fixed rootdist per PFT depending on soil depth and layer thickness
	for (unsigned int p = 0; p<pftlist.nobj; p++) {
		pft[p].pft.init_rootdist(soil_lyr_depth, soil_depth, nsoillayers);
	}

	// Create soil layers for each patch
	Gridcell::iterator gc_itr = begin();
	while (gc_itr != end()) {

		// START OF LOOP THROUGH STANDS
		Stand& stand = *gc_itr;

		stand.firstobj();
		while (stand.isobj) {

			// START OF LOOP THROUGH PATCHES

			// Get reference to this patch
			Patch& patch = stand.getobj();

			for (int i = 0; i < npft; i++) {
				for (int lyr = 0; lyr < nsoillayers; lyr++) {
					patch.pft[i].rootdist[lyr] = pft[i].pft.rootdist[lyr];
				}
			}

			// Get reference to this soil column
			Soilcolumn& soilcolumn = patch.soilcolumn;

			// Set soilcolumn depth
			soilcolumn.depth = soil_depth;

			// Set number of soil layer in soilcolumn
			soilcolumn.nsoillayers = nsoillayers;

			// Create soil layers and initalise them
			soilcolumn.create_soillayers(soil_lyr_depth);

			// Set padding layers depth
			soilcolumn.pad.set_depths(soilcolumn.soillayers[nsoillayers - 1].depth);
			
			stand.nextobj();
		}// End of loop through patches

		++gc_itr;
	}	// End of loop through stands
}

Stand& Gridcell::create_stand_lu(StandType& st, double fraction, int no_patch, bool suppress_disturbance) {

	Stand& stand = create_stand(st.landcover, no_patch);
	stand.init_stand_lu(st, fraction, suppress_disturbance);

	return stand;
}

double Gridcell::ccont() {

	double ccont = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ccont += stand.ccont() * stand.get_gridcell_fraction();
	}

	return ccont;
}

double Gridcell::ncont() {

	double ncont = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ncont += stand.ncont() * stand.get_gridcell_fraction();
	}

	return ncont;
}

double Gridcell::water_content() {

	double water_content = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		water_content += stand.water_content() * stand.get_gridcell_fraction();
	}

	return water_content;
}

double Gridcell::cflux() {

	double cflux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		cflux += stand.cflux() * stand.get_gridcell_fraction();
	}

	cflux += landcover.acflux_landuse_change;
	cflux += landcover.acflux_harvest_slow;
	cflux += landcover.acflux_wood_harvest;
	cflux += landcover.acflux_clearing;

	return cflux;
}

double Gridcell::nflux() {

	double nflux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		nflux += stand.nflux() * stand.get_gridcell_fraction();
	}

	nflux += landcover.anflux_landuse_change;
	nflux += landcover.anflux_harvest_slow;
	nflux += landcover.anflux_wood_harvest;
	nflux += landcover.anflux_clearing;

	return nflux;
}

double Gridcell::water_flux() {

	double water_flux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		water_flux += stand.water_flux() * stand.get_gridcell_fraction();
	}

	return water_flux;
}

void Gridcell::soil_layer_depths(double lyr_depths[NSOILLAYERS_MAX], int& nsoillayers) {

	if (soil_depth < SOILDEPTH_EVAP) {
		fail("Soil depth (%g) needs to be thicker than evaporation depth (%g)\n", soil_depth, SOILDEPTH_EVAP);
	}

	if (hydrology == RICHARDS) {

		// Standard soil layer thickness (mm)
		double std_depths[NSOILLAYERS_MAX] = { 100.0, 100.0, 100.0, 200.0, 200.0, 300.0, 500.0, 500.0, 1000.0, 1000.0, 1500.0, 3000.0, 4000.0, 6000.0, 8000.0 };
		
		double std_tot_depth = 0.0;
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) std_tot_depth += std_depths[lyr];

		// Calculate number of soil layers and their individual depths
		if (soil_depth > std_tot_depth) {

			// Deeper soil than total standard soil layer depths

			// Number of layers that will be adjusted
			int nadj_lyr = 6;
		
			double tot_depth = 0.0;
			// Bottom three soil layers will adjust for
			// deeper soil than standard depth
			for (int lyr = 0; lyr < NSOILLAYERS_MAX - nadj_lyr; lyr++) {
				lyr_depths[lyr] = std_depths[lyr];
				tot_depth += std_depths[lyr];
			}

			// Upscaling fraction to remaining soil layers
			double scalar = (soil_depth - tot_depth) / (std_tot_depth - tot_depth);
		
			// Adjust soil layer thickness to be off 100 mm increments
			for (int lyr = NSOILLAYERS_MAX - nadj_lyr; lyr < NSOILLAYERS_MAX; lyr++) {
				lyr_depths[lyr] = floor(std_depths[lyr] * scalar / 100.0) * 100.0;
				tot_depth += lyr_depths[lyr];
			}

			lyr_depths[NSOILLAYERS_MAX - 1] += soil_depth - tot_depth;

			nsoillayers = NSOILLAYERS_MAX;
		}
		else {

			// Shallower soil than total standard soil layer depths
			// Adjust soil layer thickness to be off 100 mm increments

			int lyr = 0;
			double depth = soil_depth;
			double tot_depth = 0.0;
			double remain_depth = depth - tot_depth;

			while (!negligible(remain_depth)) {
				double ref_depth = lyr == 0 ? std_depths[lyr] : std_depths[lyr - 1];
				if (remain_depth >= ref_depth * 2.0 && lyr < NSOILLAYERS_MAX - 1) {
					double lyr_depth = remain_depth < std_depths[lyr] * 2.0 ? remain_depth / 2.0 : std_depths[lyr];
					lyr_depths[lyr] = floor(lyr_depth / 100.0) * 100.0;
					tot_depth += lyr_depths[lyr];
					remain_depth = depth - tot_depth;
				}
				else {
					lyr_depths[lyr] = remain_depth;
					remain_depth = 0.0;
				}
				lyr++;
			}
			nsoillayers = lyr;
		}
	}
	else if (hydrology == MILLER) {
		// 100mm thick layers
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			lyr_depths[lyr] = soil_depth / NSOILLAYERS_MAX;
		}
		nsoillayers = NSOILLAYERS_MAX;
	}
	else if (hydrology == GERTEN) {
		// Upper layer 500mm thick and lower layer 1000mm
		lyr_depths[0] = SOILDEPTH_UPPER;
		lyr_depths[1] = soil_depth - lyr_depths[0];

		nsoillayers = 2;
	}
}

void Gridcell::serialize(ArchiveStream& arch) {
	arch & climate
		& landcover
		& seed
		& balance
		& nesterov_max
		& nesterov_monthly_max
		& nesterov_cur
		& fapar_recent_max;

	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}

		for (unsigned int i = 0; i < st.nobj; i++) {
			arch & st[i];
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

		st.killall();

		for (unsigned int i = 0; i < stlist.nobj; i++) {
			st.createobj(stlist[i]);
			arch & st[i];
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
	return (int) size();
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


////////////////////////////////////////////////////////////////////////////////
// Implementation of MassBalance member functions
////////////////////////////////////////////////////////////////////////////////

void MassBalance::serialize(ArchiveStream& arch) {
	arch & start_year
		& ccont_zero
		& ccont_zero_scaled
		& cflux_zero
		& ncont_zero
		& ncont_zero_scaled
		& nflux_zero
		& water_cont_zero
		& water_flux_zero
		& ccont
		& ncont
		& water_cont
		& cflux
		& nflux
		& water_flux;
}

/// Should be used together with check_indiv()
void MassBalance::init_indiv(Individual& indiv) {

	Patch& patch = indiv.vegetation.patch;
	Stand& stand = patch.stand;
	if (!stand.is_true_crop_stand())
		return;

	double scale = 1.0;
	if (patch.stand.get_gridcell().landcover.updated && (patch.nharv == 0 || date.day == 0))
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

	double ccont = indiv.ccont();
	ccont += patch.ccont(0.0);
	double cflux = patch.cflux();

	if(check_harvest && patch.isharvestday)
		ccont_zero = ccont_zero_scaled;

	if(date.year >= nyear_spinup && !negligible(ccont - ccont_zero + cflux - cflux_zero, -10)) {
		dprintf("\nStand %d Patch %d Indiv %d C balance year %d day %d: %.10f\n", patch.stand.id, patch.id, indiv.id,
			date.year, date.day, ccont - ccont_zero + cflux - cflux_zero);
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

	double ncont = indiv.ncont();
	ncont += patch.ncont(0.0);
	double nflux = patch.nflux();

	if(check_harvest && patch.isharvestday)
		ncont_zero = ncont_zero_scaled;

	if(date.year >= nyear_spinup && !negligible(ncont - ncont_zero + nflux - nflux_zero, -14)) {
		dprintf("\nStand %d Patch %d Indiv %d N balance year %d day %d: %.10f\n", patch.stand.id, patch.id, indiv.id,
			date.year, date.day, ncont - ncont_zero + nflux - nflux_zero);
		dprintf("N pool change: %.14f\n", ncont - ncont_zero);
		dprintf("N flux: %.14f\n\n",  nflux - nflux_zero);
		balance = false;
	}

	return balance;
}

/// Should be preceded by init_indiv()
/** check_harvest must be true if growth_daily() is tested
 *  canopy_exchange() and growth_daily() and functions in between cannot be tested separately
 */
bool MassBalance::check_indiv(Individual& indiv, bool check_harvest) {

	return check_indiv_C(indiv, check_harvest) && check_indiv_N(indiv, check_harvest);
}

/// Should be used together with check_patch() e.g. in framework()
void MassBalance::init_patch(Patch& patch) {

	Stand& stand = patch.stand;
	Gridcell& gridcell = stand.get_gridcell();

	double scale = 1.0;
	if (patch.stand.get_gridcell().landcover.updated && (patch.nharv == 0 || date.day == 0))
		scale = stand.scale_LC_change;

	ccont_zero = patch.ccont();
	ccont_zero_scaled = patch.ccont(scale, true);
	cflux_zero = patch.cflux();

	if (stand.get_gridcell_fraction())
		cflux_zero += gridcell.landcover.acflux_harvest_slow / stand.get_gridcell_fraction();

	ncont_zero = patch.ncont();
	ncont_zero_scaled = patch.ncont(scale, true);
	nflux_zero = patch.nflux();

	if (stand.get_gridcell_fraction())
		nflux_zero += gridcell.landcover.anflux_harvest_slow / stand.get_gridcell_fraction();

	water_cont_zero = patch.water_content();
	water_flux_zero = patch.water_flux();
}

bool MassBalance::check_patch_C(Patch& patch, bool check_harvest) {

	bool balance = true;
	Stand& stand = patch.stand;
	Gridcell& gridcell = stand.get_gridcell();
	double ccont = patch.ccont();
	double cflux = patch.cflux();

	if (stand.get_gridcell_fraction())
		cflux += gridcell.landcover.acflux_harvest_slow / stand.get_gridcell_fraction();

	if (check_harvest && patch.isharvestday)
		ccont_zero = ccont_zero_scaled;

	if (date.year >= nyear_spinup && !negligible(ccont - ccont_zero + cflux - cflux_zero, -10)) {
		dprintf("\nStand %d Patch %d C balance year %d day %d: %.10f\n", patch.stand.id, patch.id, date.year, date.day,
			ccont - ccont_zero + cflux - cflux_zero);
		dprintf("C pool change: %.10f\n", ccont - ccont_zero);
		dprintf("C flux: %.10f\n\n",  cflux - cflux_zero);
		balance = false;
	}

	return balance;
}

bool MassBalance::check_patch_N(Patch& patch, bool check_harvest) {

	bool balance = true;
	Stand& stand = patch.stand;
	//if (!stand.is_true_crop_stand())
	//	return balance;
	Gridcell& gridcell = stand.get_gridcell();
	double ncont = patch.ncont();
	double nflux = patch.nflux();

	if (stand.get_gridcell_fraction())
		nflux += gridcell.landcover.anflux_harvest_slow / stand.get_gridcell_fraction();

	if (check_harvest && patch.isharvestday)
		ncont_zero = ncont_zero_scaled;

	if (date.year >= nyear_spinup && !negligible(ncont - ncont_zero + nflux - nflux_zero, -14)) {
		dprintf("\nStand %d Patch %d N balance year %d day %d: %.14f\n", patch.stand.id, patch.id, date.year, date.day,
			ncont - ncont_zero + nflux - nflux_zero);
		dprintf("N pool change: %.14f\n", ncont - ncont_zero);
		dprintf("N flux: %.14f\n\n", nflux - nflux_zero);
		balance = false;
	}

	return balance;
}

bool MassBalance::check_patch_water(Patch& patch) {

	bool balance = true;

	double water_content = patch.water_content();
	double water_flux = patch.water_flux();

	if (date.year >= nyear_spinup && !negligible(water_content - water_cont_zero + water_flux - water_flux_zero, -8)) {
		dprintf("\nStand %d Patch %d Water balance year %d day %d: %.9f\n", patch.stand.id, patch.id, date.year, date.day, water_content - water_cont_zero + water_flux - water_flux_zero);
		dprintf("Water pool change: %.10f\n", water_content - water_cont_zero);
		dprintf("Water flux: %.10f\n\n", water_flux - water_flux_zero);
		balance = false;
	}

	return balance;
}

/// Should be preceded by init_patch() e.g. i framework()
/** check_harvest must be true if growth_daily() is tested
 *  canopy_exchange() and growth_daily() and functions in between cannot be tested separately
 *  (init_patch() must be before canopy_exchange() and check_patch() after growth_daily()
 */
bool MassBalance::check_patch(Patch& patch, bool check_harvest) {

	return check_patch_C(patch, check_harvest) && check_patch_N(patch, check_harvest) && check_patch_water(patch);
}

void MassBalance::check_year_N(Gridcell& gridcell) {

	double ncont_year = gridcell.ncont();
	double nflux_year = gridcell.nflux();

	if (date.year == start_year) {
		ncont_zero = ncont_year;
	}
	else {

		nflux += nflux_year;

		// Cropland without N-limitation is not balanced in N, fertilisation gives poorer N-balance
		// For natural vegetation or unfertilised N-limited cropland, the check can be much stricter

		// N balance check:
		double epsilon_biomass = 1.0e-9;
		if (!all_fracs_const)
			epsilon_biomass = 50 * INPUT_RESOLUTION;
		if (fabs(ncont_year - ncont + nflux_year) > epsilon_biomass) {
			dprintf("\n(%.2f, %.2f): N balance year %d: %.9f\n", gridcell.get_lon(), gridcell.get_lat(), date.year,
				ncont_year - ncont + nflux_year);
			dprintf("N pool change: %.9f\n", ncont_year - ncont);
			dprintf("N flux: %.9f\n", nflux_year);
		}
	}

	ncont = ncont_year;
}


void MassBalance::check_year_C(Gridcell& gridcell) {

	double ccont_year = gridcell.ccont();
	double cflux_year = gridcell.cflux();

	if (date.year == start_year) {
		ccont_zero = ccont_year;
	}
	else {

		cflux += cflux_year;

		// C balance check:
		double epsilon_biomass = 1.0e-9;
		if (!all_fracs_const)
			epsilon_biomass = 50 * INPUT_RESOLUTION;
		if (fabs(ccont_year - ccont + cflux_year) > epsilon_biomass) {
			dprintf("\n(%.2f, %.2f): C balance year %d: %.10f\n", gridcell.get_lon(), gridcell.get_lat(), date.year,
				ccont_year - ccont + cflux_year);
			dprintf("C pool change: %.5f\n", ccont_year - ccont);
			dprintf("C flux: %.5f\n", cflux_year);
		}
	}

	ccont = ccont_year;
}

void MassBalance::check_year_water(Gridcell& gridcell) {

	double water_cont_year = gridcell.water_content();
	double water_flux_year = gridcell.water_flux();

	if (date.year == start_year) {
		water_cont_zero = water_cont_year;
	}
	else {

		water_flux += water_flux_year;

		// Water balance check:
		if (!negligible(water_cont_year - water_cont + water_flux_year, -3)) {
			dprintf("\n(%.2f, %.2f): Water balance year %d: %.6f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, water_cont_year - water_cont + water_flux_year);
			dprintf("Water pool change: %.7f\n", water_cont_year - water_cont);
			dprintf("Water flux: %.7f\n", water_flux_year);
			dprintf("On year, day: %4d, %4d \n", date.year, date.day);
		}
	}

	water_cont = water_cont_year;
}


void MassBalance::check_year(Gridcell& gridcell) {

	if (date.year < start_year) {
		return;
	}

	check_year_C(gridcell);

	if (ifcentury) {
		check_year_N(gridcell);
	}

	check_year_water(gridcell);

}

void MassBalance::check_period(Gridcell& gridcell) {

	double epsilon_biomass = 1.0e-9;
	if(!all_fracs_const)
		epsilon_biomass = 50 * INPUT_RESOLUTION;

	// C balance check:
	if (fabs(ccont - ccont_zero + cflux) > epsilon_biomass) {
		dprintf("\nWARNING: (%.2f, %.2f): Period C balance: %.10f\n", gridcell.get_lon(), gridcell.get_lat(), ccont - ccont_zero + cflux);
		dprintf("C pool change: %.10f\n", ccont - ccont_zero);
		dprintf("C fluxes: %.10f\n",  cflux);
	}
	// Cropland without N-limitation is not balanced in N, fertilisation gives poorer N-balance
	// For natural vegetation or unfertilised N-limited cropland, the check can be much stricter
	
	// N balance check:
	if (fabs(ncont - ncont_zero + nflux) > epsilon_biomass) {
		dprintf("\nWARNING: (%.2f, %.2f): Period N balance: %.10f\n", gridcell.get_lon(), gridcell.get_lat(), ncont - ncont_zero + nflux);
		dprintf("N pool change: %.10f\n", ncont - ncont_zero);
		dprintf("N fluxes: %.10f\n", nflux);
	}

	// Water balance check:
	if (!negligible(water_cont - water_cont_zero + water_flux, -2)) { // NB! Not too strict as small errors can accumulate in time
		dprintf("\nWARNING: (%.2f, %.2f): Period water balance: %.5f\n", gridcell.get_lon(), gridcell.get_lat(), water_cont - water_cont_zero + water_flux);
		dprintf("Water pool change: %.6f\n", water_cont - water_cont_zero);
		dprintf("Water fluxes: %.6f\n", water_flux);
	}
}

void MassBalance::init(Gridcell& gridcell) {

//	start_year = date.year;
	ccont_zero = gridcell.ccont();
	cflux_zero = gridcell.cflux();
	ncont_zero = gridcell.ncont();
	nflux_zero = gridcell.nflux();
	water_cont_zero = gridcell.water_content();
	water_flux_zero = gridcell.water_flux();
}

void MassBalance::check(Gridcell& gridcell) {

	double ccont = gridcell.ccont();
	double cflux = gridcell.cflux();

	if (!negligible(ccont - ccont_zero + cflux, -5)) {
		dprintf("\n(%.2f, %.2f): C balance year %d: %.5f\n", gridcell.get_lon(), gridcell.get_lat(), date.get_calendar_year(), ccont - ccont_zero + cflux);
		dprintf("C pool change: %.5f\n", ccont - ccont_zero);
		dprintf("C flux: %.5f\n\n",  cflux);
	}

	double ncont = gridcell.ncont();
	double nflux = gridcell.nflux();

	if (!negligible(ncont - ncont_zero + nflux, -5)) {
		dprintf("\n(%.2f, %.2f): N balance year %d: %.5f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, ncont - ncont_zero + nflux);
		dprintf("N pool change: %.5f\n", ncont - ncont_zero);
	}

	double water_content = gridcell.water_content();
	double water_flux = gridcell.water_flux();

	if (!negligible(water_content - water_cont_zero + water_flux, -4)) {
		dprintf("\n(%.2f, %.2f): Water balance year %d: %.6f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, water_content - water_cont_zero + water_flux);
		dprintf("Water pool change: %.7f\n", water_content - water_cont_zero);
		dprintf("Water flux: %.7f\n\n", water_flux);
	}
}

bool ManagementType::pftinselection(const char* name) {

	if(planting_system == "")
		return false;

	return issubstring(selection, name);
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995.
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
