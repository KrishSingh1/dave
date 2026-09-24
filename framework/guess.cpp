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
#include "growth.h"
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

double SOILDEPTH_EVAP = SOILDEPTH_EVAP_TRUNK;

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
// Implementation of PhotosynthesisSubdaily member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisSubdaily::serialize(ArchiveStream& arch) {
	arch & vcmax
		 & jmax
		 & gamstar
		 & Rd
		 & An
		 & gsw
		 & ci;
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
		& daylength_max
		& daylength_min
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
		& aprec
		& rainfall_annual_avg
		& last_rainfall
		& days_since_last_rainfall
		& kbdi
		& ffdi_monthly
		& weathergenstate;
}

double Climate::get_tgrowth() const {
	return dtemp_31.periodicmean(30);
}

double Climate::get_gtemp_mean() const {
	if (date.diurnal()) {
		double sum = 0;
		for (size_t i = 0; i < date.subdaily; i++) {
			sum += gtemps[i];
		}
		return sum / date.subdaily;
	}
	return gtemp;
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
    annual_fluxes_per_pft(npft, std::vector<double>(NPERPFTFLUXTYPES)),
    daily_fluxes_per_pft(size_t(date.year_length()) * size_t(npft) * size_t(NPERPFTFLUXTYPES), 0.0) {

	reset();
}

size_t Fluxes::daily_flux_index(int day, int pft_id, PerPFTFluxType flux_type) const {
	const size_t pft_count = annual_fluxes_per_pft.size();
	return (size_t(day) * pft_count + size_t(pft_id)) * size_t(NPERPFTFLUXTYPES) + size_t(flux_type);
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
	}
	std::fill(daily_fluxes_per_pft.begin(), daily_fluxes_per_pft.end(), 0);
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
	daily_fluxes_per_pft[daily_flux_index(date.day, pft_id, flux_type)] += value;
}

void Fluxes::report_flux(PerPatchFluxType flux_type, double value) {
	monthly_fluxes_patch[date.month][flux_type] += value;
	daily_fluxes_patch[date.day][flux_type] += value;
}

double Fluxes::get_daily_flux(PerPFTFluxType flux_type, int day, int pft_id) const {
	return daily_fluxes_per_pft[daily_flux_index(day, pft_id, flux_type)];
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

double Fluxes::get_flux(PerPatchFluxType flux_type, Period period) const {
	switch (period) {
		case Period::DAY:
			return get_daily_flux(flux_type, date.day);
		case Period::MONTH:
			return get_monthly_flux(flux_type, date.month);
		case Period::YEAR:
			return get_annual_flux(flux_type);
		default:
			fail("Invalid flux period");
	}
}

double Fluxes::get_flux(PerPFTFluxType flux_type, Period period) const {
	switch (period) {
		case Period::DAY:
			return get_daily_flux(flux_type, date.day);
		case Period::MONTH:
			return get_monthly_flux(flux_type, date.month);
		case Period::YEAR:
			return get_annual_flux(flux_type);
		default:
			fail("Invalid flux period");
	}
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
		& wscal_5
		& wscal_10
		& soiltemp_10
		& daylength_10
		& reldaylength
		& ftemp_cum
		& wscal_30
		& airtemp_10
		& anetps_ff_est
		& anetps_ff_est_initial
		& wscal_mean_est
		& phen
		& mphen
		& aphen
		& establish
		& nsapling
		& nsapling_year
		& cmass_litter_leaf
		& cmass_litter_root
		& cmass_litter_sap
		& cmass_litter_heart
		& cmass_litter_repr
		& nmass_litter_repr
		& cmass_ext_seed
		& nmass_ext_seed
		& gcbase
		& gcbase_day
		& wsupply
		& wsupply_leafon
		& wstress
		& wstress_day
		& cmass_harvested_products_slow
		& nmass_litter_leaf
		& nmass_litter_root
		& nmass_litter_sap
		& nmass_litter_heart
		& nmass_litter_repr
		& nmass_harvested_products_slow
		& nmass_litter_NH4
		& swindow
		& wmass_deficit_y
		& inund_count
		& inund_stress;
	if (pft.landcover==CROPLAND)
		arch & *cropphen;

}

/*
Ensure the given wscal value lies within range [0, 1]. If not, terminate the
simulation.
*/
void validate_wscal(double wscal, const char *context, const char *pft) {
	if (wscal < 0 || wscal > 1) {
		fail("%s: wscal value %f for PFT %s out of range [0,1]", context, wscal, pft);
	}
}

void Patchpft::set_wscal(double wstress) {
	if (wscal_set_today) {
		// Debatable what we should do in this case. For now I'm going with an
		// immediate fail. Alternatively, we could update to the newer values,
		// although why you would want to do that is beyond me. If you do that,
		// make sure you take that into account wrt wscal_mean and the moving
		// averages.
		fail("Patchpft::set_wscal(): wscal has already been set today. Don't do this to yourself.");
	}

	if (pft.landcover == CROPLAND) {
		fail("Patchpft::set_wscal(): crops not supported (need to implement growing-season wscal as in canexch.cpp water_scalar())");
	}
	validate_wscal(wstress, "Patchpft::set_wscal()", (const char *)pft.name);

	// This gets reset to false in daily_accounting().
	wscal_set_today = true;

	wscal = wstress;
	wscal_5.add(wstress);
	wscal_10.add(wstress);
	wscal_30.add(wstress);

	// This is reset to 0 in Patchpft::annual_accounting().
	wscal_mean += wstress;
	if (date.islastmonth && date.islastday) {
		wscal_mean /= date.year_length();
	}
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

void Patchpft::daily_accounting(const Patch &patch) {
	if (date.day == 0) {
		annual_accounting(patch);
	}
	double total_lai = 0;
	double total_agb = 0;
	for (size_t i = 0; i < patch.vegetation.nobj; ++i) {
		const Individual& indiv = patch.vegetation[i];
		if (indiv.pft.id == pft.id) {
			total_lai += indiv.lai;
			total_agb += indiv.aboveground_biomass();
		}
	}

	lai_ymax = max(lai_ymax, total_lai);
	agb_ymax = max(agb_ymax, total_agb);
	lai_ymean += total_lai;
	agb_ymean += total_agb;

	if (date.islastday && date.islastmonth) {
		int year_length = date.year_length();
		lai_ymean /= year_length;
		agb_ymean /= year_length;
	}

	airtemp_10.add(patch.get_climate().temp);
	wscal_set_today = false;
	dcmass_mort_kill = 0.0;
	dmort_cstarve_ntree = 0.0;
}

void Patchpft::annual_accounting(const Patch &patch) {
	lai_ymax = 0;
	agb_ymax = 0;
	lai_ymean = 0;
	agb_ymean = 0;
	wscal_mean = 0;
	if (pft.phenology==CROPGREEN || pft.isintercropgrass) {
		cropphen->growingdays_y = 0;
	}
	acmass_mort_kill = 0.0;
	amort_cstarve = 0.0;
	amort_cstarve_ntree = 0.0;
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

			// Standard depth
			const double depth_standard = SOILDEPTH_UPPER + SOILDEPTH_LOWER;

			// Standard root depth
			double depth_root = depth_standard;

			// Fraction of fine roots to be homogenously distributed within root depth
			// is zero for MILLER hydrology
			double fhom_root = 0.0;

			if (hydrology == RICHARDS) {

				// Fraction of deep soil, beyond standard depth, that roots have access to
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
	id(i),stand(s),vegetation(*this),soil(*this, st),fluxes(*this) {

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

	transp = 0.0;
	intercep = 0.0;

	burned = false;
	fire_line_intensity = 0.0;
	avail_fuel = 0;
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
		& soil
		& fluxes
		& fpar_grass
		& fpar_ff
		& par_grass_mean
		& nday_growingseason
		& fpc_total
		& rpc_total
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
		& transp
		& awetland_water_added
		& apet
		& eet_net_veg
		& wdemand
		& wdemand_day
		& wdemand_leafon
		& fpc_rescale
		& ndemand
		& fire_line_intensity
		& avail_fuel
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

const int DAVE_DIAGNOSTIC_PRECISION = 16;

/// C content of patch
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 */
double Patch::ccont(double scale_indiv, bool luc) {

	double ccont = soil.ccont();

	ccont += soil.get_sum(&Soillayer::CH4) * KG_PER_G;

	ccont += soil.snow.ccont();

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
double Patch::ncont(double scale_indiv, bool luc) const {

	double ncont = soil.ncont();

	ncont += soil.snow.ncont();

	for (int i = 0; i<npft; i++) {
		Patchpft& ppft = pft[i];
		ncont += ppft.nmass_litter_leaf;
		ncont += ppft.nmass_litter_root;
		ncont += ppft.nmass_litter_sap;
		ncont += ppft.nmass_litter_heart;
		ncont += ppft.nmass_harvested_products_slow;
		ncont += ppft.nmass_litter_NH4;
		ncont += ppft.nmass_litter_repr;
	}

	for (unsigned int i = 0; i<vegetation.nobj; i++) {

		Individual& indiv = vegetation[i];
		ncont += indiv.ncont(scale_indiv, luc);
	}

	return ncont;
}

/// Water content of patch
double Patch::water_content() const {

	double water_content = soil.water_and_ice_content();

	water_content += soil.snow.water_and_ice_content();

	water_content += soil.ssw.water_and_ice_content();

	return water_content;
}

/// C flux of patch
double Patch::cflux(Period period) const {

	double cflux = 0.0;

	cflux += fluxes.get_flux(Fluxes::LEACHORGC, period);
	cflux += -fluxes.get_flux(Fluxes::NPP, period);
	cflux += fluxes.get_flux(Fluxes::REPRC, period);
	cflux += fluxes.get_flux(Fluxes::SOILC, period);
	cflux += fluxes.get_flux(Fluxes::FIREC, period);
	cflux += fluxes.get_flux(Fluxes::ESTC, period);
	cflux += fluxes.get_flux(Fluxes::SEEDC, period);
	cflux += fluxes.get_flux(Fluxes::MANUREC, period);
	cflux += fluxes.get_flux(Fluxes::HARVESTC, period);
	cflux += fluxes.get_flux(Fluxes::CH4C, period) * KG_PER_G; // Convert to kg CH4-C m-2 from g CH4-C m-2

	return cflux;
}

/// N flux of patch
double Patch::nflux(Period period) const {

	const Gridcell &gridcell = stand.get_gridcell();
	double nflux = 0.0;

	double nh4dep = period == Period::YEAR ? gridcell.aNH4dep : gridcell.dNH4dep;
	double no3dep = period == Period::YEAR ? gridcell.aNO3dep : gridcell.dNO3dep;

	nflux += -(nh4dep + no3dep);
	nflux += -fluxes.get_flux(Fluxes::NFERT, period);
	nflux += -fluxes.get_flux(Fluxes::NFIX, period);

	nflux += fluxes.get_flux(Fluxes::LEACHMINN, period);
	nflux += fluxes.get_flux(Fluxes::LEACHORGN, period);
	nflux += fluxes.get_flux(Fluxes::HARVESTN, period);
	nflux += fluxes.get_flux(Fluxes::SEEDN, period);
	nflux += fluxes.get_flux(Fluxes::NH3_FIRE, period);
	nflux += fluxes.get_flux(Fluxes::NOx_FIRE, period);
	nflux += fluxes.get_flux(Fluxes::N2O_FIRE, period);
	nflux += fluxes.get_flux(Fluxes::N2_FIRE, period);
	nflux += fluxes.get_flux(Fluxes::N2O_SOIL, period);
	nflux += fluxes.get_flux(Fluxes::N2_SOIL, period);
	nflux += fluxes.get_flux(Fluxes::NO_SOIL, period);
	nflux += fluxes.get_flux(Fluxes::NH3_SOIL, period);
	nflux += fluxes.get_flux(Fluxes::ESTN, period);

	return nflux;
}

/// water flux of patch
double Patch::water_flux(Period period) const {

	const Climate &climate = get_climate();
	double precip = period == Period::YEAR ? climate.aprec : climate.prec;

	double wflux = 0.0;

	// Influx
	wflux += -precip;
	wflux += -fluxes.get_flux(Fluxes::WIRRI, period);
	wflux += -fluxes.get_flux(Fluxes::WADDED, period);

	// Outflux
	wflux += fluxes.get_flux(Fluxes::WINTERCEP, period);
	wflux += fluxes.get_flux(Fluxes::WTRAN, period);
	wflux += fluxes.get_flux(Fluxes::WEVAP, period);
	wflux += fluxes.get_flux(Fluxes::WRUNOFF_SURF, period);
	wflux += fluxes.get_flux(Fluxes::WRUNOFF_SSW, period);
	wflux += fluxes.get_flux(Fluxes::WRUNOFF_DRAIN, period);
	wflux += fluxes.get_flux(Fluxes::WRUNOFF_BASE, period);

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

void Standpft::daily_accounting() {
	if (date.day == 0) {
		annual_accounting();
	}
}

void Standpft::annual_accounting() {
	cmass_repr = 0;
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
	animal_units = 0.0;
	days_grazed = 0.0;
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
			|| naturalgrass && pftx.landcover == NATURAL && pftx.is_grass() // Allow natural grass pft:s to grow in e.g. forests
			|| pftx.landcover == landcover && landcover == FOREST && pftx.is_grass()) { // Always allow forest grass pft:s to grow in forests

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

double Stand::water_content() const {

	double water_content = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		water_content += (*this)[p].water_content() / nobj;

	return water_content;
}

double Stand::cflux(Period period) const {

	double cflux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		cflux += (*this)[p].cflux(period) / nobj;

	return cflux;
}

double Stand::nflux(Period period) const {

	double nflux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		nflux += (*this)[p].nflux(period) / nobj;

	return nflux;
}

double Stand::water_flux(Period period) const {

	double water_flux = 0.0;

	for (unsigned int p = 0; p < nobj; p++)
		water_flux += (*this)[p].water_flux(period) / nobj;

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

size_t get_response_window(double dormancy_depth) {
	double mindorm = 0.2;
	double maxdorm = 0.7;
	int minwin = 2;
	int maxwin = 10;
	double slope = (maxwin - minwin) / (maxdorm - mindorm);
	return round(slope * (maxdorm - min(maxdorm,max(mindorm,dormancy_depth))) + minwin);
}

double calc_root_beta(double wscal_avg) {
	// min and max values from range presented in Jackson et al. 1996
	const double wscal_min = 0.5; // TODO: use existing wscal_min parameter (?)
	const double ROOT_BETA_MIN = 0.94;
	const double ROOT_BETA_MAX = 0.98;
	double x = (1.0 - max(wscal_min,wscal_avg)) / (1.0 - wscal_min);
    return ROOT_BETA_MIN + (ROOT_BETA_MAX - ROOT_BETA_MIN) * pow(x,2.0);
}

// Inits root fractions in each soil layer through a shape parameter beta (see Jackson et al., 1996)
void calc_root_distribution(const Soil &soil, double root_beta, double rootdist[NSOILLAYERS_MAX]) {

	double depth = Dz_soil * CM_PER_MM;

	rootdist[0] = 1.0 - pow(root_beta, depth); // init first layer

	double tot = rootdist[0];
	for (int i = 1; i < soil.nsoillayers; i++){
		depth += Dz_soil * CM_PER_MM;
		rootdist[i] = 1.0 - pow(root_beta, depth) - (1.0 - pow(root_beta, depth - Dz_soil * CM_PER_MM));
		tot += rootdist[i];
	}

	// Calibrated the root_beta for each PFT to match rootdist_upper from 'old' (pre LPJG 4.1) parameterisation.
	// Sometimes the rootdist goes below our maximum soildepth. When that happens, put the residual fraction in lowest soil layer
	rootdist[soil.nsoillayers-1] += 1.0 - tot;

	// Set root distribution in remaining unused soil layers to 0.
	for (size_t i = soil.nsoillayers; i < NSOILLAYERS_MAX; i++) {
		rootdist[i] = 0.0;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp              = 0.0;
	agpp			  = 0.0;
	aresp			  = 0.0;
	aet               = 0.0;
	fpc               = 0.0;
	fpc_greenbrown    = 0.0;
	fpc_daily		  = 0.0;
	rpc               = 0.0;
	densindiv         = 0.0;
	cmass_leaf        = 0.0;
	cmass_root        = 0.0;
	cmass_crown       = 0.0;
	cmass_repr        = 0.0;
	cmass_sap         = 0.0;
	cmass_heart       = 0.0;
	cmass_debt        = 0.0;
	cmass_wood_old	  = 0.0;
	cmass_veg         = 0.0;
	cmass_leaf_post_turnover      = 0.0;
	cmass_root_post_turnover      = 0.0;
	cmass_sap_yesterday = 0.0;
	cmass_sap_turnover = 0.0;
	kmax_drought_affected = pft.kmax;
	max_k_loss = 0.0;
	xylem_loss_frac = 0.0;
	closs_embolism = 0;
	nloss_embolism = 0;
	closs_starve = 0;
	nloss_starve = 0;
	closs_hydraulic = 0;
	nloss_hydraulic = 0;
	closs_blaze = 0;
	nloss_blaze = 0;
	closs_globfirm = 0;
	nloss_globfirm = 0;
	closs_mort = 0;
	nloss_mort = 0;
	closs_shade = 0;
	nloss_shade = 0;
	closs_drawdown = 0;
	cmass_tot_luc     = 0.0;
	cmass_leaf_limit  = 0.0;
	cmass_root_limit  = 0.0;
	cmass_sap_limit	  = 0.0;
	cmass_storage_limit = 0.0;
	absolute_alpha_leaf 		  = 0.0;
	absolute_alpha_root 		  = 0.0;
	absolute_alpha_sap  		  = 0.0;
	absolute_alpha_storage        = 0.0;
	cgrow             = 0.0;
	cgrow_leaf	      = 0.0;
	cgrow_root 		  = 0.0;
	cgrow_crown       = 0.0;
	cgrow_sap		  = 0.0;
	cgrow_repr        = 0.0;
	cgrow_storage     = 0.0;
	creturn           = 0.0;
	growth_demand     = pft.is_grass() ? 0.0 : 1.0;
	growth_stage      = 0.0;
	frepr_cstorage    = 0.0;
	fdev_growth       = 1.0;
	cmass_crown_flowering = 0.0;
	fdev_senesc       = 1.0;
	senesc_root       = 0.0;
	growth_stage_days = 0;
	frac_germ         = 0.0;
	cmass_germ        = 0.0;
	acmass_germ       = 0.0;
	phen              = 0.0;
	aphen             = 0.0;
	fphen_activity    = 1.0;
	dphen_activity    = 0.0;
	cexcess           = 0.0;
	nexcess           = 0.0;
	nmass_missing     = 0.0;
	deltafpc          = 0.0;
	lai_ymax					= 0.0;
	fpc_ymax					= 0.0;
	cmass_leaf_ymax				= 0.0;
	nmass_leaf_ymax				= 0.0;
	cmass_root_ymax				= 0.0;
	nmass_root_ymax				= 0.0;
	cmass_storage				= 0.0;
	cmass_storage_dynam			= 0.0;
	ftemp			     		= 0.0;
	fwater			   			= 0.0;
	fstorage					= 0.0;
	fphoto						= 0.0;
	cgrow_pot					= 0.0;
	cstorage_max                = 0.0;
	nmass_storage				= 0.0;
	nmass_uptake_no3            = 0.0;
	nmass_uptake_nh4            = 0.0;
	nretranslocated             = 0.0;
	anretranslocated            = 0.0;
	ntoavail                    = 0.0;
	antoavail                   = 0.0;
	nstorage_max				= 0.0;
	nstorage_drawdown           = 0.0;
	anstorage_drawdown          = 0.0;
	cmass_leaf_ygrowth			= 0.0;
	cmass_root_ygrowth			= 0.0;
	cton_leaf_opt				= 0.0;
	cton_leaf_dopt = (p.cton_leaf_min + p.cton_leaf_max) / 2.0;
	nday_leafon					= 0.0;
	csenesc                     = 0.0;
	g1 = pft.g1;
	dnpp_allocated = true;

	anc_frac                    = 0.0;
	anj_frac                    = 0.0;
	anp_frac                    = 0.0;

	Patchpft &ppft = patchpft();
	double dormancy_depth = ppft.dormancy_depth_historical.mean();
	double wscal_avg = ppft.wscal_historical.mean();
	size_t response_window = get_response_window(dormancy_depth);
	assert(response_window > 0);
	wscal_dynamic.set_capacity(response_window);
	soiltemp_dynamic.set_capacity(response_window);
	dphen_activity_max = pft.dphen_activity_max;  // TODO: OK to initialise like this for both ifadaptive_rootdist cases?
	max_mobil_cstorage = pft.max_mobil_cstorage;
	temp_dormancy = pft.temp_dormancy;

	// Initialise dynamic root_beta and rootdist to the pft-level values.
	if (ifadaptive_rootdist && pft.is_grass() && pft.dave_allocation()) {
		root_beta = calc_root_beta(wscal_avg);
		calc_root_distribution(v.patch.soil, root_beta, rootdist);
	} else { // Take PFT-level default
		root_beta = pft.root_beta;
		for (size_t i = 0; i < NSOILLAYERS_MAX; i++) {
			rootdist[i] = pft.input_rootdist[i];
		}
	}
	acalloc_leaf				= 0.0;
	acalloc_root				= 0.0;
	acalloc_repr				= 0.0;
	acalloc_sap					= 0.0;
	acalloc_storage				= 0.0;
	acalloc_crown				= 0.0;

	leaf_out = false;
	leaf_senesced = false;
	leaf_full_switch_summer = false;
	leaf_growth = false;
	senescence_active = false;
	leaf_out_phen_maxlimit = 0.9;

	// Nitrogen
	nmass_leaf        = 0.0;
	nmass_leaf_luc    = 0.0;
	nmass_root        = 0.0;
	nmass_crown       = 0.0;
	nmass_repr        = 0.0;
	nmass_sap         = 0.0;
	nmass_heart       = 0.0;
	cton_leaf_aopt    = 0.0;
	cton_leaf_aavr    = 0.0;
	cton_status       = 0.0;
	nmass_veg         = 0.0;
	nmass_tot_luc     = 0.0;
	ctolitter         = 0.0;
	ntolitter         = 0.0;

	nextin						= 1.0;
	nstore_longterm				= 0.0;
	nstore_labile				= 0.0;
	ndemand						= 0.0;
	fnuptake					= 1.0;
	anuptake					= 0.0;
	dnuptake                    = 0.0;
	max_n_storage				= 0.0;
	scale_n_storage				= 0.0;

	leafndemand					= 0.0;
	rootndemand					= 0.0;
	crownndemand                = 0.0;
	reprndemand					= 0.0;
	sapndemand					= 0.0;
	storendemand				= 0.0;
	hondemand                   = 0.0;
	avmaxnlim                   = 0.0;
	leaffndemand				= 0.0;
	rootfndemand				= 0.0;
	sapfndemand					= 0.0;
	storefndemand				= 0.0;
	leafndemand_store			= 0.0;
	rootndemand_store			= 0.0;
	ps_sink_lim                 = 1.0;

	n_opt_isabovelim  = false;
	nstress           = false;
	stomata_open	  = true;

	weighted_psi = 0;
	predawn_lwp = 0.0;
	lwp_previous = 0.0;
	lwp.resize(date.subdaily);
	k_loss.resize(date.subdaily);
	transpiration.resize(date.subdaily);
	hydraulic_failure = false;
	mort_age = 0;
	mort_greff = 0;
	mort_cstarve = 0.0;
	mort_dcstarve = 0.0;
	wscal = 0.0;
	wscal_mean_year = 0.0;
	lai_trend_max = 0;

	for (int sl = 0; sl<NSOILLAYERS_MAX; sl++) {
		Rs_r[sl] = 0.0;
		fwuptake[sl] = 0.0;
	}

	for (Day day; day.period != date.subdaily; day.next()) {
		lwp[day.period] = 0;
		k_loss[day.period] = 0;
	}

	size_t nleafpools;
	if (pft.phenology == CROPGREEN) {
		nleafpools = 2; // 0 - alive leaf, 1 - dead leaf (not used yet)
	}
	else if (ifdailygrass && pft.is_grass()) {
		nleafpools = 4; // 0 - growing leaves, 1 - first fully expanded leaves, 2 - second fully expanded leaves, 3 - senescing leaves
	}
	else {
		nleafpools = 1; // 0 - alive leaf (not used yet)
	}
	cmass_leafs.resize(nleafpools);
	nmass_leafs.resize(nleafpools);
	for (size_t i = 0; i < nleafpools; i++) {
		cmass_leafs[i] = 0.0;
		nmass_leafs[i] = 0.0;
	}

	alpha_leaf        = 0.0;
	alpha_root        = 0.0;
	alpha_sap         = 0.0;
	alpha_storage     = 0.0;
	cdeficit_leaf     = 0.0;
	cdeficit_root     = 0.0;
	cdeficit_sap      = 0.0;
	cdeficit_storage  = 0.0;
	cmass_leaf_limit  = 0.0;
	cmass_root_limit  = 0.0;
	cmass_sap_limit   = 0.0;
	cturnover_leaf  = 0.0;
 	cturnover_root  = 0.0;
	cturnover_sap  = 0.0;
	cturnover_repr  = 0.0;
	cturnover_crown  = 0.0;

	nturnover_leaf  = 0.0;
 	nturnover_root  = 0.0;
	nturnover_sap  = 0.0;
	nturnover_repr  = 0.0;
	nturnover_crown  = 0.0;
	
	pturnover_leaf   = 0.0;

	charvest_green = 0.0;
	charvest_brown = 0.0;
	charvest_repr  = 0.0;
	nharvest_green = 0.0;
	nharvest_brown = 0.0;
	nharvest_repr  = 0.0;

	// additional initialisation
	age               = 0.0;
	age_active        = 0.0;
	height            = 0.0;
	diam              = 0.0;
	d_height		  = 0.0;
	d_diam		      = 0.0;
	c_d_diam          = 0.0;
	basalarea    	  = 0.0;
	basalarea_inc	  = 0.0;
	crownarea		  = 0.0;
	latosa 			  = 0.0;
	boleht            = 0.0;
	sa 				  = 0.0;
	sapwood_depth     = 0.0;
	fpar              = 0.0;
	aphen_raingreen   = 0;
	intercep          = 0.0;
	phen_mean         = 0.0;
	wstress           = false;
	lambda            = 0.0;
	aaet              = 0.0;
	ltor              = 0.0;
	lai               = 0.0;
	lai_layer         = 0.0;
	lai_indiv         = 0.0;
	lai_brown_indiv   = 0.0;
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

	dgpp              = 0.0;
	dnpp              = 0.0;
	resp_maintenance  = 0.0;
	resp_growth       = 0.0;
	resp_leaf         = 0.0;
	resp_culm         = 0.0;
	resp_sap          = 0.0;
	resp_root         = 0.0;
	resp_repr         = 0.0;
	resp_crown        = 0.0;
	dgpp_no           = 0.0;
	dgpp_ns           = 0.0;
	agpp_no           = 0.0;
	agpp_ns           = 0.0;
	cropindiv         = NULL;
	last_turnover_day = -1;
	andemand_leaf     = 0.0;
	andemand_root     = 0.0;
	andemand_sap      = 0.0;
	andemand          = 0.0;

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
	nscal_5.add(1.0);
	nscal_30.add(1.0);

	if (ifdavephotosyn) {
		ps_inputs_layerwise.resize(date.subdaily);
		ps_outputs_layerwise.resize(date.subdaily);
		energy_balance_layerwise.resize(date.subdaily);
		ps_diagnostics_layerwise.resize(date.subdaily);
	}
}

void Individual::serialize(ArchiveStream& arch) {
	arch & cmass_leaf
		& cmass_root
		& cmass_crown 
		& cmass_repr
		& cmass_sap
		& cmass_heart
		& cmass_debt
		& ctolitter
		& ntolitter
		& cmass_veg
		& cmass_leaf_post_turnover
		& cmass_root_post_turnover
		& cmass_leaf_limit
		& cmass_root_limit
		& cmass_sap_limit
		& last_turnover_day
		& fpc
		& fpc_greenbrown
		& fpc_daily
		& rpc
		& fpar
		& densindiv
		& cgrow
		& cgrow_leaf
		& cgrow_root
		& cgrow_sap
		& cgrow_repr
		& growth_demand
		& growth_stage
		& phen
		& aphen
		& aphen_raingreen
		& fphen_activity
		& dphen_activity
		& cexcess
		& nmass_missing
		& anpp
		& agpp
		& aresp
		& gpp_15
		& gpptrend_15
		& cmass_leafs
		& nmass_leafs
		& cmass_storage
		& cmass_storage_dynam
		& ftemp
		& fwater
		& fstorage
		& cgrow_pot
		& lai_20
		& cmass_storage_20
		& cexcess_npp_5
		& cstorage_dynam_14
		& cstorage_max
		& nmass_storage
		& nstorage_max
		& nscal_5
		& nscal_30
		& cton_leaf_opt	 
		& aet
		& aaet
		& ltor
		& height
		& diam
		& basalarea
		& crownarea
		& deltafpc
		& boleht
		& lai
		& lai_layer
		& lai_indiv
		& lai_brown_indiv
		& lai_daily
		& lai_indiv_daily
		& lai_trend_max
		& greff_5
		& cmass_wood_inc_5
		& age
		& age_active
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
		& nmass_crown
		& nmass_repr
		& nmass_sap
		& nmass_heart
		& nmass_uptake_no3
		& nmass_uptake_nh4
		& ntoavail
		& antoavail
		& cmass_wood_old
		& nextin
		& nstore_longterm
		& nstore_labile
		& ndemand
		& fnuptake
		& anuptake
		& dnuptake
		& max_n_storage
		& scale_n_storage
		& avmaxnlim
		& cton_leaf_aopt
		& cton_leaf_aavr
		& cton_status
		& nmass_veg
		& leaf_full_switch_summer
		& senescence_active
		& leaf_senesced
		& leaf_growth
		& leaf_out
		& leaf_out_phen_maxlimit
		& photosynthesis
		& nstress
		& leafndemand
		& rootndemand
		& crownndemand
		& reprndemand
		& sapndemand
		& storendemand
		& leaffndemand
		& rootfndemand
		& sapfndemand
		& storefndemand
		& leafndemand_store
		& rootndemand_store
		& nday_leafon
		& ps_sink_lim

		& andemand_leaf
		& andemand_root
		& andemand_sap
		& andemand
		& resp_leaf
		& resp_culm
		& resp_sap
		& resp_root
		& resp_repr
		& resp_crown;

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

double Individual::get_cmass_leaf_brown() const {
	if (ifdailygrass && pft.is_grass()) {
		if (ifdavegrowth) {
			return cmass_leafs[brown];
		} else {
			return cmass_leafs[3];
		}
	}
	return 0;
}

// Receive active cmass leaf 
// Will change name to cmass_leaf() when indiv.cmass_leaf variable it replaced everywhere by cmass_leafs array
double Individual::get_cmass_leaf() const {
	if (pft.phenology == CROPGREEN) {
		return cmass_leaf;// cmass_leafs[0]; // alive leaves
	}
	else if (ifdailygrass && pft.is_grass() && alive) {
		if (ifdavegrowth) {
			return cmass_leafs[0]; // only one green leaf pool in davegrowth
		}
		else {
			return cmass_leafs[0] + cmass_leafs[1] + cmass_leafs[2]; // growing leaves and first and second fully expanded leaves
		}
	}
	else {
		return cmass_leaf;// cmass_leafs[0]; // all cmass_leaf
	}
}

/// Help function for reduce_biomass(), partitions nstore into leafs and roots
// Receive active nmass leaf 
// Will change name to nmass_leaf() when indiv.nmass_leaf variable it replaced everywhere by nmass_leafs array
double Individual::get_nmass_leaf() const {
	if (pft.phenology == CROPGREEN) {
		return nmass_leaf;// nmass_leafs[0]; // alive leaves
	}
	else if (ifdailygrass && pft.is_grass() && alive) {
		if (ifdavegrowth) {
			return nmass_leafs[0];  // only one green leaf pool in davegrowth
		} 
		else {
			return nmass_leafs[0] + nmass_leafs[1] + nmass_leafs[2]; // growing leaves and first and second fully expanded leaves
		}
	}
	else {
		return nmass_leaf;// nmass_leafs[0]; // all cmass_leaf
	}
}

/// Help function for reduce_biomass(), partitions nstore into leaves and roots
/// Stored N is partioned out to make sure tissue C:N ratios are within limits
/**
 *  Help function for reduce_biomass(). Stored N is partioned out to make sure 
 *  tissue C:N ratios are within limits
 */
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

void Individual::reduce_biomass(double mortality, double mortality_fire, double &closs, double &nloss) {

	// This function needs to be modified if a new lifeform is added,
	// specifically to deal with nstore().
	assert(pft.lifeform == TREE || pft.is_grass() || pft.lifeform == MOSS);

	if (!negligible(mortality)) {

		const double mortality_non_fire = mortality - mortality_fire;
		const double mortality_fire_wood = mortality_fire * pft.ag_woody_frac;
		const double mortality_non_fire_wood = mortality - mortality_fire_wood;

		// Transfer killed biomass to litter
		// (above-ground biomass killed by fire enters atmosphere, not litter)

		Patchpft& ppft = patchpft();

		double cmass_leaf_litter = mortality * get_cmass_leaf();
		double cmass_root_litter = mortality * cmass_root;

		if (pft.is_grass() && ifdailygrass && alive && !istruecrop_or_intercropgrass()){
			// add extra daily carbon pools to litter

			if (ifdavegrowth) {
				cmass_leaf_litter += mortality * cmass_leafs[1];
				cmass_root_litter += mortality * cmass_crown;
			} else {
				cmass_leaf_litter += mortality * cmass_leafs[3];
			}
			cmass_root_litter += mortality * cmass_storage;
		}

		if (pft.lifeform == TREE && ifdavegrowth && ifdailytree) {
			cmass_root_litter += mortality * cmass_storage;
		}

		if (pft.dave_allocation()) {
			cmass_root_litter += mortality * cmass_repr;
		}

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
		double nmass_leaf_litter = mortality * get_nmass_leaf();
		double nmass_root_litter = mortality * nmass_root;
		double nmass_sap_litter = mortality * nmass_sap;
		double nmass_heart_litter = mortality * nmass_heart;

		if (pft.is_grass() && ifdailygrass && alive && !istruecrop_or_intercropgrass()) {
			// add extra daily carbon pools to litter

			if (ifdavegrowth) {
				nmass_leaf_litter += mortality * nmass_leafs[1];
				nmass_root_litter += mortality * nmass_crown;
			} else {
				nmass_leaf_litter += mortality * nmass_leafs[3];
			}
			nmass_root_litter += mortality * nmass_storage;
		}

		if (pft.lifeform == TREE && ifdavegrowth && ifdailytree) {
			nmass_root_litter += mortality * nmass_storage;
		}

		if (pft.dave_allocation()) {
			nmass_root_litter += mortality * nmass_repr;
		}

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
		closs = mortality * ccont();
		nloss = mortality * ncont();

		const double remaining = 1.0 - mortality;

		if (!pft.is_grass() && pft.lifeform != MOSS) {
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

		if (ifdailygrass && pft.is_grass()) {
			// Leafs
			cmass_leafs[0] *= remaining;	// growing leaves
			nmass_leafs[0] *= remaining;
			cmass_leafs[1] *= remaining;	// first fully expanded leaves
			nmass_leafs[1] *= remaining;
			cmass_leafs[2] *= remaining;	// second fully expanded leaves
			nmass_leafs[2] *= remaining;
			cmass_leafs[3] *= remaining;	// senescing leaves
			nmass_leafs[3] *= remaining;
			
			// Roots
			cmass_storage *= remaining;
			nmass_storage *= remaining;
		}

		if (pft.dave_allocation()) {
			cmass_repr *= remaining;
			nmass_repr *= remaining;

			// Storage pool for grasses was handled above.
			if (pft.lifeform == TREE) {
				cmass_storage *= remaining;
				nmass_storage *= remaining;
			} else {
				cmass_crown   *= remaining;
				nmass_crown   *= remaining;
			}
		}
	} else {
		// Negligible mortality.
		closs = 0;
		nloss = 0;
	}
}

double Individual::cton_leaf(bool use_phen /* = true*/) const {

	Stand& stand = vegetation.patch.stand;

	if (stand.is_true_crop_stand() && !negligible(cmass_leaf_today()) && !negligible(get_nmass_leaf())) {
		return cmass_leaf_today() / get_nmass_leaf();
	}
	else if (!stand.is_true_crop_stand() && !negligible(get_cmass_leaf()) && !negligible(get_nmass_leaf())) {
		if (use_phen && !ifdailygrass) {
			if (!negligible(phen)) {
				return cmass_leaf_today() / get_nmass_leaf();
			}
			else {
				return pft.cton_leaf_avr;
			}
		}
		else {
			return get_cmass_leaf() / get_nmass_leaf();
		}
	}
	else {
		return pft.cton_leaf_max;
	}
}

double Individual::cton_root(bool use_phen /* = true*/) const {

	if (!negligible(cmass_root) && !negligible(nmass_root)) {
		if (pft.dave_allocation()) {
			return cmass_root / nmass_root;
		}
		if (use_phen) {
			if (!negligible(cmass_root_today())) {
				return max(pft.cton_root_min, cmass_root_today() / nmass_root);
			}
			else {
				return max(pft.cton_root_avr * pft.cton_leaf_min / pft.cton_leaf_avr, cmass_root / nmass_root);
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

double Individual::cton_crown() const {
	if (pft.is_grass() && ifdailygrass && ifdavegrowth) {
		if (!negligible(cmass_crown) && !negligible(nmass_crown)) {
			return (cmass_crown / nmass_crown);
		} 
		else {
			return pft.cton_root_max;
		}
	} 
	else {
		return 1.0;
	}
}



double Individual::cton_sap() const {

	if (pft.lifeform == TREE) {
		if (!negligible(cmass_sap) && !negligible(nmass_sap)) {
			if (ifdavegrowth) {
				return (cmass_sap / nmass_sap); // stoichiometry should be checked in growth
			} else {
				return max(pft.cton_sap_min, cmass_sap / nmass_sap);
			}
		} else {
			return pft.cton_sap_max;
		}
	} else {
		return 1.0;
	}
}

double Individual::cton_heart() const {

	if (pft.lifeform == TREE) {
		if (pft.dave_allocation()) {
			return cmass_heart / nmass_heart; // stoichiometry should be checked in growth
		} else if (!negligible(cmass_heart) && !negligible(nmass_heart))
			return max(pft.cton_heart_min, cmass_heart / nmass_heart);
		else
			return pft.cton_heart_max;
	} else {
		return 1.0;
	}
}

double Individual::sla() const {
 
	// Prescribed PFT specific SLA
	double current_sla = pft.sla;

	// Age dependent SLA for daily grass
	if (ifdailygrass && pft.is_grass() && !negligible(get_cmass_leaf())) {
		if (ifdavegrowth) {
			current_sla *= 1.0;
		} else {
			current_sla *= (cmass_leafs[0] * 1.1 + cmass_leafs[1] + cmass_leafs[2] / 1.1) / get_cmass_leaf();
		}
	}

	return current_sla;
}

double Individual::aboveground_biomass() const {
	double ccont = get_cmass_leaf();

	if (pft.is_grass()) {
		ccont += cmass_crown; // TODO: should this be included?
	} else if (pft.lifeform != MOSS) {
		ccont += cmass_wood() * pft.ag_woody_frac;
	}

	if (pft.dave_allocation()) {
		ccont += cmass_storage; // TODO: should this be included?
	}

	// Convert kgC/m2 to kg/m2 biomass.
	return 2.0 * ccont;
}

/// C content of individual
/** INPUT PARAMETERS
 *
 *  \param scale_indiv  		scaling factor for living C
 *  \param luc 					down-scales living C (used in C balance tests)
 *  \param ignore_alive			if true, ignore whether individual is alive or not
 */
double Individual::ccont(double scale_indiv, bool luc, bool ignore_alive) const {

	double ccont = 0.0;

	if (ignore_alive || alive || istruecrop_or_intercropgrass()) {

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
			ccont += (get_cmass_leaf() + cmass_root + cmass_sap + cmass_heart - cmass_debt) * scale_indiv;

			if (ifdailygrass && pft.is_grass()) {
				if (ifdavegrowth) {
					ccont += (cmass_storage + cmass_leafs[1] + cmass_crown + cmass_repr) * scale_indiv;
				} else {
					ccont += (cmass_storage + cmass_leafs[3]) * scale_indiv;
				}
			} else if (ifdavegrowth) {
				if (pft.dave_allocation()) {
					ccont += (cmass_storage + cmass_repr) * scale_indiv;
				} else if (!(date.islastday && date.islastmonth)){
					// This individual's carbon is allocated annually, so in order
					// to maintain the carbon balance throughout the year, we must
					// allow for the carbon which has been assimilated but not yet
					// allocated (ie anpp).
					//
					// This is a hack mainly added for debugging purposes and should
					// probably be removed in the fullness of time.
					ccont += anpp * scale_indiv;
				}

			}

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
 *  \param ignore_alive			if true, ignore whether individual is alive or not
 */
double Individual::ncont(double scale_indiv, bool luc, bool ignore_alive) const {

	double ncont = 0.0;

	if (alive || !(date.islastday && date.islastmonth)) {

		if (luc) {

			ncont += nmass_leaf - nmass_leaf_luc * (1.0 - scale_indiv);
			ncont += nmass_root - nmass_root_luc * (1.0 - scale_indiv);
			ncont += nmass_sap - nmass_sap_luc * (1.0 - scale_indiv);
			ncont += nmass_heart - nmass_heart_luc * (1.0 - scale_indiv);
			ncont += nstore_longterm - nstore_longterm_luc * (1.0 - scale_indiv);
			ncont += nstore_labile - nstore_labile_luc * (1.0 - scale_indiv);
			ncont += nmass_storage - nmass_storage_luc * (1.0 - scale_indiv);
			ncont += (nmass_uptake_no3 + nmass_uptake_nh4) * (1.0 - scale_indiv);

			if (ifdailygrass && pft.is_grass()) {
				ncont += nmass_leafs[0] - nmass_leafs_luc[0] * (1.0 - scale_indiv);
				ncont += nmass_leafs[1] - nmass_leafs_luc[1] * (1.0 - scale_indiv);
				ncont += nmass_leafs[2] - nmass_leafs_luc[2] * (1.0 - scale_indiv);
				ncont += nmass_leafs[3] - nmass_leafs_luc[3] * (1.0 - scale_indiv);
			}
		}
		else {

			ncont += get_nmass_leaf() * scale_indiv;
			ncont += nmass_root * scale_indiv;
			ncont += nmass_sap * scale_indiv;
			ncont += nmass_heart * scale_indiv;
			ncont += nstore_longterm * scale_indiv;
			ncont += nstore_labile * scale_indiv;
			ncont += nmass_storage * scale_indiv;
			if (!pft.store_uptake()) {
				ncont += (nmass_uptake_no3 + nmass_uptake_nh4) * scale_indiv;
			}

			if (ifdailygrass && pft.is_grass()) {
				if (ifdavegrowth) {
					ncont += (nmass_leafs[1] + nmass_crown + nmass_repr) * scale_indiv;
				} else {
					ncont += nmass_leafs[3] * scale_indiv;
				}
			}

			if (ifdailytree && pft.lifeform == TREE && ifdavegrowth) {
				ncont += nmass_repr * scale_indiv;
			}
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
	else if (ifdailygrass && pft.is_grass() && !istruecrop_or_intercropgrass()) {
		double froot = 1.0 / (1.0 + pft.ltor_max);
		storendemand = max(0.0, cmass_storage / (cton_leaf_opt * (1.0 - froot) + (cton_leaf_opt * pft.cton_root_avr / pft.cton_leaf_avr) * froot) - nmass_storage);
	}
	else
		storendemand = max(0.0, min(max(0.0, agpp_no) * scale_n_storage, max_n_storage) - nstore());

	return storendemand;
}

bool Pft::dave_allocation() const {
	if (!ifdavegrowth) {
		// Dave growth is disabled. This may be a DW daily grass but either way,
		// it's certainly not using dave allocation.
		return false;
	}

	if (is_grass() && ifdailygrass) {
		// This is a dave daily grass.
		return true;
	}

	if (lifeform == TREE && ifdailytree) {
		// This is a dave daily tree.
		return true;
	}

	// TODO: mosses??
	return false;
}

bool Pft::daily_establishment() const {
	if (is_grass() && ifdailygrass && ifdavegrowth) {
		// DAVE daily grass.
		return true;
	}

	// Everything else uses annual establishment.
	return false;
}

// b is s1
// c is s2

// returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
double Pft::proportion_of_conductivity(double psi) const {
  return exp(-pow((psi / s1), s2));
}

// set spline for proportion of conductivity
void Pft::setup_transpiration(double resolution) {
    // Integrate and accumulate results.
    auto x_psi = std::vector<double>{ 0.0 };
    auto y_cumulative_transpiration = std::vector<double>{ 0.0 };
    double interval = (s1 * pow((log(1 / 0.01)), (1 / s2)));
    double step = interval / (resolution - 1);
	double ymax = -1e36;

	// TODO: repeated additions here will cause rounding errors to accumulate.
	// Need to change to the commented version below, but this affects the
	// results quite a bit. Need to discuss with Clare.

    // for (size_t i = 1; i < resolution; i++) {
		// double psi_spline = step * i;
	for (double psi_spline = 0.0 + step; psi_spline <= interval; psi_spline += step) {
        double e_psi = step * ((proportion_of_conductivity(psi_spline - step) + proportion_of_conductivity(psi_spline)) / 2) + y_cumulative_transpiration.back();

        // x values for spline.
        x_psi.push_back(psi_spline);

        // y values for spline.
        y_cumulative_transpiration.push_back(e_psi);
		ymax = max(ymax, e_psi);
    }

    // Setup the interpolator.
    // We could do a reverse lookup if need, but it's not required.

    transpiration_from_psi.init(x_psi, y_cumulative_transpiration);
    transpiration_from_psi.set_extrapolate(false);

	psi_from_transpiration.init(y_cumulative_transpiration, x_psi);
	psi_from_transpiration.set_extrapolate(false);

	// Note: interpolation will be significantly faster if the x-values are
	// equidistant. Currently, x_psi is equidistant, but cum. transpiration (due
	// to the nonlinearity of the relationship) is not. Therefore, we perform a
	// one-time interpolation at regular intervals along the curve, and use
	// that as the basis for future interpolations during the simulation.
	std::vector<double> y{ 0.0 };
	std::vector<double> x{ 0.0 };
	for (size_t i = 1; i < resolution; i++) {

		// The i-th evenly-spaced point between 0 and ymax.
		double yi = 1.0 * i / (resolution - 1) * ymax;

		// Get the psi corresponding to this y-value.
		double xi = psi_from_transpiration.eval(yi);

        // x values for spline.
        x.push_back(xi);

        // y values for spline.
        y.push_back(yi);
    }

	// Re-initialise psi_from_transpiration.
    psi_from_transpiration.init(y, x);
    psi_from_transpiration.set_extrapolate(false);
}


bool Pft::get_psi_from_transpiration(double transpiration, double &psi) const {
	return psi_from_transpiration.try_eval(transpiration, psi);
}

bool Pft::get_transpiration_lower_bound(double psi, double &transpiration) const {
	return transpiration_from_psi.try_eval(psi, transpiration);
}

/// Checks C mass and zeroes any negative value, balancing by adding to npp and reducing respiration
double Individual::check_C_mass() {

	double negative_cmass = 0.0;

	if (ifdailygrass && pft.is_grass()) {

		if (cmass_storage < 0.0) {
			negative_cmass -= cmass_storage;
			cmass_storage = 0.0;
		}
		if (cmass_root < 0.0) {
			negative_cmass -= cmass_root;
			cmass_root = 0.0;
		}
		if (cmass_leafs[0] < 0.0) {
			negative_cmass -= cmass_leafs[0];
			cmass_leafs[0] = 0.0;
		}
		if (cmass_leafs[1] < 0.0) {
			negative_cmass -= cmass_leafs[1];
			cmass_leafs[1] = 0.0;
		}
		if (cmass_leafs[2] < 0.0) {
			negative_cmass -= cmass_leafs[2];
			cmass_leafs[2] = 0.0;
		}
		if (cmass_leafs[3] < 0.0) {
			negative_cmass -= cmass_leafs[3];
			cmass_leafs[3] = 0.0;
		}

		if (largerthanzero(negative_cmass, -14)) {
			anpp += negative_cmass;
			report_flux(Fluxes::NPP, negative_cmass);
			report_flux(Fluxes::RA, -negative_cmass);
		}

		return negative_cmass;
	}

	if (pft.landcover != CROPLAND)
		return 0;

	negative_cmass = 0.0;

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
			nstore_labile -= negative_nmass * nstore_labile / pos_nmass;
			nstore_longterm -= negative_nmass * nstore_longterm / pos_nmass;
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
			nstore_labile = 0.0;
			nstore_longterm = 0.0;
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
	nmass_storage_luc = nmass_storage;

	if (ifdailygrass && pft.is_grass()) {
		nmass_leafs_luc[0] = nmass_leafs[0];
		nmass_leafs_luc[1] = nmass_leafs[1];
		nmass_leafs_luc[2] = nmass_leafs[2];
		nmass_leafs_luc[3] = nmass_leafs[3];
	}

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
	else if (ifdailygrass && pft.is_grass())
		return get_cmass_leaf();
	else
		return cmass_leaf * phen;
}

/// Gets the individual's daily cmass_root value
double Individual::cmass_root_today() const {

	if (istruecrop_or_intercropgrass())
		return patchpft().cropphen->growingseason ? cropindiv->grs_cmass_root : 0;
	else if (ifdailygrass && pft.is_grass())
		return cmass_root;
	else
		return cmass_root * phen;
}

double Individual::get_cmass_crown() const {
	if (pft.is_grass() && pft.dave_allocation()) {
		return cmass_crown;
	} else {
		return 0.0;
	}
}

double Individual::get_cmass_repr() const {
	if (pft.dave_allocation()) {
		return cmass_repr;
	} else {
		return 0.0;
	}
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
			return crownarea * densindiv * (1.0 - lambertbeer(lai_indiv_today()));
		}
		else {
			// Grass and mosses
			return 1.0 - lambertbeer(lai_indiv_today());
		}
	}
}

/// Gets the individual's daily fpc value including green and brown leaves
double Individual::fpc_greenbrown_today() const {

	// TODO: adjust for vegtypes other than grasses, if needed
	if (pft.phenology == CROPGREEN)
		return patchpft().cropphen->growingseason ? fpc_daily : 0;
	else if (pft.is_grass() && ifdailygrass)
		return fpc_greenbrown;
	else
		return fpc * phen;  
}

/// Gets the individual's daily lai value
double Individual::lai_today() const {

	if (pft.dave_allocation()) {
		return lai;
	} else if (pft.phenology == CROPGREEN)
		return patchpft().cropphen->growingseason ? lai_daily : 0;
	else if (pft.is_grass() && ifdailygrass)
		return lai;
	else
		return lai * phen;
}

/// Gets the individual's daily lai_indiv value
double Individual::lai_indiv_today() const {

	if (pft.phenology == CROPGREEN)
		return patchpft().cropphen->growingseason ? lai_indiv_daily : 0;
	else if (pft.is_grass() && ifdailygrass)
		return lai_indiv;
	else
		return lai_indiv * phen;
}

/// Gets the Nitrigen limited LAI
double Individual::lai_nitrogen_today() const{
	if (pft.phenology==CROPGREEN) {

		double Ln = 0.0;
		if (patchpft().cropphen->growingseason && cmass_leaf_today() > 0.0) {
			const double k = 0.5;
			const double ktn = 0.52 * k + 0.01; // Yin et al 2003
			double nb = 1.0 / (pft.cton_leaf_max * sla());
			Ln = (1.0 / ktn) * log(1.0 + ktn * nmass_leaf / nb);
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
		dprintf("Year %d %s IID %d Diff %g Before %g After %g\n", date.year, (char*)pft.name, id, nbefore - nafter, nbefore, nafter);
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
	flush_unallocated_npp();
	Patchpft& ppft = patchpft();
	double charvest_flux = 0.0;
	double charvested_products_slow = 0.0;

	double nharvest_flux = 0.0;
	double nharvested_products_slow = 0.0;

	double harv_eff = 0.0;
	double harvest_slow_frac = 0.0;
	double res_outtake = 0.0;

	if (alive || istruecrop_or_intercropgrass()) {
		double cmass_mort_kill = ccont(1.0, false, true);
		ppft.dcmass_mort_kill += cmass_mort_kill;
		ppft.acmass_mort_kill += cmass_mort_kill;
	}

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

	if (pft.is_grass() && ifdailygrass && alive && !istruecrop_or_intercropgrass()){
		ppft.cmass_litter_root += cmass_storage;	// storage growth to litter
		if (ifdavegrowth) {
			ppft.cmass_litter_leaf += cmass_leafs[brown];
			flush_repr_litter_DAVE(1.0);
		} else {
			ppft.cmass_litter_leaf += cmass_leafs[3];	// add senescing leaves and growth storage to litter
		}
	}

	if (pft.lifeform == TREE && ifdailytree && ifdavegrowth) {
		ppft.cmass_litter_root += cmass_storage;
	}

	if (!pft.store_uptake()) {
		Soil &soil = vegetation.patch.soil;
		if (nmass_uptake_no3 > 0) {
			soil.soillayers[0].nmass_inc(nmass_uptake_no3, NO3);
			nmass_uptake_no3 = 0;
		}
		if (nmass_uptake_nh4 > 0) {
			soil.soillayers[0].nmass_inc(nmass_uptake_nh4, NH4);
			nmass_uptake_nh4 = 0;
		}
	}

	// C doesn't return to litter/harvest if the Individual isn't alive
	if (alive || istruecrop_or_intercropgrass()) {

		// For leaf and root, catches small, negative values too

		// Leaf: remove residue outtake and send the rest to litter
		if (has_daily_turnover() && cropindiv) {

			if (pft.is_grass() && pft.phenology != CROPGREEN) {
				charvest_flux += cropindiv->grs_cmass_leaf * harv_eff;
				cropindiv->grs_cmass_leaf *= (1 - harv_eff);
			}

			ppft.cmass_litter_leaf += cropindiv->grs_cmass_leaf * (1 - res_outtake);
			charvest_flux    += cropindiv->grs_cmass_leaf * res_outtake;
		}
		else {

			if (pft.is_grass() && pft.phenology != CROPGREEN && !ifdailygrass) {
				charvest_flux += cmass_leaf * harv_eff;
				cmass_leaf *= (1 - harv_eff);
			}
			ppft.cmass_litter_leaf += get_cmass_leaf() * (1 - res_outtake);
			charvest_flux    += get_cmass_leaf() * res_outtake;
		}
		// Root: all goes to litter
		if (has_daily_turnover() && cropindiv)
			ppft.cmass_litter_root += cropindiv->grs_cmass_root;
		else {
			ppft.cmass_litter_root += cmass_root;
			if (pft.is_grass() && ifdavegrowth && ifdailygrass) {
				ppft.cmass_litter_root += cmass_crown;
			}
		}

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

	// Nitrogen always return to soil litter
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
	else {
		// Transfer nitrogen storage to root nitrogen litter
		ppft.nmass_litter_root += nstore();
	}

	if (pft.is_grass() && ifdailygrass && alive && !istruecrop_or_intercropgrass()) {
		ppft.nmass_litter_root += nmass_storage;	// storage growth to litter
		if (ifdavegrowth) {
			ppft.nmass_litter_leaf += nmass_leafs[brown];
			ppft.nmass_litter_root += nmass_crown;
		} else {
			ppft.nmass_litter_leaf += nmass_leafs[3];	// add senescing leaves and growth storage to litter
		}
	}

	// Account for N storage in dave trees.
	if (pft.lifeform == TREE && ifdailytree && ifdavegrowth) {
		ppft.nmass_litter_root += nmass_storage;
	}

	// Leaf: remove residue outtake and send the rest to litter
	ppft.nmass_litter_leaf += get_nmass_leaf() * (1.0 - res_outtake);
	nharvest_flux          += get_nmass_leaf() * res_outtake;

	// Root: all goes to litter
	if (cmass_root > 0.0) {
		ppft.nmass_litter_root += nmass_root;
	}
	else {
		ppft.nmass_litter_NH4 += nmass_root;
	}

	if (pft.landcover == CROPLAND) {
		if (pft.aboveground_ho) {
			ppft.nmass_litter_leaf += cropindiv->nmass_ho * (1 - res_outtake);
			nharvest_flux += cropindiv->nmass_ho * res_outtake;
		}
		else {
			ppft.nmass_litter_root += cropindiv->nmass_ho;
		}

		ppft.nmass_litter_leaf += cropindiv->nmass_agpool * (1 - res_outtake);
		nharvest_flux += cropindiv->nmass_agpool * res_outtake;
		ppft.nmass_litter_leaf += cropindiv->nmass_dead_leaf * (1 - res_outtake);
		nharvest_flux += cropindiv->nmass_dead_leaf * res_outtake;
	}

	if (pft.dave_allocation()) {
		if (pft.is_grass()) {
			ppft.cmass_litter_repr += cmass_repr;
			ppft.nmass_litter_repr += nmass_repr;
		} else {
			// Tree reproductive tissue to be treated as leaf tissue, so that it
			// decays over time. (cmass_litter_repr decay is defined by loss of
			// seed viability, and is much slower than the turnover of the leaf/
			// root/woody "standing dead"/litter turnover).
			ppft.cmass_litter_leaf += cmass_repr;
			ppft.nmass_litter_leaf += nmass_repr;
		}
		cmass_repr = 0;
		nmass_repr = 0;
	}

	// Report harvest fluxes
	report_flux(Fluxes::HARVESTC, charvest_flux);
	report_flux(Fluxes::HARVESTN, nharvest_flux);

	// Add to biomass depositories for long-lived products
	ppft.cmass_harvested_products_slow += charvested_products_slow;
	ppft.nmass_harvested_products_slow += nharvested_products_slow;
}

// Flush reproductive carbon and nitrogen to ppft-level reproductive litter (=seedbank).
// Happens every day (flush_fraction = daily turnover) as well as at the end of 
// the reproductive period or after mortality (flush_fraction - 1.0).
void Individual::flush_repr_litter_DAVE(double flush_fraction) {

	Patchpft &ppft = patchpft();

	if (pft.lifeform == TREE) {
		return;
	}

	// C and N removed from the repr. pool 
	double cflush_removed = flush_fraction * cmass_repr;
	double nflush_removed = flush_fraction * nmass_repr;

	// C and N going to litter or seedbank (can be higher than the one removed from the repr.pool because of initialc and initialn)
	double cflush_export = cflush_removed;
	double nflush_export = nflush_removed;

	// only a fraction of cmass_export (actual seeds) goes to the seedbank.
	// the rest (e.g. panicles, husks) goes to the litter pool.
	double cflux_seed = pft.frac_seed * cflush_export;
	double nflux_seed = pft.frac_seed * nflush_export;

	double cflux_litter = cflush_export - cflux_seed;
	double nflux_litter = nflush_export - nflux_seed; 

	ppft.cmass_litter_repr += cflux_seed;
	ppft.nmass_litter_repr += nflux_seed;

	ppft.cmass_litter_leaf += cflux_litter;
	ppft.nmass_litter_leaf += nflux_litter;			

	assert(equal(cflush_export, cflux_seed + cflux_litter));
	assert(equal(nflush_export, nflux_seed + nflux_litter));

	ctolitter += cflush_export;
	ntolitter += nflush_export;

	cturnover_repr = cflush_export;
	nturnover_repr = nflush_export;

	cmass_repr -= cflush_removed;
	nmass_repr -= nflush_removed;
}

void Individual::daily_accounting() {
	if (date.day == 0) {
		annual_accounting();
	}
	cmass_germ = 0;
	mort_dcstarve = 0.0;
	closs_temp = 0;
	nloss_temp = 0;
	loss_temp = 0;
	xylem_loss_frac = 0.0;
	closs_embolism = 0;
	nloss_embolism = 0;
	closs_starve = 0;
	nloss_starve = 0;
	closs_hydraulic = 0;
	nloss_hydraulic = 0;
	closs_blaze = 0;
	nloss_blaze = 0;
	closs_globfirm = 0;
	nloss_globfirm = 0;
	closs_mort = 0;
	nloss_mort = 0;
	closs_shade = 0;
	nloss_shade = 0;
	closs_drawdown = 0;
	dnpp_allocated = true;
	canopy.clear();

	if (ifdavephotosyn) {
		for (size_t i = 0; i < date.subdaily; i++) {
			ps_inputs_layerwise[i].clear();
			ps_outputs_layerwise[i].clear();
			energy_balance_layerwise[i].clear();
			ps_diagnostics_layerwise[i].clear();
		}
		ps_subdaily.clear();
		photosynthesis.clear();
	}
}

void Individual::annual_accounting() {
	acmass_germ = 0;
	wscal_mean_year = 0;
	mort_age = 0;
	mort_greff = 0;
	mort_cstarve = 0;
	antoavail = 0;
	mort_cstarve = 0.0;
	aaet = 0.0;
	acloss_temp = 0;
	anloss_temp = 0;
	aloss_temp = 0;

	const Gridcell &gridcell = vegetation.patch.stand.get_gridcell();
	if (gridcell.landcover.updated && pft.dave_allocation()) {
		// Rescale indiv following landcover change.
		scale_indiv(*this, false);
	}
}

size_t Individual::get_wscal_dynamic_capacity() const {
	return wscal_dynamic.get_capacity();
}

void Individual::set_wscal_dynamic_capacity(size_t capacity) {
	wscal_dynamic.set_capacity(capacity);
}

double Individual::get_wscal_moving_average() const {
	return wscal_dynamic.size() > 0 ? wscal_dynamic.mean() : 1.0;
}

void Individual::set_wscal(double wscal) {
	validate_wscal(wscal, "Individual::set_wscal()", (const char *)pft.name);

	this->wscal = wscal;
	wscal_dynamic.add(wscal);
	// TODO: combine the 5-, 10-, and 30-day moving averages into a single
	// 30-day moving average, and have get_wscal_N() return the mean of the most
	// recent N days' values.
	wscal_5.add(wscal);
	wscal_10.add(wscal);
	wscal_30.add(wscal);
	wscal_mean_year += wscal;
	if (date.islastmonth && date.islastday) {
		wscal_mean_year /= Date::MAX_YEAR_LENGTH;
	}
}

double Individual::wscal_mean() const {
	if (use_indiv_wscal()) {
		return wscal_mean_year;
	} else {
		return patchpft().wscal_mean;
	}
}

bool Individual::use_indiv_wscal() const {
	// When using the zhou model, individual-level water stress will be
	// different to ppft-level water stress, because it's weighted by root
	// resistance (an indiv-level property).
	//
	// When using adaptive root distribution for dave grasses, individual-level
	// water stress will be different to ppft-level water stress, because the
	// individuals' root distributions can vary over time.
	//
	// When using dave allocation more generally, individual-level water stress
	// will be the same as ppft-level water stress.
	//
	// When not using dave allocation (or the zhou model, which requires dave
	// allocation), individual-level water stress will not be set, and we must
	// therefore fall back to using the ppft-level water stress.
	return ifzhou || pft.dave_allocation();
}

void Individual::flush_unallocated_npp() {
	if (dnpp_allocated) {
		return;
	}
	if (!pft.dave_allocation()) {
		return;
	}
	if (dnpp <= 0.0) {
		return;
	}

	cmass_storage += dnpp;
	dnpp_allocated = true;
}

double Individual::get_wscal() const {
	if (use_indiv_wscal()) {
		return wscal;
	} else {
		return patchpft().wscal;
	}
}

const Historic<double,5> &Individual::get_wscal_5() const {
	if (use_indiv_wscal()) {
		return wscal_5;
	} else {
		return patchpft().wscal_5;
	}
}

const Historic<double,10> &Individual::get_wscal_10() const {
	if (use_indiv_wscal()) {
		return wscal_10;
	} else {
		return patchpft().wscal_10;
	}
}

const Historic<double,30> &Individual::get_wscal_30() const {
	if (use_indiv_wscal()) {
		return wscal_30;
	} else {
		return patchpft().wscal_30;
	}
}

void Individual::validate_canopy() const {
	double lai_sum = 0;
	for (size_t j = 0; j < canopy.size(); j++) {
		lai_sum += canopy[j].lai;
		if (j > 0) {
			// fipar should decrease as we go down the canopy.
			if (canopy[j].fipar > canopy[j - 1].fipar) {
				log_fatal("Invalid canopy for individual %d (%s): fipar should decrease as we go down the canopy, but fipar in layer %d is %g and in layer %d is %g.",
						  id, (const char *)pft.name, j, canopy[j].fipar,
						  j - 1, canopy[j - 1].fipar);
			}
		}
	}

	// Ensure that entire LAI is allocated across layers.
	if (!equal(lai_sum, lai_today())) {
		log_fatal("Invalid canopy for individual %d (%s): total LAI across layers (%g) does not match today's LAI (%g).",
				  id, (const char *)pft.name, lai_sum, lai_today());
	}
}

double Individual::get_fwuptake(size_t layer) const {
	if (pft.dave_allocation()) {
		return fwuptake[layer];
	}
	return vegetation.patch.pft[pft.id].fwuptake[layer];
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

	distinterval_gc = 1.0e10;
	for (int m=0;m<12;m++) {
		monthly_burned_area[m] = 0.0;
		monthly_fire_risk[m] = 0.0;

	}
	burned_area = 0.0;
	simfire_region = 0;
	simfire_biome = 0;
	pop_density = 0;
	simfire_annual_burned_area = 0;
	nesterov_max = 0;
	annual_burned_area = 0;
	k_tun_litter = 0;
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

void Gridcell::set_soillayer_properties() {

	double soil_lyr_thick[NSOILLAYERS_MAX];
	int nsoillayers = 0;
	
	soillayer_thickness(soil_lyr_thick, nsoillayers);

	// Set fixed rootdist per PFT depending on soil depth and layer thickness
	for (unsigned int p = 0; p<pftlist.nobj; p++) {
		pft[p].pft.init_rootdist(soil_lyr_thick, soil_depth, nsoillayers);
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
					patch.pft[i].rootdist[lyr] = pft[i].pft.input_rootdist[lyr];
				}
			}

			// Get reference to this soil 
			Soil& soil = patch.soil;

			// Set soil depth
			soil.depth = soil_depth;

			// Set number of soil layer in soil
			soil.nsoillayers = nsoillayers;

			// Create soil layers and initalise them
			soil.create_soillayers(soil_lyr_thick);

			// Set padding layers depth
			soil.pad.set_depths(soil.soillayers[nsoillayers - 1].thickness);
			
			stand.nextobj();
		}// End of loop through patches

		++gc_itr;
	}	// End of loop through stands
}

void Gridcell::soillayer_thickness(double lyr_thickness[NSOILLAYERS_MAX], int& nsoillayers) {

	if (soil_depth < SOILDEPTH_EVAP) {
		fail("Soil depth (%g) needs to be thicker than evaporation depth (%g)\n", soil_depth, SOILDEPTH_EVAP);
	}

	if (hydrology == RICHARDS) {

		// Standard soil layer thickness (mm)
		double std_thick[NSOILLAYERS_MAX] = { 100.0, 100.0, 100.0, 200.0, 200.0, 300.0, 500.0, 500.0, 1000.0, 1000.0, 1500.0, 3000.0, 4000.0, 6000.0, 8000.0 };
		
		double std_tot_thick = 0.0;
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) std_tot_thick += std_thick[lyr];

		// Calculate number of soil layers and their individual depths
		if (soil_depth > std_tot_thick) {

			// Deeper soil than total standard soil layer depths

			// Number of layers that will be adjusted
			int nadj_lyr = 6;
		
			double tot_thick = 0.0;
			// Bottom three soil layers will adjust for
			// deeper soil than standard depth
			for (int lyr = 0; lyr < NSOILLAYERS_MAX - nadj_lyr; lyr++) {
				lyr_thickness[lyr] = std_thick[lyr];
				tot_thick += std_thick[lyr];
			}

			// Upscaling fraction to remaining soil layers
			double scalar = (soil_depth - tot_thick) / (std_tot_thick - tot_thick);
		
			// Adjust soil layer thickness to be off 100 mm increments
			for (int lyr = NSOILLAYERS_MAX - nadj_lyr; lyr < NSOILLAYERS_MAX; lyr++) {
				lyr_thickness[lyr] = floor(std_thick[lyr] * scalar / 100.0) * 100.0;
				tot_thick += lyr_thickness[lyr];
			}

			lyr_thickness[NSOILLAYERS_MAX - 1] += soil_depth - tot_thick;

			nsoillayers = NSOILLAYERS_MAX;
		}
		else {

			// Shallower soil than total standard soil layer depths
			// Adjust soil layer thickness to be off 100 mm increments

			int lyr = 0;
			double tot_thick = 0.0;
			double remain_depth = soil_depth;

			while (!negligible(remain_depth)) {
				double ref_depth = lyr == 0 ? std_thick[lyr] : std_thick[lyr - 1];
				if (remain_depth >= ref_depth * 2.0 && lyr < NSOILLAYERS_MAX - 1) {
					double lyr_depth = remain_depth < std_thick[lyr] * 2.0 ? remain_depth / 2.0 : std_thick[lyr];
					lyr_thickness[lyr] = floor(lyr_depth / 100.0) * 100.0;
					tot_thick += lyr_thickness[lyr];
					remain_depth = soil_depth - tot_thick;
				}
				else {
					lyr_thickness[lyr] = remain_depth;
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
			lyr_thickness[lyr] = soil_depth / NSOILLAYERS_MAX;
		}
		nsoillayers = NSOILLAYERS_MAX;
	}
	else if (hydrology == GERTEN) {
		// Upper layer 500mm thick and lower layer 1000mm
		lyr_thickness[0] = SOILDEPTH_UPPER;
		lyr_thickness[1] = soil_depth - lyr_thickness[0];

		nsoillayers = 2;
	}
}

Stand& Gridcell::create_stand_lu(StandType& st, double fraction, int no_patch, bool suppress_disturbance) {

	Stand& stand = create_stand(st.landcover, no_patch);
	stand.init_stand_lu(st, fraction, suppress_disturbance);

	return stand;
}

double Gridcell::ccont() const {

	double ccont = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ccont += stand.ccont() * stand.get_gridcell_fraction();
	}

	return ccont;
}

double Gridcell::ncont() const {

	double ncont = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		ncont += stand.ncont() * stand.get_gridcell_fraction();
	}

	return ncont;
}

double Gridcell::water_content() const {

	double water_content = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		water_content += stand.water_content() * stand.get_gridcell_fraction();
	}

	return water_content;
}

double Gridcell::cflux(Period period) const {

	double cflux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		cflux += stand.cflux(period) * stand.get_gridcell_fraction();
	}

	cflux += landcover.acflux_landuse_change;
	cflux += landcover.acflux_harvest_slow;
	cflux += landcover.acflux_wood_harvest;
	cflux += landcover.acflux_clearing;

	return cflux;
}

double Gridcell::nflux(Period period) const {

	double nflux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		nflux += stand.nflux(period) * stand.get_gridcell_fraction();
	}

	// TODO DH: these won't be accounted for via the period argument. We need
	// to create per-stand and per-gridcell flux objects to track these kinds
	// of fluxes. Everything should be a flux!
	nflux += landcover.anflux_landuse_change;
	nflux += landcover.anflux_harvest_slow;
	nflux += landcover.anflux_wood_harvest;
	nflux += landcover.anflux_clearing;

	return nflux;
}

double Gridcell::water_flux(Period period) const {

	double water_flux = 0.0;

	for (unsigned int s = 0; s < nbr_stands(); s++) {
		Stand& stand = (*this)[s];
		water_flux += stand.water_flux(period) * stand.get_gridcell_fraction();
	}

	return water_flux;
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
		& water_flux
		& ccont_yesterday
		& ncont_yesterday
		& water_cont_yesterday;
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
	cflux_zero = patch.cflux(Period::YEAR);

	ncont_zero = indiv.ncont();
	ncont_zero_scaled = indiv.ncont(scale, true);
	// Add soil N
	ncont_zero += patch.ncont(0.0);
	ncont_zero_scaled += patch.ncont(0.0);
	nflux_zero = patch.nflux(Period::YEAR);
}

bool MassBalance::check_indiv_C(Individual& indiv, bool check_harvest) {

	bool balance = true;
	Patch& patch = indiv.vegetation.patch;
	Stand& stand = patch.stand;
	if (!stand.is_true_crop_stand())
		return balance;

	double ccont = indiv.ccont();
	ccont += patch.ccont(0.0);
	double cflux = patch.cflux(Period::YEAR);

	if (check_harvest && patch.isharvestday)
		ccont_zero = ccont_zero_scaled;

	if (date.year >= nyear_spinup && !negligible(ccont - ccont_zero + cflux - cflux_zero, -10)) {
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
	if (!stand.is_true_crop_stand())
		return balance;

	double ncont = indiv.ncont();
	ncont += patch.ncont(0.0);
	double nflux = patch.nflux(Period::YEAR);

	if (check_harvest && patch.isharvestday)
		ncont_zero = ncont_zero_scaled;

	if (date.year >= nyear_spinup && !negligible(ncont - ncont_zero + nflux - nflux_zero, -14)) {
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
	cflux_zero = patch.cflux(Period::YEAR);

	if (stand.get_gridcell_fraction())
		cflux_zero += gridcell.landcover.acflux_harvest_slow / stand.get_gridcell_fraction();

	ncont_zero = patch.ncont();
	ncont_zero_scaled = patch.ncont(scale, true);
	nflux_zero = patch.nflux(Period::YEAR);

	if (stand.get_gridcell_fraction())
		nflux_zero += gridcell.landcover.anflux_harvest_slow / stand.get_gridcell_fraction();

	water_cont_zero = patch.water_content();
	water_flux_zero = patch.water_flux(Period::YEAR);
}

bool MassBalance::check_patch_C(Patch& patch, bool check_harvest) {

	bool balance = true;
	Stand& stand = patch.stand;
	Gridcell& gridcell = stand.get_gridcell();
	double ccont = patch.ccont();
	double cflux = patch.cflux(Period::YEAR);

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
	double nflux = patch.nflux(Period::YEAR);

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
	double water_flux = patch.water_flux(Period::YEAR);

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
	double nflux_year = gridcell.nflux(Period::YEAR);

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

void MassBalance::rebase_day(Gridcell &gridcell) {
	ccont_yesterday = gridcell.ccont();
	ncont_yesterday = gridcell.ncont();
	water_cont_yesterday = gridcell.water_content();
}

void MassBalance::check_day(Gridcell &gridcell) {

	if (ifdavegrowth) {
		// Sub-annual mass balance is not preserved in the annual allocation
		// scheme (ie trunk allocation). This is because carbon is assimilated
		// (and that amount reported as NPP flux) on a daily basis, but the
		// carbon is not allocated until the end of the year. My workaround is
		// to count indiv.anpp inside indiv.ccont(), but this is a hack and
		// should ultimately be removed.

		check_day_c(gridcell);

		if (ifcentury) {
			check_day_n(gridcell);
		}
	}

	check_day_water(gridcell);
}

bool is_first_day() {
	if (date.day == 0 && date.year == 0) {
		return true;
	}
	if (date.day == 0 && date.year == nyear_spinup) {
		return true;
	}
	return false;
}

bool check_mass_balance_today(const Gridcell &gridcell) {
	int solvesomcent_endyr   = get_solvesomcent_endyr();
	bool solve_somcent_today = date.islastday && date.islastmonth && date.year == solvesomcent_endyr;
	return ifdavephotosyn && !solve_somcent_today && !is_first_day();
}

void MassBalance::check_day_c(Gridcell &gridcell) {
	double ccont = gridcell.ccont();
	double cflux = gridcell.cflux(Period::DAY);

	if (check_mass_balance_today(gridcell)) {
		// double cflux_today = cflux_yesterday - cflux;
		double expected_ccont = ccont_yesterday - cflux;

		const double epsilon = 1e-9;

		double err = ccont - expected_ccont;
		if (fabs(err) > epsilon) {
			dprintf("\n(%.2f, %.2f): C balance error year %d day %d: %.10f\n"
				, gridcell.get_lon(), gridcell.get_lat(),
				date.year, date.day, err);
			dprintf("C pool change: %.10f\n", ccont - ccont_yesterday);
			dprintf("C flux: %.10f\n", cflux);
		}
	}

	ccont_yesterday = ccont;
}

void MassBalance::check_day_n(Gridcell &gridcell) {
	double ncont = gridcell.ncont();
	double nflux = gridcell.nflux(Period::DAY); // NO_SOIL, N2O_SOIL, gridcell.aNH4dep, gridcell.aNO3dep

	if (check_mass_balance_today(gridcell)) {
		// double nflux_today = nflux_yesterday - nflux;
		double expected_ncont = ncont_yesterday - nflux;

		const double epsilon = 1e-9;

		double err = ncont - expected_ncont;
		if (fabs(err) > epsilon) {
			dprintf("\n(%.2f, %.2f): N balance year %d day %d: %.10f\n"
				, gridcell.get_lon(), gridcell.get_lat(),
				date.year, date.day, err);
			dprintf("N pool change: %.10f\n", ncont - ncont_yesterday);
			dprintf("N flux: %.10f\n", nflux);
		}
	}

	ncont_yesterday = ncont;
}

void MassBalance::check_day_water(Gridcell &gridcell) {
	double wcont = gridcell.water_content();
	double wflux = gridcell.water_flux(Period::DAY);

	if (check_mass_balance_today(gridcell)) {
		// double wflux_today = water_flux_yesterday - wflux;
		double expected_wcont = water_cont_yesterday - wflux;

		// Ensure wcont matches yesterday's value + today's water flux.
		double err = wcont - expected_wcont;
		if (fabs(err) > EPS_WATER) {
			dprintf("\n(%.2f, %.2f): Water balance year %d day %d: %.10f\n"
				, gridcell.get_lon(), gridcell.get_lat(),
				date.year, date.day, err);
			dprintf("Water content change: %.10f\n", wcont - water_cont_yesterday);
			dprintf("Water flux: %.10f\n", wflux);
		}
	}

	water_cont_yesterday = wcont;
}

void MassBalance::check_year_C(Gridcell& gridcell) {

	double ccont_year = gridcell.ccont();
	double cflux_year = gridcell.cflux(Period::YEAR);

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
			dprintf("C flux: %.5f\n",  cflux_year);
			dprintf("Imbalance: %.9f\n", ccont_year - ccont + cflux_year);
		}
	}

	ccont = ccont_year;
}

void MassBalance::check_year_water(Gridcell& gridcell) {

	double water_cont_year = gridcell.water_content();
	double water_flux_year = gridcell.water_flux(Period::YEAR);

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
	if ((ncont - ncont_zero + nflux) > epsilon_biomass) {
		dprintf("\nWARNING: (%.2f, %.2f): Period N balance year %d: %.10f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, ncont - ncont_zero + nflux);
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
	cflux_zero = gridcell.cflux(Period::YEAR);
	ncont_zero = gridcell.ncont();
	nflux_zero = gridcell.nflux(Period::YEAR);
	water_cont_zero = gridcell.water_content();
	water_flux_zero = gridcell.water_flux(Period::YEAR);
}

void MassBalance::check(Gridcell& gridcell) {

	double ccont = gridcell.ccont();
	double cflux = gridcell.cflux(Period::YEAR);

	if (!negligible(ccont - ccont_zero + cflux, -5)) {
		dprintf("\n(%.2f, %.2f): C balance year %d: %.5f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, ccont - ccont_zero + cflux);
		dprintf("C pool change: %.5f\n", ccont - ccont_zero);
		dprintf("C flux: %.5f\n\n",  cflux);
	}

	double ncont = gridcell.ncont();
	double nflux = gridcell.nflux(Period::YEAR);

	if (!negligible(ncont - ncont_zero + nflux, -5)) {
		dprintf("\n(%.2f, %.2f): N balance year %d: %.5f\n", gridcell.get_lon(), gridcell.get_lat(), date.year, ncont - ncont_zero + nflux);
		dprintf("N pool change: %.5f\n", ncont - ncont_zero);
	}

	double water_content = gridcell.water_content();
	double water_flux = gridcell.water_flux(Period::YEAR);

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

inline double calc_temp_growth(double a, const Pft *pft) {
	return temp_growth(a, pft->Ha_growth, pft->T0_growth, 1.0, pft->temp_dormancy);
}

double calc_a_growth(const Pft *pft) {
	double xmin = 0;
	double xmax = 50;
	const size_t NSTEP_MAX = 100;
	const double xres = 0.001;
	double ymin = calc_temp_growth(xmin, pft);
	double ymax = calc_temp_growth(xmax, pft);
	double y;

	size_t i;
	for (i = 0; i < NSTEP_MAX; i++) {
		double x = (xmin + xmax) / 2.0;
		y = calc_temp_growth(x, pft);
		if (fabs(x - xmin) < xres) {
			break;
		}
		double ynext = calc_temp_growth(x + 0.01, pft);
		if (ynext > y) {
			xmin = x;
		} else {
			xmax = x;
		}
	}
	if (i >= NSTEP_MAX - 1) {
		// Failed to find a solution - perhaps some of the parameters were not
		// supplied in the instruction file.
		fail("Failed to optimise temp_growth curve - check pft parameters Ha_growth, T0_growth, temp_dormancy!!!");
	}
	return 1.0 / y;
}

void Canopy::clear() {
	// Reset all daily values.
	for (size_t i = 0; i < size(); i++) {
		layers[i].lai = 0.0;
		layers[i].fipar = 0.0;
	}
}

void Canopy::set(size_t layer, double lai, double fipar) {
	verify_size(layer);
	layers[layer].lai = lai;
	layers[layer].fipar = fipar;
}

CanopyLayer &Canopy::operator[] (size_t i) {
	verify_size(i);
	return layers[i];
}

const CanopyLayer &Canopy::operator[] (size_t i) const {
	// FIXME: only necessary until daveoutput refactor to improve handling of
	// cohort-level outputs.
	if (i >= size()) {
		// Accessing a layer that doesn't exist - return an empty layer.
		static CanopyLayer empty_layer;
		return empty_layer;
	}
	return layers[i];
}

void Canopy::verify_size(size_t size) {
	// Ensure that layer `size` exists.
	while (this->size() < size + 1) {
		layers.push_back(CanopyLayer{});
	}
}

size_t Canopy::size() const {
	return layers.size();
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995.
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
