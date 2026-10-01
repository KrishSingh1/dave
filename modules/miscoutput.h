///////////////////////////////////////////////////////////////////////////////////////
/// \file cropoutput.h
/// \brief Output module for the land use and management information
///
/// \author Joe Siltberg
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_MISC_OUTPUT_H
#define LPJ_GUESS_MISC_OUTPUT_H

#include "outputmodule.h"
#include "outputchannel.h"
#include "gutil.h"

namespace GuessOutput {

// Definitions for separate output files per NATURAL and FOREST stand (when
// instruction file parameter printseparatestands == true, printseparatestands
// is set to false in LandcoverInput::init() when input land cover
// fraction data file has data for > 50 gridcells)

/// Output module for the most commonly needed output files
class MiscOutput : public OutputModule {
public:

	MiscOutput();

	~MiscOutput();

	// implemented functions inherited from OutputModule
	// (see documentation in OutputModule)

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);

private:

	/// Upper limit for files in multiple stand printout
	static const int MAXNUMBER_STANDS = 1000;

	/// Upper limit for files in multiple pft printout
	static const int MAXNUMBER_PFTS = 100; // from EC-Earth

	/// Printout of first stand from first historic year
	static const bool PRINTFIRSTSTANDAFTERSPINUP = true;

	/// Defines all output tables
	void openlocalfiles(Gridcell& gridcell, int coordinates_precision);

	void define_output_tables();

	void closelocalfiles(Gridcell& gridcell);

	// Output file names ...
	xtring file_yield_active, file_yield1_active, file_yield2_active, file_sdate1, file_sdate2,
		   file_hdate1, file_hdate2, file_lgp, file_phu, file_fphu, file_fhi,
		   file_irrigation, file_seasonality, file_cflux_cropland,
		   file_cflux_pasture, file_cflux_natural, file_cflux_forest,
		   file_cpool_cropland, file_cpool_pasture, file_cpool_natural,
		   file_cpool_forest, file_nflux_cropland, file_nflux_pasture,
		   file_nflux_natural, file_nflux_forest, file_npool_cropland,
		   file_npool_pasture, file_npool_natural, file_npool_forest,
		   file_anpp_active_cropland, file_anpp_active_pasture, file_anpp_active_natural,
		   file_anpp_active_forest, file_cmass_active_cropland, file_cmass_active_pasture,
		   file_cmass_active_natural, file_cmass_active_forest, file_dens_natural,
		   file_anpp_cropland, file_anpp_pasture,
		   file_anpp_natural, file_anpp_forest, file_anpp_peatland,
		   file_cmass_cropland, file_cmass_pasture,
		   file_cmass_natural, file_cmass_forest, file_cmass_peatland,
		   file_lai_cropland, file_lai_pasture,
		   file_lai_natural, file_lai_forest, file_lai_peatland,
		   file_dens_forest, file_soil_nflux_cropland, file_soil_nflux_pasture,
		   file_soil_nflux_natural, file_soil_nflux_forest,
		   file_cmass_active_peatland, file_cflux_peatland,
		   file_cpool_peatland, file_nflux_peatland, file_npool_peatland,
		   file_anpp_active_peatland,
		   file_agestruct_natural, file_agestruct_forest, file_diamstruct_natural,
		   file_diamstruct_forest, file_diamstruct_cmass_natural, file_diamstruct_cmass_forest,
		   file_diamstruct_cmass_leaf_natural, file_diamstruct_cmass_leaf_forest,
		   file_diamstruct_cmass_wood_natural, file_diamstruct_cmass_wood_forest,
		   file_diamstruct_cmass_froot_natural, file_diamstruct_cmass_froot_forest,
		   file_diamstruct_cmass_wood_inc_natural, file_diamstruct_cmass_wood_inc_forest,
		   file_diamstruct_ba_natural, file_diamstruct_ba_forest,
		   file_agestruct_standage_natural, file_agestruct_standage_forest,
		   file_agestruct_standage_cropland, file_agestruct_standage_pasture,
		   file_aaet_active_natural, file_aaet_active_forest, file_agpp_active_natural, file_agpp_active_forest,
		   file_speciesdiam_natural, file_speciesdiam_forest,
		   file_speciesheights_natural, file_speciesheights_forest,
		   file_lai_active_natural, file_lai_active_forest, file_fpc_active_natural, file_fpc_active_forest,
		   file_forest_cmass_harv_killed_active, file_forest_vegc, file_forest_cflux_veg, file_forest_harvest,
		   file_harvest_flux_luc, file_cflux_forestry, file_cflux_regrowth, file_cflux_primary,
		   file_cpool_forestry, file_cpool_regrowth, file_cpool_primary;

	xtring file_frac_natural, file_frac_forest;

	xtring file_cmass_active_plus10cm, file_forest_cmass_harv_killed_active_plus10cm, file_dens_plus10cm, file_speciesheights_plus10cm, file_speciesdiam_plus10cm;

	// Stand type output file_cmass_harv_killed_sts, file_cmass_thin_killed_sts, file_cmass_dist_default_sts, file_cmass_tree_mort_sts
	xtring file_anpp_sts, file_agpp_sts, file_lai_sts, file_cmass_sts, file_cmass_harv_killed_sts, file_cmass_thin_killed_sts;
	xtring file_anpp_tree_sts, file_agpp_tree_sts, file_aautresp_tree_sts, file_aautresp_root_tree_sts, file_aautresp_sap_tree_sts, file_aautresp_leaf_tree_sts, file_aautresp_growth_tree_sts;

	xtring file_lai_tree_sts, file_cmass_tree_sts, file_cmass_tree_mort_sts;
	xtring file_frac_sts;

	xtring file_cmass_wood_sts, file_cmass_wood_harv_killed_sts, file_cmass_wood_harv_sts, file_cmass_wood_harv_toprod_sts,
		   file_cmass_wood_thin_sts, file_cmass_wood_clearcut_sts, file_cmass_wood_mort_sts, file_cmass_wood_inc_sts;

	xtring file_cmass_clearcut_killed_sts, file_cmass_fire_killed_sts, file_cmass_wind_killed_sts, file_cmass_barkbeetle_killed_sts;
	xtring file_cmass_dist_default_sts, file_cmass_dist_prescribed_sts, file_cmass_dist_clearcut_killed_sts, file_cmass_dist_fire_killed_sts, file_cmass_dist_wind_killed_sts, file_cmass_dist_barkbeetle_killed_sts;

	xtring file_csoil_sts, file_clitter_sts, file_clitter_course_sts, file_clitter_fine_sts, file_csink_sts, file_luc_cflux_sts, file_nsoil_sts, file_nlitter_sts, file_nlitter_course_sts, file_nlitter_fine_sts;
	xtring file_rh_sts, file_rh_clitter_course_sts, file_rh_clitter_fine_sts;

	xtring file_dens_sts, file_height_sts, file_diam_sts, file_diam_g_sts;

	xtring file_cmass_plus10cm_sts, file_cmass_harv_killed_plus10cm_sts;
	xtring file_cmass_dist_prescribed_sts_plus10cm_sts, file_cmass_dist_clearcut_killed_plus10cm_sts, file_cmass_dist_fire_killed_plus10cm_sts, file_cmass_dist_wind_killed_plus10cm_sts, file_cmass_dist_barkbeetle_killed_plus10cm_sts;

	xtring file_cmass_wood_plus10cm_sts, file_cmass_wood_harv_killed_plus10cm_sts, file_cmass_wood_mort_plus10cm_sts, file_cmass_wood_inc_plus10cm_sts;

	xtring file_dens_plus10cm_sts, file_height_plus10cm_sts, file_diam_plus10cm_sts, file_diam_g_plus10cm_sts;

	xtring file_height_sts_patch_mean, file_diam_sts_patch_mean, file_diam_g_sts_patch_mean;
	xtring file_height_plus10cm_sts_patch_mean, file_diam_plus10cm_sts_patch_mean, file_diam_g_plus10cm_sts_patch_mean;

	xtring file_cutinterval_sts, file_cutinterval_thisyear_sts, file_nstand_sts, file_npatch_sts;

	xtring file_transpiration_sts, file_runoff_sts, file_intercep_sts, file_evap_sts;

	xtring file_storm, file_storm_si_sts, file_storm_frac_sts, file_storm_frac_pot_sts, file_storm_dw_sts, file_storm_dw_pot_sts, file_insect_dw_sts;

	xtring file_cmort_storm_sts, file_cmort_insect_sts, file_cmort_fire_sts, file_cmort_clearcut_sts;

	xtring file_cflux_harvest;
	xtring file_harvest_volume;

	// soil n2o flux per landcover, info needed for estimating global n2o emissions
	// soil n2o from BARREN which is LUH2 icwtr shows up in the n2o in soil_nflux
	// TODO: allow BARREN to get excluded from soil_nflux
	xtring file_soil_n2o_flux;

	// daily
	xtring file_daily_lai, file_daily_npp, file_daily_nmass, file_daily_cmass,
		   file_daily_cton, file_daily_ndemand, file_daily_cmass_leaf,
		   file_daily_nmass_leaf, file_daily_cmass_root, file_daily_nmass_root,
		   file_daily_cmass_stem, file_daily_nmass_stem,
		   file_daily_cmass_storage, file_daily_nmass_storage,
		   file_daily_n_input_soil, file_daily_avail_nmass_soil,
		   file_daily_upper_wcont, file_daily_lower_wcont,
		   file_daily_irrigation, file_daily_climate,
		   file_daily_cmass_dead_leaf, file_daily_nmass_dead_leaf, 
		   file_daily_fphu, file_daily_nminleach,
		   file_daily_norgleach, file_daily_nuptake, file_daily_ds,
		   file_daily_stem, file_daily_leaf, file_daily_root,
		   file_daily_storage,file_daily_blaze;
	xtring file_monthly_mlambda; 

	xtring file_aaet_gridcell_natural, file_aaet_gridcell_forest, file_lai_gridcell_natural, file_lai_gridcell_forest, 
		   file_anpp_gridcell_natural, file_anpp_gridcell_forest, file_cmass_gridcell_natural, file_cmass_gridcell_forest;

	// Output tables
	Table out_yield_active, out_yield1_active, out_yield2_active, out_sdate1, out_sdate2,
		  out_hdate1, out_hdate2, out_lgp, out_phu, out_fhi, out_fphu,
		  out_irrigation, out_seasonality, out_cflux_cropland,
		  out_cflux_pasture, out_cflux_natural, out_cflux_forest,
		  out_cpool_cropland, out_cpool_pasture, out_cpool_natural,
		  out_cpool_forest, out_nflux_cropland, out_nflux_pasture,
		  out_nflux_natural, out_nflux_forest, out_npool_cropland,
		  out_npool_pasture, out_npool_natural, out_npool_forest,
		  out_anpp_active_cropland, out_anpp_active_pasture, out_anpp_active_natural,
		  out_anpp_active_forest, out_cmass_active_cropland, out_cmass_active_pasture,
		  out_cmass_active_natural, out_cmass_active_forest, out_dens_natural, out_dens_forest,
		  out_anpp_cropland, out_anpp_pasture,
		  out_anpp_natural, out_anpp_forest, out_anpp_peatland,
		  out_cmass_cropland, out_cmass_pasture,
		  out_cmass_natural, out_cmass_forest, out_cmass_peatland,
		  out_lai_cropland, out_lai_pasture,
		  out_lai_natural, out_lai_forest, out_lai_peatland,
		  out_soil_nflux_cropland, out_soil_nflux_pasture,
		  out_soil_nflux_natural, out_soil_nflux_forest,
		  out_cflux_peatland, out_cpool_peatland,
		  out_nflux_peatland, out_npool_peatland, out_cmass_active_peatland,
		  out_anpp_active_peatland,
		  out_agestruct_natural, out_agestruct_forest, out_diamstruct_natural,
		  out_diamstruct_forest, out_diamstruct_cmass_natural, out_diamstruct_cmass_forest,
		  out_diamstruct_cmass_leaf_natural, out_diamstruct_cmass_leaf_forest,
		  out_diamstruct_cmass_wood_natural, out_diamstruct_cmass_wood_forest,
		  out_diamstruct_cmass_froot_natural, out_diamstruct_cmass_froot_forest,
		  out_diamstruct_cmass_wood_inc_natural, out_diamstruct_cmass_wood_inc_forest,
		  out_diamstruct_ba_natural, out_diamstruct_ba_forest,
		  out_agestruct_standage_natural, out_agestruct_standage_forest,
		  out_agestruct_standage_cropland, out_agestruct_standage_pasture,
		  out_aaet_active_natural, out_aaet_active_forest, out_speciesdiam_natural, out_speciesdiam_forest,
		  out_speciesheights_natural, out_speciesheights_forest,
		  out_lai_active_natural, out_lai_active_forest, out_fpc_active_natural, out_fpc_active_forest,
		  out_forest_harvest, out_forest_vegc, out_forest_cflux_veg, out_forest_cmass_harv_killed_active,
		  out_harvest_flux_luc, out_cflux_forestry, out_cflux_regrowth, out_cflux_primary,
		  out_cpool_forestry, out_cpool_regrowth, out_cpool_primary,
		  out_lai_gridcell_natural, out_lai_gridcell_forest, out_anpp_gridcell_natural, out_anpp_gridcell_forest, 
		  out_cmass_gridcell_natural, out_cmass_gridcell_forest, out_aaet_gridcell_natural, out_aaet_gridcell_forest,
		  out_mlambda; // EC-Earth only;

	Table out_cflux_harvest;
	Table out_harvest_volume;

	Table out_frac_natural, out_frac_forest;

	Table out_cmass_active_plus10cm, out_forest_cmass_harv_killed_active_plus10cm, out_dens_plus10cm, out_speciesheights_plus10cm, out_speciesdiam_plus10cm;

	// Output files with stand type columns
	Table out_anpp_sts, out_agpp_sts, out_lai_sts, out_cmass_sts, out_cmass_harv_killed_sts, out_cmass_thin_killed_sts;
	Table out_anpp_tree_sts, out_agpp_tree_sts, out_aautresp_tree_sts, out_aautresp_root_tree_sts,out_aautresp_sap_tree_sts, out_aautresp_leaf_tree_sts, out_aautresp_growth_tree_sts;
	Table out_cmass_dist_default_sts;

	Table out_cmass_clearcut_killed_sts, out_cmass_fire_killed_sts, out_cmass_wind_killed_sts, out_cmass_barkbeetle_killed_sts;
	Table out_cmass_dist_prescribed_killed_sts, out_cmass_dist_clearcut_killed_sts, out_cmass_dist_fire_killed_sts, out_cmass_dist_wind_killed_sts, out_cmass_dist_barkbeetle_killed_sts;

	Table out_lai_tree_sts, out_cmass_tree_sts, out_cmass_tree_mort_sts;
	Table out_frac_sts;

	Table out_cmass_wood_sts, out_cmass_wood_harv_killed_sts, out_cmass_wood_harv_sts, out_cmass_wood_harv_toprod_sts, out_cmass_wood_thin_sts, 
		out_cmass_wood_clearcut_sts, out_cmass_wood_mort_sts, out_cmass_wood_inc_sts;

	Table out_csoil_sts, out_clitter_sts, out_clitter_course_sts, out_clitter_fine_sts, out_csink_sts, out_luc_cflux_sts, out_nsoil_sts, out_nlitter_sts, out_nlitter_course_sts, out_nlitter_fine_sts;
	Table out_rh_sts, out_rh_clitter_course_sts, out_rh_clitter_fine_sts;

	Table out_dens_sts, out_height_sts, out_diam_sts, out_diam_g_sts;

	Table out_cmass_plus10cm_sts, out_cmass_harv_killed_plus10cm_sts;
	Table out_cmass_dist_prescribed_killed_plus10cm_sts, out_cmass_dist_clearcut_killed_plus10cm_sts, out_cmass_dist_fire_killed_plus10cm_sts, out_cmass_dist_wind_killed_plus10cm_sts, out_cmass_dist_barkbeetle_killed_plus10cm_sts;
	Table out_cmass_wood_plus10cm_sts, out_cmass_wood_harv_killed_plus10cm_sts, out_cmass_wood_mort_plus10cm_sts, out_cmass_wood_inc_plus10cm_sts;

	Table out_dens_plus10cm_sts, out_height_plus10cm_sts, out_diam_plus10cm_sts, out_diam_g_plus10cm_sts;

	Table out_height_sts_patch_mean, out_diam_sts_patch_mean, out_diam_g_sts_patch_mean;
	Table out_height_plus10cm_sts_patch_mean, out_diam_plus10cm_sts_patch_mean, out_diam_g_plus10cm_sts_patch_mean;

	Table out_cutinterval_sts, out_cutinterval_thisyear_sts, out_nstand_sts, out_npatch_sts;

	Table out_transpiration_sts, out_runoff_sts, out_intercep_sts, out_evap_sts;

	Table out_storm, out_storm_si_sts, out_storm_frac_pot_sts, out_storm_frac_sts, out_storm_dw_sts, out_storm_dw_pot_sts, out_insect_dw_sts;

	Table out_cmort_storm_sts, out_cmort_insect_sts, out_cmort_fire_sts, out_cmort_clearcut_sts;

	Table out_soil_n2o_flux;

	// Separate output files for stand types with pft columns
	Table* out_cmass_pft_st;
	Table* out_cmass_harv_killed_pft_st;
	Table* out_transpiration_pft_st;

	// Separate output files for stand types with diameter-class columns
	Table* out_diamstruct_cmass_st;
	Table* out_diamstruct_cmass_leaf_st;
	Table* out_diamstruct_cmass_wood_st;
	Table* out_diamstruct_cmass_froot_st;
	Table* out_diamstruct_cmass_wood_inc_st;

	Table* out_standagestruct_st;

	// Separate output files for stands with pft columns
	Table* out_anpp_stand[MAXNUMBER_STANDS];
	Table* out_lai_stand[MAXNUMBER_STANDS];
	Table* out_cmass_stand[MAXNUMBER_STANDS];
	Table* out_diam_stand[MAXNUMBER_STANDS];
	Table* out_height_stand[MAXNUMBER_STANDS];
	Table* out_dens_stand[MAXNUMBER_STANDS];
	Table* out_cmass_mort_stand[MAXNUMBER_STANDS];
	Table* out_cmass_wood_stand[MAXNUMBER_STANDS];
	Table* out_cmass_wood_harv_stand[MAXNUMBER_STANDS];
	Table* out_storm_stand[MAXNUMBER_STANDS];

	// Separate output files for stands with age-class columns
	Table* out_agestruct_stand[MAXNUMBER_STANDS];

	// Separate output files for stands with diameter-class columns
	Table* out_diamstruct_stand[MAXNUMBER_STANDS];
	Table* out_diamstruct_cmass_stand[MAXNUMBER_STANDS];

	//daily
	Table out_daily_lai, out_daily_npp, out_daily_cton, out_daily_nmass,
		  out_daily_cmass, out_daily_ndemand, out_daily_cmass_leaf,
		  out_daily_nmass_leaf, out_daily_cmass_root, out_daily_nmass_root,
		  out_daily_cmass_stem, out_daily_nmass_stem, out_daily_cmass_storage,
		  out_daily_nmass_storage, out_daily_n_input_soil,
		  out_daily_cmass_dead_leaf, out_daily_nmass_dead_leaf, out_daily_fphu,
		  out_daily_avail_nmass_soil, out_daily_upper_wcont,
		  out_daily_lower_wcont, out_daily_irrigation, out_daily_climate,
		  out_daily_nminleach, out_daily_norgleach, out_daily_nuptake, out_daily_ds, 
		  out_daily_stem, out_daily_leaf, out_daily_root, out_daily_storage, out_daily_blaze;

	bool print_anpp_stand;
	bool print_lai_stand;
	bool print_cmass_stand;
	bool print_cmass_wood_stand;
	bool print_cmass_wood_harv_stand;
	bool print_cmass_mort_stand;
	bool print_height_stand;
	bool print_diam_stand;
	bool print_dens_stand;
	bool print_agestruct_stand;
	bool print_diamstruct_stand;
	bool print_diamstruct_cmass_stand;

	bool printstandtypes;
	bool print_cmass_pft_st;
	bool print_cmass_harv_killed_pft_st;
	bool print_transpiration_pft_st;
	bool print_diamstruct_cmass_st;
	bool print_diamstruct_cmass_leaf_st;
	bool print_diamstruct_cmass_wood_st;
	bool print_diamstruct_cmass_froot_st;
	bool print_diamstruct_cmass_wood_inc_st;

	bool print_standagestruct_st;
};

}

#endif // LPJ_GUESS_MISC_OUTPUT_H
