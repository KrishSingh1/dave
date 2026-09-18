///////////////////////////////////////////////////////////////////////////////////////
/// \file outputmodule.cpp
/// \brief Implementation of the common output module
///
/// \author Joe Siltberg
/// $Date: 2015-04-09 18:40:34 +0200 (Thu, 09 Apr 2015) $
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include <filesystem>

#include "config.h"
#include "miscoutput.h"
#include "parameters.h"
#include "guess.h"
#include "management.h"
#include <sys/stat.h>

namespace GuessOutput {

// Nitrogen output is in kgN/ha instead of kgC/m2 as for carbon

REGISTER_OUTPUT_MODULE("misc", MiscOutput)

MiscOutput::MiscOutput() {
	// Annual output variables
	declare_parameter("file_cmass_active_cropland", &file_cmass_active_cropland, 300, "Annual cropland cmass output file");
	declare_parameter("file_cmass_active_pasture", &file_cmass_active_pasture, 300, "Annual pasture cmass output file");
	declare_parameter("file_cmass_active_natural", &file_cmass_active_natural, 300, "Annual natural vegetation cmass output file");
	declare_parameter("file_cmass_active_forest", &file_cmass_active_forest, 300, "Annual managed forest cmass output file");
	declare_parameter("file_cmass_active_peatland", &file_cmass_active_peatland, 300, "Annual peatland cmass output file");
	declare_parameter("file_cmass_cropland", &file_cmass_cropland, 300, "Annual cropland cmass output file");
	declare_parameter("file_cmass_pasture", &file_cmass_pasture, 300, "Annual pasture cmass output file");
	declare_parameter("file_cmass_natural", &file_cmass_natural, 300, "Annual natural vegetation cmass output file");
	declare_parameter("file_cmass_forest", &file_cmass_forest, 300, "Annual managed forest cmass output file");
	declare_parameter("file_cmass_peatland", &file_cmass_peatland, 300, "Annual peatland cmass output file");
	declare_parameter("file_anpp_active_cropland", &file_anpp_active_cropland, 300, "Annual cropland NPP output file");
	declare_parameter("file_anpp_active_pasture", &file_anpp_active_pasture, 300, "Annual pasture NPP output file");
	declare_parameter("file_anpp_active_natural", &file_anpp_active_natural, 300, "Annual natural vegetation NPP output file");
	declare_parameter("file_anpp_active_forest", &file_anpp_active_forest, 300, "Annual managed forest NPP output file");
	declare_parameter("file_anpp_active_peatland", &file_anpp_active_peatland, 300, "Annual peatland NPP output file");
	declare_parameter("file_anpp_cropland", &file_anpp_cropland, 300, "Annual cropland NPP output file");
	declare_parameter("file_anpp_pasture", &file_anpp_pasture, 300, "Annual pasture NPP output file");
	declare_parameter("file_anpp_natural", &file_anpp_natural, 300, "Annual natural vegetation NPP output file");
	declare_parameter("file_anpp_forest", &file_anpp_forest, 300, "Annual managed forest NPP output file");
	declare_parameter("file_anpp_peatland", &file_anpp_peatland, 300, "Annual peatland NPP output file");
	declare_parameter("file_lai_cropland", &file_lai_cropland, 300, "Annual cropland LAI output file");
	declare_parameter("file_lai_pasture", &file_lai_pasture, 300, "Annual pasture LAI output file");
	declare_parameter("file_lai_natural", &file_lai_natural, 300, "Annual natural vegetation LAI output file");
	declare_parameter("file_lai_forest", &file_lai_forest, 300, "Annual managed forest LAI output file");
	declare_parameter("file_lai_peatland", &file_lai_peatland, 300, "Annual peatland LAI output file");

	declare_parameter("file_anpp_gridcell_natural", &file_anpp_gridcell_natural, 300, "Annual natural vegetation NPP per gridcell area output file");
	declare_parameter("file_anpp_gridcell_forest", &file_anpp_gridcell_forest, 300, "Annual managed forest NPP per gridcell area output file");
	declare_parameter("file_cmass_gridcell_natural", &file_cmass_gridcell_natural, 300, "Annual natural vegetation cmass per gridcell area output file");
	declare_parameter("file_cmass_gridcell_forest", &file_cmass_gridcell_forest, 300, "Annual managed forest cmass per gridcell area output file");
	declare_parameter("file_aaet_gridcell_natural", &file_aaet_gridcell_natural, 300, "Annual natural vegetation AET per gridcell area output file");
	declare_parameter("file_aaet_gridcell_forest", &file_aaet_gridcell_forest, 300, "Annual managed forest AET per gridcell area output file");
	declare_parameter("file_lai_gridcell_natural", &file_lai_gridcell_natural, 300, "Annual natural vegetation LAI per gridcell area output file");
	declare_parameter("file_lai_gridcell_forest", &file_lai_gridcell_forest, 300, "Annual managed forest LAI per gridcell area output file");

	declare_parameter("file_yield_active",&file_yield_active,300, "Crop yield output file");
	declare_parameter("file_yield1_active",&file_yield1_active,300,"Crop first yield output file");
	declare_parameter("file_yield2_active",&file_yield2_active,300,"Crop second yield output file");
	declare_parameter("file_sdate1",&file_sdate1,300,"Crop first sowing date output file");
	declare_parameter("file_sdate2",&file_sdate2,300,"Crop second sowing date output file");
	declare_parameter("file_hdate1",&file_hdate1,300,"Crop first harvest date output file");
	declare_parameter("file_hdate2",&file_hdate2,300,"Crop second harvest date output file");
	declare_parameter("file_lgp",&file_lgp,300,"Crop length of growing period output file");
	declare_parameter("file_phu",&file_phu,300,"Crop potential heat units output file");
	declare_parameter("file_fphu",&file_fphu,300,"Crop attained fraction of potential heat units output file");
	declare_parameter("file_fhi",&file_fhi,300,"Crop attained fraction of harvest index output file");
	declare_parameter("file_irrigation",&file_irrigation,300,"Crop irrigation output file");
	declare_parameter("file_seasonality",&file_seasonality,300,"Seasonality output file");
	declare_parameter("file_cflux_cropland", &file_cflux_cropland, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_pasture", &file_cflux_pasture, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_natural", &file_cflux_natural, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_forest", &file_cflux_forest, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_peatland", &file_cflux_peatland, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_forestry", &file_cflux_forestry, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_regrowth", &file_cflux_regrowth, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_cflux_primary", &file_cflux_primary, 300, "C fluxes per gridcell area output file");
	declare_parameter("file_dens_natural", &file_dens_natural, 300, "Natural vegetation tree density output file");
	declare_parameter("file_dens_forest", &file_dens_forest, 300, "Managed forest tree density output file");
	declare_parameter("file_cpool_cropland", &file_cpool_cropland, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_pasture", &file_cpool_pasture, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_natural", &file_cpool_natural, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_forest", &file_cpool_forest, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_peatland", &file_cpool_peatland, 300, "C pool per gridcell area file");
	declare_parameter("file_cpool_forestry", &file_cpool_forestry, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_regrowth", &file_cpool_regrowth, 300, "C pool per gridcell area output file");
	declare_parameter("file_cpool_primary", &file_cpool_primary, 300, "C pool per gridcell area output file");
	declare_parameter("file_nflux_cropland", &file_nflux_cropland, 300, "N fluxes per gridcell area output file");
	declare_parameter("file_nflux_pasture", &file_nflux_pasture, 300, "N fluxes per gridcell area output file");
	declare_parameter("file_nflux_natural", &file_nflux_natural, 300, "N fluxes per gridcell area output file");
	declare_parameter("file_nflux_forest", &file_nflux_forest, 300, "N fluxes per gridcell area output file");
	declare_parameter("file_nflux_peatland", &file_nflux_peatland, 300, "N fluxes per gridcell area output file");
	declare_parameter("file_npool_cropland", &file_npool_cropland, 300, "N pool per gridcell area output file");
	declare_parameter("file_npool_pasture", &file_npool_pasture, 300, "N pool per gridcell area output file");
	declare_parameter("file_npool_natural", &file_npool_natural, 300, "N pool per gridcell area output file");
	declare_parameter("file_npool_forest", &file_npool_forest, 300, "N pool per gridcell area output file");
	declare_parameter("file_npool_peatland", &file_npool_peatland, 300, "Soil N output per gridcell area file");
	declare_parameter("file_aaet_active_natural", &file_aaet_active_natural, 300, "Annual natural vegetation AET output file");
	declare_parameter("file_aaet_active_forest", &file_aaet_active_forest, 300, "Annual managed forest AET output file");
	declare_parameter("file_fpc_active_natural", &file_fpc_active_natural, 300, "Annual natural vegetation FPC output file");
	declare_parameter("file_fpc_active_forest", &file_fpc_active_forest, 300, "Annual managed forest FPC output file");
	declare_parameter("file_lai_active_natural", &file_lai_active_natural, 300, "Annual natural vegetation LAI output file");
	declare_parameter("file_lai_active_forest", &file_lai_active_forest, 300, "Annual managed forest LAI output file");
	declare_parameter("file_speciesdiam_natural", &file_speciesdiam_natural, 300, "Mean species diameter (cm)");
	declare_parameter("file_speciesdiam_forest", &file_speciesdiam_forest, 300, "Mean species diameter (cm)");
	declare_parameter("file_speciesheights_natural", &file_speciesheights_natural, 300, "Mean species height (m)");
	declare_parameter("file_speciesheights_forest", &file_speciesheights_forest, 300, "Mean species height (m)");

	declare_parameter("file_soil_nflux_cropland", &file_soil_nflux_cropland, 300, "Soil N fluxes output file");
	declare_parameter("file_soil_nflux_pasture", &file_soil_nflux_pasture, 300, "Soil N fluxes output file");
	declare_parameter("file_soil_nflux_natural", &file_soil_nflux_natural, 300, "Soil N fluxes output file");
	declare_parameter("file_soil_nflux_forest", &file_soil_nflux_forest, 300, "Soil N fluxes output file");
	declare_parameter("file_soil_n2o_flux", &file_soil_n2o_flux, 300, "Soil N2O fluxes per landcover output file, includes a Total.noBar (Total without BARREN contributions)");
	declare_parameter("file_frac_natural", &file_frac_natural, 300, "PFT area fraction output file");
	declare_parameter("file_frac_forest", &file_frac_forest, 300, "PFT area fraction output file");
	declare_parameter("file_agestruct_natural", &file_agestruct_natural, 300, "Age structure (tree density)");
	declare_parameter("file_agestruct_forest", &file_agestruct_forest, 300, "Age structure (tree density)");
	declare_parameter("file_diamstruct_natural", &file_diamstruct_natural, 300, "Diameter structure (tree density)");
	declare_parameter("file_diamstruct_forest", &file_diamstruct_forest, 300, "Diameter structure (tree density)");
	declare_parameter("file_diamstruct_cmass_natural", &file_diamstruct_cmass_natural, 300, "Diameter structure (cmass_wood)");
	declare_parameter("file_diamstruct_cmass_forest", &file_diamstruct_cmass_forest, 300, "Diameter structure (cmass_wood)");
	declare_parameter("file_diamstruct_cmass_leaf_natural", &file_diamstruct_cmass_leaf_natural, 300, "Diameter structure (cmass_wood_leaf)");
	declare_parameter("file_diamstruct_cmass_leaf_forest", &file_diamstruct_cmass_leaf_forest, 300, "Diameter structure (cmass_wood_leaf)");
	declare_parameter("file_diamstruct_cmass_wood_natural", &file_diamstruct_cmass_wood_natural, 300, "Diameter structure (cmass_wood)");
	declare_parameter("file_diamstruct_cmass_wood_forest", &file_diamstruct_cmass_wood_forest, 300, "Diameter structure (cmass_wood)");
	declare_parameter("file_diamstruct_cmass_froot_natural", &file_diamstruct_cmass_froot_natural, 300, "Diameter structure (cmass_root)");
	declare_parameter("file_diamstruct_cmass_froot_forest", &file_diamstruct_cmass_froot_forest, 300, "Diameter structure (cmass_root)");
	declare_parameter("file_diamstruct_cmass_wood_inc_natural", &file_diamstruct_cmass_wood_inc_natural, 300, "Diameter structure (cmass increment)");
	declare_parameter("file_diamstruct_cmass_wood_inc_forest", &file_diamstruct_cmass_wood_inc_forest, 300, "Diameter structure (cmass increment)");
	declare_parameter("file_diamstruct_ba_natural", &file_diamstruct_ba_natural, 300, "Diameter structure (basal area)");
	declare_parameter("file_diamstruct_ba_forest", &file_diamstruct_ba_forest, 300, "Diameter structure (basal area)");
	declare_parameter("file_agestruct_standage_natural", &file_agestruct_standage_natural, 300, "Age structure (stand age)");
	declare_parameter("file_agestruct_standage_forest", &file_agestruct_standage_forest, 300, "Age structure (stand age)");
	declare_parameter("file_agestruct_standage_cropland", &file_agestruct_standage_cropland, 300, "Age structure (stand age)");
	declare_parameter("file_agestruct_standage_pasture", &file_agestruct_standage_pasture, 300, "Age structure (stand age)");

	declare_parameter("file_anpp_sts", &file_anpp_sts, 300, "stand type anpp output file");
	declare_parameter("file_lai_sts", &file_lai_sts, 300, "stand type lai output file");
	declare_parameter("file_cmass_sts", &file_cmass_sts, 300, "stand type cmass output file");
	declare_parameter("file_cmass_harv_killed_sts", &file_cmass_harv_killed_sts, 300,
		"stand type whole tree harvest cmass output file");
	declare_parameter("file_cmass_thin_killed_sts", &file_cmass_thin_killed_sts, 300,
		"stand type whole tree thinning cmass output file");
	declare_parameter("file_cmass_clearcut_killed_sts", &file_cmass_clearcut_killed_sts, 300,
		"stand type whole tree non-satellite clearcut disturbance cmass output file");
	declare_parameter("file_cmass_fire_killed_sts", &file_cmass_fire_killed_sts, 300,
		"stand type whole tree prognostic fire disturbance cmass output file");
	declare_parameter("file_cmass_wind_killed_sts", &file_cmass_wind_killed_sts, 300,
		"stand type whole tree prognostic wind disturbance cmass output file");
	declare_parameter("file_cmass_barkbeetle_killed_sts", &file_cmass_barkbeetle_killed_sts, 300,
		"stand type whole tree prognostic bark beetle disturbance cmass output file");
	declare_parameter("file_cmass_dist_default_sts", &file_cmass_dist_default_sts, 300,
		"stand type whole tree killed in default patch-destroying disturbance cmass output file");
	declare_parameter("file_cmass_dist_prescribed_sts", &file_cmass_dist_prescribed_sts, 300,
		"stand type whole tree killed in prescribed patch-destroying disturbance cmass output file");
	declare_parameter("file_cmass_dist_clearcut_killed_sts", &file_cmass_dist_clearcut_killed_sts, 300,
		"stand type whole tree prescribed (satellite) clearcut disturbance cmass output file");
	declare_parameter("file_cmass_dist_fire_killed_sts", &file_cmass_dist_fire_killed_sts, 300,
		"stand type whole tree prescribed fire disturbance cmass output file");
	declare_parameter("file_cmass_dist_wind_killed_sts", &file_cmass_dist_wind_killed_sts, 300,
		"stand type whole tree prescribed wind disturbance cmass output file");
	declare_parameter("file_cmass_dist_barkbeetle_killed_sts", &file_cmass_dist_barkbeetle_killed_sts, 300,
		"stand type whole tree prescribed bark beetle disturbance cmass output file");
	declare_parameter("file_agpp_sts", &file_agpp_sts, 300, "stand type agpp output file");

	declare_parameter("file_anpp_tree_sts", &file_anpp_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_agpp_tree_sts", &file_agpp_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_aautresp_tree_sts", &file_aautresp_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_aautresp_root_tree_sts", &file_aautresp_root_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_aautresp_sap_tree_sts", &file_aautresp_sap_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_aautresp_leaf_tree_sts", &file_aautresp_leaf_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_aautresp_growth_tree_sts", &file_aautresp_growth_tree_sts, 300, "stand type tree lai output file");

	declare_parameter("file_lai_tree_sts", &file_lai_tree_sts, 300, "stand type tree lai output file");
	declare_parameter("file_cmass_tree_sts", &file_cmass_tree_sts, 300, "stand type tree cmass output file");
	declare_parameter("file_cmass_tree_mort_sts", &file_cmass_tree_mort_sts, 300,
		"stand type cmass of trees killed by mortality output file");

	declare_parameter("file_cmass_wood_sts", &file_cmass_wood_sts, 300, "stand type wood cmass output file");
	declare_parameter("file_cmass_wood_harv_killed_sts", &file_cmass_wood_harv_killed_sts, 300, "stand type gross wood harvest (standard plus prescribed) cmass output file");
	declare_parameter("file_cmass_wood_harv_sts", &file_cmass_wood_harv_sts, 300,
		"stand type net wood harvest (standard plus prescribed) cmass output file");
	declare_parameter("file_cmass_wood_harv_toprod_sts", &file_cmass_wood_harv_toprod_sts, 300,
		"stand type wood harvest (standard plus prescribed) product cmass output file");
	declare_parameter("file_cmass_wood_thin_sts", &file_cmass_wood_thin_sts, 300,
		"stand type thinning net wood harvest cmass output file");
	declare_parameter("file_cmass_wood_clearcut_sts", &file_cmass_wood_clearcut_sts, 300,
		"stand type clearcut (standard plus prescribed) net wood harvest cmass output file");
	declare_parameter("file_cmass_wood_mort_sts", &file_cmass_wood_mort_sts, 300, "stand type cmass of trees killed by mortality output file");
	declare_parameter("file_cmass_wood_inc_sts", &file_cmass_wood_inc_sts, 300, "stand type wood productivity output file");

	declare_parameter("file_csoil_sts", &file_csoil_sts, 300, "stand type soil C output file");
	declare_parameter("file_clitter_sts", &file_clitter_sts, 300, "stand type litter C output file");
	declare_parameter("file_clitter_course_sts", &file_clitter_course_sts, 300, "stand type course litter C output file");
	declare_parameter("file_clitter_fine_sts", &file_clitter_fine_sts, 300, "stand type fine litter C output file");
	declare_parameter("file_rh_sts", &file_rh_sts, 300, "stand type heterotrophic respiration C output file");
	declare_parameter("file_rh_clitter_course_sts", &file_rh_clitter_course_sts, 300, "stand type course wood heterotrophic respiration C output file");
	declare_parameter("file_rh_clitter_fine_sts", &file_rh_clitter_fine_sts, 300, "stand type fine wood heterotrophic respiration C output file");
	declare_parameter("file_csink_sts", &file_csink_sts, 300, "stand type carbon sink output file");
	declare_parameter("file_luc_cflux_sts", &file_luc_cflux_sts, 300, "stand type carbon lost in LUC output file");
	declare_parameter("file_nsoil_sts", &file_nsoil_sts, 300, "stand type soil N output file");
	declare_parameter("file_nlitter_sts", &file_nlitter_sts, 300, "stand type litter N output file");
	declare_parameter("file_nlitter_course_sts", &file_nlitter_course_sts, 300, "stand type course litter N output file");
	declare_parameter("file_nlitter_fine_sts", &file_nlitter_fine_sts, 300, "stand type fine litter N output file");

	declare_parameter("file_transpiration_sts", &file_transpiration_sts, 300, "stand type transpiration output file");
	declare_parameter("file_runoff_sts", &file_runoff_sts, 300, "stand type runoff output file");
	declare_parameter("file_intercep_sts", &file_intercep_sts, 300, "stand type interception output file");
	declare_parameter("file_evap_sts", &file_evap_sts, 300, "stand type evaporation file");

	declare_parameter("file_dens_sts", &file_dens_sts, 300, "stand type tree density output file");
	declare_parameter("file_height_sts", &file_height_sts, 300, "stand type height output file");
	declare_parameter("file_diam_sts", &file_diam_sts, 300, "stand type tree mean diameter output file");
	declare_parameter("file_diam_g_sts", &file_diam_g_sts, 300, "stand type tree quadratic mean diameter output file");

	declare_parameter("file_cmass_plus10cm_sts", &file_cmass_plus10cm_sts, 300, "stand type cmass output file");
	declare_parameter("file_cmass_harv_killed_plus10cm_sts", &file_cmass_harv_killed_plus10cm_sts, 300,
		"stand type whole tree harvest cmass output file");
	declare_parameter("file_cmass_dist_clearcut_killed_plus10cm_sts", &file_cmass_dist_clearcut_killed_plus10cm_sts, 300,
		"stand type whole tree clearcut disturbance cmass output file");
	declare_parameter("file_cmass_dist_fire_killed_plus10cm_sts", &file_cmass_dist_fire_killed_plus10cm_sts, 300,
		"stand type whole tree fire disturbance cmass output file");
	declare_parameter("file_cmass_dist_wind_killed_plus10cm_sts", &file_cmass_dist_wind_killed_plus10cm_sts, 300,
		"stand type whole tree wind disturbance cmass output file");
	declare_parameter("file_cmass_dist_barkbeetle_killed_plus10cm_sts", &file_cmass_dist_barkbeetle_killed_plus10cm_sts, 300,
		"stand type whole tree bark beetle disturbance cmass output file");
	declare_parameter("file_dens_plus10cm_sts", &file_dens_plus10cm_sts, 300, "stand type tree density output file");
	declare_parameter("file_height_plus10cm_sts", &file_height_plus10cm_sts, 300, "stand type height output file");
	declare_parameter("file_diam_plus10cm_sts", &file_diam_plus10cm_sts, 300, "stand type tree quadratic mean diameter output file");
	declare_parameter("file_diam_g_plus10cm_sts", &file_diam_g_plus10cm_sts, 300, "stand type tree quadratic mean diameter output file");

	declare_parameter("file_cmass_wood_plus10cm_sts", &file_cmass_wood_plus10cm_sts, 300, "stand type wood productivity output file");
	declare_parameter("file_cmass_wood_harv_killed_plus10cm_sts", &file_cmass_wood_harv_killed_plus10cm_sts, 300, "stand type wood harvest cmass output file");
	declare_parameter("file_cmass_wood_inc_plus10cm_sts", &file_cmass_wood_inc_plus10cm_sts, 300, "stand type wood productivity output file");
	declare_parameter("file_cmass_wood_mort_plus10cm_sts", &file_cmass_wood_mort_plus10cm_sts, 300, "stand type wood productivity output file");

	declare_parameter("file_height_sts_patch_mean", &file_height_sts_patch_mean, 300, "stand type height output file");
	declare_parameter("file_diam_sts_patch_mean", &file_diam_sts_patch_mean, 300, "stand type tree quadratic mean diameter output file");
	declare_parameter("file_diam_g_sts_patch_mean", &file_diam_g_sts_patch_mean, 300, "stand type tree quadratic mean diameter output file");

	declare_parameter("file_height_plus10cm_sts_patch_mean", &file_height_plus10cm_sts_patch_mean, 300, "stand type height output file");
	declare_parameter("file_diam_plus10cm_sts_patch_mean", &file_diam_plus10cm_sts_patch_mean, 300, "stand type tree quadratic mean diameter output file");
	declare_parameter("file_diam_g_plus10cm_sts_patch_mean", &file_diam_g_plus10cm_sts_patch_mean, 300, "stand type tree quadratic mean diameter output file");

	declare_parameter("file_cutinterval_sts", &file_cutinterval_sts, 1000,
		"Mean latest cutting interval (patch age at year of clearcut) output file");
	declare_parameter("file_cutinterval_thisyear_sts", &file_cutinterval_thisyear_sts, 1000,
		"Mean stand type cut interval this year output file");

	declare_parameter("file_nstand_sts", &file_nstand_sts, 300, "stand type stand number output file");
	declare_parameter("file_npatch_sts", &file_npatch_sts, 300, "stand type patch number output file");

	declare_parameter("file_frac_sts", &file_frac_sts, 300, "stand type are fraction output file");

	declare_parameter("file_forest_cmass_harv_killed_active", &file_forest_cmass_harv_killed_active, 300,
		"Killed forest C biomass during wood harvest output file");

	declare_parameter("file_cmass_active_plus10cm", &file_cmass_active_plus10cm, 300, "C biomass output file");
	declare_parameter("file_forest_cmass_harv_killed_active_plus10cm", &file_forest_cmass_harv_killed_active_plus10cm, 300,
		"Killed forest C biomass during wood harvest output file");
	declare_parameter("file_dens_plus10cm", &file_dens_plus10cm, 300, "Tree density output file");
	declare_parameter("file_speciesheights_plus10cm", &file_speciesheights_plus10cm, 300, "Mean species heights");
	declare_parameter("file_speciesdiam_plus10cm", &file_speciesdiam_plus10cm, 300, "Mean species diameter (cm)");

	declare_parameter("file_forest_vegc", &file_forest_vegc, 300, "Forest vegetation output file");
	declare_parameter("file_forest_cflux_veg", &file_forest_cflux_veg, 300, "Forest C fluxes to and from vegetation output file");
	declare_parameter("file_forest_harvest", &file_forest_harvest, 300, "Forest harvest output file");
	declare_parameter("file_harvest_flux_luc", &file_harvest_flux_luc, 300,
		"Harvest output file for simulations with wood harvest modelled as luc, using eg. LUH2 input");

	declare_parameter("file_cflux_harvest", &file_cflux_harvest, 300, "Harvest flux output file");
	declare_parameter("file_harvest_volume", &file_harvest_volume, 300, "Harvest volume output file");

	declare_parameter("file_storm", &file_storm, 300, "Storm indicators output file");
	declare_parameter("file_storm_si_sts", &file_storm_si_sts, 300, "Stand type level storm sensitivity index");
	declare_parameter("file_storm_frac_pot_sts", &file_storm_frac_pot_sts, 300, "Stand type level potentially storm damaged fraction");
	declare_parameter("file_storm_frac_sts", &file_storm_frac_sts, 300, "Stand type level storm damaged fraction");
	declare_parameter("file_storm_dw_pot_sts", &file_storm_dw_pot_sts, 300, "Potential stand type level storm damaged wood biomass C");
	declare_parameter("file_storm_dw_sts", &file_storm_dw_sts, 300, "Stand type level storm damaged wood biomass C");
	declare_parameter("file_insect_dw_sts", &file_insect_dw_sts, 300, "Stand type level insect damaged wood biomass C");

	declare_parameter("file_cmort_storm_sts", &file_cmort_storm_sts, 300, "Stand type level storm damaged biomass C");
	declare_parameter("file_cmort_insect_sts", &file_cmort_insect_sts, 300, "Stand type level insect damaged biomass C");
	declare_parameter("file_cmort_fire_sts", &file_cmort_fire_sts, 300, "Stand type level fire damaged biomass C");
	declare_parameter("file_cmort_clearcut_sts", &file_cmort_clearcut_sts, 300, "Stand type level clearcut biomass C");

	//daily
	declare_parameter("file_daily_lai",&file_daily_lai,300,"Daily output.");
	declare_parameter("file_daily_npp",&file_daily_npp,300,"Daily output.");
	declare_parameter("file_daily_nmass",&file_daily_nmass,300,"Daily output.");
	declare_parameter("file_daily_ndemand",&file_daily_ndemand,300,"Daily output.");
	declare_parameter("file_daily_cmass",&file_daily_cmass,300,"Daily output.");
	declare_parameter("file_daily_cton",&file_daily_cton,300,"Daily output.");
	declare_parameter("file_daily_cmass_leaf",&file_daily_cmass_leaf,300,"Daily output.");
	declare_parameter("file_daily_nmass_leaf",&file_daily_nmass_leaf,300,"Daily output.");
	declare_parameter("file_daily_cmass_root",&file_daily_cmass_root,300,"Daily output.");
	declare_parameter("file_daily_nmass_root",&file_daily_nmass_root,300,"Daily output.");
	declare_parameter("file_daily_cmass_stem",&file_daily_cmass_stem,300,"Daily output.");
	declare_parameter("file_daily_nmass_stem",&file_daily_nmass_stem,300,"Daily output.");
	declare_parameter("file_daily_cmass_storage",&file_daily_cmass_storage,300,"Daily output.");
	declare_parameter("file_daily_nmass_storage",&file_daily_nmass_storage,300,"Daily output.");

	declare_parameter("file_daily_cmass_dead_leaf",&file_daily_cmass_dead_leaf,300,"Daily output.");
	declare_parameter("file_daily_nmass_dead_leaf",&file_daily_nmass_dead_leaf,300,"Daily output.");

	declare_parameter("file_daily_n_input_soil",&file_daily_n_input_soil,300,"Daily output.");
	declare_parameter("file_daily_avail_nmass_soil",&file_daily_avail_nmass_soil,300,"Daily output.");
	declare_parameter("file_daily_upper_wcont",&file_daily_upper_wcont,300,"Daily output.");
	declare_parameter("file_daily_lower_wcont",&file_daily_lower_wcont,300,"Daily output.");
	declare_parameter("file_daily_irrigation",&file_daily_irrigation,300,"Daily output.");

	declare_parameter("file_daily_nminleach",&file_daily_nminleach,300,"Daily output.");
	declare_parameter("file_daily_norgleach",&file_daily_norgleach,300,"Daily output.");
	declare_parameter("file_daily_nuptake",&file_daily_nuptake,300,"Daily output.");

	declare_parameter("file_daily_climate",&file_daily_climate,300,"Daily output.");

	declare_parameter("file_daily_fphu",&file_daily_fphu,300,"Daily DS output file"); //daglig ds

	if (ifnlim) {
		declare_parameter("file_daily_ds",&file_daily_ds,300,"Daily DS output file"); //daglig ds
		declare_parameter("file_daily_stem",&file_daily_stem,300,"Daily stem allocation output file");
		declare_parameter("file_daily_leaf",&file_daily_leaf,300,"Daily leaf allocation output file");
		declare_parameter("file_daily_root",&file_daily_root,300,"Daily root allocation output file");
		declare_parameter("file_daily_storage",&file_daily_storage,300,"Daily storage allocation output file");
	}

	print_anpp_stand = false;
	print_lai_stand = false;
	print_cmass_stand = false;
	print_cmass_mort_stand = false;
	print_cmass_wood_stand = false;
	print_cmass_wood_harv_stand = false;
	print_height_stand = false;
	print_diam_stand = false;
	print_dens_stand = false;
	print_agestruct_stand = false;
	print_diamstruct_stand = false;
	print_diamstruct_cmass_stand = false;

	declare_parameter("print_anpp_stand",&print_anpp_stand,
		"Whether to print pft anpp for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_lai_stand",&print_lai_stand,
		"Whether to print pft lai for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_cmass_stand",&print_cmass_stand,
		"Whether to print pft cmass for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_cmass_wood_stand",&print_cmass_wood_stand,
		"Whether to print pft cmass_wood for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_cmass_wood_harv_stand",&print_cmass_wood_harv_stand,
		"Whether to print harvested pft cmass_wood for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_cmass_mort_stand",&print_cmass_mort_stand,
		"Whether to print C lost in mortality for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_height_stand",&print_height_stand,
		"Whether to print pft height for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_diam_stand",&print_diam_stand,
		"Whether to print pft diameter for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_dens_stand",&print_dens_stand,
		"Whether to print pft density for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_agestruct_stand",&print_agestruct_stand,
		"Whether to print tree densities in age classes for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_diamstruct_stand",&print_diamstruct_stand,
		"Whether to tree densities in diameter classes for multiple stands within a stand type (except cropland) separately");
	declare_parameter("print_diamstruct_cmass_stand",&print_diamstruct_cmass_stand,
		"Whether to print cmass_potharv in diameter classes for multiple stands within a stand type (except cropland) separately");

	printstandtypes = false;
	print_cmass_pft_st = false;
	print_diamstruct_cmass_st = false;
	print_diamstruct_cmass_leaf_st = false;
	print_diamstruct_cmass_wood_st = false;
	print_diamstruct_cmass_froot_st = false;
	print_diamstruct_cmass_wood_inc_st = false;
	print_cmass_harv_killed_pft_st = false;
	print_transpiration_pft_st = false;

	declare_parameter("printstandtypes", &printstandtypes, "Whether stand type output enabled (0,1)");
	declare_parameter("print_cmass_pft_st",&print_cmass_pft_st,"Whether to print pft cmass for stand types (except cropland) separately");	// Not sure if cropland st:s are excluded: check !
	declare_parameter("print_diamstruct_cmass_st",&print_diamstruct_cmass_st,
		"Whether to print pft cmass in diameter classes for stand types (except cropland) separately");
	declare_parameter("print_diamstruct_cmass_leaf_st",&print_diamstruct_cmass_leaf_st,
		"Whether to print pft leaf cmass in diameter classes for stand types (except cropland) separately");
	declare_parameter("print_diamstruct_cmass_wood_st",&print_diamstruct_cmass_wood_st,
		"Whether to print pft wood cmass in diameter classes for stand types (except cropland) separately");
	declare_parameter("print_diamstruct_cmass_froot_st",&print_diamstruct_cmass_froot_st,
		"Whether to print pft fine root cmass in diameter classes for stand types (except cropland) separately");
	declare_parameter("print_diamstruct_cmass_wood_inc_st",&print_diamstruct_cmass_wood_inc_st,
		"Whether to print pft wood increment cmass in diameter classes for stand types (except cropland) separately");
	declare_parameter("print_cmass_harv_killed_pft_st",&print_cmass_harv_killed_pft_st,
		"Whether to print pft cmass killed in harvest for stand types (except cropland) separately");
	declare_parameter("print_transpiration_pft_st",&print_transpiration_pft_st,"Whether to print pft transpiration for stand types (except cropland) separately");

	declare_parameter("print_standagestruct_st",&print_standagestruct_st,
		"Whether to print stand distribution in age classes for stand types separately");

	out_cmass_pft_st = NULL;
	out_cmass_harv_killed_pft_st = NULL;
	out_transpiration_pft_st = NULL;
	out_diamstruct_cmass_st = NULL;
	out_diamstruct_cmass_leaf_st = NULL;
	out_diamstruct_cmass_wood_st = NULL;
	out_diamstruct_cmass_froot_st = NULL;
	out_diamstruct_cmass_wood_inc_st = NULL;
	out_standagestruct_st = NULL;

	for(int id=0;id<MAXNUMBER_STANDS;id++) {
		out_anpp_stand[id] = NULL;
		out_lai_stand[id] = NULL;
		out_cmass_stand[id] = NULL;
		out_diam_stand[id] = NULL;
		out_height_stand[id] = NULL;
		out_dens_stand[id] = NULL;
		out_cmass_wood_stand[id] = NULL;
		out_cmass_wood_harv_stand[id] = NULL;
		out_cmass_mort_stand[id] = NULL;
		out_agestruct_stand[id] = NULL;
		out_diamstruct_stand[id] = NULL;
		out_diamstruct_cmass_stand[id] = NULL;
	}
	if (firemodel==BLAZE) {
		declare_parameter("file_daily_blaze",     &file_daily_blaze,     300,"Daily BLAZE output file");
        	//declare_parameter("file_single_fire",     &file_single_fire,     300,"Single fire output file");
	}

	// From EC-Earth
	declare_parameter("file_monthly_mlambda", &file_monthly_mlambda, 300, "Monthly lambda output file");
}

MiscOutput::~MiscOutput() {

	if(printstandtypes) {
		if(out_cmass_pft_st)
			delete[] out_cmass_pft_st;
		if(out_diamstruct_cmass_st)
			delete[] out_diamstruct_cmass_st;
		if(out_diamstruct_cmass_leaf_st)
			delete[] out_diamstruct_cmass_leaf_st;
		if(out_diamstruct_cmass_wood_st)
			delete[] out_diamstruct_cmass_wood_st;
		if(out_diamstruct_cmass_froot_st)
			delete[] out_diamstruct_cmass_froot_st;
		if(out_diamstruct_cmass_wood_inc_st)
			delete[] out_diamstruct_cmass_wood_inc_st;
		if(out_transpiration_pft_st)
			delete[] out_transpiration_pft_st;
		if(out_cmass_harv_killed_pft_st)
			delete[] out_cmass_harv_killed_pft_st;
		if(out_standagestruct_st)
			delete[] out_standagestruct_st;
	}
}

/// Help function to print forest structure header columns
std::vector<std::string> get_structure_string(const char* type) {

	char buffer2[300]={'\0'};
	std::vector<std::string> struct_vect;

	if(!strcmp(type, "age")) {

		for(int i=0;i<NFOREST_STRUCTURAL_CLASSES;i++) {

			// age structure
			int age_from, age_to;

			if(i < (NFOREST_STRUCTURAL_CLASSES - 1)) {			// 1-10,...,291-300
				age_from = i*10+1;
				age_to = i*10+10;
			}
			else if(i<NFOREST_STRUCTURAL_CLASSES) {			// >300y, not used in column header
				age_from = i*10+1;
				age_to = 1500;
			}

			if(i == (NFOREST_STRUCTURAL_CLASSES - 1)) {
				sprintf(buffer2, "gt300");
			}
			else {
				sprintf(buffer2, "%d_%d", age_from, age_to);
			}
			struct_vect.push_back(buffer2);					
		}
	}
	else if (!strcmp(type, "diam")) {

		for(int i=0;i<NFOREST_STRUCTURAL_CLASSES;i++) {

			// diameter structure
			double diam_from, diam_to;

			if(i < (NFOREST_STRUCTURAL_CLASSES - 1)) {			//  1-5,...,135-150
				diam_from = (double)(i*5);
				diam_to = (double)(i*5+5);
			}
			else if(i<NFOREST_STRUCTURAL_CLASSES) {			// >150cm	// not used in column header
				diam_from = i*5+1;
				diam_to = 1500;
			}

			if(i == (NFOREST_STRUCTURAL_CLASSES - 1)) {
				sprintf(buffer2, "gt150");				// >150cm
			}
			else {
				sprintf(buffer2, "%.0f_%.0f", diam_from, diam_to);
			}
			struct_vect.push_back(buffer2);
		}
	}

	return struct_vect;
}

/// Define all output tables and their formats
void MiscOutput::init() {
	
	define_output_tables();
}

/// Specify all columns in all output tables
/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void MiscOutput::define_output_tables() {
	// create a vector with the pft names
	std::vector<std::string> pfts;

	// create a vector with the crop pft names
	std::vector<std::string> crop_pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft=pftlist.getobj();

		pfts.push_back((char*)pft.name);

		if (pft.landcover==CROPLAND)
			crop_pfts.push_back((char*)pft.name);

		pftlist.nextobj();
	}

	// create a vector with the stand type names
	std::vector<std::string> sts;

	stlist.firstobj();
	while (stlist.isobj) {
		 StandType& st = stlist.getobj();
		 sts.push_back((char*)st.name);
		 stlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if(run_landcover) {
		const char* landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum",
				"Forest_sum", "Natural_sum", "Peatland_sum", "Barren_sum"};
		for (int i=0; i<NLANDCOVERTYPES; i++) {
			if (run[i]) {
				landcovers.push_back(landcover_string[i]);
			}
		}
	}

	// Create the month columns
	ColumnDescriptors month_columns;
	ColumnDescriptors month_columns_wide;
	xtring months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
	for (int i = 0; i < 12; i++) {
		month_columns      += ColumnDescriptor(months[i], 8,  3);
		month_columns_wide += ColumnDescriptor(months[i], 10, 3);
	}

	// Create the columns for each output file

	// CMASS
	ColumnDescriptors cmass_columns;
	cmass_columns += ColumnDescriptors(pfts,               8, 3);
	ColumnDescriptors cmass_active_columns = cmass_columns;
	cmass_columns += ColumnDescriptor("Total",             8, 3);
	ColumnDescriptors cmass_columns_lc = cmass_columns;
	cmass_columns += ColumnDescriptors(landcovers,        13, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;
	ColumnDescriptors anpp_columns_lc = cmass_columns_lc;

	// AET
	ColumnDescriptors aaet_columns_lc;
	aaet_columns_lc += ColumnDescriptors(pfts,                8, 2);
	aaet_columns_lc += ColumnDescriptor("Total",              8, 2);

	// DENS
	ColumnDescriptors dens_columns;
	dens_columns += ColumnDescriptors(pfts,                8, 4);
	dens_columns += ColumnDescriptor("Total",              8, 4);
	ColumnDescriptors dens_columns_lc = dens_columns;
	dens_columns += ColumnDescriptors(landcovers,         13, 4);

	// CFLUX
	ColumnDescriptors cflux_columns;
	cflux_columns += ColumnDescriptor("Veg",					 8, 3);
	cflux_columns += ColumnDescriptor("Repr",					 8, 3);
	cflux_columns += ColumnDescriptor("Soil",					 8, 3);
	cflux_columns += ColumnDescriptor("Fire",					10, 5);
	cflux_columns += ColumnDescriptor("Est",					 8, 3);
	cflux_columns += ColumnDescriptor("Seed",					 8, 3);
	cflux_columns += ColumnDescriptor("Harvest",				 9, 5);
	cflux_columns += ColumnDescriptor("LU_ch",					 9, 5);
	cflux_columns += ColumnDescriptor("Slow_h",					 9, 5);
	cflux_columns += ColumnDescriptor("NEE",					10, 5);

	// C FLUXES FOR FORESTS WITH LUC HISTORY (E.G. LUH2)
	ColumnDescriptors harv_columns_luc;
	harv_columns_luc += ColumnDescriptor("forharv_gross",				14, 5);
	harv_columns_luc += ColumnDescriptor("forharv_toprod",				15, 5);
	harv_columns_luc += ColumnDescriptor("clearing_gross",				15, 5);
	harv_columns_luc += ColumnDescriptor("clearing_toprod",				16, 5);
	harv_columns_luc += ColumnDescriptor("lucharv_gross",				14, 5);
	harv_columns_luc += ColumnDescriptor("lucharv_toprod",				15, 5);
	harv_columns_luc += ColumnDescriptor("fromprod",					 9, 5);
	harv_columns_luc += ColumnDescriptor("prod_balance",				13, 5);
	harv_columns_luc += ColumnDescriptor("harv_balance",				13, 5);

	ColumnDescriptors cflux_columns_for_regr;
	cflux_columns_for_regr += ColumnDescriptor("Veg",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("Repr",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("Soil",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("Fire",					 8, 5);
	cflux_columns_for_regr += ColumnDescriptor("Est",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("DOC",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("Seed",					 8, 3);
	cflux_columns_for_regr += ColumnDescriptor("Harvest",				 9, 5);
	cflux_columns_for_regr += ColumnDescriptor("NEE-LU",				10, 5);

	// FOREST WOOD HARVEST
	ColumnDescriptors harv_columns_active;
	harv_columns_active += ColumnDescriptors(pfts,								10, 5);

	ColumnDescriptors forest_harv_columns;
	forest_harv_columns += ColumnDescriptor("lucharv_totvegC",		16, 5);
	forest_harv_columns += ColumnDescriptor("lucharv_stemC",		14, 5);
	forest_harv_columns += ColumnDescriptor("lucharv_toprod",       15, 5);
	forest_harv_columns += ColumnDescriptor("lucharv_toflux",       15, 5);
	forest_harv_columns += ColumnDescriptor("lucharv_tolitt",       15, 5);

	forest_harv_columns += ColumnDescriptor("forharv_totvegC",		16, 5);
	forest_harv_columns += ColumnDescriptor("forharv_stemC",		14, 5);
	forest_harv_columns += ColumnDescriptor("forharv_toprod",       15, 5);
	forest_harv_columns += ColumnDescriptor("forharv_toflux",       15, 5);
	forest_harv_columns += ColumnDescriptor("forharv_tolitt",       15, 5);

	forest_harv_columns += ColumnDescriptor("totharv_totvegC",		16, 5);
	forest_harv_columns += ColumnDescriptor("totharv_stemC",		14, 5);
	forest_harv_columns += ColumnDescriptor("totharv_toprod",       15, 5);
	forest_harv_columns += ColumnDescriptor("totharv_toflux",       15, 5);
	forest_harv_columns += ColumnDescriptor("totharv_tolitt",       15, 5);

	// FOREST VEGC
	ColumnDescriptors forest_vegc_columns;
	forest_vegc_columns += ColumnDescriptor("nat_vegC",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("nat_stemC",				10, 5);
	forest_vegc_columns += ColumnDescriptor("nat_prod",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("nat_fuel",					 9, 5);

	forest_vegc_columns += ColumnDescriptor("for_vegC",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("for_stemC",				10, 5);
	forest_vegc_columns += ColumnDescriptor("for_prod",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("for_fuel",					 9, 5);

	forest_vegc_columns += ColumnDescriptor("tot_vegC",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("tot_stemC",				10, 5);
	forest_vegc_columns += ColumnDescriptor("tot_prod",					 9, 5);
	forest_vegc_columns += ColumnDescriptor("tot_fuel",					 9, 5);

	// FOREST CFLUX_VEG
	ColumnDescriptors forest_cflux_veg_columns;
	forest_cflux_veg_columns += ColumnDescriptor("nat_npp",				 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_harvC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_mortC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_fireC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_estC",			 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_distC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_turnoverC",		14, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_reprC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_clonedC",			12, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_cflux",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("nat_NAI",				10, 5);

	forest_cflux_veg_columns += ColumnDescriptor("for_npp",				 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_harvC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_mortC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_fireC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_estC",			 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_distC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_turnoverC",		14, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_reprC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_clonedC",			12, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_cflux",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("for_NAI",				10, 5);

	forest_cflux_veg_columns += ColumnDescriptor("tot_npp",				 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_harvC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_mortC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_fireC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_estC",			 9, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_distC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_turnoverC",		14, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_reprC",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_clonedC",			12, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_cflux",			10, 5);
	forest_cflux_veg_columns += ColumnDescriptor("tot_NAI",				10, 5);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",              8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",         8, 3);
		cpool_columns += ColumnDescriptor("SoilfC",        8, 3);
		cpool_columns += ColumnDescriptor("SoilsC",        8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LitterC",       8, 3);
		cpool_columns += ColumnDescriptor("SoilC",         8, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC",   10, 3);
	}
	cpool_columns += ColumnDescriptor("Total",            10, 3);

	// C POOLS FOR FORESTS WITH LUC HISTORY (E.G. LUH2)
	ColumnDescriptors cpool_columns_for_regr;
	cpool_columns_for_regr += ColumnDescriptor("VegC",              8, 3);
	cpool_columns_for_regr += ColumnDescriptor("LitterC",			8, 3);
	cpool_columns_for_regr += ColumnDescriptor("SoilC",				8, 3);
	cpool_columns_for_regr += ColumnDescriptor("Total-prod",	   11, 3);

	//CROP YIELD
	ColumnDescriptors crop_columns;
	crop_columns += ColumnDescriptors(crop_pfts,           8, 3);

	//CROP SDATE & HDATE
	ColumnDescriptors date_columns;
	date_columns += ColumnDescriptors(crop_pfts,           8, 0);

	//IRRIGATION
	ColumnDescriptors irrigation_columns;
	irrigation_columns += ColumnDescriptor("Total",       10, 3);

	//SEASONALITY
	ColumnDescriptors seasonality_columns;
	seasonality_columns += ColumnDescriptor("Seasonal",   10, 0);
	seasonality_columns += ColumnDescriptor("V_temp",     10, 3);
	seasonality_columns += ColumnDescriptor("V_prec",     10, 3);
	seasonality_columns += ColumnDescriptor("temp_min",   10, 1);  
	seasonality_columns += ColumnDescriptor("temp_mean",  10, 1);
	seasonality_columns += ColumnDescriptor("temp_max",   10, 1);
	seasonality_columns += ColumnDescriptor("mtemp_max",  10, 1);
	seasonality_columns += ColumnDescriptor("temp_seas",  10, 0);
	seasonality_columns += ColumnDescriptor("gdd0",       10, 0);
	seasonality_columns += ColumnDescriptor("gdd5",       10, 0);
	seasonality_columns += ColumnDescriptor("prec_min",   10, 2);
	seasonality_columns += ColumnDescriptor("prec",       10, 1);
	seasonality_columns += ColumnDescriptor("prec_range", 12, 0);
	seasonality_columns += ColumnDescriptor("insolation", 12, 4);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts,      8, 2);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",              9, 4);
	npool_columns += ColumnDescriptor("LitterN",           9, 4);
	npool_columns += ColumnDescriptor("SoilN",             9, 4);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN",    10, 4);
	}

	npool_columns += ColumnDescriptor("Total",            10, 4);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("dep",               8, 2);
	nflux_columns += ColumnDescriptor("fix",               8, 2);
	nflux_columns += ColumnDescriptor("fert",              8, 2);
	nflux_columns += ColumnDescriptor("flux",              8, 2);
	nflux_columns += ColumnDescriptor("leach",             8, 2);
	if (run_landcover) {
		nflux_columns += ColumnDescriptor("seed",          8, 2);
		nflux_columns += ColumnDescriptor("harvest",       8, 2);
		nflux_columns += ColumnDescriptor("LU_ch",         8, 3);
		nflux_columns += ColumnDescriptor("Slow_h",        8, 3);
	}
	nflux_columns += ColumnDescriptor("NEE",               8, 2);
	// SOIL N TRANSFORMATION - fluxes
	ColumnDescriptors soil_nflux_columns;
	soil_nflux_columns += ColumnDescriptor("NH3",  12, 6);
	soil_nflux_columns += ColumnDescriptor("NO",    9, 3);
	soil_nflux_columns += ColumnDescriptor("N2O",  12, 6);
	soil_nflux_columns += ColumnDescriptor("N2",    9, 3);

	ColumnDescriptors soil_n2o_flux_columns;
	soil_n2o_flux_columns += ColumnDescriptors(landcovers, 13, 5);
	soil_n2o_flux_columns += ColumnDescriptor("Total.noBar", 13, 5); // Total without BARREN contributions
	soil_n2o_flux_columns += ColumnDescriptor("Total", 10, 5);

	// Storm
	ColumnDescriptors storm_columns;
	storm_columns += ColumnDescriptor("SI", 9, 3);
	storm_columns += ColumnDescriptor("DamFracPot", 11, 3);
	storm_columns += ColumnDescriptor("DamFrac", 9, 3);
	storm_columns += ColumnDescriptor("DamWoodCPot", 12, 3);
	storm_columns += ColumnDescriptor("DamWoodC", 9, 3);
	storm_columns += ColumnDescriptor("InsDamWoodC", 12, 6);
	storm_columns += ColumnDescriptor("InsectPI", 9, 3);

	ColumnDescriptors cflux_harvest_columns;
	cflux_harvest_columns += ColumnDescriptor("cc", 9, 5);
	cflux_harvest_columns += ColumnDescriptor("thin", 9, 5);
	cflux_harvest_columns += ColumnDescriptor("cc_presc", 9, 5);
	cflux_harvest_columns += ColumnDescriptor("wind_presc", 9, 5);
	cflux_harvest_columns += ColumnDescriptor("bb_presc", 9, 5);

	// ST
	ColumnDescriptors st_columns;
	st_columns += ColumnDescriptors(sts,            15, 3);
	ColumnDescriptors st_dens_columns;
	st_dens_columns += ColumnDescriptors(sts,       15, 5);	// Increased to 5 decimals
	ColumnDescriptors st_columns_age;
	st_columns_age += ColumnDescriptors(sts,       15, 1);
	ColumnDescriptors st_columns_int;
	st_columns_int += ColumnDescriptors(sts,       15, 0);
	ColumnDescriptors st_height_columns;
	st_height_columns += ColumnDescriptors(sts,    15, 2);

	ColumnDescriptors st_ext_columns = st_columns;
	st_ext_columns += ColumnDescriptor("Total", 10, 4);
	st_ext_columns += ColumnDescriptors(landcovers, 13, 3);
	ColumnDescriptors st_ext_dens_columns = st_dens_columns;
	st_ext_dens_columns += ColumnDescriptor("Total", 10, 5);
	st_ext_dens_columns += ColumnDescriptors(landcovers, 13, 5);
	ColumnDescriptors st_ext_height_columns = st_height_columns;
	st_ext_height_columns += ColumnDescriptor("Total", 10, 2);
	st_ext_height_columns += ColumnDescriptors(landcovers, 13, 2);
	ColumnDescriptors st_ext_runoff_columns = st_columns_age;
	st_ext_runoff_columns += ColumnDescriptor("Total", 10, 1);
	st_ext_runoff_columns += ColumnDescriptors(landcovers, 13, 1);

	ColumnDescriptors st_int_ext_columns = st_columns_int;
	st_columns_int += ColumnDescriptor("Total", 10, 0);
	st_columns_int += ColumnDescriptors(landcovers, 13, 0);

	// FOREST STRUCTURE OUTPUT
	ColumnDescriptors agestruct_columns;
	agestruct_columns += ColumnDescriptors(get_structure_string("age"),			       9, 2);
	ColumnDescriptors standagestruct_columns;
	standagestruct_columns += ColumnDescriptors(get_structure_string("age"),		   9, 3);
	ColumnDescriptors diamstruct_columns;
	diamstruct_columns += ColumnDescriptors(get_structure_string("diam"),              9, 2);
	ColumnDescriptors diamstruct_cmass_columns;
	diamstruct_cmass_columns += ColumnDescriptors(get_structure_string("diam"),        9, 3);

	ColumnDescriptors daily_climate_columns;
	daily_climate_columns += ColumnDescriptor("Temp",   12, 6);
	daily_climate_columns += ColumnDescriptor("Prec",   12, 6);
	daily_climate_columns += ColumnDescriptor("Rad",    14, 3);

	ColumnDescriptors daily_columns;
	daily_columns += ColumnDescriptors(crop_pfts, 13, 3);

	// BLAZE output
	ColumnDescriptors blaze_daily_columns;
	blaze_daily_columns += ColumnDescriptor("die_r",         6, 0);
	blaze_daily_columns += ColumnDescriptor("fireProb%",    16, 12);
	blaze_daily_columns += ColumnDescriptor("effBurnedArea",15, 9);
	blaze_daily_columns += ColumnDescriptor("kbdi",          6, 1);
	blaze_daily_columns += ColumnDescriptor("ffdi",          6, 0);
	blaze_daily_columns += ColumnDescriptor("max_pfli",      9, 0);
	blaze_daily_columns += ColumnDescriptor("mean_pfli",    10, 0);

	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_forest_cmass_harv_killed_active,  file_forest_cmass_harv_killed_active,     harv_columns_active);

	create_output_table(out_cmass_active_plus10cm, file_cmass_active_plus10cm, cmass_active_columns);
	create_output_table(out_forest_cmass_harv_killed_active_plus10cm, file_forest_cmass_harv_killed_active_plus10cm, harv_columns_active);
	create_output_table(out_dens_plus10cm, file_dens_plus10cm, dens_columns);
	create_output_table(out_speciesheights_plus10cm, file_speciesheights_plus10cm, speciesheights_columns);
	create_output_table(out_speciesdiam_plus10cm, file_speciesdiam_plus10cm, speciesheights_columns);

	create_output_table(out_cmass_active_cropland, file_cmass_active_cropland, cmass_columns_lc);
	create_output_table(out_cmass_active_pasture,  file_cmass_active_pasture,  cmass_columns_lc);
	create_output_table(out_cmass_active_natural,  file_cmass_active_natural,  cmass_columns_lc);
	create_output_table(out_cmass_active_forest,   file_cmass_active_forest,   cmass_columns_lc);
	create_output_table(out_cmass_active_peatland,	file_cmass_active_peatland, cmass_columns_lc);
	create_output_table(out_anpp_active_cropland,  file_anpp_active_cropland,  anpp_columns_lc);
	create_output_table(out_anpp_active_pasture,   file_anpp_active_pasture,   anpp_columns_lc);
	create_output_table(out_anpp_active_natural,   file_anpp_active_natural,   anpp_columns_lc);
	create_output_table(out_anpp_active_forest,    file_anpp_active_forest,    anpp_columns_lc);
	create_output_table(out_anpp_active_peatland, file_anpp_active_peatland, anpp_columns_lc);

	create_output_table(out_cmass_cropland,	file_cmass_cropland, cmass_columns_lc);
	create_output_table(out_cmass_pasture,	file_cmass_pasture,  cmass_columns_lc);
	create_output_table(out_cmass_natural,	file_cmass_natural,  cmass_columns_lc);
	create_output_table(out_cmass_forest,		file_cmass_forest,   cmass_columns_lc);
	create_output_table(out_cmass_peatland,	file_cmass_peatland, cmass_columns_lc);
	create_output_table(out_anpp_cropland,	file_anpp_cropland,  anpp_columns_lc);
	create_output_table(out_anpp_pasture,		file_anpp_pasture,   anpp_columns_lc);
	create_output_table(out_anpp_natural,		file_anpp_natural,   anpp_columns_lc);
	create_output_table(out_anpp_forest,		file_anpp_forest,    anpp_columns_lc);
	create_output_table(out_anpp_peatland,	file_anpp_peatland,  anpp_columns_lc);
	create_output_table(out_lai_cropland,		file_lai_cropland,   dens_columns_lc);
	create_output_table(out_lai_pasture,		file_lai_pasture,    dens_columns_lc);
	create_output_table(out_lai_natural,		file_lai_natural,    dens_columns_lc);
	create_output_table(out_lai_forest,			file_lai_forest,     dens_columns_lc);
	create_output_table(out_lai_peatland,		file_lai_peatland,   dens_columns_lc);

	create_output_table(out_dens_natural,   file_dens_natural,   dens_columns_lc);
	create_output_table(out_dens_forest,    file_dens_forest,    dens_columns_lc);
	create_output_table(out_cflux_cropland, file_cflux_cropland, cflux_columns);
	create_output_table(out_cflux_pasture,  file_cflux_pasture,  cflux_columns);
	create_output_table(out_cflux_natural,  file_cflux_natural,  cflux_columns);
	create_output_table(out_cflux_forest,	file_cflux_forest,   cflux_columns);
	create_output_table(out_cflux_peatland, file_cflux_peatland, cflux_columns);
	create_output_table(out_cpool_cropland, file_cpool_cropland, cpool_columns);
	create_output_table(out_cpool_pasture,  file_cpool_pasture,  cpool_columns);
	create_output_table(out_cpool_natural,  file_cpool_natural,  cpool_columns);
	create_output_table(out_cpool_forest,	file_cpool_forest,   cpool_columns);
	create_output_table(out_cpool_peatland, file_cpool_peatland, cpool_columns);
	create_output_table(out_aaet_active_natural,   file_aaet_active_natural,   aaet_columns_lc);
	create_output_table(out_aaet_active_forest,    file_aaet_active_forest,    aaet_columns_lc);
	create_output_table(out_lai_active_natural,    file_lai_active_natural,    dens_columns_lc);
	create_output_table(out_lai_active_forest,     file_lai_active_forest,     dens_columns_lc);
	create_output_table(out_fpc_active_natural,    file_fpc_active_natural,    cmass_columns_lc);
	create_output_table(out_fpc_active_forest,     file_fpc_active_forest,     cmass_columns_lc);
	create_output_table(out_speciesdiam_natural,   file_speciesdiam_natural,   speciesheights_columns);
	create_output_table(out_speciesdiam_forest,    file_speciesdiam_forest,    speciesheights_columns);
	create_output_table(out_speciesheights_natural, file_speciesheights_natural, speciesheights_columns);
	create_output_table(out_speciesheights_forest,  file_speciesheights_forest,  speciesheights_columns);

	create_output_table(out_aaet_gridcell_natural,   file_aaet_gridcell_natural,   aaet_columns_lc);
	create_output_table(out_aaet_gridcell_forest,    file_aaet_gridcell_forest,    aaet_columns_lc);
	create_output_table(out_lai_gridcell_natural,    file_lai_gridcell_natural,    dens_columns_lc);
	create_output_table(out_lai_gridcell_forest,     file_lai_gridcell_forest,     dens_columns_lc);
	create_output_table(out_cmass_gridcell_natural,  file_cmass_gridcell_natural,  dens_columns_lc);
	create_output_table(out_cmass_gridcell_forest,   file_cmass_gridcell_forest,   dens_columns_lc);
	create_output_table(out_anpp_gridcell_natural,   file_anpp_gridcell_natural,   dens_columns_lc);
	create_output_table(out_anpp_gridcell_forest,    file_anpp_gridcell_forest,    dens_columns_lc);

	if (run_landcover && run[CROPLAND]) {
		create_output_table(out_yield_active,      file_yield_active,          crop_columns);
		create_output_table(out_yield1_active,     file_yield1_active,         crop_columns);
		create_output_table(out_yield2_active,     file_yield2_active,         crop_columns);
		create_output_table(out_sdate1,     file_sdate1,         date_columns);
		create_output_table(out_sdate2,     file_sdate2,         date_columns);
		create_output_table(out_hdate1,     file_hdate1,         date_columns);
		create_output_table(out_hdate2,     file_hdate2,         date_columns);
		create_output_table(out_lgp,        file_lgp,            date_columns);
		create_output_table(out_phu,        file_phu,            date_columns);
		create_output_table(out_fphu,       file_fphu,           crop_columns);
		create_output_table(out_fhi,        file_fhi,            crop_columns);
	}

    create_output_table(out_seasonality,file_seasonality,    seasonality_columns);

	if(run_landcover)
		create_output_table(out_irrigation, file_irrigation,     irrigation_columns); 

	create_output_table(out_npool_cropland, file_npool_cropland, npool_columns);
	create_output_table(out_npool_pasture,  file_npool_pasture,  npool_columns);
	create_output_table(out_npool_natural,  file_npool_natural,  npool_columns);
	create_output_table(out_npool_forest,	file_npool_forest,   npool_columns);
	create_output_table(out_npool_peatland, file_npool_peatland, npool_columns);
	create_output_table(out_nflux_cropland, file_nflux_cropland, nflux_columns);
	create_output_table(out_nflux_pasture,  file_nflux_pasture,  nflux_columns);
	create_output_table(out_nflux_natural,  file_nflux_natural,  nflux_columns);
	create_output_table(out_nflux_forest,	file_nflux_forest,   nflux_columns);
	create_output_table(out_nflux_peatland, file_nflux_peatland, nflux_columns);
	create_output_table(out_soil_nflux_cropland, file_soil_nflux_cropland, soil_nflux_columns);
	create_output_table(out_soil_nflux_pasture,  file_soil_nflux_pasture,  soil_nflux_columns);
	create_output_table(out_soil_nflux_natural,  file_soil_nflux_natural,  soil_nflux_columns);
	create_output_table(out_soil_nflux_forest,	file_soil_nflux_forest,   soil_nflux_columns);

	create_output_table(out_frac_natural,	file_frac_natural,	dens_columns_lc);
	create_output_table(out_frac_forest,	file_frac_forest,	dens_columns_lc);

	create_output_table(out_soil_n2o_flux,	file_soil_n2o_flux,   soil_n2o_flux_columns);
	// TODO		create_output_table(out_nflux_peatland, file_nflux_peatland, nflux_columns);

	create_output_table(out_anpp_sts,					file_anpp_sts,					st_ext_columns);
	create_output_table(out_agpp_sts,					file_agpp_sts,					st_ext_columns);
	create_output_table(out_anpp_tree_sts,					file_anpp_tree_sts,					st_ext_columns);
	create_output_table(out_agpp_tree_sts,					file_agpp_tree_sts,					st_ext_columns);
	create_output_table(out_aautresp_tree_sts,					file_aautresp_tree_sts,					st_ext_columns);
	create_output_table(out_aautresp_root_tree_sts,					file_aautresp_root_tree_sts,					st_ext_columns);
	create_output_table(out_aautresp_sap_tree_sts,					file_aautresp_sap_tree_sts,					st_ext_columns);
	create_output_table(out_aautresp_leaf_tree_sts,					file_aautresp_leaf_tree_sts,					st_ext_columns);
	create_output_table(out_aautresp_growth_tree_sts,					file_aautresp_growth_tree_sts,					st_ext_columns);
	create_output_table(out_lai_sts,					file_lai_sts,					st_ext_columns);
	create_output_table(out_lai_tree_sts,				file_lai_tree_sts,				st_ext_columns);
	create_output_table(out_cmass_sts,					file_cmass_sts,					st_ext_columns);
	create_output_table(out_cmass_tree_sts,				file_cmass_tree_sts,			st_ext_columns);
	create_output_table(out_cmass_tree_mort_sts,		file_cmass_tree_mort_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_harv_killed_sts,		file_cmass_harv_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_thin_killed_sts,		file_cmass_thin_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_clearcut_killed_sts,	file_cmass_clearcut_killed_sts,	st_ext_dens_columns);
	create_output_table(out_cmass_fire_killed_sts,		file_cmass_fire_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_wind_killed_sts,		file_cmass_wind_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_barkbeetle_killed_sts,		file_cmass_barkbeetle_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_clearcut_killed_sts,	file_cmass_dist_clearcut_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_fire_killed_sts,	file_cmass_dist_fire_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_default_sts,				file_cmass_dist_default_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_prescribed_killed_sts,				file_cmass_dist_prescribed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_wind_killed_sts,	file_cmass_dist_wind_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_dist_barkbeetle_killed_sts,		file_cmass_dist_barkbeetle_killed_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_wood_sts,				file_cmass_wood_sts,			st_ext_columns);
	create_output_table(out_cmass_wood_harv_killed_sts, file_cmass_wood_harv_killed_sts, st_ext_dens_columns);
	create_output_table(out_cmass_wood_harv_sts,		file_cmass_wood_harv_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_wood_harv_toprod_sts, file_cmass_wood_harv_toprod_sts,st_ext_dens_columns);
	create_output_table(out_cmass_wood_thin_sts,		file_cmass_wood_thin_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_wood_clearcut_sts,	file_cmass_wood_clearcut_sts,   st_ext_dens_columns);
	create_output_table(out_cmass_wood_inc_sts,			file_cmass_wood_inc_sts,		st_ext_dens_columns);
	create_output_table(out_cmass_wood_mort_sts,		file_cmass_wood_mort_sts,		st_ext_dens_columns);
	create_output_table(out_diam_sts,					file_diam_sts,					st_ext_height_columns);	
	create_output_table(out_diam_g_sts,					file_diam_g_sts,				st_ext_height_columns);	// 2 decimals only in this branch
	create_output_table(out_cutinterval_sts,			file_cutinterval_sts,			st_columns_age);
	create_output_table(out_cutinterval_thisyear_sts,	file_cutinterval_thisyear_sts,  st_columns_age);
	create_output_table(out_dens_sts,					file_dens_sts,					st_ext_dens_columns);
	create_output_table(out_height_sts,					file_height_sts,				st_ext_height_columns);
	create_output_table(out_csoil_sts,					file_csoil_sts,					st_ext_columns);
	create_output_table(out_clitter_sts,				file_clitter_sts,				st_ext_columns);
	create_output_table(out_clitter_course_sts,			file_clitter_course_sts,		st_ext_columns);
	create_output_table(out_clitter_fine_sts,			file_clitter_fine_sts,			st_ext_columns);
	create_output_table(out_rh_sts,						file_rh_sts,					st_ext_columns);
	create_output_table(out_rh_clitter_course_sts,		file_rh_clitter_course_sts,		st_ext_columns);
	create_output_table(out_rh_clitter_fine_sts,		file_rh_clitter_fine_sts,		st_ext_columns);
	create_output_table(out_csink_sts,					file_csink_sts,					st_ext_dens_columns);
	create_output_table(out_luc_cflux_sts,				file_luc_cflux_sts,				st_ext_columns);
	create_output_table(out_nsoil_sts,					file_nsoil_sts,					st_ext_columns);
	create_output_table(out_nlitter_sts,				file_nlitter_sts,				st_ext_dens_columns);
	create_output_table(out_nlitter_course_sts,			file_nlitter_course_sts,		st_ext_dens_columns);
	create_output_table(out_nlitter_fine_sts,			file_nlitter_fine_sts,			st_ext_dens_columns);
	create_output_table(out_nstand_sts,					file_nstand_sts,				st_columns_int);
	create_output_table(out_npatch_sts,					file_npatch_sts,				st_columns_int);
	create_output_table(out_frac_sts,					file_frac_sts,					st_ext_dens_columns);

	create_output_table(out_storm_si_sts,				file_storm_si_sts,				st_ext_columns);
	create_output_table(out_storm_frac_sts,				file_storm_frac_sts,			st_ext_columns);
	create_output_table(out_storm_frac_pot_sts,			file_storm_frac_pot_sts,		st_ext_columns);
	create_output_table(out_storm_dw_sts,				file_storm_dw_sts,				st_ext_dens_columns);
	create_output_table(out_storm_dw_pot_sts,			file_storm_dw_pot_sts,			st_ext_columns);
	create_output_table(out_insect_dw_sts,				file_insect_dw_sts,				st_ext_dens_columns);

	create_output_table(out_cmort_storm_sts,			file_cmort_storm_sts,			st_ext_dens_columns);
	create_output_table(out_cmort_insect_sts,			file_cmort_insect_sts,			st_ext_dens_columns);
	create_output_table(out_cmort_fire_sts,				file_cmort_fire_sts,			st_ext_dens_columns);
	create_output_table(out_cmort_clearcut_sts,			file_cmort_clearcut_sts,		st_ext_dens_columns);

	create_output_table(out_cmass_plus10cm_sts, file_cmass_plus10cm_sts, st_ext_columns);
	create_output_table(out_cmass_harv_killed_plus10cm_sts, file_cmass_harv_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_dist_clearcut_killed_plus10cm_sts, file_cmass_dist_clearcut_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_dist_fire_killed_plus10cm_sts, file_cmass_dist_fire_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_dist_wind_killed_plus10cm_sts, file_cmass_dist_wind_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_dist_barkbeetle_killed_plus10cm_sts, file_cmass_dist_barkbeetle_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_wood_plus10cm_sts, file_cmass_wood_plus10cm_sts, st_ext_columns);
	create_output_table(out_cmass_wood_harv_killed_plus10cm_sts, file_cmass_wood_harv_killed_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_wood_mort_plus10cm_sts, file_cmass_wood_mort_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_cmass_wood_inc_plus10cm_sts, file_cmass_wood_inc_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_dens_plus10cm_sts, file_dens_plus10cm_sts, st_ext_dens_columns);
	create_output_table(out_height_plus10cm_sts, file_height_plus10cm_sts, st_ext_height_columns);
	create_output_table(out_diam_plus10cm_sts, file_diam_plus10cm_sts, st_ext_height_columns);
	create_output_table(out_diam_g_plus10cm_sts, file_diam_g_plus10cm_sts, st_ext_height_columns);

	create_output_table(out_transpiration_sts,			file_transpiration_sts,         st_ext_runoff_columns);
	create_output_table(out_runoff_sts,					file_runoff_sts,				st_ext_runoff_columns);
	create_output_table(out_intercep_sts,				file_intercep_sts,				st_ext_runoff_columns);
	create_output_table(out_evap_sts,					file_evap_sts,					st_ext_runoff_columns);

	create_output_table(out_height_plus10cm_sts_patch_mean, file_height_plus10cm_sts_patch_mean, st_ext_height_columns);
	create_output_table(out_diam_plus10cm_sts_patch_mean, file_diam_plus10cm_sts_patch_mean, st_ext_height_columns);
	create_output_table(out_diam_g_plus10cm_sts_patch_mean, file_diam_g_plus10cm_sts_patch_mean, st_ext_height_columns);

	create_output_table(out_height_sts_patch_mean, file_height_sts_patch_mean, st_ext_height_columns);
	create_output_table(out_diam_sts_patch_mean, file_diam_sts_patch_mean, st_ext_height_columns);
	create_output_table(out_diam_g_sts_patch_mean, file_diam_g_sts_patch_mean, st_ext_height_columns);

	create_output_table(out_forest_harvest,				file_forest_harvest,			forest_harv_columns);
	create_output_table(out_forest_vegc,				file_forest_vegc,				forest_vegc_columns);
	create_output_table(out_forest_cflux_veg,			file_forest_cflux_veg,			forest_cflux_veg_columns);

	create_output_table(out_cflux_forestry,				file_cflux_forestry,			cflux_columns_for_regr);
	create_output_table(out_cflux_regrowth,				file_cflux_regrowth,			cflux_columns_for_regr);
	create_output_table(out_cflux_primary,				file_cflux_primary,				cflux_columns_for_regr);
	create_output_table(out_cpool_forestry,				file_cpool_forestry,			cpool_columns_for_regr);
	create_output_table(out_cpool_regrowth,				file_cpool_regrowth,			cpool_columns_for_regr);
	create_output_table(out_cpool_primary,				file_cpool_primary,				cpool_columns_for_regr);
	create_output_table(out_harvest_flux_luc,			file_harvest_flux_luc,			harv_columns_luc);

	create_output_table(out_agestruct_natural, file_agestruct_natural, agestruct_columns);
	create_output_table(out_agestruct_forest, file_agestruct_forest, agestruct_columns);
	create_output_table(out_diamstruct_natural, file_diamstruct_natural, diamstruct_columns);
	create_output_table(out_diamstruct_forest, file_diamstruct_forest, diamstruct_columns);
	create_output_table(out_diamstruct_cmass_natural, file_diamstruct_cmass_natural, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_forest, file_diamstruct_cmass_forest, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_leaf_natural, file_diamstruct_cmass_leaf_natural, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_leaf_forest, file_diamstruct_cmass_leaf_forest, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_wood_natural, file_diamstruct_cmass_wood_natural, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_wood_forest, file_diamstruct_cmass_wood_forest, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_froot_natural, file_diamstruct_cmass_froot_natural, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_froot_forest, file_diamstruct_cmass_froot_forest, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_wood_inc_natural, file_diamstruct_cmass_wood_inc_natural, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_cmass_wood_inc_forest, file_diamstruct_cmass_wood_inc_forest, diamstruct_cmass_columns);
	create_output_table(out_diamstruct_ba_natural, file_diamstruct_ba_natural, diamstruct_columns);
	create_output_table(out_diamstruct_ba_forest, file_diamstruct_ba_forest, diamstruct_columns);
	create_output_table(out_agestruct_standage_natural, file_agestruct_standage_natural, standagestruct_columns);
	create_output_table(out_agestruct_standage_forest, file_agestruct_standage_forest, standagestruct_columns);
	create_output_table(out_agestruct_standage_cropland, file_agestruct_standage_cropland, standagestruct_columns);
	create_output_table(out_agestruct_standage_pasture, file_agestruct_standage_pasture, standagestruct_columns);

	create_output_table(out_cflux_harvest, file_cflux_harvest, cflux_harvest_columns);
	create_output_table(out_harvest_volume, file_harvest_volume, cflux_harvest_columns);

	create_output_table(out_storm, file_storm, storm_columns);

	if(printstandtypes) {
		char dirname[200]={'\0'};
		// Put files in separate directory only on Windows for now.
#ifdef _MSC_VER
		strcpy(dirname, "st_output/");
		std::filesystem::create_directories(dirname);
#endif
		out_cmass_pft_st = new Table[nst];
		out_cmass_harv_killed_pft_st = new Table[nst];
		out_transpiration_pft_st = new Table[nst];
		out_diamstruct_cmass_st = new Table[nst];
		out_diamstruct_cmass_leaf_st = new Table[nst];
		out_diamstruct_cmass_wood_st = new Table[nst];
		out_diamstruct_cmass_froot_st = new Table[nst];
		out_diamstruct_cmass_wood_inc_st = new Table[nst];
		out_standagestruct_st = new Table[nst];

		for(int i=0;i<nst;i++) {
			StandType& st = stlist[i];

			char outfilename[100]={'\0'}, buffer[50]={'\0'};

			strcat(buffer, "_pft_st.out");

			if(print_cmass_pft_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "cmass_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, buffer);
				create_output_table(out_cmass_pft_st[i], outfilename, anpp_columns_lc);
			}
			if(print_transpiration_pft_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "transpiration_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, buffer);
				create_output_table(out_transpiration_pft_st[i], outfilename, anpp_columns_lc);
			}
			if(print_diamstruct_cmass_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "diamstruct_cmass_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, "_st.out");
				create_output_table(out_diamstruct_cmass_st[i], outfilename, diamstruct_cmass_columns);
			}
			if(print_diamstruct_cmass_leaf_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "diamstruct_cmass_leaf_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, "_st.out");
				create_output_table(out_diamstruct_cmass_leaf_st[i], outfilename, diamstruct_cmass_columns);
			}
			if(print_diamstruct_cmass_wood_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "diamstruct_cmass_wood_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, "_st.out");
				create_output_table(out_diamstruct_cmass_wood_st[i], outfilename, diamstruct_cmass_columns);
			}
			if(print_diamstruct_cmass_froot_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "diamstruct_cmass_froot_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, "_st.out");
				create_output_table(out_diamstruct_cmass_froot_st[i], outfilename, diamstruct_cmass_columns);
			}
			if(print_diamstruct_cmass_wood_inc_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "diamstruct_cmass_wood_inc_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, "_st.out");
				create_output_table(out_diamstruct_cmass_wood_inc_st[i], outfilename, diamstruct_cmass_columns);
			}
			if(print_cmass_harv_killed_pft_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "cmass_harv_killed_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, buffer);
				create_output_table(out_cmass_harv_killed_pft_st[i], outfilename, anpp_columns_lc);
			}
			if(print_standagestruct_st) {
				outfilename[0] = '\0';
				strcpy(outfilename, dirname);
				strcat(outfilename, "standagestruct_");
				strcat(outfilename, (char*)st.name);
				strcat(outfilename, buffer);
				create_output_table(out_standagestruct_st[i], outfilename, standagestruct_columns);
			}
		}
	}

	// *** DAILY OUTPUT VARIABLES ***

	create_output_table(out_daily_lai,					file_daily_lai,					daily_columns);
	create_output_table(out_daily_npp,					file_daily_npp,					daily_columns);
	create_output_table(out_daily_ndemand,				file_daily_ndemand,				daily_columns);
	create_output_table(out_daily_nmass,				file_daily_nmass,				daily_columns);
	create_output_table(out_daily_cmass,				file_daily_cmass,				daily_columns);
	create_output_table(out_daily_nmass_leaf,			file_daily_nmass_leaf,			daily_columns);
	create_output_table(out_daily_cmass_leaf,			file_daily_cmass_leaf,			daily_columns);
	create_output_table(out_daily_nmass_root,			file_daily_nmass_root,			daily_columns);
	create_output_table(out_daily_cmass_root,			file_daily_cmass_root,			daily_columns);
	create_output_table(out_daily_nmass_stem,			file_daily_nmass_stem,			daily_columns);
	create_output_table(out_daily_cmass_stem,			file_daily_cmass_stem,			daily_columns);
	create_output_table(out_daily_nmass_storage,        file_daily_nmass_storage,       daily_columns);
	create_output_table(out_daily_cmass_storage,        file_daily_cmass_storage,       daily_columns);
	create_output_table(out_daily_nmass_dead_leaf,      file_daily_nmass_dead_leaf,     daily_columns);
	create_output_table(out_daily_cmass_dead_leaf,      file_daily_cmass_dead_leaf,     daily_columns);
	create_output_table(out_daily_n_input_soil,         file_daily_n_input_soil,        daily_columns);
	create_output_table(out_daily_avail_nmass_soil,     file_daily_avail_nmass_soil,    daily_columns);

	create_output_table(out_daily_upper_wcont,			file_daily_upper_wcont,         daily_columns);
	create_output_table(out_daily_lower_wcont,			file_daily_lower_wcont,         daily_columns);
	create_output_table(out_daily_irrigation,			file_daily_irrigation,			daily_columns);

	create_output_table(out_daily_climate,				file_daily_climate,				daily_climate_columns);

	create_output_table(out_daily_cton,					file_daily_cton,				daily_columns);

	create_output_table(out_daily_nminleach,			file_daily_nminleach,			daily_columns);
	create_output_table(out_daily_norgleach,			file_daily_norgleach,			daily_columns);
	create_output_table(out_daily_nuptake,				file_daily_nuptake,				daily_columns);

	if (ifnlim) {
		create_output_table(out_daily_ds,				file_daily_ds,					daily_columns);
		create_output_table(out_daily_fphu,				file_daily_fphu,				daily_columns);
		create_output_table(out_daily_stem,				file_daily_stem,				daily_columns);
		create_output_table(out_daily_leaf,				file_daily_leaf,				daily_columns);
		create_output_table(out_daily_root,				file_daily_root,				daily_columns);
		create_output_table(out_daily_storage,			file_daily_storage,				daily_columns);
	}
	if (firemodel==BLAZE) {
		create_output_table(out_daily_blaze,      file_daily_blaze,      blaze_daily_columns);
	}

	ColumnDescriptors mpft_columns; // EC-Earth
	mpft_columns += ColumnDescriptors(pfts, 9, 4);

	create_output_table(out_mlambda, file_monthly_mlambda, mpft_columns); // EC-Earth
}

/// Local analogue of OutputRows::add_value for restricting output
/** Use to restrict output to specified range of years
  * (or other user-specified limitation)
  *
  * If only yearly output between, say 1961 and 1990 is requred, use:
  *  if (date.get_calendar_year() >= 1961 && date.get_calendar_year() <= 1990)
  *  (assuming the input module has set the first calendar year in the date object)	
  */
void outlimit_misc(OutputRows& out, const Table& table, double d) {

	if (date.year >= nyear_spinup)
		out.add_value(table, d);
}

/// Output of simulation results at the end of each year
/** Output of simulation results at the end of each year, or for specific years in
  * the simulation of each stand or grid cell. 
  * This function does not have to provide any information to the framework.
  *
  * Restrict output to specific years in the local helper function outlimit_misc().
  *
  * Changes in the structure of CommonOutput::outannual() should be mirrored here.
  */
void MiscOutput::outannual(Gridcell& gridcell) {

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();

	double flux_veg_regrowth, flux_repr_regrowth, flux_soil_regrowth, flux_fire_regrowth, flux_est_regrowth, 
		flux_seed_regrowth, flux_charvest_regrowth;
	double flux_veg_forestry, flux_repr_forestry, flux_soil_forestry, flux_fire_forestry, flux_est_forestry, 
		flux_seed_forestry, flux_charvest_forestry;
	double flux_veg_primary, flux_repr_primary, flux_soil_primary, flux_fire_primary, flux_est_primary, 
		flux_seed_primary, flux_charvest_primary;
	double surfsoillitterc_regrowth, surfsoillitterc_forestry, surfsoillitterc_primary;
	double cwdc_regrowth, cwdc_forestry, cwdc_primary;
	double centuryc_regrowth, centuryc_forestry, centuryc_primary;

	Landcover& lc = gridcell.landcover;
	// The OutputRows object manages the next row of output for each
	// output table
	OutputRows out(output_channel, lon, lat, date.get_calendar_year());

	double mean_standpft_cmass_landscape=0.0;
	double mean_standpft_anpp_landscape=0.0;
	double mean_standpft_lai_landscape=0.0;

	double landcover_cmass[NLANDCOVERTYPES]={0.0};
	double landcover_nmass[NLANDCOVERTYPES]={0.0};
	double landcover_clitter[NLANDCOVERTYPES]={0.0};
	double landcover_nlitter[NLANDCOVERTYPES]={0.0};
	double landcover_anpp[NLANDCOVERTYPES]={0.0};
	double landcover_fpc[NLANDCOVERTYPES]={0.0};
	double landcover_aaet[NLANDCOVERTYPES]={0.0};
	double landcover_lai[NLANDCOVERTYPES]={0.0};
	double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};
	double landcover_cmass_harv_killed[NLANDCOVERTYPES]={0.0};

	double landcover_cmass_plus10cm[NLANDCOVERTYPES] = { 0.0 };
	double landcover_cmass_harv_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
	double landcover_densindiv_total_plus10cm[NLANDCOVERTYPES] = { 0.0 };

	double mean_standpft_cmass_harv_killed=0.0;

	double mean_standpft_cmass_plus10cm = 0.0;
	double mean_standpft_cmass_harv_killed_plus10cm = 0.0;
	double mean_standpft_densindiv_total_plus10cm = 0.0;
	double mean_standpft_heightindiv_total_plus10cm = 0.0;
	double mean_standpft_diamindiv_total_plus10cm = 0.0;

	double mean_standpft_anpp_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_cmass_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_densindiv_total_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_anpp_landscape_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_cmass_landscape_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_lai_landscape_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_aaet_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_lai_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_fpc_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_diamindiv_total_lc[NLANDCOVERTYPES]={0.0};
	double mean_standpft_heightindiv_total_lc[NLANDCOVERTYPES]={0.0};

	double gridcellpft_anpp_lc[NLANDCOVERTYPES]={0.0};
	double gridcellpft_cmass_lc[NLANDCOVERTYPES]={0.0};
	double gridcellpft_aaet_lc[NLANDCOVERTYPES]={0.0};
	double gridcellpft_lai_lc[NLANDCOVERTYPES]={0.0};

	double gridcellpft_cmass_harv_killed=0.0;	// NOT USED ?

	double cmass_gridcell=0.0;
	double anpp_gridcell=0.0;
	double lai_gridcell=0.0;

	double cmass_gridcell_forestry=0.0;
	double cmass_gridcell_regrowth=0.0;
	double cmass_gridcell_primary=0.0;
	double clitter_gridcell_forestry=0.0;
	double clitter_gridcell_regrowth=0.0;
	double clitter_gridcell_primary=0.0;
	double c_org_leach_gridcell_regrowth=0.0;
	double c_org_leach_gridcell_forestry=0.0;
	double c_org_leach_gridcell_primary=0.0;

	double irrigation_gridcell=0.0;

	double cmass_harv_killed_gridcell=0.0;

	double cmass_gridcell_plus10cm = 0.0;
	double cmass_harv_killed_gridcell_plus10cm = 0.0;
	double dens_gridcell_plus10cm = 0.0;

	double standpft_cmass=0.0;
	double standpft_cmass_wood=0.0;
	double standpft_cmass_wood_potharv=0.0;
	double standpft_cmass_wood_potharv_products=0.0;
	double standpft_cmass_potfuel=0.0;
	double standpft_diam_g=0.0;
	double standpft_cmass_wood_harv=0.0;
	double standpft_cmass_wood_harv_toprod=0.0;
	double standpft_cmass_wood_clearcut=0.0;
	double standpft_cmass_harv_killed=0.0;
	double standpft_cmass_harv_tolitter=0.0;
	double standpft_cmass_mort=0.0;
	double standpft_cmass_fire=0.0;
	double standpft_cmass_dist=0.0;
	double standpft_cmass_leaf_root_turnover=0.0;
	double standpft_cmass_repr=0.0;
	double standpft_cmass_est=0.0;
	double standpft_nmass=0.0;
	double standpft_clitter=0.0;
	double standpft_nlitter=0.0;
	double standpft_anpp=0.0;
	double standpft_aaet=0.0;
	double standpft_lai=0.0;
	double standpft_fpc=0.0;
	double standpft_yield=0.0;
	double standpft_yield1=0.0;
	double standpft_yield2=0.0;
	double standpft_densindiv_total = 0.0;
	double standpft_heightindiv_total = 0.0;
	double standpft_diamindiv_total = 0.0;

	double standpft_cmass_plus10cm = 0.0;
	double standpft_densindiv_total_plus10cm = 0.0;
	double standpft_heightindiv_total_plus10cm = 0.0;
	double standpft_diamindiv_total_plus10cm = 0.0;
	double standpft_diam_g_plus10cm = 0.0;
	double standpft_cmass_harv_killed_plus10cm = 0.0;

	// from EC-Earth
	double pft_mlambda_gpp[12][MAXNUMBER_PFTS] = { 0.0 };
	double pft_mgpp[12][MAXNUMBER_PFTS] = { 0.0 };

	double* st_pft_cmass = new double[nst];
	double* st_total_cmass = new double[nst];
	double* st_pft_cmass_harv_killed = new double[nst];
	double* st_total_cmass_harv_killed = new double[nst];
	double* st_pft_aaet = new double[nst];
	double* st_total_aaet = new double[nst];
	
	for(int stid=0; stid<nst; stid++) {
		st_total_cmass[stid] = 0.0;
		st_total_cmass_harv_killed[stid] = 0.0;
		st_total_aaet[stid] = 0.0;
	}

	pftlist.firstobj();
	while (pftlist.isobj) {

		for(int stid=0; stid<nst; stid++) {
			st_pft_cmass[stid] = 0.0;
			st_pft_cmass_harv_killed[stid] = 0.0;
			st_pft_aaet[stid] = 0.0;
		}

		Pft& pft=pftlist.getobj();

		mean_standpft_anpp_landscape = 0.0;
		mean_standpft_cmass_landscape = 0.0;
		mean_standpft_lai_landscape = 0.0;

		mean_standpft_cmass_harv_killed = 0.0;

		mean_standpft_cmass_plus10cm = 0.0;
		mean_standpft_cmass_harv_killed_plus10cm = 0.0;
		mean_standpft_densindiv_total_plus10cm = 0.0;
		mean_standpft_heightindiv_total_plus10cm = 0.0;
		mean_standpft_diamindiv_total_plus10cm = 0.0;

		// Sum values across stands, patches and PFTs
		double mean_standpft_yield=0.0;
		double mean_standpft_yield1=0.0;
		double mean_standpft_yield2=0.0;

		for (int i=0; i<NLANDCOVERTYPES; i++) {
			mean_standpft_anpp_lc[i]=0.0;
			mean_standpft_cmass_lc[i]=0.0;
			mean_standpft_densindiv_total_lc[i]=0.0;
			mean_standpft_anpp_landscape_lc[i]=0.0;
			mean_standpft_cmass_landscape_lc[i]=0.0;
			mean_standpft_lai_landscape_lc[i]=0.0;

			mean_standpft_aaet_lc[i]=0.0;
			mean_standpft_lai_lc[i]=0.0;
			mean_standpft_fpc_lc[i]=0.0;
			mean_standpft_diamindiv_total_lc[i] = 0.0;
			mean_standpft_heightindiv_total_lc[i] = 0.0;

			gridcellpft_anpp_lc[i]=0.0;
			gridcellpft_cmass_lc[i]=0.0;
			gridcellpft_aaet_lc[i]=0.0;
			gridcellpft_lai_lc[i]=0.0;
		}

		// Determine area fraction of stands where this pft is active:
		double active_fraction = 0.0;
		double active_fraction_lc[NLANDCOVERTYPES]={0.0};

		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			if (stand.pft[pft.id].active) {
				active_fraction += stand.get_gridcell_fraction();
				active_fraction_lc[stand.landcover] += stand.get_gridcell_fraction();
			}

			++gc_itr;
		}

		// Loop through Stands
		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			Standpft& standpft=stand.pft[pft.id];

			// Sum values across patches and PFTs
			standpft_cmass=0.0;
			standpft_cmass_mort=0.0;
			standpft_cmass_fire=0.0;
			standpft_cmass_dist=0.0;
			standpft_cmass_leaf_root_turnover=0.0;
			standpft_cmass_repr=0.0;
			standpft_cmass_est=0.0;
			standpft_cmass_wood=0.0;
			standpft_cmass_wood_potharv=0.0;
			standpft_cmass_wood_potharv_products=0.0;
			standpft_cmass_potfuel=0.0;
			standpft_cmass_wood_harv=0.0;
			standpft_cmass_wood_harv_toprod=0.0;
			standpft_cmass_wood_clearcut=0.0;
			standpft_cmass_harv_killed=0.0;
			standpft_cmass_harv_tolitter=0.0;
			standpft_diam_g=0.0;
			standpft_nmass=0.0;
			standpft_clitter=0.0;
			standpft_nlitter=0.0;
			standpft_anpp=0.0;
			standpft_fpc=0.0;
			standpft_aaet=0.0;
			standpft_lai=0.0;
			standpft_yield=0.0;
			standpft_yield1=0.0;
			standpft_yield2=0.0;
			standpft_densindiv_total = 0.0;
			standpft_heightindiv_total = 0.0;
			standpft_diamindiv_total = 0.0;

			standpft_cmass_plus10cm = 0.0;
			standpft_densindiv_total_plus10cm = 0.0;
			standpft_heightindiv_total_plus10cm = 0.0;
			standpft_diamindiv_total_plus10cm = 0.0;
			standpft_diam_g_plus10cm = 0.0;
			standpft_cmass_harv_killed_plus10cm = 0.0;

			stand.firstobj();

			// Loop through Patches
			while (stand.isobj) {
				Patch& patch = stand.getobj();
				Patchpft& patchpft = patch.pft[pft.id];
				Vegetation& vegetation = patch.vegetation;

				// from EC-Earth:
				double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

				standpft_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);

				standpft_cmass_wood_harv += patchpft.cmass_wood_harv;
				standpft_cmass_wood_harv_toprod += patchpft.cmass_wood_harv_toprod;
				standpft_cmass_wood_clearcut += patchpft.cmass_wood_clearcut;
				standpft_cmass_harv_tolitter += patchpft.cmass_harv_tolitter;
				standpft_cmass_harv_killed += patchpft.cmass_harv_killed;
				standpft_cmass_mort += patchpft.cmass_mort;
				standpft_cmass_fire += patchpft.cmass_fire_killed;
				standpft_cmass_dist += patchpft.cmass_dist;
				standpft_cmass_leaf_root_turnover += patchpft.cmass_leaf_root_turnover;
				standpft_cmass_repr += patchpft.cmass_repr;
				standpft_cmass_est += patchpft.cmass_est;
				standpft_clitter += patchpft.cmass_litter_leaf + patchpft.cmass_litter_root + patchpft.cmass_litter_sap
					+ patchpft.cmass_litter_heart + patchpft.cmass_litter_repr;
				standpft_nlitter += patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap
					+ patchpft.nmass_litter_heart;

				standpft_cmass_harv_killed_plus10cm += patchpft.cmass_harv_killed_plus10cm;

				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv=vegetation.getobj();

					if (indiv.id!=-1 && indiv.alive && indiv.pft.id==pft.id) {
//					if (indiv.id!=-1 && indiv.alive && indiv.pft.id==pft.id && indiv.diam > 0.1) {	// TreeMort

						standpft_cmass += indiv.ccont();
						standpft_cmass_wood += indiv.cmass_wood();
						standpft_cmass_wood_potharv += check_harvest_cmass(indiv, true);
						standpft_cmass_wood_potharv_products += check_harvest_cmass(indiv, true, true);
						standpft_cmass_potfuel += (check_harvest_cmass(indiv) - check_harvest_cmass(indiv, true, true));
						standpft_nmass += indiv.ncont();
						standpft_fpc += indiv.fpc;
						standpft_aaet += indiv.aaet;
						standpft_lai += indiv.lai;

						if (pft.lifeform == TREE && indiv.diam >= 0.1) {	// TreeMort
							standpft_cmass_plus10cm += indiv.ccont();
							standpft_densindiv_total_plus10cm += indiv.densindiv;
							standpft_heightindiv_total_plus10cm += indiv.height * indiv.densindiv;
							standpft_diamindiv_total_plus10cm += indiv.diam * indiv.densindiv;
							standpft_diam_g_plus10cm += pow(indiv.diam, 2) * indiv.densindiv;
						}

						if (pft.landcover == CROPLAND) {
							standpft_yield += indiv.cropindiv->harv_yield;
							standpft_yield1 += indiv.cropindiv->yield_harvest[0];
							standpft_yield2 += indiv.cropindiv->yield_harvest[1];
						}
						else {

							if (vegmode==COHORT || vegmode==INDIVIDUAL) {
								if (pft.lifeform==TREE) {
									standpft_densindiv_total += indiv.densindiv; // indiv/m2
									standpft_diamindiv_total += indiv.diam * indiv.densindiv;
									standpft_heightindiv_total += indiv.height * indiv.densindiv;
									standpft_diam_g += pow(indiv.diam, 2) * indiv.densindiv;
								}
							}

							// from EC-Earth
							// Calculate monthly lambda
							for (int m = 0; m < 12; m++) {
								// pft_mlambda_gpp[m][pft.id] += max(0.0, indiv.mlambda_gpp[m]) * to_gridcell_average;
								pft_mgpp[m][pft.id] += max(0.0, indiv.mgpp[m]) * to_gridcell_average;
							}
						}
					}
					vegetation.nextobj();
				}

				stand.nextobj();
			} // end of patch loop

			standpft_cmass/=(double)stand.npatch();
			standpft_cmass_wood/=(double)stand.npatch();
			standpft_cmass_wood_potharv/=(double)stand.npatch();
			standpft_cmass_wood_potharv_products/=(double)stand.npatch();
			standpft_cmass_potfuel/=(double)stand.npatch();
			standpft_cmass_mort/=(double)stand.npatch();
			standpft_cmass_fire/=(double)stand.npatch();
			standpft_cmass_dist/=(double)stand.npatch();
			standpft_cmass_leaf_root_turnover/=(double)stand.npatch();
			standpft_cmass_repr/=(double)stand.npatch();
			standpft_cmass_est/=(double)stand.npatch();
			standpft_diam_g/=(double)stand.npatch();
			standpft_cmass_wood_harv/=(double)stand.npatch();
			standpft_cmass_wood_harv_toprod/=(double)stand.npatch();
			standpft_cmass_wood_clearcut/=(double)stand.npatch();
			standpft_cmass_harv_tolitter/=(double)stand.npatch();
			standpft_cmass_harv_killed/=(double)stand.npatch();
			standpft_nmass/=(double)stand.npatch();
			standpft_clitter/=(double)stand.npatch();
			standpft_nlitter/=(double)stand.npatch();
			standpft_anpp/=(double)stand.npatch();
			standpft_fpc/=(double)stand.npatch();
			standpft_aaet/=(double)stand.npatch();
			standpft_lai/=(double)stand.npatch();
			standpft_densindiv_total/=(double)stand.npatch();
			standpft_heightindiv_total/=(double)stand.npatch();
			standpft_diamindiv_total/=(double)stand.npatch();
			standpft_yield/=(double)stand.npatch();
			standpft_yield1/=(double)stand.npatch();
			standpft_yield2/=(double)stand.npatch();

			standpft_cmass_plus10cm /= (double)stand.npatch();
			standpft_densindiv_total_plus10cm /= (double)stand.npatch();
			standpft_heightindiv_total_plus10cm /= (double)stand.npatch();
			standpft_diamindiv_total_plus10cm /= (double)stand.npatch();
			standpft_diam_g_plus10cm /= (double)stand.npatch();
			standpft_cmass_harv_killed_plus10cm /= (double)stand.npatch();

			//Update landcover totals
			landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
			landcover_nmass[stand.landcover]+=standpft_nmass*stand.get_landcover_fraction();
			landcover_clitter[stand.landcover]+=standpft_clitter*stand.get_landcover_fraction();
			landcover_nlitter[stand.landcover]+=standpft_nlitter*stand.get_landcover_fraction();
			landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
			if(!pft.isintercropgrass) {
				landcover_fpc[stand.landcover]+=standpft_fpc*stand.get_landcover_fraction();
				landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
			}
			landcover_aaet[stand.landcover]+=standpft_aaet*stand.get_landcover_fraction();
			landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();
			landcover_cmass_harv_killed[stand.landcover]+=standpft_cmass_harv_killed*stand.get_landcover_fraction();

			landcover_cmass_plus10cm[stand.landcover] += standpft_cmass_plus10cm * stand.get_landcover_fraction();
			landcover_cmass_harv_killed_plus10cm[stand.landcover] += standpft_cmass_harv_killed_plus10cm * stand.get_landcover_fraction();
			landcover_densindiv_total_plus10cm[stand.landcover] += standpft_densindiv_total_plus10cm * stand.get_landcover_fraction();

			//Update pft means for active stands
			if(active_fraction) {
				mean_standpft_cmass_harv_killed += standpft_cmass_harv_killed * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_yield += standpft_yield * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_yield1 += standpft_yield1 * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_yield2 += standpft_yield2 * stand.get_gridcell_fraction() / active_fraction;

				mean_standpft_cmass_plus10cm += standpft_cmass_plus10cm * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_cmass_harv_killed_plus10cm += standpft_cmass_harv_killed_plus10cm * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_densindiv_total_plus10cm += standpft_densindiv_total_plus10cm * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_heightindiv_total_plus10cm += standpft_heightindiv_total_plus10cm * stand.get_gridcell_fraction() / active_fraction;
				mean_standpft_diamindiv_total_plus10cm += standpft_diamindiv_total_plus10cm * stand.get_gridcell_fraction() / active_fraction;

				//Update pft mean for active stands in landcover
				if(active_fraction_lc[stand.landcover]) {
					mean_standpft_anpp_lc[stand.landcover] +=
						standpft_anpp * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_cmass_lc[stand.landcover] +=
						standpft_cmass * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_densindiv_total_lc[stand.landcover] +=
						standpft_densindiv_total * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];

					mean_standpft_aaet_lc[stand.landcover] +=
						standpft_aaet * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_lai_lc[stand.landcover] +=
						standpft_lai * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_fpc_lc[stand.landcover] +=
						standpft_fpc * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_heightindiv_total_lc[stand.landcover] +=
						standpft_heightindiv_total * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
					mean_standpft_diamindiv_total_lc[stand.landcover]  +=
						standpft_diamindiv_total * stand.get_gridcell_fraction() / active_fraction_lc[stand.landcover];
				}
			}

			mean_standpft_cmass_landscape += standpft_cmass * stand.get_gridcell_fraction();
			mean_standpft_anpp_landscape += standpft_anpp * stand.get_gridcell_fraction();
			mean_standpft_lai_landscape += standpft_lai * stand.get_gridcell_fraction();
			mean_standpft_anpp_landscape_lc[stand.landcover] += standpft_anpp * stand.get_landcover_fraction();
			mean_standpft_cmass_landscape_lc[stand.landcover] += standpft_cmass * stand.get_landcover_fraction();
			mean_standpft_lai_landscape_lc[stand.landcover] += standpft_lai * stand.get_landcover_fraction();

			gridcellpft_cmass_harv_killed += standpft_cmass_harv_killed * stand.get_gridcell_fraction();	// NOT USED  ?

			gridcellpft_cmass_lc[stand.landcover] += standpft_cmass * stand.get_gridcell_fraction();
			gridcellpft_anpp_lc[stand.landcover] += standpft_anpp * stand.get_gridcell_fraction();
			gridcellpft_aaet_lc[stand.landcover] += standpft_aaet * stand.get_gridcell_fraction();
			gridcellpft_lai_lc[stand.landcover] += standpft_lai * stand.get_gridcell_fraction();

			//Update stand totals
			stand.anpp += standpft_anpp;
			stand.lai += standpft_lai;
			stand.cmass += standpft_cmass;
			stand.cmass_wood += standpft_cmass_wood;
			stand.cmass_wood_harv += standpft_cmass_wood_harv;
			stand.cmass_mort += standpft_cmass_mort;

			// Print per-stand pft values
			if (printseparatestands) {

				int id = stand.id;
				if(stand.id >= MAXNUMBER_STANDS)
					fail("Number of stands above limit, increase MAXNUMBER_STANDS for output of individual stands !\n\n");

				if(!out_anpp_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_anpp_stand[id][stand.stid], standpft_anpp);
				if(!out_lai_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_lai_stand[id][stand.stid], standpft_lai);
				if(!out_cmass_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_cmass_stand[id][stand.stid], standpft_cmass);
				if(!out_cmass_wood_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_cmass_wood_stand[id][stand.stid], standpft_cmass_wood);
				if(!out_cmass_wood_harv_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_cmass_wood_harv_stand[id][stand.stid], standpft_cmass_wood_harv);
				if(!out_cmass_mort_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_cmass_mort_stand[id][stand.stid], standpft_cmass_mort);

				double height = 0.0;
				double diam = 0.0;
				if (standpft_densindiv_total > 0.0) {
					height = standpft_heightindiv_total / standpft_densindiv_total;
					diam = standpft_diamindiv_total / standpft_densindiv_total;
				}
				if(!out_height_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_height_stand[id][stand.stid], height);
				if(!out_diam_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_diam_stand[id][stand.stid], diam * 100.0);	//diameter output in cm
				if(!out_dens_stand[id][stand.stid].invalid())
					outlimit_misc(out, out_dens_stand[id][stand.stid], standpft_densindiv_total);
			}

			//Update stand type totals
			StandType& st = stlist[stand.stid];
			Gridcellst& gcst = gridcell.st[stand.stid];

			st_pft_cmass[stand.stid] += standpft_cmass * stand.get_gridcell_fraction();
			st_total_cmass[stand.stid] += standpft_cmass * stand.get_gridcell_fraction();
			st_pft_cmass_harv_killed[stand.stid] += standpft_cmass_harv_killed * stand.get_gridcell_fraction();
			st_total_cmass_harv_killed[stand.stid] += standpft_cmass_harv_killed * stand.get_gridcell_fraction();
			st_pft_aaet[stand.stid] += standpft_aaet * stand.get_gridcell_fraction();
			st_total_aaet[stand.stid] += standpft_aaet * stand.get_gridcell_fraction();

			// Update gridcell totals
			double fraction_of_gridcell = stand.get_gridcell_fraction();
			cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
			anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
			lai_gridcell+=standpft_lai*fraction_of_gridcell;
			cmass_harv_killed_gridcell+=standpft_cmass_harv_killed*fraction_of_gridcell;

			cmass_gridcell_plus10cm += standpft_cmass_plus10cm * fraction_of_gridcell;
			cmass_harv_killed_gridcell_plus10cm += standpft_cmass_harv_killed_plus10cm * fraction_of_gridcell;
			dens_gridcell_plus10cm += standpft_densindiv_total_plus10cm * fraction_of_gridcell;

			if(stand.landcover == NATURAL || stand.landcover == FOREST) {
				if(stand.first_year) {
					if(stand.lc_origin == NATURAL) {
						cmass_gridcell_forestry+=standpft_cmass*fraction_of_gridcell;
						clitter_gridcell_forestry+=standpft_clitter*fraction_of_gridcell;
					}
					else {
						cmass_gridcell_regrowth+=standpft_cmass*fraction_of_gridcell;
						clitter_gridcell_regrowth+=standpft_clitter*fraction_of_gridcell;
					}
				}
				else {
					cmass_gridcell_primary+=standpft_cmass*fraction_of_gridcell;
					clitter_gridcell_primary+=standpft_clitter*fraction_of_gridcell;
				}
			}

			++gc_itr;
		}//End of loop through stands

		// Print to pft per stand type output
		if(printstandtypes) {
			for(int stid=0; stid<nst; stid++) {
				Gridcellst& gcst = gridcell.st[stid];
				if(gcst.frac) {
					if(!out_cmass_pft_st[stid].invalid())
						outlimit_misc(out, out_cmass_pft_st[stid], st_pft_cmass[stid] / gcst.frac);
					if(!out_cmass_harv_killed_pft_st[stid].invalid())
						outlimit_misc(out, out_cmass_harv_killed_pft_st[stid], st_pft_cmass_harv_killed[stid] / gcst.frac);
					if(!out_transpiration_pft_st[stid].invalid())
						outlimit_misc(out, out_transpiration_pft_st[stid], st_pft_aaet[stid] / gcst.frac);
				}
				else {
					if(!out_cmass_pft_st[stid].invalid())
						outlimit_misc(out, out_cmass_pft_st[stid], 0.0);
					if(!out_cmass_harv_killed_pft_st[stid].invalid())
						outlimit_misc(out, out_cmass_harv_killed_pft_st[stid], 0.0);
					if(!out_transpiration_pft_st[stid].invalid())
						outlimit_misc(out, out_transpiration_pft_st[stid], 0.0);
				}
			}
		}

		outlimit_misc(out, out_forest_cmass_harv_killed_active,  mean_standpft_cmass_harv_killed);

		outlimit_misc(out, out_cmass_active_plus10cm, mean_standpft_cmass_plus10cm);
		outlimit_misc(out, out_forest_cmass_harv_killed_active_plus10cm, mean_standpft_cmass_harv_killed_plus10cm);
		outlimit_misc(out, out_dens_plus10cm, mean_standpft_densindiv_total_plus10cm);

		double height = 0.0;
		double diam = 0.0;
		if (mean_standpft_densindiv_total_plus10cm > 0.0) {
			height = mean_standpft_heightindiv_total_plus10cm / mean_standpft_densindiv_total_plus10cm;
			diam = mean_standpft_diamindiv_total_plus10cm / mean_standpft_densindiv_total_plus10cm;
		}
		outlimit_misc(out, out_speciesheights_plus10cm, height);
		outlimit_misc(out, out_speciesdiam_plus10cm, diam * 100.0);	//diameter output in cm

		// Print to landcover files in case pft:s are common to several landcovers (currently only used in NATURAL and FOREST)
		if (run_landcover) {
			for (int i=0;i<NLANDCOVERTYPES;i++) {
				if (run[i]) {

					switch (i) {
					case CROPLAND:
						if (run[NATURAL]) {
							outlimit_misc(out, out_anpp_cropland,	 mean_standpft_anpp_landscape_lc[i]);
							outlimit_misc(out, out_cmass_cropland, mean_standpft_cmass_landscape_lc[i]);
							outlimit_misc(out, out_lai_cropland,	 mean_standpft_lai_landscape_lc[i]);
						}
						break;
					case PASTURE:
						if (run[NATURAL]) {
							outlimit_misc(out, out_anpp_active_pasture, mean_standpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_active_pasture, mean_standpft_cmass_lc[i]);
							outlimit_misc(out, out_anpp_pasture, mean_standpft_anpp_landscape_lc[i]);
							outlimit_misc(out, out_cmass_pasture, mean_standpft_cmass_landscape_lc[i]);
							outlimit_misc(out, out_lai_pasture, mean_standpft_lai_landscape_lc[i]);
						}
						break;
					case BARREN:
						break;
					case NATURAL:
//						if(run[FOREST] || run[PASTURE]) {
						if(run[FOREST]) {
							outlimit_misc(out, out_anpp_active_natural, mean_standpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_active_natural, mean_standpft_cmass_lc[i]);
							outlimit_misc(out, out_dens_natural, mean_standpft_densindiv_total_lc[i]);
							outlimit_misc(out, out_anpp_natural, mean_standpft_anpp_landscape_lc[i]);
							outlimit_misc(out, out_cmass_natural, mean_standpft_cmass_landscape_lc[i]);
							outlimit_misc(out, out_lai_natural, mean_standpft_lai_landscape_lc[i]);
							outlimit_misc(out, out_aaet_active_natural, mean_standpft_aaet_lc[i]);
							outlimit_misc(out, out_lai_active_natural, mean_standpft_lai_lc[i]);
							outlimit_misc(out, out_fpc_active_natural, mean_standpft_fpc_lc[i]);

							outlimit_misc(out, out_anpp_gridcell_natural,			gridcellpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_gridcell_natural,			gridcellpft_cmass_lc[i]);
							outlimit_misc(out, out_aaet_gridcell_natural,			gridcellpft_aaet_lc[i]);
							outlimit_misc(out, out_lai_gridcell_natural,				gridcellpft_lai_lc[i]);

							double height = 0.0;
							double diam = 0.0;
							if (mean_standpft_densindiv_total_lc[i] > 0.0) {
								height = mean_standpft_heightindiv_total_lc[i] / mean_standpft_densindiv_total_lc[i];
								diam = mean_standpft_diamindiv_total_lc[i] / mean_standpft_densindiv_total_lc[i];
							}
							outlimit_misc(out, out_speciesheights_natural, height);
							outlimit_misc(out, out_speciesdiam_natural, diam * 100.0);	//diameter output in cm
						}
						outlimit_misc(out, out_frac_natural, active_fraction_lc[i]);
						break;
					case FOREST:
						if (run[NATURAL]) {
							outlimit_misc(out, out_anpp_active_forest, mean_standpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_active_forest, mean_standpft_cmass_lc[i]);
							outlimit_misc(out, out_dens_forest, mean_standpft_densindiv_total_lc[i]);
							outlimit_misc(out, out_anpp_forest, mean_standpft_anpp_landscape_lc[i]);
							outlimit_misc(out, out_cmass_forest, mean_standpft_cmass_landscape_lc[i]);
							outlimit_misc(out, out_lai_forest, mean_standpft_lai_landscape_lc[i]);
							outlimit_misc(out, out_aaet_active_forest, mean_standpft_aaet_lc[i]);
							outlimit_misc(out, out_lai_active_forest, mean_standpft_lai_lc[i]);
							outlimit_misc(out, out_fpc_active_forest, mean_standpft_fpc_lc[i]);

							outlimit_misc(out, out_anpp_gridcell_forest,			gridcellpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_gridcell_forest,			gridcellpft_cmass_lc[i]);
							outlimit_misc(out, out_aaet_gridcell_forest,			gridcellpft_aaet_lc[i]);
							outlimit_misc(out, out_lai_gridcell_forest,				gridcellpft_lai_lc[i]);

							double height = 0.0;
							double diam = 0.0;
							if (mean_standpft_densindiv_total_lc[i] > 0.0) {
								height = mean_standpft_heightindiv_total_lc[i] / mean_standpft_densindiv_total_lc[i];
								diam = mean_standpft_diamindiv_total_lc[i] / mean_standpft_densindiv_total_lc[i];
							}
							outlimit_misc(out, out_speciesheights_forest, height);
							outlimit_misc(out, out_speciesdiam_forest, diam * 100.0);	//diameter output in cm
						}
						outlimit_misc(out, out_frac_forest, active_fraction_lc[i]);
						break;
					case URBAN:
						break;
					case PEATLAND:
						if (run[PEATLAND]) {
							outlimit_misc(out, out_anpp_active_peatland, mean_standpft_anpp_lc[i]);
							outlimit_misc(out, out_cmass_active_peatland, mean_standpft_cmass_lc[i]);
							outlimit_misc(out, out_anpp_peatland, mean_standpft_anpp_landscape_lc[i]);
							outlimit_misc(out, out_cmass_peatland, mean_standpft_cmass_landscape_lc[i]);
							outlimit_misc(out, out_lai_peatland, mean_standpft_lai_landscape_lc[i]);
						}
						break;
					default:
						if (date.year == nyear_spinup)
							dprintf("Modify code to deal with landcover output!\n");
					}
				}
			}
		}

		if (pft.landcover == CROPLAND) {
			outlimit_misc(out, out_yield_active,   mean_standpft_yield);
			outlimit_misc(out, out_yield1_active,  mean_standpft_yield1);
			outlimit_misc(out, out_yield2_active,  mean_standpft_yield2);

			int pft_sdate1=-1;
			int pft_sdate2=-1;
			int pft_hdate1=-1;
			int pft_hdate2=-1;
			int pft_lgp=-1;
			double pft_phu=-1;
			double pft_fphu=-1;
			double pft_fhi=-1;

			Gridcell::iterator gc_itr = gridcell.begin();
			while (gc_itr != gridcell.end()) {
				Stand& stand = *gc_itr;

				if (stlist[stand.stid].pftinrotation(pft.name) >= 0) {
					pft_sdate1=stand[0].pft[pft.id].cropphen->sdate_thisyear[0];
					pft_sdate2=stand[0].pft[pft.id].cropphen->sdate_thisyear[1];
					pft_hdate1=stand[0].pft[pft.id].cropphen->hdate_harvest[0];
					pft_hdate2=stand[0].pft[pft.id].cropphen->hdate_harvest[1];
					pft_lgp=stand[0].pft[pft.id].cropphen->lgp;
					pft_phu=stand[0].pft[pft.id].cropphen->phu;
					pft_fphu=stand[0].pft[pft.id].cropphen->fphu_harv;
					pft_fhi=stand[0].pft[pft.id].cropphen->fhi_harv;
				}

				++gc_itr;
			}
			outlimit_misc(out, out_sdate1, pft_sdate1);
			outlimit_misc(out, out_sdate2, pft_sdate2);
			outlimit_misc(out, out_hdate1, pft_hdate1);
			outlimit_misc(out, out_hdate2, pft_hdate2);
			outlimit_misc(out, out_lgp,    pft_lgp);
			outlimit_misc(out, out_phu,    pft_phu);
			outlimit_misc(out, out_fphu,   pft_fphu);
			outlimit_misc(out, out_fhi,    pft_fhi);
		}

		pftlist.nextobj();

	} // *** End of PFT loop ***

	if(printstandtypes) {
		for(int stid=0; stid<nst; stid++) {
			Gridcellst& gcst = gridcell.st[stid];
			if(gcst.frac) {
				if(!out_cmass_pft_st[stid].invalid())
					outlimit_misc(out, out_cmass_pft_st[stid], st_total_cmass[stid] / gcst.frac);
				if(!out_cmass_harv_killed_pft_st[stid].invalid())
					outlimit_misc(out, out_cmass_harv_killed_pft_st[stid], st_total_cmass_harv_killed[stid] / gcst.frac);
				if(!out_transpiration_pft_st[stid].invalid())
					outlimit_misc(out, out_transpiration_pft_st[stid], st_total_aaet[stid] / gcst.frac);
			}
			else {
				if(!out_cmass_pft_st[stid].invalid())
					outlimit_misc(out, out_cmass_pft_st[stid], 0.0);
				if(!out_cmass_harv_killed_pft_st[stid].invalid())
					outlimit_misc(out, out_cmass_harv_killed_pft_st[stid],     0.0);
				if(!out_transpiration_pft_st[stid].invalid())
					outlimit_misc(out, out_transpiration_pft_st[stid], 0.0);
			}
		}
	}

	flux_veg_regrowth = flux_repr_regrowth = flux_soil_regrowth = flux_fire_regrowth = flux_est_regrowth = flux_seed_regrowth 
		= flux_charvest_regrowth = 0.0;
	flux_veg_forestry = flux_repr_forestry = flux_soil_forestry = flux_fire_forestry = flux_est_forestry = flux_seed_forestry 
		= flux_charvest_forestry = 0.0;
	flux_veg_primary = flux_repr_primary = flux_soil_primary = flux_fire_primary = flux_est_primary = flux_seed_primary 
		= flux_charvest_primary = 0.0;
	c_org_leach_gridcell_regrowth = c_org_leach_gridcell_forestry = c_org_leach_gridcell_primary = 0.0;
	surfsoillitterc_regrowth = surfsoillitterc_forestry = surfsoillitterc_primary = 0.0;
	cwdc_regrowth = cwdc_forestry = cwdc_primary = 0.0;
	centuryc_regrowth = centuryc_forestry = centuryc_primary = 0.0;

	double flux_veg_lc[NLANDCOVERTYPES], flux_repr_lc[NLANDCOVERTYPES],
		 flux_soil_lc[NLANDCOVERTYPES], flux_fire_lc[NLANDCOVERTYPES],
		 flux_est_lc[NLANDCOVERTYPES], flux_seed_lc[NLANDCOVERTYPES];
	double flux_charvest_lc[NLANDCOVERTYPES], c_org_leach_lc[NLANDCOVERTYPES];
	double c_litter_lc[NLANDCOVERTYPES], c_fast_lc[NLANDCOVERTYPES],
		  c_slow_lc[NLANDCOVERTYPES], c_harv_slow_lc[NLANDCOVERTYPES];
	double surfsoillitterc_lc[NLANDCOVERTYPES], cwdc_lc[NLANDCOVERTYPES],
		 centuryc_lc[NLANDCOVERTYPES];

	double n_harv_slow_lc[NLANDCOVERTYPES], availn_lc[NLANDCOVERTYPES],
		   andep_lc[NLANDCOVERTYPES], anfert_lc[NLANDCOVERTYPES];
	double anmin_lc[NLANDCOVERTYPES], animm_lc[NLANDCOVERTYPES],
		   anfix_lc[NLANDCOVERTYPES], n_org_leach_lc[NLANDCOVERTYPES],
		   n_min_leach_lc[NLANDCOVERTYPES];
	double flux_ntot_lc[NLANDCOVERTYPES], flux_nharvest_lc[NLANDCOVERTYPES],
		   flux_nseed_lc[NLANDCOVERTYPES];
	double surfsoillittern_lc[NLANDCOVERTYPES], cwdn_lc[NLANDCOVERTYPES],
		   centuryn_lc[NLANDCOVERTYPES];

	//N-transform
	double flux_NH3_soil[NLANDCOVERTYPES],flux_NOx_soil[NLANDCOVERTYPES],flux_N2O_soil[NLANDCOVERTYPES],flux_N2_soil[NLANDCOVERTYPES];
	double flux_N2O_soil_total;

	for (int i=0; i<NLANDCOVERTYPES; i++) {
		flux_veg_lc[i]=0.0;
		flux_repr_lc[i]=0.0;
		flux_soil_lc[i]=0.0;
		flux_fire_lc[i]=0.0;
		flux_est_lc[i]=0.0;
		flux_seed_lc[i]=0.0;
		flux_charvest_lc[i]=0.0;
		c_org_leach_lc[i]=0.0;

		c_litter_lc[i]=0.0;
		c_fast_lc[i]=0.0;
		c_slow_lc[i]=0.0;
		c_harv_slow_lc[i]=0.0;
		surfsoillitterc_lc[i]=0.0;
		cwdc_lc[i]=0.0;
		centuryc_lc[i]=0.0;

		flux_ntot_lc[i]=0.0;
		flux_nharvest_lc[i]=0.0;
		flux_nseed_lc[i]=0.0;

		availn_lc[i]=0.0;
		andep_lc[i]=0.0;	// same value for all land covers
		anfert_lc[i]=0.0;
		anmin_lc[i]=0.0;
		animm_lc[i]=0.0;
		anfix_lc[i]=0.0;
		n_org_leach_lc[i]=0.0;
		n_min_leach_lc[i]=0.0;
		n_harv_slow_lc[i]=0.0;
		surfsoillittern_lc[i]=0.0;
		cwdn_lc[i]=0.0;
		centuryn_lc[i]=0.0;
		//N-transform
		flux_NH3_soil[i]=0.0;
		flux_NOx_soil[i]=0.0;
		flux_N2O_soil[i]=0.0;
		flux_N2_soil[i]=0.0;
	}

	flux_N2O_soil_total = 0.0;
	// Sum C fluxes, dead C pools and runoff across patches

	double forestry_frac = 0.0;
	double regrowth_frac = 0.0;
	double primary_frac = 0.0;

	Gridcell::iterator gc_itr = gridcell.begin();

	// Loop through Stands
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		StandType& st = stlist[stand.stid];
		Gridcellst& gcst = gridcell.st[st.id];
		stand.firstobj();

		if(stand.landcover == NATURAL || stand.landcover == FOREST) {
			if(stand.first_year) {
				if(stand.lc_origin == NATURAL) {
					forestry_frac+=stand.get_gridcell_fraction();
				}
				else {
					regrowth_frac+=stand.get_gridcell_fraction();
				}
			}
			else {
				primary_frac+=stand.get_gridcell_fraction();
			}
		}

		//Loop through Patches
		while (stand.isobj) {
			Patch& patch = stand.getobj();

			double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

			if(stand.landcover == NATURAL || stand.landcover == FOREST) {
				if(stand.first_year) {
					if(stand.lc_origin == NATURAL) {
						flux_veg_forestry+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
						flux_repr_forestry+=-patch.fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
						flux_soil_forestry+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
						flux_fire_forestry+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
						flux_est_forestry+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
						flux_seed_forestry+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
						flux_charvest_forestry+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;
						c_org_leach_gridcell_forestry += patch.soil.aorgCleach * to_gridcell_average;
					}
					else {
						flux_veg_regrowth+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
						flux_repr_regrowth+=-patch.fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
						flux_soil_regrowth+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
						flux_fire_regrowth+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
						flux_est_regrowth+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
						flux_seed_regrowth+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
						flux_charvest_regrowth+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;
						c_org_leach_gridcell_regrowth += patch.soil.aorgCleach * to_gridcell_average;
					}
				}
				else {
					flux_veg_primary+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
					flux_repr_primary+=-patch.fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
					flux_soil_primary+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
					flux_fire_primary+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
					flux_est_primary+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
					flux_seed_primary+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
					flux_charvest_primary+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;
					c_org_leach_gridcell_primary += patch.soil.aorgCleach * to_gridcell_average;
				}
			}

			flux_nseed_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::SEEDN)*to_gridcell_average;
			flux_nharvest_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::HARVESTN)*to_gridcell_average;
			flux_ntot_lc[stand.landcover]+=(patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE) +
					   patch.fluxes.get_annual_flux(Fluxes::NOx_FIRE) +
					   patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE) +
					   patch.fluxes.get_annual_flux(Fluxes::N2_FIRE) +
					   patch.fluxes.get_annual_flux(Fluxes::NH3_SOIL) +
					   patch.fluxes.get_annual_flux(Fluxes::NO_SOIL) +
					   patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL) +
					   patch.fluxes.get_annual_flux(Fluxes::N2_SOIL)) * to_gridcell_average;

			flux_veg_lc[stand.landcover]+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
			flux_repr_lc[stand.landcover]+=-patch.fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
			flux_soil_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
			flux_fire_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
			flux_est_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
			flux_seed_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
			flux_charvest_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;

			double cpool_fast = 0;
			double cpool_slow = 0;
			for (size_t l = 0; l < patch.soil.nsoillayers; l++) {
				const Soillayer &lyr = patch.soil.soillayers[l];
				for (int r = 0; r < NSOMPOOL; r++) {
					cpool_fast += lyr.sompool[r].cmass;
					cpool_slow += lyr.sompool[r].cmass;
				}
			}
			c_fast_lc[stand.landcover]+=cpool_fast*to_gridcell_average;
			c_slow_lc[stand.landcover]+=cpool_slow*to_gridcell_average;

			//Sum slow pools of harvested products
			if (run_landcover && ifslowharvestpool) {
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft=patch.pft[q];
					//slow pool in receiving landcover
					c_harv_slow_lc[stand.landcover]+=patchpft.cmass_harvested_products_slow*to_gridcell_average;
					n_harv_slow_lc[stand.landcover]+=patchpft.nmass_harvested_products_slow*to_gridcell_average;
				}
			}

			//Gridcell irrigation
			irrigation_gridcell += patch.irrigation_y*to_gridcell_average;

			andep_lc[stand.landcover] += (gridcell.aNH4dep + gridcell.aNO3dep) * to_gridcell_average;
			anfert_lc[stand.landcover] += patch.anfert * to_gridcell_average;
			anmin_lc[stand.landcover] += patch.soil.anmin * to_gridcell_average;
			animm_lc[stand.landcover] += patch.soil.animmob * to_gridcell_average;
			anfix_lc[stand.landcover] += patch.soil.anfix * to_gridcell_average;
			n_min_leach_lc[stand.landcover] += patch.soil.aminleach * to_gridcell_average;
			n_org_leach_lc[stand.landcover] += patch.soil.aorgNleach * to_gridcell_average;
			c_org_leach_lc[stand.landcover] += patch.soil.aorgCleach * to_gridcell_average;
			availn_lc[stand.landcover] += (patch.soil.NH4_mass + patch.soil.NO3_mass + patch.soil.snow.NH4_mass 
				+ patch.soil.snow.NO3_mass) * to_gridcell_average;

			//N-transform
			flux_NH3_soil[stand.landcover]	+=patch.fluxes.get_annual_flux(Fluxes::NH3_SOIL)*to_gridcell_average;
			flux_NOx_soil[stand.landcover]	+=patch.fluxes.get_annual_flux(Fluxes::NO_SOIL)*to_gridcell_average;
			flux_N2O_soil[stand.landcover]	+=patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;
			flux_N2_soil[stand.landcover]	+=patch.fluxes.get_annual_flux(Fluxes::N2_SOIL)*to_gridcell_average;;

			flux_N2O_soil_total +=patch.fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;

			for (size_t l = 0; l < patch.soil.nsoillayers; l++) {
				const Soillayer &lyr = patch.soil.soillayers[l];
				for (int r = 0; r < NSOMPOOL; r++) {

					if (r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
						surfsoillitterc_lc[stand.landcover] += lyr.sompool[r].cmass * to_gridcell_average;
						surfsoillittern_lc[stand.landcover] += lyr.sompool[r].nmass * to_gridcell_average;

						if(stand.landcover == NATURAL || stand.landcover == FOREST) {
							if(stand.first_year) {
								if(stand.lc_origin == NATURAL) {
									surfsoillitterc_forestry += lyr.sompool[r].cmass * to_gridcell_average;
								}
								else {
									surfsoillitterc_regrowth += lyr.sompool[r].cmass * to_gridcell_average;
								}
							}
							else {
								surfsoillitterc_primary += lyr.sompool[r].cmass * to_gridcell_average;
							}
						}
					}
					else if (r == SURFFWD || r == SURFCWD) {
						cwdc_lc[stand.landcover] += lyr.sompool[r].cmass * to_gridcell_average;
						cwdn_lc[stand.landcover] += lyr.sompool[r].nmass * to_gridcell_average;

						if(stand.landcover == NATURAL || stand.landcover == FOREST) {
							if(stand.first_year) {
								if(stand.lc_origin == NATURAL) {
									cwdc_forestry += lyr.sompool[r].cmass * to_gridcell_average;
								}
								else {	
									cwdc_regrowth += lyr.sompool[r].cmass * to_gridcell_average;
								}
							}
							else {
								cwdc_primary += lyr.sompool[r].cmass * to_gridcell_average;
							}
						}
					}
					else {
						centuryc_lc[stand.landcover] += lyr.sompool[r].cmass * to_gridcell_average;
						centuryn_lc[stand.landcover]  += lyr.sompool[r].nmass * to_gridcell_average;

						if(stand.landcover == NATURAL || stand.landcover == FOREST) {
							if(stand.first_year) {
								if(stand.lc_origin == NATURAL) {
									centuryc_forestry += lyr.sompool[r].cmass * to_gridcell_average;
								}
								else {
									centuryc_regrowth += lyr.sompool[r].cmass * to_gridcell_average;
								}
							}
							else {
								centuryc_primary += lyr.sompool[r].cmass * to_gridcell_average;
							}
						}
					}
				}
			}

			stand.nextobj();
		} // patch loop
		++gc_itr;
	} // stand loop


	outlimit_misc(out, out_dens_plus10cm, dens_gridcell_plus10cm);

	// Print per-stand totals
	if (printseparatestands) {

		Gridcell::iterator gc_itr = gridcell.begin();
		while (gc_itr != gridcell.end()) {

			Stand& stand = *gc_itr;
			int id = stand.id;;

			if(!out_anpp_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_anpp_stand[id][stand.stid], stand.anpp);
			if(!out_lai_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_lai_stand[id][stand.stid], stand.lai);
			if(!out_cmass_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_cmass_stand[id][stand.stid], stand.cmass);
			if(!out_cmass_wood_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_cmass_wood_stand[id][stand.stid],  stand.cmass_wood);
			if(!out_cmass_wood_harv_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_cmass_wood_harv_stand[id][stand.stid], stand.cmass_wood_harv);
			if(!out_cmass_mort_stand[id][stand.stid].invalid())
				outlimit_misc(out, out_cmass_mort_stand[id][stand.stid], stand.cmass_mort);

			if(!out_storm_stand[id][stand.stid].invalid() && stand.storm_init_done) { // Storm calculations for new stands are done in the beginning of next year

				double stand_sens_ind = 0.0;
				double stand_dam_wood_pot = 0.0;
				double stand_dam_wood = 0.0;
				double stand_storm_cmass_wood_save = 0.0;

				stand.firstobj();

				while (stand.isobj) {
					// START OF LOOP THROUGH PATCHES
					Patch& patch = stand.getobj();

					stand_sens_ind += patch.storm_sens_ind * patch.storm_cmass_wood_save;	// cmass_wood weight
					stand_dam_wood_pot += patch.storm_dam_wood_pot;
					stand_dam_wood += patch.storm_dam_wood;
					stand_storm_cmass_wood_save += patch.storm_cmass_wood_save;
				}

				if(stand_storm_cmass_wood_save) {
					outlimit_misc(out, out_storm_stand[id][stand.stid], stand_sens_ind / stand_storm_cmass_wood_save * 1000.0);		// Sensitivity * 1000; cmass_wood weight
					outlimit_misc(out, out_storm_stand[id][stand.stid], stand_dam_wood_pot / stand_storm_cmass_wood_save * 100.0); // Damaged fraction in %
				}
				else{
					outlimit_misc(out, out_storm_stand[id][stand.stid], 0.0);
					outlimit_misc(out, out_storm_stand[id][stand.stid], 0.0);
				}
				outlimit_misc(out, out_storm_stand[id][stand.stid], stand_dam_wood_pot);
				outlimit_misc(out, out_storm_stand[id][stand.stid], stand_dam_wood);
			}

			++gc_itr;
		}
	}

	// FOREST STRUCTURE OUTPUT
	for(int i=0;i<NFOREST_STRUCTURAL_CLASSES;i++) {

		// age structure
		int age_from, age_to;

		if(i < (NFOREST_STRUCTURAL_CLASSES - 1)) {			// 1-10,...,291-300
			age_from = i*10+1;
			age_to = i*10+10;
		}
		else if(i < NFOREST_STRUCTURAL_CLASSES) {
			age_from = i*10+1;
			age_to = 1500;
		}

		// diameter structure
		double diam_from, diam_to;

		if(i < (NFOREST_STRUCTURAL_CLASSES - 1)) {			//  1-5,...,135-150
			diam_from = (double)(i*5);
			diam_to = (double)(i*5+5);
		}
		else if(i < NFOREST_STRUCTURAL_CLASSES) {					// >150cm
			diam_from = i*5+1;
			diam_to = 1500;
		}
		else {							// no output
			diam_from = 1500;
			diam_to = 10000;
		}

		for(int lc=0; lc<NLANDCOVERTYPES; lc++) {

			double age_dens_lc = 0.0, diam_dens_lc = 0.0, diam_cmass_lc = 0.0, diam_vol_lc = 0.0;
			double diam_cmass_leaf_lc = 0.0, diam_cmass_wood_lc = 0.0, diam_cmass_froot_lc = 0.0, diam_cmass_wood_inc_lc = 0.0;
			double diam_ba_lc = 0.0, standage_lc = 0.0;

			for(int stid=0;stid<nst;stid++) {
				StandType& st = stlist[stid];
				Gridcellst& gcst = gridcell.st[stid];

				double age_dens_st = 0.0, diam_dens_st = 0.0, diam_cmass_st = 0.0, diam_cmass_leaf_st = 0.0, diam_cmass_wood_st = 0.0;
				double diam_cmass_froot_st = 0.0, diam_cmass_wood_inc_st = 0.0,  diam_vol_st = 0.0;
				double standage_st = 0.0;

				if(st.landcover == lc) {

					Gridcell::iterator gc_itr = gridcell.begin();
					while (gc_itr != gridcell.end()) {	

						Stand& stand = *gc_itr;
						int id = stand.id;

						if(stand.stid == stid) {

							double age_dens = 0.0;
							double diam_dens = 0.0;
							double diam_cmass = 0.0;
							double diam_cmass_leaf = 0.0;
							double diam_cmass_wood = 0.0;
							double diam_cmass_froot = 0.0;
							double diam_cmass_wood_inc = 0.0;
							double diam_vol = 0.0;
							double diam_ba = 0.0;
							double standage = 0.0;

							int stand_age = (date.year - stand.first_year);
							if(stand_age >= age_from && stand_age <= age_to)
								standage = 1.0;

							stand.firstobj();

							//Loop through Patches
							while (stand.isobj) {				

								Patch& patch = stand.getobj();

								Vegetation& vegetation = patch.vegetation;

								vegetation.firstobj();
								while (vegetation.isobj) {

									Individual& indiv = vegetation.getobj();
									if (indiv.pft.lifeform == TREE) {

										double trunk_area = PI * pow(indiv.diam / 2.0, 2);

										if(indiv.age >= age_from && indiv.age <= age_to) {
											age_dens += indiv.densindiv / (double)stand.npatch();
										}
										if(indiv.diam >= diam_from / 100.0 && indiv.diam < diam_to / 100.0) {
											diam_dens += indiv.densindiv / (double)stand.npatch();
											diam_cmass += indiv.ccont() / (double)stand.npatch();
											diam_cmass_leaf += indiv.cmass_leaf / (double)stand.npatch();
											diam_cmass_wood += indiv.cmass_wood() / (double)stand.npatch();
											diam_cmass_froot += indiv.cmass_root / (double)stand.npatch();
											diam_cmass_wood_inc += indiv.cmass_wood_inc_5.lastadd() / (double)stand.npatch();
											diam_ba += indiv.densindiv * trunk_area / (double)stand.npatch();
										}
									}

									vegetation.nextobj();
								}
								stand.nextobj();
							}

							diam_cmass_st += diam_cmass * stand.get_gridcell_fraction() / gcst.frac;
							diam_cmass_leaf_st += diam_cmass_leaf * stand.get_gridcell_fraction() / gcst.frac;
							diam_cmass_wood_st += diam_cmass_wood * stand.get_gridcell_fraction() / gcst.frac;
							diam_cmass_froot_st += diam_cmass_froot * stand.get_gridcell_fraction() / gcst.frac;
							diam_cmass_wood_inc_st += diam_cmass_wood_inc * stand.get_gridcell_fraction() / gcst.frac;
							standage_st += standage * stand.get_gridcell_fraction() / gcst.frac;

							age_dens_lc += age_dens * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_dens_lc += diam_dens * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_cmass_lc += diam_cmass * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_cmass_leaf_lc += diam_cmass_leaf * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_cmass_wood_lc += diam_cmass_wood * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_cmass_froot_lc += diam_cmass_froot * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_cmass_wood_inc_lc += diam_cmass_wood_inc * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							diam_ba_lc += diam_ba * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
							standage_lc += standage * stand.get_gridcell_fraction() / gridcell.landcover.frac[lc];
			
							// print to stand output
							if(printseparatestands) {	
								if(!out_agestruct_stand[id][stand.stid].invalid())
									outlimit_misc(out, out_agestruct_stand[id][stand.stid],      age_dens * M2_PER_HA);
								if(!out_diamstruct_stand[id][stand.stid].invalid())
									outlimit_misc(out, out_diamstruct_stand[id][stand.stid],      diam_dens * M2_PER_HA);
								if(!out_diamstruct_cmass_stand[id][stand.stid].invalid())
									outlimit_misc(out, out_diamstruct_cmass_stand[id][stand.stid],      diam_cmass);
							}
						}
						++gc_itr;
					}
					// print to st output
					if(printstandtypes) {
						if(!out_diamstruct_cmass_st[stid].invalid())
							outlimit_misc(out, out_diamstruct_cmass_st[stid], diam_cmass_st);
						if(!out_diamstruct_cmass_leaf_st[stid].invalid())
							outlimit_misc(out, out_diamstruct_cmass_leaf_st[stid], diam_cmass_leaf_st);
						if(!out_diamstruct_cmass_wood_st[stid].invalid())
							outlimit_misc(out, out_diamstruct_cmass_wood_st[stid], diam_cmass_wood_st);
						if(!out_diamstruct_cmass_froot_st[stid].invalid())
							outlimit_misc(out, out_diamstruct_cmass_froot_st[stid], diam_cmass_froot_st);
						if(!out_diamstruct_cmass_wood_inc_st[stid].invalid())
							outlimit_misc(out, out_diamstruct_cmass_wood_inc_st[stid], diam_cmass_wood_inc_st);
						if(!out_standagestruct_st[stid].invalid())
							outlimit_misc(out, out_standagestruct_st[stid], standage_st);
					}
				}
			}


			// print to lc output
			if(lc == NATURAL) {
				outlimit_misc(out, out_agestruct_natural, age_dens_lc * M2_PER_HA);
				outlimit_misc(out, out_diamstruct_natural, diam_dens_lc * M2_PER_HA);
				outlimit_misc(out, out_diamstruct_cmass_natural, diam_cmass_lc);
				outlimit_misc(out, out_diamstruct_cmass_leaf_natural, diam_cmass_leaf_lc);
				outlimit_misc(out, out_diamstruct_cmass_wood_natural, diam_cmass_wood_lc);
				outlimit_misc(out, out_diamstruct_cmass_froot_natural, diam_cmass_froot_lc);
				outlimit_misc(out, out_diamstruct_cmass_wood_inc_natural, diam_cmass_wood_inc_lc);
				outlimit_misc(out, out_diamstruct_ba_natural, diam_ba_lc * M2_PER_HA);
				outlimit_misc(out, out_agestruct_standage_natural, standage_lc);
			}
			else if(lc == FOREST) {
				outlimit_misc(out, out_agestruct_forest, age_dens_lc * M2_PER_HA);
				outlimit_misc(out, out_diamstruct_forest, diam_dens_lc * M2_PER_HA);
				outlimit_misc(out, out_diamstruct_cmass_forest, diam_cmass_lc);
				outlimit_misc(out, out_diamstruct_cmass_leaf_forest, diam_cmass_leaf_lc);
				outlimit_misc(out, out_diamstruct_cmass_wood_forest, diam_cmass_wood_lc);
				outlimit_misc(out, out_diamstruct_cmass_froot_forest, diam_cmass_froot_lc);
				outlimit_misc(out, out_diamstruct_cmass_wood_inc_forest, diam_cmass_wood_inc_lc);
				outlimit_misc(out, out_diamstruct_ba_forest, diam_ba_lc * M2_PER_HA);
				outlimit_misc(out, out_agestruct_standage_forest, standage_lc);
			}
			else if(lc == CROPLAND) {
				outlimit_misc(out, out_agestruct_standage_cropland, standage_lc);
			}
			else if(lc == PASTURE) {
				outlimit_misc(out, out_agestruct_standage_pasture, standage_lc);
			}
		}
	}

	// out_forest_harvest output variables
	double cmass_wood_potharv_forest = 0.0;
	double cmass_wood_potharv_products_forest = 0.0;
	double cmass_potfuel_forest = 0.0;
	double cmass_wood_harv_forest = 0.0;
	double cmass_wood_harv_toprod_forest = 0.0;
	double cmass_harv_tolitter_forest = 0.0;
	double cmass_harv_killed_forest = 0.0;
	double cmass_mort_forest = 0.0;
	double cmass_fire_forest = 0.0;
	double cmass_dist_forest = 0.0;
	double cmass_leaf_root_turnover_forest = 0.0;
	double cmass_repr_forest = 0.0;
	double cmass_est_forest = 0.0;

	double cmass_wood_potharv_natural = 0.0;
	double cmass_wood_potharv_products_natural = 0.0;
	double cmass_potfuel_natural = 0.0;
	double cmass_mort_natural = 0.0;
	double cmass_fire_natural = 0.0;
	double cmass_dist_natural = 0.0;
	double cmass_leaf_root_turnover_natural = 0.0;
	double cmass_repr_natural = 0.0;
	double cmass_est_natural = 0.0;
	double cmass_wood_harv_natural = 0.0;
	double cmass_wood_harv_toprod_natural = 0.0;
	double cmass_harv_tolitter_natural = 0.0;
	double cmass_harv_killed_natural = 0.0;

	// TreeMort stand type output (with lc totals)
	bool treemort_st_output = true;
	bool average_plus10cm_patches_only = true;

	if (treemort_st_output) {

		// Mean values among patches (_patch_mean) are identical with values weighted by stand & stand type area (_st_mean), except for individual-based variables (diam, height)
		// The old code calculated st mean values
		double gridcell_area_frac_plus10cm = 0.0;

		// Tree density-weighted means
		double gridcell_anpp = 0.0;
		double gridcell_agpp = 0.0;
		double gridcell_lai = 0.0;
		double gridcell_cmass = 0.0;	// Test 1
		double gridcell_cmass_potfuel = 0.0;
		double gridcell_cmass_harv_killed = 0.0;
		double gridcell_cmass_dist_clearcut_killed = 0.0;
		double gridcell_cmass_dist_fire_killed = 0.0;
		double gridcell_cmass_dist_wind_killed = 0.0;
		double gridcell_cmass_dist_barkbeetle_killed = 0.0;
		double gridcell_cmass_harv_tolitter = 0.0;
		double gridcell_cmass_mort = 0.0;
		double gridcell_cmass_fire = 0.0;
		double gridcell_cmass_dist = 0.0;
		double gridcell_cmass_leaf_root_turnover = 0.0;
		double gridcell_cmass_repr = 0.0;
		double gridcell_cmass_est = 0.0;

		double gridcell_anpp_tree = 0.0;
		double gridcell_agpp_tree = 0.0;
		double gridcell_aautresp_root_tree = 0.0;
		double gridcell_aautresp_sap_tree = 0.0;
		double gridcell_aautresp_leaf_tree = 0.0;
		double gridcell_aautresp_growth_tree = 0.0;
		double gridcell_lai_tree = 0.0;
		double gridcell_cmass_tree = 0.0;
		double gridcell_cmass_tree_mort = 0.0;

		double gridcell_cmass_wood = 0.0;
		double gridcell_cmass_wood_potharv = 0.0;
		double gridcell_cmass_wood_potharv_products = 0.0;
		double gridcell_cmass_wood_harv_killed = 0.0;
		double gridcell_cmass_wood_harv = 0.0;
		double gridcell_cmass_wood_harv_toprod = 0.0;
		double gridcell_cmass_wood_thin = 0.0;
		double gridcell_cmass_wood_clearcut = 0.0;
		double gridcell_cmass_wood_mort = 0.0;
		double gridcell_cmass_wood_inc = 0.0;

		double gridcell_csoil = 0.0;
		double gridcell_clitter = 0.0;
		double gridcell_clitter_course = 0.0;
		double gridcell_clitter_fine = 0.0;
		double gridcell_rh = 0.0;
		double gridcell_rh_clitter_course = 0.0;
		double gridcell_rh_clitter_fine = 0.0;
		double gridcell_csink = 0.0;
		double gridcell_luc_cflux = 0.0;
		double gridcell_nsoil = 0.0;
		double gridcell_nlitter = 0.0;
		double gridcell_nlitter_course = 0.0;
		double gridcell_nlitter_fine = 0.0;

		double gridcell_aet = 0.0;
		double gridcell_runoff = 0.0;
		double gridcell_intercep = 0.0;
		double gridcell_evap = 0.0;

		double gridcell_storm_sens_ind = 0.0;
		double gridcell_storm_dam_wood_pot = 0.0;
		double gridcell_storm_dam_wood = 0.0;
		double gridcell_storm_cmass_wood_save = 0.0;
		double gridcell_insect_dam_wood = 0.0;
		double gridcell_cmass_wind_killed = 0.0;
		double gridcell_cmass_barkbeetle_killed = 0.0;
		double gridcell_cmass_fire_killed = 0.0;
		double gridcell_cmass_clearcut_killed = 0.0;
		double gridcell_cmass_thin_killed = 0.0;

		double gridcell_dens = 0.0;
		double gridcell_height = 0.0;
		double gridcell_diam = 0.0;
		double gridcell_diam_g = 0.0;

		double gridcell_nstands = 0.0;
		double gridcell_npatches = 0.0;

		double gridcell_cmass_plus10cm = 0.0;
		double gridcell_cmass_harv_killed_plus10cm = 0.0;
		double gridcell_cmass_dist_clearcut_killed_plus10cm = 0.0;
		double gridcell_cmass_dist_fire_killed_plus10cm = 0.0;
		double gridcell_cmass_dist_wind_killed_plus10cm = 0.0;
		double gridcell_cmass_dist_barkbeetle_killed_plus10cm = 0.0;
		double gridcell_cmass_wood_plus10cm = 0.0;
		double gridcell_cmass_wood_harv_killed_plus10cm = 0.0;
		double gridcell_cmass_wood_mort_plus10cm = 0.0;
		double gridcell_cmass_wood_inc_plus10cm = 0.0;

		double gridcell_dens_plus10cm = 0.0;		
		double gridcell_height_plus10cm = 0.0;
		double gridcell_diam_plus10cm = 0.0;
		double gridcell_diam_g_plus10cm = 0.0;

//		double gridcell_dens_patch_mean = 0.0;		
		double gridcell_height_patch_mean = 0.0;
		double gridcell_diam_patch_mean = 0.0;
		double gridcell_diam_g_patch_mean = 0.0;

//		double gridcell_dens_plus10cm_patch_mean = 0.0;		
		double gridcell_height_plus10cm_patch_mean = 0.0;
		double gridcell_diam_plus10cm_patch_mean = 0.0;
		double gridcell_diam_g_plus10cm_patch_mean = 0.0;

		double landcover_area_frac_plus10cm[NLANDCOVERTYPES] = { 0.0 };

		// Tree density-weighted means
		double landcover_anpp[NLANDCOVERTYPES] = { 0.0 };
		double landcover_agpp[NLANDCOVERTYPES] = { 0.0 };
		double landcover_lai[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass[NLANDCOVERTYPES] = { 0.0 };	// Test 1
		double landcover_cmass_potfuel[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_harv_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_clearcut_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_fire_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_wind_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_barkbeetle_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_harv_tolitter[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_mort[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_fire[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_leaf_root_turnover[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_repr[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_est[NLANDCOVERTYPES] = { 0.0 };

		double landcover_anpp_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_agpp_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_aautresp_root_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_aautresp_sap_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_aautresp_leaf_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_aautresp_growth_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_lai_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_tree[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_tree_mort[NLANDCOVERTYPES] = { 0.0 };

		double landcover_cmass_wood[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_potharv[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_potharv_products[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_harv_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_harv[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_harv_toprod[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_thin[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_clearcut[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_mort[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_inc[NLANDCOVERTYPES] = { 0.0 };

		double landcover_csoil[NLANDCOVERTYPES] = { 0.0 };
		double landcover_clitter[NLANDCOVERTYPES] = { 0.0 };
		double landcover_clitter_course[NLANDCOVERTYPES] = { 0.0 };
		double landcover_clitter_fine[NLANDCOVERTYPES] = { 0.0 };
		double landcover_rh[NLANDCOVERTYPES] = { 0.0 };
		double landcover_rh_clitter_course[NLANDCOVERTYPES] = { 0.0 };
		double landcover_rh_clitter_fine[NLANDCOVERTYPES] = { 0.0 };
		double landcover_csink[NLANDCOVERTYPES] = { 0.0 };
		double landcover_luc_cflux[NLANDCOVERTYPES] = { 0.0 };
		double landcover_nsoil[NLANDCOVERTYPES] = { 0.0 };
		double landcover_nlitter[NLANDCOVERTYPES] = { 0.0 };
		double landcover_nlitter_course[NLANDCOVERTYPES] = { 0.0 };
		double landcover_nlitter_fine[NLANDCOVERTYPES] = { 0.0 };

		double landcover_aet[NLANDCOVERTYPES] = { 0.0 };
		double landcover_runoff[NLANDCOVERTYPES] = { 0.0 };
		double landcover_intercep[NLANDCOVERTYPES] = { 0.0 };
		double landcover_evap[NLANDCOVERTYPES] = { 0.0 };

		double landcover_storm_sens_ind[NLANDCOVERTYPES] = { 0.0 };
		double landcover_storm_dam_wood_pot[NLANDCOVERTYPES] = { 0.0 };
		double landcover_storm_dam_wood[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wind_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_storm_cmass_wood_save[NLANDCOVERTYPES] = { 0.0 };
		double landcover_insect_dam_wood[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_barkbeetle_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_fire_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_clearcut_killed[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_thin_killed[NLANDCOVERTYPES] = { 0.0 };

		double landcover_dens[NLANDCOVERTYPES] = { 0.0 };
		double landcover_height[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_g[NLANDCOVERTYPES] = { 0.0 };

		double landcover_nstands[NLANDCOVERTYPES] = { 0.0 };
		double landcover_npatches[NLANDCOVERTYPES] = { 0.0 };

		double landcover_cmass_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_harv_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_clearcut_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_fire_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_wind_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_dist_barkbeetle_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_harv_killed_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_mort_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_wood_inc_plus10cm[NLANDCOVERTYPES] = { 0.0 };

		double landcover_dens_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_height_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_plus10cm[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_g_plus10cm[NLANDCOVERTYPES] = { 0.0 };

		// Unweighted means of patch values
//		double landcover_dens_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_height_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_g_patch_mean[NLANDCOVERTYPES] = { 0.0 };

//		double landcover_dens_plus10cm_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_height_plus10cm_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_plus10cm_patch_mean[NLANDCOVERTYPES] = { 0.0 };
		double landcover_diam_g_plus10cm_patch_mean[NLANDCOVERTYPES] = { 0.0 };

		for (int i = 0;i < nst;i++) {
			StandType& st = stlist[i];
			Gridcellst& gcst = gridcell.st[st.id];

			double st_area_frac_plus10cm = 0.0;

			// Tree density-weighted means
			double st_anpp = 0.0;
			double st_agpp = 0.0;
			double st_lai = 0.0;
			double st_cmass = 0.0;	// Test 1
			double st_cmass_potfuel = 0.0;
			double st_cmass_harv_killed = 0.0;
			double st_cmass_dist_clearcut_killed = 0.0;
			double st_cmass_dist_fire_killed = 0.0;
			double st_cmass_dist_wind_killed = 0.0;
			double st_cmass_dist_barkbeetle_killed = 0.0;
			double st_cmass_harv_tolitter = 0.0;
			double st_cmass_mort = 0.0;
			double st_cmass_fire = 0.0;
			double st_cmass_dist = 0.0;
			double st_cmass_leaf_root_turnover = 0.0;
			double st_cmass_repr = 0.0;
			double st_cmass_est = 0.0;

			double st_anpp_tree = 0.0;
			double st_agpp_tree = 0.0;
			double st_aautresp_root_tree = 0.0;
			double st_aautresp_sap_tree = 0.0;
			double st_aautresp_leaf_tree = 0.0;
			double st_aautresp_growth_tree = 0.0;
			double st_lai_tree = 0.0;
			double st_cmass_tree = 0.0;
			double st_cmass_tree_mort = 0.0;

			double st_cmass_wood = 0.0;
			double st_cmass_wood_potharv = 0.0;
			double st_cmass_wood_potharv_products = 0.0;
			double st_cmass_wood_harv_killed = 0.0;
			double st_cmass_wood_harv = 0.0;
			double st_cmass_wood_harv_toprod = 0.0;
			double st_cmass_wood_thin = 0.0;
			double st_cmass_wood_clearcut = 0.0;
			double st_cmass_wood_mort = 0.0;
			double st_cmass_wood_inc = 0.0;

			double st_csoil = 0.0;
			double st_clitter = 0.0;
			double st_clitter_course = 0.0;
			double st_clitter_fine = 0.0;
			double st_rh = 0.0;
			double st_rh_clitter_course = 0.0;
			double st_rh_clitter_fine = 0.0;
			double st_csink = 0.0;
			double st_luc_cflux = 0.0;
			double st_nsoil = 0.0;
			double st_nlitter = 0.0;
			double st_nlitter_course = 0.0;
			double st_nlitter_fine = 0.0;

			double st_aet = 0.0;
			double st_runoff = 0.0;
			double st_intercep = 0.0;
			double st_evap = 0.0;

			double st_storm_sens_ind = 0.0;
			double st_storm_dam_frac = 0.0;
			double st_storm_dam_wood_pot = 0.0;
			double st_storm_dam_wood = 0.0;
			double st_storm_cmass_wood_save = 0.0;
			double st_insect_dam_wood = 0.0;	
			double st_cmass_wind_killed = 0.0;
			double st_cmass_barkbeetle_killed = 0.0;
			double st_cmass_fire_killed = 0.0;
			double st_cmass_clearcut_killed = 0.0;
			double st_cmass_thin_killed = 0.0;

			double st_dens = 0.0;
			double st_height = 0.0;
			double st_diam = 0.0;
			double st_diam_g = 0.0;

			double st_cmass_plus10cm = 0.0;
			double st_cmass_harv_killed_plus10cm = 0.0;
			double st_cmass_dist_clearcut_killed_plus10cm = 0.0;
			double st_cmass_dist_fire_killed_plus10cm = 0.0;
			double st_cmass_dist_wind_killed_plus10cm = 0.0;
			double st_cmass_dist_barkbeetle_killed_plus10cm = 0.0;
			double st_cmass_wood_plus10cm = 0.0;
			double st_cmass_wood_harv_killed_plus10cm = 0.0;
			double st_cmass_wood_mort_plus10cm = 0.0;
			double st_cmass_wood_inc_plus10cm = 0.0;

			double st_dens_plus10cm = 0.0;
			double st_height_plus10cm = 0.0;
			double st_diam_plus10cm = 0.0;
			double st_diam_g_plus10cm = 0.0;

			// Unweighted means of patch values
//			double st_dens_patch_mean = 0.0;
			double st_height_patch_mean = 0.0;
			double st_diam_patch_mean = 0.0;
			double st_diam_g_patch_mean = 0.0;

//			double st_dens_plus10cm_patch_mean = 0.0;
			double st_height_plus10cm_patch_mean = 0.0;
			double st_diam_plus10cm_patch_mean = 0.0;
			double st_diam_g_plus10cm_patch_mean = 0.0;

			int npatches_cc = 0;
			int npatches_cc_thisyear = 0;
			double cutinterval_mean = 0.0;
			double cutinterval_thisyear_mean = 0.0;
			int nstands_st = 0;
			int npatches_st = 0;

			// Gridcellst-based variables
			if(gcst.frac_old)
				st_luc_cflux += gcst.acflux_harvest / gcst.frac_old;	// acflux_harvest is per gridcell area
			if(gridcell.landcover.frac_old[st.landcover])
				landcover_luc_cflux[st.landcover] = gcst.acflux_harvest / gridcell.landcover.frac_old[st.landcover];
			gridcell_luc_cflux += gcst.acflux_harvest;

			Gridcell::iterator gc_itr = gridcell.begin();

			// Loop through Stands
			while (gc_itr != gridcell.end()) {
				Stand& stand = *gc_itr;

				if (stand.stid == i) {

					double stand_area_frac_plus10cm = 0.0;
					int stand_npatch_plus10cm = 0;

					// Tree density-weighted means
					double stand_anpp = 0.0;
					double stand_agpp = 0.0;
					double stand_lai = 0.0;
					double stand_cmass = 0.0;	// Test 1
					double stand_cmass_potfuel = 0.0;
					double stand_cmass_harv_killed = 0.0;
					double stand_cmass_dist_clearcut_killed = 0.0;
					double stand_cmass_dist_fire_killed = 0.0;
					double stand_cmass_dist_wind_killed = 0.0;
					double stand_cmass_dist_barkbeetle_killed = 0.0;
					double stand_cmass_harv_tolitter = 0.0;
					double stand_cmass_mort = 0.0;
					double stand_cmass_fire = 0.0;
					double stand_cmass_dist = 0.0;
					double stand_cmass_leaf_root_turnover = 0.0;
					double stand_cmass_repr = 0.0;
					double stand_cmass_est = 0.0;

					double stand_anpp_tree = 0.0;
					double stand_agpp_tree = 0.0;
					double stand_aautresp_root_tree = 0.0;
					double stand_aautresp_sap_tree = 0.0;
					double stand_aautresp_leaf_tree = 0.0;
					double stand_aautresp_growth_tree = 0.0;
					double stand_lai_tree = 0.0;
					double stand_cmass_tree = 0.0;
					double stand_cmass_tree_mort = 0.0;

					double stand_cmass_wood = 0.0;
					double stand_cmass_wood_potharv = 0.0;
					double stand_cmass_wood_potharv_products = 0.0;
					double stand_cmass_wood_harv_killed = 0.0;
					double stand_cmass_wood_harv = 0.0;
					double stand_cmass_wood_harv_toprod = 0.0;
					double stand_cmass_wood_thin = 0.0;
					double stand_cmass_wood_clearcut = 0.0;
					double stand_cmass_wood_mort = 0.0;
					double stand_cmass_wood_inc = 0.0;

					double stand_csoil = 0.0;
					double stand_clitter = 0.0;
					double stand_clitter_course = 0.0;
					double stand_clitter_fine = 0.0;
					double stand_rh = 0.0;	
					double stand_rh_clitter_course = 0.0;
					double stand_rh_clitter_fine = 0.0;				
					double stand_csink = 0.0;
					double stand_nsoil = 0.0;
					double stand_nlitter = 0.0;
					double stand_nlitter_course = 0.0;
					double stand_nlitter_fine = 0.0;

					double stand_aet = 0.0;
					double stand_runoff = 0.0;
					double stand_intercep = 0.0;
					double stand_evap = 0.0;

					double stand_storm_sens_ind = 0.0;
					double stand_storm_dam_frac = 0.0;
					double stand_storm_dam_wood_pot = 0.0;
					double stand_storm_dam_wood = 0.0;
					double stand_storm_cmass_wood_save = 0.0;
					double stand_insect_dam_wood = 0.0;
					double stand_cmass_wind_killed = 0.0;
					double stand_cmass_barkbeetle_killed = 0.0;
					double stand_cmass_fire_killed = 0.0;
					double stand_cmass_clearcut_killed = 0.0;
					double stand_cmass_thin_killed = 0.0;

					double stand_dens = 0.0;
					double stand_height = 0.0;
					double stand_diam = 0.0;
					double stand_diam_g = 0.0;

					double stand_cmass_plus10cm = 0.0;
					double stand_cmass_harv_killed_plus10cm = 0.0;
					double stand_cmass_dist_clearcut_killed_plus10cm = 0.0;
					double stand_cmass_dist_fire_killed_plus10cm = 0.0;
					double stand_cmass_dist_wind_killed_plus10cm = 0.0;
					double stand_cmass_dist_barkbeetle_killed_plus10cm = 0.0;

					double stand_cmass_wood_plus10cm = 0.0;
					double stand_cmass_wood_harv_killed_plus10cm = 0.0;
					double stand_cmass_wood_mort_plus10cm = 0.0;
					double stand_cmass_wood_inc_plus10cm = 0.0;

					double stand_dens_plus10cm = 0.0;
					double stand_height_plus10cm = 0.0;
					double stand_diam_plus10cm = 0.0;
					double stand_diam_g_plus10cm = 0.0;

					// Unweighted means of patch values
//					double stand_dens_patch_mean = 0.0;		
					double stand_height_patch_mean = 0.0;
					double stand_diam_patch_mean = 0.0;
					double stand_diam_g_patch_mean = 0.0;

//					double stand_dens_plus10cm_patch_mean = 0.0;		
					double stand_height_plus10cm_patch_mean = 0.0;
					double stand_diam_plus10cm_patch_mean = 0.0;
					double stand_diam_g_plus10cm_patch_mean = 0.0;

					nstands_st++;
					gridcell_nstands++;
					landcover_nstands[stand.landcover]++;

					//Loop through Patches
					stand.firstobj();
					while (stand.isobj) {
						Patch& patch = stand.getobj();
						Vegetation& vegetation = patch.vegetation;

						double patch_dens = 0.0;
						double patch_height = 0.0;
						double patch_diam = 0.0;
						double patch_diam_g = 0.0;

						double patch_dens_plus10cm = 0.0;
						double patch_height_plus10cm = 0.0;
						double patch_diam_plus10cm = 0.0;
						double patch_diam_g_plus10cm = 0.0;

						npatches_st++;
						gridcell_npatches++;
						landcover_npatches[stand.landcover]++;

						// Patch-based variables

						for (size_t l = 0; l < patch.soil.nsoillayers; l++) {
							const Soillayer &lyr = patch.soil.soillayers[l];

							stand_clitter_course += lyr.sompool[SURFCWD].cmass;
							stand_nlitter_course += lyr.sompool[SURFCWD].nmass;
							stand_rh_clitter_course -= lyr.sompool[SURFCWD].delta_cmass * 1000.0;	// Output in g/m2
							stand_clitter_fine += lyr.sompool[SURFFWD].cmass;
							stand_nlitter_fine += lyr.sompool[SURFFWD].nmass;
							stand_rh_clitter_fine -= lyr.sompool[SURFFWD].delta_cmass * 1000.0;	// Output in g/m2

							for (int r = 0; r < NSOMPOOL; r++) {

								if (r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
									stand_clitter += lyr.sompool[r].cmass;
									stand_nlitter += lyr.sompool[r].nmass;
								}
								else if (r == SURFFWD || r == SURFCWD) {
									stand_clitter += lyr.sompool[r].cmass;
									stand_nlitter += lyr.sompool[r].nmass;
								}
								else {
									stand_csoil += lyr.sompool[r].cmass;
									stand_nsoil += lyr.sompool[r].nmass;
								}
							}
						}
						for (int m=0;m<12;m++) {
							stand_rh += patch.fluxes.get_monthly_flux(Fluxes::SOILC, m);
						}

						stand_csink -= (-patch.fluxes.get_annual_flux(Fluxes::NPP) + patch.fluxes.get_annual_flux(Fluxes::REPRC) 
							+ patch.fluxes.get_annual_flux(Fluxes::SOILC) + patch.fluxes.get_annual_flux(Fluxes::FIREC)
							+ patch.fluxes.get_annual_flux(Fluxes::ESTC) + patch.fluxes.get_annual_flux(Fluxes::SEEDC) 
							+ patch.soil.aorgCleach);

						// stand_aet += patch.aaet;
						// stand_runoff += patch.arunoff;
						// stand_intercep += patch.aintercep;
						// stand_evap += patch.aevap;

						stand_storm_sens_ind += patch.storm_sens_ind * patch.storm_cmass_wood_save;	// cmass_wood weight
						stand_storm_dam_wood_pot += patch.storm_dam_wood_pot;
						stand_storm_dam_wood += patch.storm_dam_wood;
						stand_storm_cmass_wood_save += patch.storm_cmass_wood_save;
						stand_insect_dam_wood += patch.insect_dam_wood;

						if(patch.cutinterval_actual) {
							npatches_cc++;
							cutinterval_mean += patch.cutinterval_actual;
						}
						if(patch.cutinterval_actual_thisyear) {
							npatches_cc_thisyear++;
							cutinterval_thisyear_mean += patch.cutinterval_actual_thisyear;
						}

						// Patchpft-based variables
						pftlist.firstobj();
						while (pftlist.isobj) {
							Pft& pft = pftlist.getobj();
							Patchpft& patchpft = patch.pft[pft.id];

							if (pft.lifeform == TREE) {
								stand_cmass_tree_mort += patchpft.cmass_mort;
							}

							stand_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);
							stand_agpp += patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id);

							stand_cmass_harv_killed += patchpft.cmass_harv_killed;
							stand_cmass_dist_clearcut_killed += patchpft.cmass_dist_clearcut_killed;
							stand_cmass_dist_fire_killed += patchpft.cmass_dist_fire_killed;
							stand_cmass_dist_wind_killed += patchpft.cmass_dist_wind_killed;
							stand_cmass_dist_barkbeetle_killed += patchpft.cmass_dist_barkbeetle_killed;
							stand_cmass_harv_tolitter += patchpft.cmass_harv_tolitter;
							stand_cmass_mort += patchpft.cmass_mort;
							stand_cmass_fire += patchpft.cmass_fire_killed + patchpft.cmass_dist_fire_killed;
							stand_cmass_dist += patchpft.cmass_dist;
							stand_cmass_leaf_root_turnover += patchpft.cmass_leaf_root_turnover;
							stand_cmass_repr += patchpft.cmass_repr;
							stand_cmass_est += patchpft.cmass_est;
							stand_clitter += patchpft.cmass_litter_leaf + patchpft.cmass_litter_root + patchpft.cmass_litter_sap
								+ patchpft.cmass_litter_heart + patchpft.cmass_litter_repr;
							// This year's nmass_litter currently belongs to the VegC compartment (as in cpool.out and cmass.out).
//							stand_nlitter += patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap
//								+ patchpft.nmass_litter_heart;

							stand_cmass_wood_harv_killed += patchpft.cmass_wood_harv_killed;
							stand_cmass_wood_harv += patchpft.cmass_wood_harv;
							stand_cmass_wood_harv_toprod += patchpft.cmass_wood_harv_toprod;
							stand_cmass_wood_thin += patchpft.cmass_wood_thin;
							stand_cmass_wood_clearcut += patchpft.cmass_wood_clearcut;
							// stand_cmass_wood_mort += patchpft.cmass_mort_wood;

							stand_cmass_wind_killed += patchpft.cmass_wind_killed;
							stand_cmass_barkbeetle_killed += patchpft.cmass_barkbeetle_killed;
							stand_cmass_fire_killed += patchpft.cmass_fire_killed;
							stand_cmass_clearcut_killed += patchpft.cmass_clearcut_killed;
							stand_cmass_thin_killed += patchpft.cmass_thin_killed;

							stand_cmass_harv_killed_plus10cm += patchpft.cmass_harv_killed_plus10cm;
							stand_cmass_dist_clearcut_killed_plus10cm += patchpft.cmass_dist_clearcut_killed_plus10cm;
							stand_cmass_dist_fire_killed_plus10cm += patchpft.cmass_dist_fire_killed_plus10cm;
							stand_cmass_dist_wind_killed_plus10cm += patchpft.cmass_dist_wind_killed_plus10cm;
							stand_cmass_dist_barkbeetle_killed_plus10cm += patchpft.cmass_dist_barkbeetle_killed_plus10cm;

							// stand_cmass_wood_harv_killed_plus10cm += patchpft.cmass_wood_harv_killed_plus10cm;
							// stand_cmass_wood_mort_plus10cm += patchpft.cmass_mort_wood_plus10cm;

							pftlist.nextobj();
						}

						// Individual-based variables
						vegetation.firstobj();
						while (vegetation.isobj) {
							Individual& indiv = vegetation.getobj();

							if (indiv.id != -1 && indiv.alive) {

								stand_lai += indiv.lai;
								stand_cmass += indiv.ccont();	// Test 1
								stand_cmass_potfuel += (check_harvest_cmass(indiv) - check_harvest_cmass(indiv, true, true));

								stand_cmass_wood += indiv.cmass_wood();
								stand_cmass_wood_potharv += check_harvest_cmass(indiv, true);
								stand_cmass_wood_potharv_products += check_harvest_cmass(indiv, true, true);


//								patch_cmass += indiv.ccont();	// Test 2

								if (indiv.pft.lifeform == TREE) {

									stand_anpp_tree += indiv.anpp;
									stand_agpp_tree += indiv.agpp;
									stand_aautresp_root_tree += indiv.resp_root;
									stand_aautresp_sap_tree += indiv.resp_sap;
									stand_aautresp_leaf_tree += indiv.resp_leaf;
									stand_aautresp_growth_tree += indiv.resp_growth;
									stand_lai_tree += indiv.lai;
									stand_cmass_tree += indiv.ccont();

									stand_cmass_wood_inc += indiv.cmass_wood_inc_5.lastadd();

									stand_dens += indiv.densindiv;
									stand_height += indiv.height * indiv.densindiv;
									stand_diam += indiv.diam * indiv.densindiv;
									stand_diam_g += pow(indiv.diam, 2) * indiv.densindiv;

									patch_dens += indiv.densindiv;
									patch_height += indiv.height * indiv.densindiv;
									patch_diam += indiv.diam * indiv.densindiv;
									patch_diam_g+= pow(indiv.diam, 2) * indiv.densindiv;

									if (indiv.diam >= 0.1) {
										stand_cmass_plus10cm += indiv.ccont();

										stand_cmass_wood_plus10cm += indiv.cmass_wood();
										stand_cmass_wood_inc_plus10cm += indiv.cmass_wood_inc_5.lastadd();

										stand_dens_plus10cm += indiv.densindiv;
										stand_height_plus10cm += indiv.height * indiv.densindiv;
										stand_diam_plus10cm += indiv.diam * indiv.densindiv;
										stand_diam_g_plus10cm += pow(indiv.diam, 2) * indiv.densindiv;

										patch_dens_plus10cm += indiv.densindiv;
										patch_height_plus10cm += indiv.height * indiv.densindiv;
										patch_diam_plus10cm += indiv.diam * indiv.densindiv;
										patch_diam_g_plus10cm += pow(indiv.diam, 2) * indiv.densindiv;
									}
								}
							}
							vegetation.nextobj();
						}

						if(patch_dens) {
							patch_height /= patch_dens;
							patch_diam /= patch_dens;
							patch_diam_g = pow(patch_diam_g / patch_dens, 0.5);	// Complete Dg calculation for each patch
						}

						if (patch_dens_plus10cm) {

							patch_height_plus10cm /= patch_dens_plus10cm;
							patch_diam_plus10cm /= patch_dens_plus10cm;
							patch_diam_g_plus10cm = pow(patch_diam_g_plus10cm / patch_dens_plus10cm, 0.5);	// Complete Dg calculation for each patch
						}

//						stand_cmass_patch_mean += patch_cmass;	// Test 2	EXTRA

						stand_height_patch_mean += patch_height;
						stand_diam_patch_mean += patch_diam;
						stand_diam_g_patch_mean += patch_diam_g;

						stand_height_plus10cm_patch_mean += patch_height_plus10cm;
						stand_diam_plus10cm_patch_mean += patch_diam_plus10cm;
						stand_diam_g_plus10cm_patch_mean += patch_diam_g_plus10cm;

						if (patch_dens_plus10cm)
							stand_npatch_plus10cm++;

						stand.nextobj();
					}

					stand_anpp /= (double)stand.npatch();
					stand_agpp /= (double)stand.npatch();
					stand_lai /= (double)stand.npatch();
					stand_cmass /= (double)stand.npatch();	// Test 1
					stand_cmass_potfuel /= (double)stand.npatch();
					stand_cmass_harv_killed /= (double)stand.npatch();
					stand_cmass_dist_clearcut_killed /= (double)stand.npatch();
					stand_cmass_dist_fire_killed /= (double)stand.npatch();
					stand_cmass_dist_wind_killed /= (double)stand.npatch();
					stand_cmass_dist_barkbeetle_killed /= (double)stand.npatch();
					stand_cmass_harv_tolitter /= (double)stand.npatch();
					stand_cmass_mort /= (double)stand.npatch();
					stand_cmass_fire /= (double)stand.npatch();
					stand_cmass_dist /= (double)stand.npatch();
					stand_cmass_leaf_root_turnover /= (double)stand.npatch();
					stand_cmass_repr /= (double)stand.npatch();
					stand_cmass_est /= (double)stand.npatch();

					stand_anpp_tree /= (double)stand.npatch();
					stand_agpp_tree /= (double)stand.npatch();
					stand_aautresp_root_tree /= (double)stand.npatch();
					stand_aautresp_sap_tree /= (double)stand.npatch();
					stand_aautresp_leaf_tree /= (double)stand.npatch();
					stand_aautresp_growth_tree /= (double)stand.npatch();
					stand_lai_tree /= (double)stand.npatch();
					stand_cmass_tree /= (double)stand.npatch();
					stand_cmass_tree_mort /= (double)stand.npatch();

					stand_cmass_wood /= (double)stand.npatch();
					stand_cmass_wood_potharv /= (double)stand.npatch();
					stand_cmass_wood_potharv_products /= (double)stand.npatch();
					stand_cmass_wood_harv_killed /= (double)stand.npatch();
					stand_cmass_wood_harv /= (double)stand.npatch();
					stand_cmass_wood_harv_toprod /= (double)stand.npatch();
					stand_cmass_wood_thin /= (double)stand.npatch();
					stand_cmass_wood_clearcut /= (double)stand.npatch();
					stand_cmass_wood_mort /= (double)stand.npatch();
					stand_cmass_wood_inc /= (double)stand.npatch();

					stand_csoil /= (double)stand.npatch();
					stand_clitter /= (double)stand.npatch();
					stand_clitter_course /= (double)stand.npatch();
					stand_clitter_fine /= (double)stand.npatch();
					stand_rh /= (double)stand.npatch();
					stand_rh_clitter_course /= (double)stand.npatch();
					stand_rh_clitter_fine /= (double)stand.npatch();				
					stand_csink /= (double)stand.npatch();
					stand_nsoil /= (double)stand.npatch();
					stand_nlitter /= (double)stand.npatch();
					stand_nlitter_course /= (double)stand.npatch();
					stand_nlitter_fine /= (double)stand.npatch();

					stand_aet /= (double)stand.npatch();
					stand_runoff /= (double)stand.npatch();
					stand_intercep /= (double)stand.npatch();
					stand_evap /= (double)stand.npatch();

					stand_storm_sens_ind /= (double)stand.npatch();
					stand_storm_dam_wood_pot /= (double)stand.npatch();
					stand_storm_dam_wood /= (double)stand.npatch();
					stand_storm_cmass_wood_save /= (double)stand.npatch();

					if(stand_storm_cmass_wood_save)
						stand_storm_sens_ind /= stand_storm_cmass_wood_save;	// cmass_wood weight
					else
						stand_storm_sens_ind = 0.0;

					stand_insect_dam_wood /= (double)stand.npatch();
					stand_cmass_wind_killed /= (double)stand.npatch();
					stand_cmass_barkbeetle_killed /= (double)stand.npatch();
					stand_cmass_fire_killed /= (double)stand.npatch();
					stand_cmass_clearcut_killed /= (double)stand.npatch();
					stand_cmass_thin_killed /= (double)stand.npatch();

					stand_dens /= (double)stand.npatch();
					stand_height /= (double)stand.npatch();
					stand_diam /= (double)stand.npatch();
					stand_diam_g /= (double)stand.npatch();

//					stand_cmass_patch_mean /= (double)stand.npatch();	// Test 2

					stand_cmass_plus10cm /= (double)stand.npatch();
					stand_cmass_harv_killed_plus10cm /= (double)stand.npatch();
					stand_cmass_dist_clearcut_killed_plus10cm /= (double)stand.npatch();
					stand_cmass_dist_fire_killed_plus10cm /= (double)stand.npatch();
					stand_cmass_dist_wind_killed_plus10cm /= (double)stand.npatch();
					stand_cmass_dist_barkbeetle_killed_plus10cm /= (double)stand.npatch();

//					stand_dens_plus10cm_patch_mean /= (double)stand.npatch();

					stand_cmass_wood_plus10cm /= (double)stand.npatch();
					stand_cmass_wood_harv_killed_plus10cm /= (double)stand.npatch();
					stand_cmass_wood_mort_plus10cm /= (double)stand.npatch();
					stand_cmass_wood_inc_plus10cm /= (double)stand.npatch();

					stand_dens_plus10cm /= (double)stand.npatch();
					stand_height_plus10cm /= (double)stand.npatch();
					stand_diam_plus10cm /= (double)stand.npatch();
					stand_diam_g_plus10cm /= (double)stand.npatch();

					stand_height_patch_mean /= (double)stand.npatch();
					stand_diam_patch_mean /= (double)stand.npatch();
					stand_diam_g_patch_mean /= (double)stand.npatch();

					if (average_plus10cm_patches_only) {
						if (stand_npatch_plus10cm) {
							stand_height_plus10cm_patch_mean /= (double)stand_npatch_plus10cm;
							stand_diam_plus10cm_patch_mean /= (double)stand_npatch_plus10cm;
							stand_diam_g_plus10cm_patch_mean /= (double)stand_npatch_plus10cm;

							stand_area_frac_plus10cm = stand.get_gridcell_fraction();
						}
					}
					else {
						stand_height_plus10cm_patch_mean /= (double)stand.npatch();
						stand_diam_plus10cm_patch_mean /= (double)stand.npatch();
						stand_diam_g_plus10cm_patch_mean /= (double)stand.npatch();
					}

					st_area_frac_plus10cm += stand_area_frac_plus10cm;

					if (gcst.frac) {
						double stand_st_ratio = stand.get_gridcell_fraction() / gcst.frac;

						// Patch,Patchpft & Individual-based variables

						st_anpp += stand_anpp * stand_st_ratio;
						st_agpp += stand_agpp * stand_st_ratio;
						st_lai += stand_lai * stand_st_ratio;
						st_cmass += stand_cmass * stand_st_ratio;	// Test 1
						st_cmass_potfuel += stand_cmass_potfuel * stand_st_ratio;
						st_cmass_harv_killed += stand_cmass_harv_killed * stand_st_ratio;
						st_cmass_dist_clearcut_killed += stand_cmass_dist_clearcut_killed * stand_st_ratio;
						st_cmass_dist_fire_killed += stand_cmass_dist_fire_killed * stand_st_ratio;
						st_cmass_dist_wind_killed += stand_cmass_dist_wind_killed * stand_st_ratio;
						st_cmass_dist_barkbeetle_killed += stand_cmass_dist_barkbeetle_killed * stand_st_ratio;
						st_cmass_harv_tolitter += stand_cmass_harv_tolitter * stand_st_ratio;
						st_cmass_mort += stand_cmass_mort * stand_st_ratio;
						st_cmass_fire += stand_cmass_fire * stand_st_ratio;
						st_cmass_dist += stand_cmass_dist * stand_st_ratio;
						st_cmass_leaf_root_turnover += stand_cmass_leaf_root_turnover * stand_st_ratio;
						st_cmass_repr += stand_cmass_repr * stand_st_ratio;
						st_cmass_est += stand_cmass_est * stand_st_ratio;

						st_anpp_tree += stand_anpp_tree * stand_st_ratio;
						st_agpp_tree += stand_agpp_tree * stand_st_ratio;
						st_aautresp_root_tree += stand_aautresp_root_tree * stand_st_ratio;
						st_aautresp_sap_tree += stand_aautresp_sap_tree * stand_st_ratio;
						st_aautresp_leaf_tree += stand_aautresp_leaf_tree * stand_st_ratio;
						st_aautresp_growth_tree += stand_aautresp_growth_tree * stand_st_ratio;
						st_lai_tree += stand_lai_tree * stand_st_ratio;
						st_cmass_tree += stand_cmass_tree * stand_st_ratio;
						st_cmass_tree_mort += stand_cmass_tree_mort * stand_st_ratio;

						st_cmass_wood += stand_cmass_wood * stand_st_ratio;
						st_cmass_wood_potharv += stand_cmass_wood_potharv * stand_st_ratio;
						st_cmass_wood_potharv_products += stand_cmass_wood_potharv_products * stand_st_ratio;
						st_cmass_wood_harv_killed += stand_cmass_wood_harv_killed * stand_st_ratio;
						st_cmass_wood_harv += stand_cmass_wood_harv * stand_st_ratio;
						st_cmass_wood_harv_toprod += stand_cmass_wood_harv_toprod * stand_st_ratio;
						st_cmass_wood_thin += stand_cmass_wood_thin * stand_st_ratio;
						st_cmass_wood_clearcut += stand_cmass_wood_clearcut * stand_st_ratio;
						st_cmass_wood_mort += stand_cmass_wood_mort * stand_st_ratio;
						st_cmass_wood_inc += stand_cmass_wood_inc * stand_st_ratio;

						st_csoil += stand_csoil * stand_st_ratio;
						st_clitter += stand_clitter * stand_st_ratio;
						st_clitter_course += stand_clitter_course * stand_st_ratio;
						st_clitter_fine += stand_clitter_fine * stand_st_ratio;
						st_rh += stand_rh * stand_st_ratio;
						st_rh_clitter_course += stand_rh_clitter_course * stand_st_ratio;
						st_rh_clitter_fine += stand_rh_clitter_fine * stand_st_ratio;
						st_csink += stand_csink * stand_st_ratio;
						st_nsoil += stand_nsoil * stand_st_ratio;
						st_nlitter += stand_nlitter * stand_st_ratio;
						st_nlitter_course += stand_nlitter_course * stand_st_ratio;
						st_nlitter_fine += stand_nlitter_fine * stand_st_ratio;

						st_aet += stand_aet * stand_st_ratio;
						st_runoff += stand_runoff * stand_st_ratio;
						st_intercep += stand_intercep * stand_st_ratio;
						st_evap += stand_evap * stand_st_ratio;

						st_storm_sens_ind += stand_storm_sens_ind * stand_storm_cmass_wood_save * stand_st_ratio;	// cmass_wood weight
						st_storm_dam_wood_pot += stand_storm_dam_wood_pot * stand_st_ratio;
						st_storm_dam_wood += stand_storm_dam_wood * stand_st_ratio;
						st_cmass_wind_killed += stand_cmass_wind_killed * stand_st_ratio;
						st_storm_cmass_wood_save += stand_storm_cmass_wood_save * stand_st_ratio;
						st_insect_dam_wood += stand_insect_dam_wood * stand_st_ratio;
						st_cmass_barkbeetle_killed += stand_cmass_barkbeetle_killed * stand_st_ratio;
						st_cmass_fire_killed += stand_cmass_fire_killed * stand_st_ratio;
						st_cmass_clearcut_killed += stand_cmass_clearcut_killed * stand_st_ratio;
						st_cmass_thin_killed += stand_cmass_thin_killed * stand_st_ratio;

						st_dens += stand_dens * stand_st_ratio;
						st_height += stand_height * stand_st_ratio;
						st_diam += stand_diam * stand_st_ratio;
						st_diam_g += stand_diam_g * stand_st_ratio;

						st_cmass_plus10cm += stand_cmass_plus10cm * stand_st_ratio;
						st_cmass_harv_killed_plus10cm += stand_cmass_harv_killed_plus10cm * stand_st_ratio;
						st_cmass_dist_clearcut_killed_plus10cm += stand_cmass_dist_clearcut_killed_plus10cm * stand_st_ratio;
						st_cmass_dist_fire_killed_plus10cm += stand_cmass_dist_fire_killed_plus10cm * stand_st_ratio;
						st_cmass_dist_wind_killed_plus10cm += stand_cmass_dist_wind_killed_plus10cm * stand_st_ratio;
						st_cmass_dist_barkbeetle_killed_plus10cm += stand_cmass_dist_barkbeetle_killed_plus10cm * stand_st_ratio;

						st_cmass_wood_plus10cm += stand_cmass_wood_plus10cm * stand_st_ratio;
						st_cmass_wood_harv_killed_plus10cm += stand_cmass_wood_harv_killed_plus10cm * stand_st_ratio;
						st_cmass_wood_mort_plus10cm += stand_cmass_wood_mort_plus10cm * stand_st_ratio;
						st_cmass_wood_inc_plus10cm += stand_cmass_wood_inc_plus10cm * stand_st_ratio;

						st_dens_plus10cm += stand_dens_plus10cm * stand_st_ratio;
						st_height_plus10cm += stand_height_plus10cm * stand_st_ratio;
						st_diam_plus10cm += stand_diam_plus10cm * stand_st_ratio;
						st_diam_g_plus10cm += stand_diam_g_plus10cm * stand_st_ratio;

						st_height_patch_mean += stand_height_patch_mean * stand_st_ratio;
						st_diam_patch_mean += stand_diam_patch_mean * stand_st_ratio;
						st_diam_g_patch_mean += stand_diam_g_patch_mean * stand_st_ratio;

						if (average_plus10cm_patches_only) {
							st_height_plus10cm_patch_mean += stand_height_plus10cm_patch_mean * stand_area_frac_plus10cm;
							st_diam_plus10cm_patch_mean += stand_diam_plus10cm_patch_mean * stand_area_frac_plus10cm;
							st_diam_g_plus10cm_patch_mean += stand_diam_g_plus10cm_patch_mean * stand_area_frac_plus10cm;
						}
						else {
							st_height_plus10cm_patch_mean += stand_height_plus10cm_patch_mean * stand_st_ratio;
							st_diam_plus10cm_patch_mean += stand_diam_plus10cm_patch_mean * stand_st_ratio;
							st_diam_g_plus10cm_patch_mean += stand_diam_g_plus10cm_patch_mean * stand_st_ratio;
						}
					}

					double stand_lc_ratio = stand.get_gridcell_fraction() / gridcell.landcover.frac[stand.landcover];

					landcover_anpp[stand.landcover] += stand_anpp * stand_lc_ratio;
					landcover_agpp[stand.landcover] += stand_agpp * stand_lc_ratio;
					landcover_lai[stand.landcover] += stand_lai * stand_lc_ratio;
					landcover_cmass[stand.landcover] += stand_cmass * stand_lc_ratio;	// Test 1
					landcover_cmass_potfuel[stand.landcover] += stand_cmass_potfuel * stand_lc_ratio;
					landcover_cmass_harv_killed[stand.landcover] += stand_cmass_harv_killed * stand_lc_ratio;
					landcover_cmass_harv_tolitter[stand.landcover] += stand_cmass_harv_tolitter * stand_lc_ratio;
					landcover_cmass_mort[stand.landcover] += stand_cmass_mort * stand_lc_ratio;
					landcover_cmass_fire[stand.landcover] += stand_cmass_fire * stand_lc_ratio;
					landcover_cmass_dist[stand.landcover] += stand_cmass_dist * stand_lc_ratio;
					landcover_cmass_leaf_root_turnover[stand.landcover] += stand_cmass_leaf_root_turnover * stand_lc_ratio;
					landcover_cmass_repr[stand.landcover] += stand_cmass_repr * stand_lc_ratio;

					landcover_anpp_tree[stand.landcover] += stand_anpp_tree * stand_lc_ratio;
					landcover_agpp_tree[stand.landcover] += stand_agpp_tree * stand_lc_ratio;
					landcover_aautresp_root_tree[stand.landcover] += stand_aautresp_root_tree * stand_lc_ratio;
					landcover_aautresp_sap_tree[stand.landcover] += stand_aautresp_sap_tree * stand_lc_ratio;
					landcover_aautresp_leaf_tree[stand.landcover] += stand_aautresp_leaf_tree * stand_lc_ratio;
					landcover_aautresp_growth_tree[stand.landcover] += stand_aautresp_growth_tree * stand_lc_ratio;
					landcover_lai_tree[stand.landcover] += stand_lai_tree * stand_lc_ratio;
					landcover_cmass_tree[stand.landcover] += stand_cmass_tree * stand_lc_ratio;
					landcover_cmass_tree_mort[stand.landcover] += stand_cmass_tree_mort * stand_lc_ratio;

					landcover_cmass_wood[stand.landcover] += stand_cmass_wood * stand_lc_ratio;
					landcover_cmass_wood_potharv[stand.landcover] += stand_cmass_wood_potharv * stand_lc_ratio;
					landcover_cmass_wood_potharv_products[stand.landcover] += stand_cmass_wood_potharv_products * stand_lc_ratio;
					landcover_cmass_wood_harv_killed[stand.landcover] += stand_cmass_wood_harv_killed * stand_lc_ratio;
					landcover_cmass_dist_clearcut_killed[stand.landcover] += stand_cmass_dist_clearcut_killed * stand_lc_ratio;
					landcover_cmass_dist_fire_killed[stand.landcover] += stand_cmass_dist_fire_killed * stand_lc_ratio;
					landcover_cmass_dist_wind_killed[stand.landcover] += stand_cmass_dist_wind_killed * stand_lc_ratio;
					landcover_cmass_dist_barkbeetle_killed[stand.landcover] += stand_cmass_dist_barkbeetle_killed * stand_lc_ratio;
					landcover_cmass_wood_harv[stand.landcover] += stand_cmass_wood_harv * stand_lc_ratio;
					landcover_cmass_wood_harv_toprod[stand.landcover] += stand_cmass_wood_harv_toprod * stand_lc_ratio;
					landcover_cmass_wood_thin[stand.landcover] += stand_cmass_wood_thin * stand_lc_ratio;
					landcover_cmass_wood_clearcut[stand.landcover] += stand_cmass_wood_clearcut * stand_lc_ratio;
					landcover_cmass_wood_mort[stand.landcover] += stand_cmass_wood_mort * stand_lc_ratio;
					landcover_cmass_wood_inc[stand.landcover] += stand_cmass_wood_inc * stand_lc_ratio;

					landcover_csoil[stand.landcover] += stand_csoil * stand_lc_ratio;
					landcover_clitter[stand.landcover] += stand_clitter * stand_lc_ratio;
					landcover_clitter_course[stand.landcover] += stand_clitter_course * stand_lc_ratio;
					landcover_clitter_fine[stand.landcover] += stand_clitter_fine * stand_lc_ratio;
					landcover_rh[stand.landcover] += stand_rh * stand_lc_ratio;
					landcover_rh_clitter_course[stand.landcover] += stand_rh_clitter_course * stand_lc_ratio;
					landcover_rh_clitter_fine[stand.landcover] += stand_rh_clitter_fine * stand_lc_ratio;					
					landcover_csink[stand.landcover] += stand_csink * stand_lc_ratio;
					landcover_nsoil[stand.landcover] += stand_nsoil * stand_lc_ratio;
					landcover_nlitter[stand.landcover] += stand_nlitter * stand_lc_ratio;
					landcover_nlitter_course[stand.landcover] += stand_nlitter_course * stand_lc_ratio;
					landcover_nlitter_fine[stand.landcover] += stand_nlitter_fine * stand_lc_ratio;

					landcover_aet[stand.landcover] += stand_aet * stand_lc_ratio;
					landcover_runoff[stand.landcover] += stand_runoff * stand_lc_ratio;
					landcover_intercep[stand.landcover] += stand_intercep * stand_lc_ratio;
					landcover_evap[stand.landcover] += stand_evap * stand_lc_ratio;

					landcover_storm_sens_ind[stand.landcover] += stand_storm_sens_ind * stand_storm_cmass_wood_save * stand_lc_ratio;	// cmass_wood weight
					landcover_storm_dam_wood_pot[stand.landcover] += stand_storm_dam_wood_pot * stand_lc_ratio;
					landcover_storm_dam_wood[stand.landcover] += stand_storm_dam_wood * stand_lc_ratio;
					landcover_storm_cmass_wood_save[stand.landcover] += stand_storm_cmass_wood_save * stand_lc_ratio;
					landcover_insect_dam_wood[stand.landcover] += stand_insect_dam_wood * stand_lc_ratio;
					landcover_cmass_wind_killed[stand.landcover] += stand_cmass_wind_killed * stand_lc_ratio;
					landcover_cmass_barkbeetle_killed[stand.landcover] += stand_cmass_barkbeetle_killed * stand_lc_ratio;
					landcover_cmass_fire_killed[stand.landcover] += stand_cmass_fire_killed * stand_lc_ratio;
					landcover_cmass_clearcut_killed[stand.landcover] += stand_cmass_clearcut_killed * stand_lc_ratio;
					landcover_cmass_thin_killed[stand.landcover] += stand_cmass_thin_killed * stand_lc_ratio;

					landcover_dens[stand.landcover] += stand_dens * stand_lc_ratio;
					landcover_height[stand.landcover] += stand_height * stand_lc_ratio;
					landcover_diam[stand.landcover] += stand_diam * stand_lc_ratio;
					landcover_diam_g[stand.landcover] += stand_diam_g * stand_lc_ratio;

					landcover_cmass_plus10cm[stand.landcover] += stand_cmass_plus10cm * stand_lc_ratio;
					landcover_cmass_harv_killed_plus10cm[stand.landcover] += stand_cmass_harv_killed_plus10cm * stand_lc_ratio;
					landcover_cmass_dist_clearcut_killed_plus10cm[stand.landcover] += stand_cmass_dist_clearcut_killed_plus10cm * stand_lc_ratio;
					landcover_cmass_dist_fire_killed_plus10cm[stand.landcover] += stand_cmass_dist_fire_killed_plus10cm * stand_lc_ratio;
					landcover_cmass_dist_wind_killed_plus10cm[stand.landcover] += stand_cmass_dist_wind_killed_plus10cm * stand_lc_ratio;
					landcover_cmass_dist_barkbeetle_killed_plus10cm[stand.landcover] += stand_cmass_dist_barkbeetle_killed_plus10cm * stand_lc_ratio;

					landcover_cmass_wood_plus10cm[stand.landcover] += stand_cmass_wood_plus10cm * stand_lc_ratio;
					landcover_cmass_wood_harv_killed_plus10cm[stand.landcover] += stand_cmass_wood_harv_killed_plus10cm * stand_lc_ratio;
					landcover_cmass_wood_mort_plus10cm[stand.landcover] += stand_cmass_wood_mort_plus10cm * stand_lc_ratio;
					landcover_cmass_wood_inc_plus10cm[stand.landcover] += stand_cmass_wood_inc_plus10cm * stand_lc_ratio;

					landcover_dens_plus10cm[stand.landcover] += stand_dens_plus10cm * stand_lc_ratio;
					landcover_height_plus10cm[stand.landcover] += stand_height_plus10cm * stand_lc_ratio;
					landcover_diam_plus10cm[stand.landcover] += stand_diam_plus10cm * stand_lc_ratio;
					landcover_diam_g_plus10cm[stand.landcover] += stand_diam_g_plus10cm * stand_lc_ratio;

					landcover_height_patch_mean[stand.landcover] += stand_height_patch_mean * stand_lc_ratio;
					landcover_diam_patch_mean[stand.landcover] += stand_diam_patch_mean * stand_lc_ratio;
					landcover_diam_g_patch_mean[stand.landcover] += stand_diam_g_patch_mean * stand_lc_ratio;

					if (average_plus10cm_patches_only) {
						landcover_height_plus10cm_patch_mean[stand.landcover] += stand_height_plus10cm_patch_mean * stand_area_frac_plus10cm;
						landcover_diam_plus10cm_patch_mean[stand.landcover] += stand_diam_plus10cm_patch_mean * stand_area_frac_plus10cm;
						landcover_diam_g_plus10cm_patch_mean[stand.landcover] += stand_diam_g_plus10cm_patch_mean * stand_area_frac_plus10cm;
					}
					else {
						landcover_height_plus10cm_patch_mean[stand.landcover] += stand_height_plus10cm_patch_mean * stand_lc_ratio;
						landcover_diam_plus10cm_patch_mean[stand.landcover] += stand_diam_plus10cm_patch_mean * stand_lc_ratio;
						landcover_diam_g_plus10cm_patch_mean[stand.landcover] += stand_diam_g_plus10cm_patch_mean * stand_lc_ratio;
					}

					if (stand_npatch_plus10cm)
						landcover_area_frac_plus10cm[stand.landcover] += stand_area_frac_plus10cm;

					double stand_gc_ratio = stand.get_gridcell_fraction();

					gridcell_anpp += stand_anpp * stand_gc_ratio;
					gridcell_agpp += stand_agpp * stand_gc_ratio;
					gridcell_lai += stand_lai * stand_gc_ratio;
					gridcell_cmass += stand_cmass * stand_gc_ratio;	// Test 1
					gridcell_cmass_potfuel += stand_cmass_potfuel * stand_gc_ratio;
					gridcell_cmass_harv_killed += stand_cmass_harv_killed * stand_gc_ratio;
					gridcell_cmass_dist_clearcut_killed += stand_cmass_dist_clearcut_killed * stand_gc_ratio;
					gridcell_cmass_dist_fire_killed += stand_cmass_dist_fire_killed * stand_gc_ratio;
					gridcell_cmass_dist_wind_killed += stand_cmass_dist_wind_killed * stand_gc_ratio;
					gridcell_cmass_dist_barkbeetle_killed += stand_cmass_dist_barkbeetle_killed * stand_gc_ratio;
					gridcell_cmass_harv_tolitter += stand_cmass_harv_tolitter * stand_gc_ratio;
					gridcell_cmass_mort += stand_cmass_mort * stand_gc_ratio;
					gridcell_cmass_fire += stand_cmass_fire * stand_gc_ratio;
					gridcell_cmass_dist += stand_cmass_dist * stand_gc_ratio;
					gridcell_cmass_leaf_root_turnover += stand_cmass_leaf_root_turnover * stand_gc_ratio;
					gridcell_cmass_repr += stand_cmass_repr * stand_gc_ratio;

					gridcell_anpp_tree += stand_anpp_tree * stand_gc_ratio;
					gridcell_agpp_tree += stand_agpp_tree * stand_gc_ratio;
					gridcell_aautresp_root_tree += stand_aautresp_root_tree * stand_gc_ratio;
					gridcell_aautresp_sap_tree += stand_aautresp_sap_tree * stand_gc_ratio;
					gridcell_aautresp_leaf_tree += stand_aautresp_leaf_tree * stand_gc_ratio;
					gridcell_aautresp_growth_tree += stand_aautresp_growth_tree * stand_gc_ratio;
					gridcell_lai_tree += stand_lai_tree * stand_gc_ratio;
					gridcell_cmass_tree += stand_cmass_tree * stand_gc_ratio;
					gridcell_cmass_tree_mort += stand_cmass_tree_mort * stand_gc_ratio;

					gridcell_cmass_wood += stand_cmass_wood * stand_gc_ratio;
					gridcell_cmass_wood_potharv += stand_cmass_wood_potharv * stand_gc_ratio;
					gridcell_cmass_wood_potharv_products += stand_cmass_wood_potharv_products * stand_gc_ratio;
					gridcell_cmass_wood_harv_killed += stand_cmass_wood_harv_killed * stand_gc_ratio;
					gridcell_cmass_wood_harv += stand_cmass_wood_harv * stand_gc_ratio;
					gridcell_cmass_wood_harv_toprod += stand_cmass_wood_harv_toprod * stand_gc_ratio;
					gridcell_cmass_wood_thin += stand_cmass_wood_thin * stand_gc_ratio;
					gridcell_cmass_wood_clearcut += stand_cmass_wood_clearcut * stand_gc_ratio;
					gridcell_cmass_wood_mort += stand_cmass_wood_mort * stand_gc_ratio;
					gridcell_cmass_wood_inc += stand_cmass_wood_inc * stand_gc_ratio;

					gridcell_csoil += stand_csoil * stand_gc_ratio;
					gridcell_clitter += stand_clitter * stand_gc_ratio;
					gridcell_clitter_course += stand_clitter_course * stand_gc_ratio;
					gridcell_clitter_fine += stand_clitter_fine * stand_gc_ratio;
					gridcell_rh += stand_rh * stand_gc_ratio;
					gridcell_rh_clitter_course += stand_rh_clitter_course * stand_gc_ratio;
					gridcell_rh_clitter_fine += stand_rh_clitter_fine * stand_gc_ratio;
					gridcell_csink += stand_csink * stand_gc_ratio;
					gridcell_nsoil += stand_nsoil * stand_gc_ratio;
					gridcell_nlitter += stand_nlitter * stand_gc_ratio;
					gridcell_nlitter_course += stand_nlitter_course * stand_gc_ratio;
					gridcell_nlitter_fine += stand_nlitter_fine * stand_gc_ratio;

					gridcell_aet += stand_aet * stand_gc_ratio;
					gridcell_runoff += stand_runoff * stand_gc_ratio;
					gridcell_intercep += stand_intercep * stand_gc_ratio;
					gridcell_evap += stand_evap * stand_gc_ratio;

					gridcell_storm_sens_ind += stand_storm_sens_ind * stand_storm_cmass_wood_save * stand_gc_ratio;	// cmass_wood weight
					gridcell_storm_dam_wood_pot += stand_storm_dam_wood_pot * stand_gc_ratio;
					gridcell_storm_dam_wood += stand_storm_dam_wood * stand_gc_ratio;
					gridcell_storm_cmass_wood_save += stand_storm_cmass_wood_save * stand_gc_ratio;
					gridcell_insect_dam_wood += stand_insect_dam_wood * stand_gc_ratio;
					gridcell_cmass_wind_killed += stand_cmass_wind_killed * stand_gc_ratio;
					gridcell_cmass_barkbeetle_killed += stand_cmass_barkbeetle_killed * stand_gc_ratio;
					gridcell_cmass_fire_killed += stand_cmass_fire_killed * stand_gc_ratio;
					gridcell_cmass_clearcut_killed += stand_cmass_clearcut_killed * stand_gc_ratio;
					gridcell_cmass_thin_killed += stand_cmass_thin_killed * stand_gc_ratio;

					gridcell_dens += stand_dens * stand_gc_ratio;
					gridcell_height += stand_height * stand_gc_ratio;
					gridcell_diam += stand_diam * stand_gc_ratio;
					gridcell_diam_g += stand_diam_g * stand_gc_ratio;

					gridcell_cmass_plus10cm += stand_cmass_plus10cm * stand_gc_ratio;
					gridcell_cmass_harv_killed_plus10cm += stand_cmass_harv_killed_plus10cm * stand_gc_ratio;
					gridcell_cmass_dist_clearcut_killed_plus10cm += stand_cmass_dist_clearcut_killed_plus10cm * stand_gc_ratio;
					gridcell_cmass_dist_fire_killed_plus10cm += stand_cmass_dist_fire_killed_plus10cm * stand_gc_ratio;
					gridcell_cmass_dist_wind_killed_plus10cm += stand_cmass_dist_wind_killed_plus10cm * stand_gc_ratio;
					gridcell_cmass_dist_barkbeetle_killed_plus10cm += stand_cmass_dist_barkbeetle_killed_plus10cm * stand_gc_ratio;

					gridcell_cmass_wood_plus10cm += stand_cmass_wood_plus10cm * stand_gc_ratio;
					gridcell_cmass_wood_harv_killed_plus10cm += stand_cmass_wood_harv_killed_plus10cm * stand_gc_ratio;
					gridcell_cmass_wood_mort_plus10cm += stand_cmass_wood_mort_plus10cm * stand_gc_ratio;
					gridcell_cmass_wood_inc_plus10cm += stand_cmass_wood_inc_plus10cm * stand_gc_ratio;

					gridcell_dens_plus10cm += stand_dens_plus10cm * stand_gc_ratio;
					gridcell_height_plus10cm += stand_height_plus10cm * stand_gc_ratio;
					gridcell_diam_plus10cm += stand_diam_plus10cm * stand_gc_ratio;
					gridcell_diam_g_plus10cm += stand_diam_g_plus10cm * stand_gc_ratio;

					gridcell_height_patch_mean += stand_height_patch_mean * stand_gc_ratio;
					gridcell_diam_patch_mean += stand_diam_patch_mean * stand_gc_ratio;
					gridcell_diam_g_patch_mean += stand_diam_g_patch_mean * stand_gc_ratio;

					if (average_plus10cm_patches_only) {
						gridcell_height_plus10cm_patch_mean += stand_height_plus10cm_patch_mean * stand_area_frac_plus10cm;
						gridcell_diam_plus10cm_patch_mean += stand_diam_plus10cm_patch_mean * stand_area_frac_plus10cm;
						gridcell_diam_g_plus10cm_patch_mean += stand_diam_g_plus10cm_patch_mean * stand_area_frac_plus10cm;
					}
					else {
						gridcell_height_plus10cm_patch_mean += stand_height_plus10cm_patch_mean * stand_gc_ratio;
						gridcell_diam_plus10cm_patch_mean += stand_diam_plus10cm_patch_mean * stand_gc_ratio;
						gridcell_diam_g_plus10cm_patch_mean += stand_diam_g_plus10cm_patch_mean * stand_gc_ratio;
					}

					if (stand_npatch_plus10cm)
						gridcell_area_frac_plus10cm += stand_area_frac_plus10cm;
				}
				++gc_itr;
			}

			if(st_storm_cmass_wood_save)
				st_storm_sens_ind /= st_storm_cmass_wood_save;						// cmass_wood weight
			else
				st_storm_sens_ind = 0.0;

			outlimit_misc(out, out_anpp_sts, st_anpp);
			outlimit_misc(out, out_agpp_sts, st_agpp);
			outlimit_misc(out, out_lai_sts, st_lai);
			outlimit_misc(out, out_cmass_sts, st_cmass);	// TEST 1
//			outlimit_misc(out, out_cmass_potfuel_sts, st_cmass_potfuel);
			outlimit_misc(out, out_cmass_harv_killed_sts, st_cmass_harv_killed);
			outlimit_misc(out, out_cmass_thin_killed_sts, st_cmass_thin_killed);
			outlimit_misc(out, out_cmass_clearcut_killed_sts, st_cmass_clearcut_killed);
			outlimit_misc(out, out_cmass_fire_killed_sts, st_cmass_fire_killed);
			outlimit_misc(out, out_cmass_wind_killed_sts, st_cmass_wind_killed);
			outlimit_misc(out, out_cmass_barkbeetle_killed_sts, st_cmass_barkbeetle_killed);
			outlimit_misc(out, out_cmass_dist_clearcut_killed_sts, st_cmass_dist_clearcut_killed);
			outlimit_misc(out, out_cmass_dist_fire_killed_sts, st_cmass_dist_fire_killed);
			outlimit_misc(out, out_cmass_dist_wind_killed_sts, st_cmass_dist_wind_killed);
			outlimit_misc(out, out_cmass_dist_barkbeetle_killed_sts, st_cmass_dist_barkbeetle_killed);
			outlimit_misc(out, out_cmass_dist_prescribed_killed_sts, st_cmass_dist_clearcut_killed + 
				st_cmass_dist_fire_killed + st_cmass_dist_wind_killed + st_cmass_dist_barkbeetle_killed);
//			outlimit_misc(out, out_cmass_harv_tolitter_sts, st_cmass_harv_tolitter);
//			outlimit_misc(out, out_cmass_mort_sts, st_cmass_mort);
//			outlimit_misc(out, out_cmass_fire_sts, st_cmass_fire);
			outlimit_misc(out, out_cmass_dist_default_sts, st_cmass_dist);
//			outlimit_misc(out, out_cmass_leaf_root_turnover_sts, st_cmass_leaf_root_turnover);
//			outlimit_misc(out, out_cmass_repr_sts, st_cmass_repr);

			outlimit_misc(out, out_anpp_tree_sts, st_anpp_tree);
			outlimit_misc(out, out_agpp_tree_sts, st_agpp_tree);
			outlimit_misc(out, out_aautresp_tree_sts, (st_aautresp_root_tree + st_aautresp_sap_tree + st_aautresp_leaf_tree + st_aautresp_growth_tree));
			outlimit_misc(out, out_aautresp_root_tree_sts, st_aautresp_root_tree);
			outlimit_misc(out, out_aautresp_sap_tree_sts, st_aautresp_sap_tree);
			outlimit_misc(out, out_aautresp_leaf_tree_sts, st_aautresp_leaf_tree);
			outlimit_misc(out, out_aautresp_growth_tree_sts, st_aautresp_growth_tree);
			outlimit_misc(out, out_lai_tree_sts, st_lai_tree);
			outlimit_misc(out, out_cmass_tree_sts, st_cmass_tree);
			outlimit_misc(out, out_cmass_tree_mort_sts, st_cmass_tree_mort);

			outlimit_misc(out, out_cmass_wood_sts, st_cmass_wood);
//			outlimit_misc(out, out_cmass_wood_potharv_sts, st_cmass_wood_potharv);
//			outlimit_misc(out, out_cmass_wood_potharv_products_sts, st_cmass_wood_potharv_products);
			outlimit_misc(out, out_cmass_wood_harv_killed_sts, st_cmass_wood_harv_killed);
			outlimit_misc(out, out_cmass_wood_harv_sts, st_cmass_wood_harv);
			outlimit_misc(out, out_cmass_wood_harv_toprod_sts, st_cmass_wood_harv_toprod);
			outlimit_misc(out, out_cmass_wood_thin_sts, st_cmass_wood_thin);
			outlimit_misc(out, out_cmass_wood_clearcut_sts, st_cmass_wood_clearcut);
			outlimit_misc(out, out_cmass_wood_mort_sts, st_cmass_wood_mort);
			outlimit_misc(out, out_cmass_wood_inc_sts, st_cmass_wood_inc);

			outlimit_misc(out, out_csoil_sts, st_csoil);
			outlimit_misc(out, out_clitter_sts, st_clitter);
			outlimit_misc(out, out_clitter_course_sts, st_clitter_course);
			outlimit_misc(out, out_clitter_fine_sts, st_clitter_fine);
			outlimit_misc(out, out_rh_sts, st_rh);
			outlimit_misc(out, out_rh_clitter_course_sts, st_rh_clitter_course);
			outlimit_misc(out, out_rh_clitter_fine_sts, st_rh_clitter_fine);
			outlimit_misc(out, out_csink_sts, st_csink);
			outlimit_misc(out, out_luc_cflux_sts, st_luc_cflux);
			outlimit_misc(out, out_nsoil_sts, st_nsoil);
			outlimit_misc(out, out_nlitter_sts, st_nlitter);
			outlimit_misc(out, out_nlitter_course_sts, st_nlitter_course);
			outlimit_misc(out, out_nlitter_fine_sts, st_nlitter_fine);

			outlimit_misc(out, out_transpiration_sts, st_aet);
			outlimit_misc(out, out_runoff_sts, st_runoff);
			outlimit_misc(out, out_intercep_sts, st_intercep);
			outlimit_misc(out, out_evap_sts, st_evap);

			outlimit_misc(out, out_storm_si_sts, st_storm_sens_ind * 1000.0);	// Storm sensitivity * 1000
			if(st_storm_cmass_wood_save)
				outlimit_misc(out, out_storm_frac_pot_sts, st_storm_dam_wood_pot / st_storm_cmass_wood_save * 100.0);	// Potential damaged fraction in %
			else
				outlimit_misc(out, out_storm_frac_pot_sts, 0.0);
			if(st_storm_cmass_wood_save)
				outlimit_misc(out, out_storm_frac_sts, st_storm_dam_wood / st_storm_cmass_wood_save * 100.0);	// Actual damaged fraction in %
			else
				outlimit_misc(out, out_storm_frac_sts, 0.0);
			outlimit_misc(out, out_storm_dw_pot_sts, st_storm_dam_wood_pot);
			outlimit_misc(out, out_storm_dw_sts, st_storm_dam_wood);
			outlimit_misc(out, out_insect_dw_sts, st_insect_dam_wood);
			outlimit_misc(out, out_cmort_storm_sts, st_cmass_wind_killed + st_cmass_dist_wind_killed);
			outlimit_misc(out, out_cmort_insect_sts, st_cmass_barkbeetle_killed + st_cmass_dist_barkbeetle_killed);
			outlimit_misc(out, out_cmort_fire_sts, st_cmass_fire_killed + st_cmass_dist_fire_killed);
			outlimit_misc(out, out_cmort_clearcut_sts, st_cmass_clearcut_killed + st_cmass_dist_clearcut_killed);

			outlimit_misc(out, out_dens_sts, st_dens);

			outlimit_misc(out, out_frac_sts, gcst.frac);

			double height = 0.0;
			double diam = 0.0;
			double diam_g = 0.0;

			if (st_dens) {
				height = st_height / st_dens;
				diam = st_diam / st_dens;
				diam_g = st_diam_g / st_dens;
				diam_g = pow(diam_g, 0.5);
			}

			outlimit_misc(out, out_height_sts, height);
			outlimit_misc(out, out_diam_sts, diam * 100.0);		//diameter output in cm
			outlimit_misc(out, out_diam_g_sts, diam_g * 100.0);	//diameter output in cm;

			outlimit_misc(out, out_cmass_plus10cm_sts, st_cmass_plus10cm);
			outlimit_misc(out, out_cmass_harv_killed_plus10cm_sts, st_cmass_harv_killed_plus10cm);
			outlimit_misc(out, out_cmass_dist_clearcut_killed_plus10cm_sts, st_cmass_dist_clearcut_killed_plus10cm);
			outlimit_misc(out, out_cmass_dist_fire_killed_plus10cm_sts, st_cmass_dist_fire_killed_plus10cm);
			outlimit_misc(out, out_cmass_dist_wind_killed_plus10cm_sts, st_cmass_dist_wind_killed_plus10cm);
			outlimit_misc(out, out_cmass_dist_barkbeetle_killed_plus10cm_sts, st_cmass_dist_barkbeetle_killed_plus10cm);

			outlimit_misc(out, out_cmass_wood_plus10cm_sts, st_cmass_wood_plus10cm);
			outlimit_misc(out, out_cmass_wood_harv_killed_plus10cm_sts, st_cmass_wood_harv_killed_plus10cm);
			outlimit_misc(out, out_cmass_wood_mort_plus10cm_sts, st_cmass_wood_mort_plus10cm);
			outlimit_misc(out, out_cmass_wood_inc_plus10cm_sts, st_cmass_wood_inc_plus10cm);

			outlimit_misc(out, out_dens_plus10cm_sts, st_dens_plus10cm);

			height = 0.0;
			diam = 0.0;
			diam_g = 0.0;

			if (st_dens_plus10cm) {
				height = st_height_plus10cm / st_dens_plus10cm;
				diam = st_diam_plus10cm / st_dens_plus10cm;
				diam_g = st_diam_g_plus10cm / st_dens_plus10cm;
				diam_g = pow(diam_g, 0.5);
			}

			outlimit_misc(out, out_height_plus10cm_sts, height);
			outlimit_misc(out, out_diam_plus10cm_sts, diam * 100.0);		//diameter output in cm
			outlimit_misc(out, out_diam_g_plus10cm_sts, diam_g * 100.0);	//diameter output in cm;

			outlimit_misc(out, out_height_sts_patch_mean, st_height_patch_mean);
			outlimit_misc(out, out_diam_sts_patch_mean, st_diam_patch_mean * 100.0);	//diameter output in cm
			outlimit_misc(out, out_diam_g_sts_patch_mean, st_diam_g_patch_mean * 100.0);	//diameter output in cm

			if (average_plus10cm_patches_only) {
				if (st_area_frac_plus10cm) {
					st_height_plus10cm_patch_mean /= st_area_frac_plus10cm;
					st_diam_plus10cm_patch_mean /= st_area_frac_plus10cm;
					st_diam_g_plus10cm_patch_mean /= st_area_frac_plus10cm;
				}
			}
			
			outlimit_misc(out, out_height_plus10cm_sts_patch_mean, st_height_plus10cm_patch_mean);
			outlimit_misc(out, out_diam_plus10cm_sts_patch_mean, st_diam_plus10cm_patch_mean * 100.0);	//diameter output in cm
			outlimit_misc(out, out_diam_g_plus10cm_sts_patch_mean, st_diam_g_plus10cm_patch_mean * 100.0);	//diameter output in cm

			outlimit_misc(out, out_cutinterval_sts, cutinterval_mean / max(1, npatches_cc));
			outlimit_misc(out, out_cutinterval_thisyear_sts, cutinterval_thisyear_mean / max(1, npatches_cc_thisyear));
			outlimit_misc(out, out_nstand_sts, nstands_st);
			outlimit_misc(out, out_npatch_sts, npatches_st);

			if(st.landcover == FOREST) {

				cmass_wood_potharv_forest += st_cmass_wood_potharv * gcst.frac;
				cmass_wood_harv_forest += st_cmass_wood_harv * gcst.frac;											// Standard harvest plus prescribed clearcut
				cmass_wood_harv_toprod_forest += st_cmass_wood_harv_toprod * gcst.frac;								// Standard harvest plus prescribed clearcut
				cmass_harv_tolitter_forest += st_cmass_harv_tolitter * gcst.frac;									// Standard harvest plus prescribed clearcut
				cmass_harv_killed_forest += (st_cmass_harv_killed + st_cmass_dist_clearcut_killed + st_cmass_dist_wind_killed + st_cmass_dist_barkbeetle_killed) * gcst.frac;		// Standard harvest plus prescribed clearcut
				cmass_mort_forest += st_cmass_mort * gcst.frac;
				cmass_fire_forest += st_cmass_fire * gcst.frac;
				cmass_dist_forest += st_cmass_dist * gcst.frac;
				cmass_leaf_root_turnover_forest += st_cmass_leaf_root_turnover * gcst.frac;
				cmass_repr_forest += st_cmass_repr * gcst.frac;
				cmass_est_forest += st_cmass_est * gcst.frac;
				cmass_wood_potharv_products_forest += st_cmass_wood_potharv_products * gcst.frac;
				cmass_potfuel_forest += st_cmass_potfuel * gcst.frac;
			}
			else if(st.landcover == NATURAL) {

				cmass_wood_potharv_natural += st_cmass_wood_potharv * gcst.frac;
				cmass_wood_potharv_products_natural += st_cmass_wood_potharv_products * gcst.frac;
				cmass_potfuel_natural += st_cmass_potfuel * gcst.frac;
				cmass_mort_natural += st_cmass_mort * gcst.frac;
				cmass_fire_natural += st_cmass_fire * gcst.frac;
				cmass_dist_natural += st_cmass_dist * gcst.frac;
				cmass_leaf_root_turnover_natural += st_cmass_leaf_root_turnover * gcst.frac;
				cmass_repr_natural += st_cmass_repr * gcst.frac;
				cmass_est_natural += st_cmass_est * gcst.frac;
				cmass_wood_harv_natural += st_cmass_wood_harv * gcst.frac;
				cmass_wood_harv_toprod_natural += st_cmass_wood_harv_toprod * gcst.frac;
				cmass_harv_tolitter_natural += st_cmass_harv_tolitter * gcst.frac;
				cmass_harv_killed_natural += (st_cmass_harv_killed + st_cmass_dist_clearcut_killed + st_cmass_dist_wind_killed + st_cmass_dist_barkbeetle_killed) * gcst.frac;
			}
		}

		if(gridcell_storm_cmass_wood_save)
			gridcell_storm_sens_ind /= gridcell_storm_cmass_wood_save;			// cmass_wood weight
		else
			gridcell_storm_sens_ind = 0.0;

		outlimit_misc(out, out_anpp_sts, gridcell_anpp);
		outlimit_misc(out, out_agpp_sts, gridcell_agpp);
		outlimit_misc(out, out_lai_sts, gridcell_lai);
		outlimit_misc(out, out_cmass_sts, gridcell_cmass);	// TEST 1
//		outlimit_misc(out, out_cmass_potfuel_sts, gridcell_cmass_potfuel);				// Decide later which output is needed here
		outlimit_misc(out, out_cmass_harv_killed_sts, gridcell_cmass_harv_killed);
		outlimit_misc(out, out_cmass_thin_killed_sts, gridcell_cmass_thin_killed);
		outlimit_misc(out, out_cmass_clearcut_killed_sts, gridcell_cmass_clearcut_killed);
		outlimit_misc(out, out_cmass_fire_killed_sts, gridcell_cmass_fire_killed);
		outlimit_misc(out, out_cmass_wind_killed_sts, gridcell_cmass_wind_killed);
		outlimit_misc(out, out_cmass_barkbeetle_killed_sts, gridcell_cmass_barkbeetle_killed);
		outlimit_misc(out, out_cmass_dist_clearcut_killed_sts, gridcell_cmass_dist_clearcut_killed);
		outlimit_misc(out, out_cmass_dist_fire_killed_sts, gridcell_cmass_dist_fire_killed);
		outlimit_misc(out, out_cmass_dist_wind_killed_sts, gridcell_cmass_dist_wind_killed);
		outlimit_misc(out, out_cmass_dist_barkbeetle_killed_sts, gridcell_cmass_dist_barkbeetle_killed);
		outlimit_misc(out, out_cmass_dist_prescribed_killed_sts, gridcell_cmass_dist_clearcut_killed + 
			gridcell_cmass_dist_fire_killed + gridcell_cmass_dist_wind_killed + gridcell_cmass_dist_barkbeetle_killed);
//		outlimit_misc(out, out_cmass_harv_tolitter_sts, gridcell_cmass_harv_tolitter);
//		outlimit_misc(out, out_cmass_mort_sts, gridcell_cmass_mort);
//		outlimit_misc(out, out_cmass_fire_sts, gridcell_cmass_fire);
		outlimit_misc(out, out_cmass_dist_default_sts, gridcell_cmass_dist);
//		outlimit_misc(out, out_cmass_leaf_root_turnover_sts, gridcell_cmass_leaf_root_turnover);
//		outlimit_misc(out, out_cmass_repr_sts, gridcell_cmass_repr);

		outlimit_misc(out, out_anpp_tree_sts, gridcell_anpp_tree);
		outlimit_misc(out, out_agpp_tree_sts, gridcell_agpp_tree);
		outlimit_misc(out, out_aautresp_tree_sts, (gridcell_aautresp_root_tree + gridcell_aautresp_sap_tree + gridcell_aautresp_leaf_tree + gridcell_aautresp_growth_tree));
		outlimit_misc(out, out_aautresp_root_tree_sts, gridcell_aautresp_root_tree);
		outlimit_misc(out, out_aautresp_sap_tree_sts, gridcell_aautresp_sap_tree);
		outlimit_misc(out, out_aautresp_leaf_tree_sts, gridcell_aautresp_leaf_tree);
		outlimit_misc(out, out_aautresp_growth_tree_sts, gridcell_aautresp_growth_tree);
		outlimit_misc(out, out_lai_tree_sts, gridcell_lai_tree);
		outlimit_misc(out, out_cmass_tree_sts, gridcell_cmass_tree);
		outlimit_misc(out, out_cmass_tree_mort_sts, gridcell_cmass_tree_mort);

		outlimit_misc(out, out_cmass_wood_sts, gridcell_cmass_wood);
//		outlimit_misc(out, out_cmass_wood_potharv_sts, gridcell_cmass_wood_potharv);			// Decide later which output is needed here
//		outlimit_misc(out, out_cmass_wood_potharv_products_sts, gridcell_cmass_wood_potharv_products);
		outlimit_misc(out, out_cmass_wood_harv_killed_sts, gridcell_cmass_wood_harv_killed);
		outlimit_misc(out, out_cmass_wood_harv_sts, gridcell_cmass_wood_harv);
		outlimit_misc(out, out_cmass_wood_harv_toprod_sts, gridcell_cmass_wood_harv_toprod);
		outlimit_misc(out, out_cmass_wood_thin_sts, gridcell_cmass_wood_thin);
		outlimit_misc(out, out_cmass_wood_clearcut_sts, gridcell_cmass_wood_clearcut);
		outlimit_misc(out, out_cmass_wood_mort_sts, gridcell_cmass_wood_mort);
		outlimit_misc(out, out_cmass_wood_inc_sts, gridcell_cmass_wood_inc);

		outlimit_misc(out, out_csoil_sts, gridcell_csoil);
		outlimit_misc(out, out_clitter_sts, gridcell_clitter);
		outlimit_misc(out, out_clitter_course_sts, gridcell_clitter_course);
		outlimit_misc(out, out_clitter_fine_sts, gridcell_clitter_fine);
		outlimit_misc(out, out_rh_sts, gridcell_rh);
		outlimit_misc(out, out_rh_clitter_course_sts, gridcell_rh_clitter_course);
		outlimit_misc(out, out_rh_clitter_fine_sts, gridcell_rh_clitter_fine);
		outlimit_misc(out, out_csink_sts, gridcell_csink);
		outlimit_misc(out, out_luc_cflux_sts, gridcell_luc_cflux);
		outlimit_misc(out, out_nsoil_sts, gridcell_nsoil);
		outlimit_misc(out, out_nlitter_sts, gridcell_nlitter);
		outlimit_misc(out, out_nlitter_course_sts, gridcell_nlitter_course);
		outlimit_misc(out, out_nlitter_fine_sts, gridcell_nlitter_fine);

		outlimit_misc(out, out_runoff_sts, gridcell_runoff);
		outlimit_misc(out, out_transpiration_sts, gridcell_aet);
		outlimit_misc(out, out_intercep_sts, gridcell_intercep);
		outlimit_misc(out, out_evap_sts, gridcell_evap);

		outlimit_misc(out, out_storm_si_sts, gridcell_storm_sens_ind * 1000.0);	// Storm sensitivity * 1000
		if(gridcell_storm_cmass_wood_save)
			outlimit_misc(out, out_storm_frac_pot_sts, gridcell_storm_dam_wood_pot / gridcell_storm_cmass_wood_save * 100.0);	// Potential damaged fraction in %
		else
			outlimit_misc(out, out_storm_frac_pot_sts, 0.0);
		if(gridcell_storm_cmass_wood_save)
			outlimit_misc(out, out_storm_frac_sts, gridcell_storm_dam_wood / gridcell_storm_cmass_wood_save * 100.0);	// Actual damaged fraction in %
		else
			outlimit_misc(out, out_storm_frac_sts, 0.0);
		outlimit_misc(out, out_storm_dw_pot_sts, gridcell_storm_dam_wood_pot);
		outlimit_misc(out, out_storm_dw_sts, gridcell_storm_dam_wood);
		outlimit_misc(out, out_cmort_storm_sts, gridcell_cmass_wind_killed + gridcell_cmass_dist_wind_killed);
		outlimit_misc(out, out_insect_dw_sts, gridcell_insect_dam_wood);
		outlimit_misc(out, out_cmort_insect_sts, gridcell_cmass_barkbeetle_killed + gridcell_cmass_dist_barkbeetle_killed);
		outlimit_misc(out, out_cmort_fire_sts, gridcell_cmass_fire_killed + gridcell_cmass_dist_fire_killed);
		outlimit_misc(out, out_cmort_clearcut_sts, gridcell_cmass_clearcut_killed + gridcell_cmass_dist_clearcut_killed);

		outlimit_misc(out, out_dens_sts, gridcell_dens);

		outlimit_misc(out, out_frac_sts, 1.0);

		outlimit_misc(out, out_nstand_sts, gridcell_nstands);
		outlimit_misc(out, out_npatch_sts, gridcell_npatches);

		double height = 0.0;
		double diam = 0.0;
		double diam_g = 0.0;
		if (gridcell_dens) {
			height = gridcell_height / gridcell_dens;
			diam = gridcell_diam / gridcell_dens;
			diam_g = gridcell_diam_g / gridcell_dens;
			diam_g = pow(diam_g, 0.5);
		}

		outlimit_misc(out, out_height_sts, height);
		outlimit_misc(out, out_diam_sts, diam * 100.0);		//diameter output in cm
		outlimit_misc(out, out_diam_g_sts, diam_g * 100.0);	//diameter output in cm

		outlimit_misc(out, out_cmass_plus10cm_sts, gridcell_cmass_plus10cm);
		outlimit_misc(out, out_cmass_harv_killed_plus10cm_sts, gridcell_cmass_harv_killed_plus10cm);
		outlimit_misc(out, out_cmass_dist_clearcut_killed_plus10cm_sts, gridcell_cmass_dist_clearcut_killed_plus10cm);
		outlimit_misc(out, out_cmass_dist_fire_killed_plus10cm_sts, gridcell_cmass_dist_fire_killed_plus10cm);
		outlimit_misc(out, out_cmass_dist_wind_killed_plus10cm_sts, gridcell_cmass_dist_wind_killed_plus10cm);
		outlimit_misc(out, out_cmass_dist_barkbeetle_killed_plus10cm_sts, gridcell_cmass_dist_barkbeetle_killed_plus10cm);

		outlimit_misc(out, out_cmass_wood_plus10cm_sts, gridcell_cmass_wood_plus10cm);
		outlimit_misc(out, out_cmass_wood_harv_killed_plus10cm_sts, gridcell_cmass_wood_harv_killed_plus10cm);
		outlimit_misc(out, out_cmass_wood_mort_plus10cm_sts, gridcell_cmass_wood_mort_plus10cm);
		outlimit_misc(out, out_cmass_wood_inc_plus10cm_sts, gridcell_cmass_wood_inc_plus10cm);

		outlimit_misc(out, out_dens_plus10cm_sts, gridcell_dens_plus10cm);

		height = 0.0;
		diam = 0.0;
		diam_g = 0.0;
		if (gridcell_dens_plus10cm) {
			height = gridcell_height_plus10cm / gridcell_dens_plus10cm;
			diam = gridcell_diam_plus10cm / gridcell_dens_plus10cm;
			diam_g = gridcell_diam_g_plus10cm / gridcell_dens_plus10cm;
			diam_g = pow(diam_g, 0.5);
		}

		outlimit_misc(out, out_height_plus10cm_sts, height);
		outlimit_misc(out, out_diam_plus10cm_sts, diam * 100.0);		//diameter output in cm
		outlimit_misc(out, out_diam_g_plus10cm_sts, diam_g * 100.0);	//diameter output in cm

		outlimit_misc(out, out_height_sts_patch_mean, gridcell_height_patch_mean);
		outlimit_misc(out, out_diam_sts_patch_mean, gridcell_diam_patch_mean * 100.0);		//diameter output in cm
		outlimit_misc(out, out_diam_g_sts_patch_mean, gridcell_diam_g_patch_mean * 100.0);	//diameter output in cm

		if (average_plus10cm_patches_only) {
			if (gridcell_area_frac_plus10cm) {
				gridcell_height_plus10cm_patch_mean /= gridcell_area_frac_plus10cm;
				gridcell_diam_plus10cm_patch_mean /= gridcell_area_frac_plus10cm;
				gridcell_diam_g_plus10cm_patch_mean /= gridcell_area_frac_plus10cm;
			}
		}

		outlimit_misc(out, out_height_plus10cm_sts_patch_mean, gridcell_height_plus10cm_patch_mean);
		outlimit_misc(out, out_diam_plus10cm_sts_patch_mean, gridcell_diam_plus10cm_patch_mean * 100.0);		//diameter output in cm
		outlimit_misc(out, out_diam_g_plus10cm_sts_patch_mean, gridcell_diam_g_plus10cm_patch_mean * 100.0);	//diameter output in cm

		if (run_landcover) {
			for (int lc = 0;lc < NLANDCOVERTYPES;lc++) {
				if (run[lc]) {

					if(landcover_storm_cmass_wood_save[lc])
						landcover_storm_sens_ind[lc] /= landcover_storm_cmass_wood_save[lc];	// cmass_wood weight
					else
						landcover_storm_sens_ind[lc] = 0.0;

					outlimit_misc(out, out_anpp_sts, landcover_anpp[lc]);
					outlimit_misc(out, out_agpp_sts, landcover_agpp[lc]);
					outlimit_misc(out, out_lai_sts, landcover_lai[lc]);
					outlimit_misc(out, out_cmass_sts, landcover_cmass[lc]);
//					outlimit_misc(out, out_cmass_potfuel_sts, landcover_cmass_potfuel[lc]);			// Decide later which output is needed here
					outlimit_misc(out, out_cmass_harv_killed_sts, landcover_cmass_harv_killed[lc]);
					outlimit_misc(out, out_cmass_thin_killed_sts, landcover_cmass_thin_killed[lc]);
					outlimit_misc(out, out_cmass_clearcut_killed_sts, landcover_cmass_clearcut_killed[lc]);
					outlimit_misc(out, out_cmass_fire_killed_sts, landcover_cmass_fire_killed[lc]);
					outlimit_misc(out, out_cmass_wind_killed_sts, landcover_cmass_wind_killed[lc]);
					outlimit_misc(out, out_cmass_barkbeetle_killed_sts, landcover_cmass_barkbeetle_killed[lc]);
					outlimit_misc(out, out_cmass_dist_clearcut_killed_sts, landcover_cmass_dist_clearcut_killed[lc]);
					outlimit_misc(out, out_cmass_dist_fire_killed_sts, landcover_cmass_dist_fire_killed[lc]);
					outlimit_misc(out, out_cmass_dist_wind_killed_sts, landcover_cmass_dist_wind_killed[lc]);
					outlimit_misc(out, out_cmass_dist_barkbeetle_killed_sts, landcover_cmass_dist_barkbeetle_killed[lc]);
					outlimit_misc(out, out_cmass_dist_prescribed_killed_sts, landcover_cmass_dist_clearcut_killed[lc] + 
						landcover_cmass_dist_fire_killed[lc] + landcover_cmass_dist_wind_killed[lc] + landcover_cmass_dist_barkbeetle_killed[lc]);
//					outlimit_misc(out, out_cmass_harv_tolitter_sts, landcover_cmass_harv_tolitter[lc]);
//					outlimit_misc(out, out_cmass_mort_sts, landcover_cmass_mort[lc]);
//					outlimit_misc(out, out_cmass_fire_sts, landcover_cmass_fire[lc]);
					outlimit_misc(out, out_cmass_dist_default_sts, landcover_cmass_dist[lc]);
//					outlimit_misc(out, out_cmass_leaf_root_turnover_sts, landcover_cmass_leaf_root_turnover[lc]);
//					outlimit_misc(out, out_cmass_repr_sts, landcover_cmass_repr[lc]);

					outlimit_misc(out, out_anpp_tree_sts, landcover_anpp_tree[lc]);
					outlimit_misc(out, out_agpp_tree_sts, landcover_agpp_tree[lc]);
					outlimit_misc(out, out_aautresp_tree_sts, (landcover_aautresp_root_tree[lc] + landcover_aautresp_sap_tree[lc] + landcover_aautresp_leaf_tree[lc] + landcover_aautresp_growth_tree[lc]));
					outlimit_misc(out, out_aautresp_root_tree_sts, landcover_aautresp_root_tree[lc]);
					outlimit_misc(out, out_aautresp_sap_tree_sts, landcover_aautresp_sap_tree[lc]);
					outlimit_misc(out, out_aautresp_leaf_tree_sts, landcover_aautresp_leaf_tree[lc]);
					outlimit_misc(out, out_aautresp_growth_tree_sts, landcover_aautresp_growth_tree[lc]);
					outlimit_misc(out, out_lai_tree_sts, landcover_lai_tree[lc]);
					outlimit_misc(out, out_cmass_tree_sts, landcover_cmass_tree[lc]);
					outlimit_misc(out, out_cmass_tree_mort_sts, landcover_cmass_tree_mort[lc]);

					outlimit_misc(out, out_cmass_wood_sts, landcover_cmass_wood[lc]);
//					outlimit_misc(out, out_cmass_wood_potharv_sts, landcover_cmass_wood_potharv[lc]);		// Decide later which output is needed here
//					outlimit_misc(out, out_cmass_wood_potharv_products_sts, landcover_cmass_wood_potharv_products[lc]);
					outlimit_misc(out, out_cmass_wood_harv_killed_sts, landcover_cmass_wood_harv_killed[lc]);
					outlimit_misc(out, out_cmass_wood_harv_sts, landcover_cmass_wood_harv[lc]);
					outlimit_misc(out, out_cmass_wood_harv_toprod_sts, landcover_cmass_wood_harv_toprod[lc]);
					outlimit_misc(out, out_cmass_wood_thin_sts, landcover_cmass_wood_thin[lc]);
					outlimit_misc(out, out_cmass_wood_clearcut_sts, landcover_cmass_wood_clearcut[lc]);
					outlimit_misc(out, out_cmass_wood_mort_sts, landcover_cmass_wood_mort[lc]);
					outlimit_misc(out, out_cmass_wood_inc_sts, landcover_cmass_wood_inc[lc]);

					outlimit_misc(out, out_csoil_sts, landcover_csoil[lc]);
					outlimit_misc(out, out_clitter_sts, landcover_clitter[lc]);
					outlimit_misc(out, out_clitter_course_sts, landcover_clitter_course[lc]);
					outlimit_misc(out, out_clitter_fine_sts, landcover_clitter_fine[lc]);
					outlimit_misc(out, out_rh_sts, landcover_rh[lc]);
					outlimit_misc(out, out_rh_clitter_course_sts, landcover_rh_clitter_course[lc]);
					outlimit_misc(out, out_rh_clitter_fine_sts, landcover_rh_clitter_fine[lc]);
					outlimit_misc(out, out_csink_sts, landcover_csink[lc]);
					outlimit_misc(out, out_luc_cflux_sts, landcover_luc_cflux[lc]);
					outlimit_misc(out, out_nsoil_sts, landcover_nsoil[lc]);
					outlimit_misc(out, out_nlitter_sts, landcover_nlitter[lc]);
					outlimit_misc(out, out_nlitter_course_sts, landcover_nlitter_course[lc]);
					outlimit_misc(out, out_nlitter_fine_sts, landcover_nlitter_fine[lc]);

					outlimit_misc(out, out_runoff_sts, landcover_runoff[lc]);
					outlimit_misc(out, out_transpiration_sts, landcover_aet[lc]);
					outlimit_misc(out, out_intercep_sts, landcover_intercep[lc]);
					outlimit_misc(out, out_evap_sts, landcover_evap[lc]);

					outlimit_misc(out, out_storm_si_sts, landcover_storm_sens_ind[lc] * 1000.0);	// Storm sensitivity * 1000
					if(landcover_storm_cmass_wood_save[lc])
						outlimit_misc(out, out_storm_frac_pot_sts, landcover_storm_dam_wood_pot[lc] / landcover_storm_cmass_wood_save[lc] * 100.0);	// Prognostic damaged fraction in %
					else
						outlimit_misc(out, out_storm_frac_pot_sts, 0.0);
					if(landcover_storm_cmass_wood_save[lc])
						outlimit_misc(out, out_storm_frac_sts, landcover_storm_dam_wood[lc] / landcover_storm_cmass_wood_save[lc] * 100.0);	// Actual damaged fraction in %
					else
						outlimit_misc(out, out_storm_frac_sts, 0.0);
					outlimit_misc(out, out_storm_dw_pot_sts, landcover_storm_dam_wood_pot[lc]);
					outlimit_misc(out, out_storm_dw_sts, landcover_storm_dam_wood[lc]);
					outlimit_misc(out, out_insect_dw_sts, landcover_insect_dam_wood[lc]);
					outlimit_misc(out, out_cmort_storm_sts, landcover_cmass_wind_killed[lc] + landcover_cmass_dist_wind_killed[lc]);
					outlimit_misc(out, out_cmort_insect_sts, landcover_cmass_barkbeetle_killed[lc] + landcover_cmass_dist_barkbeetle_killed[lc]);
					outlimit_misc(out, out_cmort_fire_sts, landcover_cmass_fire_killed[lc] + landcover_cmass_dist_fire_killed[lc]);
					outlimit_misc(out, out_cmort_clearcut_sts, landcover_cmass_clearcut_killed[lc] + landcover_cmass_dist_clearcut_killed[lc]);

					outlimit_misc(out, out_dens_sts, landcover_dens[lc]);

					outlimit_misc(out, out_frac_sts, gridcell.landcover.frac[lc]);

					outlimit_misc(out, out_nstand_sts, landcover_nstands[lc]);
					outlimit_misc(out, out_npatch_sts, landcover_npatches[lc]);

					double height = 0.0;
					double diam = 0.0;
					double diam_g = 0.0;
					if (landcover_dens[lc]) {
						height = landcover_height[lc] / landcover_dens[lc];
						diam = landcover_diam[lc] / landcover_dens[lc];
						diam_g = landcover_diam_g[lc] / landcover_dens[lc];
						diam_g = pow(diam_g, 0.5);
					}

					outlimit_misc(out, out_height_sts, height);
					outlimit_misc(out, out_diam_sts, diam * 100.0);		//diameter output in cm
					outlimit_misc(out, out_diam_g_sts, diam_g * 100.0);	//diameter output in cm

					outlimit_misc(out, out_cmass_plus10cm_sts, landcover_cmass_plus10cm[lc]);
					outlimit_misc(out, out_cmass_harv_killed_plus10cm_sts, landcover_cmass_harv_killed_plus10cm[lc]);
					outlimit_misc(out, out_cmass_dist_clearcut_killed_plus10cm_sts, landcover_cmass_dist_clearcut_killed_plus10cm[lc]);
					outlimit_misc(out, out_cmass_dist_fire_killed_plus10cm_sts, landcover_cmass_dist_fire_killed_plus10cm[lc]);
					outlimit_misc(out, out_cmass_dist_wind_killed_plus10cm_sts, landcover_cmass_dist_wind_killed_plus10cm[lc]);
					outlimit_misc(out, out_cmass_dist_barkbeetle_killed_plus10cm_sts, landcover_cmass_dist_barkbeetle_killed_plus10cm[lc]);

					outlimit_misc(out, out_cmass_wood_plus10cm_sts, landcover_cmass_wood_plus10cm[lc]);
					outlimit_misc(out, out_cmass_wood_harv_killed_plus10cm_sts, landcover_cmass_wood_harv_killed_plus10cm[lc]);
					outlimit_misc(out, out_cmass_wood_mort_plus10cm_sts, landcover_cmass_wood_mort_plus10cm[lc]);
					outlimit_misc(out, out_cmass_wood_inc_plus10cm_sts, landcover_cmass_wood_inc_plus10cm[lc]);
					
					outlimit_misc(out, out_dens_plus10cm_sts, landcover_dens_plus10cm[lc]);

					height = 0.0;
					diam = 0.0;
					diam_g = 0.0;
					if (landcover_dens_plus10cm[lc]) {
						height = landcover_height_plus10cm[lc] / landcover_dens_plus10cm[lc];
						diam = landcover_diam_plus10cm[lc] / landcover_dens_plus10cm[lc];
						diam_g = landcover_diam_g_plus10cm[lc] / landcover_dens_plus10cm[lc];
						diam_g = pow(diam_g, 0.5);
					}

					outlimit_misc(out, out_height_plus10cm_sts, height);
					outlimit_misc(out, out_diam_plus10cm_sts, diam * 100.0);	//diameter output in cm
					outlimit_misc(out, out_diam_g_plus10cm_sts, diam_g * 100.0);	//diameter output in cm

					outlimit_misc(out, out_height_sts_patch_mean, landcover_height_patch_mean[lc]);
					outlimit_misc(out, out_diam_sts_patch_mean, landcover_diam_patch_mean[lc] * 100.0);		//diameter output in cm
					outlimit_misc(out, out_diam_g_sts_patch_mean, landcover_diam_g_patch_mean[lc] * 100.0);	//diameter output in cm

					if (average_plus10cm_patches_only) {
						if (landcover_area_frac_plus10cm[lc]) {
							landcover_height_plus10cm_patch_mean[lc] /= landcover_area_frac_plus10cm[lc];
							landcover_diam_plus10cm_patch_mean[lc] /= landcover_area_frac_plus10cm[lc];
							landcover_diam_g_plus10cm_patch_mean[lc] /= landcover_area_frac_plus10cm[lc];
						}
					}
					outlimit_misc(out, out_height_plus10cm_sts_patch_mean, landcover_height_plus10cm_patch_mean[lc]);
					outlimit_misc(out, out_diam_plus10cm_sts_patch_mean, landcover_diam_plus10cm_patch_mean[lc] * 100.0);		//diameter output in cm
					outlimit_misc(out, out_diam_g_plus10cm_sts_patch_mean, landcover_diam_g_plus10cm_patch_mean[lc] * 100.0);	//diameter output in cm
				}
			}
		}

		outlimit_misc(out, out_storm, gridcell_storm_sens_ind * 1000.0); // Sensitivity * 1000
		if(gridcell_storm_cmass_wood_save)
			outlimit_misc(out, out_storm, gridcell_storm_dam_wood_pot / gridcell_storm_cmass_wood_save * 100.0); // Potential damaged fraction in %
		else
			outlimit_misc(out, out_storm, 0.0);
		if(gridcell_storm_cmass_wood_save)
			outlimit_misc(out, out_storm, gridcell_storm_dam_wood / gridcell_storm_cmass_wood_save * 100.0); // Actual damaged fraction in %
		else
			outlimit_misc(out, out_storm, 0.0);
		outlimit_misc(out, out_storm, gridcell_storm_dam_wood_pot);
		outlimit_misc(out, out_storm, gridcell_storm_dam_wood);

		double insect_damwood = 0.0;
		double insect_gc_popi = 0.0;
		for (int i = 0; i < nift; i++) {
			insect_damwood += gridcell.ift[i].insect_damage;
			insect_gc_popi += gridcell.ift[i].pop_ind;
		}
		outlimit_misc(out, out_storm, insect_damwood);
		outlimit_misc(out, out_storm, insect_gc_popi);
	}

	outlimit_misc(out, out_cflux_harvest, gridcell.cflux_harv_cc);
	outlimit_misc(out, out_cflux_harvest, gridcell.cflux_harv_thin);
	outlimit_misc(out, out_cflux_harvest, gridcell.cflux_harv_cc_prescribed);
	outlimit_misc(out, out_cflux_harvest, gridcell.cflux_harv_wind_prescribed);
	outlimit_misc(out, out_cflux_harvest, gridcell.cflux_harv_barkbeetle_prescribed);
	
	outlimit_misc(out, out_harvest_volume, gridcell.cfluxvol_harv_cc);
	outlimit_misc(out, out_harvest_volume, gridcell.cfluxvol_harv_thin);
	outlimit_misc(out, out_harvest_volume, gridcell.cfluxvol_harv_cc_prescribed);
	outlimit_misc(out, out_harvest_volume, gridcell.cfluxvol_harv_wind_prescribed);
	outlimit_misc(out, out_harvest_volume, gridcell.cfluxvol_harv_barkbeetle_prescribed);

	Landcover& lcC = gridcell.landcover;

	// Print forest vegetation total, wood and AG compartments

	// Primary forest

	// Total vegetation cmass
	outlimit_misc(out, out_forest_vegc, landcover_cmass[NATURAL] * lcC.frac[NATURAL]);
	// Potential harvestable stem wood, taking harvest efficiency into account
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_natural);
	// Potential wood products, taking harvest efficiency into account
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_products_natural);
	// Potential fuel wood, taking harvest efficiency and residue outtake into account
	outlimit_misc(out, out_forest_vegc, cmass_potfuel_natural);

	// Managed forest

	outlimit_misc(out, out_forest_vegc, landcover_cmass[FOREST] * lcC.frac[FOREST]);
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_forest);
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_products_forest);
	outlimit_misc(out, out_forest_vegc, cmass_potfuel_forest);

	// Total forest

	outlimit_misc(out, out_forest_vegc, landcover_cmass[NATURAL] * lcC.frac[NATURAL] + landcover_cmass[FOREST] * lcC.frac[FOREST]);
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_natural + cmass_wood_potharv_forest);
	outlimit_misc(out, out_forest_vegc, cmass_wood_potharv_products_natural + cmass_wood_potharv_products_forest);
	outlimit_misc(out, out_forest_vegc, cmass_potfuel_natural + cmass_potfuel_forest);

	// Print wood harvest total killed and fate of compartments (All harvest events using harvest_forest(), including prescribed clearcut, storm and insect mortality)

	// Primary forest harvest

	// Killed tree cmass during harvest
	outlimit_misc(out, out_forest_harvest, lcC.cmass_harv_killed + cmass_harv_killed_natural);
	// Harvested stem wood cmass; wood products and fuelwood from stems
	outlimit_misc(out, out_forest_harvest, lcC.cmass_stem_harvest + cmass_wood_harv_natural);
	// Harvested wood product cmass
	outlimit_misc(out, out_forest_harvest, lcC.cmass_stem_toprod + cmass_wood_harv_toprod_natural);
	// Fuel wood cmass from products and residues
	outlimit_misc(out, out_forest_harvest, lcC.acflux_wood_harvest + flux_charvest_lc[NATURAL]);
	// Killed tree cmass entering litter pool
	outlimit_misc(out, out_forest_harvest, lcC.cmass_harv_tolitter + cmass_harv_tolitter_natural);

	// Managed forest harvest

	outlimit_misc(out, out_forest_harvest, cmass_harv_killed_forest);
	outlimit_misc(out, out_forest_harvest, cmass_wood_harv_forest);
	outlimit_misc(out, out_forest_harvest, cmass_wood_harv_toprod_forest);
	outlimit_misc(out, out_forest_harvest, flux_charvest_lc[FOREST]);
	outlimit_misc(out, out_forest_harvest, cmass_harv_tolitter_forest);

	// Total forest harvest

	outlimit_misc(out, out_forest_harvest, lcC.cmass_harv_killed + cmass_harv_killed_natural + cmass_harv_killed_forest);
	outlimit_misc(out, out_forest_harvest, lcC.cmass_stem_harvest + cmass_wood_harv_natural + cmass_wood_harv_forest);
	outlimit_misc(out, out_forest_harvest, lcC.cmass_stem_toprod + cmass_wood_harv_toprod_natural + cmass_wood_harv_toprod_forest);
	outlimit_misc(out, out_forest_harvest, lcC.acflux_wood_harvest + flux_charvest_lc[NATURAL] + flux_charvest_lc[FOREST]);
	outlimit_misc(out, out_forest_harvest, lcC.cmass_harv_tolitter + cmass_harv_tolitter_natural + cmass_harv_tolitter_forest);

	double anpp_natural = landcover_anpp[NATURAL] * lcC.frac[NATURAL];	// identical values to sum of st_anpp values
	double anpp_forest = landcover_anpp[FOREST] * lcC.frac[FOREST];

	double cflux_veg_natural = lcC.acflux_cloned_lc[NATURAL] - anpp_natural + lcC.cmass_harv_killed + cmass_harv_killed_natural 
		+ cmass_mort_natural + cmass_fire_natural + cmass_est_natural + cmass_dist_natural + cmass_leaf_root_turnover_natural 
		+ cmass_repr_natural;
	double cflux_veg_forest = lcC.acflux_cloned_lc[FOREST] - anpp_forest + cmass_harv_killed_forest + cmass_mort_forest 
		+ cmass_fire_forest + cmass_est_forest + cmass_dist_forest + cmass_leaf_root_turnover_forest + cmass_repr_forest;
	double cflux_veg_tot = cflux_veg_forest + cflux_veg_natural;
	double NAI_natural = -(cflux_veg_natural - lcC.acflux_cloned_lc[NATURAL] - lcC.cmass_harv_killed - cmass_harv_killed_natural);
	double NAI_forest = -(cflux_veg_forest - lcC.acflux_cloned_lc[FOREST] - cmass_harv_killed_forest);
	double NAI_tot = NAI_forest + NAI_natural;

	// Print C fluxes to and from vegetation,  (harvC includes all harvest events using harvest_forest(), including prescribed clearcut, storm and insect mortality)
	// Primary forest
	outlimit_misc(out, out_forest_cflux_veg, -anpp_natural);
	outlimit_misc(out, out_forest_cflux_veg, lcC.cmass_harv_killed + cmass_harv_killed_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_mort_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_fire_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_est_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_dist_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_leaf_root_turnover_natural);
	outlimit_misc(out, out_forest_cflux_veg, cmass_repr_natural);
	outlimit_misc(out, out_forest_cflux_veg, lcC.acflux_cloned_lc[NATURAL]);
	outlimit_misc(out, out_forest_cflux_veg, cflux_veg_natural);
	outlimit_misc(out, out_forest_cflux_veg, NAI_natural);
	// Managed forest
	outlimit_misc(out, out_forest_cflux_veg, -anpp_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_harv_killed_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_mort_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_fire_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_est_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_dist_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_leaf_root_turnover_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_repr_forest);
	outlimit_misc(out, out_forest_cflux_veg, lcC.acflux_cloned_lc[FOREST]);
	outlimit_misc(out, out_forest_cflux_veg, cflux_veg_forest);
	outlimit_misc(out, out_forest_cflux_veg, NAI_forest);
	// Total forest
	outlimit_misc(out, out_forest_cflux_veg, -anpp_natural - anpp_forest);
	outlimit_misc(out, out_forest_cflux_veg, lcC.cmass_harv_killed + cmass_harv_killed_natural + cmass_harv_killed_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_mort_natural + cmass_mort_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_fire_natural + cmass_fire_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_est_natural + cmass_est_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_dist_natural + cmass_dist_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_leaf_root_turnover_natural + cmass_leaf_root_turnover_forest);
	outlimit_misc(out, out_forest_cflux_veg, cmass_repr_natural + cmass_repr_forest);
	outlimit_misc(out, out_forest_cflux_veg, lcC.acflux_cloned_lc[NATURAL] + lcC.acflux_cloned_lc[FOREST]);
	outlimit_misc(out, out_forest_cflux_veg, cflux_veg_tot);
	outlimit_misc(out, out_forest_cflux_veg, NAI_tot);

	if (run[CROPLAND]) {
		outlimit_misc(out, out_irrigation,   irrigation_gridcell);
	}

	// Print landcover totals to files
	if (run_landcover) {
		for (int i=0;i<NLANDCOVERTYPES;i++) {
			if (run[i]) {

				outlimit_misc(out, out_dens_plus10cm, landcover_densindiv_total_plus10cm[i]);

				switch (i) {
				case CROPLAND:
					if (run[CROPLAND]) {
						outlimit_misc(out, out_anpp_cropland,			landcover_anpp[i]);
						outlimit_misc(out, out_cmass_cropland,		landcover_cmass[i]);
						outlimit_misc(out, out_lai_cropland,			landcover_lai[i]);
					}
					break;
				case PASTURE:
					if (run[NATURAL]) {
						outlimit_misc(out, out_anpp_active_pasture,					landcover_anpp[i]);
						outlimit_misc(out, out_cmass_active_pasture,					landcover_cmass[i]);
						outlimit_misc(out, out_anpp_pasture,			landcover_anpp[i]);
						outlimit_misc(out, out_cmass_pasture,			landcover_cmass[i]);
						outlimit_misc(out, out_lai_pasture,			landcover_lai[i]);
					}
					break;
				case BARREN:
					break;
				case NATURAL:
					// if(run[FOREST] || run[PASTURE]) {
					if(run[FOREST]) {
						outlimit_misc(out, out_anpp_active_natural,					landcover_anpp[i]);
						outlimit_misc(out, out_cmass_active_natural,					landcover_cmass[i]);
						outlimit_misc(out, out_dens_natural,					landcover_densindiv_total[i]);
						outlimit_misc(out, out_anpp_natural,			landcover_anpp[i]);
						outlimit_misc(out, out_cmass_natural,			landcover_cmass[i]);
						outlimit_misc(out, out_lai_natural,			landcover_lai[i]);
						outlimit_misc(out, out_aaet_active_natural,					landcover_aaet[i]);
						outlimit_misc(out, out_lai_active_natural,						landcover_lai[i]);
						outlimit_misc(out, out_fpc_active_natural,						landcover_fpc[i]);

						outlimit_misc(out, out_anpp_gridcell_natural,    landcover_anpp[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_cmass_gridcell_natural,	landcover_cmass[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_aaet_gridcell_natural,	landcover_aaet[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_lai_gridcell_natural,		landcover_lai[i] * gridcell.landcover.frac[i]);
					}
					outlimit_misc(out, out_frac_natural, gridcell.landcover.frac[i]);
					break;
				case FOREST:
					if (run[NATURAL]) {
						outlimit_misc(out, out_anpp_active_forest,						landcover_anpp[i]);
						outlimit_misc(out, out_cmass_active_forest,					landcover_cmass[i]);
						outlimit_misc(out, out_dens_forest,						landcover_densindiv_total[i]);
						outlimit_misc(out, out_anpp_forest,			landcover_anpp[i]);
						outlimit_misc(out, out_cmass_forest,			landcover_cmass[i]);
						outlimit_misc(out, out_lai_forest,			landcover_lai[i]);
						outlimit_misc(out, out_aaet_active_forest,						landcover_aaet[i]);
						outlimit_misc(out, out_lai_active_forest,						landcover_lai[i]);
						outlimit_misc(out, out_fpc_active_forest,						landcover_fpc[i]);

						outlimit_misc(out, out_anpp_gridcell_forest,     landcover_anpp[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_cmass_gridcell_forest,	landcover_cmass[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_aaet_gridcell_forest,		landcover_aaet[i] * gridcell.landcover.frac[i]);
						outlimit_misc(out, out_lai_gridcell_forest,		landcover_lai[i] * gridcell.landcover.frac[i]);
					}
					outlimit_misc(out, out_frac_forest, gridcell.landcover.frac[i]);
					break;
				case URBAN:
					break;
				case PEATLAND:
					if (run[PEATLAND]) {
						outlimit_misc(out, out_anpp_active_peatland,					landcover_anpp[i]);
						outlimit_misc(out, out_cmass_active_peatland,					landcover_cmass[i]);
						outlimit_misc(out, out_anpp_peatland,			landcover_anpp[i]);
						outlimit_misc(out, out_cmass_peatland,		landcover_cmass[i]);
						outlimit_misc(out, out_lai_peatland,			landcover_lai[i]);
					}
					break;
				default:
					if (date.year == nyear_spinup)
						dprintf("Modify code to deal with landcover output!\n");
				}

				// N2O info output per landcover
				outlimit_misc(out, out_soil_n2o_flux, flux_N2O_soil[i] * M2_PER_HA);
			}
		}
		// exclude BARREN, LUH2 BARREN is ice/water fraction, but BARREN produces n-fluxes on barren
		// TODO: are there N-fluxes from water? we really need an extra landcover ice/water.
		outlimit_misc(out, out_soil_n2o_flux, (flux_N2O_soil_total - flux_N2O_soil[BARREN]) * M2_PER_HA);
		outlimit_misc(out, out_soil_n2o_flux, flux_N2O_soil_total * M2_PER_HA);
	}

	// Print C fluxes of forests (typically PNV, so potentially grassland) with LUC history
	outlimit_misc(out, out_cflux_forestry, flux_veg_forestry);
	outlimit_misc(out, out_cflux_forestry, -flux_repr_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_soil_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_fire_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_est_forestry);
	outlimit_misc(out, out_cflux_forestry, c_org_leach_gridcell_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_seed_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_charvest_forestry);
	outlimit_misc(out, out_cflux_forestry, flux_veg_forestry - flux_repr_forestry + flux_soil_forestry + flux_fire_forestry 
		+ flux_est_forestry + c_org_leach_gridcell_forestry + flux_seed_forestry);

	outlimit_misc(out, out_cflux_regrowth, flux_veg_regrowth);
	outlimit_misc(out, out_cflux_regrowth, -flux_repr_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_soil_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_fire_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_est_regrowth);
	outlimit_misc(out, out_cflux_regrowth, c_org_leach_gridcell_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_seed_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_charvest_regrowth);
	outlimit_misc(out, out_cflux_regrowth, flux_veg_regrowth - flux_repr_regrowth + flux_soil_regrowth + flux_fire_regrowth 
		+ flux_est_regrowth + c_org_leach_gridcell_regrowth + flux_seed_regrowth);

	outlimit_misc(out, out_cflux_primary, flux_veg_primary);
	outlimit_misc(out, out_cflux_primary, -flux_repr_primary);
	outlimit_misc(out, out_cflux_primary, flux_soil_primary);
	outlimit_misc(out, out_cflux_primary, flux_fire_primary);
	outlimit_misc(out, out_cflux_primary, flux_est_primary);
	outlimit_misc(out, out_cflux_primary, c_org_leach_gridcell_primary);
	outlimit_misc(out, out_cflux_primary, flux_seed_primary);
	outlimit_misc(out, out_cflux_primary, flux_charvest_primary);
	outlimit_misc(out, out_cflux_primary, flux_veg_primary - flux_repr_primary + flux_soil_primary + flux_fire_primary 
		+ flux_est_primary + c_org_leach_gridcell_primary + flux_seed_primary);

	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_wood_harvest_orig);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_wood_harvest_orig - gridcell.landcover.acflux_wood_harvest);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_clearing_orig);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_clearing_orig - gridcell.landcover.acflux_clearing);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_landuse_change_orig);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_landuse_change_orig - gridcell.landcover.acflux_landuse_change);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_harvest_slow);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_wood_harvest_orig - gridcell.landcover.acflux_wood_harvest 
		+ gridcell.landcover.acflux_clearing_orig - gridcell.landcover.acflux_clearing + gridcell.landcover.acflux_landuse_change_orig 
		- gridcell.landcover.acflux_landuse_change + gridcell.landcover.acflux_harvest_slow);
	outlimit_misc(out, out_harvest_flux_luc, gridcell.landcover.acflux_wood_harvest + gridcell.landcover.acflux_clearing_orig 
		+ gridcell.landcover.acflux_landuse_change + gridcell.landcover.acflux_harvest_slow);

	// Print C fluxes to per-landcover files
	if (run_landcover) {
		for (int i=0;i<NLANDCOVERTYPES;i++) {
			if (run[i]) {

				GuessOutput::Table* table_p=NULL;
				GuessOutput::Table* table_p_N=NULL;

				GuessOutput::Table* table_p_N_soil=NULL;
				switch (i) {
				case CROPLAND:
					table_p=&out_cflux_cropland;
					table_p_N=&out_nflux_cropland;
					table_p_N_soil=&out_soil_nflux_cropland;
					break;
				case PASTURE:
					table_p=&out_cflux_pasture;
					table_p_N=&out_nflux_pasture;
					table_p_N_soil=&out_soil_nflux_pasture;
					break;
				case NATURAL:
					table_p=&out_cflux_natural;
					table_p_N=&out_nflux_natural;
					table_p_N_soil=&out_soil_nflux_natural;
					break;
				case FOREST:
					table_p=&out_cflux_forest;
					table_p_N=&out_nflux_forest;
					table_p_N_soil=&out_soil_nflux_forest;
					break;
				case URBAN:
					break;
				case PEATLAND:
					table_p=&out_cflux_peatland;
					table_p_N=&out_nflux_peatland;
					break;
				case BARREN:
					break;
				default:
					if (date.year == nyear_spinup)
						dprintf("Modify code to deal with landcover output!\n");
				}

				if (table_p) {
					outlimit_misc(out, *table_p, flux_veg_lc[i]);
					outlimit_misc(out, *table_p, -flux_repr_lc[i]);
					outlimit_misc(out, *table_p, flux_soil_lc[i] + c_org_leach_lc[i]);
					outlimit_misc(out, *table_p, flux_fire_lc[i]);
					outlimit_misc(out, *table_p, flux_est_lc[i]);

					if (run_landcover) {
						outlimit_misc(out, *table_p, flux_seed_lc[i]);
						outlimit_misc(out, *table_p, flux_charvest_lc[i]);
						outlimit_misc(out, *table_p, lc.acflux_wood_harvest_lc[i] + lc.acflux_clearing_lc[i] 
							+ lc.acflux_landuse_change_lc[i]);
						outlimit_misc(out, *table_p, lc.acflux_harvest_slow_lc[i]);
					}
				}

				if (table_p_N) {
					outlimit_misc(out, *table_p_N, -andep_lc[i] * M2_PER_HA);
					outlimit_misc(out, *table_p_N, -anfix_lc[i] * M2_PER_HA);
					outlimit_misc(out, *table_p_N, -anfert_lc[i] * M2_PER_HA);
					outlimit_misc(out, *table_p_N, flux_ntot_lc[i] * M2_PER_HA);
					outlimit_misc(out, *table_p_N, (n_min_leach_lc[i] + n_org_leach_lc[i]) * M2_PER_HA);

					if (run_landcover) {
						outlimit_misc(out, *table_p_N, flux_nseed_lc[i] * M2_PER_HA);
						outlimit_misc(out, *table_p_N, flux_nharvest_lc[i] * M2_PER_HA);
						outlimit_misc(out, *table_p_N, (lc.anflux_wood_harvest_lc[i] + lc.anflux_clearing_lc[i] 
						+ gridcell.landcover.anflux_landuse_change_lc[i]) * M2_PER_HA);
						outlimit_misc(out, *table_p_N, lc.anflux_harvest_slow_lc[i] * M2_PER_HA);
					}
				}

				double cflux_total = flux_veg_lc[i] - flux_repr_lc[i] + flux_soil_lc[i] + flux_fire_lc[i] 
					+ flux_est_lc[i] + c_org_leach_lc[i];
				double nflux_total = -andep_lc[i] - anfix_lc[i] - anfert_lc[i] + flux_ntot_lc[i] + n_min_leach_lc[i] 
					+ n_org_leach_lc[i];

				if (run_landcover) {
					cflux_total += flux_seed_lc[i];
					cflux_total += flux_charvest_lc[i];
					cflux_total += lc.acflux_wood_harvest_lc[i];
					cflux_total += lc.acflux_clearing_lc[i];
					cflux_total += lc.acflux_landuse_change_lc[i];
					cflux_total += lc.acflux_harvest_slow_lc[i];
					nflux_total += flux_nseed_lc[i];
					nflux_total += flux_nharvest_lc[i];
					nflux_total += lc.anflux_wood_harvest_lc[i];
					nflux_total += lc.anflux_clearing_lc[i];
					nflux_total += lc.anflux_landuse_change_lc[i];
					nflux_total += lc.anflux_harvest_slow_lc[i];
				}
				if (table_p) {
					outlimit_misc(out, *table_p,  cflux_total);
				}
				if (table_p_N) {
					outlimit_misc(out, *table_p_N,  nflux_total * M2_PER_HA);
				}
			}
		}

		// Print C pools of forests with LUC history
		outlimit_misc(out, out_cpool_forestry, cmass_gridcell_forestry);
		outlimit_misc(out, out_cpool_forestry, clitter_gridcell_forestry + surfsoillitterc_forestry + cwdc_forestry);
		outlimit_misc(out, out_cpool_forestry, centuryc_forestry);
		outlimit_misc(out, out_cpool_forestry, cmass_gridcell_forestry + clitter_gridcell_forestry + centuryc_forestry 
			+ surfsoillitterc_forestry + cwdc_forestry);

		outlimit_misc(out, out_cpool_regrowth, cmass_gridcell_regrowth);
		outlimit_misc(out, out_cpool_regrowth, clitter_gridcell_regrowth + surfsoillitterc_regrowth + cwdc_regrowth);
		outlimit_misc(out, out_cpool_regrowth, centuryc_regrowth);
		outlimit_misc(out, out_cpool_regrowth, cmass_gridcell_regrowth + clitter_gridcell_regrowth + centuryc_regrowth 
			+ surfsoillitterc_regrowth + cwdc_regrowth);

		outlimit_misc(out, out_cpool_primary, cmass_gridcell_primary);
		outlimit_misc(out, out_cpool_primary, clitter_gridcell_primary + surfsoillitterc_primary + cwdc_primary);
		outlimit_misc(out, out_cpool_primary, centuryc_primary);
		outlimit_misc(out, out_cpool_primary, cmass_gridcell_primary + clitter_gridcell_primary + centuryc_primary 
			+ surfsoillitterc_primary + cwdc_primary);

		// Print C pools to per-landcover files
		for (int i=0;i<NLANDCOVERTYPES;i++) {
			if (run[i]) {

				GuessOutput::Table* table_p=NULL;
				GuessOutput::Table* table_p_N=NULL;

				switch (i) {
				case CROPLAND:
					table_p=&out_cpool_cropland;
					table_p_N=&out_npool_cropland;
					break;
				case PASTURE:
					table_p=&out_cpool_pasture;
					table_p_N=&out_npool_pasture;
					break;
				case NATURAL:
					table_p=&out_cpool_natural;
					table_p_N=&out_npool_natural;
					break;
				case FOREST:
					table_p=&out_cpool_forest;
					table_p_N=&out_npool_forest;
					break;
				case URBAN:
					break;
				case PEATLAND:
					table_p=&out_cpool_peatland;
					table_p_N=&out_npool_peatland;
					break;
				case BARREN:
					break;
				default:
					if (date.year == nyear_spinup)
						dprintf("Modify code to deal with landcover output!\n");
				}

				if (table_p) {
					outlimit_misc(out, *table_p, landcover_cmass[i] * lc.frac[i]);
					outlimit_misc(out, *table_p_N, (landcover_nmass[i] + landcover_nlitter[i]) * lc.frac[i]);

					if (!ifcentury) {
						outlimit_misc(out, *table_p, landcover_clitter[i] * lc.frac[i]);
						outlimit_misc(out, *table_p, c_fast_lc[i]);
						outlimit_misc(out, *table_p, c_slow_lc[i]);
					}
					else {
						outlimit_misc(out, *table_p, landcover_clitter[i] * lc.frac[i] + surfsoillitterc_lc[i] + cwdc_lc[i]);
						outlimit_misc(out, *table_p, centuryc_lc[i]);
						outlimit_misc(out, *table_p_N, surfsoillittern_lc[i] + cwdn_lc[i]);
						outlimit_misc(out, *table_p_N, centuryn_lc[i] + availn_lc[i]);
					}

					if (run_landcover && ifslowharvestpool) {
						outlimit_misc(out, *table_p, c_harv_slow_lc[i]);
						outlimit_misc(out, *table_p_N, n_harv_slow_lc[i]);
					}
				}

				// Calculate total cpool, starting with cmass and litter...
				double cpool_total = (landcover_cmass[i] + landcover_clitter[i]) * lc.frac[i];
				double npool_total = (landcover_nmass[i] + landcover_nlitter[i]) * lc.frac[i];

				// Add SOM pools
				if (!ifcentury) {
					cpool_total += c_fast_lc[i] + c_slow_lc[i];
				}
				else {
					cpool_total += centuryc_lc[i] + surfsoillitterc_lc[i] + cwdc_lc[i];
					npool_total += centuryn_lc[i] + surfsoillittern_lc[i] + cwdn_lc[i] + availn_lc[i];
				}

				// Add slow harvest pool if needed
				if (run_landcover && ifslowharvestpool) {
					cpool_total += c_harv_slow_lc[i];
					npool_total += n_harv_slow_lc[i];
				}
				if (table_p) {
					outlimit_misc(out, *table_p, cpool_total);
					outlimit_misc(out, *table_p_N, npool_total);
				}
			}
		}
	}

	double gridcell_climate_agdd0_20_mean;
	if (gridcell.climate.agdd0_20.size()>0) {
		gridcell_climate_agdd0_20_mean = gridcell.climate.agdd0_20.mean();
	}
	else {
		gridcell_climate_agdd0_20_mean = -999;
	}

	//Output of seasonality variables
	outlimit_misc(out, out_seasonality,   gridcell.climate.seasonality);
	outlimit_misc(out, out_seasonality,   gridcell.climate.var_temp);
	outlimit_misc(out, out_seasonality,   gridcell.climate.var_prec);
	outlimit_misc(out, out_seasonality,   gridcell.climate.mtemp_min20);
	outlimit_misc(out, out_seasonality,   gridcell.climate.atemp_mean);
	outlimit_misc(out, out_seasonality,   gridcell.climate.mtemp_max20);
	outlimit_misc(out, out_seasonality,   gridcell.climate.mtemp_max);
	outlimit_misc(out, out_seasonality,   gridcell.climate.temp_seasonality);
	outlimit_misc(out, out_seasonality,   gridcell_climate_agdd0_20_mean);
	outlimit_misc(out, out_seasonality,   gridcell.climate.agdd5);
	outlimit_misc(out, out_seasonality,   gridcell.climate.mprec_petmin20);
	outlimit_misc(out, out_seasonality,   gridcell.climate.aprec);
	outlimit_misc(out, out_seasonality,   gridcell.climate.prec_range);
	outlimit_misc(out, out_seasonality,   gridcell.climate.ainsol/date.year_length());

	if(st_pft_cmass)
		delete[] st_pft_cmass;
	if(st_total_cmass)
		delete[] st_total_cmass;
	if(st_pft_cmass_harv_killed)
		delete[] st_pft_cmass_harv_killed;
	if(st_total_cmass_harv_killed)
		delete[] st_total_cmass_harv_killed;

	// from EC-Earth
	// Monthly output of lambda
	for (int m = 0; m < 12; m++) {

		OutputRows out(output_channel, lon, lat, date.get_calendar_year(), m + 1);

		for (int q = 0; q < npft; q++) {
			outlimit_misc(out, out_mlambda, !negligible(pft_mgpp[m][q]) ? pft_mlambda_gpp[m][q] / pft_mgpp[m][q] : 0.0);
		}
	}

}

/// Output of simulation results at the end of each day
/** This function does not have to provide any information to the framework.
  */
void MiscOutput::outdaily(Gridcell& gridcell) {

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	OutputRows out(output_channel, lon, lat, date.get_calendar_year(), date.day);

	if (date.year < nyear_spinup) {
		return;
	}

	outlimit_misc(out, out_daily_climate, gridcell.climate.temp);
	outlimit_misc(out, out_daily_climate, gridcell.climate.prec);
	outlimit_misc(out, out_daily_climate, gridcell.climate.rad);

	if (firemodel==BLAZE) {
		outlimit_misc(out, out_daily_blaze, (int)gridcell.dfire_die_roll);
		outlimit_misc(out, out_daily_blaze, gridcell.dfire_probability);
		outlimit_misc(out, out_daily_blaze, gridcell.effective_burned_area);
		outlimit_misc(out, out_daily_blaze, gridcell.climate.kbdi);
		outlimit_misc(out, out_daily_blaze, gridcell.climate.mcarthur_forest_fire_index);
		outlimit_misc(out, out_daily_blaze, gridcell.climate.max_potfli);
		outlimit_misc(out, out_daily_blaze, gridcell.climate.mean_potfli);
	}
}

void MiscOutput::openlocalfiles(Gridcell& gridcell, int coordinates_precision) {

	if(!printseparatestands)
		return;

	char dirname[200]={'\0'};
#ifdef _MSC_VER
	strcpy(dirname, "stand_output\\");
#else
	strcpy(dirname, "../stand_output/");
#endif
	std::filesystem::create_directories(dirname);

	if (!date.year || restart && date.year == state_year) {
		for(int id=0;id<MAXNUMBER_STANDS;id++) {
			out_anpp_stand[id] = new Table[nst];
			out_lai_stand[id] = new Table[nst];
			out_cmass_stand[id] = new Table[nst];
			out_diam_stand[id] = new Table[nst];
			out_height_stand[id] = new Table[nst];
			out_dens_stand[id] = new Table[nst];
			out_cmass_wood_stand[id] = new Table[nst];
			out_cmass_wood_harv_stand[id] = new Table[nst];
			out_cmass_mort_stand[id] = new Table[nst];
			out_agestruct_stand[id] = new Table[nst];
			out_diamstruct_stand[id] = new Table[nst];
			out_diamstruct_cmass_stand[id] = new Table[nst];
			out_storm_stand[id] = new Table[nst];
		}
	}

	if(date.year < nyear_spinup)
		return;

	bool open[NLANDCOVERTYPES];
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();

	for(int i=0;i<NLANDCOVERTYPES;i++)
		open[i] = false;

	Gridcell::iterator gc_itr = gridcell.begin();

	// Loop through Stands
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;

		if(stand.first_year == date.year || stand.clone_year == date.year) {
			if(stand.landcover == NATURAL) {
				open[NATURAL] = true;
			}
			else if(stand.landcover == FOREST) {
				open[FOREST] = true;
			}
		}

		++gc_itr;
	}

	if(PRINTFIRSTSTANDAFTERSPINUP && date.year == nyear_spinup) {
		open[NATURAL] = true;
		open[FOREST] = true;
	}

	if(open[NATURAL] || open[FOREST]) {

		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {

			Stand& stand = *gc_itr;
			StandType& st = stlist[stand.stid];

			int id = stand.id;
			char outfilename[100] = {'\0'}, buffer[50] = {'\0'}, format_string[50] = {'\0'};

			sprintf(format_string, "%%.%df", coordinates_precision);
			sprintf(buffer, "_");
			sprintf(buffer + strlen(buffer), format_string, lon);
			strcat(buffer, "_");
			sprintf(buffer + strlen(buffer), format_string, lat);
			sprintf(buffer + strlen(buffer), "_%d", id);
			strcat(buffer, ".out");

			// create a vector with the pft names
			std::vector<std::string> pfts;

			pftlist.firstobj();
			while (pftlist.isobj) {

				 Pft& pft=pftlist.getobj();	 
				 Standpft& standpft=stand.pft[pft.id];

				 pfts.push_back((char*)pft.name);

				 pftlist.nextobj();
			}

			ColumnDescriptors anpp_columns;
			anpp_columns += ColumnDescriptors(pfts,               8, 3);
			anpp_columns += ColumnDescriptor("Total",             8, 3);

			ColumnDescriptors diam_columns;
			diam_columns += ColumnDescriptors(pfts,               8, 3);

			ColumnDescriptors dens_columns;
			dens_columns += ColumnDescriptors(pfts,               8, 4);
			
			ColumnDescriptors storm_columns;
			storm_columns += ColumnDescriptor("SI",                9, 3);
			storm_columns += ColumnDescriptor("DamFrac",           9, 3);
			storm_columns += ColumnDescriptor("DWoodC_P", 9, 3);
			storm_columns += ColumnDescriptor("DamWoodC", 9, 3);

			ColumnDescriptors agestruct_columns;
			agestruct_columns += ColumnDescriptors(get_structure_string("age"),               9, 2);

			ColumnDescriptors diamstruct_columns;
			diamstruct_columns += ColumnDescriptors(get_structure_string("diam"),               9, 2);

			ColumnDescriptors diamstruct_cmass_columns;
			diamstruct_cmass_columns += ColumnDescriptors(get_structure_string("diam"),         9, 3);

			if(open[stand.landcover]) {

				if(print_anpp_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "anpp_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_anpp_stand[id][stand.stid].invalid())
						create_output_table(out_anpp_stand[id][stand.stid], outfilename, anpp_columns);
				}
				if(print_lai_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "lai_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_lai_stand[id][stand.stid].invalid())
						create_output_table(out_lai_stand[id][stand.stid], outfilename, anpp_columns);
				}
				if(print_cmass_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "cmass_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_cmass_stand[id][stand.stid].invalid())
						create_output_table(out_cmass_stand[id][stand.stid], outfilename, anpp_columns);
				}
				if(print_cmass_wood_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "cmass_wood_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_cmass_wood_stand[id][stand.stid].invalid())
						create_output_table(out_cmass_wood_stand[id][stand.stid], outfilename, anpp_columns);
				}
				if(print_cmass_wood_harv_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "cmass_wood_harv_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_cmass_wood_harv_stand[id][stand.stid].invalid())
						create_output_table(out_cmass_wood_harv_stand[id][stand.stid], outfilename, anpp_columns);
				}		
				if(print_cmass_mort_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "cmass_mort_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_cmass_mort_stand[id][stand.stid].invalid())
						create_output_table(out_cmass_mort_stand[id][stand.stid], outfilename, anpp_columns);
				}
				if(print_diam_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "diam_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_diam_stand[id][stand.stid].invalid())
						create_output_table(out_diam_stand[id][stand.stid], outfilename, diam_columns);
				}
				if(print_height_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "height_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_height_stand[id][stand.stid].invalid())
						create_output_table(out_height_stand[id][stand.stid], outfilename, diam_columns);
				}
				if(print_dens_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "dens_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_dens_stand[id][stand.stid].invalid())
						create_output_table(out_dens_stand[id][stand.stid], outfilename, dens_columns);
				}
				if(print_agestruct_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "agestruct_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_agestruct_stand[id][stand.stid].invalid())
						create_output_table(out_agestruct_stand[id][stand.stid], outfilename, agestruct_columns);
				}
				if(print_diamstruct_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "diamstruct_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_diamstruct_stand[id][stand.stid].invalid())
						create_output_table(out_diamstruct_stand[id][stand.stid], outfilename, diamstruct_columns);
				}
				if(print_diamstruct_cmass_stand) {
					strcpy(outfilename, dirname);
					strcat(outfilename, "diamstruct_cmass_wood_potharv_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_diamstruct_cmass_stand[id][stand.stid].invalid())
						create_output_table(out_diamstruct_cmass_stand[id][stand.stid], outfilename, diamstruct_cmass_columns);
				}
				if(print_storm_stand) {
					outfilename[0] = '\0';
					strcpy(outfilename, dirname);
					strcat(outfilename, "storm_");
					strcat(outfilename, (char*)st.name);
					strcat(outfilename, buffer);

					if(out_storm_stand[id][stand.stid].invalid())
						create_output_table(out_storm_stand[id][stand.stid], outfilename, storm_columns);
				}
			}

			++gc_itr;
		}
	}
}

void MiscOutput::closelocalfiles(Gridcell& gridcell) {

	if(!printseparatestands)
		return;
// This was removed in the crop_ncep branch (without any comments), not sure if it is needed.
	for(int id=0;id<MAXNUMBER_STANDS;id++) {

		for(int st=0;st<nst;st++) {
			if(!out_anpp_stand[id][st].invalid())
				close_output_table(out_anpp_stand[id][st]);
			if(!out_lai_stand[id][st].invalid())
				close_output_table(out_lai_stand[id][st]);
			if(!out_cmass_stand[id][st].invalid())
				close_output_table(out_cmass_stand[id][st]);
			if(!out_diam_stand[id][st].invalid())
				close_output_table(out_diam_stand[id][st]);
			if(!out_height_stand[id][st].invalid())
				close_output_table(out_height_stand[id][st]);
			if(!out_dens_stand[id][st].invalid())
				close_output_table(out_dens_stand[id][st]);
			if(!out_cmass_wood_stand[id][st].invalid())
				close_output_table(out_cmass_wood_stand[id][st]);
			if(!out_cmass_wood_harv_stand[id][st].invalid())
				close_output_table(out_cmass_wood_harv_stand[id][st]);
			if(!out_agestruct_stand[id][st].invalid())
				close_output_table(out_agestruct_stand[id][st]);
			if(!out_cmass_mort_stand[id][st].invalid())
				close_output_table(out_cmass_mort_stand[id][st]);
			if(!out_diamstruct_stand[id][st].invalid())
				close_output_table(out_diamstruct_stand[id][st]);
			if(!out_diamstruct_cmass_stand[id][st].invalid())
				close_output_table(out_diamstruct_cmass_stand[id][st]);
			if(!out_storm_stand[id][st].invalid())
				close_output_table(out_storm_stand[id][st]);
		}
	}

	for(int id=0;id<MAXNUMBER_STANDS;id++) {
		if(out_anpp_stand[id])
			delete[] out_anpp_stand[id];
		if(out_lai_stand[id])
			delete[] out_lai_stand[id];
		if(out_cmass_stand[id])
			delete[] out_cmass_stand[id];
		if(out_diam_stand[id])
			delete[] out_diam_stand[id];
		if(out_height_stand[id])
			delete[] out_height_stand[id];
		if(out_dens_stand[id])
			delete[] out_dens_stand[id];
		if(out_cmass_wood_stand[id])
			delete[] out_cmass_wood_stand[id];
		if(out_cmass_wood_harv_stand[id])
			delete[] out_cmass_wood_harv_stand[id];
		if(out_cmass_mort_stand[id])
			delete[] out_cmass_mort_stand[id];
		if(out_agestruct_stand[id])
			delete[] out_agestruct_stand[id];
		if(out_diamstruct_stand[id])
			delete[] out_diamstruct_stand[id];
		if(out_diamstruct_cmass_stand[id])
			delete[] out_diamstruct_cmass_stand[id];
		if (out_storm_stand[id])
			delete[] out_storm_stand[id];

		out_anpp_stand[id] = NULL;
		out_lai_stand[id] = NULL;
		out_cmass_stand[id] = NULL;
		out_diam_stand[id] = NULL;
		out_height_stand[id] = NULL;
		out_dens_stand[id] = NULL;
		out_cmass_wood_stand[id] = NULL;
		out_cmass_wood_harv_stand[id] = NULL;
		out_cmass_mort_stand[id] = NULL;
		out_agestruct_stand[id] = NULL;
		out_diamstruct_stand[id] = NULL;
		out_diamstruct_cmass_stand[id] = NULL;
		out_storm_stand[id] = NULL;
	}
}

} // namespace
