///////////////////////////////////////////////////////////////////////////////////////
/// \file outputmodule.cpp
/// \brief Implementation of the common output module
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "commonoutput.h"
#include "parameters.h"
#include "guess.h"

// Years between output time points in plot function output (Windows shell only)
const int PLOT_INTERVAL = 1;

// Years between updates of soil water and stand structure plot (Windows shell only)
const int PLOT_UPDATE_INTERVAL = 20;

// Years between updates of 3D vegetation view (Windows shell only)
const int VEG3D_UPDATE_INTERVAL = 5;

namespace GuessOutput {

REGISTER_OUTPUT_MODULE("common", CommonOutput)

CommonOutput::CommonOutput() {
	// Annual output variables
	declare_parameter("file_cmass", &file_cmass, 300, "C biomass output file");
	declare_parameter("file_cmass_agb", &file_cmass_agb, 300, "C above ground biomass output file");
	declare_parameter("file_anpp", &file_anpp, 300, "Annual NPP output file");
	declare_parameter("file_agpp", &file_agpp, 300, "Annual GPP with all stress output file");
	declare_parameter("file_agpp_no", &file_agpp_no, 300, "Annual GPP with no stresses output file");
	declare_parameter("file_agpp_ns", &file_agpp_ns, 300, "Annual GPP with N stress output file");
	declare_parameter("file_fpc", &file_fpc, 300, "FPC output file");
	declare_parameter("file_aaet", &file_aaet, 300, "Annual AET output file");
	declare_parameter("file_lai", &file_lai, 300, "LAI output file");
	declare_parameter("file_cflux", &file_cflux, 300, "C fluxes output file");
	declare_parameter("file_doc", &file_doc, 300, "DOC output file");
	declare_parameter("file_dens", &file_dens, 300, "Tree density output file");
	declare_parameter("file_cpool", &file_cpool, 300, "C output file");
	declare_parameter("file_clitter", &file_clitter, 300, "Litter C output file");
	declare_parameter("file_soildepth", &file_soildepth, 300, "Soil layer depths output file");
	declare_parameter("file_csoil", &file_csoil, 300, "Soil layer C output file");
	declare_parameter("file_csoil_depth", &file_csoil_depth, 300, "Soil depth C output file");
	declare_parameter("file_nsoil_depth", &file_nsoil_depth, 300, "Soil depth N output file");
	declare_parameter("file_csoil_cm", &file_csoil_cm, 300, "Soil depth with 10cm interval C output file");
	declare_parameter("file_csoilp_cm", &file_csoilp_cm, 300, "Soil depth with 10cm interval permafrost C output file");
	declare_parameter("file_netnminsoil_depth", &file_netnminsoil_depth, 300, "Soil depth net N mineralisation output file");
	declare_parameter("file_rh", &file_rh, 300, "Soil layer C respiration output file");
	declare_parameter("file_rh_depth", &file_rh_depth, 300, "Soil depth C respiration output file");
	declare_parameter("file_tau_soil", &file_tau_soil, 300, "Soil layer C residence time output file");
	declare_parameter("file_tau_soil_depth", &file_tau_soil_depth, 300, "Soil depth C residence time output file");
	declare_parameter("file_volume_org_depth", &file_volume_org_depth, 300, "Soil depths volume of organic matter output file");
	declare_parameter("file_volume_min_depth", &file_volume_min_depth, 300, "Soil depths volume of mineral matter output file");
	declare_parameter("file_volume_root_depth", &file_volume_root_depth, 300, "Soil depths volume of root matter output file");
	declare_parameter("file_temp_mod", &file_temp_mod, 300, "Soil layer temperature rate modifier output file");
	declare_parameter("file_moist_mod", &file_moist_mod, 300, "Soil layer moisture rate modifier output file");
	declare_parameter("file_tot_mod", &file_tot_mod, 300, "Soil layer temperature and moisture rate modifier output file");
	declare_parameter("file_temp_mod_depth", &file_temp_mod_depth, 300, "Soil layer temperature rate modifier output file");
	declare_parameter("file_moist_mod_depth", &file_moist_mod_depth, 300, "Soil layer moisture rate modifier output file");
	declare_parameter("file_tot_mod_depth", &file_tot_mod_depth, 300, "Soil layer temperature and moisture rate modifier output file");
	declare_parameter("file_csom", &file_csom, 300, "Soil Organic Matter C output file");
	declare_parameter("file_ctonsom", &file_ctonsom, 300, "Soil Organic Matter C:N ratio output file");
	declare_parameter("file_soildeclim", &file_soildeclim, 300, "Annual nitrogen decay rate limitation on SOM decomposition output file");
	declare_parameter("file_soilgrolim", &file_soilgrolim, 300, "Annual nitrogen growth rate limitation on SOM decomposition output file");
	declare_parameter("file_soiltotlim", &file_soiltotlim, 300, "Annual nitrogen total limitation on SOM decomposition output file");
	declare_parameter("file_soilcuelim", &file_soilcuelim, 300, "Annual carbon use efficiency on SOM decomposition output file");
	declare_parameter("file_clitterinput", &file_clitterinput, 300, "Annual carbon litter input to soil column output file");
	declare_parameter("file_clitterinput_depth", &file_clitterinput_depth, 300, "Annual carbon litter input to soil column output file");
	declare_parameter("file_csomlitterinput", &file_csomlitterinput, 300, "Annual carbon litter input to SOM pools output file");
	declare_parameter("file_soilcperma", &file_soilcperma, 300, "Fraction of soil C in permafrost output file");
	declare_parameter("file_snowdays", &file_snowdays, 300, "Snow days output file");
	declare_parameter("file_mtemp", &file_mtemp, 300, "Air, snow and soil layer temperature output file");
	declare_parameter("file_mwmasslyr", &file_mwmasslyr, 300, "Monthly soil layer water output file");
	declare_parameter("file_mimasslyr", &file_mimasslyr, 300, "Monthly soil layer ice output file");
	declare_parameter("file_mwimasslyr", &file_mwimasslyr, 300, "Monthly soil layer water+ice output file");
	declare_parameter("file_wsatslyr", &file_wsatslyr, 300, "Yearly water holding capacity at saturation output file");
	declare_parameter("file_mdrainlyr", &file_mdrainlyr, 300, "Monthly soil layer lateral drainage output file");
	declare_parameter("file_mtemp_depth", &file_mtemp_depth, 300, "Air, snow and soil layer temperature output file");
	declare_parameter("file_mwmass_depth", &file_mwmass_depth, 300, "Monthly soil layer water output file");
	declare_parameter("file_mimass_depth", &file_mimass_depth, 300, "Monthly soil layer ice output file");
	declare_parameter("file_mwimass_depth", &file_mwimass_depth, 300, "Monthly soil layer water+ice output file");
	declare_parameter("file_mdrain_depth", &file_mdrain_depth, 300, "Monthly soil layer lateral drainage output file");
	declare_parameter("file_runoff", &file_runoff, 300, "Runoff output file");
	declare_parameter("file_wetland_water_added", &file_wetland_water_added, 300, "Wetland water added output file");

	declare_parameter("file_firert", &file_firert, 300, "Fire retrun time output file");

	declare_parameter("file_nmass", &file_nmass, 300, "N biomass output file");
	declare_parameter("file_cton_leaf", &file_cton_leaf, 300, "Mean leaf C:N output file");
	declare_parameter("file_nsources", &file_nsources, 300, "Annual nitrogen sources output file");
	declare_parameter("file_npool", &file_npool, 300, "Nitrogen pools output file");
	declare_parameter("file_nlitter", &file_nlitter, 300, "Litter nitrogen output file");
	declare_parameter("file_nuptake", &file_nuptake, 300, "Annual nitrogen uptake output file");
	declare_parameter("file_vmaxnlim", &file_vmaxnlim, 300, "Annual nitrogen limitation on vm output file");
	declare_parameter("file_lambda", &file_lambda, 300, "Annual lambda limitation on photosynthesis output file");
	declare_parameter("file_nflux", &file_nflux, 300, "Annual nitrogen fluxes output file");
	declare_parameter("file_ngases", &file_ngases, 300, "Annual nitrogen gases output file");
	declare_parameter("file_soil_npool", &file_soil_npool, 300, "Annual soil N pools output file");
	declare_parameter("file_soil_nflux", &file_soil_nflux, 300, "Annual soil N fluxes output file");

	declare_parameter("file_speciesheights", &file_speciesheights, 300, "Mean species heights");

	declare_parameter("file_wflux", &file_wflux, 300, "Annual water fluxes output file");
	declare_parameter("file_wpool", &file_wpool, 300, "Annual water pools output file");
	declare_parameter("file_wsats", &file_wsats, 300, "Annual water holding capacity at saturation output file");
	declare_parameter("file_wsats_depth", &file_wsats_depth, 300, "Annual water holding capacity at saturation output file");

	// Monthly output variables
	declare_parameter("file_mnpp", &file_mnpp, 300, "Monthly NPP output file");
	declare_parameter("file_mlai", &file_mlai, 300, "Monthly LAI output file");
	declare_parameter("file_mgpp", &file_mgpp, 300, "Monthly GPP-LeafResp output file");
	declare_parameter("file_mra", &file_mra, 300, "Monthly autotrophic respiration output file");
	declare_parameter("file_maet", &file_maet, 300, "Monthly AET output file");
	declare_parameter("file_mpet", &file_mpet, 300, "Monthly PET output file");
	declare_parameter("file_mevap", &file_mevap, 300, "Monthly Evap output file");
	declare_parameter("file_mrunoff", &file_mrunoff, 300, "Monthly runoff output file");
	declare_parameter("file_mintercep", &file_mintercep, 300, "Monthly intercep output file");
	declare_parameter("file_mrh", &file_mrh, 300, "Monthly heterotrophic respiration output file");
	declare_parameter("file_mnee", &file_mnee, 300, "Monthly NEE output file");
	declare_parameter("file_mnetmin", &file_mnetmin, 300, "Monthly net mineralisation output file");
	declare_parameter("file_mcton_leaf", &file_mcton_leaf, 300, "Monthly leaf C:N ratio output file");
	declare_parameter("file_mvmaxnlim", &file_mvmaxnlim, 300, "Monthly nitrogen limitation on vm output file");
	declare_parameter("file_mlambda", &file_mlambda, 300, "Monthly lambda limitation on photosynthesis output file");
	declare_parameter("file_mwmass", &file_mwmass, 300, "Monthly liquid water output file");
	declare_parameter("file_mimass", &file_mimass, 300, "Monthly frozen water output file");
	declare_parameter("file_mwcont_upper", &file_mwcont_upper, 300, "Monthly wcont_upper output file");
	declare_parameter("file_mwcont_lower", &file_mwcont_lower, 300, "Monthly wcont_lower output file");
	// bvoc
	declare_parameter("file_aiso", &file_aiso, 300, "annual isoprene flux output file");
	declare_parameter("file_miso", &file_miso, 300, "monthly isoprene flux output file");
	declare_parameter("file_amon", &file_amon, 300, "annual monoterpene flux output file");
	declare_parameter("file_mmon", &file_mmon, 300, "monthly monoterpene flux output file");
	declare_parameter("file_amon_mt1", &file_amon_mt1, 300, "annual endocyclic monoterpene flux output file");	
	declare_parameter("file_amon_mt2", &file_amon_mt2, 300, "annual other monoterpene flux output file");
	declare_parameter("file_mmon_mt1", &file_mmon_mt1, 300, "monthly endocyclic monoterpene flux output file");	
	declare_parameter("file_mmon_mt2", &file_mmon_mt2, 300, "monthly other monoterpene flux output file");

	if ( firemodel == BLAZE ) {
		declare_parameter("file_aburned_area_out", &file_aburned_area_out, 300, "BLAZE burned area output file");
		declare_parameter("file_mburned_area_out", &file_mburned_area_out, 300, "BLAZE monthly burned area output file");
		declare_parameter("file_simfireanalysis_out", &file_simfireanalysis_out, 300, "SIMFIRE analytics output");
	}

	declare_parameter("file_mch4", &file_mch4, 300, "Monthly CH4 emissions, total");
	declare_parameter("file_mch4diff", &file_mch4diff, 300, "Monthly CH4 emissions, diffusion");
	declare_parameter("file_mch4plan", &file_mch4plan, 300, "Monthly CH4 emissions, plant-mediated");
	declare_parameter("file_mch4ebull", &file_mch4ebull, 300, "Monthly CH4 emissions, ebullition");
	declare_parameter("file_msnow", &file_msnow, 300, "Monthly snow depth");
	declare_parameter("file_msnowcover", &file_msnowcover, 300, "Monthly snow cover");
	declare_parameter("file_mwmass_stand", &file_mwmass_stand, 300, "Monthly standing water");
	declare_parameter("file_mwtp", &file_mwtp, 300, "Monthly water table depth");
	declare_parameter("file_mald", &file_mald, 300, "Monthly active layer depth");
	declare_parameter("file_msoiltempdepth25", &file_msoiltempdepth25, 300, "Soil temperature output file (25cm depth)");
}


CommonOutput::~CommonOutput() {
}

/// Define all output tables and their formats
void CommonOutput::init() {

	define_output_tables();
}

/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void CommonOutput::define_output_tables() {

	//Extra number of decimals when output for benchmarks
#ifdef RUN_BENCHMARKS	
	const int bm_extra_prec = 2;
#else
	const int bm_extra_prec = 0;
#endif

	// create a vector with the pft names
	std::vector<std::string> pfts;

	// create a vector with the crop pft names
	std::vector<std::string> crop_pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();

		pfts.push_back((char*)pft.name);

		if (pft.landcover == CROPLAND)
			crop_pfts.push_back((char*)pft.name);

		pftlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if (run_landcover) {
		const char* landcover_string[] = { "Urban_sum", "Crop_sum", "Pasture_sum",
			"Forest_sum", "Natural_sum", "Peatland_sum", "Barren_sum" };
		for (int i = 0; i < NLANDCOVERTYPES; i++) {
			if (run[i]) {
				landcovers.push_back(landcover_string[i]);
			}
		}
	}

	// Create the month columns
	ColumnDescriptors month_columns, month_columns_wide, month_columns_large;
	xtring months[] = { "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec" };
	for (int i = 0; i < 12; i++) {
		month_columns += ColumnDescriptor(months[i], 8, 3);
		month_columns_wide += ColumnDescriptor(months[i], 10, 3);
		month_columns_large += ColumnDescriptor(months[i], 7, 0);
	}

	// Create soil depth columns
	ColumnDescriptors soildepths_columns;
	soildepths_columns += ColumnDescriptor("0.3m", 9, 3);
	soildepths_columns += ColumnDescriptor("0.5m", 9, 3);
	soildepths_columns += ColumnDescriptor("1.0m", 9, 3);
	soildepths_columns += ColumnDescriptor("1.5m", 9, 3);
	soildepths_columns += ColumnDescriptor("2.0m", 9, 3);
	soildepths_columns += ColumnDescriptor("Total", 9, 3);

	// Create soil layer columns
	ColumnDescriptors soillyr_columns_xlarge, soillyr_columns_large, soillyr_columns_short, soillyr_columns_med, soillyr_columns_high, soillyr_columns, surfsoillyr_columns;
	surfsoillyr_columns += ColumnDescriptor("surf", 11, 6);
	char lyr_nr[20];
	int number_of_layers = hydrology == GERTEN ? NSOILLAYERS_SIMPLE : NSOILLAYERS_MAX;
	for (int lyr = 0; lyr < number_of_layers; lyr++) {
		sprintf(lyr_nr, "lyr%d", lyr);
		surfsoillyr_columns += ColumnDescriptor(lyr_nr, 10, 5);
		soillyr_columns_xlarge += ColumnDescriptor(lyr_nr, 12, 1);
		soillyr_columns_large += ColumnDescriptor(lyr_nr, 7, 1);
		soillyr_columns_short += ColumnDescriptor(lyr_nr, 7, 3);
		soillyr_columns_med += ColumnDescriptor(lyr_nr, 10, 5);
		soillyr_columns_high += ColumnDescriptor(lyr_nr, 13, 10);
		soillyr_columns += ColumnDescriptor(lyr_nr, 7, 4);
	}
	surfsoillyr_columns += ColumnDescriptor("Total", 10, 5);
	soillyr_columns_xlarge += ColumnDescriptor("Total", 10, 1);
	soillyr_columns_large += ColumnDescriptor("Total", 7, 1);
	soillyr_columns_short += ColumnDescriptor("Total", 7, 3);
	soillyr_columns_med += ColumnDescriptor("Total", 12, 7);
	soillyr_columns_high += ColumnDescriptor("Total", 13, 10);

	// Create soil layer columns
	ColumnDescriptors temp_columns;
	temp_columns += ColumnDescriptor("air", 8, 2);
	if (hydrology == GERTEN) {
		for (int lyr = 0; lyr < NSOILLAYERS_SIMPLE; lyr++) {
			sprintf(lyr_nr, "soil%d", lyr);
			temp_columns += ColumnDescriptor(lyr_nr, 8, 2);
		}
	}
	else {
		for (int lyr = NSNOWLAYERS_MAX - 1; lyr >= 0; lyr--) {
			sprintf(lyr_nr, "snow%d", lyr);
			temp_columns += ColumnDescriptor(lyr_nr, 8, 2);
		}
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			sprintf(lyr_nr, "soil%d", lyr);
			temp_columns += ColumnDescriptor(lyr_nr, 8, 2);
		}
		for (int lyr = 0; lyr < NPADLAYERS; lyr++) {
			sprintf(lyr_nr, "pad%d", lyr);
			temp_columns += ColumnDescriptor(lyr_nr, 8, 2);
		}
	}

	// Create soil layer columns
	ColumnDescriptors soms_columns, soms_ratio_columns;
	xtring soms[] = { "SURFMETA", "SURFSTRUCT", "SURFFWD", "SURFCWD", "SURFMICRO", "SURFHUMUS", "SOILMETA", "SOILSTRUCT", "SOILFWD", "SOILCWD", "SOILMICRO", "SLOWSOM", "PASSIVESOM" };
	for (int p = 0; p < NSOMPOOL; p++) {
		soms_columns += ColumnDescriptor(soms[p], 11, 6);
		soms_ratio_columns += ColumnDescriptor(soms[p], 11, 2);
	}
	soms_columns += ColumnDescriptor("Total", 11, 6);
	soms_ratio_columns += ColumnDescriptor("Total", 11, 2);

	// Create soil depth columns
	ColumnDescriptors surfsoildepth_columns, surfsoildepth_columns_med, surfsoildepth_columns_low,
		soildepth_columns, soildepth_columns_med, soildepth_columns_low;
	xtring depths[] = { "d0_0.3", "d0.3_0.5", "d0.5_1", "d1_1.5", "d1.5_2", "d2_3", "d3_max", "Total" };
	surfsoildepth_columns += ColumnDescriptor("Surf", 9, 3);
	surfsoildepth_columns_med += ColumnDescriptor("Surf", 11, 6);
	surfsoildepth_columns_low += ColumnDescriptor("Surf", 9, 1);
	for (int d = 0; d < 8; d++) {
		soildepth_columns += ColumnDescriptor(depths[d], 9, 3);
		soildepth_columns_med += ColumnDescriptor(depths[d], 11, 6);
		soildepth_columns_low += ColumnDescriptor(depths[d], 9, 1);
		surfsoildepth_columns += ColumnDescriptor(depths[d], 9, 3);
		surfsoildepth_columns_med += ColumnDescriptor(depths[d], 11, 6);
		surfsoildepth_columns_low += ColumnDescriptor(depths[d], 9, 1);
	}

	// Create soil temperature depth columns
	ColumnDescriptors tempdepth_columns;
	xtring temp_depths[] = { "d0_0.3", "d0.3_0.5", "d0.5_1", "d1_1.5", "d1.5_2", "d2_3", "d3_5", "d5_10", "d10_20", "d20_30", "d30_40", "d40_50" };
	tempdepth_columns += ColumnDescriptor("AirT", 9, 3);
	for (int d = 0; d < 12; d++) {
		tempdepth_columns += ColumnDescriptor(temp_depths[d], 9, 3);
	}

	// Create 10cm soil layer down to 3m columns
	ColumnDescriptors depth_cm_columns;
	double depth_cm = 0;
	char depth_in_m[20];
	for (int d = 0; d < 30; d++) {
		depth_cm += 0.1;
			sprintf(depth_in_m, "d%g", depth_cm);
			depth_cm_columns += ColumnDescriptor(depth_in_m, 8, 2);
	}

	// Create the columns for each output file

	// PFT
	ColumnDescriptors pft_columns;
	pft_columns += ColumnDescriptors(pfts,                 8, 2);

	// CMASS
	ColumnDescriptors cmass_columns;
	cmass_columns += ColumnDescriptors(pfts,               8, 3);
	cmass_columns += ColumnDescriptor("Total",             8, 3);
	cmass_columns += ColumnDescriptors(landcovers,        13, 3);

	// ALD
	ColumnDescriptors mald_columns = month_columns;
	mald_columns += ColumnDescriptor("MAXALD", 8, 3);
	mald_columns += ColumnDescriptor("PF_full", 8, 3);
	mald_columns += ColumnDescriptor("PF_3m", 8, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;

	// AGPP
	ColumnDescriptors agpp_columns = cmass_columns;

	// FPC
	ColumnDescriptors fpc_columns = cmass_columns;

	// AET
	ColumnDescriptors aaet_columns;
	aaet_columns += ColumnDescriptors(pfts,                8, 2);
	aaet_columns += ColumnDescriptor("Total",              8, 2);
	aaet_columns += ColumnDescriptors(landcovers,         13, 2);

	// DENS
	ColumnDescriptors dens_columns;
	dens_columns += ColumnDescriptors(pfts,                8, 4);
	dens_columns += ColumnDescriptor("Total",              8, 4);
	dens_columns += ColumnDescriptors(landcovers,         13, 4);

	// LAI
	ColumnDescriptors lai_columns = dens_columns;

	// CFLUX
	ColumnDescriptors cflux_columns;
	cflux_columns += ColumnDescriptor("Veg",               8, 3);
	cflux_columns += ColumnDescriptor("Repr",              8, 3);
	cflux_columns += ColumnDescriptor("Soil",              8, 3);
	cflux_columns += ColumnDescriptor("Fire",             10, 5);
	cflux_columns += ColumnDescriptor("Est",               8, 3);
	if (run_landcover) {
		 cflux_columns += ColumnDescriptor("Seed",         8, 3);
		 cflux_columns += ColumnDescriptor("Harvest",      9, 3);
		 cflux_columns += ColumnDescriptor("LU_ch",        9, 3);
		 cflux_columns += ColumnDescriptor("Slow_h",       9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",              10 + bm_extra_prec, 5 + bm_extra_prec);

	ColumnDescriptors doc_columns;
	doc_columns += ColumnDescriptor("Total",              10, 3);
	doc_columns += ColumnDescriptors(landcovers,          13, 3);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",              8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",         8, 3);
		cpool_columns += ColumnDescriptor("SoilfC",        8, 3);
		cpool_columns += ColumnDescriptor("SoilsC",        8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LitterC",       9, 3);
		cpool_columns += ColumnDescriptor("SoilC",         9, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC",   10, 3);
	}
	cpool_columns += ColumnDescriptor("Total",            10 + bm_extra_prec, 3 + bm_extra_prec);

	// CLITTER
	ColumnDescriptors clitter_columns = cmass_columns;

	// Permafrost fraction
	ColumnDescriptors soilc_perma_columns;
	soilc_perma_columns += ColumnDescriptor("Active", 8, 3);
	soilc_perma_columns += ColumnDescriptor("Frost", 8, 3);

	// Permafrost fraction
	ColumnDescriptors snow_columns;
	snow_columns += ColumnDescriptor("Snowdays", 9, 1);

	// FIRERT
	ColumnDescriptors firert_columns;
	firert_columns += ColumnDescriptor("FireRT",			8, 1);
	firert_columns += ColumnDescriptor("BurntFr",			8, 5);

	// BLAZE burnt area 
	ColumnDescriptors blaze_columns;
	blaze_columns += ColumnDescriptor("BurntFr",			9, 5);

	// SIMFIRE Analysis 
	ColumnDescriptors simfireanalysis_columns;
	simfireanalysis_columns += ColumnDescriptor("Biome",	6, 0);
	simfireanalysis_columns += ColumnDescriptor("MxNest",	7, 0);
	simfireanalysis_columns += ColumnDescriptor("PopDens",	10, 3);
	simfireanalysis_columns += ColumnDescriptor("AMxFApar",	10, 7);
	simfireanalysis_columns += ColumnDescriptor("FireProb",	10, 7);
	simfireanalysis_columns += ColumnDescriptor("Region",	7, 0);
	
	// RUNOFF
	ColumnDescriptors runoff_columns;
	runoff_columns += ColumnDescriptor("Surf",             8, 1);
	runoff_columns += ColumnDescriptor("Drain",            8, 1);
	runoff_columns += ColumnDescriptor("Base",             8, 1);
	runoff_columns += ColumnDescriptor("Total",            9, 1);

	// WETLAND WATER ADDED
	ColumnDescriptors wetland_water_added_columns;
	wetland_water_added_columns += ColumnDescriptor("H2OAdded", 10, 1);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts,      8, 2);

	// AISO
	ColumnDescriptors aiso_columns;
	aiso_columns += ColumnDescriptors(pfts,               10, 3);
	aiso_columns += ColumnDescriptor("Total",             10, 3);
	aiso_columns += ColumnDescriptors(landcovers,         13, 3);

	// MONOTERPENES
	ColumnDescriptors amon_columns = aiso_columns;

	// CTON
	ColumnDescriptors cton_columns;
	cton_columns += ColumnDescriptors(pfts,                8, 1);
	cton_columns += ColumnDescriptor("Total",              8, 1);
	cton_columns += ColumnDescriptors(landcovers,         13, 1);

	// NSOURCES
	ColumnDescriptors nsources_columns;
	nsources_columns += ColumnDescriptor("NH4dep",         8, 2);
	nsources_columns += ColumnDescriptor("NO3dep",         8, 2);
	nsources_columns += ColumnDescriptor("fix",            8, 2);
	nsources_columns += ColumnDescriptor("fert",           8, 2);
	nsources_columns += ColumnDescriptor("input",          8, 2);
	nsources_columns += ColumnDescriptor("min",            8, 2);
	nsources_columns += ColumnDescriptor("imm",            8, 2);
	nsources_columns += ColumnDescriptor("netmin",         8, 2);
	nsources_columns += ColumnDescriptor("Total",          8, 2);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",              9, 4);
	npool_columns += ColumnDescriptor("LitterN",           9, 4);
	npool_columns += ColumnDescriptor("SoilN",             9, 4);
	npool_columns += ColumnDescriptor("AvailN",           11, 6);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN",    10, 4);
	}

	npool_columns += ColumnDescriptor("Total",            10 + bm_extra_prec, 4 + bm_extra_prec);

	// NMASS
	ColumnDescriptors nmass_columns;
	nmass_columns += ColumnDescriptors(pfts,               8, 2);
	nmass_columns += ColumnDescriptor("Total",             8, 2);
	nmass_columns += ColumnDescriptors(landcovers,        13, 2);

	// NUPTAKE
	ColumnDescriptors nuptake_columns = nmass_columns;

	// NLITTER
	ColumnDescriptors nlitter_columns = nmass_columns;

	// VMAXNLIM
	ColumnDescriptors vmaxnlim_columns;
	vmaxnlim_columns += ColumnDescriptors(pfts,            8, 2);
	vmaxnlim_columns += ColumnDescriptor("Total",          8, 2);
	vmaxnlim_columns += ColumnDescriptors(landcovers,     13, 2);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("NH4dep",            8, 2);
	nflux_columns += ColumnDescriptor("NO3dep",            8, 2);
	nflux_columns += ColumnDescriptor("fix",               8, 2);
	nflux_columns += ColumnDescriptor("fert",              8, 2);
	nflux_columns += ColumnDescriptor("flux",              8, 2);
	nflux_columns += ColumnDescriptor("leach",             8, 2);
	if (run_landcover) {
		nflux_columns += ColumnDescriptor("seed",		   8, 2);
		nflux_columns += ColumnDescriptor("harvest",       8, 2);
		nflux_columns += ColumnDescriptor("LU_ch",         8, 3);
		nflux_columns += ColumnDescriptor("Slow_h",        8, 3);
	}
	nflux_columns += ColumnDescriptor("NEE",               8 + bm_extra_prec, 2 + bm_extra_prec);

	// NGASES
	ColumnDescriptors ngases_columns;
	ngases_columns += ColumnDescriptor("NH3_fire",         9, 4);
	ngases_columns += ColumnDescriptor("NH3_soil",         9, 4);
	ngases_columns += ColumnDescriptor("NOx_fire",          9, 4);
	ngases_columns += ColumnDescriptor("NOx_soil",          9, 4);
	ngases_columns += ColumnDescriptor("N2O_fire",         9, 4);
	ngases_columns += ColumnDescriptor("N2O_soil",         9, 4);
	ngases_columns += ColumnDescriptor("N2_fire",          9, 4);
	ngases_columns += ColumnDescriptor("N2_soil",          9, 4);
	ngases_columns += ColumnDescriptor("Total",            9, 4);

	// SOIL N TRANSFORMATION - pools
	ColumnDescriptors soil_npool_columns;
	soil_npool_columns += ColumnDescriptor("NH4", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO3", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO2", 11, 4);
	soil_npool_columns += ColumnDescriptor("NO",  11, 4);
	soil_npool_columns += ColumnDescriptor("N2O", 11, 4);
	soil_npool_columns += ColumnDescriptor("N2",  11, 4);


	// SOIL N TRANSFORMATION - fluxes
	ColumnDescriptors soil_nflux_columns;
	soil_nflux_columns += ColumnDescriptor("NH3",  12, 6);
	soil_nflux_columns += ColumnDescriptor("NO",   12, 6);
	soil_nflux_columns += ColumnDescriptor("N2O",  12, 6);
	soil_nflux_columns += ColumnDescriptor("N2",   12, 6);

	// WFLUX
	ColumnDescriptors wflux_columns;
	wflux_columns += ColumnDescriptor("Prec", 10, 1);
	if (run_landcover) {
		wflux_columns += ColumnDescriptor("Irri", 8, 1);
		wflux_columns += ColumnDescriptor("Added", 8, 1);
	}
	wflux_columns += ColumnDescriptor("Intcep", 8, 1);
	wflux_columns += ColumnDescriptor("AET", 8, 1);
	wflux_columns += ColumnDescriptor("Evap", 8, 1);
	wflux_columns += ColumnDescriptor("Rsurf", 8, 1);
	wflux_columns += ColumnDescriptor("Rdrain", 8, 1);
	wflux_columns += ColumnDescriptor("Rbase", 8, 1);
	wflux_columns += ColumnDescriptor("NEE", 8, 1);

	// WPOOL
	ColumnDescriptors wpool_columns;
	wpool_columns += ColumnDescriptor("SnowWater", 10, 1);
	wpool_columns += ColumnDescriptor("SnowIce", 8, 1);
	wpool_columns += ColumnDescriptor("SoilWater", 10, 1);
	wpool_columns += ColumnDescriptor("SoilIce", 9, 1);

	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_cmass,          file_cmass,          cmass_columns);
	create_output_table(out_cmass_agb,      file_cmass_agb,      cmass_columns); 
	create_output_table(out_anpp,           file_anpp,           anpp_columns);
	create_output_table(out_agpp,           file_agpp,           agpp_columns);
	create_output_table(out_agpp_no,        file_agpp_no,        agpp_columns);
	create_output_table(out_agpp_ns,        file_agpp_ns,        agpp_columns);
	create_output_table(out_fpc,            file_fpc,            fpc_columns);
	create_output_table(out_aaet,           file_aaet,           aaet_columns);
	create_output_table(out_dens,           file_dens,           dens_columns);
	create_output_table(out_lai,            file_lai,            lai_columns);
	create_output_table(out_cflux,          file_cflux,          cflux_columns);
	create_output_table(out_doc,	        file_doc,			 doc_columns);
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
	create_output_table(out_clitter,        file_clitter,        clitter_columns);
	create_output_table(out_soildepth,      file_soildepth,      soillyr_columns_large);
	create_output_table(out_csoil,          file_csoil,          soillyr_columns_med);
	create_output_table(out_csoil_depth,    file_csoil_depth,    surfsoildepth_columns);
	create_output_table(out_nsoil_depth,    file_nsoil_depth,    surfsoildepth_columns);
	create_output_table(out_csoil_cm,       file_csoil_cm,       depth_cm_columns);
	create_output_table(out_csoilp_cm,      file_csoilp_cm,      depth_cm_columns);
	create_output_table(out_netnminsoil_depth, file_netnminsoil_depth, soildepth_columns);
	create_output_table(out_rh,             file_rh,             soillyr_columns_med); 
	create_output_table(out_rh_depth,       file_rh_depth,       surfsoildepth_columns_med);
	create_output_table(out_tau_soil,		file_tau_soil,       soillyr_columns_xlarge);
	create_output_table(out_tau_soil_depth, file_tau_soil_depth, surfsoildepth_columns_low);
	create_output_table(out_temp_mod,       file_temp_mod,       soillyr_columns);
	create_output_table(out_moist_mod,      file_moist_mod,      soillyr_columns);
	create_output_table(out_tot_mod,        file_tot_mod,        soillyr_columns);
	create_output_table(out_temp_mod_depth, file_temp_mod_depth, soildepth_columns);
	create_output_table(out_moist_mod_depth,file_moist_mod_depth,soildepth_columns);
	create_output_table(out_tot_mod_depth,  file_tot_mod_depth,  soildepth_columns);
	create_output_table(out_volume_org_depth, file_volume_org_depth, soildepth_columns);
	create_output_table(out_volume_min_depth, file_volume_min_depth, soildepth_columns);
	create_output_table(out_volume_root_depth, file_volume_root_depth, soildepth_columns);
	create_output_table(out_csom,           file_csom,           soms_columns);
	create_output_table(out_ctonsom,        file_ctonsom,        soms_ratio_columns);
	create_output_table(out_soildeclim,     file_soildeclim,     soillyr_columns_short);
	create_output_table(out_soilgrolim,     file_soilgrolim,     soillyr_columns_short);
	create_output_table(out_soiltotlim,     file_soiltotlim,     soillyr_columns_short);
	create_output_table(out_soilcuelim,     file_soilcuelim,     soillyr_columns_short);
	create_output_table(out_clitterinput,   file_clitterinput,   surfsoillyr_columns);
	create_output_table(out_clitterinput_depth, file_clitterinput_depth, surfsoildepth_columns);
	create_output_table(out_csomlitterinput,file_csomlitterinput,soms_columns);
	create_output_table(out_soilcperma,     file_soilcperma,     soilc_perma_columns);
	create_output_table(out_mtemp,          file_mtemp,          temp_columns);
	create_output_table(out_mwmasslyr,      file_mwmasslyr,      soillyr_columns_short);
	create_output_table(out_mimasslyr,      file_mimasslyr,      soillyr_columns_short);
	create_output_table(out_mwimasslyr,     file_mwimasslyr,     soillyr_columns_short);
	create_output_table(out_wsatslyr,       file_wsatslyr,       soillyr_columns_large);
	create_output_table(out_mdrainlyr,      file_mdrainlyr,	     soillyr_columns_short);
	create_output_table(out_mtemp_depth,    file_mtemp_depth,    tempdepth_columns);
	create_output_table(out_mwmass_depth,   file_mwmass_depth,   soildepth_columns);
	create_output_table(out_mimass_depth,   file_mimass_depth,   soildepth_columns);
	create_output_table(out_mwimass_depth,  file_mwimass_depth,  soildepth_columns);
	create_output_table(out_mdrain_depth,   file_mdrain_depth,   soildepth_columns);
	create_output_table(out_snowdays,       file_snowdays,       snow_columns);

	if ( firemodel == BLAZE ) {
		create_output_table(out_aburned_area,		file_aburned_area_out,		blaze_columns);
		create_output_table(out_simfireanalysis,	file_simfireanalysis_out,	simfireanalysis_columns);
	} else if ( firemodel == GLOBFIRM ) {
		create_output_table(out_firert,			file_firert,			firert_columns);
	}

	create_output_table(out_runoff,			file_runoff,         runoff_columns);
	create_output_table(out_wetland_water_added, file_wetland_water_added, wetland_water_added_columns);
	create_output_table(out_speciesheights, file_speciesheights, speciesheights_columns);
	create_output_table(out_aiso,           file_aiso,           aiso_columns);
	create_output_table(out_amon,           file_amon,           amon_columns);
	create_output_table(out_amon_mt1,       file_amon_mt1,       amon_columns);
	create_output_table(out_amon_mt2,       file_amon_mt2,       amon_columns);

	create_output_table(out_nmass,          file_nmass,          nmass_columns);
	create_output_table(out_cton_leaf,      file_cton_leaf,      cton_columns);
	create_output_table(out_nsources,       file_nsources,       nsources_columns);
	create_output_table(out_npool,          file_npool,          npool_columns);
	create_output_table(out_nlitter,        file_nlitter,        nlitter_columns);
	create_output_table(out_nuptake,        file_nuptake,        nuptake_columns);
	create_output_table(out_vmaxnlim,       file_vmaxnlim,       vmaxnlim_columns);
	create_output_table(out_lambda,         file_lambda,         vmaxnlim_columns);
	create_output_table(out_nflux,          file_nflux,          nflux_columns);
	create_output_table(out_ngases,         file_ngases,         ngases_columns);
	create_output_table(out_soil_npool,		file_soil_npool,     soil_npool_columns);
	create_output_table(out_soil_nflux,	    file_soil_nflux,     soil_nflux_columns);

	create_output_table(out_wflux,          file_wflux,          wflux_columns);
	create_output_table(out_wpool,          file_wpool,          wpool_columns);
	create_output_table(out_wsats,          file_wsats,          soillyr_columns_short);
	create_output_table(out_wsats_depth,    file_wsats_depth,    soildepth_columns);

	// *** MONTHLY OUTPUT VARIABLES ***

	create_output_table(out_mnpp,           file_mnpp,           month_columns);
	create_output_table(out_mlai,           file_mlai,           month_columns);
	create_output_table(out_mgpp,           file_mgpp,           month_columns);
	create_output_table(out_mra,            file_mra,            month_columns);
	create_output_table(out_maet,           file_maet,           month_columns);
	create_output_table(out_mpet,           file_mpet,           month_columns);
	create_output_table(out_mevap,          file_mevap,          month_columns);
	create_output_table(out_mrunoff,        file_mrunoff,        month_columns_wide);
	create_output_table(out_mintercep,      file_mintercep,      month_columns);
	create_output_table(out_mrh,            file_mrh,            month_columns);
	create_output_table(out_mnee,           file_mnee,           month_columns);
	create_output_table(out_mnetmin,        file_mnetmin,        month_columns);
	create_output_table(out_mcton_leaf,     file_mcton_leaf,     pft_columns);
	create_output_table(out_mvmaxnlim,      file_mvmaxnlim,      pft_columns);
	create_output_table(out_mlambda,        file_mlambda,        pft_columns);
	create_output_table(out_mwmass,         file_mwmass,         month_columns_large);
	create_output_table(out_mimass,         file_mimass,         month_columns_large);
	create_output_table(out_mwcont_upper,   file_mwcont_upper,   month_columns);
	create_output_table(out_mwcont_lower,   file_mwcont_lower,   month_columns);
	create_output_table(out_miso,           file_miso,           month_columns_wide);
	create_output_table(out_mmon,           file_mmon,           month_columns_wide);
	create_output_table(out_mmon_mt1,       file_mmon_mt1,       month_columns_wide);
	create_output_table(out_mmon_mt2,       file_mmon_mt2,       month_columns_wide);
	create_output_table(out_mburned_area,   file_mburned_area_out, month_columns);
    
	// Methane
	create_output_table(out_mch4,           file_mch4,           month_columns);
	create_output_table(out_mch4diff,       file_mch4diff,       month_columns);
	create_output_table(out_mch4plan,       file_mch4plan,       month_columns);
	create_output_table(out_mch4ebull,      file_mch4ebull,      month_columns);
    
	// Snow
	create_output_table(out_msnow,          file_msnow,          month_columns);
	create_output_table(out_msnowcover,     file_msnowcover,     month_columns);
	create_output_table(out_mwmass_stand,   file_mwmass_stand,   month_columns_large);
	create_output_table(out_mwtp,           file_mwtp,           month_columns);
	create_output_table(out_mald,           file_mald,           mald_columns);
	create_output_table(out_msoiltempdepth25, file_msoiltempdepth25, month_columns);
}

/// Function for producing data file used to communicate information on stand structure
/** for 3D vegetation plot in Windows shell
 */
void output_vegetation(Gridcell& gridcell, Pftlist& pftlist) {
	
	// File for output of 3D vegetation structure (invoked by Windows shell only)
	plot3d_fileopen(); 
	
	if (plot3d_getfilehandle()) {

		int ival, p, npft_tree, npft_grass, npft_total;
		double grasslai;
		const bool FALSCH = false;
		const double rgb[3] = { -1, -1, -1 };
		char pftname[16];

		// Loop through Stands
		Gridcell::iterator gc_itr = gridcell.begin();
		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			npft_tree = npft_grass = 0;
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == TREE) npft_tree++;
				else if (pft.lifeform == GRASS) npft_grass++;
				pftlist.nextobj();
			}
			npft_total = npft_tree + npft_grass;
			fwrite(&npft_total, sizeof(int), 1, plot3d_getfilehandle());
			fwrite(&npft_tree, sizeof(int), 1, plot3d_getfilehandle());
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == TREE) {
					sprintf(pftname, "%s", (char*)(pft.name.left(15)));
					fwrite(pftname, sizeof(char), 16, plot3d_getfilehandle());
					fwrite(&FALSCH, sizeof(bool), 1, plot3d_getfilehandle());
					fwrite(&rgb, sizeof(double), 3, plot3d_getfilehandle());
					//fwrite(&pft.ifconifer, sizeof(bool), 1, plot3d_getfilehandle());
					//fwrite(pft.preferredrgb, sizeof(double), 3, plot3d_getfilehandle());
				}
				pftlist.nextobj();
			}
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if (pft.lifeform == GRASS) {
					sprintf(pftname, "%s", (char*)(pft.name.left(15)));
					fwrite(pftname, sizeof(char), 16, plot3d_getfilehandle());
					fwrite(&rgb, sizeof(double), 3, plot3d_getfilehandle());
					//fwrite(pft.preferredrgb, sizeof(double), 3, plot3d_getfilehandle());
				}
				pftlist.nextobj();
			}
			int npatch = stand.npatch();
			fwrite(&npatch, sizeof(int), 1, plot3d_getfilehandle());
			fwrite(&patcharea, sizeof(double), 1, plot3d_getfilehandle());
			for (p = 0; p<npatch; p++) {
				Patch& patch = stand[p];
				Vegetation& vegetation = patch.vegetation;
				grasslai = 0.0;
				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv = vegetation.getobj();
					if (indiv.pft.lifeform == TREE && indiv.alive) {
						ival = indiv.pft.id;
						fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle());
						ival = indiv.id;
						fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle());
						fwrite(&indiv.densindiv, sizeof(double), 1, plot3d_getfilehandle());
						fwrite(&indiv.height, sizeof(double), 1, plot3d_getfilehandle());
						fwrite(&indiv.crownarea, sizeof(double), 1, plot3d_getfilehandle());
					}
					else if (indiv.pft.lifeform == GRASS) grasslai += indiv.lai;
					vegetation.nextobj();
				}
				ival = -9999;
				fwrite(&ival, sizeof(int), 1, plot3d_getfilehandle()); // indicates no more cohorts in this patch
				if (grasslai<0) grasslai = 0.0;
				fwrite(&grasslai, sizeof(double), 1, plot3d_getfilehandle());
			}
			++gc_itr;
		} //while (gridcell.isobj) 

		plot3d_fileclose();

	} //if(out)

	plot3d();
}

/// Gets stand age structure to argument densindiv of dimensions [npft,nageclass]
/** First call with havedims=false to allocate memory and return nageclass,
  * then with havedims=true to get data into densindiv
  * Calling function is responsible for deallocating memory by:
  * delete[] densindiv;
*/
void get_stand_age_structure(Gridcell& gridcell,double* densindiv,int& nageclass,bool havedims) {

	int p, c;
	double active_fraction;

	if (!(vegmode == COHORT || vegmode == INDIVIDUAL)) {
		nageclass = 0;
		return;
	}

	if (havedims) {
		for (p = 0; p < npft; p++)
			for (c = 0; c < nageclass; c++)
				densindiv[p*nageclass + c] = 0.0;
	}
	else nageclass = 0;

	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();

		// Determine area fraction of stands where this pft is active:
		active_fraction = 0.0;
		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			if (stand.pft[pft.id].active) {
				active_fraction += stand.get_gridcell_fraction();
			}
			++gc_itr;
		}

		// Loop through Stands
		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			Standpft& standpft = stand.pft[pft.id];
			if (standpft.active) {

				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch = stand.getobj();
					Patchpft& patchpft = patch.pft[pft.id];
					Vegetation& vegetation = patch.vegetation;

					// Loop through individuals/cohorts

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv = vegetation.getobj();

						if (indiv.id != -1 && indiv.alive) {

							if (indiv.pft.id == pft.id) {

								// Age structure
								c = (int)(indiv.age / estinterval);

								if (havedims) {

									densindiv[pft.id*nageclass + c] += indiv.densindiv / (double)stand.npatch()* stand.get_gridcell_fraction() / active_fraction;

								}
								else if (c > nageclass) nageclass = c;
							}
						}
						vegetation.nextobj();
					}
					stand.nextobj();
				}
			}
			++gc_itr;
		}
		pftlist.nextobj();
	}
}

/// Local analogue of OutputRows::add_value for restricting output
/** Use to restrict output to specified range of years
  * (or other user-specified limitation)
  *
  * If only yearly output between, say 1961 and 1990 is requred, use:
  *  if (date.get_calendar_year() >= 1961 && date.get_calendar_year() <= 1990)
  *  (assuming the input module has set the first calendar year in the date object)	
  */
void outlimit(OutputRows& out, const Table& table, double d) {

	if (date.year >= nyear_spinup-100)
		out.add_value(table, d);
}

/// Output of simulation results at the end of each year
/** Output of simulation results at the end of each year, or for specific years in
  * the simulation of each stand or grid cell. 
  * This function does not have to provide any information to the framework.
  *
  * Restrict output to specific years in the local helper function outlimit().
  *
  * Changes in the structure of this function should be mirrored in outannual()
  * of the other output modules, e.g. MiscOutput::outannual().
  */
void CommonOutput::outannual(Gridcell& gridcell) {

	if (date.year < nyear_spinup-100)
		return;

	const double NOT_A_VALUE_SMALL = -9.999;
	const double NOT_A_VALUE_LARGE = -999.9;

	// Soil sections layers [mm]
	const int sslyr = 7;
	double soil_sections[sslyr] = { 300.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 50000.0 };

	// Soil sections temperature layers [mm]
	const int sstlyr = 12;
	double soil_temp_sections[sstlyr] = { 300.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0 };

	int c, m;
	double flux_veg, flux_repr, flux_soil, flux_fire, flux_est, flux_seed, flux_charvest;
	double c_fast, c_slow, c_harv_slow;

	double surflitterc, surflittern, surfhumusc, surfhumusn, soillitterc, soillittern, surfcwdc, surfcwdn, soilcwdc, soilcwdn, centuryc, centuryn, n_harv_slow, availn;
	double flux_nharvest, flux_nseed;
	double flux_NH3_soil, flux_NOx_soil, flux_N2O_soil, flux_N2_soil;
	double flux_NH3_fire, flux_NOx_fire, flux_N2O_fire, flux_N2_fire;
	double flux_ntot;

	double NH4_mass, NH3_mass, NO3_mass;
	double NO2_mass, NO_mass, N2O_mass, N2_mass;
	double NO_mass_inc, N2O_mass_inc, N2_mass_inc, NH3_mass_inc;
	double gross_nitrif, gross_denitrif, net_nitrif, net_denitrif;

	double soildepth[NSOILLAYERS_MAX] = { 0.0 };
	double soilc[NSOILLAYERS_MAX + 1] = { 0.0 };
	double soilc_depth[sslyr + 2] = { 0.0 };
	double soilc_cm[30] = { 0.0 };
	double soilcp_cm[30] = { 0.0 };
	double soilr[NSOILLAYERS_MAX + 1] = { 0.0 };
	double soilr_depth[sslyr + 2] = { 0.0 };
	double soiln_depth[sslyr + 2] = { 0.0 };
	double soilnetnmin_depth[sslyr + 1] = { 0.0 };
	double clitterinput_depth[sslyr + 2] = { 0.0 };
	double volume_org_depth[sslyr + 1] = { 0.0 };
	double volume_min_depth[sslyr + 1] = { 0.0 };
	double volume_root_depth[sslyr + 1] = { 0.0 };
	double temp_mod[NSOILLAYERS_MAX] = { 0.0 };
	double moist_mod[NSOILLAYERS_MAX] = { 0.0 };
	double tot_mod[NSOILLAYERS_MAX] = { 0.0 };
	double temp_mod_depth[sslyr + 1] = { 0.0 };
	double moist_mod_depth[sslyr + 1] = { 0.0 };
	double tot_mod_depth[sslyr + 1] = { 0.0 };
	double somc[NSOMPOOL][NSOILLAYERS_MAX] = { 0.0 };
	double somn[NSOMPOOL][NSOILLAYERS_MAX] = { 0.0 };
	double cdec[NSOILLAYERS_MAX][2] = { 0.0 };
	double cinc[NSOILLAYERS_MAX][2] = { 0.0 };
	double clitterinput[NSOILLAYERS_MAX + 1] = { 0.0 };
	double csomlitterinput[NSOMPOOL + 1] = { 0.0 };
	double wsats[NSOILLAYERS_MAX + 1] = { 0.0 };
	double wsats_depth[sslyr + 1] = { 0.0 };
	double soilc_tot, soilc_ald;
	double snowdays;

	// hold the monthly average across patches
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
	double mnetmin[12];
	double mcton_leaf[100][12] = { 0.0 };
	double mvmaxnlim[100][12] = { 0.0 };
	double mlambda[100][12] = { 0.0 };
	double mwmass[12];
	double mimass[12];
	double mwcont_upper[12];
	double mwcont_lower[12];
	double miso[12];
	double mmon[12];
	double mmon_mt1[12];
	double mmon_mt2[12];

	double mch4[12];
	double mch4_diff[12];
	double mch4_ebull[12];
	double mch4_plant[12];
	double msnowdepth[12];
	double msnowcover[12];
	double mwmass_stand[12];
	double mwtp[12];
	double mald[12];
	double mtemp[NSNOWLAYERS_MAX + NSOILLAYERS_MAX + NPADLAYERS][12];
	double mtemp_depth[sstlyr][12] = { 0.0 };
	double mtemp25cm[12];
	double mwmasslyr[NSOILLAYERS_MAX + 1][12];
	double mimasslyr[NSOILLAYERS_MAX + 1][12];
	double mdrainlyr[NSOILLAYERS_MAX + 1][12];
	double wsatslyr[NSOILLAYERS_MAX + 1];
	double mwmass_depth[sslyr + 1][12] = { 0.0 };
	double mimass_depth[sslyr + 1][12] = { 0.0 };
	double mwimass_depth[sslyr + 1][12] = { 0.0 };
	double mdrain_depth[sslyr + 1][12] = { 0.0 };

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();

	// The OutputRows object manages the next row of output for each
	// output table
	OutputRows out(output_channel, lon, lat, date.get_calendar_year());

	// guess2008 - reset monthly and annual sums across patches each year
	for (m = 0; m < 12; m++) {
		mnpp[m] = mlai[m] = mgpp[m] = mra[m] = maet[m] = mpet[m] = mevap[m] = mintercep[m] = mrunoff[m] = mrh[m] = mnee[m] = mnetmin[m] = mwmass[m] = mimass[m] = mwcont_upper[m] = mwcont_lower[m] = miso[m] = mmon[m] = mmon_mt1[m] = mmon_mt2[m] = 0.0;
		mch4[m] = mch4_diff[m] = mch4_ebull[m] = mch4_plant[m] = msnowdepth[m] = msnowcover[m] = mwmass_stand[m] = mwtp[m] = mald[m] = mtemp25cm[m] = 0.0;
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX + NSOILLAYERS_MAX + NPADLAYERS; lyr++) {
			mtemp[lyr][m] = 0.0;
		}
		for (int lyr = 0; lyr < NSOILLAYERS_MAX + 1; lyr++) {
			mwmasslyr[lyr][m] = 0.0;
			mimasslyr[lyr][m] = 0.0;
			mdrainlyr[lyr][m] = 0.0;
			if (m == 0) {
				wsatslyr[lyr] = 0.0;
			}
		}
	}

	double aaet, apet, aevap, arunoff, aintercep, awetland_wmass_added;
	aaet = apet = aevap = arunoff = aintercep = awetland_wmass_added = 0.0;

	double landcover_cmass[NLANDCOVERTYPES] = { 0.0 };
	double landcover_nmass[NLANDCOVERTYPES] = { 0.0 };
	double landcover_cmass_agb[NLANDCOVERTYPES] = { 0.0 };
	double landcover_cmass_leaf[NLANDCOVERTYPES] = { 0.0 };
	double landcover_nmass_leaf[NLANDCOVERTYPES] = { 0.0 };
	double landcover_cmass_veg[NLANDCOVERTYPES] = { 0.0 };
	double landcover_nmass_veg[NLANDCOVERTYPES] = { 0.0 };
	double landcover_clitter[NLANDCOVERTYPES] = { 0.0 };
	double landcover_nlitter[NLANDCOVERTYPES] = { 0.0 };
	double landcover_anpp[NLANDCOVERTYPES] = { 0.0 };
	double landcover_agpp[NLANDCOVERTYPES] = { 0.0 };
	double landcover_agpp_no[NLANDCOVERTYPES] = { 0.0 };
	double landcover_agpp_no_norm[NLANDCOVERTYPES] = { 0.0 };
	double landcover_agpp_ns[NLANDCOVERTYPES] = { 0.0 };
	double landcover_fpc[NLANDCOVERTYPES] = { 0.0 };
	double landcover_aaet[NLANDCOVERTYPES] = { 0.0 };
	double landcover_lai[NLANDCOVERTYPES] = { 0.0 };
	double landcover_densindiv_total[NLANDCOVERTYPES] = { 0.0 };
	double landcover_aiso[NLANDCOVERTYPES] = { 0.0 };
	double landcover_amon[NLANDCOVERTYPES] = { 0.0 };
	double landcover_amon_mt1[NLANDCOVERTYPES] = { 0.0 };
	double landcover_amon_mt2[NLANDCOVERTYPES] = { 0.0 };
	double landcover_nuptake[NLANDCOVERTYPES] = { 0.0 };
	double landcover_vmaxnlim[NLANDCOVERTYPES] = { 0.0 };
	double landcover_lambda[NLANDCOVERTYPES] = { 0.0 };

	double mean_standpft_cmass = 0.0;
	double mean_standpft_nmass = 0.0;
	double mean_standpft_cmass_agb = 0.0;
	double mean_standpft_cmass_leaf = 0.0;
	double mean_standpft_nmass_leaf = 0.0;
	double mean_standpft_cmass_veg = 0.0;
	double mean_standpft_nmass_veg = 0.0;
	double mean_standpft_clitter = 0.0;
	double mean_standpft_nlitter = 0.0;
	double mean_standpft_anpp = 0.0;
	double mean_standpft_agpp = 0.0;
	double mean_standpft_agpp_no = 0.0;
	double mean_standpft_agpp_no_norm = 0.0;
	double mean_standpft_agpp_ns = 0.0;
	double mean_standpft_fpc = 0.0;
	double mean_standpft_aaet = 0.0;
	double mean_standpft_lai = 0.0;
	double mean_standpft_densindiv_total = 0.0;
	double mean_standpft_heightindiv_total = 0.0;
	double mean_standpft_aiso = 0.0;
	double mean_standpft_amon = 0.0;
	double mean_standpft_amon_mt1 = 0.0;
	double mean_standpft_amon_mt2 = 0.0;
	double mean_standpft_nuptake = 0.0;
	double mean_standpft_vmaxnlim = 0.0;
	double mean_standpft_lambda = 0.0;
	double mean_standpft_mcmass_leaf[12] = { 0.0 };
	double mean_standpft_mnmass_leaf[12] = { 0.0 };
	double mean_standpft_mvmaxnlim[12] = { 0.0 };
	double mean_standpft_mlambda[12] = { 0.0 };
	double mean_standpft_mgpp_no_norm[12] = { 0.0 };

	double cmass_gridcell = 0.0;
	double nmass_gridcell = 0.0;
	double cmass_agb_gridcell = 0.0;
	double cmass_leaf_gridcell = 0.0;
	double nmass_leaf_gridcell = 0.0;
	double cmass_veg_gridcell = 0.0;
	double nmass_veg_gridcell = 0.0;
	double clitter_gridcell = 0.0;
	double nlitter_gridcell = 0.0;
	double csom_gridcell[NSOMPOOL + 1] = { 0.0 };
	double nsom_gridcell[NSOMPOOL + 1] = { 0.0 };
	double anpp_gridcell = 0.0;
	double agpp_gridcell = 0.0;
	double agpp_no_gridcell = 0.0;
	double agpp_no_norm_gridcell = 0.0;
	double agpp_ns_gridcell = 0.0;
	double fpc_gridcell = 0.0;
	double aaet_gridcell = 0.0;
	double lai_gridcell = 0.0;
	double surfrunoff_gridcell = 0.0;
	double drainrunoff_gridcell = 0.0;
	double baserunoff_gridcell = 0.0;
	double runoff_gridcell = 0.0;
	double wetland_wmass_added_gridcell = 0.0;
	double irrigation_gridcell = 0.0;
	double intercept_gridcell = 0.0;
	double evap_gridcell = 0.0;
	double aet_gridcell = 0.0;
	double dens_gridcell = 0.0;
	double firert_gridcell = 0.0;
	double burned_area_gridcell = 0.0;
	double aiso_gridcell = 0.0;
	double amon_gridcell = 0.0;
	double amon_mt1_gridcell = 0.0;
	double amon_mt2_gridcell = 0.0;
	double nuptake_gridcell = 0.0;
	double vmaxnlim_gridcell = 0.0;
	double lambda_gridcell = 0.0;
	double maxald_gridcell = 0.0;
	double permafrost_full_gridcell = 0.0;
	double permafrost_3m_gridcell = 0.0;

	double aNH4dep_gridcell = 0.0;
	double aNO3dep_gridcell = 0.0;
	double anfert_gridcell = 0.0;
	double anmin_gridcell = 0.0;
	double animm_gridcell = 0.0;
	double anfix_gridcell = 0.0;
	double n_min_leach_gridcell = 0.0;
	double n_org_leach_gridcell = 0.0;
	double c_org_leach_gridcell = 0.0;

	double wmass_snow_gridcell = 0.0;
	double imass_snow_gridcell = 0.0;
	double wmass_soil_gridcell = 0.0;
	double imass_soil_gridcell = 0.0;

	double standpft_cmass = 0.0;
	double standpft_nmass = 0.0;
	double standpft_cmass_agb = 0.0;
	double standpft_cmass_leaf = 0.0;
	double standpft_nmass_leaf = 0.0;
	double standpft_cmass_veg = 0.0;
	double standpft_nmass_veg = 0.0;
	double standpft_clitter = 0.0;
	double standpft_nlitter = 0.0;
	double standpft_anpp = 0.0;
	double standpft_agpp = 0.0;
	double standpft_agpp_no = 0.0;
	double standpft_agpp_no_norm = 0.0;
	double standpft_agpp_ns = 0.0;
	double standpft_fpc = 0.0;
	double standpft_aaet = 0.0;
	double standpft_lai = 0.0;
	double standpft_densindiv_total = 0.0;
	double standpft_heightindiv_total = 0.0;
	double standpft_aiso = 0.0;
	double standpft_amon = 0.0;
	double standpft_amon_mt1 = 0.0;
	double standpft_amon_mt2 = 0.0;
	double standpft_nuptake = 0.0;
	double standpft_vmaxnlim = 0.0;
	double standpft_lambda = 0.0;
	double standpft_mcmass_leaf[12] = { 0.0 };
	double standpft_mnmass_leaf[12] = { 0.0 };
	double standpft_mvmaxnlim[12] = { 0.0 };
	double standpft_mlambda[12] = { 0.0 };
	double standpft_mgpp_no_norm[12] = { 0.0 };

	// *** Loop through PFTs ***

	pftlist.firstobj();
	while (pftlist.isobj) {

		Pft& pft=pftlist.getobj();

		// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
		mean_standpft_cmass=0.0;
		mean_standpft_nmass=0.0;
		mean_standpft_cmass_agb=0.0;
		mean_standpft_cmass_leaf=0.0;
		mean_standpft_nmass_leaf=0.0;
		mean_standpft_cmass_veg=0.0;
		mean_standpft_nmass_veg=0.0;
		mean_standpft_clitter=0.0;
		mean_standpft_nlitter=0.0;
		mean_standpft_anpp=0.0;
		mean_standpft_agpp=0.0;
		mean_standpft_agpp_no=0.0;
		mean_standpft_agpp_no_norm=0.0;
		mean_standpft_agpp_ns=0.0;
		mean_standpft_fpc=0.0;
		mean_standpft_aaet=0.0;
		mean_standpft_lai=0.0;
		mean_standpft_densindiv_total=0.0;
		mean_standpft_aiso=0.0;
		mean_standpft_amon=0.0;
		mean_standpft_amon_mt1=0.0;
		mean_standpft_amon_mt2=0.0;
		mean_standpft_nuptake=0.0;
		mean_standpft_vmaxnlim = 0.0;
		mean_standpft_lambda = 0.0;

		for (int m = 0; m < 12; m++) {
			mean_standpft_mcmass_leaf[m] = 0.0;
			mean_standpft_mnmass_leaf[m] = 0.0;
			mean_standpft_mvmaxnlim[m] = 0.0;
			mean_standpft_mlambda[m] = 0.0;
			mean_standpft_mgpp_no_norm[m] = 0.0;
		}

		mean_standpft_heightindiv_total = 0.0;

		// Determine area fraction of stands where this pft is active:
		double active_fraction = 0.0;

		Gridcell::iterator gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			if(stand.pft[pft.id].active) {
				active_fraction += stand.get_gridcell_fraction();
			}

			++gc_itr;
		}

		// Loop through Stands
		gc_itr = gridcell.begin();

		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;

			Standpft& standpft=stand.pft[pft.id];
			if(standpft.active) {

				// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
				standpft_cmass=0.0;
				standpft_nmass=0.0;
				standpft_cmass_agb=0.0;
				standpft_cmass_leaf=0.0;
				standpft_nmass_leaf=0.0;
				standpft_cmass_veg=0.0;
				standpft_nmass_veg=0.0;
				standpft_clitter=0.0;
				standpft_nlitter=0.0;
				standpft_anpp=0.0;
				standpft_agpp=0.0;
				standpft_agpp_no=0.0;
				standpft_agpp_no_norm=0.0;
				standpft_agpp_ns=0.0;
				standpft_fpc=0.0;
				standpft_aaet=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total = 0.0;
				standpft_heightindiv_total = 0.0;
				standpft_aiso=0.0;
				standpft_amon=0.0;
				standpft_amon_mt1=0.0;
				standpft_amon_mt2=0.0;
				standpft_nuptake=0.0;
				standpft_vmaxnlim=0.0;
				standpft_lambda=0.0;

				for (int m = 0; m < 12; m++) {
					standpft_mcmass_leaf[m] = 0.0;
					standpft_mnmass_leaf[m] = 0.0;
					standpft_mvmaxnlim[m] = 0.0;
					standpft_mlambda[m] = 0.0;
					standpft_mgpp_no_norm[m] = 0.0;
				}

				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch = stand.getobj();
					Patchpft& patchpft = patch.pft[pft.id];
					Vegetation& vegetation = patch.vegetation;
					Fluxes& fluxes = patch.fluxes;

					standpft_anpp += fluxes.get_annual_flux(Fluxes::NPP, pft.id);
					standpft_agpp += fluxes.get_annual_flux(Fluxes::GPP, pft.id);
					standpft_agpp_no += fluxes.get_annual_flux(Fluxes::GPPno, pft.id);
					standpft_agpp_ns += fluxes.get_annual_flux(Fluxes::GPPns, pft.id);
					standpft_aiso += fluxes.get_annual_flux(Fluxes::ISO, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_APIN, pft.id); 					
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_LIMO, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_TRIC, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_BPIN, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_MYRC, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_SABI, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_CAMP, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_TBOC, pft.id);
					standpft_amon += fluxes.get_annual_flux(Fluxes::MT_OTHR, pft.id);
					standpft_amon_mt1 += fluxes.get_annual_flux(Fluxes::MT_APIN, pft.id);
					standpft_amon_mt1 += fluxes.get_annual_flux(Fluxes::MT_LIMO, pft.id);
					standpft_amon_mt1 += fluxes.get_annual_flux(Fluxes::MT_TRIC, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_BPIN, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_MYRC, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_SABI, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_CAMP, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_TBOC, pft.id);
					standpft_amon_mt2 += fluxes.get_annual_flux(Fluxes::MT_OTHR, pft.id);
					
					standpft_clitter += patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr;
					standpft_nlitter += patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap + patchpft.nmass_litter_heart;

					for (int m = 0; m < 12; m++) {
						standpft_mcmass_leaf[m] += patchpft.mcmass_leaf[m];
						standpft_mnmass_leaf[m] += patchpft.mnmass_leaf[m];
					}

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv=vegetation.getobj();

						if (indiv.id!=-1 && indiv.alive) {

							if (indiv.pft.id==pft.id) {

								standpft_cmass_leaf += indiv.cmass_leaf;
								standpft_cmass += indiv.ccont();
								standpft_nmass += indiv.ncont();
								standpft_nmass_leaf += indiv.cmass_leaf / indiv.cton_leaf_aavr;
								standpft_nmass_veg += indiv.nmass_veg;
								standpft_cmass_agb += indiv.cmass_leaf + indiv.cmass_wood() * indiv.pft.ag_woody_frac;
								standpft_fpc += indiv.fpc;
								standpft_aaet += indiv.aaet;
								standpft_lai += indiv.lai;
								if (pft.lifeform==TREE) {	
									standpft_densindiv_total += indiv.densindiv;
									standpft_heightindiv_total += indiv.height * indiv.densindiv;
								}
								standpft_agpp_no_norm += indiv.agpp_no;
								standpft_vmaxnlim += indiv.avmaxnlim * indiv.agpp_no;
								standpft_lambda += indiv.alambda * indiv.agpp_no;
								standpft_nuptake += indiv.anuptake;
								
								for (int m = 0; m < 12; m++) {
									standpft_mvmaxnlim[m] += indiv.mvmaxnlim[m] * indiv.mgpp_no[m];
									standpft_mlambda[m] += indiv.mlambda[m] * indiv.mgpp_no[m];
									standpft_mgpp_no_norm[m] += indiv.mgpp_no[m];
								}

								if(pft.landcover == CROPLAND) {
									standpft_cmass_veg += indiv.cmass_leaf + indiv.cmass_root;
									if(indiv.cropindiv) {
										standpft_cmass_veg += indiv.cropindiv->cmass_ho + indiv.cropindiv->cmass_agpool + indiv.cropindiv->cmass_stem;
										standpft_nmass_leaf += indiv.cropindiv->ynmass_leaf + indiv.cropindiv->ynmass_dead_leaf;
										standpft_nmass_veg += indiv.cropindiv->ynmass_leaf + indiv.cropindiv->ynmass_dead_leaf + indiv.cropindiv->ynmass_root + indiv.cropindiv->ynmass_ho + indiv.cropindiv->ynmass_agpool;
									}
								}
								else {

									standpft_cmass_veg += indiv.cmass_veg;
									standpft_nmass_veg += indiv.nmass_veg;
								}
							}

						} // alive?
						vegetation.nextobj();
					}

					stand.nextobj();
				} // end of patch loop

				// Update to stand
				double fraction_of_stand = 1.0 / (double)stand.npatch();

				standpft_cmass*=fraction_of_stand;
				standpft_nmass*=fraction_of_stand;
				standpft_cmass_agb*=fraction_of_stand;
				standpft_cmass_leaf*=fraction_of_stand;
				standpft_nmass_leaf*=fraction_of_stand;
				standpft_cmass_veg*=fraction_of_stand;
				standpft_nmass_veg*=fraction_of_stand;
				standpft_clitter*=fraction_of_stand;
				standpft_nlitter*=fraction_of_stand;
				standpft_anpp*=fraction_of_stand;
				standpft_agpp*=fraction_of_stand;
				standpft_agpp_no*=fraction_of_stand;
				standpft_agpp_no_norm*=fraction_of_stand;
				standpft_agpp_ns*=fraction_of_stand;
				standpft_fpc*=fraction_of_stand;
				standpft_aaet*=fraction_of_stand;
				standpft_lai*=fraction_of_stand;
				standpft_densindiv_total*=fraction_of_stand;
				standpft_aiso*=fraction_of_stand; // missing above!
				standpft_amon*=fraction_of_stand; // missing above!
				standpft_amon_mt1*=fraction_of_stand;
				standpft_amon_mt2*=fraction_of_stand;
				standpft_nuptake*=fraction_of_stand;
				standpft_vmaxnlim*=fraction_of_stand;
				standpft_lambda*=fraction_of_stand;
				standpft_heightindiv_total*=fraction_of_stand;

				for (int m = 0; m < 12; m++) {
					standpft_mcmass_leaf[m]*=fraction_of_stand;
					standpft_mnmass_leaf[m]*=fraction_of_stand;
				}

				if (!negligible(standpft_agpp_no_norm)) {
					standpft_vmaxnlim/=standpft_agpp_no_norm;
					standpft_lambda/=standpft_agpp_no_norm;
				}

				// Update landcover totals
				double fraction_of_landcover = stand.get_landcover_fraction();

				landcover_cmass[stand.landcover]+=standpft_cmass*fraction_of_landcover;
				landcover_nmass[stand.landcover]+=standpft_nmass*fraction_of_landcover;
				landcover_cmass_agb[stand.landcover]+=standpft_cmass_agb*fraction_of_landcover;
				landcover_cmass_leaf[stand.landcover]+=standpft_cmass_leaf*fraction_of_landcover;
				landcover_nmass_leaf[stand.landcover]+=standpft_nmass_leaf*fraction_of_landcover;
				landcover_cmass_veg[stand.landcover]+=standpft_cmass_veg*fraction_of_landcover;
				landcover_nmass_veg[stand.landcover]+=standpft_nmass_veg*fraction_of_landcover;
				landcover_clitter[stand.landcover]+=standpft_clitter*fraction_of_landcover;
				landcover_nlitter[stand.landcover]+=standpft_nlitter*fraction_of_landcover;
				landcover_anpp[stand.landcover]+=standpft_anpp*fraction_of_landcover;
				landcover_agpp[stand.landcover]+=standpft_agpp*fraction_of_landcover;
				landcover_agpp_no[stand.landcover]+=standpft_agpp_no*fraction_of_landcover;
				landcover_agpp_no_norm[stand.landcover]+=standpft_agpp_no_norm*fraction_of_landcover;
				landcover_agpp_ns[stand.landcover]+=standpft_agpp_ns*fraction_of_landcover;
				if(!pft.isintercropgrass) {
					landcover_fpc[stand.landcover]+=standpft_fpc*fraction_of_landcover;
					landcover_lai[stand.landcover]+=standpft_lai*fraction_of_landcover;
				}
				landcover_aaet[stand.landcover]+=standpft_aaet*fraction_of_landcover;
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*fraction_of_landcover;
				landcover_aiso[stand.landcover]+=standpft_aiso*fraction_of_landcover;
				landcover_amon[stand.landcover]+=standpft_amon*fraction_of_landcover;
				landcover_amon_mt1[stand.landcover]+=standpft_amon_mt1*fraction_of_landcover;
				landcover_amon_mt2[stand.landcover]+=standpft_amon_mt2*fraction_of_landcover;
				landcover_nuptake[stand.landcover]+=standpft_nuptake*fraction_of_landcover;
				landcover_vmaxnlim[stand.landcover]+=standpft_vmaxnlim*standpft_agpp_no_norm*fraction_of_landcover;
				landcover_lambda[stand.landcover]+=standpft_lambda*standpft_agpp_no_norm*fraction_of_landcover;

				// Update gridcell totals
				double fraction_of_gridcell = stand.get_gridcell_fraction();

				// Update pft means for active stands
				if(active_fraction) {
					mean_standpft_cmass += standpft_cmass * fraction_of_gridcell / active_fraction;
					mean_standpft_nmass += standpft_nmass * fraction_of_gridcell / active_fraction;
					mean_standpft_cmass_agb += standpft_cmass_agb * fraction_of_gridcell / active_fraction;
					mean_standpft_cmass_leaf += standpft_cmass_leaf * fraction_of_gridcell / active_fraction;
					mean_standpft_nmass_leaf += standpft_nmass_leaf * fraction_of_gridcell / active_fraction;
					mean_standpft_cmass_veg += standpft_cmass_veg * fraction_of_gridcell / active_fraction;
					mean_standpft_nmass_veg += standpft_nmass_veg * fraction_of_gridcell / active_fraction;
					mean_standpft_clitter += standpft_clitter * fraction_of_gridcell / active_fraction;
					mean_standpft_nlitter += standpft_nlitter * fraction_of_gridcell / active_fraction;
					mean_standpft_anpp += standpft_anpp * fraction_of_gridcell / active_fraction;
					mean_standpft_agpp += standpft_agpp * fraction_of_gridcell / active_fraction;
					mean_standpft_agpp_no += standpft_agpp_no * fraction_of_gridcell / active_fraction;
					mean_standpft_agpp_no_norm += standpft_agpp_no_norm * fraction_of_gridcell / active_fraction;
					mean_standpft_agpp_ns += standpft_agpp_ns * fraction_of_gridcell / active_fraction;
					mean_standpft_fpc += standpft_fpc * fraction_of_gridcell / active_fraction;
					mean_standpft_aaet += standpft_aaet * fraction_of_gridcell / active_fraction;
					mean_standpft_lai += standpft_lai * fraction_of_gridcell / active_fraction;
					mean_standpft_densindiv_total += standpft_densindiv_total * fraction_of_gridcell / active_fraction;
					mean_standpft_heightindiv_total += standpft_heightindiv_total * fraction_of_gridcell / active_fraction;
					mean_standpft_aiso += standpft_aiso * fraction_of_gridcell / active_fraction;
					mean_standpft_amon += standpft_amon * fraction_of_gridcell / active_fraction;
					mean_standpft_amon_mt1 += standpft_amon_mt1 * fraction_of_gridcell / active_fraction;
					mean_standpft_amon_mt2 += standpft_amon_mt2 * fraction_of_gridcell / active_fraction;
					mean_standpft_nuptake += standpft_nuptake * fraction_of_gridcell / active_fraction;
					mean_standpft_vmaxnlim += standpft_vmaxnlim * fraction_of_gridcell / active_fraction;
					mean_standpft_lambda += standpft_lambda * fraction_of_gridcell / active_fraction;
					for (int m = 0; m < 12; m++) {
						mean_standpft_mcmass_leaf[m] += standpft_mcmass_leaf[m] * fraction_of_gridcell / active_fraction;
						mean_standpft_mnmass_leaf[m] += standpft_mnmass_leaf[m] * fraction_of_gridcell / active_fraction;
						mean_standpft_mvmaxnlim[m] += standpft_mvmaxnlim[m] * fraction_of_gridcell / active_fraction;
						mean_standpft_mlambda[m] += standpft_mlambda[m] * fraction_of_gridcell / active_fraction;
						mean_standpft_mgpp_no_norm[m] += standpft_mgpp_no_norm[m] * fraction_of_gridcell / active_fraction;
					}
				}

				// Update stand totals
				stand.anpp += standpft_anpp;
				stand.cmass += standpft_cmass;

				cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
				nmass_gridcell+=standpft_nmass*fraction_of_gridcell;
				cmass_agb_gridcell+=standpft_cmass_agb*fraction_of_gridcell;
				cmass_leaf_gridcell+=standpft_cmass_leaf*fraction_of_gridcell;
				nmass_leaf_gridcell+=standpft_nmass_leaf*fraction_of_gridcell;
				cmass_veg_gridcell+=standpft_cmass_veg*fraction_of_gridcell;
				nmass_veg_gridcell+=standpft_nmass_veg*fraction_of_gridcell;
				clitter_gridcell+=standpft_clitter*fraction_of_gridcell;
				nlitter_gridcell+=standpft_nlitter*fraction_of_gridcell;
				anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
				agpp_gridcell+=standpft_agpp*fraction_of_gridcell;
				agpp_no_gridcell+=standpft_agpp_no*fraction_of_gridcell;
				agpp_no_norm_gridcell+=standpft_agpp_no_norm*fraction_of_gridcell;
				agpp_ns_gridcell+=standpft_agpp_ns*fraction_of_gridcell;
				if(!pft.isintercropgrass) {
					fpc_gridcell+=standpft_fpc*fraction_of_gridcell;
					lai_gridcell+=standpft_lai*fraction_of_gridcell;
				}
				aaet_gridcell+=standpft_aaet*fraction_of_gridcell;
				dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
				aiso_gridcell+=standpft_aiso*fraction_of_gridcell;
				amon_gridcell+=standpft_amon*fraction_of_gridcell;
				amon_mt1_gridcell+=standpft_amon_mt1*fraction_of_gridcell;
				amon_mt2_gridcell+=standpft_amon_mt2*fraction_of_gridcell;
				nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
				vmaxnlim_gridcell+=standpft_vmaxnlim*standpft_agpp_no_norm*fraction_of_gridcell;
				lambda_gridcell += standpft_lambda*standpft_agpp_no_norm*fraction_of_gridcell;
				
				// Graphical output every PLOT_INTERVAL years
				// (Windows shell only - "plot" statements have no effect otherwise)
				if (!(date.year%PLOT_INTERVAL)) {
					plot("C mass [kgC/m2]",pft.name,date.year,mean_standpft_cmass);
					plot("NPP [kgC/m2/yr]",pft.name,date.year,mean_standpft_anpp);
					plot("LAI [m2/m2]",pft.name,date.year,mean_standpft_lai);
					if (pft.lifeform == TREE) plot("Dens [indiv/ha]",pft.name,date.year,mean_standpft_densindiv_total*M2_PER_HA);
					if (mean_standpft_cmass_leaf > 0.0 && ifnlim) {
						plot("Vmax N lim",pft.name,date.year,mean_standpft_vmaxnlim);
						plot("leaf C:N [kgC/kg N]",pft.name,date.year,mean_standpft_cmass_leaf/mean_standpft_nmass_leaf);
					}
				}

			}//if(active)
			++gc_itr;
		}//End of loop through stands

		// Print PFT sums to files

		double mean_standpft_cton_leaf = limited_cton(mean_standpft_cmass_leaf, mean_standpft_nmass_leaf, pft.cton_leaf_max);

		outlimit(out,out_cmass,     mean_standpft_cmass);
		outlimit(out,out_cmass_agb, mean_standpft_cmass_agb);
		outlimit(out,out_anpp,      mean_standpft_anpp);
		outlimit(out,out_agpp,      mean_standpft_agpp);
		outlimit(out,out_agpp_no,   mean_standpft_agpp_no);
		outlimit(out,out_agpp_ns,   mean_standpft_agpp_ns);
		outlimit(out,out_fpc,       mean_standpft_fpc);
		outlimit(out,out_aaet,      mean_standpft_aaet);
		outlimit(out,out_clitter,   mean_standpft_clitter);
		outlimit(out,out_dens,      mean_standpft_densindiv_total);
		outlimit(out,out_lai,       mean_standpft_lai);
		outlimit(out,out_aiso,      mean_standpft_aiso);
		outlimit(out,out_amon,      mean_standpft_amon);
		outlimit(out,out_amon_mt1,  mean_standpft_amon_mt1);
		outlimit(out,out_amon_mt2,  mean_standpft_amon_mt2);
		outlimit(out,out_nmass,     (mean_standpft_nmass + mean_standpft_nlitter) * M2_PER_HA);
		outlimit(out,out_cton_leaf, mean_standpft_cton_leaf);
		outlimit(out,out_vmaxnlim,  mean_standpft_vmaxnlim);
		outlimit(out,out_lambda,    mean_standpft_lambda);
		outlimit(out,out_nuptake,   mean_standpft_nuptake * M2_PER_HA);
		outlimit(out,out_nlitter,   mean_standpft_nlitter * M2_PER_HA);

		// print species heights
		double height = 0.0;
		if (mean_standpft_densindiv_total > 0.0)
			height = mean_standpft_heightindiv_total / mean_standpft_densindiv_total;

		outlimit(out,out_speciesheights, height);

		for (int m = 0; m < 12; m++) {
			if (!negligible(mean_standpft_mgpp_no_norm[m])) {
				mvmaxnlim[pft.id][m] = mean_standpft_mvmaxnlim[m] / mean_standpft_mgpp_no_norm[m];
				mlambda[pft.id][m] = mean_standpft_mlambda[m] / mean_standpft_mgpp_no_norm[m];
			}

			mcton_leaf[pft.id][m] = limited_cton(mean_standpft_mcmass_leaf[m], mean_standpft_mnmass_leaf[m], pft.cton_leaf_max);
		}

		pftlist.nextobj();

	} // *** End of PFT loop ***

	flux_veg = flux_repr = flux_soil = flux_fire = flux_est = flux_seed = flux_charvest = 0.0;

	// guess2008 - carbon pools
	c_fast = c_slow = c_harv_slow = 0.0;

	surflitterc = surflittern = surfhumusc = surfhumusn = soillitterc = soillittern = surfcwdc = surfcwdn = soilcwdc = soilcwdn = centuryc = centuryn = n_harv_slow = availn = 0.0;
	aNH4dep_gridcell = aNO3dep_gridcell = anfert_gridcell = anmin_gridcell = animm_gridcell = anfix_gridcell = 0.0;
	n_org_leach_gridcell = n_min_leach_gridcell = c_org_leach_gridcell = 0.0;
	flux_NH3_soil = flux_NOx_soil = flux_N2O_soil = flux_N2_soil = 0.0;
	flux_NH3_fire = flux_NOx_fire = flux_N2O_fire = flux_N2_fire = 0.0;
	flux_ntot = flux_nharvest = flux_nseed = 0.0;
	soilc_ald = soilc_tot = 0.0;
	snowdays = 0.0;

	// Nitrogen in soil
	NH4_mass = NH3_mass = NO3_mass = 0.0;
	NO_mass = NO2_mass = N2O_mass = N2_mass = 0.0;
	NH3_mass_inc = NO_mass_inc = N2O_mass_inc = N2_mass_inc = 0.0;
	gross_nitrif = gross_denitrif = net_nitrif = net_denitrif= 0.0;

	double c_org_leach_lc[NLANDCOVERTYPES];

	for (int i=0; i<NLANDCOVERTYPES; i++) {
		c_org_leach_lc[i]=0.0;
	}

	// Sum C fluxes, dead C pools and runoff across patches

	Gridcell::iterator gc_itr = gridcell.begin();

	// Loop through Stands
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();

		//Loop through Patches
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			Fluxes& fluxes = patch.fluxes;
			Soilcolumn& soilcolumn = patch.soilcolumn;
			Snow& snow = soilcolumn.snow;

			double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

			flux_veg+=-fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
			flux_repr+=-fluxes.get_annual_flux(Fluxes::REPRC)*to_gridcell_average;
			flux_soil+=fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
			flux_fire+=fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
			flux_est+=fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
			flux_seed+=fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
			flux_charvest+=fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;

			flux_nseed+=fluxes.get_annual_flux(Fluxes::SEEDN)*to_gridcell_average;
			flux_nharvest+=fluxes.get_annual_flux(Fluxes::HARVESTN)*to_gridcell_average;
			//Nfluxes
			flux_NH3_fire+=fluxes.get_annual_flux(Fluxes::NH3_FIRE)*to_gridcell_average;
			flux_NOx_fire+=fluxes.get_annual_flux(Fluxes::NOx_FIRE)*to_gridcell_average;
			flux_N2O_fire+=fluxes.get_annual_flux(Fluxes::N2O_FIRE)*to_gridcell_average;
			flux_N2_fire+=fluxes.get_annual_flux(Fluxes::N2_FIRE)*to_gridcell_average;
			flux_NH3_soil+=fluxes.get_annual_flux(Fluxes::NH3_SOIL)*to_gridcell_average;
			flux_NOx_soil+=fluxes.get_annual_flux(Fluxes::NO_SOIL)*to_gridcell_average;
			flux_N2O_soil+=fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;
			flux_N2_soil+=fluxes.get_annual_flux(Fluxes::N2_SOIL)*to_gridcell_average;

			//Soil N flux from ntransform.cpp
			NH3_mass_inc+=fluxes.get_annual_flux(Fluxes::NH3_SOIL)*to_gridcell_average;
			NO_mass_inc+=fluxes.get_annual_flux(Fluxes::NO_SOIL)*to_gridcell_average;
			N2O_mass_inc+=fluxes.get_annual_flux(Fluxes::N2O_SOIL)*to_gridcell_average;
			N2_mass_inc+=fluxes.get_annual_flux(Fluxes::N2_SOIL)*to_gridcell_average;
			gross_nitrif+=fluxes.get_annual_flux(Fluxes::GROSS_NITRIF)*to_gridcell_average;
			net_nitrif+=fluxes.get_annual_flux(Fluxes::NET_NITRIF)*to_gridcell_average;
			gross_denitrif+=fluxes.get_annual_flux(Fluxes::GROSS_DENITRIF)*to_gridcell_average;
			net_denitrif+=fluxes.get_annual_flux(Fluxes::NET_DENITRIF)*to_gridcell_average;

			//Fire N flux
			flux_ntot+=(fluxes.get_annual_flux(Fluxes::NH3_FIRE) +
						fluxes.get_annual_flux(Fluxes::NOx_FIRE) +
						fluxes.get_annual_flux(Fluxes::N2O_FIRE) +
						fluxes.get_annual_flux(Fluxes::N2_FIRE) +
			//Soil N flux
					fluxes.get_annual_flux(Fluxes::NH3_SOIL) +
					fluxes.get_annual_flux(Fluxes::NO_SOIL) +
					fluxes.get_annual_flux(Fluxes::N2O_SOIL) +
					fluxes.get_annual_flux(Fluxes::N2_SOIL)) * to_gridcell_average;

			//Sum slow pools of harvested products
			if(run_landcover && ifslowharvestpool) {
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft=patch.pft[q];
					c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
					n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
				}
			}

			surfrunoff_gridcell += fluxes.get_annual_flux(Fluxes::WRUNOFF_SURF)*to_gridcell_average;
			drainrunoff_gridcell += fluxes.get_annual_flux(Fluxes::WRUNOFF_DRAIN)*to_gridcell_average;
			baserunoff_gridcell += fluxes.get_annual_flux(Fluxes::WRUNOFF_BASE)*to_gridcell_average;
			runoff_gridcell += (fluxes.get_annual_flux(Fluxes::WRUNOFF_SURF) + fluxes.get_annual_flux(Fluxes::WRUNOFF_DRAIN) + fluxes.get_annual_flux(Fluxes::WRUNOFF_BASE)) * to_gridcell_average;
			wetland_wmass_added_gridcell += fluxes.get_annual_flux(Fluxes::WADDED) * to_gridcell_average;
			irrigation_gridcell += fluxes.get_annual_flux(Fluxes::WIRRI) * to_gridcell_average;
			intercept_gridcell += fluxes.get_annual_flux(Fluxes::WINTERCEP) * to_gridcell_average;
			evap_gridcell += fluxes.get_annual_flux(Fluxes::WEVAP) * to_gridcell_average;
			aet_gridcell += fluxes.get_annual_flux(Fluxes::WAET) * to_gridcell_average;

			wmass_snow_gridcell += snow.wmass() * to_gridcell_average;
			imass_snow_gridcell += snow.imass() * to_gridcell_average;
			wmass_soil_gridcell += soilcolumn.get_sum(WMASS) * to_gridcell_average;
			imass_soil_gridcell += soilcolumn.get_sum(IMASS) * to_gridcell_average;

			// Fire return time
			if (!patch.has_fires() || patch.fireprob < 0.001) {
				firert_gridcell += 1000.0 * to_gridcell_average; // Set a limit of 1000 years
			}
			else {
				firert_gridcell += (1.0 / patch.fireprob) * to_gridcell_average;
				burned_area_gridcell += patch.fireprob * to_gridcell_average;
			}

			aNH4dep_gridcell += gridcell.aNH4dep * to_gridcell_average;
			aNO3dep_gridcell += gridcell.aNO3dep * to_gridcell_average;
			anfert_gridcell += patch.anfert * to_gridcell_average;
			anmin_gridcell += soilcolumn.get_sum(ANMIN) * to_gridcell_average;
			animm_gridcell += soilcolumn.get_sum(ANIMM) * to_gridcell_average;
			anfix_gridcell += soilcolumn.anfix * to_gridcell_average;
			n_min_leach_gridcell += fluxes.get_annual_flux(Fluxes::LEACHMINN) * to_gridcell_average;
			c_org_leach_gridcell += fluxes.get_annual_flux(Fluxes::LEACHORGC) * to_gridcell_average;
			n_org_leach_gridcell += fluxes.get_annual_flux(Fluxes::LEACHORGN) * to_gridcell_average;
			availn += (soilcolumn.nmass_avail() + soilcolumn.snow.nmass_avail()) * to_gridcell_average;

			c_org_leach_lc[stand.landcover] += fluxes.get_annual_flux(Fluxes::LEACHORGC) * to_gridcell_average;

			for (int p = 0; p < NSOMPOOL; p++) {
				csomlitterinput[p] += soilcolumn.csom_litter_input[p] * to_gridcell_average;
				csomlitterinput[NSOMPOOL] += soilcolumn.csom_litter_input[p] * to_gridcell_average;
			}

			clitterinput[0] += soilcolumn.asurflitterinput * to_gridcell_average;

			for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

				Soillayer& soillayer = soilcolumn.soillayers[lyr];

				NH4_mass += soillayer.NH4_mass * to_gridcell_average;
				NO3_mass += soillayer.NO3_mass * to_gridcell_average;
				NO2_mass += soillayer.NO2_mass * to_gridcell_average;
				NO_mass += soillayer.NO_mass * to_gridcell_average;
				N2O_mass += soillayer.N2O_mass * to_gridcell_average;
				N2_mass += soillayer.N2_mass * to_gridcell_average;

				// Soil layer depths
				soildepth[lyr] += soilcolumn.soillayerdepths[lyr] * to_gridcell_average;

				cdec[lyr][0] += soillayer.acdec[0] * to_gridcell_average;
				cdec[lyr][1] += soillayer.acdec[1] * to_gridcell_average;
				cinc[lyr][0] += soillayer.acinc[0] * to_gridcell_average;
				cinc[lyr][1] += soillayer.acinc[1] * to_gridcell_average;

				clitterinput[lyr + 1] += soillayer.get(ALITTERINPUT) * to_gridcell_average;

				c_fast += soillayer.cpool_fast * to_gridcell_average;
				c_slow += soillayer.cpool_slow * to_gridcell_average;

				for (int p = 0; p < NSOMPOOL; p++) {

					Sompool& sompool = soillayer.sompool[p];

					if (p == SURFMETA || p == SURFSTRUCT) {
						surflitterc += sompool.cmass * to_gridcell_average;
						surflittern += sompool.nmass * to_gridcell_average;
					}
					else if (p == SURFFWD || p == SURFCWD) {
						surfcwdc += sompool.cmass * to_gridcell_average;
						surfcwdn += sompool.nmass * to_gridcell_average;
					}
					else if (p == SURFMICRO || p == SURFHUMUS) {
						surfhumusc += sompool.cmass * to_gridcell_average;
						surfhumusn += sompool.nmass * to_gridcell_average;
					}
					else if (p == SOILMETA || p == SOILSTRUCT) {
						soillitterc += sompool.cmass * to_gridcell_average;
						soillittern += sompool.nmass * to_gridcell_average;
					}
					else if (p == SOILFWD || p == SOILCWD) {
						soilcwdc += sompool.cmass * to_gridcell_average;
						soilcwdn += sompool.nmass * to_gridcell_average;
					}
					else {
						centuryc += sompool.cmass * to_gridcell_average;
						centuryn += sompool.nmass * to_gridcell_average;
					}

					somc[p][lyr] += sompool.cmass * to_gridcell_average;
					somn[p][lyr] += sompool.nmass * to_gridcell_average;
				}
			}

			for (int lyr = soilcolumn.nsoillayers; lyr < NSOILLAYERS_MAX; lyr++) {
				// Soil layer depths
				soildepth[lyr] += soilcolumn.soillayerdepths[soilcolumn.nsoillayers - 1] * to_gridcell_average;
			}

			// Soil depths total C mass excluding surface litter C mass and heteorotrophic respiration down to depth
			double surflitc = soilcolumn.get_surflittercmass();
			double surflitn = soilcolumn.get_surflitternmass();
			double surflitr = soilcolumn.asurflitterresp;
			double surflitinput = soilcolumn.asurflitterinput;

			soilc_depth[0] += surflitc * to_gridcell_average;
			soiln_depth[0] += surflitn * to_gridcell_average;
			soilr_depth[0] += surflitr * to_gridcell_average;
			clitterinput_depth[0] += surflitinput * to_gridcell_average;

			// Soil to specific depths
			double depth_prev = 0.0;
			double soilc_prev = surflitc;
			double soiln_prev = surflitn;
			double soilr_prev = surflitr;

			for (int lyr = 0; lyr < sslyr; lyr++) {
				double depth_this = min(soilcolumn.depth, soil_sections[lyr]);
				double depth_lyr = depth_this - depth_prev;

				if (depth_lyr > 1.0e-12) {

					soilc_depth[lyr + 1] += (soilcolumn.var_depth(depth_prev, depth_this, CCONT) - soilc_prev) * to_gridcell_average;
					soiln_depth[lyr + 1] += (soilcolumn.var_depth(depth_prev, depth_this, NCONT) - soiln_prev)* to_gridcell_average;
					soilr_depth[lyr + 1] += (soilcolumn.var_depth(depth_prev, depth_this, ARH) - soilr_prev)* to_gridcell_average;
					clitterinput_depth[lyr + 1] += soilcolumn.var_depth(depth_prev, depth_this, ALITTERINPUT) * to_gridcell_average;
					soilnetnmin_depth[lyr] += (soilcolumn.var_depth(depth_prev, depth_this, ANMIN) - soilcolumn.var_depth(depth_prev, depth_this, ANIMM)) * to_gridcell_average;
					volume_org_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, VOL_ORG) * to_gridcell_average;
					volume_min_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, VOL_MIN) * to_gridcell_average;
					volume_root_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, VOL_ROOT) * to_gridcell_average;
					double wsats = soilcolumn.var_depth(depth_prev, depth_this, WSATS);
					wsats_depth[lyr] += wsats / depth_lyr * to_gridcell_average;
					temp_mod_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, MODTEMP, DEPTH) * to_gridcell_average;
					moist_mod_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, MODMOIS, DEPTH) * to_gridcell_average;
					tot_mod_depth[lyr] += soilcolumn.var_depth(depth_prev, depth_this, MODTOT, DEPTH) * to_gridcell_average;

					for (int mth = 0; mth < 12; mth++) {
						double wmass = soilcolumn.var_depth(depth_prev, depth_this, MWMASS, mth) / wsats;
						double imass = soilcolumn.var_depth(depth_prev, depth_this, MIMASS, mth) / wsats;
						mwmass_depth[lyr][mth] += wmass * to_gridcell_average;
						mimass_depth[lyr][mth] += imass * to_gridcell_average;
						mwimass_depth[lyr][mth] += (wmass + imass) * to_gridcell_average;
						mdrain_depth[lyr][mth] += soilcolumn.var_depth(depth_prev, depth_this, MDRAIN, mth) * to_gridcell_average;
					}

					soilc_prev = 0.0;
					soiln_prev = 0.0;
					soilr_prev = 0.0;
				}

				depth_prev = depth_this;
			}

			// 10 cm intervals down to 3m
			double dcm = 0.0;
			double dint = 100.0;
			for (int d = 0; d < 30; d++) {
				double frac = 1.0;
				if (soilcolumn.maxthaw >= dcm && soilcolumn.maxthaw < dcm + dint) {
					frac = (soilcolumn.maxthaw - dcm) / dint;
				} 
				else if (soilcolumn.maxthaw < dcm) {
					frac = 0.0;
				}

				double carbon = soilcolumn.var_depth(dcm, dcm + dint, CCONT);
				soilc_cm[d] += carbon * frac * to_gridcell_average;
				soilcp_cm[d] += carbon * (1.0 - frac) * to_gridcell_average;
				dcm += dint;
			}

			soilc_depth[sslyr + 1] += soilcolumn.get_sum(CCONT) * to_gridcell_average;
			soiln_depth[sslyr + 1] += soilcolumn.get_sum(NCONT) * to_gridcell_average;
			soilr_depth[sslyr + 1] += soilcolumn.get_sum(ARH) * to_gridcell_average;
			soilnetnmin_depth[sslyr] += (soilcolumn.get_sum(ANMIN) - soilcolumn.get_sum(ANIMM)) * to_gridcell_average;
			clitterinput_depth[sslyr + 1] += (surflitinput + soilcolumn.get_sum(ALITTERINPUT)) * to_gridcell_average;
			volume_org_depth[sslyr] += soilcolumn.get_sum(VOL_ORG) * to_gridcell_average;
			volume_min_depth[sslyr] += soilcolumn.get_sum(VOL_MIN) * to_gridcell_average;
			volume_root_depth[sslyr] += soilcolumn.get_sum(VOL_ROOT) * to_gridcell_average;
			double wsats_tot = soilcolumn.get_sum(WSATS);
			wsats_depth[sslyr] += wsats_tot / soilcolumn.depth * to_gridcell_average;
			temp_mod_depth[sslyr] += soilcolumn.var_depth(0.0, soilcolumn.depth, MODTEMP, DEPTH) * to_gridcell_average;
			moist_mod_depth[sslyr] += soilcolumn.var_depth(0.0, soilcolumn.depth, MODMOIS, DEPTH) * to_gridcell_average;
			tot_mod_depth[sslyr] += soilcolumn.var_depth(0.0, soilcolumn.depth, MODTOT, DEPTH) * to_gridcell_average;

			depth_prev = 0.0;
			for (int lyr = 0; lyr < sstlyr; lyr++) {
				double depth_this = soil_temp_sections[lyr];
				double depth_lyr = depth_this - depth_prev;

				for (int mth = 0; mth < 12; mth++) {
					mtemp_depth[lyr][mth] += soilcolumn.var_depth(depth_prev, depth_this, MTEMP, DEPTH, mth) * to_gridcell_average;
				}

				depth_prev = depth_this;
			}

			for (int mth = 0; mth < 12; mth++) {
				double wmass_tot = soilcolumn.get_sum(MWMASS, mth);
				double imass_tot = soilcolumn.get_sum(MIMASS, mth);

				mwmass_depth[sslyr][mth] += wmass_tot / wsats_tot * to_gridcell_average;
				mimass_depth[sslyr][mth] += imass_tot / wsats_tot * to_gridcell_average;
				mwimass_depth[sslyr][mth] += (wmass_tot + imass_tot) / wsats_tot * to_gridcell_average;
				mdrain_depth[sslyr][mth] += soilcolumn.get_sum(MDRAIN, mth) * to_gridcell_average;
			}

			for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {
				Soillayer& soillayer = soilcolumn.soillayers[lyr];
				if (lyr == 0) {
					soilc[lyr] += (soillayer.ccont() - surflitc) * to_gridcell_average;
					soilr[lyr] += (soillayer.get(ARH) - surflitr) * to_gridcell_average;
				}
				else {
					soilc[lyr] += soillayer.ccont() * to_gridcell_average;
					soilr[lyr] += soillayer.get(ARH) * to_gridcell_average;
				}
				temp_mod[lyr] += soillayer.get(MODTEMP) * to_gridcell_average;
				moist_mod[lyr] += soillayer.get(MODMOIS) * to_gridcell_average;
				tot_mod[lyr] += soillayer.get(MODTOT) * to_gridcell_average;
				wsats[lyr] += soillayer.wsats / soillayer.depth * to_gridcell_average;
				wsats[NSOILLAYERS_MAX] += soillayer.wsats / soilcolumn.depth * to_gridcell_average;
			}
			soilc[NSOILLAYERS_MAX] += (soilcolumn.get_sum(CCONT) - surflitc) * to_gridcell_average;
			soilr[NSOILLAYERS_MAX] += (soilcolumn.get_sum(ARH) - soilcolumn.asurflitterresp) * to_gridcell_average;

			soilc_ald += soilcolumn.var_depth(0.0, soilcolumn.maxthaw, CCONT) * to_gridcell_average;
			soilc_tot += soilcolumn.get_sum(CCONT) * to_gridcell_average;
			snowdays += (double)soilcolumn.snow.snowdays * to_gridcell_average;

			// Monthly output variables

			for (m=0;m<12;m++) {

				mpet[m] += fluxes.get_monthly_flux(Fluxes::WPET, m) * to_gridcell_average;
				maet[m] += fluxes.get_monthly_flux(Fluxes::WAET, m) * to_gridcell_average;
				mevap[m] += fluxes.get_monthly_flux(Fluxes::WEVAP, m) * to_gridcell_average;
				mintercep[m] += fluxes.get_monthly_flux(Fluxes::WINTERCEP, m) * to_gridcell_average;
				mrunoff[m] += (fluxes.get_monthly_flux(Fluxes::WRUNOFF_SURF, m) + fluxes.get_monthly_flux(Fluxes::WRUNOFF_DRAIN, m) + fluxes.get_monthly_flux(Fluxes::WRUNOFF_BASE, m)) * to_gridcell_average;
				mrh[m] += fluxes.get_monthly_flux(Fluxes::SOILC, m)*to_gridcell_average;

				mwmass[m] += soilcolumn.get_sum(MWMASS, m) * to_gridcell_average;
				mimass[m] += soilcolumn.get_sum(MIMASS, m) * to_gridcell_average;
				mwcont_upper[m] += soilcolumn.var_depth(0.0, SOILDEPTH_UPPER, MWCONT, AWC, m) * to_gridcell_average;
				mwcont_lower[m] += soilcolumn.var_depth(SOILDEPTH_UPPER, soilcolumn.depth, MWCONT, AWC, m) * to_gridcell_average;

				mgpp[m] += fluxes.get_monthly_flux(Fluxes::GPP, m)*to_gridcell_average;
				mra[m] += fluxes.get_monthly_flux(Fluxes::RA, m)*to_gridcell_average;

				miso[m]+=fluxes.get_monthly_flux(Fluxes::ISO, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_APIN, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_LIMO, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_TRIC, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_BPIN, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_MYRC, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_SABI, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_CAMP, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_TBOC, m)*to_gridcell_average;
				mmon[m]+=fluxes.get_monthly_flux(Fluxes::MT_OTHR, m)*to_gridcell_average;

				mmon_mt1[m]+=fluxes.get_monthly_flux(Fluxes::MT_APIN, m)*to_gridcell_average;
				mmon_mt1[m]+=fluxes.get_monthly_flux(Fluxes::MT_LIMO, m)*to_gridcell_average;
				mmon_mt1[m]+=fluxes.get_monthly_flux(Fluxes::MT_TRIC, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_BPIN, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_MYRC, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_SABI, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_CAMP, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_TBOC, m)*to_gridcell_average;
				mmon_mt2[m]+=fluxes.get_monthly_flux(Fluxes::MT_OTHR, m)*to_gridcell_average;
				
				msnowdepth[m] += patch.soilcolumn.snow.mdepth[m] * M_PER_MM * to_gridcell_average;
				msnowcover[m] += patch.soilcolumn.snow.mcover[m] * to_gridcell_average;
				mwmass_stand[m] += patch.soilcolumn.mwmass_stand[m] * to_gridcell_average;
				mald[m] += min(50.0, (hydrology == GERTEN ? patch.soilcolumn.depth : patch.soilcolumn.mthaw[m]) * M_PER_MM) * to_gridcell_average;
				mwtp[m] += min(50.0, (hydrology == GERTEN ? patch.soilcolumn.depth : patch.soilcolumn.mwtp[m]) * M_PER_MM) * to_gridcell_average;
				mtemp25cm[m] += (hydrology == GERTEN ? patch.soilcolumn.soillayers[0].mtemp[m] : (patch.soilcolumn.nsoillayers >= 3 ? patch.soilcolumn.soillayers[2].mtemp[m] : patch.soilcolumn.pad[0].mtemp[m])) * to_gridcell_average;
				
				mch4[m] += fluxes.get_monthly_flux(Fluxes::CH4C, m)*to_gridcell_average;				// g CH4-C/m2 - as CH4 is in gC, but CO2 fluxes are kgC
				mch4_diff[m] += fluxes.get_monthly_flux(Fluxes::CH4C_DIFF, m)*to_gridcell_average;	// g CH4-C/m2
				mch4_plant[m] += fluxes.get_monthly_flux(Fluxes::CH4C_PLAN, m)*to_gridcell_average;	// g CH4-C/m2
				mch4_ebull[m] += fluxes.get_monthly_flux(Fluxes::CH4C_EBUL, m)*to_gridcell_average;	// g CH4-C/m2
				
				for (int slyr = 0; slyr < NSNOWLAYERS_MAX; slyr++) {
					mtemp[slyr][m] += patch.soilcolumn.snow.mtemp[abs(slyr - (NSNOWLAYERS_MAX - 1))][m] * to_gridcell_average;
				}

				for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
					if (lyr < soilcolumn.nsoillayers) {
						Soillayer& soillayer = soilcolumn.soillayers[lyr];
						mtemp[lyr + NSNOWLAYERS_MAX][m] += soillayer.mtemp[m] * to_gridcell_average;
						mnetmin[m] += soillayer.mnetmin[m] * to_gridcell_average;
					}
					else {
						mtemp[lyr + NSNOWLAYERS_MAX][m] += NOT_A_VALUE_LARGE * to_gridcell_average;
					}
				}

				for (int plyr = 0; plyr < NPADLAYERS; plyr++) {
					mtemp[plyr + NSNOWLAYERS_MAX + NSOILLAYERS_MAX][m] += soilcolumn.pad[plyr].mtemp[m] * to_gridcell_average;
				}
				double wmass_tot = 0.0;
				double imass_tot = 0.0;
				double drain_tot = 0.0;
				for (int lyr = 0; lyr < NSOILLAYERS_MAX+1; lyr++) {
					if (lyr < soilcolumn.nsoillayers) {
						Soillayer& soillayer = soilcolumn.soillayers[lyr];
						mwmasslyr[lyr][m] += max(0.0, soillayer.mwmass[m]) / soillayer.wsats * to_gridcell_average;
						mimasslyr[lyr][m] += max(0.0, soillayer.mimass[m]) / soillayer.wsats * to_gridcell_average;
						wmass_tot += max(0.0, soillayer.mwmass[m]);
						imass_tot += max(0.0, soillayer.mimass[m]);
						mdrainlyr[lyr][m] += soillayer.mdrain[m] * to_gridcell_average;
						drain_tot += soillayer.mdrain[m];
						if (m == 0) {
							wsatslyr[lyr] += soillayer.wsats * to_gridcell_average;
						}
					}
					else if (lyr == NSOILLAYERS_MAX) {
						mwmasslyr[lyr][m] += wmass_tot / soilcolumn.get_sum(WSATS) * to_gridcell_average;
						mimasslyr[lyr][m] += imass_tot / soilcolumn.get_sum(WSATS) * to_gridcell_average;
						mdrainlyr[lyr][m] += drain_tot * to_gridcell_average;
						if (m == 0) {
							wsatslyr[lyr] += soilcolumn.get_sum(WSATS) * to_gridcell_average;
						}
					}
					else {
						mwmasslyr[lyr][m] += NOT_A_VALUE_SMALL * to_gridcell_average;
						mimasslyr[lyr][m] += NOT_A_VALUE_SMALL * to_gridcell_average;
						mdrainlyr[lyr][m] += NOT_A_VALUE_SMALL * to_gridcell_average; 
						if (m == 0) {
							wsatslyr[lyr] += NOT_A_VALUE_SMALL * to_gridcell_average;
						}
					}
				}
			}

			maxald_gridcell += min(50.0, (hydrology == GERTEN ? PAD_DEPTH : patch.soilcolumn.maxthaw) * M_PER_MM) * to_gridcell_average;
			permafrost_full_gridcell += (hydrology == GERTEN || patch.soilcolumn.maxthaw >= PAD_DEPTH ? 0.0 : 1.0) * to_gridcell_average;
			permafrost_3m_gridcell += (hydrology == GERTEN || patch.soilcolumn.maxthaw >= 3000.0 ? 0.0 : 1.0) * to_gridcell_average;
			
			// Calculate monthly NPP and LAI

			Vegetation& vegetation = patch.vegetation;

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				// guess2008 - alive check added
				if (indiv.id != -1 && indiv.alive) {

					for (m=0;m<12;m++) {
						mlai[m] += indiv.mlai[m] * to_gridcell_average;
					}

				} // alive?

				vegetation.nextobj();

			} // while/vegetation loop
			stand.nextobj();
		} // patch loop
		++gc_itr;
	} // stand loop

	// In contrast to annual NEE, monthly NEE does not include fire
	// or establishment fluxes
	for (m=0;m<12;m++) {
		mnpp[m] = mgpp[m] - mra[m];
		mnee[m] = mrh[m] - mnpp[m];
	}

	// Total gridcell soil C
	for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {

		for (int p = 0; p < NSOMPOOL; p++) {
			csom_gridcell[p] += somc[p][lyr];
			nsom_gridcell[p] += somn[p][lyr];
			csom_gridcell[NSOMPOOL] += somc[p][lyr];
			nsom_gridcell[NSOMPOOL] += somn[p][lyr];
		}
	}

	// Print gridcell totals to files

	// Determine total leaf C:N ratio
	double cton_leaf_gridcell = limited_cton(cmass_leaf_gridcell, nmass_leaf_gridcell);

	// Determine total vmax nitrogen limitation
	if (agpp_no_norm_gridcell > 0.0) {
		vmaxnlim_gridcell /= agpp_no_norm_gridcell;
		lambda_gridcell /= agpp_no_norm_gridcell;
	}

	outlimit(out,out_cmass,					cmass_gridcell);
	outlimit(out,out_cmass_agb,				cmass_gridcell); 
	outlimit(out,out_anpp,					anpp_gridcell);
	outlimit(out,out_agpp,					agpp_gridcell);
	outlimit(out,out_agpp_no,				agpp_no_gridcell);
	outlimit(out,out_agpp_ns,				agpp_ns_gridcell);
	outlimit(out,out_fpc,					fpc_gridcell);
	outlimit(out,out_aaet,					aaet_gridcell);
	outlimit(out,out_dens,					dens_gridcell);
	outlimit(out,out_lai,					lai_gridcell);
	outlimit(out,out_clitter,				clitter_gridcell);
	outlimit(out,out_aburned_area,				gridcell.annual_burned_area);
	outlimit(out,out_simfireanalysis,			gridcell.simfire_biome);
	outlimit(out,out_simfireanalysis,			gridcell.nesterov_max);
	outlimit(out,out_simfireanalysis,			gridcell.pop_density);
	outlimit(out,out_simfireanalysis,			gridcell.fapar_annual_max);
	outlimit(out,out_simfireanalysis,			gridcell.simfire_annual_burned_area);
	outlimit(out,out_simfireanalysis,			gridcell.simfire_region);
	outlimit(out,out_firert,				firert_gridcell);
	outlimit(out,out_firert,				burned_area_gridcell);
	outlimit(out,out_runoff,				surfrunoff_gridcell);
	outlimit(out,out_runoff,				drainrunoff_gridcell);
	outlimit(out,out_runoff,				baserunoff_gridcell);
	outlimit(out,out_runoff,				runoff_gridcell);
	outlimit(out,out_wetland_water_added,	wetland_wmass_added_gridcell);

	outlimit(out,out_aiso,		aiso_gridcell);
	outlimit(out,out_amon,		amon_gridcell);
	outlimit(out,out_amon_mt1,  amon_mt1_gridcell);
	outlimit(out,out_amon_mt2,  amon_mt2_gridcell);

	outlimit(out,out_nmass,    (nmass_gridcell + nlitter_gridcell) * M2_PER_HA);
	outlimit(out,out_cton_leaf, cton_leaf_gridcell);
	outlimit(out,out_vmaxnlim,  vmaxnlim_gridcell);
	outlimit(out,out_lambda,    lambda_gridcell);
	outlimit(out,out_nuptake,   nuptake_gridcell * M2_PER_HA);
	outlimit(out,out_nlitter,   nlitter_gridcell * M2_PER_HA);

	outlimit(out,out_nsources, aNH4dep_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, aNO3dep_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, anfix_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, anfert_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, (aNH4dep_gridcell + aNO3dep_gridcell + anfix_gridcell + anfert_gridcell) * M2_PER_HA);
	outlimit(out,out_nsources, anmin_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, animm_gridcell * M2_PER_HA);
	outlimit(out,out_nsources, (anmin_gridcell - animm_gridcell) * M2_PER_HA);
	outlimit(out,out_nsources, (anmin_gridcell - animm_gridcell + aNH4dep_gridcell + aNO3dep_gridcell +
				anfix_gridcell + anfert_gridcell) * M2_PER_HA);

	outlimit(out, out_wflux, -gridcell.climate.aprec);
	if (run_landcover) {
		outlimit(out, out_wflux, -irrigation_gridcell);
		outlimit(out, out_wflux, -wetland_wmass_added_gridcell);
	}
	outlimit(out, out_wflux, intercept_gridcell);
	outlimit(out, out_wflux, aet_gridcell);
	outlimit(out, out_wflux, evap_gridcell);
	outlimit(out, out_wflux, surfrunoff_gridcell);
	outlimit(out, out_wflux, drainrunoff_gridcell);
	outlimit(out, out_wflux, baserunoff_gridcell);
	if (run_landcover) {
		outlimit(out, out_wflux, intercept_gridcell + aet_gridcell + evap_gridcell + 
			surfrunoff_gridcell + drainrunoff_gridcell + baserunoff_gridcell - 
			(gridcell.climate.aprec + irrigation_gridcell + wetland_wmass_added_gridcell));
	}
	else {
		outlimit(out, out_wflux, intercept_gridcell + aet_gridcell + evap_gridcell +
			surfrunoff_gridcell + drainrunoff_gridcell + baserunoff_gridcell - gridcell.climate.aprec);
	}

	outlimit(out, out_wpool, wmass_snow_gridcell);
	outlimit(out, out_wpool, imass_snow_gridcell);
	outlimit(out, out_wpool, wmass_soil_gridcell);
	outlimit(out, out_wpool, imass_soil_gridcell);
	

	// Print landcover totals to files
	if (run_landcover) {
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			if(run[i]) {
				outlimit(out,out_cmass, landcover_cmass[i]);
				outlimit(out,out_cmass_agb, landcover_cmass_agb[i]);
				outlimit(out,out_anpp,  landcover_anpp[i]);
				outlimit(out,out_agpp, landcover_agpp[i]);
				outlimit(out,out_agpp_no, landcover_agpp_no[i]);
				outlimit(out,out_agpp_ns, landcover_agpp_ns[i]);
				outlimit(out,out_fpc,   landcover_fpc[i]);
				outlimit(out,out_aaet,  landcover_aaet[i]);
				outlimit(out,out_dens,  landcover_densindiv_total[i]);
				outlimit(out,out_lai,   landcover_lai[i]);
				outlimit(out,out_clitter, landcover_clitter[i]);
				outlimit(out,out_aiso,  landcover_aiso[i]);
				outlimit(out,out_amon,  landcover_amon[i]);
				outlimit(out,out_amon_mt1,  landcover_amon_mt1[i]);
				outlimit(out,out_amon_mt2,  landcover_amon_mt2[i]);

				double landcover_cton_leaf = limited_cton(landcover_cmass_leaf[i],
						landcover_nmass_leaf[i]);

				if (landcover_agpp_no_norm[i] > 0.0) {
					landcover_vmaxnlim[i] /= landcover_agpp_no_norm[i];
					landcover_lambda[i] /= landcover_agpp_no_norm[i];
				}

				outlimit(out,out_nmass,     (landcover_nmass[i] + landcover_nlitter[i]) * M2_PER_HA);
				outlimit(out,out_cton_leaf, landcover_cton_leaf);
				outlimit(out,out_vmaxnlim,  landcover_vmaxnlim[i]);
				outlimit(out,out_lambda,    landcover_lambda[i]);
				outlimit(out,out_nuptake,   landcover_nuptake[i] * M2_PER_HA);
				outlimit(out,out_nlitter,   landcover_nlitter[i] * M2_PER_HA);
			}
		}
	}

	// Print soil layer variables
	// Total gridcell soil C
	int number_of_layers = hydrology == GERTEN ? NSOILLAYERS_SIMPLE : NSOILLAYERS_MAX;
	double tot_cdec[2] = { 0.0 };
	double tot_cinc[2] = { 0.0 };
	double clitterinput_tot = clitterinput[0];
	outlimit(out, out_clitterinput, clitterinput[0]);
	for (int lyr = 0; lyr < number_of_layers; lyr++) {
		if (date.year == nyear_spinup)
			outlimit(out, out_soildepth, soildepth[lyr] * M_PER_MM);

		tot_cdec[0] += cdec[lyr][0];
		tot_cdec[1] += cdec[lyr][1];
		tot_cinc[0] += cinc[lyr][0];
		tot_cinc[1] += cinc[lyr][1];

		double sdl_decay = 1.0, sdl_growth = 1.0, sdl_total = 1.0, sdl_cue = 1.0;
		if (!negligible(cdec[lyr][1]) && !negligible(cinc[lyr][0]) && !negligible(cdec[lyr][0])) {
			sdl_decay = cdec[lyr][1] / cdec[lyr][0];
			sdl_growth = (cinc[lyr][1] / cdec[lyr][1]) / (cinc[lyr][0] / cdec[lyr][0]);
			sdl_total = (cinc[lyr][1] / cdec[lyr][0]) / (cinc[lyr][0] / cdec[lyr][0]);
			sdl_cue = cinc[lyr][1] / cdec[lyr][0];
		}
		outlimit(out, out_soildeclim, sdl_decay);
		outlimit(out, out_soilgrolim, sdl_growth);
		outlimit(out, out_soiltotlim, sdl_total);
		outlimit(out, out_soilcuelim, sdl_cue);

		outlimit(out, out_clitterinput, clitterinput[lyr + 1]);
		clitterinput_tot += clitterinput[lyr + 1];

		outlimit(out, out_csoil, soilc[lyr]);
		outlimit(out, out_rh, soilr[lyr]);
		outlimit(out, out_tau_soil, !negligible(soilr[lyr]) ? soilc[lyr] / soilr[lyr] : 0.0);
		outlimit(out, out_temp_mod, temp_mod[lyr]);
		outlimit(out, out_moist_mod, moist_mod[lyr]);
		outlimit(out, out_tot_mod, tot_mod[lyr]);
		outlimit(out, out_wsats, wsats[lyr]);
	}

	double sdl_decay_tot = 1.0, sdl_growth_tot = 1.0, sdl_total_tot = 1.0, sdl_cue_tot = 1.0;
	if (!negligible(tot_cdec[1]) && !negligible(tot_cinc[0]) && !negligible(tot_cdec[0])) {
		sdl_decay_tot = tot_cdec[1] / tot_cdec[0];
		sdl_growth_tot = (tot_cinc[1] / tot_cdec[1]) / (tot_cinc[0] / tot_cdec[0]);
		sdl_total_tot = (tot_cinc[1] / tot_cdec[0]) / (tot_cinc[0] / tot_cdec[0]);
		sdl_cue_tot = tot_cinc[1] / tot_cdec[0];
	}

	if (date.year == nyear_spinup)
		outlimit(out, out_soildepth, soildepth[number_of_layers - 1] * M_PER_MM);
	outlimit(out, out_soildeclim, sdl_decay_tot);
	outlimit(out, out_soilgrolim, sdl_growth_tot);
	outlimit(out, out_soiltotlim, sdl_total_tot);
	outlimit(out, out_soilcuelim, sdl_cue_tot);

	outlimit(out, out_clitterinput, clitterinput_tot);

	outlimit(out, out_csoil, soilc[NSOILLAYERS_MAX]);
	outlimit(out, out_rh, soilr[NSOILLAYERS_MAX]);
	outlimit(out, out_tau_soil, !negligible(soilr[NSOILLAYERS_MAX]) ? soilc[NSOILLAYERS_MAX] / soilr[NSOILLAYERS_MAX] : 0.0);
	outlimit(out, out_wsats, wsats[NSOILLAYERS_MAX]);

	// Surface output
	outlimit(out, out_csoil_depth, soilc_depth[0]);
	outlimit(out, out_nsoil_depth, soiln_depth[0]);
	outlimit(out, out_rh_depth, soilr_depth[0]);
	outlimit(out, out_tau_soil_depth, !negligible(soilr_depth[0]) ? soilc_depth[0] / soilr_depth[0] : 0.0);
	outlimit(out, out_clitterinput_depth, clitterinput_depth[0]);

	// Soil depth and total output
	for (int depths = 0; depths < sslyr + 1; depths++) {
		outlimit(out, out_csoil_depth, soilc_depth[depths+1]);
		outlimit(out, out_nsoil_depth, soiln_depth[depths+1]);
		outlimit(out, out_netnminsoil_depth, soilnetnmin_depth[depths] * M2_PER_HA);
		outlimit(out, out_rh_depth, soilr_depth[depths+1]);
		outlimit(out, out_tau_soil_depth, !negligible(soilr_depth[depths+1]) ? soilc_depth[depths+1] / soilr_depth[depths+1] : 0.0);
		outlimit(out, out_clitterinput_depth, clitterinput_depth[depths + 1]);
		outlimit(out, out_volume_org_depth, volume_org_depth[depths] * M_PER_MM);
		outlimit(out, out_volume_min_depth, volume_min_depth[depths] * M_PER_MM);
		outlimit(out, out_volume_root_depth, volume_root_depth[depths] * M_PER_MM);
		outlimit(out, out_wsats_depth, wsats_depth[depths]);
		outlimit(out, out_temp_mod_depth, temp_mod_depth[depths]);
		outlimit(out, out_moist_mod_depth, moist_mod_depth[depths]);
		outlimit(out, out_tot_mod_depth, tot_mod_depth[depths]);
	}

	// Soil depth and total output
	for (int d = 0; d < 30; d++) {
		outlimit(out, out_csoil_cm, soilc_cm[d]);
		outlimit(out, out_csoilp_cm, soilcp_cm[d]);
	}

	outlimit(out, out_soilcperma, soilc_ald);
	outlimit(out, out_soilcperma, soilc_tot - soilc_ald);
	outlimit(out, out_snowdays, snowdays);

	for (int p = 0; p < NSOMPOOL; p++) {
		outlimit(out, out_csom, csom_gridcell[p]);
		outlimit(out, out_csomlitterinput, csomlitterinput[p]);
		outlimit(out, out_ctonsom, !negligible(nsom_gridcell[p]) ? min(999.99, csom_gridcell[p] / nsom_gridcell[p]) : 0.0);
	}
	outlimit(out, out_csom, csom_gridcell[NSOMPOOL]);
	outlimit(out, out_csomlitterinput, csomlitterinput[NSOMPOOL]);
	outlimit(out, out_ctonsom, !negligible(nsom_gridcell[NSOMPOOL]) ? min(999.99, csom_gridcell[NSOMPOOL] / nsom_gridcell[NSOMPOOL]) : 0.0);

	// Print monthly output variables
	for (m=0;m<12;m++) {
		outlimit(out,out_mnpp,         mnpp[m]);
		outlimit(out,out_mlai,         mlai[m]);
		outlimit(out,out_mgpp,         mgpp[m]);
		outlimit(out,out_mra,          mra[m]);
		outlimit(out,out_maet,         maet[m]);
		outlimit(out,out_mpet,         mpet[m]);
		outlimit(out,out_mevap,        mevap[m]);
		outlimit(out,out_mrunoff,      mrunoff[m]);
		outlimit(out,out_mintercep,    mintercep[m]);
		outlimit(out,out_mrh,          mrh[m]);
		outlimit(out,out_mnee,         mnee[m]);
		outlimit(out,out_mnetmin,      mnetmin[m] * M2_PER_HA);
		outlimit(out,out_mwmass,       mwmass[m]);
		outlimit(out,out_mimass,       mimass[m]);
		outlimit(out,out_mwcont_upper, mwcont_upper[m]);
		outlimit(out,out_mwcont_lower, mwcont_lower[m]);
		outlimit(out,out_miso,         miso[m]);
		outlimit(out,out_mmon,         mmon[m]);
		outlimit(out,out_mmon_mt1,     mmon_mt1[m]);
		outlimit(out,out_mmon_mt2,     mmon_mt2[m]);
		outlimit(out,out_mburned_area, (float)gridcell.monthly_burned_area[m]);

		aaet += maet[m];
		apet += mpet[m];
		aevap += mevap[m];
		arunoff += mrunoff[m];
		aintercep += mintercep[m];

		outlimit(out,out_mch4, mch4[m]);
		outlimit(out,out_mch4diff, mch4_diff[m]);
		outlimit(out,out_mch4plan, mch4_plant[m]);
		outlimit(out,out_mch4ebull, mch4_ebull[m]);
		outlimit(out,out_msnow, msnowdepth[m]);
		outlimit(out,out_msnowcover, msnowcover[m]);
		outlimit(out,out_mwmass_stand, mwmass_stand[m]);
		outlimit(out,out_mwtp, mwtp[m]);
		outlimit(out,out_mald, mald[m]);
		outlimit(out,out_msoiltempdepth25,mtemp25cm[m]);
	}

	outlimit(out,out_mald, maxald_gridcell); // [m]
	outlimit(out,out_mald, permafrost_full_gridcell);
	outlimit(out,out_mald, permafrost_3m_gridcell);

	// Graphical output every PLOT_INTERVAL years
	// (Windows shell only - no effect otherwise)

	if (!(date.year%PLOT_INTERVAL)) {
		if(gridcell.nbr_stands() > 0)	//Fixed bug here if no stands were present.
		{
			Stand& stand = gridcell[0];
			plot("C flux [kgC/m2/yr]","veg",  date.year, flux_veg);
			plot("C flux [kgC/m2/yr]","repr", date.year, flux_repr);
			plot("C flux [kgC/m2/yr]","soil", date.year, flux_soil);
			plot("C flux [kgC/m2/yr]","fire", date.year, flux_fire);
			plot("C flux [kgC/m2/yr]","est",  date.year, flux_est);
			plot("C flux [kgC/m2/yr]","NEE",  date.year, flux_veg + flux_repr + flux_soil + flux_fire + flux_est);

			if (!ifcentury) {
				plot("Soil C [kgC/m2]","slow", date.year, stand[0].soilcolumn.get_sum(CSLOW));
				plot("Soil C [kgC/m2]","fast", date.year, stand[0].soilcolumn.get_sum(CFAST));
			}
			else {
				plot("N flux [kgN/ha/yr]","fix",   date.year, -anfix_gridcell * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","dep",   date.year, -(aNH4dep_gridcell - aNO3dep_gridcell) * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","fert",  date.year, -anfert_gridcell * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","leach", date.year, (n_min_leach_gridcell + n_org_leach_gridcell) * M2_PER_HA);
				plot("N flux [kgN/ha/yr]","emissions",  date.year, flux_ntot * M2_PER_HA);

				plot("N flux [kgN/ha/yr]","NEE",   date.year, (flux_ntot + n_min_leach_gridcell + n_org_leach_gridcell -
					(anfix_gridcell + aNH4dep_gridcell + aNO3dep_gridcell + anfert_gridcell)) * M2_PER_HA);

				plot("N mineralization [kgN/ha/yr]","N", date.year, (anmin_gridcell - animm_gridcell) * M2_PER_HA);

				plot("Soil C [kgC/m2]","fine litter",   date.year, surflitterc + soillitterc);
				plot("Soil C [kgC/m2]","coarse litter", date.year, surfcwdc + soilcwdc);
				plot("Soil C [kgC/m2]","soil",          date.year, centuryc + surfhumusc);
				plot("Soil C [kgC/m2]","total",         date.year, surflitterc + soillitterc + surfcwdc + soilcwdc + centuryc + surfhumusc);

				plot("Soil N [kgN/m2]","fine litter",   date.year, surflittern + soillittern + surfhumusn);
				plot("Soil N [kgN/m2]","coarse litter", date.year, surfcwdn + soilcwdn);
				plot("Soil N [kgN/m2]","soil",          date.year, centuryn);
				plot("Soil N [kgN/m2]","total",         date.year, surflittern + soillittern + surfcwdn + soilcwdn + centuryn + surfhumusn);
			}

			plot("H2O flux [mm/yr]", "transp", date.year, aaet);
			plot("H2O flux [mm/yr]", "intercep", date.year, aintercep);
			plot("H2O flux [mm/yr]", "evap", date.year, aevap);
			plot("H2O flux [mm/yr]", "runoff", date.year, arunoff);
			plot("H2O flux [mm/yr]", "PET", date.year, apet);
		}
	}

    // Write fluxes to file

	Landcover& lc = gridcell.landcover;

	outlimit(out,out_cflux, flux_veg);
	outlimit(out,out_cflux, -flux_repr);
	outlimit(out,out_cflux, flux_soil + c_org_leach_gridcell);
	outlimit(out,out_cflux, flux_fire);
	outlimit(out,out_cflux, flux_est);
	if (run_landcover) {
			outlimit(out,out_cflux, flux_seed);
			outlimit(out,out_cflux, flux_charvest);
			outlimit(out,out_cflux, lc.acflux_landuse_change);
			outlimit(out,out_cflux, lc.acflux_harvest_slow);
	}
	outlimit(out,out_cflux, flux_veg - flux_repr + flux_soil + flux_fire + flux_est + c_org_leach_gridcell +
			flux_seed + flux_charvest + lc.acflux_landuse_change + lc.acflux_harvest_slow);

	outlimit(out,out_doc, (c_org_leach_gridcell) * M2_PER_HA);

	if (run_landcover) {
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			if(run[i]) {
				outlimit(out,out_doc, c_org_leach_lc[i] * M2_PER_HA);
			}
		}
	}


	outlimit(out,out_nflux, -aNH4dep_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -aNO3dep_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -anfix_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, -anfert_gridcell * M2_PER_HA);
	outlimit(out,out_nflux, flux_ntot * M2_PER_HA);
	outlimit(out,out_nflux, (n_min_leach_gridcell + n_org_leach_gridcell) * M2_PER_HA);

	if (run_landcover) {
			outlimit(out,out_nflux, flux_nseed * M2_PER_HA);
			outlimit(out,out_nflux, flux_nharvest * M2_PER_HA);
			outlimit(out,out_nflux, lc.anflux_landuse_change * M2_PER_HA);
			outlimit(out,out_nflux, lc.anflux_harvest_slow * M2_PER_HA);
	}
	outlimit(out,out_nflux, (flux_nharvest + lc.anflux_landuse_change +
				lc.anflux_harvest_slow + flux_nseed + flux_ntot +
				n_min_leach_gridcell + n_org_leach_gridcell -
					 (aNH4dep_gridcell + aNO3dep_gridcell + anfix_gridcell + anfert_gridcell)) * M2_PER_HA);

	// CPOOL Write cpool to file

	outlimit(out,out_cpool, cmass_gridcell);
	if (!ifcentury) {
		outlimit(out,out_cpool, clitter_gridcell);
		outlimit(out,out_cpool, c_fast);
		outlimit(out,out_cpool, c_slow);
	}
	else {
		outlimit(out,out_cpool, clitter_gridcell + surflitterc + surfcwdc);
		outlimit(out,out_cpool, centuryc + surfhumusc + soillitterc + soilcwdc);
	}

	if (run_landcover && ifslowharvestpool) {
		outlimit(out,out_cpool, c_harv_slow);
	}

	// Calculate total cpool, starting with cmass and litter...
	double cpool_total = cmass_gridcell + clitter_gridcell;

	// Add SOM pools
	if (!ifcentury) {
		cpool_total += c_fast + c_slow;
	}
	else {
		cpool_total += centuryc + surflitterc + soillitterc + surfcwdc + soilcwdc + surfhumusc;
	}

	// Add slow harvest pool if needed
	if (run_landcover && ifslowharvestpool) {
		cpool_total += c_harv_slow;
	}

	outlimit(out,out_cpool, cpool_total);

	// NPOOL Write npool to file
	
	if (ifcentury) {
		outlimit(out,out_npool, nmass_gridcell + nlitter_gridcell);
		outlimit(out,out_npool, surflittern + surfcwdn + surfhumusn);
		outlimit(out,out_npool, centuryn + soillittern + soilcwdn);
		outlimit(out,out_npool, availn);

		if(run_landcover && ifslowharvestpool) {
			outlimit(out,out_npool, n_harv_slow);
			outlimit(out,out_npool, (nmass_gridcell + nlitter_gridcell + surflittern + surfhumusn + surfcwdn + soillittern + soilcwdn + centuryn + availn + n_harv_slow));
		}
		else {
			outlimit(out,out_npool, (nmass_gridcell + nlitter_gridcell + surflittern + surfhumusn + surfcwdn + soillittern + soilcwdn + centuryn + availn));
		}
	}

	outlimit(out,out_ngases, flux_NH3_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NH3_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NOx_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_NOx_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2O_fire   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2O_soil   * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2_fire    * M2_PER_HA);
	outlimit(out,out_ngases, flux_N2_soil    * M2_PER_HA);
	outlimit(out,out_ngases, flux_ntot  * M2_PER_HA);
	// Volatile nitrogen from soil
	outlimit(out,out_soil_nflux, NH3_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, NO_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, N2O_mass_inc * M2_PER_HA);
	outlimit(out,out_soil_nflux, N2_mass_inc * M2_PER_HA);
	// Nitrification/denitrification fluxes
	//outlimit(out,out_soil_nflux, gross_nitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, net_nitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, gross_denitrif * M2_PER_HA);
	//outlimit(out,out_soil_nflux, net_denitrif * M2_PER_HA);
		
	// Ackumulated nitrogen in pools		
	outlimit(out,out_soil_npool, NH4_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO3_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO2_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, NO_mass  * M2_PER_HA);
	outlimit(out,out_soil_npool, N2O_mass * M2_PER_HA);
	outlimit(out,out_soil_npool, N2_mass  * M2_PER_HA);


	// Output of tree stand age structure, monthly soil water and 3D vegetation view
	// (Windows shell only - no effect otherwise)

	if (vegmode==COHORT || vegmode==INDIVIDUAL) {

		if (!(date.year%PLOT_UPDATE_INTERVAL)) {

			double* densindiv=NULL;
			int nageclass;
			get_stand_age_structure(gridcell, densindiv, nageclass, false);

			if (nageclass) {

				densindiv = new double[(npft+1)*nageclass];
				if (densindiv) {

					resetwindow("Age structure [indiv/ha]");

					get_stand_age_structure(gridcell, densindiv, nageclass, true);

					pftlist.firstobj();
					while (pftlist.isobj) {
						Pft& pft = pftlist.getobj();

						if (pft.lifeform == TREE) {

							for (c = 0; c<nageclass; c++)
								plot("Age structure [indiv/ha]", pft.name,
								c * estinterval + estinterval*0.5,
								densindiv[pft.id*nageclass+c]*1e4); // includes conversion from /m2 --> /ha
						}

						pftlist.nextobj();
					}

					delete[] densindiv;
				}
			}
		}
	}

	if (!(date.year % PLOT_UPDATE_INTERVAL)) {

		resetwindow("Soil water [%AWC]");

		for (m = 0; m < 12; m++) {
			plot("Soil water [%AWC]", "upper", m, mwcont_upper[m]);
			plot("Soil water [%AWC]", "lower", m, mwcont_lower[m]);
		}
	}

	if (vegmode == COHORT) {
		if (!(date.year%VEG3D_UPDATE_INTERVAL)) {
			if (date.year == 0) open3d();
			output_vegetation(gridcell, pftlist);
		}
	}


	// Monthly row output
	for (int mth = 0; mth < 12; mth++) {
		OutputRows out(output_channel, lon, lat, date.get_calendar_year(), CTMONTH, mth + 1);

		outlimit(out, out_mtemp, gridcell.climate.hmtemp_20[mth].lastadd());
		if (hydrology == GERTEN) {
			for (int lyr = NSNOWLAYERS_MAX; lyr < NSNOWLAYERS_MAX + NSOILLAYERS_SIMPLE; lyr++) {
				outlimit(out, out_mtemp, mtemp[lyr][mth]);
			}
		}
		else {
			for (int lyr = 0; lyr < NSNOWLAYERS_MAX + NSOILLAYERS_MAX + NPADLAYERS; lyr++) {
				outlimit(out, out_mtemp, mtemp[lyr][mth]);
			}
		}

		outlimit(out, out_mtemp_depth, gridcell.climate.hmtemp_20[mth].lastadd());

		for (int lyr = 0; lyr < number_of_layers; lyr++) {
			outlimit(out, out_mwmasslyr, mwmasslyr[lyr][mth]);
			outlimit(out, out_mimasslyr, mimasslyr[lyr][mth]);
			outlimit(out, out_mwimasslyr, mwmasslyr[lyr][mth] >= 0.0 ? mwmasslyr[lyr][mth] + mimasslyr[lyr][mth] : NOT_A_VALUE_SMALL);
			outlimit(out, out_mdrainlyr, mdrainlyr[lyr][mth]);
			if (mth == 0) {
				outlimit(out, out_wsatslyr, wsatslyr[lyr]);
			}
		}
		outlimit(out, out_mwmasslyr, mwmasslyr[NSOILLAYERS_MAX][mth]);
		outlimit(out, out_mimasslyr, mimasslyr[NSOILLAYERS_MAX][mth]);
		outlimit(out, out_mwimasslyr, mwmasslyr[NSOILLAYERS_MAX][mth] + mimasslyr[NSOILLAYERS_MAX][mth]);
		outlimit(out, out_mdrainlyr, mdrainlyr[NSOILLAYERS_MAX][mth]);
		if (mth == 0) {
			outlimit(out, out_wsatslyr, wsatslyr[NSOILLAYERS_MAX]);
		}

		for (int depths = 0; depths < sslyr + 1; depths++) {
			outlimit(out, out_mwmass_depth, mwmass_depth[depths][mth]);
			outlimit(out, out_mimass_depth, mimass_depth[depths][mth]);
			outlimit(out, out_mwimass_depth, mwimass_depth[depths][mth]);
			outlimit(out, out_mdrain_depth, mdrain_depth[depths][mth]);
		}

		for (int depths = 0; depths < sstlyr; depths++) {
			outlimit(out, out_mtemp_depth, mtemp_depth[depths][mth]);
		}

		// Loop through PFTs
		pftlist.firstobj();
		while (pftlist.isobj) {

			Pft& pft = pftlist.getobj();

			outlimit(out, out_mcton_leaf, min(pft.cton_leaf_max, mcton_leaf[pft.id][mth]));
			outlimit(out, out_mvmaxnlim, mvmaxnlim[pft.id][mth]);
			outlimit(out, out_mlambda, mlambda[pft.id][mth]);

			pftlist.nextobj();

		} // End of PFT loop

	}
}

/// Output of simulation results at the end of each day
/** This function does not have to provide any information to the framework.
  */
void CommonOutput::outdaily(Gridcell& gridcell) {
}

} // namespace
