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

namespace GuessOutput {

REGISTER_OUTPUT_MODULE("common", CommonOutput)

CommonOutput::CommonOutput() {
	// Annual output variables
	declare_parameter("file_cmass", &file_cmass, 300, "C biomass output file");
	declare_parameter("file_anpp", &file_anpp, 300, "Annual NPP output file");
	declare_parameter("file_lai", &file_lai, 300, "LAI output file");
	declare_parameter("file_yield",&file_yield,300, "Crop yield output file");
	declare_parameter("file_yield1",&file_yield1,300,"Crop first yield output file");
	declare_parameter("file_yield2",&file_yield2,300,"Crop second yield output file");
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
	declare_parameter("file_cflux", &file_cflux, 300, "C fluxes output file");
	declare_parameter("file_cflux_cropland", &file_cflux_cropland, 300, "C fluxes output file");
	declare_parameter("file_cflux_pasture", &file_cflux_pasture, 300, "C fluxes output file");
	declare_parameter("file_cflux_natural", &file_cflux_natural, 300, "C fluxes output file");
	declare_parameter("file_cflux_forest", &file_cflux_forest, 300, "C fluxes output file");
	declare_parameter("file_dens", &file_dens, 300, "Tree density output file");
	declare_parameter("file_cpool", &file_cpool, 300, "Soil C output file");
	declare_parameter("file_cpool_cropland", &file_cpool_cropland, 300, "Soil C output file");
	declare_parameter("file_cpool_pasture", &file_cpool_pasture, 300, "Soil C output file");
	declare_parameter("file_cpool_natural", &file_cpool_natural, 300, "Soil C output file");
	declare_parameter("file_cpool_forest", &file_cpool_forest, 300, "Soil C output file");
	declare_parameter("file_runoff", &file_runoff, 300, "Runoff output file");
	declare_parameter("file_firert", &file_firert, 300, "Fire retrun time output file");

	declare_parameter("file_cton_leaf", &file_cton_leaf, 300, "Mean leaf C:N output file");
	declare_parameter("file_cton_veg", &file_cton_veg, 300, "Mean vegetation C:N output file");
	declare_parameter("file_nsources", &file_nsources, 300, "Annual nitrogen sources output file");
	declare_parameter("file_npool", &file_npool, 300, "Soil nitrogen output file");
	declare_parameter("file_nuptake", &file_nuptake, 300, "Annual nitrogen uptake output file");
	declare_parameter("file_vmaxnlim", &file_vmaxnlim, 300, "Annual nitrogen limitation on vm output file");
	declare_parameter("file_nflux", &file_nflux, 300, "Annual nitrogen fluxes output file");
	declare_parameter("file_ngases", &file_ngases, 300, "Annual nitrogen gases output file");

	declare_parameter("file_speciesheights", &file_speciesheights, 300, "Mean species heights");

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
	declare_parameter("file_mwcont_upper", &file_mwcont_upper, 300, "Monthly wcont_upper output file");
	declare_parameter("file_mwcont_lower", &file_mwcont_lower, 300, "Monthly wcont_lower output file");
	// bvoc
	declare_parameter("file_aiso", &file_aiso, 300, "annual isoprene flux output file");
	declare_parameter("file_miso", &file_miso, 300, "monthly isoprene flux output file");
	declare_parameter("file_amon", &file_amon, 300, "annual monoterpene flux output file");
	declare_parameter("file_mmon", &file_mmon, 300, "monthly monoterpene flux output file");
}


CommonOutput::~CommonOutput() {
}

void CommonOutput::init() {
	// Define all output tables and their formats
	define_output_tables();
}

/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void CommonOutput::define_output_tables() {
	// create a vector with the pft names
	std::vector<std::string> pfts;

	// create a vector with the crop pft names
	std::vector<std::string> crop_pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		 Pft& pft=pftlist.getobj();

		 pfts.push_back((char*)pft.name);

		 if(pft.landcover==CROPLAND)
			 crop_pfts.push_back((char*)pft.name);

		 pftlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if (run_landcover) {
		 const char* landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum", "Forest_sum", "Natural_sum", "Peatland_sum"};
		 for (int i=0; i<NLANDCOVERTYPES; i++) {
			  if(run[i]) {
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
	cmass_columns += ColumnDescriptor("Total",             8, 3);
	cmass_columns += ColumnDescriptors(landcovers,        13, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;

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
	cflux_columns += ColumnDescriptor("Soil",              8, 3);
	cflux_columns += ColumnDescriptor("Fire",              8, 3);
	cflux_columns += ColumnDescriptor("Est",               8, 3);
	if (run_landcover) {
		 cflux_columns += ColumnDescriptor("Seed",         8, 3);
		 cflux_columns += ColumnDescriptor("Harvest",      9, 3);
		 cflux_columns += ColumnDescriptor("LU_ch",        9, 3);
		 cflux_columns += ColumnDescriptor("Slow_h",       9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",              10, 5);

	ColumnDescriptors cflux2_columns;
	cflux2_columns += ColumnDescriptor("Veg",               8, 3);
	cflux2_columns += ColumnDescriptor("Soil",              8, 3);
	cflux2_columns += ColumnDescriptor("Fire",              8, 3);
	cflux2_columns += ColumnDescriptor("Est",               8, 3);
	if (run_landcover) {
		 cflux2_columns += ColumnDescriptor("Seed",         8, 3);
		 cflux2_columns += ColumnDescriptor("Harvest",      9, 3);
	}
	cflux2_columns += ColumnDescriptor("NEE",              10, 5);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",              8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",         8, 3);
		cpool_columns += ColumnDescriptor("SoilfC",        8, 3);
		cpool_columns += ColumnDescriptor("SoilsC",        8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LittVC",        8, 3);
		cpool_columns += ColumnDescriptor("LittSC",        8, 3);
		cpool_columns += ColumnDescriptor("CwdC",          8, 3);
		cpool_columns += ColumnDescriptor("SoilC",         8, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC",   10, 3);
	}
	cpool_columns += ColumnDescriptor("Total",            10, 3);

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
	seasonality_columns += ColumnDescriptor("temp_seas",  10, 0);
	seasonality_columns += ColumnDescriptor("prec_min",   10, 2);
	seasonality_columns += ColumnDescriptor("prec",       10, 1);
	seasonality_columns += ColumnDescriptor("prec_range", 12, 0);
//	seasonality_columns += ColumnDescriptor("biseasonal", 12, 0);

	// FIRERT
	ColumnDescriptors firert_columns;
	firert_columns += ColumnDescriptor("FireRT",           8, 1);

	// RUNOFF
	ColumnDescriptors runoff_columns;
	runoff_columns += ColumnDescriptor("Surf",             8, 1);
	runoff_columns += ColumnDescriptor("Drain",            8, 1);
	runoff_columns += ColumnDescriptor("Base",             8, 1);
	runoff_columns += ColumnDescriptor("Total",            8, 1);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts,      8, 2);

	// AISO
	ColumnDescriptors aiso_columns;
	aiso_columns += ColumnDescriptors(pfts,               10, 3);
	aiso_columns += ColumnDescriptor("Total",             10, 3);
	aiso_columns += ColumnDescriptors(landcovers,         13, 3);

	// AMON
	ColumnDescriptors amon_columns = aiso_columns;

	//TODO Fix these for landcover

	// CTON
	ColumnDescriptors cton_columns;
	cton_columns += ColumnDescriptors(pfts,                8, 1);
	cton_columns += ColumnDescriptor("Total",              8, 1);
	cton_columns += ColumnDescriptors(landcovers,         12, 1);

	// NSOURCES
	ColumnDescriptors nsources_columns;
	nsources_columns += ColumnDescriptor("dep",            8, 2);
	nsources_columns += ColumnDescriptor("fix",            8, 2);
	nsources_columns += ColumnDescriptor("fert",           8, 2);
	nsources_columns += ColumnDescriptor("input",          8, 2);
	nsources_columns += ColumnDescriptor("min",            7, 2);
	nsources_columns += ColumnDescriptor("imm",            7, 2);
	nsources_columns += ColumnDescriptor("netmin",         7, 2);
	nsources_columns += ColumnDescriptor("Total",          7, 2);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",              9, 5);
	npool_columns += ColumnDescriptor("LittVN",            9, 5);
	npool_columns += ColumnDescriptor("LittSN",            9, 5);
	npool_columns += ColumnDescriptor("CwdN",              9, 5);
	npool_columns += ColumnDescriptor("SoilN",             9, 5);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN",     9, 5);
	}

	npool_columns += ColumnDescriptor("Total",            10, 5);

	// NUPTAKE
	ColumnDescriptors nuptake_columns;
	nuptake_columns += ColumnDescriptors(pfts,             7, 2);
	nuptake_columns += ColumnDescriptor("Total",           7, 2);
	nuptake_columns += ColumnDescriptors(landcovers,      10, 2);

	// VMAXNLIM
	ColumnDescriptors vmaxnlim_columns;
	vmaxnlim_columns += ColumnDescriptors(pfts,            8, 2);
	vmaxnlim_columns += ColumnDescriptor("Total",          8, 2);
	vmaxnlim_columns += ColumnDescriptors(landcovers,     13, 2);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("dep",               8, 2);
	nflux_columns += ColumnDescriptor("fix",               8, 2);
	nflux_columns += ColumnDescriptor("fert",              8, 2);
	nflux_columns += ColumnDescriptor("flux",              8, 2);
	nflux_columns += ColumnDescriptor("leach",             8, 2);
	if (run_landcover) {
		nflux_columns += ColumnDescriptor("harvest",       8, 2);
	}
	nflux_columns += ColumnDescriptor("NEE",               8, 2);

	// NGASES
	ColumnDescriptors ngases_columns;
	ngases_columns += ColumnDescriptor("NH3",              9, 3);
	ngases_columns += ColumnDescriptor("NO",               9, 3);
	ngases_columns += ColumnDescriptor("NO2",              9, 3);
	ngases_columns += ColumnDescriptor("N2O",              9, 3);
	ngases_columns += ColumnDescriptor("NSoil",            9, 3);
	ngases_columns += ColumnDescriptor("Total",            9, 3);

	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_cmass,          file_cmass,          cmass_columns);
	create_output_table(out_anpp,           file_anpp,           anpp_columns);
	create_output_table(out_dens,           file_dens,           dens_columns);
	create_output_table(out_lai,            file_lai,            lai_columns);
	create_output_table(out_cflux,          file_cflux,          cflux_columns);
	create_output_table(out_cflux_cropland, file_cflux_cropland, cflux_columns);
	create_output_table(out_cflux_pasture,  file_cflux_pasture,  cflux_columns);
	create_output_table(out_cflux_natural,  file_cflux_natural,  cflux_columns);
	create_output_table(out_cflux_forest,	file_cflux_forest,	 cflux_columns);	
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
	create_output_table(out_cpool_cropland, file_cpool_cropland, cpool_columns);
	create_output_table(out_cpool_pasture,  file_cpool_pasture,  cpool_columns);
	create_output_table(out_cpool_natural,  file_cpool_natural,  cpool_columns);
	create_output_table(out_cpool_forest,	file_cpool_forest,	 cpool_columns);

	if(run_landcover && run[CROPLAND])
	{
		create_output_table(out_yield,      file_yield,          crop_columns);
		create_output_table(out_yield1,     file_yield1,         crop_columns);
		create_output_table(out_yield2,     file_yield2,         crop_columns);
		create_output_table(out_sdate1,     file_sdate1,         date_columns);
		create_output_table(out_sdate2,     file_sdate2,         date_columns);
		create_output_table(out_hdate1,     file_hdate1,         date_columns);
		create_output_table(out_hdate2,     file_hdate2,         date_columns);
		create_output_table(out_lgp,	    file_lgp,	         date_columns);
		create_output_table(out_phu,	    file_phu,	         date_columns);
		create_output_table(out_fphu,	    file_fphu,	         crop_columns);
		create_output_table(out_fhi,	    file_fhi,	         crop_columns);
		create_output_table(out_irrigation, file_irrigation,     irrigation_columns);
		create_output_table(out_seasonality,file_seasonality,    seasonality_columns);	
	}
	create_output_table(out_firert,         file_firert,         firert_columns);
	create_output_table(out_runoff,         file_runoff,         runoff_columns);
	create_output_table(out_speciesheights, file_speciesheights, speciesheights_columns);
	create_output_table(out_aiso,           file_aiso,           aiso_columns);
	create_output_table(out_amon,           file_amon,           amon_columns);

	create_output_table(out_cton_leaf,      file_cton_leaf,      cton_columns);
	create_output_table(out_cton_veg,       file_cton_veg,       cton_columns);
	create_output_table(out_nsources,       file_nsources,       nsources_columns);
	create_output_table(out_npool,          file_npool,          npool_columns);
	create_output_table(out_nuptake,        file_nuptake,        nuptake_columns);
	create_output_table(out_vmaxnlim,       file_vmaxnlim,       vmaxnlim_columns);
	create_output_table(out_nflux,          file_nflux,          nflux_columns);
	create_output_table(out_ngases,         file_ngases,         ngases_columns);

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
	create_output_table(out_mwcont_upper,   file_mwcont_upper,   month_columns);
	create_output_table(out_mwcont_lower,   file_mwcont_lower,   month_columns);
	create_output_table(out_miso,           file_miso,           month_columns_wide);
	create_output_table(out_mmon,           file_mmon,           month_columns_wide);
}


/// Help function to prepare C:N values for output
/** Avoids division by zero and limits the results to a maximum
 *  value to avoid inf or values large enough to ruin the alignment
 *  in the output.
 *
 *  If both cmass and nmass is 0, the function returns 0.
 */ 
double limited_cton(double cmass, double nmass) {
	const double MAX_CTON = 1000;

	if (nmass > 0.0) {
		return min(MAX_CTON, cmass / nmass);
	}
	else if (cmass > 0.0) {
		return MAX_CTON;
	}
	else {
		return 0.0;
	}
}


void CommonOutput::outannual(Gridcell& gridcell) {
	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int c, m, nclass;
	double flux_veg, flux_soil, flux_fire, flux_est, flux_seed, flux_charvest;
	double c_litter, c_fast, c_slow, c_harv_slow; 

	double surfsoillitterc,surfsoillittern,cwdc,cwdn,centuryc,centuryn,n_litter,n_harv_slow,availn;
	double flux_nh3,flux_no,flux_no2,flux_n2o,flux_nsoil,flux_ntot,flux_nharvest;

	// Nitrogen output is in kgN/ha instead of kgC/m2 as for carbon 
	double m2toha = 10000.0;

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
	double mwcont_upper[12];
	double mwcont_lower[12];
	// bvoc
	double miso[12];
	double mmon[12];

	if (vegmode == COHORT)
		nclass = min(date.year / estinterval + 1, OUTPUT_MAXAGECLASS);
	
	// yearly output after spinup
		
	// If only yearly output between, say 1961 and 1990 is requred, use: 
	//  if (date.get_calendar_year() >= 1961 && date.get_calendar_year() <= 1990) {
	//  (assuming the input module has set the first calendar year in the date object)

	if (date.year >= nyear_spinup) {

		double lon = gridcell.get_lon();
		double lat = gridcell.get_lat();

		// The OutputRows object manages the next row of output for each
		// output table
		OutputRows out(output_channel, lon, lat, date.get_calendar_year());

		// guess2008 - reset monthly average across patches each year
		for (m=0;m<12;m++)
			mnpp[m]=mlai[m]=mgpp[m]=mra[m]=maet[m]=mpet[m]=mevap[m]=mintercep[m]=mrunoff[m]=mrh[m]=mnee[m]=mwcont_upper[m]=mwcont_lower[m]=miso[m]=mmon[m]=0.0;



		double landcover_cmass[NLANDCOVERTYPES]={0.0};
		double landcover_cmass_leaf[NLANDCOVERTYPES]={0.0};
		double landcover_nmass_leaf[NLANDCOVERTYPES]={0.0};
		double landcover_cmass_veg[NLANDCOVERTYPES]={0.0};
		double landcover_nmass_veg[NLANDCOVERTYPES]={0.0};
		double landcover_anpp[NLANDCOVERTYPES]={0.0};
		double landcover_lai[NLANDCOVERTYPES]={0.0};
		double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};
		double landcover_aiso[NLANDCOVERTYPES]={0.0};
		double landcover_amon[NLANDCOVERTYPES]={0.0};
		double landcover_nuptake[NLANDCOVERTYPES]={0.0};
		double landcover_vmaxnlim[NLANDCOVERTYPES]={0.0};

		double stand_mean_cmass=0.0;
		double stand_mean_nmass=0.0;
		double stand_mean_cmass_leaf=0.0;
		double stand_mean_nmass_leaf=0.0;
		double stand_mean_cmass_veg=0.0;
		double stand_mean_nmass_veg=0.0;
		double stand_mean_anpp=0.0;
		double stand_mean_lai=0.0;
		double stand_mean_yield=0.0;
		double stand_mean_yield1=0.0;
		double stand_mean_yield2=0.0;
		double stand_mean_densindiv_total=0.0;
		double stand_mean_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double stand_mean_aiso=0.0;
		double stand_mean_amon=0.0;
		double stand_mean_nuptake=0.0;
		double stand_mean_vmaxnlim=0.0;

		double cmass_gridcell=0.0;
		double nmass_gridcell= 0.0;
		double cmass_leaf_gridcell=0.0;
		double nmass_leaf_gridcell=0.0;
		double cmass_veg_gridcell=0.0;
		double nmass_veg_gridcell=0.0;
		double anpp_gridcell=0.0;
		double lai_gridcell=0.0;
		double surfrunoff_gridcell=0.0;
		double drainrunoff_gridcell=0.0;
		double baserunoff_gridcell=0.0;
		double runoff_gridcell=0.0;
		double dens_gridcell=0.0;
		double firert_gridcell=0.0;
		double aiso_gridcell=0.0;
		double amon_gridcell=0.0;
		double nuptake_gridcell=0.0;
		double vmaxnlim_gridcell=0.0;

		double andep_gridcell=0.0;
		double anfert_gridcell=0.0;
		double anmin_gridcell=0.0;
		double animm_gridcell=0.0;
		double anfix_gridcell=0.0;
		double n_min_leach_gridcell=0.0;
		double n_org_leach_gridcell=0.0;

		double irrigation_gridcell=0.0;

		double standpft_cmass=0.0;
		double standpft_nmass=0.0;
		double standpft_cmass_leaf=0.0;
		double standpft_nmass_leaf=0.0;
		double standpft_cmass_veg=0.0;
		double standpft_nmass_veg=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_yield=0.0;	
		double standpft_yield1=0.0;	
		double standpft_yield2=0.0;	
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double standpft_aiso=0.0;
		double standpft_amon=0.0;
		double standpft_nuptake=0.0;
		double standpft_vmaxnlim=0.0;
		double standpft_ic_nuptake=0.0;
		double standpft_ic_vmaxnlim=0.0;


		// *** Loop through PFTs ***

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Gridcellpft& gridcellpft=gridcell.pft[pft.id];

			// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs		
			stand_mean_cmass=0.0;
			stand_mean_nmass=0.0;
			stand_mean_cmass_leaf=0.0;
			stand_mean_nmass_leaf=0.0;
			stand_mean_cmass_veg=0.0;
			stand_mean_nmass_veg=0.0;
			stand_mean_anpp=0.0;
			stand_mean_lai=0.0;
			stand_mean_yield=0.0;
			stand_mean_yield1=0.0;
			stand_mean_yield2=0.0;
			stand_mean_densindiv_total=0.0;		
			stand_mean_aiso=0.0;
			stand_mean_amon=0.0;
			stand_mean_nuptake=0.0;
			stand_mean_vmaxnlim=0.0;

			double heightindiv_total = 0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				if(standpft.active) {
					// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
					standpft_cmass=0.0;
					standpft_nmass=0.0;
					standpft_cmass_leaf=0.0;
					standpft_nmass_leaf=0.0;
					standpft_cmass_veg=0.0;
					standpft_nmass_veg=0.0;
					standpft_anpp=0.0;
					standpft_lai=0.0;
					standpft_yield=0.0;
					standpft_yield1=0.0;
					standpft_yield2=0.0;
					standpft_densindiv_total = 0.0;
					standpft_aiso=0.0;
					standpft_amon=0.0;
					standpft_nuptake=0.0;
					standpft_vmaxnlim=0.0;
					standpft_ic_nuptake=0.0;
					standpft_ic_vmaxnlim=0.0;

					// Initialise age structure array

					if (vegmode==COHORT || vegmode==INDIVIDUAL)
						for (c=0;c<nclass;c++){
							standpft_densindiv_ageclass[c] = 0.0;
						}			
					stand.firstobj();

					// Loop through Patches
					while (stand.isobj) {
						Patch& patch = stand.getobj();
						Vegetation& vegetation=patch.vegetation;

						standpft_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);
						standpft_aiso += patch.fluxes.get_annual_flux(Fluxes::ISO, pft.id);
						standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MON, pft.id);

						vegetation.firstobj();
						while (vegetation.isobj) {
							Individual& indiv=vegetation.getobj();
								
							// guess2008 - alive check added
							if (indiv.id!=-1 && indiv.alive) { 
								
								if (indiv.pft.id==pft.id) {

									if(indiv.has_daily_turnover() && indiv.cropindiv)
										standpft_cmass += indiv.cropindiv->grs_cmass_leaf + indiv.cropindiv->grs_cmass_root;
									else
										standpft_cmass += indiv.cmass_leaf + indiv.cmass_root +  indiv.cmass_wood();
									standpft_nmass += indiv.nmass_leaf + indiv.nmass_root + 
													  indiv.nmass_wood() + indiv.nstore();
									standpft_cmass_leaf += indiv.cmass_leaf;
									standpft_nmass_leaf += indiv.cmass_leaf / indiv.cton_leaf_aavr;
									standpft_cmass_veg += indiv.cmass_veg;
									standpft_nmass_veg += indiv.nmass_veg;
									standpft_lai+=indiv.lai;
									standpft_vmaxnlim += indiv.avmaxnlim * indiv.cmass_leaf;
									standpft_nuptake += indiv.anuptake;

									if(pft.landcover == CROPLAND) {

										if(indiv.has_daily_turnover())
											standpft_cmass += indiv.cropindiv->grs_cmass_ho + indiv.cropindiv->grs_cmass_agpool;
										else
											standpft_cmass += indiv.cropindiv->cmass_ho + indiv.cropindiv->cmass_agpool;
										standpft_nmass += indiv.cropindiv->nmass_ho + indiv.cropindiv->nmass_agpool;

										standpft_yield += indiv.cropindiv->harv_yield;
										standpft_yield1 += indiv.cropindiv->yield_harvest[0];
										standpft_yield2 += indiv.cropindiv->yield_harvest[1];
									}
									else {
										if (vegmode==COHORT || vegmode==INDIVIDUAL) {
											
											// Age structure
											
											c=(int)(indiv.age/estinterval); // guess2008
											if (c<OUTPUT_MAXAGECLASS)
												standpft_densindiv_ageclass[c]+=indiv.densindiv;

											// guess2008 - only count trees with a trunk above a certain diameter  
											if (pft.lifeform==TREE && indiv.age>0) {
												double diam=pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3);
												if (diam>0.03) {
													standpft_densindiv_total+=indiv.densindiv; // indiv/m2

													heightindiv_total+=indiv.height * indiv.densindiv;
												}
											}
										}
									}						
								}

							} // alive?
							vegetation.nextobj();
						}

						stand.nextobj();
					} // end of patch loop

					standpft_cmass/=(double)stand.npatch();
					standpft_nmass/=(double)stand.npatch();
					standpft_cmass_leaf/=(double)stand.npatch();
					standpft_nmass_leaf/=(double)stand.npatch();
					standpft_cmass_veg/=(double)stand.npatch();
					standpft_nmass_veg/=(double)stand.npatch();
					standpft_anpp/=(double)stand.npatch();
					standpft_lai/=(double)stand.npatch();
					standpft_densindiv_total/=(double)stand.npatch();
					standpft_aiso/=(double)stand.npatch();
					standpft_amon/=(double)stand.npatch();
					standpft_nuptake/=(double)stand.npatch();
					standpft_vmaxnlim/=(double)stand.npatch();
					heightindiv_total/=(double)stand.npatch();

					if (!negligible(standpft_cmass_leaf))
						standpft_vmaxnlim /= standpft_cmass_leaf;

					//Update landcover totals
					landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
					landcover_cmass_leaf[stand.landcover]+=standpft_cmass_leaf*stand.get_landcover_fraction();
					landcover_nmass_leaf[stand.landcover]+=standpft_nmass_leaf*stand.get_landcover_fraction();
					landcover_cmass_veg[stand.landcover]+=standpft_cmass_veg*stand.get_landcover_fraction();
					landcover_nmass_veg[stand.landcover]+=standpft_nmass_veg*stand.get_landcover_fraction();
					landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
					landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
					landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();
					landcover_aiso[stand.landcover]+=standpft_aiso*stand.get_landcover_fraction();
					landcover_amon[stand.landcover]+=standpft_amon*stand.get_landcover_fraction();
					landcover_nuptake[stand.landcover]+=standpft_nuptake*stand.get_landcover_fraction();
					landcover_vmaxnlim[stand.landcover]+=standpft_vmaxnlim*stand.get_landcover_fraction();

					//Update pft totals

					if(pft.landcover!=CROPLAND  || pft.isintercropgrass) {

						stand_mean_cmass+=standpft_cmass*stand.get_landcover_fraction();
						stand_mean_nmass+=standpft_nmass*stand.get_landcover_fraction();
						stand_mean_cmass_leaf+=standpft_cmass_leaf*stand.get_landcover_fraction();
						stand_mean_nmass_leaf+=standpft_nmass_leaf*stand.get_landcover_fraction();
						stand_mean_cmass_veg+=standpft_cmass_veg*stand.get_landcover_fraction();
						stand_mean_nmass_veg+=standpft_nmass_veg*stand.get_landcover_fraction();
						stand_mean_anpp+=standpft_anpp*stand.get_landcover_fraction();
						stand_mean_lai+=standpft_lai*stand.get_landcover_fraction();
						stand_mean_densindiv_total+=standpft_densindiv_total*stand.get_landcover_fraction();
						stand_mean_aiso+=standpft_aiso*stand.get_landcover_fraction();
						stand_mean_amon+=standpft_amon*stand.get_landcover_fraction();
						stand_mean_nuptake+=standpft_nuptake*stand.get_landcover_fraction();
						stand_mean_vmaxnlim+=standpft_vmaxnlim*stand.get_landcover_fraction();

						stand.anpp += standpft_anpp;

						if (vegmode==COHORT || vegmode==INDIVIDUAL)
							for (c=0;c<nclass;c++)
								stand_mean_densindiv_ageclass[c]+=standpft_densindiv_ageclass[c]*stand.get_landcover_fraction();
					}
					else {		//Not mean value for crop stands, but value for unique stand with crop pft
						stand_mean_cmass+=standpft_cmass;		
						stand_mean_nmass+=standpft_nmass;
						stand_mean_cmass_leaf+=standpft_cmass_leaf;
						stand_mean_nmass_leaf+=standpft_nmass_leaf;
						stand_mean_cmass_veg+=standpft_cmass_veg;
						stand_mean_nmass_veg+=standpft_nmass_veg;
						stand_mean_anpp+=standpft_anpp;
						stand_mean_lai+=standpft_lai;
						stand_mean_yield+=standpft_yield;
						stand_mean_yield1+=standpft_yield1;
						stand_mean_yield2+=standpft_yield2;
						stand_mean_densindiv_total+=standpft_densindiv_total;
						stand_mean_aiso+=standpft_aiso;
						stand_mean_amon+=standpft_amon;
						stand_mean_nuptake+=standpft_nuptake;
						stand_mean_vmaxnlim+=standpft_vmaxnlim;

						if (vegmode==COHORT || vegmode==INDIVIDUAL)
							for (c=0;c<nclass;c++)
								stand_mean_densindiv_ageclass[c]+=standpft_densindiv_ageclass[c];
					}


					// Update gridcell totals
					double fraction_of_gridcell = stand.get_gridcell_fraction();

					cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
					nmass_gridcell+=standpft_nmass*fraction_of_gridcell;
					cmass_leaf_gridcell+=standpft_cmass_leaf*fraction_of_gridcell;
					nmass_leaf_gridcell+=standpft_nmass_leaf*fraction_of_gridcell;
					cmass_veg_gridcell+=standpft_cmass_veg*fraction_of_gridcell;
					nmass_veg_gridcell+=standpft_nmass_veg*fraction_of_gridcell;
					anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
					lai_gridcell+=standpft_lai*fraction_of_gridcell;
					dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
					aiso_gridcell+=standpft_aiso*fraction_of_gridcell;
					amon_gridcell+=standpft_amon*fraction_of_gridcell;
					nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
					vmaxnlim_gridcell+=standpft_vmaxnlim*standpft_cmass_leaf*fraction_of_gridcell;
			
					// Graphical output every 10 years
					// (Windows shell only - "plot" statements have no effect otherwise)
					if (!(date.year%10)) {
						plot("C mass [kg C/m2]",pft.name,date.year,stand_mean_cmass);
						plot("NPP [kg C/m2/yr]",pft.name,date.year,stand_mean_anpp);
						plot("LAI [m2/m2]",pft.name,date.year,stand_mean_lai);
						plot("dens [indiv/ha]",pft.name,date.year,stand_mean_densindiv_total*m2toha);
						if (stand_mean_cmass_leaf > 0.0 && ifnlim) {
							plot("vmax nitrogen lim [dimless]",pft.name,date.year,stand_mean_vmaxnlim);
							plot("leaf C:N ratio [kg C/kg N]",pft.name,date.year,stand_mean_cmass_leaf/stand_mean_nmass_leaf);
						}
					}

					if(printseparatestands) {

						int id = stand.id;

						if(stand.landcover == NATURAL) {

							if(!out_anpp_stand_natural[id].invalid())
								out.add_value(out_anpp_stand_natural[id],      standpft_anpp);
						}
						else if(stand.landcover == FOREST) {

							if(!out_anpp_stand_forest[id].invalid())
								out.add_value(out_anpp_stand_forest[id],      standpft_anpp);
						}
					}

				}//if(active)
				gridcell.nextobj();
			}//End of loop through stands

			// Print PFT sums to files

			double stand_mean_cton_leaf = limited_cton(stand_mean_cmass_leaf, stand_mean_nmass_leaf);
			double stand_mean_cton_veg = limited_cton(stand_mean_cmass_veg, stand_mean_nmass_veg);

			out.add_value(out_cmass,     stand_mean_cmass);
			out.add_value(out_anpp,      stand_mean_anpp);
			out.add_value(out_dens,      stand_mean_densindiv_total);
			out.add_value(out_lai,       stand_mean_lai);
			out.add_value(out_aiso,      stand_mean_aiso);
			out.add_value(out_amon,      stand_mean_amon);
			out.add_value(out_cton_leaf, stand_mean_cton_leaf);
			out.add_value(out_cton_veg,  stand_mean_cton_veg);
			out.add_value(out_vmaxnlim,  stand_mean_vmaxnlim);
			out.add_value(out_nuptake,   stand_mean_nuptake * m2toha);

			if (pft.landcover == CROPLAND)
			{
				out.add_value(out_yield,   stand_mean_yield);
				out.add_value(out_yield1,  stand_mean_yield1);
				out.add_value(out_yield2,  stand_mean_yield2);
			}

			// print species heights
			double height = 0.0;
			if (stand_mean_densindiv_total > 0.0)
				height = heightindiv_total / stand_mean_densindiv_total;

			out.add_value(out_speciesheights, height);

			// print crop sowing- and harvest dates
			if (pft.landcover==CROPLAND) {
				int pft_sdate1=-1;
				int pft_sdate2=-1;
				int pft_hdate1=-1;
				int pft_hdate2=-1;
				int pft_lgp=-1;
				double pft_phu=-1;
				double pft_fphu=-1;
				double pft_fhi=-1;

				gridcell.firstobj();
				while (gridcell.isobj) {
					Stand& stand=gridcell.getobj();

					if(stand.cftid==pft.cftid) {
						pft_sdate1=stand[0].pft[pft.id].cropphen->sdate_thisyear[0];
						pft_sdate2=stand[0].pft[pft.id].cropphen->sdate_thisyear[1];
						pft_hdate1=stand[0].pft[pft.id].cropphen->hdate_harvest[0];
						pft_hdate2=stand[0].pft[pft.id].cropphen->hdate_harvest[1];
						pft_lgp=stand[0].pft[pft.id].cropphen->lgp;
						pft_phu=stand[0].pft[pft.id].cropphen->phu;
						pft_fphu=stand[0].pft[pft.id].cropphen->fphu_harv;
						pft_fhi=stand[0].pft[pft.id].cropphen->fhi_harv;
					}
					gridcell.nextobj();
				}
				out.add_value(out_sdate1, pft_sdate1);
				out.add_value(out_sdate2, pft_sdate2);
				out.add_value(out_hdate1, pft_hdate1);
				out.add_value(out_hdate2, pft_hdate2);
				out.add_value(out_lgp,	  pft_lgp);
				out.add_value(out_phu,	  pft_phu);
				out.add_value(out_fphu,	  pft_fphu);
				out.add_value(out_fhi,	  pft_fhi);
			}

			pftlist.nextobj();
		
		} // *** End of PFT loop ***

		flux_veg = flux_soil = flux_fire = flux_est = flux_seed = flux_charvest = 0.0;

		// guess2008 - carbon pools
		c_litter = c_fast = c_slow = c_harv_slow = 0.0;

		surfsoillitterc = surfsoillittern = cwdc = cwdn = centuryc = centuryn = n_litter = n_harv_slow = availn = 0.0;
		andep_gridcell = anfert_gridcell = anmin_gridcell = animm_gridcell = anfix_gridcell = 0.0;
		n_org_leach_gridcell = n_min_leach_gridcell = 0.0;
		flux_nh3 = flux_no = flux_no2 = flux_n2o = flux_nsoil = flux_ntot = flux_nharvest = 0.0;

		double flux_veg_lc[NLANDCOVERTYPES], flux_soil_lc[NLANDCOVERTYPES], flux_fire_lc[NLANDCOVERTYPES], flux_est_lc[NLANDCOVERTYPES], flux_seed_lc[NLANDCOVERTYPES];
		double flux_charvest_lc[NLANDCOVERTYPES];
		double c_litter_lc[NLANDCOVERTYPES], c_fast_lc[NLANDCOVERTYPES], c_slow_lc[NLANDCOVERTYPES], c_harv_slow_lc[NLANDCOVERTYPES];
		double surfsoillitterc_lc[NLANDCOVERTYPES], cwdc_lc[NLANDCOVERTYPES], centuryc_lc[NLANDCOVERTYPES];

		for (int i=0; i<NLANDCOVERTYPES; i++) {
			flux_veg_lc[i]=0.0;
			flux_soil_lc[i]=0.0;
			flux_fire_lc[i]=0.0;
			flux_est_lc[i]=0.0;
			flux_seed_lc[i]=0.0;
			flux_charvest_lc[i]=0.0;

			c_litter_lc[i]=0.0;
			c_fast_lc[i]=0.0;
			c_slow_lc[i]=0.0;
			c_harv_slow_lc[i]=0.0;
			surfsoillitterc_lc[i]=0.0;
			cwdc_lc[i]=0.0;
			centuryc_lc[i]=0.0;
		}

		// Sum C fluxes, dead C pools and runoff across patches

		gridcell.firstobj();

		// Loop through Stands
		while (gridcell.isobj) {
			Stand& stand = gridcell.getobj();
			stand.firstobj();

			//Loop through Patches
			while (stand.isobj) {
				Patch& patch = stand.getobj();

				double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

				flux_veg+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
				flux_soil+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
				flux_fire+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
				flux_est+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
				flux_seed+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
				flux_charvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;
				flux_nharvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTN)*to_gridcell_average;
				flux_nh3+=patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE)*to_gridcell_average;
				flux_no+=patch.fluxes.get_annual_flux(Fluxes::NO_FIRE)*to_gridcell_average;
				flux_no2+=patch.fluxes.get_annual_flux(Fluxes::NO2_FIRE)*to_gridcell_average;
				flux_n2o+=patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE)*to_gridcell_average;
				flux_nsoil+=patch.fluxes.get_annual_flux(Fluxes::N_SOIL)*to_gridcell_average;	
				flux_ntot+=(patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE) + 
				           patch.fluxes.get_annual_flux(Fluxes::NO_FIRE) + 
				           patch.fluxes.get_annual_flux(Fluxes::NO2_FIRE) +
				           patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE) +
				           patch.fluxes.get_annual_flux(Fluxes::N_SOIL)) * to_gridcell_average;

				flux_veg_lc[stand.landcover]+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
				flux_soil_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
				flux_fire_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
				flux_est_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
				flux_seed_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::SEEDC)*to_gridcell_average;
				flux_charvest_lc[stand.landcover]+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;

				c_fast+=patch.soil.cpool_fast*to_gridcell_average;
				c_slow+=patch.soil.cpool_slow*to_gridcell_average;

				c_fast_lc[stand.landcover]+=patch.soil.cpool_fast*to_gridcell_average;
				c_slow_lc[stand.landcover]+=patch.soil.cpool_slow*to_gridcell_average;

				// Sum all litter
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft = patch.pft[q];
					c_litter += (patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr)  * to_gridcell_average;
					c_litter_lc[stand.landcover] += (patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr)  * to_gridcell_average;
					n_litter += (patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap + patchpft.nmass_litter_heart) * to_gridcell_average;
				}

				//Sum slow pools of harvested products
				if(run_landcover && ifslowharvestpool) {
					for (int q=0;q<npft;q++) {
						Patchpft& patchpft=patch.pft[q];
						c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
						c_harv_slow_lc[stand.landcover]+=patchpft.harvested_products_slow*to_gridcell_average;		  //slow pool in receiving landcover (1)
//						c_harv_slow_lc[patchpft.pft.landcover]+=patchpft.harvested_products_slow*to_gridcell_average; //slow pool in donating landcover (2)
						n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
					}
				}

				//Gridcell irrigation
				irrigation_gridcell += patch.irrigation_y*to_gridcell_average;

				surfrunoff_gridcell+=patch.asurfrunoff*to_gridcell_average;
				drainrunoff_gridcell+=patch.adrainrunoff*to_gridcell_average;
				baserunoff_gridcell+=patch.abaserunoff*to_gridcell_average;
				runoff_gridcell+=patch.arunoff*to_gridcell_average;
	
				// Fire return time
				if (!iffire || patch.fireprob < 0.001)
					firert_gridcell+=1000.0 * to_gridcell_average; // Set a limit of 1000 years
				else	
					firert_gridcell+=(1.0/patch.fireprob) * to_gridcell_average;


				andep_gridcell += stand.gridcell.climate.andep * to_gridcell_average;
				anfert_gridcell += stand.gridcell.climate.anfert * to_gridcell_average;
				anmin_gridcell += patch.soil.anmin * to_gridcell_average;
				animm_gridcell += patch.soil.animmob * to_gridcell_average;
				anfix_gridcell += patch.soil.anfix * to_gridcell_average;
				n_min_leach_gridcell += patch.soil.aminleach * to_gridcell_average;
				n_org_leach_gridcell += patch.soil.aorgleach * to_gridcell_average;
				availn += (patch.soil.nmass_avail + patch.soil.snowpack_nmass) * to_gridcell_average;

				for (int r = 0; r < NSOMPOOL-1; r++) {

					if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
						surfsoillitterc += patch.soil.sompool[r].cmass * to_gridcell_average;
						surfsoillitterc_lc[stand.landcover] += patch.soil.sompool[r].cmass * to_gridcell_average;
						surfsoillittern += patch.soil.sompool[r].nmass * to_gridcell_average;
					}
					else if (r == SURFFWD || r == SURFCWD) {
						cwdc += patch.soil.sompool[r].cmass * to_gridcell_average;
						cwdc_lc[stand.landcover] += patch.soil.sompool[r].cmass * to_gridcell_average;
						cwdn += patch.soil.sompool[r].nmass * to_gridcell_average;
					}
					else {	
						centuryc += patch.soil.sompool[r].cmass * to_gridcell_average;
						centuryc_lc[stand.landcover] += patch.soil.sompool[r].cmass * to_gridcell_average;
						centuryn += patch.soil.sompool[r].nmass * to_gridcell_average;
					}
				}

				// Monthly output variables

				for (m=0;m<12;m++) {
					maet[m] += patch.maet[m]*to_gridcell_average;
					mpet[m] += patch.mpet[m]*to_gridcell_average;
					mevap[m] += patch.mevap[m]*to_gridcell_average;
					mintercep[m] += patch.mintercep[m]*to_gridcell_average;
					mrunoff[m] += patch.mrunoff[m]*to_gridcell_average;
					mrh[m] += patch.fluxes.get_monthly_flux(Fluxes::SOILC, m)*to_gridcell_average;
					mwcont_upper[m] += patch.soil.mwcont[m][0]*to_gridcell_average;
					mwcont_lower[m] += patch.soil.mwcont[m][1]*to_gridcell_average;

					mgpp[m] += patch.fluxes.get_monthly_flux(Fluxes::GPP, m)*to_gridcell_average;
					mra[m] += patch.fluxes.get_monthly_flux(Fluxes::RA, m)*to_gridcell_average;

					miso[m]+=patch.fluxes.get_monthly_flux(Fluxes::ISO, m)*to_gridcell_average;
					mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MON, m)*to_gridcell_average;
				}


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
			gridcell.nextobj();
		} // stand loop


		// In contrast to annual NEE, monthly NEE does not include fire
		// or establishment fluxes
		for (m=0;m<12;m++) {
			mnpp[m] = mgpp[m] - mra[m];
			mnee[m] = mnpp[m] - mrh[m];
		}

		// Print gridcell totals to files

		// Determine total leaf C:N ratio
		double cton_leaf_gridcell = limited_cton(cmass_leaf_gridcell, nmass_leaf_gridcell);
		
		// Determine total vegetation C:N ratio
		double cton_veg_gridcell = limited_cton(cmass_veg_gridcell, nmass_veg_gridcell);

		// Determine total vmax nitrogen limitation
		if (cmass_leaf_gridcell > 0.0) {
			vmaxnlim_gridcell /= cmass_leaf_gridcell;
		}

		out.add_value(out_cmass,  cmass_gridcell);
		out.add_value(out_anpp,   anpp_gridcell);
		out.add_value(out_dens,   dens_gridcell);
		out.add_value(out_lai,    lai_gridcell);
		out.add_value(out_firert, firert_gridcell);
		out.add_value(out_runoff, surfrunoff_gridcell);
		out.add_value(out_runoff, drainrunoff_gridcell);
		out.add_value(out_runoff, baserunoff_gridcell);
		out.add_value(out_runoff, runoff_gridcell);
		out.add_value(out_aiso,   aiso_gridcell);
		out.add_value(out_amon,   amon_gridcell);

		if(printseparatestands) {

			gridcell.firstobj();
			while (gridcell.isobj) {

				Stand& stand=gridcell.getobj();
				int id = stand.id;;

				if(stand.landcover==NATURAL) {

					if(!out_anpp_stand_natural[id].invalid())
						out.add_value(out_anpp_stand_natural[id],   stand.anpp);
				}
				else if(stand.landcover == FOREST) {

					if(!out_anpp_stand_forest[id].invalid())
						out.add_value(out_anpp_stand_forest[id],    stand.anpp);
				}

				gridcell.nextobj();
			}
		}

		out.add_value(out_cton_leaf, cton_leaf_gridcell);
		out.add_value(out_vmaxnlim,  vmaxnlim_gridcell);
		out.add_value(out_cton_veg,  cton_veg_gridcell);
		out.add_value(out_nuptake,   nuptake_gridcell * m2toha);

		out.add_value(out_nsources, andep_gridcell * m2toha);
		out.add_value(out_nsources, anfix_gridcell * m2toha);
		out.add_value(out_nsources, anfert_gridcell * m2toha);
		out.add_value(out_nsources, (andep_gridcell + anfix_gridcell + anfert_gridcell) * m2toha);
		out.add_value(out_nsources, anmin_gridcell * m2toha);
		out.add_value(out_nsources, animm_gridcell * m2toha);
		out.add_value(out_nsources, (anmin_gridcell - animm_gridcell) * m2toha);
		out.add_value(out_nsources, (anmin_gridcell - animm_gridcell + andep_gridcell + anfix_gridcell + anfert_gridcell) * m2toha);

		if (run_landcover && run[CROPLAND]) {
			out.add_value(out_irrigation,   irrigation_gridcell);
		}

		// Print landcover totals to files
		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {
					out.add_value(out_cmass, landcover_cmass[i]);
					out.add_value(out_anpp,  landcover_anpp[i]);
					out.add_value(out_dens,  landcover_densindiv_total[i]);
					out.add_value(out_lai,   landcover_lai[i]);
					out.add_value(out_aiso,  landcover_aiso[i]);
					out.add_value(out_amon,  landcover_amon[i]);

					double landcover_cton_leaf = limited_cton(landcover_cmass_leaf[i], landcover_nmass_leaf[i]);
					double landcover_cton_veg = limited_cton(landcover_cmass_veg[i], landcover_nmass_veg[i]);

					if (landcover_cmass_leaf[i] > 0.0) {
						landcover_vmaxnlim[i] /= landcover_cmass_leaf[i];
					}

					out.add_value(out_cton_leaf, landcover_cton_leaf);
					out.add_value(out_cton_veg,  landcover_cton_veg);
					out.add_value(out_vmaxnlim,  landcover_vmaxnlim[i]);
					out.add_value(out_nuptake,   landcover_nuptake[i] * m2toha);
				}
			}
		}

		// Print monthly output variables
		for (m=0;m<12;m++) {
			 out.add_value(out_mnpp,         mnpp[m]);
			 out.add_value(out_mlai,         mlai[m]);
			 out.add_value(out_mgpp,         mgpp[m]);
			 out.add_value(out_mra,          mra[m]);
			 out.add_value(out_maet,         maet[m]);
			 out.add_value(out_mpet,         mpet[m]);
			 out.add_value(out_mevap,        mevap[m]);
			 out.add_value(out_mrunoff,      mrunoff[m]);
			 out.add_value(out_mintercep,    mintercep[m]);
			 out.add_value(out_mrh,          mrh[m]);
			 out.add_value(out_mnee,         mnee[m]);
			 out.add_value(out_mwcont_upper, mwcont_upper[m]);
			 out.add_value(out_mwcont_lower, mwcont_lower[m]);
			 out.add_value(out_miso,         miso[m]);
			 out.add_value(out_mmon,         mmon[m]);
		}


		// Graphical output every 10 years
		// (Windows shell only - no effect otherwise)

		if (!(date.year%10)) {
			gridcell.firstobj();
			if(gridcell.isobj)	//Fixed bug here if no stands were present.
			{
				Stand& stand=gridcell.getobj();
				plot("C flux [kg C/m2/yr]","flux_veg",  date.year, flux_veg);
				plot("C flux [kg C/m2/yr]","flux_soil", date.year, flux_soil);
				plot("C flux [kg C/m2/yr]","flux_fire", date.year, flux_fire);
				plot("C flux [kg C/m2/yr]","flux_est",  date.year, flux_est);
				plot("C flux [kg C/m2/yr]","NEE",       date.year, flux_veg + flux_soil + flux_fire + flux_est);

				if (!ifcentury) {
					plot("Soil C [kg C/m2]","slow", date.year, stand[0].soil.cpool_slow);
					plot("Soil C [kg C/m2]","fast", date.year, stand[0].soil.cpool_fast);
				}
				else {
					plot("N flux (kg N/ha/yr)","Fix",   date.year, -anfix_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Dep",   date.year, -andep_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Fert",  date.year, -anfert_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Leach", date.year, (n_min_leach_gridcell + n_org_leach_gridcell) * m2toha);
					plot("N flux (kg N/ha/yr)","Flux",  date.year, flux_ntot * m2toha);

					plot("N mineralization [kg N/ha/yr]","N", date.year, (anmin_gridcell - animm_gridcell) * m2toha);

					plot("Soil C [kg C/m2]","fine litter",   date.year, surfsoillitterc);
					plot("Soil C [kg C/m2]","coarse litter", date.year, cwdc);
					plot("Soil C [kg C/m2]","soil",          date.year, centuryc); 
					plot("Soil C [kg C/m2]","total",         date.year, surfsoillitterc + cwdc + centuryc); 

					plot("Soil N [kg N/ha]","fine litter",   date.year, surfsoillittern);
					plot("Soil N [kg N/ha]","coarse litter", date.year, cwdn);
					plot("Soil N [kg N/ha]","soil",          date.year, centuryn); 
					plot("Soil N [kg N/ha]","total",         date.year, surfsoillittern + cwdn + centuryn); 
				}
			}
		}

		// Write fluxes to file

		out.add_value(out_cflux, flux_veg);
		out.add_value(out_cflux, flux_soil);
		out.add_value(out_cflux, flux_fire);
		out.add_value(out_cflux, flux_est);
		if (run_landcover) {
			 out.add_value(out_cflux, flux_seed);
			 out.add_value(out_cflux, flux_charvest);
			 out.add_value(out_cflux, gridcell.acflux_landuse_change);
			 out.add_value(out_cflux, gridcell.acflux_harvest_slow);
		}
		out.add_value(out_cflux, flux_veg + flux_soil + flux_fire + flux_est + flux_seed + flux_charvest + gridcell.acflux_landuse_change + gridcell.acflux_harvest_slow);

		// Print C fluxes to per-landcover files
		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {

					GuessOutput::Table* table_p=NULL;

					switch (i)
					{
					case CROPLAND:
						table_p=&out_cflux_cropland;
						break;
					case PASTURE:
						table_p=&out_cflux_pasture;
						break;
					case NATURAL:
						table_p=&out_cflux_natural;
						break;
					case FOREST:
						table_p=&out_cflux_forest;
						break;
					default:
						if(date.year == nyear_spinup)
							dprintf("Modify code to deal with landcover output!\n");
					}

					if(table_p) {
					out.add_value(*table_p, flux_veg_lc[i]);
						out.add_value(*table_p, flux_soil_lc[i]);
						out.add_value(*table_p, flux_fire_lc[i]);
						out.add_value(*table_p, flux_est_lc[i]);
						if (run_landcover) {
							 out.add_value(*table_p, flux_seed_lc[i]);
							 out.add_value(*table_p, flux_charvest_lc[i]);
							 out.add_value(*table_p, gridcell.acflux_landuse_change_lc[i]);
							 out.add_value(*table_p, gridcell.acflux_harvest_slow_lc[i]);
						}
					}

					double cflux_total = flux_veg_lc[i] + flux_soil_lc[i] + flux_fire_lc[i] + flux_est_lc[i];
					if (run_landcover) {
						cflux_total += flux_seed_lc[i];
						cflux_total += flux_charvest_lc[i];
						cflux_total += gridcell.acflux_landuse_change_lc[i];
						cflux_total += gridcell.acflux_harvest_slow_lc[i];
					}
					if(table_p)
						out.add_value(*table_p,  cflux_total);
				}
			}
		}

		out.add_value(out_nflux, -andep_gridcell * m2toha);
		out.add_value(out_nflux, -anfix_gridcell * m2toha);
		out.add_value(out_nflux, -anfert_gridcell * m2toha);
		out.add_value(out_nflux, flux_ntot * m2toha);
		out.add_value(out_nflux, (n_min_leach_gridcell + n_org_leach_gridcell) * m2toha);
		if (run_landcover) {
			 out.add_value(out_nflux, flux_nharvest * m2toha);
		}
		out.add_value(out_nflux, (flux_nharvest + flux_ntot + n_min_leach_gridcell + n_org_leach_gridcell - (andep_gridcell + anfix_gridcell + anfert_gridcell)) * m2toha);

		out.add_value(out_cpool, cmass_gridcell);
		out.add_value(out_cpool, c_litter);
		if (!ifcentury) {
			out.add_value(out_cpool, c_fast);
			out.add_value(out_cpool, c_slow);
		}
		else {
			out.add_value(out_cpool, surfsoillitterc);
			out.add_value(out_cpool, cwdc);
			out.add_value(out_cpool, centuryc);
		}

		if (run_landcover && ifslowharvestpool) {
			out.add_value(out_cpool, c_harv_slow);
		}

		// Calculate total cpool, starting with cmass and litter...
		double cpool_total = cmass_gridcell + c_litter;

		// Add SOM pools
		if (!ifcentury) {
			cpool_total += c_fast + c_slow;
		}
		else {
			cpool_total += centuryc + surfsoillitterc + cwdc;
		}

		// Add slow harvest pool if needed
		if (run_landcover && ifslowharvestpool) {
			cpool_total += c_harv_slow;
		}

		out.add_value(out_cpool, cpool_total);

		// Print C pools to per-landcover files
		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {

					GuessOutput::Table* table_p=NULL;

					switch (i)
					{
					case CROPLAND:
						table_p=&out_cpool_cropland;
						break;
					case PASTURE:
						table_p=&out_cpool_pasture;
						break;
					case NATURAL:
						table_p=&out_cpool_natural;
						break;
					case FOREST:
						table_p=&out_cpool_forest;
						break;
					default:
						if(date.year == nyear_spinup)
							dprintf("Modify code to deal with landcover output!\n");
					}

					if(table_p) {
						out.add_value(*table_p, landcover_cmass[i] * gridcell.landcoverfrac[i]);
						out.add_value(*table_p, c_litter_lc[i]);

						if (!ifcentury) {
							out.add_value(*table_p, c_fast_lc[i]);
							out.add_value(*table_p, c_slow_lc[i]);
						}
						else {
							out.add_value(*table_p, surfsoillitterc_lc[i]);
							out.add_value(*table_p, cwdc_lc[i]);
							out.add_value(*table_p, centuryc_lc[i]);
						}

						if (run_landcover && ifslowharvestpool) {
							out.add_value(*table_p, c_harv_slow_lc[i]);
						}
					}

					// Calculate total cpool, starting with cmass and litter...
					double cpool_total = landcover_cmass[i] * gridcell.landcoverfrac[i] + c_litter_lc[i];

					// Add SOM pools
					if (!ifcentury) {
						cpool_total += c_fast_lc[i] + c_slow_lc[i];
					}
					else {
						cpool_total += centuryc_lc[i] + surfsoillitterc_lc[i] + cwdc_lc[i];
					}

					// Add slow harvest pool if needed
					if (run_landcover && ifslowharvestpool) {
						cpool_total += c_harv_slow_lc[i];
					}
					if(table_p)
						out.add_value(*table_p, cpool_total);

				}
			}
		}



		if (ifcentury) {
			out.add_value(out_npool, nmass_gridcell);
			out.add_value(out_npool, n_litter);
			out.add_value(out_npool, surfsoillittern);
			out.add_value(out_npool, cwdn);
			out.add_value(out_npool, centuryn + availn);

			if(run_landcover && ifslowharvestpool) {
				out.add_value(out_npool, n_harv_slow);
				out.add_value(out_npool, (nmass_gridcell + n_litter + surfsoillittern + cwdn + centuryn + availn + n_harv_slow));
			}
			else {
				out.add_value(out_npool, (nmass_gridcell + n_litter + surfsoillittern + cwdn + centuryn + availn));
			}
		}

		out.add_value(out_ngases, flux_nh3   * m2toha);
		out.add_value(out_ngases, flux_no    * m2toha);
		out.add_value(out_ngases, flux_no2   * m2toha);
		out.add_value(out_ngases, flux_n2o   * m2toha);
		out.add_value(out_ngases, flux_nsoil * m2toha);
		out.add_value(out_ngases, flux_ntot  * m2toha);

		// Output of age structure (Windows shell only - no effect otherwise)

		if (vegmode==COHORT || vegmode==INDIVIDUAL) {

			if (!(date.year%20) && date.year<2000) {
			
				resetwindow("Age structure [yr]");

				pftlist.firstobj();
				while (pftlist.isobj) {
					Pft& pft=pftlist.getobj();

					if (pft.lifeform==TREE) {

						Gridcellpft& gridcellpft=gridcell.pft[pft.id];

						for (c=0;c<nclass;c++)
							plot("Age structure [yr]",pft.name,
								c * estinterval + estinterval / 2,
								stand_mean_densindiv_ageclass[c] / (double)npatch);
					}
					
					pftlist.nextobj();
				}
			}
		}

		//Output of seasonality variables
		if (run_landcover && run[CROPLAND]) {
			out.add_value(out_seasonality,   gridcell.climate.seasonality);
			out.add_value(out_seasonality,   gridcell.climate.var_temp);
			out.add_value(out_seasonality,   gridcell.climate.var_prec);
			out.add_value(out_seasonality,   gridcell.climate.mtemp_min20);
			out.add_value(out_seasonality,   gridcell.climate.atemp_mean);
			out.add_value(out_seasonality,   gridcell.climate.temp_seasonality);
			out.add_value(out_seasonality,   gridcell.climate.mprec_petmin20);
			out.add_value(out_seasonality,   gridcell.climate.aprec);
			out.add_value(out_seasonality,   gridcell.climate.prec_range);
//			out.add_value(out_seasonality,   gridcell.climate.biseasonal);	//Not implemented in this version
		}
	}
}

void CommonOutput::outdaily(Gridcell& gridcell) {
	
}

void CommonOutput::openlocalfiles(Gridcell& gridcell) {

	if(!printseparatestands)
		return;

	if (date.year >= nyear_spinup) {

		int nnaturalstands = 0;
		int nforeststands = 0;
		bool open_natural = false;
		bool open_forest = false;
		double lon = gridcell.get_lon();
		double lat = gridcell.get_lat();


		gridcell.firstobj();
		while (gridcell.isobj) //Loop through stands:
		{
			Stand& stand=gridcell.getobj();
			stand.anpp=0.0;

			if(stand.landcover == NATURAL) {
				nnaturalstands++;
				if(stand.first_year == date.year) {
					open_natural = true;
				}
			}
			else if(stand.landcover == FOREST) {
				nforeststands++;
				if(stand.first_year == date.year) {
					open_forest = true;
				}
			}

			gridcell.nextobj();
		}

#ifdef PRINTFIRSTSTANDFROM1901
		if(date.year == nyear_spinup) {
			open_natural = true;
			open_forest = true;
		}
#endif

		if(open_natural || open_forest) {

			gridcell.firstobj();
			while(gridcell.isobj) {

				Stand& stand=gridcell.getobj();

				int id = stand.id;
				char outfilename[100]={'\0'}, buffer[50]={'\0'};

				sprintf(buffer, "%.1f_%.1f_%d",lon, lat, id);
				strcat(buffer, ".out");

				// create a vector with the pft names
				std::vector<std::string> pfts;

				pftlist.firstobj();
				while (pftlist.isobj) {

					 Pft& pft=pftlist.getobj();	 
					 Standpft& standpft=stand.pft[pft.id];

					 if(standpft.active)
						 pfts.push_back((char*)pft.name);

					 pftlist.nextobj();
				}
				ColumnDescriptors anpp_columns;
				anpp_columns += ColumnDescriptors(pfts,               8, 3);
				anpp_columns += ColumnDescriptor("Total",             8, 3);

				if(open_natural && stand.landcover == NATURAL) {

					strcpy(outfilename, "anpp_natural_");
					strcat(outfilename, buffer);

					if(out_anpp_stand_natural[id].invalid())
						create_output_table(out_anpp_stand_natural[id], outfilename, anpp_columns);
				}
				else if(open_forest && stand.landcover == FOREST) {

					strcpy(outfilename, "anpp_forest_");
					strcat(outfilename, buffer);

					if(out_anpp_stand_forest[id].invalid())
						create_output_table(out_anpp_stand_forest[id], outfilename, anpp_columns);
				}

				gridcell.nextobj();
			}
		}
	}
}

void CommonOutput::closelocalfiles(Gridcell& gridcell) {

	if(!printseparatestands)
		return;

	for(int id=0;id<MAXNUMBER_STANDS;id++) {

		if(!out_anpp_stand_natural[id].invalid())
			close_output_table(out_anpp_stand_natural[id]);
		if(!out_anpp_stand_forest[id].invalid())
			close_output_table(out_anpp_stand_forest[id]);
	}
}

} // namespace
