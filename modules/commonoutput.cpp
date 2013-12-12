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
	declare_parameter("file_cflux", &file_cflux, 300, "C fluxes output file");
	declare_parameter("file_dens", &file_dens, 300, "Tree density output file");
	declare_parameter("file_cpool", &file_cpool, 300, "Soil C output file");
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

	pftlist.firstobj();
	while (pftlist.isobj) {
		 Pft& pft=pftlist.getobj();

		 pfts.push_back((char*)pft.name);

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
		 cflux_columns += ColumnDescriptor("Harvest",      9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",              10, 5);

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
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
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
	double flux_veg, flux_soil, flux_fire, flux_est, flux_charvest;
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

		double gcpft_cmass=0.0;
		double gcpft_nmass=0.0;
		double gcpft_cmass_leaf=0.0;
		double gcpft_nmass_leaf=0.0;
		double gcpft_cmass_veg=0.0;
		double gcpft_nmass_veg=0.0;
		double gcpft_anpp=0.0;
		double gcpft_lai=0.0;
		double gcpft_densindiv_total=0.0;
		double gcpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double gcpft_aiso=0.0;
		double gcpft_amon=0.0;
		double gcpft_nuptake=0.0;
		double gcpft_vmaxnlim=0.0;

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

		double standpft_cmass=0.0;
		double standpft_nmass=0.0;
		double standpft_cmass_leaf=0.0;
		double standpft_nmass_leaf=0.0;
		double standpft_cmass_veg=0.0;
		double standpft_nmass_veg=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double standpft_aiso=0.0;
		double standpft_amon=0.0;
		double standpft_nuptake=0.0;
		double standpft_vmaxnlim=0.0;


		// *** Loop through PFTs ***

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Gridcellpft& gridcellpft=gridcell.pft[pft.id];

			// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs		
			gcpft_cmass=0.0;
			gcpft_nmass=0.0;
			gcpft_cmass_leaf=0.0;
			gcpft_nmass_leaf=0.0;
			gcpft_cmass_veg=0.0;
			gcpft_nmass_veg=0.0;
			gcpft_anpp=0.0;
			gcpft_lai=0.0;
			gcpft_densindiv_total=0.0;		
			gcpft_aiso=0.0;
			gcpft_amon=0.0;
			gcpft_nuptake=0.0;
			gcpft_vmaxnlim=0.0;

			double heightindiv_total = 0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
				standpft_cmass=0.0;
				standpft_nmass=0.0;
				standpft_cmass_leaf=0.0;
				standpft_nmass_leaf=0.0;
				standpft_cmass_veg=0.0;
				standpft_nmass_veg=0.0;
				standpft_anpp=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total = 0.0;
				standpft_aiso=0.0;
				standpft_amon=0.0;
				standpft_nuptake=0.0;
				standpft_vmaxnlim=0.0;

				// Initialise age structure array

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++){
						standpft_densindiv_ageclass[c] = 0.0;
					}
				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch = stand.getobj();

					standpft_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);
					standpft_aiso += patch.fluxes.get_annual_flux(Fluxes::ISO, pft.id);
					standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MON, pft.id);

					Vegetation& vegetation = patch.vegetation;

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv=vegetation.getobj();
							
						// guess2008 - alive check added
						if (indiv.id!=-1 && indiv.alive) { 
							
							if (indiv.pft.id==pft.id) {
								standpft_cmass += indiv.cmass_leaf + indiv.cmass_root + 
								                  indiv.cmass_wood();
								standpft_nmass += indiv.nmass_leaf + indiv.nmass_root + 
								                  indiv.nmass_wood() + indiv.nstore();
								standpft_cmass_leaf += indiv.cmass_leaf;
								standpft_nmass_leaf += indiv.cmass_leaf / indiv.cton_leaf_aavr;
								standpft_cmass_veg += indiv.cmass_veg;
								standpft_nmass_veg += indiv.nmass_veg;
								standpft_lai += indiv.lai;
								standpft_vmaxnlim += indiv.avmaxnlim * indiv.cmass_leaf;
								standpft_nuptake += indiv.anuptake;

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
				gcpft_cmass+=standpft_cmass;
				gcpft_nmass+=standpft_nmass;
				gcpft_cmass_leaf+=standpft_cmass_leaf;
				gcpft_nmass_leaf+=standpft_nmass_leaf;
				gcpft_cmass_veg+=standpft_cmass_veg;
				gcpft_nmass_veg+=standpft_nmass_veg;
				gcpft_anpp+=standpft_anpp;
				gcpft_lai+=standpft_lai;
				gcpft_densindiv_total+=standpft_densindiv_total;
				gcpft_aiso+=standpft_aiso;
				gcpft_amon+=standpft_amon;
				gcpft_nuptake+=standpft_nuptake;
				gcpft_vmaxnlim+=standpft_vmaxnlim;

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++)
						gcpft_densindiv_ageclass[c] += standpft_densindiv_ageclass[c];

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
					plot("C mass [kg C/m2]",pft.name,date.year,gcpft_cmass);
					plot("NPP [kg C/m2/yr]",pft.name,date.year,gcpft_anpp);
					plot("LAI [m2/m2]",pft.name,date.year,gcpft_lai);
					plot("dens [indiv/ha]",pft.name,date.year,gcpft_densindiv_total*m2toha);
					if (gcpft_cmass_leaf > 0.0 && ifnlim) {
						plot("vmax nitrogen lim [dimless]",pft.name,date.year,gcpft_vmaxnlim);
						plot("leaf C:N ratio [kg C/kg N]",pft.name,date.year,gcpft_cmass_leaf/gcpft_nmass_leaf);
					}
				}
				gridcell.nextobj();
			}//End of loop through stands

			// Print PFT sums to files

			double gcpft_cton_leaf = limited_cton(gcpft_cmass_leaf, gcpft_nmass_leaf);
			double gcpft_cton_veg = limited_cton(gcpft_cmass_veg, gcpft_nmass_veg);
			
			out.add_value(out_cmass,     gcpft_cmass);
			out.add_value(out_anpp,      gcpft_anpp);
			out.add_value(out_dens,	     gcpft_densindiv_total);
			out.add_value(out_lai,       gcpft_lai);
			out.add_value(out_aiso,      gcpft_aiso);
			out.add_value(out_amon,      gcpft_amon);
			out.add_value(out_cton_leaf, gcpft_cton_leaf);
			out.add_value(out_cton_veg,  gcpft_cton_veg);
			out.add_value(out_vmaxnlim,  gcpft_vmaxnlim);
			out.add_value(out_nuptake,   gcpft_nuptake * m2toha);	

			// print species heights
			double height = 0.0;
			if (gcpft_densindiv_total > 0.0)
				height = heightindiv_total / gcpft_densindiv_total;
			
			out.add_value(out_speciesheights, height);

			pftlist.nextobj();

		} // *** End of PFT loop ***

		flux_veg = flux_soil = flux_fire = flux_est = flux_charvest = 0.0;

		// guess2008 - carbon pools
		c_litter = c_fast = c_slow = c_harv_slow = 0.0;

		surfsoillitterc = surfsoillittern = cwdc = cwdn = centuryc = centuryn = n_litter = n_harv_slow = availn = 0.0;
		andep_gridcell = anfert_gridcell = anmin_gridcell = animm_gridcell = anfix_gridcell = 0.0;
		n_org_leach_gridcell = n_min_leach_gridcell = 0.0;
		flux_nh3 = flux_no = flux_no2 = flux_n2o = flux_nsoil = flux_ntot = flux_nharvest = 0.0;

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
				
				c_fast+=patch.soil.cpool_fast*to_gridcell_average;
				c_slow+=patch.soil.cpool_slow*to_gridcell_average;

				// Sum all litter
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft = patch.pft[q];
					c_litter += (patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr)  * to_gridcell_average;
					n_litter += (patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap + patchpft.nmass_litter_heart) * to_gridcell_average;
				}

				//Sum slow pools of harvested products
				if(run_landcover && ifslowharvestpool)
				{
					for (int q=0;q<npft;q++) 
					{
						Patchpft& patchpft=patch.pft[q];
						c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
						n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
					}
				}

				surfrunoff_gridcell+=patch.asurfrunoff*to_gridcell_average;
				drainrunoff_gridcell+=patch.adrainrunoff*to_gridcell_average;
				baserunoff_gridcell+=patch.abaserunoff*to_gridcell_average;
				runoff_gridcell+=patch.arunoff*to_gridcell_average;
	
				// Fire return time
				if (!iffire || patch.fireprob < 0.001)
					firert_gridcell+=1000.0/(double)stand.npatch(); // Set a limit of 1000 years
				else	
					firert_gridcell+=(1.0/patch.fireprob)/(double)stand.npatch();


				andep_gridcell += stand.gridcell.climate.andep / (double)stand.npatch();
				anfert_gridcell += stand.gridcell.climate.anfert / (double)stand.npatch();
				anmin_gridcell += patch.soil.anmin / (double)stand.npatch();
				animm_gridcell += patch.soil.animmob / (double)stand.npatch();
				anfix_gridcell += patch.soil.anfix / (double)stand.npatch();
				n_min_leach_gridcell += patch.soil.aminleach / (double)stand.npatch();
				n_org_leach_gridcell += patch.soil.aorgleach / (double)stand.npatch();
				availn += (patch.soil.nmass_avail + patch.soil.snowpack_nmass)
				          / (double)stand.npatch();

				for (int r = 0; r < NSOMPOOL-1; r++) {
					if (patch.soil.sompool[r].nmass > 0.0) {
						if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
							surfsoillitterc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							surfsoillittern += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
						else if (r == SURFFWD || r == SURFCWD) {
							cwdc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							cwdn += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
						else {	
							centuryc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							centuryn += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
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
			 out.add_value(out_cflux, flux_charvest);
		}
		out.add_value(out_cflux, flux_veg + flux_soil + flux_fire + flux_est + flux_charvest);

		out.add_value(out_nflux, -andep_gridcell * m2toha);
		out.add_value(out_nflux, -anfix_gridcell * m2toha);
		out.add_value(out_nflux, -anfert_gridcell * m2toha);
		out.add_value(out_nflux, flux_ntot * m2toha);
		out.add_value(out_nflux, (n_min_leach_gridcell + n_org_leach_gridcell) * m2toha);
		if (run_landcover) {
			 out.add_value(out_nflux, flux_nharvest);
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
								gcpft_densindiv_ageclass[c] / (double)npatch);
					}
					
					pftlist.nextobj();
				}
			}
		}
	}
}

void CommonOutput::outdaily(Gridcell& gridcell) {
	
}

} // namespace
