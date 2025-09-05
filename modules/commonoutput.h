///////////////////////////////////////////////////////////////////////////////////////
/// \file commonoutput.h
/// \brief Output module for the most commonly needed output files
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_COMMON_OUTPUT_H
#define LPJ_GUESS_COMMON_OUTPUT_H

#include "outputmodule.h"
#include "outputchannel.h"
#include "gutil.h"

namespace GuessOutput {

/// Output module for the most commonly needed output files
class CommonOutput : public OutputModule {
public:

	CommonOutput();

	~CommonOutput();

	// implemented functions inherited from OutputModule
	// (see documentation in OutputModule)

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);

	void openlocalfiles(Gridcell& gridcell) {};

	void closelocalfiles(Gridcell& gridcell) {};

private:

	/// Defines all output tables
	void define_output_tables();

	// Output file names ...
	xtring file_cmass,file_cmass_agb,file_anpp,file_agpp,file_agpp_no,file_agpp_ns,file_fpc,file_aaet,file_dens,file_lai,file_cflux,file_doc,file_cpool,file_clitter;
	xtring file_csom, file_ctonsom, file_soildepth, file_csoil, file_csoil_depth, file_nsoil_depth, file_csoil_cm, file_csoilp_cm, file_rh_cm, file_flitter_cm, file_clitter_cm, file_csoil_cm_peat, file_csoilp_cm_peat, file_rh_cm_peat, file_netnminsoil_depth, file_rh, file_rh_depth, file_tau_soil, file_tau_soil_depth;
	xtring file_temp_mod, file_moist_mod, file_tot_mod, file_temp_mod_depth, file_moist_mod_depth, file_tot_mod_depth, file_volume_org_depth, file_volume_min_depth, file_volume_root_depth, file_runoff;
	xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
	xtring file_mnee,file_mnetmin,file_mcton_leaf,file_mvmaxnlim,file_mlambda,file_mwmass,file_mimass,file_mwcont_upper,file_mwcont_lower;
	xtring file_firert,file_speciesheights;
	xtring file_wetland_water_added;
	xtring file_soildeclim, file_soilgrolim, file_soiltotlim, file_soilcuelim, file_clitterinput, file_clitterinput_depth, file_csomlitterinput, file_soilcperma;

	// bvoc
	xtring file_aiso, file_miso, file_amon, file_mmon, file_amon_mt1, file_amon_mt2, file_mmon_mt1, file_mmon_mt2;

	// nitrogen
	xtring file_nmass, file_cton_leaf, file_nsources, file_npool, file_nlitter, file_nuptake, file_vmaxnlim, file_lambda, file_nflux, file_ngases;
	xtring file_soil_npool, file_soil_nflux;
	xtring file_csoillyr;

	// BLAZE & SIMFIRE
	xtring file_aburned_area_out, file_mburned_area_out;
	xtring file_simfireanalysis_out;
		
	// Soil temperature and water
	xtring file_mtemp, file_mairtemp, file_mwmasslyr, file_mimasslyr, file_mwimasslyr, file_wsatslyr, file_mdrainlyr;
	xtring file_mtemp_depth, file_mwmass_depth, file_mimass_depth, file_mwimass_depth, file_mdrain_depth;

	// Methane fluxes
	xtring file_mch4, file_mch4diff, file_mch4plan, file_mch4ebull; 
	
	// Snow, water table depth and active layer depth
	xtring file_msnow, file_msnowcover, file_snowdays, file_mwmass_ssw, file_mwtp, file_mald, file_mald_peat, file_msoiltempdepth25;
	xtring file_wflux, file_wpool, file_wsats, file_wsats_depth;
	xtring file_peat_frac, file_peat_depth, file_csoil_data_top, file_csoil_data_bot;

	// Output tables
	Table out_cmass, out_cmass_agb, out_anpp, out_agpp, out_agpp_no, out_agpp_ns, out_fpc, out_aaet, out_dens, out_lai, out_cflux, out_doc, out_cpool, out_clitter; 
	Table out_csom, out_ctonsom, out_soildepth, out_csoil, out_csoil_depth, out_nsoil_depth, out_csoil_cm, out_csoilp_cm, out_rh_cm, out_flitter_cm, out_clitter_cm, out_csoil_cm_peat, out_csoilp_cm_peat, out_rh_cm_peat, out_netnminsoil_depth, out_rh, out_rh_depth, out_tau_soil, out_tau_soil_depth;
	Table out_temp_mod, out_moist_mod, out_tot_mod, out_temp_mod_depth, out_moist_mod_depth, out_tot_mod_depth, out_volume_org_depth, out_volume_min_depth, out_volume_root_depth, out_firert, out_runoff, out_speciesheights;
	Table out_wetland_water_added;
	Table out_soildeclim, out_soilgrolim, out_soiltotlim, out_soilcuelim, out_clitterinput, out_clitterinput_depth, out_csomlitterinput, out_soilcperma;

	Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
	Table out_mrh, out_mnee, out_mnetmin, out_mcton_leaf, out_mvmaxnlim, out_mlambda, out_mwmass, out_mimass, out_mwcont_upper, out_mwcont_lower;
	
	// bvoc
	Table out_aiso, out_miso, out_amon, out_mmon, out_amon_mt1, out_amon_mt2, out_mmon_mt1, out_mmon_mt2;
	
	Table out_nmass, out_cton_leaf, out_nsources, out_npool, out_nlitter, out_nuptake, out_vmaxnlim, out_lambda, out_nflux, out_ngases;
	Table out_soil_npool, out_soil_nflux;
	Table out_csoillyr;

	// BLAZE && SIMFIRE
	Table out_aburned_area, out_mburned_area;
	Table out_simfireanalysis;

	// Methane, snow, water table and active layer depth
	Table out_mch4, out_mch4diff, out_mch4plan, out_mch4ebull, out_msnow, out_msnowcover, out_snowdays, out_mwmass_ssw, out_mwtp, out_mald, out_mald_peat, out_msoiltempdepth25;
	Table out_wflux, out_wpool, out_wsats, out_wsats_depth;
	Table out_peat_frac, out_peat_depth, out_csoil_data_top, out_csoil_data_bot;

	// Soil temperatures and water
	Table out_mtemp, out_mairtemp, out_mwmasslyr, out_mimasslyr, out_mwimasslyr, out_wsatslyr, out_mdrainlyr;
	Table out_mtemp_depth, out_mwmass_depth, out_mimass_depth, out_mwimass_depth, out_mdrain_depth;
};

}

#endif // LPJ_GUESS_COMMON_OUTPUT_H
