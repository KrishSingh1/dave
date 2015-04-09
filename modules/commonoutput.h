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

// Definitions for separate output files per NATURAL and FOREST stand (when instruction file parameter printseparatestands == true,
// printseparatestands is set to false in LandcoverInputModule::init() when input land cover fraction data file has data for > 50 gridcells)
#define PRINTFIRSTSTANDFROM1901		// Printout of first stand from first historic year
#define MAXNUMBER_STANDS 100		// Upper limit for files in multiple stand printout

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

	void openlocalfiles(Gridcell& gridcell);

	void closelocalfiles(Gridcell& gridcell);

private:

	/// Defines all output tables
	void define_output_tables();

	// Output file names ...
	xtring file_cmass,file_anpp,file_agpp,file_fpc,file_aaet,file_dens,file_lai,file_cflux,file_cpool,file_clitter,file_runoff;
	xtring file_yield, file_yield1, file_yield2, file_sdate1, file_sdate2, file_hdate1, file_hdate2, file_lgp, file_phu, file_fphu, file_fhi, file_irrigation, file_seasonality;
	xtring file_cflux_cropland, file_cflux_pasture, file_cflux_natural, file_cflux_forest, file_cpool_cropland, file_cpool_pasture, file_cpool_natural, file_cpool_forest;
	xtring file_nflux_cropland, file_nflux_pasture, file_nflux_natural, file_nflux_forest, file_npool_cropland, file_npool_pasture, file_npool_natural, file_npool_forest;
	xtring file_anpp_cropland, file_anpp_pasture, file_anpp_natural, file_anpp_forest, file_cmass_cropland, file_cmass_pasture, file_cmass_natural, file_cmass_forest, file_dens_natural, file_dens_forest;
	xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
	xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
	xtring file_firert,file_speciesheights;

	// bvoc
	xtring file_aiso,file_miso,file_amon,file_mmon;

	// nitrogen
	xtring file_nmass, file_cton_leaf, file_nsources, file_npool, file_nlitter, file_nuptake, file_vmaxnlim, file_nflux, file_ngases;
	
	// daily
	xtring file_daily_lai,file_daily_npp,file_daily_nmass,file_daily_cmass,file_daily_cton,file_daily_ndemand;
	xtring file_daily_cmass_leaf,file_daily_nmass_leaf,file_daily_cmass_root,file_daily_nmass_root,file_daily_cmass_stem,file_daily_nmass_stem,file_daily_cmass_storage,file_daily_nmass_storage,file_daily_n_input_soil;
	xtring file_daily_avail_nmass_soil,file_daily_upper_wcont,file_daily_lower_wcont,file_daily_irrigation;
	xtring file_daily_temp,file_daily_prec,file_daily_rad;
	xtring file_daily_cmass_dead_leaf,file_daily_nmass_dead_leaf,file_daily_fphu;
	xtring file_daily_nminleach, file_daily_norgleach, file_daily_nuptake;
	xtring file_daily_ds,file_daily_stem,file_daily_leaf,file_daily_root,file_daily_storage;


	// Output tables
	Table out_cmass, out_anpp, out_agpp, out_fpc, out_aaet, out_dens, out_lai, out_cflux, out_cpool, out_clitter, out_firert, out_runoff, out_speciesheights;
	Table out_yield, out_yield1, out_yield2, out_sdate1, out_sdate2, out_hdate1, out_hdate2, out_lgp, out_phu, out_fhi, out_fphu, out_irrigation, out_seasonality;	
	Table out_cflux_cropland, out_cflux_pasture, out_cflux_natural, out_cflux_forest, out_cpool_cropland, out_cpool_pasture, out_cpool_natural, out_cpool_forest;
	Table out_nflux_cropland, out_nflux_pasture, out_nflux_natural, out_nflux_forest, out_npool_cropland, out_npool_pasture, out_npool_natural, out_npool_forest;
	Table out_anpp_cropland, out_anpp_pasture, out_anpp_natural, out_anpp_forest, out_cmass_cropland, out_cmass_pasture, out_cmass_natural, out_cmass_forest, out_dens_natural, out_dens_forest;
	Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
	Table out_mrh, out_mnee, out_mwcont_upper, out_mwcont_lower;
	Table out_anpp_stand_natural[MAXNUMBER_STANDS];
	Table out_cmass_stand_natural[MAXNUMBER_STANDS];
	Table out_anpp_stand_forest[MAXNUMBER_STANDS];
	Table out_cmass_stand_forest[MAXNUMBER_STANDS];

	// bvoc
	Table out_aiso, out_miso, out_amon, out_mmon;
	
	Table out_nmass, out_cton_leaf, out_nsources, out_npool, out_nlitter, out_nuptake, out_vmaxnlim, out_nflux, out_ngases;

	//daily
	Table out_daily_lai,out_daily_npp,out_daily_cton,out_daily_nmass,out_daily_cmass,out_daily_ndemand;
	Table out_daily_cmass_leaf,out_daily_nmass_leaf,out_daily_cmass_root,out_daily_nmass_root,out_daily_cmass_stem,out_daily_nmass_stem,out_daily_cmass_storage,out_daily_nmass_storage,out_daily_n_input_soil;
	Table out_daily_cmass_dead_leaf,out_daily_nmass_dead_leaf,out_daily_fphu;
	Table out_daily_avail_nmass_soil,out_daily_upper_wcont,out_daily_lower_wcont,out_daily_irrigation;
	Table out_daily_temp,out_daily_prec,out_daily_rad;
	Table out_daily_nminleach, out_daily_norgleach, out_daily_nuptake;
	Table out_daily_ds,out_daily_stem,out_daily_leaf,out_daily_root,out_daily_storage;
};

}

#endif // LPJ_GUESS_COMMON_OUTPUT_H
