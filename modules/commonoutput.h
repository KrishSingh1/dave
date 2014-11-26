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

	void openlocalfiles(Gridcell& gridcell);

	void closelocalfiles(Gridcell& gridcell);

private:

	/// Defines all output tables
	void define_output_tables();

	// Output file names ...
	xtring file_cmass,file_anpp,file_agpp,file_fpc,file_aaet,file_dens,file_lai,file_cflux,file_cpool,file_clitter,file_runoff;
	xtring file_yield, file_yield1, file_yield2, file_sdate1, file_sdate2, file_hdate1, file_hdate2, file_lgp, file_phu, file_fphu, file_fhi, file_irrigation, file_seasonality;
	xtring file_cflux_cropland, file_cflux_pasture, file_cflux_natural, file_cflux_forest, file_cpool_cropland, file_cpool_pasture, file_cpool_natural, file_cpool_forest;
	xtring file_anpp_cropland, file_anpp_pasture, file_anpp_natural, file_anpp_forest, file_cmass_cropland, file_cmass_pasture, file_cmass_natural, file_cmass_forest, file_dens_natural, file_dens_forest;
	xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
	xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
	xtring file_firert,file_speciesheights;

	// bvoc
	xtring file_aiso,file_miso,file_amon,file_mmon;

	// nitrogen
	xtring file_nmass, file_cton_leaf, file_nsources, file_npool, file_nlitter, file_nuptake, file_vmaxnlim, file_nflux, file_ngases;
	
	// Output tables
	Table out_cmass, out_anpp, out_agpp, out_fpc, out_aaet, out_dens, out_lai, out_cflux, out_cpool, out_clitter, out_firert, out_runoff, out_speciesheights;
	Table out_yield, out_yield1, out_yield2, out_sdate1, out_sdate2, out_hdate1, out_hdate2, out_lgp, out_phu, out_fhi, out_fphu, out_irrigation, out_seasonality;	
	Table out_cflux_cropland, out_cflux_pasture, out_cflux_natural, out_cflux_forest, out_cpool_cropland, out_cpool_pasture, out_cpool_natural, out_cpool_forest;
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
};

}

#endif // LPJ_GUESS_COMMON_OUTPUT_H
