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

private:

	/// Defines all output tables
	void define_output_tables();

	// Output file names ...
	xtring file_cmass,file_anpp,file_agpp,file_fpc,file_aaet,file_dens,file_lai,file_cflux,file_doc,file_cpool,file_clitter,file_runoff;
	xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
	xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
	xtring file_firert,file_speciesheights;
	xtring file_wetland_water_added;

	// bvoc
	xtring file_aiso, file_miso, file_amon, file_mmon, file_amon_mt1, file_amon_mt2, file_mmon_mt1, file_mmon_mt2;

	// nitrogen
	xtring file_nmass, file_cton_leaf, file_nsources, file_npool, file_nlitter, file_nuptake, file_vmaxnlim, file_nflux, file_ngases;
		
	// Soil temperature at 25cm depth
	xtring file_msoiltempdepth5, file_msoiltempdepth15, file_msoiltempdepth25, file_msoiltempdepth35, file_msoiltempdepth45, file_msoiltempdepth55, file_msoiltempdepth65, file_msoiltempdepth75, file_msoiltempdepth85, file_msoiltempdepth95, file_msoiltempdepth105, file_msoiltempdepth115, file_msoiltempdepth125, file_msoiltempdepth135, file_msoiltempdepth145;
	
	// Methane fluxes
	xtring file_mch4, file_mch4diff, file_mch4plan, file_mch4ebull; 
	
	// Snow, water table depth and active layer depth
	xtring file_msnow, file_mwtp, file_mald;

	// outdaily files
	xtring file_dlai, file_dflux, file_dcmass_grass;

	// dave outputs
	xtring file_fpar, file_dave_lai, file_dave_cmass, file_dave_cmass_leaf, file_dave_sla, file_dave_nmass, file_dave_wscal, file_fpar_layerwise, file_lai_layerwise, file_agd_g, file_rd_g, file_je, file_vm, file_dave_assim, file_resp, file_dnpp, file_subdaily_an, file_subdaily_rd, file_subdaily_anc, file_subdaily_anj, file_subdaily_gsw, file_subdaily_ci, file_subdaily_vcmax, file_subdaily_jmax, file_temps, file_pars, file_vpds;

	// Output tables
	Table out_cmass, out_anpp, out_agpp, out_fpc, out_aaet, out_dens, out_lai, out_cflux, out_doc, out_cpool, out_clitter, out_firert, out_runoff, out_speciesheights;
	Table out_wetland_water_added;

	Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
	Table out_mrh, out_mnee, out_mwcont_upper, out_mwcont_lower;
	
	// bvoc
	Table out_aiso, out_miso, out_amon, out_mmon, out_amon_mt1, out_amon_mt2, out_mmon_mt1, out_mmon_mt2;

	// nitrogen
	Table out_nmass, out_cton_leaf, out_nsources, out_npool, out_nlitter, out_nuptake, out_vmaxnlim, out_nflux, out_ngases;

	// Methane, snow, water table and active layer depth
	Table out_mch4, out_mch4diff, out_mch4plan, out_mch4ebull, out_msnow, out_mwtp, out_mald;

	// Soil temperatures
	Table out_msoiltempdepth5, out_msoiltempdepth15, out_msoiltempdepth25, out_msoiltempdepth35, out_msoiltempdepth45, out_msoiltempdepth55, out_msoiltempdepth65, out_msoiltempdepth75, out_msoiltempdepth85, out_msoiltempdepth95,
    out_msoiltempdepth105, out_msoiltempdepth115, out_msoiltempdepth125, out_msoiltempdepth135, out_msoiltempdepth145;

	// Outdaily
	Table out_dlai, out_dflux, out_dcmass_grass;

	// Dave
	Table out_fpar, out_dave_lai, out_dave_cmass, out_dave_cmass_leaf, out_dave_sla, out_dave_nmass, out_dave_wscal, out_fpar_layerwise, out_lai_layerwise, out_agd_g, out_rd_g, out_je, out_vm, out_dave_assim, out_resp, out_dnpp, out_subdaily_an, out_subdaily_rd, out_subdaily_anc, out_subdaily_anj, out_subdaily_gsw, out_subdaily_ci, out_subdaily_vcmax, out_subdaily_jmax, out_temps, out_pars, out_vpds;
};

}

#endif // LPJ_GUESS_COMMON_OUTPUT_H
