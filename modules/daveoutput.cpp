///////////////////////////////////////////////////////////////////////////////////////
/// \file daveoutput.cpp
/// \brief Implementation of the DAVE output module
///
/// \author Drew Holzworth
/// $Date: 2023-02-17 10:58:54 +1100 (Fri, 17 Feb 2023) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <limits.h>
#include <sstream>

#include "config.h"
#include "commonoutput.h"
#include "daveoutput.h"
#include "daveoutputhelpers.h"
#include "gutil.h"
#include "outputchannel.h"
#include "outputmodule.h"
#include "parameters.h"
#include "util.h"

namespace GuessOutput {

// Register the DAVE output module.
REGISTER_OUTPUT_MODULE("dave", DaveOutput)

// #############################################################################
// DaveOutput implementation.
// #############################################################################

DaveOutput::DaveOutput() {
	double t_start = get_wall_time();

	// Annual patch-level outputs.
	add_patch_output("file_dave_patch_age"           , ANNUAL, "age", [](const Patch &p) { return 1.0 * p.age; }, "Annual patch-level age since last disturbance (years)");
	add_patch_output("file_dave_arunoff"             , ANNUAL, runoff_columns, "Annual runoff (mm)");
	add_patch_output("file_dave_acpool"              , ANNUAL, cpool_columns, "Annual C pools (kgC/m2)");
	add_patch_output("file_dave_anpool"              , ANNUAL, npool_columns, "Annual N pools (kgN/m2)");
	add_patch_output("file_dave_acflux"              , ANNUAL, cflux_columns, "Annual C flux (gC/m2)");
	add_patch_output("file_dave_globfirm"            , ANNUAL, globfirm_columns, "Annual GLOBFIRM output file");
	add_patch_output("file_dave_mwcont_upper"        , ANNUAL, mwcont_upper_columns, "Monthly wcont_upper output file (0-1)");
	add_patch_output("file_dave_mwcont_lower"        , ANNUAL, mwcont_lower_columns, "Monthly wcont_lower output file (0-1)");
	add_patch_output("file_dave_apet"                , ANNUAL, pet_columns, "Annual potential patch-level evapotranspiration (mm)");
	add_patch_output("file_dave_asimfire"            , ANNUAL, simfire_columns, "Annual simfire analysis");
	add_patch_output("file_dave_afuel"               , ANNUAL, fuel_columns, "Annual blaze fuel availability (gC/m2)");
	add_patch_output("file_dave_acoarse_woody_debris", ANNUAL, coarse_woody_debris_columns, "Annual coarse woody debris (gC/m2)");
	add_patch_output("file_dave_amet_year"           , ANNUAL, met_year_columns, "Current year of met data being used");
	add_patch_output("file_dave_aco2"                , ANNUAL, "co2", [](const Patch &p) { return p.get_climate().co2; }, "Annual atmospheric co2 concentration (ppm)");
	add_patch_output("file_dave_aminleach"           , ANNUAL, "aminleach", [](const Patch &patch) { return patch.soil.aminleach; }, "Leaching of soil mineral N (kgN/m2/yr)");
	add_patch_output("file_dave_sompool_acmass"      , ANNUAL, sompool_cmass_columns, "Daily SOM pool C mass (kgC/m2)");
	add_patch_output("file_dave_sompool_anmass"      , ANNUAL, sompool_nmass_columns, "Daily SOM pool N mass (kgN/m2)");
	add_patch_output("file_dave_andep"               , ANNUAL, aninput_columns, "Annual nitrogen deposition (kgN/m2)");

	// Annual per-patch individual-level outputs.
	add_aggregated_indiv_output("file_dave_alai"     , ANNUAL, get_lai, "Annual LAI(m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_afpc"     , ANNUAL, [](const Individual &i) { return i.fpc; }, "Annual FPC ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_acmass"   , ANNUAL, [](const Individual &i) { return i.ccont(); }, "Annual cmass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_anmass"   , ANNUAL, [](const Individual &i) { return i.ncont(); }, "Annual nmass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aheight"  , ANNUAL, [](const Individual &i) { return i.height; }, "Annual plant height (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_aaet"     , ANNUAL, [](const Individual &i) { return i.aaet; }, "Annual actual evapotranspiration (mm)", sum, true);
	add_aggregated_indiv_output("file_dave_adensity" , ANNUAL, [](const Individual &i) { return i.densindiv; }, "Annual density of individuals over patch (/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_altor"    , ANNUAL, [](const Individual &i) { return i.ltor; }, "Annual leaf to root ratio (unitless)", mean, true);

	// Annual per-patch pft-level outputs.
	add_pft_output("file_dave_anpp"                  , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id); }, "Total individual-level Annual NPP (kgC/m2/year)", true);
	add_pft_output("file_dave_agpp"                  , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id); }, "Total individual-level Annual GPP (kgC/m2/year)", true);
	add_pft_output("file_dave_aresp"                 , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::RA, pft.id); }, "Total individual-level Annual Leaf Respiration (kgC/m2/year)", true);
	add_pft_output("file_dave_acmass_mort"           , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].cmass_mort; }, "Mass of annual killed vegetation (kgC/m2/year)", true);
	add_pft_output("file_dave_aclitter"              , ANNUAL, get_pft_cmass_litter, "Total individual-level Annual GPP (kgC/m2/year)");
	add_pft_output("file_dave_ancohort"              , ANNUAL, get_nindiv, "Get the number of cohorts of this PFT currently established in this patch");

	// Annual stand-level outputs.
	add_stand_output("file_dave_stand_frac"          , ANNUAL, "fraction", [](const Stand &s) { return s.get_gridcell_fraction(); }, "Fraction of the gridcell occupied by each stand");
	add_stand_output("file_dave_stand_type"          , ANNUAL, "type", [](const Stand &s) { return (double)s.landcover; }, "Stand landcover types");

	// Daily patch-level outputs.
	add_patch_output("file_dave_daylength"           , DAILY, "daylength", [](const Patch &p) { return p.get_climate().daylength; }, "Daily patch-level day-length (h)");
	add_patch_output("file_dave_soil_nmass_avail"    , DAILY, "nmass", [](const Patch &p) { return p.soil.get_nmass_avail(); }, "Available Soil N for plant uptake");
	add_patch_output("file_dave_dsimfire"            , DAILY, simfire_columns, "Daily simfire analysis");
	add_patch_output("file_dave_dfuel"               , DAILY, fuel_columns, "Daily blaze fuel availability (kgC/m2)");
	add_patch_output("file_dave_dcoarse_woody_debris", DAILY, coarse_woody_debris_columns, "Daily coarse woody debris (gC/m2)");
	add_patch_output("file_dave_met_temp"            , DAILY, "temp", [](const Patch &patch) { return patch.stand.get_gridcell().climate.temp;    }, "Daily air temperature (°C)");
	add_patch_output("file_dave_met_par"             , DAILY, "par", [](const Patch &patch) { return patch.stand.get_gridcell().climate.par;     }, "Daily PAR (J/m2/timestep)");
	add_patch_output("file_dave_met_vpd"             , DAILY, "vpd", [](const Patch &patch) { return mean(&patch.stand.get_gridcell().climate.vpds.front(), patch.stand.get_gridcell().climate.vpds.size());     }, "Daily VPD (kPa)");
	add_patch_output("file_dave_met_insol"           , DAILY, "insol", [](const Patch &patch) { return patch.stand.get_gridcell().climate.insol;   }, "Daily insolation (units depend on instype)");
	add_patch_output("file_dave_met_precip"          , DAILY, "precip", [](const Patch &patch) { return patch.stand.get_gridcell().climate.prec;    }, "Daily total precipitation (mm)");
	add_patch_output("file_dave_met_pressure"        , DAILY, "pressure", [](const Patch &patch) { return mean(&patch.stand.get_gridcell().climate.pressure.front(), patch.stand.get_gridcell().climate.pressure.size()); }, "Daily air pressure (Pa)");
	add_patch_output("file_dave_met_co2"             , DAILY, "co2", [](const Patch &patch) { return patch.stand.get_gridcell().climate.co2;         }, "Daily atmospheric co2 concentration (ppm)");
	add_patch_output("file_dave_daylength"           , DAILY, "daylength", [](const Patch &patch) { return patch.stand.get_gridcell().climate.daylength;         }, "Daylength (h)");
	add_patch_output("file_dave_sompool_cmass"       , DAILY, sompool_cmass_columns, "Daily SOM pool C mass (kgC/m2)");
	add_patch_output("file_dave_sompool_nmass"       , DAILY, sompool_nmass_columns, "Daily SOM pool N mass (kgN/m2)");
	add_patch_output("file_dave_ninput"              , DAILY, ninput_columns, "Daily nitrogen deposition (kgN/m2)");
	add_patch_output("file_dave_fpar_ff"             , DAILY, "fpar_ff", [](const Patch &p) { return p.fpar_ff; }, "Forest floor fpar (0-1)");

	// Daily per-patch aggregated individual-level outputs.
	add_aggregated_indiv_output("file_dave_lai"                 , DAILY, [](const Individual &i) { return i.lai; }, "LAI by individual (m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_sla"                 , DAILY, [](const Individual &i) { return i.sla(); }, "SLA by individual (m2/kgC)", mean, false, true);
	add_aggregated_indiv_output("file_dave_cmass_leaf"          , DAILY, [](const Individual &i) { return i.get_cmass_leaf(); }, "Leaf C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_leaf_brown"    , DAILY, [](const Individual &i) { return i.cmass_leafs[1]; }, "Leaf brown C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_leaf"          , DAILY, [](const Individual &i) { return i.get_nmass_leaf(); }, "Leaf N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_crown"         , DAILY, [](const Individual &i) { return i.cmass_crown; }, "Crown C mass (kgC/m2)", mean, false, true);
	add_aggregated_indiv_output("file_dave_cmass_repr"          , DAILY, [](const Individual &i) { return i.cmass_repr; }, "Reproductive C mass (kgC/m2)", mean, false, true);
	add_aggregated_indiv_output("file_dave_cmass_root"          , DAILY, [](const Individual &i) { return i.cmass_root; }, "Root C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass"               , DAILY, [](const Individual &i) { return i.ccont(); }, "Total C biomass per individual per patch (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass"               , DAILY, [](const Individual &i) { return i.ncont(); }, "Total N biomass per individual per patch (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage"       , DAILY, [](const Individual &i) { return i.cmass_storage; }, "Storage C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage_max"   , DAILY, [](const Individual &i) { return i.cstorage_max; }, "Storage capacity for C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_storage"       , DAILY, [](const Individual &i) { return i.nmass_storage; }, "Storage N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_storage_max"   , DAILY, [](const Individual &i) { return i.nstorage_max; }, "Storage capacity for N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_sap"           , DAILY, [](const Individual &i) { return i.cmass_sap; }, "Sapwood C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_sap"           , DAILY, [](const Individual &i) { return i.nmass_sap; }, "Sapwood N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_heart"         , DAILY, [](const Individual &i) { return i.cmass_heart; }, "Heartwood C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_heart"         , DAILY, [](const Individual &i) { return i.nmass_heart; }, "Heartwood N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_repr"          , DAILY, [](const Individual &i) { return i.nmass_repr; }, "Reproductive N mass (kgN/m2)", mean, false, true);
	add_aggregated_indiv_output("file_dave_ndemand"             , DAILY, [](const Individual &i) { return i.ndemand; }, "N demand mass (kgN/m2)", mean, false, true);
	add_aggregated_indiv_output("file_dave_height"              , DAILY, [](const Individual &i) { return i.height; }, "Plant height (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_diameter"            , DAILY, [](const Individual &i) { return i.diam; }, "Stem diameter (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_density"             , DAILY, [](const Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_sapwood_area"        , DAILY, [](const Individual &i) { return i.sa; }, "sapwood area sum of individuals (m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_latosa"              , DAILY, [](const Individual &i) { return i.latosa; }, "leaf area to sapwood area ratio (unitless)", mean, true);
	add_aggregated_indiv_output("file_dave_fpar"                , DAILY, [](const Individual &i) { return i.fpar; }, "Daily fpar [0-1]", sum, true);
	add_aggregated_indiv_output("file_dave_indiv_gpp"           , DAILY, [](const Individual &i) { return i.dgpp * G_PER_KG; }, "Daily GPP (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_indiv_resp_auto"     , DAILY, [](const Individual &i) { return (i.dgpp - i.dnpp) * G_PER_KG; }, "Daily autotrophic respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_indiv_npp"           , DAILY, [](const Individual &i) { return i.dnpp * G_PER_KG; }, "Daily NPP (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_fpc"                 , DAILY, [](const Individual &i) { return i.fpc; }, "Daily FPC ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_crownarea"           , DAILY, [](const Individual &i) { return i.crownarea; }, "Daily Crown Area (m2)", sum, true);
	add_aggregated_indiv_output("file_dave_agd_g"               , DAILY, [](const Individual &i) { return i.photosynthesis.agd_g; }, "Daily gross photosynthesis (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_rd_g"                , DAILY, [](const Individual &i) { return i.photosynthesis.rd_g; }, "Daily leaf respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_je"                  , DAILY, [](const Individual &i) { return i.photosynthesis.je; }, "Daily PAR-limited photosynthesis rate (gC/m2/h)", sum, true);
	add_aggregated_indiv_output("file_dave_vm"                  , DAILY, [](const Individual &i) { return i.photosynthesis.vm; }, "RuBisCO capacity (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_fphen_activity"      , DAILY, [](const Individual &i) { return i.fphen_activity; }, "Dormancy downregulation (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_fdev_growth"         , DAILY, [](const Individual &i) { return i.fdev_growth; }, "Development factor for growth demand (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_frepr_cstruct"       , DAILY, [](const Individual &i) { return i.frepr_cstruct; }, "Ratio of reproductive to aboveground structural biomass", mean, false, true);
	add_aggregated_indiv_output("file_dave_growth_demand"       , DAILY, [](const Individual &i) { return i.growth_demand; }, "Growth Demand (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_et"                  , DAILY, [](const Individual &i) { return i.aet; }, "Daily individual-level evapotranspiration (mm/day)", sum, true, true);
	add_aggregated_indiv_output("file_dave_nscal"               , DAILY, [](const Individual &i) { return i.nscal_5.lastadd(); }, "N stress scalar (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_nscal_mean"          , DAILY, [](const Individual &i) { return i.nscal_5.mean(); }, "N stress scalar 5 day running mean (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_ltor"            	, DAILY, [](const Individual &i) { return i.ltor; }, "Daily leaf to root ratio (unitless)", mean, true);
	add_aggregated_indiv_output("file_dave_cue"                 , DAILY, [](const Individual &i) { return i.dnpp / i.dgpp; }, "carbon use efficiency ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_alpha_leaf"          , DAILY, [](const Individual &i) { return i.alpha_leaf; }, "Daily leaf sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_root"          , DAILY, [](const Individual &i) { return i.alpha_root; }, "Daily root sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_sap"           , DAILY, [](const Individual &i) { return i.alpha_sap; }, "Daily sap sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_repr"          , DAILY, [](const Individual &i) { return i.alpha_repr; }, "Daily repr sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cmass_leaf_limit"    , DAILY, [](const Individual &i) { return i.cmass_leaf_limit; }, "Daily optimum leaf pool size for trees (kgC/m2)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cmass_root_limit"    , DAILY, [](const Individual &i) { return i.cmass_root_limit; }, "Daily optimum root pool size for trees (kgC/m2)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cmass_sap_limit"     , DAILY, [](const Individual &i) { return i.cmass_sap; }, "Daily optimum sap pool size for trees (kgC/m2)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cmass_repr_limit"    , DAILY, [](const Individual &i) { return i.cmass_repr_limit; }, "Daily optimum reproductive pool size for trees (kgC/m2)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cmass_storage_limit" , DAILY, [](const Individual &i) { return i.cmass_storage_limit; }, "Daily optimum storage pool size (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_leaf"          , DAILY, [](const Individual &i) { return i.cgrow_leaf; }, "Daily leaf growth (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_root"          , DAILY, [](const Individual &i) { return i.cgrow_root; }, "Daily root growth (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_sap"           , DAILY, [](const Individual &i) { return i.cgrow_sap; }, "Daily sap growth (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_repr"          , DAILY, [](const Individual &i) { return i.cgrow_repr; }, "Daily repr growth (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_diameter_inc"        , DAILY, [](const Individual &i) { return i.d_diam; }, "Daily diameter growth (m/day)", mean, true, false);
	add_aggregated_indiv_output("file_dave_height_inc"          , DAILY, [](const Individual &i) { return i.d_height; }, "Daily height growth (m/day)", mean, true, false);
	add_aggregated_indiv_output("file_dave_dturnover_leaf"      , DAILY, [](const Individual &i) { return i.cturnover_leaf; }, "daily leaf C turnover flux (kgC/m2/day)", mean, false, false);
	add_aggregated_indiv_output("file_dave_dturnover_root"      , DAILY, [](const Individual &i) { return i.cturnover_root; }, "daily root C turnover flux (kgC/m2/day)", mean, false, false);
	add_aggregated_indiv_output("file_dave_dturnover_sap"       , DAILY, [](const Individual &i) { return i.cturnover_sap; }, "daily sapwood C turnover flux (kgC/m2/day)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anc_frac"            , DAILY, [](const Individual &i) { return i.anc_frac; }, "Daily anc as fraction of total daily photosynthesis (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anj_frac"            , DAILY, [](const Individual &i) { return i.anj_frac; }, "Daily anj as fraction of total daily photosynthesis (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anp_frac"            , DAILY, [](const Individual &i) { return i.anp_frac; }, "Daily anp as fraction of total daily photosynthesis (0-1)", mean, false, false);

	// Daily individual-level outputs.
	// These will produce one row per individual per patch per timestep.
	add_indiv_output("file_dave_indiv_cpool"                    , DAILY, indiv_cpool_columns, "Indiv-level carbon pools (kgC/m2)", mean, false, true);
	add_indiv_output("file_dave_indiv_npool"                    , DAILY, indiv_npool_columns, "Indiv-level nitrogen pools (kgN/m2)", mean, false, true);

	// Daily individual-level array outputs.
	add_indiv_array_output("file_dave_layerwise_fpar", DAILY, "layer", get_nlayer_fpar, [](const Individual &indiv, int i) { return indiv.fipar_layerwise[i]; }, "Daily layered fpar [0-1]");
	add_indiv_array_output("file_dave_layerwise_lai" , DAILY, "layer", get_nlayer_lai , [](const Individual &indiv, int i) { return indiv.lai_layerwise[i];  }, "Daily layered LAI (m2/m2)");

	// Daily pft-level per-patch outputs.
	add_pft_output("file_dave_wscal"                 , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].wscal; }, "Water scalar [0-1], 1 = min stress",  false, true);
	add_pft_output("file_dave_cmass_litter_repr"     , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_repr; }, "Seedbank cmass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_litter_repr"     , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_repr; }, "Seedbank nmass (kgN/m2)",  true);

	// Subdaily individual-level outputs.
	//                      param name in .ins file    timestep  column name    function which gets number of array elems     function which returns the i-th output value                               Description of output file
	add_indiv_array_output("file_dave_subdaily_an"   , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.An[i];    }, "Subdaily net photosynthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_rd"   , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Rd[i];    }, "Subdaily leaf respiration (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_anc"  , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anc[i];   }, "Subdaily RuBisCO-limited net photosynthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_anj"  , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anj[i];   }, "Subdaily RuBP regeneration-limited photoysnthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_gsw"  , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.gsw[i];   }, "Subdaily stomatal conductance to water vapour (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_ci"   , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.ci[i];    }, "Subdaily intercellular CO2 concentration (mol/mol)");
	add_indiv_array_output("file_dave_subdaily_vcmax", DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.vcmax[i]; }, "Subdaily max carboxylation rate (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_jmax" , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.jmax[i];  }, "Subdaily max electron transport rate (mol/m2/s)");

	// Daily patch-level array outputs.
	add_patch_array_output("file_dave_sw"            , DAILY, "sw_", NSOILLAYER, [](const Patch &p, int i) { return p.soil.get_layer_soil_water(i); }, "Soil water content as fraction of available water holding capacity", true, false);
	add_patch_array_output("file_dave_swmm"          , DAILY, "swmm_", NSOILLAYER, [](const Patch &p, int i) { return p.soil.get_layer_soil_water(i) * p.soil.soiltype.awc[i] + p.soil.soiltype.wp[i]; }, "Soil water content (mm)", true, false);
	add_patch_array_output("file_dave_swvol"         , DAILY, "swvol_", NSOILLAYER, [](const Patch &p, int i) { return (p.soil.get_layer_soil_water(i) * p.soil.soiltype.awc[i] + p.soil.soiltype.wp[i]) / Dz_soil; }, "Volumetric soil water content (m3/m3)", true, false);
	// add_daily_patch_level_array_output("file_dave_cfluxes_patch", get_patch_flux_names(), (int)Fluxes::NPERPATCHFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPatchFluxType)i, date.day) * G_PER_KG; }, "Daily patch-level fluxes (g/m2/day)");
	// add_daily_patch_level_array_output("file_dave_cfluxes_pft", get_pft_flux_names(), (int)Fluxes::NPERPFTFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPFTFluxType)i, date.day) * G_PER_KG; }, "Daily PFT-level fluxes (g/m2/day)");
	add_patch_array_output("file_dave_resp"          , DAILY, "resp", 1, get_er, "Daily patch-level ecosystem respiration (gC/m2/day)", true);
	add_patch_array_output("file_dave_gpp"           , DAILY, "gpp", 1, get_gpp, "Daily patch-level GPP (gC/m2/day)", true);
	add_patch_array_output("file_dave_npp"           , DAILY, "npp", 1, get_npp, "Daily patch-level NPP (gC/m2/day)", true);
	add_patch_array_output("file_dave_nee"           , DAILY, "nee", 1, get_nee, "Daily patch-level NEE (gC/m2/day)", true);

	// Subdaily patch-level outputs.
	//                                               column name              function which gets number of array elems     function which returns the i-th output value
	add_patch_array_output("file_dave_met_subdaily_temp"    , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.temps[i];    }, "Subdaily temperature (°C)", false, true);
	add_patch_array_output("file_dave_met_subdaily_par"     , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pars[i];     }, "Subdaily PAR (J/m2/timestep)", true, false);
	add_patch_array_output("file_dave_met_subdaily_vpd"     , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.vpds[i];     }, "Subdaily VPD (kPa)", false, true);
	add_patch_array_output("file_dave_met_subdaily_insol"   , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.insols[i];   }, "Subdaily insolation (units depend on instype)", true, false);
	add_patch_array_output("file_dave_met_subdaily_precip"  , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.precs[i];    }, "Subdaily precipitation (mm)", true, false);
	add_patch_array_output("file_dave_met_subdaily_pressure", DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pressure[i]; }, "Subdaily air pressure (Pa)", false, true);
	add_patch_array_output("file_dave_met_subdaily_co2"     , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.co2;         }, "Subdaily atmospheric co2 concentration (ppm)", false, true);

	// Standpft-level outputs.
	//add_spft_output("file_dave_cmass_repr"           , DAILY, [](const Standpft &spft) { return spft.cmass_repr; }, "Daily reproductive cmass at standpft-level (gC/m2/day)");

	// Register instruction file parameters.
	declare_parameters();

	total_time = get_wall_time() - t_start;
}

DaveOutput::~DaveOutput() {
	printf("Total time spent in dave output: %.2fs\n", total_time);
}

void DaveOutput::init() {
	double t_start = get_wall_time();
	open_outfiles();
	total_time += get_wall_time() - t_start;
}

void DaveOutput::write_outputs(Gridcell &gridcell, std::vector<std::unique_ptr<OutputFile>> &outputs) {
	ensure_channel_open();

	for (std::unique_ptr<OutputFile> &output : outputs) {
		if (output->enabled()) {
			output->write_rows(*output_channel, gridcell);
		}
	}
}

bool should_write_today() {
	return spinup_output || date.year >= nyear_spinup;
}

bool write_stats() {
	return patch_output_stats && npatch > 1;
}

void DaveOutput::outannual(Gridcell &gridcell) {
	double t0 = get_wall_time();
	if (should_write_today()) {
		write_outputs(gridcell, annual_outputs);
	}
	total_time += get_wall_time() - t0;
}

void DaveOutput::outdaily(Gridcell &gridcell) {
	double t0 = get_wall_time();
	if (should_write_today()) {
		write_outputs(gridcell, daily_outputs);
	}
	total_time += get_wall_time() - t0;
}

void DaveOutput::add_patch_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		double (*read_output)(const Patch&),
		std::string description,
		bool sum,
		bool mean) {
	std::vector<std::string> colnames;
	colnames.push_back(colname);
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colnames, 1, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_output(
		std::string name,
		timestep_t timestep,
		std::vector<ColumnDefinition<Patch>> (*get_columns)(),
		std::string description,
		bool sum,
		bool mean) {
	CustomPatchLevelOutput *output = new CustomPatchLevelOutput(name, get_columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_stand_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		double (*read_output)(const Stand &),
		std::string description,
		bool sum,
		bool mean) {
	ColumnDefinition<Stand> column(colname, read_output);
	std::vector<ColumnDefinition<Stand>> columns;
	columns.push_back(column);
	StandLevelOutput *output = new StandLevelOutput(name, columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_stand_output(
	std::string name,
	timestep_t timestep,
	std::vector<ColumnDefinition<Stand>> (*get_columns)(),
	std::string description,
	bool sum,
	bool mean) {
	StandLevelOutput *output = new StandLevelOutput(name, get_columns, description, sum, mean);
}

void DaveOutput::add_patch_array_output(
		std::string name,
		timestep_t timestep,
		std::vector<std::string> colnames,
		double (*read_output)(const Patch&, int),
		std::string description,
		bool sum,
		bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colnames, colnames.size(), reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_array_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		int nelem,
		double (*read_output)(const Patch&, int),
		std::string description,
		bool sum,
		bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, nelem, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_array_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		int (*get_nelem)(),
		double (*read_output)(const Patch&, int),
		std::string description,
		bool sum,
		bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, get_nelem, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_aggregated_indiv_output(std::string name, timestep_t timestep, double (*read_output)(const Individual&), std::string description, double (*aggregator)(double*, int), bool sum, bool mean) {
	OutputReader<Individual> *reader = new OutputReader<Individual>(read_output);
	auto output = new AggregatedIndivOutput(name, reader, aggregator, description, sum, mean);

	register_output(output, timestep);
}

void DaveOutput::add_indiv_output(
		std::string name,
		timestep_t timestep,
		std::vector<ColumnDefinition<Individual>> (*get_columns)(),
		std::string description,
		double (*aggregator)(double *, int),
		bool sum,
		bool mean) {
	IndividualOutput *output = new IndividualOutput(name, get_columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_pft_output(std::string name, timestep_t timestep, double (*read_output)(const Patch &, const Pft&), std::string description, bool sum, bool mean) {
	PftOutput *output = new PftOutput(name, read_output, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_array_output(
		std::string name,
		timestep_t timestep,
		std::string elem_name,
		int (*get_nelem)(const Patch &),
		double (*read_output)(const Individual&, int),
		std::string description) {
	ArrayOutputReader<Individual> *reader = new ArrayOutputReader<Individual>(read_output);
	IndivLevelArrayOutput *output = new IndivLevelArrayOutput(name, elem_name, get_nelem, reader, description);
	register_output(output, timestep);
}

void DaveOutput::add_spft_output(
		std::string name,
		timestep_t timestep,
		double (*read_output)(const Standpft&),
		std::string description,
		bool sum,
		bool mean) {
	OutputReader<Standpft> *reader = new OutputReader<Standpft>(read_output);
	StandPftOutput *output = new StandPftOutput(name, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::register_output(OutputFile *file, timestep_t timestep) {
	auto &vec = get_outputs_vector(timestep);
	// vec.push_back(file);
	vec.emplace_back(file);
	// vec.push_back(std::move(file));
}

std::vector<std::unique_ptr<DaveOutput::OutputFile>> &DaveOutput::get_outputs_vector(timestep_t timestep) {
	if (timestep == ANNUAL) {
		return annual_outputs;
	}
	return daily_outputs;
}

void DaveOutput::declare_parameters() {
	declare_parameter("spinup_output", &spinup_output, "True to write DAVE outputs during the spinup, false otherwise.");
	declare_parameter("patch_output_stats", &patch_output_stats, "True to ONLY write stats (mean/std. deviation) for patch-level outputs. False to write the full patch-level outputs.");

	declare_parameters(daily_outputs);
	declare_parameters(annual_outputs);
}

void DaveOutput::declare_parameters(std::vector<std::unique_ptr<OutputFile>> &outputs) {
	for (auto &file : outputs) {
		file->declare_parameters();
	}
}

void DaveOutput::open_outfiles() {
	open_outfiles(daily_outputs);
	open_outfiles(annual_outputs);
}

void DaveOutput::open_outfiles(std::vector<std::unique_ptr<OutputFile>> &outputs) {
	ensure_channel_open();

	for (std::unique_ptr<OutputFile> &output : outputs) {
		if (output->enabled()) {
			output->open(*output_channel);
		}
	}
}

// #############################################################################
// DaveOutput::OutputFile implementation.
// #############################################################################

DaveOutput::OutputFile::~OutputFile() {
	// table will be NULL if this output file is disabled.
	if (table) {
		delete table;
	}
	table = NULL;
}

DaveOutput::OutputFile::OutputFile(std::string param, std::string desc) {
	param_name = param;
	description = desc;
	table = NULL;
}

void DaveOutput::OutputFile::declare_parameters() {
	declare_parameter(param_name.c_str(), &file_name, PATH_MAX, description.c_str());
}

void DaveOutput::OutputFile::open(OutputChannel &channel) {
	const char *name = file_name.c_str();
	Table tbl = channel.create_table(TableDescriptor(name, get_cols()));
	table = new Table(tbl.id());
}

void DaveOutput::OutputFile::close(OutputChannel &channel) {
	assert_table_open();
	channel.close_table(*table);
}

bool DaveOutput::OutputFile::should_write(OutputRows &rows) {
	if (!should_write_today()) {
		return false;
	}
	return true;
}

void DaveOutput::OutputFile::write_output(OutputRows &rows, double value)
{
	assert_table_open();
	if (should_write(rows)) {
		rows.add_value(*table, value);
	}
}

void DaveOutput::OutputFile::write_output(OutputRows &rows, int value)
{
	assert_table_open();
	if (should_write(rows)) {
		rows.add_value(*table, value);
	}
}

void DaveOutput::OutputFile::write_output(OutputRows &rows, const char *value)
{
	assert_table_open();
	if (should_write(rows)) {
		rows.add_value(*table, value);
	}
}

void DaveOutput::OutputFile::finish_row(OutputChannel &rows, double lon, double lat, int year, int day) {
	assert_table_open();
	rows.finish_row(*table, lon, lat, year, day);
}

void DaveOutput::OutputFile::assert_table_open() {
	if (!table) {
		fail("Output table for '%s' is not open", param_name.c_str());
	}
}

bool DaveOutput::OutputFile::enabled() {
	// Output file is enabled iff output file name is not empty.
	return !file_name.empty();
}

// #############################################################################
// A few helper functions.
// #############################################################################

void add_column_with_suffix(ColumnDescriptors &cols, const char *colname, const char *sfx, int width, int precision) {
	std::stringstream ss;
	ss << colname << "_" << sfx;
	std::string sname = ss.str();
	const char *name = sname.c_str();
	width = max((size_t)width, sname.length() + 1);
	cols += ColumnDescriptor(name, width, precision);
}

void add_column(ColumnDescriptors &cols, const char *colname, int width, int precision) {
	if (write_stats()) {
		add_column_with_suffix(cols, colname, "mean", width, precision);
		add_column_with_suffix(cols, colname, "median", width, precision);
		add_column_with_suffix(cols, colname, "sd", width, precision);
	} else {
		cols += ColumnDescriptor(colname, width, precision);
	}
}

void add_columns(ColumnDescriptors &cols, std::vector<std::string> names, int width, int precision) {
	for (std::string name : names) {
		add_column(cols, name.c_str(), width, precision);
	}
}

ColumnDescriptors create_column_descriptors() {
	ColumnDescriptors columns;
	columns += ColumnDescriptor(COLNAME_STAND, COLUMN_WIDTH, 0);
	if (!write_stats()) {
		columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	}
	return columns;
}

// #############################################################################
// StandLevelOutput implementation.
// #############################################################################

DaveOutput::StandLevelOutput::~StandLevelOutput() { }

DaveOutput::StandLevelOutput::StandLevelOutput(
		std::string param,
		std::vector<ColumnDefinition<Stand>> (*get_columns)(),
		std::string desc,
		bool sum,
		bool mean) : OutputFile(param, desc) {
	get_column_definitions = get_columns;
	write_sum = sum;
	write_mean = mean;
}

DaveOutput::StandLevelOutput::StandLevelOutput(
		std::string param,
		std::vector<ColumnDefinition<Stand>> columns,
		std::string desc,
		bool sum,
		bool mean) : OutputFile(param, desc) {
	get_column_definitions = NULL;
	this->columns = columns;
	write_sum = sum;
	write_mean = mean;
}

void DaveOutput::StandLevelOutput::write_rows(
		OutputChannel &channel,
		Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	// Each vector of this vector represents a single stand. Each element in the
	// sub-vectors represents a value which will be written for column i of that
	// particular stand. This is only really required for write_stats mode.
	std::vector<std::vector<double>> data;

	for (auto gc_itr = gridcell.begin(); gc_itr != gridcell.end(); gc_itr++) {
		Stand &stand = *gc_itr;
		std::vector<double> stand_data = write_stand(stand);
		if (write_stats()) {
			data.push_back(stand_data);
		} else {
			// The code would be cleaner if we were to use the 2d data array in
			// both cases, but I suspect that performance will be better with
			// this shortcut when the number of stands is high.

			// Write stand ID.
			// todo: would it be better to write the stand name instead?
			write_output(rows, stand.id);

			// Write stand data.
			double total = 0;
			for (double x : stand_data) {
				write_output(rows, x);
				total += x;
			}
			if (write_sum) {
				write_output(rows, total);
			}
			if (write_mean) {
				write_output(rows, total / stand_data.size());
			}
			finish_row(channel, lon, lat, year, day);
		}
	}

	if (write_stats() && data.size()) {
		// Each element of data contains all columns for each stand. We want
		// to write one column at a time.

		// First, get the number of columns.
		int ncol = data[0].size();

		// Now iterate over the columns that we want to write.
		for (size_t j = 0; j < ncol; j++) {
			// Get this column's value for each stand.
			std::vector<double> column_data;
			double total = 0;
			for (size_t i = 0; i < data.size(); i++) {
				column_data.push_back(data[i][j]);
			}
			// Compute summary statistics for this column.
			double mean = 0, sd = 0, median = 0;
			compute_stats(column_data, &mean, &sd, &median);
			write_output(rows, mean);
			write_output(rows, median);
			write_output(rows, sd);
		}
		if (write_sum || write_mean) {
			// Calculate the total for each stand.
			std::vector<double> totals;
			for (size_t i = 0; i < data.size(); i++) {
				// Calculate the total for the i-th stand.
				double total = 0;
				for (size_t j = 0; j < data[i].size(); i++) {
					total += data[i][j];
				}
				totals.push_back(total);
			}
			double mean = 0, sd = 0, median = 0;
			compute_stats(totals, &mean, &sd, &median);
			if (write_sum) {
				write_output(rows, mean);
				write_output(rows, median);
				write_output(rows, sd);
			}
			if (write_mean) {
				write_output(rows, mean / ncol);
				write_output(rows, median / ncol);
				write_output(rows, sd / ncol);
			}
		}
		finish_row(channel, lon, lat, year, day);
	}
}

ColumnDescriptors DaveOutput::StandLevelOutput::get_cols() {
	ColumnDescriptors descriptors;
	if (!write_stats()) {
		add_column(descriptors, COLNAME_STAND, COLUMN_WIDTH, 0);
	}

	// Only invoke get_column_definitions if necessary.
	if (get_column_definitions) {
		columns = get_column_definitions();
	}

	for (const ColumnDefinition<Stand> &column : columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(descriptors, name.c_str(), width, NUMERIC_PRECISION);
	}

	if (write_sum) {
		add_column(descriptors, COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	if (write_mean) {
		add_column(descriptors, COLNAME_MEAN, COLUMN_WIDTH, NUMERIC_PRECISION);
	}

	return descriptors;
}

std::vector<double> DaveOutput::StandLevelOutput::write_stand(
		const Stand &stand) {
	std::vector<double> data;
	for (ColumnDefinition<Stand> &column : columns) {
		double value = column.read_value(stand);
		data.push_back(value);
	}
	return data;
}

// #############################################################################
// PatchLevelOutput implementation.
// #############################################################################

DaveOutput::PatchLevelOutput::~PatchLevelOutput() { }

DaveOutput::PatchLevelOutput::PatchLevelOutput(std::string param, std::string desc, bool sum, bool mean)
: OutputFile(param, desc) {
	write_sum = sum;
	write_mean = mean;
}

void DaveOutput::PatchLevelOutput::write_rows(OutputChannel &channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;

		stand.firstobj();

		std::vector<std::vector<double>> data;

		while (stand.isobj) {
			Patch &patch = stand.getobj();

			// Write stand ID.
			write_output(rows, stand.id);

			std::vector<double> patch_data = write_patch(patch, rows);
			if (write_stats()) {
				data.push_back(patch_data);
			} else {
				// The code would be cleaner if we were to use the data 2d array
				// in both cases, but I suspect that performance will be better
				// with this shortcut when the number of patches is high.

				// Write patch ID.
				write_output(rows, patch.id);

				// Write patch data.
				double total = 0;
				for (double x : patch_data) {
					write_output(rows, x);
					total += x;
				}
				if (write_sum) {
					write_output(rows, total);
				}
				if (write_mean) {
					write_output(rows, total / patch_data.size());
				}
			}

			if (!write_stats()) {
				finish_row(channel, lon, lat, year, day);
			}

			stand.nextobj();
		}

		if (write_stats() && data.size()) {
			// Each element of data contains all columns for each patch. We want
			// to write one column at a time.

			// First, get the number of columns.
			int ncol = data[0].size();

			// Now iterate over the columns that we want to write.
			for (size_t j = 0; j < ncol; j++) {
				// Get this column's value for each patch.
				std::vector<double> column_data;
				double total = 0;
				for (size_t i = 0; i < data.size(); i++) {
					column_data.push_back(data[i][j]);
				}
				// Compute summary statistics for this column.
				double mean = 0, sd = 0, median = 0;
				compute_stats(column_data, &mean, &sd, &median);
				write_output(rows, mean);
				write_output(rows, median);
				write_output(rows, sd);

			}
			if (write_sum || write_mean) {
				std::vector<double> totals;
				// Iterate over patches.
				for (int i = 0; i < data.size(); i++) {
					// Get the total for this patch.
					double total = 0;
					for (size_t j = 0; j < data[i].size(); j++) {
						total += data[i][j];
					}
					totals.push_back(total);
				}
				double mean = 0, sd = 0, median = 0;
				compute_stats(totals, &mean, &sd, &median);
				if (write_sum) {
					write_output(rows, mean);
					write_output(rows, median);
					write_output(rows, sd);
				}
				if (write_mean) {
					write_output(rows, mean / ncol);
					write_output(rows, median / ncol);
					write_output(rows, sd / ncol);
				}
			}
			finish_row(channel, lon, lat, year, day);
		}

		gc_itr++;
	}
}

ColumnDescriptors DaveOutput::PatchLevelOutput::get_cols() {
	ColumnDescriptors columns = get_columns();
	if (write_sum) {
		add_column(columns, COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	if (write_mean) {
		add_column(columns, COLNAME_MEAN, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	return columns;
}

// #############################################################################
// IndividualOutput implementation.
// #############################################################################

DaveOutput::IndividualOutput::IndividualOutput(
		std::string param_name,
		std::vector<ColumnDefinition<Individual>> (*get_columns)(),
		std::string description,
		bool sum,
		bool mean) : OutputFile(param_name, description) {
	get_column_definitions = get_columns;
	write_sum = sum;
	write_mean = mean;
}

DaveOutput::IndividualOutput::~IndividualOutput() { }

void DaveOutput::IndividualOutput::write_rows(OutputChannel &channel, Gridcell &gridcell) {

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;

		stand.firstobj();

		std::vector<std::vector<double>> data;

		while (stand.isobj) {
			const Patch &patch = stand.getobj();
			const Vegetation &vegetation = patch.vegetation;

			for (size_t i = 0; i < vegetation.nobj; i++) {
				const Individual &indiv = vegetation[i];

				// Write stand ID.
				write_output(rows, stand.id);

				// Write patch ID.
				write_output(rows, patch.id);

				// Write individual ID.
				write_output(rows, indiv.id);

				// Write PFT name.
				write_output(rows, indiv.pft.name);

				// Write data.
				for (size_t j = 0; j < columns.size(); j++) {
					write_output(rows, columns[j].read_value(indiv));
				}

				// End of this row.
				finish_row(channel, lon, lat, year, day);
			}

			stand.nextobj();
		}

		gc_itr++;
	}
}

ColumnDescriptors DaveOutput::IndividualOutput::get_cols() {
	columns = get_column_definitions();

	ColumnDescriptors descriptors;
	descriptors += ColumnDescriptor(COLNAME_STAND, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_INDIV_ID, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_PFT, COLUMN_WIDTH + 2, 0, TEXT);

	for (const ColumnDefinition<Individual> &column : columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(descriptors, name.c_str(), width, NUMERIC_PRECISION);
	}

	return descriptors;
}

// #############################################################################
// PftOutput implementation.
// #############################################################################

DaveOutput::PftOutput::PftOutput(
	std::string param_name,
	double (*get_output)(const Patch &, const Pft &),
	std::string description,
	bool sum,
	bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	read_output = get_output;
}

ColumnDescriptors DaveOutput::PftOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	std::vector<std::string> pfts = get_pfts(true);
	add_columns(columns, pfts, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

std::vector<double> DaveOutput::PftOutput::write_patch(Patch &patch, OutputRows &out) {
	std::vector<double> results;

	for (int i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];

		double result = get_value(patch, pft);
		results.push_back(result);
	}

	return results;
}

double DaveOutput::PftOutput::get_value(const Patch &patch, const Pft &pft) const {
	return read_output(patch, pft);
}

DaveOutput::PftOutput::~PftOutput() { }

// #############################################################################
// CustomPatchLevelOutput implementation.
// #############################################################################

DaveOutput::CustomPatchLevelOutput::CustomPatchLevelOutput(
	std::string param_name,
	std::vector<ColumnDefinition<Patch>> (*get_columns)(),
	std::string description,
	bool sum,
	bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
		get_column_definitions = get_columns;
}

std::vector<double> DaveOutput::CustomPatchLevelOutput::write_patch(Patch &patch, OutputRows &out) {
	std::vector<double> result;

	for (const ColumnDefinition<Patch> &col : columns) {
		result.push_back(col.read_value(patch));
	}

	return result;
}

ColumnDescriptors DaveOutput::CustomPatchLevelOutput::get_columns() {
	this->columns = this->get_column_definitions();

	ColumnDescriptors columns = create_column_descriptors();

	for (const ColumnDefinition<Patch> &column : this->columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(columns, name.c_str(), width, NUMERIC_PRECISION);
	}

	return columns;
}

DaveOutput::CustomPatchLevelOutput::~CustomPatchLevelOutput() { }

// #############################################################################
// StandPftOutput implementation.
// #############################################################################

DaveOutput::StandPftOutput::StandPftOutput(
		std::string param_name,
		OutputReader<Standpft> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	this->read_output = std::unique_ptr<OutputReader<Standpft>>(get_output);
}

DaveOutput::StandPftOutput::~StandPftOutput() { }

std::vector<double> DaveOutput::StandPftOutput::write_patch(Patch &patch, OutputRows &out) {
	Vegetation &vegetation = patch.vegetation;
	Stand &stand = patch.stand;

	std::vector<double> results;

	// Note that patch id is written before this function is called, so we
	// shouldn't write it ourselves (that would result in a 'too many columns'
	// error).

	// Write one column per PFT.
	for (int i = 0; i < pftlist.nobj; i++) {
		Pft &pft = pftlist[i];
		const Standpft &spft = stand.pft[pft.id];
		double result = read_output->read(spft);
		results.push_back(result);
	}

	return results;
}

ColumnDescriptors DaveOutput::StandPftOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	add_columns(columns, get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// PatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int n,
		ArrayOutputReader<Patch> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {

	nelem = n;
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
	get_nelem = NULL;

	for (int i = 0; i < n; i++) {
		if (n == 1) {
			column_names.push_back(col_name);
		} else {
			column_names.push_back(col_name + std::to_string(i));
		}
	}
}

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
	std::string param_name,
	std::vector<std::string> col_names,
	int n,
	ArrayOutputReader<Patch> *get_output,
	std::string description,
	bool sum,
	bool mean)
 : PatchLevelOutput(param_name, description, sum, mean) {
	column_names = col_names;
	nelem = n;
	get_nelem = NULL;
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
}
DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int (*get_nelem)(),
		ArrayOutputReader<Patch> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	this->get_nelem = get_nelem;
	column_names.push_back(col_name);
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
}

DaveOutput::PatchLevelArrayOutput::~PatchLevelArrayOutput() { }

std::vector<double> DaveOutput::PatchLevelArrayOutput::write_patch(
		Patch &patch, OutputRows &out) {
	std::vector<double> results;
	for (int i = 0; i < nelem; i++) {
		double value = read_output->read(patch, i);
		results.push_back(value);
	}
	return results;
}

ColumnDescriptors DaveOutput::PatchLevelArrayOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();

	if (get_nelem) {
		nelem = get_nelem();
		std::string basename = column_names[0];
		column_names.clear();
		for (int i = 0; i < nelem; i++) {
			if (nelem == 1) {
				column_names.push_back(basename);
			} else {
				column_names.push_back(basename + std::to_string(i));
			}
		}
	}

	for (int i = 0; i < nelem; i++) {
		std::string name = column_names[i];
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(columns, name.c_str(), width, NUMERIC_PRECISION);
	}
	return columns;
}

// #############################################################################
// IndivLevelArrayOutput implementation.
// #############################################################################

DaveOutput::IndivLevelArrayOutput::IndivLevelArrayOutput(
		std::string param_name,
		std::string element_name,
		int (*get_nelem)(const Patch &),
		ArrayOutputReader<Individual> *get_output,
		std::string description) : OutputFile(param_name, description)
{
	elem_name = element_name;
	get_n = get_nelem;
	read_elem = std::unique_ptr<ArrayOutputReader<Individual>>(get_output);
}

void DaveOutput::IndivLevelArrayOutput::write_rows(OutputChannel &channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	std::vector<std::string> pfts = get_pfts(false);

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;

		for (int j = 0; j < stand.nobj; j++) {
			Patch &patch = stand[j];

			int nelem = get_n(patch);
			Vegetation &vegetation = patch.vegetation;
			for (int i = 0; i < nelem; i++) {

				// Write patch ID.
				write_output(rows, patch.id);

				// Write array index.
				write_output(rows, i);

				vegetation.firstobj();
				double total = 0;

				for (int k = 0; k < pfts.size(); k++) {
					std::vector<Individual *> indivs = get_indivs(vegetation, pfts[k].c_str());
					double indiv_total = 0;
					for (Individual *indiv : indivs) {
						// Read array element for this individual.
						indiv_total += read_elem->read(*indiv, i);
					}
					total += indiv_total;

					// Write current array element for this individual.
					write_output(rows, indiv_total);
				}

				// Write total value across all individuals for this array elem.
				write_output(rows, total);

				// Finish the row, so next array element is written to next row.
				finish_row(channel, lon, lat, year, day);
			}
		}
		gc_itr++;
	}
}

ColumnDescriptors DaveOutput::IndivLevelArrayOutput::get_cols() {
	ColumnDescriptors columns;
	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	columns += ColumnDescriptor(elem_name.c_str(), COLUMN_WIDTH, 0);
	columns += ColumnDescriptors(get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// AggregatedIndivOutput implementation.
// #############################################################################

DaveOutput::AggregatedIndivOutput::AggregatedIndivOutput(
		std::string param_name,
		OutputReader<Individual> *get_output,
		double (*agg)(double *, int),
		std::string description,
		bool sum,
		bool mean) : PftOutput(param_name, NULL, description, sum, mean) {
	read_output = std::unique_ptr<OutputReader<Individual>>(get_output);
	aggregator = agg;
}

double DaveOutput::AggregatedIndivOutput::get_value(
		const Patch &patch,
		const Pft &pft) const {

	std::vector<double> data;
	const Vegetation &vegetation = patch.vegetation;
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		if (indiv.pft.id == pft.id) {
			double value = read_output->read(indiv);
			data.push_back(value);
		}
	}
	return aggregator(&data.front(), data.size());
}

bool spinup_output;
bool patch_output_stats;

const int COLUMN_WIDTH = 14;
const int COLUMN_WIDTH_PFT = 24;
const int NUMERIC_PRECISION = 8;
const char *COLNAME_TOTAL = "total";
const char *COLNAME_MEAN = "mean";
const char *COLNAME_PATCH = "patch";
const char *COLNAME_STAND = "stand";
const char *COLNAME_TIMESTEP = "timestep";

} // namespace GuessOutput
