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
#include "guessmath.h"
#include "gutil.h"
#include "outputchannel.h"
#include "outputmodule.h"
#include "parameters.h"
#include "util.h"
#include "logging.h"

namespace GuessOutput {

// Register the DAVE output module.
REGISTER_OUTPUT_MODULE("dave", DaveOutput)

// #############################################################################
// DaveOutput implementation.
// #############################################################################

DaveOutput::DaveOutput() {
	double t_start = get_wall_time();

	// Annual gridcell-level outputs.
	add_gridcell_aggregated_indiv_output("file_aagb"             , ANNUAL, [](const Individual &i) { return i.aboveground_biomass(); }, sum, "Annual gridcell-level above-ground biomass (kg/m2)", true);
	add_gridcell_aggregated_tree_output("file_aagb_tree"         , ANNUAL, [](const Individual &i) { return i.aboveground_biomass(); }, sum, "Annual gridcell-level above-ground tree biomass (kg/m2)", true);
	add_gridcell_output("file_dave_amet_prec"                    , ANNUAL, "prec", [](const Gridcell &g) { return g.climate.aprec; }, "Total annual precipitation (mm)");
	add_gridcell_output("file_dave_amet_par"                     , ANNUAL, "par", [](const Gridcell &g) { return g.climate.apar; }, "Total annual PAR (J/m2)");
	add_gridcell_aggregated_ppft_output("file_dave_lai_ymax"     , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].lai_ymax; }, "Maximum gridcell-level annual leaf area index (m2/m2)");
	add_gridcell_aggregated_ppft_output("file_dave_agb_ymax"     , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].agb_ymax; }, "Maximum gridcell-level annual above-ground biomass (kg/m2)");
	add_gridcell_aggregated_ppft_output("file_dave_lai_ymean"    , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].lai_ymean; }, "Mean gridcell-level annual leaf area index (m2/m2)");
	add_gridcell_aggregated_ppft_output("file_dave_agb_ymean"    , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].agb_ymean; }, "Mean gridcell-level annual above-ground biomass (kg/m2)");

	// Annual patch-level outputs.
	add_patch_output("file_dave_patch_age"                      , ANNUAL, "age", [](const Patch &p) { return 1.0 * p.age; }, "Annual patch-level age since last disturbance (years)");
	add_patch_output("file_dave_arunoff"                        , ANNUAL, annual_runoff_columns, "Annual runoff (mm)");
	add_patch_output("file_dave_acpool"                         , ANNUAL, cpool_columns, "Annual C pools (kgC/m2)", true);
	add_patch_output("file_dave_anpool"                         , ANNUAL, npool_columns, "Annual N pools (kgN/m2)");
	add_patch_output("file_dave_acflux"                         , ANNUAL, cflux_columns, "Annual C flux (kgC/m2)");
	add_patch_output("file_dave_anflux"                         , ANNUAL, nflux_columns, "Annual N flux (kgN/m2/year)");
	add_patch_output("file_dave_globfirm"                       , ANNUAL, globfirm_columns, "Annual GLOBFIRM output file");
	add_patch_output("file_dave_mwcont_upper"                   , ANNUAL, mwcont_upper_columns, "Monthly wcont_upper output file (0-1)");
	add_patch_output("file_dave_mwcont_lower"                   , ANNUAL, mwcont_lower_columns, "Monthly wcont_lower output file (0-1)");
	add_patch_output("file_dave_apet"                           , ANNUAL, pet_columns, "Annual potential patch-level evapotranspiration (mm)");
	add_patch_output("file_dave_asimfire"                       , ANNUAL, simfire_columns, "Annual simfire analysis");
	add_patch_output("file_dave_afuel"                          , ANNUAL, fuel_columns, "Annual blaze fuel availability (gC/m2)");
	add_patch_output("file_dave_acoarse_woody_debris"           , ANNUAL, coarse_woody_debris_columns, "Annual coarse woody debris (gC/m2)");
	add_patch_output("file_dave_amet_year"                      , ANNUAL, met_year_columns, "Current year of met data being used");
	add_patch_output("file_dave_aco2"                           , ANNUAL, "co2", [](const Patch &p) { return p.get_climate().co2; }, "Annual atmospheric co2 concentration (ppm)");
	add_patch_output("file_dave_aminleach"                      , ANNUAL, "aminleach", [](const Patch &patch) { return patch.soil.aminleach; }, "Leaching of soil mineral N (kgN/m2/yr)");
	add_patch_output("file_dave_sompool_acmass"                 , ANNUAL, sompool_cmass_columns, "Daily SOM pool C mass (kgC/m2)", true);
	add_patch_output("file_dave_sompool_anmass"                 , ANNUAL, sompool_nmass_columns, "Daily SOM pool N mass (kgN/m2)", true);
	add_patch_output("file_dave_andep"                          , ANNUAL, aninput_columns, "Annual nitrogen deposition (kgN/m2)");
	add_patch_output("file_dave_anfixation"                     , ANNUAL, "nfixation", [](const Patch &p) { return p.soil.anfix; }, "Total annual biological N fixation (kgN/m2)");
	add_patch_output("file_dave_water_budget_annual"            , ANNUAL, water_budget_cols_annual, "Annual water budget (mm)");
	add_patch_output("file_dave_nitrogen_budget_annual"         , ANNUAL, nitrogen_budget_cols_annual, "Annual nitrogen budget (kgN/m2)");
	add_patch_output("file_dave_anee"                           , ANNUAL, "nee", get_anee, "Total annual patch-level NEE (kgC/m2)");
	add_patch_output("file_dave_aclitter"                       , ANNUAL, clitter_columns, "Total End-of-Year Patch-Level C Mass in Litter (kgC/m2)");
	add_patch_output("file_dave_par_grass"                      , ANNUAL, "par_grass", [](const Patch &p) { return p.par_grass_mean / 1000000; }, "Annual mean PAR at top of grass canopy (MJ/m2/day)");

	// Annual per-patch individual-level outputs.
	add_aggregated_indiv_output("file_dave_alai"                , ANNUAL, get_lai, "Annual LAI(m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_afpc"                , ANNUAL, [](const Individual &i) { return i.fpc; }, "Annual FPC ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_acmass"              , ANNUAL, [](const Individual &i) { return i.ccont(); }, "Annual cmass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_anmass"              , ANNUAL, [](const Individual &i) { return i.ncont(); }, "Annual nmass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aheight"             , ANNUAL, [](const Individual &i) { return i.height; }, "Annual plant height (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_altor"               , ANNUAL, [](const Individual &i) { return i.ltor; }, "Annual plant leaf to root allocation ratio", mean, false, true);
	add_aggregated_indiv_output("file_dave_aaet"                , ANNUAL, [](const Individual &i) { return i.aaet; }, "Annual actual evapotranspiration (mm)", sum, true);
	add_aggregated_indiv_output("file_dave_adensity"            , ANNUAL, [](const Individual &i) { return i.densindiv; }, "Annual density of individuals over patch (/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_altor"               , ANNUAL, [](const Individual &i) { return i.ltor; }, "Annual leaf to root ratio (unitless)", mean, true);
	add_aggregated_indiv_output("file_dave_anuptake"            , ANNUAL, [](const Individual &i) { return i.anuptake; }, "Total annual nitrogen uptake (kgN/m2/year)", sum, true);
	add_aggregated_indiv_output("file_dave_a_aboveground_cmass" , ANNUAL, [](const Individual &i) { return i.get_cmass_leaf() + i.cmass_wood() * i.pft.ag_woody_frac + i.cmass_repr; }, "Annual total above-ground C biomass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_a_belowground_cmass" , ANNUAL, [](const Individual &i) { return i.cmass_root + i.cmass_crown + i.cmass_wood() * (1.0 - i.pft.ag_woody_frac); }, "Annual total below-ground C biomass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_a_aboveground_nmass" , ANNUAL, [](const Individual &i) { return i.get_nmass_leaf() + i.nmass_wood() * i.pft.ag_woody_frac + i.nmass_repr; }, "Annual total above-ground N biomass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_a_belowground_nmass" , ANNUAL, [](const Individual &i) { return i.nmass_root + i.nmass_crown + i.nmass_wood() * (1.0 - i.pft.ag_woody_frac); }, "Annual total below-ground N biomass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aresponse_window"    , ANNUAL, [](const Individual &i) { return (double)i.get_wscal_dynamic_capacity(); }, "Annual output of response window length (days)", mean);
	add_aggregated_indiv_output("file_dave_anretranslocated"    , ANNUAL, [](const Individual &i) { return i.anretranslocated; }, "Annual nitrogen retranslocated (kgN/m2)", sum);
	add_aggregated_indiv_output("file_dave_amort_greff"         , ANNUAL, [](const Individual &i) { return i.mort_greff; }, "Annual PFT-level growth-efficiency-related mortality, expressed as fraction of individuals within the cohort which died this year, summed across all indivs of a given PFT.", sum);
	add_aggregated_indiv_output("file_dave_amort_age"           , ANNUAL, [](const Individual &i) { return i.mort_age; }, "Annual PFT-level age-related mortality, expressed as fraction of individuals within the cohort which died this year, summed across all indivs of a given PFT.", sum);
	add_aggregated_indiv_output("file_dave_acmass_germ"         , ANNUAL, [](const Individual &i) { return i.cmass_germ; }, "Annual germinated C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aroot_beta"          , ANNUAL, [](const Individual &i) { return i.root_beta; }, "Annual output of pft-level root beta value (0-1)", mean);

	// Annual per-patch pft-level outputs.
	add_pft_output("file_dave_anpp"                             , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id); }, "Total individual-level Annual NPP (kgC/m2/year)", true);
	add_pft_output("file_dave_agpp"                             , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id); }, "Total individual-level Annual GPP (kgC/m2/year)", true);
	add_pft_output("file_dave_acue"                             , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id) ? patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id) / patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id) : 0; }, "Annual CUE (unitless)", true);
	add_pft_output("file_dave_aresp"                            , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::RA, pft.id); }, "Total individual-level Annual autotrophic Respiration (kgC/m2/year)", true);
	add_pft_output("file_dave_acmass_mort"                      , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].acmass_mort_kill; }, "Annual live vegetation C removed at Individual::kill() events (kgC/m2/year)", true);
	add_pft_output("file_dave_ancohort"                         , ANNUAL, get_nindiv, "Get the number of cohorts of this PFT currently established in this patch");
	add_pft_output("file_dave_anetps_ff"                        , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].anetps_ff; }, "Annual net photosynthesis at forest floor (kgC/m2)", true);
	add_pft_output("file_dave_ansapling"                        , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].nsapling; }, "Number of saplings of each PFT ready to establish (trees only)", true);
	add_pft_output("file_dave_ansapling_year"                   , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].nsapling_year; }, "Recruitment rate of each PFT (saplings per year) (trees only)", true);
	add_pft_output("file_dave_cton_leaf_max"                    , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_leaf_max; }, "Maximum C:N ratio of leaves (unitless)");
	add_pft_output("file_dave_cton_leaf_min"                    , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_leaf_min; }, "Minimum C:N ratio of leaves (unitless)");
	add_pft_output("file_dave_cton_root_max"                    , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_root_max; }, "Maximum C:N ratio of roots (unitless)");
	add_pft_output("file_dave_cton_root_min"                    , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_root_min; }, "Minimum C:N ratio of roots (unitless)");
	add_pft_output("file_dave_cton_sap_max"                     , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_sap_max; }, "Maximum C:N ratio of sapwood (unitless)");
	add_pft_output("file_dave_cton_sap_min"                     , ANNUAL, [](const Patch &p, const Pft &pft) { return pft.cton_sap_min; }, "Minimum C:N ratio of sapwood (unitless)");
	add_pft_output("file_dave_amort_cstarve"                    , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].amort_cstarve; }, "Annual PFT-level carbon starvation-related mortality, expressed as fraction of individuals within the cohort which died this year, summed across all indivs of a given PFT.");
	add_pft_output("file_dave_amort_cstarve_ntree"              , ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].amort_cstarve_ntree; }, "Annual PFT-level carbon starvation-related mortality, expressed as fraction of individuals within the cohort which died this year, summed across all indivs of a given PFT.");

	// Annual stand-level outputs.
	add_stand_output("file_dave_stand_frac"                     , ANNUAL, "fraction", [](const Stand &s) { return s.get_gridcell_fraction(); }, "Fraction of the gridcell occupied by each stand");
	add_stand_output("file_dave_stand_type"                     , ANNUAL, "type", [](const Stand &s) { return (double)s.landcover; }, "Stand landcover types");

	// Annual stand-level per-pft outputs.
	add_spft_output("file_dave_anetps_ff_max"                   , ANNUAL, [](const Standpft &spft) { return spft.anetps_ff_max; }, "Maximum recorded annual net photosynthesis at forest-floor conditions (kgC/m2)",  false, true);

	// Annual per-patch aggregated individual-level outputs.
	add_aggregated_indiv_output("file_dave_acalloc_leaf"	    , ANNUAL, [](const Individual &i) { return i.acalloc_leaf; }, "Total annual C allocation to leaf (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acalloc_root"	    , ANNUAL, [](const Individual &i) { return i.acalloc_root; }, "Total annual C allocation to root (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acalloc_repr"	    , ANNUAL, [](const Individual &i) { return i.acalloc_repr; }, "Total annual C allocation to repr (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acalloc_sap"		    , ANNUAL, [](const Individual &i) { return i.acalloc_sap; }, "Total annual C allocation to sap (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acalloc_storage"	    , ANNUAL, [](const Individual &i) { return i.acalloc_storage; }, "Total annual C allocation to storage (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acalloc_crown"	    , ANNUAL, [](const Individual &i) { return i.acalloc_crown; }, "Total annual C allocation to crown (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_anstorage_drawdown"  , ANNUAL, [](const Individual &i) { return i.anstorage_drawdown; }, "Total annual amount of N withdrawn from storage for growth (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_acloss_temp"         , ANNUAL, [](const Individual &i) { return i.acloss_temp; }, "Total annual amount of C lost to high temperature (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_anloss_temp"         , ANNUAL, [](const Individual &i) { return i.anloss_temp; }, "Total annual amount of N lost to high temperature (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aloss_temp"          , ANNUAL, [](const Individual &i) { return i.aloss_temp; }, "Total annual fraction of leaves lost to high temperature (unitless)", mean, true);

	// Annual indiv-level outputs.
	add_indiv_output("file_dave_indiv_mort"                     , ANNUAL, mortality_cols, "Annual indiv-level mortality divided by component/cause, expressed as fraction of individuals within the cohort which died this year.", true);
	add_indiv_output("file_dave_indiv_age"                      , ANNUAL, "age", [](const Individual &i) { return i.age; }, "Individual-level age (years)");

	// Daily gridcell-level outputs.
	add_gridcell_aggregated_indiv_output("file_agb"             , DAILY, [](const Individual &i) { return i.aboveground_biomass(); }, sum, "Daily gridcell-level above-ground biomass (kg/m2)", true);
	add_gridcell_aggregated_tree_output("file_agb_tree"         , DAILY, [](const Individual &i) { return i.aboveground_biomass(); }, sum, "Daily gridcell-level above-ground tree biomass (kg/m2)", true);

	// Daily patch-level outputs.
	add_patch_output("file_dave_daylength"                      , DAILY, "daylength", [](const Patch &p) { return p.get_climate().daylength; }, "Daily patch-level day-length (h)");
	add_patch_output("file_dave_cpool"                          , DAILY, cpool_columns, "Daily C pools (kgC/m2)", true);
	add_patch_output("file_dave_soil_nmass_avail"               , DAILY, "nmass", [](const Patch &p) { return sum_layers(p.soil, [](const Soillayer &layer) { return layer.NO3_mass + layer.NH4_mass; }); }, "Available Soil N for plant uptake");
	add_patch_output("file_dave_soil_nh4_mass"                  , DAILY, "nh4_mass", [](const Patch &p) { return p.soil.NH4_mass; }, "NH4 mass in soil (kgN/m2)");
	add_patch_output("file_dave_soil_no3_mass"                  , DAILY, "no3_mass", [](const Patch &p) { return p.soil.NO3_mass; }, "NO3 mass in soil (kgN/m2)");
	add_patch_output("file_dave_dsimfire"                       , DAILY, simfire_columns, "Daily simfire analysis");
	add_patch_output("file_dave_dfuel"                          , DAILY, fuel_columns, "Daily blaze fuel availability (kgC/m2)");
	add_patch_output("file_dave_dcoarse_woody_debris"           , DAILY, coarse_woody_debris_columns, "Daily coarse woody debris (gC/m2)");
	add_patch_output("file_dave_met_temp"                       , DAILY, "temp", [](const Patch &patch) { return patch.stand.get_gridcell().climate.temp;    }, "Daily air temperature (°C)");
	add_patch_output("file_dave_met_par"                        , DAILY, "par", [](const Patch &patch) { return patch.stand.get_gridcell().climate.par * 0.001;     }, "Daily PAR (kJ/m2/timestep)");
	add_patch_output("file_dave_met_par_max"                    , DAILY, "par_max", [](const Patch &patch) { return patch.stand.get_gridcell().climate.par_max/1000000;     }, "Maximum daily PAR (MJ/m2/timestep)");
	add_patch_output("file_dave_rel_par"                        , DAILY, "rel_par", [](const Patch &patch) { return patch.stand.get_gridcell().climate.rel_par;     }, "Relative PAR (unitless)");
	add_patch_output("file_dave_met_vpd"                        , DAILY, "vpd", [](const Patch &patch) { return mean(&patch.stand.get_gridcell().climate.vpds.front(), patch.stand.get_gridcell().climate.vpds.size());     }, "Daily VPD (kPa)");
	add_patch_output("file_dave_met_insol"                      , DAILY, "insol", [](const Patch &patch) { return patch.stand.get_gridcell().climate.insol;   }, "Daily insolation (units depend on instype)");
	add_patch_output("file_dave_met_precip"                     , DAILY, "precip", [](const Patch &patch) { return patch.stand.get_gridcell().climate.prec;    }, "Daily total precipitation (mm)");
	add_patch_output("file_dave_met_pressure"                   , DAILY, "pressure", [](const Patch &patch) { return mean(&patch.stand.get_gridcell().climate.pressure.front(), patch.stand.get_gridcell().climate.pressure.size()) * 0.001; }, "Daily air pressure (kPa)");
	add_patch_output("file_dave_met_co2"                        , DAILY, "co2", [](const Patch &patch) { return patch.stand.get_gridcell().climate.co2;         }, "Daily atmospheric co2 concentration (ppm)");
	add_patch_output("file_dave_daylength"                      , DAILY, "daylength", [](const Patch &patch) { return patch.stand.get_gridcell().climate.daylength;         }, "Daylength (h)");
	add_patch_output("file_dave_sompool_cmass"                  , DAILY, sompool_cmass_columns, "Daily SOM pool C mass (kgC/m2)", true);
	add_patch_output("file_dave_sompool_nmass"                  , DAILY, sompool_nmass_columns, "Daily SOM pool N mass (kgN/m2)", true);
	add_patch_output("file_dave_sompool_nimmob"                 , DAILY, sompool_nimmob_columns, "Daily SOM pool N immobilisation (kgN/m2)", true);
	add_patch_output("file_dave_sompool_nmin_gross"             , DAILY, sompool_nmin_gross_columns, "Daily SOM pool N immobilisation (kgN/m2)", true);
	add_patch_output("file_dave_sompool_decay_reduction"        , DAILY, sompool_decay_reduction_columns, "Daily SOM pool decay reduction (0-1)", true);
	add_patch_output("file_dave_sompool_frac_remaining"         , DAILY, sompool_frac_remaining_columns, "Daily SOM pool fraction remaining (0-1)", true);
	add_patch_output("file_dave_sompool_ntoc"                   , DAILY, sompool_ntoc_columns, "Daily SOM pool N to C ratio (kgN/kgC)");
	add_patch_output("file_dave_ninput"                         , DAILY, ninput_columns, "Daily nitrogen deposition (kgN/m2)");
	add_patch_output("file_dave_fpar_ff"                        , DAILY, "fpar_ff", [](const Patch &p) { return p.fpar_ff; }, "Forest floor fpar (0-1)");
	add_patch_output("file_dave_aet"                            , DAILY, daily_aet_columns, "Daily patch-level evaporation (mm)", true);
	add_patch_output("file_dave_fevap"                          , DAILY, fevap_columns, "Daily patch-level fraction of patch subject to soil evap", true);
	add_patch_output("file_dave_resp_heterotrophic"             , DAILY, "resp_heterotrophic", [](const Patch &p) { return p.fluxes.get_daily_flux(Fluxes::SOILC, date.day) * G_PER_KG; }, "Daily carbon flux to atmosphere from soil respiration (gC/m2/day)");
	add_patch_output("file_dave_resp"                           , DAILY, "resp", get_er, "Daily patch-level ecosystem respiration (gC/m2/day)");
	add_patch_output("file_dave_gpp"                            , DAILY, "gpp", get_gpp, "Daily patch-level GPP (gC/m2/day)");
	add_patch_output("file_dave_npp"                            , DAILY, "npp", get_npp, "Daily patch-level NPP (gC/m2/day)");
	add_patch_output("file_dave_nee"                            , DAILY, "nee", get_nee, "Daily patch-level NEE (gC/m2/day)");
	add_patch_output("file_dave_evaporation"                    , DAILY, daily_evap_columns, "Daily patch-level evaporation (mm/day)");
	add_patch_output("file_dave_soilc"                          , DAILY, soil_cmass_columns, "Daily soil carbon content (kgC/m2)", true);
	add_patch_output("file_dave_soiln"                          , DAILY, soil_nmass_columns, "Daily soil nitrogen content (kgN/m2)", true);
	add_patch_output("file_dave_soil_psi"                       , DAILY,"soil_psi", [](const Patch &p) { return mean_layers(p.soil, [](const Soillayer &l) { return l.psi_mpa(); }); }, "Soil psi (MPa)");
	add_patch_output("file_dave_pet"                            , DAILY, "pet", [](const Patch &patch) { return patch.stand.get_gridcell().climate.eet; }, "Daily potential patch-level evapotranspiration (mm)");
	add_patch_output("file_dave_swmm_10"                        , DAILY, "total", [](const Patch &p) { return get_soil_moisture(p.soil, 0, 100); }, "Soil Moisture at depth: 0-0.1m (mm)");
	add_patch_output("file_dave_swmm_100"                       , DAILY, "total", [](const Patch &p) { return get_soil_moisture(p.soil, 100, 1000); }, "Soil Moisture at depth: 0.1-1m (mm)");
	add_patch_output("file_dave_swmm_600"                       , DAILY, "total", [](const Patch &p) { return get_6m_soil_moisture_estimate(p.soil); }, "Soil Moisture at depth: 0.1-1m (mm)");
	add_patch_output("file_dave_swavail_100"                    , DAILY, "total", [](const Patch &p) { return get_available_soil_moisture(p.soil, 100, 1000); }, "Available Soil Moisture at depth: 0.1-1m (mm)");
	add_patch_output("file_dave_water_budget_daily"             , DAILY, water_budget_cols_daily, "Daily water budget (mm)");
	add_patch_output("file_dave_nitrogen_budget_daily"          , DAILY, nitrogen_budget_cols_daily, "Daily nitrogen budget (kgN/m2)");
	add_patch_output("file_dave_dnpool"                         , DAILY, npool_columns, "Daily N pools (kgN/m2)");
	add_patch_output("file_dave_wcont_upper"                    , DAILY, "wcont", [](const Patch &p) { return p.soil.var_depth(0, SOILDEPTH_UPPER, &Soillayer::wcont); }, "Daily soil moisture in the upper 50cm (0-1)");
	add_patch_output("file_dave_soiltemp25"                     , DAILY, "soiltemp25", [](const Patch &p) { return p.soil.temp_at_depth(25 * MM_PER_CM); }, "Daily soil temperature at 25cm (°C)");
	add_patch_output("file_dave_clitter"                        , DAILY, clitter_columns, "Daily C Mass in Litter (kgC/m2)", true);
	add_patch_output("file_dave_nlitter"                        , DAILY, nlitter_columns, "Daily N Mass in Litter (kgN/m2)", true);
	add_patch_output("file_dave_patch_ba"                       , DAILY, "total", get_patch_ba, "Total daily patch-level basal area (m2/ha)");

	// Daily per-patch aggregated individual-level outputs.
	add_aggregated_indiv_output("file_dave_lai"                 , DAILY, [](const Individual &i) { return i.lai_today(); }, "LAI by individual (m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_meanleafage"         , DAILY, [](const Individual &i) { return i.meanleafage; }, "Mean leaf age (years)", mean, false, true);
	add_aggregated_indiv_output("file_dave_meansapage"          , DAILY, [](const Individual &i) { return i.meansapage; }, "Mean sapwood age (years)", mean, false, true);
	add_aggregated_indiv_output("file_dave_meanrootage"         , DAILY, [](const Individual &i) { return i.meanrootage; }, "Mean root age (years)", mean, false, true);
	add_aggregated_indiv_output("file_dave_sla"                 , DAILY, [](const Individual &i) { return i.sla(); }, "SLA by individual (m2/kgC)", mean, false, true);
	add_aggregated_indiv_output("file_dave_cmass_leaf"          , DAILY, [](const Individual &i) { return i.get_cmass_leaf(); }, "Leaf C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_leaf_brown"    , DAILY, [](const Individual &i) { return i.pft.is_grass() ? i.cmass_leafs[1] : 0.0; }, "Leaf brown C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_leaf"          , DAILY, [](const Individual &i) { return i.get_nmass_leaf(); }, "Leaf N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_crown"         , DAILY, [](const Individual &i) { return i.cmass_crown; }, "Crown C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_crown"         , DAILY, [](const Individual &i) { return i.nmass_crown; }, "Crown N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_repr"          , DAILY, [](const Individual &i) { return i.cmass_repr; }, "Reproductive C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_root"          , DAILY, [](const Individual &i) { return i.cmass_root; }, "Root C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_sap"           , DAILY, [](const Individual &i) { return i.cmass_sap; }, "Sapwood C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cton_sap"            , DAILY, [](const Individual &i) { return min(i.pft.cton_sap_max, i.cton_sap()); }, "Sapwood C to N ratio (kgC/kgN)", sum, true);
	add_aggregated_indiv_output("file_dave_cton_leaf"           , DAILY, [](const Individual &i) { return min(i.pft.cton_leaf_max, i.cton_leaf()); }, "Leaf C to N ratio (kgC/kgN)", mean);
	add_aggregated_indiv_output("file_dave_cton_root"           , DAILY, [](const Individual &i) { return min(i.pft.cton_root_max, i.cton_root()); }, "Root C to N ratio (kgC/kgN)", mean);
	add_aggregated_indiv_output("file_dave_cton_heart"          , DAILY, [](const Individual &i) { return min(i.pft.cton_sap_max, i.cton_heart()); }, "Heartwood C to N ratio (kgC/kgN)", mean);
	add_aggregated_indiv_output("file_dave_cton_storage"        , DAILY, [](const Individual &i) { return i.cmass_storage / i.nmass_storage; }, "Storage C to N ratio (kgC/kgN)", mean);
	add_aggregated_indiv_output("file_dave_cton_crown"          , DAILY, [](const Individual &i) { return min(i.pft.cton_root_max, i.cton_crown()); }, "Crown C to N ratio (kgC/kgN)", mean);
	add_aggregated_indiv_output("file_dave_cmass_heart"         , DAILY, [](const Individual &i) { return i.cmass_heart; }, "Heartwood C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_root"          , DAILY, [](const Individual &i) { return i.nmass_root; }, "Root N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass"               , DAILY, [](const Individual &i) { return i.ccont(); }, "Total C biomass per individual per patch (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass"               , DAILY, [](const Individual &i) { return i.ncont(); }, "Total N biomass per individual per patch (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage"       , DAILY, [](const Individual &i) { return i.cmass_storage; }, "Storage C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage_max"   , DAILY, [](const Individual &i) { return i.cstorage_max; }, "Storage capacity for C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_storage"       , DAILY, [](const Individual &i) { return i.nmass_storage; }, "Storage N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_storage_max"   , DAILY, [](const Individual &i) { return i.nstorage_max; }, "Storage capacity for N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage_dynam" , DAILY, [](const Individual &i) { return i.cmass_storage_dynam; }, "C mass storage dynamics (> 0 when saving C) (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_storage_dynam_14sum" , DAILY, [](const Individual &i) { return i.cstorage_dynam_14.sum(); }, "Mean 14-day cumulative C mass storage dynamics (> 0 if saving C) (kgC/m2)", mean, true);
	add_aggregated_indiv_output("file_dave_nmass_sap"           , DAILY, [](const Individual &i) { return i.nmass_sap; }, "Sapwood N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_heart"         , DAILY, [](const Individual &i) { return i.cmass_heart; }, "Heartwood C mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_heart"         , DAILY, [](const Individual &i) { return i.nmass_heart; }, "Heartwood N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_nmass_repr"          , DAILY, [](const Individual &i) { return i.nmass_repr; }, "Reproductive N mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_ndemand"             , DAILY, [](const Individual &i) { return i.ndemand; }, "N demand mass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_height"              , DAILY, [](const Individual &i) { return i.height; }, "Plant height (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_diameter"            , DAILY, [](const Individual &i) { return i.diam; }, "Stem diameter (m)", mean, false, true);
	add_aggregated_indiv_output("file_dave_basalarea"           , DAILY, [](const Individual &i) { return i.basalarea; }, "Mean Basal Area (m2/tree)", mean, false, true);
	add_aggregated_indiv_output("file_dave_basalarea_inc"       , DAILY, [](const Individual &i) { return i.basalarea_inc * 10000; }, "Mean Basal Area increment (cm2/tree/day)", mean, false, true);
	add_aggregated_indiv_output("file_dave_density"             , DAILY, [](const Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_sapwood_area"        , DAILY, [](const Individual &i) { return i.sa; }, "sapwood area sum of individuals (m2/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_sapwood_depth"       , DAILY, [](const Individual &i) { return i.sapwood_depth; }, "Sapwood depth sum over individuals (m)", sum, true);
	add_aggregated_indiv_output("file_dave_latosa"              , DAILY, [](const Individual &i) { return i.latosa /1000; }, "leaf area to sapwood area ratio / 1000 (unitless)", mean, true);
	add_aggregated_indiv_output("file_dave_fpar"                , DAILY, [](const Individual &i) { return i.fpar; }, "Daily fpar [0-1]", sum, true);
	add_aggregated_indiv_output("file_dave_indiv_gpp"           , DAILY, [](const Individual &i) { return i.dgpp * G_PER_KG; }, "Daily GPP (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_autotrophic"    , DAILY, [](const Individual &i) { return (i.dgpp - i.dnpp) * G_PER_KG; }, "Daily autotrophic respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_maintenance"    , DAILY, [](const Individual &i) { return i.resp_maintenance * G_PER_KG; }, "Daily maintenance respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_growth"         , DAILY, [](const Individual &i) { return i.resp_growth * G_PER_KG; }, "Daily growth respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_repr"           , DAILY, [](const Individual &i) { return i.resp_repr * G_PER_KG; }, "Daily reproductive respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_crown"          , DAILY, [](const Individual &i) { return i.resp_crown * G_PER_KG; }, "Daily crown respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_root"           , DAILY, [](const Individual &i) { return i.resp_root * G_PER_KG; }, "Daily root respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_sap"            , DAILY, [](const Individual &i) { return i.resp_sap * G_PER_KG; }, "Daily sapwood respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_indiv_npp"           , DAILY, [](const Individual &i) { return i.dnpp * G_PER_KG; }, "Daily NPP (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_fpc"                 , DAILY, [](const Individual &i) { return i.fpc; }, "Daily FPC ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_crownarea"           , DAILY, [](const Individual &i) { return i.crownarea; }, "Daily Crown Area (m2)", sum, true);
	add_aggregated_indiv_output("file_dave_agd_g"               , DAILY, [](const Individual &i) { return i.photosynthesis.agd_g; }, "Daily gross photosynthesis (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_rd_g"                , DAILY, [](const Individual &i) { return i.photosynthesis.rd_g; }, "Daily leaf respiration (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_je"                  , DAILY, [](const Individual &i) { return i.photosynthesis.je; }, "Daily PAR-limited photosynthesis rate (gC/m2/h)", sum, true);
	add_aggregated_indiv_output("file_dave_vcmax"               , DAILY, [](const Individual &i) { return i.photosynthesis.vm; }, "RuBisCO capacity (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_vcmax25"             , DAILY, [](const Individual &i) { return i.ps_subdaily.vcmax25; }, "Max carboxylation rate at 25 degrees (vcmax 25) (umol/m2/s)", mean);
	add_aggregated_indiv_output("file_dave_jmax25"              , DAILY, [](const Individual &i) { return i.ps_subdaily.jmax25; }, "Maximum e- transport rate at 25 degrees (jmax) (umol/m2/s)", mean);
	add_aggregated_indiv_output("file_dave_fphen_activity"      , DAILY, [](const Individual &i) { return i.fphen_activity; }, "Dormancy downregulation (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_fdev_growth"         , DAILY, [](const Individual &i) { return i.fdev_growth; }, "Development factor for growth demand (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_frepr_cstorage"      , DAILY, [](const Individual &i) { return i.frepr_cstorage; }, "Ratio of reproductive biomass to storage capacity", mean, false, true);
	add_aggregated_indiv_output("file_dave_growth_demand"       , DAILY, [](const Individual &i) { return i.growth_demand; }, "Growth Demand (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_cmass_germ"          , DAILY, [](const Individual &i) { return i.cmass_germ; }, "Germinated C Mass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_transpiration"       , DAILY, [](const Individual &i) { return i.aet; }, "Daily individual-level transpiration (mm/day)", sum, true, true);
	add_aggregated_indiv_output("file_dave_nscal"               , DAILY, [](const Individual &i) { return i.nscal_5.lastadd(); }, "N stress scalar (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_nscal_mean"          , DAILY, [](const Individual &i) { return i.nscal_5.mean(); }, "N stress scalar 5 day running mean (0-1)", mean, false, true);
	add_aggregated_indiv_output("file_dave_ltor"            	, DAILY, [](const Individual &i) { return i.ltor; }, "Daily leaf to root ratio (unitless)", mean, true);
	add_aggregated_indiv_output("file_dave_cue"                 , DAILY, [](const Individual &i) { return i.dnpp / i.dgpp; }, "carbon use efficiency ()", mean, false, true);
	add_aggregated_indiv_output("file_dave_weighted_psi"        , DAILY, [](const Individual &i) { return i.weighted_psi; }, "Weighted psi (MPa)", mean, false, true);
	add_aggregated_indiv_output("file_dave_phen"                , DAILY, [](const Individual &i) { return i.phen; }, "Phenology multiplier", mean, false, false);
	add_aggregated_indiv_output("file_dave_lai_today"           , DAILY, [](const Individual &i) { return i.lai_today(); }, "LAI today (m2/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_kmax_drought_affected", DAILY, [](const Individual &i) { return i.kmax_drought_affected; }, "kmax accounting for drought impacts (mol/m2/s/MPa)", mean, false, false);
	add_aggregated_indiv_output("file_dave_max_k_loss"          , DAILY, [](const Individual &i) { return i.max_k_loss; }, "reduction in kmax due to xylem embolism", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_leaf"          , DAILY, [](const Individual &i) { return i.alpha_leaf; }, "Daily leaf sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_root"          , DAILY, [](const Individual &i) { return i.alpha_root; }, "Daily root sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_sap"           , DAILY, [](const Individual &i) { return i.alpha_sap; }, "Daily sap sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_alpha_storage"       , DAILY, [](const Individual &i) { return i.alpha_storage; }, "Daily storage sink strength (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cdeficit_leaf"       , DAILY, [](const Individual &i) { return i.cdeficit_leaf; }, "Daily leaf C deficit (kgC/m2)", sum, false, false);
	add_aggregated_indiv_output("file_dave_cdeficit_root"       , DAILY, [](const Individual &i) { return i.cdeficit_root; }, "Daily root C deficit (kgC/m2)", sum, false, false);
	add_aggregated_indiv_output("file_dave_cdeficit_sap"        , DAILY, [](const Individual &i) { return i.cdeficit_sap; }, "Daily sap C deficit (kgC/m2)", sum, false, false);
	add_aggregated_indiv_output("file_dave_cdeficit_storage"    , DAILY, [](const Individual &i) { return i.cdeficit_storage; }, "Daily storage C deficit (kgC/m2)", sum, false, false);
	add_aggregated_indiv_output("file_dave_cmass_leaf_limit"    , DAILY, [](const Individual &i) { return i.cmass_leaf_limit; }, "Daily optimum leaf pool size for trees (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_root_limit"    , DAILY, [](const Individual &i) { return i.cmass_root_limit; }, "Daily optimum root pool size for trees (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_sap_limit"     , DAILY, [](const Individual &i) { return i.cmass_sap_limit; }, "Daily optimum sap pool size for trees (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_cmass_storage_limit" , DAILY, [](const Individual &i) { return i.cmass_storage_limit; }, "Daily optimum storage pool size (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_absolute_alpha_leaf"   , DAILY, [](const Individual &i) { return i.absolute_alpha_leaf; }, "Relative leaf sink strength (unitless)", mean, false, false);
	add_aggregated_indiv_output("file_dave_absolute_alpha_root"   , DAILY, [](const Individual &i) { return i.absolute_alpha_root; }, "Relative root sink strength (unitless)", mean, false, false);
	add_aggregated_indiv_output("file_dave_absolute_alpha_sap"    , DAILY, [](const Individual &i) { return i.absolute_alpha_sap; }, "Relative sapwood sink strength (unitless)", mean, false, false);
	add_aggregated_indiv_output("file_dave_absolute_alpha_storage", DAILY, [](const Individual &i) { return i.absolute_alpha_storage; }, "Relative storage sink strength (unitless)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cgrow_leaf"          , DAILY, [](const Individual &i) { return i.cgrow_leaf; }, "Daily leaf growth (kgC/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_root"          , DAILY, [](const Individual &i) { return i.cgrow_root; }, "Daily root growth (kgC/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_sap"           , DAILY, [](const Individual &i) { return i.cgrow_sap; }, "Daily sap growth (kgC/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_repr"          , DAILY, [](const Individual &i) { return i.cgrow_repr; }, "Daily repr growth (kgC/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_cgrow_storage"       , DAILY, [](const Individual &i) { return i.cgrow_storage; }, "Daily storage allocation (kgC/m2)", sum, true, false);
	add_aggregated_indiv_output("file_dave_creturn"             , DAILY, [](const Individual &i) { return i.creturn; }, "daily storage growth plus excess from N limitation (kgC/m2/day)", sum, true, false);
	add_aggregated_indiv_output("file_dave_cgrow"               , DAILY, [](const Individual &i) { return i.cgrow; }, "Daily growth (kgC/m2)", mean, true, false);
	add_aggregated_indiv_output("file_dave_diameter_inc"        , DAILY, [](const Individual &i) { return i.d_diam; }, "Daily diameter growth (m/day)", mean, true, false);
	add_aggregated_indiv_output("file_dave_height_inc"          , DAILY, [](const Individual &i) { return i.d_height; }, "Daily height growth (m/day)", mean, true, false);
	add_aggregated_indiv_output("file_dave_cturnover_leaf"      , DAILY, [](const Individual &i) { return i.cturnover_leaf; }, "Daily leaf C turnover flux (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cturnover_root"      , DAILY, [](const Individual &i) { return i.cturnover_root; }, "Daily root C turnover flux (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cturnover_sap"       , DAILY, [](const Individual &i) { return i.cturnover_sap; }, "Daily sapwood C turnover flux (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cturnover_repr"      , DAILY, [](const Individual &i) { return i.cturnover_repr; }, "Daily reproductive C turnover flux (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cturnover_crown"     , DAILY, [](const Individual &i) { return i.cturnover_crown; }, "Daily crown C turnover flux (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nturnover_leaf"      , DAILY, [](const Individual &i) { return i.nturnover_leaf * G_PER_KG; }, "Daily leaf N turnover flux (gN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nturnover_root"      , DAILY, [](const Individual &i) { return i.nturnover_root * G_PER_KG; }, "Daily root N turnover flux (gN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nturnover_sap"       , DAILY, [](const Individual &i) { return i.nturnover_sap * G_PER_KG; }, "Daily sapwood N turnover flux (gN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nturnover_repr"      , DAILY, [](const Individual &i) { return i.nturnover_repr * G_PER_KG; }, "Daily reproductive N turnover flux (gN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nturnover_crown"     , DAILY, [](const Individual &i) { return i.nturnover_crown * G_PER_KG; }, "Daily crown N turnover flux (gN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_turnover_leaf"       , DAILY, [](const Individual &i) { return i.pturnover_leaf; }, "Daily leaf turnover as a proportion of leaf C biomass (unitless)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anc_frac"            , DAILY, [](const Individual &i) { return i.anc_frac; }, "Daily anc as fraction of total daily photosynthesis (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anj_frac"            , DAILY, [](const Individual &i) { return i.anj_frac; }, "Daily anj as fraction of total daily photosynthesis (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_anp_frac"            , DAILY, [](const Individual &i) { return i.anp_frac; }, "Daily anp as fraction of total daily photosynthesis (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_dnuptake"            , DAILY, [](const Individual &i) { return i.dnuptake; }, "Total daily nitrogen uptake (kgN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cexcess"             , DAILY, [](const Individual &i) { return i.cexcess; }, "Total daily carbon overflow (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_aboveground_cmass"   , DAILY, [](const Individual &i) { return i.get_cmass_leaf() + i.cmass_wood() * i.pft.ag_woody_frac + i.cmass_repr; }, "Daily above-ground C biomass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_belowground_cmass"   , DAILY, [](const Individual &i) { return i.cmass_root + i.cmass_crown + i.cmass_wood() * (1.0 - i.pft.ag_woody_frac); }, "Daily below-ground C biomass (kgC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aboveground_nmass"   , DAILY, [](const Individual &i) { return i.get_nmass_leaf() + i.nmass_wood() * i.pft.ag_woody_frac + i.nmass_repr; }, "Daily above-ground N biomass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_belowground_nmass"   , DAILY, [](const Individual &i) { return i.nmass_root + i.nmass_crown + i.nmass_wood() * (1.0 - i.pft.ag_woody_frac); }, "Daily below-ground N biomass (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_aboveground_tree_biomass", DAILY, [](const Individual &i) { return i.pft.is_grass() ? 0 : 2.0 * (i.cmass_leaf + i.cmass_wood() * i.pft.ag_woody_frac /* + i.cmass_storage*/); }, "Daily total above-ground tree biomass (kg/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_csenesc"             , DAILY, [](const Individual &i) { return G_PER_KG * i.csenesc; }, "Daily total C senesced by grasses (gC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_closs_temp"          , DAILY, [](const Individual &i) { return i.closs_temp; }, "Daily C lost to high temperature (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_nloss_temp"          , DAILY, [](const Individual &i) { return i.nloss_temp; }, "Daily N lost to high temperature (kgN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_loss_temp"           , DAILY, [](const Individual &i) { return i.loss_temp; }, "Daily fraction of leaves lost to high temperature (unitless)", mean, true);

	add_aggregated_indiv_output("file_dave_cgrazed_repr"        , DAILY, [](const Individual &i) { return i.charvest_repr; }, "Daily reproductive cmass grazed (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cgrazed_green_leaf"  , DAILY, [](const Individual &i) { return i.charvest_green; }, "Daily green leaf cmass grazed (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_cgrazed_brown_leaf"  , DAILY, [](const Individual &i) { return i.charvest_brown; }, "Daily brown leaf cmass grazed (kgC/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_ngrazed_repr"        , DAILY, [](const Individual &i) { return i.nharvest_repr; }, "Daily reproductive nmass grazed (kgN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_ngrazed_green_leaf"  , DAILY, [](const Individual &i) { return i.nharvest_green; }, "Daily green leaf nmass grazed (kgN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_ngrazed_brown_leaf"  , DAILY, [](const Individual &i) { return i.nharvest_brown; }, "Daily brown leaf nmass grazed (kgN/m2/day)", sum, true);
	add_aggregated_indiv_output("file_dave_rpc"                 , DAILY, [](const Individual &i) { return i.rpc; }, "Daily root projective cover (0-1)", sum);
	add_aggregated_indiv_output("file_dave_nretranslocated"     , DAILY, [](const Individual &i) { return i.nretranslocated; }, "Annual nitrogen retranslocated (kgN/m2)", sum);

	// Daily outputs of annually-calculated pft-level outputs.
	// (Outputted daily for plotting purposes.)
	add_aggregated_indiv_output("file_dave_dresponse_window"    , DAILY, [](const Individual &i) { return (double)i.get_wscal_dynamic_capacity(); }, "Daily output of response window length (days)", mean);
	add_aggregated_indiv_output("file_dave_root_beta"           , DAILY, [](const Individual &i) { return i.root_beta; }, "Daily output of pft-level root beta value (0-1)", mean);

	add_aggregated_indiv_output("file_dave_live_biomass"        , DAILY, [](const Individual &i) { return i.ccont() + i.ncont(); }, "Daily total live biomass (kg/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_wscal"               , DAILY, [](const Individual &i) { return i.get_wscal(); }, "Water scalar [0-1], 1 = min stress",  mean, false, true);
	add_aggregated_indiv_output("file_dave_dphen"               , DAILY, [](const Individual &i) { return i.phen; }, "daily phen (-)", mean, false, false);
	add_aggregated_indiv_output("file_dave_ftemp"               , DAILY, [](const Individual &i) { return i.ftemp; }, "coefficient that scales growth demand using air temperature (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_fwater"              , DAILY, [](const Individual &i) { return i.fwater; }, "coefficient that scales growth demand using wscal (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_fstorage"            , DAILY, [](const Individual &i) { return i.fstorage; }, "coefficient that scales growth according to actualt storage state (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_fphoto"              , DAILY, [](const Individual &i) { return i.fphoto; }, "coefficient that scales growth according to relative day length (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_cgrow_pot"           , DAILY, [](const Individual &i) { return i.cgrow_pot; }, "coefficient used to scale cgrow by taking into account ftem, fwater and fstrage (0-1)", mean, false, false);
	add_aggregated_indiv_output("file_dave_resp_leaf"           , DAILY, [](const Individual &i) { return i.resp_leaf * G_PER_KG; }, "Total Leaf respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_root"           , DAILY, [](const Individual &i) { return i.resp_root * G_PER_KG; }, "Total Root respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_sap"            , DAILY, [](const Individual &i) { return i.resp_sap * G_PER_KG; }, "Total Sapwood respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_crown"          , DAILY, [](const Individual &i) { return i.resp_crown * G_PER_KG; }, "Total Crown respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_repr"           , DAILY, [](const Individual &i) { return i.resp_repr * G_PER_KG; }, "Total Reproductive respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_maintenance"    , DAILY, [](const Individual &i) { return i.resp_maintenance * G_PER_KG; }, "Total Maintenance respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_resp_growth"         , DAILY, [](const Individual &i) { return i.resp_growth * G_PER_KG; }, "Total Growth respiration by PFT (gC/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_fnuptake"            , DAILY, [](const Individual &i) { return i.fnuptake; }, "Fraction of nitrogen demand taken up", mean);
	add_aggregated_indiv_output("file_dave_nstorage_drawdown"   , DAILY, [](const Individual &i) { return i.nstorage_drawdown; }, "Amount of N withdrawn from storage today for growth (kgN/m2)", sum, true);
	add_aggregated_indiv_output("file_dave_g1"                  , DAILY, [](const Individual &i) { return i.g1; }, "Stomatal slope parameter (Medlyn et al. 2011) (kPa^0.5)", mean, false, true);

	// Daily individual-level outputs.
	// These will produce one row per individual per patch per timestep.
	add_indiv_output("file_dave_indiv_cpool"                    , DAILY, indiv_cpool_columns, "Indiv-level carbon pools (kgC/m2)", false, true);
	add_indiv_output("file_dave_indiv_npool"                    , DAILY, indiv_npool_columns, "Indiv-level nitrogen pools (kgN/m2)", false, true);
	add_indiv_output("file_dave_hydraulics"                     , DAILY, indiv_hydraulics_columns, "kmax, height, lai, transpiration, daily minimum leaf water potential and weighted soil water potential at individual level", false, false);
	add_indiv_output("file_dave_indiv_lai"                    	, DAILY, "lai", [](const Individual &i) { return i.lai; }, "Cohort-level LAI (m2/m2)");
	add_indiv_output("file_dave_indiv_c_d_dbh"                  , ANNUAL, "d_dbh", [](const Individual &i) { return i.acalloc_sap; }, "annual diamter increment (m)");
	add_indiv_output("file_dave_indiv_respiration"              , DAILY, indiv_resp_columns, "Indiv-level respiration components (kgC/m2)"); // No auto total column
	add_indiv_output("file_dave_indiv_tree_cpool"               , DAILY, indiv_cpool_tree_columns, "Indiv-level tree carbon pools (kgC/m2)", false, true);
	add_indiv_output("file_dave_indiv_dbh"             	        , DAILY, indiv_dbhcoloumns_columns, "Cohort-level DBH and Delta DBH (m, m/yr)");
	add_indiv_output("file_dave_annual_indiv_dbh"             	, ANNUAL, "dbh", [](const Individual &i) { return i.diam; }, "Cohort-level DBH  (m)");
	add_indiv_output("file_dave_annual_indiv_c_d_dbh"           , ANNUAL, "d_dbh", [](const Individual &i) { return i.c_d_diam; }, "Annual diameter increment (m/yr)");
	add_indiv_output("file_dave_annual_lai"          		    , ANNUAL, "alai", [](const Individual &i) { return i.lai; }, "Annual LAI");
	add_indiv_output("file_dave_indiv_age"             	        , DAILY, "age", [](const Individual &i) { return i.age; }, "Cohort-level Age (years)");
	add_indiv_output("file_dave_indiv_nmass_root"             	, DAILY, "nmass_root", [](const Individual &i) { return i.nmass_root; }, "Cohort-level root N content (kgN/m2)");
	add_indiv_output("file_dave_indiv_height"                   , DAILY, "height", [](const Individual &i) { return i.height; }, "Cohort-level Height (m)");
	add_indiv_output("file_dave_indiv_crownarea"                , DAILY, "crownarea", [](const Individual &i) { return i.crownarea; }, "crownarea (m2)");
	add_indiv_output("file_dave_indiv_fpc"                      , DAILY, "fpc", [](const Individual &i) { return i.fpc; }, "Cohort-level FPC (-)");
	add_indiv_output("file_dave_indiv_fpar"                     , DAILY, "fpar", [](const Individual &i) { return i.fpar; }, "Cohort-level fpar (-)");
	add_indiv_output("file_dave_indiv_density"                  , DAILY, "density", [](const Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)");
	add_indiv_output("file_dave_indiv_alpha_sap"                , DAILY, "alpha_sap", [](const Individual &i) { return i.alpha_sap; }, "proportion allocated to sapwood (unitless)");
	add_indiv_output("file_dave_indiv_alpha_leaf"               , DAILY, "alpha_leaf", [](const Individual &i) { return i.alpha_leaf; }, "proportion allocated to leaf (unitless)");
	add_indiv_output("file_dave_indiv_alpha_root"               , DAILY, "alpha_root", [](const Individual &i) { return i.alpha_root; }, "proportion allocated to root (unitless)");
	add_indiv_output("file_dave_indiv_alpha_storage"            , DAILY, "alpha_storage", [](const Individual &i) { return i.alpha_storage; }, "proportion allocated to  storage (unitless)");
	add_indiv_output("file_dave_indiv_cgrow"          		    , DAILY, "cgrow", [](const Individual &i) { return i.cgrow; }, "daily total growth rate of individual (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_cgrow_leaf"          		, DAILY, "cgrow_leaf", [](const Individual &i) { return i.cgrow_leaf; }, "daily leaf allocation of individual (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_cgrow_root"         	    , DAILY, "cgrow_root", [](const Individual &i) { return i.cgrow_root; }, "daily root allocation of individual (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_cgrow_sap"          	    , DAILY, "cgrow_sap", [](const Individual &i) { return i.cgrow_sap; }, "daily sapwood allocation of individual (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_cgrow_repr"     		    , DAILY, "cgrow_repr", [](const Individual &i) { return i.cgrow_repr; }, "daily reproduction allocation of individual (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_cgrow_storage"     	    , DAILY, "cgrow_storage", [](const Individual &i) { return i.cgrow_storage; }, "daily storage allocation of individual (kgC/m2/d)");
	add_indiv_output("file_dave_cohort_gpp"     	            , DAILY, "gpp", [](const Individual &i) { return i.dgpp; }, "daily gpp (kgC/m2/d)");
	add_indiv_output("file_dave_cohort_npp"     	            , DAILY, "npp", [](const Individual &i) { return i.dnpp; }, "daily npp (kgC/m2/d)");
	add_indiv_output("file_dave_indiv_litter_pools"             , DAILY, get_pft_litter_pools, "Total Daily PFT-Level C Mass in Litter (kgC/m2)", true);
	add_indiv_output("file_dave_indiv_sapwood_area"             , DAILY, "sapwood_area", [](const Individual &i) { return i.sa; }, "Sapwood area (m2)");
	add_indiv_output("file_dave_indiv_sapwood_depth"            , DAILY, "sapwood_depth", [](const Individual &i) { return i.sapwood_depth; }, "Sapwood depth (m)");
	add_indiv_output("file_dave_closs_hydraulics"               , DAILY, get_closs_hydraulics_cols, "Daily C lost to hydraulic stress (kgC/m2/day)", true);
	add_indiv_output("file_dave_transpiration_check"            , DAILY, transp_check_columns, "Daily check of old vs new transpiration (mm/day)");
	add_indiv_output("file_dave_indiv_closs"                    , DAILY, indiv_closs_columns, "Daily individual-level C loss components (kgC/m2/day)", true);

	// Daily individual-level array outputs.
	add_indiv_array_output("file_dave_layerwise_fpar"           , DAILY, "layer", get_nlayer_canopy, [](const Individual &indiv, int i) { return indiv.canopy[i].fipar; }, "Daily layered fpar [0-1]", mean);
	add_indiv_array_output("file_dave_layerwise_lai"            , DAILY, "layer", get_nlayer_canopy , [](const Individual &indiv, int i) { return indiv.canopy[i].lai;  }, "Daily layered LAI (m2/m2)", sum);
	add_indiv_array_output("file_dave_Rs_r"                     , DAILY, "layer", [](const Patch &p) { return (size_t)NSOILLAYERS_MAX; }, [](const Individual &indiv, int i) { return indiv.Rs_r[i] / 1000.0;  }, "Daily soil-root resistance (mol/m/s/MPa)", sum);
	add_indiv_array_output("file_dave_water_uptake"             , DAILY, "layer", [](const Patch &p) { return (size_t)NSOILLAYERS_MAX; }, [](const Individual &indiv, int i) { return indiv.fwuptake[i] * indiv.aet;  }, "Water uptake from each layer (mm/day)", sum);

	// Daily pft-level per-patch outputs.
	add_pft_output("file_dave_standing_dead"                    , DAILY, get_pft_standing_dead_cmass, "Daily PFT-Level Total Standing Dead C Mass (kgC/m2)", true);
	add_pft_output("file_dave_cmass_standing_repr"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_repr; }, "Standing Dead Reproductive (Seedbank for grasses) C mass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_standing_repr"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_repr; }, "Standing Dead Reproductive (Seedbank for grasses) N mass (kgN/m2)",  true);
	add_pft_output("file_dave_cmass_standing_leaf"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_leaf; }, "Standing Dead Leaf C Mass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_standing_leaf"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_leaf; }, "Standing Dead Leaf N Mass (kgN/m2)",  true);
	add_pft_output("file_dave_cmass_standing_root"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_root; }, "Standing Dead Root C Mass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_standing_root"         , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_root; }, "Standing Dead Root N Mass (kgN/m2)",  true);
	add_pft_output("file_dave_cmass_standing_sap"          , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_sap; }, "Standing Dead Sap C Mass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_standing_sap"          , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_sap; }, "Standing Dead Sap N Mass (kgN/m2)",  true);
	add_pft_output("file_dave_cmass_standing_heart"        , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].cmass_litter_heart; }, "Seedbank cmass (kgC/m2)", true);
	add_pft_output("file_dave_nmass_standing_heart"        , DAILY, [](const Patch &p, const Pft &pft) { return p.pft[pft.id].nmass_litter_heart; }, "Seedbank nmass (kgN/m2)",  true);
	add_pft_output("file_dave_dresp"                       , DAILY, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::RA, pft.id); }, "Total daily pft-level autotrophic Respiration (kgC/m2)", true);
	add_pft_output("file_dave_dcmass_mort"                 , DAILY, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].dcmass_mort_kill; }, "Daily live vegetation C removed at Individual::kill() events (kgC/m2/day)", true);
	add_pft_output("file_dave_cmass_seed_ext"              , DAILY, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].cmass_ext_seed; }, "Total daily pft-level hypothetical seedbank (kgC/m2)", true);
	add_pft_output("file_dave_dncohort"                    , DAILY, get_nindiv, "Get the number of cohorts of this PFT currently established in this patch");
	add_pft_output("file_dave_dmort_cstarve_ntree"         , DAILY, [](const Patch &patch, const Pft &pft) { return patch.pft[pft.id].dmort_cstarve_ntree; }, "Daily PFT-level carbon starvation-related mortality, expressed as fraction of individuals within the cohort which died this day, summed across all indivs of a given PFT.");

	// Subdaily individual-level outputs.
	//                      param name in .ins file    timestep  column name    function which gets number of array elems     function which returns the i-th output value                               Description of output file
	add_indiv_array_output("file_dave_subdaily_an"              , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.An[i];    }, "Subdaily net photosynthesis (mol/m2/s)", sum);
	add_indiv_array_output("file_dave_subdaily_rd"              , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Rd[i];    }, "Subdaily leaf respiration (mol/m2/s)", sum);
	add_indiv_array_output("file_dave_subdaily_anc"             , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anc[i];   }, "Subdaily RuBisCO-limited net photosynthesis (mol/m2/s)", sum);
	add_indiv_array_output("file_dave_subdaily_anj"             , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anj[i];   }, "Subdaily RuBP regeneration-limited photoysnthesis (mol/m2/s)", sum);
	add_indiv_array_output("file_dave_subdaily_gsw"             , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.gsw[i];   }, "Subdaily stomatal conductance to water vapour (mol/m2/s)", sum);
	add_indiv_array_output("file_dave_subdaily_ci"              , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.ci[i];    }, "Subdaily intercellular CO2 concentration (mol/mol)", mean);
	add_indiv_array_output("file_dave_subdaily_vcmax"           , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.vcmax[i]; }, "Subdaily max carboxylation rate (mol/m2/s)", mean);
	add_indiv_array_output("file_dave_subdaily_jmax"            , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.jmax[i];  }, "Subdaily max electron transport rate (mol/m2/s)", mean);
	add_indiv_array_output("file_dave_subdaily_lwp"             , DAILY, COLNAME_TIMESTEP, [](const Patch &p) { return (size_t)date.subdaily; }, [](const Individual &indiv, int i) { return indiv.lwp[i];  }, "Subdaily leaf water potential (MPa)", get_lai /* Mean weighted by LAI*/);

	// Subdaily cohort-level per-layer photosynthesis inputs.
	add_indiv_layer_timestep_output("file_dave_ps_inputs", ps_input_cols(), "Subdaily layered cohort-level photosynthesis inputs");

	// Subdaily cohort-level per-layer photosynthesis outputs.
	add_indiv_layer_timestep_output("file_dave_ps_outputs", ps_output_cols(), "Subdaily layered cohort-level photosynthesis outputs");

	// Subdaily cohort-level per-layer energy balance.
	add_indiv_layer_timestep_output("file_dave_energy_balance", energy_balance_cols(), "Subdaily layered cohort-level energy balance outputs");

	add_indiv_layer_timestep_output("file_dave_ps_diagnostics", ps_diagnostics_cols(), "Subdaily layered cohort-level photosynthesis diagnostics");

	// Daily patch-level array outputs.
	// One row per array element per patch per timestep, multiple data columns (tall/long table).
	add_tall_patch_array_output("file_dave_soil_state"               , DAILY, "layer", [](const Patch &p) { return (size_t)p.soil.nsoillayers; }, soil_state_cols(), "Soil water content as fraction of available water holding capacity");
	add_tall_patch_array_output("file_dave_soil_properties"          , DAILY, "layer", [](const Patch &p) { return (size_t)p.soil.nsoillayers; }, soil_properties_cols(), "Soil layer thickness (m)");

	// One row per patch per timestep, one column per array element (wide table).
	add_wide_patch_array_output("file_dave_sw"                       , DAILY, "sw_", NSOILLAYERS_MAX, [](const Patch &p, int i) { return i >= p.soil.nsoillayers ? 0.0 : p.soil.soillayers[i].wcont(); }, "Soil water fraction full (0-1)", true, false);
	add_wide_patch_array_output("file_dave_swmm"                     , DAILY, "swmm_", NSOILLAYERS_MAX, [](const Patch &p, int i) { return i >= p.soil.nsoillayers ? 0.0 : p.soil.soillayers[i].wmass; }, "Soil water content (mm)", true, false);
	add_wide_patch_array_output("file_dave_swvol"                    , DAILY, "swvol_", NSOILLAYERS_MAX, [](const Patch &p, int i) { return i >= p.soil.nsoillayers ? 0.0 : p.soil.soillayers[i].wmass / p.soil.soillayers[i].thickness; }, "Volumetric soil water content (m3/m3)", false, true);
	add_wide_patch_array_output("file_dave_soil_psi_layer"           , DAILY, "soil_psi", NSOILLAYERS_MAX, [](const Patch &p, int i) { return i >= p.soil.nsoillayers ? 0.0 : p.soil.soillayers[i].psi_mpa(); }, "Soil water potential (MPa) in each layer", false, false);
	add_wide_patch_array_output("file_dave_ksoil_layer"              , DAILY, "ksoil_", NSOILLAYERS_MAX, [](const Patch &p, int i) { return i >= p.soil.nsoillayers ? 0.0 : p.soil.soillayers[i].k_actual / 2.778e-6; }, "Soil hydraulic conductivity (cm/hr)", false, true);

	// Subdaily patch-level outputs.
	add_wide_patch_array_output("file_dave_met_subdaily_temp"        , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.temps[i];    }, "Subdaily temperature (°C)", false, true);
	add_wide_patch_array_output("file_dave_met_subdaily_par"         , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pars[i] * 0.001;     }, "Subdaily PAR (kJ/m2/timestep)", true, false);
	add_wide_patch_array_output("file_dave_met_subdaily_vpd"         , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.vpds[i];     }, "Subdaily VPD (kPa)", false, true);
	add_wide_patch_array_output("file_dave_met_subdaily_insol"       , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.insols[i];   }, "Subdaily insolation (units depend on instype)", true, false);
	add_wide_patch_array_output("file_dave_met_subdaily_precip"      , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.precs[i];    }, "Subdaily precipitation (mm)", true, false);
	add_wide_patch_array_output("file_dave_met_subdaily_pressure"    , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pressure[i] * 0.001; }, "Subdaily air pressure (kPa)", false, true);
	add_wide_patch_array_output("file_dave_met_subdaily_co2"         , DAILY, COLNAME_TIMESTEP, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.co2;         }, "Subdaily atmospheric co2 concentration (ppm)", false, true);

	// Standpft-level outputs.
	add_spft_output("file_dave_spft_cmass_repr"                 , DAILY, [](const Standpft &spft) { return spft.cmass_repr; }, "Daily reproductive cmass at standpft-level (kgC/m2)");
	add_spft_output("file_dave_anetps_ff_max"                   , DAILY, [](const Standpft &spft) { return spft.anetps_ff_max; }, "Maximum recorded annual net photosynthesis at forest-floor (kgC/m2)");

	// Register instruction file parameters.
	declare_parameters();

	total_time = get_wall_time() - t_start;
}

DaveOutput::~DaveOutput() {
	log_diagnostic("Total time spent in dave output: %.2fs\n", total_time);
}

void DaveOutput::init() {
	double t_start = get_wall_time();
	// rows = new OutputRows(output_channel, -1, -1, -1, -1);
	open_outfiles();
	total_time += get_wall_time() - t_start;
}

void DaveOutput::write_outputs(Gridcell &gridcell, std::vector<std::unique_ptr<OutputFile>> &outputs) {
	ensure_channel_open();

	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	// OutputRows rows(output_channel, lon, lat, year, date.day);

	for (std::unique_ptr<OutputFile> &output : outputs) {
		if (output->enabled()) {
			output->write_rows(output_channel, gridcell);
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
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, []() { return 1; }, reader, description, sum, mean);
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
	register_output(output, timestep);
}

void DaveOutput::add_gridcell_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		double (*read_output)(const Gridcell &),
		std::string description,
		bool sum,
		bool mean) {
	ColumnDefinition<Gridcell> column(colname, read_output);
	std::vector<ColumnDefinition<Gridcell>> columns;
	columns.push_back(column);
	GridcellLevelOutput *output = new GridcellLevelOutput(name, columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_gridcell_output(
		std::string name,
		timestep_t timestep,
		std::vector<ColumnDefinition<Gridcell>> (*get_columns)(),
		std::string description,
		bool sum,
		bool mean) {
	GridcellLevelOutput *output = new GridcellLevelOutput(name, get_columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_gridcell_aggregated_indiv_output(
		std::string name,
		timestep_t timestep,
		double (*read_output)(const Individual &i),
		double (*aggregator)(const double *, size_t),
		std::string description,
		bool sum,
		bool mean) {
	GridcellLevelOutput *output = new GridcellLevelOutput(
		name,
		[read_output, aggregator]() -> std::vector<ColumnDefinition<Gridcell>> { return get_gridcell_aggregated_indiv_columns(read_output, aggregator); },
		description,
		sum,
		mean);
	register_output(output, timestep);
}

void DaveOutput::add_gridcell_aggregated_ppft_output(
		std::string name,
		timestep_t timestep,
		double (*read_output)(const Patch &, const Pft &),
		std::string description,
		bool sum,
		bool mean
) {
	GridcellLevelOutput *output = new GridcellLevelOutput(
		name,
		[read_output]() -> std::vector<ColumnDefinition<Gridcell>> { return get_gridcell_patchpft_columns(read_output); },
		description,
		sum,
		mean);
	register_output(output, timestep);
}

void DaveOutput::add_gridcell_aggregated_tree_output(
		std::string name,
		timestep_t timestep,
		double (*read_output)(const Individual &i),
		double (*aggregator)(const double *, size_t),
		std::string description,
		bool sum,
		bool mean) {
	GridcellLevelOutput *output = new GridcellLevelOutput(
		name,
		[read_output, aggregator]() -> std::vector<ColumnDefinition<Gridcell>> { return get_gridcell_aggregated_tree_columns(read_output, aggregator); },
		description,
		sum,
		mean);
	register_output(output, timestep);
}

void DaveOutput::add_wide_patch_array_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		int nelem,
		double (*read_output)(const Patch&, int),
		std::string description,
		bool sum,
		bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, [nelem]() { return nelem; }, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_wide_patch_array_output(
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

void DaveOutput::add_tall_patch_array_output(std::string name,
									         timestep_t timestep,
											 std::string index_col_name,
											 std::function<size_t(const Patch &)> get_nelem,
											 std::vector<ArrayColumnDefinition<Patch>> columns,
											 std::string description) {
    DynamicPatchLevelArrayOutput *output = new DynamicPatchLevelArrayOutput(name, index_col_name, get_nelem, columns, description);
    register_output(output, timestep);
}

void DaveOutput::add_aggregated_indiv_output(std::string name, timestep_t timestep, double (*read_output)(const Individual&), std::string description, double (*aggregator)(const double*, size_t), bool sum, bool mean) {
	OutputReader<Individual> *reader = new OutputReader<Individual>(read_output);
	auto output = new AggregatedIndivOutput(name, reader, aggregator, description, sum, mean);

	register_output(output, timestep);
}

void DaveOutput::add_indiv_output(
		std::string name,
		timestep_t timestep,
		std::vector<ColumnDefinition<Individual>> (*get_columns)(),
		std::string description,
		bool sum,
		bool mean) {
	IndividualOutput *output = new IndividualOutput(name, get_columns, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_output(
		std::string name,
		timestep_t timestep,
		std::string colname,
		double (*read_output)(const Individual &),
		std::string description) {
	std::vector<ColumnDefinition<Individual>> columns;
	columns.push_back(ColumnDefinition<Individual>(colname, read_output));
	IndividualOutput *output = new IndividualOutput(name, columns, description, false, false);
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
		size_t (*get_nelem)(const Patch &),
		double (*read_output)(const Individual&, int),
		std::string description,
		double (*aggregator)(const double *, size_t)) {
	ArrayOutputReader<Individual> *reader = new ArrayOutputReader<Individual>(read_output);
	IndivLevelArrayOutput *output = new IndivLevelArrayOutput(name, elem_name, get_nelem, reader, description, aggregator);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_array_output(
		std::string name,
		timestep_t timestep,
		std::string elem_name,
		size_t (*get_nelem)(const Patch &),
		double (*read_output)(const Individual&, int),
		std::string description,
		double (*read_weight)(const Individual&)) {
	ArrayOutputReader<Individual> *reader = new ArrayOutputReader<Individual>(read_output);
	OutputReader<Individual> *weight_reader = new OutputReader<Individual>(read_weight);
	IndivLevelArrayOutput *output = new IndivLevelArrayOutput(name, elem_name, get_nelem, reader, description, weight_reader);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_layer_timestep_output(
		std::string name,
		std::vector<IndivLayerTimestepOutput::Column> columns,
		std::string description) {
	IndivLayerTimestepOutput *output = new IndivLayerTimestepOutput(name, columns, description);
	register_output(output, DAILY);
}

void DaveOutput::add_indiv_layer_timestep_output(
		std::string name,
		std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> column_defs,
		std::string description) {
	std::vector<IndivLayerTimestepOutput::Column> columns;
	for (auto def : column_defs) {
		IndivLayerTimestepOutput::Column col;
		col.name = def.first;
		col.read = def.second;
		columns.push_back(col);
	}
	IndivLayerTimestepOutput *output = new IndivLayerTimestepOutput(name, columns, description);
	register_output(output, DAILY);
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

bool DaveOutput::OutputFile::should_write(OutputChannel *channel) {
	if (!should_write_today()) {
		return false;
	}
	return true;
}

void DaveOutput::OutputFile::write_output(OutputChannel *channel, double value)
{
	assert_table_open();
	if (should_write(channel)) {
		channel->add_value(*table, value);
	}
}

void DaveOutput::OutputFile::write_output(OutputChannel *channel, int value)
{
	assert_table_open();
	if (should_write(channel)) {
		channel->add_value(*table, value);
	}
}

void DaveOutput::OutputFile::write_output(OutputChannel *channel, const char *value)
{
	assert_table_open();
	if (should_write(channel)) {
		channel->add_value(*table, value);
	}
}

void DaveOutput::OutputFile::finish_row(OutputChannel *channel, double lon, double lat, int year, int day) {
	assert_table_open();
	channel->finish_row(*table, lon, lat, year, day);
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
// GridcellLevelOutput implementation.
// #############################################################################

DaveOutput::GridcellLevelOutput::~GridcellLevelOutput() { }

DaveOutput::GridcellLevelOutput::GridcellLevelOutput(
		std::string param,
		std::function<std::vector<ColumnDefinition<Gridcell>>()> get_columns,
		std::string desc,
		bool sum,
		bool mean) : OutputFile(param, desc) {
	get_column_definitions = get_columns;
	write_sum = sum;
	write_mean = mean;
}

DaveOutput::GridcellLevelOutput::GridcellLevelOutput(
		std::string param,
		std::vector<ColumnDefinition<Gridcell>> columns,
		std::string desc,
		bool sum,
		bool mean) : OutputFile(param, desc) {
	get_column_definitions = NULL;
	this->columns = columns;
	write_sum = sum;
	write_mean = mean;
}

void DaveOutput::GridcellLevelOutput::write_rows(
		OutputChannel *channel,
		Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

	double total = 0;
	for (ColumnDefinition<Gridcell> &column : columns) {
		double value = column.read_value(gridcell);
		write_output(channel, value);
		total += value;
	}

	if (write_sum) {
		write_output(channel, total);
	}
	if (write_mean) {
		write_output(channel, total / columns.size());
	}
	finish_row(channel, lon, lat, year, day);
}

ColumnDescriptors DaveOutput::GridcellLevelOutput::get_cols() {
	ColumnDescriptors descriptors;

	// Only invoke get_column_definitions if necessary.
	if (get_column_definitions) {
		columns = get_column_definitions();
	}

	for (const ColumnDefinition<Gridcell> &column : columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		descriptors += ColumnDescriptor(name.c_str(), width, NUMERIC_PRECISION);
	}

	if (write_sum) {
		descriptors += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	if (write_mean) {
		descriptors += ColumnDescriptor(COLNAME_MEAN, COLUMN_WIDTH, NUMERIC_PRECISION);
	}

	return descriptors;
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
		OutputChannel *channel,
		Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

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
			write_output(channel, stand.id);

			// Write stand data.
			double total = 0;
			for (double x : stand_data) {
				write_output(channel, x);
				total += x;
			}
			if (write_sum) {
				write_output(channel, total);
			}
			if (write_mean) {
				write_output(channel, total / stand_data.size());
			}
			finish_row(output_channel, lon, lat, year, day);
		}
	}

	if (write_stats() && data.size()) {
		// Each element of data contains all columns for each stand. We want
		// to write one column at a time.

		// First, get the number of columns.
		size_t ncol = data[0].size();

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
			write_output(channel, mean);
			write_output(channel, median);
			write_output(channel, sd);
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
				write_output(channel, mean);
				write_output(channel, median);
				write_output(channel, sd);
			}
			if (write_mean) {
				write_output(channel, mean / ncol);
				write_output(channel, median / ncol);
				write_output(channel, sd / ncol);
			}
		}
		finish_row(output_channel, lon, lat, year, day);
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

void DaveOutput::PatchLevelOutput::write_rows(OutputChannel *channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;

		stand.firstobj();

		std::vector<std::vector<double>> data;

		while (stand.isobj) {
			Patch &patch = stand.getobj();

			// Write stand ID; when writing "standard" outputs, we write stand
			// ID once per patch (ie once per row of output). When writing
			// stats, we only write one row per stand, so we cannot write stand
			// ID from inside a patch-level loop.
			if (!write_stats()) {
				write_output(channel, stand.id);
			}

			std::vector<double> patch_data = write_patch(patch, channel);
			if (write_stats()) {
				data.push_back(patch_data);
			} else {
				// The code would be cleaner if we were to use the data 2d array
				// in both cases, but I suspect that performance will be better
				// with this shortcut when the number of patches is high.

				// Write patch ID.
				write_output(channel, patch.id);

				// Write patch data.
				double total = 0;
				for (double x : patch_data) {
					write_output(channel, x);
					total += x;
				}
				if (write_sum) {
					write_output(channel, total);
				}
				if (write_mean) {
					write_output(channel, patch_data.size() ? total / patch_data.size() : 0);
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
			size_t ncol = data[0].size();

			// We write one row of output per stand, so we need to first output
			// the stand ID.
			write_output(channel, stand.id);

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
				write_output(channel, mean);
				write_output(channel, median);
				write_output(channel, sd);

			}
			if (write_sum || write_mean) {
				std::vector<double> totals;
				// Iterate over patches.
				for (size_t i = 0; i < data.size(); i++) {
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
					write_output(channel, mean);
					write_output(channel, median);
					write_output(channel, sd);
				}
				if (write_mean) {
					write_output(channel, mean / ncol);
					write_output(channel, median / ncol);
					write_output(channel, sd / ncol);
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

DaveOutput::IndividualOutput::IndividualOutput(
		std::string param_name,
		std::vector<ColumnDefinition<Individual>> columns,
		std::string description,
		bool sum,
		bool mean) : OutputFile(param_name, description) {
	this->columns = columns;
	get_column_definitions = NULL;
	write_sum = sum;
	write_mean = mean;
}

DaveOutput::IndividualOutput::~IndividualOutput() { }

void DaveOutput::IndividualOutput::write_rows(OutputChannel *channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

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
				write_output(channel, stand.id);

				// Write patch ID.
				write_output(channel, patch.id);

				// Write individual ID.
				write_output(channel, indiv.id);

				// Write PFT name.
				write_output(channel, indiv.pft.name);

				// Write data.
				for (size_t j = 0; j < columns.size(); j++) {
					write_output(channel, columns[j].read_value(indiv));
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
	if (get_column_definitions) {
		columns = get_column_definitions();
	}
	if (columns.size() == 0) {
		fail("No columns specified for output file");
	}

	ColumnDescriptors descriptors;
	descriptors += ColumnDescriptor(COLNAME_STAND, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_INDIV_ID, COLUMN_WIDTH, 0);
	descriptors += ColumnDescriptor(COLNAME_PFT, COLUMN_WIDTH + 2, 0, TEXT);

	for (const ColumnDefinition<Individual> &column : columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		descriptors += ColumnDescriptor(name.c_str(), width, NUMERIC_PRECISION);
	}

	return descriptors;
}

// #############################################################################
// IndivLayerTimestepOutput implementation.
// #############################################################################

DaveOutput::IndivLayerTimestepOutput::IndivLayerTimestepOutput(
		std::string param_name,
		std::vector<Column> columns,
		std::string description)
	: OutputFile(param_name, description), columns(columns) { }

void DaveOutput::IndivLayerTimestepOutput::write_rows(OutputChannel *channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand &stand = *gc_itr;
		stand.firstobj();

		while (stand.isobj) {
			const Patch &patch = stand.getobj();
			const Vegetation &vegetation = patch.vegetation;

			for (size_t i = 0; i < vegetation.nobj; i++) {
				const Individual &indiv = vegetation[i];

				int ntimesteps = (int)indiv.ps_inputs_layerwise.size();
				for (int t = 0; t < ntimesteps; t++) {
					int nlayers = (int)indiv.ps_inputs_layerwise[t].size();
					for (int l = 0; l < nlayers; l++) {
						write_output(channel, stand.id);
						write_output(channel, patch.id);
						write_output(channel, indiv.id);
						write_output(channel, indiv.pft.name);
						write_output(channel, t);
						write_output(channel, l);

						for (const Column &col : columns) {
							write_output(channel, col.read(indiv, t, l));
						}

						finish_row(channel, lon, lat, year, day);
					}
				}
			}

			stand.nextobj();
		}

		gc_itr++;
	}
}

ColumnDescriptors DaveOutput::IndivLayerTimestepOutput::get_cols() {
	ColumnDescriptors descriptors;
	descriptors += ColumnDescriptor(COLNAME_STAND,    COLUMN_WIDTH,     0);
	descriptors += ColumnDescriptor(COLNAME_PATCH,    COLUMN_WIDTH,     0);
	descriptors += ColumnDescriptor(COLNAME_INDIV_ID, COLUMN_WIDTH,     0);
	descriptors += ColumnDescriptor(COLNAME_PFT,      COLUMN_WIDTH + 2, 0, TEXT);
	descriptors += ColumnDescriptor(COLNAME_TIMESTEP, COLUMN_WIDTH,     0);
	descriptors += ColumnDescriptor(COLNAME_LAYER,    COLUMN_WIDTH,     0);

	for (const Column &col : columns) {
		int width = max((size_t)COLUMN_WIDTH, col.name.length() + 1);
		descriptors += ColumnDescriptor(col.name.c_str(), width, NUMERIC_PRECISION);
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

void add_pft_columns(ColumnDescriptors &columns) {
	for (size_t i = 0; i < pftlist.nobj; i++) {
		add_column(columns, (const char *)pftlist[i].name, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
}

ColumnDescriptors DaveOutput::PftOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	add_pft_columns(columns);
	return columns;
}

std::vector<double> DaveOutput::PftOutput::write_patch(Patch &patch, OutputChannel *channel) {
	std::vector<double> results(pftlist.nobj);

	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];

		double result = get_value(patch, pft);
		results[i] = result;
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

std::vector<double> DaveOutput::CustomPatchLevelOutput::write_patch(Patch &patch, OutputChannel *channel) {
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

std::vector<double> DaveOutput::StandPftOutput::write_patch(Patch &patch, OutputChannel *channel) {
	Vegetation &vegetation = patch.vegetation;
	Stand &stand = patch.stand;

	std::vector<double> results;

	// Note that patch id is written before this function is called, so we
	// shouldn't write it ourselves (that would result in a 'too many columns'
	// error).

	// Write one column per PFT.
	for (size_t i = 0; i < pftlist.nobj; i++) {
		Pft &pft = pftlist[i];
		const Standpft &spft = stand.pft[pft.id];
		double result = read_output->read(spft);
		results.push_back(result);
	}

	return results;
}

ColumnDescriptors DaveOutput::StandPftOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	add_pft_columns(columns);
	return columns;
}

// #############################################################################
// PatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		std::function<int()> get_n,
		ArrayOutputReader<Patch> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {

	get_nelem = get_n;
	col_name_prefix = col_name;
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
}

DaveOutput::PatchLevelArrayOutput::~PatchLevelArrayOutput() { }

std::vector<double> DaveOutput::PatchLevelArrayOutput::write_patch(
		Patch &patch, OutputChannel *channel) {
	std::vector<double> results;
	for (int i = 0; i < get_nelem(); i++) {
		double value = read_output->read(patch, i);
		results.push_back(value);
	}
	return results;
}

ColumnDescriptors DaveOutput::PatchLevelArrayOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();

	int nelem = get_nelem();
	for (int i = 0; i < nelem; i++) {
		std::string name = nelem == 1 ? col_name_prefix : col_name_prefix + std::to_string(i);
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(columns, name.c_str(), width, NUMERIC_PRECISION);
	}

	return columns;
}

// #############################################################################
// DynamicPatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::DynamicPatchLevelArrayOutput::DynamicPatchLevelArrayOutput(
		std::string param_name,
		std::string index_col_name,
		std::function<size_t(const Patch &)> get_nelem,
		std::vector<ArrayColumnDefinition<Patch>> columns,
		std::string description) : OutputFile(param_name, description) {
	this->index_col_name = index_col_name;
	this->columns = columns;
	this->get_nelem = get_nelem;
}

DaveOutput::DynamicPatchLevelArrayOutput::~DynamicPatchLevelArrayOutput() {

}

void DaveOutput::DynamicPatchLevelArrayOutput::write_rows(OutputChannel *channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	for (size_t i = 0; i < gridcell.size(); i++) {
		const Stand &stand = gridcell[i];
		for (size_t j = 0; j < stand.nobj; j++) {
			const Patch &patch = stand[j];

			int nelem = get_nelem(patch);
			for (size_t k = 0; k < nelem; k++) {
				write_output(channel, stand.id);

				if (!write_stats()) {
					write_output(channel, patch.id);
				}

				write_output(channel, (int)k);
				for (const ArrayColumnDefinition<Patch> &column : columns) {
					double value = column.read_value(patch, k);
					write_output(channel, value);
				}
				finish_row(channel, lon, lat, year, day);
			}
		}
	}
}

ColumnDescriptors DaveOutput::DynamicPatchLevelArrayOutput::get_cols() {
	ColumnDescriptors descriptors = create_column_descriptors();

	add_column(descriptors, index_col_name.c_str(), COLUMN_WIDTH, 0);

	for (const ArrayColumnDefinition<Patch> &column : columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(descriptors, name.c_str(), width, NUMERIC_PRECISION);
	}

	return descriptors;
}

// #############################################################################
// IndivLevelArrayOutput implementation.
// #############################################################################

DaveOutput::IndivLevelArrayOutput::IndivLevelArrayOutput(
		std::string param_name,
		std::string element_name,
		size_t (*get_nelem)(const Patch &),
		ArrayOutputReader<Individual> *get_output,
		std::string description,
		double (*agg)(const double*, size_t)) : OutputFile(param_name, description)
{
	elem_name = element_name;
	get_n = get_nelem;
	read_elem = std::unique_ptr<ArrayOutputReader<Individual>>(get_output);
	aggregator = agg;
	read_weight = nullptr;
}

DaveOutput::IndivLevelArrayOutput::IndivLevelArrayOutput(
		std::string param_name,
		std::string element_name,
		size_t (*get_nelem)(const Patch &),
		ArrayOutputReader<Individual> *get_output,
		std::string description,
		OutputReader<Individual> *weight_reader) : OutputFile(param_name, description)
{
	elem_name = element_name;
	get_n = get_nelem;
	read_elem = std::unique_ptr<ArrayOutputReader<Individual>>(get_output);
	aggregator = NULL;
	read_weight = std::unique_ptr<OutputReader<Individual>>(weight_reader);
}

void DaveOutput::IndivLevelArrayOutput::write_rows(OutputChannel *channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;

		for (size_t j = 0; j < stand.nobj; j++) {
			Patch &patch = stand[j];

			int nelem = get_n(patch);
			Vegetation &vegetation = patch.vegetation;
			for (int i = 0; i < nelem; i++) {

				// Write patch ID.
				write_output(channel, patch.id);

				// Write array index.
				write_output(channel, i);

				vegetation.firstobj();
				double total = 0;

				for (size_t k = 0; k < pftlist.nobj; k++) {
					std::vector<Individual *> indivs = get_indivs(vegetation, (const char *)pftlist[k].name);

					double aggregated = 0;
					if (read_weight) {
						// Compute weighted mean, using user-defined weight.
						double total_weight = 0;
						double weighted_sum = 0;
						for (Individual *indiv : indivs) {
							double value = read_elem->read(*indiv, i);
							double weight = read_weight->read(*indiv);
							weighted_sum += value * weight;
							total_weight += weight;
						}
						aggregated = total_weight > 0 ? weighted_sum / total_weight : 0;
					} else {
						// Aggregate using custom aggregator.
						std::vector<double> values;
						values.reserve(indivs.size());
						for (Individual *indiv : indivs) {
							double value = read_elem->read(*indiv, i);
							values.push_back(value);
						}
						aggregated = values.size() && aggregator ? aggregator(&values.front(), values.size()) : 0;
					}
					total += aggregated;

					// Write current array element for this individual.
					write_output(channel, aggregated);
				}

				// Write total value across all individuals for this array elem.
				write_output(channel, total);

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
	// TODO: add support for stats outputs in the indiv-level array outputs.
	// add_pft_columns(columns);
	for (size_t i = 0; i < pftlist.nobj; i++) {
		columns += ColumnDescriptor((const char *)pftlist[i].name, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// AggregatedIndivOutput implementation.
// #############################################################################

DaveOutput::AggregatedIndivOutput::AggregatedIndivOutput(
		std::string param_name,
		OutputReader<Individual> *get_output,
		double (*agg)(const double *, size_t),
		std::string description,
		bool sum,
		bool mean) : PftOutput(param_name, NULL, description, sum, mean) {
	read_output = std::unique_ptr<OutputReader<Individual>>(get_output);
	aggregator = agg;
}

double DaveOutput::AggregatedIndivOutput::get_value(
		const Patch &patch,
		const Pft &pft) const {

	const Vegetation &vegetation = patch.vegetation;
	std::vector<double> data(vegetation.nobj);
	size_t index = 0;
	for (size_t i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		if (indiv.pft.id == pft.id) {
			double value = read_output->read(indiv);
			data[index++] = value;
		}
	}
	// index now has the number of elements in data (the size and capacity of
	// data are the number of cohorts in the patch).
	return data.size() ? aggregator(&data.front(), index) : 0;
}

bool patch_output_stats;

const int COLUMN_WIDTH = 15; // Enough to fit a 6-digit number
const int COLUMN_WIDTH_PFT = 24;
const int NUMERIC_PRECISION = 8;
const char *COLNAME_TOTAL = "total";
const char *COLNAME_MEAN = "mean";
const char *COLNAME_PATCH = "patch";
const char *COLNAME_STAND = "stand";
const char *COLNAME_TIMESTEP = "timestep";

} // namespace GuessOutput
