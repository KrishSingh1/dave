#ifndef LPJ_GUESS_DAVE_OUTPUT_HELPERS_H
#define LPJ_GUESS_DAVE_OUTPUT_HELPERS_H

#include <vector>

#include "guess.h"
#include "daveoutput.h"

namespace GuessOutput {

/*
Ensure the output channel is open for writing. Abort the simualtion if not.
*/
void ensure_channel_open();

/*
This is an example of a helper function which may be used instead of a lambda
expression when defining an OutputFile.

@param indiv: An individual object..
*/
double get_lai(const Individual &indiv);

/*
Get the patch-level ecosystem ER.

@param patch: The patch.
*/
double get_er(const Patch &patch);

/*
Get the patch-level ecosystem GPP.

@param patch: The patch.
*/
double get_gpp(const Patch &patch);

/*
Get the patch-level ecosystem NPP.

@param patch: The patch.
*/
double get_npp(const Patch &patch);

/*
Get the patch-level ecosystem NEE.

@param patch: The patch.
*/
double get_nee(const Patch &patch);

/*
Get the patch-level annual NEE in kgC/m2.

@param patch: The patch.
*/
double get_anee(const Patch &patch);

/*
Get the patch-level daily evapotranspiration (mm/day).

@param patch: The patch.
*/
double get_et(const Patch &patch);

/*
Get the total basal area in the patch (m2/ha).
*/
double get_patch_ba(const Patch &patch);

/*
Get the individual in the patch with the specified name, or NULL if not found.

@param vegetation: The vegetation in the patch.
@param name: Name of a PFT.
*/
std::vector<Individual *> get_indivs(const Vegetation &vegetation, const char *name);

/*
Get the name of all known PFTs.

@param remove_whitespace: True to remove whitespace. False otherwise.
*/
std::vector<std::string> get_pfts(bool remove_whitespace);

/*
Get the names of all patch-level fluxes.
*/
std::vector<std::string> get_patch_flux_names();

/*
Calculate the sum of an array.

@param array: An array of data.
@param n: Length of the array.
*/
double sum(const double *array, size_t n);

/*
Calculate a mean.

@param array: An array of data.
@param n: Length of the array.
*/
double mean(const double *array, size_t n);

/*
Sum a quantity across soil layers.
*/
double sum_layers(const Soil &soil, double (*read_layer)(const Soillayer &));

/*
Sum a quantity across soil layers.
*/
double mean_layers(const Soil &soil, double (*read_layer)(const Soillayer &));

/*
Get the columns used in the patch-level runoff output file.
*/
std::vector<ColumnDefinition<Patch>> annual_runoff_columns();

/*
Get the columns used in the patch-level cpool file.
*/
std::vector<ColumnDefinition<Patch>> cpool_columns();

/*
Get the columns used in the vegetation C pool output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_cpool_columns();


std::vector<ColumnDefinition<Individual>> indiv_dbhcoloumns_columns();


/*
Get the columns used in the vegetation C pool output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_cpool_tree_columns();

/*
Get the columns used in the patch-level cflux file.
*/
std::vector<ColumnDefinition<Patch>> cflux_columns();

/*
Get the columns used in the patch-level nflux output file.
*/
std::vector<ColumnDefinition<Patch>> nflux_columns();

/*
Get the columns used in the patch-level npool file.
*/
std::vector<ColumnDefinition<Patch>> npool_columns();

/*
Get the columns used in the patch-level pft litter pools file.
*/
std::vector<ColumnDefinition<Individual>> get_pft_litter_pools();

/*
Get the columns used in the indiv-level transpiration debug file.
*/
std::vector<ColumnDefinition<Individual>> transp_check_columns();

/*
Get the columns used in the patch-level closs hydraulics output file.
*/
std::vector<ColumnDefinition<Individual>> get_closs_hydraulics_cols();

/*
Get the columns used in the patch-level npool veg file.
*/
std::vector<ColumnDefinition<Patch>> npool_veg_columns();

/*
Get the columns used in the patch-level globfirm output file.
*/
std::vector<ColumnDefinition<Patch>> globfirm_columns();

/*
Get the columns used in the patch-level mwcont_upper output file.
*/
std::vector<ColumnDefinition<Patch>> mwcont_upper_columns();

/*
Get the columns used in the patch-level mwcont_lower output file.
*/
std::vector<ColumnDefinition<Patch>> mwcont_lower_columns();

/*
Get the columns used in the patch-level PET output file.
*/
std::vector<ColumnDefinition<Patch>> pet_columns();

/*
Get the columns used in the patch-level simfire output file.
*/
std::vector<ColumnDefinition<Patch>> simfire_columns();

/*
Get the columns used in the patch-level fire fuel output file.
*/
std::vector<ColumnDefinition<Patch>> fuel_columns();

/*
Get the columns used in the patch-level coarse woody debris output file.
*/
std::vector<ColumnDefinition<Patch>> coarse_woody_debris_columns();

/*
Get the columns used in the patch-level SOM pool cmass output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_cmass_columns();

/*
Get the columns used in the patch-level SOM pool nmass output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_nmass_columns();

/*
Get the columns used in the patch-level SOM pool nimmob output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_nimmob_columns();

/*
Get the columns used in the patch-level SOM pool nmin_gross output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_nmin_gross_columns();

/*
Get the columns used in the patch-level SOM pool decay reduction output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_decay_reduction_columns();

/*
Get the columns used in the patch-level SOM pool fraction remaining output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_frac_remaining_columns();

/*
Get the columns used in the patch-level SOM pool N to C ratio output file.
*/
std::vector<ColumnDefinition<Patch>> sompool_ntoc_columns();

/*
Get the columns used in the patch-level C litter output file.
*/
std::vector<ColumnDefinition<Patch>> clitter_columns();

/*
Get the columns used in the patch-level N litter output file.
*/
std::vector<ColumnDefinition<Patch>> nlitter_columns();

/*
Get the columns used in the patch-level nitrogen input output file.
*/
std::vector<ColumnDefinition<Patch>> ninput_columns();

/*
Get the columns used in the patch-level evaporation output file.
*/
std::vector<ColumnDefinition<Patch>> daily_evap_columns();

/*
Get the columns used in the patch-level soil C mass output file.
*/
std::vector<ColumnDefinition<Patch>> soil_cmass_columns();

/*
Get the columns used in the patch-level soil N mass output file.
*/
std::vector<ColumnDefinition<Patch>> soil_nmass_columns();

/*
Get the columns used in the patch-level soil N flux output file.
*/
std::vector<ColumnDefinition<Patch>> soil_nflux_columns();

/*
Get the columns used in the annual patch-level nitrogen input output file.
*/
std::vector<ColumnDefinition<Patch>> aninput_columns();

/*
Get the columns used in the met year output file.
*/
std::vector<ColumnDefinition<Patch>> met_year_columns();

/*
Get the columns used in the aet output file.
*/
std::vector<ColumnDefinition<Patch>> daily_aet_columns();

/*
Get the columns used in the fevap output file.
*/
std::vector<ColumnDefinition<Patch>> fevap_columns();

/*
Get the columns used to evaluate the hydraulics model.
*/
std::vector<ColumnDefinition<Individual>> indiv_hydraulics_columns();

/*
Get the columns used to evaluate the hydraulics model.
*/
std::vector<ColumnDefinition<Individual>> indiv_resp_columns();

/*
Get the columns in the soil state output file.
*/
std::vector<ArrayColumnDefinition<Patch>> soil_state_cols();

/*
Get the columns in the soil properties output file.
*/
std::vector<ArrayColumnDefinition<Patch>> soil_properties_cols();


/*
Get the total standing dead biomass for the given PFT.
*/
double get_pft_standing_dead_cmass(const Patch &patch, const Pft &pft);

/*
Get the number of canopy layers used for light interception today in the
specified patch.

@param patch: A patch object.
*/
size_t get_nlayer_canopy(const Patch &patch);

/*
Get the litter cmass for the specified pft in the specified patch (kgC/m2).
*/
double get_pft_standing_dead_cmass(const Patch &patch, const Pft &pft);

/*
Get the soil moisture in mm between the specified top and bottom layers.

@param soil: A soil object.
@param start_layer: The index of the top layer.
@param end_layer: The index of the bottom layer.
*/
double get_soil_moisture(const Soil &soil, double start_depth, double end_depth);

/*
Get the available soil moisture in mm between the specified top and bottom
layers. This is the sum of AWC * wcont for each layer.

@param soil: A soil object.
@param start_layer: The index of the top layer.
@param end_layer: The index of the bottom layer.
*/
double get_available_soil_moisture(const Soil &soil, double start_depth, double end_depth);

/*
Get an estimate of total soil moisture in mm from 1-6m.

@param soil: A soil object.
*/
double get_6m_soil_moisture_estimate(const Soil &soil);

/*
Compute summary statistics for the given row of data.

@param data: A row of output data. This will be mutated (sorted) by this function.
@param mean: (out parameter): The mean of the data.
@param sd: (out parameter): The standard deviation of the data.
@param median: (out parameter): The median of the data.
*/
void compute_stats(std::vector<double> data, double *mean, double *sd, double *median);

/*
Get the number of cohorts of this PFT currently established in this patch.
*/
double get_nindiv(const Patch &patch, const Pft &pft);

/*
Get the columns used in the vegetation C pool output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_cpool_columns();

/*
Get the columns used in the vegetation C pool output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_cpool_tree_columns();

/*
Get the columns used in the vegetation N pool output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_npool_columns();

/*
Get the columns used in the lai output file.
*/
std::vector<ColumnDefinition<Individual>> indiv_lai_columns();

/*
Aggregate a value across all individuals in a patch.
*/
double aggregate_indivs(const Patch &patch, const Pft &pft, double (*read_value)(const Individual &), double (*aggregate)(double *, int));

std::vector<ColumnDefinition<Individual>> indiv_dbhcoloumns_columns();

int get_ncohort(const Patch &patch, const Pft &pft);
/*
Aggregate a patch-level value to gridcell-level, by reading it for all patches
in the gridcell and aggregating accordingly.
*/
double aggregate_to_gridcell(const Gridcell &gridcell, double (*read_value)(const Patch &p));

/*
Generate a vector of column defintions, which aggregate the specified value read
from an individual up to the gridcell level, for all PFTs.
*/
std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_indiv_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t)
);

/*
Generate a vector of column definitions, which aggregate the specified value
read from a pft in a given patch up to the gridcell level, for all PFTs.
*/
std::vector<ColumnDefinition<Gridcell>> get_gridcell_patchpft_columns(
	double (*read_value)(const Patch &, const Pft &)
);

/*
Generate a vector of column defintions, which aggregate the specified value read
from an individual up to the gridcell level, for all tree PFTs.
*/
std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_tree_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t)
);

/*
Get the columns used in the annual water budget output file.
*/
std::vector<ColumnDefinition<Patch>> water_budget_cols_annual();

/*
Get the total standing dead C mass in the patch (kgC/m2).
*/
double get_standing_c(const Patch &patch);

/*
Get the total standing dead N mass in the patch (kgN/m2).
*/
double get_standing_n(const Patch &patch);

/*
Get the columns used in the daily water budget output file.
*/
std::vector<ColumnDefinition<Patch>> water_budget_cols_daily();

/*
Get the columns used in the annual nitrogen budget output file.
*/
std::vector<ColumnDefinition<Patch>> nitrogen_budget_cols_annual();

/*
Get the columns used in the daily nitrogen budget output file.
*/
std::vector<ColumnDefinition<Patch>> nitrogen_budget_cols_daily();

/*
Get the columns used in the mortality output file.
*/
std::vector<ColumnDefinition<Individual>> mortality_cols();

/*
Get the columns used in the high-fidelity photosynthesis inputs file.
*/
std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> ps_input_cols();

/*
Get the columns used in the high-fidelity photosynthesis outputs file.
*/
std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> ps_output_cols();

/*
Get the columns used in the high-fidelity energy balance output file.
*/
std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> energy_balance_cols();

} // namespace GuessOutput

#endif // LPJ_GUESS_DAVE_OUTPUT_HELPERS_H
