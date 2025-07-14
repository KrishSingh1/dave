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
Get the patch-level daily evapotranspiration (mm/day).

@param patch: The patch.
*/
double get_et(const Patch &patch);

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
Get the columns used in the patch-level runoff output file.
*/
std::vector<ColumnDefinition<Patch>> runoff_columns();

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
Get the columns used in the patch-level npool file.
*/
std::vector<ColumnDefinition<Patch>> npool_columns();

/*
Get the columns used in the patch-level nflux file.
*/
std::vector<ColumnDefinition<Patch>> nflux_columns();

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
Get the columns used in the patch-level nitrogen input output file.
*/
std::vector<ColumnDefinition<Patch>> ninput_columns();

/*
Get the columns used in the patch-level evaporation output file.
*/
std::vector<ColumnDefinition<Patch>> evap_columns();

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
std::vector<ColumnDefinition<Patch>> aet_columns();

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
Get the litter cmass for the specified pft in the specified patch (kgC/m2).
*/
double get_pft_cmass_litter(const Patch &patch, const Pft &pft);

/*
Get the first individual in the patch, if any are established.

@param patch: A patch object.
*/
Individual *get_first_indiv(Patch &patch);

/*
Get the number of canopy layers used for fpar calculations today in the
specified patch.

@param patch: A patch object.
*/
int get_nlayer_fpar(const Patch &patch);

/*
Get the number of canopy layers used for light interception today in the
specified patch.

@param patch: A patch object.
*/
int get_nlayer_lai(const Patch &patch);

int get_nlayer_an(const Patch &patch);

/*
Get the litter cmass for the specified pft in the specified patch (kgC/m2).
*/
double get_pft_cmass_litter(const Patch &patch, const Pft &pft);

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
Generate a vector of column defintions, which aggregate the specified value read
from an individual up to the gridcell level, for all tree PFTs.
*/
std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_tree_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t)
);

} // namespace GuessOutput

#endif // LPJ_GUESS_DAVE_OUTPUT_HELPERS_H
