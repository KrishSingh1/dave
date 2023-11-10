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
@param i: Unused (required by API).
*/
double get_er(const Patch &patch, int i);

/*
Get the patch-level ecosystem GPP.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_gpp(const Patch &patch, int i);

/*
Get the patch-level ecosystem NPP.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_npp(const Patch &patch, int i);

/*
Get the patch-level ecosystem NEE.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_nee(const Patch &patch, int i);

/*
Get the patch-level daily evapotranspiration (mm/day).

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_et(const Patch &patch, int i);

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
double sum(double *array, int n);

/*
Calculate a mean.

@param array: An array of data.
@param n: Length of the array.
*/
double mean(double *array, int n);

/*
Get the columns used in the patch-level runoff output file.
*/
std::vector<ColumnDefinition<Patch>> runoff_columns();

/*
Get the columns used in the patch-level cpool file.
*/
std::vector<ColumnDefinition<Patch>> cpool_columns();

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
Get the columns used in the met year output file.
*/
std::vector<ColumnDefinition<Patch>> met_year_columns();

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
Get the fpar value for a given individual in the specified canopy layer.

@param indiv: The individual.
@param i: The layer index - ie the i-th layer above the ground.
*/
double get_fpar_layer(const Individual &indiv, int i);

/*
Get the lai value for a given individual in the specified canopy layer.

@param indiv: The individual.
@param i: The layer index - ie the i-th layer above the ground.
*/
double get_lai_layer(const Individual &indiv, int i);

/*
Get the number of canopy layers used for light interception today in the
specified patch.

@param patch: A patch object.
*/
int get_nlayer_lai(const Patch &patch);

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

} // namespace GuessOutput

#endif // LPJ_GUESS_DAVE_OUTPUT_HELPERS_H