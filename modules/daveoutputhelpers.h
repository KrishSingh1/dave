#ifndef LPJ_GUESS_DAVE_OUTPUT_HELPERS_H
#define LPJ_GUESS_DAVE_OUTPUT_HELPERS_H

#include <vector>

#include "guess.h"

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

} // namespace GuessOutput

#endif // LPJ_GUESS_DAVE_OUTPUT_HELPERS_H