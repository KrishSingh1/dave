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
double get_lai(Individual &indiv);

/*
Get the patch-level ecosystem ER.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_er(Patch &patch, int i);

/*
Get the patch-level ecosystem GPP.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_gpp(Patch &patch, int i);

/*
Get the patch-level ecosystem NPP.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_npp(Patch &patch, int i);

/*
Get the patch-level ecosystem NEE.

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_nee(Patch &patch, int i);

/*
Get the patch-level daily evapotranspiration (mm/day).

@param patch: The patch.
@param i: Unused (required by API).
*/
double get_et(Patch &patch, int i);

/*
Get the individual in the patch with the specified name, or NULL if not found.

@param vegetation: The vegetation in the patch.
@param name: Name of a PFT.
*/
std::vector<Individual *> get_indivs(Vegetation &vegetation, const char *name);

/*
Get the name of all known PFTs.
*/
std::vector<std::string> get_pfts();

/*
Get the names of all patch-level fluxes.
*/
std::vector<std::string> get_patch_flux_names();

} // namespace GuessOutput

#endif // LPJ_GUESS_DAVE_OUTPUT_HELPERS_H