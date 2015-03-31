////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cropallocation.h
/// \brief Crop allocation and growth				
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CROPALLOCATION_H
#define CROPALLOCATION_H

/// Updates patch.members fpc_total and fpc_rescale for crops (to be called after crop_phenology())
void update_patch_fpc(Patch& patch);
/// Handles daily crop allocation and daily lai calculation
void growth_daily(Patch& patch);
/// Turnover function for continuous grass.
//void turnover_grass(Individual& indiv);

#endif // CROPALLOCATION_H