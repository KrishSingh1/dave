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
/// Updates yearly lai and fpc for cropland
void allometry_crop(Individual& indiv);

#endif // CROPALLOCATION_H