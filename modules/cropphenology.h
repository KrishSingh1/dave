////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cropphenology.h
/// \brief Crop phenology including phu calculations					
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CROPPHENOLOGY_H
#define CROPPHENOLOGY_H

/// Calculation of down-scaling of lai during crop senescence
double senescence_curve(Pft& pft, double fphu);
void phu_init(cropphen_struct& ppftcrop, Gridcellpft& gridcellpft, Patch& patch);
//void crop_phenology(Pft& pft, Patch& patch);
void crop_phenology(Patch& patch);
/// Updates crop phen from yesterday's lai_daily
void leaf_phenology_crop(Pft& pft, Patch& patch);

#endif // CROPPHENOLOGY_H