///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H

#include "guess.h"

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell,Pftlist& pftlist);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell,Pftlist& pftlist);
void harvest_crop(double& cmass_plant,double& cmass_leaf,double& cmass_root,double& cmass_ho,double& cmass_agpool,double& litter_leaf,double& litter_root,
	double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
void harvest_natural(double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,double& cmass_debt,double& litter_leaf,double& litter_root,double& litter_wood,
	double& acflux_harvest,double& harvested_products_slow,Individual& indiv);
void harvest_pasture(double& cmass_leaf,double& cmass_root,double& litter_leaf,double& litter_root,double& acflux_harvest,double& harvested_products_slow,Individual& indiv);


#endif // LPJ_GUESS_LANDCOVER_H