////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file input.h
/// \brief Master class for all environmental input
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INPUT_H
#define INPUT_H

#include "guess.h"
#include "config.h"
#include "InData.h"
#include "inputdefinitions.h"
#include "landcoverinput.h"
#include "managementinput.h"
#include "lamarquendep.h"
#include "globalco2file.h"
#include "inputmodule.h"

/// Reads gridlist in lon-lat-description format from text unput file
void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist);
/// Parses the spatial resolution for a lon-lat Coord gridlist
double parse_gridlist_spatial_resolution(ListArray_id<Coord>& gridlist);
/// Sets simulation period accoring to instruction file settings and/or climate time period
void set_simulation_years(InputModule* input_module);

#endif // INPUT_H
