

#ifndef INPUTDEFINITIONS
#define INPUTDEFINITIONS
#include "guess.h"

/// Write land use fraction data to memory; enables efficient usage of randomised gridlists for parallell simulations
#define LUTOMEMORY

/// Default value for gridlist and text input spatial resolution.
/*  Input data will be parsed for finer resolution than the default value. For coarser resolutions, raise default value.or set manually */
const double DEFAULT_SPATIAL_RESOLUTION = 0.5;

/// Activates the conversion of coordinates to the centre of the gridcell before searching for the nearest data point.
/**	The gridlist spatial resolution is needed for finding the correct input coordinates with a finer-scaled gridlist.
 *	Make sure that the input module supports different spatial resolutions and that the gridlist is sufficiently large 
 *	to be parsed successfully or set gridlist_spatial_resolution manually in the code.
 */
const bool search_for_centre_of_gridcell = false;

namespace inputdef {

/// Type for storing grid cell longitude, latitude and description text
struct Coord {
	
	int id;
	double lon;
	double lat;
	xtring descrip;
};

}

#endif