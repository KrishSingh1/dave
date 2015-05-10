

#ifndef INPUTDEFINITIONS
#define INPUTDEFINITIONS
#include "guess.h"

/// Write land use fraction data to memory; enables efficient usage of randomised gridlists for parallell simulations
#define LUTOMEMORY

/// Default value for gridlist and text input spatial resolution.
/*  Input data will be parsed for finer resolution than the default value. For coarser resolutions, raise default value.or set manually */
const double DEFAULT_SPATIAL_RESOLUTION = 0.5;

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