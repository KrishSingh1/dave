

#ifndef INPUTDEFINITIONS
#define INPUTDEFINITIONS
#include "guess.h"

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