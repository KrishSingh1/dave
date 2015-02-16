

#ifndef INPUTDEFINITIONS
#define INPUTDEFINITIONS

#define LUTOMEMORY		// Write land use fraction data to memory; enables efficient usage of randomized gridlists for parallell runs on Simba.

//namespace inputdef {

/// Type for storing grid cell longitude, latitude and description text
struct Coord {
	
	int id;
	double lon;
	double lat;
	xtring descrip;
};
//}
#endif