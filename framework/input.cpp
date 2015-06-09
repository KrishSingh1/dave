////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file input.cpp
/// \brief Master class for all environmental input		
/// \author Mats Lindeskog
/// $Date: $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "input.h"
#include "inputdefinitions.h"

void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist) {

	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon, dlat;
	bool eof = false;
	xtring descrip;

	// Read list of grid coordinates and store in Coord object 'gridlist'

	FILE* in_grid = fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input", (char*)file_gridlist);

	while (!eof) {
		
		// Read next record in file
		eof =! readfor(in_grid, "f,f,a#", &dlon, &dlat, &descrip);

		if (!eof && !(dlon == 0.0 && dlat == 0.0)) { // ignore blank lines at end (if any)
			Coord& c = gridlist.createobj(); // add new coordinate to grid list

			c.lon = dlon;
			c.lat = dlat;
			c.descrip = descrip;
		}
	}
	fclose(in_grid);
	gridlist.firstobj();
}

double parse_gridlist_spatial_resolution(ListArray_id<Coord>& gridlist) {

	// Parse spatial resolution
	double precision = 100;
	double dif_lon;
	double dif_lat;
	const int maxnsample = 200;
	int nsample = min((unsigned)maxnsample, gridlist.nobj);

	for (int i=0; i<nsample; i++) {
		for (int j=0; j<nsample; j++) {

			dif_lon = fabs(gridlist[i].lon - gridlist[j].lon);
			dif_lat = fabs(gridlist[i].lat - gridlist[j].lat);
			if(dif_lon > 1.0e-12)
				precision = min(precision, dif_lon);
			if(dif_lat > 1.0e-12)
				precision = min(precision, dif_lat);
		}
	}
	return precision;
}

