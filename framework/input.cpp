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
