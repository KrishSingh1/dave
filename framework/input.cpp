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

/// Setting of firsthistyear_sim and nyear_hist_sim; firsthistyear_sim and lasthistyear_sim are initialised to -1, nyear_hist_sim to 0.
void set_simulation_years(InputModule* input_module) {

	// Let climate input module set firsthistyear and lasthistyear
	if(!input_module->supports_firsthistyear_in_insfile()) {
		firsthistyear_sim = input_module->getfirsthistyear();
		nyear_hist_sim = (input_module->getfirsthistyear() + input_module->getnyear_hist() - 1) - firsthistyear_sim + 1;
	}
	// First look for firsthistyear and lasthistyear in instruction file
	if(firsthistyear_sim > -1 && lasthistyear_sim > -1) {
		nyear_hist_sim = lasthistyear_sim - firsthistyear_sim + 1;
	}
	else if(firsthistyear_sim < 0 || nyear_hist_sim == 0) {

		// If firsthistyear or nyear_hist not defined in instruction file, use values in climate input
		if(firsthistyear_sim < 0)
			firsthistyear_sim = input_module->getfirsthistyear_climate();
		if(nyear_hist_sim == 0) {
			if(input_module->getnyear_hist_climate() > 0)
				nyear_hist_sim = (input_module->getfirsthistyear_climate() + input_module->getnyear_hist_climate() - 1) - firsthistyear_sim + 1;
			else
				nyear_hist_sim = 0;
		}

		// If not already found, use values from landcover input module
		if(firsthistyear_sim < 0)
			firsthistyear_sim = input_module->get_landcover_module()->getfirsthistyear();
		if(nyear_hist_sim == 0) {
			if(input_module->get_landcover_module()->getnyear_hist() > 0)
				nyear_hist_sim = (input_module->get_landcover_module()->getfirsthistyear() + input_module->get_landcover_module()->getnyear_hist() - 1) - firsthistyear_sim + 1;
			else
				nyear_hist_sim = 0;
		}

		lasthistyear_sim = firsthistyear_sim + nyear_hist_sim - 1;	// Not used further

		// firsthistyear and nyear_hist must be defined by now
		if(firsthistyear_sim < 0 || nyear_hist_sim == 0)
			fail("firsthistyear_sim or nyear_hist_sim not defined\n");
	}
}

