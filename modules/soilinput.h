/*
 * soilinput.h
 *
 *  Created on: 24 nov 2014
 *      Author: stefan
 */

#ifndef SOILINPUT_H_
#define SOILINPUT_H_

#include "guess.h"
#include <fstream>
#include <sstream>
#include <iterator>     // std::istream_iterator
#include <algorithm>
#include <vector>
class SoilInput {
public:
	SoilInput();

	~SoilInput();

	void init();

	bool getsoil(Gridcell &gridcell, const int soilcode);

	bool loaddatafromfile(const double lon, const double lat);

	bool loaddatafromfile(const int rlon, const int rlat);
private:
	// Soil variables
	// Vectors not needed in the current implementation
	// preparation for future implementations
	std::vector<double> sand;
	std::vector<double> silt;
	std::vector<double> clay;
	std::vector<double> orgc;
	std::vector<double> ph;
	std::vector<double> soil_area;

	// preparation for future implementations
	int soilmap_index;

	//
	bool use_ph_from_data;
	int sand_i, silt_i, clay_i, orgc_i, ph_i;
	int ncols;

};

#endif /* SOILINPUT_H_ */
