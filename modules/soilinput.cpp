/*
 * soilinput.cpp
 *
 *  Created on: 24 nov 2014
 *      Author: stefan
 */



#include "soilinput.h"

SoilInput::SoilInput(){};
SoilInput::~SoilInput(){};

int find_index_soilfile_header(std::vector<std::string>& header, char* s) {
	int r = -1;
	if ( std::find(header.begin(), header.end(), s) != header.end() ) {
		r = std::find(header.begin(), header.end(), s)-header.begin()-2;
	} else {
	   fail("No %s in the header of the soil file\n",s);
	}
	return r;
}

//void SoilInput::init(xtring fname, ListArray gridlist) {
//	// iterate over gridlist and create map
//	std::map<std::pair<double, double>, std::vector<double> >
//	return cache(fname, map)
//}
//
//void SoilInput::init(xtring fname, std::vector<Coord> gridlist) {
//	// iterate over vector  and
//	return cache(map);
//}
//bool cache(map) {


void SoilInput::init() {

	if (param["file_soilmap"].str != "") {

		soilmap_index = (int)param["soilmap_index"].num;
		if(param.find(xtring("use_ph_from_data"))) {
			use_ph_from_data = (bool)(int)param["use_ph_from_data"].num;
		} else {
			use_ph_from_data = false;
		}

		sand.push_back(0.0);
		silt.push_back(0.0);
		clay.push_back(0.0);
		orgc.push_back(0.0);
		ph.push_back(0.0);

		sand_i = silt_i = clay_i = orgc_i = ph_i = -1;
		xtring fname = param["file_soilmap"].str;
		std::ifstream ifs(fname, std::ifstream::in);
		bool found_soil = false;
		if (!ifs.good()) {
			fail("load_soils: could not open %s for input", (char*)fname);
		}

		std::string line;

		// reads the first line which should be a header
		getline(ifs, line);
		std::istringstream ss(line);
		std::istream_iterator<std::string> begin(ss);
		std::istream_iterator<std::string> end;
		std::vector<std::string> header(begin, end);
		ncols = header.size();
		for (int i=0; i<ncols;i++) {
			std::transform(header[i].begin(), header[i].end(), header[i].begin(),::tolower);
		}
		sand_i = find_index_soilfile_header(header, (char*)"sand");
		silt_i = find_index_soilfile_header(header, (char*)"silt");
		clay_i = find_index_soilfile_header(header, (char*)"clay");
		orgc_i = find_index_soilfile_header(header, (char*)"orgc");
		ph_i = 	 find_index_soilfile_header(header, (char*)"ph");

	}

}

bool SoilInput::loaddatafromfile(const double lon, const double lat) {
	if (param["file_soilmap"].str != "") {
		//dprintf("Soilmap coordinates lon %f lat %f\n",lon,lat);
		xtring fname = param["file_soilmap"].str;
		std::ifstream ifs(fname, std::ifstream::in);
		bool found_soil = false;
		if (!ifs.good()) {
			fail("load_soils: could not open %s for input", (char*)fname);
		}

		//std::vector<SoilClass> soil_properties;
		std::string line;

		// reads the first line which should be a header
		getline(ifs, line);

		double lon_temp, lat_temp;
		while (getline(ifs, line)) {
			std::istringstream iss(line);
			//SoilClass sc;
			// This routine searches the input file until it finds the location,
			// TODO not the most efficient implementation, but sufficient
			// since the data set is quite small
			if (iss >> lon_temp >> lat_temp) {
				if(fabs(lon_temp - lon) < 0.05 && fabs(lat_temp - lat) < 0.05){
					std::vector<double> T(ncols-2);
					for (int i=0; i<ncols-2; i++) {
						iss  >> T[i];
					}
					//iss  >> T[0] >> T[1] >> T[2] >> T[3] >> T[4] >> T[5];
					sand[0] = T[sand_i];
					silt[0] = T[silt_i];
					clay[0] = T[clay_i];
					orgc[0] = T[orgc_i];
					ph[0] = T[ph_i];
					found_soil = true;
					return found_soil;
				}
			}
		}
		ifs.close();
		fail("load_soils: could not find: lon %f, lat %f in file: %s",lon,lat, (char*)fname);
		return found_soil;
	} else {
		return false;
	}
}



bool SoilInput::getsoil(Gridcell &gridcell, const int soilstatus) {
	if (soilstatus>=0) {
		if (param["file_soilmap"].str != "") {
			int si = 0; // currently ont used
			if (si>=0 && soilmap_index>=0) {
				if (si>=soilmap_index) {
					si = soilmap_index;
				}
				gridcell.soiltype.sand_frac = sand[si];
				gridcell.soiltype.silt_frac = silt[si];
				gridcell.soiltype.clay_frac = clay[si];
				gridcell.soiltype.organic_frac = orgc[si];
				if (use_ph_from_data) {
					gridcell.soiltype.pH = ph[si];
				}
			}
			return true;
		} else {
			dprintf("Neither soil code, nor input soil data set is supplied\n");
			return false;
		}
	} else {
		return true; // standard values are used
	}
}
