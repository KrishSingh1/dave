/*
 * soilinput.h
 *
 *  Created on: 24 nov 2014
 *      Author: stefan
 */

#ifndef SOILINPUT_H
#define SOILINPUT_H

#include "guess.h"
#include <fstream>
#include <sstream>
#include <iterator>     // std::istream_iterator
#include <algorithm>
#include <vector>
#include <map>
#include <set>

typedef std::pair<double, double> coord;

class SoilInput {
public:
	SoilInput() {};

	~SoilInput() {};

	void init(const char* filename, const std::vector<coord>& gridlist=std::vector<coord>());

	void getsoil(double lon, double lat, Soiltype &type);

private:

	bool lpj;

	struct SoilProperties {
		//    0  empirical parameter in percolation equation (k1) (mm/day)
		double b;
		//    1  volumetric water holding capacity at field capacity minus vol water
		double volumetric_whc_field_capacity;
		//       holding capacity at wilting point (Hmax), as fraction of soil layer
		//       depth
		//    2  thermal diffusivity (mm2/s) at wilting point (0% WHC)
		double thermal_wilting_point;
		//    3  thermal diffusivity (mm2/s) at 15% WHC
		double thermal_15_whc;
		//    4  thermal diffusivity at field capacity (100% WHC)
		//       Thermal diffusivities follow van Duin (1963),
		//       Jury et al (1991), Fig 5.11.
		double thermal_field_capacity;
		//    5  wilting point as fraction of depth (calculation method described in
		//       Prentice et al 1992)
		double wilting_point;
		//    6  saturation capacity following Cosby (1984)
		double saturation_capacity;
		//    7  sand fraction
		double sand;
		//    8  clay fraction
		double clay;
		double bulk_density;
		double pH;
		double soil_OC;
	};

	void load_mineral_soils(const char* fname, const std::set<coord>& coords);
	void load_lpj_soilcodes(const char* fname, const std::set<coord>& coords);

	SoilProperties get_lpj(coord c);
	SoilProperties get_mineral(coord c);

	struct SoilDataMineral {
		double sand;
		double clay;
		double orgc;
		double bulkdensity;
		double pH;
	};

	std::map<coord, int> lpj_map;
	std::map<coord, SoilDataMineral> mineral_map;
};

#endif /* SOILINPUT_H */
