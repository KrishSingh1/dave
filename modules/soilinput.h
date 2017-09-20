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

const double data[9][9] = {

	//    0  empirical parameter in percolation equation (k1) (mm/day)
	//    1  volumetric water holding capacity at field capacity minus vol water
	//       holding capacity at wilting point (Hmax), as fraction of soil layer
	//       depth
	//    2  thermal diffusivity (mm2/s) at wilting point (0% WHC)
	//    3  thermal diffusivity (mm2/s) at 15% WHC
	//    4  thermal diffusivity at field capacity (100% WHC)
	//       Thermal diffusivities follow van Duin (1963),
	//       Jury et al (1991), Fig 5.11.
	//    5  wilting point as fraction of depth (calculation method described in
	//       Prentice et al 1992)
	//    6  saturation capacity following Cosby (1984)
	//    7  sand fraction
	//    8  clay fraction

	//    0      1      2      3      4      5      6       7       8         soilcode
	//  ------------------------------------------------------------------

	{ 5.0, 0.110,   0.2, 0.800,   0.4,	0.074,	0.395,	0.90,	0.05 },    // 1	Coarse
	{ 4.0, 0.150,   0.2, 0.650,   0.4,	0.184,	0.439,	0.35,	0.15 },    // 2	Medium
	{ 3.0, 0.120,   0.2, 0.500,   0.4,	0.274,	0.454,	0.30,	0.45 },    // 3	Fine
	{ 4.5, 0.130,   0.2, 0.725,   0.4,	0.129,	0.417,	0.60,	0.15 },    // 4	Medium-coarse
	{ 4.0, 0.115,   0.2, 0.650,   0.4,	0.174,	0.425,	0.60,	0.30 },    // 5	Fine-coarse
	{ 3.5, 0.135,   0.2, 0.575,   0.4,	0.229,	0.447,	0.20,	0.30 },    // 6	Fine-medium
	{ 4.0, 0.127,   0.2, 0.650,   0.4,	0.177,	0.430,	0.45,	0.25 },    // 7	Fine-medium-coarse
	{ 9.0, 0.300,   0.1, 0.100,   0.1,	0.200,	0.600,	0.28,	0.12 },    // 8	Organic (values not know for wp), sand and clay are from Parton 2010
	{ 0.2, 0.100,   0.2, 0.500,   0.4,	0.100,	0.250,	0.10,	0.80 }     // 9	Vertisols (values not know for wp)
};

/// Set soil physical properties based in soilcode
void soilparameters(Soiltype& soiltype, int soilcode);


#endif /* SOILINPUT_H */
