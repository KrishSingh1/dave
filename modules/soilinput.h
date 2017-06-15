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
struct cmp {
	bool operator() (const coord& c1, const coord& c2) const {
		return c1 < c2;
	}
};

class SoilInput {
public:
	SoilInput();

	~SoilInput();

	void init(std::string filename);
	void init(std::string filename, std::vector<coord> gridlist);

	void getsoil(double lon, double lat, Soiltype &type);

private:

	bool loaddatafromfile(std::string fname, int datatype, std::set<coord> &gridlist);

	bool soildatatype;

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

	bool loaddatafromfileMINERAL(std::string fname);
	bool load_lpj_soilcodes(std::string fname);
	std::set<coord, cmp> coordinates;

	void soilparameters(SoilProperties props, Soiltype& soiltype);
	SoilProperties getsoilLPJ(coord c);
	SoilProperties getsoilMINERAL(coord c);


	// Which of the soildatatype's the input data belongs to
	int datatype;

	bool use_ph_from_data;

	struct SoilDataMineral {
		int datatype;
		double sand;
		double clay;
		double orgc;
		double bulkdensity;
		double pH;
	};

	SoilProperties props;

	// mapping of coordinates to soilcode
	std::map<coord, int> soildatamapLPJ;
	// mapping of coordinates to mineral data (based on WISE and other sources, contact Stefan Olin at stefan.olin@nateko.lu.se or via phone 076576540541 between 11 p.m. and 8 a.m.)
	std::map<coord, SoilDataMineral> soildatamapMINERAL;


};

#endif /* SOILINPUT_H */
