/*
 * soilinput.cpp
 *
 *  Created on: 24 nov 2014
 *      Author: stefan
 */

#include "soilinput.h"

namespace {
bool fuzz(const char* fname) {
	std::ifstream ifs(fname, std::ifstream::in);
	if (!ifs.good()) {
		fail("SoilInput::init: could not open %s for input", fname);
	}

	std::string line;

	// reads the first line which should be a header
	getline(ifs, line);

	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> header(begin, end);

	int ncols = header.size();
	ifs.close();
	return ncols == 3;
}
}


void SoilInput::init(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	lpj = fuzz(fname);

	if (lpj) {
		load_lpj_soilcodes(fname, coords);
	} else {
		load_mineral_soils(fname, coords);
	}
}

void SoilInput::load_lpj_soilcodes(const char* fname, const std::set<coord>& coords) {
	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	while (getline(ifs, line)) {
		std::istringstream iss(line);

		double lon, lat;
		int classnbr;
		if (iss >> lon >> lat >> classnbr) {
			coord c(lon, lat);
			if (!coords.empty() && coords.find(c) == coords.end()) {
				continue;
			}
			if (classnbr<1 || classnbr>9) {
				fail("SoilInput::init: invalid LPJ soil code (%d) for location"
					" (%g, %g) in file %s", classnbr, lon, lat, fname);
			}
			lpj_map[c] = classnbr - 1;
			if (lpj_map.size() == coords.size()) {
				break;
			}
		}
	}
	ifs.close();
}

void SoilInput::load_mineral_soils(const char* fname, const std::set<coord>& coords) {
	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	// reads the first line which should be a header
	getline(ifs, line);

	std::transform(line.begin(), line.end(), line.begin(), ::tolower);

	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> header(begin, end);
	std::vector<std::string>::iterator it = header.begin()+2;

	// first four variable are intentionally left uninitialised, as those
	// columns are requiered. Should blow up further down in while loop.
	int sand_i, clay_i, orgc_i, ph_i, bd_i = -1;

	for (int i=0; it != header.end(); ++it, ++i) {
		if (*it == "sand") {
			sand_i = i;
		} else if (*it == "clay") {
			clay_i = i;
		} else if (*it == "orgc") {
			orgc_i = i;
		} else if (*it == "ph") {
			ph_i = i;
		} else if (*it == "bulkdensity") {
			bd_i = i;
		}
	}

	std::vector<double> T(header.size()-2);

	while (getline(ifs, line)) {
		std::istringstream iss(line);
		double lon, lat;
		if (iss >> lon >> lat) {
			coord c(lon, lat);
			if (!coords.empty() && coords.find(c) == coords.end()) {
				continue;
			}
			for (std::vector<double>::iterator it=T.begin(); it != T.end(); ++it) {
				iss  >> *it;
			}
			SoilDataMineral& soildata = mineral_map[c];
			soildata.sand = T[sand_i];
			soildata.clay = T[clay_i];
			soildata.orgc = T[orgc_i];
			soildata.pH = T[ph_i];
			if (bd_i<0) {
				soildata.bulkdensity = (double)bd_i;
			}

			if (mineral_map.size() == coords.size()) {
				break;
			}
		}
	}
	ifs.close();
}

SoilInput::SoilProperties SoilInput::get_lpj(coord c) {

	double data[9][9] = {

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

		{   5.0, 0.110,   0.2, 0.800,   0.4,	0.074,	0.395,	0.90,	0.05},    // 1	Coarse
		{   4.0, 0.150,   0.2, 0.650,   0.4,	0.184,	0.439,	0.35,	0.15},    // 2	Medium
		{   3.0, 0.120,   0.2, 0.500,   0.4,	0.274,	0.454,	0.30,	0.45},    // 3	Fine
		{   4.5, 0.130,   0.2, 0.725,   0.4,	0.129,	0.417,	0.60,	0.15},    // 4	Medium-coarse
		{   4.0, 0.115,   0.2, 0.650,   0.4,	0.174,	0.425,	0.60,	0.30},    // 5	Fine-coarse
		{   3.5, 0.135,   0.2, 0.575,   0.4,	0.229,	0.447,	0.20,	0.30},    // 6	Fine-medium
		{   4.0, 0.127,   0.2, 0.650,   0.4,	0.177,	0.430,	0.45,	0.25},    // 7	Fine-medium-coarse
		{   9.0, 0.300,   0.1, 0.100,   0.1,	0.200,	0.600,	0.28,	0.12},    // 8	Organic (values not know for wp), sand and clay are from Parton 2010
		{   0.2, 0.100,   0.2, 0.500,   0.4,	0.100,	0.250,	0.10,	0.80}     // 9	Vertisols (values not know for wp)
	};

	int soilcode = lpj_map[c];

	SoilProperties soiltype;
	soiltype.sand = data[soilcode][7];
	soiltype.clay = data[soilcode][8];

	soiltype.b = data[soilcode][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	soiltype.thermal_wilting_point = data[soilcode][2];
	soiltype.thermal_15_whc = data[soilcode][3];
	soiltype.thermal_field_capacity = data[soilcode][4];
	soiltype.wilting_point = data[soilcode][5];
	soiltype.saturation_capacity = data[soilcode][6];
	soiltype.pH = 6.5;
	soiltype.soil_OC = 0.05;
	return soiltype;
}

SoilInput::SoilProperties SoilInput::get_mineral(coord c) {

	SoilDataMineral& soil = mineral_map[c];
	double silt = 1.0 - soil.sand - soil.clay;

	// Equation 1 from Cosby 1984
	// Psi = Psi_s * (Theta/Theta_s)^b
	// Psi is the pressure head in cm
	// *_s is the values at saturation
	// Theta is the volumetric moisture content in percent
	// Re-arranged to get the Theta
	// Theta = Theta_s * (Psi/Psi_s)^(1/b)

	// from Table 4, Cosby 1984
	double b = 3.10 + 15.7*soil.clay - 0.3*soil.sand;

	//logK_s = -0.6 + 1.26 * soiltype.sand_frac - 0.64 * soiltype.clay_frac;
	double logPsi_s = 1.54 - 0.95 * soil.sand + 0.63 * silt;

	// Theta_s in Cosby expressed as %
	double Theta_s = 0.01*(50.5 - 14.2 * soil.sand - 3.7 * soil.clay);

	double Psi_s = pow(10.0, -logPsi_s);
	double Psi_wilt = pow(10.0, -4.2);
	double Psi_whc = pow(10.0, -2.0);

	double Theta_whc = Theta_s * pow((Psi_whc/Psi_s),1.0/b);
	double Theta_wilt = Theta_s * pow((Psi_wilt/Psi_s),1.0/b);

	// A linear dependence between the percolation coefficient from Haxeltine 1996a
	// and the texture dependent parameter b from Cosby 1984 was established
	// K = 5.87 - 0.29*b

	SoilProperties soiltype;
	soiltype.b = 5.87 - 0.29 * b;

	soiltype.sand = soil.sand;
	soiltype.clay = soil.clay;

	soiltype.volumetric_whc_field_capacity = (Theta_whc - Theta_wilt);
	soiltype.thermal_wilting_point = 0.2;
	soiltype.thermal_15_whc = 0.15 * b + 0.05;
	soiltype.thermal_field_capacity = 0.4;
	soiltype.wilting_point = Theta_wilt;
	soiltype.saturation_capacity = Theta_s;
	soiltype.pH = soil.pH;
	soiltype.soil_OC = soil.orgc;
	return soiltype;
}

void SoilInput::getsoil(double lon, double lat, Soiltype& soiltype) {
	coord c(lon, lat);
	SoilProperties soilprop = lpj ? get_lpj(c) : get_mineral(c);

	soiltype.sand_frac = soilprop.sand;
	soiltype.clay_frac = soilprop.clay;
	soiltype.silt_frac = 1 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = soilprop.b;
	soiltype.perc_exp = 2;
	soiltype.awc[0] = SOILDEPTH_UPPER * soilprop.volumetric_whc_field_capacity;
	soiltype.awc[1] = SOILDEPTH_LOWER * soilprop.volumetric_whc_field_capacity;
	soiltype.thermdiff_0 = soilprop.thermal_wilting_point;
	soiltype.thermdiff_15 = soilprop.thermal_15_whc;
	soiltype.thermdiff_100 = soilprop.thermal_field_capacity;
	soiltype.wp[0] = SOILDEPTH_UPPER * soilprop.wilting_point;
	soiltype.wp[1] = SOILDEPTH_LOWER * soilprop.wilting_point;
	soiltype.wsats[0] = SOILDEPTH_UPPER * soilprop.saturation_capacity;
	soiltype.wsats[1] = SOILDEPTH_LOWER *  soilprop.saturation_capacity;
	soiltype.wtot = (soilprop.volumetric_whc_field_capacity + soilprop.wilting_point) * (SOILDEPTH_UPPER + SOILDEPTH_LOWER);
	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}

