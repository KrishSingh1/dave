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

///////////////////////////////////////////////////////////////////////////////////////
// SOILPARAMETERS
// May be called from input/output module to initialise stand Soiltype objects when
// soil data supplied as LPJ soil code rather than soil physical parameter values


void soilparameters(Soiltype& soiltype, int soilcode) {

	// DESCRIPTION
	// Derivation of soil physical parameters given LPJ soil code

	// INPUT AND OUTPUT PARAMETER
	// soil = patch soil

	const double PERC_EXP = 2.0;
	// exponent in percolation equation [k2; LPJF]
	// (Eqn 31, Haxeltine & Prentice 1996)
	// Changed from 4 to 2 (Sitch, Thonicke, pers comm 26/11/01)


	if (soilcode<1 || soilcode>9)
		fail("soilparameters: invalid LPJ soil code (%d)", soilcode);


	if (textured_soil) {
		soiltype.sand_frac = data[soilcode - 1][7];
		soiltype.clay_frac = data[soilcode - 1][8];
	}
	else {
		// Using fixed values from Parton et al. (2010)
		soiltype.sand_frac = 0.28;
		soiltype.clay_frac = 0.12;
	}

	soiltype.silt_frac = 1 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = data[soilcode - 1][0];
	soiltype.perc_exp = PERC_EXP;
	soiltype.awc[0] = SOILDEPTH_UPPER * data[soilcode - 1][1];
	soiltype.awc[1] = SOILDEPTH_LOWER * data[soilcode - 1][1];
	soiltype.thermdiff_0 = data[soilcode - 1][2];
	soiltype.thermdiff_15 = data[soilcode - 1][3];
	soiltype.thermdiff_100 = data[soilcode - 1][4];
	soiltype.wp[0] = SOILDEPTH_UPPER * data[soilcode - 1][5];
	soiltype.wp[1] = SOILDEPTH_LOWER * data[soilcode - 1][5];
	soiltype.wsats[0] = SOILDEPTH_UPPER * data[soilcode - 1][6];
	soiltype.wsats[1] = SOILDEPTH_LOWER * data[soilcode - 1][6];
	soiltype.wtot = (data[soilcode - 1][1] + data[soilcode - 1][5]) * (SOILDEPTH_UPPER + SOILDEPTH_LOWER);

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}
