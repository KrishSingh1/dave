/*
 * soilinput.cpp
 *
 *  Created on: 24 nov 2014
 *      Author: stefan
 */

#include <iterator>
#include "soilinput.h"

SoilInput::SoilInput(){};
SoilInput::~SoilInput(){};

namespace {
int find_index_soilfile_header(std::vector<std::string>& header, char* s) {
	int r = -1;
	if ( std::find(header.begin(), header.end(), s) != header.end() ) {
		r = std::find(header.begin(), header.end(), s)-header.begin()-2;
	}
	return r;
}
bool fuzz(std::string fname) {
	std::ifstream ifs(fname.c_str(), std::ifstream::in);
	if (!ifs.good()) {
		fail("load_soils: could not open %s for input", fname.c_str());
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


void SoilInput::init(std::string fname) {

	datatype = fuzz(fname);

	if (datatype) {
		load_lpj_soilcodes(fname);
	} else {
		loaddatafromfileMINERAL(fname);
	}
}

void SoilInput::init(std::string fname, std::vector<coord> gridlist) {
	for (std::vector<coord>::iterator it = gridlist.begin(); it != gridlist.end(); ++it) {
		coordinates.insert(*it);
	}
	init(fname);
}

bool SoilInput::load_lpj_soilcodes(std::string fname) {
	std::ifstream ifs(fname.c_str(), std::ifstream::in);

	std::string line;

	std::map<coord, int>::iterator it = soildatamapLPJ.begin();
	while (getline(ifs, line)) {
		std::istringstream iss(line);
		//SoilClass sc;
		// This routine searches the input file until it finds the location,
		// TODO not the most efficient implementation, but sufficient
		// since the data set is quite small
		// std::set<coord>::iterator IT;

		double lon_temp, lat_temp;
		int classnbr;
		if (iss >> lon_temp >> lat_temp >> classnbr) {
			coord c(lon_temp, lat_temp);
			if(coordinates.size() == 0 || coordinates.find(c) != coordinates.end()) {
				// SoilDataLPJ soildata;
				// soildata.soilcode = classnbr;
				// soildata.datatype = LPJSOIL;
				soildatamapLPJ.insert(it,std::pair<coord, int>(c,classnbr));
				it++;
				if(soildatamapLPJ.size() == coordinates.size()) {
					dprintf("TJUPP\n");
					break;
				}
			}
		}
	}
	ifs.close();
	//dprintf("MAP %d SET %d\n",soildatamapLPJ.size(), gridlist.size());
	return true;
}

bool SoilInput::loaddatafromfileMINERAL(std::string fname) {
	std::ifstream ifs(fname.c_str(), std::ifstream::in);

	std::string line;

	// reads the first line which should be a header
	getline(ifs, line);

	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> header(begin, end);

	int ncols = header.size();
	for (int i=0; i<ncols;i++) {
		std::transform(header[i].begin(), header[i].end(), header[i].begin(),::tolower);
	}
	int sand_i = find_index_soilfile_header(header, (char*)"sand");
	int clay_i = find_index_soilfile_header(header, (char*)"clay");
	int orgc_i = find_index_soilfile_header(header, (char*)"orgc");
	int ph_i   = find_index_soilfile_header(header, (char*)"ph");
	int bd_i   = find_index_soilfile_header(header, (char*)"bulkdensity");

	double lon_temp, lat_temp;
	std::map<coord,SoilDataMineral>::iterator it = soildatamapMINERAL.begin();
	while (getline(ifs, line)) {
		std::istringstream iss(line);
		//SoilClass sc;
		// This routine searches the input file until it finds the location,
		// TODO not the most efficient implementation, but sufficient
		// since the data set is quite small
		if (iss >> lon_temp >> lat_temp) {
			coord c(lon_temp, lat_temp);
			if(coordinates.size() == 0 || coordinates.find(c) != coordinates.end()) {
				std::vector<double> T(ncols-2);
				for (int i=0; i<ncols-2; i++) {
					iss  >> T[i];
				}
				SoilDataMineral soildata;
				soildata.sand = T[sand_i];
				soildata.clay = T[clay_i];
				soildata.orgc = T[orgc_i];
				soildata.pH = T[ph_i];
				if (bd_i<0) {
					soildata.bulkdensity = (double)bd_i;
				}
				soildatamapMINERAL.insert(it,std::pair<coord,SoilDataMineral>(c,soildata));
				it++;
				if(soildatamapLPJ.size() == coordinates.size()) {
					break;
				}
			}
		}
	}
	ifs.close();
	return true;
}



///////////////////////////////////////////////////////////////////////////////////////
// SOILPARAMETERS
// May be called from input/output module to initialise stand Soiltype objects when
// soil data supplied as LPJ soil code rather than soil physical parameter values


SoilInput::SoilProperties SoilInput::getsoilLPJ(coord c) {

	int soilcode = soildatamapLPJ[c];

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

	if (soilcode<1 || soilcode>9)
		fail("soilparameters: invalid LPJ soil code (%d)",soilcode);

	SoilProperties soiltype;
	soiltype.sand = data[soilcode-1][7];
	soiltype.clay = data[soilcode-1][8];

	soiltype.b = data[soilcode-1][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode-1][1];
	soiltype.thermal_wilting_point = data[soilcode-1][2];
	soiltype.thermal_15_whc = data[soilcode-1][3];
	soiltype.thermal_field_capacity = data[soilcode-1][4];
	soiltype.wilting_point = data[soilcode-1][5];
	soiltype.saturation_capacity = data[soilcode-1][6];
	soiltype.pH = 6.5;
	soiltype.soil_OC = 0.05;
	return soiltype;
}

SoilInput::SoilProperties SoilInput::getsoilMINERAL(coord c) {

	SoilDataMineral soil = soildatamapMINERAL[c];
	double sand = soil.sand;
	double clay = soil.clay;
	double silt = 1.0 - sand - clay;

	double b = 0.0;
	double logPsi_s = 0.0;
	double Theta_s = 0.0;
	double Theta_wilt = 0.0;
	double Theta_whc = 0.0;
	// Equation 1 from Cosby 1984
	// Psi = Psi_s * (Theta/Theta_s)^b
	// Psi is the pressure head in cm
	// *_s is the values at saturation
	// Theta is the volumetric moisture content in percent
	// Re-arranged to get the Theta
	// Theta = Theta_s * (Psi/Psi_s)^(1/b)


	// from Table 4, Cosby 1984
	b = 3.10+15.7 * clay - 0.3 * sand;
	//logK_s = -0.6 + 1.26 * soiltype.sand_frac - 0.64 * soiltype.clay_frac;
	logPsi_s = 1.54 - 0.95 * sand + 0.63 * silt;
	// Theta_s in Cosby expressed as %
	Theta_s = 0.01*(50.5 - 14.2 * sand - 3.7 * clay);

	double Psi_s = pow(10.0, -logPsi_s);
	double Psi_wilt = pow(10.0, -4.2);
	double Psi_whc = pow(10.0, -2.0);

	Theta_whc = Theta_s * pow((Psi_whc/Psi_s),1.0/b);
	Theta_wilt = Theta_s * pow((Psi_wilt/Psi_s),1.0/b);

	// A linear dependence between the percolation coefficient from Haxeltine 1996a
	// and the texture dependent parameter b from Cosby 1984 was established
	// K = 5.87 - 0.29*b

	SoilProperties soiltype;
	soiltype.b = 5.87 - 0.29 * b;

	soiltype.sand = sand;
	soiltype.clay = clay;

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

void SoilInput::soilparameters(SoilProperties soilprop, Soiltype& soiltype) {
	const double PERC_EXP = 2.0;
	soiltype.sand_frac = soilprop.sand;
	soiltype.clay_frac = soilprop.clay;
	soiltype.silt_frac = 1 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = soilprop.b;
	soiltype.perc_exp = PERC_EXP;
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

void SoilInput::getsoil(double lon, double lat, Soiltype& type) {
	coord c(lon, lat);
	SoilProperties props = datatype ? getsoilLPJ(c) : getsoilMINERAL(c);
	soilparameters(props, type);
}

