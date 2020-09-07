///////////////////////////////////////////////////////////////////////////////////////
/// \file soilinput.cpp
/// \brief Implementation of soil input from text files, soil codes or soil physical properties.
///
///  Created on: 24 nov 2014
/// \author : Stefan Olin
///
/// $Date:  $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "soilinput.h"
#include <fstream>
#include <sstream>


namespace {
bool format_input_header(const char* fname) {
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

	size_t ncols = header.size();
	ifs.close();
	// Return TRUE if the soil data file has three columns (lon, lat and soilcode) and FALSE otherwise (such as when we specify sand, silt, clay, pH etc.) 
	return ncols == 3;
}

}
coord SoilInput::find_closest_point(double searchradius, coord C) {

	// First try the exact coordinate
	bool found = soil_code ? lpj_map.count(C) > 0 :
							 mineral_map.count(C) > 0;
	
	if (found) {
		return C;
	}

	double lon = C.first;
	double lat = C.second;

	if (searchradius == 0) {
		// Don't try to search
		fail("Coordinate at %f, %f could not be found in the soil map, and no search radius specified.\n",lon,lat);
		return C;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon * 2) / 2 + 0.25;
	double center_lat = floor(lat * 2) / 2 + 0.25;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy
	// sorting.
	using std::pair;
	using std::make_pair;
	std::vector<pair<double, coord> > search_points;

	
	const double EPS = 1e-15;

	for (double y = center_lon - searchradius; y <= center_lon + searchradius + EPS; y += STEP) {
		for (double x = center_lat - searchradius; x <= center_lat + searchradius + EPS; x += STEP) {
			double xdist = x - lat;
			double ydist = y - lon;
			double dist = sqrt(xdist * xdist + ydist * ydist);

			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found
	for (unsigned int i = 0; i < search_points.size(); i++) {
		coord search_point = search_points[i].second;

		found = soil_code ? lpj_map.count(search_point) > 0 :
							mineral_map.count(search_point) > 0;

		if (found) {
			return search_point;
		}
	}

	fail("Coordinate at %f, %f could not be found in the soil map.\n",lon,lat);
	return C;
}

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	soil_code = format_input_header(fname);

	if (soil_code) {
		load_lpj_soilcodes(fname, coords);
	}
	else {
		load_mineral_soils(fname, coords);
	}
}

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init_depth(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	while (getline(ifs, line)) {
		std::istringstream iss(line);

		double lon, lat;
		double depth;
		if (iss >> lon >> lat >> depth) {
			coord c(lon, lat);
			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}
			depth_map[c] = depth;
			if (depth_map.size() == coords.size()) {
				break;
			}
		}
	}
	ifs.close();
}

// Read in routine for text file with LPJ soil codes.
void SoilInput::load_lpj_soilcodes(const char* fname, const std::set<coord>& coords) {
	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	while (getline(ifs, line)) {
		std::istringstream iss(line);

		double lon, lat;
		int classnbr;
		if (iss >> lon >> lat >> classnbr) {
			coord c(lon, lat);
			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}
			if (classnbr<0 || classnbr>9) {
				fail("SoilInput::init: invalid LPJ soil code (%d) for location"
					" (%g, %g) in file %s", classnbr, lon, lat, fname);
			}
			lpj_map[c] = classnbr;
			if (lpj_map.size() == coords.size()) {
				break;
			}
		}
	}
	ifs.close();
}

// Load mineral soil data from text file with columns describing the data. They need to be named:
// lon, lat, sand, silt, clay, orgc, bulkdensity, ph, soilc, cn.
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
	std::vector<std::string>::iterator it = header.begin() + 2;

	// First four variables are intentionally left uninitialised, as those
	// columns are required.
	// int sand_i, clay_i, silt_i, orgc_i, ph_i, bd_i, cn_i, soilc_i = -1;
	int sand_i =-1;
	int clay_i = -1;
	int silt_i = -1; 
	int orgc_i = -1;
	int ph_i = -1;
	int bd_i = -1;
	int cn_i = -1; 
	int soilc_i = -1;

	for (int i = 0; it != header.end(); ++it, ++i) {
		if (*it == "sand") {
			sand_i = i;
		}
		else if (*it == "clay") {
			clay_i = i;
		}	
		else if (*it == "silt") {
			silt_i = i;
		}
		else if (*it == "orgc") {
			orgc_i = i;
		}
		else if (*it == "ph") {
			ph_i = i;
		}
		else if (*it == "bulkdensity") {
			bd_i = i;
		}	
		else if (*it == "cn") {
			cn_i = i;
		}
		else if (*it == "soilc") {
			soilc_i = i;
		}
	}
	if (soilc_i == -1 && iforganicsoilproperties == 1) {
		fail("Error! No Soil C column found in %s. \nTip: do not use iforganicsoilproperties 1 together with a soilmap file without a SoilC column\n", fname);
	}

	// Create a empty vector T with header.size-2 elements
	std::vector<double> T((unsigned int)header.size() - 2);

	while (getline(ifs, line)) {

		std::istringstream iss(line);
		double lon, lat;
		if (iss >> lon >> lat) {	// Get lon and lat from current line.
			coord c(lon, lat);

			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}

			for (std::vector<double>::iterator it = T.begin(); it != T.end(); ++it) {
				// Add the next values from the line to T
				iss >> *it;
			}

			SoilDataMineral& soildata = mineral_map[c];
			soildata.sand = T[sand_i];
			soildata.clay = T[clay_i];
			soildata.orgc = T[orgc_i];
			soildata.pH = T[ph_i];
			soildata.CN = T[cn_i];
			if (soilc_i != -1) {
				soildata.soilC = T[soilc_i];
			}
			else {
				soildata.soilC = 0.0;
			}
				
			// Not all data sets includes bulk density data, here it is set to a negative number if no column with that name. 
			// TODO, set it to value: 1.6
			if (bd_i<0) {	
				soildata.bulkdensity = (double)bd_i;
			} 
			else {
				soildata.bulkdensity = T[bd_i];
			}

			if (mineral_map.size() == coords.size()) {
				break;
			}
		}
	}

	ifs.close();
}

// Get and set soil properties based on LPJ soil code.
SoilInput::SoilProperties SoilInput::get_lpj(coord c) {

	coord C = find_closest_point(searchradius_soil, c);

	int soilcode = lpj_map[C];
	
	SoilProperties soiltype;

	// Mineral fractions
	soiltype.sand = data[soilcode][7];
	soiltype.clay = data[soilcode][8];

	// Percolation coefficient
	soiltype.b = data[soilcode][0];

	// Thermal diffusivity
	soiltype.thermal_wilting_point = data[soilcode][2];
	soiltype.thermal_15_whc = data[soilcode][3];
	soiltype.thermal_field_capacity = data[soilcode][4];

	// Water holding capacities
	soiltype.saturation_capacity = data[soilcode][6];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	soiltype.wilting_point = data[soilcode][5];

	// Soil pH
	soiltype.pH = 6.5;

	// Porosity fractions
	soiltype.soil_OC = 0.05;
	soiltype.porosity = data[soilcode][11];

	// Soil carbon
	soiltype.soilC = 0.0; // default value (kgC/m2)

	return soiltype;
}

// Get and set soil properties based on mineral soil input.
SoilInput::SoilProperties SoilInput::get_mineral(coord c) {

	coord C = find_closest_point(searchradius_soil, c);

	SoilDataMineral& soil = mineral_map[C];

	// Mineral silt fraction
	double silt = 1.0 - soil.sand - soil.clay;

	// Equation 1 from Cosby 1984
	// Psi = Psi_s * (Theta/Theta_s)^b
	// Psi is the pressure head in cm
	// *_s is the values at saturation
	// Theta is the volumetric moisture content in percent
	// Re-arranged to get the Theta
	// Theta = Theta_s * (Psi/Psi_s)^(1/b)

	// from Table 4, Cosby 1984
	double b = 3.10 + 15.7 * soil.clay - 0.3 * soil.sand;

	double logPsi_s = 1.54 - 0.95 * soil.sand + 0.63 * silt;

	// Theta_s in Cosby expressed as %
	double Theta_s = 0.01 * (50.5 - 14.2 * soil.sand - 3.7 * soil.clay);

	double Psi_s = pow(10.0, -logPsi_s);
	double Psi_wilt = pow(10.0, -4.2);
	double Psi_whc = pow(10.0, -2.0);

	double Theta_whc = Theta_s * pow((Psi_whc / Psi_s), 1.0 / b);
	double Theta_wilt = Theta_s * pow((Psi_wilt / Psi_s),1.0 / b);

	// A linear dependence between the percolation coefficient from Haxeltine 1996a
	// and the texture dependent parameter b from Cosby 1984 was established
	// K = 5.87 - 0.29*b

	SoilProperties soiltype;

	// Percolation coefficient
	soiltype.b = 5.87 - 0.29 * b;

	// Mineral fractions
	soiltype.sand = soil.sand;
	soiltype.clay = soil.clay;

	// Thermal diffusivity
	soiltype.thermal_wilting_point = 0.2;
	soiltype.thermal_15_whc = 0.15 * b + 0.05;
	soiltype.thermal_field_capacity = 0.4;

	// Water holding capacities
	soiltype.volumetric_whc_field_capacity = (Theta_whc - Theta_wilt);
	soiltype.saturation_capacity = Theta_s;
	soiltype.wilting_point = Theta_wilt;

	// Soil pH
	soiltype.pH = soil.pH;

	// Soil carbon content
	soiltype.soilC = soil.soilC;

	// Soil carbon fraction
	soiltype.soil_OC = soil.orgc;

	// Porosity fractions
	soiltype.porosity = Theta_s;
	
	return soiltype;
}

int get_texture(double clay, double sand, double silt) {

	int tc;

	if (clay > 0.4 && sand < 0.55 && silt < 0.4)
		tc = 0; // Clay
	else if(clay > 0.4 && sand < 0.2)
		tc = 1; // Silty Clay
	else if(clay > 0.27 && sand < 0.2)
		tc = 2; // Silty Clay Loam
	else if(clay > 0.27 && sand < 0.45)
		tc = 3; // Clay Loam
	else if(clay > 0.35)
		tc = 4; // Sandy Clay
	else if(clay < 0.13 && sand < 0.2 && silt > 0.8)
		tc = 5; // Silt
	else if(silt > 0.5)
		tc = 6; // Silt Loam
	else if(clay > 0.08 && sand < 0.52 && silt > 0.28)
		tc = 7; // Loam
	else if(clay > 0.2)
		tc = 8; // Sandy Clay Loam
	else if(clay < (0.1 - 0.67 * silt) && sand > 0.85 && silt < 0.15)
		tc = 9; // Sand
	else if(clay > (0.15 - 0.5 * silt) && sand > 0.7 && silt < 0.3)
		tc = 10; // Loamy Sand
	else
		tc = 11; // Sandy loam

	return tc;
}


void SoilInput::get_soil(double lon, double lat, Gridcell& gridcell) {

	Soiltype& soiltype = gridcell.soiltype;

	coord c(lon, lat);
	SoilProperties soilprop = soil_code ? get_lpj(c) : get_mineral(c);

	soiltype.soilcode = lpj_map[c];
	int organic_soilcode = 8;

	// Mineral fractions
	soiltype.sand_frac = soilprop.sand;
	soiltype.clay_frac = soilprop.clay;
	soiltype.silt_frac = 1.0 - soiltype.sand_frac - soiltype.clay_frac;

	// Porosity fractions for mineral and organic soil
	soiltype.organic_frac = soilprop.soil_OC;
	soiltype.porosity = soilprop.porosity;
	soiltype.mineral_frac = 1.0 - soiltype.organic_frac - soiltype.porosity;
	soiltype.porosity_org = data[organic_soilcode][11];

	// Percolation coefficients
	soiltype.perc_base = soilprop.b;
	soiltype.perc_base_evap = soilprop.b;
	soiltype.perc_exp = 2;
	soiltype.perc_base_org = data[organic_soilcode][0];

	// Thermal diffusivity for mineral and organic soil
	soiltype.thermdiff_0 = soilprop.thermal_wilting_point;
	soiltype.thermdiff_15 = soilprop.thermal_15_whc;
	soiltype.thermdiff_100 = soilprop.thermal_field_capacity;
	soiltype.thermdiff_0_org = data[organic_soilcode][2];
	soiltype.thermdiff_15_org = data[organic_soilcode][3];
	soiltype.thermdiff_100_org = data[organic_soilcode][4];

	// Water holding capacities for mineral and organic soil
	soiltype.wsats = soilprop.saturation_capacity;
	soiltype.awc = soilprop.volumetric_whc_field_capacity;
	soiltype.wp = soilprop.wilting_point;	
	soiltype.wsats_org = data[organic_soilcode][6];
	soiltype.awc_org = data[organic_soilcode][1];
	soiltype.wp_org = data[organic_soilcode][5];

	// Soil texture class based on clay, sand and silt fractions
	int texture = get_texture(soiltype.clay_frac, soiltype.sand_frac, soiltype.silt_frac);

	// van Genuchten soil parameters (Carsel 1988)
	soiltype.Kw_max = vG_data[texture][2];
	soiltype.alpha = vG_data[texture][3];
	soiltype.m = 1.0 - 1.0 / vG_data[texture][4];
	soiltype.l = 0.5; // Default value
	soiltype.Kw_max_org = vG_data[12][2];
	soiltype.alpha_org = vG_data[12][3];
	soiltype.m_org = 1.0 - 1.0 / vG_data[12][4];

	// Soil pH
	soiltype.pH = soilprop.pH;

	// Soil carbon content
	soiltype.soilC = soilprop.soilC;

	// Bulk density (g cm-3) (eqn 8.3 Bonan 2019 with particle density set to 2650 kg m–3
	soiltype.bulk_density = soilprop.bulk_density > 0.0 ? soilprop.bulk_density : (1.0 - soiltype.wsats) * 2.650;

	// Soil depth
	// soil_depth is defined in global.ins
	// if zero then use input map of soil depth
	if (negligible(soil_depth)) {
		soil_depth = depth_map[c] * MM_PER_M;
	}

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// SOILPARAMETERS
// May be called from input/output module to initialise stand Soiltype objects when
// soil data supplied as LPJ soil code rather than soil physical parameter values
void soil_parameters(Soiltype& soiltype, int soilcode) {

	// DESCRIPTION
	// Derivation of soil physical parameters given LPJ soil code

	// INPUT AND OUTPUT PARAMETER
	// soil = patch soil

	// exponent in percolation equation [k2; LPJF]
	// (Eqn 31, Haxeltine & Prentice 1996)
	// Changed from 4 to 2 (Sitch, Thonicke, pers comm 26/11/01)
	const double PERC_EXP = 2.0;

	if (soilcode<0 || soilcode>9) {
		fail("soil_parameters: invalid LPJ soil code (%d)", soilcode);
	}

	if (textured_soil) {
		soiltype.sand_frac = data[soilcode][7];
		soiltype.clay_frac = data[soilcode][8];
	}
	else {
		// Using fixed values from Parton et al. (2010)
		soiltype.sand_frac = 0.28;
		soiltype.clay_frac = 0.12;
	}

	soiltype.silt_frac = 1.0 - soiltype.sand_frac - soiltype.clay_frac;
	soiltype.perc_base = data[soilcode][0];
	soiltype.perc_base_evap = data[soilcode][0];
	soiltype.perc_exp = PERC_EXP;
	soiltype.thermdiff_0 = data[soilcode][2];
	soiltype.thermdiff_15 = data[soilcode][3];
	soiltype.thermdiff_100 = data[soilcode][4];

	soiltype.awc = data[soilcode][1];
	soiltype.wp = data[soilcode][5];
	soiltype.wsats =  data[soilcode][6];

	// Store the soilcode and new soil properties
	soiltype.soilcode = soilcode;
	soiltype.organic_frac = data[soilcode][9];
	soiltype.porosity = data[soilcode][11];
	soiltype.mineral_frac = 1.0 - soiltype.organic_frac - soiltype.porosity;

	if (!ifcentury) {
		// override the default SOM years with 70-80% of the spin-up period
		soiltype.updateSolveSOMvalues(nyear_spinup);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//  and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Carsel, R.F. and Parrish, R.S. 1988 Developing joint probability distributions of soil 
//   water retention characteristics. Water Resources Research, 24: 755-769.
