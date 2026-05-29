///////////////////////////////////////////////////////////////////////////////////////
/// \file soilinput.cpp
/// \brief Implementation of soil input from text files, soil codes or soil physical properties.
///
///  Created on: 24 nov 2014
/// \author : Stefan Olin
///
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

const double PERC_EXP_DAVE = 4.0;

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

	// Try a closest-match.
	const double EXACT_MATCH_EPS = 1e-3;
	double dist_min = 1e100;
	if (!soil_code) {
		// Search through the map
		for (std::map<coord, SoilDataMineral>::iterator it = mineral_map.begin(); it != mineral_map.end(); ++it) {
			coord search_point = it->first;
			double xdist = search_point.second - lat;
			double ydist = search_point.first - lon;
			double dist = sqrt(xdist * xdist + ydist * ydist);

			dist_min = min(dist, dist_min);
			if (dist <= EXACT_MATCH_EPS) {
				return search_point;
			}
		}
	}

	fail("Coordinate at %f, %f could not be found in the soil map.\n",lon,lat);
	return C;
}

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	xtring soilClay = param["file_soilclay"].str;
	xtring soilSand = param["file_soilsand"].str;

	// Have soil depth varying properties
	if (soilClay.len() && soilSand.len()) {

		soil_varying = true;

		init_depth_varying(coords);
	} 
	else {

		soil_varying = false;

		soil_code = format_input_header(fname);

		if (soil_code) {
			load_lpj_soilcodes(fname, coords);
		}
		else {
			load_mineral_soils(fname, coords);
		}
	}
}

// Helper: extract depths from header
std::vector<double> SoilInput::extract_layer_depths(std::ifstream& prop) {

	std::string line;
	getline(prop, line);
	std::transform(line.begin(), line.end(), line.begin(), ::tolower);

	const int nr_coord = 2; // lon + lat
	std::istringstream ss(line);
	std::istream_iterator<std::string> begin(ss), end;
	std::vector<std::string> header(begin, end);
	const int nr_layers = static_cast<int>(header.size()) - nr_coord;

	std::vector<double> lyr_depth(nr_layers);
	for (int h = nr_coord; h < static_cast<int>(header.size()); ++h) {
		size_t dashPos = header[h].find('-');
		std::string lastNumberStr = header[h].substr(dashPos + 1);
		lyr_depth[h - nr_coord] = std::stod(lastNumberStr);
	}

	return lyr_depth;
}

void SoilInput::get_depth_varying_prop(
	const std::set<coord>& coords,
	std::ifstream& prop,
	std::map<coord, SoilDataMineral>& mineral_map,
	std::vector<double> SoilDataMineral::* member) {
	if (!prop.is_open())
		return;

	// Read and parse header to extract layer depths
	std::vector<double> lyr_depth = extract_layer_depths(prop);
	const int nr_layers = static_cast<int>(lyr_depth.size());

	std::vector<double> T(nr_layers);
	std::string line;

	int sites = 0;

	while (getline(prop, line)) {
		std::istringstream iss(line);
		double lon, lat;
		if (!(iss >> lon >> lat))
			continue;

		coord c(lon, lat);
		if (!coords.empty() && coords.count(c) == 0)
			continue;

		for (double& val : T)
			iss >> val;

		SoilDataMineral& soildata = mineral_map[c];
		(soildata.*member).resize(nr_layers);
		for (int lyr = 0; lyr < nr_layers; ++lyr)
			(soildata.*member)[lyr] = T[lyr];

		sites++;

		if (sites == (int)coords.size())
			break;
	}
}

void SoilInput::assign_depths(
	const std::set<coord>& coords,
	std::map<coord, SoilDataMineral>& mineral_map,
	const std::vector<double>& depths)
{
	for (const auto& c : coords) {
		SoilDataMineral& soildata = mineral_map[c];
		soildata.depth_lyr = depths;
	}
}

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init_depth_varying(const std::set<coord>& coords) {

	soil_varying = true;

	// Assign depths once
	std::ifstream clayHeader(param["file_soilclay"].str);
	std::vector<double> depths = extract_layer_depths(clayHeader);
	assign_depths(coords, mineral_map, depths);

	std::ifstream clayFile(param["file_soilclay"].str);
	get_depth_varying_prop(coords, clayFile, mineral_map, &SoilDataMineral::clay_lyr);

	std::ifstream sandFile(param["file_soilsand"].str);
	get_depth_varying_prop(coords, sandFile, mineral_map, &SoilDataMineral::sand_lyr);

	std::ifstream coarseFile(param["file_soilcoarse"].str);
	if (coarseFile) {
		get_depth_varying_prop(coords, coarseFile, mineral_map, &SoilDataMineral::coarse_lyr);
	}
	else {
		dprintf("Using soil depth varying properties, but not coarse fragments!\n");
	}

	std::ifstream pHFile(param["file_soilpH"].str);
	if (pHFile) {
		get_depth_varying_prop(coords, pHFile, mineral_map, &SoilDataMineral::pH_lyr);
	}
	else {
		dprintf("Using soil depth varying properties, but not pH!\n");
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

// Initialises the soil input class and determines the format of the supplied input file.
void SoilInput::init_slope(const char* fname, const std::vector<coord>& gridlist) {

	std::set<coord> coords(gridlist.begin(), gridlist.end());

	std::ifstream ifs(fname, std::ifstream::in);

	std::string line;

	while (getline(ifs, line)) {
		std::istringstream iss(line);

		double lon, lat;
		double slope_mean, slope_median, slope_max, slope_min, slope_std;
		if (iss >> lon >> lat >> slope_mean >> slope_median >> slope_max >> slope_min >> slope_std) {
			coord c(lon, lat);
			if (!coords.empty() && coords.count(c) == 0) {
				continue;
			}
			slope_map[c] = slope_median;
			if (slope_map.size() == coords.size()) {
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

	soiltype.K = data[soilcode][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	// Percolation coefficient
	soiltype.b = data[soilcode][0];

	// Thermal diffusivity
	soiltype.thermal_wilting_point = data[soilcode][2];
	soiltype.thermal_15_whc = data[soilcode][3];
	soiltype.thermal_field_capacity = data[soilcode][4];

	// Water holding capacities
	soiltype.saturation_capacity = data[soilcode][6];
	soiltype.pH = 6.5;
	soiltype.soil_OC = 0.05;
	soiltype.soilC = 0.0; // default value (kgC/m2)
	soiltype.porosity = data[soilcode][11];
	soiltype.Psi_s = 20; // cm, just a guess
	soiltype.b = 4; // unitless, just a guess
	soiltype.Psi_wp = 5; // cm, just a guess
	soiltype.Psi_whc = 15; // cm, just a guess
	return soiltype;
}

// Stores the basic properties of organic soil. 
// Used only if scaling between pure mineral and organic soils is activated.
SoilInput::SoilProperties SoilInput::get_lpj_organic_soil() {

	int soilcode = 8;

	SoilProperties soiltype;
	soiltype.sand = data[soilcode][7];
	soiltype.clay = data[soilcode][8];

	soiltype.K = data[soilcode][0];
	soiltype.volumetric_whc_field_capacity = data[soilcode][1];
	soiltype.wilting_point = data[soilcode][5];

	// Soil pH
	soiltype.pH = 6.5;

	// Porosity fractions
	soiltype.soil_OC = 0.05;

	// Soil carbon
	soiltype.soilC = 0.0; // default value (kgC/m2)

	return soiltype;
}

// Get and set soil properties based on mineral soil input.
SoilInput::SoilProperties SoilInput::get_mineral(coord c) {

	coord C = find_closest_point(searchradius_soil, c);

	SoilDataMineral& soil = mineral_map[C];

	// soil.sand = 0.7; ! Testing sandier soil at Calperum, reduced soil evap but didn't change LAI much
	// soil.clay = 0.15;

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

	double Psi_s = 0.0;
	double Psi_wilt = 0.0;
	double Psi_whc = 0.0;
	double Theta_wilt = 0.0;
	double Theta_whc = 0.0;

	if (ifdavehydrology){

		Psi_s = -pow(10.0, logPsi_s);  /// CS changed from pow(10.0, -logPsi_s), why the negative there??? Compensates for not having '-b' in equation below but these Psi values didn't make sense (and now I'm using them elsewhere)
		Psi_wilt = -pow(10.0, 5.3);   /// CS changed from pow(10.0, -4.2). See Hohenbrink et al (2023), output is in cm. Testing parameter change from 4.2 to 5.3, to allow plants/evap to extract more water, up to ~ -20MPa.
		Psi_whc = -pow(10.0, 1.8);    /// CS changed from pow(10.0, -2). See Hohenbrink et al (2023), output is in cm, 1.8

		Theta_wilt = Theta_s * pow((Psi_wilt / Psi_s), 1.0 / -b);  // CS changed to -b, I think there's a mistake in Cosby 1984
		Theta_whc = Theta_s * pow((Psi_whc / Psi_s), 1.0 / -b); 

	} else {

		Psi_s = pow(10.0, -logPsi_s);  
		Psi_wilt = pow(10.0, -4.2);   
		Psi_whc = pow(10.0, -2.0);    

		Theta_wilt = Theta_s * pow((Psi_wilt / Psi_s), 1.0 / b);  
		Theta_whc = Theta_s * pow((Psi_whc / Psi_s), 1.0 / b); 

	}

	// A linear dependence between the percolation coefficient from Haxeltine 1996a
	// and the texture dependent parameter b from Cosby 1984 was established
	// K = 5.87 - 0.29*b

	SoilProperties soiltype;
	soiltype.K = 5.87 - 0.29 * b; // note can't find where this comes from or what the units are???

	// hydraulic conductivity (m/s) based on equation 10 from Saxton et al (1986)
	soiltype.Ksat = 2.778e-6 * exp(12.012 - 0.0755 * soil.sand*100.0 + (-3.895 + 0.03671 * soil.sand*100.0 - 0.1103 * 
					soil.clay*100.0 + 8.7546e-4 * soil.clay*100.0) * 1.0/Theta_s);

	// Mineral fractions
	soiltype.sand = soil.sand;
	soiltype.clay = soil.clay;

	// Soil pH
	soiltype.pH = soil.pH;

	// Soil carbon content
	soiltype.soilC = soil.soilC;
	soiltype.Psi_s = Psi_s;
	soiltype.b = b;
	soiltype.Theta_s = Theta_s;
	soiltype.Psi_wp = Psi_wilt;
	soiltype.Psi_whc = Psi_whc;

	// Soil carbon fraction
	soiltype.soil_OC = soil.orgc;
	
	return soiltype;
}


void SoilInput::get_soil(double lon, double lat, Gridcell& gridcell) {

	Soiltype& soiltype = gridcell.soiltype;

	coord c(lon, lat);

	// Organic properties
	int organic_soilcode = 8;

	// Percolation coefficients
	soiltype.perc_exp = ifdavehydrology ? PERC_EXP_DAVE : 2;
	soiltype.perc_base_org = data[organic_soilcode][0];

	// Thermal diffusivity for organic soil
	soiltype.thermdiff_0_org = data[organic_soilcode][2];
	soiltype.thermdiff_15_org = data[organic_soilcode][3];
	soiltype.thermdiff_100_org = data[organic_soilcode][4];

	if (!soil_varying) {
		SoilProperties soilprop = soil_code ? get_lpj(c) : get_mineral(c);

		// Homogeneous soil profile when it comes to minerals and pH
		soiltype.soil_input = soil_code ? LPJ : MINERAL;

		int nr_lyr = 1;

		soiltype.sand_frac.resize(nr_lyr);
		soiltype.clay_frac.resize(nr_lyr);
		soiltype.coarse_frac.resize(nr_lyr);
		soiltype.organic_frac.resize(nr_lyr);
		soiltype.soilC.resize(nr_lyr);
		soiltype.pH.resize(nr_lyr);
		soiltype.text_depths.resize(nr_lyr);

		// Mineral fractions
		soiltype.sand_frac[0] = soilprop.sand;
		soiltype.clay_frac[0] = soilprop.clay;

		// Organic fraction
		soiltype.organic_frac[0] = soilprop.soil_OC;

		if (soiltype.soil_input == LPJ) {

			// Percolation coefficients
			soiltype.perc_base = soilprop.K;			

			// Thermal diffusivity for mineral soil
			soiltype.thermdiff_0 = soilprop.thermal_wilting_point;
			soiltype.thermdiff_15 = soilprop.thermal_15_whc;
			soiltype.thermdiff_100 = soilprop.thermal_field_capacity;

			// Water holding capacities for mineral and organic soil
			soiltype.wsats_min = soilprop.saturation_capacity;
			soiltype.awc_min = soilprop.volumetric_whc_field_capacity;
			soiltype.wp_min = soilprop.wilting_point;
		}

		// Soil pH
		soiltype.pH[0] = soilprop.pH;

		// Soil carbon content
		soiltype.soilC[0] = soilprop.soilC;

		// Texture depth set to 100 mm
		soiltype.text_depths[0] = 100.0;
	}
	else {
		// Inhomogeneous soil profile when it comes to minerals and pH
		soiltype.soil_input = VARYING;

		SoilDataMineral& soil = mineral_map[c];

		int nr_lyr = (int)soil.depth_lyr.size();

		soiltype.sand_frac.resize(nr_lyr);
		soiltype.clay_frac.resize(nr_lyr);
		soiltype.coarse_frac.resize(nr_lyr);
		soiltype.organic_frac.resize(nr_lyr);
		soiltype.soilC.resize(nr_lyr);
		soiltype.pH.resize(nr_lyr);
		soiltype.text_depths.resize(nr_lyr);

		for (int lyr = 0; lyr < nr_lyr; lyr++) {

			// Mineral fractions
			soiltype.sand_frac[lyr] = soil.sand_lyr[lyr] * PERCENT_TO_FRACT;
			soiltype.clay_frac[lyr] = soil.clay_lyr[lyr] * PERCENT_TO_FRACT;

			// Coarse fragmnets
			soiltype.coarse_frac[lyr] = soil.coarse_lyr[lyr] * PERCENT_TO_FRACT;

			// pH
			if (!negligible(soil.pH_lyr[lyr])) {
				soiltype.pH[lyr] = soil.pH_lyr[lyr];
			}
			else {
				soiltype.pH[lyr] = 6.5;
			}

			// Depths to bottom (mm)
			soiltype.text_depths[lyr] = soil.depth_lyr[lyr] * MM_PER_CM;
		}
	}

	// Slope of soil
	if (slope < 0.0) {
		if (slope_map.size()) {
			slope = slope_map[c];
		}
		else {
			// No input data for slope. Set it to 2 degree as default
			slope = 2.0;
		}
	}

	// Slope [degree] limits
	slope = min(90.0, max(0.0, slope));

	// Soil depth
	// soil_depth is defined in global.ins
	// if zero then use input map of soil depth
	if (hydrology == RICHARDS) {
		if (negligible(depth_soil)) {
			gridcell.soil_depth = depth_map[c] * MM_PER_M;
		}
		else {
			gridcell.soil_depth = depth_soil;
		}
	}
	else {
		gridcell.soil_depth = SOILDEPTH_UPPER + SOILDEPTH_LOWER;
	}

	// Make sure soil depth is within ranges (mm)
	gridcell.soil_depth = min(SOILDEPTH_MAX, max(SOILDEPTH_MIN, gridcell.soil_depth));

	// Soil depth need to be of 100mm intervals
	double soil_depth_interval = round(gridcell.soil_depth / 100.0) * 100.0;

	if (!negligible(fabs(gridcell.soil_depth-soil_depth_interval))) {
		dprintf("Soil depth needs to be of 100mm intervals. Input %g Used %g\n", gridcell.soil_depth, soil_depth_interval);
		gridcell.soil_depth = soil_depth_interval;
	}

	fail("FIXME: figure out what K, b, and Psi_s are, and add to new type if necessary");
	// soiltype.K = soilpropmineral.K;
	// soiltype.b = soilpropmineral.b;
	// soiltype.Psi_s = soilpropmineral.Psi_s;

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

	int nr_lyr = 1;

	soiltype.sand_frac.resize(nr_lyr);
	soiltype.clay_frac.resize(nr_lyr);
	soiltype.coarse_frac.resize(nr_lyr);
	soiltype.pH.resize(nr_lyr);
	soiltype.text_depths.resize(nr_lyr);

	// exponent in percolation equation [k2; LPJF]
	// (Eqn 31, Haxeltine & Prentice 1996)
	// Changed from 4 to 2 (Sitch, Thonicke, pers comm 26/11/01)
	const double PERC_EXP = 2.0;

	int organic_soilcode = 8;

	if (soilcode<0 || soilcode>9) {
		fail("soil_parameters: invalid LPJ soil code (%d)", soilcode);
	}

	if (textured_soil) {
		soiltype.sand_frac[0] = data[soilcode][7];
		soiltype.clay_frac[0] = data[soilcode][8];
	}
	else {
		// Using fixed values from Parton et al. (2010)
		soiltype.sand_frac[0] = 0.28;
		soiltype.clay_frac[0] = 0.12;
	}

	soiltype.perc_base = data[soilcode][0];
	soiltype.perc_base_evap = data[soilcode][0];
	soiltype.perc_exp = ifdavehydrology ? PERC_EXP_DAVE : PERC_EXP;
	soiltype.perc_base_org = data[organic_soilcode][0];

	soiltype.thermdiff_0 = data[soilcode][2];
	soiltype.thermdiff_15 = data[soilcode][3];
	soiltype.thermdiff_100 = data[soilcode][4];
	soiltype.thermdiff_0_org = data[organic_soilcode][2];
	soiltype.thermdiff_15_org = data[organic_soilcode][3];
	soiltype.thermdiff_100_org = data[organic_soilcode][4];

	soiltype.awc_min = data[soilcode][1];
	soiltype.wp_min = data[soilcode][5];
	soiltype.wsats_min = data[soilcode][6];

	// Store the soilcode and new soil properties
	soiltype.soilcode = soilcode;
	soiltype.organic_frac[0] = data[soilcode][9];
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
