///////////////////////////////////////////////////////////////////////////////////////
/// \file externalinput.cpp
/// \brief Input code for land cover, management and other data from text files.
///				
/// \author Mats Lindeskog
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include <limits.h>

#include "externalinput.h"
#include "landcover.h"
#include "blaze.h"

/// Landcover area fraction input resolution used in the code to reject changes caused by rounding errors.
/// Initial value willl be modified by parsing input files or setting value in instruction file
double INPUT_RESOLUTION = 1.0e-6;
/// Reasonable guess of landcover area fraction input resolution in original data; used to take into account 
/// very small additional error added by rescaling area fractions so the sum is 1.0.
const double ORIGINAL_INPUT_RESOLUTION = 0.5e-6;

void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist) {

	/* Reads list of grid cells and (optional) description text from grid list file
	 * This file should consist of any number of one-line records in the format:
	 *   <longitude> <latitude> [<description>]
	 */

	double dlon, dlat;
	bool eof = false;
	xtring descrip;

	// Read list of grid coordinates and store in Coord object 'gridlist'

	FILE* in_grid = fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input", (char*)file_gridlist);

	while (!eof) {
		
		// Read next record in file
		eof =! readfor(in_grid, "f,f,a#", &dlon, &dlat, &descrip);

		// Deal with rounding errors associated with reading in text files
		dlon = roundoff(dlon, 6);
		dlat = roundoff(dlat, 6);

		if (!eof && !(dlon == 0.0 && dlat == 0.0)) { // ignore blank lines at end (if any)
			Coord& c = gridlist.createobj(); // add new coordinate to grid list

			c.lon = dlon;
			c.lat = dlat;
			c.descrip = descrip;
		}
	}
	fclose(in_grid);
	gridlist.firstobj();
}

MiscInput::~MiscInput() {

	delete[] forest_init_st;
}

void MiscInput::init(xtring file_gridlist) {

	ListArray_id<Coord> gridlist;
	read_gridlist(gridlist, file_gridlist);

	// See getdisturbance() for input file format and content of text files

	file_disturbance = param["file_disturbance"].str;
	if(file_disturbance != "")	{
		if(!disturbance.Open(file_disturbance, gridlist))
			fail("initio: could not open %s for input",(char*)file_disturbance);
		readdisturbance = true;
	}

	file_disturbance_st = param["file_disturbance_st"].str;
	if(file_disturbance_st != "")	{
		if(!disturbance_st.Open(file_disturbance_st, gridlist))
			fail("initio: could not open %s for input",(char*)file_disturbance_st);
		readdisturbance_st = true;
	}

	if (ifstorm) {
		file_windload = param["file_windload"].str;
		if (file_windload != "") {
			if (!windload.Open(file_windload, gridlist))
				fail("initio: could not open %s for input", (char*)file_windload);
			readwindload = true;
		}
	}

	if (param.isparam("file_elevation_gc_mean")) {
		file_elevation_gc_mean = param["file_elevation_gc_mean"].str;
		if (file_elevation_gc_mean != "")	{
			if(!elevation_gc_mean.Open(file_elevation_gc_mean, gridlist))
				fail("initio: could not open %s for input",(char*)file_elevation_gc_mean);
			readelevation_gc_mean = true;
		}
	}

	if (param.isparam("file_elevation_local")) {
		file_elevation_local = param["file_elevation_local"].str;
		if (file_elevation_local != "")	{
			if(!readelevation_gc_mean)
				fail("initio: Local elevation input not allowed without mean gridcell input in this branch!\n\n");
			if(!elevation_local.Open(file_elevation_local, gridlist))
				fail("initio: could not open %s for input",(char*)file_elevation_local);
			readelevation_local = true;
		}
	}

	if (param.isparam("file_elevation_st")) {
		file_elevation_st = param["file_elevation_st"].str;
		if(file_elevation_st != "")	{
			if(!elevation_st.Open(file_elevation_st, gridlist))
				fail("initio: could not open %s for input",(char*)file_elevation_st);
			readelevation_st = true;
		}
	}
	// Functionality not tested with run_landcover == false
	if(run_landcover) {
		file_dist_clearcut = param["file_dist_clearcut"].str;
		if(file_dist_clearcut != "")	{
			if(!dist_clearcut.Open(file_dist_clearcut, gridlist))
				fail("initio: could not open %s for input",(char*)file_dist_clearcut);
			readdist_clearcut = true;
		}

		file_dist_fire = param["file_dist_fire"].str;
		if(file_dist_fire != "")	{
			if(!dist_fire.Open(file_dist_fire, gridlist))
				fail("initio: could not open %s for input",(char*)file_dist_fire);
			readdist_fire = true;
		}

		file_dist_wind = param["file_dist_wind"].str;
		if(file_dist_wind != "")	{
			if(!dist_wind.Open(file_dist_wind, gridlist))
				fail("initio: could not open %s for input",(char*)file_dist_wind);
			readdist_wind = true;
		}

		file_dist_barkbeetle = param["file_dist_barkbeetle"].str;
		if(file_dist_barkbeetle != "")	{
			if(!dist_barkbeetle.Open(file_dist_barkbeetle, gridlist))
				fail("initio: could not open %s for input",(char*)file_dist_barkbeetle);
			readdist_barkbeetle = true;
		}

		file_storm_prob = param["file_storm_prob"].str;
		if(file_storm_prob != "")	{
			if(!storm_prob.Open(file_storm_prob, gridlist))
				fail("initio: could not open %s for input",(char*)file_storm_prob);
			readstorm_prob = true;
		}
	}

	forest_init_st = new TextInput::TimeDataD[nst];

	for (int i = 0; i < nst; i++) {
		StandType& st = stlist[i];
		if (st.file_treeinit_st != "") {
			if (!forest_init_st[i].Open(st.file_treeinit_st, gridlist))
				fail("initio: could not open %s for input", (char*)st.file_treeinit_st);
			readforest_init = true;
		}
	}
}

void MiscInput::load(double lon, double lat) {
	bool error = false;

	error |= loaddisturbance(lon, lat);
	error |= loaddist_clearcut(lon, lat);
	error |= loaddist_fire(lon, lat);
	error |= loaddist_wind(lon, lat);
	error |= loaddist_barkbeetle(lon, lat);
	error |= loadstorm_prob(lon, lat);
	error |= loadwindload(lon, lat);
	error |= loadelevation(lon, lat);
	error |= loadforest_init(lon, lat);

	if (error) {
		fail("loadmiscinput: error loading external input data for gridcell (lon, lat) = (%.6f, %.6f)",
			 lon, lat);
	}
}

bool MiscInput::loaddisturbance(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readdisturbance) { 
		if(!disturbance.Load(c)) {
			disterror = true;
		}
	}
	if(readdisturbance_st) { 
		if(!disturbance_st.Load(c)) {
			disterror = true;
		}
	}
	return disterror;
}

bool MiscInput::loaddist_clearcut(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readdist_clearcut) { 
		if(!dist_clearcut.Load(c)) {
			disterror = true;
		}
	}

	// st level input here

	return disterror;
}

bool MiscInput::loaddist_fire(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readdist_fire) { 
		if(!dist_fire.Load(c)) {
			disterror = true;
		}
	}

	// st level input here

	return disterror;
}

bool MiscInput::loaddist_wind(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readdist_wind) { 
		if(!dist_wind.Load(c)) {
			disterror = true;
		}
	}

	// st level input here

	return disterror;
}

bool MiscInput::loaddist_barkbeetle(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readdist_barkbeetle) { 
		if(!dist_barkbeetle.Load(c)) {
			disterror = true;
		}
	}

	// st level input here

	return disterror;
}

bool MiscInput::loadstorm_prob(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	// Not all gridcells have to be included in input file
	if(readstorm_prob) { 
		if(!storm_prob.Load(c)) {
			disterror = true;
		}
	}

	// st level input here

	return disterror;
}

bool MiscInput::loadwindload(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool disterror = false;

	if (readwindload) {
		if (!windload.Load(c)) {
			disterror = true;
		}
	}
	return disterror;
}

bool MiscInput::loadelevation(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool error = false;

	if(readelevation_gc_mean) { 
		if(!elevation_gc_mean.Load(c)) {
			error = true;
		}
	}
	if(readelevation_local) { 
		if(!elevation_local.Load(c)) {
			error = true;
		}
	}
	if(readelevation_st) { 
		if(!elevation_st.Load(c)) {
			error = true;
		}
	}
	return error;
}

bool MiscInput::loadforest_init(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool error = false;

	if (readforest_init) {
		for (int i = 0; i < nst; i++) {
			StandType& st = stlist[i];
			if (st.file_treeinit_st != "") {
				if (!forest_init_st[i].Load(c)) {
					error = true;	// skip this stand
				}
			}
		}
	}
	return error;
}

/// Read patch-destroying disturbance interval from file
/** This implementation uses disturbance interval input at gridcell or stand type level in standard text format input 
 *	files (see indata.h). The column for the gridcell disturbance interval has a header name of "Return" and columns
 *  for stand type disturbance have headers with the stand type names.
 */
void MiscInput::getdisturbance(Gridcell& gridcell) {

	// Retrieve disturbance for stand types from instruxtion file if text input not present,
	if(!date.year && !readdisturbance_st) {
		for(int i=0; i<nst; i++) {
			StandType& st = stlist[i];
			gridcell.st[st.id].distinterval_st = st.distinterval;
		}
	}

	int year = date.get_calendar_year();

	// Retrieve disturbance for grid cell from text input file. If not present, use global intruction file parameter.
	if(disturbance.isloaded())
		gridcell.distinterval_gc = disturbance.Get(year, "Return");
	else
		gridcell.distinterval_gc = distinterval;

	// Retrieve disturbance for stand types from text input file
	for(int i=0; i<nst; i++) {
		if(disturbance_st.isloaded()) {
			double dist = disturbance_st.Get(year, stlist[i].name, true);
			if(dist != NOTFOUND)
				gridcell.st[i].distinterval_st = dist;
		}
		// If not set at this stage, use gridcell value.
		if(gridcell.st[i].distinterval_st == 1.0e10)
			gridcell.st[i].distinterval_st = gridcell.distinterval_gc;
//		if(!date.year)
//			log_info("st %d distinterval_st =%f", i, gridcell.st[i].distinterval_st);
	}
}

/// Read patch-destroying clearcut disturbance biomass fraction from file
void MiscInput::getdist_clearcut(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(dist_clearcut.isloaded()) {
		double dist = dist_clearcut.Get(year, "clearcut", false, true);
		if(dist != NOTFOUND)
			gridcell.dist_clearcut_gc = dist;
		else
			fail("Satellite clearcut disturbance input file is missing \"clearcut\" column\n");
	}
}

/// Read patch-destroying fire disturbance biomass fraction from file
void MiscInput::getdist_fire(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(dist_fire.isloaded()) {
		double dist = dist_fire.Get(year, "fire", false, true);
		if(dist != NOTFOUND)
			gridcell.dist_fire_gc = dist;
		else
			fail("Satellite fire disturbance input file is missing \"fire\" column\n");
	}
}

/// Read patch-destroying wind disturbance biomass fraction from file
void MiscInput::getdist_wind(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(dist_wind.isloaded()) {
		double dist = dist_wind.Get(year, "wind", false, true);
		if(dist != NOTFOUND)
			gridcell.dist_wind_gc = dist;
		else
			fail("Satellite wind disturbance input file is missing \"wind\" column\n");
	}
}

/// Read patch-destroying barkbeetle disturbance biomass fraction from file
void MiscInput::getdist_barkbeetle(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(dist_barkbeetle.isloaded()) {
		double dist = dist_barkbeetle.Get(year, "barkbeetle", false, true);
		if(dist != NOTFOUND)
			gridcell.dist_barkbeetle_gc = dist;
		else
			fail("Satellite barkbeetle disturbance input file is missing \"barkbeetle\" column\n");
	}
}

/// Read patch-destroying wind disturbance probability from file
void MiscInput::getstorm_prob(Gridcell& gridcell) {
	if (!ifstorm) {
		return;
	}

	int year = date.get_calendar_year();

	if(storm_prob.isloaded()) {
		double stormprob = storm_prob.Get(year, "WDP");
		if(stormprob != NOTFOUND)
			gridcell.storm_prob_gc = stormprob;
		else
			fail("Storm probability input file is missing \"WDP\" column\n");
	}
}

void MiscInput::getwindload(Gridcell& gridcell) {

	// Retrieve windload for grid cell
	if (windload.isloaded()) {
		// Cycle the windload data before and after the data period
		int data_year;
		if (date.get_calendar_year() < windload.GetFirstyear())
			data_year = windload.GetnYears() - (((windload.GetFirstyear() - date.get_calendar_year() - 1) % windload.GetnYears()) + 1) + windload.GetFirstyear();
		else
			data_year = (date.get_calendar_year() - windload.GetFirstyear()) % windload.GetnYears() + windload.GetFirstyear();

		double wind_load = windload.Get(data_year, "Wind_load");
		if(wind_load != NOTFOUND)
			gridcell.windload = wind_load;
		else
			fail("Wind load input file is missing \"Wind_load\" column\n");
	}
	else {
		gridcell.windload = 0.0;
	}
}

void MiscInput::getelevation(Gridcell& gridcell) {

	// Year not really needed...
	int year = date.get_calendar_year();

	// Initialised to 0.0; this means missing value.
	if(elevation_gc_mean.isloaded()) {
		double elevation = elevation_gc_mean.Get(year, "elevation");
		if(elevation != NOTFOUND)
			gridcell.elevation_gc_mean = elevation;
		else
			fail("Gridcell elevation input file is missing \"elevation\" column\n");

		log_info("Mean elevation in text input file= %.0f", gridcell.elevation_gc_mean);
	}
	// Initialised to 0.0; this means missing value.
	if(elevation_local.isloaded()) {
		double elevation = elevation_local.Get(year, "elevation");
		if(elevation != NOTFOUND)
			gridcell.elevation_local = elevation;
		else
			fail("Local elevation input file is missing \"elevation\" column\n");

		log_info("Mean elevation = %.0f", gridcell.elevation_gc_mean);
		log_info("Local elevation = %.0f", gridcell.elevation_local);
	}

	// Retrieve elevation for stand types
	if(elevation_st.isloaded()) {
		for(int i=0; i<nst; i++) {	
			double elevation = elevation_st.Get(year, stlist[i].name, true);
//			log_info("st %s: elevation = %.0f", (char*)stlist[i].name, elevation);
			if(elevation != NOTFOUND)
				gridcell.st[i].elevation_st = elevation;
		}
	}
}

void MiscInput::getforest_init(Gridcell& gridcell) {

	int establishment_year_test = 2000;
	double dens_test = 0.01;
	bool test_run = false;
	bool print_init_info = true;

	for (int s=0; s<nst;s++) {
		StandType& st = stlist[s];
		Gridcellst& gcst = gridcell.st[s];

		if (forest_init_st[s].isloaded()) {

			gcst.prescribed_establishment = true;

			int ntrees_max = forest_init_st[s].GetnYears();
			if (ntrees_max > MAX_NESTABLISH_INDIV)
				fail("Input tree number = %d, increase value of MAX_NESTABLISH_INDIV (now %d)\n\n", ntrees_max, MAX_NESTABLISH_INDIV);
//			int est_year = (int)forest_init_st[s].Get(0, "census_year");	// census year should be the same for all trees in a stand type in a gricell
			int est_year = (int)forest_init_st[s].Get(0, "init_year");	// changed this for the European simulations

			if (est_year != NOTFOUND)
				gcst.establishment_list.establishment_year = est_year;
			else
				gcst.establishment_list.establishment_year = establishment_year_test;

			if(print_init_info)
				log_info("Tree cohorts initialised for stand type %s from input file year %d", (char*)st.name, gcst.establishment_list.establishment_year);

			for (int i = 0;i < ntrees_max; i++) {	// read values for all trees (for now: Year column 0, 1, ...)
				for (int p = 0;p < npft;p++) {
					Pft& pft = pftlist[p];
					char pftname[MAXNAMESIZE] = { NULL };

					forest_init_st[s].Get_string(i, "str_pft", pftname);
					if (!strcmp((char*)(pftlist[p].name), pftname)) {
						double biomass = forest_init_st[s].Get(i, "biomass");
						double dbh = forest_init_st[s].Get(i, "dbh");
						double dens = forest_init_st[s].Get(i, "dens");
						if (dens != NOTFOUND) {
							gcst.establishment_list.est_array[i].dens = dens / 10000.0;		// per ha to per m2
						}
						else {
							gcst.establishment_list.est_array[i].dens = dens_test;
						}
						int age = (int)forest_init_st[s].Get(i, "tree_age", true);	// Suppress warning if column not found.
						if(age != NOTFOUND && age > 0) {
							gcst.establishment_list.est_array[i].age = age;
						}
						else {
							 gcst.establishment_list.est_array[i].age = 0;
							/// Age set in establish_input() if column not found in input file or if value is <= 0
						}

						if(print_init_info) {
							log_info("Tree no.%d\t", i);
							log_info("%s: dbh = %.1f cm, dens = %.1f/ha, age = %d", pftname, dbh, dens, age);
						}
						strcpy(gcst.establishment_list.est_array[i].pft_name, pftname);
						gcst.establishment_list.est_array[i].biomass = biomass;
						gcst.establishment_list.est_array[i].diam = dbh;
					}
				}
			}
		}
	}
}

std::set<int> MiscInput::get_static_disturbances() const {
	if (file_disturb == "") {
		fail("file_disturb must be set when using disttype = DETERMINISTIC");
	}

	FILE* file = fopen(file_disturb, "r");
	if (!file) {
		fail("get_static_disturbances(): could not open %s for input: %s", (const char *)file_disturb, strerror(errno));
	}

	std::set<int> years;

	const size_t BUFSIZE = 1024;
	char buffer[BUFSIZE];
	while (fgets(buffer, BUFSIZE, file)) {
		// Convert to integer.
		char *end;
		errno = 0;
		unsigned long year = strtoul(buffer, &end, 10);
		if (errno != 0) {
			fail("get_static_disturbances(): failed to parse line '%s': %s", buffer, strerror(errno));
		}
		if (end == buffer) {
			fail("get_static_disturbances(): failed to parse line '%s': line cannot be converted to integer", buffer);
		}
		if (year > INT_MAX) {
			fail("get_static_disturbances(): year %lu is out of range", year);
		}
		years.insert((int)year);
	}

	int err = fclose(file);
	if (err) {
		printf("WARNING: Failed to close file_disturb: %s", (const char *)file_disturb);
	}

	return years;
}

void MiscInput::getmiscinput_static(Gridcell& gridcell) {

	getelevation(gridcell);
	getforest_init(gridcell);
	if (ifdisturb && disttype == DETERMINISTIC) {
		gridcell.static_disturbances = get_static_disturbances();
	}
}

void MiscInput::getmiscinput_yearly(Gridcell& gridcell) {

	if (date.day) {
		return;
	}
	getdisturbance(gridcell);
	getdist_clearcut(gridcell);
	getdist_fire(gridcell);
	getdist_wind(gridcell);
	getdist_barkbeetle(gridcell);
	getstorm_prob(gridcell);
	getwindload(gridcell);
}

LandcoverInput::LandcoverInput()
	: nyears_cropland_ramp(0), input_precision_force(0) {

	declare_parameter("minimizecftlist", &minimizecftlist, "Whether pfts not in crop fraction input file are removed from pftlist (0,1)");
	declare_parameter("nyears_cropland_ramp", &nyears_cropland_ramp, 0, 10000,
		"Number of years to increase cropland fraction linearly from 0 to first year's value");
	declare_parameter("frac_fixed_default_crops", &frac_fixed_default_crops,
		"Whether to use all active crop stand types (0) or only stand types with suitable rainfed crops (based on crop pft tb and gridcell latitude) (1) when using fixed crop fractions");
	declare_parameter("input_precision_force", &input_precision_force, 3, 15, "Precision of landcover area fraction input files");
}

void LandcoverInput::init(xtring file_gridlist) {

	if(!run_landcover)
		return;

	int input_precision_parsed = 0;
	int input_precision_parsed_grossLUC = 0;
	int input_precision_parsed_grossSTC = 0;
	int input_precision_use = 0;
	ListArray_id<Coord> gridlist;
	read_gridlist(gridlist, file_gridlist);

	all_fracs_const = true;	/* If any of the opened files have yearly data, all_fracs_const will be set to false and
							 * landcover_dynamics will call get_landcover() each year
							 */

	// No peatland allowed when using the two layer soil so treat the peatland fraction as natural
	if (hydrology == GERTEN && run[PEATLAND]) {
		fail("LandcoverInput::init(): do not set run_peatland to 1 in landcover.ins if hydrology is set to GERTEN");
	}

	// Active peatland and hydrology set to Richards, then check if peat duration input file exist
	if (hydrology == RICHARDS && run[PEATLAND]) {
		// Retrieve file name for peat duration and start date
		file_peat = param["file_peat"].str;

		if (file_peat != "") {
			if (!Peatdata.Open(file_peat, gridlist)) {
				fail("LandcoverInput::init(): could not open %s for input", (char*)file_peat);
			}
			else {
				peatduration = true;
			}
		}
	}


	bool LU_input_present = false;

	//Retrieve file names for landcover files and open them if static values from ins-file are not used

	bool openLUfile = false;

	for(int i=0; i<NLANDCOVERTYPES; i++) {
		if(run[i] && i != NATURAL)
			openLUfile = true;
	}

	if (openLUfile) {

		// Retrieve file name for landcover fraction file and open them unless empty string and static equal-size values are used.
		file_lu=param["file_lu"].str;

		if(file_lu != "")	{
			if(!LUdata.Open(file_lu, gridlist)) {
				fail("LandcoverInput::init(): could not open %s for input",(char*)file_lu);
			}
			else {
				lcfrac_fixed = false;

				if(LUdata.GetFormat() == TextInput::LOCAL_YEARLY || LUdata.GetFormat() == TextInput::GLOBAL_YEARLY)
					all_fracs_const = false;				//Set all_fracs_const to false if yearly data

				// Avoid large number of output files
				if(printseparatestands && LUdata.GetNCells() > 100) {
					printseparatestands = false;
					log_info("Stand output disabled for simulations with land cover input files with more than 100 gridcells to avoid a large number of output files\nLandcover input file has data for %d gridcells", LUdata.GetNCells());
				}

				input_precision_parsed = LUdata.GetPrecision();
//				if(input_precision_parsed == -1 && !input_precision_force)
//					fail("Scientific notation area fraction input only supported when input_precision_force is defined\n");

				LU_input_present = true;
			}
		}
	}

	if (!lcfrac_fixed) {
		// Open file for LUC transitions
		file_grossLUC = param["file_grossLUC"].str;
		file_grossSTC = param["file_grossSTC"].str;
		if(file_grossLUC != "") {
			if(!grossLUC.Open(file_grossLUC, gridlist))
				fail("initio: could not open %s for input",(char*)file_grossLUC);
			all_fracs_const = false;	// needed for some tests with static net land cover are fractions
			gross_land_transfer = LC_GROSS_INPUT;

			input_precision_parsed_grossLUC = grossLUC.GetPrecision();
//			if(input_precision_parsed == -1 && !input_precision_force)
//				fail("Scientific notation area fraction input only supported when input_precision_force is defined\n");
		}
		else if(file_grossSTC != "") {
			if(!grossSTC.Open(file_grossSTC, gridlist))
				fail("initio: could not open %s for input",(char*)file_grossSTC);
			all_fracs_const = false;	// needed for some tests with static net stand type area fractions
			gross_land_transfer = ST_GROSS_INPUT;

			input_precision_parsed_grossSTC = grossSTC.GetPrecision();
//			if(input_precision_parsed == -1 && !input_precision_force)
//				fail("Scientific notation area fraction input only supported when input_precision_force is defined\n");
		}
		else {
			gross_land_transfer = NO_GROSS_INPUT;
		}
	}
	else {
		gross_land_transfer = NO_GROSS_INPUT;
	}

	//Retrieve file name for stand type fraction files and open them if static equal-size values are not used.
	file_lu_st[CROPLAND] = param["file_lucrop"].str;
	file_lu_st[PASTURE] = param["file_lupasture"].str;
	file_lu_st[NATURAL] = param["file_lunatural"].str;
	file_lu_st[FOREST] = param["file_luforest"].str;

	if(file_lu_st[CROPLAND] != "" || file_lu_st[PASTURE] != "" || file_lu_st[NATURAL] != "" || file_lu_st[FOREST] != "")
		LU_input_present = true;

	int input_precision_parsed_st_max = 0;
	bool st_input = false;

	if(LU_input_present) {
		log_info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		log_info("Landcover text input information:");
	}

	for(int lc=0; lc<NLANDCOVERTYPES; lc++) {
		if(run[lc] && file_lu_st[lc] != "")	{
			int input_precision_parsed_st = 0;
			if(!st_data[lc].Open(file_lu_st[lc], gridlist)) {
				fail("could not open %s for input",(char*)file_lu_st[lc]);
			}
			else {
				st_input = true;
				frac_fixed[lc] = false;
				if(st_data[lc].GetFormat() == TextInput::LOCAL_YEARLY || st_data[lc].GetFormat() == TextInput::GLOBAL_YEARLY) 
					all_fracs_const = false;				// Set all_fracs_const to false if yearly data

				input_precision_parsed_st = st_data[lc].GetPrecision();
//				if(input_precision_parsed_st == -1 && !input_precision_force)
//					fail("Scientific notation area fraction input only supported when input_precision_force is defined\n");
				if(input_precision_parsed_st > input_precision_parsed_st_max)
					input_precision_parsed_st_max = input_precision_parsed_st;
				if(input_precision_parsed != input_precision_parsed_st) {
					log_info("Stand type fraction input file %s\n has different input precision (%d) from LC file (%d)",
						(char*)file_lu_st[lc], input_precision_parsed_st, input_precision_parsed);
				}
			}
		}
	}

	// Increase precision if stand type input file has normalised data.
	bool found = false;
	stlist.firstobj();
	while (stlist.isobj && !found) {
		StandType& st = stlist.getobj();

		if(file_lu_st[st.landcover] != "" && st_data[st.landcover].NormalisedData()) {
			input_precision_parsed_st_max = input_precision_parsed_st_max + input_precision_parsed;
			found = true;
			log_info("Normalised stand type fraction input file %s:\n input precision increased to (%d) ",
				(char*)file_lu_st[st.landcover], input_precision_parsed_st_max);
		}
		stlist.nextobj();
	}

	bool use_STC_prec = false;
	bool use_grossLUC_prec = false;
	bool use_grossSTC_prec = false;

	if(st_input) {
		if(input_precision_parsed > input_precision_parsed_st_max) {
			log_info("Using LC input precision (%d)", input_precision_parsed);
		}
		else if(input_precision_parsed < input_precision_parsed_st_max) {
			log_info("Using ST input precision (%d)", input_precision_parsed_st_max);
			input_precision_parsed = input_precision_parsed_st_max;
			use_STC_prec = true;
		}
	}

	if(input_precision_parsed_grossLUC > input_precision_parsed) {
		input_precision_parsed = input_precision_parsed_grossLUC;
		use_grossLUC_prec = true;
	}
	if(input_precision_parsed_grossSTC > input_precision_parsed) {
		input_precision_parsed = input_precision_parsed_grossSTC;
		use_grossSTC_prec = true;
	}
	if(use_grossSTC_prec)
		log_info("Using gross STC input precision (%d)", input_precision_parsed);
	else if(use_grossLUC_prec)
		log_info("Using gross LUC input precision (%d)", input_precision_parsed);
	else if(use_STC_prec)
		log_info("Using ST input precision (%d)", input_precision_parsed_st_max);
	else
		log_info("Using LC input precision (%d)", input_precision_parsed);

	input_precision_use = input_precision_parsed;

	// Overwrite parsed input precision with instruction file parameter, if set.
	if(input_precision_force)
		input_precision_use = input_precision_force;

	// Limit precision range to 3-14.
	input_precision_use = max(min(input_precision_use, 14), 3);

	INPUT_RESOLUTION = pow(10.0, -input_precision_use);

	// Take into account very small additional error added by rescaling area fractions so the sum is 1.0.
	INPUT_RESOLUTION = INPUT_RESOLUTION - INPUT_RESOLUTION * ORIGINAL_INPUT_RESOLUTION;

	if(input_precision_force)
		log_info("Manual instruction file input precision=%d", input_precision_force);
	if(LU_input_present) {
		log_info("Land cover fraction input precision parsed=%d, used=%d", input_precision_parsed, input_precision_use);
		log_info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	}

	if(!frac_fixed[CROPLAND]) {

		TextInput::TimeDataD& CFTdata = st_data[CROPLAND];
		bool do_minimize = false;

		// Remove crop stand types from stlist that always have zero area fraction in all cells in grid list
		if(minimizecftlist && gridlist.nobj < 100) {	// Reduce the risk of accidentally using minimized cft lists when using split gridlists.
			log_info("WARNING: minimizecftlist is activated in crop.ins.\n Restart will NOT WORK properly if you are running on multiple processors. ");
			CFTdata.CheckIfPresent(gridlist);
			do_minimize = true;
		}

		int n=0;
		stlist.firstobj();
		while(stlist.isobj) {	

			StandType& st = stlist.getobj();
			bool remove = false;

			if(st.landcover == CROPLAND) {
				if(do_minimize)
					remove = !CFTdata.item_has_data(st.name);
				else
					remove = !CFTdata.item_in_header(st.name);
			}

			if(remove) {
				n+=1;
				stlist.killobj();
				nst--;
				nst_lc[CROPLAND]--;
			}
			else {
				st.id-=n;
				stlist.nextobj();
			}			
		}
	}

	// Remove pft:s from pftlist that are not grown in simulated stand types
	int n = 0;
	pftlist.firstobj();
	while(pftlist.isobj) {	

		bool remove = false;
		Pft& pft = pftlist.getobj();

		if(pft.landcover == CROPLAND) {

			bool found = false;

			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();

				if(st.pftinrotation(pft.name) >= 0) {
					found = true;
					break;
				}
				stlist.nextobj();
			}

			if(!found)
				remove = true;
		}

		if(remove && !(pft.isintercropgrass && ifintercropgrass)) {
			n += 1;
			pftlist.killobj();
			npft--;
		}
		else {
			pft.id -= n;
			pftlist.nextobj();
		}			
	}

	// Remove mt:s from mtlist that are not used in simulated stand types
	int m = 0;
	mtlist.firstobj();
	while(mtlist.isobj) {	

		bool remove = false;
		ManagementType& mt = mtlist.getobj();

		bool found = false;

		stlist.firstobj();
		while(stlist.isobj) {
			StandType& st = stlist.getobj();

			if(st.mtinrotation(mt.name) >= 0) {
				found = true;
				break;
			}
			stlist.nextobj();
		}

		if(!found)
			remove = true;

		if(remove) {
			log_info("Removing mt %s", (char*)mt.name);
			m += 1;
			mtlist.killobj();
			nmt--;
		}
		else {
			// Set new id for remaining management types, corresponding to place in mtlist.
			mt.id -= m;

			// Also set new id for corresponding st.management copy (identical to mtlist[st.management.id]).
			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();
				if(st.management.name == mt.name) {
					st.management.id = mt.id;
				}
				stlist.nextobj();
			}

			mtlist.nextobj();
		}			
	}

	gridlist.killall();
}

bool LandcoverInput::loadlandcover(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool LUerror = false;


	/* Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
	 * transferred to gridcell.landcover.frac each year in getlandcover()
	 */

	if (!lcfrac_fixed) {

		bool loadLU = false;

		for(int i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				loadLU = true;
		}

		if (loadLU) {
			// Load landcover area fraction data from input file to data object
			if (!LUdata.Load(c)) {
				log_info("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.",c.lon,c.lat);
				LUerror = true;		// skip this stand
			}
		}

		//Read LUC transitions
		if(gross_land_transfer == LC_GROSS_INPUT && !LUerror) {

			if(!grossLUC.Load(c)) {
				log_info("Data for %.3f,%.3f missing in gross LUC transitions input file.",c.lon,c.lat);
				gross_LC_input_present = false;
			}
			else {
				gross_LC_input_present = true;
			}
		}
		if(gross_land_transfer == ST_GROSS_INPUT && !LUerror) {

			if(!grossSTC.Load(c)) {
				log_info("Data for %.3f,%.3f missing in gross STC transitions input file.",c.lon,c.lat);
				gross_ST_input_present = false;
			}
			else {
				gross_ST_input_present = true;
			}
		}
	}

	for(int lc=0; lc<NLANDCOVERTYPES && !LUerror; lc++) {

		if(run[lc] && !frac_fixed[lc]) {

			/* Stand type fraction data: read from stand type fraction file; dynamic, so data for all years are loaded to st_data object and 
			 * transferred to gridcell.landcover.frac[] each year in getlandcover()	
			 */

			if(!st_data[lc].Load(c)) {
				log_info("Problems with stand type fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.",c.lon,c.lat);
				LUerror = true;	// skip this stand
			}
		}
	}

	return LUerror;
}

void LandcoverInput::getlandcover(Gridcell& gridcell) {

	const char* lcnames[]={"URBAN", "CROPLAND", "PASTURE", "FOREST", "NATURAL", "PEATLAND", "BARREN"};
	double sum_tot=0.0, sum_active=0.0;
	double sum_tot_nextyear=0.0, sum_active_nextyear=0.0;
	Landcover& lc = gridcell.landcover;

	// Calender year of start of simulation (after spinup)
	int first_historic_year = date.first_calendar_year + nyear_spinup;

	/* Use values for first historic year during spinup period, unless data exist before firsthistyear
	 * Use values for last historic year during the time after that
	 * This is handled by the text input class.
	 */
	int year = date.get_calendar_year();
	int year_saved = year;
	int first_reduction_year = first_historic_year - nyear_spinup + (int)(SOLVESOMCENT_SPINEND * (nyear_spinup - freenyears) + freenyears) + 1;
	if(year < first_reduction_year && year > LUdata.GetFirstyear()) {
		if(year == first_reduction_year - 1 && !nyears_cropland_ramp)
			log_info("Land cover change before soil spinup is not allowed, first lcc year will be %d. Using lc fractions for %d earlier.",
				first_reduction_year, first_reduction_year);
		year = first_reduction_year;
	}

	//Save old fraction values:									
	for(int i=0; i<NLANDCOVERTYPES; i++)
		lc.frac_old[i] = lc.frac[i];
	for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
		Gridcellst& gcst = gridcell.st[i];
		gcst.frac_old = gcst.frac;
	}


	if(lcfrac_fixed) {		// Set land covers to equal fractions
		if(date.year == 0) {	// Year 0: called by landcover_init

			int nactive_landcovertypes = 0;

			for(int i=0; i<NLANDCOVERTYPES; i++) {
				if(run[i])
					nactive_landcovertypes++;
			}

			for(int i=0;i<NLANDCOVERTYPES;i++) {
				lc.frac[i] = 1.0 * run[i] / (double)nactive_landcovertypes;	// only set fractions that are active
				lc.frac_nextyear[i] = 1.0 * run[i] / (double)nactive_landcovertypes;	// only set fractions that are active
				sum_active += lc.frac[i];
				sum_tot = sum_active;
				sum_active_nextyear += lc.frac_nextyear[i];
				sum_tot_nextyear = sum_active_nextyear;
			}		
		}
	}
	else {	// landcover area fractions are read from input file(s)	

		bool printyear = year >= LUdata.GetFirstyear() && LUdata.GetFirstyear() >= 0;	
		bool getLU = false;
		double sum = 0.0;

		for(int i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				getLU = true;
		}

		if(getLU) {	

			if(LUdata.Get(year, 0) < 0.0) {		// Missing data (negative values)
				if(date.year == 1)
					log_info("Missing landcover fraction data for year %d, natural vegetation fraction set to 1.0", year);
				for(int i=0;i<NLANDCOVERTYPES;i++) {
					lc.frac[i] = 0.0;
					lc.frac_nextyear[i] = 0.0;
				}
				lc.frac[NATURAL] = 1.0;
				sum_active = 1.0;
				lc.frac_nextyear[NATURAL] = 1.0;
				sum_active_nextyear = 1.0;
			}
			else {

				for(int i=0; i<NLANDCOVERTYPES; i++) {

					if(run[i]) {
						double lcfrac = 0.0;
						double lcfrac_nextyear = 0.0;
						switch(i)
						{
						case URBAN:
							lcfrac = LUdata.Get(year,"URBAN");
							lcfrac_nextyear = LUdata.Get((year + 1),"URBAN");
							break;
						case CROPLAND:
							lcfrac = LUdata.Get(year,"CROPLAND");
							lcfrac_nextyear = LUdata.Get((year + 1),"CROPLAND");
							break;
						case PASTURE:
							lcfrac = LUdata.Get(year,"PASTURE");
							lcfrac_nextyear = LUdata.Get((year + 1),"PASTURE");
							break;
						case FOREST:
							lcfrac = LUdata.Get(year,"FOREST");
							lcfrac_nextyear = LUdata.Get((year + 1),"FOREST");
							break;
						case NATURAL:
							lcfrac = LUdata.Get(year,"NATURAL");
							lcfrac_nextyear = LUdata.Get((year + 1),"NATURAL");
							break;
						case PEATLAND:
							lcfrac = LUdata.Get(year,"PEATLAND");
							lcfrac_nextyear = LUdata.Get((year + 1),"PEATLAND");
							break;
						case BARREN:
							lcfrac = LUdata.Get(year,"BARREN");
							lcfrac_nextyear = LUdata.Get((year + 1),"BARREN");
							break;
						default:
							if(date.year == 0)
								log_info("Modify code to deal with landcover input!");
						}

						if(lcfrac == NOTFOUND) {	// land cover not found in input file
							lcfrac = 0.0;
						}
						else if(lcfrac < 0.0) {	// discard unreasonable values
							if(printyear)
								log_info("WARNING ! landcover fraction size out of limits, set to 0.0");
							lcfrac = 0.0;
						}
						else if(lcfrac > 1.01) {	// discard unreasonable values
							if(printyear)
								log_info("WARNING ! %d landcover %d fraction size %f out of limits, set to 1.0",
									date.get_calendar_year(), i, lcfrac);
							lcfrac = 1.0;
						}
						lc.frac[i] = lcfrac;
						sum_tot += lcfrac;
						sum_active += run[i] * lc.frac[i];
/// frac_nextyear used for tree init procedure
						if(lcfrac_nextyear == NOTFOUND) {	// land cover not found in input file
							lcfrac_nextyear = 0.0;
						} 
						else if(lcfrac_nextyear < 0.0 ) {	// discard unreasonable values
							if(printyear)
								log_info("WARNING ! %d landcover %d %s fraction size %f out of limits, set to 0.0", 
									date.get_calendar_year(), i, lcnames[i], lcfrac);
							lcfrac_nextyear = 0.0;
						} 
						else if(lcfrac_nextyear > 1.01) {	// discard unreasonable values
							if(printyear)
								log_info("WARNING ! %d landcover %d %s fraction size %f out of limits, set to 1.0", 
									date.get_calendar_year(), i, lcnames[i], lcfrac);
							lcfrac_nextyear = 1.0;
						}
						lc.frac_nextyear[i] = lcfrac_nextyear;
						sum_tot_nextyear += lcfrac_nextyear;
						sum_active_nextyear += run[i] * lc.frac_nextyear[i];
////////////////////////////
					}
				}

				if (grassforcrop) {
					lc.frac[PASTURE]+=lc.frac[CROPLAND];
					lc.frac[CROPLAND]=0.0;
				}
				if(sum_tot != 1.0 && sum_tot != 0.0) {		// Check input data, rescale if sum !=1.0	

					sum_active = 0.0;		// reset sum of active landcover fractions

					if(sum_tot < 0.99 || sum_tot > 1.01) {
						if(printyear) {
							log_info("WARNING ! landcover fraction sum is %4.2f for input year %d", sum_tot, year);
							log_info("Rescaling landcover fractions !");
						}
					}

					double sum_dynamic = sum_tot;
					double sum_dynamic_adjusted = 1.0;
					if (no_barren_frac_corr) {
						sum_dynamic -= lc.frac[BARREN];
						sum_dynamic_adjusted -= lc.frac[BARREN];
					}
					if (no_peatland_frac_corr) {
						sum_dynamic -= lc.frac[PEATLAND];
						sum_dynamic_adjusted -= lc.frac[PEATLAND];
					}
					
					// if only BARREN and/or PEATLAND are set, it can happen that sum_dynamic_adjusted is zero
					if (negligible(sum_dynamic, -14)) {
						sum_dynamic = sum_tot; // prevent div by zero
					}
					if (negligible(sum_dynamic_adjusted, -14)) {
						sum_dynamic_adjusted = 1.0; // prevent div by zero
					}

					for(int i=0; i<NLANDCOVERTYPES; i++) {
						if(no_barren_frac_corr||no_peatland_frac_corr) {
							bool flag_adjust_lc = true;
							if (i == BARREN && no_barren_frac_corr) {
								flag_adjust_lc = false;
							}
							if (i == PEATLAND && no_peatland_frac_corr) {
								flag_adjust_lc = false;
							}
							if(flag_adjust_lc) {
								lc.frac[i] /= sum_dynamic / sum_dynamic_adjusted;
							}
						}
						else {
							lc.frac[i] /= sum_tot;
						}
						sum_active += lc.frac[i];
					}
				}
/// frac_nextyear used for tree init procedure
				if (grassforcrop) {
					lc.frac_nextyear[PASTURE]+=lc.frac_nextyear[CROPLAND];
					lc.frac_nextyear[CROPLAND]=0.0;
				}
				if(sum_tot_nextyear != 1.0 && sum_tot_nextyear != 0.0) {		// Check input data, rescale if sum !=1.0	

					sum_active_nextyear = 0.0;		// reset sum of active landcover fractions

					if(sum_tot_nextyear < 0.99 || sum_tot_nextyear > 1.01) {
						if(printyear) {
							log_info("WARNING ! next year's landcover fraction sum is %4.2f for input year %d", sum_tot_nextyear, (year + 1));
							log_info("Rescaling landcover fractions !");
						}
					}

					double sum_dynamic_nextyear = sum_tot_nextyear - lc.frac_nextyear[BARREN];
					double sum_dynamic_adjusted_nextyear = 1.0 - lc.frac_nextyear[BARREN];

					for(int i=0; i<NLANDCOVERTYPES; i++) {
						if(no_barren_frac_corr) {
							if(i != BARREN) {
								lc.frac_nextyear[i] /= sum_dynamic_nextyear / sum_dynamic_adjusted_nextyear;
							}
						}
						else {
							lc.frac_nextyear[i] /= sum_tot_nextyear;
						}
						sum_active_nextyear += lc.frac_nextyear[i];
					}
				}

//////////////////////////////////////////////
			}
		}
		else {
			lc.frac[NATURAL] = 0.0;
			lc.frac_nextyear[NATURAL] = 0.0;
		}

		/* NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed.
		 * if landcover types are turned off in the instruction file, or if more landcover types
		 * are added in other input files, can be either less or more than 1.0
		 */
		if(!negligible(sum_active - 1.0, -14))	{

			if(date.year == 0)
				log_info("Landcover fraction sum not 1.0 !");

			if(run[NATURAL]) {	// Transfer landcover areas not simulated to NATURAL fraction, if simulated.		
				if(date.year == 0) {
					if(sum_active < 1.0)
						log_info("Inactive fractions (%4.3f) transferred to NATURAL fraction.", 1.0-sum_active);
					else
						log_info("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).",
							sum_active-1.0, lc.frac[NATURAL]);
				}

				lc.frac[NATURAL] += (1.0 - sum_active);
				// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					log_info("New NATURAL fraction is %4.3f.", lc.frac[NATURAL]);

				sum_active = 1.0;		// sum_active should now be 1.0

				if(lc.frac[NATURAL] < 0.0) {	// If new landcover type fraction is bigger than the natural fraction
												// (something wrong in the distribution of input file area fractions)						
					if(date.year == 0)
						log_info("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.");

					sum_active -= lc.frac[NATURAL];	// fraction not possible to transfer moved back to sum_active, which
													// will now be >1.0 again
					lc.frac[NATURAL] = 0.0;

					for(int i=0; i<NLANDCOVERTYPES; i++) {
						lc.frac[i] /= sum_active;		// fraction rescaled to unity sum
						if(run[i] && date.year == 0)
							log_info("Landcover type %d fraction is %4.3f", i, lc.frac[i]);
					}
				}
			}
			else {
				if(date.year == 0)
					log_info("Non-unity fraction sum retained.");				// let sum remain non-unity
			}
		}

/// frac_nextyear used for tree init procedure
		if(!negligible(sum_active_nextyear - 1.0, -14))	{

//			if(date.year == 0)
//				log_info("Landcover fraction sum not 1.0 !");

			if(run[NATURAL]) {	// Transfer landcover areas not simulated to NATURAL fraction, if simulated.		
//				if(date.year == 0) {
//					if(sum_active < 1.0)
//						log_info("Inactive fractions (%4.3f) transferred to NATURAL fraction.", 1.0-sum_active);
//					else
//						log_info("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).",
//							sum_active-1.0, lc.frac[NATURAL]);
//				}

				lc.frac_nextyear[NATURAL] += (1.0 - sum_active_nextyear);
				// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
//				if(date.year==0)
//					log_info("New NATURAL fraction is %4.3f.", lc.frac[NATURAL]);

				sum_active_nextyear = 1.0;		// sum_active should now be 1.0

				if(lc.frac_nextyear[NATURAL] < 0.0) {	// If new landcover type fraction is bigger than the natural fraction
												// (something wrong in the distribution of input file area fractions)						
//					if(date.year == 0)
//						log_info("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.");

					sum_active_nextyear -= lc.frac_nextyear[NATURAL];	// fraction not possible to transfer moved back to sum_active, which
													// will now be >1.0 again
					lc.frac_nextyear[NATURAL] = 0.0;

					for(int i=0; i<NLANDCOVERTYPES; i++) {
						lc.frac_nextyear[i] /= sum_active_nextyear;		// fraction rescaled to unity sum
//						if(run[i] && date.year == 0)
//							log_info("Landcover type %d fraction is %4.3f", i, lc.frac[i]);
					}
				}
			}
			else {
//				if(date.year == 0)
//					log_info("Non-unity fraction sum retained.");				// let sum remain non-unity
			}
		}

//////////////////////////////////////////////

		if(nyears_cropland_ramp) {

			bool doramp = false;
			int firstyear;
			if(LUdata.GetFirstyear() >= 0) {
				if(year_saved < LUdata.GetFirstyear()) {
					doramp = true;
					firstyear = LUdata.GetFirstyear();
				}
			}
			else {			
				if(year_saved < first_historic_year) {		
					doramp = true;
					firstyear = first_historic_year;
				}
			}

			if(doramp) {
				int max_ramp_years = firstyear - first_reduction_year;
				if(nyears_cropland_ramp > max_ramp_years && year_saved == first_reduction_year)
					log_info("Requested cropland ramp period too long for given nyear_spinup. Maximum is %d.", max_ramp_years);

				double reduce_cropland = min((double)(firstyear - year_saved) / min(nyears_cropland_ramp, max_ramp_years), 1.0) * lc.frac[CROPLAND];
				lc.frac[CROPLAND] -= reduce_cropland;
				lc.frac[NATURAL] += reduce_cropland;
				if(year_saved == firstyear -1)
					log_info("Cropland area fraction ramp from 0 to %.3f during period %d to %d",
						LUdata.Get(firstyear,"CROPLAND"), max(first_reduction_year, firstyear - nyears_cropland_ramp), firstyear-1);
			}
		}	
	}


	// Set fractions for static stand types
	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		Gridcellst& gcst = gridcell.st[st.id];

		gcst.frac = 0.0;
		if(nst_lc[st.landcover] == 1)
			gcst.frac = 1.0;
		else if(frac_fixed[st.landcover] && (st.landcover != CROPLAND || !frac_fixed_default_crops))
			gcst.frac = 1.0 / (double)nst_lc[st.landcover];
/// frac_nextyear used for tree init procedure
		gcst.frac_nextyear = 0.0;
		if(nst_lc[st.landcover] == 1)
			gcst.frac_nextyear = 1.0;
		else if(frac_fixed[st.landcover] && (st.landcover != CROPLAND || !frac_fixed_default_crops))
			gcst.frac_nextyear = 1.0 / (double)nst_lc[st.landcover];
///
		stlist.nextobj();
	}

	// Get fractions for dynamic stand types within a land cover

	for(int lc=0; lc<NLANDCOVERTYPES; lc++) {

		if(run[lc] && gridcell.landcover.frac[lc] && (!frac_fixed[lc] || lc == CROPLAND && frac_fixed_default_crops)) {

			double sum = 0.0;
			double sum_nextyear = 0.0;

			bool printyear = year >= st_data[lc].GetFirstyear() && st_data[lc].GetFirstyear() >= 0;

			if(lc == CROPLAND) {
				sum = get_crop_fractions(gridcell, year, st_data[CROPLAND], sum_tot);
			}
			else {

				if(year == st_data[lc].GetFirstyear() + st_data[lc].GetnYears())
					log_info("Last year of %s st fraction data used from year %d and onwards", lcnames[lc], year);

				for(int i=0; i<nst; i++) {
					if(stlist[i].landcover == lc)	{

						double frac = st_data[lc].Get(year,stlist[i].name);

						if(frac == NOTFOUND) {	// stand type not found in input file
							frac = 0.0;
						}
						else if(frac < 0.0)	{	// correct unreasonable values
							if(printyear && gridcell.landcover.frac[lc])
								log_info("WARNING ! %s fraction size out of limits, set to 0.0", lcnames[lc]);
							frac = 0.0;
						}
						else if(frac > 1.01) {	// correct unreasonable values
							if(printyear)
								log_info("WARNING ! %s fraction size out of limits, set to 1.0", lcnames[lc]);
							frac = 1.0;
						}
						if(!st_data[lc].NormalisedData())
							frac /= sum_tot;	// Scale by same factor as land cover fractions !
						sum += gridcell.st[i].frac = frac;

/// frac_nextyear used for tree init procedure
						double frac_nextyear = st_data[lc].Get((year + 1),stlist[i].name);

						if(frac_nextyear == NOTFOUND) {	// stand type not found in input file
							frac_nextyear = 0.0;
						}
						else if(frac_nextyear < 0.0)	{	// correct unreasonable values
							if(printyear && gridcell.landcover.frac[lc])
								log_info("WARNING ! %s fraction size out of limits, set to 0.0", lcnames[lc]);
							frac_nextyear = 0.0;
						}
						else if(frac_nextyear > 1.01) {	// correct unreasonable values
							if(printyear)
								log_info("WARNING ! %s fraction size out of limits, set to 1.0", lcnames[lc]);
							frac_nextyear = 1.0;
						}
						if(!st_data[lc].NormalisedData())
							frac_nextyear /= sum_tot_nextyear;	// Scale by same factor as land cover fractions !
						sum_nextyear += gridcell.st[i].frac_nextyear = frac_nextyear;
//						gridcell.st[i].frac_nextyear = frac_nextyear;
					}
				}

				if(sum == 0.0) {
					fail("%s st fraction sum is zero. Check input data\n", lcnames[lc]);
				}
			}

			if(sum > 0.0) {

				// rescale active fractions so sum is 1.0 (or land cover fraction value)

				double ratio = sum;

				if(!st_data[lc].NormalisedData())
					ratio /= gridcell.landcover.frac[lc];

				stlist.firstobj();
				while(stlist.isobj) {
					StandType& st = stlist.getobj();
					Gridcellst& gcst = gridcell.st[st.id];
					if(st.landcover == lc) {
						gcst.frac_old_orig = gcst.frac;
						gcst.frac /= ratio;
					}
					stlist.nextobj();
				}

				if((ratio < 0.99 || ratio > 1.01) && printyear) {	// warn if sum is significantly different from 1.0 
					log_info("WARNING ! %s st fraction sum is %7.5f for input year %d while LC fraction is %7.5f", 
						(char*)lcnames[lc], sum, year, gridcell.landcover.frac[lc]);
					log_info("Rescaling %s  fractions year %d !", lcnames[lc], date.get_calendar_year());
				}
			}
/// frac_nextyear used for tree init procedure
			if(sum_nextyear > 0.0) {

				// rescale active fractions so sum is 1.0 (or land cover fraction value)

				double ratio = sum_nextyear;

				if(!st_data[lc].NormalisedData())
					ratio /= gridcell.landcover.frac_nextyear[lc];

				stlist.firstobj();
				while(stlist.isobj) {
					StandType& st = stlist.getobj();
					Gridcellst& gcst = gridcell.st[st.id];
					if(st.landcover == lc) {
//						gcst.frac_old_orig = gcst.frac;
						gcst.frac_nextyear /= ratio;
					}
					stlist.nextobj();
				}

				if((ratio < 0.99 || ratio > 1.01) && printyear) {	// warn if sum is significantly different from 1.0 
					log_info("WARNING ! %s st fraction sum is %7.5f for input year %d while LC fraction is %7.5f", 
						(char*)lcnames[lc], sum_nextyear, (year + 1), gridcell.landcover.frac_nextyear[lc]);
					log_info("Rescaling %s fractions for next year %d !", lcnames[lc], (date.get_calendar_year() + 1));
				}
			}
		}
	}

	// Update landcover frac_change arrays
	for(int i=0; i<NLANDCOVERTYPES; i++) {
		lc.frac_change[i] = lc.frac[i] - lc.frac_old[i];
	}

	// Convert fractions from landcover-based to gridcell-based and update gridcellst frac_change arrays
	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		Gridcellst& gcst = gridcell.st[st.id];

		if(frac_fixed[st.landcover] || st_data[st.landcover].NormalisedData())
			gcst.frac = gcst.frac * lc.frac[st.landcover];
		// Ignore very small fraction changes unless frac_old is 0 and frac larger than limit or if frac smaller than limit.
		if(fabs(gcst.frac_old - gcst.frac) < INPUT_RESOLUTION * 0.1 && !(!gcst.frac_old && gcst.frac > INPUT_RESOLUTION)
				&& !(gcst.frac < INPUT_RESOLUTION))
			gcst.frac = gcst.frac_old;
		if(gcst.frac < INPUT_RESOLUTION)
			gcst.frac = 0.0; 
		gcst.frac_change = gcst.frac - gcst.frac_old;

/// frac_nextyear used for tree init procedure
		if(frac_fixed[st.landcover] || st_data[st.landcover].NormalisedData())
			gcst.frac_nextyear = gcst.frac_nextyear * lc.frac_nextyear[st.landcover];
		// Ignore very small fraction changes unless frac_old is 0 and frac larger than limit or if frac smaller than limit.
		if(fabs(gcst.frac - gcst.frac_nextyear) < INPUT_RESOLUTION * 0.1 && !(!gcst.frac && gcst.frac_nextyear > INPUT_RESOLUTION)
				&& !(gcst.frac_nextyear < INPUT_RESOLUTION))
			gcst.frac_nextyear = gcst.frac;
		if(gcst.frac_nextyear < INPUT_RESOLUTION)
			gcst.frac_nextyear = 0.0; 

		stlist.nextobj();
	}

	if (peatduration) {
		if (date.year == 0) {	// Fixed through out time // DW_COMMENT need to work with restarts!!!
			int year_fixed = Peatdata.GetFirstyear();
			lc.wet_startday = Peatdata.Get(year_fixed, "Startday");
			lc.wet_days = Peatdata.Get(year_fixed, "Days");
		}
	}
}

double LandcoverInput::get_crop_fractions(Gridcell& gridcell, int year, TimeDataD& CFTdata, double sum_tot) {

	if(!run[CROPLAND] || !gridcell.landcover.frac[CROPLAND] || frac_fixed[CROPLAND] && !frac_fixed_default_crops)
		return 0.0;

	bool printyear = year >= CFTdata.GetFirstyear() && CFTdata.GetFirstyear() >= 0;
	double sum=0.0;
	Landcover& lc = gridcell.landcover;

	if(!(frac_fixed[CROPLAND] && frac_fixed_default_crops)) {

		/* sum fractions for active crop pft:s and discard unreasonable values
		 * if crop fraction sum is 0 this year, first try last year's values, then try the following years
		 */
		int first_data_year = CFTdata.GetFirstyear();
		int last_data_year = first_data_year + CFTdata.GetnYears() - 1;
		// first_data_year is -1 for static data
		if(first_data_year == -1) {
			first_data_year = year;
			last_data_year = year;
		}
		if(first_data_year != -1 && year == last_data_year + 1)
			log_info("Last year of cropland fraction data used from year %d and onwards", year);

		for(int y=year;y<=max(year,last_data_year) + 1;y++) {

			for(int i=0; i<nst; i++) {
				if(stlist[i].landcover == CROPLAND)	{

					double cropfrac = CFTdata.Get(y,stlist[i].name);

					if(cropfrac == NOTFOUND) {	// crop not found in input file
						cropfrac = 0.0;
					}
					else if(cropfrac < 0.0) {	// correct unreasonable values
						if(printyear && gridcell.landcover.frac[CROPLAND])
							log_info("WARNING ! crop fraction size out of limits, set to 0.0");
						cropfrac = 0.0;
					}
					else if(cropfrac > 1.01)	{	// correct unreasonable values
						if(printyear)
							log_info("WARNING ! crop fraction size out of limits, set to 1.0");
						cropfrac = 1.0;
					}
					if(!st_data[CROPLAND].NormalisedData())
						cropfrac /= sum_tot;
					sum += gridcell.st[i].frac = cropfrac;
				}
			}
			if(sum) {
				if(y != year && (printyear || !date.year)) {
					log_info("WARNING ! crop fraction sum is 0.0 for year %d while LU[CROPLAND] is > 0 !", year);
					log_info("Using values for year %d.", y);
				}
				break;
			}
			else if(y == year) {
				// If no crop values for this year, first try to use last year's values
				for(int i=0; i<nst; i++) {
					if(stlist[i].landcover == CROPLAND && gridcell.landcover.frac_old[CROPLAND])
						sum += gridcell.st[i].frac = gridcell.st[i].frac_old_orig;
				}
				if(sum) {
					if(printyear) {
						log_info("WARNING ! crop fraction sum is 0.0 for year %d while LU[CROPLAND] is > 0 !", year);
						log_info("Using previous values.");
					}
					break;
				}
			}
		}
		if(printyear && !sum) {
			log_info("WARNING ! crop fraction sum is 0.0 for year %d while LU[CROPLAND] is > 0 !", year);
		}
	}

	// If crop fraction data are missing or if fixed default crops are used, try to find suitable stand types to use.
	if(sum==0.0) {

		if(lc.frac[CROPLAND]) {
			// Use equal areas of rainfed stand types with tropical or temperate crop pft:s based on base temperatures
			int nsts = 0;
			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();
				Gridcellst& gcst = gridcell.st[st.id];
				if(st.landcover == CROPLAND) {
					if(st.rotation.nmanagements == 1 && st.get_management(0).hydrology == RAINFED) {
						if(gridcell.get_lat() > 30 || gridcell.get_lat() < -30) {
							if(pftlist[pftlist.getpftid(st.get_management(0).pftname)].tb <= 5) {
								gcst.frac = 1.0;
								nsts++;
							}
						}
						else {
							if(pftlist[pftlist.getpftid(st.get_management(0).pftname)].tb > 5) {
								gcst.frac = 1.0;
								nsts++;
							}
						}
					}
				}
				stlist.nextobj();
			}
			if(nsts) {
				if(lc.frac[CROPLAND] && (printyear || !date.year) && !(frac_fixed[CROPLAND] && frac_fixed_default_crops))
					log_info("Missing crop fraction data. Using available suitable stand types for year %d", year);
			}
			else {
				fail("No suitable default stand types available for filling missing data. Check crop fraction input data\n\n");
			}
			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();
				Gridcellst& gcst = gridcell.st[st.id];
				if(st.landcover == CROPLAND) {
					gcst.frac /= nsts;
					gcst.frac_old_orig = gcst.frac;
				}
				stlist.nextobj();
			}
		}
	}

	return sum;
}

bool LandcoverInput::get_land_transitions(Gridcell& gridcell, double* st_frac_transfer) {

	bool result = false;
	if(gross_land_transfer == ST_GROSS_INPUT) {

		/* Read stand type transfer fractions from file here and put them into the st_frac_transfer array.
		 * Landcover and stand type net fractions still need to be read from file as previously.
		 * return get_st_transfer(gridcell);
		 */
		result = get_st_transfer(gridcell, gridcell.landcover.frac_transfer, st_frac_transfer);

		// Currently no code for gross stand type transfer option
	}
	else if(gross_land_transfer == LC_GROSS_INPUT) {

		/* Read landcover transfer fractions from file here and put them into the st_frac_transfer array.
		 * Landcover and stand type net fractions still need to be read from file as previously.
		 */
		result = get_lc_transfer(gridcell);
	}
	return result;
}

/// Read ST transitions from file
/** This implementation uses stand type area fraction transfer values in columns with names of the transfers according to
 *	[from_st]_[to_st] in the header line. The calendar year column value in this file needs to be the year preceding the
 *	updated area fraction value in the net LC and ST fraction files. No warnings are issued for missing columns for the ST transfer combinations.
 *  The consistency between gross ST transfer input and net LC and ST files is the responsibility of the end user and is not checked by this function.
 */
bool LandcoverInput::get_st_transfer(Gridcell& gridcell, double lc_frac_transfer[][NLANDCOVERTYPES], double* st_frac_transfer) {

	double tot_frac_ch = 0.0;

	if(!grossSTC.isloaded() || date.get_calendar_year() < getfirsthistyear() + 1)
		return false;

	// Before second year of net land cover input gross_lc_change_frac must be zero.

	int year = date.get_calendar_year() - 1;

	for(int from=0; from<nst; from++) {
		for(int to=0; to<nst; to++) {

			char st_transfer_column[100] = { '\0'};
			sprintf(st_transfer_column, "%s_%s", (char*)stlist[from].name, (char*)stlist[to].name);

			double frac_transfer;
			frac_transfer = grossSTC.Get(year, st_transfer_column, true);
			if(frac_transfer == NOTFOUND)
				frac_transfer = 0.0;
			st_frac_transfer[index(from, to)] += frac_transfer;
			lc_frac_transfer[stlist[from].landcover][stlist[to].landcover] += frac_transfer;
		}
	}

	if(tot_frac_ch > 1.0e-14)
		return true;
	else
		return false;
}

/// Read LUC transitions from file
/** This implementation uses LUH2 gross transfer standard text format input files (see guess.doc). Transitions to and from  
 *  cropland (c), pasture (p), barren land (b), urban land (u), primary PNV (v) and secondary PNV (s) are represented as
 *  columns with header names, e.g. "cp" denoting.cropland to pasture transfers.
 *  The function call to adjust_gross_transfers() tries to solve inconsistencies between net landcover and gross land
 *	transfer input and adjusts the gross transfer two-dimensional array, keeping the net landcover change array intact.
 */
bool LandcoverInput::get_lc_transfer(Gridcell& gridcell) {

	double tot_frac_change = 0.0;
	Landcover& lc = gridcell.landcover;

	if(!grossLUC.isloaded() || date.get_calendar_year() < getfirsthistyear() + 1)
		return false;

	// Before second year of net land cover input gross_lc_change_frac must be zero.

	int year = date.get_calendar_year() - 1;

	/* Assume that transitions in file are correct at end of year, therefore want to get 
	 * "last year's" transitions, as landcover_dynamics is called at the beginning of the year.
	 * Transfers from primary (v) and secondary (s) land preferentially reduces the oldest and the 
	 * youngest stands, respectively. Transitions from primary to secondary NATURAL land result 
	 * in killing of vegetation and creating a new NATURAL stand.
	 */

	// Options to not include barren and urban transfers in the input file (avoiding warnings if absent in file)
	const bool use_barren_transfers = false;
	const bool use_urban_transfers = true;
	const bool use_peatland_transfers = true;
	double frac_transfer;

	lc.frac_transfer[CROPLAND][PASTURE] += (frac_transfer = grossLUC.Get(year,"cp")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[PASTURE][CROPLAND] += (frac_transfer = grossLUC.Get(year,"pc")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[PASTURE][NATURAL] += (frac_transfer = grossLUC.Get(year,"pv")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[NATURAL][PASTURE] += (frac_transfer = grossLUC.Get(year,"vp")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[NATURAL][CROPLAND] += (frac_transfer = grossLUC.Get(year,"vc")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[CROPLAND][NATURAL] += (frac_transfer = grossLUC.Get(year,"cv")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[NATURAL][CROPLAND] += (frac_transfer = grossLUC.Get(year,"sc")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[CROPLAND][NATURAL] += (frac_transfer = grossLUC.Get(year,"cs")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[NATURAL][PASTURE] += (frac_transfer = grossLUC.Get(year,"sp")) != NOTFOUND ? frac_transfer : 0.0;
	lc.frac_transfer[PASTURE][NATURAL] += (frac_transfer = grossLUC.Get(year,"ps")) != NOTFOUND ? frac_transfer : 0.0;

	if(use_barren_transfers) {
		lc.frac_transfer[BARREN][CROPLAND] += (frac_transfer = grossLUC.Get(year,"bc")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[CROPLAND][BARREN] += (frac_transfer = grossLUC.Get(year,"cb")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[BARREN][PASTURE] += (frac_transfer = grossLUC.Get(year,"bp")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[PASTURE][BARREN] += (frac_transfer = grossLUC.Get(year,"pb")) != NOTFOUND ? frac_transfer : 0.0;

		lc.frac_transfer[BARREN][NATURAL] += (frac_transfer = grossLUC.Get(year,"bs")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][BARREN] += (frac_transfer = grossLUC.Get(year,"sb")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][BARREN] += (frac_transfer = grossLUC.Get(year,"vb")) != NOTFOUND ? frac_transfer : 0.0;
	}
	if(use_urban_transfers && run[URBAN]) {
		lc.frac_transfer[URBAN][CROPLAND] += (frac_transfer = grossLUC.Get(year,"uc")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[CROPLAND][URBAN] += (frac_transfer = grossLUC.Get(year,"cu")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[URBAN][PASTURE] += (frac_transfer = grossLUC.Get(year,"up")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[PASTURE][URBAN] += (frac_transfer = grossLUC.Get(year,"pu")) != NOTFOUND ? frac_transfer : 0.0;

		lc.frac_transfer[URBAN][NATURAL] += (frac_transfer = grossLUC.Get(year,"us")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][URBAN] += (frac_transfer = grossLUC.Get(year,"su")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][URBAN] += (frac_transfer = grossLUC.Get(year,"vu")) != NOTFOUND ? frac_transfer : 0.0;
		if (use_barren_transfers) {
			// There should be transition to/from urban and barren in case of use_barren_transfers && use_urban_transfers
			// but at present the grossLUC file does not have these anyhow
			lc.frac_transfer[URBAN][BARREN] += (frac_transfer = grossLUC.Get(year,"ub")) != NOTFOUND ? frac_transfer : 0.0;
			lc.frac_transfer[BARREN][URBAN] += (frac_transfer = grossLUC.Get(year,"bu")) != NOTFOUND ? frac_transfer : 0.0;
		}
	}
	if(use_peatland_transfers && run[PEATLAND]) {
		/** Various potential wetland approaches exists
		 * 
		 * wetland approach E:
		 *   scale down avalable LU (1 - BARREN) and make room for time invariant PEATLAND
		 *   all gross scalled down accordingly, all from and to PEATLAND should be zero
		 *
		 * wetland approach F: Note: approach F failed DO NOT USE!
		 *   split the LUH2 NATURAL into NATURAL and PEATLAND(wetland)
		 *   All net and gross lu data were adjusted, based on the wetland/(wetland+natural) ratios accordingly.
		 *
		 * A potential existing wood harvest on secondary natural (see file_woodharv_frac) is currently not split up,
		 * since wetland has no tree PFTs, hence no wood harvest possible.
		 */
		lc.frac_transfer[PEATLAND][CROPLAND] += (frac_transfer = grossLUC.Get(year,"wc")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[CROPLAND][PEATLAND] += (frac_transfer = grossLUC.Get(year,"cw")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[PEATLAND][PASTURE] += (frac_transfer = grossLUC.Get(year,"wp")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[PASTURE][PEATLAND] += (frac_transfer = grossLUC.Get(year,"pw")) != NOTFOUND ? frac_transfer : 0.0;

		lc.frac_transfer[PEATLAND][NATURAL] += (frac_transfer = grossLUC.Get(year,"ws")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][PEATLAND] += (frac_transfer = grossLUC.Get(year,"sw")) != NOTFOUND ? frac_transfer : 0.0;
		lc.frac_transfer[NATURAL][PEATLAND] += (frac_transfer = grossLUC.Get(year,"vw")) != NOTFOUND ? frac_transfer : 0.0;
		if (use_barren_transfers) {
			lc.frac_transfer[PEATLAND][BARREN] += (frac_transfer = grossLUC.Get(year,"wb")) != NOTFOUND ? frac_transfer : 0.0;
			lc.frac_transfer[BARREN][PEATLAND] += (frac_transfer = grossLUC.Get(year,"bw")) != NOTFOUND ? frac_transfer : 0.0;
		}
		if (use_urban_transfers && run[URBAN]) {
			lc.frac_transfer[PEATLAND][URBAN] += (frac_transfer = grossLUC.Get(year,"wu")) != NOTFOUND ? frac_transfer : 0.0;
			lc.frac_transfer[URBAN][PEATLAND] += (frac_transfer = grossLUC.Get(year,"uw")) != NOTFOUND ? frac_transfer : 0.0;
		}
	}

	// Distinguish between primary and secondary PNV in the input file (default true)
	if(ifprimary_lc_transfer) {
		lc.forest_lc_subset_transfer.primary[NATURAL][PASTURE] += (frac_transfer = grossLUC.Get(year,"vp")) != NOTFOUND ? frac_transfer : 0.0;
		lc.forest_lc_subset_transfer.primary[NATURAL][CROPLAND] += (frac_transfer = grossLUC.Get(year,"vc")) != NOTFOUND ? frac_transfer : 0.0;
		if(use_barren_transfers)
			lc.forest_lc_subset_transfer.primary[NATURAL][BARREN] += (frac_transfer = grossLUC.Get(year,"vb")) != NOTFOUND ? frac_transfer : 0.0;
		if(use_urban_transfers)
			lc.forest_lc_subset_transfer.primary[NATURAL][URBAN] += (frac_transfer = grossLUC.Get(year,"vu")) != NOTFOUND ? frac_transfer : 0.0;
		if(use_peatland_transfers)
			lc.forest_lc_subset_transfer.primary[NATURAL][PEATLAND] += (frac_transfer = grossLUC.Get(year, "vw")) != NOTFOUND ? frac_transfer : 0.0;

		// Use transitions from virgin to secondary natural land (default false).
		if(ifprimary_to_secondary_transfer) {
			lc.frac_transfer[NATURAL][NATURAL] += (frac_transfer = grossLUC.Get(year,"vs")) != NOTFOUND ? frac_transfer : 0.0;
			lc.forest_lc_subset_transfer.primary[NATURAL][NATURAL] += (frac_transfer = grossLUC.Get(year,"vs")) != NOTFOUND ? frac_transfer : 0.0;
		}
		use_primary_lc_transfer = true;
	}

	for(int from=0; from<NLANDCOVERTYPES; from++) {
		if(run[from]) {
			for(int to=0; to<NLANDCOVERTYPES; to++) {
				if(run[to])					
					tot_frac_change += lc.frac_transfer[to][from];
			}
		}
	}

	// Check if gross lcc input data are consistent with net lcc input file. Try to adjust if not.
	adjust_gross_transfers(gridcell, lc.frac_change, lc.frac_transfer, lc.forest_lc_subset_transfer, tot_frac_change);

	if(largerthanzero(tot_frac_change, -14))
		return true;
	else
		return false;
}

/// Help function for get_lc_transfer() to adjust inconsistencies between net land cover inout and gross land cover transitions.
void adjust_gross_transfers(Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES], 
		forest_lc_frac_transfer& forest_lc_subset_transfer, double& tot_frac_change) {

	const bool print_adjustment_info = false;
	bool error = false;
	double net_lc_change[NLANDCOVERTYPES] = {0.0};
	double gross_lc_increase[NLANDCOVERTYPES] = {0.0};
	double gross_lc_decrease[NLANDCOVERTYPES] = {0.0};
	double pos_error = 0.0;
	double neg_error = 0.0;

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		if(run[from]) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if(run[to]) {
					
					net_lc_change[from] -= lc_frac_transfer[from][to];
					net_lc_change[from] += lc_frac_transfer[to][from];
					gross_lc_decrease[from] += lc_frac_transfer[from][to];
					gross_lc_increase[from] += lc_frac_transfer[to][from];
				}
			}
			if(fabs(landcoverfrac_change[from] - net_lc_change[from])  > 1.0e-15) {
				error = true;
				if(print_adjustment_info) {
					log_info("In get_lc_transfer: lc_change_array sum not equal to landcoverfrac_change value for landcover %d", 
						     from);
					log_info("dif=%.15f", net_lc_change[from] - landcoverfrac_change[from]);
					log_info("lc_change_array sum=%.15f", net_lc_change[from]);
					log_info("landcoverfrac_change=%.15f", landcoverfrac_change[from]);
				}
			}

		}
	}

	// Save forest class percentages before correcting transitions

	double prim_ratio[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	double sec_young_ratio[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			if(lc_frac_transfer[from][to]) {
				prim_ratio[from][to] = forest_lc_subset_transfer.primary[from][to] / lc_frac_transfer[from][to];
				sec_young_ratio[from][to] = forest_lc_subset_transfer.secondary_young[from][to] / lc_frac_transfer[from][to];
			}
		}
	}

	// Try to balance overshoot; only existing transitions are altered.

	// 1
	for(int first=0; first<NLANDCOVERTYPES; first++) {

		double twoway_overshoot = min(gross_lc_increase[first] - gridcell.landcover.frac[first], gross_lc_decrease[first] 
			- gridcell.landcover.frac_old[first]);

		if(twoway_overshoot > 1.0e-15) {

			for(int second=0; second<NLANDCOVERTYPES; second++) {

				for(int third=0; third<NLANDCOVERTYPES; third++) {

					if(lc_frac_transfer[first][second] >= twoway_overshoot 
						&& lc_frac_transfer[third][second] >= twoway_overshoot
						&& lc_frac_transfer[third][first] >= twoway_overshoot
						&& first != second && first != third) {

						if(print_adjustment_info) {
							log_info("Year %d: Trying to balance two-way overshoot %.18f of lc %d.", date.year, twoway_overshoot, first);
							log_info("lc_frac_transfer[%d][%d] before: %.15f", first, second, lc_frac_transfer[first][second]);
							log_info("lc_frac_transfer[%d][%d] before: %.15f", third, second, lc_frac_transfer[third][second]);
							log_info("lc_frac_transfer[%d][%d] before: %.15f", third, first, lc_frac_transfer[third][first]);
						}

						lc_frac_transfer[first][second] -= twoway_overshoot;
						lc_frac_transfer[third][second] += twoway_overshoot;
						lc_frac_transfer[third][first] -= twoway_overshoot;
						gross_lc_decrease[first] -= twoway_overshoot;
						gross_lc_increase[second] -= twoway_overshoot;
						gross_lc_decrease[third] += twoway_overshoot;
						gross_lc_increase[second] += twoway_overshoot;
						gross_lc_decrease[third] -= twoway_overshoot;
						gross_lc_increase[first] -= twoway_overshoot;

						if(print_adjustment_info) {
							log_info("Year %d: After balancing lc %d.", date.year, first);
							log_info("lc_frac_transfer[%d][%d] after: %.15f", first, second, lc_frac_transfer[first][second]);
							log_info("lc_frac_transfer[%d][%d] after: %.15f", third, second, lc_frac_transfer[third][second]);
							log_info("lc_frac_transfer[%d][%d] after: %.15f", third, first, lc_frac_transfer[third][first]);
						}
					}
				}
			}
		}
	}

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		if(print_adjustment_info) {
			double balance = gross_lc_decrease[from] - gridcell.landcover.frac_old[from];
			if(balance > 1.0e-15)
				log_info("Year %d: remaining undershoot %.18f of lc %d.", date.year, balance, from);
			balance = gross_lc_increase[from] - gridcell.landcover.frac[from];
			if(balance > 1.0e-15)
				log_info("Year %d: remaining overshoot %.18f of lc %d.", date.year, balance, from);				}
	}


	// 2
	for(int from=0; from<NLANDCOVERTYPES; from++) {

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			if(lc_frac_transfer[from][to]) {

				if(gross_lc_increase[to] > gridcell.landcover.frac[to]) {

					double balance = gross_lc_increase[to] - gridcell.landcover.frac[to];
					if(print_adjustment_info) {
						log_info("Year %d: Trying to balance overshoot %.18f of lc %d.", date.year, balance, to);
						log_info("lc_frac_transfer[%d][%d] before: %.15f", from, to, lc_frac_transfer[from][to]);
						log_info("lc_frac_transfer[%d][%d] before: %.15f", to, from, lc_frac_transfer[to][from]);
					}
					balance = min(balance, lc_frac_transfer[from][to]);
					balance = min(balance, lc_frac_transfer[to][from]);
					lc_frac_transfer[from][to] -= balance;
					gross_lc_decrease[from] -= balance;
					gross_lc_increase[from] -= balance;
					if(from != to) {
						lc_frac_transfer[to][from] -= balance;
						gross_lc_decrease[to] -= balance;
						gross_lc_increase[to] -= balance;
					}
					if(print_adjustment_info) {
						log_info("lc_frac_transfer[%d][%d] after: %.15f", from, to, lc_frac_transfer[from][to]);
						log_info("lc_frac_transfer[%d][%d] after: %.15f", to, from, lc_frac_transfer[to][from]);
					}
				}
				if(gridcell.landcover.frac_old[from] - gross_lc_decrease[from] < 0.0) {

					double balance = gross_lc_decrease[from] - gridcell.landcover.frac_old[from];
					if(print_adjustment_info) {
						log_info("Year %d: Trying to balance overshoot %.18f of lc %d.", date.year, balance, from);
						log_info("lc_frac_transfer[%d][%d] before: %.15f", from, to, lc_frac_transfer[from][to]);
						log_info("lc_frac_transfer[%d][%d] before: %.15f", to, from, lc_frac_transfer[to][from]);
					}
					balance = min(balance, lc_frac_transfer[from][to]);
					balance = min(balance, lc_frac_transfer[to][from]);
					lc_frac_transfer[from][to] -= balance;
					gross_lc_decrease[from] -= balance;
					gross_lc_increase[from] -= balance;
					if(from != to) {
						lc_frac_transfer[to][from] -= balance;
						gross_lc_decrease[to] -= balance;
						gross_lc_increase[to] -= balance;
					}
					if(print_adjustment_info) {
						log_info("lc_frac_transfer[%d][%d] after: %.15f", from, to, lc_frac_transfer[from][to]);
						log_info("lc_frac_transfer[%d][%d] after: %.15f", to, from, lc_frac_transfer[to][from]);
					}
				}
			}
		}
	}

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		if(print_adjustment_info) {
			double balance = gross_lc_decrease[from] - gridcell.landcover.frac_old[from];
			if(balance > 1.0e-15)
				log_info("Year %d: remaining undershoot %.18f of lc %d.", date.year, balance, from);
			balance = gross_lc_increase[from] - gridcell.landcover.frac[from];
			if(balance > 1.0e-15)
				log_info("Year %d: remaining overshoot %.18f of lc %d.", date.year, balance, from);				}
	}

	// Discard obvious artefacts
	for(int from=0; from<NLANDCOVERTYPES; from++) {
		for(int to=0; to<NLANDCOVERTYPES; to++) {
			if(lc_frac_transfer[from][to]) {
				if(!gridcell.landcover.frac_old[from] || !gridcell.landcover.frac[to]) {
					if(print_adjustment_info) {
						log_info("Year %d: Rejecting transfer %.18f from lc %d to %d", 
							date.year, lc_frac_transfer[from][to], from, to);
						log_info("frac_old[%d]=%.15f, frac[%d]=%.15f", 
							from, gridcell.landcover.frac_old[from], to, gridcell.landcover.frac[to]);
					}
					gross_lc_decrease[from] -= lc_frac_transfer[from][to];
					gross_lc_increase[to] -= lc_frac_transfer[from][to];
					net_lc_change[from] += lc_frac_transfer[from][to];
					net_lc_change[to] -= lc_frac_transfer[from][to];
					tot_frac_change -= lc_frac_transfer[from][to];
					lc_frac_transfer[from][to] = 0.0;
				}
			}
		}
	}

	double frac_remain[NLANDCOVERTYPES];
	double frac_remain_old[NLANDCOVERTYPES];
	for(int lc=0; lc<NLANDCOVERTYPES; lc++) {
		frac_remain[lc] = gridcell.landcover.frac[lc];
		frac_remain_old[lc] = gridcell.landcover.frac_old[lc];
	}
	for(int from=0; from<NLANDCOVERTYPES; from++) {
		for(int to=0; to<NLANDCOVERTYPES; to++) {
			if(lc_frac_transfer[from][to]) {
				if(lc_frac_transfer[from][to] > frac_remain_old[from] || lc_frac_transfer[from][to] > frac_remain[to]) {
					if(print_adjustment_info) {
						log_info("Year %d: Adjusting transfer %.18f from lc %d to %d", 
							date.year, lc_frac_transfer[from][to], from, to);
						log_info("frac_old[%d]=%.15f, frac[%d]=%.15f", 
							from, gridcell.landcover.frac_old[from], to, gridcell.landcover.frac[to]);
					}
					double transfer = min(lc_frac_transfer[from][to], min(frac_remain_old[from], frac_remain[to]));
					gross_lc_decrease[from] -= lc_frac_transfer[from][to] - transfer;
					gross_lc_increase[to] -= lc_frac_transfer[from][to] - transfer;
					net_lc_change[from] += lc_frac_transfer[from][to] - transfer;
					net_lc_change[to] -= lc_frac_transfer[from][to] - transfer;
					tot_frac_change -= lc_frac_transfer[from][to] - transfer;
					lc_frac_transfer[from][to] = transfer;
					frac_remain_old[from] -= transfer;
					frac_remain[to] -= transfer;
				}
			}
		}
	}

	// Adjust transitions for lc:s with net changes deviating from the landcover fractions

	if(error) {

		double partition_adjustment[NLANDCOVERTYPES][NLANDCOVERTYPES];
		double original_error[NLANDCOVERTYPES] = {0.0};

		for(int from=0; from<NLANDCOVERTYPES; from++) {
			for(int to=0; to<NLANDCOVERTYPES; to++) {
				partition_adjustment[from][to] = 0.0;
			}
		}

		// Determine how much to change transfer in each direction

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if((lc_frac_transfer[from][to] + lc_frac_transfer[to][from]) > 0.0)
					partition_adjustment[from][to] = lc_frac_transfer[from][to] / (lc_frac_transfer[from][to] + lc_frac_transfer[to][from]);
			}

			if(run[from])
				original_error[from] = net_lc_change[from] - landcoverfrac_change[from];

			if(original_error[from] > 0.0)
				pos_error += original_error[from];
			else if(original_error[from] < 0.0)
				neg_error += original_error[from];
		}

		if(fabs(pos_error + neg_error) > 1.0e-15)
			fail("\nYear %d: pos_error + neg_error = %.15f\n\n", date.year, pos_error + neg_error);


		// 1. Only lc:s with existing opposing errors are altered.

		double residual_error[NLANDCOVERTYPES] = {0.0};
		for(int from=0; from<NLANDCOVERTYPES; from++)
			residual_error[from] = original_error[from];

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(run[from] && fabs(residual_error[from]) > 1.0e-15) {

				if(print_adjustment_info)
					log_info("residual_error[%d] before = %.15f", from, residual_error[from]);

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(fabs(residual_error[from]) > 1.0e-15 && run[to] && fabs(residual_error[to])  > 1.0e-15 && from != to) {

						// Errors must have opposite signs
						if(fabs(residual_error[from] + residual_error[to]) - (fabs(residual_error[from]) + fabs(residual_error[to])) < -1.0e-15) {

							if(print_adjustment_info)
								log_info("Trying with lc %d and %d", from, to);
							// Correct transfer between two lc:s
							if((lc_frac_transfer[from][to] + lc_frac_transfer[to][from]) > 0.0) {

								if(print_adjustment_info) {
									log_info("Before: transfer lc %d to %d: %.15f", from, to, lc_frac_transfer[from][to]);
									log_info("Before: transfer lc %d to %d: %.15f", to, from, lc_frac_transfer[to][from]);
								}
								double effective_corr;

								if(residual_error[from] >= 0.0) {
									effective_corr = min(residual_error[from], fabs(residual_error[to]));
									// Make sure lc_frac_transfer[to][from] not negative
									if(partition_adjustment[to][from])
										effective_corr = min(effective_corr, lc_frac_transfer[to][from] / partition_adjustment[to][from]);
									// Make sure resulting frac[to] not overshot and frac[from] not depleted
									if(partition_adjustment[from][to]) {
										effective_corr = min(effective_corr, max(0.0, gridcell.landcover.frac[to] - gross_lc_increase[to]) 
											/ partition_adjustment[from][to]);
										effective_corr = min(effective_corr, max(0.0, gridcell.landcover.frac_old[from] - gross_lc_decrease[from]) 
											/ partition_adjustment[from][to]);
									}
								}
								else {
									effective_corr = max(residual_error[from], -residual_error[to]);
									// Make sure lc_frac_transfer[from][to] not negative
									if(partition_adjustment[from][to])
										effective_corr = max(effective_corr, -lc_frac_transfer[from][to] / partition_adjustment[from][to]);
									// Make sure resulting frac[from] not overshot and frac[to] not depleted
									if(partition_adjustment[to][from]) {
										effective_corr = max(effective_corr, -max(0.0, gridcell.landcover.frac[from] - gross_lc_increase[from]) 
											/ partition_adjustment[to][from]);
										effective_corr = max(effective_corr, -max(0.0, gridcell.landcover.frac_old[to] - gross_lc_decrease[to]) 
											/ partition_adjustment[to][from]);
									}
								}

								lc_frac_transfer[from][to] += effective_corr * partition_adjustment[from][to];
								lc_frac_transfer[to][from] -= effective_corr * partition_adjustment[to][from];
								residual_error[from] -= effective_corr;
								residual_error[to] += effective_corr;
								gross_lc_decrease[from] += effective_corr * partition_adjustment[from][to];
								gross_lc_increase[from] -= effective_corr * partition_adjustment[to][from];
								gross_lc_decrease[to] -= effective_corr * partition_adjustment[to][from];
								gross_lc_increase[to] += effective_corr * partition_adjustment[from][to];
								if(print_adjustment_info) {
									log_info("After: transfer lc %d to %d: %.15f", from, to, lc_frac_transfer[from][to]);
									log_info("After: transfer lc %d to %d: %.15f", to, from, lc_frac_transfer[to][from]);
								}
							}
						}
					}
				}
				if(print_adjustment_info)
					log_info("residual_error[%d] after = %.15f", from, residual_error[from]);
			}
		}

		// 2. Using third lc

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(run[from] && fabs(residual_error[from]) > 1.0e-15) {

				if(print_adjustment_info)
					log_info("residual_error[%d] before = %.15f", from, residual_error[from]);

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(fabs(residual_error[from]) > 1.0e-15 && run[to] && fabs(residual_error[to]) > 1.0e-15 && from != to) {

						// Errors must have opposite signs
						if(fabs(residual_error[from] + residual_error[to]) - (fabs(residual_error[from]) + fabs(residual_error[to])) < -1.0e-15) {

							if(print_adjustment_info)
								log_info("Trying with lc %d and %d", from, to);

							if(print_adjustment_info) 
								log_info("Using third land cover type");
							for(int third=0; third<NLANDCOVERTYPES; third++) {

								if(lc_frac_transfer[from][third] + lc_frac_transfer[third][from] > 0.0 
									&& lc_frac_transfer[to][third] + lc_frac_transfer[third][to] && third != from && third != to) {

									if(print_adjustment_info) {
										log_info("Before: transfer lc %d to %d: %.15f", from, third, lc_frac_transfer[from][third]);
										log_info("Before: transfer lc %d to %d: %.15f", third, from, lc_frac_transfer[third][from]);
										log_info("Before: transfer lc %d to %d: %.15f", to, third, lc_frac_transfer[to][third]);
										log_info("Before: transfer lc %d to %d: %.15f", third, to, lc_frac_transfer[third][to]);
									}

									double effective_corr;

									if(residual_error[from] >= 0.0) {
										effective_corr = min(residual_error[from], fabs(residual_error[to]));
										// Make sure lc_frac_transfer[third][from] not negative
										if(partition_adjustment[third][from]) {
											effective_corr = min(effective_corr, lc_frac_transfer[third][from]
												/ partition_adjustment[third][from]);
										}
										// Make sure resulting frac[third] not overshot and frac[from] not depleted
										if(partition_adjustment[from][third]) {
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac[third] - (gross_lc_increase[third] + effective_corr * partition_adjustment[to][third])) 
													/ partition_adjustment[from][third]);
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac_old[from] - gross_lc_decrease[from]) 
													/ partition_adjustment[from][third]);
										}
										// Make sure lc_frac_transfer[to][third] not negative
										if(partition_adjustment[to][third])
											effective_corr = min(effective_corr, lc_frac_transfer[to][third] / partition_adjustment[to][third]);
										// Make sure resulting frac[to] not overshot and frac[third] not depleted
										if(partition_adjustment[third][to]) {
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac[to] - gross_lc_increase[to]) 
													/ partition_adjustment[third][to]);
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac_old[third] - (gross_lc_decrease[third] + effective_corr * partition_adjustment[third][from])) 
													/ partition_adjustment[third][to]);
										}
									}
									else {
										effective_corr = max(residual_error[from], -residual_error[to]);
										// Make sure lc_frac_transfer[from][third] not negative
										if(partition_adjustment[from][third])
											effective_corr = max(effective_corr, -lc_frac_transfer[from][third] / partition_adjustment[from][third]);
										// Make sure resulting frac[from] not overshot and frac[third] not depleted
										if(partition_adjustment[third][from]) {
											effective_corr = max(effective_corr, -max(0.0, gridcell.landcover.frac[from] - gross_lc_increase[from]) 
												/ partition_adjustment[third][from]);
											effective_corr = 
												max(effective_corr, -max(0.0, gridcell.landcover.frac_old[third] - (gross_lc_decrease[third] - effective_corr * partition_adjustment[third][to])) 
													/ partition_adjustment[third][from]);
										}
										// Make sure lc_frac_transfer[third][to] not negative
										if(partition_adjustment[third][to])
											effective_corr = max(effective_corr, -lc_frac_transfer[third][to] / partition_adjustment[third][to]);
										// Make sure resulting frac[third] not overshot and frac[to] not depleted
										if(partition_adjustment[to][third]) {
											effective_corr = 
												max(effective_corr, -max(0.0, gridcell.landcover.frac[third] - (gross_lc_increase[third] - effective_corr * partition_adjustment[from][third])) 
													/ partition_adjustment[to][third]);
											effective_corr = 
												max(effective_corr, -max(0.0, gridcell.landcover.frac_old[to] - gross_lc_decrease[to]) 
													/ partition_adjustment[to][third]);
										}
									}

									lc_frac_transfer[from][third] += effective_corr * partition_adjustment[from][third];
									lc_frac_transfer[third][from] -= effective_corr * partition_adjustment[third][from];
									residual_error[from] -= effective_corr;
									residual_error[third] += effective_corr;

									lc_frac_transfer[to][third] -= effective_corr * partition_adjustment[to][third];
									lc_frac_transfer[third][to] += effective_corr * partition_adjustment[third][to];
									residual_error[to] += effective_corr;
									residual_error[third] -= effective_corr;

									gross_lc_decrease[from] += effective_corr * partition_adjustment[from][third];
									gross_lc_increase[from] -= effective_corr * partition_adjustment[third][from];
									gross_lc_decrease[third] -= effective_corr * partition_adjustment[third][from];
									gross_lc_increase[third] += effective_corr * partition_adjustment[from][third];

									gross_lc_decrease[to] -= effective_corr * partition_adjustment[to][third];
									gross_lc_increase[to] += effective_corr * partition_adjustment[third][to];
									gross_lc_decrease[third] += effective_corr * partition_adjustment[third][to];
									gross_lc_increase[third] -= effective_corr * partition_adjustment[to][third];

									if(print_adjustment_info) {
										log_info("After: transfer lc %d to %d: %.15f", from, third, lc_frac_transfer[from][third]);
										log_info("After: transfer lc %d to %d: %.15f", third, from, lc_frac_transfer[third][from]);
										log_info("After: transfer lc %d to %d: %.15f", to, third, lc_frac_transfer[to][third]);
										log_info("After: transfer lc %d to %d: %.15f", third, to, lc_frac_transfer[third][to]);
									}
								}
							}

						}
					}
				}
				if(print_adjustment_info)
					log_info("residual_error[%d] after = %.15f", from, residual_error[from]);
			}
		}

		// 3. New direct transfer

		if(print_adjustment_info)
			log_info("Dealing with rounding artefacts with new direct transfer");
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(run[from] && fabs(residual_error[from]) > 1.0e-15) {

				if(print_adjustment_info)
					log_info("residual_error[%d] before = %.15f", from, residual_error[from]);

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(fabs(residual_error[from]) > 1.0e-15 && run[to] && fabs(residual_error[to])  > 1.0e-15 && from != to) {

						// Errors must have opposite signs
						if(fabs(residual_error[from] + residual_error[to]) - (fabs(residual_error[from]) + fabs(residual_error[to])) < -1.0e-15) {

							if(print_adjustment_info) {
								log_info("Dealing with rounding artefacts with new direct transfer");
								log_info("Before: transfer lc %d to %d: %.15f", from, to, lc_frac_transfer[from][to]);
								log_info("Before: transfer lc %d to %d: %.15f", to, from, lc_frac_transfer[to][from]);
							}
							double effective_corr;

							if(residual_error[from] >= 0.0) {
								effective_corr = min(residual_error[from], fabs(residual_error[to]));
								effective_corr = min(effective_corr, max(0.0, gridcell.landcover.frac_old[from] - gross_lc_decrease[from]));
								effective_corr = min(effective_corr, max(0.0, gridcell.landcover.frac[to] - gross_lc_increase[to]));
							}
							else {
								effective_corr = max(residual_error[from], -residual_error[to]);
								effective_corr = max(effective_corr, -max(0.0, gridcell.landcover.frac_old[to] - gross_lc_decrease[to]));
								effective_corr = max(effective_corr, -max(0.0, gridcell.landcover.frac[from] - gross_lc_increase[from]));
							}
							if(effective_corr >= 0.0) {
								lc_frac_transfer[from][to] += effective_corr;
								gross_lc_decrease[from] += effective_corr;
								gross_lc_increase[to] += effective_corr;
							}
							else {
								lc_frac_transfer[to][from] -= effective_corr;
								gross_lc_decrease[to] -= effective_corr;
								gross_lc_increase[from] -= effective_corr;
							}
							residual_error[from] -= effective_corr;
							residual_error[to] += effective_corr;

							if(print_adjustment_info) {
								log_info("After: transfer lc %d to %d: %.15f", from, to, lc_frac_transfer[from][to]);
								log_info("After: transfer lc %d to %d: %.15f", to, from, lc_frac_transfer[to][from]);
							}
						}
					}
				}
				if(print_adjustment_info)
					log_info("residual_error[%d] after = %.15f", from, residual_error[from]);
			}
		}

		// 4. Third lc used, relaxed rules

		if(print_adjustment_info) 
			log_info("Using third land cover type, relaxed rules");
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(run[from] && fabs(residual_error[from]) > 1.0e-15) {

				if(print_adjustment_info)
					log_info("residual_error[%d] before = %.15f", from, residual_error[from]);

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(fabs(residual_error[from]) > 1.0e-15 && run[to] && fabs(residual_error[to])  > 1.0e-15 && from != to) {

						// Errors must have opposite signs
						if(fabs(residual_error[from] + residual_error[to]) - (fabs(residual_error[from]) + fabs(residual_error[to])) < -1.0e-15) {

							if(print_adjustment_info) 
								log_info("Using third land cover type, relaxed rules");
							if(print_adjustment_info)
								log_info("Trying again with lc %d and %d", from, to);

							for(int third=0; third<NLANDCOVERTYPES; third++) {

								if(partition_adjustment[from][third] + partition_adjustment[third][from] > 0.0 && third != from && third != to) {
									if(print_adjustment_info) {
										log_info("Before: transfer lc %d to %d: %.15f", from, third, lc_frac_transfer[from][third]);
										log_info("Before: transfer lc %d to %d: %.15f", third, from, lc_frac_transfer[third][from]);
										log_info("Before: transfer lc %d to %d: %.15f", to, third, lc_frac_transfer[to][third]);
										log_info("Before: transfer lc %d to %d: %.15f", third, to, lc_frac_transfer[third][to]);
									}
									double effective_corr;

									// Transfer between from and third: 
									if(residual_error[from] >= 0.0) {
										effective_corr = min(residual_error[from], fabs(residual_error[to]));
										// Make sure lc_frac_transfer[third][from] not negative
										if(partition_adjustment[third][from])
											effective_corr = min(effective_corr, lc_frac_transfer[third][from] / partition_adjustment[third][from]);
										// Make sure resulting frac[third] not overshot and frac[from] not depleted
										if(partition_adjustment[from][third]) {
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac[third] - gross_lc_increase[third] + effective_corr)
													/ partition_adjustment[from][third]);
											effective_corr = 
												min(effective_corr, max(0.0, gridcell.landcover.frac_old[from] - gross_lc_decrease[from]) 
													/ partition_adjustment[from][third]);
										}
										// Make sure resulting frac[to] not overshot and frac[third] not depleted
										effective_corr = 
											min(effective_corr, max(0.0, gridcell.landcover.frac[to] - gross_lc_increase[to]));
										effective_corr = 
											min(effective_corr, max(0.0, gridcell.landcover.frac_old[third] - gross_lc_decrease[third] + effective_corr));
									}
									else {
										effective_corr = max(residual_error[from], -residual_error[to]);
										// Make sure lc_frac_transfer[from][third] not negative
										if(partition_adjustment[from][third])
											effective_corr = max(effective_corr, -lc_frac_transfer[from][third] / partition_adjustment[from][third]);
										// Make sure resulting frac[from] not overshot and frac[third] not depleted
										if(partition_adjustment[third][from]) {
											effective_corr = 
												max(effective_corr, -max(0.0, gridcell.landcover.frac[from] - gross_lc_increase[from]) 
													/ partition_adjustment[third][from]);
											effective_corr = 
												max(effective_corr, -max(0.0, gridcell.landcover.frac_old[third] - gross_lc_decrease[third] - effective_corr) 
													/ partition_adjustment[third][from]);
										}
										// Balancing [to][third] transfer: make sure resulting frac[third] not overshot and frac[to] not depleted
										effective_corr = 
											max(effective_corr, -max(0.0, gridcell.landcover.frac[third] - gross_lc_increase[third] - effective_corr));
										effective_corr = 
											max(effective_corr, -max(0.0, gridcell.landcover.frac_old[to] - gross_lc_decrease[to]));
									}

									lc_frac_transfer[from][third] += effective_corr * partition_adjustment[from][third];
									lc_frac_transfer[third][from] -= effective_corr * partition_adjustment[third][from];
									residual_error[from] -= effective_corr;
									residual_error[third] += effective_corr;

									gross_lc_decrease[from] += effective_corr * partition_adjustment[from][third];
									gross_lc_increase[from] -= effective_corr * partition_adjustment[third][from];
									gross_lc_decrease[third] -= effective_corr * partition_adjustment[third][from];
									gross_lc_increase[third] += effective_corr * partition_adjustment[from][third];


									effective_corr = -effective_corr;
									if(effective_corr >= 0.0) {
										effective_corr = min(effective_corr, fabs(residual_error[third]));
									}
									else {
										effective_corr = max(effective_corr, -residual_error[third]);
									}

									if(residual_error[to] >= 0.0) {
										lc_frac_transfer[to][third] += effective_corr;
										gross_lc_decrease[to] += effective_corr;
										gross_lc_increase[third] += effective_corr;
									}
									else {
										lc_frac_transfer[third][to] -= effective_corr;
										gross_lc_decrease[third] -= effective_corr;
										gross_lc_increase[to] -= effective_corr;
									}

									residual_error[to] -= effective_corr;
									residual_error[third] += effective_corr;

									if(print_adjustment_info) {
										log_info("After: transfer lc %d to %d: %.15f", from, third, lc_frac_transfer[from][third]);
										log_info("After: transfer lc %d to %d: %.15f", third, from, lc_frac_transfer[third][from]);
										log_info("After: transfer lc %d to %d: %.15f", to, third, lc_frac_transfer[to][third]);
										log_info("After: transfer lc %d to %d: %.15f", third, to, lc_frac_transfer[third][to]);
									}
								}
							}
						}
					}
				}
				if(print_adjustment_info)
					log_info("residual_error[%d] after = %.15f", from, residual_error[from]);
			}
		}

		// Correcting negative transfer fraction

		if(print_adjustment_info)
			log_info("Correcting negative transfer fraction:");
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if(lc_frac_transfer[from][to] < 0.0) {
					if(lc_frac_transfer[from][to] > 1.0e-30) {	// Reset "Negative zeros" from equations above silently.
						log_info("Correcting negative transfer fraction:");
						log_info("Before: transfer lc %d to %d: %.20f", from, to, lc_frac_transfer[from][to]);
						log_info("Before: transfer lc %d to %d: %.20f", to, from, lc_frac_transfer[to][from]);
					}
					lc_frac_transfer[to][from] -= lc_frac_transfer[from][to];
					lc_frac_transfer[from][to] = 0.0;
					if(lc_frac_transfer[from][to] > 1.0e-30) {
						log_info("After: transfer lc %d to %d: %.20f", from, to, lc_frac_transfer[from][to]);
						log_info("After: transfer lc %d to %d: %.20f", to, from, lc_frac_transfer[to][from]);
					}
				}
			}
		}


		bool stop = false;
		for(int from=0; from<NLANDCOVERTYPES; from++) {
			if(print_adjustment_info && fabs(original_error[from]) > 1.0e-15) {
				log_info("original_error[%d] before = %.15f", from, original_error[from]);
				log_info("residual_error[%d] after = %.15f", from, residual_error[from]);
			}
			if(residual_error[from] > 1.0e-14)
				stop = true;
		}
		if(print_adjustment_info)
			log_info("");
		if(stop)
			fail("Failing to balance lc transitions");
	}
	// Adjust primary land fractions
	for(int from=0; from<NLANDCOVERTYPES; from++) {
		for(int to=0; to<NLANDCOVERTYPES; to++) {
			forest_lc_subset_transfer.primary[from][to] = prim_ratio[from][to] * lc_frac_transfer[from][to];
			forest_lc_subset_transfer.secondary_young[from][to] = sec_young_ratio[from][to] * lc_frac_transfer[from][to];
		}
	}
}

int LandcoverInput::getfirsthistyear() {

	return LUdata.GetFirstyear();
}

ManagementInput::ManagementInput() {

	targetfrac_pft_mt = NULL;
}

ManagementInput::~ManagementInput() {

	if(targetfrac_pft_mt) {
		delete[] targetfrac_pft_mt;
	}
}

void ManagementInput::init(xtring file_gridlist) {

	if(!run_landcover)
		return;

	ListArray_id<Coord> gridlist;
	read_gridlist(gridlist, (const char *)file_gridlist);

	if(run[CROPLAND]) {

		file_sdates=param["file_sdates"].str;
		if(file_sdates != "")	{
			if(!sdates.Open(file_sdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_sdates);
			readsowingdates = true;
		}

		file_hdates=param["file_hdates"].str;
		if(file_hdates != "")	{
			if(!hdates.Open(file_hdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_hdates);
			readharvestdates = true;
		}

		file_Nfert=param["file_Nfert"].str;
		if(	file_Nfert != "")	{
			if(!Nfert.Open(file_Nfert, gridlist))
				fail("initio: could not open %s for input",(char*)file_Nfert);
			readNfert = true;
		}
		if(param.find(xtring("file_NfertMan"))) {
			file_NfertMan=param["file_NfertMan"].str;
			if(	file_NfertMan != "")	{
				if(!NfertMan.Open(file_NfertMan, gridlist))
					fail("initio: could not open %s for input",(char*)file_NfertMan);
				readNman = true;
			}
		}
	}

	if(run_landcover) {
		file_Nfert_st=param["file_Nfert_st"].str;
		if(	file_Nfert_st != "")	{
			if(!Nfert_st.Open(file_Nfert_st, gridlist))
				fail("initio: could not open %s for input",(char*)file_Nfert_st);
			readNfert_st = true;
		}
		// See getwoodharvest() for input file format and content of text files
		file_woodharv_frac = param["file_woodharv_frac"].str;
		if(	file_woodharv_frac != "")	{
			if(!woodharv_frac.Open(file_woodharv_frac, gridlist))
				fail("initio: could not open %s for input",(char*)file_woodharv_frac);
			readwoodharvest_frac = true;
		}
		file_woodharv_cmass = param["file_woodharv_cmass"].str;
		if(	file_woodharv_cmass != "")	{
			if(!woodharv_cmass.Open(file_woodharv_cmass, gridlist))
				fail("initio: could not open %s for input",(char*)file_woodharv_cmass);
			readwoodharvest_cmass = true;
		}

		file_woodharv_prob_BL = param["file_woodharv_prob_BL"].str;
		if(file_woodharv_prob_BL != "")	{
			if(!woodharv_prob_BL.Open(file_woodharv_prob_BL, gridlist))
				fail("initio: could not open %s for input",(char*)file_woodharv_prob_BL);
			readwoodharv_prob_BL = true;
		}
		file_woodharv_prob_NL = param["file_woodharv_prob_NL"].str;
		if (file_woodharv_prob_NL != "") {
			if (!woodharv_prob_NL.Open(file_woodharv_prob_NL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_prob_NL);
			readwoodharv_prob_NL = true;
		}
		file_CC_prob_BL = param["file_CC_prob_BL"].str;
		if (file_CC_prob_BL != "") {
			if (!CC_prob_BL.Open(file_CC_prob_BL, gridlist))
				fail("initio: could not open %s for input", (char*)file_CC_prob_BL);
			readCC_prob_BL = true;
		}
		file_CC_prob_NL = param["file_CC_prob_NL"].str;
		if (file_CC_prob_NL != "") {
			if (!CC_prob_NL.Open(file_CC_prob_NL, gridlist))
				fail("initio: could not open %s for input", (char*)file_CC_prob_NL);
			readCC_prob_NL = true;
		}
		file_woodharv_intens_mean_BL = param["file_woodharv_intens_mean_BL"].str;
		if (file_woodharv_intens_mean_BL != "") {
			if (!woodharv_intens_mean_BL.Open(file_woodharv_intens_mean_BL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_intens_mean_BL);
			readwoodharv_intens_mean_BL = true;
		}
		file_woodharv_intens_mean_NL = param["file_woodharv_intens_mean_NL"].str;
		if (file_woodharv_intens_mean_NL != "") {
			if (!woodharv_intens_mean_NL.Open(file_woodharv_intens_mean_NL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_intens_mean_NL);
			readwoodharv_intens_mean_NL = true;
		}
		file_woodharv_intens_sd_BL = param["file_woodharv_intens_sd_BL"].str;
		if (file_woodharv_intens_sd_BL != "") {
			if (!woodharv_intens_sd_BL.Open(file_woodharv_intens_sd_BL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_intens_sd_BL);
			readwoodharv_intens_sd_BL = true;
		}
		file_woodharv_intens_sd_NL = param["file_woodharv_intens_sd_NL"].str;
		if (file_woodharv_intens_sd_NL != "") {
			if (!woodharv_intens_sd_NL.Open(file_woodharv_intens_sd_NL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_intens_sd_NL);
			readwoodharv_intens_sd_NL = true;
		}
		file_woodharv_size_preference_BL = param["file_woodharv_size_preference_BL"].str;
		if (file_woodharv_size_preference_BL != "") {
			if (!woodharv_size_preference_BL.Open(file_woodharv_size_preference_BL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_size_preference_BL);
			readwoodharv_size_preference_BL = true;
		}
		file_woodharv_size_preference_NL = param["file_woodharv_size_preference_NL"].str;
		if (file_woodharv_size_preference_NL != "") {
			if (!woodharv_size_preference_NL.Open(file_woodharv_size_preference_NL, gridlist))
				fail("initio: could not open %s for input", (char*)file_woodharv_size_preference_NL);
			readwoodharv_size_preference_NL = true;
		}
		file_woodharv_intensity_scalar = param["file_woodharv_intensity_scalar"].str;
		if(file_woodharv_intensity_scalar != "")	{
			if(!woodharv_intensity_scalar.Open(file_woodharv_intensity_scalar, gridlist))
				fail("initio: could not open %s for input",(char*)file_woodharv_intensity_scalar);
			readwoodharv_intensity_scalar = true;
		}

		// See getcutinterval() for input file format and content of text file
		file_cutinterval_st = param["file_cutinterval_st"].str;
		if(file_cutinterval_st != "")	{
			if(!cutinterval_st.Open(file_cutinterval_st, gridlist))
				fail("initio: could not open %s for input",(char*)file_cutinterval_st);
			readcutinterval_st = true;
		}
		file_firstmanageyear_st = param["file_firstmanageyear_st"].str;
		if(file_firstmanageyear_st != "")	{
			if(!firstmanageyear_st.Open(file_firstmanageyear_st, gridlist))
				fail("initio: could not open %s for input",(char*)file_firstmanageyear_st);
			readfirstmanageyear_st = true;
		}

		targetfrac_pft_mt = new TextInput::TimeDataD[nmt];

		// See gettargetcutting() for input file format and content of text file
		for(int i=0;i<nmt;i++) {
			ManagementType& mt = mtlist[i];
			if(mt.file_targetfrac_pft_mt != "") {
				if(!targetfrac_pft_mt[i].Open(mt.file_targetfrac_pft_mt, gridlist))
					fail("initio: could not open %s for input",(char*)mt.file_targetfrac_pft_mt);
				readtargetcutting = true;
			}
		}
	}

	file_burned_area = param["file_presc_burned_area"].str;
	if(file_burned_area != "")	{
		if(!presc_burned_area.Open(file_burned_area, gridlist))
			fail("initio: could not open %s for input",(char*)file_burned_area);
		readpresc_burned_area = true;
		if(presc_burned_area_start == FAR_PREHISTORIC_YEAR || presc_burned_area_end == FAR_FUTURE_YEAR)
			fail("presc_burned_area_start and presc_burned_area_end must be defined in the instruction file when burned area input is present\n");
	}

	if (run[PASTURE]) {
		file_grazing = param["file_grazing"].str;
		if (file_grazing != "") {
			if (!grazing.init(file_grazing)) {
				fail("initio: could not open %s for input", (const char *)file_grazing);
			}
			read_grazing = true;
		}
	}

	gridlist.killall();
}

bool ManagementInput::loadmanagement(double lon, double lat) {

	Coord c;
	c.lon = lon;
	c.lat = lat;
	bool LUerror = false;

	if(readpresc_burned_area) {
		if (!presc_burned_area.Load(c)) {
			log_info("Problems with prescribed burned area input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Prescribed burned area data not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if(readsowingdates) { 
		if(!sdates.Load(c)) {
			log_info("Problems with sowing date input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Sowing date data not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readharvestdates && !LUerror) {
		if(!hdates.Load(c)) {
			log_info("Problems with harvest date input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Harvest date data not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readNfert && !LUerror) {
		if(!Nfert.Load(c)) {
			log_info("Problems with N fertilization input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("N fertilization data not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readNman && !LUerror) {
		if(!NfertMan.Load(c)) {
			log_info("Manure data not found for %.2f,%.2f.",c.lon,c.lat);
		}
	}
	if(readNfert_st && !LUerror) {
		if(!Nfert_st.Load(c)) {
			log_info("Problems with N fertilization input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("N fertilization data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if(readwoodharvest_frac && !LUerror) {
		if(!woodharv_frac.Load(c)) {
			log_info("Problems with wood harvest fraction input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Wood harvest fraction data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if(readwoodharvest_cmass && !LUerror) {
		if(!woodharv_cmass.Load(c)) {
			log_info("Problems with wood harvest volume input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Wood harvest cmass data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if (readwoodharv_prob_BL && !LUerror) {
		if (!woodharv_prob_BL.Load(c)) {
			log_info("Problems with BL wood harvest probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("BL wood harvest probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_prob_NL && !LUerror) {
		if (!woodharv_prob_NL.Load(c)) {
			log_info("Problems with NL wood harvest probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("NL wood harvest probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readCC_prob_BL && !LUerror) {
		if (!CC_prob_BL.Load(c)) {
			log_info("Problems with BL wood harvest probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("BL wood harvest probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readCC_prob_NL && !LUerror) {
		if (!CC_prob_NL.Load(c)) {
			log_info("Problems with NL wood harvest probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Wood harvest fraction data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_intens_mean_BL && !LUerror) {
		if (!woodharv_intens_mean_BL.Load(c)) {
			log_info("Problems with BL tree thinning probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info(" BL tree thinning probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_intens_mean_NL && !LUerror) {
		if (!woodharv_intens_mean_NL.Load(c)) {
			log_info("Problems with NL tree thinning probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("NL tree thinning probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_intens_sd_BL && !LUerror) {
		if (!woodharv_intens_sd_BL.Load(c)) {
			log_info("Problems with BL tree thinning SD probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("BL tree thinning SD probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_intens_sd_NL && !LUerror) {
		if (!woodharv_intens_sd_NL.Load(c)) {
			log_info("Problems with NL tree thinning SD probability input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("NL tree thinning SD probability data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_size_preference_BL && !LUerror) {
		if (!woodharv_size_preference_BL.Load(c)) {
			log_info("Problems with thinning size preference input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Thinning size preference data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if (readwoodharv_size_preference_NL && !LUerror) {
		if (!woodharv_size_preference_NL.Load(c)) {
			log_info("Problems with thinning size preference input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("Thinning size preference data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readcutinterval_st && !LUerror) { 
		if(!cutinterval_st.Load(c)) {
			LUerror = true;	// skip this stand
			log_info("cutinterval data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}
	if(readfirstmanageyear_st && !LUerror) { 
		if(!firstmanageyear_st.Load(c)) {
			log_info("Problems with first management year input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("firstmanageyear data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readwoodharv_intensity_scalar && !LUerror) { 
		if(!woodharv_intensity_scalar.Load(c)) {
			log_info("Problems with probabilistic wood harvest intensity scalar input file. EXCLUDING STAND at %.3f,%.3f from simulation.", c.lon, c.lat);
			LUerror = true;	// skip this stand
			log_info("probabilistic wood harvest intensity scalar data not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if (read_grazing && !LUerror) {
		if (!grazing.load()) {
			LUerror = true; // skip this stand
			log_info("grazing data for stand types not found in input file for %.2f,%.2f.", c.lon, c.lat);
		}
	}

	if(readtargetcutting && !LUerror) {
		for(int i=0;i<nmt;i++) {
			ManagementType& mt = mtlist[i];
			if(mt.file_targetfrac_pft_mt != "") {
				if(!targetfrac_pft_mt[i].Load(c)) {
					LUerror = true;	// skip this stand
					log_info("Target cutting data for stand types not found in input file %s for %.2f,%.2f.",
						(char*)mt.file_targetfrac_pft_mt, c.lon, c.lat);
				}
				else if(mt.targetfrac_input_mode) {
					// This input is expected to be static or with non-zero values for pft:s in the selection the first input year.
					// selection must be set before stand creation
					char tempstring[200] = {0};
					if(targetfrac_pft_mt[i].isloaded()) {
						for(int pft=0; pft<npft; pft++)	{
							double pft_target = targetfrac_pft_mt[i].Get(0, pftlist[pft].name, true);
							if(pft_target != NOTFOUND) {
								if(pft_target) {
									strcat(tempstring, pftlist[pft].name);
									strcat(tempstring, " ");
								}
							}
						}
						mt.selection = tempstring;
					}
				}
			}
		}
	}
	return LUerror;
}

void ManagementInput::getsowingdates(Gridcell& gridcell) {

	if(!sdates.isloaded())
		return;

	int year = date.get_calendar_year();

	for(int i=0; i<npft; i++) {
		if(pftlist[i].phenology == CROPGREEN) {	

			gridcell.pft[i].sdate_force = (int)sdates.Get(year,pftlist[i].name);

			// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
			for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
				Standpft& standpft = gridcell[j].pft[i];
				if(standpft.active)
					standpft.sdate_force = gridcell.pft[i].sdate_force;
			}
		}
	}
}

void ManagementInput::getharvestdates(Gridcell& gridcell) {

	if(!hdates.isloaded())
		return;

	int year = date.get_calendar_year();

	for(int i=0; i<npft; i++)	{
		if(pftlist[i].phenology == CROPGREEN) {		

			gridcell.pft[pftlist[i].id].hdate_force = (int)hdates.Get(year,pftlist[i].name);

			// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
			for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
				Standpft& standpft = gridcell[j].pft[i];
				if(standpft.active)
					standpft.hdate_force = gridcell.pft[i].hdate_force;
			}
		}
	}
}

void ManagementInput::getNfert(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(Nfert.isloaded()) {
		for(int i=0; i<npft; i++)	{
			if(pftlist[i].phenology == CROPGREEN) {		
				gridcell.pft[pftlist[i].id].Nfert_read = Nfert.Get(year,pftlist[i].name);
				if(NfertMan.isloaded()) {
					gridcell.pft[pftlist[i].id].Nfert_man_read = NfertMan.Get(year,pftlist[i].name);
				}
			}
		}
	}

	if(Nfert_st.isloaded()) {
		for(int i=0; i<nst; i++)	{
			gridcell.st[stlist[i].id].nfert = Nfert_st.Get(year,stlist[i].name);
		}
	}
}

/// Load prescribed burned area from file
void ManagementInput::getburned_area(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	if(presc_burned_area.isloaded()) {
		//CLNdprintf("CLN getba \n");
		for(int i=0; i<365; i++)	{

			//CLNdprintf("CLN getba %d ",i);
			char day_srt[50] = {0};
			sprintf(day_srt, "%d", i);
			gridcell.prescribed_burned_area[i] = presc_burned_area.Get(year, day_srt);
			//CLNdprintf(" ba %f \n",gridcell.prescribed_burned_area[i]);
		}
	}
}
/// Read wood harvest from file
/** This implementation uses LUH2 wood harvest (area fraction or C mass) standard text format input files (see indata.h).
 *  Harvest of primary forest, primary non-forested land, mature secondary and young secondary forest and secondary non-forested
 *  land are represented in columns with header names "primf_harv", "primn_harv", "secmf_harv", "secyf_harv" and "secnf_harv" in
 *	gridcell area fraction input files and "primf_bioh", "primn_bioh", "secmf_bioh", "secyf_bioh" and "secnf_bioh" in C mass (kg/m2) input files.
 */
void ManagementInput::getwoodharvest(Gridcell& gridcell, LandcoverInput& landcover_input) {

	// Ignore harvest on "non-forested" land
	const bool ignore_non_forest_harvest = false;
	// Ignore all wood harvest on gridcells with smaller average cmass_wood than limit
	const double CMASS_WOOD_LIMIT = 0.1;
	int calender_year = date.get_calendar_year();
	int firsthistyear = landcover_input.getfirsthistyear();

	gridcell.landcover.wood_harvest.zero();

	/// Consistent with LUH2 transition input data: read previous year's values.
	if((calender_year >= firsthistyear + 1)) {

		double cmass_wood = 0.0;
		int nstands = 0;
		for(unsigned int i=0;i<gridcell.nbr_stands();i++) {
			Stand& stand = gridcell[i];
			if(stand.landcover == NATURAL) {
				for(unsigned int j=0;j<stand.nobj;j++) {
					Patch& patch = stand[j];
					cmass_wood += patch.cmass_wood() / (double)stand.npatch();
				}
				nstands++;
			}
		}
		cmass_wood /= double(nstands);

		if(cmass_wood < CMASS_WOOD_LIMIT)
			return;

		double frac_transfer;
		int year = calender_year - 1;

		if(woodharv_frac.isloaded()) {
			// Avoid using same data twice
			if(!ifprimary_to_secondary_transfer) {
				gridcell.landcover.wood_harvest.prim_frac += (frac_transfer = woodharv_frac.Get(year,"primf_harv")) != NOTFOUND ? frac_transfer : 0.0;
				if(!ignore_non_forest_harvest)
					gridcell.landcover.wood_harvest.prim_frac += (frac_transfer = woodharv_frac.Get(year,"primn_harv")) != NOTFOUND ? frac_transfer : 0.0;
			}
			gridcell.landcover.wood_harvest.sec_mature_frac += (frac_transfer = woodharv_frac.Get(year,"secmf_harv")) != NOTFOUND ? frac_transfer : 0.0;
			gridcell.landcover.wood_harvest.sec_young_frac += (frac_transfer = woodharv_frac.Get(year,"secyf_harv")) != NOTFOUND ? frac_transfer : 0.0;
			if(!ignore_non_forest_harvest)
				gridcell.landcover.wood_harvest.sec_young_frac += (frac_transfer = woodharv_frac.Get(year,"secnf_harv")) != NOTFOUND ? frac_transfer : 0.0;
		}
		if(woodharv_cmass.isloaded()) {
			// kg: divide by gridcell area (m2); NB. last parameter 3 if old cru climate (llc)
			// All removed C mass (including branches and attached leaves); kg/m2
			double cellsize =  pixelsize(gridcell.get_lat(), 0.5, 0.5, 0) * 1.0e6;	// km2 to m2
			gridcell.landcover.wood_harvest.prim_cmass += ((frac_transfer = woodharv_cmass.Get(year,"primf_bioh")) != NOTFOUND ? frac_transfer : 0.0) / cellsize;
			gridcell.landcover.wood_harvest.prim_cmass += ((frac_transfer = woodharv_cmass.Get(year,"primn_bioh")) != NOTFOUND ? frac_transfer : 0.0) / cellsize;
			gridcell.landcover.wood_harvest.sec_mature_cmass += ((frac_transfer = woodharv_cmass.Get(year,"secmf_bioh")) != NOTFOUND ? frac_transfer : 0.0) / cellsize;
			gridcell.landcover.wood_harvest.sec_young_cmass += ((frac_transfer = woodharv_cmass.Get(year,"secyf_bioh")) != NOTFOUND ? frac_transfer : 0.0) / cellsize;
			gridcell.landcover.wood_harvest.sec_young_cmass += ((frac_transfer = woodharv_cmass.Get(year,"secnf_bioh")) != NOTFOUND ? frac_transfer : 0.0) / cellsize;
		}
	}
}

/// Read clear-cutting interval from file
/** May use cutinterval output files from a simulation with automated clear-cutting along with e,g, a non-automated thinning
 *	scheme or no thinnings. The preferred input in that case is a time series with zero values in years without clear-cutting
 *  and the age of the clear-cut patch in years of clear-cutting (cutinterval_actual_thisyear). The default vales of keep_nonzero_value 
 *  and use_nextvalue are set to comply with this kind of input, but should cause no problems when using input with cutting interval
 *  values every year. If input file is present and no values found for a stand type, the default cutinterval_st value of 0
 *  will result in an absence of both thinnings and clearcuts,
 */
void ManagementInput::getcutinterval(Gridcell& gridcell) {

	int year = date.get_calendar_year();
	bool keep_nonzero_value = true;
	bool use_nextvalue = true;

	// Retrieve cutinterval for stand types
	for(int i=0; i<nst; i++) {
		if(cutinterval_st.isloaded()) {
			Gridcellst& gcst = gridcell.st[i];
			if(gcst.reset_cutinterval_st)
				gcst.cutinterval_st = 0.0;

			double cutinterval = cutinterval_st.Get(year, stlist[i].name, true);
			if(cutinterval != NOTFOUND) {
				if(!keep_nonzero_value || cutinterval) {
					gcst.cutinterval_st = cutinterval;
				}
				else if(use_nextvalue && !gcst.cutinterval_st) {
					// During the time input value is 0, use the first year with value > 0
					int year = date.get_calendar_year();
					int first_data_year = cutinterval_st.GetFirstyear();
					int last_data_year = first_data_year + cutinterval_st.GetnYears() - 1;
					// first_data_year is -1 for static data
					if(first_data_year == -1) {
						first_data_year = year;
						last_data_year = year;
					}
					// if value is 0 this year, try the following years
					for(int y=year;y<=max(year,last_data_year) + 1 && !cutinterval;y++) {
						cutinterval = cutinterval_st.Get(y, stlist[i].name, true);
						if(cutinterval)
							gcst.cutinterval_st = cutinterval;
					}
				}
			}
			gcst.reset_cutinterval_st = false;
		}
	}
}

/// Read firstmanageyear from file
/** Saves value in input file for stand type in gridcell.st[stid].firstmanageyear_st and overrides mt firstmanageyear value
 */
void ManagementInput::getfirstmanageyear(Gridcell& gridcell) {

	if(firstmanageyear_st.isloaded()) {
		for(int i=0; i<nst; i++)	{
			int firstmanageyear = (int)firstmanageyear_st.Get(0,stlist[i].name);
			if(firstmanageyear != NOTFOUND)
				gridcell.st[i].firstmanageyear_st = firstmanageyear;
		}
	}
}

/// Read firstmanageyear from file
/** Saves value in input file for stand type in gridcell.st[stid].firstmanageyear_st and overrides mt firstmanageyear value
 */
void ManagementInput::get_grazing(Gridcell& gridcell) {
	if(grazing.is_initialised()) {
		GrazingDate grazing_day(date.get_calendar_year(), date.month + 1, date.dayofmonth + 1);
		double animal_units = grazing.get_animal_units(grazing_day);
		for (int i=0; i<gridcell.nbr_stands(); i++){
			Stand& stand=gridcell[i];
			if (stand.landcover == PASTURE){
				stand.animal_units = animal_units;
			}
		}
	}
}

/// Read pft target cutting biomass fractions from file
/** File name with path, file_targetfrac_pft_mt, is defined in management types (or in stand types if if they don't have
 *  a rotation scheme). Columns for PFT values in text file are as in standard text input.In years with zero target cutting 
 *  sum, no target cutting will occur. To be certain to achieve target cutting at a specific time, at least targetcutinterval 
 *  number of consecutive years will have to have non-zero values.mt.targetfrac is not updated, but standpft.targetfrac is 
 *  filled with the values in the input file directly. During a year with management change (rotation), the instruction file 
 *  values of the target fractions are used rather than values in the input file to avoid inconsistencies, since the new 
 *	management type is loaded after the call to this function.
 */
void ManagementInput::gettargetcutting(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	for(int mtid=0;mtid<nmt;mtid++) {
		ManagementType& mt = mtlist[mtid];
		if(targetfrac_pft_mt[mtid].isloaded()) {
			for(int pft=0; pft<npft; pft++)	{
				double pft_target = targetfrac_pft_mt[mtid].Get(year,pftlist[pft].name, true);
				if(pft_target != NOTFOUND) {
					for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
						Stand& stand = gridcell[j];
						ManagementType& mt_cur = stand.get_current_management();
						if(mt_cur.name == mt.name) {
							Standpft& standpft = stand.pft[pft];
							standpft.targetfrac = pft_target;
						}
					}
				}
			}
		}
	}
}

void get_forest_structure(Patch& patch, bool& conifer, char* forest_structure_column) {

	double BA_BL = 0.0;
	double BA_NL = 0.0;
	double BA_tot = 0.0;

	double dens = 0.0;
	double diam_sq = 0.0;

	// Decide if broadleaved or needleleaved patch
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		if (indiv.pft.lifeform == TREE && indiv.diam >= 0.1) {		// trees > 10 cm diam; OK for BA too ?
			double BA = indiv.diam * indiv.diam * 0.25 * PI;
			if (indiv.pft.leafphysiognomy == BROADLEAF)
				BA_BL += BA * indiv.densindiv;
			else if (indiv.pft.leafphysiognomy == NEEDLELEAF)
				BA_NL += BA * indiv.densindiv;
			dens += indiv.densindiv;
			diam_sq += indiv.diam * indiv.diam * indiv.densindiv;
		}
		vegetation.nextobj();
	}
	if (BA_NL > BA_BL)
		conifer = true;
	else
		conifer = false;

	BA_tot = BA_NL + BA_BL;

	// Convert from m-2 to ha-1
	BA_tot *= M2_PER_HA;

	if (dens)
		diam_sq = pow(diam_sq / dens, 0.5);

	// Convert diam_sq m to cm
	diam_sq *= 100;

	int ba_bin = -1;
	int diam_bin = -1;

	for (int ba = 0; ba < nBA_bins; ba++) {
		double ba_low = BA_bin_min + ba * BA_bin_size;
		double ba_high = BA_bin_min + (ba + 1) * BA_bin_size;

		if (BA_tot >= ba_low && BA_tot < ba_high || iflast_bin_inclusive && ba == (nBA_bins - 1) && BA_tot >= ba_low)
			ba_bin = ba;
	}
	for (int d = 0; d < ndiam_bins; d++) {

		double d_low = diam_bin_min + d * diam_bin_size;
		double d_high = diam_bin_min + (d + 1) * diam_bin_size;

		if (diam_sq >= d_low && diam_sq < d_high || iflast_bin_inclusive && d == (ndiam_bins - 1) && d >= d_low)
			diam_bin = d;
	}

	char buffer[100] = { '\0'};
	sprintf(buffer, "D%d_BA%d", (diam_bin_min + diam_bin * diam_bin_size), (BA_bin_min + ba_bin * BA_bin_size));
	strcpy(forest_structure_column, buffer);
}

void ManagementInput::getwoodharvest_probability(Gridcell& gridcell) {

	int year = date.get_calendar_year();

	// Loop through all forest stands and patches that have harvest_system "PROBABILITY_INPUT"
	// Set harvest and clearcut probabilitied and thinning intensity means and standard deviations
	for (unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
		Stand& stand = gridcell[i];
		ManagementType& mt = stand.get_current_management();
		if (mt.harvest_system == "PROBABILITY_INPUT") {
			for (unsigned int j = 0; j < stand.nobj;j++) {
				Patch& patch = stand[j];
				if (!patch.management_started || !patch.cmass_wood())
					continue;
				char forest_structure_column[100] = { '\0'};
				bool conifer = false;
				get_forest_structure(patch, conifer, forest_structure_column);

				// If forest structure is not within any BA-diameter groups (as when mean diameter is < 10 cm), the variables are set to 0.
				if (conifer) {
					if (woodharv_prob_NL.isloaded()) {
						patch.harvest_probability = woodharv_prob_NL.Get(year, forest_structure_column, true);
						if (patch.harvest_probability == NOTFOUND)
							patch.harvest_probability = 0.0;
					}
					if (CC_prob_NL.isloaded()) {
						patch.clearcut_probability = CC_prob_NL.Get(year, forest_structure_column, true);
						if (patch.clearcut_probability == NOTFOUND)
							patch.clearcut_probability = 0.0;
					}
					if (woodharv_intens_mean_NL.isloaded()) {
						patch.thinning_intensity_mean = woodharv_intens_mean_NL.Get(year, forest_structure_column, true);
						if (patch.thinning_intensity_mean == NOTFOUND)
							patch.thinning_intensity_mean = 0.0;
					}
					if (woodharv_intens_sd_NL.isloaded()) {
						patch.thinning_intensity_sd = woodharv_intens_sd_NL.Get(year, forest_structure_column, true);
						if (patch.thinning_intensity_sd == NOTFOUND)
							patch.thinning_intensity_sd = 0.0;
					}
					if (woodharv_size_preference_NL.isloaded()) {
						patch.thinning_size_preference = (int)woodharv_size_preference_NL.Get(year, forest_structure_column, true);
						if (patch.thinning_size_preference == NOTFOUND)
							patch.thinning_size_preference = 0;
					}
					if (woodharv_intensity_scalar.isloaded()) {
						double harvest_intensity = woodharv_intensity_scalar.Get(year, "CF", true);
						if (harvest_intensity != NOTFOUND)
							patch.harvest_probability *= harvest_intensity;
					}
				}
				else {
					if (woodharv_prob_BL.isloaded()) {
						patch.harvest_probability = woodharv_prob_BL.Get(year, forest_structure_column, true);
						if (patch.harvest_probability == NOTFOUND)
							patch.harvest_probability = 0.0;
					}
					if (CC_prob_BL.isloaded()) {
						patch.clearcut_probability = CC_prob_BL.Get(year, forest_structure_column, true);
						if (patch.clearcut_probability == NOTFOUND)
							patch.clearcut_probability = 0.0;
					}
					if (woodharv_intens_mean_BL.isloaded()) {
						patch.thinning_intensity_mean = woodharv_intens_mean_BL.Get(year, forest_structure_column, true);
						if (patch.thinning_intensity_mean == NOTFOUND)
							patch.thinning_intensity_mean = 0.0;
					}
					if (woodharv_intens_sd_BL.isloaded()) {
						patch.thinning_intensity_sd = woodharv_intens_sd_BL.Get(year, forest_structure_column, true);
						if (patch.thinning_intensity_sd == NOTFOUND)
							patch.thinning_intensity_sd = 0.0;
					}
					if (woodharv_size_preference_BL.isloaded()) {
						patch.thinning_size_preference = (int)woodharv_size_preference_BL.Get(year, forest_structure_column, true);
						if (patch.thinning_size_preference == NOTFOUND)
							patch.thinning_size_preference = 0;
					}
					if (woodharv_intensity_scalar.isloaded()) {
						double harvest_intensity = woodharv_intensity_scalar.Get(year, "BF", true);
						if (harvest_intensity != NOTFOUND)
							patch.harvest_probability *= harvest_intensity;
					}
				}
			}
		}
	}
}

void ManagementInput::getmanagement(Gridcell& gridcell, LandcoverInput& landcover_input) {

	if (!date.day && readpresc_burned_area) {
		getburned_area(gridcell);
	}

	// check for grazing every day
	if (read_grazing) {
		get_grazing(gridcell);
	}

	if (!run_landcover || date.day) {
		return;
	}

	if(run[CROPLAND]) {

		//Read sowing dates from input file, put into gridcellpft.sdate_force
		if(readsowingdates)		
			getsowingdates(gridcell);
		//Read harvest dates from input file, put into gridcellpft.hdate_force
		if(readharvestdates)		
			getharvestdates(gridcell);
		//Read N fertilization from input file, put into gridcellpft.Nfert_read
		if (readNfert || readNfert_st)
			getNfert(gridcell);
	}
	// Read wood harvest from input file, put into gridcell.landcover.wood_harvest struct
	if(readwoodharvest_frac || readwoodharvest_cmass)
		getwoodharvest(gridcell, landcover_input);
	// Read cutting interval from input file
	if(readcutinterval_st)
		getcutinterval(gridcell);
	// Read selection target fractions from input file
	if(readtargetcutting)
		gettargetcutting(gridcell);
	// Read firstmanageyear from input file
	if(readfirstmanageyear_st)
		getfirstmanageyear(gridcell);

	// Read probabilistic wood harvest from input files; allow either thinning or clearcut components to be turned off
	if (readwoodharv_prob_BL && readwoodharv_prob_NL) {	
		if((readCC_prob_BL && readCC_prob_NL) || (readwoodharv_intens_mean_BL && readwoodharv_intens_mean_NL 
			&& readwoodharv_intens_sd_BL && readwoodharv_intens_sd_NL))
			getwoodharvest_probability(gridcell);
	}
}
