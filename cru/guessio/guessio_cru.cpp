///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio_cru.cpp
/// \brief LPJ-GUESS input/output module with input from instruction script
///
/// This I/O module reads in CRU climate data in a customised binary format.
/// The binary files contain CRU half-degree global historical climate data
/// for 1901-2006.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"

#ifdef USE_CRU_IO

#include "guessio.h"

#include "driver.h"
#include "parameters.h"
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>
#include "globalco2file.h"


// guess2008 - header file for the CRU TS 3.0 data archives
#include "cru_1901_2006.h"
#include "cru_1901_2006misc.h"

// header file for reading binary data archive of global nitrogen deposition
#include "GlobalNitrogenDeposition.h"


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH FILE SCOPE

double searchradius; // search radius to use when finding CRU data

/// Landcover fractions read from ins-file (% area).
int lc_fixed_frac[NLANDCOVERTYPES]={0};

/// Whether gridcell is divided into equal active landcover fractions.
bool equal_landcover_area;


namespace {
void initsettings() {

	searchradius = 0;
}
}

void construct_io() {
	initsettings();

	declare_parameter("searchradius", &searchradius, 0, 100,
		"If specified, CRU data will be searched for in a circle");

	declare_parameter("lcfrac_fixed", &lcfrac_fixed, "Whether static landcover fractions are set in the ins-file (0,1)");
	declare_parameter("equal_landcover_area", &equal_landcover_area, "Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
	declare_parameter("lc_fixed_urban", &lc_fixed_frac[URBAN], 0, 100, "% lc_fixed_urban");
	declare_parameter("lc_fixed_cropland", &lc_fixed_frac[CROPLAND], 0, 100, "% lc_fixed_cropland");
	declare_parameter("lc_fixed_pasture", &lc_fixed_frac[PASTURE], 0, 100, "% lc_fixed_pasture");
	declare_parameter("lc_fixed_forest", &lc_fixed_frac[FOREST], 0, 100, "% lc_fixed_forest");
	declare_parameter("lc_fixed_natural", &lc_fixed_frac[NATURAL], 0, 100, "% lc_fixed_natural");
	declare_parameter("lc_fixed_peatland", &lc_fixed_frac[PEATLAND], 0, 100, "% lc_fixed_peatland");
}

///////////////////////////////////////////////////////////////////////////////////////
//
//             SECTION: INPUT OF ENVIRONMENTAL DRIVING DATA FOR SIMULATION
//                            OUTPUT OF SIMULATION RESULTS
//
// In general it is the responsibility of the user of the model to provide code for
// this section of the input/output module. The following functions are called by the
// framework at various stages of the simulation and should contain appropriate code:
//
// void initio(const xtring& insfilename)
//   Initialises input/output (e.g. opening files), sets values for the global
//   simulation parameter variables (currently vegmode, npatch, patcharea,
//   ifbgestab, ifsme, ifstochestab, ifstochmort, iffire, estinterval,
//   npft), initialises pftlist (the one and only list of PFTs and their static
//   parameters for this run of the model). Normally all of the above parameters,
//   and possibly others, are read from the ins file (see above). Function readins
//   should be called to input settings from the ins file.
//
// bool getgridcell(Gridcell& gridcell)
//   Obtains coordinates and soil static parameters for the next grid cell to
//   simulate. The function should return false if no grid cells remain to be simulated,
//   otherwise true. Currently the following member variables of gridcell should be
//   initialised: longitude, latitude and climate.instype; the following members of
//   member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
//   thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
//   code (Sitch et al 2000) by a call to function soilparameters in the driver
//   module (driver.cpp):
//
//   soilparameters(gridcell.soiltype,soilcode);
//
//   If the model is to be driven by quasi-daily values of the climate variables
//   derived from monthly means, this function may be the appropriate place to
//   perform the required interpolations. The utility functions interp_monthly_means
//   and interp_monthly_totals in driver.cpp may be called for this purpose.
//
// bool getclimate(Gridcell& gridcell)
//   Obtains climate data (including atmospheric CO2 and insolation) for this day.
//   The function should return false if the simulation is complete for this grid cell,
//   otherwise true. This will normally require querying the year and day member
//   variables of the global class object date:
//
//   if (date.day==0 && date.year==nyear_spinup) return false;
//   // else
//   return true;
//
//   Currently the following member variables of the climate member of gridcell must be
//   initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
//   values of the climate variables derived from monthly means, this day's values
//   will presumably be extracted from arrays containing the interpolated daily
//   values (see function getgridcell):
//
//   gridcell.climate.temp=dtemp[date.day];
//   gridcell.climate.prec=dprec[date.day];
//   gridcell.climate.insol=dsun[date.day];
//
//   Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
//   BVOC:
//   gridcell.climate.dtr=ddtr[date.day]; 
//
//   If model is run in diurnal mode, which requires appropriate climate forcing data, 
//   additional members of the climate must be initialised: temps, insols. Both of the
//   variables must be of type std::vector. The length of these vectors should be equal
//   to value of date.subdaily which also needs to be set either in getclimate or 
//   getgridcell functions. date.subdaily is a number of sub-daily period in a single 
//   day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
//   diurnal mode.
//
// void outannual(Gridcell& gridcell)
//   Called at the end of the last day of each simulation year to permit output of
//   model results.
//
// termio()
//   Called after simulation is complete for all gridcells to allow memory deallocation,
//   closing of files or other cleanup functions.
//
// printhelp()
//   Prints out information about all available ins file parameters. Is typically
//   called when the user starts the program with the -help option.
//
///////////////////////////////////////////////////////////////////////////////////////

struct Coord {

	// Type for storing grid cell longitude, latitude and description text

	int id;
	double lon;
	double lat;
	xtring descrip;

};


ListArray_id<Coord> gridlist;
	// Will maintain a list of Coord objectsc ontaining coordinates
	// of the grid cells to simulate

int ngridcell; // the number of grid cells to simulate

class Spinup_data {

	// Class for management of climate data for spinup
	// (derived from first few years of historical climate data)

private:
	int nyear;
	int thisyear;
	double* data;
	bool havedata;

	// guess2008 - this array holds the climatology for the spinup period
	double dataclim[12];


public:
	Spinup_data(int nyear_loc) {
		nyear=nyear_loc;
		havedata=false;
		data=new double[nyear*12];
		if (!data) fail("Spinup_data::Spinup_data: out of memory");
		thisyear=0;
		havedata=true;
		reset_clim(); // guess2008
	}

	~Spinup_data() {
		if (havedata) delete[] data;
	}

	double& operator[](int month) {

		return data[thisyear*12+month];
	}

	void nextyear() {
		if (thisyear==nyear-1) thisyear=0;
		else thisyear++;
	}

	void firstyear() {
		thisyear=0;
	}

	void get_data_from(double source[][12]) {
		
		int y,m;
		thisyear=0; // guess2008 - ML bugfix
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m]=source[y][m];
			}
		}
	}

	// guess2008 - NEW METHODS 

	void reset_clim() {
		for (int ii = 0; ii < 12; ii++) dataclim[ii] = 0.0;
	}


	void make_clim() {
		
		reset_clim(); // Always reset before calculating

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				dataclim[m] += data[y*12+m] / (double)nyear;
			}
		}
	}


	bool extract_data(double source[][12], const int& startyear, const int& endyear) {
		
		// Populate data with data from the middle of source. 
		// Condition: endyear - startyear + 1 == nyear
		// if startyear == 1 and endyear == 30 then this function is identical to get_data_from above.

		if (endyear < startyear) return false;
		if (endyear - startyear + 1 == nyear) {

			int y,m;
			for (y=startyear-1;y<endyear;y++) {
				for (m=0;m<12;m++) {
					data[(y-(startyear-1))*12+m]=source[y][m];
				}
			}

		} else return false;

		return true;
	}


	void adjust_data(double anom[12], bool additive) {
		
		// Adjust the spinup data to the conditions prevailing at a particular time, as given by 
		// the (additive or multiplicative) anomalies in anom 
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (additive)	
					data[y*12+m] += anom[m];
				else
					data[y*12+m] *= anom[m];
			}
		}

	}
	
	
	// Replace interannual data with the period's climatology.
	void use_clim_data() {
	
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m];
			}
		}
	}


	// Alter variability about the mean climatology
	void adjust_data_variability(const double& factor) {
	
		// factor == 0 gives us the climatology (i.e. generalises use_clim_data above)
		// factor == 1 leaves everything unchanged
		// Remember to check the for negative precip or cloudiness values etc. 
		// after calling this method.

		if (factor == 1.0) return;

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m] + (data[y*12+m] - dataclim[m]) * factor;
			}
		}
	}


	void limit_data(double minval, double maxval) {

		// Limit data to a range
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < minval) data[y*12+m] = minval;
				if (data[y*12+m] > maxval) data[y*12+m] = maxval;
			}
		}

	}
	
	
	void set_min_val(const double& oldval, const double& newval) {

		// Change values < oldval to newval
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < oldval) data[y*12+m] = newval;
			}
		}

	}

	// guess2008 - END OF NEW METHODS


	void detrend_data() {

		int y,m;
		double a,b,anomaly;
		double* annual_mean=new double[nyear];
		double* year_number=new double[nyear];

		if (!annual_mean || !year_number)
			fail("Spinup_driver::detrend_data: out of memory");

		for (y=0;y<nyear;y++) {
			annual_mean[y]=0.0;
			for (m=0;m<12;m++) annual_mean[y]+=data[y*12+m];
			annual_mean[y]/=12.0;
			year_number[y]=y;
		}

		regress(year_number,annual_mean,nyear,a,b);

		for (y=0;y<nyear;y++) {
			anomaly=b*(double)y;
			for (m=0;m<12;m++)
				data[y*12+m]-=anomaly;
		}
		
		// guess2008 - added [] - Clean up
		delete[] annual_mean;
		delete[] year_number;
	}
};

// Constants associated with historical climate data set

// number of years of historical climate
// CRU TS 3.0 has 106 years of data (1901-2006)
const int NYEAR_HIST=106;

// calender year corresponding to first year in CRU climate data set
const int FIRSTHISTYEAR=1901;

// calender year corresponding to first year nitrogen deposition
const int FIRSTHISTYEARNDEP=1850;

// number of years of historical nitrogen deposition 
const int NYEAR_HISTNDEP=16;

// number of years to use for temperature-detrended spinup data set
// (not to be confused with the number of years to spinup model for, which
// is read from the ins file)	
const int NYEAR_SPINUP_DATA=30;

// Stream pointer to binary CRU historical climate data file (read from ins file)
FILE *in_cru;

// Full pathname of ASCII file containing annual CO2 values (read from ins file)
xtring file_co2;


// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

/// Yearly CO2 data read from file
/**
 * This object is indexed with calendar years, so to get co2 value for
 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
 * more information.
 */
GlobalCO2File co2;

// Monthly temperature, precipitation and sunshine data for current grid cell
// and historical period
double hist_mtemp[NYEAR_HIST][12];
double hist_mprec[NYEAR_HIST][12];
double hist_msun[NYEAR_HIST][12];

// Monthly frost days, precipitation days and DTR data for current grid cell
// and historical period
double hist_mfrs[NYEAR_HIST][12];
double hist_mwet[NYEAR_HIST][12];
double hist_mdtr[NYEAR_HIST][12];

/// Monthly data on daily dry NHx deposition (kgN/m2/day)
double NHxDryDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily wet NHx deposition (kgN/m2/day)
double NHxWetDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily dry NOy deposition (kgN/m2/day)
double NOyDryDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily wet NOy deposition (kgN/m2/day)
double NOyWetDep[NYEAR_HISTNDEP][12];

// Spinup data sets for current grid cell
Spinup_data spinup_mtemp(NYEAR_SPINUP_DATA);
Spinup_data spinup_mprec(NYEAR_SPINUP_DATA);
Spinup_data spinup_msun(NYEAR_SPINUP_DATA);

// guess2008
// Spinup data sets for monthly frost days, precipitation days and DTR data for 
// current grid cell
Spinup_data spinup_mfrs(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwet(NYEAR_SPINUP_DATA);
Spinup_data spinup_mdtr(NYEAR_SPINUP_DATA);


// Daily temperature, precipitation and sunshine for one year
double dtemp[365],dprec[365],dsun[365];
// bvoc
// Daily diurnal temperature range for one year
double ddtr[365];

// Daily N deposition for one year
double dndep[365];

// guess2008 - make file_cru and file_cru_misc global variables
xtring file_cru;
xtring file_cru_misc;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double mtemp[12], double mprec[12], double msun[12], double mdtr[12],
					double dtemp[365], double dprec[365], double dsun[365], double ddtr[365]) {
	interp_monthly_means(mtemp, dtemp);
	interp_monthly_totals(mprec, dprec);
	interp_monthly_means(msun, dsun);
	interp_monthly_means(mdtr, ddtr);
}

//Landuse:

//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
//TimeDataD input code may be put here
TimeDataD LUdata(LOCAL_YEARLY);
TimeDataD Peatdata;
#endif
xtring file_lu, file_peat;
const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 
// guess2008 - new functions for reading CRU TS 3.0 binary files.
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU
// Determine temp, precip, sunshine & soilcode
 
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
	double mtemp[NYEAR_HIST][12],double mprec[NYEAR_HIST][12],
	double msun[NYEAR_HIST][12]) {

	// !!!! NEW VERSION OF THIS FUNCTION - guess2008 - NEW VERSION OF THIS FUNCTION !!!!
	// Please note the new function signature. 

	// Archive object. Definition in new header file, cru.h
	Cru_1901_2006Archive ark;

	int target_ilon=(int)(dlon*10.0);
	int target_ilat=(int)(dlat*10.0);

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		soilcode=(int)data.soilcode[0];


		for (y=0;y<NYEAR_HIST;y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[y*12+m]*0.1; // now degC
				mprec[y][m] = data.mprec[y*12+m]*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				msun[y][m]  = data.msun[y*12+m]*0.1;   // % sun 

			}
		}


		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}




///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU_MISC
// Determine elevation, frs frq, wet frq & DTR

bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
	double mfrs[NYEAR_HIST][12],double mwet[NYEAR_HIST][12],
	double mdtr[NYEAR_HIST][12]) {
	
	// Please note the new function signature. 

	// Archive object
	Cru_1901_2006miscArchive ark; 
	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006misc data;

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays.
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000

		for (y=0;y<NYEAR_HIST;y++) { 
			for (m=0;m<12;m++) {

				// guess2008 - catch rounding errors 
				mfrs[y][m] = data.mfrs[y*12+m]*0.01; // days
				if (mfrs[y][m] < 0.1) 
					mfrs[y][m] = 0.0; // Catches rounding errors

				mwet[y][m] = data.mwet[y*12+m]*0.01; // days
				if (mwet[y][m] <= 0.1) 
					mwet[y][m] = 0.0; // Catches rounding errors

				mdtr[y][m] = data.mdtr[y*12+m]*0.1;  // degC

				/*
				If vapour pressure is needed:
				mvap[y][m] = data.mvap[y*12+m]*0.01;
				*/
			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


// Utility function that returns the CRU data from the nearest cell to (lon,lat) within
// a given search radius
bool findnearestCRUdata(double searchradius, char* cruark, double& lon, double& lat, 
                        int& scode, double hist_mtemp1[NYEAR_HIST][12], 
                        double hist_mprec1[NYEAR_HIST][12], 
                        double hist_msun1[NYEAR_HIST][12]) {

	// First try the exact coordinate
	if (searchcru(cruark, lon, lat, scode, hist_mtemp1, hist_mprec1, hist_msun1)) {
		return true;
	}
	
	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon*2)/2;
	double center_lat = floor(lat*2)/2;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.5;
	const double EPS = 1e-15;

	for (double y = center_lon-searchradius; y <= center_lon+searchradius+EPS; y += STEP) {
		for (double x = center_lat-searchradius; x <= center_lat+searchradius+EPS; x += STEP) {
			double xdist = x-center_lat;
			double ydist = y-center_lon;
			double dist = sqrt(xdist*xdist + ydist*ydist);
			
			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found in CRU
	for (unsigned int i = 0; i < search_points.size(); i++) {
		point search_point = search_points[i].second;
		double search_lon = search_point.first;
		double search_lat = search_point.second;

		if (searchcru(cruark, search_lon, search_lat, scode, 
		              hist_mtemp1, hist_mprec1, hist_msun1)) {
			lon = search_lon;
			lat = search_lat;
			return true;
		}
	}

	return false;
}


///////////////////////////////////////////////////////////////////////////////////////
// INITIO
// Called by the framework at the start of the model run

void initio(const xtring& insfilename) {

	// DESCRIPTION
	// Initialises input/output (e.g. opening files), sets values for the global
	// simulation parameter variables (currently vegmode, npatch, patcharea,
	// ifbgestab, ifsme, ifstochestab, ifstochmort, iffire,
	// estinterval, npft), initialises pftlist (the one and only list of PFTs and their
	// static parameters for this run of the model). Normally all of the above
	// parameters, and possibly others, are read from the ins file (see above).
	// Function readins should be called to input settings from the ins file.

	///////////////////////////////////////////////////////////////////////////////////
	// GENERIC SECTION - DO NOT MODIFY

	xtring header;

	unixtime(header);
	header=(xtring)"[LPJ-GUESS  "+header+"]\n\n";
	dprintf((char*)header);

	// Print the title of this run
	dprintf("\n\n-----------------------------------------------\n%s\n-----------------------------------------------\n",(char*)title);

	///////////////////////////////////////////////////////////////////////////////////
	// USER-SPECIFIC SECTION (Modify as necessary or supply own code)
	//
	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon,dlat;
	bool eof=false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;

	
	ngridcell=0;
	while (!eof) {
		
		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
			ngridcell++;
		}
	}


	fclose(in_grid);

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	if (run_landcover) {
		all_fracs_const=true;	//If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year

		//Retrieve file names for landcover files and open them if static values from ins-file are not used !
		if (!lcfrac_fixed) {	//This version does not support dynamic landcover fraction data

			if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
				file_lu=param["file_lu"].str;
#if defined DYNAMIC_LANDCOVER_INPUT
				if(!LUdata.Open(file_lu))				//Open Bondeau area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_lu);
				else if(LUdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}

			if (run[PEATLAND]) {	//special case for peatland: separate fraction file
				file_peat=param["file_peat"].str;
#if defined DYNAMIC_LANDCOVER_INPUT				
				if(!Peatdata.Open(file_peat))			//Open peatland area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_peat);
				else if(Peatdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}

		}
	}

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

///	Loads landcover area fraction data from file(s) for a gridcell.
/** Called from getgridcell() if run_landcover is true. 
  */
bool loadlandcover(Gridcell& gridcell, Coord c)	{
	bool LUerror=false;

	if (!lcfrac_fixed) {
		// Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
		// transferred to gridcell.landcoverfrac each year in getlandcover()

		if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
#if defined DYNAMIC_LANDCOVER_INPUT					
			if (!LUdata.Load(c))		//Load area fraction data from Bondeau input file to data object
			{
				dprintf("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;		// skip this stand
			}
#endif
		}

		if (run[PEATLAND] && !LUerror) {
#if defined DYNAMIC_LANDCOVER_INPUT
			if(!Peatdata.Load(c))	//special case for peatland: separate fraction file
			{
				dprintf("Problems with natural fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand						
			}
#endif
		}
	}

	return LUerror;
}

/// Retrieves nitrogen deposition for a particular gridcell
/** The values are either taken from a binary archive file or when it's not
 *  provided default to pre-industrial level of 2 kgN/ha/year.
 *
 *  The binary archive files have nitrogen deposition in gN/m2 on a monthly timestep
 *  for 16 years with 10 year interval starting from 1850 (Lamarque et. al., 2011).
 *
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 */
void getndep(double lon, double lat) {
	
	const double convert = 1e-7;				// converting from gN ha-1 to kgN m-2

	xtring file_ndep = param["file_ndep"].str;

	if (file_ndep == "") {

		// pre-industrial total nitrogen depostion set to 2 kgN/ha/year [kgN m-2]
		double dailyndep = 2000.0 / (4 * 365) * convert;

		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = dailyndep;
				NHxWetDep[y][m] = dailyndep;
				NOyDryDep[y][m] = dailyndep;
				NOyWetDep[y][m] = dailyndep;
			}
		}
	}
	else {
		GlobalNitrogenDepositionArchive ark;
		if (!ark.open(file_ndep)) {
			fail("Could not open %s for input", (char*)file_ndep);
		}

		GlobalNitrogenDeposition rec;
		rec.longitude = lon;
		rec.latitude = lat;

		if (!ark.getindex(rec)) {
			ark.close();
			fail("Grid cell not found in %s", (char*)file_ndep);
		}

		// Found the record, get the values
		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = rec.NHxDry[y*12+m] * convert;
				NHxWetDep[y][m] = rec.NHxWet[y*12+m] * convert;
				NOyDryDep[y][m] = rec.NOyDry[y*12+m] * convert;
				NOyWetDep[y][m] = rec.NOyWet[y*12+m] * convert;
			}
		}
		ark.close();
	}
}

/// Called by the framework at the start of the simulation for a particular grid cell
bool getgridcell(Gridcell& gridcell) {

	// DESCRIPTION
	// Obtains coordinates and soil static parameters for the next grid cell to
	// simulate. The function should return false if no grid cells remain to be simulated,
	// otherwise true. Currently the following member variables of Gridcell should be
	// initialised: longitude, latitude and climate.instype; the following members of
	// member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	// thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	// code (Sitch et al 2000) by a call to function soilparameters in the driver
	// module (driver.cpp):
	//
	// soilparameters(gridcell.soiltype,soilcode);
	//
	// If the model is to be driven by quasi-daily values of the climate variables
	// derived from monthly means, this function may be the appropriate place to
	// perform the required interpolations. The utility functions interp_monthly_means
	// and interp_monthly_totals in driver.cpp may be called for this purpose.

	// Select coordinates for next grid cell in linked list
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	bool gridfound;
	bool LUerror=false;

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	static bool first_call = true;

	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		
		// guess2008 - New searchcru functions takee the CRU filenames as their first 
		// argument, i.e. cru_1901_2002.bin and cru_1901_2002_misc.bin

		// New code:

		double lon = gridlist.getobj().lon;
		double lat = gridlist.getobj().lat;
		gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
		                               hist_mtemp, hist_mprec, hist_msun);

		if (gridfound) // Get more historical CRU data for this grid cell
			gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);

		if (run_landcover) {
			Coord& c=gridlist.getobj();
			LUerror=loadlandcover(gridcell, c);
		}
		if (LUerror)
			gridfound=false;

		while (!gridfound) {

			if (run_landcover && LUerror)
				dprintf("\nError: could not find stand at (%g,%g) in landcover data file\n", gridlist.getobj().lon,gridlist.getobj().lat);
			else
				dprintf("\nError: could not find stand at (%g,%g) in CRU data file\n", gridlist.getobj().lon,gridlist.getobj().lat);

			gridlist.nextobj();
			if (gridlist.isobj) {
				lon = gridlist.getobj().lon;
				lat = gridlist.getobj().lat;
				gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
				                               hist_mtemp, hist_mprec, hist_msun);
			  
				if (gridfound) // Get more historical CRU data for this grid cell
					gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation,
					                           hist_mfrs, hist_mwet, hist_mdtr);

				if (run_landcover) {
					Coord& c=gridlist.getobj();
					LUerror=loadlandcover(gridcell, c);
				}
				if (LUerror)
					gridfound=false;
			}
			else return false;
		}

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);
		spinup_mdtr.detrend_data();


		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);
		
		// Get nitrogen deposition data
		getndep(lon, lat);

		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		gridcell.climate.instype=SUNSHINE;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

///	Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
void getlandcover(Gridcell& gridcell) {
	int i, year;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;

	if(date.year<nyear_spinup)					//Use values for first historic year during spinup period !
		year=0;
	else if(date.year>=nyear_spinup+NYEAR_LU)	//AR4 adaptation
	{
//		dprintf("setting LU data for scenario period\n");
		year=NYEAR_LU-1;
	}
	else
		year=date.year-nyear_spinup;

	if(lcfrac_fixed)	// If area fractions are set in the ins-file.
	{
		if(date.year==0) // called by landcover_init
		{
			int nactive_landcovertypes=0;

			if(equal_landcover_area)
			{
				for(i=0;i<NLANDCOVERTYPES;i++)
				{
					if(run[i])
						nactive_landcovertypes++;
				}
			}

			for(i=0;i<NLANDCOVERTYPES;i++)
			{
				if(equal_landcover_area)
				{
					sum_active+=gridcell.landcoverfrac[i]=1.0*run[i]/(double)nactive_landcovertypes;	// only set fractions that are active !
					sum_tot=sum_active;
				}
				else
				{
					sum_tot+=gridcell.landcoverfrac[i]=(double)lc_fixed_frac[i]/100.0;					//count sum of all fractions (should be 1.0)

					if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)					//discard unreasonable values
					{
						if(date.year==0)
							dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
						sum_tot-=gridcell.landcoverfrac[i];
						gridcell.landcoverfrac[i]=0.0;
					}

					sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];				//only set fractions that are active !
				}
			}
			
			if(sum_tot<0.99 || sum_tot>1.01)	// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions
				if(date.year==0)
					dprintf("WARNING ! landcover fixed fraction sum is %4.2f, rescaling landcover fractions !\n", sum_tot);

				for(i=0;i<NLANDCOVERTYPES;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}

			//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
			if(sum_active<0.99)	//if landcover types are turned off in the ini-file, always <=1.0 here
			{
				if(date.year==0)
					dprintf("WARNING ! landcover active fraction sum is %4.2f.\n", sum_active);

				if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
				{
					if(date.year==0)
						dprintf("Inactive fractions (%4.2f) transferred to NATURAL fraction.\n", 1.0-sum_active);

					gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference 1.0-(sum of active landcover fractions) are added to the natural fraction
				}
				else
				{
/*					if(date.year==0)
						dprintf("Rescaling landcover fractions !\n");
					for(i=0;i<NLANDCOVERTYPES;i++)
						gridcell.landcoverfrac[i]/=sum_active;			// if NATURAL not simulated, rescale active fractions to 1.0
*/					if(date.year==0)
						dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
				}
																	
			}
		}
	}
	else	//area fractions are read from input file(s);
	{
		if(run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST])
		{	

			for(i=0;i<PEATLAND;i++)		//peatland fraction data is not in this file, otherwise i<NLANDCOVERTYPES.
			{	
#if defined DYNAMIC_LANDCOVER_INPUT
				sum_tot+=gridcell.landcoverfrac[i]=LUdata.Get(year,i);					//count sum of all fractions (should be 1.0)
#endif
				if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)			//discard unreasonable values
				{		
					if(date.year==0)
						dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
					sum_tot-=gridcell.landcoverfrac[i];
					gridcell.landcoverfrac[i]=0.0;
				}

				sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];
			}

			if(sum_tot!=1.0)		// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions

				if(sum_tot<0.99 || sum_tot>1.01)
				{
					if(date.year==0)
					{
						dprintf("WARNING ! landcover fraction sum is %4.2f for year %d\n", sum_tot, year+FIRSTHISTYEAR);
						dprintf("Rescaling landcover fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);
					}
				}
				else				//added scaling to sum=1.0 (sum often !=1.0)
					dprintf("Rescaling landcover fractions year %d ! (sum is within 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);

				for(i=0;i<PEATLAND;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}
		}
		else
			gridcell.landcoverfrac[NATURAL]=0.0;

		if(run[PEATLAND])
		{
#if defined DYNAMIC_LANDCOVER_INPUT
			sum_active+=gridcell.landcoverfrac[PEATLAND]=Peatdata.Get(year,"PEATLAND");			//peatland fraction data is currently in a separate file !
#endif
		}

		//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
		if(sum_active!=1.0)		//if landcover types are turned off in the ini-file, or if more landcover types are added in other input files, can be either less or more than 1.0
		{
			if(date.year==0)
				dprintf("Landcover fraction sum not 1.0 !\n");

			if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
			{
				if(date.year==0)
				{
					if(sum_active<1.0)
						dprintf("Inactive fractions (%4.3f) transferred to NATURAL fraction.\n", 1.0-sum_active);
					else
						dprintf("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).\n", sum_active-1.0, gridcell.landcoverfrac[NATURAL]);
				}

				gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					dprintf("New NATURAL fraction is %4.3f.\n", gridcell.landcoverfrac[NATURAL]);

				sum_active=1.0;		//sum_active should now be 1.0

				if(gridcell.landcoverfrac[NATURAL]<0.0)	//If new landcover type fraction is bigger than the natural fraction (something wrong in the distribution of input file area fractions)
				{										
					if(date.year==0)
						dprintf("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.\n");

					sum_active-=gridcell.landcoverfrac[NATURAL];	//fraction not possible to transfer moved back to sum_active, which will now be >1.0 again
					gridcell.landcoverfrac[NATURAL]=0.0;

					for(i=0;i<NLANDCOVERTYPES;i++)
					{
						gridcell.landcoverfrac[i]/=sum_active;		//fraction rescaled to unity sum
						if(run[i])
							if(date.year==0)
								dprintf("Landcover type %d fraction is %4.3f\n", i, gridcell.landcoverfrac[i]);
					}
				}
			}
			else
			{
//				if(date.year==0)
//					dprintf("Rescaling landcover fractions !\n");
//				for(i=0;i<NLANDCOVERTYPES;i++)
//					gridcell.landcoverfrac[i]/=sum_active;						// if NATURAL not simulated, rescale active fractions to 1.0
				if(date.year==0)
					dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
			}
		}
	}
}


/// Called by the framework each simulation day before any process modelling is performed for this day
/** Obtains climate data (including atmospheric CO2 and insolation) for this day. */
bool getclimate(Gridcell& gridcell) {

	// DESCRIPTION
	// The function should returns false if the simulation is complete for this grid cell,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of gridcell must be
	// initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	// values of the climate variables derived from monthly means, this day's values
	// will presumably be extracted from arrays containing the interpolated daily
	// values (see function getgridcell):
	//
	// gridcell.climate.temp=dtemp[date.day];
	// gridcell.climate.prec=dprec[date.day];
	// gridcell.climate.insol=dsun[date.day];
	// 
	// Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
	// BVOC:
	// gridcell.climate.dtr=ddtr[date.day]; 
	//
	// If model is run in diurnal mode, which requires appropriate climate forcing data, 
	// additional members of the climate must be initialised: temps, insols. Both of the
	// variables must be of type std::vector. The length of these vectors should be equal
	// to value of date.subdaily which also needs to be set either in getclimate or 
	// getgridcell functions. date.subdaily is a number of sub-daily period in a single 
	// day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
	// diurnal mode.

	double progress;

	// guess2008 - changed name from mwet to mwet_all
	double mwet_all[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // number of rain days per month
	Climate& climate = gridcell.climate;

	if (date.day == 0) {

		// First day of year ...

		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		int first_ndep_year = nyear_spinup + FIRSTHISTYEARNDEP - FIRSTHISTYEAR;

		double mndrydep[12], mnwetdep[12];
		int ndep_year = 0;

		if (date.year >= first_ndep_year) {
			ndep_year = (int)((date.year - first_ndep_year)/10);
		}

		for (int m = 0; m < 12; m++) {
			mndrydep[m] = NHxDryDep[ndep_year][m] + NOyDryDep[ndep_year][m];

			mnwetdep[m] = NHxWetDep[ndep_year][m] + NOyWetDep[ndep_year][m];
		}
		
		if (date.year < nyear_spinup) {

			// During spinup period

			int m;
			double mtemp[12],mprec[12],msun[12];
			double mfrs[12],mwet[12],mdtr[12];

			for (m=0;m<12;m++) {
				mtemp[m] = spinup_mtemp[m];
				mprec[m] = spinup_mprec[m];
				msun[m]	 = spinup_msun[m];

				// guess2008
				mfrs[m] = spinup_mfrs[m];
				mwet[m] = spinup_mwet[m];
				mdtr[m] = spinup_mdtr[m];
			}

			// Interpolate monthly spinup data to quasi-daily values
			interp_climate(mtemp,mprec,msun,mdtr,dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(mprec, dprec, mwet, gridcell.seed);
			}
			
			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);

			spinup_mtemp.nextyear();
			spinup_mprec.nextyear();
			spinup_msun.nextyear();

			// guess2008
			spinup_mfrs.nextyear();
			spinup_mwet.nextyear();
			spinup_mdtr.nextyear();

		}
		else if (date.year < nyear_spinup + NYEAR_HIST) {

			// Historical period

			// Interpolate this year's monthly data to quasi-daily values
			interp_climate(hist_mtemp[date.year-nyear_spinup],
				hist_mprec[date.year-nyear_spinup],hist_msun[date.year-nyear_spinup],
					   hist_mdtr[date.year-nyear_spinup],
				       dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year-nyear_spinup], dprec, hist_mwet[date.year-nyear_spinup], gridcell.seed);
			}

			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}
	}


	// Send environmental values for today to framework

	climate.co2 = co2[FIRSTHISTYEAR + date.year - nyear_spinup];

	climate.temp  = dtemp[date.day];
	climate.prec  = dprec[date.day];
	climate.insol = dsun[date.day];

	// Nitrogen deposition
	climate.dndep = dndep[date.day];

	// Nitrogen fertilization
	climate.dnfert = 0.0;
	
	// bvoc
	if(ifbvoc){
	  climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+NYEAR_HIST)
				+date.year)/(double)(ngridcell*(nyear_spinup+NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////
// TERMIO
// Called at end of model run (i.e. following simulation of all stands)

void termio() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
}

#endif // USE_CRU_IO


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Lamarque, J.-F., Kyle, G. P., Meinshausen, M., Riahi, K., Smith, S. J., Van Vuuren, 
//   D. P., Conley, A. J. & Vitt, F. 2011. Global and regional evolution of short-lived
//   radiatively-active gases and aerosols in the Representative Concentration Pathways. 
//   Climatic Change, 109, 191-212.
// Nakai, T., Sumida, A., Kodama, Y., Hara, T., Ohta, T. (2010). A comparison between
//   various definitions of forest stand height and aerodynamic canopy height.
//   Agricultural and Forest Meteorology, 150(9), 1225-1233
