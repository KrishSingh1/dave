///////////////////////////////////////////////////////////////////////////////////////
/// \file demoinput.cpp
/// \brief LPJ-GUESS input module for a toy data set (for demonstration purposes)
///
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "demoinput.h"

#include "driver.h"
#include "outputchannel.h"
#include <stdio.h>

REGISTER_INPUT_MODULE("demo", DemoInput)

// Anonymous namespace for variables and functions with file scope
namespace {

// File names for temperature, precipitation, sunshine and soil code driver files
xtring file_temp,file_prec,file_sun,file_soil;

// LPJ soil code
int soilcode;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double* mtemp, double* mprec, double* msun, double* mdtr,
					double* dtemp, double* dprec, double* dsun, double* ddtr) {
	interp_monthly_means_conserve(mtemp, dtemp);
	interp_monthly_totals_conserve(mprec, dprec, 0);
	interp_monthly_means_conserve(msun, dsun, 0, 100);
	interp_monthly_means_conserve(mdtr, ddtr, 0);
}

} // namespace


DemoInput::DemoInput(Input& in) 
	: gridlist(in.gridlist),
	  input(in) {
#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
	// Retrieve specified N value as read from ins file
	declare_parameter("ndep_fixed", &ndep_fixed, 0.0, 2000.0,"Fixed ndep value");
#endif
}

void DemoInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files), and reads in the gridlist

	// Retrieve input file names as read from ins file

	file_temp=param["file_temp"].str;
	file_prec=param["file_prec"].str;
	file_sun=param["file_sun"].str;
#ifdef SOIL_INPUT_IN_CLIMATE_MODULE
	file_soil=param["file_soil"].str;
#endif
	if(search_for_centre_of_gridcell)
		fail("Demo input does not support searchradius\n");
}


bool DemoInput::read_from_file(Coord coord, xtring fname, const char* format,
                               double monthly[12], bool soil /* = false */) {
	double dlon, dlat;
	int elev;
	FILE* in = fopen(fname, "r");
	if (!in) {
		fail("readenv: could not open %s for input", (char*)fname);
	}

	bool foundgrid = false;
	while (!feof(in) && !foundgrid) {
		if (!soil) {
			readfor(in, format, &dlon, &dlat, &elev, monthly);
		} else {
			readfor(in, format, &dlon, &dlat, &soilcode);
		}
		foundgrid = equal(coord.lon, dlon) && equal(coord.lat, dlat);
	}


	fclose(in);
	if (!foundgrid) {
		dprintf("readenv: could not find record for (%g,%g) in %s",
										coord.lon, coord.lat, (char*)fname);
	}
	return foundgrid;
}

bool DemoInput::readenv(Coord coord, long& seed) {

	// Searches for environmental data in driver temperature, precipitation,
	// sunshine and soil code files for the grid cell whose coordinates are given by
	// 'coord'. Data are written to arrays mtemp, mprec, msun and the variable
	// soilcode, which are defined as global variables in this file

	// The temperature, precipitation and sunshine files (Cramer & Leemans,
	// unpublished) should be in ASCII text format and contain one-line records for the
	// climate (mean monthly temperature, mean monthly percentage sunshine or total
	// monthly precipitation) of a particular 0.5 x 0.5 degree grid cell. Elevation is
	// also included in each grid cell record.

	// The following sample record from the temperature file:
	//   " -4400 8300 293-298-311-316-239-105  -7  26  -3 -91-184-239-277"
	// corresponds to the following data:
	//   longitude 44 deg (-=W)
	//   latitude 83 deg (+=N)
	//   elevation 293 m
	//   mean monthly temperatures (deg C) -29.8 (Jan), -31.1 (Feb), ..., -27.7 (Dec)

	// The following sample record from the precipitation file:
	//   " 12750 -200 223 190 165 168 239 415 465 486 339 218 162 149 180"
	// corresponds to the following data:
	//   longitude 127.5 deg (+=E)
	//   latitude 20 deg (-=S)
	//   elevation 223 m
	//   monthly precipitation sum (mm) 190 (Jan), 165 (Feb), ..., 180 (Dec)

	// The following sample record from the sunshine file:
	//   "  2600 7000 293  0 20 38 37 31 28 28 25 21 17  9  7"
	// corresponds to the following data:
	//   longitude 26 deg (+=E)
	//   latitude 70 deg (+=N)
	//   elevation 293 m
	//   monthly mean %age of full sunshine 0 (Jan), 20 (Feb), ..., 7 (Dec)

	// The LPJ soil code file is in ASCII text format and contains one-line records for
	// each grid cell in the form:
	//   <lon> <lat> <soilcode>
	// where <lon>      = longitude as a floating point number (-=W, +=E)
	//       <lat>      = latitude as a floating point number (-=S, +=N)
	//       <soilcode> = integer in the range 0 (no soil) to 9 (see function
	//                    soilparameters in driver module)
	// The fields in each record are separated by spaces

	double mtemp[12];		// monthly mean temperature (deg C)
	double mprec[12];		// monthly precipitation sum (mm)
	double msun[12];		// monthly mean percentage sunshine values

	double mwet[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // number of rain days per month

	double mdtr[12];		// monthly mean diurnal temperature range (oC)
	for(int m=0; m<12; m++) {
		mdtr[m] = 0.;
		if (ifbvoc) {
			dprintf("WARNING: No data available for dtr in sample data set!\nNo daytime temperature correction for BVOC calculations applied.");
		}
	}

	bool gridfound = read_from_file(coord, file_temp, "f6.2,f5.2,i4,12f4.1", mtemp);
	if(gridfound)
		gridfound = read_from_file(coord, file_prec, "f6.2,f5.2,i4,12f4", mprec);
	if(gridfound)
		gridfound = read_from_file(coord, file_sun, "f6.2,f5.2,i4,12f3", msun);
#ifdef SOIL_INPUT_IN_CLIMATE_MODULE
	read_from_file(coord, file_soil, "f,f,i", msun, true);	// msun is not used here: just dummy
#endif
	if(gridfound) {
		// Interpolate monthly values for environmental drivers to daily values
		// (relevant daily values will be sent to the framework each simulation
		// day in function getclimate, below)
		interp_climate(mtemp, mprec, msun, mdtr, dtemp, dprec, dsun, ddtr);

		// Recalculate precipitation values using weather generator
		// (from Dieter Gerten 021121)
		prdaily(mprec, dprec, mwet, seed);
	}

	return gridfound;
}


bool DemoInput::getgridcell(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	// Select coordinates for next grid cell in linked list
	Coord& c = gridlist.getobj();

	// Load climate data for this grid cell from files
	bool gridfound = readenv(c, gridcell.seed);

	if(gridfound) {
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine	
		gridcell.climate.instype = SUNSHINE;
#ifdef SOIL_INPUT_IN_CLIMATE_MODULE
		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype, soilcode);
#endif
	}

	return gridfound;
}


int DemoInput::getfirsthistyear() {

	return -1;
}

int DemoInput::getnyear_hist() {

	return 0;
}

bool DemoInput::getclimate(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	Climate& climate = gridcell.climate;


	// Send environmental values for today to framework
#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
	climate.dndep  = ndep_fixed / (365.0 * 10000.0);
#endif
	climate.temp  = dtemp[date.day];
	climate.prec  = dprec[date.day];
	climate.insol = dsun[date.day];

	// bvoc
	climate.dtr=ddtr[date.day];

/*	if (date.day == 0) {

		// Return false if last year was the last for the simulation
		if (date.year == nyear_spinup + nyear) return false;
	}
*/
	return true;
}

bool DemoInput::getsoil(Gridcell& gridcell, const int soilmap_index){
	return true;
}
