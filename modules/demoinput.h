///////////////////////////////////////////////////////////////////////////////////////
/// \file demoinput.h
/// \brief Input module for demo data set
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_DEMOINPUT_H
#define LPJ_GUESS_DEMOINPUT_H

#include "inputmodule.h"
#include <vector>
#include "gutil.h"

class DemoInput : public InputModule {
public:

	DemoInput();

	~DemoInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);

private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	bool loadlandcover(Gridcell& gridcell, Coord c);

	void read_from_file(Coord coord, xtring fname, const char* format,
	                    double monthly[12], bool soil = false);

	bool readenv(Coord coord, long& seed);

	/// number of simulation years to run after spinup
	int nyear;

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether gridcell is divided into equal active landcover fractions.
	bool equal_landcover_area;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// The number of grid cells to simulate
	int ngridcell;

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages

	// Daily temperature, precipitation and sunshine for one year
	double dtemp[365],dprec[365],dsun[365];
	// bvoc
	// Daily diurnal temperature range for one year
	double ddtr[365];

	/// atmospheric CO2 concentration (ppmv) (read from ins file)
	double co2;

	/// atmospheric nitrogen deposition (kgN/yr/ha) (read from ins file)
	double ndep;

	//Landuse:

	//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
	//TimeDataD input code may be put here
	TimeDataD LUdata(LOCAL_YEARLY);
	TimeDataD Peatdata;
#endif
	xtring file_lu, file_peat;
	static const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise
};

#endif // LPJ_GUESS_DEMOINPUT_H
