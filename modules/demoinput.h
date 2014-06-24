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

#include "guess.h"
#include "inputmodule.h"
#include <vector>
#include "gutil.h"

/// An input module for a toy data set (for demonstration purposes)
/** This input module is provided as an example of an input module.
 *  Included together with the LPJ-GUESS source code is a small
 *  toy data set (text based) which can be used to test the model,
 *  or to learn about how to write input modules.
 *
 *  \see InputModule for more documentation about writing input modules.
 */
class DemoInput : public InputModule {
public:

	/// Constructor
	/** Declares the instruction file parameters used by the input module.
	 */
	DemoInput();

	/// Destructor, cleans up used resources
	~DemoInput();

	/// Reads in gridlist and initialises the input module
	/** Gets called after the instruction file has been read */
	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	void getlandcover(Gridcell& gridcell);

private:

	/// Type for storing grid cell longitude, latitude and description text
	struct Coord {
		
		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	///	Loads landcover area fraction data from file(s) for a gridcell.
	/** Called from getgridcell() if run_landcover is true. 
	 */
	bool loadlandcover(Gridcell& gridcell, Coord c);

	/// Help function to readenv, reads in 12 monthly values from a text file
	void read_from_file(Coord coord, xtring fname, const char* format,
	                    double monthly[12], bool soil = false);

	/// Reads in environmental data for a location
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
	double dtemp[Date::MAX_YEAR_LENGTH];
	double dprec[Date::MAX_YEAR_LENGTH];
	double dsun[Date::MAX_YEAR_LENGTH];
	// bvoc
	// Daily diurnal temperature range for one year
	double ddtr[Date::MAX_YEAR_LENGTH];

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
