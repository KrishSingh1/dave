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
#include "guess.h"
#include "input.h"

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

	// Returns reference to gridlist
	ListArray_id<Coord>& getgridlist() { return gridlist;}

	/// Returns hte spatial resolution of the gridlist
	double getgridlist_spatial_resolution() { return gridlist_spatial_resolution;}

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) {management_input_module.getmanagement(gridcell);}

	bool getsoil(Gridcell& gridcell, const int soilmap_index);

	/// Returns first historic year of climate input data
	int getfirsthistyear_climate();

	/// Returns number of years of climate input data
	int getnyear_hist_climate();

	/// Returns first historic year of simulation
	int getfirsthistyear();

	/// Returns number of historic years of simulation
	int getnyear_hist();

	bool supports_firsthistyear_in_insfile() { return true;}

	/// Returns pointer to land cover input module
	LandcoverInputModule* get_landcover_module() {return &landcover_input_module;}
	/// Returns pointer to land management input module
	ManagementInputModule* get_management_module() {return &management_input_module;}


private:

	/// Land cover input module
	LandcoverInputModule landcover_input_module;
	/// Management input module
	ManagementInputModule management_input_module;

	/// Help function to readenv, reads in 12 monthly values from a text file
	bool read_from_file(Coord coord, xtring fname, const char* format,
	                    double monthly[12], bool soil = false);

	/// Reads in environmental data for a location
	bool readenv(Coord coord, long& seed);

	/// A list of Lon-Lat Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// Spatial resolution of gridlist (degrees)
	double gridlist_spatial_resolution;

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
};

#endif // LPJ_GUESS_DEMOINPUT_H
