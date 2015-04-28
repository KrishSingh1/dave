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
	DemoInput(Input&);

	/// Reads in gridlist and initialises the input module
	/** Gets called after the instruction file has been read */
	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);

	bool getsoil(Gridcell& gridcell, const int soilmap_index);

	int getfirsthistyear();

	int getnyear_hist();

	double* getdprec() {return dprec;};

	bool supports_firsthistyear_in_insfile() { return true;}

///
private:

	/// Reference to input container object
	Input& input;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;

	/// Help function to readenv, reads in 12 monthly values from a text file
	bool read_from_file(Coord coord, xtring fname, const char* format,
	                    double monthly[12], bool soil = false);

	/// Reads in environmental data for a location
	bool readenv(Coord coord, long& seed);

	// Daily temperature, precipitation and sunshine for one year
	double dtemp[Date::MAX_YEAR_LENGTH];
	double dprec[Date::MAX_YEAR_LENGTH];
	double dsun[Date::MAX_YEAR_LENGTH];
	// bvoc
	// Daily diurnal temperature range for one year
	double ddtr[Date::MAX_YEAR_LENGTH];

	/// atmospheric nitrogen deposition (kgN/yr/ha) (read from ins file)
	double ndep_fixed;
};

#endif // LPJ_GUESS_DEMOINPUT_H
