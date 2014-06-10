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
#include "guess.h"
#if defined DYNAMIC_LANDCOVER_INPUT
#include "InData.h"
#endif

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

	/// See base class for documentation about this function's responsibilities
	void getsowingdates(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	void getharvestdates(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	void getNfert(Gridcell& gridcell);

	int getfirsthistyear();	

private:

	/// Type for storing grid cell longitude, latitude and description text
	struct Coord {

		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	/// Loads landcover and crop fractions plus sowing and harvest dates from input files
	bool loadlandcover(Gridcell& gridcell, Coord c);

#if defined DYNAMIC_LANDCOVER_INPUT
	/// Transfers coordinates from CRUInput::Coord to InData::Coord
	InData::Coord GetLonLat(Coord coord);
	/// Transfers gridlist of coordinates from CRUInput::Coord to InData::Coord
	void GetLonLatList(ListArray_id<InData::Coord>& lonlatlist, ListArray_id<Coord>& gridlist);
#endif

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

	/// Whether enforced static landcover fractions are equal-sized stands of all included landcovers
	bool equal_landcover_area;

	/// Whether pfts not in crop fraction input file are removed from pftlist (0,1)
	bool minimizecftlist;

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

	//Landuse input:

#if defined DYNAMIC_LANDCOVER_INPUT
	// Objects handling landcover fraction data input
	InData::TimeDataD LUdata;
	InData::TimeDataD Peatdata;
	InData::TimeDataD CFTdata;
	InData::TimeDataD sdates;
	InData::TimeDataD hdates;
	InData::TimeDataD Nfert;
#ifdef LUTOMEMORY
	InData::TimeDataDmem LUdata_mem;
	InData::TimeDataDmem CFTdata_mem;
	InData::TimeDataDmem sdates_mem;
	InData::TimeDataDmem hdates_mem;
	InData::TimeDataDmem Nfert_mem;
#endif

#endif
	xtring file_lu, file_lucrop, file_peat, file_sdates, file_hdates, file_Nfert;
	// Number of years of landcover fraction data in input files
	static const int NYEAR_LU=103;
	// Number of years in sowing date and harvest date input files
	static const int NYEAR_HIST=103;
	// First historical year in landcover fraction and sowing/harvest date input files
	static const int FIRSTHISTYEAR=1901;
};

#endif // LPJ_GUESS_DEMOINPUT_H
