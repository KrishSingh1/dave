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

class DemoInput : public InputModule {
public:

	DemoInput();

	~DemoInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);

	void getsowingdates(Gridcell& gridcell);

	void getharvestdates(Gridcell& gridcell);

private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	bool loadlandcover(Gridcell& gridcell, Coord c);

#if defined DYNAMIC_LANDCOVER_INPUT
	InData::Coord GetLonLat(Coord coord);

//void GetLonLatListFromCoord(ListArray_id<InData::Coord>&lonlatlist, ListArray_id<Coord>& gridlist);
	ListArray_id<InData::Coord> GetLonLatList(ListArray_id<Coord>& gridlist);
#endif

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

	//Whether enforced static landcover fractions are equal-sized stands of all included landcovers
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

	//Landuse:

#if defined DYNAMIC_LANDCOVER_INPUT

	InData::TimeDataD LUdata;
	InData::TimeDataD Peatdata;
	InData::TimeDataD CFTdata;
	InData::TimeDataD sdates;
	InData::TimeDataD hdates;

#ifdef LUTOMEMORY
	InData::TimeDataDmem LUdata_mem;
	InData::TimeDataDmem CFTdata_mem;
#endif

#endif
	xtring file_lu, file_lucrop, file_peat, file_sdates, file_hdates;
	static const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise
	static const int NYEAR_HIST=103; // Number of years in sowing date and harvest date input files
};

#endif // LPJ_GUESS_DEMOINPUT_H
