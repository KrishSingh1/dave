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

	/// Gets landcover and crop fractions for one year
	void getlandcover(Gridcell& gridcell);

	/// Gets sowing dates for one year
	void getsowingdates(Gridcell& gridcell);

	/// Gets harvest dates for one year
	void getharvestdates(Gridcell& gridcell);

	/// Gets N fertilization for one year
	void getNfert(Gridcell& gridcell);

private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
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

	void read_from_file(Coord coord, xtring fname, const char* format,
	                    double monthly[12], bool soil = false);

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
