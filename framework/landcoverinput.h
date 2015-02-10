//#define USE
#ifdef USE
#ifndef LANDCOVERINPUTMODULE_H
#define LANDCOVERINPUTMODULE_H

#include "input.h"

class LandcoverInputModule {

public:

	LandcoverInputModule(const char* input_module_name, Input& in);

	void init();

	bool getgridcell(Gridcell& gridcell);
	void getlandcover(Gridcell& gridcell);
	bool get_lc_transfer(Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES], double primary_lc_frac_transfer[][NLANDCOVERTYPES]);
	int getfirsthistyear();
	int getnyear_hist();

private:

	/// Reference to input container object
	Input& input;

	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;

#if defined DYNAMIC_LANDCOVER_INPUT
	// Objects handling landcover fraction data input
	InData::TimeDataD LUdata;
	InData::TimeDataD Peatdata;
	InData::TimeDataD grossLUC;
	InData::TimeDataD CFTdata;
#endif
	xtring file_lu, file_grossLUC, file_lucrop, file_peat;

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether enforced static landcover fractions are equal-sized stands of all included landcovers
	bool equal_landcover_area;

	/// Whether pfts not in crop fraction input file are removed from pftlist (0,1)
	bool minimizecftlist;

	/// Number of years to increase cropland fraction linearly from 0 to first year's value
	int nyears_cropland_ramp;

	/// Loads landcover and crop fractions plus sowing and harvest dates from input files
	bool loadlandcover(Gridcell& gridcell, Coord c);
};

#endif
#endif