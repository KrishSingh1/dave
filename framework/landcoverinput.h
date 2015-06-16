////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file managementinput.h
/// \brief Input code for land cover management	from text files					
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LANDCOVERINPUT_H
#define LANDCOVERINPUT_H

#include "InData.h"

/// Reads gridlist in lon-lat-description format from text unput file
void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist);

/// Class that deals with all land cover input from text files
class LandcoverInputModule {

public:

	/// Constructor
	LandcoverInputModule();
	/// Opens management data files
	void init();
	/// Loads landcover and crop fractions from input files
	bool loadlandcover(Coord c);
	/// Gets land cover data for a year
	void getlandcover(Gridcell& gridcell);
	/// Gets land cover transition data for a year
	bool get_lc_transfer(Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES], double primary_lc_frac_transfer[][NLANDCOVERTYPES]);

private:

	// Objects handling land cover fraction data input
	InData::TimeDataD LUdata;
	InData::TimeDataD Peatdata;
	InData::TimeDataD grossLUC;
	InData::TimeDataD CFTdata;

	/// Files names for land cover fraction input files
	xtring file_lu, file_grossLUC, file_lucrop, file_peat;

	/// Whether pfts not in crop fraction input file are removed from pftlist (0,1)
	bool minimizecftlist;

	/// Number of years to increase cropland fraction linearly from 0 to first year's value
	int nyears_cropland_ramp;

	/// whether to use stand types with suitable rainfed crops (based on crop pft tb and gridcell latitude) when using fixed crop fractions
	bool frac_fixed_default_crops;
};

/// Class that deals with all crop management input from text files
class ManagementInputModule {

public:

	/// Constructor
	ManagementInputModule();
	/// Opens management data files
	void init();
	/// Loads fertilisation, sowing and harvest dates from input files
	bool loadmanagement(Coord c);
	/// Gets management data for a year
	void getmanagement(Gridcell& gridcell);

private:

	/// Input objects for each management text input file
	InData::TimeDataD sdates;
	InData::TimeDataD hdates;
	InData::TimeDataD Nfert;

	/// Files names for management input file
	xtring file_sdates, file_hdates, file_Nfert;

	/// Gets sowing date data for a year
	void getsowingdates(Gridcell& gridcell);
	/// Gets harvest date data for a year
	void getharvestdates(Gridcell& gridcell);
	/// Gets nitrogen fertilisation data for a year
	void getNfert(Gridcell& gridcell);
};

#endif // LANDCOVERINPUT_H