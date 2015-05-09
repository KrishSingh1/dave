////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file managementinput.h
/// \brief Input code for land cover management	from text files					
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LANDCOVERINPUT_H
#define LANDCOVERINPUT_H

// Forward declaration
class InputModule;

/// Class that deals with all land cover input from text files
class LandcoverInputModule {

public:

	/// Constructor
	LandcoverInputModule(InputModule& in);
	/// Opens management data files
	void init();
	/// Loads landcover and crop fractions from input files
	bool loadlandcover(Gridcell& gridcell, Coord c);
	/// Gets land cover data for a year
	void getlandcover(Gridcell& gridcell);
	/// Gets land cover transition data for a year
	bool get_lc_transfer(Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES], double primary_lc_frac_transfer[][NLANDCOVERTYPES]);
	/// Gets first year of land cover data
	int getfirsthistyear();
	/// Gets number of years of land cover data
	int getnyear_hist();

private:

	/// Reference to the InputModule object
	InputModule& input;

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
};

#endif // LANDCOVERINPUT_H