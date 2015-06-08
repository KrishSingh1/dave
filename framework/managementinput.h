////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file managementinput.h
/// \brief Input code for land cover management	from text files					
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MANAGEMENTINPUT_H
#define MANAGEMENTINPUT_H

/// Class that deals with all crop management input from text files
class ManagementInputModule {

public:

	/// Constructor
	ManagementInputModule(InputModule& in);
	/// Opens management data files
	void init();
	/// Loads fertilisation, sowing and harvest dates from input files
	bool loadmanagement(Coord c);
	/// Gets management data for a year
	void getmanagement(Gridcell& gridcell);

private:

	/// Reference to the InputModule object
	InputModule& input;

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

#endif // MANAGEMENTINPUT_H