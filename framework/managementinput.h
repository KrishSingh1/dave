////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file managementinput.h
/// \brief Input code for land cover management	from text files					
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MANAGEMENTINPUT_H
#define MANAGEMENTINPUT_H

class ManagementInputModule {

public:

	ManagementInputModule(Input& in);

	void init();

	bool getgridcell(Gridcell& gridcell);
	void getmanagement(Gridcell& gridcell);
	void getsowingdates(Gridcell& gridcell);
	void getharvestdates(Gridcell& gridcell);
	void getNfert(Gridcell& gridcell);

private:

	/// Reference to input container object
	Input& input;

	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;

	InData::TimeDataD sdates;
	InData::TimeDataD hdates;
	InData::TimeDataD Nfert;

	xtring file_sdates, file_hdates, file_Nfert;

	/// Loads fertilisation, sowing and harvest dates from input files
	bool loadmanagement(Gridcell& gridcell, Coord c);
};

#endif // MANAGEMENTINPUT_H