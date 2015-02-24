
#ifndef INPUT_H
#define INPUT_H

#include "guess.h"
#include "config.h"
#ifdef DYNAMIC_LANDCOVER_INPUT
#include "InData.h"
#endif
#include "inputdefinitions.h"
#include "lamarquendep.h"
#include "globalco2file.h"
#include "inputmodule.h"
#include <memory>

#define SOIL_INPUT_IN_CLIMATE_MODULE
#define NDEP_INPUT_IN_CLIMATE_MODULE

using std::auto_ptr;
//using namespace inputdef;

class Input;

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

#if defined DYNAMIC_LANDCOVER_INPUT
	InData::TimeDataD sdates;
	InData::TimeDataD hdates;
	InData::TimeDataD Nfert;
#endif

	xtring file_sdates, file_hdates, file_Nfert;

	/// Loads fertilisation, sowing and harvest dates from input files
	bool loadmanagement(Gridcell& gridcell, Coord c);
};

class NdepInput {

public:

	NdepInput(Input& in);
	bool getgridcell(Gridcell& gridcell);
	double getndep(Gridcell& gridcell);

private:

	/// Reference to input container object
	Input& input;

	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;
	double ndep_fixed;

	/// Daily N deposition for current year
	double dndep[365];
};

class SoilInput {

public:

	SoilInput(Input& in);
	void init();
	bool getgridcell(Gridcell& gridcell);
	int getsoilcode();

private:

	/// Reference to input container object
	Input& input;

	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;

#if defined DYNAMIC_LANDCOVER_INPUT
	InData::TimeDataD soilcode;
#endif

	xtring file_soilcode;

	bool loadsoilcode(Gridcell& gridcell, Coord c);
};


class Input {

public:

	ListArray_id<Coord> gridlist;
	/// The number of grid cells to simulate
	int ngridcell;
	int firsthistyear;
	int lasthistyear;
	int nyear_hist;

	Input(const char* climate_input_module_name, const char* landcover_input_module_name);
	~Input();

	void init();
	bool getgridcell(Gridcell& gridcell);
	bool getclimate(Gridcell& gridcell);
	void getmanagement(Gridcell& gridcell);
	bool getsoil(Gridcell& gridcell, const int soilmap_index);
	double getco2(Gridcell& gridcell);
	double getndep(Gridcell& gridcell);
	int getfirsthistyear();
	InputModule* get_climate_module();
	LandcoverInputModule* get_landcover_module();
	ManagementInputModule* get_management_module();
	void read_gridlist();

private:

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages

	auto_ptr<InputModule> climate_input_module;
	auto_ptr<LandcoverInputModule> landcover_input_module;
	auto_ptr<ManagementInputModule> management_input_module;

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;
	double co2_fixed;

	NdepInput ndep_input;
	SoilInput soil_input;

};

#endif // INPUT_H
