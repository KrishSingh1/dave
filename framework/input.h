////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file input.h
/// \brief Master class for all environmental input
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INPUT_H
#define INPUT_H

#define DYNAMIC_LANDCOVER_INPUT		// Reads landcover data from text files, using the TimeDataD class.

#include "guess.h"
#include "config.h"
#ifdef DYNAMIC_LANDCOVER_INPUT
#include "InData.h"
#endif
#include "inputdefinitions.h"
#include "landcoverinput.h"
#include "managementinput.h"
#include "lamarquendep.h"
#include "globalco2file.h"
#include "inputmodule.h"
#include <memory>

#define SOIL_INPUT_IN_CLIMATE_MODULE
#define NDEP_INPUT_IN_CLIMATE_MODULE

using std::auto_ptr;

// Forward declaration of Input
class Input;

/// Class for nitrogen deposition input
class NdepInput {

public:
	/// Constructor
	NdepInput(Input& in);
	/// Obtains nitrogen deposition data for current grid cell
	bool getgridcell(Gridcell& gridcell);
	/// Returns nitrogen deposition today
	double getndep(Gridcell& gridcell);

private:

	/// Reference to input container object
	Input& input;
	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;
	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;
	/// Static nitrogen deposition value optionally set in instruction file
	double ndep_fixed;
	/// Daily nitrogen deposition for current year
	double dndep[365];
};

/// Class for soil code and auxilliary code structure input
class SoilInput {

public:
	/// Constructor
	SoilInput(Input& in);
	/// Opens soilcode input file
	void init();
	/// Obtains soil data for current grid cell
	bool getgridcell(Gridcell& gridcell);
	/// Returns soil code for current grid cell
	int getsoilcode();

private:

	/// Reference to input container object
	Input& input;
	/// Reference to the list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord>& gridlist;
#if defined DYNAMIC_LANDCOVER_INPUT
	/// Soil code input object
	InData::TimeDataD soilcode;
#endif
	/// Soil code input file name
	xtring file_soilcode;

	/// Loads soilcode for current grid cell
	bool loadsoilcode(Gridcell& gridcell, Coord c);
};

/// Class containing all different environmental inputs
class Input {

public:

	/// Gridlist used by all different input classes
	ListArray_id<Coord> gridlist;
	/// The number of grid cells to simulate
	int ngridcell;

	/// Constructor
	Input(const char* climate_input_module_name, const char* landcover_input_module_name);
	/// Deconstructor
	~Input();

	/// Initiates all input classes
	void init();
	/// Obtains coordinates and loads environmental input for tcurrent grid cell
	bool getgridcell(Gridcell& gridcell);
	/// Obtains climate data, atmospheric (and nitrogen deposition) for one day
	bool getclimate(Gridcell& gridcell);
	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) {management_input_module->getmanagement(gridcell);}
	/// Obtains soil data for the next grid cell to simulate
	bool getsoil(Gridcell& gridcell, const int soilmap_index);
	/// Obtains atmospheric CO2 for one day
	double getco2(Gridcell& gridcell);
	/// Obtains nitrogen deposition for one day
	double getndep(Gridcell& gridcell) {return ndep_input.getndep(gridcell);}

	/// Returns the first historic year of simulation
	int getfirsthistyear() {return firsthistyear;}
	/// Returns the first historic year of simulation
	int getnyear_hist() {return nyear_hist;}

	/// Returns pointer to climate input module
	InputModule* get_climate_module() {return climate_input_module.get();}
	/// Returns pointer to land cover input module
	LandcoverInputModule* get_landcover_module() {return landcover_input_module.get();}
	/// Returns pointer to land management input module
	ManagementInputModule* get_management_module() {return management_input_module.get();}
	/// Reads gridlist from file
	void read_gridlist();

private:
	/// First historic year in simulation. May be set in instruction file.
	int firsthistyear;
	/// Last historic year in simulation. May be set in instruction file.
	int lasthistyear;
	/// Number of historic years in simulation. May be set in instruction file.
	int nyear_hist;

	/// Timer for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages

	/// Pointer to climate input module
	auto_ptr<InputModule> climate_input_module;
	/// Pointer to land cover input module
	auto_ptr<LandcoverInputModule> landcover_input_module;
	/// Pointer to land management input module
	auto_ptr<ManagementInputModule> management_input_module;

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;
	/// Static CO2 value optionally set in instruction file
	double co2_fixed;
	/// Nitrogen deposition input
	NdepInput ndep_input;
	/// Soil input
	SoilInput soil_input;

};

#endif // INPUT_H
