///////////////////////////////////////////////////////////////////////////////////////
/// \file cfinput.h
/// \brief Input module for CF conforming NetCDF files
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CFINPUT_H
#define LPJ_GUESS_CFINPUT_H

#ifdef HAVE_NETCDF

#include "cruinput.h"
#include "guessnc.h"
#include <memory>
#include <limits>
#include "input.h"

class CFInput : public InputModule {
public:
	CFInput();

	~CFInput();

	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);

	/// Returns hte spatial resolution of the gridlist
	double getgridlist_spatial_resolution() { return gridlist_spatial_resolution;}

	// Creates cf gridlist from lon-lat gridlist
	void create_cf_gridlist();

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) {management_input_module.getmanagement(gridcell);}

	/// Returns first historic year of climate input data
	int getfirsthistyear_climate();

	/// Returns number of years of climate input data
	int getnyear_hist_climate();

	/// Returns first historic year of simulation
	int getfirsthistyear();

	/// Returns number of historic years of simulation
	int getnyear_hist();

	bool supports_firsthistyear_in_insfile() { return false;}

	/// Returns pointer to land cover input module
	LandcoverInputModule* get_landcover_module() {return &landcover_input_module;}
	/// Returns pointer to land management input module
	ManagementInputModule* get_management_module() {return &management_input_module;}

	static const int NYEAR_SPINUP_DATA=30;

private:

	/// Land cover input module
	LandcoverInputModule landcover_input_module;
	/// Management input module
	ManagementInputModule management_input_module;
	
	/// List of Lon-Lat Coord objects containing coordinates of the grid cells to simulate (used to initate land cover and management input)
	ListArray_id<inputdef::Coord> gridlist;

	/// Spatial resolution of gridlist (degrees)
	double gridlist_spatial_resolution;

	/// search radius to use when finding soil data
	double searchradius;

	/// cf-specific Coord defenition
	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int rlon;
		int rlat;
		int landid;
		std::string descrip;
	};

	/// List of cf Coord objects containing NetCDF indeces and description of the grid cells
	std::vector<Coord> gridlistCF;

	/// The current grid cell to simulate
	std::vector<Coord>::iterator current_gridcell;

	/// Loads data from NetCDF files for current grid cell
	/** Returns the coordinates for the current grid cell, for
	 *  the closest CRU grid cell and the soilcode for the cell.
	 *  \returns whether it was possible to load data and find nearby CRU cell */
	bool load_data_from_files(double& lon, double& lat,
	                          double& cru_lon, double& cru_lat,
	                          int& soilcode);

	/// Gets the first few years of data from cf_var and puts it into spinup_data
	void load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
	                      GenericSpinupData& spinup_data);

	/// Gets data for one year, for one variable.
	/** Returns either 12 or 365/366 values (depending on LPJ-GUESS year length, not 
	 *  data set year length). Gets the values from spinup and/or historic period. */
	void get_yearly_data(std::vector<double>& data,
	                     const GenericSpinupData& spinup,
	                     GuessNC::CF::GridcellOrderedVariable* cf_historic,
	                     int& historic_timestep);

	/// Fills one array of daily values with forcing data for the current year
	void populate_daily_array(double* daily,
	                          const GenericSpinupData& spinup,
	                          GuessNC::CF::GridcellOrderedVariable* cf_historic,
	                          int& historic_timestep,
	                          double minimum = -std::numeric_limits<double>::max(),
	                          double maximum = std::numeric_limits<double>::max());

	/// Same as populate_daily_array, but for precipitation which is special
	/** Uses number of wet days if available and handles extensive/intensive conversion */
	void populate_daily_prec_array(long& seed);
	
	/// Fills dtemp, dprec, etc. with forcing data for the current year
	void populate_daily_arrays(long& seed);

	/// \returns all (used) variables
	std::vector<GuessNC::CF::GridcellOrderedVariable*> all_variables() const;

	/// Sets the spatial resolution of the climate data
	double parse_climate_spatial_resolution();

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	// The variables

	GuessNC::CF::GridcellOrderedVariable* cf_temp;

	GuessNC::CF::GridcellOrderedVariable* cf_prec;

	GuessNC::CF::GridcellOrderedVariable* cf_insol;

	GuessNC::CF::GridcellOrderedVariable* cf_wetdays;

	GuessNC::CF::GridcellOrderedVariable* cf_min_temp;

	GuessNC::CF::GridcellOrderedVariable* cf_max_temp;

	// Spinup data for each variable

	GenericSpinupData spinup_temp;

	GenericSpinupData spinup_prec;

	GenericSpinupData spinup_insol;

	GenericSpinupData spinup_wetdays;

	GenericSpinupData spinup_min_temp;
	
	GenericSpinupData spinup_max_temp;

	/// Temperature for current gridcell and current year (deg C)
	double dtemp[Date::MAX_YEAR_LENGTH];

	/// Precipitation for current gridcell and current year (mm/day)
	double dprec[Date::MAX_YEAR_LENGTH];

	/// Insolation for current gridcell and current year (\see instype)
	double dinsol[Date::MAX_YEAR_LENGTH];

	/// Daily N deposition for one year
	double dndep[Date::MAX_YEAR_LENGTH];

	/// Minimum temperature for current gridcell and current year (deg C)
	double dmin_temp[Date::MAX_YEAR_LENGTH];

	/// Maximum temperature for current gridcell and current year (deg C)
	double dmax_temp[Date::MAX_YEAR_LENGTH];

	/// Whether the forcing data for precipitation is an extensive quantity
	/** If given as an amount (kg m-2) per timestep it is extensive, if it's
	 *  given as a mean rate (kg m-2 s-1) it is an intensive quantity */
	bool extensive_precipitation;

	// Current timestep in CF files

	int historic_timestep_temp;

	int historic_timestep_prec;

	int historic_timestep_insol;

	int historic_timestep_wetdays;

	int historic_timestep_min_temp;

	int historic_timestep_max_temp;

	/// Path to CRU binary archive
	xtring file_cru;

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Nitrogen deposition time series to use (historic,rcp26,...)
	std::string ndep_timeseries;

	/// Spatial resolution of gridlist (degrees)
	double climate_spatial_resolution;

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_CFINPUT_H
