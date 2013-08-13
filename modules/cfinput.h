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

#include "cruinput.h"
#include "guessnc.h"
#include <memory>

class CFInput : public InputModule {
public:
	CFInput();

	~CFInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);

	static const int NYEAR_SPINUP_DATA=30;

private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int id;
		int rlon;
		int rlat;
		xtring descrip;
	};

	/// The grid cells to simulate
	std::vector<Coord> gridlist;

	/// The current grid cell to simulate
	std::vector<Coord>::iterator current_gridcell;

	/// Gets the first few years of data from cf_var and puts it into spinup_data
	void load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
	                      GenericSpinupData& spinup_data);

	/// Fills dtemp, dprec, etc. with forcing data for the current year
	void populate_daily_arrays();

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	GuessNC::CF::GridcellOrderedVariable* cf_temp;

	GuessNC::CF::GridcellOrderedVariable* cf_prec;

	GuessNC::CF::GridcellOrderedVariable* cf_insol;

	GenericSpinupData spinup_temp;

	GenericSpinupData spinup_prec;

	GenericSpinupData spinup_insol;

	/// Temperature for current gridcell and current year (deg C)
	double dtemp[365];

	/// Precipitation for current gridcell and current year (mm/day)
	double dprec[365];

	/// Insolation for current gridcell and current year (\see instype)
	double dinsol[365];

	// Daily N deposition for one year
	double dndep[365];

	/// Current timestep CF files
	int historic_timestep;

	/// Path to CRU binary archive
	xtring file_cru;

	/// Monthly data on daily dry NHx deposition (kgN/m2/day)
	double NHxDryDep[Lamarque::NYEAR_HISTNDEP][12];
	/// Monthly data on daily wet NHx deposition (kgN/m2/day)
	double NHxWetDep[Lamarque::NYEAR_HISTNDEP][12];
	/// Monthly data on daily dry NOy deposition (kgN/m2/day)
	double NOyDryDep[Lamarque::NYEAR_HISTNDEP][12];
	/// Monthly data on daily wet NOy deposition (kgN/m2/day)
	double NOyWetDep[Lamarque::NYEAR_HISTNDEP][12];

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether gridcell is divided into equal active landcover fractions.
	bool equal_landcover_area;
};

#endif // LPJ_GUESS_CFINPUT_H
