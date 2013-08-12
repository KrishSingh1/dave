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

class CFInput : public CRUInput {
public:
	CFInput();

	~CFInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	static const int NYEAR_SPINUP_DATA=30;

private:

	/// Gets the first few years of data from cf_var and puts it into spinup_data
	void load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
	                      GenericSpinupData& spinup_data);

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

	int historic_timestep;
};

#endif // LPJ_GUESS_CFINPUT_H
