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

private:

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	std::auto_ptr<GuessNC::CF::GridcellOrderedVariable> cf_temp;

	std::auto_ptr<GuessNC::CF::GridcellOrderedVariable> cf_prec;

	std::auto_ptr<GuessNC::CF::GridcellOrderedVariable> cf_insol;

	int historic_timestep;
	int spinup_timestep;

	// Daily N deposition for one year
	double dndep[365];
};

#endif // LPJ_GUESS_CFINPUT_H
