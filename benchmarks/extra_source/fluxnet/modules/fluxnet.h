///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.h
/// \brief Extra code used by the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_FLUXNET_H
#define LPJ_GUESS_FLUXNET_H

#include "cruinput.h"
#include "outputmodule.h"

/// Input module for Fluxnet benchmark
/** This is a subclass of the CRU input module. The subclass will alter the CRU
 *  forcing data according to site data.
 */
class FluxnetInput : public CRUInput {
public:
	bool getgridcell(Gridcell& gridcell);
	bool getclimate(Gridcell& gridcell);

private:
	std::vector<double> rain, tair, swrad;

	Lamarque::NDepData ndep;
	/// Daily N deposition for current year
	double dndep[Date::MAX_YEAR_LENGTH];
};

/// Output module for the extra files for FLUXNET
class FluxnetOutput : public GuessOutput::OutputModule {
public:
	FluxnetOutput();
	void init();
	void outannual(Gridcell& gridcell) {}
	void outdaily(Gridcell& gridcell);

private:
	xtring file_fluxnetdaily, file_fluxnetmonthly, file_fluxnetclim, file_monthlynee,
		   file_monthlyaet, file_monthlygpp, file_fluxnetmonth;
	GuessOutput::Table out_fluxnetmonthly, out_fluxnetdaily, out_fluxnetclim,
		ut_monthlynee, out_monthlyaet, out_monthlygpp, out_fluxnetmonth;

};

#endif // LPJ_GUESS_FLUXNET_H
