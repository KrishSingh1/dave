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
#include <gutil.h>

/// Type for storing Fluxnet grid cell information
struct FluxnetData {

	xtring desc;
	int start_y;
	int end_y;


	FluxnetData() {
		// initialise Fluxnet arrays with missing values;
		desc = "";
		start_y = 0;
		end_y = 0;


	}
};

/// Input module for Fluxnet benchmark
/** This is a subclass of the CRU input module. The subclass
 *  will alter the CRU forcing data according to site data,
 *  and also update the soiltype with soildepth information
 *  after the base class has initialized it according to soil code.
 */
class FluxnetInput : public CRUInput {
public:
	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);


private:

	std::vector<double> rain,Ta,swrad;
	std::vector<int> yr;
	Lamarque::NDepData ndep;
	std::map<std::pair<double, double>, FluxnetData> Fluxnetdata;
	/// Daily N deposition for current year
	double dndep[Date::MAX_YEAR_LENGTH];
};

/// Output module for the extra files for FLUXNET
class FluxnetOutput : public GuessOutput::OutputModule {
public:
	FluxnetOutput();

	void init();

	void outannual(Gridcell& gridcell);

	void outdaily(Gridcell& gridcell);


private:
	// Files for FLUXNET output and stats
	xtring file_fluxnetdaily, file_fluxnetmonthly, file_fluxnetclim, file_monthlynee, file_monthlyaet, file_monthlygpp,file_fluxnetmonth;
	// Output tables
	GuessOutput::Table out_fluxnetmonthly, out_fluxnetdaily, out_fluxnetclim,out_monthlynee,out_monthlyaet,out_monthlygpp,out_fluxnetmonth;

};

#endif // LPJ_GUESS_Fluxnet_H
