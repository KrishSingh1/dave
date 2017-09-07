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

/// The value used for missing data in the Fluxnet files.
const double MISSING_DATA = -9999.0;

/// The number of Fluxnet years, currently 1996-2006, inclusive.
const int NFLUXYEARS=11;

/// Type for storing Fluxnet grid cell information
struct FluxnetData {

	xtring desc;
	xtring desc2;
	int start_y;
	int end_y;



	int isfluxdata[NFLUXYEARS];




	// Modelled flux data for the same site
	double modelNEE[NFLUXYEARS][12];
	double modelAET[NFLUXYEARS][12];
	double modelGPP[NFLUXYEARS][12];

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
	std::map<std::pair<double, double>, FluxnetData> Fluxnetdata;
	std::vector<double> rain,Ta,swrad;
	std::vector<int> yr;
	Lamarque::NDepData ndep;

	/// Daily N deposition for current year
	double dndep[Date::MAX_YEAR_LENGTH];
};



#endif // LPJ_GUESS_Fluxnet_H
