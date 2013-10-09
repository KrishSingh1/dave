///////////////////////////////////////////////////////////////////////////////////////
/// \file lamarquendep.cpp
/// \brief Functionality for reading the Lamarque Nitrogen deposition data set
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "lamarquendep.h"
#include "shell.h"
#include <stdio.h>
#include <string>

// header file for reading binary data archive of global nitrogen deposition
#include "GlobalNitrogenDeposition.h"

namespace Lamarque {

const double convert = 1e-7;				// converting from gN ha-1 to kgN m-2

NDepData::NDepData() {
	set_to_pre_industrial();
}

void NDepData::getndep(const char* file_ndep,
					   double lon, double lat) {

	if (std::string(file_ndep) == "") {
		set_to_pre_industrial();
	}
	else {
		GlobalNitrogenDepositionArchive ark;
		if (!ark.open(file_ndep)) {
			fail("Could not open %s for input", (char*)file_ndep);
		}

		GlobalNitrogenDeposition rec;
		rec.longitude = lon;
		rec.latitude = lat;

		if (!ark.getindex(rec)) {
			ark.close();
			fail("Grid cell not found in %s", (char*)file_ndep);
		}

		// Found the record, get the values
		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = rec.NHxDry[y*12+m] * convert;
				NHxWetDep[y][m] = rec.NHxWet[y*12+m] * convert;
				NOyDryDep[y][m] = rec.NOyDry[y*12+m] * convert;
				NOyWetDep[y][m] = rec.NOyWet[y*12+m] * convert;
			}
		}
		ark.close();
	}
}

void NDepData::get_one_calendar_year(int calendar_year,
									 double mndrydep[12],
									 double mnwetdep[12]) {
	int ndep_year = 0;

	if (calendar_year >= Lamarque::FIRSTHISTYEARNDEP) {
		ndep_year = (int)((calendar_year - Lamarque::FIRSTHISTYEARNDEP)/10);
	}

	if (ndep_year >= NYEAR_HISTNDEP) {
		fail("Tried to get ndep for year %d (not included in Lamarque ndep data set)",
		     calendar_year);
	}

	for (int m = 0; m < 12; m++) {
		mndrydep[m] = NHxDryDep[ndep_year][m] + NOyDryDep[ndep_year][m];
		
		mnwetdep[m] = NHxWetDep[ndep_year][m] + NOyWetDep[ndep_year][m];
	}	
}

void NDepData::set_to_pre_industrial() {
	// pre-industrial total nitrogen depostion set to 2 kgN/ha/year [kgN m-2]
	double dailyndep = 2000.0 / (4 * 365) * convert;

	for (int y=0; y<NYEAR_HISTNDEP; y++) {
		for (int m=0; m<12; m++) {
			NHxDryDep[y][m] = dailyndep;
			NHxWetDep[y][m] = dailyndep;
			NOyDryDep[y][m] = dailyndep;
			NOyWetDep[y][m] = dailyndep;
		}
	}
}

}
