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

void getndep(const char* file_ndep,
             double lon, double lat, 
             double NHxDryDep[NYEAR_HISTNDEP][12],
             double NHxWetDep[NYEAR_HISTNDEP][12],
             double NOyDryDep[NYEAR_HISTNDEP][12],
             double NOyWetDep[NYEAR_HISTNDEP][12]) {

	const double convert = 1e-7;				// converting from gN ha-1 to kgN m-2

	if (std::string(file_ndep) == "") {

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

}
