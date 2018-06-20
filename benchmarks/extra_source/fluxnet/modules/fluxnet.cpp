///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Input and output modules for the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "fluxnet.h"
#include "driver.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>

namespace {
    xtring file_cru;
    xtring file_cru_misc;
}


void FluxnetInput::init() {
    CRUInput::init();
    
    file_cru=param["file_cru"].str;
    file_cru_misc=param["file_cru_misc"].str;
}

bool FluxnetInput::get_fluxnet_data_from_file() {
	
	std::string line;
	
	tair.clear();
	rain.clear();
	swrad.clear();
	
	xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip.printable() + ".csv";
	std::ifstream ifs(fluxfile, std::ifstream::in);

	if (!ifs.good()) {
		dprintf("FluxnetInput::getgridcell: could not open %s for input\n", (char*)fluxfile);
		return false;
	}
	
	while (getline(ifs, line)) {
		
		std::string monthday;
		double file_temp, file_prec, file_insol;
		int year;
		
		std::istringstream iss(line);
		if (iss >> year >> monthday >> file_temp >> file_insol >> file_prec) {
			if (tair.empty()) {
				first_fluxnet_year = year;
			}

			tair.push_back(file_temp);
			swrad.push_back(file_insol);
			rain.push_back(file_prec);
		}
	}
	
	ifs.close();
	
	std::vector<double>::size_type days = rain.size();
	if (days % 365) {
		dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
		return false;
	}
	return true;
}

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
	if (!CRUInput::getgridcell(gridcell)){
		return false;
	}
	
	return FluxnetInput::get_fluxnet_data_from_file();

}

    
    

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	
	if (!CRUInput::getclimate(gridcell)) {
		return false;
	};
	
	int year = date.get_calendar_year();
	if (year >= first_fluxnet_year) {
		// Overwrite/extend climate data with site values for the period
		int id = (year - first_fluxnet_year) * 365 + date.day;

		
		if (date.day == 0) {
			// Put this years daily precipitation into an array
			int ct=0;
			for (int i=id; i < id+365; i++) {
				drain[ct] = rain[i];
				ct++;
			}
			double mndrydep[12], mnwetdep[12];
			ndep.get_one_calendar_year(year, mndrydep, mnwetdep);
			distribute_ndep(mndrydep, mnwetdep, drain, dndep_fluxnet);
		}
		
		gridcell.climate.prec = rain[id];
		gridcell.climate.temp = tair[id];
		gridcell.climate.insol = swrad[id];
		gridcell.climate.dndep = dndep_fluxnet[date.day];
	}
	return true;
}


REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)
