///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Input and output modules for the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "cruinput.h"
#include "guess.h"
#include "driver.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>

class FluxnetInput : public CRUInput {
public:
bool getgridcell(Gridcell& gridcell) {
	if (!CRUInput::getgridcell(gridcell)) {
		return false;
	}
	tair.clear();
	rain.clear();
	swrad.clear();

	xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip + ".csv";
	std::ifstream ifs(fluxfile, std::ifstream::in);

	if (!ifs.good()) {
		dprintf("FluxnetInput::getgridcell: could not open %s for input\n", (char*)fluxfile);
		return false;
	}

	std::string line;
	getline(ifs, line);
	std::istringstream iss(line);
	double lon, lat;
	iss >> lon >> lat;

	while (getline(ifs, line)) {

		std::string monthday;
		double temp, prec, insol;

		std::istringstream iss(line);

		if (iss >> end_year >> monthday >> temp >> insol >> prec) {
			if (tair.empty()) {
				dailyoutput_firstyear = end_year;
			}
			tair.push_back(temp);
			swrad.push_back(insol);
			rain.push_back(prec);
		}
	}

	ifs.close();

	std::vector<double>::size_type days = rain.size();
	if (days % 365) {
		dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
		return false;
	}

	gridcell.set_coordinates(lon, lat);

	return true;
}

bool getclimate(Gridcell& gridcell) {
	CRUInput::getclimate(gridcell);
	
	int year = date.get_calendar_year();
	if (year < dailyoutput_firstyear) {
		return true;
	} else if (year > end_year) {
		return false;
	}

	// Overwrite / extend climate data with site values for the period
	int id = (year - dailyoutput_firstyear) * 365 + date.day;

	if (date.day == 0) {
		double mndrydep[12], mnwetdep[12];
		ndep.get_one_calendar_year(year, mndrydep, mnwetdep);
		distribute_ndep(mndrydep, mnwetdep, &rain[id], dndep);
	}

	gridcell.climate.prec = rain[id];
	gridcell.climate.temp = tair[id];
	gridcell.climate.insol = swrad[id];
	gridcell.climate.dndep = dndep[date.day];
	return true;
}

private:
	std::vector<double> rain, tair, swrad;
	int end_year;
};

REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)
