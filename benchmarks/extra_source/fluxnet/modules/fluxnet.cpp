///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Input and output modules for the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "fluxnet.h"
#include "guess.h"
#include "driver.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>

REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)
REGISTER_OUTPUT_MODULE("fluxnet", FluxnetOutput)

using namespace GuessOutput;

namespace {
int start_year = 0, end_year = 0;
}

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
	if (!CRUInput::getgridcell(gridcell)) {
		return false;
	}
	tair.clear();
	rain.clear();
	swrad.clear();

	xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip + ".csv";
	std::ifstream ifs(fluxfile, std::ifstream::in);

	if (!ifs.good()) {
		dprintf("FluxnetInput::getgridcell: could not open %s for input", (char*)fluxfile);
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
				start_year = end_year;
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

	ndep.getndep(param["file_ndep"].str, gridcell.get_lon(), gridcell.get_lat(), Lamarque::RCP60);

	gridcell.set_coordinates(lon, lat);

	return true;
}

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	CRUInput::getclimate(gridcell);

	int year = date.get_calendar_year();
	if (year < start_year) {
		return true;
	} else if (year > end_year) {
		return false;
	}

	// Overwrite / extend climate data with site values for the period
	int id = (year - start_year) * 365 + date.day;

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

/// OUTPUT MODULE FOR FLUXNET
FluxnetOutput::FluxnetOutput() {
	declare_parameter("file_fluxnetdaily", &file_fluxnetdaily, 300, "FLUXNET daily output file");
	declare_parameter("file_fluxnetclim", &file_fluxnetclim, 300, "FLUXNET daily Climate output file");
}

void FluxnetOutput::init() {

	ColumnDescriptors fluxnet_columns;
	fluxnet_columns += ColumnDescriptor("GPP", 14, 6);
	fluxnet_columns += ColumnDescriptor("NEE", 14, 6);
	create_output_table(out_fluxnetdaily, file_fluxnetdaily, fluxnet_columns);

	ColumnDescriptors fluxnetclim_columns;
	fluxnetclim_columns += ColumnDescriptor("temp", 11, 3);
	fluxnetclim_columns += ColumnDescriptor("prec", 11, 3);
	fluxnetclim_columns += ColumnDescriptor("insol", 11, 3);
	create_output_table(out_fluxnetclim, file_fluxnetclim, fluxnetclim_columns);
}

void FluxnetOutput::outdaily(Gridcell& gridcell) {
	if (date.get_calendar_year() < start_year) {
		return;
	}
	OutputRows out(output_channel, gridcell.get_lon(), gridcell.get_lat(),
			date.get_calendar_year(), date.day);

	out.add_value(out_fluxnetclim, gridcell.climate.temp);
	out.add_value(out_fluxnetclim, gridcell.climate.prec);
	out.add_value(out_fluxnetclim, gridcell.climate.insol);

	double dgpp = 0, dra = 0, drh = 0;

	Gridcell::iterator gc_itr = gridcell.begin();

	// Loop through Stands
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();

		//Loop through Patches
		while (stand.isobj) {
			Patch& patch = stand.getobj();

			double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

			dgpp +=-patch.fluxes.get_daily_flux(Fluxes::GPP,date.day)*to_gridcell_average;
			dra  +=-patch.fluxes.get_daily_flux(Fluxes::RA,date.day)*to_gridcell_average;
			drh  += patch.fluxes.get_daily_flux(Fluxes::SOILC,date.day)*to_gridcell_average;

			stand.nextobj();
		} // patch loop
		++gc_itr;
	} // stand loop

	out.add_value(out_fluxnetdaily, dgpp);
	// daily NEE do not include fire, establishment as monthly NEE
	out.add_value(out_fluxnetdaily, dgpp - dra + drh);
}
