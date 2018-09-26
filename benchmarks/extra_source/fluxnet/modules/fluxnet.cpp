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
#include "fluxnet.h"
#include <fstream>
#include <sstream>

using namespace std;

namespace {

	xtring file_cru;
	xtring file_cru_misc;

}

REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)

void FluxnetInput::init() {
	CRUInput::init(); 

	file_cru = param["file_cru"].str;
	file_cru_misc = param["file_cru_misc"].str;

	for (int i = 0; i < 12; i++) {
		monthly_fluxnet_rad[i] = 0.0;
		monthly_fluxnet_rain[i] = 0.0;
		monthly_fluxnet_temp[i] = 0.0;
	}
};

void FluxnetInput::adjust_raw_forcing_data(double hist_mtemp[NYEAR_HIST][12],
			double hist_mprec[NYEAR_HIST][12], double hist_msun[NYEAR_HIST][12],
			double fluxnet_temp[12], double fluxnet_prec[12], double fluxnet_rad[12]) {
	
	double cru_temp_mean[12], cru_rain_mean[12], cru_rad_mean[12];
	double cru_temp_anom[12], cru_rain_anom[12], cru_rad_anom[12];

    // initialise
	for (int i = 0; i < 12; i++) {
		cru_rad_mean[i] = 0.0;
		cru_temp_mean[i] = 0.0;
		cru_rain_mean[i] = 0.0;

		cru_rad_anom[i] = 0.0;
		cru_temp_anom[i] = 0.0;
		cru_rain_anom[i] = 0.0;
	}

	for (int yr = first_year - FIRSTHISTYEAR; yr < last_year - FIRSTHISTYEAR; yr++) {
		for (int m = 0; m < 12; m++) {
			cru_temp_mean[m] += hist_mtemp[yr][m] / nyear;
			cru_rad_mean[m] += hist_msun[yr][m] / nyear;
			cru_rain_mean[m] += hist_mprec[yr][m] / nyear;
		};
	};

	for (int i = 0; i < 12; i++) {
		cru_rad_anom[i] = fluxnet_rad[i] / cru_rad_mean[i];
		cru_temp_anom[i] = fluxnet_temp[i] - cru_temp_mean[i];
		cru_rain_anom[i] = fluxnet_prec[i] / cru_rain_mean[i];
	}

	for (int yr = 0; yr < NYEAR_HIST; yr++) {
		for (int m = 0; m < 12; m++) {
			hist_mtemp[yr][m] += cru_temp_anom[m];
            hist_mprec[yr][m] += min(cru_rain_anom[m], 0.0);
			hist_msun[yr][m] *= cru_rad_anom[m];
		}
	}


};

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
	
	// See base class for documentation about this function's responsibilities

	int soilcode;
	int elevation;

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		bool gridfound = false;
		bool LUerror = false;
		double lon;
		double lat;

		while (!gridfound) {

			if (gridlist.isobj) {

				lon = gridlist.getobj().lon;
				lat = gridlist.getobj().lat;
				gridfound = CRU_TS30::findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
					hist_mtemp, hist_mprec, hist_msun);

				if (gridfound) // Get more historical CRU data for this grid cell
					gridfound = CRU_TS30::searchcru_misc(file_cru_misc, lon, lat, elevation,
						hist_mfrs, hist_mwet, hist_mdtr);

				if (run_landcover && gridfound) {
					LUerror = landcover_input.loadlandcover(lon, lat);
					if (!LUerror)
						LUerror = management_input.loadmanagement(lon, lat);
				}

				if (!gridfound || LUerror) {
					if (!gridfound)
						dprintf("\nError: could not find stand at (%g,%g) in climate data files\n", gridlist.getobj().lon, gridlist.getobj().lat);
					else if (LUerror)
						dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", gridlist.getobj().lon, gridlist.getobj().lat);
					gridfound = false;
					gridlist.nextobj();
				}
			}
			else return false;
		}

		///FLUXNET UNIQUE CODE:

		tair.clear();
		rain.clear();
		swrad.clear();

		// Fetch the fluxnetdata from the gridlist description
		xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip + ".csv";
		std::ifstream ifs(fluxfile, std::ifstream::in);

		if (!ifs.good()) {
			dprintf("FluxnetInput::getgridcell: could not open %s for input\n", (char*)fluxfile);
			return false;
		}

		std::string line;
		std::istringstream iss(line);

		while (getline(ifs, line)) {

			int month;
			double temp, prec, rad;

			std::istringstream iss(line);

			if (iss >> last_year >> month >> temp >> rad >> prec) {
				if (tair.empty()) {
					dailyoutput_firstyear = first_year = last_year;
				}
				tair.push_back(temp);
				swrad.push_back(rad);
				rain.push_back(prec);
			}
		}

		ifs.close();

		// Check so that input contains full years
		std::vector<double>::size_type days = rain.size();
		if (days % 12) {
			dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
			return false;
		}

		nyear = last_year - first_year + 1;

		// Make monthly averages from the input data
		int idx = 0;
		for (int i = 0; i < tair.size(); i++) {
			monthly_fluxnet_temp[idx] += tair[i] / nyear;
			monthly_fluxnet_rad[idx] += swrad[i] / nyear;
			monthly_fluxnet_rain[idx] += rain[i] / nyear;
			idx++;
			if (idx % 12 == 0) { idx = 0; }
		}

		// Give sub-classes a chance to modify the data
		adjust_raw_forcing_data(hist_mtemp, hist_mprec, hist_msun, monthly_fluxnet_temp, 
			monthly_fluxnet_rain, monthly_fluxnet_rad);

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);

		// We wont detrend dtr for now. Partly because dtr is at the moment only
		// used for BVOC, so what happens during the spinup is not affecting
		// results in the period thereafter, and partly because the detrending
		// can give negative dtr values.
		//spinup_mdtr.detrend_data();


		dprintf("\nCommencing simulation for stand at (%g,%g)", gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip != "") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");

		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Get nitrogen deposition data. 
		/* Since the historic data set does not reach decade 2010-2019,
		* we need to use the RCP data for the last decade. */
		ndep.getndep(param["file_ndep"].str, lon, lat, Lamarque::RCP60);

		// The insolation data will be sent (in function getclimate, below)
		// as incoming shortwave radiation, averages are over 24 hours

		gridcell.climate.instype = SWRAD_TS;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype, soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)

		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
	
}

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	if (!CRUInput::getclimate(gridcell)) {
		return false;
	}
	return true;
}

