///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Input and output modules for the fluxnet benchmarks
/// \this module is most of all a copy of CRUInput and do use CRU data for its climate.
/// \The CRU climate is then biascorrected with monthly fluxnet climate.
///
/// \author Niklas Boke Olén and Adrian Gustafson
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "guess.h"
#include "driver.h"
#include "guessstring.h"
#include "fluxnet.h"
#include <fstream>
#include <sstream>

using namespace std;

REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)

namespace {

	xtring file_cru;
	xtring file_cru_misc;

	/// Interpolates monthly data to quasi-daily values.
	void interp_climate(double* mtemp, double* mprec, double* msun,
		double* dtemp, double* dprec, double* dsun) {
		interp_monthly_means_conserve(mtemp, dtemp);
		interp_monthly_totals_conserve(mprec, dprec, 0);
		interp_monthly_means_conserve(msun, dsun, 0);
	}
}

FluxnetInput::FluxnetInput()
	: searchradius(0),
	spinup_mtemp(NYEAR_SPINUP_DATA),
	spinup_mprec(NYEAR_SPINUP_DATA),
	spinup_msun(NYEAR_SPINUP_DATA),
	spinup_mfrs(NYEAR_SPINUP_DATA),
	spinup_mwet(NYEAR_SPINUP_DATA),
	spinup_mdtr(NYEAR_SPINUP_DATA) {

	// Declare instruction file parameters

	declare_parameter("searchradius", &searchradius, 0, 100,
		"If specified, CRU data will be searched for in a circle");
}


void FluxnetInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files), and reads in the gridlist

	//
	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon, dlat;
	bool eof = false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist = param["file_gridlist"].str;

	FILE* in_grid = fopen(file_gridlist, "r");
	if (!in_grid) fail("initio: could not open %s for input", (char*)file_gridlist);

	file_cru = param["file_cru"].str;
	file_cru_misc = param["file_cru_misc"].str;

	gridlist.killall();
	first_call = true;

	while (!eof) {

		// Read next record in file
		eof = !readfor(in_grid, "f,f,a#", &dlon, &dlat, &descrip);

		if (!eof && !(dlon == 0.0 && dlat == 0.0)) { // ignore blank lines at end (if any)
			Coord& c = gridlist.createobj(); // add new coordinate to grid list

			c.lon = dlon;
			c.lat = dlat;
			c.descrip = descrip;
		}
	}


	fclose(in_grid);

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	// Open landcover files
	landcover_input.init();
	// Open management files
	management_input.init();

	date.set_first_calendar_year(FIRSTHISTYEAR - nyear_spinup);
	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);


};

void FluxnetInput::adjust_raw_forcing_data(double hist_mtemp[NYEAR_HIST][12],
			double hist_mprec[NYEAR_HIST][12], double hist_msun[NYEAR_HIST][12],
			double fluxnet_temp[12], double fluxnet_prec[12], double fluxnet_rad[12]) {

};

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
	/*
	This method reads in monthly CRU-data and then fluxnet station data. Using the 
	latter it bias corrects the CRU-data to the station using an anomaly appraoch. 
	*/
	
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
		int nyear, first_year, last_year;
		bool firstcall = true;
		std::vector<double> rain, tair, swrad;

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

			int year, month;
			double temp, prec, rad;

			std::istringstream iss(line);

			if (iss >> year >> month >> temp >> rad >> prec) {
				if (firstcall) {
					first_year = year;
				}

				tair.push_back(temp);
				swrad.push_back(rad);
				rain.push_back(prec);
				
				last_year = year;
				firstcall = false;
			}
		}

		ifs.close();

		// Check so that input contains full years
		std::vector<double>::size_type days = rain.size();
		if (days % 12) {
			dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
			return false;
		}

		double mtemp_fluxnet[12], mprec_fluxnet[12], mrad_fluxnet[12];

		for (int i = 0; i < 12; i++) {
			mrad_fluxnet[i] = 0.0;
			mprec_fluxnet[i] = 0.0;
			mtemp_fluxnet[i] = 0.0;
		}

		nyear = last_year - first_year + 1;
		
		// Make monthly averages from the input data
		int idx = 0;
		for (int i = 0; i < tair.size(); i++) {
			mtemp_fluxnet[idx] += tair[i] / nyear;
			mrad_fluxnet[idx] += swrad[i] / nyear;
			mprec_fluxnet[idx] += rain[i] / nyear;
			idx++;
			if (idx % 12 == 0) { idx = 0; }
		}

		// Give sub-classes a chance to modify the data
		adjust_raw_forcing_data(hist_mtemp, hist_mprec, hist_msun, mtemp_fluxnet,
			mprec_fluxnet, mrad_fluxnet);

		/* CALCULATE AND APPLY ANOMALIES*/
		double cru_temp_mean[12], cru_rain_mean[12], cru_rad_mean[12];
		double temp_anom[12], rain_anom[12], rad_anom[12];

		// initialise
		for (int i = 0; i < 12; i++) {
			cru_rad_mean[i] = 0.0;
			cru_temp_mean[i] = 0.0;
			cru_rain_mean[i] = 0.0;

			rad_anom[i] = 0.0;
			temp_anom[i] = 0.0;
			rain_anom[i] = 0.0;
		}

		for (int yr = first_year - FIRSTHISTYEAR; yr < last_year - FIRSTHISTYEAR; yr++) {
			for (int m = 0; m < 12; m++) {
				cru_temp_mean[m] += hist_mtemp[yr][m] / nyear;
				cru_rad_mean[m] += hist_msun[yr][m] / nyear;
				cru_rain_mean[m] += hist_mprec[yr][m] / nyear;
			};
		};

		for (int i = 0; i < 12; i++) {
			rad_anom[i] =  mrad_fluxnet[i] - cru_rad_mean[i];
			temp_anom[i] =  mtemp_fluxnet[i] - cru_temp_mean[i];
			if (cru_rain_mean[i] == 0.0) {
				// Protect against potential division by zero
				// can happen in very arid areas
				rain_anom[i] = 0.0;
			}
			else {
				rain_anom[i] = mprec_fluxnet[i] / cru_rain_mean[i];
			}
			
		}

		for (int yr = 0; yr < NYEAR_HIST; yr++) {
			for (int m = 0; m < 12; m++) {
				hist_mtemp[yr][m] += temp_anom[m];
				hist_mprec[yr][m] *= rain_anom[m];
				hist_msun[yr][m] += rad_anom[m];
				hist_msun[yr][m] = max(hist_msun[yr][m], 0.0); // Radiation cannot be negative
			}
		}


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

void FluxnetInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	// See base class for documentation about this function's responsibilities

	double progress;

	Climate& climate = gridcell.climate;

	if (date.day == 0) {

		// First day of year ...

		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		double mndrydep[12], mnwetdep[12];
		ndep.get_one_calendar_year(date.year - nyear_spinup + FIRSTHISTYEAR,
			mndrydep, mnwetdep);

		if (date.year < nyear_spinup) {

			// During spinup period

			if (date.year == state_year && restart) {

				int year_offset = state_year % NYEAR_SPINUP_DATA;

				for (int y = 0; y<year_offset; y++) {
					spinup_mtemp.nextyear();
					spinup_mprec.nextyear();
					spinup_msun.nextyear();
					spinup_mfrs.nextyear();
					spinup_mwet.nextyear();
					spinup_mdtr.nextyear();
				}
			}

			int m;
			double mtemp[12], mprec[12], msun[12];
			double mfrs[12], mwet[12], mdtr[12];

			for (m = 0; m<12; m++) {
				mtemp[m] = spinup_mtemp[m];
				mprec[m] = spinup_mprec[m];
				msun[m] = spinup_msun[m];

				mfrs[m] = spinup_mfrs[m];
				mwet[m] = spinup_mwet[m];
				mdtr[m] = spinup_mdtr[m];
			}

			// Interpolate monthly spinup data to quasi-daily values
			interp_climate(mtemp, mprec, msun, dtemp, dprec, dsun);

			// Only recalculate precipitation values using weather generator
			// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) {
				// (from Dieter Gerten 021121)
				prdaily(mprec, dprec, mwet, gridcell.seed);
			}

			spinup_mtemp.nextyear();
			spinup_mprec.nextyear();
			spinup_msun.nextyear();

			spinup_mfrs.nextyear();
			spinup_mwet.nextyear();
			spinup_mdtr.nextyear();

		}
		else if (date.year < nyear_spinup + NYEAR_HIST) {

			// Historical period

			// Interpolate this year's monthly data to quasi-daily values
			interp_climate(hist_mtemp[date.year - nyear_spinup],
				hist_mprec[date.year - nyear_spinup], hist_msun[date.year - nyear_spinup],
				dtemp, dprec, dsun);

			// Only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) {
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year - nyear_spinup], dprec, hist_mwet[date.year - nyear_spinup], gridcell.seed);
			}
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}

		// Distribute N deposition
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
	}

	// Send environmental values for today to framework

	climate.co2 = co2[FIRSTHISTYEAR + date.year - nyear_spinup];

	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dsun[date.day];

	// Nitrogen deposition
	climate.dndep = dndep[date.day];

	// bvoc
	if (ifbvoc) {
		climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress() >= 1.0) {
			progress = (double)(gridlist.getobj().id*(nyear_spinup + NYEAR_HIST)
				+ date.year) / (double)(gridlist.nobj*(nyear_spinup + NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n", (int)(progress*100.0),
				tprogress.elapsed.str, tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


FluxnetInput::~FluxnetInput() {
	// Performs deconstruction of FluxnetInput object
}