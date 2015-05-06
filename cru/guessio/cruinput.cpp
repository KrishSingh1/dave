///////////////////////////////////////////////////////////////////////////////////////
/// \file cruinput.cpp
/// \brief LPJ-GUESS input module for CRU TS 3.0 data set
///
/// This input module reads in CRU climate data in a customised binary format.
/// The binary files contain CRU half-degree global historical climate data
/// for 1901-2006.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "cruinput.h"

#include "driver.h"
#include "parameters.h"
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>

// Switch to keep CO2 level at first historical year
bool fixedco2_hist=0;

REGISTER_INPUT_MODULE("cru", CRUInput)

// Anonymous namespace for variables and functions with file scope
namespace {

xtring file_cru;
xtring file_cru_misc;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double* mtemp, double* mprec, double* msun, double* mdtr,
					double* dtemp, double* dprec, double* dsun, double* ddtr) {
	interp_monthly_means_conserve(mtemp, dtemp);
	interp_monthly_totals_conserve(mprec, dprec, 0);
	interp_monthly_means_conserve(msun, dsun, 0, 100);
	interp_monthly_means_conserve(mdtr, ddtr, 0);
}

} // namespace


CRUInput::CRUInput()
	: searchradius(0),
	  climate_spatial_resolution(CRU_TS30::SPATIAL_RESOLUTION),
	  gridlist_spatial_resolution(DEFAULT_SPATIAL_RESOLUTION),
	  landcover_input_module(*this),
	  management_input_module(*this),
	  spinup_mtemp(NYEAR_SPINUP_DATA),
	  spinup_mprec(NYEAR_SPINUP_DATA),
	  spinup_msun(NYEAR_SPINUP_DATA),
	  spinup_mfrs(NYEAR_SPINUP_DATA),
	  spinup_mwet(NYEAR_SPINUP_DATA),
	  spinup_mdtr(NYEAR_SPINUP_DATA),
	  extended_mtemp(NYEAR_FUTURE_DATA),
	  extended_mprec(NYEAR_FUTURE_DATA),
	  extended_msun(NYEAR_FUTURE_DATA),
	  extended_mfrs(NYEAR_FUTURE_DATA),
	  extended_mwet(NYEAR_FUTURE_DATA),
	  extended_mdtr(NYEAR_FUTURE_DATA) {

	// Declare instruction file parameters

	declare_parameter("searchradius", &searchradius, 0, 100,
		"If specified, CRU data will be searched for in a circle");
}

void CRUInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files), and reads in the gridlist

	file_cru = param["file_cru"].str;
	file_cru_misc = param["file_cru_misc"].str;

	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]
	read_gridlist(gridlist, param["file_gridlist"].str);

	// Set the spatial resolution of the gridlist (needed when using data with different spatial resolution).
	gridlist_spatial_resolution = min(parse_gridlist_spatial_resolution(gridlist), gridlist_spatial_resolution);
	// Set gridlist_spatial_resolution here manually if needed (if higher than DEFAULT_SPATIAL_RESOLUTION or if 
	// gridlist too short to be sucessfully parsed for spatial resolution)

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	// Open landcover files
	landcover_input_module.init();
	// Open management files
	management_input_module.init();

	// Set the simulation period accoring to instruction file settings and/or climate time period
	set_simulation_years(this);

//	date.set_first_calendar_year(getfirsthistyear() - nyear_spinup);

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}


void CRUInput::get_monthly_ndep(int calendar_year,
                                double* mndrydep,
                                double* mnwetdep) {

	ndep.get_one_calendar_year(calendar_year, mndrydep, mnwetdep);
}

		
void CRUInput::adjust_raw_forcing_data(double lon,
                                       double lat,
                                       double hist_mtemp[NYEAR_HIST][12],
                                       double hist_mprec[NYEAR_HIST][12],
                                       double hist_msun[NYEAR_HIST][12]) {

	// The default (base class) implementation does nothing here.
}


bool CRUInput::getgridcell(Gridcell& gridcell) {

	/// See base class for documentation about this function's responsibilities
	
	int soilcode;
	int elevation;

	if(climate_spatial_resolution != gridlist_spatial_resolution && !search_for_centre_of_gridcell)
		fail("We must use a searchradius and search for centre of a gridcell when using different spatial resolution in input data\n");

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	static bool first_call = true;

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

		while(!gridfound) {

			if(gridlist.isobj) {

				lon = gridlist.getobj().lon;
				lat = gridlist.getobj().lat;

				gridfound = CRU_TS30::findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
											   hist_mtemp, hist_mprec, hist_msun);

				if (gridfound) // Get more historical CRU data for this grid cell
					gridfound = CRU_TS30::searchcru_misc(file_cru_misc, lon, lat, elevation, 
												   hist_mfrs, hist_mwet, hist_mdtr);

				if (run_landcover && gridfound) {
					LUerror = landcover_input_module.loadlandcover(gridcell, gridlist.getobj());
					if(!LUerror)
						LUerror = management_input_module.loadmanagement(gridcell, gridlist.getobj());
				}

				if(!gridfound || LUerror) {
					if(!gridfound)
						dprintf("\nError: could not find stand at (%g,%g) in climate data files\n", gridlist.getobj().lon,gridlist.getobj().lat);
					else if(LUerror)
						dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", gridlist.getobj().lon,gridlist.getobj().lat);
					gridfound = false;
					gridlist.nextobj();
				}
			}
			else return false;
		}


		// Give sub-classes a chance to modify the data
		adjust_raw_forcing_data(gridlist.getobj().lon,
								gridlist.getobj().lat,
								hist_mtemp, hist_mprec, hist_msun);

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

		// Build extended data sets
		int lastyear = NYEAR_HIST; // Year count from 1 here.
		int firstyear = lastyear - NYEAR_FUTURE_DATA + 1;
		extended_mtemp.extract_data(hist_mtemp, firstyear, lastyear);
		extended_mprec.extract_data(hist_mprec, firstyear, lastyear);
		extended_msun.extract_data(hist_msun, firstyear, lastyear);
		extended_mfrs.extract_data(hist_mfrs, firstyear, lastyear);
		extended_mwet.extract_data(hist_mwet, firstyear, lastyear);
		extended_mdtr.extract_data(hist_mdtr, firstyear, lastyear);

		// Detrend extended temperature data
		extended_mtemp.detrend_data(true);

		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
//		double offset = gridlist_spatial_resolution / 2.0;
//		gridcell.set_coordinates(gridlist.getobj().lon + offset, gridlist.getobj().lat + offset);	// Corrects previous errror
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Get nitrogen deposition data
		ndep.getndep(param["file_ndep"].str, lon, lat);

		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine

		gridcell.climate.instype = SUNSHINE;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype, soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

int CRUInput::getfirsthistyear_climate() {

	return FIRSTHISTYEAR;
}

int CRUInput::getnyear_hist_climate() {

	return NYEAR_HIST;
}

int CRUInput::getfirsthistyear() {

	return firsthistyear_sim;
}

int CRUInput::getnyear_hist() {

	return nyear_hist_sim;
}

bool CRUInput::getclimate(Gridcell& gridcell) {

	/// See base class for documentation about this function's responsibilities

	double progress;

	int calender_year = date.year - nyear_spinup + firsthistyear_sim;
	Climate& climate = gridcell.climate;

	if (date.day == 0) {

		// First day of year ...
	
		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		double mndrydep[12], mnwetdep[12];
		ndep.get_one_calendar_year(calender_year, mndrydep, mnwetdep);

		if(date.year < nyear_spinup + nyear_hist_sim) {

			if (calender_year < FIRSTHISTYEAR) {

				// During spinup period

				int m;
				double mtemp[12],mprec[12],msun[12];
				double mfrs[12],mwet[12],mdtr[12];

				for (m=0;m<12;m++) {
					mtemp[m] = spinup_mtemp[m];
					mprec[m] = spinup_mprec[m];
					msun[m] = spinup_msun[m];

					mfrs[m] = spinup_mfrs[m];
					mwet[m] = spinup_mwet[m];
					mdtr[m] = spinup_mdtr[m];
				}

				// Interpolate monthly spinup data to quasi-daily values
				interp_climate(mtemp,mprec,msun,mdtr,dtemp,dprec,dsun,ddtr);

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
			else if (calender_year < FIRSTHISTYEAR + NYEAR_HIST) {

				// Historical period

				int data_year = calender_year - FIRSTHISTYEAR;

				// Interpolate this year's monthly data to quasi-daily values
				interp_climate(hist_mtemp[data_year], hist_mprec[data_year], hist_msun[data_year],
						   hist_mdtr[data_year], dtemp, dprec, dsun, ddtr);

				// Only recalculate precipitation values using weather generator
				// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
				if (ifrainonwetdaysonly) { 
					// (from Dieter Gerten 021121)
					prdaily(hist_mprec[data_year], dprec, hist_mwet[data_year], gridcell.seed);
				}
			}
			else {

				// Extended period

				int m;
				double mtemp[12],mprec[12],msun[12];
				double mfrs[12],mwet[12],mdtr[12];

				for (m=0;m<12;m++) {
					mtemp[m] = extended_mtemp[m];
					mprec[m] = extended_mprec[m];
					msun[m] = extended_msun[m];

					mfrs[m] = extended_mfrs[m];
					mwet[m] = extended_mwet[m];
					mdtr[m] = extended_mdtr[m];
				}

				// Interpolate monthly spinup data to quasi-daily values
				interp_climate(mtemp, mprec, msun, mdtr, dtemp, dprec, dsun, ddtr);

				// Only recalculate precipitation values using weather generator
				// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
				if (ifrainonwetdaysonly) { 
					// (from Dieter Gerten 021121)
					prdaily(mprec, dprec, mwet, gridcell.seed);
				}

				extended_mtemp.nextyear();
				extended_mprec.nextyear();
				extended_msun.nextyear();
				extended_mfrs.nextyear();
				extended_mwet.nextyear();
				extended_mdtr.nextyear();

				if(calender_year == FIRSTHISTYEAR + NYEAR_HIST)
					dprintf("Last %d years of CRU climate data used from year %d and onwards\n", NYEAR_FUTURE_DATA, calender_year);
			}
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}

		// Distribute N deposition
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
	}

	if(date.day == 0) {

		for(int m=0; m<12; m++)
			climate.mpet_year[m] = 0.0;
	}
	climate.mpet_year[date.month] += climate.eet * PRIESTLEY_TAYLOR;

	// Send environmental values for today to framework
	if(fixedco2_hist)
		climate.co2 = co2[firsthistyear_sim];
	else
		climate.co2 = co2[calender_year];

	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dsun[date.day];

	// Nitrogen deposition
	climate.dndep = dndep[date.day];

	// bvoc
	if(ifbvoc){
	  climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+getnyear_hist())
				+date.year)/(double)(gridlist.nobj*(nyear_spinup+getnyear_hist()));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}

bool CRUInput::getsoil(Gridcell& gridcell, const int soilmap_index){
	return true;
}

CRUInput::~CRUInput() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
}

// Copies gridlist to calling function's gridlist
void CRUInput::getgridlist(ListArray_id<Coord>& outlist) {

	gridlist.firstobj();
	while(gridlist.isobj) {
		Coord& c = outlist.createobj();
		c.lon = gridlist.getobj().lon;
		c.lat = gridlist.getobj().lat;
		c.descrip = gridlist.getobj().descrip;
		gridlist.nextobj();
	}
	gridlist.firstobj();
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Lamarque, J.-F., Kyle, G. P., Meinshausen, M., Riahi, K., Smith, S. J., Van Vuuren, 
//   D. P., Conley, A. J. & Vitt, F. 2011. Global and regional evolution of short-lived
//   radiatively-active gases and aerosols in the Representative Concentration Pathways. 
//   Climatic Change, 109, 191-212.
// Nakai, T., Sumida, A., Kodama, Y., Hara, T., Ohta, T. (2010). A comparison between
//   various definitions of forest stand height and aerodynamic canopy height.
//   Agricultural and Forest Meteorology, 150(9), 1225-1233
