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


CRUInput::CRUInput(Input& in)
	: searchradius(0),
	  gridlist(in.gridlist),
	  input(in),
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
	// Initialises input

	file_cru = param["file_cru"].str;
	file_cru_misc = param["file_cru_misc"].str;
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

	bool gridfound;

	double lon = gridlist.getobj().lon;
	double lat = gridlist.getobj().lat;
	gridfound = CRU_TS30::findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
	                               hist_mtemp, hist_mprec, hist_msun);

	if (gridfound) // Get more historical CRU data for this grid cell
		gridfound = CRU_TS30::searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);
	if(gridfound) {

		// Give sub-classes a chance to modify the data
		adjust_raw_forcing_data(gridlist.getobj().lon,
								gridlist.getobj().lat,
								hist_mtemp, hist_mprec, hist_msun);

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

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

		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine			
		gridcell.climate.instype = SUNSHINE;

#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
		// Get nitrogen deposition data
		ndep.getndep(param["file_ndep"].str, lon, lat);
#endif

#ifdef SOIL_INPUT_IN_CLIMATE_MODULE
		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype, soilcode);
#endif
	}
	return gridfound;
}

int CRUInput::getfirsthistyear() {

	return FIRSTHISTYEAR;
}

int CRUInput::getnyear_hist() {

	return NYEAR_HIST;
}

bool CRUInput::getclimate(Gridcell& gridcell) {

	/// See base class for documentation about this function's responsibilities

	int calender_year = date.year - nyear_spinup + input.firsthistyear;
	Climate& climate = gridcell.climate;

	if (date.day == 0) {

		// First day of year ...
	
#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		double mndrydep[12], mnwetdep[12];
		ndep.get_one_calendar_year(calender_year, mndrydep, mnwetdep);
#endif
		if (calender_year < FIRSTHISTYEAR) {

			// During spinup period

			int m;
			double mtemp[12],mprec[12],msun[12];
			double mfrs[12],mwet[12],mdtr[12];

			for (m=0;m<12;m++) {
				mtemp[m] = spinup_mtemp[m];
				mprec[m] = spinup_mprec[m];
				climate.mtemp_year[m] = spinup_mtemp[m];	
				climate.mprec_year[m] = spinup_mprec[m];	
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

			// save climate date for calculation of seasonality
			for(int m=0;m<12;m++) {
				climate.mtemp_year[m] = hist_mtemp[data_year][m];
				climate.mprec_year[m] = hist_mprec[data_year][m];
			}
		}
		else if(date.year < nyear_spinup + input.nyear_hist) {

			// Extended period

			int m;
			double mtemp[12],mprec[12],msun[12];
			double mfrs[12],mwet[12],mdtr[12];

			for (m=0;m<12;m++) {
				mtemp[m] = extended_mtemp[m];
				mprec[m] = extended_mprec[m];
				climate.mtemp_year[m] = extended_mtemp[m];	
				climate.mprec_year[m] = extended_mprec[m];	
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
		else {
			// Return false if last year was the last for the simulation
			return false;
		}

#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
		// Distribute N deposition
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
#endif
	}

	if(date.day == 0) {
		climate.aprec = 0.0;

		for(int m=0; m<12; m++) {
			climate.aprec += climate.mprec_year[m];
			climate.mpet_year[m] = 0.0;
		}
	}
	climate.mpet_year[date.month] += climate.eet * PRIESTLEY_TAYLOR;

	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dsun[date.day];

#ifdef NDEP_INPUT_IN_CLIMATE_MODULE
	// Nitrogen deposition
	climate.dndep = dndep[date.day];
#endif

	// bvoc
	if(ifbvoc){
	  climate.dtr = ddtr[date.day];
	}

	return true;
}

bool CRUInput::getsoil(Gridcell& gridcell, const int soilmap_index){
	return true;
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
