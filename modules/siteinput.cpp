///////////////////////////////////////////////////////////////////////////////////////
/// \file SiteInput.cpp
/// \brief LPJ-GUESS input module for CRU-NCEP binary data set.
///
/// This input module reads in CRU-NCEP climate data in a customised binary format.
/// The binary files contain CRU-NCEP half-degree global historical climate data
/// for 1901-2015.
///
/// \author Ben Smith
/// $Date: 2020-03-03 16:31:01 +0100 (tis, 03 mar 2020) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "siteinput.h"
#include "weathergen.h"
#include "driver.h"
#include "parameters.h"
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>
#include "indata.h"

REGISTER_INPUT_MODULE("site", SiteInput)

// Anonymous namespace for variables and functions with file scope
namespace {

xtring file_cru;
xtring file_cru_misc;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double* mtemp, double* mprec, double* msun, double* mdtr,
					double* dtemp, double* dprec, double* dsun, double* ddtr) {
	interp_monthly_means_conserve(mtemp, dtemp);
	interp_monthly_totals_conserve(mprec, dprec, 0);
	interp_monthly_means_conserve(msun, dsun, 0);
	//interp_monthly_means_conserve(mdtr, ddtr, 0);
}

/// Transfers coordinates from SiteInput::Coord to InData::Coord
//InData::Coord SiteInput::GetLonLat(SiteInput::Coord coord) {
InData::Coord GetLonLat(Coord coord) {
	InData::Coord lonlat;
	lonlat.lon = coord.lon;
	lonlat.lat = coord.lat;

	return lonlat;
}

} // namespace


std::vector<std::pair<double, double> > SiteInput::translate_gridlist_to_coord(ListArray_id<Coord>& gridlist) {
	gridlist.firstobj();
	std::vector<std::pair<double, double> > output;
	while (gridlist.isobj) {
		Coord& c = gridlist.getobj();
		double center_lon = floor(c.lon * (1/ soilinput.STEP)) * soilinput.STEP + soilinput.STEP / 2.0;
		double center_lat = floor(c.lat * (1 / soilinput.STEP)) * soilinput.STEP + soilinput.STEP / 2.0;
		output.push_back(std::make_pair(center_lon, center_lat));
		gridlist.nextobj();
	}
	return output;
}




SiteInput::SiteInput()
	: searchradius(0),
	  spinup_mtemp(NYEAR_SPINUP_DATA),
	  spinup_mprec(NYEAR_SPINUP_DATA),
	  spinup_msun(NYEAR_SPINUP_DATA),
	  spinup_mfrs(NYEAR_SPINUP_DATA),
	  spinup_mwet(NYEAR_SPINUP_DATA),
	  spinup_mdtr(NYEAR_SPINUP_DATA),
	  spinup_mwind(NYEAR_SPINUP_DATA),
	  spinup_mrhum(NYEAR_SPINUP_DATA) {

	// Declare instruction file parameters

	declare_parameter("searchradius", &searchradius, 0, 100,
		"If specified, CRU data will be searched for in a circle");
}

void SiteInput::get_daily_climate(Gridcell& gridcell,xtring variable, double data[], int year) {
/* 	int i, year;
	int first_historic_year = FIRSTHISTYEAR;

	int d;

	if(date.year < nyear_spinup)
		year = first_historic_year;
	else
		year = date.year - nyear_spinup + first_historic_year;
*/
	int d = (year-FIRSTHISTYEAR)*365; 
	for(int i=0; i<365; i++)	{
		data[i] = daily_climate.Get(d,variable);
		d++;
	}
}


void SiteInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files), and reads in the gridlist

	//
	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon,dlat;
	bool eof=false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);


	gridlist.killall();
	first_call = true;

	while (!eof) {

		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		// Deal with rounding errors associated with reading in text files
		dlon = roundoff(dlon, 6);
		dlat = roundoff(dlat, 6);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
		}
	}

	file_daily_climate=param["file_daily_input_climate"].str;
	if(!daily_climate.Open(file_daily_climate)) {
		fail("initio: could not open %s for input",(char*)file_daily_climate);
	}

	fclose(in_grid);

	// Read CO2 data from file
	// co2.load_file(param["file_co2"].str);

	// Open landcover files
	landcover_input.init(param["file_gridlist"].str);
	// Open management files
	management_input.init();

	date.set_first_calendar_year(FIRSTHISTYEAR - nyear_spinup);

	// I've made the soilmap optional for eucface runs. We're not actually using it anyway.
	xtring soil_data_file = param["file_soildata"].str;
	soilinput.init(soil_data_file, translate_gridlist_to_coord(gridlist));

	// DH: David W added a soil depth input file which is used to
	// set soil depth which is in turn used to calculate some
	// further parameters which I don't think are related to eucface.
	// soilinput.init_depth(param["file_soildepth"].str, translate_gridlist_to_coord(gridlist));

	// DH: David W added an sdor (standard deviation of orography)
	// input file which allows sdor to be set on a per-gridcell basis.
	// The default value for sdor is 300. This is used to calculate b_oro
	// (Balsamo et al 2009 eqn 10) which in turn is not used at all. I'm
	// somewhat confident that this is related to P, not eucface, so
	// I've disabled this in our branch.
	// soilinput.init_sdor(param["file_soilsdor"].str, translate_gridlist_to_coord(gridlist));

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}


void SiteInput::get_monthly_ndep(int calendar_year,
                                double* mNHxdrydep, double* mNOydrydep,
								double* mNHxwetdep, double* mNOywetdep) {

	ndep.get_one_calendar_year(calendar_year,
	                           mNHxdrydep, mNOydrydep,
							   mNHxwetdep, mNOywetdep);
}


void SiteInput::adjust_raw_forcing_data(double lon,
                                       double lat,
                                       double hist_mtemp[NYEAR_HIST][12],
                                       double hist_mprec[NYEAR_HIST][12],
                                       double hist_msun[NYEAR_HIST][12]) {

	// The default (base class) implementation does nothing here.
}


bool SiteInput::getgridcell(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	int soilcode;
	int elevation;

	int days[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	
	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else
		gridlist.nextobj();

	if (gridlist.isobj) {

		bool gridfound = false;
		bool LUerror = false;
		double lon;
		double lat;

		while (!gridfound) {

			if (gridlist.isobj) {

				lon = gridlist.getobj().lon;
				lat = gridlist.getobj().lat;

				InData::Coord c;
				//Coord& c = gridlist.createobj();
				c.lon = lon;
				c.lat = lat;
				if (!daily_climate.Load(c)) {
					dprintf("Problems with daily climate file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n", c.lon, c.lat);
				}
				else {
					gridfound = true;
				}

				int y = 0;
				for (int year = FIRSTHISTYEAR; year < FIRSTHISTYEAR + NYEAR_HIST; year++) {

					//Coord& cc=gridlist.getobj();
					//InData::Coord c = GetLonLat(cc);

					get_daily_climate(gridcell, "Tair", dtemp, year);
					get_daily_climate(gridcell, "ddtr", ddtr, year);
					get_daily_climate(gridcell, "PAR", dpar, year);
					get_daily_climate(gridcell, "SWdown", dswdow, year);
					get_daily_climate(gridcell, "LWdown", dlwdow, year);
					get_daily_climate(gridcell, "Rain", dprec, year);
					get_daily_climate(gridcell, "RH", drhum, year);
					get_daily_climate(gridcell, "Wind", dwind, year);
					get_daily_climate(gridcell, "CO2ambient", daco2, year);
					get_daily_climate(gridcell, "CO2elevated", deco2, year);
					get_daily_climate(gridcell, "Ndep", dndep, year);

					// Assuming 50% Ndep for NH4 and NO3 each.
					for (int day = 0; day < 365; day++) {

						hist_dtemp[y][day] = dtemp[day];
						hist_ddtr[y][day] = ddtr[day];
						hist_dprec[y][day] = dprec[day];
						hist_dwind[y][day] = dwind[day];
						hist_drhum[y][day] = drhum[day];
						hist_dpar[y][day] = dpar[day];
						hist_dswdown[y][day] = dswdow[day];
						hist_dlwdown[y][day] = dlwdow[day];
						hist_daco2[y][day] = daco2[day];
						hist_deco2[y][day] = deco2[day];
						hist_dndep[y][day] = dndep[day];

					}
					y++;
				}

				if (run_landcover && gridfound) {
					LUerror = landcover_input.loadlandcover(gridlist.getobj().lon, gridlist.getobj().lat);
					if (!LUerror)
						LUerror = management_input.loadmanagement(gridlist.getobj().lon, gridlist.getobj().lat);
				}

				if (!gridfound || LUerror) {
					if (!gridfound)
						dprintf("\nError: could not find stand at (%g,%g) in climate data files\n\n", gridlist.getobj().lon, gridlist.getobj().lat);
					else if (LUerror)
						dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n\n", gridlist.getobj().lon, gridlist.getobj().lat);
					gridfound = false;
					gridlist.nextobj();
				}
			}
			else {
				return false;
			}
		}

		// Give sub-classes a chance to modify the data
		//adjust_raw_forcing_data(gridlist.getobj().lon,
		//                        gridlist.getobj().lat,
		//                        hist_mtemp, hist_mprec, hist_msun);

		// Build spinup data sets
		//spinup_mtemp.get_data_from(hist_mtemp);
		//spinup_mprec.get_data_from(hist_mprec);
		//spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		//spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		/*
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);
		spinup_mwind.get_data_from(hist_mwind);
		spinup_mrhum.get_data_from(hist_mrhum);
		*/
		// We wont detrend dtr for now. Partly because dtr is at the moment only
		// used for BVOC, so what happens during the spinup is not affecting
		// results in the period thereafter, and partly because the detrending
		// can give negative dtr values.
		//spinup_mdtr.detrend_data();


		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");

		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Get nitrogen deposition data.
		/* Since the historic data set does not reach decade 2010-2019,
		 * we need to use the RCP data for the last decade. */
		//ndep.getndep(param["file_ndep"].str, lon, lat, Lamarque::RCP60);

		// The insolation data will be sent (in function getclimate, below)
		// as incoming shortwave radiation, averages are over 24 hours
		/*if (param["insoltype"].num == 0) {
			gridcell.climate.instype = SWRAD_TS;
		}
		else {
			gridcell.climate.instype = NETSWRAD_TS;

		}*/

		gridcell.climate.instype = SWRAD;
			
		soilinput.get_soil(lon, lat, gridcell);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}


void SiteInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}


bool SiteInput::getclimate(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	double progress;

	Climate& climate = gridcell.climate;

	if (date.year < nyear_spinup + NYEAR_HIST) {

		if (date.day == 0) {

			int sim_year = date.get_calendar_year() - FIRSTHISTYEAR;

			// During spinup period	
			if (date.get_calendar_year() < FIRSTHISTYEAR) {

				sim_year = date.year % NYEAR_SPINUP_DATA;
			}

			for (int day = 0; day < date.MAX_YEAR_LENGTH; day++) {
				dtemp[day] = hist_dtemp[sim_year][day];
				ddtr[day] = hist_ddtr[sim_year][day];
				dprec[day] = hist_dprec[sim_year][day];
				dpar[day] = hist_dpar[sim_year][day];
				dwind[day] = hist_dwind[sim_year][day];
				drhum[day] = hist_drhum[sim_year][day];
				dpar[day] = hist_dpar[sim_year][day];
				dlwdow[day] = hist_dlwdown[sim_year][day];
				daco2[day] = hist_daco2[sim_year][day];
				deco2[day] = hist_deco2[sim_year][day];
				dndep[day] = hist_dndep[sim_year][day];
			}
		}
	}
	else {
		return false;
	}

	// Send environmental values for today to framework

	// CO2
	if (ambient_CO2) {
		climate.co2 = daco2[date.day];
	}
	else {
		climate.co2 = deco2[date.day];
	}

	climate.temp = dtemp[date.day] - K2degC;
	climate.prec = dprec[date.day];
	climate.insol = dswdow[date.day];
	climate.par_umol = dpar[date.day];

	// Nitrogen deposition g N m-2 yr-1 to 
	gridcell.dNH4dep = dndep[date.day] * KG_PER_G / date.MAX_YEAR_LENGTH / 2.0;
	gridcell.dNO3dep = dndep[date.day] * KG_PER_G / date.MAX_YEAR_LENGTH / 2.0;

	/* dh - commenting out p-specific changes

	// Phosphorus deposition // DW_COMMENT EucFACE value (0.093 kg P ha-1 year-1) 
	gridcell.dPO4dep = pdep * HA_PER_M2 / date.year_length();

	// Phosphorus fertilisation // DW_COMMENT EucFACE value (0.093 kg P ha-1 year-1) 
	gridcell.dPO4fert = 0.0;
	if (date.get_calendar_year() >= 2020 && date.get_calendar_year() <= 2022)
		gridcell.dPO4fert = euc_pfert / date.year_length(); 
	*/
	// Tmin, Tmax for BLAZE
	// initialise first
	climate.tmin = 0.;
	climate.tmax = 0.;
	if (firemodel == BLAZE) {
		climate.tmin = climate.temp - 0.5 * ddtr[date.day];
		climate.tmax = climate.temp + 0.5 * ddtr[date.day];
	}

	// Assuming rhum and wind are wanted when GWGEN is run
	// initialise first
	climate.u10 = 0.;
	climate.relhum = 0.;
	if (weathergenerator == GWGEN) {
		climate.u10 = dwind[date.day];
		climate.relhum = drhum[date.day];
	}
	climate.u10 = dwind[date.day];
	climate.relhum = drhum[date.day];

	// bvoc
	if (ifbvoc) {
		climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress() >= 1.0) {
			progress = (double)(gridlist.getobj().id * (nyear_spinup + NYEAR_HIST)
				+ date.year) / (double)(gridlist.nobj * (nyear_spinup + NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n", (int)(progress * 100.0),
				tprogress.elapsed.str, tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;

}


SiteInput::~SiteInput() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

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
