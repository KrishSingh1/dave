///////////////////////////////////////////////////////////////////////////////////////
/// \file cfinput.cpp
/// \brief Input module for CF conforming NetCDF files
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "cfinput.h"
#include "guess.h"
#include "driver.h"
#include <fstream>
#include <sstream>

REGISTER_INPUT_MODULE("cf", CFInput)

using namespace GuessNC::CF;

namespace {

insoltype cf_standard_name_to_insoltype(const std::string& standard_name) {
	if (standard_name == "surface_downwelling_shortwave_flux_in_air" ||
	    standard_name == "surface_downwelling_shortwave_flux") {
		return SWRAD_TS;
	}
	else if (standard_name == "surface_net_downward_shortwave_flux") {
		return NETSWRAD_TS;
	}
	else {
		fail("Unknown insolation type: %s", standard_name.c_str());
	}
}

// Checks if a DateTime is at the first day of the year
bool first_day_of_year(GuessNC::CF::DateTime dt) {
	return dt.get_month() == 1 && dt.get_day() == 1;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on an earlier day
bool earlier_day(const Date& date, int calendar_year, 
                 const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = calendar_year;
	d2[0] = date_time.get_year();
	
	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 < d2;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on a later day
bool later_day(const Date& date, int calendar_year,
               const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = calendar_year;
	d2[0] = date_time.get_year();
	
	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 > d2;	
}

}

CFInput::CFInput()
	: cf_temp(0),
	  cf_prec(0),
	  cf_insol(0),
	  lc_fixed_frac(NLANDCOVERTYPES, 0),
	  equal_landcover_area(false) {

	// Not used by this input module currently, but included as parameters so
	// common ins files can be used.

	declare_parameter("equal_landcover_area", &equal_landcover_area, "Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
	declare_parameter("lc_fixed_urban", &lc_fixed_frac[URBAN], 0, 100, "% lc_fixed_urban");
	declare_parameter("lc_fixed_cropland", &lc_fixed_frac[CROPLAND], 0, 100, "% lc_fixed_cropland");
	declare_parameter("lc_fixed_pasture", &lc_fixed_frac[PASTURE], 0, 100, "% lc_fixed_pasture");
	declare_parameter("lc_fixed_forest", &lc_fixed_frac[FOREST], 0, 100, "% lc_fixed_forest");
	declare_parameter("lc_fixed_natural", &lc_fixed_frac[NATURAL], 0, 100, "% lc_fixed_natural");
	declare_parameter("lc_fixed_peatland", &lc_fixed_frac[PEATLAND], 0, 100, "% lc_fixed_peatland");

}

CFInput::~CFInput() {
	delete cf_temp;
	delete cf_prec;
	delete cf_insol;
}

void CFInput::init() {

	// A warning about this input module not being proper from a scientific
	// perspective yet. For instance we're using historical ndep values for 
	// the future (if the NetCDF data set has a timespan that reaches further 
	// than the CRU data set).
	dprintf("Please note: this input module is a draft and not meant to be used for\n");
	dprintf("anything except technical evaluation of the file format.\n");

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	file_cru = param["file_cru"].str;
	
	// Try to open the NetCDF files
	try {
		cf_temp = new GridcellOrderedVariable(param["file_temp"].str, param["variable_temp"].str);
		cf_prec = new GridcellOrderedVariable(param["file_prec"].str, param["variable_prec"].str);
		cf_insol = new GridcellOrderedVariable(param["file_insol"].str, param["variable_insol"].str);
	}
	catch (const std::runtime_error& e) {
		fail(e.what());
	}

	// Make sure they contain what we expect

	if (cf_temp->get_standard_name() != "air_temperature") {
		fail("Temperature variable doesn't seem to contain air temperature data");
	}
	if (cf_temp->get_units() != "K") {
		fail("Temperature variable doesn't seem to be in Kelvin");
	}

	if (cf_prec->get_standard_name() != "precipitation_flux") {
		fail("Precipitation variable doesn't seem to contain precipitation flux data");
	}
	if (cf_prec->get_units() != "kg m-2 s-1") {
		fail("Precipitation variable doesn't seem to be in kg m-2 s-1");
	}

	if (cf_insol->get_standard_name() != "surface_downwelling_shortwave_flux_in_air" &&
	    cf_insol->get_standard_name() != "surface_downwelling_shortwave_flux" &&
	    cf_insol->get_standard_name() != "surface_net_downward_shortwave_flux") {
		fail("Insolation variable doesn't seem to contain insolation data");
	}
	if (cf_insol->get_units() != "W m-2") {
		fail("Insolation variable doesn't seem to be in W m-2");
	}

	// TODO: check that all variables have the same timespan
	// check time resolution?
	// other checks?


	// Read list of localities and store in gridlist member variable

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	std::ifstream ifs(file_gridlist, std::ifstream::in);

	if (!ifs.good()) fail("CFInput::init: could not open %s for input",(char*)file_gridlist);

	std::string line;
	while (getline(ifs, line)) {

		// Read next record in file
		int rlat, rlon;
		int landid;
		std::string descrip;
		Coord c;

		std::istringstream iss(line);

		if (cf_temp->is_reduced()) {
			if (iss >> landid) {
				iss >> descrip;

				c.landid = landid;
				c.descrip = descrip;

				gridlist.push_back(c);
			}
		}
		else {
			if (iss >> rlat >> rlon) {
				iss >> descrip;
				
				c.rlat = rlat;
				c.rlon = rlon;
				c.descrip = descrip;

				gridlist.push_back(c);
			}
		}
	}

	current_gridcell = gridlist.begin();

	date.set_first_calendar_year(cf_temp->get_date_time(0).get_year() - nyear_spinup);
}

bool CFInput::getgridcell(Gridcell& gridcell) {
	
	if (current_gridcell == gridlist.end()) {
		return false;
	}

	int rlon = current_gridcell->rlon;
	int rlat = current_gridcell->rlat;
	int landid = current_gridcell->landid;

	if (cf_temp->is_reduced()) {
		if (!cf_temp->load_data_for(landid) ||
		    !cf_prec->load_data_for(landid) ||
		    !cf_insol->load_data_for(landid)) {
			fail("Failed to load data for (%d) from NetCDF files", landid);
		}
	}
	else {
		if (!cf_temp->load_data_for(rlon, rlat) ||
		    !cf_prec->load_data_for(rlon, rlat) ||
		    !cf_insol->load_data_for(rlon, rlat)) {
			fail("Failed to load data for (%d, %d) from NetCDF files", rlat, rlon);
		}		
	}

	load_spinup_data(cf_temp, spinup_temp);
	load_spinup_data(cf_prec, spinup_prec);
	load_spinup_data(cf_insol, spinup_insol);

	spinup_temp.detrend_data();

	gridcell.climate.instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());

	double lon, lat;

	if (cf_temp->is_reduced()) {
		cf_temp->get_coords_for(landid, lon, lat);
	}
	else {
		cf_temp->get_coords_for(rlon, rlat, lon, lat);
	}

	gridcell.set_coordinates(lon, lat);

	// Find nearest CRU grid cell in order to get the soilcode

	int soilcode;
	double cru_lon = lon, cru_lat = lat;
	double dummy[CRU::NYEAR_HIST][12];

	const double searchradius = 1;

	if (!CRU::findnearestCRUdata(searchradius, file_cru, cru_lon, cru_lat, soilcode,
	                             dummy, dummy, dummy)) {
		fail("Failed to find soil code from CRU archive, close to coordinates (%g,%g)", cru_lon, cru_lat);
	}

	// Get nitrogen deposition, using the found CRU coordinates
	Lamarque::getndep(param["file_ndep"].str, cru_lon, cru_lat,
	                  NHxDryDep, NHxWetDep,
	                  NOyDryDep, NOyWetDep);

	// Setup the soil type
	soilparameters(gridcell.soiltype, soilcode);

	historic_timestep = -1;

	dprintf("\nCommencing simulation for stand at (%g,%g)", lon, lat);
	if (current_gridcell->descrip != "") 
		dprintf(" (%s)\n\n", current_gridcell->descrip.c_str());
	else dprintf("\n\n");
	
	return true;
}

void CFInput::populate_daily_arrays() {
	// Extract daily values for all days in this year, either from
	// spinup dataset or historical dataset

	int calendar_year = date.get_calendar_year();

	Date current_day = date;

	while (current_day.year == date.year) {

		// In the spinup?
		if (earlier_day(current_day, calendar_year, cf_temp->get_date_time(0))) {
			dtemp[current_day.day]  = spinup_temp[current_day.day];
			dprec[current_day.day]  = spinup_prec[current_day.day];
			dinsol[current_day.day] = spinup_insol[current_day.day];
		}
		else {
			// Historical period

			if (historic_timestep + 1 < cf_temp->get_timesteps()) {

				++historic_timestep;
				GuessNC::CF::DateTime dt = cf_temp->get_date_time(historic_timestep);

				if (dt.get_month() == 2 && dt.get_day() == 29) {
					++historic_timestep;
				}
			}
				
			if (historic_timestep < cf_temp->get_timesteps()) {
				dtemp[current_day.day]  = cf_temp->get_value(historic_timestep);
				dprec[current_day.day]  = cf_prec->get_value(historic_timestep);
				dinsol[current_day.day] = cf_insol->get_value(historic_timestep);
			}
			else {
				// Past the end of the historical period, these days wont be simulated.
				dtemp[current_day.day] = 0;
				dprec[current_day.day] = 0;
				dinsol[current_day.day] = 0;
			}
		}

		// Convert to units the model expects
		dtemp[current_day.day] -= K2degC;
		dprec[current_day.day] *= 3600*24;

		current_day.next();
	}

	// Move to next year in spinup dataset

	spinup_temp.nextyear();
	spinup_prec.nextyear();
	spinup_insol.nextyear();


	// Get monthly ndep values and convert to daily

	double mndrydep[12];
	double mnwetdep[12];

	// The ndep data set only goes up to 2009, after that we use the 2009 data
	Lamarque::get_one_calendar_year(min(2009, calendar_year),
	                                NHxDryDep, NHxWetDep,
	                                NOyDryDep, NOyWetDep,
	                                mndrydep, mnwetdep);

	// Distribute N deposition
	distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
}

void CFInput::getlandcover(Gridcell& gridcell) {
	for (int i = 0; i < NLANDCOVERTYPES; ++i) {
		gridcell.landcoverfrac[i] = 0;
	}
	gridcell.landcoverfrac[NATURAL] = 1;
}

bool CFInput::getclimate(Gridcell& gridcell) {
	
	Climate& climate = gridcell.climate;

	int calendar_year = date.get_calendar_year();

	GuessNC::CF::DateTime last_date = cf_temp->get_date_time(cf_temp->get_timesteps()-1);

	if (later_day(date, calendar_year, last_date)) {
		++current_gridcell;
		return false;
	}

	climate.co2 = co2[calendar_year];

	if (date.day == 0) {
		populate_daily_arrays();
	}

	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dinsol[date.day];

	// Nitrogen deposition
	climate.dndep = dndep[date.day];
		
	// Nitrogen fertilization
	climate.dnfert = 0.0;

	// bvoc
	if(ifbvoc){
		//	  climate.dtr=ddtr[date.day];
		fail("bvoc not supported by this input module");
	}

	return true;

}

void CFInput::load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
                               GenericSpinupData& spinup_data) {

	GenericSpinupData::RawData source;

	int timestep = 0;

	// Skip the first year if it doesn't start on Jan 1
	while (!first_day_of_year(cf_var->get_date_time(timestep))) {
		++timestep;
	}

	// Get all the daily values for the first NYEAR_SPINUP_DATA years, 
	// and put them into source
	for (int i = 0; i < NYEAR_SPINUP_DATA; ++i) {
		std::vector<double> year(GenericSpinupData::DAYS_PER_YEAR);

		for (int d = 0; d < year.size(); ++d) {
			GuessNC::CF::DateTime dt = cf_var->get_date_time(timestep);

			if (dt.get_month() == 2 && dt.get_day() == 29) {
				++timestep;
			}

			year[d] = cf_var->get_value(timestep);
			++timestep;
		}
		
		source.push_back(year);
	}

	spinup_data.get_data_from(source);
}
