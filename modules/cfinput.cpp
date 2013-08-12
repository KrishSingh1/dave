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
	  cf_insol(0) {
}

CFInput::~CFInput() {
	delete cf_temp;
	delete cf_prec;
	delete cf_insol;
}

void CFInput::init() {

	CRUInput::init();

	// A warning about this input module not being proper from a scientific
	// perspective yet. For instance we're not doing the spinup properly
	// yet, and we're using historical ndep values for the future (if
	// the NetCDF data set has a timespan that reaches further than the
	// CRU data set). Also ndep isn't distributed correctly according
	// to wet days.
	dprintf("Please note: this input module is a draft and not meant to be used for\n");
	dprintf("anything except technical evaluation of the file format.\n");

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);
	
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
}

bool CFInput::getgridcell(Gridcell& gridcell) {
	if (!CRUInput::getgridcell(gridcell)) {
		return false;
	}

	// Somehow get these based on the gridcell's coordinates
	int rlon = 0, rlat = 0;

	// If it didn't work, call CRUInput::getgridcell until it works

	cf_temp->load_data_for(rlon, rlat);
	cf_prec->load_data_for(rlon, rlat);
	cf_insol->load_data_for(rlon, rlat);

	load_spinup_data(cf_temp, spinup_temp);
	load_spinup_data(cf_prec, spinup_prec);
	load_spinup_data(cf_insol, spinup_insol);

	spinup_temp.detrend_data();

	gridcell.climate.instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());

	historic_timestep = -1;

	return true;
}

bool CFInput::getclimate(Gridcell& gridcell) {
	
	// We won't call the base class' getclimate here since:
	// - the mapping from simulation year to calendar year might be different
	//   (so we could get incorrect values for i.e. co2 or ndep)
	// - we want to be able to continue after the last CRU year if needed

	Climate& climate = gridcell.climate;

	int calendar_year = cf_temp->get_date_time(0).get_year() + date.year - nyear_spinup;

	GuessNC::CF::DateTime last_date = cf_temp->get_date_time(cf_temp->get_timesteps()-1);

	if (later_day(date, calendar_year, last_date)) {
		return false;
	}

	climate.co2 = co2[calendar_year];

	if (date.day == 0) {

		// Extract daily values for all days in this year, either from
		// spinup dataset or historical dataset

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
		get_monthly_ndep(min(2009, calendar_year), mndrydep, mnwetdep);

		// Distribute N deposition
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
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
