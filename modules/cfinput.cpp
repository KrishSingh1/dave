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

}

CFInput::CFInput() {
}

void CFInput::init() {

	CRUInput::init();

	// A warning about this input module not being proper from a scientific
	// perspective yet. For instance we're not doing the spinup properly
	// yet, and we're using historical ndep values for the future (if
	// the NetCDF data set has a timespan that reaches further than the
	// CRU data set). Also ndep isn't distributed correctly according
	// to wet days.
	dprintf("Please note: this input module is a draft and not meant to be used for");
	dprintf("anything except technical evaluation of the file format.");

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);
	
	// Try to open the NetCDF files
	try {
		cf_temp.reset(new GridcellOrderedVariable(param["file_temp"].str, param["variable_temp"].str));
		cf_prec.reset(new GridcellOrderedVariable(param["file_prec"].str, param["variable_prec"].str));
		cf_insol.reset(new GridcellOrderedVariable(param["file_insol"].str, param["variable_insol"].str));
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

	gridcell.climate.instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());

	historic_timestep = -1;
	spinup_timestep = -1;

	return true;
}

bool CFInput::getclimate(Gridcell& gridcell) {
	
	// We won't call the base class' getclimate here since:
	// - the mapping from simulation year to calendar year might be different
	//   (so we could get incorrect values for i.e. co2 or ndep)
	// - we want to be able to continue after the last CRU year if needed

	Climate& climate = gridcell.climate;

	int calendar_year = cf_temp->get_date_time(0).get_year() + date.year - nyear_spinup;

	climate.co2 = co2[calendar_year];

	if (date.day == 0) {
		double dprec[365];
		double mndrydep[12];
		double mnwetdep[12];

		// The ndep data set only goes up to 2009, after that we use the 2009 data
		get_monthly_ndep(min(2009, calendar_year), mndrydep, mnwetdep);

		// Distribute N deposition - without rain days
		std::fill_n(dprec, 365, 0);
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
	}

	int timestep;

	if (date.year < nyear_spinup) {
		++spinup_timestep;

		GuessNC::CF::DateTime dt_start = cf_temp->get_date_time(0);
		GuessNC::CF::DateTime dt_now = cf_temp->get_date_time(spinup_timestep);
		
		if (dt_now.get_year() - dt_start.get_year() == NYEAR_SPINUP_DATA) {
			spinup_timestep = 0;
		}
		else if (dt_now.get_month() == 2 && dt_now.get_day() == 29) {
			++spinup_timestep;
		}

		timestep = spinup_timestep;
	}
	else {
		++historic_timestep;
		GuessNC::CF::DateTime dt = cf_temp->get_date_time(historic_timestep);

		if (dt.get_month() == 2 && dt.get_day() == 29) {
			++historic_timestep;
		}

		timestep = historic_timestep;
	}

	climate.temp  = cf_temp->get_value(timestep)-K2degC;
	climate.prec  = cf_prec->get_value(timestep)*3600*24;
	climate.insol = cf_insol->get_value(timestep);

	// Nitrogen deposition
	climate.dndep = dndep[date.day];
		
	// Nitrogen fertilization
	climate.dnfert = 0.0;

	// bvoc
	if(ifbvoc){
		//	  climate.dtr=ddtr[date.day];
		fail("bvoc not supported by this input module");
	}

	if (historic_timestep + 1 == cf_temp->get_timesteps()) {
		return false;
	}

	return true;

}
