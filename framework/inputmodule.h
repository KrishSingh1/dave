///////////////////////////////////////////////////////////////////////////////////////
/// \file inputmodule.h
/// \brief Base class for input modules
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_INPUT_MODULE_H
#define LPJ_GUESS_INPUT_MODULE_H

class Gridcell;

class InputModule {
public:
	/// Called after the instruction file has been read
	/** Initialises the input module (e.g. opening files). Typically
	 *  reads in a gridlist.
	 */
	virtual void init() = 0;

	/// Obtains coordinates and soil static parameters for the next grid cell to simulate 
	/** The function should return false if no grid cells remain to be simulated,
	 *  otherwise true. Currently the following member variables of gridcell should be
	 *  initialised: longitude, latitude and climate.instype; the following members of
	 *  member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	 *  thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	 *  code (Sitch et al 2000) by a call to function soilparameters in the driver
	 *  module (driver.cpp):
	 *
	 *  soilparameters(gridcell.soiltype,soilcode);
	 *
	 *  If the model is to be driven by quasi-daily values of the climate variables
	 *  derived from monthly means, this function may be the appropriate place to
	 *  perform the required interpolations. The utility functions interp_monthly_means
	 *  and interp_monthly_totals in driver.cpp may be called for this purpose.
	 */
	virtual bool getgridcell(Gridcell& gridcell) = 0;

	/// Obtains climate data (including atmospheric CO2 and insolation) for this day
	/** The function should return false if the simulation is complete for this grid cell,
	 *  otherwise true. This will normally require querying the year and day member
	 *  variables of the global class object date:
	 *
	 *  if (date.day==0 && date.year==nyear_spinup) return false;
	 *  // else
	 *  return true;
	 *
	 *  Currently the following member variables of the climate member of gridcell must be
	 *  initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	 *  values of the climate variables derived from monthly means, this day's values
	 *  will presumably be extracted from arrays containing the interpolated daily
	 *  values (see function getgridcell):
	 *
	 *  gridcell.climate.temp=dtemp[date.day];
	 *  gridcell.climate.prec=dprec[date.day];
	 *  gridcell.climate.insol=dsun[date.day];
	 *
	 *  Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
	 *  BVOC:
	 *  gridcell.climate.dtr=ddtr[date.day]; 
	 *
	 *  If model is run in diurnal mode, which requires appropriate climate forcing data, 
	 *  additional members of the climate must be initialised: temps, insols. Both of the
	 *  variables must be of type std::vector. The length of these vectors should be equal
	 *  to value of date.subdaily which also needs to be set either in getclimate or 
	 *  getgridcell functions. date.subdaily is a number of sub-daily period in a single 
	 *  day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
	 *  diurnal mode.
	 */
	virtual bool getclimate(Gridcell& gridcell) = 0;

	/// Sets land cover fractions for the gridcell for the current year
	virtual void getlandcover(Gridcell& gridcell) = 0;
};

#endif // LPJ_GUESS_INPUT_MODULE_H
