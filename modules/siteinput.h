///////////////////////////////////////////////////////////////////////////////////////
/// \file cruinput.h
/// \brief Input module for the CRU-NCEP data set
///
/// \author Joe Siltberg
/// $Date: 2019-11-06 17:45:44 +0100 (ons, 06 nov 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_SITEINPUT_H
#define LPJ_GUESS_SITEINPUT_H

#include "guess.h"
#include "inputmodule.h"
#include "siteinput.h"
#include <vector>
#include "gutil.h"
#include "globalco2file.h"
#include "soilinput.h"
#include "spinupdata.h"
#include "lamarquendep.h"
#include "externalinput.h"
#include "indata.h"


/// An input module for CRU climate data
/** This input module gets climate data from binary archives built from
 *  CRU-NCEP (1901-2015).
 */
class SiteInput : public InputModule {
public:

	/// Constructor
	/** Declares the instruction file parameters used by the input module.
	 */
	SiteInput();

	/// Destructor, cleans up used resources
	~SiteInput();

	/// Reads in gridlist and initialises the input module
	/** Gets called after the instruction file has been read */
	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	void getlandcover(Gridcell& gridcell);

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) { management_input.getmanagement(gridcell); }

	// Constants associated with historical climate data set

	/// number of years of historical climate
	static const int NYEAR_HIST = 370;// 119; // EucFACE data

	/// calendar year corresponding to first year in data set
	static const int FIRSTHISTYEAR = 1700; //  1901; // EucFACE data

	/// number of years to use for temperature-detrended spinup data set
	/** (not to be confused with the number of years to spinup model for, which
	 * is read from the ins file)
	 */
	static const int NYEAR_SPINUP_DATA = 50; // 30; // EucFACE data

	void get_daily_climate(Gridcell& gridcell, xtring variable, double data[], int year);

	//InData::Coord SiteInput::GetLonLat(SiteInput::Coord coord);

protected:

	/// Gets monthly ndep values for a given calendar year
	/** To be used by sub-classes that wish to do their own
	 *  distribution of monthly values to daily values.
	 *
	 *  The ndep values returned are for the current gridcell,
	 *  i.e. the gridcell chosen in the most recent call to
	 *  getgridcell().
	 *
	 *  \param calendar_year The calendar (not simulation!) year for which to get ndep
	 *  \param mndrydep      Pointer to array holding 12 doubles
	 *  \param mnwetdep      Pointer to array holding 12 doubles
	 */
	void get_monthly_ndep(int calendar_year,
		double* mNHxdrydep, double* mNOydrydep,
		double* mNHxwetdep, double* mNOywetdep);

	/// Gives sub-classes a chance to modify the forcing data
	/** This function will be called just after the forcing data for the historical
	 *  period has been read in for a gridcell. Sub-classes can override this function
	 *  and modify the data if needed, for instance adjusting according to site data.
	 *
	 *  Note that modifying this data will also affect the spinup period since
	 *  the spinup forcing is based on the historical period.
	 *
	 *  \param hist_mtemp  Monthly temperature values for each year
	 *  \param hist_mprec  Monthly precipitation values for each year
	 *  \param hist_msun   Monthly sunshine values for each year
	 */
	virtual void adjust_raw_forcing_data(double lon,
		double lat,
		double hist_mtemp[NYEAR_HIST][12],
		double hist_mprec[NYEAR_HIST][12],
		double hist_msun[NYEAR_HIST][12]);

private:

	/// Type for storing grid cell longitude, latitude and description text
	struct Coord {

		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	std::vector<std::pair<double, double> > translate_gridlist_to_coord(ListArray_id<Coord>& gridlist);

	SoilInput soilinput;

	/// Land cover input module
	LandcoverInput landcover_input;
	/// Management input module
	ManagementInput management_input;

	/// search radius to use when finding CRU data
	double searchradius;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// Flag for getgridcell(). True indicates that the first gridcell has not been read yet by getgridcell()
	bool first_call;

	// Timers for keeping track of progress through the simulation
	Timer tprogress, tmute;
	static const int MUTESEC = 20; // minimum number of sec to wait between progress messages

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	InData::TimeDataD daily_climate;
	xtring file_daily_climate;

	/// Daily temperature for current grid cell and historical period
	double hist_dtemp[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily precipitation for current grid cell and historical period
	double hist_dprec[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily sunshine for current grid cell and historical period
	double hist_dpar[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily frost days for current grid cell and historical period
	double hist_dfrs[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily precipitation days for current grid cell and historical period
	double hist_dwet[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily DTR (diurnal temperature range) for current grid cell and historical period
	double hist_ddtr[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean wind for current grid cell and historical period
	double hist_dwind[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_drhum[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_dlwdown[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_dswdown[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_daco2[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_deco2[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Daily mean relative humidity for current grid cell and historical period
	double hist_dndep[NYEAR_HIST][Date::MAX_YEAR_LENGTH];

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Spinup data for current grid cell - temperature
	Spinup_data spinup_mtemp;
	/// Spinup data for current grid cell - precipitation
	Spinup_data spinup_mprec;
	/// Spinup data for current grid cell - sunshine
	Spinup_data spinup_msun;

	/// Spinup data for current grid cell - frost days
	Spinup_data spinup_mfrs;
	/// Spinup data for current grid cell - precipitation days
	Spinup_data spinup_mwet;
	/// Spinup data for current grid cell - DTR (diurnal temperature range)
	Spinup_data spinup_mdtr;
	/// Spinup data for current grid cell - wind
	Spinup_data spinup_mwind;
	/// Spinup data for current grid cell - relative humidity
	Spinup_data spinup_mrhum;

	// Daily temperature for current year
	double dtemp[Date::MAX_YEAR_LENGTH];
	// Daily precipitation for current year
	double dprec[Date::MAX_YEAR_LENGTH];
	// Daily sunshine for current year
	double dpar[Date::MAX_YEAR_LENGTH];
	// Daily shortwave radiation for current year
	double dswdow[Date::MAX_YEAR_LENGTH];
	// Daily longwave radiation for current year
	double dlwdow[Date::MAX_YEAR_LENGTH];
	// Daily diurnal temperature range for current year
	double ddtr[Date::MAX_YEAR_LENGTH];
	// Daily mean wind 
	double dwind[Date::MAX_YEAR_LENGTH];
	// Daily mean relative humidity
	double drhum[Date::MAX_YEAR_LENGTH];
	// Daily mean ambient CO2
	double daco2[Date::MAX_YEAR_LENGTH];
	// Daily mean elevated CO2
	double deco2[Date::MAX_YEAR_LENGTH];
	// Daily N deposition for current year
	double dndep[Date::MAX_YEAR_LENGTH];
	double dNH4dep[Date::MAX_YEAR_LENGTH];
	double dNO3dep[Date::MAX_YEAR_LENGTH];
};

#endif // LPJ_GUESS_CRUINPUT_H
