///////////////////////////////////////////////////////////////////////////////////////
/// \file cruinput.h
/// \brief Input module for the CRU TS 3.0 data set
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CRUINPUT_H
#define LPJ_GUESS_CRUINPUT_H

#include "inputmodule.h"
#include <vector>
#include "gutil.h"
#include "globalco2file.h"
#include "spinupdata.h"
#include "cru_ts30.h"
#include "lamarquendep.h"

class CRUInput : public InputModule {
public:

	CRUInput();

	~CRUInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);


	// Constants associated with historical climate data set

	/// number of years of historical climate
	static const int NYEAR_HIST = CRU_TS30::NYEAR_HIST;

	/// calendar year corresponding to first year in data set
	static const int FIRSTHISTYEAR = CRU_TS30::FIRSTHISTYEAR;

	// number of years to use for temperature-detrended spinup data set
	// (not to be confused with the number of years to spinup model for, which
	// is read from the ins file)	
	static const int NYEAR_SPINUP_DATA=30;

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
	                      double* mndrydep,
	                      double* mnwetdep);

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

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	bool loadlandcover(Gridcell& gridcell, Coord c);

	/// search radius to use when finding CRU data
	double searchradius;

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether gridcell is divided into equal active landcover fractions.
	bool equal_landcover_area;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// The number of grid cells to simulate
	int ngridcell;

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	// Monthly temperature, precipitation and sunshine data for current grid cell
	// and historical period
	double hist_mtemp[NYEAR_HIST][12];
	double hist_mprec[NYEAR_HIST][12];
	double hist_msun[NYEAR_HIST][12];

	// Monthly frost days, precipitation days and DTR data for current grid cell
	// and historical period
	double hist_mfrs[NYEAR_HIST][12];
	double hist_mwet[NYEAR_HIST][12];
	double hist_mdtr[NYEAR_HIST][12];

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	// Spinup data sets for current grid cell
	Spinup_data spinup_mtemp;
	Spinup_data spinup_mprec;
	Spinup_data spinup_msun;

	// guess2008
	// Spinup data sets for monthly frost days, precipitation days and DTR data for 
	// current grid cell
	Spinup_data spinup_mfrs;
	Spinup_data spinup_mwet;
	Spinup_data spinup_mdtr;


	// Daily temperature, precipitation and sunshine for one year
	double dtemp[365],dprec[365],dsun[365];
	// bvoc
	// Daily diurnal temperature range for one year
	double ddtr[365];

	// Daily N deposition for one year
	double dndep[365];

	//Landuse:

	//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
	//TimeDataD input code may be put here
	TimeDataD LUdata(LOCAL_YEARLY);
	TimeDataD Peatdata;
#endif
	xtring file_lu, file_peat;
	static const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise
};

#endif // LPJ_GUESS_CRUINPUT_H
