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
#include "guess.h"
#if defined DYNAMIC_LANDCOVER_INPUT
#include "InData.h"
#endif

class CRUInput : public InputModule {
public:

	CRUInput();

	~CRUInput();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	/// Gets landcover and crop fractions for one year
	void getlandcover(Gridcell& gridcell);

	/// Gets sowing dates for one year
	void getsowingdates(Gridcell& gridcell);

	/// Gets harvest dates for one year
	void getharvestdates(Gridcell& gridcell);

	// Constants associated with historical climate data set

	// number of years of historical climate
	// CRU TS 3.0 has 106 years of data (1901-2006)
	static const int NYEAR_HIST=106;

	// calender year corresponding to first year in CRU climate data set
	static const int FIRSTHISTYEAR=1901;

	// calender year corresponding to first year nitrogen deposition
	static const int FIRSTHISTYEARNDEP=1850;

	// number of years of historical nitrogen deposition 
	static const int NYEAR_HISTNDEP=16;

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

	/// Loads landcover and crop fractions plus sowing and harvest dates from input files
	bool loadlandcover(Gridcell& gridcell, Coord c);

#if defined DYNAMIC_LANDCOVER_INPUT
	/// Transfers coordinates from CRUInput::Coord to InData::Coord
	InData::Coord GetLonLat(Coord coord);
	/// Transfers gridlist of coordinates from CRUInput::Coord to InData::Coord
	void GetLonLatList(ListArray_id<InData::Coord>& lonlatlist, ListArray_id<Coord>& gridlist);
#endif

	void getndep(double lon, double lat);

	/// search radius to use when finding CRU data
	double searchradius;

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether enforced static landcover fractions are equal-sized stands of all included landcovers
	bool equal_landcover_area;

	/// Whether pfts not in crop fraction input file are removed from pftlist (0,1)
	bool minimizecftlist;

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

	/// Monthly data on daily dry NHx deposition (kgN/m2/day)
	double NHxDryDep[NYEAR_HISTNDEP][12];
	/// Monthly data on daily wet NHx deposition (kgN/m2/day)
	double NHxWetDep[NYEAR_HISTNDEP][12];
	/// Monthly data on daily dry NOy deposition (kgN/m2/day)
	double NOyDryDep[NYEAR_HISTNDEP][12];
	/// Monthly data on daily wet NOy deposition (kgN/m2/day)
	double NOyWetDep[NYEAR_HISTNDEP][12];

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

	//Landuse input:

#if defined DYNAMIC_LANDCOVER_INPUT

	// Objects handling landcover fraction data input
	InData::TimeDataD LUdata;
	InData::TimeDataD Peatdata;
	InData::TimeDataD CFTdata;
	InData::TimeDataD sdates;
	InData::TimeDataD hdates;

#ifdef LUTOMEMORY
	InData::TimeDataDmem LUdata_mem;
	InData::TimeDataDmem CFTdata_mem;
#endif

#endif
	xtring file_lu, file_lucrop, file_peat, file_sdates, file_hdates;
	// Number of years of landcover fraction data in input files
	static const int NYEAR_LU=103;
};

#endif // LPJ_GUESS_CRUINPUT_H
