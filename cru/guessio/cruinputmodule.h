///////////////////////////////////////////////////////////////////////////////////////
/// \file cruinputmodule.h
/// \brief Input module for the CRU TS 3.0 data set
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CRUINPUTMODULE_H
#define LPJ_GUESS_CRUINPUTMODULE_H

#include "inputmodule.h"
#include <vector>
#include "gutil.h"
#include "globalco2file.h"
#include "spinupdata.h"

class CRUInputModule : public InputModule {
public:

	CRUInputModule();

	~CRUInputModule();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);


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


private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		
		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	bool loadlandcover(Gridcell& gridcell, Coord c);

	void getndep(double lon, double lat);

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

#endif // LPJ_GUESS_CRUINPUTMODULE_H
