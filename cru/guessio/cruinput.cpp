///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio_cru.cpp
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
#include "globalco2file.h"
#if defined DYNAMIC_LANDCOVER_INPUT
#include "InData.h"
#endif
// guess2008 - header file for the CRU TS 3.0 data archives
#include "cru_1901_2006.h"
#include "cru_1901_2006misc.h"

// header file for reading binary data archive of global nitrogen deposition
#include "GlobalNitrogenDeposition.h"

bool fixedtemp_hist=0;
bool fixedprec_hist=0;
bool fixedrad_hist=0;
bool fixedco2_hist=0;
bool fixedlu_hist=0;
bool fixedcrop_hist=0;

REGISTER_INPUT_MODULE("cru", CRUInput)

CRUInput::CRUInput()
	: searchradius(0),
	  lc_fixed_frac(NLANDCOVERTYPES, 0),
	  equal_landcover_area(false),
	  spinup_mtemp(NYEAR_SPINUP_DATA),
	  spinup_mprec(NYEAR_SPINUP_DATA),
	  spinup_msun(NYEAR_SPINUP_DATA),
	  spinup_mfrs(NYEAR_SPINUP_DATA),
	  spinup_mwet(NYEAR_SPINUP_DATA),
	  spinup_mdtr(NYEAR_SPINUP_DATA) {

	// Declare instruction file parameters

	declare_parameter("searchradius", &searchradius, 0, 100,
		"If specified, CRU data will be searched for in a circle");

	declare_parameter("equal_landcover_area", &equal_landcover_area, "Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
	declare_parameter("minimizecftlist", &minimizecftlist, "Whether pfts not in crop fraction input file are removed from pftlist (0,1)");
	declare_parameter("lc_fixed_urban", &lc_fixed_frac[URBAN], 0, 100, "% lc_fixed_urban");
	declare_parameter("lc_fixed_cropland", &lc_fixed_frac[CROPLAND], 0, 100, "% lc_fixed_cropland");
	declare_parameter("lc_fixed_pasture", &lc_fixed_frac[PASTURE], 0, 100, "% lc_fixed_pasture");
	declare_parameter("lc_fixed_forest", &lc_fixed_frac[FOREST], 0, 100, "% lc_fixed_forest");
	declare_parameter("lc_fixed_natural", &lc_fixed_frac[NATURAL], 0, 100, "% lc_fixed_natural");
	declare_parameter("lc_fixed_peatland", &lc_fixed_frac[PEATLAND], 0, 100, "% lc_fixed_peatland");
}


// guess2008 - make file_cru and file_cru_misc global variables
xtring file_cru;
xtring file_cru_misc;

namespace {

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double mtemp[12], double mprec[12], double msun[12], double mdtr[12],
					double dtemp[365], double dprec[365], double dsun[365], double ddtr[365]) {
	interp_monthly_means(mtemp, dtemp);
	interp_monthly_totals(mprec, dprec);
	interp_monthly_means(msun, dsun);
	interp_monthly_means(mdtr, ddtr);
}

} // namespace

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 
// guess2008 - new functions for reading CRU TS 3.0 binary files.
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU
// Determine temp, precip, sunshine & soilcode
 
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
               double mtemp[CRUInput::NYEAR_HIST][12],
               double mprec[CRUInput::NYEAR_HIST][12],
               double msun[CRUInput::NYEAR_HIST][12]) {

	// !!!! NEW VERSION OF THIS FUNCTION - guess2008 - NEW VERSION OF THIS FUNCTION !!!!
	// Please note the new function signature. 

	// Archive object. Definition in new header file, cru.h
	Cru_1901_2006Archive ark;

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		soilcode=(int)data.soilcode[0];


		for (y = 0; y < CRUInput::NYEAR_HIST; y++) {
			for (m=0;m<12;m++) {
				if(fixedtemp_hist)
					mtemp[y][m] = data.mtemp[m]*0.1;
				else
					mtemp[y][m] = data.mtemp[y*12+m]*0.1; // now degC
				if(fixedprec_hist)
					mprec[y][m] = data.mprec[m]*0.1;
				else
					mprec[y][m] = data.mprec[y*12+m]*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;

				if(fixedrad_hist)
					msun[y][m]  = data.msun[m]*0.1;
				else
					msun[y][m]  = data.msun[y*12+m]*0.1;   // % sun 

			}
		}


		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}




///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU_MISC
// Determine elevation, frs frq, wet frq & DTR

bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
                    double mfrs[CRUInput::NYEAR_HIST][12],
                    double mwet[CRUInput::NYEAR_HIST][12],
                    double mdtr[CRUInput::NYEAR_HIST][12]) {
	
	// Please note the new function signature. 

	// Archive object
	Cru_1901_2006miscArchive ark; 
	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006misc data;

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays.
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000

		for (y = 0; y < CRUInput::NYEAR_HIST; y++) { 
			for (m=0;m<12;m++) {

				// guess2008 - catch rounding errors 
				mfrs[y][m] = data.mfrs[y*12+m]*0.01; // days
				if (mfrs[y][m] < 0.1) 
					mfrs[y][m] = 0.0; // Catches rounding errors

				mwet[y][m] = data.mwet[y*12+m]*0.01; // days
				if (mwet[y][m] <= 0.1) 
					mwet[y][m] = 0.0; // Catches rounding errors

				mdtr[y][m] = data.mdtr[y*12+m]*0.1;  // degC

				/*
				If vapour pressure is needed:
				mvap[y][m] = data.mvap[y*12+m]*0.01;
				*/
			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


// Utility function that returns the CRU data from the nearest cell to (lon,lat) within
// a given search radius
bool findnearestCRUdata(double searchradius, char* cruark, double& lon, double& lat, 
                        int& scode, 
                        double hist_mtemp1[CRUInput::NYEAR_HIST][12], 
                        double hist_mprec1[CRUInput::NYEAR_HIST][12], 
                        double hist_msun1[CRUInput::NYEAR_HIST][12]) {

	// First try the exact coordinate
	if (searchcru(cruark, lon, lat, scode, hist_mtemp1, hist_mprec1, hist_msun1)) {
		return true;
	}
	
	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon*2)/2;
	double center_lat = floor(lat*2)/2;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.5;
	const double EPS = 1e-15;

	for (double y = center_lon-searchradius; y <= center_lon+searchradius+EPS; y += STEP) {
		for (double x = center_lat-searchradius; x <= center_lat+searchradius+EPS; x += STEP) {
			double xdist = x-center_lat;
			double ydist = y-center_lon;
			double dist = sqrt(xdist*xdist + ydist*ydist);
			
			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found in CRU
	for (unsigned int i = 0; i < search_points.size(); i++) {
		point search_point = search_points[i].second;
		double search_lon = search_point.first;
		double search_lat = search_point.second;

		if (searchcru(cruark, search_lon, search_lat, scode, 
		              hist_mtemp1, hist_mprec1, hist_msun1)) {
			lon = search_lon;
			lat = search_lat;
			return true;
		}
	}

	return false;
}


void CRUInput::init() {

	// DESCRIPTION
	// Initialises input (e.g. opening files), and reads in the gridlist

	///////////////////////////////////////////////////////////////////////////////////
	// GENERIC SECTION - DO NOT MODIFY

	xtring header;

	unixtime(header);
	header=(xtring)"[LPJ-GUESS  "+header+"]\n\n";
	dprintf((char*)header);

	// Print the title of this run
	dprintf("\n\n-----------------------------------------------\n%s\n-----------------------------------------------\n",(char*)title);

	///////////////////////////////////////////////////////////////////////////////////
	// USER-SPECIFIC SECTION (Modify as necessary or supply own code)
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

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;

	
	ngridcell=0;
	while (!eof) {
		
		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
			ngridcell++;
		}
	}


	fclose(in_grid);

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	if (run_landcover) {
		all_fracs_const=true;	//If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year

		//Retrieve file names for landcover files and open them if static values from ins-file are not used
		if (!lcfrac_fixed) {

			bool openLUfile = false;

			for(int i=0; i<NLANDCOVERTYPES; i++) {
				if(run[i] && i != NATURAL)
					openLUfile = true;
			}

			if (openLUfile) {
				file_lu=param["file_lu"].str;
#if defined DYNAMIC_LANDCOVER_INPUT
				// Open landcover area fraction file, return false if problem
				if(!LUdata.Open(file_lu))
					fail("initio: could not open %s for input",(char*)file_lu);
				else {
					if(LUdata.format==InData::LOCAL_YEARLY)
						all_fracs_const=false;				//Set all_fracs_const to false if yearly data

#ifdef LUTOMEMORY
					// Save all landcover area fraction data in memory
					LUdata_mem.Open(gridlist.nobj, LUdata.nRecords,LUdata.nYears);

					ListArray_id<InData::Coord> lonlatlist;
					GetLonLatList(lonlatlist, gridlist);
					LUdata_mem.CopyFromTimeDataD(LUdata, lonlatlist);
#endif
				}
#endif
			}
		}

		//Retrieve file names for crop fraction file and open them if static equal-size values are not used.
		if(run[CROPLAND] && !cftfrac_fixed)
		{
			file_lucrop=param["file_lucrop"].str;
#if defined DYNAMIC_LANDCOVER_INPUT
			// Open crop fraction file, return false if problem
			if(!CFTdata.Open(file_lucrop))
				fail("initio: could not open %s for input",(char*)file_lucrop);
			else if(minimizecftlist) {
				// remove all crop pft:s from gridlist that always have zero area fraction
				ListArray_id<InData::Coord> lonlatlist;
				GetLonLatList(lonlatlist, gridlist);
				CFTdata.CheckIfPresent(lonlatlist);
				
				int n=0;
				pftlist.firstobj();
				while(pftlist.isobj) {		
					if(pftlist.getobj().cftid>=0 && !CFTdata.CFTPresent(pftlist.getobj().cftid) && 
						!(pftlist.getobj().isintercropgrass && ifintercropgrass)) {
						n+=1;
						pftlist.killobj();
						npft--;
						ncft--;
					}
					else {
						pftlist.getobj().id-=n;
						if(pftlist.getobj().cftid>=0)
							CFTdata.active[pftlist.getobj().cftid]=1;
						pftlist.nextobj();
					}			
				}
			}
			else {
				pftlist.firstobj();
				while(pftlist.isobj) {
					if(pftlist.getobj().cftid>=0)
						CFTdata.active[pftlist.getobj().cftid]=1;
					pftlist.nextobj();
				}
			}

			if(CFTdata.format==InData::LOCAL_YEARLY) 
				all_fracs_const=false;				// Set all_fracs_const to false if yearly data

#ifdef LUTOMEMORY
				// Save all crop area fraction data in memory
				CFTdata_mem.Open(gridlist.nobj, CFTdata.nRecords,CFTdata.nYears);

				ListArray_id<InData::Coord> lonlatlist;
				GetLonLatList(lonlatlist, gridlist);
				CFTdata_mem.CopyFromTimeDataD(CFTdata, lonlatlist);
#endif

//			for(int i=0;i<CFTdata.nRecords;i++)
//				dprintf("%s:CFTdata.active=%d\n", CFTdata.GetHeader(i), CFTdata.active[i]);

			if(CFTdata.GetnRecords()!=NCROPSTANDS_MAX)
				fail("\ninitio: NCROPSTANDS_MAX is incorrectly set in guess.h !\n");
#endif
		}

		if(run[CROPLAND]) {
#if defined DYNAMIC_LANDCOVER_INPUT
			if(forcesowingdates)
			{
				file_sdates=param["file_sdates"].str;
				if(!sdates.Open(file_sdates))
					fail("initio: could not open %s for input",(char*)file_sdates);
#ifdef LUTOMEMORY
				// Save all harvest date data in memory
				sdates_mem.Open(gridlist.nobj, sdates.nRecords,sdates.nYears);

				ListArray_id<InData::Coord> lonlatlist;
				GetLonLatList(lonlatlist, gridlist);
				sdates_mem.CopyFromTimeDataD(sdates, lonlatlist);
#endif
			}
			if(forceharvestdates)
			{
				file_hdates=param["file_hdates"].str;
				if(!hdates.Open(file_hdates))
					fail("initio: could not open %s for input",(char*)file_hdates);
#ifdef LUTOMEMORY
				// Save all harvest date data in memory
				hdates_mem.Open(gridlist.nobj, hdates.nRecords,hdates.nYears);

				ListArray_id<InData::Coord> lonlatlist;
				GetLonLatList(lonlatlist, gridlist);
				hdates_mem.CopyFromTimeDataD(hdates, lonlatlist);
#endif
			}
			if(readNfert) {
				file_Nfert=param["file_Nfert"].str;
				if(!Nfert.Open(file_Nfert))
					fail("initio: could not open %s for input",(char*)file_Nfert);
#ifdef LUTOMEMORY
				// Save all N fertilization data in memory
				Nfert_mem.Open(gridlist.nobj, Nfert.nRecords,Nfert.nYears);

				ListArray_id<InData::Coord> lonlatlist;
				GetLonLatList(lonlatlist, gridlist);
				Nfert_mem.CopyFromTimeDataD(Nfert, lonlatlist);
#endif
			}

#endif
		}
	}
	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}


void CRUInput::get_monthly_ndep(int calendar_year,
                                double* mndrydep,
                                double* mnwetdep) {
	int ndep_year = 0;

	if (calendar_year >= FIRSTHISTYEARNDEP) {
		ndep_year = (int)((calendar_year - FIRSTHISTYEARNDEP)/10);
	}

	for (int m = 0; m < 12; m++) {
		mndrydep[m] = NHxDryDep[ndep_year][m] + NOyDryDep[ndep_year][m];
		
		mnwetdep[m] = NHxWetDep[ndep_year][m] + NOyWetDep[ndep_year][m];
	}
}


void CRUInput::adjust_raw_forcing_data(double lon,
                                       double lat,
                                       double hist_mtemp[NYEAR_HIST][12],
                                       double hist_mprec[NYEAR_HIST][12],
                                       double hist_msun[NYEAR_HIST][12]) {

	// The default (base class) implementation does nothing here.
}


#if defined DYNAMIC_LANDCOVER_INPUT

/// Transfers coordinates from CRUInput::Coord to InData::Coord
InData::Coord CRUInput::GetLonLat(Coord coord) {

	InData::Coord lonlat;
	lonlat.lon=coord.lon;
	lonlat.lat=coord.lat;

	return lonlat;
}

/// Transfers gridlist of coordinates from CRUInput::Coord to InData::Coord
void CRUInput::GetLonLatList(ListArray_id<InData::Coord>&lonlatlist, ListArray_id<Coord>& gridlist) {

	for(unsigned int i = 0; i < gridlist.nobj; i++) {
		InData::Coord& c= lonlatlist.createobj();
		c.lon=gridlist[i].lon;
		c.lat=gridlist[i].lat;
	}
}

#endif

///	Loads landcover area fraction data from file(s) for a gridcell.
/** Called from getgridcell() if run_landcover is true. 
  */
bool CRUInput::loadlandcover(Gridcell& gridcell, Coord cc) {

#if defined DYNAMIC_LANDCOVER_INPUT
	InData::Coord c = GetLonLat(cc);
#endif
	bool LUerror=false;

	if (!lcfrac_fixed) {

		// Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
		// transferred to gridcell.landcoverfrac each year in getlandcover()

			bool loadLU = false;

			for(int i=0; i<NLANDCOVERTYPES; i++) {
				if(run[i] && i != NATURAL)
					loadLU = true;
			}

#if defined DYNAMIC_LANDCOVER_INPUT
		if (loadLU) {
			// Load landcover area fraction data from input file to data object
#ifdef LUTOMEMORY
			if (!LUdata_mem.Load(c)) {
#else
			if (!LUdata.Load(c)) {
#endif		
				dprintf("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;		// skip this stand
			}
		}
#endif
	}

	if(run[CROPLAND] && !LUerror)
	{
#if defined DYNAMIC_LANDCOVER_INPUT
		if(!cftfrac_fixed) {
			
			// Crop fraction data: read from crop fraction file; dynamic, so data for all years are loaded to CFTdata object and 
			// transferred to gridcell.cftfrac each year in getlandcover()			

#ifdef LUTOMEMORY
			if(!CFTdata_mem.Load(c)) {
#else
			if(!CFTdata.Load(c)) {
#endif		
				dprintf("Problems with CFT fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand
			}
		}

		if(forcesowingdates && !LUerror) { 
#ifdef LUTOMEMORY
			if(!sdates_mem.Load(c)) {
#else
			if(!sdates.Load(c)) {
#endif
				dprintf("Problems with sowing date input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand
			}
		}
		if(forceharvestdates && !LUerror) {
#ifdef LUTOMEMORY
			if(!hdates_mem.Load(c)) {
#else
			if(!hdates.Load(c)) {
#endif
				dprintf("Problems with harvest date input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand
			}
		}
		if(readNfert && !LUerror) {


#ifdef LUTOMEMORY
			if(!Nfert_mem.Load(c)) {
#else
			if(!Nfert.Load(c)) {
#endif
				dprintf("Problems with N fertilization input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand
			}
		}
#endif
	}

	return LUerror;
}

/// Retrieves nitrogen deposition for a particular gridcell
/** The values are either taken from a binary archive file or when it's not
 *  provided default to pre-industrial level of 2 kgN/ha/year.
 *
 *  The binary archive files have nitrogen deposition in gN/m2 on a monthly timestep
 *  for 16 years with 10 year interval starting from 1850 (Lamarque et. al., 2011).
 *
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 */
void CRUInput::getndep(double lon, double lat) {
	
	const double convert = 1e-7;				// converting from gN ha-1 to kgN m-2

	xtring file_ndep = param["file_ndep"].str;

	if (file_ndep == "") {

		// pre-industrial total nitrogen depostion set to 2 kgN/ha/year [kgN m-2]
		double dailyndep = 2000.0 / (4 * 365) * convert;

		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = dailyndep;
				NHxWetDep[y][m] = dailyndep;
				NOyDryDep[y][m] = dailyndep;
				NOyWetDep[y][m] = dailyndep;
			}
		}
	}
	else {
		GlobalNitrogenDepositionArchive ark;
		if (!ark.open(file_ndep)) {
			fail("Could not open %s for input", (char*)file_ndep);
		}

		GlobalNitrogenDeposition rec;
		rec.longitude = lon;
		rec.latitude = lat;

		if (!ark.getindex(rec)) {
			ark.close();
			fail("Grid cell not found in %s", (char*)file_ndep);
		}

		// Found the record, get the values
		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = rec.NHxDry[y*12+m] * convert;
				NHxWetDep[y][m] = rec.NHxWet[y*12+m] * convert;
				NOyDryDep[y][m] = rec.NOyDry[y*12+m] * convert;
				NOyWetDep[y][m] = rec.NOyWet[y*12+m] * convert;
			}
		}
		ark.close();
	}
}

/// Called by the framework at the start of the simulation for a particular grid cell
bool CRUInput::getgridcell(Gridcell& gridcell) {

	// DESCRIPTION
	// Obtains coordinates and soil static parameters for the next grid cell to
	// simulate. The function should return false if no grid cells remain to be simulated,
	// otherwise true. Currently the following member variables of Gridcell should be
	// initialised: longitude, latitude and climate.instype; the following members of
	// member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	// thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	// code (Sitch et al 2000) by a call to function soilparameters in the driver
	// module (driver.cpp):
	//
	// soilparameters(gridcell.soiltype,soilcode);
	//
	// If the model is to be driven by quasi-daily values of the climate variables
	// derived from monthly means, this function may be the appropriate place to
	// perform the required interpolations. The utility functions interp_monthly_means
	// and interp_monthly_totals in driver.cpp may be called for this purpose.

	// Select coordinates for next grid cell in linked list
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	bool gridfound;
	bool LUerror=false;

	// Make sure we use the first gridcell in the first call to this function,
	// and then step through the gridlist in subsequent calls.
	static bool first_call = true;

	if (first_call) {
		gridlist.firstobj();

		// Note that first_call is static, so this assignment is remembered
		// across function calls.
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		
		// guess2008 - New searchcru functions takee the CRU filenames as their first 
		// argument, i.e. cru_1901_2002.bin and cru_1901_2002_misc.bin

		// New code:

		double lon = gridlist.getobj().lon;
		double lat = gridlist.getobj().lat;
		gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
		                               hist_mtemp, hist_mprec, hist_msun);

		if (gridfound) // Get more historical CRU data for this grid cell
			gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);

		// Load environmental data for this grid cell from files
		if (run_landcover) {
			Coord& c=gridlist.getobj();
			LUerror=loadlandcover(gridcell, c);
		}
		if (LUerror)
			gridfound=false;

		while (!gridfound) {

			if (run_landcover && LUerror)
				dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", gridlist.getobj().lon,gridlist.getobj().lat);
			else
				dprintf("\nError: could not find stand at (%g,%g) in CRU data file\n", gridlist.getobj().lon,gridlist.getobj().lat);

			gridlist.nextobj();
			if (gridlist.isobj) {
				lon = gridlist.getobj().lon;
				lat = gridlist.getobj().lat;
				gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
				                               hist_mtemp, hist_mprec, hist_msun);
			  
				if (gridfound) // Get more historical CRU data for this grid cell
					gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation,
					                           hist_mfrs, hist_mwet, hist_mdtr);

				if (run_landcover) {
					Coord& c=gridlist.getobj();
					LUerror=loadlandcover(gridcell, c);
				}
				if (LUerror)
					gridfound=false;
			}
			else return false;
		}

		// Give sub-classes a chance to modify the data
		adjust_raw_forcing_data(gridlist.getobj().lon,
		                        gridlist.getobj().lat,
		                        hist_mtemp, hist_mprec, hist_msun);

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);
		spinup_mdtr.detrend_data();


		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Get nitrogen deposition data
		getndep(lon, lat);

		// Set CFT-specific members of climate and gridcellpft: 
		if (run_landcover && run[CROPLAND]) {
			if (gridcell.get_lat()>=0) 
			{
				gridcell.climate.testday_temp=180;		//June 30(day 180)
				gridcell.climate.testday_prec=364;		//Dec.31(day 364)
				gridcell.climate.coldestday=COLDEST_DAY_NHEMISPHERE;
				gridcell.climate.adjustlat=0;
			}
			else {
				gridcell.climate.testday_temp=364;		//Dec.31(day 364)
				gridcell.climate.testday_prec=180;		//June 30(day 180)
				gridcell.climate.coldestday=COLDEST_DAY_SHEMISPHERE;
				gridcell.climate.adjustlat=181;
			}

			for(unsigned int p = 0; p < gridcell.pft.nobj; p++) {
				Gridcellpft& gcpft=gridcell.pft[p];

				if (gridcell.get_lat()>=0.0) {
					gcpft.sdate_default=gcpft.pft.sdatenh;
					gcpft.hlimitdate_default=gcpft.pft.hlimitdatenh;
				}
				else {
					gcpft.sdate_default=gcpft.pft.sdatesh;
					gcpft.hlimitdate_default=gcpft.pft.hlimitdatesh;
				}
				// double cropping in China and Japan.
				if (!strncmp(gcpft.pft.name,"TrRi", strlen("TrRi")) && gridcell.get_lon()>=60.0 && gridcell.get_lat()<=30.0)
					gcpft.singlecrop=false;
			}
		}
	
		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		gridcell.climate.instype=SUNSHINE;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)
		
		clear_all_graphs();

		return true; // simulate this stand
	}

	return false; // no more stands
}

///	Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
void CRUInput::getlandcover(Gridcell& gridcell) {
	int i, year, year_saved;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;

	if(date.year<nyear_spinup)					// Use values for first historic year during spinup period.
		year=0;
	else if(date.year>=nyear_spinup+NYEAR_LU) {	// scenario adaptation
		year=NYEAR_LU-1;
	}
	else
		year=date.year-nyear_spinup;

	if(fixedlu_hist) {
		year_saved=year;
		year=0;
	}

	if(lcfrac_fixed) {		// If landcover area fractions are set in the ins-file.
		if(date.year==0) {	// Year 0: called by landcover_init
	
			int nactive_landcovertypes=0;

			if(equal_landcover_area) {
				for(i=0;i<NLANDCOVERTYPES;i++) {
					if(run[i])
						nactive_landcovertypes++;
				}
			}

			for(i=0;i<NLANDCOVERTYPES;i++) {
				if(equal_landcover_area) {
					sum_active+=gridcell.landcoverfrac[i]=1.0*run[i]/(double)nactive_landcovertypes;	// only set fractions that are active
					sum_tot=sum_active;
				}
				else {
#if defined GRASSFORCROP
					if(i==PASTURE)
						sum_tot+=gridcell.landcoverfrac[PASTURE]=(double)lc_fixed_frac[CROPLAND]/100.0;
					else if(i==CROPLAND)
						gridcell.landcoverfrac[CROPLAND]=0.0;
					else
#endif
					sum_tot+=gridcell.landcoverfrac[i]=(double)lc_fixed_frac[i]/100.0;					// count sum of all fractions (should be 1.0)

					if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0) {				// discard unreasonable values					
						if(date.year==0)
							dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
						sum_tot-=gridcell.landcoverfrac[i];
						gridcell.landcoverfrac[i]=0.0;
					}

					sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];				// only set fractions that are active
				}
			}
			
			if(sum_tot<0.99 || sum_tot>1.01) {	// Check input data, rescale if sum !=1.0
			
				sum_active=0.0;					// reset sum of active landcover fractions
				if(date.year==0)
					dprintf("WARNING ! landcover fixed fraction sum is %4.2f, rescaling landcover fractions !\n", sum_tot);

				for(i=0;i<NLANDCOVERTYPES;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}

			// NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
			if(sum_active<0.99)	{	// if landcover types are turned off in the ini-file, always <=1.0 here		
				if(date.year==0)
					dprintf("WARNING ! landcover active fraction sum is %4.2f.\n", sum_active);

				if(run[NATURAL]) {	// Transfer landcover areas not simulated to NATURAL fraction, if simulated.				
					if(date.year==0)
						dprintf("Inactive fractions (%4.2f) transferred to NATURAL fraction.\n", 1.0-sum_active);

					gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference 1.0-(sum of active landcover fractions) are added to the natural fraction
				}
				else {
/*					if(date.year==0)
						dprintf("Rescaling landcover fractions !\n");
					for(i=0;i<NLANDCOVERTYPES;i++)
						gridcell.landcoverfrac[i]/=sum_active;					// if NATURAL not simulated, rescale active fractions to 1.0
*/					if(date.year==0)
						dprintf("Non-unity fraction sum retained.\n");			// OR let sum remain non-unity
				}																
			}
		}
	}
	else {	// landcover area fractions are read from input file(s)	
#if defined DYNAMIC_LANDCOVER_INPUT
		bool getLU = false;

		for(i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				getLU = true;
		}

		if(getLU) {	
			// To allow run without LU data in landcover file (sets NATURAL to 1.0)
#ifdef LUTOMEMORY
			if(LUdata_mem.Get(year,0)==-9.999) {
#else
			if(LUdata.Get(year,0)==-9.999) {
#endif		
				dprintf("WARNING ! missing landcover fraction data for year %d, natural vegetation fraction set to 1.0\n", year+FIRSTHISTYEAR);
				memset(gridcell.landcoverfrac, 0, sizeof(double)*NLANDCOVERTYPES);
				gridcell.landcoverfrac[NATURAL]=1.0;
				sum_active=1.0;	//fix 110316
			}
			else {
				for(i=0;i<NLANDCOVERTYPES;i++)	{				
#ifdef LUTOMEMORY
					sum_tot+=gridcell.landcoverfrac[i]=LUdata_mem.Get(year,i);					// count sum of all fractions (should be 1.0)
#else
					sum_tot+=gridcell.landcoverfrac[i]=LUdata.Get(year,i);						// count sum of all fractions (should be 1.0)
#endif
				}
#ifdef GRASSFORCROP
				gridcell.landcoverfrac[PASTURE]+=gridcell.landcoverfrac[CROPLAND];
				gridcell.landcoverfrac[CROPLAND]=0.0;
#endif
				for(i=0;i<NLANDCOVERTYPES;i++)	{				
					if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0) {		// discard unreasonable values							
						if(date.year==0)
							dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
						sum_tot-=gridcell.landcoverfrac[i];
						gridcell.landcoverfrac[i]=0.0;
					}

					sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];
				}

				if(sum_tot!=1.0) {		// Check input data, rescale if sum !=1.0	

					sum_active=0.0;		// reset sum of active landcover fractions

					if(sum_tot<0.99 || sum_tot>1.01) {
						if(date.year==0) {
							dprintf("WARNING ! landcover fraction sum is %4.2f for year %d\n", sum_tot, year+FIRSTHISTYEAR);
							dprintf("Rescaling landcover fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);
						}
					}
					else				// sum often !=1.0 in input file
						if(!SUPPRESSLARGEOUTPUT)
							dprintf("Rescaling landcover fractions year %d ! (sum is within 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);

					for(i=0;i<NLANDCOVERTYPES;i++)
						sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
				}
			}
		}
		else
			gridcell.landcoverfrac[NATURAL]=0.0;

		// NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
		if(sum_active!=1.0)	{	// if landcover types are turned off in the ini-file, or if more landcover types are added in other input files, can be either less or more than 1.0
			if(!SUPPRESSLARGEOUTPUT)
				if(date.year==0)
					dprintf("Landcover fraction sum not 1.0 !\n");

			if(run[NATURAL]) {	// Transfer landcover areas not simulated to NATURAL fraction, if simulated.		
				if(date.year==0) {
					if(sum_active<1.0)
						dprintf("Inactive fractions (%4.3f) transferred to NATURAL fraction.\n", 1.0-sum_active);
					else
						dprintf("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).\n", sum_active-1.0, gridcell.landcoverfrac[NATURAL]);
				}

				gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					dprintf("New NATURAL fraction is %4.3f.\n", gridcell.landcoverfrac[NATURAL]);

				sum_active=1.0;		// sum_active should now be 1.0

				if(gridcell.landcoverfrac[NATURAL]<0.0) {	// If new landcover type fraction is bigger than the natural fraction (something wrong in the distribution of input file area fractions)						
					if(date.year==0)
						dprintf("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.\n");

					sum_active-=gridcell.landcoverfrac[NATURAL];	// fraction not possible to transfer moved back to sum_active, which will now be >1.0 again
					gridcell.landcoverfrac[NATURAL]=0.0;

					for(i=0;i<NLANDCOVERTYPES;i++) {
						gridcell.landcoverfrac[i]/=sum_active;		// fraction rescaled to unity sum
						if(run[i])
							if(date.year==0)
								dprintf("Landcover type %d fraction is %4.3f\n", i, gridcell.landcoverfrac[i]);
					}
				}
			}
			else {
//				if(date.year==0)
//					dprintf("Rescaling landcover fractions !\n");
//				for(i=0;i<NLANDCOVERTYPES;i++)
//					gridcell.landcoverfrac[i]/=sum_active;						// 1) if NATURAL not simulated, rescale active fractions to 1.0
				if(date.year==0)
					dprintf("Non-unity fraction sum retained.\n");				// 2) let sum remain non-unity
			}
		}
#endif
	}

	if(run[CROPLAND]) {
		sum=0.0;
		if(cftfrac_fixed) {		// If static equal crop fractions
			if(date.year==0) {	// Year 0: called by landcover_init
				for(int i=0;i<npft;i++)	{
					int index=-9;

					if(pftlist[i].cftid>=0)	{
						index=pftlist[i].cftid;
						sum+=gridcell.cftfrac[index]=1.0/(double)ncft;	
					}
				}
			}
		}
		else {					// crop area fractions are read from input file(s)
#if defined DYNAMIC_LANDCOVER_INPUT
			if(fixedcrop_hist)
				year=0;
			else if(fixedlu_hist)
				year=year_saved;
#ifdef LUTOMEMORY
			if(CFTdata_mem.Get(year,0)==-9.999) {	// to cope with missing Bondeau fraction data
#else
			if(CFTdata.Get(year,0)==-9.999) {		// to cope with missing Bondeau fraction data
#endif		
				dprintf("WARNING ! missing crop fraction data  for year %d, all set to 0.0\n", year+FIRSTHISTYEAR);
				memset(gridcell.cftfrac, 0, sizeof(double)*NCROPSTANDS_MAX);
			}
			else {
				// sum fractions for active crop pft:s and discard unreasonable values
				for(i=0;i<NCROPSTANDS_MAX;i++) {
					if(CFTdata.active[i]) {		// forces rescaling of fractions of active pft:s					
#ifdef LUTOMEMORY
						sum+=gridcell.cftfrac[i]=CFTdata_mem.Get(year,i);
#else
						sum+=gridcell.cftfrac[i]=CFTdata.Get(year,i);
#endif
						if(gridcell.cftfrac[i]<0.0 || gridcell.cftfrac[i]>1.0) {
							dprintf("WARNING ! crop fraction size out of limits, set to 0.0\n");
							sum-=gridcell.cftfrac[i];
							gridcell.cftfrac[i]=0.0;
						}
					}
				}
			}
#endif
		}

		if(!cftfrac_fixed || date.year==0) {
#if defined DYNAMIC_LANDCOVER_INPUT
			if(gridcell.landcoverfrac[CROPLAND]==0.0) {
				if(sum!=0.0) {
					dprintf("WARNING ! crop landcover fraction is 0.0 for year %d while crop data exist !\n", year+FIRSTHISTYEAR);
				}
			}
			else {
				if(sum==0.0) {
					if(!SUPPRESSLARGEOUTPUT)
						dprintf("WARNING ! crop fraction sum is 0.0 for year %d while LU[CROPLAND] is > 0 !\n", year+FIRSTHISTYEAR);

					//	Set to most common crop according to Bondeau
					pftlist.firstobj();
					while(pftlist.isobj) {
						Pft& pft=pftlist.getobj();
						if(pft.landcover==CROPLAND)	{
							
							if(!strcmp(pft.name,"TeWW") && (gridcell.get_lat()>30 || gridcell.get_lat()<-30)) {
								gridcell.cftfrac[pft.cftid]=1.0;				
								dprintf("Wheat fraction set to 1.0.\n");
							}
							else if(!strcmp(pft.name,"TrMi") && (gridcell.get_lat()<=30 && gridcell.get_lat()>=-30)) {
								gridcell.cftfrac[pft.cftid]=1.0;				
								dprintf("Millet fraction set to 1.0.\n");
							}
						}
						pftlist.nextobj();	
					}
					
				}
				else {
					// rescale active crop fraction so sum is 1.0
					for(i=0;i<NCROPSTANDS_MAX;i++)
						gridcell.cftfrac[i]/=sum;

					if(sum<0.99 || sum>1.01) {	// warn if sum is significantly different from 1.0 
						if(!SUPPRESSLARGEOUTPUT) {
							dprintf("WARNING ! crop fraction sum is %5.3f for year %d\n", sum, date.year-nyear_spinup+FIRSTHISTYEAR);
							dprintf("Rescaling crop fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);
						}
					}
				}
			}
#endif
		}
	}
}

/// Get sowing dates for one year
void CRUInput::getsowingdates(Gridcell& gridcell) {
	int i, year;

	if(date.year < nyear_spinup)
		year=0;
	else
		year = date.year - nyear_spinup;

	if(date.year < nyear_spinup + NYEAR_HIST) {
		for(i=0; i<npft; i++) {
			if(pftlist[i].cftid >= 0 && pftlist[i].forcesowingdate)	{ //natural pft:s have cftid=-1	
#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
				gridcell.pft[i].sdate_force = (int)sdates_mem.Get(year,pftlist[i].name);
#else
				gridcell.pft[i].sdate_force = (int)sdates.Get(year,pftlist[i].name);
#endif
#endif
			}
		}
	}
}

/// Get harvest dates for one year
void CRUInput::getharvestdates(Gridcell& gridcell) {
	int i, year;

	if(date.year < nyear_spinup)
		year=0;
	else
		year = date.year - nyear_spinup;

	if(date.year < nyear_spinup + NYEAR_HIST) {
 		for(i=0; i<npft; i++)	{
			if(pftlist[i].cftid >= 0 && pftlist[i].forceharvestdate) {	//natural pft:s have cftid=-1			
#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
				gridcell.pft[pftlist[i].id].hdate_force = (int)hdates_mem.Get(year,pftlist[i].name);
#else
				gridcell.pft[pftlist[i].id].hdate_force = (int)hdates.Get(year,pftlist[i].name);
#endif
#endif
			}
		}
	}
}

/// Get N fertilization for one year
void CRUInput::getNfert(Gridcell& gridcell) {
	int i, year;

	if(date.year < nyear_spinup)
		year=0;
	else
		year = date.year - nyear_spinup;

	if(date.year < nyear_spinup + NYEAR_HIST) {
 		for(i=0; i<npft; i++)	{
			if(pftlist[i].cftid >= 0 && pftlist[i].readNfert) {	//natural pft:s have cftid=-1			
#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
				gridcell.pft[pftlist[i].id].Nfert_read = Nfert_mem.Get(year,pftlist[i].name);
#else
				gridcell.pft[pftlist[i].id].Nfert_read = Nfert.Get(year,pftlist[i].name);
#endif
#endif
			}
		}
	}
}

/// Called by the framework each simulation day before any process modelling is performed for this day
/** Obtains climate data (including atmospheric CO2 and insolation) for this day. */
bool CRUInput::getclimate(Gridcell& gridcell) {

	// DESCRIPTION
	// The function should returns false if the simulation is complete for this grid cell,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of gridcell must be
	// initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	// values of the climate variables derived from monthly means, this day's values
	// will presumably be extracted from arrays containing the interpolated daily
	// values (see function getgridcell):
	//
	// gridcell.climate.temp=dtemp[date.day];
	// gridcell.climate.prec=dprec[date.day];
	// gridcell.climate.insol=dsun[date.day];
	// 
	// Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
	// BVOC:
	// gridcell.climate.dtr=ddtr[date.day]; 
	//
	// If model is run in diurnal mode, which requires appropriate climate forcing data, 
	// additional members of the climate must be initialised: temps, insols. Both of the
	// variables must be of type std::vector. The length of these vectors should be equal
	// to value of date.subdaily which also needs to be set either in getclimate or 
	// getgridcell functions. date.subdaily is a number of sub-daily period in a single 
	// day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
	// diurnal mode.

	double progress;

	Climate& climate = gridcell.climate;

	if (date.day == 0) {

		// First day of year ...
	
		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		int first_ndep_year = nyear_spinup + FIRSTHISTYEARNDEP - FIRSTHISTYEAR;

		double mndrydep[12], mnwetdep[12];
		int ndep_year = 0;

		if (date.year >= first_ndep_year) {
			ndep_year = (int)((date.year - first_ndep_year)/10);
		}

		for (int m = 0; m < 12; m++) {
			mndrydep[m] = NHxDryDep[ndep_year][m] + NOyDryDep[ndep_year][m];

			mnwetdep[m] = NHxWetDep[ndep_year][m] + NOyWetDep[ndep_year][m];
		}

		if (date.year < nyear_spinup) {

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

				// guess2008
				mfrs[m] = spinup_mfrs[m];
				mwet[m] = spinup_mwet[m];
				mdtr[m] = spinup_mdtr[m];
			}

			// Interpolate monthly spinup data to quasi-daily values
			interp_climate(mtemp,mprec,msun,mdtr,dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(mprec, dprec, mwet, gridcell.seed);
			}

			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);

			spinup_mtemp.nextyear();
			spinup_mprec.nextyear();
			spinup_msun.nextyear();

			// guess2008
			spinup_mfrs.nextyear();
			spinup_mwet.nextyear();
			spinup_mdtr.nextyear();

		}
		else if (date.year < nyear_spinup + NYEAR_HIST) {

			// Historical period

			// Interpolate this year's monthly data to quasi-daily values
			interp_climate(hist_mtemp[date.year-nyear_spinup],
				hist_mprec[date.year-nyear_spinup],hist_msun[date.year-nyear_spinup],
					   hist_mdtr[date.year-nyear_spinup],
				       dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year-nyear_spinup], dprec, hist_mwet[date.year-nyear_spinup], gridcell.seed);
			}

			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);

			// save climate date for calculation of seasonality
			for(int m=0;m<12;m++) {
				climate.mtemp_year[m] = hist_mtemp[date.year-nyear_spinup][m];
				climate.mprec_year[m] = hist_mprec[date.year-nyear_spinup][m];
			}
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}
	}

	if(date.day == 0) {
		climate.aprec = 0.0;

		for(int m=0; m<12; m++) {
			climate.aprec += climate.mprec_year[m];
			climate.mpet_year[m] = 0.0;
		}
	}
	climate.mpet_year[date.month] += climate.eet * PRIESTLEY_TAYLOR;

	// Send environmental values for today to framework
	if(fixedco2_hist)
		climate.co2 = co2[FIRSTHISTYEAR];
	else
		climate.co2 = co2[FIRSTHISTYEAR + date.year - nyear_spinup];

	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dsun[date.day];

	// Nitrogen deposition
	climate.dndep = dndep[date.day];

	// Nitrogen fertilization
	climate.dnfert = 0.0;

	// bvoc
	if(ifbvoc){
	  climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+NYEAR_HIST)
				+date.year)/(double)(ngridcell*(nyear_spinup+NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


CRUInput::~CRUInput() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
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
