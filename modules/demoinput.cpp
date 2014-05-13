///////////////////////////////////////////////////////////////////////////////////////
/// \file demoinput.cpp
/// \brief LPJ-GUESS input module for a toy data set (for demonstration purposes)
///
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "demoinput.h"

#include "driver.h"
#include "outputchannel.h"
#include <plib.h>
#include <stdio.h>

REGISTER_INPUT_MODULE("demo", DemoInput)

// Anonymous namespace for variables and functions with file scope
namespace {

// File names for temperature, precipitation, sunshine and soil code driver files
xtring file_temp,file_prec,file_sun,file_soil;

// LPJ soil code
int soilcode;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double mtemp[12], double mprec[12], double msun[12], double mdtr[12],
					double dtemp[365], double dprec[365], double dsun[365], double ddtr[365]) {
	interp_monthly_means_conserve(mtemp, dtemp);
	interp_monthly_totals_conserve(mprec, dprec, 0);
	interp_monthly_means_conserve(msun, dsun, 0, 100);
	interp_monthly_means_conserve(mdtr, ddtr, 0);
}

} // namespace


DemoInput::DemoInput() 
	: nyear(1),
	  lc_fixed_frac(NLANDCOVERTYPES, 0),
	  equal_landcover_area(false) {

	// Declare instruction file parameters

	declare_parameter("nyear", &nyear, 1, 10000, "Number of simulation years to run after spinup");

	declare_parameter("equal_landcover_area", &equal_landcover_area, "Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
	declare_parameter("minimizecftlist", &minimizecftlist, "Whether pfts not in crop fraction input file are removed from pftlist (0,1)");
	declare_parameter("lc_fixed_urban", &lc_fixed_frac[URBAN], 0, 100, "% lc_fixed_urban");
	declare_parameter("lc_fixed_cropland", &lc_fixed_frac[CROPLAND], 0, 100, "% lc_fixed_cropland");
	declare_parameter("lc_fixed_pasture", &lc_fixed_frac[PASTURE], 0, 100, "% lc_fixed_pasture");
	declare_parameter("lc_fixed_forest", &lc_fixed_frac[FOREST], 0, 100, "% lc_fixed_forest");
	declare_parameter("lc_fixed_natural", &lc_fixed_frac[NATURAL], 0, 100, "% lc_fixed_natural");
	declare_parameter("lc_fixed_peatland", &lc_fixed_frac[PEATLAND], 0, 100, "% lc_fixed_peatland");
}


void DemoInput::read_from_file(Coord coord, xtring fname, const char* format,
                               double monthly[12], bool soil /* = false */) {
	double dlon, dlat;
	int elev;
	FILE* in = fopen(fname, "r");
	if (!in) {
		fail("readenv: could not open %s for input", (char*)fname);
	}

	bool foundgrid = false;
	while (!feof(in) && !foundgrid) {
		if (!soil) {
			readfor(in, format, &dlon, &dlat, &elev, monthly);
		} else {
			readfor(in, format, &dlon, &dlat, &soilcode);
		}
		foundgrid = equal(coord.lon, dlon) && equal(coord.lat, dlat);
	}


	fclose(in);
	if (!foundgrid) {
		fail("readenv: could not find record for (%g,%g) in %s",
										coord.lon, coord.lat, (char*)fname);
	}
}

bool DemoInput::readenv(Coord coord, long& seed) {

	// Searches for environmental data in driver temperature, precipitation,
	// sunshine and soil code files for the grid cell whose coordinates are given by
	// 'coord'. Data are written to arrays mtemp, mprec, msun and the variable
	// soilcode, which are defined as global variables in this file

	// The temperature, precipitation and sunshine files (Cramer & Leemans,
	// unpublished) should be in ASCII text format and contain one-line records for the
	// climate (mean monthly temperature, mean monthly percentage sunshine or total
	// monthly precipitation) of a particular 0.5 x 0.5 degree grid cell. Elevation is
	// also included in each grid cell record.

	// The following sample record from the temperature file:
	//   " -4400 8300 293-298-311-316-239-105  -7  26  -3 -91-184-239-277"
	// corresponds to the following data:
	//   longitude 44 deg (-=W)
	//   latitude 83 deg (+=N)
	//   elevation 293 m
	//   mean monthly temperatures (deg C) -29.8 (Jan), -31.1 (Feb), ..., -27.7 (Dec)

	// The following sample record from the precipitation file:
	//   " 12750 -200 223 190 165 168 239 415 465 486 339 218 162 149 180"
	// corresponds to the following data:
	//   longitude 127.5 deg (+=E)
	//   latitude 20 deg (-=S)
	//   elevation 223 m
	//   monthly precipitation sum (mm) 190 (Jan), 165 (Feb), ..., 180 (Dec)

	// The following sample record from the sunshine file:
	//   "  2600 7000 293  0 20 38 37 31 28 28 25 21 17  9  7"
	// corresponds to the following data:
	//   longitude 26 deg (+=E)
	//   latitude 70 deg (+=N)
	//   elevation 293 m
	//   monthly mean %age of full sunshine 0 (Jan), 20 (Feb), ..., 7 (Dec)

	// The LPJ soil code file is in ASCII text format and contains one-line records for
	// each grid cell in the form:
	//   <lon> <lat> <soilcode>
	// where <lon>      = longitude as a floating point number (-=W, +=E)
	//       <lat>      = latitude as a floating point number (-=S, +=N)
	//       <soilcode> = integer in the range 0 (no soil) to 9 (see function
	//                    soilparameters in driver module)
	// The fields in each record are separated by spaces

	double mtemp[12];		// monthly mean temperature (deg C)
	double mprec[12];		// monthly precipitation sum (mm)
	double msun[12];		// monthly mean percentage sunshine values

	double mwet[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // number of rain days per month

	double mdtr[12];		// monthly mean diurnal temperature range (oC)
	for(int m=0; m<12; m++) {
		mdtr[m] = 0.;
		if (ifbvoc) {
			dprintf("WARNING: No data available for dtr in sample data set!\nNo daytime temperature correction for BVOC calculations applied.");
		}
	}

	read_from_file(coord, file_temp, "f6.2,f5.2,i4,12f4.1", mtemp);
	read_from_file(coord, file_prec, "f6.2,f5.2,i4,12f4", mprec);
	read_from_file(coord, file_sun, "f6.2,f5.2,i4,12f3", msun);
	read_from_file(coord, file_soil, "f,f,i", msun, true);	// msun is not used here: just dummy

	// Interpolate monthly values for environmental drivers to daily values
	// (relevant daily values will be sent to the framework each simulation
	// day in function getclimate, below)
	interp_climate(mtemp, mprec, msun, mdtr, dtemp, dprec, dsun, ddtr);

	// Recalculate precipitation values using weather generator
	// (from Dieter Gerten 021121)
	prdaily(mprec, dprec, mwet, seed);

	return true;
}


void DemoInput::init() {

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

	// Retrieve specified CO2 value as read from ins file
	co2=param["co2"].num;

	// Retrieve specified N value as read from ins file
	ndep=param["ndep"].num;

	if (run_landcover) {
		all_fracs_const=true;	// If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year

		// Retrieve file names for landcover files and open them if static values from ins-file are not used.
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
				else if(LUdata.format==InData::LOCAL_YEARLY)
					all_fracs_const=false;				// Set all_fracs_const to false if yearly data
#ifdef LUTOMEMORY
					// Save all landcover area fraction data in memory
					LUdata_mem.Open(gridlist.nobj, LUdata.nRecords,LUdata.nYears);

					ListArray_id<InData::Coord> lonlatlist;
					GetLonLatList(lonlatlist, gridlist);
					LUdata_mem.CopyFromTimeDataD(LUdata, lonlatlist);
#endif
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
			else {

				bool do_minimize = false;

				if(minimizecftlist && CFTdata.GetNCells() < 1000) {	// Reduce the risk of accidentally using minimized cft lists when using split gridlists.
					// remove all crop pft:s from gridlist that always have zero area fraction
					ListArray_id<InData::Coord> lonlatlist;
					GetLonLatList(lonlatlist, gridlist);
					CFTdata.CheckIfPresent(lonlatlist);
					do_minimize = true;
				}
				
				int n=0;
				pftlist.firstobj();
				while(pftlist.isobj) {	

					bool remove = false;

					if(pftlist.getobj().cftid>=0) {
						if(do_minimize)
							remove = !CFTdata.item_has_data(pftlist.getobj().name);
						else
							remove = !CFTdata.item_in_header(pftlist.getobj().name);
					}

					if(remove && !(pftlist.getobj().isintercropgrass && ifintercropgrass)) {
						n+=1;
						pftlist.killobj();
						npft--;
						ncft--;
					}
					else {
						pftlist.getobj().id-=n;
						pftlist.nextobj();
					}			
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
			if(CFTdata.GetnRecords() > NCROPSTANDS_MAX)
				fail("\ninitio: NCROPSTANDS_MAX is incorrectly set in guess.h !\n");
#endif
		}
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

	// Retrieve input file names as read from ins file

	file_temp=param["file_temp"].str;
	file_prec=param["file_prec"].str;
	file_sun=param["file_sun"].str;
	file_soil=param["file_soil"].str;

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}


#if defined DYNAMIC_LANDCOVER_INPUT

/// Transfers coordinates from DemoInput::Coord to InData::Coord
InData::Coord DemoInput::GetLonLat(Coord coord) {

	InData::Coord lonlat;
	lonlat.lon=coord.lon;
	lonlat.lat=coord.lat;

	return lonlat;
}

/// Transfers gridlist of coordinates from DemoInput::Coord to InData::Coord
void DemoInput::GetLonLatList(ListArray_id<InData::Coord>&lonlatlist, ListArray_id<Coord>& gridlist) {

	for(unsigned int i = 0; i < gridlist.nobj; i++) {
		InData::Coord& c= lonlatlist.createobj();
		c.lon=gridlist[i].lon;
		c.lat=gridlist[i].lat;
	}
}

#endif

bool DemoInput::loadlandcover(Gridcell& gridcell, Coord cc)	{

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

	if(run[CROPLAND] && !LUerror) {
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
//				dprintf("Problems with N fertilization input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
//				LUerror=true;	// skip this stand
				dprintf("N fertilization data not found in input file for %.2f,%.2f.\n\n",c.lon,c.lat);
			}
		}
#endif
	}

	return LUerror;
}


bool DemoInput::getgridcell(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	// Select coordinates for next grid cell in linked list
	bool gridfound = false;

	bool LUerror = false;

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

		while(!gridfound) {

			// Retrieve coordinate of next grid cell from linked list
			Coord& c = gridlist.getobj();

			// Load environmental data for this grid cell from files
			if(run_landcover) {
				LUerror = loadlandcover(gridcell, c);
			}
			if (!LUerror) {
				gridfound = readenv(c, gridcell.seed);
			} else {
				gridlist.nextobj();
			}
		}


		dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
			gridlist.getobj().lat);
		if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
			(char*)gridlist.getobj().descrip);
		else dprintf("\n\n");
		
		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// Set CFT-specific members of climate and gridcellpft: 
		if (run_landcover && run[CROPLAND]) {
			if (gridcell.get_lat()>=0) {
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


void DemoInput::getlandcover(Gridcell& gridcell) {
	int i, year;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;

	if(date.year<nyear_spinup)					// Use values for first historic year during spinup period.
		year=0;
	else if(date.year>=nyear_spinup+NYEAR_LU) {	// scenario adaptation
		year=NYEAR_LU-1;
	}
	else
		year=date.year-nyear_spinup;

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
		memset(gridcell.cftfrac, 0, sizeof(double)*NCROPSTANDS_MAX);

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
#ifdef LUTOMEMORY
			if(CFTdata_mem.Get(year,0)==-9.999) {	// to cope with missing Bondeau fraction data
#else
			if(CFTdata.Get(year,0)==-9.999) {		// to cope with missing Bondeau fraction data
#endif		
				dprintf("WARNING ! missing crop fraction data  for year %d, all set to 0.0\n", year+FIRSTHISTYEAR);
			}
			else {
				// sum fractions for active crop pft:s and discard unreasonable values
				for(i=0; i<npft; i++) {
					if(pftlist[i].cftid >= 0)	{ //natural pft:s have cftid=-1	
#ifdef LUTOMEMORY
						sum += gridcell.cftfrac[pftlist[i].cftid] = CFTdata_mem.Get(year,pftlist[i].name);
#else
						sum += gridcell.cftfrac[pftlist[i].cftid] = CFTdata.Get(year,pftlist[i].name);
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

void DemoInput::getsowingdates(Gridcell& gridcell) {
	int i, year;

#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
	if(!sdates_mem.isloaded())
#else
	if(!sdates.isloaded())
#endif
		return;
#endif

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

void DemoInput::getharvestdates(Gridcell& gridcell) {
	int i, year;

#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
	if(!hdates_mem.isloaded())
#else
	if(!hdates.isloaded())
#endif
		return;
#endif

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
void DemoInput::getNfert(Gridcell& gridcell) {
	int i, year;

#if defined DYNAMIC_LANDCOVER_INPUT
#ifdef LUTOMEMORY
	if(!Nfert_mem.isloaded())
#else
	if(!Nfert.isloaded())
#endif
		return;
#endif

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


bool DemoInput::getclimate(Gridcell& gridcell) {

	// See base class for documentation about this function's responsibilities

	double progress;

	Climate& climate = gridcell.climate;


	// Send environmental values for today to framework

	climate.dndep  = ndep / (365.0 * 10000.0);
	climate.dnfert = 0.0;

	climate.co2=co2;

	climate.temp  = dtemp[date.day];
	climate.prec  = dprec[date.day];
	climate.insol = dsun[date.day];

	// bvoc

	climate.dtr=ddtr[date.day];


	// First day of year only ...

	if (date.day == 0) {

		// Return false if last year was the last for the simulation
		if (date.year==nyear_spinup+nyear) return false;

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+nyear)
				+date.year)/(double)(ngridcell*(nyear_spinup+nyear));


			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


DemoInput::~DemoInput() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
}
