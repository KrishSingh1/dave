///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Extra code used by the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "fluxnet.h"
#include "parameters.h"
#include "guess.h"
#include "driver.h"

REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)


using namespace GuessOutput;

// guess2008 - Fluxnet - new int to keep track of the simulation year
// Needed for management etc., used in vegetation dynamics
// century_year = 0, when date.year < nyear, i.e. during spin up.
// century_year = 1, when date.year = nyear, i.e. 1901,
// century_year = 80, when date.year = nyear+79, i.e. 1980, etc.
int century_year;

FluxnetData* current_stand_fluxdata = 0;

void FluxnetInput::init() {
	// First let base class initialize
	CRUInput::init();

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	bool eof = false;

	while (!eof) {

		FluxnetData edata;

		double dlon, dlat;

		eof=!readfor(in_grid,"f,f,a;a",&dlon,&dlat,&edata.desc,&edata.desc2);



		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Fluxnetdata[std::make_pair(dlon, dlat)] = edata;
		}
	}

	fclose(in_grid);
}

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
	if (CRUInput::getgridcell(gridcell)) {


		//current_stand_fluxdata = &Fluxnetdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];
		//
		FluxnetData& coord = Fluxnetdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];


		// *** Adjust all CRU temp and precip data to site conditions


		// Get actual temp and precip data for the site, as well as NEE and latent heat flux

		// Create some strings
		xtring fluxdirectory=param["flux_dir"].str;
		xtring fluxfilestart = "FLUXNET2015_FULLSET_";
		xtring fluxfileend = ".txt";

		xtring fluxfile =  fluxdirectory + fluxfilestart + coord.desc + fluxfileend;


		dprintf("%s \n",(char*)fluxfile);

		FILE* in_flux=fopen(fluxfile,"r");
		if (!in_flux) fail("getgridcell: could not open %s for input",(char*)fluxfile);

		bool eof = false;
		xtring header;

		// Read the header first. We don't use this.
		eof=!readfor(in_flux,"a",&header);

		coord.start_y=0;
		coord.end_y= 0;


		while (!eof) {

			double sitedata[5];
			eof=!readfor(in_flux,"f,f,f,f,f",
					&sitedata[0],&sitedata[1],&sitedata[2],&sitedata[3],&sitedata[4]);
			if (!eof) {

				if (!coord.start_y) {
					coord.start_y = sitedata[0];
				}



				// Read the relevant climate data from the site
				yr.push_back( sitedata[0]);	// Year
				Ta.push_back(sitedata[2]);	// degC
				rain.push_back(sitedata[3]); // mm day-1
				swrad.push_back(sitedata[4]); // W m-2


			}




		}


		fclose(in_flux);

		std::vector<double>::size_type days = rain.size();

		dprintf("%d \n",days);

		if (days % 365) {
			fail("Given timeseries doesn't extend for a full number of years (length: %d)\n", days);
		}


		coord.end_y = yr[days-1];

		dprintf("s=%d e=%d\n",coord.start_y,coord.end_y);

		ndep.getndep(param["file_ndep"].str, gridcell.get_lon(), gridcell.get_lat(), Lamarque::RCP60);

		return true;
	}
	else {
		return false;
	}
}

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	if (date.day == 0) {
		if (date.year < nyear_spinup) {
			century_year = 0;
		}
		else {
			century_year++;
		}
	}

	//HÄR SKA vi uppdatera daily climate om vi har fluxnet år., kanske bättre att läsa in data i init och spara till. return false om fluxnetdata är slut annars true,
	//uppdatera även ndep så som det görs i cru getclimate.
	//kan skapa vektorer dynamisk så som görs i ascii input tex i expeer branch.

	bool cru_input = CRUInput::getclimate(gridcell);


	//Update with fluxnet daily data

	FluxnetData& fluxdata = Fluxnetdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];

	if(date.get_calendar_year() >= fluxdata.start_y && date.get_calendar_year() <= fluxdata.end_y){

		int id = (date.get_calendar_year()-fluxdata.start_y) * 365 + date.day;


		if(date.day == 0){
			double mndrydep[12], mnwetdep[12];
			ndep.get_one_calendar_year(date.get_calendar_year(), mndrydep, mnwetdep);
			distribute_ndep(mndrydep, mnwetdep, &rain[id], dndep);
		}


		gridcell.climate.prec = rain[id];
		gridcell.climate.temp = Ta[id];
		gridcell.climate.insol = swrad[id];
		gridcell.climate.dndep = dndep[date.day];


		return true;
	}else if(date.get_calendar_year() > fluxdata.end_y){
		return false;
	}

	return cru_input;



	//	return CRUInput::getclimate(gridcell);
}






