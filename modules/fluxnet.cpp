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
REGISTER_OUTPUT_MODULE("fluxnet", FluxnetOutput)

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

	// Read site name description from gridlist file

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


		current_stand_fluxdata = &Fluxnetdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];
		//
		FluxnetData& fluxdata = Fluxnetdata[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];


		// *** Adjust all CRU temp and precip data to site conditions


		// Get actual temp and precip data for the site

		// Create some strings
		xtring fluxdirectory=param["flux_dir"].str;
		xtring fluxfilestart = "FLUXNET2015_FULLSET_";
		xtring fluxfileend = ".txt";

		xtring fluxfile =  fluxdirectory + fluxfilestart + fluxdata.desc + fluxfileend;


		dprintf("Fluxnet file: %s \n",(char*)fluxfile);

		FILE* in_flux=fopen(fluxfile,"r");
		if (!in_flux) fail("getgridcell: could not open %s for input",(char*)fluxfile);

		bool eof = false;
		xtring header;

		// Read the header first. We don't use this.
		eof=!readfor(in_flux,"a",&header);


		while (!eof) {

			double sitedata[5];
			eof=!readfor(in_flux,"f,f,f,f,f",
					&sitedata[0],&sitedata[1],&sitedata[2],&sitedata[3],&sitedata[4]);
			if (!eof) {

				if (fluxdata.start_y == 0) {
					fluxdata.start_y = sitedata[0];
				}



				// Read the relevant climate data from the site
				yr.push_back( sitedata[0]);	// Year
				Ta.push_back(sitedata[2]);	// air temp degC
				rain.push_back(sitedata[3]); // precip mm day-1
				swrad.push_back(sitedata[4]); // swrad W m-2


			}

		}

		fclose(in_flux);

		std::vector<double>::size_type days = rain.size();


		if (days % 365) {
			fail("Given time series doesn't extend for a full number of years (length: %d)\n", days);
		}


		fluxdata.end_y = yr[days-1];

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


/// OUTPUT MODULE FOR FLUXNET

FluxnetOutput::FluxnetOutput() {
	// Files for FLUXNET output
	declare_parameter("file_fluxnetdaily", &file_fluxnetdaily, 300, "FLUXNET daily output file");
	declare_parameter("file_fluxnetclim", &file_fluxnetclim, 300, "FLUXNET daily Climate output file");
	declare_parameter("file_fluxnetmonth", &file_fluxnetmonth, 300, "FLUXNET Monthly output file");
}

void FluxnetOutput::init() {
	// create the output tables


		ColumnDescriptors fluxnet_columns;
		fluxnet_columns += ColumnDescriptor("Veg", 14, 6);
		fluxnet_columns += ColumnDescriptor("Repr", 14, 6);
		fluxnet_columns += ColumnDescriptor("Soil", 14, 6);
		fluxnet_columns += ColumnDescriptor("Fire", 14, 6);
		fluxnet_columns += ColumnDescriptor("Est", 14, 6);
		fluxnet_columns += ColumnDescriptor("NEE", 14, 6);
		create_output_table(out_fluxnetdaily, file_fluxnetdaily, fluxnet_columns);


		ColumnDescriptors fluxnetclim_columns;
		fluxnetclim_columns += ColumnDescriptor("temp", 14, 6);
		fluxnetclim_columns += ColumnDescriptor("prec", 14, 6);
		fluxnetclim_columns += ColumnDescriptor("insol", 14, 6);
		create_output_table(out_fluxnetclim, file_fluxnetclim, fluxnetclim_columns);


		ColumnDescriptors fluxnetmonth_columns;
		fluxnetmonth_columns += ColumnDescriptor("Month", 8, 0);
		fluxnetmonth_columns += ColumnDescriptor("NEE", 10, 3);
		fluxnetmonth_columns += ColumnDescriptor("AET", 10, 3);
		fluxnetmonth_columns += ColumnDescriptor("GPP", 10, 3);
		create_output_table(out_fluxnetmonth, file_fluxnetmonth, fluxnetmonth_columns);


}


void FluxnetOutput::outannual(Gridcell& gridcell) {


}


void FluxnetOutput::outdaily(Gridcell& gridcell) {

	// DESCRIPTION
	// Output of simulation results at the end of each day
	// added by niklas

	double lon,lat;
	double flux_veg, flux_repr, flux_soil, flux_fire, flux_est, flux_seed, flux_charvest;
	lon=gridcell.get_lon();
	lat=gridcell.get_lat();
	//climate output

	if (date.get_calendar_year() >= (current_stand_fluxdata->start_y - 10)) {
		// The OutputRows object manages the next row of output for each
		// output table
		OutputRows out1(output_channel, lon, lat, date.get_calendar_year(), date.day);

		out1.add_value(out_fluxnetclim, gridcell.climate.temp);
		out1.add_value(out_fluxnetclim, gridcell.climate.prec);
		out1.add_value(out_fluxnetclim, gridcell.climate.insol);

	}


	//Flux output

	if (date.get_calendar_year() >= current_stand_fluxdata->start_y) {




		// The OutputRows object manages the next row of output for each
		// output table
		OutputRows out(output_channel, lon, lat, date.get_calendar_year(), date.day);


		flux_veg = flux_repr = flux_soil = flux_fire = flux_est = flux_seed = flux_charvest = 0.0;



		Gridcell::iterator gc_itr = gridcell.begin();

		// Loop through Stands
		while (gc_itr != gridcell.end()) {
			Stand& stand = *gc_itr;
			stand.firstobj();

			//Loop through Patches
			while (stand.isobj) {
				Patch& patch = stand.getobj();

				double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

				flux_veg+=-patch.fluxes.get_daily_flux(Fluxes::NPP,date.day)*to_gridcell_average;
				flux_repr+=-patch.fluxes.get_daily_flux(Fluxes::REPRC,date.day)*to_gridcell_average;
				flux_soil+=patch.fluxes.get_daily_flux(Fluxes::SOILC,date.day)*to_gridcell_average;
				flux_fire+=patch.fluxes.get_daily_flux(Fluxes::FIREC,date.day)*to_gridcell_average;
				flux_est+=patch.fluxes.get_daily_flux(Fluxes::ESTC,date.day)*to_gridcell_average;
				flux_seed+=patch.fluxes.get_daily_flux(Fluxes::SEEDC,date.day)*to_gridcell_average;
				flux_charvest+=patch.fluxes.get_daily_flux(Fluxes::HARVESTC,date.day)*to_gridcell_average;


				stand.nextobj();
			} // patch loop
			++gc_itr;
		} // stand loop

		out.add_value(out_fluxnetdaily, flux_veg);
		out.add_value(out_fluxnetdaily, -flux_repr);
		out.add_value(out_fluxnetdaily, flux_soil);
		out.add_value(out_fluxnetdaily, flux_fire);
		out.add_value(out_fluxnetdaily, flux_est);

		// daily NEE do not include fire, establishment as monthly NEE
		out.add_value(out_fluxnetdaily, flux_veg   +  flux_soil );


		if(date.islastday){ //output monthly values
			OutputRows out2(output_channel, gridcell.get_lon(), gridcell.get_lat(), date.get_calendar_year());

			double mnee, mgpp,maet;
			mnee = mgpp = maet = 0.0;
			Gridcell::iterator gc_itr = gridcell.begin();

			while (gc_itr != gridcell.end()) {

				Stand& stand = *gc_itr;
				stand.firstobj();

				while (stand.isobj) {
					Patch& patch = stand.getobj();

					double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();


						double gpp = patch.fluxes.get_monthly_flux(Fluxes::GPP, date.month);
						double ra = patch.fluxes.get_monthly_flux(Fluxes::RA, 	date.month);
						double rh = patch.fluxes.get_monthly_flux(Fluxes::SOILC, date.month);
						double npp = gpp - ra;
						mnee += (rh - npp)*to_gridcell_average;
						mgpp += gpp*to_gridcell_average;
						maet += patch.maet[date.month]*to_gridcell_average;

					stand.nextobj();
				}

				++gc_itr;
			}

			 out2.add_value(out_fluxnetmonth,date.month+1);
			 out2.add_value(out_fluxnetmonth, -1000.0*mnee/date.dayofmonth);// gC/m2/day (average)
			 out2.add_value(out_fluxnetmonth,  maet); // mm/month
			 out2.add_value(out_fluxnetmonth, 1000.0*mgpp/date.dayofmonth);  // gC/m2/day (average)

		}



	} // end if date year > spinup year
}

