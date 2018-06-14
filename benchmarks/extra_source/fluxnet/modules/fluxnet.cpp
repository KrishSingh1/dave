///////////////////////////////////////////////////////////////////////////////////////
/// \file fluxnet.cpp
/// \brief Input and output modules for the fluxnet benchmarks
///
/// \author Niklas Boke Olén
/// $Date: 2015-11-13 16:25:45 +0100 (Fri, 13 Nov 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "fluxnet.h"
#include "driver.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>

namespace {
    xtring file_cru;
    xtring file_cru_misc;
}
/*
bool getgridcell(Gridcell& gridcell) {
    bool success = CRUInput::getgridcell(gridcell));
    
    if (!success){
		return false;
	}
	tair.clear();
	rain.clear();
	swrad.clear();

	xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip + ".csv";
	std::ifstream ifs(fluxfile, std::ifstream::in);

	if (!ifs.good()) {
		dprintf("FluxnetInput::getgridcell: could not open %s for input\n", (char*)fluxfile);
		return false;
	}

	std::string line;
	//getline(ifs, line);
	//std::istringstream iss(line);
    double lon, lat;
    lon = gridlist.getobj().lon;
    lat = gridlist.getobj().lat;

	while (getline(ifs, line)) {

		std::string monthday;
		double temp, prec, insol;

		std::istringstream iss(line);

		if (iss >> end_year >> monthday >> temp >> insol >> prec) {
			if (tair.empty()) {
				dailyoutput_firstyear = end_year;
			}
			tair.push_back(temp);
			swrad.push_back(insol);
			rain.push_back(prec);
		}
	}

	ifs.close();

	std::vector<double>::size_type days = rain.size();
	if (days % 365) {
		dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
		return false;
	}

	//gridcell.set_coordinates(lon, lat);

	return true;
}
*/

void FluxnetInput::init() {
    CRUInput::init();
    
    file_cru=param["file_cru"].str;
    file_cru_misc=param["file_cru_misc"].str;
}

void FluxnetInput::get_fluxnet_data_from_file() {
	
	std::string line;
	
	tair.clear();
	rain.clear();
	swrad.clear();
	
	xtring fluxfile = param["flux_dir"].str + gridlist.getobj().descrip + ".csv";
	std::ifstream ifs(fluxfile, std::ifstream::in);
	
	if (!ifs.good()) {
		dprintf("FluxnetInput::getgridcell: could not open %s for input\n", (char*)fluxfile);
		fail();
	}
	
	while (getline(ifs, line)) {
		
		std::string monthday;
		double temp, prec, insol;
		
		std::istringstream iss(line);
		
		if (iss >> firstyear >> monthday >> temp >> insol >> prec) {
			if (tair.empty()) {
				first_fluxnet_year = firstyear;
			}
			tair.push_back(temp);
			swrad.push_back(insol);
			rain.push_back(prec);
		}
	}
	
	ifs.close();
	
	std::vector<double>::size_type days = rain.size();
	if (days % 365) {
		dprintf("Given time series doesn't extend for a full number of years (length: %d)\n", days);
		fail();
	}
}

bool FluxnetInput::getgridcell(Gridcell& gridcell) {
    // See base class for documentation about this function's responsibilities
        
    int soilcode;
    int elevation;
        
    // Make sure we use the first gridcell in the first call to this function,
    // and then step through the gridlist in subsequent calls.
    if (first_call) {
        gridlist.firstobj();
            
        // Note that first_call is static, so this assignment is remembered
        // across function calls.
        first_call = false;
    }
    else gridlist.nextobj();
        
    if (gridlist.isobj) {
            
            bool gridfound = false;
            bool LUerror = false;
            double lon;
            double lat;
		
			while(!gridfound) {

				if(gridlist.isobj) {
                    
					lon = gridlist.getobj().lon;
                    lat = gridlist.getobj().lat;
                    gridfound = CRU_TS30::findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode,
                                                             hist_mtemp, hist_mprec, hist_msun);
                    
                    if (gridfound) // Get more historical CRU data for this grid cell
                        gridfound = CRU_TS30::searchcru_misc(file_cru_misc, lon, lat, elevation,
                                                             hist_mfrs, hist_mwet, hist_mdtr);
                    
                    if (run_landcover && gridfound) {
                        LUerror = landcover_input.loadlandcover(lon, lat);
                        if(!LUerror)
                            LUerror = management_input.loadmanagement(lon, lat);
                    }
                    
                    if(!gridfound || LUerror) {
                        if(!gridfound)
                            dprintf("\nError: could not find stand at (%g,%g) in climate data files\n", gridlist.getobj().lon,gridlist.getobj().lat);
                        else if(LUerror)
                            dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", gridlist.getobj().lon,gridlist.getobj().lat);
                        gridfound = false;
                        gridlist.nextobj();
                    }
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
		
			// Get the extra fluxnet data
			get_fluxnet_data_from_file();
		
            // We wont detrend dtr for now. Partly because dtr is at the moment only
            // used for BVOC, so what happens during the spinup is not affecting
            // results in the period thereafter, and partly because the detrending
            // can give negative dtr values.
            //spinup_mdtr.detrend_data();
		
            dprintf("\nCommencing simulation for stand at (%g,%g)",gridlist.getobj().lon,
                    gridlist.getobj().lat);
            if (gridlist.getobj().descrip!="") dprintf(" (%s)\n\n",
                                                       (char*)gridlist.getobj().descrip);
            else dprintf("\n\n");
            
            // Tell framework the coordinates of this grid cell
            gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);
            
            // Get nitrogen deposition data.
            /* Since the historic data set does not reach decade 2010-2019,
             * we need to use the RCP data for the last decade. */
            ndep.getndep(param["file_ndep"].str, lon, lat, Lamarque::RCP60);
            
            // The insolation data will be sent (in function getclimate, below)
            // as incoming shortwave radiation, averages are over 24 hours
            
            gridcell.climate.instype = SWRAD_TS;
            
            // Tell framework the soil type of this grid cell
            soilparameters(gridcell.soiltype,soilcode);
            
            clear_all_graphs();
            
            return true; // simulate this stand
    }
        
    return false; // no more stands
}

    
    

bool FluxnetInput::getclimate(Gridcell& gridcell) {
	
	if (!CRUInput::getclimate(gridcell)) {
		return false;
	};
	
	int year = date.get_calendar_year();
	
	// Overwrite / extend climate data with site values for the period
	int id = (year - first_fluxnet_year) * 365 + date.day;

	//if (date.day == 0) {
	//	double mndrydep[12], mnwetdep[12];
	//	ndep.get_one_calendar_year(year, mndrydep, mnwetdep);
	//	distribute_ndep(mndrydep, mnwetdep, &rain[id], dndep);
	//}

	gridcell.climate.prec = rain[id];
	gridcell.climate.temp = tair[id];
	gridcell.climate.insol = swrad[id];
	gridcell.climate.dndep = dndep[date.day];
	return true;
}


REGISTER_INPUT_MODULE("fluxnet", FluxnetInput)
