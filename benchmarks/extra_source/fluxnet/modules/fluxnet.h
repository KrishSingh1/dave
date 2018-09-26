//
//  fluxnet.h
//  guess
//
//  Created by Adrian Gustafson on 2018-06-14.
//
//

#ifndef LPJ_GUESS_FLUXNET_H
#define LPJ_GUESS_FLUXNET_H

#include "guess.h"
#include "cruinput.h"
#include <gutil.h>

class FluxnetInput : public CRUInput {
public:
    void init();
    
    xtring extracted();
    
    bool getgridcell(Gridcell& gridcell);
    
    bool getclimate(Gridcell& gridcell);
	
	bool get_fluxnet_data_from_file();

	void adjust_raw_forcing_data(double hist_mtemp[NYEAR_HIST][12], 
								 double hist_mprec[NYEAR_HIST][12],
								 double hist_msun[NYEAR_HIST][12], 
								 double fluxnet_temp[12],
								 double fluxnet_prec[12],
								 double fluxnet_sun[12]);
    
private:
	std::vector<double> rain, tair, swrad;
	int first_year, last_year, nyear;
	
	double monthly_fluxnet_temp[12], monthly_fluxnet_rain[12], monthly_fluxnet_rad[12];

	// Daily precipitation during one year. Needed for the N-deposition.
	double drain[Date::MAX_YEAR_LENGTH];
	// Daily nitrogen deposition distributed with fluxnet precipitation
	double dndep_fluxnet[Date::MAX_YEAR_LENGTH];
};

#endif // LPJ_GUESS_FLUXNET_H
