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
    
    bool getgridcell(Gridcell& gridcell);
    
    bool getclimate(Gridcell& gridcell);
	
	bool get_fluxnet_data_from_file();
    
private:
	std::vector<double> rain, tair, swrad;
	int first_fluxnet_year;
	
	// Daily precipitation during one year. Needed for the N-deposition.
	double drain[Date::MAX_YEAR_LENGTH];
	// Daily nitrogen deposition distributed with fluxnet precipitation
	double dndep_fluxnet[Date::MAX_YEAR_LENGTH];
};

#endif // LPJ_GUESS_FLUXNET_H
