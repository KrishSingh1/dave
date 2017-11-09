///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.cpp
/// \brief BLAZE fire simulation and combustion
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 17:03:10 +0100 (Tue, 24 Jan 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
// 
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "guess.h"
#include "blaze.h"
#include <fenv.h>

// Internal help function for splitting up nitrogen fire fluxes into components
// copied from vegdynam.cpp
void report_fire_flux_n(Patch& patch, double nflux_fire) {
	patch.fluxes.report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::NO_FIRE,  Fluxes::NO_FIRERATIO  * nflux_fire);
	patch.fluxes.report_flux(Fluxes::NO2_FIRE, Fluxes::NO2_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::N2_FIRE,  Fluxes::N2_FIRERATIO  * nflux_fire);
}

/// Compute area of Gridcell
double pixelsize(double latpos,double longsize,double latsize,int postype) {

	// taken from aslice.cpp in the utilities trunk

        // Returns area in square km of a pixel of a given size at a given point
        // on the world.  The formula applied is the surface area of a segment of
        // a hemisphere of radius r from the equator to a parallel (circular)
        // plane h vertical units towards the pole: S=2*pi*r*h

        // latpos    latitude position (see postype)
        // longsize  longitude range in degrees
        // latsize   latitude range in degrees
        // postype   declares which part of the pixel latpos
        //           refer to:
        //           0 = centre
        //           1 = NW corner
        //           2 = NE corner
        //           3 = SW corner
        //           4 = SE corner

	double pi,r,h1,h2,lattop,latbot,s;
      
	pi=3.1415926536;
	//r=6367.425;   // mean radius of the earth
	r=6371.2213;   // mean radius of the earth (revised LN 04/2015)
	
	lattop=latpos;
	if (postype==0) lattop=latpos+latsize*0.5;
	if (postype==3 || postype==4) lattop=latpos+latsize;
	if (lattop<0.0) lattop=-lattop+latsize;
	latbot=lattop-latsize;
	h1=r*sin(lattop*pi/180.0);
	h2=r*sin(latbot*pi/180.0);
	s=2.0*pi*r*(h1-h2);  //for this latitude band
	
	return s*longsize/360.0;  //for this pixel
}

void blaze_accounting_gridcell(Climate& climate) {
	
	/* Called by:  dailyaccounting_gridcell in driver.cpp 
	   Calls    :  available_fuel (local)
	               blaze_ignition (local)
	   routine to keep track of various met-related and fire specific 
	   parameters 
	*/
	//CLNdprintf("CLN blaze accounting GO  \n");

	const int average_span = 5; // span to average annual rainfall

	// some initialisations

	if (date.year == 0 && date.day == 0 && ! restart) {
		if ( vegmode == INDIVIDUAL ) fail("INDIVDUAL MODE not ready in BLAZE!");
		climate.avg_annual_rainf = 0.0;
		climate.cur_rainf        = 0.0;
		climate.dslr             = 0  ;
		climate.last_rainfall    = 0.0;
		climate.kbdi             = 0.0;
		climate.can_burn         = 0;
	}

	// Set Australian trees to be sprouters
	if (date.year == 0 && date.day == 0 ) {
		double lat = climate.gridcell.get_lat();
		double lon = climate.gridcell.get_lon();
		if ( lat < -10. && lon > 110. && lon < 158.) {
			climate.is_sprouter = 1;
		} else {
			climate.is_sprouter = 0;
		}
	}
             

	// Update running mean of average annual rainfall
	climate.cur_rainf += climate.prec;
	if (date.islastday && date.islastmonth) {
		double wght; // used to compute running average of ann rainfall
		if (date.year < average_span) {
			wght = date.year + 1;
		} else {
			wght = average_span;
		}
		climate.avg_annual_rainf = ((wght - 1.) * climate.avg_annual_rainf 
					    + climate.cur_rainf ) / wght;
		climate.cur_rainf   = 0.0;
	}

	//CLNdprintf("CLN blaze accounting 1  \n");

	// Keep track of Days-since-last-rainfall and accumulated last ainfall
	if (climate.prec > 0.01) {
		if (climate.dslr > 0) {
			climate.last_rainfall = climate.prec;
		} 
		else {
			climate.last_rainfall += climate.prec;
		}
		climate.dslr = 0;
	}
	else climate.dslr++;

	//CLNdprintf("CLN blaze accounting dslr %f  \n",climate.dslr);
	// Update the Keetch-Byram-Drought-Index 
	double v        = climate.u10   ; // Wind speed at 10m height [km/h] (for KBDI)
	double rh       = climate.relhum; // relative humidity [%] (for KBDI)         
	double t        = climate.tmax  ; // day's max temperature [deg C] (for KBDI) 

	v *= 3.6; // m/s -> km/h
	v = ( 214.7 * pow(  v + 10. ,-1.6968)  + 1. ) * v;

	double dkbdi; // change in kbdi due to rainfall history
	if (climate.dslr == 0) {
		if (climate.last_rainfall > 5.) {
			dkbdi = 5. - climate.last_rainfall;
		} 
		else {
			dkbdi = 0.0;
		}
	} 
	else {
		dkbdi = (( 800. - climate.kbdi) * (.968 * exp(.0486 * (t * 9./5. + 32.)) 
			 - 8.3) / 1000. / (1. + 10.88 * exp(-.0441 * 
			 climate.avg_annual_rainf/25.4)) * .254);
	}
	climate.kbdi = max(0.0,climate.kbdi + dkbdi);

	// ...and McArthur-Drought-Factor D ...
	double mcarthur_d = .191 * ( climate.kbdi + 104. ) * pow( climate.dslr + 1.,1.5 ) / 
		( 3.52 * pow( climate.dslr + 1. ,1.5 ) + climate.last_rainfall - 1. );
	mcarthur_d = max(0.0,min(10.0,mcarthur_d));
	
	// ... and finally: McArthur's Fire index
	double mcarthur_fire_index = 2. * exp( -.45 + .987 * log(mcarthur_d+.001) - 
				     .03456 * rh + .0338 * t + .0234 * v );
	mcarthur_fire_index = max(0.0,mcarthur_fire_index);

        // save old index
        double old_mfi = climate.mcarthur_fire_index;

        climate.mcarthur_fire_index = max(mcarthur_fire_index,old_mfi);

	// get burned area if sufficient biomassavailable fuel
	// loop over patches
	blaze_ignition(climate);

}		     

double available_fuel (Patch& patch,int flix)  {
			
	/* Called by:  get_firelineintensity (local)
	               blaze_account_gridcell (local)
 	   Calls    :  get_combustion_rates (local)
	   compute the amount of fuel that is readily available 
	   for burning 
	*/

	get_combustion_rates(patch,flix);

	//CLN Add trans-litter here!
	double trans_litter_leaf  = 0.;
	double trans_litter_sap   = 0.;
	double trans_litter_heart = 0.;
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& patchpft = patch.pft.getobj();
		trans_litter_leaf  += patchpft.litter_leaf;
		trans_litter_sap   += patchpft.litter_sap;
		trans_litter_heart += patchpft.litter_heart;
		patch.pft.nextobj();
	}


	// compute readily available fuel load for given FLI-index in PATCH
	double available_fuel = patch.litf2atm * (patch.soil.sompool[SURFSTRUCT].cmass + 
					   patch.soil.sompool[SURFMETA].cmass + 
					   trans_litter_leaf)
		+ patch.lfwd2atm * ( patch.soil.sompool[SURFFWD].cmass + trans_litter_sap )
		+ patch.lcwd2atm * ( patch.soil.sompool[SURFCWD].cmass + trans_litter_heart);

	Vegetation& vegetation=patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
		if (indiv.pft.lifeform == GRASS) {
			available_fuel += 0.5 * indiv.cmass_leaf ;
		}
		vegetation.nextobj();
	}
	return available_fuel ;
}

int get_fli_index(double fli, bool is_sprouter) {

	/* Called by:  get_firelineintensity (local)
	               get_combustion_rates (local)
	   Calls    :  -
	   get appropriate FLI category for look-up tables 
	*/

	// determine intensity category
	int flix; // fli - index
	if ( fli > 7000. ) {
		if ( is_sprouter ) {
			flix = 3;
		} else {
			flix = 4; 
		}
	}
	else if ( fli > 3000. ) {
		flix = 2;
	}
	else if ( fli > 750. ) {
		flix = 1;
	}
	else if ( fli > 0. ) {
		flix = 0;
	} 
	else {
		flix = -1;
	}
	return flix;
}
	

void get_firelineintensity(Patch& patch, Climate climate) {
	
	/* Called by:  blaze (local)
	   Calls    :  available_fuel (local)
	               get_fli_index (local)
	   compute potential fire-line intensity under given
	   met and fuel conditions
	*/

	const double heat_yield = 20.; // Energy contents of fuel [MJ/kg]

	double w;                      // Readily available fuel  [g/m2]
	double ros;                    // rate of spread          [m/s]
	double fli;                    // fire-line intensity     [W/m]

	int flix, i;                   // fire intensity category index, counter

	flix = 0;

	for ( i=0; i<4; i++ ) {

		// get available fuel and convert kg/m2 to g/m2
		w = available_fuel(patch,flix) * kg2g;
		// check whether a fire makes sense
		if ( w < min_fuel ) { 
			fli  =  -1. ;
			break;
		}
		// Compute Rate-of-spread [m/s]
		ros = 3.3333e-05 * climate.mcarthur_fire_index * w;
		
		// flame Height  [m]
		// original Z   = 13. * ROS + 0.24 * w - 2.
		//Z   = 46.8 * ROS + 0.024 * w - 2.;
		//Z   = MAX(0.,Z);
		
		// fire line intensity[W/m]
		fli = heat_yield * w * ros;
		
		flix = get_fli_index(fli, climate.is_sprouter);
		
		if (i >= flix ) break;
	}
	// check whether a fire makes sense
	patch.fli = max(patch.fli,fli);

}
/*
int gfed31_availability() {
	
	int cyear = date.get_calendar_year(); // current year
	// bitwise availability d,m,s,a => 1,2,4,8
	int av = 0;
	// annual data available
	if ( cyear >= 1997 && cyear <= 2012  ) av += 8;
	// seasonal availability 
	if ( cyear >= 1997 && cyear <= 2012  ) av += 4;


	if ( cyear < 1996 || (cyear == 1996 && date.month <= 6) ||
	     (cyear == 2012 && date.month >= 3) || cyear > 2012 ) {
		av = 0; // none
	}
	else if ( cyear >= 1996 && cyear <= 2012  )
		av = "a"

*/ 
bool burntime() { 

	/* Called by: blaze (local)
	              blaze_ignition (local) 
	   Calls    : sendmessage (plib)
	   function to determine appropriate time for burning w.r.t. to 
	   given ignition model/data-set, desired blaze_tstep and current date
	   returns TRUE if THIS timestep is a burn-timestep or FALSE else
	*/

	int cyear = date.get_calendar_year(); // current year

	if (blaze_tstep == DAILY) {
		// GFED3.1 only!!!
		if (ignition == GFED31 && (cyear < 2003 || cyear > 2011)) {
			sendmessage("Warning","No daily burned area available for this year");
		}
		return true;
	}
	else if (blaze_tstep == MONTHLY && date.islastday) {
		if ( ignition == GFED31 && (cyear <= 1996 || cyear >= 2012 )) {
			sendmessage("Warning","No monthly burned area available for this month");
		}
		return true;
	}
	else if (blaze_tstep == ANNUAL && date.islastday && date.islastmonth) {
		if ( ignition == GFED31 && (cyear <= 1996 || cyear >= 2012 )) {
			sendmessage("Warning","No GFED31 data available for this year");
		}
		return true;
	}
	else if (blaze_tstep == SEASONAL && date.islastday &&  
		 (date.month == 2 || date.month == 5 || date.month == 8 || date.month == 11)) {
		sendmessage("Warning","Please check end of seasons!!!!");
		if ( ignition == GFED31 && (cyear <= 1996 || cyear >= 2012 )) {
			sendmessage("Warning","No monthly GFED31 data available for this year");
		}
		return true;
	}
	else if (blaze_tstep == HYBRID ) {
		// this combines GFED 3.1 and SIMFIRE ignition 
		if (ignition == SIMGFED) {
			sendmessage("Error","HYBRID+SIMGFED not yet imlemented burntime()!");
			return false;
		}
	}
	else {
		return false;
	}
}

double surv_prob_boreal(double fli) {
	// following Dalziel et al. 2008
	double surv_prob_boreal = exp(-fli/500.);
	return surv_prob_boreal;
}

double surv_prob_temp_nl(double dbh, double fli, double mass_cwd) {
	// temperate Needleleaf forest
	// following Kobziar 2006
	// 
	double frac_cwd = 0.5;
	if ( fli > 7000. ) {
		frac_cwd = 0.8;
	}
	else if ( fli > 750. ) {
		frac_cwd = 0.75;
	}

	double cdbh = dbh * 100; // in cm
	double con1000 = frac_cwd * mass_cwd * 0.1 ; // in Mg/ha
	double p750, p_surv;
	if ( fli < 750. ) {
		p750   = 1. - (1./(1.+ exp(-(1.0337 + 0.000151*750. 
						- .221*cdbh + .0219*con1000))));
		p_surv = 1. - (fli/750. * (1. - p750) );
	}
	else {
		p_surv = 1. - (1./(1.+ exp(-(1.0337 + 0.000151*fli
						- .221*cdbh + .0219*con1000))));
	}
        //# WRONG allometry (NL)
	return p_surv;
}

double surv_prob_tropics(double dbh, double fli) {
	// van Nieuwstadt et al. 2005
	// fli as a scaling factor to mortaliaty from 0 to double
	double p_surv = 1.;
	double p_surv3000 = 1. - ( 0.82 - 0.035 * pow(dbh,0.7) );
	if ( fli > 7000. ) {
		double scal_fac = 1. - log((fli/7000.)) ;
		p_surv = scal_fac * p_surv3000;
	}
	else if ( fli > 3000. ) {
		p_surv =  p_surv3000;
	}
	else {
		p_surv = exp(fli/3000. * log(p_surv3000));
	}
	p_surv   = max(min(1.,p_surv), 0.001);
    	
	return p_surv;
}

double surv_prob_Savanna(double height, double fli) {
	// following  Bond 2008
	double intensity = fli / 1000. ;
	double p_surv = 1. - 1./(1. + exp(1.5*(height - 0.5 * intensity - 1. )));
    
	return p_surv;
}

double surv_prob_OzSavanna(double height, double fli) {
	
	/* Called by: combust (local)
	   Calls    : -
	   compute survival probability depending on height of 
	   individual or cohort and current fire-line intensity
	*/

	// inputs
	// height: tree/avg. cohort height
	// fli   : fire-line intensity from BLAZE [kW/m]
	
	// height of max survival probability [m]
	// taller trees are vulnerable due to age

	// Australian Savanna 
	double const max_prob_hgt = 8.5;

	// fire-line intensity [MW/m]
	double intensity = fli / 1000.; // Conversion to MW/m 

	// minimum height for trees to survive [m]
	double min_height = 3.7 * (1.-exp(-0.19 * intensity));

	// survival probability [fract.]
	double p_survival; 
	// Empirically generated functions by Vanessa Haverd
	// based on observations from G. Cook
	

	if (height > max_prob_hgt && height > min_height) {
		p_survival = ( -.0011 * intensity - .00002) * height
			+ .0075 * intensity + 1. ;
	}
	else if (height > min_height) {
		p_survival = ( .0178 * intensity + .0144) * height
			+ ( -.1174 * intensity + 0.9158 );
	}
	else {
		p_survival = 0.001;
	}
	p_survival = max(1.e-3,min(1.,p_survival));
	return p_survival;
}

double surv_prob_temp_bl(double dbh, double fli, bool res) {
	
	// This approach follows Hikcler et al. 2004, Using a generalized
	// vegetation model to simulate vegetation dynamics in NE USA
	double R;
	if ( res ) {
		R = 0.04;
	} else {
		R = 0.07;
	}

	// following Hickler et al. 2004
	double p_surv3000 = 0.95 - 1./(1.+ pow((dbh/R),1.5)) ;
	double surv_prob_temp_bl ;
	if ( fli > 7000. ) {
		surv_prob_temp_bl = 0.001;
	}
	else if ( fli > 3000 ) { 
		// Hickler 2004
		surv_prob_temp_bl = p_surv3000 * (1. - (fli-3000.)/ 4000. );
	}
	else {
		surv_prob_temp_bl = exp(fli/3000. * log(p_surv3000)); 
	}

	return surv_prob_temp_bl;
	
}

double survival_probability(Patch& patch, Individual& indiv, Climate& climate) {

	// Depending on biome and geolocation the appropriate survival_probabilities
	// will be selected

	Gridcell& gridcell = climate.gridcell;

	double height = indiv.height;
	double fli    = patch.fli;
	double lat    = gridcell.get_lat();
	double lon    = gridcell.get_lon();

	double survival_probability = 1.;

	double dbh   = pow(height * 100. / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3) / 100.;
	//double dbh   = pow(height* / 50., 1.5);

	if ( ignition == SIMFIRE ) {
		int biome = climate.simfire_biome;

		// Temperate Needleleaf
		if ( biome == 1) { 
			double mass_cwd = patch.soil.sompool[SURFCWD].cmass   ;   
			survival_probability = surv_prob_temp_nl(dbh, fli, mass_cwd);
			//CLN survival_probability = 1. - 0.25 * (1. -survival_probability); 
		}
		// Broadleaf and mixed
		else if ( biome == 2 || biome == 3 ) {
			// tropical 
			if  (lat > -30 && lat < 30 ) {
				// moist
				//CLN	dbh *= 0.5;
				survival_probability = surv_prob_tropics(dbh,fli);
			} else if (lat< -20. && lon > 110. && lon < 158.){
				// temperate Oz
				survival_probability = surv_prob_temp_bl(dbh, fli, 1);

			} else {
				// temperate 
				survival_probability = surv_prob_temp_bl(dbh, fli, 0);
				//CLN	survival_probability = 1. - 0.5 * (1. -survival_probability); 

			}
		}
		// Savanna, shrubland and sparsely vegetated
		else if ( biome == 4 || biome == 5 || biome == 7) {
			if ( lat < -10. && lon > 110.) {
				survival_probability = surv_prob_OzSavanna(height, fli);
			}
			else {
				survival_probability = surv_prob_Savanna(height, fli);
			}
		}
		// Tundra
		else if ( biome == 6 ) {
			survival_probability = surv_prob_boreal(fli);
		} 
		else {
			return 1.;
		}
	}

	// GFED and others
	else {
		survival_probability = surv_prob_OzSavanna(height, fli);
	}

	survival_probability = min(1. ,max(survival_probability,0.0001));

	//CLN dprintf("BIOME %d HGT %f DBH %f P_surv %f \n",climate.simfire_biome,height,dbh,survival_probability);
	return survival_probability;
}

void get_combustion_rates(Patch& patch, int flix) {

	/* Called by: combust (local)
	              available_fuel (local)
	   Calls    : get_fli_index (local)
	   compute the flux rates between live, litter pools and
	   atmosphere given current fire-line intensity
	*/

	// relative fluxes from wood to atmosphere and litter pools
	patch.wood2atm = (1.-fbranch-fbark) * turnoverfract[ 0][flix] +
		         fbranch            * turnoverfract[ 1][flix] +
		         fbark              * turnoverfract[ 2][flix];
	patch.wood2str = fbark              * turnoverfract[ 6][flix];
	patch.wood2fwd = fbranch            * turnoverfract[ 5][flix];
	patch.wood2cwd = (1.-fbranch-fbark) * turnoverfract[ 4][flix];
	
	// relative fluxes from leaf to atmosphere and litter pools
	patch.leaf2atm = turnoverfract[ 3][flix];
	patch.leaf2lit = turnoverfract[ 7][flix];

	// relative fluxes from litter pools to atmosphere
	patch.litf2atm = turnoverfract[11][flix];
	patch.lfwd2atm = turnoverfract[10][flix];
	patch.lcwd2atm = turnoverfract[ 9][flix];
	//CLN hier noch den transitional litter mit ein bauen!!
	return;
}
void combust(Patch& patch, Climate& climate) {

	/* Called by: blaze (local)
	   Calls    : survival_probability (local)
	              get_combustion_rates (local)
		      indiv.blaze_reduce_biomass (local)
		      randfrac (driver.cpp)
		      vegetation.<obj-functions> (guess.h)
		      allometry (growth.cpp)
		      negligible (guessmath.h)
	   Input ab: Area burnt [fract.]
	   The combustion part of the model. Here, all fire related fluxes
	   are computed and the changes applied to the affected pools.
	   This routine handles all current (i.e. April2015) available 
	   settings for vegetation model and uses the century som-model soil 
	   soil-pools.
	 */

	const double LIGCFRAC_leaf = 0.2;
	
	double ab  = climate.areaburnt;

	// Correction for previously burnt area
	double accf= 1. / (1. - climate.acc_areaburnt);
	
	// see if it burns at all
	if (! ( randfrac(patch.stand.seed) <= ab || vegmode == POPULATION ) )
		return ;


	// get relative fluxes between pools
	int flix = get_fli_index(patch.fli, climate.is_sprouter);

	// if fuel availability is too low return 
	if ( flix < 0 ) return;

	// CHECK FOR CLOSURE WITHIN BLAZE
	double c_prev = patch.ccont();
	double n_prev = patch.ncont();
	double c_fluxsum = 0.;
	double n_fluxsum = 0.;

	// adjustment factor for fluxes
	double fab = 1.0;
	if ( vegmode == POPULATION )
		fab = ab * accf;
       
	get_combustion_rates(patch,flix);

	dprintf("Here start burn =============== %i %i \n",date.year,flix);

	// compute fluxes FROM soil litter pools to atmosphere first!
	// since they are patch-specific only and the fluxes INTO 
	// soil litter will be added in loop over INDIVIDUALS below
	double cmtb2atm = fab * patch.litf2atm * patch.soil.sompool[SURFMETA].cmass   ;
	double cstr2atm = fab * patch.litf2atm * patch.soil.sompool[SURFSTRUCT].cmass ;
	double cfwd2atm = fab * patch.lfwd2atm * patch.soil.sompool[SURFFWD].cmass    ;
	double ccwd2atm = fab * patch.lcwd2atm * patch.soil.sompool[SURFCWD].cmass    ;   
	
	// nitrogen proportional to cmass flux
	double nmtb2atm = fab * patch.litf2atm * patch.soil.sompool[SURFMETA].nmass   ;
	double nstr2atm = fab * patch.litf2atm * patch.soil.sompool[SURFSTRUCT].nmass ;
	double nfwd2atm = fab * patch.lfwd2atm * patch.soil.sompool[SURFFWD].nmass    ;
	double ncwd2atm = fab * patch.lcwd2atm * patch.soil.sompool[SURFCWD].nmass    ;   
	
	// update soil-surface-litter pools
	// carbon
	patch.soil.sompool[SURFMETA].cmass   -= cmtb2atm ;
	patch.soil.sompool[SURFSTRUCT].cmass -= cstr2atm ;
	patch.soil.sompool[SURFFWD].cmass    -= cfwd2atm ;
	patch.soil.sompool[SURFCWD].cmass    -= ccwd2atm ;   
	
	// nitrogen 
	patch.soil.sompool[SURFMETA].nmass   -= nmtb2atm ;
	patch.soil.sompool[SURFSTRUCT].nmass -= nstr2atm ;
	patch.soil.sompool[SURFFWD].nmass    -= nfwd2atm ;
	patch.soil.sompool[SURFCWD].nmass    -= ncwd2atm ;   
 	
	// report C litter -> atm fluxes
	patch.fluxes.report_flux(Fluxes::FIREC, cmtb2atm + cstr2atm + cfwd2atm + ccwd2atm);

	// report N litter -> atm fluxes
	report_fire_flux_n(patch, nmtb2atm + nstr2atm + nfwd2atm + ncwd2atm );
       
	c_fluxsum += cmtb2atm + cstr2atm + cfwd2atm + ccwd2atm; 
	n_fluxsum += nmtb2atm + nstr2atm + nfwd2atm + ncwd2atm;

	//dprintf("cflux soil %f %f \n", c_fluxsum, c_prev - patch.ccont());

	Vegetation& vegetation=patch.vegetation;

	//kill trees, grass cohorts and population fractions 
	bool killed = false;
	double frac_survive;
	if ( vegmode == POPULATION ) {
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();
			// For this individual ...
			killed=false;
			// apply fire by burnt area
			indiv.blaze_reduce_biomass(patch, (1.-fab),c_fluxsum, n_fluxsum);

			// Remove this cohort completely if all individuals killed
			// (in individual mode: removes individual if killed)
			if (negligible(indiv.densindiv)) {
				vegetation.killobj();
				killed=true;
			}
			if (!killed) vegetation.nextobj(); // ... on to next individual
		}
		fab = ab * accf;
	}
	else {
		// vegmode == Individual or Cohort
		// taken from mortality_guess (vegdynam.cpp)
		// Impose fire in this patch with probability ab as given in the top
		// of this routine

		// Loop through individuals
		
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();
			if ( !indiv.alive ) {
				vegetation.nextobj();
				continue;
			}
			// For this individual ...
			killed=false;
			
			if (indiv.pft.lifeform==GRASS) {
				dprintf("GRASSO \n");
				// Reduce individual live biomass and freshly created litter
				indiv.reduce_biomass(.99,.99);
				
				// Update allometry
				allometry(indiv);
			}
			else {
				
				// TREE PFT
				if (ifstochmort) {
					// Impose stochastic mortality
					// Each individual in cohort dies with probability 'mort_fire'
					// Number of individuals represented by 'indiv'
					// (round up to be on the safe side)
					// In individual mode densindiv is 1!
					
					int nindiv=(int)(indiv.densindiv*patcharea+0.5);
					int nindiv_prev=nindiv;
					for (int i=0;i<nindiv_prev;i++) {
						if (randfrac(patch.stand.seed)  >=
						    survival_probability(patch,indiv, climate)) nindiv--;
					}
					
					if (nindiv_prev)
						frac_survive=(double)nindiv/(double)nindiv_prev;
					else
						frac_survive=0.0;
				}
				
				// Deterministic mortality (cohort mode only)
				
				else { 
					frac_survive=survival_probability(patch, indiv, climate);
				}

				// Reduce individual biomass on patch area basis
				// to account for loss of killed individuals
				dprintf("frac_survive id %i %f \n",indiv.id,frac_survive);
				indiv.blaze_reduce_biomass(patch,frac_survive,c_fluxsum, n_fluxsum);
				
				// Remove this cohort completely if all individuals killed
				// (in individual mode: removes individual if killed)
				if (negligible(indiv.densindiv)) {
					//dprintf("===================================[ \n");
					//dprintf("killed %s (%i) at patch %i y %i doy %i \n",(char*)indiv.pft.name,indiv.id,patch.id, date.get_calendar_year(), date.day);
			//dprintf("litter_leaf    %18.14f \n",litter_leaf  );
			//dprintf("litter_sap     %18.14f \n",litter_sap   );
			//dprintf("litter_heart   %18.14f \n",litter_heart );
			//dprintf("nmass_litter_leaf    %18.14f \n",indiv.pft.nmass_litter_leaf  );
			//dprintf("nmass_litter_sap     %18.14f \n",indiv.pft.nmass_litter_sap   );
			//dprintf("nmass_litter_heart   %18.14f \n",indiv.pft.nmass_litter_heart );
			//dprintf("totcmass    %18.14f \n",patch.ccont());
			//dprintf("totnmass    %18.14f \n",patch.ncont());
			dprintf("----------- KILL ----------- \n");
			//dprintf("anpp    %18.14f \n",indiv.anpp  );
					indiv.kill();
					vegetation.killobj();
					killed=true;

					//dprintf("anppa   %18.14f \n",indiv.anpp  );
					//dprintf("-------------------------- \n" );
			//dprintf("litter_leaf    %18.14f \n",indiv.pft.litter_leaf  );
			//dprintf("litter_sap     %18.14f \n",indiv.pft.litter_sap   );
			//dprintf("litter_heart   %18.14f \n",indiv.pft.litter_heart );
			//dprintf("nmass_litter_leaf    %18.14f \n",indiv.pft.nmass_litter_leaf  );
			//dprintf("nmass_litter_sap     %18.14f \n",indiv.pft.nmass_litter_sap   );
			//dprintf("nmass_litter_heart   %18.14f \n",indiv.pft.nmass_litter_heart );
			//dprintf("totcmass    %18.14f \n",patch.ccont());
			//dprintf("totnmass    %18.14f \n",patch.ncont());
			//dprintf("===================================] \n");
				} 
				else {
					// Update allometry 
					allometry(indiv);
				}

			}
			
			if (!killed) vegetation.nextobj(); // ... on to next individual
		}
	}
		
	// transitional LITTER removed by area burned (ab)
	// total fluxes away in patch
	double cmtb2atm_t = 0.; 
	double cstr2atm_t = 0.; 
	double cfwd2atm_t = 0.; 
	double ccwd2atm_t = 0.; 
	double nmtb2atm_t = 0.; 
	double nstr2atm_t = 0.; 
	double nfwd2atm_t = 0.;
	double ncwd2atm_t = 0.;
	
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& patchpft = patch.pft.getobj();
		
		double lton = lignin_to_n_ratio(patchpft.litter_leaf, patchpft.nmass_litter_leaf,
						LIGCFRAC_leaf, patchpft.pft.cton_leaf_avr);
		double fm_leaf = metabolic_litter_fraction(lton); 
		
		// carbon transitional litter fluxes
		double cmtb2atm = fab * patch.litf2atm * patchpft.litter_leaf * fm_leaf;
		double cstr2atm = fab * patch.litf2atm * patchpft.litter_leaf * (1.-fm_leaf);
		double cfwd2atm = fab * patch.lfwd2atm * patchpft.litter_sap;
		double ccwd2atm = fab * patch.lcwd2atm * patchpft.litter_heart;

		// nitrogen transitional litter fluxes
		double nmtb2atm = fab * patch.litf2atm * patchpft.nmass_litter_leaf * fm_leaf;
		double nstr2atm = fab * patch.litf2atm * patchpft.nmass_litter_leaf * (1.-fm_leaf);
		double nfwd2atm = fab * patch.lfwd2atm * patchpft.nmass_litter_sap;
		double ncwd2atm = fab * patch.lcwd2atm * patchpft.nmass_litter_heart;
		
		// update transitional rest-of-year litter pools
		// carbon
		patchpft.litter_leaf             -= (cmtb2atm + cstr2atm);
		patchpft.litter_sap              -= cfwd2atm;
		patchpft.litter_heart            -= ccwd2atm;

		// nitrogen
		// Are dropped at begining of year
		patchpft.nmass_litter_leaf       -= (nmtb2atm + nstr2atm);
		patchpft.nmass_litter_sap        -= nfwd2atm; 
		patchpft.nmass_litter_heart      -= ncwd2atm; 
		
		// calculate total fluxes away in patch
		// carbon
		cmtb2atm_t += cmtb2atm; 
		cstr2atm_t += cstr2atm; 
		cfwd2atm_t += cfwd2atm; 
		ccwd2atm_t += ccwd2atm; 
		// nitrogen
		nmtb2atm_t += nmtb2atm; 
		nstr2atm_t += nstr2atm; 
		nfwd2atm_t += nfwd2atm;
		ncwd2atm_t += ncwd2atm;
		
		patch.pft.nextobj();
	}
	
	
	// report C litter -> atm flux from transitional pools
	patch.fluxes.report_flux(Fluxes::FIREC, cmtb2atm_t + cstr2atm_t + cfwd2atm_t + ccwd2atm_t);
	// report N litter -> atm flux from transitional pools
	report_fire_flux_n(patch, nmtb2atm_t + nstr2atm_t + nfwd2atm_t + ncwd2atm_t );
	

	// CHECK CLOSURE
	c_fluxsum += cmtb2atm_t + cstr2atm_t + cfwd2atm_t + ccwd2atm_t; 
	n_fluxsum += nmtb2atm_t + nstr2atm_t + nfwd2atm_t + ncwd2atm_t;

	double dc_sum = c_prev - patch.ccont();
	double dn_sum = n_prev - patch.ncont();
	bool is_err = false;
	if ( abs(c_fluxsum - dc_sum) > 0.000000001 ) {
		dprintf("C-IMBALANCE in BLAZE: %18.14f \n",c_fluxsum - dc_sum);
		dprintf("dc_sum              : %18.14f \n",dc_sum);
		dprintf("c_fluxsum           : %18.14f \n",c_fluxsum);
		is_err = true;
	}	
	if ( abs(n_fluxsum - dn_sum) > 0.000000001 ) {
		dprintf("N-IMBALANCE in BLAZE: %18.14f \n",n_fluxsum - dn_sum);
		dprintf("dn_sum              : %18.14f \n",dn_sum);
		dprintf("n_fluxsum           : %18.14f \n",n_fluxsum);
		is_err = true;
	}	
	//if ( is_err ) fail("in blaze.cpp");
}  //combust

void Individual::blaze_reduce_biomass(Patch& patch, double frac_survive, 
				      double& c_fluxsum, double& n_fluxsum) {

	/* Called by: combust (local)
	   Calls    : lignin_to_n_ratio (somdynam.cpp)
	              metabolic_litter_fraction (somdynam.cpp)
	   Applies the fluxes computed in blaze on the class::Individual
	   level affecting the live pools, transitional 
	   litter pools and influx to CENTURY litter pools.
	   In INDIVIDUAL and COHORT mode the actual biomass killed in
	   combust is used to scale live fluxes while in POPULATION
	   mode the burnt area is taken.
	   Input frac_survive means fraction of surviving INDIVIDUAL/COHORT
	   in respective mode or will be burned area in case of POPULATION 
	   mode.
	*/

	bool now = false;
	//if ( patch.id == 16) now = true;

	double frac_killed = 1. - frac_survive;

	// FROM transfer_litter in somdynam.cpp
	// CLN SHOULD BE CENTRALLY SAVED!!!!! (somdynam.h???)
	// Leaf, root and wood litter lignin fractions
	// Leaf and root fractions: Comins & McMurtrie 1993; Friend et al 1997
	// Not sure of wood fraction
	const double LIGCFRAC_leaf = 0.2;
	const double LIGCFRAC_root = 0.16;
	//const double LIGCFRAC_wood = 0.3;

	if ( negligible(frac_killed) ) return;

	// local copies of live fluxes
	double wood2atm = patch.wood2atm;
	double wood2str = patch.wood2str;
	double wood2fwd = patch.wood2fwd;
	double wood2cwd = patch.wood2cwd;

	double leaf2atm = patch.leaf2atm;
	double leaf2lit = patch.leaf2lit;

	double fab = 1.0;

	if ( now ) dprintf("red bio 0    \n");
	if ( vegmode == INDIVIDUAL || vegmode == COHORT ) {
		double wtotw = patch.wood2atm + patch.wood2str + patch.wood2fwd + patch.wood2cwd;
		//CLN dprintf("w2atm/wtotw %f \n",patch.wood2atm/wtotw);
		// adjust relative fluxes from wood when stochastic killing has occured
		if ( wtotw > 0.0 ) {
			wood2atm = patch.wood2atm * frac_killed ;			
			//wood2str = patch.wood2str * frac_killed / wtotw ;
			//wood2fwd = patch.wood2fwd * frac_killed / wtotw ;
			//wood2cwd = patch.wood2cwd * frac_killed / wtotw ;
			wood2str = (1. - patch.wood2atm ) * fbark   * frac_killed ;
			wood2fwd = (1. - patch.wood2atm ) * fbranch * frac_killed ;
			wood2cwd = (1. - patch.wood2atm ) * (1.-fbark-fbranch) * frac_killed ;
		}
		else {
			wood2atm = frac_killed * .48;
			wood2str = frac_killed * .04;
			wood2fwd = frac_killed * .24; // CLN CHECK for ratio
			wood2cwd = frac_killed * .24; // CLN CHECK for ratio
		}
	
		// adjust relative fluxes from leaves
		double ltotw = patch.leaf2atm + patch.leaf2lit;
		if ( ltotw > 0.0 ) {
			//leaf2atm = patch.leaf2atm * frac_killed / ltotw ;
			//leaf2lit = patch.leaf2lit * frac_killed / ltotw ;
			leaf2atm = patch.leaf2atm * frac_killed ;
			leaf2lit = (1. - patch.leaf2atm) * frac_killed ;
		}			               
		else {			               
			leaf2atm = frac_killed * 1./3.;
			leaf2lit = frac_killed * 2./3.;
		}
		fab = 1.0;
	} 
	else if ( vegmode == POPULATION ) {
		fab = frac_killed;
		fail("The several x2y factors are zero nelow!!! indiv:blaze_reduce_biomass in blaze.cpp");
	}
	if ( now )dprintf("red bio 1   \n");
		

	// ===== compute mass/area fluxes ====
	Patchpft& ppft = patchpft();
	double lton = lignin_to_n_ratio(ppft.litter_leaf, ppft.nmass_litter_leaf, LIGCFRAC_leaf, 
					ppft.pft.cton_leaf_avr);


	// Metabolic litter fraction for litter
	double fm_leaf = metabolic_litter_fraction(lton);
	// double fm_leaf = 0.;

	// LEAVES
	double cleaf2atm = fab * leaf2atm * cmass_leaf;
	double cleaf2met = fab * leaf2lit * cmass_leaf * fm_leaf;
	double cleaf2str = fab * leaf2lit * cmass_leaf * (1. - fm_leaf);
	double nleaf2atm = fab * leaf2atm * nmass_leaf;
	double nleaf2met = fab * leaf2lit * nmass_leaf * fm_leaf;
	double nleaf2str = fab * leaf2lit * nmass_leaf * (1. - fm_leaf);

	// SAP-WOOD
	double csapw2atm = fab * wood2atm * cmass_sap ; 
	double csapw2str = fab * wood2str * cmass_sap ;
	double csapw2fwd = fab * (wood2fwd + wood2cwd) * cmass_sap ;
	double nsapw2atm = fab * wood2atm * nmass_sap ; 
	double nsapw2str = fab * wood2str * nmass_sap ;
	double nsapw2fwd = fab * (wood2fwd + wood2cwd) * nmass_sap ;

	// HEART-WOOD
	double chrtw2atm = fab * wood2atm * cmass_heart; 
	double chrtw2str = fab * wood2str * cmass_heart;
	double chrtw2cwd = fab * (wood2fwd + wood2cwd) * cmass_heart;
	double nhrtw2atm = fab * wood2atm * nmass_heart; 
	double nhrtw2str = fab * wood2str * nmass_heart;
	double nhrtw2cwd = fab * (wood2fwd + wood2cwd) * nmass_heart;

	// ROOT
	// assume that same percentage of root biomass is killed as for total 
	// above ground woody biomass
	double lossratio = 0.;
	if ( cmass_sap + cmass_heart > 0. ) {
		lossratio = (csapw2atm + csapw2str + csapw2fwd + 
			    chrtw2atm + chrtw2str + chrtw2cwd) / (cmass_sap + cmass_heart);
	}

	// Root litter lignin:N ratio
	lton = lignin_to_n_ratio(ppft.litter_root, ppft.nmass_litter_root, LIGCFRAC_root, 
					ppft.pft.cton_root_avr);

	// Metabolic litter fraction for root 
	double fm_root = metabolic_litter_fraction(lton);
	
	double croot2met = fm_root        * lossratio * cmass_root;
	double croot2str = (1. - fm_root) * lossratio * cmass_root;
	double nroot2met = fm_root        * lossratio * nmass_root;
	double nroot2str = (1. - fm_root) * lossratio * nmass_root;

	// ===== UPDATE POOLS ====

	// live carbon
	cmass_leaf      -= (cleaf2atm + cleaf2met + cleaf2str);
	cmass_sap       -= (csapw2atm + csapw2str + csapw2fwd);
	cmass_heart     -= (chrtw2atm + chrtw2str + chrtw2cwd); 
	cmass_root      -= (croot2met + croot2str); 
	
	if ( now ) dprintf("red bio 3    \n");
	// Deal with c-debt before asserting to litter & atm pool
	double saploss = csapw2atm + csapw2str + csapw2fwd;
	double hrtloss = chrtw2atm + chrtw2str + chrtw2cwd;

	if ( cmass_debt > 0. ) { 
		if ( cmass_debt <= saploss + hrtloss ) {
			if ( cmass_debt <= hrtloss ) {
				double scalefac = (hrtloss - cmass_debt) / hrtloss;
				chrtw2str  *= scalefac;
				chrtw2cwd  *= scalefac; 
				chrtw2atm  *= scalefac;
				cmass_debt  = 0.;
			}
			else {
				chrtw2str  = 0.;
				chrtw2cwd  = 0.; 
				chrtw2atm  = 0.;
				double scalefac = (saploss - (cmass_debt - hrtloss)) / saploss;
				csapw2str  *= scalefac;
				csapw2fwd  *= scalefac; 
				csapw2atm  *= scalefac;
				cmass_debt  = 0.;
			}
		}
		else {
			chrtw2str  = 0.;
			chrtw2cwd  = 0.; 
			chrtw2atm  = 0.;
			csapw2str  = 0.;
			csapw2fwd  = 0.; 
			csapw2atm  = 0.;
			cmass_debt -= (saploss + hrtloss);
		}
	}

	if ( now ) dprintf("red bio 4    \n");
	// live nitrogen
	nmass_leaf      -= (nleaf2atm + nleaf2met + nleaf2str);
	nmass_sap       -= (nsapw2atm + nsapw2str + nsapw2fwd);
	nmass_heart     -= (nhrtw2atm + nhrtw2str + nhrtw2cwd);
	nmass_root      -= (nroot2met + nroot2str); 
	
       	double d_nstore = (nstore_longterm + nstore_labile) * (1. - frac_survive);
	nstore_longterm *= frac_survive; 
	nstore_labile   *= frac_survive; 

	// deal with accumulated NPP
	double anpp_loss = anpp * (1. - frac_survive);
	double tot_loss  = cleaf2atm + csapw2atm + chrtw2atm + cleaf2met + 
		cleaf2str + csapw2str + chrtw2str + csapw2fwd + chrtw2cwd +
		croot2met + croot2str;
	double anpp2atm  =  0.;
	double anpp2met  =  0.;
	double anpp2str  =  0.;
	double anpp2fwd  =  0.;
	double anpp2cwd  =  0.;
	double anpp2smtb =  0.;
	double anpp2sstr =  0.;
	if ( tot_loss > 0. ) {
		anpp2atm  = anpp_loss * ( cleaf2atm + csapw2atm + chrtw2atm ) / tot_loss;
		anpp2met  = anpp_loss *                           cleaf2met   / tot_loss;
		anpp2str  = anpp_loss * ( cleaf2str + csapw2str + chrtw2str ) / tot_loss;
		anpp2fwd  = anpp_loss *                           csapw2fwd   / tot_loss;
		anpp2cwd  = anpp_loss *                           chrtw2cwd   / tot_loss;
		anpp2smtb = anpp_loss *                           croot2met   / tot_loss;
		anpp2sstr = anpp_loss *                           croot2str   / tot_loss;
		anpp     -= anpp_loss;
	}

	// report C live -> atm flux 
	patch.fluxes.report_flux(Fluxes::FIREC, cleaf2atm + csapw2atm + chrtw2atm + anpp2atm);

	// report N live -> atm flux 
	report_fire_flux_n(patch, nleaf2atm + nsapw2atm + nhrtw2atm + d_nstore);
 
	// CHECK CLOSURE
	c_fluxsum += cleaf2atm + csapw2atm + chrtw2atm + anpp2atm; 
	n_fluxsum += nleaf2atm + nsapw2atm + nhrtw2atm + d_nstore;

        /*patch.fluxes.report_flux(Fluxes::C_leaf2atm, cleaf2atm);
	patch.fluxes.report_flux(Fluxes::C_leaf2met, cleaf2met);
	patch.fluxes.report_flux(Fluxes::C_leaf2str, cleaf2str);
	patch.fluxes.report_flux(Fluxes::C_sapw2atm, csapw2atm);
	patch.fluxes.report_flux(Fluxes::C_sapw2str, csapw2str);
	patch.fluxes.report_flux(Fluxes::C_sapw2fwd, csapw2fwd);
	patch.fluxes.report_flux(Fluxes::C_hrtw2atm, chrtw2atm);
	patch.fluxes.report_flux(Fluxes::C_hrtw2str, chrtw2str);
	patch.fluxes.report_flux(Fluxes::C_hrtw2cwd, chrtw2cwd);
	patch.fluxes.report_flux(Fluxes::C_root2met, croot2met);
	patch.fluxes.report_flux(Fluxes::C_root2str, croot2str);*/


	// soil-surface-litter carbon
	patch.soil.sompool[SURFMETA].cmass   += anpp2met  + cleaf2met;
	patch.soil.sompool[SURFSTRUCT].cmass += anpp2str  + cleaf2str + csapw2str + chrtw2str;
	patch.soil.sompool[SURFFWD].cmass    += anpp2fwd  + csapw2fwd;
	patch.soil.sompool[SURFCWD].cmass    += anpp2cwd  + chrtw2cwd;   
	//deep soil litter carbon
	patch.soil.sompool[SOILMETA].cmass   += anpp2smtb + croot2met;
	patch.soil.sompool[SOILSTRUCT].cmass += anpp2sstr + croot2str;

	// soil-surface-litter nitrogen 
	patch.soil.sompool[SURFMETA].nmass   += nleaf2met;
	patch.soil.sompool[SURFSTRUCT].nmass += nleaf2str + nsapw2str + nhrtw2str;
	patch.soil.sompool[SURFFWD].nmass    += nsapw2fwd;
	patch.soil.sompool[SURFCWD].nmass    += nhrtw2cwd;   
	//deep soil litter nitrogen
	patch.soil.sompool[SOILMETA].nmass   += nroot2met;
	patch.soil.sompool[SOILSTRUCT].nmass += nroot2str;

	//dprintf("Burno hiero! \n");
	//CLN	if (vegmode == POPULATION && pft.lifeform != GRASS) {
	if (pft.lifeform != GRASS) {
		densindiv *= frac_survive;

		if ( negligible(densindiv) && cmass_debt > 0.) {
			//excess c_debt shifted from RA to NPP
			report_flux(Fluxes::NPP, cmass_debt);
			report_flux(Fluxes::RA, -cmass_debt);
			cmass_debt = 0.;
		}

		//		if ( negligible(densindiv) ) {
			//dprintf("===================================[ \n");
			//dprintf("killed with  \n");
			//dprintf("cmass_leaf  %18.14f \n",cmass_leaf  );
			//dprintf("cmass_sap   %18.14f \n",cmass_sap   );
			//dprintf("cmass_heart %18.14f \n",cmass_heart );
			//dprintf("cmass_root  %18.14f \n",cmass_root  );
			//dprintf("nmass_leaf  %18.14f \n",nmass_leaf  );
			//dprintf("nmass_sap   %18.14f \n",nmass_sap   );
			//dprintf("nmass_heart %18.14f \n",nmass_heart );
			//dprintf("nmass_root  %18.14f \n",nmass_root  );
			//dprintf("litter_leaf    %18.14f \n",ppft.litter_leaf  );
			//dprintf("litter_sap     %18.14f \n",ppft.litter_sap   );
			//dprintf("litter_heart   %18.14f \n",ppft.litter_heart );
			//dprintf("nmass_litter_leaf    %18.14f \n",ppft.nmass_litter_leaf  );
			//dprintf("nmass_litter_sap     %18.14f \n",ppft.nmass_litter_sap   );
			//dprintf("nmass_litter_heart   %18.14f \n",ppft.nmass_litter_heart );
			//dprintf("totcmass    %18.14f \n",ccont());
			//dprintf("totnmass    %18.14f \n",ncont());
		
//		}
	}

	//dprintf("SOMPOOLS BLAZE %i %i \n",patch.id, patch.stand.nobj);
	//dprintf("patch.soil.sompool[SURFMETA].nmass   %f \n",patch.soil.sompool[SURFMETA].nmass   );
	//dprintf("patch.soil.sompool[SURFSTRUCT].nmass %f \n",patch.soil.sompool[SURFSTRUCT].nmass );
	//dprintf("patch.soil.sompool[SURFFWD].nmass    %f \n",patch.soil.sompool[SURFFWD].nmass    );
	//dprintf("patch.soil.sompool[SURFCWD].nmass    %f \n",patch.soil.sompool[SURFCWD].nmass    );
	//dprintf("patch.soil.sompool[SOILMETA].nmass   %f \n",patch.soil.sompool[SOILMETA].nmass   );
	//dprintf("patch.soil.sompool[SOILSTRUCT].nmass %f \n",patch.soil.sompool[SOILSTRUCT].nmass );


}


void blaze_ignition(Climate& climate) {

	/* Called by: blaze_accounting_gridcell (local)
	   Calls    : simfire_ba (simfire.cpp)
	              gfed3_ba (IO???) CLN
	   provides burned area at given timestep by
	   calling appropriate IO-routines or 
	   model respectively
	*/

	Gridcell& gridcell = climate.gridcell;

	// initialisation of ba
	// ba will be zeroed in blaze after fire has occurred 
	if ( date.day == 0 && date.year == 0 ) {
		climate.areaburnt = 0.0;
	}

	// reset annual accumulative values
	if (date.day == 0) {
		climate.annual_areaburnt = 0.0;
		for (int i = 0; i < 12; i++) {
			climate.monthly_areaburnt[i] = 0.0;
		}
	}

	if ( ignition == SIMFIRE && burntime()) {
		climate.areaburnt += simfire_ba(climate, gridcell);
		//climate.annual_areaburnt += simfire_ba(climate, gridcell);
		//climate.monthly_areaburnt[date.month] += simfire_ba(climate, gridcell);
	}
	else if ( ignition == GFED31 ) {
		climate.areaburnt += gfed31_ba(gridcell);
		//CLN if (!negligible(climate.areaburnt)) dprintf("GFDED AB: %d %d %f \n",date.day ,date.year, gfed31_ba(gridcell));
	} 
	else if ( ignition == PRESCRIBED && burntime() ) {
		double tfac = 1. ;
		if ( blaze_tstep == DAILY ) {
			tfac = climate.monthly_fire_risk[date.month] /
				date.ndaymonth[date.month];
		}
		else if ( blaze_tstep == MONTHLY ) {
			tfac = climate.monthly_fire_risk[date.month];
		}
		else if ( blaze_tstep == ANNUAL ) {
			tfac = 1.;
		}
		else {
			fail ("PRESCRIBED Burning only available for daily, monthly, annual!");
		}
		climate.areaburnt = climate.prescribed_ba * tfac ;	
	}

} 

void blaze(Patch& patch, Climate& climate) {

	/* Called by: simulate_day (framework.cpp)
	   Calls    : get_firelineintensity (local)
	              burntime (local)
		      combust (local)
	   Does patch-wise accounting for FLI and activates
	   combustion when it's time
	*/

	// convert from km2 to ha
	const double kmsq2ha = 100.;

	// resolution for cell area
	// CLN Check for correct lat/lon res!!!
	double lat_res = 0.5;
	double lon_res = 0.5;

	// initialise fli
	if (date.day == 0 && date.year == 0)
		patch.fli = 0.0;

	// reset accumulated area_burnt to 0 on begining of year
	if (date.day == 0 )
		climate.acc_areaburnt = 0.0;
	
	// Accounting of max episodic fli
	get_firelineintensity(patch,climate);

	// start combustion at appropriate time-step
	if (burntime()) { 
	        // get relative fluxes between pools
	        int flix = get_fli_index(patch.fli, climate.is_sprouter);
		if ( flix >= 0 ) 
		        climate.can_burn += 1;

		if (!negligible(climate.areaburnt)) {
			combust(patch, climate);
		}
                //BLAZE-OUTPUT: climate.areaburnt, climate.mcarthur_fire_index, climate.areaburnt 

		// after burning reset accumulated variables
		patch.fli = 0.0;
		if ( patch.id == patch.stand.nobj-1 ) {
			// Now add BA to output if there was enough fuel...
			if ( climate.can_burn > 0 ) {
			        climate.annual_areaburnt              += climate.areaburnt;
			        climate.monthly_areaburnt[date.month] += climate.areaburnt;
			        climate.can_burn = 0;
			}
			//CLN			climate.max_nesterov        = 0.0;
			climate.areaburnt           = 0.0;
			//CLN climate.mcarthur_fire_index = 0.0;
		}
	}
}

