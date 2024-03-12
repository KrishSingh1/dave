///////////////////////////////////////////////////////////////////////////////////////
/// \file hydraulics.cpp
/// \brief The plant hydraulics module
///
/// Water potentials calculated from soil water potential and gs
/// to inform drought mortality risk (P88)
///
/// \author Clare Stephens
/// $Date$
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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


#include "guess.h"
#include "hydraulics.h"
#include "canexch.h"

// Calculate soil water potential according to Campbell (1974)

void calc_soil_water_potenial_dave(Soil & soil) {

    for (int i=0; i<NSOILLAYER; i++){
	    double wcont = soil.get_layer_soil_water(i);
		double wilting_point_wc = soil.soiltype.wp[i]; // wilting point in mm
		double awc = soil.soiltype.awc[i];                     // awc in mm

		double tot_wcont = wcont * awc + wilting_point_wc; 

		double sat_wcont = soil.soiltype.wsats[i]; // wstats in mm

		double Psi_s = soil.soiltype.Psi_s;     // saturated pressure head in cm (Cosby 1984)
		double b = soil.soiltype.b;     // slope of soil water retention curve (Cosby 1984)

		soil.psi_layer[i] = Psi_s * pow((tot_wcont / sat_wcont),-b) / 100 * 1000 * 9.81 * 1e-6; // Cosby 1984 uses Psi_s in cm, convert to MPa. 
        // Also seems to be a mistake in this equation in Cosby 1984 as other literature (e.g. Campbell 1974 and Clap & Hornberger 1978) use '-b', no negative sign in Cosby 1984 but parameters seem to be the same.

	} 
}

double calc_weighted_psi(Soil & soil, Pft& pft){

    double weighted_psi;

    for (int j=0; j<NSOILLAYER; j++){
        double psi_layer = soil.psi_layer[j];
        double rootfrac_layer = pft.rootdist[j];
        weighted_psi = weighted_psi + psi_layer * rootfrac_layer; // weight soil water potential by root fraction in each layer
    }
    return weighted_psi;
}

void hydraulics(Day& day, const Climate& climate, Soil& soil, Individual& indiv) {

    if (!ifdavephotosyn) {
		// Only apply if using DAVE photosynthesis
		return;
	}

    Pft& pft = indiv.pft;
    Patchpft& ppft = indiv.patchpft();

    double patmos = climate.pressure[day.period];
    double vpd = climate.vpds[day.period];  

    double gsw = indiv.ps_subdaily.gsw[day.period];

    double transpiration = gsw / (patmos / (vpd * 1000.0)) * 1000.0; // mmol/m2/s 

    if (indiv.stomata_open){ // stomata are not fully closed

        double psi1 = abs(indiv.weighted_psi);  /// for the purposes of this loop, we treat water potentials as positive (rather than using abs everywhere)
        double E_guess = 0;   // mmol/m2/s
        double psi2 = 0;

        for (int x=1; x<10000000; x++){
            psi2 = psi1 + 0.0001;
            E_guess = E_guess + (pft.kmax*exp(-pow(((psi1+psi2)/2.0/pft.s1),pft.s2)))*(psi2-psi1);
            if (E_guess >= transpiration){ 
                break;
            } else {
                psi1 = psi2;
            }
            if (x == 10000000){
                printf("loop to estimate E failed/n");
            }
        }

        indiv.lwp[day.period] = -psi2;

        if (indiv.lwp[day.period] < pft.p10){
            indiv.stomata_open = false;
            printf("stomata have closed\n");
        }

    } else { // stomata are fully closed

        /// Alternative for high water stress - consuming internal storage

        // average cross-sectional xylem area per unit ground area (m2/m2)
        double Ax = PI*pow(indiv.diam/2,2)*indiv.densindiv*0.15;   // calculated trunk area per m2, 0.15 is assumed xylem:trunk area

        double Vmax = pft.nx*Ax*indiv.height; // see Paschalis et al. (2024), currently ignoring leaf water storage

        double Vl_init = Vmax*(1-indiv.lwp_previous/pft.p88); // note using P88 as approximation for LWP when Vl = 0. Equation from Paschalis et al. (2024)

        const double MOL2MMOL = 1000.0;
        const double M3_H2O_PER_KG = 0.001;

        double Vl = Vl_init - transpiration*(SECS_PER_DAY/date.subdaily) * MOL2G_H2O * KG_PER_G * M3_H2O_PER_KG / MOL2MMOL; // m3/m2

        if (Vl < 0){
            printf("plants are dead\n");
            indiv.lwp[day.period] = pft.p88;
        } else {
            indiv.lwp[day.period] =  pft.p88*(1-Vl/Vmax);
        }

        if (climate.precs[day.period] > 0.0){
            indiv.stomata_open = true;
            printf("stomata have opened\n");
        }


    }

    indiv.lwp_previous = indiv.lwp[day.period];

    if (day.period > 0){
        double par_previous   = climate.pars[day.period-1];    // W m-2; 0.5 * climate.insols
        if (par_previous < 200 && climate.pars[day.period] > 200){
            indiv.predawn_lwp = indiv.lwp[day.period]; // Note Zhou model uses actual (negative) pd lwp
        }
    }

    
    
}