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

// The cohort density (/m2) above which cohorts will lose MORT_FRAC biomass
// upon hydraulic failure, and below which they will be killed outright.
const double KILL_THRESHOLD = 0.01;

// The fraction of biomass lost due to hydraulic failure.
const double MORT_FRAC = 0.1;

// Calculate soil water potential according to Campbell (1974)

void calc_soil_water_potenial_dave(Soil & soil) {

    if (!ifdavehydraulics) {
        return;
    }

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

    if (!ifdavephotosyn || !ifdavehydraulics) {
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
        double psi2 = 0;

        double transpiration_psi_soil = 0;
        bool e_met = indiv.pft.get_transpiration_from_psi(psi1, transpiration_psi_soil);

        double transpiration_full = transpiration_psi_soil + transpiration / indiv.kmax_drought_affected; // normalise kmax to xylem/leaf area?? because it reduces based on xlyem loss, not leaf loss.
        e_met &= indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);

        if (!e_met){ // if hydraulic vulnerability curve doesn't allow E to be met, we account for transpiration from capacitance in next step (closed stomata loop). This happens often when ET is high and hydraulic parameters (kmax, P50 and P88) are too small, check values.
                psi2 = -pft.stom_clos; /// set water potential to stom_clos to close stomata
        //        printf("hydraulic vulnerability curve doesn't allow E to be met\n"); 
        }

        // calculate loss of conductance
        indiv.k_loss[day.period] = 1-exp(-pow(((psi2)/pft.s1),pft.s2));

        indiv.lwp[day.period] = -psi2;
        indiv.E_not_met = !e_met;

        if (indiv.lwp[day.period] <= pft.stom_clos){
            indiv.stomata_open = false;
            // printf("stomata have closed\n");
        }

    } else { // stomata are fully closed

        /// Alternative for high water stress - consuming internal storage

        // average cross-sectional xylem area per unit ground area (m2/m2)
        double Ax = PI*pow(indiv.diam/2,2)*indiv.densindiv*(indiv.cmass_sap/(indiv.cmass_sap+indiv.cmass_heart));   // calculated trunk area per m2, xylem area = fractional sapwood area
        double nx = 1-(pft.wooddens*2)/1200;  /// water holding capacity approximated based on wood density - but check this!
        double Vmax = nx*Ax*indiv.height; // see Paschalis et al. (2024), currently ignoring leaf water storage

        double Vl_init = Vmax*(1-indiv.lwp_previous/pft.p88); // note using P88 as approximation for LWP when Vl = 0. Equation from Paschalis et al. (2024)

        const double MOL2MMOL = 1000.0;
        const double M3_H2O_PER_KG = 0.001;

        if (indiv.E_not_met){ // if the previous day's transpiration couldn't be met, we just assume the plant starts at p12 but the previous day's transpiration is taken from the capacitance store
            Vl_init = Vl_init - indiv.transp_previous*(SECS_PER_DAY/date.subdaily) * MOL2G_H2O * KG_PER_G * M3_H2O_PER_KG / MOL2MMOL;
        }

        double Vl = Vl_init - transpiration*(SECS_PER_DAY/date.subdaily) * MOL2G_H2O * KG_PER_G * M3_H2O_PER_KG / MOL2MMOL; // m3/m2

        if (Vl < 0){
            // printf("%d-%d hydraulic failure\n",date.get_calendar_year(),date.day);
            indiv.lwp[day.period] = pft.p88;
        } else {
            indiv.lwp[day.period] =  pft.p88*(1-Vl/Vmax);
        }

        indiv.E_not_met = false; // only account for this once (can only be set to 'true' in the loop for open stomata)

        // calculate loss of conductance
        indiv.k_loss[day.period] = 1-exp(-pow(((-indiv.lwp[day.period])/pft.s1),pft.s2));

    }
    
    if (climate.prec > 0.0 & day.period == date.subdaily-1 && !indiv.stomata_open){ // if there is rain this day and stomata are closed, open stomata next day
        indiv.stomata_open = true;
        // printf("stomata have opened\n");
    }
    
    indiv.lwp_previous = indiv.lwp[day.period];
    indiv.transp_previous = transpiration;

    if (day.period > 0){
        double par_previous   = climate.pars[day.period-1];    // W m-2; 0.5 * climate.insols
        if (par_previous < 200 && climate.pars[day.period] > 200){
            indiv.predawn_lwp = indiv.lwp[day.period]; // Note Zhou model uses actual (negative) pd lwp
        }
    }

   
}

/*
Reduce the biomass of the given individual by MORT_FRAC fraction.

This is complicated by the fact that individuals using trunk allocation cannot
have their removed biomass added to "actual" (e.g. litter) pools during the
first year of their "life", as they are still considered to be "on probation" at
this point in time, and their carbon pools are not yet considered to be part of
the system.
*/
void reduce_biomass(Individual &indiv) {
    if (indiv.pft.dave_allocation() || indiv.alive) {
        // "Normal" allocation - either daily (dave) allocation, or this indiv
        // has been alive for >1 year and is considered alive. In this case, we
        // can just call the standard reduce_biomass() function, which will move
        // biomass to the appropriate litter pools.
        indiv.reduce_biomass(MORT_FRAC, 0);
        return;
    }

    // trunk allocation + non-alive indiv. This is where things get hairy!

    // This is really not ideal, but we can't have biomass entering the live
    // biomass pools because a corresponding input flux has not yet been
    // reported for this indiv (nor will such a flux be reported until it
    // survives its first year of "life").
    double remaining = 1.0 - MORT_FRAC;

    indiv.cmass_leaf *= remaining;
    indiv.cmass_root *= remaining;
    indiv.cmass_sap *= remaining;

    // Debt and heartwood will be zero because allocation() has not yet been
    // called for this indiv.
    assert(is_zero(indiv.cmass_heart));
    assert(is_zero(indiv.cmass_debt));
    // indiv.cmass_heart *= remaining;
    // indiv.cmass_debt *= remaining;

    // Nitrogen is more complicated because all indivs are given daily N uptake.
    // In this case, N has already been moved from the soil to the indiv. We now
    // need to move the N somewhere, but it's not clear where it should go. If
    // the indiv dies while on probation, the N should probably go back into the
    // soil (not the litter pools), because the indiv is considered as never
    // having been alive in the first place. On the other hand, if the indiv
    // does survive its juvenile year, an argument could be made that any N lost
    // during that time via this hydraulic failure mechanism should have been
    // sent to litter (as that's what happens after the first year).
    //
    // For now, I'm going to simply move the N back into the soil pool, and note
    // that we may need to revisit this.

    double nloss = MORT_FRAC * (indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_sap);
    indiv.vegetation.patch.soil.nmass_inc(nloss);

    indiv.nmass_leaf *= remaining;
    indiv.nmass_root *= remaining;
    indiv.nmass_sap *= remaining;

    // Heartwood will be zero in first year of life, as growth() will not yet
    // have been called.
    assert(is_zero(indiv.nmass_heart));
    // indiv.nmass_heart *= remaining;

    // Reduce density.
    indiv.densindiv *= remaining;
}

void kill_indiv(Individual &indiv) {
    if (!indiv.pft.dave_allocation()) {
        // Trunk allocation and we are midway through the year.
        // (We are always midway through the year because this
        // is called before growth.)

        // NPP, GPP, and RA fluxes are reported daily, but the
        // carbon doesn't actually enter/leave the system until
        // the allocation occurs at the end of the year. If we
        // call kill() now and then go about our lives, we will
        // not have lost the C we claimed to have respired, and
        // similarly will not have gained the C we have claimed
        // to have assimilated, which will result in the mass
        // balance checks failing at the end of the year.
        indiv.report_flux(Fluxes::NPP, -indiv.anpp);
        indiv.report_flux(Fluxes::GPP, -indiv.agpp);
        indiv.report_flux(Fluxes::RA, -indiv.aresp);
    }

    indiv.kill();
    indiv.vegetation.killobj();
}

void hydraulic_failure(Patch& patch) {

    Vegetation& vegetation = patch.vegetation;

    vegetation.firstobj();
    size_t nobj = vegetation.nobj;

    while(vegetation.isobj) {
        Individual& indiv = vegetation.getobj();
		Pft& pft = indiv.pft;

        indiv.max_k_loss = 0;

        for (Day day; day.period != date.subdaily; day.next()) {
            if (indiv.lwp[day.period] <= pft.p88){
                double densindiv = indiv.densindiv;
                if (indiv.densindiv > KILL_THRESHOLD) {
                    reduce_biomass(indiv);
                } else {
                    kill_indiv(indiv);
                }
                break;
            }
            indiv.max_k_loss = std::max(indiv.max_k_loss,indiv.k_loss[day.period]);
        }

        if (nobj != vegetation.nobj) {
            // This indiv was killed due to hydraulic failure. We cannot proceed
            // with this iteration (how could we?), and vegetation.getobj() will
            // now return the indiv which was previously after the indiv which
            // was killed, so we don't want to call indiv.nextobj(). We also
            // need to set nobj to vegetation.nobj so that this check still
            // functions correctly for the remaining indivs in vegetation.
            nobj = vegetation.nobj;
            continue;
        }

        if (indiv.max_k_loss > (1 - indiv.kmax_drought_affected/pft.kmax)){ // drought only has an impact if loss of conductance exceeds current state
            indiv.kmax_drought_affected = std::min(pft.kmax,indiv.kmax_drought_affected*(1-indiv.max_k_loss*indiv.kmax_drought_affected/pft.kmax)); // kmax decreases proportionately to loss but scaled by how affected it already is   
       // kmax2[d] <- min(kmax_init, kmax2[d-1]*(1-drought[d]/100*kmax2[d-1]/kmax_init))
    
    
        } else {
            //this line is designed to work with daily growth
            //  indiv.kmax_drought_affected = std::min(pft.kmax,(indiv.kmax_drought_affected*(1-indiv.cmass_sap_turnover) + pft.kmax*indiv.cmass_sap_growth)*indiv.cmass_sap_yesterday/indiv.cmass_sap);
            //these lines are a messy placeholder trying to mimic the above with annual growth
            double day_turnover_sap = pft.turnover_sap/365.0;
            double day_growth_sap = 0;
            if (indiv.cmass_sap_yesterday > 0.0){
                day_growth_sap = indiv.cmass_sap_growth/indiv.cmass_sap/365.0; // growth defined as a unitless rate, as per turnover - check, is cmass_sap_growth annual?
            } else {
                day_growth_sap = 0;
            }
            double cmass_sap_today = std::max(0.0,indiv.cmass_sap+indiv.cmass_sap*day_growth_sap*(date.day+1)-indiv.cmass_sap*day_turnover_sap*(date.day+1)); // indiv.cmass_sap in kg/m2
            double cmass_sap_yesterday = std::max(0.0,indiv.cmass_sap+indiv.cmass_sap*day_growth_sap*date.day-indiv.cmass_sap*day_turnover_sap*date.day);
            if (cmass_sap_today > 0.0){
                indiv.kmax_drought_affected = (indiv.kmax_drought_affected*(1-day_turnover_sap) + pft.kmax*day_growth_sap)*cmass_sap_yesterday/cmass_sap_today;
                indiv.kmax_drought_affected = std::min(indiv.kmax_drought_affected,pft.kmax);
            } else {
                indiv.kmax_drought_affected = indiv.kmax_drought_affected; // if sapwood goes negative, stop changing kmax
            }
            // kmax2[d] <- (kmax2[d-1]*(1-turnover) + kmax_init*growth)*sapwood2[d-1]/sapwood2[d] // R code
            if (date.year == 100){
                printf("");
            }

        }

        vegetation.nextobj();
        
    }

}
