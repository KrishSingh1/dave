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
const double MORT_FRAC = 0.1; // revisit how much dies based on level of stress over the year

// PAR threshold between predawn and postdawn (J/m2).
const double PAR_DAWN = 10;

double calc_weighted_psi_rootdist(Soil & soil, Pft& pft, Individual &indiv){
    
    double weighted_psi;

    weighted_psi = 0.0;
 
    for (int j=0; j<NSOILLAYER; j++){
        double psi_layer = soil.psi_layer[j];
        double rootfrac_layer = pft.rootdist[j];
        weighted_psi = weighted_psi + psi_layer * rootfrac_layer; // weight soil water potential by root fraction in each layer
    }
    return weighted_psi;
}

double calc_weighted_psi(Soil &soil, Pft &pft, Individual &indiv) { // De Kauwe et al (2015) recommends weighting by soil-to-root resistance

    double weighted_psi = 0;
    double Rs_inv_sum = 0;
    
    double rr = 0.5e-3; // fine root radius (m), from Table 2 in Williams et al (2001)
    double froot_density = 0.31e3; // specific fine root density, kg[biomass]/m3[root], from Table 3 in Bonan et al (2014)
    double root_resist =  25; // root resistivity in MPa*s*g/mmol, see Table 3 in Bonan et al (2014)

    for (int j = 0; j < NSOILLAYER; j++) {
        double thickness_lyr = soil.Dz[j+7] / 1000.0; // soil layer depth in m, note that we skip the 6 snow layers and layer of air in betwen (7 total)
        
        double rootc_layer = pft.rootdist[j] * indiv.cmass_root;    // kgC/m2
        double root_biomass = rootc_layer * 2.0; // kg/m2     !!UPDATE THIS - WHAT CONVERSION FACTOR IS SENSIBLE FOR ROOTS?????
        

        double lr = (1 / (PI * pow(rr, 2))) * root_biomass * (1 / froot_density); // length of fine roots in soil (m) 
        double rs = pow(1 / (PI * lr * thickness_lyr), 0.5); // mean distance between roots (m)
        
        double Ksoil = soil.K_actual[j] / 0.009807 * MMOL_PER_MOL / MOL2G_H2O * G_PER_KG * water_density; // convert from m/s to m2/s/MPa by dividing by head of pressure in MPa/m, then convert to mmol/m/s/MPa

        if (root_biomass > 0.0){
            indiv.Rs_r[j] = log(rs / rr) / (2.0 * PI * lr * thickness_lyr * Ksoil) + root_resist / (root_biomass);  // soil-root resistance + root resistance, see Williams et al. (2001). Note we don't multiply root biomass by D because it's already in kg/m2 (equation assumes kg/m3)
        } else {
            indiv.Rs_r[j] = pow(10.0,10.0); // conductivity is effectively zero if there are no roots
        }

        double psi_layer = soil.psi_layer[j];
        
        // Weight soil water potential by root fraction in each layer.
        weighted_psi = weighted_psi + psi_layer * (1.0 / indiv.Rs_r[j]);
        Rs_inv_sum += 1.0 / indiv.Rs_r[j];

    }

    // Divide by sum of 1/Rs to get overall weighted average.
    return weighted_psi / Rs_inv_sum;
}

// Calculate soil water potential according to Campbell (1974)
void calc_soil_water_potenial_dave(Soil &soil) {

    if (!ifzhou) {
        // ifzhou will always be true if dave hydraulics routine is enabled.
        return;
    }

    if (!ifdavehydrology) {
        fail("ifdavehydrology must be enabled to calculate soil water potential");
    }

    // Sand content, as percentage.
    double sand_perc = soil.soiltype.sand_frac * 100.0;

    // Clay content, as percentage.
    double clay_perc = soil.soiltype.clay_frac * 100.0;

    for (int i = 0; i < NSOILLAYER; i++) {

	    double wcont = soil.get_layer_soil_water(i);

        // Wilting point in mm.
		double wilting_point_wc = soil.soiltype.wp[i];

        // Fixed available water holding capacity of this layer in mm.
		double awc = soil.soiltype.awc[i];

		double tot_wcont = wcont * awc + wilting_point_wc;

        // Saturation point in mm.
		double sat_wcont = soil.soiltype.wsats[i];

        // Saturated pressure head in cm (Cosby 1984).
		double Psi_s = soil.soiltype.Psi_s;

        // Slope of soil water retention curve (Cosby 1984).
		double b = soil.soiltype.b;

        // Same approach used in De Kauwe et al (2020)
        
        // Seems to be a mistake in this equation in Cosby 1984 as other
        // literature (e.g. Campbell 1974 and Clap & Hornberger 1978) use '-b',
        // no negative sign in Cosby 1984 but parameters seem to be the same.

        // Cosby 1984 uses Psi_s in cm, convert to MPa.
		soil.psi_layer[i] = Psi_s * pow((tot_wcont / sat_wcont), -b) / 100 * 1000 * 9.81 * 1e-6;

        // double Ksoil_sat = soil.soiltype.Ksat; // saturated hydraulic conductivity in m/s
        // double k = pow(tot_wcont / sat_wcont, 2 * b + 3); // Equation 2 from Clapp & Hornberger 1978
        // soil.K_actual[i] = k * Ksoil_sat; // one way of converting Ksat to Kactual, but can also use method below

        double theta = tot_wcont / soil.Dz[i+7]; // convert from mm to m3/m3
        
        // hydraulic conductivity (m/s) based on equation 10 from Saxton et al (1986)
        soil.K_actual[i] = 2.778e-6 * exp(12.012 - 0.0755 * sand_perc + (-3.895 + 0.03671 *
                            sand_perc - 0.1103 * clay_perc + 8.7546e-4 * pow(clay_perc, 2.0)) * 1.0 / theta);
    }

    soil.soil_psi = mean(soil.psi_layer, NSOILLAYER);

    for (int i = 0; i < soil.patch.vegetation.nobj; i++) {
		Individual &indiv = soil.patch.vegetation[i];
		Pft &pft = indiv.pft;
        if (wateruptake == WR_ROOT_RESISTANCE) {
		    indiv.weighted_psi = calc_weighted_psi(soil, pft, indiv);
        } else {
            fail("Water uptake method different to root_resistance: Growth, soil, and photosynthesis modules need to be updated for compatibility with soil water potential calculated in hydraulics module");
        }
	}
}

void tree_hydraulics(const Day &day, Individual &indiv) {
    const Patch &patch = indiv.vegetation.patch;
    const Climate &climate = patch.get_climate();
    const Soil &soil = patch.soil;

    if (!ifdavephotosyn || !ifdavehydraulics) {
		// Only apply if using DAVE photosynthesis
		return;
	}

    const Pft &pft = indiv.pft;

    double patmos = climate.pressure[day.period];
    double vpd = climate.vpds[day.period];
    double gsw = indiv.ps_subdaily.gsw[day.period];

    // Transpiration in mmol/m2/s.
    double transpiration = gsw / (patmos / (vpd * 1000.0)) * 1000.0; /// shouldn't be proportional to stomatal conductance - look at how transpiration compares to flux data

    // Average cross-sectional xylem area per unit ground area (m2/m2).
    double Ax = PI * pow(indiv.diam / 2.0, 2.0) * indiv.densindiv * (indiv.cmass_sap / (indiv.cmass_sap + indiv.cmass_heart));

    if (indiv.stomata_open) {
        // Stomata are not fully closed.

        // For the purposes of this loop, we treat water potentials as positive
        // (rather than using abs everywhere).
        double psi1 = abs(indiv.weighted_psi);
        double psi2 = 0;

        double transpiration_psi_soil = 0;
        bool e_met = indiv.pft.get_transpiration_from_psi(psi1, transpiration_psi_soil);

        double transpiration_full = transpiration_psi_soil + transpiration / (indiv.kmax_drought_affected * 1000 * Ax); // convert from mol to mmol
        e_met &= indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);

        if (!e_met) {
            // If hydraulic vulnerability curve doesn't allow E to be met, we
            // account for transpiration from capacitance in next step (closed
            // stomata loop). This happens often when ET is high and hydraulic
            // parameters (kmax, P50 and P88) are too small, check values.

            // Set water potential to stom_clos to close stomata.
            psi2 = -pft.stom_clos;
            //printf("hydraulic vulnerability curve doesn't allow E to be met\n");
        }

        // Calculate loss of conductance.
        indiv.k_loss[day.period] = 1 - exp(-pow(psi2 / pft.s1, pft.s2)); // calculate loss of conductance from vulnerability curve

        indiv.lwp[day.period] = -psi2;
        indiv.E_not_met = !e_met;

        if (indiv.lwp[day.period] <= pft.stom_clos) {
            indiv.stomata_open = false;
            // printf("stomata have closed\n");
        }
    } else {
        // Stomata are fully closed.

        // Alternative for high water stress - consuming internal storage.

        // Water holding capacity approximated based on wood density:
        double nx = 1 - (pft.wooddens * 2) / 1000; // find out what paper had the wood density vs water holding capacity relationship

        // See Paschalis et al. (2024), currently ignoring leaf water storage.
        double Vmax = nx * Ax * indiv.height;

        // Note using P88 as approximation for LWP when Vl = 0. Equation from
        // Paschalis et al. (2024).
        double Vl_init = Vmax*(1-indiv.lwp_previous/pft.p88);

        const double MOL2MMOL = 1000.0;
        const double M3_H2O_PER_KG = 0.001;

        if (indiv.E_not_met) {
            // If the previous day's transpiration couldn't be met, we just
            // assume the plant starts at p12 but the previous day's
            // transpiration is taken from the capacitance store.
            Vl_init = Vl_init - indiv.transp_previous * (SECS_PER_DAY / date.subdaily) * MOL2G_H2O * KG_PER_G * M3_H2O_PER_KG / MOL2MMOL;
        }

        // m3/m2
        double Vl = Vl_init - transpiration * (SECS_PER_DAY / date.subdaily) * MOL2G_H2O * KG_PER_G * M3_H2O_PER_KG / MOL2MMOL;

        if (Vl < 0) {
            printf("%04d-%02d-%02d indiv %d (%s) hydraulic failure\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, indiv.id, (const char *)pft.name);
            indiv.lwp[day.period] = pft.p88;
            indiv.hydraulic_failure = true;
        } else {
            indiv.lwp[day.period] =  pft.p88 * (1 - Vl / Vmax);
        }

        // Only account for the extra ET once (can only be set to 'true' in the loop for
        // open stomata).
        indiv.E_not_met = false;

        // Calculate loss of conductance.
        indiv.k_loss[day.period] = 1 - exp(-pow(-indiv.lwp[day.period] / pft.s1, pft.s2));
    }

    if (climate.prec > 1.0 && day.period == (date.subdaily - 1) && !indiv.stomata_open) {
        // If there is >1mm rain this day and stomata are closed, open stomata next
        // day.
        indiv.stomata_open = true;
        // printf("stomata have opened\n");
    }

    indiv.lwp_previous = indiv.lwp[day.period];
    indiv.transp_previous = transpiration;

    if (day.period > 0) {
        // W m-2; 0.5 * climate.insols
        double par_previous = climate.pars[day.period - 1];

        if (par_previous < PAR_DAWN && climate.pars[day.period] > PAR_DAWN) {
            // Note Zhou model uses actual (negative) pd lwp.
            indiv.predawn_lwp = indiv.lwp[day.period];
        }
    }
}

// Hydraulic failure mortality applied at the daily timescale if allocation is
// daily, else at the annual timescale (can be different for different PFTs).
void hydraulic_mortality(Patch &patch) {

    Vegetation &vegetation = patch.vegetation;

    vegetation.firstobj();
    size_t nobj = vegetation.nobj;

    while(vegetation.isobj) {
        Individual &indiv = vegetation.getobj();
		Pft &pft = indiv.pft;

        if (!pft.dave_allocation() && !(date.islastday && date.islastmonth)) {
            vegetation.nextobj();
            continue;
        }

        if (indiv.hydraulic_failure) {
            double densindiv = indiv.densindiv;
            if (indiv.densindiv > KILL_THRESHOLD) {
                printf("%04d-%2d-%02d indiv %d (%s): biomass loss due to hydraulic failure\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, indiv.id, (const char *)pft.name);
                indiv.reduce_biomass(MORT_FRAC, 0);
                indiv.hydraulic_failure = false;
            } else {
                printf("%04d-%02d-%02d indiv %d (%s): mortality due to hydraulic failure\n", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1, indiv.id, (const char *)indiv.pft.name);
                indiv.kill();
                vegetation.killobj();
            }
            break;
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

        vegetation.nextobj();
    }
}

// Reduce kmax if drought damage has occurred today, otherwise allow kmax to
// recover with new growth / turnover.
void adjust_kmax(Patch &patch) {

    Vegetation &vegetation = patch.vegetation;

    vegetation.firstobj();
    size_t nobj = vegetation.nobj;

    while(vegetation.isobj) {
        Individual &indiv = vegetation.getobj();
		Pft &pft = indiv.pft;
        Patchpft &ppft = indiv.patchpft();

        if (pft.dave_allocation() || date.day == 0){ // if using annual allocation (not DAVE) then kmax loss is over whole year (not every day)
            indiv.max_k_loss = 0; // reset at beginning of year (annual allocation) or beginning of day (daily allocation)
        }

        for (Day day; day.period != date.subdaily; day.next()) {
            indiv.max_k_loss = std::max(indiv.max_k_loss, indiv.k_loss[day.period]); // take maximum k_loss for that day, will continue to maximise over all days of year if annual allocation
        }

        if (!pft.dave_allocation() && !(date.islastday && date.islastmonth)){ // if annual allocation, only move forward to xylem loss calcs at end of year
            vegetation.nextobj();
            continue;
        }

        if (indiv.max_k_loss > (1 - indiv.kmax_drought_affected/pft.kmax) && indiv.max_k_loss > 0.12){ // drought only has an impact if loss of conductance exceeds current state and also we've crossed P12
            double xylem_loss = indiv.max_k_loss * indiv.kmax_drought_affected / pft.kmax; // xylem capacity decreases proportionately to k_loss but scaled by how affected it already is

            //double closs_leaf = indiv.get_cmass_leaf() * xylem_loss; /// note dropping leaves causes problems due to fixed klatosa - consider adding back in with new allocation
            //double closs_root = indiv.cmass_root * xylem_loss;

            //double nloss_leaf = indiv.get_nmass_leaf() * xylem_loss; 
            //double nloss_root = indiv.nmass_root * xylem_loss;

            if (pft.is_grass()) {
                fail("Grasses not supported");
            }
            //indiv.cmass_leaf -= closs_leaf; // lose leaves and roots in proportion to lost xylem
            //ppft.cmass_litter_leaf += closs_leaf; // lost leaves and roots go to litter

            //indiv.cmass_root -= closs_root;
            //ppft.cmass_litter_root += closs_root;

            //indiv.nmass_leaf -= nloss_leaf;
            //ppft.nmass_litter_leaf += nloss_leaf;

            //indiv.nmass_root -= nloss_root;
            //ppft.nmass_litter_root += nloss_root;

            indiv.kmax_drought_affected = std::min(pft.kmax,indiv.kmax_drought_affected*(1-xylem_loss));
        }

        if (!pft.dave_allocation()){ // this loop allows kmax to recover as sapwood regrows
            double old_wood = std::max(indiv.cmass_sap - indiv.acalloc_sap, 0.0);
            double new_wood = std::max(std::min(indiv.acalloc_sap, indiv.cmass_sap),0.0);
            if (indiv.cmass_sap_yesterday == 0.0){ // there was no sapwood before, no drought legacy
                indiv.kmax_drought_affected = pft.kmax;
            } else { // drought-impacted wood is lost and fresh wood (with maximum kmax) regrows
                indiv.kmax_drought_affected = (indiv.kmax_drought_affected * old_wood + pft.kmax * new_wood) / (old_wood + new_wood);
            }
        } else {
            fail("Update sapwood turnover in accordance with new allocation");
        }

        assert(indiv.kmax_drought_affected >= 0.0);
        assert(less_than_or_equal(indiv.kmax_drought_affected,pft.kmax));


    //     MIGHT COME BACK TO THIS FOR DAILY ALLOCATION

    //     for (Day day; day.period != date.subdaily; day.next()) {
    //         indiv.max_k_loss = std::max(indiv.max_k_loss,indiv.k_loss[day.period]);
    //     }

    //     if (indiv.max_k_loss > (1 - indiv.kmax_drought_affected/pft.kmax)){ // drought only has an impact if loss of conductance exceeds current state
    //         indiv.kmax_drought_affected = std::min(pft.kmax,indiv.kmax_drought_affected*(1-indiv.max_k_loss*indiv.kmax_drought_affected/pft.kmax)); // kmax decreases proportionately to loss but scaled by how affected it already is   
    //    // kmax2[d] <- min(kmax_init, kmax2[d-1]*(1-drought[d]/100*kmax2[d-1]/kmax_init))
    
    
    //     } else { // kmax recovery as sapwood regrows
    //         //this line is designed to work with daily growth
    //         //  indiv.kmax_drought_affected = std::min(pft.kmax,(indiv.kmax_drought_affected*(1-indiv.cmass_sap_turnover) + pft.kmax*indiv.cmass_sap_growth)*indiv.cmass_sap_yesterday/indiv.cmass_sap);
    //      }

        vegetation.nextobj();
    }
}



void grass_hydraulics(const Day &day, Individual &indiv) {
    const Patch &patch = indiv.vegetation.patch;
    const Climate &climate = patch.get_climate();
    const Soil &soil = patch.soil;

    const Pft &pft = indiv.pft;
    Patchpft &ppft = indiv.patchpft();

    if (!ifdavephotosyn || !ifdavehydraulics) {
		// Only apply if using DAVE photosynthesis and hydraulics
		return;
	}

    if (!pft.is_grass()) {
        fail("This routine is only for grasses");
    }

    double patmos = climate.pressure[day.period];
    double vpd = climate.vpds[day.period];
    double gsw = indiv.ps_subdaily.gsw[day.period];

    // Transpiration in mmol/m2/s.
    double transpiration = gsw / (patmos / (vpd * 1000.0)) * 1000.0;

    // Stomata are not fully closed.

    // For the purposes of this loop, we treat water potentials as positive
    // (rather than using abs everywhere).
    double psi1 = abs(indiv.weighted_psi);
    double psi2 = 0;

    double transpiration_psi_soil = 0;
    bool e_met = indiv.pft.get_transpiration_from_psi(psi1, transpiration_psi_soil);

    // Normalise kmax to xylem/leaf area?? because it reduces based on xlyem
    // loss, not leaf loss.
    double transpiration_full = transpiration_psi_soil + transpiration / (pft.kmax * 1000 * indiv.lai);
    e_met &= indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);

    if (!e_met) {
        // If hydraulic vulnerability curve doesn't allow E to be met, we
        // assume stomata are fully closed which means grass will fully senesce

        // Set water potential to stom_clos to close stomata.
        psi2 = -pft.stom_clos;
        // printf("hydraulic vulnerability curve doesn't allow E to be met\n");
    }

    indiv.lwp[day.period] = -psi2;

    if (day.period > 0) {
        // W m-2; 0.5 * climate.insols
        double par_previous = climate.pars[day.period - 1];

        if (par_previous < PAR_DAWN && climate.pars[day.period] > PAR_DAWN) {
            // Note Zhou model uses actual (negative) pd lwp.
            indiv.predawn_lwp = indiv.lwp[day.period];
        }
    }

    double min_lwp_day = 0.0;

    // Senescence based on lowest water potential that day - runs at end of each day.
    if (day.period == date.subdaily - 1){

        for (Day day; day.period != date.subdaily; day.next()) {
            min_lwp_day = std::min(min_lwp_day, indiv.lwp[day.period]);
        }

        // We start to get senescence at turgor loss point, full senescense once stomata are fully closed.
        if (min_lwp_day <= pft.tlp) {

            double hyd_senesce = std::min(std::max(pft.tlp + 1.0/(pft.stom_clos-pft.tlp)*min_lwp_day, 0.0), 1.0); // fraction of lai lost to senescence increases linearly with lwp

            double closs_leaf = indiv.get_cmass_leaf() * hyd_senesce;
            double nloss_leaf = indiv.get_nmass_leaf() * hyd_senesce;

            if (ifdailygrass) {
                indiv.cmass_leafs[0] -= closs_leaf; //green leaf carbon removed
                indiv.nmass_leafs[0] -= nloss_leaf; //green leaf N removed
                indiv.cmass_leafs[1] += closs_leaf; //added to brown leaf carbon
                indiv.nmass_leafs[1] += nloss_leaf; //added to brown leaf N
            } else {
                fail("TODO: Should go to soil litter (discuss with Drew)");
                indiv.cmass_leaf -= closs_leaf;
                indiv.nmass_leaf -= nloss_leaf;
                ppft.cmass_litter_leaf += closs_leaf;
                ppft.nmass_litter_leaf += nloss_leaf;
            }
        }
    }
}
