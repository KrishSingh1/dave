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
void calc_soil_water_potenial_dave(Soil &soil) {

    if (!ifdavehydraulics) {
        return;
    }

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

        // Also seems to be a mistake in this equation in Cosby 1984 as other
        // literature (e.g. Campbell 1974 and Clap & Hornberger 1978) use '-b',
        // no negative sign in Cosby 1984 but parameters seem to be the same.

        // Cosby 1984 uses Psi_s in cm, convert to MPa.
		soil.psi_layer[i] = Psi_s * pow((tot_wcont / sat_wcont), -b) / 100 * 1000 * 9.81 * 1e-6;
	}
}

double calc_weighted_psi(Soil &soil, Pft &pft) {

    double weighted_psi = 0;
    double rootfrac_sum = 0;

    for (int j = 0; j < NSOILLAYER; j++) {
        double psi_layer = soil.psi_layer[j];
        double rootfrac_layer = pft.rootdist[j];

        // Weight soil water potential by root fraction in each layer.
        weighted_psi = weighted_psi + psi_layer * rootfrac_layer;
        rootfrac_sum += rootfrac_layer;
    }

    // Divide by sum of root distribution, in case it doesn't sum to 1.
    return weighted_psi / rootfrac_sum;
}

void hydraulics(const Day &day, Individual &indiv) {
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
    double transpiration = gsw / (patmos / (vpd * 1000.0)) * 1000.0;

    if (indiv.stomata_open) {
        // Stomata are not fully closed.

        // For the purposes of this loop, we treat water potentials as positive
        // (rather than using abs everywhere).
        double psi1 = abs(indiv.weighted_psi);
        double psi2 = 0;

        double transpiration_psi_soil = 0;
        bool e_met = indiv.pft.get_transpiration_from_psi(psi1, transpiration_psi_soil);

        // Normalise kmax to xylem/leaf area?? because it reduces based on xlyem
        // loss, not leaf loss.
        double transpiration_full = transpiration_psi_soil + transpiration / indiv.kmax_drought_affected;
        e_met &= indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);

        if (!e_met) {
            // If hydraulic vulnerability curve doesn't allow E to be met, we
            // account for transpiration from capacitance in next step (closed
            // stomata loop). This happens often when ET is high and hydraulic
            // parameters (kmax, P50 and P88) are too small, check values.

            // Set water potential to stom_clos to close stomata.
            psi2 = -pft.stom_clos;
            // printf("hydraulic vulnerability curve doesn't allow E to be met\n");
        }

        // Calculate loss of conductance.
        indiv.k_loss[day.period] = 1 - exp(-pow(psi2 / pft.s1, pft.s2));

        indiv.lwp[day.period] = -psi2;
        indiv.E_not_met = !e_met;

        if (indiv.lwp[day.period] <= pft.stom_clos) {
            indiv.stomata_open = false;
            // printf("stomata have closed\n");
        }
    } else {
        // Stomata are fully closed.

        // Alternative for high water stress - consuming internal storage.

        // Average cross-sectional xylem area per unit ground area (m2/m2).
        // Calculated trunk area per m2, xylem area = fractional sapwood area.
        double Ax = PI * pow(indiv.diam / 2, 2) * indiv.densindiv * (indiv.cmass_sap / (indiv.cmass_sap + indiv.cmass_heart));

        // Water holding capacity approximated based on wood density: check this!
        double nx = 1 - (pft.wooddens * 2) / 1200;

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
            // printf("%d-%d hydraulic failure\n",date.get_calendar_year(),date.day);
            indiv.lwp[day.period] = pft.p88;
            indiv.hydraulic_failure = true;
        } else {
            indiv.lwp[day.period] =  pft.p88 * (1 - Vl / Vmax);
        }

        // Only account for this once (can only be set to 'true' in the loop for
        // open stomata).
        indiv.E_not_met = false;

        // Calculate loss of conductance.
        indiv.k_loss[day.period] = 1 - exp(-pow(-indiv.lwp[day.period] / pft.s1, pft.s2));
    }

    if (climate.prec > 0.0 && day.period == (date.subdaily - 1) && !indiv.stomata_open) {
        // If there is rain this day and stomata are closed, open stomata next
        // day.
        indiv.stomata_open = true;
        // printf("stomata have opened\n");
    }

    indiv.lwp_previous = indiv.lwp[day.period];
    indiv.transp_previous = transpiration;

    if (day.period > 0) {
        // W m-2; 0.5 * climate.insols
        double par_previous = climate.pars[day.period - 1];

        if (par_previous < 200 && climate.pars[day.period] > 200) {
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
            printf("\n %d-%d hydraulic failure %s\n", date.get_calendar_year(), date.day, (const char *)pft.name);
            double densindiv = indiv.densindiv;
            if (indiv.densindiv > KILL_THRESHOLD) {
                indiv.reduce_biomass(MORT_FRAC, 0);
                indiv.hydraulic_failure = false;
            } else {
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

        indiv.max_k_loss = 0;

        for (Day day; day.period != date.subdaily; day.next()) {
            indiv.max_k_loss = std::max(indiv.max_k_loss, indiv.k_loss[day.period]);
        }

        // Drought only has an impact if loss of conductance exceeds current
        // state.
        if (indiv.max_k_loss > (1 - indiv.kmax_drought_affected / pft.kmax)) {
            indiv.kmax_drought_affected = std::min(pft.kmax,indiv.kmax_drought_affected*(1-indiv.max_k_loss*indiv.kmax_drought_affected/pft.kmax)); // kmax decreases proportionately to loss but scaled by how affected it already is
            // kmax2[d] <- min(kmax_init, kmax2[d - 1] * (1 - drought[d] / 100 * kmax2[d - 1] / kmax_init))
        } else {
            // kmax recovery as sapwood regrows.

            // This line is designed to work with daily growth:
            // indiv.kmax_drought_affected = std::min(pft.kmax, (indiv.kmax_drought_affected * (1 - indiv.cmass_sap_turnover) + pft.kmax * indiv.cmass_sap_growth) * indiv.cmass_sap_yesterday / indiv.cmass_sap);

            // These lines are a messy placeholder trying to mimic the above
            // with annual growth:
            double day_turnover_sap = pft.turnover_sap / 365.0;
            double day_growth_sap = 0;
            if (indiv.cmass_sap_yesterday > 0.0) {
                // Growth defined as a unitless rate, as per turnover - check,
                // is cmass_sap_growth annual?
                day_growth_sap = indiv.cmass_sap_growth / indiv.cmass_sap / 365.0;
            } else {
                day_growth_sap = 0;
            }

            // indiv.cmass_sap in kg/m2
            double cmass_sap_today = std::max(0.0, indiv.cmass_sap + indiv.cmass_sap * day_growth_sap * (date.day + 1) - indiv.cmass_sap * day_turnover_sap * (date.day + 1));

            double cmass_sap_yesterday = std::max(0.0, indiv.cmass_sap + indiv.cmass_sap * day_growth_sap * date.day - indiv.cmass_sap * day_turnover_sap * date.day);
            if (cmass_sap_today > 0.0) {
                indiv.kmax_drought_affected = (indiv.kmax_drought_affected * (1 - day_turnover_sap) + pft.kmax * day_growth_sap) * cmass_sap_yesterday / cmass_sap_today;
                indiv.kmax_drought_affected = std::min(indiv.kmax_drought_affected, pft.kmax);
            } else {
                // If sapwood goes negative, stop changing kmax.
                indiv.kmax_drought_affected = indiv.kmax_drought_affected;
            }
            // R code
            // kmax2[d] <- (kmax2[d - 1] * (1 - turnover) + kmax_init * growth) * sapwood2[d - 1] / sapwood2[d]
        }

        vegetation.nextobj();
    }
}
