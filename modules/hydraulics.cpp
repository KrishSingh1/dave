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
#include "log_scope.h"
#include "canexch.h"
#include "somdynam.h"
#include "logging.h"
#include "config.h"

// The cohort density (/m2) above which cohorts will lose MORT_FRAC biomass
// upon hydraulic failure, and below which they will be killed outright.
const double KILL_THRESHOLD = 0.01;

// The fraction of biomass lost due to hydraulic failure.
// Revisit how much dies based on level of stress over the year (or day if using dave dailytree).
const double MORT_FRAC = 0.01;

// PAR threshold between predawn and postdawn (J/m2).
const double PAR_DAWN = 10;

// Alternative formulation that uses root distribution directly - recommended is 
// the root resistance method given in calc_weighted_psi
/// water density [kg m-3]
const double WATER_DENSITY = 1000.0;

double calc_weighted_psi_rootdist(Soil & soil, Pft& pft, Individual &indiv){

    double weighted_psi = 0.0;

    for (int j = 0; j < soil.nsoillayers; j++){
        double psi_layer = soil.soillayers[j].psi_mpa();
        double rootfrac_layer = indiv.rootdist[j];
        // Weight soil water potential by root fraction in each layer.
        weighted_psi = weighted_psi + psi_layer * rootfrac_layer;
    }
    return weighted_psi;
}

// De Kauwe et al (2015) recommends weighting by soil-to-root resistance.
double calc_weighted_psi(Soil &soil, Pft &pft, Individual &indiv) {

    double weighted_psi = 0;
    double Rs_inv_sum = 0;

    // Fine root radius (m), from Table 2 in Williams et al (2001).
    double rr = 0.5e-3;

    // Specific fine root density, kg[biomass]/m3[root], from Table 3 in Bonan
    // et al (2014).
    double froot_density = 0.31e3;

    // Root resistivity in MPa*s*g/mmol, see Table 3 in Bonan et al (2014).
    double root_resist =  25;

    for (int i = 0; i < soil.nsoillayers; i++) {
        // Soil layer depth in m, note that we skip the 6 snow layers and layer
        // of air in betwen (7 total).
        double thickness_lyr = soil.soillayers[i].thickness / MM_PER_M;

        // kgC/m2
        double rootc_layer = indiv.rootdist[i] * indiv.cmass_root;

        // Root biomass density, using 2.0 as conversion between carbon and biomass
        double root_biomass = rootc_layer * 2.0 / thickness_lyr; // kg/m3

        // Length density of fine roots in soil (m/m3).
        double lr = (1 / (PI * pow(rr, 2))) * root_biomass * (1 / froot_density);

        // Volumetric root density (m)
        double rs = pow(1 / (PI * lr), 0.5);

        // Convert from m/s to m2/s/MPa by dividing by head of pressure in
        // MPa/m, then convert to mmol/m/s/MPa.
        double Ksoil = soil.soillayers[i].k_actual / 0.009807 * MMOL_PER_MOL / MOL2G_H2O * G_PER_KG * WATER_DENSITY;

        if (root_biomass > 0.0 && soil.soillayers[i].wcont() > 0.0){
            // Soil-root resistance + root resistance, see Williams et al.
            // (2001). Note we don't multiply root biomass by D because it's
            // already in kg/m2 (equation assumes kg/m3). MPa * s * m2 / mmol    MPa * s * m2 / mmol
            indiv.Rs_r[i] = log(rs / rr) / (2.0 * PI * lr * Ksoil * thickness_lyr) + root_resist / (root_biomass * G_PER_KG * thickness_lyr);
        } else {
            // Resistance is infinite (conductivity exactly zero) if there are
            // no roots, or if the layer is at or below wilting point (no
            // plant-available water to assign uptake demand to). 
            indiv.Rs_r[i] = INFINITY;
        }

        double psi_layer = soil.soillayers[i].psi_mpa();

        // Weight soil water potential by root fraction in each layer.
        weighted_psi = weighted_psi + psi_layer * (1.0 / indiv.Rs_r[i]);
        Rs_inv_sum += 1.0 / indiv.Rs_r[i];
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

    double psi_sum = 0;
    for (int i = 0; i < soil.nsoillayers; i++) {
        Soillayer &layer = soil.soillayers[i];

        psi_sum += layer.psi_mpa();

        // Hydraulic conductivity, now using the same Campbell (1974)/
        // Clapp-Hornberger formula, Kw_max and b as Richards solver, instead of
        // previous Saxton et al. (1986) formula. The two were inconsistent:
        // Saxton used theta = wmass/thickness (fraction of total volume),
        // while the Richards solver (and this new version) uses
        // theta = wmass/wsats (fraction of saturation) - these differ
        // substantially since wsats < thickness (soil isn't 100% pore
        // space). Kw_max/b are in mm/s; divide by MM_PER_M to keep
        // k_actual's existing downstream unit contract (m/s) in
        // calc_weighted_psi(). Note: unlike the Richards solver, this
        // doesn't include ice impedance (no frozen-soil handling here).
        double theta = max(0.0, min(1.0, layer.wmass / layer.wsats));
        double k_actual_mm_s;
        if (theta >= 1.0) {
            k_actual_mm_s = layer.Kw_max;
        } else if (theta < 0.01) {
            k_actual_mm_s = 0.0;
        } else {
            k_actual_mm_s = layer.Kw_max * pow(theta, 2.0 * layer.b + 3.0);
        }
        layer.k_actual = k_actual_mm_s / MM_PER_M;
    }

    soil.soil_psi = psi_sum / soil.nsoillayers;

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
		// Only apply if using DAVE photosynthesis and hydraulics.
		return;
	}

    const Pft &pft = indiv.pft;

    indiv.predawn_lwp = indiv.weighted_psi; // Assume zero transpiration overnight

    double patmos = climate.pressure[day.period];
    double vpd = climate.vpds[day.period];
    double gsw = indiv.ps_subdaily.gsw[day.period];

    // Transpiration in mmol/m2/s.
    double transpiration = 0.0;

    if (penman_monteith){

        const std::vector<EnergyBalance> &ebal = indiv.energy_balance_layerwise[day.period];
	    for (size_t j = 0; j < ebal.size(); j++) {
		    transpiration += std::max(ebal[j].get_transp() * indiv.canopy[j].lai * 1000, 0.0);
	    }
    } else {
        transpiration = gsw / (patmos / (vpd * 1000.0)) * 1000.0;
    }

    // Average cross-sectional xylem area per unit ground area (m2/m2).
    double Ax = PI * pow(indiv.diam / 2.0, 2.0) * indiv.densindiv * (indiv.cmass_sap / (indiv.cmass_sap + indiv.cmass_heart));
    // Water holding capacity approximated based on wood density. Find out what
    // paper had the wood density vs water holding capacity relationship.
    double nx = 1 - (pft.wooddens * 2) / 1000;

    double Vmax = 0.0; // CALC BASED ON AREA UNDER CURVE BETWEEN STOM_CLOS AND P88
	HydraulicDiagnostic &diag = indiv.hydraulic_diagnostics[day.period];
	diag = HydraulicDiagnostic{};
	diag.valid = true;
	diag.stomata_open_start = indiv.stomata_open;
	diag.weighted_soil_psi = indiv.weighted_psi;
	diag.transpiration = transpiration;
	diag.conducting_area = Ax;

    if (indiv.stomata_open) {
        // Stomata are not fully closed.

        // For the purposes of this loop, we treat water potentials as positive
        // (rather than using abs everywhere).
        double psi1 = abs(indiv.weighted_psi);
        double psi2 = 0;

        // We calculate the area under the hydraulic curve up to the point of psi_soil. This
        // will be subtracted from the total are under the curve up to lwp to get at transpiration
        // between psi_soil and lwp. We work in a unitless space (i.e. per unit kmax) so the 
        // actual area under the curve would be transpiration_psi_soil x kmax
        double transpiration_psi_soil = 0;
        bool e_met = indiv.pft.get_transpiration_lower_bound(psi1, transpiration_psi_soil);
		diag.psi_start = psi1;
		diag.area_start = transpiration_psi_soil;
		diag.start_lookup_ok = e_met;

        double kmax_param;

        if (pft.dave_allocation()){
            // If we're using dave_allocation, we account for xylem loss through
            // sapwood turnover.
            kmax_param = pft.kmax;
        } else {
            // If we're using annual (trunk) allocation, we don't change the
            // sapwood but instead reduce the kmax parameter.
            kmax_param = indiv.kmax_drought_affected;
        }

        // Convert kmax from mol to mmol.
        double transpiration_full = transpiration_psi_soil + transpiration / (kmax_param * 1000 * Ax);
		diag.kmax = kmax_param;
		diag.area_requested = transpiration_full - transpiration_psi_soil;
		diag.area_full = transpiration_full;
		diag.end_lookup_ok = indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);
		diag.psi_end_lookup = psi2;
		e_met &= diag.end_lookup_ok;

        if (!e_met) {
            // If hydraulic vulnerability curve doesn't allow E to be met, parameters aren't
            // coordinating properly. Either a higher kmax or lower g1 may be required.
            double transpiration_check = 0.0;
	        indiv.pft.get_transpiration_lower_bound(-indiv.pft.p88,transpiration_check);
	        log_diagnostic("total area under curve = %f/n",transpiration_check * (kmax_param * 1000 * Ax));
            log_diagnostic("transpiration = %f/n",transpiration);

            log_error("indiv %i (%s): Hydraulic vulnerability curve doesn't "
                      "allow E to be met", indiv.id, (const char *)indiv.pft.name);
            
        }

        // Calculate loss of conductance.
        // Calculate loss of conductance from vulnerability curve (note this
        // will change sapwood, not k parameter, if using daily allocation)
        indiv.k_loss[day.period] = 1 - exp(-pow(psi2 / pft.s1, pft.s2));

        indiv.lwp[day.period] = -psi2; // - 0.01 * indiv.height; // height adjustment is to account for gravitational head

        if (indiv.lwp[day.period] <= pft.stom_clos) {
            indiv.stomata_open = false;
            // log_diagnostic("stomata have closed\n");
        }

    } else {
        // Stomata are fully closed.

        // Alternative for high water stress - consuming internal storage. Remaining 'area under the curve' between soil 
        // water potential at stom_clos and p88 is the volume available for transpiration before the plant dies. Roots are 
        // assumed to be fully disconnected from the soil, no further extraction.
        double psi1 = 0.0;

        if (day.period == 0){
            psi1 = abs(indiv.lwp[date.subdaily-1]);
        } else {
            psi1 = abs(indiv.lwp[day.period-1]);
        }

        double psi2 = 0;

        // We calculate the area under the hydraulic curve up to the point of the previous lwp. This
        // will be subtracted from the total are under the curve up to the new lwp to get at transpiration
        // between last timestep and this one.
        double transpiration_psi_prev = 0;
        bool e_met = indiv.pft.get_transpiration_lower_bound(psi1, transpiration_psi_prev);
		diag.psi_start = psi1;
		diag.area_start = transpiration_psi_prev;
		diag.start_lookup_ok = e_met;

        double kmax_param;

        if (pft.dave_allocation()){
            // If we're using dave_allocation, we account for xylem loss through
            // sapwood turnover.
            kmax_param = pft.kmax;
        } else {
            // If we're using annual (trunk) allocation, we don't change the
            // sapwood but instead reduce the kmax parameter.
            kmax_param = indiv.kmax_drought_affected;
        }

        // Convert from mol to mmol.
        double transpiration_full = transpiration_psi_prev + transpiration / (kmax_param * 1000 * Ax);
		diag.kmax = kmax_param;
		diag.area_requested = transpiration_full - transpiration_psi_prev;
		diag.area_full = transpiration_full;
		diag.end_lookup_ok = indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);
		diag.psi_end_lookup = psi2;
		e_met &= diag.end_lookup_ok;

        if (!e_met) {
            psi2 = -pft.p88;
            log_warning("indiv %i (%s): Hydraulic vulnerability curve doesn't "
                        "allow E to be met in water stressed module", indiv.id,
                        (const char *)indiv.pft.name);
        }

        indiv.k_loss[day.period] = 1 - exp(-pow(psi2 / pft.s1, pft.s2));

        indiv.lwp[day.period] = -psi2; // - 0.01 * indiv.height; // height adjustment is to account for gravitational head
   
        if (indiv.lwp[day.period] < pft.p88) {
            log_info("Hydraulic failure for indiv %d (%s)", indiv.id, (const char *)pft.name);
            indiv.lwp[day.period] = pft.p88;
            indiv.hydraulic_failure = true;
        }

        // Refill tree water stores if we have rainfall
        if (climate.precs[day.period] > 0.5 * 24 / date.subdaily) { // assume 0.5 mm every 1h is lost e.g. to interception
            double cmass_root_tot = 0;
            for (size_t i = 0 ; i < indiv.vegetation.nobj; i++) {
                cmass_root_tot += indiv.vegetation[i].cmass_root;
            }
            
            double prec_avail = std::max(climate.precs[day.period] - 0.1 * 24 / date.subdaily, 0.0); // rainfall that isn't intercepted / lost
            double frac_prec_this_indiv = indiv.cmass_root / cmass_root_tot; // this indiv gets an amount based on its root carbon

            // How much rain this individual gets to refill stem, convert from mm to mmol/m2/s.
            // 1 mm water = 1e-3 m * WATER_DENSITY kg/m3 = kg/m2, then kg -> g -> mol -> mmol.
            double seconds_per_timestep = (24.0 / date.subdaily) * 3600.0;
            double stem_refill = prec_avail * frac_prec_this_indiv * (1e-3 * WATER_DENSITY) * G_PER_KG / MOL2G_H2O * MMOL_PER_MOL / seconds_per_timestep;

            // Total volume under hydraulic curve up to current lwp, then subtract rainfall refill.
            // This will tell us the 'transpiration under the curve' value we need to read off the 
            // new lwp value with some stem refilling  
            double after_rain = transpiration_full - stem_refill / (kmax_param * 1000 * Ax);

            if (after_rain < 0.0){ // If there's more than enough rain to 'refill' the whole hydraulic curve
                indiv.stomata_open = true;
                indiv.lwp[day.period] = pft.stom_clos; 
            } else {
                e_met &= indiv.pft.get_psi_from_transpiration(after_rain, psi2);
                indiv.lwp[day.period] = -psi2; // - 0.01 * indiv.height; // height adjustment is to account for gravitational head
            }

            // Stomata reopen if volume refills to the point that lwp is equal
            // to stom_clos.

            if (indiv.lwp[day.period] >= pft.stom_clos) {
                indiv.stomata_open = true;
                indiv.lwp[day.period] = pft.stom_clos;
            }

        }

    }

	diag.lwp_end = indiv.lwp[day.period];
	diag.stomata_open_end = indiv.stomata_open;
    indiv.lwp_previous = indiv.lwp[day.period];

    // Calculate total transpiration in mm/timestep.
    indiv.transpiration[day.period] = transpiration / MMOL_PER_MOL * MOL2G_H2O * KG_PER_G * SECS_PER_DAY / date.subdaily;

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
		LogScope log_scope(indiv);

        if (indiv.hydraulic_failure) {
            double densindiv = indiv.densindiv;
            if (indiv.densindiv > KILL_THRESHOLD) {
                log_info("Indiv %d (%s): biomass loss due to hydraulic failure", indiv.id, (const char *)pft.name);
                indiv.reduce_biomass(MORT_FRAC, 0, 0, 0, 0, 0, false, indiv.closs_hydraulic, indiv.nloss_hydraulic);
                indiv.hydraulic_failure = false;
            } else {
                log_info("Indiv %d (%s): mortality due to hydraulic failure", indiv.id, (const char *)indiv.pft.name);
                indiv.kill(MORT_OTHER);
                vegetation.killobj();
            }
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

// Reduce hydraulic conductance if drought damage has occurred today (or this
// year if annual growth), otherwise allow conductance to recover with new
// growth / turnover. This is done via direct adjustment to the kmax parameter
// if running with annual growth, or via loss of xylem, leaves and roots if
// running with dave growth. Note that we still calculate kmax_drought_affected
// even if we aren't actually altering the kmax parameter in the hydraulics
// routine (ie when we run with DAVE growth) because it keeps tabs on drought
// damage over time and is used to scale the level of damage accounting for
// drought legacy as applicable
void adjust_kmax(Patch &patch) {

    Vegetation &vegetation = patch.vegetation;

    vegetation.firstobj();
    size_t nobj = vegetation.nobj;

    for (size_t i = 0; i < vegetation.nobj; i++) {
        Individual &indiv = vegetation[i];
		Pft &pft = indiv.pft;
        Patchpft &ppft = indiv.patchpft();

        // Only trees are affected by adjust_kmax(). Grasses use
        // grass_hydraulics().
        if (pft.is_grass()) {
            continue;
        }

        // if using annual allocation (not DAVE) then kmax loss is over whole
        // year (not every day).
        if (pft.dave_allocation() || date.day == 0){
            // Reset at beginning of year (annual allocation) or beginning of
            // day (daily allocation).
            indiv.max_k_loss = 0;
        }

        for (Day day; day.period != date.subdaily; day.next()) {
            // Take maximum k_loss for that day, will continue to maximise over
            // all days of year if annual allocation.
            indiv.max_k_loss = std::max(indiv.max_k_loss, indiv.k_loss[day.period]);
        }

        // If annual allocation, only move forward to conductance loss calcs at
        // end of year.
        if (!pft.dave_allocation() && !(date.islastday && date.islastmonth)) {
            continue;
        }

        // Drought only has an impact if stomata are fully closed, loss of
        // conductance exceeds current state and we've crossed P12.
        if (indiv.stomata_open == false && indiv.max_k_loss > (1 - indiv.kmax_drought_affected / pft.kmax) && indiv.max_k_loss > 0.12) {

            // Xylem capacity decreases proportionately to k_loss starting with
            // zero loss at p12 (we scale from here to avoid sudden 'threshold'
            // losses of 12%) but scaled by how affected it already is.
            indiv.xylem_loss_frac = (indiv.max_k_loss - 0.12) / (1.0 - 0.12) * indiv.kmax_drought_affected / pft.kmax;

            indiv.kmax_drought_affected = indiv.kmax_drought_affected * (1.0 - indiv.xylem_loss_frac);
            indiv.kmax_drought_affected = std::min(pft.kmax,indiv.kmax_drought_affected);

            if (pft.dave_allocation()) {

                // Note dropping leaves causes problems due to fixed klatosa -
                // consider adding back in with new allocation.
                double closs_sap = indiv.cmass_sap * indiv.xylem_loss_frac;
                double nloss_sap = indiv.nmass_sap * indiv.xylem_loss_frac;

                indiv.cmass_sap -= closs_sap;
                indiv.cmass_heart += closs_sap;
                indiv.nmass_sap -= nloss_sap;
                indiv.nmass_heart += nloss_sap;

                // Note dropping leaves causes problems due to fixed klatosa -
                // consider adding back in with new allocation.
                double closs_leaf = indiv.get_cmass_leaf() * indiv.xylem_loss_frac;
                double closs_root = indiv.cmass_root * indiv.xylem_loss_frac;

                double nloss_leaf = indiv.get_nmass_leaf() * indiv.xylem_loss_frac;
                double nloss_root = indiv.nmass_root * indiv.xylem_loss_frac;

                // Lose leaves and roots in proportion to lost xylem.
                // Lost leaves and roots go to litter.
                indiv.cmass_leaf -= closs_leaf;
                indiv.nmass_leaf -= nloss_leaf;
                leaf_to_litter(patch.soil, pft, closs_leaf, nloss_leaf);

                indiv.cmass_root -= closs_root;
                indiv.nmass_root -= nloss_root;
                root_to_litter(patch.soil, pft, closs_root, nloss_root);

                indiv.closs_embolism += closs_leaf + closs_root;
                indiv.nloss_embolism += nloss_leaf + nloss_root;
            }
        }

        double old_wood, new_wood;

        // This loop allows kmax to recover as sapwood regrows.
        if (!pft.dave_allocation()) {
            old_wood = std::max(indiv.cmass_sap - indiv.acalloc_sap, 0.0);
            new_wood = std::max(std::min(indiv.acalloc_sap, indiv.cmass_sap), 0.0);
        } else {
            old_wood = std::max(indiv.cmass_sap - indiv.cgrow_sap, 0.0);
            new_wood = std::max(std::min(indiv.cgrow_sap, indiv.cmass_sap), 0.0);
        }

        // There was no sapwood before, no drought legacy.
        if (old_wood == 0.0){
            indiv.kmax_drought_affected = pft.kmax;
        } else {
            // Drought-impacted wood is lost and fresh wood (with maximum kmax)
            // regrows.
            indiv.kmax_drought_affected = std::min((indiv.kmax_drought_affected * old_wood + pft.kmax * new_wood) / (old_wood + new_wood),pft.kmax);
        }

        assert(indiv.kmax_drought_affected >= 0.0);
        assert(less_than_or_equal(indiv.kmax_drought_affected,pft.kmax));
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
    bool e_met = indiv.pft.get_transpiration_lower_bound(psi1, transpiration_psi_soil);
	HydraulicDiagnostic &diag = indiv.hydraulic_diagnostics[day.period];
	diag = HydraulicDiagnostic{};
	diag.valid = true;
	diag.stomata_open_start = indiv.stomata_open;
	diag.weighted_soil_psi = indiv.weighted_psi;
	diag.psi_start = psi1;
	diag.area_start = transpiration_psi_soil;
	diag.start_lookup_ok = e_met;
	diag.transpiration = transpiration;
	diag.kmax = pft.kmax;
	diag.conducting_area = indiv.lai;

    // Normalise kmax to xylem/leaf area?? because it reduces based on xlyem
    // loss, not leaf loss.
    double transpiration_full = transpiration_psi_soil + transpiration / (pft.kmax * 1000 * indiv.lai);
	diag.area_requested = transpiration_full - transpiration_psi_soil;
	diag.area_full = transpiration_full;
	diag.end_lookup_ok = indiv.pft.get_psi_from_transpiration(transpiration_full, psi2);
	diag.psi_end_lookup = psi2;
	e_met &= diag.end_lookup_ok;

    if (!e_met) {
        // If hydraulic vulnerability curve doesn't allow E to be met, we
        // assume stomata are fully closed which means grass will fully senesce, but this shouldn't really happen

        // Set water potential to stom_clos to close stomata.
        psi2 = -pft.stom_clos;
        log_error("indiv %i (%s): Hydraulic vulnerability curve doesn't allow E"
                  " to be met", indiv.id, (const char *)indiv.pft.name);
    }

    indiv.lwp[day.period] = -psi2;
	diag.lwp_end = indiv.lwp[day.period];
	diag.stomata_open_end = indiv.stomata_open;

    double min_lwp_day = 0.0;

    // Senescence based on lowest water potential that day - runs at end of each day.
    if (day.period == date.subdaily - 1){

        for (Day day; day.period != date.subdaily; day.next()) {
            min_lwp_day = std::min(min_lwp_day, indiv.lwp[day.period]);
        }

        // We start to get senescence at turgor loss point, full senescense once stomata are fully closed.
        if (min_lwp_day <= pft.tlp) {
            // Fraction of LAI lost rises linearly from 0 at TLP to 1 at stomatal closure.
            double hyd_senesce = 0.0;
            double lwp_span = pft.tlp - pft.stom_clos;
            if (lwp_span > 0.0) {
                hyd_senesce = std::min(std::max((pft.tlp - min_lwp_day) / lwp_span, 0.0), 1.0);
            } else {
                // Fallback for invalid/equal thresholds.
                hyd_senesce = (min_lwp_day <= pft.stom_clos) ? 1.0 : 0.0;
            }

            double closs_leaf = indiv.get_cmass_leaf() * hyd_senesce;
            double nloss_leaf = indiv.get_nmass_leaf() * hyd_senesce;

            if (ifdailygrass) {
                indiv.cmass_leafs[0] -= closs_leaf; //green leaf carbon removed
                indiv.nmass_leafs[0] -= nloss_leaf; //green leaf N removed
                indiv.cmass_leafs[1] += closs_leaf; //added to brown leaf carbon
                indiv.nmass_leafs[1] += nloss_leaf; //added to brown leaf N
            } else {
                indiv.cmass_leaf -= closs_leaf;
                indiv.nmass_leaf -= nloss_leaf;
                leaf_to_litter(indiv.vegetation.patch.soil, pft, closs_leaf, nloss_leaf);

                indiv.closs_hydraulic += closs_leaf;
                indiv.nloss_hydraulic += nloss_leaf;
            }
        }
    }
}
