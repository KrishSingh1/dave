///////////////////////////////////////////////////////////////////////////////////////
/// \file soil.cpp
/// \brief LPJ-GUESS Soil class and related definitions
///
/// Implementation of member functions of class Soil.
/// The class Soil and its member functions and variables are declared in guess.h
///
/// \author Paul Miller
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "soil.h"
#include "parameters.h"
#include "guessmath.h"
#include "lookup.h"

// Anonymous namespace for variables with file scope
namespace {

/// lookup table for natural exponential for thermal calculations in soil layers
NatExp_LUT exp_thermal_lut(Kpeat_log, Kice_log, 0.001, 1.0);

/// lookup table for natural exponential for hydrological calculations in soil layers
/// min and max are estimated from a minimum value of theta of 0.01 and min and max on b
NatExp_LUT exp_hydro_lut(-187.0, 87.0, 0.01, 1.0);

/// lookup table for natural logarith for hydrological calculations in soil layers
NatLog_LUT log_hydro_lut(0.001, 1.0, 0.0001, 1.0);

}

////////////////////////////////////////////////////////////////////////////////
// SNOWLAYER

/// Clears all snow layer members
/** Clear all snow layer members
 */
void Snowlayer::clear() {
	imass = 0.0;
	wmass = 0.0;
	dens = 0.0;
	thickness = 0.0;
	whc = 0.0;
	temp = 0.0;
	Kh = 0.0;
	Ch = 0.0;
}

/// Calculates heat capacity of snow
/** Calculates heat capacity of snow
 *  Fukusako, Eq. 2. (J m-3 K-1) ///Best et al. 2011, Eq. 22 [J m-3 K-1]
 */
double snowheatcap(double imass, double wmass, double thickness, double temp) {
	//return (Cp_ice * imass + Cp_water * wmass + Cp_air * max(0.0, (depth - wmass - imass))) / depth;

	// Snow density including water withing snow
	double tot_dens = (imass + wmass) / thickness * MM_PER_M;
	return (0.185 + 0.689 * (K2degC + temp) * 0.01) * tot_dens * J_PER_KJ; // J kg-1 K-1
}

/// Calculates thermal conductivity of snow
/** Calculates thermal conductivity of snow
 *  Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
 */
double snowcond(double imass, double wmass, double depth) {

	// Snow density including water withing snow
	double tot_dens = min(rho_ICE, (imass + wmass) / depth * MM_PER_M);

	// Thermal conductivity of snow
	return 2.22 * pow((tot_dens / rho_H2O), 1.88);
}

/// Updates snow layer temperature
/** Updates snow layer temperature due to movement of mass 
 * with different temperature
 */
void Snowlayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible(thickness * Ch + add_mass * add_Ch))
		dprintf("Year %d Day %d Error snow temp \n",date.year,date.day);

	// New snow layer but haven't updated snow properties yet 
	if (!negligible(thickness) && negligible(Ch)) {
		update_thermal_properties();
	}

	// No snow depth nor heat capacity
	if (negligible(thickness) && negligible(Ch)) {
		Ch = add_Ch;
	}

	// Add to existing snow layer with updated snow properties
	if (!negligible(thickness) && !negligible(Ch)) {
		temp = (temp * thickness * Ch + add_temp * add_mass * add_Ch) / (thickness * Ch + add_mass * add_Ch);
	}
	else { 
		// No existing snow layer, so added snow temperature is new snow layer temperature
		temp = add_temp;
	}
}

/// Update heat capacity and conductivity of snow layer
/** Update heat capacity and conductivity of snow layer
 */
void Snowlayer::update_thermal_properties() {

	// Heat capacity of snow - Fukusako, Eq. 2. (J m-3 K-1)
	Ch = snowheatcap(imass, wmass, thickness, temp);

	// Thermal conductivity of snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = snowcond(imass, wmass, thickness);
}

/// Update snow water holding capacity
/** Update snow water holding capacity
 * based on eqn 10 Wang et al. 2013
 */
void Snowlayer::update_whc() {

	// Reference density [kg m-3]
	const double dens_t = 400;	
	// Minimum snow water holding capacity
	const double rw_min = 0.03;
	// Maximum snow water holding capacity
	const double rw_max = 0.1;

	// Water holding capacity of snow layer [mm]
	whc = max(0.0, min(thickness - imass, imass * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - dens) / dens_t)));
}

/** Phase change of snow layer
 */
void Snowlayer::phase_change(bool toplayer) {

	// Initial water and ice content
	double wtotal_initial = wmass + imass;

	// Should ice thaw
	bool ifThaw = false;
	// Should water freeze
	bool ifFreeze = false;
	// Fraction of water which freezes
	double Ffreez = 0.0;
	// Fraction of ice which thaws
	double Fthaw = 0.0;
	// Fraction of water
	double Fwmass;
	// Fraction of ice
	double Fimass;
	// Energy per volume [J m-3]
	double energy;
	// Temperature at freezing
	double temp_old = 0.0;
	// Temperature difference (current - previous)
	double delta_T = temp - temp_old;
	// Maximum error
	const double EPS = 1.0e-10;
	// Addition thickness to top snow layer to accound for radiation energy
	double lyr_thick = toplayer ? thickness + top_thickness : thickness;

	// Water and ice fraction of full volume
	Fwmass = wmass / lyr_thick;
	Fimass = imass / lyr_thick;

	// Should water freeze
	if (temp < 0.0 && Fwmass > 0.0) {
		ifFreeze = true;
	}
	// Should ice thaw
	else if (temp > 0.0 && Fimass > 0.0) {
		ifThaw = true;
	}

	// No phase change, return zero
	if (!ifThaw && !ifFreeze) {
		return;
	}

	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy from full soil layer to freeze water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction water that would freeze [m3 m-3]
		Ffreez = energy / Lheat;

		if (Fwmass > Ffreez) {

			// Update water and ice content
			wmass = (Fwmass - Ffreez) * lyr_thick;
			imass = (Fimass + Ffreez) * lyr_thick;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			wmass = 0.0;
			imass = (Fimass + Fwmass) * lyr_thick;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fwmass * Lheat;

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update soil temperature
			temp = temp_old - energy_remain / Ch;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from full soil layer to thaw water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = energy / Lheat;

		if (Fimass > Fthaw) {

			// Update water and ice content
			wmass = (Fwmass + Fthaw) * lyr_thick;
			imass = (Fimass - Fthaw) * lyr_thick;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass = (Fwmass + Fimass) * lyr_thick;
			imass = 0.0;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fimass * Lheat;

			// Calculate remaining energy to increase the temperature
			double energy_remain = energy - energy_loss;

			// Update soil temperature
			temp = temp_old + energy_remain / Ch;
		}

	} // End of Thawing

	// Freezing
	if (ifFreeze) {

		// Update density as depth doesn't change
		if (!negligible(thickness)) {
			dens = imass * rho_H2O / thickness;
		}
	}
	// Thawing
	else if (ifThaw) {
		
		// Update depth as density doesn't change
		thickness = imass / (dens / rho_H2O);

		// Update water holding capacity. Water exceeding this limit will be handled later.
		update_whc();
	}

	if (temp > 0.0 && imass > 0.0)
		dprintf("ERROR: SNOW Phase Change Year %d Day %d Temp %g imass %g\n", date.year, date.day, temp, imass);

	// No negative wmass
	if (wmass < -EPS)
		fail("ERROR: Negative wmass at end of snow phase_change! Year %d Day %d wmass %g\n", date.year, date.day, wmass);

	// No negative imass
	if (imass < -EPS)
		fail("ERROR: Negative imass at end of snow phase_change! Year %d Day %d imass %g\n", date.year, date.day, imass);

	double wtotal_final = wmass + imass;

	// water balance checks:
	if ((fabs(wtotal_initial - wtotal_final)) > EPS) {
		fail("ERROR: Water mass balance during freezing (%d) or thawing (%d) in Snow::phase_change. Init %g, Final %g, Diff %g \n", 
			ifFreeze, ifThaw, wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
	}

	// Update thermal properties
	update_thermal_properties();
}

/// Serialize snow layer variables
/** Serialize snow layer variables
 */
void Snowlayer::serialize(ArchiveStream& arch) {
	arch & imass
		& wmass
		& dens
		& thickness
		& whc
		& temp;
}

////////////////////////////////////////////////////////////////////////////////
// SNOW 

/// Initial state of snow 
/** Initial state of snow 
 */
void Snow::init_state() {

	nsnowlayers = 0;
	NO3_mass = 0.0;
	NH4_mass = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mdepth[mth] = 0.0;
		mcover[mth] = 0.0;
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			mtemp[lyr][mth] = 0.0;
		}
	}

	for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
		snow_excess_frac[lyr] = (sin(PI * (2.0 * lyr - (NSNOWLAYERS_MAX - 1.0) / 2.0) / (NSNOWLAYERS_MAX - 1.0)) + 1.0) / (NSNOWLAYERS_MAX - 1.0);
	}
}

/// Function returns total available N mass in snow 
/** Function returns total available N mass in snow 
 */
double Snow::nmass_avail() {
	return NO3_mass + NH4_mass;
}

/// Function returns total C mass in snow 
/** Function returns total C mass in snow 
 */
double Snow::ccont() const {
	return 0.0;
}

/// Function returns total N mass in snow 
/** Function returns total N mass in snow 
 */
double Snow::ncont() const {
	return NO3_mass + NH4_mass;
}

/// Function returns total water and ice mass in snow 
/** Function returns total water and ice mass in snow 
 */
double Snow::water_and_ice_content() const {
	double wcontent = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		const Snowlayer& snowlayer = snowlayers[lyr];
		wcontent += snowlayer.wmass + snowlayer.imass;
	}
	return wcontent;
}

/// Function returns total water and ice mass in snow 
/** Function returns total water and ice mass in snow 
 */
double Snow::wmass() {
	double water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		water += snowlayer.wmass;
	}
	return water;
}

/// Function returns total water and ice mass in snow 
/** Function returns total water and ice mass in snow 
 */
double Snow::imass() {
	double ice = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		ice += snowlayer.imass;
	}
	return ice;
}

/// Function returns total snow depth
/** Function returns total snow depth
 */
double Snow::depth() {
	double sd = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sd += snowlayers[lyr].thickness;
	}
	return sd;
}

/// Function returns average snow depth from the three months with highest snow depth
/** Function returns average snow depth 
 * from the three months with highest snow depth
 */
double Snow::snowpack_avrmax3mth() {

	double mdepth_sort[12];

	// Copy values from mdepth to local variable mdepth_sort
	std::copy(mdepth, mdepth + 12, mdepth_sort);

	// Sort depths ascending
	std::sort(mdepth_sort, mdepth_sort + 12);

	// Return average snow depth of the three months of highest snow depth
	return (mdepth_sort[9] + mdepth_sort[10] + mdepth_sort[11]) / 3.0;
}

/// Function updates heat capacity of each snow layer
/** Function updates heat capacity of each snow layer
 */
void Snow::update_thermal_properties() {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].update_thermal_properties();
	}
}

/// Potential redistribution of ice and density after melt
/** Potential redistribution of ice and density after melt
 */
void Snow::potential_depth(double snowpack_depth, double pot_depth[NSNOWLAYERS_MAX], int& pot_nsnowlayers) {

	pot_nsnowlayers = 0;

	// All snow layer together exceed standard depth
	if (snowpack_depth >= SNOW_DEPTH_STANDARD * NSNOWLAYERS_MAX) {

		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			pot_depth[lyr] = SNOW_DEPTH_STANDARD;
		}

		double left_depth = snowpack_depth - SNOW_DEPTH_STANDARD * NSNOWLAYERS_MAX;

		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			pot_depth[lyr] += snow_excess_frac[lyr] * left_depth;
		}

		pot_nsnowlayers = NSNOWLAYERS_MAX;
	}
	else {
		// fill layers so that depths are between half SNOW_DEPTH_STANDARD and SNOW_DEPTH_STANDARD
		int lyr = 0;
		while (!negligible(snowpack_depth)) {
			if (snowpack_depth > 2.0 * SNOW_DEPTH_STANDARD) {
				pot_depth[lyr] = SNOW_DEPTH_STANDARD;
				snowpack_depth -= SNOW_DEPTH_STANDARD;
				pot_nsnowlayers++;
				lyr++;
			}
			else if (snowpack_depth <= SNOW_DEPTH_STANDARD) {
				pot_depth[lyr] = snowpack_depth;
				snowpack_depth = 0.0;
				pot_nsnowlayers++;
			}
			else if (snowpack_depth <= 2.0 * SNOW_DEPTH_STANDARD) {
				pot_depth[lyr] = snowpack_depth / 2.0;
				pot_depth[lyr + 1] = snowpack_depth / 2.0;
				snowpack_depth = 0.0;
				pot_nsnowlayers += 2;
			}
		}
	}
}

/// Re-distributing snow after freeze, melt or compaction event
/** Re-distributing snow after freeze, melt or compaction event
 */
void Snow::redist_layers(double snowpack_depth, double snowfall_depth, 
	double snowfall_imass, double snowfall_temp, double snowfall_Ch) {

	const double EPS = 1.0e-12;
	int lyr;

	// Save current snowpack status
	double snowlyr_thickness[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_imass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_wmass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_temp[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_Ch[NSNOWLAYERS_MAX] = { 0.0 };

	// Water in snow layer that has no ice after phase change, hence no depth (a full layer has melted)
	double wmass_no_snowlyr = 0.0;

	// Extract information from all snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		Snowlayer& snowlayer = snowlayers[lyr];

		snowlyr_thickness[lyr] = snowlayer.thickness;
		snowlyr_imass[lyr] = snowlayer.imass;
		if (snowlayer.thickness > 0.0) {
			snowlyr_wmass[lyr] = snowlayer.wmass;
		}
		else {
			wmass_no_snowlyr += snowlayer.wmass;
		}

		snowlyr_temp[lyr] = snowlayer.temp;
		snowlyr_Ch[lyr] = snowlayer.Ch;

		snowlayer.thickness = snowlayer.imass = snowlayer.wmass = snowlayer.temp = snowlayer.Ch = 0.0;
	}

	// If water in snow layer without depth, then add water to top snowlayer
	if (wmass_no_snowlyr > 0.0) {
		int lyr = nsnowlayers - 1;
		while (wmass_no_snowlyr > 0.0 && lyr >= 0) {
			if (snowlyr_thickness[lyr] > 0.0) {
				snowlyr_wmass[lyr] += wmass_no_snowlyr;
				wmass_no_snowlyr = 0.0;
			}
			else {
				lyr--;
			}
		}
	}

	// Potential depth array
	double pot_thickness[NSNOWLAYERS_MAX] = { 0.0 };

	// Potential number of snow layers
	int pot_nsnowlayers;

	// Check potential snow layer depth
	potential_depth(snowpack_depth + snowfall_depth, pot_thickness, pot_nsnowlayers);

	// Create new snow layer if needed
	for (lyr = nsnowlayers; nsnowlayers < pot_nsnowlayers; lyr++) {
		snowlayers.resize(nsnowlayers + 1);
		snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, 0.0);
		nsnowlayers++;
	}

	// Movement variables
	double frac_move, diff_lyr;
	int lyr_add;
	bool diff;

	// Recalculate ice & dens of layers
	for (lyr = 0; lyr < nsnowlayers; lyr++) {

		Snowlayer& snowlayer = snowlayers[lyr];

		diff_lyr = pot_thickness[lyr] - snowlayer.thickness;

		lyr_add = 0;
		diff = true;

		// Fill every snow layer to pot_ice from bottom up
		while (diff && lyr_add < nsnowlayers) {

			if (!negligible(snowlyr_thickness[lyr_add])) {

				// Frac of lyr_add to add to this snow layer
				frac_move = min(1.0, diff_lyr / snowlyr_thickness[lyr_add]);

				// Update snow layer temperature
				if (negligible(snowlayer.thickness)) {
					snowlayer.temp = snowlyr_temp[lyr_add];
				}
				else {
					snowlayer.update_temp(snowlyr_temp[lyr_add], frac_move * (snowlyr_imass[lyr_add] + snowlyr_wmass[lyr_add]), snowlyr_Ch[lyr_add]);
				}

				// Add snow
				snowlayer.thickness += frac_move * snowlyr_thickness[lyr_add];
				snowlayer.imass += frac_move * snowlyr_imass[lyr_add];
				snowlayer.wmass += frac_move * snowlyr_wmass[lyr_add];

				// Calculate new density
				snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.thickness;

				// Update thermal properties
				snowlayer.update_thermal_properties();

				// Subtract added snow
				diff_lyr -= frac_move * snowlyr_thickness[lyr_add];
				snowlyr_thickness[lyr_add] *= (1.0 - frac_move);
				snowlyr_imass[lyr_add] *= (1.0 - frac_move);
				snowlyr_wmass[lyr_add] *= (1.0 - frac_move);

				if (negligible(diff_lyr)) {
					diff = false;
				}
				else {
					lyr_add++;
				}
			}
			else {
				lyr_add++;
			}
		}
	}

	// Add snowfall
	lyr = nsnowlayers;
	while (!negligible(snowfall_depth)) {

		lyr--;

		if (lyr >= 0) {

			Snowlayer& snowlayer = snowlayers[lyr];

			// Frac of snowfall to add to this snow layer
			frac_move = min(1.0, max(0.0, pot_thickness[lyr] - snowlayer.thickness) / snowfall_depth);

			// Update snow layer temperature
			if (negligible(snowlayer.thickness)) {
				snowlayer.temp = snowfall_temp;
			}
			else {

				snowlayer.update_temp(snowfall_temp, frac_move * snowfall_imass, snowfall_Ch);
			}

			// Add snow
			snowlayer.thickness += frac_move * snowfall_depth;
			snowlayer.imass += frac_move * snowfall_imass;

			// Calculate new density
			snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.thickness;

			// Update thermal properties
			snowlayer.update_thermal_properties();

			// Subtract added snowfall
			snowfall_depth *= (1.0 - frac_move);
			snowfall_imass *= (1.0 - frac_move);
		}
		else {

			// Add remaining (negligible) snowfall to the top snow layer
			int top_lyr = nsnowlayers - 1;

			Snowlayer& snowlayer = snowlayers[top_lyr];
			
			// Update snow layer temperature
			snowlayer.update_temp(snowfall_temp, snowfall_imass, snowfall_Ch);

			// Add snow
			snowlayer.thickness += snowfall_depth;
			snowlayer.imass += snowfall_imass;

			// Calculate new density
			snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.thickness;

			// Update thermal properties
			snowlayer.update_thermal_properties();

			// Remove snowfall
			snowfall_depth = 0.0;
			snowfall_imass = 0.0;
		}
	}

	// Check if snow layers should be erased
	if (nsnowlayers > pot_nsnowlayers) {
		for (lyr = nsnowlayers - 1; lyr >= pot_nsnowlayers; lyr--) {

			Snowlayer& snowlayer_erase = snowlayers[lyr];
			Snowlayer& snowlayer_top = snowlayers[pot_nsnowlayers - 1];

			if (!negligible(snowlayer_erase.wmass + snowlayer_erase.imass)) {
				snowlayer_top.update_temp(0.0, snowlayer_erase.wmass + snowlayer_erase.imass, Cp_water);
				snowlayer_top.wmass += snowlayer_erase.wmass + snowlayer_erase.imass;
			}
			snowlayers.erase(snowlayers.begin() + lyr);
		}
		nsnowlayers = static_cast<int>(std::round(snowlayers.size()));
	}
}

/// Reduce snowpack to its maximum depth, excess goes to runoff
/** Reduce snowpack to its maximum depth, excess goes to runoff
 *  Reduce from bottom of snowpack
 */
void Snow::reduce_snowpack(Patch& patch, double depth_reduction) {

	// Bottom snow layer to be reduced
	int blyr = 0;

	// Runoff
	double runoff_water = 0.0;
	double runoff_ice = 0.0;

	while (!negligible(depth_reduction) && blyr < nsnowlayers) {

		Snowlayer& snowlayer = snowlayers[blyr];

		// Fraction to remove
		double frac = min(1.0, depth_reduction / snowlayer.thickness);

		// Reduce reduction
		depth_reduction -= frac * snowlayer.thickness;

		// Runoff water
		runoff_water += frac * snowlayer.wmass;

		// Runoff ice
		runoff_ice += frac * snowlayer.imass;

		// Reduce snow
		snowlayer.wmass *= (1.0 - frac);
		snowlayer.imass *= (1.0 - frac);
		snowlayer.thickness *= (1.0 - frac);

		blyr++;
	}

	// Add water to surface runoff
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_water + runoff_ice);
}

/// Update snow water holding capacity
/** Update snow water holding capacity
* based on eqn 10 Wang et al. 2013
*/
void Snow::update_whc() {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].update_whc();
	}
}

/// Balance ice and water in snow layers after melting and refreezing
/** Balance ice and water in snow layers after melting and refreezing
 */
void Snow::balance_snow(Patch& patch, double& rain_melt, double& rain_melt_temp) {

	// No snow layers or not multilayer snow scheme then leave function
	if (!nsnowlayers || !ifmultilayersnow)
		return;

	const double EPS = 1.0e-11;
	bool snowpack_removed = false;
	SSW& ssw = patch.soil.ssw;

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];

		snowlayer.thickness = snowlayer.imass / (snowlayer.dens / rho_H2O);
	}

	// Balance layers depending on if snowpack is over maximum or not
	double snowpack_depth = depth();

	// Too thick snowpack (max 2m), remove some ice/water as runoff to 
	// not have a runaway effect of snowpack depth
	if (snowpack_depth > MAX_SNOWPACK_DEPTH) {
		reduce_snowpack(patch, snowpack_depth - MAX_SNOWPACK_DEPTH);
		snowpack_depth = MAX_SNOWPACK_DEPTH;
	}

	// Redistribute snow layers if snow exist
	if (!negligible(snowpack_depth)) {
		redist_layers(snowpack_depth);
	}
	// Insignificant snow depth. Remove all snow layers and add it to standing water
	else if (nsnowlayers) {
		for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

			Snowlayer& snowlayer = snowlayers[lyr];

			double added_wmass = snowlayer.imass + snowlayer.wmass;
			if (!negligible(added_wmass)) {

				// Add excess water to surface standing water
				ssw.update_temp(snowlayer.temp, added_wmass);

				// Update standing water content (mm)
				ssw.wmass += added_wmass;
			}

			snowlayers.erase(snowlayers.begin() + lyr);
		}
		nsnowlayers = 0;
		snowpack_removed = true;
	}

	// Update snow water holding capacity
	update_whc();

	// Save this timesteps extra melt water either from insignificant snow layers or that snow can't hold all water
	double ts_melt = 0.0;
	double ts_melt_temp = 0.0;

	// Redistribute water within the snow pack
	// moving from top to bottom snow layers
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		Snowlayer& snowlayer = snowlayers[lyr];

		if (!negligible(ts_melt)) {
			snowlayer.update_temp(ts_melt_temp, ts_melt, Cp_water);
			snowlayer.wmass += ts_melt;
		}

		if (snowlayer.wmass > snowlayer.whc) {
			ts_melt = max(0.0, snowlayer.wmass - snowlayer.whc);
			snowlayer.wmass -= ts_melt;
			ts_melt_temp = max(0.0, snowlayer.temp);
		}
		else {
			ts_melt = 0.0;
			ts_melt_temp = 0.0;
		}

		// Update thermal properties
		snowlayer.update_thermal_properties();
	}

	// Add excess water to standing water
	ssw.update_temp(ts_melt_temp, ts_melt);

	// Update standing water content (mm)
	ssw.wmass += ts_melt;
}

/// Phase change within the snow pack
/** Phase change within the snow pack
 * Calculate water, ice and air fraction in each layer
 * and snow phase change
 */
void Snow::phase_change() {

	// No phase change in simple snow scheme
	if (!ifmultilayersnow)
		return;

	// Phase change for each snow layer
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		if (lyr == nsnowlayers - 1) {
			snowlayers[lyr].phase_change(true);
		}
		else {
			snowlayers[lyr].phase_change();
		}
	}
}

/// Soak snow pack with water
/** Soak snow pack with water
 * Water soaks snowpack from bottom up
 */
void Snow::soak_snowpack(double& wmass_ssw, double& imass_ssw, double temp_ssw) {

	// Added excess ice to the bottom of the snowpack
	if (imass_ssw > 0.0) {
		Snowlayer& snowlayer = snowlayers[0];

		// Update temperature due to added water
		snowlayer.update_temp(temp_ssw, imass_ssw, Cp_ice);

		// Add ice to snow layer
		snowlayer.imass += imass_ssw;

		// Update snow layer thickness
		snowlayer.thickness += imass_ssw;

		// Update thermal properties
		snowlayer.update_thermal_properties();

		// Set standing ice to zero
		imass_ssw = 0.0;
	}

	// If there was water exceeding whc and couldn't be infiltrated to the soil, 
	// then it will fill the snowpack from bottom up

	int lyr = 0;
	while (wmass_ssw > 0.0 && lyr < nsnowlayers) {

		Snowlayer& snowlayer = snowlayers[lyr];

		// Water to add within water holding capacity
		double add_wmass = max(0.0, min(wmass_ssw, snowlayer.thickness - (snowlayer.wmass + snowlayer.imass)));

		// Update temperature due to added water
		snowlayer.update_temp(temp_ssw, add_wmass, Cp_water);

		// Add water to snow layer
		snowlayer.wmass += add_wmass;

		// Update thermal properties
		snowlayer.update_thermal_properties();

		// Remove added water from standing water
		wmass_ssw -= add_wmass;

		// Next layer
		lyr++;
	}

	// Added excess water to the top of the snowpack
	if (wmass_ssw > 0.0) {
		Snowlayer& snowlayer = snowlayers[nsnowlayers - 1];

		// Update temperature due to added water
		snowlayer.update_temp(temp_ssw, wmass_ssw, Cp_water);

		// Add water to snow layer
		snowlayer.wmass += wmass_ssw;

		// Update snow layer thickness
		snowlayer.thickness += wmass_ssw;

		// Update thermal properties
		snowlayer.update_thermal_properties();

		// Set standing water to zero
		wmass_ssw = 0.0;
	}
}

/// Snow water over water holding capacity
/** Snow water over water holding capacity
 *  Sum up water above water holding capacity of snow
 *  that will be allowed to infiltrate or runoff 
 */
double Snow::wmass_excess() {

	double wmass = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		wmass += max(0.0, snowlayers[lyr].wmass - snowlayers[lyr].whc);
	}
	return wmass;
}


/// Snow compaction and update density and depth
/** Snow compaction and update density and depth
 */
void Snow::compaction(Patch& patch, double& rain_melt, double& rain_melt_temp) {

	// No compaction in simple snow scheme
	if (!ifmultilayersnow)
		return;

	// Increase in layer density due to compaction (Best et. al. 2011, Eq.21)
	double k_s = 4000.0;
	// Compactive viscosity factor for snow
	double eta_0 = 1E7;

	// Snow mass above current layer
	double I_above = 0.0;

	// Mechanical compaction scheme, following Best et al. 2011
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		Snowlayer& snowlayer = snowlayers[lyr];

		// Current layers ice and water mass
		double I_sum = snowlayer.imass + snowlayer.wmass;

		// Cumulative snow mass affecting compaction
		double M_k = 0.5 * I_sum + I_above;

		// Layer temperature in Kelvin (never above zero degree Celcius)
		double T_k = min(K2degC, snowlayer.temp + K2degC);

		if (snowlayer.dens < rho_ICE) {
			double delta_dens = snowlayer.dens * gravity * M_k / eta_0 *
				exp((k_s / K2degC) - (k_s / T_k) - (snowlayer.dens / rho_SNOW_0)) * SECS_PER_DAY;

			// Increase density
			snowlayer.dens += delta_dens;

			// Snow density can not exceed pure ice density
			snowlayer.dens = min(rho_ICE, snowlayer.dens);
		}

		// Add this layer to the above layers for future calculations
		I_above += I_sum;
	}

	// Calculate snow layers thickness after compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		snowlayer.thickness = snowlayer.imass / (snowlayer.dens / rho_H2O);
	}
}

/// Handle rain on snow event and rainwater distribution within the snowpack
/** Handle rain on snow event and rainwater distribution within the snowpack
 * Calculate water, ice and air fraction in each layer
 * and snow phase change in layer based on the available energy
 */
void Snow::rainwater_melt(double airT, double& prec, double& rain_melt, double& rain_melt_temp) {

	// Water transferred between snow layer and in the end to the soil
	double water_transferred = prec;

	// Temperature difference between air temperture and zero
	double delta_T = airT - 0.0;
	
	// Total energy from rain [J]
	double rain_energy = Cp_water * delta_T * prec * M_PER_MM;

	// Top snow lyr index
	int lyr = nsnowlayers - 1;

	// Loop through layers as long as water still has a temperture above zero (still energy to melt ice) 
	while (rain_energy > 0.0 && !negligible(snowlayers[lyr].imass)) {

		Snowlayer& snowlayer = snowlayers[lyr];

		// Heat capcity of snow layer [J m-3 K-1]
		double heat_cap = snowlayer.Ch;
		// Volume of snow layer [m-3]
		double snow_volume = snowlayer.thickness * M_PER_MM;
		// Energy needed to bring snow layer to 0C [J]
		double energy_to_zero = fabs(snowlayer.temp - 0.0) * heat_cap * snow_volume;

		// Enough energy to melt some snow
		if (rain_energy > energy_to_zero) {

			// Snow layer temp is set to zero and erase energy needed to bring up temperature
			rain_energy -= energy_to_zero;
			snowlayer.temp = 0.0;

			// ice melt [mm]
			double imass_melt = rain_energy / Lheat * MM_PER_M;

			// Enough energy to melt all ice (snow layer will be removed) 
			if (imass_melt > snowlayer.imass) {
				imass_melt -= snowlayer.imass;
				snowlayer.wmass += snowlayer.imass;
				snowlayer.imass = 0.0;
				rain_energy = imass_melt * Lheat / MM_PER_M;
			}
			// Only melt a frac of ice
			else {
				snowlayer.wmass += imass_melt;
				snowlayer.imass -= imass_melt;
				rain_energy = 0.0;
			}
		}
		// Only enough energy to change temperature
		else {
			snowlayer.temp += rain_energy / (snow_volume * heat_cap);
			rain_energy = 0.0;
		}

		// If layer still exist update thickness, water content and temperature
		if (!negligible(snowlayer.imass)) {

			// Recalculate thickness
			snowlayer.thickness = snowlayer.imass / (snowlayer.dens / rho_H2O);

			// Recalculate water holding capacity of snow layer
			snowlayer.update_whc();

			// Transferred water added to this layer
			double water_added = min(water_transferred, snowlayer.whc - snowlayer.wmass);

			water_transferred -= water_added;

			// Update temperature assuming water temperature is 0C
			if (water_added > 0.0) {
				snowlayer.update_temp(0.0, water_added, Cp_water);
			}

			snowlayer.wmass += water_added;			
		}
		else {

			// Snow layer will be erased, transfer all water to next snow layer or the soil
			water_transferred += snowlayer.wmass;

			// Erase this empty snow layer
			snowlayers.erase(snowlayers.end()-1);
			nsnowlayers--;

			// If this was the last snow layer, then exit loop
			if (nsnowlayers == 0)
				break;
		}

		// Go thorugh lyrs while rain_energy not 0
		lyr--;

		// This was the last snow layer
		if (lyr < 0) {
			break;
		}
	} // End of while loop

	// Still energy left, update melted water temperature
	double water_transferred_temp = 0.0;
	if (!negligible(rain_energy)) {
		water_transferred_temp = rain_energy / (Cp_water * water_transferred * M_PER_MM);
	}

	// If water to be transferred to the soil
	if (!negligible(water_transferred)) {

		// Update rain_melt temperture
		rain_melt_temp = (rain_melt_temp * rain_melt + water_transferred_temp * water_transferred) / (rain_melt + water_transferred);

		// Save melt water above snow whc which will be pass downwards
		rain_melt += water_transferred;
	}
}

/// Serialize snow variables
/** Serialize snow variables
 */
void Snow::serialize(ArchiveStream& arch) {
	arch & snowlayers
		& nsnowlayers
		& NH4_mass
		& NO3_mass;
}

/// Calculate ice impedance factor
/** Calculate ice impedance factor (Swenson et al. 2012)
 */
double ice_imp(double theta_i) {

	 return min(1.0, pow(10.0, -6.0 * theta_i));
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soillayer member functions

/// Soil layer constructor
/** Soil layer constructor: initialises soil, depth, profile depth
 * and if soil layer is an acrotelm or catotelm
 */
Soillayer::Soillayer(Soil* sc, double d, double dp, bool at, bool ct) : soil(sc), thickness(d), depth_mid(dp), acrotelm(at), catotelm(ct) {

	// Initialises certain member variables

	// Energy
	temp = 0.0;
	temp_old = 0.0;
	Kh = 0.0;
	Ch = 0.0;
	thermdiff_0 = 0.0;
	thermdiff_15 = 0.0;
	thermdiff_100 = 0.0;
	Cp_min = 0.0;
	gtemp = 0.0;

	// Water
	wmass = 0.0;
	imass = 0.0;
	wmass_sc = 0.0;
	wsats = 0.0;
	awc = 0.0;
	wp = 0.0;
	Kw = 0.0;
	dKwdw = 0.0;
	Kw_max = 0.0;
	Psi = 0.0;
	dPsidw = 0.0;
	Psi_max = 0.0;
	b = 0.0;
	perc_mass = 0.0;
	runoff_drain = 0.0;
	runoff_baseflow = 0.0;
	mwcontabove = 0.0;
	awcont = 0.0;

	root_volume = 0.0;
	mineral_volume = thickness;
	coarse_volume = 0.0;
	organic_volume = 0.0;
	effective_volume = thickness;
	perc_base = 0.0;
	depth_bot = depth_mid + thickness / 2.0;

	rho_org = rho_high;

	// Leaching
	orgleachfrac = 0.0;
	aorgCleach = 0.0;
	aorgNleach = 0.0;
	aminleach = 0.0;

	// Old SOM module
	alag = 0.0;
	exp_alag = 1.0;
	cpool_slow = 0.0;
	cpool_fast = 0.0;
	decomp_litter_mean = 0.0;
	k_soilfast_mean = 0.0;
	k_soilslow_mean = 0.0;

	// CENTURY SOM pools
	for (int p = 0; p < NSOMPOOL; p++) {
		sompool[p].cmass = 0.0;
		sompool[p].nmass = 0.0;
		sompool[p].ligcfrac = 0.0;
		sompool[p].delta_cmass = 0.0;
		sompool[p].delta_nmass = 0.0;
		sompool[p].fracremain = 0.0;
		sompool[p].litterme = 0.0;
		sompool[p].fireresist = 0.0;
		sompool[p].ntoc = 0.0;
		sompool[p].cdec = 0.0;
		sompool[p].cinc = 0.0;
		sompool[p].ndec = 0.0;
		sompool[p].ninc = 0.0;
		sompool[p].respsum = 0.0;
		sompool[p].respsum_surf = 0.0;
		sompool[p].cmass_init = 0.0;
		sompool[p].multi = 0.0;
		sompool[p].sign_trend = 0;
		sompool[p].same_trend = 0;
		sompool[p].changed_trend = 0;

		for (int m = 0; m < 12; m++) {
			sompool[p].mfracremain_mean[m] = 0.0;
		}
	}

	// Initial CENTURY pool N:C ratios
	// Parton et al 1993, Fig 4

	sompool[SOILMICRO].ntoc = 1.0 / 15.0;
	sompool[SURFHUMUS].ntoc = 1.0 / 15.0;
	sompool[SLOWSOM].ntoc = 1.0 / 20.0;
	sompool[SURFMICRO].ntoc = 1.0 / 20.0;

	// Passive SOM pool has a fixed value
	sompool[PASSIVESOM].ntoc = 1.0 / 9.0;

	// Decay rate modifier depending on depth profile (Koven et al 2013)
	// Reference depth set to 10m
	depth_mod = ifdepthdecaymod ? exp(-depth_mid / 10000.0) : 1.0;
	a = 0.0;

	// N pools
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	NO2_mass = 0.0;
	NO_mass = 0.0;
	N2O_mass = 0.0;
	N2_mass = 0.0;
	anmin = 0.0;
	animmob = 0.0;

	// Gases
	CO2 = 0.0;
	CH4 = 0.0;
	CH4_ebull_ind = 0.0;
	CH4_gas = 0.0;
	CH4_diss = 0.0;
	O2 = 0.0;

	// Tiller area
	tiller_area = 0.0001;

	// Heterotrophic respiration
	respsum = 0.0;
	arespsum = 0.0;

	// Decay rates
	acdec[0] = 0.0;
	acdec[1] = 0.0;
	acinc[0] = 0.0;
	acinc[1] = 0.0;
	atemp_mod = 0.0;
	amoist_mod = 0.0;
	atot_mod = 0.0;

	// Soil labile carbon availability
	labile_carbon = 0.0;

	// Soil organic matter age
	socage = 0.0;

	peat_frac = 0.0;
	bulk_density = 0.0;

	// Monthly variables
	for (int mth = 0; mth < 12; mth++) {
		mwmass[mth] = 0.0;
		mimass[mth] = 0.0;
		mwcont[mth] = 0.0;
		mdrain[mth] = 0.0;
		mnetmin[mth] = 0.0;
		mtemp[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
	}
}

/// Calculates the effective organic density
/** Based on Chadburn et al. 2022
 *	Coarse root litter is not affecting soil properties, just taking up space
 */ 
double Soillayer::calc_organic_density(grouppooltype group) {

	double volume_tot = 0.0;
	double cmass_tot = 0.0;

	for (int p = 0; p < NSOMPOOL; p++) {

		if (group == GALL || group == GSURFL) {
			if (p == SURFMETA || p == SURFSTRUCT || p == SURFFWD || p == SURFCWD) {
				volume_tot += sompool[p].cmass / rho_low; // DW_COMMENT need to set rho better!!!
				cmass_tot += sompool[p].cmass;
			}
		}
		
		if (group == GALL || group == GSOM || group == GHUMUSSOM) {
			if (p == SOILMETA || p == SOILSTRUCT) {
				volume_tot += sompool[p].cmass / rho_low;
				cmass_tot += sompool[p].cmass;
			}
			else if (p == SOILMICRO || p == SLOWSOM) {
				volume_tot += sompool[p].cmass / rho_mean;
				cmass_tot += sompool[p].cmass;
			}
			else if (p == SOILFWD || p == SOILCWD || p == PASSIVESOM) {
				volume_tot += sompool[p].cmass / rho_high;
				cmass_tot += sompool[p].cmass;
			}
		}

		if (group == GALL || group == GHUMUS || group == GHUMUSSOM) {
			if (p == SURFMICRO) {
				volume_tot += sompool[p].cmass / rho_mean;
				cmass_tot += sompool[p].cmass;
			}
			else if (p == SURFHUMUS) {
				volume_tot += sompool[p].cmass / rho_high;
				cmass_tot += sompool[p].cmass;
			}
		}
	}

	return !negligible(volume_tot) ? cmass_tot / volume_tot : rho_high;
}

// Mineral exponential in Campbell (1974) parameterization for soil 
// water retention curves from Table 4, Cosby (1984)
double Soillayer::mineral_exp(double clay_frac, double sand_frac) {
	return 3.10 + 15.7 * clay_frac - 0.3 * sand_frac;
}

// Set soil layer properties
/** Set soil layer properties
 * Sets the following properties for this individual layer
 * clay, silt and sand fraction
 * pH
 * hydrologic parameters (wsats, awc, wp, Psi_max, b and Kw_max)
 * thermal properties for analytic temperature calculation
 * and two layer soil perculation parameters
*/
void Soillayer::set_properties() {

	// Get soil type
	Soiltype soiltype = (*soil).soiltype;

	if (soiltype.soil_input == LPJ) {

		// Water holding capacities for mineral
		wsats_min = soiltype.wsats_min;
		awc_min = soiltype.awc_min;
		wp_min = soiltype.wp_min;

		// Percolation coefficient
		perc_base_min = soiltype.perc_base;

		// Mineral exponential in Campbell (1974) parameterization for soil 
		// water retention curves from Table 4, Cosby (1984)
		b_min = mineral_exp(clay_frac, sand_frac);
	}
	else if (soiltype.soil_input == MINERAL || soiltype.soil_input == VARYING) {

		// Equation 1 from Cosby 1984
		// Psi = Psi_s * (Theta/Theta_s)^b
		// Psi is the pressure head in cm
		// *_s is the values at saturation
		// Theta is the volumetric moisture content
		// Re-arranged to get the Theta
		// Theta = Theta_s * (Psi/Psi_s)^(-1/b)

		// Mineral exponential in Campbell (1974) parameterization for soil 
		// water retention curves from Table 4, Cosby (1984)
		b_min = mineral_exp(clay_frac, sand_frac);

		// Saturated soil water suction (cm)
		double Psi_s = -pow(10.0, 1.54 - 0.95 * sand_frac + 0.63 * silt_frac);

		// Volumetric soil water concentration at saturation point (porosity)
		double Theta_s = 0.505 - 0.142 * sand_frac - 0.037 * clay_frac;

		// Volumetric soil moisture concentration at Field Capacity
		// Field capacity generally occurs at a suction of 100 cm
		double Theta_whc = Theta_s * pow(-100.0 / Psi_s, -1.0 / b_min);

		// Volumetric soil moisture concentration at wilting point
		// Wilting point typically occurs at a suction of 15000 cm
		double Theta_wilt = Theta_s * pow(-15000.0 / Psi_s, -1.0 / b_min);

		// Water holding capacities
		wsats_min = Theta_s;
		awc_min = (Theta_whc - Theta_wilt);
		wp_min = Theta_wilt;

		// Percolation coefficient
		perc_base_min = 5.87 - 0.29 * b_min;

		// Analytic thermal diffusivity parameters
		thermdiff_0_min = 0.2;
		thermdiff_15_min = 0.15 * b_min + 0.05;
		thermdiff_100_min = 0.4;
	}

	// Heat capacity of mineral soil [J m-3 K-1] - Bonan (2002)
	Cp_min = (2.128 * sand_frac + 2.385 * clay_frac) / (sand_frac + clay_frac) * 1.0e6;// 2380000;

	// Saturated mineral hydraulic conductivity (from inch h-1 to mm s-1) from Table 4, Cosby 1984
	Kw_max_min = pow(10.0, -0.6 + 1.26 * sand_frac - 0.64 * clay_frac) * MM_PER_INCH / SECS_PER_HOUR;
	
	// The saturated hydraulic conductivity decreases exponentially with depth (Niu et al. 2003)
	if (depth_mid > soiltype.text_depths.back()) {
		
		// Global exponential decay factor, f_min, of the saturated hydraulic conductivity (Niu et al. 2003)
		double f_min = 0.002;
		// Minimum saturated hydraulic conductivity
		double Kw_min_min = 1.0e-8;
		// Updated saturated hydraulic conductivity due to depth from last measured texture of the soil
		Kw_max_min = Kw_min_min + (Kw_max_min - Kw_min_min) * exp(-f_min * (depth_mid - soiltype.text_depths.back()));
	}

	// Saturated soil water suction (from cm to mm) from Table 4, Cosby 1984
	Psi_max_min = -pow(10.0, 1.54 - 0.95 * sand_frac + 0.63 * silt_frac) * MM_PER_CM;

	// Super cooled water setting 2
	// Soil specific surface area (m2 g-1) from eqn 11 in Sepaskhah et al. (2011)
	// eqn 11
	double D = 3.0 - 0.118 * (-log(clay_frac) + silt_frac + sand_frac);
	// eqn 13
	ssa = exp(-25.06 + 10.66 * D);

	// Soil organic C microorganism collocation factor for decomposition (Yan 2018)
	a = !(acrotelm || catotelm) ? (clay_frac < 0.0165 ? 0.0 : (clay_frac > 0.373 ? 1.0 : 2.8 * clay_frac - 0.046)) : 0.0;


	// Organic

	// Organic exponential in Campbell (1974) parameterization for soil 
	// water retention curves from CLM5.0 documentation 
	b_org = min(2.7 + 9.3 * depth_mid / sapric_depth, 12.0);

	// Organic soil matric potential at saturation 
	Psi_max_org = -max(10.3 - 0.2 * depth_mid / sapric_depth, 10.1);

	// Organic volumetric soil water concentration at saturation point (porosity)
	wsats_org = max(0.93 - 0.1 * depth_mid / sapric_depth, 0.83);

	// Volumetric soil moisture concentration at Field Capacity
	// Field capacity generally occurs at a suction of 100 cm
	double Theta_whc = wsats_org * pow(-100.0 / Psi_max_org, -1.0 / b_org);

	// Volumetric soil moisture concentration at wilting point
	// Wilting point typically occurs at a suction of 15000 cm
	wp_org = wsats_org * pow(-15000.0 / Psi_max_org, -1.0 / b_org);

	// Organic available water holding capacity of soil layer (-)
	awc_org = Theta_whc - wp_org;

	// Hydraulic conductivity at saturation for organic matter (mm s-1) (CLM5.0 - Documentation; Hoag and Price 1995)
	Kw_max_org = max(0.28 - 0.2799 * depth_mid / sapric_depth, 1.0e-6);


	// Coarse fragment

	/// Coarse fragment (rock) available water holding fraction of soil layer (-)
	awc_rock = 0.0;
	/// Coarse fragment (rock) wilting point of soil layer (-)
	wp_rock = 0.0;
	/// Coarse fragment (rock) saturation point (-)
	wsats_rock = 0.0;

	// Coarse fragment volume
	coarse_volume = !(acrotelm || catotelm) ? coarse_frac * thickness : 0.0;
}

static double calc_weighted_avg(
	const std::vector<double>& prop_lyr,
	const std::vector<double>& depth_lyr,
	double z_top, double z_bottom) {

	double weighted_sum = 0.0;
	double overlap_sum = 0.0;
	double layer_top = 0.0;
	const int nr_lyr = (int)prop_lyr.size();

	for (int lyr = 0; lyr < nr_lyr; lyr++) {
		double layer_bottom = depth_lyr[lyr];
		double prop = prop_lyr[lyr];

		double overlap_top = max(layer_top, z_top);
		double overlap_bottom = min(layer_bottom, z_bottom);
		double overlap = overlap_bottom - overlap_top;

		if (overlap > 0.0) {
			weighted_sum += prop * overlap;
			overlap_sum += overlap;
		}

		layer_top = layer_bottom;

		if (layer_top >= z_bottom)
			break;
	}

	// If the target layer extends below the last soiltype layer,
	// add the remaining thickness weighted by the last property value
	if (z_bottom > layer_top && nr_lyr > 0) {
		double remaining_thickness = z_bottom - max(layer_top, z_top);
		if (remaining_thickness > 0.0) {
			weighted_sum += prop_lyr.back() * remaining_thickness;
			overlap_sum += remaining_thickness;
		}
	}

	return (overlap_sum > 0.0) ? (weighted_sum / overlap_sum) : 0.0;
}


// Initiate soil layer properties
/** Initiate soil layer properties
 * Sets the following properties for this individual layer
 * clay, silt and sand fraction
 * pH
 * hydrologic parameters (wsats, awc, wp, Psi_max, b and Kw_max)
 * thermal properties for analytic temperature calculation
 * and two layer soil perculation parameters
*/
void Soillayer::initiate_properties() {

		// Get soil type
	Soiltype soiltype = (*soil).soiltype;

	// Analytic temp calc parameters
	thermdiff_0_org = soiltype.thermdiff_0_org;
	thermdiff_15_org = soiltype.thermdiff_15_org;
	thermdiff_100_org = soiltype.thermdiff_100_org;
	thermdiff_0_min = soiltype.thermdiff_0;
	thermdiff_15_min = soiltype.thermdiff_15;
	thermdiff_100_min = soiltype.thermdiff_100;

	// Percolation coefficients
	perc_exp = soiltype.perc_exp;
	perc_base_org = soiltype.perc_base_org;

	if (soiltype.soil_input == LPJ || soiltype.soil_input == MINERAL) {

		clay_frac = soiltype.clay_frac.front();
		sand_frac = soiltype.sand_frac.front();
		silt_frac = 1.0 - (clay_frac + sand_frac);
		coarse_frac = 0.0;
		pH = soiltype.pH.front();
	}
	else if (soiltype.soil_input == VARYING) {

		double z_top = max(0.0, depth_bot - thickness);
		double z_bottom = depth_bot;

		clay_frac = calc_weighted_avg(soiltype.clay_frac, soiltype.text_depths, z_top, z_bottom);
		sand_frac = calc_weighted_avg(soiltype.sand_frac, soiltype.text_depths, z_top, z_bottom);
		silt_frac = 1.0 - (clay_frac + sand_frac);
		coarse_frac = !(acrotelm || catotelm) ? calc_weighted_avg(soiltype.coarse_frac, soiltype.text_depths, z_top, z_bottom) : 0.0;
		pH = calc_weighted_avg(soiltype.pH, soiltype.text_depths, z_top, z_bottom);
	}

	set_properties();

	// Mineral volume
	mineral_volume = thickness - coarse_volume;
}


// Set soil layer properties
/** Set soil layer properties
 * Sets the following properties for this individual layer
 * clay, silt and sand fraction
 * pH
 * hydrologic parameters (wsats, awc, wp, por, Psi_max, b and Kw_max)
 * mineral, organic and peat fractions
 * thermal properties for analytic temperature calculation
 * and two layer soil perculation parameters
*/
void Soillayer::update_properties(double& wmass_move, double& imass_move, bool must_set) {

	// Only enter first year except if dynamic soil organic C influence soil properties
	if (!must_set && !ifeffvolume && !iforganicsoilproperties && !(*soil).patch.stand.first_year == date.year) { // DW_COMMENT restart years needs to update it then as well!
		return;
	}

	// Mineral particle density set to 2.65 g cm�3 (Bonan 2019)
	double min_density = Mindensity * G_PER_KG * M3_PER_CM3;

	// Effective soil volume when excluding volume roots and rock take up 
	if (ifeffvolume) {
		effective_volume = thickness - root_volume - coarse_volume;
	}
	else {
		effective_volume = thickness;
	}

	if (!iforganicsoilproperties) {

		// Acrotelm
		if (acrotelm) {
			wsats = acrotelm_por * effective_volume;
			// Field capacity is set at 75% of saturation (Chaudhary 2017)
			awc = (acrotelm_por * 0.75 - peat_wp) * effective_volume;
			wp = peat_wp * effective_volume;
			peat_frac = 1.0 - acrotelm_por;
			bulk_density = peat_frac * rho_org * G_PER_KG * M3_PER_CM3;

			Kw_max = Kw_max_org;
			b = b_org;
			Psi_max = Psi_max_org;

			// Initially saturated
			if ((*soil).patch.stand.first_year == date.year) {
				double wadded = wsats - (wmass + imass);
				wmass += wadded;
				(*soil).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
			}
		}
		// Catotelm
		else if (catotelm) {
			wsats = catotelm_por * effective_volume;
			// Field capacity is set at 75% of saturation (Chaudhary 2017)
			awc = (catotelm_por * 0.75 - peat_wp) * effective_volume;
			wp = peat_wp * effective_volume;
			peat_frac = 1.0 - catotelm_por;
			bulk_density = peat_frac * rho_org * G_PER_KG * M3_PER_CM3;

			Kw_max = Kw_max_org;
			b = b_org;
			Psi_max = Psi_max_org;

			// Initially saturated
			if ((*soil).patch.stand.first_year == date.year) {
				double wadded = wsats - (wmass + imass);
				wmass += wadded;
				(*soil).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
			}
		}
		// Mineral soil and no influence of soil organic C
		else {

			// added moved water and ice due to shrinking wsats
			wmass += wmass_move;
			imass += imass_move;
			wmass_move = imass_move = 0.0;

			wsats = wsats_min * effective_volume;
			awc = awc_min * effective_volume;
			wp = wp_min * effective_volume;
			peat_frac = 0.0;
			bulk_density = min_density * (effective_volume - wsats) / effective_volume;

			Kw_max = Kw_max_min;
			b = b_min;
			Psi_max = Psi_max_min;

			// Move water downwards if water content exceeds maximum water content
			double wmass_max_lyr = wmass_max();
			if (wmass_max_lyr < wmass + imass) {
				double frac = (wmass + imass - wmass_max_lyr) / (wmass + imass);
				wmass_move = wmass * frac;
				imass_move = imass * frac;
				wmass -= wmass_move;
				imass -= imass_move;
			}
		}
	}
	// Soil organic content influences soil properties
	// Hudson (1994) showed that for each 1-percent increase in soil organic matter, the available water holding 
	// capacity in the soil increased by 3.7 percent. Soil water is held by adhesive and cohesive forces within 
	// the soil and an increase in pore space will lead to an increase in water holding capacity of the soil.
	else {

		// added moved water and ice due to shrinking maximum water content
		wmass += wmass_move;
		imass += imass_move;
		wmass_move = imass_move = 0.0;

		double org_frac = organic_volume / thickness;
		double min_frac = mineral_volume / thickness;
		double rock_frac = coarse_volume / thickness;

		if (fabs(organic_volume + mineral_volume + coarse_volume + root_volume - thickness) > 1.0e-12)
			dprintf("Year %d Day %d ERROR in soil volume (%g) compared to thickness (%g)\n", 
				date.year,date.day, organic_volume + mineral_volume + coarse_volume + root_volume, thickness);

		wsats = (min_frac * wsats_min + org_frac * wsats_org + rock_frac * wsats_rock) * thickness;
		awc = (min_frac * awc_min + org_frac * awc_org + rock_frac * awc_rock) * thickness;
		wp = (min_frac * wp_min + org_frac * wp_org + rock_frac * wp_rock) * thickness;
		peat_frac = 0.0;

		double density = (min_frac + rock_frac) * min_density + org_frac * rho_org * G_PER_KG * M3_PER_CM3;

		// Bulk density of mineral and organic matter fraction of soillayer excluding roots
		bulk_density = density * (min_frac + rock_frac + org_frac - wsats / thickness);

		// Move water downwards if water content exceeds maximum water content
		double wmass_max_lyr = wmass_max();
		if (wmass_max_lyr < wmass + imass) {
			double frac = (wmass + imass - wmass_max_lyr) / (wmass + imass);
			wmass_move = wmass * frac;
			imass_move = imass * frac;
			wmass -= wmass_move;
			imass -= imass_move;
		}

		// The saturated hydraulic conductivity
		// The soil organic matter fraction is further subdivided into �connected� and �unconnected� 
		// fractions using a result from percolation theory (Stauffer and Aharony 1994, Berkowitz and 
		// Balberg 1992). Assuming that the organic and mineral fractions are randomly distributed 
		// throughout a soil layer, percolation theory predicts that above a threshold value 
		// org_frac = f_threshold, connected flow pathways consisting of organic material only exist
		// and span the soil space. Flow through these pathways interacts only with organic material, 
		// and thus can be described by Kw_max_org (CLM5.0 documentation).
		if (!negligible(org_frac)) {
			double perc_norm = 1.101;
			double perc_frac = org_frac >= threshold_frac ? perc_norm * pow(org_frac - threshold_frac, beta_perc) * org_frac : 0.0;
			double uncon_frac = 1.0 - perc_frac - rock_frac;
			double Kw_max_uncon = org_frac < 1.0 ? uncon_frac / ((1.0 - org_frac - rock_frac) / Kw_max_min + (org_frac - perc_frac) / Kw_max_org) : Kw_max_org;
			Kw_max = uncon_frac * Kw_max_uncon + (1.0 - uncon_frac - rock_frac) * Kw_max_org;
		}
		else {
			Kw_max = Kw_max_min;
		}

		// Exponential in Campbell (1974) parameterization for soil water retention curves
		b = (min_frac * b_min + org_frac * b_org);

		// Saturated soil water suction (mm)
		Psi_max = (min_frac * Psi_max_min + org_frac * Psi_max_org);

		// Soil organic C microorganism collocation factor for decomposition (Yan 2018)
		a = fclay() < 0.016 ? 0.0 : (fclay() > 0.37 ? 1.0 : 2.8 * fclay() - 0.046);

		// Analytic temp calc parameters
		thermdiff_0 = min_frac * thermdiff_0_min + org_frac * thermdiff_0_org;
		thermdiff_15 = min_frac * thermdiff_15_min + org_frac * thermdiff_15_org;
		thermdiff_100 = min_frac * thermdiff_100_min + org_frac * thermdiff_100_org;

		// Two-layer hydrology parameter
		perc_base = min_frac * perc_base_min + org_frac * perc_base_org;
	}

	// For Gerten and Miller hydrology schemes wmass should never be below wp
	if (hydrology != RICHARDS) {
		double wadded = max(0.0, wp - (wmass+imass));
		wmass += wadded;
		(*soil).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
	}
}

double Soillayer::fclay(mincalctype type) {

	// Absolute value of mineral fraction
	if (type == ABS) {
		return clay_frac;
	}
	// Absolute value of mineral fraction when accounting for effective volume
	else if (type == ABSEFF) {
		return clay_frac * effective_volume / thickness;
	}
	// Value of mineral fraction when only considering the mineral volume of the soil
	else if (type == MINEFF) {
		return clay_frac * mineral_volume / thickness;
	}
	else {
		fail("ERROR! Can't calculate clay fraction\n");
		return clay_frac;
	}
}

double Soillayer::fsand(mincalctype type) {

	// Absolute value of mineral fraction
	if (type == ABS) {
		return sand_frac;
	}
	// Absolute value of mineral fraction when accounting for effective volume
	else if (type == ABSEFF) {
		return sand_frac * effective_volume / thickness;
	}
	// Value of mineral fraction when only considering the mineral volume of the soil
	else if (type == MINEFF) {
		return sand_frac * mineral_volume / thickness;
	}
	else {
		fail("ERROR! Can't calculate sand fraction\n");
		return clay_frac;
	}
}

double Soillayer::fsilt(mincalctype type) {

	// Absolute value of mineral fraction
	if (type == ABS) {
		return silt_frac;
	}
	// Absolute value of mineral fraction when accounting for effective volume
	else if (type == ABSEFF) {
		return silt_frac * effective_volume / thickness;
	}
	// Value of mineral fraction when only considering the mineral volume of the soil
	else if (type == MINEFF) {
		return silt_frac * mineral_volume / thickness;
	}
	else {
		fail("ERROR! Can't calculate silt fraction\n");
		return silt_frac;
	}
}

/// Check if all water should be frozen
/** Check if all water should be frozen
 *  If maximum soil temperature of last year without phase change is below zero, 
 *  then freeze all water as an assumption
 */
void Soillayer::freeze_lyr() {
	int freeze_mths = 0;
	for (int mth = 0; mth < 12; mth++) {
		freeze_mths += mtemp[mth] < 0.0 ? 1 : 0;
	}

	// If monthly mean soil temperature is below freezing half of the year
	// then fill soil layer with half of what is left to saturation and then 
	// freeze all water.
	if (freeze_mths >= 6) {
		double imass_max = hydrology == RICHARDS ? wsats : fc();
		imass = wmass + imass + 0.5 * max(0.0, (imass_max - (wmass + imass)));
		wmass = 0.0;
	}
}

/// Calculate phase change of water and ice for this soil layer
/** Calculate phase change of water and ice for this soil layer
 */
void Soillayer::phase_change() {

	// Initial water and ice content
	double wtotal_initial = wmass + imass;

	// Should ice thaw
	bool ifThaw = false;
	// Should water freeze
	bool ifFreeze = false;
	// Fraction of water which freezes
	double Ffreez = 0.0;
	// Fraction of ice which thaws
	double Fthaw = 0.0;
	// Fraction of water
	double Fwmass = wmass / thickness;
	// Fraction of water super cooled
	double Fwmass_sc = 0.0;
	// Fraction of ice
	double Fimass = imass / thickness;
	// Energy per volume [J m-3]
	double energy;
	// Temperature difference (current - previous)
	double delta_T;
	// Maximum error
	const double EPS = 1.0e-10;

	// Calculate fraction supercooled soil layer water
	// Two options
	// (1) - Niu et al. 2006 eqn 3 (CLM scheme) or
	// (2) - Jin et al. 2020 eqn 28 
	if (ifsupercooled) {

		// Maximum water capacity of soil layer
		double wmass_max_lyr = wmass_max();

		// Temperature below freezing
		if (temp < 0.0) {

			// Maximum super cooled water fraction
			const double max_sc_frac = 0.9;

			// Temperature close to zero, maximum water is super cooled
			if (temp > -0.002) {
				wmass_sc = max_sc_frac * wmass_max_lyr;
			}
			// Niu scheme
			else if (ifsupercooled == 1) {
				wmass_sc = wmass_max_lyr * min(max_sc_frac, pow(Lheat * temp / (gravity * (K2degC + temp) * Psi_max), -1.0 / b));
			}
			// Jin scheme
			else if (ifsupercooled == 2) {
				double w_m = 0.05 * ssa;
				double S = 0.37;
				double C_1 = 0.041;
				double C = ssa * S * C_1 * pow(K2degC - temp, 0.5);
				wmass_sc = wmass_max_lyr * min(max_sc_frac, (C * pow(-temp, -0.5) + w_m - 0.25 * C) / 100.0);
			}

			// Fraction super cooled water
			Fwmass_sc = wmass_sc / thickness;

			// Old temperature can't be above freezing
			temp_old = min(0.0, temp_old);

			// Should water freeze
			if (Fwmass > Fwmass_sc && temp < temp_old) {
				ifFreeze = true;
			}
			// Should ice thaw
			else if (Fwmass < Fwmass_sc && temp > temp_old) {
				ifThaw = true;
			}
		}
		// Should ice thaw
		else if (temp > 0.0 && Fimass > 0.0) {
			ifThaw = true;
			temp_old = 0.0;
			wmass_sc = wmass_max_lyr;
			Fwmass_sc = wmass_sc / thickness;
		}
	}
	else {
		// Should water freeze
		if (temp < 0.0 && Fwmass > 0.0) {
			ifFreeze = true;
			temp_old = 0.0;
		}
		// Should ice thaw
		else if (temp > 0.0 && Fimass > 0.0) {
			ifThaw = true;
			temp_old = 0.0;
		}
	}

	// No phase change, return without doing phase change
	if (!ifThaw && !ifFreeze) {
		temp_old = temp;
		return;
	}

	// Temperature difference from previous phase change
	delta_T = temp - temp_old;

	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy from full soil layer to freeze water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / Lheat;

		// Super cooled water that could freeze
		double Fwmass_freeze = Fwmass - Fwmass_sc;

		// Can't freeze all
		if (Fwmass_freeze > Ffreez) {

			// Update water and ice content
			wmass = (Fwmass - Ffreez) * thickness;
			imass = (Fimass + Ffreez) * thickness;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			wmass = (Fwmass - Fwmass_freeze) * thickness;
			imass = (Fimass + Fwmass_freeze) * thickness;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fwmass_freeze * Lheat;

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update soil temperature
			temp = temp_old - energy_remain / Ch;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from full soil layer to thaw water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = energy / Lheat;

		// Maximum ice that could be turned into liquid form
		double Fimass_thaw = temp > 0.0 ? Fimass : min(Fimass, Fwmass_sc - Fwmass);

		// Can't thaw all
		if (Fimass_thaw > Fthaw) {

			// Update water and ice content
			wmass = (Fwmass + Fthaw) * thickness;
			imass = (Fimass - Fthaw) * thickness;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass = (Fwmass + Fimass_thaw) * thickness;
			imass = (Fimass - Fimass_thaw) * thickness;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fimass_thaw * Lheat;

			// Calculate remaining energy to lower the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update soil temperature in this layer
			temp = temp_old + energy_remain / Ch;
		}
	} // End of Thawing

	// Set new old temperature
	temp_old = temp;

	// Final water
	double wtotal_final = wmass + imass;
	
	// Water balance checks
	if ((fabs(wtotal_initial - wtotal_final)) > EPS) {
		fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n",
			wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
	}

	// Ice despite positive temperature
	if (temp > 0.0 && imass > 0.0) {
		dprintf("ERROR: Phase Change Year %d Day %d Depth %g Temp %g imass %g\n",
			date.year, date.day, depth_mid, temp, imass);
	}

	// No negative wmass
	if (wmass < -EPS) {
		fail("ERROR: Negative wmass at end of phase_change! Year %d Day %d depth_mid %g thickness %g delta_T %g wmass %g imass %g\n",
			date.year, date.day, depth_mid, thickness, delta_T, wmass, imass);
	}

	// No negative imass
	if (imass < -EPS) {
		fail("ERROR: Negative imass at end of phase_change! Year %d Day %d depth_mid %g thickness %g delta_T %g wmass %g imass %g\n",
			date.year, date.day, depth_mid, thickness, delta_T, wmass, imass);
	}
}

/** Water content of soil layer
 *  wcont never reports over field capacity even if there is more water
 */
double Soillayer::wcont() const {
	return wcont(0);
}

/// Water content of soil layer
/** Water content of soil layer
 *  wcont never reports over field capacity even if there is more water
 */
double Soillayer::wcont(double wmass_added) const{
	return min(1.0, max(0.0, wmass_added + wmass - wp) / awc);
}

/// Field capacity of soil layer
/** Field capacity of soil layer
 */
double Soillayer::fc() {
	return awc + wp;
}

/// Maximum water capacity of soil layer
/** Maximum water capacity of soil layer
 */
double Soillayer::wmass_max() {
	if (hydrology == RICHARDS || 
		(*soil).patch.stand.is_highlatitude_peatland_stand() || 
		(*soil).patch.stand.is_true_wetland_stand()) {
		return wsats;
	}
	else {
		return fc();
	}
}

/// Fraction of air in soil layer
/** Fraction of air in soil layer
 */
double Soillayer::fracair() {
	return max(0.0, (wsats - (max(0.0, wmass) + imass)) / thickness);
}

/// Get the water filled pore space of the layer
/** Get the water filled pore space of the layer
 */
double Soillayer::wfps() {
	return max(0.0, min(1.0, (wmass + imass) / wsats));
}

/// C content of soil layer
/** C content of soil layer
 * Could also include/exclude surface litter and woody soil litter
 */
double Soillayer::ccont(grouppooltype group) const {
	double ccont = 0.0;

	ccont += cpool_fast;
	ccont += cpool_slow;

	for (int p = 0; p < NSOMPOOL; p++) {
		if (group == GALL || group == GSURFL) {
			if (p == SURFMETA || p == SURFSTRUCT || p == SURFFWD || p == SURFCWD) {
				ccont += sompool[p].cmass;
			}
		}

		if (group == GALL || group == GHUMUS || group == GHUMUSSOM) {
			if (p == SURFMICRO || p == SURFHUMUS) {
				ccont += sompool[p].cmass;
			}
		}
		
		if (group == GALL || group == GSOM || group == GHUMUSSOM) {
			if (p == SOILMETA || p == SOILSTRUCT || p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
				ccont += sompool[p].cmass;
			}
		}
	}

	return ccont;
}

/// N content of soil layer
/** N content of soil layer
 * Could also include/exclude litter
 */
double Soillayer::ncont(grouppooltype group) const {
	double ncont = NH4_mass + NO3_mass + NO2_mass + NO_mass + N2O_mass + N2_mass;

	for (int p = 0; p < NSOMPOOL; p++) {
		if (group == GALL || group == GSURFL) {
			if (p == SURFMETA || p == SURFSTRUCT || p == SURFFWD || p == SURFCWD) {
				ncont += sompool[p].nmass;
			}
		}

		if (group == GALL || group == GHUMUS || group == GHUMUSSOM) {
			if (p == SURFMICRO || p == SURFHUMUS) {
				ncont += sompool[p].nmass;
			}
		}

		if (group == GALL || group == GSOM || group == GHUMUSSOM) {
			if (p == SOILMETA || p == SOILSTRUCT || p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
				ncont += sompool[p].nmass;
			}
		}
	}

	return ncont;
}

/// Available N content of soil layer
/** Available N content of soil layer
 * Could be type specific. Could also be water limited.
 */
double Soillayer::nmass_avail(n_pref_type pref, bool wmass_limited) {

	double nmass = 0.0;

	// If N transformation is turned off, then all available N should be of the type NH4
	if (!ifntransform) {
		pref = NH4;
		if (NO3_mass > 0.0) {
			NH4_mass += NO3_mass;
			NO3_mass = 0.0;
		}
	}

	if (pref == NO) {
		nmass = NH4_mass + NO3_mass;
	}
	else if (pref == NH4) {
		nmass = NH4_mass;
	}
	else if (pref == NO3) {
		nmass = NO3_mass;
	}

	// return N value with water scalar if water limiting
	return nmass * (wmass_limited ? wcont() : 1.0);
}

/// N substraction from soil layer
/** N substraction from soil layer
 */
void Soillayer::nmass_subtract(double nmass, n_pref_type pref) {

	if (!ifntransform) {
		pref = NH4;
	}
	double residual = 0.0;
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		residual = nmass_tot - nmass;

		// N is available for subtraction
		if (residual >= 0.0) {
			nmass_subtract(nmass * NH4_mass / nmass_tot, NH4);
			nmass_subtract(nmass * NO3_mass / nmass_tot, NO3);
		}
		// Very small negative residual, add missing N as NO3 flux
		else if (negligible(residual, -10) || date.year < freenyears) {
			if (date.year >= freenyears && residual < -1.0e-12) {
				dprintf("ERROR year %d and day %d! NH4 and NO3 mass %f is smaller than N substr %f with the residual %g. Adding residual as NO flux\n",
					date.year, date.day, nmass_tot, nmass, residual);
			}
			NO3_mass = 0.0;
			NH4_mass = 0.0;
			(*soil).patch.fluxes.report_flux(Fluxes::NH3_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NO) than available (residual %g)\n",
				date.year, date.day, (*soil).patch.stand.id, (*soil).patch.id, depth_mid, residual);
		}
	}
	else if (pref == NH4) {
		residual = NH4_mass - nmass;
		// N is available for subtraction
		if (residual >= 0.0) {
			NH4_mass -= nmass;
		}
		// Very small negative residual, add missing N as NO3 flux
		else if (negligible(residual, -10) || date.year < freenyears) {
			if (date.year >= freenyears && residual < -1.0e-12) {
				dprintf("ERROR year %d and day %d! NH4 mass %f is smaller than N substr %f with the residual %g. Adding residual as NH3 flux\n",
					date.year, date.day, NH4_mass, nmass, residual);
			}
			NH4_mass = 0.0;
			(*soil).patch.fluxes.report_flux(Fluxes::NH3_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NH4) than available (residual %g)\n",
				date.year, date.day, (*soil).patch.stand.id, (*soil).patch.id, depth_mid, residual);
		}
	}
	else if (pref == NO3) {
		residual = NO3_mass - nmass;
		// N is available for subtraction
		if (residual >= 0.0) {
			NO3_mass -= nmass;
		}
		// Very small negative residual, add missing N as NO flux
		else if (negligible(residual, -10) || date.year < freenyears) {
			if (date.year >= freenyears && residual < -1.0e-12) {
				dprintf("ERROR year %d and day %d! NO3 mass %f is smaller than N substr %f with the residual %g. Adding residual as NO flux\n",
					date.year, date.day, NO3_mass, nmass, residual);
			}
			NO3_mass = 0.0;
			(*soil).patch.fluxes.report_flux(Fluxes::NO_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NO3) than available (residual %g)\n",
				date.year, date.day, (*soil).patch.stand.id, (*soil).patch.id, depth_mid, residual);
		}
	}
}

/// N increment to soil layer
/** N increment to soil layer
 */
void Soillayer::nmass_inc(double nmass, n_pref_type pref) {

	// If negative increment the subtract nmass
	if (nmass < 0.0) {
		nmass_subtract(fabs(nmass), pref);
	}
	else {
		// If N transformation is turned off, then all available N should be of the type NH4
		if (!ifntransform) {
			pref = NH4;
		}

		if (pref == NO) {
			double nmass_tot = NH4_mass + NO3_mass;
			NH4_mass += nmass * NH4_mass / nmass_tot;
			NO3_mass += nmass * NO3_mass / nmass_tot;
		}
		else if (pref == NH4) {
			NH4_mass += nmass;
		}
		else if (pref == NO3) {
			NO3_mass += nmass;
		}
	}
}

/// Soil layer nitrogen saturation level
/** Soil layer nitrogen saturation level
 * Corresponds to the amount of soil available nitrogen where SOM C:N ratio reach
 * their minimum (nitrogen saturation) [kg N m-3]
 */
double Soillayer::nmass_saturation(double NMASS_SAT) {
	return NMASS_SAT * thickness * M_PER_MM;
}

/// Update soil layer temperature due to mass movement
/** Update soil layer temperature due to mass movement
 *  assuming it replaces air
 */
void Soillayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible((thickness * Ch - add_mass * Cp_air) + add_mass * add_Ch))
		dprintf("ERROR: Year %d Day %d Updating soil temperature\n", date.year, date.day);

	// Update temperature
	temp = (temp * (thickness * Ch - add_mass * Cp_air) + add_temp * add_mass * add_Ch) / ((thickness * Ch - add_mass * Cp_air) + add_mass * add_Ch);

	// If temperature becomes positive (negative) with ice (water) present
	// then update thermal properties and perform a phase change for this layer
	if (date.year >= FIRST_FREEZE_YEAR && ((temp > 0.0 && imass > 0.0) || (temp < 0.0 && wmass > 0.0))) {
		// Update thermal properties of layer
		update_thermal_properties(); 

		// Phase change between water and ice and adjust soil layer temperatures accordingly
		phase_change();
	}
}

/**
 * Help function for soil water matric potential calculation
 *
 * This is an optimisation which improves simulation runtime by ~2.5%.
 *
 * Equivalent to theta^k.
 */
double hydro_pow(double theta, double k) {

	double z = k * log_hydro_lut(max(0.001, theta));

	return exp_hydro_lut(z);
}

/// Updates soil water matric potential
/** Campbell (1974) soil water retention function,
 * which predicts the soil water matric potential
 */
void Soillayer::update_matric_potential() {
	double theta = min(1.0, max(0.001, wmass) / wsats);

	Psi = max(-1.0e8, Psi_max * hydro_pow(theta, -b));

	// Derivative in matric potential
	dPsidw = -b * Psi / theta;
}

/// Updates heat capacity and thermal conductivity
/** Updates heat capacity and thermal conductivity 
 *  of the soil layer
 */
void Soillayer::update_thermal_properties() {

	// Log terms of
	// Water fraction
	double lKfwater;
	// Ice fraction
	double lKfice;
	// Organic fraction
	double lKforg;
	// Peat fraction
	double lKfpeat;
	// Mineral fraction
	double lKfmin;
	// Coarse root fraction
	double lKfcoarse_root;
	// Total solids in soil
	double Stot;
	// Organic water holding capacities at saturation
	double organic_wsats = iforganicsoilproperties ? wsats_org : wsats_min;

	// Depth each quantity occupies
	double combined_depth = thickness;
	double ice = imass;
	double water = wmass;
	double mineral = mineral_volume * (1.0 - wsats_min) + coarse_volume;
	double coarse_root = root_volume;
	double peat, organic;

	if (acrotelm || catotelm) {
		peat = organic_volume * (1.0 - organic_wsats);
		organic = 0.0;
	}
	else {
		peat = 0.0;
		organic = organic_volume * (1.0 - organic_wsats);
	}

	// Soil layer air (mm)
	double air = max(0.0, wsats - (wmass + imass));

	// Update heat capacity (J m-3 K-1)
	Ch = (mineral * Cp_min + organic * Cp_org + ice * Cp_ice + water * Cp_water + peat * Cp_peat + air * Cp_air + coarse_root * Cp_wood) / combined_depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	// Total solids in soil
	Stot = mineral + organic + ice + water + peat + coarse_root;

	if (fabs(1.0 - (Stot + air) / combined_depth) > 1.0e-6)
		fail("Year %d Day %d SID %d PID %d Profile depth %g Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g\n", 
			date.year, date.day, (*soil).patch.stand.id, (*soil).patch.id, depth_mid,
			Stot / combined_depth, air / combined_depth, (Stot + air) / combined_depth);

	// Calculate log terms
	lKforg = (organic / Stot) * Korg_log;
	lKfpeat = (peat / Stot) * Kpeat_log;
	lKfmin = (mineral / Stot) * Kmin_log;
	lKfwater = (water / Stot) * Kwater_log;
	lKfice = (ice / Stot) * Kice_log;
	lKfcoarse_root = (coarse_root / Stot) * Kwood_log;

	// Total log terms
	double lKtot = lKfwater + lKfice + lKforg + lKfmin + lKfpeat + lKfcoarse_root;

	// Update thermal conductivity
	Kh = (air * Kair + Stot * exp_thermal_lut(lKtot)) / combined_depth;
}

/// Update organic matter age of soil layer
/** Update organic matter age of soil layer
 */
void Soillayer::update_socage(double cmass, double age) {

	if (!negligible(cmass)) {

		double ctot = ccont();
		socage = (ctot * socage + cmass * age) / (ctot + cmass);
	}
}

/// Return organic fraction of soil layer
/** Return organic fraction of soil layer
 *  as org / (org + min)
 */
double Soillayer::org_frac() {

	// Calculate soil layer organic carbon masses
	double ccont_org = ccont(GHUMUSSOM);

	// Calculate soil layer organic carbon densities
	double rho_org = calc_organic_density(GHUMUSSOM);

	// Update organic volume (mm)
	double org_volume = ccont_org / rho_org * MM_PER_M;

	return org_volume / (org_volume + mineral_volume + coarse_volume);
}

double Soillayer::psi_mpa() const {
	// 1 mm water head = rho * g * 0.001 m = 9.81 Pa
	return Psi * 9.81e-6;
}

/// Serialize soil layer variables
/** Serialize soil layer variables
 */
void Soillayer::serialize(ArchiveStream& arch) {
	arch & thickness 
		& depth_mid
		& depth_bot
		& wmass
		& imass
		& temp
		& temp_old
		& rho_org
		& awc_org
		& wp_org
		& wsats_org
		& b_org
		& Psi_max_org
		& Kw_max_org
		& acrotelm
		& catotelm
		& Cp_min;

	for (int i = 0; i < NSOMPOOL; i++) {
		arch & sompool[i];
	}

	arch & orgleachfrac
		& NH4_mass
		& NO3_mass
		& NO2_mass
		& N2O_mass
		& NO_mass
		& N2_mass
		& CH4
		& O2
		& cpool_slow
		& cpool_fast
		& depth_mod
		& a

		// Equilsom variables
		& solvesom
		& morgleach_mean
		& mminleach_mean
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean;	
}

void Soillayer::daily_accounting() {
	if (date.day == 0) {
		annual_accounting();
	}
	nvol = 0;
	nloss_gas = 0;
}

void Soillayer::annual_accounting() {
	anvol = 0;
	anloss_gas = 0;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soil member functions
int get_solvesomcent_endyr() {
	return freenyears + 10;
}

/// Initilize soil variables
/** Initilize soil variables
 */
void Soil::initdrivers() {

	depth = 0.0;
	nsoillayers = 0;

	wcont_evap = 0.0;
	NH4_input = 0.0;
	NO3_input = 0.0;
	anfix = 0.0;
	anfix_mean = 0.0;
	anfix_calc = 0.0;
	awcont_upper = 0.0;
	rain_melt = 0.0;
	rain_melt_temp = 0.0;
	max_rain_melt = 0.0;
	wmass_wetland_added = 0.0;
	percolate = false;

	// Peatland hydrology variables
	thaw = 0.0;
	maxthaw = 0.0;
	awtp = 0.0;
	acro_co2 = 934.0;
	wtd = 0.0;
	wtp = 0.0;
	soil_psi = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mthaw[mth] = 0.0;
		mwtp[mth] = 0.0;
		mald_mean[mth] = 0.0;

		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
	}

	for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
		soillayerthicks[lyr] = 0.0;
		wfluxes[lyr] = 0.0;
	}
	wfluxes[NSOILLAYERS_MAX] = 0.0;

	dmoss_wtp_limit = 1.0;
	dgraminoid_wtp_limit = 1.0;

	nacrotelm = 0;
	ncatotelm = 0;

	k_O2 = 0.0;
	k_CO2 = 0.0;
	k_CH4 = 0.0;
	Ceq_O2 = 0.0;
	Ceq_CO2 = 0.0;
	Ceq_CH4 = 0.0;

	slope_grade_tan = 0.0;
	slope_grade_cos = 0.0;

	firstTempCalc = true;

	// Solve LPJ SOM pools
	solvesom_end = static_cast<int>(0.8 * nyear_spinup);
	solvesom_begin = static_cast<int>(0.7 * nyear_spinup);

	// Solve CENTURY SOM pools
	solvesomcent_beginyr = static_cast<int>(std::round(0.5*freenyears));
	solvesomcent_endyr = get_solvesomcent_endyr();
	savesomcent = false;
	NO3_mass = 0;
	NH4_mass = 0;
}

/// Create soil layers
/** Create soil layers with specific soil depths and thickness
 */
void Soil::create_soillayers(double soil_lyr_thicks[NSOILLAYERS_MAX]) {

	double acc_depth = 0.0;
	double depth_mid = 0.0;

	bool acrotelm = false;
	bool catotelm = false;
	nacrotelm = 0;
	ncatotelm = 0;
	thaw = depth;
	maxthaw_10.add(thaw);

	soillayers.reserve(nsoillayers);

	// Loop over number of soil layers to create
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		soillayerthicks[lyr] = soil_lyr_thicks[lyr];
		acc_depth += soil_lyr_thicks[lyr];
		depth_mid += lyr == 0 ? soil_lyr_thicks[lyr] / 2.0 : (soil_lyr_thicks[lyr - 1] + soil_lyr_thicks[lyr]) / 2.0;
		if (patch.stand.is_highlatitude_peatland_stand()) {
			if (acc_depth <= ACROTELM_DEPTH) {
				acrotelm = true;
				catotelm = false;
				nacrotelm++;
			}
			else {
				acrotelm = false;
				catotelm = true;
				ncatotelm++;
			}
		}

		// Add layer
		add_soillayer(soil_lyr_thicks[lyr], depth_mid, acrotelm, catotelm);
	}
}

/// Add a single soil layer
/** Add a single soil layer to the bottom with a specific depth
 */
void Soil::add_soillayer(double soil_layer_thickness, double soil_layer_depth_mid, bool acrotelm, bool catotelm) {

	// Create soil layer directly inside the vector
	soillayers.emplace_back(this, soil_layer_thickness, soil_layer_depth_mid, acrotelm, catotelm);
	
	// Access the new layer via reference
	Soillayer& soillayer = soillayers.back();

	// Initiate fixed soil layer properties
	soillayer.initiate_properties();

	// Update soil layer properties
	double dval = 0.0;
	soillayer.update_properties(dval, dval, true);
}

const Climate& Soil::get_climate() const {
	// All patches within a stand share the same climate
	return patch.stand.get_climate();
}

/// Set slope properties
/** Set slpoe properties
 */
void Soil::set_slope_properties() {
	// Set fixed slope properties
	// To be able to distinguish between and flat
	// wetland (slope 0) and other landcovers (slope as input)
	// this has to belong to the soil
	// Slope in radians

	double slope_rad = slope * PI / 180.0;

	// High latitude wetland are on flat ground. Mineral wetlands are seasonal and still 
	// have a slope when Richards hydrology is used during none wet season
	//if (patch.stand.is_highlatitude_peatland_stand()) {
	//	slope_rad = 0.0;
	//}

	// Slope gradient [m m-1]
	slope_grade_tan = tan(slope_rad);
	slope_grade_cos = cos(slope_rad);

	// Beta_0 as a function of maximum topographical distribution (Ekcic et al. 2019; eqn 3)
	double beta_0 = 1.357209;

	// Microtopography distribution (mm) (Ekcic et al. 2019; eqn 2)
	ssw.slope_sigma = pow(slope + beta_0, -3.0) * MM_PER_M;

	// Constants
	ssw.sigma_sqrt2 = 1.0 / (ssw.slope_sigma * sqrt(2.0));
	ssw.sigma_sqrt2pi = ssw.slope_sigma / sqrt(2.0 * PI);
	ssw.sigma2 = 1.0 / (2.0 * pow(ssw.slope_sigma, 2));

	// Fraction of the area that is inundated threshold below which no 
	// single connected inundated area spans the grid cell
	ssw.wmass_frac_thresh = 0.4;

	// Calculate surface water threshold below which no single connected  
	// inundated area spans the grid cell (Ekcic et al. 2019 and CLM 5.0 code)
	ssw.wmass_thresh = 0.0;

	if (ssw.slope_sigma > 1.0e-6) {

		// Surface water level (mm)
		double d = 0.0;

		// Determine surface water threshold
		for (int i = 0; i < 4; i++) {
			double fd = 0.5 * (1.0 + erf(d * ssw.sigma_sqrt2)) - ssw.wmass_frac_thresh;
			double dfdd = exp(-pow(d, 2) * ssw.sigma2) / (ssw.slope_sigma * sqrt(2.0 * PI));
			d -= fd / dfdd;
		}

		// Surface water threshold below which no single connected  
		// inundated area spans the grid cell (mm)
		ssw.wmass_thresh = 0.5 * d * (1.0 + erf(d / (ssw.slope_sigma * sqrt(2.0)))) +
			ssw.slope_sigma / sqrt(2.0 * PI) * exp(-pow(d, 2) / (2.0 * pow(ssw.slope_sigma, 2)));
	}
}

/// Initiate soil water content
/** Initiate soil water content
 */
void Soil::init_water() {

	double dz[NSOILLAYERS_MAX], vws_liq[NSOILLAYERS_MAX], wflux_root[NSOILLAYERS_MAX], dt_dz[NSOILLAYERS_MAX],
		Kw[NSOILLAYERS_MAX], dKwdw[NSOILLAYERS_MAX], Psi[NSOILLAYERS_MAX], dPsidw[NSOILLAYERS_MAX],
		qin[NSOILLAYERS_MAX], qout[NSOILLAYERS_MAX], dqidw0[NSOILLAYERS_MAX], dqidw1[NSOILLAYERS_MAX],
		dqodw1[NSOILLAYERS_MAX], dqodw2[NSOILLAYERS_MAX], rmx[NSOILLAYERS_MAX], amx[NSOILLAYERS_MAX],
		bmx[NSOILLAYERS_MAX], cmx[NSOILLAYERS_MAX], dwat[NSOILLAYERS_MAX];

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		soillayer.wmass = soillayer.wsats;
	}

	double dtsub = 3600.0;

	for (int i = 0; i < 100000; i++) {

		// Calculate matrix potential
		update_matric_potential();
		 
		// Calculate hydraulic conductivity and diffusivity
		update_hydraulic_properties(nsoillayers);

		// calculate commonly used variables
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			// saturated volumetric water content
			vws_liq[lyr] = soillayer.wsats / soillayer.thickness;
			// time per mm
			dt_dz[lyr] = dtsub / soillayer.thickness;
			Kw[lyr] = soillayer.Kw;
			Psi[lyr] = soillayer.Psi;
			dKwdw[lyr] = soillayer.dKwdw;
			dPsidw[lyr] = soillayer.dPsidw;
			dz[lyr] = soillayer.thickness;
			wflux_root[lyr] = 0.0;
		}

		// Soil moisture fluxes and their derivatives
		compute_moisture_fluxes_and_derivs(0.0, Kw, Psi, dKwdw, dPsidw, qin, qout,
			dqidw0, dqidw1, dqodw1, dqodw2, vws_liq, dz, nsoillayers);

		// RHS of system of equations
		compute_RHS_moisture_form(wflux_root, qin, qout, dt_dz, rmx, nsoillayers);

		// LHS of system of equations
		compute_LHS_moisture_form(dt_dz, dqidw0, dqidw1, dqodw1, dqodw2,
			amx, bmx, cmx, nsoillayers);

		// solve tridiagonal matrix A.X = B
		TridiagonalCol(amx, bmx, cmx, rmx, dwat, nsoillayers);

		// Update the mass of liquid water
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].wmass += dwat[lyr] * soillayers[lyr].thickness;
		}
	}

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wmass = min(soillayers[lyr].wmass, soillayers[lyr].wsats);
	}
}

/// Returns surface litter C content
/** Returns surface litter C content
 */
double Soil::get_surflittercmass() {
	Soillayer& soillayer = soillayers[0];
	return soillayer.sompool[SURFMETA].cmass + soillayer.sompool[SURFSTRUCT].cmass +
		soillayer.sompool[SURFFWD].cmass + soillayer.sompool[SURFCWD].cmass;
}

/// Returns surface litter N content
/** Returns surface litter N content
 */
double Soil::get_surflitternmass() {
	Soillayer& soillayer = soillayers[0];
	return soillayer.sompool[SURFMETA].nmass + soillayer.sompool[SURFSTRUCT].nmass +
		soillayer.sompool[SURFFWD].nmass + soillayer.sompool[SURFCWD].nmass;
}

/// Returns soil C content
/** Returns soil C content
 */
double Soil::ccont(grouppooltype group) const {

	double ccont = 0.0;
	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ccont += soillayers[lyr].ccont(group);
	}
	return ccont;
}

/// Returns soil N content
/** Returns soil N content
 */
double Soil::ncont(grouppooltype group) const {

	double ncont = 0.0;
	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ncont += soillayers[lyr].ncont(group);
	}
	return ncont;
}

/// Returns surface standing water and ice content
/** Returns surface standing water and ice content
 */
double SSW::water_and_ice_content() const {
	return wmass + imass;
}

/// Returns soil water and ice content
/** Returns soil water and ice content
 */
double Soil::water_and_ice_content() const {
	return get_sum(&Soillayer::wmass) + get_sum(&Soillayer::imass);
}

/// Returns the total value for a certain property
/** Returns the total value for a certain property
 */
double Soil::get_sum(double Soillayer::* var_member) const {
	// Variable to return
	double prop = 0.0;

	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		prop += (soillayers[lyr].*var_member);
	}

	// Return vector
	return prop;
}

/// Returns the total value for a certain property
/** Returns the total value for a certain property
 */
double Soil::get_sum(double(Soillayer::* var_arr)[12], int m) const {
	// Variable to return
	double prop = 0.0;

	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		prop += (soillayers[lyr].*var_arr)[m];
	}

	// Return vector
	return prop;
}

/// Returns a vector with soillayer N fractions
/** Returns a vector with soillayer N fractions
 */
void Soil::get_nfrac(n_pref_type pref, double nfrac[NSOILLAYERS_MAX], double nmin, bool wmass_limited) {

	double ntot = nmass_avail(pref, wmass_limited);
	double nfracsum = 0.0;

	if (ntot > nmin) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			nfrac[lyr] = soillayers[lyr].nmass_avail(pref, wmass_limited) / ntot;
			nfracsum += nfrac[lyr];
		}
		nfrac[0] += 1.0 - min(1.0, nfracsum);
	}
	else {
		nfrac[0] = 1.0;
	}
}

/// Calculate soillayer C root fractions
/** Calculate soillayer C root fractions
 */
void Soil::calc_rootfrac(double rootfrac[NSOILLAYERS_MAX]) {
	
	double cmass_root[NSOILLAYERS_MAX] = { 0.0 };

	// Get root mass per soil layer
	calc_rootmass(cmass_root);

	// Calculate total root mass
	double cmass_root_tot = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		cmass_root_tot += cmass_root[lyr];
	}

	if (!negligible(cmass_root_tot)) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			rootfrac[lyr] += cmass_root[lyr] / cmass_root_tot;
		}
	}
	else {
		rootfrac[0] = 1.0;
	}
}

/// Calculateh soillayer C root mass
/** Calculate soillayer C root mass
 */
void Soil::calc_rootmass(double cmass_root[NSOILLAYERS_MAX]) {

	Vegetation& vegetation = patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			cmass_root[lyr] += indiv.cmass_root * patch.pft[indiv.pft.id].rootdist[lyr];
		}

		vegetation.nextobj();

	} // indiviual loop
}

/// Returns available N for soil 
/** Returns available N for soil 
 */
double Soil::nmass_avail(n_pref_type pref, bool wmass_limited) {

	double nmass_avail = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		nmass_avail += soillayers[lyr].nmass_avail(pref, wmass_limited);
	}
	return nmass_avail;
}

/// Extracts N from soillayers
/** Extracts N from soillayers
 *  This should be done more advanced depending on root distribution and competition in each soil layer
 */
void Soil::nmass_avail_uptake(double uptake, n_pref_type pref) {
	double tot_nmass_avail = nmass_avail(pref);

	if (!negligible(tot_nmass_avail)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			Soillayer& soillayer = soillayers[lyr];
			soillayer.nmass_subtract(uptake * soillayer.nmass_avail(pref) / tot_nmass_avail, pref);
		}
	}
}

/// Extracts N from soillayers
/** Extracts N from soillayers
 *  This should be done more advanced depending on root distribution and competition in each soil layer
 */
void Soil::nmass_avail_inc(double inc, n_pref_type pref) {
	double tot_nmass_avail = nmass_avail(pref);

	if (!negligible(tot_nmass_avail)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			Soillayer& soillayer = soillayers[lyr];
			soillayer.nmass_inc(inc * soillayer.nmass_avail(pref) / tot_nmass_avail, pref);
		}
	}
}

/// Calculate water content for this month
/** Calculate water content for this month
 */
void Soil::calc_mwcont() {
	double acc_wcont = 0.0;
	double acc_awc = 0.0;

	// Monthly standing water and ice
	ssw.mwimass[date.month] += ssw.water_and_ice_content() / date.ndaymonth[date.month];

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Monthly wmass
		soillayer.mwmass[date.month] += soillayer.wmass / date.ndaymonth[date.month];

		// Monthly imass
		soillayer.mimass[date.month] += soillayer.imass / date.ndaymonth[date.month];

		// Monthly wcont
		soillayer.mwcont[date.month] += soillayer.wcont() / date.ndaymonth[date.month];

		// Monthly lateral drainage
		soillayer.mdrain[date.month] += soillayer.runoff_drain;


		if (date.islastday) {

			acc_wcont += lyr == 0 ? soillayer.mwcont[date.month] * soillayer.awc / 2.0 :
				(soillayers[lyr - 1].mwcont[date.month] * soillayers[lyr - 1].awc + soillayer.mwcont[date.month] * soillayer.awc) / 2.0;

			acc_awc += lyr == 0 ? soillayer.awc / 2.0 : (soillayers[lyr - 1].awc + soillayer.awc) / 2.0;
			soillayer.mwcontabove = acc_wcont / acc_awc;
		}
	}
}

/// Calculate this months average soil layer temperature
/** Calculate this months average soil layer temperature
 */
void Soil::calc_mtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		soillayer.mtemp[date.month] += soillayer.temp / date.ndaymonth[date.month];
	}

	for (int plyr = 0; plyr < NPADLAYERS; plyr++) {
		Padlayer& pad_lyr = pad.padlayers[plyr];
		pad_lyr.mtemp[date.month] += pad_lyr.temp / date.ndaymonth[date.month];
	}
}

// 
/// Function to get temperature at a specific depth
/** Function to get temperature at a specific depth
 */
double Soil::temp_at_depth(double depth) const {
	
	// Depth above first layer
	if (depth <= soillayers.front().depth_mid) {
		return soillayers.front().temp;
	}

	// Depth below last layer
	if (depth >= soillayers.back().depth_mid) {
		return soillayers.back().temp;
	}

	// Find the two layers surrounding the depth
	for (int i = 0; i < (int)soillayers.size() - 1; i++) {

		double z1 = soillayers[i].depth_mid;
		double z2 = soillayers[i + 1].depth_mid;

		if (depth >= z1 && depth <= z2) {
			double t1 = soillayers[i].temp;
			double t2 = soillayers[i + 1].temp;

			// Linear interpolation
			return t1 + (t2 - t1) * (depth - z1) / (z2 - z1);
		}
	}

	// Should not reach here
	fail("Depth (%g) not within soil depth (%g)!", depth, soillayers.back().depth_bot);
	return 0.0;
}

/// Average temperature down to a depth
/** Average temperature down to a depth
 */
double Soil::gtemp_rootdist(double rootdist[NSOILLAYERS_MAX]) {

	// Calculate root weighted gtemp
	double avr_gtemp = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		avr_gtemp += soillayers[lyr].gtemp * rootdist[lyr];
	}
	return avr_gtemp;
}

/// Performs daily updates to soil layer temperatures using an analytic function
/**
 * Call each simulation day (if hydrology set to GERTEN) following update of daily air temperature
 * prior to canopy exchange and SOM dynamics
 * Calculation of mid soil layer temperature.
 * Soil temperatures are assumed to follow surface temperatures according to an
 * annual sinusoidal cycle with damped oscillation about a common mean, and a
 * temporal lag.
 *
 * For a sinusoidal cycle, soil temperature at depth z and time t from beginning
 * of cycle given by (Carslaw & Jaeger 1959; Eqn 52; Jury et al 1991):
 *
 *   (1) t(z,t) = t_av + a*exp(-z/d)*sin(omega*t - z/d)
 *
 *   where
 *     t_av      = average (base) air/soil temperature
 *     a         = amplitude of air temp fluctuation
 *     exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z,
 *                 relative to surface temp fluctuation
 *     z/d       = oscillation lag in angular units at soil depth z
 *     z         = soil depth
 *     d         = sqrt(2*k/omega), damping depth
 *     k         = soil thermal diffusivity
 *     omega     = 2*PI/tau, angular frequency of oscillation (radians)
 *     tau       = oscillation period (365 days)
 *
 * Here we assume a sinusoidal cycle, but estimate soil temperatures based on
 * a lag (z/d, converted from angular units to days) relative to air temperature
 * and damping ( exp(-z/d) ) of soil temperature amplitude relative to air
 * temperature amplitude. A linear model for change in air temperature with time
 * during the last 31 days is used to estimate 'lagged' air temperature. Soil
 * temperature today is thus given by:
 *
 *   (2) temp_soil = atemp_mean + exp( -alag ) * ( temp_lag - temp_mean )
 *
 *   where
 *     atemp_mean = mean of monthly mean temperatures for the last year (deg C)
 *     alag       = oscillation lag in angular units at depth 0.25 m
 *                  (corresponds to z/d in Eqn 1)
 *     temp_lag   = air temperature 'lag' days ago (estimated from linear model)
 *                  where 'lag' = 'alag' converted from angular units to days
 *
 * Soil thermal diffusivity (k) is sensitive to soil water content and is estimated
 * monthly based on mean daily soil water content for the past month, interpolating
 * between estimates for 0, 15% and 100% AWHC (Van Duin 1963; Jury et al 1991,
 * Fig 5.11).
 */
void Soil::calc_temp_analytic() {

	Climate climate = patch.get_climate();

	// Conversion factor for soil thermal diffusivity from mm2/s to m2/day
	const double DIFFUS_CONV = 0.0864;
	// Corresponds to omega/2 = pi/365 (Eqn 1)	
	const double HALF_OMEGA = 8.607E-3; 
	// Conversion factor for oscillation lag from angular units to days (=365/(2*PI))
	const double LAG_CONV = 58.09;
	// Regression parameters
	double a, b;
	// Soil thermal diffusivity (m2/day)
	double k;
	// Air temperature 'lag' days ago (see above; deg C)
	double temp_lag;
	// Day numbers
	double day[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };

	// First month of simulation, use air temperature for soil temperature
	if ((date.year == 0 || date.year == patch.stand.first_year) && date.month == 0 && !date.islastday) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].temp = climate.temp;
		}
	}
	else {

		double above_k, current_k;

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			if (date.islastday) {

				// Linearly interpolate soil thermal diffusivity given mean
				// soil water content from this layer and above

				// Only top layer influencing soil thermal diffusivity
				if (lyr == 0) {
					if (soillayer.mwcontabove < 0.15)
						current_k = ((soillayer.thermdiff_15 - soillayer.thermdiff_0) / 0.15 *
							soillayer.mwcontabove + soillayer.thermdiff_0) * DIFFUS_CONV;
					else
						current_k = ((soillayer.thermdiff_100 - soillayer.thermdiff_15) / 0.85 *
						(soillayer.mwcontabove - 0.15) + soillayer.thermdiff_15) * DIFFUS_CONV;

					k = current_k;
				}
				else {
					// Current and all above layers influence soil thermal diffusivity

					Soillayer& soillayer_above = soillayers[lyr - 1];

					// Above layer
					if (soillayer_above.mwcontabove < 0.15)
						above_k = ((soillayer_above.thermdiff_15 - soillayer_above.thermdiff_0) / 0.15 *
							soillayer_above.mwcontabove + soillayer_above.thermdiff_0) * DIFFUS_CONV;
					else
						above_k = ((soillayer_above.thermdiff_100 - soillayer_above.thermdiff_15) / 0.85 *
						(soillayer_above.mwcontabove - 0.15) + soillayer_above.thermdiff_15) * DIFFUS_CONV;

					// Current layer
					if (soillayer.mwcontabove < 0.15)
						current_k = ((soillayer.thermdiff_15 - soillayer.thermdiff_0) / 0.15 *
							soillayer.mwcontabove + soillayer.thermdiff_0) * DIFFUS_CONV;
					else
						current_k = ((soillayer.thermdiff_100 - soillayer.thermdiff_15) / 0.85 *
						(soillayer.mwcontabove - 0.15) + soillayer.thermdiff_15) * DIFFUS_CONV;

					k = (k * (soillayer.depth_mid - (soillayer_above.thickness + soillayer.thickness) / 2.0) +
						above_k * soillayer_above.thickness / 2.0 + current_k * soillayer.thickness / 2.0) / soillayer.depth_mid;
				}

				// Calculate parameters alag and exp(-alag) from Eqn 2
				soillayer.alag = soillayer.depth_mid * M_PER_MM / sqrt(k / HALF_OMEGA); // from Eqn 1
				soillayer.exp_alag = exp(-soillayer.alag);
			}

			// Every day, calculate linear model for trend in daily air
			// temperatures for the last 31 days: temp_day = a + b * day

			double buffer[31];
			climate.dtemp_31.to_array(buffer);
			regress(day, buffer, 31, a, b);

			// Calculate soil temperature lag
			temp_lag = a + b * (30.0 - soillayer.alag * LAG_CONV);

			// Calculate analytic solution for soil temperature
			soillayer.temp = climate.atemp_mean + soillayer.exp_alag * (temp_lag - climate.atemp_mean);
		}
	}
}

/// Returns the value of some vector for a certain depth interval
/** Returns the value of some vector for a certain depth interval
 */
double Soil::temp_depth(double depth_start, double depth_end) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil depth
	double depth_soil_end = min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			acc_var += soillayer.temp * layer_depth;
			acc_scalar += layer_depth;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	// If monthly temperature then padding layer can also be used
	if (depth_end > depth) {
		int lyr = 0;

		while (acc_depth < depth_end && lyr < NPADLAYERS) {

			Padlayer& padlyr = pad.padlayers[lyr];

			lyr_start = min(max(acc_depth, depth_start), depth_end);
			lyr_end = min(max(acc_depth + padlyr.thickness, depth_start), depth_end);
			layer_depth = lyr_end - lyr_start;

			if (layer_depth > 0.0) {
				acc_var += padlyr.temp * layer_depth;
				acc_scalar += layer_depth;
			}

			acc_depth += padlyr.thickness;
			lyr++;
		}
	}

	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

/// Returns the month average temperature for a certain depth interval
/** Returns the month average temperature for a certain depth interval
 */
double Soil::temp_depth(double depth_start, double depth_end, int mth) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil depth
	double depth_soil_end = min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			acc_var += soillayer.mtemp[mth] * layer_depth;
			acc_scalar += layer_depth;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	// If monthly temperature then padding layer can also be used
	if (depth_end > depth) {
		int lyr = 0;

		while (acc_depth < depth_end && lyr < NPADLAYERS) {

			Padlayer& padlyr = pad.padlayers[lyr];

			lyr_start = min(max(acc_depth, depth_start), depth_end);
			lyr_end = min(max(acc_depth + padlyr.thickness, depth_start), depth_end);
			layer_depth = lyr_end - lyr_start;

			if (layer_depth > 0.0) {
				acc_var += padlyr.mtemp[mth] * layer_depth;
				acc_scalar += layer_depth;
			}

			acc_depth += padlyr.thickness;
			lyr++;
		}
	}

	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

double Soil::var_depth(double depth_start, double depth_end,
	double Soillayer::* var_member,	double Soillayer::* scal_member) const {

	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;

	int lyr = 0;

	double depth_soil_end = std::min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		const Soillayer& soillayer = soillayers[lyr];

		lyr_start = std::min(std::max(acc_depth, depth_start), depth_soil_end);
		lyr_end = std::min(std::max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {

			const double scalar = (soillayer.*scal_member) * layer_depth / soillayer.thickness;

			acc_var += (soillayer.*var_member) * scalar;
			acc_scalar += scalar;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	// No padding logic needed here; this overload is for scalar members only
	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

double Soil::var_depth(double depth_start, double depth_end,
	double(Soillayer::* var_arr)[12], double Soillayer::* scal_member, int m) const {
	
	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;

	int lyr = 0;

	double depth_soil_end = std::min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		const Soillayer& soillayer = soillayers[lyr];

		lyr_start = std::min(std::max(acc_depth, depth_start), depth_soil_end);
		lyr_end = std::min(std::max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {

			double scalar = (soillayer.*scal_member) * layer_depth / soillayer.thickness;

			acc_var += (soillayer.*var_arr)[m] * scalar;
			acc_scalar += scalar;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

double Soil::var_depth(double depth_start, double depth_end,
	double(Soillayer::* var_arr)[12], int m) const {

	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	double depth_soil_end = std::min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		const Soillayer& soillayer = soillayers[lyr];

		lyr_start = std::min(std::max(acc_depth, depth_start), depth_soil_end);
		lyr_end = std::min(std::max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);

		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0)
			acc_var += (soillayer.*var_arr)[m] * layer_depth / soillayer.thickness;

		acc_depth += soillayer.thickness;
		lyr++;
	}
	return acc_var;
}

double Soil::var_depth(double depth_start, double depth_end, double (Soillayer::*getter)() const) const {

	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	double depth_soil_end = std::min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		const Soillayer& soillayer = soillayers[lyr];

		lyr_start = std::min(std::max(acc_depth, depth_start), depth_soil_end);
		lyr_end = std::min(std::max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);

		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0)
			acc_var += (soillayer.*getter)() * layer_depth / soillayer.thickness;

		acc_depth += soillayer.thickness;
		lyr++;
	}
	return acc_var;
}

double Soil::var_depth(double depth_start, double depth_end,
	double Soillayer::* member) const {

	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	double depth_soil_end = std::min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		const Soillayer& soillayer = soillayers[lyr];

		lyr_start = std::min(std::max(acc_depth, depth_start), depth_soil_end);
		lyr_end = std::min(std::max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);

		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0)
			acc_var += soillayer.*member * layer_depth / soillayer.thickness;

		acc_depth += soillayer.thickness;
		lyr++;
	}
	return acc_var;
}


/// Returns the value of soil C for a certain depth interval
/** Returns the value of soil C for a certain depth interval
 */
double Soil::var_depth_ccont(double depth_start, double depth_end, grouppooltype group) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil depth
	//double depth_soil_end = min(depth, depth_end);
	double depth_soil_end = depth_end;

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			scalar = layer_depth / soillayer.thickness;
			acc_var += soillayer.ccont(group) * scalar;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	return acc_var;
}

/// Returns the value of soil N for a certain depth interval
/** Returns the value of soil N for a certain depth interval
 */
double Soil::var_depth_ncont(double depth_start, double depth_end, grouppooltype group) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil depth
	//double depth_soil_end = min(depth, depth_end);
	double depth_soil_end = depth_end;

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.thickness, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			scalar = layer_depth / soillayer.thickness;
			acc_var += soillayer.ncont(group) * scalar;
		}

		acc_depth += soillayer.thickness;
		lyr++;
	}

	return acc_var;
}

/// Snowfall
/** Snowfall
*  Calculate snowfall density and redistribute snow layers
*/
void Soil::snowfall(double airT, double& prec, double u10) {

	if (negligible(ssw.water_and_ice_content(), -10) || snow.nsnowlayers) {

		// Fresh snow density calculation, following eqn 1 in Vionnet et al (2012)
		// Constant [kg m-3]
		int a = 109;
		// Constant [kg m-3 K-1]
		int b = 6;
		// Constant [kg m-7/2 s-1/2]
		int c = 26;
		// Temperature of the melting point for water
		double Tmelt = 0.0;

		// Snow fall density
		double snowfall_dens = max(rho_SNOW_0, a + b * (airT - Tmelt) + c * pow(u10, 0.5));

		// Depth of the newly fallen snow based on precipitation & density (Oleason 2013, Eq. 6.88)
		double snowfall_depth = prec / (snowfall_dens / rho_H2O);

		// Snow fall heat capacity
		double snowfall_Ch = snowheatcap(prec, 0.0, snowfall_depth, airT);

		// Snowpack depth
		double snowpack_depth = snow.depth();

		// Add snowfall and redistribute snowpack
		snow.redist_layers(snowpack_depth, snowfall_depth, prec, airT, snowfall_Ch);
	}
	else {
		// Update surface standing water
		ssw.update_temp(airT, 0.0, prec);

		// Update standing ice content (mm)
		ssw.imass += prec;

		// Do phase change of standing suface water
		ssw.phase_change();

		// All standing surface water has turned into ice. Create a snow layer
		// consisting of only ice.
		if (negligible(ssw.wmass)) {
			snow.snowlayers.resize(snow.nsnowlayers + 1);
			snow.snowlayers[snow.nsnowlayers] = Snowlayer(ssw.imass, 0.0, rho_ICE, ssw.temp);
			snow.snowlayers[snow.nsnowlayers].thickness = snow.snowlayers[snow.nsnowlayers].imass * rho_H2O / snow.snowlayers[snow.nsnowlayers].dens;
			snow.nsnowlayers++;
			ssw.imass = 0.0;
			ssw.temp = 0.0;
		}
	}

	// Reset precipitation
	prec = 0.0;
}

/// Snow accumulation
/** Snow accumulation
 * Rain or snowfall, and melt if simple snow scheme is used
 */
void Soil::snow_or_rainfall(double airT, double& prec, double u10, double& rain_melt, double& rain_melt_temp) {

	// No precipitation and freezing temperatures, then return from function
	if (prec == 0 && airT < 0.0)
		return;

	double w_begin = ssw.water_and_ice_content() + snow.water_and_ice_content() + prec + rain_melt;

	// Multi layer snow scheme
	if (ifmultilayersnow) {

		// Precipitation
		if (!negligible(prec)) {

			// Snowfall
			if (airT < 0.0) {

				// Add precipitation as snow and update depth, density
				snowfall(airT, prec, u10);
			}
			// Precipitation falls on snow
			else if (airT > 0.0 && !negligible(snow.depth())) {

				// Phase changes and temperature change due to rainwater
				snow.rainwater_melt(airT, prec, rain_melt, rain_melt_temp);
			}
			else {

				// Precipitation falls as rain

				// Update temperature
				rain_melt_temp = (rain_melt * rain_melt_temp + prec * airT) / (rain_melt + prec);

				// Add precipitation
				rain_melt += prec;
			}
		}
	}

	// Simple single snow layer scheme (v4.0)
	else {

		// Daily calculation of snowfall and rainfall from precipitation and snow melt from
		// snow pack; update of snow pack with new snow and snow melt and snow melt

		// Maximum size of snowpack (mm) (S. Sitch, pers. comm. 2001-11-28)
		const double SNOWPACK_MAX = 10000.0;

		// Fixed snow density (kg m-3) JSBACH value from Ekici et al (2015)
		const double SNOWDENSITY = 250.0;

		// Snowing today
		if (airT < 0.0) {

			// Snowfall up to max snowpack size
			double snowfall;

			// Snowpack exist
			if (snow.nsnowlayers) {
				snowfall = min(prec, SNOWPACK_MAX - snow.snowlayers[0].imass);
			}
			else {
				snowfall = min(prec, SNOWPACK_MAX);
			}

			rain_melt += max(0.0, prec - snowfall);

			// If no snow layers
			if (snow.nsnowlayers == 0) {
				snow.snowlayers.resize(snow.nsnowlayers + 1);
				snow.snowlayers[snow.nsnowlayers] = Snowlayer(snowfall, 0.0, SNOWDENSITY, airT);
				snow.nsnowlayers++;
			}
			else {
				snow.snowlayers[0].imass += snowfall;
			}
		}
		// Snow melting today
		else if (airT > 0.0 && snow.nsnowlayers) {

			Snowlayer& snowlayer = snow.snowlayers[0];

			// Snow melt formulation (Choudhury et al 1998)
			double snowmelt = min((1.5 + 0.007 * prec) * (airT - 0.0), snowlayer.imass);

			snowlayer.imass -= snowmelt;

			rain_melt += prec + snowmelt;
		}
		else {

			// Add precipitation
			rain_melt += prec;
		}

		// Update single snowlayer with new snow or melt
		if (snow.nsnowlayers) {

			Snowlayer& snowlayer = snow.snowlayers[0];

			snowlayer.wmass = 0.0;
			snowlayer.thickness = snowlayer.imass / (snowlayer.dens / rho_H2O);

			if (negligible(snowlayer.imass)) {

				// Erase this empty snow layer
				snow.snowlayers.erase(snow.snowlayers.begin());
				snow.nsnowlayers--;
			}
		}
	}

	// Water balance might not meet mass balance in very low productive and saturated sites 
	double wdiff = w_begin - (ssw.water_and_ice_content() + snow.water_and_ice_content() + rain_melt);
	if (wdiff > 0.0) {
		rain_melt += wdiff;
	}
	else if (wdiff < 0.0) {
		rain_melt += wdiff;

		// Still missing water
		if (rain_melt < 0.0) {

			// Extract first from standing water
			if (ssw.water_and_ice_content() > 0.0) {
				double frac = min(1.0, -rain_melt / ssw.water_and_ice_content());
				rain_melt += frac * ssw.water_and_ice_content();
				ssw.wmass *= 1.0 - frac;
				ssw.imass *= 1.0 - frac;
			}

			// Still missing water. Take from snow
			if (rain_melt < 0.0) {
				if (snow.water_and_ice_content() > 0.0) {
					int lyr = snow.nsnowlayers - 1;
					while (rain_melt < 0.0 && lyr >= 0) {

						Snowlayer& snowlayer = snow.snowlayers[lyr];

						snowlayer.wmass += rain_melt;
						rain_melt = 0.0;
						if (snowlayer.wmass < 0.0) {
							rain_melt = snowlayer.wmass;
							snowlayer.wmass = 0.0;
						}

						if (rain_melt < 0.0) {
							snowlayer.imass += rain_melt;
							rain_melt = 0.0;

							if (snowlayer.imass < 0.0) {
								rain_melt = snowlayer.imass;
								snowlayer.imass = 0.0;
								snow.snowlayers.erase(snow.snowlayers.begin() + lyr);
								snow.nsnowlayers = static_cast<int>(std::round(snow.snowlayers.size()));
							}
							else {
								snowlayer.dens = snowlayer.imass / snowlayer.thickness * rho_H2O;
							}
						}
						lyr--;
					}
				}
				else {
					rain_melt = 0.0;
				}
			}
		}
	}

	double w_end = ssw.water_and_ice_content() + snow.water_and_ice_content() + rain_melt;
	if (fabs(w_begin - w_end) > 1.0e-10) {
		dprintf("Year %d Day %d ERROR: Snow snow_or_rainfall diff %g\n", date.year, date.day, w_begin - w_end);
	}
}

/// Saturation vapor pressure (Pa)
/** Computes the saturation vapor pressure of water in air as a function of
 * temperature using the Tetens approximation
 */
double esat(double temp) {
	return 611.2 * exp(17.67 * temp / (temp + 243.5));
}

/// Specific humidity from vapor pressure (kg/kg)
/** Converts vapor pressure to specific humidity using total air pressure
 */
double spec_hum(double e, double air_press) {
	return 0.622 * e / (air_press - 0.378 * e);
}

/// Routine that calculate surface evaporation, infiltration and runoff
/** Routine that calculate surface evaporation, infiltration and runoff
 */
void Soil::compute_evap_infilt_and_runoff(double dtsub, double wflux_rain_melt, double depth_infilt, 
	double imass_infilt, double wsats_infilt, double imp_ice, double Kw_max, double fevap, double& wflux_evap, 
	double& wflux_runoff_surf, double& wflux_runoff_ssw,
	double& wflux_ssw, double& wflux_soil, double& wflux_ssw_drain) {

	// Reset outputs
	
	// Total evaporation (mm s-1)
	wflux_evap = 0.0;
	// Total surface runoff (mm s-1)
	wflux_runoff_surf = wflux_runoff_ssw = 0.0;
	// Change in surface standing water (mm s-1)
	wflux_ssw = 0.0;
	// Total soil infiltration (mm s-1)
	wflux_soil = 0.0;
	// Water drained from the surface standing water into the soil (mm s-1) 
	wflux_ssw_drain = 0.0;

	// Constants

	// Minimum water considerred (equals to a 0.1 mm3)
	const double EPS = 1.0e-9;
	// Maximum saturated fraction (Wang et al. 2008)
	const double fsat_max = 0.3;
	// Minimum surface standing water before all is removed through runoff
	const double wmass_ssw_min = 0.01;

	// Fractional area with water table at surface (-)
	double fsat;
	// Water directly going to surface runoff (mm s-1)
	double wflux_surf_runoff = 0.0;
	// Water entering the top of the soil soil (mm s-1)
	double wflux_in = 0.0;
	// Runoff from surface standing water (mm s-1)
	double wflux_ssw_runoff = 0.0;
	// Update the interconnected fraction of surface standing water (-)
	double wmass_ssw_con_frac;
	// Evaporation from the soil fraction (mm s-1)
	double wflux_soil_evap;
	// Evaporation from the surface standing water fraction (mm s-1)
	double wflux_ssw_evap;
	// Water that couldn't enter the soil and is added to the surface standing water (mm s-1)
	double wflux_soil_excess;
	// Maximum soil infiltration capacity (mm s-1)
	double wflux_infilt_max;


	// Calculate water content for infiltration layers
	double wmass_infilt = var_depth(0.0, depth_infilt, &Soillayer::wmass);

	// Fraction pore spaces that are filled with water and ice
	double totwcont_infilt = min(1.0, (wmass_infilt + imass_infilt) / wsats_infilt);

	// Fractional area with water table at surface (Wood et al. 1992)
	fsat = 1.0 - pow(1.0 - totwcont_infilt, b_infilt / (1.0 + b_infilt));

	// Calculate surface water fraction and runoff
	if (ssw.wmass > wmass_ssw_min) {

		// Surface water level(m)
		double d = 0.0;

		double dfdd, fd = 1.0;
		int i = 0;
		while (fd > 1.0e-8 && i < 10) {
			fd = 0.5 * d * (1.0 + erf(d * ssw.sigma_sqrt2)) + ssw.sigma_sqrt2pi * exp(-pow(d, 2) * ssw.sigma2) - ssw.wmass * M_PER_MM;
			dfdd = 0.5 * (1.0 + erf(d * ssw.sigma_sqrt2));
			d -= fd / dfdd;
			i++;
		}

		// Update the submerged area fraction using the new d value (Ekici et al. 2019; eqn 1)
		ssw.wmass_frac = 0.5 * (1.0 + erf(d * ssw.sigma_sqrt2));

		// Update the interconnected fraction of surface water
		wmass_ssw_con_frac = 0.0;
		if (ssw.wmass_frac > ssw.wmass_frac_thresh) {
			wmass_ssw_con_frac = pow(ssw.wmass_frac - ssw.wmass_frac_thresh, 0.14);
		}

		// Calculate surface water lateral runoff
		wflux_ssw_runoff = min(1.0, C_dw * max(0.001, slope_grade_tan) * wmass_ssw_con_frac) * max(0.0, ssw.wmass - ssw.wmass_thresh) / dtsub;
	}
	else {
		// Very small amount of surface water -> all lost through surface runoff
		wmass_ssw_con_frac = 0.0;
		ssw.wmass_frac = 0.0;
		wflux_ssw_runoff = ssw.wmass / dtsub;
	}

	// Fast runoff directly to streams
	if (wflux_rain_melt > 0.0) {

		// Surface runoff (mm s-1)
		wflux_surf_runoff = fsat * wflux_rain_melt;

		// Water staying at surface / soil (mm s-1)
		wflux_in = wflux_rain_melt - wflux_surf_runoff;
	}

	// Partition surface inputs between soil and standing water
	wflux_ssw = ssw.wmass_frac * wflux_in;
	wflux_soil = wflux_in - wflux_ssw;

	// Remove surface water runoff
	wflux_ssw -= wflux_ssw_runoff;

	// Evaporation

	// Soil evaporation
	wflux_soil_evap = 0.0;
	// Standing water evaporation
	wflux_ssw_evap = 0.0;

	// Evaporation only occurs if snow depth is below 10 mm and 
	// there is a fraction of ground which isn't coverred by vegetation
	if (snow.depth() < 10.0 && !negligible(fevap)) {

		// Evaporation is calculated on a sub daily time step as water content 
		// in top soil layer changes with each time step

		// Evaporation for high-latitude peatland soils
		if (patch.stand.is_highlatitude_peatland_stand()) {

			// Evaporation from peat soil surface (Wania et al 2009)
			if (wtd < SOILDEPTH_EVAP) {
				wflux_soil_evap = fevap * get_climate().eet * PRIESTLEY_TAYLOR * (0.99 / (1.0 + exp(-1.0 * (-wtd + 98.7) / 22.6) + 0.02));
				wflux_soil_evap *= (1.0 - ssw.wmass_frac) / SECS_PER_DAY;
			}
		}
		else {
			// Evaporation for mineral soils (mm s-1)

			// Top soil layer
			Soillayer& surflayer = soillayers[0];

			// Soil moisture stress factor on evaporation
			if (hydrology != RICHARDS) {
				// Standard wcont for the top soil layer
				wcont_evap = surflayer.wcont();
			}
			else {
				// For Richards hydrology we allow evaporation to dry soil below wp
				double low_lim = surflayer.wp / 2.0;
				wcont_evap = min(1.0, max(0.0, surflayer.wmass - low_lim) / (surflayer.fc() - low_lim));
			}

			// Mineral soil evaporation (mm s-1)
			// Note: Wind speed doesn't affect soil evaporation
			if (ifdavehydrology) {
				// Only multiply by wcont once: can't find justification for
				// squared term. This first appears in Rost 2008, Eq A5, but no
				// squared term in Gerten 2004.
				wflux_soil_evap = get_climate().eet * PRIESTLEY_TAYLOR * wcont_evap * fevap *
					(1.0 - ssw.wmass_frac) / SECS_PER_DAY;
			} else {
				wflux_soil_evap = get_climate().eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap *
					(1.0 - ssw.wmass_frac) / SECS_PER_DAY;
			}
		}

		// Evaporation for standing water part of the soil (mm s-1)

		// Calculate evaporation only if there is a standing water fraction
		if (!negligible(ssw.wmass_frac)) {

			// Evaporation from standing surface water uses a bulk aerodynamic formulation
			//
			//   E = rho_air * CE * u10 * (q_water - q_air)
			//
			// The water surface is assumed saturated, so humidity at the surface is q_water
			// Evaporation is driven by the humidity gradient (q_water - q_air)
			// Wind enhances turbulent transport of vapor away from the surface
			// Negative fluxes (q_air > q_water) are suppressed using max(0.0, ...),
			// so condensation (dew formation) is not represented

			Climate climate = get_climate();

			// Air pressure
			const double atm_press = 1.0e5;

			// Specific gas constant for dry air (J kg^-1 K^-1)
			const double Rd = 287.05;

			// Bulk transfer (exchange) coefficient for water vapor over open water (-)
			// Typical value for 10 m wind: ~1.1e-3 to 1.5e-3
			const double CE = 1.3e-3;

			// Saturation vapor pressure at air temperature (Pa)
			double es_air = esat(climate.temp);

			// Actual vapor pressure of air (Pa)
			double e = climate.relhum * es_air;

			// Specific humidity of air (kg/kg), derived from vapor pressure and air pressure
			double q_air = spec_hum(e, atm_press);

			// Saturation vapor pressure at water surface temperature (Pa)
			double es_water = esat(ssw.temp);

			// Saturation specific humidity at the water surface (kg/kg)
			// Air in contact with water is assumed saturated
			double q_water = spec_hum(es_water, atm_press);

			// Air density (kg m^-3) using ideal gas law: rho = P / (R_d * T)
			double rho_air = atm_press / (Rd * (climate.temp + K2degC));

			// Evaporation from standing water scaled by fractional coverage (mm s-1)
			// No condensation allowed
			wflux_ssw_evap = max(0.0, rho_air * CE * climate.u10 * (q_water - q_air)) * fevap * ssw.wmass_frac;
		}
	}

	// Remove evaporation
	wflux_ssw -= wflux_ssw_evap;
	wflux_soil -= wflux_soil_evap;

	// Determine the maximum infiltration rate (mm s-1)
	wflux_infilt_max = (1.0 - fsat) * imp_ice * Kw_max * slope_grade_cos;

	// Excess water that can't infiltrate (mm s-1)
	wflux_soil_excess = max(0.0, wflux_soil - (1.0 - ssw.wmass_frac) * wflux_infilt_max);

	// Soil infiltration and surface water run-on
	wflux_ssw += wflux_soil_excess;
	wflux_soil -= wflux_soil_excess;

	// Temporary storage prior to calculating bottom drainage from surface water
	double wmass_ssw_tmp = ssw.wmass + wflux_ssw * dtsub;

	// Negative surface water -> reduce soil infiltration
	if (wmass_ssw_tmp < 0.0) {
		wflux_soil += wmass_ssw_tmp / dtsub;
		wflux_ssw = -ssw.wmass / dtsub;
		wflux_ssw_drain = 0.0;
	}
	else {
		// Drainage from surface water to soil 
		wflux_ssw_drain = min(ssw.wmass_frac * wflux_infilt_max, wmass_ssw_tmp / dtsub);
	}

	// Reduce surface standing water with drainage (mm s-1)
	wflux_ssw -= wflux_ssw_drain;

	// Total evaporation (mm s-1)
	wflux_evap = wflux_soil_evap + wflux_ssw_evap;

	// Total surface runoff this timestep (mm s-1)
	wflux_runoff_surf = wflux_surf_runoff;
	wflux_runoff_ssw = wflux_ssw_runoff;

	double wbudget = wflux_evap + wflux_soil + wflux_ssw_drain + wflux_runoff_surf + wflux_runoff_ssw + wflux_ssw - wflux_rain_melt;

	if (fabs(wbudget) > EPS) {
		dprintf("Year %d Day %d ERROR: Infiltration water budget off: %g\n", 
			date.year, date.day, wbudget);
	}
}

/// Update soil matrix potential for the soil 
/** Update soil matrix potential for the soil 
 */
void Soil::update_matric_potential() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_matric_potential();
	}
}

 /// Updates hydraulic conductivity for the soil 
 /** Updates hydraulic conductivity for the soil 
  * Relate the hydraulic conductivity to the soil water 
  * content following Campbell (1974)
  */
void Soil::update_hydraulic_properties(int nsoillayers_sim) {

	bool full = false;
	int full_lyr = nsoillayers_sim - 1;

	for (int lyr = nsoillayers_sim - 1; lyr >= 0; lyr--) {

		// Layer below
		int nlyr = lyr + 1 < nsoillayers_sim - 1 ? lyr + 1 : lyr;

		Soillayer& soillayer = soillayers[lyr];
		Soillayer& soillayer_below = soillayers[nlyr];

		// Water and ice theta for the intersect
		double theta_w = max(0.001, min(1.0, 0.5 * (soillayer.wmass / soillayer.wsats +
			soillayer_below.wmass / soillayer_below.wsats)));
		double theta_i = max(0.0, min(1.0, 0.5 * (soillayer.imass / soillayer.wsats +
			soillayer_below.imass / soillayer_below.wsats)));

		// Ice impedance factor for the intersect (Swenson et al. 2012)
		double imp_ice = 1.0;
		if (theta_i > 0.0) {
			imp_ice = ice_imp(theta_i);
		}

		// Hydraulic conductivity (Campbell 1974; mm s-1)
		if (theta_w >= 1.0) {
			soillayer.Kw = soillayer.Kw_max;
		}
		else if (theta_w < 0.01) {
			soillayer.Kw = 0.0;
		}
		else {
			soillayer.Kw = imp_ice * soillayer.Kw_max * hydro_pow(theta_w, 2.0 * soillayer.b + 3.0);
		}

		if (soillayer.Kw < 1.0e-20) {
			soillayer.Kw = 0.0;
		}

		// Derivative in hydraulic conductivity w.r.t relative saturation
		soillayer.dKwdw = imp_ice * (2.0 * soillayer.b + 3.0) * soillayer.Kw / theta_w;
	}
}

/// Update soil thermal capacity and conductivity for the soil 
/** Update soil thermal capacity and conductivity for the soil 
 */
void Soil::update_thermal_properties() {

	// Update heat capacity and conductivity for each snow layer
	snow.update_thermal_properties();

	// Calculate heat capacity and conductivity for surface standing water
	ssw.update_thermal_properties();

	// Update heat capacity and conductivity for each soil layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		soillayers[lyr].update_thermal_properties();
	}

	// Calculate heat capacity and conductivity for each padding layer
	pad.update_thermal_properties(soillayers[nsoillayers - 1].Kh, soillayers[nsoillayers - 1].Ch);
}

/// Update surface standing waters thermal capacity and conductivity
/** Update surface standing waters thermal capacity and conductivity
 */
void SSW::update_thermal_properties() {

	double depth = wmass + imass;

	// No standing water then return
	if (negligible(depth, -10)) {
		return;
	}

	// Standing water heat capacity (J m-3 K-1)
	Ch = (wmass * Cp_water + imass * Cp_ice) / depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	// Calculate log terms
	double lKfwater = (wmass / depth) * Kwater_log;
	double lKfice = (imass / depth) * Kice_log;

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh = exp_thermal_lut(lKfwater + lKfice);
}

/// Serialize surface standing water variables
/** Serialize surface standing water variables
 */
void SSW::serialize(ArchiveStream& arch) {
	arch& wmass
	    & imass
	    & temp
	    & Kh
	    & Ch
		& slope_sigma
		& sigma_sqrt2
		& sigma_sqrt2pi
		& sigma2;
}

/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
void cnstep_full(double surf_temp, double dt, int nlayers, double temp[NLAYERS_MAX + 1], 
	double dz[NLAYERS_MAX + 1], double ki[NLAYERS_MAX + 1], double ci[NLAYERS_MAX + 1]) {

	// DESCRIPTION
	//     Crank-Nicholson timestepper for temperature diffusion equation.
	//     This routine performs a single timestep (of length dt) of the heat
	//     diffusion equation
	//
	//                        6 T    6  /      6 T \
	//                   c(z) --- = --- | k(z) --- |
	//                        6 t   6 z \      6 z /
	//
	//     where here, T(z, t) is the temperature, the '6's represent partial
	//     differentiation, and c(z)/k(z) is the (depth-dependent) heat capacities and thermal conductivities
	//     The boundary conditions are a prescribed surface
	//     temperature, and no heat flow at the bottom of the solution domain.
	//
	//     INPUT PARAMETERS:
	//
	//      surf_temp    The surface temperature boundary condition.
	//
	//      dt           The timestep.
	//
	//      nlayers      Total number of layers represented in temp, dz, d, k and
	//                   c arrays.
	//
	//      temp         The temperatures in each layer at the start of the
	//                   timestep.
	//
	//      dz           The thickness of each layer.
	//
	//      d            The thermal diffusivities in each layer.
	//
	//      k            The thermal conductivities in each layer.
	//
	//      c            The heat capacities in each layer.
	//
	//     NOTE ON UNITS: no units are specified here for any of these
	//     inputs; the only constraint on units is that the layer
	//     thicknesses, the diffusion constant and the timestep should be
	//     measured in consistent units. For instance, if the layer
	//     thicknesses are measured in metres, and the timestep in days, {
	//     the diffusion constants should have units of m2 / day.
	//
	//     OUTPUT PARAMETERS
	//
	//      temp         The temperatures in each layer at the end of the
	//                   timestep. 

	// As cnstep, but now we feed in ki and ci

		// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	double diag[NLAYERS_MAX + 1] = { 0.0 };
	double left[NLAYERS_MAX + 1] = { 0.0 };
	double right[NLAYERS_MAX + 1] = { 0.0 };
	double rhs[NLAYERS_MAX + 1] = { 0.0 };
	double solution[NLAYERS_MAX + 1] = { 0.0 };

	// BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE 
	for (int lyr = 0; lyr < nlayers; lyr++) {

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.

		// D+
		if (lyr == nlayers - 1) {
			dplus = 0.0;
		}
		else if (lyr >= nlayers - NPADLAYERS) {
			dplus = ki[lyr] / ci[lyr];
		}
		else {
			dplus = ((ki[lyr] * dz[lyr] + ki[lyr + 1] * dz[lyr + 1]) / (dz[lyr] + dz[lyr + 1])) / ci[lyr];
		}

		// D-
		if (lyr == 0 || lyr >= nlayers - NPADLAYERS) {
			dminus = ki[lyr] / ci[lyr];
		}
		else {
			dminus = ((ki[lyr] * dz[lyr] + ki[lyr - 1] * dz[lyr - 1]) / (dz[lyr] + dz[lyr - 1])) / ci[lyr];
		}

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// and end cases.

		dzhere = dz[lyr];
		there = temp[lyr];

		// --- PLUS ---
		if (lyr < nlayers - 1) {
			// all soil layers apart from the bottom padding layer
			dzplus = dz[lyr + 1];
			tplus = temp[lyr + 1];
		}
		else {
			// bottom padding layer
			dzplus = dz[lyr];
			tplus = temp[lyr];
		}

		// --- MINUS ---
		if (lyr == 0) {
			// top soil layer
			dzminus = dz[lyr];
			tminus = temp[lyr];
		}
		else {
			// all layers except top soil layer
			dzminus = dz[lyr - 1];
			tminus = temp[lyr - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lyr == 0) {
			diag[lyr] = 1.0;
		}
		else {
			diag[lyr] = 1.0 + Cplus + Cp_minus;
		}

		// LEFT
		if (lyr == 0) {
			left[lyr] = 0.0;
		}
		else {
			left[lyr] = -Cp_minus;
		}

		// RIGHT
		if (lyr == 0 || lyr == nlayers - 1) {
			right[lyr] = 0.0;
		}
		else {
			right[lyr] = -Cplus;
		}

		// RHS
		// Calculate right hand side vector values.
		if (lyr == 0) {
			rhs[lyr] = surf_temp;
		}
		else {
			rhs[lyr] = (1.0 - Cplus - Cp_minus) * there +
				Cplus * tplus + Cp_minus * tminus;
		}
	}

	//   SOLVE TRIDIAGONAL SYSTEM

	tridiag(nlayers, left, diag, right, rhs, solution);

	// Update temperature array
	for (int lyr = 0; lyr < nlayers; lyr++) {
		if (verysmall(solution[lyr])) {
			solution[lyr] = 0.0;
		}
		temp[lyr] = solution[lyr];
	}
}

void cnstep(double surf_temp, double dt, int nlayers, double temp[NLAYERS_MAX + 1],
	double dz[NLAYERS_MAX + 1], double di[NLAYERS_MAX + 1]) {
	
	// DESCRIPTION
	//     Crank-Nicholson timestepper for temperature diffusion equation.
	//     This routine performs a single timestep (of length dt) of the heat
	//     diffusion equation
	//
	//                        6 T    6  /      6 T \
	//                   c(z) --- = --- | k(z) --- |
	//                        6 t   6 z \      6 z /
	//
	//     where here, T(z, t) is the temperature, the '6's represent partial
	//     differentiation, and c(z)/k(z) is the (depth-dependent) heat capacities and thermal conductivities
	//     The boundary conditions are a prescribed surface
	//     temperature, and no heat flow at the bottom of the solution domain.
	//
	//     INPUT PARAMETERS:
	//
	//      surf_temp    The surface temperature boundary condition.
	//
	//      dt           The timestep.
	//
	//      nlayers      Total number of layers represented in temp, dz, d, k and
	//                   c arrays.
	//
	//      temp         The temperatures in each layer at the start of the
	//                   timestep.
	//
	//      dz           The thickness of each layer.
	//
	//      d            The thermal diffusivities in each layer.
	//
	//      k            The thermal conductivities in each layer.
	//
	//      c            The heat capacities in each layer.
	//
	//     NOTE ON UNITS: no units are specified here for any of these
	//     inputs; the only constraint on units is that the layer
	//     thicknesses, the diffusion constant and the timestep should be
	//     measured in consistent units. For instance, if the layer
	//     thicknesses are measured in metres, and the timestep in days, {
	//     the diffusion constants should have units of m2 / day.
	//
	//     OUTPUT PARAMETERS
	//
	//      temp         The temperatures in each layer at the end of the
	//                   timestep. 

	// As cnstep, but now we feed in ki and ci

		// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	double diag[NLAYERS_MAX + 1] = { 0.0 };
	double left[NLAYERS_MAX + 1] = { 0.0 };
	double right[NLAYERS_MAX + 1] = { 0.0 };
	double rhs[NLAYERS_MAX + 1] = { 0.0 };
	double solution[NLAYERS_MAX + 1] = { 0.0 };

	// BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE 
	for (int lyr = 0; lyr < nlayers; lyr++) {

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.

		// D+
		if (lyr == nlayers - 1) {
			dplus = 0.0;
		}
		else if (lyr >= nlayers - NPADLAYERS) {
			dplus = di[lyr];
		}
		else {
			dplus = 0.5 * (di[lyr] + di[lyr + 1]);
		}

		// D-
		if (lyr == 0 || lyr >= nlayers - NPADLAYERS) {
			dminus = di[lyr];
		}
		else {
			dminus = 0.5 * (di[lyr] + di[lyr - 1]);
		}

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// and end cases.

		dzhere = dz[lyr];
		there = temp[lyr];

		// --- PLUS ---
		if (lyr < nlayers - 1) {
			// all soil layers apart from the bottom padding layer
			dzplus = dz[lyr + 1];
			tplus = temp[lyr + 1];
		}
		else {
			// bottom padding layer
			dzplus = dz[lyr];
			tplus = temp[lyr];
		}

		// --- MINUS ---
		if (lyr == 0) {
			// top soil layer
			dzminus = dz[lyr];
			tminus = temp[lyr];
		}
		else {
			// all layers except top soil layer
			dzminus = dz[lyr - 1];
			tminus = temp[lyr - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lyr == 0) {
			diag[lyr] = 1.0;
		}
		else {
			diag[lyr] = 1.0 + Cplus + Cp_minus;
		}

		// LEFT
		if (lyr == 0) {
			left[lyr] = 0.0;
		}
		else {
			left[lyr] = -Cp_minus;
		}

		// RIGHT
		if (lyr == 0 || lyr == nlayers - 1) {
			right[lyr] = 0.0;
		}
		else {
			right[lyr] = -Cplus;
		}

		// RHS
		// Calculate right hand side vector values.
		if (lyr == 0) {
			rhs[lyr] = surf_temp;
		}
		else {
			rhs[lyr] = (1.0 - Cplus - Cp_minus) * there +
				Cplus * tplus + Cp_minus * tminus;
		}
	}

	//   SOLVE TRIDIAGONAL SYSTEM

	tridiag(nlayers, left, diag, right, rhs, solution);

	// Update temperature array
	for (int lyr = 0; lyr < nlayers; lyr++) {
		if (verysmall(solution[lyr])) {
			solution[lyr] = 0.0;
		}
		temp[lyr] = solution[lyr];
	}
}

void build_diags_rhs(double dt, double f_top, double f_bot, double temp[NLAYERS_MAX], double ki[NLAYERS_MAX], double ci[NLAYERS_MAX],
	double di[NLAYERS_MAX], double dz[NLAYERS_MAX], double z[NLAYERS_MAX], double sink[NLAYERS_MAX], int nlayers,
	double lower[NLAYERS_MAX], double upper[NLAYERS_MAX], double main[NLAYERS_MAX], double rhs[NLAYERS_MAX]) {

	// cap - Ci
	// diff - Ki
	// cond - di
	// psi - temp
	// dz - in m
	// z middle of each layer

	// Upper boundary conditions
	double l_p12 = ki[0] + (ki[1] - ki[0]) / (z[1] - z[0]) * dz[0] / 2.0; // Interpolate lambda_+1/2
	double g_p12 = di[0] + (di[1] - di[0]) / (z[1] - z[0]) * dz[0] / 2.0; // Interpolate gamma_+1/2
	double beta_p = dt / 2.0 * l_p12 / dz[0] / (z[1] - z[0]) / ci[0];

	lower[0] = 1.0e6; // Not used, initialize to known value
	main[0] = 1.0 + beta_p;
	upper[0] = -beta_p;

	rhs[0] = (1.0 - beta_p) * temp[0] + beta_p * temp[1] - dt / ci[0] / dz[0] * (g_p12 - f_top) + sink[0] * dt;

	// Middle terms
	double l_m12, g_m12, beta_m;
	for (int i = 1; i < nlayers - 1; i++) {
		l_m12 = l_p12;
		l_p12 = ki[i] + (ki[i + 1] - ki[i]) / (z[i + 1] - z[i]) * dz[i] / 2.0;
		beta_m = beta_p * dz[i - 1] / dz[i] * ci[i - 1] / ci[i];
		beta_p = dt / 2.0 * l_p12 / dz[i] / (z[i + 1] - z[i]) / ci[i];

		g_m12 = g_p12;
		g_p12 = di[i] + (di[i + 1] - di[i]) / (z[i + 1] - z[i]) * dz[i] / 2.0;

		lower[i] = -beta_m;
		upper[i] = -beta_p;
		main[i] = 1.0 + beta_m + beta_p;
		rhs[i] = beta_m * temp[i - 1] + (1.0 - beta_m - beta_p) * temp[i] + beta_p * temp[i + 1]
			- dt / ci[i] / dz[i] * (g_p12 - g_m12) + sink[i] * dt;
	}

	beta_m = beta_p * dz[nlayers - 2] / dz[nlayers - 1] * ci[nlayers - 2] / ci[nlayers - 1];
	g_m12 = g_p12;

	upper[nlayers - 1] = 1.0e6; // Dummy value
	main[nlayers - 1] = 1.0 + beta_m;
	lower[nlayers - 1] = -beta_m;

	rhs[nlayers - 1] = beta_m * temp[nlayers - 2] + (1.0 - beta_m) * temp[nlayers - 1]
		- dt / ci[nlayers - 1] / dz[nlayers - 1] * (f_bot - g_m12)
		+ sink[nlayers - 1] * dt;
}

/// Solve tridiagonal system
/** Solve tridiagonal system
 */
void solve_tridiagonal(double lower[NLAYERS_MAX], double upper[NLAYERS_MAX], double main[NLAYERS_MAX], double rhs[NLAYERS_MAX],
	double psi[NLAYERS_MAX], int nlayers) {


	for (int lyr = 1; lyr < nlayers; lyr++) {
		double w = lower[lyr] / main[lyr - 1];
		main[lyr] -= w * upper[lyr - 1];
		rhs[lyr] -= w * rhs[lyr - 1];
	}

	psi[nlayers - 1] = rhs[nlayers - 1] / main[nlayers - 1];

	for (int lyr = nlayers - 2; lyr >= 0; lyr--) {
		psi[lyr] = (rhs[lyr] - upper[lyr] * psi[lyr + 1]) / main[lyr];
	}
}

/// Belda eqns
void belda(double surf_temp, double dt, int nlayers_all, double temp_all[NLAYERS_MAX + 1],
	double dz_all[NLAYERS_MAX + 1], double ki_all[NLAYERS_MAX + 1], double ci_all[NLAYERS_MAX + 1]) {

	// cap - Ci
	// diff - Ki
	// cond - zeros
	// psi - temp
	// dz - in m
	// z middle of each layer
	// 
	// first place in arrays are air and will not be used

	double temp[NLAYERS_MAX] = { 0.0 };
	double ki[NLAYERS_MAX] = { 0.0 };
	double ci[NLAYERS_MAX] = { 0.0 };
	double empty[NLAYERS_MAX] = { 0.0 };
	double dz[NLAYERS_MAX] = { 0.0 };
	double z[NLAYERS_MAX] = { 0.0 };
	double sink[NLAYERS_MAX] = { 0.0 };
	double lower[NLAYERS_MAX] = { 0.0 };
	double upper[NLAYERS_MAX] = { 0.0 };
	double main[NLAYERS_MAX] = { 0.0 };
	double rhs[NLAYERS_MAX] = { 0.0 };

	for (int lyr = 1; lyr < nlayers_all; lyr++) {
		temp[lyr - 1] = temp_all[lyr];
		ki[lyr - 1] = ki_all[lyr] / SECS_PER_DAY * MM_PER_M; // day to s, mm-1 to m-1
		ci[lyr - 1] = ci_all[lyr] / M3_PER_MM3; // mm-3 to m-3
		dz[lyr - 1] = dz_all[lyr] * M_PER_MM;
	}
	int nlayers = nlayers_all - 1;

	double depth = dz[0] / 2.0;
	z[0] = depth;
	for (int lyr = 1; lyr < nlayers; lyr++) {
		depth += (dz[lyr - 1] + dz[lyr]) / 2.0;
		z[lyr] = depth;
	}

	// Boundary conditions for soil temperature (Shouldn't use air temperature)
	// Upper boundary condition [J m-2 s-1]
	double f_top = -2.0 * ki[0] * (temp[0] - surf_temp) / z[0];

	// Lower boundary condition
	double f_bot = 0.0;

	// Initial energy mass balance
	double init_cum_energy = f_top * SECS_PER_DAY * dt;
	for (int lyr = 0; lyr < nlayers; lyr++) {
		init_cum_energy += (temp[lyr] + K2degC) * ci[lyr] * dz[lyr];
	}

	// Build diagonals and RHS vector
	build_diags_rhs(SECS_PER_DAY * dt, f_top, f_bot, temp, ki, ci, empty, dz, z, sink, nlayers, lower, upper, main, rhs);

	// Solve tridiagonal system
	solve_tridiagonal(lower, upper, main, rhs, temp, nlayers);

	// Final energy mass balance
	double final_cum_energy = 0.0;
	for (int lyr = 0; lyr < nlayers; lyr++) {
		final_cum_energy += (temp[lyr] + K2degC) * ci[lyr] * dz[lyr];
	}

	if (fabs(final_cum_energy - init_cum_energy) > 1.0e-12 * final_cum_energy)
		dprintf("Year %d Day %d TS %g Diff %g Final %g Initial %g\n", date.year, date.day, dt, fabs(final_cum_energy - init_cum_energy), final_cum_energy, init_cum_energy);

	// Update temp_all
	for (int lyr = 1; lyr < nlayers_all; lyr++) {
		temp_all[lyr] = temp[lyr - 1];
	}
}

/// Update surface standing water temperatures due to water addition
/** Update surface standing water temperatures due to water addition
 */
void SSW::update_temp(double temp_add, double wmass_add, double imass_add) {

	// Only update standing water temperature if the added water is positive. 
	// Negative would be removal of water and would not change temperature.
	if (wmass_add > 1.0e-6 || imass_add > 1.0e-6) {
		temp = (temp * (wmass + imass) * Ch + temp_add * wmass_add * Cp_water +
			temp_add * imass_add * Cp_ice) /
			((wmass + imass) * Ch + wmass_add * Cp_water + imass_add * Cp_ice);
	}
}

/// Update soil layer temperatures due to water movement
/** Update soil layer temperatures due to water movement
 */
void Soil::wflow_temp(double temp_infilt) {

	const double EPS = 1.0e-12;

	// Update top soil layer for infiltarted rain or melt water
	if (wfluxes[0] > EPS) {
		soillayers[0].update_temp(temp_infilt, wfluxes[0], Cp_water);
	}
	
	// Update soil layer temperature depending on moving water
	for (int lyr = 1; lyr < nsoillayers; lyr++) {
		if (wfluxes[lyr] > EPS) {
			soillayers[lyr].update_temp(max(0.0, soillayers[lyr - 1].temp), wfluxes[lyr], Cp_water);
		}
		else if (wfluxes[lyr] < -EPS) {
			soillayers[lyr - 1].update_temp(max(0.0, soillayers[lyr].temp), -wfluxes[lyr], Cp_water);
		}
	}
}

/// Update padding layers thermal capacity and conductivity
/** Update padding layers thermal capacity and conductivity
 */
void Padlayer::update_thermal_properties(double Kh_last_lyr, double Ch_last_lyr) {

	// Thermal conductivity [W m-1 K-1 s-1]
	Kh = Kh_last_lyr;
	// Heat capacity [J m-3 K-1]
	Ch = Ch_last_lyr;
}

/// Serialize soil padlayer variables
/** Serialize soil padlayer variables
 */
void Padlayer::serialize(ArchiveStream& arch) {
	arch& temp
		& thickness;
}

/// Set soil padding depths
/** Calculate padding layer depths for Crank-Nicholson bottom boundary
 * condition setup. The padding layer thicknesses are calculated as
 * a geometric sequence with a given total padding thickness and
 * number of layers. This gives a smooth increase in layer depth
 * away from the calculation domain.
 */
void Padding::set_depths(double depth_last_lyr) {

	// Iterative improvements of geometric factor
	double k_guess, k_last;

	// Statement functions for function and derivative for Newton
	// solver for geometric sequence factor.

	// Initial state - don't choose 1.0 as the initial value for
	// k_guess, since it's a trivial solution
	k_guess = 10.0;
	k_last = 0.0;

	// Function, derivative 
	double func, dfunc;
	// Thickness of last active soil layer above the padding layers
	double dz = depth_last_lyr;

	// If bottom soil layer multiplied with number of padding layers are less
	// than total padding layer depth, then calculate each padding layer depth
	if (dz * NPADLAYERS < PAD_DEPTH) {

		// Do Newton-Raphson solve to find geometric sequence factor.
		while (fabs(k_guess - k_last) > 1E-2) {
			k_last = k_guess;

			func = pow(k_guess, (NPADLAYERS + 1)) - k_guess * (1.0 + PAD_DEPTH / dz) + PAD_DEPTH / dz;
			dfunc = (NPADLAYERS + 1) * pow(k_guess, NPADLAYERS) - (1.0 + PAD_DEPTH / dz);
			k_guess = k_guess - func / dfunc;
		}

		double tot_depth = 0.0;
		double lyr_depth;

		// Build geometric sequence of padding layer depths.
		for (int lyr = 0; lyr < NPADLAYERS - 1; lyr++) {
			lyr_depth = dz * pow(k_guess, lyr + 1);
			lyr_depth = floor(dz * pow(k_guess, lyr + 1) / 100.0) * 100.0;
			padlayers[lyr].thickness = lyr_depth;
			tot_depth += lyr_depth;
		}

		padlayers[NPADLAYERS - 1].thickness += PAD_DEPTH - tot_depth;
	}
	// Else each padding layer decrease linear to a minimum of 2m for the bottom padding layer
	else {
		double first_pad = min(16000.0, dz);
		double diff = floor((5.0 * first_pad - PAD_DEPTH) / 10.0 / 100.0) * 100.0;
		int scalar = 0;
		double tot_depth = 0.0;

		for (int lyr = 0; lyr < NPADLAYERS - 1; lyr++) {
			padlayers[lyr].thickness = first_pad - diff * (double)scalar;
			tot_depth += padlayers[lyr].thickness;
			scalar++;
		}

		padlayers[NPADLAYERS - 1].thickness = PAD_DEPTH - tot_depth;
	}
}

/// Update thermal properties for padding layers
/** Update thermal properties for padding layers
 */
void Padding::update_thermal_properties(double Kh_last_lyr, double Ch_last_lyr) {

	// Use bottom soil layer properties for the first 2 padding layers
	const int first_bedrock_padding_layer = 2;

	for (int lyr = 0; lyr < first_bedrock_padding_layer; lyr++) {
		padlayers[lyr].update_thermal_properties(Kh_last_lyr, Ch_last_lyr);
	}
}

/// Serialize soil padding variables
/** Serialize soil padding variables
 */
void Padding::serialize(ArchiveStream& arch) {

	for (int i = 0; i < NPADLAYERS; i++) {
		arch& padlayers[i];
	};
}

/// Update the CO2 concentration in the acrotelm water
/** Update the CO2 concentration in the acrotelm water, after Wania et al. (2009)
 * Varies between atmospheric CO2 level when awtp = -300mm, to PORE_WATER_CO2 when awtp = 0, 
 * and greater if there's standing water.
 * Called once daily for wetlands, after the water table has been updated
 *
 * LPJ-WHy comments from Rita Wania:
 *
 * LPJ-WHy uses a different atmospheric CO2 concentration for mosses
 * as they can access the CO2 dissolved in acrotelm water. If the
 * water table is high, mosses can access all of the acrotelm CO2,
 * but as the water table drops the CO2 concentrations available to
 * mosses is a mixture between atmospheric CO2 and acrotelm CO2.
 * NOTE: for ideal gases, ppmv = not just mimol L-1 but also
 * mimol mol-1.
 * The mean annual wtd of the previous year is used here, as the
 * water table is only calculated after the co2 concentration is
 * needed
 */
void Soil::update_acrotelm_co2(double atmo_co2) {
	acro_co2 = min(PORE_WATER_CO2, PORE_WATER_CO2 + ((PORE_WATER_CO2 - atmo_co2) / ACROTELM_DEPTH) * awtp);
}

/// Updates the water content for each soil layer in peatland soils
/**
* Call each simulation day (if stand is of peat type).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
* Uses peatland hydrology scheme from Granberg et al. (1999) and Wania et al. (2009a, 2009b)
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_peat(const Climate& climate, double fevap) {

	// Total AET (mm)
	double aet_total = 0.0;
	// Transpiration from the acrotelm (mm)
	double aet_acrotelm = 0.0;
	// Runoff from upper soil layer (mm)
	double runoff_surf;
	// Runoff from standing water (mm)
	double runoff_ssw = 0.0;
	// Runoff (drainage) from lower soil layers (mm)
	double runoff_drain;
	// Evaporation from soil surface (mm)
	double evap;
	// Evapotranspiration (mm/day)
	double evapotranspiration;
	// An empirical evapotranspiration parameterisation from Kim & Verma
	bool empirical_evapotranspiration = false;
	// Wetland runon set externally (mm)
	double runon = 0.0;
	// Water added to keep water mass balance (mm)
	double waddon = 0.0;
	// Total fpc for mosses. Not used at the moment
	double fpc_moss_total = 0.0;
	// Total water in acrotelm layers (mm)
	double wtot_acrotelm = var_depth(0.0, ACROTELM_DEPTH, &Soillayer::wmass);
	// Total ice in acrotelm layers (mm)
	double itot_acrotelm = var_depth(0.0, ACROTELM_DEPTH, &Soillayer::imass);
	// Initial water mass in the acrotelm
	double wtot_acrotelm_init = rain_melt;
	// Minimum water mass in the acrotelm
	double wtot_acrotelm_min = 0.0;
	// Maximum water mass in the acrotelm
	double wsats_acrotelm = 0.0;

	// Calculate minimum water in acrotelm
	for (int lyr = 0; lyr < nacrotelm; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		if (lyr == 0) {
			wtot_acrotelm_min += minvtot * soillayer.thickness;
		}
		else {
			wtot_acrotelm_min += soillayer.wp;
		}
		wtot_acrotelm_init += soillayer.wmass;
		wsats_acrotelm += soillayer.wsats;
	}
	
	// Transpiration

	// Retrieve Vegetation for this patch
	Vegetation& vegetation = patch.vegetation;

	// Sum AET for across all vegetation individuals in this patch

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		for (int lyr = 0; lyr < nacrotelm; lyr++) {
			pft_acro_root_frac += patch.pft[indiv.pft.id].pft.input_rootdist[lyr];
		}

		// Fraction of the aet from the acrotelm only.
		aet_acrotelm += indiv.aet * pft_acro_root_frac;

		// Total aet, assumed from acrotelm and catotelm combined 
		aet_total += indiv.aet;

		if (indiv.pft.ismoss())
			fpc_moss_total += indiv.fpc;

		vegetation.nextobj();
	}

	// AET for the acrotelm can't make water go below minimum
	aet_acrotelm = max(0.0, min(wtot_acrotelm_init - wtot_acrotelm_min, aet_acrotelm));

	// Water extracted from catotelm
	waddon += aet_total - aet_acrotelm;

	// Water interaction with the surface
	bool wania_approach = false;
	if (wania_approach) {

		// Evaporation

		// Evaporation from soil surface (Wania et al 2009)
		if (snow.depth() < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
			evap = fevap * climate.eet * PRIESTLEY_TAYLOR * (0.99 / (1.0 + exp(-1.0 * (-wtd + 98.7) / 22.6) + 0.02));
		}
		else {
			evap = 0.0;
		}

		// Evapotranspiration

		double evap_excess = 0.0;

		if (empirical_evapotranspiration) {

			// An empirical evapotranspiration parameterisation from Kim & Verma

			// From LPJ-WHyMe:
			// Calculate transpiration and evaporation from wetlands together.
			// dpet is the equilibrium ET but is the same as the potential ET in
			// Kim and Verma (1996).  (BAD choice of variable names in original
			// LPJ!)  Hence, we can use dpet as approximation for potential
			// evapotranspiration and scale it by the water table depth.  daept
			// is allowed to happen even if the wtp is at -300mm because plants
			// will be able to extract water from the unsaturated acrotelm or the
			// saturated catotelm.  This however makes it impossible to close the
			// hydrological cycle in LPJ-WHy.
			// Kim and Verma (1996): ET = PET * (1.02 + 0.00075 * WTP)

			// Note 1. that this is now the TOTAL evapotranspiration, replacing aet_total

			// Note 2: An alternative parameterisation was given by Yurova et al. (2006)
			// evapotranspiration = eet * exp[0.02 * wtp_cm], where // wtp_cm [100,-30]
			// This would give almost identical results to the Kim & Verma formula below 
			// when the water table is at the surface

			if (snow.depth() < 10.0) { // i.e. evap only if snow depth < 1cm
				evapotranspiration = climate.eet * (1.02 - 0.00075 * wtd) * fevap; // wtd [-100,+300]
			}
			else {
				evapotranspiration = 0.0;
			}

			// Evapotranspiration can't exceed available water in the acrotelm
			evapotranspiration = min(wtot_acrotelm, evapotranspiration);
		}
		else {
			// This is JUST evapotranspiration from the acrotelm
			evapotranspiration = aet_acrotelm + evap;
		}

		// Drainage - could possibly be read in for site-specific studies
		runoff_drain = 0;

		// Check so evapotranspiration doesn't exceed available water in the acrotelm
		// If so, then decrease evap accordingly (if negative, water will be added as runon) 
		if (evapotranspiration + runoff_drain > wtot_acrotelm + rain_melt) {
			evap_excess += evapotranspiration + runoff_drain - (wtot_acrotelm + rain_melt);
		}

		// Add water that has been taken above availability
		waddon += evap_excess;

		// Update acrotelm water mass

		// Water fluxes in and out of the acrotelm 
		wtot_acrotelm += rain_melt + evap_excess - evapotranspiration - runoff_drain;
		rain_melt = 0.0;

		// Water runoff (Granberg et al. 1999 eqn 10, table 1)
		double ideal_runoff = exp(-0.01 * wtd);
		if (wtot_acrotelm > 0.0 && soillayers[0].imass * 0.7 < soillayers[0].wmass) {
			runoff_surf = min(ideal_runoff, wtot_acrotelm);
		}
		else {
			runoff_surf = 0.0;
		}
		runoff_ssw = 0.0;

		// Water runon
		// External run on or off - but only when there's flowing water
		if (runoff_surf > 0.0) {

			// e.g. 3 mm/day for bog-like/hummock conditions
			if (wetland_runon < 0)
				runoff_surf -= wetland_runon;

			// 3-5 mm/day for fen-like conditions
			if (wetland_runon > 0)
				runon += wetland_runon;
		}

		// Save initial water mass in the acrotelm
		wtot_acrotelm_init = wtot_acrotelm;

		// Update acrotelm water mass
		wtot_acrotelm += runon - runoff_surf + ssw.wmass;

		// If negative total water, then all initial water has gone to runoff and water is set to zero
		if (wtot_acrotelm < 0.0) {
			wtot_acrotelm = 0.0;
			runoff_surf = wtot_acrotelm_init;
			runon = 0.0;
			ssw.wmass = 0.0;
		}
	}
	else {
		// Water fluxes after calculating the interaction with the surface (mm s-1)
		double wflux_evap, wflux_runoff_surf, wflux_runoff_ssw, wflux_surf, wflux_infilt, wflux_ssw_drain;

		// Timeset (sec)
		double dtsub = SECS_PER_DAY;

		// Water flux into the surfcae (mm s-1)
		double wflux_rain_melt = rain_melt / dtsub;

		// Total saturation level for acrotelm layers (mm)
		double wsats_acrotelm = var_depth(0.0, ACROTELM_DEPTH, &Soillayer::wsats);

		// Calculate ice impedance factor
		double imp_ice = ice_imp(itot_acrotelm / wsats_acrotelm);

		// Average hydraulic conductivity at saturation for infiltation depth (top layer)
		double Kw_max = var_depth(0.0, ACROTELM_DEPTH, &Soillayer::Kw_max, &Soillayer::thickness);

		// Calculate infiltration and surface runoff
		compute_evap_infilt_and_runoff(dtsub, wflux_rain_melt, ACROTELM_DEPTH, itot_acrotelm,
			wsats_acrotelm, imp_ice, Kw_max, fevap, wflux_evap, wflux_runoff_surf, wflux_runoff_ssw,
			wflux_surf, wflux_infilt, wflux_ssw_drain);

		// Water added to standing water
		double wmass_ssw_added = wflux_surf * dtsub;

		// Update standing water temperature
		ssw.update_temp(rain_melt_temp, wmass_ssw_added);

		// Update standing water content (mm)
		ssw.wmass += wmass_ssw_added;

		// Accumulate evaporation
		evap = wflux_evap * dtsub;

		// Accumulate surface runoff
		runoff_surf = wflux_runoff_surf * dtsub;
		runoff_ssw = wflux_runoff_ssw * dtsub;

		// Water fluxes to the acrotelm 
		wtot_acrotelm += (wflux_infilt + wflux_ssw_drain) * dtsub - aet_acrotelm;

		// If negative total water, then take water from standing, runoff, and evap until it is zero 
		if (wtot_acrotelm < 0.0) {
			dprintf("Year %d Day %d ERROR: Peat surface water balance negative (%g)\n", 
				date.year, date.day, wtot_acrotelm);
		}

		rain_melt = 0.0;
	}

	// Very low surface standing water -> send to runoff
	if (negligible(ssw.water_and_ice_content(), -10)) {
		runoff_ssw += ssw.water_and_ice_content();
		ssw.wmass = 0.0;
		ssw.imass = 0.0;
	}

	// Update water table depth

	// LPJ-WHy comments:
	// If total water volume is greater than the porosity then the
	// WTP is above the surface. The maximum height of standing water
	// is limited by maxh. If the watertable is above that height, it
	// will run off.  Snow melt season: 1. If there is a standing
	// water layer, which is obviously frozen in winter, then the
	// snow melt will run off and only when all the ice of the stand
	// layer is thawed, the water will be added to the soil
	// layers. 2. If there is no standing water layer, then the
	// snowmelt water...should penetrate into the soil.
	if (wtot_acrotelm + itot_acrotelm > ACROTELM_DEPTH * acrotelm_por) {
	//if (wtot_acrotelm + itot_acrotelm > wsats_acrotelm) { // DW_COMMENT

		// Saturated acrotelm

		wtp = wtot_acrotelm + itot_acrotelm - ACROTELM_DEPTH * acrotelm_por;
		//wtp = wtot_acrotelm + itot_acrotelm - wsats_acrotelm;
		
		// WTP > 0 above the surface
		if (wtp > maxh) {
			// Update runoff from the surface
			if (wania_approach) {
				runoff_surf += wtp - maxh;
			}
			else {
				double wmass_added = wtp - maxh;

				// Update standing water temperature // DW_COMMENT don't know where the excess water comes from???
				//ssw.update_temp(soillayers[0].temp, wmass_added);

				// Update standing water content (mm)
				ssw.wmass += wmass_added;
			}
			wtp = maxh; // DW_COMMENT THIS IS WRONG!!!!
		}

		// WTD defined as < 0 above the surface
		wtd = -wtp;

		wtot_acrotelm = ACROTELM_DEPTH * acrotelm_por - itot_acrotelm;
		//wtot_acrotelm = wsats_acrotelm - itot_acrotelm;
	}
	else {
		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * az));
		//wtd = sqrt(3.0 * (wsats_acrotelm - (wtot_acrotelm + itot_acrotelm)) / (2.0 * az));

		if (wtd > zmin) {
			wtd = 3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * (acrotelm_por - minvtot));
			//wtd = 3.0 * (wsats_acrotelm - (wtot_acrotelm + itot_acrotelm)) / (2.0 * (wsats_acrotelm / ACROTELM_DEPTH - minvtot));
		}
		
		wtd = min(ACROTELM_DEPTH, wtd);

		wtp = -wtd;
	}

	// Water content in each acrotelm soil layer depending on wtd

	// Save acrotelm water mass before update
	wtot_acrotelm_init = wtot_acrotelm;

	// Depth of acrotelm sublayers - 10 mm
	double Dz_sub = 10.0;

	// Number of sublayers in acrotelm
	int nsublayers_acro = static_cast<int>(std::round(min(ACROTELM_DEPTH, depth) / Dz_sub));

	// Water content in each sublayer
	double surfw_sat = 0.0;
	double value[NSUBLAYERS_ACRO];
	double sub_water[NSUBLAYERS_ACRO];

	// Water table below surface
	if (wtd > 0.0) {

		// Determine the subsurface layer (0 to nsublayers_acro) within which wtd lies
		// Possible values:
		// = 0 if wtd < 10 mm 
		// = nsublayers_acro if wtd = ACROTELM_DEPTH
		int wtd_layer = (int)(wtd / Dz_sub);

		// Water content at the surface (Granberg et al. 1999 eqn 3)
		surfw_sat = max(minvtot, acrotelm_por - az * wtd);

		if (wtd_layer == 0) {
			// Total saturation throughout the acrotelm (from 10mm down)
			for (int lyr = 0; lyr < nsublayers_acro; lyr++) {
				value[lyr] = acrotelm_por;
			}
		}
		else {
			// Unsaturated zone in the acrotelm (Granberg et al. 1999 eqn 1)			
			for (int lyr = 0; lyr < wtd_layer; lyr++) {
				double z = (lyr + 1) * Dz_sub;
				value[lyr] = min(acrotelm_por, surfw_sat + (acrotelm_por - surfw_sat) * pow(z / wtd, 2));
			}

			// Saturated zone in the acrotelm
			for (int lyr = wtd_layer; lyr < nsublayers_acro; lyr++) {
				value[lyr] = acrotelm_por;
			}
		}

		// Now calculate the volumetric water content throughout the acrotelm
		sub_water[0] = 0.5 * (surfw_sat + value[0]);
		for (int lyr = 1; lyr < nsublayers_acro; lyr++) {
			sub_water[lyr] = 0.5 * (value[lyr] + value[lyr - 1]);
		}

		// Water content in each acrotelm layer
		wtot_acrotelm = 0.0;
		int sublyr = 0;
		for (int lyr = 0; lyr < nacrotelm; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			// Number of sublayers in acrotelm layer
			int sublayers_per_layer = static_cast<int>(std::round(soillayer.thickness / Dz_sub));

			double layerwater = 0.0;
			for (int sl = 0; sl < sublayers_per_layer; sl++) {
				layerwater += sub_water[sublyr] * Dz_sub;
				sublyr++;
			}

			// Water should not exceed saturation
			soillayer.wmass = min(soillayer.wsats - soillayer.imass, max(0.0, layerwater - soillayer.imass));
			wtot_acrotelm += soillayer.wmass;
		}
	}
	else {
		// WTD <= 0.0, at or above the surface - standing water
		wtot_acrotelm = 0.0;
		// Acrotelm saturated (subtract ice content)
		for (int lyr = 0; lyr < nacrotelm; lyr++) {
			Soillayer& soillayer = soillayers[lyr];
			soillayer.wmass = soillayer.wsats - soillayer.imass;
			wtot_acrotelm += soillayer.wmass;
		}
	}

	// Catotelm assumed to be always saturated
	for (int lyr = nacrotelm; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		double wadd = max(0.0, soillayer.wsats - (soillayer.wmass + soillayer.imass));
		waddon += wadd;
		soillayer.wmass += wadd;
	}

	// Check for error in acrotelm water balance
	double wtot_acrotelm_error = wtot_acrotelm - wtot_acrotelm_init;
	if (!negligible(wtot_acrotelm_error, -12)) {

		// Too much water
		if (wtot_acrotelm_error > 0.0) {

			// Water available for reduction
			double water_avail = 0.0;
			for (int lyr = 0; lyr < nacrotelm; lyr++) {
				Soillayer& soillayer = soillayers[lyr];
				if (lyr == 0) {
					water_avail += max(0.0, soillayer.wmass - minvtot * soillayer.thickness);
				}
				else {
					water_avail += max(0.0, soillayer.wmass - soillayer.wp);
				}
			}

			// Remaining water fraction per layer
			double frac = water_avail > 0.0 ? min(1.0, wtot_acrotelm_error / water_avail) : 0.0;

			// Subtract water
			for (int lyr = 0; lyr < nacrotelm; lyr++) {
				Soillayer& soillayer = soillayers[lyr];
				double wsub;
				if (lyr == 0) {
					wsub = frac * max(0.0, soillayer.wmass - minvtot * soillayer.thickness);
				}
				else {
					wsub = frac * max(0.0, soillayer.wmass - soillayer.wp);
				}
				soillayer.wmass -= wsub;
				wtot_acrotelm_error -= wsub;
			}

			// If still too much water, then account it as added water
			waddon += wtot_acrotelm_error;
		}
		// Water is missing. Fill from bottom
		else {
			for (int lyr = nacrotelm - 1; lyr >= 0; lyr--) {

				Soillayer& soillayer = soillayers[lyr];

				// Space for adding water
				double avail = max(0.0, soillayer.wsats - (soillayer.wmass + soillayer.imass));

				// Water to add to this layer
				double wmass_to_add = max(0.0, min(avail, -wtot_acrotelm_error));

				// Update fraction available water in layer and reduce error
				soillayer.wmass += wmass_to_add;
				wtot_acrotelm_error += wmass_to_add;
			}

			// If still too little water, then account it as leached water
			runoff_surf -= wtot_acrotelm_error;
		}
	}

	// Remove missing water from runoff
	double wmass_moved = min(runoff_surf, waddon);
	runoff_surf -= wmass_moved;
	waddon -= wmass_moved;
	wmass_moved = min(runoff_ssw, waddon);
	runoff_ssw -= wmass_moved;
	waddon -= wmass_moved;

	// Report water fluxes
	Fluxes& fluxes = patch.fluxes;

	// Evapotranspiration data
	if (empirical_evapotranspiration) {
		fluxes.report_flux(Fluxes::WEVAP, fevap * evapotranspiration);
		fluxes.report_flux(Fluxes::WTRAN, (1.0 - fevap) * evapotranspiration);
	}
	else {
		fluxes.report_flux(Fluxes::WEVAP, evap);
		fluxes.report_flux(Fluxes::WTRAN, aet_total);
	}

	// Surface runoff
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, runoff_ssw);

	// Runon
	fluxes.report_flux(Fluxes::WADDED, runon + waddon);

	// Update CO2 level in acrotelm
	update_acrotelm_co2(climate.co2);

	// Update limits on moss photosynthesis

	// The cap is kept at 100% while WTP is above wtp_moss_upper ((Wania et al 2009b) have -150 mm), decreases linearly
	// between wtp_moss_upper and -280mm  This has to do with the water
	// content of the mosses as they are less able to photosynthesize
	// as they dry out. NO MOSS DESSICATION

	// Mosses can be too competitive when wtp_moss_upper set to -150 mm, as in Wania et al. 2009b 
	const double wtp_moss_upper = 0.0;
	const double wtp_moss_lower = -280.0;
	const double moss_wtp_limit_lowerlimit = 0.3;

	// No limit
	if (wtp > wtp_moss_upper) {
		dmoss_wtp_limit = 1.0;
	}
	// Middle
	else if (wtp > wtp_moss_lower) {
		dmoss_wtp_limit = 1.0 - (-wtp + wtp_moss_upper) * (1.0 - moss_wtp_limit_lowerlimit) / (wtp_moss_upper - wtp_moss_lower);
	}
	// Lower limit for wtp between -300mm and -280mm
	else {
		dmoss_wtp_limit = moss_wtp_limit_lowerlimit; 
	}


	// Update limits on graminoid photosynthesis

	// LPJ-WHy code & comment:

	// The following lines will reduce the mgpp of
	// flood-tolerant grasses when the water table drops
	// below -10cm. This should reflect that these grasses
	// are usually found where the water table stays close
	// to the surface for most of the time (i.e. in fens).
	// Note that mwtp is negative!

	// mgpp(m,pft) = mgpp(m,pft) + mgpp(m,pft) *
	//               (mwtp(m) + 100) / 200

	// Upper limit (mm - Wania et al. 2009b have -100, could be -50.0)
	const double wtp_graminoid_upper = -10;
	// Lower limit (mm - Wania et al. 2009b have -300, could be -200.0)
	const double wtp_graminoid_lower = -100;
	// Lowest limit (mm)
	const double wtp_graminoid_lowestlimit = 0.0;

	// More restrictive than Wania et al., otherwise graminoids become too dominant as they never suffer from inundation stress
	// The cap drops from 1 as soon as wtp goes below wtp_graminoid_upper
	// No limit
	if (wtp >= wtp_graminoid_upper) {
		dgraminoid_wtp_limit = 1.0;
	}
	// Lower limit
	else if (wtp > wtp_graminoid_lower) {
		dgraminoid_wtp_limit = 1.0 - (-wtp + wtp_graminoid_upper) * (1.0 - wtp_graminoid_lowestlimit) / (wtp_graminoid_upper - wtp_graminoid_lower);
	}
	// Lowest limit
	else {
		dgraminoid_wtp_limit = wtp_graminoid_lowestlimit;
	}

	// Update output variables

	// Reset monthly wtp and mmoss_wtp_limit averages on Jan 1
	if (date.day == 0) {
		awtp = 0.0;
		for (int mth = 0; mth < 12; mth++) {
			mwtp[mth] = 0.0;
		}
	}

	mwtp[date.month] += wtp / (double)date.ndaymonth[date.month];
	awtp += wtp / (double)date.year_length();

	return; // ...if no problems
}

bool Soil::is_wet_season() {

	Landcover& lc = patch.stand.get_gridcell().landcover;

	// Are we in wet period
	bool is_saturated;
	if (lc.wet_startday + lc.wet_days > 365) {
		is_saturated = date.day >= lc.wet_startday || date.day <= lc.wet_startday + lc.wet_days - 365;
	}
	else {
		is_saturated = date.day >= lc.wet_startday && date.day <= lc.wet_startday + lc.wet_days;
	}

	// If no wet duration is read, then constant wet
	return !peatduration || is_saturated;
}

/// Updates the water content for each soil layer in mineral wetland soils
/**
* Call each simulation day (if is_true_wetland_stand true).
* Soil class method used to calculate water need to keep saturation or 
* surface runoff if water exceeds saturation depending on given
* snow melt, rainfall, evapotranspiration from vegetation (AET) 
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_wetland(const Climate& climate, double fevap) {

	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET), evaporation and runoff
	// Evaporation depth is set to the top soil layer

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)
	// evap       = evaporation (mm/day)
	// wmass_wetland_added = water added to soil to keep saturation (mm)

	// Minimum water considerred (equals to a 0.1 mm3)
	const double EPS = 1.0e-9;

	double wmass_begin = get_sum(&Soillayer::wmass) + rain_melt - wmass_wetland_added;

	double evap = 0.0;

	double runoff_surf = 0.0;

	// AET for a particular layer and individual (mm/day)
	double aet;
	// total AET (mm/day)
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals and extract it from soil lalyers
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			//double fwu = patch.pft[indiv.pft.id].fwuptake[lyr]; // DW_COMMENT not needed but crashed otherwise on Simba???????
			aet = indiv.get_fwuptake(lyr) * indiv.aet;
			aet_total += aet;
		}
		vegetation.nextobj();
	}

	// Calculate wcont for top layer (evaporation layer)
	wcont_evap = soillayers[0].wcont();

	// Calculate total evaporation only if snow depth < 10mm
	if (snow.depth() < 10.0 && !negligible(fevap)) {
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}
	else {
		evap = 0.0;
	}

	// Water budget
	double wbudget = rain_melt - evap - aet_total;
	rain_melt = 0.0;

	// Missing water
	if (wbudget <= 0.0) {
		wmass_wetland_added -= wbudget;
	}
	// Excess was is lost through surface runoff
	else {
		runoff_surf = wbudget;
	}

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, 0.0);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, 0.0);
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, 0.0);
	fluxes.report_flux(Fluxes::WEVAP, evap);
	fluxes.report_flux(Fluxes::WTRAN, aet_total);
	fluxes.report_flux(Fluxes::WADDED, wmass_wetland_added);

	// Calculate annual wcont for growing season (temperature > 5C)
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Reset the awcont array on the first day of every year
		if (date.day == 0) {
			soillayer.awcont = 0.0;
			if (lyr == 0) patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			soillayer.awcont += soillayer.wcont();
			if (lyr == 0) patch.growingseasondays++;
		}

		// Do the averaging on the last day of every year
		if (date.islastday && date.islastmonth) {

			// In case it's never warm enough:
			if (!patch.growingseasondays) {
				soillayer.awcont = 1.0;
			}
			else {
				soillayer.awcont /= (double)patch.growingseasondays;
			}
		}
	}

	// Water balance check
	double wmass_end = get_sum(&Soillayer::wmass) + runoff_surf + evap + aet_total - wmass_wetland_added;

	if (fabs(wmass_begin - wmass_end) > 1.0e-10)
		dprintf("Year %d Day %d WMASS ERROR: Diff %g Begin %g End %g\n", date.year, date.day, wmass_begin - wmass_end, wmass_begin, wmass_end);
}

/// Updates the water content for each soil layer
/**
* Call each simulation day (if hydrology set to GERTEN).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_gerten(const Climate& climate, double fevap) {

	// DESCRIPTION
	// This version replicates the hydrology_lpjf function in LPJ-GUESS v4.0, so 
	// it uses 2 soil layers, and does not consider impedance or ice formation. 

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wmass      = water content of soil layers [0=upper layer, nsoillayers-1=bottom layer] as
	//              water mass (kg m-2 or mm)
	// awcont_upper = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)
	// runoff     = total daily runoff from all soil layers (mm/day)

	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET) and percolation between layers;
	// calculation of runoff

	// INPUT PARAMETERS
	// rain_melt  = inward water flux to soil today (rain + snowmelt) (mm)
	// perc_base  = coefficient in percolation calculation (K in Eqn 31, Haxeltine
	//              & Prentice 1996)
	// perc_exp   = exponent in percolation calculation (=4 in Eqn 31, Haxeltine &
	//              Prentice 1996)
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface
	// snowpack   = depth of snow (mm)

	// INPUT AND OUTPUT PARAMETERS
	// wmass      = water content of soil layers [0=upper layer, nsoillayers-1=bottom layer] as
	//              water mass (kg m-2 or mm)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont_upper = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double BASEFLOW_FRAC = 0.5;
	// Fraction of total (vegetation) AET from upper soil layer that is derived
	// from the top K_DEPTH (fraction) of the upper soil layer
	// (parameters for calculating K_AET_DEPTH below)
	const double K_DEPTH = 0.4;
	const double K_AET = 0.52;
	// Weighting coefficient for AET flux from evaporation layer, assuming active
	//   root density decreases with soil depth
	// Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
	//   K_DEPTH=0.4, K_AET=0.52
	const double K_AET_DEPTH = (max(1.0, soillayers[0].thickness / SOILDEPTH_EVAP) - 1.0) *
		(K_AET / K_DEPTH - 1.0) / (1.0 / K_DEPTH - 1.0) + 1.0;

	// Water mass balance check
	double wbefore = get_sum(&Soillayer::wmass) + rain_melt;

	// Initial water masses of each layer
	double wmass_org[NSOILLAYERS_MAX] = { 0.0 };
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wmass_org[lyr] = soillayers[lyr].wmass;
	}

	// Actual Evapotranspiration from soil 

	// AET for a particular layer and individual (mm)
	double aet;
	// total AET for each Gerten-type soil layer (mm)
	double aet_layer[NSOILLAYERS_SIMPLE] = { 0.0 };
	// total AET (mm)
	double aet_total = 0.0;

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// Loop over all layers
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			aet = indiv.get_fwuptake(lyr) * indiv.aet;
			aet_total += aet;
			aet_layer[lyr] += aet;
		}
		vegetation.nextobj();
	}

	// Evaporation from soil surface (Bondeau et al. (2007),  Rost et al. (2008))
	
	Soillayer& soillayer_top = soillayers[0];

	double evap = 0.0;

	// Evaporation only if snow equivalent water < 10mm
	if (snow.water_and_ice_content() < 10.0) {	
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}

	// Implement in- and outgoing fluxes to upper soil layer
	soillayer_top.wmass += rain_melt - aet_layer[0] - evap;

	// Save water influx for water movement calculation
	double wmass_infilt_tot = rain_melt - evap;

	// Two-layer hydrology always have water up to wilting point
	if (soillayer_top.wmass <= soillayer_top.wp) { 
		evap += soillayer_top.wmass - soillayer_top.wp;
		soillayer_top.wmass = soillayer_top.wp;
	}

	// Surface runoff if surface soil layer has water over field capacity
	double runoff_surf = 0.0;
	if (soillayer_top.wmass > soillayer_top.fc()) {
		runoff_surf = soillayer_top.wmass - soillayer_top.fc();
		soillayer_top.wmass = soillayer_top.fc();
	}

	// Update water content in evaporation layer for tomorrow
	wcont_evap += (rain_melt - aet_layer[0] * K_AET_DEPTH * min(1.0, SOILDEPTH_EVAP / soillayer_top.thickness) - evap)
		/ soillayer_top.awc;

	if (wcont_evap > soillayer_top.wcont()) {
		wcont_evap = soillayer_top.wcont();
	}
	else if (wcont_evap < 0.0) {
		wcont_evap = 0.0;
	}

	// Percolation from evaporation layer
	double perc = 0.0;
	if (percolate) {
		perc = min(min(1.0, SOILDEPTH_EVAP / soillayer_top.thickness) * soillayer_top.perc_base * pow(wcont_evap, soillayer_top.perc_exp),
			max_rain_melt);
	}
	wcont_evap -= min(wcont_evap, perc / soillayer_top.awc);

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	for (int lyr = 1; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];
		Soillayer& soillayer_above = soillayers[lyr - 1];

		// Percolation

		// Allow only on days with rain or snowmelt
		if (percolate) {
			perc = min(soillayer_above.perc_base * pow(soillayer_above.wcont(), soillayer_above.perc_exp), max_rain_melt);
		}
		else {
			perc = 0.0;
		}

		// Percolation water from above layer (mm)
		soillayer_above.perc_mass = min(perc, soillayer_above.wmass - soillayer_above.wp);

		// Remove and transfer water (mm)
		soillayer_above.wmass -= soillayer_above.perc_mass;
		soillayer.wmass += soillayer_above.perc_mass;
		if (soillayer.wmass > soillayer.fc()) {
			soillayer.runoff_drain = soillayer.wmass - soillayer.fc();
			soillayer.wmass = soillayer.fc();
		}
		else {
			soillayer.runoff_drain = 0.0;
		}

		// Deduct AET from this soil layer
		soillayer.wmass -= min(soillayer.wmass - soillayer.wp, aet_layer[lyr]);
	}

	// Baseflow runoff (rain or snowmelt days only)
	Soillayer& soillayer_bottom = soillayers[nsoillayers - 1];
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC * soillayer_bottom.perc_base * 
			pow(soillayer_bottom.wcont(), soillayer_bottom.perc_exp);
		
		// Check so baseflow doesn't exceed inflow of water
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer
		soillayer_bottom.perc_mass = min(perc_baseflow, soillayer_bottom.wmass - soillayer_bottom.wp);
		soillayer_bottom.wmass -= soillayer_bottom.perc_mass;
		soillayer_bottom.runoff_baseflow = soillayer_bottom.perc_mass;
	}
	else {
		soillayer_bottom.runoff_baseflow = 0.0;
	}

	// Resolve todays water movements
	wfluxes[0] = wmass_infilt_tot;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wfluxes[lyr + 1] = wfluxes[lyr] + wmass_org[lyr] - (soillayer.wmass + soillayer.runoff_drain + soillayer.runoff_baseflow + aet_layer[lyr]);
	}

	// Water balance check
	double wafter = get_sum(&Soillayer::wmass) + runoff_surf + get_sum(&Soillayer::runoff_drain) + get_sum(&Soillayer::runoff_baseflow) + aet_total + evap;
	if (fabs(wbefore-wafter) > 1.0e-10)
		dprintf("Year %d Day %d BUCKETWATER Before %g After %g Diff %g Rain_melt %g\n",date.year,date.day, wbefore, wafter, wbefore - wafter, rain_melt);

	// Reset rain_melt
	rain_melt = 0.0;

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, 0.0);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(&Soillayer::runoff_drain));
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(&Soillayer::runoff_baseflow));
	fluxes.report_flux(Fluxes::WEVAP, evap);
	fluxes.report_flux(Fluxes::WTRAN, aet_total);

	// Update awcont_upper

	// Reset awcont_upper on the first day of every year
	if (date.day == 0) {
		awcont_upper = 0.0;
		patch.growingseasondays = 0;
	}

	// If it's warm enough for growth, update awcont_upper with this day's wcont
	if (climate.temp > 5.0) {
		awcont_upper += soillayer_top.wcont();
		patch.growingseasondays++;
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}
}

/// Updates the water content for each soil layer
/**
* Call each simulation day (if hydrology set to MILLER).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_miller(const Climate& climate, double fevap) {

	// DESCRIPTION
	// Soil class method used to update the water content for each soil layer 
	// given snow melt, rainfall, evapotranspiration from vegetation (AET) and 
	// percolation between layers. Runoff is also calculated.
	// Updated to deal with 15 soil layers and ice formation.  

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wmass      = water mass of soil layers [mm]
	// wcont_evap = water mass of evaporation sublayers (20 cm) at top of soil 
	// awcont     = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)
	// runoff     = total daily runoff from all soil layers (mm/day)

	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double SOILDEPTH_EVAP = 200.0;
	// Number of layers in the evaporation depth
	const int NEVAPLAYERS = 2;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double BASEFLOW_FRAC = 0.5;
	// Minimum water amount for percolation to occur (mm)
	const double MIN_WATER_PERC = 0.0000001;
	// max error allowed (mm)
	const double maxerr = 0.0001;

	// Initial water masses of each layer
	double wmass_org[NSOILLAYERS_MAX] = { 0.0 };
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wmass_org[lyr] = soillayers[lyr].wmass;
	}

	// Initial estimate of evaporation

	double evap_init = 0;
	double wmass_evap_init = 0.0;
	double awc_evap = 0.0;
	double wcont_evap_init = 0.0;

	for (int lyr = 0; lyr < NEVAPLAYERS; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wmass_evap_init += soillayer.wcont() * soillayer.awc;
		awc_evap += soillayer.awc;
	}

	wcont_evap_init = wmass_evap_init / awc_evap;

	// Only potential evaporation if snow depth < 10mm
	if (snow.depth() < 10.0) {
		evap_init = climate.eet * PRIESTLEY_TAYLOR * wcont_evap_init * wcont_evap_init * fevap;
	}

	// Transpiration

	// AET for a particular layer and individual (mm)
	double aet = 0.0;
	// Total AET (mm)
	double aet_total = 0.0;
	// Total AET for each soil layer (mm)
	double aet_layer[NSOILLAYERS_MAX] = { 0.0 };

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			aet = indiv.get_fwuptake(lyr) * indiv.aet;
			aet_layer[lyr] += aet;
		}
		vegetation.nextobj();
	}

	// Water that can still be added to each soil layer (mm)
	double wmass_pot_layer[NSOILLAYERS_MAX] = { 0.0 };

	// Top and bottom GERTEN layers, initial water and potential water that can be added
	double wmass_layer1 = 0.0;
	double wmass_layer2 = 0.0;
	double wmass_pot_layer1 = 0.0;
	double wmass_pot_layer2 = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Remove AET from each soil layer
		double aet_lyr = min(aet_layer[lyr], soillayer.wmass);
		soillayer.wmass -= aet_lyr;
		aet_total += aet_lyr;

		// Determine potential for each soil layer
		wmass_pot_layer[lyr] = max(0.0, soillayer.fc() - (soillayer.wmass + soillayer.imass));

		// Determine potential for GERTEN layers
		if (lyr < NSOILLAYERS_UPPER) {
			wmass_pot_layer1 += wmass_pot_layer[lyr];
		}
		else {
			wmass_pot_layer2 += wmass_pot_layer[lyr];
		}
	}

	// Evaporation

	// Update water mass available for evaporation
	double wmass_evap = 0.0;
	for (int lyr = 0; lyr < NEVAPLAYERS; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wmass_evap += soillayer.wcont() * soillayer.awc;
	}

	// Evaporation from bare surface (after AET removed)
	double evap = min(evap_init, max(0.0, wmass_evap));

	// Removal of water from the evaporation layers of the soil
	if (evap > 0.0) {

		for (int lyr = 0; lyr < NEVAPLAYERS; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			// Remove the water in proportional to the available water remaining in each layer
			double proportion = evap * (soillayer.wcont() * soillayer.awc / wmass_evap);

			// Remove water
			soillayer.wmass -= proportion;

			// And increae potential to add water
			wmass_pot_layer[lyr] += proportion;
		}
	}
	else {
		// No evaporation in this case
		evap = 0.0;
	}

	// Implement in- and outgoing fluxes to the top 50cm soil layers, taking into account AET in those layers and evaporation
	double wmass_flux_in = rain_melt;
	rain_melt = 0.0;

	// Initialise runoff variables

	// Surface runoff (mm) 
	double runoff_surf = 0.0;
	// Drainage runoff (mm)
	double runoff_drain = 0.0;
	// Baseflow runoff (mm) (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_baseflow = 0.0;

	// Input to top GERTEN layer
	if (wmass_flux_in > 0) {

		// Potential in the top layer
		double wmass_pot_top_layer = 0.0;
		for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {
			wmass_pot_top_layer += wmass_pot_layer[lyr];
		}

		// Addition of water to the top layers (rain_melt > 0)
		if (wmass_flux_in < wmass_pot_top_layer) {

			// The upper soil layers can absorb all of this water today
			// so add the water to layers in proportion to their capacity
			for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				double water_input = wmass_flux_in * (wmass_pot_layer[lyr] / wmass_pot_top_layer);

				soillayer.wmass += water_input;
				wmass_pot_layer[lyr] -= water_input;
			}

			// No surface runoff
			runoff_surf = 0.0;
		}
		else {

			// The upper soil layers cannot absorb all of this water today, so take up what can be absorbed 
			// (i.e. wmass_pot_top_layer) and add the rest to surface runoff
			for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {
				Soillayer& soillayer = soillayers[lyr];

				soillayer.wmass = soillayer.fc() - soillayer.imass;
				wmass_pot_layer[lyr] = 0.0;
			}

			// Update surface runoff
			runoff_surf = wmass_flux_in - wmass_pot_top_layer;
		}
	}

	// Total water infiltration for water flux calculations
	double wmass_infilt_tot = wmass_flux_in - (runoff_surf + evap);

	// Update wmass_layer1

	double awc1 = 0.0;
	wmass_layer1 = 0.0;
	for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wmass_layer1 += soillayer.wcont() * soillayer.awc;
		awc1 += soillayer.awc;
	}

	double wcont_layer1 = wmass_layer1 / awc1;

	// Percolation
	// We assume there are two percolation layers. 
	// 1: the top 50cm soil layers considered in the Gerten et al. scheme
	// 2: Lower 1.0 m layers considered in the Gerten et al. scheme

	if (percolate) {

		// Percolation from top layer if there is enough water (at least MIN_WATER_PERC)
		double perc = 0.0;

		if (wmass_layer1 > MIN_WATER_PERC) {

			// Percolation from the top layer (mm; limited to available liquid water)
			perc = min(soillayers[0].perc_base * pow(wcont_layer1, soillayers[0].perc_exp), wmass_layer1);

			for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {
				Soillayer& soillayer = soillayers[lyr];

				// Take water proportionally from the layers.
				double perc_layer = perc * (soillayer.wcont() * soillayer.awc / wmass_layer1);
				soillayer.wmass -= perc_layer;
				wmass_pot_layer[lyr] += perc_layer;
			}
		}

		// Percolation and fluxes to and from the soil layers below the top 50cm

		// Excess water transferred to runoff
		// Eqns 26, 27, 31, Haxeltine & Prentice 1996

		double water_input = 0.0;

		// First we need to add the water percolated from above (perc) to this layer
		if (perc > wmass_pot_layer2) {
			// Fill the layer and reduce water_input
			water_input = wmass_pot_layer2;
			// The excess is run off and not percolated further (added to runoff_drain)
			runoff_drain = perc - wmass_pot_layer2;
		}
		else {
			water_input = perc;
		}

		double wmass_perc_layer_2 = 0.0;
		double awc_perc_layer_2 = 0.0;
		double wcont_perc_layer_2 = 0.0;

		if (water_input > 0.0) {

			// Add any water (water_input) to layers 5-14 in proportion to the capacity
			for (int lyr = NSOILLAYERS_UPPER; lyr < nsoillayers; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				wmass_layer2 += soillayer.wcont() * soillayer.awc;

				double water_input_lyr = water_input * (wmass_pot_layer[lyr] / wmass_pot_layer2);

				soillayer.wmass += water_input_lyr;
				wmass_pot_layer[lyr] -= water_input_lyr;

				wmass_perc_layer_2 += soillayer.wcont() * soillayer.awc;
				awc_perc_layer_2 += soillayer.awc;
			}

			// Update wcont for these layers before percolation
			wcont_perc_layer_2 = wmass_perc_layer_2 / awc_perc_layer_2;
		}

		// Percolation from base layer
		double perc_from_base = 0.0;
		if (percolate && wmass_perc_layer_2 > MIN_WATER_PERC) {

			perc_from_base = min(wmass_perc_layer_2, min(BASEFLOW_FRAC * soillayers[NSOILLAYERS_UPPER].perc_base * 
				pow(wcont_perc_layer_2, soillayers[NSOILLAYERS_UPPER].perc_exp), max_rain_melt));

			// As in LPJ-GUESS v4.0
			if (perc_from_base > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
				perc_from_base = rain_melt - runoff_surf;
			}
		}
		else {
			// No percolation
			perc_from_base = 0.0;
		}

		// Now remove the water (i.e. perc_from_base) from layers 5-14 in proportion to the available water
		if (perc_from_base > 0) {

			for (int lyr = NSOILLAYERS_UPPER; lyr < nsoillayers; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				// Remove the water (perc_from_base) in proportional to the available water remaining in each layer
				double proportion = perc_from_base * (soillayer.wcont() * soillayer.awc / wmass_perc_layer_2); // mm
				// Subtract water
				soillayer.wmass -= proportion;
				// Increase potential to add water
				wmass_pot_layer[lyr] += proportion;
			}
		}

		// Baseflow
		runoff_baseflow = perc_from_base;
	}
	else {
		// No percolation at all, so reset runoff to 0 
		runoff_baseflow = 0.0;
		runoff_drain = 0.0;
	}

	double runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// water added when patch.stand.is_true_wetland_stand() should be be subtracted from runoff
	// in proportion to its components
	if (patch.stand.is_true_wetland_stand() && ifsaturatewetlands) {

		if (runoff <= wmass_wetland_added && runoff > 0.0) {

			// Not enough runoff to balance the water added to the wetland, i.e. we added more than was run off.

			// Update wetland_water_added by taking back the amount used for runoff.
			wmass_wetland_added -= runoff;

			// All runoff removed 
			runoff_surf = 0.0;
			runoff_baseflow = 0.0;
			runoff_drain = 0.0;
			runoff = 0.0;
		}
		else if (runoff > wmass_wetland_added && runoff > 0.0) {

			// Enough runoff to balance the water added to the wetland so we take it back.
			runoff_surf -= wmass_wetland_added * runoff_surf / runoff;
			runoff_drain -= wmass_wetland_added * runoff_drain / runoff;
			runoff_baseflow -= wmass_wetland_added * runoff_baseflow / runoff;
			// Recalculate total runoff 
			runoff = runoff_surf + runoff_drain + runoff_baseflow;

			// Reset wetland_water_added
			wmass_wetland_added = 0.0;
		}
	}

	// Assign runoff to layers
	if (wmass_layer2 > 0.0) {
		for (int lyr = NSOILLAYERS_UPPER; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];
			soillayer.runoff_drain = runoff_drain * soillayer.wcont() * soillayer.awc / wmass_layer2;
		}
	}
	soillayers[nsoillayers - 1].runoff_baseflow = runoff_baseflow;

	// Resolve todays water movements
	wfluxes[0] = wmass_infilt_tot;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wfluxes[lyr + 1] = wfluxes[lyr] + wmass_org[lyr] - (soillayer.wmass + soillayer.runoff_drain + soillayer.runoff_baseflow + aet_layer[lyr]);
	}

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, 0.0);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, runoff_drain);
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, runoff_baseflow);
	fluxes.report_flux(Fluxes::WEVAP, evap);
	fluxes.report_flux(Fluxes::WTRAN, aet_total);
	fluxes.report_flux(Fluxes::WADDED, wmass_wetland_added);

	// Drought limited establishment - update awcont_upper
	// Original algorithm by Thomas Hickler
	for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {

		// Reset on the first day of every year
		if (date.day == 0 && lyr == 0) {
			awcont_upper = 0.0;
			patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			awcont_upper += soillayers[lyr].wcont() / (double)NSOILLAYERS_UPPER;
			if (lyr == 0) {
				patch.growingseasondays++;
			}
		}
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}
}

/// Updates the water content for each soil layer in non-peatland soils
/**
* Call each simulation day (if hydrology set to RICHARDS).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soil::hydrology_richards(const Climate& climate, double fevap) {

	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET), evaporation and percolation between
	// layers and runoff
	// Infiltration and evaporation depth is set to the top soil layer

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wmass      = water mass content of soil layers

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)
	// evap       = evaporation (mm/day)

	// Minimum water considerred (equals to a 0.1 mm3)
	const double EPS = 1.0e-9;

	// Maximum effective depth of surface runoff - 0.1m
	const double DEPTH_INFILT = 100.0;

	// Evaposation (mm)
	double evap = 0.0;
	// Surface runoff (mm)
	double runoff_surf = 0.0;
	// Standing water runoff (mm)
	double runoff_ssw = 0.0;
	// Infiltration (mm)
	double wmass_infilt_tot = 0.0;
	// Infiltration water temperature
	double temp_infilt_tot = 0.0;
	// AET for a particular layer and individual (mm/day)
	double aet;
	// total AET for each soil layer (mm/day)
	double aet_layer[NSOILLAYERS_MAX] = { 0.0 };
	// total AET (mm/day)
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals and extract it from soil lalyers
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			//double fwu = patch.pft[indiv.pft.id].fwuptake[lyr]; // DW_COMMENT not needed but crashed otherwise on Simba???????
			aet = indiv.get_fwuptake(lyr) * indiv.aet;
			aet_layer[lyr] += aet;
			aet_total += aet;
		}
		vegetation.nextobj();
	}

	// Mixed-form of the Richards equation hydrology

	double dz[NSOILLAYERS_MAX], vws_liq[NSOILLAYERS_MAX], wflux_root[NSOILLAYERS_MAX], dt_dz[NSOILLAYERS_MAX],
		Kw[NSOILLAYERS_MAX], dKwdw[NSOILLAYERS_MAX], Psi[NSOILLAYERS_MAX], dPsidw[NSOILLAYERS_MAX],
		qin[NSOILLAYERS_MAX], qout[NSOILLAYERS_MAX], dqidw0[NSOILLAYERS_MAX], dqidw1[NSOILLAYERS_MAX],
		dqodw1[NSOILLAYERS_MAX], dqodw2[NSOILLAYERS_MAX], rmx[NSOILLAYERS_MAX], amx[NSOILLAYERS_MAX],
		bmx[NSOILLAYERS_MAX], cmx[NSOILLAYERS_MAX], dwat[NSOILLAYERS_MAX], wmass_org[NSOILLAYERS_MAX];
	
	int icy_lyr[NSOILLAYERS_MAX], icalc_lyr[NSOILLAYERS_MAX], calc_lyr[NSOILLAYERS_MAX];
	const double WFPS_SAT = 0.99;
	const int ICY_TOP_LYRS_MAX = 2;

	// Initial sub daily timesteps on days with precipitation
	double nsubdays_prec = 12.0;
	// Initial sub daily timesteps on days with no precipitation
	double nsubdays_noprec = 6.0;
	// Minimum timestep (15 min)
	double dtmin = 900.0;
	// Model timestep length in seconds (daily timestep)
	double dtime = SECS_PER_DAY;
	// Tolerance (mm) to halve length of substep
	double xTolerUpper = 0.1;
	// Tolerance (mm) to double length of substep
	double xTolerLower = 0.01;
	// A very small number used to check for sub step completion
	double verySmall = 1.0e-8;
	// Maximum error in timestep calculation
	double errorMax;
	// Precipitation and snow melt (mm s-1)
	double wflux_rain_melt = rain_melt / SECS_PER_DAY;
	// Total evaporation (mm s-1)
	double wflux_evap;
	// Total surface runoff (mm s-1)
	double wflux_runoff_surf;
	// Total surface standing water runoff (mm s-1)
	double wflux_runoff_ssw;
	// Total standing surface water flux (mm s-1)
	double wflux_ssw;
	// Total soil infiltration (mm s-1)
	double wflux_infilt;
	// Water drained from the standing water into the soil (mm s-1) 
	double wflux_ssw_drain;

	// Initialize the number of substeps
	int nsubstep = 0;

	// Fixed duration of rain event (6 hours)
	const double dtrain = 6.0 * SECS_PER_HOUR;

	// Substep completed
	double dtdone = 0.0;

	// Estimate length of the substep (seconds)
	// If day with water influx then shorten initial length
	double dtsub = SECS_PER_DAY / nsubdays_noprec;
	if (wflux_rain_melt > 0.0 || ssw.wmass > 0.0) {
		dtsub = SECS_PER_DAY / nsubdays_prec;
	}

	// Finished day
	bool finished = false;

	// Calculate static ice mass and saturation level for infiltation depth (top layer)
	double imass_infilt = var_depth(0.0, DEPTH_INFILT, &Soillayer::imass);
	double wsats_infilt = var_depth(0.0, DEPTH_INFILT, &Soillayer::wsats);

	// Calculate ice impedance factor for infiltation depth (top layer)
	double imp_ice = ice_imp(imass_infilt / wsats_infilt);

	// Hydraulic conductivity at saturation for infiltation depth (top layer)
	double Kw_max = var_depth(0.0, DEPTH_INFILT, &Soillayer::Kw_max, &Soillayer::thickness);

	// Reset daily infiltation flux, evaporation and surface and drainage runoff 
	// and set initial water fraction to be able to calculate daily water fluxes
	// Also determine if a layer is icy and the amounts icy layers
	int nr_icy_lyrs = 0;
	int nr_icy_top_lyrs = 0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		soillayer.runoff_drain = 0.0;
		soillayer.runoff_baseflow = 0.0;

		// Save initial wmass
		wmass_org[lyr] = soillayer.wmass;

		// Determine layers without ice starting from the top. 
		// These are the layers to be included in the 
		// hydrological calculations
		if (ice_imp(soillayer.imass / soillayer.wsats) < 0.01) {

			// Frozen layer
			icy_lyr[lyr] = 1;

			// Allow top layers to be frozen and still include 
			// them in calculations
			if (nr_icy_lyrs == lyr) {
				icalc_lyr[lyr] = 1;
			}
			else {
				// Not consisten frozen top layer
				icalc_lyr[lyr] = 0;
			}

			// Total number of icy layers
			nr_icy_lyrs++;

			// Number of icy top layers
			if (lyr < ICY_TOP_LYRS_MAX) {
				nr_icy_top_lyrs++;
			}
		}
		else {
			// Not frozen
			icy_lyr[lyr] = 0;
			// Top layer
			if (lyr == 0) {
				icalc_lyr[lyr] = 1;
			}
			else {
				// Only include if all above layers have been included
				icalc_lyr[lyr] = icalc_lyr[lyr-1];
			}
		}
	}

	// Check if there is too much ice to do hydrology calculations
	bool calc_hyd = true;

	// Too many frozen layers to do any hydrology calculations
	if (nr_icy_top_lyrs == ICY_TOP_LYRS_MAX || nr_icy_top_lyrs == nsoillayers || nr_icy_lyrs >= nsoillayers - 1) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			icalc_lyr[lyr] = 0;
		}
		calc_hyd = false;
	}

	// Adaptive substep loop
	while (!finished) {

		// Calculate substep
		bool calcsub = true;

		// Number of soil layers to simulate
		int nsoillayers_sim = 0;

		// Calculate number of soil layers that should be included
		// in hydrology calculation. Depends on frozen and saturated 
		// layers. If saturated and/or frozen layers are included, then
		// the calculations become very slow. Saturated layers are handled
		// after lateral drainage is calculated
		if (calc_hyd) {

			// Start from bottom to find where standing liquid water is
			for (int lyr = nsoillayers - 1; lyr >= 0; lyr--) {
				bool sat_lyr = soillayers[lyr].wfps() >= WFPS_SAT;

				// Bottom layer
				if (lyr == nsoillayers - 1) {
					bool sat_lyr_above = soillayers[lyr - 1].wfps() >= WFPS_SAT;

					calc_lyr[lyr] = (!sat_lyr || !sat_lyr_above) * icalc_lyr[lyr];
				}
				// Top layer
				else if (lyr == 0) {
					calc_lyr[lyr] = (!sat_lyr || calc_lyr[lyr + 1]) * icalc_lyr[lyr];
				}
				// Middle layers
				else {
					bool sat_lyr_above = soillayers[lyr - 1].wfps() >= WFPS_SAT;
					calc_lyr[lyr] = (!sat_lyr || !sat_lyr_above || calc_lyr[lyr + 1]) * icalc_lyr[lyr];
				}

				// Add up how many layers that will be included in 
				// hydrology calculations
				nsoillayers_sim += calc_lyr[lyr];

				// Reset change in water
				dwat[lyr] = 0.0;
			}
		}

		// Minimum one layer needs to be simulated
		nsoillayers_sim = max(1, nsoillayers_sim);

		while (calcsub) {

			// Keep track of the number of trials
			nsubstep++;

			// If only one layer to simulate, do the rest of this timestep in one go
			if (nsoillayers_sim == 1) {
				dtsub = min(SECS_PER_DAY / nsubdays_noprec, dtime - dtdone);
			}

			// Calculate infiltration and surface runoff
			compute_evap_infilt_and_runoff(dtsub, wflux_rain_melt, DEPTH_INFILT, imass_infilt,
				wsats_infilt, imp_ice, Kw_max, fevap, wflux_evap, wflux_runoff_surf, wflux_runoff_ssw, 
				wflux_ssw, wflux_infilt, wflux_ssw_drain);

			// If only one layer to simulate, then solve directly
			if (nsoillayers_sim == 1) {

				// Extract transpiration 
				for (int lyr = 0; lyr < nsoillayers; lyr++) {
					soillayers[lyr].wmass -= aet_layer[lyr] * dtsub / SECS_PER_DAY;
				}

				// Influx of water
				double wmass_influx = (wflux_infilt + wflux_ssw_drain) * dtsub;

				// Add water to top two layers
				for (int lyr = 0; lyr < 2; lyr++) {
					double wmass_add = min(wmass_influx, soillayers[lyr].wsats - (soillayers[lyr].wmass + soillayers[lyr].imass));
					dwat[lyr] = wmass_add / soillayers[lyr].thickness;
					wmass_influx -= wmass_add;
				}

				// Reset lower layers influx
				for (int lyr = 2; lyr < nsoillayers; lyr++) {
					dwat[lyr] = 0.0;
				}

				// Excees influx water added to surface standing water
				wflux_ssw += wmass_influx / dtsub;

				// No error in mass balance
				errorMax = 0.0;

				// No need to recalculate sub timestep
				calcsub = false;
			}
			else {
				// Calculate matrix potential
				update_matric_potential();

				// Calculate hydraulic conductivity and diffusivity
				update_hydraulic_properties(nsoillayers_sim);

				// calculate commonly used variables
				for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {

					Soillayer& soillayer = soillayers[lyr];

					// saturated volumetric water content
					vws_liq[lyr] = soillayer.wsats / soillayer.thickness;
					// time per mm
					dt_dz[lyr] = dtsub / soillayer.thickness;
					Kw[lyr] = soillayer.Kw;
					Psi[lyr] = soillayer.Psi;
					dKwdw[lyr] = soillayer.dKwdw;
					dPsidw[lyr] = soillayer.dPsidw;
					dz[lyr] = soillayer.thickness;
					wflux_root[lyr] = aet_layer[lyr] / SECS_PER_DAY;
				}

				// Soil moisture fluxes and their derivatives
				compute_moisture_fluxes_and_derivs(wflux_infilt + wflux_ssw_drain, Kw, Psi, dKwdw, dPsidw, qin, qout,
					dqidw0, dqidw1, dqodw1, dqodw2, vws_liq, dz, nsoillayers_sim);

				// RHS of system of equations
				compute_RHS_moisture_form(wflux_root, qin, qout, dt_dz, rmx, nsoillayers_sim);

				// LHS of system of equations
				compute_LHS_moisture_form(dt_dz, dqidw0, dqidw1, dqodw1, dqodw2,
					amx, bmx, cmx, nsoillayers_sim);

				// solve tridiagonal matrix A.X = B
				TridiagonalCol(amx, bmx, cmx, rmx, dwat, nsoillayers_sim);

				// check if substep can be accepted
				check_error(dwat, dz, dt_dz, qin, qout, wflux_root, nsoillayers_sim, dtmin, xTolerUpper, dtsub, errorMax, calcsub);
			}
		}

		// Water added to standing water
		double wmass_ssw_added = wflux_ssw * dtsub;

		// Update standing water temperature // Not totally correct
		ssw.update_temp(rain_melt_temp, wmass_ssw_added);

		// Update standing water content (mm)
		ssw.wmass += wmass_ssw_added;

		// Very low surface standing water -> send to runoff
		if (negligible(ssw.water_and_ice_content(), -10)) {
			runoff_ssw += ssw.water_and_ice_content();
			ssw.wmass = 0.0;
			ssw.imass = 0.0;
		}

		// Accumulate evaporation
		evap += wflux_evap * dtsub;

		// Accumulate surface runoff
		runoff_surf += wflux_runoff_surf * dtsub;
		runoff_ssw += wflux_runoff_ssw * dtsub;

		// Temperature infiltration
		if (!negligible(wmass_infilt_tot + (wflux_infilt + wflux_ssw_drain) * dtsub)) {
			temp_infilt_tot = (temp_infilt_tot * wmass_infilt_tot + rain_melt_temp * wflux_infilt * dtsub +
				ssw.temp * wflux_ssw_drain * dtsub) / (wmass_infilt_tot + (wflux_infilt + wflux_ssw_drain) * dtsub);
		}

		// Accumulater infiltration
		wmass_infilt_tot += (wflux_infilt + wflux_ssw_drain) * dtsub;

		// Update the mass of liquid water
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].wmass += dwat[lyr] * soillayers[lyr].thickness;
		}

		// Extract root uptake from soil layer not accounted for
		if (nsoillayers_sim != 1 && nsoillayers_sim != nsoillayers) {
			for (int lyr = nsoillayers_sim; lyr < nsoillayers; lyr++) {
				soillayers[lyr].wmass -= aet_layer[lyr] * dtsub / SECS_PER_DAY;
			}
		}

		// drainage from the bottom of the soil 
		double baserunoff = 0.0;

		// compute drainage if not bedrock (mm s-1)
		if (!ifbedrock) {
			baserunoff = Kw[nsoillayers - 1] + dKwdw[nsoillayers - 1] * dwat[nsoillayers - 1];
		}

		// increment todays base runoff flux
		soillayers[nsoillayers - 1].runoff_baseflow += baserunoff * dtsub;

		// increment substep and check for completion of the substep
		dtdone += dtsub;

		// The excess water greater than soil moisture capacity (wsats) 
		// is recharged to the unsaturated layers above the current layer
		// (Seok et al. 2011)
		for (int lyr = nsoillayers - 1; lyr >= 0; lyr--) {

			Soillayer& soillayer = soillayers[lyr];

			if (soillayer.wmass + soillayer.imass > soillayer.wsats) {

				double wmass_excess = max(0.0, min(soillayer.wmass, soillayer.wmass - (soillayer.wsats - soillayer.imass)));

				if (lyr > 0) {
					soillayers[lyr - 1].wmass += wmass_excess;
				}
				else {
					// Update standing water content (mm)
					ssw.wmass += wmass_excess;
				}
				soillayer.wmass -= wmass_excess;
			}
		}

		// Lateral drainage
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			double wmass_avail = soillayer.wmass - max(0.0, soillayer.fc() - soillayer.imass);

			if (wmass_avail > 0.0) {

				// Slope dependant drainage of liquid water above fixed fraction of wsats (CLM5.0 eqn 7.106 and 107; mm s-1)
				// Depends on the side area of the soil layer in m -> soillayer.depth * M_PER_MM
				double drained_wmass = min(wmass_avail, C_dw * soillayer.Kw * ice_imp(soillayer.imass / soillayer.wsats) * slope_grade_tan * soillayer.thickness * M_PER_MM * dtsub);
				soillayer.runoff_drain += drained_wmass;
				soillayer.wmass -= drained_wmass;
			}
		}

		// For saturated zone move water downwards to always have saturated layer from the bottom
		if (nsoillayers_sim < nsoillayers) {
			for (int lyr = nsoillayers - 1; lyr >= nsoillayers_sim; lyr--) {
				double wdemand = min(soillayers[lyr - 1].wmass, soillayers[lyr].wsats - (soillayers[lyr].wmass + soillayers[lyr].imass));
				soillayers[lyr].wmass += wdemand;
				soillayers[lyr - 1].wmass -= wdemand;
			}
		}

		if (abs(dtime - dtdone) < verySmall) {
			// time step completed so exit the loop
			finished = true;
		}
		else {

			// check if the error is below the lower tolerance
			if (errorMax < xTolerLower * 2.0) {
				dtsub *= 2.0;
			}

			// ensure that the substep does not exceed the time remaining
			if (dtsub > dtime - dtdone) {
				dtsub = dtime - dtdone;
			}
		}
	}  // substep loop

	// Resolve todays water movements
	wfluxes[0] = wmass_infilt_tot;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		wfluxes[lyr + 1] = wfluxes[lyr] + wmass_org[lyr] - (soillayer.wmass + soillayer.runoff_drain + soillayer.runoff_baseflow + aet_layer[lyr]);
	}

	// Determine water table position
	update_wtp();

	// Update soil layer temperatures due to soil water fluxes
	wflow_temp(temp_infilt_tot);

	double drain = get_sum(&Soillayer::runoff_drain);
	double base = get_sum(&Soillayer::runoff_baseflow);

	double wmass_end = get_sum(&Soillayer::wmass) - rain_melt + (runoff_surf + 
		runoff_ssw + drain + base + evap + aet_total);

	// Reset ran_melt as water has entered the soil
	rain_melt = 0.0;

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, runoff_ssw);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(&Soillayer::runoff_drain));
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(&Soillayer::runoff_baseflow));
	fluxes.report_flux(Fluxes::WEVAP, evap);
	fluxes.report_flux(Fluxes::WTRAN, aet_total);

	// Calculate annual wcont for growing season (temperature > 5C)
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Reset the awcont array on the first day of every year
		if (date.day == 0) {
			soillayer.awcont = 0.0;
			if (lyr == 0) patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			soillayer.awcont += soillayer.wcont();
			if (lyr == 0) patch.growingseasondays++;
		}

		// Do the averaging on the last day of every year
		if (date.islastday && date.islastmonth) {

			// In case it's never warm enough:
			if (!patch.growingseasondays) {
				soillayer.awcont = 1.0;
			}
			else {
				soillayer.awcont /= (double)patch.growingseasondays;
			}
		}
	}
}

void Soil::compute_moisture_fluxes_and_derivs(double wflux_infl, double Kw[NSOILLAYERS_MAX], 
	double Psi[NSOILLAYERS_MAX], double dKwdw[NSOILLAYERS_MAX], double dPsidw[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX], 
	double qout[NSOILLAYERS_MAX], double dqidw0[NSOILLAYERS_MAX], double dqidw1[NSOILLAYERS_MAX], double dqodw1[NSOILLAYERS_MAX], 
	double dqodw2[NSOILLAYERS_MAX], double vws_liq[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], 
	int nsoillayers_sim) {

	// Calculate fluxes at the boundary of each layerand derivatives w.r.t.
	// relevant state variables

	/*
		dz !Input : layer thickness(m)
		wflux_infl !Input : infiltration(mm H2O / s)
		vws_liq !Input : volumetric soil water at saturation(porosity)
	*/

	double dKwdw1, dKwdw2, num, den;

	// Soil moisture fluxesand their derivatives

	// Top layer
	int tlyr = 0;

	// specify infiltration (mm H2O / s)
	qin[tlyr] = wflux_infl;

	// compute the flux derivative w.r.t influx
	dqidw1[tlyr] = 0.0;

	// compute derivatives in hydraulic conductivity at the interface w.r.t.
	// volumetric liquid water content in the layer aboveand the layer  below
	// NOTE : dKwds is the derivative in hydrualic conductivity at the
	// layer interface w.r.t relative saturation at the interface
	dKwdw1 = 0.5 * dKwdw[tlyr] / vws_liq[tlyr];   // derivative w.r.t.volumetric liquid water in the upper layer
	dKwdw2 = 0.5 * dKwdw[tlyr] / vws_liq[tlyr + 1]; // derivative w.r.t.volumetric liquid water in the lower layer

	// compute flux at the bottom of the lyr - th layer
	// NOTE: Kw[lyr] is hydraulic conductivity at the bottom of the lyr - th layer
	num = Psi[tlyr + 1] - Psi[tlyr];
	den = dz[tlyr];
	qout[tlyr] = -Kw[tlyr] * num / den + Kw[tlyr];

	// compute flux derivatives
	dqodw1[tlyr] = (Kw[tlyr] * dPsidw[tlyr] - dKwdw1 * num) / den + dKwdw1;
	dqodw2[tlyr] = (-Kw[tlyr] * dPsidw[tlyr + 1] - dKwdw2 * num) / den + dKwdw2;

	// Interior nodes

	for (int lyr = 1; lyr < nsoillayers_sim - 1; lyr++) {

		// get the flux from above
		qin[lyr] = qout[lyr - 1];

		// get the derivatives from above
		dqidw0[lyr] = dqodw1[lyr - 1];
		dqidw1[lyr] = dqodw2[lyr - 1];

		// compute derivatives in hydraulic conductivity at the interface w.r.t.
		// volumetric liquid water content in the layer above and the layer below
		// NOTE : dKwdw[lyr] is the derivative in hydrualic conductivity at the
		// layer interface w.r.t relative saturation at the interface
		dKwdw1 = 0.5 * dKwdw[lyr] / vws_liq[lyr];   // derivative w.r.t.volumetric liquid water in the upper layer
		dKwdw2 = 0.5 * dKwdw[lyr] / vws_liq[lyr + 1]; // derivative w.r.t.volumetric liquid water in the lower layer

		// compute flux at the bottom of the lyr - th layer
		// NOTE: Kw[lyr] is hydraulic conductivity at the bottom of the j - th  layer
		num = Psi[lyr + 1] - Psi[lyr];
		den = dz[lyr];
		qout[lyr] = -Kw[lyr] * num / den + Kw[lyr];

		// compute flux derivatives
		dqodw1[lyr] = (Kw[lyr] * dPsidw[lyr] - dKwdw1 * num) / den + dKwdw1;
		dqodw2[lyr] = (-Kw[lyr] * dPsidw[lyr + 1] - dKwdw2 * num) / den + dKwdw2;
	}

	// bottom layer

	int blyr = nsoillayers_sim - 1;

	// get the flux from above
	qin[blyr] = qout[blyr - 1];

	// get the derivatives from above
	dqidw0[blyr] = dqodw1[blyr - 1];
	dqidw1[blyr] = dqodw2[blyr - 1];

	// select type of boundary condition

	// Bedrock with no drainage
	qout[blyr] = 0.0;

	// compute the flux derivative
	dqodw1[blyr] = 0.0;

	// free drainage flux boundary condition
	if (!ifbedrock) {

		// free drainage
		qout[blyr] = Kw[blyr];

		// compute the flux derivative
		dqodw1[blyr] = dKwdw[blyr] / vws_liq[blyr];
	}
}

void Soil::compute_RHS_moisture_form(double wflux_root[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX],
	double qout[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX], double rmx[NSOILLAYERS_MAX], int nsoillayers_sim) {


	// Calculate RHS of moisture - based form of Richards equation

	// wflux_root - vertically distributed transpiration sink (mm H2O / s) (+= to atm)

	// Compute RHS of moisture - based form of Richards equation
	for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {

		double fluxNet = qin[lyr] - qout[lyr] - wflux_root[lyr];
		rmx[lyr] = -fluxNet * dt_dz[lyr];
	}
}

void Soil::compute_LHS_moisture_form(double dt_dz[NSOILLAYERS_MAX], double dqidw0[NSOILLAYERS_MAX], 
	double dqidw1[NSOILLAYERS_MAX], double dqodw1[NSOILLAYERS_MAX], double dqodw2[NSOILLAYERS_MAX],
	double amx[NSOILLAYERS_MAX], double bmx[NSOILLAYERS_MAX], double cmx[NSOILLAYERS_MAX], int nsoillayers_sim) {

	// Calculate LHS of moisture - based form of Richards equation
	// LHS is a tridiagonal matrix

	// amx - left off diagonal of tridiagonal matrix
	// bmx - diagonal column for tridiagonal matrix
	// cmx - right off diagonal tridiagonal matrix


	// Compute LHS of moisture - based form of Richards equation

	// Top soil layer
	int tlyr = 0;
	amx[tlyr] = 0.0;
	bmx[tlyr] = -1.0 - (-dqidw1[tlyr] + dqodw1[tlyr]) * dt_dz[tlyr];
	cmx[tlyr] = -dqodw2[tlyr] * dt_dz[tlyr];

	// Interior soil layers
	for (int lyr = 1; lyr < nsoillayers_sim - 1; lyr++) {
		amx[lyr] = dqidw0[lyr] * dt_dz[lyr];
		bmx[lyr] = -1.0 - (-dqidw1[lyr] + dqodw1[lyr]) * dt_dz[lyr];
		cmx[lyr] = -dqodw2[lyr] * dt_dz[lyr];
	}

	int blyr = nsoillayers_sim - 1;
	amx[blyr] = dqidw0[blyr] * dt_dz[blyr];
	bmx[blyr] = -1.0 - (-dqidw1[blyr] + dqodw1[blyr]) * dt_dz[blyr];
	cmx[blyr] = 0.0;
}

void Soil::TridiagonalCol(double a[NSOILLAYERS_MAX], double b[NSOILLAYERS_MAX], double c[NSOILLAYERS_MAX], 
	double r[NSOILLAYERS_MAX], double u[NSOILLAYERS_MAX], int nsoillayers_sim) {

	// Tridiagonal matrix solution

	// a - left off diagonal of tridiagonal matrix
	// b - diagonal column for tridiagonal matrix
	// c - right off diagonal tridiagonal matrix
	// r - forcing term of tridiagonal matrix
	// u - solution

	// Solve the matrix

	int tlyr = 0;

	double bet = b[tlyr];
	double gam[NSOILLAYERS_MAX];

	for (int lyr = tlyr; lyr < nsoillayers_sim; lyr++) {
		if (lyr > tlyr) {
			gam[lyr] = c[lyr - 1] / bet;
			bet = b[lyr] - a[lyr] * gam[lyr];
			u[lyr] = (r[lyr] - a[lyr] * u[lyr - 1]) / bet;
		}
		else {
			u[lyr] = r[lyr] / bet;
		}
	}

	for (int lyr = nsoillayers_sim - 2; lyr >= tlyr; lyr--) {
		u[lyr] = u[lyr] - gam[lyr + 1] * u[lyr + 1];
	}
}

void Soil::check_error(double dwat[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX],
	double qin[NSOILLAYERS_MAX], double qout[NSOILLAYERS_MAX], double wflux_root[NSOILLAYERS_MAX], int nsoillayers_sim,
	double dtmin, double xTolerUpper, double& dtsub, double& errorMax, bool& calcsub) {

	calcsub = false;
	double neg_wmass = false;

	double fluxNet0[NSOILLAYERS_MAX];
	double fluxNet1[NSOILLAYERS_MAX];

	// Estimate errors for each layer
	for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {

		// Convert iteration increment to a net flux (mm s-1)
		fluxNet0[lyr] = dwat[lyr] / dt_dz[lyr];

		// Compute the net flux at the start of the sub - step
		fluxNet1[lyr] = qin[lyr] - qout[lyr] - wflux_root[lyr];

		// Check for negative water
		if (soillayers[lyr].wmass + dwat[lyr] * dz[lyr] < 0.0) {
			neg_wmass = true;
		}
	}

	errorMax = 0.0;

	// Compute absolute errors
	for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {
		if (fabs(fluxNet1[lyr] - fluxNet0[lyr]) * dtsub * 0.5 > errorMax)
			errorMax = fabs(fluxNet1[lyr] - fluxNet0[lyr]) * dtsub * 0.5;
	}

	// Check if the error is above the upper tolerance 
	// or we have a soil layer with negative wmass
	if ((errorMax > xTolerUpper || neg_wmass)) {

		if (dtsub > dtmin) {
			// Halve the length of the sub - step
			if (dtsub / 2.0 > dtmin) {
				dtsub /= 2.0;
			}
			// Go to the smallest timestep
			else {
				dtsub = dtmin;
			}
			// Substep rejected; try again
			calcsub = true;
		}
		// We can't solve the issue with a smaller allowed timestep
		else {

			// Total water (currently, positiv, and after)
			double wmass_cur_total = 0.0;
			double wmass_pos_total = 0.0;
			double wmass_aft_total = 0.0;

			for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {
				wmass_cur_total += soillayers[lyr].wmass;
				wmass_pos_total += max(0.0, soillayers[lyr].wmass + dwat[lyr] * dz[lyr]);
				wmass_aft_total += soillayers[lyr].wmass + dwat[lyr] * dz[lyr];
			}

			if (!negligible(wmass_cur_total) && !negligible(wmass_pos_total)) {

				// Fraction water left after movement and loss through evapotranspiration
				double frac_wmass = wmass_aft_total / wmass_cur_total;

				// Distribute change
				for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {
					// Fraction in each layer if excuding negative values
					double frac_pos = max(0.0, soillayers[lyr].wmass + dwat[lyr] * dz[lyr]) / wmass_pos_total;
					
					// Current fraction in each layer
					double frac_cur = soillayers[lyr].wmass / wmass_cur_total;

					// Average fraction
					double frac = (frac_cur + frac_pos) / 2.0;

					// Update water change
					dwat[lyr] = (frac * wmass_aft_total - soillayers[lyr].wmass) / dz[lyr];
				}
			}
			else {
				for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {
					if (soillayers[lyr].wmass + dwat[lyr] * dz[lyr] < 0.0) {
						// Negative water mass
						dprintf("ERROR!!! Negative wmass! Year %d Day %d lyr %d wmass %g\n",
							date.year, date.day, lyr, soillayers[lyr].wmass + dwat[lyr] * dz[lyr]);

						dwat[lyr] = -soillayers[lyr].wmass / dz[lyr];
					}
				}
			}

			// Solve water moved as good as time step resolution allows
			calcsub = false;
		}
	}
	else {
		// Everything is fine
		calcsub = false;
	}
}

void Soil::update_wtp() {

	wtp = -depth;

	bool thres = true;

	for (int lyr = nsoillayers - 1; lyr >= 0; lyr--) {

		if (thres) {

			Soillayer& soillayer = soillayers[lyr];

			// Fraction of maximum water capacity
			double frac = min(1.0, max(0.0, (soillayer.wmass + soillayer.imass) / soillayer.wmass_max()));

			if (frac < (soillayer.wp + soillayer.awc / 2.0) / soillayer.wmass_max()) {
				thres = false;
			}
			else {
				wtp += pow(frac, 4) * soillayer.thickness;
			}			
		}
	}

	wtd = -wtp;
}


// Initial percolation
/** Initial percolation for two layer soil scheme
 * and for peat and wetlands
 */
void Soil::initial_percolation() {

	wmass_wetland_added = 0.0;

	if (percolate) {

		if (hydrology == GERTEN) {

			// Update water content in upper soil layer with rain and snowmelt
			// and adjust wcont_evap

			Soillayer& surflayer = soillayers[0];

			surflayer.wmass += rain_melt;

			if (surflayer.wmass > surflayer.fc()) {
				rain_melt = surflayer.wmass - surflayer.fc();
				surflayer.wmass = surflayer.fc();
			}
			else {
				rain_melt = 0.0;
			}

			wcont_evap = surflayer.wcont();
		}
		else if (patch.stand.is_highlatitude_peatland_stand()) {

			// Peatland soils
			// Distribute the water in the acrotelm in proportion to the capacity

			// water that can still be added to each soil layer (kg m-2)
			double wmass_pot_layer[NSOILLAYERS_MAX];

			double imass_fraction = 0.0;
			double imass_impedance = 1.0;

			double imass_tot = 0.0;
			double wsats_tot = 0.0;

			// Total potential of the soil to store water, from the current content up to saturation
			double wmass_pot_total = 0.0;

			for (int lyr = 0; lyr < nacrotelm; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				wmass_pot_layer[lyr] = max(0.0, soillayer.wsats - (soillayer.wmass + soillayer.imass));

				wsats_tot += soillayer.wsats;
				imass_tot += soillayer.imass;
				wmass_pot_total += wmass_pot_layer[lyr];
			}

			// Potential for infiltration exist
			if (wmass_pot_total > 0.0) {

				imass_fraction = imass_tot / wsats_tot;

				// Impedance - this is the fraction of water that can infiltrate initially 
				if (IMASS_IMPEDANCE && can_freeze()) {
					imass_impedance = ice_imp(imass_fraction);
				}
				else {
					imass_impedance = 1.0;
				}

				double wmass_in = 0.0;

				// Assume the amount of water available for initial infiltration is limited by the ice content in the top layer
				double wmass_for_infiltration = imass_impedance * rain_melt;

				// Only infiltrate what we can. The rest remains in rain_melt.
				if (wmass_for_infiltration >= wmass_pot_total) {
					wmass_in = wmass_pot_total;
					rain_melt -= wmass_pot_total;
				}
				else {
					// All water can infiltrate
					wmass_in = wmass_for_infiltration;
					rain_melt -= wmass_for_infiltration;
				}

				// Now add the rain_melt to each layer in proportion to the capacity if total_potential > 0.0 mm
				for (int lyr = 0; lyr < nacrotelm; lyr++) {
					soillayers[lyr].wmass += wmass_in * (wmass_pot_layer[lyr] / wmass_pot_total);
				}
			}
		}
		else if (hydrology == MILLER && !patch.stand.is_true_wetland_stand()) {

			// NEITHER PEATLAND NOR WETLAND SOILS
			// Distribute the water in the upper 50cm in proportion to the capacity:

			// Water that can still be added to each soil layer (mm)
			double wmass_pot_layer[NSOILLAYERS_UPPER] = { 0.0 };
			double wmass_pot_total = 0.0;

			for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				wmass_pot_layer[lyr] = max(0.0, soillayer.fc() - (soillayer.wmass + soillayer.imass));
				wmass_pot_total += wmass_pot_layer[lyr];
			}

			double water_in = 0.0;

			// Only infiltrate what we can. The rest remains in rain_melt.
			if (rain_melt >= wmass_pot_total) {
				water_in = wmass_pot_total;
				rain_melt -= wmass_pot_total;
			}
			else {
				// Because ALL the water can infiltrate
				water_in = rain_melt;
				rain_melt = 0.0;
			}

			// There is a potential to add water
			if (wmass_pot_total > 0.0) {

				// Now add the rain_melt in proportion to the capacity if total_potential > 0.0 mm
				for (int lyr = 0; lyr < NSOILLAYERS_UPPER; lyr++) {

					soillayers[lyr].wmass += water_in * (wmass_pot_layer[lyr] / wmass_pot_total);
				}
			}
		}
	}

	if (patch.stand.is_true_wetland_stand()) {

		// Mineral wetlands

		// Water mass added to saturate soil
		double wmass_added = 0.0;

		// Saturate the soil if it is the wet season
		if (is_wet_season()) {
			for (int lyr = 0; lyr < nsoillayers; lyr++) {

				Soillayer& soillayer = soillayers[lyr];

				double wmass_input = soillayer.wsats - (soillayer.wmass + soillayer.imass);

				// Add water to the layer
				if (soillayer.temp >= 0.0) {
					soillayer.wmass += wmass_input;
				}
				else {
					soillayer.imass += wmass_input;
				}
				wmass_added += wmass_input;
			}
		}

		// Record the water added to this wetland today
		wmass_wetland_added = max(0.0, wmass_added - rain_melt);

		// Subtract added water from rain_melt
		rain_melt -= min(rain_melt, wmass_added);
	}
}

/// Calculate phase change of surface standing water and ice 
/** Calculate phase change of surface standing water and ice
 */
void SSW::phase_change() {

	// Initial water and ice content
	double wtotal_initial = wmass + imass;

	// Surface standing water depth
	double depth = wmass + imass;

	// Update thermal properties for standing water
	update_thermal_properties();

	// Should ice thaw
	bool ifThaw = false;
	// Should water freeze
	bool ifFreeze = false;
	// Fraction of water which freezes
	double Ffreez = 0.0;
	// Fraction of ice which thaws
	double Fthaw = 0.0;
	// Addition thickness to top snow layer to accound for radiation energy
	double ssw_depth = depth + top_thickness;
	// Fraction of water
	double Fwmass = wmass / ssw_depth;
	// Fraction of ice
	double Fimass = imass / ssw_depth;
	// Energy per volume [J m-3]
	double energy;
	// Temperature difference (current - previous)
	double delta_T;
	// Warning error
	const double EPS_warning = 1.0e-10;
	// Fail error
	const double EPS_fail = 1.0e-6;

	// Should water freeze
	if (temp < 0.0 && Fwmass > 0.0) {
		ifFreeze = true;
	}
	// Should ice thaw
	else if (temp > 0.0 && Fimass > 0.0) {
		ifThaw = true;
	}

	double temp_old = 0.0;

	// Temperature difference from previous phase change
	delta_T = temp - temp_old;

	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy from full standing water to freeze water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / Lheat;

		// Can't freeze all
		if (Fwmass > Ffreez) {

			// Update water and ice content
			wmass = (Fwmass - Ffreez) * ssw_depth;
			imass = (Fimass + Ffreez) * ssw_depth;

			// Update standing water temperature
			temp = temp_old;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			imass += wmass;
			wmass = 0.0;

			// Energy used for phase change and to keep standing water at initial temperature [J m-3]
			double energy_loss = Fwmass * Lheat;

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update standing water temperature
			temp = temp_old - energy_remain / Cp_ice;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from standing water to thaw water [J m-3]
		energy = Ch * fabs(delta_T);

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = energy / Lheat;

		// Can't thaw all
		if (Fimass > Fthaw) {

			// Update water and ice content
			wmass = (Fwmass + Fthaw) * ssw_depth;
			imass = (Fimass - Fthaw) * ssw_depth;

			// Standing water temperature stays at previous timesteps temperature
			temp = temp_old;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass += imass;
			imass = 0.0;

			// Energy used for phase change and to keep standing water at initial temperature [J m-3]
			double energy_loss = Fimass * Lheat;

			// Calculate remaining energy to increase the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update standing water temperature
			temp = temp_old + energy_remain / Cp_water;
		}
	} // End of Thawing

	// Final water
	double wtotal_final = wmass + imass;

	// Water balance checks
	if ((fabs(wtotal_initial - wtotal_final)) > EPS_warning) {
		if ((fabs(wtotal_initial - wtotal_final)) > EPS_fail) {
			fail("Fail! Year %d Day %d. Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n",
				date.year, date.day, wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
		}
		else {
			dprintf("WARNING! Year %d Day %d. Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n",
				date.year, date.day, wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
		}
	}

	// Update thermal properties for standing water
	update_thermal_properties();

	// Ice despite positive temperature
	if (temp > 0.0 && imass > 0.0) {
		dprintf("ERROR: Phase Change Stand Year %d Day %d Temp %g imass %g\n",
			date.year, date.day, temp, imass);
	}

	// No negative wmass
	if (wmass < -EPS_warning) {
		fail("ERROR: Negative wmass at end of ssw phase_change! Year %d Day %d depth %g delta_T %g wmass %g imass %g\n",
			date.year, date.day, depth, delta_T, wmass, imass);
	}

	// No negative imass
	if (imass < -EPS_warning) {
		fail("ERROR: Negative imass at end of ssw phase_change! Year %d Day %d depth %g delta_T %g wmass %g imass %g\n",
			date.year, date.day, depth, delta_T, wmass, imass);
	}
}

/// Calculate phase change of water and ice for each layer
/** Calculate phase change of water and ice for each layer
 *
 * To know how much ice is formed or melted at each time step,
 * The temperature change from the previous to this time step is first 
 * calculated. If not all of the available water is frozen, the
 * temperature will be held at 0 deg until everything is frozen (and
 * vice versa for the melting). The amount of heat released when 
 * freezing or absorbed when thawing is to be found. To do this
 * the total heat capacity of the soil is multiplied with delta_T. Since
 * we are working in m-3 here, we don't need to multiply by the volume
 * or mass, because we want to get the fraction of water freezing or
 * ice thawing, not the absolute value. 
 *
 * Determine where to start:
 * if there's standing water, include that layer for the phase change 
 * calculation, if not, start with the top soil layer
 */
void Soil::phase_change() {

	// Snow phase change

	// Phase change between water and ice and adjust snow layer temperatures accordingly
	snow.phase_change();

	// Standing water phase change
	ssw.phase_change();

	// Soil phase change

	// No phase change for years before first freeze year
	if (date.year >= FIRST_FREEZE_YEAR) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].phase_change();
		}
	}
}

/// Calculate thaw depth
/** Calculate thaw depth
 * Define the thaw depth as the 0 degree C isotherm (Muller 1947, Burn 1998), 
 * which equals the depth, until which liquid water exists continuously from the top down.
 */
void Soil::thaw_depth() {

	// Daily thawing depth, where at least some water is found, but ice can be present [mm]
	thaw = 0.0;

	// Add soil layer where all ice is thawn
	int lyr = 0;
	while (lyr < nsoillayers && soillayers[lyr].temp > 0.0) {
		thaw += soillayers[lyr].thickness;
		lyr++;
	}
	// If all soil layer where thawn then add padding layer with temperature above 0 degree C
	if (lyr == nsoillayers) {
		while (lyr < nsoillayers + NPADLAYERS && pad.padlayers[lyr - nsoillayers].temp > 0.0) {
			thaw += pad.padlayers[lyr - nsoillayers].thickness;
			lyr++;
		}
	}

	// If still among soil layers then assuming that water in next soil layer is at the top of the layer and ice below
	if (lyr < nsoillayers) {
		Soillayer& soillayer = soillayers[lyr];
		if (!negligible(soillayer.wmass + soillayer.imass)) {
			thaw += soillayer.wmass / (soillayer.wmass + soillayer.imass) * soillayer.thickness;
		}
	}

	// Remove very small negative values
	thaw = max(0.0, thaw);

	// Set maximum thaw depth this year
	if (thaw > maxthaw || date.day == 0)
		maxthaw = thaw;

	// Calculate mean monthly thaw depth
	mthaw[date.month] += thaw / date.ndaymonth[date.month];

	// Equilsom
	if (savesomcent && date.islastday) {
		mald_mean[date.month] = min(depth, max(0.0, mthaw[date.month])) * M_PER_MM;
	}

	if (date.islastmonth && date.islastday) {
		maxthaw_10.add(maxthaw);
	}
}

/// True if phase changes allowed, false otherwise
/** True if phase changes allowed, false otherwise
 */
bool Soil::can_freeze() const {
	return ifallowphasechanges;
}

/// Calculate temperature profile through snow, soil and padding layers
/** Calculate temperature profile through snow, soil and padding layers
 */
void Soil::calc_temp() {

	double temp[NLAYERS_MAX + 1], 
		dz[NLAYERS_MAX + 1],
		ki[NLAYERS_MAX + 1],
		ci[NLAYERS_MAX + 1],
		di[NLAYERS_MAX + 1]; // diffusivity of the soil layers [mm2 / day]

	// Air temperature 
	double airT = patch.get_climate().temp;

	// Number of sub timesteps to be calculated
	int stimesteps = 6;

	if (energy == BELDA) {
		stimesteps = 24;
	}

	for (int ts = 0; ts < stimesteps; ts++) {

		int nsurflayers = 0;
		bool combine_surf_soil = false;

		// Layer above soil exists
		if (!negligible(snow.depth() + ssw.water_and_ice_content())) {

			// If a snowpack and standing water is present then soak snowpack
			if (!negligible(snow.depth()) && !negligible(ssw.water_and_ice_content())) {
				snow.soak_snowpack(ssw.wmass, ssw.imass, ssw.temp);
			}

			// Layer above soil is thinner that 50% of top soil layer
			// Should be combined with top soil layer for computational stability
			if (snow.depth() + ssw.water_and_ice_content() <= soillayers[0].thickness * 0.5) {

				// Combine additional surface layers with top soil layer
				combine_surf_soil = true;

				// We have no surface layer
				nsurflayers = 0;
			}
			else {

				// We have a minimum of one surface layer
				nsurflayers = max(1, snow.nsnowlayers);
			}
		}

		// Update heat capacity and thermal conductivity for snow, standing water, soil and padding layers
		update_thermal_properties();

		int nairlayers = 1;
		int nairsurflayers = nairlayers + nsurflayers;
		int nairsurfsoillayers = nairlayers + nsurflayers + nsoillayers;
		int nairsurfsoilpadlayers = nairlayers + nsurflayers + nsoillayers + NPADLAYERS;

		// Air layer
		if (energy == CNSTEP_FULL) {
			dz[0] = AIR_THICKNESS_CNSTEP_FULL;
		}
		else if (energy == CNSTEP) {
			dz[0] = AIR_THICKNESS_CNSTEP;
		}
		else if (energy == BELDA) {
			// Not used
			dz[0] = 0.0;
		}
		
		ci[0] = Cp_air * M3_PER_MM3;
		di[0] = Kair / Cp_air * SECS_PER_DAY * MM2_PER_M2;	// mm2 day-1
		ki[0] = di[0] * ci[0];
		temp[0] = airT;

		// Add additional layers
		if (nsurflayers) {
			if (snow.depth()) {
				int slyr = snow.nsnowlayers - 1;
				for (int lyr = nairlayers; lyr < nairsurflayers; lyr++) {

					// Snow layer
					Snowlayer& snowlayer = snow.snowlayers[slyr];

					dz[lyr] = snowlayer.thickness;
					ci[lyr] = snowlayer.Ch / MM3_PER_M3; // m-3 to mm-3
					ki[lyr] = snowlayer.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
					di[lyr] = ki[lyr] / ci[lyr]; // mm2 day-1
					temp[lyr] = snowlayer.temp;
					
					slyr--;
				}
			}
			else {

				// Standing water layers
				int lyr = nairsurflayers - 1;

				dz[lyr] = ssw.water_and_ice_content();
				ci[lyr] = ssw.Ch / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = ssw.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				di[lyr] = ki[lyr] / ci[lyr]; // mm2 day-1
				temp[lyr] = ssw.temp;
			}
		}

		// Soil layers
		for (int lyr = nairsurflayers; lyr < nairsurfsoillayers; lyr++) {

			// Soil layer
			Soillayer& soillayer = soillayers[lyr - nairsurflayers];

			// Combine snow and standing water with top soil layer
			if (lyr == nairsurflayers && combine_surf_soil) {
				
				// Additional layer to add to top soil layer
				double thickness_surf, Ch_surf, Kh_surf, temp_surf;
				if (snow.depth()) {
					Snowlayer& snowlayer = snow.snowlayers[0];
					thickness_surf = snowlayer.thickness + ssw.water_and_ice_content();
					Ch_surf = (snowlayer.Ch * snowlayer.thickness + ssw.Ch * ssw.water_and_ice_content()) / thickness_surf;
					Kh_surf = (snowlayer.Kh * snowlayer.thickness + ssw.Kh * ssw.water_and_ice_content()) / thickness_surf;
					temp_surf = (snowlayer.temp * snowlayer.thickness * snowlayer.Ch + ssw.temp * ssw.water_and_ice_content() * ssw.Ch) /
						(snowlayer.thickness * snowlayer.Ch + ssw.water_and_ice_content() * ssw.Ch);
				}
				else {
					thickness_surf = ssw.water_and_ice_content();
					Ch_surf = ssw.Ch;
					Kh_surf = ssw.Kh;
					temp_surf = ssw.temp;
				}
				
				dz[lyr] = (soillayer.thickness + thickness_surf);
				ci[lyr] = (soillayer.Ch * soillayer.thickness + Ch_surf * thickness_surf) / dz[lyr] / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = (soillayer.Kh * soillayer.thickness + Kh_surf * thickness_surf) / dz[lyr] * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				di[lyr] = ki[lyr] / ci[lyr]; // mm2 day-1
				temp[lyr] = (soillayer.temp * soillayer.thickness * soillayer.Ch + temp_surf * thickness_surf * Ch_surf) /
					(soillayer.thickness * soillayer.Ch + thickness_surf * Ch_surf);
			}
			else {
				dz[lyr] = soillayer.thickness;
				ci[lyr] = soillayer.Ch / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = soillayer.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1 
				di[lyr] = ki[lyr] / ci[lyr]; // mm2 day-1
				temp[lyr] = soillayer.temp;
			}
		}

		// Padding layers
		for (int lyr = nairsurfsoillayers; lyr < nairsurfsoilpadlayers; lyr++) {

			// Padding layer
			Padlayer& pad_lyr = pad.padlayers[lyr - nairsurfsoillayers];

			dz[lyr] = pad_lyr.thickness;
			ci[lyr] = pad_lyr.Ch / MM3_PER_M3; // m-3 to mm-3
			ki[lyr] = pad_lyr.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
			di[lyr] = ki[lyr] / ci[lyr]; // mm2 day-1
			temp[lyr] = pad_lyr.temp;
		}

		// Calculated temperatures using the Crank-Nicholson timestepper algorithm 
		// for temperature diffusion equation
		if (energy == CNSTEP_FULL) {
			cnstep_full(airT, 1.0 / (double)stimesteps, nairsurfsoilpadlayers, temp, dz, ki, ci);
		}
		else if (energy == CNSTEP) {
			cnstep(airT, 1.0 / (double)stimesteps, nairsurfsoilpadlayers, temp, dz, di);
		}
		else if (energy == BELDA) {
			belda(airT, 1.0 / (double)stimesteps, nairsurfsoilpadlayers, temp, dz, ki, ci);
		}

		// Update temperatures

		// Set new snow layer temperatures (0 == bottom layer, nsnowlayers - 1 == top snow layer)
		if (!negligible(snow.depth())) {
			if (combine_surf_soil) {
				snow.snowlayers[0].temp = temp[1];
			}
			else {
				int slyr = snow.nsnowlayers - 1;
				for (int lyr = nairlayers; lyr < nairsurflayers; lyr++) {
					snow.snowlayers[slyr].temp = temp[lyr];
					slyr--;
				}
			}
		}

		// Set new standing water temperatures
		if (!negligible(ssw.water_and_ice_content(), -10)) {
			if (combine_surf_soil) {
				ssw.temp = temp[1];
			}
			else {
				int lyr = max(1, snow.nsnowlayers);
				ssw.temp = temp[lyr];
			}
		}
		else {
			// Very low surface standing water -> send to runoff
			patch.fluxes.report_flux(Fluxes::WRUNOFF_SSW, ssw.water_and_ice_content());
			ssw.wmass = 0.0;
			ssw.imass = 0.0;
		}

		// Set new soil layer temperatures
		for (int lyr = nairsurflayers; lyr < nairsurfsoillayers; lyr++) {
			soillayers[lyr - nairsurflayers].temp = temp[lyr];
		}

		// Set new padding layer temperatures
		for (int lyr = nairsurfsoillayers; lyr < nairsurfsoilpadlayers; lyr++) {
			pad.padlayers[lyr - nairsurfsoillayers].temp = temp[lyr];
		}

		// Phase change between water and ice and adjust snow, standing water, 
		// and soil layer temperatures accordingly
		phase_change();

		// Update snow properties
		snow.balance_snow(patch, rain_melt, rain_melt_temp);
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	// Thaw depth calculation
	thaw_depth();

	// Update bottom padding layers temperature so they faster comes to an equilibrium
	if (date.year < nyear_spinup * 0.5 && date.islastmonth && date.islastday) {
		double avr_temp = 0.0;

		for (int mth = 0; mth < 12; mth++) {
			avr_temp += soillayers[nsoillayers - 1].mtemp[mth];
		}

		avr_temp /= 12.0;

		for (int lyr = NPADLAYERS-2; lyr < NPADLAYERS; lyr++) {
			pad.padlayers[lyr].temp = avr_temp;
		}
	}
}


/// Calculate snow, soil and padding layer temperatures
/** Calculate snow, soil and padding layer temperatures
*/
void Soil::temp() {

	if (hydrology == GERTEN) {
		// Calculate analytic soil temperatures and update dtemp_31 (trend in daily air temperatures for the last 31 days)
		// Analytic soil temperatures can be used to test the accuracy of the numerical scheme in Soil::calc_temp below.
		calc_temp_analytic();
	}
	else {

		// Calculate the soil temperature in snow, soil and padding layers.
		calc_temp();
	}
}

/// Calculate monthly snow properties
/** Calculate monthly snow properties
*/
void Soil::calc_msnow() {

	// Calculate mean monthly snow depth
	snow.mdepth[date.month] += snow.depth() / date.ndaymonth[date.month];

	// Calculate mean monthly snow cover (Dutra et al. 2010 eqn 8)
	snow.mcover[date.month] += min(1.0, snow.depth() / 100.0) / date.ndaymonth[date.month];
	
	// Calculate mean monthly snow temperature
	for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
		if (lyr < snow.nsnowlayers) {
			snow.mtemp[lyr][date.month] += snow.snowlayers[lyr].temp / date.ndaymonth[date.month];
		}
		else {
			snow.mtemp[lyr][date.month] += patch.get_climate().temp / date.ndaymonth[date.month];
		}
	}

	// Sum number of snow days
	if (snow.depth() > 10.0) {
		snow.snowdays++;
	}
}

/// Calculate soillayer root volume
/**
 *
 */
void Soil::calc_root_volume() {

	if (!ifeffvolume) {
		return;
	}

	// Calculate total coarse root volume
	double volume_coarse_root[NSOILLAYERS_MAX] = { 0.0 };

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		if (indiv.pft.lifeform == TREE) {
			double volume_coarse_root_tot = indiv.cmass_wood() * (1.0 - indiv.pft.ag_woody_frac) /
				indiv.pft.wooddens * MM_PER_M;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				volume_coarse_root[lyr] += volume_coarse_root_tot * patch.pft[indiv.pft.id].rootdist[lyr];
			}
		}
		vegetation.nextobj();
	}

	// Add patch pft coarse root litter
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();
		Patchpft& ppft = patch.pft[pft.id];
		if (pft.lifeform == TREE) {
			double volume_coarse_root_tot = (ppft.cmass_litter_sap + ppft.cmass_litter_heart) * (1.0 - pft.ag_woody_frac) /
				pft.wooddens * MM_PER_M;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				volume_coarse_root[lyr] += volume_coarse_root_tot * patch.pft[pft.id].rootdist[lyr];
			}
		}
		pftlist.nextobj();
	}

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Soillayers root volume (maximum value 50% of soil layer volume)

		// Root volume
		soillayer.root_volume = min(0.5 * soillayer.thickness,
			max(0.0, volume_coarse_root[lyr]) / rho_high * MM_PER_M);

		// If organic matter doesn't affect soil properties, then they 
		// need to be updated here if effective volume is updated
		if (!iforganicsoilproperties) {
			double eff_volume_frac = 1.0 - soillayer.root_volume / soillayer.thickness;
			soillayer.wp = soillayer.wp_min * eff_volume_frac * soillayer.thickness;
			soillayer.awc = soillayer.awc_min * eff_volume_frac * soillayer.thickness;
			soillayer.wsats = soillayer.wsats_min * eff_volume_frac * soillayer.thickness;
		}
	}
}

void interp_single_lyr(double preceding_mean, double this_mean, double succeeding_mean,
	int start_lyr, int layers, std::vector<double>& result, double minimum, double maximum) {

	// The values for the beginning and the end of the layer are determined
	// from the average of the two adjacent layer means
	const double first_value = mean(this_mean, preceding_mean);
	const double last_value = mean(this_mean, succeeding_mean);

	// The mid-point value is computed as offset from the mean, so that the
	// average deviation from the mean of first_value and last_value
	// is compensated for.
	// E.g., if the two values at beginning and end of the layer are on average
	// a fraction lower than the layer mean, the mid-layer value is
	// determined as layer mean + the fractio, so that the layer mean is
	// conserved.
	const double average_deviation = mean(first_value - this_mean, last_value - this_mean);

	const double middle_value = this_mean - average_deviation;
	const double half_time = layers / 2.0;

	const double first_slope = (middle_value - first_value) / half_time;
	const double second_slope = (last_value - middle_value) / half_time;

	double sum = 0;
	int lyr = 0;

	// Interpolate the first half
	for (; lyr < layers / 2; lyr++) {
		double current_time = lyr + 0.5; // middle of layer lyr
		result[start_lyr + lyr] = first_value + first_slope * current_time;
		sum += result[start_lyr + lyr];
	}

	// Special case for dealing with the middle layer if layers is odd
	if (layers % 2 == 1) {
		// In this case we can't use the value corresponding to the middle
		// of the layer. We'll simply skip it and calculate it based on
		// whatever the other layer sum up to.
		lyr++;
	}

	// Interpolate the other half
	for (; lyr < layers; lyr++) {
		double current_time = lyr + 0.5; // middle of layer lyr
		result[start_lyr + lyr] = middle_value + second_slope * (current_time - half_time);
		sum += result[start_lyr + lyr];
	}

	if (layers % 2 == 1) {
		// Go back and set the middle value to whatever is needed to
		// conserve the mean
		result[start_lyr + layers / 2] = layers * this_mean - sum;
	}

	// Go through all values and make sure they're all above the minimum
	double added = 0;
	double sum_above = 0;

	for (int lyr = 0; lyr < layers; lyr++) {
		if (result[start_lyr + lyr] < minimum) {
			added += minimum - result[start_lyr + lyr];
			result[start_lyr + lyr] = minimum;
		}
		else {
			sum_above += result[start_lyr + lyr] - minimum;
		}
	}

	double fraction_to_remove = sum_above > 0 ? added / sum_above : 0;

	for (int lyr = 0; lyr < layers; lyr++) {
		if (result[start_lyr + lyr] > minimum) {
			result[start_lyr + lyr] -= fraction_to_remove * (result[start_lyr + lyr] - minimum);

			// Needed (only) due to limited precision in floating point arithmetic
			result[start_lyr + lyr] = max(result[start_lyr + lyr], minimum);
		}
	}

	// Go through all values and make sure they're all below the maximum
	double removed = 0;
	double sum_below = 0;

	for (int lyr = 0; lyr < layers; lyr++) {
		if (result[start_lyr + lyr] > maximum) {
			removed += result[start_lyr + lyr] - maximum;
			result[start_lyr + lyr] = maximum;
		}
		else {
			sum_below += maximum - result[start_lyr + lyr];
		}
	}

	double fraction_to_add = sum_below > 0 ? removed / sum_below : 0;

	for (int lyr = 0; lyr < layers; lyr++) {
		if (result[start_lyr + lyr] < maximum) {
			result[start_lyr + lyr] += fraction_to_add * (maximum - result[start_lyr + lyr]);

			// Needed (only) due to limited precision in floating point arithmetic
			result[start_lyr + lyr] = min(result[start_lyr + lyr], maximum);
		}
	}
}

std::vector<double> Soil::interp_soil_conserve(double fineDz) {

	std::vector<double> thicknesses(nsoillayers);
	std::vector<double> csoilfrac(nsoillayers);
	int nlyr_org = static_cast<int>(std::round(csoilfrac.size()));
	int nlyr = static_cast<int>(std::round(depth / fineDz));
	std::vector<double> result(nlyr);
	double minimum = 0.0;
	double maximum = 1.0;
	int start_lyr = 0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];

		thicknesses[lyr] = soillayer.thickness;
		csoilfrac[lyr] = soillayer.org_frac();
	}
	
	// Loop through the rest of the layers
	for (int lyr = 0; lyr < nlyr_org; lyr++) {

		double c_prev = lyr == 0 ? c_prev = csoilfrac[lyr] : c_prev = csoilfrac[lyr - 1];
		double c_this = csoilfrac[lyr];
		double c_next = lyr == nlyr_org - 1 ? 0.0 : csoilfrac[lyr + 1];

		int layers = static_cast<int>(std::round(thicknesses[lyr] / fineDz));

		interp_single_lyr(c_prev, c_this, c_next,
			start_lyr, layers, result,
			minimum, maximum);

		start_lyr += layers;
	}

	return result;
}

/// Vertical movement of organic matter if C density of soil layer exceed maximum
/** Vertical movement of organic matter if C density of soil layer exceed maximum
*  Downwards movement of SOM if soil layer C exceeds maximum
*/
void Soil::vertical_movement() {

	// Resolution of finersoil layering (mm)
	//const double fineDz = 10.0;
	
	// Conservative interpolation of soil organic matter fractions down the soil depth 
	//std::vector<double> csoilfinefrac = interp_soil_conserve(fineDz);

	double ccont_before = ccont();
	double wcont_before = water_and_ice_content();
	double wadd = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Initial mineral volume
		double min_volume = soillayer.mineral_volume;

		// Calculate soil layer organic carbon masses
		double ccont_org = soillayer.ccont(GHUMUSSOM);

		// Calculate soil layer organic carbon densities
		double rho_org = soillayer.calc_organic_density(GHUMUSSOM);

		// Update organic volume (mm)
		double org_volume = ccont_org / rho_org * MM_PER_M;
		soillayer.organic_volume = org_volume;

		// Effective soil volume
		double eff_volume = soillayer.thickness - soillayer.root_volume - soillayer.coarse_volume;

		// Difference in material and effetive volume
		double diff_volume = min_volume + org_volume - eff_volume;

		// If a difference, solve it by moving material either down or upwards
		if (!negligible(diff_volume)) {

			// For the bottom soil layer only mineral volume is updated
			if (lyr == nsoillayers - 1) {

				soillayer.mineral_volume -= diff_volume;

				// If organic volume exceeds volume of bottom layer, then mineral volume is set to zero
				if (soillayer.mineral_volume < 0.0) {
					soillayer.mineral_volume = 0.0;
				}
				else if (soillayer.mineral_volume > eff_volume - org_volume) {
					soillayer.mineral_volume = eff_volume - org_volume;
				}
			}
			else {

				if (diff_volume > 0.0) {

					// Volume will be moved to layer below
					Soillayer& soillayer_below = soillayers[lyr + 1];

					// Lyr at the interface between coarse layer in the fine scale organic fraction of soil layers // DW_COMMENT
					//int lyr_above = static_cast<int>(std::round(soillayer.depth_bot / fineDz)) - 1;

					// Organic fraction at interface between layers
					//double organic_frac_int = (csoilfinefrac[lyr_above] + csoilfinefrac[lyr_above + 1]) / 2.0; // DW_COMMENT
					double organic_frac_int = (soillayer.org_frac() * soillayer.thickness + 
						receiver_frac * soillayer_below.org_frac() * soillayer_below.thickness) / 
						(soillayer.thickness + receiver_frac * soillayer_below.thickness);

					// Fraction of matter to be moved
					double frac = min(1.0, diff_volume / (min_volume + org_volume));

					// Move mineral nitrogen
					soillayer_below.NH4_mass += soillayer.NH4_mass * frac;
					soillayer_below.NO3_mass += soillayer.NO3_mass * frac;
					soillayer_below.N2O_mass += soillayer.N2O_mass * frac;
					soillayer_below.NO2_mass += soillayer.NO2_mass * frac;
					soillayer_below.NO_mass += soillayer.NO_mass * frac;
					soillayer_below.N2_mass += soillayer.N2_mass * frac;
					soillayer.NH4_mass *= 1.0 - frac;
					soillayer.NO3_mass *= 1.0 - frac;
					soillayer.N2O_mass *= 1.0 - frac;
					soillayer.NO2_mass *= 1.0 - frac;
					soillayer.NO_mass *= 1.0 - frac;
					soillayer.N2_mass *= 1.0 - frac;

					// Check if mineral and som matter is enough to fulfill move
					double frac_min, frac_org;

					// Only organic matter to move
					if (negligible(min_volume)) {
						frac_min = 0.0;
						frac_org = min(1.0, diff_volume / org_volume);
					}
					// Only mineral matter to move
					else if (negligible(org_volume)) {
						frac_min = min(1.0, diff_volume / min_volume);
						frac_org = 0.0;
					}
					// All mineral and organic matter needs to move. Only humus left
					else if (diff_volume > min_volume + org_volume) {
						frac_min = 1.0;
						frac_org = 1.0;
					}
					else {
						// Calculate mineral volume to be moved with interface fraction
						double diff_volume_org = diff_volume * organic_frac_int;
						double diff_volume_min = diff_volume - diff_volume_org;

						if (diff_volume_min <= min_volume && diff_volume_org <= org_volume) {
							frac_min = diff_volume_min / min_volume;
							frac_org = diff_volume_org / org_volume;
						}
						else if (diff_volume_org > org_volume) {
							frac_org = 1.0;
							frac_min = (diff_volume - org_volume) / min_volume;
						}
						else {
							frac_org = (diff_volume - min_volume) / org_volume;
							frac_min = 1.0;
						}
					}

					double cmass_moved = 0.0;

					// Move som organic matter
					for (int p = 0; p < NSOMPOOL; p++) {
						if (p == SURFHUMUS || p == SURFMICRO || p == SOILMETA || p == SOILSTRUCT ||
							p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
							soillayer_below.sompool[p].cmass += soillayer.sompool[p].cmass * frac_org;
							soillayer_below.sompool[p].nmass += soillayer.sompool[p].nmass * frac_org;
							soillayer.sompool[p].cmass *= 1.0 - frac_org;
							soillayer.sompool[p].nmass *= 1.0 - frac_org;
							cmass_moved += soillayer.sompool[p].cmass * frac_org;
						}
					}

					// Update organic matter age
					soillayer_below.update_socage(cmass_moved, soillayer.socage);

					// Move mineral matter
					soillayer_below.mineral_volume += min_volume * frac_min;
					soillayer.mineral_volume *= 1.0 - frac_min;

					// Subtract moved matter
					diff_volume -= min_volume * frac_min + org_volume * frac_org;

					if (diff_volume > 1.0e-12) {
						dprintf("Year %d Day %d Lyr %d ERROR in vertical movement! No more volume to move (%g).\n", date.year, date.day, lyr, diff_volume);
					}

					// Update organic carbon density in layer below
					soillayer_below.rho_org = soillayer_below.calc_organic_density(GHUMUSSOM);

					// Update organic volume in layer below
					soillayer_below.organic_volume = soillayer_below.ccont(GHUMUSSOM) / soillayer_below.rho_org * MM_PER_M;
				}
				// Movement upwards
				else {

					// Layer below
					int lyr_below = lyr + 1;

					while (diff_volume < 0.0 && lyr_below < nsoillayers) {

						// Volume will be exchanged with layer below
						Soillayer& soillayer_below = soillayers[lyr_below];

						// Calculate soil layer organic carbon masses
						double ccont_org = soillayer_below.ccont(GHUMUSSOM);

						// Calculate soil layer organic carbon densities
						double rho_org = soillayer_below.calc_organic_density(GHUMUSSOM);

						// Update organic volume (mm)
						double org_volume = ccont_org / rho_org * MM_PER_M;
						soillayer_below.organic_volume = org_volume;

						// Move mineral, organic and humus matter upwards from the layer below
						double up_volume;

						// Organic and mineral volume below
						double volume_below = soillayer_below.mineral_volume + org_volume;

						if (volume_below > fabs(diff_volume)) {
							up_volume = fabs(diff_volume);
							diff_volume = 0.0;
						}
						else {
							up_volume = volume_below;
							diff_volume += volume_below;
						}

						// Lyr at the interface between coarse layer in the fine scale organic fraction of soil layers
						//int lyr_above = static_cast<int>(std::round(soillayer.depth_bot / fineDz)) - 1;

						// Organic fraction at interface between layers
						//double organic_frac_int = (csoilfinefrac[lyr_above] + csoilfinefrac[lyr_above + 1]) / 2.0;
						//double organic_frac_int = (receiver_frac * soillayer.org_frac() + soillayer_below.org_frac()) / (1.0 + receiver_frac);
						double organic_frac_int = (receiver_frac * soillayer.org_frac() * soillayer.thickness +
							soillayer_below.org_frac() * soillayer_below.thickness) /
							(receiver_frac * soillayer.thickness + soillayer_below.thickness);

						// Fraction of matter to be moved
						double frac = min(1.0, up_volume / volume_below);

						// Move mineral nitrogen
						soillayer.NH4_mass += soillayer_below.NH4_mass * frac;
						soillayer.NO3_mass += soillayer_below.NO3_mass * frac;
						soillayer.N2O_mass += soillayer_below.N2O_mass * frac;
						soillayer.NO2_mass += soillayer_below.NO2_mass * frac;
						soillayer.NO_mass += soillayer_below.NO_mass * frac;
						soillayer.N2_mass += soillayer_below.N2_mass * frac;
						soillayer_below.NH4_mass *= 1.0 - frac;
						soillayer_below.NO3_mass *= 1.0 - frac;
						soillayer_below.N2O_mass *= 1.0 - frac;
						soillayer_below.NO2_mass *= 1.0 - frac;
						soillayer_below.NO_mass *= 1.0 - frac;
						soillayer_below.N2_mass *= 1.0 - frac;

						// Still volume needed to move upwards
						if (!negligible(up_volume)) {

							// Check if mineral and organic matter is enough to fulfill move
							double frac_min, frac_org;
							if (up_volume > soillayer_below.mineral_volume + org_volume) {
								frac_min = 1.0;
								frac_org = 1.0;
							}
							else if (negligible(org_volume)) {
								frac_min = min(1.0, up_volume / soillayer_below.mineral_volume);
								frac_org = 0.0;
							}
							else if (negligible(soillayer_below.mineral_volume)) {
								frac_min = 0.0;
								frac_org = min(1.0, up_volume / org_volume);
							}
							else {
								// Calculate mineral volume to be moved with interface fraction
								double up_volume_org = up_volume * organic_frac_int;
								double up_volume_min = up_volume - up_volume_org;

								if (up_volume_min <= soillayer_below.mineral_volume &&
									up_volume_org <= org_volume) {
									frac_min = up_volume_min / soillayer_below.mineral_volume;
									frac_org = up_volume_org / org_volume;
								}
								else if (up_volume_org > org_volume) {
									frac_org = 1.0;
									frac_min = (up_volume - org_volume) / soillayer_below.mineral_volume;
								}
								else {
									frac_org = (up_volume - soillayer_below.mineral_volume) / org_volume;
									frac_min = 1.0;
								}
							}

							double cmass_moved = 0.0;

							// Move organic matter
							for (int p = 0; p < NSOMPOOL; p++) {
								if (p == SURFHUMUS || p == SURFMICRO || p == SOILMETA || p == SOILSTRUCT || p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
									soillayer.sompool[p].cmass += soillayer_below.sompool[p].cmass * frac_org;
									soillayer.sompool[p].nmass += soillayer_below.sompool[p].nmass * frac_org;
									soillayer_below.sompool[p].cmass *= 1.0 - frac_org;
									soillayer_below.sompool[p].nmass *= 1.0 - frac_org;
									cmass_moved += soillayer_below.sompool[p].cmass * frac_org;
								}
							}

							// Update organic matter age
							soillayer.update_socage(cmass_moved, soillayer_below.socage);

							// Move mineral matter
							soillayer.mineral_volume += soillayer_below.mineral_volume * frac_min;
							soillayer_below.mineral_volume *= 1.0 - frac_min;
						}

						// Update organic carbon density in layer below
						soillayer_below.rho_org = soillayer_below.calc_organic_density(GHUMUSSOM);

						// Update organic volume in layer below
						soillayer_below.organic_volume = soillayer_below.ccont(GHUMUSSOM) / soillayer_below.rho_org * MM_PER_M;

						lyr_below++;
					}

					// If still missing volume add mineral
					if (diff_volume < 0.0) {
						soillayer.mineral_volume += fabs(diff_volume);
					}
				}
			}
		}

		// Update soil layer organic carbon density
		soillayer.rho_org = soillayer.calc_organic_density(GHUMUSSOM);

		// Update organic volume
		soillayer.organic_volume = soillayer.ccont(GHUMUSSOM) / soillayer.rho_org * MM_PER_M;

		// If bottom layer organic volume exceed the total volume of layer
		// then increase depth of layer by 100mm and fill it with mineral soil 
		if (lyr == nsoillayers - 1 && soillayer.mineral_volume == 0 && soillayer.organic_volume > eff_volume) {

			if (hydrology == RICHARDS) {

				// Depth to add
				double add_depth = 100.0;

				// Increase effective volume
				eff_volume += add_depth;

				// Increase soil layer depth
				soillayer.thickness += add_depth;

				// Update soil layer depth middle
				soillayer.depth_mid += add_depth / 2.0;

				// Update soil layer depth bottom
				soillayer.depth_bot += add_depth;

				// Update array of soil layer thicknesses
				soillayerthicks[lyr] += add_depth;

				// Update total soil depth
				depth += add_depth;

				// Coarse fragment volume
				soillayer.coarse_volume = soillayer.coarse_frac * soillayer.thickness;

				// Update soil layer mineral volume
				soillayer.mineral_volume = soillayer.thickness - (soillayer.root_volume + soillayer.coarse_volume + soillayer.organic_volume);

				double wsats_new;
				if (soillayer.acrotelm) {
					wsats_new = soillayer.mineral_volume * acrotelm_por;
				}
				else if (soillayer.catotelm) {
					wsats_new = soillayer.mineral_volume * catotelm_por;
				}
				else {
					wsats_new = soillayer.mineral_volume * soillayer.wsats_min;
				}

				soillayer.wsats += wsats_new;

				if (soillayer.temp > 0.0) {
					soillayer.wmass += wsats_new;
				}
				else {
					soillayer.imass += wsats_new;
				}

				wadd += wsats_new;

				patch.fluxes.report_flux(Fluxes::WADDED, wsats_new);
			}
			else {
				soillayer.organic_volume = eff_volume;
			}
		}

		if (fabs(eff_volume - (soillayer.mineral_volume + soillayer.organic_volume)) > 1.0e-10) {
			dprintf("Year %d Day %d Lyr %d VERTMOVE VOLUME ERROR %g Eff %g Volume %g\n",
				date.year, date.day, lyr, fabs(eff_volume - (soillayer.mineral_volume + soillayer.organic_volume)),
				eff_volume, soillayer.mineral_volume + soillayer.organic_volume);
		}

		soillayer.mineral_volume = max(0.0, eff_volume - soillayer.organic_volume);
		soillayer.organic_volume = max(0.0, eff_volume - soillayer.mineral_volume);
	}

	double ccont_after = ccont();
	double wcont_after = water_and_ice_content();

	if (fabs(ccont_after - ccont_before) > 1.0e-10) {
		dprintf("Year %d ERROR VERTMOVE C before %g after %g Diff %g\n", date.year, ccont_before, ccont_after, ccont_before - ccont_after);
	}

	if (fabs((wcont_after - wadd) - wcont_before) > 1.0e-10) {
		dprintf("Year %d ERROR VERTMOVE W before %g after %g Diff %g\n", date.year, wcont_before, wcont_after, wcont_before - wcont_after);
	}
}

/// Serialize soil variables
/** Serialize soil variables
 */
void Soil::serialize(ArchiveStream& arch) {
	arch& depth
		& nsoillayers
		& nacrotelm
		& ncatotelm
		& snow
		& ssw
		& pad
		& maxthaw_10
		& slope_grade_tan
		& slope_grade_cos;

	for (int i = 0; i < NSOILLAYERS_MAX; i++) {
		arch& soillayerthicks[i];
	}

	for (int i = 0; i < nsoillayers; i++) {
		arch & soillayers[i];
	}
}

void Soil::daily_accounting() {
	if (date.day == 0) {
		annual_accounting();
	}
	for (size_t i = 0; i < nsoillayers; i++) {
		Soillayer &layer = soillayers[i];
		layer.daily_accounting();
	}
}

void Soil::annual_accounting() {

}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// References not related to Arctic and wetland code:
//
// Campbell, G.S.: Simple method for determining unsaturated conductivity from moisture 
//  retention data. Soil Sci. 117: 311�314. 
//  https://doi.org/10.1097/00010694-197406000-00001, 1974
// Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., 
//   and Friedlingstein, P.: An improved representation of physical permafrost dynamics 
//   in the JULES land - surface model, Geosci.Model Dev., 8, 1493 - 1508, 
//   https ://doi.org/10.5194/gmd-8-1493-2015, 2015
// Clapp, R.B., and Hornberger, G.M. Empirical equations for some soil hydraulic properties. 
//  Water Resour. Res. 14: 601�604. https://doi.org/10.1029/WR014i004p00601, 1978
// Fukusako, S. Int J Thermophys (1990) 11: 353. doi:10.1007/BF01133567
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//   and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Lawrence, D. M. and Slater, A.: Incorporating organic soil into a global climate 
//  model, Clim.Dynam., 30, 145�160, doi:10.1007/s00382-007-0278-1, 2008.
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
//
///////////////////////////////////////////////////////////////////////////////////////
//  Arctic and wetland references:
//
// Aerts, R., Verhoeven, J.T.A., & Whigham, D.F. (1999) Plant-mediated controls on nutrient cycling 
//   in temperate fenns and bogs, Ecology, 80 (7), 2170-2181. 
// Best et al. (2011) Geosci. Model Dev., 4, 677-699. www.geosci-model-dev.net/4/677/2011/
//   doi:10.5194/gmd-4-677-2011
// Bonan, G.B. 2002 Ecological Climatology, Cambridge University Press 
// Clein, J. S., and J. P. Schimel, Microbial activity of tundra and taiga soils at sub-zero 
//   temperatures, Soil Biol. Biochem., 27(9), 1231-1234, 1995.
// Cronk, J. K. and Fennessy, M. S.: Wetland Plants: Biology and Ecology, CRC Press LLC, 2001.
// Ekici, A., Chadburn, S., Chaudhary, N., Hajdu, L. H., Marmy, A., Peng, S., Boike, J., Burke, E., 
//   Friend, A. D., Hauck, C., Krinner, G., Langer, M., Miller, P. A., and Beer, C.: Site-level model 
//   intercomparison of high latitude and high altitude soil thermal dynamics in tundra and barren 
//   landscapes, The Cryosphere, 9, 1343-1361, https://doi.org/10.5194/tc-9-1343-2015, 2015.
// Frolking, S., Roulet, N. T., Moore, T. R., Richard, P. J. H., Lavoie, M., and Muller, S. D.: Modeling 
//   northern peatland decomposition and peat accumulation, Ecosystems, 4, 479-498, 2001.
// Frolking, S., Roulet, N. T., Tuittila, E., Bubier, J. L., Quillet, A., Talbot, J., and Richard, 
//   P. J. H.: A new model of Holocene peatland net primary production, decomposition, water balance, 
//   and peat accumulation, Earth Syst. Dynam., 1, 1-21, https://doi.org/10.5194/esd-1-1-2010, 2010.
// Fukusako, S. (1990) Thermophysical properties of ice, snow, and sea ice. Int. J. Thermophys., 11(2):
//	 353-372. doi:10.1007/BF01133567
// Granberg, et al. 1999 A simple model for simulation of water content, soil frost, 
//   and soil temperatures in boreal mixed mires. Water Resour. Res., 35(12), 3771-3782.
// Hillel D., (1982) Introduction to soil physics. Academic Press, San Diego, CA, USA.
// Ise, T., Dunn, A.L, Wofsy, S.C. & Moorcroft, P.R.: High sensitivity of peat decomposition to climate 
//   change through water-table feedback. Nature Geoscience volume 1, pages 763-766 (2008)
// Jaehne, B., Heinz, G., and Dietrich, W.: Measurement of the diffusion coefficients of sparingly 
//   soluble gases in water, J. Geophys. Res., 92, 10767-10776, 1987.
// Koven, C.D, Ringeval, B., Friedlingstein, P., Ciais, P., Cadule, P., Khvorostyanov, D.,
//   Krinner, G., and Tarnocai C. (2011) Permafrost carbon-climate feedbacks accelerate global warming, 
//   PNAS, vol. 108 no. 36, 14769-14774. 
// Lawrence, D. M., and A. G. Slater, 2008: Incorporating organic soil into a global climate model. 
//   Climate Dynamics, 30, 145-160, doi:10.1007/s00382-007-0278-1.
// Ling, F., and Zhang, T. (2006) Sensitivity of ground thermal regime and surface energy fluxes to
//   tundra snow density in northern Alaska. Cold Regions Science and Technology 44 (2006) 121-130
// McGuire, A. D., Christensen, T. R., Hayes, D., Heroult, A., Euskirchen, E., Kimball, J. S., Koven, C., 
//   Lafleur, P., Miller, P. A., Oechel, W., Peylin, P., Williams, M., and Yi, Y.: An assessment of 
//   the carbon balance of Arctic tundra: comparisons among observations, process models, and atmospheric 
//   inversions, Biogeosciences, 9, 3185-3204, https://doi.org/10.5194/bg-9-3185-2012, 2012.
// Potter, C. S., Davidson, E. A., and Verchot, L. V.: Estimation of global biogeochemical controls and 
//   seasonality in soil methane consumption, Chemosphere, 32, 2219-2246, 1996.
// Riera, J. L., Schindler, J. E., and Kratz, T. K.: Seasonal dynamics of carbon dioxide and methane in 
//   two clear-water lakes and two bog lakes in northern Wisconsin, USA, 
//   Can. J. Fish. Aquat. Sci., 56, 265-274, 1999.
// Sander, R.: Compilation of Henry's Law Constants for inorganic and organic species of potential 
//   importance in environmental chemistry, Tech. Rep. Version 3, MPI Mainz, Air Chemistry Department,
//   Max-Planck Institute of Chemistry, 1999.
// Schimel, J. P.: Plant transport and methane production as controls on methane flux from arctic wet 
//   meadow tundra, Biogeochem., 28, 183-200, 1995.
// Smolders, A. J. P., H. B. M. Tomassen, H. W. Pijnappel, L. P. M. Lamers, and J. G. M. Roelofs (2001), 
//   Substrate-derived CO2 is important in the development of Sphagnum spp., New Phytol., 152(2), 325- 332.
// Spahni, R., Wania, R., Neef, L., van Weele, M., Pison, I., Bousquet, P., Frankenberg, C., Foster, 
//   P. N., Joos, F., Prentice, I. C., and van Velthoven, P.: Constraining global methane emissions and 
//   uptake by ecosystems, Biogeosciences, 8, 1643-1665, https://doi.org/10.5194/bg-8-1643-2011, 2011.
// Sturm, M., Holmgren, J., Konig, M., and Morris, K. (1997) The thermal conductivity of seasonal snow. 
//   Journal of Glaciology 43 (143), 26-41.
// Swenson, S.C., Lawrence, D.M., and Lee, H. 2012. Improved Simulation of the Terrestrial Hydrological 
//   Cycle in Permafrost Regions by the Community Land Model. JAMES, 4, M08002. DOI:10.1029/2012MS000165.
// Tang, J., Miller, P. A., Persson, A., Olefeldt, D., Pilesjo, P., Heliasz, M., Jackowicz-Korczynski, 
//   M., Yang, Z., Smith, B., Callaghan, T. V., and Christensen, T. R.: Carbon budget estimation of a 
//   subarctic catchment using a dynamic ecosystem model at high spatial resolution, 
//   Biogeosciences, 12, 2791-2808, doi:10.5194/bg-12-2791-2015, 2015.
// Wania, R., Ross, I., & Prentice, I.C. (2009a) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: I. Evaluation and sensitivity of physical 
//   land surface processes. Global Biogeochemical Cycles, 23, GB3014, doi:10.1029/2008GB003412
// Wania, R., Ross, I., & Prentice, I.C. (2009b) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: II. Evaluation and sensitivity of vegetation 
//   and carbon cycle processes. Global Biogeochemical Cycles, 23, GB015, doi:10.1029/2008GB003413
// Wania, R., Ross, I., & Prentice, I.C. (2010) Implementation and evaluation of a new methane 
//   model within a dynamic global vegetation model: LPJ-WHyMe v1.3.1, Geosci. Model Dev., 3, 565-584.
// Wisser, D., Marchenko, S., Talbot, J., Treat, C., and Frolking, S. (2011) Soil temperature response 
//   to 21st century global warming: the role of and some implications for peat carbon in thawing 
//   permafrost soils in North America, Earth Syst. Dynam., 2, 121-138
// Wolf, A., Callaghan T.V., & Larson K. (2008) Future changes in vegetation and ecosystem 
//   function of the Barents Region. Climatic Change, 87:51-73 DOI 10.1007/s10584-007-9342-4
// Yurova, A., Wolf, A., Sagerfors, J., & Nilsson, M. (2007) Variations in net ecosystem 
//   exchange of carbon dioxide in a boreal mire: Modeling mechanisms linked to water table 
//   position, Journal of Geophysical Research, 112, art. no. G02025, doi:10.1029/2006JG000342.
// Zhang, W., Miller, P.A., Smith, B., Wania, R., Koenigk, T. & Doscher, R., 2013, 
//   Tundra shrubification and tree-line advance amplify arctic climate warming: results from an 
//   individual-based dynamic vegetation model. Environmental Research Letters 8: 034023.
