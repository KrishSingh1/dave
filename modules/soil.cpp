///////////////////////////////////////////////////////////////////////////////////////
/// \file soil.cpp
/// \brief LPJ-GUESS Soil class and related definitions
/// Implementation of member functions of class Soil.
/// The class Soil and its member functions and variables are declared in guess.h
///
/// \author Paul Miller
/// $Date: 2015-12-22 15:22:04 +0100 (Tue, 22 Dec 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "soil.h"
#include "parameters.h"
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
	depth = 0.0;
	whc = 0.0;
	temp = 0.0;
	Kh = 0.0;
	Ch = 0.0;
}

/// Calculates heat capacity of snow
/** Calculates heat capacity of snow
 *  Best et al. 2011, Eq. 22 [J m-3 K-1]
 */
double snowheatcap(double imass, double wmass, double depth) {
	return (Cp_ice * imass + Cp_water * wmass + Cp_air * max(0.0, (depth - wmass - imass))) / depth;
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

	if (negligible(depth * Ch + add_mass * add_Ch))
		dprintf("Year %d Day %d Error snow temp \n",date.year,date.day);

	// New snow layer but haven't updated snow properties yet 
	if (!negligible(depth) && negligible(Ch)) {
		update_thermal_properties();
	}

	// No snow depth nor heat capacity
	if (negligible(depth) && negligible(Ch)) {
		Ch = add_Ch;
	}

	// Add to existing snow layer with updated snow properties
	if (!negligible(depth) && !negligible(Ch)) {
		temp = (temp * depth * Ch + add_temp * add_mass * add_Ch) / (depth * Ch + add_mass * add_Ch);
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

	// Heat capacity of snow - Best et al. 2011, Eq. 22 [J m-3 K-1]
	Ch = snowheatcap(imass, wmass, depth);

	// Thermal conductivity of snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = snowcond(imass, wmass, depth);
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
	whc = max(0.0, min(depth - imass, imass * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - dens) / dens_t)));
}

/** Phase change of snow layer
 */
void Snowlayer::phase_change(bool toplayer) {

	// Initial water and ice content
	double wtotal_initial = wmass + imass;

	// Initial energy
	double energy_initial = Ch * (temp + K2degC);

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
	// Phase change energy consumption [J]
	double phase_change_energy = 0.0;
	// Maximum error
	const double EPS = 1.0e-10;
	// Addition depth to top snow layer to accound for radiation energy
	double lyr_depth = toplayer ? depth + top_depth : depth;

	// Water and ice fraction of full volume
	Fwmass = wmass / lyr_depth;
	Fimass = imass / lyr_depth;

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
		energy = Ch * delta_T;

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		if (Fwmass > Ffreez) {

			// Update water and ice content
			wmass = (Fwmass - Ffreez) * lyr_depth;
			imass = (Fimass + Ffreez) * lyr_depth;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Ffreez * Lheat;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			wmass = (Fwmass - Fwmass) * lyr_depth;
			imass = (Fimass + Fwmass) * lyr_depth;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fwmass * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update soil temperature
			temp = temp_old + energy_remain / Ch;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Fwmass * Lheat;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from full soil layer to thaw water [J m-3]
		energy = Ch * delta_T;

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = -energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		if (Fimass > Fthaw) {

			// Update water and ice content
			wmass = (Fwmass + Fthaw) * lyr_depth;
			imass = (Fimass - Fthaw) * lyr_depth;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fthaw * Lheat;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass = (Fwmass + Fimass) * lyr_depth;
			imass = 0.0;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = -Fimass * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update soil temperature in this layer
			temp = temp_old + energy_remain / Ch;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fimass * Lheat;
		}

	} // End of Thawing

	// Freezing
	if (ifFreeze) {

		// Update density as depth doesn't change
		if (!negligible(depth)) {
			dens = imass * rho_H2O / depth;
		}
	}
	// Thawing
	else if (ifThaw) {
		
		// Update depth as density doesn't change
		depth = imass / (dens / rho_H2O);

		// Update water holding capacity
		update_whc();

		// Update depth ifmore water than the snow can hold.
		// Correct depth will be calculated after temperature calculations are done
		depth += max(0.0, wmass - whc);
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
		fail("ERROR: Water mass balance during freezing or thawing in Snow::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
	}

	// Final energy
	double energy_final = Ch * (temp + K2degC) - phase_change_energy;

	// water balance checks:
	if ((fabs(energy_initial - energy_final)) > energy_initial * EPS && !toplayer) {
		dprintf("Year %d Day %d Energy balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", date.year, date.day, energy_initial, energy_final, energy_initial - energy_final);
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
		& depth
		& whc
		& temp;
}

////////////////////////////////////////////////////////////////////////////////
// SNOW COLUMN

/// Initial state of snow column
/** Initial state of snow column
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

/// Function returns total available N mass in snow column
/** Function returns total available N mass in snow column
 */
double Snow::nmass_avail() {
	return NO3_mass + NH4_mass;
}

/// Function returns total C mass in snow column
/** Function returns total C mass in snow column
 */
double Snow::ccont() {
	return 0.0;
}

/// Function returns total N mass in snow column
/** Function returns total N mass in snow column
 */
double Snow::ncont() {
	return NO3_mass + NH4_mass;
}

/// Function returns total water and ice mass in snow column
/** Function returns total water and ice mass in snow column
 */
double Snow::wcontent() {
	double wcontent = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		wcontent += snowlayer.wmass + snowlayer.imass;
	}
	return wcontent;
}

/// Function returns total water and ice mass in snow column
/** Function returns total water and ice mass in snow column
 */
double Snow::wmass() {
	double water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		water += snowlayer.wmass;
	}
	return water;
}

/// Function returns total water and ice mass in snow column
/** Function returns total water and ice mass in snow column
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
		sd += snowlayers[lyr].depth;
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
	double snowlyr_depth[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_imass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_wmass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_temp[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_Ch[NSNOWLAYERS_MAX] = { 0.0 };

	// Water in snow layer that has no ice after phase change, hence no depth (a full layer has melted)
	double wmass_no_snowlyr = 0.0;

	// Extract information from all snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		Snowlayer& snowlayer = snowlayers[lyr];

		snowlyr_depth[lyr] = snowlayer.depth;
		snowlyr_imass[lyr] = snowlayer.imass;
		if (snowlayer.depth > 0.0) {
			snowlyr_wmass[lyr] = snowlayer.wmass;
		}
		else {
			wmass_no_snowlyr += snowlayer.wmass;
		}

		snowlyr_temp[lyr] = snowlayer.temp;
		snowlyr_Ch[lyr] = snowlayer.Ch;

		snowlayer.depth = snowlayer.imass = snowlayer.wmass = snowlayer.temp = snowlayer.Ch = 0.0;
	}

	// If water in snow layer without depth, then add water to top snowlayer
	if (wmass_no_snowlyr > 0.0) {
		int lyr = nsnowlayers - 1;
		while (wmass_no_snowlyr > 0.0 && lyr >= 0) {
			if (snowlyr_depth[lyr] > 0.0) {
				snowlyr_wmass[lyr] += wmass_no_snowlyr;
				wmass_no_snowlyr = 0.0;
			}
			else {
				lyr--;
			}
		}
	}

	// Potential depth array
	double pot_depth[NSNOWLAYERS_MAX] = { 0.0 };

	// Potential number of snow layers
	int pot_nsnowlayers;

	// Check potential snow layer depth
	potential_depth(snowpack_depth + snowfall_depth, pot_depth, pot_nsnowlayers);

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

		diff_lyr = pot_depth[lyr] - snowlayer.depth;

		lyr_add = 0;
		diff = true;

		// Fill every snow layer to pot_ice from bottom up
		while (diff && lyr_add < nsnowlayers) {

			if (!negligible(snowlyr_depth[lyr_add])) {

				// Frac of lyr_add to add to this snow layer
				frac_move = min(1.0, diff_lyr / snowlyr_depth[lyr_add]);

				// Update snow layer temperature
				if (negligible(snowlayer.depth)) {
					snowlayer.temp = snowlyr_temp[lyr_add];
				}
				else {
					snowlayer.update_temp(snowlyr_temp[lyr_add], frac_move * (snowlyr_imass[lyr_add] + snowlyr_wmass[lyr_add]), snowlyr_Ch[lyr_add]);
				}

				// Add snow
				snowlayer.depth += frac_move * snowlyr_depth[lyr_add];
				snowlayer.imass += frac_move * snowlyr_imass[lyr_add];
				snowlayer.wmass += frac_move * snowlyr_wmass[lyr_add];

				// Calculate new density
				snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.depth;

				// Update heat capacity
				snowlayer.Ch = snowheatcap(snowlayer.imass, snowlayer.wmass, snowlayer.depth);

				// Subtract added snow
				diff_lyr -= frac_move * snowlyr_depth[lyr_add];
				snowlyr_depth[lyr_add] *= (1.0 - frac_move);
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
			frac_move = min(1.0, max(0.0, pot_depth[lyr] - snowlayer.depth) / snowfall_depth);

			// Update snow layer temperature
			if (negligible(snowlayer.depth)) {
				snowlayer.temp = snowfall_temp;
			}
			else {

				snowlayer.update_temp(snowfall_temp, frac_move * snowfall_imass, snowfall_Ch);
			}

			// Add snow
			snowlayer.depth += frac_move * snowfall_depth;
			snowlayer.imass += frac_move * snowfall_imass;

			// Calculate new density
			snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.depth;

			// Update heat capacity
			snowlayer.Ch = snowheatcap(snowlayer.imass, snowlayer.wmass, snowlayer.depth);

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
			snowlayer.depth += snowfall_depth;
			snowlayer.imass += snowfall_imass;

			// Calculate new density
			snowlayer.dens = snowlayer.imass * rho_H2O / snowlayer.depth;

			// Update heat capacity
			snowlayer.Ch = snowheatcap(snowlayer.imass, snowlayer.wmass, snowlayer.depth);

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
		nsnowlayers = (int)snowlayers.size();
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
		double frac = min(1.0, depth_reduction / snowlayer.depth);

		// Reduce reduction
		depth_reduction -= frac * snowlayer.depth;

		// Runoff water
		runoff_water += frac * snowlayer.wmass;

		// Runoff ice
		runoff_ice += frac * snowlayer.imass;

		// Reduce snow
		snowlayer.wmass *= (1.0 - frac);
		snowlayer.imass *= (1.0 - frac);
		snowlayer.depth *= (1.0 - frac);

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

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];

		snowlayer.depth = snowlayer.imass / (snowlayer.dens / rho_H2O);
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
				patch.soilcolumn.update_ssw_temp(snowlayer.temp, added_wmass);

				// Update standing water content (mm)
				patch.soilcolumn.update_ssw(patch.soilcolumn.wmass_ssw, added_wmass);
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
	}

	// Add excess water to standing water
	patch.soilcolumn.update_ssw_temp(ts_melt_temp, ts_melt);

	// Update standing water content (mm)
	patch.soilcolumn.update_ssw(patch.soilcolumn.wmass_ssw, ts_melt);
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

	// If there was water exceeding whc and couldn't be infiltrated to the soil, 
	// then it will fill the snowpack from bottom up

	int lyr = 0;
	while (wmass_ssw > 0.0 && lyr < nsnowlayers) {

		Snowlayer& snowlayer = snowlayers[lyr];

		// Water to add within water holding capacity
		double add_wmass = max(0.0, min(wmass_ssw, snowlayer.depth - (snowlayer.wmass + snowlayer.imass)));

		// Update temperature due to added water
		snowlayer.update_temp(temp_ssw, add_wmass, Cp_water);

		// Add water to snow layer
		snowlayer.wmass += add_wmass;

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

		// Update snow layer depth
		snowlayer.depth += wmass_ssw;

		// Set standing water to zero
		wmass_ssw = 0.0;
	}

	// Added excess ice to the bottom of the snowpack
	if (imass_ssw > 0.0) {
		Snowlayer& snowlayer = snowlayers[0];

		// Update temperature due to added water
		snowlayer.update_temp(temp_ssw, imass_ssw, Cp_ice);

		// Add ice to snow layer
		snowlayer.imass += imass_ssw;

		// Update snow layer depth
		snowlayer.depth += imass_ssw;

		// Set standing ice to zero
		imass_ssw = 0.0;
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
	// Reference snow density [kg m-3]
	double rho_0 = 80.0;

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
				exp((k_s / K2degC) - (k_s / T_k) - (snowlayer.dens / rho_0)) * SECS_PER_DAY;

			// Increase density
			snowlayer.dens += delta_dens;

			// Snow density can not exceed pure ice density
			snowlayer.dens = min(rho_ICE, snowlayer.dens);
		}

		// Add this layer to the above layers for future calculations
		I_above += I_sum;
	}

	// Calculate snow depth after compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		Snowlayer& snowlayer = snowlayers[lyr];
		snowlayer.depth = snowlayer.imass / (snowlayer.dens / rho_H2O);
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
		double snow_volume = snowlayer.depth * M_PER_MM;
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

		// If layer still exist update depth, water content and temperature
		if (!negligible(snowlayer.imass)) {

			// Recalculate depth
			snowlayer.depth = snowlayer.imass / (snowlayer.dens / rho_H2O);

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
/** Soil layer constructor: initialises soilcolumn, depth, profile depth
 * and if soil layer is an acrotelm or catotelm
 */
Soillayer::Soillayer(Soilcolumn* sc, double d, double dp, bool at, bool ct) : soilcolumn(sc), depth(d), depth_profile(dp), acrotelm(at), catotelm(ct) {

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
	mineral_volume = depth;
	organic_volume = 0.0;
	perc_base = 0.0;

	rho_org = rho_high;

	// Organic exponential in Campbell (1974) parameterization for soil 
	// water retention curves from CLM5.0 documentation 
	b_org = min(2.7 + 9.3 * depth_profile / sapric_depth, 12.0);
	// Organic soil matric potential at saturation 
	Psi_max_org = -max(10.3 - 0.2 * depth_profile / sapric_depth, 10.1);
	// Organic volumetric soil water concentration at saturation point (porosity)
	wsats_org = max(0.93 - 0.1 * depth_profile / sapric_depth, 0.83);
	// Volumetric soil moisture concentration at Field Capacity
	// Field capacity generally occurs at a suction of 100 cm
	double Theta_whc = wsats_org * pow(-100.0 / Psi_max_org, -1.0 / b_org);
	// Volumetric soil moisture concentration at wilting point
	// Wilting point typically occurs at a suction of 15000 cm
	wp_org = wsats_org * pow(-15000.0 / Psi_max_org, -1.0 / b_org);
	// Organic available water holding capacity of soil layer (-)
	awc_org = Theta_whc - wp_org;
	// Hydraulic conductivity at saturation for organic matter (mm s-1) (CLM5.0 - Documentation; Hoag and Price 1995)
	Kw_max_org = max(0.28 - 0.2799 * depth_profile / sapric_depth, (*soilcolumn).soiltype.Kw_max_min);

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
	depth_mod = ifdepthdecaymod ? exp(-depth_profile / 10000.0) : 1.0;
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

/// Returns certain properties of soil layer
/* Returns certain properties of soil layer
*/
double Soillayer::get(soilproptype SPT, int m) {

	// Return property
	double prop;

	switch (SPT) {
	case DEPTH:
		prop = depth;
		break;
	case TEMP:
		prop = temp;
		break;
	case MTEMP:
		prop = mtemp[m];
		break;
	case WMASS:
		prop = wmass;
		break;
	case MWMASS:
		prop = mwmass[m];
		break;
	case IMASS:
		prop = imass;
		break;
	case MIMASS:
		prop = mimass[m];
		break;
	case WCONT:
		prop = wcont();
		break;
	case MWCONT:
		prop = mwcont[m];
		break;
	case AWCONT:
		prop = awcont;
		break;
	case WSATS:
		prop = wsats;
		break;
	case AWC:
		prop = awc;
		break;
	case WP:
		prop = wp;
		break;
	case KW_MAX:
		prop = Kw_max;
		break;
	case MDRAIN:
		prop = mdrain[m];
		break;
	case NH_4:
		prop = NH4_mass;
		break;
	case NO_3:
		prop = NO3_mass;
		break;
	case CH_4:
		prop = CH4;
		break;
	case CO_2:
		prop = CO2;
		break;
	case RESPSUM:
		prop = respsum;
		break;
	case WCONTENT:
		prop = wmass + imass;
		break;
	case RUNOFFDRAIN:
		prop = runoff_drain;
		break;
	case RUNOFFBASE:
		prop = runoff_baseflow;
		break;
	case CSLOW:
		prop = cpool_slow;
		break;
	case CFAST:
		prop = cpool_fast;
		break;
	case ARH:
		prop = arespsum;
		break;
	case ANMIN:
		prop = anmin;
		break;
	case ANIMM:
		prop = animmob;
		break;
	case AORGCLEACH:
		prop = aorgCleach;
		break;
	case AORGNLEACH:
		prop = aorgNleach;
		break;
	case AMINLEACH:
		prop = aminleach;
		break;
	case AFLITTERINPUT:
		prop = acmassfinput;
		break;
	case ACLITTERINPUT:
		prop = acmasscinput;
		break;
	case VOL_ORG:
		prop = organic_volume;
		break;
	case VOL_MIN:
		prop = mineral_volume;
		break;
	case VOL_ROOT:
		prop = root_volume;
		break;
	case MODTEMP:
		prop = atemp_mod;
		break;
	case MODMOIS:
		prop = amoist_mod;
		break;
	case MODTOT:
		prop = atot_mod;
		break;
	default:
		// Should never happen
		fail("Unsupported soil property type (%d)\n", SPT);
	}
	return prop;
}


/// Calculates the effective organic density
/** Based on Chadburn et al. 2022
 *	Coarse root litter is not affecting soil properties, just taking up space
 */ 
double Soillayer::calc_organic_density(grouppooltype group) {

	double volume_tot = 0.0;
	double cmass_tot = 0.0;

	for (int p = 0; p < NSOMPOOL; p++) {
		
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
void Soillayer::set_properties(double& wmass_move, double& imass_move, bool must_set) {

	// Only enter first year except if dynamic soil organic C influence soil properties
	if (!must_set && !ifeffvolume && !((*soilcolumn).patch.stand.first_year == date.year || (iforganicsoilproperties && !acrotelm && !catotelm))) { // DW_COMMENT restart years needs to update it then as well!
		return;
	}
	
	// Get soil type
	Soiltype soiltype = (*soilcolumn).soiltype;

	// Mineral particle density set to 2.65 g cm–3 (Bonan 2019)
	double min_density = Mindensity * G_PER_KG * M3_PER_CM3;

	// Analytic temp calc parameters
	thermdiff_0 = soiltype.thermdiff_0;
	thermdiff_15 = soiltype.thermdiff_15;
	thermdiff_100 = soiltype.thermdiff_100;

	// Heat capacity of mineral soil [J m-3 K-1; Bonan 2019]
	Cp_min = soiltype.Cp_min;

	double eff_volume_frac = 1.0 - root_volume / depth;

	// Acrotelm
	if (acrotelm) {
		wsats = acrotelm_por * depth;
		// Field capacity is set at 75% of saturation (Chaudhary 2017)
		awc = (acrotelm_por * 0.75 - peat_wp) * depth;
		wp = peat_wp * depth;
		peat_frac = 1.0 - acrotelm_por;
		bulk_density = peat_frac * rho_org;

		Kw_max = Kw_max_org;
		b = b_org;
		Psi_max = Psi_max_org;
		a = 0.0;

		// Initially saturated
		if ((*soilcolumn).patch.stand.first_year == date.year) {
			double wadded = wsats - (wmass + imass);
			wmass += wadded;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
		}
	}
	// Catotelm
	else if (catotelm) {
		wsats = catotelm_por * depth;
		// Field capacity is set at 75% of saturation (Chaudhary 2017)
		awc = (catotelm_por * 0.75 - peat_wp) * depth;
		wp = peat_wp * depth;
		peat_frac = 1.0 - catotelm_por;
		bulk_density = peat_frac * rho_org;

		Kw_max = Kw_max_org;
		b = b_org;
		Psi_max = Psi_max_org;
		a = 0.0;

		// Initially saturated
		if ((*soilcolumn).patch.stand.first_year == date.year) {
			double wadded = wsats - (wmass + imass);
			wmass += wadded;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
		}
	}
	// Mineral soil and no influence of soil organic C
	else if (!iforganicsoilproperties) {

		// added moved water and ice due to shrinking wsats
		wmass += wmass_move;
		imass += imass_move;
		wmass_move = imass_move = 0.0;

		wsats = soiltype.wsats_min * depth * eff_volume_frac;
		awc = soiltype.awc_min * depth * eff_volume_frac;
		wp = soiltype.wp_min * depth * eff_volume_frac;
		peat_frac = 0.0;
		bulk_density = min_density * (depth - wsats) / depth;

		Kw_max = soiltype.Kw_max_min;
		b = soiltype.b_min;
		Psi_max = soiltype.Psi_max_min;
		a = soiltype.a;

		// Two layer hydrology parameter
		perc_base = soiltype.perc_base;

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
	// Soil organic content influences soil properties
	// Hudson (1994) showed that for each 1-percent increase in soil organic matter, the available water holding 
	// capacity in the soil increased by 3.7 percent. Soil water is held by adhesive and cohesive forces within 
	// the soil and an increase in pore space will lead to an increase in water holding capacity of the soil.
	else {

		// added moved water and ice due to shrinking maximum water content
		wmass += wmass_move;
		imass += imass_move;
		wmass_move = imass_move = 0.0;

		double org_frac = organic_volume / depth;
		double min_frac = mineral_volume / depth;

		if (fabs(organic_volume + mineral_volume + root_volume - depth) > 1.0e-12)
			dprintf("Year %d Day %d ERROR in soil volume (%g) compared to depth (%g)\n", 
				date.year,date.day, organic_volume + mineral_volume + root_volume, depth);

		wsats = (min_frac * soiltype.wsats_min + org_frac * wsats_org) * depth;
		awc = (min_frac * soiltype.awc_min + org_frac * awc_org) * depth;
		wp = (min_frac * soiltype.wp_min + org_frac * wp_org) * depth;
		peat_frac = 0.0;
		double density = min_frac * min_density + org_frac * rho_org * G_PER_KG * M3_PER_CM3;
		bulk_density = density * (depth - wsats) / depth;

		// Move water downwards if water content exceeds maximum water content
		double wmass_max_lyr = wmass_max();
		if (wmass_max_lyr < wmass + imass) {
			double frac = (wmass + imass - wmass_max_lyr) / (wmass + imass);
			wmass_move = wmass * frac;
			imass_move = imass * frac;
			wmass -= wmass_move;
			imass -= imass_move;
		}

		// Effective mineral and organic fractions
		double min_frac_eff = min_frac / eff_volume_frac;
		double org_frac_eff = org_frac / eff_volume_frac;

		// The saturated hydraulic conductivity
		// The soil organic matter fraction is further subdivided into “connected” and “unconnected” 
		// fractions using a result from percolation theory (Stauffer and Aharony 1994, Berkowitz and 
		// Balberg 1992). Assuming that the organic and mineral fractions are randomly distributed 
		// throughout a soil layer, percolation theory predicts that above a threshold value 
		// org_frac = f_threshold, connected flow pathways consisting of organic material only exist
		// and span the soil space. Flow through these pathways interacts only with organic material, 
		// and thus can be described by Kw_max_org (CLM5.0 documentation).
		double perc_norm = 1.101;
		double perc_frac = org_frac_eff >= threshold_frac ? perc_norm * pow(org_frac_eff - threshold_frac, beta_perc) * org_frac_eff : 0.0;
		double uncon_frac = 1.0 - perc_frac;
		double Kw_max_uncon = org_frac_eff < 1.0 ? uncon_frac / ((1.0 - org_frac_eff) / soiltype.Kw_max_min + (org_frac_eff - perc_frac) / Kw_max_org) : 0.0;
		Kw_max = uncon_frac * Kw_max_uncon + (1.0 - uncon_frac) * Kw_max_org;

		// Exponential in Campbell (1974) parameterization for soil water retention curves
		b = (min_frac_eff * soiltype.b_min + org_frac_eff * b_org);

		// Saturated soil water suction (mm)
		Psi_max = (min_frac_eff * soiltype.Psi_max_min + org_frac_eff * Psi_max_org);

		// Soil organic C microorganism collocation factor for decomposition (Yan 2018)
		a = min_frac_eff * soiltype.clay_frac < 0.016 ? 0.0 : (min_frac_eff * soiltype.clay_frac > 0.37 ? 1.0 : 2.8 * min_frac_eff * soiltype.clay_frac - 0.046);

		// Analytic temp calc parameters
		thermdiff_0 = min_frac_eff * soiltype.thermdiff_0 + org_frac_eff * soiltype.thermdiff_0_org;
		thermdiff_15 = min_frac_eff * soiltype.thermdiff_15 + org_frac_eff * soiltype.thermdiff_15_org;
		thermdiff_100 = min_frac_eff * soiltype.thermdiff_100 + org_frac_eff * soiltype.thermdiff_100_org;

		// Two layer hydrology parameter
		perc_base = min_frac_eff * soiltype.perc_base + org_frac_eff * soiltype.perc_base_org;
	}

	// For Gerten and Miller hydrology schemes wmass should never be below wp
	if (hydrology != RICHARDS) {
		double wadded = max(0.0, wp - (wmass+imass));
		wmass += wadded;
		(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
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

	// Initial energy
	double energy_initial = Ch * (temp + K2degC);

	// Should ice thaw
	bool ifThaw = false;
	// Should water freeze
	bool ifFreeze = false;
	// Fraction of water which freezes
	double Ffreez = 0.0;
	// Fraction of ice which thaws
	double Fthaw = 0.0;
	// Fraction of water
	double Fwmass = wmass / depth;
	// Fraction of water super cooled
	double Fwmass_sc = 0.0;
	// Fraction of ice
	double Fimass = imass / depth;
	// Energy per volume [J m-3]
	double energy;
	// Temperature difference (current - previous)
	double delta_T;
	// Phase change energy consumption [J]
	double phase_change_energy = 0.0;
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
				double w_m = 0.05 * (*soilcolumn).soiltype.ssa;
				double S = 0.37;
				double C_1 = 0.041;
				double C = (*soilcolumn).soiltype.ssa * S * C_1 * pow(K2degC - temp, 0.5);
				wmass_sc = wmass_max_lyr * min(max_sc_frac, (C * pow(-temp, -0.5) + w_m - 0.25 * C) / 100.0);
			}

			// Fraction super cooled water
			Fwmass_sc = wmass_sc / depth;

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
			Fwmass_sc = wmass_sc / depth;
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
		energy = Ch * delta_T;

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		// Super cooled water that could freeze
		double Fwmass_freeze = Fwmass - Fwmass_sc;

		// Can't freeze all
		if (Fwmass_freeze > Ffreez) {

			// Update water and ice content
			wmass = (Fwmass - Ffreez) * depth;
			imass = (Fimass + Ffreez) * depth;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Ffreez * Lheat;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			wmass = (Fwmass - Fwmass_freeze) * depth;
			imass = (Fimass + Fwmass_freeze) * depth;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = Fwmass_freeze * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update soil temperature
			temp = temp_old + energy_remain / Ch;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Fwmass_freeze * Lheat;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from full soil layer to thaw water [J m-3]
		energy = Ch * delta_T;

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = -energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		// Maximum ice that could be turned into liquid form
		double Fimass_thaw = temp > 0.0 ? Fimass : min(Fimass, Fwmass_sc - Fwmass);

		// Can't thaw all
		if (Fimass_thaw > Fthaw) {

			// Update water and ice content
			wmass = (Fwmass + Fthaw) * depth;
			imass = (Fimass - Fthaw) * depth;

			// Update thermal properties
			update_thermal_properties();

			// Soil temperature stays at previous timesteps temperature
			temp = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fthaw * Lheat;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass = (Fwmass + Fimass_thaw) * depth;
			imass = (Fimass - Fimass_thaw) * depth;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = -Fimass_thaw * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update soil temperature in this layer
			temp = temp_old + energy_remain / Ch;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fimass_thaw * Lheat;
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

	// Final energy
	double energy_final = Ch * (temp + K2degC) - phase_change_energy;

	// Energy balance checks
	if ((fabs(energy_initial - energy_final)) > energy_initial * EPS) {
		dprintf("Year %d Day %d Depth %g Energy balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n",
			date.year, date.day, depth_profile, energy_initial, energy_final, energy_initial - energy_final);
	}

	// Ice despite positive temperature
	if (temp > 0.0 && imass > 0.0) {
		dprintf("ERROR: Phase Change Year %d Day %d Depth %g Temp %g imass %g\n",
			date.year, date.day, depth_profile, temp, imass);
	}

	// No negative wmass
	if (wmass < -EPS) {
		fail("ERROR: Negative wmass at end of phase_change! Year %d Day %d depth_profile %g depth %g delta_T %g Phase_energy %g wmass %g imass %g\n",
			date.year, date.day, depth_profile, depth, delta_T, phase_change_energy, wmass, imass);
	}

	// No negative imass
	if (imass < -EPS) {
		fail("ERROR: Negative imass at end of phase_change! Year %d Day %d depth_profile %g depth %g delta_T %g Phase_energy %g wmass %g imass %g\n",
			date.year, date.day, depth_profile, depth, delta_T, phase_change_energy, wmass, imass);
	}
}


/// Water content of soil layer
/** Water content of soil layer
 *  wcont never reports over field capacity even if there is more water
 */
double Soillayer::wcont(double wmass_added) {
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
		(*soilcolumn).patch.stand.is_highlatitude_peatland_stand() || 
		(*soilcolumn).patch.stand.is_true_wetland_stand()) {
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
	return max(0.0, (wsats - (max(0.0, wmass) + imass)) / depth);
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
double Soillayer::ccont(grouppooltype group) {
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
double Soillayer::ncont(grouppooltype group) {
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
			if (date.year >= freenyears && residual < -1.0e-14) {
				dprintf("ERROR year %d and day %d! NH4 and NO3 mass %f is smaller than N substr %f with the residual %g. Adding residual as NO flux\n",
					date.year, date.day, nmass_tot, nmass, residual);
			}
			NO3_mass = 0.0;
			NH4_mass = 0.0;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::NH3_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NO) than available (residual %g)\n",
				date.year, date.day, (*soilcolumn).patch.stand.id, (*soilcolumn).patch.id, depth_profile, residual);
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
			if (date.year >= freenyears && residual < -1.0e-14) {
				dprintf("ERROR year %d and day %d! NH4 mass %f is smaller than N substr %f with the residual %g. Adding residual as NH3 flux\n",
					date.year, date.day, NH4_mass, nmass, residual);
			}
			NH4_mass = 0.0;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::NH3_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NH4) than available (residual %g)\n",
				date.year, date.day, (*soilcolumn).patch.stand.id, (*soilcolumn).patch.id, depth_profile, residual);
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
			if (date.year >= freenyears && residual < -1.0e-14) {
				dprintf("ERROR year %d and day %d! NO3 mass %f is smaller than N substr %f with the residual %g. Adding residual as NO flux\n",
					date.year, date.day, NO3_mass, nmass, residual);
			}
			NO3_mass = 0.0;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::NO_SOIL, -residual);
		}
		// Large negative residual, fail model
		else {
			fail("Fail in year %d and day %d for stand %d and patch %d at depth %g!\nTried to subtract more N (NO3) than available (residual %g)\n",
				date.year, date.day, (*soilcolumn).patch.stand.id, (*soilcolumn).patch.id, depth_profile, residual);
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
	return NMASS_SAT * depth * M_PER_MM;
}

/// Update soil layer temperature due to mass movement
/** Update soil layer temperature due to mass movement
 *  assuming it replaces air
 */
void Soillayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible((depth * Ch - add_mass * Cp_air) + add_mass * add_Ch))
		dprintf("ERROR: Year %d Day %d Updating soil temperature\n", date.year, date.day);

	// Update temperature
	temp = (temp * (depth * Ch - add_mass * Cp_air) + add_temp * add_mass * add_Ch) / ((depth * Ch - add_mass * Cp_air) + add_mass * add_Ch);

	// If temperature becomes positive (negative) with ice (water) present
	// then update thermal properties and perform a phase change for this layer
	if (date.year >= FIRST_FREEZE_YEAR && ((temp > 0.0 && imass > 0.0) || (temp < 0.0 && wmass > 0.0))) {
		// Update thermal properties of layer
		update_thermal_properties(); 

		// Phase change between water and ice and adjust soil layer temperatures accordingly
		phase_change();
	}
}

/// Help function for soil water matric potential calculation
/** Help function for soil water matric potential calculation
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
	double organic_wsats = iforganicsoilproperties ? wsats_org : (*soilcolumn).soiltype.wsats_min;

	// Depths each quantity occupies
	double combined_depth = depth;
	double ice = imass;
	double water = wmass;
	double peat = peat_frac * depth;
	double mineral = (!acrotelm && !catotelm) ? mineral_volume * (1.0 - (*soilcolumn).soiltype.wsats_min) : 0.0;
	double organic = (!acrotelm && !catotelm) ? organic_volume * (1.0 - organic_wsats) : 0.0;
	double coarse_root = (!acrotelm && !catotelm) ? root_volume : 0.0;

	// Soil layer air (mm)
	double air = max(0.0, wsats - (wmass + imass));

	// Update heat capacity (J m-3 K-1)
	Ch = (mineral * Cp_min + organic * Cp_org + ice * Cp_ice + water * Cp_water + peat * Cp_peat + air * Cp_air + coarse_root * Cp_wood) / combined_depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	// Total solids in soil
	Stot = mineral + organic + ice + water + peat + coarse_root;

	if (fabs(1.0 - (Stot + air) / combined_depth) > 1.0e-6)
		fail("Year %d Day %d Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g\n", 
			date.year, date.day, Stot / combined_depth, air / combined_depth, (Stot + air) / combined_depth);

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

/// Serialize soil layer variables
/** Serialize soil layer variables
 */
void Soillayer::serialize(ArchiveStream& arch) {
	arch & depth_profile
		& depth
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

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soilcolumn member functions


/// Initilize soil column variables
/** Initilize soil column variables
 */
void Soilcolumn::initdrivers() {

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
	wmass_ssw = 0.0;
	imass_ssw = 0.0;
	Ch_ssw = 0.0;
	Kh_ssw = 0.0;
	temp_ssw = 0.0;
	wmass_ssw_frac = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mthaw[mth] = 0.0;
		mwtp[mth] = 0.0;
		mald_mean[mth] = 0.0;

		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
	}

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

	firstTempCalc = true;

	// Solve LPJ SOM pools
	solvesom_end = static_cast<int>(0.8 * nyear_spinup);
	solvesom_begin = static_cast<int>(0.7 * nyear_spinup);

	// Solve CENTURY SOM pools
	solvesomcent_beginyr = (int)0.5*freenyears;
	solvesomcent_endyr = (int)freenyears+10;
	savesomcent = false;
}

/// Create soil layers
/** Create soil layers with specific soil depths and thickness
 */
void Soilcolumn::create_soillayers(double soil_layer_depths[NSOILLAYERS_MAX]) {

	double acc_depth = 0.0;
	double depth_profile = 0.0;

	bool acrotelm = false;
	bool catotelm = false;
	nacrotelm = 0;
	ncatotelm = 0;
	thaw = depth;
	maxthaw_10.add(thaw);

	// Loop over number of soil layers to create
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		soillayerthicks[lyr] = soil_layer_depths[lyr];
		acc_depth += soil_layer_depths[lyr];
		soillayerdepths[lyr] = acc_depth;
		depth_profile += lyr == 0 ? soil_layer_depths[lyr] / 2.0 : (soil_layer_depths[lyr - 1] + soil_layer_depths[lyr]) / 2.0;
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
		add_soillayer(soil_layer_depths[lyr], depth_profile, acrotelm, catotelm);
	}
}

/// Add a single soil layer
/** Add a single soil layer to the bottom with a specific depth
 */
void Soilcolumn::add_soillayer(double soil_layer_depth, double soil_layer_depth_profile, bool acrotelm, bool catotelm) {

	// Create soil layer
	Soillayer* soillayer = new Soillayer(this, soil_layer_depth, soil_layer_depth_profile, acrotelm, catotelm);

	// Add soil layer to soillayer vector
	soillayers.push_back((*soillayer));

	delete soillayer;
}

const Climate& Soilcolumn::get_climate() const {
	// All patches within a stand share the same climate
	return patch.stand.get_climate();
}

/// Set slope properties
/** Set slpoe properties
 */
void Soilcolumn::set_slope_properties() {
	// Set fixed slope properties
	// To be able to distinguish between and flat
	// wetland (slope 0) and other landcovers (slope as input)
	// this has to belong to the soilcolumn
	// Slope in radians

	double slope_rad = slope * PI / 180.0;

	// High latitude wetland are on flat ground. Mineral wetlands are seasonal and still 
	// have a slope when Richards hydrology is used during none wet season
	//if (patch.stand.is_highlatitude_peatland_stand()) {
	//	slope_rad = 0.0;
	//}

	// Slope gradient [m m-1]
	slope_grade_sin = sin(slope_rad);
	slope_grade_cos = cos(slope_rad);

	// Beta_0 as a function of maximum topographical distribution (Ekcic et al. 2019; eqn 3)
	double beta_0 = 1.357209;

	// Microtopography distribution (mm) (Ekcic et al. 2019; eqn 2)
	slope_sigma = pow(slope_rad + beta_0, -3.0) * MM_PER_M;

	// Constants
	sigma_sqrt2 = 1.0 / (slope_sigma * sqrt(2.0));
	sigma_sqrt2pi = slope_sigma / sqrt(2.0 * PI);
	sigma2 = 1.0 / (2.0 * pow(slope_sigma, 2));

	// Fraction of the area that is inundated threshold below which no 
	// single connected inundated area spans the grid cell
	wmass_ssw_frac_thresh = 0.4;

	// Calculate surface water threshold below which no single connected  
	// inundated area spans the grid cell (Ekcic et al. 2019 and CLM 5.0 code)
	wmass_ssw_thresh = 0.0;

	if (slope_sigma > 1.0e-6) {

		// Surface water level (mm)
		double d = 0.0;

		// Determine surface water threshold
		for (int i = 0; i < 4; i++) {
			double fd = 0.5 * (1.0 + erf(d * sigma_sqrt2)) - wmass_ssw_frac_thresh;
			double dfdd = exp(-pow(d, 2) * sigma2) / (slope_sigma * sqrt(2.0 * PI));
			d -= fd / dfdd;
		}

		// Surface water threshold below which no single connected  
		// inundated area spans the grid cell (mm)
		wmass_ssw_thresh = 0.5 * d * (1.0 + erf(d / (slope_sigma * sqrt(2.0)))) +
			slope_sigma / sqrt(2.0 * PI) * exp(-pow(d, 2) / (2.0 * pow(slope_sigma, 2)));
	}
}

/// Initiate soil water content
/** Initiate soil water content
 */
void Soilcolumn::init_water() {

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
			vws_liq[lyr] = soillayer.wsats / soillayer.depth;
			// time per mm
			dt_dz[lyr] = dtsub / soillayer.depth;
			Kw[lyr] = soillayer.Kw;
			Psi[lyr] = soillayer.Psi;
			dKwdw[lyr] = soillayer.dKwdw;
			dPsidw[lyr] = soillayer.dPsidw;
			dz[lyr] = soillayer.depth;
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
			soillayers[lyr].wmass += dwat[lyr] * soillayers[lyr].depth;
		}
	}

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wmass = min(soillayers[lyr].wmass, soillayers[lyr].wsats);
	}
}

/// Returns surface litter C content
/** Returns surface litter C content
 */
double Soilcolumn::get_surflittercmass() {
	Soillayer& soillayer = soillayers[0];
	return soillayer.sompool[SURFMETA].cmass + soillayer.sompool[SURFSTRUCT].cmass +
		soillayer.sompool[SURFFWD].cmass + soillayer.sompool[SURFCWD].cmass;
}

/// Returns surface litter N content
/** Returns surface litter N content
 */
double Soilcolumn::get_surflitternmass() {
	Soillayer& soillayer = soillayers[0];
	return soillayer.sompool[SURFMETA].nmass + soillayer.sompool[SURFSTRUCT].nmass +
		soillayer.sompool[SURFFWD].nmass + soillayer.sompool[SURFCWD].nmass;
}

/// Returns soil C content
/** Returns soil C content
 */
double Soilcolumn::ccont(grouppooltype group) {

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
double Soilcolumn::ncont(grouppooltype group) {

	double ncont = 0.0;
	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ncont += soillayers[lyr].ncont(group);
	}
	return ncont;
}

/// Returns surface standing water content
/** Returns surface standing water content
 */
double Soilcolumn::wcont_ssw() {
	return wmass_ssw + imass_ssw;
}

/// Returns a vector with soillayer values for a certain property
/** Returns a vector with soillayer values for a certain property
 */
void Soilcolumn::get(soilproptype SPT, double prop[NSOILLAYERS_MAX], int m) {

	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		switch (SPT) {
		case DEPTH:
			prop[lyr] = soillayers[lyr].depth;
			break;
		case TEMP:
			prop[lyr] = soillayers[lyr].temp;
			break;
		case MTEMP:
			prop[lyr] = soillayers[lyr].mtemp[m];
			break;
		case WMASS:
			prop[lyr] = soillayers[lyr].wmass;
			break;
		case MWMASS:
			prop[lyr] = soillayers[lyr].mwmass[m];
			break;
		case IMASS:
			prop[lyr] = soillayers[lyr].imass;
			break;
		case MIMASS:
			prop[lyr] = soillayers[lyr].mimass[m];
			break;
		case WCONT:
			prop[lyr] = soillayers[lyr].wcont();
			break;
		case MWCONT:
			prop[lyr] = soillayers[lyr].mwcont[m];
			break;
		case AWCONT:
			prop[lyr] = soillayers[lyr].awcont;
			break;
		case WSATS:
			prop[lyr] = soillayers[lyr].wsats;
			break;
		case AWC:
			prop[lyr] = soillayers[lyr].awc;
			break;
		case WP:
			prop[lyr] = soillayers[lyr].wp;
			break;
		case MDRAIN:
			prop[lyr] = soillayers[lyr].mdrain[m];
			break;
		case NH_4:
			prop[lyr] = soillayers[lyr].NH4_mass;
			break;
		case NO_3:
			prop[lyr] = soillayers[lyr].NO3_mass;
			break;
		case CH_4:
			prop[lyr] = soillayers[lyr].CH4;
			break;
		case CO_2:
			prop[lyr] = soillayers[lyr].CO2;
			break;
		case RESPSUM:
			prop[lyr] = soillayers[lyr].respsum;
			break;
		case WCONTENT:
			prop[lyr] = soillayers[lyr].wmass + soillayers[lyr].imass;
			break;
		case RUNOFFDRAIN:
			prop[lyr] = soillayers[lyr].runoff_drain;
			break;
		case RUNOFFBASE:
			prop[lyr] = soillayers[lyr].runoff_baseflow;
			break;
		case CSLOW:
			prop[lyr] = soillayers[lyr].cpool_slow;
			break;
		case CFAST:
			prop[lyr] = soillayers[lyr].cpool_fast;
			break;
		case ARH:
			prop[lyr] = soillayers[lyr].arespsum;
			break;
		case ANMIN:
			prop[lyr] = soillayers[lyr].anmin;
			break;
		case ANIMM:
			prop[lyr] = soillayers[lyr].animmob;
			break;
		case AORGCLEACH:
			prop[lyr] = soillayers[lyr].aorgCleach;
			break;
		case AORGNLEACH:
			prop[lyr] = soillayers[lyr].aorgNleach;
			break;
		case AMINLEACH:
			prop[lyr] = soillayers[lyr].aminleach;
			break;
		case AFLITTERINPUT:
			prop[lyr] = soillayers[lyr].acmassfinput;
			break;
		case ACLITTERINPUT:
			prop[lyr] = soillayers[lyr].acmasscinput;
			break;
		case VOL_ORG:
			prop[lyr] = soillayers[lyr].organic_volume;
			break;
		case VOL_MIN:
			prop[lyr] = soillayers[lyr].mineral_volume;
			break;
		case VOL_ROOT:
			prop[lyr] = soillayers[lyr].root_volume;
			break;
		case MODTEMP:
			prop[lyr] = soillayers[lyr].atemp_mod;
			break;
		case MODMOIS:
			prop[lyr] = soillayers[lyr].amoist_mod;
			break;
		case MODTOT:
			prop[lyr] = soillayers[lyr].atot_mod;
			break;
		default:
			// Should never happen
			fail("Unsupported soil property type (%d)\n", SPT);
		}
	}
}

/// Returns a vector with soillayer values for a certain property
/** Returns a vector with soillayer values for a certain property
 */
double Soilcolumn::get_sum(soilproptype SPT, int mth) {

	// Variable to return
	double prop = 0.0;

	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		prop += soillayers[lyr].get(SPT, mth);
	}

	// Return vector
	return prop;
}

/// Returns a vector with soillayer N fractions
/** Returns a vector with soillayer N fractions
 */
void Soilcolumn::get_nfrac(n_pref_type pref, double nfrac[NSOILLAYERS_MAX], double nmin, bool wmass_limited) {

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
void Soilcolumn::calc_rootfrac(double rootfrac[NSOILLAYERS_MAX]) {
	
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
void Soilcolumn::calc_rootmass(double cmass_root[NSOILLAYERS_MAX]) {

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

/// Returns available N for soil column
/** Returns available N for soil column
 */
double Soilcolumn::nmass_avail(n_pref_type pref, bool wmass_limited) {

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
void Soilcolumn::nmass_avail_uptake(double uptake, n_pref_type pref) {
	double tot_nmass_avail = nmass_avail(pref);

	if (!negligible(tot_nmass_avail)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			Soillayer& soillayer = soillayers[lyr];
			soillayer.nmass_subtract(uptake * soillayer.nmass_avail(pref) / tot_nmass_avail, pref);
		}
	}
}

/// Calculate water content for this month
/** Calculate water content for this month
 */
void Soilcolumn::calc_mwcont() {
	double acc_wcont = 0.0;
	double acc_awc = 0.0;

	// Monthly standing water and ice
	mwmass_ssw[date.month] += (wmass_ssw + imass_ssw) / date.ndaymonth[date.month];

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
void Soilcolumn::calc_mtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		soillayer.mtemp[date.month] += soillayer.temp / date.ndaymonth[date.month];
	}

	for (int plyr = 0; plyr < NPADLAYERS; plyr++) {
		padding& pad_lyr = pad[plyr];
		pad_lyr.mtemp[date.month] += pad_lyr.temp / date.ndaymonth[date.month];
	}
}

/// Average temperature down to a depth
/** Average temperature down to a depth
 */
double Soilcolumn::gtemp_rootdist(double rootdist[NSOILLAYERS_MAX]) {

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
void Soilcolumn::calc_temp_analytic() {

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

					k = (k * (soillayer.depth_profile - (soillayer_above.depth + soillayer.depth) / 2.0) +
						above_k * soillayer_above.depth / 2.0 + current_k * soillayer.depth / 2.0) / soillayer.depth_profile;
				}

				// Calculate parameters alag and exp(-alag) from Eqn 2
				soillayer.alag = soillayer.depth_profile * M_PER_MM / sqrt(k / HALF_OMEGA); // from Eqn 1
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
double Soilcolumn::var_depth(double depth_start, double depth_end, soilproptype SPT_var, soilproptype SPT_scal, int mth) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil column depth
	double depth_soil_end = min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.depth, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			scalar = soillayer.get(SPT_scal) * layer_depth / soillayer.depth;
			acc_var += soillayer.get(SPT_var, mth) * scalar;
			acc_scalar += scalar;
		}

		acc_depth += soillayer.depth;
		lyr++;
	}

	// If monthly temperature then padding layer can also be used
	if (depth_end > depth && SPT_var == MTEMP && SPT_scal == DEPTH) {
		int lyr = 0;

		while (acc_depth < depth_end && lyr < NPADLAYERS) {

			padding& padlyr = pad[lyr];

			lyr_start = min(max(acc_depth, depth_start), depth_end);
			lyr_end = min(max(acc_depth + padlyr.depth, depth_start), depth_end);
			layer_depth = lyr_end - lyr_start;

			if (layer_depth > 0.0) {
				acc_var += padlyr.mtemp[mth] * layer_depth;
				acc_scalar += layer_depth;
			}

			acc_depth += padlyr.depth;
			lyr++;
		}
	}

	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

/// Returns the value of some vector for a certain depth interval
/** Returns the value of some vector for a certain depth interval
 */
double Soilcolumn::var_depth(double depth_start, double depth_end, soilproptype SPT_var, int mth) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil column depth
	double depth_soil_end = min(depth, depth_end);

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.depth, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			acc_var += soillayer.get(SPT_var, mth) * layer_depth / soillayer.depth;
		}

		acc_depth += soillayer.depth;
		lyr++;
	}

	return acc_var;
}

/// Returns the value of soil C for a certain depth interval
/** Returns the value of soil C for a certain depth interval
 */
double Soilcolumn::var_depth_ccont(double depth_start, double depth_end, grouppooltype group) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil column depth
	//double depth_soil_end = min(depth, depth_end);
	double depth_soil_end = depth_end;

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.depth, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			scalar = layer_depth / soillayer.depth;
			acc_var += soillayer.ccont(group) * scalar;
		}

		acc_depth += soillayer.depth;
		lyr++;
	}

	return acc_var;
}

/// Returns the value of soil N for a certain depth interval
/** Returns the value of soil N for a certain depth interval
 */
double Soilcolumn::var_depth_ncont(double depth_start, double depth_end, grouppooltype group) {

	// Calculate variable values for layers within depth range
	double lyr_start, lyr_end;
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;

	// depth_end can't be larger than soil column depth
	//double depth_soil_end = min(depth, depth_end);
	double depth_soil_end = depth_end;

	while (acc_depth < depth_soil_end && lyr < nsoillayers) {

		Soillayer& soillayer = soillayers[lyr];

		lyr_start = min(max(acc_depth, depth_start), depth_soil_end);
		lyr_end = min(max(acc_depth + soillayer.depth, depth_start), depth_soil_end);
		layer_depth = lyr_end - lyr_start;

		if (layer_depth > 0.0) {
			scalar = layer_depth / soillayer.depth;
			acc_var += soillayer.ncont(group) * scalar;
		}

		acc_depth += soillayer.depth;
		lyr++;
	}

	return acc_var;
}

/// Snowfall
/** Snowfall
*  Calculate snowfall density and redistribute snow layers
*/
void Soilcolumn::snowfall(double airT, double& prec, double u10) {

	if (negligible(wmass_ssw + imass_ssw, -10) || snow.nsnowlayers) {

		// Fresh snow density calculation, following eqn 1 in Vionnet et al (2012)
		// Constant [kg m-3]
		int a = 109;
		// Constant [kg m-3 K-1]
		int b = 6;
		// Constant [kg m-7/2 s-1/2]
		int c = 26;
		// Minimum snow density 
		double dens_min = 80;
		// Temperature of the melting point for water
		double Tmelt = 0;

		// Snow fall density
		double snowfall_dens = max(dens_min, a + b * (airT - Tmelt) + c * pow(u10, 0.5));

		// Depth of the newly fallen snow based on precipitation & density (Oleason 2013, Eq. 6.88)
		double snowfall_depth = prec / (snowfall_dens / rho_H2O);

		// Snow fall heat capacity
		double snowfall_Ch = snowheatcap(prec, 0.0, snowfall_depth);

		// Snowpack depth
		double snowpack_depth = snow.depth();

		// Add snowfall and redistribute snowpack
		snow.redist_layers(snowpack_depth, snowfall_depth, prec, airT, snowfall_Ch);
	}
	else {
		// Update surface standing water
		update_ssw_temp(airT, 0.0, prec);

		// Update standing ice content (mm)
		update_ssw(imass_ssw, prec);

		// Do phase change of standing suface water
		phase_change_ssw();

		// All standing surface water has turned into ice. Create a snow layer
		// consisting of only ice.
		if (negligible(wmass_ssw)) {
			snow.snowlayers.resize(snow.nsnowlayers + 1);
			snow.snowlayers[snow.nsnowlayers] = Snowlayer(imass_ssw, 0.0, rho_ICE, temp_ssw);
			snow.snowlayers[snow.nsnowlayers].depth = snow.snowlayers[snow.nsnowlayers].imass * rho_H2O / snow.snowlayers[snow.nsnowlayers].dens;
			snow.nsnowlayers++;
			imass_ssw = 0.0;
			temp_ssw = 0.0;
		}
	}

	// Reset precipitation
	prec = 0.0;
}

/// Snow accumulation
/** Snow accumulation
 * Rain or snowfall, and melt if simple snow scheme is used
 */
void Soilcolumn::snow_or_rainfall(double airT, double& prec, double u10, double& rain_melt, double& rain_melt_temp) {

	// No precipitation and freezing temperatures, then return from function
	if (prec == 0 && airT < 0.0)
		return;

	double w_begin = wcont_ssw() + snow.wcontent() + prec + rain_melt;

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
			snowlayer.depth = snowlayer.imass / (snowlayer.dens / rho_H2O);

			if (negligible(snowlayer.imass)) {

				// Erase this empty snow layer
				snow.snowlayers.erase(snow.snowlayers.begin());
				snow.nsnowlayers--;
			}
		}
	}

	// Water balance might not meet mass balance in very low productive and saturated sites 
	double wdiff = w_begin - (wcont_ssw() + snow.wcontent() + rain_melt);
	if (wdiff > 0.0) {
		rain_melt += wdiff;
	}
	else if (wdiff < 0.0) {
		rain_melt += wdiff;

		// Still missing water
		if (rain_melt < 0.0) {

			// Extract first from standing water
			if (wcont_ssw() > 0.0) {
				double frac = min(1.0, -rain_melt / wcont_ssw());
				rain_melt += frac * wcont_ssw();
				wmass_ssw *= 1.0 - frac;
				imass_ssw *= 1.0 - frac;
			}

			// Still missing water. Take from snow
			if (rain_melt < 0.0) {
				if (snow.wcontent() > 0.0) {
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
								snow.nsnowlayers = (int)snow.snowlayers.size();
							}
							else {
								snowlayer.dens = snowlayer.imass / snowlayer.depth * rho_H2O;
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

	double w_end = wcont_ssw() + snow.wcontent() + rain_melt;
	if (fabs(w_begin - w_end) > 1.0e-10) {
		dprintf("Year %d Day %d ERROR: Snow snow_or_rainfall diff %g\n", date.year, date.day, w_begin - w_end);
	}
}

/// Routine that calculate surface evaporation, infiltration and runoff
/** Routine that calculate surface evaporation, infiltration and runoff
 */
void Soilcolumn::compute_evap_infilt_and_runoff(double dtsub, double wflux_rain_melt, double depth_infilt, 
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
	double wmass_infilt = var_depth(0.0, depth_infilt, WMASS);

	// Fraction pore spaces that are filled with water and ice
	double totwcont_infilt = min(1.0, (wmass_infilt + imass_infilt) / wsats_infilt);

	// Fractional area with water table at surface (Wood et al. 1992)
	fsat = 1.0 - pow(1.0 - totwcont_infilt, b_infilt / (1.0 + b_infilt));

	// Calculate surface water fraction and runoff
	if (wmass_ssw > wmass_ssw_min) {

		// Surface water level(m)
		double d = 0.0;

		double dfdd, fd = 1.0;
		int i = 0;
		while (fd > 1.0e-8 && i < 10) {
			fd = 0.5 * d * (1.0 + erf(d * sigma_sqrt2)) + sigma_sqrt2pi * exp(-pow(d, 2) * sigma2) - wmass_ssw;
			dfdd = 0.5 * (1.0 + erf(d * sigma_sqrt2));
			d -= fd / dfdd;
			i++;
		}

		// Update the submerged area fraction using the new d value (Ekici et al. 2019; eqn 1)
		wmass_ssw_frac = 0.5 * (1.0 + erf(d * sigma_sqrt2));

		// Update the interconnected fraction of surface water
		wmass_ssw_con_frac = 0.0;
		if (wmass_ssw_frac > wmass_ssw_frac_thresh) {
			wmass_ssw_con_frac = pow(wmass_ssw_frac - wmass_ssw_frac_thresh, 0.14);
		}

		// Calculate surface water lateral runoff
		wflux_ssw_runoff = 0.0;
		if (wmass_ssw > wmass_ssw_thresh) {
			wflux_ssw_runoff = min(1.0 / dtsub, C_dw * max(0.001, slope_grade_sin) * wmass_ssw_con_frac) * (wmass_ssw - wmass_ssw_thresh);
		}
	}
	else {
		// Very small amount of surface water -> all lost through surface runoff
		wmass_ssw_con_frac = 0.0;
		wmass_ssw_frac = 0.0;
		wflux_ssw_runoff = wmass_ssw / dtsub;
	}

	// Fast runoff directly to streams
	if (wflux_rain_melt > 0.0) {

		double fsat_surf = min(fsat_max, fsat);

		// Fractional area with water table at surface is assumed to be under standing water
		if (fsat_surf > wmass_ssw_frac) {

			// Surface runoff (mm s-1)
			wflux_surf_runoff = (fsat_surf - wmass_ssw_frac) * wflux_rain_melt;

			// Water staying at surface / soil (mm s-1)
			wflux_in = wflux_rain_melt - wflux_surf_runoff;
		}
		else {
			// No water to runoff, all water staying at surface / soil (mm s-1)
			wflux_in = wflux_rain_melt;
		}
	}

	// Partition surface inputs between soil and standing water
	wflux_ssw = wmass_ssw_frac * wflux_in;
	wflux_soil = wflux_in - wflux_ssw;

	// Remove surface water runoff
	wflux_ssw -= wflux_ssw_runoff;

	// Evaporation

	// Evaporation is calculated on a sub daily time step as water content 
	// in top soil layer changes with each time step

	// Evaporation for high-latitude peatland soils
	if (patch.stand.is_highlatitude_peatland_stand()) {

		// Evaporation from peat soil surface (Wania et al 2009)
		if (snow.depth() < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
			wflux_soil_evap = fevap * get_climate().eet * PRIESTLEY_TAYLOR * (0.99 / (1.0 + exp(-1.0 * (-wtd + 98.7) / 22.6) + 0.02));
			wflux_soil_evap *= (1.0 - wmass_ssw_frac) / SECS_PER_DAY;
		}
		else {
			wflux_soil_evap = 0.0;
		}
	}
	else {
		// Evaporation for mineral soils

		wcont_evap = soillayers[0].wcont();

		// Calculate total evaporation only if snow depth < 10mm (mm day-1)
		wflux_soil_evap = 0.0;
		if (snow.depth() < 10.0 && !negligible(fevap)) {
			wflux_soil_evap = get_climate().eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
			
			// Remove small evaporations
			if (wflux_soil_evap < EPS) {
				wflux_soil_evap = 0.0;
			}
			wflux_soil_evap *= (1.0 - wmass_ssw_frac) / SECS_PER_DAY;
		}
	}

	// Evaporation for standing water part of the soil
	double wcont_ssw_evap = soillayers[0].wcont(wmass_ssw);

	// Calculate total evaporation only if snow depth < 10mm (mm day-1)
	wflux_ssw_evap = 0.0;
	if (snow.depth() < 10.0 && !negligible(fevap)) {
		wflux_ssw_evap = get_climate().eet * PRIESTLEY_TAYLOR * wcont_ssw_evap * wcont_ssw_evap * fevap;

		if (wflux_ssw_evap < EPS) {
			wflux_ssw_evap = 0.0;
		}
		wflux_ssw_evap *= wmass_ssw_frac / SECS_PER_DAY;
	}

	// Remove evaporation
	wflux_ssw -= wflux_ssw_evap;
	wflux_soil -= wflux_soil_evap;

	// Determine the maximum infiltration rate (mm s-1)
	wflux_infilt_max = (1.0 - fsat) * imp_ice * Kw_max * slope_grade_cos;

	// Excess water that can't infiltrate (mm s-1)
	wflux_soil_excess = max(0.0, wflux_soil - (1.0 - wmass_ssw_frac) * wflux_infilt_max);

	// Soil infiltration and surface water run-on
	wflux_ssw += wflux_soil_excess;
	wflux_soil -= wflux_soil_excess;

	// Temporary storage prior to calculating bottom drainage from surface water
	double wmass_ssw_tmp = wmass_ssw + wflux_ssw * dtsub;

	// Negative surface water -> reduce soil infiltration
	if (wmass_ssw_tmp < 0.0) {
		wflux_soil += wmass_ssw_tmp / dtsub;
		wflux_ssw = -wmass_ssw / dtsub;
		wflux_ssw_drain = 0.0;
	}
	else {
		// Drainage from surface water to soil column
		wflux_ssw_drain = min(wmass_ssw_frac * wflux_infilt_max, wmass_ssw_tmp / dtsub);
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

/// Update soil matrix potential for the soil column
/** Update soil matrix potential for the soil column
 */
void Soilcolumn::update_matric_potential() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_matric_potential();
	}
}

 /// Updates hydraulic conductivity for the soil column
 /** Updates hydraulic conductivity for the soil column
  * Relate the hydraulic conductivity to the soil water 
  * content following Campbell (1974)
  */
void Soilcolumn::update_hydraulic_properties(int nsoillayers_sim) {

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

/// Update soil thermal capacity and conductivity for the soil column
/** Update soil thermal capacity and conductivity for the soil column
 */
void Soilcolumn::update_thermal_properties() {

	// Update heat capacity and conductivity for each snow layer
	snow.update_thermal_properties();

	// Calculate heat capacity and conductivity for surface standing water
	update_thermal_properties_for_ssw();

	// Update heat capacity and conductivity for each soil layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		soillayers[lyr].update_thermal_properties();
	}

	// Calculate heat capacity and conductivity for each padding layer
	update_thermal_properties_for_pad();
}


/// Update surface standing waters thermal capacity and conductivity
/** Update surface standing waters thermal capacity and conductivity
 */
void Soilcolumn::update_thermal_properties_for_ssw() {

	double depth_ssw = wmass_ssw + imass_ssw;

	// No standing water then return
	if (negligible(depth_ssw, -10)) {
		return;
	}

	// Standing water heat capacity (J m-3 K-1)
	Ch_ssw = (wmass_ssw * Cp_water + imass_ssw * Cp_ice) / depth_ssw;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	// Calculate log terms
	double lKfwater = (wmass_ssw / depth_ssw) * Kwater_log;
	double lKfice = (imass_ssw / depth_ssw) * Kice_log;

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh_ssw = exp_thermal_lut(lKfwater + lKfice);
}

/// Update padding layers thermal capacity and conductivity
/** Update padding layers thermal capacity and conductivity
 */
void Soilcolumn::update_thermal_properties_for_pad() {

	// Bedrock parameters from Chadburn et al. (2015)

	// Thermal conductivity [W m-1 K-1 s-1]
	double Kbedrock = 8.6;
	// Heat capacity [J m-3 K-1]
	double Cbedrock = 2100000.0;	

	// Use bottom soil layer properties for the first 5 padding layers
	const int first_bedrock_padding_layer = 5;

	for (int lyr = 0; lyr < NPADLAYERS; lyr++) {

		padding& pad_lyr = pad[lyr];

		// First padding layers use bottom soil layers properties
		if (lyr < first_bedrock_padding_layer) {
			Soillayer& soillayer = soillayers[nsoillayers - 1];
			pad_lyr.Kh = soillayer.Kh;
			pad_lyr.Ch = soillayer.Ch;
		}
		else {
			// Bedrock
			pad_lyr.Kh = Kbedrock;
			pad_lyr.Ch = Cbedrock;
		}
	}
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
		else {
			dplus = ((ki[lyr] * dz[lyr] + ki[lyr + 1] * dz[lyr + 1]) / (dz[lyr] + dz[lyr + 1])) / ci[lyr];
		}

		// D-
		if (lyr == 0) {
			dminus = ki[lyr] / ci[lyr];
		}
		else if (lyr == nlayers - 1) {
			dminus = Dbedrock;
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

/// Update surface standing water temperatures due to water addition
/** Update surface standing water temperatures due to water addition
 */
void Soilcolumn::update_ssw_temp(double temp_add, double wmass_add, double imass_add) {

	// Only update standing water temperature if the added water is positive. 
	// Negative would be removal of water and would not change temperature.
	if (wmass_add > 1.0e-6 || imass_add > 1.0e-6) {
		temp_ssw = (temp_ssw * (wmass_ssw + imass_ssw) * Ch_ssw + temp_add * wmass_add * Cp_water +
			temp_add * imass_add * Cp_ice) /
			((wmass_ssw + imass_ssw) * Ch_ssw + wmass_add * Cp_water + imass_add * Cp_ice);
	}
}

/// Update surface standing water content
/** Update surface standing water content
 */
void Soilcolumn::update_ssw(double& mass, double mass_add) {

	mass += mass_add;

	// Very small amount of surface standing water
	// then it is lost as suface runoff
	if (negligible(mass, -10)) {
		patch.fluxes.report_flux(Fluxes::WRUNOFF_SSW, mass);
		mass = 0.0;
	}
}


/// Update soil layer temperatures due to water movement
/** Update soil layer temperatures due to water movement
 */
void Soilcolumn::wflow_temp(double temp_infilt) {

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

/// Serialize soil padding variables
/** Serialize soil padding variables
 */
void padding::serialize(ArchiveStream& arch) {
	arch & temp
		& depth;
}

/// Set soil padding depths
/** Calculate padding layer depths for Crank-Nicholson bottom boundary
 * condition setup. The padding layer thicknesses are calculated as
 * a geometric sequence with a given total padding thickness and
 * number of layers. This gives a smooth increase in layer depth
 * away from the calculation domain.
 */
void Soilcolumn::set_padding_depth() {

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
	double dz = soillayers[nsoillayers - 1].depth;

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
			pad[lyr].depth = lyr_depth;
			tot_depth += lyr_depth;
		}

		pad[NPADLAYERS - 1].depth += PAD_DEPTH - tot_depth;
	}
	// Else each padding layer decrease linear to a minimum of 2m for the bottom padding layer
	else {
		double first_pad = min(16000.0, soillayers[nsoillayers - 1].depth);
		double diff = floor((5.0 * first_pad - PAD_DEPTH) / 10.0 / 100.0) * 100.0;
		int scalar = 0;
		double tot_depth = 0.0;

		for (int lyr = 0; lyr < NPADLAYERS - 1; lyr++) {
			pad[lyr].depth = first_pad - diff * (double)scalar;
			tot_depth += pad[lyr].depth;
			scalar++;
		}

		pad[NPADLAYERS - 1].depth = PAD_DEPTH - tot_depth;
	}
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
void Soilcolumn::update_acrotelm_co2(double atmo_co2) {
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
void Soilcolumn::hydrology_peat(const Climate& climate, double fevap) {

	// Total AET (mm)
	double aet_total = 0.0;
	// Transpiration from the acrotelm (mm)
	double aet_acrotelm = 0.0;
	// Runoff from upper soil layer (mm)
	double runoff_surf;
	// Runoff from standing water (mm)
	double runoff_ssw;
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
	double wtot_acrotelm = var_depth(0.0, ACROTELM_DEPTH, WMASS);
	// Total ice in acrotelm layers (mm)
	double itot_acrotelm = var_depth(0.0, ACROTELM_DEPTH, IMASS);
	// Initial water mass in the acrotelm
	double wtot_acrotelm_init = rain_melt;
	// Minimum water mass in the acrotelm
	double wtot_acrotelm_min = 0.0;

	// Calculate minimum water in acrotelm
	for (int lyr = 0; lyr < nacrotelm; lyr++) {
		Soillayer& soillayer = soillayers[lyr];
		if (lyr == 0) {
			wtot_acrotelm_min += minvtot * soillayer.depth;
		}
		else {
			wtot_acrotelm_min += soillayer.wp;
		}
		wtot_acrotelm_init += soillayer.wmass;
	}
	
	// Transpiration

	// Retrieve Vegetation for this patch
	Vegetation& vegetation = patch.vegetation;

	// Sum AET for across all vegetation individuals in this patch

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		int lyr = 0;
		while (soillayers[lyr].acrotelm) {
			pft_acro_root_frac += patch.pft[indiv.pft.id].pft.rootdist[lyr];
			lyr++;
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

		// Add water that has been taken about available
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
		wtot_acrotelm += runon - runoff_surf + wmass_ssw;

		// If negative total water, then all initial water has gone to runoff and water is set to zero
		if (wtot_acrotelm < 0.0) {
			wtot_acrotelm = 0.0;
			runoff_surf = wtot_acrotelm_init;
			runon = 0.0;
			wmass_ssw = 0.0;
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
		double wsats_acrotelm = var_depth(0.0, ACROTELM_DEPTH, WSATS);

		// Calculate ice impedance factor
		double imp_ice = ice_imp(itot_acrotelm / wsats_acrotelm);

		// Average hydraulic conductivity at saturation for infiltation depth (top layer)
		double Kw_max = var_depth(0.0, ACROTELM_DEPTH, KW_MAX, DEPTH);

		// Calculate infiltration and surface runoff
		compute_evap_infilt_and_runoff(dtsub, wflux_rain_melt, ACROTELM_DEPTH, itot_acrotelm,
			wsats_acrotelm, imp_ice, Kw_max, fevap, wflux_evap, wflux_runoff_surf, wflux_runoff_ssw,
			wflux_surf, wflux_infilt, wflux_ssw_drain);

		// Water added to standing water
		double wmass_ssw_added = wflux_surf * dtsub;

		// Update standing water temperature
		update_ssw_temp(rain_melt_temp, wmass_ssw_added);

		// Update standing water content (mm)
		update_ssw(wmass_ssw, wmass_ssw_added);

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

		// Saturated acrotelm

		wtp = wtot_acrotelm + itot_acrotelm - ACROTELM_DEPTH * acrotelm_por;
		
		// WTP > 0 above the surface
		if (wtp > maxh) {
			// Update runoff from the surface
			if (wania_approach) {
				runoff_surf += wtp - maxh;
			}
			else {
				double wmass_added = wtp - maxh;

				// Update standing water temperature // DW_COMMENT don't know where the excess water comes from???
				//update_ssw_temp(soillayers[0].temp, wmass_added);

				// Update standing water content (mm)
				update_ssw(wmass_ssw, wmass_added);
			}
			wtp = maxh;
		}

		// WTD defined as < 0 above the surface
		wtd = -wtp;

		wtot_acrotelm = ACROTELM_DEPTH * acrotelm_por - itot_acrotelm;
	}
	else {
		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * az));

		if (wtd > zmin)
			wtd = 3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * (acrotelm_por - minvtot));
		
		wtd = min(ACROTELM_DEPTH, wtd);

		wtp = -wtd;
	}

	// Water content in each acrotelm soil layer depending on wtd

	// Save acrotelm water mass before update
	wtot_acrotelm_init = wtot_acrotelm;

	// Depth of acrotelm sublayers - 10 mm
	double Dz_sub = 10.0;

	// Number of sublayers in acrotelm
	int nsublayers_acro = (int)(min(ACROTELM_DEPTH, depth) / Dz_sub);

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
			int sublayers_per_layer = (int)(soillayer.depth / Dz_sub);

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
					water_avail += max(0.0, soillayer.wmass - minvtot * soillayer.depth);
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
					wsub = frac * max(0.0, soillayer.wmass - minvtot * soillayer.depth);
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

bool Soilcolumn::is_wet_season() {

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
void Soilcolumn::hydrology_wetland(const Climate& climate, double fevap) {

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

	double wmass_begin = get_sum(WMASS) + rain_melt - wmass_wetland_added;

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
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
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
	double wmass_end = get_sum(WMASS) + runoff_surf + evap + aet_total - wmass_wetland_added;

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
void Soilcolumn::hydrology_gerten(const Climate& climate, double fevap) {

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
	const double K_AET_DEPTH = (max(1.0, soillayers[0].depth / SOILDEPTH_EVAP) - 1.0) *
		(K_AET / K_DEPTH - 1.0) / (1.0 / K_DEPTH - 1.0) + 1.0;

	// Water mass balance check
	double wbefore = get_sum(WMASS) + rain_melt;

	// Initial water masses of each layer
	double wmass_org[NSOILLAYERS_MAX] = { 0.0 };
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wmass_org[lyr] = soillayers[lyr].wmass;
	}

	// Actual Evapotranspiration from soil column

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
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
			aet_total += aet;
			aet_layer[lyr] += aet;
		}
		vegetation.nextobj();
	}

	// Evaporation from soil surface (Bondeau et al. (2007),  Rost et al. (2008))
	
	Soillayer& soillayer_top = soillayers[0];

	double evap = 0.0;

	// Evaporation only if snow equivalent water < 10mm
	if (snow.wcontent() < 10.0) {	
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
	wcont_evap += (rain_melt - aet_layer[0] * K_AET_DEPTH * min(1.0, SOILDEPTH_EVAP / soillayer_top.depth) - evap)
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
		perc = min(min(1.0, SOILDEPTH_EVAP / soillayer_top.depth) * soillayer_top.perc_base * pow(wcont_evap, soiltype.perc_exp),
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
			perc = min(soillayer_above.perc_base * pow(soillayer_above.wcont(), soiltype.perc_exp), max_rain_melt);
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
			pow(soillayer_bottom.wcont(), soiltype.perc_exp);
		
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
	double wafter = get_sum(WMASS) + runoff_surf + get_sum(RUNOFFDRAIN) + get_sum(RUNOFFBASE) + aet_total + evap;
	if (fabs(wbefore-wafter) > 1.0e-10)
		dprintf("Year %d Day %d BUCKETWATER Before %g After %g Diff %g Rain_melt %g\n",date.year,date.day, wbefore, wafter, wbefore - wafter, rain_melt);

	// Reset rain_melt
	rain_melt = 0.0;

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, 0.0);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(RUNOFFDRAIN));
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(RUNOFFBASE));
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
void Soilcolumn::hydrology_miller(const Climate& climate, double fevap) {

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
	// wcont_evap = water mass of evaporation sublayers (20 cm) at top of soil column
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
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
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
			perc = min(soiltype.perc_base * pow(wcont_layer1, soiltype.perc_exp), wmass_layer1);

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

			perc_from_base = min(wmass_perc_layer_2, min(BASEFLOW_FRAC * soiltype.perc_base * pow(wcont_perc_layer_2, soiltype.perc_exp), max_rain_melt));

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
	if (patch.stand.is_true_wetland_stand() && wmass_wetland_added > 0.0) {

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
void Soilcolumn::hydrology_richards(const Climate& climate, double fevap) {

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
	const double DEPTH_INFILT_MAX = 100.0;

	// Effective depth of surface runoff [mm]
	double depth_infilt = min(thaw, DEPTH_INFILT_MAX);

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
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
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
	const double WFPS_SAT = 0.98;
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
	if (wflux_rain_melt > 0.0 || wmass_ssw > 0.0) {
		dtsub = SECS_PER_DAY / nsubdays_prec;
	}

	// Finished day
	bool finished = false;

	// Calculate static ice mass and saturation level for infiltation depth (top layer)
	double imass_infilt = var_depth(0.0, depth_infilt, IMASS);
	double wsats_infilt = var_depth(0.0, depth_infilt, WSATS);

	// Calculate ice impedance factor for infiltation depth (top layer)
	double imp_ice = ice_imp(imass_infilt / wsats_infilt);

	// Hydraulic conductivity at saturation for infiltation depth (top layer)
	double Kw_max = var_depth(0.0, depth_infilt, KW_MAX, DEPTH);

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
		// layers. 
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
			compute_evap_infilt_and_runoff(dtsub, wflux_rain_melt, depth_infilt, imass_infilt, 
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
					dwat[lyr] = wmass_add / soillayers[lyr].depth;
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
					vws_liq[lyr] = soillayer.wsats / soillayer.depth;
					// time per mm
					dt_dz[lyr] = dtsub / soillayer.depth;
					Kw[lyr] = soillayer.Kw;
					Psi[lyr] = soillayer.Psi;
					dKwdw[lyr] = soillayer.dKwdw;
					dPsidw[lyr] = soillayer.dPsidw;
					dz[lyr] = soillayer.depth;
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
		update_ssw_temp(rain_melt_temp, wmass_ssw_added);

		// Update standing water content (mm)
		update_ssw(wmass_ssw, wmass_ssw_added);

		// Accumulate evaporation
		evap += wflux_evap * dtsub;

		// Accumulate surface runoff
		runoff_surf += wflux_runoff_surf * dtsub;
		runoff_ssw += wflux_runoff_ssw * dtsub;

		// Temperature infiltration
		if (!negligible(wmass_infilt_tot + (wflux_infilt + wflux_ssw_drain) * dtsub)) {
			temp_infilt_tot = (temp_infilt_tot * wmass_infilt_tot + rain_melt_temp * wflux_infilt * dtsub +
				temp_ssw * wflux_ssw_drain * dtsub) / (wmass_infilt_tot + (wflux_infilt + wflux_ssw_drain) * dtsub);
		}

		// Accumulater infiltration
		wmass_infilt_tot += (wflux_infilt + wflux_ssw_drain) * dtsub;

		// Update the mass of liquid water
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].wmass += dwat[lyr] * soillayers[lyr].depth;
		}

		// Extract root uptake from soil layer not accounted for
		if (nsoillayers_sim != 1 && nsoillayers_sim != nsoillayers) {
			for (int lyr = nsoillayers_sim; lyr < nsoillayers; lyr++) {
				soillayers[lyr].wmass -= aet_layer[lyr] * dtsub / SECS_PER_DAY;
			}
		}

		// drainage from the bottom of the soil column
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
					update_ssw(wmass_ssw, wmass_excess);
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
				double drained_wmass = min(wmass_avail, C_dw * soillayer.Kw * ice_imp(soillayer.imass / soillayer.wsats) * slope_grade_sin * soillayer.depth * M_PER_MM * dtsub);
				soillayer.runoff_drain += drained_wmass;
				soillayer.wmass -= drained_wmass;
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

	// Reset ran_melt as water has entered the soil
	rain_melt = 0.0;

	// Save fluxes
	Fluxes& fluxes = patch.fluxes;
	fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	fluxes.report_flux(Fluxes::WRUNOFF_SSW, runoff_ssw);
	fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(RUNOFFDRAIN));
	fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(RUNOFFBASE));
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

void Soilcolumn::compute_moisture_fluxes_and_derivs(double wflux_infl, double Kw[NSOILLAYERS_MAX], 
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

void Soilcolumn::compute_RHS_moisture_form(double wflux_root[NSOILLAYERS_MAX], double qin[NSOILLAYERS_MAX],
	double qout[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX], double rmx[NSOILLAYERS_MAX], int nsoillayers_sim) {


	// Calculate RHS of moisture - based form of Richards equation

	// wflux_root - vertically distributed transpiration sink (mm H2O / s) (+= to atm)

	// Compute RHS of moisture - based form of Richards equation
	for (int lyr = 0; lyr < nsoillayers_sim; lyr++) {

		double fluxNet = qin[lyr] - qout[lyr] - wflux_root[lyr];
		rmx[lyr] = -fluxNet * dt_dz[lyr];
	}
}

void Soilcolumn::compute_LHS_moisture_form(double dt_dz[NSOILLAYERS_MAX], double dqidw0[NSOILLAYERS_MAX], 
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

void Soilcolumn::TridiagonalCol(double a[NSOILLAYERS_MAX], double b[NSOILLAYERS_MAX], double c[NSOILLAYERS_MAX], 
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

void Soilcolumn::check_error(double dwat[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], double dt_dz[NSOILLAYERS_MAX],
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
					double frac= (frac_cur + frac_pos) / 2.0;

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

void Soilcolumn::update_wtp() {

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
				wtp += pow(frac, 4) * soillayer.depth;
			}			
		}
	}

	wtd = -wtp;
}


// Initial percolation
/** Initial percolation for two layer soil scheme
 * and for peat and wetlands
 */
void Soilcolumn::initial_percolation() {

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
void Soilcolumn::phase_change_ssw() {

	// Initial water and ice content
	double wtotal_initial = wmass_ssw + imass_ssw;

	// Surface standing water depth
	double depth_ssw = wmass_ssw + imass_ssw;

	if (negligible(depth_ssw, -10)) {
		patch.fluxes.report_flux(Fluxes::WRUNOFF_SSW, depth_ssw);
		wmass_ssw = 0.0;
		imass_ssw = 0.0;
		return;
	}

	// Update thermal properties for standing water
	update_thermal_properties_for_ssw();

	// Initial energy
	double energy_initial = Ch_ssw * (temp_ssw + K2degC);

	// Should ice thaw
	bool ifThaw = false;
	// Should water freeze
	bool ifFreeze = false;
	// Fraction of water which freezes
	double Ffreez = 0.0;
	// Fraction of ice which thaws
	double Fthaw = 0.0;
	// Fraction of water
	double Fwmass = wmass_ssw / depth_ssw;
	// Fraction of ice
	double Fimass = imass_ssw / depth_ssw;
	// Energy per volume [J m-3]
	double energy;
	// Temperature difference (current - previous)
	double delta_T;
	// Phase change energy consumption [J]
	double phase_change_energy = 0.0;
	// Warning error
	const double EPS_warning = 1.0e-10;
	// Fail error
	const double EPS_fail = 1.0e-6;

	// Should water freeze
	if (temp_ssw < 0.0 && Fwmass > 0.0) {
		ifFreeze = true;
	}
	// Should ice thaw
	else if (temp_ssw > 0.0 && Fimass > 0.0) {
		ifThaw = true;
	}

	double temp_old = 0.0;

	// Temperature difference from previous phase change
	delta_T = temp_ssw - temp_old;

	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy from full standing water to freeze water [J m-3]
		energy = Ch_ssw * delta_T;

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		// Can't freeze all
		if (Fwmass > Ffreez) {

			// Update water and ice content
			wmass_ssw = (Fwmass - Ffreez) * depth_ssw;
			imass_ssw = (Fimass + Ffreez) * depth_ssw;

			// Update standing water temperature
			temp_ssw = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Ffreez * Lheat;
		}
		else {
			// It is cold enough to freeze all available water for freezing
			// and decrease the temperature

			// Update water and ice content
			imass_ssw += wmass_ssw;
			wmass_ssw = 0.0;

			// Energy used for phase change and to keep standing water at initial temperature [J m-3]
			double energy_loss = Fwmass * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature
			double energy_remain = energy - energy_loss;

			// Update standing water temperature
			temp_ssw = temp_old + energy_remain / Cp_ice;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Fwmass * Lheat;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		// Energy from standing water to thaw water [J m-3]
		energy = Ch_ssw * delta_T;

		// Fraction ice that can be thawed [m3 m-3]
		Fthaw = -energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		// Can't thaw all
		if (Fimass > Fthaw) {

			// Update water and ice content
			wmass_ssw = (Fwmass + Fthaw) * depth_ssw;
			imass_ssw = (Fimass - Fthaw) * depth_ssw;

			// Standing water temperature stays at previous timesteps temperature
			temp_ssw = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fthaw * Lheat;
		}
		else {

			// It is warm enough to thaw all available ice for thawing
			// and increase the temperature

			// Update water and ice content
			wmass_ssw += imass_ssw;
			imass_ssw = 0.0;

			// Energy used for phase change and to keep standing water at initial temperature [J m-3]
			double energy_loss = -Fimass * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update standing water temperature
			temp_ssw = temp_old + energy_remain / Cp_water;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fimass * Lheat;
		}
	} // End of Thawing


	// Final water
	double wtotal_final = wmass_ssw + imass_ssw;

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
	update_thermal_properties_for_ssw();

	// Final energy
	double energy_final = Ch_ssw * (temp_ssw + K2degC) - phase_change_energy;

	// Energy balance checks
	if ((fabs(energy_initial - energy_final)) > energy_initial * EPS_warning) {
		dprintf("WARNING! Year %d Day %d Energy balance during freezing or thawing in Soilcolumn::phase_change_ssw. Init %g, Final %g, Diff %g \n",
			date.year, date.day, energy_initial, energy_final, energy_initial - energy_final);
	}

	// Ice despite positive temperature
	if (temp_ssw > 0.0 && imass_ssw > 0.0) {
		dprintf("ERROR: Phase Change Stand Year %d Day %d Temp %g imass %g\n",
			date.year, date.day, temp_ssw, imass_ssw);
	}

	// No negative wmass
	if (wmass_ssw < -EPS_warning) {
		fail("ERROR: Negative wmass at end of phase_change_ssw! Year %d Day %d depth %g delta_T %g Phase_energy %g wmass %g imass %g\n",
			date.year, date.day, depth_ssw, delta_T, phase_change_energy, wmass_ssw, imass_ssw);
	}

	// No negative imass
	if (imass_ssw < -EPS_warning) {
		fail("ERROR: Negative imass at end of phase_change_ssw! Year %d Day %d depth %g delta_T %g Phase_energy %g wmass %g imass %g\n",
			date.year, date.day, depth_ssw, delta_T, phase_change_energy, wmass_ssw, imass_ssw);
	}

	// Very small amount of surface standing water
	// then it is lost as surface runoff
	if (negligible(wmass_ssw, -10)) {
		patch.fluxes.report_flux(Fluxes::WRUNOFF_SSW, wmass_ssw);
		wmass_ssw = 0.0;
	}

	// Very small amount of surface standing ice
	// then it is lost as surface runoff
	if (negligible(imass_ssw, -10)) {
		patch.fluxes.report_flux(Fluxes::WRUNOFF_SSW, imass_ssw);
		imass_ssw = 0.0;
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
void Soilcolumn::phase_change() {

	// Snow phase change

	// Phase change between water and ice and adjust snow layer temperatures accordingly
	snow.phase_change();

	// Standing water phase change
	phase_change_ssw();

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
void Soilcolumn::thaw_depth() {

	// Daily thawing depth, where at least some water is found, but ice can be present [mm]
	thaw = 0.0;

	// Add soil layer where all ice is thawn
	int lyr = 0;
	while (lyr < nsoillayers && soillayers[lyr].temp > 0.0) {
		thaw += soillayers[lyr].depth;
		lyr++;
	}
	// If all soil layer where thawn then add padding layer with temperature above 0 degree C
	if (lyr == nsoillayers) {
		while (lyr < nsoillayers + NPADLAYERS && pad[lyr - nsoillayers].temp > 0.0) {
			thaw += pad[lyr - nsoillayers].depth;
			lyr++;
		}
	}

	// If still among soil layers then assuming that water in next soil layer is at the top of the layer and ice below
	if (lyr < nsoillayers) {
		Soillayer& soillayer = soillayers[lyr];
		if (!negligible(soillayer.wmass + soillayer.imass)) {
			thaw += soillayer.wmass / (soillayer.wmass + soillayer.imass) * soillayer.depth;
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
bool Soilcolumn::can_freeze() const {
	return ifallowphasechanges;
}

/// Calculate temperature profile through snow, soil and padding layers
/** Calculate temperature profile through snow, soil and padding layers
 */
void Soilcolumn::calc_temp() {

	double temp[NLAYERS_MAX + 1], 
		dz[NLAYERS_MAX + 1],
		ki[NLAYERS_MAX + 1],
		ci[NLAYERS_MAX + 1];

	// Air temperature 
	double airT = patch.get_climate().temp;

	// Number of sub timesteps to be calculated
	int stimesteps = 6;

	// If a snowpack and standing water is present then soak snowpack
	if (!negligible(wcont_ssw()) && !negligible(snow.depth())) {
		snow.soak_snowpack(wmass_ssw, imass_ssw, temp_ssw);
	}

	// Update heat capacity and thermal conductivity for snow, standing water, soil and padding layers
	update_thermal_properties();

	for (int ts = 0; ts < stimesteps; ts++) {

		int nsurflayers = 0;
		bool combine_surf_soil = false;

		// Layer above soil exists
		if (!negligible(snow.depth() + wcont_ssw())) {

			// Layer above soil is thinner that 50% of top soil layer
			// Should be combined with top soil layer for computational stability
			if (snow.depth() + wcont_ssw() <= soillayers[0].depth * 0.5 && nsurflayers < snow.nsnowlayers) {

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

		int nairlayers = 1;
		int nairsurflayers = nairlayers + nsurflayers;
		int nairsurfsoillayers = nairlayers + nsurflayers + nsoillayers;
		int nairsurfsoilpadlayers = nairlayers + nsurflayers + nsoillayers + NPADLAYERS;

		// Air layer
		dz[0] = AIR_THICKNESS;
		ci[0] = Cp_air * M3_PER_MM3;
		ki[0] = Kair * SECS_PER_DAY / MM_PER_M;
		temp[0] = airT;

		// Add additional layers
		if (nsurflayers) {
			if (snow.depth()) {
				int slyr = snow.nsnowlayers - 1;
				for (int lyr = nairlayers; lyr < nairsurflayers; lyr++) {

					// Snow layer
					Snowlayer& snowlayer = snow.snowlayers[slyr];

					dz[lyr] = snowlayer.depth;
					ci[lyr] = snowlayer.Ch / MM3_PER_M3; // m-3 to mm-3
					ki[lyr] = snowlayer.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
					temp[lyr] = snowlayer.temp;
					
					slyr--;
				}
			}
			else {

				// Standing water layers
				int lyr = nairsurflayers - 1;

				dz[lyr] = wcont_ssw();
				ci[lyr] = Ch_ssw / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = Kh_ssw * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				temp[lyr] = temp_ssw;
			}
		}

		// Soil layers
		for (int lyr = nairsurflayers; lyr < nairsurfsoillayers; lyr++) {

			// Soil layer
			Soillayer& soillayer = soillayers[lyr - nairsurflayers];

			// Combine snow and standing water with top soil layer
			if (lyr == nairsurflayers && combine_surf_soil) {
				
				// Additional layer to add to top soil layer
				double depth_surf, Ch_surf, Kh_surf, temp_surf;
				if (snow.depth()) {
					Snowlayer& snowlayer = snow.snowlayers[0];
					depth_surf = snowlayer.depth;
					Ch_surf = snowlayer.Ch;
					Kh_surf = snowlayer.Kh;
					temp_surf = snowlayer.temp;
				}
				else {
					depth_surf = wcont_ssw();
					Ch_surf = Ch_ssw;
					Kh_surf = Kh_ssw;
					temp_surf = temp_ssw;
				}
				

				dz[lyr] = (soillayer.depth + depth_surf);
				ci[lyr] = (soillayer.Ch * soillayer.depth + Ch_surf * depth_surf) / dz[lyr] / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = (soillayer.Kh * soillayer.depth + Kh_surf * depth_surf) / dz[lyr] * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				temp[lyr] = (soillayer.temp * soillayer.depth * soillayer.Ch + temp_surf * depth_surf * Ch_surf) / 
					(soillayer.depth * soillayer.Ch + depth_surf * Ch_surf);
			}
			else {
				dz[lyr] = soillayer.depth;
				ci[lyr] = soillayer.Ch / MM3_PER_M3; // m-3 to mm-3
				ki[lyr] = soillayer.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1 
				temp[lyr] = soillayer.temp;
			}
		}

		// Padding layers
		for (int lyr = nairsurfsoillayers; lyr < nairsurfsoilpadlayers; lyr++) {

			// Padding layer
			padding& pad_lyr = pad[lyr - nairsurfsoillayers];

			dz[lyr] = pad_lyr.depth;
			ci[lyr] = pad_lyr.Ch / MM3_PER_M3; // m-3 to mm-3
			ki[lyr] = pad_lyr.Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
			temp[lyr] = pad_lyr.temp;
		}

		// Calculated temperatures using the Crank-Nicholson timestepper algorithm 
		// for temperature diffusion equation
		cnstep_full(airT, 1.0 / (double)stimesteps, nairsurfsoilpadlayers, temp, dz, ki, ci);

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
		if (!negligible(wcont_ssw())) {
			temp_ssw = temp[1];
		}

		// Set new soil layer temperatures
		for (int lyr = nairsurflayers; lyr < nairsurfsoillayers; lyr++) {
			soillayers[lyr - nairsurflayers].temp = temp[lyr];
		}

		// Set new padding layer temperatures
		for (int lyr = nairsurfsoillayers; lyr < nairsurfsoilpadlayers; lyr++) {
			pad[lyr - nairsurfsoillayers].temp = temp[lyr];
		}

		// Phase change between water and ice and adjust snow, standing water, 
		// and soil layer temperatures accordingly
		phase_change();
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	// Update snow properties
	snow.balance_snow(patch, rain_melt, rain_melt_temp);

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
			pad[lyr].temp = avr_temp;
		}
	}
}


/// Calculate snow, soil and padding layer temperatures
/** Calculate snow, soil and padding layer temperatures
*/
void Soilcolumn::temp() {

	if (hydrology == GERTEN) {
		// Calculate analytic soil temperatures and update dtemp_31 (trend in daily air temperatures for the last 31 days)
		// Analytic soil temperatures can be used to test the accuracy of the numerical scheme in Soilcolumn::calc_temp below.
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
void Soilcolumn::calc_msnow() {

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
void Soilcolumn::calc_root_volume() {

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
			double volume_coarse_root_tot = (ppft.litter_sap + ppft.litter_heart) * (1.0 - pft.ag_woody_frac) /
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
		soillayer.root_volume = min(0.5 * soillayer.depth,
			max(0.0, volume_coarse_root[lyr]) / rho_high * MM_PER_M);

		// If organic matter doesn't affect soil properties, then they 
		// need to be updated here if effective volume is updated
		if (!iforganicsoilproperties) {
			double eff_volume_frac = 1.0 - soillayer.root_volume / (soillayer.depth * M_PER_MM);
			soillayer.wp = soiltype.wp_min * eff_volume_frac * soillayer.depth;
			soillayer.awc = soiltype.awc_min * eff_volume_frac * soillayer.depth;
			soillayer.wsats = soiltype.wsats_min * eff_volume_frac * soillayer.depth;
		}
	}
}

/// Vertical movement of organic matter if C density of soil layer exceed maximum
/** Vertical movement of organic matter if C density of soil layer exceed maximum
*  Downwards movement of SOM if soil layer C exceeds maximum
*/
void Soilcolumn::vertical_movement() {

	double ccont_before = ccont();
	double wcont_before = get_sum(WCONTENT);
	double wadd = 0.0;

	double organic_frac_move_thres = 0.001;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Initial mineral volume
		double min_volume = soillayer.mineral_volume;

		// Calculate soil layer organic carbon masses
		double ccont_humus = soillayer.ccont(GHUMUS);
		double ccont_som = soillayer.ccont(GSOM);

		// Calculate soil layer organic carbon densities
		double rho_humus = soillayer.calc_organic_density(GHUMUS);
		double rho_som = soillayer.calc_organic_density(GSOM);

		// Update organic volume (mm)
		double humus_volume = ccont_humus / rho_humus * MM_PER_M;
		double som_volume = ccont_som / rho_som * MM_PER_M;
		double org_volume = humus_volume + som_volume;
		soillayer.organic_volume = org_volume;

		// Effective soil volume
		double eff_volume = soillayer.depth - soillayer.root_volume;

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

					double organic_frac_int = org_volume / (org_volume + min_volume);

					// For peat move mineral soil first
					if (patch.stand.landcover == PEATLAND) {
						organic_frac_int = 0.0;
					}

					// Fraction of matter to be moved
					double frac;

					// Move mineral nitrogen
					frac = min(1.0, diff_volume / (min_volume + org_volume));
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

					// At very low organic matter fraction at the interface 
					// between the layers, then only mineral volume should be moved
					if (organic_frac_int < organic_frac_move_thres) {

						// Check if mineral matter is enough to fulfill move
						if (diff_volume <= min_volume) {
							frac = diff_volume / min_volume;
						}
						else {
							// Else all mineral volume will be moved
							frac = 1.0;
						}

						// Move mineral matter
						double min_move = min_volume * frac;
						diff_volume -= min_move;
						min_volume -= min_move;

						soillayer_below.mineral_volume += min_move;
						soillayer.mineral_volume -= min_move;
					}

					// Still volume to be moved
					if (diff_volume > 1.0e-15) {

						// Check if mineral and som matter is enough to fulfill move
						double frac_min, frac_som;
						if (!negligible(min_volume + som_volume)) {

							// Only organic matter to move
							if (negligible(min_volume)) {
								frac_min = 0.0;
								frac_som = min(1.0, diff_volume / som_volume);
							}
							// Only mineral matter to move
							else if (negligible(som_volume)) {
								frac_min = min(1.0, diff_volume / min_volume);
								frac_som = 0.0;
							}
							// All mineral and organic matter needs to move. Only humus left
							else if (diff_volume > min_volume + som_volume) {
								frac_min = 1.0;
								frac_som = 1.0;
							}
							else {
								// Calculate mineral volume to be moved with interface fraction
								double diff_volume_som = diff_volume * organic_frac_int;
								double diff_volume_min = diff_volume - diff_volume_som;

								if (diff_volume_min <= min_volume && diff_volume_som <= som_volume) {
									frac_min = diff_volume_min / min_volume;
									frac_som = diff_volume_som / som_volume;
								}
								else if (diff_volume_som > som_volume) {
									frac_som = 1.0;
									frac_min = (diff_volume - som_volume) / min_volume;
								}
								else {
									frac_som = (diff_volume - min_volume) / som_volume;
									frac_min = 1.0;
								}
							}

							// Move som organic matter
							for (int p = 0; p < NSOMPOOL; p++) {
								if (p == SOILMETA || p == SOILSTRUCT || p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
									soillayer_below.sompool[p].cmass += soillayer.sompool[p].cmass * frac_som;
									soillayer_below.sompool[p].nmass += soillayer.sompool[p].nmass * frac_som;
									soillayer.sompool[p].cmass *= 1.0 - frac_som;
									soillayer.sompool[p].nmass *= 1.0 - frac_som;
								}
							}

							// Move mineral matter
							soillayer_below.mineral_volume += min_volume * frac_min;
							soillayer.mineral_volume *= 1.0 - frac_min;

							// Subtract moved matter
							diff_volume -= min_volume * frac_min + som_volume * frac_som;
						}

						// Still volume to be moved
						if (diff_volume > 1.0e-15 && !negligible(humus_volume)) {

							double frac_humus;
							if (diff_volume <= humus_volume) {
								frac_humus = diff_volume / humus_volume;
							}
							else {
								frac_humus = 1.0;
							}

							// Move som organic matter
							for (int p = 0; p < NSOMPOOL; p++) {
								if (p == SURFMICRO || p == SURFHUMUS) {
									soillayer_below.sompool[p].cmass += soillayer.sompool[p].cmass * frac_humus;
									soillayer_below.sompool[p].nmass += soillayer.sompool[p].nmass * frac_humus;
									soillayer.sompool[p].cmass *= 1.0 - frac_humus;
									soillayer.sompool[p].nmass *= 1.0 - frac_humus;
								}
							}

							diff_volume -= humus_volume * frac_humus;
						}

						if (diff_volume > 1.0e-12) {
							dprintf("Year %d Day %d Lyr %d ERROR in vertical movement! No more volume to move (%g).\n", date.year, date.day, lyr, diff_volume);
						}
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
						double ccont_humus = soillayer_below.ccont(GHUMUS);
						double ccont_som = soillayer_below.ccont(GSOM);

						// Calculate soil layer organic carbon densities
						double rho_humus = soillayer_below.calc_organic_density(GHUMUS);
						double rho_som = soillayer_below.calc_organic_density(GSOM);

						// Update organic volume (mm)
						double humus_volume = ccont_humus / rho_humus * MM_PER_M;
						double som_volume = ccont_som / rho_som * MM_PER_M;
						double org_volume = humus_volume + som_volume;
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

						// Organic matter fraction at the interface between the layers
						double organic_frac_int = 0.0;

						if (!negligible(org_volume + soillayer_below.mineral_volume)) {
							organic_frac_int = org_volume / (org_volume + soillayer_below.mineral_volume);
						}

						// For peat move organic matter up first
						if (patch.stand.landcover == PEATLAND) {
							organic_frac_int = 1.0;
						}

						// Move mineral nitrogen
						double frac = min(1.0, up_volume / volume_below);
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

						if (!negligible(org_volume) || !negligible(soillayer_below.mineral_volume)) {

							// Humus should be moved upwards first
							if (!negligible(humus_volume)) {
								double frac_humus;
								if (up_volume <= humus_volume) {
									frac_humus = up_volume / humus_volume;
								}
								else {
									frac_humus = 1.0;
								}

								// Move som organic matter
								for (int p = 0; p < NSOMPOOL; p++) {
									if (p == SURFMICRO || p == SURFHUMUS) {
										soillayer.sompool[p].cmass += soillayer_below.sompool[p].cmass * frac_humus;
										soillayer.sompool[p].nmass += soillayer_below.sompool[p].nmass * frac_humus;
										soillayer_below.sompool[p].cmass *= 1.0 - frac_humus;
										soillayer_below.sompool[p].nmass *= 1.0 - frac_humus;
									}
								}

								up_volume -= humus_volume * frac_humus;
							}

							// Still volume needed to move upwards
							if (!negligible(up_volume)) {

								// Check if mineral and organic matter is enough to fulfill move
								double frac_min, frac_org;
								if (up_volume > soillayer_below.mineral_volume + som_volume) {
									frac_min = 1.0;
									frac_org = 1.0;
								}
								else if (negligible(som_volume)) {
									frac_min = min(1.0, up_volume / soillayer_below.mineral_volume);
									frac_org = 0.0;
								}
								else if (negligible(soillayer_below.mineral_volume)) {
									frac_min = 0.0;
									frac_org = min(1.0, up_volume / som_volume);
								}
								else {
									// Calculate mineral volume to be moved with interface fraction
									double up_volume_org = up_volume * organic_frac_int;
									double up_volume_min = up_volume - up_volume_org;

									if (up_volume_min <= soillayer_below.mineral_volume &&
										up_volume_org <= som_volume) {
										frac_min = up_volume_min / soillayer_below.mineral_volume;
										frac_org = up_volume_org / som_volume;
									}
									else if (up_volume_org > som_volume) {
										frac_org = 1.0;
										frac_min = (up_volume - som_volume) / soillayer_below.mineral_volume;
									}
									else {
										frac_org = (up_volume - soillayer_below.mineral_volume) / som_volume;
										frac_min = 1.0;
									}
								}

								// Move organic matter
								for (int p = 0; p < NSOMPOOL; p++) {
									if (p == SOILMETA || p == SOILSTRUCT || p == SOILFWD || p == SOILCWD || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
										soillayer.sompool[p].cmass += soillayer_below.sompool[p].cmass * frac_org;
										soillayer.sompool[p].nmass += soillayer_below.sompool[p].nmass * frac_org;
										soillayer_below.sompool[p].cmass *= 1.0 - frac_org;
										soillayer_below.sompool[p].nmass *= 1.0 - frac_org;
									}
								}

								// Move mineral matter
								soillayer.mineral_volume += soillayer_below.mineral_volume * frac_min;
								soillayer_below.mineral_volume *= 1.0 - frac_min;
							}
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
				soillayer.depth += add_depth;

				// Update soil layer depth profile
				soillayer.depth_profile += add_depth / 2.0;

				// Update soil columns array of soil layer depths
				soillayerdepths[lyr] += add_depth;

				// Update soil columns array of soil layer thicknesses
				soillayerthicks[lyr] += add_depth;

				// Update total soil column depth
				depth += add_depth;

				// Update soil layer mineral volume
				soillayer.mineral_volume = eff_volume - soillayer.organic_volume;

				double wsats_new;
				if (soillayer.acrotelm) {
					wsats_new = soillayer.mineral_volume * acrotelm_por;
				}
				else if (soillayer.catotelm) {
					wsats_new = soillayer.mineral_volume * catotelm_por;
				}
				else {
					wsats_new = soillayer.mineral_volume * soiltype.wsats_min;
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
	double wcont_after = get_sum(WCONTENT);

	if (fabs(ccont_after - ccont_before) > 1.0e-10) {
		dprintf("Year %d ERROR VERTMOVE C before %g after %g Diff %g\n", date.year, ccont_before, ccont_after, ccont_before - ccont_after);
	}

	if (fabs((wcont_after - wadd) - wcont_before) > 1.0e-10) {
		dprintf("Year %d ERROR VERTMOVE W before %g after %g Diff %g\n", date.year, wcont_before, wcont_after, wcont_before - wcont_after);
	}
}

/// Serialize soilcolumn variables
/** Serialize soilcolumn variables
 */
void Soilcolumn::serialize(ArchiveStream& arch) {
	arch & depth
		& nsoillayers
		& wmass_ssw
		& imass_ssw
		& temp_ssw
		& nacrotelm
		& ncatotelm
		& snow
		& maxthaw_10
		& slope_grade_sin
		& slope_grade_cos
		& slope_sigma
		& sigma_sqrt2
		& sigma_sqrt2pi
		& sigma2
		& wmass_ssw_frac_thresh
		& wmass_ssw_thresh;

	for (int i = 0; i < NSOILLAYERS_MAX; i++) {
		arch& soillayerdepths[i];
		arch& soillayerthicks[i];
	}

	for (int i = 0; i < nsoillayers; i++) {
		arch & soillayers[i];
	}

	for (int i = 0; i < NPADLAYERS; i++) {
		arch & pad[i];
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// References not related to Arctic and wetland code:
//
// Campbell, G.S.: Simple method for determining unsaturated conductivity from moisture 
//  retention data. Soil Sci. 117: 311–314. 
//  https://doi.org/10.1097/00010694-197406000-00001, 1974
// Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., 
//   and Friedlingstein, P.: An improved representation of physical permafrost dynamics 
//   in the JULES land - surface model, Geosci.Model Dev., 8, 1493 - 1508, 
//   https ://doi.org/10.5194/gmd-8-1493-2015, 2015
// Clapp, R.B., and Hornberger, G.M. Empirical equations for some soil hydraulic properties. 
//  Water Resour. Res. 14: 601–604. https://doi.org/10.1029/WR014i004p00601, 1978
// Fukusako, S. Int J Thermophys (1990) 11: 353. doi:10.1007/BF01133567
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//   and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Lawrence, D. M. and Slater, A.: Incorporating organic soil into a global climate 
//  model, Clim.Dynam., 30, 145–160, doi:10.1007/s00382-007-0278-1, 2008.
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
