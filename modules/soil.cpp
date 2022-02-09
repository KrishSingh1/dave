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
	Dh = 0.0;
}

/// Calculates heat capacity of snow
/** Calculates heat capacity of snow
 *  Best et al. 2011, Eq. 22 [J m-3 K-1]
 */
double snowheatcap(double imass, double wmass, double depth) {
	return (Cp_ice * imass + Cp_water * wmass) / depth;
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

	// Snow density including water withing snow
	double tot_dens = min(rho_ICE, (imass + wmass) / depth * MM_PER_M);

	// Heat capacity of snow - Best et al. 2011, Eq. 22 [J m-3 K-1]
	Ch = snowheatcap(imass, wmass, depth);

	// Thermal conductivity of snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = 2.22 * pow((tot_dens / rho_H2O), 1.88);

	// Diffusivity of snow [m2 s-1]
	Dh = Kh / Ch;
}

/// Update snow water holding capacity
/** Update snow water holding capacity
 * based on eqn 10 Wang et al. 2013
 */
void Snowlayer::update_whc() {

	// Reference density [kg m-3]
	double dens_t = 400;	
	// Minimum snow water holding capacity
	double rw_min = 0.03;
	// Maximum snow water holding capacity
	double rw_max = 0.1;

	// Water holding capacity of snow layer [mm]
	whc = max(0.0, min(depth - imass, imass * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - dens) / dens_t)));
}

/// Weigthed avg of new and old snow density
/** Weigthed avg of new and old snow density
 */
void Snowlayer::mix_dens(double prec, double snowfall_dens) { // DW_COMMENT not used!

	// Total water/ice mass
	double sum = imass + wmass + prec;

	// New water/ice mass fraction
	double w_new = prec / sum;

	// Old water/ice mass fraction
	double w_old = 1.0 - w_new;

	// New density
	dens = min(rho_ICE, w_old * dens + w_new * snowfall_dens);
}

/** Phase change of snow layer
 */
void Snowlayer::phase_change(bool single_layer) {

	// Fraction of water which freezes
	double Ffreez;
	// Fraction of ice which thaws
	double Fthaw;
	// Energy per volume [J m-3]
	double energy;
	// Temperature difference (current - previous)
	double delta_T;
	// Phase change energy consumption [J]
	double phase_change_energy = 0.0;

	bool iffreezed = false;
	bool ifthawed = false;

	// The total amount of ice and water should not change in this routine
	double initial_wimass = imass + wmass;

	double lyr_depth = single_layer ? max(depth, SNOW_DEPTH_STANDARD) : depth;

	// Water and ice fraction of full volume
	double Fwmass = wmass / lyr_depth;
	double Fimass = imass / lyr_depth;

	// Freezing temperature and water present
	if (temp < 0.0 && wmass > 0.0) {

		iffreezed = true;

		// Temperature change from zero
		delta_T = fabs(temp);

		// Energy to freeze water [J m-3]
		energy = Cp_water * delta_T;

		// Fraction water that would freeze [m3 m-3]
		Ffreez = energy / Lheat;

		// More water in a volume than could freeze
		if (Fwmass > Ffreez) {

			// There is some water left, even after freeze
			Fimass += Ffreez;
			Fwmass -= Ffreez;

			// Energy from daily temperature change is used for phase change
			temp = 0.0;
		}
		else {

			// It is cold enough to freeze all the water and decrease the temperature below 0,
			// since the rest of the water gets frozen now. However, it should never result in
			// temperatures colder than the values coming from the Crank-Nicholson scheme.

			// Snow layer temperature rise resulting from freezing water.
			double temp_rise = Fwmass * Lheat / Cp_water;

			// Update snow temperature in this layer
			temp += temp_rise;
			temp = min(0.0, temp);

			// All water is frozen
			Fimass += Fwmass;
			Fwmass = 0.0;
		}
	} // End of freezing

	// Thawing temperatures and ice present
	else if (temp >= 0.0 && imass > 0.0) {

		ifthawed = true;

		// Temperature change from zero
		delta_T = fabs(temp);

		// Energy to thaw ice [J m-3]
		energy = Cp_water * delta_T;

		// Fraction ice that would thaw [m3 m-3]
		Fthaw = energy / Lheat;

		// More ice in a volume than could thaw
		if (Fimass > Fthaw) {

			// There is some ice left, even after melt
			Fwmass += Fthaw;
			Fimass -= Fthaw;

			// Energy from daily temperature change is used for phase change
			temp = 0.0;
		}
		else {
			// There is more energy available than required for
			// melting the rest of the ice. Therefore no need to
			// calculate new temperature as snow layer will be removed

			// All ice thaws
			Fwmass += Fimass;

			// No ice left in snow layer
			Fimass = 0.0;
		}

	} // End of thawing

	// New ice fraction in snow layer
	imass = Fimass * lyr_depth;

	// New water fraction in snow layer
	wmass = Fwmass * lyr_depth;

	// Freezing - Update density as depth doesn't change
	if (iffreezed) {
		if (!negligible(depth)) {
			dens = imass * rho_H2O / depth;
		}
	}
	// Thawing - Update depth as density doesn't change
	else if (ifthawed) {
		depth = imass / (dens / rho_H2O);
	}
}

/** Phase change of snow layer
 */
void Snowlayer::phase_change_new(bool& updated) {

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

	// Water and ice fraction of full volume
	Fwmass = wmass / depth;
	Fimass = imass / depth;

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
	else {
		updated = true;
	}

	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy from full soil layer to freeze water [J m-3]
		energy = Ch * delta_T;

		// Fraction ice that can be frozen [m3 m-3]
		Ffreez = energy / ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

		if (Fwmass > Ffreez) {

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
			wmass = (Fwmass - Fwmass) * depth;
			imass = (Fimass + Fwmass) * depth;

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
			wmass = (Fwmass + Fimass) * depth;
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

	// Freezing - Update density as depth doesn't change
	if (ifFreeze) {
		if (!negligible(depth)) {
			dens = imass * rho_H2O / depth;
		}
	}
	// Thawing - Update depth as density doesn't change
	else if (ifThaw) {
		depth = imass / (dens / rho_H2O);
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
	if ((fabs(energy_initial - energy_final)) > energy_initial * EPS) {
		dprintf("Year %d Day %d Energy balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", date.year, date.day, energy_initial, energy_final, energy_initial - energy_final);
	}
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
	for (int lyr = 0; lyr < nsnowlayers; lyr++)
		wcontent += snowlayers[lyr].wmass + snowlayers[lyr].imass;
	return wcontent;
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
	double snow_max_mths[3];
	snow_max_mths[0] = snow_max_mths[1] = snow_max_mths[2] = 0.0;

	for (int mth = 0; mth < 12; mth++) {

		// If current element is greater than first
		if (mdepth[mth] > snow_max_mths[0]) {
			snow_max_mths[2] = snow_max_mths[1];
			snow_max_mths[1] = snow_max_mths[0];
			snow_max_mths[0] = mdepth[mth];
		}
		// If mdepth[mth] is in between first and second then update second
		else if (mdepth[mth] > snow_max_mths[1]) {
			snow_max_mths[2] = snow_max_mths[1];
			snow_max_mths[1] = mdepth[mth];
		}
		// If mdepth[mth] is larger than the third then update third
		else if (mdepth[mth] > snow_max_mths[2])
			snow_max_mths[2] = mdepth[mth];
	}

	// Return avergae snow depth of the three months of highest snow depth
	return (snow_max_mths[0] + snow_max_mths[1] + snow_max_mths[2]) / 3.0;
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
void Snow::redist_layers(double snowpack_depth, double snowfall_depth, double snowfall_imass, double snowfall_temp, double snowfall_Ch) {

	const double EPS = 1.0e-12;
	int lyr;

	// Save current snowpack status
	double snowlyr_depth[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_imass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_wmass[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_temp[NSNOWLAYERS_MAX] = { 0.0 };
	double snowlyr_Ch[NSNOWLAYERS_MAX] = { 0.0 };

	// Water in snowlayer that has no ice after phase change, hence no depth (a full layer has melted)
	double wmass_no_snowlyr = 0.0;

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlyr_depth[lyr] = snowlayers[lyr].depth;
		snowlyr_imass[lyr] = snowlayers[lyr].imass;
		if (snowlayers[lyr].depth > 0.0) {
			snowlyr_wmass[lyr] = snowlayers[lyr].wmass;
		}
		else {
			wmass_no_snowlyr += snowlayers[lyr].wmass;
		}

		snowlyr_temp[lyr] = snowlayers[lyr].temp;
		snowlyr_Ch[lyr] = snowlayers[lyr].Ch;

		snowlayers[lyr].depth = snowlayers[lyr].imass = snowlayers[lyr].wmass = snowlayers[lyr].temp = snowlayers[lyr].Ch = 0.0;
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

	// Re-calc ice & dens of layers
	for (lyr = 0; lyr < nsnowlayers; lyr++) {

		diff_lyr = pot_depth[lyr] - snowlayers[lyr].depth;

		lyr_add = 0;
		diff = true;

		// Fill every snow layer to pot_ice from bottom up
		while (diff && lyr_add < nsnowlayers) {

			if (!negligible(snowlyr_depth[lyr_add])) {

				// Frac of lyr_add to add to this snow layer
				frac_move = min(1.0, diff_lyr / snowlyr_depth[lyr_add]);

				// Update snow layer temperature
				if (negligible(snowlayers[lyr].depth)) {
					snowlayers[lyr].temp = snowlyr_temp[lyr_add];
				}
				else {
					snowlayers[lyr].update_temp(snowlyr_temp[lyr_add], frac_move * (snowlyr_imass[lyr_add] + snowlyr_wmass[lyr_add]), snowlyr_Ch[lyr_add]);
				}

				// Add snow
				snowlayers[lyr].depth += frac_move * snowlyr_depth[lyr_add];
				snowlayers[lyr].imass += frac_move * snowlyr_imass[lyr_add];
				snowlayers[lyr].wmass += frac_move * snowlyr_wmass[lyr_add];

				// Calculate new density
				snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;

				// Update heat capacity
				snowlayers[lyr].Ch = snowheatcap(snowlayers[lyr].imass, snowlayers[lyr].wmass, snowlayers[lyr].depth);

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

			// Frac of snowfall to add to this snow layer
			frac_move = min(1.0, max(0.0, pot_depth[lyr] - snowlayers[lyr].depth) / snowfall_depth);

			// Update snow layer temperature
			if (negligible(snowlayers[lyr].depth)) {
				snowlayers[lyr].temp = snowfall_temp;
			}
			else {

				snowlayers[lyr].update_temp(snowfall_temp, frac_move * snowfall_imass, snowfall_Ch);
			}

			// Add snow
			snowlayers[lyr].depth += frac_move * snowfall_depth;
			snowlayers[lyr].imass += frac_move * snowfall_imass;

			// Calculate new density
			snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;

			// Update heat capacity
			snowlayers[lyr].Ch = snowheatcap(snowlayers[lyr].imass, snowlayers[lyr].wmass, snowlayers[lyr].depth);

			// Subtract added snowfall
			snowfall_depth *= (1.0 - frac_move);
			snowfall_imass *= (1.0 - frac_move);
		}
		else {

			// Add remaining (negligible) snowfall to the top snow layer
			int top_lyr = nsnowlayers - 1;
			
			// Update snow layer temperature
			snowlayers[top_lyr].update_temp(snowfall_temp, snowfall_imass, snowfall_Ch);

			// Add snow
			snowlayers[top_lyr].depth += snowfall_depth;
			snowlayers[top_lyr].imass += snowfall_imass;

			// Calculate new density
			snowlayers[top_lyr].dens = snowlayers[top_lyr].imass * rho_H2O / snowlayers[top_lyr].depth;

			// Update heat capacity
			snowlayers[top_lyr].Ch = snowheatcap(snowlayers[top_lyr].imass, snowlayers[top_lyr].wmass, snowlayers[top_lyr].depth);

			// Remove snowfall
			snowfall_depth = 0.0;
			snowfall_imass = 0.0;
		}
	}

	// Check if snow layers should be erased
	if (nsnowlayers > pot_nsnowlayers) {
		for (lyr = nsnowlayers - 1; lyr >= pot_nsnowlayers; lyr--) {
			if (!negligible(snowlayers[lyr].wmass + snowlayers[lyr].imass)) {
				snowlayers[pot_nsnowlayers - 1].update_temp(0.0, snowlayers[lyr].wmass + snowlayers[lyr].imass, Cp_water);
				snowlayers[pot_nsnowlayers - 1].wmass += snowlayers[lyr].wmass + snowlayers[lyr].imass;
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

		// Fraction to remove
		double frac = min(1.0, depth_reduction / snowlayers[blyr].depth);

		// Reduce reduction
		depth_reduction -= frac * snowlayers[blyr].depth;

		// Runoff water
		runoff_water += frac * snowlayers[blyr].wmass;

		// Runoff ice
		runoff_ice += frac * snowlayers[blyr].imass;

		// Reduce snow
		snowlayers[blyr].wmass *= (1.0 - frac);
		snowlayers[blyr].imass *= (1.0 - frac);
		snowlayers[blyr].depth *= (1.0 - frac);

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

	// Balance layers depending on if snowpack is over maximum or not
	double snowpack_depth = depth();

	// Too thick snowpack (max 20m), remove some ice/water as runoff to 
	// not have a runaway effect of snowpack depth
	double max_snowpack_depth = 20000.0;
	if (snowpack_depth > max_snowpack_depth) {
		reduce_snowpack(patch, snowpack_depth - max_snowpack_depth);
		snowpack_depth = max_snowpack_depth;
	}

	// Redistribute snow layers if snow exist
	if (!negligible(snowpack_depth)) {
		redist_layers(snowpack_depth);
	}
	// Insignificant snow depth. Remove all snow layers and turn it into melt water
	else if (nsnowlayers) {
		for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {
			double added_wmass = snowlayers[lyr].imass + snowlayers[lyr].wmass;
			if (!negligible(added_wmass)) {
				rain_melt_temp = (rain_melt_temp * rain_melt + snowlayers[lyr].temp * added_wmass) / (rain_melt + added_wmass);
				rain_melt += added_wmass;
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

	if (ifsoaksnow) {

		// Soak snowpack with melted water
		// Water fill full snow volume. Water above 
		// total volume is added to rain_melt
		soak_snowpack(ts_melt, ts_melt_temp);
	}
	else {

		// Redistribute water within the snow pack
		// moving from top to bottom snow layers
		// Water above whc is added to rain_melt
		for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

			if (!negligible(ts_melt)) {
				snowlayers[lyr].update_temp(ts_melt_temp, ts_melt, Cp_water);
				snowlayers[lyr].wmass += ts_melt;
			}

			if (snowlayers[lyr].wmass > snowlayers[lyr].whc) {
				ts_melt = max(0.0, snowlayers[lyr].wmass - snowlayers[lyr].whc);
				snowlayers[lyr].wmass -= ts_melt;
				ts_melt_temp = max(0.0, snowlayers[lyr].temp);
			}
			else {
				ts_melt = 0.0;
				ts_melt_temp = 0.0;
			}
		}
	}

	// Add excess water this timesteps melt to surface melt-water to pass to soil on next timestep
	if (!negligible(ts_melt)) {
		rain_melt_temp = (rain_melt_temp * rain_melt + ts_melt_temp * ts_melt) / (rain_melt + ts_melt);
		rain_melt += ts_melt;
	}
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

	bool single_layer = nsnowlayers == 1 ? true : false;

	// Phase change for each snow layer
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].phase_change(single_layer);
	}
}
void Snow::phase_change_new(bool& updated) {

	// No phase change in simple snow scheme
	if (!ifmultilayersnow)
		return;

	// Phase change for each snow layer
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].phase_change_new(updated);
	}
}

/// Soak snow pack with water
/** Soak snow pack with water
 * Water soaks snowpack from bottom up
 */
void Snow::soak_snowpack(double& ts_melt, double ts_melt_temp) {

	if (!ifsoaksnow || nsnowlayers == 0)
		return;

	// Add rain water on top snow layer
	if (!negligible(ts_melt)) {
		snowlayers[nsnowlayers - 1].update_temp(ts_melt_temp, ts_melt, Cp_water);
		snowlayers[nsnowlayers - 1].wmass += ts_melt;
		ts_melt = 0.0;
		ts_melt_temp = 0.0;
	}

	// Final water content in each snow layer
	std::vector<double> wmass_lyr(nsnowlayers);

	// If there is water movement within snowpack
	bool wmass_movement = false;

	// Water moving between snow layers
	double excess_wmass = 0.0;

	// Loop through layers and fill layer up to whc
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {
		wmass_lyr[lyr] = snowlayers[lyr].wmass;

		if (wmass_lyr[lyr] > snowlayers[lyr].whc) {
			excess_wmass += wmass_lyr[lyr] - snowlayers[lyr].whc;
			wmass_lyr[lyr] = snowlayers[lyr].whc;
			wmass_movement = true;
		}
		else if (excess_wmass > 0.0) {
			double add_wmass = min(excess_wmass, snowlayers[lyr].whc - wmass_lyr[lyr]);
			wmass_lyr[lyr] += add_wmass;
			excess_wmass -= add_wmass;
		}
	}

	// If there has been water movement, then we need to resolve this
	if (wmass_movement) {

		// If there was water exceeding whc, then it will fill the snowpack from bottom up
		if (excess_wmass > 0.0) {

			// Maximum excess water the snowpack can hold when it is soaked
			double EXCESS_WMASS_MAX = 500.0;
			excess_wmass = min(EXCESS_WMASS_MAX, excess_wmass);

			int lyr = 0;
			while (excess_wmass > 0.0 && lyr < nsnowlayers) {
				double add_wmass = max(0.0, min(excess_wmass, snowlayers[lyr].depth - (snowlayers[lyr].whc + snowlayers[lyr].imass)));
				wmass_lyr[lyr] += add_wmass;
				excess_wmass -= add_wmass;
				lyr++;
			}
		}

		// Recalculate newwater conent in each snow layerand adjust temperature accordingly
		for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

			// Add water that has moved from snow layer above
			if (!negligible(ts_melt)) {
				snowlayers[lyr].update_temp(ts_melt_temp, ts_melt, Cp_water);
				snowlayers[lyr].wmass += ts_melt;
				ts_melt = 0.0;
				ts_melt_temp = 0.0;
			}

			// Check if water in snow layer is above final content
			if (snowlayers[lyr].wmass - wmass_lyr[lyr] > 1.0e-12) {
				ts_melt = snowlayers[lyr].wmass - wmass_lyr[lyr];
				snowlayers[lyr].wmass -= ts_melt;
				ts_melt_temp = max(0.0, snowlayers[lyr].temp);
			}
			else {
				// Check if water in snow layer is below final content
				int lyr_below = lyr - 1;
				while (snowlayers[lyr].wmass - wmass_lyr[lyr] < -1.0e-12 && lyr_below >= 0) {
					double add_wmass = min(snowlayers[lyr_below].wmass, wmass_lyr[lyr] - snowlayers[lyr].wmass);
					snowlayers[lyr].wmass += add_wmass;
					snowlayers[lyr_below].wmass -= add_wmass;
					lyr_below--;
				}
			}

			if (fabs(snowlayers[lyr].wmass - wmass_lyr[lyr]) > 1.0e-10)
				dprintf("ERROR snow layer water: Year %d Day %d Lyr %d Diff %g Current %g Should %g\n",date.year, date.day, lyr, snowlayers[lyr].wmass - wmass_lyr[lyr], snowlayers[lyr].wmass, wmass_lyr[lyr]);
		}
	}
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
	double rho_0 = 50.0;

	// Snow mass above current layer
	double I_above = 0.0;

	// Mechanical compaction scheme, following Best et al. 2011
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		// Current layers ice and water mass
		double I_sum = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		// Cumulative snow mass affecting compaction
		double M_k = 0.5 * I_sum + I_above;

		// Layer temperature in Kelvin (never above zero degree Celcius)
		double T_k = min(K2degC, snowlayers[lyr].temp + K2degC);

		if (snowlayers[lyr].dens < rho_ICE) {
			double delta_dens = snowlayers[lyr].dens * gravity * M_k / eta_0 *
				exp((k_s / K2degC) - (k_s / T_k) - (snowlayers[lyr].dens / rho_0)) * SECS_PER_DAY;

			// Increase density
			snowlayers[lyr].dens += delta_dens;

			// Snow density can not exceed pure ice density
			snowlayers[lyr].dens = min(rho_ICE, snowlayers[lyr].dens);
		}

		// Add this layer to the above layers for future calculations
		I_above += I_sum;
	}

	// Calculate snow depth after compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);
	}
}

/// Snowfall
/** Snowfall
*  Calculate snowfall density and redistribute snow layers
*/
void Snow::snowfall(double airT, double& prec, double u10) {

	// Fresh snow density calculation, following eqn 1 in Vionnet et al (2012)
	// Constant [kg m-3]
	int a = 109;
	// Constant [kg m-3 K-1]
	int b = 6;
	// Constant [kg m-7/2 s-1/2]
	int c = 26;			 
	// Minimum snow density 
	double dens_min = 100;
	// Temperature of the melting point for water
	double Tmelt = 0;

	// Snow fall density
	double snowfall_dens = max(dens_min, a + b * (airT - Tmelt) + c * pow(u10, 0.5));

	// Depth of the newly fallen snow based on precipitation & density (Oleason 2013, Eq. 6.88)
	double snowfall_depth = prec / (snowfall_dens / rho_H2O);

	// Snow fall heat capacity
	double snowfall_Ch = snowheatcap(prec, 0.0, snowfall_depth);

	// Snowpack depth
	double snowpack_depth = depth();

	// Add snowfall and redistribute snowpack
	redist_layers(snowpack_depth, snowfall_depth, prec, airT, snowfall_Ch);

	// Reset precipitation
	prec = 0.0;
}

/// Handle rain on snow event and rainwater distribution within the snowpack
/** Handle rain on snow event and rainwater distribution within the snowpack
 * Calculate water, ice and air fraction in each layer
 * and snow phase change in layer based on the available energy
 */
void Snow::rainwater_melt(double airT, double& prec, double& rain_melt, double& rain_melt_temp) {

	// Water balance check // DW_COMMENT to be removed!
	double w_total_initial = prec + rain_melt +  wcontent();

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

		// Heat capcity of snow layer [J m-3 K-1]
		double heat_cap = snowlayers[lyr].Ch;
		// Volume of snow layer [m-3]
		double snow_volume = snowlayers[lyr].depth * M_PER_MM;
		// Energy needed to bring snow layer to 0C [J]
		double energy_to_zero = fabs(snowlayers[lyr].temp - 0.0) * heat_cap * snow_volume;

		// Enough energy to melt some snow
		if (rain_energy > energy_to_zero) {

			// Snow layer temp is set to zero and erase energy needed to bring up temperature
			rain_energy -= energy_to_zero;
			snowlayers[lyr].temp = 0.0;

			// ice melt [mm]
			double imass_melt = rain_energy / Lheat * MM_PER_M;

			// Enough energy to melt all ice (snow layer will be removed) 
			if (imass_melt > snowlayers[lyr].imass) {
				imass_melt -= snowlayers[lyr].imass;
				snowlayers[lyr].wmass += snowlayers[lyr].imass;
				snowlayers[lyr].imass = 0.0;
				rain_energy = imass_melt * Lheat / MM_PER_M;
			}
			// Only melt a frac of ice
			else {
				snowlayers[lyr].wmass += imass_melt;
				snowlayers[lyr].imass -= imass_melt;
				rain_energy = 0.0;
			}
		}
		// Only enough energy to change temperature
		else {
			snowlayers[lyr].temp += rain_energy / (snow_volume * heat_cap);
			rain_energy = 0.0;
		}

		// If now layer still exist update depth, water content and temperature
		if (!negligible(snowlayers[lyr].imass)) {

			// Recalculate depth
			snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);

			// Recalculate water holding capacity of snow layer
			snowlayers[lyr].update_whc();

			// Transferred water added to this layer
			double water_added = min(water_transferred, snowlayers[lyr].whc - snowlayers[lyr].wmass);

			water_transferred -= water_added;

			// Update temperature assuming water temperature is 0C
			if (water_added > 0.0) {
				snowlayers[lyr].update_temp(0.0, water_added, Cp_water);
			}

			snowlayers[lyr].wmass += water_added;			
		}
		else {

			// Snow layer will be erased, transfer all water to next snow layer or the soil
			water_transferred += snowlayers[lyr].wmass;

			// Erase this empty snow layer
			snowlayers.erase(snowlayers.end()-1);
			nsnowlayers--;

			// If this was the last snow layer, then exit loop
			if (nsnowlayers == 0)
				break;
		}

		// Go thorugh lyrs while rain_energy not 0
		lyr--;
	} // End of while loop

	// Still energy left, update melted water temperature
	double water_transferred_temp = 0.0;
	if (!negligible(rain_energy)) {
		water_transferred_temp = rain_energy / (Cp_water * water_transferred * M_PER_MM);
	}

	// If a snowpack exist, no water influx until snowpack is fully soaked
	if (!negligible(depth()) && !negligible(water_transferred))
		soak_snowpack(water_transferred, water_transferred_temp);

	// If water to be transferred to the soil
	if (!negligible(water_transferred)) {

		// Update rain_melt temperture
		rain_melt_temp = (rain_melt_temp * rain_melt + water_transferred_temp * water_transferred) / (rain_melt + water_transferred);

		// Save melt water above snow whc which will be pass downwards
		rain_melt += water_transferred;
	}

	// Check final water content
	double w_total_final = rain_melt + wcontent();
	const double EPS = 1.0e-12;
	//assert(fabs(w_total_initial - w_total_final) < EPS);
	if (fabs(w_total_final - w_total_initial) > 1.0e-10)
		dprintf("Year %d Day %d ERROR: Snow rainwater_melt diff %g\n", date.year, date.day, w_total_final - w_total_initial);
}

/// Snow accumulation
/** Snow accumulation
 * Rain or snowfall, and melt if simple snow scheme is used
 */
void Snow::snow_or_rainfall(double airT, double& prec, double u10, double& rain_melt, double& rain_melt_temp) {

	// No precipitation, then return from function
	if (prec == 0)
		return;

	double w_begin = wcontent() + prec + rain_melt;

	// Multi layer snow scheme
	if (ifmultilayersnow) {
		
		// Snowfall
		if (airT < 0.0) {

			// Add precipitation as snow and update depth, density
			snowfall(airT, prec, u10);
		}
		// Precipitation falls on snow
		else if (airT > 0.0 && !negligible(depth())) {

			// Phase changes and temperature change due to rainwater
			rainwater_melt(airT, prec, rain_melt, rain_melt_temp);
		}
		else {

			// Precipitation falls as rain

			// Update temperature
			rain_melt_temp = (rain_melt * rain_melt_temp + prec * airT) / (rain_melt + prec);

			// Add precipitation
			rain_melt += prec;
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
			if (nsnowlayers) {
				snowfall = min(prec, SNOWPACK_MAX - snowlayers[0].imass);
			}
			else {
				snowfall = min(prec, SNOWPACK_MAX);
			}

			rain_melt += max(0.0, prec - snowfall);

			// If no snow layers
			if (nsnowlayers == 0) {
				snowlayers.resize(nsnowlayers + 1);
				snowlayers[nsnowlayers] = Snowlayer(snowfall, 0.0, SNOWDENSITY, airT);
				nsnowlayers++;
			}
			else {
				snowlayers[0].imass += snowfall;
			}			
		}
		// Snow melting today
		else if (airT > 0.0 && nsnowlayers) {

			// Snow melt formulation (Choudhury et al 1998)
			double snowmelt = min((1.5 + 0.007 * prec) * (airT - 0.0), snowlayers[0].imass);

			snowlayers[0].imass -= snowmelt;

			rain_melt += prec + snowmelt;
		}
		else {

			// Add precipitation
			rain_melt += prec;
		}

		// Update single snowlayer with new snow or melt
		if (nsnowlayers) {
			
			snowlayers[0].wmass = 0.0;
			snowlayers[0].depth = snowlayers[0].imass / (snowlayers[0].dens / rho_H2O);

			if (negligible(snowlayers[0].imass)) {

				// Erase this empty snow layer
				snowlayers.erase(snowlayers.begin());
				nsnowlayers--;
			}
		}
	}

	// Water balance might not meet mass balance in very low productive and saturated sites 
	double wdiff = w_begin - (wcontent() + rain_melt);
	if (wdiff > 0.0) {
		rain_melt += wdiff;
	}
	else if (wdiff < 0.0) {
		rain_melt += wdiff;
		if (rain_melt < 0.0) {
			if (wcontent() > 0.0) {
				int lyr = nsnowlayers - 1;
				while (rain_melt < 0.0 && lyr >= 0) {
					snowlayers[lyr].wmass += rain_melt;
					rain_melt = 0.0;
					if (snowlayers[lyr].wmass < 0.0) {
						rain_melt = snowlayers[lyr].wmass;
						snowlayers[lyr].wmass = 0.0;
					}

					if (rain_melt < 0.0) {
						snowlayers[lyr].imass += rain_melt;
						rain_melt = 0.0;
						
						if (snowlayers[lyr].imass < 0.0) {
							rain_melt = snowlayers[lyr].imass;
							snowlayers[lyr].imass = 0.0;
							snowlayers.erase(snowlayers.begin() + lyr);
							nsnowlayers = (int)snowlayers.size();
						}
						else {
							snowlayers[lyr].dens = snowlayers[lyr].imass / snowlayers[lyr].depth * rho_H2O;
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

	double w_end = rain_melt + wcontent();
	if (fabs(w_begin - w_end) > 1.0e-10)
		dprintf("Year %d Day %d ERROR: Snow snow_or_rainfall diff %g\n", date.year, date.day, w_begin - w_end);
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

	// Water (initiate water mass at half of water holding capacity plus wilting point)
	wmass = ((*soilcolumn).soiltype.wp_min + (*soilcolumn).soiltype.awc_min * 0.5) * depth;
	imass = 0.0;
	wmass_sc = 0.0;
	por = 0.0;
	Dw = 0.0;
	Kw = 0.0;
	b = 0.0;
	Kw_max = 0.0;
	a = 0.0;
	imp_ice = 0.0;
	perc_mass = 0.0;
	runoff_drain = 0.0;
	runoff_baseflow = 0.0;
	psi = 0.0;
	mwcontabove = 0.0;
	awcont = 0.0;
	wsats = 0.0;
	awc = 0.0;
	wp = 0.0;
	perc_base = 0.0;
	Psi_max = 0.0;

	// Organic exponential in Campbell (1974) parameterization for soil 
	// water retention curves from CLM5.0 documentation 
	b_org = min(2.7 + 9.3 * depth_profile / sapric_depth, 12.0);
	// Organic soil matric potential at saturation 
	Psi_max_org = -min(10.3 - 0.2 * depth_profile / sapric_depth, 10.1);
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
	// Hydraulic conductivity at saturation for organic matter (m s-1) (CLM5.0 - Documentation; Hoag and Price 1995)
	Kw_max_org = max(0.28 - 0.2799 * depth_profile / sapric_depth * M_PER_MM, (*soilcolumn).soiltype.Kw_max_min);

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

	// N pools
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	NO2_mass = 0.0;
	NO_mass = 0.0;
	N2O_mass = 0.0;
	N2_mass = 0.0;
	nmass_sat = 0.0;
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

	// Daily heterotrophic respiration
	respsum = 0.0;
	arespsum = 0.0;
	acdec[0] = 0.0;
	acdec[1] = 0.0;
	acinc[0] = 0.0;
	acinc[1] = 0.0;

	// Soil labile carbon availability
	labile_carbon = 0.0;

	organic_frac = 0.0;
	mineral_frac = 0.0;
	peat_frac = 0.0;
	soc_frac = 0.0;
	bulk_density = 0.0;

	// Monthly variables
	for (int mth = 0; mth < 12; mth++) {
		mwmass[mth] = 0.0;
		mimass[mth] = 0.0;
		mwcont[mth] = 0.0;
		mKw[mth] = 0.0;
		mDw[mth] = 0.0;
		mdrain[mth] = 0.0;
		mnetmin[mth] = 0.0;
		mtemp[mth] = 0.0;
		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
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
	case WCONTMAX:
		prop = wcont_max();
		break;
	case MKW:
		prop = mKw[m];
		break;
	case MDW:
		prop = mDw[m];
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
	case CCONT:
		prop = ccont();
		break;
	case NCONT:
		prop = ncont();
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
	case ANIMMOB:
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
	default:
		// Should never happen
		fail("Unsupported soil property type (%d)\n", SPT);
	}
	return prop;
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
	if (!must_set && !((date.day == 0 && (*soilcolumn).patch.stand.first_year == date.year) || (iforganicsoilproperties == 2 && !acrotelm && !catotelm))) // DW_COMMENT restart years needs to update it then as well!
		return;

	// Get soil type
	Soiltype soiltype = (*soilcolumn).soiltype;

	// Mineral particle density set to 2.65 g cm–3 (Bonan 2019)
	double min_density = Mindensity * G_PER_KG * M3_PER_CM3;
	// Soil organic carbon density, standard bulk density of peat (g cm-3; Farouki 1981)
	double org_density = Orgdensity * G_PER_KG * M3_PER_CM3;

	// Analytic temp calc parameters
	thermdiff_0 = soiltype.thermdiff_0;
	thermdiff_15 = soiltype.thermdiff_15;
	thermdiff_100 = soiltype.thermdiff_100;

	// Heat capacity of mineral soil [J m-3 K-1; Bonan 2019]
	Cp_min = soiltype.Cp_min;

	// Acrotelm
	if (acrotelm) {
		wsats = acrotelm_por * depth;
		// Field capacity is set at 75% of saturation (Chaudhary 2017)
		awc = (acrotelm_por * 0.75 - peat_wp) * depth;
		wp = peat_wp * depth;
		por = acrotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - acrotelm_por;
		bulk_density = peat_frac * org_density;
		
		Kw_max = Kw_max_org;
		b = b_org;
		Psi_max = Psi_max_org;
		a = 0.0;

		// Initially saturated
		double wadded = wsats - (wmass + imass);
		wmass += wadded;
		(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
	}
	// Catotelm
	else if (catotelm) {
		wsats = catotelm_por * depth;
		// Field capacity is set at 75% of saturation (Chaudhary 2017)
		awc = (catotelm_por * 0.75 - peat_wp) * depth;
		wp = peat_wp * depth;
		por = catotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - catotelm_por;
		bulk_density = peat_frac * org_density;

		Kw_max = Kw_max_org;
		b = b_org;
		Psi_max = Psi_max_org;
		a = 0.0;

		// Initially saturated
		double wadded = wsats - (wmass + imass);
		wmass += wadded;
		(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
	}
	// Mineral soil and no influence of soil organic C
	else if (!iforganicsoilproperties) {

		wsats = soiltype.wsats_min * depth;
		awc = soiltype.awc_min * depth;
		wp = soiltype.wp_min * depth;
		por = soiltype.porosity * depth;
		organic_frac = soiltype.organic_frac;
		mineral_frac = soiltype.mineral_frac;
		peat_frac = 0.0;
		bulk_density = mineral_frac * min_density + organic_frac * org_density;

		Kw_max = soiltype.Kw_max_min;
		b = soiltype.b_min;
		Psi_max = soiltype.Psi_max_min;
		a = soiltype.a;

		// Two layer hydrology parameter
		perc_base = soiltype.perc_base;
	}
	// Soil organic content influences soil properties
	// Hudson (1994) showed that for each 1-percent increase in soil organic matter, the available water holding 
	// capacity in the soil increased by 3.7 percent. Soil water is held by adhesive and cohesive forces within 
	// the soil and an increase in pore space will lead to an increase in water holding capacity of the soil.
	else {

		// added moved water and ice due to shrinking wsats
		wmass += wmass_move;
		imass += imass_move;
		wmass_move = imass_move = 0.0;

		// Soil layer C content (kg C m-2)
		double soil_c = iforganicsoilproperties == 1 ? soc_frac * soiltype.soilC : ccont(true);

		// Soil C density in this soil layer (kg C m-3) 
		double soil_density = soil_c / (depth * M_PER_MM);

		// Calculate organic fraction for soil layer
		// This is the fraction of the material (i.e. 1 - porosity) that is organic
		double org_frac = min(1.0, soil_density / Orgdensity);
		double min_frac = 1.0 - org_frac;

		double porosity = min_frac * soiltype.wsats_min + org_frac * wsats_org;

		wsats = porosity * depth;
		awc = (min_frac * soiltype.awc_min + org_frac * awc_org) * depth;
		wp = (min_frac * soiltype.wp_min + org_frac * wp_org) * depth;
		por = wsats;
		organic_frac = (1.0 - porosity) * org_frac;
		mineral_frac = (1.0 - porosity) * min_frac;
		peat_frac = 0.0;
		bulk_density = mineral_frac * min_density + organic_frac * org_density;

		// Move water downwards if water content exceeds wsats
		if (wsats < wmass + imass) {
			double wmass_frac = wmass / (wmass + imass);
			double move = wmass + imass - wsats;
			wmass_move = wmass_frac * move;
			imass_move = (1 - wmass_frac) * move;
			wmass -= wmass_move;
			imass -= imass_move;
		}

		// The saturated hydraulic conductivity
		// The soil organic matter fraction is further subdivided into “connected” and “unconnected” 
		// fractions using a result from percolation theory (Stauffer and Aharony 1994, Berkowitz and 
		// Balberg 1992). Assuming that the organic and mineral fractions are randomly distributed 
		// throughout a soil layer, percolation theory predicts that above a threshold value 
		// org_frac = f_threshold, connected flow pathways consisting of organic material only exist
		// and span the soil space. Flow through these pathways interacts only with organic material, 
		// and thus can be described by Kw_max_org (CLM5.0 documentation).
		double perc_frac = org_frac >= threshold_frac ? 1.1 * pow(org_frac - threshold_frac, beta_perc) * org_frac : 0.0; // DW_COMMENT
		double uncon_frac = 1.0 - perc_frac;
		double Kw_max_uncon = uncon_frac / ((1.0 - org_frac) / soiltype.Kw_max_min + (org_frac - uncon_frac) / Kw_max_org);
		Kw_max = uncon_frac * Kw_max_uncon + (1.0 - uncon_frac) * Kw_max_org;

		// Exponential in Campbell (1974) parameterization for soil water retention curves
		b = min_frac * soiltype.b_min + org_frac * b_org;

		// Saturated soil water suction (mm)
		Psi_max = min_frac * soiltype.Psi_max_min + org_frac * Psi_max_org;

		// Soil organic C microorganism collocation factor for decomposition (Yan 2018)
		soiltype.a = min_frac * soiltype.clay_frac < 0.016 ? 0.0 : (min_frac * soiltype.clay_frac > 0.37 ? 1.0 : 2.8 * min_frac * soiltype.clay_frac - 0.046);

		// Analytic temp calc parameters
		thermdiff_0 = min_frac * soiltype.thermdiff_0 + org_frac * soiltype.thermdiff_0_org;
		thermdiff_15 = min_frac * soiltype.thermdiff_15 + org_frac * soiltype.thermdiff_15_org;
		thermdiff_100 = min_frac * soiltype.thermdiff_100 + org_frac * soiltype.thermdiff_100_org;

		// Two layer hydrology parameter
		perc_base = min_frac * soiltype.perc_base + org_frac * soiltype.perc_base_org;
	}

	// If two layer soil scheme then wmass is never below wp
	if (iftwolayersoil) {
		wmass = max(wmass, wp);
	}
}


/// Check if all water should be frozen
/** Check if all water should be frozen
 *  If maximum soil temperature of last year without phase change is below zero, 
 *  then freeze all water as an assumption
 */
void Soillayer::freeze_lyr() {
	double max_temp = -999.9;
	for (int mth = 0; mth < 12; mth++) {
		max_temp = max(max_temp, mtemp[mth]);
	}

	// If maximum temperature during the year is below zero, then freeze all water
	if (max_temp < 0.0) {
		imass += wmass;
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

		if (temp < 0.0) {

			// Maximum super cooled water fraction
			const double max_sc_frac = 0.5;

			// Temperature close to zero, maximum water is super cooled
			if (temp > -0.1) {
				wmass_sc = max_sc_frac * wsats;
			}
			// Niu scheme
			else if (ifsupercooled == 1) {
				wmass_sc = wsats * min(max_sc_frac, pow(Lheat * temp / (gravity * (K2degC + temp) * Psi_max), -1.0 / b));
			}
			// Jin scheme
			else if (ifsupercooled == 2) {
				double w_m = 0.05 * (*soilcolumn).soiltype.ssa;
				double S = 0.37;
				double C_1 = 0.041;
				double C = (*soilcolumn).soiltype.ssa * S * C_1 * pow(K2degC - temp, 0.5);
				wmass_sc = wsats * min(max_sc_frac, (C * pow(-temp, -0.5) + w_m - 0.25 * C) / 100.0);
			}

			Fwmass_sc = wmass_sc / depth;

			temp_old = min(0.0, temp_old);
			
			if (Fwmass > Fwmass_sc && temp < temp_old) {
				ifFreeze = true;
			} 
			else if (Fwmass < Fwmass_sc && temp > temp_old){
				ifThaw = true;
			}
		}
		else if (temp > 0.0 && Fimass > 0.0) {
			ifThaw = true;
			temp_old = 0.0;
			wmass_sc = wsats;
			Fwmass_sc = wmass_sc / depth;
		}
	}
	else {

		if (temp < 0.0 && Fwmass > 0.0) {
			ifFreeze = true;
			temp_old = 0.0;
		}
		else if (temp > 0.0 && Fimass > 0.0) {
			ifThaw = true;
			temp_old = 0.0;
		}
	}

	// Temperature difference from previous phase change
	delta_T = temp - temp_old;
	
	// Freezing temperatures and water present
	if (ifFreeze) {

		// Energy to freeze water [J m-3]
		energy = Cp_water * fabs(delta_T);

		// Fraction water that would freeze [m3 m-3]
		Ffreez = energy / Lheat;

		// Super cooled water that could freeze
		double Fwmass_freeze = Fwmass - Fwmass_sc;

		if (Fwmass_freeze > Ffreez) {
			Fimass += Ffreez;
			Fwmass -= Ffreez;

			temp = temp_old;

			// Energy used for phase change [J m-3]
			phase_change_energy = -energy;
		}
		else {
			// It is cold enough to freeze all the water and decrease the temperature below 0, 
			// since the rest of the water gets frozen now. However, it should never result in 
			// temperatures colder than the values coming from the Crank-Nicholson scheme.

			// Soil layer temperature rise resulting from freezing water.
			double temp_rise = Fwmass_freeze * Lheat / Cp_water;

			// Update soil temperature in this layer
			temp += temp_rise;
			temp = min(temp_old, temp);

			// All water is frozen
			Fimass += Fwmass_freeze;
			Fwmass -= Fwmass_freeze;

			// Energy used for phase change [J m-3]
			phase_change_energy = -Fwmass_freeze * Lheat;
		}
	} // End of Freezing

	// Thawing temperatures and ice present
	else if (ifThaw) {

		if (temp > 0.0) {

			// Energy from full soil layer to freeze water [J m-3]
			energy = Cp_water * fabs(delta_T);

			// Fraction ice that would thaw [m3 m-3]
			Fthaw = energy / Lheat;

			if (Fimass > Fthaw) {

				// There is some ice left, even after melt
				Fwmass += Fthaw;
				Fimass -= Fthaw;
				temp = temp_old;

				// Energy used for phase change [J m-3]
				phase_change_energy = -energy;
			}
			else {
				// There is more energy available than required for
				// melting the rest of the ice. Therefore the
				// left-over energy will increase the temperature.

				// Soil layer temperature reduction resulting from thawing ice.
				double temp_reduction = Fimass * Lheat / Cp_water;

				// Update soil temperature in this layer
				temp -= temp_reduction;
				temp = max(temp_old, temp);

				// All water is unfrozen
				Fwmass += Fimass;
				Fimass = 0.0;

				// Energy used for phase change [J m-3]
				phase_change_energy = -Fimass * Lheat;
			}
		}
		else if (Fwmass < Fwmass_sc) {

			// Energy from full soil layer to thaw water [J m-3]
			energy = Cp_water * fabs(delta_T);

			// Fraction ice that would thaw [m3 m-3]
			Fthaw = energy / Lheat;

			// Missing super cooled water that should be in liquid form
			double Fwmass_thaw = min(Fimass, Fwmass_sc - Fwmass);

			if (Fwmass_thaw > Fthaw) {

				// There is some ice left, even after melt
				Fwmass += Fthaw;
				Fimass -= Fthaw;
				temp = temp_old;

				// Energy used for phase change [J m-3]
				phase_change_energy = -energy;
			}
			else {
				// There is more energy available than required for
				// melting the super cooled water level. The enery used for 
				// melting the ice will decrease the water temperarture

				// Soil layer temperature reduction resulting from thawing ice.
				double temp_reduction = Fwmass_thaw * Lheat / Cp_water;

				// Update soil temperature in this layer
				temp -= temp_reduction;
				temp = max(temp_old, temp);

				// Either all ice is unfrozen or super cooled fraction is meet
				Fwmass += Fwmass_thaw;
				Fimass -= Fwmass_thaw;

				// Energy used for phase change [J m-3]
				phase_change_energy = -Fwmass_thaw * Lheat;
			}
		}
	} // End of thawing

	// Water and ice mass
	wmass = Fwmass * depth;
	imass = Fimass * depth;

	// Energy used for phase change [J]
	phase_change_energy *= depth * M_PER_MM;

	// Set new old temperature
	temp_old = temp;

	// If phase change has happend then update soil layers thermal properties 
	if (ifThaw || ifFreeze) {
		update_thermal_properties();
	}

	if (temp > 0.0 && imass > 0.0)
		dprintf("ERROR: Phase Change Year %d Day %d Depth %g Temp %g imass %g\n", date.year, date.day, depth_profile, temp, imass);

	// No negative wmass
	if (wmass < -EPS)
		fail("ERROR: Negative wmass at end of phase_change! Year %d Day %d depth_profile %g wmass %g\n", date.year, date.day, depth_profile, wmass);

	// No negative imass
	if (imass < -EPS)
		fail("ERROR: Negative imass at end of phase_change! Year %d Day %d depth_profile %g imass %g\n", date.year, date.day, depth_profile, imass);

	double wtotal_final = wmass + imass;

	// water balance checks:
	if ((fabs(wtotal_initial - wtotal_final)) > EPS) {
		fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
	}
}

/// Calculate phase change of water and ice for this soil layer
/** Calculate phase change of water and ice for this soil layer
 */
void Soillayer::phase_change_new() {

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

		if (temp < 0.0) {

			// Maximum super cooled water fraction
			const double max_sc_frac = 0.9;

			// Temperature close to zero, maximum water is super cooled
			if (temp > -0.02) {
				wmass_sc = max_sc_frac * wsats;
			}
			// Niu scheme
			else if (ifsupercooled == 1) {
				wmass_sc = wsats * min(max_sc_frac, pow(Lheat * temp / (gravity * (K2degC + temp) * Psi_max), -1.0 / b));
			}
			// Jin scheme
			else if (ifsupercooled == 2) {
				double w_m = 0.05 * (*soilcolumn).soiltype.ssa;
				double S = 0.37;
				double C_1 = 0.041;
				double C = (*soilcolumn).soiltype.ssa * S * C_1 * pow(K2degC - temp, 0.5);
				wmass_sc = wsats * min(max_sc_frac, (C * pow(-temp, -0.5) + w_m - 0.25 * C) / 100.0);
			}

			Fwmass_sc = wmass_sc / depth;

			temp_old = min(0.0, temp_old);

			if (Fwmass > Fwmass_sc && temp < temp_old) {
				ifFreeze = true;
			}
			else if (Fwmass < Fwmass_sc && temp > temp_old) {
				ifThaw = true;
			}
		}
		else if (temp > 0.0 && Fimass > 0.0) {
			ifThaw = true;
			temp_old = 0.0;
			wmass_sc = wsats;
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

	// No phase change, return zero
	if (!ifThaw && !ifFreeze)
		return;

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
		double Fwmass_thaw = min(Fimass, Fwmass_sc - Fwmass);

		// Can't thaw all
		if (Fwmass_thaw > Fthaw) {

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
			wmass = (Fwmass + Fwmass_thaw) * depth;
			imass = (Fimass - Fwmass_thaw) * depth;

			// Update thermal properties
			update_thermal_properties();

			// Energy used for phase change and to keep layer at initial temperature [J m-3]
			double energy_loss = -Fwmass_thaw * ((Cp_ice - Cp_water) * (temp_old + K2degC) + Lheat);

			// Calculate remaining energy to lower the temperature [J m-3]
			double energy_remain = energy - energy_loss;

			// Update soil temperature in this layer
			temp = temp_old + energy_remain / Ch;

			// Energy used for phase change [J m-3]
			phase_change_energy = Fwmass_thaw * Lheat;
		}
	} // End of Thawing

	// Set new old temperature
	temp_old = temp;

	if (temp > 0.0 && imass > 0.0)
		dprintf("ERROR: Phase Change Year %d Day %d Depth %g Temp %g imass %g\n", date.year, date.day, depth_profile, temp, imass);

	// No negative wmass
	if (wmass < -EPS)
		fail("ERROR: Negative wmass at end of phase_change! Year %d Day %d depth_profile %g wmass %g\n", date.year, date.day, depth_profile, wmass);

	// No negative imass
	if (imass < -EPS)
		fail("ERROR: Negative imass at end of phase_change! Year %d Day %d depth_profile %g imass %g\n", date.year, date.day, depth_profile, imass);

	double wtotal_final = wmass + imass;

	// water balance checks:
	if ((fabs(wtotal_initial - wtotal_final)) > EPS) {
		fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);
	}

	// Final energy
	double energy_final = Ch * (temp + K2degC) - phase_change_energy;

	// water balance checks:
	if ((fabs(energy_initial - energy_final)) > energy_initial * EPS) {
		dprintf("Year %d Day %d Depth %g Energy balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", date.year, date.day, depth_profile, energy_initial, energy_final, energy_initial - energy_final);
	}
}


/// Water content of soil layer
/** Water content of soil layer
 *  wcont never reports over field capacity even if there is more water
 */
double Soillayer::wcont() {
	return min(1.0, max(0.0, wmass - wp) / awc);
}

/// Field capacity of soil layer
/** Field capacity of soil layer
 */
double Soillayer::fc() {
	return awc + wp;
}

/// Fraction of air in soil layer
/** Fraction of air in soil layer
 */
double Soillayer::fracair() {
	return max(0.0, (por - (max(0.0, wmass) + imass)) / depth);
}

/// Get the water filled pore space of the layer
/** Get the water filled pore space of the layer
 */
double Soillayer::wfps() {
	return max(0.0, min(1.0, (wmass + imass) / wsats));
}

/// Get the maximum wcont possible
/** Get the maximum wcont possible
 */
double Soillayer::wcont_max() {
	return (wsats - wp) / awc;
}

/// C content of soil layer
/** C content of soil layer
 * Could also include/exclude litter
 */
double Soillayer::ccont(bool exclude_litter) {
	double ccont = 0.0;

	ccont += cpool_fast;
	ccont += cpool_slow;

	for (int i = 0; i < NSOMPOOL; i++) {
		if ( !(exclude_litter && (i == SURFMETA || i == SURFSTRUCT || i == SURFFWD || i == SURFCWD || i == SOILMETA || i == SOILSTRUCT || i == SOILFWD || i == SOILCWD)))
			ccont += sompool[i].cmass;
	}

	return ccont;
}

/// N content of soil layer
/** N content of soil layer
 * Could also include/exclude litter
 */
double Soillayer::ncont(bool exclude_litter) {
	double ncont = NH4_mass + NO3_mass + NO2_mass + NO_mass + N2O_mass + N2_mass;

	for (int i = 0; i < NSOMPOOL; i++) {
		if (!(exclude_litter && (i == SURFMETA || i == SURFSTRUCT || i == SURFFWD || i == SURFCWD || i == SOILMETA || i == SOILSTRUCT || i == SOILFWD || i == SOILCWD)))
			ncont += sompool[i].nmass;
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

	// If N transformation is turned off, then all available N should be of the type NH4
	if (!ifntransform) {
		pref = NH4;
	}

	double residual = 0.0;
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		if (nmass_tot > 0.0) {
			residual = nmass_tot - nmass;
			if (residual >= 0.0 || negligible(residual, 1.0e-10)) {
				nmass_subtract(nmass * NH4_mass / nmass_tot, NH4);
				nmass_subtract(nmass * NO3_mass / nmass_tot, NO3);
			}
			else {
				fail("ERROR: Year %d Day %d Tried to subtract more N (%f) than available (%f) Diff (%f)\n", date.year, date.day, nmass, nmass_tot, residual);
			}
		}
	}
	else if (pref == NH4) {
		if (nmass <= NH4_mass || negligible(nmass - NH4_mass, 1.0e-10)) {
			NH4_mass = max(0.0, NH4_mass - nmass);
		}
		else {
			double r = nmass - NH4_mass;
			fail("ERROR: Year %d Day %d Tried to subtract more NH4 (%f) than available (%f) Diff (%f)\n", date.year, date.day, nmass, NH4_mass, r);
		}
	}
	else if (pref == NO3) {
		if (nmass <= NO3_mass || negligible(nmass - NO3_mass, 1.0e-10)) {
			NO3_mass = max(0.0, NO3_mass - nmass);
		}
		else {
			double r = nmass - NO3_mass;
			fail("ERROR: Year %d Day %d Tried to subtract more NO3 (%f) than available (%f) Diff (%f)\n", date.year, date.day, nmass, NO3_mass, r);
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

/// Update soil layer temperature due to mass movement
/** Update soil layer temperature due to mass movement
 *  assuming it replaces air
 */
void Soillayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible((depth * Ch - add_mass * Cp_air) + add_mass * add_Ch))
		dprintf("ERROR: Year %d Day %d Updating soil temperature\n", date.year, date.day);

	double old_temp = temp;

	temp = (temp * (depth * Ch - add_mass * Cp_air) + add_temp * add_mass * add_Ch) / ((depth * Ch - add_mass * Cp_air) + add_mass * add_Ch);

	// If temperature becomes positive (negative) with ice (water) present
	// then perfomr a phase_change for this layer
	if (date.year >= FIRST_FREEZE_YEAR && ((temp > 0.0 && imass > 0.0) || (temp < 0.0 && wmass > 0.0))) {
		phase_change();
	}
}

/// 
/**
 */
double hydro_pow(double theta, double k) {

	double z = k * log_hydro_lut(max(0.001, theta));

	return exp_hydro_lut(z);
}

/// Updates soil water matric potential
/** Campbell (1974) soil water retention function,
 * which predicts the soil water matric potential
 */
void Soillayer::update_matric_potential() { // DW_COMMENT not needed yet
	double theta = max(0.0, wmass) / wsats;

	// psi = Psi_max * pow(theta, -b);
	psi = Psi_max * hydro_pow(theta, -b);
}

/// Updates hydraulic conductivity and diffusivity
/** Updates hydraulic conductivity and diffusivity
 * Relate the hydraulic conductivity and diffusivity 
 * to the soil water content following Campbell	(1974) 
 * and Clapp and Hornberger (1978)
 */
void Soillayer::update_hydraulic_properties() {

	double theta = max(0.0, wmass) / wsats;

	// Hydraulic conductivity (Campbell 1974; m s-1)
	if (wmass > wsats) {
		Kw = Kw_max;
	}
	else {
		// Kw = Kw_max * pow(theta, 2.0 * b + 3.0);
		Kw = Kw_max * hydro_pow(theta, 2.0 * b + 3.0);
	}

	// Hydraulic diffusivity (Clapp and Hornberger 1978; m2 s-1)
	// Dw = b * Kw_max * (-Psi_max) / (wsats / depth) * pow(theta, b + 2.0);
	Dw = b * Kw_max * (-Psi_max) / (wsats / depth) * hydro_pow(theta, b + 2.0);

	// The water transport in frozen soil is limited in the case of a partially frozen soil, by considering the effective 
	// hydraulic conductivity and diffusivity to be a weighted average of the values for total soil water and the value 
	// at the permanent wilting point for frozen water (for convenience set to zero)
	if (imass > 0.0) {
		// Ice impedance factor (Swenson et al. 2012)
		imp_ice = pow(10.0, -6.0 * min(1.0, imass / wsats));

		Dw *= imp_ice;
		Kw *= imp_ice;
	}
	else {
		imp_ice = 1.0;
	}
}

/// Updates heat capacity and thermal conductivity
/** Updates heat capacity and thermal conductivity of the soil layer
 */
void Soillayer::update_thermal_properties(bool combine_snow, double snow_depth, double snow_imass, double snow_wmass) {

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
	// Total solids in soil
	double Stot;

	double combined_depth = depth;
	double ice = imass;
	double water = wmass;
	double peat = peat_frac * depth;
	double mineral = mineral_frac * depth;
	double organic = organic_frac * depth;

	// Soil layer air (mm)
	double air = max(0.0, por - (wmass + imass));

	// Add snow properties if to combine thin snow layer to top soil layer
	if (combine_snow) {
		combined_depth += snow_depth;
		if (snow_imass + snow_wmass > snow_depth) {
			ice += snow_imass * snow_depth / (snow_imass + snow_wmass);
			water += snow_wmass * snow_depth / (snow_imass + snow_wmass);
			// No air
		}
		else {
			air += snow_depth - (snow_imass + snow_wmass);
			ice += snow_imass;
			water += snow_wmass;
		}
	}

	// Update heat capacity (J m-3 K-1)
	Ch = (mineral * Cp_min + organic * Cp_org + ice * Cp_ice + water * Cp_water + peat * Cp_peat + air * Cp_air) / combined_depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	// Total solids in soil
	Stot = mineral + organic + ice + water + peat;

	if ((Stot + air) / combined_depth > 1.000001)
		fail("Year %d Day %d Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g\nCombined depth %g Snow depth %g Air %g Ice %g Air %g", 
			date.year, date.day, Stot / combined_depth, air / combined_depth, (Stot + air) / combined_depth, combined_depth, snow_depth, snow_depth-snow_imass-snow_wmass, snow_imass, snow_wmass);

	// Calculate log terms
	lKforg = (organic / Stot) * Korg_log;
	lKfpeat = (peat / Stot) * Kpeat_log;
	lKfmin = (mineral / Stot) * Kmin_log;
	lKfwater = (water / Stot) * Kwater_log;
	lKfice = (ice / Stot) * Kice_log;

	// Total log terms
	double lKtot = lKfwater + lKfice + lKforg + lKfmin + lKfpeat;

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh = (air * Kair + Stot * exp_thermal_lut(lKtot)) / combined_depth;
}

/// Serialize soil layer variables
/** Serialize soil layer variables
 */
void Soillayer::serialize(ArchiveStream& arch) {
	arch & depth
		& wmass
		& imass
		& temp;

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
		& cpool_slow
		& cpool_fast

		// Equilsom variables
		& solvesom
		& fNH4uptake_mean
		& fNO3uptake_mean
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

	// Padding layers
	pad.resize(PAD_LAYERS);
	soillayerdepths.resize(1);
	cmass_litter_input.resize(1);
	csom_litter_input.resize(1);
	wfluxes.resize(1);

	depth = 0.0;
	nsoillayers = 0;

	wtimesteps = 12;
	ttimesteps = 24;

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
	percolate = false;

	// Peatland hydrology variables
	thaw = 0.0;
	maxthaw = 0.0;
	awtp = 0.0;
	acro_co2 = 934.0;
	wtd = 0.0;
	wtp = 0.0;
	stand_wmass = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mthaw[mth] = 0.0;
		mwtp[mth] = 0.0;
		mald_mean[mth] = 0.0;
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
	//solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	//solvesomcent_endyr = (int)(SOLVESOMCENT_SPINEND * (nyear_spinup - freenyears) + freenyears);

	solvesomcent_beginyr = (int)0.5*freenyears; // DW_COMMENT
	solvesomcent_endyr = (int)freenyears+10;
	savesomcent = false;
}

/// Create soil layers
/** Create soil layers with specific soil depths and thickness
 */
void Soilcolumn::create_soillayers(std::vector<double> soil_layer_depths) {

	// Set number of soil layers
	nsoillayers = (int)soil_layer_depths.size();

	// Update size
	soillayerdepths.resize(nsoillayers);

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

		add_soillayer(soil_layer_depths[lyr], depth_profile, acrotelm, catotelm);
	}

	// Set size of litter input vector (surface as well)
	cmass_litter_input.resize(nsoillayers + 1, 0.0);

	// Set size of SOM litter input vector
	csom_litter_input.resize(NSOMPOOL, 0.0);

	// Set sizes of water flow calculation parameters
	wpara.set_sizes(nsoillayers);

	// Set sizes of temperature change calculation parameters
	tpara.set_sizes(nsoillayers + PAD_LAYERS);

	// Set size of water flow vector
	wfluxes.resize(nsoillayers + 1, 0.0);
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

/// Returns surface litter C content
/** Returns surface litter C content
 */
double Soilcolumn::get_surflitterC() {
	return soillayers[0].sompool[SURFMETA].cmass + soillayers[0].sompool[SURFSTRUCT].cmass +
		soillayers[0].sompool[SURFFWD].cmass + soillayers[0].sompool[SURFCWD].cmass +
		soillayers[0].sompool[SURFMICRO].cmass + soillayers[0].sompool[SURFHUMUS].cmass;
}

/// Returns a vector with soillayer values for a certain property
/** Returns a vector with soillayer values for a certain property
 */
std::vector<double> Soilcolumn::get(soilproptype SPT, int m) {

	// Create vector to return
	std::vector<double> prop(nsoillayers);

	// Loop over soil layers and extract property
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		switch (SPT) {
		case DEPTH:
			prop[lyr] = soillayers[lyr].depth;
			break;
		case TEMP:
			prop[lyr] = soillayers[lyr].temp;
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
		case WCONTMAX:
			prop[lyr] = soillayers[lyr].wcont_max();
			break;
		case MKW:
			prop[lyr] = soillayers[lyr].mKw[m];
			break;
		case MDW:
			prop[lyr] = soillayers[lyr].mDw[m];
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
		case CCONT:
			prop[lyr] = soillayers[lyr].ccont();
			break;
		case NCONT:
			prop[lyr] = soillayers[lyr].ncont();
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
		case ANIMMOB:
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
		default:
			// Should never happen
			fail("Unsupported soil property type (%d)\n", SPT);
		}
	}

	// Return vector
	return prop;
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
std::vector<double> Soilcolumn::get_nfrac(n_pref_type pref, bool wmass_limited) {

	std::vector<double> nfrac(nsoillayers);

	double ntot = nmass_avail(pref, wmass_limited);
	double nfracsum = 0.0;

	if (ntot > 1.0e-12) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			nfrac[lyr] = soillayers[lyr].nmass_avail(pref, wmass_limited) / ntot;
			nfracsum += nfrac[lyr];
		}
		nfrac[0] += 1.0 - min(1.0, nfracsum);
	}
	else {
		nfrac[0] = 1.0;
	}

	return nfrac;
}

/// Returns a vector with soillayer C root fractions
/** Returns a vector with soillayer C root fractions
 */
std::vector<double> Soilcolumn::get_rootfrac() { // DW_COMMENT needed?
	std::vector<double> cmass_root(nsoillayers);
	double cmass_root_tot = 0.0;

	Vegetation& vegetation = patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			cmass_root[lyr] += indiv.cmass_root * patch.pft[indiv.pft.id].rootdist[lyr];
		}
		cmass_root_tot += indiv.cmass_root;

		vegetation.nextobj();

	} // indiviual loop

	std::vector<double> rootfrac(nsoillayers);

	if (!negligible(cmass_root_tot)) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			rootfrac[lyr] += cmass_root[lyr] / cmass_root_tot;
		}
	}
	else {
		rootfrac[0] = 1.0;
	}

	return rootfrac;
}

/// Returns a vector with soillayer C root mass
/** Returns a vector with soillayer C root mass
 */
std::vector<double> Soilcolumn::get_rootmass() { // DW_COMMENT needed
	std::vector<double> cmass_root(nsoillayers);

	Vegetation& vegetation = patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			cmass_root[lyr] += indiv.cmass_root * patch.pft[indiv.pft.id].rootdist[lyr];
		}

		vegetation.nextobj();

	} // indiviual loop

	return cmass_root;
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
 *  This should be done more advanced depending on root distribution and competition in each soil layer // DW_COMMENT 
 */
void Soilcolumn::nmass_avail_uptake(double uptake, n_pref_type pref) {
	double tot_nmass_avail = nmass_avail(pref);

	if (!negligible(tot_nmass_avail)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].nmass_subtract(uptake * soillayers[lyr].nmass_avail(pref) / tot_nmass_avail, pref);
		}
	}
}

/// Calculate water content for this month
/** Calculate water content for this month
 */
void Soilcolumn::calc_mwcont() {
	double acc_wcont = 0.0;
	double acc_awc = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Monthly wmass
		soillayers[lyr].mwmass[date.month] += soillayers[lyr].wmass / date.ndaymonth[date.month];

		// Monthly imass
		soillayers[lyr].mimass[date.month] += soillayers[lyr].imass / date.ndaymonth[date.month];

		// Monthly wcont
		soillayers[lyr].mwcont[date.month] += soillayers[lyr].wcont() / date.ndaymonth[date.month];

		// Monthly moisture conductivity
		soillayers[lyr].mKw[date.month] += soillayers[lyr].Kw / date.ndaymonth[date.month];

		// Monthly moisture diffusivity
		soillayers[lyr].mDw[date.month] += soillayers[lyr].Dw / date.ndaymonth[date.month];

		// Monthly lateral drainage
		soillayers[lyr].mdrain[date.month] += soillayers[lyr].runoff_drain;


		if (date.islastday) {

			acc_wcont += lyr == 0 ? soillayers[lyr].mwcont[date.month] * soillayers[lyr].awc / 2.0 :
				(soillayers[lyr - 1].mwcont[date.month] * soillayers[lyr - 1].awc + soillayers[lyr].mwcont[date.month] * soillayers[lyr].awc) / 2.0;

			acc_awc += lyr == 0 ? soillayers[lyr].awc / 2.0 : (soillayers[lyr - 1].awc + soillayers[lyr].awc) / 2.0;
			soillayers[lyr].mwcontabove = acc_wcont / acc_awc;
		}
	}
}

/// Calculate this months average soil layer temperature
/** Calculate this months average soil layer temperature
 */
void Soilcolumn::calc_mtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].mtemp[date.month] += soillayers[lyr].temp / date.ndaymonth[date.month];
	}

	for (int plyr = 0; plyr < PAD_LAYERS; plyr++) {
		pad[plyr].mtemp[date.month] += pad[plyr].temp / date.ndaymonth[date.month];
	}
}
 
/// Determines static soil layer organic fractions 
/** Determines static soil layer organic fractions 
 */
void Soilcolumn::set_static_socfrac() {

	// Only enter first year except if dynamic soil organic C influence soil properties
	if (!(date.day == 0 && patch.stand.first_year == date.year && iforganicsoilproperties == 1))
		return;

	// Use a beta value to reproduce the global vertical SOC profile shown in Fig 4 of Jobbagy & Jackson (2010)
	const double beta_SOC = 0.976;
	// This value gives SOC fractions as a fraction of the total in the top 100 cm in this Table:
	// Depth range	:	Calculated fraction	:	Jobbagy & Jackson (2000), Fig 4
	// 00 - 20cm	:	0.42				:	0.42
	// 20 - 40cm	:	0.26				:	0.22
	// 40 - 60cm	:	0.16				:	0.16
	// 60 - 80cm	:	0.10				:	0.11
	// 80 - 100cm	:	0.06				:	0.09

	double depth = soillayers[0].depth * CM_PER_MM;

	// Init first layer
	soillayers[0].soc_frac = 1.0 - pow(beta_SOC, depth); 
	double soc_frac_total = soillayers[0].soc_frac;

	for (int lyr = 1; lyr < nsoillayers; lyr++) {
		depth += soillayers[lyr].depth * CM_PER_MM;
		soillayers[lyr].soc_frac = 1.0 - pow(beta_SOC, depth) - (1.0 - pow(beta_SOC, depth - soillayers[lyr].depth * CM_PER_MM));
		soc_frac_total += soillayers[lyr].soc_frac;
	}

	// Put the residual SOC fraction in lowest soil layer
	soillayers[nsoillayers - 1].soc_frac += 1.0 - min(soc_frac_total, 1.0);
}

/// Average temperature down to a depth
/** Average temperature down to a depth
 */
double Soilcolumn::gtemp_rootdist(std::vector<double> rootdist) {

	// Calculate root weighted gtemp
	double avr_gtemp = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		avr_gtemp += soillayers[lyr].gtemp * rootdist[lyr];
	}
	return avr_gtemp;
}

/// Performs daily updates to soil layer temperatures using an analytic function
/**
 * Call each simulation day (if iftwolayersoil true) following update of daily air temperature
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

			if (date.islastday) {

				// Linearly interpolate soil thermal diffusivity given mean
				// soil water content from this layer and above

				// Only top layer influencing soil thermal diffusivity
				if (lyr == 0) {
					if (soillayers[lyr].mwcontabove < 0.15)
						current_k = ((soillayers[lyr].thermdiff_15 - soillayers[lyr].thermdiff_0) / 0.15 *
							soillayers[lyr].mwcontabove + soillayers[lyr].thermdiff_0) * DIFFUS_CONV;
					else
						current_k = ((soillayers[lyr].thermdiff_100 - soillayers[lyr].thermdiff_15) / 0.85 *
						(soillayers[lyr].mwcontabove - 0.15) + soillayers[lyr].thermdiff_15) * DIFFUS_CONV;

					k = current_k;
				}
				else {
					// Current and all above layers influence soil thermal diffusivity

					// Above layer
					if (soillayers[lyr - 1].mwcontabove < 0.15)
						above_k = ((soillayers[lyr - 1].thermdiff_15 - soillayers[lyr - 1].thermdiff_0) / 0.15 *
							soillayers[lyr - 1].mwcontabove + soillayers[lyr - 1].thermdiff_0) * DIFFUS_CONV;
					else
						above_k = ((soillayers[lyr - 1].thermdiff_100 - soillayers[lyr - 1].thermdiff_15) / 0.85 *
						(soillayers[lyr - 1].mwcontabove - 0.15) + soillayers[lyr - 1].thermdiff_15) * DIFFUS_CONV;

					// Current layer
					if (soillayers[lyr].mwcontabove < 0.15)
						current_k = ((soillayers[lyr].thermdiff_15 - soillayers[lyr].thermdiff_0) / 0.15 *
							soillayers[lyr].mwcontabove + soillayers[lyr].thermdiff_0) * DIFFUS_CONV;
					else
						current_k = ((soillayers[lyr].thermdiff_100 - soillayers[lyr].thermdiff_15) / 0.85 *
						(soillayers[lyr].mwcontabove - 0.15) + soillayers[lyr].thermdiff_15) * DIFFUS_CONV;

					k = (k * (soillayers[lyr].depth_profile - (soillayers[lyr - 1].depth + soillayers[lyr].depth) / 2.0) +
						above_k * soillayers[lyr - 1].depth / 2.0 + current_k * soillayers[lyr].depth / 2.0) / soillayers[lyr].depth_profile;
				}

				// Calculate parameters alag and exp(-alag) from Eqn 2
				soillayers[lyr].alag = soillayers[lyr].depth_profile * M_PER_MM / sqrt(k / HALF_OMEGA); // from Eqn 1
				soillayers[lyr].exp_alag = exp(-soillayers[lyr].alag);
			}

			// Every day, calculate linear model for trend in daily air
			// temperatures for the last 31 days: temp_day = a + b * day

			double buffer[31];
			climate.dtemp_31.to_array(buffer);
			regress(day, buffer, 31, a, b);

			// Calculate soil temperature lag
			temp_lag = a + b * (30.0 - soillayers[lyr].alag * LAG_CONV);

			// Calculate analytic solution for soil temperature
			soillayers[lyr].temp = climate.atemp_mean + soillayers[lyr].exp_alag * (temp_lag - climate.atemp_mean);
		}
	}
}

/// Returns the value of some vector down to a depth with a scalar
/** Returns the value of some vector down to a depth with a scalar
 */
double Soilcolumn::to_depth(double depth, soilproptype SPT_var, soilproptype SPT_scal, int mth) {

	// Calculate variable values for layers down to depth
	double layer_depth, scalar;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		scalar = soillayers[lyr].get(SPT_scal);
		acc_var += soillayers[lyr].get(SPT_var, mth) * scalar * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		acc_scalar += scalar * layer_depth / soillayers[lyr].depth;
		lyr++;
	}
	return !negligible(acc_scalar) ? acc_var / acc_scalar : acc_var;
}

/// Returns the value of some vector down to a depth
/** Returns the value of some vector down to a depth
 */
double Soilcolumn::to_depth(double depth, soilproptype SPT_var) {

	// Calculate variable values for layers down to depth
	double layer_depth;
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_var += soillayers[lyr].get(SPT_var) * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_var;
}

/// Add water to a certain soil layer and return the water that cannot be added due to limited capacity
/** Add water to a certain soil layer and return the water that cannot be added due to limited capacity
 */
double Soilcolumn::add_wmass(int layer, double extrawmass) {

	Soillayer& soillayer = soillayers[layer];

	// Add water
	soillayer.wmass += extrawmass;

	double overflow = 0.0;
	
	// Check if water exceed saturation
	if (soillayer.wmass + soillayer.imass > soillayer.wsats) {
		overflow = soillayer.wmass + soillayer.imass - soillayer.wsats;
		soillayer.wmass -= overflow;
	}
	else {
		overflow = 0.0;
	}

	return overflow;
}

/// Update soil matrix potential for the soil column
/** Update soil matrix potential for the soil column
 */
void Soilcolumn::update_matric_potential() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_matric_potential();
	}
}

/// Update soil hydraulic conductivity and diffusivity for the soil column
/** Update soil hydraulic conductivity and diffusivity for the soil column
 */
void Soilcolumn::update_hydraulic_properties() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_hydraulic_properties();
	}
}

/// Update soil thermal capacity and conductivity for the soil column
/** Update soil thermal capacity and conductivity for the soil column
 */
void Soilcolumn::update_thermal_properties(bool combine_snow) {

	// Update heat capacity and conductivity for each snow layer
	if (!combine_snow)
		snow.update_thermal_properties();

	// Update heat capacity and conductivity for each soil layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		if (lyr == 0 && combine_snow) {
			soillayers[lyr].update_thermal_properties(combine_snow, snow.snowlayers[0].depth, snow.snowlayers[0].imass, snow.snowlayers[0].wmass);
		}
		else {
			soillayers[lyr].update_thermal_properties();
		}
	}

	// Calculate heat capacity and conductivity for each padding layer
	update_thermal_properties_for_pad();
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

	// Use bottom soil layer properties for the first 3 padding layers
	const int first_bedrock_padding_layer = 3;

	for (int lyr = 0; lyr < PAD_LAYERS; lyr++) {

		// First padding layers use bottom soil layers properties
		if (lyr < first_bedrock_padding_layer) {
			pad[lyr].Kh = soillayers[nsoillayers - 1].Kh;
			pad[lyr].Ch = soillayers[nsoillayers - 1].Ch;
		}
		else {
			// Bedrock
			pad[lyr].Kh = Kbedrock;
			pad[lyr].Ch = Cbedrock;
		}
	}
}

/// Build diagonals and right hand side vector
/** Build diagonals and right hand side vector
 */
void build_diags_rhs(double dt, double f_top, double f_bot, CNpara& para, int n_soil_layers) {

	// Upper boundary conditions
	double l_p12 = para.diff[0] + (para.diff[1] - para.diff[0]) / (para.z[1] - para.z[0]) * para.dz[0] / 2.0; // Interpolate lambda_+1/2
	double g_p12 = para.cond[0] + (para.cond[1] - para.cond[0]) / (para.z[1] - para.z[0]) * para.dz[0] / 2.0; // Interpolate gamma_+1/2
	double beta_p = dt / 2.0 * l_p12 / para.dz[0] / (para.z[1] - para.z[0]) / para.cap[0];

	para.lower[0] = 1.0e6; // Not used, initialize to known value
	para.main[0] = 1.0 + beta_p;
	para.upper[0] = -beta_p;

	para.rhs[0] = (1. - beta_p) * para.psi[0] + beta_p * para.psi[1] - dt / para.cap[0] / para.dz[0] * (g_p12 - f_top) + para.sink[0] * dt;

	// Middle terms
	double l_m12, g_m12, beta_m;
	for (int i = 1; i < n_soil_layers - 1; i++) {
		l_m12 = l_p12;
		l_p12 = para.diff[i] + (para.diff[i + 1] - para.diff[i]) / (para.z[i + 1] - para.z[i]) * para.dz[i] / 2.0;
		beta_m = beta_p * para.dz[i - 1] / para.dz[i] * para.cap[i - 1] / para.cap[i];
		beta_p = dt / 2.0 * l_p12 / para.dz[i] / (para.z[i + 1] - para.z[i]) / para.cap[i];

		g_m12 = g_p12;
		g_p12 = para.cond[i] + (para.cond[i + 1] - para.cond[i]) / (para.z[i + 1] - para.z[i]) * para.dz[i] / 2.;

		para.lower[i] = -beta_m;
		para.upper[i] = -beta_p;
		para.main[i] = 1.0 + beta_m + beta_p;
		para.rhs[i] = beta_m * para.psi[i - 1] + (1.0 - beta_m - beta_p) * para.psi[i] + beta_p * para.psi[i + 1]
			- dt / para.cap[i] / para.dz[i] * (g_p12 - g_m12) + para.sink[i] * dt;
	}

	beta_m = beta_p * para.dz[n_soil_layers - 2] / para.dz[n_soil_layers - 1] * para.cap[n_soil_layers - 2] / para.cap[n_soil_layers - 1];
	g_m12 = g_p12;

	para.upper[n_soil_layers - 1] = 1.0e6; // Dummy value
	para.main[n_soil_layers - 1] = 1.0 + beta_m;
	para.lower[n_soil_layers - 1] = -beta_m;

	para.rhs[n_soil_layers - 1] = beta_m * para.psi[n_soil_layers - 2] + (1.0 - beta_m) * para.psi[n_soil_layers - 1]
		- dt / para.cap[n_soil_layers - 1] / para.dz[n_soil_layers - 1] * (f_bot - g_m12)
		+ para.sink[n_soil_layers - 1] * dt;
}

/// Solve tridiagonal system
/** Solve tridiagonal system
 */
void solve_tridiagonal(CNpara& para) {

	int n = (int)para.main.size();
	assert(n == para.psi1.size());

	for (int i = 1; i < n; i++) {
		double w = para.lower[i] / para.main[i - 1];
		para.main[i] -= w * para.upper[i - 1];
		para.rhs[i] -= w * para.rhs[i - 1];
	}

	para.psi1[n - 1] = para.rhs[n - 1] / para.main[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		para.psi1[i] = (para.rhs[i] - para.upper[i] * para.psi1[i + 1]) / para.main[i];
	}
}


/// Solve tridiagonal system
/** Solve tridiagonal system
 */
void tridiag(CNpara& para) {

	// Tridiagonal system solver from Numerical Recipes.

	int n = (int)para.dz.size();
	std::vector<double> gam(n);
	double bet;

	bet = para.main[0];

	para.psi1[0] = para.rhs[0] / bet;

	for (int j = 1; j < n; j++) {
		gam[j] = para.upper[j - 1] / bet;
		bet = para.main[j] - para.lower[j] * gam[j];

		para.psi1[j] = (para.rhs[j] - para.lower[j] * para.psi1[j - 1]) / bet;
	}

	for (int j = (n - 2); j >= 0; j--) {
		para.psi1[j] -= gam[j + 1] * para.psi1[j + 1];
	}
}


/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
void cnstep_full(double surf_temp, double dt, CNpara& para) {

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
	//      NLAYERS      Total number of layers represented in temp, d and
	//                   dz arrays. Not all of the layers have to be active
	//                   at one time, but this parameter is used to constrain
	//                   the sizes of the depth-based arrays.
	//
	//      layer0       Index of the first active layer. The diffusion
	//                   equation will be solved for all layers in the range
	//                   layer0 - NLAYERS inclusive, and the surface temperature
	//                   boundary condition will be applied to the top of layer
	//                   index layer0.
	//
	//      temp         The temperatures in each layer at the start of the
	//                   timestep.
	//
	//      d            The thermal diffusivities in each layer.
	//
	//      k            The thermal conductivities in each layer.
	//
	//      c            The heat capacities in each layer.
	//
	//      dz           The thickness of each layer.
	//
	//      surf_temp    The surface temperature boundary condition.
	//
	//      dt           The timestep.
	//
	//      PAD_LAYERS   Number of layers to use below computation domain to
	//                   deal with bottom boundary condition.
	//
	//      pad_dz       Padding layer thicknesses (derived from calc_padding).
	//
	//      pad_temp     Temperature in padding layers.
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
	//                   timestep. Only layers layer0 - NLAYERS will have
	//                   valid temperature values after a to cnstep. All
	//                   layers with index less than layer0 will have missing
	//                   value temperatures set (a value of -1.0E35).
	//
	//      pad_temp     Temperature in padding layers at the end of the
	//                   timestep.

	// As cnstep, but now we feed in ki and ci

	// Retrieve number of layers
	int nlayers = para.dz.size();

	// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	// BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	// lyr runs from 0 to active_layers - 1 
	for (int lyr = 0; lyr < nlayers; lyr++) {

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.

		// D+
		if (lyr == nlayers - 1) {
			dplus = 0.0;
		}
		else {
			dplus = ((para.diff[lyr] * para.dz[lyr] + para.diff[lyr + 1] * para.dz[lyr + 1]) / (para.dz[lyr] + para.dz[lyr + 1])) / para.cap[lyr];
		}

		// D-
		if (lyr == 0) {
			dminus = para.cond[lyr];
		}
		else {
			dminus = ((para.diff[lyr] * para.dz[lyr] + para.diff[lyr - 1] * para.dz[lyr - 1]) / (para.dz[lyr] + para.dz[lyr - 1])) / para.cap[lyr];
		}

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// and end cases.

		dzhere = para.dz[lyr];
		there = para.psi[lyr];

		// --- PLUS ---
		if (lyr < nlayers - 1) {
			// all soil layers apart from the bottom padding layer
			dzplus = para.dz[lyr + 1];
			tplus = para.psi[lyr + 1];
		}
		else {
			// bottom padding layer
			dzplus = para.dz[lyr];
			tplus = para.psi[lyr];
		}

		// --- MINUS ---
		if (lyr == 0) {
			// top soil layer
			dzminus = para.dz[lyr];
			tminus = para.psi[lyr];
		}
		else {
			// all layers except top soil layer
			dzminus = para.dz[lyr - 1];
			tminus = para.psi[lyr - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lyr == 0) {
			para.main[lyr] = 1.0;
		}
		else {
			para.main[lyr] = 1.0 + Cplus + Cp_minus;
		}

		// LEFT
		if (lyr == 0) {
			para.lower[lyr] = 0.0;
		}
		else {
			para.lower[lyr] = -Cp_minus;
		}

		// RIGHT
		if (lyr == 0 || lyr == nlayers - 1) {
			para.upper[lyr] = 0.0;
		}
		else {
			para.upper[lyr] = -Cplus;
		}

		// RHS
		// Calculate right hand side vector values.
		if (lyr == 0) {
			para.rhs[lyr] = surf_temp;
		}
		else {
			para.rhs[lyr] = (1.0 - Cplus - Cp_minus) * there +
				Cplus * tplus + Cp_minus * tminus;
		}
	}

	//   SOLVE TRIDIAGONAL SYSTEM
	tridiag(para);
}


/// Update soil layer temperatures due to water movement
/** Update soil layer temperatures due to water movement
 */
void Soilcolumn::wflow_temp() {

	const double EPS = 1.0e-12;

	// Update top soil layer for infiltarted rain or melt water
	if (wfluxes[0] > EPS) {
		soillayers[0].update_temp(rain_melt_temp, wfluxes[0], Cp_water);
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

/// Calculate soil column water flow
/** Calculate soil column water flow
 */
void Soilcolumn::calc_wflow(double winfilt, std::vector<double> aet) {

	bool negwmass = false;
	double drained_wmass;

	// Calculate matrix potential (not used yet)
	update_matric_potential();

	// Calculate hydraulic conductivity and diffusivity
	update_hydraulic_properties();

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		wpara.dz[lyr] = soillayers[lyr].depth * M_PER_MM;
		wpara.z[lyr] = soillayers[lyr].depth_profile * M_PER_MM;
		wpara.cap[lyr] = 1.0;
		wpara.diff[lyr] = soillayers[lyr].Dw;
		wpara.cond[lyr] = soillayers[lyr].Kw;
		wpara.sink[lyr] = -(aet[lyr] / soillayers[lyr].depth) / SECS_PER_DAY;
		wpara.psi[lyr] = soillayers[lyr].wmass / soillayers[lyr].depth;
	}

	// Boundary conditions for soil water
	// Upper boundary condition [m s-1]
	double f_top = winfilt / rho_H2O / SECS_PER_DAY;

	// Lower boundary condition
	double f_bot = ifbedrock ? 0.0 : wpara.cond[nsoillayers - 1];

	// Build diagonals and RHS vector
	build_diags_rhs(SECS_PER_DAY / wtimesteps, f_top, f_bot, wpara, nsoillayers);

	// Solve tridiagonal system
	solve_tridiagonal(wpara);

	// Set new wmass
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Water mass in soil layer
		soillayer.wmass = wpara.psi1[lyr] * soillayer.depth;

		// Calculate lateral drainage

		// wmass can have a small negative value (resolved later when mass balance can be kept)
		// if so, no lateral drainage is calculated
		if (soillayer.wmass > 0.0) {

			// Liquid water above saturation is drained laterally
			drained_wmass = min(soillayer.wmass, max(0.0, soillayer.wmass + soillayer.imass - soillayer.wsats));

			// Sum todays lateral drainage
			soillayer.runoff_drain += drained_wmass;

			// Update total water after drainage 
			// remove drained water from wmass. Water and ice will now maximum fill up soil layer to saturation
			soillayer.wmass -= drained_wmass;

			// If orography affects sub surface lateral drainage (Steenhuis 2019 eqn 14)
			if (iforography && (soillayer.wmass > soillayer.fc())) {

				// Excees liquid water above field capacity
				double wmass_excees = soillayer.wmass - soillayer.fc();

				// Slope dependant drainage of liquid water above field capacity
				drained_wmass = wmass_excees * min(1.0, soillayer.Kw * soiltype.slope_grade * MM_PER_M * SECS_PER_DAY / wtimesteps * soillayer.imp_ice);

				// Sum todays lateral drainage
				soillayer.runoff_drain += drained_wmass;

				// Update total water after drainage 
				// remove drained water from wmass. Water and ice will now maximum fill up soil layer to saturation
				soillayer.wmass -= drained_wmass;
			}
		}
	}
}

/// Resolves negative soil layer water mass
/** Resolves negative soil layer water mass
*/
void Soilcolumn::resolve_neg_wmass(std::vector<int> neglyr, double& evap, double& runoff_surf) {

	// Resolve negative wmass
	for (int lyr = 0; lyr < neglyr.size(); lyr++) {
		int lyr_neg = neglyr[lyr];

		// Start with checking fluxes from layer

		// If drained water then take negative wmass first from it
		if (soillayers[lyr_neg].runoff_drain > 0.0) {
			double wmass_move = min(soillayers[lyr_neg].runoff_drain, -soillayers[lyr_neg].wmass);
			soillayers[lyr_neg].runoff_drain -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
		}

		// Still negative, take from imass
		if (soillayers[lyr_neg].wmass < 0.0 && soillayers[lyr_neg].imass > 0.0) {
			double wmass_move = min(soillayers[lyr_neg].imass, -soillayers[lyr_neg].wmass);
			soillayers[lyr_neg].imass -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
		}

		// Soil layer below have taken wmass from negative layer
		if (wfluxes[lyr_neg + 1] > 0.0 && soillayers[lyr_neg].wmass < 0.0) {

			double max_move;
			if (lyr_neg < nsoillayers - 1) {
				max_move = max(0.0, soillayers[lyr_neg + 1].wmass);
			}
			else {
				max_move = max(0.0, soillayers[lyr_neg].runoff_baseflow);
			}

			double wmass_move = min(max_move, min(wfluxes[lyr_neg + 1], -soillayers[lyr_neg].wmass));
			wfluxes[lyr_neg + 1] -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;

			if (lyr_neg < nsoillayers - 1) {
				soillayers[lyr_neg + 1].wmass -= wmass_move;
			}
			else {
				soillayers[lyr_neg].runoff_baseflow -= wmass_move;
			}
		}

		// Soil layer above have taken wmass from negative layer
		if (wfluxes[lyr_neg] < 0.0 && soillayers[lyr_neg].wmass < 0.0) {

			double max_move;
			if (lyr_neg > 0) {
				max_move = max(0.0, soillayers[lyr_neg - 1].wmass);
			}
			else {
				max_move = max(0.0, evap);
			}

			double wmass_move = min(max_move, min(-wfluxes[lyr_neg], -soillayers[lyr_neg].wmass));
			wfluxes[lyr_neg] += wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
			if (lyr_neg > 0) {
				soillayers[lyr_neg - 1].wmass -= wmass_move;
			}
			else {
				evap -= wmass_move;
			}
		}

		// Very abnormal error in water movment. Happens on a few occations in a global simulation.
		// Extreme measures to be taken to keep soil layer wmass above zero 

		// Still negative, take from evaporation
		if (soillayers[lyr_neg].wmass < 0.0 && evap > 0.0) {
			double wmass_move = min(evap, -soillayers[lyr_neg].wmass);
			evap -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
		}

		// Still negative, take from surface runoff
		if (soillayers[lyr_neg].wmass < 0.0 && runoff_surf > 0.0) {
			double wmass_move = min(runoff_surf, -soillayers[lyr_neg].wmass);
			runoff_surf -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
		}

		// Still negative, take from wmass from layers below
		if (soillayers[lyr_neg].wmass < 0.0) {
			for (int lyr = lyr_neg + 1; lyr < nsoillayers; lyr++) {
				
				if (soillayers[lyr_neg].wmass < 0.0 && soillayers[lyr].wmass > 0.0) {
					double wmass_move = min(soillayers[lyr].wmass, -soillayers[lyr_neg].wmass);
					soillayers[lyr].wmass -= wmass_move;
					soillayers[lyr_neg].wmass += wmass_move;
				}
			}
		}

		// Still negative, take from wmass from layers above
		if (soillayers[lyr_neg].wmass < 0.0) {
			for (int lyr = lyr_neg - 1; lyr >= 0; lyr--) {

				if (soillayers[lyr_neg].wmass < 0.0 && soillayers[lyr].wmass > 0.0) {
					double wmass_move = min(soillayers[lyr].wmass, -soillayers[lyr_neg].wmass);
					soillayers[lyr].wmass -= wmass_move;
					soillayers[lyr_neg].wmass += wmass_move;
				}
			}
		}

		// Still negative, take from drainage from layers below
		if (soillayers[lyr_neg].wmass < 0.0) {
			for (int lyr = lyr_neg + 1; lyr < nsoillayers; lyr++) {

				if (soillayers[lyr_neg].wmass < 0.0 && soillayers[lyr].runoff_drain > 0.0) {
					double wmass_move = min(soillayers[lyr].runoff_drain, -soillayers[lyr_neg].wmass);
					soillayers[lyr].runoff_drain -= wmass_move;
					soillayers[lyr_neg].wmass += wmass_move;
				}
			}
		}

		// Still negative, take from drainage from layers above
		if (soillayers[lyr_neg].wmass < 0.0) {
			for (int lyr = lyr_neg - 1; lyr >= 0; lyr--) {

				if (soillayers[lyr_neg].wmass < 0.0 && soillayers[lyr].runoff_drain > 0.0) {
					double wmass_move = min(soillayers[lyr].runoff_drain, -soillayers[lyr_neg].wmass);
					soillayers[lyr].runoff_drain -= wmass_move;
					soillayers[lyr_neg].wmass += wmass_move;
				}
			}
		}

		// Still negative, take from baseflow
		if (soillayers[lyr_neg].wmass < 0.0 && soillayers[nsoillayers - 1].runoff_baseflow > 0.0) {
			double wmass_move = min(soillayers[nsoillayers - 1].runoff_baseflow, -soillayers[lyr_neg].wmass);
			soillayers[nsoillayers - 1].runoff_baseflow -= wmass_move;
			soillayers[lyr_neg].wmass += wmass_move;
		}

		if (soillayers[lyr_neg].wmass < -1.0e-12)
			dprintf("ERROR: Negative soil layer water! Year %d Day %d Lyr %d wmass %g\n", date.year, date.day, lyr_neg, soillayers[lyr_neg].wmass);
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

	// Layer counter.
	int layer;

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
	if (dz * PAD_LAYERS < PAD_DEPTH) {

		// Do Newton-Raphson solve to find geometric sequence factor.
		while (fabs(k_guess - k_last) > 1E-2) {
			k_last = k_guess;

			func = pow(k_guess, (PAD_LAYERS + 1)) - k_guess * (1.0 + PAD_DEPTH / dz) + PAD_DEPTH / dz;
			dfunc = (PAD_LAYERS + 1) * pow(k_guess, PAD_LAYERS) - (1.0 + PAD_DEPTH / dz);
			k_guess = k_guess - func / dfunc;
		}

		double tot_depth = 0.0;
		double lyr_depth;

		// Build geometric sequence of padding layer depths.
		for (layer = 0; layer < PAD_LAYERS; layer++) {
			lyr_depth = dz * pow(k_guess, layer + 1);
			pad[layer].depth = lyr_depth;
			tot_depth += lyr_depth;
		}

		if (!negligible(PAD_DEPTH - tot_depth)) {
			pad[PAD_LAYERS - 1].depth += PAD_DEPTH - tot_depth;
		}

		if (pad[PAD_LAYERS - 1].depth < -1.0e-10)
			fail("ERROR: Padding layer negative\n");
	}
	// Else each padding layer decrease linear to a minimum of 2m for the bottom padding layer
	else {
		double first_pad = min(17200.0, soillayers[nsoillayers - 1].depth);
		double diff = (5.0 * first_pad - PAD_DEPTH) / 10.0;
		int scalar = 0;

		for (layer = 0; layer < PAD_LAYERS; layer++) {
			pad[layer].depth = first_pad - diff * (double)scalar;
			scalar++;
		}
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
* Call each simulation day (if iftwolayersoil false).
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
	// Runoff (drainage) from lower soil layers (mm)
	double runoff_drain;
	// Evaporation from soil surface (mm)
	double evap;
	// Evapotranspiration (mm/day)
	double evapotranspiration;
	// Wetland runon set externally (mm)
	double runon = 0.0;
	// Water added to keep water mass balance (mm)
	double waddon = 0.0;
	// Total fpc for mosses. Not used at the moment
	double fpc_moss_total = 0.0;
	// Total water in acrotelm layers (mm)
	double wtot_acrotelm = to_depth(ACROTELM_DEPTH, WMASS);
	// Total ice in acrotelm layers (mm)
	double itot_acrotelm = to_depth(ACROTELM_DEPTH, IMASS);
	/// Minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// The soil water characteristics in the upper peat are considered
	/// to be linear in the suction interval 0 - 0.1 m
	/// Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22
	/// Unit here being mm, m in Granberg (1999)
	double az = (acrotelm_por - minvtot) / zmin;
	
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

	// Water extracted from catotelm
	waddon += aet_total - aet_acrotelm;

	
	// Evaporation

	// Evaporation from soil surface (Wania et al 2009)
	if (snow.depth() < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
		evap = fevap * climate.eet * PRIESTLEY_TAYLOR * (0.99 / (1.0 + exp(-1.0 * (-wtd + 98.7) / 22.6) + 0.02));
	}
	else {
		evap = 0.0;
	}
		

	// Evapotranspiration

	bool empirical_evapotranspiration = false;
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
	if (wtot_acrotelm > 0.0 && soillayers[0].imass / soillayers[0].wsats < 0.7)
		runoff_surf = min(ideal_runoff, wtot_acrotelm);
	else
		runoff_surf = 0.0;	

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
	double wtot_acrotelm_init = wtot_acrotelm;

	// Update acrotelm water mass
	wtot_acrotelm += runon - runoff_surf + stand_wmass;

	// If negative total water, then all initial water has gone to runoff and water is set to zero
	if (wtot_acrotelm < 0.0) {
		wtot_acrotelm = 0.0;
		runoff_surf = wtot_acrotelm_init;
		runon = 0.0;
		stand_wmass = 0.0;
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
			runoff_surf += wtp - maxh;
			wtp = maxh;
		}

		// WTD defined as < 0 above the surface
		wtd = -wtp;

		wtot_acrotelm = ACROTELM_DEPTH * acrotelm_por - itot_acrotelm;

		// Wetlands does NOT include standing water
		// we assume no standing water 
		stand_wmass = 0.0;
	}
	else {

		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * az));

		if (wtd > zmin)
			wtd = 3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * (acrotelm_por - minvtot));
		
		wtd = min(ACROTELM_DEPTH, wtd);

		wtp = -wtd;

		// Wetlands does NOT include standing water
		// we assume no standing water 
		stand_wmass = 0.0;
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
	std::vector<double> value(nsublayers_acro);
	std::vector<double> sub_water(nsublayers_acro);

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

			// Number of sublayers in acrotelm layer
			int sublayers_per_layer = (int)(soillayers[lyr].depth / Dz_sub);

			double layerwater = 0.0;
			for (int sl = 0; sl < sublayers_per_layer; sl++) {
				layerwater += sub_water[sublyr] * Dz_sub;
				sublyr++;
			}

			// Water should exceed saturation
			soillayers[lyr].wmass = min(soillayers[lyr].wsats - soillayers[lyr].imass, max(0.0, layerwater - soillayers[lyr].imass));
			wtot_acrotelm += soillayers[lyr].wmass;
		}
	}
	else {

		// WTD <= 0.0, at or above the surface - standing water
		wtot_acrotelm = 0.0;
		// Acrotelm saturated (subtract ice content)
		for (int lyr = 0; lyr < nacrotelm; lyr++) {
			soillayers[lyr].wmass = soillayers[lyr].wsats - soillayers[lyr].imass;
			wtot_acrotelm += soillayers[lyr].wmass;
		}
	}

	// Catotelm assumed to be always saturated
	for (int lyr = nacrotelm; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wmass = soillayers[lyr].wsats - soillayers[lyr].imass;
	}

	// Total water change during calculation of water content based on WTD
	waddon += wtot_acrotelm - wtot_acrotelm_init;

	// Remove missing water from runoff
	double wmass_moved = min(runoff_surf, waddon);
	runoff_surf -= wmass_moved;
	waddon -= wmass_moved;


	// Report water fluxes

	// Evapotranspiration data
	if (empirical_evapotranspiration) {
		patch.fluxes.report_flux(Fluxes::WEVAP, fevap * evapotranspiration);
		patch.fluxes.report_flux(Fluxes::WAET, (1.0 - fevap) * evapotranspiration);
	}
	else {
		patch.fluxes.report_flux(Fluxes::WEVAP, evap);
		patch.fluxes.report_flux(Fluxes::WAET, aet_total);
	}

	// Surface runoff
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);

	// Runon
	patch.fluxes.report_flux(Fluxes::WADDED, runon + waddon);


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
	// wetland_wmass_added_today = water added to soil to keep saturation (mm)

	// Minimum water considerred (equals to a 0.1 mm3)
	const double EPS = 1.0e-9;

	double wmass_begin = get_sum(WMASS) + rain_melt - patch.wetland_wmass_added_today;

	double evap = 0.0;

	double runoff_surf = 0.0;

	// AET for a particular layer and individual (mm/day)
	double aet;
	// total AET for each soil layer (mm/day)
	std::vector<double> aet_layer(nsoillayers);
	// total AET (mm/day)
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals and extract it from soil lalyers
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			double fwu = patch.pft[indiv.pft.id].fwuptake[lyr]; // DW_COMMENT not needed but crashed otherwise on Simba???????
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
		patch.wetland_wmass_added_today -= wbudget;
	}
	// Excess was is lost through surface runoff
	else {
		runoff_surf = wbudget;
	}

	// Save fluxes
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, 0.0);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_BASE, 0.0);
	patch.fluxes.report_flux(Fluxes::WEVAP, evap);
	patch.fluxes.report_flux(Fluxes::WAET, aet_total);
	patch.fluxes.report_flux(Fluxes::WADDED, patch.wetland_wmass_added_today);

	// Calculate annual wcont for growing season (temperature > 5C)
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Reset the awcont array on the first day of every year
		if (date.day == 0) {
			soillayers[lyr].awcont = 0.0;
			if (lyr == 0) patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			soillayers[lyr].awcont += soillayers[lyr].wcont();
			if (lyr == 0) patch.growingseasondays++;
		}

		// Do the averaging on the last day of every year
		if (date.islastday && date.islastmonth) {

			// In case it's never warm enough:
			if (!patch.growingseasondays) {
				soillayers[lyr].awcont = 1.0;
			}
			else {
				soillayers[lyr].awcont /= (double)patch.growingseasondays;
			}
		}
	}

	// Water balance check
	double wmass_end = get_sum(WMASS) + runoff_surf + evap + aet_total - patch.wetland_wmass_added_today;

	if (fabs(wmass_begin - wmass_end) > 1.0e-10)
		dprintf("Year %d Day %d WMASS ERROR: Diff %g Begin %g End %g\n", date.year, date.day, wmass_begin - wmass_end, wmass_begin, wmass_end);
}

/// Updates the water content for each soil layer
/**
* Call each simulation day (if iftwolayersoil true).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soilcolumn::hydrology_lpjf_twolayer(const Climate& climate, double fevap) {

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

	// Actual Evapotranspiration from soil column

	// AET for a particular layer and individual (mm)
	double aet;
	// total AET for each Gerten-type soil layer (mm)
	std::vector<double> aet_layer(nsoillayers);	
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

	// Evaporation from soil surface

	double evap = 0.0;

	// Evaporation only if snow equivalent water < 10mm
	if (snow.wcontent() < 10.0) {	
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}

	// Implement in- and outgoing fluxes to upper soil layer
	soillayers[0].wmass += rain_melt - aet_layer[0] - evap;

	// Two-layer hydrology always have water up to wilting point
	if (soillayers[0].wmass <= soillayers[0].wp) { 
		evap += soillayers[0].wmass - soillayers[0].wp;
		soillayers[0].wmass = soillayers[0].wp;
	}

	// Surface runoff if surface soil layer has water over field capacity
	double runoff_surf = 0.0;
	if (soillayers[0].wmass > soillayers[0].fc()) {
		runoff_surf = soillayers[0].wmass - soillayers[0].fc();
		soillayers[0].wmass = soillayers[0].fc();
	}

	// Update water content in evaporation layer for tomorrow
	wcont_evap += min(wcont_evap,
		(rain_melt - aet_layer[0] * K_AET_DEPTH * min(1.0, SOILDEPTH_EVAP / soillayers[0].depth) - evap)
		/ soillayers[0].awc);

	if (wcont_evap > soillayers[0].wcont()) {
		wcont_evap = soillayers[0].wcont();
	}

	// Percolation from evaporation layer
	double perc = 0.0;
	if (percolate) {
		perc = min(min(1.0, SOILDEPTH_EVAP / soillayers[0].depth) * soillayers[0].perc_base * pow(wcont_evap, soiltype.perc_exp),
			max_rain_melt);
	}
	wcont_evap -= perc / (min(1.0, SOILDEPTH_EVAP / soillayers[0].depth) * soillayers[0].awc);

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	for (int lyr = 1; lyr < nsoillayers; lyr++) {

		// Percolation

		// Allow only on days with rain or snowmelt
		if (percolate) {
			perc = min(soillayers[lyr - 1].perc_base * pow(soillayers[lyr - 1].wcont(), soiltype.perc_exp), max_rain_melt);
		}
		else {
			perc = 0.0;
		}

		// Percolation water from above layer (mm)
		soillayers[lyr - 1].perc_mass = min(perc, soillayers[lyr - 1].wmass - soillayers[lyr - 1].wp);

		// Remove and transfer water (mm)
		soillayers[lyr - 1].wmass -= soillayers[lyr - 1].perc_mass;
		soillayers[lyr].wmass += soillayers[lyr - 1].perc_mass;
		if (soillayers[lyr].wmass > soillayers[lyr].fc()) {
			soillayers[lyr].runoff_drain = soillayers[lyr].wmass - soillayers[lyr].fc();
			soillayers[lyr].wmass = soillayers[lyr].fc();
		}
		else {
			soillayers[lyr].runoff_drain = 0.0;
		}

		// Deduct AET from this soil layer
		soillayers[lyr].wmass -= min(soillayers[lyr].wmass - soillayers[lyr].wp, aet_layer[lyr]);
	}

	// Baseflow runoff (rain or snowmelt days only)
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC * soillayers[nsoillayers - 1].perc_base * 
			pow(soillayers[nsoillayers - 1].wcont(), soiltype.perc_exp);
		
		// Check so baseflow doesn't exceed inflow of water
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer
		soillayers[nsoillayers - 1].perc_mass = min(perc_baseflow, soillayers[nsoillayers - 1].wmass - soillayers[nsoillayers - 1].wp);
		soillayers[nsoillayers - 1].wmass -= soillayers[nsoillayers - 1].perc_mass;
		soillayers[nsoillayers - 1].runoff_baseflow = soillayers[nsoillayers - 1].perc_mass;
	}
	else {
		soillayers[nsoillayers - 1].runoff_baseflow = 0.0;
	}

	// Water balance check
	double wafter = get_sum(WMASS) + runoff_surf + get_sum(RUNOFFDRAIN) + get_sum(RUNOFFBASE) + aet_total + evap;
	if (fabs(wbefore-wafter) > 1.0e-10)
		dprintf("Year %d Day %d BUCKETWATER Before %g After %g Diff %g Rain_melt %g\n",date.year,date.day, wbefore, wafter, wbefore - wafter, rain_melt);

	// Reset rain_melt
	rain_melt = 0.0;

	// Save fluxes
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(RUNOFFDRAIN));
	patch.fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(RUNOFFBASE));
	patch.fluxes.report_flux(Fluxes::WEVAP, evap);
	patch.fluxes.report_flux(Fluxes::WAET, aet_total);

	// Update awcont_upper

	// Reset awcont_upper on the first day of every year
	if (date.day == 0) {
		awcont_upper = 0.0;
		patch.growingseasondays = 0;
	}

	// If it's warm enough for growth, update awcont_upper with this day's wcont
	if (climate.temp > 5.0) {
		awcont_upper += soillayers[0].wcont();
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

/// Updates the water content for each soil layer in non-peatland soils
/**
* Call each simulation day (if iftwolayersoil false).
* Soil class method used to update the water content for each soil layer
* given snow melt, rainfall, evapotranspiration from vegetation (AET) and
* percolation between layers. Runoff is also calculated.
*
* \param climate   The climate to use to update soil water
* \param fevap     fraction of modelled area (grid cell or patch) subject to
*                  evaporation from soil surface
*/
void Soilcolumn::hydrology_lpjf(const Climate& climate, double fevap) {

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

	// Maximum effective depth of surface runoff - 0.5m
	const double DEPTH_INFILT_MAX = 500.0;

	// Effective depth of surface runoff [mm]
	double depth_infilt = min(thaw, DEPTH_INFILT_MAX);

	double initial_wmass = snow.wcontent() + get_sum(WMASS) + rain_melt;

	// Evaposation (mm)
	double evap = 0.0;
	// Surface runoff (mm)
	double runoff_surf = 0.0;
	// AET for a particular layer and individual (mm/day)
	double aet;
	// total AET for each soil layer (mm/day)
	std::vector<double> aet_layer(nsoillayers);
	// total AET (mm/day)
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals and extract it from soil lalyers
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			double fwu = patch.pft[indiv.pft.id].fwuptake[lyr]; // DW_COMMENT not needed but crashed otherwise on Simba???????
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
			aet_layer[lyr] += aet;
			aet_total += aet;
		}
		vegetation.nextobj();
	}

	// Calculate static ice mass and saturation level for infiltation depth (top layer)
	double imass_infilt = to_depth(depth_infilt, IMASS);
	double wsats_infilt = to_depth(depth_infilt, WSATS);

	// Reset daily infiltation flux, evaporation and surface and drainage runoff 
	// and set initial water fraction to be able to calculate daily water fluxes
	double winfilt = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].runoff_drain = 0.0;
		wpara.psi0[lyr] = soillayers[lyr].wmass / soillayers[lyr].depth;
	}

	// Loop through sub daily timesteps to be able to solve 
	// the vertical flow of water with the Crank - Nicolson method
	for (int ts = 0; ts < wtimesteps; ts++) {

		// Evaporation

		// Evaporation is calculated on a sub daily time step as water content 
		// in top soil layer changes with each time step

		// Calculate wcont for top layer (evaporation layer)
		wcont_evap = soillayers[0].wcont();

		// Calculate total evaporation only if snow depth < 10mm
		double evap_sd = 0.0;
		if (snow.depth() < 10.0 && !negligible(fevap)) {
			evap_sd = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;

			// Remove small evaporations
			if (evap_sd < EPS)
				evap_sd = 0.0;

			evap += evap_sd / wtimesteps;
		}

		// Calculate water influx after rain_melt initially meets the evaporation demand, can be negative (mm/day)
		double winflux_sd = rain_melt - evap_sd;

		// Surface infiltration and runoff

		// Calculate changing water content for infiltration layers
		double wmass_infilt = to_depth(depth_infilt, WMASS);

		// Fraction pore spaces that are filled with water and ice
		double totwcont_infilt = (wmass_infilt + imass_infilt) / wsats_infilt;

		// Calculate runoff if there is any throughfall or snow melt 
		double runoff_surf_sd = 0.0;

		// Surface runoff depends on orography (Balsamo 2009, eqn 6 and 11)
		/*if (iforography) {
			if (winflux_sd > EPS) {
				runoff_surf_sd = max(0.0, winflux_sd - max(0.0, wsats_infilt * (1.0 - totwcont_infilt)) + wsats_infilt *
					pow(max(0.0, pow(1.0 - totwcont_infilt, 1.0 / (1.0 + soiltype.b_oro)) -
						winflux_sd / (wsats_infilt * (1.0 + soiltype.b_oro))), soiltype.b_oro + 1.0));
			}
		}
		else {*/
			// Surface runoff is equal to excees water above infiltation capacity
			runoff_surf_sd = max(0.0, min(winflux_sd, (winflux_sd + wmass_infilt + imass_infilt) - wsats_infilt));
		//}
		runoff_surf += runoff_surf_sd / wtimesteps;

		// Calculate infiltration (mm/day)
		double winfilt_sd = winflux_sd - runoff_surf_sd;
		winfilt += winfilt_sd / wtimesteps;

		// Calculate water flows and content after knowing surface influx and AET of each soil layer
		calc_wflow(winfilt_sd, aet_layer);
	}

	// Calculate water fluxes to, inbetween and out of soil layers and
	// check for negative wmass

	// Soil layers with negative water mass
	std::vector<int> neglyr;

	wfluxes[0] = winfilt;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wfluxes[lyr + 1] = wfluxes[lyr] + wpara.psi0[lyr] * soillayers[lyr].depth - (soillayers[lyr].wmass + soillayers[lyr].runoff_drain + aet_layer[lyr]);

		if (soillayers[lyr].wmass < 0.0)
			neglyr.push_back(lyr);
	}

	// Calculate baseflow
	soillayers[nsoillayers - 1].runoff_baseflow = max(0.0, wfluxes[nsoillayers]);

	// Negative soil layer wmass exists
	if (neglyr.size()) {
		resolve_neg_wmass(neglyr, evap, runoff_surf);
	}

	// Water balance check
	double final_wmass = snow.wcontent() + get_sum(WMASS) + runoff_surf + get_sum(RUNOFFDRAIN) + get_sum(RUNOFFBASE) + evap + aet_total;

	if (fabs(initial_wmass - final_wmass) > 0.0) {

		double diff = final_wmass - initial_wmass;

		// Too little water
		if (diff < 0.0) {

			// Add excess water to top soil layer
			soillayers[0].wmass += diff;
		}
		// Too much water
		else {

			// Take from bottom soil layer upwards
			int lyr = nsoillayers - 1;
			while (diff > 0.0 && lyr >= 0) {
				
				double water = min(soillayers[lyr].wmass, diff);
				diff -= water;
				soillayers[lyr].wmass -= water;
				lyr--;
			}
		}
	}

	final_wmass = snow.wcontent() + get_sum(WMASS) + runoff_surf + get_sum(RUNOFFDRAIN) + get_sum(RUNOFFBASE) + evap + aet_total;
	if (fabs(initial_wmass - final_wmass) > 1.0e-10 && date.year > FIRST_FREEZE_YEAR)
		dprintf("Year %d Day %d WMASS ERROR: Diff %g Begin %g End %g\n", date.year, date.day, initial_wmass - final_wmass, initial_wmass, final_wmass);

	// Update soil layer temperatures due to soil water fluxes
	wflow_temp();

	// Reset ran_melt as water has entered the soil
	rain_melt = 0.0;

	// Save fluxes
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, get_sum(RUNOFFDRAIN));
	patch.fluxes.report_flux(Fluxes::WRUNOFF_BASE, get_sum(RUNOFFBASE));
	patch.fluxes.report_flux(Fluxes::WEVAP, evap);
	patch.fluxes.report_flux(Fluxes::WAET, aet_total);

	// Calculate annual wcont for growing season (temperature > 5C)
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Reset the awcont array on the first day of every year
		if (date.day == 0) {
			soillayers[lyr].awcont = 0.0;
			if (lyr == 0) patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			soillayers[lyr].awcont += soillayers[lyr].wcont();
			if (lyr == 0) patch.growingseasondays++;
		}

		// Do the averaging on the last day of every year
		if (date.islastday && date.islastmonth) {

			// In case it's never warm enough:
			if (!patch.growingseasondays) {
				soillayers[lyr].awcont = 1.0;
			}
			else {
				soillayers[lyr].awcont /= (double)patch.growingseasondays;
			}
		}
	}
}


// Initial percolation
/** Initial percolation for two layer soil scheme
 * and for peat and wetlands
 */
void Soilcolumn::initial_percolation() {

	patch.wetland_wmass_added_today = 0.0;

	if (percolate) {

		if (iftwolayersoil) {

			// As in LPJ-GUESS v4.0, where soil.wcont[0] and soil.wcont_evap were adjusted

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
		else {

			if (patch.stand.is_highlatitude_peatland_stand()) {

				// Peatland soils
				// Distribute the water in the acrotelm in proportion to the capacity

				// water that can still be added to each soil layer (kg m-2)
				std::vector<double> potential_layer(nacrotelm);

				double imass_fraction = 0.0;
				double imass_impedance = 1.0;

				double imass_tot = 0.0;
				double wsats_tot = 0.0;

				// Total potential of the soil to store water, from the current content up to saturation
				double total_potential = 0.0;

				for (int lyr = 0; lyr < nacrotelm; lyr++) {

					Soillayer& soillayer = soillayers[lyr];

					potential_layer[lyr] = max(0.0, soillayer.wsats - (soillayer.wmass + soillayer.imass));

					wsats_tot += soillayer.wsats;
					imass_tot += soillayer.imass;
					total_potential += potential_layer[lyr];
				}

				// Potential for infiltration exist
				if (total_potential > 0.0) {

					imass_fraction = imass_tot / wsats_tot;

					// Impedance - this is the fraction of water that can infiltrate initially 
					// imass_fraction = 0 (e.g. summer, no phase change etc.)
					// imass_fraction = 1 (e.g. winter, saturated soil to field capacity, all soil water frozen)
					// CLM4.5 - see Swenson et al (2012) Eqn. (8)
					if (IMASS_IMPEDANCE && can_freeze()) {

						imass_impedance = max(0.0, 1.0 - imass_fraction);
						// CLM method: imass_impedance = pow(10.0, -6.0 * imass_fraction);
					}
					else {
						imass_impedance = 1.0;
					}

					double wmass_in = 0.0;
					double rain_melt_orig = rain_melt;

					// Assume the amount of water available for initial infiltration is limited by the ice content in the top layer
					double wmass_for_infiltration = imass_impedance * rain_melt;

					// Only infiltrate what we can. The rest remains in rain_melt.
					if (wmass_for_infiltration >= total_potential) {
						wmass_in = total_potential;
						rain_melt -= total_potential;
					}
					else {
						// All water can infiltrate
						wmass_in = wmass_for_infiltration;
						rain_melt -= wmass_for_infiltration; // i.e. 0.0 when there is no ice
					}

					// Now add the rain_melt to each layer in proportion to the capacity if total_potential > 0.0 mm
					for (int lyr = 0; lyr < nacrotelm; lyr++) {

						soillayers[lyr].wmass += wmass_in * (potential_layer[lyr] / total_potential);

						if (soillayers[lyr].wmass - soillayers[lyr].wsats > 1.0e-12)
							fail("ERROR: Year %d Day %d Lyr %d initial_infiltration - TOO much water added! - is_highlatitude_peatland_stand\n", date.year,date.day,lyr);
					}
				}
			}
		}
	}

	if (patch.stand.is_true_wetland_stand()) {

		// Mineral wetlands

		// Water mass added to saturate soil
		double wmass_added = 0.0;

		// Saturate the soil
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			double wmass_input = soillayer.wsats - (soillayer.wmass + soillayer.imass);

			// Add water to the layer
			soillayer.wmass += wmass_input;
			wmass_added += wmass_input;
		}

		// Record the water added to this wetland today
		patch.wetland_wmass_added_today = max(0.0, wmass_added - rain_melt);

		rain_melt -= min(rain_melt, wmass_added);

		// End mineral wetlands
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

	// Balance ice and water in snow layers and snow layer depths after phase change
	snow.balance_snow(patch, rain_melt, rain_melt_temp);

	// Update snow thermal properties
	snow.update_thermal_properties();


	// Soil phase change

	// No phase change for years before first freeze year
	if (date.year >= FIRST_FREEZE_YEAR) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].phase_change();
		}
	}
}

void Soilcolumn::phase_change_new() {

	// Snow phase change

	// Any updates to the snowpack
	bool updated = false;

	// Phase change between water and ice and adjust snow layer temperatures accordingly
	snow.phase_change_new(updated);

	// Phase change has happened and we need to update the snowpack
	if (updated) {

		// Balance ice and water in snow layers and snow layer depths after phase change
		snow.balance_snow(patch, rain_melt, rain_melt_temp);

		// Update snow thermal properties
		snow.update_thermal_properties();
	}

	// Soil phase change

	// No phase change for years before first freeze year
	if (date.year >= FIRST_FREEZE_YEAR) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].phase_change_new();
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
		while (lyr < nsoillayers + PAD_LAYERS && pad[lyr - nsoillayers].temp > 0.0) {
			thaw += pad[lyr - nsoillayers].depth;
			lyr++;
		}
	}

	// If still among soil layers then assuming that water in next soil layer is at the top of the layer and ice below
	if (lyr < nsoillayers) {
		if (!negligible(soillayers[lyr].wmass + soillayers[lyr].imass)) {
			thaw += soillayers[lyr].wmass / (soillayers[lyr].wmass + soillayers[lyr].imass) * soillayers[lyr].depth;
		}
	}

	// Remove very small negative values
	thaw = max(0.0, thaw);

	// Set maximum thaw depth this year
	if (thaw > maxthaw || date.day == 0)
		maxthaw = thaw;

	// Calculate mean monthly thaw depth
	mthaw[date.month] += thaw / date.ndaymonth[date.month];

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

	double airT = patch.get_climate().temp;
	double top_temp = snow.nsnowlayers ? snow.snowlayers[snow.nsnowlayers - 1].temp : soillayers[0].temp;

	if (fabs(airT - top_temp) > 2.0) {
		ttimesteps = 48;
	}
	else {
		ttimesteps = 24;
	}

	// Update heat capacity and thermal conductivity for snow, soil and padding layers
	update_thermal_properties();

	/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
		if (date.year == FIRST_FREEZE_YEAR && date.day == 0) {
			dprintf("Year Day TS AirT ");
			for (int lyr = 0; lyr < nsoillayers; lyr++)
				dprintf("Lyr%d ", lyr);
			dprintf("ENERGY\n");
		}
		dprintf("%d %d %d %g ",date.year, date.day, 0, airT);
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			dprintf("%g ", soillayers[lyr].temp);
			cum_energy += soillayers[lyr].Ch * (soillayers[lyr].temp + K2degC) * soillayers[lyr].depth * M_PER_MM;
		}
		dprintf("%g\n", cum_energy);
	}*/


	for (int ts = 0; ts < ttimesteps; ts++) {

		int nsnowlayers = snow.nsnowlayers;
		bool combine_snow = false;

		// If snow layer is too thin add its depth and properties to the top soil layer (Best 2011 section 3.3)
		if (snow.nsnowlayers == 1) {
			if (snow.depth() < SNOW_DEPTH_STANDARD) {
				combine_snow = true;
				nsnowlayers = 0;

				// Calculate the combined top soil layer and single snow layers thermal properties
				soillayers[0].update_thermal_properties(combine_snow, snow.snowlayers[0].depth, snow.snowlayers[0].imass, snow.snowlayers[0].wmass);
			}
		}

		if (false) {

			int nsnowsoillayers = nsnowlayers + nsoillayers;
			int nsnowsoilpadlayers = nsnowlayers + nsoillayers + PAD_LAYERS;

			if (nsnowsoilpadlayers != tpara.dz.size()) {
				tpara.set_sizes(nsnowsoilpadlayers);
			}

			// Snow layers (0 == bottom layer, nsnowlayers - 1 == top snow layer)
			int slyr = nsnowlayers - 1;
			for (int lyr = 0; lyr < nsnowlayers; lyr++) {
				tpara.dz[lyr] = snow.snowlayers[slyr].depth * M_PER_MM;
				tpara.cap[lyr] = snow.snowlayers[slyr].Ch;
				tpara.diff[lyr] = snow.snowlayers[slyr].Kh;
				tpara.psi[lyr] = snow.snowlayers[slyr].temp;

				slyr--;
			}

			// Soil layers
			for (int lyr = nsnowlayers; lyr < nsnowsoillayers; lyr++) {
				if (lyr == nsnowlayers && combine_snow) {
					tpara.dz[lyr] = (soillayers[lyr - nsnowlayers].depth + snow.snowlayers[0].depth) * M_PER_MM;
				}
				else {
					tpara.dz[lyr] = soillayers[lyr - nsnowlayers].depth * M_PER_MM;
				}
				tpara.cap[lyr] = soillayers[lyr - nsnowlayers].Ch;
				tpara.diff[lyr] = soillayers[lyr - nsnowlayers].Kh;
				tpara.psi[lyr] = soillayers[lyr - nsnowlayers].temp;
			}

			// Padding layers
			for (int lyr = nsnowsoillayers; lyr < nsnowsoilpadlayers; lyr++) {
				tpara.dz[lyr] = pad[lyr - nsnowsoillayers].depth * M_PER_MM;
				tpara.cap[lyr] = pad[lyr - nsnowsoillayers].Ch;
				tpara.diff[lyr] = pad[lyr - nsnowsoillayers].Kh;
				tpara.psi[lyr] = pad[lyr - nsnowsoillayers].temp;
			}

			double depth = tpara.dz[0] / 2.0;
			tpara.z[0] = depth;
			for (int lyr = 1; lyr < nsnowsoilpadlayers; lyr++) {
				depth += (tpara.dz[lyr - 1] + tpara.dz[lyr]) / 2.0;
				tpara.z[lyr] = depth;
			}


			// Boundary conditions for soil temperature (Shouldn't use air temperature)
			// Upper boundary condition [J m-2 s-1]
			double f_top = -2.0 * tpara.diff[0] * (tpara.psi[0] - airT) / tpara.z[0];

			// Lower boundary condition
			double f_bot = 0.0;

			// Initial energy mass balance
			double init_cum_energy = f_top * SECS_PER_DAY / ttimesteps;
			for (int lyr = 0; lyr < nsnowsoilpadlayers; lyr++) {
				init_cum_energy += (tpara.psi[lyr] + K2degC) * tpara.cap[lyr] * tpara.dz[lyr];
			}

			// Build diagonals and RHS vector
			build_diags_rhs(SECS_PER_DAY / ttimesteps, f_top, f_bot, tpara, nsnowsoilpadlayers);

			// Solve tridiagonal system
			solve_tridiagonal(tpara);

			// Final energy mass balance
			double final_cum_energy = 0.0;
			for (int lyr = 0; lyr < nsnowsoilpadlayers; lyr++) {
				final_cum_energy += (tpara.psi1[lyr] + K2degC) * tpara.cap[lyr] * tpara.dz[lyr];
			}

			if (fabs(final_cum_energy - init_cum_energy) > 1.0e-12 * final_cum_energy)
				dprintf("Year %d Day %d TS %d Diff %g Final %g Initial %g\n", date.year, date.day, ts, fabs(final_cum_energy - init_cum_energy), final_cum_energy, init_cum_energy);

			// Update temperatures

			// Set new snow layer temperatures (0 == bottom layer, nsnowlayers - 1 == top snow layer)
			slyr = nsnowlayers - 1;
			for (int lyr = 0; lyr < nsnowlayers; lyr++) {
				snow.snowlayers[slyr].temp = tpara.psi1[lyr];
				slyr--;
			}

			// Set new soil layer temperatures
			for (int lyr = nsnowlayers; lyr < nsnowsoillayers; lyr++) {
				soillayers[lyr - nsnowlayers].temp = tpara.psi1[lyr];

				if (lyr == nsnowlayers && combine_snow)
					snow.snowlayers[0].temp = tpara.psi1[lyr];
			}

			// Set new padding layer temperatures
			for (int lyr = nsnowsoillayers; lyr < nsnowsoilpadlayers; lyr++) {
				pad[lyr - nsnowsoillayers].temp = tpara.psi1[lyr];
			}
		}
		else {

			ttimesteps = 4;
			int nairlayers = 1;
			int nairsnowlayers = nairlayers + nsnowlayers;
			int nairsnowsoillayers = nairlayers + nsnowlayers + nsoillayers;
			int nairsnowsoilpadlayers = nairlayers + nsnowlayers + nsoillayers + PAD_LAYERS;

			if (nairsnowsoilpadlayers != tpara.dz.size()) {
				tpara.set_sizes(nairsnowsoilpadlayers);
			}

			// Air layer
			tpara.dz[0] = AIR_THICKNESS;
			tpara.cap[0] = Cp_air * M3_PER_MM3;
			tpara.diff[0] = Kair * SECS_PER_DAY / MM_PER_M;
			tpara.cond[0] = tpara.diff[0] / tpara.cap[0];
			tpara.psi[0] = airT;

			// Snow layers (0 == bottom layer, nsnowlayers - 1 == top snow layer)
			int slyr = nsnowlayers - 1;
			for (int lyr = nairlayers; lyr < nairsnowlayers; lyr++) {
				tpara.dz[lyr] = snow.snowlayers[slyr].depth;
				tpara.cap[lyr] = snow.snowlayers[slyr].Ch / MM3_PER_M3; // m-3 to mm-3
				tpara.diff[lyr] = snow.snowlayers[slyr].Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
				tpara.psi[lyr] = snow.snowlayers[slyr].temp;

				slyr--;
			}

			// Soil layers
			for (int lyr = nairsnowlayers; lyr < nairsnowsoillayers; lyr++) {
				if (lyr == nairsnowlayers && combine_snow) {
					tpara.dz[lyr] = (soillayers[lyr - nairsnowlayers].depth + snow.snowlayers[0].depth);
				}
				else {
					tpara.dz[lyr] = soillayers[lyr - nairsnowlayers].depth;
				}
				tpara.cap[lyr] = soillayers[lyr - nairsnowlayers].Ch / MM3_PER_M3; // m-3 to mm-3
				tpara.diff[lyr] = soillayers[lyr - nairsnowlayers].Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
				tpara.psi[lyr] = soillayers[lyr - nairsnowlayers].temp;
			}

			// Padding layers
			for (int lyr = nairsnowsoillayers; lyr < nairsnowsoilpadlayers; lyr++) {
				tpara.dz[lyr] = pad[lyr - nairsnowsoillayers].depth;
				tpara.cap[lyr] = pad[lyr - nairsnowsoillayers].Ch / MM3_PER_M3; // m-3 to mm-3
				tpara.diff[lyr] = pad[lyr - nairsnowsoillayers].Kh * SECS_PER_DAY / MM_PER_M; // s to day, m-1 to mm-1
				tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
				tpara.psi[lyr] = pad[lyr - nairsnowsoillayers].temp;
			}

			cnstep_full(airT, 1.0 / ttimesteps, tpara);


			// Update temperatures

			// Set new snow layer temperatures (0 == bottom layer, nsnowlayers - 1 == top snow layer)
			slyr = nsnowlayers - 1;
			for (int lyr = nairlayers; lyr < nairsnowlayers; lyr++) {
				snow.snowlayers[slyr].temp = tpara.psi1[lyr];
				slyr--;
			}

			// Set new soil layer temperatures
			for (int lyr = nairsnowlayers; lyr < nairsnowsoillayers; lyr++) {
				soillayers[lyr - nairsnowlayers].temp = tpara.psi1[lyr];

				if (lyr == nairsnowlayers && combine_snow)
					snow.snowlayers[0].temp = tpara.psi1[lyr];
			}

			// Set new padding layer temperatures
			for (int lyr = nairsnowsoillayers; lyr < nairsnowsoilpadlayers; lyr++) {
				pad[lyr - nairsnowsoillayers].temp = tpara.psi1[lyr];
			}
		}

		/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
			dprintf("%d %d %d %g ", date.year, date.day, 2*ts + 1, airT);
			cum_energy = 0.0;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				dprintf("%g ", soillayers[lyr].temp);
				cum_energy += soillayers[lyr].Ch * (soillayers[lyr].temp + K2degC) * soillayers[lyr].depth * M_PER_MM;
			}
			dprintf("%g\n", cum_energy);
		}*/


		/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
			dprintf("%d %d %d %g ", date.year, date.day, ts, airT);

			for (int lyr = 0; lyr < NSNOWLAYERS_MAX - snow.nsnowlayers; lyr++) {
				dprintf("%g ", airT);
			}

			for (int lyr = snow.nsnowlayers - 1; lyr >= 0; lyr--) {
				dprintf("%g ", snow.snowlayers[lyr].temp);
			}

			for (int lyr = 0; lyr < 5; lyr++) {
				dprintf("%g %g %g %g %g ", soillayers[lyr].temp, soillayers[lyr].temp_old, soillayers[lyr].wmass, soillayers[lyr].imass, soillayers[lyr].wp);
			}

			dprintf("\n");
		}*/

		// Recalculate the top soil layer thermal properties after it been combined with the single snow layer for temperature calculations
		if (combine_snow) {
			soillayers[0].update_thermal_properties();
		}

		// Phase change between water and ice and adjust snow and soil layer temperatures accordingly
		phase_change_new();


		/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
			dprintf("%d %d %d %g ", date.year, date.day, 2 * ts + 2, airT);
			cum_energy = 0.0;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				dprintf("%g ", soillayers[lyr].temp);
				cum_energy += soillayers[lyr].Ch * (soillayers[lyr].temp + K2degC) * soillayers[lyr].depth * M_PER_MM;
			}
			dprintf("%g\n", cum_energy);
		}*/

		/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
			dprintf("%d %d %d %g ", date.year, date.day, ts, airT);

			for (int lyr = 0; lyr < NSNOWLAYERS_MAX - snow.nsnowlayers; lyr++) {
				dprintf("%g ", airT);
			}

			for (int lyr = snow.nsnowlayers - 1; lyr >= 0; lyr--) {
				dprintf("%g ", snow.snowlayers[lyr].temp);
			}

			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				dprintf("%g ", soillayers[lyr].temp);
			}

			for (int lyr = 0; lyr < PAD_LAYERS; lyr++) {
				dprintf("%g ", pad[lyr].temp);
			}

			dprintf("\n");
		}*/

		/*dprintf("%d %d %d %g ", date.year, date.day, ts, airT);

		for (int lyr = 0; lyr < NSNOWLAYERS_MAX - snow.nsnowlayers; lyr++) {
			dprintf("%g ", airT);
		}

		for (int lyr = snow.nsnowlayers - 1; lyr >= 0; lyr--) {
			dprintf("%g ", snow.snowlayers[lyr].temp);
		}

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			dprintf("%g ", soillayers[lyr].temp);
		}

		for (int lyr = 0; lyr < PAD_LAYERS; lyr++) {
			dprintf("%g ", pad[lyr].temp);
		}

		dprintf("\n");*/
	}


	/*if (date.year >= FIRST_FREEZE_YEAR && date.year <= FIRST_FREEZE_YEAR + 1) {
		dprintf("%d %d %d %g ", date.year, date.day, 99, airT);

		for (int lyr = 0; lyr < NSNOWLAYERS_MAX - snow.nsnowlayers; lyr++) {
			dprintf("%g ", airT);
		}

		for (int lyr = snow.nsnowlayers - 1; lyr >= 0; lyr--) {
			dprintf("%g ", snow.snowlayers[lyr].temp);
		}

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			dprintf("%g ", soillayers[lyr].temp);
		}

		for (int lyr = 0; lyr < PAD_LAYERS; lyr++) {
			dprintf("%g ", pad[lyr].temp);
		}

		dprintf("\n");
	}*/

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

		for (int lyr = PAD_LAYERS-2; lyr < PAD_LAYERS; lyr++) {
			pad[lyr].temp = avr_temp;
		}
	}
}


/// Calculate snow, soil and padding layer temperatures
/** Calculate snow, soil and padding layer temperatures
*/
void Soilcolumn::temp() {

	if (iftwolayersoil) {
		// Calculate analytic soil temperatures and update dtemp_31 (trend in daily air temperatures for the last 31 days)
		// Analytic soil temperatures can be used to test the accuracy of the numerical scheme in Soilcolumn::calc_temp below.
		calc_temp_analytic();
	}
	else {

		// Calculate the soil temperature in snow, soil and padding layers.
		calc_temp();
	}

	// Calculate mean monthly soil layer temperatures
	calc_mtemp();
}

/// Saves snow properties
/** Saves snow properties
*/
void Soilcolumn::save_snow() {

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
}

/// Vertical movement of organic matter if C density of soil layer exceed maximum
/** Vertical movement of organic matter if C density of soil layer exceed maximum
*  Downwards movement of SOM if soil layer C exceeds maximum
*/
void Soilcolumn::vertical_movement() {

	// Root C biomass in soil layers (currently set to zero) // DW_COMMENT future development: include root in soil layers as influence to total volume, thermal properties, hydrology etc
	std::vector<double> ccont_roots(nsoillayers);

	int donor, receiver, p_receiver;;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		double ccont_soil = soillayers[lyr].ccont(true);

		// If soil layer C exceed maximum, then move SOM downwards
		if (ccont_soil + ccont_roots[lyr] > Orgdensity * soillayers[lyr].depth * M_PER_MM) {
			donor = lyr;
			receiver = lyr + 1;

			if (receiver < nsoillayers) {

				double ccont_donor = ccont_soil + ccont_roots[lyr];
				double ccont_receiver = soillayers[receiver].ccont(true) + ccont_roots[receiver];

				if (!negligible(ccont_soil)) {

					double move_frac = 0.0;

					// Downwards movement of SOM
					move_frac = (ccont_donor - Orgdensity * soillayers[donor].depth * M_PER_MM) / ccont_soil;

					for (int p = 0; p < NSOMPOOL; p++) {

						if (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
							double cmove = soillayers[donor].sompool[p].cmass * move_frac;
							double nmove = soillayers[donor].sompool[p].nmass * move_frac;

							soillayers[donor].sompool[p].cmass -= cmove;
							soillayers[donor].sompool[p].nmass -= nmove;

							if (p == SURFHUMUS)
								p_receiver = SLOWSOM;
							else if (SURFMICRO)
								p_receiver = SOILMICRO;
							else
								p_receiver = p;

							soillayers[receiver].sompool[p_receiver].cmass += cmove;
							soillayers[receiver].sompool[p_receiver].nmass += nmove;
						}
					}
				}
			}
			else {
				// Last layer has too much soil C 
				// Instead of creating a new layer, the existing last layer
				// is expanded by 100mm. If new layer was created then it would
				// create issues for output files, the need to update all 
				// individuals and with LUC
				double depth_increase = 100.0;
				soillayers[donor].depth += depth_increase;
				soillayers[donor].depth_profile += depth_increase / 2.0;
				soillayerdepths[donor] += depth_increase;

				// Update affected soil properties
				double wmass_add = 0.0;
				double imass_add = 0.0;
				soillayers[donor].set_properties(wmass_add, imass_add, true);
			}
		}
		else if (patch.stand.landcover == PEATLAND && ccont_soil + ccont_roots[lyr] < Orgdensity * soillayers[lyr].depth * M_PER_MM) {
			// Check if SOM can be moved upwards

			donor = lyr + 1;
			receiver = lyr;

			if (donor < nsoillayers) {

				double ccont_donor = soillayers[donor].ccont(true);
				double ccont_receiver = ccont_soil + ccont_roots[lyr];

				if (!negligible(ccont_donor)) {

					double move_frac = min(Orgdensity * soillayers[receiver].depth * M_PER_MM - ccont_receiver, ccont_donor) / ccont_donor;

					for (int p = 0; p < NSOMPOOL; p++) {

						if (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
							double cmove = soillayers[donor].sompool[p].cmass * move_frac;
							double nmove = soillayers[donor].sompool[p].nmass * move_frac;

							soillayers[donor].sompool[p].cmass -= cmove;
							soillayers[donor].sompool[p].nmass -= nmove;

							soillayers[receiver].sompool[p].cmass += cmove;
							soillayers[receiver].sompool[p].nmass += nmove;
						}
					}
				}
			}
		}
	}
}

/// Vertical mixing (diffusion) of soil organic matter
/** Vertical mixing (diffusion) of soil organic matter
*  Cryoturbation processes on the vertical profile of SOM stocks are captured
*  by considering a simple diffusion scheme.
*/
void Soilcolumn::vertical_mixing(pooltype p, double ald, double timestep) {

	// Fixed multiplier to the active layer depth for maximum mixing depth
	double k = 3.0;

	// active layer depth is not deep enough to create any vertical mixing then return
	if (k * ald < soillayers[0].depth_profile * M_PER_MM)
		return;	

	// Diffusion constant resulting from the mixing due to cryoturbation (m2 yr-1)
	double Dc_0 = 0.0001;
	
	// Top and bottom boundary conditions are zero for soil mixing
	double f_top = 0.0;
	double f_bot = 0.0;

	CNpara VMpara;

	VMpara.set_sizes(nsoillayers);

	// Populate Crank-Nicholson struct
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		VMpara.dz[lyr] = soillayers[lyr].depth * M_PER_MM;
		VMpara.z[lyr] = soillayers[lyr].depth_profile * M_PER_MM;

		// SOM pool C concentration
		VMpara.psi[lyr] = soillayers[lyr].sompool[p].cmass / VMpara.dz[lyr];

		// Depth of mixing is a function of the active layer depth, with mixing rates 
		// constant through the active layer and tapering linearly to zero at a depth of some 
		// constant (k) times the active layer depth (Koven 2009 eqn 3)
		if (VMpara.z[lyr] < ald) {
			VMpara.diff[lyr] = Dc_0;
		}
		else if (VMpara.z[lyr] > k * ald) {
			VMpara.diff[lyr] = 0.0;
		}
		else {
			VMpara.diff[lyr] = Dc_0 * (1.0 - min(1.0, (VMpara.z[lyr] - ald) / ((k - 1.0) * ald)));
		}
	}

	// Build diagonals and RHS vector
	build_diags_rhs(timestep, f_top, f_bot, VMpara, nsoillayers);

	// Solve tridiagonal system
	solve_tridiagonal(VMpara);

	// Calculate carbon fluxes between soil layers
	// Positive flow to layer below
	// Negative flow from layer below
	std::vector<double> cfluxes(nsoillayers - 1);
	cfluxes[0] = (VMpara.psi[0] - VMpara.psi1[0]) * VMpara.dz[0];
	cfluxes[0] = negligible(cfluxes[0], 1.0e-12) ? 0.0 : cfluxes[0];
	for (int lyr = 1; lyr < nsoillayers - 1; lyr++) {
		cfluxes[lyr] = cfluxes[lyr - 1] + (VMpara.psi[lyr] - VMpara.psi1[lyr]) * VMpara.dz[lyr];
		cfluxes[lyr] = negligible(cfluxes[lyr], 1.0e-12) ? 0.0 : cfluxes[lyr];
	}

	// Check for inconsistencies in downwards flux, no negative C mass
	for (int lyr = 0; lyr < nsoillayers - 1; lyr++) {
		if (lyr == 0) {
			if (soillayers[lyr].sompool[p].cmass - cfluxes[lyr] < 0.0) {
				cfluxes[lyr] = soillayers[lyr].sompool[p].cmass;
			}
		}
		else {
			if (soillayers[lyr].sompool[p].cmass + cfluxes[lyr - 1] - cfluxes[lyr] < 0.0 && cfluxes[lyr] > 0.0) {
				double neg_c = soillayers[lyr].sompool[p].cmass + cfluxes[lyr - 1] - cfluxes[lyr];
				cfluxes[lyr] += max(neg_c, -cfluxes[lyr]);
			}
		}
	}

	// Check for inconsistencies in upwards flux, no negative C mass
	for (int lyr = 1; lyr < nsoillayers - 1; lyr++) {
		if (soillayers[lyr].sompool[p].cmass + cfluxes[lyr - 1] - cfluxes[lyr] < 0.0) {
			if (cfluxes[lyr - 1] < 0.0) {
				double neg_c = soillayers[lyr].sompool[p].cmass + cfluxes[lyr - 1] - cfluxes[lyr];
				cfluxes[lyr - 1] -= neg_c;
			}
			if (cfluxes[lyr - 1] > 0.0)
				dprintf("ERROR: VERTICAL_MIX lyr %d flux should be negative or zero %g\n", lyr, cfluxes[lyr - 1]);
		}
	}

	std::vector<double> delta_cmass(nsoillayers);
	std::vector<double> delta_nmass(nsoillayers);

	// Calculate C and N mass changes
	double cmove, nmove, frac;
	for (int lyr = 0; lyr < nsoillayers - 1; lyr++) {
		if (cfluxes[lyr] > 0.0 && !negligible(soillayers[lyr].sompool[p].cmass)) {
			frac = max(0.0, cfluxes[lyr] / soillayers[lyr].sompool[p].cmass);
			cmove = frac * soillayers[lyr].sompool[p].cmass;
			nmove = frac * soillayers[lyr].sompool[p].nmass;

			delta_cmass[lyr] -= cmove;
			delta_nmass[lyr] -= nmove;
			delta_cmass[lyr + 1] += cmove;
			delta_nmass[lyr + 1] += nmove;
		}
		else if (cfluxes[lyr] < 0.0 && !negligible(soillayers[lyr + 1].sompool[p].cmass)) {
			frac = max(0.0, -cfluxes[lyr] / soillayers[lyr + 1].sompool[p].cmass);
			cmove = frac * soillayers[lyr + 1].sompool[p].cmass;
			nmove = frac * soillayers[lyr + 1].sompool[p].nmass;

			delta_cmass[lyr + 1] -= cmove;
			delta_nmass[lyr + 1] -= nmove;
			delta_cmass[lyr] += cmove;
			delta_nmass[lyr] += nmove;
		}
	}

	// Make C and N changes
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].sompool[p].cmass += delta_cmass[lyr];
		soillayers[lyr].sompool[p].nmass += delta_nmass[lyr];
	}
}

/// Serialize soilcolumn variables
/** Serialize soilcolumn variables
 */
void Soilcolumn::serialize(ArchiveStream& arch) {
	arch & soillayers
		& anfix
		& anfix_mean
		& mald_mean
		& awtp
		& anfix_calc
		& pad
		& solvesomcent_beginyr
		& solvesomcent_endyr;
}



///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Campbell, G.S.: Simple method for determining unsaturated conductivity from moisture 
//  retention data. Soil Sci. 117: 311–314. 
//  https://doi.org/10.1097/00010694-197406000-00001, 1974
// Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., 
//  and Friedlingstein, P.: An improved representation of physical permafrost dynamics 
//  in the JULES land - surface model, Geosci.Model Dev., 8, 1493 - 1508, 
//  https ://doi.org/10.5194/gmd-8-1493-2015, 2015
// Clapp, R.B., and Hornberger, G.M. Empirical equations for some soil hydraulic properties. 
//  Water Resour. Res. 14: 601–604. https://doi.org/10.1029/WR014i004p00601, 1978
// Fukusako, S. Int J Thermophys (1990) 11: 353. doi:10.1007/BF01133567
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//  and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Lawrence, D. M. and Slater, A.: Incorporating organic soil into a global climate 
//  model, Clim.Dynam., 30, 145–160, doi:10.1007/s00382-007-0278-1, 2008.
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
// Wang, T., C. Ottle´ , A. Boone, P. Ciais, E. Brun, S. Morin, G. Krinner, S. Piao, 
//   and S. Peng (2013), Evaluation of an improved intermediate complexity snow scheme 
//   in the ORCHIDEE land surface model, J. Geophys. Res. Atmos., 118, 6064–6079, 
//   doi:10.1002/jgrd.50395
// For Arctic and wetland references see soil.h
