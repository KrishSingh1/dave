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

////////////////////////////////////////////////////////////////////////////////
// Implementation of Snowlayer member functions
////////////////////////////////////////////////////////////////////////////////

/// Clears all members
void Snowlayer::clear() {
	imass = 0.0;
	wmass = 0.0;
	dens = 0.0;
	depth = 0.0;
	whc = 0.0;
	temp = 0.0;
	temp_old = 0.0;
	Kh = 0.0;
	Ch = 0.0;
	Dh = 0.0;
}

double snowheatcap(double temp, double dens) {
	// C Snow - Fukusako, Eq. 2. [J m-3 K-1]
	return (0.185 + 0.689 * (K2degC + temp) * 0.01) * dens * J_PER_KJ;
}

///
/*
 *
 */
void Snowlayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible((imass + wmass) * Ch + add_mass * add_Ch))
		dprintf("Year %d Day %d ERRRRRRRRRROOOORRRRRRRRR \n",date.year,date.day);

	if (negligible(add_mass))
		int sch = 0;

	// New snow layer but haven't updated snow properties yet 
	if (!negligible(imass) && negligible(Ch)) {
		update_heatcap_and_cond();
	}

	if (negligible(imass) && negligible(Ch)) {
		Ch = add_Ch;
	}

	// Add to existing snow layer with updated snow properties
	if (!negligible(Ch) && !negligible(imass)) {
		temp = (temp * (imass + wmass) * Ch + add_temp * add_mass * add_Ch) / ((imass + wmass) * Ch + add_mass * add_Ch);
	}
	else { 
		// No existing snow layer, so added snow temperature is new snow layer temperature
		temp = add_temp;
	}
}

///
/*
 *
 */
void Snowlayer::update_heatcap_and_cond() {

	double max_temp = 50.0; // DW_COMMENT
	double min_temp = -50.0;
	double snow_temp = min(max_temp, max(min_temp, temp));

	// Snow density including water withing snow
	double tot_dens = min(rho_ICE, (imass + wmass) / depth * MM_PER_M);

	// C Snow - Fukusako, Eq. 2. [J m-3 K-1]
	Ch = snowheatcap(snow_temp, tot_dens);

	// K Snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = 2.22 * pow((tot_dens / rho_H2O), 1.88);

	// D Snow
	Dh = Kh / Ch;
}

// Update snow water holding capacity
void Snowlayer::update_whc() {
	// Calculate layer water holding capacity, based on eqn 10 Wang et al. 2013	
	double dens_t = 400;	// kg m-3, reference density
	double rw_min = 0.03;
	double rw_max = 0.1;

	whc = imass * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - dens) / dens_t);
}

/// weigthed avg of new and old snow density
void Snowlayer::mix_dens(double prec, double snowfall_dens) {
	double sum = imass + wmass + prec;
	double w_new = prec / sum;
	double w_old = 1.0 - w_new;
	dens = min(rho_ICE, w_old * dens + w_new * snowfall_dens);
}

void Snowlayer::serialize(ArchiveStream& arch) {
	arch & imass
		& wmass
		& dens
		& depth
		& whc
		& temp
		& temp_old;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Snow member functions
////////////////////////////////////////////////////////////////////////////////

void Snow::init_state() {

	nsnowlayers = 0;
	NO3_mass = 0.0;
	NH4_mass = 0.0;
	snow_depth_standard = 100.0;
	snow_depth_limit = 2.0; // Twice the standard snow layer depth

	for (int mth = 0; mth < 12; mth++) {
		mdepth[mth] = 0.0;
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			mtemp[lyr][mth] = 0.0;
		}
	}
}


double Snow::nmass_avail() {
	return NO3_mass + NH4_mass;
}

double Snow::ccont() {
	return 0.0;
}

double Snow::ncont() {
	return NO3_mass + NH4_mass;
}

/// Water mass in snow layers
double Snow::wcontent() {
	double wcontent = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++)
		wcontent += snowlayers[lyr].wmass + snowlayers[lyr].imass;
	return wcontent;
}

double Snow::depth() {
	double sd = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sd += snowlayers[lyr].depth;
	}
	return sd;
}

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

		else if (mdepth[mth] > snow_max_mths[2])
			snow_max_mths[2] = mdepth[mth];
	}

	return (snow_max_mths[0] + snow_max_mths[1] + snow_max_mths[2]) / 3.0;
}

double Snow::snowpack() {
	double sp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sp += snowlayers[lyr].imass + snowlayers[lyr].wmass;
	}
	return sp;
}

void Snow::update_heatcap_and_cond() {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].update_heatcap_and_cond();
	}
}

/// Potential redistribution of ice/dens after melt
void Snow::potential_depth(double snowpack_depth, std::vector<double>& pot_depth, int& pot_nsnowlayers) {

	pot_nsnowlayers = 0;

	// All snow layer exceed standard depth. Bottom layer is thickest
	if (snowpack_depth >= snow_depth_standard * NSNOWLAYERS_MAX) {
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			pot_depth[lyr] = snow_depth_standard;
		}
		pot_depth[0] += snowpack_depth - snow_depth_standard * NSNOWLAYERS_MAX;

		pot_nsnowlayers = NSNOWLAYERS_MAX;

		//double add_snow_depth = snowpack_depth - snow_depth_standard * NSNOWLAYERS_MAX; // DW_COMMENT future snow development 2nd and 3rd snow layer should be the thick ones
		//pot_depth[1] += add_snow_depth / 2.0;
		//pot_depth[2] += add_snow_depth / 2.0;
	}
	else {
		// fill layers so that depths are between snow_depth_standard and snow_depth_standard * snow_depth_limit
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			if (snowpack_depth > snow_depth_standard * snow_depth_limit) {
				pot_depth[lyr] = snow_depth_standard;
				snowpack_depth -= snow_depth_standard;
				pot_nsnowlayers++;
			}
			else if (snowpack_depth <= snow_depth_standard * snow_depth_limit) {
				pot_depth[lyr] = snowpack_depth;
				snowpack_depth = 0.0;

				if (!negligible(pot_depth[lyr])) {
					pot_nsnowlayers++;
				}
			}
		}
	}
}

/// Re-distributing snow after freeze/melt or compaction event
void Snow::redist_layers(double snowpack_depth) {

	const double EPS = 1.0e-10;

	double initial_depth = 0.0;
	double initial_ice = 0.0;
	double initial_water = 0.0;
	double initial_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_depth += snowlayers[lyr].depth;
		initial_ice += snowlayers[lyr].imass;
		initial_water += snowlayers[lyr].wmass;
		initial_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	initial_temp /= initial_depth;


	// Potential depth
	std::vector<double> pot_depth(NSNOWLAYERS_MAX);

	// Potential number of snow layers
	int pot_nsnowlayers;

	// Check potential snow layer depth
	potential_depth(snowpack_depth, pot_depth, pot_nsnowlayers);

	// Re-calc ice & dens of layers
	for (int lyr = 0; lyr < nsnowlayers - 1; lyr++) {
		int lyr_above = lyr;

		double diff_lyr = fabs(pot_depth[lyr] - snowlayers[lyr].depth);

		// Move snow downwards in snowpack to fill every layer to pot_ice
		while (!negligible(diff_lyr) && lyr_above < nsnowlayers - 1) {
			lyr_above++;

			// layer shrinks so that depth needs to be moved from thicker bottom layer
			if (snowlayers[lyr].depth > pot_depth[lyr]) {

				double depth_move = snowlayers[lyr].depth - pot_depth[lyr];
				double frac_move = depth_move / snowlayers[lyr].depth;

				// Recalculate temperature
				snowlayers[lyr_above].update_temp(snowlayers[lyr].temp, frac_move * (snowlayers[lyr].imass + snowlayers[lyr].wmass), snowlayers[lyr].Ch);

				// Update depths
				snowlayers[lyr].depth -= depth_move;
				snowlayers[lyr_above].depth += depth_move;

				// Move ice accordingly
				snowlayers[lyr_above].imass += frac_move * snowlayers[lyr].imass;
				snowlayers[lyr].imass *= 1.0 - frac_move;

				// Move water accordingly
				snowlayers[lyr_above].wmass += frac_move * snowlayers[lyr].wmass;
				snowlayers[lyr].wmass *= 1.0 - frac_move;

				// Update density
				snowlayers[lyr_above].dens = snowlayers[lyr_above].imass * rho_H2O / snowlayers[lyr_above].depth;
			}
			else {

				double depth_move = max(0.0, min(snowlayers[lyr_above].depth, pot_depth[lyr] - snowlayers[lyr].depth));

				// Add depth to be moved if it isn't negligible
				if (!negligible(depth_move) && depth_move > 0.0) {

					double frac_move = depth_move / snowlayers[lyr_above].depth;

					// Recalculate temperature
					snowlayers[lyr].update_temp(snowlayers[lyr_above].temp, frac_move * (snowlayers[lyr_above].imass + snowlayers[lyr_above].wmass), snowlayers[lyr_above].Ch);

					// Update depths
					snowlayers[lyr_above].depth -= depth_move;
					snowlayers[lyr].depth += depth_move;

					// Move ice accordingly
					snowlayers[lyr].imass += frac_move * snowlayers[lyr_above].imass;
					snowlayers[lyr_above].imass *= 1.0 - frac_move;

					// Move water accordingly
					snowlayers[lyr].wmass += frac_move * snowlayers[lyr_above].wmass;
					snowlayers[lyr_above].wmass *= 1.0 - frac_move;

					// Update density
					snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;
				}
			}
			diff_lyr = fabs(pot_depth[lyr] - snowlayers[lyr].depth);
		}
	}

	// Check if snow layers should be erased
	if (nsnowlayers > pot_nsnowlayers) {
		for (int lyr = nsnowlayers - 1; lyr >= pot_nsnowlayers; lyr--) {
			if (!negligible(snowlayers[lyr].wmass + snowlayers[lyr].imass)) {
				snowlayers[pot_nsnowlayers - 1].update_temp(0.0, snowlayers[lyr].wmass + snowlayers[lyr].imass, Cp_water);
				snowlayers[pot_nsnowlayers - 1].wmass += snowlayers[lyr].wmass + snowlayers[lyr].imass;
			}
			snowlayers.erase(snowlayers.begin() + lyr);
		}
		nsnowlayers = (int)snowlayers.size();
	}

	double final_depth = 0.0;
	double final_ice = 0.0;
	double final_water = 0.0;
	double final_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_depth += snowlayers[lyr].depth;
		final_ice += snowlayers[lyr].imass;
		final_water += snowlayers[lyr].wmass;
		final_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	final_temp /= final_depth;

	assert(fabs(initial_depth - final_depth) < EPS);
	assert(fabs(initial_ice - final_ice) < EPS);
	assert(fabs(initial_water - final_water) < EPS);
//	assert(fabs(initial_temp - final_temp) < EPS);

	if (fabs(initial_temp - final_temp) >= EPS)
		int sch = 0;
}

/// Redistribute snow when all snow layers are full due to melt or compaction
void Snow::redist_layers_full() {

	const double EPS = 1.0e-10;

	double initial_depth = 0.0;
	double initial_ice = 0.0;
	double initial_water = 0.0;
	double initial_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_depth += snowlayers[lyr].depth;
		initial_ice += snowlayers[lyr].imass;
		initial_water += snowlayers[lyr].wmass;
		initial_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	initial_temp /= initial_depth;

	// filling each layer to snow_depth_standard
	for (int lyr = nsnowlayers - 1; lyr > 0; lyr--) {
		int lyr_move = lyr;

		// Move snow upwards in snowpack to fill every layer to snow_depth_standard
		double lyr_diff = snow_depth_standard - snowlayers[lyr].depth;
		
		while ((lyr_diff > EPS)) {
			lyr_move--;

			if (!negligible(snowlayers[lyr_move].imass)) {
				// Move upwards maximum of what is in snow layer lyr_below
				double depth_move = min(snowlayers[lyr_move].depth, lyr_diff);
				lyr_diff -= depth_move;

				double frac_move = !negligible(snowlayers[lyr_move].depth) ? depth_move / snowlayers[lyr_move].depth : 0.0;

				// Recalculate temperature
				snowlayers[lyr].update_temp(snowlayers[lyr_move].temp, frac_move * (snowlayers[lyr_move].imass + snowlayers[lyr_move].wmass), snowlayers[lyr_move].Ch);

				// Update depths
				snowlayers[lyr].depth += depth_move;
				snowlayers[lyr_move].depth -= depth_move;

				// Move ice accordingly
				snowlayers[lyr].imass += frac_move * snowlayers[lyr_move].imass;
				snowlayers[lyr_move].imass *= 1.0 - frac_move;

				// Move water accordingly
				snowlayers[lyr].wmass += frac_move * snowlayers[lyr_move].wmass;
				snowlayers[lyr_move].wmass *= 1.0 - frac_move;

				// Update density
				snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;
			}
		}
	}

	// Check for small depths
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		if (snowlayers[lyr].depth < 1.0e-12) {
			snowlayers[lyr].depth = 0.0;
			snowlayers[lyr].imass = 0.0;
			snowlayers[lyr].wmass = 0.0;
		}
	}

	double final_depth = 0.0;
	double final_ice = 0.0;
	double final_water = 0.0;
	double final_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_depth += snowlayers[lyr].depth;
		final_ice += snowlayers[lyr].imass;
		final_water += snowlayers[lyr].wmass;
		final_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	final_temp /= final_depth;

	assert(fabs(initial_depth - final_depth) < EPS);
	assert(fabs(initial_ice - final_ice) < EPS);
	assert(fabs(initial_water - final_water) < EPS);
	//assert(fabs(initial_temp - final_temp) < EPS);
}

// Update snow water holding capacity
void Snow::update_whc() {
	// Calculate layer water holding capacity, based on eqn 10 Wang et al. 2013

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].update_whc();
	}
}

/// Balance ice and water in snow layers after melting and refreezing
void Snow::balance_snow(double& rain_melt, double& rain_melt_temp) {

	// No snow layers or not multilayer snow scheme then leave function
	if (!nsnowlayers || !ifmultilayersnow)
		return;

	const double EPS = 1.0e-10;

	double initial_depth = 0.0;
	double initial_ice = 0.0;
	double initial_water = rain_melt;
	double initial_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_depth += snowlayers[lyr].depth;
		initial_ice += snowlayers[lyr].imass;
		initial_water += snowlayers[lyr].wmass;
		initial_temp += snowlayers[lyr].temp * (snowlayers[lyr].imass + snowlayers[lyr].wmass);
	}
	initial_temp /= (initial_ice + initial_water);

	// Re-distribute layers
	double max_snow = NSNOWLAYERS_MAX * snow_depth_standard;

	// Balance layers depending on if snowpack is over maxumum or not
	double snowpack_depth = depth();
	
	if (snowpack_depth > max_snow) {
		redist_layers_full();
	}
	else if (!negligible(snowpack_depth)) {
		redist_layers(snowpack_depth);
	}	
	else if (nsnowlayers) {
		// Insignificant snow depth. Remove all snow layers
		for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {
			double added_wmass = snowlayers[lyr].imass + snowlayers[lyr].wmass;
			if (!negligible(added_wmass)) {
				rain_melt_temp = (rain_melt_temp * rain_melt + snowlayers[lyr].temp * added_wmass) / (rain_melt + added_wmass);
				rain_melt += added_wmass;
			}

			snowlayers.erase(snowlayers.begin() + lyr);
		}
		nsnowlayers = 0;
	}

	// Update snow water holding capacity
	update_whc();

	// Save this timesteps extra melt water either from insignificant snow layers or that snow can't hold all water
	double ts_melt = 0.0;
	double ts_melt_temp = 0.0;

	// Redistribute liquid water within the snow pack
	// moving from top to bottom snow layers
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

	// Add excess water this timesteps melt to surface melt-water to pass to soil on next timestep
	if (!negligible(ts_melt)) {
		rain_melt_temp = (rain_melt_temp * rain_melt + ts_melt_temp * ts_melt) / (rain_melt + ts_melt);
		rain_melt += ts_melt;
	}

	double final_depth = 0.0;
	double final_ice = 0.0;
	double final_water = rain_melt;
	double final_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_depth += snowlayers[lyr].depth;
		final_ice += snowlayers[lyr].imass;
		final_water += snowlayers[lyr].wmass;
		final_temp += snowlayers[lyr].temp * (snowlayers[lyr].imass + snowlayers[lyr].wmass);
	}
	final_temp /= (final_ice + final_water);

	assert(fabs(initial_depth - final_depth) < EPS);
	assert(fabs(initial_ice - final_ice) < EPS);
	//assert(fabs(initial_water - final_water) < EPS);
	//assert(fabs(initial_temp - final_temp) < EPS);

	//if (fabs(initial_temp - final_temp) >= EPS) // DW_COMMENT error the whole time!!!!!!!!!!!
	//	dprintf("Year %d Day %d ERROR: Snow temp difference %g\n",date.year, date.day, fabs(initial_temp - final_temp));
}

void Snow::phase_change() {

	// DESCRIPTION
	// Calculate water, ice and air fraction in each layer
	// and snow phase change

	double Ffreez;		// Fraction of water which freezes
	double Fthaw;		// Fraction of ice which thaws
	double energy;		// Energy per volume [J m-3]
	double delta_T;		// Temperature difference (current - previous)
	double snow_depth;	// Snow layer depth [mm]
	double snow_temp;	// Snow layer temperature [C]
	double melt_factor;	// factor governing the rate of snowmelt (min depth allowed)
	melt_factor = 1.0; // DW_COMMENT should check what is a ok value prefereable 1
	const double EPS = 1.0e-10;

	if (!ifmultilayersnow)
		return;

	// Minimum air volume that will influence snow layer (only affect top snow layer)
	double min_snow_depth = snow_depth_standard * (snow_depth_limit - 1.0) * melt_factor;

	// Going through snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		bool iffreezed = false;
		bool ifthawed = false;

		// The total amount of ice and water should not change in this routine
		double wtotal_initial = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		snow_depth = snowlayers[lyr].depth;
		snow_temp = snowlayers[lyr].temp;

		// Water and ice fraction of full volume
		double Fwmass = snowlayers[lyr].wmass / snow_depth;
		double Fimass = snowlayers[lyr].imass / snow_depth;

		// Volume of snow layer [m-3]
	//	snow_volume = snowlayers[lyr].depth * M_PER_MM;

		// Temperature change from yesterday
		delta_T = snow_temp - snowlayers[lyr].temp_old;

		// Freezing
		if (snow_temp < 0.0 && snow_temp < snowlayers[lyr].temp_old) { // Wania conditions 

			iffreezed = true;
			Fthaw = 0.0;

			// Water to freeze
			if (!negligible(Fwmass)) {

				// Energy to freeze water [J m-3]
				//energy = snowlayers[lyr].Ch * fabs(delta_T);
				energy = Cp_water * fabs(delta_T);

				// Fraction water that would freeze [m3 m-3]
				Ffreez = energy / Lheat;

				// More water in a volume than could freeze
				if (Fwmass > Ffreez) {

					// There is some water left, even after freeze
					Fimass += Ffreez;
					Fwmass -= Ffreez;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
				}
				else {

					// It is cold enough to freeze all the water and
					// decrease the temperature below 0, so only the 
					// energy need for phase change is removed from 
					// the already calculated temperature 

					// All water freezes
					Fimass += Fwmass;

					// Update snow temperature in this layer
					// How much energy (temperature) the phase change needed
					snowlayers[lyr].temp += Fwmass * Lheat / Cp_water;

					if (snowlayers[lyr].temp > 0.0)
						dprintf("ERROR: Year %d Day %d Snow layer %d temperature is above zero AFTER freezing %g\n",date.year,date.day,lyr, snowlayers[lyr].temp);

					snowlayers[lyr].temp = min(0.0, snowlayers[lyr].temp);

					// No liquid water left in snow layer
					Fwmass = 0.0;
				}
			}
			else {
				Ffreez = 0.0;
				Fimass += Fwmass;
				Fwmass = 0.0;
			}
		} // end of freezing

		// Thawing  
		else if (snow_temp >= 0.0  && snow_temp > snowlayers[lyr].temp_old) { // Wania conditions

			ifthawed = true;
			Ffreez = 0.0;

			// Ice to thaw
			if (!negligible(Fimass)) {

				// Energy to thaw ice [J m-3]
				//energy = snowlayers[lyr].Ch  * fabs(delta_T);
				energy = Cp_water  * fabs(delta_T);

				// Fraction ice that would thaw [m3 m-3]
				Fthaw = energy / Lheat;

				// More ice in a volume than could thaw
				if (Fimass > Fthaw) {

					// There is some ice left, even after melt
					Fwmass += Fthaw;
					Fimass -= Fthaw;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
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
			}
			else {
				Fthaw = 0.0;
				Fwmass += Fimass;
				Fimass = 0.0;
			}
		} // end of thawing

		  // New ice fraction in snow layer
		snowlayers[lyr].imass = Fimass * snow_depth;

		// New water fraction in snow layer
		snowlayers[lyr].wmass = Fwmass * snow_depth;

		// Freezing - Update density as depth doesn't change
		if (iffreezed) {
			if (!negligible(snowlayers[lyr].depth)) {
				snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;
			}
		}
		// Thawing - Update depth as density doesn't change
		else if (ifthawed) {
			snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);
		}

		// The total amount of ice and water should not change in this routine
		double wtotal_final = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		assert(fabs(wtotal_initial - wtotal_final) < EPS);

	} // end of loop through layers
}


/// Snow compaction and update density and depth
void Snow::compaction(double& rain_melt, double& rain_melt_temp) {

	if (!ifmultilayersnow)
		return;

	// Mechanical compaction scheme, following Best et al. 2011
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// Current layers ice and water mass
		double I_sum = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		// Snow mass above this layer
		double I_above = 0.0;
		for (int ii = lyr + 1; ii < nsnowlayers; ii++) {
			I_above += snowlayers[ii - 1].imass + snowlayers[ii - 1].wmass;
		}

		// Cumulative snow mass affecting compaction
		double M_k = 0.5 * I_sum + I_above;

		// Layer temperature in Kelvin (never above zero degree Celcius)
		double T_k = min(K2degC, snowlayers[lyr].temp + K2degC);

		// Increase in layer density due to compaction (Best et. al. 2011, Eq.21)
		double k_s = 4000.0;
		double eta_0 = 1E7; // Compactive viscosity factor for snow
		double rho_0 = 50.0; // Reference snow density [kg m-3]

		if (snowlayers[lyr].dens < rho_ICE) {
			double delta_dens = snowlayers[lyr].dens * gravity * M_k / eta_0 *
				exp((k_s / K2degC) - (k_s / T_k) - (snowlayers[lyr].dens / rho_0)) * SECS_PER_DAY;
			snowlayers[lyr].dens += delta_dens;

			// Snow density can not exceed maximum density
			snowlayers[lyr].dens = min(rho_ICE, snowlayers[lyr].dens);
		}
	}

	// Calculate snow depth after compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);
	}

	// Redistribute snow layers, update water holding capacity and calculate water that should be added to rain_melt
	balance_snow(rain_melt, rain_melt_temp);
}

/// Accumulate snow in bottom layer if all other layers are already full
void Snow::update_layers(double snowfall_depth, double snowfall_dens, double airT) {

	// Check initial snow and water content
	const double EPS = 1.0e-10;

	double initial_depth = snowfall_depth;
	double initial_ice = snowfall_depth * snowfall_dens / rho_H2O;
	double initial_water = 0.0;
	double initial_temp = airT * snowfall_depth;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_depth += snowlayers[lyr].depth;
		initial_ice += snowlayers[lyr].imass;
		initial_water += snowlayers[lyr].wmass;
		initial_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	initial_temp /= initial_depth;

	// Loop through snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// Total snow to move down to this layer
		double tot_depth_move = lyr == 0 ? snowfall_depth : snow_depth_standard - snowlayers[lyr].depth; // DW_COMMENT always the same for all layers?

		// Layer to move snow from
		int lyr_move = lyr + 1;

		// Move snow until enough is moved (could be from more than one layer above)
		while (!negligible(tot_depth_move)) {

			double depth_move;

			// Determine snow to move from snow layer above

			// Snow layer above exist
			if (lyr_move < nsnowlayers) {

				// Move maximum what exist in layer to move from
				depth_move = min(tot_depth_move, snowlayers[lyr_move].depth);
				double frac_move = !negligible(snowlayers[lyr_move].depth) ? depth_move / snowlayers[lyr_move].depth : 0.0;

				// Recalculate temperature
				snowlayers[lyr].update_temp(snowlayers[lyr_move].temp, frac_move * (snowlayers[lyr_move].imass + snowlayers[lyr_move].wmass), snowlayers[lyr_move].Ch);

				// Update depths
				snowlayers[lyr].depth += depth_move;
				snowlayers[lyr_move].depth -= depth_move;

				// Move ice accordingly
				snowlayers[lyr].imass += frac_move * snowlayers[lyr_move].imass;
				snowlayers[lyr_move].imass *= 1.0 - frac_move;

				// Move water accordingly
				snowlayers[lyr].wmass += frac_move * snowlayers[lyr_move].wmass;
				snowlayers[lyr_move].wmass *= 1.0 - frac_move;

				// Update density
				snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;
			}
			else { // No layer above so take from snowfall to fill layer
				depth_move = min(tot_depth_move, snow_depth_standard - snowlayers[lyr].depth);

				// Recalculate temperature
				snowlayers[lyr].update_temp(airT, depth_move * snowfall_dens / rho_H2O, snowheatcap(airT, snowfall_dens));

				// Update depths
				snowlayers[lyr].depth += depth_move;

				// Move ice accordingly
				snowlayers[lyr].imass += depth_move * snowfall_dens / rho_H2O;

				// Update density
				snowlayers[lyr].dens = snowlayers[lyr].imass * rho_H2O / snowlayers[lyr].depth;
			}

			// Remove added snow from this layer from the toal that needs to be moved
			tot_depth_move -= depth_move;

			// Go to the next layer above (if tot_snow_move > 0)
			lyr_move++;
		}
	}

	// Check final snow and water content
	double final_depth = 0.0;
	double final_ice = 0.0;
	double final_water = 0.0;
	double final_temp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_depth += snowlayers[lyr].depth;
		final_ice += snowlayers[lyr].imass;
		final_water += snowlayers[lyr].wmass;
		final_temp += snowlayers[lyr].temp * snowlayers[lyr].depth;
	}
	final_temp /= final_depth;

	assert(fabs(initial_depth - final_depth) < EPS);
	assert(fabs(initial_ice - final_ice) < EPS);
	assert(fabs(initial_water - final_water) < EPS);
	//assert(fabs(initial_temp - final_temp) < EPS);
}

///
/*
*
*/
void Snow::snowfall(double airT, double& prec, double u10) {

	// Fresh snow density calculation, following eqn 1 in Vionnet et al (2012)
	int a = 109; // kg m-3, 
	int b = 6;  // kg m-3 K-1
	int c = 26;  // kg m-7/2 s -1/2
	double dens_min = 100; // minimum snow density 
	double Tmelt = 0; // temperature of the melting point for water

	double snowfall_dens = max(dens_min, a + b * (airT - Tmelt) + c * pow(u10, 0.5));
	double snowfall_Ch = snowheatcap(airT, snowfall_dens);

	// Calculate the depth of the newly fallen snow based on precipitation & density (Oleason 2013, Eq. 6.88)
	double snowfall_depth = prec / (snowfall_dens / rho_H2O);

	/////////////////////////////////////////////////////////
	// DW_COMMENT remove
	// Check initial depth, snow, water and temperature
	double EPS = 1.0e-12;
	double initial_depth = snowfall_depth;
	double initial_ice = prec;
	double initial_water = 0.0;
	if (nsnowlayers) {
		for (int lyr = 0; lyr < nsnowlayers; lyr++) {
			initial_depth += snowlayers[lyr].depth;
			initial_ice += snowlayers[lyr].imass;
			initial_water += snowlayers[lyr].wmass;
		}
	}

	/////////////////////////////////////////////////////////

	prec = 0.0;

	// Transferred snow depth and density from another layer
	double transfer_depth = 0.0;
	double transfer_dens = 0.0;
	double transfer_imass = 0.0;
	double transfer_wmass = 0.0;
	double transfer_temp = 0.0;
	double transfer_Ch = 0.0;

	// If no snow layers
	if (nsnowlayers == 0) {
		snowlayers.resize(nsnowlayers + 1);
		snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT); // DW_COMMENT new and soilcolumn added to constructor
		nsnowlayers++;
	}

	int top_snow_lyr = nsnowlayers - 1;

	// Fill snow layers from top_snow_lyr
	while (!negligible(transfer_depth + snowfall_depth) && top_snow_lyr < NSNOWLAYERS_MAX) {

		double snow_depth_max;
		// If not the topmost layer and total layer depth is lower than max depth limit after snowfall and transfer of snow then fill layer to depth limit
		if (top_snow_lyr < NSNOWLAYERS_MAX - 1 && (snowlayers[top_snow_lyr].depth + transfer_depth + snowfall_depth) <= (snow_depth_standard * snow_depth_limit)) {
			snow_depth_max = snow_depth_standard * snow_depth_limit;
		}
		else { // else add snow to standard depth and add additional snow to a new layer or start moving snow downwards in snowpack if snowpack is full
			snow_depth_max = snow_depth_standard;
		}

		// Snow layer depth exceeds maximum snow depth -> snow needs to be transferred to layer above
		if (snow_depth_max < snowlayers[top_snow_lyr].depth) {
			transfer_depth = snowlayers[top_snow_lyr].depth - snow_depth_max;
			transfer_dens = snowlayers[top_snow_lyr].dens;
			transfer_wmass = snowlayers[top_snow_lyr].wmass * (transfer_depth / snowlayers[top_snow_lyr].depth);
			transfer_imass = snowlayers[top_snow_lyr].imass * (transfer_depth / snowlayers[top_snow_lyr].depth);
			transfer_temp = snowlayers[top_snow_lyr].temp;
			transfer_Ch = snowlayers[top_snow_lyr].Ch;

			snowlayers[top_snow_lyr].depth = snow_depth_max;
			snowlayers[top_snow_lyr].wmass -= transfer_wmass;
			snowlayers[top_snow_lyr].imass -= transfer_imass;
		}
		// else add snow to this snow layer
		else {

			// Add transferred snow to a new empty layer
			if (!negligible(transfer_depth)) {

				if (snowlayers[top_snow_lyr].depth > 0.0 || snowlayers[top_snow_lyr].imass > 0.0)
					dprintf("ERROR: Year %d Day %d New snowlayer isn't empty: Depth %g Ice %g\n", date.year, date.day, snowlayers[top_snow_lyr].depth, snowlayers[top_snow_lyr].imass);

				snowlayers[top_snow_lyr].depth = transfer_depth;
				snowlayers[top_snow_lyr].temp = transfer_temp;
				snowlayers[top_snow_lyr].imass = transfer_imass;
				snowlayers[top_snow_lyr].wmass = transfer_wmass;
				snowlayers[top_snow_lyr].Ch = transfer_Ch;

				// Reset transferred snow
				transfer_depth = 0.0;
				transfer_dens = 0.0;
				transfer_imass = 0.0;
				transfer_wmass = 0.0;
				transfer_temp = 0.0;
				transfer_Ch = 0.0;
			}

			// Add snowfall
			double snowdepthtoadd = min(snowfall_depth, snow_depth_max - snowlayers[top_snow_lyr].depth);
			double snowfall_imass = snowdepthtoadd * (snowfall_dens / rho_H2O);

			// Update snow temp
			if (!negligible(snowfall_imass)) {
				snowlayers[top_snow_lyr].update_temp(airT, snowfall_imass, snowfall_Ch);
			}

			snowlayers[top_snow_lyr].depth += snowdepthtoadd;
			snowlayers[top_snow_lyr].imass += snowfall_imass;
			snowfall_depth -= snowdepthtoadd;

			// Update snow layer density
			snowlayers[top_snow_lyr].dens = snowlayers[top_snow_lyr].imass * rho_H2O / snowlayers[top_snow_lyr].depth;
		}

		// Go to next layer above
		top_snow_lyr++;

		// If top_snow_lyr exceeds the number of snow layers and there is additional snow to be added
		// then either all snow layers are full up to standard depths and snow needs to be reshuffled
		// or a new snow layer needs to be created
		if (snowfall_depth > 0.0) {
			if (top_snow_lyr == NSNOWLAYERS_MAX) {

				if (transfer_depth > 0.0)
					dprintf("ERROR: Year %d Day %d Transfer_Snow_depth not zero: %g\n", date.year, date.day, transfer_depth);

				while (snowfall_depth > 0.0) {

					double snow_depth_to_add = min(snowfall_depth, snow_depth_standard);
					update_layers(snow_depth_to_add, snowfall_dens, airT);  // DW_COMMENT test by adding insane amount of snowfall 
					snowfall_depth -= snow_depth_to_add;
				}
			}
			else if (top_snow_lyr == nsnowlayers) {
				snowlayers.resize(nsnowlayers + 1);
				snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT); // DW_COMMENT new and soilcolumn added to constructor
				nsnowlayers++;
			}
		}
	}

	/////////////////////////////////////////////////////////
	// DW_COMMENT remove
	double final_depth = 0.0;
	double final_ice = 0.0;
	double final_water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_depth += snowlayers[lyr].depth;
		final_ice += snowlayers[lyr].imass;
		final_water += snowlayers[lyr].wmass;
	}

	assert(fabs(initial_depth - final_depth) < EPS);
	assert(fabs(initial_ice - final_ice) < EPS);
	assert(fabs(initial_water - final_water) < EPS);
	/////////////////////////////////////////////////////////
}

/// Handld ROS events and rainwater distribution within the snowpack
/* Calculate water, ice and air fraction in each layer
* and snow phase changein layer based on the available energy
*/
void Snow::rainwater_melt(double airT, double& prec, double& rain_melt, double& rain_melt_temp) {

	// water balance check
	double w_total_initial = prec;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		w_total_initial += snowlayers[lyr].wmass + snowlayers[lyr].imass;
	}

	// Water transferred between snow layer and in the end to the soil
	double water_transferred = prec;

	// temperature difference between air temperture and zero
	double delta_T = airT - 0.0;
	
	// total energy from rain [J]
	double rain_energy = Cp_water * delta_T * prec * M_PER_MM;

	// top snow lyr index
	int lyr = nsnowlayers - 1;

	// loop through layers as long as water still has a temperture above zero (still energy to melt ice) 
	while (rain_energy > 0.0 && !negligible(snowlayers[lyr].imass)) {

		// Fraction ice in snow layer
		double Fice = snowlayers[lyr].imass / snowlayers[lyr].depth;
		// Fraction water in snow layer
		double Fwater = snowlayers[lyr].wmass / snowlayers[lyr].depth;
		// Heat capcity of snow layer [J m-3 K-1]
		double heat_cap = snowlayers[lyr].Ch;
		// Volume of snow layer [m-3]
		double snow_volume = snowlayers[lyr].depth * M_PER_MM;
		// Eenergy needed to bring snow layer to 0C [J]
		double energy_to_zero = fabs(snowlayers[lyr].temp - 0.0) * heat_cap / snow_volume;

		// Enough energy to melt some snow
		if (rain_energy > energy_to_zero) {

			// Snow layer temp is set to zero and erase energy needed to bring up temperature
			rain_energy -= energy_to_zero;
			snowlayers[lyr].temp = 0.0;

			// Fraction that will melt
			double frac_change = rain_energy / Lheat;

			// Enough energy to melt all ice (snow layer will be removed) 
			if (frac_change > Fice) {
				frac_change -= Fice;
				Fwater += Fice;
				Fice = 0.0;
				rain_energy = frac_change * Lheat;
			}
			// only melt a frac of ice
			else {
				Fwater += frac_change;
				Fice -= frac_change;
				rain_energy = 0.0;
			}
		}
		// only enough energy to change temperature
		else {
			snowlayers[lyr].temp += rain_energy * snow_volume / heat_cap;
			rain_energy = 0.0;
		}

		// Convert from fraction to mm depth 
		snowlayers[lyr].wmass = Fwater * snowlayers[lyr].depth;
		snowlayers[lyr].imass = Fice * snowlayers[lyr].depth;

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
			snowlayers.erase(snowlayers.begin());
			nsnowlayers--;

			// If this was the last snow layer, then exit loop
			if (nsnowlayers == 0)
				break;
		}

		// go thorugh lyrs while rain_energy not 0
		lyr--;
	} // end of while loop

	// If water to be transferred to the soil
	if (!negligible(water_transferred)) {

		// Still energy left, update rain_melt_temp with melted water 
		double water_transferred_temp = 0.0;
		if (!negligible(rain_energy)) {
			water_transferred_temp = rain_energy / (Cp_water * water_transferred * M_PER_MM);
		}
		rain_melt_temp = (rain_melt_temp * rain_melt + water_transferred_temp * water_transferred) / (rain_melt + water_transferred);

		// Save melt water above snow whc which will be pass downwards
		rain_melt += water_transferred;
	}

	// Check final water content
	double w_total_final = water_transferred;
	for (int lyr = 0; lyr <  nsnowlayers; lyr++) {
		w_total_final += snowlayers[lyr].wmass + snowlayers[lyr].imass;
	}
	const double EPS = 1.0e-10;
	assert(fabs(w_total_initial - w_total_final) < EPS);
}

/// Snow accumulation
void Snow::snow_or_rainfall(double airT, double& prec, double u10, double& rain_melt, double& rain_melt_temp) {

	// Multi layer snow scheme
	if (ifmultilayersnow) {
		
		// Snowfall
		if (airT < 0.0 && prec > 0.0) {

			// Add precipitation as snow and update depth, density
			snowfall(airT, prec, u10);
		}
		// Precipitation falls on snow
		else if (prec > 0.0 && airT > 0.0 && !negligible(depth())) {

			// Phase changes and temperature change due to rainwater
			rainwater_melt(airT, prec, rain_melt, rain_melt_temp);
		}
		else if (prec > 0.0) {

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
		if (airT < 0.0 && prec > 0.0) {

			// If no snow layers
			if (nsnowlayers == 0) {
				snowlayers.resize(nsnowlayers + 1);
				snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT); // DW_COMMENT new and soilcolumn added to constructor
				snowlayers[0].dens = SNOWDENSITY;
				nsnowlayers++;
			}

			// Snowfall up to max snowpack size
			double snowfall = min(prec, SNOWPACK_MAX - snowlayers[0].imass);

			snowlayers[0].imass += snowfall;

			rain_melt += max(0.0, prec - snowfall);
		}
		// Snow melting today
		else if (airT > 0.0 && nsnowlayers) {

			// Snow melt formulation (Choudhury et al 1998)
			double snowmelt = min((1.5 + 0.007 * prec) * (airT - 0.0), snowlayers[0].imass);

			snowlayers[0].imass -= snowmelt;

			rain_melt += prec + snowmelt;
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
}

void Snow::serialize(ArchiveStream& arch) {
	arch & snowlayers
		& nsnowlayers
		& snow_depth_standard
		& snow_depth_limit
		& NH4_mass
		& NO3_mass;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soillayer member functions
////////////////////////////////////////////////////////////////////////////////

//Soillayer::Soillayer(Soilcolumn* sc) : soilcolumn(sc) {
//}

Soillayer::Soillayer(Soilcolumn* sc, double d, double dp, bool a, bool c) : soilcolumn(sc), depth(d), depth_profile(dp), acrotelm(a), catotelm(c) {
	init();
}


void Soillayer::init() {

	// Initialises certain member variables

	// Energy
	temp = 0.0;
	Kh = 0.0;
	Ch = 0.0;

	// Water
	wmass = 0.0;
	imass = 0.0;
	por = 0.0;
	Dw = 0.0;
	Kw = 0.0;
	perc_mass = 0.0;
	runoff_drain = 0.0;
	runoff_baseflow = 0.0;

	// Leaching
	orgleachfrac = 0.0;
	aorgCleach = 0.0;
	aorgNleach = 0.0;
	aminleach = 0.0;
	orgmovefrac = 0.0;

	// Old SOM module
	alag = 0.0;
	exp_alag = 1.0;
	cpool_slow = 0.0;
	cpool_fast = 0.0;
	decomp_litter_mean = 0.0;
	k_soilfast_mean = 0.0;
	k_soilslow_mean = 0.0;

	/////////////////////////////////////////////////////
	// Initialise CENTURY pools

	// Set initial CENTURY pool N:C ratios
	// Parton et al 1993, Fig 4

	sompool[SOILMICRO].ntoc = 1.0 / 15.0;
	sompool[SURFHUMUS].ntoc = 1.0 / 15.0;
	sompool[SLOWSOM].ntoc = 1.0 / 20.0;
	sompool[SURFMICRO].ntoc = 1.0 / 20.0;

	// passive has a fixed value
	sompool[PASSIVESOM].ntoc = 1.0 / 9.0;

	// Decay rate modifier depending on depth profile (Koven et al 2013)
	// reference depth 10m
	depth_mod = exp(-depth_profile / 10000.0);

	// N pools
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	NO2_mass = 0.0;
	NO_mass = 0.0;
	N2O_mass = 0.0;
	N2_mass = 0.0;
	anmin = 0.0;
	animmob = 0.0;
	aminleach = 0.0;
	aorgNleach = 0.0;
	aorgCleach = 0.0;

	// Gases
	CO2_soil = 0.0;
	CO2_soil_yesterday = 0.0;
	CO2_soil_prod = 0.0;
	CH4 = 0.0;
	CH4_yesterday = 0.0;
	CH4_prod = 0.0;
	CH4_oxid = 0.0;
	CH4_ebull_ind = 0.0;
	CH4_ebull_vol = 0.0;
	CH4_gas = 0.0;
	CH4_diss = 0.0;
	CH4_gas_yesterday = 0.0;
	CH4_diss_yesterday = 0.0;
	CH4_gas_vol = 0.0;
	CH4_vgc = 0.0;
	O2 = 0.0; // DW_COMMENT Same as air above

	// tiller area
	tiller_area = 0.0001;

	// Daily heterotrophic respiration
	dcflux_soil = 0.0;

	// soil labile carbon availability
	labile_carbon = 0.0;

	// Monthly variables
	for (int mth = 0; mth < 12; mth++) {
		mwcont[mth] = 0.0;
		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		morgmove_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
	}

	// Daily variables
	std::fill_n(dwcont, Date::MAX_YEAR_LENGTH, 0.0);
}


// Set soil layer properties
/* Set soil layer properties
 * Sets the following properties for this individual layer
 * clay, silt and sand fraction
 * pH
 * hydrologic parameters (wsats, awc, wp and por)
 * mineral, organic and peat fractions
 * thermal properties for analytic temperature calculation
 * and two layer soil perculation parameters
*/
void Soillayer::set_properties(double& wmass_move, double& imass_move) {

	// Get soil type
	Soiltype soiltype = (*soilcolumn).soiltype;

	// Only enter first year except if dynamic soil organic C influence soil properties (iforganicsoilproperties == 2)
	if (!((date.day == 0 && (*soilcolumn).patch.stand.first_year == date.year) || (iforganicsoilproperties == 2 && !acrotelm && !catotelm))) // DW_COMMENT restart years needs to update it then as well!
		return;

	// Analytic temp calc parameters
	thermdiff_0 = soiltype.thermdiff_0;
	thermdiff_15 = soiltype.thermdiff_15;
	thermdiff_100 = soiltype.thermdiff_100;

	if (acrotelm) {
		wsats = acrotelm_por * depth;
		awc = acrotelm_por * depth;
		wp = peat_wp * depth;
		por = acrotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - acrotelm_por;

		// Initially saturated
		double wadded = wsats - (wmass + imass);
		wmass += wadded;
		(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
	}
	else if (catotelm) {
		wsats = catotelm_por * depth;
		awc = catotelm_por * depth;
		wp = peat_wp * depth;
		por = catotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - catotelm_por;

		// Initially saturated
		double wadded = wsats - (wmass + imass);
		wmass += wadded;
		(*soilcolumn).patch.fluxes.report_flux(Fluxes::WADDED, wadded);
	}
	// Mineral soil and no influence of soil organic C
	else if (!iforganicsoilproperties) {

		wsats = soiltype.wsats * depth;
		awc = soiltype.awc * depth;
		wp = soiltype.wp * depth;
		por = soiltype.porosity * depth;
		organic_frac = soiltype.organic_frac;
		mineral_frac = soiltype.mineral_frac;
		peat_frac = 0.0;

		// Two layer hydrology parameter
		perc_base = soiltype.perc_base;

		dprintf("Year %d Day %d set_properties depth %g org_frac %g min_frec %g por_frac %g WSATS %g AWC %g WP %g\n", 
			date.year, date.day, depth_profile, organic_frac, mineral_frac, por / depth, wsats / depth, awc / depth, wp / depth);

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

		// Maximum soil layer C content (kg C m-2)
		double soil_c_max = maxSOCdensity * (depth / MM_PER_M);

		if (soil_c > soil_c_max)
			dprintf("Year %d Day %d Depth %g Soil C %g MAX %g\n",date.year,date.day,depth_profile,soil_c,soil_c_max); // DW_COMMENT move C downwards!!!!!!!!!!

		// Soil C density in this soil layer (kg C m-3) 
		double soil_density = soil_c / (depth / MM_PER_M);

		double mineral_porosity = 0.489 - 0.00126 * 100.0 * soiltype.sand_frac;	// Eqn 2 (Lawrence & Slater 2008)

		// Calculate organic fraction for soil layer using the Lawrence and Slater (2008) approach
		// This is the fraction of the material (i.e. 1 - porosity) that is organic
		double org_frac = min(1.0, soil_density / maxSOCdensity);
		double min_frac = 1.0 - org_frac;
		
		double porosity = min_frac * mineral_porosity + org_frac * organic_porosity;

		wsats = (min_frac * soiltype.wsats + org_frac * soiltype.wsats_org) * depth;
		awc = (min_frac * soiltype.awc + org_frac * soiltype.awc_org) * depth;
		wp = (min_frac * soiltype.wp + org_frac * soiltype.wp_org) * depth;
		por = porosity * depth;
		organic_frac = (1.0 - porosity) * org_frac;
		mineral_frac = (1.0 - porosity) * min_frac;
		peat_frac = 0.0;

		// Move water downwards if water content exceeds wsats
		if (wsats < wmass + imass) {
			double wmass_frac = wmass / (wmass + imass);
			double move = wmass + imass - wsats;
			wmass_move = wmass_frac * move;
			imass_move = (1 - wmass_frac) * move;
			wmass -= wmass_move;
			imass -= imass_move;
		}

		// Analytic temp calc parameters
		thermdiff_0 = min_frac * soiltype.thermdiff_0 + org_frac * soiltype.thermdiff_0_org;
		thermdiff_15 = min_frac * soiltype.thermdiff_15 + org_frac * soiltype.thermdiff_15_org;
		thermdiff_100 = min_frac * soiltype.thermdiff_100 + org_frac * soiltype.thermdiff_100_org;

		// Two layer hydrology parameter
		perc_base = min_frac * soiltype.perc_base + org_frac * soiltype.perc_base_org;

		if (date.year == 1 || date.year == 100 || date.year == 200 || date.year == 300 || date.year == 400 || date.year == 500 || date.year == 600)
			dprintf("Year %d Day %d set_properties depth %g org_frac %g min_frac %g por_frac %g WSATS %g AWC %g WP %g\n",
				date.year, date.day, depth_profile, organic_frac, mineral_frac, por / depth, wsats / depth, awc / depth, wp / depth);
	}

	// If two layer soil scheme then wmass is never below wp
	if (iftwolayersoil) {
		wmass = max(wmass, wp);
	}
}

/// C content of soil layer
double Soillayer::ccont(bool exclude_surf) {
	double ccont = 0.0;

	ccont += cpool_fast;
	ccont += cpool_slow;

	for (int i = 0; i < NSOMPOOL; i++) {
		if ( !(exclude_surf && (i == SURFMETA || i == SURFSTRUCT || i == SURFFWD || i == SURFCWD)) )
			ccont += sompool[i].cmass;
	}
	return ccont;
}
 
/// N content of soil layer
double Soillayer::ncont() {
	double ncont = NH4_mass + NO3_mass + NO2_mass + NO_mass + N2O_mass + N2_mass;

	for (int i = 0; i < NSOMPOOL; i++)
		ncont += sompool[i].nmass;

	return ncont;
}

/// Water content of soil layer
double Soillayer::wcont() {

	return max(0.0, wmass - wp) / awc;
}

/// Field capacity of soil layer
double Soillayer::fc() {

	return awc + wp;
}

double Soillayer::fracair() {
	return max(0.0, (por - (wmass + imass)) / depth);
}


double Soillayer::nmass_avail(n_pref_type pref, bool wmass_limited) {

	double nmass = 0.0;
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
	return nmass * (wmass_limited ? min(1.0, wcont()) : 1.0);
}

void Soillayer::nmass_subtract(double nmass, n_pref_type pref) {

	if (!ifntransform) {
		pref = NH4;
	}
	double residual = 0.0;
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		if (nmass_tot > 0.0) {
			residual = nmass_tot - nmass;
			if (residual >= 0.0 || negligible(residual, -10)) {
				nmass_subtract(nmass * NH4_mass / nmass_tot, NH4);
				nmass_subtract(nmass * NO3_mass / nmass_tot, NO3);
			}
			else {
				fail("tried to subtract more N than available");
			}
		}
	}
	else if (pref == NH4) {
		if (nmass<NH4_mass || negligible(nmass - NH4_mass, -10)) {
			NH4_mass -= nmass;
		}
		else {
			double r = nmass - NH4_mass;
			NH4_mass = 0.0;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::NH3_SOIL, r);
			dprintf("Year %d Day %d NH4 mass %f N substr %f Diff %f\n", date.year ,date.day, NH4_mass, nmass, r);
		}
	}
	else if (pref == NO3) {
		if (nmass<NO3_mass || negligible(nmass - NO3_mass, -10)) {
			NO3_mass -= nmass;
		}
		else {
			double r = nmass - NO3_mass;
			NO3_mass = 0.0;
			(*soilcolumn).patch.fluxes.report_flux(Fluxes::NO_SOIL, r);
			dprintf("Year %d Day %d NO3 mass %f N substr %f Diff %f\n", date.year, date.day, NO3_mass, nmass, r);
		}
	}
}

void Soillayer::nmass_inc(double nmass, n_pref_type pref) {
	if (nmass < 0.0) {
		nmass_subtract(fabs(nmass), pref);
	}
	else {
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

void Soillayer::nmass_multiplic_inc(double inc, n_pref_type pref) {
	if (!ifntransform) {
		pref = NH4;
	}
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		NH4_mass *= inc * NH4_mass / nmass_tot;
		NO3_mass *= inc * NO3_mass / nmass_tot;
	}
	else if (pref == NH4) {
		NH4_mass *= inc;
	}
	else if (pref == NO3) {
		NO3_mass *= inc;
	}
}

void Soillayer::nmass_remove_small_values(n_pref_type pref) {
	double EPS = -1.0e-18;

	if (!ifntransform) {
		pref = NH4;
	}
	if (pref == NO) {
		if (NH4_mass < 0.0 && NH4_mass > EPS)
			NH4_mass = 0.0;
		if (NO3_mass < 0.0 && NO3_mass > EPS)
			NO3_mass = 0.0;
	}
	else if (pref == NH4) {
		if (NH4_mass < 0.0 && NH4_mass > EPS)
			NH4_mass = 0.0;
	}
	else if (pref == NO3) {
		if (NO3_mass < 0.0 && NO3_mass > EPS)
			NO3_mass = 0.0;
	}
}

void Soillayer::update_temp(double add_temp, double add_mass, double add_Ch) {

	if (negligible((imass + wmass) * Ch + add_mass * add_Ch))
		dprintf("Year %d Day %d ERRRRRRRRRROOOORRRRRRRRR \n", date.year, date.day);

	temp = (temp * (imass + wmass) * Ch + add_temp * add_mass * add_Ch) / ((imass + wmass) * Ch + add_mass * add_Ch);
}

void Soillayer::update_moisture_diff_and_cond() {

	// Wösten 1999 table 5 has eqns to determin wsats, Kw_max, alpha, lambda, n from clay, sand, silt, org, bd

	// van Genuchten soil parameters
	Soiltype soiltype = (*soilcolumn).soiltype;
	const double Kw_max = soiltype.Kw_max;
	const double alpha = soiltype.alpha;
	const double lambda = soiltype.l;
	const double m = soiltype.m; 

	double W_res = wp / depth;

	double W_sat, W_tot, I_frac;

	// Account for ice filling up space
	W_sat = wsats / depth;
	W_tot = wmass / depth;
	I_frac = min(1.0, imass / wsats);

	// Moisture limits to be considered in D and K calculations
	double theta;
	if (W_tot <= 1.001 * W_res) {
		theta = 1.001 * W_res;
		Kw = 0.0;
	}
	else if (W_tot > 1.001 * W_res && W_tot <= 0.999 * W_sat) {
		theta = W_tot;
	}
	else {
		theta = W_sat * 0.999;
		Kw = Kw_max;
	}

	double the = (theta - W_res) / (W_sat - W_res);
	double po = pow(the, (1 / m));

	// Diffusivity (m2 s-1; van Genuchten 1980 eqn 11, adopted from HTESSEL code)
	Dw = ((1.0 - m) * Kw_max / (alpha * m * (W_sat - W_res))) * pow(the, (lambda - 1.0 / m)) * (pow(1.0 - po, m) + pow(1.0 - po, -m) - 2); // DW_COMMENT lookup table to speed up code?

	// Conductivity (m s-1; van Genuchten 1980 eqn 8, adopted from HTESSEL code)
	if (W_tot > 1.001 * W_res && W_tot <= 0.999 * W_sat) {
		Kw = Kw_max * pow(the, lambda) * pow(1.0 - pow(1.0 - pow(the, 1.0 / m), m), 2); // DW_COMMENT lookup table to speed up code?
	}


	// If organic content affects soil properties
	if (iforganicsoilproperties) {

		// Calculate organic and mineral fraction
		// This is the fraction of the material (i.e. 1 - porosity) that is organic
		double org_frac = !negligible(organic_frac + mineral_frac) ? organic_frac / (organic_frac + mineral_frac) : 0.0;
		double min_frac = 1.0 - org_frac;

		if (!negligible(org_frac)) {

			// Organis soil parameters
			const double Kw_max_org = soiltype.Kw_max_org;
			const double alpha_org = soiltype.alpha_org;
			const double lambda_org = soiltype.l;
			const double m_org = soiltype.m_org;

			double po_org = pow(the, (1 / m_org));

			// organic soil diffusivity
			double Dw_org = ((1.0 - m_org) * Kw_max_org / (alpha_org * m_org * (W_sat - W_res))) * pow(the, (lambda_org - 1.0 / m_org)) * (pow(1.0 - po_org, m_org) + pow(1.0 - po_org, -m_org) - 2);

			// Organic soil conductivity
			double Kw_org;
			if (W_tot <= 1.001 * W_res) {
				Kw_org = 0.0;
			}
			else if (W_tot > 0.999 * W_sat) {
				Kw_org = Kw_max_org;
			}
			else {
				Kw_org = Kw_max_org * pow(the, lambda_org) * pow(1.0 - pow(1.0 - pow(the, 1.0 / m_org), m_org), 2);
			}

			// Scale diffusivity and conductivity depending on organic fraction
			Dw = min_frac * Dw + org_frac * Dw_org;
			Kw = min_frac * Kw + org_frac * Kw_org;
		}
	}

	// The water transport in frozen soil is limited in the case of a partially frozen soil, by considering the effective 
	// hydraulic conductivity and diffusivity to be a weighted average of the values for total soil water and the value 
	// at the permanent wilting point for frozen water (for convenience set to zero)
	Dw *= (1.0 - I_frac);
	Kw *= (1.0 - I_frac);
}

void Soillayer::update_heatcap_and_cond(bool combine_snow, double snow_depth, double snow_imass, double snow_wmass) {

	// DESCRIPTION
	// Called daily to update heat capacities and thermal conductivities of the 
	// soil layers below the mixed/moss and snow layers 

	double lKfwater;			// log of therm. cond. for the water fraction
	double lKfice;				// log of therm. cond. for the ice fraction
	double lKforg;				// log of therm. cond. for the org. fraction
	double lKfpeat;				// log of therm. cond. for the peat fraction
	double lKfmin;				// log of therm. cond. for the min. fraction
	double Stot;					// total of solids in soil

	double combined_depth = depth;
	double ice = imass;
	double water = wmass;
	double peat = peat_frac * depth;
	double mineral = mineral_frac * depth;
	double organic = organic_frac * depth;

	//acrotelm catotelm DW_COMMENT



	// Soil layer air (mm)
	double air = por - (wmass + imass);

	if (air < -1.0e-8)
		fail("Illegal air mass value in update_heat_diff_and_cond(): %g", air);

	// Catch rounding errors for air
	air = max(0.0, air);

	// Add snow properties if to combine thin snow layer to top soil layer
	if (combine_snow) {
		combined_depth += snow_depth;
		air += snow_depth - (snow_imass + snow_wmass);
		ice += snow_imass;
		water += snow_wmass;
	}

	// Update heat capacity (J m-3 K-1)
	Ch = (mineral * Cp_min + organic * Cp_org + ice * Cp_ice + water * Cp_water + peat * Cp_peat + air * Cp_air) / combined_depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	Stot = mineral + organic + ice + water + peat;

	if ((Stot + air) / combined_depth > 1.000001)
		fail("Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g", Stot / combined_depth, air / combined_depth, (Stot + air) / combined_depth);

	lKforg = (organic / Stot) * log(Korg);
	lKfpeat = (peat / Stot) * log(Kpeat);
	lKfmin = (mineral / Stot) * log(Kmin);
	lKfwater = (water / Stot) * log(Kwater);
	lKfice = (ice / Stot) * log(Kice);

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh = (air * Kair + Stot * exp(lKfwater + lKfice + lKforg + lKfmin + lKfpeat)) / combined_depth;
}

// Get the Water Filled Pore Space of the layer
double Soillayer::wfps() {
	return wmass / wsats;
}

double Soillayer::wcont_max() {
	return (wsats - wp) / awc;
}

// DW_COMMENT need to fill when finished!
void Soillayer::serialize(ArchiveStream& arch) {
	arch & depth
		& wmass
		& imass
		& temp
		& temp_old;

	for (int i = 0; i < NSOMPOOL; i++) {
		arch & sompool[i];
	}

	arch & orgleachfrac
		& orgmovefrac
		& NH4_mass
		& NO3_mass
		& anmin
		& animmob
		& aminleach
		& aorgNleach
		& aorgCleach
		& solvesom
		& fNH4uptake_mean
		& fNO3uptake_mean
		& morgleach_mean
		& morgmove_mean
		& mminleach_mean

		& cpool_slow
		& cpool_fast
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soilcolumn member functions
////////////////////////////////////////////////////////////////////////////////

void Soilcolumn::initdrivers() {

	// Padding layers
	pad.resize(PAD_LAYERS);

	// Set padding layers depth
	set_padding_depth();

	b = set_b(300.0); // standard deviation of orography (300) // DW_COMMENT should be read in!!!
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
	percolate = 0.0;

	// Peatland hydrology variables
	thaw = 0.0;
	maxthawdepththisyear = 0.0;
	awtp = 0.0;
	acro_co2 = 934.0; // [mimol L-1]
	wtd = 0.0; // [-100, +300] mm
	stand_wmass = 0.0;
	ch4_store = 0.0;
	co2_store = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mthaw[mth] = 0.0;
		mwtp[mth] = 0.0;
	}

	for (int d = 0; d < Date::MAX_YEAR_LENGTH; d++) {
		wtp[d] = 0.0;
		dthaw[d] = 0.0;
	}

	dmoss_wtp_limit = 1.0;
	dgraminoid_wtp_limit = 1.0;

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
	solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	solvesomcent_endyr = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);
	savesomcent = false;
}


void Soilcolumn::create_soillayers(std::vector<double> soil_layer_depths) {

	nsoillayers = (int)soil_layer_depths.size();

	double acc_depth = 0.0;
	double depth_profile = 0.0;

	bool acrotelm = false;
	bool catotelm = false;
	nacrotelm = 0;
	ncatotelm = 0;
	thaw = soil_depth;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		acc_depth += soil_layer_depths[lyr];
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

		Soillayer* soillayer = new Soillayer(this, soil_layer_depths[lyr], depth_profile, acrotelm, catotelm);
		soillayers.push_back((*soillayer));
	}

	// Set sizes of water flow calculation parameters
	wpara.set_sizes(nsoillayers);

	// Set sizes of temperature change calculation parameters
	tpara.set_sizes(nsoillayers + PAD_LAYERS);

	wfluxes.resize(nsoillayers + 1, 0.0);
}



std::vector<double> Soilcolumn::get_prop(soilproptype SPT, int m) {
	
	std::vector<double> prop(nsoillayers);

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		if (SPT == DEPTH)
			prop[lyr] = soillayers[lyr].depth; 
		else if (SPT == TEMP)
			prop[lyr] = soillayers[lyr].temp; 
		else if (SPT == WMASS)
			prop[lyr] = soillayers[lyr].wmass;
		else if (SPT == IMASS)
			prop[lyr] = soillayers[lyr].imass;
		else if (SPT == WCONT)
			prop[lyr] = soillayers[lyr].wcont();
		else if (SPT == MWCONT)
			prop[lyr] = soillayers[lyr].mwcont[m];
		else if (SPT == AWCONT)
			prop[lyr] = soillayers[lyr].awcont;
		else if (SPT == WSATS)
			prop[lyr] = soillayers[lyr].wsats;
		else if (SPT == AWC)
			prop[lyr] = soillayers[lyr].awc;
		else if (SPT == WP)
			prop[lyr] = soillayers[lyr].wp;
		else if (SPT == WCONTMAX)
			prop[lyr] = soillayers[lyr].wcont_max();
		else if (SPT == NH_4)
			prop[lyr] = soillayers[lyr].NH4_mass;
		else if (SPT == NO_3)
			prop[lyr] = soillayers[lyr].NO3_mass;
	}
	return prop;
}

std::vector<double> Soilcolumn::get_ccontfrac() {

	std::vector<double> ccontfrac(nsoillayers);

	double cconttot = ccont();

	if (!negligible(cconttot)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			ccontfrac[lyr] += soillayers[lyr].ccont() / cconttot;
		}
	}
	else {
		ccontfrac[0] = 1.0;
	}

	return ccontfrac;
}

std::vector<double> Soilcolumn::get_nfrac(n_pref_type pref) {

	std::vector<double> nfrac(nsoillayers);

	double ntot = nmass_avail(pref);
	double nfracsum = 0.0;

	if (ntot > 1.0e-10) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			nfrac[lyr] += soillayers[lyr].nmass_avail(pref) / ntot;
			nfracsum += soillayers[lyr].nmass_avail(pref) / ntot;
		}
		nfrac[0] += 1.0 - min(1.0, nfracsum);
	}
	else {
		nfrac[0] = 1.0;
	}

	return nfrac;
}

std::vector<double> Soilcolumn::get_rootfrac() {
	std::vector<double> cmass_root(nsoillayers);
	double cmass_root_tot = 0.0;

	Vegetation& vegetation = patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			cmass_root[lyr] += indiv.cmass_root * indiv.pft.rootdist[lyr];
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

std::vector<double> Soilcolumn::get_peatland_rootfrac() {

	// Initialise the root fractions in each layer of the wetland, assuming an exponential decrease of root biomass with depth
	// in accordance with observations from peat cores.

	// Called on the first day of the simulaton.

	// See Wania et al. (2010) - Section 2.2.2

	const double b = 25.1695;
	const double corr = 2.49994026;

	std::vector<double> rootfrac(nsoillayers);

	double layermiddepth = -soillayers[0].depth / 2.0 * CM_PER_MM; // cm

	double sumrootfrac = 0.0;  // DW_COMMENT CHECK rootfrac!!!!!!!

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		if (lyr == nsoillayers - 1) {
			rootfrac[lyr] = 1.0 - sumrootfrac;
		}
		else {
			rootfrac[lyr] = exp(layermiddepth / b) / corr;
			layermiddepth -= (soillayers[lyr].depth + soillayers[lyr + 1].depth) / 2.0 * CM_PER_MM;
		}

		sumrootfrac += rootfrac[lyr];
	}

	return rootfrac;
}

double Soilcolumn::nmass_avail(n_pref_type pref, bool wmass_limited) {

	double nmass_avail = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		nmass_avail += soillayers[lyr].nmass_avail(pref, wmass_limited);
	}
	return nmass_avail;
}

void Soilcolumn::nmass_avail_uptake(double uptake, n_pref_type pref) {
	double tot_nmass_avail = nmass_avail(pref, true);

	if (!negligible(tot_nmass_avail)) {
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].nmass_subtract(uptake * soillayers[lyr].nmass_avail(pref, true) / tot_nmass_avail, pref);
			soillayers[lyr].nmass_remove_small_values(pref);
		}
	}
}

double Soilcolumn::aorgCleach() {
	double aorgC_leach = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		aorgC_leach += soillayers[lyr].aorgCleach;
	}
	return aorgC_leach;
}

double Soilcolumn::aorgNleach() {
	double aorgN_leach = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		aorgN_leach += soillayers[lyr].aorgNleach;
	}
	return aorgN_leach;
}

double Soilcolumn::aminleach() {
	double amin_leach = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		amin_leach += soillayers[lyr].aminleach;
	}
	return amin_leach;
}

double Soilcolumn::anmin() {
	double an_min = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		an_min += soillayers[lyr].anmin;
	}
	return an_min;
}

double Soilcolumn::animmob() {
	double an_immob = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		an_immob += soillayers[lyr].animmob;
	}
	return an_immob;
}

double Soilcolumn::runoff_drain() {
	double runoff_drainage = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		runoff_drainage += max(0.0, soillayers[lyr].runoff_drain);
	}
	return runoff_drainage;
}

double Soilcolumn::runoff_baseflow() {
	return max(0.0, soillayers[nsoillayers - 1].runoff_baseflow);
}

/// C content of soil column
double Soilcolumn::ccont() {
	double ccont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ccont += soillayers[lyr].ccont();
	}
	return ccont;
}

/// N content of soil column
double Soilcolumn::ncont() {
	double ncont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ncont += soillayers[lyr].ncont();
	}
	return ncont;
}

/// Water content of soil column
double Soilcolumn::wcontent() {
	double wcontent = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wcontent += soillayers[lyr].wmass + soillayers[lyr].imass;
	}
	return wcontent;
}

/// Water mass of soil column
double Soilcolumn::wmass() {
	double acc_wmass = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		acc_wmass += soillayers[lyr].wmass;
	}
	return acc_wmass;
}

/// Total cpool slow C content
double Soilcolumn::cpool_slow() {
	double ccont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ccont += soillayers[lyr].cpool_slow;
	}
	return ccont;
}

/// Total cpool slow C content
double Soilcolumn::cpool_fast() {
	double ccont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ccont += soillayers[lyr].cpool_fast;
	}
	return ccont;
}

/// Total cpool slow C content
double Soilcolumn::cmass_sompool(pooltype pool) {
	double ccont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ccont += soillayers[lyr].sompool[pool].cmass;
	}
	return ccont;
}

/// Total cpool slow C content
double Soilcolumn::nmass_sompool(pooltype pool) {
	double ncont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		ncont += soillayers[lyr].sompool[pool].nmass;
	}
	return ncont;
}


void Soilcolumn::set_dwcont() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].dwcont[date.day] = soillayers[lyr].wcont();
	}
}

// Set daily soil layer temperature for this month
void Soilcolumn::set_dtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].dtemp[date.dayofmonth] = soillayers[lyr].temp;
	}

	for (int plyr = 0; plyr < PAD_LAYERS; plyr++) {
		pad[plyr].dtemp[date.dayofmonth] = pad[plyr].temp;
	}
}

void Soilcolumn::set_mwcont() {
	double acc_wcont = 0.0;
	double acc_awc = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].mwcont[date.month] = mean(soillayers[lyr].dwcont + date.day - date.ndaymonth[date.month] + 1,
			date.ndaymonth[date.month]);

		acc_wcont += lyr == 0 ? soillayers[lyr].mwcont[date.month] * soillayers[lyr].awc / 2.0 : 
			(soillayers[lyr-1].mwcont[date.month] * soillayers[lyr-1].awc + soillayers[lyr].mwcont[date.month] * soillayers[lyr].awc) / 2.0;

		acc_awc += lyr == 0 ? soillayers[lyr].awc / 2.0 : (soillayers[lyr - 1].awc + soillayers[lyr].awc) / 2.0; 
		soillayers[lyr].mwcontabove = acc_wcont / acc_awc;
	}
}

// Set this months average soil layer temperature
void Soilcolumn::set_mtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].mtemp[date.month] = mean(soillayers[lyr].dtemp, date.ndaymonth[date.month]);
	}

	for (int plyr = 0; plyr < PAD_LAYERS; plyr++) {
		pad[plyr].mtemp[date.month] = mean(pad[plyr].dtemp, date.ndaymonth[date.month]);
	}
}

double Soilcolumn::set_b(double sigma) {
	// Accounts for the sub-grid variability related to orography
	// Min and max standard deviation of orography
	// (IFS DOCUMENTATION – Cy45r1 - PART IV: PHYSICAL PROCESSES - eqn 8.84)
	double sigma_min = 100.0;
	double sigma_max = 1000.0;
	return (sigma - sigma_min) / (sigma + sigma_max);
}

// Determines static soil layer organic fractions 
void Soilcolumn::set_static_socfrac() {

	// Use a beta value to reproduce the global vertical SOC profile shown in Fig 4 of Jobbagy & Jackson (2010)
	const double beta_SOC = 0.976;
	// This value gives SOC fractions as a fraction of the total in the top 100 cm in this Table:
	// Depth range	:	Calculated fraction	:	Jobbagy & Jackson (2000), Fig 4
	// 00 - 20cm	:	0.42				:	0.42
	// 20 - 40cm	:	0.26				:	0.22
	// 40 - 60cm	:	0.16				:	0.16
	// 60 - 80cm	:	0.10				:	0.11
	// 80 - 100cm	:	0.06				:	0.09

	double depth = soillayers[0].depth * CM_PER_MM; // [mm]

	soillayers[0].soc_frac = 1.0 - pow(beta_SOC, depth); // init first layer
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
double Soilcolumn::gtemp_rootdist(std::vector<double> rootdist) {

	// Calculate root weighted gtemp
	double avr_gtemp = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		avr_gtemp += soillayers[lyr].gtemp * rootdist[lyr];
	}
	return avr_gtemp;
}

/// Performs daily updates to soil temperatures using an analytic function
/**
* Call each simulation day (if iftwolayersoil true) following update of daily air temperature
* prior to canopy exchange and SOM dynamics
* \param climate   The climate to use to update soil water
* \param depth     The depth at which the soil temperature is calculated
*/
void Soilcolumn::calc_temp_analytic() { // DW_COMMENT

	// DESCRIPTION
	// Calculation of soil temperature at depth (m) (usually middle of upper soil layer).
	// Soil temperatures are assumed to follow surface temperatures according to an
	// annual sinusoidal cycle with damped oscillation about a common mean, and a
	// temporal lag.

	// For a sinusoidal cycle, soil temperature at depth z and time t from beginning
	// of cycle given by (Carslaw & Jaeger 1959; Eqn 52; Jury et al 1991):
	//
	//   (1) t(z,t) = t_av + a*exp(-z/d)*sin(omega*t - z/d)
	//
	//   where
	//     t_av      = average (base) air/soil temperature
	//     a         = amplitude of air temp fluctuation
	//     exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z,
	//                 relative to surface temp fluctuation
	//     z/d       = oscillation lag in angular units at soil depth z
	//     z         = soil depth
	//     d         = sqrt(2*k/omega), damping depth
	//     k         = soil thermal diffusivity
	//     omega     = 2*PI/tau, angular frequency of oscillation (radians)
	//     tau       = oscillation period (365 days)
	//
	// Here we assume a sinusoidal cycle, but estimate soil temperatures based on
	// a lag (z/d, converted from angular units to days) relative to air temperature
	// and damping ( exp(-z/d) ) of soil temperature amplitude relative to air
	// temperature amplitude. A linear model for change in air temperature with time
	// during the last 31 days is used to estimate 'lagged' air temperature. Soil
	// temperature today is thus given by:
	//
	//   (2) temp_soil = atemp_mean + exp( -alag ) * ( temp_lag - temp_mean )
	//
	//   where
	//     atemp_mean = mean of monthly mean temperatures for the last year (deg C)
	//     alag       = oscillation lag in angular units at depth 0.25 m
	//                  (corresponds to z/d in Eqn 1)
	//     temp_lag   = air temperature 'lag' days ago (estimated from linear model)
	//                  where 'lag' = 'alag' converted from angular units to days
	//
	// Soil thermal diffusivity (k) is sensitive to soil water content and is estimated
	// monthly based on mean daily soil water content for the past month, interpolating
	// between estimates for 0, 15% and 100% AWHC (Van Duin 1963; Jury et al 1991,
	// Fig 5.11).

	Climate climate = patch.get_climate();

	// conversion factor for soil thermal diffusivity from mm2/s to m2/day
	const double DIFFUS_CONV = 0.0864;

	const double HALF_OMEGA = 8.607E-3; // corresponds to omega/2 = pi/365 (Eqn 1)
										// soil depth at which to estimate temperature (m)
										// const double DEPTH = SOILDEPTH_UPPER * 0.0005;
										// conversion factor for oscillation lag from angular units to days (=365/(2*PI))		
	const double LAG_CONV = 58.09;

	double a, b; // regression parameters
	double k; // soil thermal diffusivity (m2/day)
	double temp_lag; // air temperature 'lag' days ago (see above; deg C)
	double day[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };

	if ((date.year == 0 || date.year == patch.stand.first_year) && date.month == 0 && !date.islastday) {

		// First month of simulation, use air temperature for soil temperature
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			soillayers[lyr].temp = climate.temp;
		}
	}
	else {

		double acc_depth = 0;
		double depth, above_k, current_k;

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			if (date.islastday) {

				acc_depth += (lyr == 0 ? soillayers[lyr].depth / 2.0 : (soillayers[lyr - 1].depth + soillayers[lyr].depth) / 2.0);

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

					k = (k * (acc_depth - (soillayers[lyr - 1].depth + soillayers[lyr].depth) / 2.0) +
						above_k * soillayers[lyr - 1].depth / 2.0 + current_k * soillayers[lyr].depth / 2.0) / acc_depth;
				}

				// Convert accumulated depth to meter
				depth = acc_depth * M_PER_MM;

				// Calculate parameters alag and exp(-alag) from Eqn 2

				soillayers[lyr].alag = depth / sqrt(k / HALF_OMEGA); // from Eqn 1
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
			soillayers[lyr].temp = climate.atemp_mean + soillayers[lyr].exp_alag*(temp_lag - climate.atemp_mean);
		}
	}
}

double Soilcolumn::to_depth(double depth, std::vector<double> var, std::vector<double> scal) {

	// Calculate variable values for layers down to depth
	double acc_depth = 0.0;
	double acc_var = 0.0;
	double acc_scalar = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < var.size()) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		double scalar = !negligible((double)scal.size()) ? scal[lyr] : 1.0;
		acc_var += var[lyr] * scalar * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		acc_scalar += scalar * layer_depth / soillayers[lyr].depth;
		lyr++;
	}
	return !negligible((double)scal.size()) ? acc_var / acc_scalar : acc_var;
}

double Soilcolumn::to_depth(double depth, std::vector<double> var) { // DW_COMMENT Don't want this one!

	// Calculate variable values for layers down to depth
	double acc_depth = 0.0;
	double acc_var = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < var.size()) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_var += var[lyr] * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_var;
}

// add water to a certain soil layer and return the water that cannot be added due to limited capacity
double Soilcolumn::add_wmass(int layer, double extrawmass) {

	Soillayer& soillayer = soillayers[layer];

	double overflow = 0.0;
	
	if (soillayer.wmass + soillayer.imass > soillayer.wsats) {
		overflow = soillayer.wmass + soillayer.imass - soillayer.wsats;
		soillayer.wmass -= overflow;
	}
	else {
		overflow = 0.0;
	}

	return overflow;
}

void Soilcolumn::update_moisture_diff_and_cond() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_moisture_diff_and_cond();
	}
}

// Update heat capacity and conductivity
void Soilcolumn::update_heatcap_and_cond(bool combine_snow) {

	// Update heat capacity and conductivity for each snow layer
	if (!combine_snow)
		snow.update_heatcap_and_cond();

	// Update heat capacity and conductivity for each soil layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		if (lyr == 0 && combine_snow) {
			soillayers[lyr].update_heatcap_and_cond(combine_snow, snow.snowlayers[0].depth, snow.snowlayers[0].imass, snow.snowlayers[0].wmass);
		}
		else {
			soillayers[lyr].update_heatcap_and_cond(); // DW_COMMENT NEED TO CHECK THESE!!!!!!!!!!!
		}
	}

	// Calculate heat capacity and conductivity for each padding layer
	update_heatcap_and_cond_for_pad();
}

void Soilcolumn::update_heatcap_and_cond_for_pad() {

	// Bedrock parameters from Chadburn et al. (2015)
	double Kbedrock = 8.6;			// W m-1 K-1 s-1
	double Cbedrock = 2100000.0;	// J m-3 K-1

	const int first_bedrock_padding_layer = 3;

	for (int lyr = 0; lyr < PAD_LAYERS; lyr++) {
		if (lyr < first_bedrock_padding_layer) {
			pad[lyr].Kh = soillayers[nsoillayers - 1].Kh;
			pad[lyr].Ch = soillayers[nsoillayers - 1].Ch;
		}
		else {
			pad[lyr].Kh = Kbedrock;
			pad[lyr].Ch = Cbedrock;
		}
	}
}

void build_diags_rhs(double dt, double f_top, double f_bot,
	std::vector<double> cap, std::vector<double> diff, std::vector<double> cond,
	std::vector<double>& upper, std::vector<double>& main, std::vector<double>& lower, std::vector<double>& rhs,
	std::vector<double> z, std::vector<double> dz, std::vector<double> psi, std::vector<double> sink,
	int n_soil_layers) {

	// Upper boundary conditions
	double l_p12 = diff[0] + (diff[1] - diff[0]) / (z[1] - z[0]) * dz[0] / 2.0; // Interpolate lambda_+1/2
	double g_p12 = cond[0] + (cond[1] - cond[0]) / (z[1] - z[0]) * dz[0] / 2.0; // Interpolate gamma_+1/2
	double beta_p = dt / 2.0 * l_p12 / dz[0] / (z[1] - z[0]) / cap[0];

	lower[0] = 1.0e6; // Not used, initialize to known value
	main[0] = 1.0 + beta_p;
	upper[0] = -beta_p;

	rhs[0] = (1. - beta_p) * psi[0] + beta_p * psi[1] - dt / cap[0] / dz[0] * (g_p12 - f_top) + sink[0] * dt;

	// Middle terms
	double l_m12, g_m12, beta_m;
	for (int i = 1; i < n_soil_layers - 1; i++) {
		l_m12 = l_p12;
		l_p12 = diff[i] + (diff[i + 1] - diff[i]) / (z[i + 1] - z[i]) * dz[i] / 2.0;
		beta_m = beta_p * dz[i - 1] / dz[i] * cap[i - 1] / cap[i];
		beta_p = dt / 2.0 * l_p12 / dz[i] / (z[i + 1] - z[i]) / cap[i];

		g_m12 = g_p12;
		g_p12 = cond[i] + (cond[i + 1] - cond[i]) / (z[i + 1] - z[i]) * dz[i] / 2.;

		lower[i] = -beta_m;
		upper[i] = -beta_p;
		main[i] = 1.0 + beta_m + beta_p;
		rhs[i] = beta_m * psi[i - 1] + (1.0 - beta_m - beta_p) * psi[i] + beta_p * psi[i + 1]
			- dt / cap[i] / dz[i] * (g_p12 - g_m12) + sink[i] * dt;
	}

	beta_m = beta_p * dz[n_soil_layers - 2] / dz[n_soil_layers - 1] * cap[n_soil_layers - 2] / cap[n_soil_layers - 1];
	g_m12 = g_p12;

	upper[n_soil_layers - 1] = 1.0e6; // Dummy value
	main[n_soil_layers - 1] = 1.0 + beta_m;
	lower[n_soil_layers - 1] = -beta_m;

	rhs[n_soil_layers - 1] = beta_m * psi[n_soil_layers - 2] + (1.0 - beta_m) * psi[n_soil_layers - 1]
		- dt / cap[n_soil_layers - 1] / dz[n_soil_layers - 1] * (f_bot - g_m12)
		+ sink[n_soil_layers - 1] * dt;
}

void solve_tridiagonal(std::vector<double> upper, std::vector<double> main, std::vector<double> lower, std::vector<double> rhs, std::vector<double>& solution) {

	int n = (int)main.size();
	assert(n == solution.size());

	for (int i = 1; i < n; i++) {
		double w = lower[i] / main[i - 1];
		main[i] -= w * upper[i - 1];
		rhs[i] -= w * rhs[i - 1];
	}

	solution[n - 1] = rhs[n - 1] / main[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		solution[i] = (rhs[i] - upper[i] * solution[i + 1]) / main[i];
	}
}

//
/*
 *
 */
void Soilcolumn::wflow_temp() {

	const double EPS = 1.0e-12;

	// Update top soil layer for infiltarted rain or melt water
	if (wfluxes[0] > EPS) {
		soillayers[0].update_temp(rain_melt_temp, wfluxes[0], Cp_water);
	}
	
	// Update soil layer temperature depending on moving water
	for (int lyr = 1; lyr < nsoillayers + 1; lyr++) {
		if (wfluxes[lyr] > EPS && lyr < nsoillayers) {
			soillayers[lyr].update_temp(max(0.0, soillayers[lyr - 1].temp), wfluxes[lyr], Cp_water);
		}
		else if (wfluxes[lyr] < -EPS) {
			soillayers[lyr - 1].update_temp(max(0.0, soillayers[lyr].temp), -wfluxes[lyr], Cp_water);
		}
	}
}

/// Calculate soil column water flow
/* Calculate soil column water flow
 *
 */
void Soilcolumn::calc_wflow(double winfilt, std::vector<double> aet) {

	bool ifbedrock = false;

	// Calculate moisture diffusivity and conductivity
	update_moisture_diff_and_cond();

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wpara.dz[lyr] = soillayers[lyr].depth * M_PER_MM;
		wpara.diff[lyr] = soillayers[lyr].Dw;
		wpara.cond[lyr] = soillayers[lyr].Kw;
		wpara.sink[lyr] = -(aet[lyr] / soillayers[lyr].depth) / SECS_PER_DAY;
		wpara.psi[lyr] = soillayers[lyr].wmass / soillayers[lyr].depth;
	}

	double depth = wpara.dz[0] / 2.0;
	wpara.z[0] = depth;
	for (int lyr = 1; lyr < nsoillayers; lyr++) {
		depth += (wpara.dz[lyr - 1] + wpara.dz[lyr]) / 2.0;
		wpara.z[lyr] = depth;
	}

	// Boundary conditions for soil water
	// Upper boundary condition
	double f_top = winfilt / rho_H2O / SECS_PER_DAY;

	// Lower boundary condition
	double f_bot = ifbedrock ? 0.0 : wpara.cond[nsoillayers - 1];

	// Build diagonals and RHS vector
	build_diags_rhs(SECS_PER_DAY / date.timesteps, f_top, f_bot,
		wpara.cap, wpara.diff, wpara.cond,
		wpara.upper, wpara.main, wpara.lower, wpara.rhs,
		wpara.z, wpara.dz, wpara.psi, wpara.sink, nsoillayers);

	// Solve tridiagonal system
	solve_tridiagonal(wpara.upper, wpara.main, wpara.lower, wpara.rhs, wpara.psi1);

	// Set new wmass
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Remove small negative numbers
		if (wpara.psi1[lyr] > -1.0e-14 && wpara.psi1[lyr] < 0.0) {
			if (lyr < nsoillayers - 1)
				wpara.psi1[lyr + 1] += wpara.psi1[lyr];
			wpara.psi1[lyr] = 0.0;
		}

		soillayers[lyr].wmass = wpara.psi1[lyr] * soillayers[lyr].depth;

		if (soillayers[lyr].wmass < -1.0e-10)
			dprintf("Year %d Day %d Lyr %d FIRST %g\n",date.year,date.day,lyr, soillayers[lyr].wmass);

		// Calculate drainage

		// Drainage of water above saturation including ice
		double drained_wmass = max(0.0, soillayers[lyr].wmass + soillayers[lyr].imass - soillayers[lyr].wsats);

		soillayers[lyr].runoff_drain += drained_wmass;

		// Update total water after drainage 
		// remove drained water from wmass. Water and ice will now maximum fill up soil layer to saturation
		soillayers[lyr].wmass -= drained_wmass;

		if (soillayers[lyr].wmass < -1.0e-10)
			dprintf("Year %d Day %d Lyr %d DRAIN %g Depth %g %g\n", date.year, date.day, lyr, soillayers[lyr].wmass, soillayers[lyr].depth, soillayers[lyr].depth_profile);
	}
}

void padding::serialize(ArchiveStream& arch) {
	arch & temp
		& depth;
}

void Soilcolumn::set_padding_depth() {

	// DESCRIPTION
	// Calculate padding layer depths for Crank-Nicholson bottom boundary
	// condition setup. The padding layer thicknesses are calculated as
	// a geometric sequence with a given total padding thickness and
	// number of layers. This gives a smooth increase in layer depth
	// away from the calculation domain.

	// INPUT PARAMETERS
	// dz;					// Thickness of last double layer.
	// pad_dz[PAD_LAYERS]	// Output padding layer thicknesses [mm].

	// Function, derivative and parameter for geometric factor
	// calculation.
	// double func, dfunc, k

	/// Total depth [mm] of padding layers
	const double PAD_DEPTH = 48000;

	// use the global variables PAD_DEPTH and PAD_LAYERS
	double pad_depth = PAD_DEPTH;

	// Iterative improvements of geometric factor
	double k_guess, k_last;

	// Layer counter.
	int layer;

	// Statement functions for function and derivative for Newton
	// solver for geometric sequence factor.

	// Initial state - don't choose 1.0 as the initial value for
	// k_guess, since it's a trivial solution//
	k_guess = 10.0;
	k_last = 0.0;

	double func, dfunc;
	double dz = 10000.0;// soillayers[nsoillayers - 1].depth; // DW_COMMENT

	// Do Newton-Raphson solve to find geometric sequence factor.
	while (fabs(k_guess - k_last) > 1E-2) {
		k_last = k_guess;

		func = pow(k_guess, (PAD_LAYERS + 1)) - k_guess * (1.0 + pad_depth / dz) + pad_depth / dz;
		dfunc = (PAD_LAYERS + 1) * pow(k_guess, PAD_LAYERS) - (1.0 + pad_depth / dz);
		k_guess = k_guess - func / dfunc;
	}

	// Build geometric sequence of padding layer depths.
	for (layer = 0; layer<PAD_LAYERS; layer++) {
		pad[layer].depth = dz * pow(k_guess, layer + 1);
	}
}

void Soilcolumn::update_acrotelm_co2(double atmo_co2) {

	// DESCRIPTION
	// Update the CO2 concentration in the acrotelm water, after Wania et al. (2009)
	// Varies between atmospheric CO2 level when awtp = -300mm, to PORE_WATER_CO2 when awtp = 0, 
	// and greater if there's standing water.
	// Called once daily for wetlands, after the water table has been updated

	// LPJ-WHy comments from Rita Wania:

	// LPJ-WHy uses a different atmospheric CO2 concentration for mosses
	// as they can access the CO2 dissolved in acrotelm water. If the
	// water table is high, mosses can access all of the acrotelm CO2,
	// but as the water table drops the CO2 concentrations available to
	// mosses is a mixture between atmospheric CO2 and acrotelm CO2.
	// NOTE: for ideal gases, ppmv = not just mimol L-1 but also
	// mimol mol-1.
	// The mean annual wtd of the previous year is used here, as the
	// water table is only calculated after the co2 concentration is
	// needed

	acro_co2 = min(PORE_WATER_CO2, PORE_WATER_CO2 + ((PORE_WATER_CO2 - atmo_co2) / ACROTELM_DEPTH) * awtp); // DW_COMMENT check so it is correct
}

/// Updates the water content for each soil layer in peatland soils
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
void Soilcolumn::hydrology_peat(const Climate& climate, double fevap) {

	// DESCRIPTION
	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET) and percolation between layers;
	// calculation of runoff
	// Uses peatland hydrology scheme from Granberg et al. (1999) and Wania et al. (2009a, 2009b)

	// INPUT PARAMETERS
	// climate    = climate today (mm)
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface


	double aet_total = 0.0;		// total AET (mm)
	double runoff_surf;			// runoff from upper soil layer (mm)
	double runoff_drain;		// runoff (drainage) from lower soil layers (mm)
	double evap;				// evaporation from soil surface (mm)
	double evapotranspiration;	// evapotranspiration (mm/day)
	double fpc_moss_total = 0.0;// Total fpc for mosses

	// total soilcolumn water (mm)
	double wtot = to_depth(soil_depth, get_prop(WMASS));
	// total water in acrotelm layers (mm)
	double wtot_acrotelm = to_depth(ACROTELM_DEPTH, get_prop(WMASS));
	// total soilcolumn ice (mm)
	double itot = to_depth(soil_depth, get_prop(IMASS));
	// total ice in acrotelm layers (mm)
	double itot_acrotelm = to_depth(ACROTELM_DEPTH, get_prop(IMASS));

	/// surface water content when WTD is at 20 mm
	const double surfw = 0.34;
	/// minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// gradient limit of upper gradient [mm]
	const double grad_lim = 50.0;
	/// The soil water characteristics in the upper peat are considered
	/// to be linear in the suction interval 0 - 0.1 m
	/// Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22
	double az = (acrotelm_por - minvtot) / zmin;

	
	// Transpiration

	// Retrieve Vegetation for this patch
	Vegetation& vegetation = patch.vegetation;

	// Sum AET for across all vegetation individuals in this patch

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		// Total aet, assumed from acrotelm and catotelm combined 
		aet_total += indiv.aet;

		if (indiv.pft.ismoss())
			fpc_moss_total += indiv.fpc;

		vegetation.nextobj();
	}

	
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

		// An assumption is that the catotelm should be constantly saturated
		// So AET affecting the catotelm should be replenished with water from 
		// the acrotelm. Hence all AET could directly be taken from the acrotelm.

		// Check so evapotranspiration doesn't exceed available water in the acrotelm
		// If so, then decrease evap accordingly (if negative, water will be added as runon)
		if (aet_total + evap > wtot_acrotelm) {
			evap -= aet_total + evap - wtot_acrotelm;
		}
		evapotranspiration = aet_total + evap;
	}


	// Update acrotelm and total water volume

	// Drainage - could possibly be read in for site-specific studies
	runoff_drain = 0;

	// Water fluxes in and out of the acrotelm 
	wtot_acrotelm += rain_melt - evapotranspiration - runoff_drain;
	wtot += rain_melt - evapotranspiration - runoff_drain;
	rain_melt = 0.0;


	// Water runoff (Granberg et al. 1999 eqn 10, table 1)
	double ideal_runoff = exp(-0.01 * wtd);
	if (wtot_acrotelm > 0.0 && soillayers[0].imass / soillayers[0].depth < 0.7)
		runoff_surf = min(ideal_runoff, wtot_acrotelm);
	else
		runoff_surf = 0.0;

	// .ins file option: 
	// run on or off - but only when there's flowing water

	// Water runon
	// wetland_runon set in global.ins
	double runon = 0.0;
	if (runoff_surf > 0.0) {

		// e.g. 3 mm/day for bog-like/hummock conditions
		if (wetland_runon < 0)
			runoff_surf -= wetland_runon;	

		// 3-5 mm/day for fen-like conditions
		if (wetland_runon > 0)
			runon += wetland_runon;		
	}

	double wtot_acrotelm_init = wtot_acrotelm;

	wtot_acrotelm += runon - runoff_surf + stand_wmass;
	wtot += runon - runoff_surf + stand_wmass;

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

		wtp[date.day] = wtot_acrotelm + itot_acrotelm - ACROTELM_DEPTH * acrotelm_por;
		
		// wtp > 0 above the surface
		if (wtp[date.day] > maxh) {
			// Update runoff from the surface
			runoff_surf += wtp[date.day] - maxh;
			wtp[date.day] = maxh;
		}

		// wtd defined as < 0 above the surface
		wtd = -wtp[date.day];

		wtot_acrotelm = ACROTELM_DEPTH * acrotelm_por - itot_acrotelm;
		wtot = ACROTELM_DEPTH * acrotelm_por + (soil_depth - ACROTELM_DEPTH) * catotelm_por - itot;

		// wtot does NOT include standing water
		// we assume no standing water 
		stand_wmass = 0.0;
	}
	else {

		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * az));

		if (wtd > zmin)
			wtd = 3.0 * (ACROTELM_DEPTH * acrotelm_por - (wtot_acrotelm + itot_acrotelm)) / (2.0 * (acrotelm_por - minvtot));
		
		wtd = min(ACROTELM_DEPTH, wtd);

		wtp[date.day] = -wtd;

		// wtot does NOT include standing water
		// we assume no standing water 
		stand_wmass = 0.0;
	}


	// Water content in each acrotelm soil layer depending on wtd

	// Depth of acrotelm sublayers - 10 mm
	double Dz_sub = 10.0;

	// Number of sublayers in acrotelm
	int nsublayers_acro = (int)(min(ACROTELM_DEPTH, soil_depth) / Dz_sub);

	// Water content in each sublayer
	double surfw_sat = 0.0;
	std::vector<double> value(nsublayers_acro);
	std::vector<double> sub_water(nsublayers_acro);

	// Water table below surface
	if (wtd > 0.0) {

		// Determine the subsurface layer (0 to nsublayers_acro) within which wtd lies. 
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
				double z = (lyr + 1) * Dz_sub; // mm
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
		double wtot_acrotelm_new = 0.0;
		int sublyr = 0;
		for (int lyr = 0; lyr < nacrotelm; lyr++) {

			// Number of sublayers in acrotelm layer
			int sublayers_per_layer = (int)(soillayers[lyr].depth / Dz_sub);

			double layerwater = 0.0;
			for (int sl = 0; sl < sublayers_per_layer; sl++) {
				layerwater += sub_water[sublyr] * Dz_sub; // add mm of water for this sublayer 
				sublyr++;
			}

			if (soillayers[lyr].wsats < layerwater)
				int sch = 0;

			// Water should exceed saturation
			soillayers[lyr].wmass = min(soillayers[lyr].wsats, layerwater) - soillayers[lyr].imass;
			wtot_acrotelm_new += min(soillayers[lyr].wsats, layerwater) - soillayers[lyr].imass;
		}

		// If new total acrotelm water is different than before, adjust water content accordingly
		if (!negligible(fabs(wtot_acrotelm_new - wtot_acrotelm))) {
			double werror = wtot_acrotelm_new - wtot_acrotelm;

			if (werror < 0.0) {
				// Too little water, add missing water from bottom acrotelm layer
				double lyr = nacrotelm - 1;
				while (!negligible(werror)) {

					// Add water no more than to saturation
					double wadd = min(soillayers[lyr].wsats - (soillayers[lyr].wmass + soillayers[lyr].imass), -werror);
					soillayers[lyr].wmass += wadd;
					werror += wadd;
					lyr--;
				}
			}
			else {
				// Too much water, remove excees from top acrotelm layer
				double lyr = 0;
				while (!negligible(werror)) {

					// Remove no more than water in layer
					double wremove = min(soillayers[lyr].wmass, werror);
					soillayers[lyr].wmass -= wremove;
					werror -= wremove;
					lyr++;
				}
			}
		}
	}
	else {
		// WTD <= 0.0, at or above the surface - standing water

		// acrotelm saturated (subtract ice content)
		for (int lyr = 0; lyr < nacrotelm; lyr++) {
			soillayers[lyr].wmass = soillayers[lyr].wsats - soillayers[lyr].imass;
		}
	}

	// Catotelm assumed to be always saturated
	for (int lyr = nacrotelm; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wmass = soillayers[lyr].wsats - soillayers[lyr].imass;
	}

	// Total water change during calculation of water content for each soil layer
	double wafter = to_depth(soil_depth, get_prop(WMASS));
	double wchange = wafter - wtot;

	if (fabs(wchange) > 1.0e-10)
		dprintf("Year %d Day %d ERROR: Peat hydrology water inbalance %g\n",date.year, date.day, wchange);


	// Report water fluxes

	// Evapotranspiration data
	if (empirical_evapotranspiration) {

		patch.fluxes.report_flux(Fluxes::WEVAP, fevap * evapotranspiration);
		patch.fluxes.report_flux(Fluxes::WAET, (1.0 - fevap) * evapotranspiration);
	}
	else {

		if (evap >= 0.0) {
			patch.fluxes.report_flux(Fluxes::WEVAP, evap);
		}
		else {
			// If acrotelm water couldn't meet catotelm water decrease, 
			// then evap would have become negative and is then assumed 
			// to have been added to the water as runon
			patch.fluxes.report_flux(Fluxes::WADDED, -evap);
		}

		patch.fluxes.report_flux(Fluxes::WAET, aet_total);
	}

	// Surface runoff
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);

	// Runon
	patch.fluxes.report_flux(Fluxes::WADDED, runon);


	/*if (date.year >= 100 && date.year <= 10) {
		dprintf("%d %d %g %g %g %g ", date.year, date.day, aet_total, evap, runoff_surf, wtd);
		for (int lyr = 0; lyr < nsoillayers; lyr++)
			dprintf("%g ", soillayers[lyr].wcont());
		dprintf("\n");
	}*/


	// UPDATE CO2 LEVEL IN ACROTELM
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

	if (wtp[date.day] > wtp_moss_upper) {
		dmoss_wtp_limit = 1.0; // No limit
	}
	else if (wtp[date.day] > wtp_moss_lower) {
		dmoss_wtp_limit = 1.0 - (-wtp[date.day] + wtp_moss_upper) * (1.0 - moss_wtp_limit_lowerlimit) / (wtp_moss_upper - wtp_moss_lower);
	}
	else {
		dmoss_wtp_limit = moss_wtp_limit_lowerlimit; // between -300mm and -280mm
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

	const double wtp_graminoid_upper = -10; // -50.0;		// mm - Wania et al. 2009b have -100
	const double wtp_graminoid_lower = -100; // -200.0;	// mm - Wania et al. 2009b have -300 mm
	const double graminoid_wtp_limit_lowerlimit = 0.0;

	// More restrictive than Wania et al., otherwise graminoids become too dominant as they never suffer from inundation stress
	// The cap drops from 1 as soon as wtp goes below wtp_graminoid_upper
	if (wtp[date.day] >= wtp_graminoid_upper) {
		dgraminoid_wtp_limit = 1.0; // i.e. no limit!
	}
	else if (wtp[date.day] > wtp_graminoid_lower) {
		dgraminoid_wtp_limit = 1.0 - (-wtp[date.day] + wtp_graminoid_upper) * (1.0 - graminoid_wtp_limit_lowerlimit) / (wtp_graminoid_upper - wtp_graminoid_lower);
	}
	else {
		dgraminoid_wtp_limit = graminoid_wtp_limit_lowerlimit; // between -300mm and -200mm
	}


	// Update output variables

	// Reset monthly wtp and mmoss_wtp_limit averages on Jan 1
	if (date.day == 0) {
		for (int mth = 0; mth < 12; mth++) {
			mwtp[mth] = 0.0;
		}
	}

	int mdays = date.ndaymonth[date.month];
	mwtp[date.month] += wtp[date.day] / (double)mdays;

	if (date.day == Date::MAX_YEAR_LENGTH - 1) {
		awtp = 0.0;
		for (int d = 0; d < Date::MAX_YEAR_LENGTH; d++)
			awtp += wtp[d] / (double)Date::MAX_YEAR_LENGTH;
	}

	return; // ...if no problems
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

	// Depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double SOILDEPTH_EVAP = 200.0;
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
	double wbefore = wmass() + rain_melt;

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
	if (snow.snowpack() < 10.0) {	
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
	double wafter = wmass() + runoff_surf + runoff_drain() + runoff_baseflow() + aet_total + evap;
	if (fabs(wbefore-wafter) > 1.0e-10)
		dprintf("Year %d Day %d BUCKETWATER Before %g After %g Diff %g Rain_melt %g\n",date.year,date.day, wbefore, wafter, wbefore - wafter, rain_melt);

	// Reset rain_melt
	rain_melt = 0.0;

	// Save fluxes
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, runoff_drain());
	patch.fluxes.report_flux(Fluxes::WRUNOFF_BASE, runoff_baseflow()); 
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

	// Soil depth for soil infiltration
	const double SOILDEPTH_INFILT = 500.0;

	double wmass_begin = wmass() + rain_melt;

	double evap, runoff_surf;

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

	// Calculate static ice mass and saturation level for infiltation depth
	double imass_infilt = to_depth(SOILDEPTH_INFILT, get_prop(IMASS));
	double wsats_infilt = to_depth(SOILDEPTH_INFILT, get_prop(WSATS));

	// Reset daily infiltation flux, evaporation and surface and drainage runoff 
	// and set initial water fraction to be able to calculate daily water fluxes
	double winfilt = 0.0;
	evap = 0.0;
	runoff_surf = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].runoff_drain = 0.0;
		wpara.psi0[lyr] = soillayers[lyr].wmass / soillayers[lyr].depth;
	}

	// Loop through sub daily timesteps to be able to solve 
	// the vertical flow of water with the Crank - Nicolson method
	for (int ts = 0; ts < date.timesteps; ts++) {

		// Evaporation

		// Evaporation is calculated on a sub daily time step as water content 
		// in top soil layer changes with each time step

		// Calculate wcont for top layer (evaporation layer)
		wcont_evap = min(1.0, soillayers[0].wcont());

		// Calculate total evaporation only if snow depth < 10mm
		double evap_sd = 0.0;
		if (snow.depth() < 10.0 && !negligible(fevap)) {
			evap_sd = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;

			// Remove small evaporations
			if (evap_sd < EPS)
				evap_sd = 0.0;

			evap += evap_sd / date.timesteps;
		}

		// Calculate water influx after rain_melt initially meets the evaporation demand, can be negative (mm/day)
		double winflux_sd = rain_melt - evap_sd;

		// Surface infiltration and runoff

		// Calculate changing water content for infiltration layers
		double wmass_infilt = to_depth(SOILDEPTH_INFILT, get_prop(WMASS));

		// Fraction pore spaces that are filled with water and ice
		double totwcont_infilt = (wmass_infilt + imass_infilt) / wsats_infilt;

		// Calculate runoff if there is any throughfall or snow melt (Balsamo 2009, eqn 6 and 11)
		double runoff_surf_sd = 0.0;
		if (winflux_sd > EPS) {
			runoff_surf_sd = max(0.0, winflux_sd - max(0.0, wsats_infilt * (1.0 - totwcont_infilt)) + wsats_infilt *
				pow(max(0.0, pow(1.0 - totwcont_infilt, 1.0 / (1.0 + b)) -
					winflux_sd / (wsats_infilt * (1.0 + b))), b + 1.0));
		}
		runoff_surf += runoff_surf_sd / date.timesteps;

		// Calculate infiltration (mm/day)
		double winfilt_sd = winflux_sd - runoff_surf_sd;
		winfilt += winfilt_sd / date.timesteps;

		// Calculate water flows and content after knowing surface influx and AET of each soil layer
		calc_wflow(winfilt_sd, aet_layer);
	}

	// Calculate water fluxes to, inbetween and out of soil layers
	wfluxes[0] = winfilt;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wfluxes[lyr + 1] = wfluxes[lyr] + wpara.psi0[lyr] * soillayers[lyr].depth - (soillayers[lyr].wmass + soillayers[lyr].runoff_drain + aet_layer[lyr]);
	}

	// Calculate baseflow
	soillayers[nsoillayers - 1].runoff_baseflow = max(0.0, wfluxes[nsoillayers]);

	// Update soil layer temperatures due to soil water fluxes
	wflow_temp();

	// Reset ran_melt as water has entered the soil
	rain_melt = 0.0;

	// Total runoff
	double runoff = runoff_surf + runoff_drain() + runoff_baseflow();

	// water added when patch.stand.is_true_wetland_stand() should be be subtracted from runoff
	// in proportion to its components
	if (patch.stand.is_true_wetland_stand() && ifsaturatewetlands) {

		if (runoff <= patch.wetland_wmass_added_today && runoff > 0.0) {

			// Not enough runoff to balance the water added to the wetland, i.e. we added more than was run off.

			// Update wetland_water_added by taking back the amount used for runoff.
			patch.wetland_wmass_added_today -= runoff;

			// All runoff removed 
			runoff_surf = 0.0;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				soillayers[lyr].runoff_drain = 0.0;
				soillayers[lyr].runoff_baseflow = 0.0;
			}
			runoff = 0.0;
		}
		else if (runoff > patch.wetland_wmass_added_today  && runoff > 0.0) {

			// Enough runoff to balance the water added to the wetland so we take it back.
			// Recalculate runoff 
			runoff_surf -= patch.wetland_wmass_added_today * runoff_surf / runoff;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				soillayers[lyr].runoff_drain = patch.wetland_wmass_added_today * soillayers[lyr].runoff_drain / runoff;
				soillayers[lyr].runoff_baseflow = patch.wetland_wmass_added_today * soillayers[lyr].runoff_baseflow / runoff;
			}
			
			runoff = runoff_surf + runoff_drain() + runoff_baseflow();

			// Reset wetland_water_added
			patch.wetland_wmass_added_today = 0.0;
		}
	}

	// Save fluxes
	patch.fluxes.report_flux(Fluxes::WRUNOFF_SURF, runoff_surf);
	patch.fluxes.report_flux(Fluxes::WRUNOFF_DRAIN, runoff_drain());
	patch.fluxes.report_flux(Fluxes::WRUNOFF_BASE, runoff_baseflow());
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
	double wmass_end = wmass() + runoff_surf + runoff_drain() + runoff_baseflow() + evap + aet_total + patch.wetland_wmass_added_today;
	
	if (fabs(wmass_begin - wmass_end) > 1.0e-10)
		dprintf("Year %d Day %d WMASS ERROR: Diff %g Begin %g End %g\n", date.year, date.day, wmass_begin - wmass_end, wmass_begin, wmass_end);
}


// Initial percolation
/* Initial percolation
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

				// PEATLAND SOILS
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
					// imass_fraction = 0 (e.g. summer, no phase change etc.) = 1
					// imass_fraction = 1 (e.g. winter, saturated soil to field capacity, all soil water frozen) = 0.000000
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

						if (soillayers[lyr].wmass > soillayers[lyr].wsats)
							fail("initial_infiltration - error too much water added! - is_highlatitude_peatland_stand\n");
					}
				}
			}
		}

		if (patch.stand.is_true_wetland_stand()) {

			// MINERAL WETLANDS

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
			if (ifsaturatewetlands)
				patch.wetland_wmass_added_today = max(0.0, wmass_added - rain_melt);

			rain_melt -= min(rain_melt, wmass_added);

			// END MINERAL WETLANDS
		}
	}
}


void Soilcolumn::phase_change() {

	// DESCRIPTION
	// Calculate water and ice for each layer

	// To know how much ice is formed or melted at each time step,
	// The temperature change from the previous to this time step is first 
	// calculated. If not all of the available water is frozen, the
	// temperature will be held at 0 deg until everything is frozen (and
	// vice versa for the melting). The amount of heat released when 
	// freezing or absorbed when thawing is to be found. To do this
	// the total heat capacity of the soil is multiplied with delta_T. Since
	// we are working in m-3 here, we don't need to multiply by the volume
	// or mass, because we want to get the fraction of water freezing or
	// ice thawing, not the absolute value. 

	// Determine where to start:
	// if there's standing water, include that layer for the phase change 
	// calculation, if not, start with the top soil layer

	double Ffreez;		// fraction of water which freezes
	double Fthaw;		// fraction of ice which thaws
	double Fwmass;		// fraction of water
	double Fimass;		// fraction of ice
	double energy;		// Energy per volume [J m-3]
	double delta_T;		// Temperature difference (current - previous)
	double soil_depth;	// Soil layer depth [mm]
	double soil_temp;	// Soil layer temperature [C]
	const double EPS = 1.0e-10;

	// If snow depth is below one standard depth, then add snow layer to top soil layer
	bool combine_snow = snow.depth() < snow.snow_depth_standard && !negligible(snow.depth());

	// Phase change between water and ice and adjust snow layer temperatures accordingly
	if (!combine_snow) {
		snow.phase_change();
	}

	if (date.year >= FIRST_FREEZE_YEAR) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			// Initial water and ice content
			double wtotal_initial = soillayers[lyr].wmass + soillayers[lyr].imass;

			Ffreez = 0.0;
			Fthaw = 0.0;

			// Soil layer depth
			soil_depth = soillayers[lyr].depth;

			// Soil layer temperature
			soil_temp = soillayers[lyr].temp;

			// Water and ice fraction of full volume
			if (combine_snow && lyr == 0) {
				soil_depth += snow.snowlayers[0].depth;
				Fwmass = (soillayers[lyr].wmass + snow.snowlayers[0].wmass) / soil_depth;
				Fimass = (soillayers[lyr].imass + snow.snowlayers[0].imass) / soil_depth;

				// NOTE: Heat capacity and temperature of top soil layer is already adjusted for thin 
				// snow layer in temperature calculations
			}
			else {
				Fwmass = soillayers[lyr].wmass / soil_depth;
				Fimass = soillayers[lyr].imass / soil_depth;
			}

			// temperature difference
			delta_T = soil_temp - soillayers[lyr].temp_old;

			// FREEZING
			if (soil_temp < 0.0 && soil_temp < soillayers[lyr].temp_old) { // Wania conditions 

				if (Fwmass > 0.0) {

					// Energy to freeze water [J m-3]
					//energy = soillayers[lyr].Ch * fabs(delta_T); // DW_COMMENT need Pauls explanation here!!!!!!!!!
					energy = Cp_water * fabs(delta_T);

					// Fraction water that would freeze [m3 m-3]
					Ffreez = energy / Lheat;

					if (Fwmass > Ffreez) {
						Fimass += Ffreez;
						Fwmass -= Ffreez;

						soillayers[lyr].temp = 0.0;
					}
					else {

						// It is cold enough to freeze all the water and
						// decrease the temperature below 0, since the rest
						// of the water gets frozen now.

						Fimass += Fwmass;

						// Update soil temperature in this layer
						// How much energy (temperature) the phase change needed
						soillayers[lyr].temp += Fwmass * Lheat / Cp_water;

						if (soillayers[lyr].temp > 0.0)
							dprintf("ERROR: Year %d Day %d Soil layer %d temperature is above zero AFTER freezing %g\n", date.year, date.day, lyr, soillayers[lyr].temp);

						soillayers[lyr].temp = min(0.0, soillayers[lyr].temp);

						Fwmass = 0.0;
					}
				}
			} // end of FREEZING

			  // THAWING
			else if (soil_temp >= 0.0  && soil_temp > soillayers[lyr].temp_old) { // Wania conditions

				if (Fimass > 0.0) {

					// Energy to thaw ice [J m-3]
					//energy = soillayers[lyr].Ch * fabs(delta_T);
					energy = Cp_water * fabs(delta_T);

					// Fraction ice that would thaw [m3 m-3]
					Fthaw = energy / Lheat;

					if (Fimass > Fthaw) {

						// There is some ice left, even after melt
						Fwmass += Fthaw;
						Fimass -= Fthaw;
						soillayers[lyr].temp = 0.0;
					}
					else {
						// There is more energy available than required for
						// melting the rest of the ice. Therefore the
						// left-over energy will increase the temperature.

						Fwmass += Fimass;

						// Update soil temperature in this layer
						// How much energy (temperature) the phase change needed
						soillayers[lyr].temp -= Fimass * Lheat / Cp_water;

						if (soillayers[lyr].temp < 0.0)
							dprintf("ERROR: Year %d Day %d Soil layer %d temperature is below zero AFTER thawing %g\n", date.year, date.day, lyr, soillayers[lyr].temp);

						soillayers[lyr].temp = max(0.0, soillayers[lyr].temp);

						Fimass = 0.0;
					}
				}
			} // end of THAWING

			// If top soil layer has been combined with a thin snow layer
			// then the same fraction of water (ice) is freezed (thawed)
			// in both the soil and snow layer
			if (combine_snow && lyr == 0) {
				if (!negligible(Ffreez)) {

					// Fraction of water left after freezing
					double frac = Fwmass * soil_depth / (soillayers[lyr].wmass + snow.snowlayers[0].wmass);

					// Update soil layer water and ice mass
					double imass = (1.0 - frac) * soillayers[lyr].wmass;
					soillayers[lyr].wmass -= imass;
					soillayers[lyr].imass += imass;

					// Update snow layer water and ice mass
					imass = (1.0 - frac) * snow.snowlayers[0].wmass;
					snow.snowlayers[0].wmass -= imass;
					snow.snowlayers[0].imass += imass;
				}
				else if (!negligible(Fthaw)) {

					// Fraction of ice left after thawing
					double frac = Fimass * soil_depth / (soillayers[lyr].imass + snow.snowlayers[0].imass);

					// Update soil layer water and ice mass
					double wmass = (1.0 - frac) * soillayers[lyr].imass;
					soillayers[lyr].imass -= wmass;
					soillayers[lyr].wmass += wmass;

					// Update snow layer water and ice mass
					wmass = (1.0 - frac) * snow.snowlayers[0].imass;
					snow.snowlayers[0].imass -= wmass;
					snow.snowlayers[0].wmass += wmass;
				}

				// Update snow temperature
				snow.snowlayers[0].temp = soillayers[lyr].temp;
			}
			else {
				// Water and ice fraction of full volume
				soillayers[lyr].wmass = Fwmass * soil_depth;
				soillayers[lyr].imass = Fimass * soil_depth;
			}

			double wtotal_final = soillayers[lyr].wmass + soillayers[lyr].imass;

			// water balance checks:
			if ((fabs(wtotal_initial - wtotal_final)) > EPS)
				fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);

		} // end of loop through layers
	}

	// Balance ice and water in snow layers after phase change
	snow.balance_snow(rain_melt, rain_melt_temp);
}

/// Thaw depth
void Soilcolumn::thaw_depth() {

	// Define the thaw depth as the 0 degree C isotherm (Muller 1947, Burn 1998), 
	// which equals the depth, until which liquid water exists continuously from the top down.

	// daily thawing depth, where at least some water is found, but ice can be present [mm]
	thaw = 0.0;

	// Add soil layer where all ice is thawn
	int lyr = 0;
	while (lyr < nsoillayers && soillayers[lyr].temp > 0.0) {
		thaw += soillayers[lyr].depth;
		lyr++;
	}
	// If all soil layer where thawn then add padding layer with temperature above 0 degree C
	if (lyr == nsoillayers) {
		while (lyr >= nsoillayers && lyr < nsoillayers + PAD_LAYERS && pad[lyr - nsoillayers].temp > 0.0) {
			thaw += pad[lyr - nsoillayers].depth;
			lyr++;
		}
	}

	// If still among soil layers then assuming that water in next soil layer is at the top of the layer and ice below
	if (lyr < nsoillayers) {
		if (!negligible(soillayers[lyr].wmass + soillayers[lyr].imass) && lyr < nsoillayers) {
			thaw += soillayers[lyr].wmass / (soillayers[lyr].wmass + soillayers[lyr].imass) * soillayers[lyr].depth;
		}
	}

	// Set maximum thaw depth this year
	if (thaw > maxthawdepththisyear || date.day == 0)
		maxthawdepththisyear = thaw;

	// Save daily thaw depth
	dthaw[date.day] = thaw;

	// Calculate mean monthly thaw depth
	if (date.islastday)
		mthaw[date.month] = mean(dthaw + date.day - date.ndaymonth[date.month] + 1,
			date.ndaymonth[date.month]);
}

// True if phase changes allowed, false otherwise
bool Soilcolumn::can_freeze() const {
	return ifallowphasechanges;
}

void Soilcolumn::calc_temp() {

	// Minimum dz for top snow layer [m] (should be calculated from Von Neumann Stability Analysis k=gamma deltaT/(deltaX)^2) DW_COMMENT

	double airT = patch.get_climate().temp;

	for (int ts = 0; ts < date.timesteps; ts++) {

		int nsnowlayers = snow.nsnowlayers;
		bool combine_snow = false;

		// Save temperature before temperature update, needed for phase changes
		for (int lyr = 0; lyr < snow.nsnowlayers; lyr++) {
			snow.snowlayers[lyr].temp_old = snow.snowlayers[lyr].temp;
		}
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].temp_old = soillayers[lyr].temp;
		}
		
		// If snow layer is too thin add its depth and properties to the top soil layer (Best 2011 section 3.3)
		if (snow.nsnowlayers == 1) {
			if (snow.depth() < snow.snow_depth_standard) {
				combine_snow = true;
				nsnowlayers = 0;
			}
		}

		int nsnowsoillayers = nsnowlayers + nsoillayers;
		int nsnowsoilpadlayers = nsnowlayers + nsoillayers + PAD_LAYERS;

		if (nsnowsoilpadlayers != tpara.dz.size()) {
			tpara.set_sizes(nsnowsoilpadlayers);
		}

		// Update heat capacity and thermal conductivity for snow, soil and padding layers
		update_heatcap_and_cond(combine_snow);

		// Snow layers (0 == bottom layer, nsnowlayers - 1 == top snow layer)
		int slyr = nsnowlayers - 1;
		for (int lyr = 0; lyr < nsnowlayers; lyr++) {
			tpara.dz[lyr] = snow.snowlayers[slyr].depth * M_PER_MM;
			tpara.cap[lyr] = snow.snowlayers[slyr].Ch;
			tpara.diff[lyr] = snow.snowlayers[slyr].Kh;
			tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
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
			tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
			tpara.psi[lyr] = soillayers[lyr - nsnowlayers].temp;
		}

		// Padding layers
		for (int lyr = nsnowsoillayers; lyr < nsnowsoilpadlayers; lyr++) {
			tpara.dz[lyr] = pad[lyr - nsnowsoillayers].depth * M_PER_MM;
			tpara.cap[lyr] = pad[lyr - nsnowsoillayers].Ch;
			tpara.diff[lyr] = pad[lyr - nsnowsoillayers].Kh;
			tpara.cond[lyr] = tpara.diff[lyr] / tpara.cap[lyr];
			tpara.psi[lyr] = pad[lyr - nsnowsoillayers].temp;
		}

		double depth = tpara.dz[0] / 2.0;
		tpara.z[0] = depth;
		for (int lyr = 1; lyr < nsnowsoilpadlayers; lyr++) {
			depth += (tpara.dz[lyr - 1] + tpara.dz[lyr]) / 2.0;
			tpara.z[lyr] = depth;
		}

		// Boundary conditions for soil temperature (Shouldn't use air temperature)
		// Upper boundary condition
		double f_top = -2.0 * tpara.diff[0] * (tpara.psi[0] - airT) / tpara.z[0] / date.timesteps;

		// Lower boundary condition
		double f_bot = 0.0;

		// Build diagonals and RHS vector
		build_diags_rhs(SECS_PER_DAY / date.timesteps, f_top, f_bot,
			tpara.cap, tpara.diff, tpara.cond,
			tpara.upper, tpara.main, tpara.lower, tpara.rhs,
			tpara.z, tpara.dz, tpara.psi, tpara.sink, nsnowsoilpadlayers);

		// Solve tridiagonal system
		solve_tridiagonal(tpara.upper, tpara.main, tpara.lower, tpara.rhs, tpara.psi1);

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

		// Phase change between water and ice and adjust snow and soil layer temperatures accordingly
		phase_change();

	/*	if (snow.nsnowlayers) {
			dprintf("Year %d Day %d SNOW ", date.year, date.day);
			for (int lyr = 0; lyr < snow.nsnowlayers; lyr++) {
				dprintf("%g ", snow.snowlayers[lyr].depth);
			}
			dprintf("\n");
		}*/

		/*if (date.year == 0) {
		dprintf("%d %d %g %g ", date.day, ts, patch.get_climate().temp, f_top);
		for (int lyr = 0; lyr < snow.nsnowlayers + nsoillayers + PAD_LAYERS; lyr++) {
		dprintf("%g ", tpara.psi1[lyr]);
		}
		for (int lyr = 0; lyr < snow.nsnowlayers + nsoillayers + PAD_LAYERS; lyr++) {
		dprintf("%g ", tpara.z[lyr]);
		}
		for (int lyr = 0; lyr < snow.nsnowlayers + nsoillayers + PAD_LAYERS; lyr++) {
		dprintf("%g ", tpara.dz[lyr]);
		}
		for (int lyr = 0; lyr < snow.nsnowlayers + nsoillayers + PAD_LAYERS; lyr++) {
		dprintf("%g ", tpara.cap[lyr]);
		}
		for (int lyr = 0; lyr < snow.nsnowlayers + nsoillayers + PAD_LAYERS; lyr++) {
		dprintf("%g ", tpara.diff[lyr]);
		}
		dprintf("\n");
		}*/

		/*//PRINTING OUT DAILY SNOW AND SOIL PROPERTIES
		if (date.year == 0 && date.day == 0 && ts == 0) {
			dprintf("Year Day cnstep dailyair f_top Snow1 Snow2 Snow3 Snow4 Snow5 Soil1 Soil3 Soil3 Soil4 Soil5 Soil6 snowdz1 snowdz2 snowdz3 snowdz4 snowdz5 soildz1 soildz2 soildz3 cap_C1 cap_C2 cap_C3 cap_C4 cap_C5 cap_C6 cap_C7 cap_C8 diff_C1 diff_C2 diff_C3 diff_C4 diff_C5 diff_C6 diff_C7 diff_C8 \n");
		}
		//
		if (date.year < 5){// && date.day < 100) {

			dprintf("%d %d %d %g %g", date.year, date.day, ts, patch.get_climate().temp, f_top);

			// Temp
			for (int i = 4; i >= 0; i--) {
				if (snow.nsnowlayers - 1 >= i)
					dprintf(" %g", snow.snowlayers[i].temp);
				else
					dprintf(" %g", patch.get_climate().temp);
			}
			dprintf(" %g %g %g %g %g %g", soillayers[0].temp, soillayers[1].temp, soillayers[2].temp, soillayers[3].temp, soillayers[4].temp, soillayers[5].temp);

			// Snow depth
			for (int i = 4; i >= 0; i--) {
				if (snow.nsnowlayers - 1 >= i)
					dprintf(" %g", snow.snowlayers[i].depth);
				else
					dprintf(" %g", 100.0);
			}
			dprintf(" %g %g %g", tpara.dz[tmp_snowlayers], tpara.dz[tmp_snowlayers + 1], tpara.dz[tmp_snowlayers + 2]);

			// cap_C
			for (int i = 4; i >= 0; i--) {
				if (snow.nsnowlayers - 1 >= i)
					dprintf(" %g", snow.snowlayers[i].Ch);
				else
					dprintf(" %g", 0.0);
			}
			dprintf(" %g %g %g", tpara.cap[tmp_snowlayers], tpara.cap[tmp_snowlayers +1], tpara.cap[tmp_snowlayers +2]);

			// diff_K
			for (int i = 4; i >= 0; i--) {
				if (snow.nsnowlayers - 1 >= i)
					dprintf(" %g", snow.snowlayers[i].Kh);
				else
					dprintf(" %g", 0.0);
			}
			dprintf(" %g %g %g", tpara.diff[tmp_snowlayers], tpara.diff[tmp_snowlayers + 1], tpara.diff[tmp_snowlayers + 2]);

			dprintf("\n");
		}*/
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	// Thaw depth calculation
	thaw_depth();
}


// Calculate soil layer temperatures
/*
 *
*/
void Soilcolumn::temp() {

	if (iftwolayersoil) { // DW_COMMENT analytic should work for multilayer but then with old snow!!!!!!!
		// Calculate analytic soil temperatures and update dtemp_31 (trend in daily air temperatures for the last 31 days)
		// Analytic soil temperatures can be used to test the accuracy of the numerical scheme in Soilcolumn::calc_temp below.
		calc_temp_analytic();
	}
	else {

		// Calculate the soil temperature in snow, soil and padding layers.
		calc_temp();
	}

	// Save daily soil temperature
	set_dtemp();

	// On last day of month, calculate mean soil temperature for last month
	if (date.islastday) {
		set_mtemp();
	}
}

void Soilcolumn::save_snow() {

	// Save daily snow depth
	snow.ddepth[date.day] = snow.depth();

	// Calculate mean monthly snow depth
	if (date.islastday) {
		snow.mdepth[date.month] = mean(snow.ddepth + date.day - date.ndaymonth[date.month] + 1,
			date.ndaymonth[date.month]);
	}

	for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
		if (lyr < snow.nsnowlayers) {
			snow.mtemp[lyr][date.month] += snow.snowlayers[lyr].temp / date.ndaymonth[date.month];
		}
		else {
			snow.mtemp[lyr][date.month] += patch.get_climate().temp / date.ndaymonth[date.month];
		}
	}
}

void Soilcolumn::serialize(ArchiveStream& arch) {
	arch & soillayers
		& anfix
		& anfix_mean
		& anfix_calc
		& pad
		& solvesomcent_beginyr
		& solvesomcent_endyr;
}



///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Fukusako, S. Int J Thermophys (1990) 11: 353. doi:10.1007/BF01133567
// Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., 
//  and Friedlingstein, P.: An improved representation of physical permafrost dynamics 
//  in the JULES land - surface model, Geosci.Model Dev., 8, 1493 - 1508, 
//  https ://doi.org/10.5194/gmd-8-1493-2015, 2015
// Jobbagy, E. G. & Jackson, R. B. (2000) The vertical distribution of soil carbon
//  and its relation to climate and vegetation. Ecological Applications 10(2): 423-436
// Lawrence, D. M. and Slater, A.: Incorporating organic soil into a global climate 
//  model, Clim.Dynam., 30, 145–160, doi:10.1007/s00382-007-0278-1, 2008.
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
// For Arctic and wetland references see soil.h
