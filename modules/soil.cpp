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
	ice = 0.0;
	water = 0.0;
	dens = 0.0;
	depth = 0.0;
	whc = 0.0;
	temp = 0.0;
	temp_old = 0.0;
	Kh = 0.0;
	Ch = 0.0;
	Dh = 0.0;
}

void Snowlayer::update_heatcap_and_cond() {

	// Snow density including water withing snow
	double tot_dens = min(rho_ICE, (ice + water) / depth * MM_PER_M);

	// C Snow - Fukusako, Eq. 2. [J m-3 K-1]
	Ch = (0.185 + 0.689 * (K2degC + temp) * 0.01) * tot_dens * J_PER_KJ;

	// K Snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = 2.22 * pow((tot_dens / rho_H2O), 1.88);

	// D Snow
	Dh = Kh / Ch;
}

/// weigthed avg of new and old snow density
void Snowlayer::mix_dens(double prec, double snowfall_dens) {
	double sum = ice + water + prec;
	double w_new = prec / sum;
	double w_old = 1.0 - w_new;
	dens = min(rho_ICE, w_old * dens + w_new * snowfall_dens);
}

void Snowlayer::serialize(ArchiveStream& arch) {
	arch & ice
		& water
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
	melt = 0.0;
	NO3_mass = 0.0;
	NH4_mass = 0.0;
	snow_depth_standard = 50.0;
	snow_depth_limit = 1.2;
	dec_depth = 0.0;

	for (int mth = 0; mth < 12; mth++) {
		mdepth[mth] = 0.0;
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

double Snow::depth() {
	double sd = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sd += snowlayers[lyr].depth;
	}
	return sd;
}

double Snow::snowpack() {
	double sp = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sp += snowlayers[lyr].ice + snowlayers[lyr].water;
	}
	return sp;
}

void Snow::update_heatcap_and_cond() {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].update_heatcap_and_cond();
	}
}

/// Potential redistribution of ice/dens after melt
void Snow::potential_ice(double current_ice, std::vector<double>& pot_ice) {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		if (current_ice > snow_depth_standard * snow_depth_limit) {
			pot_ice[lyr] = snow_depth_standard;
			current_ice -= snow_depth_standard;
		}
		else if (current_ice <= snow_depth_standard * snow_depth_limit) {
			pot_ice[lyr] = current_ice;
			current_ice = 0.0;
		}
	}
}

/// Re-distributing snow after freeze/melt event
void Snow::melt_dens(double current_ice) {
	double initial_water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_water += snowlayers[lyr].ice + snowlayers[lyr].water;
	}

	// Potential ice
	std::vector<double> pot_ice(nsnowlayers);

	// Check potential ice content
	potential_ice(current_ice, pot_ice);

	// Re-calc ice & dens of layers
	for (int lyr = 0; lyr < nsnowlayers - 1; lyr++) {
		int lyr_above = lyr;

		// Move snow downwards in snowpack to fill every layer to pot_ice
		while (snowlayers[lyr].ice != pot_ice[lyr] && lyr_above < nsnowlayers - 1) {
			lyr_above++;

			if (!negligible(pot_ice[lyr])) {

				// Either bottom layer shrinking from being larger than standard_depth
				// or water freezing making layer exceed potential
				if (snowlayers[lyr].ice > pot_ice[lyr]) {
					double ice_move = snowlayers[lyr].ice - pot_ice[lyr];

					// Calculate new snow density
					double ice_sum = ice_move + snowlayers[lyr_above].ice;
					double w_new = ice_move / ice_sum;
					double w_old = 1.0 - w_new;
					snowlayers[lyr_above].dens = min(rho_ICE, w_old * snowlayers[lyr_above].dens + w_new * snowlayers[lyr].dens);

					// Add and remove water
					double water_move = ice_move / snowlayers[lyr].ice * snowlayers[lyr].water;
					snowlayers[lyr_above].water += water_move;
					snowlayers[lyr].water -= water_move;

					// Add and remove ice
					snowlayers[lyr_above].ice += ice_move;
					snowlayers[lyr].ice -= ice_move;

					// If new layer created
					if (snowlayers[lyr_above].Dh == 0.0) {
						snowlayers[lyr_above].Kh = snowlayers[lyr].Kh;
						snowlayers[lyr_above].Ch = snowlayers[lyr].Ch;
						snowlayers[lyr_above].Dh = snowlayers[lyr].Dh;
						snowlayers[lyr_above].dens = snowlayers[lyr].dens;
					}
				}
				else {

					double ice_move = max(0.0, min(snowlayers[lyr_above].ice, pot_ice[lyr] - snowlayers[lyr].ice));

					// Add ice to be moved if it isn't negligible
					if (!negligible(ice_move) && ice_move > 0.0) {

						// Calculate new snow density
						double ice_sum = ice_move + snowlayers[lyr].ice;
						double w_new = ice_move / ice_sum;
						double w_old = 1.0 - w_new;
						snowlayers[lyr].dens = min(rho_ICE, w_old * snowlayers[lyr].dens + w_new * snowlayers[lyr_above].dens);

						// Add and remove water
						double water_move = ice_move / snowlayers[lyr_above].ice * snowlayers[lyr_above].water;
						snowlayers[lyr].water += water_move;
						snowlayers[lyr_above].water -= water_move;

						// Add and remove ice

						snowlayers[lyr].ice += ice_move;
						snowlayers[lyr_above].ice -= ice_move;
					}
				}
			}
			else {
				// Erase empty snow layer
				snowlayers.erase(snowlayers.begin() + lyr);
			}
		}
	}

	double final_water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_water += snowlayers[lyr].ice + snowlayers[lyr].water;
	}

	const double EPS = 1.0e-8;
	double diff = abs(initial_water - final_water);


	//dprintf(" balance %d %d  %g %g %g \n",date.year, date.day, initial_water, final_water, diff);
	assert(abs(initial_water - final_water) < EPS);
}

/// Re-distributing snow after freeze/melt event, when snowpack > max potential snowdepth
void Snow::melt_dens_full() {

	// filling each layer to standard_depth
	for (int lyr = NLAYERS_SNOW_MAX - 1; lyr > 0; lyr--) {
		int lyr_below = lyr;

		// Move snow upwards in snowpack to fill every layer to standard_depth
		while (snowlayers[lyr].ice != snow_depth_standard) {
			lyr_below--;

			// Move upwards maximum of what is in snow layer lyr_below
			double ice_move = min(snowlayers[lyr_below].ice, snow_depth_standard - snowlayers[lyr].ice);

			// New snow depth with added snow
			double new_depth = min(snow_depth_standard, snowlayers[lyr].ice + ice_move);

			// Calculate new snow density
			double w_new = ice_move / new_depth;
			double w_old = 1.0 - w_new;
			snowlayers[lyr].dens = min(rho_ICE, w_old * snowlayers[lyr].dens + w_new * snowlayers[lyr].dens);

			// Add and remove water
			double water_move = ice_move / snowlayers[lyr_below].ice * snowlayers[lyr_below].water;
			snowlayers[lyr].water += water_move;
			snowlayers[lyr_below].water -= water_move;

			// Add and remove ice
			snowlayers[lyr].ice += ice_move;
			snowlayers[lyr_below].ice -= ice_move;
		}
	}
}

// Update snow water holding capacity
void Snow::update_whc() {
	// Calculate layer water holding capacity, based on eqn 10 Wang et al. 2013	
	double dens_t = 400;	// kg m-3, reference density
	double rw_min = 0.03;
	double rw_max = 0.1;

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].whc = snowlayers[lyr].ice * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - snowlayers[lyr].dens) / dens_t);
	}
}

/// Balance ice and water in snow layers after melting and refreezing
void Snow::ice_water_snow_balance() {

	// No snow layers then leave function
	if (!nsnowlayers)
		return;

	double initial_snow = 0.0;
	double sum_ice = 0.0;
	initial_snow = melt;

	// Sum up ice in snowpack
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sum_ice += snowlayers[lyr].ice;
		initial_snow += snowlayers[lyr].ice + snowlayers[lyr].water;
	}

	// Re-distribute layers (only based on ICE!)
	double max_snow = NLAYERS_SNOW_MAX * snow_depth_standard;

	// Balance layers depending on if snowpack is over maxumum or not
	if (sum_ice != 0) {
		if (sum_ice > max_snow) {
			melt_dens_full();
		}
		else {
			melt_dens(sum_ice);
		}
	}

	// Check for insignificant top layers	
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		// If layer is too thin then add it to layer below
		if (snowlayers[lyr].ice < snow_depth_standard * (snow_depth_limit - 1.0)) {

			// Check so this layer isn't the bottom layer
			if (lyr != 0) {

				// Calculate old and new snow fractions
				double ice_sum = snowlayers[lyr - 1].ice + snowlayers[lyr].ice;
				double ice_new = !negligible(ice_sum) ? snowlayers[lyr].ice / ice_sum : 0.0;
				double ice_old = 1.0 - ice_new;

				// Update ice, water and snow density in layer below
				snowlayers[lyr - 1].ice += snowlayers[lyr].ice;
				snowlayers[lyr - 1].water += snowlayers[lyr].water;
				snowlayers[lyr - 1].dens = ice_new * snowlayers[lyr].dens + ice_old * snowlayers[lyr - 1].dens;

				// Erase this empty snow layer
				snowlayers.erase(snowlayers.begin() + lyr);
				nsnowlayers = snowlayers.size();
			}
		}
	}

	// Update snow water holding capacity
	update_whc();

	// Redistribute liquid water within the snow pack
	// moving from top to bottom snow layers

	// Save this timesteps melt water above whc which will be pass downwards
	double ts_melt = 0.0;

	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		// No ice
		if (snowlayers[lyr].ice == 0.0) {
			snowlayers[lyr].whc = 0.0;
			ts_melt += snowlayers[lyr].water;
			snowlayers[lyr].water = 0.0;

			// Erase this empty snow layer
			snowlayers.erase(snowlayers.begin() + lyr);
		}
		else {
			// Fill snow layer with all melt water from above
			snowlayers[lyr].water += ts_melt;
			ts_melt = 0.0;

			// Pass water above holding capacity downwards
			if (snowlayers[lyr].water > snowlayers[lyr].whc) {
				ts_melt = snowlayers[lyr].water - snowlayers[lyr].whc;
				snowlayers[lyr].water = snowlayers[lyr].whc;
			}
		}
	}

	// Update number of snow layers
	nsnowlayers = snowlayers.size();

	// Add excess water this timesteps melt to surface melt-water to pass to soil on next timestep 
	melt += ts_melt;

	double final_snow = melt;

	// Update snow depth and calculate total snow
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].depth =  snowlayers[lyr].ice / (snowlayers[lyr].dens / rho_H2O);
		final_snow += snowlayers[lyr].ice + snowlayers[lyr].water;
	}

	const double EPS = 1.0e-8;
	assert(abs(initial_snow - final_snow) < EPS);
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
	double snow_volume;	// Volume of snow layer [m-3]

	// Minimum air volume that will influence snow layer (only affect top snow layer)
	double min_snow_depth = snow_depth_standard * (snow_depth_limit - 1.0) * 10.0;

	// Going through snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// Snow layer depth
		snow_depth = lyr == nsnowlayers - 1 ? max(min_snow_depth, snowlayers[lyr].depth) : snowlayers[lyr].depth;

		// Water and ice fraction of full volume
		double Frac_water = snowlayers[lyr].water / snow_depth;
		double Frac_ice = snowlayers[lyr].ice / snow_depth;

		double init_water = snowlayers[lyr].water;
		double init_ice = snowlayers[lyr].ice;
		double init_temp = snowlayers[lyr].temp;

		// Volume of snow layer [m-3]
		snow_volume = snowlayers[lyr].depth * M_PER_MM;

		// The total amount of ice and water should not change in this routine
		double wtotal_initial = snowlayers[lyr].ice + snowlayers[lyr].water;

		// Temperature change from yesterday
		delta_T = snowlayers[lyr].temp - snowlayers[lyr].temp_old;

		// Freezing
		if (snowlayers[lyr].temp < 0.0 && snowlayers[lyr].temp < snowlayers[lyr].temp_old) { // Wania conditions 
			Fthaw = 0.0;

			// Water to freeze
			if (Frac_water > 0.0) {

				// Available energy to freeze water[J m-3]
				energy = Cp_water * fabs(delta_T);

				// Fraction of a volume of water that could freeze [-]
				Ffreez = energy / Lheat;

				// More water in a volume than could freeze
				if (Frac_water > Ffreez) {

					// There is some water left, even after freeze
					Frac_ice += Ffreez;
					Frac_water -= Ffreez;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
				}
				else {

					// It is cold enough to freeze all the water and
					// decrease the temperature below 0, so only the 
					// energy need for phase change is removed from 
					// the already calculated temperature 

					// All water freezes
					Frac_ice += Frac_water;

					// Rescale Fwater if snow layer depth is lower than min_snow_depth for energy calculation
					Frac_water *= snow_depth / snowlayers[lyr].depth;

					// Calculate new snow temperature
					snowlayers[lyr].temp += Frac_water * snow_volume * Lheat / snowlayers[lyr].Ch;

					// No liquid water left in snow layer
					Frac_water = 0.0;
				}
			}
			else {
				Ffreez = 0.0;
			}
		} // end of freezing

		  // Thawing  
		else if (snowlayers[lyr].temp >= 0.0  && snowlayers[lyr].temp > snowlayers[lyr].temp_old) { // Wania conditions

			Ffreez = 0.0;

			// Ice to thaw
			if (Frac_ice > 0.0) {

				// Energy to thaw ice [J m-3]
				energy = Cp_water * fabs(delta_T);

				// Fraction of a volume of ice that would thaw [-]
				Fthaw = energy / Lheat;

				// More ice in a volume than could thaw
				if (Frac_ice > Fthaw) {

					// There is some ice left, even after melt
					Frac_water += Fthaw;
					Frac_ice -= Fthaw;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
				}
				else {
					// There is more energy available than required for
					// melting the rest of the ice. Therefore no need to 
					// calculate new temperature as snow layer will be removed

					// All ice thaws
					Frac_water += Frac_ice;

					// No ice left in snow layer
					Frac_ice = 0.0;
				}
			}
			else {
				Fthaw = 0.0;
			}
		} // end of thawing

		  // New ice fraction in snow layer
		snowlayers[lyr].ice = Frac_ice * snow_depth;

		// New water fraction in snow layer
		snowlayers[lyr].water = Frac_water * snow_depth;

		// The total amount of ice and water should not change in this routine
		double wtotal_final = snowlayers[lyr].ice + snowlayers[lyr].water;

		//dprintf("%d %d %d AirT %g Temp %g %g Ice %g %g Water %g %g Tot %g %g\n", date.year, date.day, lyr, airT, snowlayers[lyr].temp, init_temp, snowlayers[lyr].ice, init_ice, snowlayers[lyr].water, init_water, wtotal_initial, wtotal_final);
		double diff = abs(wtotal_initial - wtotal_final);

		const double EPS = 1.0e-5;

		assert(abs(wtotal_initial - wtotal_final) < EPS);

		// Update T_old after possible phase change
		snowlayers[lyr].temp_old = snowlayers[lyr].temp;

	} // end of loop through layers
}


/// Update snow density
void Snow::update_density() {

	// Calculate snow depth ice content before compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// snow depth based only on ice
		snowlayers[lyr].depth = snowlayers[lyr].ice / (snowlayers[lyr].dens / rho_H2O);

		// snow density, only defined by ice fraction, Oleason 2013, Eq. 6.88. can not exceed maximum density
		snowlayers[lyr].dens = min(rho_ICE, snowlayers[lyr].ice / (snowlayers[lyr].depth * M_PER_MM));
	}

	// Mechanical compaction scheme, following Best et al. 2011
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// Current layers ice and water mass
		double I_sum = snowlayers[lyr].ice + snowlayers[lyr].water;

		// Snow mass above this layer
		double I_above = 0.0;
		for (int ii = lyr + 1; ii < nsnowlayers; ii++) {
			I_above += snowlayers[ii - 1].ice + snowlayers[ii - 1].water;
		}

		// Cumulative snow mass affecting compaction
		double M_k = 0.5 * I_sum + I_above;

		// Layer temperature in Kelvin
		double T_k = max(K2degC, snowlayers[lyr].temp + K2degC);

		// Increase in layer density due to compaction (Best et. al. 2011, Eq.21)
		double k_s = 4000.0;
		double eta_0 = 1E7; // Compactive viscosity factor for snow
		double rho_0 = 50.0; // Reference snow density [kg m-3]
		double delta_dens = snowlayers[lyr].dens * gravity * M_k / eta_0 *
			exp((k_s / K2degC) - (k_s / T_k) - (snowlayers[lyr].dens / rho_0)) * SECS_PER_DAY;
		snowlayers[lyr].dens += delta_dens;

		// Snow density can not exceed maximum density
		snowlayers[lyr].dens = min(rho_ICE, snowlayers[lyr].dens);
	}

	// Calculate snow depth after compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].depth = snowlayers[lyr].ice / (snowlayers[lyr].dens / rho_H2O);
	}

	// Update snow water holding capacity
	update_whc();
}

/// Accumulate snow in bottom layer if all other layers are already full
void Snow::update_layers(double prec, double snowfall_dens) {

	// Check initial snow and water content
	double initial_snow = prec * snowfall_dens;
	double initial_water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_snow += snowlayers[lyr].ice * snowlayers[lyr].dens;
		initial_water += snowlayers[lyr].water;
	}

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		double ice_new, ice_old, ice_sum, ice_frac_move, water_move;

		// Total snow to move down to this layer
		double tot_snow_move = lyr == 0 ? prec : snow_depth_standard - snowlayers[lyr].ice;

		// Layer to move snow from
		int lyr_move = lyr + 1;

		// Move snow until enough is moved (could be from more than one layer above)
		while (!negligible(tot_snow_move)) {

			double snow_move;

			// Determine snow to move from snow layer above

			// Snow layer above exist
			if (lyr_move < nsnowlayers) {
				snow_move = min(tot_snow_move, snowlayers[lyr_move].ice);

				// Move water as well
				ice_frac_move = !negligible(snowlayers[lyr_move].ice) ? snow_move / snowlayers[lyr_move].ice : 0.0;
				water_move = snowlayers[lyr_move].water * ice_frac_move;
				snowlayers[lyr_move].water -= water_move;
				snowlayers[lyr].water += water_move;
			}
			else { // No layer above so take from prec to fill layer
				snow_move = min(tot_snow_move, snow_depth_standard - snowlayers[lyr].ice);
			}

			// Calculate old and new snow fractions
			ice_sum = snowlayers[lyr].ice + snow_move;
			ice_old = !negligible(ice_sum) ? snowlayers[lyr].ice / ice_sum : 0.0;
			ice_new = 1.0 - ice_old;

			snowlayers[lyr].ice += snow_move;

			// Update snow density
			if (lyr_move < nsnowlayers) {
				snowlayers[lyr].dens = ice_old * snowlayers[lyr].dens + ice_new * snowlayers[lyr_move].dens;
				snowlayers[lyr_move].ice -= snow_move;
			}
			else {
				snowlayers[lyr].dens = ice_old * snowlayers[lyr].dens + ice_new * snowfall_dens;
			}

			// Remove added snow from this layer from the toal that needs to be moved
			tot_snow_move -= snow_move;

			// Go to the next layer above (if tot_snow_move > 0)
			lyr_move++;
		}
	}

	// Check final snow and water content
	double final_snow = 0.0;
	double final_water = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_snow += snowlayers[lyr].ice * snowlayers[lyr].dens;
		final_water += snowlayers[lyr].water;
	}

	const double EPS = 1.0e-8;
	assert(abs(initial_snow - final_snow) < EPS);
	assert(abs(initial_water - final_water) < EPS);
}

/// Snow accumulation
void Snow::snowfall(double airT, double prec) {

	// If snowfall
	if (airT < 0.0 && prec > 0.0) {

		// Fresh snow density calculation, following eqn 1 in Vionnet et al (2012)
		int a = 109; // kg m-3, 
		int b = 6;  // kg m-3 K-1
		int c = 26;  // kg m-7/2 s -1/2
		double U = 1.0;	// reference wind speed [m s-1]
		double dens_min = 50; // minimum snow density 
		double Tmelt = 0; // temperature of the melting point for water

		double snowfall_dens = max(dens_min, a + b * (airT - Tmelt) + c * pow(U, 0.5));

		// If no snow layers
		if (nsnowlayers == 0) {
			snowlayers.resize(nsnowlayers + 1);
			snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT);
			nsnowlayers++;
		}

		int top_snow_lyr = nsnowlayers - 1;

		// Fill snow layers from top snow layer
		while (!negligible(prec) && top_snow_lyr < NLAYERS_SNOW_MAX) {

			// If new snowfall doesn't make snow layer exceed maximum snow layer depth, then add it to layer 
			if ((top_snow_lyr < NLAYERS_SNOW_MAX - 1 && (snowlayers[top_snow_lyr].ice + prec) <= snow_depth_standard * snow_depth_limit) ||
				(top_snow_lyr == NLAYERS_SNOW_MAX - 1 && (snowlayers[top_snow_lyr].ice + prec) <= snow_depth_standard)) {

				// update dens
				snowlayers[top_snow_lyr].mix_dens(prec, snowfall_dens);

				// update temperature
				snowlayers[top_snow_lyr].temp = snowlayers[top_snow_lyr].temp * (snowlayers[top_snow_lyr].ice / (snowlayers[top_snow_lyr].ice + prec)) +
					airT * (prec / (snowlayers[top_snow_lyr].ice + prec));

				// update ice
				snowlayers[top_snow_lyr].ice += prec;
				prec = 0.0;
			}
			else {

				// Snowfall will make snow layer exceed maximum snow layer depth

				// Reduce snow layer dept as it is already higher than snow_depth_standard
				// and then add prec to next layer (NLAYERS_SNOW_MAX layer can't get above snow_depth_standard)
				if (snowlayers[top_snow_lyr].ice > snow_depth_standard) {

					double ice_move = snowlayers[top_snow_lyr].ice - snow_depth_standard;
					double frac_move = ice_move / snowlayers[top_snow_lyr].ice;

					snowlayers[top_snow_lyr].ice = snow_depth_standard;

					// Need to move water as well
					double water_move = snowlayers[top_snow_lyr].water * frac_move;
					snowlayers[top_snow_lyr].water -= water_move;

					// Create new snow layer with same density as original snow layer
					snowlayers.resize(nsnowlayers + 1);
					snowlayers[nsnowlayers] = Snowlayer(ice_move, water_move, snowlayers[top_snow_lyr].dens, snowlayers[top_snow_lyr].temp);
					nsnowlayers++;
				}
				else {
					// Fill up to standard_depth and then put rest in next layer
					double prec_lyr = snow_depth_standard - snowlayers[top_snow_lyr].ice;

					// update dens
					snowlayers[top_snow_lyr].mix_dens(prec_lyr, snowfall_dens);

					// update ice
					snowlayers[top_snow_lyr].ice += prec_lyr;
					prec -= prec_lyr;

					// Create a new empty snow layer
					if (nsnowlayers < NLAYERS_SNOW_MAX) {
						snowlayers.resize(nsnowlayers + 1);
						snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT);
						nsnowlayers++;
					}
				}

				top_snow_lyr++;
			}

			// snow left when last layer full -> special case
			if (prec > 0.0 && top_snow_lyr == NLAYERS_SNOW_MAX) {
				update_layers(prec, snowfall_dens);
			}
		}
	}
}

void Snow::serialize(ArchiveStream& arch) {
	arch & snowlayers
		& nsnowlayers
		& snow_depth_standard
		& snow_depth_limit
		& dec_depth
		& NH4_mass
		& NO3_mass;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soillayer member functions
////////////////////////////////////////////////////////////////////////////////

Soillayer::Soillayer(double d) : depth(d) {

	// Initialises certain member variables


	alag = 0.0;
	exp_alag = 1.0;
	cpool_slow = 0.0;
	cpool_fast = 0.0;
	decomp_litter_mean = 0.0;
	k_soilfast_mean = 0.0;
	k_soilslow_mean = 0.0;

	orgleachfrac = 0.0;
	aorgCleach = 0.0;
	aorgNleach = 0.0;
	aminleach = 0.0;

	wcont = 0.0;
	wpcont = 0.0;
	perc_frac = 0.0;
	orgleachfrac = 0.0;
	orgmovefrac = 0.0;
	Dw = 0.0;
	Kw = 0.0;
	runoff_drain = 0.0;
	runoff_baseflow = 0.0;
	water = 0.0;
	ice = 0.0;
	air = 0.0;
	for (int mth = 0; mth < 12; mth++) {
		mwcont[mth] = 0.0;
		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		morgmove_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
		mminmove_mean[mth] = 0.0;
	}
	temp = 0.0;
	Kh = 0.0;
	Ch = 0.0;

	std::fill_n(dwcont, Date::MAX_YEAR_LENGTH, 0.0);

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

	NO2_mass = 0.0;
	NO2_mass_w = 0.0;
	NO2_mass_d = 0.0;
	NO_mass = 0.0;
	NO_mass_w = 0.0;
	NO_mass_d = 0.0;
	N2O_mass = 0.0;
	N2O_mass_w = 0.0;
	N2O_mass_d = 0.0;
	N2_mass = 0.0;
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	anmin = 0.0;
	animmob = 0.0;
	aminleach = 0.0;
	aorgNleach = 0.0;
	aorgCleach = 0.0;

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
	O2 = 0.0;
	volume_liquid_water = 0.0;
	total_volume_water = 0.0;
	tiller_area = 0.0001;
	dcflux_soil = 0.0;

	Frac_water = 0.0;
	Frac_org = 0.0;
	Frac_peat = 0.0;
	Frac_min = 0.0;
	Fpwp_ref = 0.0;
	Frac_water_belowpwp = 0.0;
	Frac_ice = 0.0;
	Frac_ice_yesterday = 0.0;
	Frac_air = 1.0;
	rootfrac = 0.0;
	por = 0.0;
	awc_peat = 0.0;

	pH = 6.5;

	labile_carbon = 0.0;
	labile_carbon_w = 0.0;
	labile_carbon_d = 0.0;
}

// Set soil layer properties
void Soillayer::set_properties(Soiltype soiltype) {

	awc = soiltype.awc * depth;
	wp = soiltype.wp * depth;
	wsats = soiltype.wsats * depth;
	por = soiltype.porosity * depth;
	organic_frac = soiltype.organic_frac;
	mineral_frac = soiltype.mineral_frac;
	
	// soiltype variable for Medium - fine // DW_COMMENT needs to be changed
	alpha = 0.83; // DW_COMMENT should come from input
	Kw_max = 2.63E-07; // DW_COMMENT should come from input
	m = 1.0 - 1.0 / 1.2539; // DW_COMMENT should come from input
	lambda = -0.5884; // DW_COMMENT should come from input
	/////////////////////////
	
	clay_frac = soiltype.clay_frac;
	sand_frac = soiltype.sand_frac;
	silt_frac = soiltype.silt_frac;
}

/// C content of soil layer
double Soillayer::ccont() {
	double ccont = 0.0;

	ccont += cpool_fast;
	ccont += cpool_slow;

	for (int i = 0; i < NSOMPOOL; i++) {
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
double Soillayer::wcontent() {

	return wcont * awc + wpcont * wp;
}

double Soillayer::nmass_avail(int pref, bool water_limited) {

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
	return nmass * (water_limited ? min(1.0, wcont) : 1.0);
}

void Soillayer::nmass_subtract(double nmass, int pref) {

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
			//patch.fluxes.report_flux(Fluxes::NH3_SOIL, r); // DW_COMMENT
			dprintf("Year %d Day %d NH4 mass %f N substr %f \n", date.year ,date.day, NH4_mass, nmass);
			//fail("tried to subtract more NH4 than available");
		}
	}
	else if (pref == NO3) {
		if (nmass<NO3_mass || negligible(nmass - NO3_mass, -10)) {
			NO3_mass -= nmass;
		}
		else {
			double r = nmass - NO3_mass;
			NO3_mass = 0.0;
			//patch.fluxes.report_flux(Fluxes::NO_SOIL, r); // DW_COMMENT
			dprintf("Year %d Day %d NO3 mass %f N substr %f \n", date.year, date.day, NO3_mass, nmass);
			//fail("tried to subtract more NO3 than available");
		}
	}
}

void Soillayer::nmass_inc(double nmass, int pref) {
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

void Soillayer::nmass_multiplic_inc(double inc, int pref) {
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



void Soillayer::nmass_remove_small_values(int pref) {
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

void Soillayer::update_moisture_diff_and_cond() {

	const double W_res = 0.01;

	if ((wsats - ice) / depth <= 1.001 * W_res) {
		Dw = 0.0;
		Kw = 0.0;
		return;
	}

	double W_sat, W_tot;

	// Account for ice filling up space
	W_sat = max(1.001 * W_res, wsats / depth);
	W_tot = (awc * wcont + wp * wpcont) / depth;

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
	Dw = ((1.0 - m) * Kw_max / (alpha * m * (W_sat - W_res))) * pow(the, (lambda - 1.0 / m)) * (pow(1.0 - po, m) + pow(1.0 - po, -m) - 2);

	// Conductivity (m s-1; van Genuchten 1980 eqn 8, adopted from HTESSEL code)
	if (W_tot > 1.001 * W_res && W_tot <= 0.999 * W_sat) {
		Kw = Kw_max * pow(the, lambda) * pow(1.0 - pow(1.0 - pow(the, 1.0 / m), m), 2);
	}
}

void Soillayer::update_heatcap_and_cond() {

	// DESCRIPTION
	// Called daily to update heat capacities and thermal conductivities of the 
	// soil layers below the mixed/moss and snow layers 

	double lKfwater;			// log of therm. cond. for the water fraction
	double lKfice;				// log of therm. cond. for the ice fraction
	double lKforg;				// log of therm. cond. for the org. fraction
	double lKfpeat;				// log of therm. cond. for the peat fraction
	double lKfmin;				// log of therm. cond. for the min. fraction
	double Ftot;				// total of solid fraction in soil

								// DW_COMMENT peat will be included when merge with 4.1
	double peat = 0.0;

	// Update air (mm)
	air = por - (water + ice);

	if (air < -1.0e-8)
		fail("Illegal air mass value in update_heat_diff_and_cond(): %g", air);

	// Catch rounding errors for air
	air = max(0.0, air);

	Frac_air = air / depth;

	// Update heat capacity (J m-3 K-1)
	Ch = mineral_frac * Cp_min + organic_frac * Cp_org +
		(ice * Cp_ice + water * Cp_water + peat * Cp_peat + air * Cp_air) / depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	Ftot = mineral_frac + organic_frac + (ice + water + peat) / depth;

	if (Ftot + air / depth > 1.000001)
		fail("Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g", Ftot, air / depth, Ftot + air / depth);

	lKforg = (organic_frac / Ftot) * log(Korg);
	lKfpeat = (peat / depth / Ftot) * log(Kpeat);
	lKfmin = (mineral_frac / Ftot) * log(Kmin);
	lKfwater = (water / depth / Ftot) * log(Kwater);
	lKfice = (ice / depth / Ftot) * log(Kice);

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh = air / depth * Kair + Ftot * exp(lKfwater + lKfice + lKforg + lKfmin + lKfpeat);
}

// Get the Water Filled Pore Space of the layer
double Soillayer::wfps() const {
	return (wcont * awc + wpcont * wp) / wsats;
}

double Soillayer::wcont_max() {
	return (wsats - wp) / awc;
}

void Soillayer::serialize(ArchiveStream& arch) {
	arch & depth
		& wcont
		& wpcont
		& awcont
		& dwcont
		& temp
		& temp_old
		& temp_analytic
		& dtemp
		& mtemp
		& gtemp
		& mwcont
		& awc
		& wp
		& wsats;

	for (int i = 0; i<NSOMPOOL; i++) {
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
		& mminmove_mean

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

	// Set soil layer depths (needs to be the same in init_rootdist()!!!)
	std::vector<double> soil_layer_depths = patch.stand.get_gridcell().soil_layer_depths();

	nsoillayers = soil_layer_depths.size();
	nsoillayersupper = 4; // DW_COMMENT check so it is right!

	// Set number of soil layers
	soillayers.resize(nsoillayers);

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		soillayers[lyr] = Soillayer(soil_layer_depths[lyr]);
	}

	// Set sizes of water flow calculation parameters
	wpara.set_sizes(nsoillayers);

	// Set sizes of temperature change calculation parameters
	tpara.set_sizes(nsoillayers + PAD_LAYERS);

	// Padding layers
	pad.resize(PAD_LAYERS);

	// Set padding layers depth
	set_padding_depth();

	wfluxes.resize(nsoillayers + 1, 0.0);
	b = set_b(300.0); // standard deviation of orography (300) 
	NH4_input = 0.0;
	NO3_input = 0.0;
	anfix = 0.0;
	anfix_mean = 0.0;
	anfix_calc = 0.0;
	awcont_upper = 0.0;
	wcont_evap = 0.0;
	rain_melt = 0.0;
	max_rain_melt = 0.0;
	dperc = 0.0;
	percolate = 0.0;
	runoff = 0.0;

	// Peatland hydrology variables
	awtp = 0.0;
	acro_depth = 0.0;
	cato_depth = 0.0;
	acro_co2 = 934.0; // [mimol L-1]
	wtd = 0.0; // [-100, +300] mm
	acro_por = acrotelm_por - Fgas; // DW_COMMENT MAYBE SHOULD BE IN SOILLAYER
	cato_por = catotelm_por - Fgas;
	Wtot = 0.0;
	stand_water = 0.0;
	ch4_store = 0.0;
	co2_store = 0.0;

	for (int lyr = 0; lyr < NSUBLAYERS_ACRO; lyr++)
		sub_water[lyr] = 0.0; // DW_COMMENT what to do???

	dmoss_wtp_limit = 1.0;
	dgraminoid_wtp_limit = 1.0;

	k_O2 = 0.0;
	k_CO2 = 0.0;
	k_CH4 = 0.0;
	Ceq_O2 = 0.0;
	Ceq_CO2 = 0.0;
	Ceq_CH4 = 0.0;

	firstTempCalc = true;
	firstHydrologyCalc = true;
	num_evaplayers = 2; // DW_COMMENT should be calculated

	// Depths of the acrotelm & catotelm 
	acro_depth = NACROTELM * Dz_acro;
	cato_depth = soil_depth - acro_depth;

	// Solve LPJ SOM pools
	solvesom_end = static_cast<int>(0.8 * nyear_spinup);
	solvesom_begin = static_cast<int>(0.7 * nyear_spinup);

	// Solve CENTURY SOM pools
	solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	solvesomcent_endyr = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);
	savesomcent = false;
}

std::vector<double> Soilcolumn::get_wcont() {
	std::vector<double> wconts(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wconts[lyr] = soillayers[lyr].wcont;
	return wconts;
}

std::vector<double> Soilcolumn::get_awc() {
	std::vector<double> awcs(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		awcs[lyr] = soillayers[lyr].awc;
	return awcs;
}

double Soilcolumn::nmass_avail(int pref, bool water_limited) {

	double nmass_avail = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		nmass_avail += soillayers[lyr].nmass_avail(pref, water_limited);
	}
	return nmass_avail;
}

void Soilcolumn::nmass_avail_uptake(double uptake, int pref) {
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
		runoff_drainage += soillayers[lyr].runoff_drain;
	}
	return runoff_drainage;
}

double Soilcolumn::runoff_baseflow() {
	double runoff_base_flow = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		runoff_base_flow += soillayers[lyr].runoff_baseflow;
	}
	return runoff_base_flow;
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
double Soilcolumn::wcont() {
	double wcont = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wcont += soillayers[lyr].wcontent();
	}
	return wcont;
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
		soillayers[lyr].dwcont[date.day] = soillayers[lyr].wcont;
	}
}

// Set daily soil layer temperature for this month
void Soilcolumn::set_dtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].dtemp[date.dayofmonth] = soillayers[lyr].temp;
	}
}

void Soilcolumn::set_mwcont() {
	double acc_wcont = 0.0;
	double acc_depth = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].mwcont[date.month] = mean(soillayers[lyr].dwcont + date.day - date.ndaymonth[date.month] + 1,
			date.ndaymonth[date.month]);

		acc_depth += lyr == 0 ? soillayers[lyr].depth / 2.0 : (soillayers[lyr - 1].depth + soillayers[lyr].depth) / 2.0;
		acc_wcont += lyr == 0 ? soillayers[lyr].mwcont[date.month] * soillayers[lyr].depth / 2.0 : (soillayers[lyr-1].mwcont[date.month] * soillayers[lyr-1].depth + soillayers[lyr].mwcont[date.month] * soillayers[lyr].depth) / 2.0;
		soillayers[lyr].mwcontabove = acc_wcont / acc_depth;
	}
}

// Set this months average soil layer temperature
void Soilcolumn::set_mtemp() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].mtemp[date.month] = mean(soillayers[lyr].dtemp, date.ndaymonth[date.month]);
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

			soillayers[lyr].temp_analytic = climate.temp;
		}
	}
	else {

		double depth = 0;

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			depth += lyr == 0 ? soillayers[lyr].depth / 2.0 : (soillayers[lyr-1].depth + soillayers[lyr].depth) / 2.0;

			if (date.islastday) {

				// Linearly interpolate soil thermal diffusivity given mean
				// soil water content from this layer and above



				if (soillayers[lyr].mwcontabove < 0.15)
					k = ((soiltype.thermdiff_15 - soiltype.thermdiff_0) / 0.15 *
						soillayers[lyr].mwcontabove + soiltype.thermdiff_0) * DIFFUS_CONV;
				else
					k = ((soiltype.thermdiff_100 - soiltype.thermdiff_15) / 0.85 *
					(soillayers[lyr].mwcontabove - 0.15) + soiltype.thermdiff_15) * DIFFUS_CONV;

				// Calculate parameters alag and exp(-alag) from Eqn 2

				soillayers[lyr].alag = depth / sqrt(k / HALF_OMEGA); // from Eqn 1
				soillayers[lyr].exp_alag = exp(-soillayers[lyr].alag);

			}

			// Every day, calculate linear model for trend in daily air
			// temperatures for the last 31 days: temp_day = a + b * day

			double buffer[31];
			climate.dtemp_31.to_array(buffer);
			regress(day, buffer, 31, a, b);

			// Calculate soil temperature

			temp_lag = a + b * (30.0 - soillayers[lyr].alag * LAG_CONV);

			// Calculate temperatures the old way
			soillayers[lyr].temp_analytic = climate.atemp_mean + soillayers[lyr].exp_alag*(temp_lag - climate.atemp_mean);
		}
	}
}


/// Average temperature down to a depth
double Soilcolumn::temp_depth(double depth) {

	// Calculate temperature for layers down to depth
	double acc_depth = 0.0;
	double acc_temp = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_temp += soillayers[lyr].temp * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return !negligible(acc_depth) ? acc_temp / acc_depth : 0.0;
}

// Average current wcont down to a depth
double Soilcolumn::wcont_depth(double depth) {

	// Calculate wcont for layers down to depth
	double acc_depth = 0.0;
	double acc_wcont = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_wcont += soillayers[lyr].wcont * layer_depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return !negligible(acc_depth) ? acc_wcont / acc_depth : 0.0;
}

// Average monthly wcont down to a depth
double Soilcolumn::mwcont_depth(double depth, int mth) {

	// Calculate wcont for layers down to depth
	double acc_depth = 0.0;
	double acc_mwcont = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_mwcont += soillayers[lyr].mwcont[mth] * layer_depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return !negligible(acc_depth) ? acc_mwcont / acc_depth : 0.0;
}

// Average yearly wcont down to a depth
double Soilcolumn::awcont_depth(double depth) {

	// Calculate wcont for layers down to depth
	double acc_depth = 0.0;
	double acc_awcont = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_awcont += soillayers[lyr].awcont * layer_depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return !negligible(acc_depth) ? acc_awcont / acc_depth : 0.0;
}

// Average current wpcont down to a depth
double Soilcolumn::wpcont_depth(double depth) {

	// Calculate wcont for layers down to depth
	double acc_depth = 0.0;
	double acc_wpcont = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_wpcont += soillayers[lyr].wpcont * layer_depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return !negligible(acc_depth) ? acc_wpcont / acc_depth : 0.0;
}

double Soilcolumn::wsats_depth(double depth) {

	// Calculate awc for layers down to depth
	double acc_depth = 0.0;
	double acc_wsats = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_wsats += soillayers[lyr].wsats * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_wsats;
}

/*double Soilcolumn::water_depth(double depth) { // DW_COMMENT not used

	// Calculate awc for layers down to depth
	double acc_depth = 0.0;
	double acc_water = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_water += soillayers[lyr].water * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_water;
}*/

double Soilcolumn::awc_depth(double depth) {

	// Calculate awc for layers down to depth
	double acc_depth = 0.0;
	double acc_awc = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_awc += soillayers[lyr].awc * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_awc;
}

double Soilcolumn::awc_peat_depth(double depth) {

	// Calculate awc for layers down to depth
	double acc_depth = 0.0;
	double acc_awc = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_awc += soillayers[lyr].awc_peat * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_awc;
}

double Soilcolumn::wp_depth(double depth) {

	// Calculate wp for layers down to depth
	double acc_depth = 0.0;
	double acc_wp = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_wp += soillayers[lyr].wp * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_wp;
}

double Soilcolumn::wcont_max_depth(double depth) {

	// Calculate wp for layers down to depth
	double acc_depth = 0.0;
	double acc_wcont_max = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_wcont_max += soillayers[lyr].wcont_max() * layer_depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_wcont_max / acc_depth;
}

double Soilcolumn::NH4_depth(double depth) {

	// Calculate wp for layers down to depth
	double acc_depth = 0.0;
	double acc_NH4 = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_NH4 += soillayers[lyr].NH4_mass * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_NH4;
}

double Soilcolumn::NO3_depth(double depth) {

	// Calculate wp for layers down to depth
	double acc_depth = 0.0;
	double acc_NO3 = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_NO3 += soillayers[lyr].NO3_mass * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_NO3;
}

double Soilcolumn::ice_depth(double depth) {

	// Calculate ice (mm) for layers down to depth
	double acc_depth = 0.0;
	double acc_ice = 0.0;
	int lyr = 0;
	while (acc_depth < depth && lyr < nsoillayers) {
		double layer_depth = acc_depth + soillayers[lyr].depth < depth ? soillayers[lyr].depth : depth - acc_depth;
		acc_ice += soillayers[lyr].ice * layer_depth / soillayers[lyr].depth;
		acc_depth += layer_depth;
		lyr++;
	}
	return acc_ice;
}

// add water to a certain soil layer and return the water that cannot be added due to limited capacity
double Soilcolumn::add_layer_soil_water(int layer, double extrawater) {

	Soillayer& soillayer = soillayers[layer];

	double overflow = 0.0;
	double Faw_layer = 0.0; // mm
	double ice_layer = 0.0; // mm
	double layerwater_max = 0.0; // mm

	if (patch.stand.is_highlatitude_peatland_stand()) {
		Faw_layer = soillayer.wcont * soillayer.awc_peat; // mm
		ice_layer = soillayer.Frac_ice * soillayer.depth; // mm
		layerwater_max = (acro_por - peat_wp) * soillayer.depth;

	}
	else {
		Faw_layer = soillayer.wcont * soillayer.awc; // mm
		ice_layer = soillayer.Frac_ice * soillayer.depth; // mm
		layerwater_max = soillayer.wsats;
	}

	double layerwater = Faw_layer + ice_layer; // mm
	double potential_layer = layerwater_max - layerwater;

	if (extrawater < potential_layer)
		overflow = 0.0;
	else {
		overflow = extrawater - potential_layer;
		extrawater = potential_layer;
	}

	// Now add as much water to this layer as it can hold
	Faw_layer += extrawater;
	potential_layer -= extrawater;

	// wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
	if (patch.stand.is_highlatitude_peatland_stand())
		soillayer.wcont = Faw_layer / soillayer.awc_peat;
	else
		soillayer.wcont = Faw_layer / soillayer.awc;

	oob_check_wcont(soillayer.wcont, soillayer.wcont_max());
	soillayer.Frac_water = Faw_layer / soillayer.depth;

	return overflow;
}

void Soilcolumn::update_moisture_diff_and_cond() {
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_moisture_diff_and_cond();
	}
}

// Update heat capacity and conductivity
void Soilcolumn::update_heatcap_and_cond() {

	// Update heat capacity and conductivity for each snow layer
	snow.update_heatcap_and_cond();

	// Update heat capacity and conductivity for each soil layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].update_heatcap_and_cond();
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

	int n = main.size();
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

void Soilcolumn::calc_wflow(double infil) {

	bool ifbedrock = false; // 1 (Bedrock) or 0 (free drainage)

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wpara.psi0[lyr] = (soillayers[lyr].wcont * soillayers[lyr].awc + soillayers[lyr].wpcont * soillayers[lyr].wp) / soillayers[lyr].depth;
	}

	for (int ts = 0; ts < date.timesteps; ts++) {

		// Calculate moisture diffusivity and conductivity
		update_moisture_diff_and_cond();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			wpara.dz[lyr] = soillayers[lyr].depth * M_PER_MM;
			wpara.diff[lyr] = soillayers[lyr].Dw;
			wpara.cond[lyr] = soillayers[lyr].Kw;
			wpara.psi[lyr] = (soillayers[lyr].wcont * soillayers[lyr].awc + soillayers[lyr].wpcont * soillayers[lyr].wp) / soillayers[lyr].depth;
		}

		double depth = wpara.dz[0] / 2.0;
		wpara.z[0] = depth;
		for (int lyr = 1; lyr < nsoillayers; lyr++) {
			depth += (wpara.dz[lyr - 1] + wpara.dz[lyr]) / 2.0;
			wpara.z[lyr] = depth;
		}

		// Boundary conditions for soil water
		// Upper boundary condition
		double f_top = infil / rho_H2O / date.dt;

		// Lower boundary condition
		double f_bot = ifbedrock ? 0.0 : wpara.cond[nsoillayers - 1];

		// Build diagonals and RHS vector
		build_diags_rhs(date.dt / date.timesteps, f_top, f_bot,
			wpara.cap, wpara.diff, wpara.cond,
			wpara.upper, wpara.main, wpara.lower, wpara.rhs,
			wpara.z, wpara.dz, wpara.psi, wpara.sink, nsoillayers);

		// Solve tridiagonal system
		solve_tridiagonal(wpara.upper, wpara.main, wpara.lower, wpara.rhs, wpara.psi1);

		double water_end = 0.0;
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			water_end += wpara.psi1[lyr] * soillayers[lyr].depth;
		}

		// Set new wcont
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			if (wpara.psi1[lyr] * soillayers[lyr].depth > soillayers[lyr].wp) {
				soillayers[lyr].wpcont = 1.0;
				soillayers[lyr].wcont = (wpara.psi1[lyr] * soillayers[lyr].depth - soillayers[lyr].wp) / soillayers[lyr].awc;  // DW_COMMENT PEAT awc_peat MAYBE NOT ALREADY TAKEN CARE OF
			}
			else {
				soillayers[lyr].wpcont = wpara.psi1[lyr] * soillayers[lyr].depth / soillayers[lyr].wp;
				soillayers[lyr].wcont = 0.0;
			}

			// Calculate drainage

			// Total water
			double wtot = soillayers[lyr].wcont * soillayers[lyr].awc + soillayers[lyr].wpcont * soillayers[lyr].wp;

			// Drainage of water above saturation including ice
			double drained_water = max(0.0, wtot + soillayers[lyr].ice - soillayers[lyr].wsats);

			soillayers[lyr].runoff_drain += drained_water;

			// Update total water after drainage 
			// wtot and ice will now max fill up soil layer to saturation
			wtot -= drained_water;

			//if (drained_water > 0.0)
			//	dprintf("%d %d %d %g %g\n",date.year,date.day,lyr, drained_water, soillayers[lyr].temp);

			// Water will first fill up to wilting point
			soillayers[lyr].wpcont = min(1.0, wtot / soillayers[lyr].wp);

			// Left over water fills up above wilting point
			soillayers[lyr].wcont = (wtot - soillayers[lyr].wpcont * soillayers[lyr].wp) / soillayers[lyr].awc;

			// Remove small values of wcont and wpcont // DW_COMMENT what is considerred small when it comes to wcont and wpcont
			if (soillayers[lyr].wcont < 1.0e-12) soillayers[lyr].wcont = 0.0;
			if (soillayers[lyr].wpcont < 1.0e-12) soillayers[lyr].wpcont = 0.0;

			soillayers[lyr].water = soillayers[lyr].wcont * soillayers[lyr].awc + soillayers[lyr].wpcont * soillayers[lyr].wp;
		}
	}

	// Calculate water fluxes to, inbetween and out of soil layers
	wfluxes[0] = infil;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		wfluxes[lyr + 1] = wfluxes[lyr] + wpara.psi0[lyr] * soillayers[lyr].depth - (soillayers[lyr].water + soillayers[lyr].runoff_drain);

		// Remove very small water flux values
		if (negligible(abs(wfluxes[lyr + 1]), -10))
			wfluxes[lyr + 1] = 0.0;
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
	double dz = soillayers[nsoillayers - 1].depth;

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

void Soilcolumn::update_layer_fractions() {

	// DESCRIPTION
	// Inializes and updates fractions of water, mineral content, peat content, organic content, as well as porosity.

	if (patch.stand.is_highlatitude_peatland_stand()) {

		if (firstTempCalc) {

			for (int lyr = 0; lyr < NACROTELM; lyr++) {

				// ACROTELM

				soillayers[lyr].Frac_por = acro_por;
				soillayers[lyr].Frac_org = 0.0;
				soillayers[lyr].Frac_min = 0.0;
				soillayers[lyr].Frac_peat = 1.0 - (acro_por + Fgas);

				soillayers[lyr].Fpwp_ref = peat_wp;  // Wisser et al. 2011 Very low in peatlands
				soillayers[lyr].Frac_water_belowpwp = peat_wp;
				soillayers[lyr].Frac_water = acro_por - peat_wp; // i.e. initially saturated

				soillayers[lyr].awc_peat = (soillayers[lyr].Frac_por - peat_wp) * soillayers[lyr].depth;
			}

			for (int lyr = NACROTELM; lyr < nsoillayers; lyr++) {

				// CATOTELM

				soillayers[lyr].Frac_por = cato_por; // as in RW code
				soillayers[lyr].Frac_org = 0.0;
				soillayers[lyr].Frac_min = 0.0;
				soillayers[lyr].Frac_peat = 1.0 - (cato_por + Fgas);

				soillayers[lyr].Fpwp_ref = peat_wp;  // Wisser et al. 2011 Very low in peatlands
				soillayers[lyr].Frac_water_belowpwp = peat_wp;
				soillayers[lyr].Frac_water = cato_por - peat_wp; // i.e. initially saturated

				soillayers[lyr].awc_peat = (soillayers[lyr].Frac_por - peat_wp) * soillayers[lyr].depth;
			}
		} // firstTempCalc
	}
	else {

		// NON-WETLAND SOILS

		if (firstTempCalc && patch.stand.first_year == date.year) {

			for (int lyr = 0; lyr < nsoillayers; lyr++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 7) {

					// No updates to the standard values for LPJ-GUESS soils
					soillayers[lyr].Frac_por = soiltype.porosity;
					soillayers[lyr].Frac_org = soiltype.organic_frac;
					soillayers[lyr].Frac_min = soiltype.mineral_frac;
					soillayers[lyr].Fpwp_ref = soiltype.wp;
					soillayers[lyr].Frac_water_belowpwp = soiltype.wp;
				}
				else {
					// Use values updated using input of soil carbon
					soillayers[lyr].Frac_por = soiltype.porosity_gridcell[lyr];
					soillayers[lyr].Frac_org = soiltype.org_frac_gridcell[lyr];
					soillayers[lyr].Frac_min = soiltype.min_frac_gridcell[lyr];
					soillayers[lyr].Fpwp_ref = soiltype.wp;
					soillayers[lyr].Frac_water_belowpwp = soiltype.wp; // Water content initially set to the wilting point
				}

				soillayers[lyr].Frac_water = 0.0;
				soillayers[lyr].Frac_peat = 0.0;
			}
		}
		else if (firstTempCalc && (restart || patch.stand.clone_year == date.year)) {
			// initialize after restarts
			for (int lyr = 0; lyr < nsoillayers; lyr++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 7) {

					// No updates to the standard values for LPJ-GUESS soils
					soillayers[lyr].Frac_por = soiltype.porosity;
					soillayers[lyr].Frac_org = soiltype.organic_frac;
					soillayers[lyr].Frac_min = soiltype.mineral_frac;
					soillayers[lyr].Fpwp_ref = soiltype.wp;
				}
				else {
					// Use values updated using input of soil carbon
					soillayers[lyr].Frac_por = soiltype.porosity_gridcell[lyr];
					soillayers[lyr].Frac_org = soiltype.org_frac_gridcell[lyr];
					soillayers[lyr].Frac_min = soiltype.min_frac_gridcell[lyr];
					soillayers[lyr].Fpwp_ref = soiltype.wp;
				}
				soillayers[lyr].Frac_peat = 0.0;
			}
		}// firstTempCalc
	}
}
bool Soilcolumn::update_layer_water_content() {

	// DESCRIPTION
	// Private helper method to update the soil water contents in the soil layers. 
	// Called each day.

	// whc[] is water above the pwp which hasn't frozen [mm]
	// Its maximum value is awc[] (mm), i.e. whc[later] <= awc[layer], for all layers

	const double maxerr = 0.000001; // Max error allowed

	// Local storage
	std::vector <double> Faw_layer(nsoillayers);
	// available water for each soil layer (mm)
	std::vector <double> ice_layer(nsoillayers);
	// mm of ice in each soil layer (mm)

	// Peatland total soil water content, incl ice and water below the pwp
	Wtot = 0.0; // mm

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Update the total available volumetric water content (mm) of the layers
		Faw_layer[lyr] = soillayers[lyr].Frac_water * soillayers[lyr].depth; // mm
		ice_layer[lyr] = soillayers[lyr].Frac_ice * soillayers[lyr].depth; // mm

		// Check balance
		double layerwater = Faw_layer[lyr] + ice_layer[lyr]; // mm

		if (soillayers[lyr].aw_max - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_layer_water_content - error in a soil layer's water balance!\n");
			return false;
		}

		// Subtract the ice/frozen fraction from the initial value of the capacity for this layer to hold water
		soillayers[lyr].alwhc = soillayers[lyr].alwhc_init - soillayers[lyr].Frac_ice;

		// For wetlands:
		// total ice in layer = Frac_ice[i] + Fpwp_ref[i] - Frac_water_belowpwp[i];
		// total water in layer = Frac_water[i] + Frac_water_belowpwp[i];
		// So, the total ice + water = Frac_ice[i] + Fpwp_ref[i] + Frac_water[i] 
		if (lyr < NACROTELM)
			Wtot += (soillayers[lyr].Frac_water + soillayers[lyr].Frac_ice + soillayers[lyr].Fpwp_ref) * soillayers[lyr].depth;

		if (patch.stand.is_highlatitude_peatland_stand())
			soillayers[lyr].whc = max(0.0, soillayers[lyr].aw_max - ice_layer[lyr]); // mm 
		else
			soillayers[lyr].whc = max(0.0, soillayers[lyr].awc - ice_layer[lyr]); // mm 

		// Avoid dividing by small numbers later
		if (fabs(soillayers[lyr].whc) < maxerr)
			soillayers[lyr].whc = 0.0;

		// wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			soillayers[lyr].wcont = Faw_layer[lyr] / soillayers[lyr].awc_peat;
		else
			soillayers[lyr].wcont = Faw_layer[lyr] / soillayers[lyr].awc;

		oob_check_wcont(soillayers[lyr].wcont, soillayers[lyr].wcont_max());

		if (soillayers[lyr].wcont < 0.0 || soillayers[lyr].wcont > 1.0) { // DW_COMMENT allows wcont > 1
			dprintf("Soil::update_layer_water_content - error in the wcont calculation!\n");
			return false;
		}
	} // for loop (lyr)

	  // wcont_evap is the available, nonfrozen water in the top two layers (20 cm), divided by the capacity of those layers to hold water
	wcont_evap = wcont_depth(SOILDEPTH_EVAP);

	oob_check_wcont(wcont_evap, wcont_max_depth(SOILDEPTH_EVAP));

	if (wcont_evap < 0.0 || wcont_evap > 1.0) {
		dprintf("Soil::update_layer_water_content - error in the wcont_evap calculation!\n");
		return false;
	}

	return true; // ...if no problems
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

	acro_co2 = min(PORE_WATER_CO2, PORE_WATER_CO2 + ((PORE_WATER_CO2 - atmo_co2) / acro_depth) * awtp);
}

void Soilcolumn::update_soil_water() {

	// DESCRIPTION
	// Update Frac_water and wcont_evap when there has been a change to wcont

	const double maxerr = 0.000001; // Max error allowed

	// Local storage
	// available water for each soil layer (mm)
	std::vector<double> Faw_layer(nsoillayers);
	// ice each soil layer (mm)
	std::vector<double> ice_layer(nsoillayers);

	Wtot = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		oob_check_wcont(soillayers[lyr].wcont, soillayers[lyr].wcont_max());

		if (soillayers[lyr].wcont < 0.0 || soillayers[lyr].wcont > soillayers[lyr].wcont_max()) {
			dprintf("Soil::update_soil_water - bad wcont! %g\n", soillayers[lyr].wcont);
			return;
		}

		// For the following, note that, for each layer: wcont = Faw_layer / soiltype.awc
		// Also, wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			Faw_layer[lyr] = soillayers[lyr].wcont * soillayers[lyr].awc_peat; // mm // DW_COMMENT DET AWC AND ALL TO PEAT IF IT IS A PEAT STAND
		else
			Faw_layer[lyr] = soillayers[lyr].wcont * soillayers[lyr].awc; // mm

		ice_layer[lyr] = soillayers[lyr].Frac_ice * soillayers[lyr].depth; // mm

		// Check balance
		double layerwater = Faw_layer[lyr] + ice_layer[lyr]; // mm

		if (soillayers[lyr].aw_max - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_soil_water - error in a soil layer's water balance!\n");
			return;
		}

		// Now update Frac_water
		soillayers[lyr].Frac_water = Faw_layer[lyr] / soillayers[lyr].depth;

		// For wetlands:
		if (lyr < NACROTELM)
			Wtot += (soillayers[lyr].Frac_water + soillayers[lyr].Frac_ice + soillayers[lyr].Fpwp_ref) * soillayers[lyr].depth;

	} // for loop (lyr)

	// Update wcont_evap
	// wcont_evap is the available, nonfrozen water in the top two layers (20 cm), divided by the capacity of those layers to hold water
	if (patch.stand.is_highlatitude_peatland_stand())
		wcont_evap = wcont_depth(SOILDEPTH_EVAP);
	else
		wcont_evap = wcont_depth(SOILDEPTH_EVAP);

	if (wcont_evap < 0.0 || wcont_evap > wcont_max_depth(SOILDEPTH_EVAP)) {
		dprintf("Soil::update_soil_water - error in the wcont_evap calculation! Year %d Day %d wcont_evap %g max %g\n", date.year, date.day, wcont_evap, wcont_max_depth(SOILDEPTH_EVAP));
		return;
	}

	oob_check_wcont(wcont_evap, wcont_max_depth(SOILDEPTH_EVAP));

	return; // ...if no problems
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


	// PARAMETERS and variables for the calculation of the WTP
	// See: Granberg et al. (1999), Eqns. 1-5.

	/// depth of each sphagnum layer [mm]
	const double Dz_sph = 10.0;
	/// surface water content when WTD is at 20 mm
	const double surfw = 0.34;
	/// minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// gradient limit of upper gradient [mm]
	const double grad_lim = 50.0;
	/// Granberg et al. (1999) parameter, eqn 2.
	double az;
	/// Depth of sublayers in the acrotelm
	double Dz_sub;

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	// Porosity updates in the presence of ice? 
	bool subtractIceFromWtot = false;
	double acro_por_icy = acro_por; // the pore space minus the ice fraction
	double avgIceFrac = 0.0;

	if (subtractIceFromWtot) {
		for (int lyr = 0; lyr < NACROTELM; lyr++)
			avgIceFrac += (soillayers[lyr].Frac_ice + soillayers[lyr].Fpwp_ref - soillayers[lyr].Frac_water_belowpwp) / (double)NACROTELM;
	}

	acro_por_icy -= avgIceFrac;


	// *** TRANSPIRATION *** 

	// Retrieve Vegetation for this patch
	Vegetation& vegetation = patch.vegetation;

	// Transpiration from the acrotelm  
	double aet_acrotelm = 0.0;
	double aet = 0.0; // AET for a particular layer and individual (mm)

	// Sum AET for across all vegetation individuals in this patch

	 // Calculate fpc_moss
	double fpc_moss_total = 0.0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		for (int lyr = 0; lyr < NACROTELM; lyr++)
			pft_acro_root_frac += patch.pft[indiv.pft.id].pft.rootdist[lyr];

		// this includes mosses, which have root_frac = 1 in the acrotelm
		// Fraction of the aet from the acrotelm only.
		aet_acrotelm += indiv.aet * pft_acro_root_frac;

		// Total aet, assumed from acrotelm and catotelm combined 
		aet_total += indiv.aet;

		if (indiv.pft.ismoss())
			fpc_moss_total += indiv.fpc;

		vegetation.nextobj();
	}


	// *** EVAPORATION *** 
	// Evaporation from soil surface

	// Wania et al (2009)
	if (snow.depth() < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
		evap = fevap * climate.eet * PRIESTLEY_TAYLOR * (0.99 / (1.0 + exp(-1.0 * (-wtd + 98.7) / 22.6) + 0.02));
	}
	else {
		evap = 0.0;
	}

	// *** EVAPOTRANSPIRATION *** 

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
		// evapotranspiration=eet*exp[0.02*wtp_cm], where // wtp_cm [100,-30]
		// This would give almost identical results to the Kim & Verma formula below 
		// when the water table is at the surface

		if (snow.depth() < 10.0) { // i.e. evap only if snow depth < 1cm
			evapotranspiration = climate.eet * (1.02 - 0.00075 * wtd) * fevap; // wtd [-100,+300]
		}
		else {
			evapotranspiration = 0.0;
		}

	}
	else {

		// This is JUST evapotranspiration from the acrotelm
		evapotranspiration = aet_acrotelm + evap;
	}

	// *** ACROTELM WATER VOLUME ***

	az = (acro_por - minvtot) / zmin; // Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22

	// Drainage - could possibly be read in for site-specific studies
	runoff_drain = 0; // DW_COMMENT CHECK HYDROLOGY PEAT AND DRAINAGE

	// Update available water 
	Wtot += rain_melt - evapotranspiration - runoff_drain;

	// *** RUNOFF AND RUNON***

	double acrowater = 0.0; // Frac_water[MIDX] * Dz[MIDX];

	for (int lyr = 0; lyr < NACROTELM; lyr++)
		acrowater += (soillayers[lyr].Frac_water + soillayers[lyr].Frac_water_belowpwp) * soillayers[lyr].depth;

	double ideal_runoff = exp(-0.01 * wtd);
	if (acrowater > 0.0 && soillayers[0].Frac_ice < 0.7) // Granberg et al. (1999) use the same ice condition - see f_icestop in their Table 1
		runoff_surf = min(ideal_runoff, acrowater);
	else
		runoff_surf = 0.0;

	// .ins file option: 
	// run on or off - but only when there's flowing water

	double runon = 0.0;

	// runon set in global.ins
	soiltype.runon = wetland_runon;

	if (runoff_surf > 0.0) {

		double runofforon = soiltype.runon;

		if (runofforon < 0)
			runoff_surf -= runofforon;	// e.g. 3 mm/day for bog-like/hummock conditions

		if (runofforon > 0)
			runon += runofforon;		// 3-5 mm/day for fen-like conditions
	}

	double Wtot_init = Wtot;

	Wtot += runon - runoff_surf;

	// Add standing water, if any
	Wtot += stand_water;

	// Not in LPJ-WHy
	if (Wtot < 0.0) {
		Wtot = 0.0;
		runoff_surf = Wtot_init;
	}

	// Save evapotranspiration data
	if (empirical_evapotranspiration) {

		patch.aaet += (1.0 - fevap) * evapotranspiration;
		patch.aevap += fevap * evapotranspiration;
		patch.maet[date.month] += (1.0 - fevap) * evapotranspiration;
		patch.mevap[date.month] += fevap * evapotranspiration;
	}
	else {

		patch.aevap += evap;
		patch.mevap[date.month] += evap;
		patch.aaet += aet_total;
		patch.maet[date.month] += aet_total;
	}


	// *** UPDATE WATER TABLE DEPTH ***

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

	// NB 
	// e.g. Wtot = 270mm =>
	// Frac_water = 270/300 = 0.9 in the three acrotelm layers

	// Wtot <= 140.03611111 =>
	// Frac_water = 0.27419, 0.4186, 0.7075 in the top three layers (if no ice)

	if (Wtot > acro_depth * acro_por) {

		// Saturated acrotelm

		wtp[date.day] = Wtot - acro_depth * acro_por;
		// i.e. wtp[date.day] > 0 above the surface

		if (wtp[date.day] > maxh) {
			runoff_surf += wtp[date.day] - maxh;
			// Update runoff from the surface
			wtp[date.day] = maxh;
		}

		wtd = -wtp[date.day]; // i.e. wtd defined as < 0 above the surface
		Wtot = acro_depth * acro_por;
		// I.e. Wtot does NOT include standing water
		stand_water = 0.0; // = wtp[date.day]; // we assume no standing water 
	}
	else {

		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (acro_por * acro_depth - Wtot) / (2.0 * az));

		if (wtd > zmin)
			wtd = 3.0 * (acro_por * acro_depth - Wtot) / (2.0 * (acro_por - minvtot));

		if (wtd > acro_depth)
			wtd = acro_depth;

		wtp[date.day] = -wtd;
		stand_water = 0.0;
	}


	// *** RUNOFF ***

	// Save today's runoff in Soil::runoff
	runoff = runoff_surf;
	patch.arunoff += runoff;
	patch.mrunoff[date.month] += runoff;
	patch.asurfrunoff += runoff;


	// Depth of acrotelm sublayers - usually 10 mm
	Dz_sub = NACROTELM * Dz_acro / NSUBLAYERS_ACRO;


	// *** WATER CONTENT IN EACH LAYER - THE PROFILE ***

	double surfw_sat = 0.0;
	double value[NSUBLAYERS_ACRO];
	double acro_aw = 0.0; // mm of plant available water in acrotelm
	double cato_aw = 0.0; // mm of plant available water in catotelm

	if (wtd > 0.0) {

		// Below the surface

		// Determine the subsurface layer (0 to NSUBLAYERS_ACRO) within which wtd lies. 
		// Possible values:
		// = 0 if wtd < 10 mm 
		// = NSUBLAYERS_ACRO iff wtd = Dz_acro
		int wtd_layer = (int)(wtd / Dz_sub); // Takes the integer part only

		surfw_sat = max(minvtot, acro_por - az * wtd); // Granberg - eqn 3.

		if (wtd_layer == 0) {

			// Total saturation throughout the acrotelm (from 10mm down)
			for (int j = 0; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por;
			}
		}
		else {

			// Unsaturated zone in the acrotelm - use Granberg's eqn 1.
			for (int i = 0; i < wtd_layer; i++) {
				double z = (i + 1.0) * Dz_sub; // mm
				value[i] = min(acro_por, surfw_sat + (acro_por_icy - surfw_sat) * pow(z / wtd, 2));
			}

			// Saturated zone in the acrotelm
			for (int j = wtd_layer; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por;
			}
		}

		// Now calculate the volumetric water content throughout the acrotelm
		sub_water[0] = 0.5 * (surfw_sat + value[0]);
		for (int lyr = 1; lyr < NSUBLAYERS_ACRO; lyr++) {
			sub_water[lyr] = 0.5 * (value[lyr] + value[lyr - 1]);
		}

		// Water content in each acrotelm layer
		double Ftotal[NACROTELM]; // Volumetric fraction in the acrotelm layers
		for (int lyr = 0; lyr < NACROTELM; lyr++)
			Ftotal[lyr] = 0.0;

		int sublayers_per_layer = (int)Dz_acro / (int)Dz_sub;

		int count = 0;
		for (int lyr = 0; lyr < NACROTELM; lyr++) {

			double layerwater = 0.0;
			for (int subly = count * sublayers_per_layer; subly < (count + 1) * sublayers_per_layer; subly++) {
				layerwater += sub_water[subly] * Dz_sub; // add mm of water for this sublayer 
			}

			// Divide by the depth of the acrotelm layers to get the volumetric fraction
			Ftotal[lyr] = layerwater / Dz_acro;
			count++;
		}

		// Subtract ice content, water content below the permanent wilting point

		// ACROTELM

		for (int lyr = 0; lyr < NACROTELM; lyr++) {
			soillayers[lyr].Frac_water = max(Ftotal[lyr] - soillayers[lyr].Frac_ice - soillayers[lyr].Fpwp_ref, 0.0);
			acro_aw += soillayers[lyr].Frac_water * soillayers[lyr].depth;
		}
	}
	else {

		// WTD <= 0.0
		// at or above the surface - standing water

		// Subtract ice content, water content below the permanent wilting point
		// For the calculations below, note that:
		// totaliceinlayer = Frac_ice + Fpwp_ref - Frac_water_belowpwp;
		// totalwaterinlayer = Frac_water + Frac_water_belowpwp;

		// ACROTELM - saturated
		for (int lyr = 0; lyr < NACROTELM; lyr++) {
			soillayers[lyr].Frac_water = max(soillayers[lyr].Frac_por - soillayers[lyr].Frac_ice - soillayers[lyr].Fpwp_ref, 0.0);
			acro_aw += soillayers[lyr].Frac_water * soillayers[lyr].depth;
		}
	}

	// CATOTELM - assumed to be always saturated
	for (int lyr = NACROTELM; lyr < nsoillayers; lyr++) {
		soillayers[lyr].Frac_water = max(soillayers[lyr].Frac_por - soillayers[lyr].Frac_ice- soillayers[lyr].Fpwp_ref, 0.0);
		cato_aw += soillayers[lyr].Frac_water * soillayers[lyr].depth;
	}

	// Tidy up Frac_water
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		if (soillayers[lyr].Frac_water < soillayers[lyr].Fpwp_ref)
			soillayers[lyr].Frac_water = 0.0;
	}

	// *** WATER CONTENT IN EACH SOIL LATER (FOR PLANT GROWTH) ***
	if (!update_layer_water_content())
		return;

	// UPDATE CO2 LEVEL IN ACROTELM
	update_acrotelm_co2(climate.co2);


	// *** UPDATE LIMITS ON MOSS PHOTOSYNTHESIS ***

	// The cap is kept at 100% while WTP is above wtp_moss_upper ((Wania et al 2009b) have -150 mm), decreases linearly
	// between wtp_moss_upper and -280mm  This has to do with the water
	// content of the mosses as they are less able to photosynthesize
	// as they dry out. NO MOSS DESSICATION

	double wtp_moss_upper = 0.0;	// Mosses can be too competitive when set to -150 mm, as in Wania et al. 2009b 
	double wtp_moss_lower = -280.0; // Wania et al. 2009b
	double moss_wtp_limit_lowerlimit = 0.3;

	if (wtp[date.day] > wtp_moss_upper)
		dmoss_wtp_limit = 1.0; // No limit
	else if (wtp[date.day] > wtp_moss_lower)
		dmoss_wtp_limit = 1.0 - (-wtp[date.day] + wtp_moss_upper) * (1.0 - moss_wtp_limit_lowerlimit) / (wtp_moss_upper - wtp_moss_lower);
	else
		dmoss_wtp_limit = moss_wtp_limit_lowerlimit; // between -300mm and -280mm


	// *** UPDATE LIMITS ON GRAMINOID PHOTOSYNTHESIS ***

	// LPJ-WHy code & comment:

	// IF (pft .EQ. 11 .AND. mwtp(m) .LT. -100) THEN
	// The following lines will reduce the mgpp of
	// flood-tolerant grasses when the water table drops
	// below -10cm. This should reflect that these grasses
	// are usually found where the water table stays close
	// to the surface for most of the time (i.e. in fens).
	// Note that mwtp is negative!

	// mgpp(m,pft) = mgpp(m,pft) + mgpp(m,pft) *
	// $                (mwtp(m) + 100) / 200
	// ENDIF

	double wtp_graminoid_upper = -10; // -50.0;		// mm - Wania et al. 2009b have -100
	double wtp_graminoid_lower = -100; // -200.0;	// mm - Wania et al. 2009b have -300 mm
	double graminoid_wtp_limit_lowerlimit = 0.0;

	// More restrictive than Wania et al., otherwise graminoids become too dominant as they never suffer from inundation stress
	// The cap drops from 1 as soon as wtp goes below wtp_graminoid_upper
	if (wtp[date.day] >= wtp_graminoid_upper)
		dgraminoid_wtp_limit = 1.0; // i.e. no limit!
	else if (wtp[date.day] > wtp_graminoid_lower)
		dgraminoid_wtp_limit = 1.0 - (-wtp[date.day] + wtp_graminoid_upper) * (1.0 - graminoid_wtp_limit_lowerlimit) / (wtp_graminoid_upper - wtp_graminoid_lower);
	else
		dgraminoid_wtp_limit = graminoid_wtp_limit_lowerlimit; // between -300mm and -200mm


	// *** UPDATE OUTPUT VARIABLES ***

	// Reset monthly wtp and mmoss_wtp_limit averages on Jan 1
	if (date.day == 0) {
		for (int mth = 0; mth < 12; mth++) {
			mwtp[mth] = 0.0;
		}
	}

	int mdays = date.ndaymonth[date.month];
	mwtp[date.month] += wtp[date.day] / (double)mdays;

	if (firstHydrologyCalc)
		firstHydrologyCalc = false;

	if (date.day == Date::MAX_YEAR_LENGTH) {

		awtp = 0.0;
		for (int d = 0; d < 365; d++)
			awtp += wtp[d] / 365.0;
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
	// wcont      = array containing water content of soil layers [0=upper layer, 14=bottom layer] as
	//              fraction of available water holding capacity (AWC)
	//				This method ensures that layers 2-4 have identical wcont values, as do layers 5-14
	// wcont_evap = water content of evaporation sublayers (20 cm) at top of soil column
	//              as fraction of available water holding capacity (AWC)
	//				This method ensures that layers 0-1 have identical wcont_evap values
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
	// awc        = array containing available water holding capacity of soil
	//              layers (mm rainfall) [0=upper layer]
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface
	// snowpack   = depth of snow (mm)

	// INPUT AND OUTPUT PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont_upper = wcont averaged over the growing season - guess2008
	// dperc      = daily percolation beyond system (mm)

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	const double SOILDEPTH_EVAP = 200.0;
	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double BASEFLOW_FRAC = 0.5;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double K_DEPTH = 0.4;
	const double K_AET = 0.52;
	// Fraction of total (vegetation) AET from upper soil layer that is derived
	// from the top K_DEPTH (fraction) of the upper soil layer
	// (parameters for calculating K_AET_DEPTH below)
	const double K_AET_DEPTH = (soillayers[0].depth / SOILDEPTH_EVAP - 1.0) *
		(K_AET / K_DEPTH - 1.0) / (1.0 / K_DEPTH - 1.0) + 1.0;
	// Weighting coefficient for AET flux from evaporation layer, assuming active
	//   root density decreases with soil depth
	// Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
	//   K_DEPTH=0.4, K_AET=0.52

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	double aet;										// AET for a particular layer and individual (mm)
	std::vector<double> aet_lyr(nsoillayers);	// total AET for each Gerten-type soil layer (mm)
	double aet_total = 0.0;

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int s = 0; s<nsoillayers; s++) {
			aet = patch.pft[indiv.pft.id].fwuptake[s] * indiv.aet;
			aet_total += aet;
			aet_lyr[s] += aet;
		}
		vegetation.nextobj();
	}

	// Get the two-layer versions of the wcont values
	std::vector<double> wcont_lyr(nsoillayers);	// wcont for each Gerten-type soil layer
	wcont_lyr = get_wcont();

	// Evaporation from soil surface

	// guess2008 - changed to wcont_evap**2, as in LPJ-mL
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	double evap = 0.0;
	if (snow.depth() < 10.0) {	// evap only if snow depth < 10mm
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}

	// Implement in- and outgoing fluxes to upper soil layer
	// BLARP: water content can become negative, though apparently only very slightly
	//    - quick fix implemented here, should be done better later

	wcont_lyr[0] += (rain_melt - aet_lyr[0] - evap) / soillayers[0].awc;
	if (wcont_lyr[0] != 0.0 && wcont_lyr[0] < 0.0001) { // guess2008 - bugfix
		wcont_lyr[0] = 0.0;
	}

	// Surface runoff
	double runoff_surf = 0.0;
	if (wcont_lyr[0] > 1.0) {
		runoff_surf = (wcont_lyr[0] - 1.0) * soillayers[0].awc;
		wcont_lyr[0] = 1.0;
	}

	// Update water content in evaporation layer for tomorrow

	wcont_evap += (rain_melt - aet_lyr[0] * SOILDEPTH_EVAP * K_AET_DEPTH / soillayers[0].depth - evap)
		/ soillayers[0].awc;

	wcont_evap = max(0.0, min(wcont_lyr[0], wcont_evap));

	// Percolation from evaporation layer
	double perc = 0.0;
	if (percolate) {
		perc = min(SOILDEPTH_EVAP / soillayers[0].depth * soiltype.perc_base * pow(wcont_evap, soiltype.perc_exp),
			max_rain_melt);
	}
	wcont_evap -= perc / soillayers[0].awc;

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	double runoff_drain = 0.0;

	for (int lyr = 1; lyr < nsoillayers; lyr++) {

		// Percolation
		// Allow only on days with rain or snowmelt (Dieter Gerten, 021216)

		if (percolate) {
			perc = min(soiltype.perc_base * pow(wcont_lyr[lyr - 1], soiltype.perc_exp), max_rain_melt);
		}
		else {
			perc = 0.0;
		}
		soillayers[lyr - 1].perc_frac = min(perc / soillayers[lyr - 1].awc, wcont_lyr[lyr - 1]);

		wcont_lyr[lyr - 1] -= soillayers[lyr - 1].perc_frac;
		wcont_lyr[lyr] += soillayers[lyr - 1].perc_frac * soillayers[lyr - 1].awc / soillayers[lyr].awc;
		if (wcont_lyr[lyr] > 1.0) {
			runoff_drain += (wcont_lyr[lyr] - 1.0) * soillayers[lyr].awc;
			wcont_lyr[lyr] = 1.0;
		}

		// Deduct AET from this soil layer
		// BLARP! Quick fix here to prevent negative soil water

		wcont_lyr[lyr] -= min(wcont_lyr[lyr], aet_lyr[lyr] / soillayers[lyr].awc);
	}

	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_baseflow = 0.0;
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC * soiltype.perc_base * pow(wcont_lyr[nsoillayers - 1], soiltype.perc_exp);
		// guess2008 - Added "&& rain_melt >= runoff_surf" to guarantee nonnegative baseflow.
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer

		soillayers[nsoillayers - 1].perc_frac = min(perc_baseflow / soillayers[nsoillayers - 1].awc, wcont_lyr[nsoillayers - 1]);
		wcont_lyr[nsoillayers - 1] -= soillayers[nsoillayers - 1].perc_frac;
		runoff_baseflow = soillayers[nsoillayers - 1].perc_frac * soillayers[nsoillayers - 1].awc;
	}

	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;


	// update awcont_upper
	// Original algorithm by Thomas Hickler

	// Reset awcont_upper on the first day of every year
	if (date.day == 0) {
		awcont_upper = 0.0;
		patch.growingseasondays = 0;
	}

	// If it's warm enough for growth, update awcont_upper with this day's wcont
	if (climate.temp > 5.0) {
		awcont_upper += wcont_lyr[0];
		patch.growingseasondays++;
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}


	// Use the two-layer values calculated above to update the soil properties needed in all NSOILLAYER layers

	// Update wcont
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wcont = wcont_lyr[lyr];
	}

	// Finally, update wcont_evap, Frac_water etc. in this patch.soil object
	update_soil_water();
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
void Soilcolumn::hydrology_lpjf(const Climate& climate, double fevap) { // DW_COMMENT peat soils!!!!!!!!!! and other things that exist in v4.1 hydrology ICEEEEEEEE

	// Daily update of water content for each soil layer given snow melt, rainfall,
	// evapotranspiration from vegetation (AET), evaporation and percolation between
	// layers and runoff

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// INPUT AND OUTPUT PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont     = wcont averaged over the growing season - guess2008

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	/// Soil depth for soil evaporation
	const double SOILDEPTH_EVAP = 200.0;
	/// Soil depth for soil infiltration
	const double SOILDEPTH_INFILT = 500.0;

	if (date.year == 5 && date.day == 32)
		int sch = 0;

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	double wcont_begin = wcont() + rain_melt;

	// AET for a particular layer and individual (mm)
	double aet;
	// total AET for each soil layer (mm)
	std::vector<double> aet_layer(nsoillayers);
	// total AET
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			double fwu = patch.pft[indiv.pft.id].fwuptake[lyr];
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
			aet_layer[lyr] += aet;
			aet_total += aet;
		}
		vegetation.nextobj();
	}

	// Evaporation

	// Calculate wcont and awc for evaporation layers (only up to awc)
	double wcont_evap = min(1.0, wcont_depth(SOILDEPTH_EVAP));

	// Calculate total evaporation
	// guess2008 - changed to wcont_evap**2, as in LPJ-mL
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	// evap only if snow depth < 10mm
	double evap = 0.0;
	if (snow.depth() < 10.0) {
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}

	// Save evaportion in a temporary variable to check if enough water exist to evaporate
	double evap_check = evap;

	// Extract evap first from incoming precipitation and snow melt
	if (rain_melt >= evap_check) {
		rain_melt -= evap_check;
		evap_check = 0.0;
	}
	else {
		evap_check -= rain_melt;
		rain_melt = 0.0;
	}

	// Evaporation fraction of available water down to evaporation depth after aet is taken
	double acc_depth = 0.0;
	int lyr = 0;
	double avail_evap_water = 0.0;
	double evap_frac = 0.0;
	if (!negligible(evap_check)) {
		while (acc_depth < SOILDEPTH_EVAP) {
			avail_evap_water += soillayers[lyr].wcont * soillayers[lyr].awc - aet_layer[lyr];
			acc_depth += soillayers[lyr].depth;
			lyr++;
		}
		if (avail_evap_water >= evap_check) {
			evap_frac = evap_check / avail_evap_water;
			evap_check = 0.0;
		}
		else {
			evap_check -= avail_evap_water;
			evap_frac = 1.0;
		}
	}

	// Update evaporation if precipitation, snow melt and soil water couldn't meet evaporation demand
	evap -= evap_check;

	// Extract evaporation and aet from soil column
	acc_depth = 0.0;
	for (lyr = 0; lyr < nsoillayers; lyr++) {

		// Evaporation included in this sink layer
		if (acc_depth < SOILDEPTH_EVAP) {
			soillayers[lyr].wcont -= (aet_layer[lyr] + (soillayers[lyr].wcont * soillayers[lyr].awc - aet_layer[lyr]) * evap_frac) / soillayers[lyr].awc;
		}
		// Below evaporation depth so only AET in sink term
		else {
			soillayers[lyr].wcont -= aet_layer[lyr] / soillayers[lyr].awc;
		}

		// Remove small values of WCONT
		if (soillayers[lyr].wcont < 1.0e-12) {
			soillayers[lyr].wcont = 0.0;
		}

		// Accomolated soil depth to know if we still are within evaporation depth
		acc_depth += soillayers[lyr].depth;

		// Reset runoff_darin for calc_wflow
		soillayers[lyr].runoff_drain = 0.0;
	}

	// DW_COMMENT CHECK SO THAT AET AND EVAP TAKES THE RIGHT AMOUNT OF WATER

	// Surface infiltration and runoff

	// Calculate water and ice content for infiltration layers
	double wcont_infilt = wcont_depth(SOILDEPTH_INFILT);
	double wpcont_infilt = wpcont_depth(SOILDEPTH_INFILT);
	double wsats_infilt = wsats_depth(SOILDEPTH_INFILT);
	double ice_infilt = ice_depth(SOILDEPTH_INFILT);
	double awc_infilt = awc_depth(SOILDEPTH_INFILT);
	double wp_infilt = wp_depth(SOILDEPTH_INFILT);
	double totwcont_infilt = (wcont_infilt * awc_infilt + wpcont_infilt * wp_infilt + ice_infilt) / wsats_infilt;

	// Calculate runoff if there is any throughfall or snow melt (Balsamo 2009, eqn 6 and 11)
	double runoff_surf = 0.0;
	if (!negligible(rain_melt)) {
		runoff_surf = max(0.0, rain_melt - max(0.0, wsats_infilt * (1.0 - totwcont_infilt)) + wsats_infilt *
			pow(max(0.0, pow(1.0 - totwcont_infilt, 1.0 / (1.0 + b)) -
				rain_melt / (wsats_infilt * (1.0 + b))), b + 1.0));
	}

	// Calculate infiltration (mm)
	double infilt = max(0.0, rain_melt - runoff_surf);

	// Calculate water profile and fluxes
	calc_wflow(infilt);

	// Calculate baseflow
	soillayers[nsoillayers - 1].runoff_baseflow = max(0.0, wfluxes[nsoillayers]); // DW_COMMENT should be resolved if negative!!!

	if (wfluxes[nsoillayers] < 0.0)
		dprintf("Year %d Day %d Baseflow %g\n",date.year,date.day, wfluxes[nsoillayers]);

	// Total runoff
	double runoff = runoff_surf + runoff_drain() + runoff_baseflow();

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain();
	patch.abaserunoff += runoff_baseflow();
	patch.arunoff += runoff;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;

	// Calculate annual wcont for growing season (temperature > 5C)
	// guess2008 - DLE - update awcont
	// Original algorithm by Thomas Hickler
	for (lyr = 0; lyr < nsoillayers; lyr++) {

		// Reset the awcont array on the first day of every year
		if (date.day == 0) {
			soillayers[lyr].awcont = 0.0;
			if (lyr == 0) patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			soillayers[lyr].awcont += soillayers[lyr].wcont;
			if (lyr == 0) patch.growingseasondays++;
		}

		// Do the averaging on the last day of every year
		if (date.islastday && date.islastmonth) {
			// In case it's never warm enough:
			if (patch.growingseasondays < 1) {
				soillayers[lyr].awcont = 1.0;
			}
			else {
				soillayers[lyr].awcont /= (double)patch.growingseasondays;
			}
		}
	}

	double wcont_end = wcont() + runoff + evap + aet_total;

	double EPS = 1.0e-5;
	
	if (fabs(wcont_begin - wcont_end) > EPS)
		dprintf("Year %d Day %d WCONT ERROR: Diff %g Begin %g End %g\n", date.year, date.day, wcont_begin - wcont_end, wcont_begin, wcont_end);

	// Finally, recalculate Frac_water based on these calculations and updates to wcont[] 
	update_soil_water();
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

	double Ffreez;	// fraction of water which freezes
	double Fthaw;	// fraction of ice which thaws
	double energy;	// Energy per volume [J m-3]

	const double SMALLVOLFRAC = 0.0000001;	// Minimum volume fraction [m-3] allowed

	double delta_T;	// temperature difference (current - previous)

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Initial water and ice content
		double wtotal_initial = soillayers[lyr].water + soillayers[lyr].ice;

		// Water and ice fraction of full volume
		double Fwater = soillayers[lyr].water / soillayers[lyr].depth;
		double Fice = soillayers[lyr].ice / soillayers[lyr].depth;

		// temperature difference
		delta_T = soillayers[lyr].temp - soillayers[lyr].temp_old;

		// FREEZING
		if (soillayers[lyr].temp < 0.0 && soillayers[lyr].temp < soillayers[lyr].temp_old) { // Wania conditions 
			Fthaw = 0.0;
			if (Fwater > 0.0) {
				energy = Cp_water * fabs(delta_T); // energy in J m-3 
				Ffreez = energy / Lheat; // unitless fraction

				if (Fwater > Ffreez) {
					Fice += Ffreez;
					Fwater -= Ffreez;

					soillayers[lyr].temp = 0.0;
				}
				else {

					// It is cold enough to freeze all the water and
					// decrease the temperature below 0, since the rest
					// of the water gets frozen now.

					Fice += Fwater;
					soillayers[lyr].temp += Fwater * Lheat / soillayers[lyr].Ch;
					soillayers[lyr].temp = min(0.0, soillayers[lyr].temp);

					Fwater = 0.0;
				}
			}
			else {
				Ffreez = 0.0;
			}
		} // end of FREEZING

		  // THAWING
		else if (soillayers[lyr].temp >= 0.0  && soillayers[lyr].temp > soillayers[lyr].temp_old) { // Wania conditions
			Ffreez = 0.0;
			if (Fice > 0.0) {
				energy = Cp_water * fabs(delta_T); // energy in J m-3
				Fthaw = energy / Lheat; // unitless fraction
				if (Fice > Fthaw) {

					// There is some ice left, even after melt
					Fwater += Fthaw;
					Fice -= Fthaw;
					soillayers[lyr].temp = 0.0;
				}
				else {
					// There is more energy available than required for
					// melting the rest of the ice. Therefore the
					// left-over energy will increase the temperature.

					Fwater += Fice;

					// Update soil temperature in this layer (ensures that we don't go back under 0 degrees)
					//soillayers[lyr].temp = max(0.0, (Fthaw - Frac_ice) * Lheat / soillayers[lyr].Ch);
					soillayers[lyr].temp -= Fice * Lheat / soillayers[lyr].Ch;
					soillayers[lyr].temp = max(0.0, soillayers[lyr].temp);

					Fice = 0.0;
				}
			}
			else {
				Fthaw = 0.0;
			}
		} // end of THAWING


		if (Fice < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid ice volume content %g\n", Fice);

		if (Fwater < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid volumetric water content %g\n", Fwater);

		// Water and ice fraction of full volume
		soillayers[lyr].water = Fwater * soillayers[lyr].depth;
		soillayers[lyr].ice = Fice * soillayers[lyr].depth;

		double wtotal_final = soillayers[lyr].water + soillayers[lyr].ice;

		// water balance checks:
		if ((fabs(wtotal_initial - wtotal_final)) > SMALLVOLFRAC)
			fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);

		// Update wcont and wpcont
		if (soillayers[lyr].water > soillayers[lyr].wp) {
			soillayers[lyr].wpcont = 1.0;
			soillayers[lyr].wcont = (soillayers[lyr].water - soillayers[lyr].wp) / soillayers[lyr].awc;
		}
		else {
			soillayers[lyr].wpcont = soillayers[lyr].water / soillayers[lyr].wp;
			soillayers[lyr].wcont = 0.0;
		}

		// Update old temperature after possible phase change
		soillayers[lyr].temp_old = soillayers[lyr].temp;

	} // end of loop through layers
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

	// If still among soil layers then assuming that water in next soil layer is at the top of the layer
	if (lyr < nsoillayers) {
		if (!negligible(soillayers[lyr].water + soillayers[lyr].ice) && lyr < nsoillayers) {
			thaw += soillayers[lyr].water / (soillayers[lyr].water + soillayers[lyr].ice) * soillayers[lyr].depth;
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

	// Minimum dz for top snow layer [m] (should be calculated from Von Neumann Stability Analysis k=gamma deltaT/(deltaX)^2)
	double min_slyr_dz = 0.1;

	for (int ts = 0; ts < date.timesteps; ts++) {

		int nsnowsoillayers = snow.nsnowlayers + nsoillayers;
		int nsnowsoilpadlayers = snow.nsnowlayers + nsoillayers + PAD_LAYERS;
		int tmp_snowlayers = snow.nsnowlayers; //////////////////////////

		if (nsnowsoilpadlayers != tpara.dz.size()) {
			tpara.set_sizes(nsnowsoilpadlayers);
		}

		// Update heat capacity and thermal conductivity for snow, soil and padding layers
		update_heatcap_and_cond();

		// Snow layers (0 == bottom layer, nsnowlayers - 1 == top snow layer)
		int slyr = snow.nsnowlayers - 1;
		for (int lyr = 0; lyr < snow.nsnowlayers; lyr++) {
			if (lyr == 0 && snow.snowlayers[slyr].depth * M_PER_MM < min_slyr_dz) {
				double frac_snow = snow.snowlayers[slyr].depth / 100.0;
				tpara.dz[lyr] = min_slyr_dz;
				tpara.cap[lyr] = snow.snowlayers[slyr].Ch * frac_snow + Cp_air * (1.0 - frac_snow);
				tpara.diff[lyr] = snow.snowlayers[slyr].Kh * frac_snow + Kair * (1.0 - frac_snow);
			}
			else {
				tpara.dz[lyr] = snow.snowlayers[slyr].depth * M_PER_MM;
				tpara.cap[lyr] = snow.snowlayers[slyr].Ch;
				tpara.diff[lyr] = snow.snowlayers[slyr].Kh;
			}
			tpara.psi[lyr] = snow.snowlayers[slyr].temp;

			slyr--;
		}

		// Soil layers
		for (int lyr = snow.nsnowlayers; lyr < nsnowsoillayers; lyr++) {
			tpara.dz[lyr] = soillayers[lyr - snow.nsnowlayers].depth * M_PER_MM;
			tpara.cap[lyr] = soillayers[lyr - snow.nsnowlayers].Ch;
			tpara.diff[lyr] = soillayers[lyr - snow.nsnowlayers].Kh;
			tpara.psi[lyr] = soillayers[lyr - snow.nsnowlayers].temp;
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
		// Upper boundary condition
		double f_top = -2.0 * tpara.diff[0] * (tpara.psi[0] - patch.get_climate().temp) / tpara.z[0] / date.timesteps;

		// Lower boundary condition
		double f_bot = 0.0;

		// Build diagonals and RHS vector
		build_diags_rhs(date.dt / date.timesteps, f_top, f_bot,
			tpara.cap, tpara.diff, tpara.cond,
			tpara.upper, tpara.main, tpara.lower, tpara.rhs,
			tpara.z, tpara.dz, tpara.psi, tpara.sink, nsnowsoilpadlayers);

		// Solve tridiagonal system
		solve_tridiagonal(tpara.upper, tpara.main, tpara.lower, tpara.rhs, tpara.psi1);

		// Update temperatures

		// Set new snow layer temperatures (0 == bottom layer, nsnowlayers - 1 == top snow layer)
		slyr = snow.nsnowlayers - 1;
		for (int lyr = 0; lyr < snow.nsnowlayers; lyr++) {
			snow.snowlayers[slyr].temp = tpara.psi1[lyr];
			slyr--;
		}

		// Set new soil layer temperatures
		for (int lyr = snow.nsnowlayers; lyr < nsnowsoillayers; lyr++) {
			soillayers[lyr - snow.nsnowlayers].temp = tpara.psi1[lyr];
		}

		// Set new padding layer temperatures
		for (int lyr = nsnowsoillayers; lyr < nsnowsoilpadlayers; lyr++) {
			pad[lyr - nsnowsoillayers].temp = tpara.psi1[lyr];
		}

		// Phase change between water and ice and adjust snow layer temperatures accordingly
		snow.phase_change();

		// Phase change between water and ice and adjust soil layer temperatures accordingly
		phase_change();

		// Balance ice and water in snow layers after melt
		snow.ice_water_snow_balance();

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

		//PRINTING OUT DAILY SNOW AND SOIL PROPERTIES
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
		}
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	// Thaw depth calculation
	thaw_depth();

	// Save daily snow depth
	snow.ddepth[date.day] = snow.depth();

	// Calculate mean monthly snow depth
	if (date.islastday)
		snow.mdepth[date.month] = mean(snow.ddepth + date.day - date.ndaymonth[date.month] + 1,
			date.ndaymonth[date.month]);
}

void Soilcolumn::serialize(ArchiveStream& arch) {
	arch & soillayers
		& NH4_input
		& NO3_input
		& anfix
		& anfix_mean
		& anfix_calc
		& wcont_evap
		& rain_melt
		& max_rain_melt
		& pad
		& solvesomcent_beginyr
		& solvesomcent_endyr;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Soil member functions
////////////////////////////////////////////////////////////////////////////////

/*void Soil::init_states() {

	// Initialises certain member variables

	alag = 0.0;
	exp_alag = 1.0;
	cpool_slow = 0.0;
	cpool_fast = 0.0;
	decomp_litter_mean = 0.0;
	k_soilfast_mean = 0.0;
	k_soilslow_mean = 0.0;
	wcont_evap = 0.0;
	snowpack = 0.0;
	orgleachfrac = 0.0;

	// Extra initialisation
	aorgCleach = 0.0;
	aorgNleach = 0.0;
	anfix = 0.0;
	aminleach = 0.0;

	lKorg = lKpeat = lKmin = lKwater = lKice = lKair = 0.0;

	mwcontupper = 0.0;
	mwcontlower = 0.0;

	for (int mth=0; mth<12; mth++) {
		mwcont[mth][0] = 0.0;
		mwcont[mth][1] = 0.0;
		fnuptake_mean[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
	}

	std::fill_n(dwcontupper, Date::MAX_YEAR_LENGTH, 0.0);
	std::fill_n(dwcontlower, Date::MAX_YEAR_LENGTH, 0.0);
	// for fire
	std::fill_n(dthaw, Date::MAX_YEAR_LENGTH, 0.0);


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
	
	NO2_mass = 0.0;
	NO2_mass_w = 0.0;
	NO2_mass_d = 0.0;
	NO_mass = 0.0;
	NO_mass_w = 0.0;
	NO_mass_d = 0.0;
	N2O_mass = 0.0;
	N2O_mass_w = 0.0;
	N2O_mass_d = 0.0;
	N2_mass = 0.0;
	NH4_mass = 0.0;
	NO3_mass = 0.0;
	NH4_input = 0.0;
	NO3_input = 0.0;
	anmin = 0.0;
	animmob = 0.0;
	aminleach = 0.0;
	aorgNleach = 0.0;
	aorgCleach = 0.0;
	anfix = 0.0;
	anfix_calc = 0.0;
	anfix_mean = 0.0;
	snowpack_NH4_mass = 0.0;
	snowpack_NO3_mass = 0.0;
	
	pH = 6.5;
	
	dperc = 0.0;

	solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	solvesomcent_endyr   = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);

	temp25 = 10;

	/////////////////////////////////////////////////////
	// Arctic and wetland initialisation

	// Daily heterotrophic respiration. Only ever nonzero on peatland stands.
	dcflux_soil = 0.0;

	// Indices
	IDX = 0;
	SIDX = 0;
	MIDX = 0;
	SIDX_old = 0;

	firstTempCalc = true;
	firstHydrologyCalc = true;
	SIDX_old = 0;
	IDX = 0;
	nsublayer1 = 0;
	nsublayer2 = 0;
	num_evaplayers = 0;

	// initialise snow variables
	snow_active = false;
	snow_active_layers = 0;
	snow_days = 0;
	snow_days_prev = 365;
	dec_snowdepth = 0.0;

	// Peatland hydrology variables
	awtp = 0.0;
	acro_depth = 0.0;
	cato_depth = 0.0;

	acro_co2 = 934.0; // [mimol L-1]

	wtd = 0.0; // [-100, +300] mm

	acro_por = acrotelm_por - Fgas;
	cato_por = catotelm_por - Fgas;

	Wtot = 0.0;
	stand_water = 0.0;

	dmoss_wtp_limit = 1.0;
	dgraminoid_wtp_limit = 1.0;

	k_O2 = 0.0;
	k_CO2 = 0.0;
	k_CH4 = 0.0;
	Ceq_O2 = 0.0;
	Ceq_CO2 = 0.0;
	Ceq_CH4 = 0.0;

	ch4_store = 0.0;
	co2_store = 0.0;

	maxthawdepththisyear = 0.0;
	thaw = 0.0;
	snowdens = snowdens_start; // kg/m3
	dsnowdepth = 0.0;

	for (int mth = 0; mth<12; mth++) {
		mthaw[mth] = 0.0;
		msnowdepth[mth] = 0.0;
		mwtp[mth] = 0.0;
	}

	for (int sly = 0; sly<NSUBLAYERS_ACRO; sly++)
		sub_water[sly] = 0.0;

	for (int d = 0; d<365; d++) {

		wtp[d] = 0.0;

		for (int ly = 0; ly<NLAYERS; ly++) {

			Frac_ice[ly] = 0.0;
			T_soil[ly] = 0.0;

			if (d == 0) {
				// First day only:
				Frac_water[ly] = 0.0;
				T[ly] = 0.0;
				T_old[ly] = 0.0;
				Frac_org[ly] = 0.0;
				Frac_peat[ly] = 0.0;
				Frac_min[ly] = 0.0;
				Fpwp_ref[ly] = 0.0;
				Frac_water_belowpwp[ly] = 0.0;
				Frac_ice_yesterday[ly] = 0.0;
				Dz[ly] = 0.0;
				por[ly] = 0.0;
				rootfrac[ly] = 0.0;
				CH4_yesterday[ly] = 0.0;
				CO2_soil_yesterday[ly] = 0.0;
				volume_liquid_water[ly] = 0.0;
				total_volume_water[ly] = 0.0;
				tiller_area[ly] = 0.0001; // The max value
				CH4_ebull_ind[ly] = 0.0;
				CH4_ebull_vol[ly] = 0.0;
				CH4_gas_yesterday[ly] = 0.0;
				CH4_diss_yesterday[ly] = 0.0;
				CH4_gas_vol[ly] = 0.0;
			}

			CH4[ly] = 0.0;
			CO2_soil[ly] = 0.0;
			O2[ly] = 0.0;
			CO2_soil_prod[ly] = 0.0;
			CH4_prod[ly] = 0.0;
			CH4_gas[ly] = 0.0;
			CH4_diss[ly] = 0.0;
			Frac_air[ly] = 0.0;
		}
	}

	for (int sl = 0; sl<NSOILLAYER; sl++) {
		wcont[sl] = 0.0;
		alwhc[sl] = 0.0;
		alwhc_init[sl] = 0.0;
	}

	awcont_upper = 0.0;

	nsublayer1 = (int)(SOILDEPTH_UPPER / Dz_soil);		// 5, typically
	nsublayer2 = (int)(SOILDEPTH_LOWER / Dz_soil);		// 10, typically
	num_evaplayers = (int)(SOILDEPTH_EVAP / Dz_soil);	// 2, typically
	whc_evap = soiltype.awc[0] + soiltype.awc[1];		// mm

	// Depths of the acrotelm & catotelm 
	acro_depth = NACROTELM * Dz_acro;
	cato_depth = SOILDEPTH_UPPER + SOILDEPTH_LOWER - acro_depth;

	// end of initialisation
}*/

/// Get soil clay fraction
/*double Soil::get_clayfrac() {

	// Return soil clay fractionh
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.clay_frac_peat;
	else
		return soiltype.clay_frac;
}*/

/// Get soil sand fraction
/*double Soil::get_sandfrac() {

	// Return soil clay fraction
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.sand_frac_peat;
	else
		return soiltype.sand_frac;
}*/

/// Get soil silt fraction
/*double Soil::get_siltfrac() {

	// Return soil silt fraction
	if (patch.stand.is_highlatitude_peatland_stand())
		return soiltype.silt_frac_peat;
	else
		return soiltype.silt_frac;
}*/

// Get the Water Filled Pore Space of the layer. Currently when Century is restricted to the upper of the old 2 layered LPJ soil column, layer defaults to 0.
/*double Soil::wfps(int layer = 0) const {
	double wcont = layer == 0 ? get_soil_water_upper() : get_soil_water_lower();

    return (wcont * soiltype.gawc[layer] + soiltype.gwp[layer]) / soiltype.gwsats[layer];
}*/

/*double Soil::get_soil_temp_25() const {

	// Return soil temperature at 25cm depth

	if (iftwolayersoil)
		return temp_analyticsoln;
	else
		return temp25;
}*/


/// Performs daily updates to soil temperatures using an analytic function
/**
 * Call each simulation day (if iftwolayersoil true) following update of daily air temperature 
 * prior to canopy exchange and SOM dynamics
 * \param climate   The climate to use to update soil water
 * \param depth     The depth at which the soil temperature is calculated
 */
/*void Soil::soil_temp_analytic(const Climate& climate, double depth) {

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
	double day[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

	if ((date.year == 0 || date.year == patch.stand.first_year) && date.month == 0 && !date.islastday) {

		// First month of simulation, use air temperature for soil temperature
		temp_analyticsoln = climate.temp;
	}
	else {

		if (date.islastday) {

			// Linearly interpolate soil thermal diffusivity given mean
			// soil water content

			if (mwcontupper < 0.15)
				k = ((soiltype.thermdiff_15 - soiltype.thermdiff_0) / 0.15 *
					mwcontupper + soiltype.thermdiff_0) * DIFFUS_CONV;
			else
				k = ((soiltype.thermdiff_100 - soiltype.thermdiff_15) / 0.85 *
					(mwcontupper - 0.15) + soiltype.thermdiff_15) * DIFFUS_CONV;

			// Calculate parameters alag and exp(-alag) from Eqn 2

			alag = depth / sqrt(k / HALF_OMEGA); // from Eqn 1
			exp_alag = exp(-alag);

		}

		// Every day, calculate linear model for trend in daily air
		// temperatures for the last 31 days: temp_day = a + b * day

		double buffer[31];
		climate.dtemp_31.to_array(buffer);
		regress(day, buffer, 31, a, b);

		// Calculate soil temperature

		temp_lag = a + b * (30.0 - alag * LAG_CONV);

		// Calculate temperatures the old way
		temp_analyticsoln = climate.atemp_mean + exp_alag*(temp_lag - climate.atemp_mean);
	}
}*/


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
/*void Soil::hydrology_lpjf_twolayer(const Climate& climate, double fevap) {

	// DESCRIPTION
	// This version replicates the hydrology_lpjf function in LPJ-GUESS v4.0, so 
	// it uses 2 soil layers, and does not consider impedance or ice formation. 

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer, 14=bottom layer] as
	//              fraction of available water holding capacity (AWC)
	//				This method ensures that layers 2-4 have identical wcont values, as do layers 5-14
	// wcont_evap = water content of evaporation sublayers (20 cm) at top of soil column
	//              as fraction of available water holding capacity (AWC)
	//				This method ensures that layers 0-1 have identical wcont_evap values
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
	// awc        = array containing available water holding capacity of soil
	//              layers (mm rainfall) [0=upper layer]
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface
	// snowpack   = depth of snow (mm)

	// INPUT AND OUTPUT PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayer at top of upper soil layer
	//              as fraction of available water holding capacity (AWC)
	// awcont_upper = wcont averaged over the growing season - guess2008
	// dperc      = daily percolation beyond system (mm)

	// OUTPUT PARAMETER
	// runoff     = total daily runoff from all soil layers (mm/day)

	const double SOILDEPTH_EVAP = 200.0;
	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double BASEFLOW_FRAC = 0.5;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double K_DEPTH = 0.4;
	const double K_AET = 0.52;
	// Fraction of total (vegetation) AET from upper soil layer that is derived
	// from the top K_DEPTH (fraction) of the upper soil layer
	// (parameters for calculating K_AET_DEPTH below)
	const double K_AET_DEPTH = (SOILDEPTH_UPPER / SOILDEPTH_EVAP - 1.0) *
		(K_AET / K_DEPTH - 1.0) / (1.0 / K_DEPTH - 1.0) + 1.0;
	// Weighting coefficient for AET flux from evaporation layer, assuming active
	//   root density decreases with soil depth
	// Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
	//   K_DEPTH=0.4, K_AET=0.52

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	double aet;								// AET for a particular layer and individual (mm)
	double aet_twolayer[2] = {0.0, 0.0};	// total AET for each Gerten-type soil layer (mm)
	double perc_frac;
	double aet_total = 0.0;

	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int s = 0; s<NSOILLAYER; s++) {
			aet = patch.pft[indiv.pft.id].fwuptake[s] * indiv.aet;
			aet_total += aet;
			if (s < NSOILLAYER_UPPER)
				aet_twolayer[0] += aet;
			else
				aet_twolayer[1] += aet;
		}
		vegetation.nextobj();
	}

	// Get the two-layer versions of the wcont values
	double wcont_twolayer[2] = { 0.0, 0.0 };	// wcont for each Gerten-type soil layer
	wcont_twolayer[0] = get_soil_water_upper();
	wcont_twolayer[1] = get_soil_water_lower();
	double twolayerwcont_evap = get_layer_soil_water_evap();

	// Evaporation from soil surface

	// guess2008 - changed to wcont_evap**2, as in LPJ-mL
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	double evap = 0.0;
	if (snowpack < 10.0) {					// evap only if snow depth < 10mm
		evap = climate.eet * PRIESTLEY_TAYLOR * twolayerwcont_evap * twolayerwcont_evap * fevap;
	}

	// Implement in- and outgoing fluxes to upper soil layer
	// BLARP: water content can become negative, though apparently only very slightly
	//    - quick fix implemented here, should be done better later

	wcont_twolayer[0] += (rain_melt - aet_twolayer[0] - evap) / soiltype.gawc[0];
	if (wcont_twolayer[0] != 0.0 && wcont_twolayer[0] < 0.0001) { // guess2008 - bugfix
		wcont_twolayer[0] = 0.0;
	}

	// Surface runoff
	double runoff_surf = 0.0;
	if (wcont_twolayer[0] > 1.0) {
		runoff_surf = (wcont_twolayer[0] - 1.0) * soiltype.gawc[0];
		wcont_twolayer[0] = 1.0;
	}

	// Update water content in evaporation layer for tomorrow

	twolayerwcont_evap += (rain_melt - aet_twolayer[0] * SOILDEPTH_EVAP*K_AET_DEPTH / SOILDEPTH_UPPER - evap)
		/ soiltype.gawc[0];

	if (twolayerwcont_evap < 0) {
		twolayerwcont_evap = 0;
	}

	if (twolayerwcont_evap > wcont_twolayer[0]) {
		twolayerwcont_evap = wcont_twolayer[0];
	}

	// Percolation from evaporation layer
	double perc = 0.0;
	if (percolate) {
		perc = min(SOILDEPTH_EVAP / SOILDEPTH_UPPER*soiltype.perc_base*pow(twolayerwcont_evap, soiltype.perc_exp),
			max_rain_melt);
	}
	twolayerwcont_evap -= perc / soiltype.gawc[0];

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	double runoff_drain = 0.0;

	for (int s = 1; s<NSOILLAYER_SIMPLE; s++) {

		// Percolation
		// Allow only on days with rain or snowmelt (Dieter Gerten, 021216)

		if (percolate) {
			perc = min(soiltype.perc_base*pow(wcont_twolayer[s - 1], soiltype.perc_exp), max_rain_melt);
		}
		else {
			perc = 0.0;
		}
		perc_frac = min(perc / soiltype.gawc[s - 1], wcont_twolayer[s - 1]);

		wcont_twolayer[s - 1] -= perc_frac;
		wcont_twolayer[s] += perc_frac * soiltype.gawc[s - 1] / soiltype.gawc[s];
		if (wcont_twolayer[s] > 1.0) {
			runoff_drain += (wcont_twolayer[s] - 1.0)*soiltype.gawc[s];
			wcont_twolayer[s] = 1.0;
		}

		// Deduct AET from this soil layer
		// BLARP! Quick fix here to prevent negative soil water

		wcont_twolayer[s] -= aet_twolayer[s] / soiltype.gawc[s];
		if (wcont_twolayer[s] < 0.0) {
			wcont_twolayer[s] = 0.0;
		}
	}

	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_baseflow = 0.0;
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC*soiltype.perc_base*pow(wcont_twolayer[NSOILLAYER_SIMPLE - 1], soiltype.perc_exp);
		// guess2008 - Added "&& rain_melt >= runoff_surf" to guarantee nonnegative baseflow.
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer

		perc_frac = min(perc_baseflow / soiltype.gawc[NSOILLAYER_SIMPLE - 1], wcont_twolayer[NSOILLAYER_SIMPLE - 1]);
		wcont_twolayer[NSOILLAYER_SIMPLE - 1] -= perc_frac;
		runoff_baseflow = perc_frac * soiltype.gawc[NSOILLAYER_SIMPLE - 1];
	}

	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;


	// update awcont_upper
	// Original algorithm by Thomas Hickler

	// Reset awcont_upper on the first day of every year
	if (date.day == 0) {
		awcont_upper = 0.0;
		patch.growingseasondays = 0;
	}
	
	// If it's warm enough for growth, update awcont_upper with this day's wcont
	if (climate.temp > 5.0) {
		awcont_upper += wcont_twolayer[0];
		patch.growingseasondays++;
	}

	// Average awcont_upper on the last day of every year
	if (date.islastday && date.islastmonth) {
		if (patch.growingseasondays > 1)
			awcont_upper /= (double)patch.growingseasondays;
		else
			awcont_upper = 0.0; // This will lead to no establishment for this PFT/species if ifdroughtlimitedestab is TRUE
	}


	// Use the two-layer values calculated above to update the soil properties needed in all NSOILLAYER layers

	// First update wcont for the evaporation layers
	for (int s = 0; s<2; s++) {
		set_layer_soil_water(s, twolayerwcont_evap);
	}

	// Now update wcont for the remaining NSOILLAYER_UPPER layers
	// NB: we want to set the wcont[] value in layers 2-4 to be consistent with twolayerwcont_evap and
	//	wcont_twolayer[0] calculated above. Let wcont_new be this value. Thus:
	//	(2 * twolayerwcont_evap + 3 * wcont_new) / 5 = wcont_twolayer[0], 
	//	giving:
	//	wcont_new = (5 * wcont_twolayer[0] - 2 * twolayerwcont_evap) / 3
	double wcont_new = (5.0 * wcont_twolayer[0] - 2.0 * twolayerwcont_evap) / 3.0;
	wcont_new = min(max(wcont_new, 0.0), 1.0); // Limit to valid values before updating wcont[]

	for (int s = 2; s<NSOILLAYER_UPPER; s++) {
		set_layer_soil_water(s, wcont_new);
	}

	// Next update wcont for the deepest NSOILLAYER_LOWER-NSOILLAYER_UPPER layers
	for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {
		set_layer_soil_water(s, wcont_twolayer[1]);
	}

	// Finally, update wcont_evap, Frac_water etc. in this patch.soil object
	update_soil_water();

}*/


/*bool Soil::ice_in_top_layer() {

	const double max_ice_fraction = 0.05; // 5% avoids irrigation restrictions when there are tiny amounts of ice in the soil
	for (int ly = 0; ly < NSOILLAYER_UPPER; ly++) {
		if (Frac_ice[ly + IDX] > max_ice_fraction) return true;
	} 

	return false;
}*/


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
/*void Soil::hydrology_lpjf(const Climate& climate, double fevap) {

	// DESCRIPTION
	// Soil class method used to update the water content for each soil layer 
	// given snow melt, rainfall, evapotranspiration from vegetation (AET) and 
	// percolation between layers. Runoff is also calculated.
	// Updated to deal with 15 soil layers and ice formation.  

	// INPUT PARAMETERS
	// fevap      = fraction of modelled area (grid cell or patch) subject to
	//              evaporation from soil surface

	// UPDATED PARAMETERS
	// wcont      = array containing water content of soil layers [0=upper layer, 14=bottom layer] as
	//              fraction of available water holding capacity (AWC)
	// wcont_evap = water content of evaporation sublayers (20 cm) at top of soil column
	//              as fraction of available water holding capacity (AWC)
	// awcont     = wcont averaged over the growing season
	// dperc      = daily percolation beyond system (mm)
	// runoff     = total daily runoff from all soil layers (mm/day)

	// depth of sublayer at top of upper soil layer, from which evaporation is
	// possible (NB: must not exceed value of global constant SOILDEPTH_UPPER)
	const double SOILDEPTH_EVAP = 200.0;
	// Fraction of standard percolation amount from lower soil layer that is
	// diverted to baseflow runoff
	const double BASEFLOW_FRAC = 0.5;

	const double MIN_WATER_PERC = 0.0000001;	// mm - Minimum water amount for percolation to occur
	const double maxerr = 0.0001;			// mm - max error allowed

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}


	// *** EVAP - based on initial wcont_evap, as in v4.0 *** 

	double evap_init = 0;
	double Faw_evap_init = 0.0;		// mm
	double awc_init = 0.0;

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_evap_init += wcont[s] * soiltype.awc[s]; // mm;
		awc_init += soiltype.awc[s];
	}

	double wcont_evap_init = Faw_evap_init / awc_init;

	if (snowpack < 10.0) {	// evap only if snow depth < 10mm
							// The potential evaporation:
		evap_init = climate.eet * PRIESTLEY_TAYLOR * wcont_evap_init * wcont_evap_init * fevap;
		// Below, we will limit evaporation to a value <= the available water in the evaporation layer
	}


	// *** AET ***

	// Reset
	double aet = 0.0; // AET for a particular layer and individual (mm)
	double aet_total = 0.0;
	double aet_evap_layers = 0.0;
	double aet_top_layer = 0.0;
	// total AET for each soil layer (mm)
	double aet_layer[NSOILLAYER];
	for (int ly = 0; ly < NSOILLAYER; ly++) {
		aet_layer[ly] = 0.0; // Initialise
	}


	// Sum AET for across all vegetation individuals
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		for (int ly = 0; ly<NSOILLAYER; ly++) {
			aet = patch.pft[indiv.pft.id].fwuptake[ly] * indiv.aet;
			aet_layer[ly] += aet;
			aet_total += aet;

			// AET from evaporation layers
			if (ly < num_evaplayers)
				aet_evap_layers += aet;

			// AET from top (50cm) layer
			if (ly < NSOILLAYER_UPPER)
				aet_top_layer += aet;
		}
		vegetation.nextobj();
	}


	// *** INITIALISE ***

	// First determine how much water and ice is currently held in each soil layer

	// available water for each soil layer (mm)
	double Faw_layer[NSOILLAYER];
	// ice each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// water that can still be added to each soil layer (mm)
	double potential_layer[NSOILLAYER];
	// Prepare for water balance tests
	double initial_water_in_column = 0.0;


	// Top and bottom layers, initial Faw and potential
	double Faw_layer1 = 0.0;
	double Faw_layer2 = 0.0;
	double potential_layer1 = 0.0;
	double potential_layer2 = 0.0;


	for (int ly = 0; ly < NSOILLAYER; ly++) {

		if (wcont[ly] < 0.0 || wcont[ly] > 1.0) {
			dprintf("Soil::hydrology_lpjf - bad wcont!\n");
			return;
		}

		Faw_layer[ly] = wcont[ly] * soiltype.awc[ly]; // mm

		ice_layer[ly] = Frac_ice[ly + IDX] * Dz[ly + IDX]; // mm

		// Inital water in this layer, before AET
		double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

		// for water balance checks
		initial_water_in_column += layerwater;

		// Shouldn't happen
		if (aet_layer[ly] > Faw_layer[ly]) {
			int test = 1;
		}

		// AW
		// Now remove AET as it's already accounted for
		Faw_layer[ly] -= min(aet_layer[ly], Faw_layer[ly]);
		// Update layer water after subtracting AET
		layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

		// Update wcont
		wcont[ly] = Faw_layer[ly] / soiltype.awc[ly]; // mm

		if (ly < NSOILLAYER_UPPER)
			Faw_layer1 += Faw_layer[ly]; // mm
		else
			Faw_layer2 += Faw_layer[ly]; // mm


		// POTENTIAL
		potential_layer[ly] = aw_max[ly] - layerwater;

		if (ly < NSOILLAYER_UPPER)
			potential_layer1 += potential_layer[ly]; // mm
		else
			potential_layer2 += potential_layer[ly]; // mm


		// Check balance
		if (potential_layer[ly] < -1.0 * maxerr) {
			// fail("Soil::hydrology_lpjf - error in a soil layer's water balance!\n");	
			dprintf("ERROR: Soil::hydrology_lpjf - error in a soil layer's water balance! potential_layer#: %g\n", ly);
		}
		else {
			potential_layer[ly] = max(0.0, potential_layer[ly]); // catches small negative values
		}
	} // for loop (ly) 


	// *** EVAPORATION and TOP LAYER ***

	// Update wcont_evap
	// wcont_evap is the available, nonfrozen water down to SOILDEPTH_EVAP, divided by the capacity of those layers to hold water 
	wcont_evap = 0.0;
	double awc_count = 0.0;
	double Faw_layer_evap = 0.0;		// mm
	double potential_evap_layer = 0.0;	// mm
	double potential_top_layer = 0.0;	// mm

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_layer_evap += Faw_layer[s];
		awc_count += soiltype.awc[s];
		potential_evap_layer += potential_layer[s];
	}

	wcont_evap = Faw_layer_evap / awc_count;


	// *** EVAPORATION FROM BARE SURFACE (after AET removed) ***
	double evap = min(evap_init, max(0.0, Faw_layer_evap));

	// REMOVAL of water from the evaporation layers of the soil (evap > 0).
	if (evap > 0.0 && Faw_layer_evap > 0.0) {

		// Specifying evap > 0 ensures that there is at least some water in the evap layers before trying to remove water

		for (int s = 0; s<num_evaplayers; s++) {

			// Remove the water in proportional to the available water remaining in each layer
			double proportion = evap * (Faw_layer[s] / Faw_layer_evap); // mm
			Faw_layer[s] -= proportion;

			potential_layer[s] += proportion; // potential INCREASES
			wcont[s] = Faw_layer[s] / soiltype.awc[s];
			oob_check_wcont(wcont[s]);
		}
	} // > 0.0?
	else {
		evap = 0.0; // No evaporation in this case
	}


	// Implement in- and outgoing fluxes to the top 50cm soil layers, taking into account AET in those layers and evaporation
	// ice_impedance[0] reduces the input when there is ice present in the top 20cm.
	double water_flux_in = rain_melt;	// mm

	// Initialise runoff variables

	// Surface runoff (mm) 
	double runoff_surf = 0.0;
	// Baseflow runoff (Dieter Gerten 021216) (rain or snowmelt days only)
	double runoff_drain = 0.0;
	double runoff_baseflow = 0.0;


	// Potential in the top layer
	potential_top_layer = 0.0;
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		potential_top_layer += potential_layer[s];
	}

	// *** INPUT TO TOP LAYER ***

	if (water_flux_in > 0) {

		// ADDITION of water to the top layers (rain_melt > 0)
		if (water_flux_in < potential_top_layer) {

			// The upper soil layers can absorb all of this water today
			// so add the water to layers in proportion to their capacity
			for (int s = 0; s<NSOILLAYER_UPPER; s++) {

				double water_input_ly = 0.0;
				if (potential_top_layer > 0.0)
					water_input_ly = water_flux_in * (potential_layer[s] / potential_top_layer);

				Faw_layer[s] += water_input_ly;
				potential_layer[s] -= water_input_ly;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
			}

			// No surface runoff
			runoff_surf = 0.0;
		}
		else {

			// The upper soil layers cannot absorb all of this water today, so take up what can be absorbed 
			// (i.e. potential_top_layer) and add the rest to surface runoff

			for (int s = 0; s<NSOILLAYER_UPPER; s++) {
				Faw_layer[s] += potential_layer[s];
				potential_layer[s] = 0.0;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
			}

			// Update surface runoff
			runoff_surf = water_flux_in - potential_top_layer;
		}

		for (int s = 0; s<NSOILLAYER_UPPER; s++) {
			oob_check_wcont(wcont[s]);
		}
	}


	// Update wcont_evap again, after the input of rain_melt
	wcont_evap = 0.0;
	Faw_layer_evap = 0.0;	// mm

	for (int s = 0; s < num_evaplayers; s++) {
		Faw_layer_evap += Faw_layer[s];
	}

	wcont_evap = Faw_layer_evap / awc_count;
	oob_check_wcont(wcont_evap);


	// Update Faw_layer1 and wcont1

	double awc1 = 0.0;
	Faw_layer1 = 0.0;
	for (int s = 0; s < NSOILLAYER_UPPER; s++) {
		Faw_layer1 += Faw_layer[s];
		awc1 += soiltype.awc[s];
	}

	double wcont_layer1 = Faw_layer1 / awc1;


	// *** PERCOLATION ***
	// We assume there are two percolation layers. 
	// 1: the top 50cm soil layers considered in the Gerten et al. scheme
	// 2: Lower 1.0 m layer considered in the Gerten et al. scheme

	if (percolate) {

		// Percolation from top layer if there is enough water (at least MIN_WATER_PERC)
		double perc = 0.0;

		if (Faw_layer1 > MIN_WATER_PERC) {

			// Percolation from the top  layer (limited to available liquid water)
			perc = min(soiltype.perc_base * pow(wcont_layer1, soiltype.perc_exp), Faw_layer1);

			for (int s = 0; s<NSOILLAYER_UPPER; s++) {

				// Take water proportionally from the layers.
				double perc_layer = perc * (Faw_layer[s] / Faw_layer1);
				Faw_layer[s] -= perc_layer;
				potential_layer[s] += perc_layer;

				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
			}
		}


		// This is the percolation (mm) from the top layers
		double perc_from_above = perc;

		// Percolation and fluxes to and from the soil layers below the top 50cm

		// Excess water transferred to runoff
		// Eqns 26, 27, 31, Haxeltine & Prentice 1996

		double water_input = perc_from_above;

		// First we need to add the water percolated from above (perc_from_above) to this layer
		if (perc_from_above > potential_layer2) {
			water_input = potential_layer2; // Fill the layer and reduce water_input
			// Update runoff_drain
			runoff_drain = perc_from_above - potential_layer2; // The excess is run off and not percolated further
		}

		perc_from_above = 0.0;

		double Faw_perc_layer_2 = 0.0;
		double awc_perc_layer_2 = 0.0;
		double wcont_perc_layer_2 = 0.0;

		if (water_input >= 0) {

			// Now add any water (water_input) to layers 5-14 in proportion to the capacity, and recalculate Faw, wcont and awc
			for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {

				double water_input_ly = 0.0;
				if (potential_layer2 > 0.0)
					water_input_ly = water_input * (potential_layer[s] / potential_layer2);

				Faw_layer[s] += water_input_ly;
				potential_layer[s] -= water_input_ly;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);

 				Faw_perc_layer_2 += Faw_layer[s];
				awc_perc_layer_2 += soiltype.awc[s];
			}

			// Update wcont for these layers before percolation
			wcont_perc_layer_2 = Faw_perc_layer_2 / awc_perc_layer_2;
		}

		// Percolation 3 -> base

		double perc_from_base = 0.0;
		if (percolate && Faw_perc_layer_2 > MIN_WATER_PERC) {

			// ice_impedance[2] reduces the base flow rate when there is ice present in the 50-150 cm layers.
			perc_from_base = min(BASEFLOW_FRAC * soiltype.perc_base * pow(wcont_perc_layer_2, soiltype.perc_exp), max_rain_melt);
			perc_from_base = min(perc_from_base, Faw_perc_layer_2); // Only available water in these layers can percolate

			// As in LPJ-GUESS v4.0
			if (perc_from_base > rain_melt - runoff_surf && rain_melt >= runoff_surf)
				perc_from_base = rain_melt - runoff_surf;
		}
		else {
			// No percolation
			perc_from_base = 0.0;
		}

		// Now remove the water (i.e. perc_from_base) from layers 5-14 in proportion to the available water
		if (perc_from_base > 0) {

			for (int s = NSOILLAYER_UPPER; s<NSOILLAYER; s++) {

				// Remove the water (perc_from_base) in proportional to the available water remaining in each layer
				double proportion = perc_from_base * (Faw_layer[s] / Faw_perc_layer_2); // mm

				Faw_layer[s] -= proportion;

				potential_layer[s] += proportion;
				wcont[s] = Faw_layer[s] / soiltype.awc[s];
				oob_check_wcont(wcont[s]);
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


	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// water added when patch.stand.is_true_wetland_stand() should be be subtracted from runoff
	// in proportion to its components
	if (patch.stand.is_true_wetland_stand() && ifsaturatewetlands) {

		if (runoff <= patch.wetland_water_added_today && runoff > 0.0) {

			// Not enough runoff to balance the water added to the wetland, i.e. we added more than was run off.

			// Update wetland_water_added by taking back the amount used for runoff.
			patch.wetland_water_added_today -= runoff;

			// All runoff removed 
			runoff_surf = 0.0;
			runoff_baseflow = 0.0;
			runoff_drain = 0.0;
			runoff = 0.0;

		}
		else if (runoff > patch.wetland_water_added_today  && runoff > 0.0) {

			// Enough runoff to balance the water added to the wetland so we take it back.
			runoff_surf -= patch.wetland_water_added_today * runoff_surf / runoff;
			runoff_drain -= patch.wetland_water_added_today * runoff_drain / runoff;
			runoff_baseflow -= patch.wetland_water_added_today * runoff_baseflow / runoff;
			// Recalculate total runoff 
			runoff = runoff_surf + runoff_drain + runoff_baseflow;

			// Reset wetland_water_added
			patch.wetland_water_added_today = 0.0;

		}
	}

	// save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.awetland_water_added += patch.wetland_water_added_today;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;

	// *** WATER IN BALANCE? ***
	if (DEBUG_SOIL_WATER) {

		double final_water_in_column = 0.0;

		for (int ly = 0; ly < NSOILLAYER; ly++) {

			Faw_layer[ly] = wcont[ly] * soiltype.awc[ly]; // mm
			ice_layer[ly] = Frac_ice[ly + IDX] * Dz[ly + IDX]; // mm

			// Water in this layer
			double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

			// for water balance checks
			final_water_in_column += layerwater;

		} // for loop (ly)

		  // is water in + initial storage = water out + final storage?
		double water_in_storage_in = initial_water_in_column + rain_melt;
		double water_out_storage_out = final_water_in_column + evap + aet_total + runoff;

		if (fabs(water_in_storage_in - water_out_storage_out) > maxerr) {
			dprintf("Soil::hydrology_lpjf - error in the TOTAL water balance!\n");
			return;
		}
	}

	// Drought limited establishment - update awcont_upper
	// Original algorithm by Thomas Hickler
	for (int s = 0; s<NSOILLAYER_UPPER; s++) {

		// Reset on the first day of every year
		if (date.day == 0 && s == 0) {
			awcont_upper = 0.0;
			patch.growingseasondays = 0;
		}

		// If it's warm enough for growth, update awcont with this day's wcont
		if (climate.temp > 5.0) {
			awcont_upper += wcont[s] / (double)NSOILLAYER_UPPER;
			if (s == 0) {
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

	// Finally, recalculate Frac_water based on these calculations and updates to wcont[] 
	update_soil_water();
}*/


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
/*void Soil::hydrology_peat(const Climate& climate, double fevap) {

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


	// PARAMETERS and variables for the calculation of the WTP
    // See: Granberg et al. (1999), Eqns. 1-5.
	
	/// depth of each sphagnum layer [mm]
	const double Dz_sph = 10.0;
	/// surface water content when WTD is at 20 mm
	const double surfw = 0.34;
	/// minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;              
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// gradient limit of upper gradient [mm]
	const double grad_lim = 50.0;
	/// Granberg et al. (1999) parameter, eqn 2.
	double az; 
	/// Depth of sublayers in the acrotelm
	double Dz_sub; 	

	/// Day of year
	int daynum = date.day;

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_water_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	// Update Wtot, whc[], awhc[], and wcont[]
	if (!update_layer_water_content(daynum)) 
		fail();

	// Porosity updates in the presence of ice? 
	bool subtractIceFromWtot = false;
	double acro_por_icy = acro_por; // the pore space minus the ice fraction
	double avgIceFrac = 0.0;
	
	if (subtractIceFromWtot) {	
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) 
			avgIceFrac += (Frac_ice[ly] + Fpwp_ref[ly] - Frac_water_belowpwp[ly])/(double)NACROTELM;
	}

	acro_por_icy -= avgIceFrac;


	// *** TRANSPIRATION *** 
	
	// Retrieve Vegetation for this patch
	Vegetation& vegetation=patch.vegetation;

	// Transpiration from the acrotelm  
	double aet_acrotelm = 0.0;
	double aet = 0.0; // AET for a particular layer and individual (mm)

	// Sum AET for across all vegetation individuals in this patch

	// Calculate fpc_moss
	double fpc_moss_total = 0.0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		for (int ly = IDX; ly < IDX + NACROTELM; ly++)
			pft_acro_root_frac += patch.pft[indiv.pft.id].fwuptake[ly-IDX];

		// this includes mosses, which have root_frac = 1 in the acrotelm
		// Fraction of the aet from the acrotelm only.
		aet_acrotelm += indiv.aet * pft_acro_root_frac;
		
		// Total aet, assumed from acrotelm and catotelm combined 
		aet_total += indiv.aet;

		if (indiv.pft.ismoss())
			fpc_moss_total += indiv.fpc;

		vegetation.nextobj();
	}


	// *** EVAPORATION *** 
	// Evaporation from soil surface

	// Wania et al (2009)
	if (snowpack < 10.0 && wtd < SOILDEPTH_EVAP && !negligible(fevap)) { // i.e. evap only if snow depth < 1cm
		evap = fevap*climate.eet*PRIESTLEY_TAYLOR*(0.99 / (1+exp(-1.0*(-wtd + 98.7)/22.6) + 0.02));
	} 
	else {
		evap = 0.0;
	}
	
	// *** EVAPOTRANSPIRATION *** 

	bool empirical_evapotranspiration = false;
	if (empirical_evapotranspiration) {

		// An empirical evapotranspiration parameterisation from Kim & Verma

		
		// From LPJ-WHyMe:
		Calculate transpiration and evaporation from wetlands together.
		dpet is the equilibrium ET but is the same as the potential ET in
		Kim and Verma (1996).  (BAD choice of variable names in original
		LPJ!)  Hence, we can use dpet as approximation for potential
		evapotranspiration and scale it by the water table depth.  daept
		is allowed to happen even if the wtp is at -300mm because plants
		will be able to extract water from the unsaturated acrotelm or the
		saturated catotelm.  This however makes it impossible to close the
		hydrological cycle in LPJ-WHy. 
		Kim and Verma (1996): ET = PET * (1.02 + 0.00075 * WTP)
		

		// Note 1. that this is now the TOTAL evapotranspiration, replacing aet_total

		// Note 2: An alternative parameterisation was given by Yurova et al. (2006)
		// evapotranspiration=eet*exp[0.02*wtp_cm], where // wtp_cm [100,-30]
		// This would give almost identical results to the Kim & Verma formula below 
		// when the water table is at the surface

		if (snowpack < 10.0) { // i.e. evap only if snow depth < 1cm
			evapotranspiration = climate.eet*(1.02 - 0.00075 * wtd)*fevap; // wtd [-100,+300]
		}
		else {
			evapotranspiration = 0.0;
		}

	} 
	else {

		// This is JUST evapotranspiration from the acrotelm
		evapotranspiration = aet_acrotelm + evap;
	}

	// *** ACROTELM WATER VOLUME ***
	
	az = (acro_por - minvtot) / zmin; // Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22
	
	// Drainage - could possibly be read in for site-specific studies
	runoff_drain = 0; 

	// Update available water 
	Wtot += rain_melt-evapotranspiration-runoff_drain;

	// *** RUNOFF AND RUNON***

	double acrowater = 0.0; // Frac_water[MIDX] * Dz[MIDX];

	for (int ly = IDX; ly < IDX + NACROTELM; ly++)
		acrowater += (Frac_water[ly] + Frac_water_belowpwp[ly]) * Dz[ly];

	double ideal_runoff = exp(-0.01 * wtd);
	if (acrowater > 0.0 && Frac_ice[IDX] < 0.7) // Granberg et al. (1999) use the same ice condition - see f_icestop in their Table 1
		runoff_surf = min(ideal_runoff,acrowater);
	else
		runoff_surf = 0.0;

	// .ins file option: 
	// run on or off - but only when there's flowing water

	double runon = 0.0;

	// runon set in global.ins
	soiltype.runon = wetland_runon;

	if (runoff_surf > 0.0) {

		double runofforon = soiltype.runon;

		if (runofforon < 0) 
			runoff_surf -= runofforon;	// e.g. 3 mm/day for bog-like/hummock conditions

		if (runofforon > 0)
			runon += runofforon;		// 3-5 mm/day for fen-like conditions
	}

	double Wtot_init = Wtot;

	Wtot += runon - runoff_surf;

	// Add standing water, if any
	Wtot += stand_water;

	// Not in LPJ-WHy
	if (Wtot < 0.0) {
		Wtot = 0.0;
		runoff_surf = Wtot_init;
	}

	// Save evapotranspiration data
	if (empirical_evapotranspiration) {

		patch.aaet+=(1.0-fevap)*evapotranspiration;
		patch.aevap+=fevap*evapotranspiration;
		patch.maet[date.month]+=(1.0-fevap)*evapotranspiration;
		patch.mevap[date.month]+=fevap*evapotranspiration;
	} 
	else {
	
		patch.aevap+=evap;
		patch.mevap[date.month]+=evap;
		patch.aaet += aet_total;
		patch.maet[date.month] += aet_total;
	}


	// *** UPDATE WATER TABLE DEPTH ***

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
	
	// NB 
	// e.g. Wtot = 270mm =>
	// Frac_water = 270/300 = 0.9 in the three acrotelm layers
	
	// Wtot <= 140.03611111 =>
	// Frac_water = 0.27419, 0.4186, 0.7075 in the top three layers (if no ice)

	if (Wtot > acro_depth * acro_por) {
		
		// Saturated acrotelm

		wtp[daynum] = Wtot - acro_depth * acro_por; 
		// i.e. wtp[daynum] > 0 above the surface
		
		if (wtp[daynum] > maxh) {
			runoff_surf += wtp[daynum] - maxh; 
			// Update runoff from the surface
			wtp[daynum] = maxh;
		}

		wtd = -wtp[daynum]; // i.e. wtd defined as < 0 above the surface
		Wtot = acro_depth * acro_por; 
		// I.e. Wtot does NOT include standing water
		stand_water = 0.0; // = wtp[daynum]; // we assume no standing water 
	} 
	else {
		
		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (acro_por*acro_depth - Wtot)/(2.0 * az)); 
		
		if (wtd > zmin)
			wtd = 3.0 * (acro_por * acro_depth - Wtot) / (2.0 * (acro_por - minvtot));

		if (wtd > acro_depth) 
			wtd = acro_depth; 

		wtp[daynum] = -wtd;
		stand_water = 0.0;
	}


	// *** RUNOFF ***

	// Save today's runoff in Soil::runoff
	runoff = runoff_surf;
	patch.arunoff+=runoff;
	patch.mrunoff[date.month]+=runoff;
	patch.asurfrunoff += runoff;


	// Depth of acrotelm sublayers - usually 10 mm
	Dz_sub = NACROTELM * Dz_acro / NSUBLAYERS_ACRO;


	// *** WATER CONTENT IN EACH LAYER - THE PROFILE ***

	double surfw_sat = 0.0;
	double value[NSUBLAYERS_ACRO];
	double acro_aw = 0.0; // mm of plant available water in acrotelm
	double cato_aw = 0.0; // mm of plant available water in catotelm

	if (wtd > 0.0) {

		// Below the surface

		// Determine the subsurface layer (0 to NSUBLAYERS_ACRO) within which wtd lies. 
		// Possible values:
		// = 0 if wtd < 10 mm 
		// = NSUBLAYERS_ACRO iff wtd = Dz_acro
		int wtd_layer = (int)(wtd / Dz_sub); // Takes the integer part only

		surfw_sat = max(minvtot, acro_por - az * wtd); // Granberg - eqn 3.

		if (wtd_layer == 0) {
			
			// Total saturation throughout the acrotelm (from 10mm down)
			for (int j = 0; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por; 
			}
		} 
		else {
			
			// Unsaturated zone in the acrotelm - use Granberg's eqn 1.
			for (int i = 0; i < wtd_layer; i++) {
				double z = (i+1)*Dz_sub; // mm
				value[i] = min(acro_por, surfw_sat + (acro_por_icy - surfw_sat) * pow(z/wtd,2));
			}

			// Saturated zone in the acrotelm
			for (int j = wtd_layer; j < NSUBLAYERS_ACRO; j++) {
				value[j] = acro_por; 
			}
		}

		// Now calculate the volumetric water content throughout the acrotelm
		sub_water[0] = 0.5 * (surfw_sat + value[0]);
		for (int j = 1; j < NSUBLAYERS_ACRO; j++) {
			sub_water[j] = 0.5 * (value[j] + value[j-1]); 
		}

		// Water content in each acrotelm layer
		double Ftotal[NLAYERS]; // Volumetric fraction in the acrotelm layers
		for (int ly = IDX; ly < IDX + NACROTELM; ly++)
			Ftotal[ly] = 0.0;

		int sublayers_per_layer = (int)Dz_acro/(int)Dz_sub; 

		int count = 0;
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
		
			double layerwater = 0.0;
			for (int subly = count*sublayers_per_layer; subly < (count+1)*sublayers_per_layer; subly++) {
				layerwater += sub_water[subly] * Dz_sub; // add mm of water for this sublayer 
			}
			
			// Divide by the depth of the acrotelm layers to get the volumetric fraction
			Ftotal[ly] = layerwater / Dz_acro; 
			count++;
		}

		// Subtract ice content, water content below the permanent wilting point
		
		// ACROTELM

		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
			Frac_water[ly] = max(Ftotal[ly] - Frac_ice[ly] - Fpwp_ref[ly], 0.0);
			acro_aw += Frac_water[ly] * Dz[ly];
		}
	} 
	else { 

		// WTD <= 0.0
		// at or above the surface - standing water

		// Subtract ice content, water content below the permanent wilting point
		// For the calculations below, note that:
		// totaliceinlayer = Frac_ice + Fpwp_ref - Frac_water_belowpwp;
		// totalwaterinlayer = Frac_water + Frac_water_belowpwp;

		// ACROTELM - saturated
		for (int ly = IDX; ly < IDX + NACROTELM; ly++) {
			Frac_water[ly] = max(por[ly]-Frac_ice[ly]-Fpwp_ref[ly], 0.0);
			acro_aw += Frac_water[ly] * Dz[ly];
		}
	}

	// CATOTELM - assumed to be always saturated
	for (int ly = IDX + NACROTELM; ly < NLAYERS; ly++) {
		Frac_water[ly] = max(por[ly] - Frac_ice[ly] - Fpwp_ref[ly], 0.0);
		cato_aw += Frac_water[ly] * Dz[ly];
	}

	// Tidy up Frac_water
	for (int ly = 0; ly < NLAYERS; ly++) {
		if (Frac_water[ly] < Fpwp_ref[ly]) 
			Frac_water[ly] = 0.0;
	}

	// *** WATER CONTENT IN EACH SOIL LATER (FOR PLANT GROWTH) ***
	if (!update_layer_water_content(daynum)) 
		return;

	// UPDATE CO2 LEVEL IN ACROTELM
	update_acrotelm_co2(climate.co2);


	// *** UPDATE LIMITS ON MOSS PHOTOSYNTHESIS ***

	
	The cap is kept at 100% while WTP is above wtp_moss_upper ((Wania et al 2009b) have -150 mm), decreases linearly
	between wtp_moss_upper and -280mm  This has to do with the water
	content of the mosses as they are less able to photosynthesize
	as they dry out. NO MOSS DESSICATION
	
	
	double wtp_moss_upper = 0.0;	// Mosses can be too competitive when set to -150 mm, as in Wania et al. 2009b 
	double wtp_moss_lower = -280.0; // Wania et al. 2009b
	double moss_wtp_limit_lowerlimit = 0.3; 

	if (wtp[daynum] > wtp_moss_upper)
		dmoss_wtp_limit = 1.0; // No limit
	else if (wtp[daynum] > wtp_moss_lower)
		dmoss_wtp_limit = 1.0 - (-wtp[daynum] + wtp_moss_upper) * (1.0 - moss_wtp_limit_lowerlimit)/(wtp_moss_upper - wtp_moss_lower);
	else 
		dmoss_wtp_limit = moss_wtp_limit_lowerlimit; // between -300mm and -280mm
	

	// *** UPDATE LIMITS ON GRAMINOID PHOTOSYNTHESIS ***

	// LPJ-WHy code & comment:

	
	    IF (pft .EQ. 11 .AND. mwtp(m) .LT. -100) THEN
			! The following lines will reduce the mgpp of
			! flood-tolerant grasses when the water table drops
			! below -10cm. This should reflect that these grasses
			! are usually found where the water table stays close
			! to the surface for most of the time (i.e. in fens).
			! Note that mwtp is negative!
	
			mgpp(m,pft) = mgpp(m,pft) + mgpp(m,pft) *
			$                (mwtp(m) + 100) / 200 
	     ENDIF
	

	double wtp_graminoid_upper = -10; // -50.0;		// mm - Wania et al. 2009b have -100
	double wtp_graminoid_lower = -100; // -200.0;	// mm - Wania et al. 2009b have -300 mm
	double graminoid_wtp_limit_lowerlimit = 0.0; 

	// More restrictive than Wania et al., otherwise graminoids become too dominant as they never suffer from inundation stress
	// The cap drops from 1 as soon as wtp goes below wtp_graminoid_upper
	if (wtp[daynum] >= wtp_graminoid_upper)
		dgraminoid_wtp_limit = 1.0; // i.e. no limit!
	else if (wtp[daynum] > wtp_graminoid_lower)
		dgraminoid_wtp_limit = 1.0 - (-wtp[daynum] + wtp_graminoid_upper) * (1.0 - graminoid_wtp_limit_lowerlimit)/(wtp_graminoid_upper - wtp_graminoid_lower);
	else 
		dgraminoid_wtp_limit = graminoid_wtp_limit_lowerlimit; // between -300mm and -200mm


	// *** UPDATE OUTPUT VARIABLES ***

	// Reset monthly wtp and mmoss_wtp_limit averages on Jan 1
	if (daynum == 0) {
		for (int mth = 0; mth < 12; mth++) {	
			mwtp[mth] = 0.0;
		}
	}

	int mdays = date.ndaymonth[date.month];
	mwtp[date.month] += wtp[date.day] / (double)mdays;

	if (firstHydrologyCalc) 
		firstHydrologyCalc = false;

	if (date.day == Date::MAX_YEAR_LENGTH) {

		awtp = 0.0;
		for (int d = 0; d < 365; d++)
			awtp += wtp[d]/365.0; 
	}

	return; // ...if no problems
}*/


// return wcont for a certain layer
/*double Soil::get_layer_soil_water(int layer) const {

	if (layer < 0 || layer > NSOILLAYER-1 || layer != int(layer))
		return -999.0;
	else
		return wcont[layer];
}*/


// return wcont for a subset of soil layers
/*double Soil::get_soil_water(int layer1, int layer2) const {

	double total_available_water = 0.0; // available water for the thick soil layer (mm)
	double total_capacity = 0.0;

	if (layer1 < 0 || layer1 > NSOILLAYER || layer1 != int(layer1)) {
		fail("Soil::get_soil_water - bad layer1!\n");
		return -999;
	}

	if (layer2 < 0 || layer2 > NSOILLAYER || layer2 != int(layer2)) {
		fail("Soil::get_soil_water - bad layer2!\n");
		return -999;
	}

	if (layer2 <= layer1) {
		fail("Soil::get_soil_water - layer2 <= layer1\n");
		return -999;
	}

	for (int ly = layer1; ly < layer2; ly++) {
		
		if (wcont[ly] < 0.0 || (wcont[ly] > 1.0 && !negligible(wcont[ly] - 1.0, -12))) {
			fail("Soil::get_soil_water - bad wcont!\n");
			return -999;
		}

		total_available_water += wcont[ly] * soiltype.awc[ly]; // mm
		total_capacity += soiltype.awc[ly]; // mm
	}

	return total_available_water/ total_capacity;
}*/


// return wcont for the upper 50cm (upper true) or lower (100cm) soil layer (upper false)
/*double Soil::get_soil_water_upper() const {
	return get_soil_water(0, nsublayer1);
}*/

// return wcont for the lower (100cm) soil layer
/*double Soil::get_soil_water_lower() const {

	return get_soil_water(nsublayer1, NSOILLAYER);
}*/

// return copy of wcont
/*void Soil::copy_layer_soil_water_array(double wconttoreturn[NSOILLAYER]) {
		
	for (int ly = 0; ly < NSOILLAYER; ly++)
		wconttoreturn[ly] = wcont[ly];
}*/

// return wcont_evap
/*double Soil::get_layer_soil_water_evap() const {
		return wcont_evap;
}*/




// set wcont for a certain layer
/*void Soil::set_layer_soil_water(int layer, double newlayerwater) {

	if (layer < 0 || layer > NSOILLAYER-1 || layer != int(layer))
		fail("Soil::set_layer_soil_water - bad layer\n");
	else {
		wcont[layer] = newlayerwater;
	}
}*/

// set wcont_evap
/*void Soil::set_layer_soil_water_evap(double newevaplayerwater) {
		wcont_evap = newevaplayerwater;
}*/


/*void Soil::update_soil_water() {

	// DESCRIPTION
	// Update Frac_water and wcont_evap when there has been a change to wcont

	const double maxerr = 0.000001; // Max error allowed
															
	// Local storage
	double Faw_layer[NSOILLAYER];
	// available water for each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// ice each soil layer (mm)

	Wtot = 0.0;

	for (int ly = IDX; ly < IDX + nsublayer1 + nsublayer2; ly++) {

		oob_check_wcont(wcont[ly - IDX]);

		if (wcont[ly - IDX] < 0.0 || wcont[ly - IDX] > 1.0) {
			dprintf("Soil::update_soil_water - bad wcont!\n");
			return;
		}

		// For the following, note that, for each layer: wcont = Faw_layer / soiltype.awc
		// Also, wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			Faw_layer[ly - IDX] = wcont[ly - IDX] * soiltype.awc_peat[ly - IDX]; // mm
		else
			Faw_layer[ly - IDX] = wcont[ly - IDX] * soiltype.awc[ly - IDX]; // mm

 		ice_layer[ly - IDX] = Frac_ice[ly] * Dz[ly]; // mm

		// Check balance
		double layerwater = Faw_layer[ly - IDX] + ice_layer[ly - IDX]; // mm

		if (aw_max[ly - IDX] - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_soil_water - error in a soil layer's water balance!\n");
			return;
		}

		// Now update Frac_water
		Frac_water[ly] = Faw_layer[ly - IDX] / Dz[ly];

		// For wetlands:
		if (ly < IDX + NACROTELM)
			Wtot += (Frac_water[ly] + Frac_ice[ly] + Fpwp_ref[ly]) * Dz[ly]; // Wtot includes ice

	} // for loop (ly)

	// Update wcont_evap
	// wcont_evap is the available, nonfrozen water in the top two layers (20 cm), divided by the capacity of those layers to hold water
	if (patch.stand.is_highlatitude_peatland_stand())
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc_peat[0] + soiltype.awc_peat[1]);
	else
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc[0] + soiltype.awc[1]);

	if (wcont_evap < 0.0 || wcont_evap > 1.0) {
		dprintf("Soil::update_soil_water - error in the wcont_evap calculation!\n");
		return;
	}

	oob_check_wcont(wcont_evap);

	return; // ...if no problems
}*/


/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
/*void cnstep_full(int layer0, double Di[NLAYERS], double dz[NLAYERS], double surf_temp, double dt,
	double pad_dz[PAD_LAYERS], double temp[NLAYERS], double pad_temp[PAD_LAYERS],
	double ki[NLAYERS], double ci[NLAYERS]) {

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

	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Bedrock parameters from Chadburn et al. (2015)
	double Kbedrock = 8.6;			// W m-1 K-1
	double Cbedrock = 2100000.0;	// J m-3 K-1
	double Dbedrock = Kbedrock / Cbedrock * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)
	int layer, lidx;
	int active_layers;

	// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson matrix.
	long double diag[active_layersmax];
	long double left[active_layersmax];
	long double right[active_layersmax];

	// Right hand side vector for Crank-Nicholson scheme equations.
	long double rhs[active_layersmax];

	// Solution vector for Crank-Nicholson scheme equations.
	long double solution[active_layersmax];

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	for (int l = 0; l < active_layersmax; l++) {
		diag[l] = 0.0;
		left[l] = 0.0;
		right[l] = 0.0;
		rhs[l] = 0.0;
		solution[l] = 0.0;
	}

	// --- CODE STARTS HERE ---

	// Not all layers in the input arrays are actually used.  Calculate
	// how many layers we have to deal with here.
	MISSING_VALUE = -9999.0;

	active_layers = NLAYERS - layer0 + PAD_LAYERS;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers - 1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	// lidx runs from 1 to active_layers == 
	for (lidx = 1; lidx <= active_layers; lidx++) {

		// Deal with different layer counting schemes.
		layer = lidx + layer0 - 1;
		// Minimum is layer0
		// Maximum is NLAYERS + PAD_LAYERS-1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  Diffusion
		// constants in the padding layers are calculated in a sensible way.

		// D+
		double small = 0.000001;

		// D+
		const int first_bedrock_padding_layer = PAD_LAYERS; // 1 bedrock layer

		if (layer == NLAYERS + PAD_LAYERS - 1)
			dplus = 0.0;				// BC2 - Bottom layer diffusion clamped to 0
		else if (layer >= NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Dbedrock;			// The lowest padding layers are treated as bedrock
	 	else
			dplus = (1.0 / ci[layer]) * (ki[layer] * dz[layer] + ki[layer + 1] * dz[layer + 1]) / (dz[layer] + dz[layer + 1]);

		// D-

		if (layer == layer0)
			dminus = Di[layer];			// top layer
		else if (layer > NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dminus = (1.0 / ci[layer]) * (ki[layer] * dz[layer] + ki[layer - 1] * dz[layer - 1]) / (dz[layer] + dz[layer - 1]);

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// taking account of padding layers and end cases.

		// --- HERE ---
		if (layer < NLAYERS) {
			 // all soil layers
			dzhere = dz[layer];
			there = temp[layer];
		}
		else {
			// all padding layers
			dzhere = pad_dz[layer - NLAYERS];
			there = pad_temp[layer - NLAYERS];
		}

		// --- PLUS ---
		if (layer < NLAYERS - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			tplus = temp[layer + 1];
		}
		else if (layer < NLAYERS + PAD_LAYERS - 1) {
			// bottom soil layer and all 
			// padding layers apart from the very bottom padding layer 
			dzplus = pad_dz[layer - NLAYERS + 1];
			tplus = pad_temp[layer - NLAYERS + 1];
		}
		else {
			// bottom padding layer
			dzplus = pad_dz[PAD_LAYERS - 1];
			tplus = pad_temp[PAD_LAYERS - 1];
		}

		// --- MINUS ---
		if (layer == layer0) {
			// top layer
			dzminus = dz[layer];
			tminus = temp[layer];
		}
		else if (layer <= NLAYERS) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			tminus = temp[layer - 1];
		}
		else {
			// all padding layers apart from the top padding layer
			dzminus = pad_dz[layer - NLAYERS - 1];
			tminus = pad_temp[layer - NLAYERS - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lidx == 1)
			diag[0] = 1.0; // BC1 - top layer should be (1,0,...,0)
		else
			diag[lidx - 1] = 1.0 + Cplus + Cp_minus;

		// LEFT & RIGHT
		if (lidx < active_layers) {

			if (lidx > 1)
				right[lidx - 1] = -Cplus;
			else
				right[lidx - 1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx - 1] = -Cp_minus;
		}

		if (lidx == active_layers) left[lidx - 1] = -Cp_minus;

		// RHS
		// Calculate right hand side vector values.
		if (lidx == 1)
			rhs[0] = surf_temp;
		else if (lidx == active_layers) // Cplus == 0 here anyway
			rhs[lidx - 1] = (1.0 - Cp_minus) * there + Cp_minus * tminus;
		else
			rhs[lidx - 1] = (1.0 - Cplus - Cp_minus) * there +
			Cplus * tplus + Cp_minus * tminus;
	} // end for

	//   SOLVE TRIDIAGONAL SYSTEM

	tridiag(active_layers, left, diag, right, rhs, solution);

	// Numerical tests

	if (DEBUG_SOIL_TEMPERATURE) {

		int testrow = 10;
		double checksum = left[testrow] * solution[testrow - 1] + diag[testrow] * solution[testrow] +
			right[testrow] * solution[testrow + 1] - rhs[testrow];

		if (fabs(checksum) > 0.0000000001)
			dprintf("%s\n", "Bad checksum after tridiag - test1");
	}

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the temperature array.  Also update
	// the padding layer temperatures so they can be used for the
	// next timestep.
	for (lidx = 0; lidx<NLAYERS; lidx++) {

		// remove temps < 1e-12.
		if (verysmall(solution[lidx])) solution[lidx] = 0.0;

		if (lidx < layer0)
			temp[lidx] = MISSING_VALUE;
		else
			temp[lidx] = solution[lidx - layer0];
	}

	for (lidx = 0; lidx<PAD_LAYERS; lidx++) {
		pad_temp[lidx] = solution[lidx + NLAYERS - layer0];
	}
}*/

/// Crank-Nicholson timestepper algorithm for temperature diffusion equation.
/*void cnstep(int layer0, double Di[NLAYERS], double dz[NLAYERS], double surf_temp, double dt,
	double pad_dz[PAD_LAYERS], double temp[NLAYERS], double pad_temp[PAD_LAYERS]) {
	// DESCRIPTION
	// Crank-Nicholson timestepper for temperature diffusion equation.
	// This routine performs a single timestep (of length dt) of the heat
	// diffusion equation
	//
	//                        6 T    6  /      6 T \
	//                        --- = --- | D(z) --- |
	// 						  6 t   6 z \      6 z /
	//
	// where here, T(z, t) is the temperature, the '6's represent partial
	// differentiation, and D(z) is the (depth-dependent) diffusion
	// constant.  The boundary conditions are a prescribed surface
	// temperature, and no heat flow at the bottom of the solution domain.
	//
	// INPUT PARAMETERS:
	//
	// NLAYERS      Total number of layers represented in temp, d and
	//              dz arrays.  Not all of the layers have to be active
	//              at one time, but this parameter is used to constrain
	//              the sizes of the depth-based arrays.
	//
	// layer0       Index of the first active layer.  The diffusion
	//              equation will be solved for all layers in the range
	//              layer0 - NLAYERS inclusive, and the surface temperature
	//              boundary condition will be applied to the top of layer
	//              index layer0.
	//
	// temp         The temperatures in each layer at the start of the
	//              timestep.
	//
	// Di           The thermal diffusivities in each layer. (m2/day)
	//
	// dz           The thickness of each layer.
	//
	// surf_temp    The surface temperature boundary condition.
	//
	// dt           The timestep.
	//
	// PAD_LAYERS   Number of layers to use below computation domain to
	//              deal with bottom boundary condition.
	//
	// pad_dz       Padding layer thicknesses (derived from calc_padding).
	//
	// pad_temp     Temperature in padding layers.
	//
	// NOTE ON UNITS: no units are specified here for any of these
	// inputs; the only constraint on units is that the layer
	// thicknesses, the diffusion constant and the timestep should be
	// measured in consistent units.  For instance, if the layer
	// thicknesses are measured in metres, and the timestep in days, {
	// the diffusion constants should have units of m2 / day.
	//
	// OUTPUT PARAMETERS
	//
	// temp         The temperatures in each layer at the end of the
	//              timestep.  Only layers layer0 - NLAYERS will have
	//              valid temperature values after a  to cnstep.  All
	//              layers with index less than layer0 will have missing
	//              value temperatures set (a value of -1.0E35).
	//
	// pad_temp     Temperature in padding layers at the end of the
	//              timestep.


	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Bedrock parameters from Chadburn et al. (2015)
	double Kbedrock = 8.6;			// W m-1 K-1
	double Cbedrock = 2100000.0;	// J m-3 K-1
	double Dbedrock = Kbedrock / Cbedrock * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)
	int layer, lidx;
	int active_layers;

	// Diffusion constants averaged over adjacent layers.
	double dplus;
	double dminus;

	// Layer-dependent weighting factors in Crank-Nicholson scheme.
	double dz_factor;
	double Cplus;
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double there, tminus, tplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson
	// matrix.
	long double diag[active_layersmax];
	long double left[active_layersmax];
	long double right[active_layersmax];

	// Right hand side vector for Crank-Nicholson scheme equations.
	long double rhs[active_layersmax];

	// Solution vector for Crank-Nicholson scheme equations.
	long double solution[active_layersmax];

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = 0.0;
	there = tminus = tplus = 0.0;

	// --- CODE STARTS HERE ---

	// Not all layers in the input arrays are actually used.  Calculate
	// how many layers we have to deal with here.
	MISSING_VALUE = -9999.0;

	active_layers = NLAYERS - layer0 + PAD_LAYERS;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers - 1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	for (lidx = 1; lidx <= active_layers; lidx++) {

		if (lidx == 1) { // moved from above
			diag[lidx-1] = 0.0;
			left[lidx-1] = 0.0;
			right[lidx-1] = 0.0;
			rhs[lidx-1] = 0.0;
			solution[lidx-1] = 0.0;
		}

		// Deal with different layer counting schemes.
		layer = lidx + layer0 - 1;
		// Minimum is layer0
		// Maximum is NLAYERS + PAD_LAYERS-1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  Diffusion
		// constants in the padding layers are calculated in a sensible way.

		// D+
		const int first_bedrock_padding_layer = 3; // 3 layers

		if (layer == NLAYERS + PAD_LAYERS - 1)
			dplus = 0.0;				// BC2 - Bottom layer diffusion clamped to 0
		else if (layer >= NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer. 
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dplus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dplus = 0.5 * (Di[layer] + Di[layer + 1]);

		// D-

		if (layer == layer0)
			dminus = Di[layer]; // top layer
		else if (layer > NLAYERS - 1 && layer < NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Di[NLAYERS - 1];	// The top padding layers have the same Di as the bottom soil layer.
		else if (layer >= NLAYERS - 1 + first_bedrock_padding_layer)
			dminus = Dbedrock;			// The lowest padding layers are treated as bedrock
		else
			dminus = 0.5 * (Di[layer] + Di[layer - 1]); // soil layers

		// Extract sensible values to use for temperature and diffusivity
		// of the current layer, the layer above and the layer below,
		// taking account of padding layers and end cases.

		// --- HERE ---
		if (layer < NLAYERS) {
			// all soil layers
			dzhere = dz[layer];
			there = temp[layer];
		}
		else {
			// all padding layers
			dzhere = pad_dz[layer - NLAYERS];
			there = pad_temp[layer - NLAYERS];
		}

		// --- PLUS ---
		if (layer < NLAYERS - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			tplus = temp[layer + 1];
		}
		else if (layer < NLAYERS + PAD_LAYERS - 1) {
			// bottom soil layer and all 
			// padding layers apart from the very bottom padding layer 
			dzplus = pad_dz[layer - NLAYERS + 1];
			tplus = pad_temp[layer - NLAYERS + 1];
		}
		else {
			// bottom padding layer
			dzplus = pad_dz[PAD_LAYERS - 1]; // PaulM added -1
			tplus = pad_temp[PAD_LAYERS - 1]; // PaulM added -1
		}

		// --- MINUS ---
		if (layer == layer0) {
			// top layer
			dzminus = dz[layer];
			tminus = temp[layer];
		}
		else if (layer <= NLAYERS) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			tminus = temp[layer - 1];
		}
		else {
			// all padding layers apart from the top padding layer
			dzminus = pad_dz[layer - NLAYERS - 1];
			tminus = pad_temp[layer - NLAYERS - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.

		// DIAG
		if (lidx == 1)
			diag[0] = 1.0; // BC1 - top layer should be (1,0,...,0)
		else
			diag[lidx - 1] = 1.0 + Cplus + Cp_minus;

		// LEFT & RIGHT
		if (lidx < active_layers) {

			if (lidx > 1)
				right[lidx - 1] = -Cplus;
			else
				right[lidx - 1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx - 1] = -Cp_minus;
		}

		if (lidx == active_layers) left[lidx - 1] = -Cp_minus;

		// RHS
		// Calculate right hand side vector values.
		if (lidx == 1)
			rhs[0] = surf_temp;
		else if (lidx == active_layers) // Cplus == 0 here anyway
			rhs[lidx - 1] = (1.0 - Cp_minus) * there + Cp_minus * tminus;
		else
			rhs[lidx - 1] = (1.0 - Cplus - Cp_minus) * there +
			Cplus * tplus + Cp_minus * tminus;
	} // end for

	// SOLVE TRIDIAGONAL SYSTEM

	tridiag(active_layers, left, diag, right, rhs, solution);

	// Numerical test

	if (DEBUG_SOIL_TEMPERATURE) {

		int testrow = 10;
		double checksum = left[testrow] * solution[testrow - 1] + diag[testrow] * solution[testrow] +
			right[testrow] * solution[testrow + 1] - rhs[testrow];

		if (fabs(checksum) > 0.0000000001)
			dprintf("%s\n", "Bad checksum after tridiag - test1");
	}

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the temperature array.  Also update
	// the padding layer temperatures so they can be used for the
	// next timestep.
	for (lidx = 0; lidx<NLAYERS; lidx++) {

		if (lidx < layer0)
			temp[lidx] = MISSING_VALUE;
		else
			temp[lidx] = solution[lidx - layer0];

		if (lidx < PAD_LAYERS) // was in the loop below
			pad_temp[lidx] = solution[lidx + NLAYERS - layer0];
	}
}*/


/*void calc_padding(double dz, double pad_dz[PAD_LAYERS]) {

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

	// Do Newton-Raphson solve to find geometric sequence factor.
	while (fabs(k_guess - k_last) > 1E-2) {
		k_last = k_guess;

		func = pow(k_guess, (PAD_LAYERS + 1)) - k_guess * (1.0 + pad_depth / dz) + pad_depth / dz;
		dfunc = (PAD_LAYERS + 1) * pow(k_guess, PAD_LAYERS) - (1.0 + pad_depth / dz);
		k_guess = k_guess - func / dfunc;
	}

	// Build geometric sequence of padding layer depths.
	for (layer = 0; layer<PAD_LAYERS; layer++) {
		pad_dz[layer] = dz * pow(k_guess, layer + 1);
	}
}*/


void Soilcolumn::init_hydrology_variables() {

	// DESCRIPTION
	// Initialisation of hydrology. 
	// Called ONCE, at the beginning of every stand's simulation
	// The initialisation relies on the following definitions:

	// sub_water[NSUBLAYERS_ACRO] - the volumetric water content in the NSUBLAYERS_ACRO of the acrotelm
	// whc[NSOILLAYER] - available water holding capacity of soil layers [0=upper layer] [mm], taking into
	// account the unavailability of frozen water. Default value: soiltype.awc[]
	// whc_evap - available water holding capacity of evap soil layers [mm], taking into
	// account the unavailability of frozen water. Default value: (2/5) * soiltype.awc[]
	// aw_max[NSOILLAYER] - Max water (mm) that can be held in each layer
	// alwhc[NLAYERS] - the volumetric liquid water content. A fraction. 
	// Considers the entire (awc + Fpwp) volumetric water content MINUS the ice fraction. Updated daily.
	// alwhc_init[NLAYERS] - the initial volumetric liquid water content. A fraction. Considers the entire
	// (awc + Fpwp) volumetric water content MINUS the ice fraction.

	// Procedure:
	// Initialise whc[] and the total volumetric liquid water content, alwhc[]
	// Save the initial values - alwhc_init 
	// Initialise aw_max[], the max water (mm) that can be held in each layer

	for (int ly = 0; ly < nsoillayers; ly++) {
		
		Soillayer& soillayer = soillayers[ly];

		if (patch.stand.is_highlatitude_peatland_stand()) {

			// Override hydrology properties for peatland cells
			if (ly < NACROTELM) {
				soillayer.alwhc = acro_por - peat_wp;
				soillayer.awc_peat = soillayer.depth * (acro_por - peat_wp); // mm
				soillayer.wsats_peat = acro_por * soillayer.depth; // DW_COMMENT SHOULD FIX WSATS AND REST FOR PEAT
			}
			else {
				soillayer.alwhc = cato_por - peat_wp;
				soillayer.awc_peat = soillayer.depth * (cato_por - peat_wp); // mm
				soillayer.wsats_peat = cato_por * soillayer.depth;
			}
			
			soillayer.whc= soillayer.awc_peat; // mm
			
			// Max water (mm) that can be held in each layer, above the pwp
			soillayer.aw_max = soillayer.wsats_peat * soillayer.depth; // mm
		}
		else {

			// Not a high latitude peatland
			soillayer.alwhc = soillayer.awc / soillayer.depth;
			soillayer.whc = soillayer.awc; // mm

			// Max water (mm) that can be held in each layer, above the pwp
			soillayer.aw_max = soillayer.wsats * soillayer.depth; // mm
		}
		
		soillayer.alwhc_init = soillayer.alwhc;
	}
}


/*void Soil::update_acrotelm_co2(double atmo_co2) {

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

	acro_co2 = min(PORE_WATER_CO2, PORE_WATER_CO2 + ((PORE_WATER_CO2 - atmo_co2) / acro_depth) * awtp);
}*/


/*bool Soil::update_layer_water_content(int day) {

	// DESCRIPTION
	// Private helper method to update the soil water contents in the soil layers. 
	// Called each day.

	// whc[] is water above the pwp which hasn't frozen [mm]
	// Its maximum value is awc[] (mm), i.e. whc[later] <= awc[layer], for all layers

	const double maxerr = 0.000001; // Max error allowed

	// Local storage
	double Faw_layer[NSOILLAYER];
	// available water for each soil layer (mm)
	double ice_layer[NSOILLAYER];
	// mm of ice in each soil layer (mm)

	// Peatland total soil water content, incl ice and water below the pwp
	Wtot = 0.0; // mm

	for (int ly = IDX; ly < IDX + nsublayer1 + nsublayer2; ly++) {

		// Update the total available volumetric water content (mm) of the layers
		Faw_layer[ly-IDX] = Frac_water[ly] * Dz[ly]; // mm
		ice_layer[ly-IDX] = Frac_ice[ly] * Dz[ly]; // mm

		// Check balance
		double layerwater = Faw_layer[ly-IDX] + ice_layer[ly-IDX]; // mm

		if (aw_max[ly-IDX] - layerwater < -1.0 * maxerr) {
			dprintf("Soil::update_layer_water_content - error in a soil layer's water balance!\n");
			return false;
		}

		// Subtract the ice/frozen fraction from the initial value of the capacity for this layer to hold water
		alwhc[ly-IDX] = alwhc_init[ly-IDX] - Frac_ice[ly]; 

		// For wetlands:
		// total ice in layer = Frac_ice[i] + Fpwp_ref[i] - Frac_water_belowpwp[i];
		// total water in layer = Frac_water[i] + Frac_water_belowpwp[i];
		// So, the total ice + water = Frac_ice[i] + Fpwp_ref[i] + Frac_water[i] 
		if (ly < IDX + NACROTELM)
			Wtot += (Frac_water[ly] + Frac_ice[ly] + Fpwp_ref[ly]) * Dz[ly]; // Wtot includes ice

		if (patch.stand.is_highlatitude_peatland_stand())
			whc[ly-IDX] = max(0.0, aw_max[ly-IDX] - ice_layer[ly-IDX]); // mm
		else
			whc[ly-IDX] = max(0.0, soiltype.awc[ly-IDX] - ice_layer[ly-IDX]); // mm 
		
		// Avoid dividing by small numbers later
		if (fabs(whc[ly-IDX]) < maxerr)
			whc[ly-IDX] = 0.0;

		// wcont is the available, nonfrozen water in the layer, divided by the capacity of that layer to hold water. 
		if (patch.stand.is_highlatitude_peatland_stand())
			wcont[ly-IDX] = Faw_layer[ly-IDX] / soiltype.awc_peat[ly-IDX];
		else
			wcont[ly-IDX] = Faw_layer[ly-IDX] / soiltype.awc[ly-IDX];

		oob_check_wcont(wcont[ly - IDX]);

		if (wcont[ly-IDX] < 0.0 || wcont[ly-IDX] > 1.0) {
			dprintf("Soil::update_layer_water_content - error in the wcont calculation!\n");
			return false;
		}
	} // for loop (ly)

	// wcont_evap is the available, nonfrozen water in the tope two layers (20 cm), divided by the capacity of those layers to hold water
	if (patch.stand.is_highlatitude_peatland_stand())
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc_peat[0] + soiltype.awc_peat[1]);
	else 
		wcont_evap = (Faw_layer[0] + Faw_layer[1]) / (soiltype.awc[0] + soiltype.awc[1]);
	
	oob_check_wcont(wcont_evap);

	if (wcont_evap < 0.0 || wcont_evap > 1.0) {
		dprintf("Soil::update_layer_water_content - error in the wcont_evap calculation!\n");
		return false;
	}

	return true; // ...if no problems
}*/


// Private helper function
/*void Soil::update_from_yesterday() {

	// DESCRIPTION
	// Update T_soil and Frac_ice using this function. 
	// More variables could be added here.

	// Frac_ice_yesterday and T_soil_yesterday record yesterday's values, not just those on Dec 31.
	for (int i = 0; i<NLAYERS; i++) {
		T_soil[i] = T_soil_yesterday[i];
		Frac_ice[i] = Frac_ice_yesterday[i];
	}
}*/


/*void Soil::update_snow_properties(const int& daynum, const double& dailyairtemp, double& Dsnow, double& Csnow, double& Ksnow) {

	// DESCRIPTION
	// Hillel (1982) gives values for porosity, Csnow and Ksnow.       
	// Here we update snow density and use that to estimate Csnow and Ksnow 

	// Csnow:		heat capacity of snow [J m-3 K-1]
	// Ksnow:		thermal conductivity of snow [W m-1 K-1]
	// Dsnow:		thermal conductivity of snow [mm2 day-1]

	bool snow_active_old = snow_active;

	// Are there snow layers?
	if (snowpack > 1) // mm water
		snow_active = true;
	else
		snow_active = false;

	// No snow, so return
	if (!snow_active) {
		snow_active_layers = 0;
		return;
	}

	// new snow?
	if (!snow_active_old)
		snowdens = snowdens_start;

	// Update snowdens using the Wania et al. simple snow compaction scheme

	// We assume fresh snow unless snow_days > 0.75 * 365;
	if (date.year >= 1 && snow_days >= 0.75 * snow_days_prev) {
		snowdens = snowdens_start + (snowdens_end - snowdens_start) / (snow_days_prev * 0.25) * (snow_days - 0.75 * snow_days_prev);
		 
		// Restrict density to the range [snowdens_end,snowdens_start]
		snowdens = max(min(snowdens, snowdens_end), snowdens_start);
	}
	else
		snowdens = snowdens_start;

	// Override the above with a fixed snow density?
	if (snowdensityconstant) {
		snowdens = 250.0; // kg m-3, JSBACH value from Ekici et al (2015)
	}

	// Note that K/C = D has units mm2/day, as above
	// Wania values
	// Csnow from Fukusako, Eq. 2. It's the same as for ice
	static const double ice_density = 917.0; // kg m-3
	Csnow = (0.185 + 0.689 * (K2degC + dailyairtemp) * 0.01) * J_PER_KJ; // J kg-1 K-1
	Csnow *= ice_density; // conversion to volumetric heat capacity (J m-3 K-1) - should be 1,900,000 J m-3 K-1 approx.

	double snowdens_scaled = snowdens / rho_H2O;
	// From Sturm et al. 1997, and in Ling and Zhang, 2006
	if (snowdens > 156.0)
		Ksnow = 0.138 - 1.01 * snowdens_scaled + 3.233 * snowdens_scaled * snowdens_scaled;
	else
		Ksnow = 0.023 + 0.234 * snowdens_scaled;
	 
	Dsnow = Ksnow / Csnow * SECS_PER_DAY * MM2_PER_M2; // mm2 day-1

	// Unit conversions - needed in cnstep_full:
	// Conversion from (W m-1 K-1) TO (J day-1 mm-1 K-1) 
	Ksnow *= SECS_PER_DAY * M_PER_MM;
	// Conversion from (J m-3 K-1) TO (J mm-3 K-1)
	Csnow *= M3_PER_MM3;
}*/


// A valid layer number?
/*bool Soil::valid_layer_num(const int& ngr) const {

	if (ngr != (SOILDEPTH_UPPER + SOILDEPTH_LOWER) / Dz_soil)
		return false;
	else
		return true;
}*/


/*void Soil::update_soil_diffusivities(const int& daynum, bool ansoln) {

	// DESCRIPTION
	// Called daily to update diffusivities, heat capacities and thermal conductivities of the 
	// soil layers below the mixed/moss and snow layers 

	double lKfwater;					// log of therm. cond. for the water fraction
	double lKfice;						// log of therm. cond. for the ice fraction
	double lKforg;						// log of therm. cond. for the org. fraction
	double lKfpeat;						// log of therm. cond. for the peat fraction
	double lKfmin;						// log of therm. cond. for the min. fraction
	double ThermCond_total[NLAYERS];	// total thermal condunctivity
	double HeatCapacity_total[NLAYERS];	// total heat capacity
	double Ftot;						// total of solid fraction in soil

	for (int i = IDX; i<NLAYERS; i++) {

		if (ansoln) Frac_water[i] = 0.0; // No water if analytic soln. sought

		// take water and ice below pwp into account
		double totaliceinlayer = Frac_ice[i] + Fpwp_ref[i] - Frac_water_belowpwp[i];
		double totalwaterinlayer = Frac_water[i] + Frac_water_belowpwp[i];

		Frac_air[i] = por[i] - totalwaterinlayer - totaliceinlayer;

		if (Frac_air[i] < -0.00001)
			fail("Illegal Frac_air value in Soil::update_soil_diffusivities: %g", Frac_air[i]);

		// Catch rounding errors for Frac_air
		if (Frac_air[i] < 0.0 && Frac_air[i] > -0.00001)
			Frac_air[i] = 0.0;

		if (ansoln) Frac_air[i] = 0.0; // No air if analytic soln. sought

		HeatCapacity_total[i] = Frac_min[i] * Cp_min + Frac_org[i] * Cp_org + totaliceinlayer * Cp_ice +
			 totalwaterinlayer * Cp_water + Frac_peat[i] * Cp_peat + Frac_air[i] * Cp_air;

		// Update the THERMAL CONDUCTIVITY K - which follows
		// Granberg et al. 1999, who refers to Farouki, 1986

		Ftot = totalwaterinlayer + totaliceinlayer + Frac_org[i] + Frac_min[i] + Frac_peat[i]; 

		if (Ftot + Frac_air[i] > 1.000001)
			fail("Illegal Ftot (%g) + Fair (%g) value in Soil::update_soil_diffusivities: %g", Ftot, Frac_air[i], Ftot +Frac_air[i]);

		lKforg = (Frac_org[i] / Ftot) * lKorg;
		lKfpeat = (Frac_peat[i] / Ftot) * lKpeat;
		lKfmin = (Frac_min[i] / Ftot) * lKmin;
		lKfwater = (totalwaterinlayer / Ftot) * lKwater;
		lKfice = (totaliceinlayer / Ftot) * lKice;

		ThermCond_total[i] = Frac_air[i] * Kair + Ftot * exp(lKfwater + lKfice + lKforg + lKfmin + lKfpeat);

		// Update THERMAL DIFFUSIVITY for each SOIL layer
		// multiply by 86400 to get from m2 s-1 to m2 d-1
		// multiply by 1E6 to get from m2 d-1 to mm2 d-1
		Di[i] = ThermCond_total[i] / HeatCapacity_total[i] * SECS_PER_DAY * MM2_PER_M2;	

		// Conversion from (W m-1 K-1) TO (J day-1 mm-1 K-1)
		Ki[i] = ThermCond_total[i] * SECS_PER_DAY * M_PER_MM;
		// Conversion from (J m-3  K-1) TO (J mm-3 K-1)
		Ci[i] = HeatCapacity_total[i] * M3_PER_MM3;

		// Note that K/C = D has units mm2/day, as above
	} // i 
}*/


/*void Soil::update_ice_fraction(const int& daynum, const int& MIDX) {

	// DESCRIPTION
	// Calculate WATER AND ICE fraction in each layer

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


	double Ffreez;	// fraction of water which freezes
	double Fthaw;	// fraction of ice which thaws
	double energy;	// Energy per volume [J m-3]

	const double SMALLVOLFRAC = 0.0000001;	// Minimum volume fraction [m-3] allowed

	double delta_T;	// temperature difference (current - previous)

	for (int i = MIDX; i<NLAYERS; i++) {

		if (firstTempCalc)
			delta_T = 0.0;

		double heat_capacity_layer = Ci[i] * MM3_PER_M3; // J m-3 K-1

		double waterabovepwp = Frac_ice[i] + Frac_water[i]; 
		double ice_thislayer = Frac_ice[i] + (Fpwp_ref[i]-Frac_water_belowpwp[i]);
		double water_thislayer = Frac_water[i] + Frac_water_belowpwp[i];
		// The total amount of ice and water should not change in this routine
		double wtotal_initial = ice_thislayer + water_thislayer;

		// Add the water (not ice) below wp before calculating phase change
		Frac_water[i] += Frac_water_belowpwp[i];
		Frac_ice[i] += (Fpwp_ref[i]-Frac_water_belowpwp[i]);

		// Initialise temperatures on first day
		if (daynum == 0 && date.year == FIRST_FREEZE_YEAR)
			T_old[i] = T_soil[i];

		delta_T = T_soil[i] - T_old[i];

		// FREEZING
		if (T_soil[i] < 0.0 && T_soil[i] < T_old[i]) { // Wania conditions 
			Fthaw = 0.0;
			if (Frac_water[i] > 0.0){ //&& T_old[i] < 1.0) {
				energy = Cp_water * fabs(delta_T); //energy in J m-3 
				Ffreez = energy / Lheat; //unitless fraction

				if (Frac_water[i] > Ffreez) {
					Frac_ice[i] += Ffreez;
					Frac_water[i] -= Ffreez;

					// How much of the water is below pwp?
					if (Frac_water[i] < Fpwp_ref[i]) {
						Frac_water_belowpwp[i] = Frac_water[i];  
						Frac_water[i] = 0.0; // So there is liquid water, but it's all below the pwp
					}
					else {
						Frac_water_belowpwp[i] = Fpwp_ref[i];  
						Frac_water[i] -= Fpwp_ref[i]; // So there is liquid water above the pwp
					}

					T_soil[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_water[i] < SMALLVOLFRAC)
							Frac_water[i] = 0.0;
						if (Frac_ice[i] < SMALLVOLFRAC) 
							Frac_ice[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}
				}
				else {

					// It is cold enough to freeze all the water and decrease the temperature below 0, 
					// since the rest of the water gets frozen now. However, it should never result in 
					// temperatures colder than the values coming from the Crank-Nicholson scheme.

					Frac_ice[i] += Frac_water[i];
					
					// Update soil temperature in this layer
					T_soil[i] = -(Ffreez - Frac_water[i]) * Lheat / Cp_water;

					// Note that this gives the same result as the following method:
					// double temp_rise = Frac_water[i] * Lheat / Cp_water; // temp rise resulting from freezing water.
					// T_soil[i] += temp_rise;
					
					T_soil[i] = min(0.0, T_soil[i]); // Ensures we don't go back over 0 degrees

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_ice[i] < SMALLVOLFRAC)
							Frac_ice[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}

					Frac_water[i] = 0.0;
					Frac_water_belowpwp[i] = 0.0; // All water is frozen, even that below the pwp
				}
			}
			else {
				Ffreez = 0.0;
				Frac_water[i] -= Frac_water_belowpwp[i];	
			}	
		} // end of FREEZING

		// THAWING
		else if (T_soil[i] >= 0.0  && T_soil[i] > T_old[i]) { // Wania conditions
			Ffreez = 0.0;
			if (Frac_ice[i] > 0.0) {
				energy = Cp_water * fabs(delta_T); //energy in J m-3
				Fthaw = energy / Lheat; //unitless fraction
				if (Frac_ice[i] > Fthaw) {

					// There is some ice left, even after melt
					Frac_water[i] += Fthaw;

					// How much of the water is below the pwp now?
					if (Frac_water[i] < Fpwp_ref[i]) {
						Frac_water_belowpwp[i] = Frac_water[i];  
						Frac_water[i] = 0.0; // So there is liquid water, but it's all below the pwp
					}
					else {
						Frac_water_belowpwp[i] = Fpwp_ref[i];  
						Frac_water[i] -= Fpwp_ref[i]; // So there is liquid water above the pwp
					}

					Frac_ice[i] -= Fthaw;
					T_soil[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_ice[i] < SMALLVOLFRAC) 
							Frac_ice[i] = 0.0;
						if (Frac_water[i] < SMALLVOLFRAC) 
							Frac_water[i] = 0.0;
						if (Frac_water_belowpwp[i] < SMALLVOLFRAC) 
							Frac_water_belowpwp[i] = 0.0;
					}
				}
				else {
					// There is more energy available than required for
					// melting the rest of the ice. Therefore the
					// left-over energy will increase the temperature.

					Frac_water[i] += Frac_ice[i];
					Frac_water_belowpwp[i] = Fpwp_ref[i]; // All water is unfrozen
					Frac_water[i] -= Frac_water_belowpwp[i]; // Only store liquid water above the pwp

					// Update soil temperature in this layer
					T_soil[i] = (Fthaw - Frac_ice[i]) * Lheat / Cp_water;					

					// Note that this gives the same result as the following method:
					// double temp_reduction = Frac_ice[i] * Lheat / Cp_water; // temp reduction resulting from thawing ice.
					// T_soil[i] -= temp_reduction;

					T_soil[i] = max(0.0, T_soil[i]); // Ensures we don't go back under 0 degrees

					Frac_ice[i] = 0.0;

					// Remove tiny amounts of water for numerical stability
					if (DEBUG_SOIL_WATER) {
						if (Frac_water[i] < SMALLVOLFRAC) 
							Frac_water[i] = 0.0;
					}
				}
			}
			else {
				Fthaw = 0.0;
				Frac_ice[i] = 0.0;
				Frac_water[i] -= Frac_water_belowpwp[i];	
			}
		} // end of THAWING
		else {
			// Neither freezing nor thawing, so revert to Frac_water
			Frac_water[i] -= Frac_water_belowpwp[i];	
		}

		// Subtract ice below pwp from Frac_ice 
		Frac_ice[i] -= (Fpwp_ref[i]-Frac_water_belowpwp[i]);

		if (Frac_ice[i] < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid ice volume content %g\n", Frac_ice[i]);

		if (Frac_water[i] < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid volumetric water content %g\n", Frac_water[i]);

		// Remove tiny amounts of ice and water for numerical stability
		if (DEBUG_SOIL_WATER) {
			// Tidy up due to rounding errors
			if (Frac_ice[i] < 0.0 && Frac_ice[i] > -1.0 * SMALLVOLFRAC)
				Frac_ice[i] = 0.0;

			if (Frac_water[i] < 0.0 && Frac_water[i] > -1.0 * SMALLVOLFRAC)
				Frac_water[i] = 0.0;
		}

		// water balance checks:
		double ice_thislayer_new = Frac_ice[i] + (Fpwp_ref[i]-Frac_water_belowpwp[i]);
		double water_thislayer_new = Frac_water[i] + Frac_water_belowpwp[i];
		double wtotal_final = ice_thislayer_new + water_thislayer_new;

		if ((fabs(wtotal_initial - wtotal_final)) > SMALLVOLFRAC)
			fail("Water mass balance during freezing or thawing in Soil::update_ice_fraction. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);

		// Update T_old after possible phase change
		T_old[i] = T_soil[i];

	} // end of loop through layers
}*/


/*void Soil::update_layer_fractions(const int& daynum, const int& mixedl, const int& MIDX) {

	// DESCRIPTION
	// Inializes and updates fractions of water, mineral content, peat content, organic content, as well as porosity.

	if (patch.stand.is_highlatitude_peatland_stand()) {

		if (firstTempCalc) {

			for (int i = IDX; i < IDX + NACROTELM; i++) {

				// ACROTELM

				Dz[i] = Dz_acro;
				por[i] = acro_por;
				Frac_org[i] = 0.0;
				Frac_min[i] = 0.0;	
				Frac_peat[i] = 1.0 - (acro_por + Fgas);

				Fpwp_ref[i] = peat_wp;  // Wisser et al. 2011 Very low in peatlands
				Frac_water_belowpwp[i] = peat_wp;
				Frac_water[i] = acro_por - peat_wp; // i.e. initially saturated

				soiltype.awc_peat[i-IDX] = (por[i] - peat_wp) * Dz[i];
			}

			for (int i = IDX + NACROTELM; i < NLAYERS; i++) {

				// CATOTELM

				Dz[i] = Dz_cato;
				por[i] = cato_por; // as in RW code
				Frac_org[i] = 0.0;
				Frac_min[i] = 0.0;
				Frac_peat[i] = 1.0 - (cato_por + Fgas);

				Fpwp_ref[i] = peat_wp;  // Wisser et al. 2011 Very low in peatlands
				Frac_water_belowpwp[i] = peat_wp;
				Frac_water[i] = cato_por-peat_wp; // i.e. initially saturated

				soiltype.awc_peat[i-IDX] = (por[i] - peat_wp) * Dz[i];
			}

			if (mixedl > 0) { 

				Dz[MIDX] = Frac_water[MIDX] * maxh + Frac_ice[MIDX] * maxh; // MIN value: 50mm
				Frac_org[MIDX] = 0.0;
				por[MIDX] = 0.0; // As it's filled with litter water or ice.
				Frac_min[MIDX] = 0.0;
				Frac_peat[MIDX] = 0.0;
				Frac_water[MIDX] = 0.0;

				// So, Frac_air and Frac_water are NOT defined here
			}
		} // firstTempCalc
		else {

			if (mixedl > 0) {

				// Mixed layer when soilcode == organic
				Dz[MIDX] = Frac_water[MIDX] * maxh + Frac_ice[MIDX] * maxh;
				Frac_org[MIDX] = 0.0;
			}
		} // not firstTempCalc
	}
	else {

		// NON-WETLAND SOILS

		if (firstTempCalc && patch.stand.first_year == date.year) {

			for (int ii = IDX; ii<NLAYERS; ii++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 7) {

					// No updates to the standard values for LPJ-GUESS soils
					por[ii] = soiltype.porosity;
					Frac_org[ii] = soiltype.organic_frac;
					Frac_min[ii] = soiltype.mineral_frac;
					Fpwp_ref[ii] = soiltype.water_below_wp;
					Frac_water_belowpwp[ii] = soiltype.water_below_wp;
				} 
				else {			
					// Use values updated using input of soil carbon
					por[ii] = soiltype.porosity_gridcell[ii-IDX];
					Frac_org[ii] = soiltype.org_frac_gridcell[ii-IDX];
					Frac_min[ii] = soiltype.min_frac_gridcell[ii-IDX];

					double layerdepth_upper = SOILDEPTH_UPPER / NSOILLAYER_UPPER;
					double layerdepth_lower = SOILDEPTH_LOWER / NSOILLAYER_LOWER;

					if (ii < IDX+NSOILLAYER_UPPER) {
						// before organic soils: Fpwp_ref[ii] = soiltype.water_below_wp;
						Fpwp_ref[ii] = soiltype.wp[ii-IDX] / layerdepth_upper; 
					}
					else {
						Fpwp_ref[ii] = soiltype.wp[ii-IDX] / layerdepth_lower; 
					}

					Frac_water_belowpwp[ii] = Fpwp_ref[ii]; // Water content initially set to the wilting point
				}

				Dz[ii] = Dz_soil;
				Frac_water[ii] = 0.0;
				Frac_peat[ii] = 0.0;
			}
		} 
		else if (firstTempCalc && (restart || patch.stand.clone_year == date.year)) {
			// initialize after restarts
			for (int ii = IDX; ii<NLAYERS; ii++) {

				if (!iforganicsoilproperties || soiltype.soilcode == 7) {

					// No updates to the standard values for LPJ-GUESS soils
					por[ii] = soiltype.porosity;
					Frac_org[ii] = soiltype.organic_frac;
					Frac_min[ii] = soiltype.mineral_frac;
					Fpwp_ref[ii] = soiltype.water_below_wp;
				}
				else {
					// Use values updated using input of soil carbon
					por[ii] = soiltype.porosity_gridcell[ii - IDX];
					Frac_org[ii] = soiltype.org_frac_gridcell[ii - IDX];
					Frac_min[ii] = soiltype.min_frac_gridcell[ii - IDX];

					double layerdepth_upper = SOILDEPTH_UPPER / NSOILLAYER_UPPER;
					double layerdepth_lower = SOILDEPTH_LOWER / NSOILLAYER_LOWER;

					if (ii < IDX + NSOILLAYER_UPPER) {
						// before organic soils: Fpwp_ref[ii] = soiltype.water_below_wp;
						Fpwp_ref[ii] = soiltype.wp[ii - IDX] / layerdepth_upper;
					}
					else {
						Fpwp_ref[ii] = soiltype.wp[ii - IDX] / layerdepth_lower;
					}
				}
				Frac_peat[ii] = 0.0;
			}
		}// firstTempCalc
	}
 }*/

/*void Soil::snowpack_dynamics(const double &snowdepth, const int& soilsurfaceindex, int& snow_active_layers) {

	// DESCRIPTION
	// Prvate function called daily from main soil temperature routine.
	// Update the depth of each snow layer, and updates the number of active snow layers.

	// INPUT:
	// snowdepth		- the daily snow depth (mm)
	// soilsurfaceindex	- index in the Dz array

	double toplayerdepth = 50.0; // mm (as in JSBACH)
	double fulldepth = NLAYERS_SNOW * toplayerdepth;
	
	// Simple snow layer scheme. The snow layer nearest the soil surface has variable depth.
	// We add/remove toplayerdepth (e.g. 50mm) layers as the snow deepens/melts. 
	if (snowdepth > fulldepth) {
		snow_active_layers = NLAYERS_SNOW;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 3] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 4] = snowdepth-(snow_active_layers -1)*toplayerdepth;
	}
	else if (snowdepth <= fulldepth && snowdepth > (NLAYERS_SNOW - 1)*toplayerdepth) { // 200-250
		snow_active_layers = NLAYERS_SNOW - 1;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 3] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else if (snowdepth <= (NLAYERS_SNOW - 1)*toplayerdepth && snowdepth > (NLAYERS_SNOW - 2)*toplayerdepth) { // 150-200
		snow_active_layers = NLAYERS_SNOW - 2;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 2] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else if (snowdepth <= (NLAYERS_SNOW - 2)*toplayerdepth && snowdepth > (NLAYERS_SNOW - 3)*toplayerdepth) { // 100-150
		snow_active_layers = NLAYERS_SNOW - 3;
		Dz[soilsurfaceindex - snow_active_layers] = toplayerdepth;
		Dz[soilsurfaceindex - snow_active_layers + 1] = snowdepth - (snow_active_layers - 1)*toplayerdepth;
	}
	else {
		snow_active_layers = NLAYERS_SNOW - 4; // 0-100 - usually a single layer
		Dz[soilsurfaceindex - snow_active_layers] = snowdepth;
	}

}*/


/*bool Soil::soil_temp_multilayer(const double &dailyairtemp) {

	// DESCRIPTION
	// Main soil temperature routine. Called daily from dailyaccounting_patch
	// To update the temperature of each snow, soil and padding layer

	// INPUT:
	// dailyairtemp - the daily air temperature

	int layer0, layer;

	// Layer numbers
	double surf_temp;   // surface temperature [deg. C]
	double Dsnow;		// diffusivity of snow [mm2 d-1]
	double Csnow;		// heat capacity of snow [J mm-3 K-1]
	double Ksnow;		// thermal conductivity of snow [J d-1 mm-1 K-1]

	bool analyticalSolutionTest = false;

	const double ERROR_TEMP = 80.0; // degrees C


	// INITIALISE
	// ------------------------------------------

	int daynum = date.day;

	if (firstTempCalc && patch.stand.first_year == date.year) {

		for (int i = 0; i<NLAYERS; i++) {
			Frac_ice[i] = 0.0;
			// Set the soil temperature to 0 degC initally.
			// Could also initialise T_soil with Tm (MAAT) * depth(m) * 0.015 as in Wisser et al. (2011)
			T_soil[i] = 0.0;
			T_soil_yesterday[i] = 0.0;
			T_old[i] = 0.0;
			Di[i] = 0.04 * MM2_PER_M2;
			Ci[i] = Cp_air * M3_PER_MM3;
			Ki[i] = Di[i] * Ci[i];
			Dz[i] = 100.0;
			T[i] = 0.0;
		}

		// Set up padding - use calc_padding to calculate the thickness
		// of the padding layers, initially using 100mm as the lower layer, 
		// and initialise the padding layer temperatures to zero.

		calc_padding(Dz_cato, pad_dz);

		for (layer = 0; layer<PAD_LAYERS; layer++) {
			pad_temp[layer] = 0.0;
		}

		// The log is used to speed things up. Without the log you need
		// to use an exponential, with it you just multiply.
		lKorg	= log(Korg);
		lKpeat	= log(Kpeat);
		lKmin	= log(Kmin);
		lKwater	= log(Kwater);
		lKice	= log(Kice);
		lKair	= log(Kair);

		if (patch.stand.is_highlatitude_peatland_stand())
			ngroundl = NACROTELM + NCATOTELM; // Usually 15
		else
			ngroundl = NSOILLAYER;

	} // firstTempCalc
	else if (firstTempCalc && (restart || patch.stand.clone_year == date.year)) {

		// The log is used to speed things up. Without the log you need
		// to use an exponential, with it you just multiply.
		lKorg = log(Korg);
		lKpeat = log(Kpeat);
		lKmin = log(Kmin);
		lKwater = log(Kwater);
		lKice = log(Kice);
		lKair = log(Kair);

		if (patch.stand.is_highlatitude_peatland_stand())
			ngroundl = NACROTELM + NCATOTELM; // Usually 15
		else
			ngroundl = NSOILLAYER;
	}


	// Error check
	if (!valid_layer_num(ngroundl)) {
		dprintf("Soil::soil_temp_multilayer - INVALID LAYER ERROR!!!\n");
		return false;
	}

	// Reset ice_frac
	for (int i = 0; i<NLAYERS; i++) {
		Di[i] = 0.04 * MM2_PER_M2;
		Ci[i] = Cp_air * M3_PER_MM3;
		Ki[i] = Di[i] * Ci[i];
	}

	// Allocate yesterday's temperature and ice fraction to today's 
	// ------------------------------------------

	if (!firstTempCalc || restart || patch.stand.clone_year == date.year) update_from_yesterday();

	// SNOW DENSITY and THERMAL PROPRTIES
	// ------------------------------------------

	bool snow_active_old = snow_active;
	update_snow_properties(daynum, dailyairtemp, Dsnow, Csnow, Ksnow);

	if (analyticalSolutionTest) snow_active = 0;

	// MIXED LAYER
	// On wetlands there can be a layer that is filled with water, depending on the WTP.
	// Set to 1 if this is activated
	int mixedl = 0;

	// Set indexes for first active layers:
	IDX = NLAYERS - ngroundl;							// usually 7
	MIDX = NLAYERS - ngroundl - mixedl;					// usually 7 or 6

	// SNOW LAYERS
	// ------------------------------------------

	// Determine how many snow layers are active, based on snow_active, snow depth and snow_active_layers (<=NLAYERS_SNOW);
	int snow_active_layers_old = snow_active_layers;

	double snowdepth = snowpack / (snowdens / water_density); // mm
	int soilsurfaceindex = NLAYERS - ngroundl - mixedl;

	if (ifmultilayersnow && !iftwolayersoil) {

		// Allow multiple layers?
		if (snow_active) {
			snowpack_dynamics(snowdepth, soilsurfaceindex, snow_active_layers);
		}
		else {
			snow_active_layers = 0;
		}
	}
	else {

		// ONE snow layer
		if (snow_active) {

			snow_active_layers = 1; // a single layer
			Dz[soilsurfaceindex - snow_active_layers] = snowdepth;
		}
		else {
			snow_active_layers = 0;
		}
	}

	// SNOW INDEX
	SIDX = NLAYERS - ngroundl - mixedl - snow_active_layers;

	// Update the thermal properties in the snow layers
	for (int j = SIDX; j< SIDX + snow_active_layers; j++) {
		Di[j] = Dsnow;
		// correct units
		Ci[j] = Csnow;
		Ki[j] = Ksnow;
	}

	// Set fractions of soil parameters for each soil layer & litter
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) {
		init_hydrology_variables(); // Initialise whc[], alwhc[], and alwhc_init[] for this patch
	}

	update_layer_fractions(daynum, mixedl, MIDX);

	// Update surface (air) temperature 
	surf_temp = dailyairtemp;

	// Update temperature in SNOW layers 
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) {
		SIDX_old = SIDX;
		snow_days_prev = 365;
		snow_days = 0;
		snow_active = false;
		snow_active_layers = 0;
	}

	// Counting snow days for the Wania et al. compaction scheme
	if (snow_active) {

		// If there's still snow...
		snow_days++;

		if (!snow_active_old) {
			// Reset
			snow_days = 1;
		}
	}
	else {
		if (snow_active_old)
			snow_days_prev = snow_days; // The day after all the snow has melted

		snow_days = 0;
	}

	SIDX_old = SIDX;

	// Now update snow layer temperatures if snow is active && the number of snow layers has changed
	if (snow_active && snow_active_layers != snow_active_layers_old) {

		if (snow_active_layers > snow_active_layers_old) {
			// new snow layer(s) - set new layer temperatures to air temp
			for (int j = SIDX; j<SIDX + snow_active_layers - snow_active_layers_old; j++) {
				T_soil[j] = surf_temp;
			}
			// The existing snow layers keep their temperatures
		}
		else if (snow_active_layers < snow_active_layers_old) {
			T_soil[SIDX - 1] = surf_temp;
		}
	}

	// AIR layer
	// ------------------------------------------

	Dz[SIDX - 1] = AIR_THICKNESS; // mm
	Frac_air[SIDX - 1] = 1.0;
	por[SIDX - 1] = 1.0;
	Frac_org[SIDX - 1] = 0.0;
	Frac_min[SIDX - 1] = 0.0;
	Frac_peat[SIDX - 1] = 0.0;
	Di[SIDX - 1] = Kair / Cp_air * SECS_PER_DAY * MM2_PER_M2;	// mm2 day-1
	Ci[SIDX - 1] = Cp_air * M3_PER_MM3;					// unit conversion J mm-3 K-1
	Ki[SIDX - 1] = Di[SIDX - 1] * Ci[SIDX - 1];			// units: (mm2 day-1) * (J mm-3 K-1) = J day-1 mm-1 K-1 , as required

	// Surface air temperature
	T_soil[SIDX - 1] = surf_temp;

	// NB - layer0 defined here
	layer0 = SIDX - 1;

	// Initialise unused T_soil values.
	for (int i = 0; i<SIDX - 1; i++) T_soil[i] = MISSING_VALUE;

	// Update the heat capacities, thermal conductivities and diffusivities (Di[i]) in the soil layers
	if (!iftwolayersoil)
		update_soil_diffusivities(daynum, analyticalSolutionTest);

	// Do one timestep of the Crank-Nicholson method.
	// ------------------------------------------

	if (firstTempCalc && patch.stand.first_year == date.year) calc_padding(Dz[NLAYERS - 1], pad_dz);

	// Initialise T array for the layers
	for (int i = SIDX - 1; i<NLAYERS; i++) T[i] = T_soil[i];

	if (!iftwolayersoil) {

		// Calculate the temperature - added the timesteps loop, just in case
		for (int timestep = 0; timestep < TIMESTEPS; timestep++) {

			double dayfrac = Dt / (double)TIMESTEPS;

			// Wania at al. (2009a) algorithm. Assumes vertical homogeneity in Ki and Ci.
			cnstep(layer0, Di, Dz, surf_temp, dayfrac, pad_dz, T, pad_temp);
			// Alternative algorithm that does not assume vertical homogeneity in Ki and Ci:
			// cnstep_full(layer0, Di, Dz, surf_temp, dayfrac, pad_dz, T, pad_temp, Ki, Ci);

		}
	}

	// ERROR Checking!
	for (int i = SIDX; i<NLAYERS; i++) {

		if (Di[i] <= 0.0)
			dprintf("%g%s%g\n", (double)i, ", BAD Di[i]: ", Di[i]);

		if (T[i] < -ERROR_TEMP || T[i] > ERROR_TEMP) { // -70 exceeded in Siberia sometimes

			dprintf("Soil::soil_temp_multilayer - SOIL TEMP. ERROR!!!\n");
			dprintf("%g%s%g\n", (double)date.year, "  ", (double)date.day);
			dprintf("%g%s%g%s%g\n", (double)i, ":  T[i]: ", T[i], ", Di[i]: ", Di[i]);

			return false;
		} 
		else 
			T_soil[i] = T[i]; // Store temperature values
	}

	// PHASE CHANGES
	// ------------------------------------------

	// - Global Boolean to regulate phase changes
	// - Added 100-year (approx) condition to get the water balance right during spin-up
	if (!iftwolayersoil && ifallowphasechanges && date.year >= FIRST_FREEZE_YEAR) {
		update_ice_fraction(daynum, MIDX);
	} // allowPhaseChanges?

	// Update whc[], awhc[], and wcont[]
	if (!update_layer_water_content(daynum))
		return false;

	if (analyticalSolutionTest) {
		for (int i = SIDX - 1; i < NLAYERS; i++) Frac_water[i] = 0.0; // No water if analytic soln. sought;
	}

	// THAW DEPTH
	// ------------------------------------------

	if (daynum == 0)
		maxthawdepththisyear = 0.0;

	// Define the thaw depth as the 0degC isotherm (Muller 1947, Burn 1998), 
	// which equals the depth, until which liquid water exists continuously from the top down.

	double thaw_depth_today;	// daily thawing depth, where at least some water is found, but ice can be present [mm]
	double thaw_top;	// thawing depth, where all the ice has melted [mm]

	thaw_depth_today = 0.0;

	if (!iftwolayersoil) {

		// Include padding layers
		int i = IDX;
		while (i < NLAYERS && T_soil[i] > 0.0) {
			thaw_depth_today += Dz[i];
			i++;
		}
	} 
	else {
		thaw_depth_today = SOILDEPTH_UPPER + SOILDEPTH_LOWER; // no freezing
	}

	// Record the daily thawdepth. NB! Needed for fire. 
	thaw = thaw_depth_today;

	if (thaw > maxthawdepththisyear)
		maxthawdepththisyear = thaw;

	// Calculate also at which depth there is no ice anymore. This is
	// needed in the waterbalance SR.

	thaw_top = 0.0;
	int j = IDX;
	while (j < NLAYERS && Frac_ice[j]<0.001) {
		thaw_top += Dz[j];
		j++;
	}

	if (!iftwolayersoil) {

		// Record Frac_ice and T_soil EVERY day now, not just Dec 31
		for (int i = 0; i<NLAYERS; i++) {
			Frac_ice_yesterday[i] = Frac_ice[i];
			T_soil_yesterday[i] = T_soil[i];
		}

		// Finally, update soil.temp25 with the calculated value:
		int soil25_ix = IDX + 2;
		if (Dz_soil == 100.0) {
			temp25 = T_soil[soil25_ix];
		}

		// Update monthly, later averages
		int l_ix = IDX;
		for (int sl = 0; sl < SOILTEMPOUT; sl++) {
			T_soil_monthly[date.month][sl] += T_soil[l_ix + sl] / (double)date.ndaymonth[date.month];
		}
	}

	// Not the first time this method is called	
	firstTempCalc = false;

	return true;

} // Soil::calcsoiltemp*/

/*double Soil::nmass_avail(int pref) {
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
	} else if (pref == NH4) {
		nmass = NH4_mass;
	} else if (pref == NO3) {
		nmass = NO3_mass;
	}
	return nmass;
}*/

/*void Soil::nmass_subtract(double nmass, int pref) {

	if (!ifntransform) {
		pref = NH4;
	}
	double residual = 0.0;
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		if (nmass_tot > 0.0) {
			residual = nmass_tot - nmass;
			if (residual >= 0.0 || negligible(residual,10)) {
				nmass_subtract(nmass * NH4_mass / nmass_tot,NH4);
				nmass_subtract(nmass * NO3_mass / nmass_tot,NO3);
			} else {
				fail("tried to subtract more N than available");
			}
		}
	} else if (pref == NH4) {
		if (nmass<NH4_mass || negligible(nmass-NH4_mass,10)) {
			NH4_mass -= nmass;
		} else {
			double r = nmass - NH4_mass;
			NH4_mass = 0.0;
			patch.fluxes.report_flux(Fluxes::NH3_SOIL, r);
			dprintf("NH4 mass %f N substr %f \n",NH4_mass,nmass);
			//fail("tried to subtract more NH4 than available");
		}
	} else if (pref == NO3) {
		if (nmass<NO3_mass || negligible(nmass-NO3_mass,10)) {
			NO3_mass -= nmass;
		} else {
			double r = nmass - NO3_mass;
			NO3_mass = 0.0;
			patch.fluxes.report_flux(Fluxes::NO_SOIL, r);
			dprintf("NO3 mass %f N substr %f \n",NO3_mass,nmass);
			//fail("tried to subtract more NO3 than available");
		}
	}

}*/

/*void Soil::nmass_inc(double nmass, int pref) {
	if (nmass < 0.0) {
		nmass_subtract(fabs(nmass), pref);
	} else {
		if (!ifntransform) {
			pref = NH4;
		}
		if (pref == NO) {
			double nmass_tot = NH4_mass + NO3_mass;
			NH4_mass+=nmass * NH4_mass / nmass_tot;
			NO3_mass+=nmass * NO3_mass / nmass_tot;
		} else if (pref == NH4) {
			NH4_mass += nmass;
		} else if (pref == NO3) {
			NO3_mass += nmass;
		}
	}
}*/

/*void Soil::nmass_multiplic_inc(double inc, int pref) {
	if (!ifntransform) {
		pref = NH4;
	}
	if (pref == NO) {
		double nmass_tot = NH4_mass + NO3_mass;
		NH4_mass *= inc * NH4_mass / nmass_tot;
		NO3_mass *=inc * NO3_mass / nmass_tot;
	} else if (pref == NH4) {
		NH4_mass *= inc;
	} else if (pref == NO3) {
		NO3_mass *= inc;
	}
}*/

// serialize new soil variables, when finalised
/*void Soil::serialize(ArchiveStream& arch) {
	arch & wcont
		& wcont_evap
		& awcont_upper
		& dwcontupper
		& dwcontlower
		& mwcontupper
		& mwcontlower
		& mwcont
		& snowpack
		& snow_active
		& snow_days
		& snow_days_prev
		& snow_active_layers
		& snow_water
		& snow_ice
		& msnowdepth
		& dsnowdepth
		& thaw
		& runoff
		& temp25
		& Dz // optimise through init after restart?
		& T_old
		& T_soil_yesterday
		& T_soil_monthly
		& dtemp
		& mtemp
		& gtemp
		& pad_temp
		& pad_dz
		& cpool_slow
		& cpool_fast
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean
		& alag
		& exp_alag
		& Frac_ice_yesterday
		& Frac_water_belowpwp
		& Frac_water
		& Frac_air
		& whc
		& whc_evap
		& alwhc
		& aw_max // optimise through init after restart?
		// Methane parameters
		// Some of these could possibly be removed to 
		// optimise for memory in state files
		& wtp
		& ch4_store
		& co2_store
		& CO2_soil_yesterday
		& CH4_yesterday
		& CH4_diss_yesterday
		& CH4_gas_yesterday
		& CH4_gas_vol
		& O2;
		
	for (int i = 0; i<NSOMPOOL; i++) {
		arch & sompool[i];
	}

	arch & dperc
		& orgleachfrac
		& NO2_mass
		& NO_mass
		& N2O_mass
		& N2_mass
		& NH4_mass
		& NO3_mass
		& NH4_input
		& NO3_input
		& anmin
		& animmob
		& aminleach
		& aorgNleach
		& aorgCleach
		& anfix
		& anfix_calc
		& anfix_mean
		& snowpack_NH4_mass
		& snowpack_NO3_mass
		& solvesomcent_beginyr
		& solvesomcent_endyr
		& solvesom
		& fnuptake_mean
		& morgleach_mean
		& mminleach_mean;
}*/


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
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
// For Arctic and wetland references see soil.h
