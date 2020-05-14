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

void Snowlayer::update_heatcap_and_cond() {

	double max_temp = 50.0;
	double min_temp = -50.0;
	double snow_temp = min(max_temp, max(min_temp, temp));

	// Snow density including water withing snow
	double tot_dens = min(rho_ICE, (imass + wmass) / depth * MM_PER_M);

	// C Snow - Fukusako, Eq. 2. [J m-3 K-1]
	Ch = (0.185 + 0.689 * (K2degC + snow_temp) * 0.01) * tot_dens * J_PER_KJ;

	// K Snow - Best et al. 2011, Eq. 23 [J m-1 K-1 s-1]
	Kh = 2.22 * pow((tot_dens / rho_H2O), 1.88);

	// D Snow
	Dh = Kh / Ch;
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
void Snow::potential_ice(double current_imass, std::vector<double>& pot_imass) {
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		if (current_imass > snow_depth_standard * snow_depth_limit) {
			pot_imass[lyr] = snow_depth_standard;
			current_imass -= snow_depth_standard;
		}
		else if (current_imass <= snow_depth_standard * snow_depth_limit) {
			pot_imass[lyr] = current_imass;
			current_imass = 0.0;
		}
	}
}

/// Re-distributing snow after freeze/melt event
void Snow::melt_dens(double current_imass) {
	double initial_wmass = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_wmass += snowlayers[lyr].imass + snowlayers[lyr].wmass;
	}

	// Potential ice
	std::vector<double> pot_imass(nsnowlayers);

	// Check potential ice content
	potential_ice(current_imass, pot_imass);

	// Re-calc ice & dens of layers
	for (int lyr = 0; lyr < nsnowlayers - 1; lyr++) {
		int lyr_above = lyr;

		// Move snow downwards in snowpack to fill every layer to pot_ice
		while (snowlayers[lyr].imass != pot_imass[lyr] && lyr_above < nsnowlayers - 1) {
			lyr_above++;

			if (!negligible(pot_imass[lyr])) {

				// Either bottom layer shrinking from being larger than standard_depth
				// or water freezing making layer exceed potential
				if (snowlayers[lyr].imass > pot_imass[lyr]) {
					double imass_move = snowlayers[lyr].imass - pot_imass[lyr];

					// Calculate new snow density
					double imass_sum = imass_move + snowlayers[lyr_above].imass;
					double w_new = imass_move / imass_sum;
					double w_old = 1.0 - w_new;
					snowlayers[lyr_above].dens = min(rho_ICE, w_old * snowlayers[lyr_above].dens + w_new * snowlayers[lyr].dens);

					// Add and remove water
					double wmass_move = imass_move / snowlayers[lyr].imass * snowlayers[lyr].wmass;
					snowlayers[lyr_above].wmass += wmass_move;
					snowlayers[lyr].wmass -= wmass_move;

					// Add and remove ice
					snowlayers[lyr_above].imass += imass_move;
					snowlayers[lyr].imass -= imass_move;

					// If new layer created
					if (snowlayers[lyr_above].Dh == 0.0) {
						snowlayers[lyr_above].Kh = snowlayers[lyr].Kh;
						snowlayers[lyr_above].Ch = snowlayers[lyr].Ch;
						snowlayers[lyr_above].Dh = snowlayers[lyr].Dh;
						snowlayers[lyr_above].dens = snowlayers[lyr].dens;
					}
				}
				else {

					double imass_move = max(0.0, min(snowlayers[lyr_above].imass, pot_imass[lyr] - snowlayers[lyr].imass));

					// Add ice to be moved if it isn't negligible
					if (!negligible(imass_move) && imass_move > 0.0) {

						// Calculate new snow density
						double imass_sum = imass_move + snowlayers[lyr].imass;
						double w_new = imass_move / imass_sum;
						double w_old = 1.0 - w_new;
						snowlayers[lyr].dens = min(rho_ICE, w_old * snowlayers[lyr].dens + w_new * snowlayers[lyr_above].dens);

						// Add and remove water
						double wmass_move = imass_move / snowlayers[lyr_above].imass * snowlayers[lyr_above].wmass;
						snowlayers[lyr].wmass += wmass_move;
						snowlayers[lyr_above].wmass -= wmass_move;

						// Add and remove ice

						snowlayers[lyr].imass += imass_move;
						snowlayers[lyr_above].imass -= imass_move;
					}
				}
			}
			else {
				// Erase empty snow layer
				snowlayers.erase(snowlayers.begin() + lyr);
			}
		}
	}

	double final_wmass = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_wmass += snowlayers[lyr].imass + snowlayers[lyr].wmass;
	}

	const double EPS = 1.0e-8;
	double diff = abs(initial_wmass - final_wmass);


	//dprintf(" balance %d %d  %g %g %g \n",date.year, date.day, initial_wmass, final_wmass, diff);
	assert(abs(initial_wmass - final_wmass) < EPS);
}

/// Re-distributing snow after freeze/melt event, when snowpack > max potential snowdepth
void Snow::melt_dens_full() {

	// filling each layer to standard_depth
	for (int lyr = NLAYERS_SNOW_MAX - 1; lyr > 0; lyr--) {
		int lyr_below = lyr;

		// Move snow upwards in snowpack to fill every layer to standard_depth
		while (snowlayers[lyr].imass != snow_depth_standard) {
			lyr_below--;

			// Move upwards maximum of what is in snow layer lyr_below
			double imass_move = min(snowlayers[lyr_below].imass, snow_depth_standard - snowlayers[lyr].imass);

			// New snow depth with added snow
			double new_depth = min(snow_depth_standard, snowlayers[lyr].imass + imass_move);

			// Calculate new snow density
			double w_new = imass_move / new_depth;
			double w_old = 1.0 - w_new;
			snowlayers[lyr].dens = min(rho_ICE, w_old * snowlayers[lyr].dens + w_new * snowlayers[lyr].dens);

			// Add and remove water
			double wmass_move = imass_move / snowlayers[lyr_below].imass * snowlayers[lyr_below].wmass;
			snowlayers[lyr].wmass += wmass_move;
			snowlayers[lyr_below].wmass -= wmass_move;

			// Add and remove ice
			snowlayers[lyr].imass += imass_move;
			snowlayers[lyr_below].imass -= imass_move;
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
		snowlayers[lyr].whc = snowlayers[lyr].imass * (rw_min + (rw_max - rw_min) * max(0.0, dens_t - snowlayers[lyr].dens) / dens_t);
	}
}

/// Balance ice and water in snow layers after melting and refreezing
void Snow::ice_water_snow_balance() {

	// No snow layers or not multilayer snow scheme then leave function
	if (!nsnowlayers || !ifmultilayersnow)
		return;

	double initial_snow = 0.0;
	double sum_imass = 0.0;
	initial_snow = melt;

	// Sum up ice in snowpack
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		sum_imass += snowlayers[lyr].imass;
		initial_snow += snowlayers[lyr].imass + snowlayers[lyr].wmass;
	}

	// Re-distribute layers (only based on ICE!)
	double max_snow = NLAYERS_SNOW_MAX * snow_depth_standard;

	// Balance layers depending on if snowpack is over maxumum or not
	if (sum_imass != 0) {
		if (sum_imass > max_snow) {
			melt_dens_full();
		}
		else {
			melt_dens(sum_imass);
		}
	}

	// Check for insignificant top layers	
	for (int lyr = nsnowlayers - 1; lyr >= 0; lyr--) {

		// If layer is too thin then add it to layer below
		if (snowlayers[lyr].imass < snow_depth_standard * (snow_depth_limit - 1.0)) {

			// Check so this layer isn't the bottom layer
			if (lyr != 0) {

				// Calculate old and new snow fractions
				double imass_sum = snowlayers[lyr - 1].imass + snowlayers[lyr].imass;
				double imass_new = !negligible(imass_sum) ? snowlayers[lyr].imass / imass_sum : 0.0;
				double imass_old = 1.0 - imass_new;

				// Update ice, water and snow density in layer below
				snowlayers[lyr - 1].imass += snowlayers[lyr].imass;
				snowlayers[lyr - 1].wmass += snowlayers[lyr].wmass;
				snowlayers[lyr - 1].dens = imass_new * snowlayers[lyr].dens + imass_old * snowlayers[lyr - 1].dens;

				// Erase this empty snow layer
				snowlayers.erase(snowlayers.begin() + lyr);
				nsnowlayers = (int)snowlayers.size();
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
		if (snowlayers[lyr].imass == 0.0) {
			snowlayers[lyr].whc = 0.0;
			ts_melt += snowlayers[lyr].wmass;
			snowlayers[lyr].wmass = 0.0;

			// Erase this empty snow layer
			snowlayers.erase(snowlayers.begin() + lyr);
		}
		else {
			// Fill snow layer with all melt water from above
			snowlayers[lyr].wmass += ts_melt;
			ts_melt = 0.0;

			// Pass water above holding capacity downwards
			if (snowlayers[lyr].wmass > snowlayers[lyr].whc) {
				ts_melt = snowlayers[lyr].wmass - snowlayers[lyr].whc;
				snowlayers[lyr].wmass = snowlayers[lyr].whc;
			}
		}
	}

	// Update number of snow layers
	nsnowlayers = (int)snowlayers.size();

	// Add excess water this timesteps melt to surface melt-water to pass to soil on next timestep 
	melt += ts_melt;

	double final_snow = melt;

	// Update snow depth and calculate total snow
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);
		final_snow += snowlayers[lyr].imass + snowlayers[lyr].wmass;
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

	if (!ifmultilayersnow)
		return;

	// Minimum air volume that will influence snow layer (only affect top snow layer)
	double min_snow_depth = snow_depth_standard * (snow_depth_limit - 1.0) * 10.0;

	// Going through snow layers
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// Snow layer depth
		snow_depth = lyr == nsnowlayers - 1 ? max(min_snow_depth, snowlayers[lyr].depth) : snowlayers[lyr].depth;

		// Water and ice fraction of full volume
		double Frac_wmass = snowlayers[lyr].wmass / snow_depth;
		double Frac_imass = snowlayers[lyr].imass / snow_depth;

		double init_wmass = snowlayers[lyr].wmass;
		double init_imass = snowlayers[lyr].imass;
		double init_temp = snowlayers[lyr].temp;

		// Volume of snow layer [m-3]
		snow_volume = snowlayers[lyr].depth * M_PER_MM;

		// The total amount of ice and water should not change in this routine
		double wtotal_initial = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		// Temperature change from yesterday
		delta_T = snowlayers[lyr].temp - snowlayers[lyr].temp_old;

		// Freezing
		if (snowlayers[lyr].temp < 0.0 && snowlayers[lyr].temp < snowlayers[lyr].temp_old) { // Wania conditions 
			Fthaw = 0.0;

			// Water to freeze
			if (Frac_wmass > 0.0) {

				// Available energy to freeze water[J m-3]
				energy = Cp_water * fabs(delta_T);

				// Fraction of a volume of water that could freeze [-]
				Ffreez = energy / Lheat;

				// More water in a volume than could freeze
				if (Frac_wmass > Ffreez) {

					// There is some water left, even after freeze
					Frac_imass += Ffreez;
					Frac_wmass -= Ffreez;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
				}
				else {

					// It is cold enough to freeze all the water and
					// decrease the temperature below 0, so only the 
					// energy need for phase change is removed from 
					// the already calculated temperature 

					// All water freezes
					Frac_imass += Frac_wmass;

					// Rescale Fwmass if snow layer depth is lower than min_snow_depth for energy calculation
					Frac_wmass *= snow_depth / snowlayers[lyr].depth;

					// Calculate new snow temperature
					snowlayers[lyr].temp += Frac_wmass * snow_volume * Lheat / snowlayers[lyr].Ch;

					// No liquid water left in snow layer
					Frac_wmass = 0.0;
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
			if (Frac_imass > 0.0) {

				// Energy to thaw ice [J m-3]
				energy = Cp_water * fabs(delta_T);

				// Fraction of a volume of ice that would thaw [-]
				Fthaw = energy / Lheat;

				// More ice in a volume than could thaw
				if (Frac_imass > Fthaw) {

					// There is some ice left, even after melt
					Frac_wmass += Fthaw;
					Frac_imass -= Fthaw;

					// Energy from daily temperature change is used for phase change 
					snowlayers[lyr].temp = 0.0;
				}
				else {
					// There is more energy available than required for
					// melting the rest of the ice. Therefore no need to 
					// calculate new temperature as snow layer will be removed

					// All ice thaws
					Frac_wmass += Frac_imass;

					// No ice left in snow layer
					Frac_imass = 0.0;
				}
			}
			else {
				Fthaw = 0.0;
			}
		} // end of thawing

		  // New ice fraction in snow layer
		snowlayers[lyr].imass = Frac_imass * snow_depth;

		// New water fraction in snow layer
		snowlayers[lyr].wmass = Frac_wmass * snow_depth;

		// The total amount of ice and water should not change in this routine
		double wtotal_final = snowlayers[lyr].imass + snowlayers[lyr].wmass;

		//dprintf("%d %d %d AirT %g Temp %g %g Ice %g %g Water %g %g Tot %g %g\n", date.year, date.day, lyr, airT, snowlayers[lyr].temp, init_temp, snowlayers[lyr].imass, init_imass, snowlayers[lyr].wmass, init_wmass, wtotal_initial, wtotal_final);
		double diff = abs(wtotal_initial - wtotal_final);

		const double EPS = 1.0e-5;

		assert(abs(wtotal_initial - wtotal_final) < EPS);

		// Update T_old after possible phase change
		snowlayers[lyr].temp_old = snowlayers[lyr].temp;

	} // end of loop through layers
}


/// Update snow density
void Snow::update_density() {

	if (!ifmultilayersnow)
		return;

	// Calculate snow depth ice content before compaction
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		// snow depth based only on ice
		snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);

		// snow density, only defined by ice fraction, Oleason 2013, Eq. 6.88. can not exceed maximum density
		snowlayers[lyr].dens = min(rho_ICE, snowlayers[lyr].imass / (snowlayers[lyr].depth * M_PER_MM));
	}

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
		snowlayers[lyr].depth = snowlayers[lyr].imass / (snowlayers[lyr].dens / rho_H2O);
	}

	// Update snow water holding capacity
	update_whc();
}

/// Accumulate snow in bottom layer if all other layers are already full
void Snow::update_layers(double prec, double snowfall_dens) {

	// Check initial snow and water content
	double initial_snow = prec * snowfall_dens;
	double initial_wmass = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		initial_snow += snowlayers[lyr].imass * snowlayers[lyr].dens;
		initial_wmass += snowlayers[lyr].wmass;
	}

	for (int lyr = 0; lyr < nsnowlayers; lyr++) {

		double imass_new, imass_old, imass_sum, imass_frac_move, wmass_move;

		// Total snow to move down to this layer
		double tot_snow_move = lyr == 0 ? prec : snow_depth_standard - snowlayers[lyr].imass;

		// Layer to move snow from
		int lyr_move = lyr + 1;

		// Move snow until enough is moved (could be from more than one layer above)
		while (!negligible(tot_snow_move)) {

			double snow_move;

			// Determine snow to move from snow layer above

			// Snow layer above exist
			if (lyr_move < nsnowlayers) {
				snow_move = min(tot_snow_move, snowlayers[lyr_move].imass);

				// Move water as well
				imass_frac_move = !negligible(snowlayers[lyr_move].imass) ? snow_move / snowlayers[lyr_move].imass : 0.0;
				wmass_move = snowlayers[lyr_move].wmass * imass_frac_move;
				snowlayers[lyr_move].wmass -= wmass_move;
				snowlayers[lyr].wmass += wmass_move;
			}
			else { // No layer above so take from prec to fill layer
				snow_move = min(tot_snow_move, snow_depth_standard - snowlayers[lyr].imass);
			}

			// Calculate old and new snow fractions
			imass_sum = snowlayers[lyr].imass + snow_move;
			imass_old = !negligible(imass_sum) ? snowlayers[lyr].imass / imass_sum : 0.0;
			imass_new = 1.0 - imass_old;

			snowlayers[lyr].imass += snow_move;

			// Update snow density
			if (lyr_move < nsnowlayers) {
				snowlayers[lyr].dens = imass_old * snowlayers[lyr].dens + imass_new * snowlayers[lyr_move].dens;
				snowlayers[lyr_move].imass -= snow_move;
			}
			else {
				snowlayers[lyr].dens = imass_old * snowlayers[lyr].dens + imass_new * snowfall_dens;
			}

			// Remove added snow from this layer from the toal that needs to be moved
			tot_snow_move -= snow_move;

			// Go to the next layer above (if tot_snow_move > 0)
			lyr_move++;
		}
	}

	// Check final snow and water content
	double final_snow = 0.0;
	double final_wmass = 0.0;
	for (int lyr = 0; lyr < nsnowlayers; lyr++) {
		final_snow += snowlayers[lyr].imass * snowlayers[lyr].dens;
		final_wmass += snowlayers[lyr].wmass;
	}

	const double EPS = 1.0e-8;
	assert(abs(initial_snow - final_snow) < EPS);
	assert(abs(initial_wmass - final_wmass) < EPS);
}

/// Snow accumulation
void Snow::snowfall(double airT, double& prec) {

	// If snowfall
	if (airT < 0.0 && prec > 0.0) {

		if (ifmultilayersnow) {

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
				if ((top_snow_lyr < NLAYERS_SNOW_MAX - 1 && (snowlayers[top_snow_lyr].imass + prec) <= snow_depth_standard * snow_depth_limit) ||
					(top_snow_lyr == NLAYERS_SNOW_MAX - 1 && (snowlayers[top_snow_lyr].imass + prec) <= snow_depth_standard)) {

					// update dens
					snowlayers[top_snow_lyr].mix_dens(prec, snowfall_dens);

					// update temperature
					snowlayers[top_snow_lyr].temp = snowlayers[top_snow_lyr].temp * (snowlayers[top_snow_lyr].imass / (snowlayers[top_snow_lyr].imass + prec)) +
						airT * (prec / (snowlayers[top_snow_lyr].imass + prec));

					// update ice
					snowlayers[top_snow_lyr].imass += prec;
					prec = 0.0;
				}
				else {

					// Snowfall will make snow layer exceed maximum snow layer depth

					// Reduce snow layer dept as it is already higher than snow_depth_standard
					// and then add prec to next layer (NLAYERS_SNOW_MAX layer can't get above snow_depth_standard)
					if (snowlayers[top_snow_lyr].imass > snow_depth_standard) {

						double imass_move = snowlayers[top_snow_lyr].imass - snow_depth_standard;
						double frac_move = imass_move / snowlayers[top_snow_lyr].imass;

						snowlayers[top_snow_lyr].imass = snow_depth_standard;

						// Need to move water as well
						double wmass_move = snowlayers[top_snow_lyr].wmass * frac_move;
						snowlayers[top_snow_lyr].wmass -= wmass_move;

						// Create new snow layer with same density as original snow layer
						snowlayers.resize(nsnowlayers + 1);
						snowlayers[nsnowlayers] = Snowlayer(imass_move, wmass_move, snowlayers[top_snow_lyr].dens, snowlayers[top_snow_lyr].temp);
						nsnowlayers++;
					}
					else {
						// Fill up to standard_depth and then put rest in next layer
						double prec_lyr = snow_depth_standard - snowlayers[top_snow_lyr].imass;

						// update dens
						snowlayers[top_snow_lyr].mix_dens(prec_lyr, snowfall_dens);

						// update ice
						snowlayers[top_snow_lyr].imass += prec_lyr;
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

	// Simple single snow layer scheme (v4.0)
	if (!ifmultilayersnow) {

		// Daily calculation of snowfall and rainfall from precipitation and snow melt from
		// snow pack; update of snow pack with new snow and snow melt and snow melt

		// Maximum temperature for precipitation as snow (deg C)
		const double TSNOW = 0.0;
		
		// Maximum size of snowpack (mm) (S. Sitch, pers. comm. 2001-11-28)
		const double SNOWPACK_MAX = 10000.0;

		// Fixed snow density (kg m-3) JSBACH value from Ekici et al (2015)
		const double SNOWDENSITY = 250.0;

		// Snow melt
		double snowmelt = 0.0;;

		// New snow
		double snowfall = 0.0;

		// Snowing today
		if (airT < TSNOW && prec > 0.0) {

			// If no snow layers
			if (nsnowlayers == 0) {
				snowlayers.resize(nsnowlayers + 1);
				snowlayers[nsnowlayers] = Snowlayer(0.0, 0.0, 0.0, airT);
				nsnowlayers++;
			}

			// Snowfall up to max snowpack size
			snowfall = min(prec, SNOWPACK_MAX - snowlayers[0].imass);
		}
		// Raining today
		else if (airT > TSNOW && nsnowlayers) {
			// Snow melt formulation (Choudhury et al 1998)
			snowmelt = min((1.5 + 0.007 * prec) * (airT - TSNOW), snowlayers[0].imass);
		}

		// Update single snowlayer with new snow or melt
		if (nsnowlayers) {
			snowlayers[0].imass += snowfall - snowmelt;
			snowlayers[0].dens = SNOWDENSITY;
			snowlayers[0].wmass = 0.0;
			snowlayers[0].depth = snowlayers[0].imass / (snowlayers[0].dens / rho_H2O);

			if (negligible(snowlayers[0].imass)) {

				// Erase this empty snow layer
				snowlayers.erase(snowlayers.begin());
				nsnowlayers--;
			}
		}

		// Update precipitation and snow melt
		prec -= snowfall;
		melt = snowmelt;
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

Soillayer::Soillayer(double d, bool a, bool c) : depth(d), acrotelm(a), catotelm(c) {

	// Initialises certain member variables

	// Energy
	temp = 0.0;
	Kh = 0.0;
	Ch = 0.0;

	// pH
	pH = 6.5;

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

	// N pools
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
	labile_carbon_w = 0.0;
	labile_carbon_d = 0.0;

	// Monthly variables
	for (int mth = 0; mth < 12; mth++) {
		mwcont[mth] = 0.0;
		fNH4uptake_mean[mth] = 0.0;
		fNO3uptake_mean[mth] = 0.0;
		morgleach_mean[mth] = 0.0;
		morgmove_mean[mth] = 0.0;
		mminleach_mean[mth] = 0.0;
		mminmove_mean[mth] = 0.0;
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
void Soillayer::set_properties(Soiltype soiltype, double socfrac) {

	// Only enter first year except if dynamic soil organic C influence soil properties (iforganicsoilproperties == 2)
	if (date.year && (iforganicsoilproperties != 2 || acrotelm || catotelm)) // DW_COMMENT restart years needs to update it then as well!
		return;

	clay_frac = soiltype.clay_frac;
	sand_frac = soiltype.sand_frac;
	silt_frac = soiltype.silt_frac;

	// Soil pH
	pH = soiltype.pH;

	/////////////////////////
	// soiltype variable for Medium - fine // DW_COMMENT needs to be changed
	alpha = 0.83; // DW_COMMENT should come from input
	Kw_max = 2.63E-07; // DW_COMMENT should come from input
	m = 1.0 - 1.0 / 1.2539; // DW_COMMENT should come from input
	lambda = -0.5884; // DW_COMMENT should come from input
	/////////////////////////

	if (acrotelm) {
		wsats = acrotelm_por * depth;
		awc = acrotelm_por * depth;
		wp = peat_wp * depth;
		por = acrotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - acrotelm_por;
	}
	else if (catotelm) {
		wsats = catotelm_por * depth;
		awc = catotelm_por * depth;
		wp = peat_wp * depth;
		por = catotelm_por * depth;
		organic_frac = 0.0;
		mineral_frac = 0.0;
		peat_frac = 1.0 - catotelm_por;
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

		// Analytic temp calc parameters
		thermdiff_0 = soiltype.thermdiff_0;
		thermdiff_15 = soiltype.thermdiff_15;
		thermdiff_100 = soiltype.thermdiff_100;

		// Two layer hydrology parameters
		perc_base = soiltype.perc_base;
		perc_exp = soiltype.perc_exp;
	}
	// Soil organic content influences soil properties
	else {

		// Soil layer C content (kg C m-2)
		double soil_c = iforganicsoilproperties == 1 ? socfrac * soiltype.soilC : ccont();

		// Maximum soil layer C content (kg C m-2)
		double soil_c_max = maxSOCdensity * (depth / MM_PER_M);

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

		// Analytic temp calc parameters
		thermdiff_0 = min_frac * soiltype.thermdiff_0 + org_frac * soiltype.thermdiff_0_org;
		thermdiff_15 = min_frac * soiltype.thermdiff_15 + org_frac * soiltype.thermdiff_15_org;
		thermdiff_100 = min_frac * soiltype.thermdiff_100 + org_frac * soiltype.thermdiff_100_org;

		// Two layer hydrology parameters
		perc_base = min_frac * soiltype.perc_base + org_frac * soiltype.perc_base_org;
		perc_exp = soiltype.perc_exp;
	}
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
double Soillayer::wcont() {

	return max(0.0, wmass - wp) / awc;
}

/// Field capacity of soil layer
double Soillayer::fc() {

	return awc + wp;
}

double Soillayer::fracair() {
	return (por - (wmass + imass)) / depth;
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

void Soillayer::update_moisture_diff_and_cond() {

	const double W_res = 0.01;

	if ((wsats - imass) / depth <= 1.001 * W_res) {
		Dw = 0.0;
		Kw = 0.0;
		return;
	}

	double W_sat, W_tot;

	// Account for ice filling up space
	W_sat = max(1.001 * W_res, (wsats - imass) / depth);
	W_tot = wmass / depth;

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

	double peat = 0.0;			// DW_COMMENT peat will be included when merge with 4.1

	// Soil layer air (mm)
	double air = por - (wmass + imass);

	if (air < -1.0e-8)
		fail("Illegal air mass value in update_heat_diff_and_cond(): %g", air);

	// Catch rounding errors for air
	air = max(0.0, air);

	// Update heat capacity (J m-3 K-1)
	Ch = mineral_frac * Cp_min + organic_frac * Cp_org +
		(imass * Cp_ice + wmass * Cp_water + peat * Cp_peat + air * Cp_air) / depth;

	// Update the thermal conductivity - which follows
	// Granberg et al. 1999, who refers to Farouki, 1986

	Ftot = mineral_frac + organic_frac + (imass + wmass + peat) / depth;

	if (Ftot + air / depth > 1.000001)
		fail("Illegal total fraction (%g) + air fraction (%g) value in update_heat_diff_and_cond(): %g", Ftot, air / depth, Ftot + air / depth);

	lKforg = (organic_frac / Ftot) * log(Korg);
	lKfpeat = (peat / depth / Ftot) * log(Kpeat);
	lKfmin = (mineral_frac / Ftot) * log(Kmin);
	lKfwater = (wmass / depth / Ftot) * log(Kwater);
	lKfice = (imass / depth / Ftot) * log(Kice);

	// Update thermal conductivity (J m-1 K-1 s-1)
	Kh = air / depth * Kair + Ftot * exp(lKfwater + lKfice + lKforg + lKfmin + lKfpeat);
}

// Get the Water Filled Pore Space of the layer
double Soillayer::wfps() {
	return wmass / wsats;
}

double Soillayer::wcont_max() {
	return (wsats - wp) / awc;
}

void Soillayer::serialize(ArchiveStream& arch) {
	arch & depth
		& wmass
		& imass
		& temp
		& temp_old;

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

	// Set soil layer depths (needs to be the same in init_rootdist()!!!) // DW_COMMENT
	std::vector<double> soil_layer_depths = patch.stand.get_gridcell().soil_layer_depths();

	nsoillayers = (int)soil_layer_depths.size();
	nsoillayersupper = 4; // DW_COMMENT check so it is right!

	// Set number of soil layers
	soillayers.resize(nsoillayers);

	double acc_depth = 0.0;

	bool acrotelm = false;
	bool catotelm = false;
	nacrotelm = 0;
	ncatotelm = 0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		acc_depth += soil_layer_depths[lyr];
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

		soillayers[lyr] = Soillayer(soil_layer_depths[lyr], acrotelm, catotelm);
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
	wcont_evap = 0.0;
	NH4_input = 0.0;
	NO3_input = 0.0;
	anfix = 0.0;
	anfix_mean = 0.0;
	anfix_calc = 0.0;
	awcont_upper = 0.0;
	rain_melt = 0.0;
	max_rain_melt = 0.0;
	dperc = 0.0;
	percolate = 0.0;
	runoff = 0.0;

	// Peatland hydrology variables
	awtp = 0.0;
	acro_co2 = 934.0; // [mimol L-1]
	wtd = 0.0; // [-100, +300] mm
	stand_wmass = 0.0;
	ch4_store = 0.0;
	co2_store = 0.0;

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

	// Solve LPJ SOM pools
	solvesom_end = static_cast<int>(0.8 * nyear_spinup);
	solvesom_begin = static_cast<int>(0.7 * nyear_spinup);

	// Solve CENTURY SOM pools
	solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
	solvesomcent_endyr = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);
	savesomcent = false;
}

std::vector<double> Soilcolumn::get_temp() {
	std::vector<double> temp(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		temp[lyr] = soillayers[lyr].temp;
	return temp;
}

std::vector<double> Soilcolumn::get_wmass() {
	std::vector<double> wmass(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wmass[lyr] = soillayers[lyr].wmass;
	return wmass;
}

std::vector<double> Soilcolumn::get_imass() {
	std::vector<double> imass(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		imass[lyr] = soillayers[lyr].imass;
	return imass;
}

std::vector<double> Soilcolumn::get_wcont() {
	std::vector<double> wconts(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wconts[lyr] = soillayers[lyr].wcont();
	return wconts;
}

std::vector<double> Soilcolumn::get_mwcont(int m) {
	std::vector<double> mwconts(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		mwconts[lyr] = soillayers[lyr].mwcont[m];
	return mwconts;
}

std::vector<double> Soilcolumn::get_awcont() {
	std::vector<double> awconts(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		awconts[lyr] = soillayers[lyr].awcont;
	return awconts;
}

std::vector<double> Soilcolumn::get_wsat() {
	std::vector<double> wsats(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wsats[lyr] = soillayers[lyr].wsats;
	return wsats;
}

std::vector<double> Soilcolumn::get_awc() {
	std::vector<double> awcs(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		awcs[lyr] = soillayers[lyr].awc;
	return awcs;
}

std::vector<double> Soilcolumn::get_wp() {
	std::vector<double> wps(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wps[lyr] = soillayers[lyr].wp;
	return wps;
}

std::vector<double> Soilcolumn::get_depth() {
	std::vector<double> depths(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		depths[lyr] = soillayers[lyr].depth;
	return depths;
}

std::vector<double> Soilcolumn::get_wcont_max() {
	std::vector<double> wcont_max(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		wcont_max[lyr] = soillayers[lyr].wcont_max();
	return wcont_max;
}

std::vector<double> Soilcolumn::get_NH4() {
	std::vector<double> NH4(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		NH4[lyr] = soillayers[lyr].NH4_mass;
	return NH4;
}

std::vector<double> Soilcolumn::get_NO3() {
	std::vector<double> NO3(nsoillayers);
	for (int lyr = 0; lyr < nsoillayers; lyr++)
		NO3[lyr] = soillayers[lyr].NO3_mass;
	return NO3;
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
	double acc_water = 0.0;
	double acc_awc = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		acc_water += soillayers[lyr].wcont() * soillayers[lyr].awc;
		acc_awc += soillayers[lyr].awc;
	}
	return acc_awc > 0.0 ? acc_water / acc_awc : 0.0;
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
void Soilcolumn::get_static_socfrac(std::vector<double> socfrac) {

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

	socfrac[0] = 1.0 - pow(beta_SOC, depth); // init first layer
	double soc_frac_total = socfrac[0];

	for (int lyr = 1; lyr < nsoillayers; lyr++) {
		depth += soillayers[lyr].depth * CM_PER_MM;
		socfrac[lyr] = 1.0 - pow(beta_SOC, depth) - (1.0 - pow(beta_SOC, depth - soillayers[lyr].depth * CM_PER_MM));
		soc_frac_total += socfrac[lyr];
	}

	// Put the residual SOC fraction in lowest soil layer
	socfrac[nsoillayers - 1] += 1.0 - min(soc_frac_total, 1.0);
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

	if (!iftwolayersoil)
		return;

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
		double depth, k, above_k, current_k;

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

void Soilcolumn::calc_wflow(double infil) {

	bool ifbedrock = false; // 1 (Bedrock) or 0 (free drainage)

	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		wpara.psi0[lyr] = soillayers[lyr].wmass / soillayers[lyr].depth;
	}

	for (int ts = 0; ts < date.timesteps; ts++) {

		// Calculate moisture diffusivity and conductivity
		update_moisture_diff_and_cond();

		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			wpara.dz[lyr] = soillayers[lyr].depth * M_PER_MM;
			wpara.diff[lyr] = soillayers[lyr].Dw;
			wpara.cond[lyr] = soillayers[lyr].Kw;
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

		double wmass_end = 0.0;
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			wmass_end += wpara.psi1[lyr] * soillayers[lyr].depth;
		}

		// Set new wmass
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].wmass = wpara.psi1[lyr] * soillayers[lyr].depth;

			// Calculate drainage

			// Drainage of water above saturation including ice
			double drained_wmass = max(0.0, soillayers[lyr].wmass + soillayers[lyr].imass - soillayers[lyr].wsats);

			soillayers[lyr].runoff_drain += drained_wmass;

			// Update total water after drainage 
			// remove drained water from wmass. Water and ice will now maximum fill up soil layer to saturation
			soillayers[lyr].wmass -= drained_wmass;
		}
	}

	// Calculate water fluxes to, inbetween and out of soil layers
	wfluxes[0] = infil;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		wfluxes[lyr + 1] = wfluxes[lyr] + wpara.psi0[lyr] * soillayers[lyr].depth - (soillayers[lyr].wmass + soillayers[lyr].runoff_drain);

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
	double aet_acrotelm = 0.0;  // transpiration from the acrotelm layers (mm)
	double runoff_surf;			// runoff from upper soil layer (mm)
	double runoff_drain;		// runoff (drainage) from lower soil layers (mm)
	double evap;				// evaporation from soil surface (mm)
	double evapotranspiration;	// evapotranspiration (mm/day)
	double fpc_moss_total = 0.0;// Total fpc for mosses

	// total water in acrotelm layers (mm)
	double wtot_acrotelm = to_depth(ACROTELM_DEPTH, get_wmass());
	// total ice in acrotelm layers (mm)
	double itot_acrotelm = to_depth(ACROTELM_DEPTH, get_imass());

	/// surface water content when WTD is at 20 mm
	const double surfw = 0.34;
	/// minimum fractional water content at surface in mm3/mm3 - see Wania et al (2009a), Eqn 22. Granberg et al. (1999), Eqns. 1-3.
	const double minvtot = 0.25;
	/// The soil water characteristics in the upper peat are considered to be linear in the suction interval 0-zmin and constant 
	/// in the interval zmin - acrotelm depth (normally 300mm). See: Granberg et al. (1999), Eqns. 1-3.
	const double zmin = 100.0;
	/// gradient limit of upper gradient [mm]
	const double grad_lim = 50.0;
	/// Granberg (1999) - Eqn 2 & Wania et al. (2009a), Eqn 22
	double az = (acrotelm_por - minvtot) / zmin;
	

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_wmass_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

	// *** TRANSPIRATION *** 

	// Retrieve Vegetation for this patch
	Vegetation& vegetation = patch.vegetation;

	// Sum AET for across all vegetation individuals in this patch

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		double pft_acro_root_frac = 0.0;

		for (int lyr = 0; lyr < nacrotelm; lyr++)
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
	// DW_COMMENT wtd updated?
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

	// Drainage - could possibly be read in for site-specific studies
	runoff_drain = 0; // DW_COMMENT CHECK HYDROLOGY PEAT AND DRAINAGE

	// Water in the acrotelm 
	wtot_acrotelm += rain_melt - evapotranspiration - runoff_drain;


	// *** RUNOFF AND RUNON***

	double ideal_runoff = exp(-0.01 * wtd);
	if (wtot_acrotelm > 0.0 && soillayers[0].imass / soillayers[0].depth < 0.7) // Granberg et al. (1999) use the same ice condition - see f_icestop in their Table 1
		runoff_surf = min(ideal_runoff, wtot_acrotelm);
	else
		runoff_surf = 0.0;

	// .ins file option: 
	// run on or off - but only when there's flowing water

	double runon = 0.0;

	if (runoff_surf > 0.0) {

		if (wetland_runon < 0)	// wetland_runon set in global.ins
			runoff_surf -= wetland_runon;	// e.g. 3 mm/day for bog-like/hummock conditions

		if (wetland_runon > 0)
			runon += wetland_runon;		// 3-5 mm/day for fen-like conditions
	}

	double wtot_init = wtot_acrotelm;

	wtot_acrotelm += runon - runoff_surf;

	// Add standing water, if any
	wtot_acrotelm += stand_wmass;

	// Id negative total water, then all initial water has gone to runoff and water is set to zero
	if (wtot_acrotelm < 0.0) {
		wtot_acrotelm = 0.0;
		runoff_surf = wtot_init;
	}

	// Save evapotranspiration data
	if (empirical_evapotranspiration) {

		patch.aaet += (1.0 - fevap) * evapotranspiration;
		patch.aevap += fevap * evapotranspiration;
		patch.maet[date.month] += (1.0 - fevap) * evapotranspiration;
		patch.mevap[date.month] += fevap * evapotranspiration;
	}
	else {

		patch.aaet += aet_total;
		patch.aevap += evap;
		patch.maet[date.month] += aet_total;
		patch.mevap[date.month] += evap;
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

	if (wtot_acrotelm > ACROTELM_DEPTH * acrotelm_por) {

		// Saturated acrotelm

		wtp[date.day] = wtot_acrotelm - ACROTELM_DEPTH * acrotelm_por;
		// i.e. wtp[date.day] > 0 above the surface

		if (wtp[date.day] > maxh) {
			runoff_surf += wtp[date.day] - maxh;
			// Update runoff from the surface
			wtp[date.day] = maxh;
		}

		wtd = -wtp[date.day]; // i.e. wtd defined as < 0 above the surface
		wtot_acrotelm = ACROTELM_DEPTH * acrotelm_por;
		// I.e. Wtot does NOT include standing water
		stand_wmass = 0.0; // = wtp[date.day]; // we assume no standing water 
	}
	else {

		// Non-saturated acrotelm - Granberg (1999), Eqns 1-5
		wtd = sqrt(3.0 * (ACROTELM_DEPTH * acrotelm_por - wtot_acrotelm) / (2.0 * az));

		if (wtd > zmin)
			wtd = 3.0 * (ACROTELM_DEPTH * acrotelm_por - wtot_acrotelm) / (2.0 * (acrotelm_por - minvtot));

		if (wtd > ACROTELM_DEPTH)
			wtd = ACROTELM_DEPTH;

		wtp[date.day] = -wtd;
		stand_wmass = 0.0;
	}


	// *** RUNOFF ***

	// Save today's runoff in Soil::runoff
	runoff = runoff_surf;
	patch.arunoff += runoff;
	patch.mrunoff[date.month] += runoff;
	patch.asurfrunoff += runoff;


	// *** WATER CONTENT IN EACH LAYER - THE PROFILE ***

	// Porosity updates in the presence of ice?
	bool subtractIceFromWtot = true;
	double avgIceFrac = subtractIceFromWtot ? to_depth(ACROTELM_DEPTH, get_imass(), get_depth()) : 0.0;

	// the pore space minus the ice fraction
	double acrotelm_por_icy = acrotelm_por - avgIceFrac;

	if (wtd > 0.0) {

		// Below the surface

		// Determine the subsurface layer (0 to nacrotelm) within which wtd lies. 
		
		int wtd_layer = 0;
		double acc_depth = soillayers[0].depth;
		while (acc_depth < wtd && wtd_layer < nacrotelm) {
			wtd_layer++;
			acc_depth += soillayers[wtd_layer].depth;
		}
		
		double surfw_sat = max(minvtot, acrotelm_por_icy - az * wtd); // Granberg - eqn 3. 
		std::vector<double> value(nacrotelm);

		if (wtd_layer == 0) {

			// Total saturation throughout the acrotelm (from 10mm down)
			for (int lyr = 0; lyr < nacrotelm; lyr++) {
				value[lyr] = soillayers[lyr].por;
			}
		}
		else {

			// Unsaturated zone in the acrotelm - use Granberg's eqn 1.
			double z = 0.0;
			for (int lyr = 0; lyr < wtd_layer; lyr++) {
				double z = soillayers[lyr].depth; // mm
				value[lyr] = min(soillayers[lyr].por, surfw_sat + (acrotelm_por_icy - surfw_sat) * pow(z / wtd, 2));
			}

			// Saturated zone in the acrotelm
			for (int lyr = wtd_layer; lyr < nacrotelm; lyr++) {
				value[lyr] = soillayers[lyr].por;
			}
		}

		// Now calculate the water content throughout the acrotelm
		soillayers[0].wmass = min(max(0.5 * (surfw_sat + value[0]) * soillayers[0].depth, 0.0), soillayers[0].por * soillayers[0].depth - soillayers[0].imass);
		for (int lyr = 1; lyr < nacrotelm; lyr++) {
			soillayers[lyr].wmass = min(max(0.5 * (value[lyr] + value[lyr - 1]) * soillayers[lyr].depth, 0.0), soillayers[lyr].por * soillayers[lyr].depth - soillayers[lyr].imass);
		}
	}
	// WTD <= 0.0, at or above the surface - standing water
	else {
		
		// acrotelm saturated (subtract ice content)
		for (int lyr = 0; lyr < nacrotelm; lyr++) {
			soillayers[lyr].wmass = max(soillayers[lyr].por * soillayers[lyr].depth - soillayers[lyr].imass, 0.0);
		}
	}

	// CATOTELM - assumed to be always saturated
	for (int lyr = nacrotelm; lyr < nsoillayers; lyr++) {
		soillayers[lyr].wmass = max(soillayers[lyr].por * soillayers[lyr].depth - soillayers[lyr].imass, 0.0);
	}


	// UPDATE CO2 LEVEL IN ACROTELM
	update_acrotelm_co2(climate.co2);


	// *** UPDATE LIMITS ON MOSS PHOTOSYNTHESIS ***

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


	// *** UPDATE LIMITS ON GRAMINOID PHOTOSYNTHESIS ***

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

	if (date.day == Date::MAX_YEAR_LENGTH - 1) {
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

	// Reset annuals
	if (date.day == 0) {
		patch.aevap = 0.0;
		patch.asurfrunoff = 0.0;
		patch.adrainrunoff = 0.0;
		patch.abaserunoff = 0.0;
		patch.arunoff = 0.0;
		patch.awetland_wmass_added = 0.0;
	}

	// Reset monthlys
	if (date.dayofmonth == 0) {
		patch.mevap[date.month] = 0.0;
		patch.mrunoff[date.month] = 0.0;
	}

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
		perc = min(min(1.0, SOILDEPTH_EVAP / soillayers[0].depth) * soillayers[0].perc_base * pow(wcont_evap, soillayers[0].perc_exp),
			max_rain_melt);
	}
	wcont_evap -= perc / (min(1.0, SOILDEPTH_EVAP / soillayers[0].depth) * soillayers[0].awc);

	// Percolation and fluxes to and from lower soil layer(s)

	// Transfer percolation between soil layers
	// Excess water transferred to runoff
	// Eqns 26, 27, 31, Haxeltine & Prentice 1996

	double runoff_drain = 0.0;

	for (int lyr = 1; lyr < nsoillayers; lyr++) {

		// Percolation

		// Allow only on days with rain or snowmelt
		if (percolate) {
			perc = min(soillayers[lyr - 1].perc_base * pow(soillayers[lyr - 1].wcont(), soillayers[lyr - 1].perc_exp), max_rain_melt);
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
			runoff_drain += soillayers[lyr].wmass - soillayers[lyr].fc();
			soillayers[lyr].wmass = soillayers[lyr].fc();
		}

		// Deduct AET from this soil layer
		soillayers[lyr].wmass -= min(soillayers[lyr].wmass - soillayers[lyr].wp, aet_layer[lyr]);
	}

	// Baseflow runoff (rain or snowmelt days only)
	double runoff_baseflow = 0.0;
	if (percolate) {
		double perc_baseflow = BASEFLOW_FRAC * soillayers[nsoillayers - 1].perc_base * 
			pow(soillayers[nsoillayers - 1].wcont(), soillayers[nsoillayers - 1].perc_exp);
		
		// Check so baseflow doesn't exceed inflow of water
		if (perc_baseflow > rain_melt - runoff_surf && rain_melt >= runoff_surf) {
			perc_baseflow = rain_melt - runoff_surf;
		}

		// Deduct from water content of bottom soil layer
		soillayers[nsoillayers - 1].perc_mass = min(perc_baseflow, soillayers[nsoillayers - 1].wmass - soillayers[nsoillayers - 1].wp);
		soillayers[nsoillayers - 1].wmass -= soillayers[nsoillayers - 1].perc_mass;
		runoff_baseflow = soillayers[nsoillayers - 1].perc_mass;
	}

	// Total runoff
	runoff = runoff_surf + runoff_drain + runoff_baseflow;

	// Water balance check
	double wafter = wmass() + rain_melt + runoff + aet_total + evap;
	if (fabs(wbefore-wafter) > 1.0e-10)
		dprintf("Year %d Day %d BUCKETWATER Before %g After %g Diff %g\n",date.year,date.day, wbefore, wafter, wbefore - wafter);

	// Save percolation from system (needed in leaching())
	dperc = runoff_baseflow + runoff_drain;

	// Save monthly and annaul values
	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain;
	patch.abaserunoff += runoff_baseflow;
	patch.arunoff += runoff;
	patch.aaet += aet_total;
	patch.aevap += evap;

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

	double wmass_begin = wmass() + rain_melt;

	// AET for a particular layer and individual (mm)
	double aet;
	// total AET for each soil layer (mm)
	std::vector<double> aet_layer(nsoillayers);
	// total AET
	double aet_total = 0.0;

	// Actual Evapotranspiration

	// Sum AET for across all vegetation individuals and extract it from soil lalyers
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			double fwu = patch.pft[indiv.pft.id].fwuptake[lyr]; // DW_COMMENT not needed but crashed otherwise???????
			aet = patch.pft[indiv.pft.id].fwuptake[lyr] * indiv.aet;
			soillayers[lyr].wmass -= aet;
			aet_layer[lyr] += aet;
			aet_total += aet;
		}
		vegetation.nextobj();
	}

	// Evaporation

	// Calculate wcont and awc for evaporation layers (only up to awc)
	double wcont_evap = min(1.0, to_depth(SOILDEPTH_EVAP, get_wcont(), get_awc()));

	// Calculate total evaporation
	// guess2008 - changed to wcont_evap**2, as in LPJ-mL
	// - see Bondeau et al. (2007),  Rost et al. (2008)
	// Added the snowdepth restriction too.
	// evap only if snow depth < 10mm
	double evap = 0.0;
	if (snow.depth() < 10.0) {
		evap = climate.eet * PRIESTLEY_TAYLOR * wcont_evap * wcont_evap * fevap;
	}

	// Evaporation cannot exceed available water down to evaporation depth
	double avail_evap = to_depth(SOILDEPTH_EVAP, get_wmass()) - to_depth(SOILDEPTH_EVAP, get_wp());
	evap = min(avail_evap, evap);

	// Extract evaporation from soil layers as a fraction of available water
	if (!negligible(evap)) {
		double acc_depth = 0.0;
		int lyr = 0;

		while (acc_depth < SOILDEPTH_EVAP) {
			// Only use fraction of layer if SOILDEPTH_EVAP is in the middle of a layer
			double lyr_depth_frac = (acc_depth + soillayers[lyr].depth < SOILDEPTH_EVAP ? soillayers[lyr].depth : SOILDEPTH_EVAP - acc_depth) / soillayers[lyr].depth;
			soillayers[lyr].wmass -= evap * lyr_depth_frac * soillayers[lyr].wmass / avail_evap;
			acc_depth += soillayers[lyr].depth;
			lyr++;
		}
	}

	// Surface infiltration and runoff

	// Calculate water and ice content for infiltration layers
	double wmass_infilt = to_depth(SOILDEPTH_INFILT, get_wmass());
	double imass_infilt = to_depth(SOILDEPTH_INFILT, get_imass());
	double wsats_infilt = to_depth(SOILDEPTH_INFILT, get_wsat());
	double totwcont_infilt = (wmass_infilt + imass_infilt) / wsats_infilt;

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

	patch.asurfrunoff += runoff_surf;
	patch.adrainrunoff += runoff_drain();
	patch.abaserunoff += runoff_baseflow();
	patch.arunoff += runoff;
	patch.awetland_wmass_added += patch.wetland_wmass_added_today;
	patch.aaet += aet_total;
	patch.aevap += evap;

	patch.maet[date.month] += aet_total;
	patch.mevap[date.month] += evap;
	patch.mrunoff[date.month] += runoff;

	// Calculate annual wcont for growing season (temperature > 5C)
	// guess2008 - DLE - update awcont
	// Original algorithm by Thomas Hickler
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
			if (patch.growingseasondays < 1) {
				soillayers[lyr].awcont = 1.0;
			}
			else {
				soillayers[lyr].awcont /= (double)patch.growingseasondays;
			}
		}
	}

	double wmass_end = wmass() + runoff + evap + aet_total;

	double EPS = 1.0e-5;
	
	if (fabs(wmass_begin - wmass_end) > EPS)
		dprintf("Year %d Day %d WMASS ERROR: Diff %g Begin %g End %g\n", date.year, date.day, wmass_begin - wmass_end, wmass_begin, wmass_end);
}


// Initial percolation
/* Initial percolation
 */
void Soilcolumn::initial_percolation() {

	// Reset annuals
	if (date.day == 0) {
		patch.awetland_wmass_added = 0.0;
	}

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

	double Ffreez;	// fraction of water which freezes
	double Fthaw;	// fraction of ice which thaws
	double energy;	// Energy per volume [J m-3]

	const double SMALLVOLFRAC = 0.0000001;	// Minimum volume fraction [m-3] allowed

	double delta_T;	// temperature difference (current - previous)

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Initial water and ice content
		double wtotal_initial = soillayers[lyr].wmass + soillayers[lyr].imass;

		// Water and ice fraction of full volume
		double Fwmass = soillayers[lyr].wmass / soillayers[lyr].depth;
		double Fimass = soillayers[lyr].imass / soillayers[lyr].depth;

		// temperature difference
		delta_T = soillayers[lyr].temp - soillayers[lyr].temp_old;

		// FREEZING
		if (soillayers[lyr].temp < 0.0 && soillayers[lyr].temp < soillayers[lyr].temp_old) { // Wania conditions 
			Fthaw = 0.0;
			if (Fwmass > 0.0) {
				energy = Cp_water * fabs(delta_T); // energy in J m-3 
				Ffreez = energy / Lheat; // unitless fraction

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
					soillayers[lyr].temp += Fwmass * Lheat / soillayers[lyr].Ch;
					soillayers[lyr].temp = min(0.0, soillayers[lyr].temp);

					Fwmass = 0.0;
				}
			}
			else {
				Ffreez = 0.0;
			}
		} // end of FREEZING

		  // THAWING
		else if (soillayers[lyr].temp >= 0.0  && soillayers[lyr].temp > soillayers[lyr].temp_old) { // Wania conditions
			Ffreez = 0.0;
			if (Fimass > 0.0) {
				energy = Cp_water * fabs(delta_T); // energy in J m-3
				Fthaw = energy / Lheat; // unitless fraction
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

					// Update soil temperature in this layer (ensures that we don't go back under 0 degrees)
					//soillayers[lyr].temp = max(0.0, (Fthaw - Frac_imass) * Lheat / soillayers[lyr].Ch);
					soillayers[lyr].temp -= Fimass * Lheat / soillayers[lyr].Ch;
					soillayers[lyr].temp = max(0.0, soillayers[lyr].temp);

					Fimass = 0.0;
				}
			}
			else {
				Fthaw = 0.0;
			}
		} // end of THAWING


		if (Fimass < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid ice volume content %g\n", Fimass);

		if (Fwmass < -1.0 * SMALLVOLFRAC)
			fail("Soil::update_ice_fraction - invalid volumetric water content %g\n", Fwmass);

		// Water and ice fraction of full volume
		soillayers[lyr].wmass = Fwmass * soillayers[lyr].depth;
		soillayers[lyr].imass = Fimass * soillayers[lyr].depth;

		double wtotal_final = soillayers[lyr].wmass + soillayers[lyr].imass;

		// water balance checks:
		if ((fabs(wtotal_initial - wtotal_final)) > SMALLVOLFRAC)
			fail("Water mass balance during freezing or thawing in Soil::phase_change. Init %g, Final %g, Diff %g \n", wtotal_initial, wtotal_final, wtotal_initial - wtotal_final);

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

	if (!iftwolayersoil) {

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
				if (lyr == 0 && snow.snowlayers[slyr].depth * M_PER_MM < min_slyr_dz && !iftwolayersoil) {
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
		& rain_melt
		& max_rain_melt
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
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Levine, J. S. (1996) Biomass Burning and Global Change. Remote Sensing, Modeling
//   and Inventory Development, and Biomass Burning in Africa, 1J. S. Levine,
//   XXXV-XLIII, MIT Press, Mass.
// For Arctic and wetland references see soil.h
