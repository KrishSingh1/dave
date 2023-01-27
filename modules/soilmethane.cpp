///////////////////////////////////////////////////////////////////////////////////////
/// \file methane.cpp
/// \brief Methane calculations and transport.
/// Implementation of member functions of class Soil.
/// The class Soil and its member functions and variables are declared in guess.h
///
/// \author Paul Miller
/// $Date: 2017-05-203 13:32:28 +0200 (Wed, 03 May 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include "soilmethane.h"
#include "guess.h"
#include "soil.h"


//////////////////////////////////////////////////////////////////////////////////////////////////
// METHANE
//////////////////////////////////////////////////////////////////////////////////////////////////


/// Calculate methane emissions and update methane concentrations in soil layers. 
/**
	Emissions if and only if ifmethane true (in global.ins) and this is a PEATLAND stand
	and we have passed the initial years of the spinup.

	Called daily from framework.

*/
void methane_dynamics(Patch& patch) {

	if (ifmethane && patch.stand.landcover == PEATLAND && date.year > nyear_spinup-100) {
		patch.soilcolumn.methane(true);
	} 
	else {
		patch.soilcolumn.methane(false);
	}
}

/// Crank-Nicholson timestepper algorithm for gas diffusion equation
/** Crank-Nicholson timestepper algorithm for gas diffusion equation
 * UNITS:
 * Di			m2 d-1
 * dz			m
 * surf_conc	mmol m-3
 * dt			d
 * conc			mmol m-3
*/
void cnstepgas(double Di[NSOILLAYERS_MAX], double dz[NSOILLAYERS_MAX], double surf_conc, double dt,
			double conc[NSOILLAYERS_MAX], int nsoillayers) {

	double dplus;
	double dminus;
	double dz_factor;
	double Cplus; 
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double cohere, cominus, coplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson matrix.
	double diag[NSOILLAYERS_MAX];
	double left[NSOILLAYERS_MAX];
	double right[NSOILLAYERS_MAX];

	// Right hand side vector for Crank-Nicholson scheme equations.
	double rhs[NSOILLAYERS_MAX];

	// Solution vector for Crank-Nicholson scheme equations.
	double solution[NSOILLAYERS_MAX];

	// Initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = cohere = cominus = coplus = 0.0;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[nsoillayers - 1] = 0.0;

	// Process the active layers.

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  

		// D+

		if (lyr == nsoillayers - 1) {
			// BC2 - Bottom layer diffusion clamped to 0
			dplus = 0.0;
		}
		else {
			dplus = 0.5 * (Di[lyr] + Di[lyr + 1]);
		}

		// D-

		if (lyr == 0) {
			// Top layer
			dminus = Di[lyr];
		}
		else {
			// Soil layers
			dminus = 0.5 * (Di[lyr] + Di[lyr - 1]);
		}
 
		// --- HERE ---
		if (lyr < nsoillayers) {
			// All soil layers
			dzhere = dz[lyr];
			cohere = conc[lyr];
		}

		// --- PLUS ---
		if (lyr < nsoillayers - 1) {
			// All soil layers apart from the bottom soil layer
			dzplus = dz[lyr + 1];
			coplus = conc[lyr + 1];
		}

		// --- MINUS ---
		if (lyr == 0) {
			// Top layer
			dzminus = dz[lyr];
			cominus = conc[lyr];
		}
		else if (lyr < nsoillayers) {
			// All soil layers and the first padding layer
			dzminus = dz[lyr - 1];
			cominus = conc[lyr - 1];
		}

		dz_factor = 0.25 * (dzplus + 2.0 * dzhere + dzminus);
		Cplus = dplus * dt / dz_factor / (dzplus + dzhere);
		Cp_minus = dminus * dt / dz_factor / (dzhere + dzminus);

		// Fill in matrix diagonal and off-diagonal elements.
 
		// DIAG
		if (lyr == 0) {
			// BC1 - top layer should be (1,0,...,0)
			diag[0] = 1.0;
		}
		else {
			diag[lyr] = 1.0 + Cplus + Cp_minus;
		}
 
		// LEFT & RIGHT
		if (lyr < nsoillayers - 1) {

			if (lyr > 0) {
				right[lyr] = -Cplus;
			}
			else {
				// i.e. BC1, where top layer == (1,0,..,0)
				right[lyr] = 0.0;
			}

			// left[0] is set above.
			if (lyr > 0) {
				left[lyr] = -Cp_minus;
			}
		}

		if (lyr == nsoillayers - 1) {
			left[lyr] = -Cp_minus;
		}
		 
		// RHS
		// Calculate right hand side vector values.
		if (lyr == 0) {
			rhs[0] = surf_conc;
		}
		else if (lyr == nsoillayers - 1) {
			// Cplus == 0 here anyway
			rhs[lyr] = (1.0 - Cp_minus) * cohere + Cp_minus * cominus;
		}
		else {
			rhs[lyr] = (1.0 - Cplus - Cp_minus) * cohere +
				Cplus * coplus + Cp_minus * cominus;
		}
	} // end for

	//   SOLVE TRIDIAGONAL SYSTEM
	tridiag(nsoillayers, left, diag, right, rhs, solution);

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the concentration array.
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		conc[lyr] = solution[lyr];
	}
}

/// Calculate and return diffusivities
/** Calculate and return diffusivities, in units of m2 d-1
 *	Called each day
 *	See Wania et al. (2010) - Sec 2.5
 */
void Soilcolumn::calculate_gas_diffusivities(double dCH4[NSOILLAYERS_MAX], double dCO2[NSOILLAYERS_MAX], 
	double dO2[NSOILLAYERS_MAX]) {

	double D_CH4_air, D_CO2_air, D_O2_air;

	const double scale = 0.00001;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		dCH4[lyr] = 0.0;
		dCO2[lyr] = 0.0;
		dO2[lyr] = 0.0;
		
		// Soil layer temperature
		double temp = soillayer.temp;

		// Fraction air in soil layer
		double frac_air = soillayer.fracair();
	
		if (frac_air > 0.05 && soillayer.acrotelm) {

			// AIR - UNITS: 10-4 m2 s-1 == cm2 s-1
			// Lerman (1979)

			// Wania et al. (2010), Eqn. 10
			D_CH4_air = 0.1875 + 0.00130 * temp;
			D_CO2_air = 0.1325 + 0.00090 * temp;
			D_O2_air = 0.1759 + 0.00117 * temp;

			// Wania et al. (2010), Eqn. 11 
			double airpow = pow(frac_air, 10 / 3) / pow(acrotelm_por, 2.0);
			dCH4[lyr] = airpow * D_CH4_air;
			dCO2[lyr] = airpow * D_CO2_air;
			dO2[lyr] = airpow * D_O2_air;
		}
		else {

			// Water diffusivities dominate
			// WATER - UNITS: 10-4 m2 s-1 == cm2 s-1
			// Wania et al. (2010), Eqn. 9 - where the diffusivities are in units of 10**-9 m2 s-1. 
			// * by scale to get 10-4 m2 s-1
			dCH4[lyr] = (0.9798 + 0.02986 * temp + 0.0004381 * temp * temp) * scale;
			dCO2[lyr] = (0.9390 + 0.02671 * temp + 0.0004095 * temp * temp) * scale;
			dO2[lyr] = (1.1720 + 0.03443 * temp + 0.0005048 * temp * temp) * scale;
		}
		
		// Convert from cm2 s-1 to m2 d-1
		dCH4[lyr] *= SECS_PER_DAY / CM2_PER_M2;
		dCO2[lyr] *= SECS_PER_DAY / CM2_PER_M2;
		dO2[lyr] *= SECS_PER_DAY / CM2_PER_M2;
	}
}

void Soilcolumn::calculate_ch4_and_co2_production() {

	// Conservation check - C in methane produced
	double c_input = 0.0;
	double initial_CH4 = get_sum(CH_4);
	double initial_CO2 = get_sum(CO_2) + get_sum(RESPSUM) * G_PER_KG;

	double CH4_prod, CO2_prod;

	// Loop through soil layers
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Daily decomposition, converted to gC/m2 from kgC/m2
		double drh = soillayer.respsum * G_PER_KG;

		double anoxic = 1.0 - soillayer.fracair() - Fgas;

		// CH4 production in the layer (gC/m2)
		if (soillayer.wmass / soillayer.depth < Fwmass_min)
			CH4_prod = 0.0;
		else
			CH4_prod = anoxic * CH4toCO2_peat * drh;

		// Add CH4 production to CH4 pool 
		soillayer.CH4 += CH4_prod;

		// CO2 production
		CO2_prod = drh - CH4_prod;

		// CO2 pool
		// Add CO2 production to CO2 pool 
		soillayer.CO2 += CO2_prod;
	}

	double total_CH4 = get_sum(CH_4);
	double total_CO2 = get_sum(CO_2);
	double respsum_tot = get_sum(RESPSUM);

	if (fabs((initial_CH4 + initial_CO2) - (total_CH4 + total_CO2)) > 1.0e-12)
		dprintf("Year %d Day %d ERROR: CH4 production %g\n",date.year,date.day, (initial_CH4 + initial_CO2) - (total_CH4 + total_CO2));
}

/// Update gas transport velocities and equilibrium gas concentrations at surface interface
/** Update gas transport velocities [m d-1] and equilibrium gas concentrations [mmol m-3]
 * at the water - air interface boundary at the surface interface
 * for CH4, CO2 and O2
 * Called daily by Soil::methane()
 * See Wania et al. (2010) - Sec 2.5
 */
void Soilcolumn::update_daily_gas_parameters() {

	double surfT = soillayers[0].temp;
	double U10 = 0.0; // patch.get_climate().u10

	// Gas transport velocities

	// Gas transport velocity of SF6.
	// Wania et al. (2010) - Eqn. 6 (cm h-1)
	double k_600 = 2.07 + 0.215 * pow(U10, 1.7);

	// Schmidt number of O2
	// Wania et al. (2010) - Eqn. 7
	double ScO2 = 1800.6 - 120.1 * surfT + 3.7818 * pow(surfT, 2) - 0.047608 * pow(surfT, 3);

	// Gas transport velocity of O2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// Converted from cm h-1 to m d-1
	k_O2 = k_600 * pow(ScO2 / 600.0, n_coeff) * HOUR_PER_DAY * M_PER_CM;

	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCH4 = 1898.0 - 110.1 * surfT + 2.834 * pow(surfT, 2) - 0.02791 * pow(surfT, 3);

	// Gas transport velocity of CH4 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// Converted from cm h-1 to m d-1
	k_CH4 = k_600 * pow(ScCH4 / 600.0, n_coeff) * HOUR_PER_DAY * M_PER_CM;

	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCO2 = 1911.0 - 113.7 * surfT + 2.967 * pow(surfT, 2) - 0.02943 * pow(surfT, 3);

	// Gas transport velocity of CO2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// Converted from cm h-1 to m d-1
	k_CO2 = k_600 * pow(ScCO2 / 600.0, n_coeff) * HOUR_PER_DAY * M_PER_CM;

	// Equilibrium gas concentration
	// See Wania et al. (2010) - Eqn. 8

	// 25 degrees C
	double deg25 = K2degC + 25.0;

	// Henry's Law [L atm mol-1]
	double henry_coeff_O2 = henry_k_O2 * exp(-1.0 * henry_C_O2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for O2
	// pp_gas/MM2_PER_M2 (mol L-1) converts to atm units
	// mol L-1 to mmol m-3
	Ceq_O2 = pp_O2 / henry_coeff_O2 * MMOL_PER_MOL * M3_PER_DM3;

	// Henry's Law [L atm mol-1]
	double henry_coeff_CO2 = henry_k_CO2 * exp(-1.0 * henry_C_CO2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for CO2
	// Use this gridcell's CO2 concentration, not a fixed value as in Wania et al. (2010).
	// mol L-1 to mmol m-3
	Ceq_CO2 = patch.get_climate().co2 / henry_coeff_CO2 * MMOL_PER_MOL * M3_PER_DM3;

	// Henry's Law [L atm mol-1]
	double henry_coeff_CH4 = henry_k_CH4 * exp(-1.0 * henry_C_CH4 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for CH4
	// mol L-1 to mmol m-3
	Ceq_CH4 = pp_CH4 / henry_coeff_CH4 * MMOL_PER_MOL * M3_PER_DM3;
}


void Soilcolumn::diffuse_gas(double D[NSOILLAYERS_MAX], gastype thisgastype, double Ceq,
					  double kgas, double &dailyDiff) {

	// Generic gas diffusion method that works with CO2, CH4 and O2 (as specified with gastype) 
	// See Wania et al. (2010) - Sec 2.5 - for a full description

	// Called like this (e.g. for O2): 
	// diffuse_gas(O2, D_O2, O2gas, Ceq_O2, k_O2, dailyO2diffusion);

	// No diffusion until snow depth < 50mm 
	// Top boundary condition should be zero, but diffusion within the unfrozen soil
	if (snow.depth() >= 50.0)
		return;
	
	// Atomic mass of this gas [gC/mol]
	double atomic_mass;

	// Concentration of the dissolved gas in question [mmol m-3]
	double C[NSOILLAYERS_MAX];

	// Soil layer thickness in meters
	double Dz_m[NSOILLAYERS_MAX];

	// layer water + ice volume [m3]
	double total_volume_water[NSOILLAYERS_MAX];
	
	if (thisgastype != O2gas)
		atomic_mass = atomiccmass; // i.e. CO2 or CH4 - 12 gC/mol
	else
		atomic_mass = 1.0; // O2 already in mol layer-1
	
	double initialAmount = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		Dz_m[lyr] = soillayer.depth * M_PER_MM;

		total_volume_water[lyr] = (soillayer.wmass + soillayer.imass) * M_PER_MM;
		if (total_volume_water[lyr] > 1.0e-6) {
			if (thisgastype == CO2gas) {
				initialAmount += soillayer.CO2;
				// g layer-1 to mmol m-3
				C[lyr] = soillayer.CO2 / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
			}
			else if (thisgastype == CH4gas) {
				initialAmount += soillayer.CH4;
				// g layer-1 to mmol m-3
				C[lyr] = soillayer.CH4 / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
			}
			else {
				initialAmount += soillayer.O2;
				// mol layer-1 to mmol m-3
				C[lyr] = soillayer.O2 / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
			}
		}
		else {
			C[lyr] = 0.0;
		}
	}

	// Set the BC, i.e. equilibrium gas concentrations in the top layer depending on atmospheric concentrations
	// and using Henry's law. See Wania et al. (2010), Eqs. 4-8
	
	// New surface concentration [mmol m-3]
	double Cnew;
	if (thisgastype == O2gas || thisgastype == CH4gas) { // Could also run for CH4 here below.

		Soillayer& soillayer_surf = soillayers[0];

		// No diffusion if there's too little liquid water in the top layer
		// Could add a snow restriction
		if ((soillayer_surf.wmass / soillayer_surf.depth) < Fwmass_min) {
			// Unchanged surface concentration
			Cnew = C[0]; 
			dailyDiff = 0.0;
		} 
		else {
			double volume_liquid_water = soillayer_surf.wmass * M_PER_MM;

			// Analytical solution to determine Csurf - see Wania et al. Sec 2.5 (mmol m-3)
			Cnew = Ceq + (C[0] - Ceq) * exp(-kgas / volume_liquid_water);
			// mol layer-1 d-1 (O2)
			dailyDiff = (C[0] - Cnew) * atomic_mass * volume_liquid_water / MMOL_PER_MOL;
		}
		// The new surface concentration (mmol m-3) 
		C[0] = Cnew; 
	}

	double surf_conc = C[0];

	// Diffusion of gas today?
	if (dailyDiff != 0) {

		// Determine the timestep to use

		double maxD = -0.01;
		// Maximum diffusivity today below the surface layer [m2 d-1]
		// Excluding top layer as we've already considered in the top layer above
		for (int lyr = 1; lyr < nsoillayers; lyr++) { 
			if (D[lyr] > maxD) {
				maxD = D[lyr];
			}
		}

		// Maximum timestep [day]
		double max_timestep = (0.1 * 0.1) / (2.0 * maxD);

		if (max_timestep < Dt_gas && verbosity >= INFO)
			dprintf("%s%8.4f\n","Max gas diffusion timestep exceeded: ",maxD);

		// Assume a tiny diffusivity if there is very little liquid water in the layer
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			Soillayer& soillayer = soillayers[lyr];
			if ((soillayer.wmass / soillayer.depth) < Fwmass_min) {
				D[lyr] = 1e-9; // [m2 d-1]
			}
		}

		double C_init[NSOILLAYERS_MAX];
		double C_last[NSOILLAYERS_MAX];

		// Determine initial and total concentrations [mmol m-3]
		double totalConc = 0.0; 
		double total_diff = 0.0;
		double ratio = 0.0;
		int cncount = 0;

		// Loop until stable concentrations are found, or until 100 iterations have been performed
		int gasdiffix = int(1 / Dt_gas);
		
		bool stable = false;

		do {

			totalConc = 0;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				C_init[lyr] = C[lyr];
				totalConc += C_init[lyr];
			}
			
			// Use the CN scheme
			cnstepgas(D, Dz_m, surf_conc, Dt_gas, C, nsoillayers);

			total_diff = 0.0;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				C_last[lyr] = C[lyr];
				total_diff += fabs(C_last[lyr] - C_init[lyr]);
			}

			cncount++;

			ratio = total_diff / totalConc;
			
			// Stable if the ratio is < 0.01 and we've passed a third of the max iterations 
			if (ratio < 0.01 && cncount > int(gasdiffix / 3))
				stable = true;

		} while (cncount < gasdiffix && !stable); // until max 1% variation between iterations, or 100 loops
	}

	// Unit conversions - allocate concentrations to g layer-1 or mole layer-1
	double finalAmount = 0.0;
	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		if (thisgastype == CO2gas) {
			soillayer.CO2 = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayer.CO2;
		}
		else if (thisgastype == CH4gas) {
			soillayer.CH4 = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayer.CH4;
		}
		else {
			soillayer.O2 = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayer.O2;
		}
	}

	// Calculate the amount of diffused diffused to the atmosphere
	dailyDiff = initialAmount - finalAmount;
}

/// Calculate the area and density of tillers
/** Calculate the area and density of tillers
 * Loop through the individuals present and sum the leaf carbon mass for C3-graminoids [m-2] 
 * See Wania et al. (2010) - Sec 2.6
*/
void Soilcolumn::calculate_tiller_areas(double& tiller_density) {

	// Tiller density [number m-2]
	tiller_density = 0.0; 

	double graminoid_leafcmass = 0.0;
	double graminoid_dphen = 0.0;
	double graminoid_anpp_red = 0.0;
	double root_mass = 0.0;
	double rootfrac[NSOILLAYERS_MAX];

	// Loop through this patch object's vegetation and pick out graminoids

	Vegetation& vegetation=patch.vegetation;

	vegetation.firstobj();
	while (vegetation.isobj) {
	
		Individual& indiv=vegetation.getobj();
	
		// Only allow plant-mediated transport from graminoids that are alive
		if (indiv.alive && indiv.pft.has_aerenchyma) {
			// e.g. C3 graminoids on wetlands
			graminoid_leafcmass += indiv.cmass_leaf;
			graminoid_dphen += indiv.phen;
			graminoid_anpp_red += indiv.anpp;
			root_mass += indiv.cmass_root;
			for (int lyr = 0; lyr < nsoillayers; lyr++) {
				rootfrac[lyr] += indiv.cmass_root * indiv.pft.rootdist[lyr];
			}
		}
		vegetation.nextobj();
	}

	// Root fractions per layer
	for (int lyr = 0; lyr < nsoillayers; lyr++) {
		rootfrac[lyr] /= root_mass;
	}

	// Convert leaf C mass from kg to gram (kgC m-2 to gC m-2)
	graminoid_leafcmass *= G_PER_KG;

	// Tiller density [number m-2]
	tiller_density = graminoid_leafcmass * graminoid_dphen / tiller_weight;

	double tiller_frac = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		if (!negligible(tiller_density, 1.0e-15)) {
			tiller_frac = tiller_density * rootfrac[lyr];

			// Consider this as the fraction of the area (1m2) occupied by tillers in this layer
			soillayer.tiller_area = PI * pow(tiller_radius, 2.0) * tiller_frac * tiller_por; // [unitless m2 m-2]

			// Neglect small values
			if (soillayer.tiller_area < 0.000000001)
				soillayer.tiller_area = 0.0;
		}
		else {
			soillayer.tiller_area = 0.0;
		}
	}
}

/// Plant transport of O2 or CH4
/** Plant transport of O2 or CH4
 * The limiting factor to plant transport are the number of tillers as
 * they represent the cross-sectional area available to gas transport.
 * The biomass is related to tiller density based on Schimel 1995. The
 * average biomass was 185g/m2 and the average number of tillers was 380/m2.
 * This gives us 2.05 tillers per g biomass and 0.48 g per tiller, which
 * corresponds to 0.22 g C per tiller
 * See Wania et al. (2010) - Sec 2.6
 */
void Soilcolumn::plant_gas_transport(double Ceq, double kgas, gastype thisgastype, 
							 double tillers, double& plantTransportToday) {

	// Only allow plant transport if there are graminoid tillers present
	if (negligible(tillers))
		return;

	double Cgas;
	double Cnew = 0.0; // mmol m-3
	double atomic_mass;

	if (thisgastype != O2gas)
		// i.e. CO2 or CH4 - 12 g/mol
		atomic_mass = atomiccmass;
	else
		atomic_mass = 1.0; 

	// to return
	plantTransportToday = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Weight plant transport by the area of porous root cross-sections. 
		// plant_trans(i) [mmol layer-1 d-1]
        // the 2.07 are from the equation for k_CH4, when U10 is zero.

		double total_volume_water = (soillayer.wmass + soillayer.imass) * M_PER_MM;
		double volume_liquid_water = soillayer.wmass * M_PER_MM;

		// CH4 and CO2: gC layer-1 to mmol m-3
		// O2: mol layer-1 to mmol m-3
		if (thisgastype == O2gas) {
			Cgas = soillayer.O2 / atomic_mass / total_volume_water * MMOL_PER_MOL;
		}
		else {
			Cgas = soillayer.CH4 / atomic_mass / total_volume_water * MMOL_PER_MOL;
		}

		// No plant transport to/from this layer
		if (soillayer.wmass / soillayer.depth < Fwmass_min || negligible(soillayer.tiller_area) || volume_liquid_water < 0.0001) {
			Cnew = Cgas;
		} 
		else {

			double exponent = kgas / volume_liquid_water * soillayer.tiller_area;
			Cnew = Ceq + (Cgas - Ceq) * exp(-exponent); // mmol m-3

			// mol layer-1 d-1 (O2) or gC m-2 d-1 (CH4 & CO2)
			plantTransportToday += (Cgas - Cnew) 
				* atomic_mass * total_volume_water / MMOL_PER_MOL;

			// Restrict TINY negative plant transport
			if (plantTransportToday < 0.0 && fabs(plantTransportToday) < 1e-8) {
				plantTransportToday = 0.0;
				Cnew = Cgas;
			}
		}

		// mmol m-3 to gC layer-1 (CO2 or CH4), or mol layer-1 (O2) 
		if (thisgastype != O2gas) {
			soillayer.CH4 = Cnew * atomic_mass * total_volume_water / MMOL_PER_MOL;
		}
		else {
			soillayer.O2 = Cnew * atomic_mass * total_volume_water / MMOL_PER_MOL;
		}
	}
}

/// CH4 oxidation
/** CH4 oxidation
 * CH4 units: gC layer-1 
 * O2 units: mol layer-1
 * Assume that 1/2 of the O2 is utilized by other electron
 * acceptors and 2 moles of O2 is needed to oxidise 1 mole of CH4.
 * O2 is in [mol layer-1]
 */
void Soilcolumn::oxidation_ch4() {

	double CH4_mole = 0.0;
	double CH4_oxid = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		// Root use (1-oxid_frac) of O2 [mol layer-1] 
		soillayer.O2 *= oxid_frac;
		// Convert CH4 from gC to mol layer-1 
		CH4_mole = soillayer.CH4 / atomiccmass;
		// Calculate CH4 oxidation (usually ~75%) [g CH4-C layer-1 d-1]
		CH4_oxid = min(CH4_mole, 0.5 * soillayer.O2);
		// Calculate new CH4 level in gC
		soillayer.CH4 -= CH4_oxid * atomiccmass;
		// Subtract the moles used in oxidation
		soillayer.O2 -= 2.0 * CH4_oxid;
		// Oxidised CH4 becomes CO2 (and water!)
		soillayer.CO2 += CH4_oxid * atomiccmass;

		//dprintf("Watermass added %g\n",CH4_oxid * 2.0 * molecularmass_H2O * KG_PER_G);
	}
}

/** Gas ebullition
 * Gas ebullition
 * Returns the gas ebullition today
 * CH4 directly emitted to the atmosphere
 * See Wania et al. (2010) - Sec 2.7
*/
void Soilcolumn::calculate_gas_ebullition(double& ebull_today) {

	// Initial mass
	double initial_CH4 = get_sum(CH_4);

	// Ebullition to emit today [g m-2 d-1]
	ebull_today = 0.0;
	
	double waterheight = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		double TsoilK = soillayer.temp + K2degC;

		double total_volume_water = (soillayer.wmass + soillayer.imass) * M_PER_MM;
		double volume_liquid_water = soillayer.wmass * M_PER_MM;

		// Water height in meters
		waterheight += total_volume_water;

		// Max CH4 that can be dissolved (Wania et al. (2010), Eqn 15) [ml CH4 ml-1 H2O]
		double CH4_diss_max = 0.05708 - 0.001545 * max(0.0, soillayer.temp) +
			0.00002069 * max(0.0, pow(soillayer.temp, 2));

		// Water pressure [N, or kg m s-2]
		double hydro_press = rho_H2O * gravity * waterheight;

		// Maximum volume that can be dissolved [m3 CH4 layer-1}
		double CH4_diss_max_m3 = CH4_diss_max * volume_liquid_water;

		// Use ideal gas law to convert to mol CH4 layer-1
		double CH4_diss_max_mol = CH4_diss_max_m3 * (hydro_press + atm_press) / (R_gas * TsoilK);

		// Convert to g CH4-C layer-1
		double CH4_diss_max_g = CH4_diss_max_mol * atomiccmass;

		//double henry_k_cc_CH4 = TsoilK / (12.2 * henry_k_CH4);
		//soillayer.CH4_diss = min(CH4_diss_max_g, henry_k_cc_CH4 * soillayer.CH4);

		// Override (Wania et al. 2010) with this new, simpler ebullition
		soillayer.CH4_diss = min(CH4_diss_max_g, soillayer.CH4);

		soillayer.CH4_gas = soillayer.CH4 - soillayer.CH4_diss; // g CH4-C layer-1 - gas

		// Reset CH4 which bubbles out
		soillayer.CH4_ebull_ind = 0.0;

		// Restrict ebullition to cases when there is enough liquid water and when soil T > 0. 
		if (soillayer.wmass / soillayer.depth > Fwmass_min && soillayer.temp > 0.0) {

			// Gaseous CH4 volume in each layer [m3]
			double CH4_gas_vol = soillayer.CH4_gas / atomiccmass * R_gas * TsoilK / (hydro_press + atm_press);

			/// volumetric CH4 content [unitless]
			double CH4_vgc = CH4_gas_vol / (soillayer.depth / MM_PER_M * SQ_M); // m3/m3

			// Ebullition/bubble formation if the volumetric gas content (VGC) exceeds vgc_high * bubble_CH4_frac, 
			// where bubble_CH4_frac id the CH4 fraction of gas bubbles
			if (CH4_vgc > vgc_high * bubble_CH4_frac) {

				/// Volume of CH4 which bubbles out [m3]
				double CH4_ebull_vol = (CH4_vgc - vgc_low * bubble_CH4_frac) * soillayer.depth / MM_PER_M * SQ_M; // m3/m3

				// CH4 volumetric content after ebullition has occured
				CH4_vgc = vgc_low * bubble_CH4_frac;

				// CH4 which bubbles out[g CH4 - C layer - 1]
				soillayer.CH4_ebull_ind = CH4_ebull_vol * (hydro_press + atm_press) / (R_gas * TsoilK) * atomiccmass;

				// Gaseous CH4 concentration in each layer [g CH4-C layer-1 d-1]
				soillayer.CH4_gas = CH4_vgc * soillayer.depth / MM_PER_M * SQ_M * (hydro_press + atm_press) / (R_gas * TsoilK) * atomiccmass;

				// Update the amount of CH4 to emit today [g m-2 d-1], from this layer
				ebull_today += soillayer.CH4_ebull_ind;
			}
		}

		// Now recalculate the total CH4 amount in this layer, both dissolved and gaseous
		soillayer.CH4 = soillayer.CH4_diss + soillayer.CH4_gas;
	}

	// C mass at end
	double end_CH4 = get_sum(CH_4);

	if (fabs(initial_CH4 - (end_CH4 + ebull_today)) > 1.0e-12)
		dprintf("Year %d Day %d ERROR: gas_ebullition %g\n",date.year,date.day, initial_CH4 - (end_CH4 + ebull_today));
}


void Soilcolumn::methane(bool generatemethane) {

	// Main methane routine, called daily. 
	// Calculates daily methane fluxes and updates methane concentrations in peatland soil layers.
	// Algorithm etc. from Wania et al. 2010
	// This version coded by Paul Miller, based on Rita Wania's F90 code
	// Layer depths in m, and volume of water in m3

	Fluxes& fluxes = patch.fluxes;

	if (!generatemethane) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			// No generation of methane for what ever reason, so report the heterotrophic respiration 
			// flux calculated from the call to som_dynamics(), and 0 for the 4 CH4 fluxes.
			fluxes.report_flux(Fluxes::SOILC, soillayers[lyr].respsum);
			fluxes.report_flux(Fluxes::CH4C, 0.0);
			fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
			fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
			fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);
		}
	} 
	else if (patch.get_climate().lat < PEATLAND_WETLAND_LATITUDE_LIMIT) {
		
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			// Simple Spahni et al. (2011) approach in which a fraction (CH4toCO2_inundated) of the daily 
			// heterotrophic respiration is assumed to be in the form of CH4
			// g CH4-C/m2/day
			double inundated_CH4_flux_today = soillayer.respsum * G_PER_KG * CH4toCO2_inundated;

			// Report the heterotrophic respiration and the 4 CH4 fluxes.
			fluxes.report_flux(Fluxes::SOILC, soillayer.respsum - inundated_CH4_flux_today * KG_PER_G);
			fluxes.report_flux(Fluxes::CH4C, inundated_CH4_flux_today);
			fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
			fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
			fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);
		}
	} 
	else {

		// Layer diffusivities
		double D_CH4[NSOILLAYERS_MAX];
		double D_CO2[NSOILLAYERS_MAX];
		double D_O2[NSOILLAYERS_MAX];

		// Daily diffusion
		// mol O2 into the soil (and then diffused downwards). Should be negative, i.e. O2 diffuses INTO the soil
		double O2_diff_today = 0.0;
		double CH4_diff_today = 0.0;
	
		// Daily plant transport
		double O2_plant_today = 0.0;
		double CH4_plant_today = 0.0;
	
		// Total daily ebullition [g CH4-C d-1]
		double CH4_ebull_today = 0.0;

		// Tillers for C3-graminoids [m-2]
		double tillers;

		// Step 0 

		// C budget before today's methane calculations

		double ch4_budget_init = get_sum(CH_4) * KG_PER_G;
		double co2_budget_init = get_sum(RESPSUM);

		// Step 1

		// Calculate diffusivities and  gas constants

		// Set layer gas diffusivities (m2 d-1)
		calculate_gas_diffusivities(D_CH4, D_CO2, D_O2);

		// Update gas transport velocities and equilibrium gas concentrations at surface interface
		// Units: k_O2 etc: [m d-1], Ceq_O2: [mmol m-3]
		update_daily_gas_parameters();

		// Step 2

		// Calculate CH4 & CO2 production
		calculate_ch4_and_co2_production();

		// Step 3

		// Diffusion of O2
		diffuse_gas(D_O2, O2gas, Ceq_O2, k_O2, O2_diff_today);
		
		// Step 4

		// Plant transport of oxygen

		// Calculate the area and density of tillers
		calculate_tiller_areas(tillers);
	
		// O2 transport into soil through pants
		plant_gas_transport(Ceq_O2, k_O2, O2gas, tillers, O2_plant_today);
	
		// Step 5

		// Diffusion of CH4
		diffuse_gas(D_CH4, CH4gas, Ceq_CH4, k_CH4, CH4_diff_today);

		// Step 6
		
		// CH4 oxidation
		oxidation_ch4();

		// Step 7
		
		// Plant transport of CH4 into soil 
		plant_gas_transport(Ceq_CH4, k_CH4, CH4gas, tillers, CH4_plant_today);

		// Step 8

		// Ebullition of CH4
		// Reduces CH4 and updates CH4_ebull_today
		calculate_gas_ebullition(CH4_ebull_today);

		// Steps 9
		// Instead of calling diffusion and plant transport for CO2, I simply diffuse ALL the CO2. 
		// This should improve C conservation
		// In future one could look at how much IS emitted, especially as CO2 could be saved and emitted in 
		// bursts during spring thaw

		// Step 10

		// Calculate daily CH4 flux [gCH4-C m-2 d-1]
		double CH4_flux_today = CH4_diff_today + CH4_plant_today + CH4_ebull_today;

		// Calculate total heterotrophic respiration [kg C m-2]
		double tot_respsum = get_sum(CO_2) * KG_PER_G;

		// Report the heterotrophic respiration [kg C m-2] and the 4 CH4 fluxes [g C m-2].
		fluxes.report_flux(Fluxes::SOILC, tot_respsum);
		fluxes.report_flux(Fluxes::CH4C, CH4_flux_today);
		fluxes.report_flux(Fluxes::CH4C_DIFF, CH4_diff_today);
		fluxes.report_flux(Fluxes::CH4C_PLAN, CH4_plant_today);
		fluxes.report_flux(Fluxes::CH4C_EBUL, CH4_ebull_today);

		// All CO2 was respired today
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].CO2 = 0.0;
		}

		// Check C budget
		double ch4_budget_end = get_sum(CH_4) * KG_PER_G;
		double co2_budget_end = tot_respsum;
		double ch4_flux_budget_end = CH4_flux_today * KG_PER_G;

		if (fabs((ch4_budget_init + co2_budget_init) - (ch4_budget_end + ch4_flux_budget_end + co2_budget_end)) > 1.0e-12)
			dprintf("Year %d Day %d ERROR: Methane %g INIT CH4 %g CO2 %g END CH4 %g CH4_flux %g CO2 %g\n",date.year,date.day, 
				(ch4_budget_init + co2_budget_init) - (ch4_budget_end + ch4_flux_budget_end + co2_budget_end),
				ch4_budget_init, co2_budget_init, ch4_budget_end, ch4_flux_budget_end, co2_budget_end);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Wania, R., Ross, I., and Prentice, I. C.: Integrating peatlands and permafrost into a dynamic
//   global vegetation model; 1, Evaluation and sensitivity of physical land surface processes, Global
//   Biogeochemical Cycles, 23, doi:10.1029/2008gb003412, 2009a.
// Wania, R., Ross, I., and Prentice, I. C.: Integrating peatlands and permafrost into a dynamic
//   global vegetation model; 2, Evaluation and sensitivity of vegetation and carbon cycle processes,
//   Global Biogeochemical Cycles, 23, doi:10.1029/2008gb003413, 2009b.
// Wania, R., Ross, I., and Prentice, I. C.: Implementation and evaluation of a new methane model
//   within a dynamic global vegetation model: LPJ-WHyMe v1.3.1, Geosci Model Dev, 3, 565-584,
//   doi:DOI 10.5194/gmd-3-565-2010, 2010. 
