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

/// Crank-Nicholson timestepper algorithm for gas diffusion equation.
void cnstepgas(std::vector<double> Di, std::vector<double> dz, double surf_conc, double dt, 
			std::vector<double> conc, int nsoillayers) {

	// UNITS:
	// Di			m2 d-1
	// dz			m
	// surf_conc	mmol m-3
	// dt			d
	// conc			mmol m-3

	// A fill-in for layers above layer0.
	double MISSING_VALUE;

	// Layer counters: note that there are two different layer counting
	// schemes used, one for the input and output parameters (vectors of
	// length NLAYERS) and one for the values used in the Crank-Nicholson
	// solver (vectors of length active_layers)

	int layer, lidx;
	int active_layers;
	double dplus;
	double dminus;
	double dz_factor;
	double Cplus; 
	double Cp_minus;
	double dzhere, dzminus, dzplus;
	double cohere, cominus, coplus;

	// Leading diagonal, left and right subdiagonals for Crank-Nicholson
	// matrix.
	std::vector<double> diag(nsoillayers);
	std::vector<double> left(nsoillayers);
	std::vector<double> right(nsoillayers);

	// Right hand side vector for Crank-Nicholson scheme equations.
	std::vector<double> rhs(nsoillayers);

	// Solution vector for Crank-Nicholson scheme equations.
	std::vector<double> solution(nsoillayers);

	// initialise
	dplus = dminus = dz_factor = Cplus = Cp_minus = 0.0;
	dzhere = dzminus = dzplus = cohere = cominus = coplus = 0.0;

    // --- CODE STARTS HERE ---

	MISSING_VALUE = -9999.0;

	active_layers = nsoillayers;

	//   BUILD TRIDIAGONAL MATRIX AND KNOWN RIGHT HAND SIDE

	// End members for off-diagonal elements.
	left[0] = 0.0;
	right[active_layers - 1] = 0.0;

	// Process the active layers. 
	// The first time (lidx and layer = layer0) corresponds to the surface layer

	for (lidx = 1; lidx <= active_layers; lidx++) {

		// Deal with different layer counting schemes.
		layer = lidx - 1; 
		// Minimum is layer0
		// Maximum is NLAYERS - 1

		// Calculate diffusion constants averaged over adjacent layers.
		// The diffusion constant at the bottom layer is clamped to zero
		// to enforce the no heat flow boundary condition there.  

		// D+

		if (layer == nsoillayers - 1) 
			dplus = 0.0; // BC2 - Bottom layer diffusion clamped to 0
		else 
			dplus = 0.5 * (Di[layer] + Di[layer + 1]);

		// D-

		if (layer == 0) 
			dminus = Di[layer]; // top layer
		else
			dminus = 0.5 * (Di[layer] + Di[layer - 1]); // soil layers	
 
		// --- HERE ---
		if (layer < nsoillayers) {
			// all soil layers
			dzhere = dz[layer];
			cohere = conc[layer];
		}

		// --- PLUS ---
		if (layer < nsoillayers - 1) {
			// all soil layers apart from the bottom soil layer
			dzplus = dz[layer + 1];
			coplus = conc[layer + 1];
		}

		// --- MINUS ---
		if (layer == 0) {
			// top layer
			dzminus = dz[layer];
			cominus = conc[layer];
		}
		else if (layer < nsoillayers) {
			// all soil layers and the first padding layer
			dzminus = dz[layer - 1];
			cominus = conc[layer - 1];
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
				right[lidx-1] = 0.0; // i.e. BC1, where top layer == (1,0,..,0)

			// left[0] is set above.
			if (lidx > 1)
				left[lidx - 1] = -Cp_minus;			
		}

		if (lidx == active_layers) 
			left[lidx - 1] = -Cp_minus;
		 
		// RHS
		// Calculate right hand side vector values.
		if (lidx == 1) 
			rhs[0] = surf_conc;
		else if (lidx == active_layers) // Cplus == 0 here anyway
			rhs[lidx - 1] = (1.0 - Cp_minus) * cohere + Cp_minus * cominus;
		else
			rhs[lidx - 1] = (1.0 - Cplus - Cp_minus) * cohere +
			Cplus * coplus + Cp_minus * cominus;
	} // end for

	//   SOLVE TRIDIAGONAL SYSTEM
	tridiag(active_layers, left, diag, right, rhs, solution);

	//   FORMAT OUTPUT PARAMETERS

	// Transfer the solution to the concentration array.
	for (lidx = 0; lidx < nsoillayers; lidx++) {

		conc[lidx] = solution[lidx];
	}
}


bool Soilcolumn::calculate_gas_diffusivities(std::vector<double>& dCH4, std::vector<double>& dCO2, std::vector<double>& dO2) {

	// Calculate and return diffusivities, in units of m2 d-1
	// Called each day
	// See Wania et al. (2010) - Sec 2.5

	double D_CH4_water;
	double D_CO2_water;
	double D_O2_water;

	double D_CH4_air;
	double D_CO2_air;
	double D_O2_air;

	const double scale = 0.00001;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		dCH4[lyr] = 0.0;
		dCO2[lyr] = 0.0;
		dO2[lyr] = 0.0;

		// AIR - UNITS: 10-4 m2 s-1 == cm2 s-1
		// Lerman (1979)

		double layerT = soillayers[lyr].temp;

		// Wania et al. (2010), Eqn. 10
		D_CH4_air = 0.1875 + 0.00130 * layerT;
		D_CO2_air = 0.1325 + 0.00090 * layerT;
		D_O2_air = 0.1759 + 0.00117 * layerT;

		// WATER - UNITS: 10-4 m2 s-1 == cm2 s-1
		// Wania et al. (2010), Eqn. 9 - where the diffusivities are in units of 10**-9 m2 s-1. 
		// * by scale to get 10-4 m2 s-1
		D_CH4_water = (0.9798 + 0.02986 * layerT + 0.0004381 * layerT * layerT) * scale;
		D_CO2_water = (0.9390 + 0.02671 * layerT + 0.0004095 * layerT * layerT) * scale;
		D_O2_water = (1.1720 + 0.03443 * layerT + 0.0005048 * layerT * layerT) * scale;

		// Fraction air in soil layer
		double frac_air = soillayers[lyr].fracair();
	
		if (frac_air > 0.05 && soillayers[lyr].acrotelm) {
			// Wania et al. (2010), Eqn. 11 
			double airpow = pow(frac_air, 10 / 3) / pow(acrotelm_por, 2.0);
			dCH4[lyr] = airpow * D_CH4_air;
			dCO2[lyr] = airpow * D_CO2_air;
			dO2[lyr] = airpow * D_O2_air;
		}
		else {
			// Water diffusivities dominate
			dCH4[lyr] = D_CH4_water;
			dCO2[lyr] = D_CO2_water;
			dO2[lyr] = D_O2_water;
		}
		
		// Convert from cm2 s-1 to m2 d-1
		dCH4[lyr] *= SECS_PER_DAY / CM2_PER_M2;
		dCO2[lyr] *= SECS_PER_DAY / CM2_PER_M2;
		dO2[lyr] *= SECS_PER_DAY / CM2_PER_M2;
	}

	return true;
}


bool Soilcolumn::update_daily_gas_parameters() {

	// Update gas transport velocities [m d-1] and equilibrium gas concentrations [mmol m-3]
	// for CH4, CO2 and O2
	// Called daily by Soil::methane()
	// See Wania et al. (2010) - Sec 2.5

	double surfT = soillayers[0].temp;

	// Gas transport velocities

	// gas transport velocity of SF6.
	// Wania et al. (2010) - Eqn. 6 (cm h-1)
	double k_600 = 2.07 + 0.215 * pow(patch.get_climate().u10, 1.7);

	// Schmidt number of O2
	// Wania et al. (2010) - Eqn. 7
	double ScO2 = 1800.6 - 120.1 * surfT + 3.7818 * pow(surfT, 2) - 0.047608 * pow(surfT, 3);

	// gas transport velocity of O2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// converted from cm h-1 to m d-1
	k_O2 = k_600 * pow(ScO2 / 600.0, n_coeff);
	k_O2 *= 24.0 * M_PER_CM;

	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCH4 = 1898.0 - 110.1 * surfT + 2.834 * pow(surfT, 2) - 0.02791 * pow(surfT, 3);

	// gas transport velocity of CH4 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// converted from cm h-1 to m d-1
	k_CH4 = k_600 * pow(ScCH4 / 600.0, n_coeff);
	k_CH4 *= 24.0 * M_PER_CM;

	// Schmidt number of CH4
	// Wania et al. (2010) - Eqn. 7
	double ScCO2 = 1911.0 - 113.7 * surfT + 2.967 * pow(surfT, 2) - 0.02943 * pow(surfT, 3);

	// gas transport velocity of CO2 [cm h-1]
	// Wania et al. (2010) - Eqn. 5
	// converted from cm h-1 to m d-1
	k_CO2 = k_600 * pow(ScCO2 / 600.0, n_coeff);
	k_CO2 *= 24.0 * M_PER_CM;

	// Equilibrium gas concentration
	// See Wania et al. (2010) - Eqn. 8

	// 25 degrees C
	double deg25 = K2degC + 25.0;

	// [L atm mol-1] DW_COMMENT description
	double henry_coeff_O2 = henry_k_O2 * exp(-1.0 * henry_C_O2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for O2
	// pp_gas/MM2_PER_M2 (mol L-1) converts to atm units DW_COMMENT description
	// mol L-1 to mmol m-3
	Ceq_O2 = pp_O2 / MM2_PER_M2 / henry_coeff_O2;
	Ceq_O2 *= MM2_PER_M2; 

	// [L atm mol-1] DW_COMMENT description
	double henry_coeff_CO2 = henry_k_CO2 * exp(-1.0 * henry_C_CO2 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for CO2
	// Use this gridcell's CO2 concentration, not a fixed value as in Wania et al. (2010).
	// mol L-1 to mmol m-3
	Ceq_CO2 = patch.get_climate().co2 / MM2_PER_M2 / henry_coeff_CO2;
	Ceq_CO2 *= MM2_PER_M2;

	// [L atm mol-1] DW_COMMENT description
	double henry_coeff_CH4 = henry_k_CH4 * exp(-1.0 * henry_C_CH4 * (1.0 / (surfT + K2degC) - 1.0 / deg25));

	// Equilibrium gas concentrations (mmol m - 3) for CH4
	// mol L-1 to mmol m-3
	Ceq_CH4 = pp_CH4 / MM2_PER_M2 / henry_coeff_CH4;
	Ceq_CH4 *= MM2_PER_M2;

	return true;
}


double Soilcolumn::diffuse_gas(std::vector<double> D, gastype thisgastype, double Ceq,
					  double kgas, std::vector<double> Dz_m, double &dailyDiff) {

	// Generic gas diffusion method that works with CO2, CH4 and O2 (as specified with gastype) 
	// See Wania et al. (2010) - Sec 2.5 - for a full description

	// Called like this (e.g. for O2): 
	// diffuse_gas(O2, D_O2, O2gas, Ceq_O2, k_O2, Dz_metre, dailyO2diffusion);
	
	// Atomic mass of this gas [gC/mol]
	double atomic_mass;

	// Concentration of the dissolved gas in question [mmol m-3]
	std::vector<double> C(nsoillayers);

	// layer water + ice volume [m3]
	std::vector<double> total_volume_water(nsoillayers);
	
	if (thisgastype != O2gas)
		atomic_mass = atomiccmass; // i.e. CO2 or CH4 - 12 gC/mol
	else
		atomic_mass = 1.0; // O2 already in mol layer-1
	
	double initialAmount = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		total_volume_water[lyr] = (soillayers[lyr].wmass + soillayers[lyr].imass) * M_PER_MM;

		if (thisgastype == CO2gas) {
			initialAmount += soillayers[lyr].CO2_soil;
			// g layer - 1 to mmol m - 3
			C[lyr] = soillayers[lyr].CO2_soil / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
		}
		else if (thisgastype == CH4gas) {
			initialAmount += soillayers[lyr].CH4;
			// g layer - 1 to mmol m - 3
			C[lyr] = soillayers[lyr].CH4 / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
		}
		else {
			initialAmount += soillayers[lyr].O2;
			// mol layer-1 to mmol m-3
			C[lyr] = soillayers[lyr].O2 / atomic_mass / total_volume_water[lyr] * MMOL_PER_MOL;
		}
	}

	// Set the BC, i.e. equilibrium gas concentrations in the top layer depending on atmospheric concentrations
	// and using Henry's law. See Wania et al. (2010), Eqs. 4-8
	
	// New surface concentration [mmol m-3]
	double Cnew;
	if (thisgastype == O2gas || thisgastype == CH4gas) { // Could also run for CH4 here below.

		// No diffusion if there's too little liquid water in the top layer
		// Could add a snow restriction
		if ((soillayers[0].wmass / soillayers[0]. depth) < Fwmass_min) {
			// Unchanged surface concentration
			Cnew = C[0]; 
			dailyDiff = 0.0;
		} 
		else {
			double volume_liquid_water = soillayers[0].wmass * M_PER_MM;

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
			if ((soillayers[lyr].wmass / soillayers[lyr].depth) < Fwmass_min) {
				D[lyr] = 1e-9; // [m2 d-1]
			}
		}

		std::vector<double> C_init(nsoillayers);
		std::vector<double> C_last(nsoillayers);

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

		if (thisgastype == CO2gas) {
			soillayers[lyr].CO2_soil = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayers[lyr].CO2_soil;
		}
		else if (thisgastype == CH4gas) {
			soillayers[lyr].CH4 = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayers[lyr].CH4;
		}
		else {
			soillayers[lyr].O2 = C[lyr] * atomic_mass * total_volume_water[lyr] / MMOL_PER_MOL;
			finalAmount += soillayers[lyr].O2;
		}
	}

	// Return the amount of gas that has diffused INTO the soil
	// CO2 & CH4 - g
	// O2 - mol
	double gasIntoSoil = finalAmount - initialAmount;

	if (thisgastype == CH4gas)
		dailyDiff = gasIntoSoil; // gC/m2/d

	return gasIntoSoil;
}


double Soilcolumn::calculate_tiller_areas() {

	// Calculate the area of tillers
	// See Wania et al. (2010) - Sec 2.6

	double tiller_density = 0.0; // to return [number m-2]

	double graminoid_leafcmass = 0.0; // kgC/m2
	double graminoid_dphen = 0.0;
	double graminoid_anpp_red = 0.0;

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
		}
		vegetation.nextobj();
	}

	graminoid_leafcmass *= G_PER_KG; // kgC m-2 to gC m-2

	tiller_density = graminoid_leafcmass * graminoid_dphen / tiller_weight; // m-2

	double tiller_frac = 0.0;

	std::vector<double> rootfrac = get_peatland_rootfrac(); // DW_COMMENT saved in soilcolumn in trunk and only updated on 1st Jan

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		if (!negligible(tiller_density)) {
			tiller_frac = tiller_density * rootfrac[lyr];

			// Consider this as the fraction of the area (1m2) occupied by tillers in this later
			soillayers[lyr].tiller_area = PI * pow(tiller_radius, 2.0) * tiller_frac * tiller_por; // [unitless m2 m-2]

			if (soillayers[lyr].tiller_area < 0.000000001)
				soillayers[lyr].tiller_area = 0.0;
		}
		else {
			soillayers[lyr].tiller_area = 0.0;
		}
	}

	return tiller_density;
}


bool Soilcolumn::plant_gas_transport(double Ceq, double kgas, gastype thisgastype, 
							 double& plantTransportToday) {

	// Plant transport of O2 or CH4
	// See Wania et al. (2010) - Sec 2.6

	// The limiting factor to plant transport are the number of tillers as
    // they represent the cross-sectional area available to gas transport.
	// The biomass is related to tiller density based on Schimel 1995. The
	// average biomass was 185g/m2 and the average number of tillers was 380/m2.
    // This gives us 2.05 tillers per g biomass and 0.48 g per tiller, which
    // corresponds to 0.22 g C per tiller

	std::vector<double> Cgas(nsoillayers);
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

		// Weight plant transport by the area of porous root cross-sections. 
		// plant_trans(i) [mmol layer-1 d-1]
        // the 2.07 are from the equation for k_CH4, when U10 is zero.

		double total_volume_water = (soillayers[lyr].wmass + soillayers[lyr].imass) * M_PER_MM;
		double volume_liquid_water = soillayers[lyr].wmass * M_PER_MM;

		// CH4 and CO2: gC layer-1 to mmol m-3
		// O2: mol layer-1 to mmol m-3
		if (thisgastype == O2gas)
			Cgas[lyr] = soillayers[lyr].O2 / atomic_mass / total_volume_water * MMOL_PER_MOL;
		else
			Cgas[lyr] = soillayers[lyr].CH4 / atomic_mass / total_volume_water * MMOL_PER_MOL;

		if (soillayers[lyr].wmass / soillayers[lyr].depth < Fwmass_min || negligible(soillayers[lyr].tiller_area) || volume_liquid_water < 0.0001) {
			Cnew = Cgas[lyr]; // no plant transport to/from this layer
		} 
		else {
			double exponent = kgas / volume_liquid_water;
			exponent *= soillayers[lyr].tiller_area;
			double fac = exp(-1.0 * exponent);
			double change = Cnew + (Cgas[lyr] - Ceq) * fac;

			Cnew = Ceq + (Cgas[lyr] - Ceq) * exp(-kgas / (volume_liquid_water / soillayers[lyr].tiller_area)); // mmol m-3
			plantTransportToday += (Cgas[lyr] - Cnew) 
				* atomic_mass * total_volume_water / MMOL_PER_MOL; // mol layer-1 d-1 (O2) or // gC m-2 d-1 (CH4 & CO2)

			// Restrict TINY negative plant transport
			if (plantTransportToday < 0.0 && fabs(plantTransportToday) < 1e-8) {
				plantTransportToday = 0.0;
				Cnew = Cgas[lyr];
			}
		}

		// mmol m-3 to gC layer-1 (CO2 or CH4), or mol layer-1 (O2) 
		if (thisgastype != O2gas)
			soillayers[lyr].CH4 = Cnew * atomic_mass * total_volume_water / MMOL_PER_MOL;
		else
			soillayers[lyr].O2 = Cnew * atomic_mass * total_volume_water / MMOL_PER_MOL;
	}
	return true;
}


bool Soilcolumn::calculate_gas_ebullition(double& ebull_today) {

	// Gas ebullition
	// Returns the gas ebullition today
	// See Wania et al. (2010) - Sec 2.7

	double waterheight = 0.0;
	ebull_today = 0.0;

	// Today's water table position (-300 to +100) mm
	double wtpp = wtp[date.day];

	// Whether to emit CH4 directly to the atmosphere (true by default)
	bool emitToAtmosphere = true;

	// The soil layer to bubble to.
	int bubbleToLayer = 0;

	if (wtpp <= 0) {
		if (wtpp >= -100.0)
			bubbleToLayer = 0;
		else if (wtpp >= -200.0 && wtpp < -100.0)
			bubbleToLayer = 1;
		else
			bubbleToLayer = 2;
	}
	else { // Standing water
		bubbleToLayer = -1;
	}

	// Ensure ebullition from all layers
	bubbleToLayer = -1;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[lyr];

		double TsoilK = soillayer.temp + K2degC;

		double total_volume_water = (soillayer.wmass + soillayer.imass) * M_PER_MM;
		double volume_liquid_water = soillayer.wmass * M_PER_MM;

		waterheight += total_volume_water; // m

		// Bubble formation below the WT only
		if (lyr > bubbleToLayer) {

			// Max CH4 that can be dissolved (Wania et al. (2010), Eqn 15) [ml CH4 ml-1 H2O]
			double CH4_diss_max = 0.05708 - 0.001545 * max(0.0, soillayer.temp) +
				0.00002069 * max(0.0, pow(soillayer.temp,2.0));

			// Water pressure [N, or kg m s-2]
			double hydro_press = rho_H2O * gravity * waterheight;

			// Maximum volume that can be dissolved [m3 CH4 layer-1}
			double CH4_diss_max_m3 = CH4_diss_max * volume_liquid_water;

			// Use ideal gas law to convert to mol CH4 layer-1
			double CH4_diss_max_mol = CH4_diss_max_m3 * (hydro_press + atm_press) / (R_gas * TsoilK);

			// Convert to g CH4-C layer-1
			double CH4_diss_max_g = CH4_diss_max_mol * mr_C;
			
			soillayer.CH4_ebull_ind = 0.0;
			
			// Restrict ebullition to cases when there is enough liquid water and when soil T > 0. 
			if (soillayer.wmass / soillayer.depth > Fwmass_min && soillayer.temp > 0.0) {

				double henry_k_cc_CH4 = TsoilK / (12.2 * henry_k_CH4);
				soillayer.CH4_diss = min(CH4_diss_max_g, henry_k_cc_CH4 * soillayer.CH4);
				
				// Override (a la Wania et al. 2010) with this new, simpler ebullition
				soillayer.CH4_diss = min(CH4_diss_max_g, soillayer.CH4);

				soillayer.CH4_gas = soillayer.CH4 - soillayer.CH4_diss; // g CH4-C layer-1 - gas
				soillayer.CH4_gas_vol = soillayer.CH4_gas / atomiccmass * R_gas * TsoilK / (hydro_press + atm_press); // [m3]
				soillayer.CH4_vgc = soillayer.CH4_gas_vol / (soillayer.depth / MM_PER_M * SQ_M); // m3/m3
			
				// Ebullition/bubble formation if the volumetric gas content (VGC) exceeds vgc_high * bubble_CH4_frac, 
				// where bubble_CH4_frac id the CH4 fraction of gas bubbles
				if (soillayer.CH4_vgc > vgc_high * bubble_CH4_frac) {
				
					soillayer.CH4_ebull_vol = (soillayer.CH4_vgc - vgc_low * bubble_CH4_frac) * soillayer.depth / MM_PER_M * SQ_M; // m3/m3
					soillayer.CH4_vgc = vgc_low * bubble_CH4_frac;
					soillayer.CH4_ebull_ind = soillayer.CH4_ebull_vol * (hydro_press + atm_press) / (R_gas * TsoilK); // mol
					soillayer.CH4_gas = soillayer.CH4_vgc * soillayer.depth / MM_PER_M * SQ_M * (hydro_press + atm_press) / (R_gas * TsoilK); // mol
					soillayer.CH4_ebull_ind *= atomiccmass; // mol to g
					soillayer.CH4_gas *= atomiccmass; // mol to g 

					// Update the amount of CH4 to emit today, from this layer
					ebull_today += soillayer.CH4_ebull_ind; // g m-2 d-1
				}
				// Now recalculate the total CH4 amount in this layer, both dissolved and gaseous
				soillayer.CH4 = soillayer.CH4_diss + soillayer.CH4_gas;
			} 
			else {
				// No change in the dissolved CH4 amount
				soillayer.CH4_diss = soillayer.CH4_diss_yesterday;
			}
		} 
		else {
			// No change in the dissolved CH4 amount
			soillayer.CH4_diss = soillayer.CH4_diss_yesterday;
		}
	} 

	// If the water table is below the surface, then put the bubbled methane into the first unsaturated layer
	if (!emitToAtmosphere) {

		Soillayer& bubblelayer = soillayers[bubbleToLayer];

		bubblelayer.CH4 += ebull_today;
		ebull_today = 0.0; // No emission today.

		// O2 is in [mol layer-1]
		bubblelayer.O2 *= oxid_frac; // since 25% of O2 used by roots themselves...  
		bubblelayer.CH4 /= atomiccmass; // mol layer-1
		bubblelayer.CH4_oxid = min(bubblelayer.CH4, 0.5 * bubblelayer.O2); // usually 75%
		bubblelayer.CH4 = (bubblelayer.CH4 - bubblelayer.CH4_oxid) * atomiccmass; // gC layer again	
		
		// Subtract the moles used in oxidation
		bubblelayer.O2 -= 2.0 * bubblelayer.CH4_oxid;

		// Oxidised CH4 becomes CO2 (and water!)
		bubblelayer.CO2_soil += bubblelayer.CH4_oxid * atomiccmass;
	}

	return true;
}

void Soilcolumn::calculate_carbon_store(int dy, bool today) {

	// Carbon accounting routine - updates co2_store and ch4_store 
	// See Wania et al. (2010)

	int day = 0;

	// Reset the stores
	co2_store = 0.0;
	ch4_store = 0.0;

	for (int lyr = 0; lyr < nsoillayers; lyr++) {

		Soillayer& soillayer = soillayers[0];

		if (today)
			day = dy;
		else
			day = dy - 1;

		if (day > 0 || today) {
			co2_store += soillayer.CO2_soil;
			ch4_store += soillayer.CH4;
		}
		else {
			// Jan 1 or yesterday
			co2_store += soillayer.CO2_soil_yesterday;
			ch4_store += soillayer.CH4_yesterday;
		}
	}
}


bool Soilcolumn::methane(bool generatemethane) {

	// Main methane routine, called daily. 
	// Calculates daily methane fluxes and updates methane concentrations in peatland soil layers.
	// Algorithm etc. from Wania et al. 2010
	// This version coded by Paul Miller, based on Rita Wania's F90 code
	// Layer depths in m, and volume of water in m3

	if (!generatemethane) {

		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			// No generation of methane for what ever reason, so report the heterotrophic respiration 
			// flux calculated from the call to som_dynamics(), and 0 for the 4 CH4 fluxes.
			patch.fluxes.report_flux(Fluxes::SOILC, soillayers[lyr].dcflux_soil); // kgC/m2
			patch.fluxes.report_flux(Fluxes::CH4C, 0.0);
			patch.fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
			patch.fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
			patch.fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);
		}
	} 
	else if (patch.get_climate().lat < PEATLAND_WETLAND_LATITUDE_LIMIT) {
		
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			// Simple Spahni et al. (2011) approach in which a fraction (CH4toCO2_inundated) of the daily 
			// heterotrophic respiration is assumed to be in the form of CH4

			double inundated_CH4_flux_today = soillayers[lyr].dcflux_soil * G_PER_KG * CH4toCO2_inundated; // g CH4-C/m2/day

			// Report the heterotrophic respiration and the 4 CH4 fluxes.
			patch.fluxes.report_flux(Fluxes::SOILC, soillayers[lyr].dcflux_soil - inundated_CH4_flux_today * KG_PER_G); // kgC/m2
			patch.fluxes.report_flux(Fluxes::CH4C, inundated_CH4_flux_today); // g CH4-C/m2/day
			patch.fluxes.report_flux(Fluxes::CH4C_DIFF, 0.0);
			patch.fluxes.report_flux(Fluxes::CH4C_PLAN, 0.0);
			patch.fluxes.report_flux(Fluxes::CH4C_EBUL, 0.0);
		}
	} 
	else {

		if (date.day == 200)
			int sch = 0;


        // Max allowed error in checks
        const double MAX_ERR = 0.000001;
		const double MAX_ERR_BALANCE = 0.0001;
		const double LARGE_ERR = 0.01;

		// Variables for debugging
		bool allow_planttransport = true;
		bool allow_ebullition = true;
		bool allow_ch4diffusion = true;
		bool allow_o2diffusion = true;

		// Layer diffusivities
		std::vector<double> D_CH4(nsoillayers);
		std::vector<double> D_CO2(nsoillayers);
		std::vector<double> D_O2(nsoillayers);

		std::vector<double> Dz_metre(nsoillayers); // layer depths [m]

		// Daily diffusion 
		double O2_diff_today;
		double CO2_diff_today;
		double CH4_diff_today;
	
		// Daily plant transport
		double O2_plant_today;
		double CH4_plant_today;
		double CO2_plant_today;
	
		// Total daily ebullition [g CH4-C d-1]
		double CH4_ebull_today;
		
		int daynum = date.day;

		// Budget/conservation variables [gC/m2]
		double ch4_c_in = 0.0;
		double co2_c_in = 0.0;

		double ch4_c_store_init;
		double co2_c_store_init;

		double ch4_c_store_now;
		double co2_c_store_now;

		// Initialise components of total_C_flux, CO2_flux and CH4_flux
		O2_diff_today = 0.0;
		CO2_diff_today = 0.0;
		CH4_diff_today = 0.0;
		CH4_ebull_today = 0.0;
		CH4_plant_today = 0.0;
		CO2_plant_today = 0.0;
		O2_plant_today = 0.0;

		// Step 0 

		// C budget before today's methane calculations

		calculate_carbon_store(daynum, false);
		ch4_c_store_init = ch4_store;
		co2_c_store_init = co2_store;

		// Step 1

		// Diffusivities and layer depths, and update of gas constants

		// Set layer gas diffusivities (m2 d-1)
		calculate_gas_diffusivities(D_CH4, D_CO2, D_O2);

		std::vector<double> volume_liquid_water(nsoillayers);
		std::vector<double> total_volume_water(nsoillayers);

		// Layer water and ice content (m3)
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			volume_liquid_water[lyr] = soillayers[lyr].wmass * M_PER_MM;
			total_volume_water[lyr] = (soillayers[lyr].wmass + soillayers[lyr].imass) * M_PER_MM; 
		}

		// If we have standing water add it to the top soil layer
		if (wtp[daynum] > 0.0) {
			volume_liquid_water[0] += wtp[daynum] * M_PER_MM;
			total_volume_water[0] += wtp[daynum] * M_PER_MM;
		}


		// Update the temperature-dependent gas parameters
		// Units: 
		// k_O2 etc: [m d-1]
		// Ceq_O2: [mmol m-3]
		bool gasParamsOK = update_daily_gas_parameters();

		if (!gasParamsOK)
			return false;

		// Step 2

		// Calculate CH4 & CO2 production

		// Conservation check - C in methane produced
		double c_input = 0.0;
		double total_ch4 = 0.0;
		double total_drh = 0.0;

		// Loop through soil layers
		for (int lyr = 0; lyr < nsoillayers; lyr++) {

			Soillayer& soillayer = soillayers[lyr];

			// Daily decomposition, converted to gC/m2 from kgC/m2
			double drh = soillayer.dcflux_soil * G_PER_KG;
			total_drh += drh;

			double anoxic = 1.0 - soillayer.fracair() - Fgas;

			// CH4 production in the layer (gC/m2)
			if (soillayer.wmass / soillayer.depth < Fwmass_min)
				soillayer.CH4_prod = 0.0;
			else
				soillayer.CH4_prod = anoxic * CH4toCO2_peat * drh;

			c_input += soillayer.CH4_prod;

			// Add CH4 production to CH4 pool 
			soillayer.CH4 = soillayer.CH4_yesterday + soillayer.CH4_prod;
			total_ch4 += soillayer.CH4;

			// CO2 production
			soillayer.CO2_soil_prod = drh - soillayer.CH4_prod;
			
			c_input += soillayer.CO2_soil_prod;

			// CO2 pool
			// Add CO2 production to CO2 pool 
			soillayer.CO2_soil = soillayer.CO2_soil_yesterday + soillayer.CO2_soil_prod;
		}

		// Debugging:
		if ((total_ch4 < 0.0000000 || total_ch4 > 10000000) && verbosity >= WARNING) {
			dprintf("%s%8.5f\n","Bad total_ch4 at the start of Soil::methane()",total_ch4);	
		}

		// C conservation test:
		calculate_carbon_store(daynum, true);
		ch4_c_store_now = ch4_store;
		co2_c_store_now = co2_store;
		double check = ch4_c_store_now + co2_c_store_now - ch4_c_store_init - co2_c_store_init - total_drh;
		if (fabs(check) > MAX_ERR && verbosity >= WARNING)
			dprintf("C imbalance in Soil::methane() : After C in: %g\n",check);

		// Step 3

		// Diffusion of O2

		double dailyO2diffusion = 0.0;
		double molesO2IntoSoil = 0.0; 

		// No O2 diffusion until snow depth < 50mm 
		if (allow_o2diffusion && snow.depth() < 50.0) {
			molesO2IntoSoil = diffuse_gas(D_O2, O2gas, Ceq_O2, k_O2, Dz_metre, dailyO2diffusion);
		}

		// mol O2 into the soil (and then diffused downwards)
		// Should be negative, i.e. O2 diffuses INTO the soil
		O2_diff_today = dailyO2diffusion;

		// Step 4

		// Plant transport of oxygen

		// Tiller set-up
		// Loop through the individuals present and sum the leaf carbon mass for C3-graminoids [m-2] 
		double tillers = calculate_tiller_areas();
	
		// O2 plant transport
		double plantTransportToday_O2 = 0.0;
	
		// Only allow plant transport if there are graminoid tillers present
		if (allow_planttransport && tillers > 0)
			plant_gas_transport(Ceq_O2, k_O2, O2gas, plantTransportToday_O2);

		O2_plant_today = plantTransportToday_O2; // mol O2 into soil through plants
	
		// Step 5

		// Diffusion of CH4

		double ch4_before_diffusion = ch4_store;

		double dailyCH4diffusion = 0.0;
		double gramCH4IntoSoil = 0.0; 
		// No CH4 diffusion until snow depth < 50mm 
		if (allow_ch4diffusion && snow.depth() < 50.0) {
			gramCH4IntoSoil = diffuse_gas(D_CH4, CH4gas, Ceq_CH4, k_CH4, Dz_metre, dailyCH4diffusion);
		}
		// Should be positive, i.e. upward flux [gC m-2 d-1] - should be ZERO?
		CH4_diff_today = -dailyCH4diffusion;

		// Remove tiny values
		if (fabs(CH4_diff_today) < MAX_ERR)
			CH4_diff_today = 0.0;

		// C conservation test:
		if ((CH4_diff_today < -0.1 || CH4_diff_today > 10000000) && verbosity >= WARNING) {
			dprintf("%s%8.5f\n","Bad CH4 diffusion in Soil::methane()",CH4_diff_today);	
		}

		// C conserved?
		calculate_carbon_store(daynum, true);
		check = ch4_store -ch4_before_diffusion+CH4_diff_today;
		if (fabs(check) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after diffusion of CH4: ",check);

		// Step 6
		
		// CH4 oxidation

		// CH4 units: gC layer-1 
		// O2 units: mol layer-1

		// Assume that 1/2 of the O2 is utilized by other electron
		// acceptors and 2 moles of O2 is needed to oxidise 1 mole of CH4.
		// O2 is in [mol layer-1]
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			// Since ((1-oxid_frac)*100)% of O2 used by roots themselves... [mol layer-1] 
			soillayers[lyr].O2 *= oxid_frac;
			soillayers[lyr].CH4 /= atomiccmass; // mol layer-1
			soillayers[lyr].CH4_oxid = min(soillayers[lyr].CH4, 0.5 * soillayers[lyr].O2); // usually 75%
			soillayers[lyr].CH4 = (soillayers[lyr].CH4 - soillayers[lyr].CH4_oxid) * atomiccmass; // gC layer again

			// Subtract the moles used in oxidation
			soillayers[lyr].O2 -= 2.0 * soillayers[lyr].CH4_oxid;

			// Oxidised CH4 becomes CO2 (and water!)
			soillayers[lyr].CO2_soil += soillayers[lyr].CH4_oxid * atomiccmass;
		}

		// Step 6b

		// CH4 diffusion from top layer AFTER diffusion and oxidation
		// Now done after CN diffusion

		double ch4_diff_before_BC = CH4_diff_today;

		// Check C budget
		double ch4_before_oxid_and_diffusion = ch4_store;
		double co2_before_oxid_and_diffusion = co2_store;
	
		double checkbefore = ch4_before_oxid_and_diffusion + co2_before_oxid_and_diffusion + ch4_diff_before_BC;

		calculate_carbon_store(daynum, true);
		check = ch4_store + co2_store + CH4_diff_today;

		double checkafter1 = check - checkbefore;

		if (fabs(checkafter1) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after oxidation and diffusion of CH4: ", checkafter1);

		checkbefore = check;

		// Step 7
		
		// Plant transport of CH4

		double plantTransportToday_CH4 = 0.0;

		// Only allow plant transport of there are gramoinoid tillers present
		if (allow_planttransport && tillers > 0)
			plant_gas_transport(Ceq_CH4, k_CH4, CH4gas, plantTransportToday_CH4);
		
		// CH4 into soil through plants (gC m-2)
		CH4_plant_today = plantTransportToday_CH4;

		// Check C budget
		calculate_carbon_store(daynum, true);
		check = ch4_store + co2_store + CH4_diff_today + CH4_plant_today;
		double checkafter2 = check - checkbefore;
		if (fabs(checkafter2) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after plant transport of CH4: ",checkafter2);
		checkbefore = check;

		// Step 8

		// Ebullition of CH4
		// Reduces CH4 and updates CH4_ebull_today
		if (allow_ebullition)
			calculate_gas_ebullition(CH4_ebull_today);

		// Check C budget
		calculate_carbon_store(daynum, true);
		check = ch4_store + co2_store + CH4_diff_today + CH4_plant_today + CH4_ebull_today;

		double checkafter3 = check - checkbefore;

		if (fabs(checkafter3) > MAX_ERR && verbosity >= WARNING)
			dprintf("%s%g\n","Soil::methane() - after ebullition of CH4: ",checkafter3);

		// Steps 9
		// Instead of calling diffusion and plant transport for CO2, I simply diffuse ALL the CO2. 
		// This should improve C conservation
		// In future one could look at how much IS emitted, especially as CO2 could be saved and emitted in 
		// bursts during spring thaw // DW_COMMENT something A. Pongracz could look at

		CO2_diff_today = co2_store;
		double tot_dcflux_soil = 0.0;

		// Correct & record Dec 31 data
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].CO2_soil = 0.0;
			tot_dcflux_soil += soillayers[lyr].dcflux_soil;
		}
		
		// Recalculate ch4_store and co2_store now, after CO2 diffusion
		calculate_carbon_store(daynum, true);

		CO2_plant_today = 0.0;

		// Step 9

		// Conservation steps

		// Total C flux out of the soil today [gC m-2 d-1] 
		double total_C_flux = CO2_diff_today + CH4_diff_today + CH4_ebull_today + CH4_plant_today + CO2_plant_today;

		// CO2_flux_today = CO2_diff_today + CO2_plant_today;
		// Calculate daily CH4 flux [gCH4-C m-2 d-1]
		double CH4_flux_today = CH4_diff_today + CH4_plant_today + CH4_ebull_today;

		// Reduce heterotrophic respiration by this CH4-C amount

		// Report the heterotrophic respiration and the 4 CH4 fluxes.
		patch.fluxes.report_flux(Fluxes::SOILC, tot_dcflux_soil-CH4_flux_today/G_PER_KG); // Units for Fluxes::SOILC are kgC/m2/time
		patch.fluxes.report_flux(Fluxes::CH4C, CH4_flux_today); // gC/m2
		patch.fluxes.report_flux(Fluxes::CH4C_DIFF, CH4_diff_today);
		patch.fluxes.report_flux(Fluxes::CH4C_PLAN, CH4_plant_today);
		patch.fluxes.report_flux(Fluxes::CH4C_EBUL, CH4_ebull_today);

		// Total C content of soil
		double c_soil = 0.0;
		double co2_soil = 0.0;

		calculate_carbon_store(daynum, true);
		c_soil += ch4_store + co2_store;

		if (total_C_flux < -LARGE_ERR && verbosity >= WARNING) {
			dprintf("%s%10.5f\n","Large negative C flux (total_C_flux) in Soil::methane()",total_C_flux);	
		}

		double final_C_budget = (ch4_c_store_init + co2_c_store_init + c_input) - (c_soil + total_C_flux);
		if (fabs(final_C_budget) > MAX_ERR_BALANCE && verbosity >= WARNING) {
			dprintf("%s%10.5f\n","Unbalanced C budget in Soil::methane()",total_C_flux);	
		}

		// Record today's values
		for (int lyr = 0; lyr < nsoillayers; lyr++) {
			soillayers[lyr].CO2_soil_yesterday = soillayers[lyr].CO2_soil;
			soillayers[lyr].CH4_yesterday = soillayers[lyr].CH4;
			soillayers[lyr].CH4_diss_yesterday = soillayers[lyr].CH4_diss;
			soillayers[lyr].CH4_gas_yesterday = soillayers[lyr].CH4_gas;
		}
	} // generatemethane?

	return true;
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
