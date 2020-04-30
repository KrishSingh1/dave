///////////////////////////////////////////////////////////////////////////////////////
/// \file soilwater.cpp
/// \brief Soil hydrology and snow
///
/// Original version including evaporation from soil surface, based on work by Dieter Gerten,
/// Sibyll Schaphoff and Wolfgang Lucht, Potsdam
///
/// Includes baseflow runoff
///
/// Updates include initial_infiltration, irrigation and calls to new hydrology scheme that takes into 
/// account soil water content in a greater number of layers. The Boolean iftwolayersoil (set in .ins file) 
/// determines if the original (Gerten et al) scheme is used or the newer scheme.
///
/// \author Ben Smith
/// $Date$
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

#include "config.h"
#include "soilwater.h"

void snow(double prec, double temp, Soilcolumn& soilcolumn) {

	// Daily calculation of snowfall and rainfall from precipitation and snow melt from
	// snow pack; update of snow pack with new snow and snow melt and snow melt

	// INPUT PARAMETERS
	// prec = precipitation today (mm)
	// temp = air temperature today (deg C)

	// INPUT AND OUTPUT PARAMETER
	// solicolumn = stored snow (rainfall mm equivalents) and rainfall and snow melt today (mm)

	// Add snow if snowfall (temp < 0)
	soilcolumn.snow.snowfall(temp, prec);

	// Update snow density after mechanical compaction
	soilcolumn.snow.update_density();

	// Add snow melt to precipitation
	soilcolumn.rain_melt = prec + soilcolumn.snow.melt;

	// Reset snow melt
	soilcolumn.snow.melt = 0.0;
}

/// SNOW_NINPUT
/** Nitrogen deposition and fertilization on a snowpack stays in snowpack
 *  until it starts melting. If no snowpack daily nitrogen deposition and
 *  fertilization goes to the soil available mineral nitrogen pool.
 */
void snow_ninput(double prec, double snowpack, double rain_melt,
	           double dNH4dep, double dNO3dep, double dnfert, 
			   double& snowpack_NH4_mass, double& snowpack_NO3_mass,
			   double& NH4_input, double& NO3_input) {

	// calculates this day melt and original snowpack size
	double melt = max(0.0, rain_melt - prec);
	double snowpack_total = melt + snowpack;

	// does snow exist?
	if (!negligible(snowpack_total)) {

		// if some snow is melted, fraction of nitrogen in snowpack
		// will go to soil available nitrogen pool
		if (melt > 0.0) {
			double frac_melt  = melt / snowpack_total;
			double melt_NH4_mass = frac_melt * snowpack_NH4_mass;
			NH4_input            = melt_NH4_mass + dNH4dep + dnfert / 2.0;
			snowpack_NH4_mass   -= melt_NH4_mass;

			double melt_NO3_mass = frac_melt * snowpack_NO3_mass;
			NO3_input            = melt_NO3_mass + dNO3dep + dnfert / 2.0;
			snowpack_NO3_mass   -= melt_NO3_mass;
		}
		// if no snow melts, then add daily nitrogen deposition
		// and fertilization to snowpack nitrogen pool
		else {
			snowpack_NH4_mass += dNH4dep + dnfert / 2.0;
			NH4_input = 0.0;

			snowpack_NO3_mass += dNO3dep + dnfert / 2.0;
			NO3_input = 0.0;
		}
	}
	else {
		NH4_input = dNH4dep + dnfert / 2.0;
		NO3_input = dNO3dep + dnfert / 2.0;
	}
}


/// Derive and re-distribute available rain-melt for today
/** Function to be called after interception and before canopy_exchange
 *  Calculate snowmelt
 *  If there's any rain-melt available, refill top layer, leaving any excessive
 *  rainmelt to be re-distributed later in hydrology_lpjf
 */
void initial_infiltration(Patch& patch, Climate& climate) {

	Gridcell& gridcell = climate.gridcell;
	Soilcolumn& soilcolumn = patch.soilcolumn;

	// Update snowpack
	snow(climate.prec - patch.intercep, climate.temp, soilcolumn);

	// Update snowpack N content
	snow_ninput(climate.prec - patch.intercep, soilcolumn.snow.snowpack(), soilcolumn.rain_melt, 
		        gridcell.dNH4dep, gridcell.dNO3dep, patch.dnfert,
				soilcolumn.snow.NH4_mass, soilcolumn.snow.NO3_mass, 
				soilcolumn.NH4_input, soilcolumn.NO3_input);

	// Soil water percolation
	soilcolumn.percolate = soilcolumn.rain_melt >= 0.1;

	// Save maximum rain and melt water
	soilcolumn.max_rain_melt = soilcolumn.rain_melt;

	// Reset annuals
	if (date.day == 0) {
		patch.awetland_water_added = 0.0;
	}
 
	patch.wetland_water_added_today = 0.0;

	if (soilcolumn.percolate) {

		if (iftwolayersoil) {

			// As in LPJ-GUESS v4.0, where soil.wcont[0] and soil.wcont_evap were adjusted

			Soillayer& surflayer = soilcolumn.soillayers[0];

			double twolayerwcont0 = surflayer.wcont;
			twolayerwcont0 += soilcolumn.rain_melt / surflayer.awc;

			if (twolayerwcont0 > surflayer.wsats / surflayer.wcont) {
				soilcolumn.rain_melt = (twolayerwcont0 - surflayer.wsats / surflayer.wcont) * surflayer.awc;
				twolayerwcont0 = surflayer.wsats / surflayer.wcont;
			}
			else {
				soilcolumn.rain_melt = 0.0;
			}

			// Update wcont_evap, Frac_water etc. too
			soilcolumn.update_soil_water();
		}
		else {

			if (patch.stand.is_highlatitude_peatland_stand()) {

				// PEATLAND SOILS
				// Distribute the water in the acrotelm in proportion to the capacity

				std::vector<double> Faw_layer(NACROTELM);
				// available water for each soil layer (mm)
				std::vector<double> ice_layer(NACROTELM);
				// ice each soil layer (mm)
				std::vector<double> potential_layer(NACROTELM);
				// water that can still be added to each soil layer (mm)

				double ice_fraction = 0.0; // [0-1]
				double ice_impedance = 1.0; // [0-1]

				// initialise to 0 mm
				for (int sl = 0; sl<NACROTELM; sl++) {
					Faw_layer[sl] = 0.0;
					ice_layer[sl] = 0.0;
					potential_layer[sl] = 0.0;
				}

				// Total potential of the soil to store water, from the current content up to field capacity
				double total_potential = 0.0;

				// Fraction of acrotelm pore space filled by ice
				double ice_fraction_of_pore_space = 0.0;

				for (int ly = 0; ly < NACROTELM; ly++) {

					Soillayer& soillayer = soilcolumn.soillayers[ly];

					Faw_layer[ly] = soillayer.Frac_water * soillayer.depth;
					ice_layer[ly] = soillayer.Frac_ice * soillayer.depth;
					ice_fraction += ice_layer[ly] / soillayer.aw_max / (double)NACROTELM;

					ice_fraction_of_pore_space += 1.0 / (double)NACROTELM *
						(soillayer.Frac_ice + soillayer.Fpwp_ref - soillayer.Frac_water_belowpwp) / soilcolumn.acro_por; // DW_COMMENT soillayer.Frac_por

					// Water + ice in this layer, above the wp
					double layerwater = Faw_layer[ly] + ice_layer[ly];

					potential_layer[ly] = (soilcolumn.acro_por - peat_wp) * soillayer.depth - layerwater;
					total_potential += potential_layer[ly];

					// Check balance
					if (potential_layer[ly] < -0.00001) {
						fail("initial_infiltration - error in a soil layer's water balance!\n");
					}

				} // for loop (ly)


				// Impedance - this is the fraction of water that can infiltrate initially 
				// ice_fraction = 0 (e.g. summer, no phase change etc.) = 1
				// ice_fraction = 1 (e.g. winter, daturated soil to field capacity, all soil water frozen) = 0.000000
				// CLM4.5 - see Swenson et al (2012) Eqn. (8)
				if (ICE_IMPEDANCE && patch.soilcolumn.can_freeze()) {

					ice_impedance = max(0.0, 1.0 - ice_fraction_of_pore_space);
					// CLM method: ice_impedance = pow(10.0, -6.0 * ice_fraction_of_pore_space);
				}
				else {
					ice_impedance = 1.0;
				}

				double water_in = 0.0;
				double rain_melt_orig = soilcolumn.rain_melt;

				// Assume the amount of water available for initial infiltration is limited by the ice content in the top layer
				double water_for_infiltration = ice_impedance * soilcolumn.rain_melt;

				// Only infiltrate what we can. The rest remains in rain_melt.
				if (water_for_infiltration >= total_potential) {
					water_in = total_potential;
					soilcolumn.rain_melt -= total_potential;
				}
				else {
					// Because ALL the water can infiltrate
					water_in = water_for_infiltration; // was: soil.rain_melt; before water_for_infiltration
													   //soil.rain_melt = 0.0; 
					soilcolumn.rain_melt -= water_for_infiltration; // i.e. 0.0 when there is no ice
				}

				if (total_potential > 0.0) {

					// Now add the rain_melt to each layer in proportion to the capacity if total_potential > 0.0 mm
					for (int ly = 0; ly<NACROTELM; ly++) {

						double water_input_ly = water_in * (potential_layer[ly] / total_potential);

						// Add water to the layer, and update wcont and Frac_water for this layer:
						soilcolumn.add_layer_soil_water(ly, water_input_ly); 
					}
				}
				else {
					soilcolumn.rain_melt = rain_melt_orig;
				}

			} 
			/*else if (!patch.stand.is_true_wetland_stand()) {
			
				// NEITHER PEATLAND NOR WETLAND SOILS
				// Distribute the water in the upper 50cm in proportion to the capacity:
				// The calculations rely on the fact that, in each layer: wcont = Faw_layer / soiltype.awc;

				int nsoillayers_upper = soilcolumn.nsoillayers >= 5 ? 5 : soilcolumn.nsoillayers; // DW_COMMENT could do better

				// available water for each soil layer (mm)
				std::vector<double> Faw_layer(nsoillayers_upper);
				// ice each soil layer (mm)
				std::vector<double> ice_layer(nsoillayers_upper);
				// water that can still be added to each soil layer (mm)
				std::vector<double> potential_layer(nsoillayers_upper);
		
				double ice_fraction = 0.0; // [0-1]
				double ice_impedance = 1.0; // [0-1]
				double total_potential = 0.0;

				// Average ice fraction as a fraction of pore space.

				for (int ly = 0; ly < nsoillayers_upper; ly++) {

					Soillayer& soillayer = soilcolumn.soillayers[ly];

					Faw_layer[ly] = soillayer.wcont * soillayer.awc; // mm
					ice_layer[ly] = soillayer.Frac_ice * soillayer.depth; // mm

					ice_fraction += ice_layer[ly] / soillayer.awc / (double)nsoillayers_upper;

					// Water in this layer
					double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

					potential_layer[ly] = soillayer.aw_max - layerwater;
					total_potential += potential_layer[ly];
					
					// Check balance
					if (potential_layer[ly] < -0.00001) {
						fail("initial_infiltration - error in a soil layer's water balance!\n");
					}

				} // for loop (ly)

				// Impedance - this is the fraction of water that can infiltrate initially 
				// ice_fraction = 0 (e.g. summer, no phase change etc.) = 1
				// ice_fraction = 1 (e.g. winter, daturated soil to field capacity, all soil water frozen) = 0.000000
				// CLM4.5 - see Swenson et al (2012) Eqn. (8)

				if (ICE_IMPEDANCE && patch.soilcolumn.can_freeze()) {
				
					double ice_fraction_of_pore_space = 0.0;

					for (int ly = 0; ly < soilcolumn.num_evaplayers; ly++) { // DW_COMMENT

						Soillayer& soillayer = soilcolumn.soillayers[ly];
						ice_fraction_of_pore_space += 1.0 / soilcolumn.num_evaplayers *
							(soillayer.Frac_ice + soillayer.Fpwp_ref - soillayer.Frac_water_belowpwp) / soillayer.Frac_por;
					}

					ice_impedance = max(0.0, 1.0 - ice_fraction_of_pore_space);
					// CLM method: ice_impedance = pow(10.0, -6.0 * ice_fraction_of_pore_space);
				}
				else {
					ice_impedance = 1.0;
				}

				double water_in = 0.0;
				double rain_melt_orig = soilcolumn.rain_melt;

				// Assume the amount of water available for initial infiltration is limited by the ice content in the top layer
				double water_for_infiltration = ice_impedance * soilcolumn.rain_melt;

				// Only infiltrate what we can. The rest remains in rain_melt.
				if (water_for_infiltration >= total_potential) {
					water_in = total_potential;
					soilcolumn.rain_melt -= total_potential;
				} 
				else {
					// Because ALL the water can infiltrate
					water_in = water_for_infiltration; // was: soil.rain_melt; before water_for_infiltration
					soilcolumn.rain_melt -= water_for_infiltration; // i.e. 0.0 when there is no ice
				}

				if (total_potential > 0.0) {

					// Now add the rain_melt in proportion to the capacity if total_potential > 0.0 mm
					for (int ly = 0; ly<nsoillayers_upper; ly++) {
					
						double water_input_ly = water_in * (potential_layer[ly] / total_potential);

						// Add water to the layer, and update wcont and Frac_water for this layer:
						soilcolumn.add_layer_soil_water(ly, water_input_ly);
					}
				} 
				else {
					soilcolumn.rain_melt = rain_melt_orig;
				}

			} // PEATLAND or not*/

			soilcolumn.update_soil_water(); // update wcont_evap, whc[], Frac_water etc. based on wcont (as the first layer's water content has changed) 
		}

		if (patch.stand.is_true_wetland_stand()) {

			// MINERAL WETLANDS

			// Saturate the soil
			// The calculations rely on the fact that, in each layer: wcont = Faw_layer / soiltype.awc;

			// available water for each soil layer (mm)
			std::vector<double> Faw_layer(soilcolumn.nsoillayers);
			// ice each soil layer (mm)
			std::vector<double> ice_layer(soilcolumn.nsoillayers);
			// water that can still be added to each soil layer (mm)
			std::vector<double> potential_layer(soilcolumn.nsoillayers);

			double ice_fraction = 0.0; // [0-1]
			double total_potential = 0.0;

			for (int ly = 0; ly < soilcolumn.nsoillayers; ly++) {

				Soillayer& soillayer = soilcolumn.soillayers[ly];
				Faw_layer[ly] = soillayer.wcont * soillayer.awc; // mm
				ice_layer[ly] = soillayer.Frac_ice * soillayer.depth; // mm
				ice_fraction += ice_layer[ly] / soillayer.awc * (soillayer.depth / soil_depth);

				// Water in this layer
				double layerwater = Faw_layer[ly] + ice_layer[ly]; // mm

				potential_layer[ly] = soillayer.aw_max - layerwater;
				total_potential += potential_layer[ly];

				// Check balance
				if (potential_layer[ly] < -0.00001) {
					dprintf("initial_infiltration - error in a soil layer's water balance!\n");
					return;
				}

			} // for loop (ly)

			if (soilcolumn.rain_melt < total_potential)
				soilcolumn.rain_melt = 0.0;
			else
				soilcolumn.rain_melt -= total_potential;

			if (total_potential > 0.0 && ice_fraction <= 0.0001) {

				// Now add the rain_melt in proportion to the capacity if total_potential > 0.0 mm
				for (int ly = 0; ly<soilcolumn.nsoillayers; ly++) {

					double water_input_ly = potential_layer[ly];

					// Add water to the layer, and update wcont and Frac_water for this layer:
					soilcolumn.add_layer_soil_water(ly, water_input_ly);
				}

				// Record the water added to this wetland today
				if (ifsaturatewetlands)
					patch.wetland_water_added_today = total_potential;
			}

			soilcolumn.update_soil_water(); // update wcont_evap, whc[], Frac_water etc. based on wcont

			// END MINERAL WETLANDS
		}

	}
}

/// Calculate required irrigation according to water deficiency.
/** Function to be called after canopy_exchange and before soilwater.
 */
void irrigation(Patch& patch) {

	Soilcolumn& soilcolumn = patch.soilcolumn;

	patch.irrigation_d = 0.0;
	if (date.day == 0) {
		patch.irrigation_y = 0.0;
	}

	if (!patch.stand.isirrigated) {
		return;
	}

	for (int i = 0; i < npft; i++) {

		Patchpft& ppft = patch.pft[i];
		if (patch.stand.pft[i].irrigated && ppft.growingseason()) {
			if (ppft.water_deficit_d < 0.0) {
				fail("irrigation: Negative water deficit for PFT %s!\n", (char*)ppft.pft.name);
			}
			patch.irrigation_d += ppft.water_deficit_d;
		}
	}
	patch.irrigation_y += patch.irrigation_d;
	soilcolumn.rain_melt += patch.irrigation_d;
	soilcolumn.max_rain_melt += patch.irrigation_d;
}


/// Performs daily accounting of soil water
/** 
 * Call this function each simulation day for each modelled area or patch, following
 * calculation of vegetation production and evapotranspiration and before soil organic
 * matter and litter dynamics
 * \param patch      The patch to simulate
 * \param climate   The climate to use to update soil water
 */
void soilwater(Patch& patch, Climate& climate) {

	// DESCRIPTION
	// Performs daily accounting of soil water

	// update the daily snow depth
	// by converting from mm water to snow depth
	Soilcolumn& soilcolumn = patch.soilcolumn;

	// Sum vegetation phenology-weighted FPC
	// Fraction of grid cell subject to evaporation from soil surface is
	// complement of summed vegetation projective cover (FPC)
	double fpc_phen_total = 0.0;	// phenology-weighted FPC sum for patch

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		fpc_phen_total += indiv.fpc_today();

		vegetation.nextobj();
	}
	
	double bare_ground = max(1.0 - fpc_phen_total, 0.0);

	// HYDROLOGY FOR HIGH-LATITUDE PEATLAND 

	if (patch.stand.is_highlatitude_peatland_stand()) {
        
		soilcolumn.hydrology_peat(climate, bare_ground);

		// Update inundation stress variables for this patch
		for (int p=0;p<npft;p++) {

			double wtpp = patch.soilcolumn.wtp[date.day]; // [-300,100] mm
			double wtpm = patch.pft[p].pft.wtp_max; // mm

			if (date.day == 0) 
				patch.pft[p].inund_count = 0; // Reset on Jan 1

			bool wania_graminoid_inundation = false; // Follows Wania et al. (2009b) if true

			if (wania_graminoid_inundation) { 

				// Following Wania et al. (2009b)

				if (date.dayofmonth == 0) 
					patch.pft[p].inund_count = 0; // Reset on the 1st day of every month too

				if (wtpp > wtpm) {
					patch.pft[p].inund_count++; // Days this month with stress
				}
			} 
			else {

				// New inundation algorithm that doesn't depend on months 

				// Inundation restrictions only applied when phen > 0
				if (patch.pft[p].phen > 0.0) {
					if (wtpp > wtpm) {
						patch.pft[p].inund_count++; // days
					} 
					else {
						patch.pft[p].inund_count--; 
					}
				} 
				else {
					patch.pft[p].inund_count = 0;
				}

				const int inundation_delay = 3; // days
				// Alternative: could restrict inund_count to be between 0 and this PFT's upper limit + 3 days.
				// Ensures that the wetland PFTs benefit from a drop in the water table after a delay of 3 days
				patch.pft[p].inund_count = max(min(patch.pft[p].inund_count,patch.pft[p].pft.inund_duration + inundation_delay),0);
			}

			// Inundation stress is updated daily, between 0 (full stress) to 1 (no stress), and used in 
			// photosynthesis to reduce daily GPP.
			
			double pft_inund_dur = (double)patch.pft[p].pft.inund_duration;
			double pft_inund_count = (double)patch.pft[p].inund_count;

			if (!ifinundationstress) {
				patch.pft[p].inund_stress = 1.0; // i.e. never any stress if this is 0 in the .ins file
			}
			else {
				if (pft_inund_dur <= 0) { // constant from from .ins file
					// FULL stress, so NO photosynthesis and this PFT will not survive on peatlands 
					patch.pft[p].inund_stress = 0.0;
				}
				else {
					if (pft_inund_count <= 0.0)
						patch.pft[p].inund_stress = 1.0; // i.e. no stress
					else
						patch.pft[p].inund_stress = 1.0 - min(1.0, pft_inund_count / pft_inund_dur); // partial stress, 0 (full) to 1 (none) 
				}
			}
		} // pft loop
	} 
	else {	
		
		// HYDROLOGY FOR MINERAL SOILS 
		if (iftwolayersoil)
			soilcolumn.hydrology_lpjf_twolayer(climate, bare_ground);
		else
			soilcolumn.hydrology_lpjf(climate, bare_ground); // DW_COMMENT NEW LSM HYDROLOGY
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Bondeau, A., Smith, P.C., Zaehle, S., Schaphoff, S., Lucht, W., Cramer, W.,
//   Gerten, D., Lotze-Campen, H., Müller, C., Reichstein, M. and Smith, B. (2007),
//   Modelling the role of agriculture for the 20th century global terrestrial carbon balance.
//   Global Change Biology, 13: 679-706. doi: 10.1111/j.1365-2486.2006.01305.x
// Rost, S., D. Gerten, A. Bondeau, W. Luncht, J. Rohwer, and S. Schaphoff (2008),
//   Agricultural green and blue water consumption and its influence on the global
//   water system, Water Resour. Res., 44, W09405, doi:10.1029/2007WR006331
// Swenson, S. C., D. M. Lawrence, and H. Lee (2012), Improved simulation of the terrestrial 
//   hydrological cycle in permafrost regions by the Community Land Model, 
//   J.Adv.Model.Earth Syst., 4, M08002, doi:10.1029 / 2012MS000165
// Wania, R., Ross, I., & Prentice, I.C. (2009a) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: I. Evaluation and sensitivity of physical 
//   land surface processes. Global Biogeochemical Cycles, 23, GB3014, doi:10.1029/2008GB003412
// Wania, R., Ross, I., & Prentice, I.C. (2009b) Integrating peatlands and permafrost 
//   into a dynamic global vegetation model: II. Evaluation and sensitivity of vegetation 
//   and carbon cycle processes. Global Biogeochemical Cycles, 23, GB015, doi:10.1029/2008GB003413
