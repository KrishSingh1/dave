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

void snow(Soilcolumn& soilcolumn, double prec, double temp) {

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

	// Soil water percolation
	soilcolumn.percolate = soilcolumn.rain_melt >= 0.1;

	// Save maximum rain and melt water
	soilcolumn.max_rain_melt = soilcolumn.rain_melt;
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

	// Water mass balance check
	double wbefore = soilcolumn.wmass() + soilcolumn.snow.snowpack() + climate.prec - patch.intercep;

	// Update snowpack
	snow(soilcolumn, climate.prec - patch.intercep, climate.temp);

	// Update snowpack N content
	snow_ninput(climate.prec - patch.intercep, soilcolumn.snow.snowpack(), soilcolumn.rain_melt, 
		        gridcell.dNH4dep, gridcell.dNO3dep, patch.dnfert,
				soilcolumn.snow.NH4_mass, soilcolumn.snow.NO3_mass, 
				soilcolumn.NH4_input, soilcolumn.NO3_input);

	// Percolate initial water
	soilcolumn.initial_percolation();

	// Water mass balance check
	double wafter = soilcolumn.wmass() + soilcolumn.snow.snowpack() + soilcolumn.rain_melt;
	if (fabs(wbefore - wafter) > 1.0e-10 && date.day == 0 && date.year == 0)
		dprintf("Year %d Day %d INITIALINFILT Before %g After %g Diff %g\n", date.year, date.day, wbefore, wafter, wbefore - wafter);
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
			if (ppft.wmass_deficit_d < 0.0) {
				fail("irrigation: Negative water deficit for PFT %s!\n", (char*)ppft.pft.name);
			}
			patch.irrigation_d += ppft.wmass_deficit_d;
		}
	}
	patch.irrigation_y += patch.irrigation_d;
	soilcolumn.rain_melt += patch.irrigation_d;
	soilcolumn.max_rain_melt += patch.irrigation_d;
}

// Inundation stress
/* Inundation stress is updated daily, between 0 (full stress) to 1 (no stress), and used in 
 * photosynthesis to reduce daily GPP
*/
void inundation_stress(Patch& patch) {
	for (int p = 0; p<npft; p++) {

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
			patch.pft[p].inund_count = max(min(patch.pft[p].inund_count, patch.pft[p].pft.inund_duration + inundation_delay), 0);
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

	Soilcolumn& soilcolumn = patch.soilcolumn;

	// Sum vegetation phenology-weighted FPC
	// Fraction of grid cell subject to evaporation from soil surface is
	// complement of summed vegetation projective cover (FPC)

	// Phenology-weighted FPC sum for patch
	double fpc_phen_total = 0.0;	

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		fpc_phen_total += indiv.fpc_today();

		vegetation.nextobj();
	}
	
	// Fraction of ground which isn't coverred by vegetation
	double bare_ground = max(1.0 - fpc_phen_total, 0.0);

	// HYDROLOGY FOR HIGH-LATITUDE PEATLAND 

	if (patch.stand.is_highlatitude_peatland_stand()) {
        
		// Hydrology for peatland soils
		soilcolumn.hydrology_peat(climate, bare_ground);

		// Update inundation stress variables for this patch
		inundation_stress(patch);
	} 
	else {	
		
		// HYDROLOGY FOR MINERAL SOILS 
		if (iftwolayersoil)
			// Standard two layer LPJ-GUESS hydrology
			soilcolumn.hydrology_lpjf_twolayer(climate, bare_ground);
		else
			// New LSM hydrology based on HTESSEL
			soilcolumn.hydrology_lpjf(climate, bare_ground);
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
