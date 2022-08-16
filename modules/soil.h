///////////////////////////////////////////////////////////////////////////////////////
/// \file soil.h
/// \brief Constants and parameters used in Arctic and wetland code, with references.
/// NB: The class Soil and its member functions and variables are declared in guess.h, 
/// while its member functions are implemented in soil.cpp and in soilmethane.cpp.
/// 
/// \author Paul Miller
/// $Date: 2019-03-10 16:34:14 +0100 (Sun, 10 Mar 2019) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_SOIL_H
#define LPJ_GUESS_SOIL_H

// DEBUGGING BOOLEANS

const bool DEBUG_SOIL_WATER = false;
const bool DEBUG_SOIL_TEMPERATURE = false;
const bool DEBUG_METHANE = false;


// CONSTANTS

/// Number of soil layers used in LPJ-GUESS v4.0 for soil water calculations. Typically 2 layers, 500mm + 10000mm, making up the 150cm-deep soil column
const int NSOILLAYERS_SIMPLE = 2;

/// maximum number of soil layers modelled
const int NSOILLAYERS_MAX = 15;

/// number of padding layers in the soil
const int NPADLAYERS = 5;

/// depth of the acrotelm (mm), catotelm rest of soil column
const double ACROTELM_DEPTH = 300.0;

/// Total depth [mm] of padding layers - set to 8000 when analyticalSolutionTest = true
const double PAD_DEPTH = 50000;

/// Depth of each mineral soil layer [mm]
const double Dz_soil = 100.0;

/// Depth of each acrotelm soil layer [mm]
const double Dz_acro = 100.0;

/// Depth of each catotelm soil layer [mm]
const double Dz_cato = 100.0;

/// Max height of standing water [mm]
const double maxh = 0.0;

/// Density of mineral particles [kg m-3] (Bonan 2019)
const double Mindensity = 265.0;

// C fraction of organic matter (56%) (Chadburn et al. 2022)
const double fc = 0.56;

// For dpm and rpm [kg C m-3] (Chadburn et al. 2022) DW_COMMENT
const double rho_low = 35.0 * fc;

// For bio and hum [kg C m-3] (Chadburn et al. 2022) DW_COMMENT
const double rho_high = 210.0 * fc;

// SNOW PARAMETERS

/// snow density at start of the snow season [kg m-3] (Wania et al. (2009a) have 150, Best et al. (2012) have 50)
const double snowdens_start = 275.0;

/// snow density at end of the snow season [kg m-3] (Wania et al. 2009a)
const double snowdens_end = 500.0;

/// maximum number of snow layers allowed
const int NSNOWLAYERS_MAX = 5;

// Standard depth of a snow layer
const double SNOW_DEPTH_STANDARD = 100.0;

// Maximum snowpack depth allowed [mm]
const double MAX_SNOWPACK_DEPTH = 2000.0;

// Maximum number of layers allowed
const int NLAYERS_MAX = NSNOWLAYERS_MAX + NSOILLAYERS_MAX + NPADLAYERS;

// WATER PARAMETERS

/// density of ice [kg m-3]
const double rho_ICE = 917.0;

/// density of water [kg m-3]
const double rho_H2O = 1000.0;


// POROSITIES

/// catotelm porosity
const double catotelm_por = 0.92;

/// acrotelm porosity
const double acrotelm_por = 0.98;

/// Gas fraction in peat (Possible to reintroduce - was 0.08 in Wania et al (2010))
const double Fgas = 0.00;

/// Wilting point in peat
const double peat_wp = 0.066;

/// Threshold for randomly distributed organic and mineral fractions for connecting flow pathways
const double threshold_frac = 0.5;

/// Variable for power function of organic fraction for flow pathways
const double beta_perc = 0.139;

/// Depth that organic matter takes on the characteristics of sapric peat
const double sapric_depth = 500.0;

/// First year when phase change is allowed
const int FIRST_FREEZE_YEAR = 40;

/// time step [day]
const double Dt = 1;

// HEAT CAPACITIES

/// heat capacity of air [J m-3 K-1] - Bonan (2002)
const double Cp_air = 1230;

/// heat capacity of water [J m-3 K-1] - Bonan (2002)
const double Cp_water = 4180000;
 
/// heat capacity of ice [J m-3 K-1] = 2117.27 [J kg-1 K-1] * 917 [kg m-3 (ice density)] - CLM
const double Cp_ice = 1941537;		 

/// heat capacity of organic matter [J m-3 K-1]
const double Cp_org = 2500000;

/// heat capacity of dry peat [J m-3 K-1] - Bonan (2002), 0% water
const double Cp_peat = 580000;

/// heat capacity of moss [J m-3 K-1] Ekici et al. (2015)
const double Cp_moss = 2500000;

/// heat capacity of wood (spruce, 0% water) [J m-3 K-1] Božiková et al. (2021)
const double Cp_wood = 480000;

// NOTE: using the Cp_org values for Cp_peat and Korg for Kpeat does not
// seem to influence the upper soil layer Ts, but increases the range
// in the lower layers (2m)

// THERMAL CONDUCTIVITIES
// Values are from Hillel (1982) unless otherwise stated

/// thermal conductivity of air [W m-1 K-1]
const double Kair = 0.025;
const double Kair_log = log(Kair);

/// thermal conductivity of water [W m-1 K-1]
const double Kwater = 0.57;
const double Kwater_log = log(Kwater);

/// thermal conductivity of ice [W m-1 K-1]
const double Kice = 2.2;
const double Kice_log = log(Kice);

/// thermal conductivity of organic matter [W m-1 K-1]
const double Korg = 0.25;
const double Korg_log = log(Korg);

/// thermal conductivity of dry peat [W m-1 K-1] - Bonan (2002)
const double Kpeat = 0.06;
const double Kpeat_log = log(Kpeat);

/// thermal conductivity of mineral soil [W m-1 K-1] - Wania et al. (2009a)
const double Kmin = 2.0;
const double Kmin_log = log(Kmin);

/// thermal conductivity of moss [W m-1 K-1]
const double Kmoss = 0.25;
const double Kmoss_log = log(Kmoss);

/// thermal conductivity of wood (spruce, 0% water) [W m-1 K-1] Božiková et al. (2021)
const double Kwood = 0.075;
const double Kwood_log = log(Kmoss);

/// latent heat of fusion [J m-3]
const double Lheat = 3.0587E8;	

/// A FILL-IN for layers above layer0.
const double MISSING_VALUE = -9999.0;      

/// Thickness of the topmost air layer [mm]
// Note - keep air thickness low when using cnstep_full
const double AIR_THICKNESS = 100.0;


////////////////////////////////////////////////////////////////////////

// SOIL CONSTANTS

/// Increased soil decomposition due to tillage (crops)
// (Inverse of "tillage factor" in Chatskikh et al. 2009, Value selected by T.Pugh)
static const double TILLAGE_FACTOR = 33.0 / 17.0;

/// Increased soil decoposition due to tillage only down to a certain depth  
static const double TILLAGE_DEPTH = 300.0;

/// reduce infiltration and percolation rates sharply when there is ice present in the soil
// See CLM4.5 - Swenson et al. (2012)
const bool IMASS_IMPEDANCE = false;
	
/// min temp [deg C] for heterotrophic decomposition
const double MIN_DECOMP_TEMP = -20.0;		
	
/// max CO2 [mimol/L] available to mosses in the acrotelm - see Wania et al (2009b)
// value taken from Smolders et al. 2001 which give an average CO2 conc. of 70 sites as 934 mimol L-1
const double PORE_WATER_CO2 = 934.0; 

// METHANE CONSTANTS 

/// optimised moisture response under inundation (Wania et al. 2010, Table 5)
const double RMOIST = 0.4;

/// Frolking et al (2001, 2010), Ise et al. (2008)
const double RMOIST_ANAEROBIC=0.025; 
	
/// time step for gas diffusion calculations [day]
const double Dt_gas = 0.01;
	
/// CH4:CO2 ratio for peatland soils (> PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Wania et al (2010)
// Wania et al. (2010) optimal value: 0.1 (see Table 4). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.25, after optimisation
const double CH4toCO2_peat = 0.085; 
	
/// CH4:CO2 ratio for inundated soils (< PEATLAND_WETLAND_LATITUDE_LIMIT N) - See Spahni et al. (2011)
//const double CH4toCO2_inundated = 0.024; // SC1 value in Spahni et al. SC2 is 0.0415
//const double CH4toCO2_inundated = 0.027; // Updated from Spahni et al. (2011) to match global emissions   
const double CH4toCO2_inundated = 0.0415; // Spahni et al. (2011)           
	
/// acceleration due to gravity [m s-2]
const double gravity = 9.81;

/// Molecular mass of CH4 [g mol-1]
const double molecularmass_CH4 = 16.0;

/// molecular weight of water [g mol-1]
const double molecularmass_H2O = 18.0;

/// universal gas constant [J mol-1 K-1]
const double R_gas = 8.314472;

/// standard atmospheric pressure [Pa]
const double atm_press = 101325.0;

/// coefficient for the calculation of the gas transport velocity, given in Riera et al. 1999
const double n_coeff = -0.5;
      
///  Henry's Law constants [L atm mol-1] at 298.15K. Wania et al. (2010), Table 2
const double henry_k_CO2 = 29.41;
const double henry_k_CH4 = 714.29;
const double henry_k_O2 = 769.23;
     
/// Constants [K] for CO2, CH4 and O2 for calculation of Henry's coefficient cited by Sander (1999). Wania et al. (2010), Table 2
const double henry_C_CO2 = 2400.0;
const double henry_C_CH4 = 1600.0;
const double henry_C_O2 = 1500.0;

/// partial pressure of CH4 above water
const double pp_CH4 = 1.7; // micro atm

/// partial pressure of O2 above water (value consistent with PO2 in canexch.h)
const double pp_O2 = 209000; // micro atm

/// when ebullition occurs, the volumetric gas content (VGC) will drop to this level [unitless]
const double vgc_low = 0.145;

///	ebullition occurs, when the volumetric gas content (VGC) exceeds this level [unitless, m3/m3]
const double vgc_high = 0.15;          

/// CH4 fraction of gas bubbles [unitless]
const double bubble_CH4_frac = 0.57;

/// Fraction of oxygen used to oxidise CH4
// Wania et al. (2010) optimal value: 0.5 (see Table 5). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.9, after optimisation
const double oxid_frac = 0.5;

/// Fraction of ANPP used to calculate number of tillers
const double ag_frac = 0.4;

/// Radius of an average tiller [m]
// (tiller_radius = 0.004)  ! Schimel (1995) - Average over E. angustifolium (diam=7.9mm)
// and C. aquatilis (diam=3.8mm)
// Wania et al. (2010) optimal value: 0.003mm (see Table 5). 
// McGuire et al (2012), Tang et al (2015) and Zhang et al (2013) use 0.0035, after optimisation
const double tiller_radius = 0.0035;       


/// Tiller porosity
// (tiller_por = 0.6) ! Wetland plants book, eds. Cronk and Fennessy, p.90, values for 2 Erioph. spp.
const double tiller_por = 0.7; // Wania et al. (2010) optimal value: 0.7 (see Table 5)
   
/// C content of biomass
const double c_content = 0.45;

/// atomic mass of carbon [g/mol]
const double atomiccmass = 12.0;

/// Individual tiller weight [g C]
const double tiller_weight = 0.22;
     
/// a threshold factor for a minimum water content in the layer [unitless]
const double Fwmass_min = 0.1;

/// Latitude (N). North of this and PEATLAND stands are treated as peatland as in Wania et al. (2009a, 2009b, 2010)
// But south of this, then the PEATLAND stands are irrigated to avoid water stress, and can be a source of methane
const double PEATLAND_WETLAND_LATITUDE_LIMIT = 40.0;

// INLINE FUNCTIONS

/*inline void tridiag(int n, double a[NSOILLAYER_MAX], double b[NSOILLAYER_MAX], double c[NSOILLAYER_MAX],
	double r[NSOILLAYER_MAX], double u[NSOILLAYER_MAX]) {

	// Tridiagonal system solver from Numerical Recipes.

	double gam[NSOILLAYER_MAX];
	double bet;

	bet = b[0];

	u[0] = r[0] / bet;

	for (int j = 1; j<n; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];

		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}

	for (int j = (n - 2); j >= 0; j--) {
		u[j] -= gam[j + 1] * u[j + 1];
	}
}*/

inline void tridiag(int n, double a[], double b[], double c[], double r[], double u[]) {

	// Tridiagonal system solver from Numerical Recipes.

	double gam[NLAYERS_MAX + 1];
	double bet;

	bet = b[0];

	u[0] = r[0] / bet;

	for (int j = 1; j < n; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];

		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}

	for (int j = (n - 2); j >= 0; j--) {
		u[j] -= gam[j + 1] * u[j + 1];
	}
}

#endif //LPJ_GUESS_SOIL_H
