///////////////////////////////////////////////////////////////////////////////////////
/// \file canexch.h
/// \brief The canopy exchange module header file
///
/// Vegetation-atmosphere exchange of H2O and CO2 via
/// production, respiration and evapotranspiration.
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_CANEXCH_H
#define LPJ_GUESS_CANEXCH_H

#include "guess.h"

void interception(Patch& patch, Climate& climate);
void canopy_exchange(Patch& patch, Climate& climate);
void photosynthesis(const PhotosynthesisEnvironment& ps_env, 
					const PhotosynthesisStresses& ps_stresses,
					const Pft& pft,
					double lambda, 
					double nactive, 
					double vm,
					PhotosynthesisResult& result);


/// Nitrogen- and landuse specific alpha a
double alphaa(const Pft& pft);

// Constants for photosynthesis calculations

/// conversion factor for solar radiation at 550 nm from J/m2 to mol_quanta/m2 (E=mol quanta); mol J-1
const double CQ = 4.6e-6;

/// intrinsic quantum efficiency of CO2 uptake, C3 plants
const double ALPHA_C3 = 0.08;

/// intrinsic quantum efficiency of CO2 uptake, C4 plants
const double ALPHA_C4 = 0.053;

/// O2 partial pressure (Pa)
const double PO2 = 2.09e4;

/// colimitation (shape) parameter
const double THETA = 0.7;

/// 'saturation' ratio of intercellular to ambient CO2 partial pressure for C4 plants
const double LAMBDA_SC4 = 0.4;

/// leaf respiration as fraction of maximum rubisco, C3 plants
const double BC3 = 0.015;

/// leaf respiration as fraction of maximum rubisco, C4 plants
const double BC4 = 0.02;

/// Radiation conversion from J m-2 s-1 (=W m-2) to mol (quanta) m-2 s-1
const double JTOMOL = 4.6E-6;

/// Ratio of the molecular diffusivities for water vapour and CO2
const double DWDC = 1.6;

/// Molecular mass of Carbon (g mol-1)
const double CMOL2G = 12.011;

/// Michaelis-Menten constants at 25 degC (mol mol-1)
const double KO25      = 278.4E-3; // MM constant for O2 (Bernacchi et al. 2001)
const double KC25      = 404.9E-6; // MM constant for CO2 (Bernacchi et al. 2001)
const double GAMSTAR25 = 42.75E-6; // photorespiratory CO2 comp. point (Bernacchi et al. 2001)

/// Activation energies (kJ mol-1)
const double HA_V  = 59.70;  // Vcmax (Kumarathunge et al. 2019)
const double HA_J  = 40.71;  // Jmax (Kumarathunge et al. 2019)
const double HA_R  = 25.49;  // Rd (Kumarathunge et al. 2019) TODO: check!
const double HA_KO = 36.38;  // Ko (Bernacchi et al. 2001)
const double HA_KC = 79.43;  // Kc (Bernacchi et al. 2001)
const double HA_GS = 37.83;  // Gamma* (Bernacchi et al. 2001)

/// Entropy terms (J mol-1 K-1)
const double DS_V = 639.43;  // Vcmax (Kumarathunge et al. 2019)
const double DS_J = 642.97;  // Jmax (Kumarathunge et al. 2019)

// Deactivation energy for Vcmax, Jmax (kJ mol-1)
const double HD = 200.0;


/// -----------------------------------------------
/// C4 photosynthesis parameters
/// -----------------------------------------------
// Maximum Rubisco carboxylation rate (mol m-2 s-1)
const double Vcmax25_C4 = 30.0E-6;

// Activation and deactivation energy and entropy term for Vcmax (Massad et al. 2007)
const double HA_V4 = 67.29;   // kJ mol-1
const double HD_V4 = 144.57;  // kJ mol-1
const double DS_V4 = 472.0;   // J mol-1 K-1

// Maximum electron transport rate (mol m-2 s-1)
const double Jmax25_C4 = 200.0E-6; // assumed same Vcmax:Jmax ratio as in von Caemmerer 2000

// Activation and deactivation energy and entropy term for Jmax (Massad et al. 2007)
const double HA_J4 = 77.90;   // kJ mol-1
const double HD_J4 = 191.93;  // kJ mol-1
const double DS_J4 = 627.0;   // J mol-1 K-1

// Maximum PEP carboxylation rate at 25 degC (mol m-2 s-1)
const double Vpmax25 = 120.0E-6; 

// Activation and deactivation energy and entropy term for Vpmax (Massad et al. 2007) 
const double HA_VP = 70.37;  // kJ mol-1
const double HD_VP = 117.91; // kJ mol-1
const double DS_VP = 376.0;  // J mol-1 K-1 

// PEP regeneration rate (mol m-2 s-1)
const double Vpr = 80.0E-6;

// Michaelis-Menten constant of PEP carboxylase 
const double KP25 = 80.0E-6;  // (mol mol-1) (von Caemmerer 2000)
const double HA_KP = 68.1;    // activation energy for KP25 (kJ mol-1)
							  // Yin et al. 2016

// Partitioning factor of electron transport in mesophyll and bundle-sheath cells
const double xj = 0.4;

// Bundle-sheath conductance to CO2 (mol m-2 s-1)
const double gbs = 0.003;  // von Caemmerer 2000
						   // 0.00287 in Yin et al. 2016 







/// leaf respiration as fraction of maximum rubisco, mosses
// see Wania et al. (2009b)
const double BC_moss = 0.03;

const double CMASS = 12.0;		// atomic mass of carbon
const double ALPHAA = 0.45;		// value chosen to give global carbon pool and flux values that
								// agree with published estimates.
								// scaling factor for PAR absorption from leaf to plant projective area level
								// alias "twigloss". Should normally be in the range 0-1

const double ALPHAA_NLIM = 0.6; // Same as ALPHAA above but chosen to give pools and flux values
								// that agree with published estimates when Nitrogen limitation is
								// switched on.

const double ALPHAA_CROP = 0.65;		// Value for crops without N limitation.
const double ALPHAA_CROP_NLIM = 0.85;	// Value for crops with N limitation

/// Lambert-Beer extinction law (Prentice et al 1993; Monsi & Saeki 1953)
inline double lambertbeer(double lai) {
	return exp(-.5 * lai);
}

/// Alternative parameterisations of the convective boundary layer
/**
 *	AET_MONTEITH_HYPERBOLIC = hyperbolic parameterisation (Huntington & Monteith 1998)
 *	AET_MONTEITH_EXPONENTIAL = exponential parameterisation (Monteith 1995)
 *	aet_monteith		Returns AET given equilibrium evapotranspiration and canopy conductance
 *	gc_monteith			Returns canopy conductance given AET and equilibrium evapotranspiration
 */

// Comment out one of the following two lines:
#define AET_MONTEITH_HYPERBOLIC
//#define AET_MONTEITH_EXPONENTIAL

// Check:
#if defined(AET_MONTEITH_HYPERBOLIC) && defined(AET_MONTEITH_EXPONENTIAL)
#error Only one of AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL should be #defined
#elif !defined(AET_MONTEITH_HYPERBOLIC) && !defined(AET_MONTEITH_EXPONENTIAL)
#error One of AET_MONTEITH_HYPERBOLIC and AET_MONTEITH_EXPONENTIAL must be #defined
#endif

#if defined(AET_MONTEITH_EXPONENTIAL)

const double ALPHAM = 1.4;
const double GM = 5.0;

inline double aet_monteith(double& eet, double& gc) {
	return negligible(gc) ? 0.0 : eet*ALPHAM*(1.0-exp(-gc/GM));
}

inline double gc_monteith(double& aet, double& eet) {
	if (negligible(eet)) return 0.0;
	double t = aet/eet/ALPHAM;
	if (t >= 1.0) fail("gc_monteith: invalid value for aet/eet/ALPHAM");
	return -GM * log(1.0 - t);
}
#elif defined(AET_MONTEITH_HYPERBOLIC)

const double ALPHAM = 1.391;
const double GM = 3.26;

inline double aet_monteith(double& eet, double& gc) {
	return eet*ALPHAM*gc/(gc+GM);
}

inline double gc_monteith(double& aet, double& eet) {
	return (aet*GM) / (eet*ALPHAM-aet);
}
#endif

#endif // LPJ_GUESS_CANEXCH_H
