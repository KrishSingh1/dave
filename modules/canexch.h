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

// Temperature response function following a modified Arrhenius equation (Johnson et al. 1942) 
double Arrhenius_temp_response(double k25, double Tl, double Ha, double Hd, double dS);

// Temperature response function for leaf mitochondrial respiration (Heskel et al. 2016, PNAS)
double respiration_temp_response(double Tl);

// Light inhibition of leaf mitochondrial respiration.
double light_inhibition(double PAR);

// System of equations to solve coupled photosynthesis-stomatal conductance scheme.
// Wang & Leuning 1998
double photosyn_solver(double g0, double X, double gamma, double beta, 
					   double gamstar, double Rd, double Ca);

// Calculate electron transport rate from PAR
// e.g. von Caemmerer 2013
double electron_transport_rate(double PAR,double Jmax);

/// Nitrogen- and landuse specific alpha a
double alphaa(const Pft& pft);

/// Calculates individual fnuptake based on surface of fine root
/** Calculates individual fraction nitrogen uptake based on surface of fine root
 *  Roots are cone formed with height == radius.
 *  V = PI * r^3 / 3
 *  A = (2^1/2 + 1) * PI * r^2
 *  -> A = const * cmass_root^2/3
 */
double nitrogen_uptake_strength(const Individual& indiv);

/// Use nitrogen storage to limit stress
void nstore_usage(Vegetation& vegetation);

// Constants for photosynthesis calculations

// atmospheric pressure (Pa)
const double PATMOS = 1e5;

/// conversion factor for solar radiation at 550 nm from J/m2 to mol_quanta/m2 (E=mol quanta); mol J-1
const double CQ = 4.6e-6;

/// intrinsic quantum efficiency of CO2 uptake, C3 plants
const double ALPHA_C3 = 0.08;

/// intrinsic quantum efficiency of CO2 uptake, C4 plants
const double ALPHA_C4 = 0.053;

/// O2 partial pressure (Pa)
const double PO2 = 2.09e4;

/// atmospheric O2 concentration (mixing ratio in mol mol-1)
const double O2atm = 0.2095;

/// colimitation (shape) parameter
const double THETA = 0.7;

/// 'saturation' ratio of intercellular to ambient CO2 partial pressure for C4 plants
const double LAMBDA_SC4 = 0.4;

/// leaf respiration as fraction of maximum carboxylation rate, C3 plants
const double BC3 = 0.015;

/// leaf respiration as fraction of maximum carboxylation rate, C4 plants
const double BC4 = 0.025;

/// Radiation conversion from J m-2 s-1 (=W m-2) to mol (quanta) m-2 s-1
const double JTOMOL = 4.6E-6;

/// Ratio of the molecular diffusivities for water vapour and CO2
const double DWDC = 1.6;

/// Molecular mass of Carbon (g mol-1)
const double MOL2G_C = 12.011;

/// Molecular mass of Water (g mol-1)
const double MOL2G_H2O = 18.01528;

/// Michaelis-Menten constants at 25 degC (mol mol-1)
const double KO25      = 278.4E-3; // MM constant for O2 (Bernacchi et al. 2001)
const double KC25      = 404.9E-6; // MM constant for CO2 (Bernacchi et al. 2001)
const double GAMSTAR25 = 42.75E-6; // photorespiratory CO2 comp. point (Bernacchi et al. 2001)

/// Activation energies (kJ mol-1)
const double HA_V  = 59.70;  // Vcmax (Kumarathunge et al. 2019)
const double HA_J  = 40.71;  // Jmax (Kumarathunge et al. 2019)
const double HA_KO = 36.38;  // Ko (Bernacchi et al. 2001)
const double HA_KC = 79.43;  // Kc (Bernacchi et al. 2001)
const double HA_GS = 37.83;  // Gamma* (Bernacchi et al. 2001)

/// Entropy terms (J mol-1 K-1)
const double DS_V = 639.43;  // Vcmax (Kumarathunge et al. 2019)
const double DS_J = 642.97;  // Jmax (Kumarathunge et al. 2019)

// Deactivation energy for Vcmax, Jmax, Vpmax (C4) (kJ mol-1)
const double HD = 200.0;

// Exponents in soil water stress function (wscal^q)
const double qs = 0.5;   // stomatal conductance
const double qb = 1.0;   // biochemistry, including Rd

/// -----------------------------------------------
/// C4 photosynthesis parameters
/// -----------------------------------------------
// Maximum Rubisco carboxylation rate at 25 degC relative to Vpmax25 (Fan et al. dataset)
const double rVC4 = 0.56;
// Maximum electron transport rate at 25 degC relative to Vcmax25 (von Caemmerer 2021)
const double rJC4 = 6.2;
// PEP regeneration rate relative to Vcmax25 (von Caemmerer 2021)
const double rVPr = 2.0;

// Activation and deactivation energy and entropy term for Vpmax (NAD-ME in Sonawane et al. 2017) 
const double HA_VP = 58.0;   // kJ mol-1
const double DS_VP = 640.0;  // J mol-1 K-1

// Activation and deactivation energy and entropy term for Vcmax (NAD-ME in Sonawane et al. 2017)
const double HA_V4 = 52.0;   // kJ mol-1 
const double DS_V4 = 630.0;  // J mol-1 K-1

// Activation and deactivation energy and entropy term for Jmax (fitted to data shown in Yamori et al. 2010 (C3 plant!!))
const double HA_J4 = 37.90;  // kJ mol-1
const double DS_J4 = 630.4;  // J mol-1 K-1

// Michaelis-Menten constant for Kp (Sonawane et al. 2017)
const double KP25    = 115.0E-6;  // (mol mol-1) // fitted assuming gm of 1.1 mol m-2 s-1.
const double HA_KP   = 38.3;      // activation energy for KP25 (kJ mol-1), DiMario & Cousins 2019

// Partitioning factor of electron transport in mesophyll and bundle-sheath cells
const double xj = 0.4;

// Fraction of cyclic electron transport (von Caemmerer 2021)
const double fcyc = 0.3;

// Bundle-sheath conductance to CO2 (mol m-2 s-1)
const double gbs = 0.003;  // von Caemmerer 2000
						   // 0.00287 in Yin et al. 2016 

// Average leaf N:P ratio for C3 grasses from Thum et al. 2019 (SI Table 8).
const double C3_LEAF_NTOP = 10.7;

// The following parameters are only used for the complex C4 photosynthesis formulation
// Michaelis-Menten constant for CO2 and O2 (Boyd et al. 2015)
const double KC25_C4 = 1210E-6;  // (mol mol-1)
const double KO25_C4 = 292E-3;   // (mol mol-1)

// Activation energies for Michaelis-Menten constants
const double HA_KC_C4 = 64.2;   // activation energy for KC_C4 (kJ mol-1), Boyd et al. 2015
const double HA_KO_C4 = 10.5;   // activation energy for KO_C4 (kJ mol-1), Boyd et al. 2015

// Half the reciprocal of Rubisco specificity
const double lgamstar = 0.0003817;  // Boyd et al. 2015

// Ratio of solubility and diffusivity of O2 to CO2
const double a0 = 0.047;

// Fraction of PSII activity in the bundle sheath
const double alpha_PSII = 0.2;  // assumed 0 in von Caemmerer 2000, but so is Rd...



/// -----------------------------------------------
/// Canopy-level parameters
/// -----------------------------------------------

/// leaf respiration as fraction of maximum rubisco, mosses
// see Wania et al. (2009b)
const double BC_moss = 0.03;

const double CMASS = 12.0;		// atomic mass of carbon
const double ALPHAA = 0.5;		// value chosen to give global carbon pool and flux values that
								// agree with published estimates.
								// scaling factor for PAR absorption from leaf to plant projective area level
								// alias "twigloss". Should normally be in the range 0-1

const double ALPHAA_NLIM = 0.6;	// Same as ALPHAA above but chosen to give pools and flux values
								// that agree with published estimates when Nitrogen limitation is
								// switched on.

const double ALPHAA_CROP = 0.7;			// Value for crops without N limitation.
const double ALPHAA_CROP_NLIM = 0.9;	// Value for crops with N limitation

/// Lambert-Beer extinction law (Prentice et al 1993; Monsi & Saeki 1953)
inline double lambertbeer(double lai, double k=0.5) {
	return exp(-k * lai);
}

// Leaf transmission coefficient (Sands, 1995)
const double m_trans = 0.0;

// Light extinction coefficient used for grasses in ifdavephotosyn
const double k_beer_grass = 0.5;

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
