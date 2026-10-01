#include "davephotosynthesis.h"
#include "hydraulics.h"
#include "guess.h"
#include "ncutil.h"
#include "canexch.h" // Many constants defined here
#include "logging.h"
#include "log_scope.h"
#include "guessmath.h"
#include "driver.h"
#include "somdynam.h"

// Average leaf N:P ratio for C3 grasses from Thum et al. 2019 (SI Table 8).
const double C3_LEAF_NTOP = 10.7;

// Depth (mm) to which forest-floor photosynthesis has access to water.
const double FOREST_FLOOR_WSCAL_DEPTH = 500;

// Scaling parameter for gcutic (Cochard 2021 Eq 13 & 14).
const double Q10A = 1.2;

// Scaling parameter for gcutic (Cochard 2021 Eq 13 & 14).
const double Q10B = 4.8;

// Transition temperature degC (to discuss with Erick) (Cochard 2021 Eq 13 & 14)
const double TP = 40;

struct EnergyBalanceResult {
	// Leaf temperature (degC).
	double tleaf;
	// Transpiration (mol H2O m-2 s-1).
	double transp;
	// Net radiation (W m-2).
	double rnet;
};

class ClimateInputs {
public:
	ClimateInputs(const Climate &climate, const Day &day) :
		tair(climate.temps[day.period]),
		// May not have subdaily wind speed. In not, fall back to daily mean.
		wind(climate.winds.size() > day.period ? climate.winds[day.period] : climate.u10), // convert from km/h to m/s
		patm(climate.pressure[day.period] / 1000.0), // convert from Pa to kPa
		vpd(climate.vpds[day.period]) {
	}
	double get_tair() const { return tair; }
	double get_wind() const { return wind; }
	double get_patm() const { return patm; }
	double get_vpd() const { return vpd; }
private:
	// Air temperature (degC).
	double tair;
	// Wind speed (m/s).
	double wind;
	// Air pressure (kPa).
	double patm;
	// Vapour pressure deficit (kPa).
	double vpd;
};

/*
Describes the state of a cohort relevant for photosynthesis calculations.
*/
class CohortState {
public:
	CohortState(double lai, double wscal, double lwp, double predawn_lwp,
				double stom_clos, double fipar, bool stomata_open,
				double fphen_activity, double ps_sink_lim, double vcmax25,
				double jmax25, double vpmax25) :
		lai(lai), wscal(wscal), lwp(lwp), predawn_lwp(predawn_lwp),
		stom_clos(stom_clos),fipar(fipar),
		stomata_open(stomata_open), fphen_activity(fphen_activity),
		ps_sink_lim(ps_sink_lim), vcmax25(vcmax25), jmax25(jmax25),
		vpmax25(vpmax25) {
	}

	CohortState(const Individual &indiv, size_t layer) : CohortState(
		indiv.canopy[layer].lai,
		indiv.get_wscal(),
		ifdavehydraulics ? indiv.lwp_previous : indiv.weighted_psi,
		indiv.predawn_lwp,
		indiv.pft.stom_clos,
		indiv.canopy[layer].fipar,
		indiv.stomata_open,
		indiv.fphen_activity,
		indiv.ps_sink_lim,
		indiv.ps_subdaily.vcmax25,
		indiv.ps_subdaily.jmax25,
		indiv.ps_subdaily.vpmax25
	) { }
	double get_lai() const { return lai; }
	double get_wscal() const { return wscal; }
	double get_lwp() const { return lwp; }
	double get_predawn_lwp() const { return predawn_lwp; }
	double get_stom_clos() const { return stom_clos; }
	double get_fipar() const { return fipar; }
	bool get_stomata_open() const { return stomata_open; }
	double get_fphen_activity() const { return fphen_activity; }
	double get_ps_sink_lim() const { return ps_sink_lim; }
	double get_vcmax25() const { return vcmax25; }
	double get_jmax25() const { return jmax25; }
	double get_vpmax25() const { return vpmax25; }
private:
	// LAI (m2/m2).
	double lai;
	// Water stress (0-1) (0=max stress, 1=no stress).
	double wscal;
	// Leaf water potential (MPa).
	double lwp;
	// Predawn leaf water potential (MPa).
	double predawn_lwp;
	// Leaf water potential (MPa) at stomatal closure.
	double stom_clos;
	// Fraction of PAR incident to the leaf area (0-1).
	double fipar;
	// Are the stomata open?
	bool stomata_open;
	// fphen_activity (phenological state for grass. 0=dormant, 1=fully active).
	double fphen_activity;
	// Sink limitation scalar for photosynthesis (0-1).
	double ps_sink_lim;
	// Maximum carboxylation rate at 25 degC (μmol m-2 s-1).
	double vcmax25;
	// Maximum electron transport rate at 25 degC (μmol m-2 s-1).
	double jmax25;
	// Maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1).
	double vpmax25;
};

/*
Get vcmax25, jmax25, and vpmax25 from PFT parameters and store them in
indiv.ps_subdaily.
*/
void vcmax_from_parameters(Individual &indiv) {
	if (indiv.pft.pathway == C3) {
		indiv.ps_subdaily.vcmax25 = indiv.pft.vcmax_25;
		indiv.ps_subdaily.jmax25 = indiv.pft.jmax_25;
		indiv.ps_subdaily.vpmax25 = 0;
	} else {
		indiv.ps_subdaily.vcmax25 = 0;
		indiv.ps_subdaily.jmax25 = 0;
		indiv.ps_subdaily.vpmax25 = indiv.pft.vpmax_25;
	}
}

/*
Calculate leaf N content at the top of the canopy (kgN m-2 per green leaf area).
*/
double calc_nmass_leaf0(const Pft &pft, double lai, double nmass_leaf) {
	// extkn cannot be 0 (but LAI can be).
	if (lai == 0) {
		return 0;
	}

	// We know Ntot is the integral of N0 * exp(-extkn * LAI) over LAI. Solving
	// for N0 gives:
	return nmass_leaf / ((1.0 - exp(-pft.extkn * lai)) / pft.extkn);
}

/*
Calculate maximum carboxylation rate at 25 degC (μmol m-2 s-1).

@param nmass_leaf: Leaf N content at the top of the canopy (kgN m-2 per green leaf area).
@param pmass_leaf: Leaf P content at the top of the canopy (kgP m-2 per green leaf area).
TODO: check default for pmass_leaf!
*/
double calc_vcmax25(double nmass_leaf, double pmass_leaf, const Pft &pft) {
	// TODO: parameterise these intercepts and slopes.
	return exp(pft.vcmax25_intercept + pft.vcmax25p_slope * log(pmass_leaf * G_PER_KG) + pft.vcmax25n_slope *
							   log(nmass_leaf * G_PER_KG));
}

/*
Calculate maximum electron transport rate at 25 degC (μmol m-2 s-1).

@param nmass_leaf: Leaf N content at the top of the canopy (kgN m-2 per green leaf area).
@param pmass_leaf: Leaf P content at the top of the canopy (kgP m-2 per green leaf area).
*/
double calc_jmax25(double nmass_leaf, double pmass_leaf, const Pft &pft) {
	return exp(pft.jmax25_intercept + pft.jmax25p_slope * log(pmass_leaf * G_PER_KG) + pft.jmax25n_slope *
							   log(nmass_leaf * G_PER_KG));
}

/*
Calculate maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1).

@param nmass_leaf: Leaf N content at the top of the canopy (kgN m-2 per green leaf area).
*/
double calc_vpmax25(double nmass_leaf, const Pft &pft) {
	return exp(pft.vpmax25_intercept + pft.vpmax25n_slope * log(nmass_leaf * G_PER_KG));
}

/*
Calculate Vcmax and Jmax from leaf N following Ellsworth et al. 2022 (page 8) for C3 grasses
and preliminary data from Yuzhen Fan for C4 grasses.
*/
void vcmax_from_N(Individual &indiv) {

	const Pft &pft = indiv.pft;
	double nmass_leaf0 = 0.0;

	// Avoid log(0)
	if (!is_zero(indiv.get_nmass_leaf())) {

		// Calculate leaf N content at the top of the canopy (kgN m-2 per green leaf area)
		nmass_leaf0 = calc_nmass_leaf0(pft, indiv.lai_today(), indiv.get_nmass_leaf());

		assert(positive(nmass_leaf0));

		if (indiv.pft.pathway == C3) {

			// leaf P in first fully expanded leaves per unit leaf area (kgP m-2 leaf) 
			// pmass_leaf is in kgP m-2 leaf area
			// nmass_leaf is in kgN m-2 leaf area 
			double pmass_leaf0 = nmass_leaf0 / C3_LEAF_NTOP;

			indiv.ps_subdaily.vcmax25 = calc_vcmax25(nmass_leaf0, pmass_leaf0, pft);
			indiv.ps_subdaily.jmax25 = calc_jmax25(nmass_leaf0, pmass_leaf0, pft);
		} 
		else if (indiv.pft.pathway == C4) {

			indiv.ps_subdaily.vpmax25 = indiv.pft.vpmax_scal * calc_vpmax25(nmass_leaf0, pft);
		} else {
			fail("Unknown photosynthetic pathway");
		}
	}
	else {
		// This may happen on first day after establishment when using the
		// annual (trunk) allocation, because individuals are given no N at
		// establishment. Nitrogen uptake happens daily though, so on the
		// next day, nmass_leaf should be nonzero.
		indiv.ps_subdaily.vcmax25 = 0.0;
		indiv.ps_subdaily.jmax25  = 0.0;
		indiv.ps_subdaily.vpmax25 = 0.0;
	}
}

/*
Store vcmax25, jmax25, and vpmax25 values for all individuals. The values will
be taken either from PFT parameters, or calculated from leaf N and P content,
based on simulation configuration.
*/
void store_vcmax_values(Vegetation &vegetation) {
	for (size_t i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		if (indiv.pft.ifcalcvcmax) {
			// Calculate photosynthetic capacity from leaf N
			vcmax_from_N(indiv);
		} else {
			// Use parameter value.
			vcmax_from_parameters(indiv);
		}
	}
}

/**
* Calculate water stress factor [0-1] in the root zone. This is calculated as a
* weighted mean of wcont / wcont_crit (clamped to [0, 1]) across soil layers.
*
* @param soil: Soil object containing soil layer information.
* @param weight: Array of weights for each soil layer, which could be based on root distribution or inverse root resistance depending on the configuration.
* @param max_depth: Maximum depth (mm).
*/
double calc_water_stress_dave(const Soil &soil, const double weight[NSOILLAYERS_MAX], double max_depth = -1) {
	if (max_depth < 0) {
		max_depth = soil.depth;
	}

	if (max_depth > soil.depth) {
		fail("calc_water_stress_dave(): max_depth (%d) cannot be greater than soil depth (%d)",
			 max_depth, soil.depth);
	}

	double wstress = 0;
	double total_weight = 0;
	double depth_cum = 0;
	for (size_t i = 0; i < soil.nsoillayers && depth_cum < max_depth; i++) {
		double stress_base = soil.soillayers[i].wcont() / wcont_crit;

		// Clamp to [0, 1].
		stress_base = max(min(stress_base, 1.0), 0.0);

		wstress += stress_base * weight[i];

		// Root distribution should add up to 1 but summed here just in case.
		total_weight += weight[i];

		depth_cum += soil.soillayers[i].thickness;
	}

	// If none of the weighted layers have any weight (e.g. all are dry/
	// rootless), there's no water available to this individual from any of
	// them - maximally stressed.
	if (is_zero(total_weight)) {
		return 0.0;
	}

	// Normalise by sum of weights.
	wstress /= total_weight;
	assert(wstress >= 0.0);
	assert(wstress <= 1.0);
	return wstress;
}

/*
DAVE water stress scalar routine.

This calculates a water stress scalar for each individual and each PFT.
*/
void water_stress_dave(Patch &patch) {
	const Soil &soil = patch.soil;

	// Iterate over individuals in patch.
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		Individual &indiv = patch.vegetation[i];

		// Water stress factor [0-1] in the root zone.
		// With zhou enabled, weight by inverse of root resistance.
		// With zhou disabled, weight by root distribution.
		double weight[NSOILLAYERS_MAX];
		for (size_t j = 0; j < soil.nsoillayers; j++) {
			weight[j] = ifzhou ? 1.0 / indiv.Rs_r[j] : indiv.rootdist[j];
		}
		double wstress = calc_water_stress_dave(soil, weight);

		// Mean daily water stress value.
		indiv.set_wscal(wstress);
	}

	// Calculate a PFT-level water stress. This is only really used for the
	// forest floor photosynthesis, which is used for establishment. Since
	// seedlings don't really have access to deep soil layers, we only use the
	// top few layers' soil moisture.
	for (int i = 0; i < pftlist.nobj; i++) {
		Patchpft &ppft = patch.pft[i];
		double wstress = calc_water_stress_dave(soil, ppft.pft.input_rootdist,
												FOREST_FLOOR_WSCAL_DEPTH);
		ppft.set_wscal(wstress);
	}
}

/*
Acclimation of vcmax/jmax to temperature (Kumarathunge et al, 2018).
*/
void temp_acclimation_kummarathunge(double tgrowth, double thome, double &ha_v,
									double &ha_j, double &ds_v, double &ds_j) {

	// eav = ha_v, eaj = ha_j, Δsv = ds_v, Δsj = ds_j
	ha_v = 42.6 + 1.14 * tgrowth;
	ha_j = 40.71;
	ds_v = 645.13 - 0.38 * tgrowth;
	ds_j = 658.77 - 0.84 * thome - 0.52 * (tgrowth - thome);
	// jvr = 2.56 - 0.0375 * thome - 0.0202 * (tgrowth - thome)
	// jmax25 = vcmax25 * jvr
}

/*
Acclimation of vcmax/jmax to temperature (Kumarathunge et al, 2018).
*/
void get_temp_acclimation_parameters(const Climate &climate, const Pft &pft,
									 double &ha_v, double &ha_j, double &ds_v,
									 double &ds_j) {

	if (ifvcmaxtempacclimation) {
		// Tgrowth, mean air temperature of preceding 30 d.
		double tgrowth = climate.get_tgrowth();
		temp_acclimation_kummarathunge(tgrowth, climate.thome, ha_v, ha_j, ds_v, ds_j);
	} else {
		// Use PFT constants.
		ha_v = pft.ha_v;
		ha_j = pft.ha_j;
		ds_v = pft.ds_v;
		ds_j = pft.ds_j;
	}
}

/*
Calculate stomatal term (Medlyn et al. 2011).
X needs to represent gs for CO2, hence no 1.6 factor.
*/
double stomatal_conductance_term(double wscal, double g1, double qs, double vpd, double ca) {
	return divide(1.0 + divide(pow(wscal, qs) * g1, sqrt(max(vpd, 0.1))), ca);
}

/*
Residual stomatal conductance term.
Note: g0 is for H20, here conversion to conductances to CO2.
TODO: do we need a different qs here?
*/
double stomatal_conductance_residual(double g0, double qs, double wscal) {
	return divide(max(0.0, g0 * pow(wscal, qs)), DWDC);
}

/*
System of equations to solve coupled photosynthesis-stomatal conductance scheme.
Wang & Leuning 1998
*/
double photosyn_solver(double g0, double X, double gamma, double beta, 
					   double gamstar, double Rd, double Ca) {
	// km = beta, gamma = vcmax
	double b2 = g0 + X * (gamma - Rd);
	double b1 = (1.0 - Ca * X) * (gamma - Rd) + g0 * (beta - Ca) - X * (gamma * gamstar + beta * Rd);
	double b0 = -(1.0 - Ca * X) * (gamma * gamstar + beta * Rd) - g0 * beta * Ca;

	double Ci;
	if (!solve_quadratic(b2, b1, b0, Ci, false)) {
		log_error("photosyn_solver: No solution for quadratic Ci. Returning An = -Rd to avoid crash. g0=%g, X=%g, gamma=%g, beta=%g, gamstar=%g, Rd=%g, Ca=%g",
				  g0, X, gamma, beta, gamstar, Rd, Ca);
		return -Rd;
	}

	if (is_zero(Ci + beta)) {
		log_warning("Ci + beta is zero in photosyn_solver. Returning An = -Rd to avoid crash.");
		return -Rd;
	}

	double solution = gamma * divide(Ci - gamstar, Ci + beta);

	if (solution < 0.0) {
		// Logging as a recoverable error because this should never happen.
		log_error("Negative solution (gpp) in photosyn_solver: %g. A gross value of 0 will be used instead.", solution);
		solution = 0;
	}

	double An = solution - Rd;

	return An;  // net photosynthesis (mol m-2 s-1)
}

/*
Temperature response function following a modified Arrhenius equation (Johnson
et al. 1942).

@param k25 process rate at 25 degC
@param Tl leaf temperature (deg C)
@param Ha activation energy (kJ mol-1)
@param Hd deactivation energy (kJ mol-1)
@param dS entropy term (J mol-1 K-1)
*/
double Arrhenius_temp_response(double k25, double Tl, double Ha, double Hd, double dS) {
	// Define reference temperature and convert units
	double Tref = 25.0 + K2degC; // Reference Temperature (25 degC in K)
    Tl += K2degC;
    Ha *= J_PER_KJ; 
	Hd *= J_PER_KJ; 

	double kTl = k25 * exp(Ha * (Tl - Tref)/(Tref * R_gas * Tl)) * 
				      (1.0 + exp((Tref * dS - Hd)/(Tref * R_gas))) / 
				      (1.0 + exp((Tl * dS - Hd)/(Tl * R_gas)));

	return kTl;  // process rate at leaf Temperature Tl
}

/*
Temperature response of leaf mitochondrial respiration.

Eq. 2 in Heskel et al. 2016, PNAS.

Value of parameter a was calculated to give a value of 1 for Tl = 25degC.

@param tl leaf temperature (deg C)
*/
double respiration_temp_response(double tl) {
	double kRl = exp(-2.2175 + 0.1012 * tl - 0.0005 * pow(tl, 2));
	return kRl; // value of leaf respiration relative to its value at 25 deg C
}

/*
Light inhibition of leaf mitochondrial respiration. Assumes 30% reduction under
light as in Huntingford et al. 2017, Nature Communications.

@param par incident photosynthetically active radiation (W m-2)
*/
double light_inhibition(double par) {
	double inhib = par > 2.0 ? 0.7 : 1.0;
	return inhib; // inhibition of leaf respiration by light (0-1)
}

/*
Calculate electron transport rate from PAR.

von Caemmerer 2013 (Plant, Cell and Environment)

@param apar  = Absorbed photosynthetically active radiation (W m-2)
@param jmax = max. electron transport rate (mol m-2 s-1)
*/
double electron_transport_rate(double apar, double alpha, double jmax) {
	// Photosynthetic photon flux density (mol (photons) m-2 s-1).
	double ppfd = apar * alpha * JTOMOL;

	double a = THETA;
	double b = -(ppfd + jmax);
	double c = ppfd * jmax;

	double j;
	if (!solve_quadratic(a, b, c, j, true)) {
		log_error("electron_transport_rate: No solution for quadratic J. Returning J=0 to avoid crash. apar=%g, alpha=%g, jmax=%g", apar, alpha, jmax);
		return 0;
	}

	return j;   // electron transport rate in mol (electrons) m-2 s-1
}

/*
Scales Vcmax and Jmax under water stress based on Zhou (2014)

@param s_f         Sensitivity parameter indicating steepness of decline in
                   Vcmax with water stress (Zhou 2014)
@param psi_f       Reference value at which Vcmax or Jmax is half its non-water-
                   stresed value (Zhou 2014)
@param predawn_lwp Pre-dawn leaf water potential (MPa)
*/
double pslim_lwp_scalar(double s_f, double psi_f, double predawn_lwp) {
	return (1 + exp(s_f * psi_f)) / (1 + exp(s_f * (psi_f - predawn_lwp)));
}

/*
Get the temperature scalar encapsulating the effect of temperature acclimation
on respiration.

@param tgrowth The growth temperature for the current day, calculated as the
               mean of the previous 30 days' mean daily air temperatures.
@param racclim The acclimation sensitivity coefficient for respiration (1/degC).
*/
double resp_temp_acclimation_scalar(double tgrowth, double racclim) {
	if (!ifresptempacclimation) {
		return 1.0;
	}

    // TODO: add references, document, etc.
	return exp(-racclim * (tgrowth - 15.0));
}

/*
Get the water stress effect on respiration.
*/
double resp_water_scalar(const Individual &indiv) {
	if (indiv.pft.lifeform != TREE) {
		return 1.0;
	}

	double wscal = indiv.get_wscal();
	double scalar = 0.5 * (1.0 + wscal);
	scalar = max(min(scalar, 1.0), 0.0);
	return scalar;
}

/*
No photosynthesis (stomata closed - night time, water stress, etc).

Sets gross assimilation to 0, net assimilation to -Rd, and Ci to Ca.
*/
void no_photosynthesis(const PhotosynthesisInputs &inputs, PhotosynthesisOutputs &outputs, double rd) {
	// Gross assimilation = 0, net assimilation = -Rd, and Ci = Ca.

	outputs.rd = rd;
	outputs.an = -rd;
	outputs.anc = -rd;
	outputs.anj = -rd;
	outputs.vcmax = 0;
	outputs.jmax = 0;
	outputs.anp = 0;

	// Hydraulic stress can take gsc to zero.
	outputs.gsc = inputs.stomata_open ? inputs.gdark : 0;

	// Unused.
	outputs.ci = inputs.ca;
}

/*
C3 photosynthesis (largely Farquhar et al. 1980).

Calculates photosynthesis for a single C3 indiv in a single layer at a given
timestep.

@param inputs: Inputs to this function.
@param outputs: Outputs from this function (passed by reference).
*/
void photosynthesis_dave_c3(const PhotosynthesisInputs &inputs, PhotosynthesisOutputs &outputs) {

	// vcmax25/jmax25 are input in μmol m-2 s-1, but the equations expect them
	// in mol m-2 s-1.
	double vcmax25_mol = inputs.vcmax25 / UMOL_PER_MOL;
	double jmax25_mol  = inputs.jmax25  / UMOL_PER_MOL;
	double resp_scalar = resp_temp_acclimation_scalar(inputs.tgrowth, inputs.racclim);

	double rd = pow(inputs.wscal, qb) * BC3 * vcmax25_mol * inputs.nscal *
			respiration_temp_response(inputs.tleaf) * light_inhibition(inputs.ipar()) * resp_scalar;

	if (is_zero(inputs.ipar())) {
		no_photosynthesis(inputs, outputs, rd);
		return;
	}

	// Michaelis-Menten parameters and photorespiratory CO2 compensation point (mol mol-1).
	double ko = Arrhenius_temp_response(KO25, inputs.tleaf, HA_KO, 0.0, 0.0);
	double kc = Arrhenius_temp_response(KC25, inputs.tleaf, HA_KC, 0.0, 0.0);
	double gamstar = Arrhenius_temp_response(GAMSTAR25, inputs.tleaf, HA_GS, 0.0, 0.0);

	// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax.
	double pslim = inputs.ps_sink_lim * inputs.nscal;

	// apply water stress scaling to vcmax / jmax
	if (ifzhou) {
	// function to scale vcmax and jmax by pre-dawn water potential, see Zhou (2014)
	//	double wp_vmax_scale = pslim_lwp_scalar(inputs.s_f, inputs.psi_f, inputs.predawn_lwp);
	// alternative formula from Xu et al. (2016), New Phyt, https://doi.org/10.1111/nph.14009
	// note we considered a stronger effect (using pg50 rather than tlp as in the paper), but review
	// of literature indicated that about 50% loss of vcmax at stomatal closure was about right.
		double a1 = 8.0;
		//double pg50 = inputs.stom_clos * pow(9.0, -1.0 / a1);
		double wp_vcmax_scale = 1.0 / (1.0 + pow(inputs.lwp / inputs.stom_clos, a1));
		pslim = pslim * wp_vcmax_scale;
	} else {
		pslim = pslim 
		* pow(inputs.wscal, qb);
	}

	double vcmax = pslim * Arrhenius_temp_response(vcmax25_mol, inputs.tleaf, inputs.ha_v, HD, inputs.ds_v);
	double jmax  = pslim * Arrhenius_temp_response(jmax25_mol, inputs.tleaf, inputs.ha_j, HD, inputs.ds_j);

	// 1) Rubisco-limited gross photosynthesis
	double gamma = vcmax;
	double beta  = kc * (1.0 + divide(O2atm, ko));
	double anc   = photosyn_solver(0.0, inputs.X, gamma, beta, gamstar, rd, inputs.ca);

	if (!ifdavehydraulics) {
		// TODO: Discuss how to allow water stress to close stomata further
		if (inputs.gdark > (anc * inputs.X)) {
			// If gdark is greater than gs from the 'normal' gs model,
			// repeat calculation with gdark only.
			anc = photosyn_solver(inputs.gdark, 0.0, gamma, beta, gamstar, rd, inputs.ca);
		}
	}

	// 2) RuBP regeneration-limited gross photosynthesis
	double je = electron_transport_rate(inputs.apar(), inputs.alpha, jmax);

	// assuming demand of 4 electrons per CO2 molecule
	gamma = 0.25 * je;

	beta  = 2.0 * gamstar;
	double anj = photosyn_solver(0.0, inputs.X, gamma, beta, gamstar, rd, inputs.ca);
	if (!ifdavehydraulics) {
		// TODO: Discuss how to allow water stress to close stomata further
		if (inputs.gdark > (anj * inputs.X)) {
			anj = photosyn_solver(inputs.gdark, 0.0, gamma, beta, gamstar, rd, inputs.ca);
		}
	}

	// 3) TPU-limited gross photosynthesis
	// Ap = 3TPU; TPU = 0.167Vcmax => Ap = 0.5 Vcmax
	double anp = 0.5 * vcmax - rd;

	// Net photosynthesis as minimum of 1), 2), and 3) (mol m-2 s-1)
	double an = min(anc, min(anj, anp));
	//assert(an >= 0);

	// Stomatal conductance to CO2 (mol m-2 s-1)
	//gsc = g0 + An * X;
	//double gsc = max(inputs.gdark, an * inputs.X);
	double floor = ifdavehydraulics ? 0.0 : inputs.gdark;
	double gsc = max(floor, an * inputs.X);	// CS removed floor at gdark (since water stress can take gsc to zero)

	// Re-calculate intercellular CO2 concentration (mol mol-1)
	double ci;

	if (!inputs.stomata_open) {
		gsc = 0.0;
		an = -rd;
	}

	if (is_zero(gsc)) {
		ci = inputs.ca;
	} else {
		ci = inputs.ca - divide(an, gsc);
	}

	outputs.rd = rd;
	outputs.an = an;
	outputs.anc = anc;
	outputs.anj = anj;
	outputs.anp = anp;
	outputs.gsc = gsc;
	outputs.ci = ci;
	outputs.vcmax = vcmax;
	outputs.jmax = jmax;
}

/*
Approximate solution of Ac for C4 photosynthesis, following the equations in von
Caemmerer (2021).
*/
double calc_c4_ac_approx(const C4Biochemistry &biochemistry,
						 const PhotosynthesisInputs &inputs,
						 const C4IterationState &state) {
	double vcmax = biochemistry.vcmax;
	double rd = biochemistry.rd;
	double rm = biochemistry.rm;
	double vp = state.vp;
	double ci = state.ci;

	// Equation 25
	return min(vp + gbs * ci - rm, vcmax - rd);
}

/*
Full solution for Ac for C4 photosynthesis, following the equations in von
Caemmerer (2021).
*/
double calc_c4_ac_full(const C4Biochemistry &biochemistry,
	                   const PhotosynthesisInputs &inputs,
				       const C4IterationState &state) {
    double kc = biochemistry.kc;
	double ko = biochemistry.ko;
	double vp = state.vp;
	double rm = biochemistry.rm;
	double ci = state.ci;
	double vcmax = biochemistry.vcmax;
	double rd = biochemistry.rd;
	double a0 = biochemistry.a0;
	double lgamstar = biochemistry.lgamstar;

	// kc * (1 + om / ko)
	double km = biochemistry.km;

	// Equations 22-24
	double a_c = alpha_PSII / a0 * kc / ko;
	double b_c = -((vp - rm + gbs * ci) + (vcmax - rd) + gbs * km + alpha_PSII / a0
			* (lgamstar * vcmax + rd * kc / ko));
	double c_c = (vcmax - rd) * (vp - rm + gbs * ci) -
			(vcmax * gbs * lgamstar * O2atm + rd * gbs * km);

	double anc;
	if (solve_quadratic(a_c, b_c, c_c, anc, true)) {
		if (anc < -rd || anc > 2.0 * vcmax) {
			log_warning("Implausible Anc value with the lower root! anc=%f, rd=%f, vcmax=%f",
				 anc, rd, vcmax);
		}
		return anc;
	}
	log_warning("No solution to C4 photosynthesis quadratic equation for Ac. A value of -Rd (%g) will be used for Ac. a_c: %g, b_c: %g, c_c: %g",
				rd, a_c, b_c, c_c);
	return -rd;
}

/*
Solve C4 photosynthesis for Ac, using either the approximate or full solution
based on the simulation configuration.
*/
double calc_c4_ac(const C4Biochemistry &biochemistry,
	              const PhotosynthesisInputs &inputs,
				  const C4IterationState &state) {
	if (ifapproximate_C4) {
		return calc_c4_ac_approx(biochemistry, inputs, state);
	}
	return calc_c4_ac_full(biochemistry, inputs, state);
}

/*
Approximate solution of Aj for C4 photosynthesis, following the equations in von
Caemmerer (2021).
*/
double calc_c4_aj_approx(const C4Biochemistry &biochemistry,
						 const PhotosynthesisInputs &inputs,
						 const C4IterationState &state) {
	// electron transport rate (same as C3) (mol (electrons) m-2 s-1)
	double je = electron_transport_rate(inputs.apar(), inputs.alpha, biochemistry.jmax);

	double rm = biochemistry.rm;
	double ci = state.ci;
	double rd = biochemistry.rd;

	// Equation 44
	return min(((divide(zj,2.0) * xj * je) - rm + gbs * ci),
					((divide(zj,3.0) * (1.0 - xj) * je) - rd));
}

/*
Full solution for Aj for C4 photosynthesis, following the equations in von
Caemmerer (2021).
*/
double calc_c4_aj_full(const C4Biochemistry &biochemistry,
					   const PhotosynthesisInputs &inputs,
					   const C4IterationState &state) {

	// electron transport rate (same as C3) (mol (electrons) m-2 s-1)
	double je = electron_transport_rate(inputs.apar(), inputs.alpha, biochemistry.jmax);

	double rm = biochemistry.rm;
	double ci = state.ci;
	double rd = biochemistry.rd;
	double a0 = biochemistry.a0;
	double lgamstar = biochemistry.lgamstar;

	// Equations 41-43
	double a_j = 1.0 - (7.0 * lgamstar * alpha_PSII) / (3.0 * a0);
	double b_j = -((divide(zj,2.0) * xj * je - rm + gbs * ci) + (divide(zj,3.0) * (1.0 - xj) * je - rd) +
			gbs * (divide((7.0 * lgamstar * O2atm),3.0)) + divide((alpha_PSII * lgamstar),a0) *
			(divide(zj,3.0) * (1.0 - xj) * je + divide((7.0 * rd),3.0)));
	double c_j = ((divide(zj,2.0) * xj * je) - rm + gbs * ci) * (divide(zj,3.0) * (1.0 - xj) * je - rd) -
			gbs * lgamstar * O2atm * (divide(zj,3.0) * (1.0 - xj) * je + divide((7.0 * rd),3.0));

	double anj;
	if (solve_quadratic(a_j, b_j, c_j, anj, false)) {
		if (anj < -rd || anj > biochemistry.jmax) {
			log_warning("Implausible Anj value with the higher root! anj=%f, rd=%f, jmax=%f",
						anj, rd, biochemistry.jmax);
		}
		return anj;
	}

	log_warning("No solution to C4 photosynthesis quadratic equation for Aj. A value of -Rd (%g) will be used for Aj. a_j: %g, b_j: %g, c_j: %g",
				rd, a_j, b_j, c_j);
	return -rd;
}

/*
Solve C4 photosynthesis for Aj, using either the approximate or full solution
based on the simulation configuration.
*/
double calc_c4_aj(const C4Biochemistry &biochemistry,
				  const PhotosynthesisInputs &inputs,
				  const C4IterationState &state) {
	if (ifapproximate_C4) {
		return calc_c4_aj_approx(biochemistry, inputs, state);
	}
	return calc_c4_aj_full(biochemistry, inputs, state);
}

/*
Solve C4 photosynthesis for a single C4 indiv in a single layer at a given
timestep, following the equations in von Caemmerer (2021).
*/
C4SolverResult solve_c4_photosynthesis(const C4Biochemistry &rates,
									   const PhotosynthesisInputs &inputs) {

	// C4/day time
	double anc = 0;
	double anj = 0;
	double an  = 0;
	double gsc = 0;

	double vpmax = rates.vpmax;
	double vpr = rates.vpr;
	double kp = rates.kp;

	// Iteration to find An, gs, and ci
	double ci = 0.5 * inputs.ca;

	double ci_diff;
	int it;
	for (it = 0; it < NIT_MAX; it++) {
		if (it == NIT_MAX - 1) {
			// No solution found. Use initial ci guess on last iteration.
			ci = 0.5 * inputs.ca;
		}

		// PEP carboxylation rate
		double vp = min(divide(ci * vpmax, ci + kp), vpr);

		C4IterationState state;
		state.ci = ci;
		state.vp = vp;

		// Enzyme-limited net photosynthesis rate (mol m-2 s-1)
		anc = calc_c4_ac(rates, inputs, state);

		// Light-limited net photosynthesis rate (mol m-2 s-1)
		anj = calc_c4_aj(rates, inputs, state);

		// Net photosynthesis rate (mol m-2 s-1)
		an = min(anc, anj);

		// Stomatal conductance to CO2 (mol m-2 s-1)
		// gsc = g0 + An * X;
		// gsc = max(inputs.gdark, an * inputs.X);
		double floor = ifdavehydraulics ? 0.0 : inputs.gdark;

		// CS removed floor at gdark since water stress can take gsc to zero.
		gsc = max(floor, an * inputs.X);

		// Intercellular CO2 concentration (mol mol-1)
		if (!inputs.stomata_open) {
			gsc = 0.0;
		}

		double ci_new = is_zero(gsc) ? inputs.ca : inputs.ca - divide(an, gsc);

		if (it != NIT_MAX - 1) {
			// Don't update ci_diff on final iteration, because its value will
			// be misleading (as ci was reset to its initial value for this
			// iteration).
			ci_diff = ci - ci_new;

			// If we can't find a solution, use initial guess on last iteration.
			ci -= ci_diff / 2.0;
		}

		// Cease iteration if ci has converged to within CI_DIFF_MAX.
		if (it > 0 && abs(ci_diff) <= CI_DIFF_MAX) {
			break;
		}
	} // end solver

	C4SolverResult state;
	state.an = an;
	state.ac = anc;
	state.aj = anj;
	state.gsc = gsc;
	state.ci = ci;
	state.iterations = it;
	state.ci_diff = ci_diff;
	state.converged = (it < NIT_MAX - 1);
	return state;
}

/*
Get the C4 biochemical parameters for a given PFT/timestep/layer.
*/
C4Biochemistry c4_biochemistry(const PhotosynthesisInputs &inputs) {
	C4Biochemistry result;

	// Parameters at 25 °C (mol m-2 s-1).
	double vpmax25_mol = inputs.vpmax25 / UMOL_PER_MOL;
	double vcmax25_c4 = rVC4 * vpmax25_mol;
	double jmax25_c4  = rJC4 * vcmax25_c4;
	result.vpr        = rVPr * vcmax25_c4;
	double resp_scalar = resp_temp_acclimation_scalar(inputs.tgrowth, inputs.racclim);

	// Daytime (mitochondrial) leaf respiration (mol m-2 s-1).
	result.rd = pow(inputs.wscal, qb) * BC4 * vcmax25_c4 * inputs.nscal *
			    respiration_temp_response(inputs.tleaf) * light_inhibition(inputs.ipar()) * resp_scalar;

	// Rd in the mesophyll cells.
	result.rm = 0.5 * result.rd;

	// Rd in the bundle-sheath cells.
	result.rs = result.rd - result.rm;

	// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax
	double pslim = inputs.ps_sink_lim * inputs.nscal;

	// apply water stress scaling to vcmax / jmax
	if (ifzhou) {
		// function to scale vcmax and jmax by pre-dawn water potential, see Zhou (2014)
		//	double wp_vmax_scale = pslim_lwp_scalar(inputs.s_f, inputs.psi_f, inputs.predawn_lwp);
	// alternative formula from Xu et al. (2016), New Phyt, https://doi.org/10.1111/nph.14009
		double a1 = 10.0;
		// double pg50 = inputs.stom_clos * pow(9.0, -1.0 / a1);
		double wp_vcmax_scale = 1.0 / (1.0 + pow(inputs.lwp / inputs.stom_clos, a1));
		pslim = pslim * wp_vcmax_scale;
	} else {
		pslim = pslim * pow(inputs.wscal, qb);
	}

	result.vcmax = pslim * Arrhenius_temp_response(vcmax25_c4, inputs.tleaf, HA_V4, HD, DS_V4);

	// Note: Jmax temperature response has been refitted to an Arrhenius equation
	// Jmax  = pslim * (Jmax25_C4 * exp(-(pow((temp - 43.0)/26.0,2))));
	result.jmax  = pslim * Arrhenius_temp_response(jmax25_c4, inputs.tleaf, HA_J4, HD, DS_J4);

	// Temperature dependency of Vpmax and Michaelis-Menten constants
	result.vpmax = pslim * Arrhenius_temp_response(vpmax25_mol, inputs.tleaf, HA_VP, HD, DS_VP);

	result.kp = Arrhenius_temp_response(KP25, inputs.tleaf, HA_KP, 0.0, 0.0);
	result.kc = Arrhenius_temp_response(KC25_C4, inputs.tleaf, HA_KC_C4, 0.0, 0.0);
	result.ko = Arrhenius_temp_response(KO25_C4, inputs.tleaf, HA_KO_C4, 0.0, 0.0);
	result.km = result.kc * (1.0 + divide(O2atm, result.ko));
	result.a0 = Arrhenius_temp_response(a025, inputs.tleaf, HA_a0, 0.0, 0.0);
	result.lgamstar = Arrhenius_temp_response(lgamstar25, inputs.tleaf, HA_lgamstar, 0.0, 0.0);

	return result;
}

/*
C4 photosynthesis (von Caemmerer 2000).

Calculates photosynthesis for a single C4 indiv in a single layer at a given
timestep.

@param inputs: Inputs to this function.
@param outputs: Outputs from this function (passed by reference).
*/
void photosynthesis_dave_c4(const PhotosynthesisInputs &inputs, PhotosynthesisOutputs &outputs) {
	C4Biochemistry biochemistry = c4_biochemistry(inputs);

	if (is_zero(inputs.ipar()) || !inputs.stomata_open) {
		no_photosynthesis(inputs, outputs, biochemistry.rd);
		return;
	}

	C4SolverResult result = solve_c4_photosynthesis(biochemistry, inputs);

	if (result.converged) {
		log_diagnostic("C4 photosynthesis solver converged after %d iterations (ci_diff = %g).",
					   result.iterations, result.ci_diff);
	} else {
		log_warning("C4 photosynthesis solver did not converge after %d iterations (ci_diff = %g).",
					NIT_MAX, result.ci_diff);
	}

	outputs.rd = biochemistry.rd;
	outputs.an = result.an;
	outputs.anc = result.ac;
	outputs.anj = result.aj;
	outputs.gsc = result.gsc;
	outputs.ci = result.ci;
	outputs.vcmax = biochemistry.vcmax;
	outputs.jmax = biochemistry.jmax;
}

// Populate photosynthesis inputs struct and return the g1 value to be used.
PhotosynthesisDiagnostics get_photosynthesis_inputs(const CohortState &cohort,
		const Pft &pft, const Climate &climate, const Patch &patch,
		const Day &day, double lai_cum, PhotosynthesisInputs &inputs) {

	// Atmospheric CO2 concentration (mol mol-1)
	// TODO: check if annual global means would make more sense 
	inputs.ca = get_co2(patch, climate, pft) / UMOL_PER_MOL;  // mol mol-1
	inputs.fphen_activity = cohort.get_fphen_activity();
	inputs.ps_sink_lim = cohort.get_ps_sink_lim();

	inputs.vcmax25 = cohort.get_vcmax25();
	inputs.vpmax25 = cohort.get_vpmax25();
	inputs.jmax25 = cohort.get_jmax25();
	// wscal based on root resistance if using hydraulics, otherwise wcont.
	inputs.wscal = cohort.get_wscal();
	inputs.s_f = pft.s_f;
	inputs.psi_f = pft.psi_f;
	inputs.tgrowth = climate.get_tgrowth();
	inputs.racclim = pft.racclim;

	get_temp_acclimation_parameters(climate, pft, inputs.ha_v, inputs.ha_j,
									inputs.ds_v, inputs.ds_j);

	// Air temperature in °C.
	inputs.tleaf  = climate.temps[day.period];

	// VPD in kPa.
	inputs.vpd   = climate.vpds[day.period];

	inputs.stomata_open = cohort.get_stomata_open();

	// lwp for Xu (2016) adjustment to vcmax, predawn_lwp for Zhou (2014) stomatal closure model
	inputs.lwp = cohort.get_lwp();
	inputs.predawn_lwp = cohort.get_predawn_lwp();

	inputs.stom_clos = cohort.get_stom_clos();
	inputs.tlp = pft.tlp;

	// CS, Testing adjustment to gcutic with temperature - based on
	// Cochard 2021 Eq 13 & 14,
	// https://peercommunityjournal.org/item/10.24072/pcjournal.45.pdf

	// Calculate gcutic at phase change T for use in equation below.
	double gcutic_Tp = pft.gcutic * pow(Q10A, (TP - 25) / 10);

	// This should all be based on leaf temperature, but using air
	// temperature as an approximation.
	double gcutic_adj = 0.0;
	if (climate.temps[day.period] < TP) {
		// Erick's measurements start at 25, so changed 20 to 25 in eqns
		gcutic_adj = pft.gcutic * pow(Q10A, (climate.temps[day.period] - 25) / 10);
	} else {
		gcutic_adj = gcutic_Tp * pow(Q10B, (climate.temps[day.period] - TP) / 10);
	}

	PhotosynthesisDiagnostics diagnostics;
	diagnostics.gcutic_adj = gcutic_adj;
	diagnostics.g0 = pft.g0;
	diagnostics.g1 = pft.g1;

	if (ifzhou) {
		// Stomatal term (Medlyn et al. 2011), with g1 as a declining
		// exponential function of leaf water potential (Zhou et al. 2014,
		// eqn 3: g1 = a*exp(b*psi)).

		// We set parameter b1z to the value that results in 90% reduction in
		// g1 at stom_clos
		double b1z = log(0.1) / (pft.stom_clos + 0.3);

		// for water
		double rescale = exp(b1z * (inputs.predawn_lwp + 0.3));

		double g1 = pft.g1x * rescale;

		// for water
		// rescale similar to g1, but don't allow values bigger than assigned g0
		double g0 = min(pft.g0, pft.g0 * rescale);

		// Don't allow g0 to fall below cuticular conductance.
		g0 = max(g0, gcutic_adj);

		// X needs to represent gs for CO2, hence no 1.6 factor
		inputs.X = divide(1.0 + divide(g1, sqrt(max(inputs.vpd, 0.1))), inputs.ca);

		// Residual stomatal conductance term
		// Note: g0 is for H20, here conversion to conductances to CO2
		// TODO: do we need a different qs here?
		inputs.gdark = divide(max(0.0, g0), DWDC);

		diagnostics.g1 = g1;
		diagnostics.g0 = g0;
	} else {
		// Stomatal term (Medlyn et al. 2011)
		// X needs to represent gs for CO2, hence no 1.6 factor
		inputs.X = stomatal_conductance_term(inputs.wscal, pft.g1, pft.qs, inputs.vpd, inputs.ca);

		// Residual stomatal conductance term
		// Note: g0 is for H20, here conversion to conductances to CO2
		// TODO: do we need a different qs here?
		inputs.gdark = stomatal_conductance_residual(pft.g0, pft.qs, inputs.wscal);
	}

	// LAI of this individual in this layer.
	// Note: this takes phen into account for indivs running trunk
	// allocation, because fpar() looks at indiv.lai_today().
	inputs.lai = cohort.get_lai();
	inputs.lai_cum = lai_cum;

	// We can actually have a canopy layer with 0 indiv-level LAI
	// but non-zero fipar. This is because fipar is not a truly
	// indiv-level property, as it is calculated using canopy-level
	// layered LAI.
	//
	// Question: can we skip this layer if lai_layer is 0?
	assert(inputs.lai >= 0);

	inputs.par = max(0.0, climate.pars[day.period]);
	inputs.fipar = cohort.get_fipar();
	// TODO: alpha only used in photosynthesis
	inputs.alpha = pft.alpha;
	inputs.absorptance_nir = pft.absorptance_nir;
	inputs.absorptance_par = pft.absorptance_par;
	inputs.leaf_width = pft.leaf_width;

	// Account for the nitrogen gradient within the canopy.
	inputs.nscal = max(1.0e-12, exp(-pft.extkn * lai_cum));

	return diagnostics;
}

/*
Calculate leaf temperature based on energy balance.

@param tleaf     Leaf temperature (°C).
@param tair      Air temperature (°C).
@param gs		 Stomatal conductance to water vapour (mol m-2 s-1).
@param ppfd		 Photosynthetic photon flux density (aka incident PAR) (umol m-2 s-1).
@param vpd		 Vapour pressure deficit (kPa).
@param patm		 Atmospheric pressure (kPa).
@param wind		 Wind speed in this layer (m s-1)
@param wleaf	 Leaf width (m) - PFT specific.
@param stomatal_ratio Captures the effect of leaves having stomata on one vs both sides.
@param abs_par	 Fraction of PAR absorbed by leaf area (0-1).
@param abs_nir	 Fraction of near-infrared radiation absorbed by leaf area (0-1).
@param fipar     Fraction of PAR incident on the leaf (0-1).

@return Leaf temperature (°C).
*/
EnergyBalanceResult energy_balance_dave(double tleaf, double tair, double gs,
	                                    double ppfd, double vpd, double patm,
										double wind, double wleaf,
					                    double stomatal_ratio, double abs_par,
										double abs_nir, double fipar) {
	if (gs == 0) {
		fail("energy_balance_dave(): gs=0");
	}

	// TODO: consolidate with daveinput into a single global parameter?
	const esat_method_t m = BUCK_1981;

	double tair_k = tair + K2degC;

	// Short calculations.
	double air_density = patm * 1000.0 / (287.058 * tair_k); // kg m-3 ? 
	double lhv = (H2OLV0 - 2365.0 * tair) * MOL2G_H2O * KG_PER_G; // J mol-1
	double slope = 1000.0 * (esat(tair + 0.1, m) - esat(tair, m)) / 0.1; // Pa K-1

	// Units: 1 J/m3 = 1 Pa. So Pa / (J mol-1 K-1) / K = Pa mol / J = J m-3 mol /J = mol m-3
	double cmolar = patm * 1000.0 / (R_gas * tair_k); // mol m-3

	// Conductances - all in mol m-2 s-1.
	double gradiation = 4.0 * BOLTZ * pow(tair_k, 3.0) * EMISSIVITY / (CP_AIR * AIRMA) * fipar;
	double gbhforced = 0.003 * sqrt(wind / wleaf) * cmolar ;
	double grashof = 1.6e+08 * fabs(tleaf - tair) * pow(wleaf, 3.0);
	double gbhfree = 0.5 * DHEAT * pow(grashof, 0.25) / wleaf * cmolar;
	double gbh = 2.0 * (gbhfree + gbhforced);
	// double gbhr = gbh + 2.0 * gradiation;
	double gbw = stomatal_ratio * 1.075 * gbh;
	double gw = gs * gbw / (gs + gbw);

	// Radiation.
	// double Rlongup = EMISSIVITY * BOLTZ * pow(Tk(tleaf), 4);
	double rsol = (abs_par + abs_nir) * ppfd / UMOL_PER_J; // J m-2 s-1 - absorbed SW radiation = 2 * APAR
	// double Rnet = rsol - Rlongup;
	double es = esat(tair, m);
	double ea = 1000.0 * (es - vpd); // In Pa. If this is negative we have a problem
	if (ea < 0) {
		log_error("energy_balance_dave(): ea < 0, this means VPD exceeds saturation vapour pressure, which shouldn't happen. esat=%f kPa, vpd=%f kPa, tair=%f °C",
			        es, vpd, tair);
		ea = 0.0;
	}

	double ema = 0.642 * pow(ea / tair_k, 1.0 / 7.0);
	double rnetiso = rsol - (EMISSIVITY - ema) * BOLTZ * pow(tair_k, 4.0) * fipar;

	// Penman-Monteith equation
	double gamma = CP_AIR * AIRMA * patm * 1000.0 / lhv;

	double y = 1.0 / (1.0 + gradiation / gbh);

	// Transpiration rate (mol H2O m-2 s-1) : in place of gs * VPD
	double et = (1.0 / lhv) * (slope * y * rnetiso + 1000.0 * vpd * gbh *
							   CP_AIR * AIRMA) / (slope * y + gamma * gbh / gw);

	// Calculate leaf temperature.
	double lambda_et = lhv * et;
	double h2 = y * (rnetiso - lambda_et);
	// double H = -CP_AIR * air_density * (gbh / cmolar) * (tair - tleaf);
	double tleaf2 = tair + h2 / (CP_AIR * air_density * (gbh / cmolar));
	assert(tleaf2 == tleaf2); // check for NaNs

	EnergyBalanceResult result;
	result.tleaf = tleaf2;
	result.transp = et;
	result.rnet = rnetiso;
	return result;
}

void photosynthesis_dave(const PhotosynthesisInputs &inputs,
						 PhotosynthesisOutputs &outputs,
						 pathwaytype pathway) {
	if (pathway == C3) {
		photosynthesis_dave_c3(inputs, outputs);
	} else {
		photosynthesis_dave_c4(inputs, outputs);
	}
}

double get_gsw(double gsc, double gcutic) {
	double gsw = gsc * DWDC; // convert from gs for CO2 to gs for H2O
	if (ifdavehydraulics && gsw < gcutic) {
		// Don't allow stomatal conductance to water vapour to fall below
		// cuticular conductance.
		gsw = gcutic;
	}
	return gsw;
}

EnergyBalance solve_energy_balance_dave(PhotosynthesisInputs &inputs,
							   const ClimateInputs &climate_inputs,
						       PhotosynthesisOutputs &outputs,
							   const PhotosynthesisDiagnostics &diagnostics,
						       pathwaytype pathway,
							   size_t niter_max,
							   double epsilon) {
    if (niter_max == 0) {
		fail("solve_energy_balance_dave(): niter_max must be greater than 0.");
	}

    // TODO: might experiment with changing this at some point.
    const double stomatal_ratio = 1.0;

	// Initialise guess: leaf temperature = air temperature.
	inputs.tleaf = climate_inputs.get_tair();

	bool converged = false;

	// Guaranteed to be set.
	size_t i;
	double tleaf_diff;
	EnergyBalanceResult ebal;
    for (i = 0; i < niter_max && !converged; i++) {
		outputs.clear();
		photosynthesis_dave(inputs, outputs, pathway);

		double tleaf_old = inputs.tleaf;
		double gsw = get_gsw(outputs.gsc, diagnostics.gcutic_adj);
		double wind = max(climate_inputs.get_wind(), 0.01);
		wind *= lambertbeer(inputs.lai_cum);
		double ppfd = inputs.ipar() * JTOMOL * UMOL_PER_MOL;

		ebal = energy_balance_dave(inputs.tleaf,
										   climate_inputs.get_tair(),
									       gsw,
										   ppfd,
									       climate_inputs.get_vpd(),
									       climate_inputs.get_patm(),
									       wind,
									       inputs.leaf_width,
										   stomatal_ratio,
									       inputs.absorptance_par,
										   inputs.absorptance_nir,
										   inputs.fipar);
		inputs.tleaf = ebal.tleaf;
        tleaf_diff = fabs(inputs.tleaf - tleaf_old);
		converged = tleaf_diff <= epsilon;
	}

	// This is not the same as the inputs.tleaf used to generate these
	// photosynthesis outputs, but rather it's the tleaf that results from them.
	// This reproduces old no-energy-balance code when niter_max = 1.
	outputs.tleaf = inputs.tleaf;
	EnergyBalance result(inputs.tleaf, i, tleaf_diff, converged, ebal.transp,
		                 ebal.rnet);
	return result;
}

/*
Coupled model of photosynthesis and stomatal conductance for C3 and C4 plants.
*/
void photosynthesis_gs(Patch& patch) {
	const Vegetation &vegetation = patch.vegetation;
	const Climate &climate = patch.get_climate();

	PhotosynthesisInputs inputs;
	PhotosynthesisOutputs outputs;

	// Loop through cohorts (=individuals)
	for (int i = 0; i < vegetation.nobj; i++) {

		// Retrieve information about cohort
		Individual& indiv = vegetation[i];
		LogScope log_scope(indiv);
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if ((indiv.pft.ifcalcvcmax && is_zero(indiv.get_nmass_leaf())) ||
		    is_zero(indiv.get_wscal())) {
			continue;
		}
		// TODO: find another home for this check. It doesn't belong here.
		double lai_layer_sum = 0;
		for (size_t j = 0; j < indiv.canopy.size(); j++) {
			lai_layer_sum += indiv.canopy[j].lai;
		}
		if (!equal(lai_layer_sum, indiv.lai_today())) {
			log_warning("%s: layerwise LAI sum (%g) does not equal total LAI (%g) (phen=%g)",
						(const char *)pft.name, lai_layer_sum,
						indiv.lai_today(), indiv.phen);
		}

		indiv.photosynthesis.agd_g = 0.0;
		indiv.photosynthesis.rd_g = 0.0;
		indiv.photosynthesis.je = 0.0;

		double an_total = 0;
		indiv.anc_frac = 0;
		indiv.anj_frac = 0;
		indiv.anp_frac = 0;

		// Loop over time steps per day
        for (Day day; day.period != date.subdaily; day.next()) {
			// Initialise timestep-level variables.
			double Rd_canopy = 0.0;
			double An_canopy = 0.0;
			double Ag_canopy = 0.0;
			double Anc_canopy = 0.0;
			double Anj_canopy = 0.0;
			double Anp_canopy = 0.0;
			double gsw_canopy = 0.0;
			double ci_canopy = 0.0;
			double vcmax_canopy = 0.0;
			double jmax_canopy = 0.0;

			// Initialise cumulative LAI.
			double lai_cum_bottom = 0.0;

			// Loop over canopy layers (from top layer to bottom layer)
			size_t nlayer = indiv.canopy.size();

            // Initialise the layerwise inputs and outputs for this timestep.
            // This can't happen earlier because we need to know the number of
            // layers, which isn't known until we call fpar().
			indiv.ps_inputs_layerwise[day.period].assign(nlayer, PhotosynthesisInputs{});
			indiv.ps_outputs_layerwise[day.period].assign(nlayer, PhotosynthesisOutputs{});
			indiv.energy_balance_layerwise[day.period].assign(nlayer, EnergyBalance(0, 0, 0, false, 0, 0));
			indiv.ps_diagnostics_layerwise[day.period].assign(nlayer, PhotosynthesisDiagnostics{});

		    for (size_t j = 0; j < nlayer; j++) {

				// No LAI means we can skip this layer as an optimisation. This
				// can happen because the number of layers is the total number
				// of canopy layers. So if this tree is shorter than the tallest
				// tree in the patch, it may not have any LAI in this layer.
				if (is_zero(indiv.canopy[j].lai)) {
					continue;
				}

				// Determine, for this indiv, LAI and cumulative LAI of the layer.

				// Cumulative LAI at the bottom of current layer.
				lai_cum_bottom += indiv.canopy[j].lai;

				// Cumulative LAI at the middle of current layer.
				double lai_cum = lai_cum_bottom - 0.5 * indiv.canopy[j].lai;

				CohortState cohort(indiv, j);

				PhotosynthesisDiagnostics diagnostics = get_photosynthesis_inputs(
					cohort, pft, climate, patch, day, lai_cum, inputs);
				indiv.g1 = diagnostics.g1;

				ClimateInputs clim(climate, day);

				EnergyBalance balance = solve_energy_balance_dave(inputs,
					clim, outputs, diagnostics, pft.pathway,
					energy_balance_maxiter, energy_balance_eps);

				// If the user requested 1 iteration max, they're obviously not
				// interested in whether the energy balance converged, so we can
				// refrain from spamming them with log messages. In this case,
				// it doesn't really matter whether it converged anyway - tleaf
				// is effectively only an output.
				if (energy_balance_maxiter > 1) {
					if (balance.get_converged()) {
						log_debug("Energy balance solver converged for %s in layer %zu at day period %d in %d iterations.",
									(const char *)pft.name, j, day.period, balance.get_niter());
					} else {
						log_diagnostic("Energy balance solver did not converge for %s in layer %zu at day period %d after %d iterations (tleaf_diff = %g).",
									(const char *)pft.name, j, day.period,
									balance.get_niter(), balance.get_delta_tleaf());
					}
				}

				indiv.ps_inputs_layerwise[day.period][j] = inputs;
				indiv.ps_outputs_layerwise[day.period][j] = outputs;
				indiv.energy_balance_layerwise[day.period][j] = balance;
				indiv.ps_diagnostics_layerwise[day.period][j] = diagnostics;

				// TODO drew: multiply by LAI only once, after canopy loop.

				// sum variables over canopy layers
				double lai_layer = inputs.lai;

				Rd_canopy  += outputs.rd * lai_layer;
				An_canopy  += outputs.an * lai_layer;
				Ag_canopy  += (outputs.an + outputs.rd) * lai_layer;
				Anc_canopy += outputs.anc * lai_layer;
				Anj_canopy += outputs.anj * lai_layer;
				Anp_canopy += outputs.anp * lai_layer;
				double gsw = get_gsw(outputs.gsc, diagnostics.gcutic_adj);
				gsw_canopy += gsw * lai_layer;
				ci_canopy  += outputs.ci * lai_layer;
				vcmax_canopy += outputs.vcmax * lai_layer;
				jmax_canopy  += outputs.jmax * lai_layer;

				assert(positive(outputs.gsc));

				double anc = outputs.anc;
				double anj = outputs.anj;
				double anp = outputs.anp;
				if (anc < anj && (pft.pathway == C4 || anc < anp) ) {
					// Anc-limited
					anc = max(0.0, anc);
					double anc_layer = anc * lai_layer;
					indiv.anc_frac += anc_layer;
					an_total += anc_layer;
				} else if (anj < anc && (pft.pathway == C4 || anj < anp) ) {
					// anj-limited
					anj = max(0.0, anj);
					double anj_layer = anj * lai_layer;
					indiv.anj_frac += anj_layer;
					an_total += anj_layer;
				} else if (pft.pathway == C3 && anp < anc && anp < anj) {
					anp = max(0.0, anp);
					double anp_layer = anp * lai_layer;
					indiv.anp_frac += anp_layer;
					an_total += anp_layer;
				}
			} // end canopy layer loop

			// lai_layerwise should sum to lai_today(), which is lai * phen. For
			// deciduous trees not using the DAVE allocation pathway, this may
			// not be the same as indiv.lai.

			// Calculate mean for some variables.
			if (lai_layer_sum == 0) {
				log_warning("%s: Total layerwise LAI is zero (%g) (indiv.lai = %g)",
					        (const char *)pft.name, lai_layer_sum,
							indiv.lai_today());
			} else {
				ci_canopy /= lai_layer_sum;
				vcmax_canopy /= lai_layer_sum;
				jmax_canopy /= lai_layer_sum;
			}

			assert(indiv.lai_today() < 20);

			// Save key variables to PhotosynthesisSubdaily Structure.
			indiv.ps_subdaily.Rd[day.period]    = Rd_canopy;
			indiv.ps_subdaily.An[day.period]    = An_canopy;
			indiv.ps_subdaily.Anc[day.period]   = Anc_canopy;
			indiv.ps_subdaily.Anj[day.period]   = Anj_canopy;
			indiv.ps_subdaily.gsw[day.period]   = gsw_canopy;
			indiv.ps_subdaily.ci[day.period]    = ci_canopy;
			indiv.ps_subdaily.vcmax[day.period] = vcmax_canopy;
			indiv.ps_subdaily.jmax[day.period]  = jmax_canopy;

			// sum up over the day
			indiv.photosynthesis.agd_g += Ag_canopy;
			indiv.photosynthesis.rd_g  += Rd_canopy;
			indiv.photosynthesis.je    += indiv.ps_subdaily.Anj[day.period];

			// Run hydraulics module and check for stomatal closure in next
			// timestep.
			if (ifdavehydraulics) {
				if (indiv.pft.is_grass()) {
					grass_hydraulics(day, indiv);
				} else {
					tree_hydraulics(day, indiv);
				}
			}
		} // end subdaily loop

		if (an_total > 0) {
			indiv.anc_frac /= an_total;
			indiv.anj_frac /= an_total;
			indiv.anp_frac /= an_total;
		}

		assert(indiv.anc_frac <= 1.0);
		assert(indiv.anj_frac <= 1.0);
		assert(indiv.anp_frac <= 1.0);
		assert(indiv.anc_frac >= 0.0);
		assert(indiv.anj_frac >= 0.0);
		assert(indiv.anp_frac >= 0.0);

		// TODO Drew: make assert statement work
		//if (an_total > 0) {	
			//assert(equal(indiv.anc_frac + indiv.anj_frac + indiv.anp_frac, 1.0));
		//}

		// Calculate daily means for the PhotosynthesisResult Structure that goes to output
		indiv.photosynthesis.agd_g /= date.subdaily;
		indiv.photosynthesis.rd_g  /= date.subdaily;
		indiv.photosynthesis.je    /= date.subdaily;

		// Convert units from molar to mass units and from seconds to day
		indiv.photosynthesis.agd_g *= MOL2G_C * SECS_PER_DAY;  // gC m-2 d-1
		indiv.photosynthesis.rd_g *= MOL2G_C * SECS_PER_DAY;   // gC m-2 d-1
		// Give values all variables even they aren't used (yet)
		// indiv.photosynthesis.vm          = indiv.ps_subdaily.vcmax25 / UMOL_PER_MOL * MOL2G_C * SECS_PER_DAY;  // umol m-2 s-1; not used
		indiv.photosynthesis.vm          = max(indiv.ps_subdaily.vcmax) * MOL2G_C * SECS_PER_DAY;
		indiv.photosynthesis.adtmm       = 0.0;    // not used
		indiv.photosynthesis.vmaxnlim    = 0.0;    // not used 
	} // end cohort loop
}

double calc_water_uptake_dave(Individual &indiv) {
    Patch &patch = indiv.vegetation.patch;
	const Soil &soil = patch.soil;
    Pft &pft = indiv.pft;
    const Standpft& spft = patch.stand.pft[indiv.pft.id];

    if (spft.irrigated && pft.id == patch.stand.pftid) {
        Day day;
        return irrigated_water_uptake(patch, pft, day, indiv.rootdist, indiv.fwuptake);
    }

    // Available water in each layer, as a fraction of water holding
    // capacity.
    double wcont_local[NSOILLAYERS_MAX];
	double awc[NSOILLAYERS_MAX];
    for (int i = 0; i<soil.nsoillayers; i++) {
        wcont_local[i] = soil.soillayers[i].wcont();
		awc[i] = soil.soillayers[i].awc;
    }

    return water_uptake(wcont_local, awc, indiv.rootdist, pft.emax,
		                patch.fpc_rescale, indiv.fwuptake,
                        pft.lifeform == TREE, pft.drought_tolerance, indiv.Rs_r,
				        soil.nsoillayers);
}

double calc_transp_old(const Individual &indiv) {
	const Climate &climate = indiv.vegetation.patch.get_climate();
	double transp = 0;

	// Loop over time steps per day
	for (Day day; day.period != date.subdaily; day.next()) {	

		// T = Gc / (pressure/VPD)
		// Gc in mol m-2 s-1
		// VPD and pressure in Pa
		double patmos = climate.pressure[day.period];
		double vpd = climate.vpds[day.period];
		double gsw = indiv.ps_subdaily.gsw[day.period];
		if (vpd > 0 && patmos > 0) {
			transp += gsw / (patmos / (vpd * 1000.0));
		}

	} // end sub-daily loop
	
	// daily transpiration rate in mol m-2 s-1
	transp /= date.subdaily;

	// daily transpiration rate in mm d-1
	// Limiting transpiration to the actual supply is handled later.
	return transp * MOL2G_H2O * KG_PER_G * SECS_PER_DAY;
}

/*
Aggregate transpiration across canopy layers and timesteps to get a total daily
transpiration in mm day-1 for this cohort.
*/
double transpiration_dave(const Individual &indiv) {
	// Loop over time steps per day
	double transp = 0;
	for (Day day; day.period != date.subdaily; day.next()) {	
		// Iterate through canopy layers for this timestep.
		// const& to avoid copy - important!
		const std::vector<EnergyBalance> &ebal = indiv.energy_balance_layerwise[day.period];
		for (size_t j = 0; j < ebal.size(); j++) {
			// Transpiration in mol H2O m-2 s-1.
			// TODO: should transpiration in each layer/timestep be clamped to
			// zero if it's negative? Or should we do this once per day?
			transp += ebal[j].get_transp() * indiv.canopy[j].lai;
		}
	}

	// Daily transpiration rate in mol m-2 s-1.
	transp /= date.subdaily;

	// Daily transpiration rate in mm d-1.
	// Limiting transpiration to the actual supply is handled later.
	// assuming 1kg m-2 = 1l m-2 == 1mm
	transp = transp * MOL2G_H2O * KG_PER_G * SECS_PER_DAY;

	// Don't allow negative daily transpiration.
	return max(0.0, transp);
}

/*
Function to redistribute extra AET if any one layer's water content is exceeded.
*/
void prevent_wuptake_oversubscription(Patch &patch) {
	Vegetation &vegetation = patch.vegetation;
	Soil &soil = patch.soil;

	// Allocate uptake explicitly so that demand in a dry layer can be moved to
	// other layers without reducing transpiration from the whole soil profile.

	// Amount of water taken up from each layer by each individual (mm).
	std::vector<std::vector<double> > uptake(vegetation.nobj, std::vector<double>(soil.nsoillayers, 0.0));

	std::vector<std::vector<double> > uptake_weights(vegetation.nobj, std::vector<double>(soil.nsoillayers, 0.0));
	std::vector<double> original_aet(vegetation.nobj, 0.0);
	std::vector<double> supply(soil.nsoillayers, 0.0);
	std::vector<bool> depleted(soil.nsoillayers, false);
	double total_demand = 0.0;
	double total_supply = 0.0;

	for (int i = 0; i < soil.nsoillayers; i++) {
		Soillayer &layer = soil.soillayers[i];

		supply[i] = max(0.0, layer.wmass - layer.wp);
		depleted[i] = negligible(supply[i]);
		total_supply += supply[i];
	}

	for (int j = 0; j < vegetation.nobj; j++) {
		const Individual &indiv = vegetation[j];
		original_aet[j] = indiv.aet;
		double indiv_uptake = 0.0;
		for (int i = 0; i < soil.nsoillayers; i++) {
			uptake_weights[j][i] = indiv.fwuptake[i];
			uptake[j][i] = indiv.aet * indiv.fwuptake[i];
			indiv_uptake += uptake[j][i];
		}

		// A hydraulically isolated individual deliberately has no soil uptake.
		assert(negligible(indiv_uptake) || equal(indiv_uptake, original_aet[j]));
		total_demand += indiv_uptake;
	}

	if (total_demand > total_supply) {
		double supply_fraction = total_supply / total_demand;

		for (size_t i = 0; i < soil.nsoillayers; i++) {
			const Soillayer &layer = soil.soillayers[i];
			log_diagnostic("Layer %zu %s: supply=%fmm, psi=%fMPa, wp=%fmm, wcont=%fmm, awc=%fmm",
			          i, depleted[i] ? " (depleted)" : "", supply[i],
					  layer.psi_mpa(), layer.wp, layer.wmass, layer.awc);
		}
		log_error("water_uptake_dave(): transpiration demand (%fmm) exceeds the "
			 "water available (%fmm) in the soil profile",
			 total_demand, total_supply);

		for (int j = 0; j < vegetation.nobj; j++) {
			for (int i = 0; i < soil.nsoillayers; i++)
				uptake[j][i] *= supply_fraction;
		}
	}

	bool oversubscribed;
	int redistribution_count = 0;
	do {
		oversubscribed = false;

		// The amount of water which couldn't be taken up by each individual due
		// to oversubscription of a layer.
		std::vector<double> rejected(vegetation.nobj, 0.0);

		// Cap oversubscribed layers. Individuals share the available water in
		// proportion to their uptake requested from that layer.
		for (int i = 0; i < soil.nsoillayers; i++) {
			double layer_uptake = 0.0;
			for (int j = 0; j < vegetation.nobj; j++)
				layer_uptake += uptake[j][i];

			if (layer_uptake > supply[i]) {
				oversubscribed = true;
				depleted[i] = true;
				double available_fraction = supply[i] / layer_uptake;
				for (int j = 0; j < vegetation.nobj; j++) {
					double accepted = uptake[j][i] * available_fraction;
					rejected[j] += uptake[j][i] - accepted;
					uptake[j][i] = accepted;
				}
			}
		}

		// Redistribute each individual's rejected uptake according to its
		// original preferences among layers which have not been depleted.
		if (oversubscribed) {
			for (int j = 0; j < vegetation.nobj; j++) {
				if (negligible(rejected[j]))
					continue;

				double weight_sum = 0.0;
				for (int i = 0; i < soil.nsoillayers; i++) {
					if (!depleted[i])
						weight_sum += uptake_weights[j][i];
				}

				if (negligible(weight_sum) && total_demand <= total_supply) {
					fail("water_uptake_dave(): unable to redistribute water uptake "
						"to a non-depleted layer");
				}

				for (int i = 0; i < soil.nsoillayers; i++) {
					if (!depleted[i])
						uptake[j][i] += rejected[j] *
							uptake_weights[j][i] / weight_sum;
				}
			}
		}

		redistribution_count++;
		if (redistribution_count > soil.nsoillayers && total_demand <= total_supply) {
			fail("water_uptake_dave(): water uptake redistribution did not converge");
		}
	} while (oversubscribed);

	for (int i = 0; i < soil.nsoillayers; i++) {
		double layer_uptake = 0.0;
		for (int j = 0; j < vegetation.nobj; j++)
			layer_uptake += uptake[j][i];
		assert(layer_uptake < supply[i] || equal(layer_uptake, supply[i]));
	}

	// Convert the allocation back to the representation used by the soil code.
	// Redistribution must conserve each active individual's original AET.
	for (int j = 0; j < vegetation.nobj; j++) {
		Individual &indiv = vegetation[j];
		double indiv_uptake = 0.0;
		for (int i = 0; i < soil.nsoillayers; i++)
			indiv_uptake += uptake[j][i];

		// Ensure that total uptake remains unchanged under normal conditions.
		if (total_demand <= total_supply && !negligible(indiv_uptake)) {
			assert(equal(indiv_uptake, original_aet[j]));
		}

		indiv.aet = indiv_uptake;
		for (int i = 0; i < soil.nsoillayers; i++) {
			indiv.fwuptake[i] = negligible(indiv.aet) ?
				0.0 : uptake[j][i] / indiv.aet;
		}
	}
}

/*
Calculate daily water uptake by layer (but don't take it out of the soil).
*/
void water_uptake_dave(Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;

	// Loop through cohorts (=individuals)
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		if (penman_monteith) {
			indiv.aet = transpiration_dave(indiv);
		} else {
			indiv.aet = calc_transp_old(indiv);
		}

		indiv.aaet += indiv.aet;

		// Note that in standard LPJ-GUESS transpiration is described as
		// evapotranspiration (AET), which is incorrect as soil evaporation is
		// calculated separately.

		// Limit water uptake.
		double wr = calc_water_uptake_dave(indiv);

        if (ifdavehydraulics && !indiv.stomata_open) {
            for (int i = 0; i < NSOILLAYERS_MAX; i++) {
                // If the stomata are closed, we assume roots are hydraulically
                // isolated from the soil.
                indiv.fwuptake[i] = 0;
            }
        }
	} // end cohort loop

	// Ensure water uptake doesn't exceed availability in any layer.
	prevent_wuptake_oversubscription(patch);
}

struct CtonRatios {
	double sap;
	double leaf;
	double root;
	double crown;
};

/*
Get the effective C:N ratios for sap, root, and crown to use in respiration
calculations, based on the individual's C:N ratios and whether N limitation is
enabled.
*/
CtonRatios cton_effective(const Individual &indiv) {
	CtonRatios cton{};

	// Static root and sap wood C:N ratio if no N limitation to not let N affect
	// respiration for C only version of model Also not allow lower C:N ratios
	// than minimum which would affect respiration severely (happens for
	// deciduous trees in the beginning of the season)

	if (ifnlim) {
		cton.sap   = indiv.cton_sap();
		cton.leaf  = indiv.cton_leaf();
		cton.root  = indiv.cton_root();
		cton.crown = indiv.cton_crown();
	} else {
		cton.sap   = indiv.pft.cton_sap_avr;
		cton.root  = indiv.pft.cton_root_avr;
		cton.leaf  = indiv.pft.cton_leaf_avr;
		// Crowns use root stoichiometric limits.
		cton.crown = indiv.pft.cton_root_avr;
	}

	if (!indiv.pft.dave_allocation()) {
		// Set C:N ratios if not using DAVE allocation (since trunk
		// allocation + DAVE photosynthesis can give unrealistic values).
		cton.leaf  = min(cton.leaf, indiv.pft.cton_leaf_max);
		cton.sap   = min(cton.sap, indiv.pft.cton_sap_max);
		cton.root  = min(cton.root, indiv.pft.cton_root_max);
		// Non-DAVE allocated PFTs don't have crowns.
	}

	return cton;
}

/*
This is a modification of the original npp routine of standard LPJ-GUESS.

Modifications include:

- called once a day (no sub-daily option, day structure not an input)
- no more calculation of photosynthesis
- bvoc function not implemented (could be added later)
*/
void npp_dave(Patch& patch) {

	// Response of respiration rate to temperature
	double gtemp = patch.get_climate().get_gtemp_mean();

	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		Individual& indiv = patch.vegetation[i];
		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[pft.id];

		// Don't do calculations for crops outside their growing season.
		if (!indiv.growingseason()) {
			indiv.dnpp=0.0;
			indiv.dnpp_allocated = true;
			indiv.lambda = 0;
			continue;
		}

		// get daily gross assimilation (kgC/m2)
		double assim = indiv.photosynthesis.gross_assimilation();

		// Calculate autotrophic respiration

		double cmass_root = indiv.cmass_root_today();

		// 0 for trees
		double cmass_culm = indiv.get_cmass_leaf() * pft.culmfrac;

		// 0 for trees
		double cmass_crown = indiv.get_cmass_crown();

		// 0 for PFTs using annual (trunk) allocation
		double cmass_repr = indiv.get_cmass_repr();

		// Get effective C:N ratios for sap, root, and crown.
		CtonRatios cton = cton_effective(indiv);

		// Calculate respiration scalar due to dormancy
		double fphen_resp = 1;
		if (pft.is_grass()) {
			fphen_resp = indiv.fphen_activity + (1.0 - indiv.fphen_activity) *
		                    max(0.0,(indiv.get_wscal_moving_average() - indiv.pft.wscal_min_dormant) / (1.0 - indiv.pft.wscal_min_dormant));
		}


		// Set leaf respiration.
		indiv.resp_leaf = indiv.photosynthesis.rd_g * KG_PER_G;

		// 30 day moving average of air temperature.
		double tgrowth = patch.get_climate().get_tgrowth();

        // Get scalar encapsulating the effect of temperature acclimation. This
        // will be 1 if acclimation is disabled.
        double scalar = resp_temp_acclimation_scalar(tgrowth, pft.racclim);
        double wscalar = resp_water_scalar(indiv);

		// Respiration coefficient is reduced based on temperature acclimation.
        double respcoeff = pft.respcoeff * scalar * wscalar;

		double gtemp = patch.soil.gtemp_rootdist(ppft.rootdist);
		double gtemp_air = mean(patch.get_climate().gtemps);

		double resp = 0;
		respiration(gtemp_air, gtemp, indiv.pft.lifeform, respcoeff, cton.sap,
			cton.leaf, cton.root, cton.crown, indiv.cmass_sap, cmass_culm,
			cmass_root, cmass_crown, cmass_repr, indiv.cgrow, fphen_resp, assim,
			indiv.resp_leaf, resp, indiv.resp_maintenance, indiv.resp_growth,
			indiv.resp_repr, indiv.resp_culm, indiv.resp_crown, indiv.resp_root,
			indiv.resp_sap);

		// Update accumulated annual NPP and daily vegetation-atmosphere flux
		// Note that resp does not include leaf day respiration.
		indiv.dgpp = assim;
		indiv.dnpp = assim - resp;
		indiv.anpp += indiv.dnpp;
		indiv.dnpp_allocated = !pft.dave_allocation();
		indiv.agpp	+= assim;
		indiv.aresp	+= resp;

		// Save historical GPP for daily growth of grass
		indiv.gpp_15.add(assim);
		indiv.mgpp[date.month] += max(0.0, assim);

		indiv.report_flux(Fluxes::NPP, indiv.dnpp);
		indiv.report_flux(Fluxes::GPP, assim);
		indiv.report_flux(Fluxes::RA, resp);
	}
}

/*
Calculate potential forest-floor photosynthesis for all PFTs not using daily
establishment.
*/
void forest_floor_dave(Patch &patch) {
	// These calculations are only required for pfts using the annual (trunk)
	// establishment.

	Vegetation &vegetation = patch.vegetation;
	const Climate &climate = patch.get_climate();

	for (Day day; day.period < date.subdaily; day.next()) {
		// Assume LAI of 1 m2 m-2.
		const double LAI_HYPOTHETICAL = 1.0;

		for (size_t i = 0; i < pftlist.nobj; i++) {
			const Pft &pft = pftlist[i];
			Patchpft &ppft = patch.pft[i];

			// Skip this calculation for PFTs using daily establishment.
			if (pft.daily_establishment()) {
				continue;
			}

			// Since the seedlings don't have access to deep layers, only look
			// at the top few layers' SWP.
			double lwp = patch.soil.var_depth(0, FOREST_FLOOR_WSCAL_DEPTH, &Soillayer::psi_mpa);

			// No soil water on first day of simulation - this results in wscal
			// of 0, which results in a division by zero. We can just skip this
			// calculation on the first day of the simulation.
			if (is_zero(ppft.wscal)) {
				continue;
			}

			CohortState cohort(LAI_HYPOTHETICAL, ppft.wscal, lwp, lwp, pft.stom_clos, patch.fpar_ff,
							   true, 1.0, 1.0, pft.vcmax_25, pft.jmax_25,
							   pft.vpmax_25);
			PhotosynthesisInputs inputs;
			get_photosynthesis_inputs(cohort, pft, climate, patch, day,
									  LAI_HYPOTHETICAL / 2.0, inputs);

			// Assume a leaf N content of 0.01 kgN/m2.
			PhotosynthesisOutputs outputs;
			if (pft.pathway == C3) {
				photosynthesis_dave_c3(inputs, outputs);
			} else {
				photosynthesis_dave_c4(inputs, outputs);
			}

			// Convert from mol/m2/s to kgC/m2/day.
			double secs_per_timestep = SECS_PER_DAY / date.subdaily;
			ppft.anetps_ff += outputs.an * MOL2G_C / 1000 * secs_per_timestep;
		}
	}

	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];
		Patchpft &ppft = patch.pft[i];
		Standpft &spft = patch.stand.pft[pft.id];
		if (date.islastday && date.islastmonth) {
			spft.anetps_ff_max = max(spft.anetps_ff_max, ppft.anetps_ff);
			ppft.anetps_ff = max(0.0, ppft.anetps_ff);
		}
	}
}

/*
Calculate the LAI in a given layer for an individual.

@param indiv Individual for which to calculate LAI.
@param lowbound Lower height of layer (m).
@param highbound Upper height of layer (m).
@param vstep Width of layer (m).
*/
double calc_lai_layer(const Individual &indiv, double lowbound,
					  double highbound, double vstep) {
	if (!(indiv.pft.lifeform == TREE &&
		  indiv.height > lowbound &&
		  indiv.boleht < highbound &&
		  !negligible(indiv.height - indiv.boleht))) {
		return 0;
	}

	// Calculate vertical fraction of current layer occupied by
	// crown cylinders of this cohort

	double frac = 1.0;
	if (indiv.height < highbound)
		frac -= (highbound - indiv.height) / vstep;
	if (indiv.boleht > lowbound)
		frac -= (indiv.boleht - lowbound) / vstep;

	// Calculate summed LAI of this cohort in this layer

	double atoh = indiv.lai_today() / (indiv.height - indiv.boleht);
	return atoh * frac * vstep;
}

/*
Calculate total tree LAI in a layer.

@param lowbound Lower height of layer (m).
@param upbound Upper height of layer (m).
@param vstep Width of layer (m).
*/
double calc_lai_layer(const Vegetation &vegetation, double lowbound,
					  double upbound, double vstep) {
	double lai = 0;
	for (size_t i = 0; i < vegetation.nobj; i++) {
		lai += calc_lai_layer(vegetation[i], lowbound, upbound, vstep);
	}
	return lai;
}

/*
Calculate the fraction of incoming PAR incident to a layer (0-1).

@param lai Cumulative LAI from top of canopy to middle of layer (m2 m-2).
@param k Extinction coefficient for PAR.
*/
double calc_fipar(double lai, double k = K_BEER_GRASS) {
	return k / (1.0 - M_TRANS) * lambertbeer(lai, k);
}

void partition_tree_canopy_dave(Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;

	// Canopy height in meters.
	double height_veg = max(vegetation, [](const Individual &indiv) { return indiv.height; });

	// Width of vertical layers for canopy-area integration (m).
	double vstep = nlayer_canopy ? height_veg / nlayer_canopy : layer_thickness_canopy;

	// Calculate number of layers (minus 1) from ground surface to top of canopy
	size_t toplayer = (size_t)(height_veg / vstep);

	if (toplayer > 0 && height_veg <= toplayer * vstep) {
		// Top layer results in an empty layer. This can happen if the
		// vegetation height falls exactly on a layer boundary. In this case
		// we simply reduce the layer count by one.
		toplayer--;
	}

	size_t nlayer = toplayer + 1;

	double lai_cum = 0;
	for (size_t i = 0; i < nlayer; i++) {
		// i = 0 is top layer.
		size_t layer = toplayer - i;
		double lowbound = (double)layer * vstep;
		double upbound = (double)(layer + 1) * vstep;

		// Total LAI of all cohorts in this layer.
		double lai_layer = calc_lai_layer(vegetation, lowbound, upbound, vstep);
		lai_cum += lai_layer;

		// see e.g. Sands 1995, Australian Journal of Plant Physiology
		double fipar_layer = calc_fipar(lai_cum - lai_layer / 2, K_BEER_TREE);

		// Store layer-wise fpar and LAI in tree cohorts.
		for (size_t j = 0; j < vegetation.nobj; j++) {
			Individual &indiv = vegetation[j];
			if (indiv.pft.lifeform != TREE) {
				continue;
			}
			double lai = calc_lai_layer(indiv, lowbound, upbound, vstep);
			indiv.canopy.set(i, lai, fipar_layer);
		}
	}
}

// LAI per layer for each grass individual.
const double LAI_PER_LAYER = 0.5;

void partition_grass_canopy_even_dave(const Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;

	// Get total grass LAI (m2 m-2).
	double lai_grass = sum(vegetation, [](const Individual &indiv) {
			return indiv.pft.is_grass() ? indiv.lai_today() : 0; });

	// Total number of canopy layers for grasses.
	size_t nlayer = ceil(lai_grass / LAI_PER_LAYER);

	// Iterate through canopy layers top to bottom. Top layer is i=0.
	double lai_cum = 0;
	for (size_t i = 0; i < nlayer; i++) {
		// In bottom layer we allocate all remaining LAI. All other layers get
		// LAI_PER_LAYER.
		double lai_layer = LAI_PER_LAYER;
		if (i == nlayer - 1) {
			lai_layer = lai_grass - lai_cum;
		}
		lai_cum += lai_layer;

		// Calculate fipar in middle of layer.
		double fipar_layer = calc_fipar(lai_cum - lai_layer / 2, K_BEER_GRASS);

		// All grasses get the same fraction of LAI in this layer.
		double lai_frac = lai_layer / lai_grass;

		for (size_t j = 0; j < vegetation.nobj; j++) {
			Individual &indiv = vegetation[j];
			if (!indiv.pft.is_grass()) {
				continue;
			}

			// Evenly distribute LAI across layers.
			double lai = indiv.lai_today() * lai_frac;

			// Add LAI in to lai_layerwise vector.
			indiv.canopy.set(i, lai, fipar_layer);
		}
	}
}

std::vector<int> get_sorted_grass_indices(const Vegetation &vegetation) {
	std::vector<int> indices;
	indices.reserve(vegetation.nobj);

	// Collect in vegetation order.
	for (size_t i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		bool is_grass = indiv.pft.is_grass() || indiv.pft.lifeform == MOSS;
		if (is_grass && indiv.lai_today() > 0) {
			indices.push_back(i);
		}
	}

	// Sort ascending by LAI.
	std::stable_sort(indices.begin(), indices.end(),
		[&vegetation](int a, int b) {
			return vegetation[a].lai_today() < vegetation[b].lai_today();
		});
	return indices;
}

void partition_grass_canopy_uneven_dave(const Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;

	// Get indices of grasses in ascending order of LAI.
	std::vector<int> grass_indices = get_sorted_grass_indices(vegetation);

	// Total grass LAI (m2 m-2).
	double lai_grass = sum(vegetation, [](const Individual &indiv) {
		return indiv.pft.is_grass() ? indiv.lai_today() : 0; });

	size_t nlayer = static_cast<size_t>(ceil(lai_grass / LAI_PER_LAYER));

	// The total LAI allocated so far for each grass.
	std::vector<double> allocated(grass_indices.size(), 0);

	// Iterate over grass canopy layers, starting at bottom of canopy.
	// Note: i=0 is top layer, which is the indexing used by indiv.canopy.
	for (int i = nlayer - 1; i >= 0; i--) {
		// Total lai in this layer across all grass PFTs.
		double lai_layer = LAI_PER_LAYER;

		// Bottom LAI gets remaining LAI which is <= LAI_PER_LAYER.
		if (i == nlayer - 1) {
			lai_layer = lai_grass - (nlayer - 1) * LAI_PER_LAYER;
		}

		// i is the number of layers above this one. Every layer except the
		// bottom one is guaranteed to have exactly LAI_PER_LAYER.
		double lai_above = i * LAI_PER_LAYER;
		double fipar = calc_fipar(lai_above + 0.5 * lai_layer, K_BEER_GRASS);

		// Each individual's entitlement in this layer. All grasses have this
		// much LAI in this layer. If some shortfall (ie remaining LAI is less
		// than their entitlement), entitlement for remaining grasses is
		// increased. Because grasses are ordered by LAI ascending, we ensure
		// that lai_layer_total is allocated to the maximum extent possible.
		double entitlement = lai_layer / grass_indices.size();

		// The LAI not yet allocated in this layer.
		double remaining = lai_layer;

		// Loop through grasses, smallest to largest.
		for (size_t j = 0; j < grass_indices.size(); j++) {
			Individual &indiv = vegetation[grass_indices[j]];

			// Calculate LAI contributed by this grass in this layer.
			double lai = min(indiv.lai_today() - allocated[j], entitlement);
			indiv.canopy.set(i, lai, fipar);
			allocated[j] += lai;

			// Subtract this grass's LAI contribution from the remaining LAI
			// required in this layer.
			remaining -= lai;

			if (lai < entitlement) {
				// This grass is too small to use its full entitlement.
				// Increase entitlement for remaining grasses.
				int nindiv_remaining = grass_indices.size() - j - 1;
				if (nindiv_remaining > 0) {
					entitlement = remaining / nindiv_remaining;
				}
			}
		} // grasses loop
	} // layer loop
}

void partition_grass_canopy_dave(Patch &patch) {
	// Short circuit if no grasses are established.
	double ngrass = sum(patch.vegetation, [](const Individual &indiv) { return indiv.pft.is_grass() ? 1.0 : 0.0; });
	if (ngrass == 0) {
		return;
	}

	if (ifevenlaidist) {
		partition_grass_canopy_even_dave(patch);
	} else {
		partition_grass_canopy_uneven_dave(patch);
	}
}

/*
Partition canopy into discrete layers, and calculate LAI and fipar in each layer
for each cohort.
*/
void partition_canopy_dave(Patch &patch) {
	// Short-circuit empty patches.
	if (patch.vegetation.nobj == 0) {
		return;
	}

	partition_tree_canopy_dave(patch);
	partition_grass_canopy_dave(patch);

	// Verify canopy vs canopy2.
	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		const Individual &indiv = patch.vegetation[i];
		LogScope log_scope(indiv);
		indiv.validate_canopy();
	}
}

/*
Damage caused by sitting at temprature tleaf for mins minutes.

@param tleaf   Leaf temperature (℃).
@param t50_par Reference temperature (°C) at which 50% damage occurs
               after 1 minute of exposure.
@param z       Temperature sensitivity coefficient (°C). Controls how
               rapidly damage changes with temperature; smaller values
               produce a steeper damage response curve.
@param mins    Length of exposure (minutes).

@return Fraction of LAI lost (0-1).
*/
double calc_leaf_damage(double tleaf, double t50_par, double z, double mins) {
	double dam = 0.5 * mins * pow(10.0, -1.0 / z * (t50_par - tleaf));
	dam = std::min(1.0, dam);
	return dam;
}

void leaf_damage_dave(Patch &patch) {
	if (!ifleafdamage) {
		return;
	}

	Soil &soil = patch.soil;

	// Timestep length in minutes.
	double timestep = SECS_PER_DAY / date.subdaily / SECONDS_PER_MINUTE;

	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		Individual &indiv = patch.vegetation[i];
		LogScope log_scope(indiv);
		const Pft &pft = indiv.pft;

		if (is_zero(indiv.lai_today())) {
			continue;
		}
		double lai_loss = 0;
		for (size_t j = 0; j < indiv.energy_balance_layerwise.size(); j++) {
			for (size_t k = 0; k < indiv.energy_balance_layerwise[j].size(); k++) {
				double tleaf = indiv.energy_balance_layerwise[j][k].get_tleaf();
				double dam = calc_leaf_damage(tleaf, pft.t50, pft.z, timestep);
				lai_loss += dam * indiv.canopy[k].lai;
			}
		}

		if (lai_loss == 0) {
			continue;
		}

		double frac_loss = lai_loss / indiv.lai_today();
		assert(frac_loss >= 0);

		// Sum of indiv.canopy[k].lai across layers can't exceed
		// indiv.lai_today(), but over the course of the day it could in theory,
		// so we need to guard against that.
		if (frac_loss >= 1.0) {
			log_warning("%s: LAI loss (%g) exceeds total LAI (%g), setting LAI loss to total LAI. This isn't mathematically impossible, and probably indicates sustained high temperatures over the course of the day.",
						(const char *)pft.name, lai_loss, indiv.lai_today());
			frac_loss = 1.0;
		}

		double frac_remain = 1.0 - frac_loss;

		double closs = frac_loss * indiv.cmass_leaf_today();
		double nloss = frac_loss * indiv.get_nmass_leaf();

		// Move biomass from leaf to litter.
		if (pft.is_grass() && pft.dave_allocation()) {
			// get_cmass_leaf() and get_nmass_leaf() only include green leaf
			// mass, so that's what we reduce. Should we also reduce brown leaf
			// mass? And should green leaf mass be moved to litter or to brown
			// leaf mass?
			indiv.cmass_leafs[green] *= frac_remain;
			indiv.nmass_leafs[green] *= frac_remain;
		} else {
			indiv.cmass_leaf *= frac_remain;
			indiv.nmass_leaf *= frac_remain;
		}
		leaf_to_litter(soil, pft, closs, nloss);

		// TODO: rationalise this. We don't really want to calculate LAI in lots
		// of different places.
		indiv.lai *= frac_remain;

		// Accounting.
		indiv.closs_temp = closs;
		indiv.nloss_temp = nloss;
		indiv.loss_temp = frac_loss;
		indiv.acloss_temp += closs;
		indiv.anloss_temp += nloss;
		indiv.aloss_temp += frac_loss;
	}
}

void canopy_exchange_dave(Patch &patch) {

	partition_canopy_dave(patch);

    // Get vcmax25/jmax25/vpmax25 values.
    store_vcmax_values(patch.vegetation);

    // Calculate soil water potential
    calc_soil_water_potenial_dave(patch.soil);

    // Calculate water stress factor
	water_stress_dave(patch);

    // Calculate subdaily photosynthesis and stomatal conductance.
    photosynthesis_gs(patch);

    // Calculation daily transpiration and fractional uptake by layer.
    water_uptake_dave(patch);

    // Calculate respiration and NPP on a daily basis.
    npp_dave(patch);

	// Calculate leaf damage due to high temperature.
	leaf_damage_dave(patch);

    // Leaf senescence for crops.
    leaf_senescence(patch.vegetation);

    // Nitrogen demand for cohorts NOT using DAVE allocation.
    ndemand(patch, patch.vegetation);

    if (!ifdavegrowth) {
        // Nitrogen competition and ~uptake~. Well, actually one could argue
        // that this isn't the uptake, because the nitrogen is removed from
        // the soil later (in somdynam::vegetation_n_uptake()). But this
        // function certainly calculates how much each indiv will uptake,
        // taking into account factors such as competition, availability,
        // stress, etc.
        nitrogen_uptake(patch.vegetation);
    }

    forest_floor_dave(patch);
}
