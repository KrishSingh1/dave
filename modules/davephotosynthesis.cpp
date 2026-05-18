#include "davephotosynthesis.h"
#include "hydraulics.h"
#include "guess.h"
#include "ncutil.h"
#include "canexch.h" // Many constants defined here
#include "logging.h"
#include "guessmath.h"
#include "driver.h"

// Average leaf N:P ratio for C3 grasses from Thum et al. 2019 (SI Table 8).
const double C3_LEAF_NTOP = 10.7;

// Number of soil layers to use for calculating PFT-level water stress scalar.
const int NLAYER_FOREST_FLOOR = 5;

// Scaling parameter for gcutic (Cochard 2021 Eq 13 & 14).
const double Q10A = 1.2;

// Scaling parameter for gcutic (Cochard 2021 Eq 13 & 14).
const double Q10B = 4.8;

// Transition temperature degC (to discuss with Erick) (Cochard 2021 Eq 13 & 14)
const double TP = 40;

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

class PhotosynthesisDiagnostics {
public:
	double g1;
	double g0;
	double gcutic_adj;
};

/*
Describes the state of a cohort relevant for photosynthesis calculations.
*/
class CohortState {
public:
	CohortState(double lai, double wscal, double lwp,
				double fipar, bool stomata_open, double fphen_activity,
				double ps_sink_lim, double vcmax25, double jmax25,
				double vpmax25) :
		lai(lai), wscal(wscal), lwp(lwp), fipar(fipar),
		stomata_open(stomata_open), fphen_activity(fphen_activity),
		ps_sink_lim(ps_sink_lim), vcmax25(vcmax25), jmax25(jmax25),
		vpmax25(vpmax25) {
	}

	CohortState(const Individual &indiv, size_t layer) : CohortState(
		indiv.lai_layerwise[layer],
		indiv.get_wscal(),
		ifdavehydraulics ? indiv.predawn_lwp : indiv.weighted_psi,
		indiv.fipar_layerwise[layer],
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
double calc_vcmax25(double nmass_leaf, double pmass_leaf) {
	return exp(4.308 + 0.298 * log(pmass_leaf * G_PER_KG) + 0.197 *
							   log(nmass_leaf * G_PER_KG));
}

/*
Calculate maximum electron transport rate at 25 degC (μmol m-2 s-1).

@param nmass_leaf: Leaf N content at the top of the canopy (kgN m-2 per green leaf area).
@param pmass_leaf: Leaf P content at the top of the canopy (kgP m-2 per green leaf area).
*/
double calc_jmax25(double nmass_leaf, double pmass_leaf) {
	return exp(5.139 + 0.325 * log(pmass_leaf * G_PER_KG) + 0.112 *
							   log(nmass_leaf * G_PER_KG));
}

/*
Calculate maximum PEP carboxylation rate at 25 degC (μmol m-2 s-1).

@param nmass_leaf: Leaf N content at the top of the canopy (kgN m-2 per green leaf area).
*/
double calc_vpmax25(double nmass_leaf) {
	return exp(3.118 + 0.751 * log(nmass_leaf * G_PER_KG));
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

			indiv.ps_subdaily.vcmax25 = calc_vcmax25(nmass_leaf0, pmass_leaf0);
			indiv.ps_subdaily.jmax25 = calc_jmax25(nmass_leaf0, pmass_leaf0);
		} 
		else if (indiv.pft.pathway == C4) {

			indiv.ps_subdaily.vpmax25 = calc_vpmax25(nmass_leaf0);

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

/*
Calculate water stress factor [0-1] in the root zone. This is calculated as a
weighted mean of wcont / wcont_crit (clamped to [0, 1]) across soil layers.
*/
double calc_water_stress_dave(const Soil &soil, const double weight[NSOILLAYER], size_t max_layers = NSOILLAYER) {
	if (max_layers > NSOILLAYER) {
		fail("calc_water_stress_dave(): max_layers cannot be greater than NSOILLAYER");
	}

	double wstress = 0;
	double total_weight = 0;
	for (size_t i = 0; i < max_layers; i++) {
		double stress_base = soil.get_layer_soil_water(i) / wcont_crit;

		// Clamp to [0, 1].
		stress_base = max(min(stress_base, 1.0), 0.0);

		wstress += stress_base * weight[i];

		// Root distribution should add up to 1 but summed here just in case.
		total_weight += weight[i];
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
		double weight[NSOILLAYER];
		for (size_t j = 0; j < NSOILLAYER; j++) {
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
												NLAYER_FOREST_FLOOR);
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
double stomatal_conductance_term(double wscal, double g1, double vpd, double ca) {
	return divide(divide(1.0 + pow(wscal, qs) * g1, sqrt(max(vpd, 0.1))), ca);
}

/*
Residual stomatal conductance term.
Note: g0 is for H20, here conversion to conductances to CO2.
TODO: do we need a different qs here?
*/
double stomatal_conductance_residual(double g0, double wscal) {
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
double electron_transport_rate(double apar,double jmax) {
	// Photosynthetic photon flux density (mol (photons) m-2 s-1).
	double ppfd = apar * JTOMOL;

	double a = THETA;
	double b = -(ppfd + jmax);
	double c = ppfd * jmax;

	double j;
	if (!solve_quadratic(a, b, c, j, true)) {
		log_error("electron_transport_rate: No solution for quadratic J. Returning J=0 to avoid crash. apar=%g, jmax=%g", apar, jmax);
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

	// Daytime (mitochondrial) leaf respiration (mol m-2 s-1)
	// Note: fphen_activity is always 1 for trees.

	// vcmax25/jmax25 are input in μmol m-2 s-1, but the equations expect them
	// in mol m-2 s-1.
	double vcmax25_mol = inputs.vcmax25 / UMOL_PER_MOL;
	double jmax25_mol  = inputs.jmax25  / UMOL_PER_MOL;
	double resp_scalar = resp_temp_acclimation_scalar(inputs.tgrowth, inputs.racclim);

	double rd = inputs.fphen_activity * pow(inputs.wscal, qb) * BC3 * vcmax25_mol * inputs.nscal *
			respiration_temp_response(inputs.tleaf) * light_inhibition(inputs.ipar) * resp_scalar;

	if (is_zero(inputs.ipar)) {
		no_photosynthesis(inputs, outputs, rd);
		return;
	}

	// Michaelis-Menten parameters and photorespiratory CO2 compensation point (mol mol-1).
	double ko = Arrhenius_temp_response(KO25, inputs.tleaf, HA_KO, 0.0, 0.0);
	double kc = Arrhenius_temp_response(KC25, inputs.tleaf, HA_KC, 0.0, 0.0);
	double gamstar = Arrhenius_temp_response(GAMSTAR25, inputs.tleaf, HA_GS, 0.0, 0.0);

	// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax.
	double pslim = inputs.fphen_activity * inputs.ps_sink_lim * inputs.nscal;

	// apply water stress scaling to vcmax / jmax
	if (ifzhou) {
	// function to scale vcmax and jmax by pre-dawn water potential, see Zhou (2014)
		double wp_vmax_scale = pslim_lwp_scalar(inputs.s_f, inputs.psi_f, inputs.predawn_lwp);
		pslim = pslim * wp_vmax_scale;
	} else {
		pslim = pslim * pow(inputs.wscal, qb);
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
	double je = electron_transport_rate(inputs.apar, jmax);

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

	// kc * (1 + om / ko)
	double km = biochemistry.km;

	// TODO: check value of alpha_PSII.

	// Equations 22-24
	double a_c = alpha_PSII / a0 * kc / ko;
	double b_c = -((vp - rm + gbs * ci) + (vcmax - rd) + gbs * km + alpha_PSII / a0
			* (lgamstar * vcmax + rd * kc / ko));
	double c_c = (vcmax - rd) * (vp - rm + gbs * ci) -
			(vcmax * gbs * lgamstar * O2atm + rd * gbs * km);

	double anc;
	if (solve_quadratic(a_c, b_c, c_c, anc, true)) {
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
	double je = electron_transport_rate(inputs.apar, biochemistry.jmax);

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
	double je = electron_transport_rate(inputs.apar, biochemistry.jmax);

	double rm = biochemistry.rm;
	double ci = state.ci;
	double rd = biochemistry.rd;

	// Equations 41-43
	double a_j = 1.0 - (7.0 * lgamstar * alpha_PSII) / (3.0 * a0);
	double b_j = -((divide(zj,2.0) * xj * je - rm + gbs * ci) + (divide(zj,3.0) * (1.0 - xj) * je - rd) +
			gbs * (divide((7.0 * lgamstar * O2atm),3.0)) + divide((alpha_PSII * lgamstar),a0) *
			(divide(zj,3.0) * (1.0 - xj) * je + divide((7.0 * rd),3.0)));
	double c_j = ((divide(zj,2.0) * xj * je) - rm + gbs * ci) * (divide(zj,3.0) * (1.0 - xj) * je - rd) -
			gbs * lgamstar * O2atm * (divide(zj,3.0) * (1.0 - xj) * je + divide((7.0 * rd),3.0));

	double anj;
	if (solve_quadratic(a_j, b_j, c_j, anj, true)) {
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
	result.rd = inputs.fphen_activity * pow(inputs.wscal, qb) * BC4 * vcmax25_c4 * inputs.nscal *
			respiration_temp_response(inputs.tleaf) * light_inhibition(inputs.ipar) * resp_scalar;

	// Rd in the mesophyll cells.
	result.rm = 0.5 * result.rd;

	// Rd in the bundle-sheath cells.
	result.rs = result.rd - result.rm;

	// Nitrogen and water scaling and temp. dependencies of Vcmax and Jmax
	double pslim = inputs.fphen_activity * inputs.ps_sink_lim * inputs.nscal;

	// apply water stress scaling to vcmax / jmax
	if (ifzhou) {
		// function to scale vcmax and jmax by pre-dawn water potential, see Zhou (2014)
		double wp_vmax_scale = pslim_lwp_scalar(inputs.s_f, inputs.psi_f, inputs.predawn_lwp);
		pslim = pslim * wp_vmax_scale;
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

	if (is_zero(inputs.ipar) || !inputs.stomata_open) {
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

	// Predawn LWP will hold yesterday's value for the first few time
	// steps, until we hit the dawn-transition timestep. Therefore it's
	// important to set this inside the subdaily loop.
	inputs.predawn_lwp = cohort.get_lwp();

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
		// Stomatal term (Medlyn et al. 2011), with g1 defined by
		// pre-dawn soil water potential (Zhou et al. 2014).

		// We set parameter b1z to the value that results in 90% reduction in g1
		// at stom_clos (for compatibility with hydraulic module).
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
		inputs.X = stomatal_conductance_term(inputs.wscal, pft.g1, inputs.vpd, inputs.ca);

		// Residual stomatal conductance term
		// Note: g0 is for H20, here conversion to conductances to CO2
		// TODO: do we need a different qs here?
		inputs.gdark = stomatal_conductance_residual(pft.g0, inputs.wscal);
	}

	// LAI of this individual in this layer.
	// Note: this takes phen into account for indivs running trunk
	// allocation, because fpar() looks at indiv.lai_today().
	inputs.lai = cohort.get_lai();

	// We can actually have a canopy layer with 0 indiv-level LAI
	// but non-zero fipar. This is because fipar is not a truly
	// indiv-level property, as it is calculated using canopy-level
	// layered LAI.
	//
	// Question: can we skip this layer if lai_layer is 0?
	assert(inputs.lai >= 0);

	double par = max(0.0, climate.pars[day.period]);
	inputs.ipar = cohort.get_fipar() * par;
	inputs.apar = cohort.get_fipar() * par * pft.alpha * pft.absorptance;

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
@param patm		 Atmospheric pressure (kPa) - default 101.5.
@param wind		 Windspeed (m s-1) - default 2.5.
@param wleaf	 Leaf width (m) - PFT specific.
@param stomatal_ratio ???
@param leaf_abs	 Leaf absorptance (fraction).
@return Leaf temperature (°C).
*/
double energy_balance_dave(double tleaf, double tair, double gs, double ppfd,
				           double vpd, double patm, double wind, double wleaf,
					       double stomatal_ratio, double leaf_abs) {
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
	double gradiation = 4.0 * BOLTZ * pow(tair_k, 3.0) * EMISSIVITY / (CP_AIR * AIRMA);
	double gbhforced = 0.003 * sqrt(wind / wleaf) * cmolar ;
	double grashof = 1.6e+08 * fabs(tleaf - tair) * pow(wleaf, 3.0);
	double gbhfree = 0.5 * DHEAT * pow(grashof, 0.25) / wleaf * cmolar;
	double gbh = 2.0 * (gbhfree + gbhforced);
	double gbhr = gbh + 2.0 * gradiation;
	double gbw = stomatal_ratio * 1.075 * gbh;
	double gw = gs * gbw / (gs + gbw);

	// Radiation.
	// double Rlongup = EMISSIVITY * BOLTZ * pow(Tk(tleaf), 4);
	double rsol = 2.0 * leaf_abs * ppfd / UMOL_PER_J; // J m-2 s-1 - absorbed SW radiation = 2 * APAR
	// double Rnet = rsol - Rlongup;
	double ea = 1000.0 * (esat(tair, m) - vpd); // In Pa. If this is negative we have a problem
	if (ea < 0) {
		fail("energy_balance_dave(): ea < 0, this means saturation vapour pressure is less than actual vapour pressure, which shouldn't happen. east=%f kPa, vpd=%f kPa, tair=%f °C",
			 esat(tair, m), vpd, tair);
	}

	double ema = 0.642 * pow(ea / tair_k, 1.0 / 7.0);
	double rnetiso = rsol - (1.0 - ema) * BOLTZ * pow(tair_k, 4.0);

	// # Penman-Monteith equation
	double gamma = CP_AIR * AIRMA * patm * 1000.0 / lhv;

	// Transpiration rate (mol H2O m-2 s-1) : in place of gs * VPD
	double et = (1.0 / lhv) * (slope * rnetiso + 1000.0 * vpd * gbh * CP_AIR *
					AIRMA) / (slope + gamma * gbhr / gw);

	// Calculate leaf temperature.
	double lambda_et = lhv * et;
	double y = 1.0 / (1.0 + gradiation / gbh);
	double h2 = y * (rnetiso - lambda_et);
	// double H = -CP_AIR * air_density * (gbh / cmolar) * (tair - tleaf);
	double tleaf2 = tair + h2 / (CP_AIR * air_density * (gbh / cmolar));
	assert(tleaf2 == tleaf2); // check for NaNs
	return tleaf2;
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
							   const Pft &pft,
						       PhotosynthesisOutputs &outputs,
							   const PhotosynthesisDiagnostics &diagnostics,
						       pathwaytype pathway,
							   size_t niter_max,
							   double epsilon) {
    if (niter_max == 0) {
		fail("solve_energy_balance_dave(): niter_max must be greater than 0.");
	}

    // TODO: figure out what this is.
    const double stomatal_ratio = 1.0;

	// Initialise guess: leaf temperature = air temperature.
	inputs.tleaf = climate_inputs.get_tair();

	bool converged = false;

	// Guaranteed to be set.
	size_t i;
	double tleaf_diff;

    for (i = 0; i < niter_max && !converged; i++) {
		outputs.clear();
		photosynthesis_dave(inputs, outputs, pathway);

		double tleaf_old = inputs.tleaf;
		double gsw = get_gsw(outputs.gsc, diagnostics.gcutic_adj);
		double wind = max(climate_inputs.get_wind(), 0.01);

		inputs.tleaf = energy_balance_dave(inputs.tleaf,
										   climate_inputs.get_tair(),
									       gsw,
										   inputs.ipar,
									       climate_inputs.get_vpd(),
									       climate_inputs.get_patm(),
									       wind,
									       pft.leaf_width,
										   stomatal_ratio,
									       pft.absorptance);

        tleaf_diff = fabs(inputs.tleaf - tleaf_old);
		converged = tleaf_diff <= epsilon;
	}

	// This is not the same as the inputs.tleaf used to generate these
	// photosynthesis outputs, but rather it's the tleaf that results from them.
	// This reproduces old no-energy-balance code when niter_max = 1.
	outputs.tleaf = inputs.tleaf;
	EnergyBalance result(inputs.tleaf, i, tleaf_diff, converged);
	return result;
}

/*
Coupled model of photosynthesis and stomatal conductance for C3 and C4 plants.
*/
void photosynthesis_gs(Patch& patch, Climate& climate, Vegetation& vegetation) {
	PhotosynthesisInputs inputs;
	PhotosynthesisOutputs outputs;

	// Loop through cohorts (=individuals)
	for (int i = 0; i < vegetation.nobj; i++) {

		// Retrieve information about cohort
		Individual& indiv = vegetation[i];
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (is_zero(indiv.get_nmass_leaf()) || is_zero(indiv.get_wscal())) {
			continue;
		}

		// TODO: find another home for this check. It doesn't belong here.
		double lai_layer_sum = sum(indiv.lai_layerwise);
		if (!equal(lai_layer_sum, indiv.lai_today())) {
			log_warning("Individual %d (%s): layerwise LAI sum (%g) does not equal total LAI (%g) (phen=%g)",
						i, (const char *)pft.name, lai_layer_sum,
						indiv.lai_today(), indiv.phen);
		}

		indiv.photosynthesis.agd_g = 0.0;
		indiv.photosynthesis.rd_g = 0.0;
		indiv.photosynthesis.je = 0.0;

		double an_total = 0;
		indiv.anc_frac = 0;
		indiv.anj_frac = 0;
		indiv.anp_frac = 0;

		indiv.an_layerwise.clear();
		indiv.an_layerwise.resize(indiv.lai_layerwise.size());

		for (int j = 0; j < indiv.lai_layerwise.size(); j++) {
			indiv.an_layerwise[j] = 0;
		}

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
			size_t nlayer = indiv.fipar_layerwise.size();

            // Initialise the layerwise inputs and outputs for this timestep.
            // This can't happen earlier because we need to know the number of
            // layers, which isn't known until we call fpar().
			indiv.ps_inputs_layerwise[day.period].assign(nlayer, PhotosynthesisInputs{});
			indiv.ps_outputs_layerwise[day.period].assign(nlayer, PhotosynthesisOutputs{});
			indiv.energy_balance_layerwise[day.period].assign(nlayer, EnergyBalance(0, 0, 0, false));

		    for (size_t j = 0; j < nlayer; j++) {

				// No LAI means we can skip this layer as an optimisation. This
				// can happen because the number of layers is the total number
				// of canopy layers. So if this tree is shorter than the tallest
				// tree in the patch, it may not have any LAI in this layer.
				if (is_zero(indiv.lai_layerwise[j])) {
					continue;
				}

				// Determine, for this indiv, LAI and cumulative LAI of the layer.

				// Cumulative LAI at the bottom of current layer.
				lai_cum_bottom += indiv.lai_layerwise[j];

				// Cumulative LAI at the middle of current layer.
				double lai_cum = lai_cum_bottom - 0.5 * indiv.lai_layerwise[j];

				CohortState cohort(indiv, j);

				PhotosynthesisDiagnostics diagnostics = get_photosynthesis_inputs(
					cohort, pft, climate, patch, day, lai_cum, inputs);
				indiv.g1 = diagnostics.g1;

				ClimateInputs clim(climate, day);

				EnergyBalance balance = solve_energy_balance_dave(inputs,
					clim, pft, outputs, diagnostics, pft.pathway,
					energy_balance_maxiter, energy_balance_eps);

				if (balance.get_converged()) {
					log_diagnostic("Energy balance solver converged for individual %d (%s) in layer %d at day period %d in %d iterations.",
								   indiv.id, (const char *)pft.name, j, day.period, balance.get_niter());
				} else {
					log_warning("Energy balance solver did not converge for individual %d (%s) in layer %d at day period %d after %d iterations (tleaf_diff = %g).",
								indiv.id, (const char *)pft.name, j, day.period,
								balance.get_niter(), balance.get_delta_tleaf());
				}

				indiv.ps_inputs_layerwise[day.period][j] = inputs;
				indiv.ps_outputs_layerwise[day.period][j] = outputs;
				indiv.energy_balance_layerwise[day.period][j] = balance;

				// TODO drew: multiply by LAI only once, after canopy loop.

				// sum variables over canopy layers
				double lai_layer = inputs.lai;

				indiv.an_layerwise[j] += outputs.an;

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
			double lai_tot = sum(indiv.lai_layerwise);

			// Calculate mean for some variables.
			if (lai_tot == 0) {
				log_warning("Individual %d (%s): Total layerwise LAI is zero (%g) (indiv.lai = %g)",
					        indiv.id, (const char *)pft.name, lai_tot,
							indiv.lai);
			} else {
				ci_canopy /= lai_tot;
				vcmax_canopy /= lai_tot;
				jmax_canopy /= lai_tot;
			}

			assert(indiv.lai < 20);

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
		indiv.photosynthesis.nactive_opt = 0.0;    // not used
		indiv.photosynthesis.vmaxnlim    = 0.0;    // not used 
	} // end cohort loop
}

double calc_water_uptake_dave(Individual &indiv) {
    Patch &patch = indiv.vegetation.patch;
    Pft &pft = indiv.pft;
    const Standpft& spft = patch.stand.pft[indiv.pft.id];

    if (spft.irrigated && pft.id == patch.stand.pftid) {
        Day day;
        return irrigated_water_uptake(patch, pft, day, indiv);
    }

    // Available water in each layer, as a fraction of water holding
    // capacity.
    double wcont_local[NSOILLAYER];
    for (int i = 0; i<NSOILLAYER; i++) {
        wcont_local[i] = patch.soil.get_layer_soil_water(i);
    }

    if (iftwolayersoil) {
        return water_uptake_twolayer(wcont_local, patch.soil.soiltype.awc,
                    indiv.rootdist, pft.emax, patch.fpc_rescale, indiv.fwuptake,
                    pft.lifeform == TREE, pft.drought_tolerance);
    }

    return water_uptake(wcont_local, patch.soil.soiltype.awc,
                indiv.rootdist, pft.emax, patch.fpc_rescale, indiv.fwuptake,
                pft.lifeform == TREE, pft.drought_tolerance, indiv.Rs_r);
}

/*
Calculate canopy transpiration.
*/
void transpiration_dave(Patch &patch, Climate& climate, Vegetation& vegetation) {

	// this formulation assumes infinite Ga
	// TODO: include Ga --> Penman-Monteith (?)

	double wsupply = 0;

	// Available water, by layer, in mm.
	double wcont[NSOILLAYER];
	for (int i = 0; i < NSOILLAYER; i++) {
		wcont[i] = patch.soil.get_layer_soil_water(i) * patch.soil.soiltype.awc[i];
		wsupply += wcont[i];
	}

	// Loop through cohorts (=individuals)
	for (int i = 0; i < vegetation.nobj; i++) {

		// Retrieve information about cohort
		Individual& indiv = vegetation[i];
		Pft& pft = indiv.pft;
		Patchpft &ppft = patch.pft[indiv.pft.id];

		if (date.day == 0)
			indiv.aaet = 0.0;

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
		indiv.aet = transp * MOL2G_H2O * KG_PER_G * SECS_PER_DAY; // assuming 1kg m-2 = 1l m-2 == 1mm

		// Note: aet may be further limited by the water supply (see below).
		// This also means that we can't update aaet here (also see below).

		// Note that in standard LPJ-GUESS transpiration is described as evapotranspiration (AET),
		// which is incorrect as soil evaporation is calculated separately.

		// Limit water uptake.
		double wr = calc_water_uptake_dave(indiv);

        if (ifdavehydraulics && !indiv.stomata_open) {
            for (int i = 0; i<NSOILLAYER; i++) {
                // If the stomata are closed, we assume roots are hydraulically
                // isolated from the soil.
                indiv.fwuptake[i] = 0;
            }
        }
	} // end cohort loop

	// Ensure water uptake doesn't exceed availability in any layer.
	for (int i = 0; i < NSOILLAYER; i++) {
		double uptake = 0;
		for (int j = 0; j < vegetation.nobj; j++) {
			const Individual &indiv = vegetation[j];
			const Patchpft &ppft = indiv.patchpft();
			uptake += indiv.aet * indiv.fwuptake[i]; 
		}
		if (uptake > wcont[i]) {
			if (ifdavehydraulics) {
				log_error("Transpiration exceeds soil water content, this shouldn't happen\n");
			}

			// Demand in this layer exceeds supply.
			double fwavail = 0.999 * wcont[i] / uptake;
			for (int j = 0; j < vegetation.nobj; j++) {
				Individual &indiv = vegetation[j];
				indiv.aet *= fwavail;
			}
		}
	}

	// aaet updated here to allow for reduction in aet due to source limitation.
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];
		indiv.aaet += indiv.aet;
	}
}

struct CtonRatios {
	double sap;
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
		cton.root  = indiv.cton_root();
		cton.crown = indiv.cton_crown();
	} else {
		cton.sap   = indiv.pft.cton_sap_avr;
		cton.root  = indiv.pft.cton_root_avr;
		// Crowns use root stoichiometric limits.
		cton.crown = indiv.pft.cton_root_avr;
	}

	if (!indiv.pft.dave_allocation()) {
		// Set C:N ratios if not using DAVE allocation (since trunk
		// allocation + DAVE photosynthesis can give unrealistic values).
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
void npp_dave(Patch& patch, Climate& climate, Vegetation& vegetation) {

	// Response of respiration rate to temperature
	double gtemp = climate.get_gtemp_mean();

	for (size_t i = 0; i < vegetation.nobj; i++) {
		Individual& indiv = vegetation[i];
		Pft& pft = indiv.pft;
		Patchpft& ppft = patch.pft[pft.id];

		// Don't do calculations for crops outside their growing season.
		if (!indiv.growingseason()) {
			indiv.dnpp=0.0;
			indiv.lambda = 0;
			continue;
		}

		// get daily gross assimilation (kgC/m2)
		double assim = indiv.photosynthesis.gross_assimilation();

		// Calculate autotrophic respiration

		double cmass_root = indiv.cmass_root_today();

		// 0 for trees
		double cmass_crown = indiv.get_cmass_crown();

		// 0 for PFTs using annual (trunk) allocation
		double cmass_repr = indiv.get_cmass_repr();

		// Get effective C:N ratios for sap, root, and crown.
		CtonRatios cton = cton_effective(indiv);

		// Set leaf respiration.
		indiv.resp_leaf = indiv.photosynthesis.rd_g * KG_PER_G;

		// 30 day moving average of air temperature.
		double tgrowth = patch.get_climate().get_tgrowth();

        // Get scalar encapsulating the effect of temperature acclimation. This
        // will be 1 if acclimation is disabled.
        double scalar = resp_temp_acclimation_scalar(tgrowth, pft.racclim);

		// Respiration coefficient is reduced based on temperature acclimation.
        double respcoeff = pft.respcoeff * scalar;

		double resp = 0;
		respiration(gtemp, patch.soil.gtemp, indiv.pft.lifeform,
			respcoeff, cton.sap, cton.root, cton.crown, indiv.cmass_sap,
			cmass_root, cmass_crown, cmass_repr, indiv.cgrow, indiv.fphen_activity,
			assim, indiv.resp_leaf, resp, indiv.resp_maintenance,
			indiv.resp_growth, indiv.resp_repr, indiv.resp_crown,
			indiv.resp_root, indiv.resp_sap);

		// Update accumulated annual NPP and daily vegetation-atmosphere flux
		// Note that resp does not include leaf day respiration.
		indiv.dgpp = assim;
		indiv.dnpp = assim - resp;
		indiv.anpp += indiv.dnpp;
		indiv.agpp	+= assim;
		indiv.aresp	+= resp;

		// Save historical GPP for daily growth of grass
		indiv.gpp_15.add(assim);

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
			double lwp = mean(patch.soil.psi_layer, NLAYER_FOREST_FLOOR);

			// No soil water on first day of simulation - this results in wscal
			// of 0, which results in a division by zero. We can just skip this
			// calculation on the first day of the simulation.
			if (is_zero(ppft.wscal)) {
				continue;
			}

			CohortState cohort(LAI_HYPOTHETICAL, ppft.wscal, lwp, patch.fpar_ff,
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

			Standpft &spft = patch.stand.pft[pft.id];
			if (day.isend && date.islastday && date.islastmonth) {
				spft.anetps_ff_max = max(spft.anetps_ff_max, ppft.anetps_ff);
				ppft.anetps_ff = max(0.0, ppft.anetps_ff);
			}
		}
	}
}

void canopy_exchange_dave(Patch &patch, Climate &climate) {

    // Get vcmax25/jmax25/vpmax25 values.
    store_vcmax_values(patch.vegetation);

    // Calculate soil water potential
    calc_soil_water_potenial_dave(patch.soil);

    // Calculate water stress factor
	water_stress_dave(patch);

    // Calculate subdaily photosynthesis and stomatal conductance 
    photosynthesis_gs(patch, climate, patch.vegetation);

    // Transpiration from VPD and canopy conductance
    transpiration_dave(patch, climate, patch.vegetation);

    // Calculate respiration and NPP on a daily basis
    npp_dave(patch, climate, patch.vegetation);

    // needed here?
    leaf_senescence(patch.vegetation);

    // Nitrogen demand
    // This will only affect PFTs:
    // - Grasses if using trunk allocation
    // - Trees always
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
