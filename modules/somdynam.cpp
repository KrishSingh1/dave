///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.cpp
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith (LPJ SOM dynamics, CENTURY), David Wårlind (CENTURY)
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
#include "somdynam.h"
#include "ntransform.h"
#include "driver.h"
#include <assert.h>
#include <bitset>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Turnover times (in years, approximate) for litter and SOM fractions at 10 deg C with
// ample moisture (Meentemeyer 1978; Foley 1995)

static const double TAU_LITTER=2.85; // Thonicke, Sitch, pers comm, 26/11/01
static const double TAU_SOILFAST=33.0;
static const double TAU_SOILSLOW=1000.0;	

static const double FASTFRAC=0.985;
	// fraction of litter decomposition entering fast SOM pool
static const double ATMFRAC=0.7;
	// fraction of litter decomposition entering atmosphere

// Corresponds to the amount of soil available nitrogen where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
// Comment: NMASS_SAT is too high when considering BNF - Zaehle
static const double NMASS_SAT = 0.002 * 0.05;
// Corresponds to the nitrogen concentration in litter where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
static const double NCONC_SAT = 0.02;

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

// Exponential decay constants for litter and SOM fractions
// Values set from turnover times (constants above) on first call to decayrates

static double k_litter10;
static double k_soilfast10;
static double k_soilslow10;

static bool firsttime=true;
	// indicates whether function decayrates has been called before


///////////////////////////////////////////////////////////////////////////////////////
// SETCONSTANTS
// Internal function (do not call directly from framework)

void setconstants() {

	// DESCRIPTION
	// Calculate exponential decay constants (annual basis) for litter and
	// SOM fractions first time function decayrates is called

	k_litter10=1.0/TAU_LITTER;
	k_soilfast10=1.0/TAU_SOILFAST;
	k_soilslow10=1.0/TAU_SOILSLOW;
	firsttime=false;
}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Internal function (do not call directly from framework)
// used by som_dynamic_lpj()

void decayrates(double wcont, double gtemp_soil, double& k_soilfast, double& k_soilslow,
	double& fr_litter, double& fr_soilfast, double& fr_soilslow, bool tillage) {

	// DESCRIPTION
	// Calculation of fractional decay amounts for litter and fast and slow SOM
	// fractions given current soil moisture and temperature

	// INPUT PARAMETERS
	// wcont       = water content of upper soil layer (fraction of AWC)
	// gtemp_soil  = respiration temperature response incorporating damping of Q10
	//               response due to temperature acclimation (Eqn 11, Lloyd & Taylor
	//               1994)

	// OUTPUT PARAMETERS
	// k_soilfast  = adjusted daily decay constant for fast SOM fraction
	// k_soilslow  = adjusted daily decay constant for slow SOM fraction
	// fr_litter   = litter fraction remaining following today's decomposition
	// fr_soilfast = fast SOM fraction remaining following today's decomposition
	// fr_soilslow = slow SOM fraction remaining following today's decomposition

	double moist_response; // moisture modifier of decomposition rate

	// On first call only: set exponential decay constants

	if (firsttime) setconstants();

	// Calculate response of soil respiration rate to moisture content of upper soil layer
	// Foley 1995 Eqn 19

	moist_response = 0.25 + 0.75 * wcont;

	// Calculate litter and SOM fractions remaining following today's decomposition
	// (Sitch et al 2000 Eqn 71) adjusting exponential decay constants by moisture and
	// temperature responses and converting from annual to daily basis
	// NB: Temperature response (gtemp; Lloyd & Taylor 1994) set by framework

	k_soilfast = k_soilfast10 * gtemp_soil * moist_response / (double)date.year_length();
	if (tillage) {
		k_soilfast *= TILLAGE_FACTOR; // Increased HR for crops (tillage)
	}
	k_soilslow = k_soilslow10 * gtemp_soil * moist_response / (double)date.year_length();

	fr_litter = exp(-k_litter10 * gtemp_soil * moist_response / (double)date.year_length());
	fr_soilfast = exp(-k_soilfast);
	fr_soilslow = exp(-k_soilslow);
}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Should be called by framework on last day of simulation year, following call to
// som_dynamics, once annual litter production and vegetation PFT composition are close
// to their long term equilibrium (typically 500-1000 simulation years).
// NB: should be called ONCE ONLY during simulation for a particular grid cell

void equilsom_lpj(Soilcolumn& soilcolumn) {

	// DESCRIPTION
	// Analytically solves differential flux equations for fast and slow SOM pools
	// assuming annual litter inputs close to long term equilibrium

	// INPUT PARAMETER (class defined in framework header file)
	// soil = current soil status

	// number of years over which decay constants and litter inputs averaged
	double nyear;
		
	nyear = soilcolumn.soiltype.solvesom_end - soilcolumn.soiltype.solvesom_begin + 1;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		soillayer.decomp_litter_mean /= nyear;
		soillayer.k_soilfast_mean /= nyear;
		soillayer.k_soilslow_mean /= nyear;

		soillayer.cpool_fast = (1.0 - ATMFRAC) * FASTFRAC * soillayer.decomp_litter_mean /
			soillayer.k_soilfast_mean;
		soillayer.cpool_slow = (1.0 - ATMFRAC) * (1.0 - FASTFRAC) * soillayer.decomp_litter_mean /
			soillayer.k_soilslow_mean;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// SOM DYNAMICS
// To be called each simulation day for each modelled area or patch, following update
// of soil temperature and soil water.

void som_dynamics_lpj(Patch& patch, bool tillage) {

	// DESCRIPTION
	// Calculation of soil decomposition and transfer of C between litter and soil
	// organic matter pools.

	double k_soilfast; // adjusted daily decay constant for fast SOM fraction
	double k_soilslow; // adjusted daily decay constant for slow SOM fraction
	double fr_litter;
		// litter fraction remaining following one day's/one month's decomposition
	double fr_soilfast;
		// fast SOM fraction remaining following one day's/one month's decomposition
	double fr_soilslow;
		// slow SOM fraction remaining following one day's/one month's decomposition
	double decomp_litter; // litter decomposition today/this month (kgC/m2)
	double cflux; // accumulated C flux to atmosphere today/this month (kgC/m2)
	int p;
	bool tillage_lyr;

	// Obtain reference to Soil object
	Soilcolumn& soilcolumn = patch.soilcolumn;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Tillage in this layer
		if (lyr == 0 && tillage)
			tillage_lyr = true;
		else
			tillage_lyr = false;

		// Calculate decay constants and rates given today's soil moisture and
		// temperature

		decayrates(soillayer.wcont(), soillayer.gtemp, k_soilfast, k_soilslow, fr_litter,
			fr_soilfast, fr_soilslow, tillage_lyr);

		// From year soilcolumn.solvesom_begin, update running means for later solution
		// (at year soilcolumn.solvesom_end) of equilibrium SOM pool sizes

		if (date.year >= soilcolumn.solvesom_begin) {
			soillayer.k_soilfast_mean += k_soilfast;
			soillayer.k_soilslow_mean += k_soilslow;
		}

		// Reduce litter and SOM pools, sum C flux to atmosphere from decomposition
		// and transfer correct proportions of litter decomposition to fast and slow
		// SOM pools

		// Reduce individual litter pools and calculate total litter decomposition
		// for today/this month

		decomp_litter = 0.0;

		// Loop through PFTs

		for (p = 0; p < npft; p++) {	//NB. also inactive pft's

			// For this PFT ...

			if (lyr == 0) {

				decomp_litter += (patch.pft[p].litter_leaf +
					patch.pft[p].litter_sap +
					patch.pft[p].litter_heart +
					patch.pft[p].litter_repr) * (1.0 - fr_litter);

				patch.pft[p].litter_leaf *= fr_litter;
				patch.pft[p].litter_sap *= fr_litter;
				patch.pft[p].litter_heart *= fr_litter;
				patch.pft[p].litter_repr *= fr_litter;
			}

			decomp_litter += patch.pft[p].litter_root * (1.0 - fr_litter) * patch.pft[p].rootdist[lyr];

			if (lyr == soilcolumn.nsoillayers - 1)
				patch.pft[p].litter_root *= fr_litter;
		}

		if (date.year >= soilcolumn.soiltype.solvesom_begin)
			soillayer.decomp_litter_mean += decomp_litter;

		// Partition litter decomposition among fast and slow SOM pools
		// and flux to atmosphere

		// flux to atmosphere
		cflux = decomp_litter * ATMFRAC;

		// remaining decomposition - goes to ...
		decomp_litter -= cflux;

		// ... fast SOM pool ...
		soillayer.cpool_fast += decomp_litter * FASTFRAC;

		// ... and slow SOM pool
		soillayer.cpool_slow += decomp_litter * (1.0 - FASTFRAC);

		// Increment C flux to atmosphere by SOM decomposition
		cflux += soillayer.cpool_fast * (1.0 - fr_soilfast) + soillayer.cpool_slow * (1.0 - fr_soilslow);

		// Reduce SOM pools
		soillayer.cpool_fast *= fr_soilfast;
		soillayer.cpool_slow *= fr_soilslow;

		// Updated soil fluxes. In wetlands, some portion of respiration can be emitted as CH4, so save respiration as respsum until Soil::methane() is called.  
		if (patch.stand.landcover != PEATLAND) {
			patch.fluxes.report_flux(Fluxes::SOILC, cflux);
			soillayer.respsum = 0.0;
		}
		else {
			soillayer.respsum = cflux;
		}
	}

	// Solve SOM pool sizes at end of year given by soilcolumn.solvesom_end

	if (date.year == soilcolumn.soiltype.solvesom_end && date.islastmonth && date.islastday)
		equilsom_lpj(soilcolumn);

}

/////////////////////////////////////////////////
// CENTURY SOM DYNAMICS

/// Reduce decay rates and microbial carbon use efficiency to keep the daily nitrogen and phosphorus balance in the soil
/**
*  Reduce decay rates and microbial carbon use efficiency to keep the daily nitrogen and phosphorus balance in the soil
*/
void reduce_decay_growth_rates(double decay_reduction[NSOMPOOL], double growth_reduction[NSOMPOOL],
	double cdec[NSOMPOOL], double ndec[NSOMPOOL], double ninc[NSOMPOOL],
	double respsum[NSOMPOOL], double tot_net_n_min) {

	// Keeping track on which pools that are negative
	bool neg_pools[NSOMPOOL] = { 0 };

	// Total negative N from pools with N total immobilization
	double tot_imm_n = 0.0;

	// Add up immobilisation for considered pools
	for (int p = 0; p < NSOMPOOL; p++) {
		if (ndec[p] - ninc[p] < 0.0) {
			tot_imm_n += ndec[p] - ninc[p];
			neg_pools[p] = true;
		}
	}

	// Calculate difference between total immobilsation and total negative net mineralisation
	double diff_imm_to_netmin = tot_net_n_min / tot_imm_n;

	// With lower immobilisation difference to net negative mineralisation (almost no mineralisation - imm == netmin) 
	// decay rates is reduced more as a fraction of the total solution to the nutrient mass balance missmatch. 
	// At higher difference most of the solution is through changes in CUE

	// Decay reduction can not go over 50% of the solution to the mass balance missmatch
	//double decay_rate_reduction = 1.0 - 0.5 * diff_imm_to_netmin * diff_imm_to_netmin;
	double decay_rate_reduction = 1.0 - 0.5 * diff_imm_to_netmin;

	double netmin_before = 0.0;
	double netmin_after = 0.0;

	// Reduce decay rate for immobilising pools
	for (int p = 0; p < NSOMPOOL; p++) {
		if (neg_pools[p]) {
			decay_reduction[p] = decay_rate_reduction;

			// Net min before
			netmin_before += ndec[p] - ninc[p];

			// Adjust corresponding C, N and P movements
			cdec[p] *= decay_rate_reduction;
			ndec[p] *= decay_rate_reduction;
			respsum[p] *= decay_rate_reduction;
			ninc[p] *= decay_rate_reduction;

			// Net min before
			netmin_after += ndec[p] - ninc[p];
		}
	}

	// Add reduction in immobilisation to total net mineralisation
	tot_net_n_min += netmin_after - netmin_before;

	// Update sum of total C and N movements
	double tot_cdec, tot_ninc, tot_resp;
	tot_cdec = tot_ninc = tot_resp = 0.0;

	for (int p = 0; p < NSOMPOOL; p++) {
		if (neg_pools[p]) {
			tot_cdec += cdec[p];
			tot_ninc += ninc[p];
			tot_resp += respsum[p];
		}
	}

	// Calculate growth fraction decrease to satisfy N deficit

	// Calculate total growth fraction after reduced decay rates
	double tot_cinc = tot_cdec - tot_resp;
	double tot_growthfrac = tot_cinc / tot_cdec;

	// Calculate needed growth fraction to satisfy a neutral net mineralisation
	// and the C:N ratio of receiver pools

	// Receiver total C:N ratio
	double cton = tot_cinc / tot_ninc;

	// C increment allowed to meet a neutral net mineralisation rate
	double tot_cinc_new = (tot_ninc + tot_net_n_min) * cton;

	// New total growth fraction
	double tot_growthfrac_new = tot_cinc_new / tot_cdec;

	// Growth fraction decrease
	double gr = tot_growthfrac_new / tot_growthfrac;

	for (int p = 0; p < NSOMPOOL; p++) {
		if (neg_pools[p]) {
			growth_reduction[p] = gr;
		}
	}
}

/// Set target N:C ratios for SOM pools
/** Set target N:C ratios for slow, passive, humus and soil microbial pools
 *  based on mineral nitrogen pool or litter nitrogen fraction (Parton et al 1993, Fig 4)
 */
void Soillayer::setntoc(double fac, pooltype pool, double cton_max, double cton_min,
	double fmin, double fmax) {

	if (fac <= fmin)
		sompool[pool].ntoc = 1.0 / cton_max;
	else if (fac >= fmax)
		sompool[pool].ntoc = 1.0 / cton_min;
	else {
		sompool[pool].ntoc = 1.0 / (cton_min + (cton_max - cton_min) *
			(fmax - fac) / (fmax - fmin));
	}
}


/// Temperature modifier for decomposition
/** Calculate decomposition temperature modifier (in range 0-1)
 *  [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
 *  effect of soil temperature, Parton et al 1993, Fig 2)
 *
 *  \param temp_soil  Soil temperature at 25 cm depth
 */
double temperature_modifier(double temp_soil) {

	double temp_mod = temp_soil > 0.0 ? max(0.0, 0.0326 + 0.00351 * pow(temp_soil, 1.652) - pow(temp_soil / 41.748, 7.19)) : 0.0;

	// Include as an overide option when: MIN_DECOMP_TEMP < temp < 0 degC
	// This increases the respiration (from 0) between MIN_DECOMP_TEMP and 0 degC - cf Koven et al. 2011
	if (ifcarbonfreeze && temp_soil <= 0.0 && temp_soil >= MIN_DECOMP_TEMP && !iftwolayersoil) {

		double decomp_at_freezing_point = 0.0326; // temp_mod above when temp_soil = 0;
		bool linear_decrease_below_freezing = false;

		if (linear_decrease_below_freezing) {
			// Alternative: Linear approach (Koven et al. 2011)
			double slope = decomp_at_freezing_point / fabs(MIN_DECOMP_TEMP);
			temp_mod = slope * temp_soil + decomp_at_freezing_point; // i.e. a linear decrease from decomp_at_freezing_point at 0C to 0 at MIN_DECOMP_TEMP (-4C).
		} 
		else {
			// Default: Q10 relationship (Schaefer & Jafarov, 2016)
			double q10_freeze = 2.9; // (Natali et al., 2019)
			temp_mod = decomp_at_freezing_point * pow(q10_freeze, temp_soil / 10.0);
		}
	}

	return temp_mod;
}


/// Water Modifier for Decomposition
/** Calculate decomposition moisture modifier (in range 0-1)
 *  Friend et al 1997, Eqn 53 (Parton et al 1993, Fig 2)
 *
 *  \param wfps  Water-filled pore space
 *  \param theta  Total water mass (ice + water)
 */
double Soillayer::moisture_modifier() {

	// double moist_mod = wfps < 60.0 ? exp((wfps - 60.0) * (wfps - 60.0) / -800.0) : 0.000371 * wfps * wfps - 0.0748 * wfps + 4.13;

	// Yan et al 2018 - A moisture function of soil heterotrophic
	// respiration that incorporates microscale processes

	// Minimum moist modifier for dry soils
	double moist_mod_dry = 0.064;

	// Maximum moist modifier at optimal moisture
	double moist_mod_max = (*soilcolumn).patch.stand.landcover == PEATLAND ? RMOIST : 1.0;

	// Total water mass
	double theta = wmass + imass;


	double moist_mod;

	if (theta < fc()) {
		moist_mod = max(moist_mod_dry, moist_mod_max * (fc() + 0.1) / (theta + 0.1) * pow(theta / fc(), 1.0 + a * 2.0));
	}
	else if (wsats - theta > 1.0e-12) {
		moist_mod = max(RMOIST_ANAEROBIC, moist_mod_max * pow((wsats - theta) / (wsats - fc()), 0.75));
	}
	else {
		moist_mod = RMOIST_ANAEROBIC;
	}
	
	return moist_mod;
}

/// Calculates CENTURY instantaneous decay rates
/** Calculates CENTURY instantaneous decay rates given soil layer temperature,
*  water content, texture etc
*/
void Soillayer::decayrates_century(bool tillage, bool savesomcent) {

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris decay
	// pools SURFMETA,SURFSTRUCT,SURFFWD,SURFCWD,SURFMICRO,SURFHUMUS,SOILMETA,SOILSTRUCT,SOILFWD,SOILCWD,SOILMICRO,SLOWSOM,PASSIVESOM
	//const double K_MAX[] = { 3.8e-2, 9.5e-3, 1.1e-2, 2.2e-3, 2.7e-2, 4.8e-4, 7.0e-2, 1.9e-2, 1.1e-2, 2.2e-3, 4.2e-2, 1.7e-3, 1.9e-6 };
	const double K_MAX[] = { 8.0e-2, 2.2e-2, 1.1e-2, 2.2e-3, 4.0e-2, 2.7e-4, 1.0e-1, 2.7e-2, 1.4e-2, 2.7e-3, 5.0e-2, 2.7e-4, 4.2e-6 };

	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	const double texture_mod = 1.0 - 0.75 * ((*soilcolumn).soiltype.clay_frac + (*soilcolumn).soiltype.silt_frac);

	// Calculate decomposition temperature modifier (in range 0-1)
	double temp_mod = temperature_modifier(temp);

	// Calculate decomposition moisture modifier
	double moist_mod = moisture_modifier();	

	// Combined moisture and temperature modifier

	for (int p = 0; p < NSOMPOOL; p++) {

		// Calculate decay constant
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		double scalar = 1.0; // DW_COMMENT
		if (p == SURFHUMUS || p == SLOWSOM || p == PASSIVESOM)
			scalar = 1.0;

		double k = scalar * K_MAX[p] * temp_mod * moist_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2
		// Kirschbaum et al 2001 changed the exponential term
		// from 3 to 5.

		if (p == SURFSTRUCT || p == SOILSTRUCT || p == SURFFWD || p == SOILFWD || p == SURFCWD || p == SOILCWD) {
			k *= exp(-5.0 * sompool[p].ligcfrac);
		}
		else if (p == SOILMICRO && !acrotelm && !catotelm) {
			k *= texture_mod;
		}

		// Decay rate modifier depending on depth only after non wetland soils
		if ((*soilcolumn).patch.stand.landcover != PEATLAND) {
			k *= depth_mod;
		}

		// Increased decay rates for crops (tillage)
		if (tillage && (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM) && depth_profile < TILLAGE_DEPTH) {
			k *= TILLAGE_FACTOR;
		}

		// Calculate fraction of carbon pool remaining after today's decomposition
		sompool[p].fracremain = exp(-k);

		if (savesomcent) {
			sompool[p].mfracremain_mean[date.month] += sompool[p].fracremain / date.ndaymonth[date.month];
		}
	}
}

/// Transfers specified fraction (frac) of today's decomposition
/** Transfers specified fraction (frac) of today's decomposition in donor pool type
  *  to receiver pool, transferring fraction respfrac of this to the accumulated CO2
 *  flux respsum (representing total microbial respiration today)
 */
void transferdecomp(Soillayer& soillayer, pooltype donor, pooltype receiver,
	double frac, double respfrac, double growth_reduction, double& nmin_actual,
	double& nimmob) {

	respfrac = 1.0 - (1.0 - respfrac) * growth_reduction;

	// decrement in donor carbon and nitrogen pools
	double cdec = soillayer.sompool[donor].cdec * frac;
	double ndec = soillayer.sompool[donor].ndec * frac;

	// associated nitrogen increment in receiver pool (Friend et al 1997, Eqn 49)
	double cinc = cdec * (1.0 - respfrac);
	double ninc = cinc * soillayer.sompool[receiver].ntoc;

	soillayer.sompool[donor].cinc += cinc;
	soillayer.sompool[donor].ninc += ninc;

	// if increase in receiver nitrogen greater than decrease in donor nitrogen,
	// balance must be immobilisation from mineral nitrogen pool
	// otherwise balance is nitrogen mineralisation
	if (ninc > ndec) {
		nimmob += ninc - ndec;
	}
	else {
		nmin_actual += ndec - ninc;
	}

	// "Transfer" carbon and nitrogen to receiver
	soillayer.sompool[receiver].delta_cmass += cdec * (1.0 - respfrac);
	soillayer.sompool[receiver].delta_nmass += ninc;

	// Transfer microbial respiration
	soillayer.sompool[donor].respsum += cdec * respfrac;
}

/// Fluxes between the CENTURY pools, and CO2 release to the atmosphere
/** Daily or monthly fluxes between the ten CENTURY pools, and CO2 release to the atmosphere
 *  Parton et al 1993, Fig 1; Comins & McMurtrie 1993, Appendix A
 *
 *  \param ifequilsom Whether the function is called during calculation om SOM pool equilibrium,
 *                    \see equilsom(). During this stage, somfluxes shouldn't calculate decayrates
 *                    itself, and shouldn't produce output like fluxes etc.
 */
void somfluxes(Patch& patch, bool ifequilsom, bool tillage) {

	Soilcolumn& soilcolumn = patch.soilcolumn;
	Soiltype soiltype = soilcolumn.soiltype;

	double leachsum_cmass = 0.0;
	double leachsum_nmass = 0.0;
	// Actual (not net) nitrogen mineralisation
	double nmin_actual = 0.0;
	// Nitrogen immobilisation
	double nimmob = 0.0;		

	const double EPS = 1.0e-16;

	// Check mass balances
	// double ncont_before = soilcolumn.get_sum(NCONT);
	// double ccont_before = soilcolumn.get_sum(CCONT);
	// double nextra = 0.0;
	// double cextra = 0.0;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Mass balance check
		// double cmass_before = soillayer.ccont();
		// double nmass_before = soillayer.ncont();

		// mineral nitrogen mass available
		const double nmin_mass = soillayer.nmass_avail(NH4);
		if (date.day == 0) {
			soillayer.anmin = 0.0;
			soillayer.animmob = 0.0;
			soillayer.acinc[0] = 0.0;
			soillayer.acinc[1] = 0.0;
			soillayer.acdec[0] = 0.0;
			soillayer.acdec[1] = 0.0;
		}

		// Warning if soil available nitrogen is negative (if happens once or so no problem, but if it propagates through time then it is)
		if (ifnlim) {
			assert(soillayer.NH4_mass > -EPS);
		}

		// Set N:C ratios for humus, soil microbial, slow and passive pool based on estimated mineral nitrogen pool
		// (Parton et al 1993, Fig 4)
		// ForCent (Parton 2010) values
		soillayer.setntoc(nmin_mass, SURFHUMUS, 30.0, 10.0, 0.0, NMASS_SAT);
		soillayer.setntoc(nmin_mass, SOILMICRO, 15.0, 6.0, 0.0, NMASS_SAT);
		soillayer.setntoc(nmin_mass, SLOWSOM, 30.0, 10.0, 0.0, NMASS_SAT);

		if (!ifequilsom) {

			// Calculate potential fraction remaining following decay today for all pools
			// (assumes no nitrogen limitation)
			soillayer.decayrates_century(tillage, soilcolumn.savesomcent);
		}

		// Calculate decomposition in all pools assuming these decay rates

		// Save delta carbon and nitrogen mass

		bool net_mineralization = false;
		int times = 0;
		double decay_reduction[NSOMPOOL];
		double growth_reduction[NSOMPOOL];
		for (int p = 0; p < NSOMPOOL; p++) {
			decay_reduction[p] = 1.0;
			growth_reduction[p] = 1.0;

			if (soillayer.sompool[p].cmass < 0.0) {
				dprintf("Year %d Day %d Lyr %d Pool %d CMASS %g NMASS %g\n", date.year, date.day, lyr, p, soillayer.sompool[p].cmass, soillayer.sompool[p].nmass);
				soillayer.sompool[p].cmass = 0.0;
				soillayer.sompool[p].nmass = 0.0;
			}
		}
		double init_negative_nmass;
		double init_ntoc_reduction = 1.0;
		double ntoc_reduction = 0.8;
		double cdec_save[2] = { 0.0 }, cinc_save[2] = { 0.0 };

		// If mineralization together with soil available nitrogen is negative then decay rates are decreased
		// The SOM system have five try to get a positive result, after that all pools decay rate has been
		// affected by nitrogen limitation
		while (!net_mineralization && times < 2) {

			leachsum_cmass = leachsum_nmass = 0.0;
			nmin_actual = nimmob = 0.0;

			// Calculate decomposition in all pools assuming these decay rates
			for (int p = 0; p < NSOMPOOL; p++) {
				soillayer.sompool[p].cdec = soillayer.sompool[p].cmass * (1.0 - soillayer.sompool[p].fracremain) * decay_reduction[p];
				soillayer.sompool[p].ndec = soillayer.sompool[p].nmass * (1.0 - soillayer.sompool[p].fracremain) * decay_reduction[p];

				soillayer.sompool[p].delta_cmass = -soillayer.sompool[p].cdec;
				soillayer.sompool[p].delta_nmass = -soillayer.sompool[p].ndec;

				soillayer.sompool[p].cinc = 0.0;
				soillayer.sompool[p].ninc = 0.0;

				soillayer.sompool[p].respsum = 0.0;
			}

			// Partition potential decomposition among receiver pools

			// Donor pool SURFACE STRUCTURAL

			transferdecomp(soillayer, SURFSTRUCT, SURFMICRO, 1.0 - soillayer.sompool[SURFSTRUCT].ligcfrac,
				0.6, growth_reduction[SURFSTRUCT], nmin_actual, nimmob);

			transferdecomp(soillayer, SURFSTRUCT, SURFHUMUS, soillayer.sompool[SURFSTRUCT].ligcfrac, 0.3,
				growth_reduction[SURFSTRUCT], nmin_actual, nimmob);

			// Donor pool SURFACE METABOLIC

			transferdecomp(soillayer, SURFMETA, SURFMICRO, 1.0, 0.6, growth_reduction[SURFMETA], nmin_actual, nimmob);

			// Donor pool SOIL STRUCTURAL

			transferdecomp(soillayer, SOILSTRUCT, SOILMICRO, 1.0 - soillayer.sompool[SOILSTRUCT].ligcfrac,
				0.55, growth_reduction[SOILSTRUCT], nmin_actual, nimmob);

			transferdecomp(soillayer, SOILSTRUCT, SLOWSOM, soillayer.sompool[SOILSTRUCT].ligcfrac, 0.3,
				growth_reduction[SOILSTRUCT], nmin_actual, nimmob);

			// Donor pool SOIL METABOLIC

			transferdecomp(soillayer, SOILMETA, SOILMICRO, 1.0, 0.55, growth_reduction[SOILMETA], nmin_actual, nimmob);

			// Donor pool SURFACE FINE WOODY DEBRIS

			transferdecomp(soillayer, SURFFWD, SURFMICRO, 1.0 - soillayer.sompool[SURFFWD].ligcfrac,
				0.76, growth_reduction[SURFFWD], nmin_actual, nimmob);

			transferdecomp(soillayer, SURFFWD, SURFHUMUS, soillayer.sompool[SURFFWD].ligcfrac, 0.4,
				growth_reduction[SURFFWD], nmin_actual, nimmob);

			// Donor pool SURFACE COARSE WOODY DEBRIS

			transferdecomp(soillayer, SURFCWD, SURFMICRO, 1.0 - soillayer.sompool[SURFCWD].ligcfrac,
				0.9, growth_reduction[SURFCWD], nmin_actual, nimmob);

			transferdecomp(soillayer, SURFCWD, SURFHUMUS, soillayer.sompool[SURFCWD].ligcfrac, 0.5,
				growth_reduction[SURFCWD], nmin_actual, nimmob);

			// Donor pool SURFACE MICROBE

			transferdecomp(soillayer, SURFMICRO, SURFHUMUS, 1.0, 0.6, growth_reduction[SURFMICRO], nmin_actual, nimmob);

			// Donor pool SURFACE HUMUS

			transferdecomp(soillayer, SURFHUMUS, SLOWSOM, 1.0, 0.6, growth_reduction[SURFHUMUS], nmin_actual, nimmob);

			// Donor pool SOIL FINE WOODY DEBRIS

			transferdecomp(soillayer, SOILFWD, SOILMICRO, 1.0 - soillayer.sompool[SOILFWD].ligcfrac,
				0.76, growth_reduction[SOILFWD], nmin_actual, nimmob);

			transferdecomp(soillayer, SOILFWD, SLOWSOM, soillayer.sompool[SOILFWD].ligcfrac, 0.4,
				growth_reduction[SOILFWD], nmin_actual, nimmob);

			// Donor pool SOIL COARSE WOODY DEBRIS

			transferdecomp(soillayer, SOILCWD, SOILMICRO, 1.0 - soillayer.sompool[SOILCWD].ligcfrac,
				0.9, growth_reduction[SOILCWD], nmin_actual, nimmob);

			transferdecomp(soillayer, SOILCWD, SLOWSOM, soillayer.sompool[SOILCWD].ligcfrac, 0.5,
				growth_reduction[SOILCWD], nmin_actual, nimmob);

			// Donor pool SLOW SOM

			// Calculate partitioning coefficients (Fig 1, Parton et al 1993)
			double respfrac = 0.55;
			double csp = max(0.0, 0.003 - 0.009 * soiltype.clay_frac);
			double csa = 1.0 - csp - respfrac;
			double frac_csp = csp / (csp + csa);
			double frac_csa = 1.0 - frac_csp;

			transferdecomp(soillayer, SLOWSOM, SOILMICRO, frac_csa, respfrac, growth_reduction[SLOWSOM], nmin_actual, nimmob);

			transferdecomp(soillayer, SLOWSOM, PASSIVESOM, frac_csp, respfrac, growth_reduction[SLOWSOM], nmin_actual, nimmob);

			// Donor pool SOIL MICROBE

			// Calculate partitioning coefficients (Fig 1, Parton et al 1993)
			// Fraction lost to  microbial respiration (F_t, Parton et al 1993 Eqn 7)
			respfrac = max(0.0, 0.85 - 0.68 * (soiltype.clay_frac + soiltype.silt_frac));

			// Fraction entering passive SOM pool (Parton et al 1993, Eqn 9)
			double cap = 0.003 + 0.032 * soiltype.clay_frac;

			// Fraction entering slow SOM pool
			double cas = 1.0 - respfrac - soillayer.orgleachfrac - cap;

			double frac_cap = cap / (cap + cas) * (1.0 - soillayer.orgleachfrac);
			double frac_cas = cas / (cap + cas) * (1.0 - soillayer.orgleachfrac);
			respfrac /= (1.0 - soillayer.orgleachfrac);


			transferdecomp(soillayer, SOILMICRO, PASSIVESOM, frac_cap, respfrac, growth_reduction[SOILMICRO], nmin_actual, nimmob);

			transferdecomp(soillayer, SOILMICRO, SLOWSOM, frac_cas, respfrac, growth_reduction[SOILMICRO], nmin_actual, nimmob);

			// Account for organic carbon and nitrogen leaching loss
			leachsum_cmass = soillayer.orgleachfrac * soillayer.sompool[SOILMICRO].cdec;

			if (!negligible(soillayer.sompool[SOILMICRO].cmass)) {
				leachsum_nmass = soillayer.orgleachfrac * soillayer.sompool[SOILMICRO].cdec * soillayer.sompool[SOILMICRO].nmass / soillayer.sompool[SOILMICRO].cmass;
			}

			// Donor pool PASSIVE SOM

			transferdecomp(soillayer, PASSIVESOM, SOILMICRO, 1.0, 0.55, growth_reduction[PASSIVESOM], nmin_actual, nimmob);

			// Total net mineralization
			double tot_net_min = nmin_actual - nimmob;

			// Estimate daily soil mineral nitrogen pool after decomposition
			// (negative value = immobilisation)
			if (tot_net_min + EPS >= 0.0) {

				net_mineralization = true;
			}
			else if (ifnlim && date.year > freenyears) {

				// Net N and/or P mineralisation is negative -> reduce decay rates and microbial CUE
				if (times < 1) {

					double cdec[NSOMPOOL], ndec[NSOMPOOL], ninc[NSOMPOOL], respsum[NSOMPOOL];
					for (int p = 0; p < NSOMPOOL; p++) {
						cdec[p] = soillayer.sompool[p].cdec;
						ndec[p] = soillayer.sompool[p].ndec;
						ninc[p] = soillayer.sompool[p].ninc;
						respsum[p] = soillayer.sompool[p].respsum;
					}
					reduce_decay_growth_rates(decay_reduction, growth_reduction, cdec, ndec, ninc, respsum, tot_net_min);
				}
				else {
					if (!ifequilsom)
						dprintf("Year %d Day %d FAIL: Lyr %d Times %d Min %g Imm %g NetMin %g\n", date.year, date.day, lyr, times, nmin_actual, nimmob, tot_net_min);
				}
				net_mineralization = false;
			}
			else {

				// Not minding immobilisation higher than nmass_avail during free nitrogen years
				if (date.year > freenyears) {

					// Immobilization larger than soil available nitrogen -> reduce targeted N concentration in SOM pool with flexible N:C ratios
					if (times == 0) {
						// initial reduction
						init_negative_nmass = tot_net_min + nmin_mass;
						init_ntoc_reduction = ntoc_reduction;
					}
					else {
						// trying to match needed N:C reduction
						ntoc_reduction = min(init_ntoc_reduction, pow(init_ntoc_reduction, 1.0 / (1.0 - (tot_net_min + nmin_mass) / init_negative_nmass) + 1.0));
					}

					soillayer.sompool[SLOWSOM].ntoc *= ntoc_reduction;
					soillayer.sompool[SOILMICRO].ntoc *= ntoc_reduction;
					soillayer.sompool[SURFHUMUS].ntoc *= ntoc_reduction;

					net_mineralization = false;
				}
				else {
					net_mineralization = true;
				}
			}

			// Save initial and final decrease and increase in C pools
			// Used for annual output files of nutrient limitation on 
			// decomposition and growth
			if (!ifequilsom) {
				if (times > 0) {
					cdec_save[1] = 0.0;
					cinc_save[1] = 0.0;
				}

				for (int p = 0; p < NSOMPOOL; p++) {
					if (times == 0) {
						cdec_save[0] += soillayer.sompool[p].cdec;
						cinc_save[0] += soillayer.sompool[p].cinc;
					}
					cdec_save[1] += soillayer.sompool[p].cdec;
					cinc_save[1] += soillayer.sompool[p].cinc;
				}
			}

			times++;
		}

		// Sum annual C increase and decrease to be able to
		// calculate annual nutient limitation on decomposition
		soillayer.acdec[0] += cdec_save[0];
		soillayer.acdec[1] += cdec_save[1];
		soillayer.acinc[0] += cinc_save[0];
		soillayer.acinc[1] += cinc_save[1];

		double respsum_tot = 0.0;

		// Update pool sizes

		for (int p = 0; p < NSOMPOOL; p++) {
			soillayer.sompool[p].cmass += soillayer.sompool[p].delta_cmass;
			soillayer.sompool[p].nmass += soillayer.sompool[p].delta_nmass;

			respsum_tot += soillayer.sompool[p].respsum;
		}

		if (!ifequilsom) {

			// Updated soil fluxes. In wetlands, some portion of respiration can be emitted as CH4, so save respiration as respsum until Soil::methane() is called.  
			if (patch.stand.landcover != PEATLAND) {
				patch.fluxes.report_flux(Fluxes::SOILC, respsum_tot);
				soillayer.respsum = 0.0;
			}
			else {
				soillayer.respsum = respsum_tot;
			}

			// Mass balance check
			// cextra += respsum_tot;

			soillayer.arespsum += respsum_tot;

			// Sum annual organic nitrogen leaching
			soillayer.aorgNleach += leachsum_nmass;

			// Mass balance check
			// nextra += leachsum_nmass;

			// Sum annual organic carbon leaching
			soillayer.aorgCleach += leachsum_cmass;

			// Mass balance check
			// cextra += leachsum_cmass;

			// Sum annuals
			soillayer.anmin += nmin_actual;
			soillayer.animmob += nimmob;
		}

		// Fraction of microbial resp. is assumed to produce labile carbon
		soillayer.labile_carbon = respsum_tot * frac_labile_carbon;

		// Adding mineral nitrogen to soil available pool
		double nmin_inc = nmin_actual - nimmob;

		// If no nitrogen limitation or during free nitrogen years
		// net mineralisation could be negative. Set to zero
		// in that case
		if (!ifnlim || date.year <= freenyears)
			nmin_inc = max(0.0, nmin_inc);

		if (!ifequilsom) {
			soillayer.mnetmin[date.month] += nmin_inc;
		}
		double nflux_tmp = 0.0;
		// Estimate of N flux from soil (simple CLM-CN approach)
		if (!ifntransform) {
			double nflux = nmin_inc > 0.0 ? (nmin_inc) * 0.01 : 0.0;
			nmin_inc -= nflux;
			// Mass balance check
			// nextra += nflux;
			nflux_tmp += nflux;

			if (!ifequilsom) {
				patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);
			}
		}
		soillayer.nmass_inc(nmin_inc, NH4);

		// If no nitrogen limitation or during free nitrogen years set soil
		// available nitrogen to its saturation level.
		if (!ifnlim || date.year <= freenyears) {
			if (ifntransform) {
				soillayer.NH4_mass = NMASS_SAT / 2.0;
				soillayer.NO3_mass = NMASS_SAT / 2.0;
			}
			else {
				soillayer.NH4_mass = NMASS_SAT;
			}
		}

		// Mass balance check
		// double cmass_after = soillayer.ccont() + respsum_tot + leachsum_cmass;
		// double nmass_after = soillayer.ncont() + nflux_tmp + leachsum_nmass;

		// if (date.year > freenyears && fabs(cmass_before - cmass_after) > 1.0e-12 && !ifequilsom)
		//	dprintf("Year %d Day %d Lyr %d C SOM diff %g %g %g Leach %g\n", date.year, date.day, lyr, fabs(cmass_before - cmass_after), cmass_before, cmass_after, leachsum_cmass);

		// if (date.year > freenyears && fabs(nmass_before - nmass_after) > 1.0e-12 && !ifequilsom)
		//	dprintf("Year %d Day %d Lyr %d N SOM diff %g %g %g Leach %g\n", date.year, date.day, lyr, fabs(nmass_before - nmass_after), nmass_before, nmass_after, leachsum_nmass);

	}

	// Mass balance check
	// double ncont_after = soilcolumn.get_sum(NCONT) + nextra;
	// double ccont_after = soilcolumn.get_sum(CCONT) + cextra;

	//if (date.year > freenyears && fabs(ncont_before - ncont_after) > 1.0e-12 && !ifequilsom)
	//	dprintf("Year %d Day %d NITROGEN SOM diff %g %g %g\n", date.year, date.day, fabs(ncont_before - ncont_after), ncont_before, ncont_after);

	//if (date.year > freenyears && fabs(ccont_before - ccont_after) > 1.0e-12 && !ifequilsom)
	//	dprintf("Year %d Day %d CARBON SOM diff %g %g %g\n", date.year, date.day, fabs(ccont_before - ccont_after), ccont_before, ccont_after);
}

/// Litter lignin to N ratio (for leaf and root litter)
/** Specific lignin fractions for leaf and root
 *  are specified in transfer_litter()
 */
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr) {

	if (!negligible(nmass_litter) && ifnlim) {
		return max(0.0, LIGCFRAC * cmass_litter / nmass_litter);
	}
	else {
		return max(0.0, LIGCFRAC * cton_avr / (1.0 - nrelocfrac));
	}
}

/// Metabolic litter fraction (for leaf and root litter)
/** Fm, Parton et al 1993, Eqn 1:
 *  NB: incorrect/out-of-date intercept and slope given in Eqn 1; values used in
 *  code of CENTURY 4.0 used instead (also correct in Parton et al. 1993, figure 1)
 *
 * \param lton  Litter lignin:N ratio
 */
double metabolic_litter_fraction(double lton) {
	return max(0.0, 0.85 - lton * 0.013);
}

/// Transfers litter from growth, mortality and fire
/** Called monthly to transfer last year's litter from vegetation
 *  (turnover, mortality and fire) to soil litter pools.
 *  Alternatively, with daily allocation and harvest/turnover,
 *  the litter produced a certain day.
 */
void transfer_litter(Patch& patch) {

	// For CROPLAND only enter on harvest/turnover day otherwise every day
	if (patch.stand.landcover != CROPLAND || patch.is_litter_day) {

		// Soil column of all soil objects
		Soilcolumn& soilcolumn = patch.soilcolumn;

		double lat = patch.get_climate().lat;

		double EPS = -1.0e-16;

		double ligcmass_new, ligcmass_old;

		// Fire (GlobFIRM)
		double litterme[NSOMPOOL];
		double fireresist[NSOMPOOL];

		patch.pft.firstobj();
		while (patch.pft.isobj) {
			Patchpft& ppft = patch.pft.getobj();

			// For stands with yearly growth, drop leaf and root litter on first month of the year for northern hemisphere
			// and first month of the second half of the year for southern hemisphere for summergreen trees. For evergreens
			// as a fraction every day and for raingreens on the drierst month from last year. 
			// For stands with daily growth, harvest and/or turnover, do this when patch.is_litter_day is true.

			double cmass_litter_leaf, nmass_litter_leaf;
			double cmass_litter_root, nmass_litter_root;
			double cmass_litter_sap, nmass_litter_sap;
			double cmass_litter_heart, nmass_litter_heart;

			//  Is a litter day, drop all leaf litter (crop)
			if (patch.is_litter_day) {
				cmass_litter_leaf = ppft.litter_leaf;
				nmass_litter_leaf = ppft.nmass_litter_leaf;
				cmass_litter_root = ppft.litter_root;
				nmass_litter_root = ppft.nmass_litter_root;

				ppft.litter_leaf = 0.0;
				ppft.nmass_litter_leaf = 0.0;
				ppft.litter_root = 0.0;
				ppft.nmass_litter_root = 0.0;
			}
			// For summergreens drop leaf litter over all days during Jan in NH and July SH
			else if ((ppft.pft.phenology == SUMMERGREEN && ((lat >= 0.0 && date.month == 0) || (lat < 0.0 && date.month == 6))) ||
				(ppft.pft.phenology == RAINGREEN && date.month == ppft.driest_mth)) {
				double frac = 1.0 / (date.ndaymonth[date.month] - date.dayofmonth);
				cmass_litter_leaf = ppft.litter_leaf * frac;
				nmass_litter_leaf = ppft.nmass_litter_leaf * frac;
				cmass_litter_root = ppft.litter_root * frac;
				nmass_litter_root = ppft.nmass_litter_root * frac;

				ppft.litter_leaf -= cmass_litter_leaf;
				ppft.nmass_litter_leaf -= nmass_litter_leaf;
				ppft.litter_root -= cmass_litter_root;
				ppft.nmass_litter_root -= nmass_litter_root;
			}
			// Drops leaf litter every day
			else {
				double frac = 1.0 / (date.year_length() - date.day);
				cmass_litter_leaf = ppft.litter_leaf * frac;
				nmass_litter_leaf = ppft.nmass_litter_leaf * frac;
				cmass_litter_root = ppft.litter_root * frac;
				nmass_litter_root = ppft.nmass_litter_root * frac;

				ppft.litter_leaf -= cmass_litter_leaf;
				ppft.nmass_litter_leaf -= nmass_litter_leaf;
				ppft.litter_root -= cmass_litter_root;
				ppft.nmass_litter_root -= nmass_litter_root;
			}

			if (ppft.pft.lifeform == TREE) {

				// Woody litter is dropped as a fraction every day
				double frac = 1.0 / (date.year_length() - date.day);
				cmass_litter_sap = ppft.litter_sap * frac;
				nmass_litter_sap = ppft.nmass_litter_sap * frac;
				cmass_litter_heart = ppft.litter_heart * frac;
				nmass_litter_heart = ppft.nmass_litter_heart * frac;

				ppft.litter_sap -= cmass_litter_sap;
				ppft.nmass_litter_sap -= nmass_litter_sap;
				ppft.litter_heart -= cmass_litter_heart;
				ppft.nmass_litter_heart -= nmass_litter_heart;
			}
			else {
				cmass_litter_sap = 0.0;
				nmass_litter_sap = 0.0;
				cmass_litter_heart = 0.0;
				nmass_litter_heart = 0.0;
			}

			if (!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf) || !negligible(cmass_litter_root) || !negligible(nmass_litter_root) || !negligible(cmass_litter_sap) || !negligible(nmass_litter_sap) || !negligible(cmass_litter_heart) || !negligible(nmass_litter_heart)) {

				for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

					// Soil object for this layer
					Soillayer& soillayer = soilcolumn.soillayers[lyr];

					double root_lyr_frac = ppft.rootdist[lyr];

					// For peatland all litter enter the top layer
					if (soilcolumn.patch.stand.landcover == PEATLAND) {
						root_lyr_frac = !lyr ? 1.0 : 0.0;
					}

					if (firemodel == GLOBFIRM && lyr == 0) {
						litterme[SURFSTRUCT] = soillayer.sompool[SURFSTRUCT].cmass * soillayer.sompool[SURFSTRUCT].litterme;
						litterme[SURFMETA] = soillayer.sompool[SURFMETA].cmass * soillayer.sompool[SURFMETA].litterme;
						litterme[SURFFWD] = soillayer.sompool[SURFFWD].cmass * soillayer.sompool[SURFFWD].litterme;
						litterme[SURFCWD] = soillayer.sompool[SURFCWD].cmass * soillayer.sompool[SURFCWD].litterme;

						fireresist[SURFSTRUCT] = soillayer.sompool[SURFSTRUCT].cmass * soillayer.sompool[SURFSTRUCT].fireresist;
						fireresist[SURFMETA] = soillayer.sompool[SURFMETA].cmass * soillayer.sompool[SURFMETA].fireresist;
						fireresist[SURFFWD] = soillayer.sompool[SURFFWD].cmass * soillayer.sompool[SURFFWD].fireresist;
						fireresist[SURFCWD] = soillayer.sompool[SURFCWD].cmass * soillayer.sompool[SURFCWD].fireresist;
					}

					// Transferred leaf litter to surface soil layer
					if ((!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf)) && lyr == 0) {

						// LEAF

						// Calculate inputs to surface structural and metabolic litter

						// Leaf litter lignin:N ratio
						double leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, ppft.pft.ligcfrac_leaf, ppft.pft.cton_leaf_avr);

						// Metabolic litter fraction for leaf litter
						double fm = metabolic_litter_fraction(leaf_lton);

						ligcmass_new = cmass_litter_leaf * (1.0 - fm) * ppft.pft.ligcfrac_leaf;
						ligcmass_old = soillayer.sompool[SURFSTRUCT].cmass * soillayer.sompool[SURFSTRUCT].ligcfrac;

						// Add to pools
						soillayer.sompool[SURFSTRUCT].cmass += cmass_litter_leaf * (1.0 - fm);
						soillayer.sompool[SURFSTRUCT].nmass += nmass_litter_leaf * (1.0 - fm);
						soillayer.sompool[SURFMETA].cmass += cmass_litter_leaf * fm;
						soillayer.sompool[SURFMETA].nmass += nmass_litter_leaf * fm;

						soilcolumn.cmass_litter_input[0] += cmass_litter_leaf;
						soilcolumn.csom_litter_input[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm);
						soilcolumn.csom_litter_input[SURFMETA] += cmass_litter_leaf * fm;

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_leaf * (1.0 - fm), nmass_litter_leaf * (1.0 - fm), SURFSTRUCT);
							soillayer.litterSolveSOM.add_litter(cmass_litter_leaf * fm, nmass_litter_leaf * fm, SURFMETA);
						}

						// Fire
						if (firemodel == GLOBFIRM) {

							litterme[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * ppft.pft.litterme;
							fireresist[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * ppft.pft.fireresist;

							litterme[SURFMETA] += cmass_litter_leaf * fm * ppft.pft.litterme;
							fireresist[SURFMETA] += cmass_litter_leaf * fm * ppft.pft.fireresist;
						}

						if (negligible(soillayer.sompool[SURFSTRUCT].cmass, -16)) {
							soillayer.sompool[SURFSTRUCT].ligcfrac = 0.0;
						}
						else {
							soillayer.sompool[SURFSTRUCT].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SURFSTRUCT].cmass));
						}
					}

					// ROOT

					if (!negligible(cmass_litter_root) || !negligible(nmass_litter_root)) {

						// Calculate inputs to soil structural and metabolic litter

						// Root litter lignin:N ratio
						double root_lton = lignin_to_n_ratio(cmass_litter_root, nmass_litter_root, ppft.pft.ligcfrac_root, ppft.pft.cton_root_avr);

						// Metabolic litter fraction for root litter
						double fm = metabolic_litter_fraction(root_lton);

						ligcmass_new = cmass_litter_root * root_lyr_frac * (1.0 - fm) * ppft.pft.ligcfrac_root;
						ligcmass_old = soillayer.sompool[SOILSTRUCT].cmass * soillayer.sompool[SOILSTRUCT].ligcfrac;

						// Add to pools and update lignin fraction in structural pool
						soillayer.sompool[SOILSTRUCT].cmass += cmass_litter_root * root_lyr_frac * (1.0 - fm);
						soillayer.sompool[SOILSTRUCT].nmass += nmass_litter_root * root_lyr_frac * (1.0 - fm);
						soillayer.sompool[SOILMETA].cmass += cmass_litter_root * root_lyr_frac * fm;
						soillayer.sompool[SOILMETA].nmass += nmass_litter_root * root_lyr_frac * fm;

						soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_root * root_lyr_frac;
						soilcolumn.csom_litter_input[SOILSTRUCT] += cmass_litter_root * root_lyr_frac * (1.0 - fm);
						soilcolumn.csom_litter_input[SOILMETA] += cmass_litter_root * root_lyr_frac * fm;

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_root * root_lyr_frac * (1.0 - fm), nmass_litter_root * root_lyr_frac * (1.0 - fm), SOILSTRUCT);
							soillayer.litterSolveSOM.add_litter(cmass_litter_root * root_lyr_frac * fm, nmass_litter_root * root_lyr_frac * fm, SOILMETA);
						}

						if (negligible(soillayer.sompool[SOILSTRUCT].cmass, -16)) {
							soillayer.sompool[SOILSTRUCT].ligcfrac = 0.0;
						}
						else {
							soillayer.sompool[SOILSTRUCT].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SOILSTRUCT].cmass));
						}
					}

					// WOOD

					// Woody debris enters two woody litter pools as described in
					// Kirschbaum and Paul (2002).

					// SAP WOOD

					if (!negligible(cmass_litter_sap) || !negligible(nmass_litter_sap)) {

						// Surface sap wood
						if (lyr == 0) {

							ligcmass_new = cmass_litter_sap * (1.0 - FCWD_ROOT) * ppft.pft.ligcfrac_wood;
							ligcmass_old = soillayer.sompool[SURFFWD].cmass * soillayer.sompool[SURFFWD].ligcfrac;

							// Add to fine woody pool and update lignin fraction in pool
							soillayer.sompool[SURFFWD].cmass += cmass_litter_sap * (1.0 - FCWD_ROOT);
							soillayer.sompool[SURFFWD].nmass += nmass_litter_sap * (1.0 - FCWD_ROOT);

							soilcolumn.cmass_litter_input[0] += cmass_litter_sap * (1.0 - FCWD_ROOT);
							soilcolumn.csom_litter_input[SURFFWD] += cmass_litter_sap * (1.0 - FCWD_ROOT);

							// Save litter input for equilsom()
							if (soilcolumn.savesomcent) {
								soillayer.litterSolveSOM.add_litter(cmass_litter_sap * (1.0 - FCWD_ROOT), nmass_litter_sap * (1.0 - FCWD_ROOT), SURFFWD);
							}

							if (negligible(soillayer.sompool[SURFFWD].cmass, -16)) {
								soillayer.sompool[SURFFWD].ligcfrac = 0.0;
							}
							else {
								soillayer.sompool[SURFFWD].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
									soillayer.sompool[SURFFWD].cmass));
							}

							// Fire
							if (firemodel == GLOBFIRM) {
								litterme[SURFFWD] += cmass_litter_sap * (1.0 - FCWD_ROOT) * ppft.pft.litterme;
								fireresist[SURFFWD] += cmass_litter_sap * (1.0 - FCWD_ROOT) * ppft.pft.fireresist;
							}
						}

						// Soil sap wood (fine woody roots)

						ligcmass_new = cmass_litter_sap * FCWD_ROOT * root_lyr_frac * ppft.pft.ligcfrac_wood;
						ligcmass_old = soillayer.sompool[SOILFWD].cmass * soillayer.sompool[SOILFWD].ligcfrac;

						// Add to structural pool and update lignin fraction in pool
						soillayer.sompool[SOILFWD].cmass += cmass_litter_sap * FCWD_ROOT * root_lyr_frac;
						soillayer.sompool[SOILFWD].nmass += nmass_litter_sap * FCWD_ROOT * root_lyr_frac;

						soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_sap * FCWD_ROOT * root_lyr_frac;
						soilcolumn.csom_litter_input[SOILFWD] += cmass_litter_sap * FCWD_ROOT * root_lyr_frac;

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_sap * FCWD_ROOT * root_lyr_frac, nmass_litter_sap * FCWD_ROOT * root_lyr_frac, SOILFWD);
						}

						if (negligible(soillayer.sompool[SOILFWD].cmass, -16)) {
							soillayer.sompool[SOILFWD].ligcfrac = 0.0;
						}
						else {
							soillayer.sompool[SOILFWD].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SOILFWD].cmass));
						}
					}

					// HEART WOOD

					if (!negligible(cmass_litter_heart) || !negligible(nmass_litter_heart)) {

						// Surface heart wood
						if (lyr == 0) {

							ligcmass_new = cmass_litter_heart * (1.0 - FCWD_ROOT) * ppft.pft.ligcfrac_wood;
							ligcmass_old = soillayer.sompool[SURFCWD].cmass * soillayer.sompool[SURFCWD].ligcfrac;

							// Add to fine woody pool and update lignin fraction in pool
							soillayer.sompool[SURFCWD].cmass += cmass_litter_heart * (1.0 - FCWD_ROOT);
							soillayer.sompool[SURFCWD].nmass += nmass_litter_heart * (1.0 - FCWD_ROOT);

							soilcolumn.cmass_litter_input[0] += cmass_litter_heart * (1.0 - FCWD_ROOT);
							soilcolumn.csom_litter_input[SURFCWD] += cmass_litter_heart * (1.0 - FCWD_ROOT);

							// Save litter input for equilsom()
							if (soilcolumn.savesomcent) {
								soillayer.litterSolveSOM.add_litter(cmass_litter_heart * (1.0 - FCWD_ROOT), nmass_litter_heart * (1.0 - FCWD_ROOT), SURFCWD);
							}

							if (negligible(soillayer.sompool[SURFCWD].cmass, -16)) {
								soillayer.sompool[SURFCWD].ligcfrac = 0.0;
							}
							else {
								soillayer.sompool[SURFCWD].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
									soillayer.sompool[SURFCWD].cmass));
							}

							// Fire
							if (firemodel == GLOBFIRM) {
								litterme[SURFCWD] += cmass_litter_heart * (1.0 - FCWD_ROOT) * ppft.pft.litterme;
								fireresist[SURFCWD] += cmass_litter_heart * (1.0 - FCWD_ROOT) * ppft.pft.fireresist;
							}
						}

						// Soil heart wood (fine woody roots)

						ligcmass_new = cmass_litter_heart * FCWD_ROOT * root_lyr_frac * ppft.pft.ligcfrac_wood;
						ligcmass_old = soillayer.sompool[SOILCWD].cmass * soillayer.sompool[SOILCWD].ligcfrac;

						// Add to structural pool and update lignin fraction in pool
						soillayer.sompool[SOILCWD].cmass += cmass_litter_heart * FCWD_ROOT * root_lyr_frac;
						soillayer.sompool[SOILCWD].nmass += nmass_litter_heart * FCWD_ROOT * root_lyr_frac;

						soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_heart * FCWD_ROOT * root_lyr_frac;
						soilcolumn.csom_litter_input[SOILCWD] += cmass_litter_heart * FCWD_ROOT * root_lyr_frac;

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_heart * FCWD_ROOT * root_lyr_frac, nmass_litter_heart * FCWD_ROOT * root_lyr_frac, SOILCWD);
						}

						if (negligible(soillayer.sompool[SOILCWD].cmass, -16)) {
							soillayer.sompool[SOILCWD].ligcfrac = 0.0;
						}
						else {
							soillayer.sompool[SOILCWD].ligcfrac = max(0.01, min(0.99, (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SOILCWD].cmass));
						}
					}
				}
			}

			patch.pft.nextobj();
		}

		// Fire only in the surface soil layer

		// Soil object for surface layer
		Soillayer& soillayer = soilcolumn.soillayers[0];

		// FIRE
		if (firemodel == GLOBFIRM) {
			if (soillayer.sompool[SURFSTRUCT].cmass > 0.0) {
				soillayer.sompool[SURFSTRUCT].litterme = litterme[SURFSTRUCT] / soillayer.sompool[SURFSTRUCT].cmass;
				soillayer.sompool[SURFSTRUCT].fireresist = fireresist[SURFSTRUCT] / soillayer.sompool[SURFSTRUCT].cmass;
			}
			if (soillayer.sompool[SURFMETA].cmass > 0.0) {
				soillayer.sompool[SURFMETA].litterme = litterme[SURFMETA] / soillayer.sompool[SURFMETA].cmass;
				soillayer.sompool[SURFMETA].fireresist = fireresist[SURFMETA] / soillayer.sompool[SURFMETA].cmass;
			}
			if (soillayer.sompool[SURFFWD].cmass > 0.0) {
				soillayer.sompool[SURFFWD].litterme = litterme[SURFFWD] / soillayer.sompool[SURFFWD].cmass;
				soillayer.sompool[SURFFWD].fireresist = fireresist[SURFFWD] / soillayer.sompool[SURFFWD].cmass;
			}
			if (soillayer.sompool[SURFCWD].cmass > 0.0) {
				soillayer.sompool[SURFCWD].litterme = litterme[SURFCWD] / soillayer.sompool[SURFCWD].cmass;
				soillayer.sompool[SURFCWD].fireresist = fireresist[SURFCWD] / soillayer.sompool[SURFCWD].cmass;
			}
		}

		// Calculate total litter carbon and nitrogen mass for set N:C ratio of surface microbial pool
		double litter_cmass = soillayer.sompool[SURFSTRUCT].cmass + soillayer.sompool[SURFMETA].cmass +
			soillayer.sompool[SURFFWD].cmass + soillayer.sompool[SURFCWD].cmass;
		double litter_nmass = soillayer.sompool[SURFSTRUCT].nmass + soillayer.sompool[SURFMETA].nmass +
			soillayer.sompool[SURFFWD].nmass + soillayer.sompool[SURFCWD].nmass;

		// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
		// Parton et al 1993 Fig 4. Dry mass litter == cmass litter * 2
		if (!negligible(litter_cmass)) {
			soillayer.setntoc(litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
		}
	}

	// Add litter to solvesom array every month
	if (patch.soilcolumn.savesomcent && date.islastday) {

		//dprintf("Year %d Day %d SURFMETA C %g N %g SURFCWD C %g N %g\n",date.year,date.day,
		//	patch.soilcolumn.soillayers[0].litterSolveSOM.get_clitter(SURFMETA), patch.soilcolumn.soillayers[0].litterSolveSOM.get_nlitter(SURFMETA),
		//	patch.soilcolumn.soillayers[0].litterSolveSOM.get_clitter(SURFCWD), patch.soilcolumn.soillayers[0].litterSolveSOM.get_nlitter(SURFCWD));

		for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
			patch.soilcolumn.soillayers[lyr].solvesom.push_back(patch.soilcolumn.soillayers[lyr].litterSolveSOM);
			patch.soilcolumn.soillayers[lyr].litterSolveSOM.clear();
		}
	}
}

/// LEACHING
/** Leaching fractions for both organic and mineral leaching
 *  Should be called every day in both daily and monthly mode
 */
void leaching_LPJ(Soilcolumn& soilcolumn) {

	// Mass balance check
	// double ncont_before = soilcolumn.get_sum(NCONT);
	// double nextra = 0.0;

	double minleachfrac, minmovefrac, runoff;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		if (lyr == 0) {
			runoff = soillayer.runoff_drain + soilcolumn.patch.fluxes.get_daily_flux(Fluxes::WRUNOFF_SURF, date.day);
		}
		else {
			runoff = soillayer.runoff_drain;
		}

		if (!negligible(soillayer.perc_mass + runoff)) {

			// Leaching from available nitrogen mineral pool
			// in proportion to amount of water precolated and drained
			minleachfrac = runoff / (soillayer.perc_mass + runoff + soillayer.wmass);
			minmovefrac = soillayer.perc_mass / (soillayer.perc_mass + runoff + soillayer.wmass);

			// Leaching from decayed organic carbon/nitrogen
			// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
			// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
			// water_leaching was originally in cm/month
			double OMLECH_1 = 0.03;
			double OMLECH_2 = 0.12;
			double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H2O/month)
			double cmpermonth_to_mmperday = 10.0 * 12.0 / 365.0;

			// precolation and drainage results in direct loss of organic C
			soillayer.orgleachfrac = min(1.0, (soillayer.perc_mass + runoff) / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soilcolumn.soiltype.sand_frac);
		}
		else {
			minleachfrac = 0.0;
			minmovefrac = 0.0;
			soillayer.orgleachfrac = 0.0;
		}

		// Leaching of soil mineral nitrogen
		// Allowed on days with residual nitrogen following vegetation uptake
		// in proportion to amount of water drainage
		const double nmin_avail = soillayer.nmass_avail(NO3);
		if (nmin_avail > 0.0) {

			double leaching = nmin_avail * minleachfrac;

			soillayer.nmass_subtract(leaching, NO3);
			soillayer.aminleach += leaching;

			// Mass balance 
			// nextra += leaching;

			leaching = nmin_avail * minmovefrac;

			soillayer.nmass_subtract(leaching, NO3);

			// Bottom layer then N is lost
			if (lyr == soilcolumn.nsoillayers - 1) {
				soillayer.aminleach += leaching;
				// Mass balance
				//nextra += leaching;
			}
			// else if is transferred to layer below
			else {
				soilcolumn.soillayers[lyr + 1].nmass_inc(leaching, NO3);
			}
		}

		if (soilcolumn.savesomcent) {
			soillayer.morgleach_mean[date.month] += soillayer.orgleachfrac / date.ndaymonth[date.month];
			soillayer.mminleach_mean[date.month] += minleachfrac / date.ndaymonth[date.month];
		}
	}

	// Mass balance
	// double ncont_after = soilcolumn.get_sum(NCONT) + nextra;

	// if (date.year > freenyears && fabs(ncont_before - ncont_after) > 1.0e-14)
	//	dprintf("Year %d Day %d NITROGEN LEACH diff %g %g %g\n", date.year, date.day, fabs(ncont_before - ncont_after), ncont_before, ncont_after);
}

/// LEACHING
/** Leaching fractions for both organic and mineral leaching
*  Should be called every day in both daily and monthly mode
*/
void leaching_LSM(Soilcolumn& soilcolumn) {

	// Water fluxes between layers during (mm s-1)
	// First element is infiltration from suface
	// Last element is baseflow out of the bottom layer
	// Negative value means upwards movement of water to 
	// layer above (lyr - 1)
	std::vector<double> wfluxes = soilcolumn.wfluxes;

	// Nitrogen moved to other layers
	std::vector<double> NO3_delta(soilcolumn.nsoillayers);

	double runoff;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Water runoff for this time step
		if (lyr == 0) {
			runoff = soillayer.runoff_drain + soilcolumn.patch.fluxes.get_daily_flux(Fluxes::WRUNOFF_SURF, date.day);
		}
		else {
			runoff = soillayer.runoff_drain;
		}

		double minleachfrac, minmovefrac;

		// Water percolation for this time step
		double dperc = wfluxes[lyr + 1];

		if (!negligible(fabs(dperc) + runoff, -10)) {

			// Leaching from available nitrogen mineral pool
			// in proportion to amount of water percolated or drainage
			minleachfrac = runoff / (fabs(dperc) + runoff + soillayer.wmass + soillayer.imass);
			minmovefrac = dperc / (fabs(dperc) + runoff + soillayer.wmass + soillayer.imass);

			// Leaching from decayed organic carbon/nitrogen
			// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
			// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
			// water_leaching was originally in cm/month
			double OMLECH_1 = 0.03;
			double OMLECH_2 = 0.12;
			double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H20/month)
			double cmpermonth_to_mmperday = MM_PER_CM * MONTH_PER_YEAR / date.year_length();

			soillayer.orgleachfrac = min(1.0, runoff / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soilcolumn.soiltype.sand_frac);
		}
		else {
			minleachfrac = 0.0;
			minmovefrac = 0.0;
			soillayer.orgleachfrac = 0.0;
		}

		// Leaching of soil mineral nitrogen
		// Allowed on days with residual nitrogen following vegetation uptake
		// in proportion to amount of water drainage
		const double nmin_avail = soillayer.nmass_avail(NO3);
		if (nmin_avail > 0.0) {

			double leach = nmin_avail * minleachfrac;
			double move = nmin_avail * fabs(minmovefrac);
			NO3_delta[lyr] -= (leach + move);
			soillayer.aminleach += leach;
			int move_lyr = min(soilcolumn.nsoillayers, max(0, (minmovefrac > 0.0) ? lyr + 1 : lyr - 1));
			if (move_lyr == soilcolumn.nsoillayers) {
				soillayer.aminleach += move;
			}
			else {
				NO3_delta[move_lyr] += move;
			}
		}

		if (soilcolumn.savesomcent) {
			soillayer.morgleach_mean[date.month] += soillayer.orgleachfrac / date.ndaymonth[date.month];
			soillayer.mminleach_mean[date.month] += fabs(minleachfrac) / date.ndaymonth[date.month];
		}
	}

	// Add nmass_avail that has been moved around in soil column
	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++)
		soilcolumn.soillayers[lyr].nmass_inc(NO3_delta[lyr], NO3);
}

void leaching(Patch& patch) {

	if (iftwolayersoil)
		leaching_LPJ(patch.soilcolumn);
	else
		leaching_LSM(patch.soilcolumn);
}


/// Nitrogen addition to the soil
/** Daily nitrogen addition to the soil
 *  from deposition and fixation
 */
void soilnadd(Patch& patch) {

	Soilcolumn& soilcolumn = patch.soilcolumn;

	if (!ifntransform) {
		// Add an additional N flux from soil
		double nflux = 0.01 * (soilcolumn.NH4_input + soilcolumn.NO3_input);
		soilcolumn.NH4_input -= 0.01 * soilcolumn.NH4_input;
		soilcolumn.NO3_input -= 0.01 * soilcolumn.NO3_input;

		patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);
	}

	// Nitrogen deposition and fertilization input to the soil (calculated in snow_ninput())
	soilcolumn.soillayers[0].NH4_mass += soilcolumn.NH4_input;
	soilcolumn.soillayers[0].NO3_mass += soilcolumn.NO3_input;
	
	// Nitrogen fixation
	// If soil available nitrogen is above the value for minimum SOM C:N ratio, then
	// nitrogen fixation is reduced (nitrogen rich soils)

	// Get soil layers total C mass
	double cconttot = soilcolumn.get_sum(CCONT);

	if (!negligible(cconttot)) {
		// Nitrogen fixation in each soil layer
		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			Soillayer& soillayer = soilcolumn.soillayers[lyr];

			// Nitrogen mass saturation level
			double NMASS_SAT_LYR = NMASS_SAT;

			// NH4 in soil layer
			const double nmin_avail = soillayer.get(NAVAIL);

			// Nitrogen fixation
			// If soil available nitrogen is above the value for minimum SOM C:N ratio, then
			// nitrogen fixation is reduced (nitrogen rich soils)
			if (nmin_avail < NMASS_SAT_LYR) {

				double daily_nfix = soilcolumn.anfix_calc * (soillayer.get(CCONT) / cconttot) / date.year_length();

				if (nmin_avail + daily_nfix > NMASS_SAT_LYR) {
					daily_nfix = NMASS_SAT_LYR - nmin_avail;
				}

				soillayer.nmass_inc(daily_nfix, NH4);
				soilcolumn.anfix += daily_nfix;
			}
		}
	}

	// Calculate nitrogen fixation (Cleveland et al. 1999)
	// by using five year average aaet
	if (date.islastmonth && date.islastday) {

		// Add this year's AET to aaet_5 which keeps track of the last 5 years
		patch.aaet_5.add(patch.fluxes.get_annual_flux(Fluxes::WAET) + patch.fluxes.get_annual_flux(Fluxes::WEVAP) + patch.fluxes.get_annual_flux(Fluxes::WINTERCEP));

		// Calculate estimated maximum nitrogen fixation (aaet should be in cm/yr, eqn is in nitrogen/ha/yr)
		double cmtomm = 0.1;
		double hatom2 = 0.0001;
		soilcolumn.anfix_calc = max((nfix_a * patch.aaet_5.mean() * cmtomm + nfix_b) * hatom2, 0.0);

		if (soilcolumn.savesomcent) {
			soilcolumn.anfix_mean += soilcolumn.anfix;
		}
	}
}

/// Vertical mixing (diffusion) or movement of soil organic matter
/** Vertical mixing (diffusion) or movement of soil organic matter
*  Cryoturbation processes on the vertical profile of SOM stocks are captured
*  by considering a simple diffusion scheme.
*  And if organic C exceed maximum in soil layer, SOM is moved downwards
*/
void vertical_movement(Patch& patch) {

	// Vertical movement of organic matter if C density of soil layer exceed maximum
	if (date.islastday && date.islastmonth)
		patch.soilcolumn.vertical_movement();


	if (patch.stand.landcover == PEATLAND || !date.islastday || !(ifverticalmixing || date.year < patch.soilcolumn.solvesomcent_beginyr))
		return;

	// Low carbon in soil -> no vertical mixing
	if (date.year <= patch.soilcolumn.solvesomcent_endyr || patch.soilcolumn.get_sum(CCONT) < 1.0e-3)
		return;

	// Average monthly active layer depth
	double mald = min(patch.soilcolumn.depth, max(0.0, patch.soilcolumn.mthaw[date.month])) * M_PER_MM;

	if (patch.soilcolumn.savesomcent) {
		patch.soilcolumn.mald_mean[date.month] += mald;
	}

	// Store pool sizes before mixing
	std::vector<double> cpools(2);
	std::vector<double> npools(2);

	for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
		cpools[0] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;
		npools[0] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass;
		cpools[1] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;
		npools[1] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass;
	}

	// Timestep in Crank-Nicholson method 
	// Diffusion is on a yearly timestep and as vertical mixing is done monthly
	// the timestep will be 1 / 12
	double timestep = 1.0 / 12.0;

	//patch.soilcolumn.vertical_mixing(SOILFWD, mald, timestep);
	//patch.soilcolumn.vertical_mixing(SOILCWD, mald, timestep);
	patch.soilcolumn.vertical_mixing(SLOWSOM, mald, timestep);
	patch.soilcolumn.vertical_mixing(PASSIVESOM, mald, timestep);


	// Check pool sizes after mixing
	std::vector<double> cpools_after(2);
	std::vector<double> npools_after(2);

	for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
		cpools_after[0] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;
		npools_after[0] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass;
		cpools_after[1] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;
		npools_after[1] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass;
	}

	// Might be small rounding errors
	// Add to/subtract from top pool
	patch.soilcolumn.soillayers[0].sompool[SLOWSOM].cmass += cpools[0] - cpools_after[0];
	patch.soilcolumn.soillayers[0].sompool[SLOWSOM].nmass += npools[0] - npools_after[0];
	patch.soilcolumn.soillayers[0].sompool[PASSIVESOM].cmass += cpools[1] - cpools_after[1];
	patch.soilcolumn.soillayers[0].sompool[PASSIVESOM].nmass += npools[1] - npools_after[1];
}

/// Add litter to equilsom()
void add_litter(Soillayer& soillayer, int timestep, int pool) {

	// If LUC have occurred in between solvesomcent_beginyr and solvesomcent_endyr then array might be too smalll
	if (timestep >= (int)soillayer.solvesom.size())
		timestep = timestep % (int)soillayer.solvesom.size();

	// Add to litter pools
	soillayer.sompool[pool].cmass += soillayer.solvesom[timestep].get_clitter(pool);
	soillayer.sompool[pool].nmass += soillayer.solvesom[timestep].get_nlitter(pool);
}

/// Iteratively solving differential flux equations for century SOM pools
/** Iteratively solving differential flux equations for century SOM pools
 *  assuming annual litter inputs, nitrogen uptake and leaching is close
 *  to long term equilibrium
 */
void equilsom(Patch& patch) {

	// Do large SOM size prediction on year solvesomcent_endyr
	// after solvesomcent_endyr every 10 year until spinup-100 extrapolate trends
	if (!((date.year == patch.soilcolumn.solvesomcent_endyr || (date.year >= patch.soilcolumn.solvesomcent_endyr + 5 && date.year <= nyear_spinup - 100 && !(date.year % 5))) && date.islastmonth && date.islastday))
		return;

	Soilcolumn& soilcolumn = patch.soilcolumn;
	const Climate& climate = soilcolumn.patch.get_climate();
	const Gridcell& gridcell = patch.stand.get_gridcell();

	// Main equilsom that happens once
	if (date.year == patch.soilcolumn.solvesomcent_endyr) {

		std::vector<double> save_NH4_avail(soilcolumn.nsoillayers);
		std::vector<double> save_NO3_avail(soilcolumn.nsoillayers);

		// Number of years with mean input data
		int EQUILSOM_YEARS = soilcolumn.solvesomcent_endyr - soilcolumn.solvesomcent_beginyr + 1;

		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			Soillayer& soillayer = soilcolumn.soillayers[lyr];

			// Save soil mineral nitrogen status
			save_NH4_avail[lyr] = soillayer.nmass_avail(NH4);
			save_NO3_avail[lyr] = soillayer.nmass_avail(NO3);

			for (int m = 0; m < 12; m++) {

				// Monthly average decay rates
				for (int p = 0; p < NSOMPOOL; p++) {
					soillayer.sompool[p].mfracremain_mean[m] = pow(min(1.0, soillayer.sompool[p].mfracremain_mean[m] / EQUILSOM_YEARS), (double)date.ndaymonth[m]);
				}

				// Monthly average mineral nitrogen uptake
				soillayer.fNH4uptake_mean[m] /= EQUILSOM_YEARS;
				soillayer.fNO3uptake_mean[m] /= EQUILSOM_YEARS;

				// Monthly average organic carbon and nitrogen leaching
				soillayer.morgleach_mean[m] /= EQUILSOM_YEARS;

				// Monthly average mineral nitrogen leaching
				soillayer.mminleach_mean[m] /= EQUILSOM_YEARS;

				// Monthly average active layer depth
				if (lyr == 0)
					soilcolumn.mald_mean[m] /= EQUILSOM_YEARS;
			}

			// Annual average nitrogen fixation
			if (lyr == 0)
				soilcolumn.anfix_mean /= EQUILSOM_YEARS;
		}

		// Vertical mixing time step (every 4 year)
		const int VMTIMESTEP = 4;
		// Maximum number of iterations the trend needs to be the same for the multiplier to increase
		const int SAME_TREND_MAX = 5;
		// Maximum multiplier
		double MAX_MULTI = pow(2.0, 15);
		// Minimum multiplier
		const double MIN_MULTI = 0.01;
		// Multiplier modifier
		const double MULTI_MOD = 2.0;
		// Minimum tries of iterations
		const int MIN_TRIES = 200;
		// Maximum tries of iterations
		const int MAX_TRIES = 1000;
		// Maximum trend / std
		const double EPS_trend = 1.0e-9;
		const double EPS_trend_passive = 1.0e-11;
		// Number of trends that doesn't meet the criteria
		int bad_trend = 1;
		// Iteration tries
		int tries = 0;
		// SOM pool N:C ratio
		double ntoc;
		// C mass increase
		double cmass_inc;
		// Max negative C mass change
		double max_neg_change;
		// Max C change for SOM pools
		double max_change;
		// Total SOM
		double cconttot;

		// Initial SOM pool size per layer
		std::vector< std::vector<double> > init(soilcolumn.nsoillayers, std::vector<double>(NSOMPOOL));
		// Multiplier for SOM pools per layer
		std::vector< std::vector<double> > multi(soilcolumn.nsoillayers, std::vector<double>(NSOMPOOL));
		// Current sign of trend for SOM pools per layer
		std::vector< std::vector<int> > sign(soilcolumn.nsoillayers, std::vector<int>(NSOMPOOL));
		// Number of iterations with same trend sign for SOM pools per layer
		std::vector< std::vector<int> > same_trend(soilcolumn.nsoillayers, std::vector<int>(NSOMPOOL));

		// Historical values Soilc CWD, slow and passive SOM, and total soillayer C
		std::vector< Historic<double, 20> > hist_slow(soilcolumn.nsoillayers);
		std::vector< Historic<double, 20> > hist_passive(soilcolumn.nsoillayers);
		std::vector< Historic<double, 20> > hist_ccont(soilcolumn.nsoillayers);

		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			for (int p = 0; p < NSOMPOOL; p++) {
				multi[lyr][p] = MAX_MULTI;
				sign[lyr][p] = 1;
				same_trend[lyr][p] = 0;
			}
		}

		int peat_years = 0;

		// Spinup SOM pools until trends of slow and passive SOM pools 
		// are within limits or maximum tries have been done
		while ((bad_trend && tries < MAX_TRIES) || tries < MIN_TRIES) {

			// Update soil layers total C mass fractions
			cconttot = soilcolumn.get_sum(CCONT);

			for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {
				for (int p = 0; p < NSOMPOOL; p++) {
					if (soilcolumn.soillayers[lyr].sompool[p].cmass < 0.0)
						dprintf("CMASS neg lyr %d pool %d C %g\n", lyr, p, soilcolumn.soillayers[lyr].sompool[p].cmass);
					if (soilcolumn.soillayers[lyr].sompool[p].nmass < 0.0)
						dprintf("NMASS neg lyr %d pool %d N %g\n", lyr, p, soilcolumn.soillayers[lyr].sompool[p].nmass);
				}
			}

			// Spin SOM pools with saved litter input, nitrogen addition and fractions of
			// nitrogen uptake and leaching for NYEARS years with monthly timesteps
			for (int yr = 0; yr < EQUILSOM_YEARS; yr++) {

				// To do vertical mixing this year
				int vertmix = yr%VMTIMESTEP;

				// Monthly time steps
				for (int m = 0; m < 12; m++) {

					// Go through each soil layer
					for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

						Soillayer& soillayer = soilcolumn.soillayers[lyr];

						// Initial C for SOM pools in this iteration
						if (yr == 0 && m == 0) {
							for (int p = 0; p < NSOMPOOL; p++) {
								init[lyr][p] = soillayer.sompool[p].cmass;
							}
						}

						// Transfer monthly mean litter (added a maximum litter pool size so a runaway pool doesn't occur)
						int timestep = yr * 12 + m;
						if (lyr == 0) {
							add_litter(soillayer, timestep, SURFSTRUCT);
							add_litter(soillayer, timestep, SURFMETA);
							add_litter(soillayer, timestep, SURFFWD);
							add_litter(soillayer, timestep, SURFCWD);
						}
						add_litter(soillayer, timestep, SOILSTRUCT);
						add_litter(soillayer, timestep, SOILMETA);
						add_litter(soillayer, timestep, SOILFWD);
						add_litter(soillayer, timestep, SOILCWD);

						// Calculate total litter carbon and nitrogen mass for set N:C ratio of surface microbial pool
						if (lyr == 0) {
							double litter_cmass = soillayer.sompool[SURFSTRUCT].cmass + soillayer.sompool[SURFMETA].cmass +
								soillayer.sompool[SURFFWD].cmass + soillayer.sompool[SURFCWD].cmass;
							double litter_nmass = soillayer.sompool[SURFSTRUCT].nmass + soillayer.sompool[SURFMETA].nmass +
								soillayer.sompool[SURFFWD].nmass + soillayer.sompool[SURFCWD].nmass;

							// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
							if (!negligible(litter_cmass)) {
								soillayer.setntoc(litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
							}
						}

						// Monthly nitrogen uptake
						soillayer.NH4_mass *= (1.0 - soillayer.fNH4uptake_mean[m]);
						soillayer.NO3_mass *= (1.0 - soillayer.fNO3uptake_mean[m]);

						// Monthly mineral nitrogen leaching
						double nmass = soillayer.nmass_avail(NO3);
						double nleach = nmass * (1.0 - soillayer.mminleach_mean[m]);
						soillayer.nmass_subtract(nleach, NO3);

						// Transferred to layer below if not bottom layer then N is lost
						if (lyr < soilcolumn.nsoillayers - 1) {
							soilcolumn.soillayers[lyr + 1].nmass_inc(nleach, NO3);
						}

						// Monthly nitrogen addition to the system
						if (lyr == 0) {
							soillayer.nmass_inc(gridcell.aNH4dep / 12.0, NH4);
							soillayer.nmass_inc(gridcell.aNO3dep / 12.0, NO3);
						}

						// Monthly nitrogen addition from fixation to the soil layer
						if (!negligible(cconttot, -12)) // DW_COMMENT
							soillayer.nmass_inc(soilcolumn.anfix_mean / 12.0 * soillayer.get(CCONT) / cconttot, NH4);

						// Monthly decomposition and fluxes between SOM pools

						// Set this months decay rates
						for (int p = 0; p < NSOMPOOL; p++) {
							soillayer.sompool[p].fracremain = soillayer.sompool[p].mfracremain_mean[m];
						}

						// Set this months organic matter leaching
						soillayer.orgleachfrac = soillayer.morgleach_mean[m];
					}

					// Monthly decomposition and fluxes between SOM pools
					// and nitrogen flux from soil
					somfluxes(patch, true, false);

					if (patch.stand.landcover != PEATLAND && !vertmix && ifverticalmixing) {

						soilcolumn.vertical_mixing(SLOWSOM, soilcolumn.mald_mean[m], (double)VMTIMESTEP / 12.0);
						soilcolumn.vertical_mixing(PASSIVESOM, soilcolumn.mald_mean[m], (double)VMTIMESTEP / 12.0);
					}
				}

				// Vertical movement of organic matter if C density of soil layer exceed maximum
				patch.soilcolumn.vertical_movement();
			}

			// Reset number of trends that doesn't meet the criteria
			bad_trend = 0;

			// Increase number of iterations
			tries++;

			// If peat age is defined, then only do soil spinup for peat_age minus nyear_spinup
			if (peat_age && patch.stand.landcover == PEATLAND) {
				peat_years += EQUILSOM_YEARS;

				// If peat age is to short, then set bad_trend to 1 so while loop continues
				if (peat_years < peat_age - nyear_spinup)
					bad_trend = 1;
				else
					bad_trend = 0;
			}
			else {

				// Loop through layers and try and to solve equilibrium by appling the trend of pools
				// by several order of magnitude
				for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

					for (int p = 0; p < NSOMPOOL; p++) {

						if (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM || p == PASSIVESOM) {
							if (tries >= 100) {
								max_change = soilcolumn.soillayers[lyr].sompool[p].cmass * 0.5;
								max_neg_change = -max_change;
							}
							else {
								max_change = 0.5 * soilcolumn.soillayers[lyr].depth * M_PER_MM;
								max_neg_change = -soilcolumn.soillayers[lyr].sompool[p].cmass * 0.5;
							}

							if (soilcolumn.soillayers[lyr].sompool[p].cmass >= 1.0e-16) {

								// Trend over NYEARS years
								double trend = soilcolumn.soillayers[lyr].sompool[p].cmass - init[lyr][p];

								multi[lyr][p] = min(multi[lyr][p], MAX_MULTI);

								// Trend changed sign
								if (trend * sign[lyr][p] < 0.0) {

									// Save trend sign and set same trend to zero
									sign[lyr][p] *= -1;

									// Flipped sign maximum times then reduce multiplier
									if (multi[lyr][p] > MIN_MULTI) {
										multi[lyr][p] /= MULTI_MOD;
									}

									// Flipping so no same trend
									same_trend[lyr][p] = 0;
								}
								// Same sign trend
								else {

									// Add to number of times same trend has happen in a raw
									same_trend[lyr][p]++;

									// Reached maximum same trend then increase multiplier
									if (same_trend[lyr][p] > SAME_TREND_MAX && multi[lyr][p] < MAX_MULTI) {
										multi[lyr][p] *= MULTI_MOD;
										same_trend[lyr][p] = 0;
									}
								}

								// C change
								cmass_inc = min(max_change, max(max_neg_change, trend * multi[lyr][p]));

								// N:C ratio of SOM pool
								ntoc = soilcolumn.soillayers[lyr].sompool[p].nmass / soilcolumn.soillayers[lyr].sompool[p].cmass;

								// Add C change
								soilcolumn.soillayers[lyr].sompool[p].cmass += cmass_inc;

								// Update N pool
								soilcolumn.soillayers[lyr].sompool[p].nmass = ntoc * soilcolumn.soillayers[lyr].sompool[p].cmass;
							}
						}
					}

					// Save C pools in historical arrays					
					hist_slow[lyr].add(soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass);
					hist_passive[lyr].add(soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass);
					hist_ccont[lyr].add(soilcolumn.soillayers[lyr].ccont(true));

					// Calculate historical mean
					double hslow_mean = hist_slow[lyr].mean();
					double hpas_mean = hist_passive[lyr].mean();
					double hccont_mean = hist_ccont[lyr].mean();

					// Calculate trend normalised against historical mean
					double hsslow = !negligible(hslow_mean, -12) ? hist_slow[lyr].slope() / hslow_mean : 0.0;
					double hspas = !negligible(hpas_mean, -12) ? hist_passive[lyr].slope() / hpas_mean : 0.0;
					double hsccont = !negligible(hccont_mean, -12) ? hist_ccont[lyr].slope() / hccont_mean : 0.0;

					// Calculate std normalised against historical mean
					double hstdslow = !negligible(hslow_mean, -12) ? hist_slow[lyr].std() / hslow_mean : 0.0;
					double hstdpas = !negligible(hpas_mean, -12) ? hist_passive[lyr].std() / hpas_mean : 0.0;
					double hstdccont = !negligible(hccont_mean, -12) ? hist_ccont[lyr].std() / hccont_mean : 0.0;

					// Calculte how many trends and stds that are too big
					if ((iftwolayersoil || lyr > 0) && (fabs(hsslow) > EPS_trend || fabs(hspas) > EPS_trend_passive || fabs(hsccont) > EPS_trend) &&
						(fabs(hstdslow) > EPS_trend || fabs(hstdpas) > EPS_trend_passive || fabs(hstdccont) > EPS_trend)) {

						bad_trend++;
					}
				}

				// Maximum multiplier should go down over time.
				if (!(tries % 20)) {
					MAX_MULTI = max(8.0, MAX_MULTI / MULTI_MOD);
				}
			}
		}

		dprintf("EquilSOM Tries %d Landcover %d Soil C %g\n", tries, patch.stand.landcover, soilcolumn.get_sum(CCONT));		

		// Reset mineral nitrogen status
		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			Soillayer& soillayer = soilcolumn.soillayers[lyr];

			soillayer.NH4_mass = save_NH4_avail[lyr];
			soillayer.NO3_mass = save_NO3_avail[lyr];

			// Reset C and N litter status
			for (int p = 0; p < NSOMPOOL; p++) {

				soillayer.sompool[p].cmass = max(0.0, soillayer.sompool[p].cmass);
				soillayer.sompool[p].nmass = max(0.0, soillayer.sompool[p].nmass);

				// Save cmass for trend based equilsom for every 10 year
				soillayer.sompool[p].cmass_init = soillayer.sompool[p].cmass;
			}

			// Reset variables for next possible equilsom()
			for (int m = 0; m < 12; m++) {

				for (int p = 0; p < NSOMPOOL; p++) {
					soillayer.sompool[p].mfracremain_mean[m] = 0.0;
				}

				soillayer.fNH4uptake_mean[m] = 0.0;
				soillayer.fNO3uptake_mean[m] = 0.0;
				soillayer.morgleach_mean[m] = 0.0;
				soillayer.mminleach_mean[m] = 0.0;
			}

			soillayer.solvesom.clear();
			std::vector<LitterSolveSOM>().swap(soillayer.solvesom); // clear array memory
		}

		// Reset soil column annual N fixation
		soilcolumn.anfix_mean = 0.0;
	}
	// Trend based equilsom (not if peat age is defined)
	else if (!(peat_age && patch.stand.landcover == PEATLAND)) {

		// Maximum number of iterations the trend needs to be the same for the multiplier to increase
		const int SAME_TREND_MAX = 3;
		// Maximum number of iterations the trend needs to flip for the multiplier to decrease
		const int FLIP_MAX = 2;
		// Maximum multiplier
		const double MAX_MULTI = 8.0;
		double MAX_MULTI_passive = date.year <= patch.soilcolumn.solvesomcent_endyr + 100 ? 1024.0 :
								   date.year <= patch.soilcolumn.solvesomcent_endyr + 200 ? 128 : 16;
		// Minimum multiplier
		const double MIN_MULTI = 0.125;
		// Multiplier modifier
		const double MULTI_MOD = 2.0;
		// SOM pool N:C ratio
		double ntoc;
		// C mass increase
		double cmass_inc;
		// Maximum positive C change for SOM pool
		double max_pos_change;
		// Maximum negative C change for SOM pool
		double max_neg_change;

		// Loop through layers and try and to solve equilibrium by appling the trend of pools
		// by several order of magnitude
		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			Soillayer& soillayer = soilcolumn.soillayers[lyr];

			double ccont = soillayer.ccont(true);

			for (int p = 0; p < NSOMPOOL; p++) {
				// SOM pool needs to be of some size
				if (soillayer.sompool[p].cmass > 1.0e-12 && (p == SURFHUMUS || p == SLOWSOM || p == PASSIVESOM)) {

					// Initial C and N mass
					double cmass_init = soillayer.sompool[p].cmass;
					double nmass_init = soillayer.sompool[p].nmass;

					// Maximum positive and negative change
					max_pos_change = soillayer.sompool[p].cmass_init * 0.1;
					max_neg_change = -soillayer.sompool[p].cmass * 0.1;

					// C mass trend
					double trend = soillayer.sompool[p].cmass - soillayer.sompool[p].cmass_init;

					// If trend changes sign a certain number of times then reduce multiplier
					// If trend has the same sign a certain number of times then increase multiplier

					// Trend changed sign
					if (trend * soillayer.sompool[p].sign_trend < 0.0) {

						// Save trend sign and set same trend to zero
						soillayer.sompool[p].sign_trend *= -1;
						soillayer.sompool[p].changed_trend++;

						// Flipped sign maximum times then reduce multiplier
						if (soillayer.sompool[p].changed_trend >= FLIP_MAX && soillayer.sompool[p].multi > MIN_MULTI) {
							soillayer.sompool[p].multi /= MULTI_MOD;
						}

						// Flipping so no same trend
						soillayer.sompool[p].same_trend = 0;
					}
					// Same sign trend
					else {

						// Add to number of times same trend has happen in a raw
						soillayer.sompool[p].same_trend++;

						// Reached maximum same trend then increase multiplier
						if (soillayer.sompool[p].same_trend >= SAME_TREND_MAX && ((soillayer.sompool[p].multi < MAX_MULTI && p != PASSIVESOM) || (soillayer.sompool[p].multi < MAX_MULTI_passive && p == PASSIVESOM))) {
							soillayer.sompool[p].multi *= MULTI_MOD;
						}

						// Same trend so no flipping
						soillayer.sompool[p].changed_trend = 0;
					}

					// C change
					cmass_inc = min(max_pos_change, max(max_neg_change, trend * soillayer.sompool[p].multi));

					// N:C ratio of SOM pool
					ntoc = soillayer.sompool[p].nmass / soillayer.sompool[p].cmass;

					// Add C change
					soillayer.sompool[p].cmass += cmass_inc;

					// Update N pool
					soillayer.sompool[p].nmass = ntoc * soillayer.sompool[p].cmass;

					if (soillayer.sompool[p].cmass < 1.0e-13 || soillayer.sompool[p].nmass < 1.0e-14) {
						soillayer.sompool[p].cmass = max(0.0, cmass_init);
						soillayer.sompool[p].nmass = max(0.0, nmass_init);
					}

					// Save C mass as initial value for next equilsom trend
					soillayer.sompool[p].cmass_init = soillayer.sompool[p].cmass;
				}
				else if (soillayer.sompool[p].cmass < 1.0e-13 || soillayer.sompool[p].nmass < 1.0e-14) {
					soillayer.sompool[p].cmass = 0.0;
					soillayer.sompool[p].nmass = 0.0;
				}
			}
		}
	}
}


/// SOM CENTURY DYNAMICS
/** To be called each simulation day for each modelled patch, following update
 *  of soil temperature and soil water.
 *  Transfers litter, performes nitrogen uptake and addition, leaching and decomposition.
 */
void som_dynamics_century(Patch& patch, Climate& climate, bool tillage) {

	// Transfer litter to SOM pools
	transfer_litter(patch);

	// Daily mineral and organic nitrogen leaching
	leaching(patch);

	// Daily nitrogen addition to the soil
	soilnadd(patch);

	// Daily or monthly decomposition and fluxes between SOM pools
	somfluxes(patch, false, tillage);

	// Nitrogen transformation in soil
	ntransform(patch, climate);

	// Vertical movement of soil organic matter
	vertical_movement(patch);

	// Solve SOM pool sizes at end of year given by soilcolumn.solvesomcent_endyr
	equilsom(patch);
}

/// Choose between CENTURY or standard LPJ SOM dynamics
/**
*/
void som_dynamics(Patch& patch, Climate& climate) {

	bool tillage = iftillage && patch.stand.landcover == CROPLAND;
	if (ifcentury) {
		som_dynamics_century(patch, climate, tillage);
	}
	else {
		som_dynamics_lpj(patch, tillage);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Chatskikh, D., Hansen, S., Olesen, J.E. & Petersen, B.M. 2009. A simplified modelling approach
//	 for quantifying tillage effects on soil carbon stocks. Eur.J.Soil.Sci., 60:924-934.
// CENTURY Soil Organic Matter Model Version 5; Century User's Guide and Reference.
//  http://www.nrel.colostate.edu/projects/century5/reference/index.htm; accessed Oct.7, 2016.
// Chaudhary, N., Miller, P.A., and Smith, B.: Modelling Holocene peatland dynamics with an 
//   individual - based dynamic vegetation model, Biogeosciences, 14, 2571–2596, 
//   https://doi.org/10.5194/bg-14-2571-2017, 2017.
// Comins, H. N. & McMurtrie, R. E. 1993. Long-Term Response of Nutrient-Limited
//   Forests to CO2 Enrichment - Equilibrium Behavior of Plant-Soil Models.
//   Ecological Applications, 3, 666-681.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Cleveland C C (1999) Global patterns of terrestrial biological nitrogen (N2) fixation
//   in natural ecosystems. GBC 13: 623-645
// Foley J A 1995 An equilibrium model of the terrestrial carbon budget
//   Tellus (1995), 47B, 310-319
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A
//   process-based, terrestrial biosphere model of ecosystem dynamics
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Kirschbaum, M. U. F. and K. I. Paul (2002). "Modelling C and N dynamics in forest soils
//   with a modified version of the CENTURY model." Soil Biology & Biochemistry 34(3): 341-354.
// Ise, T., Dunn, A., Wofsy, S. et al. High sensitivity of peat decomposition to climate change 
//   through water-table feedback. Nature Geosci 1, 763–766 (2008). https://doi.org/10.1038/ngeo331
// Meentemeyer, V. (1978) Macroclimate and lignin control of litter decomposition
//   rates. Ecology 59: 465-472.
// Mikan, C. J., Schimel, J. P., and Doyle, A. P. 2002. Temperature controls of microbial respiration 
//   in arctic tundra soils above and below freezing, Soil Biol.Biochem., 34, 1785–1795, 2002.
// Parton, W. J., Scurlock, J. M. O., Ojima, D. S., Gilmanov, T. G., Scholes, R. J., Schimel, D. S.,
//   Kirchner, T., Menaut, J. C., Seastedt, T., Moya, E. G., Kamnalrut, A. & Kinyamario, J. I. 1993.
//   Observations and Modeling of Biomass and Soil Organic-Matter Dynamics for the Grassland Biome
//   Worldwide. Global Biogeochemical Cycles, 7, 785-809.
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W. & Kelly, R. 2010.
//   ForCent model development and testing using the Enriched Background Isotope Study experiment.
//   Journal of Geophysical Research-Biogeosciences, 115.
// Schaefer, K. and Jafarov, E. 2016. 
//   A parameterization of respiration in frozen soils based on substrate availability
//   Biogeosciences, 13, 1991 - 2001, https ://doi.org/10.5194/bg-13-1991-2016

