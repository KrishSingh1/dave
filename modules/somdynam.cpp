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

	moist_response = 0.25 + 0.75 * min(1.0, wcont);

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

	// Obtain reference to Soil object
	Soilcolumn& soilcolumn = patch.soilcolumn;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Calculate decay constants and rates given today's soil moisture and
		// temperature

		decayrates(soillayer.wcont(), soillayer.gtemp, k_soilfast, k_soilslow, fr_litter,
			fr_soilfast, fr_soilslow, tillage); // DW_COMMENT tillage depth

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

			decomp_litter += patch.pft[p].litter_root * (1.0 - fr_litter) * patch.pft[p].pft.rootdist[lyr];

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

		// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
		if (patch.stand.landcover != PEATLAND)
			patch.fluxes.report_flux(Fluxes::SOILC, cflux);
		else
			soillayer.dcflux_soil = cflux;
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

	// Calculate growth fraction decrease to satify N deficit

	// Calculate total growth fraction after reduced decay rates
	double tot_cinc = tot_cdec - tot_resp;
	double tot_growthfrac = tot_cinc / tot_cdec;

	// Calculate needed growth fraction to satify a neutral net mineralisation
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

/// Calculates CENTURY instantaneous decay rates
/** Calculates CENTURY instantaneous decay rates given soil layer temperature,
*  water content, texture etc
*/
void Soillayer::decayrates_century(bool tillage, bool ispeatland, bool ismineralwetland, bool savesomcent) {

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris decay
	// pools SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFFWD,SURFCWD,SOILFWD,SOILCWD,SOILMETA,SLOWSOM,PASSIVESOM
	const double K_MAX[] = { 9.5e-3, 1.9e-2, 4.2e-2, 4.8e-4, 2.7e-2, 3.8e-2, 1.1e-2, 2.2e-3, 1.1e-2, 2.2e-3, 7.0e-2, 1.7e-3, 1.9e-6 };
	
	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	const double texture_mod = 1.0 - 0.75 * ((*soilcolumn).soiltype.clay_frac + (*soilcolumn).soiltype.silt_frac);

	// Calculate decomposition temperature modifier (in range 0-1)
	// [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
	// effect of soil temperature, Parton et al 1993, Fig 2)

	double temp_mod = 0.0;

	if (temp > 0.0) {
		temp_mod = max(0.0,
			0.0326 + 0.00351 * pow(temp, 1.652) - pow(temp / 41.748, 7.19));
	}

	// Include as an overide option when: MIN_DECOMP_TEMP < temp < 0 degC
	// This increases the respiration (from 0) between MIN_DECOMP_TEMP and 0 degC - cf Koven et al. 2011
	if (ifcarbonfreeze && temp <= 0.0 && temp >= MIN_DECOMP_TEMP && !iftwolayersoil) {

		double decomp_at_freezing_point = 0.0326; // temp_mod above when temp_soil = 0;
		bool linear_decrease_below_freezing = false;

		if (linear_decrease_below_freezing) {
			// Linear approach (Koven et al. 2011)
			double slope = decomp_at_freezing_point / fabs(MIN_DECOMP_TEMP);
			temp_mod = slope * temp + decomp_at_freezing_point; // i.e. a linear decrease from decomp_at_freezing_point at 0C to 0 at MIN_DECOMP_TEMP (-4C).
		}
		else {
			// Alternative Q10 relationship (Schaefer & Jafarov, 2016)
			double q10_freeze = 2.9; // (Natali et al., 2019)
			temp_mod = decomp_at_freezing_point * pow(q10_freeze, temp / 10.0);
		}
	}

	// Calculate decomposition moisture modifier (in range 0-1)
	// Friend et al 1997, Eqn 53 and Parton et al 1993, Fig 2

	// Water Filled Pore Spaces (wfps) %
	// water holding capacity at wilting point (wp) and saturation capacity (wsats)
	// is calculated with the help of Cosby et al 1984;
	// use Gerten equivalents here, but wfps COULD be made depth equivalent
	const double wfps100 = wfps() * 100.0;

	double moist_mod;

	if (wfps100 < 60.0)
		moist_mod = exp((wfps100 - 60.0) * (wfps100 - 60.0) / -800.0);
	else
		moist_mod = 0.000371 * wfps100 * wfps100 - 0.0748 * wfps100 + 4.13;

	// Additional moist modifier for saturated soils
	// no effect unless this is a peat soil
	double moist_mod_saturated = 1.0; 

	// Moisture modifiers for peat soils
	if (ispeatland) {

		moist_mod = RMOIST;

		double cmass_total = 0.0;
		for (int p = 0; p < NSOMPOOL; p++) {
			cmass_total += sompool[p].cmass;
		}

		double acrotelm_climit = 2.5 * depth * DM_PER_MM; // kgC/m2 - 7.5 kgC m-2 is max C content in a 30cm-deep acrotelm - see Wania et al. (2009b) so 2.5 kg C m-2 per 10cm
		if (cmass_total > acrotelm_climit) { // take a weighted average of the aerobic and anaerobic moisture modifiers 
			moist_mod = (acrotelm_climit * RMOIST + (cmass_total - acrotelm_climit) * RMOIST_ANAEROBIC) / cmass_total;
			// moist_mod approaches a value of RMOIST_ANAEROBIC asymptotically as cmass_total increases
		}

		moist_mod_saturated = RMOIST_ANAEROBIC / moist_mod;
	}

	// Moisture modifiers for wetlands on mineral soils
	if (ismineralwetland) {
		// 100% WFPS for wetlands on mineral soils, 0.36 approx
		moist_mod = 0.000371 * 100.0 * 100.0 - 0.0748 * 100 + 4.13;
	}

	// Combined moisture and temperature modifier

	for (int p = 0; p < NSOMPOOL; p++) {

		// Calculate decay constant
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		double k = K_MAX[p] * temp_mod * moist_mod * depth_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2
		// Kirschbaum et al 2001 changed the exponential term
		// from 3 to 5.

		if (p == SURFSTRUCT || p == SOILSTRUCT || p == SURFFWD || p == SURFCWD) {
			k *= exp(-5.0 * sompool[p].ligcfrac);
		}
		else if (p == SOILMICRO && !ispeatland) {
			k *= texture_mod;
		}

		// Increased decay rates for crops (tillage)
		if (tillage && (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM) && !ispeatland) {
			k *= TILLAGE_FACTOR;
		}

		// Reduced decomposition for the passive and slow pools in peatlands, as they are assumed to be in the catotelm
		if (p == PASSIVESOM || p == SLOWSOM) {
			k *= moist_mod_saturated; // ensures that a modifier of RMOIST_ANAEROBIC is used.
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

	double leachsum_cmass, leachsum_nmass;
	// Actual (not net) nitrogen mineralisation
	double nmin_actual;
	// Nitrogen immobilisation
	double nimmob;		

	const double EPS = 1.0e-16;

	double ncont_before = soilcolumn.ncont();
	double ccont_before = soilcolumn.ccont();
	double nextra = 0.0;
	double cextra = 0.0;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		double cmass_before = soillayer.ccont();
		double nmass_before = soillayer.ncont();

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

		// Set N:C ratios for humus, soil microbial, passive and slow pool based on estimated mineral nitrogen pool
		// (Parton et al 1993, Fig 4)

		// ForCent (Parton 2010) values
		soillayer.setntoc(nmin_mass, SLOWSOM, 30.0, 15.0, 0.0, NMASS_SAT);
		soillayer.setntoc(nmin_mass, SOILMICRO, 15.0, 6.0, 0.0, NMASS_SAT);
		soillayer.setntoc(nmin_mass, SURFHUMUS, 30.0, 15.0, 0.0, NMASS_SAT);

		if (!ifequilsom) {

			// Calculate potential fraction remaining following decay today for all pools
			// (assumes no nitrogen limitation)
			soillayer.decayrates_century(tillage, patch.stand.is_highlatitude_peatland_stand(), patch.stand.is_true_wetland_stand(), soilcolumn.savesomcent);
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
		}
		double init_negative_nmass, init_ntoc_reduction;
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

			// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
			if (patch.stand.landcover != PEATLAND)
				patch.fluxes.report_flux(Fluxes::SOILC, respsum_tot);
			else
				soillayer.dcflux_soil = respsum_tot;

			cextra += respsum_tot;

			soillayer.arespsum += respsum_tot;

			// Sum annual organic nitrogen leaching
			soillayer.aorgNleach += leachsum_nmass;
			nextra += leachsum_nmass;

			// Sum annual organic carbon leaching
			soillayer.aorgCleach += leachsum_cmass;
			cextra += leachsum_cmass;

			// Sum annuals
			soillayer.anmin += nmin_actual;
			soillayer.animmob += nimmob;
		}

		// Fraction of microbial resp. is assumed to produce labile carbon
		soillayer.labile_carbon = respsum_tot * frac_labile_carbon;

		// Adding mineral nitrogen to soil available pool
		double nmin_inc = nmin_actual - nimmob;

		if (!ifequilsom) {
			soillayer.mnetmin[date.month] += nmin_inc;
		}
		double nflux_tmp = 0.0;
		// Estimate of N flux from soil (simple CLM-CN approach)
		if (!ifntransform) {
			double nflux = nmin_inc > 0.0 ? (nmin_inc) * 0.01 : 0.0;
			nmin_inc -= nflux;
			nextra += nflux;
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

		double cmass_after = soillayer.ccont() + respsum_tot + leachsum_cmass;
		double nmass_after = soillayer.ncont() + nflux_tmp + leachsum_nmass;

		if (date.year > freenyears && fabs(cmass_before - cmass_after) > 1.0e-12 && !ifequilsom) // DW_COMMENT remove
			dprintf("Year %d Day %d Lyr %d C SOM diff %g %g %g Leach %g\n", date.year, date.day, lyr, fabs(cmass_before - cmass_after), cmass_before, cmass_after, leachsum_cmass);

		if (date.year > freenyears && fabs(nmass_before - nmass_after) > 1.0e-12 && !ifequilsom) // DW_COMMENT remove
			dprintf("Year %d Day %d Lyr %d N SOM diff %g %g %g Leach %g\n", date.year, date.day, lyr, fabs(nmass_before - nmass_after), nmass_before, nmass_after, leachsum_nmass);
	}

	double ncont_after = soilcolumn.ncont() + nextra;
	double ccont_after = soilcolumn.ccont() + cextra;

	if (date.year > freenyears && fabs(ncont_before - ncont_after) > 1.0e-12 && !ifequilsom) // DW_COMMENT remove
		dprintf("Year %d Day %d NITROGEN SOM diff %g %g %g\n", date.year, date.day, fabs(ncont_before - ncont_after), ncont_before, ncont_after);

	if (date.year > freenyears && fabs(ccont_before - ccont_after) > 1.0e-12 && !ifequilsom) // DW_COMMENT remove
		dprintf("Year %d Day %d CARBON SOM diff %g %g %g\n", date.year, date.day, fabs(ccont_before - ccont_after), ccont_before, ccont_after);

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

	// For CROPLAND only enter on harvest/turnover day
	if (patch.stand.landcover == CROPLAND && !patch.is_litter_day)
		return;

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
		Patchpft& pft = patch.pft.getobj();

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
			cmass_litter_leaf = pft.litter_leaf;
			nmass_litter_leaf = pft.nmass_litter_leaf;
			cmass_litter_root = pft.litter_root;
			nmass_litter_root = pft.nmass_litter_root;

			pft.litter_leaf = 0.0;
			pft.nmass_litter_leaf = 0.0;
			pft.litter_root = 0.0;
			pft.nmass_litter_root = 0.0;
		}
		// For summergreens drop leaf litter over all days during Jan in NH and July SH
		else if ((pft.pft.phenology == SUMMERGREEN && ((lat >= 0.0 && date.month == 0) || (lat < 0.0 && date.month == 6))) ||
			(pft.pft.phenology == RAINGREEN && date.month == pft.driest_mth)) {
			double frac = 1.0 / (date.ndaymonth[date.month] - date.dayofmonth);
			cmass_litter_leaf = pft.litter_leaf * frac;
			nmass_litter_leaf = pft.nmass_litter_leaf * frac;
			cmass_litter_root = pft.litter_root * frac;
			nmass_litter_root = pft.nmass_litter_root * frac;

			pft.litter_leaf -= cmass_litter_leaf;
			pft.nmass_litter_leaf -= nmass_litter_leaf;
			pft.litter_root -= cmass_litter_root;
			pft.nmass_litter_root -= nmass_litter_root;
		}
		// Drops leaf litter every day
		else {
			double frac = 1.0 / (date.year_length() - date.day);
			cmass_litter_leaf = pft.litter_leaf * frac;
			nmass_litter_leaf = pft.nmass_litter_leaf * frac;
			cmass_litter_root = pft.litter_root * frac;
			nmass_litter_root = pft.nmass_litter_root * frac;

			pft.litter_leaf -= cmass_litter_leaf;
			pft.nmass_litter_leaf -= nmass_litter_leaf;
			pft.litter_root -= cmass_litter_root;
			pft.nmass_litter_root -= nmass_litter_root;
		}

		if (pft.pft.lifeform == TREE) {

			// Woody litter is dropped as a fraction every day
			double frac = 1.0 / (date.year_length() - date.day);
			cmass_litter_sap = pft.litter_sap * frac;
			nmass_litter_sap = pft.nmass_litter_sap * frac;
			cmass_litter_heart = pft.litter_heart * frac;
			nmass_litter_heart = pft.nmass_litter_heart * frac;

			pft.litter_sap -= cmass_litter_sap;
			pft.nmass_litter_sap -= nmass_litter_sap;
			pft.litter_heart -= cmass_litter_heart;
			pft.nmass_litter_heart -= nmass_litter_heart;
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
					double leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, LIGCFRAC_LEAF, pft.pft.cton_leaf_avr);

					// Metabolic litter fraction for leaf litter
					double fm = metabolic_litter_fraction(leaf_lton);

					ligcmass_old = soillayer.sompool[SURFSTRUCT].cmass * soillayer.sompool[SURFSTRUCT].ligcfrac;

					// Add to pools
					soillayer.sompool[SURFSTRUCT].cmass += cmass_litter_leaf * (1.0 - fm);
					soillayer.sompool[SURFSTRUCT].nmass += nmass_litter_leaf * (1.0 - fm);
					soillayer.sompool[SURFMETA].cmass += cmass_litter_leaf * fm;
					soillayer.sompool[SURFMETA].nmass += nmass_litter_leaf * fm;

					soilcolumn.cmass_litter_input[0] += cmass_litter_leaf;

					// Save litter input for equilsom()
					if (soilcolumn.savesomcent) {
						soillayer.litterSolveSOM.add_litter(cmass_litter_leaf * (1.0 - fm), nmass_litter_leaf * (1.0 - fm), SURFSTRUCT);
						soillayer.litterSolveSOM.add_litter(cmass_litter_leaf * fm, nmass_litter_leaf * fm, SURFMETA);
					}

					// Fire
					if (firemodel == GLOBFIRM) {

						litterme[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.litterme;
						fireresist[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.fireresist;

						litterme[SURFMETA] += cmass_litter_leaf * fm * pft.pft.litterme;
						fireresist[SURFMETA] += cmass_litter_leaf * fm * pft.pft.fireresist;
					}
					// NB: reproduction litter cannot contain nitrogen!!

					ligcmass_new = cmass_litter_leaf * (1.0 - fm) * LIGCFRAC_LEAF;

					if (negligible(soillayer.sompool[SURFSTRUCT].cmass)) {
						soillayer.sompool[SURFSTRUCT].ligcfrac = 0.0;
					}
					else {
						soillayer.sompool[SURFSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
							soillayer.sompool[SURFSTRUCT].cmass;
					}
				}

				// ROOT

				if (!negligible(cmass_litter_root) || !negligible(nmass_litter_root)) {

					// Calculate inputs to soil structural and metabolic litter

					// Root litter lignin:N ratio
					double root_lton = lignin_to_n_ratio(cmass_litter_root, nmass_litter_root, LIGCFRAC_ROOT, pft.pft.cton_root_avr);

					// Metabolic litter fraction for root litter
					double fm = metabolic_litter_fraction(root_lton);

					ligcmass_new = cmass_litter_root * (1.0 - fm) * LIGCFRAC_ROOT;
					ligcmass_old = soillayer.sompool[SOILSTRUCT].cmass * soillayer.sompool[SOILSTRUCT].ligcfrac;

					// Add to pools and update lignin fraction in structural pool
					soillayer.sompool[SOILSTRUCT].cmass += cmass_litter_root * pft.pft.rootdist[lyr] * (1.0 - fm);
					soillayer.sompool[SOILSTRUCT].nmass += nmass_litter_root * pft.pft.rootdist[lyr] * (1.0 - fm);
					soillayer.sompool[SOILMETA].cmass += cmass_litter_root * pft.pft.rootdist[lyr] * fm;
					soillayer.sompool[SOILMETA].nmass += nmass_litter_root * pft.pft.rootdist[lyr] * fm;

					soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_root * pft.pft.rootdist[lyr];

					// Save litter input for equilsom()
					if (soilcolumn.savesomcent) {
						soillayer.litterSolveSOM.add_litter(cmass_litter_root * pft.pft.rootdist[lyr] * (1.0 - fm), nmass_litter_root * pft.pft.rootdist[lyr] * (1.0 - fm), SOILSTRUCT);
						soillayer.litterSolveSOM.add_litter(cmass_litter_root * pft.pft.rootdist[lyr] * fm, nmass_litter_root * pft.pft.rootdist[lyr] * fm, SOILMETA);
					}

					if (negligible(soillayer.sompool[SOILSTRUCT].cmass)) {
						soillayer.sompool[SOILSTRUCT].ligcfrac = 0.0;
					}
					else {
						soillayer.sompool[SOILSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
							soillayer.sompool[SOILSTRUCT].cmass;
					}
				}

				// WOOD

				// Woody debris enters two woody litter pools as described in
				// Kirschbaum and Paul (2002).

				// SAP WOOD

				if (!negligible(cmass_litter_sap) || !negligible(nmass_litter_sap)) {

					// Surface sap wood
					if (lyr == 0) {

						ligcmass_new = cmass_litter_sap * (1.0 - FCWD_ROOT) * LIGCFRAC_WOOD;
						ligcmass_old = soillayer.sompool[SURFFWD].cmass * soillayer.sompool[SURFFWD].ligcfrac;

						// Add to fine woody pool and update lignin fraction in pool
						soillayer.sompool[SURFFWD].cmass += cmass_litter_sap * (1.0 - FCWD_ROOT);
						soillayer.sompool[SURFFWD].nmass += nmass_litter_sap * (1.0 - FCWD_ROOT);

						soilcolumn.cmass_litter_input[0] += cmass_litter_sap * (1.0 - FCWD_ROOT);

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_sap * (1.0 - FCWD_ROOT), nmass_litter_sap * (1.0 - FCWD_ROOT), SURFFWD);
						}

						if (negligible(soillayer.sompool[SURFFWD].cmass)) {
							soillayer.sompool[SURFFWD].ligcfrac = 0.0;
						}
						else {
							double ligcfrac = (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SURFFWD].cmass;
							soillayer.sompool[SURFFWD].ligcfrac = ligcfrac;
						}

						// Fire
						if (firemodel == GLOBFIRM) {
							litterme[SURFFWD] += cmass_litter_sap * (1.0 - FCWD_ROOT) * pft.pft.litterme;
							fireresist[SURFFWD] += cmass_litter_sap * (1.0 - FCWD_ROOT) * pft.pft.fireresist;
						}
					}

					// Soil sap wood (fine woody roots)

					ligcmass_new = cmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr] * LIGCFRAC_WOOD;
					ligcmass_old = soillayer.sompool[SOILFWD].cmass * soillayer.sompool[SOILFWD].ligcfrac;

					// Add to structural pool and update lignin fraction in pool
					soillayer.sompool[SOILFWD].nmass += nmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr];
					soillayer.sompool[SOILFWD].cmass += cmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr];

					soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr];

					// Save litter input for equilsom()
					if (soilcolumn.savesomcent) {
						soillayer.litterSolveSOM.add_litter(cmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr], nmass_litter_sap * FCWD_ROOT * pft.pft.rootdist[lyr], SOILFWD);
					}

					if (negligible(soillayer.sompool[SOILFWD].cmass)) {
						soillayer.sompool[SOILFWD].ligcfrac = 0.0;
					}
					else {
						double ligcfrac = (ligcmass_new + ligcmass_old) /
							soillayer.sompool[SOILFWD].cmass;
						soillayer.sompool[SOILFWD].ligcfrac = ligcfrac;
					}
				}

				// HEART WOOD

				if (!negligible(cmass_litter_heart) || !negligible(nmass_litter_heart)) {

					// Surface heart wood
					if (lyr == 0) {

						ligcmass_new = cmass_litter_heart * (1.0 - FCWD_ROOT) * LIGCFRAC_WOOD;
						ligcmass_old = soillayer.sompool[SURFCWD].cmass * soillayer.sompool[SURFCWD].ligcfrac;

						// Add to fine woody pool and update lignin fraction in pool
						soillayer.sompool[SURFCWD].cmass += cmass_litter_heart * (1.0 - FCWD_ROOT);
						soillayer.sompool[SURFCWD].nmass += nmass_litter_heart * (1.0 - FCWD_ROOT);

						soilcolumn.cmass_litter_input[0] += cmass_litter_heart * (1.0 - FCWD_ROOT);

						// Save litter input for equilsom()
						if (soilcolumn.savesomcent) {
							soillayer.litterSolveSOM.add_litter(cmass_litter_heart * (1.0 - FCWD_ROOT), nmass_litter_heart * (1.0 - FCWD_ROOT), SURFCWD);
						}

						if (negligible(soillayer.sompool[SURFCWD].cmass)) {
							soillayer.sompool[SURFCWD].ligcfrac = 0.0;
						}
						else {
							double ligcfrac = (ligcmass_new + ligcmass_old) /
								soillayer.sompool[SURFCWD].cmass;
							soillayer.sompool[SURFCWD].ligcfrac = ligcfrac;
						}

						// Fire
						if (firemodel == GLOBFIRM) {
							litterme[SURFCWD] += cmass_litter_heart * (1.0 - FCWD_ROOT) * pft.pft.litterme;
							fireresist[SURFCWD] += cmass_litter_heart * (1.0 - FCWD_ROOT) * pft.pft.fireresist;
						}
					}

					// Soil heart wood (fine woody roots)

					ligcmass_new = cmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr] * LIGCFRAC_WOOD;
					ligcmass_old = soillayer.sompool[SOILCWD].cmass * soillayer.sompool[SOILCWD].ligcfrac;

					// Add to structural pool and update lignin fraction in pool
					soillayer.sompool[SOILCWD].nmass += nmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr];
					soillayer.sompool[SOILCWD].cmass += cmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr];

					soilcolumn.cmass_litter_input[lyr + 1] += cmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr];

					// Save litter input for equilsom()
					if (soilcolumn.savesomcent) {
						soillayer.litterSolveSOM.add_litter(cmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr], nmass_litter_heart * FCWD_ROOT * pft.pft.rootdist[lyr], SOILCWD);
					}

					if (negligible(soillayer.sompool[SOILCWD].cmass)) {
						soillayer.sompool[SOILCWD].ligcfrac = 0.0;
					}
					else {
						double ligcfrac = (ligcmass_new + ligcmass_old) /
							soillayer.sompool[SOILCWD].cmass;
						soillayer.sompool[SOILCWD].ligcfrac = ligcfrac;
					}
				}
			}
		}

		patch.pft.nextobj();
	}

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		// Soil object for this layer
		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Fire only in the surface soil layer
		if (lyr == 0) {

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
		if (soilcolumn.savesomcent) {
			if (date.dayofmonth == 0) {
				soillayer.solvesom.push_back(soillayer.litterSolveSOM);
				soillayer.litterSolveSOM.clear();
			}
		}
	}
}

/// LEACHING
/** Leaching fractions for both organic and mineral leaching
 *  Should be called every day in both daily and monthly mode
 */
void leaching_LPJ(Soilcolumn& soilcolumn) {

	double ncont_before = soilcolumn.ncont();
	double nextra = 0.0;

	double minleachfrac, minmovefrac;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		if (!negligible(soillayer.perc_mass + soillayer.runoff_drain)) {

			// Leaching from available nitrogen mineral pool
			// in proportion to amount of water precolated and drained
			minleachfrac = soillayer.runoff_drain / (soillayer.perc_mass + soillayer.runoff_drain + soillayer.wmass);
			minmovefrac = soillayer.perc_mass / (soillayer.perc_mass + soillayer.runoff_drain + soillayer.wmass);

			// Leaching from decayed organic carbon/nitrogen
			// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
			// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
			// water_leaching was originally in cm/month
			double OMLECH_1 = 0.03;
			double OMLECH_2 = 0.12;
			double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H2O/month)
			double cmpermonth_to_mmperday = 10.0 * 12.0 / 365.0;

			// precolation and drainage results in direct loss of organic C
			soillayer.orgleachfrac = min(1.0, (soillayer.perc_mass + soillayer.runoff_drain) / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soilcolumn.soiltype.sand_frac);
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
			nextra += leaching;

			leaching = nmin_avail * minmovefrac;

			soillayer.nmass_subtract(leaching, NO3);

			// Bottom layer then N is lost
			if (lyr == soilcolumn.nsoillayers - 1) {
				soillayer.aminleach += leaching;
				nextra += leaching;
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

	double ncont_after = soilcolumn.ncont() + nextra;

	if (date.year > freenyears && fabs(ncont_before - ncont_after) > 1.0e-14) // DW_COMMENT remove
		dprintf("Year %d Day %d NITROGEN LEACH diff %g %g %g\n", date.year, date.day, fabs(ncont_before - ncont_after), ncont_before, ncont_after);
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
	std::vector<double> NO3_delta(soilcolumn.nsoillayers); // DW_COMMENT

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		double minleachfrac, minmovefrac;

		// Water percolation and drainage for this time step
		double dperc = wfluxes[lyr + 1];
		double drain = soillayer.runoff_drain;

		if (!negligible(fabs(dperc) + drain, -10)) {
		//if (!negligible(fabs(dperc) + drain + soillayer.wmass, -10)) { // DW_COMMENT

			// Leaching from available nitrogen mineral pool
			// in proportion to amount of water percolated or drainage
			minleachfrac = drain / (fabs(dperc) + drain + soillayer.wmass);
			minmovefrac = dperc / (fabs(dperc) + drain + soillayer.wmass);

			// Leaching from decayed organic carbon/nitrogen
			// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
			// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
			// water_leaching was originally in cm/month
			double OMLECH_1 = 0.03;
			double OMLECH_2 = 0.12;
			double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H20/month)
			double cmpermonth_to_mmperday = 10.0 * 12.0 / 365.0;

			soillayer.orgleachfrac = min(1.0, drain / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soilcolumn.soiltype.sand_frac);
			soillayer.orgmovefrac = 0.0;// min(1.0, dperc / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soillayer.sand_frac);
		}
		else {
			minleachfrac = 0.0;
			minmovefrac = 0.0;
			soillayer.orgleachfrac = 0.0;
			soillayer.orgmovefrac = 0.0;
		}

		// Leaching of soil mineral nitrogen
		// Allowed on days with residual nitrogen following vegetation uptake
		// in proportion to amount of water drainage
		const double nmin_avail = soillayer.nmass_avail(NO3, true);
		if (nmin_avail > 0.0) {

			double leach = nmin_avail * minleachfrac;
			double move = nmin_avail * fabs(minmovefrac);
			NO3_delta[lyr] -= (leach + move);
			soillayer.aminleach += leach;
			int move_lyr = min(soilcolumn.nsoillayers - 1, max(0, (minmovefrac > 0.0) ? lyr + 1 : lyr - 1));
			NO3_delta[move_lyr] += move;
		}

		if (soilcolumn.savesomcent) {
			soillayer.morgleach_mean[date.month] += soillayer.orgleachfrac / date.ndaymonth[date.month];
			soillayer.morgmove_mean[date.month] += soillayer.orgmovefrac / date.ndaymonth[date.month];
			soillayer.mminleach_mean[date.month] += fabs(minleachfrac) / date.ndaymonth[date.month];
		}
	}

	// Add nmass_avail that has been moved around in soil column
	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++)
		soilcolumn.soillayers[lyr].nmass_inc(NO3_delta[lyr], NO3);
}

void leaching(Patch& patch) { // DW_COMMENT should be able to be the same function!!!!!

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

	// Get soil layers total C mass fractions
	std::vector<double> ccontfrac = soilcolumn.get_ccontfrac();

	// Nitrogen fixation in each soil layer
	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Nitrogen mass saturation level
		double NMASS_SAT_LYR = NMASS_SAT * ccontfrac[lyr];

		// NH4 in soil layer
		const double nmin_avail = soillayer.nmass_avail(NH4);

		// Nitrogen fixation
		// If soil available nitrogen is above the value for minimum SOM C:N ratio, then
		// nitrogen fixation is reduced (nitrogen rich soils)
		if (nmin_avail < NMASS_SAT_LYR) {

			double daily_nfix = soilcolumn.anfix_calc * ccontfrac[lyr] / date.year_length();

			if (nmin_avail + daily_nfix > NMASS_SAT_LYR) {
				daily_nfix = NMASS_SAT_LYR - nmin_avail;
			}

			soillayer.nmass_inc(daily_nfix, NH4);
			soilcolumn.anfix += daily_nfix;
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

/// Vertical mixing (diffusion) of soil organic matter
/** Vertical mixing (diffusion) of soil organic matter
*  Cryoturbation processes on the vertical profile of SOM stocks are captured
*  by considering a simple diffusion scheme.
*/
void vertical_mixing(Patch& patch) {

	if (!date.islastday || !ifverticalmixing || date.year < patch.soilcolumn.solvesomcent_beginyr)
		return;

	// Average monthly active layer depth
	double mald = min(soil_depth, max(0.0, patch.soilcolumn.mthaw[date.month])) * M_PER_MM;

	if (patch.soilcolumn.savesomcent) {
		patch.soilcolumn.mald_mean[date.month] += mald;
	}

	// Low carbon in soil -> no vertical mixing
	if (date.year <= patch.soilcolumn.solvesomcent_endyr || patch.soilcolumn.ccont() < 1.0e-3)
		return;

	// Store pool sizes before mixing
	std::vector<double> cpools(3);
	std::vector<double> npools(3);

	for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
		cpools[0] += patch.soilcolumn.soillayers[lyr].sompool[SOILMICRO].cmass;
		npools[0] += patch.soilcolumn.soillayers[lyr].sompool[SOILMICRO].nmass;
		cpools[1] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;
		npools[1] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass;
		cpools[2] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;
		npools[2] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass;
	}

	// Timestep in Crank-Nicholson method 
	// Diffusion is on a yearly timestep and as vertical mixing is done monthly
	// the timestep will be 1 / 12
	double timestep = 1.0 / 12.0;

	patch.soilcolumn.vertical_mixing(SOILMICRO, mald, timestep);
	patch.soilcolumn.vertical_mixing(SLOWSOM, mald, timestep);
	patch.soilcolumn.vertical_mixing(PASSIVESOM, mald, timestep);


	// Check pool sizes after mixing
	std::vector<double> cpools_after(3);
	std::vector<double> npools_after(3);

	for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
		cpools_after[0] += patch.soilcolumn.soillayers[lyr].sompool[SOILMICRO].cmass;
		npools_after[0] += patch.soilcolumn.soillayers[lyr].sompool[SOILMICRO].nmass;
		cpools_after[1] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;
		npools_after[1] += patch.soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass;
		cpools_after[2] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;
		npools_after[2] += patch.soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass;
	}

	// Might be small rounding errors
	// Add to/subtract from top pool
	patch.soilcolumn.soillayers[0].sompool[SOILMICRO].cmass += cpools[0] - cpools_after[0];
	patch.soilcolumn.soillayers[0].sompool[SOILMICRO].nmass += npools[0] - npools_after[0];
	patch.soilcolumn.soillayers[0].sompool[SLOWSOM].cmass += cpools[1] - cpools_after[1];
	patch.soilcolumn.soillayers[0].sompool[SLOWSOM].nmass += npools[1] - npools_after[1];
	patch.soilcolumn.soillayers[0].sompool[PASSIVESOM].cmass += cpools[2] - cpools_after[2];
	patch.soilcolumn.soillayers[0].sompool[PASSIVESOM].nmass += npools[2] - npools_after[2];
}

/// Add litter to equilsom()
void add_litter(Soillayer& soillayer, int timestep, int pool, double init_cpool) {

	// If LUC have occurred in between solvesomcent_beginyr and solvesomcent_endyr then array might be too smalll
	if (timestep >= (int)soillayer.solvesom.size())
		timestep = timestep % (int)soillayer.solvesom.size();

	// How big the litter pool can grow comparred to the initial size (20%)
	const double litter_max_size_frac = 5.0;
	double frac_added = 1.0;
	double clitter_to_add = soillayer.solvesom[timestep].get_clitter(pool);

	if (soillayer.sompool[pool].cmass + clitter_to_add > litter_max_size_frac * init_cpool)
		frac_added = (litter_max_size_frac * init_cpool - soillayer.sompool[pool].cmass) / clitter_to_add;

	// Add to litter pools
	soillayer.sompool[pool].cmass += frac_added * clitter_to_add;
	soillayer.sompool[pool].nmass += frac_added * soillayer.solvesom[timestep].get_nlitter(pool);
}

/// Iteratively solving differential flux equations for century SOM pools
/** Iteratively solving differential flux equations for century SOM pools
 *  assuming annual litter inputs, nitrogen uptake and leaching is close
 *  to long term equilibrium
 */
void equilsom(Patch& patch) {

	if (!(date.year == patch.soilcolumn.solvesomcent_endyr && date.islastmonth && date.islastday))
		return;	

	Soilcolumn& soilcolumn = patch.soilcolumn;
	const Climate& climate = soilcolumn.patch.get_climate();
	const Gridcell& gridcell = patch.stand.get_gridcell();

	std::vector<double> save_NH4_avail(soilcolumn.nsoillayers);
	std::vector<double> save_NO3_avail(soilcolumn.nsoillayers);
	std::vector< std::vector<double> > save_litter_cmass(soilcolumn.nsoillayers, std::vector<double>(NSOMPOOL));
	std::vector< std::vector<double> > save_litter_nmass(soilcolumn.nsoillayers, std::vector<double>(NSOMPOOL));

	// Number of years with mean input data
	int EQUILSOM_YEARS = soilcolumn.solvesomcent_endyr - soilcolumn.solvesomcent_beginyr + 1;

	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		// Save soil mineral nitrogen status
		save_NH4_avail[lyr] = soillayer.nmass_avail(NH4);
		save_NO3_avail[lyr] = soillayer.nmass_avail(NO3);

		// Save C and N soil status
		for (int p = 0; p < NSOMPOOL; p++) {
			save_litter_cmass[lyr][p] = soillayer.sompool[p].cmass;
			save_litter_nmass[lyr][p] = soillayer.sompool[p].nmass;
		}

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

			// Monthly average organic carbon and nitrogen moved within soil column
			soillayer.morgmove_mean[m] /= EQUILSOM_YEARS;

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
	// Maximum number of iterations the trend needs to flip for the multiplier to decrease
	const int FLIP_MAX = 2;
	// Maximum multiplier
	const double MAX_MULTI = pow(2, 20);
	// Minimum multiplier
	const double MIN_MULTI = 0.2;
	// Multiplier modifier
	const double MULTI_MOD = 2.0;
	// Maximum tries of iterations
	const int MAX_TRIES = 10000;;
	// Number of trends that doesn't meet the criteria
	int bad_trend = 1;
	// Iteration tries
	int tries = 0;
	// SOM pool N:C ratio
	double ntoc;
	// C mass increase
	double cmass_inc;

	// Soil layers total C mass fractions
	std::vector<double> ccontfrac(soilcolumn.nsoillayers);
	// Initial SOM pool size per layer
	std::vector<double> init_slow(soilcolumn.nsoillayers);
	std::vector<double> init_passive(soilcolumn.nsoillayers);
	// Maximum C change for SOM pools per layer
	std::vector<double> max_change_slow(soilcolumn.nsoillayers);
	std::vector<double> max_change_passive(soilcolumn.nsoillayers);
	// Multiplier for SOM pools per layer
	std::vector<double> multi_slow(soilcolumn.nsoillayers);
	std::vector<double> multi_passive(soilcolumn.nsoillayers);
	// Number of trend sign flipping for SOM pools per layer
	std::vector<int> flip_slow(soilcolumn.nsoillayers);
	std::vector<int> flip_passive(soilcolumn.nsoillayers);
	// Current sign of trend for SOM pools per layer
	std::vector<int> sign_slow_trend(soilcolumn.nsoillayers);
	std::vector<int> sign_passive_trend(soilcolumn.nsoillayers);
	// Number of iterations with same trend sign for SOM pools per layer
	std::vector<int> same_trend_slow(soilcolumn.nsoillayers);
	std::vector<int> same_trend_passive(soilcolumn.nsoillayers);
	
	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {
		multi_slow[lyr] = MAX_MULTI;
		multi_passive[lyr] = MAX_MULTI;
		sign_slow_trend[lyr] = 1;
		sign_passive_trend[lyr] = 1;
		max_change_slow[lyr] = 0.5 * soilcolumn.soillayers[lyr].depth * M_PER_MM;
		max_change_passive[lyr] = max_change_slow[lyr];
	}

	// Spinup SOM pools until trends of slow and passive SOM pools 
	// are within limits or maximum tries have been done
	while (bad_trend && tries < MAX_TRIES) {

		// Update soil layers total C mass fractions
		ccontfrac = soilcolumn.get_ccontfrac();

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

					// Initial C for slow and passive SOM pools in this iteration
					if (yr == 0 && m == 0) {
						init_slow[lyr] = soillayer.sompool[SLOWSOM].cmass;
						init_passive[lyr] = soillayer.sompool[PASSIVESOM].cmass;
					}

					// Transfer monthly mean litter (added a maximum litter pool size so a runaway pool doesn't occur)
					int timestep = yr * 12 + m;
					if (lyr == 0) {
						add_litter(soillayer, timestep, SURFSTRUCT, save_litter_cmass[lyr][SURFSTRUCT]);
						add_litter(soillayer, timestep, SURFMETA, save_litter_cmass[lyr][SURFMETA]);
						add_litter(soillayer, timestep, SURFFWD, save_litter_cmass[lyr][SURFFWD]);
						add_litter(soillayer, timestep, SURFCWD, save_litter_cmass[lyr][SURFCWD]);
					}
					add_litter(soillayer, timestep, SOILSTRUCT, save_litter_cmass[lyr][SOILSTRUCT]);
					add_litter(soillayer, timestep, SOILMETA, save_litter_cmass[lyr][SOILMETA]);
					add_litter(soillayer, timestep, SOILFWD, save_litter_cmass[lyr][SOILFWD]);
					add_litter(soillayer, timestep, SOILCWD, save_litter_cmass[lyr][SOILCWD]);

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
					soillayer.nmass_inc(soilcolumn.anfix_mean / 12.0 * ccontfrac[lyr], NH4);

					// Monthly decomposition and fluxes between SOM pools

					// Set this months decay rates
					for (int p = 0; p < NSOMPOOL; p++) {
						soillayer.sompool[p].fracremain = soillayer.sompool[p].mfracremain_mean[m];
					}

					// Set this months organic matter leaching and move fraction
					soillayer.orgleachfrac = soillayer.morgleach_mean[m];
					soillayer.orgmovefrac = soillayer.morgmove_mean[m];
				}

				// Monthly decomposition and fluxes between SOM pools
				// and nitrogen flux from soil
				somfluxes(patch, true, false);

				if (!vertmix && ifverticalmixing) {
					soilcolumn.vertical_mixing(SLOWSOM, soilcolumn.mald_mean[m], (double)VMTIMESTEP / 12.0);
					soilcolumn.vertical_mixing(PASSIVESOM, soilcolumn.mald_mean[m], (double)VMTIMESTEP / 12.0);
				}
			}
		}
		
		// Reset number of trends that doesn't meet the criteria
		bad_trend = 0;

		// Increase number of iterations
		tries++;

		//dprintf("Tries %d CCONT ", tries);

		// Loop through layers and try and to solve equilibrium by appling the trend of pools
		// by several order of magnitude
		for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

			// After 100 iterations, then maximum increase is allowed to be 10% of pool size
			if (tries == 100) {
				max_change_slow[lyr] = soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass * 0.1;
				max_change_passive[lyr] = soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass * 0.1;
			}

			// Both slow and passive pool needs to be of some size
			if (!negligible(soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass, -15) &&
				!negligible(soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass, -15)) {

				// Trend over NYEARS years
				double trend_slow = soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass - init_slow[lyr];
				double trend_passive = soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass - init_passive[lyr];
				
				// If trend changes sign a certain number of times then reduce multiplier
				// If trend has the same sign a certain number of times then increase multiplier

				// SLOW SOM

				// Trend changed sign
				if (trend_slow * sign_slow_trend[lyr] < 0.0) {

					// Add to number of times trend has changed sign in a raw
					flip_slow[lyr]++;

					// Save trend sign and set same trend to zero
					sign_slow_trend[lyr] *= -1;

					// Flipped sign maximum times then reduce multiplier
					if (flip_slow[lyr] == FLIP_MAX && multi_slow[lyr] > MIN_MULTI) {
						multi_slow[lyr] /= MULTI_MOD;
						flip_slow[lyr] = 0;
					}

					// Flipping so no same trend
					same_trend_slow[lyr] = 0;
				}
				// Same sign trend
				else {

					// Add to number of times same trend has happen in a raw
					same_trend_slow[lyr]++;

					// Reached maximum same trend then increase multiplier
					if (same_trend_slow[lyr] > SAME_TREND_MAX && multi_slow[lyr] < MAX_MULTI) {
						multi_slow[lyr] *= MULTI_MOD;
						same_trend_slow[lyr] = 0;
					}

					// Same trend so no flipping
					flip_slow[lyr] = 0;
				}

				// N:C ratio of SOM pool
				ntoc = soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass / soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;

				// C change
				cmass_inc = min(max_change_slow[lyr], max(-max_change_slow[lyr], max(-soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass * 0.1, trend_slow * multi_slow[lyr])));

				// Add C change
				soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass += cmass_inc;

				// Update N pool
				soilcolumn.soillayers[lyr].sompool[SLOWSOM].nmass = ntoc * soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass;


				// PASSIVE SOM

				// Trend changed sign
				if (trend_passive * sign_passive_trend[lyr] < 0.0) {

					// Add to number of times trend has changed sign in a raw
					flip_passive[lyr]++;

					// Save trend sign and set same trend to zero
					sign_passive_trend[lyr] *= -1;

					// Flipped sign maximum times then reduce multiplier
					if (flip_passive[lyr] == FLIP_MAX && multi_passive[lyr] > MIN_MULTI) {
						multi_passive[lyr] /= MULTI_MOD;
						flip_passive[lyr] = 0;
					}

					// Flipping so no same trend
					same_trend_passive[lyr] = 0;
				}
				// Same sign trend
				else {
					// Add to number of times same trend has happen in a raw
					same_trend_passive[lyr]++;

					// Reached maximum same trend then increase multiplier
					if (same_trend_passive[lyr] > SAME_TREND_MAX && multi_passive[lyr] < MAX_MULTI) {
						multi_passive[lyr] *= MULTI_MOD;
						same_trend_passive[lyr] = 0;
					}

					// Same trend so no flipping
					flip_passive[lyr] = 0;
				}

				// N:C ratio of SOM pool
				ntoc = soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass / soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;

				// C change
				cmass_inc = min(max_change_passive[lyr], max(-max_change_passive[lyr], max(-soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass * 0.1, trend_passive * multi_passive[lyr])));

				// Add C change
				soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass += cmass_inc;
				
				// Update N pool
				soilcolumn.soillayers[lyr].sompool[PASSIVESOM].nmass = ntoc * soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass;


				// Calculte how many trends that are too big
				if (iftwolayersoil || (lyr > 1 && fabs(trend_passive) / soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass > 1.0e-6))
					bad_trend++;

				//dprintf("%.5f %.5f %.5f %.10f %.10f %.10f %.0f ", (soilcolumn.soillayers[lyr].ccont() / (soilcolumn.soillayers[lyr].depth * M_PER_MM)) / maxSOCdensity, soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass, soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass, trend_slow / soilcolumn.soillayers[lyr].sompool[SLOWSOM].cmass, trend_passive / soilcolumn.soillayers[lyr].sompool[PASSIVESOM].cmass, cmass_inc, multi_passive[lyr]);
			}
		}

		//dprintf("Bad %d \n", bad_trend);
	}

	// Reset mineral nitrogen status
	for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {

		Soillayer& soillayer = soilcolumn.soillayers[lyr];

		soillayer.NH4_mass = save_NH4_avail[lyr];
		soillayer.NO3_mass = save_NO3_avail[lyr];

		// Reset C and N litter status
		for (int p = 0; p < NSOMPOOL; p++) {
			if (p != SURFMICRO && p != SURFHUMUS && p != SOILMICRO && p != SLOWSOM && p != PASSIVESOM) {
				soillayer.sompool[p].cmass = save_litter_cmass[lyr][p];
				soillayer.sompool[p].nmass = save_litter_nmass[lyr][p];
			}
		}

		// Reset variables for next possible equilsom()
		for (int m = 0; m < 12; m++) {

			for (int p = 0; p < NSOMPOOL; p++) {
				soillayer.sompool[p].mfracremain_mean[m] = 0.0;
			}

			soillayer.fNH4uptake_mean[m] = 0.0;
			soillayer.fNO3uptake_mean[m] = 0.0;
			soillayer.morgleach_mean[m] = 0.0;
			soillayer.morgmove_mean[m] = 0.0;
			soillayer.mminleach_mean[m] = 0.0;
		}

		soillayer.solvesom.clear();
		std::vector<LitterSolveSOM>().swap(soillayer.solvesom); // clear array memory
	}

	// Reset soil column annual N fixation
	soilcolumn.anfix_mean = 0.0;
}


/// SOM CENTURY DYNAMICS
/** To be called each simulation day for each modelled patch, following update
 *  of soil temperature and soil water.
 *  Transfers litter, performes nitrogen uptake and addition, leaching and decomposition.
 */
void som_dynamics_century(Patch& patch, Climate& climate, bool tillage) {

	// Transfer litter to SOM pools
	transfer_litter(patch);

	// Daily nitrogen addition to the soil
	soilnadd(patch);

	// Daily mineral and organic nitrogen leaching
	leaching(patch);

	// Daily or monthly decomposition and fluxes between SOM pools
	somfluxes(patch, false, tillage);

	// Nitrogen transformation in soil
	ntransform(patch, climate);

	// Vertical mixing(diffusion) of soil organic matter
	vertical_mixing(patch);

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

