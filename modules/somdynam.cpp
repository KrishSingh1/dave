///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.cpp
/// \brief Soil organic matter dynamics
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
#include "somdynam.h"

#include "driver.h"
#include "growth.h"
#include "vegdynam.h"


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

// Corresponds to minimum soil available N where SOM C:N ratio reach
// their minimum (Parton et al 1993, Fig. 4)
static const double nmin_balance_max = 0.002;	

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

void decayrates(double wcont,double gtemp_soil,double& k_soilfast,double& k_soilslow,
	double& fr_litter,double& fr_soilfast,double& fr_soilslow) {

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

	moist_response=0.25+0.75*wcont;

	// Calculate litter and SOM fractions remaining following today's decomposition
	// (Sitch et al 2000 Eqn 71) adjusting exponential decay constants by moisture and
	// temperature responses and converting from annual to daily basis
	// NB: Temperature response (gtemp; Lloyd & Taylor 1994) set by framework

	k_soilfast=k_soilfast10*gtemp_soil*moist_response/365.0;
	k_soilslow=k_soilslow10*gtemp_soil*moist_response/365.0;

	fr_litter=exp(-k_litter10*gtemp_soil*moist_response/365.0);
	fr_soilfast=exp(-k_soilfast);
	fr_soilslow=exp(-k_soilslow);
}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Should be called by framework on last day of simulation year, following call to
// som_dynamics, once annual litter production and vegetation PFT composition are close
// to their long term equilibrium (typically 500-1000 simulation years).
// NB: should be called ONCE ONLY during simulation for a particular grid cell

void equilsom(Soil& soil) {

	// DESCRIPTION
	// Analytically solves differential flux equations for fast and slow SOM pools
	// assuming annual litter inputs close to long term equilibrium

	// INPUT PARAMETER (class defined in framework header file)
	// soil = current soil status

	double nyear;
		// number of years over which decay constants and litter inputs averaged

	nyear=soil.soiltype.solvesom_end-soil.soiltype.solvesom_begin+1;

	soil.decomp_litter_mean/=nyear;
	soil.k_soilfast_mean/=nyear;
	soil.k_soilslow_mean/=nyear;

	soil.cpool_fast=(1.0-ATMFRAC)*FASTFRAC*soil.decomp_litter_mean/
		soil.k_soilfast_mean;
	soil.cpool_slow=(1.0-ATMFRAC)*(1.0-FASTFRAC)*soil.decomp_litter_mean/
		soil.k_soilslow_mean;
}


///////////////////////////////////////////////////////////////////////////////////////
// SOM DYNAMICS
// To be called each simulation day for each modelled area or patch, following update
// of soil temperature and soil water.

void som_dynamics_lpj(Patch& patch) {

	// DESCRIPTION
	// Calculation of soil decomposition and transfer of C between litter and soil
	// organic matter pools.
	//
	// NB: The global variable 'ifdailydecomp' determines whether soil decomposition
	// calculations are performed every day, or on the last day of each month, based on
	// average conditions for the month (the latter mode is much faster). Daily flux
	// values are never valid in monthly mode. If you require daily output, use daily
	// mode

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
	Soil& soil=patch.soil;

	if (ifdailydecomp) {

		// "DAILY" MODE

		// Calculate decay constants and rates given today's soil moisture and
		// temperature

		decayrates(soil.wcont[0],soil.gtemp,k_soilfast,k_soilslow,fr_litter,
			fr_soilfast,fr_soilslow);

		// From year soil.solvesom_begin, update running means for later solution
		// (at year soil.solvesom_end) of equilibrium SOM pool sizes

		if (date.year>=soil.soiltype.solvesom_begin) {
			soil.k_soilfast_mean+=k_soilfast;
			soil.k_soilslow_mean+=k_soilslow;
		}
	}
	else if (date.islastday) {

		// "MONTHLY" MODE (last day of month only)

		// Calculate respiration temperature response if not yet done for this month

		if (soil.last_mgtemp!=date.month) {
			respiration_temperature_response(soil.mtemp,soil.mgtemp);
			soil.last_mgtemp=date.month;
		}

		// Calculate decay constants and rates given monthly means

		decayrates(soil.mwcontupper,soil.mgtemp,k_soilfast,k_soilslow,fr_litter,
			fr_soilfast,fr_soilslow);

		// From year soil.solvesom_begin, update running means for later solution
		// (at year soil.solvesom_end) of equilibrium SOM pool sizes

		if (date.year>=soil.soiltype.solvesom_begin) {
			soil.k_soilfast_mean+=k_soilfast*(double)date.ndaymonth[date.month];
			soil.k_soilslow_mean+=k_soilslow*(double)date.ndaymonth[date.month];
		}

		// Convert fractional scalars from daily to monthly basis

		fr_litter=pow(fr_litter,date.ndaymonth[date.month]);
		fr_soilfast=pow(fr_soilfast,date.ndaymonth[date.month]);
		fr_soilslow=pow(fr_soilslow,date.ndaymonth[date.month]);
	}

	// DAILY AND MONTHLY MODES
	// Only on last day of month if monthly mode

	// Reduce litter and SOM pools, sum C flux to atmosphere from decomposition
	// and transfer correct proportions of litter decomposition to fast and slow
	// SOM pools

	if (ifdailydecomp || date.islastday) {

		// Reduce individual litter pools and calculate total litter decomposition
		// for today/this month

		decomp_litter=0.0;

		// Loop through PFTs

		for (p=0;p<npft;p++) {
			
			// For this PFT ...

			decomp_litter+=(patch.pft[p].litter_leaf+
				patch.pft[p].litter_root+
				patch.pft[p].litter_wood+
				patch.pft[p].litter_repr)*(1.0-fr_litter);

			patch.pft[p].litter_leaf*=fr_litter;
			patch.pft[p].litter_root*=fr_litter;
			patch.pft[p].litter_wood*=fr_litter;
			patch.pft[p].litter_repr*=fr_litter;
		}

		if (date.year>=soil.soiltype.solvesom_begin)
			soil.decomp_litter_mean+=decomp_litter;

		// Partition litter decomposition among fast and slow SOM pools
		// and flux to atmosphere

		// flux to atmosphere
		cflux=decomp_litter*ATMFRAC;

		// remaining decomposition - goes to ...
		decomp_litter-=cflux;

		// ... fast SOM pool ...
		soil.cpool_fast+=decomp_litter*FASTFRAC;

		// ... and slow SOM pool
		soil.cpool_slow+=decomp_litter*(1.0-FASTFRAC);

		// Increment C flux to atmosphere by SOM decomposition
		cflux+=soil.cpool_fast*(1.0-fr_soilfast)+soil.cpool_slow*(1.0-fr_soilslow);

		// Monthly C flux

		if (ifdailydecomp) {
			if (date.dayofmonth==0)
				patch.fluxes.mcflux_soil[date.month]=cflux;
			else
				patch.fluxes.mcflux_soil[date.month]+=cflux;
		}
		else
			patch.fluxes.mcflux_soil[date.month]=cflux;

		// Reduce SOM pools 
		soil.cpool_fast*=fr_soilfast;
		soil.cpool_slow*=fr_soilslow;

		// Updated daily and annual fluxes
		patch.fluxes.dcflux_soil[date.day]=cflux;
		patch.fluxes.acflux_soil+=cflux;

		// Solve SOM pool sizes at end of year given by soil.solvesom_end

		if (date.year==soil.soiltype.solvesom_end && date.islastmonth && date.islastday)
			equilsom(soil);
	}
}

/////////////////////////////////////////////////
// CENTURY SOM DYNAMICS


/// Estimates the daily mineral nitrogen available for plant uptake  
/** Instead of using the soil.nmass_avail at day==0 as in Parton et al 1993,
 *  nmass_avail is "updated" each day (soil.nmin_balance) depending on mineralization, 
 *  immobilization, deposition, fixation, leaching and plant uptake (no N limitation  
 *  and with last years growth C:N ratio). Then the ntoc ratios for SOM pools are able to be 
 *  calculated each day. 
 *  Mineral nitrogen leaching is also performed in this function.
 *  In future development denitrification and nitrfication will be called from this function 
 */
void est_nmin_balance(Patch& patch, Soil& soil, Climate& climate) {

	// N fixation (using last year as an estimate as it is calculated on last day of year)
	double dnfix;

	// If disturbance then no aaet -> can't use last years estimate as leaching might then exceed
	// available N in end of year
	if (patch.age)
		dnfix = soil.anfix / 365.0;
	else {
		if (nfix_b > 0.0)
			dnfix = nfix_b / 100000.0 / 365.0;
		else
			dnfix = 0.0;
	}

	// First day of year
	if (date.day == 0) {
		soil.aminleach = 0.0;
		soil.nmin_balance = soil.nmass_avail + climate.dndep[date.day] + dnfix;
	}
	else {
		// Update "daily" nmass available 
		soil.nmin_balance += soil.nmin_daily[date.day-1] - soil.nimmob_daily[date.day-1] + climate.dndep[date.day] + dnfix;
	}

	// Loop through individuals to determine N demand
	double N_demand = 0.0;
		
	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) 
	{
		Individual& indiv = vegetation.getobj();

		// For this individual ...

		double NPPp = indiv.assim - indiv.resp;

		if (NPPp > 0.0)
			N_demand += NPPp / indiv.cton_growth;		

		// reset values for next next days photosynthesis and respiration
		indiv.assim = 0.0;
		indiv.resp = 0.0;

		vegetation.nextobj();
	}

	double N_availability = max(0.0, soil.nmin_balance);
	double N_uptake = min(N_demand, N_availability);
	N_uptake = max(0.0, N_uptake);

	soil.nmin_balance -= N_uptake;

	// LEACHING OF SOIL MINERAL N
	// Allowed on days with residual N following estimated vegetation N uptake
	// in proportion to percolation following Parton et al. 1993 eqn 13
	if (soil.nmin_balance > 0.0) {
		double leaching = soil.nmin_balance * soil.minleachfrac_daily[date.day];
		soil.nmin_balance -= leaching;
		soil.aminleach += leaching;
		soil.sompool[LEACHED].nmass += leaching;
	}
}

/// Decreases decay rates to keep the daily N balance in the soil  
/** Function goes through the different SOM pools in a specific order
 *  depending on the order they feed into eachother.
 */
void neg_mineralization(double decay_reduction[NSOMPOOL], double net_min[NSOMPOOL], int start, int end, double daily_nmass, double nmin_balance) {

	int p;
	double tot_neg_min;
	int neg_min[NSOMPOOL] = {0};	// Keeping track on which pools that are negative
	int order[NSOMPOOL] = {0,5,1,7,6,4,3,8,2,9,10}; // Order to go through soil pools when decomposition is N limited
	// SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFCWD,SOILMETA,SLOWSOM,PASSIVESOM,LEACHED,NSOMPOOL

	for (p=0; p<NSOMPOOL; p++)
		neg_min[p] = 0;

	double decay_red = 0.0;	// decay reduction

	tot_neg_min = 0.0;
	for (p=start; p<end; p++) {
		if (net_min[p] < 0.0) {
			tot_neg_min += net_min[p];
			neg_min[p] = 1;
		}
	}
	if (tot_neg_min < daily_nmass)
		decay_red = 1.0 - (tot_neg_min - (daily_nmass + nmin_balance)) / tot_neg_min;
	else
		decay_red = 1.0;
		
	for (p=start;p<end;p++)
		decay_reduction[order[p]] = decay_red * neg_min[p];
}

/// Set N:C ratios for SOM pools  
/** Set N:C ratios for slow, passive, humus and soil microbial pools
 *  based on mineral N pool or litter N fraction (Parton et al 1993, Fig 4)
 */
void setntoc(Soil& soil, double fac, pooltype pool, double cton_max, double cton_min,
	double fmin, double fmax) {

	if (fac <= fmin) 
		soil.sompool[pool].ntoc = 1.0 / cton_max;
	else if (fac >= fmax) 
		soil.sompool[pool].ntoc = 1.0 / cton_min;
	else {
		soil.sompool[pool].ntoc = 1.0 / (cton_min + (cton_max - cton_min) *
			(fmax - fac) / (fmax - fmin));
	}
}

/// Calculates CENTURY instantaneous decay rates   
/** Calculates CENTURY instantaneous decay rates given soil temperature,   
 *  water content of upper soil layer
 */
void decayrates(Soil& soil, double temp_soil, double wcont_soil) {

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris decay	
	const double K_MAX[] = {9.2e-3, 1.8e-2, 4.0e-2, 4.8e-4, 2.6e-2, 3.7e-2, 2.2e-3, 6.8e-2, 1.7e-3, 6.9e-5};
	// SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFCWD,SOILMETA,SLOWSOM,PASSIVESOM

	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	double texture_mod = 1.0 - 0.75 * (soil.soiltype.clay_frac + soil.soiltype.silt_frac);

	double temp_mod;
	double moist_mod;
	double wfps;
	double k;
	int p;

	// Calculate decomposition temperature modifier (in range 0-1)
	// [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
	// effect of soil temperature, Parton et al 1993, Fig 2)

	if (temp_soil > 0.0)
		temp_mod = max(0.0,
			0.0326 + 0.00351 * pow(temp_soil, 1.652) - pow(temp_soil / 41.748, 7.19));
	else
		temp_mod = 0.0;

	// Calculate decomposition moisture modifier (in range 0-1)
	// Friend et al 1997, Eqn 53
	// (Parton et al 1993, Fig 2)

	// Water Filled Pore Spaces (wfps)
	// water holding capacity at wilting point (wp) and saturation capacity (wsats) 
	// is calculated with the help of Cosby et al 1984;
	wfps = (wcont_soil * soil.soiltype.awc[0] + soil.soiltype.wp[0]) * 100.0 / soil.soiltype.wsats[0];			

	if (wfps < 60.0)
		moist_mod = exp((wfps - 60.0) * (wfps - 60.0) / -800.0);
	else
		moist_mod = 0.000371 * wfps * wfps - 0.0748 * wfps + 4.13;

	for (p=0;p<NSOMPOOL-1;p++) {

		// Calculate decay constant (annual basis)
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		k = K_MAX[p] * temp_mod * moist_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2

		if (p == SURFSTRUCT || p == SOILSTRUCT)// || p == SURFCWD)
			k *= exp(-3.0 * soil.sompool[p].ligcfrac);
		else if (p == SOILMICRO)
			k *= texture_mod;

		// Calculate fraction of C pool remaining after today's decomposition

		soil.sompool[p].frc = exp(-k);	
	}
}

/// Transfers specified fraction (frac) of today's decomposition   
/** Transfers specified fraction (frac) of today's decomposition in donor pool type   
 *  to receiver pool, transferring fraction respfrac of this to the accumulated CO2
 *  flux respsum (representing total microbial respiration today)
 */
void transferdecomp(Soil& soil, pooltype donor, pooltype receiver,
	double frac, double respfrac, double& respsum, double& nmin_actual,
	double& nimmob, double& net_min) {

	// decrement in donor C pool and N pools
	double cdec = soil.sompool[donor].cdec * frac;
	double ndec = soil.sompool[donor].ndec * frac;

	// associated N increment in receiver pool (Friend et al 1997, Eqn 49)
	double ninc = cdec * (1.0 - respfrac) * soil.sompool[receiver].ntoc;

	// if increase in receiver N greater than decrease in donor N,
	// balance must be immobilisation from mineral N pool
	// otherwise balance is N mineralisation
	if (ninc > ndec) {
		nimmob += ninc - ndec;
		net_min += ndec - ninc;
	}
	else {
		nmin_actual += ndec - ninc;
		net_min += ndec - ninc;
	}

	// "Transfer" C and N to receiver
	soil.sompool[receiver].delta_cmass += cdec * (1.0 - respfrac);
	soil.sompool[receiver].delta_nmass += ninc;

	// Transfer microbial respiration
	respsum += cdec * respfrac;
}

/// Fluxes between the CENTURY pools, and CO2 release to the atmosphere   
/** Daily or monthly fluxes between the ten CENTURY pools, and CO2 release to the atmosphere   
 *  Parton et al 1993, Fig 1; Comins & McMurtrie 1993, Appendix A
 */
void somfluxes(Patch& patch, Soil& soil,Fluxes& fluxes) {	

	int p, d;
	double csp, csa, respfrac, cap;
	double respsum = 0.0;
	double leachsum_cmass, leachsum_nmass;
	double nmin_actual = 0.0;	// actual (not net) N mineralisation
	double nimmob = 0.0;		// N immobilisation

	const double EPS = 1.0e-16;

	// Set N:C ratios for humus, soil microbial, passive and slow pool based on estimated mineral N pool
	// (Parton et al 1993, Fig 4)

	est_nmin_balance(patch, soil, patch.stand.gridcell.climate);

	if (soil.nmin_balance < -EPS && date.year >= freenyears)
		dprintf("Year %d Day %d WRONG nmin_balance %g \n",date.year, date.day, soil.nmin_balance);

	// ForCent values
	setntoc(soil, soil.nmin_balance, SLOWSOM, 30.0, 15.0, 0.0, nmin_balance_max);
		
	setntoc(soil, soil.nmin_balance, PASSIVESOM, 10.0, 6.0, 0.0, nmin_balance_max);

	setntoc(soil, soil.nmin_balance, SOILMICRO, 15.0, 6.0, 0.0, nmin_balance_max);

	setntoc(soil, soil.nmin_balance, SURFHUMUS, 30.0, 15.0, 0.0, nmin_balance_max);

	if (ifdailydecomp || ifnlim) {

		// DAILY MODE

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no N limitation)
		decayrates(soil, soil.temp, soil.wcont[0]); 

	}
	else if (date.islastday) {

		// MONTHLY MODE (last day of month only)

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no N limitation)
		decayrates(soil, soil.temp, soil.mwcontupper);

		// Convert fractional scalars from daily to monthly basis

		for (p=0; p<NSOMPOOL; p++) {
			soil.sompool[p].frc = pow(soil.sompool[p].frc, date.ndaymonth[date.month]);
		}
	}
	else return;  // Nothing to do if monthly mode but not last day of month

	// Calculate decomposition in all pools assuming these decay rates

	// Save delta C and N mass

	bool net_mineralization = false;
	int times = 0;
	double decay_reduction[NSOMPOOL] = {0.0};

	// If mineralization together with soil available N is negative then pools decay rates are decreased 
	// The SOM system gives five try to get a positive result
	while(!net_mineralization && times < 5) {

		respsum = 0.0;
		nmin_actual = 0.0;
		nimmob = 0.0;
		leachsum_cmass = 0.0;
		leachsum_nmass = 0.0;

		// Calculate decomposition in all pools assuming these decay rates
		for (p=0;p<NSOMPOOL;p++) {
			soil.sompool[p].cdec = soil.sompool[p].cmass * (1.0 - soil.sompool[p].frc) * (1.0 - decay_reduction[p]);
			soil.sompool[p].ndec = soil.sompool[p].nmass * (1.0 - soil.sompool[p].frc) * (1.0 - decay_reduction[p]);
			
			soil.sompool[p].delta_cmass = 0.0;
			soil.sompool[p].delta_nmass = 0.0;
			soil.sompool[p].delta_cmass -= soil.sompool[p].cdec;
			soil.sompool[p].delta_nmass -= soil.sompool[p].ndec;
		}

		double net_min[NSOMPOOL] = {0};
		
		// Partition potential decomposition among receiver pools

		// Donor pool SURFACE STRUCTURAL

		transferdecomp(soil, SURFSTRUCT, SURFMICRO, 1.0 - soil.sompool[SURFSTRUCT].ligcfrac,
			0.6, respsum, nmin_actual, nimmob, net_min[0]);	

		transferdecomp(soil, SURFSTRUCT, SURFHUMUS, soil.sompool[SURFSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, nimmob, net_min[0]);

		// Donor pool SURFACE METABOLIC

		transferdecomp(soil, SURFMETA, SURFMICRO, 1.0, 0.6, respsum, nmin_actual, nimmob, net_min[1]);

		// Donor pool SOIL STRUCTURAL

		transferdecomp(soil, SOILSTRUCT, SOILMICRO, 1.0 - soil.sompool[SOILSTRUCT].ligcfrac,
			0.55, respsum, nmin_actual, nimmob, net_min[2]);

		transferdecomp(soil, SOILSTRUCT, SLOWSOM, soil.sompool[SOILSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, nimmob, net_min[2]);

		// Donor pool SOIL METABOLIC

		transferdecomp(soil, SOILMETA, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, nimmob, net_min[3]);

		// Donor pool SURFACE COARSE WOODY DEBRIS

		transferdecomp(soil, SURFCWD, SURFMICRO, 1.0 - soil.sompool[SURFCWD].ligcfrac,	
			0.76, respsum, nmin_actual, nimmob, net_min[4]);

		transferdecomp(soil, SURFCWD, SURFHUMUS, soil.sompool[SURFCWD].ligcfrac, 0.76,
			respsum, nmin_actual, nimmob, net_min[4]);
	
		// Donor pool SURFACE MICROBE

		transferdecomp(soil, SURFMICRO, SURFHUMUS, 1.0, 0.6, respsum, nmin_actual, nimmob, net_min[5]);

		// Donor pool SURFACE HUMUS

		transferdecomp(soil, SURFHUMUS, SLOWSOM, 1.0, 0.6, respsum, nmin_actual, nimmob, net_min[6]);

		// Donor pool SLOW SOM
	
		// First work out partitioning coefficients (Fig 1, Parton et al 1993)

		csp = 0.003 - 0.009 * soil.soiltype.clay_frac;
		respfrac = 0.55;
		csa= 1.0 - csp - respfrac;

		transferdecomp(soil, SLOWSOM, SOILMICRO, csa, 0.0, respsum, nmin_actual, nimmob, net_min[7]);

		transferdecomp(soil, SLOWSOM, PASSIVESOM, csp, 0.0, respsum, nmin_actual, nimmob, net_min[7]);

		// Account for respiration flux
		// N associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SLOWSOM].cdec;

		if(!negligible(soil.sompool[SLOWSOM].cmass))
			nmin_actual += respfrac * soil.sompool[SLOWSOM].cdec * soil.sompool[SLOWSOM].nmass / soil.sompool[SLOWSOM].cmass;	

		// Donor pool SOIL MICROBE

		// Fraction lost to microbial respiration (F_t, Parton et al 1993 Eqn 7)
		respfrac = max(0.0, 0.85 - 0.68 * (soil.soiltype.clay_frac + soil.soiltype.silt_frac));

		// Fraction entering passive SOM pool (Parton et al 1993, Eqn 9)
		cap = 0.003 + 0.032 * soil.soiltype.clay_frac;

		transferdecomp(soil, SOILMICRO, PASSIVESOM, cap, 0.0, respsum, nmin_actual, nimmob, net_min[8]);

		// Fraction entering slow SOM pool
		csp = 1.0 - respfrac - soil.orgleachfrac_daily[date.day] - cap;

		transferdecomp(soil, SOILMICRO, SLOWSOM, csp, 0.0, respsum, nmin_actual, nimmob, net_min[8]);

		// Account for respiration flux
		// N associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SOILMICRO].cdec;

		// Account for organic carbon leaching loss
		leachsum_cmass = soil.orgleachfrac_daily[date.day] * soil.sompool[SOILMICRO].cdec;
		
		if(!negligible(soil.sompool[SOILMICRO].cmass)) {
			nmin_actual += respfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
			
			// Account for organic nitrogen leaching loss
			leachsum_nmass = soil.orgleachfrac_daily[date.day] * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
		}

		// Donor pool PASSIVE SOM

		transferdecomp(soil, PASSIVESOM, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, nimmob, net_min[9]);

		// Estimate daily soil mineral N pool after decomposition
		// (negative value = immobilisation) 

		double daily_nmass = nmin_actual - nimmob;

		if ((daily_nmass + soil.nmin_balance + EPS >= 0.0) || date.year < freenyears) {

			net_mineralization = true;
		}
		else {

			// Immobilization larger than soil available N -> decrease decay rates
			if (times == 0)
				neg_mineralization(decay_reduction, net_min, 0, 5, daily_nmass, soil.nmin_balance);
			else if (times == 1)
				neg_mineralization(decay_reduction, net_min, 5, 6, daily_nmass, soil.nmin_balance);
			else if (times == 2)
				neg_mineralization(decay_reduction, net_min, 6, 7, daily_nmass, soil.nmin_balance);
			else if (times == 3)
				neg_mineralization(decay_reduction, net_min, 7, 10, daily_nmass, soil.nmin_balance);

			net_mineralization = false;
		}

		times++;
	}

	// Update pool sizes

	for (p=0;p<NSOMPOOL-1;p++) {
		soil.sompool[p].cmass += soil.sompool[p].delta_cmass;
		soil.sompool[p].nmass += soil.sompool[p].delta_nmass;
	}

	// Transfer respiration sum to fluxes

	fluxes.dcflux_soil[date.day] = respsum;
	fluxes.mcflux_soil[date.month] += respsum;
	fluxes.acflux_soil += respsum;

	// Transfer organic leaching to pool

	soil.sompool[LEACHED].cmass += leachsum_cmass;
	soil.sompool[LEACHED].nmass += leachsum_nmass;

	// Sum annual organic nitrogen leaching
	if (date.day == 0)
		soil.aorgleach = 0.0;

	soil.aorgleach += leachsum_nmass;

	// Store daily mineralisation and immobilisation to permit calculation of daily
	// mineral nitrogen balance at end of year

	if (ifdailydecomp) {
		soil.nmin_daily[date.day] = nmin_actual;
		soil.nimmob_daily[date.day] = nimmob;
	}
	else { // monthly mode - distribute current values evenly through the current month
		for (d=0; d<date.ndaymonth[date.month]; d++) {
			soil.nmin_daily[date.day-d] = nmin_actual / (double)date.ndaymonth[date.month];
			soil.nimmob_daily[date.day-d] = nimmob / (double)date.ndaymonth[date.month];
		}
	}
}

/// Transfers litter from this year's growth, mortality and fire   
/** Call annually after growth, mortality and fire to transfer this year's litter   
 *  from vegetation to soil litter pools
 */
void transfer_litter(Patch& patch, Soil& soil) {

	// Leaf, root and wood litter lignin fractions
	// Leaf and root fractions: Comins & McMurtrie 1993; Friend et al 1997
	// Not sure of wood fraction
	const double LIGCFRAC_LEAF = 0.2;
	const double LIGCFRAC_ROOT = 0.16;
	const double LIGCFRAC_WOOD = 0.3;

	double lton;	// Leaf litter ligning to N ratio
	double fm;
	double ligcmass_old, ligcmass_new;
	double litter_leaf_n;

	double litter_nmass = 0.0;	
	double litter_cmass = 0.0;	

	// Fire
	double litterme[3];
	litterme[0] = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].litterme;
	litterme[1] = soil.sompool[SURFMETA].cmass * soil.sompool[SURFMETA].litterme;
	litterme[2] = soil.sompool[SURFCWD].cmass * soil.sompool[SURFCWD].litterme;

	double fireresist[3];
	fireresist[0] = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].fireresist;
	fireresist[1] = soil.sompool[SURFMETA].cmass * soil.sompool[SURFMETA].fireresist;
	fireresist[2] = soil.sompool[SURFCWD].cmass * soil.sompool[SURFCWD].fireresist;

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft=patch.pft.getobj();

		// Calculate total litter C and N mass for set N:C ratio of surface microbial pool
		litter_nmass += (pft.nmass_litter_leaf + pft.nmass_litter_root + pft.nmass_litter_wood);
		litter_cmass += (pft.litter_leaf + pft.litter_root + pft.litter_wood);

		// LEAF

		// Calculate inputs to surface structural and metabolic litter

		litter_leaf_n = pft.nmass_litter_leaf; 

		// Leaf litter lignin:N ratio
		if (!negligible(litter_leaf_n))
			lton = max(0.0, LIGCFRAC_LEAF * pft.litter_leaf / litter_leaf_n);
		else
			lton = max(0.0, LIGCFRAC_LEAF * pft.pft.cton_leaf_avr / (1.0 - nrelocfrac));

		// Metabolic litter fraction for leaf litter (Fm, Parton et al 1993, Eqn 1:
		// NB: incorrect/out-of-date intercept and slope given in Eqn 1; values used in
		// code of CENTURY 4.0 used instead)
		fm = max(0.0, 0.85 - lton * 0.013);

		ligcmass_old = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].ligcfrac;

		if (fm < 0.0 || fm > 1.0) 
			dprintf("Year %d LEAF fm %g pft %s\n", date.year, fm, (char*)pft.pft.name);

		// Add to pools
		soil.sompool[SURFSTRUCT].cmass += pft.litter_leaf * (1.0 - fm);
		soil.sompool[SURFSTRUCT].nmass += litter_leaf_n * (1.0 - fm);
		soil.sompool[SURFMETA].cmass += pft.litter_leaf * fm;
		soil.sompool[SURFMETA].nmass += litter_leaf_n * fm;

		// Fire
		litterme[0] += pft.litter_leaf * (1.0 - fm) * pft.pft.litterme;
		litterme[1] += pft.litter_leaf * fm * pft.pft.litterme;

		fireresist[0] += pft.litter_leaf * (1.0 - fm) * pft.pft.fireresist;
		fireresist[1] += pft.litter_leaf * fm * pft.pft.fireresist;

		// NB: reproduction litter cannot contain nitrogen!!

		ligcmass_new = pft.litter_leaf * (1.0 - fm) * LIGCFRAC_LEAF;

		if (negligible(soil.sompool[SURFSTRUCT].cmass))
			soil.sompool[SURFSTRUCT].ligcfrac = 0.0;
		else
			soil.sompool[SURFSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old)/
				soil.sompool[SURFSTRUCT].cmass;

		// Remove association with vegetation
		pft.litter_leaf = 0.0;
		pft.nmass_litter_leaf = 0.0;
		pft.litter_repr = 0.0;

		// ROOT

		// Calculate inputs to soil structural and metabolic litter

		// Root litter lignin:N ratio
		if (!negligible(pft.nmass_litter_root))
			lton = max(0.0, LIGCFRAC_ROOT * pft.litter_root / pft.nmass_litter_root);
		else
			lton = max(0.0, LIGCFRAC_ROOT * pft.pft.cton_root_avr / (1.0 - nrelocfrac));

		// Metabolic litter fraction for root litter (Fm, Parton et al 1993, Eqn 1)
		fm = max(0.0, 0.85 - lton * 0.013);

		if (fm < 0.0 || fm > 1.0) 
			dprintf("Year %d ROOT fm %g pft %s\n", date.year, fm, (char*)pft.pft.name);

		ligcmass_new = pft.litter_root * (1.0 - fm) * LIGCFRAC_ROOT;
		ligcmass_old = soil.sompool[SOILSTRUCT].cmass * soil.sompool[SOILSTRUCT].ligcfrac;

		// Add to pools and update lignin fraction in structural pool
		soil.sompool[SOILSTRUCT].cmass += pft.litter_root * (1.0 - fm);
		soil.sompool[SOILSTRUCT].nmass += pft.nmass_litter_root * (1.0 - fm);
		if (negligible(soil.sompool[SOILSTRUCT].cmass))
			soil.sompool[SOILSTRUCT].ligcfrac = 0.0;
		else
			soil.sompool[SOILSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
				soil.sompool[SOILSTRUCT].cmass;
		soil.sompool[SOILMETA].cmass += pft.litter_root * fm;
		soil.sompool[SOILMETA].nmass += pft.nmass_litter_root * fm;
		
		// Remove association with vegetation
		pft.litter_root = 0.0;
		pft.nmass_litter_root = 0.0;

		// WOOD

		if (pft.pft.lifeform == TREE && !negligible(pft.litter_wood)) {
			// Woody debris enters a woody litter pool as described in
			// Kirschbaum and Paul (2002).

			// Coarse woody debris

			ligcmass_new = max(0.0, pft.litter_wood) * LIGCFRAC_WOOD;
			ligcmass_old = soil.sompool[SURFCWD].cmass * soil.sompool[SURFCWD].ligcfrac;

			if (pft.litter_wood < 0.0)
				dprintf("Year %d pft %s Negative litter wood %g \n", date.year, (char*)pft.pft.name, pft.litter_wood);

			// Add to structural pool and update lignin fraction in pool
			soil.sompool[SURFCWD].cmass += pft.litter_wood;
			soil.sompool[SURFCWD].nmass += pft.nmass_litter_wood;
			if (negligible(soil.sompool[SURFCWD].cmass))
				soil.sompool[SURFCWD].ligcfrac = 0.0;
			else {
				double ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFCWD].cmass;
				soil.sompool[SURFCWD].ligcfrac = ligcfrac;
			}

			// Fire
			litterme[2] += pft.litter_wood * pft.pft.litterme;
			fireresist[2] += pft.litter_wood * pft.pft.fireresist;
		
			// Update vegetation
			pft.litter_wood = 0.0;
			pft.nmass_litter_wood = 0.0;
		}

		patch.pft.nextobj();
	}

	// FIRE
	if (soil.sompool[SURFSTRUCT].cmass > 0.0) {
		soil.sompool[SURFSTRUCT].litterme = litterme[0] / soil.sompool[SURFSTRUCT].cmass;
		soil.sompool[SURFSTRUCT].fireresist = fireresist[0] / soil.sompool[SURFSTRUCT].cmass;
	}
	if (soil.sompool[SURFMETA].cmass > 0.0) {
		soil.sompool[SURFMETA].litterme = litterme[1] / soil.sompool[SURFMETA].cmass;
		soil.sompool[SURFMETA].fireresist = fireresist[1] / soil.sompool[SURFMETA].cmass;
	}
	if (soil.sompool[SURFCWD].cmass > 0.0) {
		soil.sompool[SURFCWD].litterme = litterme[2] / soil.sompool[SURFCWD].cmass;
		soil.sompool[SURFCWD].fireresist = fireresist[2] / soil.sompool[SURFCWD].cmass;
	}

	// Set N:C ratio of surface microbial pool based on C:N ratio of litter from all PFTs
	// Parton et al 1993 Fig 4. Dry mass litter == cmass litter * 2
	if (!negligible(litter_cmass))
		setntoc(soil, litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, 0.02);
}


/// LEACHING   
/** Leaching fractions for both organic and mineral leaching
 *  Should be called every day in both daily and monthly mode
 */
void leaching(Soil& soil) {
	
	if (!negligible(soil.dperc) && ifleachn) {

		// Leaching from available N mineral pool
		// using Parton et al. eqn. 13 
		soil.minleachfrac_daily[date.day] = min(1.0, soil.dperc * 0.1 / 18.0 * (0.2 + 0.7 * soil.soiltype.sand_frac) * 365.0 / 12.0);		

		// Leaching from decayed organic N
		// using Parton et al. eqn. 8
		soil.orgleachfrac_daily[date.day] = min(1.0, soil.dperc * 0.1 / 18.0 * (0.01 + 0.04 * soil.soiltype.sand_frac) * 365.0 / 12.0);
	}
	else { 
		soil.minleachfrac_daily[date.day] = 0.0;
		soil.orgleachfrac_daily[date.day] = 0.0;
	}
}

/// VEGETATION N UPTAKE  
/** Daily vegetation uptake of mineral N
 *  Partitioned among individuals according to this year's N demand
 *  Distributed through the year according to individual daily assimilation
 */
void vegetation_n_uptake(Patch& patch,Pftlist& pftlist) {

	// Daily N demand given by:
	//	 For individual:
	//     (1)  ndemand_day = dassim/aassim*ndemand_indiv
    //          where
	//          dassim = daily assimilation above base value 0
	//          aassim = annual sum of daily assimilation above base value 0
	//          ndemand_indiv = this year's N demand for growth by this individual 
	//							without any N limitation
	//
	//	 For patch:
	//          ndemand_patch_day = sum of ndemand_day over all individuals
	//                        
	// Actual N uptake for each day and individual given by:
	//     (3)  nuptake_day = ndemand_day * fnuptake
	//     where
	//     (4)  fnuptake = min(patch.nsupply/patch.ndemand,1.0)
	//	   or fnuptake is determined per individual 
	//
	// N fixation is done on patch basis using Cleveland 1999 approach
	//
	// To be called on last day of year following SOM dynamics

	const double EPS = 1e-12;

	double ndemand_day, nuptake_day;
	double excessn;
	double ndemand = 0.0;

	Vegetation& vegetation=patch.vegetation;
	Soil& soil = patch.soil;	

	// ANNUAL N SUPPLY
	// Potential N supply is remaining pool from last year
	// PLUS annual deposition
	// PLUS estimate of annual N fixation
	// PLUS sum of daily mineralisation MINUS sum of daily immobilisation

	// N deposition
	soil.andep = patch.stand.gridcell.climate.andep;

	// N fixation
	soil.anfix = max((nfix_a * patch.aaet + nfix_b) / 100000.0, 0.0);	

	// N mineralisation and immobilisation
	soil.anmin = 0.0;
	soil.animmob = 0.0;

	for (int d=0;d<365;d++) {
		soil.anmin += soil.nmin_daily[d];
		soil.animmob += soil.nimmob_daily[d];
	}

	// Total N supply in patch
	patch.nsupply = soil.nmass_avail + soil.andep + soil.anfix+
		soil.anmin - soil.animmob - soil.aminleach;

	// ANNUAL N DEMAND FOR PATCH

	// Loop through individuals

	patch.ndemand = 0.0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		indiv_ndemand(patch.stand.gridcell, patch.stand.landcover, patch.fluxes, indiv);

		//	store N in individual reserve
		if (date.year > freenyears && !negligible(indiv.ndemand)){

			double max_n_reserve_uptake;

			if (!negligible(indiv.ndemand))
				max_n_reserve_uptake = min(1.0, max(0.0, (indiv.max_n_reserve - indiv.nmass_reserve) / indiv.ndemand));
			else
				max_n_reserve_uptake = 0.0;

			// if N storage is larger than what can be stored then don't store more
			if (indiv.nmass_reserve > indiv.max_n_reserve)
				indiv.n_reserve_uptake = 0.0;
			// fill up storage to max
			else if (indiv.max_n_reserve - indiv.nmass_reserve < max_n_reserve_uptake * indiv.ndemand)
				indiv.n_reserve_uptake = (indiv.max_n_reserve - indiv.nmass_reserve) / indiv.ndemand;
			// store as much as possible 
			else
				indiv.n_reserve_uptake = max_n_reserve_uptake;

			indiv.ndemand *= (1.0 + indiv.n_reserve_uptake);
		}

		// Sum nitrogen demand of individuals with positive assimilation

		indiv.aassim = 0.0;
		for (int d=0; d<365; d++)
			indiv.aassim += max(0.0, indiv.dassim[d]);
		
		if (!negligible(indiv.aassim)) 
			patch.ndemand += indiv.ndemand;

		 vegetation.nextobj();
	}

	// Create individuals that determines amount of N that each pft has for establishment
	if (date.year >= freenyears)
		ndemand_new_est(patch, pftlist, patch.ndemand);

	// Rescale demand to not exceed supply (Eqn 4)
	if (patch.nsupply <= 0.0) 
		patch.fnuptake = 0.0; 
	else if (patch.ndemand > patch.nsupply) 
		patch.fnuptake = patch.nsupply / patch.ndemand;	
	else 
		patch.fnuptake = 1.0;

	// If soil available N is above the value for minimum SOM C:N ratio, then
	// N fixation is reduced (N rich soils)
	if (patch.fnuptake == 1.0 && patch.nsupply > patch.ndemand + nmin_balance_max) {
		if (soil.anfix <= patch.nsupply - (patch.ndemand + nmin_balance_max)) {
			patch.nsupply -= soil.anfix;
			soil.anfix = 0.0;
		}
		else {
			soil.anfix -= patch.nsupply - (patch.ndemand + nmin_balance_max);
			patch.nsupply = (patch.ndemand + nmin_balance_max);
		}
	}

	// Individual fnuptake
	if (patch.fnuptake < 1.0 && patch.fnuptake > 0.0) {
		if(ifindiv_fnuptake)
			indiv_fnuptake(vegetation, patch.nsupply, patch.fnuptake);
		else
			// Resolve Raingreen nitrogen demand
			raingreen_ndemand(vegetation, patch.ndemand, patch.nsupply, patch.fnuptake);
	}

	// VEGETATION N UPTAKE
	// Uptake in excess of daily supply permitted

	if (!negligible(patch.ndemand)) { // (some vegetation N uptake this year)
		
		// Loop through days of year

		for (int d=0;d<365;d++) {

			// Loop through individuals

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();

				if (d == 0)
					indiv.nuptake = 0.0;

				if (!ifindiv_fnuptake || patch.fnuptake == 1.0 || patch.fnuptake == 0.0)
					indiv.fnuptake = patch.fnuptake;
				
				if ((!negligible(indiv.aassim) && indiv.dassim[d] > 0.0)) {

					// Daily N demand by this individual (Eqn 1)
					
					ndemand_day = indiv.dassim[d] / indiv.aassim * indiv.ndemand;

					// Daily N uptake by this individual (Eqn 3)
					nuptake_day = ndemand_day * indiv.fnuptake;

					// Add to individual's nitrogen stores
					indiv.nstore += nuptake_day;
					
					// Add to yearly N uptake
					indiv.nuptake += nuptake_day;
				}
				// ... on to next individual
				vegetation.nextobj();
			}
		}
	}

	// Store N uptake by establishment individuals and then kill them
	if (date.year >= freenyears) {
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			if (indiv.age == -9999) {
				patch.pft[indiv.pft.id].nstore_est += indiv.nstore;
				vegetation.killobj();
			}
			else				
				vegetation.nextobj();	// ... on to next individual
		}
	}

	// EXCESS MINERAL N
	// Return remaining N to soil store for next year

	excessn = patch.nsupply - patch.ndemand * patch.fnuptake;

	// Should never be negative! (allow it for very small values for now ...)
	if (excessn < -EPS && ifnlim && date.year > freenyears)
		dprintf("Year %d vegetation_n_uptake: patch %d age %d Unexpected NEGATIVE value (%g) for annual excess mineral N before leach (%g)\n",
			date.year, patch.id, patch.age, excessn, patch.nsupply - patch.ndemand * patch.fnuptake);

	if (date.year > freenyears)
		soil.nmass_avail = excessn;
	else
		soil.nmass_avail = 0.0;
}


// Variables for checking N Balance! only works with one patch
double old_total = 0.0;
double old_vegn = 0.0;
double old_vegstore = 0.0;
double old_estn = 0.0;
double old_centuryn = 0.0;
double old_nmass_avail = 0.0;
double old_littern = 0.0;
double old_leachn = 0.0;
double nadded = 0.0;


/// Checking nitrogen balance  
/** Function to check if nitrogen is in balance
 */
void check_nbalance(Patch& patch, bool print) {

	double vegn, centuryn, littern, estn, vegstore;
	int p;

	Soil& soil=patch.soil;
	Vegetation& vegetation=patch.vegetation;

	if (patch.id == 0) {

		// Work out total ecosystem N for checking
		vegn = 0.0;
		vegstore = 0.0;
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();
			if (indiv.alive) {	
				vegn += indiv.nmass_leaf + indiv.nmass_root + indiv.nmass_sap + indiv.nmass_heart;
				vegstore += indiv.nstore + indiv.nmass_reserve; 
			}

			vegetation.nextobj();
		}

		centuryn = 0.0;
		for (p=0;p<NSOMPOOL-1;p++) {
			centuryn += soil.sompool[p].nmass;
		}	

		littern = 0.0;
		estn = 0.0;
		patch.pft.firstobj();
		while (patch.pft.isobj) {
			Patchpft& pft=patch.pft.getobj();

			littern += pft.nmass_litter_leaf +
				pft.nmass_litter_root +
				pft.nmass_litter_wood;

			estn += pft.nstore_est;

			patch.pft.nextobj();
		}

		nadded += soil.anfix + soil.andep;

		if (print && date.year > nyear_spinup) 
			dprintf("N BALANCE - difference over %d years: %g\n",
				date.year - nyear_spinup, old_total + nadded - (vegn + centuryn + soil.nmass_avail + vegstore + estn + littern + soil.sompool[LEACHED].nmass));

		if (date.year == nyear_spinup) {
			old_vegn = vegn;
			old_vegstore = vegstore;
			old_estn = estn;
			old_centuryn = centuryn;
			old_nmass_avail = soil.nmass_avail;
			old_littern = littern;
			old_leachn = soil.sompool[LEACHED].nmass;
			nadded = 0.0;

			old_total = vegn + centuryn + soil.nmass_avail + vegstore + estn + littern + soil.sompool[LEACHED].nmass + nadded;
		}
	}
}

/// SOM CENTURY DYNAMICS
/** To be called each simulation day for each modelled patch, following update
 *  of soil temperature and soil water.
 *  Transfers litter on first day, performes leaching and decomposition, and 
 *  on last day nitrogen uptake.
 */
void som_dynamics_century(Patch& patch,Pftlist& pftlist) {

	// First day of year only
	if (date.day == 0) { 

		// Function to check the N balance
		check_nbalance(patch, false);	

		// Transfer last year's litter to SOM pools
		transfer_litter(patch, patch.soil);
	}

	if (date.dayofmonth == 0) 
		patch.fluxes.mcflux_soil[date.month] = 0.0;

	// Potential daily leaching fraction for mineral N

	leaching(patch.soil);

	// Daily or monthly decomposition and fluxes between SOM pools
	
	somfluxes(patch, patch.soil, patch.fluxes);	

	if (date.islastmonth && date.islastday) {

		// Last day of year

		// Distribute plant N uptake and leaching of mineral nitrogen throughout the past year
		// Calculate mineral N pool at end of year
		
		vegetation_n_uptake(patch, pftlist);	
	}

	// Plotting soil pools
	Soil& soil=patch.soil;
	if (date.day == 0 && (date.year%10 == 0) && patch.id == 0 && ifcentury) {
		
		double total_cpool = 0.0;
		double total_npool = 0.0;
		for (int p=0;p<NSOMPOOL-1;p++) {
			total_cpool += soil.sompool[p].cmass;
			total_npool += soil.sompool[p].nmass;
		}
		
		plot("century C","surfstruct", date.year, soil.sompool[SURFSTRUCT].cmass);
		plot("century C","surfmeta", date.year, soil.sompool[SURFMETA].cmass);
		plot("century C","surfcwd", date.year, soil.sompool[SURFCWD].cmass);
		plot("century C","surfmicro", date.year, soil.sompool[SURFMICRO].cmass);
		plot("century C","soilstruct", date.year, soil.sompool[SOILSTRUCT].cmass);
		plot("century C","soilmeta", date.year, soil.sompool[SOILMETA].cmass);		
		plot("century C","soilmicro", date.year, soil.sompool[SOILMICRO].cmass);
		plot("century C","humussom", date.year, soil.sompool[SURFHUMUS].cmass);
		plot("century C","slowsom", date.year, soil.sompool[SLOWSOM].cmass);
		plot("century C","passivesom", date.year, soil.sompool[PASSIVESOM].cmass); 
		plot("century C","total", date.year, total_cpool); 
 
		plot("century N","surfstruct", date.year, soil.sompool[SURFSTRUCT].nmass);
		plot("century N","surfmeta", date.year, soil.sompool[SURFMETA].nmass);
		plot("century N","surfcwd", date.year, soil.sompool[SURFCWD].nmass);
		plot("century N","surfmicro", date.year, soil.sompool[SURFMICRO].nmass);
		plot("century N","soilstruct", date.year, soil.sompool[SOILSTRUCT].nmass);
		plot("century N","soilmeta", date.year, soil.sompool[SOILMETA].nmass);
		plot("century N","soilmicro", date.year, soil.sompool[SOILMICRO].nmass);
		plot("century N","humussom", date.year, soil.sompool[SURFHUMUS].nmass);
		plot("century N","slowsom", date.year, soil.sompool[SLOWSOM].nmass);
		plot("century N","passivesom", date.year, soil.sompool[PASSIVESOM].nmass);
		plot("century N","total", date.year, total_npool); 
	}

}

/// Choose between CENTURY or standard LPJ SOM dynamics
/**
*/
void som_dynamics(Patch& patch,Pftlist& pftlist) {

	if (ifcentury) som_dynamics_century(patch,pftlist);
	else som_dynamics_lpj(patch);
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Brunner, A. and J. P. Kimmins (2003). "Nitrogen fixation in coarse woody debris of Thuja 
//   plicata and Tsuga heterophylla forests on northern Vancouver Island." Canadian Journal 
//   of Forest Research-Revue Canadienne De Recherche Forestiere 33(9): 1670-1682.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Cleveland C C (1999) Global patterns of terrestrial biological nitrogen (N2) fixation
//   in natural ecosystems. GBC 13: 623-645
// Foley J A 1995 An equilibrium model of the terrestrial carbon budget
//   Tellus (1995), 47B, 310-319
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A process-based, 
//   terrestrial biosphere model of ecosystem dynamics (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Kirschbaum, M. U. F. and K. I. Paul (2002). "Modelling C and N dynamics in forest soils 
//   with a modified version of the CENTURY model." Soil Biology & Biochemistry 34(3): 341-354.
// Meentemeyer, V. (1978) Macroclimate and lignin control of litter decomposition
//   rates. Ecology 59: 465-472.
// Parton (2010) ForCent model development and testing using the Enriched Background 
//	 Isotope Study experiment JoGR 115: 
// Zaehle, S. & Friend, A. D. 2010. Carbon and nitrogen cycle dynamics in the O-CN land surface 
//   model: 1. Model description, site-scale evaluation, and sensitivity to parameter estimates. 
//   Global Biogeochemical Cycles, 24.
