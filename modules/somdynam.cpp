///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Soil organic matter dynamics
// Header file name:      somdynam.h
// Source code file name: somdynam.cpp
// Written by:            Ben Smith
// Version dated:         2002-09-22
// Updated:               2010-11-22


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

		// Calculate respiration temperature response if not yet done for this day

		if (soil.last_gtemp!=date.day) {
			respiration_temperature_response(soil.temp,soil.gtemp);
			soil.last_gtemp=date.day;
		}

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
		patch.fluxes.dcflux_soil=cflux;
		patch.fluxes.acflux_soil+=cflux;

		// Solve SOM pool sizes at end of year given by soil.solvesom_end

		if (date.year==soil.soiltype.solvesom_end && date.islastmonth && date.islastday)
			equilsom(soil);
	}
}

/////////////////////////////////////////////////
// CENTURY SOM DYNAMICS

void setntoc(Soil& soil,double fac,pooltype pool,double cton_max,double cton_min,
	double fmin,double fmax) {

	// Set N:C ratios for active, passive and SOM pools based on mineral N pool
	// or litter N fraction
	// (Parton et al 1993, Fig 4)

	// Crucial pool for N limitation is the active pool.
	if (fac<=fmin) 
		soil.sompool[pool].ntoc=1.0/cton_max;
	else if (fac>=fmax) 
		soil.sompool[pool].ntoc=1.0/cton_min;
	else {
		soil.sompool[pool].ntoc=1.0/(cton_min+(cton_max-cton_min)*
			(fmax-fac)/(fmax-fmin));
	}
}

// void decayrates(Soil& soil,double temp_soil,double wcont_soil) {
void decayrates(Soil& soil,double temp_soil,double wcont_soil,double net_nmass) {	// GUESSNFIX wood

	// Calculates CENTURY instantaneous decay rates given soil temperature, water
	// content of upper soil layer

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 1993, Eqn 2-4, K_I/365)
	//const double K_MAX[]={1.1e-2,1.3e-2,2.0e-2,1.6e-2,4.1e-2,5.1e-2,5.5e-4,1.2e-5}; // GUESSN DayCent

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	//const double K_MAX[]={6.8e-3,1.3e-2,3.0e-2,4.4e-4,1.9e-2,2.7e-2,5.1e-2,1.3e-3,6.8e-6};	// GUESSN ForCent

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris  decay	
	const double K_MAX[]={6.8e-3,1.3e-2,3.0e-2,4.4e-4,1.9e-2,2.7e-2,3.3e-3,5.1e-2,1.3e-3,6.8e-6};	// GUESSN ForCent + Kirschbaum

	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	double texture_mod=1.0-0.75*(soil.soiltype.clay_frac+soil.soiltype.silt_frac);

	double temp_mod;
	double moist_mod;
	double nmass_mod;	// GUESSN
	double wfps;
	double k;
	int p;

	// Calculate decomposition temperature modifier (in range 0-1)
	// [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
	// effect of soil temperature, Parton et al 1993, Fig 2)

	if (temp_soil>0.0)
		temp_mod=max(0.0,
			0.0326+0.00351*pow(temp_soil,1.652)-pow(temp_soil/41.748,7.19));
	else
		temp_mod=0.0;

	// Calculate decomposition moisture modifier (in range 0-1)
	// Friend et al 1997, Eqn 53
	// (Parton et al 1993, Fig 2)

	// wfps=wcont_soil*100.0/1.72;

	// Updated
	// water holding capacity at wilting point (wp) and ratio between saturation capacity and field capacity (f_FC) 
	// is calculated with the help of Cosby et al 1984;
	wfps=(wcont_soil*soil.soiltype.awc[0]+soil.soiltype.wp[0])*100.0/(soil.soiltype.f_FC[0]*soil.soiltype.awc[0]);
			

	if (wfps<60.0)
		moist_mod=exp((wfps-60.0)*(wfps-60.0)/-800.0);
	else
		moist_mod=0.000371*wfps*wfps-0.0748*wfps+4.13;

	
	// GUESSNFIX wood
	// N limitation modifier for decomposition
	//if (ifnlim)
	//	nmass_mod=max(1.0-exp(-pow(max((net_nmass+2.0e-5)/(2.0e-5),0.0),4.0)),0.5);
	//else
		nmass_mod=1.0;

	for (p=0;p<NSOMPOOL;p++) {

		// Calculate decay constant (annual basis)
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		k=K_MAX[p]*temp_mod*moist_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2

		if (p==SURFSTRUCT || p==SOILSTRUCT)
			k*=exp(-3.0*soil.sompool[p].ligcfrac);
		else if (p==SOILMICRO)
			k*=texture_mod;

		// GUESSNFIX wood
		if (p==SURFCWD || p==SOILSTRUCT || p==SURFSTRUCT || p==SOILMETA || p==SURFMETA)
			k*=nmass_mod;

		// Calculate fraction of C pool remaining after today's decomposition

		soil.sompool[p].frc=exp(-k);	
	}
}

void transferdecomp(Soil& soil,pooltype donor,pooltype receiver,
	double frac,double respfrac,double& respsum,double& nmin_actual,
	double& nimmob) {

	// Transfers specified fraction (frac) of today's decomposition in donor pool type
	// to receiver pool, transferring fraction respfrac of this to the accumulated CO2
	// flux respsum (representing total microbial respiration today)

	// decrement in donor C pool and N pools
	double cdec=soil.sompool[donor].cdec*frac;
	double ndec;
	
	if (!negligible(soil.sompool[donor].cmass))
		ndec=cdec*soil.sompool[donor].nmass/soil.sompool[donor].cmass;
	else ndec=0.0;

	// associated N increment in receiver pool (Friend et al 1997, Eqn 49)
	double ninc=cdec*(1.0-respfrac)*soil.sompool[receiver].ntoc;

	// if increase in receiver N greater than decrease in donor N,
	// balance must be immobilisation from mineral N pool
	// otherwise balance is N mineralisation
	if (ninc>ndec) 
		nimmob+=ninc-ndec;
	else 
		nmin_actual+=ndec-ninc;

	// "Transfer" C and N to receiver

	soil.sompool[receiver].delta_cmass+=cdec*(1.0-respfrac);
	soil.sompool[receiver].delta_nmass+=ninc;

	// Transfer microbial respiration
		
	respsum+=cdec*respfrac;

	// Decrease N mineralisation sum relative to donor pool
	soil.sompool[donor].ndec-=ninc;	
}

// GUESSN
void somfluxes(Patch& patch, Soil& soil,Fluxes& fluxes) {	

	// Daily or monthly fluxes between the eight CENTURY pools, and CO2 release to the atmosphere
	// Parton et al 1993, Fig 1; Comins & McMurtrie 1993, Appendix A

	int p,d;
	double csp,csa,respfrac,cap;
	double respsum=0.0;
	double nmin_actual=0.0; // actual (not net) N mineralisation
	double nimmob=0.0;		// N immobilisation
	double nmin_balance;
	double N_demand;

	if (date.day==0 && (!ifdailysetntoc || !ifdailydecomp)) {

		// First day of year

		// Set N:C ratios for active, passive and SOM pools based on mean mineral N pool for past year
		// (Parton et al 1993, Fig 4)

		if (date.year==0) 
			nmin_balance=0.0;
		else 
			nmin_balance=soil.nmin_annual+soil.ndep_annual-soil.nimmob_annual-soil.nleach_annual;

		/*// DayCent values 
		setntoc(soil,nmin_balance,SLOWSOM,20.0,12.0,0.0,nmass_avail_max);
		
		setntoc(soil,nmin_balance,PASSIVESOM,10.0,7.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SOILMICRO,15.0,3.0,0.0,nmass_avail_max);*/

		// ForCent values
		setntoc(soil,nmin_balance,SLOWSOM,30.0,15,0.0,nmass_avail_max);
		
		setntoc(soil,nmin_balance,PASSIVESOM,10.0,6.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SOILMICRO,15.0,6.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SURFHUMUS,30.0,15,0.0,nmass_avail_max);
	}
	else if(ifdailysetntoc)	
	{		
		// Instead of using the soil.nmass_avail at day==0, nmass_avail is
		// "updated" each day depending on mineralization, immobilization
		// N deposition, and plant uptake (no N limitation). 
		// Then ntoc ratios is calculated each day
							
		// First day of year

		if (date.day == 0)
			soil.setntoc_nmass_avail = soil.nmass_avail; 
		else
			// Update "daily" nmass available 
			soil.setntoc_nmass_avail += soil.daily_minimmndep;
		
		// Loop through individuals

		N_demand = 0.0;
		
		Vegetation& vegetation = patch.vegetation;
		vegetation.firstobj();
		while (vegetation.isobj) 
		{
			Individual& indiv=vegetation.getobj();

			// For this individual ...

			double NPPp;
			NPPp = indiv.assim-indiv.resp;

			if (NPPp > 0.0) {
				double ref_cton = (indiv.cmass_leaf+indiv.cmass_root+indiv.cmass_sap)/(indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap);
					//reference value of the C:N ratio for plant production

				N_demand += NPPp/ref_cton;
			}		
			vegetation.nextobj();
		}

		double N_availability = max(0.0,soil.setntoc_nmass_avail);
		double N_uptake = min(N_demand,N_availability);
		N_uptake = max(0.0,N_uptake);

		soil.setntoc_nmass_avail-=N_uptake;

		if (soil.setntoc_nmass_avail > 0.0) {// Leaching
			double leaching=soil.setntoc_nmass_avail*(soil.dperc/18.0*(0.2+0.7*soil.soiltype.sand_frac));
			soil.setntoc_nmass_avail-=leaching;
		}

		// Set N:C ratios for active, passive and SOM pools based on mean mineral N pool for past year
		// (Parton et al 1993, Fig 4)

		nmin_balance = soil.setntoc_nmass_avail;

		/*// DayCent values 
		setntoc(soil,nmin_balance,SLOWSOM,20.0,12.0,0.0,nmass_avail_max);
		
		setntoc(soil,nmin_balance,PASSIVESOM,10.0,7.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SOILMICRO,15.0,3.0,0.0,nmass_avail_max);*/

		// ForCent values
		setntoc(soil,nmin_balance,SLOWSOM,30.0,15,0.0,nmass_avail_max);
		
		setntoc(soil,nmin_balance,PASSIVESOM,10.0,6.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SOILMICRO,15.0,6.0,0.0,nmass_avail_max);

		setntoc(soil,nmin_balance,SURFHUMUS,30.0,15,0.0,nmass_avail_max);
	}

	if (ifdailydecomp) {

		// DAILY MODE

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no N limitation)

		//decayrates(soil,soil.temp,soil.wcont[0]);
		decayrates(soil,soil.temp,soil.wcont[0],soil.daily_minimmndep); // GUESSNFIX wood

	}
	else if (date.islastday) {

		// MONTHLY MODE (last day of month only)

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no N limitation)

		//decayrates(soil,soil.temp,soil.mwcontupper);
		decayrates(soil,soil.temp,soil.mwcontupper,soil.daily_minimmndep);	// GUESSNFIX wood

		// Convert fractional scalars from daily to monthly basis

		for (p=0;p<NSOMPOOL;p++) {
			soil.sompool[p].frc=pow(soil.sompool[p].frc,date.ndaymonth[date.month]);
		}
	}
	else return;  // Nothing to do if monthly mode but not last day of month

	// Calculate decomposition in all pools assuming these decay rates

	for (p=0;p<NSOMPOOL;p++) {
		soil.sompool[p].cdec=soil.sompool[p].cmass*(1.0-soil.sompool[p].frc);
		soil.sompool[p].ndec=soil.sompool[p].nmass*(1.0-soil.sompool[p].frc); 
		soil.sompool[p].delta_cmass=-soil.sompool[p].cdec;
		soil.sompool[p].delta_nmass=-soil.sompool[p].ndec;
	}

	// Partition potential decomposition among receiver pools

	// Donor pool SURFACE STRUCTURAL

	transferdecomp(soil,SURFSTRUCT,SURFMICRO,1.0-soil.sompool[SURFSTRUCT].ligcfrac,
		0.6,respsum,nmin_actual,nimmob);	// 

	transferdecomp(soil,SURFSTRUCT,SURFHUMUS,soil.sompool[SURFSTRUCT].ligcfrac,0.3,
		respsum,nmin_actual,nimmob);

	// Donor pool SURFACE METABOLIC

	transferdecomp(soil,SURFMETA,SURFMICRO,1.0,0.6,respsum,nmin_actual,nimmob);

	// Donor pool SOIL STRUCTURAL

	transferdecomp(soil,SOILSTRUCT,SOILMICRO,1.0-soil.sompool[SOILSTRUCT].ligcfrac,
		0.55,respsum,nmin_actual,nimmob);

	transferdecomp(soil,SOILSTRUCT,SLOWSOM,soil.sompool[SOILSTRUCT].ligcfrac,0.3,
		respsum,nmin_actual,nimmob);

	// Donor pool SOIL METABOLIC

	transferdecomp(soil,SOILMETA,SOILMICRO,1.0,0.55,respsum,nmin_actual,nimmob);

	// Donor pool SURFACE COARSE WOODY DEBRIS

	transferdecomp(soil,SURFCWD,SURFMICRO,1.0-soil.sompool[SURFCWD].ligcfrac,	// GUESSNFIX wood
		0.76,respsum,nmin_actual,nimmob);					// 0.55

	transferdecomp(soil,SURFCWD,SURFHUMUS,soil.sompool[SURFCWD].ligcfrac,0.76, // 0.3	// GUESSNFIX wood
		respsum,nmin_actual,nimmob);
	
	// Donor pool SURFACE MICROBE

	transferdecomp(soil,SURFMICRO,SURFHUMUS,1.0,0.6,respsum,nmin_actual,nimmob);

	// Donor pool SURFACE HUMUS

	transferdecomp(soil,SURFHUMUS,SLOWSOM,1.0,0.6,respsum,nmin_actual,nimmob);

	// Donor pool SLOW SOM
	
	// First work out partitioning coefficients (Fig 1, Parton et al 1993)

	csp=0.003-0.009*soil.soiltype.clay_frac;
	respfrac=0.55;
	csa=1.0-csp-respfrac;

	transferdecomp(soil,SLOWSOM,SOILMICRO,csa,0.0,respsum,nmin_actual,nimmob);

	transferdecomp(soil,SLOWSOM,PASSIVESOM,csp,0.0,respsum,nmin_actual,nimmob);

	// Account for respiration flux
	// N associated with this respiration is mineralised (Parton et al 1993, p 791)
	respsum+=respfrac*soil.sompool[SLOWSOM].cdec;

	if(!negligible(soil.sompool[SLOWSOM].cmass))
		nmin_actual+=respfrac*soil.sompool[SLOWSOM].cdec*soil.sompool[SLOWSOM].nmass/soil.sompool[SLOWSOM].cmass;	

	// Donor pool SOIL MICROBE

	// Fraction lost to microbial respiration (F_t, Parton et al 1993 Eqn 7)
	respfrac=max(0.0,0.85-0.68*(soil.soiltype.clay_frac+soil.soiltype.silt_frac));

	// Fraction entering passive SOM pool (Parton et al 1993, Eqn 9)
	cap=0.003+0.032*soil.soiltype.clay_frac;

	transferdecomp(soil,SOILMICRO,PASSIVESOM,cap,0.0,respsum,nmin_actual,nimmob);

	transferdecomp(soil,SOILMICRO,SLOWSOM,max(0.0,1.0-respfrac-cap),0.0,respsum,
		nmin_actual,nimmob);

	// Account for respiration flux
	// N associated with this respiration is mineralised (Parton et al 1993, p 791)
	respsum+=respfrac*soil.sompool[SOILMICRO].cdec;

	if(!negligible(soil.sompool[SLOWSOM].cmass))
		nmin_actual+=respfrac*soil.sompool[SOILMICRO].cdec*soil.sompool[SOILMICRO].nmass/soil.sompool[SOILMICRO].cmass;

	// Donor pool PASSIVE SOM

	transferdecomp(soil,PASSIVESOM,SOILMICRO,1.0,0.55,respsum,nmin_actual,nimmob);

	// Update pool sizes

	double nnmass = 0.0;

	for (p=0;p<NSOMPOOL;p++) {
		soil.sompool[p].cmass+=soil.sompool[p].delta_cmass;
		soil.sompool[p].nmass+=soil.sompool[p].delta_nmass;
		nnmass+=soil.sompool[p].delta_nmass;
	}

	// calculate the daily result of min, imm, and ndep
	soil.daily_minimmndep = nmin_actual-nimmob+soil.ndep_annual/365.0;

	// Transfer respiration sum to fluxes

	fluxes.dcflux_soil=respsum;
	fluxes.mcflux_soil[date.month]+=respsum;
	fluxes.acflux_soil+=respsum;

	// Store daily mineralisation and immobilisation to permit calculation of daily
	// mineral nitrogen balance at end of year

	if (ifdailydecomp) {
		soil.nmin_daily[date.day]=nmin_actual;
		soil.nimmob_daily[date.day]=nimmob;
	}
	else { // monthly mode - distribute current values evenly through the current month
		for (d=0;d<date.ndaymonth[date.month];d++) {
			soil.nmin_daily[date.day-d]=nmin_actual/(double)date.ndaymonth[date.month];
			soil.nimmob_daily[date.day-d]=nimmob/(double)date.ndaymonth[date.month];
		}
	}
}
// end GUESSN

// GUESSN
void transfer_litter(Patch& patch,Soil& soil) {

	// Call annually after growth, mortality and fire to transfer this year's litter
	// from vegetation to soil litter pools

	// Leaf, root and wood litter lignin fractions
	// Leaf and root fractions: Comins & McMurtrie 1993; Friend et al 1997
	// Wood fraction made up for now ... Thomas comment: check this!
	const double LIGCFRAC_LEAF=0.2;
	const double LIGCFRAC_ROOT=0.16;
	const double LIGCFRAC_WOOD=0.4;

	double lton;	// Leaf litter ligning to N ratio
	double fm;
	double ligcmass_old,ligcmass_new;
	double litter_leaf_n;

	double litter_nmass = 0.0;	
	double litter_cmass = 0.0;	

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft=patch.pft.getobj();

		// Calculate total litter C and N mass for set N:C ratio of surface microbial pool
		litter_nmass += (pft.nmass_litter_leaf + pft.nmass_litter_root + pft.nmass_litter_wood);
		litter_cmass += (pft.litter_leaf + pft.litter_root + pft.litter_wood);

		// LEAF

		// Calculate inputs to surface structural and metabolic litter

		litter_leaf_n=pft.nmass_litter_leaf; 

		// Leaf litter lignin:N ratio
		if (negligible(litter_leaf_n))
			lton=LIGCFRAC_LEAF*pft.pft.cton_leaf/(1.0-nrelocfrac);
		else
			lton=LIGCFRAC_LEAF*pft.litter_leaf/litter_leaf_n;

		// Metabolic litter fraction for leaf litter (Fm, Parton et al 1993, Eqn 1:
		// NB: incorrect/out-of-date intercept and slope given in Eqn 1; values used in
		// code of CENTURY 4.0 used instead)
		fm=max(0.0,0.85-lton*0.013);

		ligcmass_old=soil.sompool[SURFSTRUCT].cmass*soil.sompool[SURFSTRUCT].ligcfrac;

		if (fm < 0.0 || fm > 1.0)
			dprintf("Year %d LEAF fm %g\n",date.year,fm);
		if (litter_leaf_n < 0.0)
			dprintf("Year %d LEAF litter_leaf_n %g\n",date.year,litter_leaf_n);

		// Add to pools
		soil.sompool[SURFSTRUCT].cmass+=pft.litter_leaf*(1.0-fm);
		soil.sompool[SURFSTRUCT].nmass+=litter_leaf_n*(1.0-fm);
		soil.sompool[SURFMETA].cmass+=pft.litter_leaf*fm;
		soil.sompool[SURFMETA].nmass+=litter_leaf_n*fm;

		// NB: reproduction litter cannot contain nitrogen!!

		ligcmass_new=pft.litter_leaf*(1.0-fm)*LIGCFRAC_LEAF;

		if (negligible(soil.sompool[SURFSTRUCT].cmass))
			soil.sompool[SURFSTRUCT].ligcfrac=0.0;
		else
			soil.sompool[SURFSTRUCT].ligcfrac=(ligcmass_new+ligcmass_old)/
				soil.sompool[SURFSTRUCT].cmass;

		// Remove association with vegetation
		pft.litter_leaf=0.0;
		pft.nmass_litter_leaf=0.0;
		pft.litter_repr=0.0;

		// ROOT

		// Calculate inputs to soil structural and metabolic litter

		// Root litter lignin:N ratio
		if (negligible(pft.litter_root))
			lton=LIGCFRAC_ROOT*pft.pft.cton_root/(1-nrelocfrac);
		else
			lton=LIGCFRAC_ROOT*pft.litter_root/pft.nmass_litter_root;

		// Metabolic litter fraction for root litter (Fm, Parton et al 1993, Eqn 1)
		fm=max(0.0,0.85-lton*0.013);

		if (fm < 0.0 || fm > 1.0)
			dprintf("Year %d ROOT fm %g\n",date.year,fm);

		ligcmass_new=pft.litter_root*(1.0-fm)*LIGCFRAC_ROOT;
		ligcmass_old=soil.sompool[SOILSTRUCT].cmass*soil.sompool[SOILSTRUCT].ligcfrac;

		// Add to pools and update lignin fraction in structural pool
		soil.sompool[SOILSTRUCT].cmass+=pft.litter_root*(1.0-fm);
		soil.sompool[SOILSTRUCT].nmass+=pft.nmass_litter_root*(1.0-fm);
		if (negligible(soil.sompool[SOILSTRUCT].cmass))
			soil.sompool[SOILSTRUCT].ligcfrac=0.0;
		else
			soil.sompool[SOILSTRUCT].ligcfrac=(ligcmass_new+ligcmass_old)/
				soil.sompool[SOILSTRUCT].cmass;
		soil.sompool[SOILMETA].cmass+=pft.litter_root*fm;
		soil.sompool[SOILMETA].nmass+=pft.nmass_litter_root*fm;
		
		// Remove association with vegetation
		pft.litter_root=0.0;
		pft.nmass_litter_root=0.0;

		// WOOD

		if (pft.pft.lifeform==TREE && !negligible(pft.litter_wood)) {
			// Woody debris enters a woody litter pool as described in
			// Kirschbaum and Paul (2002).

			// Coarse woody debris

			double nmass_mod;
			if (ifnlim)
				nmass_mod = min(pow(patch.fuptake_patch,0.7),1.0);	// GUESSNFIX wood
			else
				nmass_mod = 1.0;

			ligcmass_new=max(0.0,pft.litter_wood*cwdtransfer*nmass_mod)*LIGCFRAC_WOOD;
			ligcmass_old=soil.sompool[SURFCWD].cmass*soil.sompool[SURFCWD].ligcfrac;

			if (pft.litter_wood < 0.0)
				dprintf("Year %d Negative litter wood %g \n",date.year,pft.litter_wood);

			// Add to structural pool and update lignin fraction in pool
			soil.sompool[SURFCWD].cmass+=pft.litter_wood*cwdtransfer*nmass_mod;
			soil.sompool[SURFCWD].nmass+=pft.nmass_litter_wood*cwdtransfer*nmass_mod;
			if (negligible(soil.sompool[SURFCWD].cmass))
				soil.sompool[SURFCWD].ligcfrac=0.0;
			else {
				double ligcfrac=(ligcmass_new+ligcmass_old)/
					soil.sompool[SURFCWD].cmass;
				soil.sompool[SURFCWD].ligcfrac=ligcfrac;
			}
		
			// Update vegetation
			pft.litter_wood*=(1.0-cwdtransfer*nmass_mod);
			pft.nmass_litter_wood*=(1.0-cwdtransfer*nmass_mod);
		}		

		patch.pft.nextobj();
	}

	// Set N:C ratio of surface microbial pool based on C:N ratio of litter from all PFTs
	// Parton et al 1993 Fig 4
	if (!negligible(litter_cmass))
		setntoc(soil,litter_nmass/(litter_cmass*2.0),SURFMICRO,20.0,10.0,0,0.02);
}
// end GUESSN

void leaching(Soil& soil) {

	// LEACHING
	// Should be called every day in both daily and monthly mode

	double leachfrac;

	// Leaching of organics from active pool (Parton et al 1993, Eqn 8)

	if (ifleachn)
		leachfrac=soil.dperc*0.1/18.0*(0.01+0.04*soil.soiltype.sand_frac);
	else
		leachfrac=0.0;

	soil.sompool[LEACHED].cmass+=soil.sompool[SOILMICRO].cmass*leachfrac;
	soil.sompool[LEACHED].nmass+=soil.sompool[SOILMICRO].nmass*leachfrac;
	
	soil.sompool[SOILMICRO].cmass*=(1.0-leachfrac);
	soil.sompool[SOILMICRO].nmass*=(1.0-leachfrac);

	// Leaching from mineral pool
	// Assume this affects daily mineral N excess after vegetation uptake
	// in proportion to baseflow as a fraction of total soil water

	//if (!negligible(soil.wcontmm_yesterday) && ifleachn) 
	//	soil.leachfrac_daily[date.day]=soil.dbaseflow/soil.wcontmm_yesterday;
	if (!negligible(soil.dperc) && ifleachn) 
		// using Parton et al. eqn. 13 instead
		soil.leachfrac_daily[date.day]=soil.dperc/18.0*(0.2+0.7*soil.soiltype.sand_frac);		
	else 
		soil.leachfrac_daily[date.day]=0.0;
}

void turnover_oecd_ndemand(double turnover_leaf,double turnover_root,double turnover_sap,lifeformtype lifeform,
		double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,
		double& nmass_leaf,double& nmass_root,double& nmass_sap,double& nstore_turnover,bool alive) {

	// DESCRIPTION
	// Transfers carbon from leaves and roots to litter, and from sapwood to heartwood
	// Version for OECD experiment:
	// For crops (specially labelled grass type) 50% of above-ground biomass transferred
	// to litter, remainder stored as a flux to the atmosphere (i.e. increments Rh)
	// (equal amount for each month)

	// guess2008 - new (indiv.)alive boolean throughout. Also, only turnover from 'alive' 
	// individuals is transferred to litter

	double turnover;

	// TREES AND GRASSES:

	// Leaf turnover
	turnover=turnover_leaf*cmass_leaf;
	cmass_leaf-=turnover;

	turnover=turnover_leaf*nmass_leaf;
	nmass_leaf-=turnover;
	if (alive) 
		nstore_turnover+=turnover*nrelocfrac;	

	// Root turnover
	turnover=turnover_root*cmass_root;
	cmass_root-=turnover;

	turnover=turnover_root*nmass_root;
	nmass_root-=turnover;
	if (alive) 
		nstore_turnover+=turnover*nrelocfrac;	

	if (lifeform==TREE) {
			
		// TREES ONLY:

		// Sapwood turnover by conversion to heartwood
		// NB: assumes N is translocated from sapwood prior to conversion to
		//     heartwood and that this is the same fraction that is conserved
		//     in conjunction with leaf and root shedding
		
		turnover=turnover_sap*cmass_sap;
		cmass_sap-=turnover;
		cmass_heart+=turnover;

		turnover=turnover_sap*nmass_sap;
		nmass_sap-=turnover;
		nstore_turnover+=turnover*nrelocfrac;	
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// ALLOCATION
// Function allocation is an internal function (do not call directly from framework);

// File scope global variables: used by function f below (see function allocation)

static double k1,k2,k3,b;
static double ltor_g;
static double cmass_heart_g;
static double cmass_leaf_g;

inline double f2(double& cmass_leaf_inc) {
	// Thomas this year's ndemand: changed name of this function also in allocation_ndemand() to f2 (from f)

	// Returns value of f2(cmass_leaf_inc), given by:
	//
	// f2(cmass_leaf_inc) = 0 =
	//   k1 * (b - cmass_leaf_inc - cmass_leaf_inc/ltor + cmass_heart) -
	//   [ (b - cmass_leaf_inc - cmass_leaf_inc/ltor)
	//   / (cmass_leaf + cmass_leaf_inc )*k3 ] ** k2
	//
	// See function allocation (below), Eqn (13)

	return k1*(b-cmass_leaf_inc-cmass_leaf_inc/ltor_g+cmass_heart_g)
		-pow((b-cmass_leaf_inc-cmass_leaf_inc/ltor_g)/(cmass_leaf_g+cmass_leaf_inc)*k3,
		k2);
}

void allocation_ndemand(double bminc,double cmass_leaf,double cmass_root,double cmass_sap,
	double cmass_debt,double cmass_heart,double ltor,double height,double sla,
	double wooddens,lifeformtype lifeform,double k_latosa,double k_allom2,
	double k_allom3,double& cmass_leaf_inc,double& cmass_root_inc,
	double& cmass_sap_inc,
	double& cmass_debt_inc,
	double& cmass_heart_inc) {

	// DESCRIPTION
	// Calculates changes in C compartment sizes (leaves, roots, sapwood, heartwood)
	// and litter for a plant individual as a result of allocation of biomass increment.
	// Assumed allometric relationships are given in function allometry below.

	// INPUT PARAMETERS
	// bminc       = biomass increment this time period on individual basis (kgC)
	// cmass_leaf  = leaf C biomass for last time period on individual basis (kgC)
	// cmass_root  = root C biomass for last time period on individual basis (kgC)
	// cmass_sap   = sapwood C biomass for last time period on individual basis (kgC)
	// cmass_heart = heartwood C biomass for last time period on individual basis (kgC)
	// ltor        = leaf to root mass ratio following allocation
	// height      = individual height (m)
	// sla         = specific leaf area (PFT-specific constant) (m2/kgC)
	// wooddens    = wood density (PFT-specific constant) (kgC/m3)
	// lifeform    = life form class (TREE or GRASS)
	// k_latosa    = ratio of leaf area to sapwood cross-sectional area (PFT-specific
	//               constant)
	// k_allom2    = constant in allometry equations
	// k_allom3    = constant in allometry equations

	// OUTPUT PARAMETERS
	// cmass_leaf_inc  = increment (may be negative) in leaf C biomass following
	//                   allocation (kgC)
	// cmass_root_inc  = increment (may be negative) in root C biomass following
	//                   allocation (kgC)
	// cmass_sap_inc   = increment (may be negative) in sapwood C biomass following
	//                   allocation (kgC)
	// cmass_heart_inc = increment in heartwood C biomass following allocation (kgC)
	// litter_leaf_inc = increment in leaf litter following allocation, on individual 
	//                   basis (kgC)
	// litter_root_inc = increment in root litter following allocation, on individual
	//                   basis (kgC)

	// MATHEMATICAL DERIVATION FOR TREE ALLOCATION
	// Allocation attempts to distribute biomass increment (bminc) among the living
	// tissue compartments, i.e.
	//   (1) bminc = cmass_leaf_inc + cmass_root_inc + cmass_sap_inc
	// while satisfying the allometric relationships (Shinozaki et al. 1964a,b; Waring
	// et al 1982, Huang et al 1992; see also function allometry, below) [** =
	// raised to the power of]:
	//   (2) (leaf area) = k_latosa * (sapwood xs area)
	//   (3) cmass_leaf = ltor * cmass_root
	//   (4) height = k_allom2 * (stem diameter) ** k_allom3
	// From (1) and (3),
	//   (5) cmass_sap_inc = bminc - cmass_leaf_inc -
	//         (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root
	// Let diam_new and height_new be stem diameter and height following allocation.
	// Then (see allometry),
	//   (6) diam_new = 2 * [ ( cmass_sap + cmass_sap_inc + cmass_heart )
	//         / wooddens / height_new / PI ]**(1/2)
	// From (4), (6) and (5),
	//   (7) height_new**(1+2/k_allom3) = 
	//         k_allom2**(2/k_allom3) * 4 * [cmass_sap + bminc - cmass_leaf_inc
	//         - (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root + cmass_heart]
	//         / wooddens / PI
	// Now,
	//   (8) wooddens = cmass_sap / height / (sapwood xs area)
	// From (8) and (2),
	//   (9) wooddens = cmass_sap / height / sla / cmass_leaf * k_latosa
	// From (9) and (1),
	//  (10) wooddens = (cmass_sap + bminc - cmass_leaf_inc -
	//         (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root)
	//          / height_new / sla / (cmass_leaf + cmass_leaf_inc) * k_latosa
	// From (10),
	//  (11) height_new**(1+2/k_allom3) =
	//         [ (cmass_sap + bminc - cmass_leaf_inc - (cmass_leaf + cmass_leaf_inc)
	//           / ltor + cmass_root) / wooddens / sla
	//           / (cmass_leaf + cmass_leaf_inc ) * k_latosa ] ** (1+2/k_allom3)
	//
	// Combining (7) and (11) gives a function of the unknown cmass_leaf_inc:
	//
	//  (12) f(cmass_leaf_inc) = 0 =
	//         k_allom2**(2/k_allom3) * 4/PI * [cmass_sap + bminc - cmass_leaf_inc
	//         - (cmass_leaf + cmass_leaf_inc) / ltor + cmass_root + cmass_heart]
	//         / wooddens -
	//         [ (cmass_sap + bminc - cmass_leaf_inc - (cmass_leaf + cmass_leaf_inc)
	//           / ltor + cmass_root) / (cmass_leaf + cmass_leaf_inc)
	//           / wooddens / sla * k_latosa] ** (1+2/k_allom3)
	//
	// Let k1 = k_allom2**(2/k_allom3) * 4/PI / wooddens
	//     k2 = 1+2/k_allom3
	//     k3 = k_latosa / wooddens / sla
	//     b  = cmass_sap + bminc - cmass_leaf/ltor + cmass_root
	//
	// Then,
	//  (13) f(cmass_leaf_inc) = 0 =
	//         k1 * (b - cmass_leaf_inc - cmass_leaf_inc/ltor + cmass_heart) -
	//         [ (b - cmass_leaf_inc - cmass_leaf_inc/ltor)
	//         / (cmass_leaf + cmass_leaf_inc )*k3 ] ** k2
	//
	// Numerical methods are used to solve Eqn (13) for cmass_leaf_inc

	const int NSEG=20; // number of segments (parameter in numerical methods)
	const int JMAX=40; // maximum number of iterations (in numerical methods)
	const double XACC=0.0001; // threshold x-axis precision of allocation solution
	const double YACC=1.0e-10; // threshold y-axis precision of allocation solution
	const double PI=3.14159265;
	const double CDEBT_MAXLOAN_DEFICIT=0.8; // maximum loan as a fraction of deficit
	const double CDEBT_MAXLOAN_MASS=0.2; // maximum loan as a fraction of (sapwood-cdebt)

	double cmass_leaf_inc_min;
	double cmass_root_inc_min;
	double x1,x2,dx,xmid,fx1,fmid,rtbis,sign;
	int j;
	double cmass_deficit,cmass_loan;

	if (ltor<1.0e-10) {
		
		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_leaf_inc=0.0;
		cmass_root_inc=bminc;

		if (lifeform==TREE) {
			cmass_sap_inc=-cmass_sap;
			cmass_heart_inc=-cmass_sap_inc;
		}

		return;
	}

	if (lifeform==TREE) {

		// TREE ALLOCATION

		cmass_heart_inc=0.0;
		cmass_sap_inc=0.0;	

		// Calculate minimum leaf increment to maintain current sapwood biomass
		// Given Eqn (2)

		if (height>0.0)
			cmass_leaf_inc_min=k_latosa*cmass_sap/(wooddens*height*sla)-cmass_leaf;
		else
			cmass_leaf_inc_min=0.0;

		// Calculate minimum root increment to support minimum resulting leaf biomass
		// Eqn (3)

		if (height>0.0)
			cmass_root_inc_min=k_latosa*cmass_sap/(wooddens*height*sla*ltor)-
				cmass_root;
		else
			cmass_root_inc_min=0.0;

		if (cmass_root_inc_min<0.0) { // some roots would have to be killed

			cmass_leaf_inc_min=cmass_root*ltor-cmass_leaf;
			cmass_root_inc_min=0.0;
		}

		// BLARP! C debt stuff
		if (ifcdebt) {
			cmass_deficit=cmass_leaf_inc_min+cmass_root_inc_min-bminc;
			if (cmass_deficit>0.0) {
				cmass_loan=max(min(cmass_deficit*CDEBT_MAXLOAN_DEFICIT,
					(cmass_sap-cmass_debt)*CDEBT_MAXLOAN_MASS),0.0);
				bminc+=cmass_loan;
				cmass_debt_inc=cmass_loan;
			}
			else cmass_debt_inc=0.0;
		}
		else cmass_debt_inc=0.0;

		if (cmass_root_inc_min>=0.0 && cmass_leaf_inc_min>=0.0 &&
			cmass_root_inc_min+cmass_leaf_inc_min<=bminc || bminc<=0.0) {

			// Normal allocation (positive increment to all living C compartments)
			// NOTE: includes allocation of zero or negative NPP, c.f. LPJF

			// Calculation of leaf mass increment (lminc_ind) satisfying Eqn (13)
			// using bisection method (Press et al 1986)

			// Set values for global variables for reuse by function f

			k1=pow(k_allom2,2.0/k_allom3)*4.0/PI/wooddens;
			k2=1.0+2/k_allom3;
			k3=k_latosa/wooddens/sla;
			b=cmass_sap+bminc-cmass_leaf/ltor+cmass_root;
			ltor_g=ltor;
			cmass_leaf_g=cmass_leaf;
			cmass_heart_g=cmass_heart;

			x1=0.0;
			x2=(bminc-(cmass_leaf/ltor-cmass_root))/(1.0+1.0/ltor);
			dx=(x2-x1)/(double)NSEG;

			if (cmass_leaf<1.0e-10) x1+=dx; // to avoid division by zero

			// Evaluate f(x1), i.e. Eqn (13) at cmass_leaf_inc = x1

			fx1=f2(x1);

			// Find approximate location of leftmost root on the interval
			// (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
			// change in sign of f(xmid) relative to f(x1).

			fmid=f2(x1);

			xmid=x1;

			while (fmid*fx1>0.0 && xmid<x2) {

				xmid+=dx;
				fmid=f2(xmid);
			}

			x1=xmid-dx;
			x2=xmid;

			// Apply bisection to find root on new interval (x1,x2)

			if (f2(x1)>=0.0) sign=-1.0;
			else sign=1.0;

			rtbis=x1;
			dx=x2-x1;

			// Bisection loop
			// Search iterates on value of xmid until xmid lies within
			// xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

			fmid=1.0; // dummy value to guarantee entry into loop
			j=0; // number of iterations so far

			while (dx>=XACC && fabs(fmid)>YACC && j<=JMAX) {

				dx*=0.5;
				xmid=rtbis+dx;

				fmid=f2(xmid);

				if (fmid*sign<=0.0) rtbis=xmid;
				j++;
			}

			// Now rtbis contains numerical solution for cmass_leaf_inc given Eqn (13)

			cmass_leaf_inc=rtbis;

			// Calculate increments in other compartments

			cmass_root_inc=(cmass_leaf_inc+cmass_leaf)/ltor-cmass_root; // Eqn (3)
			cmass_sap_inc=bminc-cmass_leaf_inc-cmass_root_inc; // Eqn (1)

			// guess2008 - extra check - abnormal allocation can still happen if ltor is very small
			if ((cmass_root_inc > 50 || cmass_root_inc < -50) && ltor < 0.0001) {
				cmass_leaf_inc=0.0;
				cmass_root_inc=bminc;

				if (lifeform==TREE) {
					cmass_sap_inc=-cmass_sap;
					cmass_heart_inc=-cmass_sap_inc;
				}

				return;			
			}

			// Convert killed sapwood to heartwood
			if(cmass_sap_inc<0.0)	
				if (-cmass_sap_inc > cmass_sap)
					cmass_heart_inc += cmass_sap;
				else
					cmass_heart_inc =- cmass_sap_inc;
		}
		else {

			// Abnormal allocation: reduction in some biomass compartment(s) to
			// satisfy allometry

			// Attempt to distribute this year's production among leaves and roots only
			// Eqn (3)

			cmass_leaf_inc=(bminc-cmass_leaf/ltor+cmass_root)/(1.0+1.0/ltor);

			if (cmass_leaf_inc>0.0) {

				// Positive allocation to leaves

				cmass_root_inc=bminc-cmass_leaf_inc; // Eqn (1)

				// Add killed roots (if any) to litter

				// guess2008 - back to LPJF way in this case
				// if (cmass_root_inc<0.0) litter_root_inc=-cmass_root_inc;
				if (cmass_root_inc<0.0) {
					cmass_leaf_inc = bminc;
					cmass_root_inc=(cmass_leaf_inc+cmass_leaf)/ltor-cmass_root; // Eqn (3)
				}

			}
			else {

				// Negative or zero allocation to leaves
				// Eqns (1), (3)

				cmass_root_inc=bminc;
				cmass_leaf_inc=(cmass_root+cmass_root_inc)*ltor-cmass_leaf;
			}

			// Calculate increase in sapwood mass (which must be negative)
			// Eqn (2)

			cmass_sap_inc=(cmass_leaf_inc+cmass_leaf)*wooddens*height*sla/k_latosa-
				cmass_sap;

			// Convert killed sapwood to heartwood
			if(cmass_sap_inc<0.0)	
				if (-cmass_sap_inc > cmass_sap)
					cmass_heart_inc += cmass_sap;
				else
					cmass_heart_inc =- cmass_sap_inc;
		}
	}
	else if (lifeform==GRASS || lifeform==CROP) {

		// GRASS ALLOCATION
		// Allocation attempts to distribute biomass increment (bminc) among leaf
		// and root compartments, i.e.
		//   (14) bminc = cmass_leaf_inc + cmass_root_inc
		// while satisfying Eqn(3)

		cmass_leaf_inc=(bminc-cmass_leaf/ltor+cmass_root)/(1.0+1.0/ltor);
		cmass_root_inc=bminc-cmass_leaf_inc;

		if (cmass_leaf_inc<0.0) {

			// Negative allocation to leaves

			cmass_root_inc=bminc;
			cmass_leaf_inc=(cmass_root+cmass_root_inc)*ltor-cmass_leaf; // Eqn (3)			
		}
		else if (cmass_root_inc<0.0) {

			// Negative allocation to roots

			cmass_leaf_inc=bminc;
			cmass_root_inc=(cmass_leaf+bminc)/ltor-cmass_root;
		}
	}
}

// GUESSN
double this_years_ndemand(double cmass_leaf,double cmass_root,double cmass_sap,double cmass_heart,double cmass_debt,
		double nmass_leaf,double nmass_root,double nmass_sap,double nmass_heart,double leafn_mean,
		double cton_leaf,double cton_root,double cton_sap,double anpp,double reprfrac,double wscal_mean,double ltor_max,
		double height,double sla,double wooddens,double k_latosa,double k_allom2,double k_allom3,
		phenologytype phenology,double aphen_raingreen, double leaflong,
		double turnover_leaf,double turnover_root,double turnover_sap,
		double nstore,double densindiv,bool alive,lifeformtype lifeform) {

	// DESCRIPTION
	// Calculates this year's N demand by doing a fake growth with tissue turnover
	// and allocation of fixed carbon to reproduction and new biomass, and by doing 
	// so determining the N demand. 
	// Accumulated NPP (assimilation minus maintenance and growth respiration) on
	// patch or modelled area basis assumed to be given by 'anpp' member variable for
	// each individual.

	double CDEBT_PAYBACK_RATE=0.2;
	double cmass_payback;
	double cmass_excess;

	double ndemand_uptake=0.0;	// This individuals year N demand
 
	double nstore_turnover=0.0;	// N retranslocated in turnover
	double ndemand_new_cton;		// N demand/gained from change in cton ratio
	double raingreen_ndemand=0.0;

	// Thomas this year's ndemand: declaration
	double bminc;
	double dval;

	// present N mass
	double before=nmass_leaf+nmass_root+nmass_sap+nmass_heart;
	double after;

	// C:N ratios
	double cton_leaf_new,cton_root_new,cton_sap_new,cton_sap_old;

	// Increases
	double cmass_leaf_inc,cmass_root_inc,cmass_sap_inc,cmass_heart_inc;
	double cmass_debt_inc;
	double litter_leaf_inc,litter_root_inc;

	double ltor;
	
	// C:N ratio for new and current biomass		
	if (ifvarycn) {
		if (!negligible(leafn_mean) && !negligible(cmass_leaf))
			cton_leaf_new=cmass_leaf/leafn_mean; 
				// actual mean leafN based on Vmax from photosynthesis
		else
			cton_leaf_new=cton_leaf;

		cton_root_new=
			cton_leaf_new*(cton_root/cton_leaf);
		cton_sap_new=
			cton_leaf_new*(cton_sap/cton_leaf);
	}
	else {
		cton_leaf_new=cton_leaf;
		cton_root_new=cton_root;
		cton_sap_new=cton_sap;
	}

	if (!negligible(nmass_sap))
			cton_sap_old=cmass_sap/nmass_sap;
		else
			cton_sap_old=cton_sap;

	// reproduction
	if (anpp>=0.0) 
		bminc=anpp*(1.0-reprfrac);
	else
		bminc=anpp;

	if (bminc >= 0 && phenology==RAINGREEN) {

		// Raingreen PFTs: reduce biomass increment to account for NPP
		// allocated to extra leaves during the past year.
		// Excess allocation to leaves given by:
		//   aphen_raingreen / ( leaf_longevity * 365) * cmass_leaf -
		//   cmass_leaf

		// BLARP! excess allocation to roots now also included (assumes leaf longevity = root longevity)

		cmass_excess=max((double)aphen_raingreen/
			(leaflong*365.0)*(cmass_leaf+cmass_root)-
			cmass_leaf-cmass_root,0.0);

		if (cmass_excess>bminc) cmass_excess=bminc;

		// Deduct from this year's C biomass increment
		if (alive) {
			bminc-=cmass_excess;
			raingreen_ndemand = cmass_excess/cton_leaf_new;
		}
	}

	// Set leaf:root mass ratio based on water stress parameter
	ltor=wscal_mean*ltor_max;

	turnover_oecd_ndemand(turnover_leaf,turnover_root,turnover_sap,lifeform,
						cmass_leaf,cmass_root,cmass_sap,cmass_heart,
						nmass_leaf,nmass_root,nmass_sap,nstore_turnover,alive);
	
	// ndemand_new_cton is the amount of N demanded after update of the cton ratios
	ndemand_new_cton = (cmass_leaf/cton_leaf_new-nmass_leaf)+	
		(cmass_root/cton_root_new-nmass_root);
	
	if (lifeform == TREE)
		ndemand_new_cton += (cmass_sap/cton_sap_new-nmass_sap); 	

	// N demand not associated with growth
	ndemand_uptake -= nstore_turnover - ndemand_new_cton - raingreen_ndemand;

	if (lifeform==TREE) { 

		// pay back part of cdebt
		if (ifcdebt && bminc>0.0) {
			cmass_payback=min(cmass_debt*CDEBT_PAYBACK_RATE,bminc);
			bminc-=cmass_payback;
			cmass_debt-=cmass_payback;
		}

		allocation_ndemand(bminc/densindiv,cmass_leaf/densindiv,
			cmass_root/densindiv,cmass_sap/densindiv,
			cmass_debt/densindiv,cmass_heart/densindiv,
			ltor,height,sla,wooddens,TREE,
			k_latosa,k_allom2,k_allom3,
			cmass_leaf_inc,cmass_root_inc,cmass_sap_inc,cmass_debt_inc,
			cmass_heart_inc);

		// Calculate N needed for this new biomass
		double nmass_inc =	
			(max(0.0,cmass_leaf_inc)*densindiv/cton_leaf_new +
			max(0.0,cmass_root_inc)*densindiv/cton_root_new +
			max(0.0,cmass_sap_inc)*densindiv/cton_sap_new);

		ndemand_uptake += nmass_inc;

		if (ndemand_uptake<0.0)
			ndemand_uptake=0.0;

		after = (cmass_leaf+cmass_leaf_inc*densindiv)/cton_leaf+(cmass_root+cmass_root_inc*densindiv)/cton_root+
						(cmass_sap+cmass_sap_inc*densindiv)/cton_sap+nmass_heart+cmass_heart_inc*densindiv/cton_sap*nrelocfrac;

		return ndemand_uptake;
	}
	else {
		allocation_ndemand(bminc,cmass_leaf,cmass_root,
			0.0,0.0,0.0,ltor,0.0,0.0,0.0,GRASS,0.0,
			0.0,0.0,cmass_leaf_inc,cmass_root_inc,dval,dval,dval);

		double nmass_inc =
			(max(0.0,cmass_leaf_inc)*densindiv/cton_leaf_new+
			max(0.0,cmass_root_inc)*densindiv/cton_root_new);

		ndemand_uptake += nmass_inc;

		if (ndemand_uptake<0.0)
			ndemand_uptake=0.0;

		return ndemand_uptake;
	}
}

// end GUESSN

//////////////////////////////////////////////////////////////////////////////////
// GUESSN Calculates individual fuptake based on nmass_root and crownarea
//

double nitrogen_uptake_strength(const Individual& indiv) {
      const double min_crownarea = 0.1;
      const double crownarea_power = 0.2;

      return indiv.nmass_root;///(max(min_crownarea,pow(indiv.crownarea,crownarea_power)));//*indiv.densindiv;
}

void indiv_fuptake(Vegetation& vegetation, double nsupply_patch, double fuptake) {

	double GRASS_part = 0.05;			// Grass should at least get 5% of total available N
	double GRASS_ndemand = 0.0;			// Grass total N demand
	bool GRASS_100 = false;				// if grass gets what it demands from its part of the total N supply 
	bool not_more_grass = false;		// keeping track of if GRASS can compite with TREEs for more N than what is
										// espacially assigned for GRASS (GRASS_part)
	double grass_uptake_decider = 0.0;
	double total_uptake_decider = 0.0;
	double temp_nsupply_patch = nsupply_patch;
	double ratio_uptake;				// how much N taken up per uptake strength
	bool full_uptake = true;			// if indiv.fuptake should be updated

	// GRASS
	// determine strength and demand of grasses
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		if (indiv.pft.lifeform == GRASS && !negligible(indiv.ndemand_uptake)) {
			GRASS_ndemand += indiv.ndemand_uptake;
			grass_uptake_decider += nitrogen_uptake_strength(indiv);
		}
		vegetation.nextobj();
	}

	// GRASS
	// Does grass get enough N from its part of the total
	if (GRASS_ndemand < GRASS_part*nsupply_patch) 
		GRASS_100 = true;
	else 
		GRASS_100 = false;


	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		indiv.fuptake = fuptake;

		// GRASS
		if (indiv.pft.lifeform == GRASS && GRASS_100 && !negligible(indiv.ndemand_uptake)) {

			// when grass part of total N is enough, then subtract it from total
			temp_nsupply_patch-=indiv.ndemand_uptake;
			// set uptake to meet demand
			indiv.fuptake = 1.0;
			// and subtract uptake strength as it will be added further down
			total_uptake_decider -= nitrogen_uptake_strength(indiv);
		}
		
		// TREE
		// Sum up uptake strengths
		if (!negligible(indiv.ndemand_uptake)) {
			double checkk=nitrogen_uptake_strength(indiv);
			total_uptake_decider += checkk;
		}

		vegetation.nextobj();
	}
	
	// Loop through indiv and decide their fuptake
	while (full_uptake){

		full_uptake = false;	
		
		// restore N supply and uptake decider if not_more_grass == true
		// (which can happen after the first round if there is a full_uptake)
		// so that it can be calculated if they might be able to take up more 
		// than just the GRASS part
		if (not_more_grass) {
			temp_nsupply_patch += GRASS_part*nsupply_patch;
			total_uptake_decider += grass_uptake_decider;
			not_more_grass=false;
		}

		// decide how much N that will be taken up by each uptake strength 
		if (total_uptake_decider > 0.0 && temp_nsupply_patch > 0.0)
			ratio_uptake = temp_nsupply_patch / total_uptake_decider;
		else
			ratio_uptake = 0.0;

		// GRASS
		// Grass part of avail N is not enough
		if (!GRASS_100 && !not_more_grass) {

			// See if grass can't get more than the 5%
			if (GRASS_part*nsupply_patch>ratio_uptake*grass_uptake_decider) {

				not_more_grass=true;
				// then grass takes GRASS_part of total N supply
				temp_nsupply_patch -= GRASS_part*nsupply_patch;
				// and GRASS strength is subtracted from totaluptake strength
				total_uptake_decider -= grass_uptake_decider;
				// and a new ratio uptake is calculated for TREEs 
				ratio_uptake = temp_nsupply_patch / total_uptake_decider; 
			}
			else {
				// GRASS can compite for more than 5%
				not_more_grass=false;
			}
		}

		vegetation.firstobj();
		while (vegetation.isobj && !full_uptake) {
			Individual& indiv=vegetation.getobj();

			// if lifeform is GRASS and they can't compite with TREEs for more than their part of the total N supply
			if (indiv.pft.lifeform == GRASS && not_more_grass && indiv.fuptake != 1.0) {
				if (!negligible(indiv.ndemand_uptake)) {

					indiv.fuptake = GRASS_part*nsupply_patch*(nitrogen_uptake_strength(indiv)
						/grass_uptake_decider)/indiv.ndemand_uptake;
					if (indiv.fuptake > 1.0)
						indiv.fuptake = 1.0;
				}
				else
					indiv.fuptake = 0.0;
			}

			// if lifeform is TREE and GRASS if it can compete with TREEs
			else {

				// if fuptake does't meet its N demand, then calculate a new value for fuptake
				if (indiv.fuptake != 1.0) {

					// if indiv has the strenght to take up more than N demand
					if (ratio_uptake * nitrogen_uptake_strength(indiv) > indiv.ndemand_uptake && !negligible(indiv.ndemand_uptake)){
						
						indiv.fuptake = 1.0;
						
						// subtract N demand from N supply
						temp_nsupply_patch -= indiv.ndemand_uptake;

						// and take away this indiv uptake strength from total
						total_uptake_decider -= nitrogen_uptake_strength(indiv);

						// and redo indiv fuptake calc for the rest of the indiv as this indiv probably could
						// take up more than its N demand -> more available for the rest of the indiv
						full_uptake = true;
					}
					// normal N limited uptake (0.0 < fuptake < 1.0)
					else if (indiv.ndemand_uptake > 0.0)
						indiv.fuptake = (ratio_uptake * nitrogen_uptake_strength(indiv)) / indiv.ndemand_uptake;
					else
						indiv.fuptake = 0.0;
				}
			}

			vegetation.nextobj();
		}
	}

	// Check so that N uptake matches available N 

	double EPS = 1.0e-15;
	double check_N_uptake = 0.0;

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
		
		check_N_uptake += indiv.fuptake*indiv.ndemand_uptake;

		// Individual fuptake average over 5 years
		if (indiv.alive) {
			int years = 0;
			indiv.fuptake_avr=0.0;

			for (int i=0;i<5;i++){
				if (indiv.fuptake_hist[i]>0.0) {
					indiv.fuptake_avr+=indiv.fuptake_hist[i];
					years++;
				}
				if (i < 4)
					indiv.fuptake_hist[i]=indiv.fuptake_hist[i+1];
				else
					indiv.fuptake_hist[i]=indiv.fuptake*(1.0+indiv.nstorage_uptake);
			}
			indiv.fuptake_avr/=(double)years;
		}
		else {
			indiv.fuptake_avr=indiv.fuptake*(1.0+indiv.nstorage_uptake);
			indiv.fuptake_hist[4]=indiv.fuptake*(1.0+indiv.nstorage_uptake);
		}

		vegetation.nextobj();
	}
}
// end GUESSN

///////////////////////////////////////////////////////////////////////////////
// GUESSN VEGETATION N UPTAKE
///////////////////////////////////////////////////////////////////////////////
void vegetation_n_uptake(Patch& patch) {

	// Daily vegetation uptake of mineral N
	// Partitioned among individuals according to this year's N demand
	// Distributed through the year according to individual daily assimilation

	// Daily N demand given by:
	//	 For individual:
	//     (1)  ndemand_day = dassim/aassim*ndemand_indiv
    //          where
	//          dassim = daily assimilation above base value 0
	//          aassim = annual sum of daily assimilation above base value 0
	//          ndemand_indiv = this year's N demand for growth by this individual 
	//							without any N limitation
	//
	//	 If Extra Nitrogen (EN) is used (ifstoreextran == true), then eq(1) looks like this
	//     (1)  ndemand_day = dassim/aassim*ndemand_indiv*(1+EN_actual_uptake)
	//          where
	//			EN_actual_uptake = the extra nitrogen that will be stored
	//
	//	 For patch:
	//          ndemand_patch_day = sum of ndemand_day over all individuals
	//                        
	// Actual N uptake for each day and individual given by:
	//     (3)  nuptake_day = ndemand_day*fuptake
	//     where
	//     (4)  fuptake = min(nsupply_patch/ndemand_patch,1.0)
	//	   there is also an option for individual fuptake 
	//
	// N deposition and leaching of mineral N are also performed by this function
	// Leaching of organic N is done separately by function leaching() above
	// Uptake to / usage of extra N storage is also done by this function
	//
	// N fixation is done on patch basis using Cleveland 1999 approach
	// Conservative N fixation
	//			N_fix (kgN/ha/yr) = 0.102*ET(cm/yr)+0.524
	//			N_fix (kgN/m2/yr) = (0.102*patch.aaet/10.0+0.524)/10000.0 
	//		(5)	N_fix = 0.00000102*patch.aaet-0.0000524
	//
	// To be called on last day of year following SOM dynamics

	const double EPS=1e-12;

	double ndemand_patch,nsupply_patch,fuptake,fuptake_avr=0.0;
	double dndep,dnmass_avail,dnfix,ndemand_day,nuptake_day;
	double leachn,fleach,excessn;
	double nmass_avail[365]; // daily soil N pool

	Vegetation& vegetation=patch.vegetation;
	Soil& soil=patch.soil;

	ndemand_patch=0.0;

	// ANNUAL N DEMAND FOR PATCH

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();

		indiv.ndemand_uptake = this_years_ndemand(indiv.cmass_leaf,indiv.cmass_root,indiv.cmass_sap,indiv.cmass_heart,
			indiv.cmass_debt,indiv.nmass_leaf,indiv.nmass_root,indiv.nmass_sap,indiv.nmass_heart,indiv.leafn_mean,
			indiv.pft.cton_leaf,indiv.pft.cton_root,indiv.pft.cton_sap,indiv.anpp,indiv.pft.reprfrac,indiv.wscal_mean,indiv.pft.ltor_max,
			indiv.height,indiv.pft.sla,indiv.pft.wooddens,indiv.pft.k_latosa,indiv.pft.k_allom2,indiv.pft.k_allom3,
			indiv.pft.phenology,indiv.aphen_raingreen,indiv.pft.leaflong,
			indiv.pft.turnover_leaf,indiv.pft.turnover_root,indiv.pft.turnover_sap,
			indiv.nstore,indiv.densindiv,indiv.alive,indiv.pft.lifeform);

		//	If to use longterm N storage within an individual
		if (ifnstorage && date.year > freenyears){
			// if N storage is larger than what can be stored then don't store more
			if (indiv.nmass_store > max_nstorage*(indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap))
				indiv.nstorage_uptake=0.0;
			// fill up storage to max
			else if (max_nstorage*(indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap)-indiv.nmass_store < max_nstorage_uptake*indiv.ndemand_uptake)
				indiv.nstorage_uptake=(max_nstorage*(indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap)-indiv.nmass_store)/indiv.ndemand_uptake;
			// store as much as possible 
			else
				indiv.nstorage_uptake=max_nstorage_uptake;

			indiv.ndemand_uptake*=(1.0+indiv.nstorage_uptake);
		}

		// Sum assimilation over period of positive assimilation

		indiv.aassim=0.0;
		for (int d=0;d<365;d++)
			if (indiv.dassim[d]>0.0) indiv.aassim+=indiv.dassim[d];

		if (!negligible(indiv.aassim)) 
			ndemand_patch+=indiv.ndemand_uptake;

		// Check that individual's N store is zero (N in excess of demand should have
		// been returned to soil following allocation last year)

		 vegetation.nextobj();
	}

	// ANNUAL N SUPPLY
	// Potential N supply is remaining pool from last year
	// PLUS annual deposition
	// PLUS estimate of annual N fixation
	// PLUS sum of daily mineralisation MINUS sum of daily immobilisation

	// N deposition
	soil.ndep_annual=patch.stand.climate.andep;

	// N fixation
	if (ifnfix==1)
		soil.N_fix = max(0.00000102*patch.aaet+0.0000524,0.0);	
			// Conservative N fixation (Cleveland 1999 fig. 1)
	else if (ifnfix==2)
		soil.N_fix = max(0.00000234*patch.aaet-0.0000172,0.0);
			// Central N fixation (Cleveland 1999 fig. 1)
	else if (ifnfix==3)
		soil.N_fix = max(0.00000367*patch.aaet-0.0000754,0.0);
			// Upper N fixation (Cleveland 1999 fig. 1)
	else
		soil.N_fix = 0.0;

	// Coarse Woody Debris N fixation	// GUESSNFIX wood
	double cwd_litter = 0.0;

	for (int q=0;q<npft;q++) {
		Patchpft& pft=patch.stand[patch.id].pft[q];

		cwd_litter+=pft.litter_wood;
	}

	cwd_litter+=patch.soil.sompool[SURFCWD].cmass;

	// CWD N Fixation (Brunner and Kimmins 2003)
	soil.cwd_N_fix=max(0.0000165*cwd_litter,0.0);

	//soil.N_fix+=soil.cwd_N_fix;

	// N mineralisation and immobilisation
	soil.nmin_annual=0.0;
	soil.nimmob_annual=0.0;

	for (int d=0;d<365;d++) {
		soil.nmin_annual+=soil.nmin_daily[d];
		soil.nimmob_annual+=soil.nimmob_daily[d];
	}

	// Total N supply in patch
	nsupply_patch=soil.nmass_avail+soil.ndep_annual+soil.N_fix+
		+soil.nmin_annual-soil.nimmob_annual;

	if (patch.id==0){// && !(date.year%5)){
		plot("N fixation (kgN/ha/yr)","Soil N fix",date.year,soil.N_fix*10000.0);
		plot("N fixation (kgN/ha/yr)","CWD N fix",date.year,soil.cwd_N_fix*10000.0);
		//plot("N deposition (kgN/ha/yr)","Ndep",date.year,soil.ndep_annual*10000.0);
		plot("N min-immob (kgN/ha/yr)","N",date.year,(soil.nmin_annual-soil.nimmob_annual)*10000.0);
		//plot("mineral N avail (kgN/ha/yr)","N",date.year,soil.nmass_avail*10000.0);
		plot("N demand/supply (kgN/ha/yr)","N supply",date.year,nsupply_patch*10000.0);
		plot("N demand/supply (kgN/ha/yr)","N demand",date.year,ndemand_patch*10000.0);
		//plot("new establishment N demand (kgN/ha/yr)","N demand",date.year,patch.new_est_ndemand*10000.0);
	}

	// DAILY N SUPPLY

	// Daily N deposition (distributed evenly through the year and doubled to 
	// account for NH4 deposition)
	dndep=soil.ndep_annual/365.0;
	dnfix=soil.N_fix/365.0;
	dnmass_avail=soil.nmass_avail/365.0;

	for (int day=0;day<365;day++)	// Loop through days
		nmass_avail[day]=dnmass_avail+dnfix+dndep+
			soil.nmin_daily[day]-soil.nimmob_daily[day];

	// Rescale demand to not exceed supply (Eqn 4)

	if (nsupply_patch<=0.0) {
		fuptake=0.0;
		patch.nlim=false; 
	}
	else if (ndemand_patch>nsupply_patch) 
	{
		fuptake=nsupply_patch/ndemand_patch;
		patch.nlim=true;	
	}
	else 
	{	
		fuptake=1.0;
		patch.nlim=false;
	}

	patch.fuptake_patch=fuptake;

	// Plot yearly patch N limitation
	if (patch.id==0 && !(date.year%5))
		plot("fuptake","fuptake",date.year,fuptake);

	// Individual fuptake
	if (ifindiv_fuptake && patch.nlim) 
		indiv_fuptake(vegetation,nsupply_patch,fuptake);

	// VEGETATION N UPTAKE
	// Uptake in excess of daily supply permitted

	if (!negligible(ndemand_patch)) { // (some vegetation N uptake this year)
		
		// Loop through days of year

		for (int d=0;d<365;d++) {

			// Loop through individuals

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();

				if (!ifindiv_fuptake || !patch.nlim)
					indiv.fuptake = fuptake;
				
				if (!negligible(indiv.aassim) && indiv.dassim[d]>0.0) {

					// Daily N demand by this individual (Eqn 1)
					
					ndemand_day=indiv.dassim[d]/indiv.aassim*indiv.ndemand_uptake;

					// Daily N uptake by this individual (Eqn 3)
					nuptake_day=ndemand_day*indiv.fuptake;

					// Add to individual's nitrogen stores
					indiv.nstore+=nuptake_day;
					
					// Deduct from soil N pool (negative result allowed)
					nmass_avail[d]-=nuptake_day;
				}

				// ... on to next individual
				vegetation.nextobj();
			}
		}
	}

	// LONG-TERM N STORAGE USAGE/STORING
	if (ifnstorage && date.year > freenyears) {
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();
			
			// if n uptake is larger than N demand (can save N)
			if (indiv.nstore > indiv.ndemand_uptake/(1.0+indiv.nstorage_uptake)) {
					indiv.nmass_store+=indiv.nstore-indiv.ndemand_uptake/(1.0+indiv.nstorage_uptake);
					indiv.nstore-=indiv.nstore-indiv.ndemand_uptake/(1.0+indiv.nstorage_uptake);
			}

			// how this year compares to the average of the last five years
			double diff = indiv.fuptake*(1.0+indiv.nstorage_uptake) - indiv.fuptake_avr;

			if (diff > 0.0 && indiv.fuptake_avr < 1.0) { // Save some N as this year is better than the average (50%)
				if (indiv.fuptake*(1.0+indiv.nstorage_uptake) > 1.0)
					diff -= (indiv.fuptake*(1.0+indiv.nstorage_uptake)-1.0);

				indiv.nmass_store += indiv.nstore*0.5*diff;
				indiv.nstore -= indiv.nstore*0.5*diff;
			}
			else if (diff < 0.0 && !negligible(indiv.nmass_store) && indiv.fuptake*(1.0+indiv.nstorage_uptake) < 1.0){	
				// Try to fill up nstore with 70% of the difference from the average with N from nmass_store
				indiv.nstore += min(0.7,(indiv.ndemand_uptake*-diff)/indiv.nmass_store)*indiv.nmass_store;
				indiv.nmass_store -= min(0.7,(indiv.ndemand_uptake*-diff)/indiv.nmass_store)*indiv.nmass_store;
			}

			// ... on to next individual
			vegetation.nextobj();
		}
	}

	// LEACHING OF SOIL MINERAL N
	// Allowed on days with residual N following vegetation uptake
	// Daily leaching fractions were pre-computed by function leaching() above

	soil.nleach_annual=0.0;

	excessn=nsupply_patch-ndemand_patch*fuptake;
	double save_excessn=excessn;

	if (excessn>0.0) {

		double nmass_sum=0.0;

		for (int d=0;d<365;d++) {
			nmass_sum+=nmass_avail[d];
			if (nmass_sum>0.0)  {
				leachn=nmass_sum*soil.leachfrac_daily[d];

				if (soil.nleach_annual+leachn <= excessn) {
					nmass_sum-=leachn;
					soil.sompool[LEACHED].nmass+=leachn;
					soil.nleach_annual+=leachn;
					nmass_avail[d]-=leachn;
				}
			}
		}	
	}

	// EXCESS MINERAL N
	// Return remaining N to soil store for next year

	excessn=0.0;
	for (int days=0;days<365;days++)
		excessn+=nmass_avail[days];

	// Should never be negative! (allow it for very small values for now ...)
	//if (excessn<-EPS)
	//	dprintf("Year %d vegetation_n_uptake: patch %d Unexpected NEGATIVE value (%g) for annual excess mineral N before leach (%g)\n",
	//		date.year,patch.id,excessn,nsupply_patch-ndemand_patch*fuptake);

	soil.nmass_avail=excessn;
}
// end GUESSN


// GUESSN
// Temporary (for checking) NB! only works with one patch
double old_total=0.0;
double old_vegn=0.0;
double old_vegstore=0.0;
double old_vegnmass_store=0.0;
double old_centuryn=0.0;
double old_nmass_avail=0.0;
double old_littern=0.0;
double old_leachn=0.0;

void som_dynamics_century(Patch& patch) {

	double vegn,centuryn,littern,totaln,vegstore;
	double sumassim;
	double vegnmass_store;	
	int p;

	if (date.day==0) { // First day of year only

		Soil& soil=patch.soil;
		Vegetation& vegetation=patch.vegetation;

		if (patch.id==0) {

			// Work out total ecosystem N for checking

			vegn=0.0;
			vegstore=0.0;
			vegnmass_store=0.0;	
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv=vegetation.getobj();
				if (indiv.alive) {	
					vegn+=indiv.nmass_leaf+indiv.nmass_root+indiv.nmass_sap+indiv.nmass_heart;
					vegstore+=indiv.nstore;
					vegnmass_store+=indiv.nmass_store; 
				}

				vegetation.nextobj();
			}

			centuryn=0.0;
			for (p=0;p<NSOMPOOL;p++) {
				centuryn+=soil.sompool[p].nmass;
			}	

			littern=0.0;
			patch.pft.firstobj();
			while (patch.pft.isobj) {
				Patchpft& pft=patch.pft.getobj();

				littern+=pft.nmass_litter_leaf+
					pft.nmass_litter_root+
					pft.nmass_litter_wood;

				patch.pft.nextobj();
			}
			// Thomas monitoring
			if (date.year>500) {

				plot("Total N","vegn",date.year,vegn);
				plot("Total N","vegstore",date.year,vegstore);
				plot("Total N","vegnmass_store",date.year,vegnmass_store);
				plot("Total N","SOM",date.year,centuryn);
				plot("Total N","mineral N",date.year,soil.nmass_avail);
			//	plot("Actual mineral N after excessn","mineral N",date.year,soil.nmass_avail);
			//	plot("Actual leached N after excessn","leached N",date.year,soil.nleach_annual);
				plot("Total N","litter N",date.year,littern);
				plot("Total N","leached N",date.year,soil.sompool[LEACHED].nmass);
				plot("Total N","total",date.year,
					(vegn+vegstore+vegnmass_store+centuryn+soil.nmass_avail+littern+soil.sompool[LEACHED].nmass));

				plot("deltaN","vegn",date.year,vegn-old_vegn);
				plot("deltaN","vegstore",date.year,vegstore-old_vegstore);
				plot("deltaN","vegnmass_store",date.year,vegnmass_store-old_vegnmass_store);
				plot("deltaN","SOM",date.year,centuryn-old_centuryn);
				plot("deltaN","mineral N",date.year,soil.nmass_avail-old_nmass_avail);
				plot("deltaN","litter N",date.year,littern-old_littern);
				plot("deltaN","leached N",date.year,soil.sompool[LEACHED].nmass-old_leachn);
				plot("deltaN","total",date.year,(vegn+centuryn+soil.nmass_avail+vegstore+littern+soil.sompool[LEACHED].nmass)-old_total);
			}

			old_vegn=vegn;
			old_vegstore=vegstore;
			old_vegnmass_store=vegnmass_store;
			old_centuryn=centuryn;
			old_nmass_avail=soil.nmass_avail;
			old_littern=littern;
			old_leachn=soil.sompool[LEACHED].nmass;

			old_total=vegn+centuryn+soil.nmass_avail+vegstore+littern+soil.sompool[LEACHED].nmass;
		}

		// Transfer last year's litter to SOM pools
		transfer_litter(patch,patch.soil);
	}

	if (date.dayofmonth==0) patch.fluxes.mcflux_soil[date.month]=0.0;

	// Daily or monthly decomposition and fluxes between SOM pools
	
	somfluxes(patch,patch.soil,patch.fluxes);	

	// Leaching of organic N/C and potential daily leaching fraction for mineral N

	leaching(patch.soil);

	if (date.islastmonth && date.islastday) {

		// Last day of year

		// Distribute plant N uptake and leaching of mineral nitrogen throughout the past year
		// Calculate mineral N pool at end of year
		
		vegetation_n_uptake(patch);	
	}

	Soil& soil=patch.soil;
	if (date.day==0 && (date.year%10==0) && patch.id==0 && ifcentury) {
		
		plot("century C","surfstruct",date.year,soil.sompool[SURFSTRUCT].cmass);
		plot("century C","surfmeta",date.year,soil.sompool[SURFMETA].cmass);
		plot("century C","surfcwd",date.year,soil.sompool[SURFCWD].cmass);
		plot("century C","surfmicro",date.year,soil.sompool[SURFMICRO].cmass);
		plot("century C","soilstruct",date.year,soil.sompool[SOILSTRUCT].cmass);
		plot("century C","soilmeta",date.year,soil.sompool[SOILMETA].cmass);		
		plot("century C","soilmicro",date.year,soil.sompool[SOILMICRO].cmass);
		plot("century C","humussom",date.year,soil.sompool[SURFHUMUS].cmass);
		plot("century C","slowsom",date.year,soil.sompool[SLOWSOM].cmass);
		plot("century C","passivesom",date.year,soil.sompool[PASSIVESOM].cmass); 
 
		plot("century N","surfstruct",date.year,soil.sompool[SURFSTRUCT].nmass);
		plot("century N","surfmeta",date.year,soil.sompool[SURFMETA].nmass);
		plot("century N","surfcwd",date.year,soil.sompool[SURFCWD].nmass);
		plot("century N","surfmicro",date.year,soil.sompool[SURFMICRO].nmass);
		plot("century N","soilstruct",date.year,soil.sompool[SOILSTRUCT].nmass);
		plot("century N","soilmeta",date.year,soil.sompool[SOILMETA].nmass);
		plot("century N","soilmicro",date.year,soil.sompool[SOILMICRO].nmass);
		plot("century N","humussom",date.year,soil.sompool[SURFHUMUS].nmass);
		plot("century N","slowsom",date.year,soil.sompool[SLOWSOM].nmass);
		plot("century N","passivesom",date.year,soil.sompool[PASSIVESOM].nmass);
	}

}
// end GUESSN

void som_dynamics(Patch& patch) {

	// Choose between CENTURY or standard LPJ SOM dynamics

	if (ifcentury) som_dynamics_century(patch);
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
// Kirschbaum, M. U. F. and K. I. Paul (2002). "Modelling C and N dynamics in forest soils 
//   with a modified version of the CENTURY model." Soil Biology & Biochemistry 34(3): 341-354.
// Meentemeyer, V. (1978) Macroclimate and lignin control of litter decomposition
//   rates. Ecology 59: 465-472.
// Parton (2010) ForCent model development and testing using the Enriched Background 
//	 Isotope Study experiment JoGR 115: 
