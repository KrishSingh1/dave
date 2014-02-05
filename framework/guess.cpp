///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.cpp
/// \brief LPJ-GUESS Combined Modular Framework
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
int npft; // number of possible PFTs

Pftlist pftlist;

// emission ratios from fire (NH3, NO, NO2, N2O) Delmas et al. 1995

const double Fluxes::NH3_FIRERATIO = 0.014;
const double Fluxes::NO_FIRERATIO  = 0.531;
const double Fluxes::NO2_FIRERATIO = 0.379;
const double Fluxes::N2O_FIRERATIO = 0.076;


////////////////////////////////////////////////////////////////////////////////
// Implementation of PhotosynthesisResult member functions
////////////////////////////////////////////////////////////////////////////////


void PhotosynthesisResult::serialize(ArchiveStream& arch) {
	arch & agd_g
		& adtmm
		& rd_g
		& vm
		& je
		& nactive_opt
		& vmaxnlim;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Climate member functions
////////////////////////////////////////////////////////////////////////////////


void Climate::serialize(ArchiveStream& arch) {
	arch & temp
		& rad
		& par
		& prec
		& daylength
		& co2
		& lat
		& insol
		& instype
		& eet
		& mtemp
		& mtemp_min20
		& mtemp_max20
		& mtemp_max
		& gdd5
		& agdd5 
		& chilldays
		& ifsensechill
		& gtemp
		& dtemp_31
		& mtemp_min_20
		& mtemp_max_20
		& mtemp_min
		& atemp_mean
		& sinelat
		& cosinelat
		& qo & u & v & hh & sinehh
		& daylength_save
		& doneday
		& andep
		& dndep
		& anfert
		& dnfert;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Fluxes member functions
////////////////////////////////////////////////////////////////////////////////

Fluxes::Fluxes(Patch& p) 		
  : patch(p), 
    annual_fluxes_per_pft(npft, std::vector<double>(NPERPFTFLUXTYPES)) {
	
	reset();
}

void Fluxes::reset() {
	for (size_t i = 0; i < annual_fluxes_per_pft.size(); ++i) {
		std::fill_n(annual_fluxes_per_pft[i].begin(), int(NPERPFTFLUXTYPES), 0);
	}

	for (int m = 0; m < 12; ++m) {
		std::fill_n(monthly_fluxes_pft[m], int(NPERPFTFLUXTYPES), 0);

		std::fill_n(monthly_fluxes_patch[m], int(NPERPATCHFLUXTYPES), 0);
	}
}

void Fluxes::serialize(ArchiveStream& arch) {
	arch & annual_fluxes_per_pft 
		& monthly_fluxes_patch
		& monthly_fluxes_pft;
}

void Fluxes::report_flux(PerPFTFluxType flux_type, int pft_id, double value) {
	annual_fluxes_per_pft[pft_id][flux_type] += value;
	monthly_fluxes_pft[date.month][flux_type] += value;
}

void Fluxes::report_flux(PerPatchFluxType flux_type, double value) {
	monthly_fluxes_patch[date.month][flux_type] += value;
}

double Fluxes::get_monthly_flux(PerPFTFluxType flux_type, int month) const {
	return monthly_fluxes_pft[month][flux_type];
}

double Fluxes::get_monthly_flux(PerPatchFluxType flux_type, int month) const {
	return monthly_fluxes_patch[month][flux_type];
}

double Fluxes::get_annual_flux(PerPFTFluxType flux_type, int pft_id) const {
	return annual_fluxes_per_pft[pft_id][flux_type];
}

double Fluxes::get_annual_flux(PerPFTFluxType flux_type) const {
	double sum = 0;
	for (size_t i = 0; i < annual_fluxes_per_pft.size(); ++i) {
		sum += annual_fluxes_per_pft[i][flux_type];
	}
	return sum;
}

double Fluxes::get_annual_flux(PerPatchFluxType flux_type) const {
	double sum = 0;
	for (int m = 0; m < 12; ++m) {
		sum += monthly_fluxes_patch[m][flux_type];
	}
	return sum;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Vegetation member functions
////////////////////////////////////////////////////////////////////////////////


void Vegetation::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		arch & nobj;

		for (unsigned int i = 0; i < nobj; i++) {
			Individual& indiv = (*this)[i];
			arch & indiv.pft.id
				& indiv;
		}
	}
	else {
		killall();
		unsigned int number_of_individuals;
		arch & number_of_individuals;

		for (unsigned int i = 0; i < number_of_individuals; i++) {
			int pft_id;
			arch & pft_id;
			Individual& indiv = createobj(pftlist[pft_id], *this);
			arch & indiv;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of LitterSolveSOM member functions
////////////////////////////////////////////////////////////////////////////////


void LitterSolveSOM::serialize(ArchiveStream& arch) {
	arch & clitter
		& nlitter;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Soil member functions
////////////////////////////////////////////////////////////////////////////////


void Soil::serialize(ArchiveStream& arch) {
	arch & wcont
		& awcont
		& wcont_evap
		& dwcontupper
		& mwcontupper
		& snowpack
		& runoff
		& temp
		& dtemp
		& mtemp
		& gtemp
		& cpool_slow
		& cpool_fast
		& decomp_litter_mean
		& k_soilfast_mean
		& k_soilslow_mean
		& alag
		& exp_alag
		& mwcont
		& dwcontlower
		& mwcontlower
		// probably shouldn't need to serialize these
		& rain_melt
		& max_rain_melt
		& percolate;

	for (int i = 0; i<NSOMPOOL; i++) {
		arch & sompool[i];
	} 

	arch & dperc		
		& orgleachfrac
		& nmass_avail	
		& ninput
		& anmin			
		& animmob			
		& aminleach		
		& aorgleach					
		& anfix
		& anfix_calc
		& anfix_mean
		& snowpack_nmass
		& solvesomcent_beginyr
		& solvesomcent_endyr
		& solvesom
		& fnuptake_mean
		& morgleach_mean
		& mminleach_mean; 
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patchpft member functions
////////////////////////////////////////////////////////////////////////////////


void Patchpft::serialize(ArchiveStream& arch) {
	arch & anetps_ff
		& wscal
		& wscal_mean
		& anetps_ff_est
		& anetps_ff_est_initial
		& wscal_mean_est
		& phen
		& aphen
		& establish
		& nsapling
		& litter_leaf
		& litter_root
		& litter_sap
		& litter_heart
		& litter_repr
		& gcbase
		& gcbase_day
		& wsupply
		& wsupply_leafon
		& fwuptake
		& wstress
		& wstress_day
		& harvested_products_slow
		& nmass_litter_leaf
		& nmass_litter_root
		& nmass_litter_sap
		& nmass_litter_heart
		& harvested_products_slow_nmass;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Patch member functions
////////////////////////////////////////////////////////////////////////////////


void Patch::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}
	}
	else {
		pft.killall();
				
		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			pft.createobj(pftlist[i]);
			arch & pft[i];
		}
	}

	arch & vegetation
		& soil
		& fluxes
		& fpar_grass
		& fpar_ff
		& par_grass_mean
		& nday_growingseason
		& fpc_total
		& disturbed
		& age
		& fireprob
		& growingseasondays
		& intercep
		& aaet
		& aaet_5
		& aevap
		& aintercep
		& arunoff
		& apet
		& eet_net_veg
		& wdemand
		& wdemand_day
		& wdemand_leafon
		& fpc_rescale
		& maet
		& mevap
		& mintercep
		& mrunoff
		& mpet
		& ndemand;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Standpft member functions
////////////////////////////////////////////////////////////////////////////////


void Standpft::serialize(ArchiveStream& arch) {
	arch & cmass_repr
		& anetps_ff_max
		& gpterm
		& fpc_total
		& active;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Stand member functions
////////////////////////////////////////////////////////////////////////////////

Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	unsigned int p;
	unsigned int npatchL = 1;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}


	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND) {
		npatchL=1;
	}
	else if(landcover==NATURAL || landcover==FOREST) {
		npatchL=::npatch; // use the global variable npatch (not Stand::npatch)
	}
	else {
		// Someone has added a new landcover type, the code above needs to be updated and
		// npatchL needs to be set properly.
		fail("Unrecognized landcover type");
	}

	for (p=0;p<npatchL;p++) {
		createobj(*this,gc.soiltype);
	}

	first_year=date.year;
	seed = 12345678;
}

double Stand::get_gridcell_fraction() const {
	return frac*gridcell.landcoverfrac[landcover];
}

double Stand::get_landcover_fraction() const {
	return frac;
}

void Stand::set_landcover_fraction(double fraction) {
	frac = fraction;
}

void Stand::serialize(ArchiveStream& arch) {
	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}

		arch & nobj;
		for (unsigned int k = 0; k < nobj; k++) {
			arch & (*this)[k];
		}
	}
	else {
		pft.killall();
		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			Standpft& standpft = pft.createobj(pftlist[i]);
			arch & standpft;
		}

		killall();
		unsigned int npatch;
		arch & npatch;
		for (unsigned int k = 0; k < npatch; k++) {
			Patch& patch = createobj(*this, gridcell.soiltype);
			arch & patch;
		}
	}

	arch & first_year
		& frac
		& seed;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Individual member functions
////////////////////////////////////////////////////////////////////////////////

Individual::Individual(int i,Pft& p,Vegetation& v):pft(p),vegetation(v),id(i) {

	anpp              = 0.0;
	fpc               = 0.0;
	densindiv         = 0.0;
	cmass_leaf        = 0.0;
	cmass_root        = 0.0;
	cmass_sap         = 0.0;
	cmass_heart       = 0.0;
	cmass_debt        = 0.0;
	phen              = 0.0;
	aphen             = 0.0;
	deltafpc          = 0.0;

	nmass_leaf        = 0.0;
	nmass_root        = 0.0;
	nmass_sap         = 0.0;
	nmass_heart       = 0.0;
	cton_leaf_aopt    = 0.0;
	cton_leaf_aavr    = 0.0;
	cton_status       = 0.0;
	cmass_veg         = 0.0;
	nmass_veg         = 0.0;

	nactive           = 0.0;
	nextin            = 1.0;
	nstore_longterm   = 0.0;
	nstore_labile     = 0.0;
	ndemand           = 0.0;
	fnuptake          = 1.0;
	anuptake          = 0.0;
	max_n_storage     = 0.0;
	scale_n_storage   = 0.0;

	leafndemand       = 0.0;
	rootndemand       = 0.0;
	sapndemand        = 0.0;
	storendemand      = 0.0;
	leaffndemand      = 0.0;
	rootfndemand      = 0.0;
	sapfndemand       = 0.0;
	storefndemand     = 0.0;
	leafndemand_store = 0.0;
	rootndemand_store = 0.0;

	nstress           = false;

	// additional initialisation
	age               = 0.0;
	fpar              = 0.0;
	aphen_raingreen   = 0;
	intercep          = 0.0;
	phen_mean         = 0.0;
	wstress           = false;
	lai               = 0.0;
	lai_layer         = 0.0;
	lai_indiv         = 0.0;
	alive             = false;

	int m;
	for (m=0; m<12; m++) {
		mlai[m] = 0.0;
	}

	// bvoc
	monstor           = 0.;
	iso               = 0.;
	mon               = 0.;
	fvocseas          = 1.;
}

void Individual::serialize(ArchiveStream& arch) {
	arch & cmass_leaf
		& cmass_root
		& cmass_sap 
		& cmass_heart
		& cmass_debt
		& fpc
		& fpar
		& densindiv
		& phen
		& aphen
		& aphen_raingreen
		& anpp
		& aet
		& ltor
		& height
		& crownarea
		& deltafpc
		& wscal_mean
		& boleht
		& lai
		& lai_layer
		& lai_indiv
		& greff_5
		& age
		& mlai
		& fpar_leafon
		& lai_leafon_layer
		& intercep
		& phen_mean
		& wstress 
		& alive 
		& iso 
		& mon 
		& monstor 
		& fvocseas 
		& nmass_leaf
		& nmass_root
		& nmass_sap
		& nmass_heart
		& nactive
		& nextin
		& nstore_longterm
		& nstore_labile
		& ndemand
		& fnuptake
		& anuptake
		& max_n_storage
		& scale_n_storage
		& avmaxnlim
		& cton_leaf_aopt
		& cton_leaf_aavr
		& cton_status
		& cmass_veg
		& nmass_veg

		& nstress
		& leafndemand
		& rootndemand
		& sapndemand
		& storendemand
		& leaffndemand
		& rootfndemand
		& sapfndemand
		& storefndemand
		& leafndemand_store
		& rootndemand_store
		
		& nday_leafon;
}

void Individual::report_flux(Fluxes::PerPFTFluxType flux_type, double value) {
	if (alive) {
		vegetation.patch.fluxes.report_flux(flux_type, pft.id, value);
	}
}

void Individual::report_flux(Fluxes::PerPatchFluxType flux_type, double value) {
	if (alive) {
		vegetation.patch.fluxes.report_flux(flux_type, value);
	}
}

void Individual::reduce_biomass(double mortality, double mortality_fire) {

	// This function needs to be modified if a new lifeform is added,
	// specifically to deal with nstore().
	assert(pft.lifeform == TREE || pft.lifeform == GRASS);

	const double mortality_non_fire = mortality - mortality_fire;

	// Transfer killed biomass to litter
	// (above-ground biomass killed by fire enters atmosphere, not litter)

	Patchpft& ppft = patchpft();

	ppft.litter_leaf  += mortality_non_fire * cmass_leaf;
	ppft.litter_root  += mortality * cmass_root;

	if (cmass_debt <= cmass_heart + cmass_sap) {
		if (cmass_debt <= cmass_heart) {
			ppft.litter_sap   += mortality_non_fire * cmass_sap;
			ppft.litter_heart += mortality_non_fire * (cmass_heart - cmass_debt);
		}
		else {
			ppft.litter_sap   += mortality_non_fire * (cmass_sap + cmass_heart - cmass_debt);
		}
	}
	else {
		double debt_excess = mortality_non_fire * (cmass_debt - (cmass_sap + cmass_heart));
		report_flux(Fluxes::NPP, debt_excess);
		report_flux(Fluxes::RA, -debt_excess);
	}
		
	ppft.nmass_litter_leaf  += mortality_non_fire * nmass_leaf;
	ppft.nmass_litter_root  += mortality * nmass_root;
	ppft.nmass_litter_sap   += mortality_non_fire * nmass_sap;
	ppft.nmass_litter_heart += mortality_non_fire * nmass_heart;

	if (pft.lifeform == TREE) {				
		// Transfer nitrogen storage to wood nitrogen litter for now 	
		ppft.nmass_litter_sap += mortality_non_fire * nstore();
	}
	else { // GRASS
		// Transfer nitrogen storage to root nitrogen litter for now
		ppft.nmass_litter_root += mortality * nstore();			
	}


	// Flux to atmosphere from burnt above-ground biomass

	double cflux_fire = mortality_fire * (cmass_leaf + cmass_wood());
	double nflux_fire = mortality_fire * (nmass_leaf + nmass_wood());

	if (pft.lifeform == TREE) {
		nflux_fire += mortality_fire * nstore();
	}

	report_flux(Fluxes::FIREC,    cflux_fire);

	report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire); 
	report_flux(Fluxes::NO_FIRE,  Fluxes::NO_FIRERATIO  * nflux_fire);
	report_flux(Fluxes::NO2_FIRE, Fluxes::NO2_FIRERATIO * nflux_fire);
	report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);

	// Reduce this Individual's biomass values

	const double remaining = 1.0 - mortality;

	if (pft.lifeform != GRASS) {
		densindiv *= remaining;
	}

	cmass_leaf      *= remaining;
	cmass_root      *= remaining;
	cmass_sap       *= remaining;
	cmass_heart     *= remaining;
	cmass_debt      *= remaining;
	nmass_leaf      *= remaining;
	nmass_root      *= remaining;
	nmass_sap       *= remaining;
	nmass_heart     *= remaining;
	nstore_longterm *= remaining;
	nstore_labile   *= remaining;
}

double Individual::cton_leaf(bool use_phen /* = true*/) const {
	if (ifnlim) {
		if (!negligible(cmass_leaf) && !negligible(nmass_leaf)) {
			if (use_phen) {
				if (!negligible(phen)) {
					return cmass_leaf * phen / nmass_leaf;
				}
				else {
					return pft.cton_leaf_avr;
				}
			}
			else {
				return cmass_leaf / nmass_leaf;
			}
		}
		else {
			return pft.cton_leaf_max;
		}
	}
	else {
		return pft.cton_leaf_avr;
	}
}

double Individual::cton_root(bool use_phen /* = true*/) const {
	if (ifnlim) {
		if (!negligible(cmass_root) && !negligible(nmass_root)) { 
			if (use_phen) {
				if (!negligible(phen)) {
					return cmass_root * phen / nmass_root;
				}
				else {
					return pft.cton_root_avr;
				}
			}
			else {
				return cmass_root / nmass_root;
			}
		}
		else {
			return pft.cton_root_max;
		}
	}
	else {
		return pft.cton_root_avr;
	}
}

double Individual::cton_sap() const {
	if (pft.lifeform == TREE) {
		if (ifnlim) {
			if (!negligible(cmass_sap) && !negligible(nmass_sap))
				return cmass_sap / nmass_sap;
			else
				return pft.cton_sap_max;
		}
		else {
			return pft.cton_sap_avr;
		}
	}
	else {
		return 1.0;
	}
}


Patchpft& Individual::patchpft() const {
	return vegetation.patch.pft[pft.id];
}

/// Help function for kill(), partitions wood biomass into litter and harvest
/** 
 *  Wood biomass (either C or N) is partitioned into litter pools and
 *  harvest, according to PFT specific harvest fractions.
 *
 *  Biomass is sent in as sap and heart, any debt should already have been
 *  subtracted from these before calling this function.
 *
 *  \param mass_sap          Sapwood
 *  \param mass_heart        Heartwood
 *  \param harv_eff          Harvest efficiency (fraction of biomass harvested)
 *  \param harvest_slow_frac Fraction of harvested products that goes into slow depository
 *  \param res_outtake       Fraction of residue outtake at harvest
 *  \param litter_sap        Biomass going to sapwood litter pool
 *  \param litter_heart      Biomass going to heartwood litter pool
 *  \param fast_harvest      Biomass going to harvest flux
 *  \param slow_harvest      Biomass going to slow depository
 */
void partition_wood_biomass(double mass_sap, double mass_heart,
                            double harv_eff, double harvest_slow_frac, double res_outtake,
                            double& litter_sap, double& litter_heart,
                            double& fast_harvest, double& slow_harvest) {

	double sap_left = mass_sap;
	double heart_left = mass_heart;

	// Remove harvest
	double total_wood_harvest = harv_eff * (sap_left + heart_left);

	sap_left   *= 1 - harv_eff;
	heart_left *= 1 - harv_eff;

	// Partition wood harvest into slow and fast
	slow_harvest = total_wood_harvest * harvest_slow_frac;
	fast_harvest = total_wood_harvest * (1 - harvest_slow_frac);

	// Remove residue outtake
	fast_harvest += res_outtake * (sap_left + heart_left);
				
	sap_left   *= 1 - res_outtake;
	heart_left *= 1 - res_outtake;

	// The rest goes to litter
	litter_sap   = sap_left;
	litter_heart = heart_left;
}


void Individual::kill(bool harvest /* = false */) {
	Patchpft& ppft = patchpft();

	double charvest_flux = 0.0;
	double charvested_products_slow = 0.0;

	double nharvest_flux = 0.0;
	double nharvested_products_slow = 0.0;

	double harv_eff = 0.0;
	double harvest_slow_frac = 0.0;
	double res_outtake = 0.0;

	// The function always deals with harvest, but the harvest
	// fractions are zero when there is no harvest.
	if (harvest) {
		harv_eff = pft.harv_eff;

		if (ifslowharvestpool) {
			harvest_slow_frac = pft.harvest_slow_frac;
		}
		
		res_outtake = pft.res_outtake;
	}

	// C doesn't return to litter/harvest if the Individual isn't alive
	if (alive) {

		// For leaf and root, catches small, negative values too

		// Leaf: remove residue outtake and send the rest to litter
		ppft.litter_leaf += cmass_leaf * (1 - res_outtake);
		charvest_flux    += cmass_leaf * res_outtake;

		// Root: all goes to litter
		ppft.litter_root += cmass_root;

		// Deal with the wood biomass and carbon debt for trees
		if (pft.lifeform == TREE) {

			// debt smaller than existing wood biomass
			if (cmass_debt <= cmass_sap + cmass_heart) {

				// before partitioning the biomass into litter and harvest,
				// first get rid of the debt so we're left with only
				// sap and heart
				double to_partition_sap   = 0.0;
				double to_partition_heart = 0.0;

				if (cmass_heart >= cmass_debt) {
					to_partition_sap   = cmass_sap;
					to_partition_heart = cmass_heart - cmass_debt;
				}
				else {
					to_partition_sap   = cmass_sap + cmass_heart - cmass_debt;
				}

				double clitter_sap, clitter_heart, cwood_harvest;

				partition_wood_biomass(to_partition_sap, to_partition_heart,
				                       harv_eff, harvest_slow_frac, res_outtake,
				                       clitter_sap, clitter_heart,
				                       cwood_harvest, charvested_products_slow);
				
				ppft.litter_sap   += clitter_sap;
				ppft.litter_heart += clitter_heart;

				charvest_flux += cwood_harvest;
			}
			// debt larger than existing wood biomass
			else {
				double debt_excess = cmass_debt - (cmass_sap + cmass_heart);
				report_flux(Fluxes::NPP, debt_excess);
				report_flux(Fluxes::RA, -debt_excess);
			}
		}
	}

	// Nitrogen always return to soil litter
	if (pft.lifeform == TREE) {

		double nlitter_sap, nlitter_heart, nwood_harvest;

		// Transfer nitrogen storage to sapwood nitrogen litter/harvest
		partition_wood_biomass(nmass_sap + nstore(), nmass_heart,
		                       harv_eff, harvest_slow_frac, res_outtake,
		                       nlitter_sap, nlitter_heart,
		                       nwood_harvest, nharvested_products_slow);

		ppft.nmass_litter_sap   += nlitter_sap;
		ppft.nmass_litter_heart += nlitter_heart;

		nharvest_flux += nwood_harvest;
	}
	else {
		// Transfer nitrogen storage to root nitrogen litter
		ppft.nmass_litter_root += nstore();
	}

	// Leaf: remove residue outtake and send the rest to litter
	ppft.nmass_litter_leaf += nmass_leaf * (1 - res_outtake);
	nharvest_flux          += nmass_leaf * res_outtake;

	// Root: all goes to litter
	ppft.nmass_litter_root += nmass_root;

	// Report harvest fluxes
	report_flux(Fluxes::HARVESTC, charvest_flux);
	report_flux(Fluxes::HARVESTN, nharvest_flux);

	// Add to biomass depositories for long-lived products
	ppft.harvested_products_slow += charvested_products_slow;
	ppft.harvested_products_slow_nmass += nharvested_products_slow;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcellpft member functions
////////////////////////////////////////////////////////////////////////////////


void Gridcellpft::serialize(ArchiveStream& arch) {
	arch & addtw
		 & Km;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of Gridcell member functions
////////////////////////////////////////////////////////////////////////////////

double Gridcell::get_lon() const {
	return lon;
}

double Gridcell::get_lat() const {
	return lat;
}

void Gridcell::set_coordinates(double longitude, double latitude) {
	lon = longitude;
	lat = latitude;
}

void Gridcell::serialize(ArchiveStream& arch) {
	arch & climate
		& landcoverfrac
		& landcoverfrac_old
		& LC_updated
		& seed;

	if (arch.save()) {
		for (unsigned int i = 0; i < pft.nobj; i++) {
			arch & pft[i];
		}

		arch & nobj;
		for (unsigned int s = 0; s < nobj; s++) {
			arch & (*this)[s].landcover
				& (*this)[s];
		}
	}
	else {
		pft.killall();

		for (unsigned int i = 0; i < pftlist.nobj; i++) {
			pft.createobj(pftlist[i]);
			arch & pft[i];
		}

		killall();
		unsigned int number_of_stands;
		arch & number_of_stands;
				
		for (unsigned int s = 0; s < number_of_stands; s++) {
			landcovertype landcover;
			arch & landcover;
			createobj(*this, landcover);
			arch & (*this)[s];
		}
	}
}

void Sompool::serialize(ArchiveStream& arch) {
	arch & cmass
		& nmass
		& cdec 
		& ndec 
		& delta_cmass
		& delta_nmass
		& ligcfrac
		& fracremain
		& ntoc
		& litterme
		& fireresist
		& mfracremain_mean;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995. 
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991 
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
