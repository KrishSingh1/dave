///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// \author Mats Lindeskog
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H
#include "guess.h"
#include "growth.h"
#include "canexch.h"
#include "cropsowing.h"
#include "cropphenology.h"
#include "cropallocation.h"
#include "management.h"
#include "externalinput.h"
#include "inputmodule.h"

const int def_verndate_ndays_after_last_springdate = 60;

///	Creates stands for landcovers present in the gridcell first year of the simulation
void landcover_init(Gridcell& gridcell, InputModule* input_module);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module);

/// Query whether a date is within a period spanned by two dates.
bool dayinperiod(int day, int start, int end);
/// Step n days from a date.
int stepfromdate(int day, int step);

/// struct storing carbon, nitrogen and water (plus aaet_5 and anfix_calc) during landcover change
struct landcover_change_transfer {

	double *transfer_litter_leaf;
	double *transfer_litter_sap;
	double *transfer_litter_heart;
	double *transfer_litter_root;
	double *transfer_litter_repr;
	double *transfer_harvested_products_slow;
	double *transfer_nmass_litter_leaf;
	double *transfer_nmass_litter_sap;
	double *transfer_nmass_litter_heart;
	double *transfer_nmass_litter_root;
	double *transfer_harvested_products_slow_nmass;

	double transfer_acflux_harvest;
	double transfer_anflux_harvest;

	double transfer_cpool_fast[NSOILLAYER_MAX];
	double transfer_cpool_slow[NSOILLAYER_MAX];
	double transfer_wcont[NSOILLAYER_MAX];
	double transfer_wcont_evap;
	double transfer_decomp_litter_mean[NSOILLAYER_MAX];
	double transfer_k_soilfast_mean[NSOILLAYER_MAX];
	double transfer_k_soilslow_mean[NSOILLAYER_MAX];
	Sompool transfer_sompool[NSOMPOOL][NSOILLAYER_MAX];
	double transfer_NH4_mass[NSOILLAYER_MAX];
	double transfer_NO3_mass[NSOILLAYER_MAX];
	double transfer_snowpack[NLAYERS_SNOW_MAX];
	double transfer_snowpack_NH4_mass[NLAYERS_SNOW_MAX];
	double transfer_snowpack_NO3_mass[NLAYERS_SNOW_MAX];

	double transfer_N2_mass[NSOILLAYER_MAX];
	double transfer_N2O_mass[NSOILLAYER_MAX];
	double transfer_NO_mass[NSOILLAYER_MAX];
	double transfer_NO2_mass[NSOILLAYER_MAX];

	double transfer_NH4_mass_w[NSOILLAYER_MAX];
	double transfer_NO3_mass_w[NSOILLAYER_MAX];
	double transfer_NH4_mass_d[NSOILLAYER_MAX];
	double transfer_NO3_mass_d[NSOILLAYER_MAX];
	double transfer_N2O_mass_w[NSOILLAYER_MAX];
	double transfer_NO_mass_w[NSOILLAYER_MAX];
	double transfer_NO2_mass_w[NSOILLAYER_MAX];
	double transfer_N2O_mass_d[NSOILLAYER_MAX];
	double transfer_NO_mass_d[NSOILLAYER_MAX];
	double transfer_NO2_mass_d[NSOILLAYER_MAX];
	Historic<double, NYEARAAET> transfer_aaet_5;
	double transfer_anfix_calc;

	landcover_change_transfer();
	~landcover_change_transfer();
	void allocate();

	double ccont() {

		double ccont = transfer_acflux_harvest;

		for(int i=0; i<npft; i++) {
			ccont += transfer_litter_leaf[i];
			ccont += transfer_litter_sap[i];
			ccont += transfer_litter_heart[i];
			ccont += transfer_litter_root[i];
			ccont += transfer_litter_repr[i];
			ccont += transfer_harvested_products_slow[i];
		}
		for (int lyr = 0; lyr < NSOILLAYER_MAX; lyr++) {
			for (int i = 0; i < NSOMPOOL; i++)
				ccont += transfer_sompool[lyr][i].cmass;
		}

		return ccont;
	}

	double ncont() {

		double ncont = transfer_anflux_harvest;

		for(int i=0; i<npft; i++) {
			ncont += transfer_nmass_litter_leaf[i];
			ncont += transfer_nmass_litter_sap[i];
			ncont += transfer_nmass_litter_heart[i];
			ncont += transfer_nmass_litter_root[i];
			ncont += transfer_harvested_products_slow_nmass[i];
		}

		for (int lyr = 0; lyr < NSOILLAYER_MAX; lyr++) {
			for (int i = 0; i < NSOMPOOL; i++)
				ncont += transfer_sompool[lyr][i].nmass;

			ncont += transfer_N2_mass[lyr] + transfer_NH4_mass[lyr] + transfer_NO3_mass[lyr];
			ncont += transfer_N2O_mass[lyr] + transfer_NO_mass[lyr] + transfer_NO2_mass[lyr];
			ncont += transfer_N2O_mass_w[lyr] + transfer_NO_mass_w[lyr] + transfer_NO2_mass_w[lyr];
			ncont += transfer_N2O_mass_d[lyr] + transfer_NO_mass_d[lyr] + transfer_NO2_mass_d[lyr];
			ncont += transfer_NH4_mass_w[lyr] + transfer_NO3_mass_w[lyr];
			ncont += transfer_NH4_mass_d[lyr] + transfer_NO3_mass_d[lyr];
		}

		for (int slyr = 0; slyr < NLAYERS_SNOW_MAX; slyr++) {
			ncont += transfer_snowpack_NH4_mass[slyr] + transfer_snowpack_NO3_mass[slyr];
		}

		return ncont;
	}

	// Adds non-living C, N and water from a stand
	void add_from_stand(Stand& stand, double scale) {

		stand.firstobj();
		while(stand.isobj) {
			Patch& patch = stand.getobj();

			// sum original litter C & N:
			for(int n=0; n<npft; n++)
			{
				transfer_litter_leaf[n] += patch.pft[n].litter_leaf * scale;
				transfer_litter_root[n] += patch.pft[n].litter_root * scale;
				transfer_litter_sap[n] += patch.pft[n].litter_sap * scale;
				transfer_litter_heart[n] += patch.pft[n].litter_heart * scale;
				transfer_litter_repr[n] += patch.pft[n].litter_repr * scale;

				transfer_nmass_litter_leaf[n] += patch.pft[n].nmass_litter_leaf * scale;
				transfer_nmass_litter_root[n] += patch.pft[n].nmass_litter_root*scale;
				transfer_nmass_litter_sap[n] += patch.pft[n].nmass_litter_sap * scale;
				transfer_nmass_litter_heart[n] += patch.pft[n].nmass_litter_heart * scale;

				if(ifslowharvestpool) {
					transfer_harvested_products_slow[n] += patch.pft[n].harvested_products_slow * scale;
					transfer_harvested_products_slow_nmass[n] += patch.pft[n].harvested_products_slow_nmass * scale;
				}
			}

			for (int lyr = 0; lyr < patch.soilcolumn.nsoillayers; lyr++) {
				// sum litter C & N:
				transfer_cpool_fast[lyr] += patch.soilcolumn.soillayers[lyr].cpool_fast * scale;
				transfer_cpool_slow[lyr] += patch.soilcolumn.soillayers[lyr].cpool_slow * scale;

				// sum soil C & N:
				for (int i = 0; i < NSOMPOOL; i++) {
					transfer_sompool[lyr][i].cmass += patch.soilcolumn.soillayers[lyr].sompool[i].cmass * scale;
					transfer_sompool[lyr][i].fireresist += patch.soilcolumn.soillayers[lyr].sompool[i].fireresist * scale;
					transfer_sompool[lyr][i].fracremain += patch.soilcolumn.soillayers[lyr].sompool[i].fracremain * scale;
					transfer_sompool[lyr][i].ligcfrac += patch.soilcolumn.soillayers[lyr].sompool[i].ligcfrac * scale;
					transfer_sompool[lyr][i].litterme += patch.soilcolumn.soillayers[lyr].sompool[i].litterme * scale;
					transfer_sompool[lyr][i].nmass += patch.soilcolumn.soillayers[lyr].sompool[i].nmass * scale;
					transfer_sompool[lyr][i].ntoc += patch.soilcolumn.soillayers[lyr].sompool[i].ntoc * scale;
				}

				transfer_NH4_mass[lyr] += patch.soilcolumn.soillayers[lyr].NH4_mass * scale;
				transfer_NO3_mass[lyr] += patch.soilcolumn.soillayers[lyr].NO3_mass * scale;
				transfer_N2_mass[lyr] += patch.soilcolumn.soillayers[lyr].N2_mass * scale;
				transfer_NO_mass[lyr] += patch.soilcolumn.soillayers[lyr].NO_mass * scale;
				transfer_NO2_mass[lyr] += patch.soilcolumn.soillayers[lyr].NO2_mass * scale;
				transfer_N2O_mass[lyr] += patch.soilcolumn.soillayers[lyr].N2O_mass * scale;

				transfer_NH4_mass_w[lyr] += patch.soilcolumn.soillayers[lyr].NH4_mass_w * scale;
				transfer_NO3_mass_w[lyr] += patch.soilcolumn.soillayers[lyr].NO3_mass_w * scale;
				transfer_NO_mass_w[lyr] += patch.soilcolumn.soillayers[lyr].NO_mass_w * scale;
				transfer_NO2_mass_w[lyr] += patch.soilcolumn.soillayers[lyr].NO2_mass_w * scale;
				transfer_N2O_mass_w[lyr] += patch.soilcolumn.soillayers[lyr].N2O_mass_w * scale;

				transfer_NH4_mass_d[lyr] += patch.soilcolumn.soillayers[lyr].NH4_mass_d * scale;
				transfer_NO3_mass_d[lyr] += patch.soilcolumn.soillayers[lyr].NO3_mass_d * scale;
				transfer_NO_mass_d[lyr] += patch.soilcolumn.soillayers[lyr].NO_mass_d * scale;
				transfer_NO2_mass_d[lyr] += patch.soilcolumn.soillayers[lyr].NO2_mass_d * scale;
				transfer_N2O_mass_d[lyr] += patch.soilcolumn.soillayers[lyr].N2O_mass_d * scale;

				transfer_wcont[lyr] += patch.soilcolumn.soillayers[lyr].wcont * scale;

				transfer_decomp_litter_mean[lyr] += patch.soilcolumn.soillayers[lyr].decomp_litter_mean * scale;
				transfer_k_soilfast_mean[lyr] += patch.soilcolumn.soillayers[lyr].k_soilfast_mean * scale;
				transfer_k_soilslow_mean[lyr] += patch.soilcolumn.soillayers[lyr].k_soilslow_mean * scale;
			}

			for (int lyr = 0; lyr < patch.soilcolumn.snow.nsnowlayers; lyr++) {
				transfer_snowpack_NH4_mass[lyr] += patch.soilcolumn.snow.NH4_mass * scale; // DW_COMMENT SNOW
				transfer_snowpack_NO3_mass[lyr] += patch.soilcolumn.snow.NO3_mass * scale;
			}

			transfer_wcont_evap += patch.soilcolumn.wcont_evap * scale;
			double aaet_5_cp[NYEARAAET] = {0.0};
			transfer_aaet_5.to_array(aaet_5_cp);
			for(unsigned int i=0;i<patch.aaet_5.size();i++)
				transfer_aaet_5.add(aaet_5_cp[i] + patch.aaet_5[i] * scale);

			transfer_anfix_calc += patch.soilcolumn.anfix_calc * scale;

			stand.nextobj();
		}

		return;
	}

	void copy(landcover_change_transfer& from) {

		if(transfer_litter_leaf)
			memcpy(transfer_litter_leaf, from.transfer_litter_leaf, sizeof(double) * npft);
		if(transfer_litter_sap)
			memcpy(transfer_litter_sap, from.transfer_litter_sap, sizeof(double) * npft);
		if(transfer_litter_heart)
			memcpy(transfer_litter_heart, from.transfer_litter_heart, sizeof(double) * npft);
		if(transfer_litter_root)
			memcpy(transfer_litter_root, from.transfer_litter_root, sizeof(double) * npft);
		if(transfer_litter_repr)
			memcpy(transfer_litter_repr, from.transfer_litter_repr, sizeof(double) * npft);
		if(transfer_harvested_products_slow)
			memcpy(transfer_harvested_products_slow, from.transfer_harvested_products_slow, sizeof(double) * npft);
		if(transfer_nmass_litter_leaf)
			memcpy(transfer_nmass_litter_leaf, from.transfer_nmass_litter_leaf, sizeof(double) * npft);
		if(transfer_nmass_litter_sap)
			memcpy(transfer_nmass_litter_sap, from.transfer_nmass_litter_sap, sizeof(double) * npft);
		if(transfer_nmass_litter_heart)
			memcpy(transfer_nmass_litter_heart, from.transfer_nmass_litter_heart, sizeof(double) * npft);
		if(transfer_nmass_litter_root)
			memcpy(transfer_nmass_litter_root, from.transfer_nmass_litter_root, sizeof(double) * npft);
		if(transfer_harvested_products_slow_nmass)
			memcpy(transfer_harvested_products_slow_nmass, from.transfer_harvested_products_slow_nmass, sizeof(double) * npft);

		transfer_acflux_harvest = from.transfer_acflux_harvest;
		transfer_anflux_harvest = from.transfer_anflux_harvest;
		transfer_wcont_evap = from.transfer_wcont_evap;

		for (int lyr = 0; lyr < NSOILLAYER_MAX; lyr++) {
			transfer_cpool_fast[lyr] = from.transfer_cpool_fast[lyr];
			transfer_cpool_slow[lyr] = from.transfer_cpool_slow[lyr];
			transfer_wcont[lyr] = from.transfer_wcont[lyr];
			transfer_decomp_litter_mean[lyr] = from.transfer_decomp_litter_mean[lyr];
			transfer_k_soilfast_mean[lyr] = from.transfer_k_soilfast_mean[lyr];
			for (int i = 0; i < NSOMPOOL; i++) {
				transfer_sompool[lyr][i].cmass = from.transfer_sompool[lyr][i].cmass;
				transfer_sompool[lyr][i].fireresist = from.transfer_sompool[lyr][i].fireresist;
				transfer_sompool[lyr][i].fracremain = from.transfer_sompool[lyr][i].fracremain;
				transfer_sompool[lyr][i].ligcfrac = from.transfer_sompool[lyr][i].ligcfrac;
				transfer_sompool[lyr][i].litterme = from.transfer_sompool[lyr][i].litterme;
				transfer_sompool[lyr][i].nmass = from.transfer_sompool[lyr][i].nmass;
				transfer_sompool[lyr][i].ntoc = from.transfer_sompool[lyr][i].ntoc;
			}
			transfer_NH4_mass[lyr] = from.transfer_NH4_mass[lyr];
			transfer_NO3_mass[lyr] = from.transfer_NO3_mass[lyr];
			transfer_NO_mass[lyr] = from.transfer_NO_mass[lyr];
			transfer_NO2_mass[lyr] = from.transfer_NO2_mass[lyr];
			transfer_N2O_mass[lyr] = from.transfer_N2O_mass[lyr];
			transfer_N2_mass[lyr] = from.transfer_N2_mass[lyr];

			transfer_NH4_mass_w[lyr] = from.transfer_NH4_mass_w[lyr];
			transfer_NO3_mass_w[lyr] = from.transfer_NO3_mass_w[lyr];
			transfer_NO_mass_w[lyr] = from.transfer_NO_mass_w[lyr];
			transfer_NO2_mass_w[lyr] = from.transfer_NO2_mass_w[lyr];
			transfer_N2O_mass_w[lyr] = from.transfer_N2O_mass_w[lyr];

			transfer_NH4_mass_d[lyr] = from.transfer_NH4_mass_d[lyr];
			transfer_NO3_mass_d[lyr] = from.transfer_NO3_mass_d[lyr];
			transfer_NO_mass_d[lyr] = from.transfer_NO_mass_d[lyr];
			transfer_NO2_mass_d[lyr] = from.transfer_NO2_mass_d[lyr];
			transfer_N2O_mass_d[lyr] = from.transfer_N2O_mass_d[lyr];
		}

		for (int slyr = 0; slyr < NLAYERS_SNOW_MAX; slyr++) {
			transfer_snowpack_NH4_mass[slyr] = from.transfer_snowpack_NH4_mass[slyr];
			transfer_snowpack_NO3_mass[slyr] = from.transfer_snowpack_NO3_mass[slyr];
		}

		for(unsigned int i=0;i<from.transfer_aaet_5.size();i++)
			transfer_aaet_5.add(from.transfer_aaet_5[i]);

		transfer_anfix_calc = from.transfer_anfix_calc;
	}

	void add(landcover_change_transfer& from, double multiplier = 1.0) {

		for(int i=0; i<npft; i++) {
			if(transfer_litter_leaf)
				transfer_litter_leaf[i] += from.transfer_litter_leaf[i] * multiplier;
			if(transfer_litter_sap)
				transfer_litter_sap[i] += from.transfer_litter_sap[i] * multiplier;
			if(transfer_litter_heart)
				transfer_litter_heart[i] += from.transfer_litter_heart[i] * multiplier;
			if(transfer_litter_root)
				transfer_litter_root[i] += from.transfer_litter_root[i] * multiplier;
			if(transfer_litter_repr)
				transfer_litter_repr[i] += from.transfer_litter_repr[i] * multiplier;
			if(transfer_harvested_products_slow)
				transfer_harvested_products_slow[i] += from.transfer_harvested_products_slow[i] * multiplier;
			if(transfer_nmass_litter_leaf)
				transfer_nmass_litter_leaf[i] += from.transfer_nmass_litter_leaf[i] * multiplier;
			if(transfer_nmass_litter_sap)
				transfer_nmass_litter_sap[i] += from.transfer_nmass_litter_sap[i] * multiplier;
			if(transfer_nmass_litter_heart)
				transfer_nmass_litter_heart[i] += from.transfer_nmass_litter_heart[i] * multiplier;
			if(transfer_nmass_litter_root)
				transfer_nmass_litter_root[i] += from.transfer_nmass_litter_root[i] * multiplier;
			if(transfer_harvested_products_slow_nmass)
				transfer_harvested_products_slow_nmass[i] += from.transfer_harvested_products_slow_nmass[i] * multiplier;
		}
		transfer_acflux_harvest += from.transfer_acflux_harvest * multiplier;
		transfer_anflux_harvest += from.transfer_anflux_harvest * multiplier;
		transfer_wcont_evap += from.transfer_wcont_evap * multiplier;

		for (int lyr = 0; lyr < NSOILLAYER_MAX; lyr++) {
			transfer_cpool_fast[lyr] += from.transfer_cpool_fast[lyr] * multiplier;
			transfer_cpool_slow[lyr] += from.transfer_cpool_slow[lyr] * multiplier;
			transfer_wcont[lyr] += from.transfer_wcont[lyr] * multiplier;
			transfer_decomp_litter_mean[lyr] += from.transfer_decomp_litter_mean[lyr] * multiplier;
			transfer_k_soilfast_mean[lyr] += from.transfer_k_soilfast_mean[lyr] * multiplier;
			for (int i = 0; i < NSOMPOOL; i++) {
				transfer_sompool[lyr][i].cmass += from.transfer_sompool[lyr][i].cmass * multiplier;
				transfer_sompool[lyr][i].fireresist += from.transfer_sompool[lyr][i].fireresist * multiplier;
				transfer_sompool[lyr][i].fracremain += from.transfer_sompool[lyr][i].fracremain * multiplier;
				transfer_sompool[lyr][i].ligcfrac += from.transfer_sompool[lyr][i].ligcfrac * multiplier;
				transfer_sompool[lyr][i].litterme += from.transfer_sompool[lyr][i].litterme * multiplier;
				transfer_sompool[lyr][i].nmass += from.transfer_sompool[lyr][i].nmass * multiplier;
				transfer_sompool[lyr][i].ntoc += from.transfer_sompool[lyr][i].ntoc * multiplier;
			}
			transfer_NH4_mass[lyr] += from.transfer_NH4_mass[lyr] * multiplier;
			transfer_NO3_mass[lyr] += from.transfer_NO3_mass[lyr] * multiplier;
			transfer_N2O_mass[lyr] += from.transfer_N2O_mass[lyr] * multiplier;
			transfer_N2_mass[lyr] += from.transfer_N2_mass[lyr] * multiplier;
			transfer_NO_mass[lyr] += from.transfer_NO_mass[lyr] * multiplier;
			transfer_NO2_mass[lyr] += from.transfer_NO2_mass[lyr] * multiplier;

			transfer_NH4_mass_w[lyr] += from.transfer_NH4_mass_w[lyr] * multiplier;
			transfer_NO3_mass_w[lyr] += from.transfer_NO3_mass_w[lyr] * multiplier;
			transfer_N2O_mass_w[lyr] += from.transfer_N2O_mass_w[lyr] * multiplier;
			transfer_NO_mass_w[lyr] += from.transfer_NO_mass_w[lyr] * multiplier;
			transfer_NO2_mass_w[lyr] += from.transfer_NO2_mass_w[lyr] * multiplier;

			transfer_NH4_mass_d[lyr] += from.transfer_NH4_mass_d[lyr] * multiplier;
			transfer_NO3_mass_d[lyr] += from.transfer_NO3_mass_d[lyr] * multiplier;
			transfer_N2O_mass_d[lyr] += from.transfer_N2O_mass_d[lyr] * multiplier;
			transfer_NO_mass_d[lyr] += from.transfer_NO_mass_d[lyr] * multiplier;
			transfer_NO2_mass_d[lyr] += from.transfer_NO2_mass_d[lyr] * multiplier;
		}

		for (int slyr = 0; slyr < NLAYERS_SNOW_MAX; slyr++) {
			transfer_snowpack_NH4_mass[slyr] += from.transfer_snowpack_NH4_mass[slyr] * multiplier;
			transfer_snowpack_NO3_mass[slyr] += from.transfer_snowpack_NO3_mass[slyr] * multiplier;
		}

		double aaet_5_cp[NYEARAAET] = {0.0};
		for(unsigned int i=0;i<transfer_aaet_5.size();i++)
			aaet_5_cp[i] = transfer_aaet_5[i];
		for(unsigned int i=0;i<from.transfer_aaet_5.size();i++)
			transfer_aaet_5.add(aaet_5_cp[i] + from.transfer_aaet_5[i] * multiplier);

		transfer_anfix_calc += from.transfer_anfix_calc * multiplier;
	}
};

struct lc_change_harvest_params {

	double harv_eff;
	double res_outtake_twig;
	double res_outtake_coarse_root;
	bool burn;	// to be developed

	lc_change_harvest_params() {

		res_outtake_twig = 0.0;
		res_outtake_coarse_root = 0.0;
		burn = false;
	}
};

#endif // LPJ_GUESS_LANDCOVER_H
