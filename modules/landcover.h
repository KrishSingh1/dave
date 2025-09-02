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

	// Patchpft litter
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
	double *transfer_nmass_litter_NH4;
	double *transfer_harvested_products_slow_nmass;

	double transfer_acflux_harvest;
	double transfer_anflux_harvest;

	// Soil
	double transfer_depth;
	double transfer_wcont_evap;
	Historic<double, NYEARAAET> transfer_aaet_5;
	Historic<double, 10> transfer_maxthaw_10;
	double transfer_anfix_calc;
	double transfer_wmass_ssw;
	double transfer_imass_ssw;
	double transfer_temp_ssw;
	int transfer_soil_nsoillayers;

	// Soil layers
	double transfer_soil_soillayerdepths[NSOILLAYERS_MAX];
	double transfer_soil_soillayerthicks[NSOILLAYERS_MAX];
	double transfer_soil_depth[NSOILLAYERS_MAX];
	double transfer_soil_depth_profile[NSOILLAYERS_MAX];
	double transfer_soil_wmass[NSOILLAYERS_MAX];
	double transfer_soil_imass[NSOILLAYERS_MAX];
	double transfer_soil_mineral_volume[NSOILLAYERS_MAX];
	double transfer_soil_organic_volume[NSOILLAYERS_MAX];
	double transfer_soil_root_volume[NSOILLAYERS_MAX];
	double transfer_soil_temp[NSOILLAYERS_MAX];
	double transfer_soil_CO2[NSOILLAYERS_MAX];
	double transfer_soil_CH4[NSOILLAYERS_MAX];
	double transfer_soil_O2[NSOILLAYERS_MAX];
	double transfer_soil_NH4[NSOILLAYERS_MAX];
	double transfer_soil_NO3[NSOILLAYERS_MAX];
	double transfer_soil_NO2[NSOILLAYERS_MAX];
	double transfer_soil_NO[NSOILLAYERS_MAX];
	double transfer_soil_N2O[NSOILLAYERS_MAX];
	double transfer_soil_N2[NSOILLAYERS_MAX];
	double transfer_soil_depth_mod[NSOILLAYERS_MAX];
	double transfer_soil_cpool_slow[NSOILLAYERS_MAX];
	double transfer_soil_cpool_fast[NSOILLAYERS_MAX];
	double transfer_soil_sompool_cmass[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_nmass[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_ntoc[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_fireresist[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_fracremain[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_ligcfrac[NSOILLAYERS_MAX][NSOMPOOL];
	double transfer_soil_sompool_litterme[NSOILLAYERS_MAX][NSOMPOOL];

	// Soil padding
	double transfer_soilpad_depth[NPADLAYERS];
	double transfer_soilpad_temp[NPADLAYERS];

	// Snow
	int transfer_snow_nsnowlayers;
	double transfer_snow_NH4_mass;
	double transfer_snow_NO3_mass;
	double transfer_snow_wmass[NSNOWLAYERS_MAX];
	double transfer_snow_imass[NSNOWLAYERS_MAX];
	double transfer_snow_depth[NSNOWLAYERS_MAX];
	double transfer_snow_temp[NSNOWLAYERS_MAX];

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
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			for (int i = 0; i < NSOMPOOL; i++) {
				ccont += transfer_soil_sompool_cmass[lyr][i];
			}
			ccont += transfer_soil_CH4[lyr] + transfer_soil_CO2[lyr];
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
			ncont += transfer_nmass_litter_NH4[i];
			ncont += transfer_harvested_products_slow_nmass[i];
		}

		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			for (int i = 0; i < NSOMPOOL; i++)
				ncont += transfer_soil_sompool_nmass[lyr][i];

			ncont += transfer_soil_N2[lyr] + transfer_soil_NH4[lyr] + transfer_soil_NO3[lyr];
			ncont += transfer_soil_N2O[lyr] + transfer_soil_NO[lyr] + transfer_soil_NO2[lyr];
		}

		ncont += transfer_snow_NH4_mass + transfer_snow_NO3_mass;

		return ncont;
	}

	double wcont() {

		double wcont = 0.0;
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			wcont += transfer_soil_wmass[lyr] + transfer_soil_imass[lyr];
		}
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			wcont += transfer_snow_wmass[lyr] + transfer_snow_imass[lyr];
		}

		return wcont;
	}

	// Adds non-living C, N and water from a stand
	void add_from_stand(Stand& stand, double scale) {

		stand.firstobj();
		while(stand.isobj) {
			Patch& patch = stand.getobj();

			// Sum patchpft litter
			for(int n=0; n<npft; n++) {

				Patchpft& ppft = patch.pft[n];

				transfer_litter_leaf[n] += ppft.litter_leaf * scale;
				transfer_litter_root[n] += ppft.litter_root * scale;
				transfer_litter_sap[n] += ppft.litter_sap * scale;
				transfer_litter_heart[n] += ppft.litter_heart * scale;
				transfer_litter_repr[n] += ppft.litter_repr * scale;

				transfer_nmass_litter_leaf[n] += ppft.nmass_litter_leaf * scale;
				transfer_nmass_litter_root[n] += ppft.nmass_litter_root * scale;
				transfer_nmass_litter_sap[n] += ppft.nmass_litter_sap * scale;
				transfer_nmass_litter_heart[n] += ppft.nmass_litter_heart * scale;
				transfer_nmass_litter_NH4[n] += ppft.nmass_litter_NH4 * scale;

				if(ifslowharvestpool) {
					transfer_harvested_products_slow[n] += ppft.harvested_products_slow * scale;
					transfer_harvested_products_slow_nmass[n] += ppft.harvested_products_slow_nmass * scale;
				}
			}

			// Sum soil
			Soilcolumn& soilcolumn = patch.soilcolumn;
			transfer_depth += soilcolumn.depth * scale;
			transfer_wcont_evap += soilcolumn.wcont_evap * scale;
			double aaet_5_cp[NYEARAAET] = { 0.0 };
			transfer_aaet_5.to_array(aaet_5_cp);
			for (unsigned int i = 0; i < patch.aaet_5.size(); i++) {
				transfer_aaet_5.add(aaet_5_cp[i] + patch.aaet_5[i] * scale);
			}
			double maxthaw_10_cp[10] = { 0.0 };
			transfer_maxthaw_10.to_array(maxthaw_10_cp);
			for (unsigned int i = 0; i < soilcolumn.maxthaw_10.size(); i++) {
				transfer_maxthaw_10.add(maxthaw_10_cp[i] + soilcolumn.maxthaw_10[i] * scale);
			}
			transfer_anfix_calc += soilcolumn.anfix_calc * scale;
			transfer_wmass_ssw += soilcolumn.ssw.wmass * scale;
			transfer_imass_ssw += soilcolumn.ssw.imass * scale;
			transfer_temp_ssw += soilcolumn.ssw.temp * scale;
			transfer_soil_nsoillayers = max(transfer_soil_nsoillayers, soilcolumn.nsoillayers);

			// Sum soil layers
			for (int lyr = 0; lyr < soilcolumn.nsoillayers; lyr++) {
				Soillayer& soillayer = soilcolumn.soillayers[lyr];
				transfer_soil_soillayerdepths[lyr] += soilcolumn.soillayerdepths[lyr] * scale;
				transfer_soil_soillayerthicks[lyr] += soilcolumn.soillayerthicks[lyr] * scale;
				transfer_soil_depth[lyr] += soillayer.depth * scale;
				transfer_soil_depth_profile[lyr] += soillayer.depth_profile * scale;
				transfer_soil_wmass[lyr] += soillayer.wmass * scale;
				transfer_soil_imass[lyr] += soillayer.imass * scale;
				transfer_soil_mineral_volume[lyr] += soillayer.mineral_volume * scale;
				transfer_soil_organic_volume[lyr] += soillayer.organic_volume * scale;
				transfer_soil_root_volume[lyr] += soillayer.root_volume * scale;
				transfer_soil_temp[lyr] += soillayer.temp * scale;
				transfer_soil_CO2[lyr] += soillayer.CO2 * scale;
				transfer_soil_CH4[lyr] += soillayer.CH4 * scale;
				transfer_soil_O2[lyr] += soillayer.O2 * scale;
				transfer_soil_NH4[lyr] += soillayer.NH4_mass * scale;
				transfer_soil_NO3[lyr] += soillayer.NO3_mass * scale;
				transfer_soil_NO2[lyr] += soillayer.NO2_mass * scale;
				transfer_soil_NO[lyr] += soillayer.NO_mass * scale;
				transfer_soil_N2O[lyr] += soillayer.N2O_mass * scale;
				transfer_soil_N2[lyr] += soillayer.N2_mass * scale;
				transfer_soil_depth_mod[lyr] += soillayer.depth_mod * scale;
				transfer_soil_cpool_slow[lyr] += soillayer.cpool_slow * scale;
				transfer_soil_cpool_fast[lyr] += soillayer.cpool_fast * scale;
				for (int p = 0; p < NSOMPOOL; p++) {
					Sompool& sompool = soillayer.sompool[p];
					transfer_soil_sompool_cmass[lyr][p] += sompool.cmass * scale;
					transfer_soil_sompool_nmass[lyr][p] += sompool.nmass * scale;
					transfer_soil_sompool_ntoc[lyr][p] += sompool.ntoc * scale;
					transfer_soil_sompool_fireresist[lyr][p] += sompool.fireresist * scale;
					transfer_soil_sompool_fracremain[lyr][p] += sompool.fracremain * scale;
					transfer_soil_sompool_ligcfrac[lyr][p] += sompool.ligcfrac * scale;
					transfer_soil_sompool_litterme[lyr][p] += sompool.litterme * scale;
				}
			}

			// Soil padding
			for (int lyr = 0; lyr < NPADLAYERS; lyr++) {
				padding& pad = soilcolumn.pad[lyr];
				transfer_soilpad_depth[lyr] += pad.depth * scale;
				transfer_soilpad_temp[lyr] += pad.temp * scale;
			}

			// Sum snow
			Snow& snow = soilcolumn.snow;
			transfer_snow_nsnowlayers = max(transfer_snow_nsnowlayers, snow.nsnowlayers);
			transfer_snow_NH4_mass += snow.NH4_mass * scale;
			transfer_snow_NO3_mass += snow.NO3_mass * scale;
			for (int lyr = 0; lyr < snow.nsnowlayers; lyr++) {
				Snowlayer& snowlayer = snow.snowlayers[lyr];
				transfer_snow_wmass[lyr] += snowlayer.wmass * scale;
				transfer_snow_imass[lyr] += snowlayer.imass * scale;
				transfer_snow_depth[lyr] += snowlayer.depth * scale;
				transfer_snow_temp[lyr] += snowlayer.temp * scale;
			}

			stand.nextobj();
		}

		return;
	}

	void copy(landcover_change_transfer& from) {

		if (transfer_litter_leaf)
			memcpy(transfer_litter_leaf, from.transfer_litter_leaf, sizeof(double) * npft);
		if (transfer_litter_sap)
			memcpy(transfer_litter_sap, from.transfer_litter_sap, sizeof(double) * npft);
		if (transfer_litter_heart)
			memcpy(transfer_litter_heart, from.transfer_litter_heart, sizeof(double) * npft);
		if (transfer_litter_root)
			memcpy(transfer_litter_root, from.transfer_litter_root, sizeof(double) * npft);
		if (transfer_litter_repr)
			memcpy(transfer_litter_repr, from.transfer_litter_repr, sizeof(double) * npft);
		if (transfer_harvested_products_slow)
			memcpy(transfer_harvested_products_slow, from.transfer_harvested_products_slow, sizeof(double) * npft);
		if (transfer_nmass_litter_leaf)
			memcpy(transfer_nmass_litter_leaf, from.transfer_nmass_litter_leaf, sizeof(double) * npft);
		if (transfer_nmass_litter_sap)
			memcpy(transfer_nmass_litter_sap, from.transfer_nmass_litter_sap, sizeof(double) * npft);
		if (transfer_nmass_litter_heart)
			memcpy(transfer_nmass_litter_heart, from.transfer_nmass_litter_heart, sizeof(double) * npft);
		if (transfer_nmass_litter_root)
			memcpy(transfer_nmass_litter_root, from.transfer_nmass_litter_root, sizeof(double) * npft);
		if (transfer_nmass_litter_NH4)
			memcpy(transfer_nmass_litter_NH4, from.transfer_nmass_litter_NH4, sizeof(double) * npft);
		if (transfer_harvested_products_slow_nmass)
			memcpy(transfer_harvested_products_slow_nmass, from.transfer_harvested_products_slow_nmass, sizeof(double) * npft);

		transfer_acflux_harvest = from.transfer_acflux_harvest;
		transfer_anflux_harvest = from.transfer_anflux_harvest;

		// Soil
		transfer_depth = from.transfer_depth;
		transfer_wcont_evap = from.transfer_wcont_evap;
		for (unsigned int i = 0; i < from.transfer_aaet_5.size(); i++) {
			transfer_aaet_5.add(from.transfer_aaet_5[i]);
		}
		for (unsigned int i = 0; i < from.transfer_maxthaw_10.size(); i++) {
			transfer_maxthaw_10.add(from.transfer_maxthaw_10[i]);
		}
		transfer_anfix_calc = from.transfer_anfix_calc;
		transfer_wmass_ssw = from.transfer_wmass_ssw;
		transfer_imass_ssw = from.transfer_imass_ssw;
		transfer_temp_ssw = from.transfer_temp_ssw;
		transfer_soil_nsoillayers = from.transfer_soil_nsoillayers;

		// Sum soil layers
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			transfer_soil_soillayerdepths[lyr] = from.transfer_soil_soillayerdepths[lyr];
			transfer_soil_soillayerthicks[lyr] = from.transfer_soil_soillayerthicks[lyr];
			transfer_soil_depth[lyr] = from.transfer_soil_depth[lyr];
			transfer_soil_depth_profile[lyr] = from.transfer_soil_depth_profile[lyr];
			transfer_soil_wmass[lyr] = from.transfer_soil_wmass[lyr];
			transfer_soil_imass[lyr] = from.transfer_soil_imass[lyr];
			transfer_soil_mineral_volume[lyr] = from.transfer_soil_mineral_volume[lyr];
			transfer_soil_organic_volume[lyr] = from.transfer_soil_organic_volume[lyr];
			transfer_soil_root_volume[lyr] = from.transfer_soil_root_volume[lyr];
			transfer_soil_temp[lyr] = from.transfer_soil_temp[lyr];
			transfer_soil_CO2[lyr] = from.transfer_soil_CO2[lyr];
			transfer_soil_CH4[lyr] = from.transfer_soil_CH4[lyr];
			transfer_soil_O2[lyr] = from.transfer_soil_O2[lyr];
			transfer_soil_NH4[lyr] = from.transfer_soil_NH4[lyr];
			transfer_soil_NO3[lyr] = from.transfer_soil_NO3[lyr];
			transfer_soil_NO2[lyr] = from.transfer_soil_NO2[lyr];
			transfer_soil_NO[lyr] = from.transfer_soil_NO[lyr];
			transfer_soil_N2O[lyr] = from.transfer_soil_N2O[lyr];
			transfer_soil_N2[lyr] = from.transfer_soil_N2[lyr];
			transfer_soil_depth_mod[lyr] = from.transfer_soil_depth_mod[lyr];
			transfer_soil_cpool_slow[lyr] = from.transfer_soil_cpool_slow[lyr];
			transfer_soil_cpool_fast[lyr] = from.transfer_soil_cpool_fast[lyr];
			for (int p = 0; p < NSOMPOOL; p++) {
				transfer_soil_sompool_cmass[lyr][p] = from.transfer_soil_sompool_cmass[lyr][p];
				transfer_soil_sompool_nmass[lyr][p] = from.transfer_soil_sompool_nmass[lyr][p];
				transfer_soil_sompool_ntoc[lyr][p] = from.transfer_soil_sompool_ntoc[lyr][p];
				transfer_soil_sompool_fireresist[lyr][p] = from.transfer_soil_sompool_fireresist[lyr][p];
				transfer_soil_sompool_fracremain[lyr][p] = from.transfer_soil_sompool_fracremain[lyr][p];
				transfer_soil_sompool_ligcfrac[lyr][p] = from.transfer_soil_sompool_ligcfrac[lyr][p];
				transfer_soil_sompool_litterme[lyr][p] = from.transfer_soil_sompool_litterme[lyr][p];
			}
		}

		// Soil padding
		for (int lyr = 0; lyr < NPADLAYERS; lyr++) {
			transfer_soilpad_depth[lyr] = from.transfer_soilpad_depth[lyr];
			transfer_soilpad_temp[lyr] = from.transfer_soilpad_temp[lyr];
		}

		// Sum snow
		transfer_snow_nsnowlayers = from.transfer_snow_nsnowlayers;
		transfer_snow_NH4_mass = from.transfer_snow_NH4_mass;
		transfer_snow_NO3_mass = from.transfer_snow_NO3_mass;
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			transfer_snow_wmass[lyr] = from.transfer_snow_wmass[lyr];
			transfer_snow_imass[lyr] = from.transfer_snow_imass[lyr];
			transfer_snow_depth[lyr] = from.transfer_snow_depth[lyr];
			transfer_snow_temp[lyr] = from.transfer_snow_temp[lyr];
		}
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
			if (transfer_nmass_litter_NH4)
				transfer_nmass_litter_NH4[i] += from.transfer_nmass_litter_NH4[i] * multiplier;
			if(transfer_harvested_products_slow_nmass)
				transfer_harvested_products_slow_nmass[i] += from.transfer_harvested_products_slow_nmass[i] * multiplier;
		}
		transfer_acflux_harvest += from.transfer_acflux_harvest * multiplier;
		transfer_anflux_harvest += from.transfer_anflux_harvest * multiplier;

		// Soil
		transfer_depth += from.transfer_depth * multiplier;
		transfer_wcont_evap += from.transfer_wcont_evap * multiplier;
		double aaet_5_cp[NYEARAAET] = { 0.0 };
		for (unsigned int i = 0; i < transfer_aaet_5.size(); i++)
			aaet_5_cp[i] = transfer_aaet_5[i];
		for (unsigned int i = 0; i < from.transfer_aaet_5.size(); i++)
			transfer_aaet_5.add(aaet_5_cp[i] + from.transfer_aaet_5[i] * multiplier);
		double maxthaw_10_cp[10] = { 0.0 };
		for (unsigned int i = 0; i < transfer_maxthaw_10.size(); i++)
			maxthaw_10_cp[i] = transfer_maxthaw_10[i];
		for (unsigned int i = 0; i < from.transfer_maxthaw_10.size(); i++)
			transfer_maxthaw_10.add(maxthaw_10_cp[i] + from.transfer_maxthaw_10[i] * multiplier);
		transfer_anfix_calc += from.transfer_anfix_calc * multiplier;
		transfer_wmass_ssw += from.transfer_wmass_ssw * multiplier;
		transfer_imass_ssw += from.transfer_imass_ssw * multiplier;
		transfer_temp_ssw += from.transfer_temp_ssw * multiplier;
		transfer_soil_nsoillayers = max(transfer_soil_nsoillayers, from.transfer_soil_nsoillayers);

		// Soil layers
		for (int lyr = 0; lyr < NSOILLAYERS_MAX; lyr++) {
			transfer_soil_soillayerdepths[lyr] += from.transfer_soil_soillayerdepths[lyr] * multiplier;
			transfer_soil_soillayerthicks[lyr] += from.transfer_soil_soillayerthicks[lyr] * multiplier;
			transfer_soil_depth[lyr] += from.transfer_soil_depth[lyr] * multiplier;
			transfer_soil_depth_profile[lyr] += from.transfer_soil_depth_profile[lyr] * multiplier;
			transfer_soil_wmass[lyr] += from.transfer_soil_wmass[lyr] * multiplier;
			transfer_soil_imass[lyr] += from.transfer_soil_imass[lyr] * multiplier;
			transfer_soil_mineral_volume[lyr] += from.transfer_soil_mineral_volume[lyr] * multiplier;
			transfer_soil_organic_volume[lyr] += from.transfer_soil_organic_volume[lyr] * multiplier;
			transfer_soil_root_volume[lyr] += from.transfer_soil_root_volume[lyr] * multiplier;
			transfer_soil_temp[lyr] += from.transfer_soil_temp[lyr] * multiplier;
			transfer_soil_CO2[lyr] += from.transfer_soil_CO2[lyr] * multiplier;
			transfer_soil_CH4[lyr] += from.transfer_soil_CH4[lyr] * multiplier;
			transfer_soil_O2[lyr] += from.transfer_soil_O2[lyr] * multiplier;
			transfer_soil_NH4[lyr] += from.transfer_soil_NH4[lyr] * multiplier;
			transfer_soil_NO3[lyr] += from.transfer_soil_NO3[lyr] * multiplier;
			transfer_soil_NO2[lyr] += from.transfer_soil_NO2[lyr] * multiplier;
			transfer_soil_NO[lyr] += from.transfer_soil_NO[lyr] * multiplier;
			transfer_soil_N2O[lyr] += from.transfer_soil_N2O[lyr] * multiplier;
			transfer_soil_N2[lyr] += from.transfer_soil_N2[lyr] * multiplier;
			transfer_soil_depth_mod[lyr] += from.transfer_soil_depth_mod[lyr] * multiplier;
			transfer_soil_cpool_slow[lyr] += from.transfer_soil_cpool_slow[lyr] * multiplier;
			transfer_soil_cpool_fast[lyr] += from.transfer_soil_cpool_fast[lyr] * multiplier;
			for (int p = 0; p < NSOMPOOL; p++) {
				transfer_soil_sompool_cmass[lyr][p] += from.transfer_soil_sompool_cmass[lyr][p] * multiplier;
				transfer_soil_sompool_nmass[lyr][p] += from.transfer_soil_sompool_nmass[lyr][p] * multiplier;
				transfer_soil_sompool_ntoc[lyr][p] += from.transfer_soil_sompool_ntoc[lyr][p] * multiplier;
				transfer_soil_sompool_fireresist[lyr][p] += from.transfer_soil_sompool_fireresist[lyr][p] * multiplier;
				transfer_soil_sompool_fracremain[lyr][p] += from.transfer_soil_sompool_fracremain[lyr][p] * multiplier;
				transfer_soil_sompool_ligcfrac[lyr][p] += from.transfer_soil_sompool_ligcfrac[lyr][p] * multiplier;
				transfer_soil_sompool_litterme[lyr][p] += from.transfer_soil_sompool_litterme[lyr][p] * multiplier;
			}
		}

		// Soil padding
		for (int lyr = 0; lyr < NPADLAYERS; lyr++) {
			transfer_soilpad_depth[lyr] = from.transfer_soilpad_depth[lyr] * multiplier;
			transfer_soilpad_temp[lyr] = from.transfer_soilpad_temp[lyr] * multiplier;
		}

		// Snow
		transfer_snow_nsnowlayers = max(transfer_snow_nsnowlayers, from.transfer_snow_nsnowlayers);
		transfer_snow_NH4_mass += from.transfer_snow_NH4_mass * multiplier;
		transfer_snow_NO3_mass += from.transfer_snow_NO3_mass * multiplier;
		for (int lyr = 0; lyr < NSNOWLAYERS_MAX; lyr++) {
			transfer_snow_wmass[lyr] += from.transfer_snow_wmass[lyr] * multiplier;
			transfer_snow_imass[lyr] += from.transfer_snow_imass[lyr] * multiplier;
			transfer_snow_depth[lyr] += from.transfer_snow_depth[lyr] * multiplier;
			transfer_snow_temp[lyr] += from.transfer_snow_temp[lyr] * multiplier;
		}		
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
