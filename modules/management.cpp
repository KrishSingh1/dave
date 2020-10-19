////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file management.cpp
/// \brief Harvest functions for cropland, managed forest and pasture
/// \author Mats Lindeskog
/// $Date:  $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "landcover.h"
#include "management.h"
#include "driver.h"

/// Functions to check available wood for harvest at individual, patch and stand levels
double check_harvest_cmass(Individual& indiv, bool wood_cmass_only) {

	double cmass_harvest;
	Stand& stand = indiv.vegetation.patch.stand;
	Harvest_CN cp;

	cp.copy_from_indiv(indiv, indiv.has_daily_turnover(), false);

	// Forestry:
	double harv_eff_wood_harvest = indiv.pft.harv_eff;				// 0.9
	double res_outtake_twig_wood_harvest = indiv.pft.res_outtake;	// 0.4
	double res_outtake_coarse_root_wood_harvest = 0.1;

	// Harvest of transferred areas:
	harvest_wood(cp, indiv.pft, indiv.alive, 1.0, harv_eff_wood_harvest, res_outtake_twig_wood_harvest, res_outtake_coarse_root_wood_harvest);
	if(wood_cmass_only)
		cmass_harvest = (cp.acflux_harvest_wood) * stand.get_gridcell_fraction() / (double)stand.nobj;
	else
		cmass_harvest = (cp.acflux_harvest + cp.harvested_products_slow) * stand.get_gridcell_fraction() / (double)stand.nobj;

	return cmass_harvest;
}

double check_harvest_cmass(Patch& patch, bool wood_cmass_only, bool check_selection) {

	double cmass_harvest = 0.0;
	ManagementType& mt = patch.stand.get_current_management();
	StandType& st = stlist[patch.stand.stid];
	Vegetation& vegetation = patch.vegetation;

	for(unsigned int i=0;i<vegetation.nobj;i++) {
		Individual& indiv = vegetation[i];

		if(!check_selection || mt.pftinselection((const char*)indiv.pft.name) || mt.planting_system == "MONOCULTURE" && mt.pftname == indiv.pft.name || !st.restrictpfts)
			cmass_harvest += check_harvest_cmass(indiv, wood_cmass_only);
	}
	return cmass_harvest;
}

double check_harvest_cmass(Stand& stand, bool wood_cmass_only, bool check_selection) {

	double cmass_harvest = 0.0;

	for(unsigned int i=0;i<stand.nobj;i++) {
		Patch& patch = stand[i];
		cmass_harvest += check_harvest_cmass(patch, wood_cmass_only, check_selection);
	}
	return cmass_harvest;
}

/// Harvest function used for managed forest and for clearing natural vegetation at land use change
/** A fraction of trees is cut down (frac_cut)
 *  A fraction of wood is harvested (pft.harv_eff) and returned as acflux_harvest
 *  A fraction of harvested wood (pft.harvest_slow_frac) is returned as harvested_products_slow
 *  The rest, including leaves and roots, is returned as litter, unless a fraction of twigs or roots removed.
 *  Called from landcover_dynamics() first day of the year if any natural vegetation is transferred to another land use.
 *  INPUT PARAMETER
 *  \param frac_cut					fraction of trees cut
 *  \param harv_eff					harvest efficiency
 *  \param res_outtake_twig			removed twig fraction
 *  \param res_outtake_coarse_root	removed course root fraction
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgC/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - litter_sap   				new sapwood C litter (kgC/m2)
 *   - litter_heart   				new heartwood C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *   - nmass_litter_sap 			new sapwood nitrogen litter (kgN/m2)
 *   - nmass_litter_heart        	new heartwood nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_wood(Harvest_CN& i, Pft& pft, bool alive, double frac_cut, double harv_eff, double res_outtake_twig, double res_outtake_coarse_root) {

	double harvest = 0.0;
	double residue_outtake = 0.0;
	/// Fraction of wood cmass that are stems
	double stem_frac = 0.65;	// Temporary values, should be pft-specific
	/// Fraction of wood cmass that are twigs
	double twig_frac = 0.13;
	/// Fraction of wood cmass that are coarse roots
	double coarse_root_frac = 1.0 - stem_frac - twig_frac;	// 0.22 with default stem_frac and twig_frac values
	/// Fraction of leaves adhering to twigs at the time of removal
	double adhering_leaf_frac = 0.75;

	// only harvest trees
	if (pft.lifeform == GRASS)
		return;

	// all root carbon and nitrogen goes to litter
	if (alive) {

		i.litter_root += i.cmass_root * frac_cut;
		i.cmass_root *= (1.0 - frac_cut);
	}

	i.nmass_litter_root += i.nmass_root * frac_cut;
	i.nmass_litter_root += (i.nstore_labile + i.nstore_longterm) * frac_cut;
	i.nmass_root *= (1.0 - frac_cut);
	i.nstore_labile *= (1.0 - frac_cut);
	i.nstore_longterm *= (1.0 - frac_cut);

	if (alive) {

		// Carbon:

		if (i.cmass_debt <= i.cmass_sap + i.cmass_heart) {

			// harvested stem wood
			harvest += harv_eff * stem_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;
			i.acflux_harvest_wood += harvest;

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if (ifslowharvestpool) {
				i.harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1.0 - pft.harvest_slow_frac);
			}

			// harvested products consumed (oxidised) this year put into acflux_harvest
			i.acflux_harvest += harvest;

			// removed leaves adhering to twigs
			residue_outtake += res_outtake_twig * adhering_leaf_frac * i.cmass_leaf * frac_cut;

			// removed twigs
			residue_outtake += res_outtake_twig * twig_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

			// removed coarse roots
			residue_outtake += res_outtake_coarse_root * coarse_root_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

			// removed residues are oxidised
			i.acflux_harvest += residue_outtake;

			// not removed residues are put into litter
			i.litter_leaf += i.cmass_leaf * (1.0 - res_outtake_twig * adhering_leaf_frac) * frac_cut;

			double to_partition_sap   = 0.0;
			double to_partition_heart = 0.0;

			if (i.cmass_heart >= i.cmass_debt) {
				to_partition_sap   = i.cmass_sap;
				to_partition_heart = i.cmass_heart - i.cmass_debt;
			}
			else {
				to_partition_sap   = i.cmass_sap + i.cmass_heart - i.cmass_debt;
//				dprintf("ATTENTION: pft %s: cmass_debt > cmass_heart; difference=%f\n", (char*)pft.name, i.cmass_debt-i.cmass_heart);
			}
			i.litter_sap += to_partition_sap * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac - harv_eff * stem_frac) * frac_cut;
			i.litter_heart += to_partition_heart * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac - harv_eff * stem_frac) * frac_cut;
		}
		// debt larger than existing wood biomass
		else {
			double debt_excess = i.cmass_debt - (i.cmass_sap + i.cmass_heart);
			dprintf("ATTENTION: cmass_debt > i.cmass_sap + i.cmass_heart; debt_excess=%f\n", debt_excess);
//			i.debt_excess += debt_excess * frac_cut;	// debt_excess currently not dealt with during wood harvest
		}

		// unharvested trees:
		i.cmass_leaf *= (1.0 - frac_cut);
		i.cmass_sap *= (1.0 - frac_cut);
		i.cmass_heart *= (1.0 - frac_cut);
		i.cmass_debt *= (1.0 - frac_cut);

		//Nitrogen:

		harvest = 0.0;

		// harvested products
		harvest += harv_eff * stem_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// harvested products not consumed this year put into harvested_products_slow_nmass
		if (ifslowharvestpool) {
			i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
			harvest = harvest * (1.0 - pft.harvest_slow_frac);
		}

		// harvested products consumed this year put into anflux_harvest
		i.anflux_harvest += harvest;

		residue_outtake = 0.0;

		// removed leaves adhering to twigs
		residue_outtake += res_outtake_twig * adhering_leaf_frac * i.nmass_leaf * frac_cut;

		// removed twigs
		residue_outtake += res_outtake_twig * twig_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// removed coarse roots
		residue_outtake += res_outtake_coarse_root * coarse_root_frac * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// removed residues are oxidised
		i.anflux_harvest += residue_outtake;

		// not removed residues are put into litter
		i.nmass_litter_leaf += i.nmass_leaf * (1.0 - res_outtake_twig * adhering_leaf_frac) * frac_cut;
		i.nmass_litter_sap += i.nmass_sap * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac - harv_eff * stem_frac) * frac_cut;
		i.nmass_litter_heart += i.nmass_heart * (1.0 - res_outtake_twig * twig_frac - res_outtake_coarse_root * coarse_root_frac - harv_eff * stem_frac) * frac_cut;

		// unharvested trees:
		i.nmass_leaf *= (1.0 - frac_cut);
		i.nmass_sap *= (1.0 - frac_cut);
		i.nmass_heart *= (1.0 - frac_cut);
	}
}

/// Harvest function used for managed forest and for clearing natural vegetation at land use change
/** A fraction of trees is cut down (frac_cut)
 *  A fraction of wood is harvested (pft.harv_eff) and returned as acflux_harvest
 *  A fraction of harvested wood (pft.harvest_slow_frac) is returned as harvested_products_slow
 *  The rest, including leaves and roots, is returned as litter.
 *  Called from landcover_dynamics() first day of the year if any natural vegetation is transferred to another land use.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT PARAMETER
 *  \param frac_cut					fraction of trees cut
 *  \param harv_eff					harvest efficiency
 *  \param res_outtake_twig			removed twig fraction
 *  \param res_outtake_coarse_root	removed course root fraction
 *  \param lc_change				whether to save harvest in gridcell-level luc variable
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgC/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - litter_sap   				new sapwood C litter (kgC/m2)
 *   - litter_heart   				new heartwood C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *   - nmass_litter_sap 			new sapwood nitrogen litter (kgN/m2)
 *   - nmass_litter_heart        	new heartwood nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_wood(Individual& indiv, double frac_cut, double harv_eff, double res_outtake_twig, 
		double res_outtake_coarse_root, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_wood(indiv_cp, indiv.pft, indiv.alive, frac_cut, harv_eff, res_outtake_twig, res_outtake_coarse_root);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if (!lc_change) {
		return;
	}

	Stand& stand = indiv.vegetation.patch.stand;
	Landcover& lc = stand.get_gridcell().landcover;
	lc.acflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
	lc.anflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
	if(stand.lc_origin < NLANDCOVERTYPES) {
		lc.acflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		lc.anflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
	}
}

/// Use for normal forest management in calls from growth(). For clearcut during landcover change, use harvest_wood() and kill_remaining_vegetation()
void clearcut(Individual& indiv, double anpp, bool& killed) {

	Patch& patch = indiv.vegetation.patch;
	Patchpft& ppft = patch.pft[indiv.pft.id];

	if (indiv.pft.lifeform == TREE) {

		if(indiv.alive) {
			if(anpp > 0.0)
				ppft.litter_sap += anpp;
			else
				patch.fluxes.report_flux(Fluxes::HARVESTC, anpp);
		}
		harvest_wood(indiv, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake, 0.1); // frac_cut=1, harv_eff=pft.harv_eff, res_outtake_twig=pft.res_outtake, res_outtake_coarse_root=0.1
		indiv.kill();
		indiv.vegetation.killobj();
		killed = true;
	}
}

/// Applies diameter rules if mt.diam_limit is set (default 0) and returns maximum man_strength for this individual. If not set, 1 is returned
double diameter_rules(Individual& indiv) {

	Patch& patch = indiv.vegetation.patch;
	Gridcellst& gst = patch.stand.get_gridcell().st[patch.stand.stid];
	double diam_limit = gst.diam_limit;
	ManagementType& mt = patch.stand.get_current_management();
	if(!diam_limit || mt.secondintervalstart == -1 || patch.age < mt.secondintervalstart)	// Only use diameter limit for continuous cover
		return 1.0;
	double man_strength;
	double diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);
	if(!indiv.height)
		diam = 0.0;

	double diam_max = diam_limit * 2.0;	// Fredrik's continuous management

	if (diam > diam_limit) {
		man_strength = patch.man_strength;
		if(diam > diam_max)
			man_strength = 0.9;	// Fredrik's continuous management
	}
	else {
		man_strength = 0.0;
	}

	return man_strength;
}

/// Distributes man_strength from patch level to individual level rules
/** The management strength * cmass_wood demand is distributed to the individuals according to
 *  in this function and options in the function parameter list. The amount is conserved unless diameter_rules() returns a
 *  maximum man_strength value (individuals with diameters below limit left) and not enough wood cmass is available.
 *  In this case, the diameter limit is lowered by 1% each year until demand fulfilled.
 *
 *  INPUT PARAMETERS
 *  \param select_diam				Whether small (1) or large (2) diameter individuals are preferentially cut, trees above diam_limit only (3).or no preference (0)
 *  \param select_age				Whether young (1) or old (2) individuals are preferentially cut, or no preference (0)
 *  \param select_species			Whether non-selected (1) or selected (2) pft:s are preferentially cut, unselected and selected cutting strengths specified separately (3) or no preference (0)
 *  \param str_unsel				Cutting strength for unselected pft:s if select_species = 3
 *  \param str_sel					String of cutting strengths for selected pft:s if select_species = 3
 */
void distribute_cutting(Patch& patch, int select_diam = 0, int select_age = 0, int select_species = 0, double str_unsel = 0.0, double* str_sel = NULL) {

	if(patch.distributed_cutting)
		return;

	Rank_individuals indiv_class(patch);
	indiv_class.sort_diameter();
//	for(unsigned int i=0;i<patch.vegetation.nobj;i++)
//		dprintf("%f\t", indiv_class.get_diam(i));
//	dprintf("\n");
	ManagementType& mt = patch.stand.get_current_management();
	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	const bool wood_cmass_only = true;
	double cmass_harvest_patch = check_harvest_cmass(patch, wood_cmass_only);
	double cmass_harvest_patch_selection = check_harvest_cmass(patch, wood_cmass_only, true);
	double cmass_harvest_remain = cmass_harvest_patch * patch.man_strength;
	double cmass_harvest_remain_init = cmass_harvest_remain;
	double cmass_harvest_patch_unsel = cmass_harvest_patch - cmass_harvest_patch_selection;

	double* cmass_pft = new double[stand.npftsinselection];
	memset(cmass_pft, 0, sizeof(double)*stand.npftsinselection);
	double* cmass_harvest_remain_pft = new double[stand.npftsinselection];
	memset(cmass_harvest_remain_pft, 0, sizeof(double)*stand.npftsinselection);

	// Determine cmass:harvest for each cohort first for select_age > 0
	double cmass_harvest_ageclass[500] = {0.0};
	double cmass_harvest_ageclass_selection[500] = {0.0};
	for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {
		Individual& indiv = patch.vegetation[i];
		bool pft_selection = mt.pftinselection((const char*)indiv.pft.name) || mt.planting_system == "MONOCULTURE" && mt.pftname == indiv.pft.name || !st.restrictpfts;
		cmass_harvest_ageclass[(int)indiv.age] += check_harvest_cmass(indiv, wood_cmass_only);
		if(pft_selection) {
			cmass_harvest_ageclass_selection[(int)indiv.age] += check_harvest_cmass(indiv, wood_cmass_only);
			if(str_sel) {
				cmass_pft[stand.pft[indiv.pft.id].selection] += check_harvest_cmass(indiv, wood_cmass_only);
				cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] += str_sel[stand.pft[indiv.pft.id].selection] * check_harvest_cmass(indiv, wood_cmass_only);
			}
		}
	}

	// Two laps if selected/unselected species selectively cut

	int nlaps = 1;
	if(select_species)
		nlaps = 2;
	for(int n=0; n<nlaps; n++) {

		double cmass_harvest_remain_init_selection;
		double cmass_harvest_remain_init_ageclass;
		int age_save = -1;
		bool cut_selection;

		if(select_species == 1) {
			if(!n)
				cut_selection = false;
			else
				cut_selection = true;
		}
		else {
			if(!n)
				cut_selection = true;
			else
				cut_selection = false;

			if(select_species == 3) {
				if(!n)
					cmass_harvest_remain = cmass_harvest_patch_selection * patch.man_strength;	// Value used for selection in this case
				else
					cmass_harvest_remain = cmass_harvest_patch_unsel * str_unsel;
			}
		}

		for(unsigned int i = 0; i < patch.vegetation.nobj; i++) {

			if(cmass_harvest_remain < 1e-15 && !(cut_selection && select_species == 3 && str_sel))
				break;

			int index;

			if(select_age == 1)
				index = patch.vegetation.nobj - 1 - i;
			else
				index = i;

			// select_diam overrides select_age
			if(select_diam == 1)
				index = indiv_class.get_index(i);
			else if(select_diam == 2)
				index = indiv_class.get_index(patch.vegetation.nobj -1 - i);

			if(!i)
				cmass_harvest_remain_init_selection = cmass_harvest_remain;

			Individual& indiv = patch.vegetation[index];

			if(indiv.pft.lifeform != TREE)
				continue;

			bool pft_selection = mt.pftinselection((const char*)indiv.pft.name) || mt.planting_system == "MONOCULTURE" && mt.pftname == indiv.pft.name || !st.restrictpfts;

			if(select_species) {
				if(!cut_selection && pft_selection)
					continue;
				if(cut_selection && !pft_selection)
					continue;
			}

			if(!i || (int)indiv.age != age_save)
				cmass_harvest_remain_init_ageclass = cmass_harvest_remain;
			age_save = (int)indiv.age;

			// Harvestable cmass_wood for this individual/cohort
			double cmass_harvest_cohort = check_harvest_cmass(indiv, wood_cmass_only);

			// select_diam overrides select_age
			if(select_diam) {

				// Determine if diameter rules restricts harvestable amount
				double max_cut = diameter_rules(indiv);

				// Satisfy cutting demand by cutting down each cohort starting with thinnest trees first (alternatively thickest trees, select_diam=2);
				// select_diam=3: cut patch.man_strength of individuals with diameter > diam_limit (only in continuous period)
				if(cmass_harvest_cohort) {
					if(select_diam == 3) {
						indiv.man_strength = max_cut;
					}
					else {
						if(cut_selection && select_species == 3 && str_sel) {
							indiv.man_strength = min(1.0, cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] / cmass_harvest_cohort);
							cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] -= indiv.man_strength * cmass_harvest_cohort;
						}
						else {
							indiv.man_strength = min(max_cut, cmass_harvest_remain / cmass_harvest_cohort);	// should use 1.0 instead of max_cut ????
						}
					}
				}
				else {
					indiv.man_strength = 0.0;
				}
			}
			else if(select_age) {

				// Using equal cutting within an age-class
				if(select_species) {
					if(!cut_selection && (cmass_harvest_ageclass[(int)indiv.age] - cmass_harvest_ageclass_selection[(int)indiv.age])) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass / (cmass_harvest_ageclass[(int)indiv.age] - cmass_harvest_ageclass_selection[(int)indiv.age]));
					}
					else if(cut_selection && cmass_harvest_ageclass_selection[(int)indiv.age]) {
						if(select_species == 3 && str_sel) {
							// ony one individual per pft per age !
							indiv.man_strength = min(1.0, cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] / cmass_harvest_cohort);
							cmass_harvest_remain_pft[stand.pft[indiv.pft.id].selection] -= indiv.man_strength * cmass_harvest_cohort;
						}
						else {
							indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass / cmass_harvest_ageclass_selection[(int)indiv.age]);
						}
					}
					else {
						indiv.man_strength = 0.0;
					}
				}
				else {
					if(cmass_harvest_ageclass[(int)indiv.age])
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_ageclass / cmass_harvest_ageclass[(int)indiv.age]);
					else 
						indiv.man_strength = 0.0;
				}
			}
			else {

				if(select_species) {

					// First cut unwanted species by an equal amount
					if(!cut_selection && (cmass_harvest_patch - cmass_harvest_patch_selection)) {
						indiv.man_strength = min(1.0, cmass_harvest_remain_init_selection / (cmass_harvest_patch - cmass_harvest_patch_selection));
					}
					// then cut an equal amount of the pft:s in selection
					else if(cut_selection && cmass_harvest_patch_selection) {
						if(select_species == 3 && str_sel) {
							indiv.man_strength = str_sel[stand.pft[indiv.pft.id].selection];	// Special case: not equal amounts !
						}
						else {
							indiv.man_strength = min(1.0, cmass_harvest_remain_init_selection / cmass_harvest_patch_selection);
						}
					}
					else {
						indiv.man_strength = 0.0;
					}
				}
				else {
					// Cut an equal amount of each cohort
					indiv.man_strength = patch.man_strength;
				}
			}
			cmass_harvest_remain -= indiv.man_strength * cmass_harvest_cohort;
		}
	}
	// If diameter rules used, prescribed cutting may not be acheived (especially in young stands). Try to solve demand by reducing diameter limit by 1% each year when this happens.
	if(mt.secondcutinterval && ((select_diam == 3 && cmass_harvest_remain_init && (cmass_harvest_remain_init - cmass_harvest_remain) < 1e-15)
		|| (select_diam == 1 || select_diam == 2) && cmass_harvest_remain > 1e-15)) {
		dprintf("Year %d: Warning: cmass_harvest_remain = %f, inital cmass_harvest_remain = %f; age = %d\n", date.get_calendar_year(), cmass_harvest_remain,cmass_harvest_remain_init, patch.age);
		// Approximate full rotation age:
		double str_sum = 0.0;

		// Only in continuous period
		if(mt.diam_limit && mt.secondintervalstart > -1 && patch.age >= mt.secondintervalstart) {
			for(int t=0;t<NTHINNINGS;t++) {
				str_sum += mt.thinning_strength[1][t];
			}
			// Take into account number of patches that are cut each year:
			double harvests_per_year = (1.0 * patch.stand.npatch()) / mt.secondcutinterval;

			Gridcellst& gst = patch.stand.get_gridcell().st[patch.stand.stid];
			gst.diam_limit *= (1.0 - (0.01 / harvests_per_year));
			dprintf("New diam_limit = %f\n", gst.diam_limit);
		}
	}

	patch.distributed_cutting = true;

	return;
}

int split_string(char* str) {

	char *p = strtok(str, "\t\n ");
	int count = 0;
	while(p) {
		count++;
		p = strtok(NULL, "\t\n ");
	}

	return count;
}

/// Set forest management intensity for all stands this year
void manage_forests(Gridcell& gridcell) {

	if (!run_landcover || date.day) {
		return;
	}

	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
	
		stand.firstobj();
		while (stand.isobj && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			Patch& patch = stand.getobj();
			if(harvest_secondary_to_new_stand)	{	// avoid when sending harv_cmass to man_frac
				manage_forest(patch);
			}
			stand.nextobj();
		}
		++gc_itr;
	}

	Gridcell::iterator gc_itr2 = gridcell.begin();
	while (gc_itr2 != gridcell.end()) {
		Stand& stand = *gc_itr2;

		stand.firstobj();
		while (stand.isobj && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			Patch& patch = stand.getobj();
			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				bool killed = false;
				harvest_forest(indiv, indiv.pft, indiv.alive, 0.0, killed);

				if(!killed)
					vegetation.nextobj();
			}
			stand.nextobj();
		}
		++gc_itr2;
	}
}

// The following two functions are simplified adaptations (continous cutting) from Swedish forest management code by Fredrik Lagergren and should be
// developed further. Specifically, the productivity values, which should ideally be observed values for each gridcell, are set to a static value.
// Also, the calculated diameter limits and rotation times (which are dependent on productivity) are for Swedish forests.

/// Determines whether this patch should be cut this year.
double manage_forest(Patch& patch) {

	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	ManagementType& mt = stand.get_current_management();

	int first_cutyear = nyear_spinup; // Simulation year when forestry harvesting starts; default is directly after spinup.

	if(st.firstmanageyear < 100000)	// Initialised to 1000000; other values set in instruction file.
		first_cutyear = st.firstmanageyear - date.first_calendar_year;

	if(date.year < first_cutyear || !mt.is_managed())
		return 0.0;

	if(stand.get_current_management().is_managed())
		patch.managed = true;
	else
		return 0.0;

	if(patch.stand.first_year == date.year) {
		patch.plant_this_year = true;	// Forces establishment first stand year to behave like after clearcut
		patch.has_been_cut = true;
	}

	const double minbon = 2.351;	// The minimum average "bonitet" for a county in Sweden
	const double maxbon = 11.311;	// The maximum average "bonitet" for a county in Sweden
	const double bonitet = 10.0;	// Temporary static value (gives cut_int=17)
	double cut_fraction = 0.0;
	double cut_fraction_unsel = 0.0;
	int cut_interval = mt.cutinterval;
	bool clearcut_now = false;
	if(st.cutfirstyear && date.year == first_cutyear)
		clearcut_now = true;

	if(mt.harvest_system == "CLEARCUT") {

		int patch_order = (int)(patch.id * cut_interval * 1.0 / (1.0 * stand.npatch())); //Which year in a cutting interval the patch belongs to

		// clearcut interval set in stand/management type
		if(cut_interval) {
			// Cut according to number of patches and patch id
//			if(!((date.year - first_cutyear - patch_order) % cut_interval))
			// Randomised cut years
//			if(randfrac(stand.seed) < ((double)patch.age / (cut_interval * (cut_interval + 1.0) / 2.0)))
			// Use original patch ages
			if(!(patch.age % cut_interval)) {
				if(patch.age)
					clearcut_now = true;
			}
			else {
				for(int t=0;t<NTHINNINGS;t++) {
					if((mt.thinning_strength[0][t] || mt.thinning_strength_unsel[0][t]) && (patch.age == (int)(cut_interval * mt.thinning_time[0][t]))) {
						cut_fraction = mt.thinning_strength[0][t];
						cut_fraction_unsel = mt.thinning_strength_unsel[0][t];
						patch.man_strength = cut_fraction;
						// Pre-commercial thinning: harvested biomass to litter, youngest cohorts cut first
						if(!t) {
							patch.harvest_to_litter = true;
						}
						distribute_cutting(patch, mt.thinning_select_diam[0][t], mt.thinning_select_age[0][t], mt.thinning_select_pft[0][t]);
					}
				}
			}
		}
		// Use optimal rotation age (mt.cutinterval = 0)
		else {
			// First attempt to calculate optimum rotation age for clearcut 
			if(patch.cmass_wood(true) / max(1,patch.age) > patch.get_tree_cmass_wood_inc_5() && patch.age > 20) {
				clearcut_now = true;
			}
		}
	}

	if(clearcut_now) {
		cut_fraction = 1.0;
		patch.man_strength = cut_fraction;
		patch.age = 0;
		patch.plant_this_year = true;
		patch.clearcut_this_year = true;
	}
	else if(mt.harvest_system == "CONTINUOUS") {

		if(!cut_interval) {
//			cut_interval=30-(int)(15.0*(stand.bonitet-minbon)/(maxbon-minbon));
			cut_interval=30-(int)(15.0*(bonitet-minbon)/(maxbon-minbon));
		}

		int n = 0;	// thinningloop
		int age = patch.age;
		if(mt.secondintervalstart > -1 && patch.age >= mt.secondintervalstart) {
			n = 1;
			cut_interval = mt.secondcutinterval;
			age = patch.age - mt.secondintervalstart;
		}

		int patch_order = (int)(patch.id * cut_interval * 1.0 / (1.0 * stand.npatch()));	// Which year in a cutting interval the patch belongs to

		for(int t=0;t<NTHINNINGS;t++) {
			// Cut according to number of patches and patch id
//			if ((mt.thinning_strength[n][t] || mt.thinning_strength_unsel[n][t]) && (((date.year - first_cutyear - patch_order) % cut_interval) == (int)(cut_interval * mt.thinning_time[n][t]))) {
			// Use original patch ages
			if((mt.thinning_strength[n][t] || mt.thinning_strength_unsel[n][t]) && (age % cut_interval) == (int)(cut_interval * mt.thinning_time[n][t])) {
				cut_fraction = mt.thinning_strength[n][t];
				cut_fraction_unsel = mt.thinning_strength_unsel[n][t];
				patch.man_strength = cut_fraction;
				distribute_cutting(patch, mt.thinning_select_diam[n][t], mt.thinning_select_age[n][t], mt.thinning_select_pft[n][t]);
			}
		}
	}

	return cut_fraction;
}

/// Harvest of tree individuals by an amount man_strength
/*	If clearcut is selected (depending on result from cut_fraction()), individual is killed
 */
void harvest_forest(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed) {

	Patch& patch = indiv.vegetation.patch;
	Patchpft& ppft = patch.pft[indiv.pft.id];

	double man_strength = patch.man_strength;
	if(patch.distributed_cutting)
		man_strength = indiv.man_strength;
	ManagementType& mt = indiv.vegetation.patch.stand.get_current_management();

	if (pft.lifeform==TREE && man_strength > 0.00) {

		double diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);

		if (man_strength == 1.00) {
//			dprintf("Year %d: Clearcut in %s stand, patch %d: %s inidvidual killed\n", date.get_calendar_year(), (char*)stlist[patch.stand.stid].name, patch.id, (char*)indiv.pft.name);
			clearcut(indiv, anpp, killed);
		}
		else {

			double harv_eff_wood_harvest = indiv.pft.harv_eff;				// 0.9
			double res_outtake_twig_wood_harvest = indiv.pft.res_outtake;	// 0.4
			double res_outtake_coarse_root_wood_harvest = 0.1;

			if(patch.harvest_to_litter) {
				harv_eff_wood_harvest = 0.0;
				res_outtake_twig_wood_harvest = 0.0;
				res_outtake_coarse_root_wood_harvest = 0.0;
			}

			if(man_strength) {
				harvest_wood(indiv, man_strength, harv_eff_wood_harvest, res_outtake_twig_wood_harvest, res_outtake_coarse_root_wood_harvest);
				indiv.densindiv *= (1.0 - man_strength);
				if (negligible(indiv.densindiv)) {
					indiv.vegetation.killobj();
					killed=true;
				}
				else {
					allometry(indiv);
				}
			}
		}
		// Will tell the program to skip mortality if management has been performed on this patch,
		patch.managed_this_year = true;		
		patch.has_been_cut = true;
	}
}


/// Harvest function for pasture, representing grazing (previous year).
/*  Function for balancing carbon and nitrogen fluxes from last year's growth
 *  A fraction of leaves is harvested (pft.harv_eff) and returned as acflux_harvest
 *  This represents grazing minus return as manure.
 *  The rest is handled like natural grass in turnover().
 *  Called from growth() last day of the year for normal harvest/grazing.
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_pasture(Harvest_CN& i, Pft& pft, bool alive) {

	double harvest;

	// harvest of leaves (grazing)

	// Carbon:
	harvest = pft.harv_eff * i.cmass_leaf;

	if (ifslowharvestpool) {
		i.harvested_products_slow += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	if (alive)
		i.acflux_harvest += harvest;
	i.cmass_leaf -= harvest;

	// Nitrogen
	// Reduced removal of N relative to C during grazing.
	double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
	harvest = pft.harv_eff * i.nmass_leaf * N_harvest_scale;

	if (ifslowharvestpool) {
		i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	i.anflux_harvest += harvest;
	i.nmass_leaf -= harvest;

	if (grassforcrop && alive) {
		// Carbon:
		double residue_outtake = pft.res_outtake * i.cmass_leaf;	// res_outtake currently set to 0.0,
		i.acflux_harvest += residue_outtake;				// could be used for burning
		i.cmass_leaf -= residue_outtake;

		// Nitrogen:
		residue_outtake = pft.res_outtake * i.nmass_leaf;
		i.anflux_harvest += residue_outtake;
		i.nmass_leaf -= residue_outtake;
	}
}

/// Harvest function for pasture, representing grazing (previous year).
/*  Function for balancing carbon and nitrogen fluxes from last year's growth
 *  A fraction of leaves is harvested (pft.harv_eff) and returned as acflux_harvest
 *  This represents grazing minus return as manure.
 *  The rest is handled like natural grass in turnover().
 *  Called from growth() last day of the year for normal harvest/grazing.
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *  OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_pasture(Individual& indiv, Pft& pft, bool alive, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_pasture(indiv_cp, pft, alive);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if(lc_change) {
		Stand& stand = indiv.vegetation.patch.stand;
		stand.get_gridcell().landcover.acflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		stand.get_gridcell().landcover.acflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		stand.get_gridcell().landcover.anflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		stand.get_gridcell().landcover.anflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
	}
}

/// Harvest function for cropland, including true crops, intercrop grass
/**   and pasture grass grown in cropland.
 *  Function for balancing carbon and nitrogen fluxes from this year's harvested carbon and nitrogen.
 *  A fraction of harvestable organs (grass:leaves) is harvested (pft.harv_eff) and returned as acflux_harvest.
 *  A fraction of leaves is removed (pft.res_outtake) and returned as acflux_harvest
 *  The rest, including roots, is returned as litter, leaving NO carbon or nitrogen in living tissue.
 *  Called from growth() last day of the year for old-style harvest/grazing or, alternatively, from crop_growth_daily() at harvest day
 *	(hdate) or last intercrop day (eicdate).
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function takes a Harvest_CN struct as an input parameter, copied from an individual and it's associated patchpft and patch.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - param nmass_ho				harvestable organ nitrogen biomass (kgC/m2)
 *   - param nmass_agpool			above-ground pool nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_crop(Harvest_CN& i, Pft& pft, bool alive, bool isintercropgrass) {

	double residue_outtake, harvest;

	if (pft.phenology==CROPGREEN) {

	// all root carbon and nitrogen goes to litter
		if (i.cmass_root > 0.0)
			i.litter_root += i.cmass_root;
		i.cmass_root = 0.0;

		if (i.nmass_root > 0.0)
			i.nmass_litter_root += i.nmass_root;
		if (i.nstore_labile > 0.0)
			i.nmass_litter_root += i.nstore_labile;
		if (i.nstore_longterm > 0.0)
			i.nmass_litter_root += i.nstore_longterm;
		i.nmass_root = 0.0;
		i.nstore_labile = 0.0;
		i.nstore_longterm = 0.0;

		// harvest of harvestable organs
		// Carbon:
		if (i.cmass_ho > 0.0) {
			// harvested products
			harvest = pft.harv_eff * i.cmass_ho;

			// not removed harvestable organs are put into litter
			if (pft.aboveground_ho)
				i.litter_leaf += (i.cmass_ho - harvest);
			else
				i.litter_root += (i.cmass_ho - harvest);

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if (ifslowharvestpool) {
				i.harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed (oxidised) this year put into acflux_harvest
			i.acflux_harvest += harvest;
		}
		i.cmass_ho = 0.0;

		// Nitrogen:
		if (i.nmass_ho > 0.0) {

			// harvested products
			harvest = pft.harv_eff * i.nmass_ho;

			// not removed harvestable organs are put into litter
			if (pft.aboveground_ho)
				i.nmass_litter_leaf += (i.nmass_ho - harvest);
			else
				i.nmass_litter_root += (i.nmass_ho - harvest);

			// harvested products not consumed this year put into harvested_products_slow_nmass
			if (ifslowharvestpool) {
				i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed this year put into anflux_harvest
			i.anflux_harvest += harvest;
		}
		i.nmass_ho = 0.0;

		// residues
		// Carbon
		if ((i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem);
			i.acflux_harvest += residue_outtake;

			// not removed residues are put into litter
			i.litter_leaf += i.cmass_leaf + i.cmass_agpool + i.cmass_dead_leaf + i.cmass_stem - residue_outtake;
		}
		i.cmass_leaf = 0.0;
		i.cmass_agpool = 0.0;
		i.cmass_dead_leaf = 0.0;
		i.cmass_stem = 0.0;

		// Nitrogen:
		if ((i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf);
			i.nmass_litter_leaf += i.nmass_leaf + i.nmass_agpool + i.nmass_dead_leaf - residue_outtake;

			// not removed residues are put into litter
			i.anflux_harvest += residue_outtake;
		}
		i.nmass_leaf = 0.0;
		i.nmass_agpool = 0.0;
		i.nmass_dead_leaf = 0.0;
	}
	else if (pft.phenology == ANY) {

		// Intercrop grass
		if (isintercropgrass) {

			// roots

			// all root carbon and nitrogen goes to litter
			if (i.cmass_root > 0.0)
				i.litter_root += i.cmass_root;
			if (i.nmass_root > 0.0)
				i.nmass_litter_root += i.nmass_root;
			if (i.nstore_labile > 0.0)
				i.nmass_litter_root += i.nstore_labile;
			if (i.nstore_longterm > 0.0)
				i.nmass_litter_root += i.nstore_longterm;

			i.cmass_root = 0.0;
			i.nmass_root = 0.0;
			i.nstore_labile = 0.0;
			i.nstore_longterm = 0.0;


			// leaves

			// Carbon:
			if (i.cmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.cmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.litter_leaf += i.cmass_leaf - harvest;

				if (ifslowharvestpool) {
					i.harvested_products_slow += harvest * pft.harvest_slow_frac; // no slow harvest for grass
					harvest = harvest * (1 - pft.harvest_slow_frac);
				}

				i.acflux_harvest += harvest;
			}
			i.cmass_leaf = 0.0;
			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			if (i.nmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.nmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.nmass_litter_leaf += i.nmass_leaf - harvest;

				if (ifslowharvestpool) {
					i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac; // no slow harvest for grass
					harvest = harvest * (1 - pft.harvest_slow_frac);
				}

				i.anflux_harvest += harvest;
			}
			i.nmass_leaf = 0.0;
			i.nmass_ho = 0.0;
			i.nmass_agpool = 0.0;

		}
		else {	// pasture grass

			// harvest of leaves (grazing)

			// Carbon:
			harvest = pft.harv_eff * i.cmass_leaf;

			if (ifslowharvestpool) {
				i.harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}
			if (alive)
				i.acflux_harvest += harvest;
			i.cmass_leaf -= harvest;

			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			// Reduced removal of N relative to C during grazing.
			double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
			harvest = pft.harv_eff * i.nmass_leaf * N_harvest_scale;

			if (ifslowharvestpool) {
				i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}
			i.anflux_harvest += harvest;
			i.nmass_leaf -= harvest;

			i.nmass_ho=0.0;
			i.nmass_agpool=0.0;
		}
	}
}

/// Harvest function for cropland, including true crops, intercrop grass
/**   and pasture grass grown in cropland.
 *  Function for balancing carbon and nitrogen fluxes from this year's harvested carbon and nitrogen.
 *  A fraction of harvestable organs (grass:leaves) is harvested (pft.harv_eff) and returned as acflux_harvest.
 *  A fraction of leaves is removed (pft.res_outtake) and returned as acflux_harvest
 *  The rest, including roots, is returned as litter, leaving NO carbon or nitrogen in living tissue.
 *  Called from growth() last day of the year for old-style harvest/grazing or, alternatively, from crop_growth_daily() at harvest day
 *	(hdate) or last intercrop day (eicdate).
 *  Also called from landcover_dynamics() first day of the year if any natural vegetation
 *    is transferred to another land use.
 *  This calls for a scaling factor, when the pasture area has increased.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop() function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param harvest_grsC				whether harvest daily carbon values are harvested
 *
 *  INPUT/OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - param nmass_ho				harvestable organ nitrogen biomass (kgC/m2)
 *   - param nmass_agpool			above-ground pool nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param indiv					reference to an Individual containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - harvested_products_slow		harvest products to slow pool (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 *   - harvested_products_slow_nmass harvest nitrogen products to slow pool (kgC/m2)
 */
void harvest_crop(Individual& indiv, Pft& pft, bool alive, bool isintercropgrass, bool harvest_grsC) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv, harvest_grsC);

	harvest_crop(indiv_cp, pft, alive, isintercropgrass);

	indiv_cp.copy_to_indiv(indiv, harvest_grsC);

}


/// Transfers all carbon and nitrogen from living tissue to litter
/** Mainly used at land cover change when remaining vegetation after harvest (grass) is
 *   killed by tillage, following an optional burning.
 *
 *  This function takes a Harvest_CN struct as an input parameter, copied from an individual and it's associated patchpft and patch.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param burn						whether above-ground vegetation C & N is sent to the atmosphere
 *								     rather than to litter
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgC/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgC/m2)
 *   - param nmass_ho				harvestable organ nitrogen biomass (kgC/m2)
 *   - param nmass_agpool			above-ground pool nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 */
void kill_remaining_vegetation(Harvest_CN& cp, Pft& pft, bool alive, bool istruecrop_or_intercropgrass, bool burn) {


	if (alive || istruecrop_or_intercropgrass)  {
		cp.litter_root += cp.cmass_root;

		if (burn) {
			cp.acflux_harvest += cp.cmass_leaf;
			cp.acflux_harvest += cp.cmass_sap;
			cp.acflux_harvest += cp.cmass_heart - cp.cmass_debt;
		}
		else {
			cp.litter_leaf += cp.cmass_leaf;
			cp.litter_sap += cp.cmass_sap;
			cp.litter_heart += cp.cmass_heart - cp.cmass_debt;
		}
	}

	cp.nmass_litter_root += cp.nmass_root;
	cp.nmass_litter_root += cp.nstore_longterm;
	cp.nmass_litter_root += cp.nstore_labile;

	if (burn) {
		cp.anflux_harvest += cp.nmass_leaf;
		cp.anflux_harvest += cp.nmass_sap;
		cp.anflux_harvest += cp.nmass_heart;
	}
	else {
		cp.nmass_litter_leaf += cp.nmass_leaf;
		cp.nmass_litter_sap += cp.nmass_sap;
		cp.nmass_litter_heart += cp.nmass_heart;
	}

	if (pft.landcover == CROPLAND) {
		if (pft.aboveground_ho) {
			if (burn) {
				cp.acflux_harvest += cp.cmass_ho;
				cp.anflux_harvest += cp.nmass_ho;
			}
			else {
				cp.litter_leaf += cp.cmass_ho;
				cp.nmass_litter_leaf += cp.nmass_ho;
			}
		}
		else {
			cp.litter_root += cp.cmass_ho;
			cp.nmass_litter_root += cp.nmass_ho;
		}

		if (burn) {
			cp.acflux_harvest += cp.cmass_agpool;
			cp.anflux_harvest += cp.nmass_agpool;
		}
		else {
			cp.litter_leaf += cp.cmass_agpool;
			cp.nmass_litter_leaf += cp.nmass_agpool;
		}
	}

	cp.cmass_leaf = 0.0;
	cp.cmass_root = 0.0;
	cp.cmass_sap = 0.0;
	cp.cmass_heart = 0.0;
	cp.cmass_debt = 0.0;
	cp.cmass_ho = 0.0;
	cp.cmass_agpool = 0.0;
	cp.nmass_leaf = 0.0;
	cp.nmass_root = 0.0;
	cp.nstore_longterm = 0.0;
	cp.nstore_labile = 0.0;
	cp.nmass_sap = 0.0;
	cp.nmass_heart = 0.0;
	cp.nmass_ho = 0.0;
	cp.nmass_agpool = 0.0;

}

/// Transfers all carbon and nitrogen from living tissue to litter
/** Mainly used at land cover change when remaining vegetation after harvest (grass) is
 *   killed by tillage, following an optional burning.
 *
 *  This function copies variables from an individual and it's associated patchpft and patch to
 *  a Harvest_CN struct, which is then passed on to the main harvest_crop() function.
 *  After the execution of the main harvest_crop function, the output variables are copied
 *  back to the individual and patchpft and the patch-level fluxes are updated.
 *
 *  INPUT PARAMETERS
 *  \param alive					whether individual has survived the first year
 *  \param isintercropgrass			whether individual is cover crop grass
 *  \param burn						whether above-ground vegetation C & N is sent to the atmosphere
 *								     rather than to litter
 *  INPUT/OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following indiv-specific public members:
 *   - cmass_leaf 					leaf C biomass (kgC/m2)
 *   - cmass_root					fine root C biomass (kgC/m2)
 *   - cmass_ho						harvestable organ C biomass (kgC/m2)
 *   - cmass_agpool					above-ground pool C biomass (kgC/m2)
 *   - cmass_sap					sapwood C biomass (kgC/m2)
 *   - cmass_heart   				heartwood C biomass (kgC/m2)
 *   - cmass_debt					C "debt" (retrospective storage) (kgC/m2)
 *   - nmass_leaf 					leaf nitrogen biomass (kgN/m2)
 *   - nmass_root 					fine root nitrogen biomass (kgN/m2)
 *   - nmass_sap   					sapwood nitrogen biomass (kgC/m2)
 *   - nmass_heart    				heartwood nitrogen biomass (kgC/m2)
 *   - param nmass_ho				harvestable organ nitrogen biomass (kgC/m2)
 *   - param nmass_agpool			above-ground pool nitrogen biomass (kgC/m2)
 *   - nstore_labile    			labile nitrogen storage (kgC/m2)
 *   - nstore_longterm    			longterm nitrogen storage (kgC/m2)
 *  OUTPUT PARAMETERS
 *  \param Harvest_CN& i			struct containing the following patchpft-specific public members:
 *   - litter_leaf    				new leaf C litter (kgC/m2)
 *   - litter_root 					new root C litter (kgC/m2)
 *   - nmass_litter_leaf 			new leaf nitrogen litter (kgN/m2)
 *   - nmass_litter_root			new root nitrogen litter (kgN/m2)
 *									,and the following patch-level public members:
 *   - acflux_harvest				harvest flux to atmosphere (kgC/m2)
 *   - anflux_harvest   			harvest nitrogen flux out of system (kgC/m2)
 */
void kill_remaining_vegetation(Individual& indiv, bool burn, bool lc_change) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	kill_remaining_vegetation(indiv_cp, indiv.pft, indiv.alive, indiv.istruecrop_or_intercropgrass(), burn);

	indiv_cp.copy_to_indiv(indiv, false, lc_change);

	if (burn && lc_change) {
		Stand& stand = indiv.vegetation.patch.stand;
		Landcover& lc = stand.get_gridcell().landcover;
		lc.acflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		lc.acflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.acflux_harvest / (double)stand.nobj;
		lc.anflux_landuse_change += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
		lc.anflux_landuse_change_lc[stand.lc_origin] += stand.get_gridcell_fraction() * indiv_cp.anflux_harvest / (double)stand.nobj;
	}

}

/// Scaling of last year's or harvest day individual carbon and nitrogen member values in stands that have increased their area fraction this year.
/** Called immediately before harvest functions in growth() or allocation_crop_daily().
 */
void scale_indiv(Individual& indiv, bool scale_grsC) {

	Stand& stand = indiv.vegetation.patch.stand;
	Gridcell& gridcell = stand.get_gridcell();

	if (stand.scale_LC_change >= 1.0) {
		return;
	}

	// Scale individual's C and N mass in stands that have increased in area
	// this year by (old area/new area):
	double scale = stand.scale_LC_change;

	if (scale_grsC) {

		if (indiv.pft.landcover == CROPLAND) {

			if (indiv.has_daily_turnover()) {

				indiv.cropindiv->grs_cmass_leaf -= indiv.cropindiv->grs_cmass_leaf_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_root -= indiv.cropindiv->grs_cmass_root_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_ho -= indiv.cropindiv->grs_cmass_ho_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_agpool -= indiv.cropindiv->grs_cmass_agpool_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_dead_leaf -= indiv.cropindiv->grs_cmass_dead_leaf_luc * (1.0 - scale);
				indiv.cropindiv->grs_cmass_stem -= indiv.cropindiv->grs_cmass_stem_luc * (1.0 - scale);

				indiv.check_C_mass();
			}
			else {
				indiv.cropindiv->grs_cmass_leaf *= scale;
				indiv.cropindiv->grs_cmass_root *= scale;
				indiv.cropindiv->grs_cmass_ho *= scale;
				indiv.cropindiv->grs_cmass_agpool *= scale;
				indiv.cropindiv->grs_cmass_plant *= scale;	//grs_cmass_plant not used
				indiv.cropindiv->grs_cmass_dead_leaf *= scale;
				indiv.cropindiv->grs_cmass_stem *= scale;
			}
		}
	}
	else {

		indiv.cmass_root *= scale;
		indiv.cmass_leaf *= scale;
		indiv.cmass_heart *= scale;
		indiv.cmass_sap *= scale;
		indiv.cmass_debt *= scale;

		if (indiv.pft.landcover == CROPLAND) {
			indiv.cropindiv->cmass_agpool *= scale;
			indiv.cropindiv->cmass_ho *= scale;
		}
	}

	// Deduct individual N present day 0 this year in stands that have increased in area this year, scaled by (1 - old area/new area):
	indiv.nmass_root = indiv.nmass_root - indiv.nmass_root_luc * (1.0 - scale);
	indiv.nmass_leaf = indiv.nmass_leaf - indiv.nmass_leaf_luc * (1.0 - scale);
	indiv.nmass_heart = indiv.nmass_heart - indiv.nmass_heart_luc * (1.0 - scale);
	indiv.nmass_sap = indiv.nmass_sap - indiv.nmass_sap_luc * (1.0 - scale);

	if (indiv.pft.landcover == CROPLAND) {
		indiv.cropindiv->nmass_agpool = indiv.cropindiv->nmass_agpool - indiv.cropindiv->nmass_agpool_luc * (1.0 - scale);
		indiv.cropindiv->nmass_ho = indiv.cropindiv->nmass_ho - indiv.cropindiv->nmass_ho_luc * (1.0 - scale);
		indiv.cropindiv->nmass_dead_leaf =indiv.cropindiv->nmass_dead_leaf - indiv.cropindiv->nmass_dead_leaf_luc * (1.0 - scale);
	}

	if (indiv.nstore_labile > indiv.nstore_labile_luc * (1.0 - scale))
		indiv.nstore_labile -= indiv.nstore_labile_luc * (1.0 - scale);
	else
		indiv.nstore_longterm -= indiv.nstore_labile_luc * (1.0 - scale);
	indiv.nstore_longterm = indiv.nstore_longterm - indiv.nstore_longterm_luc * (1.0 - scale);

	indiv.check_N_mass();
}

/// Yearly function for harvest of all land covers that have yearly allocation, turnover and gridcell.expand_to_new_stand[lc] = false.
/** Should only be called from growth().
//  Harvest functions are preceded by rescaling of living C.
//  Only affects natural stands if gridcell.expand_to_new_stand[NATURAL] is false.
 */
bool harvest_year(Individual& indiv, double anpp) {

	Stand& stand = indiv.vegetation.patch.stand;
	Landcover& landcover = stand.get_gridcell().landcover;
	bool killed = false;

	// Reduce individual's C and N mass in stands that have increased in area this year:
	if (landcover.updated && !indiv.has_daily_turnover()) {
		scale_indiv(indiv, false);
	}

	if (stand.landcover == CROPLAND) {
		if (!indiv.has_daily_turnover())
			harvest_crop(indiv, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass, false);
	}
	else if (stand.landcover == PASTURE) {
		harvest_pasture(indiv, indiv.pft, indiv.alive);
	}


	return killed;
}

/// Yield function for true crops and intercrop grass.
void yield_crop(Individual& indiv) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());

	if (indiv.pft.phenology == ANY) {			// grass intercrop yield

		// Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
		if (cropindiv.ycmass_leaf > 0.0)
			cropindiv.yield = cropindiv.ycmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		else
			cropindiv.yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		if (cropindiv.harv_cmass_leaf > 0.0)
			cropindiv.harv_yield = cropindiv.harv_cmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		else
			cropindiv.harv_yield = 0.0;
	}
	else if (indiv.pft.phenology == CROPGREEN) {		//true crop yield

		// Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
		if (cropindiv.ycmass_ho > 0.0)
			cropindiv.yield = cropindiv.ycmass_ho * indiv.pft.harv_eff * 2.0;// Should be /0.446 instead
		else
			cropindiv.yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		if (cropindiv.harv_cmass_ho > 0.0)
			cropindiv.harv_yield=cropindiv.harv_cmass_ho * indiv.pft.harv_eff * 2.0;
		else
			cropindiv.harv_yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		for (int i=0;i<2;i++) {
			if (cropindiv.cmass_ho_harvest[i] > 0.0)
				cropindiv.yield_harvest[i] = cropindiv.cmass_ho_harvest[i] * indiv.pft.harv_eff * 2.0;
			else
				cropindiv.yield_harvest[i]=0.0;
		}
	}

	return;
}

/// Yield function for pasture grass grown in cropland landcover
void yield_pasture(Individual& indiv, double cmass_leaf_inc) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());

	// Normal CC3G/CC4G stand growth (Pasture)

	// OK if turnover_leaf==1.0, else (cmass_leaf+cmass_leaf_inc)*indiv.pft.harv_eff*2.0
	cropindiv.yield = max(0.0, cmass_leaf_inc) * indiv.pft.harv_eff * 2.0;
	// Although no specified harvest date, harv_yield is set for compatibility.
	cropindiv.harv_yield = cropindiv.yield;
}

/// Function that determines amount of nitrogen applied today. Crop-specific, pft-based.
void nfert_crop(Patch& patch) {

	Gridcell& gridcell = patch.stand.get_gridcell();

	patch.dnfert = 0.0;

	pftlist.firstobj();
	// Loop through PFTs
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if (patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
			if(!ppftcrop.growingseason) {
				pftlist.nextobj();
				continue;
			}

			double nfert = pft.N_appfert;
			double mineral = 1.0;
			if (gridcellpft.Nfert_read >= 0.0) {
				nfert = gridcellpft.Nfert_read;
				if (gridcellpft.Nfert_man_read > 0.0) {
					mineral = 1.0 - gridcellpft.Nfert_man_read;
				}

			}
			if (!ppftcrop.fertilised[0] && ppftcrop.dev_stage > 0.0) {
				// Fertiliser application at dev_stage = 0, sowing.
				patch.dnfert = nfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * (1.0 - pft.fertrate[0] - pft.fertrate[1]));
				ppftcrop.fertilised[0] = true;
				if (mineral<1.0) {
					patch.soil.sompool[SOILMETA].nmass+= nfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILMETA].cmass+=nfert * (1.0 - mineral) * 30.0 * 0.25;
					patch.soil.sompool[SOILSTRUCT].nmass+= nfert * (1.0 - mineral) * 0.5;
					patch.soil.sompool[SOILSTRUCT].cmass+= nfert * (1.0 - mineral) * 30.0 * 0.75;
					patch.fluxes.report_flux(Fluxes::MANUREC,-nfert * (1.0 - mineral) * 30.0 );
					patch.anfert += nfert * (1.0 - mineral);
					patch.fluxes.report_flux(Fluxes::MANUREN,nfert * (1.0 - mineral));
				}
			}
			else if (!ppftcrop.fertilised[1] && ppftcrop.dev_stage > pft.fert_stages[0]) {
				patch.dnfert = nfert * mineral * pft.fertrate[0];
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * pft.fertrate[0]);
				ppftcrop.fertilised[1] = true;
			}
			else if (!ppftcrop.fertilised[2] && ppftcrop.dev_stage > pft.fert_stages[1]) {
				patch.dnfert = nfert * mineral * (pft.fertrate[1]);
				patch.fluxes.report_flux(Fluxes::NFERT,nfert * mineral * pft.fertrate[1]);
				ppftcrop.fertilised[2] = true;
			}
		}
		pftlist.nextobj();
	}
	patch.anfert += patch.dnfert;
}

/// Function that determines amount of nitrogen applied today.
void nfert(Patch& patch) {

	Stand& stand = patch.stand;
	StandType& st = stlist[stand.stid];
	Gridcell& gridcell = stand.get_gridcell();

	if(stand.landcover == CROPLAND) {
		nfert_crop(patch);
		return;
	}

	// General code for applying nitrogen to other land cover types, an equal amount every day.
	double nfert;
	if(gridcell.st[st.id].nfert >= 0.0) {	// todo: management type variable (mt.nfert)
		nfert = gridcell.st[st.id].nfert;
	}
	else {
		nfert = 0.0;
	}
	patch.dnfert = nfert / date.year_length();
	patch.anfert += patch.dnfert;
}

// Updates forest rotation status
/** Sets new forest management variables by calling stand.rotate() on st.mtstartyear[m]
 */
void forest_rotation(Stand& stand) {

	StandType& st = stlist[stand.stid];

	if(stand.landcover != FOREST || st.rotation.ncrops < 2)
		return;

	stand.nyears_inrotation++;

	for(int m=0;m<NROTATIONPERIODS_MAX;m++) {
		if(st.mtstartyear[m] == date.get_calendar_year()) {
			stand.rotate(m);
			break;
		}
	}
}

/// Updates crop rotation status
/** Sets new crop management variables, typically on harvest day
 */
void crop_rotation(Stand& stand) {

	if (stand.landcover != CROPLAND) {
		return;
	}

	CropRotation& rotation = stlist[stand.stid].rotation;

	stand.ndays_inrotation++;

	if (rotation.ncrops < 2 || !stand.isrotationday) {
		return;
	}

	int firstrotyear = rotation.firstrotyear - date.first_calendar_year;
	bool postpone_rotation = false;

	// Alternative uses of firstrotyear:
	// 1. Before firstrotyear, grow only crop1:
//	if(date.year < firstrotyear)
//		postpone_rotation = true;

	// 2. Synchronise rotation with firstrotyear:

	// A. At the creation of the stand:
	if (date.year < stand.first_year + 3)
	// B. At firstrotyear
//		if(date.year == firstrotyear - 1)
	// C. Continuously:
	{
		if ((abs(firstrotyear - date.year) % rotation.ncrops) != stand.current_rot)
			postpone_rotation = true;
	}

	if (!postpone_rotation) {

		if(stand.infallow) {
			stand.infallow = false;
			stand.get_gridcell().pft[stand.pftid].sowing_restriction = false;
		}

		int old_pftid = stand.pftid;

		stand.rotate();

		for (unsigned int p=0; p<stand.nobj; p++) {

			cropphen_struct& previous = *(stand[p].pft[old_pftid].get_cropphen());
			cropphen_struct& current = *(stand[p].pft[stand.pftid].get_cropphen());

			previous.bicdate = -1;
			if(!previous.intercropseason)
				current.bicdate = stepfromdate(date.day, 15);
			previous.eicdate = -1;
			current.eicdate = -1;
			previous.hdate = -1;
			current.intercropseason = previous.intercropseason;
		}

		// Adds sowing and harvest dates for the second crop in a double cropping system
		if (rotation.multicrop && rotation.ncrops == 2 && stand.current_rot == 1) {
			if (stand.pft[stand.pftid].sdate_force < 0)
				stand.pft[stand.pftid].sdate_force = stepfromdate(date.day, 10);
			if (stand.pft[stand.pftid].hdate_force < 0) {
				stand.pft[stand.pftid].hdate_force = stepfromdate(stand.pft[old_pftid].sdate_force, -10);
				}
		}

		if(stand.get_current_management().fallow) {
			stand.infallow = true;
			stand.get_gridcell().pft[stand.pftid].sowing_restriction = true;
		}
	}

	stand.isrotationday = false;
}
