////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file storm.cpp
/// \Calculates storm sensitivity and damaged fraction of individuals, which is also scaled to stand and gridcell
/// \Equations are taken from Lagergren et al 2012, Ecological Modelling 247:71-82
/// \The calculation of the "Frozen soil index" is placed in "dailyaccounting_patch" in driver.cpp
/// \author Fredrik Lagergren
/// $Date:  $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "landcover.h"
#include "storm.h"
#include "driver.h"

void storm(Gridcell& gridcell) {

	//if (date.day) {
	if (!(date.islastday && date.islastmonth)) {
		return;
	}

	const double a = 0.3; //Parameter a (Eq 10), sensitivity of a unthinned stand relative to a recently thinned
	const double b = 0.35; //Parameter b (Eq 10), thinning strength at which RSI = 1 directly after thinning

	const double kh = 5.0; // Height (m) below which no wind damage occurs (Eq 13)
	const double kd = 12.0; // Parameter (cm) that reduce the dependency of diameter (Eq 13)

	const double eip_min = 0.25; // Minimum level of the random factor for the patch level exposure to different wind direction (Eq 8)
	// Set to 1.0 for no random variation in exposure to different wind directions

	const double h_ref_p = 30.0; // Reference height in height index calculation at patch level (Eq 9)
	const double h_add_p = 0.5; // Additive term in height index calculation at patch level (Eq 9)

	const double h_ref_c = 30.0; // Reference height in height index calculation at cohort level (Eq 12)
	const double h_add_c = 1.0; // Additive term in height index calculation at cohort level (Eq 12)
	
	// Constant wind load (Eq 6 in Lagergren et al)
	double wind_load[6] = {0.0};
	double wind_load_pot[6];

	// Evenly distributed wind, gives less random results when few patches are simulated
	wind_load_pot[0] = 0.05;
	wind_load_pot[1] = 0.05;
	wind_load_pot[2] = 0.05;
	wind_load_pot[3] = 0.05;
	wind_load_pot[4] = 0.05;
	wind_load_pot[5] = 0.05;
	/*
		// All wind from one sector, gives a more realistic distribution of damage among patches
		wind_load[0] = 0.3;
		wind_load[1] = 0.0;
		wind_load[2] = 0.0;
		wind_load[3] = 0.0;
		wind_load[4] = 0.0;
		wind_load[5] = 0.0;
	*/

	gridcell.storm_prob_force = false;

	if (readstorm_prob) {
		if(gridcell.storm_prob_gc > randfrac(gridcell.storm_seed)) {
			for (int w = 0; w < 6; w++)
				wind_load[w] = wind_load_pot[w];
		}
	}
	else {
		for (int w = 0; w < 6; w++)
			wind_load[w] = gridcell.windload / 6.0;	
	}

	double d2hsum, d2sum; // Sums for calculation of BA-weighted patch height, section 2.6 in Lagergren et al
	double diam; // Diameter (cm)
	double forest_frac_gridcell = gridcell.landcover.frac[NATURAL] + gridcell.landcover.frac[FOREST];

	// Mean height over stands and patches
	double height = 0.0;
	gridcell.storm_dam_wood_prev = gridcell.storm_dam_wood;
	gridcell.storm_dam_wood = 0.0;

	// Set default values
	// Mean treatment area in gridcell (ha), retrived from ins-file, default = 3 ha
	double obj_size = 3.0;
	// Prescribed average height in gridcell, retrived from insfile, if set to 0.0 calculated height is used
	double height_p = 0.0; // 0.0 means that calculated mean height is used instead
	// Todo, enable reading of the parameter by gridcell from external file

	// Exposure index at gridcell level, Eq. 7 in Lagergren et al
	const double k_forfrac = 0.8; // Factor to reduce sensitivity with 100% cover of forest in landscape
	const double k_factor = 2.0; //Scaling factor
	// Exposure index for the forest in the gridcell
	double exp_ind = (1 - k_forfrac * forest_frac_gridcell) * k_factor / sqrt(obj_size);

	// FL: NATURAL is included in the calculation 20 lines above, but the code could be simplified

	// First loops, used to calculate average height for gridcells and patches
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		// START OF LOOP THROUGH STANDS
		Stand& stand = *gc_itr;

		stand.firstobj();
		while (stand.isobj && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			// START OF LOOP THROUGH PATCHES
			Patch& patch = stand.getobj();

			// Storm. Updating the thinning information used for root stability index before resetting
			if (patch.managed_this_year && patch.man_strength < 1.0) {//Was if(patch.managed_previous_year)
				patch.thinned = true;
				patch.root_m_th = patch.root_mass;
			}

			// Called the first year or for a new stand, for calculating fixed exposure indices at patch level
			if (date.year == 0 || !stand.storm_init_done) {
				// Exposure index at patch level from 6 different wind sectors, Eq. 8 in Lagergren et al
				for (int wind_sector = 0; wind_sector < 6; wind_sector++)
//					patch.exp_ind[wind_sector] = exp_ind * (eip_min + (1-eip_min) * randfrac(patch.stand.seed));
					patch.exp_ind[wind_sector] = exp_ind * (eip_min + (1-eip_min) * 0.5);	// Stochastic exposure indices turned off (temporarily ?)
			}

			d2hsum = 0.0;
			d2sum = 0.0;

			double old_root_mass = patch.root_mass;
			patch.root_mass = 0.0;

			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			// START OF LOOP THROUGH INDIVIDUALS
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				if (indiv.pft.lifeform == TREE) {
					diam = indiv.diam * 100.0;
					d2sum += diam*diam*indiv.densindiv;
					d2hsum += diam*diam*indiv.height*indiv.pft.eff_h*indiv.densindiv;
					patch.root_mass += indiv.cmass_root;
				}

				vegetation.nextobj();
			}	// End of loop through individuals
			
			if (patch.managed_this_year && patch.man_strength < 1.0) {//Was if(patch.managed_previous_year)
				// Calculation of the most recent thinning strength from the reduction in root biomass
				patch.strength_b = min(1.0,1.0 - patch.root_mass/old_root_mass);
			}

			if (patch.strength_b < 0.1) {
				patch.thinned = false;
				patch.strength_b = 0.0;
			}
			if (patch.root_mass >= patch.root_m_th)
				patch.thinned = false;


			if (d2sum > 0.0)
				patch.height = d2hsum/d2sum; //BA-weighted mean height
			else
				patch.height = 0.0;

			height += patch.height / (double)stand.npatch() * stand.get_gridcell_fraction();
			
			stand.nextobj();
		}	// End of loop through patches
		stand.storm_init_done = true;

		++gc_itr;
	}	// End of loop through stands

	height /= forest_frac_gridcell;

	// Use of prescribed or calculated average height in gridcell
	double gridcellh;
	if (height_p == 0.0)
		gridcellh = height;
	else
		gridcellh = height_p;

	Gridcell::iterator gc_itr2 = gridcell.begin();
	while (gc_itr2 != gridcell.end()) {
		// START OF LOOP THROUGH STANDS
		Stand& stand = *gc_itr2;

		stand.firstobj();

		while (stand.isobj && (stand.landcover == FOREST || stand.landcover == NATURAL)) {
			// START OF LOOP THROUGH PATCHES
			Patch& patch = stand.getobj();

			patch.height_ind = max(0.0,(patch.height - gridcellh)/h_ref_p + h_add_p); // Eq. 9 in Lagergren et al

			// Root stability index. Eq. 10 in Lagergren et al
			if (patch.thinned) {
				patch.root_stab_ind = a + (1 - a)/b * (1 - patch.root_mass/patch.root_m_th);
			}
			else {
				patch.root_stab_ind = a;
			}

			Vegetation& vegetation = patch.vegetation;
			vegetation.firstobj();
			// START OF LOOP THROUGH INDIVIDUALS
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();
				if (indiv.pft.lifeform == TREE) {
					indiv.height_ind = max(0.0,(indiv.height - patch.height)/h_ref_c + h_add_c); // Eq. 12 in Lagergren et al
						
					// Allometric index according to Eq. 13 of Lagergren el al
					diam = indiv.diam * 100.0;
					indiv.allom_ind = max(0.0,(indiv.height - kh)/(diam + kd) * indiv.pft.storm_sensiv);
					indiv.dam_frac = 0.0;
					indiv.dam_frac_pot = 0.0;
					indiv.sens_ind_tot = 0.0;
						
					for (int wind_sector = 0;wind_sector<6;wind_sector++) {

						// Total sensitivity for a cohort, Eq 14 in Lagergren et al
						indiv.sens_ind[wind_sector] = patch.exp_ind[wind_sector] * patch.height_ind * patch.root_stab_ind *
							indiv.height_ind * indiv.allom_ind;
						indiv.sens_ind_tot += indiv.sens_ind[wind_sector];

						// Damaged fraction of a cohort, Eq 15 in Lagergren et al
						indiv.dam_frac += indiv.sens_ind[wind_sector] * wind_load[wind_sector] * patch.soil.frozen_soil_ind * m_stormfactor;
						indiv.dam_frac_pot += indiv.sens_ind[wind_sector] * wind_load_pot[wind_sector] * patch.soil.frozen_soil_ind * m_stormfactor;

						//if (date.get_calendar_year() > 1799)
							//dprintf("Year %d: sector si = %f, cmwood = %f, df = %f cmwood_p = %f\n", date.get_calendar_year(), indiv.sens_ind[wind_sector], indiv.cmass_wood(), indiv.dam_frac, patch.cmass_wood());
					}

					indiv.dam_frac = min(1.0, indiv.dam_frac);			// Making sure that not more than 100% is damaged
					indiv.dam_frac_pot = min(1.0, indiv.dam_frac_pot);	// Making sure that not more than 100% is damaged
				}

				vegetation.nextobj();
			}	// End of loop through individuals

			stand.nextobj();
		}	// End of loop through patches

		++gc_itr2;
	}	// End of loop through stands
}
