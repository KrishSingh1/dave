///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.cpp
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// Landcover change, crop and pasture definitions.
///
/// \author Mats Lindeskog,
/// \based on LPJ-mL C++ code received from Alberte Bondeau in 2008.
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "landcover.h"
#include "canexch.h"
#include "guessmath.h"

#define MAXHUTEMP					//30 degree limit for heat unit summation
#define SD_TEMP_WINDOW				//Uses sowing window for temperature-dependent sowing.
#define IRRIGATED_USE_TEMP_SDATE	//Use temperature-dependent sowing date for irrigated crops at site with PRECTEMP seasonality.
//#define LOW_SOWING_TEMPERATURE_LIMIT	//Sowing not allowed when temperature is always below sowing limit. Intercrop grass groen instead year through.
#define HIGH_SOWING_TEMPERATURE_LIMIT	//Sowing not allowed when mean temperature is above limit (TeWW).
//#define DELAYED_SEEDCARBON		//Seed carbon allocation to leaves and roots are done over a 10-day period.
//#define GRASS_SEED_CMASS	// Carbon allocated to grass on bicdate or the day after turnover.

/////////////////////////// Functions facilitating handling time periods spanning newyear //////////////////////////////////////

/// Query whether a date is within a period spanned by two dates.
bool dayinperiod(int day, int start, int end) {

	bool acrossnewyear = false;

	if(start  <0 || end < 0)	// a negative value should not be a valid day
		return false;

	if(start > end)
		acrossnewyear = true;

	if(day >= start && day <= end && !acrossnewyear || (day >= start || day <= end) && acrossnewyear)
		return true;
	else
		return false;
}

/// Step n days from a date.
int stepfromdate(int day, int step) {

	if(day + step > 0)
		return (day + step) % 365;
	else if(day + step < 0)
		return day + step + 365;
	else
		return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////  Landcover stand dynamics and C-partitioning  /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creation of stands when run_landcover==true
void landcover_init(Gridcell& gridcell, InputModule* input_module) {

	// get landcover and crop area fractions from landcover input file(s) or ins-file.
	input_module->getlandcover(gridcell);

	// create stands for landcovers with only one (initial) stand
	for(int i=0; i<NLANDCOVERTYPES; i++) { //For all landcover types without subclasses
		if(i != CROPLAND) {				
			if(run[i]) {
				if(gridcell.landcoverfrac[i] > 0.0) {
					gridcell.create_stand_lu((landcovertype)i, gridcell.landcoverfrac[i]);
				}
			}
		}
	}

	if(run[CROPLAND]) {
		// create crop stands
		if(gridcell.landcoverfrac[CROPLAND] > 0.0)
		{
			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft = pftlist.getobj();
				if(pft.landcover == CROPLAND && pft.cftid >= 0) {
					if(gridcell.cftfrac[pft.cftid] > 0.0) {
						gridcell.create_stand_lu(CROPLAND, gridcell.cftfrac[pft.cftid] * gridcell.landcoverfrac[CROPLAND], pft.cftid);
					}
				}
				pftlist.nextobj();
			}
		}
	}
}

/// landcover_change_transfer constructor
landcover_change_transfer::landcover_change_transfer() {

	transfer_litter_leaf = transfer_litter_sap = transfer_litter_heart = transfer_litter_root = transfer_litter_repr = transfer_harvested_products_slow = NULL;
	transfer_nmass_litter_leaf = transfer_nmass_litter_sap = transfer_nmass_litter_heart = transfer_nmass_litter_root = transfer_harvested_products_slow_nmass = NULL;

	transfer_acflux_harvest = transfer_anflux_harvest = transfer_cpool_fast = transfer_cpool_slow = transfer_wcont_evap = transfer_decomp_litter_mean = 0.0;
	transfer_k_soilfast_mean = transfer_k_soilslow_mean = transfer_nmass_avail = transfer_snowpack = transfer_snowpack_nmass = 0.0;

	memset(transfer_wcont,0,NSOILLAYER*sizeof(double));

	for(int i=0; i<NSOMPOOL; i++)
		transfer_sompool[i].ntoc = 0.0;

}

/// landcover_change_transfer deconstructor
landcover_change_transfer::~landcover_change_transfer() {

	if(transfer_litter_leaf) delete[] transfer_litter_leaf;
	if(transfer_litter_sap) delete[] transfer_litter_sap;
	if(transfer_litter_heart) delete[] transfer_litter_heart;
	if(transfer_litter_root) delete[] transfer_litter_root;
	if(transfer_litter_repr) delete[] transfer_litter_repr;
	if(transfer_harvested_products_slow) delete[] transfer_harvested_products_slow;
	if(transfer_nmass_litter_leaf) delete[] transfer_nmass_litter_leaf;
	if(transfer_nmass_litter_sap) delete[] transfer_nmass_litter_sap;
	if(transfer_nmass_litter_heart) delete[] transfer_nmass_litter_heart;
	if(transfer_nmass_litter_root) delete[] transfer_nmass_litter_root;
	if(transfer_harvested_products_slow_nmass) delete[] transfer_harvested_products_slow_nmass;
}

/// allocates memory for landcover_change_transfer object
void landcover_change_transfer::allocate() {

		transfer_litter_leaf = new double[npft];
		transfer_litter_sap = new double[npft];
		transfer_litter_heart = new double[npft];
		transfer_litter_root = new double[npft];
		transfer_litter_repr = new double[npft];
		transfer_harvested_products_slow = new double[npft];

		transfer_nmass_litter_leaf = new double[npft];
		transfer_nmass_litter_sap = new double[npft];
		transfer_nmass_litter_heart = new double[npft];
		transfer_nmass_litter_root = new double[npft];
		transfer_harvested_products_slow_nmass = new double[npft];

		memset(transfer_litter_leaf, 0, sizeof(double) * npft);
		memset(transfer_litter_sap, 0, sizeof(double) * npft);
		memset(transfer_litter_heart, 0, sizeof(double) * npft);
		memset(transfer_litter_root, 0, sizeof(double) * npft);
		memset(transfer_litter_repr, 0, sizeof(double) * npft);
		memset(transfer_harvested_products_slow, 0, sizeof(double) * npft);

		memset(transfer_nmass_litter_leaf,0,sizeof(double)*npft);
		memset(transfer_nmass_litter_sap,0,sizeof(double)*npft);
		memset(transfer_nmass_litter_heart,0,sizeof(double)*npft);
		memset(transfer_nmass_litter_root,0,sizeof(double)*npft);
		memset(transfer_harvested_products_slow_nmass,0,sizeof(double)*npft);
}


/// Gets this year's landcover and crop area fractions, checks that area changes are significant and that the net changes are zero.
/** Stores changes in area fractions for the landcovers and crops.
 *
 *  OUTPUT PARAMETERS
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param cropstand_change					array with this year's difference in area fractions of the different crop stands
 *  \param changeLC							sum of all stands' absolute changes
 *  \param change_crop						sum of all crop stands' absolute changes
 *  \param receiving_fraction				sum of added area to expanding stands
 *  \param LCchangeCtransfer				whether to transfer carbon, nitrogen and water of reduced stands to expanding stands
 */
bool checkLCchange(Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], double cropstand_change[NCROPSTANDS_MAX], 
				   double& changeLC, double& change_crop, double& receiving_fraction, bool& LCchangeCtransfer, InputModule* input_module) {

	double cropfrac_change[NCROPSTANDS_MAX];
	double cropfrac_sum_old = 0.0;
	double change_stand = 0.0;
	double transferred_fraction = 0.0;

	memset(cropfrac_change, 0, NCROPSTANDS_MAX * sizeof(double));

	//Save old fraction values:									
	for(int i=0; i<NLANDCOVERTYPES; i++)
		gridcell.landcoverfrac_old[i] = gridcell.landcoverfrac[i];
	for(int i=0; i<NCROPSTANDS_MAX; i++)
		cropfrac_sum_old += gridcell.cftfrac_old[i] = gridcell.cftfrac[i];

	//Get new gridcell.landcoverfrac and/or gridcell.cftfrac from LUdata and CFTdata.		
	input_module->getlandcover(gridcell);	

	if(!lcfrac_fixed) {
		for(int i=0; i<NLANDCOVERTYPES; i++) {
			landcoverfrac_change[i] = gridcell.landcoverfrac[i] - gridcell.landcoverfrac_old[i];
			changeLC += fabs(landcoverfrac_change[i]) / 2.0;
			if(i != CROPLAND) {
				if(landcoverfrac_change[i] < 0.0)
					transferred_fraction -= landcoverfrac_change[i];
				if(landcoverfrac_change[i] > 0.0)
					receiving_fraction += landcoverfrac_change[i];
				change_stand += fabs(landcoverfrac_change[i]) / 2.0;
			}
		}
	}

	if(run[CROPLAND] && (!cftfrac_fixed || !lcfrac_fixed)) {
		for(int i=0;i<NCROPSTANDS_MAX;i++) {
			cropfrac_change[i] = gridcell.cftfrac[i] - gridcell.cftfrac_old[i];
			cropstand_change[i] = gridcell.cftfrac[i] * gridcell.landcoverfrac[CROPLAND] - gridcell.cftfrac_old[i] * gridcell.landcoverfrac_old[CROPLAND];

			if(cropstand_change[i] < 0.0)
				transferred_fraction -= cropstand_change[i];
			if(cropstand_change[i] > 0.0)
				receiving_fraction += cropstand_change[i];

			if(cropfrac_sum_old != 0.0) {
				change_crop += fabs(cropfrac_change[i]) / 2.0;
			}
			else {
				change_crop += fabs(cropfrac_change[i]);
			}
			change_stand += fabs(cropstand_change[i]) / 2.0;	//cropfrac_sum_old+gridcell.landcoverfrac[NATURAL] should never be 0.0
		}
	}

	// if no changes, do nothing.
	if(changeLC < 0.00001 && change_crop < 0.00001) {
		return false;
	}
	// check for balance of reduced and increased stand fractions
	else {
		if(fabs(transferred_fraction - receiving_fraction) > 0.0001 || fabs(change_stand-receiving_fraction) > 0.0001) {
			if(run[CROPLAND] && run[NATURAL] && run[PASTURE]) {
				// end program if balance is expected (no landcovers inactivated)
				fail("Transferred landcover fractions not balanced !\n");
			}
			else {
				// allow program to continue, but inactivate landcover change mass transfer 
				LCchangeCtransfer = false;
				if(!SUPPRESSLARGEOUTPUT)
					dprintf("Transferred landcover fractions not balanced !\nLandcover change carbon flux not calculated.\n");
			}
		}
		return true;
	}
}


/// identifies which natural stands to reduce in area, sets new fractions for these and sets nnaturalstands
/** Updates stand.frac. for reduced natural stands
 *  Stores carbon, nitrogen and water of harvested area in a temporary struct.
 *  Should be followed by a call to stand_dynamics() to kill stands with a new area of 0 
 *  Do not call from loop with call to gridcell.nextobj.
 *
 *  OUTPUT PARAMETERS
 *
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param nnaturalstands					number of natural stands in the gridcell
 */
void reduce_natural_stands(Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], int& nnaturalstands) {

	gridcell.firstobj();
	while (gridcell.isobj) {
		Stand& stand=gridcell.getobj();

		if(stand.landcover == NATURAL) {
			stand.natural_frac_change = 0.0;
			nnaturalstands++;
		}

		gridcell.nextobj();
	}

	if(nnaturalstands > 1  && landcoverfrac_change[NATURAL] < 0.0) {
		double natural_change_remain = landcoverfrac_change[NATURAL];

		bool reduce_all_stands = false;	//convert equal percentage of area from all stands
		bool young_stands_first = true;	//convert area from youngest stands first

// Remove this section and always use the default values ?
		// convert equal percentage of area from all stands if both managed forest and other managed land expands
		if((landcoverfrac_change[CROPLAND] > 0.0 || landcoverfrac_change[PASTURE] > 0.0 || landcoverfrac_change[URBAN] > 0.0 || landcoverfrac_change[PEATLAND] > 0.0) 
				&& landcoverfrac_change[FOREST] > 0.0) {
			reduce_all_stands=true;
			young_stands_first=false;
		}
		// convert area from youngest stands first if managedforest does not expand
		else if(landcoverfrac_change[CROPLAND]>0.0 || landcoverfrac_change[PASTURE]>0.0 || landcoverfrac_change[URBAN]>0.0 || landcoverfrac_change[PEATLAND]>0.0)
			young_stands_first=true;
		// convert area from oldest stands first if only managed forest expands
		else if(landcoverfrac_change[FOREST]>0.0)	
			young_stands_first=false;
//////////////////////

		for(unsigned int i = 0; i < gridcell.nobj; i++) {
			int index;

			if(young_stands_first)
				index = gridcell.nobj - 1 - i;
			else
				index = i;

			Stand& stand = gridcell[index];	

			if(stand.landcover == NATURAL) {
				// convert equal areas from all stands
				if(reduce_all_stands) {
					stand.natural_frac_change = landcoverfrac_change[NATURAL] * stand.get_gridcell_fraction() / gridcell.landcoverfrac_old[NATURAL];
					stand.set_gridcell_fraction(stand.get_gridcell_fraction() + stand.natural_frac_change);
				}
				else {		
					if(stand.get_gridcell_fraction() > 0.0) {
						//all natural landcover decrease is taken from this stand
						if(stand.get_gridcell_fraction() >= -natural_change_remain) {
							stand.natural_frac_change = natural_change_remain;
							stand.set_gridcell_fraction(stand.get_gridcell_fraction() + stand.natural_frac_change);
							natural_change_remain = 0.0;
							break;
						}
						//more stands will have to be reduced
						else {						
							stand.natural_frac_change = -stand.get_gridcell_fraction();
							natural_change_remain += stand.get_gridcell_fraction();
							stand.set_gridcell_fraction(0.0);	//will be killed below
						}				
					}
				}
			}
		}
	}		
}


 /// Handles harvest and turnover of reduced stands at landcover change.
/** Updates stand.frac.
 *  Sets LC_updated to true
 *  Stores carbon, nitrogen and water of harvested area in a temporary struct.
 *  Should be followed by a call to stand_dynamics() to kill stands with a new area of 0 
 *  Do not call from loop with call to gridcell.nextobj.
 *
 *  INPUT PARAMETERS
 *
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param cropstand_change					array with this year's difference in area fractions of the different crop stands
 *  \param receiving_fraction				sum of added area to expanding stands
 *  \param nnaturalstands					number of natural stands in the gridcell
 *
 *  OUTPUT PARAMETERS
 *  \param landcover_change_transfer        struct containing the following pft-specific public members:
 *   - transfer_litter_leaf         
 *   - transfer_litter_sap      
 *   - transfer_litter_heart      
 *   - transfer_litter_root    
 *   - transfer_litter_repr      
 *   - transfer_nmass_litter_leaf
 *   - transfer_nmass_litter_root  
 *   - transfer_nmass_litter_sap  
 *   - transfer_nmass_litter_heart  
 *   - transfer_harvested_products_slow  
 *   - transfer_harvested_products_slow_nmass 
 *											,the following patch-level public members:
 *   - transfer_cpool_fast  
 *   - transfer_cpool_slow
 *   - transfer_nmass_avail  
 *   - transfer_wcont_evap  
 *   - transfer_snowpack  
 *   - transfer_decomp_litter_mean  
 *   - transfer_k_soilfast_mean  
 *   - transfer_acflux_harvest  
 *   - transfer_anflux_harvest 
 *											,the following water soil layer-specific public member:
 *   - transfer_wcont 
 *											and the following century soil pool-specific public members:
 *   - transfer_sompool.cmass 
 *   - transfer_sompool.fireresist 
 *   - transfer_sompool.fracremain 
 *   - transfer_sompool.ligcfrac 
 *   - transfer_sompool.nmass 
 *   - transfer_sompool.ntoc 
 */
void donor_stand_change (Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], double cropstand_change[NCROPSTANDS_MAX], double& receiving_fraction, 
					int nnaturalstands, landcover_change_transfer& to) {

	gridcell.firstobj();
	while (gridcell.isobj) {
		double scale;

		Stand& stand = gridcell.getobj();

		if(stand.landcover != CROPLAND && stand.landcover != NATURAL && landcoverfrac_change[stand.landcover] < 0.0						
			|| stand.landcover == NATURAL && landcoverfrac_change[NATURAL] < 0.0 && (nnaturalstands == 1 || stand.natural_frac_change < 0.0)
			|| stand.landcover == CROPLAND && cropstand_change[stand.cftid] < 0.0) {

			// All landcovers that only have one stand:
			if(stand.landcover != CROPLAND && stand.landcover != NATURAL) {
				scale = -landcoverfrac_change[stand.landcover] / receiving_fraction / (double)stand.nobj;			
				stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
			}
			// Landcovers that may have several stands:
			else if(stand.landcover == NATURAL) {
				if(nnaturalstands > 1) {
					scale = -stand.natural_frac_change / receiving_fraction / (double)stand.nobj;
					//stand.frac already set in reduce_natural_stands() for these new stands
				}
				else {
					scale = -landcoverfrac_change[stand.landcover] / receiving_fraction / (double)stand.nobj;
					stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
				}
			}
			else if(stand.landcover == CROPLAND) {
				scale = -cropstand_change[stand.cftid] / receiving_fraction / (double)stand.nobj;
				stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid] * gridcell.landcoverfrac[CROPLAND]);
			}	
	
			stand.firstobj();
			while(stand.isobj) {
				Patch& patch = stand.getobj();

				// sum original litter C & N:
				for(int n=0; n<npft; n++)
				{
					to.transfer_litter_leaf[n] += patch.pft[n].litter_leaf * scale;
					to.transfer_litter_root[n] += patch.pft[n].litter_root * scale;
					to.transfer_litter_sap[n] += patch.pft[n].litter_sap * scale;
					to.transfer_litter_heart[n] += patch.pft[n].litter_heart * scale;
					to.transfer_litter_repr[n] += patch.pft[n].litter_repr * scale;

					to.transfer_nmass_litter_leaf[n] += patch.pft[n].nmass_litter_leaf * scale;
					to.transfer_nmass_litter_root[n] += patch.pft[n].nmass_litter_root*scale;
					to.transfer_nmass_litter_sap[n] += patch.pft[n].nmass_litter_sap * scale;
					to.transfer_nmass_litter_heart[n] += patch.pft[n].nmass_litter_heart * scale;

					if(ifslowharvestpool) {
						to.transfer_harvested_products_slow[n] += patch.pft[n].harvested_products_slow * scale;
						to.transfer_harvested_products_slow_nmass[n] += patch.pft[n].harvested_products_slow_nmass * scale;
					}
				}

				Vegetation& vegetation = patch.vegetation;
				vegetation.firstobj();
				while(vegetation.isobj) {

					Harvest_CN cp;

					Individual& indiv = vegetation.getobj();
					Patchpft& patchpft = patch.pft[indiv.pft.id];

					if(indiv.has_daily_turnover())
						cp.copy_from_indiv(indiv, true, false);
					else
						cp.copy_from_indiv(indiv, false, false);

					// Harvest of transferred areas:
					switch (indiv.pft.landcover)
					{
					case CROPLAND:
						harvest_crop(cp, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass);
						break;
					case PASTURE:
						harvest_pasture(cp, indiv.pft, indiv.alive);
						break;
					case NATURAL:
					case FOREST:
						harvest_wood(cp, indiv.pft, indiv.alive, 1.0);
						break;
					default:
						fail("Modify code to deal with landcover harvest at landcover change!\n");
					}

					turnover(indiv.pft.turnover_leaf, indiv.pft.turnover_root,
						indiv.pft.turnover_sap, indiv.pft.lifeform, indiv.pft.landcover,
						cp.cmass_leaf, cp.cmass_root, cp.cmass_sap, cp.cmass_heart,
						cp.nmass_leaf, cp.nmass_root, cp.nmass_sap, cp.nmass_heart,
						cp.litter_leaf,
						cp.litter_root,
						cp.nmass_litter_leaf,
						cp.nmass_litter_root,
						cp.nstore_longterm, 
						indiv.alive);

					gridcell.LC_updated = true;

					// In case any vegetation left (eg. cmass_root in pasture or grass in woodland):
					kill_remaining_vegetation(cp, indiv.pft, indiv.alive, indiv.istruecrop_or_intercropgrass(), false);

					//Sum added litter C & N:
					to.transfer_litter_leaf[indiv.pft.id] += cp.litter_leaf * scale;
					to.transfer_litter_root[indiv.pft.id] += cp.litter_root * scale;
					to.transfer_litter_sap[indiv.pft.id] += cp.litter_sap * scale;
					to.transfer_litter_heart[indiv.pft.id] += cp.litter_heart * scale;

					to.transfer_nmass_litter_leaf[indiv.pft.id] += cp.nmass_litter_leaf * scale;
					to.transfer_nmass_litter_root[indiv.pft.id] += cp.nmass_litter_root * scale;
					to.transfer_nmass_litter_sap[indiv.pft.id] += cp.nmass_litter_sap * scale;
					to.transfer_nmass_litter_heart[indiv.pft.id] += cp.nmass_litter_heart * scale;

					to.transfer_anflux_harvest += cp.anflux_harvest * scale;


					double change_frac;
					if(stand.landcover == NATURAL) {
						if(nnaturalstands > 1)
							change_frac = stand.natural_frac_change;
						else
							change_frac = landcoverfrac_change[NATURAL];
					}
					else if(stand.landcover == CROPLAND)
						change_frac = cropstand_change[stand.cftid];
					else
						change_frac = landcoverfrac_change[stand.landcover];

					gridcell.acflux_landuse_change += -cp.acflux_harvest * change_frac / (double)stand.nobj;
					gridcell.acflux_landuse_change_lc[stand.landcover] += -cp.acflux_harvest * change_frac / (double)stand.nobj;

					if(ifslowharvestpool) {
						to.transfer_harvested_products_slow[indiv.pft.id] += cp.harvested_products_slow * scale;
						to.transfer_harvested_products_slow_nmass[indiv.pft.id] += cp.harvested_products_slow_nmass * scale;
					}

					vegetation.nextobj();
				}

				// sum litter C & N:
				to.transfer_cpool_fast += patch.soil.cpool_fast * scale;
				to.transfer_cpool_slow += patch.soil.cpool_slow * scale;

				// sum soil C & N:
				for(int i=0; i<NSOMPOOL; i++) {
					to.transfer_sompool[i].cmass += patch.soil.sompool[i].cmass * scale;
					to.transfer_sompool[i].fireresist += patch.soil.sompool[i].fireresist * scale;
					to.transfer_sompool[i].fracremain += patch.soil.sompool[i].fracremain * scale;
					to.transfer_sompool[i].ligcfrac += patch.soil.sompool[i].ligcfrac * scale;
					to.transfer_sompool[i].litterme += patch.soil.sompool[i].litterme * scale;
					to.transfer_sompool[i].nmass += patch.soil.sompool[i].nmass * scale;
					to.transfer_sompool[i].ntoc += patch.soil.sompool[i].ntoc * scale;
				}

				to.transfer_nmass_avail += patch.soil.nmass_avail * scale;

				// sum wcont:
				for(int i=0; i<NSOILLAYER; i++) {
					to.transfer_wcont[i] += patch.soil.wcont[i] * scale;
				}
				to.transfer_wcont_evap += patch.soil.wcont_evap * scale;

				to.transfer_snowpack += patch.soil.snowpack * scale;
				to.transfer_snowpack_nmass += patch.soil.snowpack_nmass * scale;

				to.transfer_decomp_litter_mean += patch.soil.decomp_litter_mean * scale;
				to.transfer_k_soilfast_mean += patch.soil.k_soilfast_mean * scale;
				to.transfer_k_soilslow_mean += patch.soil.k_soilslow_mean * scale;

				stand.nextobj();
			}
		}
		gridcell.nextobj();
	}
}


/// Creates and kills stands at landcover change.
/** Harvest of reduced stands need to be done before with donor_stand_change().
 *  Should be followed by a call to receiving_stand_change() for transfer of 
 *  carbon, nitrogen and water.
 *  Do not call from loop with call to gridcell.nextobj.
 * 
 *  INPUT PARAMETERS
 *
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param changeLC							sum of all stands' absolute changes
 *  \param change_crop						sum of all crop stands' absolute changes		
 *  \param nnaturalstands					number of natural stands in the gridcell
 */
void stand_dynamics(Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], double changeLC, double change_crop, int nnaturalstands) {

	// dynamics for stands other than cropland (from updated landcoverfrac):
	if(!lcfrac_fixed && changeLC > 0.0) {
		for(int i=0; i<NLANDCOVERTYPES; i++) {
			if(i != CROPLAND) {
				if(run[i]) {
					// stand created
					if(gridcell.landcoverfrac_old[i] == 0.0 && gridcell.landcoverfrac[i] > 0.0) {
						gridcell.create_stand_lu((landcovertype)i, gridcell.landcoverfrac[i]);
					}
					// stand killed
					else if(gridcell.landcoverfrac_old[i] > 0.0 && gridcell.landcoverfrac[i] == 0.0) {
						gridcell.firstobj();
						while (gridcell.isobj) {
							Stand& stand = gridcell.getobj();
							if(stand.landcover == i) {
								gridcell.killobj();
							}
							else
								gridcell.nextobj();
						}
					}
					// new NATURAL stand created from other landcover type
					else if(i == NATURAL && landcoverfrac_change[i] > 0.0) {
						gridcell.create_stand_lu((landcovertype)i, landcoverfrac_change[i]);
					}
					// secondary natural stand killed if all of its area converted to other landcover type
					else if(i == NATURAL && nnaturalstands > 1 && landcoverfrac_change[i] < 0.0) {
						gridcell.firstobj();
						while (gridcell.isobj) {
							Stand& stand = gridcell.getobj();
							if(stand.landcover == i && stand.get_gridcell_fraction() == 0) {
								gridcell.killobj();
							}
							else
								gridcell.nextobj();
						}
					}
				}
			}
		}
	}

	// crop stand dynamics (from updated cftfrac):
	if(run[CROPLAND] && (change_crop>0.0 || landcoverfrac_change[CROPLAND]!=0.0)) {
		if(gridcell.landcoverfrac[CROPLAND]>0.0) {

			pftlist.firstobj();
			while (pftlist.isobj) {
				Pft& pft=pftlist.getobj();
				if(pft.landcover == CROPLAND && pft.cftid >= 0) {

					// Is this PFT already present in a crop stand ?
					bool present = false;
					gridcell.firstobj();
					while (gridcell.isobj && !present) {
						Stand& stand = gridcell.getobj();
						if (stand.landcover == CROPLAND && stand.pftid == pft.id)
							present = true;
						else
							gridcell.nextobj();
					}

					// Should this crop PFT be present in the gridcell this year ? 
					if(gridcell.cftfrac[pft.cftid]>0.0) {
						// if so, and not already present, create new crop stand
						if(!present) {
							gridcell.create_stand_lu(CROPLAND, gridcell.cftfrac[pft.cftid] * gridcell.landcoverfrac[CROPLAND], pft.cftid);
						}
					}
					else {
						// if not, and is present, kill stand
						if(present) {
							Stand& stand=gridcell.getobj();
							gridcell.killobj();
						}
					}
				}
				pftlist.nextobj(); // ... on to next PFT
			}
		}
		else if(gridcell.landcoverfrac_old[CROPLAND]>0.0) {	//(if !(gridcell.landcoverfrac[CROPLAND]>0.0))	
			gridcell.firstobj();
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();
				if(stand.landcover==CROPLAND) {
					stand.firstobj();
					gridcell.killobj();
				}
				else
					gridcell.nextobj();
			}
		}
	}
}

 /// Transfers litter etc. of reduced stands to expanding stands at landcover change.
/** Updates stand.frac. of expanded stands.
 *  Transfers carbon, nitrogen and water of harvested area to expanded areas from a temporary struct.
 *  Do not call from loop with call to gridcell.nextobj.
 *
 *  INPUT PARAMETERS
 *
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param cropstand_change					array with this year's difference in area fractions of the different crop stands			
 *  \param landcover_change_transfer& from  struct containing the following pft-specific public members:
 *   - transfer_litter_leaf         
 *   - transfer_litter_sap      
 *   - transfer_litter_heart      
 *   - transfer_litter_root    
 *   - transfer_litter_repr      
 *   - transfer_nmass_litter_leaf
 *   - transfer_nmass_litter_root  
 *   - transfer_nmass_litter_sap  
 *   - transfer_nmass_litter_heart  
 *   - transfer_harvested_products_slow  
 *   - transfer_harvested_products_slow_nmass 
 *											,the following patch-level public members:
 *   - transfer_cpool_fast  
 *   - transfer_cpool_slow
 *   - transfer_nmass_avail  
 *   - transfer_cpool_slow  
 *   - transfer_wcont_evap  
 *   - transfer_snowpack  
 *   - transfer_decomp_litter_mean  
 *   - transfer_k_soilfast_mean  
 *   - transfer_acflux_harvest  
 *   - transfer_anflux_harvest 
 *											,the following water soil layer-specific public member:
 *   - transfer_wcont 
 *											and the following century soil pool-specific public members:
 *   - transfer_sompool.cmass 
 *   - transfer_sompool.fireresist 
 *   - transfer_sompool.fracremain 
 *   - transfer_sompool.ligcfrac 
 *   - transfer_sompool.nmass 
 *   - transfer_sompool.ntoc 
 *  \param LCchangeCtransfer				whether to transfer carbon, nitrogen and water of reduced stands to expanding stands
 */
void receiving_stand_change (Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], double cropstand_change[NCROPSTANDS_MAX], 
							 landcover_change_transfer& from, bool LCchangeCtransfer) {

	gridcell.firstobj();
	while (gridcell.isobj) {
		Stand& stand = gridcell.getobj();
		if(stand.landcover != CROPLAND && landcoverfrac_change[stand.landcover] > 0.0 || stand.landcover == CROPLAND && cropstand_change[stand.cftid] > 0.0)
		{
			double old_frac, added_frac, new_frac;

			// define 
			if(stand.landcover != CROPLAND && stand.landcover != NATURAL) {
				old_frac = gridcell.landcoverfrac_old[stand.landcover];
				added_frac = landcoverfrac_change[stand.landcover];
				new_frac = gridcell.landcoverfrac[stand.landcover];
				stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
			}
			else if(stand.landcover == NATURAL) {
				if(stand.first_year == date.year) {	// expanding natural area always results in a new stand
					old_frac = 0.0;
					added_frac = landcoverfrac_change[stand.landcover];
					new_frac = landcoverfrac_change[stand.landcover];
					// stand.frac already set for new natural stands in stand_dynamics()
				}
				else {
					gridcell.nextobj();
					continue;
				}					
			}
			else if(stand.landcover == CROPLAND) {
				old_frac = gridcell.landcoverfrac_old[CROPLAND] * gridcell.cftfrac_old[stand.cftid];
				added_frac = cropstand_change[stand.cftid];
				new_frac = gridcell.landcoverfrac[CROPLAND] * gridcell.cftfrac[stand.cftid];
				stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid] * gridcell.landcoverfrac[CROPLAND]);
			}	

			if(LCchangeCtransfer) {
				stand.firstobj();
				while(stand.isobj) {
					Patch& patch=stand.getobj();

					// add litter C & N:
					for (int i=0; i<npft; i++) {
						Patchpft& patchpft = patch.pft[i];

						patchpft.litter_leaf = (patchpft.litter_leaf * old_frac + from.transfer_litter_leaf[i] * added_frac) / new_frac;
						patchpft.litter_sap = (patchpft.litter_sap * old_frac + from.transfer_litter_sap[i] * added_frac) / new_frac;
						patchpft.litter_heart = (patchpft.litter_heart * old_frac + from.transfer_litter_heart[i] * added_frac) / new_frac;
						patchpft.litter_root = (patchpft.litter_root * old_frac + from.transfer_litter_root[i] * added_frac) / new_frac;
						patchpft.litter_repr = (patchpft.litter_repr * old_frac + from.transfer_litter_repr[i] * added_frac) / new_frac;

						patchpft.nmass_litter_leaf = (patchpft.nmass_litter_leaf * old_frac + from.transfer_nmass_litter_leaf[i] * added_frac) / new_frac;
						patchpft.nmass_litter_root = (patchpft.nmass_litter_root * old_frac + from.transfer_nmass_litter_root[i] * added_frac) / new_frac;
						patchpft.nmass_litter_sap = (patchpft.nmass_litter_sap * old_frac + from.transfer_nmass_litter_sap[i] * added_frac) / new_frac;
						patchpft.nmass_litter_heart = (patchpft.nmass_litter_heart * old_frac + from.transfer_nmass_litter_heart[i] * added_frac) / new_frac;

						if(ifslowharvestpool)
							patchpft.harvested_products_slow = (patchpft.harvested_products_slow * old_frac + from.transfer_harvested_products_slow[i] * added_frac) / new_frac;
							patchpft.harvested_products_slow_nmass = (patchpft.harvested_products_slow_nmass * old_frac + from.transfer_harvested_products_slow_nmass[i] * added_frac) / new_frac;
					}

					// add soil C & N:
					patch.soil.cpool_fast = (patch.soil.cpool_fast * old_frac + from.transfer_cpool_fast * added_frac) / new_frac;
					patch.soil.cpool_slow = (patch.soil.cpool_slow * old_frac + from.transfer_cpool_slow * added_frac) / new_frac;

					for(int i=0; i<NSOMPOOL; i++) {
						patch.soil.sompool[i].cmass = (patch.soil.sompool[i].cmass * old_frac + from.transfer_sompool[i].cmass * added_frac) / new_frac;
						patch.soil.sompool[i].fireresist = (patch.soil.sompool[i].fireresist * old_frac + from.transfer_sompool[i].fireresist * added_frac) / new_frac;
						patch.soil.sompool[i].fracremain = (patch.soil.sompool[i].fracremain * old_frac + from.transfer_sompool[i].fracremain * added_frac) / new_frac;
						patch.soil.sompool[i].ligcfrac = (patch.soil.sompool[i].ligcfrac *old_frac + from.transfer_sompool[i].ligcfrac * added_frac) / new_frac;
						patch.soil.sompool[i].litterme = (patch.soil.sompool[i].litterme * old_frac + from.transfer_sompool[i].litterme * added_frac) / new_frac;
						patch.soil.sompool[i].nmass = (patch.soil.sompool[i].nmass * old_frac + from.transfer_sompool[i].nmass * added_frac) / new_frac;
						patch.soil.sompool[i].ntoc = (patch.soil.sompool[i].ntoc * old_frac + from.transfer_sompool[i].ntoc * added_frac) / new_frac;
					}

					patch.soil.nmass_avail = (patch.soil.nmass_avail * old_frac + from.transfer_nmass_avail * added_frac) / new_frac;
					

					// add other soil stuff:
					for(int i=0; i<NSOILLAYER; i++) {
						patch.soil.wcont[i] = (patch.soil.wcont[i] * old_frac + from.transfer_wcont[i] * added_frac) / new_frac;
					}
					patch.soil.wcont_evap = (patch.soil.wcont_evap * old_frac + from.transfer_wcont_evap * added_frac) / new_frac;

					patch.soil.snowpack = (patch.soil.snowpack * old_frac + from.transfer_snowpack * added_frac) / new_frac;
					patch.soil.snowpack_nmass = (patch.soil.snowpack_nmass * old_frac + from.transfer_snowpack_nmass * added_frac) / new_frac;

					patch.soil.decomp_litter_mean = (patch.soil.decomp_litter_mean * old_frac + from.transfer_decomp_litter_mean * added_frac) / new_frac;
					patch.soil.k_soilfast_mean = (patch.soil.k_soilfast_mean * old_frac + from.transfer_k_soilfast_mean * added_frac) / new_frac;
					patch.soil.k_soilslow_mean = (patch.soil.k_soilslow_mean * old_frac + from.transfer_k_soilslow_mean * added_frac) / new_frac;

					// add fluxes:
//					patch.fluxes.report_flux(Fluxes::HARVESTC, from.transfer_acflux_harvest * added_frac / new_frac); // no harvest C here anymore, goes to gridcell.acflux_harvest instead
					patch.fluxes.report_flux(Fluxes::HARVESTN, from.transfer_anflux_harvest * added_frac / new_frac);
		
					// set scaling factor to be used in growth() for scaling vegetation C and N:
					stand.scale_LC_change = old_frac / new_frac;

					// save individual C and N content for use in scale_indiv()
					for(unsigned int i=0; i<patch.vegetation.nobj ;i++) {
						Individual& indiv = patch.vegetation[i];

						indiv.save_cmass_luc();
						indiv.save_nmass_luc();
					}

					stand.nextobj();
				}
			}
		}
		gridcell.nextobj();
	}
}

/// Updates all landcover and crop stand area fractions each year, possibly resulting in the creation and killing of stands.
/** Harvests transferred areas and transfers litter etc. of reduced stands to expanding stands.
 *  Transfers litter etc. of reduced stands to expanding stands at landcover change.
 *  Do not call from loop with call to gridcell.nextobj.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module) {

	double landcoverfrac_change[NLANDCOVERTYPES];
	double cropstand_change[NCROPSTANDS_MAX];
	double changeLC=0.0;
	double change_crop=0.0;
	double receiving_fraction=0.0;
	int nnaturalstands=0;
	bool LCchangeCtransfer=true;
	landcover_change_transfer transfer;

	memset(landcoverfrac_change,0,NLANDCOVERTYPES*sizeof(double));
	memset(cropstand_change,0,NCROPSTANDS_MAX*sizeof(double));

	gridcell.LC_updated=false;

	gridcell.firstobj();
	while (gridcell.isobj) {
		Stand& stand=gridcell.getobj();
		stand.scale_LC_change = 1.0;
		gridcell.nextobj();
	}

	// get new landcover and crop stand area fractions from input files
	if(!all_fracs_const) {
		// this call returns 0, causing this function to return, if no significant landcover changes this year, 
		// sets LCchangeCtransfer to 0 if unbalanced landcover changes (if some landcovers are inactivated), thus inactivating transfer of C and N
		if(!checkLCchange(gridcell, landcoverfrac_change, cropstand_change, changeLC, change_crop, receiving_fraction, LCchangeCtransfer, input_module))
			return;
	}
	else return;

	transfer.allocate();

	// check how many natural stands exist
	// if necessary, identify which natural stands to reduce in area
	reduce_natural_stands(gridcell, landcoverfrac_change, nnaturalstands);

	 // handle harvest and turnover of reduced stands at landcover change
	donor_stand_change(gridcell, landcoverfrac_change, cropstand_change, receiving_fraction, nnaturalstands, transfer);

	// create and kill stands at landcover change
	stand_dynamics(gridcell, landcoverfrac_change, changeLC, change_crop, nnaturalstands);

	// transfer litter etc. of reduced stands to expanding stands at landcover change
	receiving_stand_change(gridcell, landcoverfrac_change, cropstand_change, transfer, LCchangeCtransfer);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  End of Landcover stand dynamics and C&N-partitioning  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////  Sowing date algorithm. ////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Monitors whether temperature limits have been attained for this pft this day
/** Called from crop_sowing_gridcell() each day 
 */
void check_crop_temp_limits(Climate& climate, Gridcellpft& gridcellpft) {

	Pft& pft = gridcellpft.pft;

	// check if spring conditions are present this day:
	if(pft.ifsdspring && climate.temp > pft.tempspring && climate.dtemp_31[29] <= pft.tempspring) {	// NB. after updating dtemp_31 with today's value
		// TeWW,TeCo,TeSf,TeRa: 12,14,13,12 (NB 5,14,15,5 in Bondeau 2007);
		if (climate.lat >= 0.0 && date.day > 300)
			gridcellpft.last_springdate=date.day-365;
		else 
			gridcellpft.last_springdate = date.day;
		gridcellpft.springoccurred = true;
	}

	// check if autumn conditions are present this day:
	if(pft.ifsdautumn) {

		if(climate.temp<pft.tempautumn && climate.dtemp_31[29]>=pft.tempautumn && !gridcellpft.autumnoccurred) { //TeWW,TeRa: 12,17
			if (climate.lat >= 0.0 && date.day<100)		
				gridcellpft.first_autumndate = date.day + 365;
			else
				gridcellpft.first_autumndate = date.day;
			gridcellpft.autumnoccurred = true;
		}

		// check if vernilisation conditions are present this day:

		if(climate.temp < pft.trg && climate.dtemp_31[29] >= pft.trg && !gridcellpft.vernstartoccurred) //TeWW,TeRa: 12,12
			gridcellpft.vernstartoccurred = true;

		if(climate.temp>pft.trg && climate.dtemp_31[29] <= pft.trg) { //TeWW,TeRa: 12,12
			if (climate.lat >= 0.0 && date.day > 300)
				gridcellpft.last_verndate = date.day - 365;	
			else 
				gridcellpft.last_verndate = date.day;
			gridcellpft.vernendoccurred = true;
		}
	}
}

/// Updates 20-year mean of dates when temperature limits obtained, used for sowing date calculation
/** Called from crop_sowing_gridcell() once a year
 */
void calc_crop_dates_20y_mean(Climate& climate, Gridcellpft& gridcellpft) {

	int y,startyear;
	Pft& pft = gridcellpft.pft;

	/////////////////////////////////////////////////////////////////////////////////
	// Check if spring and frost conditions occurred during the past year.		   //	
	// If not, set this year's date to either sdate_default or climate.coldestday: //
	/////////////////////////////////////////////////////////////////////////////////

	// if no spring occured during last year 
	if (pft.ifsdspring && !gridcellpft.springoccurred) {	//TeWW,TeCo,TeSf,TeRa
	
		if(climate.temp <= pft.tempspring)
			gridcellpft.last_springdate = date.day;	
		else
			gridcellpft.last_springdate = climate.coldestday;
	}

	// if no autumn occured during last year
	if(pft.ifsdautumn && !gridcellpft.autumnoccurred) {		//TeWW,TeRa
	
		if(climate.maxtemp < pft.tempautumn || climate.maxtemp >= pft.tempautumn && climate.mtemp_min < pft.tempautumn) {

			gridcellpft.first_autumndate = date.day;		
		}
		else {										//too warm
		
			gridcellpft.first_autumndate = climate.coldestday;

			if(climate.lat >= 0.0 && gridcellpft.first_autumndate < 180)
				gridcellpft.first_autumndate += 365;	
		}
	}

	// if temperature did not pass above the vernalisation temperature last year
	if(pft.ifsdautumn && !gridcellpft.vernendoccurred)		//TeWW,TeRa
			gridcellpft.last_verndate=gridcellpft.last_springdate+60; // to avoid last_verndate20 to precede last_springdate20 (Bondeau used coldest day)
																	  // 60 days is the maximum number of vernalization days

	///////////////////////////////////////////////////////////////////////////////////////
	// Update spring and frost date 20-year arrays and calculate 20 years average means: //
	///////////////////////////////////////////////////////////////////////////////////////

	// 1) this year
	if(pft.ifsdspring)									// TeWW,TeCo,TeSf,TeRa
		gridcellpft.last_springdate20=gridcellpft.last_springdate;

	if(pft.ifsdautumn) {								// TeWW,TeRa
	
		gridcellpft.first_autumndate20 = gridcellpft.first_autumndate;
		if(date.year == 1 && climate.lat >= 0.0)		// No autumn first half of first year, set value to same as for second year
			gridcellpft.first_autumndate_20[19] = gridcellpft.first_autumndate;

		gridcellpft.last_verndate20 = gridcellpft.last_verndate;
	}

	// 2) starting year (1st of 20 or less)
	startyear = 20 - (int)min(19, date.year);

	// 3) past 20 years or less
	for (y=startyear; y<20; y++) {
		if(pft.ifsdspring) {										// TeWW,TeCo,TeSf,TeRa

			gridcellpft.last_springdate_20[y-1] = gridcellpft.last_springdate_20[y];
			gridcellpft.last_springdate20 += gridcellpft.last_springdate_20[y];
		}
		if (pft.ifsdautumn)	{										// TeWW,TeRa
		
			gridcellpft.first_autumndate_20[y-1] = gridcellpft.first_autumndate_20[y];
			gridcellpft.first_autumndate20 += gridcellpft.first_autumndate_20[y];

			gridcellpft.last_verndate_20[y-1] = gridcellpft.last_verndate_20[y];
			gridcellpft.last_verndate20 += gridcellpft.last_verndate_20[y];
		}
	}

	// 4) 20 years average means:
	if(pft.ifsdspring) {											// TeWW,TeCo,TeSf,TeRa
		gridcellpft.last_springdate20 /= (int)min(20,date.year + 1);
		if (gridcellpft.last_springdate20 < 0)
			gridcellpft.last_springdate20 += 365;	
		gridcellpft.last_springdate_20[19] = gridcellpft.last_springdate;
	}

	if(pft.ifsdautumn) {											// TeWW,TeRa	
		gridcellpft.first_autumndate20 /= (int)min(20,date.year + 1);
		if (gridcellpft.first_autumndate20 > 364)
			gridcellpft.first_autumndate20 -= 365;			
		gridcellpft.first_autumndate_20[19] = gridcellpft.first_autumndate;

		gridcellpft.last_verndate20 /= (int)min(20,date.year + 1);
		if (gridcellpft.last_verndate20 < 0)
			gridcellpft.last_verndate20 += 365;
		gridcellpft.last_verndate_20[19] = gridcellpft.last_verndate;
	}
}

/// Calculates sdatecalc_temp
/** Called from crop_sowing_gridcell() once a year
 */
void set_sdatecalc_temp(Climate& climate, Gridcellpft& gridcellpft)
{
	Pft& pft = gridcellpft.pft;

	if(pft.ifsdautumn) {							// TeWW,TeRa:
	
		// Use autumn sowing if first_autumndate20 is set (autumn conditions met during the past 20 years):
		if(!((gridcellpft.first_autumndate20 == climate.testday_temp || gridcellpft.first_autumndate20 == climate.coldestday) && 
				gridcellpft.first_autumndate % 365 == gridcellpft.first_autumndate20)) {

			gridcellpft.sdatecalc_temp = gridcellpft.first_autumndate20;
			gridcellpft.wintertype = true;
		}
		// if not, use spring sowing
		else {	// if(gridcellpft.first_autumndate20==climate.coldestday)
		
			if(!((gridcellpft.last_springdate20 == climate.testday_temp || gridcellpft.last_springdate20 == climate.coldestday) && 
					gridcellpft.last_springdate == gridcellpft.last_springdate20)) {

				gridcellpft.sdatecalc_temp = gridcellpft.last_springdate20;
				gridcellpft.wintertype = false;
			}
			else {	// If neither spring nor autumn occurred during the past 20 years.
			
				if(climate.maxtemp < pft.tempspring)	// Too cold to sow at all.
					gridcellpft.sdatecalc_temp = -1;
				else {									// Too warm; avoid warmest period.
				
					gridcellpft.sdatecalc_temp = climate.coldestday;
					gridcellpft.wintertype = true;					
				}
			}
		}

		// If autumn first_autumndate20 is earlier than hlimitdate, use last_springdate20 (winter is too long):
		if (climate.lat >= 0.0 && gridcellpft.sdatecalc_temp <= pft.hlimitdatenh && gridcellpft.sdatecalc_temp > 180 
			|| climate.lat < 0.0 && gridcellpft.sdatecalc_temp <= pft.hlimitdatesh) {

			gridcellpft.sdatecalc_temp = gridcellpft.last_springdate20;	// use last_springdate20 disregarding earlier choices	
			gridcellpft.wintertype = false;
		}

		// Forced sowing date read from input file.
		// Calculated value used if value for pft not found in file.
		if(forcesowingdates && pft.forcesowingdate && gridcellpft.sdate_force >= 0) {

			if((abs(gridcellpft.sdate_force - gridcellpft.first_autumndate20) <= abs(gridcellpft.sdate_force - gridcellpft.last_springdate20)))
				gridcellpft.wintertype = true;
			else
				gridcellpft.wintertype = false;
			gridcellpft.sdatecalc_temp = gridcellpft.sdate_force;
		}
	}
	else if(pft.ifsdspring)	{							// TeCo,TeSf
	
		if(!strncmp(pft.name,"TeCo", strlen("TeCo")))
			gridcellpft.sdatecalc_temp = (int)(60.0 / 85.0 * (gridcellpft.last_springdate20 - climate.adjustlat) + 29.5 + climate.adjustlat);
		else if(!strncmp(pft.name,"TeSf", strlen("TeSf")))
			gridcellpft.sdatecalc_temp = gridcellpft.last_springdate20;
#if defined NEWSOWINGDATE
		else
			gridcellpft.sdatecalc_temp = gridcellpft.last_springdate20;
#endif

#ifdef LOW_SOWING_TEMPERATURE_LIMIT
		//Same lower temperature limit for sowing as for winter crops above.
		if((gridcellpft.last_springdate20 == climate.testday_temp && 
			gridcellpft.last_springdate == gridcellpft.last_springdate20)) {
			// If spring has not occurred during the past 20 years.		
			if(climate.maxtemp < pft.tempspring)	// Too cold to sow.
				gridcellpft.sdatecalc_temp = -1;
		}
#endif
	}
#ifdef HIGH_SOWING_TEMPERATURE_LIMIT
	// Climatic limits for TeWW growth:	
	if(!strncmp(pft.name,"TeWW", strlen("TeWW")) && climate.mtemp_min20 > 15.0)
		gridcellpft.sdatecalc_temp = -1;
#endif

	gridcellpft.springoccurred = false;	
	gridcellpft.vernstartoccurred = false;	
	gridcellpft.vernendoccurred = false;	
	gridcellpft.autumnoccurred = false;
}

/// Sets sdatecalc_prec first day of the rain period if NEWSOWINGDATE is defined
/** Called from crop_sowing_gridcell() each day 
 */
void set_sdatecalc_prec(Climate& climate, Gridcellpft& gridcellpft)
{
	Pft& pft = gridcellpft.pft;

	if(pft.hydrology == IRRIGATED) {
		gridcellpft.sdatecalc_prec = gridcellpft.sdate_default;
		gridcellpft.precoccurred = true;
	}
	else {

		bool SOAsia;

		if(climate.lat > -15.0 && climate.lat < 20.0 && climate.lon > 90.0)
			SOAsia = true;
		else
			SOAsia = false;

		if(!gridcellpft.precoccurred && (SOAsia && climate.sprec_2[1] >= 110.0 && climate.sprec_2[0] < 110.0
										|| !SOAsia && climate.sprec_2[1] >= 40.0 && climate.sprec_2[0] < 40.0)) {

			gridcellpft.first_precdate = date.day;
			gridcellpft.sdatecalc_prec = gridcellpft.first_precdate;
			gridcellpft.precoccurred = true;
		}

		if(date.day == climate.testday_prec) {	// December 31(364) north, June 30(180) in the south; just resets precoccurred and fcalc_prec.
		
			gridcellpft.precoccurred = false;
			gridcellpft.sdatecalc_prec = -1;
		}
	}
}

/// Calculates the Julian start and end day of a month.
/** January 1st is set to 0.
 */
void monthdates(int& start, int& end, int month) {
	int months[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	start = 0;
	int m = 0;
	while (m < month){
		start += months[m];
		m++;
	}
	end = start + months[m] - 1;
}

/// Calculates sowing window for each crop pft
/** Called from crop_sowing_gridcell() once a year
 */
void calc_sowing_windows(Gridcell& gridcell)
{
	Climate& climate = gridcell.climate;
	seasonality_type seasonality = climate.seasonality;
	int sow_month = 0;


	// Find the wettest month
	if(climate.seasonality == SEASONALITY_PREC || climate.seasonality == SEASONALITY_PRECTEMP) {

		double max = 0.0;
		double sum = 0.0;

		for (int m=0; m<12; m++) {
			sum = 0.0;
			for (int i=0; i<4; i++) {
				int mm = m + i;
				if(mm >= 12)
					mm -= 12;

				// Implement a check later to see if it makes any difference to to use the precipitation only
				if(true) {
					if (gridcell.climate.mpet20[mm] > 0.0) 
						sum += gridcell.climate.mprec20[mm] / gridcell.climate.mpet20[mm];
				}
				else
					sum += gridcell.climate.mprec20[mm];
			}

			if (sum > max) {
				max = sum;
				// Months are stored as, 0-11
				sow_month=m;
			}
		}
	}

	// Set gridcell-level sowing windows for crop pft:s
	pftlist.firstobj();
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if(pft.phenology == CROPGREEN) {

			bool temp_sdate = false, prec_sdate = false, def_sdate = false;

			// Different sowing date options for irrigated crops at sites with climate.seasonality == SEASONALITY_PRECTEMP:
			// 1. use temperature-dependent sowing limits (define IRRIGATED_USE_TEMP_SDATE)
			// 2. use precipitation-triggered sowing (IRRIGATED_USE_TEMP_SDATE undefined)

			// Determine, based upon site climate seasonality and sowing preferences for irrigated crops, if sowing should be triggered by 
			// temperature or precipitation, or whether to use a default sowing date.
#if defined IRRIGATED_USE_TEMP_SDATE
			if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC || seasonality == SEASONALITY_PRECTEMP && pft.hydrology == IRRIGATED)
				temp_sdate = true;
			else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP && pft.hydrology != IRRIGATED) && climate.prec_range != WET)
				prec_sdate = true;
#else
			if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC)
				temp_sdate = true;
			else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range != WET)
				prec_sdate = true;
#endif
			else // if(seasonality == SEASONALITY_NO) || (seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range == WET)
				def_sdate = true;


			if(temp_sdate && gridcellpft.sdatecalc_temp != -1) {

				// Set sowing window around sdatecalc_temp
				gridcellpft.swindow[0] = stepfromdate(gridcellpft.sdatecalc_temp, -15);
				gridcellpft.swindow[1] = stepfromdate(gridcellpft.sdatecalc_temp, 15);

				if(!gridcellpft.wintertype && dayinperiod(gridcellpft.swindow[0], stepfromdate(climate.coldestday, -100), climate.coldestday)) {

					gridcellpft.swindow[0] = climate.coldestday;
//					gridcellpft.swindow[0] = gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(gridcellpft.swindow[1], stepfromdate(climate.coldestday, -100), climate.coldestday))
						gridcellpft.swindow[1] = climate.coldestday;
				}

				if(gridcellpft.wintertype && dayinperiod(gridcellpft.swindow[1], climate.coldestday, stepfromdate(climate.coldestday, 100))) {

					gridcellpft.swindow[1] = climate.coldestday;
//					gridcellpft.swindow[1] = gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(gridcellpft.swindow[0], climate.coldestday, stepfromdate(climate.coldestday, 100)))
						gridcellpft.swindow[0] = climate.coldestday;
				}
			}

			if(prec_sdate) {

					monthdates(gridcellpft.swindow[0], gridcellpft.swindow[1],sow_month);

					// A conservative choice to expand the month  that is used for the search
					gridcellpft.swindow[0] = stepfromdate(gridcellpft.swindow[0], -15);
			}
			else if(def_sdate) {
				gridcellpft.swindow[0] = gridcellpft.sdate_default;
				gridcellpft.swindow[1] = stepfromdate(gridcellpft.sdate_default, 15);
			}

			// Includes all temperature limits for sowing set in set_sdatecalc_temp() also for sites with any type of temperature seasonality
			if(gridcellpft.sdatecalc_temp == -1) {
				gridcellpft.swindow[0] = -1;
				gridcellpft.swindow[1] = -1;
			}

		}
		pftlist.nextobj();
	}
}

/// Updates various climate 20-year means, used for sowing date calculation
/** Called from crop_sowing_gridcell() once a year
 */
void calc_m_climate_20y_mean(Climate& climate)
{
	int m, y;
	int startyear = 20 - (int)min(19, date.year);
	double var_temp = 0, var_prec = 0;
	double mtemp20kelvin[12], prec_pet_ratio20[12];
	double mprec_petmin_thisyear = 1.0;
	double mprec_petmax_thisyear = 0.0;

	memset(mtemp20kelvin, 0, 12 * sizeof(double));
	memset(prec_pet_ratio20, 0, 12 * sizeof(double));

	for(m=0; m<12; m++) {

		// 1) this year
		climate.mtemp20[m] = climate.mtemp_year[m];
		climate.mprec20[m] = climate.mprec_year[m];
		climate.mpet20[m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
			climate.mprec_pet20[m] = climate.mprec_year[m] / climate.mpet_year[m];
		else
			climate.mprec_pet20[m] = 0.0;


		if(climate.mprec_year[m] / climate.mpet_year[m] < mprec_petmin_thisyear)
			mprec_petmin_thisyear = climate.mprec_year[m] / climate.mpet_year[m];
		if(climate.mprec_year[m] / climate.mpet_year[m] > mprec_petmax_thisyear)
			mprec_petmax_thisyear = climate.mprec_year[m] / climate.mpet_year[m];

		// 2) past 20 years or less
		for (y=startyear; y<20; y++) {
			climate.mtemp_20[y-1][m] = climate.mtemp_20[y][m];
			climate.mtemp20[m] += climate.mtemp_20[y][m];

			climate.mprec_20[y-1][m] = climate.mprec_20[y][m];
			climate.mprec20[m] += climate.mprec_20[y][m];

			climate.mpet_20[y-1][m] = climate.mpet_20[y][m];
			climate.mpet20[m] += climate.mpet_20[y][m];

			climate.mprec_pet_20[y-1][m] = climate.mprec_pet_20[y][m];
			climate.mprec_pet20[m] += climate.mprec_pet_20[y][m];
		}
		// 3) 20 years average means:
		climate.mtemp20[m] /= min(20, date.year + 1);
		climate.mprec20[m] /= min(20, date.year + 1);
		climate.mpet20[m] /= min(20, date.year + 1);
		climate.mprec_pet20[m] /= min(20, date.year + 1);

		climate.mtemp_20[19][m] = climate.mtemp_year[m];
		climate.mprec_20[19][m] = climate.mprec_year[m];
		climate.mpet_20[19][m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
			climate.mprec_pet_20[19][m] = climate.mprec_year[m] / climate.mpet_year[m];
		else
			climate.mprec_pet_20[19][m] = 0.0;
	}

	climate.mprec_petmin20 = mprec_petmin_thisyear;
	climate.mprec_petmax20 = mprec_petmax_thisyear;
	for (y=startyear; y<20; y++) {
		climate.mprec_petmin_20[y-1] = climate.mprec_petmin_20[y];
		climate.mprec_petmin20 += climate.mprec_petmin_20[y];
		climate.mprec_petmax_20[y-1] = climate.mprec_petmax_20[y];
		climate.mprec_petmax20 += climate.mprec_petmax_20[y];
	}
	climate.mprec_petmin20 /= min(20, date.year + 1);
	climate.mprec_petmin_20[19] = mprec_petmin_thisyear;
	climate.mprec_petmax20 /= min(20, date.year + 1);
	climate.mprec_petmax_20[19] = mprec_petmax_thisyear;
}

/// Determines climate seasonality of gridcell
/** Called from crop_sowing_gridcell() once a year
 */
void calc_seasonality(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;
	double var_temp = 0, var_prec = 0;
	double TEMPMIN = 10.0; // temperature limit of coldest month used to determine type of temperature seasonality
	const int NMONTH = 12;												
	double mtempKelvin[NMONTH], prec_pet_ratio20[12];
	double maxprec_pet20 = 0.0;
	double minprec_pet20 = 1000;

	memset(mtempKelvin, 0, NMONTH * sizeof(double));
	memset(prec_pet_ratio20, 0, NMONTH * sizeof(double));


	// calculate absolute temperature and prec/pet ratio for each month this year
	for(int i=0; i < NMONTH; ++i) {
		 // The temperature has got to be in Kelvin, the limit 0.010 is based on that.
		mtempKelvin[i] = gridcell.climate.mtemp20[i] + 273.15;
		// Calculate precipitation/PET ratio if monthly PET is above zero
		prec_pet_ratio20[i] = (gridcell.climate.mpet20[i] > 0) ? gridcell.climate.mprec20[i] / gridcell.climate.mpet20[i] : 0;
	}
	
	// calculate variation coeffecients of temperature and prec/pet ratio for this year
	var_temp = variation_coefficient(mtempKelvin, NMONTH);
	var_prec = variation_coefficient(prec_pet_ratio20, NMONTH);

	gridcell.climate.var_prec = var_prec;
	gridcell.climate.var_temp = var_temp;

	if (var_prec <= 0.4 && var_temp <= 0.010)				// no seasonality
		climate.seasonality_lastyear = SEASONALITY_NO;				// 0
	else if (var_prec > 0.4) {

		if(var_temp <= 0.010)								// precipitation seasonality only
			climate.seasonality_lastyear = SEASONALITY_PREC;			// 1
		else if(var_temp > 0.010) {
		
			if(gridcell.climate.mtemp_min20 > TEMPMIN)		// both seasonalities, but "weak" temperature seasonality (coldest month > 10degC)
				climate.seasonality_lastyear = SEASONALITY_PRECTEMP;	// 2
			else if(gridcell.climate.mtemp_min20 < TEMPMIN)	// both seasonalities, but temperature most important
				climate.seasonality_lastyear = SEASONALITY_TEMPPREC;	// 4
		}
	}
	else if(var_prec <= 0.4) {

		if (var_temp > 0.010)								// Temperature seasonality only
			climate.seasonality_lastyear = SEASONALITY_TEMP;			// 3
		 
/*															// SEASONALITY_TEMPWARM currently not used, default sdate value is coldest day anyway when always above PFT limit.
			if(gridcell.climate.mtemp_min20 < TEMPMIN)		// Temperature seasonality only
				climate.seasonality_lastyear = SEASONALITY_TEMP;		// 3
			else if(gridcell.climate.mtemp_min20 >= TEMPMIN))	// Temperature seasonality, always above 10 degrees
				climate.seasonality_lastyear = SEASONALITY_TEMPWARM;	// 5
*/
	}

	for(int m=0; m<12; m++) {
		if(climate.mprec_pet20[m] > maxprec_pet20)
			maxprec_pet20 = climate.mprec_pet20[m];
		if(climate.mprec_pet20[m] < minprec_pet20)
			minprec_pet20 = climate.mprec_pet20[m];
	}

	if(minprec_pet20 <= 0.5 && maxprec_pet20 <= 0.5)							//Extremes of monthly means
		climate.prec_seasonality_lastyear = DRY;						// 0
	else if(minprec_pet20 <= 0.5 && maxprec_pet20>0.5 && maxprec_pet20 <= 1.0)		
		climate.prec_seasonality_lastyear = DRY_INTERMEDIATE;		// 1
	else if(minprec_pet20 <= 0.5 && maxprec_pet20 > 1.0)
		climate.prec_seasonality_lastyear = DRY_WET;					// 2
	else if(minprec_pet20 > 0.5 && minprec_pet20 <= 1.0 && maxprec_pet20 > 0.5 && maxprec_pet20 <= 1.0)
		climate.prec_seasonality_lastyear = INTERMEDIATE;			// 3
	else if(minprec_pet20 > 1.0 && maxprec_pet20 > 1.0)
		climate.prec_seasonality_lastyear = WET;						// 5
	else if(minprec_pet20 > 0.5 && minprec_pet20 <= 1.0 && maxprec_pet20 > 1.0)		
		climate.prec_seasonality_lastyear = INTERMEDIATE_WET;		// 4
	else
		dprintf("Problem with calculating precipitation seasonality !\n");

	if(climate.mprec_petmin20 <= 0.5 && climate.mprec_petmax20 <= 0.5)			//Average of extremes
		climate.prec_range_lastyear = DRY;							//0
	else if(climate.mprec_petmin20 <= 0.5 && climate.mprec_petmax20 > 0.5 && climate.mprec_petmax20 <= 1.0)
		climate.prec_range_lastyear = DRY_INTERMEDIATE;				//1
	else if(climate.mprec_petmin20 <= 0.5 && climate.mprec_petmax20 > 1.0)			
		climate.prec_range_lastyear = DRY_WET;						//2
	else if(climate.mprec_petmin20 > 0.5 && climate.mprec_petmin20 <= 1.0 && climate.mprec_petmax20 > 0.5 && climate.mprec_petmax20 <= 1.0)
		climate.prec_range_lastyear = INTERMEDIATE;					//3
	else if(climate.mprec_petmin20 > 1.0 && climate.mprec_petmax20 > 1.0)
		climate.prec_range_lastyear = WET;							//5
	else if(climate.mprec_petmin20 > 0.5 && climate.mprec_petmin20 <= 1.0 && climate.mprec_petmax20 > 1.0)
		climate.prec_range_lastyear = INTERMEDIATE_WET;				//4
	else
		dprintf("Problem with calculating precipitation range !\n");

	if(climate.mtemp_max20 <= 10)
		climate.temp_seasonality_lastyear = COLD;					//0
	else if(climate.mtemp_min20 <= 10 && climate.mtemp_max20 > 10 && climate.mtemp_max20 <= 30)
		climate.temp_seasonality_lastyear = COLD_WARM;				//1
	else if(climate.mtemp_min20 <= 10 && climate.mtemp_max20 > 30)
		climate.temp_seasonality_lastyear = COLD_HOT;					//2
	else if(climate.mtemp_min20 > 10 && climate.mtemp_max20 <= 30)
		climate.temp_seasonality_lastyear = WARM;					//3
	else if(climate.mtemp_min20 > 30)
		climate.temp_seasonality_lastyear = HOT;						//5
	else if(climate.mtemp_min20 > 10 && climate.mtemp_max20 > 30)	
		climate.temp_seasonality_lastyear = WARM_HOT;				//4
	else
		dprintf("Problem with calculating temperature seasonality !\n");
}

void update_seasonality(Climate& climate) {

	climate.seasonality = climate.seasonality_lastyear;
	climate.temp_seasonality = climate.temp_seasonality_lastyear;
	climate.prec_seasonality = climate.prec_seasonality_lastyear;
	climate.prec_range = climate.prec_range_lastyear;
}

/// Monitors climate history relevant for sowing date calculation. Calculates initial sowing dates/windows.
void crop_sowing_gridcell(Gridcell& gridcell) {

	int d;
	Climate& climate = gridcell.climate;

	if (date.year==0 && date.day == 0) {

		for (d=0; d<10; d++)
			climate.dprec_10[d] = climate.prec;
		for (d=0; d<2; d++)
			climate.sprec_2[d] = climate.prec;
	}

	climate.sprec_2[0] = climate.sprec_2[1];
	climate.sprec_2[1] = climate.prec;
	for (d=0; d<9; d++) {
		climate.dprec_10[d] = climate.dprec_10[d+1];
		climate.sprec_2[1] += climate.dprec_10[d];
	}
	climate.dprec_10[9] = climate.prec;

	if(climate.temp > climate.maxtemp)	//To know if temperature rises over vernalization limit
		climate.maxtemp = climate.temp;

	// Loop through PFTs
	pftlist.firstobj();
	while (pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if(pft.landcover == CROPLAND) {

			// If NEWSOWINGDATE is defined, all pft:s enter this code
			if(pft.ifsdcalc) {		// TeWW,TrRi,TeCo,TrMi,TrMa,TeSf,TrPe,TeRa; sdate set in getgridcell() kept for the rest;  (no code here for TrRi)
			
				if(pft.ifsdtemp) {	// TeWW,TeCo,TeSf,TeRa
	
					// Check whether temperature limits have been attained today
					check_crop_temp_limits(climate, gridcellpft);

					if(date.day == climate.testday_temp) {	// June 30(180) in the north, December 31(364) in the south
					
						// Update 20-year mean of dates when temperature limits obtained
						calc_crop_dates_20y_mean(climate, gridcellpft);

						// Determine sdatecalc_temp:						
						set_sdatecalc_temp(climate, gridcellpft);

						// Reset maxtemp to today's value
						climate.maxtemp = climate.temp;
					}
				}
#ifndef NEWSOWINGDATE
				// Determine sdatecalc_prec:
				if(pft.ifsdprec)	//TeCo,TrMi,TrMa,TrPe:
					set_sdatecalc_prec(climate, gridcellpft);
#endif
			}
		}

		pftlist.nextobj();
	}

	if(date.islastmonth && date.islastday) {
		// Update various climate 20-year means
		calc_m_climate_20y_mean(climate);
		// Determines climate seasonality of gridcell
		calc_seasonality(gridcell);
	}

	if(date.day == climate.testday_temp) {	// day 180/364

		update_seasonality(climate);

		// Calculate sowing window for each crop pft
		calc_sowing_windows(gridcell);
	}
}

/// old temperature-dependent sowing date method (Bondeau et al. 2007)
void Crop_sowing_date_temp(Patch& patch, Pft& pft) {

	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];

	patchpft.set_cropphen()->sdate = gridcellpft.sdatecalc_temp;
}

/// old precipitation-dependent sowing date method (Bondeau et al. 2007)
void Crop_sowing_date_prec(Patch& patch, Pft& pft) {

	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];
	int first_sowdate;
	int last_sowdate;

	// Limits here (dates, latitudes) derive from Bondeau code.

	if(!strncmp(pft.name,"TeCo", strlen("TeCo"))) {	// Sådd mellan sdatecalc_temp och sdate_default (140)
	
		first_sowdate = gridcellpft.sdatecalc_temp;
		last_sowdate = gridcellpft.sdate_default;

#if defined ALLOW_TWOSEASONSPREC
		// allows for two growing seasons; may sow outside the window ! eg. sow 140 - harvest 360 - sow 361
		if(climate.lat < 45.0) {	// Precipitation only determines sdate at latitudes < 45.0
		
			if(date.day == ppftcrop.sdate) {

				if(gridcellpft.precoccurred)
					return;							// use gridcellpft.sdatecalc_temp as sdate
				else
					patchpft.cropphen->sdate = -1;	// wait for rain period to begin
			}
			else if(dayinperiod(date.day,first_sowdate,last_sowdate) && gridcellpft.precoccurred)	// sow when rain period begins
				ppftcrop.sdate = date.day;
			else if(date.day == last_sowdate && !gridcellpft.precoccurred)							// if rainperiod has not begun at last_sowdate, sow anyway
				ppftcrop.sdate = date.day;
		}
		else
			return;									//use gridcellpft.sdatecalc_temp as sdate
#else
		// Sowing window starts here at sdatecalc_temp:
		if (date.day == ppftcrop.sdate && !gridcellpft.precoccurred && climate.lat < 45.0)		// sdatecalc_temp: if rain period not yet started, wait for rain
			ppftcrop.sdate = -1;																
		// Sowing occurred when rain has triggered sdatecalc_prec to be set:
		if (date.day == gridcellpft.sdatecalc_prec && ppftcrop.sdate == -1 && climate.lat < 45.0) {	// sdatecalc_temp not set (sdate set to -1 at hdate)
		
			if(dayinperiod(date.day,first_sowdate,last_sowdate))
				ppftcrop.sdate = date.day;
		}
		// If sdatecalc_prec has not been set by rain at sdate_default, sow anyway:
		if(date.day == gridcellpft.sdate_default && !gridcellpft.precoccurred && climate.lat < 45.0)	// if no rain at occurred before day 140, sow anyway 
			ppftcrop.sdate = date.day;
#endif

	}
	else {			//TrMi,TrMa,TrPe:					// sowing has to occur between firstsowdatenh/sh and day 210/30
	
		if(climate.lat >= 0.0) {
			first_sowdate = pft.firstsowdatenh_prec;
			last_sowdate = 210;
		}
		else {
			first_sowdate = pft.firstsowdatesh_prec;
			last_sowdate = 30;
		}

#if defined ALLOW_TWOSEASONSPREC
		// allows for two growing seasons
		if (date.day == first_sowdate) {	
			// if rain period has started before first_sowdate, sow immediately
			if(gridcellpft.precoccurred)
				ppftcrop.sdate = first_sowdate;
		}
		else if(dayinperiod(date.day,first_sowdate,last_sowdate) && gridcellpft.precoccurred)	// two seasons if hdate < last_sowdate
			ppftcrop.sdate = date.day;
		else if(date.day == last_sowdate && !gridcellpft.precoccurred)
			ppftcrop.sdate = date.day;
#else
		if (date.day == gridcellpft.sdatecalc_prec)	{ //firstsowdatenh_prec = sdatenh-(20 to 40)
		
			if(date.day <= pft.firstsowdatenh_prec && climate.lat >= 0.0 || date.day <= pft.firstsowdatesh_prec && date.day > 180 && climate.lat < 0.0)
				//if calculated sowing date is earlier than 20-40 days before the default sdate, use the latter
				ppftcrop.sdate=first_sowdate;
			else if (date.day <= 210 && climate.lat >= 0.0 || (date.day <= 30 || date.day > 180) && climate.lat < 0.0)
				ppftcrop.sdate = gridcellpft.sdatecalc_prec;
		}

		if (gridcellpft.sdatecalc_prec == -1 && (date.day == 210 && climate.lat >= 0.0 || date.day == 30 && climate.lat < 0.0))
			ppftcrop.sdate = date.day;
#endif
	}
}

/// sowing date calculatation with two growing seasons
void Crop_sowing_date_rice(Patch& patch, Pft& pft) {
	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];

	if(!gridcellpft.singlecrop)	{ // Rice with only one season has fixed sdate and hlimitdate.
	
		if(date.day == ppftcrop.sdate) {

			if(ppftcrop.maincrop) {
				if(climate.lat >= 0.0)
					ppftcrop.hlimitdate = pft.hlimitdatesh;
				else
					ppftcrop.hlimitdate = pft.hlimitdatenh;
			}
			else
				ppftcrop.hlimitdate = gridcellpft.hlimitdate_default;
		}
		else if(date.day == stepfromdate(ppftcrop.hdate, 1) && ppftcrop.hdate != -1) {

			if(ppftcrop.maincrop) {
				ppftcrop.maincrop = false;
				ppftcrop.sdate = (ppftcrop.hdate + 30) % 365;	
			}
			else {
				ppftcrop.maincrop = true;
				ppftcrop.sdate = gridcellpft.sdate_default;
			}
		}
	}
}

/// Forced sowing date read from input file.
/** Calculated value used if value for pft not found in file.
 */
void Crop_sowing_date_forced(Patch& patch, Pft& pft) {
	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];

	if(pft.forcesowingdate && pft.phenology != ANY && !(!strncmp(pft.name,"TrRi", strlen("TrRi")) && !gridcellpft.singlecrop)) {
														//not applicable for rice with several growingseasons
		if(gridcellpft.sdate_force >= 0) {

			patchpft.cropphen->sdate = gridcellpft.sdate_force;

			// avoid hlimitdate to cut growing period too early
			if (climate.lat >= 0.0 && gridcellpft.sdate_force <= pft.hlimitdatenh && gridcellpft.sdate_force > 180 
					|| climate.lat < 0.0 && gridcellpft.sdate_force <= pft.hlimitdatesh) {

				if(gridcellpft.sdate_force >= 0.0)
					patchpft.cropphen->hlimitdate = gridcellpft.sdate_force - 1;
				else
					patchpft.cropphen->hlimitdate = 364;
			}
			else {
				// reset hlimitdate to default (in case changed by Crop_sowing_date_new())
				if (climate.lat >= 0.0)
					patchpft.cropphen->hlimitdate = pft.hlimitdatenh;
				else
					patchpft.cropphen->hlimitdate = pft.hlimitdatesh;
			}
		}
	}
	else if(pft.forcesowingdate && (!strncmp(pft.name,"TrRi", strlen("TrRi")) && !gridcellpft.singlecrop))
		dprintf("No sowing date data available for rice with several growingseasons, using calculated values\n");
}

/// Sowing date method from Bondeau et al. 2007
/** Enters here every day when growingseason==false
 */
void Crop_sowing_date(Patch& patch, Pft& pft) {
	Patchpft& patchpft=patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

	// for crop pft:s with calculated sowing dates; all others use pft default values
	if(pft.ifsdcalc) {

		if(pft.ifsdtemp) {	//TeWW,TeCo,TeSf,TeRa		
			if(date.day == stepfromdate(ppftcrop.hdate, 1) && ppftcrop.hdate != -1 || date.day == stepfromdate(ppftcrop.hlimitdate, 1) && ppftcrop.hlimitdate != -1)
				Crop_sowing_date_temp(patch, pft);
		}

		if(pft.ifsdprec)	//TeCo,TrMi,TrMa,TrPe:
			Crop_sowing_date_prec(patch, pft);

		if(!strncmp(pft.name,"TrRi", strlen("TrRi")))
			Crop_sowing_date_rice(patch, pft);
	}
}

/// Sowing date method from Waha et al. 2010
/** Enters here every day when growingseason==false
 */
void Crop_sowing_date_new(Patch& patch, Pft& pft) {

	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Gridcell& gridcell = patch.stand.gridcell;
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];
	Climate& climate = gridcell.climate;
	seasonality_type seasonality = climate.seasonality;
	int length_growseas_def;
	bool temp_sdate = false, prec_sdate = false, def_sdate = false;

	// Different sowing date options for irrigated crops at sites with climate.seasonality == SEASONALITY_PRECTEMP:
	// 1. use temperature-dependent sowing limits (define IRRIGATED_USE_TEMP_SDATE)
	// 2. use precipitation-triggered sowing (IRRIGATED_USE_TEMP_SDATE undefined)
	// Option not to constrain sowing to the sowing date window (as in the old sowing date method); SD_TEMP_WINDOW undefied

	// Determine, based upon site climate seasonality and sowing preferences for irrigated crops, if sowing should be triggered by 
	// temperature or precipitation, or whether to use a default sowing date.
#if defined IRRIGATED_USE_TEMP_SDATE
	if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC || seasonality == SEASONALITY_PRECTEMP && pft.hydrology == IRRIGATED)
		temp_sdate = true;
	else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP && pft.hydrology != IRRIGATED) && climate.prec_range != WET)
		prec_sdate = true;
#else
	if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC)
		temp_sdate = true;
	else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range != WET)
		prec_sdate = true;
#endif
	else // if(seasonality == SEASONALITY_NO) || (seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range == WET)
		def_sdate = true;


	// adjust sowing window if too close to hlimitdate (before)
	if(temp_sdate && patchpft.swindow[0] != -1) {

		if(dayinperiod(patchpft.swindow[0], stepfromdate(ppftcrop.hlimitdate, -100), ppftcrop.hlimitdate)) {

			patchpft.swindow[0] = ppftcrop.hlimitdate + 1;
			if(dayinperiod(patchpft.swindow[1], stepfromdate(ppftcrop.hlimitdate, -100), ppftcrop.hlimitdate))
				patchpft.swindow[1] = ppftcrop.hlimitdate + 1;
		}
	}

	// option not to constrain sowing to the sowing date window (as in the old sowing date method)
#ifndef SD_TEMP_WINDOW
	if(temp_sdate) {
		if(date.day == stepfromdate(ppftcrop.hdate, 1) && ppftcrop.hdate != -1 || date.day == stepfromdate(ppftcrop.hlimitdate, 1) && ppftcrop.hlimitdate > -1)
			Crop_sowing_date_temp(patch, pft);
		return;
	}
#endif

	// monitor climate triggers within the sowing window
	if(dayinperiod(date.day, patchpft.swindow[0],patchpft.swindow[1])) {

		if(date.day != patchpft.swindow[1]) {

			if(temp_sdate) {
#if defined SD_TEMP_WINDOW
				if(gridcellpft.wintertype && climate.temp < pft.tempautumn || !gridcellpft.wintertype && climate.temp > pft.tempspring)
					ppftcrop.sdate = date.day;
#endif
			}
			else if(prec_sdate) {	
				if(climate.prec > 0.1 || pft.hydrology == IRRIGATED)
					ppftcrop.sdate = date.day;
			}
			else // if(def_sdate)
				ppftcrop.sdate = date.day; // first day of sowing window
		}
		else	 // last day of sowing window	
			ppftcrop.sdate = date.day;

		// calculation of hucountend (last day of heat unit sampling period); sdate is first day
		// NB. currently the sampling periods (which roughly correspond to growing periods) are unrealistically long,
		// Since shorter growing periods result in lower yield, a revision of this section will also make a revision of crop productivity necessary
		if(date.day == ppftcrop.sdate) {

			if(gridcellpft.sdate_default > gridcellpft.hlimitdate_default)
				length_growseas_def = gridcellpft.hlimitdate_default + 365 - gridcellpft.sdate_default;
			else
				length_growseas_def = gridcellpft.hlimitdate_default - gridcellpft.sdate_default;

			if(prec_sdate)
				ppftcrop.hlimitdate = stepfromdate(date.day, length_growseas_def);
			else
				ppftcrop.hlimitdate = gridcellpft.hlimitdate_default;

			length_growseas_def = min(length_growseas_def, 245);				// set an upper limit of 245 for the growing season

			if(pft.ifsdautumn) { // winter crops

				if(temp_sdate)
					ppftcrop.hucountend = stepfromdate(ppftcrop.hlimitdate, -20);
				else if(prec_sdate && climate.prec_seasonality <= DRY_WET) { // dry some time during the year
					if(pft.hydrology == IRRIGATED)
						ppftcrop.hucountend = stepfromdate(date.day, 230);
					else
						ppftcrop.hucountend = stepfromdate(date.day, 210);		 // shorter growing period when risk for water stress.
				}
				else if(def_sdate)
					ppftcrop.hucountend = stepfromdate(date.day, 230);
			}
			else if(!strncmp(pft.name,"TrRi", strlen("TrRi")) && gridcellpft.singlecrop) // rice
				ppftcrop.hucountend = stepfromdate(date.day, 230);
			else { // all other crops
				if(prec_sdate && climate.prec_seasonality <= DRY_WET) {	// dry some time during the year

					if(pft.hydrology == IRRIGATED)
						ppftcrop.hucountend = stepfromdate(date.day, length_growseas_def);
					else
						ppftcrop.hucountend = stepfromdate(date.day, min(length_growseas_def, 210)); // shorter growing period when risk for water stress.
				}
				else
					ppftcrop.hucountend = stepfromdate(date.day, length_growseas_def);
			}			
		}
	}

	if(!strncmp(pft.name,"TrRi", strlen("TrRi")))	// no changes are made if singlecrop == true
		// sets hlimitdate for maincrop and sdate (30 days after hdate) and hlimitdate (default) for second crop
		Crop_sowing_date_rice(patch, pft);
}

/// handles sowing date calculations for crop pft:s on patch level
void crop_sowing_patch(Patch& patch) {

	pftlist.firstobj();
	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;

	// Loop through PFTs
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if(patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

			if(date.day == climate.testday_temp) {

				patchpft.swindow[0] = gridcellpft.swindow[0];
				patchpft.swindow[1] = gridcellpft.swindow[1];
			}

			if(!ppftcrop.growingseason) {

				// copy sowing window from gridcellpft
				if(date.day == stepfromdate(ppftcrop.hdate, 1) && ppftcrop.hdate != -1 || date.day == climate.testday_temp) {

					if(gridcellpft.swindow[0] == -1) {
						gridcellpft.sowing_restriction = true;
						ppftcrop.hdate = -1;
						ppftcrop.eicdate = -1;	//redundant
					}
					else {
						gridcellpft.sowing_restriction = false;
					}
				}

				if(!gridcellpft.sowing_restriction)	{

#if defined NEWSOWINGDATE
 					// new sowing date method (Waha et al. 2010)
					Crop_sowing_date_new(patch, pft);
#else				// old sowing date method (Bondeau et al. 2007)
					Crop_sowing_date(patch, pft);
#endif
					// use sowing date read from input file if (pft.forcesowingdate==true)
					if(forcesowingdates)
						Crop_sowing_date_forced(patch, pft);
				}
			}

			// set eicdate (last intercrop day)
			if(ppftcrop.sdate != -1) {

				if(!ppftcrop.growingseason)
					ppftcrop.eicdate = stepfromdate(ppftcrop.sdate, - 15);

				if(dayinperiod(date.day, ppftcrop.eicdate, ppftcrop.sdate)) {

					if(ppftcrop.intercropseason)
						ppftcrop.eicdate = date.day;
					else if(ppftcrop.bicdate == ppftcrop.sdate) 
						ppftcrop.bicdate = -1;
				}
			}
		}
		pftlist.nextobj();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////  End of Sowing date algorithm. /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////  Crop phenology  /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculation of down-scaling of lai during crop senescence
/** Follows Bondeau et al. 2007.
 */
double senescence_curve(Pft& pft, double fphu) {

	double senfactor;

	if (pft.shapesenescencenorm)
		senfactor = pow(1-fphu,2) / pow(1 - pft.fphusen, 2) * (1 - pft.flaimaxharvest) + pft.flaimaxharvest;
	else
		senfactor = pow(1 - fphu, 0.5) / pow(1 - pft.fphusen, 0.5) * (1 - pft.flaimaxharvest) + pft.flaimaxharvest;

	return senfactor;
}

/// Initiation of potential heat unit calculation
/** Calculates pvd (required vernalising days), tb (base temperature)
 *  and phu (potential heat units) based on Bondeau et al. 2007.
 *  Dynamic phu calculation based on Lindeskog et al. 2013.
 *  Called on sowing date.
 */
void phu_init(cropphen_struct& ppftcrop, Gridcellpft& gridcellpft, Patch& patch) {

	Pft& pft = gridcellpft.pft;
	Climate& climate = patch.stand.gridcell.climate;
	double phu_last_year = ppftcrop.phu;

	ppftcrop.husum = 0.0;	
	ppftcrop.vrf = 1.0;
	ppftcrop.vdsum = 0;
	ppftcrop.prf = 1.0;

	ppftcrop.pvd = pft.pvd;		// default; kept for TrMi, TePu, TeSb, TrMa, TeSo, TrPe
	ppftcrop.phu = pft.phu;
	ppftcrop.tb = pft.tb;


	if(pft.ifsdautumn) {	// TeWW,TeRa
	
		if(gridcellpft.wintertype) {	// Autumn sowing
			// if neither spring or winter conditions for the past 20 years
			if((gridcellpft.first_autumndate20 == climate.testday_temp || gridcellpft.first_autumndate20 == climate.coldestday) 
				&& gridcellpft.first_autumndate % 365 == gridcellpft.first_autumndate20
				&& (gridcellpft.last_springdate20 == climate.testday_temp || gridcellpft.last_springdate20 == climate.coldestday) 
				&& gridcellpft.last_springdate == gridcellpft.last_springdate20) {

				ppftcrop.pvd = pft.pvd;
				ppftcrop.phu = pft.phu;
			}
			// not all past 20 years without vernendoccurred: too cold
			else if(!(gridcellpft.last_verndate20 == gridcellpft.last_springdate20 + 60 && gridcellpft.last_verndate == gridcellpft.last_verndate20)) {

				// pvd:
				// vernalization (below 12 degrees) is supposed to occur directly at sowing (trg=tempautumn) for TeWW, for TeRa a 20-day lag (tempautumn=17) ??
				// first_autumndate20 occurred before last_verndate20
				if((ppftcrop.sdate < 180 || gridcellpft.last_verndate20 >= 180) && climate.lat >= 0.0 || climate.lat < 0.0) {
					if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
						ppftcrop.pvd = (int)min(60, gridcellpft.last_verndate20 - ppftcrop.sdate);
					else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
						ppftcrop.pvd = (int)max(0, min(60, gridcellpft.last_verndate20 - ppftcrop.sdate - 20));
				}
				// first_autumndate20 occurred after last_verndate20
				else if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
					ppftcrop.pvd = (int)min(60, gridcellpft.last_verndate20 + 365 - ppftcrop.sdate);
				// first_autumndate20 occurred after last_verndate20
				else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
					ppftcrop.pvd = (int)max(0, min(60, gridcellpft.last_verndate20 + 365 - ppftcrop.sdate - 20));

				// phu:
				if (ppftcrop.sdate < 184 + climate.adjustlat) {
					if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
						ppftcrop.phu = max(1700.0, -0.1081 * pow((double)(ppftcrop.sdate - climate.adjustlat),2) + 3.1633 * ((double)(ppftcrop.sdate - climate.adjustlat)) + 2876.9);
					else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
						ppftcrop.phu = max(2100.0, -0.1081 * pow((double)(ppftcrop.sdate - climate.adjustlat),2) + 3.1633 * ((double)(ppftcrop.sdate - climate.adjustlat)) + 3279.7);
				}
				else {
					if(!strncmp(pft.name,"TeWW", strlen("TeWW"))) {
						ppftcrop.phu = max(1700.0, -0.1081 * pow((double)ppftcrop.sdate - 365,2) + 3.1633 * ((double)ppftcrop.sdate - 365) + 2876.9);
						ppftcrop.phu *= 0.8;
					}
					else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
						ppftcrop.phu = max(2100.0, -0.1081 * pow((double)ppftcrop.sdate - 365,2) + 3.1633 * ((double)ppftcrop.sdate - 365) + 3279.7);
				}
			}
		}
		else if (!gridcellpft.wintertype) {	// spring sowing
			// If last_verndate has occurred during the past 20 year (or too warm):
			if(!(gridcellpft.last_verndate20 == ppftcrop.sdate + 60 && gridcellpft.last_verndate == gridcellpft.last_verndate20)) {
				ppftcrop.pvd = (int)min(60, gridcellpft.last_verndate20 - ppftcrop.sdate);
				ppftcrop.phu = 1300.0;
			}
			// If no last_verndate occurred during the past 20 year (too cold):
			else {
				ppftcrop.pvd = 30;
				ppftcrop.phu = 1500.0;
			}
		}
	}
	else if(!strncmp(pft.name,"TrRi", strlen("TrRi")) || !strncmp(pft.name,"TeSf", strlen("TeSf"))) {

		if(!strncmp(pft.name,"TeSf", strlen("TeSf")))
			ppftcrop.phu = min(2000.0, max(1300.0, -700.0 / 90.0 * (ppftcrop.sdate - climate.adjustlat) + 2460.0));

		if (!strncmp(pft.name,"TrRi", strlen("TrRi")) && date.year <= 1) {
			if (patch.stand.gridcell.get_lon() < 60.0 || patch.stand.gridcell.get_lat() > 30.0)
				ppftcrop.phu = 1600.0;
		}
	}

	// Calculation of potential heat units according to local climate.
	if(ifcalcdynamic_phu) {

		// add last year's husum_max to running mean
		if(ppftcrop.husum_max) {
			ppftcrop.nyears_hu_sample++;
			int years = min(ppftcrop.nyears_hu_sample, 10);
			ppftcrop.husum_max_10 = (ppftcrop.husum_max_10 * (years - 1) + ppftcrop.husum_max) / years;
		}
		ppftcrop.husum_max = 0.0;

		ppftcrop.phu_old = ppftcrop.phu;		// phu_old for printout

		// set phu according to running mean
		if(ppftcrop.nyears_hu_sample)
			ppftcrop.phu = max(900.0, 0.9 * ppftcrop.husum_max_10);

		// ... unless past specified time period
		if(ifdyn_phu_limit && date.year >= patch.stand.first_year + 20 && date.year >= nyear_spinup + nyear_dyn_phu)
			ppftcrop.phu = phu_last_year;
	}
}

/// Calculation of harvest index
/** Based on fphu. Restricted by water stress.
 *  Equations are from Neitsch et al. 2002.
 */ 
void calc_hi(Patch& patch, Pft& pft) {

	double wdf, fwdf, hi_save;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

	ppftcrop.hi = pft.hiopt * 100 * ppftcrop.fphu / (100 * ppftcrop.fphu + exp(11.1 - 10.0 * ppftcrop.fphu));	// SWAT 5:2.4.1
	ppftcrop.fhi_phen = ppftcrop.hi / pft.hiopt;

	// Correction of HI according to water stress:
	ppftcrop.demandsum_crop += patch.wdemand;
	if (patchpft.wsupply > patch.wdemand) 
		ppftcrop.supplysum_crop += patch.wdemand; 
	else
		ppftcrop.supplysum_crop += patchpft.wsupply;

	if(ppftcrop.demandsum_crop > 0.0)							
		wdf = 100.0 * ppftcrop.supplysum_crop / ppftcrop.demandsum_crop;	// SWAT 5:3.3.2	: aetsum/petsum
	else
		wdf = 100.0;

	fwdf = wdf / (wdf + exp(6.13 - 0.0883 * wdf));		// (SWAT 5:3.3.1) 

	hi_save = ppftcrop.hi;

	ppftcrop.hi = ppftcrop.fhi_phen * ((pft.hiopt - pft.himin) * fwdf + pft.himin);											

	if(ppftcrop.hi > 0.0 && hi_save)
		ppftcrop.fhi_water = ppftcrop.hi / hi_save;

	ppftcrop.fhi = ppftcrop.fhi_phen * ppftcrop.fhi_water;
}


/// Calculation of accumulated of heat units
/** Accumulation of heat units during sampling period used for calculation of dynamic phu if DYNAMIC_PHU defined.
 *  Equation is from Neitsch et al. 2002.
 */ 
void calc_hu(Patch& patch, Pft& pft) {

	double hu;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Climate& climate = patch.stand.gridcell.climate;

	// calculation av fphu:
#if defined MAXHUTEMP
	hu = max(0.0,min(climate.temp, 30.0) - ppftcrop.tb);
#else
	hu = max(0,climate.temp - ppftcrop.tb);
#endif

	// account for vernalization if needs for vernalization not yet satisfied	//trg=tb for crops other than TeWW and TeRa and don't enter here
	if (climate.temp < pft.trg && ppftcrop.vdsum < ppftcrop.pvd) {				//trg=12 for TeWW & TeRa
	
		ppftcrop.vdsum++;
		ppftcrop.vrf = min(1.0, (double)ppftcrop.vdsum / (double)ppftcrop.pvd);

		// vernalisation reduction factor has no effect once temp>trg even if needs are not satisfied
		// no effect as well if temp>trg at the beginning of the growing season...
		hu *= ppftcrop.vrf;																		
	}

	// account for response to photoperiod
	ppftcrop.prf = (1 - pft.psens) * min(1.0, max(0.0, (climate.daylength_save[date.day] - pft.pb) / (pft.ps - pft.pb))) + pft.psens;
	hu *= ppftcrop.prf;

	if(date.day == ppftcrop.sdate)
		ppftcrop.husum = 0.0;

	// Accumulate heat units during growing period
	if (ppftcrop.growingseason) {
		// daily effective temperature sum (degree-days)
		ppftcrop.husum += hu;

		// phenological scale (fraction of growing season)
		ppftcrop.fphu = min(1.0, ppftcrop.husum / ppftcrop.phu);		// SWAT 5:2.1.11
	}
	
	// Sample heat units for dynamic phu calculation
	if(ifcalcdynamic_phu) {

		if(date.day == ppftcrop.sdate) {
			ppftcrop.hu_samplingperiod = true;
			ppftcrop.hu_samplingdays = 0;
			ppftcrop.husum_sampled = 0;
		}

		if(ppftcrop.hu_samplingperiod) {

			ppftcrop.husum_sampled += hu;
			ppftcrop.hu_samplingdays++;
		
			if(date.day == ppftcrop.hucountend) {

				ppftcrop.husum_sampled -= hu;					// Don't count the hu's on last day
				ppftcrop.husum_max = ppftcrop.husum_sampled;

				ppftcrop.hu_samplingperiod=false;
			}
		}
	}
}

/// Handles heat unit and harvest index calculation and identifies harvest, senescence and intercrop events.
/** Accumulation of heat units during sampling period used for calculation of dynamic phu if DYNAMIC_PHU defined.
 *  Sets patchpft.cropphen variables growingseason, hdate, intercropseason and senescence
 */ 
void crop_phenology(Patch& patch) 
{
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& patchpft = patch.pft.getobj();
		Pft& pft=patchpft.pft;
		Gridcell& gridcell = patch.stand.gridcell;
		Climate& climate = gridcell.climate;
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];
		double hu = 0.0;

		if(patch.stand.pft[pft.id].active) {
		if(pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
			ppftcrop.growingseason_ystd = ppftcrop.growingseason;

			// resets on first day of the year:	
			if(date.day == 0) {
				ppftcrop.fphu_harv = -1.0;
				ppftcrop.fhi_harv = -1.0;
				ppftcrop.sdate_harv = -1;
				ppftcrop.nsow = 0;
				ppftcrop.sendate = -1;
				ppftcrop.nharv = 0;

				ppftcrop.sownlastyear=false;

				for(int i=0;i<2;i++) {
					ppftcrop.sdate_harvest[i] = -1;
					ppftcrop.hdate_harvest[i] = -1;
					//ppftcrop.fphu_harvest[i] = -1.0;
					//ppftcrop.fhi_harvest[i] = -1.0;
					ppftcrop.sdate_thisyear[i] = -1;
				}
			}

			// initiations on sowing day:
			if(date.day == ppftcrop.sdate) {
				ppftcrop.fphu = 0.0;
				ppftcrop.fhi = 0.0;
				ppftcrop.fhi_phen = 0.0;
				ppftcrop.fhi_water = 1.0;
				ppftcrop.hdate = -1;
				ppftcrop.bicdate = -1;

				ppftcrop.growingseason = true;
				ppftcrop.growingdays = 0;
				ppftcrop.nsow++;

				if(ppftcrop.nsow == 1)
					ppftcrop.sdate_thisyear[0] = ppftcrop.sdate;	
				else if(ppftcrop.nsow == 2)
					ppftcrop.sdate_thisyear[1] = ppftcrop.sdate;

				// calculate pvd, phu & tb
				phu_init(ppftcrop, gridcellpft, patch);
			}

			// Calculation of accumulated heat units and harvest index from sowing to maturity
			if (ppftcrop.growingseason) 
			{
				ppftcrop.senescence_ystd = ppftcrop.senescence;
				ppftcrop.hi_ystd = ppftcrop.hi;
				ppftcrop.intercropseason = false;
				ppftcrop.growingdays++;

				// check if harvest is prescribed
				bool force_harvest = forceharvestdates && pft.forceharvestdate && gridcellpft.hdate_force != -1 && date.day == gridcellpft.hdate_force;

				// before maturity is reached
				if(ppftcrop.husum < ppftcrop.phu && dayinperiod(date.day, ppftcrop.sdate, stepfromdate(ppftcrop.hlimitdate, -1)) && !force_harvest) {

					// count accumulated heat units after sowing date
					calc_hu(patch, pft);

					//  test for senescence
					if (ppftcrop.fphu >= pft.fphusen) {

						if(ppftcrop.senescence_ystd == false)
							ppftcrop.sendate = date.day;
						ppftcrop.senescence = true;
					}

					// calculated harvest index
					calc_hi(patch, pft);

				}
				else {	// harvest

					// save today as harvest day
					ppftcrop.hdate = date.day;

					ppftcrop.growingseason = false;
					ppftcrop.intercropseason = false;		
					ppftcrop.senescence = false;

					// set start of intercrop grass growth
					ppftcrop.bicdate = stepfromdate(ppftcrop.hdate, 15);

					// count number of harvest events this year
					ppftcrop.nharv++;

					// Save phenological values and dates at harvest:
					ppftcrop.fphu_harv = ppftcrop.fphu;
					ppftcrop.fhi_harv = ppftcrop.fhi;
					ppftcrop.sdate_harv = ppftcrop.sdate;
					ppftcrop.lgp = ppftcrop.growingdays;

					// allowing saving at two harvests per year
					if(ppftcrop.nharv == 1) {
						ppftcrop.sdate_harvest[0] = ppftcrop.sdate;
						ppftcrop.hdate_harvest[0] = date.day;
						//ppftcrop.fphu_harvest[0] = ppftcrop.fphu;
						//ppftcrop.fhi_harvest[0] = ppftcrop.fhi;
						if(ppftcrop.sdate > date.day)							
							ppftcrop.sownlastyear = true;
					}
					else if(ppftcrop.nharv == 2) {
						ppftcrop.sdate_harvest[1] = ppftcrop.sdate;
						ppftcrop.hdate_harvest[1] = date.day;
						//ppftcrop.fphu_harvest[1] = ppftcrop.fphu;
						//ppftcrop.fhi_harvest[1] = ppftcrop.fhi;
					}

					ppftcrop.demandsum_crop = 0.0;
					ppftcrop.supplysum_crop = 0.0;

					if(pft.ifsdprec) {	// TeCo,TrMi,TrMa,TrPe; with NEWSOWINGDATE: all crops			
						ppftcrop.sdate = -1;
						ppftcrop.eicdate = -1;
					}

				} //end harvest
			}  //from sowing has taken place until harvest day

			// continue sampling heat units from hdate until last sampling date
			if(ifcalcdynamic_phu && ppftcrop.growingseason == false && ppftcrop.hu_samplingperiod)
				calc_hu(patch, pft);

			if(pft.intercrop == NATURALGRASS) {

				if(!ppftcrop.intercropseason && date.day == ppftcrop.bicdate)
					ppftcrop.intercropseason = true;

				if(date.day == ppftcrop.eicdate) {
					ppftcrop.intercropseason = false;
				}
			}
		}
		else if(pft.phenology == ANY) { // crop grasses using standard guess phenology calculation
		
			if(patch.stand.pftid != pft.id) {
				cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

				if(date.day == 0)
					ppftcrop.nharv = 0;

				if(date.day == patch.pft[patch.stand.pftid].cropphen->bicdate)
					ppftcrop.growingseason = true;
				else if(date.day == patch.pft[patch.stand.pftid].cropphen->eicdate)
					ppftcrop.growingseason = false;
			}
		}
		}

		patch.pft.nextobj();
	}
}


/// Updates crop phen from yesterday's lai_daily
/** True crops derive phen from yesterday's fpc_daily,
 *   assuming only one crop individual.
 *  Intercrop grass and pasture grass grown in crop stands use gdd5.
 *   and is treated similar to pasture grass.
 */ 
void leaf_phenology_crop(Pft& pft, Patch& patch) 
{
	Gridcell& gridcell = patch.stand.gridcell;
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

	if (pft.phenology == CROPGREEN && patch.stand.pft[pft.id].active) {
		if (ppftcrop.growingseason) {

			Vegetation& vegetation = patch.vegetation;

			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();

				if(indiv.pft.id == pft.id) {
					if(indiv.fpc)
						patchpft.phen = indiv.fpc_daily / indiv.fpc;
					else
						patchpft.phen = 0.0;
				}

				vegetation.nextobj();
			}
		}
		else if (date.day == ppftcrop.hdate)
			patchpft.phen = 0.0;
	}
	else if(pft.phenology == ANY && patch.stand.pft[pft.id].active) { // crop grasses using standard guess phenology calculation

		if(patch.stand.pftid == pft.id // pasture grass in crop stand
			|| patch.stand.hasgrassintercrop && (patch.pft[patch.stand.pftid].cropphen->intercropseason // intercrop grass
			|| date.day == patch.pft[patch.stand.pftid].cropphen->eicdate)) {							// intercrop grass

			if(patch.stand.pftid != pft.id) {

				if(date.day == patch.pft[patch.stand.pftid].cropphen->bicdate) {
					patch.stand.gdd0_intercrop = climate.gdd5;
				}
				if(date.day == patch.pft[patch.stand.pftid].cropphen->eicdate) {
					patch.stand.gdd0_intercrop = 0.0;
					patchpft.phen = 0.0;
				}
			}

			// reset stand.gdd0_intercrop same day as gdd5
			if (climate.lat >= 0.0 && date.day == COLDEST_DAY_NHEMISPHERE || climate.lat < 0.0 && date.day == COLDEST_DAY_SHEMISPHERE
					|| climate.gdd5 == 0.0)
				patch.stand.gdd0_intercrop = 0.0;

			if(ppftcrop.growingseason) {	// includes bicdate, not eicdate
			
				if(patch.stand.pftid == pft.id || gridcell.pft[patch.stand.pftid].sowing_restriction)	// Normal grass growth: gives identical result to natural stands.
					patchpft.phen = min(1.0, climate.gdd5 / pft.phengdd5ramp);
				else if(patch.stand.gdd0_intercrop > 0.0)
					patchpft.phen = min(1.0, (climate.gdd5 - patch.stand.gdd0_intercrop) / (pft.phengdd5ramp * 0.9)); // intercrop grass
				else
					patchpft.phen = min(1.0, (climate.gdd5 - patch.stand.gdd0_intercrop) / pft.phengdd5ramp); // intercrop grass

				if(patchpft.phen < 0.0)
					patchpft.phen = 0.0;

				// raingreen phenology
				if (patchpft.wscal < pft.wscal_min) {

					patchpft.phen = 0.0;
					patch.stand.gdd0_intercrop = climate.gdd5;
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////  End of crop phenology  ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////  Crop allocation  ///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Updates patch.members fpc_total and fpc_rescale for crops (to be called after crop_phenology())
void update_patch_fpc(Patch& patch) {

	if(patch.stand.landcover == CROPLAND) {
		Vegetation& vegetation = patch.vegetation;
		patch.fpc_total = 0.0;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();

			if(patch.pft[indiv.pft.id].cropphen->growingseason == true)
				patch.fpc_total += indiv.fpc;
			vegetation.nextobj();
		}
		// Calculate rescaling factor to account for overlap between populations/
		// cohorts/individuals (i.e. total FPC > 1)
		// necessary to undate here after growingseason updated
		patch.fpc_rescale = 1.0 / max(patch.fpc_total, 1.0);
	}
}


/// Updates lai_daily and fpc_daily from daily grs_cmass_leaf-value
/** lai during senescence declines according the function senescence_curve()
 */
void lai_crop(Patch& patch) {

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();

	while (vegetation.isobj)
	{
		Individual& indiv = vegetation.getobj();
		cropindiv_struct& cropindiv = *(indiv.get_cropindiv());
		Patchpft& patchpft = patch.pft[indiv.pft.id];
		cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

		if(indiv.pft.phenology == CROPGREEN) {

			if(ppftcrop.growingseason) {

				if(!ppftcrop.senescence)
					indiv.lai_daily = cropindiv.grs_cmass_leaf * indiv.pft.sla;
				else
					// Follow the senescence curve from leaf cmass at senescence (cmass_leaf_sen):
					indiv.lai_daily = cropindiv.cmass_leaf_sen * indiv.pft.sla * senescence_curve(indiv.pft, ppftcrop.fphu);

				if(indiv.lai_daily < 0.0)
					indiv.lai_daily = 0.0;

				indiv.fpc_daily = 1.0 - lambertbeer(indiv.lai_daily);

				indiv.lai_indiv_daily = indiv.lai_daily;
			}
			else if(date.day == ppftcrop.hdate) {
				indiv.lai_daily = 0.0;
				indiv.lai_indiv_daily = 0.0;
				indiv.fpc_daily = 0.0;
			}

			if(!(indiv.lai_daily>=0.0 && indiv.lai_daily<=20.0))//Test for unrealistically high lai.
if(!SUPPRESSLARGEOUTPUT)	
				dprintf("In lai_crop() stand %d pft %d year %d day %d: senescence=%d, grs_cmass_leaf=%f, grs_cmass_ho=%f, lai_daily=%f, out of bounds !\n", patch.stand.id, indiv.pft.id, date.year-nyear_spinup+1901, date.day, ppftcrop.senescence, cropindiv.grs_cmass_leaf, cropindiv.grs_cmass_ho, indiv.lai_daily);
		}
		vegetation.nextobj();
	}
}

/// Turnover function for continuous grass, to be called from any day of the year from allocation_crop_daily().
void turnover_grass(Individual& indiv) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());
	Patchpft& patchpft = indiv.patchpft();

	double cmass_leaf_inc = cropindiv.grs_cmass_leaf - indiv.cmass_leaf_post_turnover;
	double cmass_root_inc = cropindiv.grs_cmass_root - indiv.cmass_root_post_turnover;

	double grs_npp = cmass_leaf_inc + cmass_root_inc;
#ifdef GRASS_SEED_CMASS
	grs_npp -= CMASS_SEED;
#endif
	double cmass_leaf_pre_turnover = cropindiv.grs_cmass_leaf;
	double cmass_root_pre_turnover = cropindiv.grs_cmass_root;
	double cton_leaf_bg = indiv.cton_leaf(false);
	double cton_root_bg = indiv.cton_root(false);

	indiv.nstore_longterm += indiv.nstore_labile;
	indiv.nstore_labile = 0.0;

	turnover(indiv.pft.turnover_leaf, indiv.pft.turnover_root,
		indiv.pft.turnover_sap, indiv.pft.lifeform, indiv.pft.landcover,
		indiv.cropindiv->grs_cmass_leaf, indiv.cropindiv->grs_cmass_root, indiv.cmass_sap, indiv.cmass_heart,
		indiv.nmass_leaf, indiv.nmass_root, indiv.nmass_sap, indiv.nmass_heart,
		patchpft.litter_leaf,
		patchpft.litter_root,
		patchpft.nmass_litter_leaf,
		patchpft.nmass_litter_root,
		indiv.nstore_longterm, 
		true);

	indiv.cmass_leaf_post_turnover = cropindiv.grs_cmass_leaf;
	indiv.cmass_root_post_turnover = cropindiv.grs_cmass_root;

	if (indiv.nstore_longterm > indiv.max_n_storage) {
						
		// Nitrogen stored above maximum will be returned to litter
		double nsurplus = indiv.nstore_longterm - indiv.max_n_storage;
		indiv.nstore_longterm -= nsurplus;

		// Return surplus nitrogen to litter
		patchpft.nmass_litter_leaf += nsurplus * (indiv.pft.turnover_leaf / (indiv.pft.turnover_leaf + indiv.pft.turnover_root));
		patchpft.nmass_litter_root += nsurplus * (indiv.pft.turnover_root / (indiv.pft.turnover_leaf + indiv.pft.turnover_root));
	}

	// Nitrogen longtime storage
	// Nitrogen approx retranslocated next season
	double retransn_nextyear = cmass_leaf_pre_turnover * indiv.pft.turnover_leaf / cton_leaf_bg * nrelocfrac +
		cmass_root_pre_turnover * indiv.pft.turnover_root / cton_root_bg * nrelocfrac;

	// Max longterm nitrogen storage
	indiv.max_n_storage = min(cmass_root_pre_turnover * indiv.pft.fnstorage, 
		(max(0.0, cmass_leaf_inc) + max(0.0, cmass_root_inc)) * indiv.densindiv) / cton_leaf_bg;

	// Scale this year productivity to max storage
	if (grs_npp > 0.0) {
		indiv.scale_n_storage = max(0.5 * indiv.max_n_storage, indiv.max_n_storage - retransn_nextyear) * cton_leaf_bg / grs_npp;
	}

	indiv.nstore_labile = indiv.nstore_longterm;
	indiv.nstore_longterm = 0.0;
}

/// Daily growth routine for crops
/** Allocates daily npp to leaf, roots and harvestable organs
 *  Requires updated value of fphu and hi.
 *  Equations are from Neitsch et al. 2002.
 */
void allocation_crop_daily(Patch& patch) {

	double froot, fleaf;
	double grs_cmass_root_old;
	double grs_cmass_leaf_old;
	double grs_cmass_ho_old;
	double grs_cmass_ag;

	Vegetation& vegetation = patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {

		Individual& indiv = vegetation.getobj();
		cropindiv_struct& cropindiv = *(indiv.get_cropindiv());
		Patchpft& patchpft = patch.pft[indiv.pft.id];
		cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

		if(date.day == 0) {

			cropindiv.ycmass_plant = 0.0;
			cropindiv.ycmass_leaf = 0.0;
			cropindiv.ycmass_root = 0.0;
			cropindiv.ycmass_ho = 0.0;
			cropindiv.ycmass_agpool = 0.0;	

			cropindiv.harv_cmass_plant = 0.0;
			cropindiv.harv_cmass_root = 0.0;
			cropindiv.harv_cmass_leaf = 0.0;
			cropindiv.harv_cmass_ho = 0.0;
			cropindiv.harv_cmass_agpool = 0.0;

			cropindiv.cmass_ho_harvest[0] = 0.0;
			cropindiv.cmass_ho_harvest[1] = 0.0;

			cropindiv.cmass_leaf_max = 0.0;

			if(indiv.pft.phenology == ANY && !cropindiv.isintercropgrass) {		//zero of normal cc3g/cc4g-growth arbitrarily at new year
			
				cropindiv.grs_cmass_plant = 0.0;
				cropindiv.grs_cmass_root = 0.0;
				cropindiv.grs_cmass_ho = 0.0;
				cropindiv.grs_cmass_leaf = 0.0;
				cropindiv.grs_cmass_agpool = 0.0;
			}
		}

		// true crop allocation
		if(indiv.pft.phenology == CROPGREEN) {

			if(ppftcrop.growingseason) {

				cropindiv.dcmass_plant = 0.0;

#ifdef DELAYED_SEEDCARBON
				// Seed carbon; portion the seed carbon over a 10-day period.
				if(dayinperiod(date.day, pppftcrop.sdate, (patchpft.cropphen->sdate + 9)) % 365 ) {			
					cropindiv.grs_cmass_plant += 0.1 * CMASS_SEED;
					cropindiv.ycmass_plant += 0.1 * CMASS_SEED;
				}
#else
				// add seed carbon on sowing date
				if(date.day == ppftcrop.sdate) {

					cropindiv.grs_cmass_plant += CMASS_SEED;
					cropindiv.ycmass_plant += CMASS_SEED;
					cropindiv.dcmass_plant += CMASS_SEED;

					// This flux will be balancing litter fluxes for the NEXT year.
					patch.fluxes.report_flux(Fluxes::SEEDC, -CMASS_SEED);
				}
#endif
				// add today's npp
				cropindiv.dcmass_plant += indiv.dnpp;
				cropindiv.grs_cmass_plant += indiv.dnpp;
				cropindiv.ycmass_plant += indiv.dnpp;

				// allocation to roots
				froot = indiv.pft.frootstart -(indiv.pft.frootstart - indiv.pft.frootend) * ppftcrop.fphu;	// SWAT 5:2,1,21	
				grs_cmass_root_old = cropindiv.grs_cmass_root;
				cropindiv.grs_cmass_root = froot * cropindiv.grs_cmass_plant;
				cropindiv.dcmass_root = cropindiv.grs_cmass_root - grs_cmass_root_old;
				cropindiv.ycmass_root += cropindiv.dcmass_root;

				// allocation to harvestable organs
				grs_cmass_ag = (1.0-froot) * cropindiv.grs_cmass_plant;
				grs_cmass_ho_old = cropindiv.grs_cmass_ho;

				if(indiv.pft.hiopt <= 1.0)
					cropindiv.grs_cmass_ho = ppftcrop.hi * grs_cmass_ag;									// SWAT 5:2.4.2, 5:2.4.4
				else	// below-ground harvestable organs
					cropindiv.grs_cmass_ho = (1.0 - 1.0 / (1.0 + ppftcrop.hi)) * cropindiv.grs_cmass_plant;	// SWAT 5:2.4.3 8 

				cropindiv.dcmass_ho = cropindiv.grs_cmass_ho - grs_cmass_ho_old;	
				cropindiv.ycmass_ho += cropindiv.dcmass_ho;	

				// allocation to leaves
				grs_cmass_leaf_old = cropindiv.grs_cmass_leaf;	
				cropindiv.grs_cmass_leaf = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root - cropindiv.grs_cmass_ho;
 
				if(cropindiv.grs_cmass_plant > 0.0)
					fleaf = cropindiv.grs_cmass_leaf / cropindiv.grs_cmass_plant;
				cropindiv.dcmass_leaf = cropindiv.grs_cmass_leaf - grs_cmass_leaf_old;
				cropindiv.ycmass_leaf += cropindiv.dcmass_leaf;

				// allocation to above-ground pool (currently not used)
				cropindiv.dcmass_agpool = cropindiv.dcmass_plant - cropindiv.dcmass_root - cropindiv.dcmass_leaf - cropindiv.dcmass_ho;		
				cropindiv.grs_cmass_agpool = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root - cropindiv.grs_cmass_leaf - cropindiv.grs_cmass_ho;
				cropindiv.ycmass_agpool = cropindiv.ycmass_plant - cropindiv.ycmass_root - cropindiv.ycmass_leaf - cropindiv.ycmass_ho;

				if(cropindiv.grs_cmass_agpool < 10e-10)
					cropindiv.grs_cmass_agpool = 0,0;
				if(cropindiv.ycmass_agpool < 10e-10)
					cropindiv.ycmass_agpool = 0,0;

				// save this year's maximum leaf carbon mass
				if(cropindiv.grs_cmass_leaf > cropindiv.cmass_leaf_max)	
					cropindiv.cmass_leaf_max = cropindiv.grs_cmass_leaf;

				// save leaf carbon mass at the beginning of senescence
				if(date.day == ppftcrop.sendate)
					cropindiv.cmass_leaf_sen = cropindiv.grs_cmass_leaf;

				// Check that no plant cmass is negative, if so, zero cmass and correct C fluxes
				indiv.check_C_mass();
			}
			else if(date.day == ppftcrop.hdate) {

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;

				if(ppftcrop.nharv == 1)
					cropindiv.cmass_ho_harvest[0] = cropindiv.grs_cmass_ho;
				else if(ppftcrop.nharv == 2)
					cropindiv.cmass_ho_harvest[1] = cropindiv.grs_cmass_ho;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.gridcell.LC_updated && patchpft.cropphen->nharv == 1)
						scale_indiv(indiv, true);
					harvest_crop(indiv, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass, true);
					patch.is_litter_day = true;
				}

				cropindiv.grs_cmass_plant = 0.0;
				cropindiv.grs_cmass_root = 0.0;
				cropindiv.grs_cmass_ho = 0.0;
				cropindiv.grs_cmass_leaf = 0.0;
				cropindiv.grs_cmass_agpool = 0.0;
				cropindiv.cmass_leaf_sen = 0.0;	

				cropindiv.dcmass_plant = 0.0;
				cropindiv.dcmass_root = 0.0;
				cropindiv.dcmass_ho = 0.0;
				cropindiv.dcmass_leaf = 0.0;
				cropindiv.dcmass_agpool = 0.0;

			}
		}
		// crop grass allocation
		// NB: Only intercrop grass enters here ! normal cc3g/cc4g grass is treated just like natural grass 
		else if(indiv.pft.phenology == ANY && indiv.pft.id != patch.stand.pftid) {

			if(ppftcrop.growingseason) {

				cropindiv.dcmass_plant = 0.0;
#ifdef GRASS_SEED_CMASS
				// add seed carbon
				if(!indiv.continous_grass() && date.day == patch.pft[patch.stand.pftid].cropphen->bicdate 
					|| indiv.continous_grass() && date.day == stepfromdate(indiv.last_turnover_day, 1)) {

					cropindiv.grs_cmass_plant += CMASS_SEED;
					cropindiv.ycmass_plant += CMASS_SEED;
					cropindiv.dcmass_plant += CMASS_SEED;

					// This flux will be balancing litter fluxes for the NEXT year.
					patch.fluxes.report_flux(Fluxes::SEEDC, -CMASS_SEED);

					indiv.last_turnover_day = -1;
				}
#endif									
				cropindiv.dcmass_plant += indiv.dnpp;
				cropindiv.grs_cmass_plant += indiv.dnpp;		
				cropindiv.ycmass_plant += indiv.dnpp;

				indiv.ltor = indiv.wscal_mean * indiv.pft.ltor_max;

				// allocation to roots
				froot = 1.0 / (1.0 + indiv.ltor);
				grs_cmass_root_old = cropindiv.grs_cmass_root;	

				//Cumulative wscal-dependent root increase						
				cropindiv.grs_cmass_root = froot * cropindiv.grs_cmass_plant;
				cropindiv.dcmass_root = cropindiv.grs_cmass_root - grs_cmass_root_old;
				cropindiv.ycmass_root += cropindiv.dcmass_root;

				// allocation to leaves
				fleaf = 1.0 - froot;
				grs_cmass_leaf_old = cropindiv.grs_cmass_leaf;
				cropindiv.grs_cmass_leaf = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root;
				cropindiv.dcmass_leaf = cropindiv.grs_cmass_leaf - grs_cmass_leaf_old;
				cropindiv.ycmass_leaf += cropindiv.dcmass_leaf;

				// Check that no plant cmass is negative, if so, zero cmass and correct C fluxes
				indiv.check_C_mass();
			}
			else if(date.day == patch.pft[patch.stand.pftid].get_cropphen()->eicdate) {

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;	
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;	
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;	
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;		
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;

				ppftcrop.nharv++;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.gridcell.LC_updated && patchpft.cropphen->nharv == 1)
						scale_indiv(indiv, true);
					harvest_crop(indiv, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass, true);
					patch.is_litter_day = true;
				}
				else {
					cropindiv.grs_cmass_root = 0.0;
					cropindiv.grs_cmass_ho = 0.0;
					cropindiv.grs_cmass_leaf = 0.0;
					cropindiv.grs_cmass_agpool = 0.0;
				}

				cropindiv.grs_cmass_plant = cropindiv.grs_cmass_root + cropindiv.grs_cmass_leaf;

				cropindiv.dcmass_plant = 0.0;
				cropindiv.dcmass_root = 0.0;
				cropindiv.dcmass_ho = 0.0;
				cropindiv.dcmass_leaf = 0.0;
				cropindiv.dcmass_agpool = 0.0;
			}
			
			if(indiv.continous_grass() && indiv.is_turnover_day()) {

				indiv.last_turnover_day = date.day;

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;	
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;	
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;	
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;		
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;

				ppftcrop.nharv++;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.gridcell.LC_updated && patchpft.cropphen->nharv == 1)
						scale_indiv(indiv, true);

					turnover_grass(indiv);
					patch.is_litter_day = true;
				}
				else {
					cropindiv.grs_cmass_root = 0.0;
					cropindiv.grs_cmass_ho = 0.0;
					cropindiv.grs_cmass_leaf = 0.0;
				}

				cropindiv.grs_cmass_plant = cropindiv.grs_cmass_root + cropindiv.grs_cmass_leaf;
				cropindiv.grs_cmass_agpool = 0.0;

				cropindiv.dcmass_plant = 0.0;
				cropindiv.dcmass_root = 0.0;
				cropindiv.dcmass_ho = 0.0;
				cropindiv.dcmass_leaf = 0.0;
				cropindiv.dcmass_agpool = 0.0;
			}
		}
		vegetation.nextobj();
	}

	return;
}

/// Handles daily crop allocation and daily lai calculation
/** Simple allocation based on heat unit accumulation.
 *  LAI is set directly after allocation from leaf carbon mass.
 */
void crop_growth_daily(Patch& patch) {

	// allocate daily npp to leaf, roots and harvestable organs
	allocation_crop_daily(patch);

	// update patchpft.lai_daily and fpc_daily
	lai_crop(patch);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////  End of crop allocation  //////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////  Landcover harvest functions  ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Harvest function used for managed forest and for clearing natural vegetation at land use change
/** A fraction of trees is cut down (frac_cut)
 *  A fraction of wood is harvested (pft.harv_eff) and returned as acflux_harvest
 *  A fraction of harvested wood (pft.harvest_slow_frac) is returned as harvested_products_slow
 *  The rest, including leaves and roots, is returned as litter.
 *  Called from landcover_dynamics() first day of the year if any natural vegetation is transferred to another land use.
 *  INPUT PARAMETER
 *  \param frac_cut					fraction of trees cut   
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
void harvest_wood(Harvest_CN& i, Pft& pft, bool alive, double frac_cut) {

	double harvest = 0.0;
	double residue_outtake = 0.0;

	// only harvest trees
	if(pft.lifeform == GRASS)
		return;

	// all root carbon and nitrogen goes to litter
	if(alive) {

		i.litter_root += i.cmass_root * frac_cut;
		i.nmass_litter_root += i.nmass_root * frac_cut;
		i.nmass_litter_root += (i.nstore_labile + i.nstore_longterm) * frac_cut;
	}
	i.cmass_root *= (1.0 - frac_cut);
	i.nmass_root *= (1.0 - frac_cut);
	i.nstore_labile *= (1.0 - frac_cut);
	i.nstore_longterm *= (1.0 - frac_cut);

	if(alive) {	

		// Carbon:

		// harvested wood
		harvest = pft.harv_eff * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

		// harvested products not consumed (oxidised) this year put into harvested_products_slow
		if(ifslowharvestpool) {
			i.harvested_products_slow += harvest * pft.harvest_slow_frac;
			harvest = harvest * (1 - pft.harvest_slow_frac);
		}

		// harvested products consumed (oxidised) this year put into acflux_harvest
		i.acflux_harvest += harvest;				

		// removed residues are oxidised
		residue_outtake += pft.res_outtake * i.cmass_leaf * frac_cut;
		residue_outtake += pft.res_outtake * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * (1 - pft.harv_eff) * frac_cut;
		i.acflux_harvest += residue_outtake;															

		// not removed residues are put into litter
		i.litter_leaf += i.cmass_leaf * (1-pft.res_outtake) * frac_cut;												
		i.litter_sap += i.cmass_sap * (1-pft.res_outtake) * (1 - pft.harv_eff) * frac_cut;
		i.litter_heart += (i.cmass_heart-i.cmass_debt) * (1-pft.res_outtake) * (1 - pft.harv_eff) * frac_cut;

		// unharvested trees:
		i.cmass_leaf *= (1.0 - frac_cut);
		i.cmass_sap *= (1.0 - frac_cut);
		i.cmass_heart *= (1.0 - frac_cut);
		i.cmass_debt *= (1.0 - frac_cut);		

		//Nitrogen:

		// harvested products
		harvest = pft.harv_eff * (i.nmass_sap + i.nmass_heart) * frac_cut;

		// harvested products not consumed this year put into harvested_products_slow_nmass
		if(ifslowharvestpool) {
			i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
			harvest = harvest * (1 - pft.harvest_slow_frac);
		}

		// harvested products consumed this year put into anflux_harvest
		i.anflux_harvest += harvest;

		// removed residues are oxidised
		residue_outtake = 0.0;
		residue_outtake += pft.res_outtake * i.nmass_leaf * frac_cut;
		residue_outtake += pft.res_outtake * (i.nmass_sap + i.nmass_heart) * (1 - pft.harv_eff) * frac_cut;
		i.anflux_harvest += residue_outtake;															

		// not removed residues are put into litter
		i.nmass_litter_leaf += i.nmass_leaf * (1 - pft.res_outtake) * frac_cut;												
		i.nmass_litter_sap += i.nmass_sap * (1 - pft.res_outtake) * (1 - pft.harv_eff) * frac_cut;
		i.nmass_litter_heart += i.nmass_heart * (1 - pft.res_outtake) * (1 - pft.harv_eff) * frac_cut;

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
void harvest_wood(Individual& indiv,Pft& pft, bool alive, double frac_cut) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_wood(indiv_cp, pft, alive, frac_cut);

	indiv_cp.copy_to_indiv(indiv);
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

	if(ifslowharvestpool) {
		i.harvested_products_slow += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	if(alive)
		i.acflux_harvest += harvest;
	i.cmass_leaf -= harvest;

	// Nitrogen:
	// Reduced removal of N relative to C during grazing.
	double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
	harvest = pft.harv_eff * i.nmass_leaf * N_harvest_scale;

	if(ifslowharvestpool) {
		i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
		harvest = harvest * (1 - pft.harvest_slow_frac);
	}
	i.anflux_harvest += harvest;
	i.nmass_leaf -= harvest;

#if defined GRASSFORCROP
	if (alive) {
		// Carbon:
		residue_outtake = pft.res_outtake * i.cmass_leaf;	// res_outtake currently set to 0.0, 
		i.acflux_harvest += residue_outtake;						// could be used for burning		
		i.cmass_leaf -= residue_outtake;

		// Nitrogen:
		residue_outtake = pft.res_outtake * i.nmass_leaf;
		i.anflux_harvest += residue_outtake;								
		i.nmass_leaf -= residue_outtake;
	}
#endif
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
void harvest_pasture(Individual& indiv, Pft& pft, bool alive) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_pasture(indiv_cp, pft, alive);

	indiv_cp.copy_to_indiv(indiv);

}

/// Harvest function for cropland, including true crops, intercrop grass 
/**   and pasture grass grown in cropland.
 *  Function for balancing carbon and nitrogen fluxes from last year's growth if old-style harvest is selected (HARVEST_GRSC defined),
 *  or, alternatively, this years harvested carbon and nitrogen.
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

	if(pft.phenology==CROPGREEN) {

	// all root carbon and nitrogen goes to litter
		if(i.cmass_root > 0.0)
			i.litter_root += i.cmass_root;
		i.cmass_root = 0.0;

		if(i.nmass_root > 0.0)
			i.nmass_litter_root += i.nmass_root;
		if(i.nstore_labile > 0.0)
			i.nmass_litter_root += i.nstore_labile;
		if(i.nstore_longterm > 0.0)
			i.nmass_litter_root += i.nstore_longterm;
		i.nmass_root = 0.0;
		i.nstore_labile = 0.0;
		i.nstore_longterm = 0.0;

		// harvest of harvestable organs
		// Carbon:
		if(i.cmass_ho > 0.0) {
			// harvested products
			harvest = pft.harv_eff * i.cmass_ho;

			// not removed harvestable organs are put into litter
			if(pft.aboveground_ho)
				i.litter_leaf += (i.cmass_ho - harvest);
			else
				i.litter_root += (i.cmass_ho - harvest);

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if(ifslowharvestpool) {
				i.harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed (oxidised) this year put into acflux_harvest
			i.acflux_harvest += harvest;
		}
		i.cmass_ho = 0.0;

		// Nitrogen:
		if(i.nmass_ho > 0.0) {

			// harvested products
			harvest = pft.harv_eff * i.nmass_ho;

			// not removed harvestable organs are put into litter
			if(pft.aboveground_ho)
				i.nmass_litter_leaf += (i.nmass_ho - harvest);
			else
				i.nmass_litter_root += (i.nmass_ho - harvest);			

			// harvested products not consumed this year put into harvested_products_slow_nmass
			if(ifslowharvestpool) {
				i.harvested_products_slow_nmass += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}

			// harvested products consumed this year put into anflux_harvest
			i.anflux_harvest += harvest;
		}
		i.nmass_ho = 0.0;

		// residues
		// Carbon
		if ((i.cmass_leaf + i.cmass_agpool) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.cmass_leaf + i.cmass_agpool);
			i.acflux_harvest += residue_outtake;

			// not removed residues are put into litter
			i.litter_leaf += i.cmass_leaf + i.cmass_agpool - residue_outtake;
		}
		i.cmass_leaf = 0.0;
		i.cmass_agpool = 0.0;

		// Nitrogen:
		if ((i.nmass_leaf + i.nmass_agpool) > 0.0) {

			// removed residues are oxidised
			residue_outtake = pft.res_outtake * (i.nmass_leaf + i.nmass_agpool);
			i.nmass_litter_leaf += i.nmass_leaf + i.nmass_agpool - residue_outtake;

			// not removed residues are put into litter
			i.anflux_harvest += residue_outtake;
		}
		i.nmass_leaf = 0.0;
		i.nmass_agpool = 0.0;
	}
	else if(pft.phenology == ANY) {

		// Intercrop grass
		if(isintercropgrass) {

			// roots

			// all root carbon and nitrogen goes to litter
			if(i.cmass_root > 0.0)
				i.litter_root += i.cmass_root;
			if(i.nmass_root > 0.0)
				i.nmass_litter_root += i.nmass_root;
			if(i.nstore_labile > 0.0)
				i.nmass_litter_root += i.nstore_labile;
			if(i.nstore_longterm > 0.0)
				i.nmass_litter_root += i.nstore_longterm;

			i.cmass_root = 0.0;
			i.nmass_root = 0.0;
			i.nstore_labile = 0.0;
			i.nstore_longterm = 0.0;


			// leaves

			// Carbon:
			if(i.cmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.cmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.litter_leaf += i.cmass_leaf - harvest;

				if(ifslowharvestpool) {
					i.harvested_products_slow += harvest * pft.harvest_slow_frac; // no slow harvest for grass
					harvest = harvest * (1 - pft.harvest_slow_frac);
				}

				i.acflux_harvest += harvest;
			}
			i.cmass_leaf = 0.0;
			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			if(i.nmass_leaf > 0.0) {

				// Harvest/Grazing of leaves:
				harvest = pft.harv_eff_ic * i.nmass_leaf;	// currently no harvest of intercrtop grass

				// not removed grass is put into litter
				i.nmass_litter_leaf += i.nmass_leaf - harvest;

				if(ifslowharvestpool) {
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

			if(ifslowharvestpool) {
				i.harvested_products_slow += harvest * pft.harvest_slow_frac;
				harvest = harvest * (1 - pft.harvest_slow_frac);
			}
			if(alive)
				i.acflux_harvest += harvest;
			i.cmass_leaf -= harvest;

			i.cmass_ho = 0.0;
			i.cmass_agpool = 0.0;

			// Nitrogen:
			// Reduced removal of N relative to C during grazing.
			double N_harvest_scale = 0.25; // Value that works. Needs to be verified in literature.
			harvest = pft.harv_eff * i.nmass_leaf * N_harvest_scale;

			if(ifslowharvestpool) {
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
 *  Function for balancing carbon and nitrogen fluxes from last year's growth if old-style harvest is selected (HARVEST_GRSC defined),
 *  or, alternatively, this years harvested carbon and nitrogen.
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
 *  INPUT PARAMETER
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


	if(alive || istruecrop_or_intercropgrass)  {
		cp.litter_root += cp.cmass_root;

		if(burn) {
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

	if(burn) {
		cp.anflux_harvest += cp.nmass_leaf;
		cp.anflux_harvest += cp.nmass_sap;
		cp.anflux_harvest += cp.nmass_heart;
	}
	else {
		cp.nmass_litter_leaf += cp.nmass_leaf;
		cp.nmass_litter_sap += cp.nmass_sap;
		cp.nmass_litter_heart += cp.nmass_heart;
	}

	if(pft.landcover == CROPLAND) {
		if(pft.aboveground_ho) {
			if(burn) {
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

		if(burn) {
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

/// Transfer of this year's growth (ycmass_xxx) to cmass_xxx_inc
/**   and pasture grass grown in cropland.
 *  INPUT PARAMETERS 
 *  \param cmass_leaf					leaf C biomass (kgC/m2)
 *  \param cmass_root					fine root C biomass (kgC/m2)
 *  \param cmass_ho						harvestable organ C biomass (kgC/m2)
 *  \param cmass_agpool					above-ground pool C biomass (kgC/m2)
 *  OUTPUT PARAMETERS 
 *  \param cmass_leaf_inc				leaf C biomass (kgC/m2)
 *  \param cmass_root_inc				fine root C biomass (kgC/m2)
 *  \param cmass_ho_inc					harvestable organ C biomass (kgC/m2)
 *  \param cmass_agpool_inc 			above-ground pool C biomass (kgC/m2)  
 */ 
void growth_crop_year(double cmass_leaf, double cmass_root, double cmass_ho, double cmass_agpool,
	double& cmass_leaf_inc, double& cmass_root_inc, double& cmass_ho_inc, double& cmass_agpool_inc) {

	// true crop growth and grass intercrop growth; NB: bminit (cmass_repr & cmass_excess subtracted) not used !

	cmass_leaf_inc = cmass_leaf;
	cmass_root_inc = cmass_root;
	cmass_ho_inc = cmass_ho;
	cmass_agpool_inc = cmass_agpool;

	return;
}

void growth_crop_year(Individual& indiv, double& cmass_leaf_inc, double& cmass_root_inc, double& cmass_ho_inc, double& cmass_agpool_inc) {

	// true crop growth and grass intercrop growth; NB: bminit (cmass_repr & cmass_excess subtracted) not used !

	double cmass_leaf = indiv.cropindiv->ycmass_leaf;
	double cmass_root = indiv.cropindiv->ycmass_root;
	double cmass_ho = indiv.cropindiv->ycmass_ho;
	double cmass_agpool = indiv.cropindiv->ycmass_agpool;

	if(indiv.has_daily_turnover()) {

		indiv.cmass_leaf = 0.0;
		indiv.cmass_root = 0.0;
		indiv.cropindiv->cmass_ho = 0.0;
		indiv.cropindiv->cmass_agpool = 0.0;

		// Not completely accurate here when comparing this year's cmass after turnover with cmass increase (ycmass),
		// which could be from the preceding season, but probably OK, since values are not used for C balance.
		if(indiv.continous_grass()) {
			indiv.cmass_leaf = indiv.cmass_leaf_post_turnover;
			indiv.cmass_root = indiv.cmass_root_post_turnover;
		}
	}

	cmass_leaf_inc = indiv.cropindiv->ycmass_leaf;
	cmass_root_inc = indiv.cropindiv->ycmass_root;
	cmass_ho_inc = indiv.cropindiv->ycmass_ho;
	cmass_agpool_inc = indiv.cropindiv->ycmass_agpool;

	return;
}

/// Yield function for true crops and intercrop grass.
void yield_crop(Individual& indiv) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());

	if(indiv.pft.phenology == ANY) {			// grass intercrop yield
	
		// Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
		if(cropindiv.ycmass_leaf > 0.0)									
			cropindiv.yield = cropindiv.ycmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		else
			cropindiv.yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		if(cropindiv.harv_cmass_leaf > 0.0)			
			cropindiv.harv_yield = cropindiv.harv_cmass_leaf * indiv.pft.harv_eff_ic * 2.0;
		else
			cropindiv.harv_yield = 0.0;
	}
	else if(indiv.pft.phenology == CROPGREEN) {		//true crop yield
	
		// Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
		if(cropindiv.ycmass_ho > 0.0)									
			cropindiv.yield = cropindiv.ycmass_ho * indiv.pft.harv_eff * 2.0;
		else
			cropindiv.yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		if(cropindiv.harv_cmass_ho > 0.0)									
			cropindiv.harv_yield=cropindiv.harv_cmass_ho * indiv.pft.harv_eff * 2.0;
		else
			cropindiv.harv_yield = 0.0;

		// Yield dry wieght of actually harvest products this year; NB as above
		for(int i=0;i<2;i++) {
			if(cropindiv.cmass_ho_harvest[i] > 0.0)								
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
	if(cmass_leaf_inc > 0.0)									
		cropindiv.yield = cmass_leaf_inc * indiv.pft.harv_eff * 2.0;
	else
		cropindiv.yield = 0.0;
	cropindiv.harv_yield = cropindiv.yield;	// Although no specified harvest date, harv_yield is set for compatibility.

	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////  End of landcover harvest functions  ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Bondeau A, Smith PC, Zaehle S, Schaphoff S, Lucht W, Cramer W, Gerten D, Lotze-Campen H,
//   Müller C, Reichstein M & Smith B 2007. Modelling the role of agriculture for the 
//   20th century global terrestrial carbon balance. Global Change Biology, 13:679-706.
// Lindeskog M, Arneth A, Bondeau A, Waha K, Seaquist J, Olin S, & Smith B 2013. 
//   Implications of accounting for land use in simulations of ecosystem services and  
//   carbon cycling in Africa. Earth Syst Dynam Discuss 4:235-278.
// Neitsch SL, Arnold JG, Kiniry JR et al.2002 Soil and Water Assessment Tool, Theorethical 
//   Documentation + User's Manual. USDA_ARS-SR Grassland, Soil and Water Research Laboratory.
//   Agricultural Reasearch Service, Temple,Tx, US.
// Waha K, van Bussel LGJ, Müller C, and Bondeau A.2012. Climate-driven simulation of global 
//   crop sowing dates, Global Ecol Biogeogr 21:247-259
