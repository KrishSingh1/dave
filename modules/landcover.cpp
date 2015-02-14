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
#define GRASS_SEED_CMASS	// Carbon allocated to cover-crop grass on bicdate or the day after turnover.
//#define PRINT_GROSS_LC_CHANGE_INFO

/// Autumn sowing types for crops
enum {NOFORCING, AUTUMNSOWING, SPRINGSOWING};

/////////////////////////// Functions facilitating handling time periods spanning newyear //////////////////////////////////////

/// Query whether a date is within a period spanned by two dates.
bool dayinperiod(int day, int start, int end) {

	bool acrossnewyear = false;

	if(day < 0 || start < 0 || end < 0)	// a negative value should not be a valid day
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

	if(day < 0)							// a negative value should not be a valid day
		return -1;
	else if(day + step > 0)
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
void landcover_init(Gridcell& gridcell, LandcoverInputModule* landcover_input_module) {

	// Set CFT-specific members of gridcellpft:
	for(unsigned int p = 0; p < gridcell.pft.nobj; p++) {
		Gridcellpft& gcpft = gridcell.pft[p];

		if (gridcell.get_lat() >= 0.0) {
			gcpft.sdate_default = gcpft.pft.sdatenh;
			gcpft.hlimitdate_default = gcpft.pft.hlimitdatenh;
		}
		else {
			gcpft.sdate_default = gcpft.pft.sdatesh;
			gcpft.hlimitdate_default = gcpft.pft.hlimitdatesh;
		}
		// double cropping in China and Japan.
		if (!strncmp((char*)gcpft.pft.name, "TrRi", strlen("TrRi")) && gridcell.get_lon() >= 60.0 && gridcell.get_lat() <= 30.0)
			gcpft.multicrop = true;
	}

	// get landcover and crop area fractions from landcover input file(s) or ins-file.
	landcover_input_module->getlandcover(gridcell);

	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		st.frac_old = st.frac;

		if(st.frac > 0.0) {
			gridcell.create_stand_lu(st, st.frac);
		}

		stlist.nextobj();
	}
}

/// landcover_change_transfer constructor
landcover_change_transfer::landcover_change_transfer() {

	transfer_litter_leaf = transfer_litter_sap = transfer_litter_heart = transfer_litter_root = transfer_litter_repr = transfer_harvested_products_slow = NULL;
	transfer_nmass_litter_leaf = transfer_nmass_litter_sap = transfer_nmass_litter_heart = transfer_nmass_litter_root = transfer_harvested_products_slow_nmass = NULL;

	transfer_acflux_harvest = transfer_anflux_harvest = transfer_cpool_fast = transfer_cpool_slow = transfer_wcont_evap = transfer_decomp_litter_mean = 0.0;
	transfer_k_soilfast_mean = transfer_k_soilslow_mean = transfer_nmass_avail = transfer_snowpack = transfer_snowpack_nmass = transfer_anfix_calc = 0.0;

	memset(transfer_wcont,0,NSOILLAYER*sizeof(double));

	for(int i=0; i<NSOMPOOL; i++)
		transfer_sompool[i].ntoc = 0.0;

	allocate();

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

int index(int from, int to, int ncols = nst) {

	return from * ncols + to;
}

/// Gets this year's landcover and crop area fractions, checks that area changes are significant and that the net changes are zero.
/** Stores changes in area fractions for the stand types
 *
 *  OUTPUT PARAMETERS
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param LCchangeCtransfer				whether to transfer carbon, nitrogen and water of reduced stands to expanding stands
 */
bool checkLCchange(Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], bool& LCchangeCtransfer, LandcoverInputModule* input_module) {

	double cropfrac_sum_old = 0.0;
	double change_stand = 0.0;
	double changeLC = 0.0;
	double change_crop = 0.0;
	double transferred_fraction = 0.0;
	double receiving_fraction = 0.0;
	bool change = true;


	//Save old fraction values:									
	for(int i=0; i<NLANDCOVERTYPES; i++)
		gridcell.landcoverfrac_old[i] = gridcell.landcoverfrac[i];
	for(int i=0; i<nst; i++) {
		StandType& st = stlist[i];

		st.frac_old = st.frac;
		if(st.landcover == CROPLAND)
			cropfrac_sum_old += st.frac_old;
	}

	//Get new gridcell.landcoverfrac and/or standtype.frac from LUdata and CFTdata.		
	input_module->getlandcover(gridcell);	

	for(int i=0; i<nst; i++) {
		StandType& st = stlist[i];

		st.frac_change = st.frac - st.frac_old ;
	}

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

	if(run[CROPLAND] && (!frac_fixed[CROPLAND] || !lcfrac_fixed)) {
		for(int i=0;i<nst;i++) {

			if(stlist[i].landcover == CROPLAND) {

				double stfrac_change = stlist[i].frac_change;

				if(stfrac_change < 0.0)
					transferred_fraction -= stfrac_change;
				if(stfrac_change > 0.0)
					receiving_fraction += stfrac_change;

				if(cropfrac_sum_old != 0.0) {
					change_crop += fabs(stfrac_change) / 2.0;
				}
				else {
					change_crop += fabs(stfrac_change);
				}
				change_stand += fabs(stfrac_change) / 2.0;
			}
		}
	}

	// if no changes, do nothing.
	if(changeLC < 1.0e-15 && change_crop < 1.0e-15) {
		change = false;
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
		change = true;
	}

	return change;
}

enum {NONEWSTAND, CLONESTAND, CLONESTAND_KILLTREES, NEWSTAND_KILLALL};

/// contains rules for creation of new stands at land cover change
/** Options CLONESTAND and CLONESTAND_KILLTREES require that the receptor landcover
 *  allows growth of natural grass and/or tree PFTs
 *
 *  INTPUT PARAMETERS
 *
 *  \param landcover_donor				landcover type of donor stand
 *  \param landcover_receptor			landcover type of rexeptor stand
 */
int copy_stand_type(int landcover_donor, int landcover_receptor) {

	int copy_type = NONEWSTAND;

	if(landcover_donor == NATURAL) {
		
		if(landcover_receptor == FOREST)
			copy_type = CLONESTAND;

//		if(landcover_receptor == PASTURE)
//			copy_type = CLONESTAND_KILLTREES;

//		if(landcover_receptor == CROPLAND)
//			copy_type = NEWSTAND_KILLALL;

//		if(landcover_receptor == PASTURE)
//			copy_type = NEWSTAND_KILLALL;

	}
	else if(landcover_donor == FOREST) {
		if(landcover_receptor == NATURAL)
			copy_type = CLONESTAND;
	}

	return copy_type;
}

/// identifies which stands to reduce in area and sets standtype.nstands
/** Updates frac, frac_change, frac_old and gross_frac_decrease for reduced stands
 *  Updates frac_old and frac_temp for all stands
 *
 *  INTPUT PARAMETERS
 *
 *  \param st_frac_transfer				array with this year's transferred area fractions between the different stand types
 */
void reduce_stands(Gridcell& gridcell, double* st_frac_transfer, double* primary_st_frac_transfer) {

	for(int i=0; i<nst; i++)
		stlist[i].nstands = 0;

	for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
		Stand& stand = gridcell[i];

		stand.frac_old = stand.get_gridcell_fraction();
		stand.frac_temp = stand.frac_old;
		stand.frac_change = 0.0;
		stand.gross_frac_increase = 0.0;
		stand.gross_frac_decrease = 0.0;
		stand.cloned_fraction = 0.0;
		if(stand.transfer_area_st)
			memset(stand.transfer_area_st, 0, nst * sizeof(double));
		stlist[stand.stid].nstands++;
	}


	for(int from=0; from<nst; from++) {

		StandType& st = stlist[from];

		if(st.gross_frac_decrease > 0.0) {

			if(st.nstands > 1) {

				int nlaps = 1;
				if(ifprimary_lc_transfer)
					nlaps = 2;

				for(int n=0; n<nlaps; n++) {

					for(int to=0; to<nst; to++) {

						if(st_frac_transfer[index(from, to)] > 0.0) {

							double st_change_remain = -st_frac_transfer[index(from, to)];

							bool young_stands_first = true;	//convert area from youngest stands first

							if(st.landcover == NATURAL || st.landcover == FOREST) {

								if(stlist[to].landcover == NATURAL || stlist[to].landcover == FOREST)
									young_stands_first = false;
								else
									young_stands_first = true;
							}

							if(n == 0) {
								if(ifprimary_lc_transfer)
									st_change_remain += primary_st_frac_transfer[index(from, to)];
							}
							else {
								st_change_remain = -primary_st_frac_transfer[index(from, to)];
								young_stands_first = false;
							}

							bool restart = false;

							while(st_change_remain < 0.0) {

								int count_st = 0;

								for(unsigned int i = 0; i < gridcell.size(); i++) {
									int index;

									double stands_frac_sum = 0.0;

									for(unsigned int j = 0; j < gridcell.nbr_stands(); j++) {
										Stand& stand = gridcell[j];

										if(stand.stid == st.id)
											stands_frac_sum += stand.get_gridcell_fraction();
									}

									if(st_change_remain > -1.0e-15 || stands_frac_sum == 0.0) {
										if(st_change_remain <= -1.0e-15 && stands_frac_sum == 0.0)
											dprintf("\nWarning: no more stand area left of stand type %d ! Residual reduction demand %.15f ignored.\n", st.id, st_change_remain);
										st_change_remain = 0.0;
										break;
									}

									if(young_stands_first)
										index = gridcell.size() - 1 - i;
									else
										index = i;

									Stand& stand = gridcell[index];	

									if(stand.stid == st.id) {
										count_st++;

										// Don't reduce stands younger than the age limit, unless this is the last stand in the loop or the initial stand has been killed
	//									if(gross_land_transfer && date.year - stand.first_year < age_limit_reduce && count_st != st.nstands && !reduce_all_stands && !restart) continue;
										if(date.year - stand.first_year < age_limit_reduce && count_st != st.nstands && !reduce_all_stands && !restart) continue;
							
										// convert equal percentage of areas from all stands
										if(reduce_all_stands) {
											stand.frac_change += st_change_remain * stand.get_gridcell_fraction() / st.frac;
											stand.gross_frac_decrease -= st_change_remain * stand.get_gridcell_fraction() / st.frac;
											stand.transfer_area_st[to] = -st_change_remain * stand.get_gridcell_fraction() / st.frac;
											stand.set_gridcell_fraction(stand.get_gridcell_fraction() + st_change_remain * stand.get_gridcell_fraction() / st.frac);
										}
										else {
											if(stand.get_gridcell_fraction() > 0.0) {

												//all natural landcover decrease is taken from this stand
												if(stand.get_gridcell_fraction() >= -st_change_remain) {
													stand.frac_change += st_change_remain;
													stand.gross_frac_decrease -= st_change_remain;
													stand.transfer_area_st[to] -= st_change_remain;
													stand.set_gridcell_fraction(stand.get_gridcell_fraction() + st_change_remain);

													if(stand.get_gridcell_fraction() < 1.0e-15 || st.frac == 0.0 && stand.get_gridcell_fraction() < 1.0e-14)
														stand.set_gridcell_fraction(0.0);

													st_change_remain = 0.0;
													break;
												}
												//more stands will have to be reduced
												else {						
													stand.frac_change -= stand.get_gridcell_fraction();
													stand.gross_frac_decrease += stand.get_gridcell_fraction();
													stand.transfer_area_st[to] += stand.get_gridcell_fraction();
													st_change_remain += stand.get_gridcell_fraction();
													stand.set_gridcell_fraction(0.0);	//will be killed below
												}				
											}
										}
									}
								}

								// Restart loop and reduce oldest stand first if the rules above precluded reduction of the required area.
								if(st_change_remain < 0.0) {
									if(reduce_all_stands) {
										st_change_remain = 0.0;
									}
									else {
										young_stands_first = false;
										restart = true;
									}
								}
							}		
						}
					}
				}
			}
			else if(st.nstands == 1 ) {

				for(unsigned int i = 0; i < gridcell.size(); i++) {

					Stand& stand = gridcell[i];

					if(stand.stid == st.id) {

						stand.frac_change = -st.gross_frac_decrease;
						stand.gross_frac_decrease = st.gross_frac_decrease;
						stand.set_gridcell_fraction(stand.get_gridcell_fraction() + stand.frac_change);

						if(stand.get_gridcell_fraction() < 1.0e-15 || st.frac == 0.0 && stand.get_gridcell_fraction() < 1.0e-14)
							stand.set_gridcell_fraction(0.0);

						for(int to=0; to<nst; to++) {

							if(st_frac_transfer[index(st.id, to)])
								stand.transfer_area_st[to] = st_frac_transfer[index(st.id, to)];
						}
					}
				}
			}
		}
	}
}

/// identifies which stands to expand in area
/** Updates frac, frac_change, frac_old and gross_frac_increase for expanded stands
 * Should be preceded by a call to reduce_stands() and, optionally, transfer_to_new_stand()
 *
 *  INTPUT PARAMETERS
 *
 *  \param st_frac_transfer				array with this year's transferred area fractions between the different stand types
 */
void expand_stands(Gridcell& gridcell, double* st_frac_transfer) {

	// Updated variables are zeroed in reduce_stands()

	for(int i=0; i<nst; i++) {

		StandType& st = stlist[i];
		landcovertype lc = st.landcover;
		bool expand_to_new_stand = ifexpand_to_new_stand && gridcell.expand_to_new_stand[lc];

		if(st.gross_frac_increase > 0.0) {	// Not cloned stands

			if(!expand_to_new_stand) {

				for(unsigned int i = 0; i < gridcell.size(); i++) {

					Stand& stand = gridcell[i];

					if(stand.stid == st.id) {

						// Identify which stands to expand (not secondary stands)
 
						if(!stand.cloned) {

							stand.frac_change += st.gross_frac_increase;
							stand.gross_frac_increase = st.gross_frac_increase;
							stand.set_gridcell_fraction(stand.get_gridcell_fraction() + st.gross_frac_increase);
							if(fabs(stand.frac_change) < 1.0e-15)
								stand.frac_change = 0.0;
						}
					}
				}
			}
		}
	}
}

/// sets land cover transfer matrix from gross land cover change data when no input is available
/** Uses rules to select preferred transfers between land covers
 *
 *  INPUT PARAMETERS
 *  \param landcoverfrac_change			array with this year's difference in area fractions of the different landcovers
 *
 *  OUTPUT PARAMETERS
 *
 *  \param lc_frac_transfer				array with this year's transferred area fractions between the different land covers
 */
void set_lc_change_array(double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES]) {

	const int NRANK = 3;
	int target_preference[NLANDCOVERTYPES][NLANDCOVERTYPES];
	int origin_preference[NLANDCOVERTYPES][NLANDCOVERTYPES];
	double receptor_remain[NLANDCOVERTYPES];
	double donor_remain[NLANDCOVERTYPES];
	int ndonor_lc = 0;
	int nreceptor_lc = 0;

	memset(target_preference, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(int));
	memset(origin_preference, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(int));
	memset(donor_remain, 0, NLANDCOVERTYPES * sizeof(double));
	memset(receptor_remain, 0, NLANDCOVERTYPES * sizeof(double));

	target_preference[CROPLAND][PASTURE] = 2;
	target_preference[CROPLAND][NATURAL] = 3;
	target_preference[CROPLAND][FOREST] = 1;

	target_preference[PASTURE][CROPLAND] = 1;
	target_preference[PASTURE][NATURAL] = 3;
	target_preference[PASTURE][FOREST] = 2;

	target_preference[FOREST][CROPLAND] = 1;
	target_preference[FOREST][PASTURE] = 2;
	target_preference[FOREST][NATURAL] = 3;

	target_preference[NATURAL][CROPLAND] = 1;
	target_preference[NATURAL][PASTURE] = 2;
	target_preference[NATURAL][FOREST] = 3;

	origin_preference[PASTURE][CROPLAND] = 2;
	origin_preference[NATURAL][CROPLAND] = 3;
	origin_preference[FOREST][CROPLAND] = 1;

	origin_preference[CROPLAND][PASTURE] = 2;
	origin_preference[NATURAL][PASTURE] = 3;
	origin_preference[FOREST][PASTURE] = 1;

	origin_preference[CROPLAND][NATURAL] = 2;
	origin_preference[PASTURE][NATURAL] = 3;
	origin_preference[FOREST][NATURAL] = 1;

	origin_preference[CROPLAND][FOREST] = 1;
	origin_preference[PASTURE][FOREST] = 2;
	origin_preference[NATURAL][FOREST] = 3;

	for(int i=0; i<NLANDCOVERTYPES; i++) {

		if(landcoverfrac_change[i] > 0.0) {
			receptor_remain[i] = landcoverfrac_change[i];
			nreceptor_lc++;
		}
		else if(landcoverfrac_change[i] < 0.0){
			donor_remain[i] = -landcoverfrac_change[i];
			ndonor_lc++;
		}
	}


	// Simplest cases: no ambiguities
	if(ndonor_lc == 1 || nreceptor_lc == 1) {

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(landcoverfrac_change[from] < 0.0) {

				for(int to=0; to<NLANDCOVERTYPES; to++) {
				
					if(ndonor_lc == 1) {

						if(landcoverfrac_change[to] > 0.0)
							lc_frac_transfer[from][to] = landcoverfrac_change[to];
					}
					else if(nreceptor_lc == 1) {

						if(landcoverfrac_change[to] > 0.0)
							lc_frac_transfer[from][to] = - landcoverfrac_change[from];
					}
				}
			}
		}
	}
	else {

		for(int score=NRANK*2; score>0; score--) {

			for(int from=0; from<NLANDCOVERTYPES; from++) {

				if(donor_remain[from] > 1.0e-14) {

					for(int to=0; to<NLANDCOVERTYPES; to++) {

						// Identify receiving land covers:	
						if(target_preference[from][to] + origin_preference[from][to] == score && receptor_remain[to] > 1.0e-14 && donor_remain[from] > 1.0e-14) {

							// all donor land is put into this land cover
							if(receptor_remain[to] >= donor_remain[from]) {

								lc_frac_transfer[from][to] = donor_remain[from];
								receptor_remain[to] -= donor_remain[from];
								donor_remain[from] = 0.0;
								break;
							}
							// transfer to more land cover types
							else {

								lc_frac_transfer[from][to] += receptor_remain[to];
								donor_remain[from] -= receptor_remain[to];
								receptor_remain[to] = 0.0;
							}
						}
					}
				}
			}
		}
	}
}

/// sets stand type transfer matrix from land cover transfer matrix when no stand type transfer input is available
/** Distributes land cover transfers equally between stand types within a land cover but may use
 *  rules to select preferred transfers between stand types/land covers (see set_lc_change_array()) if
 *  equal_distribution is set to false.
 *
 *  INPUT PARAMETERS
 *  \param lc_frac_transfer				array with this year's transferred area fractions between the different land covers
 *
 *  OUTPUT PARAMETERS
 *
 *  \param st_frac_transfer				array with this year's transferred area fractions between the different stand types
 */
void set_st_change_array(Gridcell& gridcell, double lc_frac_transfer[][NLANDCOVERTYPES], double* st_frac_transfer, double primary_lc_frac_transfer[][NLANDCOVERTYPES], double* primary_st_frac_transfer) {

	double* net_donor_remain;
	double* net_receptor_remain;
	double* recip_donor_remain;
	double* recip_receptor_remain;
	double recip_lc_change[NLANDCOVERTYPES] = {0.0};
	double recip_lc_frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	double recip_transfer_remain[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	double net_lc_frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	double net_transfer_remain[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	double net_lc_increase[NLANDCOVERTYPES] = {0.0};
	double net_lc_decrease[NLANDCOVERTYPES] = {0.0};
	int ndonor_st = 0;
	int nreceptor_st = 0;
	bool equal_distribution = false;

	net_donor_remain = new double[nst];
	net_receptor_remain = new double[nst];
	memset(net_donor_remain, 0, nst * sizeof(double));
	memset(net_receptor_remain, 0, nst * sizeof(double));
	recip_donor_remain = new double[nst];
	recip_receptor_remain = new double[nst];
	memset(recip_donor_remain, 0, nst * sizeof(double));
	memset(recip_receptor_remain, 0, nst * sizeof(double));


	// Quantify "reciprocal" lc change (gross lc change - net lc change)

	// The special case when donor and recipient lc are the same (primary to secondary land) is treated as a "reciprocal" transfer.

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			recip_lc_frac_transfer[from][to] = min(lc_frac_transfer[from][to], lc_frac_transfer[to][from]);
			recip_transfer_remain[from][to] = recip_lc_frac_transfer[from][to];
			recip_lc_change[from] += recip_lc_frac_transfer[from][to];
			net_lc_frac_transfer[from][to] = lc_frac_transfer[from][to] - recip_lc_frac_transfer[from][to];
			net_transfer_remain[from][to] = net_lc_frac_transfer[from][to];
			net_lc_increase[to] += net_lc_frac_transfer[from][to];
			net_lc_decrease[from] += net_lc_frac_transfer[from][to];
		}
	}

	for(int i=0; i<nst; i++) {

		StandType& st = stlist[i];

		if(st.frac_change > 0.0) {
			net_receptor_remain[i] = st.frac_change;
			nreceptor_st++;
		}
		else if(st.frac_change < 0.0){
			net_donor_remain[i] = -st.frac_change;
			ndonor_st++;
		}

		if(recip_lc_change[st.landcover] > 0.0) {
			recip_receptor_remain[i] = recip_lc_change[st.landcover] * st.frac_old / gridcell.landcoverfrac_old[st.landcover];
			recip_donor_remain[i] = recip_lc_change[st.landcover] * st.frac_old / gridcell.landcoverfrac_old[st.landcover];
		}
	}


	// Net land cover change

	// Simplest cases: no ambiguities
	if(ndonor_st == 1 || nreceptor_st == 1) {

		for(int from=0; from<nst; from++) {

			StandType& st_donor = stlist[from];

			if(st_donor.frac_change < 0.0) {

				for(int to=0; to<nst; to++) {

					StandType& st_receptor = stlist[to];
		
					if(ndonor_st == 1) {

						if(st_receptor.frac_change > 0.0)
							st_frac_transfer[index(from, to)] =  st_receptor.frac_change;
					}
					else if(nreceptor_st == 1) {

						if(st_receptor.frac_change > 0.0)
							st_frac_transfer[index(from, to)] = - st_donor.frac_change;
					}
				}
			}
		}
	}
	else {	// The assignment of receptor stand types when ambiguities exist is arbitrary (following position in stand type list), 
			// so pooling of all receptor stand types within a land cover is recommended in this case: pool_to_all_standtypes[lc_receptor] = true

		for(int i=0; i<NLANDCOVERTYPES; i++) {

			double abs_frac_change_sum = 0.0;
			double frac_change_sum = 0.0;
			double receptor_sum = 0.0;
			double donor_sum = 0.0;
			for(int from=0; from<nst; from++) {

				StandType& st = stlist[from];
				if(st.landcover == i) {
					abs_frac_change_sum += fabs(st.frac_change);
					frac_change_sum += st.frac_change;
					if(st.frac_change < 0.0)
						donor_sum -= st.frac_change;
					else if(st.frac_change > 0.0)
						receptor_sum += st.frac_change;
				}
			}
			if(abs_frac_change_sum != fabs(frac_change_sum)) {

				double intraLCtransferx = min(donor_sum, receptor_sum);
				net_lc_frac_transfer[i][i] += intraLCtransferx;
				net_transfer_remain[i][i] = net_lc_frac_transfer[i][i];
				net_lc_decrease[i] += intraLCtransferx;
				net_lc_increase[i] += intraLCtransferx;
			}

		}

		for(int from=0; from<nst; from++) {

			StandType& st_donor = stlist[from];

			if(net_donor_remain[from] > 1.0e-14) {

				for(int to=0; to<nst; to++) {

					StandType& st_receptor = stlist[to];

					if((net_transfer_remain[st_donor.landcover][st_receptor.landcover]  > 1.0e-14)
						&& net_receptor_remain[to] > 1.0e-14 && net_donor_remain[from] > 1.0e-14) {

						double donor_effective = min(net_donor_remain[from], net_transfer_remain[st_donor.landcover][st_receptor.landcover]);
						double receptor_effective = min(net_receptor_remain[to], net_transfer_remain[st_donor.landcover][st_receptor.landcover]);

						if(equal_distribution) {
							if(net_lc_decrease[st_donor.landcover] && net_lc_increase[st_receptor.landcover])
								st_frac_transfer[index(from, to)] = net_lc_frac_transfer[st_donor.landcover][st_receptor.landcover] * net_donor_remain[from] / net_lc_decrease[st_donor.landcover] * net_receptor_remain[to] / net_lc_increase[st_receptor.landcover];
						}
						else {

							// all donor land going to this land cover is put into this stand type
							if(receptor_effective >= donor_effective) {

								st_frac_transfer[index(from, to)] += donor_effective;
								net_receptor_remain[to] -= donor_effective;
								net_donor_remain[from] -= donor_effective;
								net_transfer_remain[st_donor.landcover][st_receptor.landcover] -= donor_effective;
							}
							// transfer to more stand types within this land cover
							else {

								st_frac_transfer[index(from, to)] += receptor_effective;
								net_donor_remain[from] -= receptor_effective;
								net_receptor_remain[to] -= receptor_effective;
								net_transfer_remain[st_donor.landcover][st_receptor.landcover] -= receptor_effective;
							}
						}
					}
				}
			}
		}

	}

	// Add "reciprocal" lc change (gross-net)
	for(int from=0; from<nst; from++) {

		StandType& st_donor = stlist[from];

		if(recip_donor_remain[from] > 1.0e-14) {

			for(int to=0; to<nst; to++) {

				StandType& st_receptor = stlist[to];

				if(recip_transfer_remain[st_donor.landcover][st_receptor.landcover]  > 1.0e-14 
					&& recip_receptor_remain[to] > 1.0e-14 && recip_donor_remain[from] > 1.0e-14) {

					if(equal_distribution) {
						st_frac_transfer[index(from, to)] += recip_lc_frac_transfer[st_donor.landcover][st_receptor.landcover] * st_donor.frac_old / gridcell.landcoverfrac_old[st_donor.landcover] * st_receptor.frac_old / gridcell.landcoverfrac_old[st_receptor.landcover];
					}
					else {

						double donor_effective = min(recip_donor_remain[from], recip_transfer_remain[st_donor.landcover][st_receptor.landcover]);
						double receptor_effective = min(recip_receptor_remain[to], recip_transfer_remain[st_donor.landcover][st_receptor.landcover]);

						// all donor land is put into this stand type
						if(receptor_effective >= donor_effective) {

							st_frac_transfer[index(from, to)] += donor_effective;
							recip_receptor_remain[to] -= donor_effective;
							recip_donor_remain[from] -= donor_effective;
							recip_transfer_remain[st_donor.landcover][st_receptor.landcover] -= donor_effective;
						}
						// transfer to more stand types
						else {

							st_frac_transfer[index(from, to)] += receptor_effective;
							recip_donor_remain[from] -= receptor_effective;
							recip_receptor_remain[to] -= receptor_effective;
							recip_transfer_remain[st_donor.landcover][st_receptor.landcover] -= receptor_effective;
						}
					}
				}
			}
		}
	}

	// Set the transfer fraction from primary stands within a stand type (simple case: equal primary/secondary ratio for all stand types in a lc->lc transfer)
	if(ifprimary_lc_transfer) {
		for(int from=0; from<nst; from++) {

			StandType& st_donor = stlist[from];

			for(int to=0; to<nst; to++) {

				StandType& st_receptor = stlist[to];

				double primary_frac = 0.0;

				if(lc_frac_transfer[st_donor.landcover][st_receptor.landcover] > 0.0)
					primary_frac = primary_lc_frac_transfer[st_donor.landcover][st_receptor.landcover] / lc_frac_transfer[st_donor.landcover][st_receptor.landcover];

				primary_st_frac_transfer[index(from, to)] = primary_frac * st_frac_transfer[index(from, to)];
			}
		}
	}

	if(net_donor_remain)
		delete[] net_donor_remain;
	if(net_receptor_remain)
		delete[] net_receptor_remain;
	if(recip_donor_remain)
		delete[] recip_donor_remain;
	if(recip_receptor_remain)
		delete[] recip_receptor_remain;
}

/// Handles harvest and turnover of reduced stands at landcover change.
/** Sets LC_updated to true
 *  Stores carbon, nitrogen and water of harvested area in a temporary struct.
 *  Should be followed by a call to stand_dynamics() to kill stands with a new area of 0 
 *
 *  INPUT PARAMETERS
 *
 *  \param receiving_fraction				sum of added area to expanding stands
 *  \param landcover_receptor				restriction of donor stands transferring to specified receptor landcover
 *  \param landcover_donor 					restriction of donor stands according to landcover
 *  \param stid_receptor 					restriction of donor stands transferring to specified receptor stand type
 *  \param stid_donor  						restriction of stands according to stand type
 *  \param standid							restriction of donor stand
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
void donor_stand_change(Gridcell& gridcell, double& receiving_fraction, landcover_change_transfer& to, int landcover_receptor = -1, int landcover_donor = -1, int stid_receptor = -1, int stid_donor = -1, int standid = -1, bool killgrass = true) {

	double ccont_pre_orig_tot = 0.0;
	double dif_tot = 0.0;
	double acflux_landuse_change_tot = 0.0;
	double ncont_pre_orig_tot = 0.0;
	double ndif_tot = 0.0;
	double anflux_landuse_change_tot = 0.0;
	int count = 0;

	for(unsigned int i=0; i<gridcell.nbr_stands(); i++) {

		double scale;

		Stand& stand = gridcell[i];

		double donor_area = 0.0;

		if(stid_receptor >= 0)
			donor_area = stand.transfer_area_st[stid_receptor];
		else if(landcover_receptor >= 0)
			donor_area = stand.transfer_area_lc(landcover_receptor);
		else
			donor_area = stand.gross_frac_decrease;

		if(donor_area > 0.0 && (stand.id == standid || standid < 0) && (stand.landcover == landcover_donor || landcover_donor < 0) && (stand.stid == stid_donor || stid_donor < 0)) {

			double zero_gridcell_luc_cflux = gridcell.acflux_landuse_change;
			double ccont_stand_pre_orig = 0.0;
			double zero_to_ccont_pre = to.ccont();
			double to_ccont_pre = 0.0;
			double zero_gridcell_luc_nflux = gridcell.anflux_landuse_change;
			double ncont_stand_pre_orig = 0.0;
			double zero_to_ncont_pre = to.ncont();
			double to_ncont_pre = 0.0;
			bool single_stand = true;

			if(fabs(donor_area - receiving_fraction) > 1.0e-14)
				single_stand = false;
			else
				donor_area = receiving_fraction;	// Correct for rounding errors: small differences in results
			count++;

			stand.frac_temp -= donor_area;
			scale = donor_area / receiving_fraction / (double)stand.nobj;

			// Add non-living C, N and water of stand to transfer struct:
			to.add_from_stand(stand, scale);

			// C accounting:
			ccont_stand_pre_orig = stand.ccont();
			to_ccont_pre = to.ccont();
			double dif = ccont_stand_pre_orig - stand.ccont(0);
			ncont_stand_pre_orig = stand.ncont();
			to_ncont_pre = to.ncont();
			double ndif = ncont_stand_pre_orig - stand.ncont(0);
if(transfer_mode != 0) {
			// Harvest and turnover of copies of individuals, add to transfer copy:
			stand.firstobj();
			while(stand.isobj) {
				Patch& patch = stand.getobj();
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
if(transfer_mode != 1) {
					// Harvest of transferred areas:
					switch (stand.landcover)
					{
					case CROPLAND:
						harvest_crop(cp, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass); 
						break;
					case PASTURE:
						harvest_pasture(cp, indiv.pft, indiv.alive);
						break;
					case NATURAL:
					case FOREST:
						harvest_wood(cp, indiv.pft, indiv.alive, 1.0, 1.0, 0.95, 0.9);	// frac_cut=1, harv_eff=1, res_outtake_twig=0.95, res_outtake_coarse_root=0.9
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
						cp.nstore_longterm, cp.max_n_storage,
						indiv.alive);


					// In case any vegetation left (eg. cmass_root in pasture or grass in woodland):
					if(killgrass)
						kill_remaining_vegetation(cp, indiv.pft, indiv.alive, indiv.istruecrop_or_intercropgrass(), false);
}
					//Sum added litter C & N:
					to.transfer_litter_leaf[indiv.pft.id] += cp.litter_leaf * scale;
					to.transfer_litter_root[indiv.pft.id] += cp.litter_root * scale;
					to.transfer_litter_sap[indiv.pft.id] += cp.litter_sap * scale;
					to.transfer_litter_heart[indiv.pft.id] += cp.litter_heart * scale;

					to.transfer_nmass_litter_leaf[indiv.pft.id] += cp.nmass_litter_leaf * scale;
					to.transfer_nmass_litter_root[indiv.pft.id] += cp.nmass_litter_root * scale;
					to.transfer_nmass_litter_sap[indiv.pft.id] += cp.nmass_litter_sap * scale;
					to.transfer_nmass_litter_heart[indiv.pft.id] += cp.nmass_litter_heart * scale;

					gridcell.acflux_landuse_change += cp.acflux_harvest * donor_area / (double)stand.nobj;
					gridcell.acflux_landuse_change_lc[stand.landcover] += cp.acflux_harvest * donor_area / (double)stand.nobj;
					gridcell.anflux_landuse_change += cp.anflux_harvest * donor_area / (double)stand.nobj;

//					gridcell.acflux_landuse_change += -cp.debt_excess * donor_area / (double)stand.nobj;

					if(ifslowharvestpool) {
						to.transfer_harvested_products_slow[indiv.pft.id] += cp.harvested_products_slow * scale;
						to.transfer_harvested_products_slow_nmass[indiv.pft.id] += cp.harvested_products_slow_nmass * scale;
					}

					vegetation.nextobj();
				}

				stand.nextobj();
			}
}
if(transfer_mode == 0) {	// No transfer
	//Test mode 1:
	gridcell.acflux_landuse_change += stand.ccont() * donor_area;
	gridcell.acflux_landuse_change_lc[stand.landcover] += stand.ccont() * donor_area;
	gridcell.anflux_landuse_change += stand.ncont() * donor_area;
}
else if(transfer_mode == 1) {	// Transfer of soil only
	//Test mode 2 (no harvest):
	gridcell.acflux_landuse_change += (stand.ccont() - stand.ccont(0)) * donor_area;
	gridcell.acflux_landuse_change_lc[stand.landcover] += (stand.ccont() - stand.ccont(0)) * donor_area;
	gridcell.anflux_landuse_change += (stand.ncont() - stand.ncont(0)) * donor_area;
}

#ifdef PRINT_GROSS_LC_CHANGE_INFO
/*
			// Rounding errors in reduce_stands can cause differences between donor_area and receiving_fraction, when they should be identical
			// (stand.transfer_area_st[st] and st_transfer[from][to]), creating small imbalances here (ca. 1.0e-12)
			dprintf("\nCcont of donor stand %d before = %.15f\n", stand.id, ccont_stand_pre_orig);
			if(single_stand)
				dprintf("Ccont of donor stand copy before harvest = %.15f\n", to_ccont_pre);
			else
				dprintf("Ccont of donor stand copy before harvest = %.15f\n", (to_ccont_pre - zero_to_ccont_pre) / (donor_area / receiving_fraction));
			dprintf("Living C of stand %d = %.15f\n", stand.id, dif);
			if(single_stand)
				dprintf("Ccont of donor stand copy AFTER = %.15f\n", to.ccont());
			else
				dprintf("Ccont of donor stand copy AFTER = %.15f\n", (to.ccont() - zero_to_ccont_pre) / (donor_area / receiving_fraction));
			dprintf("C lost to atmosphere = %.15f\n", (gridcell.acflux_landuse_change - zero_gridcell_luc_cflux) / donor_area);
			if(single_stand) {
				dprintf("C balance = %.15f\n", to.ccont() - ccont_stand_pre_orig + (gridcell.acflux_landuse_change - zero_gridcell_luc_cflux) / donor_area);
			}
			else {
				dprintf("C balance = %.15f\n", (to.ccont() - zero_to_ccont_pre) / (donor_area / receiving_fraction) - ccont_stand_pre_orig + (gridcell.acflux_landuse_change - zero_gridcell_luc_cflux) / donor_area);
			}

			dprintf("\nNcont of donor stand %d before = %.15f\n", stand.id, ncont_stand_pre_orig);
			if(single_stand)
				dprintf("Ncont of donor stand copy before harvest = %.15f\n", to_ncont_pre);
			else
				dprintf("Ncont of donor stand copy before harvest = %.15f\n", (to_ncont_pre - zero_to_ncont_pre) / (donor_area / receiving_fraction));
			dprintf("Living N of stand %d = %.15f\n", stand.id, ndif);
			if(single_stand)
				dprintf("Ncont of donor stand copy AFTER = %.15f\n", to.ncont());
			else
				dprintf("Ncont of donor stand copy AFTER = %.15f\n", (to.ncont() - zero_to_ncont_pre) / (donor_area / receiving_fraction));
			dprintf("N lost to atmosphere = %.15f\n", (gridcell.anflux_landuse_change - zero_gridcell_luc_nflux) / donor_area);
			if(single_stand) {
				dprintf("N balance = %.15f\n", to.ncont() - ncont_stand_pre_orig + (gridcell.anflux_landuse_change - zero_gridcell_luc_nflux) / donor_area);
			}
			else {
				dprintf("N balance = %.15f\n", (to.ncont() - zero_to_ncont_pre) / (donor_area / receiving_fraction) - ncont_stand_pre_orig + (gridcell.anflux_landuse_change - zero_gridcell_luc_nflux) / donor_area);
			}
*/
#endif
			ccont_pre_orig_tot += ccont_stand_pre_orig * donor_area / receiving_fraction;
			dif_tot += dif;
			acflux_landuse_change_tot += (gridcell.acflux_landuse_change - zero_gridcell_luc_cflux) / receiving_fraction;
			ncont_pre_orig_tot += ncont_stand_pre_orig * donor_area / receiving_fraction;
			ndif_tot += ndif;
			anflux_landuse_change_tot += (gridcell.anflux_landuse_change - zero_gridcell_luc_nflux) / receiving_fraction;

		}
	}
#ifdef PRINT_GROSS_LC_CHANGE_INFO
/*
	if(count > 1) {
		dprintf("\nMultiple stands to transfer copy:\n");
		dprintf("Ccont of donor stands before = %.15f\n", ccont_pre_orig_tot);
		dprintf("Living C of stands = %.15f\n", dif_tot);
		dprintf("Ccont of donor transfer copy = %.15f\n", to.ccont());
		dprintf("C lost to atmosphere = %.15f\n", acflux_landuse_change_tot);
		dprintf("C balance = %.15f\n", to.ccont() - ccont_pre_orig_tot + acflux_landuse_change_tot);
	}

	if(count > 1) {
		dprintf("\nMultiple stands to transfer copy:\n");
		dprintf("Ncont of donor stands before = %.15f\n", ncont_pre_orig_tot);
		dprintf("Living N of stands = %.15f\n", ndif_tot);
		dprintf("Ncont of donor transfer copy = %.15f\n", to.ccont());
		dprintf("N lost = %.15f\n", anflux_landuse_change_tot);
		dprintf("N balance = %.15f\n", to.ncont() - ncont_pre_orig_tot + anflux_landuse_change_tot);
	}
*/
#endif
	if(count) {
		if(fabs(to.ccont() - ccont_pre_orig_tot + acflux_landuse_change_tot) > 1.0e-10)
			dprintf("WARNING: C balance in donor_stand_change() = %.15f\n", to.ccont() - ccont_pre_orig_tot + acflux_landuse_change_tot);
		if(fabs(to.ncont() - ncont_pre_orig_tot + anflux_landuse_change_tot) > 1.0e-10)
			dprintf("WARNING: N balance in donor_stand_change() = %.15f\n", to.ncont() - ncont_pre_orig_tot + anflux_landuse_change_tot);
	}

}


/// Creates and kills stands at landcover change.
/** Harvest of reduced stands need to be done before with donor_stand_change().
 *  Should be followed by a call to receiving_stand_change() for transfer of 
 *  carbon, nitrogen and water.
 * 
 */
void stand_dynamics(Gridcell& gridcell) {

	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st=stlist.getobj();
		landcovertype lc = st.landcover;

		bool expand_to_new_stand = ifexpand_to_new_stand && gridcell.expand_to_new_stand[lc];

		if(st.gross_frac_increase || st.gross_frac_decrease || st.frac == 0.0) {
			// first stand created
			if(st.frac_old == 0.0 && st.frac > 0.0) {
				Stand& stand = gridcell.create_stand_lu(st, st.frac);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("\nYear %d, stand %d, st %d created: initial area=%.20f\n\n", date.year, stand.id, stand.stid, st.frac);
#endif
			}
			// last stand killed
			else if(st.frac_old > 0.0 && st.frac == 0.0) {
				Gridcell::iterator gc_itr = gridcell.begin();
				while (gc_itr != gridcell.end()) {
					Stand& stand = *gc_itr;
					if(stand.stid == st.id) {
#ifdef PRINT_GROSS_LC_CHANGE_INFO
						dprintf("\nYear %d, stand %d, st %d killed: age=%d, area=%.20f\n\n", date.year, stand.id, stand.stid, date.year - stand.first_year, stand.frac_old);
#endif
						gc_itr = gridcell.delete_stand(gc_itr);
					}
					else
						++gc_itr;
				}
			}
			else {
				if(expand_to_new_stand) {
					// new stand created from other landcover types (pooled option, not created in transfer_to_new_stand())
					if(st.gross_frac_increase > 0.0) {

						// Fewer patches in secondary stands if npatch_secondarystand < npatch:
						Stand& stand = gridcell.create_stand_lu(st, st.gross_frac_increase, npatch_secondarystand);

#ifdef PRINT_GROSS_LC_CHANGE_INFO
						dprintf("\nYear %d, stand %d, st %d created: initial area=%.20f\n\n", date.year, stand.id, stand.stid, stand.frac_change);
#endif
					}
				}
				// secondary stand killed if all of its area converted to other landcover types
				if(st.nstands > 1 && st.gross_frac_decrease > 0.0) {
					Gridcell::iterator gc_itr = gridcell.begin();
					while (gc_itr != gridcell.end()) {
						Stand& stand = *gc_itr;
						if(stand.stid == st.id && stand.get_gridcell_fraction() == 0) {
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("\nYear %d, stand %d, st %d killed: age=%d, area=%.20f\n\n", date.year, stand.id, stand.stid, date.year - stand.first_year, stand.frac_old);
#endif
							gc_itr = gridcell.delete_stand(gc_itr);
						}
						else
							++gc_itr;
					}
				}
			}
		}

		stlist.nextobj();
	}
}


 /// Transfers litter etc. of reduced stands to expanding stands at landcover change.
/** Transfers carbon, nitrogen and water of harvested area to expanded areas from a temporary struct.
 *  Is additive - can be called several times to the same landcover or standtype, if the parameter donorfrac_rel is used
 *
 *  INPUT PARAMETERS
 *
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
 *	
 *  \param LCchangeCtransfer				whether to transfer carbon, nitrogen and water of reduced stands to expanding stands
 *  \param landcover						restricts receiving stands to a certain landcover, default no restriction
 *  \param stid								restricts receiving stands to a certain standtype, default no restriction
 *  \param donorfrac_rel					relative part of fraction (to a receiving stand) from a particular donor
 *  \param standid							restricts receiving stands to a certain stand, default no restriction
 */
void receiving_stand_change(Gridcell& gridcell, landcover_change_transfer& from, bool LCchangeCtransfer, int landcover = -1, int stid = -1, double donorfrac_rel = 1.0, int standid = -1) {

	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;

		if(stand.gross_frac_increase > 0.0 && (stand.landcover == landcover || landcover < 0) && (stand.stid == stid || stid < 0 && (stand.id == standid || standid < 0))) {

			double old_frac, added_frac, new_frac;
			double ccont_stand_pre = stand.ccont();
			double ccont_stand_post = 0.0;

			old_frac = stand.frac_temp;
			added_frac = donorfrac_rel * stand.gross_frac_increase;
			new_frac = old_frac + added_frac;
			stand.frac_temp += added_frac;
if(transfer_mode != 0)
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

						if(ifslowharvestpool) {
							patchpft.harvested_products_slow = (patchpft.harvested_products_slow * old_frac + from.transfer_harvested_products_slow[i] * added_frac) / new_frac;
							patchpft.harvested_products_slow_nmass = (patchpft.harvested_products_slow_nmass * old_frac + from.transfer_harvested_products_slow_nmass[i] * added_frac) / new_frac;
						}
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

					double aaet_5_cp[NYEARAAET] = {0.0};
					patch.aaet_5.to_array(aaet_5_cp);
					for(unsigned int i=0;i<NYEARAAET;i++)
						patch.aaet_5.add((aaet_5_cp[i] * old_frac + from.transfer_aaet_5[i] * added_frac) / new_frac);

					patch.soil.anfix_calc = (patch.soil.anfix_calc * old_frac + from.transfer_anfix_calc * added_frac) / new_frac;

					// save individual C and N content for use in scale_indiv()
					for(unsigned int i=0; i<patch.vegetation.nobj ;i++) {
						Individual& indiv = patch.vegetation[i];

						indiv.save_cmass_luc();
						indiv.save_nmass_luc();
					}

					stand.nextobj();
				}
				
				ccont_stand_post += stand.ccont(old_frac / new_frac);

				// set scaling factor to be used in growth() for scaling vegetation C and N:
				stand.scale_LC_change = (stand.frac_old - stand.gross_frac_decrease) / stand.get_gridcell_fraction();
				gridcell.LC_updated = true;

if(transfer_mode == 1 || transfer_mode == 2 && scaling_mode == 0) {
	gridcell.acflux_landuse_change -= (stand.ccont() - stand.ccont(0)) * added_frac; // living C
	gridcell.acflux_landuse_change_lc[stand.landcover] -= (stand.ccont() - stand.ccont(0)) * added_frac;
	gridcell.anflux_landuse_change -= (stand.ncont() - stand.ncont(0)) * added_frac; // living N
//	water not balanced in these tests !
}
			}
if(transfer_mode == 0) {
	gridcell.acflux_landuse_change -= stand.ccont() * added_frac;
	gridcell.acflux_landuse_change_lc[stand.landcover] -= stand.ccont() * added_frac;
	gridcell.anflux_landuse_change -= stand.ncont() * added_frac;
//	water not balanced in these tests !
}
#ifdef PRINT_GROSS_LC_CHANGE_INFO
			if(fabs(ccont_stand_post - (ccont_stand_pre * old_frac + from.ccont() * added_frac) / new_frac) > 1.0e-12)
				dprintf("WARNING: C balance in receiving_stand_change() = %.10f\n", ccont_stand_post - (ccont_stand_pre * old_frac + from.ccont() * added_frac) / new_frac);
/*
			dprintf("Ccont of donor copy = %.15f\n", from.ccont());
			dprintf("Ccont of receptor stand %d before = %.15f\n", stand.id, ccont_stand_pre);
			dprintf("Ccont of receptor stand AFTER = %.15f\n", ccont_stand_post);
			dprintf("C balance = %.15f\n", ccont_stand_post - (ccont_stand_pre * old_frac + from.ccont() * added_frac) / new_frac);
			dprintf("\n");
*/
#endif
		}
		++gc_itr;
	}
}

/// Creates unique stands from transfer events if copy_stand_type() returns >= 1 for landcovers combination
/** Either clones donor stand or creates new stand from scratch
 *
 *  INPUT PARAMETERS
 *
 *  \param stid_donor  						donor stand type
 *  \param stid_receptor 					receptor stand type
 */
double transfer_to_new_stand(Gridcell& gridcell, int stid_donor = -1, int stid_receptor = -1) {

	double cloned_area = 0.0;

	for(unsigned int i=0; i<gridcell.nbr_stands(); i++) {

		Stand& stand = gridcell[i];
		double transfer_area = stand.transfer_area_st[stid_receptor];

		if(transfer_area > 0.0 && (stand.stid == stid_donor || stid_donor < 0)) {

			int copy_type = copy_stand_type(stand.landcover, stlist[stid_receptor].landcover);	// NONEWSTAND, CLONESTAND, CLONESTAND_KILLTREES, NEWSTAND_KILLALL

			if(copy_type) {

				if(copy_type == CLONESTAND || copy_type == CLONESTAND_KILLTREES) {

					Stand& new_stand = stand.clone(stlist[stid_receptor], transfer_area);

					if(stand.landcover == NATURAL && (!stlist[stid_receptor].naturalveg || !stlist[stid_receptor].naturalgrass))
						dprintf("WARNING: cloning natural stand without allowing natural pft:s to grow in the new stand. Is this intended ?\n");

#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Year %d: stand %d (st %d) cloned from stand %d (st %d): ccont=%.15f; frac=%f\n", date.year, new_stand.id, new_stand.stid, stand.id, stand.stid, new_stand.ccont(), new_stand.get_gridcell_fraction());
#endif
					landcover_change_transfer transfer;

					if(copy_type == CLONESTAND_KILLTREES) {
						donor_stand_change(gridcell, transfer_area, transfer, -1,-1, stid_receptor, -1, stand.id, false);
						receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand.id);
					}
	
					new_stand.firstobj();
					while(new_stand.isobj) {
						Patch& patch = new_stand.getobj();
						Vegetation& vegetation = patch.vegetation;
						vegetation.firstobj();
						while(vegetation.isobj) {

							Individual& indiv = vegetation.getobj();
							Standpft& standpft = new_stand.pft[indiv.pft.id];

							if(!standpft.active) {
								indiv.kill(true);
								vegetation.killobj();
							}
							else
								vegetation.nextobj();
						}
						new_stand.nextobj();
					}

					new_stand.cloned = true;
					new_stand.gross_frac_increase = 0.0;
					new_stand.frac_change = 0.0;
					new_stand.cloned_fraction = transfer_area;

#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Year %d: cloned stand after harvest: %d ccont=%.15f\n", date.year, new_stand.id, new_stand.ccont());
//					dprintf("Year %d: cloned stand after harvest: %d ncont=%.15f\n", date.year, new_stand.id, new_stand.ncont());
#endif
				}
				else if(copy_type == NEWSTAND_KILLALL) {

					landcover_change_transfer transfer;

					donor_stand_change(gridcell, transfer_area, transfer, -1,-1, stid_receptor, -1, stand.id);
					Stand& new_stand = gridcell.create_stand_lu(stlist[stid_receptor], transfer_area, npatch_secondarystand);
					receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand.id);

#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Year %d: new stand %d st %d after transfer from st %d: ccont=%.15f\n\tfraction=%f\n", date.year, new_stand.id, new_stand.stid, stand.stid, new_stand.ccont(), transfer_area);
//					dprintf("Year %d: new stand %d st %d after transfer from st %d: ncont=%.15f\n\tfraction=%f\n", date.year, new_stand.id, new_stand.stid, stand.stid, new_stand.ncont(), transfer_area);
#endif

					new_stand.cloned = true;
					new_stand.gross_frac_increase = 0.0;
					new_stand.frac_change = 0.0;
					new_stand.cloned_fraction = transfer_area;
				}

				cloned_area += transfer_area;
				stand.cloned_fraction -= transfer_area;
				stlist[stid_donor].gross_frac_decrease -= transfer_area;
				stlist[stid_donor].frac_change += transfer_area;
				stand.gross_frac_decrease -= transfer_area;
				stand.frac_change += transfer_area;
				stand.transfer_area_st[stid_receptor] -= transfer_area;
				stlist[stid_receptor].gross_frac_increase -= transfer_area;
				stlist[stid_receptor].frac_change -= transfer_area;

				if(stlist[stid_donor].gross_frac_decrease < 1.0e-15)
					stlist[stid_donor].gross_frac_decrease = 0.0;
				if(stand.gross_frac_decrease < 1.0e-15)
					stand.gross_frac_decrease = 0.0;
				if(stand.transfer_area_st[stid_receptor] < 1.0e-15)
					stand.transfer_area_st[stid_receptor] = 0.0;
				if(stlist[stid_receptor].gross_frac_increase < 1.0e-15)
					stlist[stid_receptor].gross_frac_increase = 0.0;
			}
		}
	}

	return cloned_area;
}

/// Test function. Prints landcover change and transfers, called from landcover_dynamics()
void print_fractions( Gridcell& gridcell, double landcoverfrac_change[], double lc_frac_transfer[][NLANDCOVERTYPES],double* st_frac_transfer, double primary_lc_frac_transfer[][NLANDCOVERTYPES],double* primary_st_frac_transfer) {

	char landcovernames[NLANDCOVERTYPES][10] = {"URBAN", "CROPLAND", "PASTURE", "FOREST", "NATURAL", "PEATLAND", "BARREN"};

	double gross_landcoverfrac_increase[NLANDCOVERTYPES] = {0.0};
	double gross_landcoverfrac_decrease[NLANDCOVERTYPES] = {0.0};
	double gross_landcoverfrac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};

	for(int from=0; from<nst; from++) {

		for(int to=0; to<nst; to++) {

			if(st_frac_transfer[index(from, to)] > 0.0) {
				gross_landcoverfrac_increase[stlist[to].landcover] += st_frac_transfer[index(from, to)];
				gross_landcoverfrac_decrease[stlist[from].landcover] += st_frac_transfer[index(from, to)];
				gross_landcoverfrac_transfer[stlist[from].landcover][stlist[to].landcover] += st_frac_transfer[index(from, to)];
			}
		}
	}

	dprintf("landcoverfrac_change year %d:\n", date.year);
	dprintf("%10s","");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10s", landcovernames[i]);
	dprintf("\n");
	dprintf("%10s","Net");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10.5f", landcoverfrac_change[i]);
	dprintf("\n%10s","Decrease");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10.5f", gross_landcoverfrac_decrease[i]);
	dprintf("\n%10s","Increase");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10.5f", gross_landcoverfrac_increase[i]);
	dprintf("\n");
	dprintf("transfer_array: \n");
	dprintf("%10s","");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10s", landcovernames[i]);
	dprintf("\n");
	for(int i=0;i<NLANDCOVERTYPES;i++) {
		dprintf("%10s",landcovernames[i]);
		for(int j=0;j<NLANDCOVERTYPES;j++) {
			dprintf("%10.5f", gross_landcoverfrac_transfer[i][j]);
		}
		dprintf("\n");
	}

	dprintf("landcover fractions:\n");
	dprintf("%10s","");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10s", landcovernames[i]);
	dprintf("\n");
	dprintf("%10s","Old");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10.5f", gridcell.landcoverfrac_old[i]);
	dprintf("\n");
	dprintf("%10s","New");
	for(int i=0;i<NLANDCOVERTYPES;i++)
		dprintf("%10.5f", gridcell.landcoverfrac[i]);
	dprintf("\n");

	if(ifprimary_lc_transfer) {
		dprintf("primary transfer_array: \n");
		dprintf("%10s","");
		for(int i=0;i<NLANDCOVERTYPES;i++)
			dprintf("%10s", landcovernames[i]);
		dprintf("\n");
		for(int i=0;i<NLANDCOVERTYPES;i++) {
			dprintf("%10s",landcovernames[i]);
			for(int j=0;j<NLANDCOVERTYPES;j++) {
				dprintf("%10.5f", primary_lc_frac_transfer[i][j]);
			}
			dprintf("\n");
		}
	}
	dprintf("\n");

	dprintf("standtype_change year %d:\n", date.year);
	dprintf("%10s","");
	for(int i=0;i<nst;i++)
		dprintf("%10s", (char*)stlist[i].name);
	dprintf("\n");
	dprintf("%10s","Net");
	for(int i=0;i<nst;i++)
		dprintf("%10.5f", stlist[i].frac_change);
	dprintf("\n%10s","Decrease");
	for(int i=0;i<nst;i++)
		dprintf("%10.5f", stlist[i].gross_frac_decrease);
	dprintf("\n%10s","Increase");
	for(int i=0;i<nst;i++)
		dprintf("%10.5f", stlist[i].gross_frac_increase);
	dprintf("\n");
	dprintf("transfer_array: \n");
	dprintf("%10s","");
	for(int i=0;i<nst;i++)
		dprintf("%10s", (char*)stlist[i].name);
	dprintf("\n");
	for(int i=0;i<nst;i++) {
		dprintf("%10s",(char*)stlist[i].name);
		for(int j=0;j<nst;j++)
			dprintf("%10.5f", st_frac_transfer[index(i, j)]);
		dprintf("\n");
	}

	dprintf("stand type fractions:\n");
	dprintf("%10s","");
	for(int i=0;i<nst;i++)
		dprintf("%10s", (char*)stlist[i].name);
	dprintf("\n");
	dprintf("%10s","Old");
	for(int i=0;i<nst;i++)
		dprintf("%10.5f", stlist[i].frac_old);
	dprintf("\n");
	dprintf("%10s","New");
	for(int i=0;i<nst;i++)
		dprintf("%10.5f", stlist[i].frac);
	dprintf("\n");

	if(ifprimary_lc_transfer) {
		dprintf("primary transfer_array: \n");
		dprintf("%10s","");
		for(int i=0;i<nst;i++)
			dprintf("%10s", (char*)stlist[i].name);
		dprintf("\n");
		for(int i=0;i<nst;i++) {
			dprintf("%10s",(char*)stlist[i].name);
			for(int j=0;j<nst;j++)
				dprintf("%10.5f", primary_st_frac_transfer[index(i, j)]);
			dprintf("\n");
		}
	}
	dprintf("\n");
}



/// Simulates gross landcover change by adding a specified fraction to the lc_frac_transfer array
void simulate_gross_lc_transfer(Gridcell& gridcell, double lc_frac_transfer[][NLANDCOVERTYPES]) {

	// Added gross fraction transfer as percentage of the lesser of two stand types belonging to certain landcover types
	double gross_lc_change_frac[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	gross_lc_change_frac[CROPLAND][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][CROPLAND] = 0.05;
	gross_lc_change_frac[CROPLAND][PASTURE] = 0.05;
	gross_lc_change_frac[PASTURE][CROPLAND] = 0.05;
	gross_lc_change_frac[PASTURE][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][PASTURE] = 0.05;
	gross_lc_change_frac[FOREST][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][FOREST] = 0.05;
	gross_lc_change_frac[FOREST][PASTURE] = 0.05;
	gross_lc_change_frac[PASTURE][FOREST] = 0.05;

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			lc_frac_transfer[from][to] += gross_lc_change_frac[from][to] * min(gridcell.landcoverfrac_old[from], gridcell.landcoverfrac_old[to]);
		}
	}
}


/// Simulates gross landcover change by adding a specified fraction to the st_frac_transfer array
void simulate_gross_st_transfer(double* st_frac_transfer) {

	// Added gross fraction transfer as percentage of the lesser of two stand types belonging to certain landcover types
	double gross_lc_change_frac[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};
	gross_lc_change_frac[CROPLAND][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][CROPLAND] = 0.05;
	gross_lc_change_frac[CROPLAND][PASTURE] = 0.05;
	gross_lc_change_frac[PASTURE][CROPLAND] = 0.05;
	gross_lc_change_frac[PASTURE][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][PASTURE] = 0.05;
	gross_lc_change_frac[FOREST][NATURAL] = 0.05;
	gross_lc_change_frac[NATURAL][FOREST] = 0.05;
	gross_lc_change_frac[FOREST][PASTURE] = 0.05;
	gross_lc_change_frac[PASTURE][FOREST] = 0.05;

	for(int from=0; from<nst; from++) {

		StandType& st_donor = stlist[from];

		for(int to=0; to<nst; to++) {

			StandType& st_receptor = stlist[to];

			if(gross_lc_change_frac[st_donor.landcover][st_receptor.landcover]) {

				// All stand types with an area:
				st_frac_transfer[index(from, to)] += gross_lc_change_frac[st_donor.landcover][st_receptor.landcover] * min(min(st_donor.frac_old, st_donor.frac), min(st_receptor.frac_old, st_receptor.frac));
			}
		}
	}
}

/// Called after update of st fraction and transfer values
bool check_fractions(Gridcell& gridcell, double landcoverfrac_change[], double lc_change_array[][NLANDCOVERTYPES], double* st_change_array, bool check_lc_st_transfer = false) {

	bool error = false;

	// Check that landcover sum is 1:
	double lc_frac_sum = 0.0;
	for(int i=0; i<NLANDCOVERTYPES; i++)
		lc_frac_sum += gridcell.landcoverfrac[i];

	if(fabs(lc_frac_sum - 1.0)  > 1.0e-14) {
		dprintf("\nCheck 1: Year %d: landcover fraction sum: %.15f", date.year, lc_frac_sum);
		error = true;
	}
	/////////

	// Check that stand type sum is 1:
	double st_frac_sum = 0.0;
	for(int i=0; i<nst; i++)
		st_frac_sum += stlist[i].frac;
	if(fabs(st_frac_sum - 1.0)  > 1.0e-14) {
		dprintf("\nCheck 2: Year %d: stand type fraction sum: %.15f", date.year, st_frac_sum);
		error = true;
	}
	/////////

	// Test if the sum of gross lcc for a landcover is the same as net lcc:
	double test_lc_change[NLANDCOVERTYPES] = {0.0};

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			test_lc_change[from] -= lc_change_array[from][to];
			test_lc_change[to] += lc_change_array[from][to];
		}
	}

	for(int i=0; i<NLANDCOVERTYPES; i++) {

		if(fabs(test_lc_change[i] - landcoverfrac_change[i]) > 1.0e-14) {

			dprintf("\nCheck 3: Year %d: lc_change_array sum not equal to landcoverfrac_change value for landcover %d\n", date.year, i);
			dprintf("dif=%.15f", test_lc_change[i] - landcoverfrac_change[i]);
			error = true;
		}
	}
	////////////

	// Test if the sum of gross lcc for a stand type is the same as net lcc:
	double *test_st_change;
	test_st_change = new double[nst];
	memset(test_st_change,0,nst * sizeof(double));

	for(int from=0; from<nst; from++) {

		for(int to=0; to<nst; to++) {

			test_st_change[from] -= st_change_array[index(from, to)];
			test_st_change[to] += st_change_array[index(from, to)];
		}
	}

	for(int i=0; i<nst; i++) {

		StandType& st = stlist[i];

		if(fabs(test_st_change[i] - st.frac_change) > 1.0e-14) {
			dprintf("\nCheck 4: Year %d: st_change_array sum not equal to st.frac_change value for stand type %d\n", date.year, i);
			dprintf("dif=%.15f", fabs(test_st_change[i] - st.frac_change));
			error = true;
		}
	}
	if(test_st_change)
		delete[] test_st_change;
	////////////

	// Test if st_change_array is compatible with st.frac_old/frac values
	for(int from=0; from<nst; from++) {

		StandType& st_donor = stlist[from];

		for(int to=0; to<nst; to++) {

			StandType& st_receptor = stlist[to];
			if(st_change_array[index(from, to)] > (st_donor.frac_old + 1.0e-14) || st_change_array[index(from, to)] > (st_receptor.frac + 1.0e-14)) {
				dprintf("\nCheck 5: Year %d: st_change_array sum not compatible with st.frac_old/frac values for stand types %d and %d\n", date.year, from, to);
				dprintf("st_change_array=%.15f, st_donor.frac_old=%.15f, st_receptor.frac=%.15f", st_change_array[index(from, to)], st_donor.frac_old, st_receptor.frac);
				error = true;
			}
		}
	}

	// Test if stand type and landcover transfer matrices match each other:
	if(check_lc_st_transfer) {
		double lc_change[NLANDCOVERTYPES] = {0.0};
		double lc_change_arr[NLANDCOVERTYPES][NLANDCOVERTYPES];
		memset(lc_change_arr, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(double));

		for(int from=0; from<nst; from++) {

			StandType& st = stlist[from];

			for(int to=0; to<nst; to++) {

				StandType& st_dest = stlist[to];

				if(st.landcover != st_dest.landcover) {
					lc_change[st.landcover] -= st_change_array[index(from, to)];
					lc_change[st_dest.landcover] += st_change_array[index(from, to)];
					lc_change_arr[st.landcover][st_dest.landcover] += st_change_array[index(from, to)];
				}
			}
		}

		for(int i=0; i<NLANDCOVERTYPES; i++) {
		
			if(fabs(lc_change[i] - landcoverfrac_change[i]) > 1.0e-14) {
				dprintf("\nCheck 6: Year %d: st_change_array LC sum not equal to LC value for %d\n", date.year, i);
				dprintf("dif=%.15f", fabs(lc_change[i] - landcoverfrac_change[i]));
				error = true;
			}
		}

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if(fabs(lc_change_arr[from][to] - lc_change_array[from][to]) > 1.0e-14) {
					dprintf("\nCheck 7: Year %d: lc_change_arr sum not equal to lc_change_array value for %d, %d\n", date.year, from, to);
					dprintf("dif=%.15f", fabs(lc_change_arr[from][to] - lc_change_array[from][to]));
					error = true;
				}
			}
		}
	}
	if(error)
		dprintf("\n\n");

	return error;
}

/// Called before updating stand.frac (reduce_stands)
bool check_fractions1(Gridcell& gridcell) {

	bool error = false;

	// Test if the stand type reduction is smaller than the remaining stand sum:
	for(int s=0; s<nst; s++) {

		StandType& st = stlist[s];

		if(st.frac_change < 0.0){

			double stands_frac_sum = 0.0;

			for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
				Stand& stand = gridcell[i];

				if(stand.stid == st.id)
					stands_frac_sum += stand.get_gridcell_fraction();
			}

			if(-st.frac_change - stands_frac_sum > 1.0e-14) {
				dprintf("\nCheck 8: Year %d: stand type %d fraction reduction bigger than sum of stands\n", date.year, s);
				dprintf("dif=%.15f", -st.frac_change - stands_frac_sum);
				error = true;
			}
		}
	}
	if(error)
		dprintf("\n\n");

	return error;
}

/// Called after updating stand.frac and transfer values (reduce_stands)
bool check_fractions2(Gridcell& gridcell, double* st_change_array) {

	bool error = false;

	// Test if stand.transfer_area_st[to] sum for a stand type is equal to st_change_array[from][to)]
	for(int from=0; from<nst; from++) {

		StandType& st = stlist[from];

		for(int to=0; to<nst; to++) {

			double *test_st_change;
			test_st_change = new double[nst];
			memset(test_st_change,0,nst * sizeof(double));
			
			for(unsigned int i=0; i<gridcell.nbr_stands(); i++) {

				Stand& stand = gridcell[i];

				if(stand.stid == st.id)
					test_st_change[to] += stand.transfer_area_st[to];
			}

			if(fabs(test_st_change[to] - st_change_array[index(from, to)]) > 1.0e-14) {
				dprintf("\nCheck 9: Year %d: stand transfer area sum not equal to stand type value for stand types %d and %d\n", date.year, from, to);
				dprintf("dif=%.15f", fabs(test_st_change[to] - st_change_array[index(from, to)]));
				error = true;
			}

			if(test_st_change)
				delete[] test_st_change;
		}
	}
	////////////

	// Test if stand.frac_change for a stand is equal to stand.gross_frac_increase + stand.gross_frac_decrease:
	for(unsigned int i=0; i<gridcell.nbr_stands(); i++) {

		Stand& stand = gridcell[i];

		if(fabs(stand.frac_change - (stand.gross_frac_increase - stand.gross_frac_decrease)) > 1.0e-14) {
			dprintf("\nCheck 10: Year %d: frac_change is not equal to gross_frac_increase + gross_frac_decrease for stand %d\n", date.year, stand.id);
			dprintf("dif=%.15f\n", fabs(stand.frac_change - (stand.gross_frac_increase - stand.gross_frac_decrease)));
			dprintf("frac_change=%.15f, gross_frac_increase=%.15f, gross_frac_decrease=%.15f", stand.frac_change, stand.gross_frac_increase, stand.gross_frac_decrease);
			error = true;
		}
	}
	////////////

	// Test that no stands remain when stand type fraction is 0 (after reduce_stands, so stands may remain, but with frac 0):
	for(int s=0; s<nst; s++) {

		StandType& st = stlist[s];

		if(st.frac == 0.0) {

			double stands_frac_sum = 0.0;

			for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
				Stand& stand = gridcell[i];

				if(stand.stid == st.id)
					stands_frac_sum += stand.get_gridcell_fraction();
			}

			if(stands_frac_sum) {
				dprintf("\nCheck 11: Year %d: remaining stand when stand type %d fraction is 0\n", date.year, s);
				dprintf("remain=%.20f", stands_frac_sum);
				error = true;
			}
		}
	}

	if(error)
		dprintf("\n\n");

	return error;
}

/// Called after updating stand.frac and transfer values of reduced stands (reduce_stands)
bool check_fractions3(Gridcell& gridcell) {

	bool error = false;

	// Test if sum of stands not equal to stand type value for stand type for reduced stands
	for(int s=0; s<nst; s++) {

		StandType& st = stlist[s];

		double stands_frac_sum = 0.0;

		for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
			Stand& stand = gridcell[i];

			if(stand.stid == st.id)
				stands_frac_sum += stand.get_gridcell_fraction();
		}

		if(st.frac_change < 0.0) {
			if(fabs(st.frac - stands_frac_sum - st.gross_frac_increase) > 1.0e-14) {
				dprintf("\nCheck 12: Year %d: fraction sum of stands not equal to stand type value for stand type %d\n", date.year, s);
				dprintf("dif=%.15f", fabs(st.frac - stands_frac_sum - st.gross_frac_increase));
				error = true;
			}
		}
	}
	if(error)
		dprintf("\n\n");

	return error;
}

/// Called after creating new stands (stand_dynamics)
bool check_fractions4(Gridcell& gridcell) {

	bool error = false;

	// Test if sum of stands not equal to stand type value for stand type for increased or new stands
	for(int s=0; s<nst; s++) {

		StandType& st = stlist[s];

		double stands_frac_sum = 0.0;

		for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
			Stand& stand = gridcell[i];

			if(stand.stid == st.id)
				stands_frac_sum += stand.get_gridcell_fraction();
		}

		if(st.frac_change >= 0.0)
		if(fabs(st.frac - stands_frac_sum) > 1.0e-14) {
			dprintf("\nCheck 13: Year %d: fraction sum of stands not equal to stand type value for stand type %d\n", date.year, s);
			dprintf("dif=%.15f", fabs(st.frac - stands_frac_sum));
			error = true;
		}
	}
	if(error)
		dprintf("\n\n");

	return error;
}

/// Updates all landcover and crop stand area fractions each year, possibly resulting in the creation and killing of stands.
/** Harvests transferred areas and transfers litter etc. of reduced stands to expanding stands.
 *  Transfers litter etc. of reduced stands to expanding stands at landcover change.
 */
void landcover_dynamics(Gridcell& gridcell, LandcoverInputModule* landcover_input_module) {

	double landcoverfrac_change[NLANDCOVERTYPES];
	double lc_frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES];
	double primary_lc_frac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES];
	double* st_frac_transfer = NULL;
	double* primary_st_frac_transfer = NULL;
	bool LCchangeCtransfer = true;

	st_frac_transfer = new double[nst * nst];
	primary_st_frac_transfer = new double[nst * nst];

	memset(landcoverfrac_change, 0, NLANDCOVERTYPES * sizeof(double));
	memset(lc_frac_transfer, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(double));
	memset(primary_lc_frac_transfer, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(double));
	memset(st_frac_transfer, 0, nst * nst * sizeof(double));
	memset(primary_st_frac_transfer, 0, nst * nst * sizeof(double));

	gridcell.LC_updated = false;

	for(unsigned int i=0; i<gridcell.nbr_stands(); ++i)
		gridcell[i].scale_LC_change = 1.0;

	bool no_changes = true;

	// get new landcover and stand type area fractions from input files, set standtype frac_change
	if(!all_fracs_const) {
		// this call returns 0, causing this function to return, if no significant landcover changes this year, 
		// sets LCchangeCtransfer to 0 if unbalanced landcover changes (if some landcovers are inactivated), thus inactivating transfer of C and N
		if(checkLCchange(gridcell, landcoverfrac_change, LCchangeCtransfer, landcover_input_module))
			no_changes = false;
	}

	if(no_changes && !gross_land_transfer) {
		delete[] st_frac_transfer;
		delete[] primary_st_frac_transfer;
		return;
	}

	if(gross_land_transfer == 3) {

		// Read stand type transfer fractions from file here (or in checkLCchange) and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.

		// input_module->get_st_transfer();
		dprintf("Currently no code for option gross_land_transfer==3\n");
	}
	else if(gross_land_transfer == 2) {

		// Read landcover transfer fractions from file here (or in checkLCchange) and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.

		if(landcover_input_module->get_lc_transfer(gridcell, landcoverfrac_change, lc_frac_transfer, primary_lc_frac_transfer)) {
			no_changes = false;
			set_st_change_array(gridcell, lc_frac_transfer, st_frac_transfer, primary_lc_frac_transfer, primary_st_frac_transfer);
		}
	}
	else {

		const bool simulate_st = true;	// gcc simulation at land cover level (false) or stand type level (true)

		set_lc_change_array(landcoverfrac_change, lc_frac_transfer); // the lc_frac_transfer-array is only used in set_st_change_array()

		if(gross_land_transfer && !simulate_st)
			simulate_gross_lc_transfer(gridcell, lc_frac_transfer);

		set_st_change_array(gridcell, lc_frac_transfer, st_frac_transfer, primary_lc_frac_transfer, primary_st_frac_transfer);

		check_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer, true);

		if(gross_land_transfer && simulate_st)
			simulate_gross_st_transfer(st_frac_transfer);
	}

	if(no_changes) {
		delete[] st_frac_transfer;
		delete[] primary_st_frac_transfer;
		return;
	}

	check_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer);
	check_fractions1(gridcell);

	for(int i=0; i<nst; i++) {

		stlist[i].gross_frac_increase = 0.0;
		stlist[i].gross_frac_decrease = 0.0;
	}

	for(int from=0; from<nst; from++) {

		for(int to=0; to<nst; to++) {

			if(st_frac_transfer[index(from, to)] > 0.0) {
				stlist[to].gross_frac_increase += st_frac_transfer[index(from, to)];
				stlist[from].gross_frac_decrease += st_frac_transfer[index(from, to)];
			}
		}
	}

#ifdef PRINT_GROSS_LC_CHANGE_INFO
	print_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer, primary_lc_frac_transfer, primary_st_frac_transfer);
#endif
	double ccont_tot_1 = gridcell.ccont();
	double cflux_tot_1 = gridcell.cflux();
	double ncont_tot_1 = gridcell.ncont();
	double nflux_tot_1 = gridcell.nflux();
#ifdef PRINT_GROSS_LC_CHANGE_INFO
	dprintf("\nYear %d: ccont_tot before luc=%.15f\n", date.year, ccont_tot_1);
	dprintf("Year %d: ncont_tot before luc=%.15f\n\n", date.year, ncont_tot_1);
#endif
	// check how many stands of each stand type exist
	// identify which stands to reduce in area
	// set stand variables frac, frac_old, frac_change and gross_frac_decrease for reduced stands and stand types
	reduce_stands(gridcell, st_frac_transfer, primary_st_frac_transfer);

	int error = 0;
	error += check_fractions2(gridcell, st_frac_transfer);
	error += check_fractions3(gridcell);
	if(error)
//		dprintf("Fraction error after reduce_stands()\n\n");
		fail("Fraction error after reduce_stands()\n\n");

	// Create new stands for land cover transitions when natural vegetation remains or when each transition requires a new stand:
	if(iftransfer_to_new_stand) {
		bool new_stand = false;
		for(int from=0; from<nst; from++) {
			for(int to=0; to<nst; to++) {
				if(st_frac_transfer[index(from, to)] > 0.0) {
					double new_stand_frac = 0.0;
					new_stand_frac = transfer_to_new_stand(gridcell, from, to);
					if(new_stand_frac) {
						st_frac_transfer[index(from, to)] -= new_stand_frac;
						if(st_frac_transfer[index(from, to)] < 1.0e-15)
							st_frac_transfer[index(from, to)] = 0.0;
						new_stand = true;
					}
				}
			}
		}

		if(new_stand) {
			error = 0;
			error += check_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer);
			error += check_fractions2(gridcell, st_frac_transfer);
			if(error)
				dprintf("Fraction error after transfer_to_new_stand()\n\n");
		}
	}

	// set stand variables frac, frac_change and gross_frac_increase for expanding stands and stand types
	expand_stands(gridcell, st_frac_transfer);

	error += check_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer);
	error += check_fractions2(gridcell, st_frac_transfer);
	if(error)
		fail("Fraction error after expand_stands()\n");


	// Pooling options
	// transfer_level: 0: one big pool; 1: land cover-level; 2: stand type-level

	if(transfer_level == 0) {	// One big pool for all transfers (as in old code)

		landcover_change_transfer transfer;
		double receiving_fraction = 0.0;

		for(int from=0; from<nst; from++) {

			for(int to=0; to<nst; to++) {

				if(st_frac_transfer[index(from, to)] > 0.0)
					receiving_fraction += st_frac_transfer[index(from, to)];
			}
		}

		 // handle harvest and turnover of reduced stands at landcover change
		donor_stand_change(gridcell, receiving_fraction, transfer);

		// create and kill stands at landcover change
		stand_dynamics(gridcell);
		error += check_fractions4(gridcell);
		if(error)
			fail("Fraction error after expand_stands()\n");

		// transfer litter etc. of reduced stands to expanding stands at landcover change
		receiving_stand_change(gridcell, transfer, LCchangeCtransfer);

#ifdef PRINT_GROSS_LC_CHANGE_INFO
		dprintf("Year %d: C cont in pooled transfer=%.15f\n", date.year, transfer.ccont());
		dprintf("\treceiving_fraction=%f\n", receiving_fraction);
#endif
	}
	else if(transfer_level == 1) {	// Landcover-level pools, larger pools available by the gridcell variables pool_from_all_landcovers and pool_to_all_landcovers (set in gridcell constructor)

		landcover_change_transfer transfer_lc_2d[NLANDCOVERTYPES][NLANDCOVERTYPES];
		landcover_change_transfer transfer_lc[NLANDCOVERTYPES];
		landcover_change_transfer transfer_lc_from[NLANDCOVERTYPES];
		double gross_landcoverfrac_increase[NLANDCOVERTYPES] = {0.0};
		double gross_landcoverfrac_decrease[NLANDCOVERTYPES] = {0.0};
		double gross_landcoverfrac_transfer[NLANDCOVERTYPES][NLANDCOVERTYPES] = {0.0};

		for(int from=0; from<nst; from++) {

			for(int to=0; to<nst; to++) {

				if(st_frac_transfer[index(from, to)] > 0.0) {

					gross_landcoverfrac_increase[stlist[to].landcover] += st_frac_transfer[index(from, to)];
					gross_landcoverfrac_decrease[stlist[from].landcover] += st_frac_transfer[index(from, to)];
					gross_landcoverfrac_transfer[stlist[from].landcover][stlist[to].landcover] += st_frac_transfer[index(from, to)];
				}
			}
		}

		// handle harvest and turnover of reduced stands at landcover change

		// pool all transferred area from one landcover:
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(gridcell.pool_to_all_landcovers[from]) { // alt.c
				if(gross_landcoverfrac_decrease[from] > 0.0) {

					donor_stand_change(gridcell, gross_landcoverfrac_decrease[from], transfer_lc_from[from], -1, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Year %d after donor_stand_change(): C cont in transfer from landcover %d=%.15f\n", date.year, from, transfer_lc_from[from].ccont());//.15f
					dprintf("\treceiving_fraction=%f\n", gross_landcoverfrac_decrease[from]);
#endif
				}
			}
		}

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			double receiving_fraction = gross_landcoverfrac_increase[to];

			for(int from=0; from<NLANDCOVERTYPES; from++) {

				double receiving_fraction_2d = gross_landcoverfrac_transfer[from][to];

				if(receiving_fraction_2d > 0.0) {

					if(gridcell.pool_from_all_landcovers[to]) {
						if(!gridcell.pool_to_all_landcovers[from]) {	// alt.a
							// Add from->to transfer to to-pool:
							donor_stand_change(gridcell, receiving_fraction, transfer_lc[to], to, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after donor_stand_change(): C cont in transfer from landcover %d to %d=%.15f\n", date.year, from, to, transfer_lc[to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction);
#endif
						}
						else {	// alt.a+c
							// Add from-pool value to to-pool:
							transfer_lc[to].add(transfer_lc_from[from], receiving_fraction_2d / receiving_fraction);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after transfer from donor pool: C cont in transfer from landcover %d to %d=%.15f\n", date.year, from, to, transfer_lc[to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction);
#endif
						}
					}
					else {
						if(!gridcell.pool_to_all_landcovers[from]) {	// alt.b
							// Add from->to transfer to [from][to] place in 2d-array:
							donor_stand_change(gridcell, receiving_fraction_2d, transfer_lc_2d[from][to], to, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after donor_stand_change(): C cont in transfer from landcover %d to %d=%.15f\n", date.year, from, to, transfer_lc_2d[from][to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction_2d);
#endif
						}
						else {	// alt.c
							// Copy from-pool value to [from][to] place in 2d-array:
							transfer_lc_2d[from][to].copy(transfer_lc_from[from]);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after transfer from donor pool: C cont in transfer from landcover %d to %d=%.15f\n", date.year, from, to, transfer_lc_2d[from][to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction_2d);
#endif
						}
					}
				}	
			}
		}

		// create and kill stands at landcover change
		stand_dynamics(gridcell);
		error += check_fractions4(gridcell);
		if(error)
			fail("Fraction error after expand_stands()\n");

		// transfer litter etc. of reduced stands to expanding stands at landcover change
		for(int to=0; to<NLANDCOVERTYPES; to++) {

			if(gridcell.pool_from_all_landcovers[to]) {
				// Alt.a: All donor lc:s are pooled into receiving lc:s
				if(gross_landcoverfrac_increase[to] > 0.0) {
#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Before receiving_stand_change(): C cont in transfer to landcover %d=%.15f\n", to, transfer_lc[to].ccont());//.15f
#endif
					receiving_stand_change(gridcell, transfer_lc[to], LCchangeCtransfer, to);
				}
			}
			else {

				for(int from=0; from<NLANDCOVERTYPES; from++) {

					if(gross_landcoverfrac_transfer[from][to] > 0.0) {
						if(!gridcell.pool_to_all_landcovers[from]) {
						// Alt.b: Transfers from donors are independent:
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Before receiving_stand_change(): C cont in transfer from landcover %d to %d=%.15f\n", from, to, transfer_lc_2d[from][to].ccont());//.15f
#endif
							receiving_stand_change(gridcell, transfer_lc_2d[from][to], LCchangeCtransfer, to, -1, gross_landcoverfrac_transfer[from][to] / gross_landcoverfrac_increase[to]);
						}
						else {
						// Alt.c: Donor lc:s are pooled before transfer to recipients:
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Before receiving_stand_change(): C cont in transfer from landcover %d to %d=%.15f\n", from, to, transfer_lc_from[from].ccont());//.15f
#endif
							receiving_stand_change(gridcell, transfer_lc_from[from], LCchangeCtransfer, to, -1, gross_landcoverfrac_transfer[from][to] / gross_landcoverfrac_increase[to]);
						}
					}
				}
			}
		}
	}
	else if(transfer_level == 2) { // Unique transfers between stand types. Stand type-level pools available by the gridcell variables pool_from_all_landcovers and pool_to_all_landcovers (set in gridcell constructor)

		landcover_change_transfer* transfer_st_2d;
		landcover_change_transfer* transfer_st;
		landcover_change_transfer* transfer_st_from;

		transfer_st_2d = new landcover_change_transfer[nst * nst];
		transfer_st = new landcover_change_transfer[nst];
		transfer_st_from = new landcover_change_transfer[nst];

		// handle harvest and turnover of reduced stands at landcover change

		// pool all transferred area from one landcover:
		for(int from=0; from<nst; from++) {

			StandType& st = stlist[from];

			// Pooling of donor stands within a stand type
//			if(gridcell.pool_to_all_landcovers[st.landcover]) { // alt.c
			if(gridcell.pool_to_all_landcovers[st.landcover] || st.nstands == 1) { // alt.c
				if(st.gross_frac_decrease > 0.0) {

					donor_stand_change(gridcell, st.gross_frac_decrease, transfer_st_from[from], -1, -1, -1, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Year %d after donor_stand_change(): C cont in transfer pool from st %d=%.15f\n", date.year, from, transfer_st_from[from].ccont());//.15f
					dprintf("\treceiving_fraction=%f\n", st.gross_frac_decrease);
#endif
				}
			}
		}

		for(int to=0;to<nst;to++) {

			StandType& st_to = stlist[to];
			double receiving_fraction = st_to.gross_frac_increase;

			for(int from=0; from<nst; from++) {

				StandType& st_from = stlist[from];
				double receiving_fraction_2d = st_frac_transfer[index(from, to)];

				if(receiving_fraction_2d > 0.0) {

					if(gridcell.pool_from_all_landcovers[st_to.landcover]) {
//						if(!gridcell.pool_to_all_landcovers[st_from.landcover]) {	// alt.a
						if(!(gridcell.pool_to_all_landcovers[st_from.landcover] || st_from.nstands == 1)) {	// alt.a
							// Add from->to transfer to to-pool:
							donor_stand_change(gridcell, receiving_fraction, transfer_st[to], -1, -1, to, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after donor_stand_change(): C cont in transfer from st %d to %d=%.15f\n", date.year, from, to, transfer_st[to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction);
#endif
						}
						else {	// alt.a+c
							// Add from-pool value to to-pool:
							transfer_st[to].add(transfer_st_from[from], receiving_fraction_2d / receiving_fraction);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after transfer from donor pool: C cont in transfer from st %d to %d=%.15f\n", date.year, from, to, transfer_st[to].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction);
#endif
						}
					}
					else {
//						if(!gridcell.pool_to_all_landcovers[st_from.landcover]) {	// alt.b
						if(!(gridcell.pool_to_all_landcovers[st_from.landcover] || st_from.nstands == 1)) {	// alt.b
							// Add from->to transfer to [from][to] place in 2d-array:
							donor_stand_change(gridcell, receiving_fraction_2d, transfer_st_2d[index(from, to)], -1, -1, to, from);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after donor_stand_change(): C cont in transfer from st %d to %d=%.15f\n", date.year, from, to, transfer_st_2d[index(from, to)].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction_2d);
#endif
						}
/*						else {	// alt.c
							// Copy from-pool value to [from][to] place in 2d-array:
							transfer_st_2d[index(from, to)].copy(transfer_st_from[from]);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Year %d after transfer from donor pool: C cont in transfer from st %d to %d=%f\n", date.year, from, to, transfer_st_2d[index(from, to)].ccont());//.15f
							dprintf("\treceiving_fraction=%f\n", receiving_fraction_2d);
#endif
						}
*/					}
				}	
			}
		}

		// create and kill stands at landcover change
		stand_dynamics(gridcell);
		error += check_fractions4(gridcell);
		if(error)
			fail("Fraction error after expand_stands()\n");

		// transfer litter etc. of reduced stands to expanding stands at landcover change
		for(int to=0; to<nst; to++) {

			StandType& st = stlist[to];

			if(gridcell.pool_from_all_landcovers[st.landcover]) {
				// Alt.a: All donor st:s are pooled into receiving st:s
				if(st.gross_frac_increase > 0.0) {
#ifdef PRINT_GROSS_LC_CHANGE_INFO
					dprintf("Before receiving_stand_change(): C cont in transfer to st %d=%.15f\n", to, transfer_st[to].ccont());
#endif
					receiving_stand_change(gridcell, transfer_st[to], LCchangeCtransfer, -1, to);
				}
			}
			else {

				for(int from=0; from<nst; from++) {

					if(st_frac_transfer[index(from, to)] > 0.0) {
						StandType& st_from = stlist[from];
//						if(!gridcell.pool_to_all_landcovers[st_from.landcover]) {
						if(!(gridcell.pool_to_all_landcovers[st_from.landcover] || st_from.nstands == 1)) {
						// Alt.b: Transfers between donor and receptor st:s are independent:
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Before receiving_stand_change(): C cont in transfer from st %d to %d=%.15f\n", from, to, transfer_st_2d[index(from, to)].ccont());
#endif
							receiving_stand_change(gridcell, transfer_st_2d[index(from, to)], LCchangeCtransfer, st.landcover, to, st_frac_transfer[index(from, to)] / st.gross_frac_increase);
						}
						else {
						// Alt.c: Donor stands in a stand type are pooled before transfer to recipients:
#ifdef PRINT_GROSS_LC_CHANGE_INFO
							dprintf("Before receiving_stand_change(): C cont in transfer from st %d to %d=%.15f\n", from, to, transfer_st_from[from].ccont());
#endif
							receiving_stand_change(gridcell, transfer_st_from[from], LCchangeCtransfer, st.landcover, to, st_frac_transfer[index(from, to)] / st.gross_frac_increase);
						}
					}
				}
			}
		}

		if(transfer_st_2d)
			delete[] transfer_st_2d;
		if(transfer_st)
			delete[] transfer_st;
		if(transfer_st_from)
			delete[] transfer_st_from;
	}

#ifdef PRINT_GROSS_LC_CHANGE_INFO
	dprintf("\n");
#endif

	double ccont_tot = gridcell.ccont();
	double cflux_tot = gridcell.cflux();

	for(unsigned int i=0; i<gridcell.size(); i++) {
		Stand& stand = gridcell[i];

		double ccont_stand = stand.ccont(stand.scale_LC_change);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
		if(stand.landcover == NATURAL) 
			dprintf("Year %d: natural stand %d ccont=%.15f; age=%d, frac=%f", date.year, stand.id, ccont_stand, date.year - stand.first_year, stand.get_gridcell_fraction());
		else
			dprintf("Year %d: stand %d st %d, ccont=%.15f; age=%d, frac=%f", date.year, stand.id, stand.stid, ccont_stand, date.year - stand.first_year, stand.get_gridcell_fraction());
		if(stand.gross_frac_decrease)
			dprintf(", red %f", stand.gross_frac_decrease - stand.cloned_fraction);
		else
			dprintf("              ");
		if(stand.gross_frac_increase)
			dprintf(", inc %f", stand.gross_frac_increase + stand.cloned_fraction);
		dprintf("\n");
#endif
	}
#ifdef PRINT_GROSS_LC_CHANGE_INFO
	dprintf("\nYear %d: ccont_tot after luc=%.15f\n", date.year, ccont_tot);
	dprintf("Year %d: cflux_tot after luc=%.15f\n", date.year, cflux_tot);
	dprintf("Year %d: total C after luc=%.15f\n", date.year, ccont_tot + cflux_tot);
	dprintf("Year %d: c balance after luc=%.15f\n\n", date.year, cflux_tot + ccont_tot - ccont_tot_1);
#endif
	if(fabs(cflux_tot - cflux_tot_1 + ccont_tot - ccont_tot_1) > 1.0e-12)
		dprintf("WARNING ! C balance after lcc off\n");

	double ncont_tot = gridcell.ncont();
	double nflux_tot = gridcell.nflux();;

	for(unsigned int i=0; i<gridcell.size(); i++) {
		Stand& stand = gridcell[i];

		double ncont_stand = stand.ncont(stand.scale_LC_change);
#ifdef PRINT_GROSS_LC_CHANGE_INFO
		if(stand.landcover == NATURAL) 
			dprintf("Year %d: natural stand %d ncont=%.15f\n", date.year, stand.id, ncont_stand);
		else
			dprintf("Year %d: stand %d st %d, ncont=%.15f\n", date.year, stand.id, stand.stid, ncont_stand);
#endif
	}

#ifdef PRINT_GROSS_LC_CHANGE_INFO
	dprintf("\nYear %d: ncont_tot after luc=%.15f\n", date.year, ncont_tot);
	dprintf("Year %d: nflux from luc=%.15f\n", date.year, gridcell.anflux_landuse_change + gridcell.anflux_harvest_slow);
	dprintf("Year %d: total N after luc=%.15f\n", date.year, ncont_tot + nflux_tot);
	dprintf("Year %d: N balance after luc=%.15f\n\n", date.year, nflux_tot - nflux_tot_1 + ncont_tot - ncont_tot_1);
#endif
	if(fabs(nflux_tot - nflux_tot_1 + ncont_tot - ncont_tot_1) > 1.0e-12)
		dprintf("WARNING ! N balance after lcc off\n");

	if(st_frac_transfer)
		delete[] st_frac_transfer;
	if(primary_st_frac_transfer)
		delete[] primary_st_frac_transfer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  End of Landcover stand dynamics and C&N-partitioning  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void crop_nfert(Patch& patch) {

	Gridcell& gridcell = patch.stand.get_gridcell();

	pftlist.firstobj();
	// Loop through PFTs
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if(patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

			double nfert = pft.N_appfert;
			if(gridcellpft.Nfert_read >= 0.0) {
				nfert = gridcellpft.Nfert_read;
			}
			if(!ppftcrop.fertilised[0] && ppftcrop.dev_stage > 0.0 && ppftcrop.growingseason){
				patch.dnfert = nfert * (1.0 - pft.fertrate[0] - pft.fertrate[1]);
				ppftcrop.fertilised[0] = true;
			}
			else if(!ppftcrop.fertilised[1] && ppftcrop.dev_stage > pft.fert_stages[0] && ppftcrop.growingseason){
				patch.dnfert = nfert * pft.fertrate[0];
				ppftcrop.fertilised[1] = true;
			}
			else if(!ppftcrop.fertilised[2] && ppftcrop.dev_stage > pft.fert_stages[1] && ppftcrop.growingseason ){
				patch.dnfert = nfert * (pft.fertrate[1]);
				ppftcrop.fertilised[2] = true;
			}
			else {
				patch.dnfert = 0.0;
			}
//			if(date.day == ppftcrop.bicdate)
//				patch.dnfert = 0.003;
			patch.anfert += patch.dnfert;
		}
		pftlist.nextobj();
	}
}

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

	if(pft.ifsdautumn  && pft.forceautumnsowing != SPRINGSOWING) {							// TeWW,TeRa:
	
		// Use autumn sowing if first_autumndate20 is set (autumn conditions met during the past 20 years):
		if(!((gridcellpft.first_autumndate20 == climate.testday_temp || gridcellpft.first_autumndate20 == climate.coldestday) && 
				gridcellpft.first_autumndate % 365 == gridcellpft.first_autumndate20)) {

			gridcellpft.sdatecalc_temp = gridcellpft.first_autumndate20;
			gridcellpft.wintertype = true;
		}
		// if not, use spring sowing
		else {	// if(gridcellpft.first_autumndate20==climate.coldestday)
		
			if(!((gridcellpft.last_springdate20 == climate.testday_temp || gridcellpft.last_springdate20 == climate.coldestday) && 
					gridcellpft.last_springdate == gridcellpft.last_springdate20)
					 && pft.forceautumnsowing != AUTUMNSOWING) {

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
		if(dayinperiod(gridcellpft.sdatecalc_temp, climate.testday_temp, gridcellpft.hlimitdate_default)
			 && pft.forceautumnsowing != AUTUMNSOWING) {

			gridcellpft.sdatecalc_temp = gridcellpft.last_springdate20;	// use last_springdate20 disregarding earlier choices	
			gridcellpft.wintertype = false;
		}

		// Forced sowing date read from input file.
		// Calculated value used if value for pft not found in file.
		if(readsowingdates && pft.readsowingdate && gridcellpft.sdate_force >= 0) {

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
//	if(!strncmp(pft.name,"TeWW", strlen("TeWW")) && climate.mtemp_min20 > -20.0)	//Test continuous grass
		gridcellpft.sdatecalc_temp = -1;
#endif

	gridcellpft.springoccurred = false;	
	gridcellpft.vernstartoccurred = false;	
	gridcellpft.vernendoccurred = false;	
	gridcellpft.autumnoccurred = false;
}

/// Sets sdatecalc_prec first day of the rain period
/** Called from crop_sowing_gridcell() each day if NEWSOWINGDATE is not defined
 */
void set_sdatecalc_prec(Climate& climate, Gridcellpft& gridcellpft)
{
	Pft& pft = gridcellpft.pft;

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

			int swindow_temp[2];
			int swindow_prec[2];

			// Calculate temperature-dependent sowing windows
			if(gridcellpft.sdatecalc_temp != -1) {

				// Set sowing window around sdatecalc_temp
				swindow_temp[0] = stepfromdate(gridcellpft.sdatecalc_temp, -15);
				swindow_temp[1] = stepfromdate(gridcellpft.sdatecalc_temp, 15);

				if(!gridcellpft.wintertype && dayinperiod(swindow_temp[0], stepfromdate(climate.coldestday, -100), climate.coldestday)) {

					swindow_temp[0] = climate.coldestday;
//					swindow_temp[0] = gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(swindow_temp[1], stepfromdate(climate.coldestday, -100), climate.coldestday))
						swindow_temp[1] = climate.coldestday;
				}

				if(gridcellpft.wintertype && dayinperiod(swindow_temp[1], climate.coldestday, stepfromdate(climate.coldestday, 100))) {

					swindow_temp[1] = climate.coldestday;
//					swindow_temp[1] = gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(swindow_temp[0], climate.coldestday, stepfromdate(climate.coldestday, 100)))
						swindow_temp[0] = climate.coldestday;
				}
			}

			// Calculate precipitation-dependent sowing windows
			monthdates(swindow_prec[0],swindow_prec[1],sow_month);
			// A conservative choice to expand the sowing window
			swindow_prec[0] = stepfromdate(swindow_prec[0], -15);


			// Determine, based upon site climate seasonality, if sowing in rainfed stands should be triggered by 
			// temperature or precipitation, or whether to use a default sowing date.

			bool temp_sdate = false, prec_sdate = false, def_sdate = false;

			if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC)
				temp_sdate = true;
			else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range != WET)
				prec_sdate = true;
			else // if(seasonality == SEASONALITY_NO) || (seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP) && climate.prec_range == WET)
				def_sdate = true;

			if(temp_sdate) {
				gridcellpft.swindow[0] = swindow_temp[0];
				gridcellpft.swindow[1] = swindow_temp[1];
			}
			else if(prec_sdate) {
				gridcellpft.swindow[0] = swindow_prec[0];
				gridcellpft.swindow[1] = swindow_prec[1];
			}
			else if(def_sdate) {
				gridcellpft.swindow[0] = gridcellpft.sdate_default;
				gridcellpft.swindow[1] = stepfromdate(gridcellpft.sdate_default, 15);
			}

			// Rules for irrigated crops:

			// Different sowing date options for irrigated crops at sites with climate.seasonality == SEASONALITY_PRECTEMP:
			// 1. use temperature-dependent sowing limits (define IRRIGATED_USE_TEMP_SDATE)
			// 2. use precipitation-triggered sowing (IRRIGATED_USE_TEMP_SDATE undefined)

			bool irr_use_temp_sdate = false;

#if defined IRRIGATED_USE_TEMP_SDATE
			if(seasonality == SEASONALITY_PRECTEMP)
				irr_use_temp_sdate = true;
#endif

			if(irr_use_temp_sdate) {
				gridcellpft.swindow_irr[0] = swindow_temp[0];
				gridcellpft.swindow_irr[1] = swindow_temp[1];
			}
			else {
				gridcellpft.swindow_irr[0] = gridcellpft.swindow[0];
				gridcellpft.swindow_irr[1] = gridcellpft.swindow[1];
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

	climate.aprec = 0.0;

	for(m=0; m<12; m++) {

		// 1) this year
		climate.mtemp20[m] = climate.mtemp_year[m];
		climate.mprec20[m] = climate.mprec_year[m];
		//historic
		climate.mtemp20[m] = climate.hmtemp_20[m].lastadd();
		climate.mprec20[m] = climate.hmprec_20[m].lastadd();
		climate.aprec += climate.hmprec_20[m].lastadd();
		climate.mpet_year[m] = climate.hmeet_20[m].lastadd()*PRIESTLEY_TAYLOR;
		//
		climate.mpet20[m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
			climate.mprec_pet20[m] = climate.mprec_year[m] / climate.mpet_year[m];
		else
			climate.mprec_pet20[m] = 0.0;


/*		if(climate.mprec_year[m] / climate.mpet_year[m] < mprec_petmin_thisyear)
			mprec_petmin_thisyear = climate.mprec_year[m] / climate.mpet_year[m];
		if(climate.mprec_year[m] / climate.mpet_year[m] > mprec_petmax_thisyear)
			mprec_petmax_thisyear = climate.mprec_year[m] / climate.mpet_year[m];
*/
		if(climate.hmprec_20[m].lastadd() / climate.mpet_year[m] < mprec_petmin_thisyear)
			mprec_petmin_thisyear = climate.hmprec_20[m].lastadd() / climate.mpet_year[m];
		if(climate.hmprec_20[m].lastadd() / climate.mpet_year[m] > mprec_petmax_thisyear)
			mprec_petmax_thisyear = climate.hmprec_20[m].lastadd() / climate.mpet_year[m];

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

/*		climate.mtemp_20[19][m] = climate.mtemp_year[m];
		climate.mprec_20[19][m] = climate.mprec_year[m];
*/
		climate.mtemp_20[19][m] = climate.hmtemp_20[m].lastadd();
		climate.mprec_20[19][m] = climate.hmprec_20[m].lastadd();

		climate.mpet_20[19][m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
//			climate.mprec_pet_20[19][m] = climate.mprec_year[m] / climate.mpet_year[m];
			climate.mprec_pet_20[19][m] = climate.hmprec_20[m].lastadd() / climate.mpet_year[m];
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

	Gridcell& gridcell = patch.stand.get_gridcell();
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];

	patchpft.set_cropphen()->sdate = gridcellpft.sdatecalc_temp;

	if(dayinperiod(gridcellpft.sdatecalc_temp, climate.testday_temp, gridcellpft.hlimitdate_default) && pft.forceautumnsowing == AUTUMNSOWING)
		patchpft.cropphen->hlimitdate = stepfromdate(patchpft.cropphen->sdate, - 1);
}

/// old precipitation-dependent sowing date method (Bondeau et al. 2007)
void Crop_sowing_date_prec(Patch& patch, Pft& pft) {

	Gridcell& gridcell = patch.stand.get_gridcell();
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];
	int first_sowdate;
	int last_sowdate;

	if(patch.stand.pft[pft.id].irrigated)  {
		gridcellpft.sdatecalc_prec = gridcellpft.sdate_default;
		gridcellpft.precoccurred = true;
	}

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
	}
}

/// Sowing date method from Waha et al. 2010
/** Enters here every day when growingseason==false
 */
void Crop_sowing_date_new(Patch& patch, Pft& pft) {

	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	Gridcell& gridcell = patch.stand.get_gridcell();
	Gridcellpft& gridcellpft = gridcell.pft[pft.id];
	Standpft& standpft = patch.stand.pft[pft.id];
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
	if(seasonality == SEASONALITY_TEMP || seasonality == SEASONALITY_TEMPPREC || seasonality == SEASONALITY_PRECTEMP && standpft.irrigated)
		temp_sdate = true;
	else if((seasonality == SEASONALITY_PREC || seasonality == SEASONALITY_PRECTEMP && !standpft.irrigated) && climate.prec_range != WET)
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
				if(climate.prec > 0.1 || standpft.irrigated)
					ppftcrop.sdate = date.day;
			}
			else // if(def_sdate)
				ppftcrop.sdate = date.day; // first day of sowing window
		}
		else	 // last day of sowing window	
			ppftcrop.sdate = date.day;
	}

	// 
	if(standpft.sdate_force >= 0)
		patchpft.cropphen->sdate = standpft.sdate_force;

	// calculation of hucountend (last day of heat unit sampling period); sdate is first day
	// NB. currently the sampling periods (which roughly correspond to growing periods) are unrealistically long,
	// Since shorter growing periods result in lower yield, a revision of this section will also make a revision of crop productivity necessary
	if(date.day == ppftcrop.sdate) {

		if(gridcellpft.sdate_default > gridcellpft.hlimitdate_default)
			length_growseas_def = gridcellpft.hlimitdate_default + 365 - gridcellpft.sdate_default;
		else
			length_growseas_def = gridcellpft.hlimitdate_default - gridcellpft.sdate_default;

//		if(stlist[patch.stand.stid].rotation.multicrop && gridcellpft.multicrop)
		if(stlist[patch.stand.stid].rotation.multicrop)
			length_growseas_def = 150;

		if(prec_sdate)
			ppftcrop.hlimitdate = stepfromdate(date.day, length_growseas_def);
		else
			ppftcrop.hlimitdate = gridcellpft.hlimitdate_default;

		length_growseas_def = min(length_growseas_def, 245);				// set an upper limit of 245 for the growing season

		if(pft.ifsdautumn) { // winter crops
//		if(pft.ifsdautumn && gridcellpft.wintertype) { // winter crops: try this

			if(temp_sdate)
				ppftcrop.hucountend = stepfromdate(ppftcrop.hlimitdate, -20);
			else if(prec_sdate && climate.prec_seasonality <= DRY_WET) { // dry some time during the year
				if(standpft.irrigated)
					ppftcrop.hucountend = stepfromdate(date.day, min(length_growseas_def, 230));
				else
					ppftcrop.hucountend = stepfromdate(date.day, min(length_growseas_def, 210));		 // shorter growing period when risk for water stress.
			}
			else if(def_sdate)
				ppftcrop.hucountend = stepfromdate(date.day, 230);
		}
		else if(!strncmp(pft.name,"TrRi", strlen("TrRi"))) // rice
			ppftcrop.hucountend = stepfromdate(date.day, min(length_growseas_def, 230));
		else { // all other crops
			if(prec_sdate && climate.prec_seasonality <= DRY_WET) {	// dry some time during the year

			if(standpft.irrigated)
					ppftcrop.hucountend = stepfromdate(date.day, length_growseas_def);
				else
					ppftcrop.hucountend = stepfromdate(date.day, min(length_growseas_def, 210)); // shorter growing period when risk for water stress.
			}
			else
				ppftcrop.hucountend = stepfromdate(date.day, length_growseas_def);
		}	

		if(standpft.sdate_force >= 0) {
			if(standpft.hdate_force >= 0) {

				patchpft.cropphen->hlimitdate = standpft.hdate_force;
				patchpft.cropphen->hucountend = standpft.hdate_force;
			}
		}
	}
}

/// handles sowing date calculations for crop pft:s on patch level
void crop_sowing_patch(Patch& patch) {

	pftlist.firstobj();
	Gridcell& gridcell = patch.stand.get_gridcell();
	Climate& climate = gridcell.climate;

	// Loop through PFTs
	while(pftlist.isobj) {

		Pft& pft = pftlist.getobj();
		Patchpft& patchpft = patch.pft[pft.id];
		Gridcellpft& gridcellpft = gridcell.pft[pft.id];

		if(patch.stand.pft[pft.id].active && pft.phenology == CROPGREEN) {

			cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

			if(date.day == climate.testday_temp) {

				if(patch.stand.pft[pft.id].irrigated) {
					patchpft.swindow[0] = gridcellpft.swindow_irr[0];
					patchpft.swindow[1] = gridcellpft.swindow_irr[1];
				}
				else {
					patchpft.swindow[0] = gridcellpft.swindow[0];
					patchpft.swindow[1] = gridcellpft.swindow[1];
				}
			}

			if(patch.stand.pftid == pft.id && !ppftcrop.growingseason) {

				// copy sowing window from gridcellpft
				if(date.day == stepfromdate(ppftcrop.hdate, 1) && ppftcrop.hdate != -1 || date.day == climate.testday_temp) {

					if(gridcellpft.swindow[0] == -1 && date.year) {
						gridcellpft.sowing_restriction = true;
						ppftcrop.hdate = -1;
						ppftcrop.eicdate = -1;	//redundant
						patch.stand.isrotationday = true;
					}
					else {
						if(!patch.stand.infallow)
							gridcellpft.sowing_restriction = false;
						else if(patch.stand.ndays_inrotation > 180)
							patch.stand.isrotationday = true;
					}
				}

				if(!gridcellpft.sowing_restriction)	{

#if defined NEWSOWINGDATE
 					// new sowing date method (Waha et al. 2010)
					Crop_sowing_date_new(patch, pft);
#else				// old sowing date method (Bondeau et al. 2007)
					Crop_sowing_date(patch, pft);
#endif
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
	const Climate& climate = patch.get_climate();
	double phu_last_year = ppftcrop.phu;

	ppftcrop.husum = 0.0;	
	ppftcrop.vrf = 1.0;
	ppftcrop.vdsum = 0;
	ppftcrop.prf = 1.0;

	ppftcrop.pvd = pft.pvd;		// default; kept for TrMi, TePu, TeSb, TrMa, TeSo, TrPe
	ppftcrop.phu = pft.phu;
	ppftcrop.tb = pft.tb;

	ppftcrop.vdsum_alloc=0.0;
	ppftcrop.vd=0.0;
	ppftcrop.dev_stage=0.0;

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
		else {	// spring sowing
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
			if (patch.stand.get_gridcell().get_lon() < 60.0 || patch.stand.get_gridcell().get_lat() > 30.0)
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
	const Climate& climate = patch.get_climate();

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

/// Calculation of development stage
/** Accumulation of development during sampling period. TODO Add reference
 */ 
void calc_ds(Patch& patch, Pft& pft) {

	Patchpft& patchpft = patch.pft[pft.id];

	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	const Climate& climate = patch.get_climate();
	double T = climate.temp;

	// account for vernalization if needs for vernalization not yet satisfied	//trg=tb for crops other than TeWW and TeRa and don't enter here
	if (ppftcrop.vdsum_alloc < 1)	{				

		if (T > pft.T_vn_min && T < pft.T_vn_max) {
			double alpha_v = log(2.0) / (log((pft.T_vn_max - pft.T_vn_min) / (pft.T_vn_opt - pft.T_vn_min)));
			double fT_v = (2.0 * pow((T - pft.T_vn_min),alpha_v) * pow((pft.T_vn_opt - pft.T_vn_min), alpha_v) - pow((T - pft.T_vn_min), 2.0 * alpha_v)) / pow((pft.T_vn_opt - pft.T_vn_min),2.0 * alpha_v);
			ppftcrop.vd = ppftcrop.vd + fT_v;
			ppftcrop.vdsum_alloc = min(1.0, pow((double)ppftcrop.vd, 5.0) / (pow(22.5, 5.0) + pow((double)ppftcrop.vd, 5.0)));
		}																	
	}

	double e = 2.71828183;
	double P = climate.daylength_save[date.day];
	double fP = 0;

	if(pft.photo[2] > 0) //short day plant 
	{
		if(P < pft.photo[0])
			fP = 1;
		else
			fP = min(1.0, pow(e,(-pft.photo[1] * (P - pft.photo[0]))));
	} 
	else //long day plant
	{
		if(P < pft.photo[0])
			fP = 0;
		else
			fP = min(1.0, 1.0 - pow(e,(-pft.photo[1] * (P - pft.photo[0]))));
	}
	double fT = 0.0;
	double T_min = pft.T_veg_min;
	double T_opt = pft.T_veg_opt;
	double T_max = pft.T_veg_max;

	if(ppftcrop.dev_stage >= 1) {
		T_min = pft.T_rep_min;
		T_opt = pft.T_rep_opt;
		T_max = pft.T_rep_max;
	}

	double alpha = log(2.0) / (log((T_max - T_min) / (T_opt - T_min)));

	if(T > T_min && T < T_max)
		fT = min(1.0, (2.0 * pow((T - T_min), alpha) * pow((T_opt - T_min), alpha) - pow((T - T_min), 2.0 * alpha)) / pow((T_opt - T_min), 2.0 * alpha));

	double DR = 0.0;
	if (ppftcrop.dev_stage < 1.0)
		DR = pft.dev_rate_veg * ppftcrop.vdsum_alloc * fP * fT;
	else
		DR = pft.dev_rate_rep * fT;

	ppftcrop.dev_stage = min(2.0, ppftcrop.dev_stage + DR);			
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
		Standpft& standpft = patch.stand.pft[pft.id];
		Gridcell& gridcell = patch.stand.get_gridcell();
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
				bool force_harvest = date.day == standpft.hdate_force;

				// before maturity is reached
				bool pre_maturity = (ifnlim_lc[CROPLAND]) ? ppftcrop.dev_stage < 2.0 : ppftcrop.husum < ppftcrop.phu;

				if(pre_maturity && dayinperiod(date.day, ppftcrop.sdate, stepfromdate(ppftcrop.hlimitdate, -1)) && !force_harvest) {

					// count accumulated heat units after sowing date
					calc_hu(patch, pft);

					if(ifnlim_lc[CROPLAND])
						calc_ds(patch, pft);

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

					ppftcrop.fertilised[0] = false;
					ppftcrop.fertilised[1] = false;
					ppftcrop.fertilised[2] = false;

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

			if(patch.stand.pftid == pft.id && stlist[patch.stand.stid].intercrop == NATURALGRASS) {

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
	Gridcell& gridcell = patch.stand.get_gridcell();
	Climate& climate = gridcell.climate;
	Patchpft& patchpft = patch.pft[pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

	if (pft.phenology == CROPGREEN) {
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
	else if(pft.phenology == ANY) { // crop grasses using standard guess phenology calculation

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

			if(indiv.growingseason())
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

				if(!ppftcrop.senescence || ifnlim_lc[CROPLAND])
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
		indiv.nstore_longterm, indiv.max_n_storage,
		true);

	indiv.cmass_leaf_post_turnover = cropindiv.grs_cmass_leaf;
	indiv.cmass_root_post_turnover = cropindiv.grs_cmass_root;

	// Nitrogen longtime storage
	// Nitrogen approx retranslocated next season
	double retransn_nextyear = cmass_leaf_pre_turnover * indiv.pft.turnover_leaf / cton_leaf_bg * nrelocfrac +
		cmass_root_pre_turnover * indiv.pft.turnover_root / cton_root_bg * nrelocfrac;

	// Max longterm nitrogen storage
	indiv.max_n_storage = max(0.0, min(cmass_root_pre_turnover * indiv.pft.fnstorage / cton_leaf_bg, retransn_nextyear));

	// Scale this year productivity to max storage
	if (grs_npp > 0.0) {
		indiv.scale_n_storage = max(indiv.max_n_storage * 0.1, indiv.max_n_storage - retransn_nextyear) * cton_leaf_bg / grs_npp;
	}

	indiv.nstore_labile = indiv.nstore_longterm;
	indiv.nstore_longterm = 0.0;
}

void crop_allocation_WE(cropphen_struct& ppftcrop, Individual& indiv) {

	ppftcrop.dev_stage = max(0.0,min(2.0, -0.595 * pow(ppftcrop.fphu, 2.0) + 2.595 * ppftcrop.fphu));

	double t = 0.0;
	if(ppftcrop.fphu < 0.4367) {
		t = -0.07 + 2.45 * ppftcrop.fphu;
	} 
	else {
		t = 0.06 + 2.0 * ppftcrop.fphu;
		t = 0.2247 + 1.7753 * ppftcrop.fphu;
	}
	ppftcrop.dev_stage = max(0.0,min(2.0,t));

	double f1 = min(1.0, max(0.0, richards_curve(indiv.pft.a1, indiv.pft.b1, indiv.pft.c1, indiv.pft.d1, ppftcrop.dev_stage)));
	double f2 = min(1.0, max(0.0, richards_curve(indiv.pft.a2, indiv.pft.b2, indiv.pft.c2, indiv.pft.d2, ppftcrop.dev_stage)));
	double f3 = min(1.0, max(0.0, richards_curve(indiv.pft.a3, indiv.pft.b3, indiv.pft.c3, indiv.pft.d3, ppftcrop.dev_stage)));

	if(indiv.daily_cmass_leafloss > 0.0)
		f2 *= f2 * f2;

	ppftcrop.f_alloc_root = f1 * (1-f3);
	ppftcrop.f_alloc_leaf = f2 * (1-f1)*(1-f3);
	ppftcrop.f_alloc_stem = (1.0 - f2)*(1.0 - f1)*(1.0 - f3);
	ppftcrop.f_alloc_horg = f3;
}

/// Daily allocation routine for crops with nitrogen limitation
/** Allocates daily npp to leaf, roots and harvestable organs
 *  Equations are from Neitsch et al. 2002.
 */
void allocation_crop_nlim(Individual& indiv, double cmass_seed, double nmass_seed) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());
	Patch& patch = indiv.vegetation.patch;
	Patchpft& patchpft = patch.pft[indiv.pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());
	double cmass_extra = 0.0;

	if(ppftcrop.growingseason) {

		// report seed fluxes
		indiv.report_flux(Fluxes::SEEDC, -cmass_seed);
		indiv.report_flux(Fluxes::SEEDN, -nmass_seed);

		// add seed carbon
		cmass_extra += cmass_seed;

		// add seed nitrogen
		indiv.nmass_leaf += nmass_seed / 2.0;
		indiv.nmass_root += nmass_seed / 2.0;

		crop_allocation_WE(ppftcrop, indiv);

		if(indiv.dnpp < 0.0){
#define STEMSUGAR
			//#undef STEMSUGAR
#ifdef STEMSUGAR
			if(-indiv.dnpp < cropindiv.grs_cmass_agpool) {
				cropindiv.grs_cmass_agpool -= -indiv.dnpp;
				cropindiv.ycmass_agpool -= -indiv.dnpp;
				indiv.dnpp = 0.0;
			} 
			else {
				indiv.report_flux(Fluxes::NPP, (-indiv.dnpp - cropindiv.grs_cmass_agpool));
				cropindiv.ycmass_agpool -= cropindiv.grs_cmass_agpool;
				//TODO Kill the individual if ag pool is zero.
				cropindiv.grs_cmass_agpool = 0.0;
				indiv.dnpp = 0.0;
			}
#else
			indiv.report_flux(Fluxes::NPP, -indiv.dnpp);
			indiv.dnpp = 0.0;
#endif
		}

#ifdef STEMSUGAR
#define USESTEMSUGAR
#ifdef USESTEMSUGAR

		if (cropindiv.grs_cmass_agpool > 0.0 && patchpft.cropphen->f_alloc_horg > 0.95) {
			cmass_extra += 0.1 * cropindiv.grs_cmass_agpool;
			cropindiv.ycmass_agpool -= 0.1 * cropindiv.grs_cmass_agpool;
			cropindiv.grs_cmass_agpool *= 0.9;
		}

#endif
#endif
		indiv.daily_cmass_rootloss = 0.0;
		indiv.daily_nmass_rootloss = 0.0;

		if (indiv.daily_cmass_leafloss > 0.0) {

			cropindiv.dcmass_leaf = (indiv.dnpp + cmass_extra) * patchpft.cropphen->f_alloc_leaf - indiv.daily_cmass_leafloss;
			cropindiv.grs_cmass_dead_leaf += indiv.daily_cmass_leafloss;
			cropindiv.ycmass_dead_leaf += indiv.daily_cmass_leafloss;
			if (indiv.daily_cmass_leafloss / 100.0<indiv.nmass_leaf) {
				cropindiv.nmass_dead_leaf += indiv.daily_cmass_leafloss / 100.0; //TODO super low C:N in the dead leaf
				cropindiv.ynmass_dead_leaf += indiv.daily_cmass_leafloss / 100.0;
				indiv.nmass_leaf -= indiv.daily_cmass_leafloss / 100.0;
			}
			//cropindiv.grs_cmass_leaf -= indiv.daily_cmass_leafloss;
			double new_CN = (cropindiv.grs_cmass_leaf + cropindiv.dcmass_leaf) / indiv.nmass_leaf;
			// If the result is smaller (higher [N]) than the min C:N then that N is
			// put in to the ag N pool
			if( new_CN < indiv.pft.cton_leaf_min ) {

#define THREEQUARTER
#ifdef THREEQUARTER
				indiv.daily_nmass_leafloss = max(0.0, indiv.nmass_leaf - (cropindiv.grs_cmass_leaf + cropindiv.dcmass_leaf) / (1.33 * indiv.pft.cton_leaf_min));
#else
				indiv.daily_nmass_leafloss = max(0.0, indiv.nmass_leaf - (cropindiv.grs_cmass_leaf + cropindiv.dcmass_leaf) / indiv.pft.cton_leaf_min);
#endif
				if(indiv.daily_nmass_leafloss > indiv.nmass_leaf) {
					indiv.daily_nmass_leafloss = 0.0;
				}
			} else {
				indiv.daily_nmass_leafloss = 0.0;
			}
			// Very experimental root senescence
			// N and C loss when root senescence is allowed f_HO > 0.5
#define ROOTLOSS
			//#undef ROOTLOSS
#ifdef ROOTLOSS
			//d3, the DS after which more than half of the daily assimilates are going to the grains.
			if(patchpft.cropphen->dev_stage > indiv.pft.d3) {
				//only have root senescence when leaf scenescence har occured
				if (indiv.daily_nmass_leafloss > 0.0) {
					double kC = 0.0;
					double kN = 0.0;
					//The root senescence is proportional to that of the leaves
					if(indiv.nmass_leaf > 0.0) {
						kN = indiv.daily_nmass_leafloss / indiv.nmass_leaf;
					}
					if (indiv.cmass_leaf_today() > 0.0) {
						kC = indiv.daily_cmass_leafloss / indiv.cmass_leaf_today();
					}
					indiv.daily_cmass_rootloss = indiv.cmass_root_today() * kC;
					indiv.daily_nmass_rootloss = indiv.nmass_root * kN;
				}
			}
#endif
			indiv.nmass_leaf -= indiv.daily_nmass_leafloss;
			cropindiv.nmass_agpool += indiv.daily_nmass_leafloss;
		} 
		else {
			cropindiv.dcmass_leaf = (indiv.dnpp + cmass_extra) * patchpft.cropphen->f_alloc_leaf;
		}
		cropindiv.dcmass_stem = (indiv.dnpp + cmass_extra) * patchpft.cropphen->f_alloc_stem;

		if (indiv.daily_cmass_rootloss > indiv.cmass_root_today())
			indiv.daily_cmass_rootloss = 0.0;

		cropindiv.dcmass_root = (indiv.dnpp + cmass_extra) * patchpft.cropphen->f_alloc_root - indiv.daily_cmass_rootloss;

		//TODO
		patch.soil.sompool[SOILMETA].cmass += indiv.daily_cmass_rootloss;

		if (indiv.daily_nmass_rootloss < indiv.nmass_root) {
			indiv.nmass_root -= indiv.daily_nmass_rootloss;
			cropindiv.nmass_agpool += indiv.daily_nmass_rootloss * 0.5; // 50% of the N in the lost root is retranslocated.
			patch.soil.sompool[SOILMETA].nmass += indiv.daily_nmass_rootloss * 0.5;//The rest is going in to litter
		}
		if (indiv.daily_cmass_rootloss > 0.0){
			patch.is_litter_day = true;
		}

		cropindiv.dcmass_ho = (indiv.dnpp + cmass_extra) * patchpft.cropphen->f_alloc_horg;
		cropindiv.dcmass_plant = cropindiv.dcmass_ho + cropindiv.dcmass_root + cropindiv.dcmass_stem + cropindiv.dcmass_leaf;

		cropindiv.ycmass_leaf += cropindiv.dcmass_leaf;
		cropindiv.ycmass_root += cropindiv.dcmass_root;
		cropindiv.ycmass_ho += cropindiv.dcmass_ho;
		cropindiv.ycmass_plant += cropindiv.dcmass_plant;

		cropindiv.grs_cmass_leaf += cropindiv.dcmass_leaf;
#ifdef STEMSUGAR
		cropindiv.grs_cmass_stem += (1.0 - 0.4) * cropindiv.dcmass_stem;
		cropindiv.ycmass_stem += (1.0 - 0.4) * cropindiv.dcmass_stem;
		cropindiv.grs_cmass_agpool += 0.4 * cropindiv.dcmass_stem;
		cropindiv.ycmass_agpool += 0.4 * cropindiv.dcmass_stem;
#else
		cropindiv.grs_cmass_stem += cropindiv.dcmass_stem;
		cropindiv.ycmass_stem += cropindiv.dcmass_stem;
#endif
		cropindiv.grs_cmass_root += cropindiv.dcmass_root;
		cropindiv.grs_cmass_ho += cropindiv.dcmass_ho;
		cropindiv.grs_cmass_plant += cropindiv.dcmass_plant;

		double ndemand_ho = 0.0;
		double avail_leaf_N = max(0.0, (1.0 / indiv.cton_leaf(false) - 1.0 / indiv.pft.cton_leaf_max) * indiv.cmass_leaf_today());
		double avail_root_N = max(0.0, (1.0 / indiv.cton_root(false) - 1.0 / indiv.pft.cton_root_max) * indiv.cmass_root_today());
		double avail_stem_N = max(0.0,cropindiv.nmass_agpool - 1.0 / indiv.pft.cton_stem_max * cropindiv.grs_cmass_stem);
		double avail_N = avail_leaf_N + avail_root_N + avail_stem_N;
		if (avail_N > 0.0 && cropindiv.dcmass_ho > 0.0) {
			ndemand_ho = cropindiv.dcmass_ho / indiv.pft.cton_leaf_avr;
		}
		//N mass to be translocated from leaves and roots
		double trans_leaf_N = 0.0;
		double trans_root_N = 0.0;
		if (ndemand_ho > 0.0) {
			if(avail_stem_N > 0.0) {
				if (ndemand_ho > avail_stem_N) {
					ndemand_ho -= avail_stem_N;
					cropindiv.dnmass_ho += avail_stem_N;
					cropindiv.nmass_agpool -= avail_stem_N;
				} 
				else {
					cropindiv.nmass_agpool -= ndemand_ho;
					cropindiv.dnmass_ho += ndemand_ho;
					ndemand_ho = 0.0;
				}
			}
			//Seligman 1975
			//"willingness" to let go of the N in the organ to meet the demand from the storage organ
			double w = 0.0;
			double w_r = 0.0;
			double w_l = 0.0;
			double w_s = 0.0;
			double y0 = (1.0 / indiv.pft.cton_leaf_min + 1.0 / indiv.pft.cton_leaf_avr) / 2.0;
			double y = 1.0 / indiv.cton_leaf(false);
			double y2 = 1.0 / (1.0 * indiv.pft.cton_leaf_max);
			double z = (y0 - y)/(y0 - y2);
			w_l = 1.0 - max(0.0, min(1.0, pow(1.0 - z, 2.0)));
			y0 = 1,0 / indiv.pft.cton_root_avr;
			y = 1.0/ indiv.cton_root(false);
			y2 = 1.0 / (1.0 * indiv.pft.cton_root_max);
			z = (y0 - y) / (y0 - y2);
			w_r = 1.0 - max(0.0, min(1.0, pow(1.0 - z, 2.0)));
			w_s = w_r + w_l;
			w = min(1.0, w_s);
			if(w_s > 0.0) {
				trans_leaf_N = max(0.0, w_l * w * ndemand_ho / w_s);
				trans_root_N = max(0.0, w_r * w * ndemand_ho / w_s);
				if(trans_leaf_N > avail_leaf_N) {
					trans_leaf_N = avail_leaf_N;
				}
				if(trans_root_N > avail_root_N) {
					trans_root_N = avail_root_N;
				}
				cropindiv.dnmass_ho += trans_leaf_N;
				cropindiv.dnmass_ho += trans_root_N;
			}
		}
		indiv.nmass_leaf -= trans_leaf_N;
		indiv.nmass_root -= trans_root_N;
		cropindiv.nmass_ho += cropindiv.dnmass_ho;
	}
	return;
}

void allocation_crop(Individual& indiv, double cmass_seed, double nmass_seed) {

	cropindiv_struct& cropindiv = *(indiv.get_cropindiv());
	Patch& patch = indiv.vegetation.patch;
	Patchpft& patchpft = patch.pft[indiv.pft.id];
	cropphen_struct& ppftcrop = *(patchpft.get_cropphen());

	nmass_seed = 0.0;	//temporary ?

	// report seed flux
	indiv.report_flux(Fluxes::SEEDC, -cmass_seed);
	indiv.report_flux(Fluxes::SEEDN, -nmass_seed);

	// add seed carbon
	cropindiv.grs_cmass_plant += cmass_seed;
	cropindiv.ycmass_plant += cmass_seed;
	cropindiv.dcmass_plant += cmass_seed;

	// add seed nitrogen
	indiv.nmass_leaf += nmass_seed / 2.0;
	indiv.nmass_root += nmass_seed / 2.0;

	// add today's npp
	cropindiv.dcmass_plant += indiv.dnpp;
	cropindiv.grs_cmass_plant += indiv.dnpp;
	cropindiv.ycmass_plant += indiv.dnpp;

	// allocation to roots
	double froot = indiv.pft.frootstart -(indiv.pft.frootstart - indiv.pft.frootend) * ppftcrop.fphu;	// SWAT 5:2,1,21	
	double grs_cmass_root_old = cropindiv.grs_cmass_root;
	cropindiv.grs_cmass_root = froot * cropindiv.grs_cmass_plant;
	cropindiv.dcmass_root = cropindiv.grs_cmass_root - grs_cmass_root_old;
	cropindiv.ycmass_root += cropindiv.dcmass_root;

	// allocation to harvestable organs
	double grs_cmass_ag = (1.0-froot) * cropindiv.grs_cmass_plant;
	double grs_cmass_ho_old = cropindiv.grs_cmass_ho;

	if(indiv.pft.hiopt <= 1.0)
		cropindiv.grs_cmass_ho = ppftcrop.hi * grs_cmass_ag;									// SWAT 5:2.4.2, 5:2.4.4
	else	// below-ground harvestable organs
		cropindiv.grs_cmass_ho = (1.0 - 1.0 / (1.0 + ppftcrop.hi)) * cropindiv.grs_cmass_plant;	// SWAT 5:2.4.3 8 

	cropindiv.dcmass_ho = cropindiv.grs_cmass_ho - grs_cmass_ho_old;	
	cropindiv.ycmass_ho += cropindiv.dcmass_ho;	

	// allocation to leaves
	double grs_cmass_leaf_old = cropindiv.grs_cmass_leaf;	
	cropindiv.grs_cmass_leaf = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root - cropindiv.grs_cmass_ho;
 
	cropindiv.dcmass_leaf = cropindiv.grs_cmass_leaf - grs_cmass_leaf_old;
	cropindiv.ycmass_leaf += cropindiv.dcmass_leaf;

	// allocation to above-ground pool (currently not used)
	cropindiv.dcmass_agpool = cropindiv.dcmass_plant - cropindiv.dcmass_root - cropindiv.dcmass_leaf - cropindiv.dcmass_ho;		
	cropindiv.grs_cmass_agpool = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root - cropindiv.grs_cmass_leaf - cropindiv.grs_cmass_ho;
	cropindiv.ycmass_agpool = cropindiv.ycmass_plant - cropindiv.ycmass_root - cropindiv.ycmass_leaf - cropindiv.ycmass_ho;

	if(cropindiv.grs_cmass_agpool < 1.0e-9)
		cropindiv.grs_cmass_agpool = 0,0;
	if(cropindiv.ycmass_agpool < 1.0e-9)
		cropindiv.ycmass_agpool = 0,0;

	return;
}

/// Daily growth routine for crops
/** Allocates daily npp to leaf, roots and harvestable organs
 *  Requires updated value of fphu and hi.
 *  Equations are from Neitsch et al. 2002.
 */
void growth_crop_daily(Patch& patch) {

	if(date.day == 0)
		patch.nharv = 0;
	patch.isharvestday = false;
	double nharv_today = 0;

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
			cropindiv.ycmass_dead_leaf = 0.0;	
			cropindiv.ycmass_stem = 0.0;	

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

		cropindiv.dcmass_plant = 0.0;
		cropindiv.dcmass_leaf = 0.0;
		cropindiv.dcmass_root = 0.0;
		cropindiv.dcmass_ho = 0.0;
		cropindiv.dcmass_agpool = 0.0;
		cropindiv.dcmass_stem = 0.0;
		cropindiv.dnmass_ho = 0.0;

		// true crop allocation
		if(indiv.pft.phenology == CROPGREEN) {

			if(ppftcrop.growingseason) {

				double cmass_seed = 0.0;
				double nmass_seed = 0.0;

#ifdef DELAYED_SEEDCARBON
				// Seed carbon; portion the seed carbon over a 10-day period.
				if(dayinperiod(date.day, pppftcrop.sdate, (patchpft.cropphen->sdate + 9)) % 365 ) {			

					cmass_seed = 0.1 * CMASS_SEED;
					nmass_seed = 0.1 * CMASS_SEED / indiv.pft.cton_leaf_min;
				}
#else
				// add seed carbon on sowing date
				if(date.day == ppftcrop.sdate) {

					cmass_seed = CMASS_SEED;
					nmass_seed = CMASS_SEED / indiv.pft.cton_leaf_min;
				}
#endif

				if(ifnlim_lc[CROPLAND])
					allocation_crop_nlim(indiv, cmass_seed, nmass_seed);
				else
					allocation_crop(indiv, cmass_seed, nmass_seed);

				// save this year's maximum leaf carbon mass
				if(cropindiv.grs_cmass_leaf > cropindiv.cmass_leaf_max)	
					cropindiv.cmass_leaf_max = cropindiv.grs_cmass_leaf;

				// save leaf carbon mass at the beginning of senescence
				if(date.day == ppftcrop.sendate)
					cropindiv.cmass_leaf_sen = cropindiv.grs_cmass_leaf;

				// Check that no plant cmass or nmass is negative, if so, and correct fluxes
				double negative_cmass = indiv.check_C_mass();
				if(negative_cmass > 1.0e-14)
					dprintf("Year %d day %d Stand %d indiv %d: Negative main crop C mass in growth_crop_daily: %.15f\n", date.year, date.day, indiv.vegetation.patch.stand.id, indiv.id, -negative_cmass);
				double negative_nmass = indiv.check_N_mass();
				if(negative_nmass > 1.0e-14)
					dprintf("Year %d day %d Stand %d indiv %d: Negative main crop N mass in growth_crop_daily: %.15f\n", date.year, date.day, indiv.vegetation.patch.stand.id, indiv.id, -negative_nmass);
			}
			else if(date.day == ppftcrop.hdate) {

				patch.stand.isrotationday = true;

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;
				cropindiv.harv_cmass_stem += cropindiv.grs_cmass_stem;

				cropindiv.harv_nmass_root += indiv.nmass_root;
				cropindiv.harv_nmass_ho += cropindiv.nmass_ho;
				cropindiv.harv_nmass_leaf += indiv.nmass_leaf;
				cropindiv.harv_nmass_agpool += cropindiv.nmass_agpool;
				// dead_leaf to be addad

				if(ppftcrop.nharv == 1)
					cropindiv.cmass_ho_harvest[0] = cropindiv.grs_cmass_ho;
				else if(ppftcrop.nharv == 2)
					cropindiv.cmass_ho_harvest[1] = cropindiv.grs_cmass_ho;

				patch.isharvestday = true;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.get_gridcell().LC_updated && patchpft.cropphen->nharv == 1)
						scale_indiv(indiv, true);
					harvest_crop(indiv, indiv.pft, indiv.alive, indiv.cropindiv->isintercropgrass, true);
					patch.is_litter_day = true;
				}
				patch.nharv++;

				cropindiv.grs_cmass_plant = 0.0;
				cropindiv.grs_cmass_root = 0.0;
				cropindiv.grs_cmass_ho = 0.0;
				cropindiv.grs_cmass_leaf = 0.0;
				cropindiv.grs_cmass_agpool = 0.0;
				cropindiv.cmass_leaf_sen = 0.0;	

				cropindiv.grs_cmass_stem = 0.0;
				cropindiv.grs_cmass_dead_leaf = 0.0;
				cropindiv.nmass_dead_leaf = 0.0;
				cropindiv.nmass_agpool = 0.0;
				indiv.nmass_leaf = 0.0;
				indiv.nmass_root = 0.0;
				cropindiv.nmass_ho = 0.0;
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

					double cmass_seed = CMASS_SEED;
					cropindiv.grs_cmass_plant += cmass_seed;
					cropindiv.ycmass_plant += cmass_seed;
					cropindiv.dcmass_plant += cmass_seed;
					double nmass_seed = cmass_seed / indiv.pft.cton_leaf_min;
					indiv.nmass_leaf += nmass_seed / 2.0;
					indiv.nmass_root += nmass_seed / 2.0;

					indiv.report_flux(Fluxes::SEEDC, -cmass_seed);
					indiv.report_flux(Fluxes::SEEDN, -nmass_seed);

					indiv.last_turnover_day = -1;
				}
#endif
				cropindiv.dcmass_plant += indiv.dnpp;
				cropindiv.grs_cmass_plant += indiv.dnpp;		
				cropindiv.ycmass_plant += indiv.dnpp;

				indiv.ltor = indiv.wscal_mean * indiv.pft.ltor_max;

				// allocation to roots
				double froot = 1.0 / (1.0 + indiv.ltor);
				double grs_cmass_root_old = cropindiv.grs_cmass_root;	

				//Cumulative wscal-dependent root increase						
				cropindiv.grs_cmass_root = froot * cropindiv.grs_cmass_plant;
				cropindiv.dcmass_root = cropindiv.grs_cmass_root - grs_cmass_root_old;
				cropindiv.ycmass_root += cropindiv.dcmass_root;

				// allocation to leaves
				double fleaf = 1.0 - froot;
				double grs_cmass_leaf_old = cropindiv.grs_cmass_leaf;
				cropindiv.grs_cmass_leaf = cropindiv.grs_cmass_plant - cropindiv.grs_cmass_root;
				cropindiv.dcmass_leaf = cropindiv.grs_cmass_leaf - grs_cmass_leaf_old;
				cropindiv.ycmass_leaf += cropindiv.dcmass_leaf;

				// Check that no plant cmass is negative, if so, zero cmass and correct C fluxes
				double negative_cmass = indiv.check_C_mass();
				if(!SUPPRESSLARGEOUTPUT) {
					if(negative_cmass > 1.0e-14)
						dprintf("Year %d day %d Stand %d indiv %d: Negative intercrop C mass in growth_crop_daily: %.15f\n", date.year, date.day, indiv.vegetation.patch.stand.id, indiv.id, -negative_cmass);
					double negative_nmass = indiv.check_N_mass();
					if(negative_nmass > 1.0e-14)
						dprintf("Year %d day %d Stand %d indiv %d: Negative intercrop N mass in growth_crop_daily: %.15f\n", date.year, date.day, indiv.vegetation.patch.stand.id, indiv.id, -negative_nmass);
				}

				// save this year's maximum leaf carbon mass
				if(cropindiv.grs_cmass_leaf > cropindiv.cmass_leaf_max)	
					cropindiv.cmass_leaf_max = cropindiv.grs_cmass_leaf;
			}
			else if(date.day == patch.pft[patch.stand.pftid].get_cropphen()->eicdate) {

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;	
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;	
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;	
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;		
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;

				cropindiv.harv_nmass_root += indiv.nmass_root;
				cropindiv.harv_nmass_ho += cropindiv.nmass_ho;
				cropindiv.harv_nmass_leaf += indiv.nmass_leaf;
				cropindiv.harv_nmass_agpool += cropindiv.nmass_agpool;

				ppftcrop.nharv++;
				patch.isharvestday = true;
				nharv_today++;
				if(nharv_today > 1)	// In case of both C3 and C4 growing
					patch.nharv--;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.get_gridcell().LC_updated && patchpft.cropphen->nharv == 1)
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
				patch.nharv++;

				cropindiv.grs_cmass_plant = cropindiv.grs_cmass_root + cropindiv.grs_cmass_leaf;
			}
			
			if(indiv.continous_grass() && indiv.is_turnover_day()) {

				indiv.last_turnover_day = date.day;

				cropindiv.harv_cmass_plant += cropindiv.grs_cmass_plant;	
				cropindiv.harv_cmass_root += cropindiv.grs_cmass_root;	
				cropindiv.harv_cmass_leaf += cropindiv.grs_cmass_leaf;	
				cropindiv.harv_cmass_ho += cropindiv.grs_cmass_ho;		
				cropindiv.harv_cmass_agpool += cropindiv.grs_cmass_agpool;

				cropindiv.harv_nmass_root += indiv.nmass_root;
				cropindiv.harv_nmass_ho += cropindiv.nmass_ho;
				cropindiv.harv_nmass_leaf += indiv.nmass_leaf;
				cropindiv.harv_nmass_agpool += cropindiv.nmass_agpool;

				ppftcrop.nharv++;
				patch.isharvestday = true;
				nharv_today++;
				if(nharv_today > 1)	// In case of both C3 and C4 growing
					patch.nharv--;

				if(indiv.has_daily_turnover()) {
					if(patch.stand.get_gridcell().LC_updated && patchpft.cropphen->nharv == 1)
						scale_indiv(indiv, true);

					turnover_grass(indiv);
					patch.is_litter_day = true;
				}
				else {
					cropindiv.grs_cmass_root = 0.0;
					cropindiv.grs_cmass_ho = 0.0;
					cropindiv.grs_cmass_leaf = 0.0;
				}
				patch.nharv++;

				cropindiv.grs_cmass_plant = cropindiv.grs_cmass_root + cropindiv.grs_cmass_leaf;
				cropindiv.grs_cmass_agpool = 0.0;
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
void growth_daily(Patch& patch) {

	if(patch.stand.landcover == CROPLAND) {

		// allocate daily npp to leaf, roots and harvestable organs
		growth_crop_daily(patch);

		// update patchpft.lai_daily and fpc_daily
		lai_crop(patch);
	}
}


// Updates crop rotation status
/** Sets new crop management variables, typically on harvest day
 */
void crop_rotation(Stand& stand, int firsthistyear) {

	if(stand.landcover == CROPLAND) {

		CropRotation& rotation = stlist[stand.stid].rotation;
		bool postpone_rotation = false;

		stand.ndays_inrotation++;

		if(rotation.ncrops > 1 && stand.isrotationday) {

			int firstrotyear = rotation.firstrotyear + nyear_spinup - firsthistyear;

			// Alternative uses of firstrotyear:
/*			// 1. Before firstrotyear, grow only crop1:
			if(date.year < firstrotyear)
				postpone_rotation = true;
*/
			// 2. Synchronise rotation with firstrotyear:

			// A. At the creation of the stand:
			if(date.year < stand.first_year + 3)
			// B. At firstrotyear
//			if(date.year == firstrotyear - 1)
			// C. Continuously:
			{
				if((abs(firstrotyear - date.year) % rotation.ncrops) != stand.current_rot)
					postpone_rotation = true;
			}

			if(!postpone_rotation) {

				if(stand.infallow) {
					stand.infallow = false;
					stand.get_gridcell().pft[stand.pftid].sowing_restriction = false;
				}

				int old_pftid = stand.pftid;

				stand.rotate();

				for(unsigned int p=0; p<stand.nobj; p++) {

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
//				if((rotation.multicrop && stand.get_gridcell().pft[stand.pftid].multicrop) && rotation.ncrops == 2 && stand.current_rot == 1) {
				if(rotation.multicrop && rotation.ncrops == 2 && stand.current_rot == 1) {
					if(stand.pft[stand.pftid].sdate_force < 0)
						stand.pft[stand.pftid].sdate_force = stepfromdate(date.day, 10);
					if(stand.pft[stand.pftid].hdate_force < 0) {
						stand.pft[stand.pftid].hdate_force = stepfromdate(stand.pft[old_pftid].sdate_force, -10);
					}
				}

				if(stlist[stand.stid].management[stand.current_rot].fallow) {
					stand.infallow = true;
					stand.get_gridcell().pft[stand.pftid].sowing_restriction = true;
				}
			}

			stand.isrotationday = false;
		}
	}
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
void harvest_wood(Harvest_CN& i, Pft& pft, bool alive, double frac_cut, double harv_eff, double res_outtake_twig, double res_outtake_coarse_root) {

	double harvest = 0.0;
	double residue_outtake = 0.0;
	double stem_frac = 0.65;	// Temporary values, should be pft-specific
	double twig_frac = 0.13;
//	double coarse_root_frac = 0.22;
	double coarse_root_frac = 1.0 - stem_frac - twig_frac;
	double adhering_leaf_frac = 0.75;

	// only harvest trees
	if(pft.lifeform == GRASS)
		return;

	// all root carbon and nitrogen goes to litter
	if(alive) {

		i.litter_root += i.cmass_root * frac_cut;
		i.cmass_root *= (1.0 - frac_cut);
	}

	i.nmass_litter_root += i.nmass_root * frac_cut;
	i.nmass_litter_root += (i.nstore_labile + i.nstore_longterm) * frac_cut;
	i.nmass_root *= (1.0 - frac_cut);
	i.nstore_labile *= (1.0 - frac_cut);
	i.nstore_longterm *= (1.0 - frac_cut);

	if(alive) {	

		// Carbon:

		if (i.cmass_debt <= i.cmass_sap + i.cmass_heart) {

			// harvested stem wood
			harvest += harv_eff * stem_frac * (i.cmass_sap + i.cmass_heart - i.cmass_debt) * frac_cut;

			// harvested products not consumed (oxidised) this year put into harvested_products_slow
			if(ifslowharvestpool) {
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
		// debt larger than existing wood biomass
		}
		else {
			double debt_excess = i.cmass_debt - (i.cmass_sap + i.cmass_heart);
			dprintf("ATTENTION: cmass_debt > i.cmass_sap + i.cmass_heart; debt_excess=%f\n", debt_excess);
//			i.debt_excess += debt_excess * frac_cut;
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
		if(ifslowharvestpool) {
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
void harvest_wood(Individual& indiv, Pft& pft, bool alive, double frac_cut, double harv_eff, double res_outtake_twig, double res_outtake_coarse_root) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	harvest_wood(indiv_cp, pft, alive, frac_cut, harv_eff, res_outtake_twig, res_outtake_coarse_root);

	indiv_cp.copy_to_indiv(indiv);
}

void clearcut(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed) {

	Patch& patch = indiv.vegetation.patch;
	Patchpft& ppft = patch.pft[indiv.pft.id];

	if (indiv.pft.lifeform==TREE) {

		ppft.litter_sap += anpp;
		harvest_wood(indiv, indiv.pft, indiv.alive, 1.0, indiv.pft.harv_eff, indiv.pft.res_outtake); // frac_cut=1, harv_eff=pft.harv_eff, res_outtake_twig=pft.res_outtake, res_outtake_coarse_root=0
//		indiv.kill(true);
		indiv.vegetation.killobj();
		killed = true;
	}

//	patch.age=0;	//important for results
	patch.managed=true;
}

double forest_management(Patch& patch,bool age_class_run, int age_class) {

	Stand& stand = patch.stand;
	const double minbon=2.351; //The minimum average "bonitet" for a county in Sweden
	const double maxbon=11.311; //The maximum average "bonitet" for a county in Sweden
	const double bonitet = 10.0;	// Temporary atatic value

		// Code used for contineous forestry
	const int first_cutyear = nyear_spinup; //Simulation year when continues forestry harvesting starts
	int cut_int; //Interval between cuttings
	int patch_order; //Which year in a cutting interval the patch belongs to
	div_t cut_check;
//	cut_int=30-(int)(15.0*(stand.bonitet-minbon)/(maxbon-minbon));
	cut_int=30-(int)(15.0*(bonitet-minbon)/(maxbon-minbon));
	patch_order = (int)(patch.id * cut_int * 1.0 / (1.0 * stand.npatch()));
	cut_check = div(date.year - first_cutyear - patch_order, cut_int);

	if (date.year>=first_cutyear && cut_check.rem==0)
		return 0.40;
	else 
		return 0.00;

}

void harvest_forest(Individual& indiv, Pft& pft, bool alive, double anpp, bool& killed) {

	Patch& patch = indiv.vegetation.patch;
	Patchpft& ppft = patch.pft[indiv.pft.id];
	const double minbon=2.351; //The minimum average "bonitet" for a county in Sweden
	const double maxbon=11.311; //The maximum average "bonitet" for a county in Sweden
	const double bonitet = 10.0;	// Temporary static value
		
	int age_class = 0;
	double man_strength = 0.0;
	if(date.year > nyear_spinup && indiv.pft.lifeform == TREE)		
		man_strength = forest_management(patch, true, age_class);
	bool management_done=false;
		// Will tell the program to skip establishment and mortality if management has been
		// performed on this patch, Management add, FL 081127 (not implemented in this code yet ML, needs to be at patch-level)

	if (pft.lifeform==TREE && man_strength>0.00) {

		if (man_strength == 1.00) {
			clearcut(indiv, pft, alive, anpp, killed);
//			planting(stand,patch,pftlist);
		}
		else {

			double diam=pow(indiv.height/indiv.pft.k_allom2, 1.0/indiv.pft.k_allom3);
//			double diam_limit=0.13+0.07*(stand.bonitet-minbon)/(maxbon-minbon); //Harvest of trees > 13-20 cm
			double diam_limit=0.13+0.07*(bonitet-minbon)/(maxbon-minbon); //Harvest of trees > 13-20 cm
			double diam_max = diam_limit * 2.0;

			if (diam>diam_limit) {
				if(diam > diam_max)
					man_strength = 0.9;
				harvest_wood(indiv, pft, alive, man_strength, indiv.pft.harv_eff, indiv.pft.res_outtake); // frac_cut=man_strength, harv_eff=pft.harv_eff, res_outtake_twig=pft.res_outtake, res_outtake_coarse_root=0
				indiv.densindiv *= (1.0 - man_strength);
			}
		}
		management_done = true;	
		patch.managed = true;
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
		if ((i.nmass_leaf + i.nmass_agpool) > 0.0) {

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

void kill_remaining_vegetation(Individual& indiv, Pft& pft, bool alive, bool istruecrop_or_intercropgrass, bool burn) {

	Harvest_CN indiv_cp;

	indiv_cp.copy_from_indiv(indiv);

	kill_remaining_vegetation(indiv_cp, pft, alive, istruecrop_or_intercropgrass, burn);

	indiv_cp.copy_to_indiv(indiv);

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

void growth_crop_year(Individual& indiv, double& cmass_leaf_inc, double& cmass_root_inc, double& cmass_ho_inc, double& cmass_agpool_inc, double& cmass_stem_inc) {

	// true crop growth and grass intercrop growth; NB: bminit (cmass_repr & cmass_excess subtracted) not used !

	if(indiv.has_daily_turnover()) {

		indiv.cmass_leaf = 0.0;
		indiv.cmass_root = 0.0;
		indiv.cropindiv->cmass_ho = 0.0;
		indiv.cropindiv->cmass_agpool = 0.0;
		indiv.cropindiv->cmass_stem = 0.0;

		// Not completely accurate here when comparing this year's cmass after turnover with cmass increase (ycmass),
		// which could be from the preceding season, but probably OK, since values are not used for C balance.
		if(indiv.continous_grass()) {
			indiv.cmass_leaf = indiv.cmass_leaf_post_turnover;
			indiv.cmass_root = indiv.cmass_root_post_turnover;
		}
	}

	cmass_leaf_inc = indiv.cropindiv->ycmass_leaf + indiv.cropindiv->ycmass_dead_leaf;
	cmass_root_inc = indiv.cropindiv->ycmass_root;
	cmass_ho_inc = indiv.cropindiv->ycmass_ho;
	cmass_agpool_inc = indiv.cropindiv->ycmass_agpool;
	cmass_stem_inc = indiv.cropindiv->ycmass_stem;

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
			cropindiv.yield = cropindiv.ycmass_ho * indiv.pft.harv_eff * 2.0;// Should be /0.446 instead
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
