///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.cpp
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// Landcover change.
///
/// \author Mats Lindeskog,
/// \Part of code in this file as well as in cropphenology.cpp, cropallocation.cpp and 
/// \management.cpp based on LPJ-mL C++ code received from Alberte Bondeau in 2008.
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "landcover.h"
#include "guessmath.h"

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

/// Help function to access two-dimentional arrays that have been created dynamically
int index(int from, int to, int ncols = nst) {

	return from * ncols + to;
}


/// Creation of stands when run_landcover==true
void landcover_init(Gridcell& gridcell, InputModule* input_module) {

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
	}

	// get landcover and crop area fractions from landcover input file(s) or ins-file.
	LandcoverInput *landcover_input = input_module->get_landcover_input();
	landcover_input->getlandcover(gridcell);

	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		Gridcellst& gcst = gridcell.st[st.id];

		gcst.frac_old = gcst.frac;

		if(gcst.frac > 0.0) {
			gridcell.create_stand_lu(st, gcst.frac);
		}

		stlist.nextobj();
	}
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

	for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
		Gridcellst& gcst = gridcell.st[i];
		gcst.nstands = 0;
	}

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
		gridcell.st[stand.stid].nstands++;
	}


	for(int from=0; from<nst; from++) {

		StandType& st = stlist[from];
		Gridcellst& gcst = gridcell.st[from];

		if(gcst.gross_frac_decrease > 0.0) {

			if(gcst.nstands > 1) {

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
										int first_year = max(stand.first_year, stand.clone_year);

										// Don't reduce stands younger than the age limit, unless this is the last stand in the loop or the initial stand has been killed
//										if(gross_land_transfer && date.year - first_year < age_limit_reduce && count_st != st.nstands && !reduce_all_stands && !restart) continue;
										if(date.year - first_year < age_limit_reduce && count_st != gcst.nstands && !reduce_all_stands && !restart) continue;
							
										// convert equal percentage of areas from all stands
										if(reduce_all_stands) {
											stand.frac_change += st_change_remain * stand.get_gridcell_fraction() / gcst.frac;
											stand.gross_frac_decrease -= st_change_remain * stand.get_gridcell_fraction() / gcst.frac;
											stand.transfer_area_st[to] = -st_change_remain * stand.get_gridcell_fraction() / gcst.frac;
											stand.set_gridcell_fraction(stand.get_gridcell_fraction() + st_change_remain * stand.get_gridcell_fraction() / gcst.frac);
										}
										else {
											if(stand.get_gridcell_fraction() > 0.0) {

												//all natural landcover decrease is taken from this stand
												if(stand.get_gridcell_fraction() >= -st_change_remain) {
													stand.frac_change += st_change_remain;
													stand.gross_frac_decrease -= st_change_remain;
													stand.transfer_area_st[to] -= st_change_remain;
													stand.set_gridcell_fraction(stand.get_gridcell_fraction() + st_change_remain);

													if(stand.get_gridcell_fraction() < 1.0e-15 || gcst.frac == 0.0 && stand.get_gridcell_fraction() < 1.0e-14)
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
			else if(gcst.nstands == 1 ) {

				for(unsigned int i = 0; i < gridcell.size(); i++) {

					Stand& stand = gridcell[i];

					if(stand.stid == gcst.id) {

						stand.frac_change = -gcst.gross_frac_decrease;
						stand.gross_frac_decrease = gcst.gross_frac_decrease;
						stand.set_gridcell_fraction(stand.get_gridcell_fraction() + stand.frac_change);

						if(stand.get_gridcell_fraction() < 1.0e-15 || gcst.frac == 0.0 && stand.get_gridcell_fraction() < 1.0e-14)
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
		Gridcellst& gcst = gridcell.st[i];
		landcovertype lc = st.landcover;
		bool expand_to_new_stand = gridcell.landcover.expand_to_new_stand[lc];

		if(gcst.gross_frac_increase > 0.0) {	// Not cloned stands

			if(!expand_to_new_stand) {

				for(unsigned int i = 0; i < gridcell.size(); i++) {

					Stand& stand = gridcell[i];

					if(stand.stid == st.id) {

						// Identify which stands to expand (not secondary stands)
 
						if(!stand.cloned) {

							stand.frac_change += gcst.gross_frac_increase;
							stand.gross_frac_increase = gcst.gross_frac_increase;
							stand.set_gridcell_fraction(stand.get_gridcell_fraction() + gcst.gross_frac_increase);
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
 *  NB. New land cover types must be included in the preference arrays !
 *  Also, PEATLAND, URBAN and BARREN needs to be included when using dynamic fractions for these land cover types.
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
		Gridcellst& gcst = gridcell.st[i];

		if(gcst.frac_change > 0.0) {
			net_receptor_remain[i] = gcst.frac_change;
			nreceptor_st++;
		}
		else if(gcst.frac_change < 0.0){
			net_donor_remain[i] = -gcst.frac_change;
			ndonor_st++;
		}

		if(recip_lc_change[st.landcover] > 0.0) {
			recip_receptor_remain[i] = recip_lc_change[st.landcover] * gcst.frac_old / gridcell.landcover.frac_old[st.landcover];
			recip_donor_remain[i] = recip_lc_change[st.landcover] * gcst.frac_old / gridcell.landcover.frac_old[st.landcover];
		}
	}


	// Net land cover change

	// Simplest cases: no ambiguities
	if(ndonor_st == 1 || nreceptor_st == 1) {

		for(int from=0; from<nst; from++) {

			StandType& st_donor = stlist[from];
			Gridcellst& gcst_donor = gridcell.st[from];

			if(gcst_donor.frac_change < 0.0) {

				for(int to=0; to<nst; to++) {

					StandType& st_receptor = stlist[to];
					Gridcellst& gcst_receptor = gridcell.st[to];
		
					if(ndonor_st == 1) {

						if(gcst_receptor.frac_change > 0.0)
							st_frac_transfer[index(from, to)] =  gcst_receptor.frac_change;
					}
					else if(nreceptor_st == 1) {

						if(gcst_receptor.frac_change > 0.0)
							st_frac_transfer[index(from, to)] = - gcst_donor.frac_change;
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
				Gridcellst& gcst = gridcell.st[from];
				if(st.landcover == i) {
					abs_frac_change_sum += fabs(gcst.frac_change);
					frac_change_sum += gcst.frac_change;
					if(gcst.frac_change < 0.0)
						donor_sum -= gcst.frac_change;
					else if(gcst.frac_change > 0.0)
						receptor_sum += gcst.frac_change;
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
		Gridcellst& gcst_donor = gridcell.st[from];

		if(recip_donor_remain[from] > 1.0e-14) {

			for(int to=0; to<nst; to++) {

				StandType& st_receptor = stlist[to];
				Gridcellst& gcst_receptor = gridcell.st[to];

				if(recip_transfer_remain[st_donor.landcover][st_receptor.landcover]  > 1.0e-14 
					&& recip_receptor_remain[to] > 1.0e-14 && recip_donor_remain[from] > 1.0e-14) {

					if(equal_distribution) {
						st_frac_transfer[index(from, to)] += recip_lc_frac_transfer[st_donor.landcover][st_receptor.landcover] * gcst_donor.frac_old / gridcell.landcover.frac_old[st_donor.landcover] * gcst_receptor.frac_old / gridcell.landcover.frac_old[st_receptor.landcover];
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
			donor_area = stand.transfer_area_lc((landcovertype)landcover_receptor);
		else
			donor_area = stand.gross_frac_decrease;

		if(donor_area > 0.0 && (stand.id == standid || standid < 0) && (stand.landcover == landcover_donor || landcover_donor < 0) && (stand.stid == stid_donor || stid_donor < 0)) {

			double zero_gridcell_luc_cflux = gridcell.landcover.acflux_landuse_change;
			double ccont_stand_pre_orig = 0.0;
			double zero_to_ccont_pre = to.ccont();
			double to_ccont_pre = 0.0;
			double zero_gridcell_luc_nflux = gridcell.landcover.anflux_landuse_change;
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

					//Sum added litter C & N:
					to.transfer_litter_leaf[indiv.pft.id] += cp.litter_leaf * scale;
					to.transfer_litter_root[indiv.pft.id] += cp.litter_root * scale;
					to.transfer_litter_sap[indiv.pft.id] += cp.litter_sap * scale;
					to.transfer_litter_heart[indiv.pft.id] += cp.litter_heart * scale;

					to.transfer_nmass_litter_leaf[indiv.pft.id] += cp.nmass_litter_leaf * scale;
					to.transfer_nmass_litter_root[indiv.pft.id] += cp.nmass_litter_root * scale;
					to.transfer_nmass_litter_sap[indiv.pft.id] += cp.nmass_litter_sap * scale;
					to.transfer_nmass_litter_heart[indiv.pft.id] += cp.nmass_litter_heart * scale;

					gridcell.landcover.acflux_landuse_change += cp.acflux_harvest * donor_area / (double)stand.nobj;
					gridcell.landcover.acflux_landuse_change_lc[stand.landcover] += cp.acflux_harvest * donor_area / (double)stand.nobj;
					gridcell.landcover.anflux_landuse_change += cp.anflux_harvest * donor_area / (double)stand.nobj;
					gridcell.landcover.anflux_landuse_change_lc[stand.landcover] += cp.anflux_harvest * donor_area / (double)stand.nobj;

//					gridcell.acflux_landuse_change += -cp.debt_excess * donor_area / (double)stand.nobj;

					if(ifslowharvestpool) {
						to.transfer_harvested_products_slow[indiv.pft.id] += cp.harvested_products_slow * scale;
						to.transfer_harvested_products_slow_nmass[indiv.pft.id] += cp.harvested_products_slow_nmass * scale;
					}

					vegetation.nextobj();
				}

				stand.nextobj();
			}

			ccont_pre_orig_tot += ccont_stand_pre_orig * donor_area / receiving_fraction;
			dif_tot += dif;
			acflux_landuse_change_tot += (gridcell.landcover.acflux_landuse_change - zero_gridcell_luc_cflux) / receiving_fraction;
			ncont_pre_orig_tot += ncont_stand_pre_orig * donor_area / receiving_fraction;
			ndif_tot += ndif;
			anflux_landuse_change_tot += (gridcell.landcover.anflux_landuse_change - zero_gridcell_luc_nflux) / receiving_fraction;

		}
	}

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
		Gridcellst& gcst = gridcell.st[st.id];
		landcovertype lc = st.landcover;

		bool expand_to_new_stand = gridcell.landcover.expand_to_new_stand[lc];

		if(gcst.gross_frac_increase || gcst.gross_frac_decrease || gcst.frac == 0.0) {
			// first stand created
			if(gcst.frac_old == 0.0 && gcst.frac > 0.0)
				Stand& stand = gridcell.create_stand_lu(st, gcst.frac);

			// last stand killed
			else if(gcst.frac_old > 0.0 && gcst.frac == 0.0) {
				Gridcell::iterator gc_itr = gridcell.begin();
				while (gc_itr != gridcell.end()) {
					Stand& stand = *gc_itr;
					if(stand.stid == st.id)
						gc_itr = gridcell.delete_stand(gc_itr);
					else
						++gc_itr;
				}
			}
			else {
				if(expand_to_new_stand) {
					// new stand created from other landcover types (pooled option, not created in transfer_to_new_stand())
					if(gcst.gross_frac_increase > 0.0) {
						// Fewer patches in secondary stands if npatch_secondarystand < npatch:
						Stand& stand = gridcell.create_stand_lu(st, gcst.gross_frac_increase, npatch_secondarystand);
					}
				}
				// secondary stand killed if all of its area converted to other landcover types
				if(gcst.nstands > 1 && gcst.gross_frac_decrease > 0.0) {
					Gridcell::iterator gc_itr = gridcell.begin();
					while (gc_itr != gridcell.end()) {
						Stand& stand = *gc_itr;
						if(stand.stid == st.id && stand.get_gridcell_fraction() == 0)
							gc_itr = gridcell.delete_stand(gc_itr);
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
				gridcell.landcover.LC_updated = true;
			}
		}
		++gc_itr;
	}
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
	else if(landcover_donor == PASTURE) {
//		if(landcover_receptor == NATURAL)
//			copy_type = CLONESTAND;
	}

	return copy_type;
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
						dprintf("WARNING: cloning natural stand without allowing natural pft:s to grow in the new stand. Was this intended ?\n");

					landcover_change_transfer transfer;

					if(copy_type == CLONESTAND_KILLTREES) {
						donor_stand_change(gridcell, transfer_area, transfer, -1,-1, stid_receptor, -1, stand.id, false);
						receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand.id);
					}
	

					new_stand.cloned = true;
					new_stand.gross_frac_increase = 0.0;
					new_stand.frac_change = 0.0;
					new_stand.cloned_fraction = transfer_area;
					new_stand.origin = stand.landcover;

					new_stand.firstobj();
					while(new_stand.isobj) {
						Patch& patch = new_stand.getobj();
						Vegetation& vegetation = patch.vegetation;
						vegetation.firstobj();
						while(vegetation.isobj) {

							Individual& indiv = vegetation.getobj();
							Standpft& standpft = new_stand.pft[indiv.pft.id];

							if(!standpft.active) {
								// Treatment of pft individuals not allowed to grow anymore in the new stand:
								// Clearcut
								harvest_wood(indiv, 1.0, 1.0, 0.95, 0.9, true);	// frac_cut=1, harv_eff=1, res_outtake_twig=0.95, res_outtake_coarse_root=0.9
								// Grass killed, C+N goes to soil
								kill_remaining_vegetation(indiv);
								vegetation.killobj();
							}
							else
								vegetation.nextobj();
						}
						new_stand.nextobj();
					}
				}
				else if(copy_type == NEWSTAND_KILLALL) {

					landcover_change_transfer transfer;

					donor_stand_change(gridcell, transfer_area, transfer, -1,-1, stid_receptor, -1, stand.id);
					Stand& new_stand = gridcell.create_stand_lu(stlist[stid_receptor], transfer_area, npatch_secondarystand);
					receiving_stand_change(gridcell, transfer, true, -1, -1, 1.0, new_stand.id);

					new_stand.cloned = true;
					new_stand.gross_frac_increase = 0.0;
					new_stand.frac_change = 0.0;
					new_stand.cloned_fraction = transfer_area;
				}

				cloned_area += transfer_area;
				stand.cloned_fraction -= transfer_area;
				gridcell.st[stid_donor].gross_frac_decrease -= transfer_area;
				gridcell.st[stid_donor].frac_change += transfer_area;
				stand.gross_frac_decrease -= transfer_area;
				stand.frac_change += transfer_area;
				stand.transfer_area_st[stid_receptor] -= transfer_area;
				gridcell.st[stid_receptor].gross_frac_increase -= transfer_area;
				gridcell.st[stid_receptor].frac_change -= transfer_area;

				if(gridcell.st[stid_donor].gross_frac_decrease < 1.0e-15)
					gridcell.st[stid_donor].gross_frac_decrease = 0.0;
				if(stand.gross_frac_decrease < 1.0e-15)
					stand.gross_frac_decrease = 0.0;
				if(stand.transfer_area_st[stid_receptor] < 1.0e-15)
					stand.transfer_area_st[stid_receptor] = 0.0;
				if(gridcell.st[stid_receptor].gross_frac_increase < 1.0e-15)
					gridcell.st[stid_receptor].gross_frac_increase = 0.0;
			}
		}
	}

	return cloned_area;
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

			lc_frac_transfer[from][to] += gross_lc_change_frac[from][to] * min(gridcell.landcover.frac_old[from], gridcell.landcover.frac_old[to]);
		}
	}
}


/// Simulates gross landcover change by adding a specified fraction to the st_frac_transfer array
void simulate_gross_st_transfer(Gridcell& gridcell, double* st_frac_transfer) {

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
		Gridcellst& gcst_donor = gridcell.st[from];

		for(int to=0; to<nst; to++) {

			StandType& st_receptor = stlist[to];
			Gridcellst& gcst_receptor = gridcell.st[to];

			if(gross_lc_change_frac[st_donor.landcover][st_receptor.landcover]) {

				// All stand types with an area:
				st_frac_transfer[index(from, to)] += gross_lc_change_frac[st_donor.landcover][st_receptor.landcover] * min(min(gcst_donor.frac_old, gcst_donor.frac), min(gcst_receptor.frac_old, gcst_receptor.frac));
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
		lc_frac_sum += gridcell.landcover.frac[i];

	if(fabs(lc_frac_sum - 1.0)  > 1.0e-14) {
		dprintf("\nCheck 1: Year %d: landcover fraction sum: %.15f", date.year, lc_frac_sum);
		error = true;
	}
	/////////

	// Check that stand type sum is 1:
	double st_frac_sum = 0.0;
	for(int i=0; i<nst; i++)
		st_frac_sum += gridcell.st[i].frac;
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
		Gridcellst& gcst = gridcell.st[i];

		if(fabs(test_st_change[i] - gcst.frac_change) > 1.0e-14) {
			dprintf("\nCheck 4: Year %d: st_change_array sum not equal to st.frac_change value for stand type %d\n", date.year, i);
			dprintf("dif=%.15f", fabs(test_st_change[i] - gcst.frac_change));
			error = true;
		}
	}
	if(test_st_change)
		delete[] test_st_change;
	////////////

	// Test if st_change_array is compatible with st.frac_old/frac values
	for(int from=0; from<nst; from++) {

		StandType& st_donor = stlist[from];
		Gridcellst& gcst_donor = gridcell.st[from];

		for(int to=0; to<nst; to++) {

			StandType& st_receptor = stlist[to];
			Gridcellst& gcst_receptor = gridcell.st[to];
			if(st_change_array[index(from, to)] > (gcst_donor.frac_old + 1.0e-14) || st_change_array[index(from, to)] > (gcst_receptor.frac + 1.0e-14)) {
				dprintf("\nCheck 5: Year %d: st_change_array sum not compatible with st.frac_old/frac values for stand types %d and %d\n", date.year, from, to);
				dprintf("st_change_array=%.15f, st_donor.frac_old=%.15f, st_receptor.frac=%.15f", st_change_array[index(from, to)], gcst_donor.frac_old, gcst_receptor.frac);
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
		Gridcellst& gcst = gridcell.st[s];

		if(gcst.frac_change < 0.0){

			double stands_frac_sum = 0.0;

			for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
				Stand& stand = gridcell[i];

				if(stand.stid == st.id)
					stands_frac_sum += stand.get_gridcell_fraction();
			}

			if(-gcst.frac_change - stands_frac_sum > 1.0e-14) {
				dprintf("\nCheck 8: Year %d: stand type %d fraction reduction bigger than sum of stands\n", date.year, s);
				dprintf("dif=%.15f", -gcst.frac_change - stands_frac_sum);
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
		Gridcellst& gcst = gridcell.st[s];

		if(gcst.frac == 0.0) {

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
		Gridcellst& gcst = gridcell.st[s];

		double stands_frac_sum = 0.0;

		for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
			Stand& stand = gridcell[i];

			if(stand.stid == st.id)
				stands_frac_sum += stand.get_gridcell_fraction();
		}

		if(gcst.frac_change < 0.0) {
			if(fabs(gcst.frac - stands_frac_sum - gcst.gross_frac_increase) > 1.0e-14) {
				dprintf("\nCheck 12: Year %d: fraction sum of stands not equal to stand type value for stand type %d\n", date.year, s);
				dprintf("dif=%.15f", fabs(gcst.frac - stands_frac_sum - gcst.gross_frac_increase));
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
		Gridcellst& gcst = gridcell.st[s];

		double stands_frac_sum = 0.0;

		for(unsigned int i = 0; i < gridcell.nbr_stands(); i++) {
			Stand& stand = gridcell[i];

			if(stand.stid == st.id)
				stands_frac_sum += stand.get_gridcell_fraction();
		}

		if(gcst.frac_change >= 0.0)
		if(fabs(gcst.frac - stands_frac_sum) > 1.0e-14) {
			dprintf("\nCheck 13: Year %d: fraction sum of stands not equal to stand type value for stand type %d\n", date.year, s);
			dprintf("dif=%.15f", fabs(gcst.frac - stands_frac_sum));
			error = true;
		}
	}
	if(error)
		dprintf("\n\n");

	return error;
}

/// Gets this year's landcover and crop area fractions, landcover transitions and checks that area changes are significant.
/** Stores changes in area fractions for the stand types
 *
 *  OUTPUT PARAMETERS
 *  \param landcoverfrac_change				array with this year's difference in area fractions of the different landcovers
 *  \param lc_frac_transfer					array with this year's transitions in area fractions between the different landcovers
 *  \param st_frac_transfer					array with this year's transitions in area fractions between the different landcovers
 *  \param primary_lc_frac_transfer			array with this year's transitions in area fractions from primary landcovers
 *  \param primary_st_frac_transfer			array with this year's transitions in area fractions from primary stand types
 *  \param LCchangeCtransfer				whether to transfer carbon, nitrogen and water of reduced stands to expanding stands
 */
bool checkLCchange(Gridcell& gridcell, double landcoverfrac_change[NLANDCOVERTYPES], double lc_frac_transfer[][NLANDCOVERTYPES], 
				   double* st_frac_transfer, double primary_lc_frac_transfer[][NLANDCOVERTYPES], double* primary_st_frac_transfer, 
				   bool& LCchangeCtransfer, InputModule* input_module) {

	double cropfrac_sum_old = 0.0;
	double change_stand = 0.0;
	double changeLC = 0.0;
	double change_crop = 0.0;
	double transferred_fraction = 0.0;
	double receiving_fraction = 0.0;
	bool change = true, gross_LCC = false;
	LandcoverInput *landcover_input = input_module->get_landcover_input();


	//Save old fraction values:									
	for(int i=0; i<NLANDCOVERTYPES; i++)
		gridcell.landcover.frac_old[i] = gridcell.landcover.frac[i];
	for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
		Gridcellst& gcst = gridcell.st[i];

		gcst.frac_old = gcst.frac;
		if(gcst.st.landcover == CROPLAND)
			cropfrac_sum_old += gcst.frac_old;
	}

	//Get new gridcell.landcoverfrac and/or standtype.frac from LUdata and CFTdata.
	landcover_input->getlandcover(gridcell);	

	for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
		Gridcellst& gcst = gridcell.st[i];
		gcst.frac_change = gcst.frac - gcst.frac_old ;
	}

	if(!lcfrac_fixed) {
		for(int i=0; i<NLANDCOVERTYPES; i++) {
			landcoverfrac_change[i] = gridcell.landcover.frac[i] - gridcell.landcover.frac_old[i];
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
		for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
			Gridcellst& gcst = gridcell.st[i];

			if(gcst.st.landcover == CROPLAND) {

				double stfrac_change = gcst.frac_change;

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

	if(gross_land_transfer == 3) {

		// Read stand type transfer fractions from file here and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.

		// input_module->get_st_transfer();
		dprintf("Currently no code for option gross_land_transfer==3\n");
	}
	else if(gross_land_transfer == 2) {

		// Read landcover transfer fractions from file here and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.

		if(landcover_input->get_lc_transfer(gridcell, landcoverfrac_change, lc_frac_transfer, primary_lc_frac_transfer)) {
			gross_LCC = false;
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
			simulate_gross_st_transfer(gridcell, st_frac_transfer);
	}

	check_fractions(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer);
	check_fractions1(gridcell);

	// if no changes, do nothing.
	if(changeLC < 1.0e-15 && change_crop < 1.0e-15 && !gross_LCC) {
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
				dprintf("Transferred landcover fractions not balanced !\nLandcover change carbon flux not calculated.\n");
			}
		}
		change = true;
	}

	return change;
}

/// Updates all landcover, stand type and stand area fractions each year, possibly resulting in the creation and killing of stands.
/** Harvests transferred areas and transfers litter etc. of reduced stands to expanding stands and harvested matter to fluxes
 *  and (in the case of wood) to long-lived pools.
 *
 *  The instruction file parameters lcfrac_fixed and cftfrac_fixed will determine if land fractions are read from input files 
 *  (gridcell static or dynamic) or from the instruction file (global static) and if any land cover change is possible.
 *
 *  The instruction file parameter gross_land_transfer determines whether gross land transfer fractions are to be read from an input  
 *  file (2,3) or simulated by an amount set in simulate_gross_lc_transfer() or simulate_gross_st_transfer() (1) or be turned off (0).
 *
 *  Land cover fraction files need to be read together with gross land transition files. Net and gross land cover change is 
 *  expected to be compatible, although rounding errors are handled by the input code.
 *
 *  If instruction file parameter iftransfer_to_new_stand is true, new stands may be created for separate land transfers  
 *  in transfer_to_new_stand(), according to rules in copy_stand_type().
 *
 *  Otherwise, new transferred areas are pooled according to instruction file parameter transfer_level (0: one big pool; 1: land cover-level; 
 *  2: stand type-level) and rules in the gridcell arrays pool_from_all_landcovers[to] and pool_to_all_landcovers[from], set
 *  in the Gridcell constructor.
 *  New stands may then be created in stand_dynamics() from the pooled land if the gridcell array expand_to_new_stand[lc] value for the 
 *  receiving land cover, set in the Gridcell constructor, is true (natural and forest).
 *  If not, soil and litter carbon and nitrogen as well as water is pooled with the receiving stands (cropland, pasture).
 *
 *  Transferred land with living plant mass after land cover change should be put in new stands, using the stand cloning mode in 
 *  transfer_to_new_stand(). Stands with living plant mass should not be pooled with newly transferred land, unconditionally so 
 *  if they contain trees (pasture can be expanded without too much problems). New stands should instead be created, either in stand_dynamics()
 *  or transfer_to_new_stand() as described above.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module) {

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

	gridcell.landcover.LC_updated = false;

	for(unsigned int i=0; i<gridcell.nbr_stands(); ++i)
		gridcell[i].scale_LC_change = 1.0;

	// get new landcover and stand type area fractions from input files, set transition arrays
	if(!all_fracs_const) {
		// this call returns 0, causing this function to return, if no significant landcover changes this year, 
		// sets LCchangeCtransfer to 0 if unbalanced landcover changes (if some landcovers are inactivated), thus inactivating transfer of C and N
		if(!checkLCchange(gridcell, landcoverfrac_change, lc_frac_transfer, st_frac_transfer, primary_lc_frac_transfer, 
				primary_st_frac_transfer, LCchangeCtransfer, input_module)) {
			delete[] st_frac_transfer;
			delete[] primary_st_frac_transfer;
			return;
		}
	}

	for(int i=0; i<nst; i++) {

		gridcell.st[i].gross_frac_increase = 0.0;
		gridcell.st[i].gross_frac_decrease = 0.0;
	}

	for(int from=0; from<nst; from++) {

		for(int to=0; to<nst; to++) {

			if(st_frac_transfer[index(from, to)] > 0.0) {
				gridcell.st[to].gross_frac_increase += st_frac_transfer[index(from, to)];
				gridcell.st[from].gross_frac_decrease += st_frac_transfer[index(from, to)];
			}
		}
	}

	double ccont_tot_1 = gridcell.ccont();
	double cflux_tot_1 = gridcell.cflux();
	double ncont_tot_1 = gridcell.ncont();
	double nflux_tot_1 = gridcell.nflux();

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

			if(gridcell.landcover.pool_to_all_landcovers[from]) { // alt.c
				if(gross_landcoverfrac_decrease[from] > 0.0) {
					donor_stand_change(gridcell, gross_landcoverfrac_decrease[from], transfer_lc_from[from], -1, from);
				}
			}
		}

		for(int to=0; to<NLANDCOVERTYPES; to++) {

			double receiving_fraction = gross_landcoverfrac_increase[to];

			for(int from=0; from<NLANDCOVERTYPES; from++) {

				double receiving_fraction_2d = gross_landcoverfrac_transfer[from][to];

				if(receiving_fraction_2d > 0.0) {

					if(gridcell.landcover.pool_from_all_landcovers[to]) {
						if(!gridcell.landcover.pool_to_all_landcovers[from]) {	// alt.a
							// Add from->to transfer to to-pool:
							donor_stand_change(gridcell, receiving_fraction, transfer_lc[to], to, from);
						}
						else {	// alt.a+c
							// Add from-pool value to to-pool:
							transfer_lc[to].add(transfer_lc_from[from], receiving_fraction_2d / receiving_fraction);
						}
					}
					else {
						if(!gridcell.landcover.pool_to_all_landcovers[from]) {	// alt.b
							// Add from->to transfer to [from][to] place in 2d-array:
							donor_stand_change(gridcell, receiving_fraction_2d, transfer_lc_2d[from][to], to, from);
						}
						else {	// alt.c
							// Copy from-pool value to [from][to] place in 2d-array:
							transfer_lc_2d[from][to].copy(transfer_lc_from[from]);
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

			if(gridcell.landcover.pool_from_all_landcovers[to]) {
				// Alt.a: All donor lc:s are pooled into receiving lc:s
				if(gross_landcoverfrac_increase[to] > 0.0) {
					receiving_stand_change(gridcell, transfer_lc[to], LCchangeCtransfer, to);
				}
			}
			else {

				for(int from=0; from<NLANDCOVERTYPES; from++) {

					if(gross_landcoverfrac_transfer[from][to] > 0.0) {
						if(!gridcell.landcover.pool_to_all_landcovers[from]) {
							// Alt.b: Transfers from donors are independent:
							receiving_stand_change(gridcell, transfer_lc_2d[from][to], LCchangeCtransfer, to, -1, gross_landcoverfrac_transfer[from][to] / gross_landcoverfrac_increase[to]);
						}
						else {
							// Alt.c: Donor lc:s are pooled before transfer to recipients:
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
			Gridcellst& gcst = gridcell.st[from];

			// Pooling of donor stands within a stand type
			if(gridcell.landcover.pool_to_all_landcovers[st.landcover] || gcst.nstands == 1) { // alt.c
				if(gcst.gross_frac_decrease > 0.0) {
					donor_stand_change(gridcell, gcst.gross_frac_decrease, transfer_st_from[from], -1, -1, -1, from);
				}
			}
		}

		for(int to=0;to<nst;to++) {

			StandType& st_to = stlist[to];
			Gridcellst& gcst_to = gridcell.st[to];
			double receiving_fraction = gcst_to.gross_frac_increase;

			for(int from=0; from<nst; from++) {

				StandType& st_from = stlist[from];
				Gridcellst& gcst_from = gridcell.st[from];
				double receiving_fraction_2d = st_frac_transfer[index(from, to)];

				if(receiving_fraction_2d > 0.0) {

					if(gridcell.landcover.pool_from_all_landcovers[st_to.landcover]) {
						if(!(gridcell.landcover.pool_to_all_landcovers[st_from.landcover] || gcst_from.nstands == 1)) {	// alt.a
							// Add from->to transfer to to-pool:
							donor_stand_change(gridcell, receiving_fraction, transfer_st[to], -1, -1, to, from);
						}
						else {	// alt.a+c
							// Add from-pool value to to-pool:
							transfer_st[to].add(transfer_st_from[from], receiving_fraction_2d / receiving_fraction);
						}
					}
					else {
						if(!(gridcell.landcover.pool_to_all_landcovers[st_from.landcover] || gcst_from.nstands == 1)) {	// alt.b
							// Add from->to transfer to [from][to] place in 2d-array:
							donor_stand_change(gridcell, receiving_fraction_2d, transfer_st_2d[index(from, to)], -1, -1, to, from);
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
		for(int to=0; to<nst; to++) {

			StandType& st = stlist[to];
			Gridcellst& gcst = gridcell.st[to];

			if(gridcell.landcover.pool_from_all_landcovers[st.landcover]) {
				// Alt.a: All donor st:s are pooled into receiving st:s
				if(gcst.gross_frac_increase > 0.0) {
					receiving_stand_change(gridcell, transfer_st[to], LCchangeCtransfer, -1, to);
				}
			}
			else {

				for(int from=0; from<nst; from++) {

					if(st_frac_transfer[index(from, to)] > 0.0) {
						StandType& st_from = stlist[from];
						Gridcellst& gcst_from = gridcell.st[from];
						if(!(gridcell.landcover.pool_to_all_landcovers[st_from.landcover] || gcst_from.nstands == 1)) {
							// Alt.b: Transfers between donor and receptor st:s are independent:
							receiving_stand_change(gridcell, transfer_st_2d[index(from, to)], LCchangeCtransfer, st.landcover, to, st_frac_transfer[index(from, to)] / gcst.gross_frac_increase);
						}
						else {
							// Alt.c: Donor stands in a stand type are pooled before transfer to recipients:
							receiving_stand_change(gridcell, transfer_st_from[from], LCchangeCtransfer, st.landcover, to, st_frac_transfer[index(from, to)] / gcst.gross_frac_increase);
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

	double ccont_tot = 0.0;
	double cflux_tot = gridcell.cflux();

	for(unsigned int i=0; i<gridcell.size(); i++) {
		Stand& stand = gridcell[i];

		double ccont_stand = stand.ccont(stand.scale_LC_change);
		ccont_tot += ccont_stand * stand.get_gridcell_fraction();
	}

	if(fabs(cflux_tot - cflux_tot_1 + ccont_tot - ccont_tot_1) > 1.0e-12)
		dprintf("WARNING ! C balance after lcc off\n");

	double ncont_tot = 0.0;
	double nflux_tot = 0.0;

	for(unsigned int i=0; i<gridcell.size(); i++) {
		Stand& stand = gridcell[i];

		double ncont_stand = stand.ncont(stand.scale_LC_change);
		ncont_tot += ncont_stand * stand.get_gridcell_fraction();
		nflux_tot += stand.nflux() * stand.get_gridcell_fraction();
	}

	nflux_tot_1 = nflux_tot;	// Disregard fluxes not already reset this year and the associated scaling problems 
	nflux_tot += gridcell.landcover.anflux_landuse_change + gridcell.landcover.anflux_harvest_slow;

	if(fabs(nflux_tot - nflux_tot_1 + ncont_tot - ncont_tot_1) > 1.0e-12)
		dprintf("WARNING ! N balance after lcc off\n");

	if(st_frac_transfer)
		delete[] st_frac_transfer;
	if(primary_st_frac_transfer)
		delete[] primary_st_frac_transfer;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of landcover_change_transfer member functions
////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Bondeau A, Smith PC, Zaehle S, Schaphoff S, Lucht W, Cramer W, Gerten D, Lotze-Campen H,
//   Müller C, Reichstein M & Smith B 2007. Modelling the role of agriculture for the 
//   20th century global terrestrial carbon balance. Global Change Biology, 13:679-706.
// Lindeskog M, Arneth A, Bondeau A, Waha K, Seaquist J, Olin S, & Smith B 2013. 
//   Implications of accounting for land use in simulations of ecosystem services and  
//   carbon cycling in Africa. Earth Syst Dynam Discuss 4:235-278.
