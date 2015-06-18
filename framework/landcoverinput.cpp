////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file landcoverinput.cpp
/// \brief Input code for land cover area fractions	from text files					
/// \author Mats Lindeskog
/// $Date: $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "landcoverinput.h"


// Switches to keep landcover area fractions at first historical year
bool fixedlu_histX = 0;
bool fixedcrop_histX = 0;

void read_gridlist(ListArray_id<Coord>& gridlist, const char* file_gridlist) {

	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	double dlon, dlat;
	bool eof = false;
	xtring descrip;

	// Read list of grid coordinates and store in Coord object 'gridlist'

	FILE* in_grid = fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input", (char*)file_gridlist);

	while (!eof) {
		
		// Read next record in file
		eof =! readfor(in_grid, "f,f,a#", &dlon, &dlat, &descrip);

		if (!eof && !(dlon == 0.0 && dlat == 0.0)) { // ignore blank lines at end (if any)
			Coord& c = gridlist.createobj(); // add new coordinate to grid list

			c.lon = dlon;
			c.lat = dlat;
			c.descrip = descrip;
		}
	}
	fclose(in_grid);
	gridlist.firstobj();
}

LandcoverInput::LandcoverInput()
	: nyears_cropland_ramp(0) {

	declare_parameter("minimizecftlist", &minimizecftlist, "Whether pfts not in crop fraction input file are removed from pftlist (0,1)");
	declare_parameter("nyears_cropland_ramp", &nyears_cropland_ramp, 1, 10000, "Number of years to increase cropland fraction linearly from 0 to first year's value");
	declare_parameter("frac_fixed_default_crops", &frac_fixed_default_crops, " whether to use all active crop stand types (0) or only stand types with suitable rainfed crops (based on crop pft tb and gridcell latitude) (1) when using fixed crop fractions");
}

void LandcoverInput::init() {

	if(!run_landcover)
		return;

	ListArray_id<Coord> gridlist;
	read_gridlist(gridlist, param["file_gridlist"].str);

	all_fracs_const=true;	//If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year

	//Retrieve file names for landcover files and open them if static values from ins-file are not used
	if (!lcfrac_fixed) {

		bool openLUfile = false;

		for(int i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				openLUfile = true;
		}

		if (openLUfile) {
			file_lu=param["file_lu"].str;

			// Open landcover area fraction file, return false if problem
			if(!LUdata.Open(file_lu, gridlist))
				fail("initio: could not open %s for input",(char*)file_lu);
			else {
				if(LUdata.GetFormat()==InData::LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data

				// Avoid large number of output files
				if(LUdata.GetNCells() > 50)
					printseparatestands = false;
			}
		}

		//Read LUC transitions
		if(gross_land_transfer == 2) {
			file_grossLUC=param["file_grossLUC"].str;
			if(!grossLUC.Open(file_grossLUC, gridlist))
				fail("initio: could not open %s for input",(char*)file_grossLUC);
		}
	}

	//Retrieve file names for crop fraction file and open them if static equal-size values are not used.
	if(run[CROPLAND] && !frac_fixed[CROPLAND])
	{
		file_lucrop=param["file_lucrop"].str;

		// Open crop fraction file, return false if problem
		if(!CFTdata.Open(file_lucrop, gridlist))
			fail("initio: could not open %s for input",(char*)file_lucrop);
		else {

			bool do_minimize = false;

			// Remove crop stand types from stlist that always have zero area fraction in all cells in grid list

			if(minimizecftlist && gridlist.nobj < 100) {	// Reduce the risk of accidentally using minimized cft lists when using split gridlists.

				CFTdata.CheckIfPresent(gridlist);
				do_minimize = true;
			}

			int n=0;
			stlist.firstobj();
			while(stlist.isobj) {	

				StandType& st = stlist.getobj();
				bool remove = false;

				if(st.landcover == CROPLAND) {
					if(do_minimize)
						remove = !CFTdata.item_has_data(st.name);
					else
						remove = !CFTdata.item_in_header(st.name);
				}

				if(remove) {
					n+=1;
					stlist.killobj();
					nst--;
				}
				else {
					st.id-=n;
					stlist.nextobj();
				}			
			}
		}

		if(CFTdata.GetFormat()==InData::LOCAL_YEARLY) 
			all_fracs_const=false;				// Set all_fracs_const to false if yearly data

		}

	// Remove pft:s from pftlist that are not grown in simulated stand types
	int n = 0;
	pftlist.firstobj();
	while(pftlist.isobj) {	

		bool remove = false;
		Pft& pft = pftlist.getobj();

		if(pft.landcover == CROPLAND) {

			bool found = false;

			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();

				if(st.pftinrotation(pft.name) >= 0) {
					found = true;
					break;
				}
				stlist.nextobj();
			}

			if(!found)
				remove = true;
		}

		if(remove && !(pft.isintercropgrass && ifintercropgrass)) {
			n += 1;
			pftlist.killobj();
			npft--;
		}
		else {
			pft.id -= n;
			pftlist.nextobj();
		}			
	}
	gridlist.killall();
}

bool LandcoverInput::loadlandcover(Coord c) {

	bool LUerror = false;

	if (!lcfrac_fixed) {

		// Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
		// transferred to gridcell.landcoverfrac each year in getlandcover()

		bool loadLU = false;

		for(int i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				loadLU = true;
		}

		if (loadLU) {
			// Load landcover area fraction data from input file to data object
			if (!LUdata.Load(c)) {
				dprintf("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror = true;		// skip this stand
			}
		}

		//Read LUC transitions
		if(gross_land_transfer == 2 && !LUerror) {

			if(!grossLUC.Load(c)) {
				dprintf("Problems with gross LUC transitions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror = true;	// skip this stand
			}
		}
	}

	if(run[CROPLAND] && !LUerror)
	{
		if(!frac_fixed[CROPLAND]) {
			
			// Crop fraction data: read from crop fraction file; dynamic, so data for all years are loaded to CFTdata object and 
			// transferred to gridcell.cftfrac each year in getlandcover()			

			if(!CFTdata.Load(c)) {
				dprintf("Problems with CFT fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror = true;	// skip this stand
			}
		}
	}

	return LUerror;
}

void LandcoverInput::getlandcover(Gridcell& gridcell) {

	int i, year_saved;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;
	Landcover& lc = gridcell.landcover;

	// Calender year of start of simulation (after spinup)
	int first_historic_year = date.first_calendar_year + nyear_spinup;

	// Use values for first historic year during spinup period, unless data exist before firsthistyear
	// Use values for last historic year during the time after that
	// This is handled by the text input class.
	int year = date.get_calendar_year();

	if(fixedlu_histX) {
		year_saved = year;
		year = first_historic_year;
	}

	//Save old fraction values:									
	for(int i=0; i<NLANDCOVERTYPES; i++)
		lc.frac_old[i] = lc.frac[i];
	for(unsigned int i = 0; i < gridcell.st.nobj; i++) {
		Gridcellst& gcst = gridcell.st[i];
		gcst.frac_old = gcst.frac;
	}

	if(lcfrac_fixed) {		// If landcover area fractions are set in the ins-file.
		if(date.year == 0) {	// Year 0: called by landcover_init
	
			int nactive_landcovertypes = 0;

			for(i=0; i<NLANDCOVERTYPES; i++) {
				if(run[i])
					nactive_landcovertypes++;
			}

			for(i=0;i<NLANDCOVERTYPES;i++) {
				lc.frac[i] = 1.0 * run[i] / (double)nactive_landcovertypes;	// only set fractions that are active
				sum_active += gridcell.landcover.frac[i];
				sum_tot = sum_active;
			}		
		}
	}
	else {	// landcover area fractions are read from input file(s)	

		bool printyear = year >= LUdata.GetFirstyear() && LUdata.GetFirstyear() >= 0;	
		bool getLU = false;

		for(i=0; i<NLANDCOVERTYPES; i++) {
			if(run[i] && i != NATURAL)
				getLU = true;
		}

		if(getLU) {	

			if(LUdata.Get(year, 0) < 0.0) {		// Missing data (negative values)
				if(date.year == 1)
					dprintf("Missing landcover fraction data for year %d, natural vegetation fraction set to 1.0\n", year);
				memset(lc.frac, 0, sizeof(double) * NLANDCOVERTYPES);
				gridcell.landcover.frac[NATURAL] = 1.0;
				sum_active = 1.0;
			}
			else {

				if(year == LUdata.GetFirstyear() + LUdata.GetnYears() + 1)
					dprintf("Last year of landcover fraction data used from year %d and onwards\n", year);

				for(i=0; i<NLANDCOVERTYPES; i++) {

					if(run[i]) {
						double lcfrac = 0.0;
						switch(i)
						{
						case URBAN:
							lcfrac = LUdata.Get(year,"URBAN");
							break;
						case CROPLAND:
							lcfrac = LUdata.Get(year,"CROPLAND");
							break;
						case PASTURE:
							lcfrac = LUdata.Get(year,"PASTURE");
							break;
						case FOREST:
							lcfrac = LUdata.Get(year,"FOREST");
							break;
						case NATURAL:
							lcfrac = LUdata.Get(year,"NATURAL");
							break;
						case PEATLAND:
	//						lcfrac = LUdata.Get(year,"PEATLAND");
							break;
						case BARREN:
							lcfrac = LUdata.Get(year,"BARREN");
							break;
						default:
							if(date.year == 0)
								dprintf("Modify code to deal with landcover input!\n");
						}
							
						if(lcfrac == NOTFOUND)	// land cover not found in input file
							lcfrac = 0.0;
						else if(lc.frac[i] < 0.0 || gridcell.landcover.frac[i] > 1.0) {	// discard unreasonable values	
							if(printyear)
								dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
							lcfrac = 0.0;
						}
						lc.frac[i] = lcfrac;
						sum_tot += lcfrac;
						sum_active += run[i] * gridcell.landcover.frac[i];
					}
				}

				if (grassforcrop) {
					lc.frac[PASTURE]+=lc.frac[CROPLAND];
					lc.frac[CROPLAND]=0.0;
				}
				if(sum_tot != 1.0 && sum_tot != 0.0) {		// Check input data, rescale if sum !=1.0	

					sum_active = 0.0;		// reset sum of active landcover fractions

					if(sum_tot < 0.99 || sum_tot > 1.01) {
						if(printyear) {
							dprintf("WARNING ! landcover fraction sum is %4.2f for input year %d\n", sum_tot, year);
							dprintf("Rescaling landcover fractions year %d !\n", date.get_calendar_year());
						}
					}

					for(i=0; i<NLANDCOVERTYPES; i++) {
						lc.frac[i] /= sum_tot;
						sum_active += lc.frac[i];
					}
				}
			}
		}
		else {
			lc.frac[NATURAL] = 0.0;
		}

		// NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
		if(fabs(sum_active - 1.0) > 1.0e-14)	{	// if landcover types are turned off in the instruction file, or if more landcover types are added in other input files, can be either less or more than 1.0

			if(date.year == 0)
				dprintf("Landcover fraction sum not 1.0 !\n");

			if(run[NATURAL]) {	// Transfer landcover areas not simulated to NATURAL fraction, if simulated.		
				if(date.year == 0) {
					if(sum_active < 1.0)
						dprintf("Inactive fractions (%4.3f) transferred to NATURAL fraction.\n", 1.0-sum_active);
					else
						dprintf("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).\n", sum_active-1.0, lc.frac[NATURAL]);
				}

				lc.frac[NATURAL] += (1.0 - sum_active);	// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					dprintf("New NATURAL fraction is %4.3f.\n", gridcell.landcover.frac[NATURAL]);

				sum_active = 1.0;		// sum_active should now be 1.0

				if(lc.frac[NATURAL] < 0.0) {	// If new landcover type fraction is bigger than the natural fraction (something wrong in the distribution of input file area fractions)						
					if(date.year == 0)
						dprintf("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.\n");

					sum_active -= lc.frac[NATURAL];	// fraction not possible to transfer moved back to sum_active, which will now be >1.0 again
					lc.frac[NATURAL] = 0.0;

					for(i=0; i<NLANDCOVERTYPES; i++) {
						lc.frac[i] /= sum_active;		// fraction rescaled to unity sum
						if(run[i])
							if(date.year == 0)
								dprintf("Landcover type %d fraction is %4.3f\n", i, gridcell.landcover.frac[i]);
					}
				}
			}
			else {
//				if(date.year == 0)
//					dprintf("Rescaling landcover fractions !\n");
//				for(i=0; i<NLANDCOVERTYPES; i++)
//					lc.frac[i] /= sum_active;					// 1) if NATURAL not simulated, rescale active fractions to 1.0
				if(date.year == 0)
					dprintf("Non-unity fraction sum retained.\n");				// 2) let sum remain non-unity
			}
		}

		if(nyears_cropland_ramp) {

			bool doramp = false;
			int firstyear;
			if(LUdata.GetFirstyear() >= 0) {
				if(year < LUdata.GetFirstyear()) {
					doramp = true;
					firstyear = LUdata.GetFirstyear();
				}
			}
			else {			
				if(year < first_historic_year) {		
					doramp = true;
					firstyear = first_historic_year;
				}
			}

			if(doramp) {
				int first_reduction_year = first_historic_year - nyear_spinup + (int)(SOLVESOMCENT_SPINEND * (nyear_spinup - freenyears) + freenyears) + 1;
				int max_ramp_years = firstyear - first_reduction_year;
				if(nyears_cropland_ramp > max_ramp_years)
					dprintf("Requested cropland ramp period too long for given nyear_spinup. Maximum is %d.\n", max_ramp_years);

				double reduce_cropland = min((double)(firstyear - year) / min(nyears_cropland_ramp, max_ramp_years), 1.0) * lc.frac[CROPLAND];
				lc.frac[CROPLAND] -= reduce_cropland;
				lc.frac[NATURAL] += reduce_cropland;
//				if(lc.frac[CROPLAND])
//					dprintf("Year %d Cropland frac=%.3f\n", year, lc.frac[CROPLAND]);
			}
		}
	}

	// Set fractions for static stand types
	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		Gridcellst& gcst = gridcell.st[st.id];

		gcst.frac = 0.0;
		if(nst_lc[st.landcover] == 1)
			gcst.frac = 1.0;
		else if(frac_fixed[st.landcover] && !frac_fixed_default_crops)
			gcst.frac = 1.0 / (double)nst_lc[st.landcover];

		stlist.nextobj();
	}

	if(run[CROPLAND] && (!frac_fixed[CROPLAND] || frac_fixed[CROPLAND] && frac_fixed_default_crops)) {

		sum=0.0;

		bool printyear = year >= CFTdata.GetFirstyear() && CFTdata.GetFirstyear() >= 0;

		// Get fractions for dynamic stand types within a land cover
		if(!(frac_fixed[CROPLAND] && frac_fixed_default_crops)) {

			if(fixedcrop_histX)
				year=first_historic_year;
			else if(fixedlu_histX)
				year=year_saved;

			if(year == CFTdata.GetFirstyear() + CFTdata.GetnYears())
				dprintf("Last year of cropland fraction data used from year %d and onwards\n", year);

			// sum fractions for active crop pft:s and discard unreasonable values
			for(i=0; i<nst; i++) {
				if(stlist[i].landcover == CROPLAND)	{

					double cropfrac = CFTdata.Get(year,stlist[i].name);

					if(cropfrac == NOTFOUND)	// crop not found in input file
						cropfrac = 0.0;
					else if(cropfrac < 0.0 || cropfrac > 1.0)	{	// discard unreasonable values
						if(!(!gridcell.landcover.frac[CROPLAND] && cropfrac < 0.0) && printyear)
							dprintf("WARNING ! crop fraction size out of limits, set to 0.0\n");
						cropfrac = 0.0;
					}
					sum += gridcell.st[i].frac = cropfrac;
				}
			}

			if(printyear) {
				if(lc.frac[CROPLAND]==0.0) {
					if(sum!=0.0) {
	//					dprintf("WARNING ! crop landcover fraction is 0.0 for year %d while crop data exist !\n", year);
					}
				}
				else if(sum==0.0) {
					dprintf("WARNING ! crop fraction sum is 0.0 for year %d while LU[CROPLAND] is > 0 !\n", year);
				}
			}
		}

		// If crop fraction data are missing or if fixed default crops are used, try to find suitable stand types to use.
		if(sum==0.0) {

			if(lc.frac[CROPLAND]) {
				// Use equal areas of rainfed stand types with tropical or temperate crop pft:s based on base temperatures
				int nsts = 0;
				stlist.firstobj();
				while(stlist.isobj) {
					StandType& st = stlist.getobj();
					Gridcellst& gcst = gridcell.st[st.id];
					if(st.landcover == CROPLAND) {
						if(st.management[0].hydrology == RAINFED) {
							if(gridcell.get_lat() > 30 || gridcell.get_lat() < -30) {
								if(pftlist[pftlist.getpftid(st.management[0].pftname)].tb <= 5) {
									gcst.frac = 1.0;
									nsts++;
								}
							}
							else {
								if(pftlist[pftlist.getpftid(st.management[0].pftname)].tb > 5) {
									gcst.frac = 1.0;
									nsts++;
								}
							}
						}
					}
					stlist.nextobj();
				}
				if(nsts) {
					if(lc.frac[CROPLAND] && printyear && !(frac_fixed[CROPLAND] && frac_fixed_default_crops))
						dprintf("Missing crop fraction data. Using available suitable stand types for year %d\n", year);
				}
				else {
					fail("No suitable default stand types available for filling missing data. Check crop fraction input data\n");
				}
				stlist.firstobj();
				while(stlist.isobj) {
					StandType& st = stlist.getobj();
					Gridcellst& gcst = gridcell.st[st.id];
					if(st.landcover == CROPLAND)
						gcst.frac /= nsts;
					stlist.nextobj();
				}
			}
		}
		else {
			// rescale active crop fraction so sum is 1.0

			stlist.firstobj();
			while(stlist.isobj) {
				StandType& st = stlist.getobj();
				Gridcellst& gcst = gridcell.st[st.id];
				if(st.landcover == CROPLAND)
					gcst.frac /= sum;
				stlist.nextobj();
			}
			if((sum < 0.99 || sum > 1.01) && printyear) {	// warn if sum is significantly different from 1.0 
				dprintf("WARNING ! crop fraction sum is %5.3f for input year %d\n", sum, year);
				dprintf("Rescaling crop fractions year %d !\n", date.get_calendar_year());
			}
		}
	}

	// Update landcover frac_change arrays
	for(int i=0; i<NLANDCOVERTYPES; i++) {
		lc.frac_change[i] = lc.frac[i] - lc.frac_old[i];
	}

	// Convert fractions from landcover-based to gridcell-based and update gridcellst frac_change arrays
	stlist.firstobj();
	while (stlist.isobj) {
		StandType& st = stlist.getobj();
		Gridcellst& gcst = gridcell.st[st.id];

		gcst.frac = gcst.frac * gridcell.landcover.frac[st.landcover];
		if(fabs(gcst.frac_old - gcst.frac) < 1.0e-14)
			gcst.frac = gcst.frac_old;
		gcst.frac_change = gcst.frac - gcst.frac_old ;
		stlist.nextobj();
	}
}

bool LandcoverInput::get_land_transitions(Gridcell& gridcell) {

	bool result = false;
	if(gross_land_transfer == 3) {

		// Read stand type transfer fractions from file here and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.
		// return get_st_transfer(gridcell);
		dprintf("Currently no code for option gross_land_transfer==3\n");
	}
	else if(gross_land_transfer == 2) {

		// Read landcover transfer fractions from file here and put them into the st_frac_transfer array.
		// Landcover and stand type net fractions still need to be read from file as previously.
		result = get_lc_transfer(gridcell);
	}
	return result;
}


/// Read LUC transitions
bool LandcoverInput::get_lc_transfer(Gridcell& gridcell) {

	int year;
	double tot_frac_ch = 0.0;
	const bool print_adjustment_info = false;
	Landcover& lc = gridcell.landcover;

	if(!grossLUC.isloaded() || date.year < nyear_spinup + 1)
		return false;

	//If have not reached second year of simulation (after spin-up) then gross_lc_change_frac must be zero (no land-use change in spin-up).

	year = date.get_calendar_year() - 1;

	// Assume that transitions in file are correct at end of year, therefore want to get 
	// "last year's" transitions, as landcover_dynamics is called at the beginning of the year.
	// Transfers from primary (v) and secondary (s) land preferentially reduces the oldest and the 
	// youngest stands, respectively. Transitions from primary to secondary NATURAL land result 
	// in killing of vegetation and creating a new NATURAL stand.

	const bool primary_to_secondary = false;

	lc.frac_transfer[CROPLAND][PASTURE] += grossLUC.Get(year,"cp");
	lc.frac_transfer[PASTURE][CROPLAND] += grossLUC.Get(year,"pc");
	lc.frac_transfer[PASTURE][NATURAL] += grossLUC.Get(year,"pv");
	lc.frac_transfer[NATURAL][PASTURE] += grossLUC.Get(year,"vp");
	lc.frac_transfer[NATURAL][CROPLAND] += grossLUC.Get(year,"vc");
	lc.frac_transfer[CROPLAND][NATURAL] += grossLUC.Get(year,"cv");
	lc.frac_transfer[NATURAL][CROPLAND] += grossLUC.Get(year,"sc");
	lc.frac_transfer[CROPLAND][NATURAL] += grossLUC.Get(year,"cs");
	lc.frac_transfer[NATURAL][PASTURE] += grossLUC.Get(year,"sp");
	lc.frac_transfer[PASTURE][NATURAL] += grossLUC.Get(year,"ps");

	if(ifprimary_lc_transfer) {
		lc.primary_frac_transfer[NATURAL][PASTURE] += grossLUC.Get(year,"vp");
		lc.primary_frac_transfer[NATURAL][CROPLAND] += grossLUC.Get(year,"vc");
		if(primary_to_secondary) {
			lc.frac_transfer[NATURAL][NATURAL] += grossLUC.Get(year,"vs");
			lc.primary_frac_transfer[NATURAL][NATURAL] += grossLUC.Get(year,"vs");
		}
	}

	// Check if gross lcc input data are consistent with net lcc input file. Try to adjust if not.
	bool error = false;
	double net_lc_change[NLANDCOVERTYPES] = {0.0};
	double gross_lc_increase[NLANDCOVERTYPES] = {0.0};
	double gross_lc_decrease[NLANDCOVERTYPES] = {0.0};
	double pos_error = 0.0;
	double neg_error = 0.0;

	for(int from=0; from<NLANDCOVERTYPES; from++) {

		if(run[from]) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if(run[to]) {
					
					net_lc_change[from] -= lc.frac_transfer[from][to];
					net_lc_change[from] += lc.frac_transfer[to][from];
					gross_lc_decrease[from] += lc.frac_transfer[from][to];
					gross_lc_increase[from] += lc.frac_transfer[to][from];
					tot_frac_ch += lc.frac_transfer[to][from];
				}
			}
			if(fabs(lc.frac_change[from] - net_lc_change[from])  > 1.0e-14) {
				error = true;
				if(print_adjustment_info) {
					dprintf("\nIn get_lc_transfer: Year %d: lc_change_array sum not equal to landcoverfrac_change value for landcover %d\n", date.year, from);
					dprintf("dif=%.15f\n", net_lc_change[from] - lc.frac_change[from]);
					dprintf("lc_change_array sum=%.15f\n", net_lc_change[from]);
					dprintf("landcoverfrac_change=%.15f\n", lc.frac_change[from]);
				}
			}

		}
	}

	if(error) {

		int n_adjust_from[NLANDCOVERTYPES] = {0};
		int n_adjust_to[NLANDCOVERTYPES] = {0};
		double partition_adjustment[NLANDCOVERTYPES][NLANDCOVERTYPES];
		double prim_sec_ratio[NLANDCOVERTYPES][NLANDCOVERTYPES];
		double original_error[NLANDCOVERTYPES] = {0.0};

		memset(partition_adjustment, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(double));
		memset(prim_sec_ratio, 0, NLANDCOVERTYPES * NLANDCOVERTYPES * sizeof(double));

		for(int from=0; from<NLANDCOVERTYPES; from++) {

			for(int to=0; to<NLANDCOVERTYPES; to++) {

				if(lc.frac_transfer[from][to])
					prim_sec_ratio[from][to] = lc.primary_frac_transfer[from][to] / lc.frac_transfer[from][to];

				if((lc.frac_transfer[from][to] + lc.frac_transfer[to][from]) > 0.0) {

					if(gridcell.landcover.frac_old[from] - gross_lc_decrease[from] < 0.0 || gridcell.landcover.frac_old[from] + gross_lc_increase[from] > 1.0) {
						partition_adjustment[from][to] = 1.0;
						partition_adjustment[to][from] = 0.0;
					}
					else if(gridcell.landcover.frac_old[to] - gross_lc_decrease[to] < 0.0 || gridcell.landcover.frac_old[to] + gross_lc_increase[to] > 1.0) {
						partition_adjustment[from][to] = 0.0;
						partition_adjustment[to][from] = 1.0;
					}
					else {
						partition_adjustment[from][to] = lc.frac_transfer[from][to] / (lc.frac_transfer[from][to] + lc.frac_transfer[to][from]);
					}
				}
			}

			if(run[from])
				original_error[from] = net_lc_change[from] - lc.frac_change[from];

			if(original_error[from] > 0.0)
				pos_error += original_error[from];
			else if(original_error[from] < 0.0)
				neg_error += original_error[from];
		}

		if(fabs(pos_error + neg_error) > 1.0e-15)
			fail("\nYear %d: pos_error + neg_error = %.15f\n\n", date.year, pos_error + neg_error);

		// Determine number of possible adjustments to and from a lc
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(fabs(original_error[from]) > 1.0e-14) {

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(fabs(original_error[to]) > 1.0e-14) {

						// Errors must have opposite signs
						if(fabs(original_error[from] + original_error[to]) - (fabs(original_error[from]) + fabs(original_error[to])) < -1.0e-14) {

							if((lc.frac_transfer[from][to] + lc.frac_transfer[to][from]) > 0.0) {

								if(lc.frac_transfer[from][to])
									n_adjust_from[from]++;
								if(lc.frac_transfer[to][from])
									n_adjust_to[to]++;
							}
							else {

								for(int third=0; third<NLANDCOVERTYPES; third++) {

									if(lc.frac_transfer[from][third] + lc.frac_transfer[third][from] > 0.0 && lc.frac_transfer[to][third] + lc.frac_transfer[third][to]) {
										n_adjust_from[from]++;
										n_adjust_to[to]++;
									}
								}
							}
						}
					}
				}
			}
		}

		// Adjust transfers for lc:s with transition input rounding errors
		for(int from=0; from<NLANDCOVERTYPES; from++) {

			if(run[from] && fabs(lc.frac_change[from] - net_lc_change[from]) > 1.0e-14) {

				for(int to=0; to<NLANDCOVERTYPES; to++) {

					if(run[to] && fabs(lc.frac_change[to] - net_lc_change[to])  > 1.0e-14 && from != to) {

						// Errors must have opposite signs
						if(fabs(net_lc_change[from] - lc.frac_change[from] + net_lc_change[to] - lc.frac_change[to]) - (fabs(net_lc_change[from] - lc.frac_change[from]) + fabs(net_lc_change[to] - lc.frac_change[to])) < -1.0e-14)	{

							// Correct transfer between two lc:s
							if((lc.frac_transfer[from][to] + lc.frac_transfer[to][from]) > 0.0) {

								if(print_adjustment_info) {
									dprintf("Before: transfer lc %d to %d: %.15f\n", from, to, lc.frac_transfer[from][to]);
									dprintf("Before: transfer lc %d to %d: %.15f\n", to, from, lc.frac_transfer[to][from]);
								}
								if(n_adjust_from[from]) {
									lc.frac_transfer[from][to] += (net_lc_change[from] - lc.frac_change[from]) / (double)n_adjust_from[from]	* partition_adjustment[from][to];
									n_adjust_from[from]--;
								}
								if(n_adjust_to[to]) {
									lc.frac_transfer[to][from] += (net_lc_change[to] - lc.frac_change[to]) / (double)n_adjust_to[to]	* partition_adjustment[to][from];
									n_adjust_to[to]--;
								}
								if(print_adjustment_info) {
									dprintf("After: transfer lc %d to %d: %.15f\n", from, to, lc.frac_transfer[from][to]);
									dprintf("After: transfer lc %d to %d: %.15f\n", to, from, lc.frac_transfer[to][from]);
								}
							}
							// Third lc used
							else {

								for(int third=0; third<NLANDCOVERTYPES; third++) {

									if(lc.frac_transfer[from][third] + lc.frac_transfer[third][from] > 0.0 && lc.frac_transfer[to][third] + lc.frac_transfer[third][to]) {

										if(print_adjustment_info) {
											dprintf("Before: transfer lc %d to %d: %.15f\n", from, third, lc.frac_transfer[from][third]);
											dprintf("Before: transfer lc %d to %d: %.15f\n", third, from, lc.frac_transfer[third][from]);
											dprintf("Before: transfer lc %d to %d: %.15f\n", to, third, lc.frac_transfer[to][third]);
											dprintf("Before: transfer lc %d to %d: %.15f\n", third, to, lc.frac_transfer[third][to]);
										}

										lc.frac_transfer[from][third] += (net_lc_change[from] - lc.frac_change[from]) * partition_adjustment[from][third];
										lc.frac_transfer[third][from] -= (net_lc_change[from] - lc.frac_change[from]) * partition_adjustment[third][from];
										lc.frac_transfer[to][third] += (net_lc_change[to] - lc.frac_change[to]) * partition_adjustment[to][third];
										lc.frac_transfer[third][to] -= (net_lc_change[to] - lc.frac_change[to]) * partition_adjustment[third][to];

										if(print_adjustment_info) {
											dprintf("After: transfer lc %d to %d: %.15f\n", from, third, lc.frac_transfer[from][third]);
											dprintf("After: transfer lc %d to %d: %.15f\n", third, from, lc.frac_transfer[third][from]);
											dprintf("After: transfer lc %d to %d: %.15f\n", to, third, lc.frac_transfer[to][third]);
											dprintf("After: transfer lc %d to %d: %.15f\n", third, to, lc.frac_transfer[third][to]);
										}
									}
								}
							}
						}
					}
				}

				// Update net_lc_change-values:
				for(int from=0; from<NLANDCOVERTYPES; from++) {

					if(run[from]) {

						net_lc_change[from] = 0.0;

						for(int to=0; to<NLANDCOVERTYPES; to++) {

							if(run[to]) {
											
								net_lc_change[from] -= lc.frac_transfer[from][to];
								net_lc_change[from] += lc.frac_transfer[to][from];
							}
						}
					}
				}
			}
		}

		// Adjust primary land fractions
		for(int from=0; from<NLANDCOVERTYPES; from++) {
			for(int to=0; to<NLANDCOVERTYPES; to++)
				lc.primary_frac_transfer[from][to] = prim_sec_ratio[from][to] * lc.frac_transfer[from][to];
		}
	}

	if(tot_frac_ch > 1.0e-14)
		return true;
	else
		return false;
}

ManagementInput::ManagementInput() {
}

void ManagementInput::init() {

	if(!run_landcover)
		return;

	ListArray_id<Coord> gridlist;
	read_gridlist(gridlist, param["file_gridlist"].str);

	if(run_landcover && run[CROPLAND]) {

		file_sdates=param["file_sdates"].str;
		if(file_sdates != "")	{
			if(!sdates.Open(file_sdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_sdates);
			readsowingdates = true;
		}

		file_hdates=param["file_hdates"].str;
		if(file_hdates != "")	{
			if(!hdates.Open(file_hdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_hdates);
			readharvestdates = true;
		}

		file_Nfert=param["file_Nfert"].str;
		if(file_Nfert != "")	{
			if(!Nfert.Open(file_Nfert, gridlist))
				fail("initio: could not open %s for input",(char*)file_Nfert);
			readNfert = true;
		}
	}
	gridlist.killall();
}

bool ManagementInput::loadmanagement(Coord c) {

	bool LUerror = false;

	if(readsowingdates) { 
		if(!sdates.Load(c)) {
			dprintf("Problems with sowing date input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n", c.lon, c.lat);
			LUerror = true;	// skip this stand
		}
	}
	if(readharvestdates && !LUerror) {
		if(!hdates.Load(c)) {
			dprintf("Problems with harvest date input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n", c.lon, c.lat);
			LUerror = true;	// skip this stand
		}
	}
	if(readNfert && !LUerror) {
		if(!Nfert.Load(c)) {
//				dprintf("Problems with N fertilization input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n", c.lon, c.lat);
//				LUerror = true;	// skip this stand
			dprintf("N fertilization data not found in input file for %.2f,%.2f.\n\n", c.lon, c.lat);
		}
	}
	return LUerror;
}

void ManagementInput::getsowingdates(Gridcell& gridcell) {

	if(!sdates.isloaded())
		return;

	int year = date.get_calendar_year();

	for(int i=0; i<npft; i++) {
		if(pftlist[i].phenology == CROPGREEN) {	

			gridcell.pft[i].sdate_force = (int)sdates.Get(year,pftlist[i].name);

			// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
			for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
				Standpft& standpft = gridcell[j].pft[i];
				if(standpft.active)
					standpft.sdate_force = gridcell.pft[i].sdate_force;
			}
		}
	}
}

void ManagementInput::getharvestdates(Gridcell& gridcell) {

	if(!hdates.isloaded())
		return;

	int year = date.get_calendar_year();

	for(int i=0; i<npft; i++)	{
		if(pftlist[i].phenology == CROPGREEN) {		

			gridcell.pft[pftlist[i].id].hdate_force = (int)hdates.Get(year,pftlist[i].name);

			// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
			for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
				Standpft& standpft = gridcell[j].pft[i];
				if(standpft.active)
					standpft.hdate_force = gridcell.pft[i].hdate_force;
			}
		}
	}
}

void ManagementInput::getNfert(Gridcell& gridcell) {

	if(!Nfert.isloaded())
		return;

	int year = date.get_calendar_year();

	for(int i=0; i<npft; i++)	{
		if(pftlist[i].phenology == CROPGREEN) {		
			gridcell.pft[pftlist[i].id].Nfert_read = Nfert.Get(year,pftlist[i].name);
		}
	}
}

void ManagementInput::getmanagement(Gridcell& gridcell) {

	if(run[CROPLAND]) {

		//Read sowing dates from input file, put into gridcellpft.sdate_force
		if(readsowingdates)		
			getsowingdates(gridcell);
		//Read harvest dates from input file, put into gridcellpft.hdate_force
		if(readharvestdates)		
			getharvestdates(gridcell);
		//Read N fertilization from input file, put into gridcellpft.Nfert_read
		if(readNfert)		
			getNfert(gridcell);
	}
}