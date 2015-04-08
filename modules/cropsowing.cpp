////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cropsowing.h
/// \brief Seasonality and sowing date calculations				
/// \author Mats Lindeskog, Stefan Olin
/// $Date:  $
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "landcover.h"
#include "cropsowing.h"

#define SD_TEMP_WINDOW				//Uses sowing window for temperature-dependent sowing.
#define IRRIGATED_USE_TEMP_SDATE	//Use temperature-dependent sowing date for irrigated crops at site with PRECTEMP seasonality.
//#define LOW_SOWING_TEMPERATURE_LIMIT	//Sowing not allowed when temperature is always below sowing limit. Intercrop grass groen instead year through.
#define HIGH_SOWING_TEMPERATURE_LIMIT	//Sowing not allowed when mean temperature is above limit (TeWW).

/// Autumn sowing types for crops
enum {NOFORCING, AUTUMNSOWING, SPRINGSOWING};

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
/** Called from crop_sowing_gridcell() once a year (June 30(180) in the north, December 31(364) in the south)
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

	if(climate.lat > -15.0 && climate.lat < 20.0 && climate.gridcell.get_lon() > 90.0)
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
		climate.mtemp20[m] = climate.hmtemp_20[m].lastadd();
		climate.mprec20[m] = climate.hmprec_20[m].lastadd();
		climate.aprec += climate.hmprec_20[m].lastadd();
		climate.mpet_year[m] = climate.hmeet_20[m].lastadd()*PRIESTLEY_TAYLOR;
		//
		climate.mpet20[m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
			climate.mprec_pet20[m] = climate.hmprec_20[m].lastadd() / climate.mpet_year[m];
		else
			climate.mprec_pet20[m] = 0.0;

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

		climate.mtemp_20[19][m] = climate.hmtemp_20[m].lastadd();
		climate.mprec_20[19][m] = climate.hmprec_20[m].lastadd();

		climate.mpet_20[19][m] = climate.mpet_year[m];
		if(climate.mpet_year[m] > 0.0)
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
/** Called from crop_sowing_gridcell() last day of the year
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


//////////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Bondeau A, Smith PC, Zaehle S, Schaphoff S, Lucht W, Cramer W, Gerten D, Lotze-Campen H,
//   Müller C, Reichstein M & Smith B 2007. Modelling the role of agriculture for the 
//   20th century global terrestrial carbon balance. Global Change Biology, 13:679-706.
// Waha K, van Bussel LGJ, Müller C, and Bondeau A.2012. Climate-driven simulation of global 
//   crop sowing dates, Global Ecol Biogeogr 21:247-259