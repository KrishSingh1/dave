///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.cpp
//WK SIMFIRE doesn't actually do ignitions but assumes an ignition saturated regime
//WK please explain
//RLN Is that better (reference in the bottom of file)
/// \brief SIMFIRE burned area simulation by W. Knorr
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 16:02:51 +0100 (Tue, 24 Jan 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

//WK Is this the framework header file as remarked in (1)?
//RLN yes
#include "config.h"
#include "simfire.h"
#include "SimfireInput.h"

//WK state a purpose here: return SIMFIRE biome from time average vegetation
//WK characteristics, also refer to where these biomes are defined in the code
//WK state meaning of 'return 0'; here it seems to mean 'unvegetated',
//WK but in the IGBP-SIMFIRE mapping routine 'biome=0' means 'cropland/natural/urban'
//WK This latter biome could also be determined here using land use information
//WK in case the land use version of LPJ-GUESS is run
//WK in general I think this needs explaining how to handle SIMFIRE biomes
//WK using different versions of LPJ-GUESS (potential natural vegetation, land use,
//WK running with observed climate or future scenarios)
//RLN Please see comments in routine-headers below.
#define NFIREBIOMES 9

/// Get simfire data for a gridcell
void getsimfiredata(Gridcell& gridcell) {
		
	/* Called by: framework (framework.cpp)
	   Calls    : simfire_biome_mapping (local)
	   Reads SIMFIRE relevant info from SimfireInput.bin:
	   Hyde 3.1 population density
	   Monthly fire climatology
	*/

	Climate& climate = gridcell.climate;
	
	/// Paths to SIMFIRE binaries
	xtring file_simfire = param["file_simfire"].str;
	
	// open file, fill podp, monthly_burned_area and igbp_class for a gridcell
	// Fill static arrays/variables here
	SimfireInputArchive ark;
	
	if (!ark.open(file_simfire)) {
		fail("Could not open %s for input \n", (char*)file_simfire);
	}
	
	SimfireInput rec;
	rec.lon = gridcell.get_lon();
	rec.lat = gridcell.get_lat();
	
	if (!ark.getindex(rec)) {
		ark.close();
		fail("Grid cell not found in %s \n", (char*)file_simfire);
	}
	
	// Found the record, get the values
	
//CRM	// IGBP Land-Cover-Classification
//CRM	gridcell.igbp_class = (int)rec.igbp_class[0];
	
	// convert IGBP into simfire internal biomes
	simfire_biome_mapping(gridcell);
	
	// Monthly fire risk (W.Knorr)
	for (int m=0; m<12; m++) {
		climate.monthly_fire_risk[m] = rec.monthly_ba[m];
	}
	// Population density from HYDE 3.1
	for (int t=0; t<57; t++) {
		gridcell.hyde31_pop_density[t] = rec.pop_density[t];
	}		
	
	ark.close();
}

int update_fire_biome (Patch& patch, double lat) {

	/* Called by: simfire_biome_mapping (local)
	   Calls    : -
	   Computes current SIMFIRE biome for this 
	   gridcell depending on the last <n_year_biomeavg> years of
	   vegetation. 
	*/

	double fgrass=0.0; // grass fraction of all vegetation
	double fndlt=0.0;  // fraction of needle-leaf tress
	double fbrlt=0.0;  // fraction of broad-leaf trees
	double fshrb=0.0;  // fraction of woody vegetation that is shrubs
	double ftot=0.0;   // total FPAR of all individuals
	int biome=0;       // biome number
	int count[NFIREBIOMES]; // incidence count of biome in previous years;
	int count_max=0;   // maximum of 'count'

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation=patch.vegetation;

//CRM	// initialise fapar averaging array
//CRM	if ( date.year == 0 && date.day == 0 && ! restart ) {
//CRM		for (int i = 0; i<n_year_biomeavg; i++) {
//CRM			patch.avg_ftot[i]   = 0. ;
//CRM			patch.avg_fgrass[i] = 0. ;
//CRM			patch.avg_fndlt[i]  = 0. ;
//CRM			patch.avg_fbrlt[i]  = 0. ;
//CRM			patch.avg_fshrb[i]  = 0. ;
//CRM		}
//CRM	}

	// Loop through individuals of this patch
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
	
		if (indiv.id!=-1 && indiv.alive) { 
	
			if (indiv.pft.lifeform==GRASS) {
				fgrass += indiv.fpar_leafon;
			}
			else { // tree or shrub
				if (indiv.height>=2.0) { // tree
					if (indiv.pft.leafphysiognomy==NEEDLELEAF) {
						fndlt += indiv.fpar_leafon;
					}
					else { // broadleaf tree
						fbrlt += indiv.fpar_leafon;
					}
				}
				else { // shrub
					fshrb += indiv.fpar_leafon;
				}
			}
			ftot += indiv.fpar_leafon;
		}
		vegetation.nextobj(); // ... on to next individual
	}
	
	if ( ftot < 0.00000001 ) {
//WK maybe it should be -1, consistent with simfire_biome_mapping?
//RLN indeed
		return -1;
	}
	// re-normalize
	fgrass /= ftot;
	fndlt  /= ftot;
	fbrlt  /= ftot;
	fshrb  /= (1.00001-fgrass);

	// save current 
	int idx = date.year % n_year_biomeavg;  
	patch.avg_ftot  [idx] = ftot   ;
	patch.avg_fgrass[idx] = fgrass ;
	patch.avg_fndlt [idx] = fndlt  ;
	patch.avg_fbrlt [idx] = fbrlt  ;
	patch.avg_fshrb [idx] = fshrb  ;
//CRM	if ( date.year == 0 ) {
//CRM		for (int i = 0; i<n_year_biomeavg; i++) {
//CRM			patch.avg_ftot  [i] = ftot   ;
//CRM			patch.avg_fgrass[i] = fgrass ;
//CRM			patch.avg_fndlt [i] = fndlt  ;
//CRM			patch.avg_fbrlt [i] = fbrlt  ;
//CRM			patch.avg_fshrb [i] = fshrb  ;
//CRM		}
//CRM	}
	
	// generate running avereage 
	ftot   = 0.;
	fgrass = 0.;
	fndlt  = 0.;
	fbrlt  = 0.;
	fshrb  = 0.;
	for (int i = 0; i<n_year_biomeavg; i++) {
		ftot   += patch.avg_ftot  [i];
		fgrass += patch.avg_fgrass[i];
		fndlt  += patch.avg_fndlt [i];
		fbrlt  += patch.avg_fbrlt [i];
		fshrb  += patch.avg_fshrb [i];
	}
	ftot   /=  (double)n_year_biomeavg;
	fgrass /=  (double)n_year_biomeavg;
	fndlt  /=  (double)n_year_biomeavg;
	fbrlt  /=  (double)n_year_biomeavg;
	fshrb  /=  (double)n_year_biomeavg;

//WK: remove 'and' in next line; also note above about land use
//RLN I think, I removed what might be confusing, no?
	if (ftot<0.1 && fabs(lat)<50.0) {
		biome=8; } // barren or sparsely vegetated
	else if (ftot<0.1   && fabs(lat)>=50.0) {
		biome=7; } // tundra
	else if (patch.stand.landcover==CROPLAND) {
		biome=1; } // cropland
	else if (fshrb>=0.8 && fabs(lat)<50.0) {
		biome=5; } // shrubland
	else if (fshrb>=0.8 && fabs(lat)>=50.0) {
		biome=7; } // tundra
	else if (fgrass>=0.4) {
		biome=6; } // savanna or grassland
	else if (fndlt>=0.6) {
		biome=2; } // needle-leaf forest
	else if (fbrlt>=0.6) {
		biome=3; } // broad-leaf forest
	else {
		biome=4;   // mixed forest
		} 

	return biome;
}

//WK state a purpose and explain why is this mapping needed given that the
//WK routine above uses FPAR of different vegetation types to determine
//WK the SIMFIRE biome
//WK Find it confusing that there seem to be two ways of generating biome
//WK information - do they exist side by side, is there a switch in the ins
//WK file to choose between them etc., etc. - please explain
//RLN we will only use the biome-updating from LPJ-GUESS parameters. The rest was legacy code.
void simfire_biome_mapping(Gridcell& gridcell) {

	/* Called by: simfire_accounting_gridcell (local)
	              getsimfiredata (local)
	   Calls    : update_fire_biome(local)
	   Computes current SIMFIRE biome for this 
	   gridcell depending on the last <n_year_biomeavg> years of
	   vegetation. 
	*/

//CRM	// IGBP:
//CRM	//  0 Water bodies
//CRM	//  1 Evergreen Needleleaf Forest % 1: >60% cover, height>2m
//CRM	//  2 Evergreen Broadleaf Forest % 2: >60% cover, height>2m
//CRM	//  3 Deciduous Needleleaf Forest % 3: >60% cover, height>2m
//CRM	//  4 Deciduous Broadleaf Forest % 4: >60% cover, height>2m
//CRM	//  5 Mixed Forest % 5: >60% cover, height>2m, no forest type>60% cover
//CRM	//  6 Closed Shrubland % 6: >60% woody cover, height<2m
//CRM	//  7 Open Shrubland % 7: 10-60% woody cover, height<2m
//CRM	//  8 Woody Savanna % 8: 30-60% tree cover, height>2m, herbaceous or other understory
//CRM	//  9 Savanna % 9: 10-30% tree cover, height>2m, herbaceous or other understory
//CRM	// 10 Grassland % 10: <10% tree and shrub cover
//CRM	// 11 Permanent Wetland % 11: mixture of water and herbaceous or woody vegetation
//CRM	// 12 Cropland % 12
//CRM	// 13 Urban and Built-Up % 13
//CRM	// 14 Cropland/Natural Vegetation Mosaic % 14 none of forest, shrubland, cropland, 
//CRM	//    grassland >60% cover
//CRM	// 15 Permanent Snow and Ice % 15
//CRM	// 16 Barren or Sparsely Vegetated % 16: <10% vegetation cover all year
//CRM	// 17 Unclassified / No data
//CRM	
//CRM	// IGBP2BIOME MAPPING:
//CRM	// vegetation formation (should be derived from inputs to function)
//CRM	// -1 No vegetation
//CRM	//  0 Cropland/Urban/Natural Vegetation Mosaic (IGBP 12-14)
//CRM	//  1 Needleleaf forest (IGBP 1,3): >60% cover, height>2m
//CRM	//  2 Broadleaf forest (IGBP 2,4): >60% cover, height>2m
//CRM	//  3 Mixed forest (IGBP 4): >60% cover, height>2m, none >60%
//CRM	//  4 Shrubland (IGBP 6,7 and latitude<50): >10% woody cover, height<2m
//CRM	//  5 Savanna or Grassland (IGBP 8-10): herbaceous component present, <60% tree cover
//CRM	//  6 Tundra (IGBP 6,7,16 and latitude>=50): height<2m
//CRM	//  7 Barren or Sparsely Vegetated (IGBP 16 and latitude<50): <10% vegetation cover
//CRM	//
//CRM
//CRM	// mapping IGBP -> simfire-biomes
//CRM	const double igbp2simfirebiome[2][18] = {
//CRM		{ -1, 1, 2, 1, 2, 3, 4, 4, 5, 5, 5,-1, 0, 0, 0,-1, 7, -1 },  //  ! LAT  < 50 
//CRM		{ -1, 1, 2, 1, 2, 3, 6, 6, 5, 5, 5,-1, 0, 0, 0,-1, 6, -1 }}; //  ! LAT >= 50

	Climate& climate = gridcell.climate;
//CRM
//CRM	// Regions with dedicated parameter optimisation for SIMFIRE
//CRM	// 0: global
//CRM	// 1: Europe
//CRM	// 2: AUS-NZ
//CRM

	// determine gridcell's simfire biome
	// At start of spinup use IGBP
	//WK where does the IGBP class come from?
	//RLN Removed. 
//CRM	double lat = gridcell.get_lat();
//CRM	if ( date.year == 0 && date.day == 0 ) {
//CRM		if (abs(lat) >= 50.) {
//CRM			climate.simfire_biome  = igbp2simfirebiome[1][gridcell.igbp_class];
//CRM		}
//CRM		else {
//CRM			climate.simfire_biome  = igbp2simfirebiome[0][gridcell.igbp_class];
//CRM		}
//CRM	}
//CRM	//WK explain what is meant by "biome-shift" and by "later"
//CRM	// later use biome-shift
//CRM	else {
//CRM		int nobjs=0;
	std::vector<int> biomes;
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			
			biomes.push_back(update_fire_biome(patch, gridcell.get_lat()));
			stand.nextobj();
		}
		//CRM			nobjs+=stand.nobj;
		++gc_itr;
	}
		
	int count[NFIREBIOMES];
	int biome;
	int count_max=0; // maximum of 'count'
	// find and save most common biome number
	for (biome=0;biome<NFIREBIOMES;biome++) count[biome]=0;
	for (int idx = 0; idx < biomes.size(); idx++) {
		count[biomes[idx]]++;
	}
	for (biome=0;biome<NFIREBIOMES;biome++) {
		count_max=max(count_max,count[biome]);
		//if(date.year>500) dprintf("biome%d %d ", biome,count[biome] );
	}

	for (biome=0;biome<NFIREBIOMES && count[biome]<count_max;biome++) {
	}
	//if(date.year>500) dprintf("count a1 %d tot %d biome %d \n",count[1],nobjs,biome);
//CRM	if ( ((double)count[1]/(double)nobjs)>0.33) {
//CRM		biome = 1;
//CRM		//if(date.year>500) dprintf("count b1 %d tot %d biome %d \n",count[1],nobjs,biome);
//CRM	}
	// BLAZE-biome is -1 of SIMFIRE-biome classifcation
	climate.simfire_biome  = biome - 1;
//CRM}
}

//==============================================================================
//WK state that before the earliest time point ('poptime') we use constant values
//WK and after the last a linear extrapolation using the last two entries
//WK Maybe also make clear that these are historical data, I have
//WK also gridded scenario fields every ten years until 2100
//WK (but this would require an option for the choice of scenario)
// INTERPOLATE HYDE 3.1 POPULATION DENSITY BETWEEN TIME-STEPS
//==============================================================================
//RLN We will not provide future scenarios. 

void simfire_update_pop_density(Gridcell& gridcell) {

	/* Called by: simfire_accounting_gridcell (local)
	   Calls    : -
	   Computes population density from the Hyde 3.1 dataset. 
	   Annual data is computed by linearly interpolating between the existing values.
	   Before 10000 BC the 10000 BC value is used, after 2005 linear extrapolation 
	   using the change between the last two values is performed 
	*/

	const int npopt = 57;
	// years at which pop data is available in HYDE3.1
	const int poptime[npopt]  = {-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,0,
				    100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
				    1500,1600,1700,1710,1720,1730,1740,1750,1760,1780,1790,1810,
				    1820,1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,
				    1940,1950,1960,1970,1980,1990,2000,2005};
	// get calendar-year
	int cyear = date.get_calendar_year();

	// find start and end year index of pop interpolation
	int idx = 0 ;
      	while (poptime[idx] < cyear) idx++;
	double popd;
	if ( cyear <= poptime[0] ) {
		// use first year's value (10000 BC) for earlier years.
		popd = gridcell.hyde31_pop_density[0];
	}
	else if ( cyear >= poptime[npopt-1] ) {
		// linearly extrapolate latest growth/decline
		popd = gridcell.hyde31_pop_density[npopt-1] + 
			(gridcell.hyde31_pop_density[npopt-1]-gridcell.hyde31_pop_density[npopt-2]) /
			(double)(poptime[npopt-1] - poptime[npopt-2]) * (double)(cyear-poptime[npopt-1]);
	}
	else {
		double interpf = (double)(cyear-poptime[idx-1]) /
			(double)( poptime[idx]-poptime[idx-1] );
		popd = (1. - interpf) * gridcell.hyde31_pop_density[idx-1] + 
			interpf * gridcell.hyde31_pop_density[idx];
	}

	gridcell.pop_density = max(0.,popd);
}
	

/// Called each day from dailyaccounting
//WK I think this has been implemented very well here!
//RLN Thanks :)
void simfire_accounting_gridcell(Gridcell& gridcell) {
	
	/* Called by: dailyaccounting_gridcell   (driver.cpp)
	   Calls    : simfire_biome_mapping      (local)
	              simfire_update_pop_density (local)
	   Updates SIMFIRE's Max Annual Mesterov Index
	   and running mean of max annual FPAR (from canexch.cpp)
	   Updates fire biome 
	*/

	Climate& climate = gridcell.climate;
	// absolute upper boundary for the accumulative nesterov index
	const double maximum_nesterov = 1000000; //150000.;

	// check whether this day is the first day of simulation 
	// (i.e. start of spinup or first day after restart
	bool is_first_day = ( date.day == 0 && ( date.year == 0 || 
		( restart && date.year == state_year ) ) );

	if (date.day == 0 ) {
		// Set global simfire region as fixed: Global=0
		gridcell.simfire_region = 0;

		// Determine SIMFIRE biome for this year
		simfire_biome_mapping(climate.gridcell);

		// update population density
		simfire_update_pop_density(climate.gridcell);
		

		// initialise averaging array 
		if ( date.year == 0 ) {
			for(int i=0;i<avg_interv_fapar;i++) { 
				climate.recent_max_fapar[i] = 0.5;
			}
			climate.ann_max_fapar = 0.5;	
		} 
		else {
//CRM			climate.recent_max_fapar[a] = climate.cur_max_fapar;
			double avg = 0.;
			for(int i=0;i<avg_interv_fapar;i++) { 
				avg += climate.recent_max_fapar[i];
			}
			climate.ann_max_fapar = avg / (double) avg_interv_fapar;
		}
		// finally (re)set this years max fapar
		climate.cur_max_fapar = 0.0;

		// set Max annual Nesterov Index on first day of simulation
		if ( is_first_day ) {
			for ( int i=0; i<12; i++) 
				climate.monthly_max_nesterov[i] = 0.;
			if ( restart )
				climate.monthly_max_nesterov[11] = climate.max_nesterov;
			else
				climate.monthly_max_nesterov[11] = 10000. * cos(gridcell.get_lat());
		}
	} 	
	// multi-year accounting of maximum annual fapar	
	else if ( date.islastday && date.islastmonth ) {
		int a = date.year % avg_interv_fapar;
		climate.recent_max_fapar[a] = climate.cur_max_fapar;
	}

        if ( date.dayofmonth == 0 ) {
		double mnest = 0.;
		for ( int i=0; i<12; i++) 
			if ( climate.monthly_max_nesterov[i] > mnest )
				mnest = climate.monthly_max_nesterov[i];
		climate.max_nesterov = mnest;
		climate.monthly_max_nesterov[date.month] = 0. ;
	}

        if (  climate.monthly_max_nesterov[date.month] < climate.cur_nesterov )
		climate.monthly_max_nesterov[date.month] = climate.cur_nesterov;

	// PATCHLOOP FOR fpar
	int cnt= 0;
	double run_fapar = 0.;
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			run_fapar += (1. - patch.fpar_ff);
			//initialise averaging array
			if ( date.year == 0 && date.day == 0 ) {
				for (int i = 0; i<n_year_biomeavg; i++) {
					patch.avg_ftot  [i] = 0. ;
					patch.avg_fgrass[i] = 0. ;
					patch.avg_fndlt [i] = 0. ;
					patch.avg_fbrlt [i] = 0. ;
					patch.avg_fshrb [i] = 0. ;
				}
			}
			cnt += 1;
			stand.nextobj();
		}		
		++gc_itr;
	}

	// average over each patch 
	run_fapar /= (double) cnt;
	// awkward solution but no idea....
	if ( date.year == 0 && date.day == 0 && run_fapar > 0.99999 ) {
		run_fapar = 0.;
	}

	// update the this years maximum
	climate.cur_max_fapar = max(run_fapar, climate.cur_max_fapar);

	// compute running Nesterov index
	if ( climate.prec >= 3. || climate.tmax - climate.tmin < 4. ) {
		climate.cur_nesterov = 0.0; 
	}
	else {
		climate.cur_nesterov += ( climate.tmax - climate.tmin + 4. ) * climate.tmax ;
	}
	climate.cur_nesterov = min(climate.cur_nesterov,maximum_nesterov) ;

	// finally update Max Annual Mesterov Index
	if (climate.cur_nesterov > climate.max_nesterov ) 
		climate.max_nesterov = climate.cur_nesterov ;
}

double simfire_ba(Climate& climate, Gridcell& gridcell) {

	/* Called by:  blaze_burned_area(blaze.cpp)
	   Calls    :  -
	   Calculate burned area in ha following Knorr 2014. 
	*/
//WK See comments on ignitions vs. burned area above
//WK Make comment here on how the SIMFIRE region is set
//WK Maybe also explain here that SIMFIRE only predicts
//WK the annual mean burned area, and that the seasonal cycle
//WK at monthly time steps (?) is set from observations
//WK obtained from ?

//CLN I thought I got this from you. Will follow uop on this. 
//RLN I will remove the two non-global cases but at this stage I can't make changes to the branch. 
	// Regions with dedicated parameter optimisation for SIMFIRE
	// 0: global
	// 1: Europe
	// 2: AUS-NZ

	const double a[3][8] = {
		{ 0.110,  0.095    ,0.092  ,0.127  ,0.470  ,0.889 ,0.059  ,0.113  }, //	GLOBAL
		{ 0.02589,0.0008087,0.04896,0.06248,0.01966,0.1191,0.01872,0.08873}, //	EUR   
		{ 0.06974,0.6535   ,0.6341 ,0.6438 ,2.209  ,1.710 ,    0  ,2.572  }};//	Australia-NZ
	
	const double b[3] = {
		0.905,   // GLOBAL
		0.9164,  // EUR
		1.297 }; // Australia-NZ   
	
	const double c[3] = {
		0.860,   // GLOBAL
		0.4876,  // EUR
		1.038 }; // Australia-NZ   
	
	const double e[3] = {
		-0.0168, // GLOBAL
		-0.017,  // EUR
		-0.2131 }; // Australia-NZ original *corrected LN
		//	 -0.05 }; // Australia-NZ   

	const double scalar = 1.0e-5;

	int ri = gridcell.simfire_region;

	// return if improper biome-type
	if (climate.simfire_biome == -1) return 0.;

	// fPAR correction Knorr
//WK Needs some explanation: the FPAR correction is used
//WK when GUESS generates the FPAR, but not when observed
//WK FPAR is used to compute burned area
//RLN This should have been used all along...
	const double fpar_corr1 = 0.428;
	const double fpar_corr2 = 0.148;

	double fpar_cor = fpar_corr1 * climate.ann_max_fapar + fpar_corr2 * climate.ann_max_fapar * 
	  climate.ann_max_fapar;

	// compute annual burned area
	double ba = a[ri][climate.simfire_biome] * 
		pow(fpar_cor, b[ri]) *
		pow((scalar * climate.max_nesterov), c[ri]) *
		exp(e[ri] * gridcell.pop_density);

	//CLNdprintf("sim ba %f fapar %f  popd %f biome %d \n",ba ,fpar_cor, gridcell.pop_density,climate.simfire_biome);
	

//CRM	// compute annual burned area
//CRM	double ba = a[ri][climate.simfire_biome] * 
//CRM		pow(climate.ann_max_fapar, b[ri]) *
//CRM		pow((scalar * climate.max_nesterov), c[ri]) *
//CRM		exp(e[ri] * gridcell.pop_density);
//CRM	

	// compute daily burnt_area
	ba *= climate.monthly_fire_risk[date.month] /
		(double)date.ndaymonth[date.month];

//CRM	// disaggregate annual burnt area into sub-annual timescales
//CRM	if (blaze_tstep == DAILY) {
//CRM		ba *= climate.monthly_fire_risk[date.month] /
//CRM			date.ndaymonth[date.month];
//CRM	} 
//CRM	else if (blaze_tstep == MONTHLY) {
//CRM		ba *= climate.monthly_fire_risk[date.month];
//CRM	} 
//CRM	else if (blaze_tstep == SEASONAL) {
//CRM		int s = (double)date.month / 3.0;
//CRM		ba *=   climate.monthly_fire_risk[s*3  ] +
//CRM			climate.monthly_fire_risk[s*3+1] +
//CRM			climate.monthly_fire_risk[s*3+2];
//CRM	}
	
	// keep track of area burnt so far this year
	climate.acc_areaburnt += ba;

	return ba;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Knorr, W. et al., Impact of human population density on fire frequency at the 
//  global scale, BIOGEOSCIENCES, 11, 4, 2014, DOI: 10.5194/bg-11-1085-2014
