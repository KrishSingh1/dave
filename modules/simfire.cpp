///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.cpp
/// \brief SIMFIRE ignition simulation by W. Knorr
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

#include "config.h"
#include "simfire.h"

/// SIMFIRE biome mapping
#define NFIREBIOMES 9
int update_fire_biome (Patch& patch, double lat) {

	double fgrass=0.0; // grass fraction of all vegetation
	double fndlt=0.0; // fraction of needle-leaf tress
	double fbrlt=0.0; // fraction of broad-leaf trees
	double fshrb=0.0; // fraction of woody vegetation that is shrubs
	double ftot=0.0; // total FPAR of all individuals
	int biome=0; // biome number
	int count[NFIREBIOMES]; // incidence count of biome in previous years;
	int count_max=0; // maximum of 'count'

	// Obtain reference to Vegetation object for this patch
	Vegetation& vegetation=patch.vegetation;

	// initialise fapar averaging array
	if ( date.year == 0 && date.day == 0 && ! restart ) {
		for (int i = 0; i<n_year_biomeavg; i++) {
			patch.avg_ftot[i]   = 0. ;
			patch.avg_fgrass[i] = 0. ;
			patch.avg_fndlt[i]  = 0. ;
			patch.avg_fbrlt[i]  = 0. ;
			patch.avg_fshrb[i]  = 0. ;
		}
	}

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

	if ( ftot < 0.00000001 )   
		return 0;
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
	ftot   /=  n_year_biomeavg;
	fgrass /=  n_year_biomeavg;
	fndlt  /=  n_year_biomeavg;
	fbrlt  /=  n_year_biomeavg;
	fshrb  /=  n_year_biomeavg;

	// assign biome (neglecting and agricultural land use)
	if (ftot<0.1 && fabs(lat)<50.0) {
		biome=8; } // barren or sparsely vegetated
	else if (ftot<0.1) {
		biome=7; } // tundra
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
		biome=4;  // mixed forest
	}
	return biome;
}

void simfire_biome_mapping(Gridcell& gridcell) {
	
	// IGBP:
	//  0 Water bodies
	//  1 Evergreen Needleleaf Forest % 1: >60% cover, height>2m
	//  2 Evergreen Broadleaf Forest % 2: >60% cover, height>2m
	//  3 Deciduous Needleleaf Forest % 3: >60% cover, height>2m
	//  4 Deciduous Broadleaf Forest % 4: >60% cover, height>2m
	//  5 Mixed Forest % 5: >60% cover, height>2m, no forest type>60% cover
	//  6 Closed Shrubland % 6: >60% woody cover, height<2m
	//  7 Open Shrubland % 7: 10-60% woody cover, height<2m
	//  8 Woody Savanna % 8: 30-60% tree cover, height>2m, herbaceous or other understory
	//  9 Savanna % 9: 10-30% tree cover, height>2m, herbaceous or other understory
	// 10 Grassland % 10: <10% tree and shrub cover
	// 11 Permanent Wetland % 11: mixture of water and herbaceous or woody vegetation
	// 12 Cropland % 12
	// 13 Urban and Built-Up % 13
	// 14 Cropland/Natural Vegetation Mosaic % 14 none of forest, shrubland, cropland, 
	//    grassland >60% cover
	// 15 Permanent Snow and Ice % 15
	// 16 Barren or Sparsely Vegetated % 16: <10% vegetation cover all year
	// 17 Unclassified / No data
	
	// IGBP2BIOME MAPPING:
	// vegetation formation (should be derived from inputs to function)
	// -1 No vegetation
	//  0 Cropland/Urban/Natural Vegetation Mosaic (IGBP 12-14)
	//  1 Needleleaf forest (IGBP 1,3): >60% cover, height>2m
	//  2 Broadleaf forest (IGBP 2,4): >60% cover, height>2m
	//  3 Mixed forest (IGBP 4): >60% cover, height>2m, none >60%
	//  4 Shrubland (IGBP 6,7 and latitude<50): >10% woody cover, height<2m
	//  5 Savanna or Grassland (IGBP 8-10): herbaceous component present, <60% tree cover
	//  6 Tundra (IGBP 6,7,16 and latitude>=50): height<2m
	//  7 Barren or Sparsely Vegetated (IGBP 16 and latitude<50): <10% vegetation cover
	//

	const double igbp2simfirebiome[2][18] = {
		{ -1, 1, 2, 1, 2, 3, 4, 4, 5, 5, 5,-1, 0, 0, 0,-1, 7, -1 },  //  ! LAT  < 50 
		{ -1, 1, 2, 1, 2, 3, 6, 6, 5, 5, 5,-1, 0, 0, 0,-1, 6, -1 }}; //  ! LAT >= 50

	//CLN here auswahlkriterien!
	Climate& climate = gridcell.climate;

	// Regions with dedicated parameter optimisation for SIMFIRE
	// 0: global
	// 1: Europe
	// 2: AUS-NZ
 	// come up with a EUrope, ANZ definition.
	gridcell.simfire_region = 0; // Global

	double lat = gridcell.get_lat();

	if ( date.year == 0 && date.day == 0 ) {
		if (abs(lat) >= 50.) {
			climate.simfire_biome  = igbp2simfirebiome[1][gridcell.igbp_class];
		}
		else {
			climate.simfire_biome  = igbp2simfirebiome[0][gridcell.igbp_class];
		}
	}
	else {
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
		for (biome=0;biome<NFIREBIOMES;biome++) count_max=max(count_max,count[biome]);
		for (biome=0;biome<NFIREBIOMES && count[biome]<count_max;biome++) {
		}
		// Lars biome is -1 of Wolfgangs old biome classifcation
		climate.simfire_biome  = biome - 1;
	}	

}

//==============================================================================
// INTERPOLATE HYDE 3.1 POPULATION DENSITY BETWEEN TIME-STEPS
//==============================================================================
void simfire_update_pop_density(Gridcell& gridcell) {

	const int npopt = 57;
	const int poptime[npopt]  = {-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,0,
				    100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
				    1500,1600,1700,1710,1720,1730,1740,1750,1760,1780,1790,1810,
				    1820,1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,
				    1940,1950,1960,1970,1980,1990,2000,2005};

	int cyear = date.get_calendar_year();

	// start and end year index of pop interpolation
	int idx = 0 ;
      	while (poptime[idx] < cyear) idx++;
	double popd;
	if ( cyear <= poptime[0] ) {
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
			(double)(poptime[idx-1] - poptime[idx]);
		popd = (1. - interpf) * gridcell.hyde31_pop_density[idx-1] + 
			interpf * gridcell.hyde31_pop_density[idx];
	}

	gridcell.pop_density = max(0.,popd);
}
	

/// Called each day from dailyaccounting
void simfire_accounting_gridcell(Gridcell& gridcell) {
	
	// DESCRIPTION
	// Updates SIMFIRE's Max Annual Mesterov Index
	// and running mean of max annual FPAR (from canexch.cpp)
	// as well as checks for biome shifting
	// L. Nieradzik 03/2015
	Climate& climate = gridcell.climate;
	// absolute upper boundary for the accumulative nesterov index
	const double maximum_nesterov = 150000.;

	if (date.day == 0 ) {
		// set SIMFIRE biomes based on IGBP classification
		// NOW DONE IN getgridcell
		simfire_biome_mapping(climate.gridcell);

		// update population density
		simfire_update_pop_density(climate.gridcell);
		
		// reset Max annual Nesterov Index
		// CLN is tjhat correct????? Reset regarding to burntime!!!!!
		climate.max_nesterov = 0.0;
		
		// initialise averaging array (CLN MOVE TO restartvalues!)
		if ( date.year == 0 && ! restart ) {
			for(int i=0;i<avg_interv_fapar;i++) { 
				climate.recent_max_fapar[i] = 0.5;
			}
			climate.ann_max_fapar = 0.5;	
		} else {
			int a = date.year % avg_interv_fapar;
			climate.recent_max_fapar[a] = climate.cur_max_fapar;
			double avg = 0.;
			for(int i=0;i<avg_interv_fapar;i++) { 
				avg += climate.recent_max_fapar[i];
			}
			climate.ann_max_fapar = avg / (double) avg_interv_fapar;
		}
		// finally (re)set this years max fapar
		climate.cur_max_fapar = 0.0;
	}

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
	// finally update Max Annual Mesterov Index
	if (climate.cur_nesterov > climate.max_nesterov ) 
		climate.max_nesterov = min(climate.cur_nesterov,maximum_nesterov) ;
}

double simfire_ba(Climate& climate, Gridcell& gridcell) {

	/* Called by:  blaze_ignition(blaze.cpp)
	   Calls    :  -
	   calculate annual (for now) burned area in ha
	*/

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

	//	Gridcell gridcell = climate.gridcell;
	
	int ri = gridcell.simfire_region;

	/*
	  Lines below need to go into init part of LPJ-GUESS
	//reading in the data:
	// Retrieve grid information from burned area annual cycle data base
	init_sql_acba(param["file_sql_ac_burn"].str);

	// Get annual cycle of burned area and set month marking start of fire season
	read_sql_acba(lon,lat,param["file_sql_ac_burn"].str);
	*/

	//	dprintf("CLN simf ba reg: %d \n"  ,gridcell.simfire_region);
	
	/*
	dprintf("CLN simf ba biome: %d \n",climate.simfire_biome);
	dprintf("CLN simf a: %f \n",a[ri][climate.simfire_biome]);
	dprintf("CLN simf ann_fapar: %f \n",climate.ann_max_fapar);
	dprintf("CLN simf b[ri]: %f \n",b[ri]);
	dprintf("CLN simf climate.max_nesterov: %f \n",climate.max_nesterov);
	dprintf("CLN simf c[ri]: %f \n",c[ri]);
	dprintf("CLN simf e[ri]: %f \n",e[ri]);
	dprintf("CLN simf gridcell.pop_density: %f \n",gridcell.pop_density);
	*/

	// return if improper biome-type
	if (climate.simfire_biome == -1) return 0.;

	// fPAR correction Knorr
	/*const double fpar_corr1 = 0.428;
	const double fpar_corr2 = 0.148;

	double fpar_cor = fpar_corr1 * climate.ann_max_fapar + fpar_corr2 * climate.ann_max_fapar * 
	  climate.ann_max_fapar;

	double ba = a[ri][climate.simfire_biome] * 
		pow(fpar_cor, b[ri]) *
		pow((scalar * climate.max_nesterov), c[ri]) *
		exp(e[ri] * gridcell.pop_density);
	*/

	//  annual burned area
	double ba = a[ri][climate.simfire_biome] * 
		pow(climate.ann_max_fapar, b[ri]) *
		pow((scalar * climate.max_nesterov), c[ri]) *
		exp(e[ri] * gridcell.pop_density);

	if (blaze_tstep == DAILY) {
		ba *= climate.monthly_fire_risk[date.month] /
			date.ndaymonth[date.month];
	} 
	else if (blaze_tstep == MONTHLY) {
		ba *= climate.monthly_fire_risk[date.month];
	} 
	else if (blaze_tstep == SEASONAL) {
		// CLN CHECK FOR SEASON SETTING !!!
		int s = (double)date.month / 3.0;
		ba *=   climate.monthly_fire_risk[s*3  ] +
			climate.monthly_fire_risk[s*3+1] +
			climate.monthly_fire_risk[s*3+2];
	}

	climate.acc_areaburnt += ba;
	
	/*dprintf("CLN ba1 : %f \n",a[ri][climate.simfire_biome] * pow(climate.ann_max_fapar, b[ri]));
	dprintf("CLN ba2 : %f \n",pow((scalar * climate.max_nesterov), c[ri]) );
	dprintf("CLN ba3 : %f \n",exp(e[ri] * gridcell.pop_density));
	dprintf("CLN ba : %f \n",ba);
	*/
	// this is for debugging
	/*	char command[ 1024 ];
	long pid = ( long ) ::getpid(); // use long to ensure correct format specifier
	sprintf( command, "pmap %ld", pid );
	system( command );	*/
	// this is for debugging until here

	return ba;
}
