///////////////////////////////////////////////////////////////////////////////////////
/// \file simfire.cpp
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

#include "config.h"
#include "simfire.h"
#include "simfire_input.h"

#define NFIREBIOMES 9

int update_fire_biome(Patch& patch, double lat) {

	/* Called by: simfire_biome_mapping (local)
	   Calls    : -
	   Computes current SIMFIRE biome for this 
	   gridcell depending on the last <n_year_biomeavg> years of
	   vegetation. 
	   SIMFIRE BIOMES:
	   0 no veg/no data
	   1 Cropland/Urban/Natural Vegetation Mosaic (IGBP 12-14)
	   2 Needleleaf forest (IGBP 1,3): >60% cover, height>2m
	   3 Broadleaf forest (IGBP 2,4): >60% cover, height>2m
	   4 Mixed forest (IGBP 4): >60% cover, height>2m, none >60%
	   5 Shrubland (IGBP 6,7 and latitude<50): >10% woody cover, height<2m
	   6 Savanna or Grassland (IGBP 8-10): herbaceous component present, <60% tree cover
	   7 Tundra (IGBP 6,7,16 and latitude>=50): height<2m
	   8 Barren or Sparsely Vegetated (IGBP 16 and latitude<50): <10% vegetation cover
	*/

	double fgrass=0.0; // grass fraction of all vegetation
	double fndlt=0.0;  // fraction of needle-leaf tress
	double fbrlt=0.0;  // fraction of broad-leaf trees
	double fshrb=0.0;  // fraction of woody vegetation that is shrubs
	double ftot=0.0;   // total FPAR of all individuals
	int biome=0;       // biome number

    // Obtain reference to Vegetation object for this patch
	Vegetation& vegetation=patch.vegetation;

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
		return 0; // no data
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

	if (ftot<0.1 && fabs(lat)<50.0) {
		biome=SF_BARREN; } // barren or sparsely vegetated
	else if (ftot<0.1   && fabs(lat)>=50.0) {
		biome=SF_TUNDRA; } // tundra
	else if (patch.stand.landcover==CROPLAND) {
		biome=SF_CROP; } // cropland
	else if (fshrb>=0.8 && fabs(lat)<50.0) {
		biome=SF_SHRUBS; } // shrubland
	else if (fshrb>=0.8 && fabs(lat)>=50.0) {
		biome=SF_TUNDRA; } // tundra
	else if (fgrass>=0.4) {
		biome=SF_SAVANNA; } // savanna or grassland
	else if (fndlt>=0.6) {
		biome=SF_NEEDLELEAF; } // needle-leaf forest
	else if (fbrlt>=0.6) {
		biome=SF_BROADLEAF; } // broad-leaf forest
	else {
		biome=SF_MIXED_FOREST;   // mixed forest
	} 

	return biome;
}

void simfire_biome_mapping(Gridcell& gridcell) {

	/* Called by: simfire_accounting_gridcell (local)
		      getsimfiredata (local)
	   Calls    : update_fire_biome(local)
	   Computes current SIMFIRE biome for this 
	   gridcell depending on the last <n_year_biomeavg> years of
	   vegetation. 
	*/

	Climate& climate = gridcell.climate;

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
	for (biome=0;biome<NFIREBIOMES;biome++) {
		count_max=max(count_max, count[biome]);
	}

	for (biome = 0; biome<NFIREBIOMES && count[biome] < count_max; biome++) {
	}
	climate.simfire_biome = biome ;
}

// Get simfire data for a gridcell
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
	
	// convert IGBP into simfire internal biomes
	simfire_biome_mapping(gridcell);
	
	// Monthly fire risk (W.Knorr)
	for (int m=0; m<12; m++) {
		climate.monthly_fire_risk[m] = rec.monthly_burned_area[m];
	}
	// Population density from HYDE 3.1
	for (int t=0; t<57; t++) {
		gridcell.hyde31_pop_density[t] = rec.pop_density[t];
	}		
	
	ark.close();
}

void simfire_update_pop_density(Gridcell& gridcell) {

	/* Called by: simfire_accounting_gridcell (local)
	   Calls    : -
	   Computes population density from the Hyde 3.1 dataset. 
	   Annual data is computed by linearly interpolating between the existing values.
	   Before 10000 BC the 10000 BC value is used, after 2005 linear extrapolation 
	   using the change between the last two values is performed 
	*/

    // number of entries for Population data
	const int NPOPT = 57;
	// years at which population-data is available in HYDE3.1
	const int POPTIME[NPOPT]  = {-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,0,
				    100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
				    1500,1600,1700,1710,1720,1730,1740,1750,1760,1780,1790,1810,
				    1820,1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,
				    1940,1950,1960,1970,1980,1990,2000,2005};
	// get calendar-year
	int cyear = date.get_calendar_year();

	// find start and end year index of pop interpolation
	int idx = 0 ;
      	while (POPTIME[idx] < cyear) idx++;
	double popd;
	if ( cyear <= POPTIME[0] ) {
		// use first year's value (10000 BC) for earlier years.
		popd = gridcell.hyde31_pop_density[0];
	}
	else if ( cyear >= POPTIME[NPOPT-1] ) {
		// linearly extrapolate latest growth/decline
		popd = gridcell.hyde31_pop_density[NPOPT-1] + 
			(gridcell.hyde31_pop_density[NPOPT-1]-gridcell.hyde31_pop_density[NPOPT-2]) /
			(double)(POPTIME[NPOPT-1] - POPTIME[NPOPT-2]) * (double)(cyear-POPTIME[NPOPT-1]);
	}
	else {
        // interpolate between two entries
		double interpf = (double)(cyear-POPTIME[idx-1]) /
			(double)( POPTIME[idx]-POPTIME[idx-1] );
		popd = (1. - interpf) * gridcell.hyde31_pop_density[idx-1] + 
			interpf * gridcell.hyde31_pop_density[idx];
	}

	gridcell.pop_density = max(0.,popd);
}
	

/// Called each day from dailyaccounting
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
	const double MAXIMUM_NESTEROV = 1000000; //150000.;

	if (date.day == 0 ) {
		// Set global simfire region as fixed: Global=0
		gridcell.simfire_region = 0;

		// Determine SIMFIRE biome for this year
		simfire_biome_mapping(gridcell);

		// update population density
		simfire_update_pop_density(gridcell);
		

		// initialise averaging array 
		if ( date.year == 0 ) {
			for(int i=0;i<avg_interv_fapar;i++) { 
				climate.recent_max_fapar[i] = 0.5;
			}
			climate.ann_max_fapar = 0.5;	

			// initialize Max annual Nesterov Index on first day of simulation
			for ( int i=0; i<12; i++) {
				climate.monthly_max_nesterov[i] = 0.;
			}
			climate.cur_nesterov = 0.;
		} 
		else {
			double avg = 0.;
			for(int i=0;i<avg_interv_fapar;i++) { 
				avg += climate.recent_max_fapar[i];
			}
			climate.ann_max_fapar = avg / (double) avg_interv_fapar;
		}
		// finally (re)set this years max fapar
		climate.cur_max_fapar = 0.0;

	} 	
	// multi-year accounting of maximum annual fapar	
	else if ( date.islastday && date.islastmonth ) {
		int a = date.year % avg_interv_fapar;
		climate.recent_max_fapar[a] = climate.cur_max_fapar;
	}

	// update running Maximum Nesterov index array at beginning of month  
	if ( date.dayofmonth == 0 ) {
		double mnest = 0.;
		for ( int i=0; i<12; i++) 
			if ( climate.monthly_max_nesterov[i] > mnest )
				mnest = climate.monthly_max_nesterov[i];
		climate.max_nesterov = mnest;
		climate.monthly_max_nesterov[date.month] = 0. ;
	}

	// update current month's Maximum Nesterov index
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
			if ( ! ( date.year == 0 && date.day == 0 ) ) {
				run_fapar += (1. - patch.fpar_ff);
			}
			//initialise averaging array
			if (date.year == 0 && date.day == 0) {
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

	// update the this years maximum
	climate.cur_max_fapar = fmax(run_fapar, climate.cur_max_fapar);
        double bleu = pow(run_fapar, 2.);

	// compute running Nesterov index
    if ( climate.prec >= 3. || climate.tmax - climate.tmin < 4. ) {
		climate.cur_nesterov = 0.0; 
	}
	else {
		climate.cur_nesterov += ( climate.tmax - climate.tmin + 4. ) * climate.tmax ;
	}
	climate.cur_nesterov = min(climate.cur_nesterov,MAXIMUM_NESTEROV) ;

	// finally update Max Annual Mesterov Index
	if (climate.cur_nesterov > climate.max_nesterov ) 
		climate.max_nesterov = climate.cur_nesterov ;
}

/// Calculate burned area in ha following Knorr 2014.
double simfire_burned_area(Climate& climate) {

	// globally trained parameters 
	const double A[8] = { 0.110,  0.095    ,0.092  ,0.127  ,0.470  ,0.889 ,0.059  ,0.113  }; 
	const double B = 0.905;  
	const double C = 0.860; 
	const double E = -0.0168; 
	const double SCALAR = 1.0e-5;

	// return if improper biome-type
	if (climate.simfire_biome == 0) return 0.;

	// fPAR correction Knorr
	const double FPAR_CORR1 = 0.428;
	const double FPAR_CORR2 = 0.148;
	double fpar_cor = FPAR_CORR1 * climate.ann_max_fapar + FPAR_CORR2 * climate.ann_max_fapar * 
	  climate.ann_max_fapar;

	Gridcell& gridcell = climate.gridcell;

	// compute annual burned area
	double burned_area = A[climate.simfire_biome-1] * 
		pow(fpar_cor, B) *
		pow((SCALAR * climate.max_nesterov), C) *
		exp(E * gridcell.pop_density);

    // compute daily burnt_area
	burned_area *= climate.monthly_fire_risk[date.month] /
		(double)date.ndaymonth[date.month];

	// keep track of area burnt so far this year
	climate.acc_areaburnt += burned_area;

	return burned_area;
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Knorr, W. et al., Impact of human population density on fire frequency at the 
//  global scale, BIOGEOSCIENCES, 11, 4, 2014, DOI: 10.5194/bg-11-1085-2014
