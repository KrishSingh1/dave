///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK SOURCE CODE FILE
//
// Framework:             LPJ-GUESS Combined Modular Framework
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
// Header file name:      guess.h
// Source code file name: guess.cpp
// Written by:            Ben Smith
// Version dated:         2002-12-16
// Updated:               2010-11-22


#include "config.h"
#include "guess.h"

#include "guessio.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are declared in the framework header file, and defined here.
// They are accessible throughout the model code.

Date date; // object describing timing stage of simulation
vegmodetype vegmode; // vegetation mode (population, cohort or individual)
int npatch; // number of patches in each stand (should always be 1 in population mode); cropland stands always have 1 patch
double patcharea; // patch area (m2) (individual and cohort mode only)
bool ifdailynpp; // whether NPP calculations performed daily (alt: monthly)
bool ifdailydecomp;
	// whether soil decomposition calculations performed daily (alt: monthly)
bool ifbgestab; // whether background establishment enabled (individual, cohort mode)
bool ifsme;
	// whether spatial mass effect enabled for establishment (individual, cohort mode)
bool ifstochestab; // whether establishment stochastic (individual, cohort mode)
bool ifstochmort; // whether mortality stochastic (individual, cohort mode)
bool iffire; // whether fire enabled
bool ifdisturb;
	// whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
bool ifcalcsla; // whether SLA calculated from leaf longevity (alt: prescribed)
int estinterval; // establishment interval in cohort mode (years)
double distinterval;
	// generic patch-destroying disturbance interval (individual, cohort mode)
int npft; // number of possible PFTs
bool iffast;
bool ifcdebt;

// guess2008 - new inputs from the .ins file
bool ifsmoothgreffmort;				// smooth growth efficiency mortality
bool ifdroughtlimitedestab;			// whether establishment affected by growing season drought
bool ifrainonwetdaysonly;			// rain on wet days only (1, true), or a little every day (0, false); 
bool ifspeciesspecificwateruptake;	// water uptake is species specific 

bool run_landcover;
bool run[NLANDCOVERTYPES];
bool lcfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;				// If a slow harvested product pool is included in patchpft.
int nyear_spinup;		

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell,Pftlist& pftlist) {
	landcovertype landcover;

	getlandcover(gridcell,pftlist);		//Gets gridcell.landcoverfrac from landcover input file(s) or ins-file.

	for(int i=0;i<NLANDCOVERTYPES;i++) { //For all landcover types without subclasses
//		if(i!=CROPLAND) {					// cropland subclasses turned off in this version
			if(gridcell.landcoverfrac[i]>0.0) {
				if(run[i]) {
					landcover=(landcovertype)i;
					Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);

					pftlist.firstobj();
					while (pftlist.isobj) {
						Pft& pft=pftlist.getobj();
						if(pft.landcover==i) {
							stand.pft[pft.id].active=true;
						}
						pftlist.nextobj();
					}
				}
			}
//		}
	}
}

void harvest_natural(double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,double& cmass_debt,
	double& litter_leaf,double& litter_root,double& litter_wood,double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{
	double harvest=0.0;
	double residue_outtake=0.0;
	bool alive=indiv.alive;

	if(alive && cmass_root>0.0)
		litter_root+=cmass_root;			//all root carbon goes to litter
	cmass_root=0.0;

	if(alive && (cmass_sap+cmass_heart-cmass_debt)>0.0)						// Only wood currently harvested in this function !
	{	//indiv.pft.harv_eff = Bondeau's removal = 0.7 for tree wood, 0 for the rest
		harvest=indiv.pft.harv_eff*(cmass_sap+cmass_heart-cmass_debt);		//harvested products

		if(ifslowharvestpool)
		{
			harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;	//harvested products not consumed (oxidized) this year put into patchpft.harvested_products_slow
			harvest=harvest*(1-indiv.pft.harvest_slow_frac);
		}

		acflux_harvest+=harvest;							//harvested products consumed (oxidized) this year put into patch.fluxes.acflux_harvest, not litter pool !

		cmass_sap=(1-indiv.pft.harv_eff)*cmass_sap;			//unharvested parts of the plant
		cmass_heart=(1-indiv.pft.harv_eff)*cmass_heart;
		cmass_debt=(1-indiv.pft.harv_eff)*cmass_debt;		// ????
	}

	residue_outtake=indiv.pft.res_outtake*(cmass_sap+cmass_heart-cmass_debt+cmass_leaf);
	acflux_harvest+=residue_outtake;																//removed residues

	litter_leaf+=cmass_leaf*(1-indiv.pft.res_outtake);												//not removed residues
	litter_wood+=(cmass_sap+cmass_heart-cmass_debt)*(1-indiv.pft.res_outtake);						//not removed residues

	cmass_sap=cmass_heart=cmass_debt=cmass_leaf=0.0;
}

void landcover_dynamics(Gridcell& gridcell,Pftlist& pftlist)
{	// Called first day of the year if run_landcover is set.
	bool present;
	int i, j;	
	landcovertype landcover;
	double landcoverfrac_change[NLANDCOVERTYPES]={0.0};
//	double cropfrac_change[NCROPSTANDS_MAX]={0.0};	
	double cropfrac_sum_old=0.0;
//	double cropstand_change[NCROPSTANDS_MAX]={0.0};

	gridcell.LC_updated=false;

/////////////////////////////////////////////////////////////
//Landcover and cft fraction update (from updated landcoverfrac):

//Save old fraction values:
	for(i=0;i<NLANDCOVERTYPES;i++)
		gridcell.landcoverfrac_old[i]=gridcell.landcoverfrac[i];
//	for(i=0;i<NCROPSTANDS_MAX;i++)
//		cropfrac_sum_old+=gridcell.cftfrac_old[i]=gridcell.cftfrac[i];

//Get new gridcell.landcoverfrac and/or gridcell.cftfrac from LUdata and CFTdata.
	if(!all_fracs_const)
		getlandcover(gridcell,pftlist);	
	else return;

	double changeLC=0.0;
	double change_crop=0.0;
	double change_stand=0.0;
	double transferred_fraction=0.0;
	double receiving_fraction=0.0;

	if(!lcfrac_fixed)
	{
		for(i=0;i<NLANDCOVERTYPES;i++)
		{
			landcoverfrac_change[i]=gridcell.landcoverfrac[i]-gridcell.landcoverfrac_old[i];
			changeLC+=fabs(landcoverfrac_change[i])/2.0;
//			if(i!=CROPLAND)											//Landcovers with only one stand.
			{
				if(landcoverfrac_change[i]<0.0)
					transferred_fraction-=landcoverfrac_change[i];
				if(landcoverfrac_change[i]>0.0)
					receiving_fraction+=landcoverfrac_change[i];
				change_stand+=fabs(landcoverfrac_change[i])/2.0;
			}
		}
	}
/*	if(run[CROPLAND] &&!cftfrac_fixed)
	{
		for(i=0;i<NCROPSTANDS_MAX;i++)
		{
			cropfrac_change[i]=gridcell.cftfrac[i]-gridcell.cftfrac_old[i];
			cropstand_change[i]=gridcell.cftfrac[i]*gridcell.landcoverfrac[CROPLAND]-gridcell.cftfrac_old[i]*gridcell.landcoverfrac_old[CROPLAND];

			if(cropstand_change[i]<0.0)
				transferred_fraction-=cropstand_change[i];
			if(cropstand_change[i]>0.0)
				receiving_fraction+=cropstand_change[i];

			if(cropfrac_sum_old!=0.0)
			{
				change_crop+=fabs(cropfrac_change[i])/2.0;
			}
			else
			{
				change_crop+=fabs(cropfrac_change[i]);
			}
			change_stand+=fabs(cropstand_change[i])/2.0;	//cropfrac_sum_old+gridcell.landcoverfrac[NATURAL] should never be 0.0
		}
	}
*/
// If no changes, do nothing.
	if(changeLC<0.00001 && change_crop<0.00001)
		return;
	else if(fabs(transferred_fraction-receiving_fraction)>0.0001 || fabs(change_stand-receiving_fraction)>0.0001)
			fail("Transferred landcover fractions not balanced !\n");
////////////////////////////////////////////////////////////

#define cropLUchangeCtransfer
#if defined cropLUchangeCtransfer

	double *transfer_litter_leaf, *transfer_litter_wood, *transfer_litter_root, *transfer_litter_repr, *transfer_harvested_products_slow;

	transfer_litter_leaf=transfer_litter_wood=transfer_litter_root=transfer_litter_repr=transfer_harvested_products_slow=NULL;

	transfer_litter_leaf=new double[npft];
	transfer_litter_wood=new double[npft];
	transfer_litter_root=new double[npft];
	transfer_litter_repr=new double[npft];

	transfer_harvested_products_slow=new double[npft];

	double transfer_acflux_harvest=0.0;

	double transfer_cpool_fast=0.0;
	double transfer_cpool_slow=0.0;
	double transfer_wcont[NSOILLAYER]={0.0};
	double transfer_decomp_litter_mean=0.0;
	double transfer_k_soilfast_mean=0.0;
	double transfer_k_soilslow_mean=0.0;

	memset(transfer_litter_leaf,0,sizeof(double)*npft);
	memset(transfer_litter_wood,0,sizeof(double)*npft);
	memset(transfer_litter_root,0,sizeof(double)*npft);
	memset(transfer_litter_repr,0,sizeof(double)*npft);
	memset(transfer_harvested_products_slow,0,sizeof(double)*npft);

//Keep track of carbon and water in lost areas.
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		double scale;

		Stand& stand=gridcell.getobj();

		// Reset fluxes. NB. landcover_dynamics() is called before dailyaccounting_patch() on date.day==0 && date.year>=nyear_spinup !
		stand.firstobj();
		while(stand.isobj) //Loop through Patches
		{
			Patch& patch=stand.getobj();
		
			patch.fluxes.acflux_harvest=0.0;

			stand.nextobj();
		}

//		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]<0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)
		if(landcoverfrac_change[stand.landcover]<0.0)
		{
//			if(stand.landcover!=CROPLAND)														
			{
				scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
			}
/*			else if(stand.landcover==CROPLAND)
			{
				scale=-cropstand_change[stand.cftid]/receiving_fraction/(double)stand.nobj;
			}	
*/
			stand.firstobj();
			while(stand.isobj) //Loop through Patches
			{
				Patch& patch=stand.getobj();

				Vegetation& vegetation=patch.vegetation;
				vegetation.firstobj();
				while(vegetation.isobj)
				{
					double cmass_leaf_cp=0.0, cmass_root_cp=0.0, cmass_sap_cp=0.0, cmass_heart_cp=0.0, cmass_debt_cp=0.0, cmass_ho_cp=0.0, cmass_agpool_cp=0.0, cmass_plant_cp=0.0;//bugfix 101103
					double litter_leaf_cp, litter_root_cp, litter_wood_cp, litter_repr_cp;
					double acflux_harvest_cp;
					double harvested_products_slow_cp;

					Individual& indiv=vegetation.getobj();
					Patchpft& patchpft=patch.pft[indiv.pft.id];

					cmass_leaf_cp=indiv.cmass_leaf;
					cmass_root_cp=indiv.cmass_root;
					cmass_sap_cp=indiv.cmass_sap;
					cmass_heart_cp=indiv.cmass_heart;
					cmass_debt_cp=indiv.cmass_debt;

/*					if(indiv.pft.landcover==CROPLAND)
					{
						cmass_ho_cp=indiv.cropindiv->cmass_ho;
						cmass_agpool_cp=indiv.cropindiv->cmass_agpool;
						cmass_plant_cp=indiv.cropindiv->cmass_plant;
					}
*/
					litter_leaf_cp=patchpft.litter_leaf;
					litter_root_cp=patchpft.litter_root;
					litter_wood_cp=patchpft.litter_wood;
					litter_repr_cp=patchpft.litter_repr;

					acflux_harvest_cp=patch.fluxes.acflux_harvest;	//flux är nollställd
					harvested_products_slow_cp=patch.pft[indiv.pft.id].harvested_products_slow;

//Harvest of transferred areas:
/*					if(indiv.pft.landcover==CROPLAND)
						harvest_crop(cmass_plant_cp,cmass_leaf_cp,cmass_root_cp,cmass_ho_cp,cmass_agpool_cp,
						litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);
					else if(patch.stand.landcover!=CROPLAND)												
*/						harvest_natural(cmass_leaf_cp,cmass_root_cp,cmass_sap_cp,cmass_heart_cp,cmass_debt_cp,	//kolla vad som händer här, både för träd och gräs !
						litter_leaf_cp,litter_root_cp,litter_wood_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);

					gridcell.LC_updated=true;

					//In case any vegetation carbon left: (eg. cmass_root for CC3G/CC4G)
					if((cmass_leaf_cp+cmass_root_cp+cmass_sap_cp+cmass_heart_cp-cmass_debt_cp+cmass_ho_cp)!=0.0)
					{
						litter_leaf_cp+=cmass_leaf_cp;
						litter_root_cp+=cmass_root_cp;
						litter_wood_cp+=cmass_sap_cp+cmass_heart_cp-cmass_debt_cp;
/*
						if(indiv.pft.aboveground_ho)
							litter_leaf_cp+=cmass_ho_cp;
						else
							litter_root_cp+=cmass_ho_cp;
*/					}

					transfer_litter_leaf[indiv.pft.id]+=litter_leaf_cp*scale;
					transfer_litter_root[indiv.pft.id]+=litter_root_cp*scale;
					transfer_litter_wood[indiv.pft.id]+=litter_wood_cp*scale;
					transfer_litter_repr[indiv.pft.id]+=litter_repr_cp*scale;

					transfer_acflux_harvest+=acflux_harvest_cp*scale;

					if(ifslowharvestpool)
						transfer_harvested_products_slow[indiv.pft.id]+=harvested_products_slow_cp*scale;

					vegetation.nextobj();
				}

//sum litter C:
				transfer_cpool_fast+=patch.soil.cpool_fast*scale;
				transfer_cpool_slow+=patch.soil.cpool_slow*scale;

//sum wcont:
				for(i=0;i<NSOILLAYER;i++)
				{
					transfer_wcont[i]+=patch.soil.wcont[i]*scale;
				}

				transfer_decomp_litter_mean+=patch.soil.decomp_litter_mean*scale;
				transfer_k_soilfast_mean+=patch.soil.k_soilfast_mean*scale;
				transfer_k_soilslow_mean+=patch.soil.k_soilslow_mean*scale;

				stand.nextobj();
			}
		}
		gridcell.nextobj();
	}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create and kill stands:

// landcover dynamics (from updated landcoverfrac):	
	if(!lcfrac_fixed && changeLC>0.0)
	{
		for(int i=0;i<NLANDCOVERTYPES;i++)	//For all landcover types without subclasses
		{
//			if(i!=CROPLAND)
			{
				if(run[i])
				{
					if(gridcell.landcoverfrac_old[i]==0.0 && gridcell.landcoverfrac[i]>0.0)
					{
						landcover=(landcovertype)i;
						Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);

						pftlist.firstobj();
						while (pftlist.isobj) 
						{
							Pft& pft=pftlist.getobj();
							if(pft.landcover==i)
							{
								stand.pft[pft.id].active=true;
							}
							pftlist.nextobj();
						}
					}
					else if(gridcell.landcoverfrac_old[i]>0.0 && gridcell.landcoverfrac[i]==0.0)
					{
						gridcell.firstobj();
						while (gridcell.isobj) //Loop through stands:
						{
							Stand& stand=gridcell.getobj();
							if(stand.landcover==i)
							{
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

//Crop stand dynamics to be put here.

//update C-pools for receiving stands:
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		Stand& stand=gridcell.getobj();

//		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]>0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]>0.0)
		if(landcoverfrac_change[stand.landcover]>0.0)
		{
			double old_frac, added_frac, new_frac;

//			if(stand.landcover!=CROPLAND)					
			{
				old_frac=gridcell.landcoverfrac_old[stand.landcover];
				added_frac=landcoverfrac_change[stand.landcover];
				new_frac=gridcell.landcoverfrac[stand.landcover];
			}
/*			else if(stand.landcover==CROPLAND)
			{
				old_frac=gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.cftid];
				added_frac=cropstand_change[stand.cftid];
				new_frac=gridcell.landcoverfrac[CROPLAND]*gridcell.cftfrac[stand.cftid];
			}	
*/
#ifdef cropLUchangeCtransfer
			stand.firstobj();
			while(stand.isobj) //Loop through Patches
			{
				Patch& patch=stand.getobj();
//add litter C:
				for (i=0;i<npft;i++) 
				{
					Patchpft& patchpft=patch.pft[i];

					patchpft.litter_leaf=(patchpft.litter_leaf*old_frac+transfer_litter_leaf[i]*added_frac)/new_frac;
					patchpft.litter_wood=(patchpft.litter_wood*old_frac+transfer_litter_wood[i]*added_frac)/new_frac;
					patchpft.litter_root=(patchpft.litter_root*old_frac+transfer_litter_root[i]*added_frac)/new_frac;
					patchpft.litter_repr=(patchpft.litter_repr*old_frac+transfer_litter_repr[i]*added_frac)/new_frac;

					if(ifslowharvestpool)
						patchpft.harvested_products_slow=(patchpft.harvested_products_slow*old_frac+transfer_harvested_products_slow[i]*added_frac)/new_frac;
				}

//add soil C:
				patch.soil.cpool_fast=(patch.soil.cpool_fast*old_frac+transfer_cpool_fast*added_frac)/new_frac;
				patch.soil.cpool_slow=(patch.soil.cpool_slow*old_frac+transfer_cpool_slow*added_frac)/new_frac;

//other soil stuff:
				for(i=0;i<NSOILLAYER;i++)
					patch.soil.wcont[i]=(patch.soil.wcont[i]*old_frac+transfer_wcont[i]*added_frac)/new_frac;

				patch.soil.decomp_litter_mean=(patch.soil.decomp_litter_mean*old_frac+transfer_decomp_litter_mean*added_frac)/new_frac;
				patch.soil.k_soilfast_mean=(patch.soil.k_soilfast_mean*old_frac+transfer_k_soilfast_mean*added_frac)/new_frac;
				patch.soil.k_soilslow_mean=(patch.soil.k_soilslow_mean*old_frac+transfer_k_soilslow_mean*added_frac)/new_frac;
//add fluxes:
				patch.fluxes.acflux_harvest=(patch.fluxes.acflux_harvest*old_frac+transfer_acflux_harvest*added_frac)/new_frac;

				stand.nextobj();
			}
#endif
		}
		gridcell.nextobj();
	}
#if defined cropLUchangeCtransfer
	if(transfer_litter_leaf) delete[] transfer_litter_leaf;
	if(transfer_litter_wood) delete[] transfer_litter_wood;
	if(transfer_litter_root) delete[] transfer_litter_root;
	if(transfer_litter_repr) delete[] transfer_litter_repr;
	if(transfer_harvested_products_slow) delete[] transfer_harvested_products_slow;
#endif
}

Stand::Stand(int i, Gridcell& gc,landcovertype landcoverX,Pftlist& pftlist):id(i),gridcell(gc),landcover(landcoverX),frac(1.0) {

		// Constructor: initialises reference member of climate and
		// builds list array of Standpft objects
		
	int p;
	int npatchL;

	for(p=0;p<pftlist.nobj;p++) {
		pft.createobj(pftlist[p]);
	}


	if(landcover==CROPLAND || landcover==PASTURE || landcover==URBAN || landcover==PEATLAND) {
		npatchL=1;
	}
	else if(landcover==NATURAL || landcover==FOREST) {
		npatchL=npatch;
	}

	for (p=0;p<npatchL;p++) {
		createobj(*this,pftlist,gc.soiltype);
	}

	first_year=date.year;
}

double Stand::get_gridcell_fraction() const {
	return frac*gridcell.landcoverfrac[landcover];
}

double Stand::get_landcover_fraction() const {
	return frac;
}

void Stand::set_landcover_fraction(double fraction) {
	frac = fraction;
}

Individual::Individual(int i,Pft& p,Vegetation& v):id(i),pft(p),vegetation(v) {

	anpp=0.0;
	fpc=0.0;
	densindiv=0.0;
	cmass_leaf=0.0;
	cmass_root=0.0;
	cmass_sap=0.0;
	cmass_heart=0.0;
	cmass_debt=0.0;
	wscal=1.0;
	phen=0.0;
	aphen=0.0;
	deltafpc=0.0;
	fpar_wstress=0.0;
	assim=0.0;

	// guess2008 - additional initialisation
	age=0.0;
	fpar=0.0;
	aphen_raingreen=0;
	demand=0.0;
	supply=0.0;
	intercep=0.0;
	phen_mean=0.0;
	temp_wstress = 0.0;
	par_wstress = 0.0;
	daylength_wstress = 0.0;
	co2_wstress = 0.0; 
	nday_wstress = 0; 
	ifwstress = false;
	lai = 0.0;
	lai_layer = 0.0;
	lai_indiv = 0.0;
	alive = false;

	int m;
	for (m=0;m<12;m++) {
		mnpp[m]=mlai[m]=mgpp[m]=mra[m]=0.0;
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// THE FRAMEWORK
// The 'mission control' of the model, responsible for maintaining the primary model
// data structures and containing all explicit loops through space (grid cells/stands)
// and time (days and years).

int framework(int argc,char* argv[]) {

	bool dogridcell;
	int p;

	// The one and only linked list of Pft objects	
	Pftlist pftlist;

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(argc,argv,pftlist);

	// Assume there is at least one grid cell to simulate
	dogridcell=true;

	while (dogridcell) {

		// START OF LOOP THROUGH GRID CELLS

		// Initialise global variable date
		// (argument nyear not used in this implementation)
		date.init(1);

		// Create and initialise a new Gridcell object for each locality
		Gridcell gridcell(pftlist);	

		// Call input/output to obtain latitude and soil driver data for this grid cell.
		// Function getgridcell returns false if no further grid cells remain to be simulated

		if (getgridcell(gridcell)) {

			// Initialise certain climate and soil drivers
			gridcell.climate.initdrivers(gridcell.climate.lat);

			if(run_landcover) {
				//Read static landcover and cft fraction data from ins-file and/or from data files for the spinup peroid and create stands.
				landcover_init(gridcell,pftlist);
			}
			
			// Call input/output to obtain climate, insolation and CO2 for this
			// day of the simulation. Function getclimate returns false if last year
			// has already been simulated for this grid cell

			while (getclimate(gridcell)) {

				// START OF LOOP THROUGH SIMULATION DAYS

				// Update daily climate drivers etc
				dailyaccounting_gridcell(gridcell,pftlist);

				// Calculate daylength, insolation and potential evapotranspiration
				daylengthinsoleet(gridcell.climate);

				if(run_landcover && date.day==0) {
				// Update dynamic landcover and crop fraction data during historical period and create/kill stands.
					if(date.year>=nyear_spinup)
						landcover_dynamics(gridcell,pftlist);
				}

				gridcell.firstobj();
				while (gridcell.isobj) {

					// START OF LOOP THROUGH STANDS

					Stand& stand=gridcell.getobj();

					dailyaccounting_stand(stand,pftlist);

					stand.firstobj();
					while (stand.isobj) {
						// START OF LOOP THROUGH PATCHES

						// Get reference to this patch
						Patch& patch=stand.getobj();
						// Update daily soil drivers including soil temperature
						dailyaccounting_patch(patch,pftlist);
						// Leaf phenology for PFTs and individuals
						leaf_phenology(patch,gridcell.climate);
						// Photosynthesis, respiration, evapotranspiration
						canopy_exchange(patch);
						// Soil water accounting, snow pack accounting
						soilwater(gridcell.climate,patch);
						// Soil organic matter and litter dynamics
						som_dynamics(patch);

						if (date.islastday && date.islastmonth) {

							// LAST DAY OF YEAR
							// Tissue turnover, allocation to new biomass and reproduction,
							// updated allometry
							growth(stand,patch);
						}
						stand.nextobj();
					}// End of loop through patches

					if (date.islastday && date.islastmonth)
					{
						// LAST DAY OF YEAR
						stand.firstobj();
						while (stand.isobj) {

							// For each patch ...
							Patch& patch=stand.getobj();
							// Establishment, mortality and disturbance by fire
							vegetation_dynamics(stand,patch,pftlist);
							stand.nextobj();
						}
					}

					gridcell.nextobj();			
				}	// End of loop through stands

				if (date.islastday && date.islastmonth) {
					// LAST DAY OF YEAR
					// Call input/output module to output results for end of year
					// or end of simulation for this grid cell
					outannual(gridcell,pftlist);

					// Check whether to abort
					if (abort_request_received()) {
						termio();
						return 99;
					}
				}

				// Advance timer to next simulation day
				date.next();

				// End of loop through simulation days
			}//while (getclimate())
		}//if getgridcell()
		else dogridcell=false; // no more grid cells to simulate

		int test = 0;

		// End of loop through grid cells
	}

	// Call to input/output module to perform any necessary clean up
	termio();

	// END OF SIMULATION

	return 0;
}
