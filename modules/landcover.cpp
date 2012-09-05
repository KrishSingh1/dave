///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.cpp
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "landcover.h"
#include "guessio.h"
#include "canexch.h"

#define DYNAMIC_PHU					//Calculation of potential heat units according to local climate.
#define MAXHUTEMP					//30 degree limit for heat unit summation
#define SD_TEMP_WINDOW				//Uses sowing window for temperature-dependent sowing.
#define FAO_PREC_SEASONS			//If prec/pet always above 1.0, default sowing date is always used.
#define IRRIGATED_USE_TEMP_SDATE	//Use temperature-dependent sowing date for irrigated crops at site with PRECTEMP seasonality.
//#define DELAYED_SEEDCARBON		//Seed carbon allocation to leaves and roots are done over a 10-day period.


//Functions facilitating handling time periods spanning newyear:
bool dayinperiod(int day, int start, int end)	//Intended for use with short windows (start-end, eg. 30 days), but can be used with longer with caution.
{
	bool acrossnewyear=false;

	if(start<0 || end<0)	//111117: a negative value should not be a valid day
		return false;

	if(start>end)
		acrossnewyear=true;

	if(day>=start && day<=end && !acrossnewyear || (day>=start || day<=end) && acrossnewyear)
		return true;
	else
		return false;
}

int stepfromdate(int day, int step)
{
	bool acrossnewyear=false;

	if(day+step>0)
		return (day+step)%365;
	else if(day+step<0)
		return day+step+365;
	else
		return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////  Landcover stand dynamics and C-partitioning  /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void landcover_init(Gridcell& gridcell) {
	landcovertype landcover;

	getlandcover(gridcell);		//Gets gridcell.landcoverfrac from landcover input file(s) or ins-file.

	for(int i=0;i<NLANDCOVERTYPES;i++) { //For all landcover types without subclasses
		if(i!=CROPLAND) {				
			if(run[i]) {
				if(gridcell.landcoverfrac[i]>0.0) {
					landcover=(landcovertype)i;
					Stand& stand=gridcell.createobj(gridcell,landcover);
					stand.set_gridcell_fraction(gridcell.landcoverfrac[i]);

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
		}
	}

	if(run[CROPLAND])
	{
		if(gridcell.landcoverfrac[CROPLAND]>0.0)	//landcoverfrac[CROP] overrides cftfrac[]
		{
			pftlist.firstobj();
			while (pftlist.isobj) 
			{
				Pft& pft=pftlist.getobj();
				if(pft.landcover==CROPLAND)
				{
					//Is this crop PFT present in the gridcell ? 

					if(gridcell.cftfrac[pft.cftid]>0.0)
					{
						landcover=CROPLAND;
						Stand& stand=gridcell.createobj(gridcell,landcover);
						stand.pftid=pft.id;
						stand.cftid=pft.cftid;
						stand.set_gridcell_fraction(gridcell.cftfrac[pft.cftid]*gridcell.landcoverfrac[CROPLAND]);

						stand.pft[pft.id].active=true;	//101213

						if(pft.hydrology==IRRIGATED)
						{
							stand.isirrigated=true;
						}

						if(pft.intercrop==NATURALGRASS && ifintercropgrass)	
						{
							stand.hasgrassintercrop=true;

							for(int i=0;i<pftlist.nobj;i++)		//101213
							{
								if(pftlist[i].isintercropgrass)
									stand.pft[pftlist[i].id].active=true;
							}
						}

					// Set crop cycle dates to default values.
						for(int i=0;i<stand.nobj;i++)
						{
							stand[i].pft[pft.id].set_cropphen()->sdate=stand.gridcell.pft[pft.id].sdate_default;
							stand[i].pft[pft.id].set_cropphen()->hlimitdate=stand.gridcell.pft[pft.id].hlimitdate_default;
	
							if(pft.phenology==ANY)
								stand[i].pft[stand.pftid].set_cropphen()->growingseason=true;
							else if(pft.phenology==CROPGREEN)
							{
								stand[i].pft[pft.id].set_cropphen()->eicdate=stand[i].pft[pft.id].get_cropphen()->sdate-15;
								if(stand[i].pft[pft.id].get_cropphen()->eicdate<0)
									stand[i].pft[pft.id].set_cropphen()->eicdate=365+stand[i].pft[pft.id].get_cropphen()->sdate-15;
							}
						}

					}
				}
				pftlist.nextobj(); // ... on to next PFT
			}
		}
	}
}

void landcover_dynamics(Gridcell& gridcell)
{	// Called first day of the year if run_landcover is set.
	bool present;
	int i, j;	
	landcovertype landcover;
	double landcoverfrac_change[NLANDCOVERTYPES];
	double cropfrac_change[NCROPSTANDS_MAX];	
	double cropfrac_sum_old=0.0;
	double cropstand_change[NCROPSTANDS_MAX];
	int nnaturalstands=0;

	bool LCchangeCtransfer=true;

	memset(landcoverfrac_change,0,NLANDCOVERTYPES*sizeof(double));
	memset(cropfrac_change,0,NCROPSTANDS_MAX*sizeof(double));
	memset(cropstand_change,0,NCROPSTANDS_MAX*sizeof(double));

	gridcell.LC_updated=false;

/////////////////////////////////////////////////////////////
//Landcover and cft fraction update (from updated landcoverfrac):

//Save old fraction values:									
	for(i=0;i<NLANDCOVERTYPES;i++)
		gridcell.landcoverfrac_old[i]=gridcell.landcoverfrac[i];
	for(i=0;i<NCROPSTANDS_MAX;i++)
		cropfrac_sum_old+=gridcell.cftfrac_old[i]=gridcell.cftfrac[i];

//Get new gridcell.landcoverfrac and/or gridcell.cftfrac from LUdata and CFTdata.
	if(!all_fracs_const)					
		getlandcover(gridcell);	
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
			if(i!=CROPLAND)									
			{
				if(landcoverfrac_change[i]<0.0)
					transferred_fraction-=landcoverfrac_change[i];
				if(landcoverfrac_change[i]>0.0)
					receiving_fraction+=landcoverfrac_change[i];
				change_stand+=fabs(landcoverfrac_change[i])/2.0;
			}
		}
	}

	if(run[CROPLAND] && (!cftfrac_fixed || !lcfrac_fixed))	
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

// If no changes, do nothing.
	if(changeLC<0.00001 && change_crop<0.00001)
		return;
	else 
	{
		if(fabs(transferred_fraction-receiving_fraction)>0.0001 || fabs(change_stand-receiving_fraction)>0.0001)
		{
			if(run[CROPLAND] && run[NATURAL])
				fail("Transferred landcover fractions not balanced !\n");
			else
			{
				LCchangeCtransfer=false;
if(!SUPPRESSLARGEOUTPUT)
					dprintf("Transferred landcover fractions not balanced !\nLandcover change carbon flux not calculated.\n");
			}
		}
	}
////////////////////////////////////////////////////////////

	double *transfer_litter_leaf, *transfer_litter_wood, *transfer_litter_root, *transfer_litter_repr, *transfer_harvested_products_slow;

	transfer_litter_leaf=transfer_litter_wood=transfer_litter_root=transfer_litter_repr=transfer_harvested_products_slow=NULL;

	double transfer_acflux_harvest=0.0;

	double transfer_cpool_fast=0.0;
	double transfer_cpool_slow=0.0;
	double transfer_wcont[NSOILLAYER];
	double transfer_decomp_litter_mean=0.0;
	double transfer_k_soilfast_mean=0.0;
	double transfer_k_soilslow_mean=0.0;

	memset(transfer_wcont,0,NSOILLAYER*sizeof(double));

	if(LCchangeCtransfer)	//stand.frac not updated if Ctransfer code not read !
	{
		transfer_litter_leaf=new double[npft];
		transfer_litter_wood=new double[npft];
		transfer_litter_root=new double[npft];
		transfer_litter_repr=new double[npft];

		transfer_harvested_products_slow=new double[npft];

		memset(transfer_litter_leaf,0,sizeof(double)*npft);
		memset(transfer_litter_wood,0,sizeof(double)*npft);
		memset(transfer_litter_root,0,sizeof(double)*npft);
		memset(transfer_litter_repr,0,sizeof(double)*npft);
		memset(transfer_harvested_products_slow,0,sizeof(double)*npft);

//Keep track of carbon and water in lost areas.

#ifdef multiple_natural_stands
		if(landcoverfrac_change[NATURAL]<0.0)	// Find the age order of natural stands.
		{

			gridcell.firstobj();
			while (gridcell.isobj) //Loop through stands:
			{
				Stand& stand=gridcell.getobj();

				if(stand.landcover==NATURAL)
				{
					stand.natural_frac_change=0.0;
					nnaturalstands++;
				}

				gridcell.nextobj();
			}

			if(nnaturalstands>1)
			{
				double natural_change_remain=landcoverfrac_change[NATURAL];

				int i=0, index;

				for(i=0;i<gridcell.nobj;i++)
				{
					int j;

					//convert equal percentage of area from all stands
					if((landcoverfrac_change[CROPLAND]>0.0 || landcoverfrac_change[PASTURE]>0.0 || landcoverfrac_change[URBAN]>0.0 || landcoverfrac_change[PEATLAND]>0.0) && landcoverfrac_change[FOREST]>0.0)
						index=i;
					//convert area from youngest stands first
  					else if(landcoverfrac_change[CROPLAND]>0.0 || landcoverfrac_change[PASTURE]>0.0 || landcoverfrac_change[URBAN]>0.0 || landcoverfrac_change[PEATLAND]>0.0)
						index=gridcell.nobj-1-i;
					//convert area from oldest stands first
					else if(landcoverfrac_change[FOREST]>0.0)	
						index=i;


					Stand& stand=gridcell[index];	

					if(stand.landcover==NATURAL)
					{
						//convert equal areas from all stands
						if((landcoverfrac_change[CROPLAND]>0.0 || landcoverfrac_change[PASTURE]>0.0 || landcoverfrac_change[URBAN]>0.0 || landcoverfrac_change[PEATLAND]>0.0) && landcoverfrac_change[FOREST]>0.0)
						{
							stand.natural_frac_change=landcoverfrac_change[NATURAL]*stand.get_gridcell_fraction()/gridcell.landcoverfrac[NATURAL];
							stand.set_gridcell_fraction(stand.get_gridcell_fraction()+stand.natural_frac_change);
							break;
						}
						else
						{		
							if(stand.get_gridcell_fraction()>0.0)
							{		
								if(stand.get_gridcell_fraction()>=-natural_change_remain)	//all natural landcover decrease is taken from this stand
								{
									stand.natural_frac_change=natural_change_remain;
									stand.set_gridcell_fraction(stand.get_gridcell_fraction()+stand.natural_frac_change);
									natural_change_remain=0.0;
									break;
								}
								else									//more stands will have to be reduced
								{
									stand.natural_frac_change=-stand.get_gridcell_fraction();
									natural_change_remain+=stand.get_gridcell_fraction();
									stand.set_gridcell_fraction(0.0);	//will be killed below
								}				
							}
						}
					}
				}
			}		
		}
#endif

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

#ifdef multiple_natural_stands
			if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL && landcoverfrac_change[stand.landcover]<0.0						
				|| stand.landcover==NATURAL && landcoverfrac_change[NATURAL]<0.0 && (nnaturalstands==1 || stand.natural_frac_change<0.0)
				|| stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)
#else
			if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]<0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)																					//2b
#endif
			{
//All landcovers that only have one stand:
#ifdef multiple_natural_stands
				if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL)
#else
				if(stand.landcover!=CROPLAND)
#endif
				{
					scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;			
					stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
				}
//Landcovers that may have several stands:
#ifdef multiple_natural_stands
				else if(stand.landcover==NATURAL)
				{
					if(nnaturalstands>1)
					{
						scale=-stand.natural_frac_change/receiving_fraction/(double)stand.nobj;
					}
					else
					{
						scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
						stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
					}
				}
#endif			
				else if(stand.landcover==CROPLAND)
				{
					scale=-cropstand_change[stand.cftid]/receiving_fraction/(double)stand.nobj;
					stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
				}	
	
				stand.firstobj();
				while(stand.isobj) //Loop through Patches
				{
					Patch& patch=stand.getobj();

//sum original litter:
					for(int n=0;n<npft;n++)
					{
						transfer_litter_leaf[n]+=patch.pft[n].litter_leaf*scale;
						transfer_litter_root[n]+=patch.pft[n].litter_root*scale;
						transfer_litter_wood[n]+=patch.pft[n].litter_wood*scale;
						transfer_litter_repr[n]+=patch.pft[n].litter_repr*scale;
						if(ifslowharvestpool)
							transfer_harvested_products_slow[n]+=patch.pft[n].harvested_products_slow*scale;
					}
					transfer_acflux_harvest+=patch.fluxes.acflux_harvest*scale;

					Vegetation& vegetation=patch.vegetation;
					vegetation.firstobj();
					while(vegetation.isobj)
					{
						double cmass_leaf_cp=0.0, cmass_root_cp=0.0, cmass_sap_cp=0.0, cmass_heart_cp=0.0, cmass_debt_cp=0.0, cmass_ho_cp=0.0, cmass_agpool_cp=0.0, cmass_plant_cp=0.0;//bugfix 101103
						double litter_leaf_cp=0.0, litter_root_cp=0.0, litter_wood_cp=0.0, litter_repr_cp=0.0;
						double acflux_harvest_cp=0.0;
						double harvested_products_slow_cp=0.0;

						Individual& indiv=vegetation.getobj();
						Patchpft& patchpft=patch.pft[indiv.pft.id];

						cmass_leaf_cp=indiv.cmass_leaf;
						cmass_root_cp=indiv.cmass_root;
						cmass_sap_cp=indiv.cmass_sap;
						cmass_heart_cp=indiv.cmass_heart;
						cmass_debt_cp=indiv.cmass_debt;

						if(indiv.pft.landcover==CROPLAND)
						{
							cmass_ho_cp=indiv.cropindiv->cmass_ho;
							cmass_agpool_cp=indiv.cropindiv->cmass_agpool;
							cmass_plant_cp=indiv.cropindiv->cmass_plant;
						}
	
	//Harvest of transferred areas:
						if(indiv.pft.landcover==CROPLAND)
							harvest_crop(cmass_plant_cp,cmass_leaf_cp,cmass_root_cp,cmass_ho_cp,cmass_agpool_cp,
							litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);
						else if(indiv.pft.landcover==PASTURE)
						{
							harvest_pasture(cmass_leaf_cp,cmass_root_cp,
								litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp, indiv);
						}
						else												
							harvest_natural(cmass_leaf_cp,cmass_root_cp,cmass_sap_cp,cmass_heart_cp,cmass_debt_cp,
							litter_leaf_cp,litter_root_cp,litter_wood_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);

						gridcell.LC_updated=true;

						//In case any vegetation carbon left: (eg. cmass_root for CC3G/CC4G)
						if((cmass_leaf_cp+cmass_root_cp+cmass_sap_cp+cmass_heart_cp-cmass_debt_cp+cmass_ho_cp)!=0.0)
						{
							litter_leaf_cp+=cmass_leaf_cp;
							litter_root_cp+=cmass_root_cp;
							litter_wood_cp+=cmass_sap_cp+cmass_heart_cp-cmass_debt_cp;
	
							if(indiv.pft.landcover==CROPLAND)
							{
								if(indiv.pft.aboveground_ho)
									litter_leaf_cp+=cmass_ho_cp;
								else
									litter_root_cp+=cmass_ho_cp;
							}
						}

						transfer_litter_leaf[indiv.pft.id]+=litter_leaf_cp*scale;
						transfer_litter_root[indiv.pft.id]+=litter_root_cp*scale;
						transfer_litter_wood[indiv.pft.id]+=litter_wood_cp*scale;
						transfer_litter_repr[indiv.pft.id]+=litter_repr_cp*scale;

						transfer_acflux_harvest+=acflux_harvest_cp*scale;

						if(stand.landcover==NATURAL)
						{
							double change_frac;
							if(nnaturalstands>1)
								change_frac=stand.natural_frac_change;
							else
								change_frac=landcoverfrac_change[NATURAL];
//							gridcell.acflux_landuse_change+=-acflux_harvest_cp*change_frac/(double)stand.nobj;
						}

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
	}//if(LCchangeCtransfer)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create and kill stands:

// landcover dynamics (from updated landcoverfrac):	
	if(!lcfrac_fixed && changeLC>0.0)
	{
		for(int i=0;i<NLANDCOVERTYPES;i++)
		{
			if(i!=CROPLAND)
			{
				if(run[i])
				{
					if(gridcell.landcoverfrac_old[i]==0.0 && gridcell.landcoverfrac[i]>0.0)
					{
						landcover=(landcovertype)i;
						Stand& stand=gridcell.createobj(gridcell,landcover);
						stand.set_gridcell_fraction(gridcell.landcoverfrac[i]);

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
#if defined multiple_natural_stands
					else if(i==NATURAL && landcoverfrac_change[i]>0.0)	//New NATURAL stand created from other landcover type.
					{
						landcover=(landcovertype)i;
						Stand& stand=gridcell.createobj(gridcell,landcover);
						stand.set_gridcell_fraction(landcoverfrac_change[i]);

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
					else if(i==NATURAL && nnaturalstands>1 && landcoverfrac_change[i]<0.0)	// Natural stand killed.
					{
						gridcell.firstobj();
						while (gridcell.isobj)
						{
							Stand& stand=gridcell.getobj();
							if(stand.landcover==i && stand.get_gridcell_fraction()==0)
							{
								gridcell.killobj();
							}
							else
								gridcell.nextobj();
						}
					}
#endif
				}
			}
		}
	}

//CFT stand dynamics (from updated cftfrac):
	if(run[CROPLAND] && (change_crop>0.0 || landcoverfrac_change[CROPLAND]!=0.0))
	{
		if(gridcell.landcoverfrac[CROPLAND]>0.0)	//landcoverfrac[CROP] overrides cftfrac[]
		{
			pftlist.firstobj();
			while (pftlist.isobj) 
			{
				Pft& pft=pftlist.getobj();
				if(pft.landcover==CROPLAND)
				{
					// Is this PFT already represented in a crop stand ?
					present=false;
					gridcell.firstobj();
					while (gridcell.isobj && !present) 
					{
						Stand& stand=gridcell.getobj();
						if (stand.landcover==CROPLAND && stand.pftid==pft.id)
							present=true;
						else
							gridcell.nextobj();
					}

					//Should this crop PFT be present in the gridcell this year ? 
					if(gridcell.cftfrac[pft.cftid]>0.0)
					{
						if(present)
						{
							Stand& stand=gridcell.getobj();
							stand.set_gridcell_fraction(gridcell.cftfrac[pft.cftid]*gridcell.landcoverfrac[CROPLAND]);
						}
						else
						{
							landcover=CROPLAND;
							Stand& stand=gridcell.createobj(gridcell,landcover);

							stand.pftid=pft.id;
							stand.cftid=pft.cftid;
							stand.set_gridcell_fraction(gridcell.cftfrac[pft.cftid]*gridcell.landcoverfrac[CROPLAND]);

							stand.pft[pft.id].active=true;

							if(pft.hydrology==IRRIGATED)
							{
								stand.isirrigated=true;
							}

							if(pft.intercrop==NATURALGRASS && ifintercropgrass)
							{
								stand.hasgrassintercrop=true;

								for(int i=0;i<pftlist.nobj;i++)	
								{
									if(pftlist[i].isintercropgrass)
										stand.pft[pftlist[i].id].active=true;
								}
							}

							// Set crop cycle dates to default values.
							for(int i=0;i<stand.nobj;i++)
							{
								stand[i].pft[pft.id].set_cropphen()->sdate=stand.gridcell.pft[pft.id].sdate_default;
								stand[i].pft[pft.id].set_cropphen()->hlimitdate=stand.gridcell.pft[pft.id].hlimitdate_default;

								if(pft.phenology==ANY)								
									stand[i].pft[stand.pftid].set_cropphen()->growingseason=true;
								else if(pft.phenology==CROPGREEN)
								{
									stand[i].pft[pft.id].set_cropphen()->eicdate=stand[i].pft[pft.id].get_cropphen()->sdate-15;
									if(stand[i].pft[pft.id].get_cropphen()->eicdate<0)
										stand[i].pft[pft.id].set_cropphen()->eicdate=365+stand[i].pft[pft.id].get_cropphen()->sdate-15;
								}
							}
						}
					}
					else
					{
						if(present)
						{
							Stand& stand=gridcell.getobj();
							gridcell.killobj();
						}
					}
				}
				pftlist.nextobj(); // ... on to next PFT
			}
		}
		else if(gridcell.landcoverfrac_old[CROPLAND]>0.0)	//(if !(gridcell.landcoverfrac[CROPLAND]>0.0))
		{
			gridcell.firstobj();
			while (gridcell.isobj) //Loop through stands:
			{
				Stand& stand=gridcell.getobj();
				if(stand.landcover==CROPLAND)
				{
					stand.firstobj();
					gridcell.killobj();
				}
				else
					gridcell.nextobj();
			}
		}
	}
//CFT stand dynamics end

//update C-pools for receiving stands:
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		Stand& stand=gridcell.getobj();
		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]>0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]>0.0)
		{
			double old_frac, added_frac, new_frac;
#ifdef multiple_natural_stands
			if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL)
#else
			if(stand.landcover!=CROPLAND)	
#endif					
			{
				old_frac=gridcell.landcoverfrac_old[stand.landcover];
				added_frac=landcoverfrac_change[stand.landcover];
				new_frac=gridcell.landcoverfrac[stand.landcover];
				stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
			}
#ifdef multiple_natural_stands
			else if(stand.landcover==NATURAL)
			{
				if(stand.first_year==date.year)
				{
					old_frac=0.0;
					added_frac=landcoverfrac_change[stand.landcover];
					new_frac=landcoverfrac_change[stand.landcover];
				}
				else
				{
					gridcell.nextobj();
					continue;
				}					
			}
#endif
			else if(stand.landcover==CROPLAND)
			{
				old_frac=gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.cftid];
				added_frac=cropstand_change[stand.cftid];
				new_frac=gridcell.landcoverfrac[CROPLAND]*gridcell.cftfrac[stand.cftid];
				stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
			}	

			if(LCchangeCtransfer)
			{
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
			}
		}
		gridcell.nextobj();
	}

#ifndef multiple_natural_stands
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		Stand& stand=gridcell.getobj();
		if(stand.landcover!=CROPLAND)
			stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
		else
			stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[stand.landcover]);

		gridcell.nextobj();
	}
#endif

	if(transfer_litter_leaf) delete[] transfer_litter_leaf;
	if(transfer_litter_wood) delete[] transfer_litter_wood;
	if(transfer_litter_root) delete[] transfer_litter_root;
	if(transfer_litter_repr) delete[] transfer_litter_repr;
	if(transfer_harvested_products_slow) delete[] transfer_harvested_products_slow;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////  End of Landcover stand dynamics and C-partitioning  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////  Sowing date algorithm. ////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void check_crop_temp_limits(Climate& climate, Gridcellpft& gridcellpft)
{
	int y,startyear;
	Pft& pft=gridcellpft.pft;

	//Check if spring and autumn conditions are present on this day:
		if(pft.ifsdspring && climate.temp>pft.tempspring && climate.dtemp_31[29]<=pft.tempspring)	//NB. after updating dtemp_31 with today's value !
		{//TeWW,TeCo,TeSf,TeRa: 12,14,13,12 (NB 5,14,15,5 in Bondeau publication);
			if (climate.lat>=0.0 && date.day>300)
				gridcellpft.last_springdate=date.day-365;
			else 
				gridcellpft.last_springdate=date.day;
			gridcellpft.springoccurred=true;
		}

		if(pft.ifsdautumn)
		{
			if(climate.temp<pft.tempautumn && climate.dtemp_31[29]>=pft.tempautumn && !gridcellpft.autumnoccurred)
			{//TeWW,TeRa: 12,17
				if (climate.lat>=0.0 && date.day<100)		
					gridcellpft.first_autumndate=date.day+365;
				else
					gridcellpft.first_autumndate=date.day;
				gridcellpft.autumnoccurred=true;
			}

			if(climate.temp<pft.trg && climate.dtemp_31[29]>=pft.trg && !gridcellpft.vernstartoccurred)
			{//TeWW,TeRa: 12,12
				gridcellpft.vernstartoccurred=true;
			}

			if(climate.temp>pft.trg && climate.dtemp_31[29]<=pft.trg)
			{//TeWW,TeRa: 12,12
				if (climate.lat>=0.0 && date.day>300)
					gridcellpft.last_verndate=date.day-365;	
				else 
					gridcellpft.last_verndate=date.day;
				gridcellpft.vernendoccurred=true;
			}
		}
}

void calc_crop_dates_20y_mean(Climate& climate, Gridcellpft& gridcellpft)
{
	int y,startyear;
	Pft& pft=gridcellpft.pft;

//Add past sowing season's dates to 20-year mean: TEST DAY !
//June 30(180) in the north, December 31(364) in the south ; calc_crop_dates_20y_mean is called from crop_sowing_gridcell() on testday_temp

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Check if spring and frost conditions occurred during the past year. If not, set this year's date to either sdate_default or climate.coldestday:

	// if no spring occured during last year 
	if (pft.ifsdspring && !gridcellpft.springoccurred)	//TeWW,TeCo,TeSf,TeRa
	{
		if(climate.temp<=pft.tempspring)		
		{
			gridcellpft.last_springdate=date.day;	
		}
		else
		{
			gridcellpft.last_springdate=climate.coldestday;
		}
	}

	// if no autumn occured during last year
	if(pft.ifsdautumn && !gridcellpft.autumnoccurred)		//TeWW,TeRa
	{
		if(climate.maxtemp<pft.tempautumn || climate.maxtemp>=pft.tempautumn && climate.mtemp_min<pft.tempautumn)
		{
			gridcellpft.first_autumndate=date.day;		
		}
		else										//too warm
		{
			gridcellpft.first_autumndate=climate.coldestday;	// day 14

			if(climate.lat>=0.0 && gridcellpft.first_autumndate<180)
				gridcellpft.first_autumndate+=365;	
		}
///////////
	}

	if(pft.ifsdautumn && !gridcellpft.vernendoccurred)				//TeWW,TeRa		
	{
			gridcellpft.last_verndate=gridcellpft.last_springdate+60;			// Bondeau sets this to coldest day, but I don't want last_verndate20 to be earlier than last_springdate20 (60 days is the maximum number of vernalization days)
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Update spring and frost date 20-year arrays and calculate 20 years average means:

	//1) this year
	if(pft.ifsdspring)									//TeWW,TeCo,TeSf,TeRa
		gridcellpft.last_springdate20=gridcellpft.last_springdate;

	if(pft.ifsdautumn)									//TeWW,TeRa
	{
		gridcellpft.first_autumndate20=gridcellpft.first_autumndate;
		if(date.year==1 && climate.lat>=0.0)				//No autumn first half of first year, set value to same as for second year !
			gridcellpft.first_autumndate_20[19]=gridcellpft.first_autumndate;

		gridcellpft.last_verndate20=gridcellpft.last_verndate;
	}

	//2) starting year (1st of 20 or less)
	startyear=20-(int)min(19,date.year);

	//3) past 20 years or less

	for (y=startyear;y<20;y++) 
	{
		if(pft.ifsdspring)
		{
			gridcellpft.last_springdate_20[y-1]=gridcellpft.last_springdate_20[y];//TeWW,TeCo,TeSf,TeRa
			gridcellpft.last_springdate20+=gridcellpft.last_springdate_20[y];
		}
		if (pft.ifsdautumn)										//TeWW,TeRa
		{
			gridcellpft.first_autumndate_20[y-1]=gridcellpft.first_autumndate_20[y];
			gridcellpft.first_autumndate20+=gridcellpft.first_autumndate_20[y];

			gridcellpft.last_verndate_20[y-1]=gridcellpft.last_verndate_20[y];
			gridcellpft.last_verndate20+=gridcellpft.last_verndate_20[y];
		}
	}

	//4) 20 years average means:
	if(pft.ifsdspring)
	{
		gridcellpft.last_springdate20/=(int)min(20,date.year+1);//TeWW,TeCo,TeSf,TeRa
		if (gridcellpft.last_springdate20<0)
			gridcellpft.last_springdate20+=365;	
		gridcellpft.last_springdate_20[19]=gridcellpft.last_springdate;
	}

	if(pft.ifsdautumn)								//TeWW,TeRa
	{
		gridcellpft.first_autumndate20/=(int)min(20,date.year+1);
		if (gridcellpft.first_autumndate20>364)
			gridcellpft.first_autumndate20-=365;			
		gridcellpft.first_autumndate_20[19]=gridcellpft.first_autumndate;

		gridcellpft.last_verndate20/=(int)min(20,date.year+1);
		if (gridcellpft.last_verndate20<0)
			gridcellpft.last_verndate20+=365;
		gridcellpft.last_verndate_20[19]=gridcellpft.last_verndate;
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void set_sdatecalc_temp(Climate& climate, Gridcellpft& gridcellpft)
{
	Pft& pft=gridcellpft.pft;

	if(pft.ifsdautumn)							//TeWW,TeRa:
	{
		//Use first_autumndate20 as sowing date if first_autumndate20 is set (autumn conditions met during the past 20 years):
		if(!((gridcellpft.first_autumndate20==climate.testday_temp || gridcellpft.first_autumndate20==climate.coldestday) && gridcellpft.first_autumndate%365==gridcellpft.first_autumndate20))
		{
			gridcellpft.sdatecalc_temp=gridcellpft.first_autumndate20;
			gridcellpft.wintertype=true;
		}
		// decision d'avoir du ble de printemps
		else				//if (gridcellpft.first_autumndate20==climate.coldestday)
		{
			if(!((gridcellpft.last_springdate20==climate.testday_temp || gridcellpft.last_springdate20==climate.coldestday) && gridcellpft.last_springdate==gridcellpft.last_springdate20))
			{
				gridcellpft.sdatecalc_temp=gridcellpft.last_springdate20;
				gridcellpft.wintertype=false;
			}
			else	// If neither spring or autumn occurred for the past 20 years.
			{
				if(climate.maxtemp<pft.tempspring)	//Too cold to sow at all.
					gridcellpft.sdatecalc_temp=-1;
				else								//Too warm; avoid warmest period.
				{
					gridcellpft.sdatecalc_temp=climate.coldestday;
					gridcellpft.wintertype=true;					
				}
			}
		}

		// If autumn first_autumndate20 is earlier than hlimitdate, use last_springdate20 (winter is too long):
		if (climate.lat>=0.0 && gridcellpft.sdatecalc_temp<=pft.hlimitdatenh && gridcellpft.sdatecalc_temp>180 
			|| climate.lat<0.0 && gridcellpft.sdatecalc_temp<=pft.hlimitdatesh)
		{
			gridcellpft.sdatecalc_temp=gridcellpft.last_springdate20;	// use last_springdate20 disregarding earlier choices	
			gridcellpft.wintertype=false;
		}

/*
		if(forcesowingdates && pft.forcesowingdate && gridcellpft.sdate_force>0)
		{
			if((abs(gridcellpft.sdate_force-gridcellpft.first_autumndate20)<=abs(gridcellpft.sdate_force-gridcellpft.last_springdate20)))
				gridcellpft.wintertype=true;
			else
				gridcellpft.wintertype=false;
			gridcellpft.sdatecalc_temp=gridcellpft.sdate_force;
		}
*/

		// Climatic limits for TeWW growth:	
		if(!strncmp(pft.name,"TeWW", strlen("TeWW")) && climate.mtemp_min20>15.0)
		{
			gridcellpft.sdatecalc_temp=-1;
		}
	}
	else if(pft.ifsdspring)								//TeCo,TeSf
	{
		if(!strncmp(pft.name,"TeCo", strlen("TeCo")))
			gridcellpft.sdatecalc_temp=(int)(60.0/85.0*(gridcellpft.last_springdate20-climate.adjustlat)+29.5+climate.adjustlat);
		else if(!strncmp(pft.name,"TeSf", strlen("TeSf")))
			gridcellpft.sdatecalc_temp=gridcellpft.last_springdate20;
#if defined NEWSOWINGDATE
		else
			gridcellpft.sdatecalc_temp=gridcellpft.last_springdate20;
#endif
	}

	gridcellpft.springoccurred=false;	
	gridcellpft.vernstartoccurred=false;	
	gridcellpft.vernendoccurred=false;	
	gridcellpft.autumnoccurred=false;
}

void set_sdatecalc_prec(Climate& climate, Gridcellpft& gridcellpft)
{
	Pft& pft=gridcellpft.pft;

	if(pft.hydrology==IRRIGATED) 
	{
		gridcellpft.sdatecalc_prec=gridcellpft.sdate_default;
		gridcellpft.precoccurred=true;
	}
	else
	{
		if(!gridcellpft.precoccurred && (climate.SOAsia && climate.sprec_2[1]>=110.0 && climate.sprec_2[0]<110.0
										|| !climate.SOAsia && climate.sprec_2[1]>=40.0 && climate.sprec_2[0]<40.0))
		{
			gridcellpft.first_precdate=date.day;
			gridcellpft.sdatecalc_prec=gridcellpft.first_precdate;
			gridcellpft.precoccurred=true;
		}

		if(date.day==climate.testday_prec)	// December 31(364) north, June 30(180) in the south; just resets precoccurred and fcalc_prec.
		{
			gridcellpft.precoccurred=false;
			gridcellpft.sdatecalc_prec=-1;
		}
	}
}

 /**
 * Calculates the Julian start and end day of a month.
 * January 1st is set to 0.
 */
void monthdates(int& start, int& end, int month){
	int months[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	start=0;
	int m = 0;
	while (m<month){
		start+=months[m];
		m++;
	}
	end=start+months[m]-1;
}

void calc_sowing_windows(Gridcell& gridcell)	// called on climate.testday_temp
{
	Climate& climate=gridcell.climate;
	int sow_month=0;

	if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP)
	{
		double max = 0.0;
		double sum=0.0;

		for (int m=0;m<12;m++)
		{
			sum=0.0;
			for (int i=0;i<4;i++){
				int mm = m+i;
				if(mm>=12){
					mm-=12;
				}
				//Implement a check later to see if it makes any difference to to use the precipitation only
				if(true){
					if (gridcell.climate.mpet20[mm] > 0.0){
						sum+= gridcell.climate.mprec20[mm]/gridcell.climate.mpet20[mm];
					}
				}else{
					sum+= gridcell.climate.mprec20[mm];
				}
			}

			if (sum>max){
				max=sum;
				//Months are stored as, 0-11
				sow_month=m;
			}
		}
	}

	pftlist.firstobj();
	while(pftlist.isobj)
	{
		Pft& pft=pftlist.getobj();
		Gridcellpft& gridcellpft=gridcell.pft[pft.id];

		if(pft.phenology==CROPGREEN)
		{
#if defined IRRIGATED_USE_TEMP_SDATE
			if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology==IRRIGATED)
#else
			if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC)
#endif
			{
//Set sowing window around sdatecalc_temp
				gridcellpft.swindow[0]=stepfromdate(gridcellpft.sdatecalc_temp, -15);
				gridcellpft.swindow[1]=stepfromdate(gridcellpft.sdatecalc_temp, 15);

				if(!gridcellpft.wintertype && dayinperiod(gridcellpft.swindow[0], stepfromdate(climate.coldestday, -100), climate.coldestday))
				{
					gridcellpft.swindow[0]=climate.coldestday;
//					gridcellpft.swindow[0]=gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(gridcellpft.swindow[1], stepfromdate(climate.coldestday, -100), climate.coldestday))
						gridcellpft.swindow[1]=climate.coldestday;
				}

				if(gridcellpft.wintertype && dayinperiod(gridcellpft.swindow[1], climate.coldestday, stepfromdate(climate.coldestday, 100)))
				{
					gridcellpft.swindow[1]=climate.coldestday;
//					gridcellpft.swindow[1]=gridcellpft.sdatecalc_temp;	//gives better yields, but sdate transition not smooth
					if(dayinperiod(gridcellpft.swindow[0], climate.coldestday, stepfromdate(climate.coldestday, 100)))
						gridcellpft.swindow[0]=climate.coldestday;
				}

				if(gridcellpft.sdatecalc_temp==-1)
				{
					gridcellpft.swindow[0]=-1;
					gridcellpft.swindow[1]=-1;
				}
			}
#if defined IRRIGATED_USE_TEMP_SDATE
			else if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology!=IRRIGATED)
#else
			else if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP)
#endif
			{
#if defined FAO_PREC_SEASONS
				if(climate.prec_range==WET)	//Variability of prec, but prec/pet always above 1.0: same sowing as for SEASONALITY_NO
				{
					gridcellpft.swindow[0]=gridcellpft.sdate_default;
					gridcellpft.swindow[1]=stepfromdate(gridcellpft.sdate_default, 15);
				}
				else
#endif
				{
					monthdates(gridcellpft.swindow[0],gridcellpft.swindow[1],sow_month);

					//A conservative choice to expand the month
					//that is used for the search
					gridcellpft.swindow[0]=stepfromdate(gridcellpft.swindow[0], -15);
				}
			}
			else if(climate.seasonality==SEASONALITY_NO)
			{
				gridcellpft.swindow[0]=gridcellpft.sdate_default;
				gridcellpft.swindow[1]=stepfromdate(gridcellpft.sdate_default, 15);
			}

		// Climatic limits for TeWW growth:	
			if(!strncmp(pft.name,"TeWW", strlen("TeWW")) && climate.mtemp_min20>15.0)
			{
				gridcellpft.swindow[0]=-1;
				gridcellpft.swindow[1]=-1;
			}
		}
		pftlist.nextobj();
	}
}

void calc_m_climate_20y_mean(Climate& climate)	//called last day of the year from crop_sowing_gridcell()
{
	int i, m, y;
	int startyear=20-(int)min(19,date.year);
	double var_temp=0, var_prec=0;
	double mtemp20kelvin[12], prec_pet_ratio20[12];
	double mprec_petmin_thisyear=1.0;
	double mprec_petmax_thisyear=0.0;

	memset(mtemp20kelvin,0,12*sizeof(double));
	memset(prec_pet_ratio20,0,12*sizeof(double));

	for(m=0;m<12;m++)
	{
		//1) this year
		climate.mtemp20[m]=climate.mtemp_year[m];
		climate.mprec20[m]=climate.mprec_year[m];
		climate.mpet20[m]=climate.mpet_year[m];
		if(climate.mpet_year[m]>0.0)
			climate.mprec_pet20[m]=climate.mprec_year[m]/climate.mpet_year[m];
		else
			climate.mprec_pet20[m]=0.0;


		if(climate.mprec_year[m]/climate.mpet_year[m]<mprec_petmin_thisyear)
			mprec_petmin_thisyear=climate.mprec_year[m]/climate.mpet_year[m];
		if(climate.mprec_year[m]/climate.mpet_year[m]>mprec_petmax_thisyear)
			mprec_petmax_thisyear=climate.mprec_year[m]/climate.mpet_year[m];

		//3) past 20 years or less
		for (y=startyear;y<20;y++) 
		{
			climate.mtemp_20[y-1][m]=climate.mtemp_20[y][m];
			climate.mtemp20[m]+=climate.mtemp_20[y][m];

			climate.mprec_20[y-1][m]=climate.mprec_20[y][m];
			climate.mprec20[m]+=climate.mprec_20[y][m];

			climate.mpet_20[y-1][m]=climate.mpet_20[y][m];
			climate.mpet20[m]+=climate.mpet_20[y][m];

			climate.mprec_pet_20[y-1][m]=climate.mprec_pet_20[y][m];
			climate.mprec_pet20[m]+=climate.mprec_pet_20[y][m];
		}
		//4) 20 years average means:
		climate.mtemp20[m]/=min(20,date.year+1);
		climate.mprec20[m]/=min(20,date.year+1);
		climate.mpet20[m]/=min(20,date.year+1);
		climate.mprec_pet20[m]/=min(20,date.year+1);

		climate.mtemp_20[19][m]=climate.mtemp_year[m];
		climate.mprec_20[19][m]=climate.mprec_year[m];
		climate.mpet_20[19][m]=climate.mpet_year[m];
		if(climate.mpet_year[m]>0.0)
			climate.mprec_pet_20[19][m]=climate.mprec_year[m]/climate.mpet_year[m];
		else
			climate.mprec_pet_20[19][m]=0.0;
	}

	climate.mprec_petmin20=mprec_petmin_thisyear;
	climate.mprec_petmax20=mprec_petmax_thisyear;
	for (y=startyear;y<20;y++) 
	{
		climate.mprec_petmin_20[y-1]=climate.mprec_petmin_20[y];
		climate.mprec_petmin20+=climate.mprec_petmin_20[y];
		climate.mprec_petmax_20[y-1]=climate.mprec_petmax_20[y];
		climate.mprec_petmax20+=climate.mprec_petmax_20[y];
	}
	climate.mprec_petmin20/=min(20,date.year+1);
	climate.mprec_petmin_20[19]=mprec_petmin_thisyear;
	climate.mprec_petmax20/=min(20,date.year+1);
	climate.mprec_petmax_20[19]=mprec_petmax_thisyear;
}

static double variation_coefficient(double data[],int n)
{
	//0 and 1 will give division with zero.
	if(n>1){
		  double avg,dev=0,varcoe=0,sum=0;
		  int i;
		  double std=0;
		  for (i=0;i<n;i++){
			  sum+=data[i];
		  }
		  avg=fabs(sum/n);
		  for (i=0;i<n;i++){
			  dev+=(data[i]-avg)*(data[i]-avg);
		  }
		  std=sqrt(fabs(dev/(n-1)));
		  if (std>0 && avg>0){
			  varcoe=std/avg;//just if temperature and precipitation data appear in the cell
		  }
		  return varcoe;
	}else{
		return -1.0;
	}
}

void calc_seasonality(Gridcell& gridcell){	//called last day of the year from crop_sowing_gridcell()

	Climate& climate=gridcell.climate;
	double var_temp=0, var_prec=0;
	double TEMPMIN = 10.0; //minimum temperature of coldest month
	const int NMONTH = 12;												
	double mtempKelvin[NMONTH], prec_pet_ratio20[12];
	double maxprec_pet20=0.0;
	double minprec_pet20=1000;

	memset(mtempKelvin,0,NMONTH*sizeof(double));
	memset(prec_pet_ratio20,0,12*sizeof(double));


	 //The temperature has got to be in Kelvin,
	 //the limit 0.010 is based on that.

	for(int i=0;i<NMONTH;++i){
		mtempKelvin[i]=gridcell.climate.mtemp20[i]+273.15;
		prec_pet_ratio20[i]=(gridcell.climate.mpet20[i] > 0) ? gridcell.climate.mprec20[i]/gridcell.climate.mpet20[i] : 0; //calculate P/PET ratio if monthly PET is above zero
	}
	
	var_temp=variation_coefficient(mtempKelvin,NMONTH);
	var_prec=variation_coefficient(prec_pet_ratio20,NMONTH);

	gridcell.climate.var_prec = var_prec;
	gridcell.climate.var_temp = var_temp;

	if (var_prec<=0.4 && var_temp<=0.010)					//no seasonality
	{
		climate.seasonality=SEASONALITY_NO;					//0
	}
	else if (var_prec>0.4)
	{
		if(var_temp<=0.010)									//precipitation seasonality only
			climate.seasonality=SEASONALITY_PREC;			//1
		else if(var_temp>0.010)
		{
			if(gridcell.climate.mtemp_min20>TEMPMIN)		//both seasonalities, but "weak" temperature seasonality (coldest month > 10degC)
				climate.seasonality=SEASONALITY_PRECTEMP;	//2
			else if(gridcell.climate.mtemp_min20<TEMPMIN)	//both seasonalities, but temperature most important
				climate.seasonality=SEASONALITY_TEMPPREC;	//4
		}
	}
	else if(var_prec<=0.4)
	{ 
		if (var_temp>0.010)
		{													//Temperature seasonality only
			climate.seasonality=SEASONALITY_TEMP;			//3
		}
	}

	for(int m=0;m<12;m++)
	{
		if(climate.mprec_pet20[m]>maxprec_pet20)
			maxprec_pet20=climate.mprec_pet20[m];
		if(climate.mprec_pet20[m]<minprec_pet20)
			minprec_pet20=climate.mprec_pet20[m];
	}

	if(minprec_pet20<=0.5 && maxprec_pet20<=0.5)								//Extremes av monthly means.
		climate.prec_seasonality=DRY;				//0
	else if(minprec_pet20<=0.5 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)		
		climate.prec_seasonality=DRY_INTERMEDIATE;	//1
	else if(minprec_pet20<=0.5 && maxprec_pet20>1.0)
		climate.prec_seasonality=DRY_WET;			//2
	else if(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>0.5 && maxprec_pet20<=1.0)
		climate.prec_seasonality=INTERMEDIATE;		//3
	else if(minprec_pet20>1.0 && maxprec_pet20>1.0)
		climate.prec_seasonality=WET;				//5
	else if(minprec_pet20>0.5 && minprec_pet20<=1.0 && maxprec_pet20>1.0)		
		climate.prec_seasonality=INTERMEDIATE_WET;	//4
	else
		dprintf("Problem with calculating precipitation seasonality !\n");

	if(climate.mprec_petmin20<=0.5 && climate.mprec_petmax20<=0.5)				//Average of extremes.
		climate.prec_range=DRY;						//0
	else if(climate.mprec_petmin20<=0.5 && climate.mprec_petmax20>0.5 && climate.mprec_petmax20<=1.0)
		climate.prec_range=DRY_INTERMEDIATE;		//1
	else if(climate.mprec_petmin20<=0.5 && climate.mprec_petmax20>1.0)			
		climate.prec_range=DRY_WET;					//2
	else if(climate.mprec_petmin20>0.5 && climate.mprec_petmin20<=1.0 && climate.mprec_petmax20>0.5 && climate.mprec_petmax20<=1.0)
		climate.prec_range=INTERMEDIATE;			//3
	else if(climate.mprec_petmin20>1.0 && climate.mprec_petmax20>1.0)
		climate.prec_range=WET;						//5
	else if(climate.mprec_petmin20>0.5 && climate.mprec_petmin20<=1.0 && climate.mprec_petmax20>1.0)
		climate.prec_range=INTERMEDIATE_WET;		//4
	else
		dprintf("Problem with calculating precipitation range !\n");

	if(climate.mtemp_max20<=10)
		climate.temp_seasonality=COLD;				//0
	else if(climate.mtemp_min20<=10 && climate.mtemp_max20>10 && climate.mtemp_max20<=30)
		climate.temp_seasonality=COLD_WARM;			//1
	else if(climate.mtemp_min20<=10 && climate.mtemp_max20>30)
		climate.temp_seasonality=COLD_HOT;			//2
	else if(climate.mtemp_min20>10 && climate.mtemp_max20<=30)
		climate.temp_seasonality=WARM;				//3
	else if(climate.mtemp_min20>30)
		climate.temp_seasonality=HOT;				//5
	else if(climate.mtemp_min20>10 && climate.mtemp_max20>30)	
		climate.temp_seasonality=WARM_HOT;			//4
	else
		dprintf("Problem with calculating temperature seasonality !\n");
}

void crop_sowing_gridcell(Gridcell& gridcell)
{
	int d,y,startyear;
	Climate& climate=gridcell.climate;

	if (date.year==0 && date.day==0) 
	{
		for (d=0;d<10;d++)
			climate.dprec_10[d]=climate.prec;
		for (d=0;d<2;d++)
			climate.sprec_2[d]=climate.prec;
	}

	climate.sprec_2[0]=climate.sprec_2[1];
	climate.sprec_2[1]=climate.prec;
	for (d=0;d<9;d++) 
	{
		climate.dprec_10[d]=climate.dprec_10[d+1];
		climate.sprec_2[1]+=climate.dprec_10[d];
	}
	climate.dprec_10[9]=climate.prec;

	if(climate.temp>climate.maxtemp)	//To know if temperature rises over vernalization limit
		climate.maxtemp=climate.temp;

	// Loop through PFTs
	pftlist.firstobj();
	while (pftlist.isobj) 
	{
		Pft& pft=pftlist.getobj();
		Gridcellpft& gridcellpft=gridcell.pft[pft.id];

		if(pft.landcover==CROPLAND)
		{
			if(pft.ifsdcalc)		//TeWW,TrRi,TeCo,TrMi,TrMa,TeSf,TrPe,TeRa; sdate set in getgridcell() kept for the rest;  (no code here for TrRi)
			{
				if(pft.ifsdtemp)	//TeWW,TeCo,TeSf,TeRa
				{
					check_crop_temp_limits(climate, gridcellpft);

					//Add past sowing season's dates to 20-year mean: TEST DAY !
					if(date.day==climate.testday_temp)	// June 30(180) in the north, December 31(364) in the south
					{
						calc_crop_dates_20y_mean(climate, gridcellpft);

						// Determine sdatecalc_temp:						
						set_sdatecalc_temp(climate, gridcellpft);

						climate.maxtemp=climate.temp;	//Rese maxtemp to today's value
					}
				}

				// Determine sdatecalc_prec:
				if(pft.ifsdprec)	//TeCo,TrMi,TrMa,TrPe:
				{
					set_sdatecalc_prec(climate, gridcellpft);
				}
			}
		}
		// ... on to next PFT
		pftlist.nextobj();
	}

#if defined NEWSOWINGDATE
	if(date.day==climate.testday_temp)	//day 180/364
		calc_sowing_windows(gridcell);
#endif

	if(date.islastmonth && date.islastday)
	{
		calc_m_climate_20y_mean(climate);
		calc_seasonality(gridcell);
	}
}


void Crop_sowing_date_temp(Patch& patch, Pft& pft)	
{
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;
	Patchpft& patchpft=patch.pft[pft.id];
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];

	patchpft.set_cropphen()->sdate=gridcellpft.sdatecalc_temp;
}

void Crop_sowing_date_prec(Patch& patch, Pft& pft)
{
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;
	Patchpft& patchpft=patch.pft[pft.id];
	cropphen_struct& ppftcrop=*(patchpft.get_cropphen());
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];

	int first_sowdate;
	int last_sowdate;

	if(!strncmp(pft.name,"TeCo", strlen("TeCo")))	//Sådd mellan sdatecalc_temp och sdate_default (140)
	{
		first_sowdate=gridcellpft.sdatecalc_temp;
		last_sowdate=gridcellpft.sdate_default;

#if defined ALLOW_TWOSEASONSPREC
// allows for two growing seasons; may sow outside the window ! eg. sow 140 - harvest 360 - sow 361
		if(climate.lat<45.0)	//Precipitation only determines sdate at latitudes < 45.0
		{
			if(date.day==ppftcrop.sdate)
			{
				if(gridcellpft.precoccurred)
					return;							//use gridcellpft.sdatecalc_temp as sdate
				else
					patchpft.cropphen->sdate=-1;	//wait for rain period to begin
			}
			else if(dayinperiod(date.day,first_sowdate,last_sowdate) && gridcellpft.precoccurred)	//sow when rain period begins
				ppftcrop.sdate=date.day;
			else if(date.day==last_sowdate && !gridcellpft.precoccurred)							//if rainperiod has not begun at last_sowdate, sow anyway
				ppftcrop.sdate=date.day;
		}
		else
			return;									//use gridcellpft.sdatecalc_temp as sdate
#else
		//Sowing window starts here at sdatecalc_temp:
		if (date.day==ppftcrop.sdate && !gridcellpft.precoccurred && climate.lat<45.0)		//sdatecalc_temp: om regnperiod ej börjat, vänta till regn !
			ppftcrop.sdate=-1;																
		//Sowing occurred when rain has triggered sdatecalc_prec to be set:
		if (date.day==gridcellpft.sdatecalc_prec && ppftcrop.sdate==-1 && climate.lat<45.0)	//sdatecalc_temp ej satt (sdate sätts till -1 på hdate)
		{
			if(dayinperiod(date.day,first_sowdate,last_sowdate))
				ppftcrop.sdate=date.day;
		}
		//If sdatecalc_prec has not been set by rain at sdate_default, sow anyway:
		if(date.day==gridcellpft.sdate_default && !gridcellpft.precoccurred && climate.lat<45.0)	//har inte regnat här dag 140, så ändå ! 
			ppftcrop.sdate=date.day;
#endif

	}
	else			//TrMi,TrMa,TrPe:					//sådd måste ligga mellan firstsowdatenh/sh och dag 210/30
	{

		if(climate.lat>=0.0)
		{
			first_sowdate=pft.firstsowdatenh_prec;
			last_sowdate=210;
		}
		else
		{
			first_sowdate=pft.firstsowdatesh_prec;
			last_sowdate=30;
		}

#if defined ALLOW_TWOSEASONSPREC
//new code:	allows for two growing seasons
		if (date.day==first_sowdate)	//Har regnperioden börjat innan first_sowdate ? Så i så fall denna dag !
		{
			if(gridcellpft.precoccurred)
			{
				ppftcrop.sdate=first_sowdate;
			}
		}
		else if(dayinperiod(date.day,first_sowdate,last_sowdate) && gridcellpft.precoccurred)	//Two seasons if hdate < last_sowdate !
			ppftcrop.sdate=date.day;
		else if(date.day==last_sowdate && !gridcellpft.precoccurred)
			ppftcrop.sdate=date.day;
#else
		if (date.day==gridcellpft.sdatecalc_prec)	//firstsowdatenh_prec = sdatenh-(20 to 40)
		{
			if(date.day<=pft.firstsowdatenh_prec && climate.lat>=0.0 || date.day<=pft.firstsowdatesh_prec && date.day>180 && climate.lat<0.0)
			{	//if calculated sowing date is earlier than 20-40 days before the default sdate, use the latter
				ppftcrop.sdate=first_sowdate;
			}
			else if (date.day<=210 && climate.lat>=0.0 || (date.day<=30 || date.day>180) && climate.lat<0.0)
				ppftcrop.sdate=gridcellpft.sdatecalc_prec;
		}

		if (gridcellpft.sdatecalc_prec==-1 && (date.day==210 && climate.lat>=0.0 || date.day==30 && climate.lat<0.0))
			ppftcrop.sdate=date.day;
#endif
	}
}

void Crop_sowing_date_rice(Patch& patch, Pft& pft)
{
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;
	Patchpft& patchpft=patch.pft[pft.id];
	cropphen_struct& ppftcrop=*(patchpft.get_cropphen());
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];

	if(!gridcellpft.singlecrop)	// Rice with only one season has fixed sdate and hlimitdate.
	{
		if(date.day==ppftcrop.sdate)
		{
			if(ppftcrop.maincrop)
			{
				if(climate.lat>=0.0)
					ppftcrop.hlimitdate=pft.hlimitdatesh; //Achtung! on purpose...
				else
					ppftcrop.hlimitdate=pft.hlimitdatenh;
			}
			else
				ppftcrop.hlimitdate=gridcellpft.hlimitdate_default;
		}
		else if(date.day==ppftcrop.hdate+1 && ppftcrop.hdate>0)	
		{
			if(ppftcrop.maincrop)
			{
				ppftcrop.maincrop=false;
				ppftcrop.sdate=(ppftcrop.hdate+30)%365;	
			}
			else
			{
				ppftcrop.maincrop=true;
				ppftcrop.sdate=gridcellpft.sdate_default;
			}
		}
	}
}

void Crop_sowing_date_forced(Patch& patch, Pft& pft)
{
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;
	Patchpft& patchpft=patch.pft[pft.id];
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];

	if(pft.forcesowingdate && pft.phenology!=ANY && !(!strncmp(pft.name,"TrRi", strlen("TrRi")) && !gridcellpft.singlecrop))	//not applicable for rice with several growingseasons
	{
		if(gridcellpft.sdate_force>0)
		{
			patchpft.cropphen->sdate=gridcellpft.sdate_force;

			if (climate.lat>=0.0 && gridcellpft.sdate_force<=pft.hlimitdatenh && gridcellpft.sdate_force>180 
				|| climate.lat<0.0 && gridcellpft.sdate_force<=pft.hlimitdatesh)
			{
				patchpft.cropphen->hlimitdate=gridcellpft.sdate_force-1;
			}
		}

	}
	else if(pft.forcesowingdate && (!strncmp(pft.name,"TrRi", strlen("TrRi")) && !gridcellpft.singlecrop))
	{
		dprintf("No sowing date data available for rice with several growingseasons, using calculated values\n");
	}
}


void Crop_sowing_date(Patch& patch, Pft& pft)
{
	Patchpft& patchpft=patch.pft[pft.id];

	if(pft.ifsdcalc)
	{
		//inidiv.sdate och hlimitdate set to default gridcellpft value in establishment, i.e. pft.sdatenh/sdatesh och pft.hlimitdatenh/sh
		//ppftcrop.hlimitdate==gridcellpft.hlimitdate_default (==pft.hlimitdatenh/sh) except for TrRi, when it may vary if !singlecrop.
		if(pft.ifsdtemp)	//TeWW,TeCo,TeSf,TeRa
		{
			if(date.day==patchpft.get_cropphen()->hdate+1 && patchpft.get_cropphen()->hdate>0 || date.day==patchpft.get_cropphen()->hlimitdate+1 && patchpft.get_cropphen()->hlimitdate>0)
				Crop_sowing_date_temp(patch, pft);
		}

		if(pft.ifsdprec)	//TeCo,TrMi,TrMa,TrPe:
			Crop_sowing_date_prec(patch, pft);

		if(!strncmp(pft.name,"TrRi", strlen("TrRi")))
			Crop_sowing_date_rice(patch, pft);
	}
}

void Crop_sowing_date_new(Patch& patch, Pft& pft)	// Enters here every day when growingseason==false
{
	Patchpft& patchpft=patch.pft[pft.id];
	cropphen_struct& ppftcrop=*(patchpft.get_cropphen());
	Gridcell& gridcell=patch.stand.gridcell;
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];
	Climate& climate=gridcell.climate;
	double length_growseas_def;
	bool temp_sdate=false, prec_sdate=false, def_sdate=false;

#if defined IRRIGATED_USE_TEMP_SDATE
	if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology==IRRIGATED)
		temp_sdate=true;
	else if((climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology!=IRRIGATED) && climate.prec_range!=WET)
		prec_sdate=true;
#else
	if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC)
		temp_sdate=true;
	else if((climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP) && climate.prec_range!=WET)
		prec_sdate=true;
#endif
	else //if(climate.seasonality==SEASONALITY_NO) || (climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP) && climate.prec_range==WET)
		def_sdate=true;

#ifdef IRRIGATED_USE_TEMP_SDATE
	if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology==IRRIGATED)
#else
	if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC)
#endif
	{
		if(dayinperiod(patchpft.swindow[0], stepfromdate(ppftcrop.hlimitdate, -100), ppftcrop.hlimitdate))
		{
			patchpft.swindow[0]=ppftcrop.hlimitdate+1;
			if(dayinperiod(patchpft.swindow[1], stepfromdate(ppftcrop.hlimitdate, -100), ppftcrop.hlimitdate))
				patchpft.swindow[1]=ppftcrop.hlimitdate+1;
		}
	}

#ifndef SD_TEMP_WINDOW
	if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology==IRRIGATED)
	{
		if(date.day==ppftcrop.hdate+1 && ppftcrop.hdate>0 || date.day==ppftcrop.hlimitdate+1 && ppftcrop.hlimitdate>0)
			Crop_sowing_date_temp(patch, pft);
		return;
	}
#endif

	if(dayinperiod(date.day,patchpft.swindow[0],patchpft.swindow[1]))
	{
		if(date.day!=patchpft.swindow[1])
		{
#if defined IRRIGATED_USE_TEMP_SDATE
			if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology==IRRIGATED)
#else
			if(climate.seasonality==SEASONALITY_TEMP || climate.seasonality==SEASONALITY_TEMPPREC)
#endif
			{
#if defined SD_TEMP_WINDOW
				if(gridcellpft.wintertype && climate.temp<pft.tempautumn || !gridcellpft.wintertype && climate.temp>pft.tempspring)
				{
					ppftcrop.sdate=date.day;
					ppftcrop.hlimitdate=gridcellpft.hlimitdate_default;	//in case of a reversion from precipitation seasonality
				}
#endif
			}
#if defined IRRIGATED_USE_TEMP_SDATE
			else if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology!=IRRIGATED)
#else
			else if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP)
#endif
			{
#if defined FAO_PREC_SEASONS
				if(climate.prec_range==WET)	//Variability of prec, but prec/pet always above 1.0: same sowing as for SEASONALITY_NO; not applicable to Africa for historic period.
				{
					patchpft.cropphen->sdate=date.day;
					patchpft.cropphen->hlimitdate=gridcellpft.hlimitdate_default;	//in case of a reversion from precipitation seasonality; NB. rice set to sdate+200 in Crop_sowing_date_rice()-default is 304 days.
				}
				else
#endif
				if(climate.prec>0.1 || pft.hydrology==IRRIGATED)
				{
					ppftcrop.sdate=date.day;

					if(gridcellpft.sdate_default>gridcellpft.hlimitdate_default)
						length_growseas_def=gridcellpft.hlimitdate_default+365-gridcellpft.sdate_default;
					else
						length_growseas_def=gridcellpft.hlimitdate_default-gridcellpft.sdate_default;

					ppftcrop.hlimitdate=stepfromdate(date.day, length_growseas_def);
				}
			}
			else if(climate.seasonality==SEASONALITY_NO)
			{
				ppftcrop.sdate=date.day;	//sowing window starts with sdate_default
				ppftcrop.hlimitdate=gridcellpft.hlimitdate_default;	//in case of a reversion from precipitation seasonality; NB. rice set to sdate+200 in Crop_sowing_date_rice()-default is 304 days.
			}
		}
		else
		{
			ppftcrop.sdate=date.day;

#ifdef IRRIGATED_USE_TEMP_SDATE
			if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP && pft.hydrology!=IRRIGATED)
#else
			if(climate.seasonality==SEASONALITY_PREC || climate.seasonality==SEASONALITY_PRECTEMP)
#endif
			{
				if(gridcellpft.sdate_default>gridcellpft.hlimitdate_default)
					length_growseas_def=gridcellpft.hlimitdate_default+365-gridcellpft.sdate_default;
				else
					length_growseas_def=gridcellpft.hlimitdate_default-gridcellpft.sdate_default;

				ppftcrop.hlimitdate=stepfromdate(date.day, length_growseas_def);
			}
			else
				ppftcrop.hlimitdate=gridcellpft.hlimitdate_default;	//in case of a reversion from precipitation seasonality
		}
		if(date.day==ppftcrop.sdate)
		{
			if(gridcellpft.sdate_default>gridcellpft.hlimitdate_default)
				length_growseas_def=gridcellpft.hlimitdate_default+365-gridcellpft.sdate_default;
			else
				length_growseas_def=gridcellpft.hlimitdate_default-gridcellpft.sdate_default;

			length_growseas_def=min(length_growseas_def, 245.0);

			if(pft.ifsdautumn)
			{
				if(temp_sdate==true)
					ppftcrop.hucountend=stepfromdate(ppftcrop.hlimitdate, -20);
				else if(prec_sdate==true && climate.prec_seasonality<=DRY_WET)
				{
					if(pft.hydrology==IRRIGATED)
						ppftcrop.hucountend=stepfromdate(date.day, 230);
					else
						ppftcrop.hucountend=stepfromdate(date.day, 210);
				}
				else if(def_sdate==true)
					ppftcrop.hucountend=stepfromdate(date.day, 230);
			}
			else if(!strncmp(pft.name,"TrRi", strlen("TrRi")) && gridcellpft.singlecrop)
				ppftcrop.hucountend=stepfromdate(date.day, 230);
			else
			{
				if(prec_sdate && climate.prec_seasonality<=DRY_WET)
				{
					if(pft.hydrology==IRRIGATED)
						ppftcrop.hucountend=stepfromdate(date.day, length_growseas_def);
					else
						ppftcrop.hucountend=stepfromdate(date.day, min(length_growseas_def,210.0)); //Shorter growing period when risk for water stress.
				}
				else
					ppftcrop.hucountend=stepfromdate(date.day, length_growseas_def);
			}			
		}
	}

	if(!strncmp(pft.name,"TrRi", strlen("TrRi")))	
	{											//TeRiirr enters here every day when growingseason==false, but does nothing when singlecrop==true
		Crop_sowing_date_rice(patch, pft);		//Sets hlimitdate (alternative default) for maincrop and sdate (30 days after hdate) and hlimitdate (default) for second crop
	}
}

void crop_sowing_patch(Patch& patch)
{
	// Loop through PFTs
	pftlist.firstobj();
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;

	while(pftlist.isobj) 
	{
		Pft& pft=pftlist.getobj();
		Patchpft& patchpft=patch.pft[pft.id];
		Gridcellpft& gridcellpft=gridcell.pft[pft.id];

		if(patch.stand.pft[pft.id].active && pft.phenology==CROPGREEN)
		{
//			cropphen_struct& ppftcrop=*(patchpft.cropphen);
			cropphen_struct& ppftcrop=*(patchpft.get_cropphen());

			if(!ppftcrop.growingseason)
			{
#if defined NEWSOWINGDATE
				if(date.day==ppftcrop.hdate+1 && ppftcrop.hdate>0 || date.day==0)
				{
					patchpft.swindow[0]=gridcellpft.swindow[0];
					patchpft.swindow[1]=gridcellpft.swindow[1];
				}
 
				Crop_sowing_date_new(patch, pft);
#else			//Old sowing date method:
				Crop_sowing_date(patch, pft);
#endif
				// Use sowing date read from input file if (pft.forcesowingdate==true)
				if(forcesowingdates)
					Crop_sowing_date_forced(patch, pft);
			}

//Calculation of eicdate:
			if(ppftcrop.sdate!=-1)
			{
				ppftcrop.eicdate=ppftcrop.sdate-15;

				if(ppftcrop.eicdate<0)
					ppftcrop.eicdate=365+ppftcrop.sdate-15;

				if(date.day>ppftcrop.eicdate && date.day<=ppftcrop.sdate
					|| ppftcrop.sdate<ppftcrop.eicdate && date.day+365>ppftcrop.eicdate && date.day<=ppftcrop.sdate)
				{
					if(ppftcrop.intercropseason)
						ppftcrop.eicdate=date.day;
					else if(ppftcrop.bicdate==ppftcrop.sdate) 
						ppftcrop.bicdate=-1;
				}
			}
		}
		// ... on to next PFT
		pftlist.nextobj();
	}//while(pftlist.isobj) 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////  End of Sowing date algorithm. /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////  Crop phenology  /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double senescence_curve(Pft& pft, double fphu)
{
	double senfactor;

	if (pft.shapesenescencenorm)
		senfactor=pow(1-fphu,2)/pow(1-pft.fphusen,2)*(1-pft.flaimaxharvest)+pft.flaimaxharvest;
	else
		senfactor=pow(1-fphu,0.5)/pow(1-pft.fphusen,0.5)*(1-pft.flaimaxharvest)+pft.flaimaxharvest;

	return senfactor;
}

void phu_init(cropphen_struct& ppftcrop, Gridcellpft& gridcellpft, Patch& patch)
{
	Pft& pft=gridcellpft.pft;
	Climate& climate=patch.stand.gridcell.climate;

			ppftcrop.husum=0.0;	
			ppftcrop.vrf=1.0;
			ppftcrop.vdsum=0;
			ppftcrop.prf=1.0;
			ppftcrop.fphu=0.0;
			ppftcrop.fhi=0.0;
			ppftcrop.fhi_phen=0.0;
			ppftcrop.fhi_water=1.0;
			ppftcrop.hdate=-1;
			ppftcrop.bicdate=-1;	

			ppftcrop.growingseason=true;
			ppftcrop.nsow++;

			if(ppftcrop.nsow==1)
				ppftcrop.sdate_thisyear[0]=ppftcrop.sdate;	
			else if(ppftcrop.nsow==2)
				ppftcrop.sdate_thisyear[1]=ppftcrop.sdate;

			ppftcrop.pvd=pft.pvd;		//default; kept for TrMi, TePu, TeSb, TrMa, TeSo, TrPe
			ppftcrop.phu=pft.phu;
			ppftcrop.tb=pft.tb;

#if defined DYNAMIC_PHU
			int years=min(date.year-patch.stand.first_year-1, 9);

			if(patch.stand.first_year!=date.year)
				ppftcrop.husum_max_10=(ppftcrop.husum_max_10*years+ppftcrop.husum_max)/(years+1);
			ppftcrop.husum_max=0.0;
#endif

			if(pft.ifsdautumn)	//TeWW,TeRa
			{
				if(gridcellpft.wintertype)	// Autumn sowing
				{
					if((gridcellpft.first_autumndate20==climate.testday_temp || gridcellpft.first_autumndate20==climate.coldestday) && gridcellpft.first_autumndate%365==gridcellpft.first_autumndate20
						&& (gridcellpft.last_springdate20==climate.testday_temp || gridcellpft.last_springdate20==climate.coldestday) && gridcellpft.last_springdate==gridcellpft.last_springdate20)	//wintertype if neither spring or winter conditions for the past 20 years
					{
						ppftcrop.pvd=pft.pvd;
						ppftcrop.phu=pft.phu;
					}
					else if(!(gridcellpft.last_verndate20==gridcellpft.last_springdate20+60 && gridcellpft.last_verndate==gridcellpft.last_verndate20))	//not all past 20 years without vernendoccurred: too cold
					{
						//pvd:

						//NB: vernalization (below 12 degrees) is supposed to occur directly at sowing (trg=tempautumn) for TeWW, for TeRa a 20-day lag (tempautumn=17) ??
						if((ppftcrop.sdate<180 || gridcellpft.last_verndate20>=180) && climate.lat>=0.0 || climate.lat<0.0)	// first_autumndate20 occurred before last_verndate20
						{
							if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
							{
								ppftcrop.pvd=(int)min(60,gridcellpft.last_verndate20-ppftcrop.sdate);
							}
							else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
								ppftcrop.pvd=(int)max(0,min(60,gridcellpft.last_verndate20-ppftcrop.sdate-20));
						}
						else if(!strncmp(pft.name,"TeWW", strlen("TeWW")))													// first_autumndate20 occurred after last_verndate20
							ppftcrop.pvd=(int)min(60,gridcellpft.last_verndate20+365-ppftcrop.sdate);
						else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))													// first_autumndate20 occurred after last_verndate20
							ppftcrop.pvd=(int)max(0,min(60,gridcellpft.last_verndate20+365-ppftcrop.sdate-20));
						//phu:
						if (ppftcrop.sdate<184+climate.adjustlat) 
						{
							if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
								ppftcrop.phu=max(1700.0,-0.1081*pow((double)(ppftcrop.sdate-climate.adjustlat),2)+3.1633*((double)(ppftcrop.sdate-climate.adjustlat))+2876.9);
							else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
								ppftcrop.phu=max(2100.0,-0.1081*pow((double)(ppftcrop.sdate-climate.adjustlat),2)+3.1633*((double)(ppftcrop.sdate-climate.adjustlat))+3279.7);
						}
						else
						{
							if(!strncmp(pft.name,"TeWW", strlen("TeWW")))
							{

								ppftcrop.phu=max(1700.0,-0.1081*pow((double)ppftcrop.sdate-365,2)+3.1633*((double)ppftcrop.sdate-365)+2876.9);
								ppftcrop.phu*=0.8;
							}
							else if(!strncmp(pft.name,"TeRa", strlen("TeRa")))
								ppftcrop.phu=max(2100.0,-0.1081*pow((double)ppftcrop.sdate-365,2)+3.1633*((double)ppftcrop.sdate-365)+3279.7);
						}
					}
				}
				else if (!gridcellpft.wintertype)
				{
					//If last_verndate has occurred during the past 20 year (or too warm):
					if(!(gridcellpft.last_verndate20==ppftcrop.sdate+60 && gridcellpft.last_verndate==gridcellpft.last_verndate20))
					{
						ppftcrop.pvd=(int)min(60,gridcellpft.last_verndate20-ppftcrop.sdate);
						ppftcrop.phu=1300.0;
					}
////				If no last_verndate occurred during the past 20 year (too cold, what about too warm ?):
					else 
					{
						ppftcrop.pvd=30;
						ppftcrop.phu=1500.0;
					}
				}
			}
			else if(!strncmp(pft.name,"TrRi", strlen("TrRi")) || !strncmp(pft.name,"TeSf", strlen("TeSf")))
			{
				if(!strncmp(pft.name,"TeSf", strlen("TeSf")))
				{
					ppftcrop.phu=min(2000.0,max(1300.0,-700.0/90.0*(ppftcrop.sdate-climate.adjustlat)+2460.0));
				}
				if (!strncmp(pft.name,"TrRi", strlen("TrRi")) && date.year<=1)
				{
					if (patch.stand.gridcell.get_lon()<60.0 || patch.stand.gridcell.get_lat()>30.0)
						ppftcrop.phu=1600.0;
				}
			}

#if defined DYNAMIC_PHU
			ppftcrop.phu_old=ppftcrop.phu;							//phu_old mainly for printout

			if(patch.stand.first_year!=date.year)					//Insert condition here to use dynamic phu for a limited time
				ppftcrop.phu=max(900.0, 0.9*ppftcrop.husum_max_10);
#endif
}


void leaf_phenology_crop(Pft& pft, Patch& patch) 
{

	double hu=0.0,k,c;
	Gridcell& gridcell=patch.stand.gridcell;
	Climate& climate=gridcell.climate;
	Patchpft& patchpft=patch.pft[pft.id];
	Gridcellpft& gridcellpft=gridcell.pft[pft.id];

	cropphen_struct& ppftcrop=*(patchpft.get_cropphen());
	ppftcrop.growingseason_ystd=ppftcrop.growingseason;

	if(pft.phenology==CROPGREEN)	//excludes CC3G and CC4G
	{

		if(date.day==0)
		{
			ppftcrop.fphu_harv=-1.0;
			ppftcrop.fhi_harv=-1.0;
			ppftcrop.sdate_harv=-1;
			ppftcrop.nsow=0;
			ppftcrop.sendate=-1;
			ppftcrop.nharv=0;

			ppftcrop.sownlastyear=false;
			for(int i=0;i<2;i++)
			{
				ppftcrop.sdate_harvest[i]=-1;
				ppftcrop.hdate_harvest[i]=-1;
//				ppftcrop.fphu_harvest[i]=-1.0;
//				ppftcrop.fhi_harvest[i]=-1.0;
				ppftcrop.sdate_thisyear[i]=-1;
			}
		}

//Calculation of PVD, PHU and TB:
		if(date.day==ppftcrop.sdate)
		{
			phu_init(ppftcrop, gridcellpft, patch);
		}

//Calculation of FPHU, PHEN och HI
		// loop on days from sowing to maturity, calculates daily fraction of plant maximal LAI after sowing has taken place
		if (date.day==ppftcrop.sdate || ppftcrop.growingseason) 
		{
			ppftcrop.senescence_ystd=ppftcrop.senescence;
			ppftcrop.hi_ystd=ppftcrop.hi;
			ppftcrop.intercropseason=false;			

			// before maturity is reached

			if (!forceharvestdates && ppftcrop.husum<ppftcrop.phu && (date.day<ppftcrop.hlimitdate || ppftcrop.sdate>ppftcrop.hlimitdate && date.day>=ppftcrop.sdate) 
				|| (forceharvestdates && pft.forceharvestdate && gridcellpft.hdate_force>0 && (date.day<gridcellpft.hdate_force || ppftcrop.sdate>gridcellpft.hdate_force && date.day>=ppftcrop.sdate)) )
			{
// Uträkning av fphu:
#if defined MAXHUTEMP
				hu=max(0.0,min(climate.temp, 30.0)-ppftcrop.tb);
#else
				hu=max(0,climate.temp-ppftcrop.tb);
#endif

				// accounting for vernalization if needs for vernalization not yet satisfied	//trg=tb for crops other than TeWW and TeRa and don't enter here
				if (climate.temp<pft.trg && ppftcrop.vdsum<ppftcrop.pvd)						//trg=12 for TeWW & TeRa
				{
					ppftcrop.vdsum=ppftcrop.vdsum+1;
					ppftcrop.vrf=min(1.0,(double)ppftcrop.vdsum/(double)ppftcrop.pvd);

					// vernalisation reduction factor has no effect once temp>trg even if needs are not satisfied
					// no effect as well if temp>trg at the beginning of the growing season...
					hu=hu*ppftcrop.vrf;																		
				}

				// accounting for response to photoperiod
				ppftcrop.prf=(1-pft.psens)*min(1.0,max(0.0,(climate.daylength_save[date.day]-pft.pb)/(pft.ps-pft.pb)))+pft.psens;
				// Achtung difference daylength/photoperiod !
				hu=hu*ppftcrop.prf;																				

				// daily effective temperature sum (ï¿½Cd) (daily heat units accumulation assuming high temperatures do not influence leaf phenology, like in SWAT)
				ppftcrop.husum=ppftcrop.husum+hu;

				// phenological scale (fraction of growing season)
				ppftcrop.fphu=min(1.0,ppftcrop.husum/ppftcrop.phu);						//SWAT 5:2.1.11

				patchpft.phen=1.0;

				if (ppftcrop.fphu>=pft.fphusen) 
				{
					if(ppftcrop.senescence_ystd==false)
					{
						ppftcrop.sendate=date.day;
					}
					ppftcrop.senescence=true;
				}


// HARVEST INDEX CALCULATION:

				ppftcrop.hi=pft.hiopt*100*ppftcrop.fphu/(100*ppftcrop.fphu+exp(11.1-10.0*ppftcrop.fphu));	//SWAT 5:2.4.1
				ppftcrop.fhi_phen=ppftcrop.hi/pft.hiopt;

// Correction of HI according to water stress:

				double wdf;
				double fwdf;
				double hi_save;

				ppftcrop.demandsum_crop+=climate.eet*PRIESTLEY_TAYLOR;		//demamdsum_crop==petsum, supplysum_crop==aetsum
				if (patchpft.supply>patch.demand) 
					ppftcrop.supplysum_crop+=patch.demand; 
				else
					ppftcrop.supplysum_crop+=patchpft.supply;

				if(ppftcrop.demandsum_crop>0.0)							
					wdf=100.0*ppftcrop.supplysum_crop/ppftcrop.demandsum_crop;	//SWAT 5:3.3.2	: aetsum/petsum
				else
					wdf=100.0;

				fwdf=wdf/(wdf+exp(6.13-0.0883*wdf));		// (SWAT 5:3.3.1) 

				hi_save=ppftcrop.hi;

				ppftcrop.hi=ppftcrop.fhi_phen* ((pft.hiopt-pft.himin)*fwdf +pft.himin);											

				if(ppftcrop.hi>0.0)
					ppftcrop.fhi_water=ppftcrop.hi/hi_save;

				ppftcrop.fhi=ppftcrop.fhi_phen*ppftcrop.fhi_water;

			} //before maturity is reached
//	HARVEST ! //
			else	//simulates harvest (a grain drying phase can follow maturity before harvest)
			{
				// Save phenological values at harvest:
				ppftcrop.fphu_harv=ppftcrop.fphu;
				ppftcrop.fhi_harv=ppftcrop.fhi;

				patchpft.phen=0.0;
				ppftcrop.lai_crop_actual=0.0;
				ppftcrop.demandsum_crop=0.0;
				ppftcrop.supplysum_crop=0.0;
				ppftcrop.hdate=date.day;
				ppftcrop.sdate_harv=ppftcrop.sdate;

				ppftcrop.nharv++;

				if(ppftcrop.nharv==1)
				{
					ppftcrop.sdate_harvest[0]=ppftcrop.sdate;
					ppftcrop.hdate_harvest[0]=date.day;
//					ppftcrop.fphu_harvest[0]=ppftcrop.fphu;
//					ppftcrop.fhi_harvest[0]=ppftcrop.fhi;
					if(ppftcrop.sdate>date.day)							
						ppftcrop.sownlastyear=true;
				}
				else if(ppftcrop.nharv==2)
				{
					ppftcrop.sdate_harvest[1]=ppftcrop.sdate;
					ppftcrop.hdate_harvest[1]=date.day;
//					ppftcrop.fphu_harvest[1]=ppftcrop.fphu;
//					ppftcrop.fhi_harvest[1]=ppftcrop.fhi;
				}

				ppftcrop.bicdate=ppftcrop.hdate+15;

				if(ppftcrop.bicdate>365)
					ppftcrop.bicdate=ppftcrop.hdate-365+15;

				if(pft.ifsdprec)					//TeCo,TrMi,TrMa,TrPe; with NEWSOWINGDATE: all crops
				{
					ppftcrop.sdate=-1;
					ppftcrop.eicdate=-1;
				}

				ppftcrop.growingseason=false;
				ppftcrop.intercropseason=false;		
				ppftcrop.senescence=false;
			} //end harvest

			ppftcrop.lai=ppftcrop.lai_crop_actual;
			ppftcrop.fpc=1.0-lambertbeer(ppftcrop.lai);
		}  //from sowing has taken place until harvest day

#if defined DYNAMIC_PHU	// Every day
		if(ppftcrop.growingseason==false && dayinperiod(date.day, ppftcrop.hdate, ppftcrop.hucountend) && ppftcrop.hdate>=0)
		{
#if defined MAXHUTEMP
				hu=max(0.0,min(climate.temp, 30.0)-ppftcrop.tb);
#else
				hu=max(0,climate.temp-ppftcrop.tb);
#endif

			// accounting for vernalization if needs for vernalization not yet satisfied
			if (ppftcrop.vdsum<ppftcrop.pvd && climate.temp<pft.trg)						//trg=12 for TeWW & TeRa 
			{
				ppftcrop.vdsum=ppftcrop.vdsum+1;
				ppftcrop.vrf=min(1.0,(double)ppftcrop.vdsum/(double)ppftcrop.pvd);

				// vernalisation reduction factor has no effect once temp>trg even if needs are not satisfied
				// no effect as well if temp>trg at the beginning of the growing season...
				hu=hu*ppftcrop.vrf;																		
			}

			// accounting for response to photoperiod
			ppftcrop.prf=(1-pft.psens)*min(1.0,max(0.0,(climate.daylength_save[date.day]-pft.pb)/(pft.ps-pft.pb)))+pft.psens;
			// Achtung difference daylength/photoperiod !
			hu=hu*ppftcrop.prf;	

			if(date.day==ppftcrop.hdate)
				ppftcrop.husum_max_postharv=0.0;
			ppftcrop.husum_max_postharv+=hu;		//Local hu

			if(date.day==ppftcrop.hucountend)
			{
				ppftcrop.husum_max_postharv-=hu;						//Don't count the hu's on harvest day !
				ppftcrop.husum_max=ppftcrop.husum_max_postharv+ppftcrop.husum;
				ppftcrop.husum_max_hlim=ppftcrop.husum_max;
				ppftcrop.husum_h=ppftcrop.husum;
			}
		}
#endif
		if(pft.intercrop==NATURALGRASS && !ppftcrop.growingseason)
		{
			if(!ppftcrop.intercropseason && date.day==ppftcrop.bicdate)
			{
				ppftcrop.intercropseason=true;
			}

			if(ppftcrop.intercropseason && date.day==ppftcrop.eicdate)
			{
				ppftcrop.intercropseason=false;
				ppftcrop.demandsum_crop=0.0;
				ppftcrop.supplysum_crop=0.0;
			}
		}
	}// if (pft.phenology==CROPGREEN)
	else if(pft.phenology==ANY) // crop grasses using standard guess phenology calculation
	{
		if(patch.stand.pftid==pft.id 
			|| patch.stand.hasgrassintercrop && (patch.pft[patch.stand.pftid].cropphen->intercropseason 
			|| date.day==patch.pft[patch.stand.pftid].cropphen->bicdate || date.day==patch.pft[patch.stand.pftid].cropphen->eicdate))//normal CC3G & CC4G (+ irrigated) growth
		{																	
			if(patch.stand.pftid!=pft.id)
			{									
				if(date.day==patch.pft[patch.stand.pftid].cropphen->bicdate)
				{
					ppftcrop.growingseason=true;
					patch.stand.gdd0_intercrop=climate.gdd5_pasture;
				}
				if(date.day==patch.pft[patch.stand.pftid].cropphen->eicdate)
				{
					ppftcrop.growingseason=false;
					patch.stand.gdd0_intercrop=0.0;
					patchpft.phen=0.0;
				}
			}

			if (climate.lat>=0.0 && date.day==COLDEST_DAY_NHEMISPHERE || climate.lat<0.0 && date.day==COLDEST_DAY_SHEMISPHERE)
				patch.stand.gdd0_intercrop=0.0;

			if(ppftcrop.growingseason)	// NB on bicdate, not on eicdate
			{
				if(patch.stand.pftid==pft.id)	//Normal grass growth: gives identical result to natural stands.
					patchpft.phen=min(1.0,climate.gdd5_pasture/pft.phengdd5ramp);
				else if(patch.stand.gdd0_intercrop>0.0)
					patchpft.phen=min(1.0,(climate.gdd5_pasture-patch.stand.gdd0_intercrop)/(pft.phengdd5ramp*0.9));
				else
					patchpft.phen=min(1.0,(climate.gdd5_pasture-patch.stand.gdd0_intercrop)/pft.phengdd5ramp);

				if(patchpft.phen<0.0)
					patchpft.phen=0.0;

				if (patchpft.wscal<pft.wscal_min)
				{
					patchpft.phen=0.0;
					patch.stand.gdd0_intercrop=climate.gdd5_pasture;
				}
			}
		}
	}
}

void fpar_crop(Patch& patch) {

	double plai_grass; // summed LAI for grasses
	double plai_leafon_grass;
		// summed LAI for grasses assuming full leaf cover for all individuals
	double flai=0.0; // fraction of total grass LAI represented by a particular grass
	double fpar_grass; // FPAR at top of grass canopy
	double fpar_ff; // FPAR at forest floor (beneath grass canopy)
	double fpar_leafon_ff;
		// FPAR at forest floor assuming full leaf cover for all individuals
	double fpar_min; // minimum FPAR required for grass growth


	// Obtain reference to Vegetation object
	Vegetation& vegetation=patch.vegetation;

	// And to Climate object
	Climate& climate=patch.stand.gridcell.climate;

	if (vegmode==POPULATION) 
	{
		// POPULATION MODE

		// Loop through individuals

		vegetation.firstobj();
		while (vegetation.isobj) 
		{
			Individual& indiv=vegetation.getobj();
		
			// For this individual ...

			indiv.fpar=indiv.fpc*indiv.phen; // Eqn 1
			indiv.fpar_leafon=indiv.fpc; // Eqn 2

			vegetation.nextobj(); // ... on to next individual
		}
	}
	else 
	{
		// INDIVIDUAL OR COHORT MODE

		// Initialise individual FPAR, find maximum height of vegetation, calculate
		// individual LAI given current phenology, calculate summed LAI for grasses

		plai_grass=0.0;
		plai_leafon_grass=0.0;

		double highest_grass_lai=0.0;

		vegetation.firstobj();
		while (vegetation.isobj) 
		{
			Individual& indiv=vegetation.getobj();

			// For this individual ...
			indiv.fpar=0.0;
			indiv.fpar_leafon=0.0;

			if(patch.pft[indiv.pft.id].cropphen->growingseason==true)
			{
				if (indiv.pft.lifeform==GRASS)
				{
					if(indiv.lai>highest_grass_lai)
						highest_grass_lai=indiv.lai;
					plai_leafon_grass=highest_grass_lai;	// avoids double lai count for intercrop grass (c3 and c4 grass competing, lai is for monocultures)
					plai_grass+=indiv.lai*indiv.phen;	
				}
			}
			vegetation.nextobj(); // ... on to next individual
		}

		// FPAR reaching grass canopy
		fpar_grass=1.0;


		// Add grass LAI to calculate PAR reaching forest floor
		// BLARP: Order changed Ben 050301 to overcome optimisation bug in pgCC

		fpar_ff=lambertbeer(plai_grass);

		// Save this
		patch.fpar_ff=fpar_ff;	//patch.fpar_ff not used further

		fpar_leafon_ff=lambertbeer(plai_leafon_grass);


		// FPAR for grass PFTs is difference between relative PAR at top of grass canopy
		// canopy and at forest floor, or lower if FPAR at forest floor below threshold
		// for grass growth. PAR reaching the grass canopy is partitioned among grasses
		// in proportion to their LAI (a somewhat simplified assumption)

		// Loop through individuals

		vegetation.firstobj();
		while (vegetation.isobj) 
		{
			Individual& indiv=vegetation.getobj();
			Patchpft& patchpft=patch.pft[indiv.pft.id];

			// For this individual ...

			if (indiv.pft.lifeform==GRASS && patch.pft[indiv.pft.id].cropphen->growingseason==true)
			{
				// Calculate minimum FPAR for growth of this grass

				// Fraction of total grass LAI represented by this grass

				if (!negligible(climate.par))
					fpar_min=min(indiv.pft.parff_min/climate.par,1.0);	//set parff_min to 0 for crops !
				else
					fpar_min=1.0;

				if(indiv.pft.phenology==CROPGREEN)
					indiv.fpar=1-lambertbeer(max(0.0,indiv.phen*indiv.lai));	//phen is 1.0 during growingseason here
				else
				{
					if(indiv.cropindiv->isintercropgrass)	//may contain both c3 and c4 grass
					{
						if (!negligible(plai_grass))
							flai=indiv.lai*indiv.phen/plai_grass;
						else
							flai=1.0;
					}
					else	//monoculture
						flai=1.0;

					indiv.fpar=max(0.0,flai-max(fpar_ff*flai,fpar_min));
				}


				// Repeat assuming full leaf cover for all individuals

				if(indiv.pft.phenology==CROPGREEN)
					indiv.fpar_leafon=1-lambertbeer(max(0.0,indiv.lai));
				else
				{
					if(indiv.cropindiv->isintercropgrass)	//may contain both c3 and c4 grass
					{
						if (!negligible(plai_leafon_grass))
							flai=indiv.lai/plai_leafon_grass;
						else
							flai=1.0;					
					}
					else	//monoculture
						flai=1.0;
						
					indiv.fpar_leafon=max(0.0,flai-max(fpar_leafon_ff*flai,fpar_min));
				}
			}

			vegetation.nextobj();
		}

		// Save grass canopy FPAR
		patch.fpar_grass=fpar_grass;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////  End of crop phenology  ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////  Crop allocation  ///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void growth_crop_daily(Patch& patch)
{
/////////  Daily growth routine for crops /ML

	double froot,fleaf,fho;
	double grs_cmass_plant_old;
	double grs_cmass_root_old;
	double grs_cmass_leaf_old;
	double grs_cmass_ho_old;
	double grs_cmass_ag;

	Vegetation& vegetation=patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj)
	{
		Individual& indiv=vegetation.getobj();
		cropindiv_struct& cropindiv=*(indiv.get_cropindiv());
		Patchpft& patchpft=patch.pft[indiv.pft.id];
		cropphen_struct& ppftcrop=*(patchpft.get_cropphen());

		if(date.day==0)
		{
			cropindiv.ycmass_plant=0.0;
			cropindiv.ycmass_leaf=0.0;
			cropindiv.ycmass_root=0.0;
			cropindiv.ycmass_ho=0.0;
			cropindiv.ycmass_agpool=0.0;	

			cropindiv.harv_cmass_plant=0.0;
			cropindiv.harv_cmass_root=0.0;
			cropindiv.harv_cmass_leaf=0.0;
			cropindiv.harv_cmass_ho=0.0;
			cropindiv.harv_cmass_agpool=0.0;

			cropindiv.cmass_ho_harvest[0]=0.0;
			cropindiv.cmass_ho_harvest[1]=0.0;

			cropindiv.cmass_leaf_max=0.0;

			if(indiv.pft.phenology==ANY && !cropindiv.isintercropgrass)				//zero of normal cc3g/cc4g-growth arbitrarily at new year
			{
				cropindiv.grs_cmass_plant=0.0;
				cropindiv.grs_cmass_root=0.0;
				cropindiv.grs_cmass_ho=0.0;
				cropindiv.grs_cmass_leaf=0.0;
				cropindiv.grs_cmass_agpool=0.0;
			}
		}
/////////////////////////////////////////////////////////  TRUE CROP ALLOCATION  ///////////////////////////////////////////////////////////////
		if(indiv.pft.phenology==CROPGREEN)
		{
			if(ppftcrop.growingseason || date.day==ppftcrop.hdate)
			{

#define CMASS_SEED 0.01	// 10g/m2;
#ifdef DELAYED_SEEDCARBON
				if(dayinperiod(date.day, pppftcrop.sdate, (patchpft.cropphen->sdate+9))%365 )	// Seed carbon 110310; portion the seed carbon over a 10-day period.
				{
					cropindiv.grs_cmass_plant+=0.1*CMASS_SEED;
					cropindiv.ycmass_plant+=0.1*CMASS_SEED;
				}

#else
				if(date.day==ppftcrop.sdate)	// Seed carbon
				{
					cropindiv.grs_cmass_plant+=CMASS_SEED;
					cropindiv.ycmass_plant+=CMASS_SEED;

					patch.fluxes.acflux_seed-=CMASS_SEED;	// This flux will be balancing litter fluxes for the NEXT year, but the amount should be OK.
				}
#endif
				cropindiv.dcmass_plant=indiv.dnpp;
				cropindiv.grs_cmass_plant+=indiv.dnpp;
				cropindiv.ycmass_plant+=indiv.dnpp;

/////////// ROOT GROWTH																								
				froot=indiv.pft.frootstart-(indiv.pft.frootstart-indiv.pft.frootend)*ppftcrop.fphu;				//SWAT 5:2,1,21	
				grs_cmass_root_old=cropindiv.grs_cmass_root;
				cropindiv.grs_cmass_root=froot*cropindiv.grs_cmass_plant;
				cropindiv.dcmass_root=cropindiv.grs_cmass_root-grs_cmass_root_old;
				cropindiv.ycmass_root+=cropindiv.dcmass_root;

/////////// STORAGE ORGANS GROWTH
				grs_cmass_ag=(1.0-froot)*cropindiv.grs_cmass_plant;
				fho=ppftcrop.hi*(1.0-froot);			//SWAT
				grs_cmass_ho_old=cropindiv.grs_cmass_ho;

				if(indiv.pft.hiopt<=1.0)
				{
					cropindiv.grs_cmass_ho=ppftcrop.hi*grs_cmass_ag;								//SWAT 5:2.4.2, 5:2.4.4	(used by Bondeau)
				}
				else							//harvested roots (SWAT)
				{
					cropindiv.grs_cmass_ho=(1.0-1.0/(1.0+ppftcrop.hi))*cropindiv.grs_cmass_plant;	//SWAT 5:2.4.3 8 
				}
				cropindiv.dcmass_ho=cropindiv.grs_cmass_ho-grs_cmass_ho_old;	
				cropindiv.ycmass_ho+=cropindiv.dcmass_ho;	

/////////// LEAF GROWTH																
				grs_cmass_leaf_old=cropindiv.grs_cmass_leaf;	
				cropindiv.grs_cmass_leaf=cropindiv.grs_cmass_plant-cropindiv.grs_cmass_root-cropindiv.grs_cmass_ho;

				if(cropindiv.grs_cmass_plant>0.0)
					fleaf=cropindiv.grs_cmass_leaf/cropindiv.grs_cmass_plant;
				cropindiv.dcmass_leaf=cropindiv.grs_cmass_leaf-grs_cmass_leaf_old;
				cropindiv.ycmass_leaf+=cropindiv.dcmass_leaf;

//////////// ABOVE GROUND POOL (Bondeau)
				cropindiv.dcmass_agpool=cropindiv.dcmass_plant-cropindiv.dcmass_root-cropindiv.dcmass_leaf-cropindiv.dcmass_ho;		
				cropindiv.grs_cmass_agpool=cropindiv.grs_cmass_plant-cropindiv.grs_cmass_root-cropindiv.grs_cmass_leaf-cropindiv.grs_cmass_ho;
				cropindiv.ycmass_agpool=cropindiv.ycmass_plant-cropindiv.ycmass_root-cropindiv.ycmass_leaf-cropindiv.ycmass_ho;

				if(cropindiv.grs_cmass_leaf>cropindiv.cmass_leaf_max)	
				{
					cropindiv.cmass_leaf_max=cropindiv.grs_cmass_leaf;
				}

				if(date.day==ppftcrop.sendate)
				{
					cropindiv.cmass_leaf_sen=cropindiv.grs_cmass_leaf;
				}
				else if(date.day==ppftcrop.hdate)
				{
					cropindiv.dcmass_plant=0.0;
					cropindiv.dcmass_root=0.0;
					cropindiv.dcmass_ho=0.0;
					cropindiv.dcmass_leaf=0.0;
					cropindiv.dcmass_agpool=0.0;
				}

				if(ppftcrop.growingseason)
				{
					if(!ppftcrop.senescence)
					{
						ppftcrop.lai_crop_actual=cropindiv.grs_cmass_leaf*indiv.pft.sla;	
					}
					else
					{

						//Follow the Potsdam senescence curve from leaf cmass at senescence (cmass_leaf_sen):
						ppftcrop.lai_crop_actual=cropindiv.cmass_leaf_sen*indiv.pft.sla*senescence_curve(indiv.pft, ppftcrop.fphu);
					}

					if(ppftcrop.lai_crop_actual<0.0)
						ppftcrop.lai_crop_actual=0.0;
				}
				if(!(ppftcrop.lai_crop_actual>=0.0 && ppftcrop.lai_crop_actual<=20.0))//Test for unrealistically high lai.
if(!SUPPRESSLARGEOUTPUT)	
					dprintf("In growth_crop_daily() stand %d pft %d year %d day %d: senescence=%d, grs_cmass_leaf=%f, grs_cmass_ho=%f, lai_crop_actual=%f, out of bounds !\n", patch.stand.id, indiv.pft.id, date.year-nyear_spinup+1901, date.day, ppftcrop.senescence, cropindiv.grs_cmass_leaf, cropindiv.grs_cmass_ho, ppftcrop.lai_crop_actual);
			}

			// CROPGREEN COMMON TASKS AT HARVEST DAY
			if(date.day==ppftcrop.hdate)
			{
				cropindiv.harv_cmass_plant+=cropindiv.grs_cmass_plant;
				cropindiv.harv_cmass_root+=cropindiv.grs_cmass_root;
				cropindiv.harv_cmass_ho+=cropindiv.grs_cmass_ho;
				cropindiv.harv_cmass_leaf+=cropindiv.grs_cmass_leaf;
				cropindiv.harv_cmass_agpool+=cropindiv.grs_cmass_agpool;

				if(ppftcrop.nharv==1)
					cropindiv.cmass_ho_harvest[0]=cropindiv.grs_cmass_ho;
				else if(ppftcrop.nharv==2)
					cropindiv.cmass_ho_harvest[1]=cropindiv.grs_cmass_ho;

				cropindiv.grs_cmass_plant=0.0;
				cropindiv.grs_cmass_root=0.0;
				cropindiv.grs_cmass_ho=0.0;
				cropindiv.grs_cmass_leaf=0.0;
				cropindiv.grs_cmass_agpool=0.0;
				cropindiv.cmass_leaf_sen=0.0;	
			}
		}//END CROPGREEN
/////////////////////////////////////////////////////////  CROP GRASS ALLOCATION  ///////////////////////////////////////////////////////////////
		//NB: Only intercrop grass enters here ! normal cc3g/cc4g grass is treated just like natural grass except for functions called from growth() (harvest_crop() & allocation_crop()) 
		else if(indiv.pft.phenology==ANY && indiv.pft.id!=patch.stand.pftid)
		{
			if(ppftcrop.growingseason || ppftcrop.growingseason_ystd)	
			{																										
				indiv.ltor=indiv.wscal_mean*indiv.pft.ltor_max;	
										
				cropindiv.dcmass_plant=indiv.dnpp;
				cropindiv.grs_cmass_plant+=indiv.dnpp;		
				cropindiv.ycmass_plant+=indiv.dnpp;

/////////// ROOT GROWTH
				froot=1.0/(1.0+indiv.ltor);
				grs_cmass_root_old=cropindiv.grs_cmass_root;	

				//Cumulative wscal-dependent root increase						
				cropindiv.grs_cmass_root=froot*cropindiv.grs_cmass_plant;
				cropindiv.dcmass_root=cropindiv.grs_cmass_root-grs_cmass_root_old;
				cropindiv.ycmass_root+=cropindiv.dcmass_root;	// alt 2

/////////// LEAF GROWTH
				fleaf=1.0-froot;
				grs_cmass_leaf_old=cropindiv.grs_cmass_leaf;
				cropindiv.grs_cmass_leaf=cropindiv.grs_cmass_plant-cropindiv.grs_cmass_root;
				cropindiv.dcmass_leaf=cropindiv.grs_cmass_leaf-grs_cmass_leaf_old;
				cropindiv.ycmass_leaf+=cropindiv.dcmass_leaf;

/////////// HARVESTABLE ORGANS GROWTH
				if(date.day==patch.pft[patch.stand.pftid].get_cropphen()->eicdate)
				{
					cropindiv.dcmass_plant=0.0;
					cropindiv.dcmass_root=0.0;
					cropindiv.dcmass_ho=0.0;
					cropindiv.dcmass_leaf=0.0;
					cropindiv.dcmass_agpool=0.0;
				}
			}

			// "ANY" COMMON TASKS AT EICDATE
			if(date.day==patch.pft[patch.stand.pftid].get_cropphen()->eicdate)
			{
				cropindiv.harv_cmass_plant+=cropindiv.grs_cmass_plant;	
				cropindiv.harv_cmass_root+=cropindiv.grs_cmass_root;	
				cropindiv.harv_cmass_leaf+=cropindiv.grs_cmass_leaf;	
				cropindiv.harv_cmass_ho+=cropindiv.grs_cmass_ho;		
				cropindiv.harv_cmass_agpool+=cropindiv.grs_cmass_agpool;

				cropindiv.grs_cmass_plant=0.0;
				cropindiv.grs_cmass_root=0.0;
				cropindiv.grs_cmass_ho=0.0;
				cropindiv.grs_cmass_leaf=0.0;
				cropindiv.grs_cmass_agpool=0.0;

			}
		}// END ANY
		vegetation.nextobj();
	}

	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////  End of crop allocation  //////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////  Landcover harvest functions  ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
	{	
		harvest=indiv.pft.harv_eff*(cmass_sap+cmass_heart-cmass_debt);		//harvested products

		if(ifslowharvestpool)
		{
			harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;	//harvested products not consumed (oxidized) this year put into patchpft.harvested_products_slow
			harvest=harvest*(1-indiv.pft.harvest_slow_frac);
		}

		acflux_harvest+=harvest;							//harvested products consumed (oxidized) this year put into patch.fluxes.acflux_harvest, not litter pool !

		cmass_sap=(1-indiv.pft.harv_eff)*cmass_sap;			//unharvested parts of the plant
		cmass_heart=(1-indiv.pft.harv_eff)*cmass_heart;
		cmass_debt=(1-indiv.pft.harv_eff)*cmass_debt;		
	}

	residue_outtake=indiv.pft.res_outtake*(cmass_sap+cmass_heart-cmass_debt+cmass_leaf);
	acflux_harvest+=residue_outtake;																//removed residues

	litter_leaf+=cmass_leaf*(1-indiv.pft.res_outtake);												//not removed residues
	litter_wood+=(cmass_sap+cmass_heart-cmass_debt)*(1-indiv.pft.res_outtake);						//not removed residues

	cmass_sap=cmass_heart=cmass_debt=cmass_leaf=0.0;
}

void harvest_pasture(double& cmass_leaf,double& cmass_root,double& litter_leaf,double& litter_root,double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{
	double turnover, residue_outtake, harvest;
	double scale=1.0;
	int m;
	bool alive=indiv.alive;

//pasture version of turnover adapted from harvest_crop 110525 ; presently just represents grass being harvested (for use with GRASSFORCROP); only yearly harvest !
//sparar inte skörden i någon variabel än !
//kan fås genom cmass_leaf*indiv.pft.harv_eff*2.0
//!!!!!!!!!!!!!!!!!!!!!!!!!!! NB last year's C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//scale harvest products of stands with increased area by (old area/new area) if landcover change has occurred:	
	Stand& stand=indiv.vegetation.patch.stand;
	Gridcell& gridcell=stand.gridcell;

	if(gridcell.LC_updated)
	{
		scale=gridcell.landcoverfrac_old[PASTURE]/gridcell.landcoverfrac[PASTURE];	

		if(scale>=1.0)
			scale=1.0;
	}
///////////////////////////

//turnover and harvest of last year's carbon:

//NB. cmass_x can be negative here only if individuals with negative cmass-x are not killed last year.

	// Root turnover
	//Bondeau: turnover_root=0.5
	cmass_root*=scale;	
	cmass_leaf*=scale;	

	turnover=indiv.pft.turnover_root*cmass_root;	//turnover_root är normalt 0.7 för gräs
	if(alive && turnover>0.0) 
		litter_root+=turnover;
	cmass_root-=turnover;

	//OBS ! skörd före turnover !!!!
	//Harvest/Grazing:					
	//Bondeau: harv_eff=0.9 i kod, 0.5 i artikel plus 0.05 till litter (faeces)
	harvest=indiv.pft.harv_eff*cmass_leaf;						//använd 0.5 (this year's yield is set in allocation_crop)

	if(ifslowharvestpool)
	{
		harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;
		harvest=harvest*(1-indiv.pft.harvest_slow_frac);
	}
	acflux_harvest+=harvest;										//skördat gräs
	cmass_leaf-=harvest;

	// Leaf turnover
	turnover=indiv.pft.turnover_leaf*cmass_leaf;	//turnover_leaf är normalt 1.0 för gräs
	if(alive && turnover>0.0) 
		litter_leaf+=turnover;
	cmass_leaf-=turnover;

}

void harvest_crop(double& cmass_plant,double& cmass_leaf,double& cmass_root,double& cmass_ho,double& cmass_agpool,double& litter_leaf,double& litter_root,
				  double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{	//NB. this function is for balancing carbon fluxes based on last year's cmass, not for calculating this year's yield. This is done in allocation_crop().
	double turnover, residue_outtake, harvest;
	double scale=1.0;	
	int m;
	bool alive=indiv.alive;

//!!!!!!!!!!!!!!!!!!!!!!!!!!! NB last year's C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//scale harvest products of stands with increased area by (old area/new area) if landcover change has occurred:
	Stand& stand=indiv.vegetation.patch.stand;
	Gridcell& gridcell=stand.gridcell;

	if(gridcell.LC_updated)	
	{
		scale=gridcell.cftfrac_old[stand.cftid]*gridcell.landcoverfrac_old[CROPLAND]/(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
		if(scale>1.0)
			scale=1.0;
	}

	cmass_root*=scale;	
	cmass_leaf*=scale;	
	cmass_agpool*=scale;
	cmass_ho*=scale;	

//turnover and harvest (and acflux_harvest) of last year's carbon :

//NB. cmass_x can be negative here only if individuals with negative cmass-x are not killed last year.
	if(indiv.pft.phenology==CROPGREEN)
	{		
		if(alive && cmass_root>0.0)
			litter_root+=cmass_root;

		cmass_root=0.0;

		//Bondeau: harv_eff=1.0		
		if(alive && cmass_ho>0.0)										// (this year's yield is set in allocation_crop)
		{
			harvest=indiv.pft.harv_eff*cmass_ho;			//skördade produkter	

			if(indiv.pft.aboveground_ho)
				litter_leaf+=(cmass_ho-harvest);			// ej skördade produkter
			else
				litter_root+=(cmass_ho-harvest);			

			if(ifslowharvestpool)	
			{
				harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;	//patchpft.harvested_products_slow
				harvest=harvest*(1-indiv.pft.harvest_slow_frac);
			}
			acflux_harvest+=harvest;											//patch.fluxes.acflux_harvest
		}
		cmass_ho=0.0;

		//Bondeau: res_outtake=0.9 or 0.0
		if (alive && (cmass_leaf+cmass_agpool)>0.0)
		{
			residue_outtake=indiv.pft.res_outtake*(cmass_leaf+cmass_agpool);
			litter_leaf+=cmass_leaf+cmass_agpool-residue_outtake;						//ej uttagna rester

			acflux_harvest+=residue_outtake;											//uttagna rester
		}
		cmass_leaf=0.0;
		cmass_agpool=0.0;
		cmass_plant=0.0;

		//No turnover (no remaining live plant tissue after harvest) for real crops.
	}
	else if(indiv.pft.phenology==ANY)
	{
		if(indiv.cropindiv->isintercropgrass)			//Intercrop growth
		{
			if(alive && cmass_root>0.0)
				litter_root+=cmass_root;

			cmass_root=0.0;

			//Harvest/Grazing:			
			//Bondeau: harv_eff=0.9 i kod, 0 i artikel
			if(alive && cmass_leaf>0.0)
			{
				harvest=indiv.pft.harv_eff_ic*cmass_leaf;
				litter_leaf+=cmass_leaf-harvest;											//ej skördat gräs

				if(ifslowharvestpool)	
				{
					harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;
					harvest=harvest*(1-indiv.pft.harvest_slow_frac);
				}

				acflux_harvest+=harvest;										//skördat gräs	(inget för närvarande)
			}
			cmass_leaf=0.0;
			cmass_ho=0.0;														//cmass_ho används ej för gräs
			cmass_agpool=0.0;													//cmass_agpool används ej för gräs
			cmass_plant=0.0;

		}
		else								//Normal CC3G/CC4G stand growth (ej kollat om cmass>0.0 behövs än)
		{
			// Root turnover
			//Bondeau: turnover_root=0.5
			turnover=indiv.pft.turnover_root*cmass_root;	//turnover_root är normalt 0.7 för gräs

			if(alive && turnover>0.0) 
				litter_root+=turnover;

			cmass_root-=turnover;

			//OBS ! skörd före turnover !!!!
			//Harvest/Grazing:					
			//Bondeau: harv_eff=0.9 i kod, 0.5 i artikel plus 0.05 till flux (faeces)

			harvest=indiv.pft.harv_eff*cmass_leaf;						//använd 0.5;	(this year's yield is set in allocation_crop)
			cmass_leaf-=harvest;

			if(ifslowharvestpool)
			{
				harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;
				harvest=harvest*(1-indiv.pft.harvest_slow_frac);
			}

			acflux_harvest+=harvest;										//skördat gräs

			// Leaf turnover
			turnover=indiv.pft.turnover_leaf*cmass_leaf;	//turnover_leaf är normalt 1.0 för gräs

			if(alive && turnover>0.0)
				litter_leaf+=turnover;

			cmass_leaf-=turnover;
			cmass_plant=cmass_leaf+cmass_root;
		}
	}
}

void allocation_crop(double bminc,double cmass_leaf,double cmass_root,double cmass_ho,double ltor,
	double& cmass_plant_inc,double& cmass_leaf_inc,double& cmass_root_inc,double& cmass_ho_inc,double& cmass_agpool_inc, double& litter_leaf_inc,double& litter_root_inc, Individual& indiv)
{

	double froot, fho;
	double hi;

	cropindiv_struct& cropindiv=*(indiv.get_cropindiv());

	litter_leaf_inc=0.0;
	litter_root_inc=0.0;


	if (ltor<1.0e-10)	// normal cc3g/cc4g-growth: yearly mean, actual crops and intercrops: growingseason mean
	{
		// No leaf production possible - put all biomass into roots
		// (Individual will die next time period)

		cmass_leaf_inc=0.0;
		if(indiv.pft.phenology==ANY && indiv.pft.id==indiv.vegetation.patch.stand.pftid)	
			cmass_root_inc=bminc;
		else
			cmass_root_inc=cropindiv.ycmass_plant;	
		cmass_ho_inc=0.0;		
		cmass_agpool_inc=0.0;	

		return;
	}

	if(indiv.pft.phenology==ANY && indiv.pft.id==indiv.vegetation.patch.stand.pftid)			//Normal CC3G/CC4G stand growth
	{
		cmass_plant_inc=bminc;	//reproduction and cmass_excess reducement; daily and yearly values will NOT be compatible !

		cmass_leaf_inc=(cmass_plant_inc-cmass_leaf/ltor+cmass_root)/(1.0+1.0/ltor);
		cmass_root_inc=cmass_plant_inc-cmass_leaf_inc;

		if(cmass_leaf_inc<0.0) 
		{
			// Negative allocation to leaves
			cmass_root_inc=cmass_plant_inc;
			cmass_leaf_inc=(cmass_root+cmass_root_inc)*ltor-cmass_leaf; // Eqn (3)

			// Add killed leaves to litter
			litter_leaf_inc=-cmass_leaf_inc;
		}
		else if(cmass_root_inc<0.0) 
		{
			// Negative allocation to roots
			cmass_leaf_inc=cmass_plant_inc;
			cmass_root_inc=(cmass_leaf+cmass_plant_inc)/ltor-cmass_root;

			// Add killed roots to litter
			litter_root_inc=-cmass_root_inc;
		}

		if(cmass_leaf_inc>0.0)									
			cropindiv.yield=cmass_leaf_inc*indiv.pft.harv_eff*2.0;	// OK om turnover_leaf=1.0, annars (cmass_leaf+cmass_leaf_inc)*indiv.pft.harv_eff*2.0 !
		else
			cropindiv.yield=0.0;
		cropindiv.harv_yield=cropindiv.yield;							//Although no specified harvest date, harv_yield is set for compatibility.
	}
	else	// true crop growth and grass intercrop growth; NB: bminit (cmass_repr & cmass_excess subtracted) not used !
	{
		cmass_plant_inc=cropindiv.ycmass_plant;
		cmass_leaf_inc=cropindiv.ycmass_leaf;
		cmass_root_inc=cropindiv.ycmass_root;
		cmass_ho_inc=cropindiv.ycmass_ho;
		cmass_agpool_inc=cropindiv.ycmass_agpool;

		if(indiv.pft.phenology==ANY)				// grass intercrop growth
		{
			if(cropindiv.ycmass_leaf>0.0)									
				cropindiv.yield=cropindiv.ycmass_leaf*indiv.pft.harv_eff_ic*2.0;	//Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
			else
				cropindiv.yield=0.0;

			if(cropindiv.harv_cmass_leaf>0.0)			
				cropindiv.harv_yield=cropindiv.harv_cmass_leaf*indiv.pft.harv_eff_ic*2.0;	//Yield dry wieght of actually harvest products this year; NB as above
			else
				cropindiv.harv_yield=0.0;
		}
		else if(indiv.pft.phenology==CROPGREEN)		//true crop growth
		{
			if(cropindiv.ycmass_ho>0.0)									
				cropindiv.yield=cropindiv.ycmass_ho*indiv.pft.harv_eff*2.0;	//Yield dry wieght of allocated harvestable organs this year; NB independent from harvest calculation in harvest_crop (different years)
			else
				cropindiv.yield=0.0;

			if(cropindiv.harv_cmass_ho>0.0)									
				cropindiv.harv_yield=cropindiv.harv_cmass_ho*indiv.pft.harv_eff*2.0;	//Yield dry wieght of actually harvest products this year; NB as above
			else
				cropindiv.harv_yield=0.0;

			for(int i=0;i<2;i++)
			{
				if(cropindiv.cmass_ho_harvest[i]>0.0)								
					cropindiv.yield_harvest[i]=cropindiv.cmass_ho_harvest[i]*indiv.pft.harv_eff*2.0;	//Yield dry wieght of actually harvest products this year; NB as above
				else
					cropindiv.yield_harvest[i]=0.0;	
			}
		}
	}

	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////  End of landcover harvest functions  ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////