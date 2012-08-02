///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.cpp
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "landcover.h"
#include "guessio.h"

void landcover_init(Gridcell& gridcell,Pftlist& pftlist) {
	landcovertype landcover;

	getlandcover(gridcell,pftlist);		//Gets gridcell.landcoverfrac from landcover input file(s) or ins-file.

	for(int i=0;i<NLANDCOVERTYPES;i++) { //For all landcover types without subclasses
		if(i!=CROPLAND) {					// cropland subclasses turned off in this version
			if(run[i]) {
				if(gridcell.landcoverfrac[i]>0.0) {
					landcover=(landcovertype)i;
					Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
					stand.set_gridcell_fraction(gridcell.landcoverfrac[i]);
#ifdef MATS_TEST
					dprintf("Stand %d, landcover type %d created year %d. Initial fraction = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac[i]);
#endif
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
						Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
						stand.pftid=pft.id;
						stand.cftid=pft.cftid;
						stand.set_gridcell_fraction(gridcell.cftfrac[pft.cftid]*gridcell.landcoverfrac[CROPLAND]);

						stand.pft[pft.id].active=true;	//101213

						if(pft.hydrology==IRRIGATED)
						{
							stand.isirrigated=true;
#ifdef MATS_TEST			
							dprintf("Irrigated crop stand created, pft=%s\n", (char*)pft.name);
#endif
						}
						else
#ifdef MATS_TEST
						dprintf("Rainfed crop stand created, pft=%s\n", (char*)pft.name);
#endif
						if(pft.intercrop==NATURALGRASS && ifintercropgrass)	// && ifintercropgrass 101213
						{
							stand.hasgrassintercrop=true;

							for(int i=0;i<pftlist.nobj;i++)		//101213
							{
								if(pftlist[i].isintercropgrass)
									stand.pft[pftlist[i].id].active=true;
							}

#ifdef MATS_TEST
							dprintf("Crop stand with intercrop growth created, pft=%s\n", (char*)pft.name);
#endif
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

void landcover_dynamics(Gridcell& gridcell,Pftlist& pftlist)
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

//Save old fraction values:										//Flytta till getlandcover?
	for(i=0;i<NLANDCOVERTYPES;i++)
		gridcell.landcoverfrac_old[i]=gridcell.landcoverfrac[i];
	for(i=0;i<NCROPSTANDS_MAX;i++)
		cropfrac_sum_old+=gridcell.cftfrac_old[i]=gridcell.cftfrac[i];

//Get new gridcell.landcoverfrac and/or gridcell.cftfrac from LUdata and CFTdata.
	if(!all_fracs_const)					//ta bort?
		getlandcover(gridcell,pftlist);		//ta bort?
	else return;							//ta bort?

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
			if(i!=CROPLAND)											//Landcovers with only one stand.
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
#ifdef MATS_TEST
		dprintf("\nYear %d: changeLC=%f\tchangeCFT=%f\n", date.year-nyear_spinup+1901, changeLC, change_crop);	
		dprintf("Year %d: transferred_fraction=%f\n", date.year-nyear_spinup+1901, transferred_fraction);
		dprintf("Year %d: receiving_fraction=%f\n", date.year-nyear_spinup+1901, receiving_fraction);
		dprintf("Year %d: change_stand=%f\n", date.year-nyear_spinup+1901, change_stand);
#endif
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

#if defined cropLUchangeCtransfer

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
#ifdef MATS_TEST
					dprintf("Stand n:o %d, type %d coverage decreased year %d ! Difference = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);	
					dprintf("Stand n:o %d, type %d updated fraction year %d: %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);	
#endif
				}
//Landcovers that may have several stands:
#ifdef multiple_natural_stands
				else if(stand.landcover==NATURAL)
				{
					if(nnaturalstands>1)
					{
						// alt. 1: remove equal percentage from all natural stands:
//						scale=-stand.frac/gridcell.landcoverfrac[stand.landcover]*landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
//						stand.frac+=stand.frac/gridcell.landcoverfrac[stand.landcover]*landcoverfrac_change[stand.landcover];
						// alt. 2: remove equal amounts from all natural stands:
//						scale=-landcoverfrac_change[stand.landcover]/nnaturalstands/receiving_fraction/(double)stand.nobj;
//						stand.frac+=landcoverfrac_change[stand.landcover]/nnaturalstands;
						// alt.3:remove from stand according to receiving stand landcovertype:
						scale=-stand.natural_frac_change/receiving_fraction/(double)stand.nobj;

#ifdef MATS_TEST
						dprintf("Natural stand n:o %d coverage decreased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, stand.natural_frac_change);
						dprintf("Natural stand n:o %d updated absolute fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, stand.get_gridcell_fraction());
#endif
					}
					else
					{
						scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
						stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);

#ifdef MATS_TEST
						dprintf("Natural stand n:o %d coverage decreased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
						dprintf("Natural stand n:o %d updated absolute fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);
#endif
					}
				}
#endif			
				else if(stand.landcover==CROPLAND)
				{
					scale=-cropstand_change[stand.cftid]/receiving_fraction/(double)stand.nobj;
					stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
#ifdef MATS_TEST
					dprintf("Crop stand n:o %d coverage decreased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, cropstand_change[stand.cftid]);
					dprintf("Crop stand n:o %d updated absolute fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
#endif
				}	
	
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

						if(indiv.pft.landcover==CROPLAND)
						{
							cmass_ho_cp=indiv.get_cropinfo()->cmass_ho;
							cmass_agpool_cp=indiv.get_cropinfo()->cmass_agpool;
							cmass_plant_cp=indiv.get_cropinfo()->cmass_plant;
						}
	
						litter_leaf_cp=patchpft.litter_leaf;
						litter_root_cp=patchpft.litter_root;
						litter_wood_cp=patchpft.litter_wood;
						litter_repr_cp=patchpft.litter_repr;

						acflux_harvest_cp=patch.fluxes.acflux_harvest;	//flux är nollställd
						harvested_products_slow_cp=patch.pft[indiv.pft.id].harvested_products_slow;

	//Harvest of transferred areas:
						if(indiv.pft.landcover==CROPLAND)
							harvest_crop(cmass_plant_cp,cmass_leaf_cp,cmass_root_cp,cmass_ho_cp,cmass_agpool_cp,
							litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);
/*						else if(indiv.pft.landcover==PASTURE)
						{
//This code needs testing before use:
							harvest_pasture(cmass_leaf_cp,cmass_root_cp,
								litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp, indiv);
						}
*/						else												
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
#ifdef MATS_TEST
						dprintf("Stand %d pft %d: transfer_acflux_harvest=%f\n",stand.id, indiv.pft.id, transfer_acflux_harvest);
#endif
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
#endif

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
						Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
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
#ifdef MATS_TEST
						dprintf("Stand %d, landcover type %d created year %d. Initial fraction = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac[i]);
#endif
					}
					else if(gridcell.landcoverfrac_old[i]>0.0 && gridcell.landcoverfrac[i]==0.0)
					{
						gridcell.firstobj();
						while (gridcell.isobj) //Loop through stands:
						{
							Stand& stand=gridcell.getobj();
							if(stand.landcover==i)
							{
#ifdef MATS_TEST		
								dprintf("Stand %d, landcover type %d killed year %d. Lost fraction = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac_old[i]);	//flyttat hit 090928
#endif
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
						Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
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
#ifdef MATS_TEST
						dprintf("NEW natural stand no. %d created year %d. Initial fraction = %f\n", stand.id, date.year-nyear_spinup+1901, stand.get_gridcell_fraction());
#endif
					}
					else if(i==NATURAL && nnaturalstands>1 && landcoverfrac_change[i]<0.0)	// Natural stand killed.
					{
						gridcell.firstobj();
						while (gridcell.isobj)
						{
							Stand& stand=gridcell.getobj();
							if(stand.landcover==i && stand.get_gridcell_fraction()==0)
							{
#ifdef MATS_TEST		
								dprintf("Natural stand %d killed year %d. Lost fraction = %f\n", stand.id, date.year-nyear_spinup+1901, -(stand.natural_frac_change));
#endif
								gridcell.killobj();
							}
							else	//NB: else necessary after killobj() !
								gridcell.nextobj();
						}
					}
#endif
				}
			}
		}
	}

//CFT stand dynamics (from updated cftfrac):
	if(run[CROPLAND] && (change_crop>0.0 || landcoverfrac_change[CROPLAND]!=0.0))	// Fix 110323
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
#ifdef MATS_TEST
							dprintf("Crop stand %d fraction updated year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.cftfrac[pft.cftid]);
#endif
						}
						else
						{
							landcover=CROPLAND;
							Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
#ifdef MATS_TEST
							dprintf("Crop stand %d created year %d. Initial fraction = %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.cftfrac[pft.cftid]);
							dprintf("In landcover_dynamics() 1: stand %d year %d day %d: soil.wcont[0]=%f\n", stand.id, date.year-nyear_spinup+1901, date.day, stand[0].soil.wcont[0]);
#endif
							stand.pftid=pft.id;
							stand.cftid=pft.cftid;
							stand.set_gridcell_fraction(gridcell.cftfrac[pft.cftid]*gridcell.landcoverfrac[CROPLAND]);

							stand.pft[pft.id].active=true;

							if(pft.hydrology==IRRIGATED)
							{
								stand.isirrigated=true;
#ifdef MATS_TEST
								dprintf("Irrigated crop stand created, pft=%s\n", (char*)pft.name);
#endif
							}
							else
#ifdef MATS_TEST
								dprintf("Rainfed crop stand created, pft=%s\n", (char*)pft.name);
#endif
							if(pft.intercrop==NATURALGRASS && ifintercropgrass)	// && ifintercropgrass 101213
							{
								stand.hasgrassintercrop=true;

								for(int i=0;i<pftlist.nobj;i++)	
								{
									if(pftlist[i].isintercropgrass)
										stand.pft[pftlist[i].id].active=true;
								}

#ifdef MATS_TEST
								dprintf("Crop stand with intercrop growth created, pft=%s\n", (char*)pft.name);
#endif
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
#ifdef MATS_TEST
							dprintf("Crop stand %d killed year %d. Lost fraction = %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.pft[stand.pftid].pft.cftid]);	//090928 lagt till //100414
#endif
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
#ifdef MATS_TEST
					dprintf("Crop stand %d killed year %d. Lost fraction = %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.pft[stand.pftid].pft.cftid]);	//090928 flyttat hit //100414
#endif
					gridcell.killobj();
				}
				else	//NB: else necessary after killobj() !
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
#ifdef MATS_TEST
				dprintf("Stand n:o %d, type %d coverage increased year %d ! Difference = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
				dprintf("Stand n:o %d, type %d updated fraction year %d: %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);
				dprintf("Stand n:o %d, type %d fraction of receiving area year %d: %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]/receiving_fraction);
#endif
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
#ifdef MATS_TEST	
				dprintf("Natural stand n:o %d coverage increased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
				dprintf("Natural stand n:o %d updated fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, stand.get_gridcell_fraction());
				dprintf("Natural stand n:o %d fraction of receiving area year %d: %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]/receiving_fraction);
#endif
			}
#endif
			else if(stand.landcover==CROPLAND)
			{
				old_frac=gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.cftid];
				added_frac=cropstand_change[stand.cftid];
				new_frac=gridcell.landcoverfrac[CROPLAND]*gridcell.cftfrac[stand.cftid];
				stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
			}	
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
#ifdef MATS_TEST
				dprintf("In landcover_dynamics() 2: stand %d year %d day %d: soil.wcont[0]=%f, soil.temp=%f\n", stand.id, date.year-nyear_spinup+1901, date.day, stand[0].soil.wcont[0], stand[0].soil.temp);
				dprintf("transfer_wcont[0]=%f,old_frac=%f, added_frac=%f, new_frac=%f\n", transfer_wcont[0], old_frac, added_frac, new_frac);
#endif

//other soil stuff:
				for(i=0;i<NSOILLAYER;i++)
					patch.soil.wcont[i]=(patch.soil.wcont[i]*old_frac+transfer_wcont[i]*added_frac)/new_frac;
#ifdef MATS_TEST
				dprintf("In landcover_dynamics() 3: stand %d year %d day %d: soil.wcont[0]=%f, soil.temp=%f\n\n", stand.id, date.year-nyear_spinup+1901, date.day, stand[0].soil.wcont[0], stand[0].soil.temp);
#endif
				patch.soil.decomp_litter_mean=(patch.soil.decomp_litter_mean*old_frac+transfer_decomp_litter_mean*added_frac)/new_frac;
				patch.soil.k_soilfast_mean=(patch.soil.k_soilfast_mean*old_frac+transfer_k_soilfast_mean*added_frac)/new_frac;
				patch.soil.k_soilslow_mean=(patch.soil.k_soilslow_mean*old_frac+transfer_k_soilslow_mean*added_frac)/new_frac;
//add fluxes:
#ifdef MATS_TEST
				dprintf("Uppdatera C-pooler: Stand %d pft %d: transfer_acflux_harvest=%f\n",stand.id, stand.pftid, transfer_acflux_harvest);
				dprintf("Uppdatera C-pooler före: Stand %d pft %d: patch.fluxes.acflux_harvest=%f\n",stand.id, stand.pftid, patch.fluxes.acflux_harvest);
#endif
				patch.fluxes.acflux_harvest=(patch.fluxes.acflux_harvest*old_frac+transfer_acflux_harvest*added_frac)/new_frac;
#ifdef MATS_TEST
				dprintf("Uppdatera C-pooler efter: Stand %d pft %d: patch.fluxes.acflux_harvest=%f\n",stand.id, stand.pftid, patch.fluxes.acflux_harvest);
#endif
				stand.nextobj();
			}
#endif
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

#if defined cropLUchangeCtransfer
	if(transfer_litter_leaf) delete[] transfer_litter_leaf;
	if(transfer_litter_wood) delete[] transfer_litter_wood;
	if(transfer_litter_root) delete[] transfer_litter_root;
	if(transfer_litter_repr) delete[] transfer_litter_repr;
	if(transfer_harvested_products_slow) delete[] transfer_harvested_products_slow;
#endif
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

void harvest_pasture(double& cmass_leaf,double& cmass_root,double& litter_leaf,double& litter_root,double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{
	double turnover, residue_outtake, harvest;
	double scale=1.0;	//091005
	int m;
	bool alive=indiv.alive;

//	dprintf("In harvest_pasture before resetting carbon Year %d: cmass_root=%f, cmass_leaf=%f\n", date.year-nyear_spinup+1901,cmass_root,cmass_leaf);	//Last year's C
//pasture version of turnover adapted from harvest_crop 110525 ; presently just represents grass being harvested (for use with GRASSFORCROP); only yearly harvest !
//sparar inte skörden i någon variabel än !
//kan fås genom cmass_leaf*indiv.pft.harv_eff*2.0
//!!!!!!!!!!!!!!!!!!!!!!!!!!! NB last year's C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//scale harvest products of stands with increased area by (old area/new area) if landcover change has occurred:			//Review scaling !!!!!! NOT OK ?
	Stand& stand=indiv.vegetation.patch.stand;
	Gridcell& gridcell=stand.gridcell;

	if(gridcell.LC_updated)
	{
		scale=gridcell.landcoverfrac_old[PASTURE]/gridcell.landcoverfrac[PASTURE];		//fix 110102

		if(scale>=1.0)
			scale=1.0;

//		dprintf("In harvest_pasture(): year %d: scale=%f\n", date.year-nyear_spinup+1901,scale);
	}
///////////////////////////

//1.turnover and harvest of last year's carbon:

//NB. cmass_x can be negative here only if individuals with negative cmass-x are not killed last year.

	// Root turnover
	//Bondeau: turnover_root=0.5
	cmass_root*=scale;	//1108018
	cmass_leaf*=scale;	//1108018

//	turnover=indiv.pft.turnover_root*cmass_root*scale;	//turnover_root är normalt 0.7 för gräs
	turnover=indiv.pft.turnover_root*cmass_root;	//turnover_root är normalt 0.7 för gräs
	if(alive && turnover>0.0) 
		litter_root+=turnover;
	cmass_root-=turnover;

	//OBS ! skörd före turnover !!!!
	//Harvest/Grazing:					
	//Bondeau: harv_eff=0.9 i kod, 0.5 i artikel plus 0.05 till litter (faeces)
//	harvest=indiv.pft.harv_eff*cmass_leaf*scale;				//använd 0.5 (this year's yield is set in allocation_crop)
	harvest=indiv.pft.harv_eff*cmass_leaf;						//använd 0.5 (this year's yield is set in allocation_crop)

	if(ifslowharvestpool)	//Added 091007
	{
		harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;
		harvest=harvest*(1-indiv.pft.harvest_slow_frac);
	}
	acflux_harvest+=harvest;										//skördat gräs
	cmass_leaf-=harvest;

#if defined GRASSFORCROP
	if (alive && cmass_leaf>0.0)
	{
		residue_outtake=indiv.pft.res_outtake*cmass_leaf;				//res_outtake 0.75 for grassforcrop
		acflux_harvest+=residue_outtake;											//uttagna rester
		cmass_leaf-=residue_outtake;
	}
#endif
	// Leaf turnover
	turnover=indiv.pft.turnover_leaf*cmass_leaf;	//turnover_leaf är normalt 1.0 för gräs
	if(alive && turnover>0.0) 
		litter_leaf+=turnover;
	cmass_leaf-=turnover;
//	cmass_plant=cmass_leaf+cmass_root;

//2. distribution of transferred carbon calculated in landcover_dynamics():

}

void harvest_crop(double& cmass_plant,double& cmass_leaf,double& cmass_root,double& cmass_ho,double& cmass_agpool,double& litter_leaf,double& litter_root,
				  double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{	//NB. this function is for balancing carbon fluxes based on last year's cmass, not for calculating this year's yield. This is done in allocation_crop().
	double turnover, residue_outtake, harvest;
	double scale=1.0;	//091005
	int m;
	bool alive=indiv.alive;

//	dprintf("In harvest_crop before resetting carbon Year %d: cmass_plant=%f, cmass_root=%f, cmass_leaf=%f, cmass_ho=%f\n", date.year-nyear_spinup+1901,cmass_plant,cmass_root,cmass_leaf,cmass_ho);	//Last year's C
//CFT version of turnover 090602

//!!!!!!!!!!!!!!!!!!!!!!!!!!! NB last year's C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//scale harvest products of stands with increased area by (old area/new area) if landcover change has occurred:
	Stand& stand=indiv.vegetation.patch.stand;
	Gridcell& gridcell=stand.gridcell;

	if(gridcell.LC_updated)		//091005
	{
		scale=gridcell.cftfrac_old[stand.cftid]*gridcell.landcoverfrac_old[CROPLAND]/(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);	//100415
		if(scale>1.0)
			scale=1.0;

//		dprintf("In harvest_crop(): stand %d pftid=%d year %d: scale=%f\n", stand.id, indiv.pft.id, date.year-nyear_spinup+1901,scale);
	}

	cmass_root*=scale;	//121108	//Fix ger ca. 0.05 kg/m3 lägre flux på 100 år i Algeriet.
	cmass_leaf*=scale;	//121108
	cmass_agpool*=scale;//121108
	cmass_ho*=scale;	//121108

//1.turnover and harvest (and acflux_harvest) of last year's carbon :

//NB. cmass_x can be negative here only if individuals with negative cmass-x are not killed last year.
	if(indiv.pft.phenology==CROPGREEN)
	{		
		if(alive && cmass_root>0.0)
			litter_root+=cmass_root;

		cmass_root=0.0;

		//Bondeau: harv_eff=1.0		
		if(alive && cmass_ho>0.0)										// (this year's yield is set in allocation_crop)
		{
			harvest=indiv.pft.harv_eff*cmass_ho;						//skördade produkter	(möjlighet att använda överföring till gridcellpft.harvested_products_slow)

			if(indiv.pft.aboveground_ho)
				litter_leaf+=(cmass_ho-harvest);			//eller lägg i separat litter_ho ?	: ej skördade produkter
			else
				litter_root+=(cmass_ho-harvest);			//corrected 100609

			if(ifslowharvestpool)	//Added 091007
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
		if(indiv.get_cropinfo()->isintercropgrass)			//Intercrop growth
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

				if(ifslowharvestpool)	//Added 091007
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

			if(ifslowharvestpool)	//Added 091007
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
//2. distribution of transferred carbon calculated in landcover_dynamics():

}
