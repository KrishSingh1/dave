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
//		if(i!=CROPLAND) {					// cropland subclasses turned off in this version
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
	double landcoverfrac_change[NLANDCOVERTYPES];
//	double cropfrac_change[NCROPSTANDS_MAX];	
	double cropfrac_sum_old=0.0;
//	double cropstand_change[NCROPSTANDS_MAX];
	int nnaturalstands=0;

	bool LCchangeCtransfer=true;

	memset(landcoverfrac_change,0,NLANDCOVERTYPES*sizeof(double));
//	memset(cropfrac_change,0,NCROPSTANDS_MAX*sizeof(double));
//	memset(cropstand_change,0,NCROPSTANDS_MAX*sizeof(double));

	gridcell.LC_updated=false;

/////////////////////////////////////////////////////////////
//Landcover and cft fraction update (from updated landcoverfrac):

//Save old fraction values:										//Flytta till getlandcover?
	for(i=0;i<NLANDCOVERTYPES;i++)
		gridcell.landcoverfrac_old[i]=gridcell.landcoverfrac[i];
//	for(i=0;i<NCROPSTANDS_MAX;i++)
//		cropfrac_sum_old+=gridcell.cftfrac_old[i]=gridcell.cftfrac[i];

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
/*	if(run[CROPLAND] && (!cftfrac_fixed || !lcfrac_fixed))	
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
	else 
	{
#ifdef MATS_TEST
//		dprintf("\nYear %d: changeLC=%6.2f\tchangeCFT=%6.2f\n", date.year-nyear_spinup+1901, changeLC, change_crop);	
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
//orig						stand.natural_frac_change=landcoverfrac_change[NATURAL]*stand.frac/gridcell.landcoverfrac[NATURAL];
							stand.natural_frac_change=landcoverfrac_change[NATURAL]*stand.get_gridcell_fraction()/gridcell.landcoverfrac[NATURAL];
//							stand.natural_frac_change=landcoverfrac_change[NATURAL]/nnaturalstands;	//unsafe: not correct if the stand is smaller than the allotted change
//orig						stand.frac+=stand.natural_frac_change;
							stand.set_gridcell_fraction(stand.get_gridcell_fraction()+stand.natural_frac_change);
							break;
						}
						else
						{		
//orig						if(stand.frac>0.0)
							if(stand.get_gridcell_fraction()>0.0)
							{		
//orig							if(stand.frac>=-natural_change_remain)	//all natural landcover decrease is taken from this stand
								if(stand.get_gridcell_fraction()>=-natural_change_remain)	//all natural landcover decrease is taken from this stand
								{
									stand.natural_frac_change=natural_change_remain;
//orig								stand.frac+=stand.natural_frac_change;
									stand.set_gridcell_fraction(stand.get_gridcell_fraction()+stand.natural_frac_change);
									natural_change_remain=0.0;
									break;
								}
								else									//more stands will have to be reduced
								{
//orig								stand.natural_frac_change=-stand.frac;
									stand.natural_frac_change=-stand.get_gridcell_fraction();
//orig								natural_change_remain+=stand.frac;
									natural_change_remain+=stand.get_gridcell_fraction();
//orig								stand.frac=0.0;	//will be killed below
									stand.set_gridcell_fraction(0.0);	//will be killed below
								}				
							}
						}
					}

//					stand.natural_frac_change=landcoverfrac_change[NATURAL]*stand.frac/gridcell.landcoverfrac[NATURAL];		//To try different options
//					stand.natural_frac_change=landcoverfrac_change[NATURAL]/nnaturalstands;	
//					stand.frac+=stand.natural_frac_change;
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
//			if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL && landcoverfrac_change[stand.landcover]<0.0 
//				|| stand.landcover==NATURAL && landcoverfrac_change[NATURAL]<0.0 && (nnaturalstands==1 || stand.natural_frac_change<0.0)
//				|| stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)
			if(stand.landcover!=NATURAL && landcoverfrac_change[stand.landcover]<0.0 
				|| stand.landcover==NATURAL && landcoverfrac_change[NATURAL]<0.0 && (nnaturalstands==1 || stand.natural_frac_change<0.0))
#else
	//		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]<0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)
			if(landcoverfrac_change[stand.landcover]<0.0)
#endif
			{
//All landcovers that only have one stand:
#ifdef multiple_natural_stands
//				if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL)
				if(stand.landcover!=NATURAL)
#else
	//			if(stand.landcover!=CROPLAND)
#endif
				{
					scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
//orig				stand.frac=gridcell.landcoverfrac[stand.landcover];	//set frac here now !					
					stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);	//set frac here now !
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
//orig					stand.frac=gridcell.landcoverfrac[stand.landcover];	//set frac here now !
						stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);	//set frac here now !

#ifdef MATS_TEST
						dprintf("Natural stand n:o %d coverage decreased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
						dprintf("Natural stand n:o %d updated absolute fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);
#endif
					}
				}
#endif			
	/*			else if(stand.landcover==CROPLAND)
				{
					scale=-cropstand_change[stand.cftid]/receiving_fraction/(double)stand.nobj;
//orig				stand.frac=gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND];	//set frac here now !
					stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);	//set frac here now !
#ifdef MATS_TEST
					dprintf("Crop stand n:o %d coverage decreased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, cropstand_change[stand.cftid]);
					dprintf("Crop stand n:o %d updated absolute fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND]);
#endif
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
//INSERT HARVEST_PASTURE() HERE !
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

						if(stand.landcover==NATURAL)
						{
							double change_frac;
							if(nnaturalstands>1)
								change_frac=stand.natural_frac_change;
							else
								change_frac=landcoverfrac_change[NATURAL];
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
						Stand& stand=gridcell.createobj(gridcell,landcover,pftlist);
						stand.set_gridcell_fraction(landcoverfrac_change[i]);
//						newnaturalstand=true;

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
#ifdef multiple_natural_stands
//			if(stand.landcover!=CROPLAND && stand.landcover!=NATURAL)
			if(stand.landcover!=NATURAL)
#else
//			if(stand.landcover!=CROPLAND)	
#endif					
			{
				old_frac=gridcell.landcoverfrac_old[stand.landcover];
				added_frac=landcoverfrac_change[stand.landcover];
				new_frac=gridcell.landcoverfrac[stand.landcover];
//orig			stand.frac=gridcell.landcoverfrac[stand.landcover];	//set frac here now !
				stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);	//set frac here now !
#ifdef MATS_TEST
				dprintf("Stand n:o %d, type %d coverage increased year %d ! Difference = %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
				dprintf("Stand n:o %d, type %d updated fraction year %d: %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);
				dprintf("Stand n:o %d, type %d fraction of receiving area year %d: %f\n", stand.id, stand.landcover, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]/receiving_fraction);
#endif
			}
#ifdef multiple_natural_stands
			else if(stand.landcover==NATURAL)
			{
//				if(newnaturalstand)	//Put transferred C in new natural stand only.			Utkommenterat onödig kod.
				{
					if(stand.first_year==date.year)
					{
						old_frac=0.0;
						added_frac=landcoverfrac_change[stand.landcover];
						new_frac=landcoverfrac_change[stand.landcover];
					}
					else
					{
//						old_frac=0.0;
//						added_frac=0.0;
//						new_frac=1.0;	//To avoid dividing with zero; result is 0 anyway.
						gridcell.nextobj();
						continue;
					}					
				}
/*				else	//Borde hända endast vid skapandet av första beståndet !
				{
					if(nnaturalstands>1)
						dprintf("CODE TO BE ADDED HERE: if(nnaturalstands>1)\n");
					else
					{
						old_frac=gridcell.landcoverfrac_old[stand.landcover];
						added_frac=landcoverfrac_change[stand.landcover];
						new_frac=gridcell.landcoverfrac[stand.landcover];
					}
				}
*/
#ifdef MATS_TEST	
				dprintf("Natural stand n:o %d coverage increased year %d ! Difference = %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]);
//				dprintf("Natural stand n:o %d updated fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, gridcell.landcoverfrac[stand.landcover]);
				dprintf("Natural stand n:o %d updated fraction year %d: %f\n", stand.id, date.year-nyear_spinup+1901, stand.get_gridcell_fraction());
				dprintf("Natural stand n:o %d fraction of receiving area year %d: %f\n", stand.id, date.year-nyear_spinup+1901, landcoverfrac_change[stand.landcover]/receiving_fraction);
#endif
			}
#endif
/*			else if(stand.landcover==CROPLAND)
			{
				old_frac=gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.cftid];
				added_frac=cropstand_change[stand.cftid];
				new_frac=gridcell.landcoverfrac[CROPLAND]*gridcell.cftfrac[stand.cftid];
				stand.set_gridcell_fraction(gridcell.cftfrac[stand.cftid]*gridcell.landcoverfrac[CROPLAND])
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

#ifndef multiple_natural_stands
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		Stand& stand=gridcell.getobj();
			stand.set_gridcell_fraction(gridcell.landcoverfrac[stand.landcover]);
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
