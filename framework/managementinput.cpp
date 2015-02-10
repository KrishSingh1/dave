

#include "input.h"

ManagementInputModule::ManagementInputModule(Input& in)
	: input(in), 
	  gridlist(in.gridlist) {

}

void ManagementInputModule::init() {

	if(run[CROPLAND]) {
#if defined DYNAMIC_LANDCOVER_INPUT

		if(readsowingdates)	{
			file_sdates=param["file_sdates"].str;
			if(!sdates.Open(file_sdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_sdates);
		}

		if(readharvestdates) {
			file_hdates=param["file_hdates"].str;
			if(!hdates.Open(file_hdates, gridlist))
				fail("initio: could not open %s for input",(char*)file_hdates);
		}

		if(readNfert) {
			file_Nfert=param["file_Nfert"].str;
			if(!Nfert.Open(file_Nfert, gridlist))
				fail("initio: could not open %s for input",(char*)file_Nfert);
		}

#endif
	}
}

bool ManagementInputModule::loadmanagement(Gridcell& gridcell, Coord c) {

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

bool ManagementInputModule::getgridcell(Gridcell& gridcell) {

	Coord& c=gridlist.getobj();
	bool LUerror = loadmanagement(gridcell, c);

	return LUerror;
}

void ManagementInputModule::getsowingdates(Gridcell& gridcell) {

#if defined DYNAMIC_LANDCOVER_INPUT
	if(!sdates.isloaded())
		return;
#endif

	int year = date.year - nyear_spinup + input.firsthistyear;

	if(date.year < nyear_spinup + input.nyear_hist) {
		for(int i=0; i<npft; i++) {
			if(pftlist[i].landcover == CROPLAND && pftlist[i].readsowingdate)	{
#if defined DYNAMIC_LANDCOVER_INPUT
				gridcell.pft[i].sdate_force = (int)sdates.Get(year,pftlist[i].name);
#endif
				// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
				for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
					Standpft& standpft = gridcell[j].pft[i];
					if(standpft.active)
						standpft.sdate_force = gridcell.pft[i].sdate_force;
				}
			}
		}
	}
}

void ManagementInputModule::getharvestdates(Gridcell& gridcell) {

#if defined DYNAMIC_LANDCOVER_INPUT
	if(!hdates.isloaded())
		return;
#endif

	int year = date.year - nyear_spinup + input.firsthistyear;

	if(date.year < nyear_spinup + input.nyear_hist) {
 		for(int i=0; i<npft; i++)	{
			if(pftlist[i].landcover == CROPLAND && pftlist[i].readharvestdate) {		
#if defined DYNAMIC_LANDCOVER_INPUT
				gridcell.pft[pftlist[i].id].hdate_force = (int)hdates.Get(year,pftlist[i].name);
#endif
				// Copy gridcellpft-value to standpft-value. If standtype values are required, modify code and input files.
				for(unsigned int j=0; j<gridcell.nbr_stands(); j++) {
					Standpft& standpft = gridcell[j].pft[i];
					if(standpft.active)
						standpft.hdate_force = gridcell.pft[i].hdate_force;
				}
			}
		}
	}
}

void ManagementInputModule::getNfert(Gridcell& gridcell) {

#if defined DYNAMIC_LANDCOVER_INPUT
	if(!Nfert.isloaded())
		return;
#endif

	int year = date.year - nyear_spinup + input.firsthistyear;

	if(date.year < nyear_spinup + input.nyear_hist) {
 		for(int i=0; i<npft; i++)	{
			if(pftlist[i].landcover == CROPLAND && pftlist[i].readNfert) {		
#if defined DYNAMIC_LANDCOVER_INPUT
				gridcell.pft[pftlist[i].id].Nfert_read = Nfert.Get(year,pftlist[i].name);
#endif
			}
		}
	}
}

void ManagementInputModule::getmanagement(Gridcell& gridcell) {

	if(run[CROPLAND]) {

		//Read sowing dates from input file, put into gridcellpft.sdate_force
		if(readsowingdates)		
			getsowingdates(gridcell);
		//Read harvest dates from input file, put into gridcellpft.hdate_force
		if(readharvestdates)		
			getharvestdates(gridcell);
		//Read N fertilization from input file, put into xxx
		if(readNfert)		
			getNfert(gridcell);
	}
}
