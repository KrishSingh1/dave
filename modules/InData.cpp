//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Classes for input data (used for landcover/crop fractions). File format can be either line 1:lon lat, line 2: year data OR line 1: header, line 2: lon lat year data   //
// Mats Lindeskog 101019                                                                                                                                                //  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Error checking of input files revised 110601: Both formats should work now. Missing lines will cause fail() to be called at the opening of the file if MAXLINESPARSE is big enough, stopping the program. 
// If not, the site with the missing line will be omitted from the simulation.
// The program will inform about blank lines in some cases, but should handle them well, continuing the simulation of the site.
// NB. The header format is much slower when searching for coordinates in global files.

#include "config.h"
#include "guess.h" //??
//#include "gutil.h"
#include "InData.h"
//#include "math.h"

using namespace InData;

#ifndef GUESS_VERSION
bool ascendinglongitudes=false;
#else

bool ascendinglongitudes=false;	//Not true for randomised gridlists; set to false for now 130902.

bool TimeDataD::item_has_data(char* name){

	int column = GetColumn(name);

	if(column == -1)
		return false;
	else
		return checkdata[column];
}

bool TimeDataD::item_in_header(char* name) {

	if(GetColumn(name) == -1)
		return false;
	else
		return true;
};

void TimeDataD::CheckIfPresent(ListArray_id<Coord>& gridlist)	//Requires gutil.h
{
	if(checkdata)
	{
		delete[] checkdata;
		year=NULL;
	}

	checkdata=new bool[nRecords];
	memset(checkdata, 0, nRecords*sizeof(bool));
	ischeckingdata=true;

	gridlist.firstobj();
	while(gridlist.isobj)
	{
		Coord& c=gridlist.getobj();
		if(Load(c))
		{
			for(int i=0;i<nYears;i++)
			{
				for(int j=0;j<nRecords;j++)
				{
					if(data[i*nRecords+j]>0.0)
						checkdata[j]=1;
				}
			}
if(!SUPPRESSLARGEOUTPUT)
dprintf("Done checking for crop pft:s in input file at %.2f,%.2f.\n", c.lon, c.lat);
		}
		gridlist.nextobj();
	}
/*
	dprintf("\n");
	for(int j=0;j<nRecords;j++)
	{
		dprintf("%s:%d\n", header_arr[j], checkdata[j]);
	}
*/
	gridlist.firstobj();
	rewind(ifp);
	ischeckingdata=false;
}
#endif

int TimeDataD::GetHeader(char cropnames[][MAXNAMESIZE]) const
{
	if(ifheader && header_arr)
	{
		for(int i=0; i<nRecords; i++)
			strncpy(cropnames[i], header_arr[i], MAXNAMESIZE*sizeof(char));
		return 1;
	}
	else
		return 0;
}

int TimeDataD::GetHeaderFull(char *header_line) const	
{
	if(ifheader && header_arr)
	{
		if(format==LOCAL_YEARLY)
			strcpy(header_line, "   Lon\t   Lat\t  Year");
		else if(format==GLOBAL_YEARLY)
			strcpy(header_line, "    lon\t    lat");
		else if(format==LOCAL_STATIC)
			strcpy(header_line, "   year");

		for(int i=0; i<nRecords; i++)
		{
			char buffer[MAXNAMESIZE];
			sprintf(buffer, "\t%8s", header_arr[i]);
			strncat(header_line, buffer, strlen(buffer));
		}
		return 1;
	}
	else
		return 0;
}

char* TimeDataD::GetHeader(int record) const
{
	if(ifheader && header_arr)
	{
//		strcpy(cropname, header_arr[record]);
		return (char*)header_arr[record];
	}
	else
		return 0;	
}

int TimeDataD::Get(double* dataX) const
{
	if(format==LOCAL_YEARLY || format==LOCAL_STATIC ||format==GLOBAL_YEARLY || format==GLOBAL_STATIC)
		memcpy(dataX, data, nYears*nRecords * sizeof(double));
	else
	{
		printf("Wrong usage of TimeDataD::Get(int year, double* dataX).\n");
		return 0;
	}

	return 1;
}

int TimeDataD::CalenderYearToPosition(int calender_year) const {

	int year = calender_year - firstyear;

	// Use first or last year's data if calender year is not within data period.
	if(year < 0)
		year = 0;
	else if(year >= nYears)
		year = nYears -1;
	else
		year = calender_year - firstyear;

	return year;
}

int TimeDataD::Get(int calender_year, double* dataX) const
{
	int yearX = CalenderYearToPosition(calender_year);

	memcpy(dataX, &data[yearX*nRecords], nRecords * sizeof(double));

	return 1;
}

double TimeDataD::Get(int calender_year, int column) const
{
	double dataX=0.0;
	int yearX = CalenderYearToPosition(calender_year);

	if(column>=nRecords)
	{
//		if(yearX==1)	//Set to 1 in crop branch, was 0.
//			printf("WARNING: Trying to retreive more columns than available in %s. Value set to 0.0 \n", fileName);
		return 0.0;
	}

	dataX=data[nRecords*yearX+column];

//printf("yearX=%d, nYears=%d\n",yearX, nYears);	//test 100719
	return dataX;
}

double TimeDataD::Get(int calender_year, const char* name) const		//Returns a single value for column with header string name
{
	int column=-1;
	double dataX=-999;

	for(int i=0;i<nRecords;i++)
	{
		if(!strcmp(name, header_arr[i]))
		{
			column=i;
			break;
		}
	}

	if(column==-1)
	{
		if(calender_year == firstyear + 1)	//Set to 1 in crop branch, was 0.
		printf("WARNING: Value for %s not found in %s. Value set to 0.0\n", name, fileName);
	}
	else
		dataX=Get(calender_year, column);

	return dataX;
}

int TimeDataD::GetColumn(const char* name) const
{
	int column = -1;

	for(int i=0; i<nRecords; i++)
	{
		if(!strcmp(name, header_arr[i]))
		{
			column = i;
			break;
		}
	}

	if(column == -1)
	{
		printf("WARNING: Data for %s not found in %s.\n", name, fileName);
		return -1;
	}
	else
		return column;
}

int TimeDataD::Open(char* name)
{
	int format_parsed=EMPTY;

	if(ifp)
	{
		fclose(ifp);
		ifp=NULL;
		fileopened=false;
	}
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
	}

	if(name)
		ifp=fopen(name, "r");
	else
		return 0;

	if(ifp)
	{
//		dprintf("\nOpened input file %s\n", name);	// Test
		fileName=new char[strlen(name)+1];
		fileopened=true;

		if(!fileName)
		{
			printf("Cannot allocate memory for file name string !\n");
			return 0;
		}
		else
			strcpy(fileName,name);

		format_parsed=ParseFormat();

		if(format!=format_parsed)	// Checks format (sets it if header), sets nRecords, ifheader and header_arr[]
		{
			printf("Wrong format in file %s (failing ParseFormat()!\n", name);
			return 0;
		}
		else if(format==GLOBAL_YEARLY || format==LOCAL_YEARLY)
		{
			nYears=ParseNYears();	// Parse numbers of years in input file
//			printf("nYears:%d\n", nYears);
			if(nYears==0)
			{
				printf("Wrong format in file %s (nYears=0)!\n", name);
				return 0;
			}
		}
		else if(format==GLOBAL_STATIC || format==LOCAL_STATIC)
		{
			nYears=1;
		}
		else if(format==EMPTY)	//should be set by now
		{
			printf("Please set data format at initialization !\n");
			return 0;
		}

		if(!Allocate())				//Allocate memory for dynamic data
		{
			printf("Could not allocate memory for data from file %s!\n", name);
			return 0;
		}
	}
	else
	{
		printf("TimeDataD::Open: File %s could not be opened for input !\n\n", name);
		return 0;
	}

	return 1;
}

int TimeDataD::OpenSpatial(char* name, bool replace_original_file)
{
	int format_parsed=EMPTY;

	if(ifp)
	{
		fclose(ifp);
		ifp=NULL;
		fileopened=false;
	}
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
	}

	if(name)
		ifp=fopen(name, "r");
	else
		return 0;

	if(ifp)
	{
		dprintf("\nOpened input file %s\n", name);	// Test
		fileName=new char[strlen(name)+1];
		fileopened=true;

		if(!fileName)
		{
			printf("Cannot allocate memory for file name string !\n");
			return 0;
		}
		else
			strcpy(fileName,name);

		format_parsed=ParseFormat();

		if(format!=format_parsed)	// Checks format (sets it if header), sets nRecords, ifheader and header_arr[]
		{
			printf("Wrong format in file %s (failing ParseFormat()!\n", name);
			return 0;
		}
		else if(format==GLOBAL_YEARLY || format==LOCAL_YEARLY)
		{
			nCells=ParseNCellsSpatial();
			nYears=ParseNYearsSpatial();	// Parse numbers of years in input file
			printf("nCells:%d\n", nCells);
			printf("nYears:%d\n", nYears);
			if(nYears==0 || nCells==0)
			{
				printf("Wrong format in file %s (nYears=0)!\n", name);
				return 0;
			}
		}
		else if(format==GLOBAL_STATIC || format==LOCAL_STATIC)
		{
			printf("Static datasets should not be opened by  OpenSpatial(), quitting !\n");
			return 0;
		}
		else if(format==EMPTY)	//should be set by now
		{
			printf("Please set data format at initialization !\n");
			return 0;
		}

		if(!Allocate())				//Allocate memory for dynamic data
		{
			printf("Could not allocate memory for data from file %s!\n", name);
			return 0;
		}
	}
	else
	{
		printf("TimeDataD::Open: File %s could not be opened for input !\n\n", name);
		return 0;
	}

	int size=0;
	char *chp=NULL, tag[]="_converted", ext[]=".txt", *outFileName=NULL;


	size=strlen(name)+strlen(tag)+1;
	outFileName=new char[size];
	memset(outFileName, 0, size); 
	strncpy(outFileName, name, strlen(name));
	if(chp=strrchr(outFileName, '.'))
	{
		*chp='\0';
	}
	strncat(outFileName, tag, strlen(tag));
	strncat(outFileName, ext, strlen(ext));

	if(!OutputConvertedSpatial(outFileName))
		return 0;

	Close();
	if(replace_original_file)
	{
		remove(name);
		rename(outFileName, name);
		Open(name);
	}
	else
		Open(outFileName);

	return 1;
}

int TimeDataD::ParseFormat()	//Checks format, sets nRecords, ifheader and header_arr[].
{ //Desired format must be set beforehand by program at initiation of class TimeDataD objects (if no header) !

	char line[MAXLINE], *p=NULL, s1[MAXRECORDS][MAXNAMESIZE]={'\0'}, s2[MAXRECORDS][MAXNAMESIZE]={'\0'};
	int count1=0, count2=0, i=0, k=0;
	int format_local=EMPTY, offset=0;
	float d[MAXRECORDS]={0.0};

//	dprintf("Parsing format for file %s \n", fileName);

//First line: 
	do			// Just in case there is a blank line at the beginning...
	{
		if(fgets(line,sizeof(line),ifp))
		{
			p=strtok(line, "\t\n ");

			if(!p)				//Fix for blank line 110531
				continue;

			strncpy(s1[count1], p, MAXNAMESIZE-1);
			count1++;
			do
			{
				p=strtok(NULL, "\t\n ");
				if(p)
				{
					strncpy(s1[count1], p, MAXNAMESIZE-1);
					count1++;
				}
				k++;
			}
			while(p);
/*
			dprintf("Line 1 count:%d\n", count1);
			for(i=0;i<count1;i++)
				dprintf("%s\t", s1[i]);
			dprintf("\n");
*/
			p=NULL;
		}
		else return 0;
	}
	while(!(count1>0));

	if(!strcmp(s1[0], "lon") || !strcmp(s1[0], "Lon") || !strcmp(s1[0], "LON"))
	{
		if(!strcmp(s1[2], "year") || !strcmp(s1[2], "Year"))
		{
			format_local=LOCAL_YEARLY;
			offset=2;
			for(i=3;i<count1;i++)
				strncpy(header_arr[i-3],s1[i], MAXNAMESIZE-1);
		}
		else
		{
			format_local=LOCAL_STATIC;
//			offset=2;
			for(i=2;i<count1;i++)
				strncpy(header_arr[i-2],s1[i], MAXNAMESIZE-1);
		}
	}
	else if(!strcmp(s1[0], "year") || !strcmp(s1[0], "Year"))
	{
			format_local=GLOBAL_YEARLY;
			for(i=1;i<count1;i++)
				strncpy(header_arr[i-1],s1[i], MAXNAMESIZE-1);
	}
	else if(!strcmp(s1[0], "static"))
	{
			format_local=GLOBAL_STATIC;
			for(i=1;i<count1;i++)
				strncpy(header_arr[i-1],s1[i], MAXNAMESIZE-1);
	}
	else
		ifheader=false;

/*	if(ifheader)
	{
		dprintf("header:\n");
		for(i=0;i<count1 && *(header_arr[i])!='\0';i++)
			dprintf("%s\t", header_arr[i]);
		dprintf("\n");
	}
*/

//Second line:
	do			// Just in case there is a blank line at the beginning...
	{
		if(fgets(line,sizeof(line),ifp))
		{
			p=strtok(line, "\t\n ");

			if(!p)				//Fix for blank line 110531
				continue;

			strncpy(s2[count2], p, 9);
			count2++;
			do
			{
				p=strtok(NULL, "\t\n ");
				if(p)
				{
					strncpy(s2[count2], p, 9);
					count2++;
				}
			}
			while(p);
/*
			dprintf("Line 2 count:%d\n", count2);
			for(i=0;i<count2;i++)
				dprintf("%s\t", s2[i]);
			dprintf("\n");
*/
		}
		else return 0;
	}
	while(!(count2>0));

	rewind(ifp);

	if(format==EMPTY)
	{
		if(ifheader)
			format=format_local;
		else
			printf("Please set data format at initialization !\n");
	}

	switch (format)
	{
	case GLOBAL_YEARLY:
		if(format_local==GLOBAL_YEARLY || count1>1 && count1==count2)
		{
			nRecords=count2-1;
//			printf("Format in input file is compatible with GLOBAL_YEARLY flag\n");
			dprintf("nRecords:%d\n", nRecords);
			return GLOBAL_YEARLY;
		}
		else
		{
			printf("Format in input file is incompatible with GLOBAL_YEARLY flag\n");
			return 0;
		}
		break;
	case LOCAL_STATIC:
		if(format_local==LOCAL_STATIC || count1>2 && count1==count2)
		{
			nRecords=count2-2;
//			printf("Format in input file is compatible with LOCAL_STATIC flag\n");
//			dprintf("nRecords:%d\n", nRecords);
			return LOCAL_STATIC;
		}
		else
		{
			printf("Format in input file is incompatible with LOCAL_STATIC flag\n");
			return 0;
		}
		break;
	case LOCAL_YEARLY:
		if(format_local==LOCAL_YEARLY || count1==2 && count2>1)
		{
			nRecords=count2-1-offset;
//			printf("Format in input file is compatible with LOCAL_YEARLY flag\n");
//			dprintf("nRecords:%d\n", nRecords);
			return LOCAL_YEARLY;
		}
		else
		{
			printf("Format in input file is incompatible with LOCAL_YEARLY flag\n");
			return 0;
		}
		break;
	case GLOBAL_STATIC:
		if(format_local==GLOBAL_STATIC || count1>1 && count1==count2)
		{
			nRecords=count2-1;
//			printf("Format in input file is compatible with GLOBAL_STATIC flag\n");
			dprintf("nRecords:%d\n", nRecords);
			return GLOBAL_STATIC;
		}
	default:	// format EMPTY
		printf("Format is not set correctly in file %s !\n", fileName);
		return 0;
	}
}

int TimeDataD::ParseNYears()
{
	int n_yearsX=0;

	switch (format)
	{
	case GLOBAL_YEARLY:
//		dprintf("Parsing nYears for GLOBAL_YEARLY format file\n");	// Test
		n_yearsX=ParseNYearsGlobal();
		break;
	case LOCAL_YEARLY:
//		printf("Parsing nYears for LOCAL_YEARLY format in file %s\n", fileName);	// Test
		n_yearsX=ParseNYearsLocal();
		break;
	default:
		printf("Format in is uncorrectly set by program for file %s !\n", fileName);
		return 0;
	}

	return n_yearsX;
}

int TimeDataD::GetNCells()
{
	if(!nCells)
		ParseNCells();

	return
		nCells;
}

int TimeDataD::GetFirstyear()
{
	return
		firstyear;
}

void TimeDataD::ParseNCells()
{
	float d1;
	long int oldpos;
	int i=0, count=0;
	char line[MAXLINE];
	bool error=false;

	oldpos=ftell(ifp);
	if(oldpos!=0)
		rewind(ifp);

	if(ifheader)
	{
		fgets(line,sizeof(line),ifp);	//ignore header line
	}

	while(!feof(ifp))
	{
//		count=0;
		line[0]=0;

		fgets(line,sizeof(line),ifp);
		count=sscanf(line,"%f", &d1);

		if(count>0)
			i++;
	}

	if(ifheader)
	{
		nCells=i/nYears;
		if(i%nYears)
			error=true;
	}
	else
	{
		nCells=i/(nYears+1);
		if(i%(nYears+1))
			error=true;
	}

	if(error)
		dprintf("Unexpected number of lines ! No.lines=%d, No.cells=%d, No.years=%d\n", i,nCells,nYears);

	fseek(ifp, oldpos, 0);

}

int TimeDataD::ParseNYearsSpatial()	//Adapted from ParseNCells()
{
	float d1;
	long int oldpos;
	int i=0, count=0;
	char line[MAXLINE];
	bool error=false;

	if(!ifheader)
	{
		printf("ParseNYearsSpatial() not functional for files without header !\n");
		return 0;
	}

	oldpos=ftell(ifp);
	if(oldpos!=0)
		rewind(ifp);

	if(ifheader)
	{
		fgets(line,sizeof(line),ifp);	//ignore header line
	}

	while(!feof(ifp))
	{
//		count=0;
		line[0]=0;

		fgets(line,sizeof(line),ifp);
		count=sscanf(line,"%f", &d1);

		if(count>0)
			i++;
	}

	if(ifheader)
	{
//		nCells=i/nYears;
		nYears=i/nCells;
		if(i%nCells)
			error=true;
	}

	fseek(ifp, oldpos, 0);

	if(error)
	{
		dprintf("Unexpected number of lines ! No.lines=%d, No.cells=%d, No.years=%d\n", i,nCells,nYears);
		return 0;
	}
	else
		return nYears;
}

int TimeDataD::ParseNCellsSpatial()	//Adapted from ParseNYearsLocal()
{
	int i=0, count1=0, prevLine=0, ncells1=0, ncells2=0, n=0;
	char line[MAXLINE];
	bool new_year=false;
	float d1=0,d2=0,d3=0, d1_prevLine=0, d2_prevLine=0, d3_prevLine=0;

//	printf("Inside ParseNCellsSpatial()\n");	// Test

	if(!ifheader)
	{
		printf("ParseNCellsSpatial() not functional for files without header !\n");
		return 0;
	}

	for(i=0;i<150000 && !feof(ifp);)
	{
		if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
		{
			count1=sscanf(line,"%f%f%f", &d1,&d2,&d3);		//does not count header strings !
			if(count1>0)			// avoids blank lines
			{
				if(ifheader && (d3!=d3_prevLine))	// First line of new year
				{
					ncells2=i-prevLine;
					new_year=true;
					if(i<2)
						firstyear=(int)d3;	//110601
				}

				if(new_year)
				{
					if((ncells1!=ncells2) && n>1)
					{
						printf("FORMAT ERROR in input file %s !\n", fileName);
						return 0;
					}
					ncells1=ncells2;		//NB. not set if input file has data for only one coordinate !
					prevLine=i;
					n++;

					new_year=false;
				}

				i++;

				d1_prevLine=d1;
				d2_prevLine=d2;
				d3_prevLine=d3;
			}
		}
	}

	if(feof(ifp))	//Sista lokalen !
	{
		if(ifheader)
			ncells2=i-prevLine;
//		else
//			ncells2=i-prevLine-1;
		if((ncells1!=ncells2) && n>1)
		{
			printf("FORMAT ERROR in input file %s !\n", fileName);
			return 0;
		}
	}
//dprintf("nyears=%d\n", nyears1);
	rewind(ifp);
	return ncells2;			//fix 101125
}

int TimeDataD::ParseNYearsLocal()
{
	int i=0, count1=0, prevLine=0, nyears1=0, nyears2=0, n=0;
	char line[MAXLINE];
	bool new_coord=false;
	float d1=0,d2=0,d3=0, d1_prevLine=0, d2_prevLine=0;

//	printf("Inside ParseNYearsLocal()\n");	// Test

	for(i=0;i<MAXLINESPARSE && !feof(ifp);)
	{
		if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
		{
			count1=sscanf(line,"%f%f%f", &d1,&d2,&d3);		//does not count header strings !
			if(count1>0)			// avoids blank lines
			{
				if(ifheader && (d1!=d1_prevLine || d2!=d2_prevLine))	// First line of new coordinate
				{
					nyears2=i-prevLine;
					new_coord=true;
					firstyear=(int)d3;	//110601
				}
//				if(count1==2)
				else if(count1==2 && d1<=180.0)		//line with coordinates	; added new condition to be able to use data files with only one value per year 100721
				{
					nyears2=i-prevLine-1;
					new_coord=true;
				}

				if(new_coord)
				{
					if((nyears1!=nyears2) && n>1)
					{
						printf("FORMAT ERROR in input file %s !\n", fileName);
						return 0;
					}
					nyears1=nyears2;		//NB. not set if input file has data for only one coordinate !
					prevLine=i;
					n++;

					new_coord=false;
//dprintf("d1=%.2f, d1=%.2f\n", d1, d2);
				}
				else if(!ifheader && i==prevLine+1)
					firstyear=(int)d1;	//110601

				i++;

				d1_prevLine=d1;
				d2_prevLine=d2;
			}
		}
	}

	if(feof(ifp))	//Sista lokalen !
	{
		if(ifheader)
			nyears2=i-prevLine;
		else
			nyears2=i-prevLine-1;
		if((nyears1!=nyears2) && n>1)
		{
			printf("FORMAT ERROR in input file %s !\n", fileName);
			return 0;
		}
	}
//dprintf("nyears=%d\n", nyears1);
	rewind(ifp);
	return nyears2;			//fix 101125
}

int TimeDataD::ParseNYearsGlobal()
{
	int count=0, nyears=0;
	char line[MAXLINE];
	float d1=0,d2=0,d3=0;

//	dprintf("Inside ParseNYearsGlobal()\n");	// Test

	if(ifheader)
		fgets(line,sizeof(line),ifp);

	while(!feof(ifp))
	{
		if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
		{
			count=sscanf(line,"%f%f%f", &d1,&d2,&d3);
			if(count>0)	// Ignore blank lines
			{
				if(count>=2)
				{
					nyears++;

					if(nyears==1)				//110601
						firstyear=(int)d1;
				}
				else
				{
					printf("FORMAT ERROR in input file %s !\n", fileName);
					nyears=0;
					break;
				}
			}
		}
	}
	rewind(ifp);
	return nyears;
}

int TimeDataD::Allocate()	// Allocates memory for dynamic data: format & nYears must be set before !
{
	if(year)
	{
		delete[] year;
		year=NULL;
	}
	if(data)
	{
		delete[] data;
		data=NULL;
	}
//	printf("\nAllocating memory for data in TimeDataD::Allocate()\n\n");

	switch(format)
	{
	case EMPTY:
		break;
	case GLOBAL_STATIC:
		year=new int;
		data=new double[nRecords];
		if(year)
			*year=0;
		if(data)
			*data=0;
		break;
	case GLOBAL_YEARLY:
		year=new int[nYears];
		data=new double[nRecords*nYears];
		if(year)
			memset(year, 0, nYears*sizeof(int));
		if(data)
			memset(data, 0, nRecords*nYears*sizeof(double));
		break;
	case LOCAL_STATIC:
		year=new int;
		data=new double[nRecords];
		if(year)
			*year=0;
		if(data)
			memset(data, 0, nRecords*sizeof(double));
		break;
	case LOCAL_YEARLY:
		year=new int[nYears];
		data=new double[nRecords*nYears];
		if(year)
			memset(year, 0, nYears*sizeof(int));
		if(data)
			memset(data, 0, nRecords*nYears*sizeof(double));
		break;
	default:
		;
	}
	if(year && data)
		return 1;
	else
		return 0;
}

int TimeDataD::Load()	// for GLOBAL_YEARLY and GLOBAL_STATIC data
{
	int i=0, count=0, yearX=0, yearX_previous, k=0;
	char line[MAXLINE], *p=NULL;
	double d1=0.0;
	double d[MAXRECORDS]={0.0};
	float extra=0.0;
	bool error=0;

//	dprintf("Loading all data from %s into memory\n", fileName);	// Test

	if(ifp)
	{
		if(format==GLOBAL_YEARLY)
		{
			if(year)
				memset(year, 0, nYears*sizeof(int));
			if(data)
				memset(data, 0, nRecords*nYears*sizeof(double));

			yearX_previous=firstyear-1;

			if(ifheader)
				fgets(line, sizeof(line), ifp);

			for(i=0;i<nYears;)
			{
				k=0;
				if(fgets(line, sizeof(line), ifp))
				{
					memset(d, 0, nRecords*sizeof(double));
					count=0;

					p=strtok(line, "\t\n ");	//year
					if(!p)				//Fix for blank line 110531
						continue;

					sscanf(p, "%d", &yearX);
//printf("count=%d\n",count);
//printf("yearX=%d\n",yearX);
					if(yearX!=yearX_previous+1)		//110607
					{
						printf("FORMAT ERROR in input file %s: Load(). Wrong year in data file ! Missing line ?\n", fileName);
						error=1;
						break;
					}
					else
						yearX_previous=yearX;

					do
					{
						p=strtok(NULL, "\t\n ");
						if(p)
						{
							count+=sscanf(p, "%lf", &d[k]);
//printf("count=%d\n",count);
//printf("d[%d]=%f\n",k, d[k]);
						}
						k++;
					}
					while(p);

					if(count>0)
					{
						if(count==nRecords)
						{
							year[i]=yearX;
							for(int j=0;j<nRecords;j++)
							{
								data[nRecords*i+j]=d[j];
							}
						}
						else
						{
							printf("FORMAT ERROR in input file %s: Load(), count!=%d, year %d\n", fileName,i+1);
							error=1;
						}
						i++;	// only count lines with something on them
					}
				}
				else
				{
					printf("An ERROR occurred reading file %s\n", fileName);
					error=1;
				}
			}
		}
		else if(format==GLOBAL_STATIC)
		{
			if(fgets(line, sizeof(line), ifp))
			{
				if(ifheader)
				{
					if(fgets(line, sizeof(line), ifp))
						p=strtok(line," \t");	//"static"
					else
					{
						printf("An ERROR occurred reading file %s\n", fileName);
						error=1;
					}
				}
				do
				{
					p=strtok(NULL, "\t\n ");
					if(p)
					{
						count+=sscanf(p, "%lf", &d[k]);
//printf("count=%d\n",count);
//printf("d[%d]=%f\n",k, d[k]);
					}
					k++;
				}
				while(p);

				if(count==nRecords)
				{
					for(i=0;i<nRecords;i++)
						data[i]=d[i];
				}
				else
				{
					printf("FORMAT ERROR in input file %sf: Load(), count!=%d\n", fileName,nRecords+1);
					error=1;
				}
			}
			else
			{
				printf("An ERROR occurred reading file %s\n", fileName);
				error=1;
			}
		}
		else
		{
			printf("Wrong usage of Load(void)\n");
			error=1;
		}
	}
	else
	{
		printf("Cannot load from unopened file !\n");
		error=1;
	}

	if(ifp)
	{
		fclose(ifp);
		ifp=NULL;
	}

	if(error) {
		loaded = false;
		return 0;
	}
	else {
		loaded = true;
		return 1;
	}
}

int TimeDataD::Load(Coord c)
{
	char line[MAXLINE], *p=NULL;
	int i=0, j=0, k=0, count1=0, nyears=0, yearX=0, yearX_previous;
	float lonX=0.0, latX=0.0;
	double d[MAXRECORDS]={0.0};
	bool error=0;

//	dprintf("Inside Load(Coord)\n");	// Test

	if(ifp)
	{
		if(format==LOCAL_YEARLY)
		{
			if(FindRecord(c))
			{
				if(year)
					memset(year, 0, nYears*sizeof(int));
				if(data)
					memset(data, 0, nRecords*nYears*sizeof(double));

				yearX_previous=firstyear-1;

				currentStand.lon=c.lon;
				currentStand.lat=c.lat;

//				for(i=0;i<nYears ;)
				while(i<nYears && yearX<firstyear+nYears-1)
				{
					k=0;
					if(fgets(line, sizeof(line), ifp))
					{
						memset(d, 0, nRecords*sizeof(double));
						count1=0;

						if(ifheader)
						{
							p=strtok(line," \t");	//lon
							sscanf(p, "%f", &lonX);			//110607
							p=strtok(NULL, " \t");	//lat
							sscanf(p, "%f", &latX);			//110607
							p=strtok(NULL, " \t");	//year

							if(lonX!=c.lon || latX!=c.lat)	//110607
							{
								printf("FORMAT ERROR in input file %s for stand at Coordinate %.2f,%.2f: Load(). Wrong coordinates in data file !\n", fileName,c.lon,c.lat);
								error=1;
								break;
							}
						}
						else
							p=strtok(line, "\t\n ");	//year

						if(!p)				//Fix for blank line 110531
							continue;

						sscanf(p, "%d", &yearX);
//printf("count1=%d\n",count1);
//printf("yearX=%d\n",yearX);
						if(yearX!=yearX_previous+1)		//110607
						{
							printf("FORMAT ERROR in input file %s for stand at Coordinate %.2f,%.2f: Load(). Wrong year in data file ! Missing line ?\n", fileName,c.lon,c.lat);
							error=1;
//							break;	//Don't break if missing line. Next search will begin at the start of the next coordinate.
						}

						yearX_previous=yearX;

						do
						{
							p=strtok(NULL, "\t\n ");
							if(p)
							{
									count1+=sscanf(p, "%lf", &d[k]);
//printf("count1=%d\n",count1);
//printf("d[%d]=%f\n",k, d[k]);
							}
							k++;
						}
						while(p);
//printf("count1=%d\n",count1);
						if(count1>0)
						{
							if(count1==nRecords)
							{
								year[i]=yearX;
								for(j=0;j<nRecords;j++)
								{
									data[nRecords*i+j]=d[j];
								}
							}
							else
							{
								printf("FORMAT ERROR in input file %s for stand at Coordinate %.2f,%.2f: Load(), count!=%d, year %d\n", fileName,c.lon,c.lat,nRecords+1,i+1);
								error=1;
								break;
							}
							i++;	// only count lines with something on them
						}
					}
					else
					{
						printf("An ERROR occurred reading file %s\n", fileName);
						error=1;
						break;
					}
				}
			}
			else
			{
				printf("COULD NOT FIND DATA for %.2f, %.2f in file %s\n",c.lon,c.lat,fileName);
				error=1;
			}
		}
		else if(format==LOCAL_STATIC)
		{
			if(FindRecord(c))
			{
				if(data)
					memset(data, 0, nRecords*sizeof(double));
				currentStand.lon=c.lon;
				currentStand.lat=c.lat;
				if(fgets(line, sizeof(line), ifp))
				{
					memset(d, 0, nRecords*sizeof(double));
					p=strtok(line," \t");	//lon
					sscanf(p, "%f", &lonX);
					p=strtok(NULL, " \t");	//lat
					sscanf(p, "%f", &latX);

					do
					{
						p=strtok(NULL, "\t\n ");
						if(p)
						{
							count1+=sscanf(p, "%lf", &d[k]);
//printf("count1=%d\n",count1);
//printf("d[%d]=%f\n",k, d[k]);
						}
						k++;
					}
					while(p);

					if(count1>0)
					{
						if(count1==nRecords)
						{
							for(j=0;j<nRecords;j++)
							{
								data[j]=d[j];
							}
						}
						else
						{
							printf("FORMAT ERROR in input file %s for stand at Coordinate %.2f,%.2f: Load(), count!=%d, year %d\n", fileName,c.lon,c.lat,nRecords+1,i+1);
							error=1;
						}
//						i++;	// only count lines with something on them
					}
				}
				else
				{
					printf("An ERROR occurred reading file %s\n", fileName);
					error=1;
				}
			}
			else
			{
				printf("COULD NOT FIND DATA for %.2f, %.2f in file %s\n",c.lon,c.lat,fileName);
				error=1;
			}
		}
		else
		{
			printf("Wrong usage of Load(Coord)\n");
			error=1;
		}
	}
	else
	{
		printf("Cannot load from unopened file !\n");
		error=1;
	}

	if(error)
		return 0;
	else
	{
if(!SUPPRESSLARGEOUTPUT)
		dprintf("Loading all data for %.2f,%.2f in %s into memory\n", c.lon, c.lat,fileName);
		return 1;
	}
}

int TimeDataD::LoadNext()	//Only implemented for LOCAL_YEARLY (100106) and LOCAL_STATIC (121016)	; Needs to be modified to handle missing lines in data files with header ! (see Load)
{

	char line[MAXLINE], *p=NULL;
	int count=0, i=0, j=0, k=0, count1=0, nyears=0, yearX=0, n=0;
	double d1, d2, d3, d[MAXRECORDS]={0.0};
	bool error=0, firstyear=true;
	long int fpos;

	if(ifp && !feof(ifp))
	{
		if(format==LOCAL_YEARLY)
		{
			if(year)
				memset(year, 0, nYears*sizeof(int));
			if(data)
				memset(data, 0, nRecords*nYears*sizeof(double));

			if(ifheader)
			{
				fpos=ftell(ifp);
				if(fpos==0)
					fgets(line,sizeof(line),ifp);	//ignore header line
			}

			if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
			{
				count=sscanf(line,"%lf%lf%lf", &d1, &d2, &d3);
				if(count>0)	// Avoid blank lines at the end of the file
				{
					if(count==2 || count>2 && (format==LOCAL_STATIC || ifheader))	//added LOCAL_STATIC compatibility 091207 (not used, use SoilData instead for soilcode)
					{																//Bugfix 110607 (was =2)
						currentStand.lon=d1;
						currentStand.lat=d2;				
					}
					else
					{
						printf("FORMAT ERROR in input file %s: LoadNext(), count!=2, line %d\n",fileName,i);
						error=1;
					}
				}
				else
				{
					printf("WARNING: blank line in file %s: LoadNext(), count==0\n",fileName);
				}	
			}
			else
				error=1;

			for(i=0;i<nYears && error==0;)
			{
				k=0;
				count1=0;

				if(ifheader && firstyear)
					firstyear=false;
				else
					fgets(line, sizeof(line), ifp);

				if(line)
				{
					memset(d, 0, nRecords*sizeof(double));

					if(ifheader)
					{
						p=strtok(line," \t");	//lon
						p=strtok(NULL, " \t");	//lat
						p=strtok(NULL, " \t");	//year

						//Kolla att koordinaten är samma här !
					}
					else
						p=strtok(line, "\t\n ");	//year
					if(!p)							//Fix for blank line 110531
						continue;
					sscanf(p, "%d", &yearX);
//printf("count1=%d\n",count1);
//printf("yearX=%d\n",yearX);
					//Kolla att årtalet är rätt här !

					do
					{
						p=strtok(NULL, "\t\n ");
						if(p)
						{
							count1+=sscanf(p, "%lf", &d[k]);
//printf("count1=%d\n",count1);
//printf("d[%d]=%f\n",k, d[k]);
						}
						k++;
					}
					while(p);
//printf("count1=%d\n",count1);

					if(count1>0)
					{
						if(count1==nRecords)
						{
							year[i]=yearX;
							for(j=0;j<nRecords;j++)
							{
								data[nRecords*i+j]=d[j];
							}
						}
						else
						{
							printf("FORMAT ERROR in input file %s: LoadNext(), count!=%d, year %d\n", fileName,nRecords+1,i+1);
							error=1;
							break;
						}
						i++;	// only count lines with something on them
					}
				}
				else
				{
					printf("An ERROR occurred reading file %s\n", fileName);		//dailytomonthly Fastnar här !
					error=1;
					break;
				}
			}
		}
		else if(format==LOCAL_STATIC)
		{
			if(data)
				memset(data, 0, nRecords*nYears*sizeof(double));

			if(ifheader)
			{
				fpos=ftell(ifp);
				if(fpos==0)
					fgets(line,sizeof(line),ifp);	//ignore header line
			}

			if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
			{
				count=sscanf(line,"%lf%lf%lf", &d1, &d2, &d3);
				if(count>0)	// Avoid blank lines at the end of the file
				{
					if(count==2 || count>2 && (format==LOCAL_STATIC || ifheader))	//added LOCAL_STATIC compatibility 091207 (not used, use SoilData instead for soilcode)
					{																//Bugfix 110607 (was =2)
						currentStand.lon=d1;
						currentStand.lat=d2;				
					}
					else
					{
						printf("FORMAT ERROR in input file %s: LoadNext(), count!=2, line %d\n",fileName,i);
						error=1;
					}
				}
				else
				{
					printf("WARNING: blank line in file %s: LoadNext(), count==0\n",fileName);
				}	
			}
			else
				error=1;

			for(i=0;i<nYears && error==0;)
			{
				k=0;
				count1=0;

				if(ifheader && firstyear)
					firstyear=false;
				else
					fgets(line, sizeof(line), ifp);

				if(line)
				{
					memset(d, 0, nRecords*sizeof(double));

					if(ifheader)
					{
						p=strtok(line," \t");	//lon
						p=strtok(NULL, " \t");	//lat
					}

					if(!p)							//Fix for blank line 110531
						continue;

//printf("count1=%d\n",count1);
//printf("yearX=%d\n",yearX);
					//Kolla att årtalet är rätt här !

					do
					{
						p=strtok(NULL, "\t\n ");
						if(p)
						{
							count1+=sscanf(p, "%lf", &d[k]);
//printf("count1=%d\n",count1);
//printf("d[%d]=%f\n",k, d[k]);
						}
						k++;
					}
					while(p);
//printf("count1=%d\n",count1);

					if(count1>0)
					{
						if(count1==nRecords)
						{
							for(j=0;j<nRecords;j++)
							{
								data[nRecords*i+j]=d[j];
							}
						}
						else
						{
							printf("FORMAT ERROR in input file %s: LoadNext(), count!=%d, year %d\n", fileName,nRecords+1,i+1);
							error=1;
							break;
						}
						i++;	// only count lines with something on them
					}
				}
				else
				{
					printf("An ERROR occurred reading file %s\n", fileName);		//dailytomonthly Fastnar här !
					error=1;
					break;
				}
			}
		}
	}
	else
		error=1;

	if(error)
	{
//		if(feof(ifp))
//			printf("End of file reached for file %s\n", fileName);
		return 0;
	}
	else
	{
//		dprintf("Data from %s loaded for coordinate <%.2f,%.2f>\n",fileName,c.lon,c.lat);
		return 1;
	}
}

//Fast version. Can not handle blank lines in some cases, will call FindRecord2() in those cases.
int TimeDataD::FindRecord(Coord c) const
{

	int i=0, count=0, n=0, lap=0, line_no=0;
	char line[MAXLINE], *p=NULL;
	double d1=0.0,d2=0.0,d3=0.0;
	bool found=0, error=0, start=true;
	long int oldpos;

//	printf("Inside FindRecord(), looking for coordinate <%.3f, %.3f> in %s\n", c.lon, c.lat, fileName);	// Test

	do
	{	
		i=0;
		start=true;							//110531

		while(!feof(ifp))
		{
			if(ifheader)
			{
//				oldpos=ftell(ifp);				// slow !

				if(!(i%nYears))					// much quicker, only 3 times per coordinate.
					oldpos=ftell(ifp);
			}


			if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
			{
				if(ifheader)
				{
					if(start==true)
					{
						start=false;

						if(oldpos==0)
							continue;
					}
				}	

				if(!(i%(nYears+1)) && !ifheader || !(i%nYears) && ifheader)		//fix 101111
				{
					count=sscanf(line,"%lf%lf%lf", &d1, &d2, &d3);
					if(count>0)	// Avoid blank line at the end of the file
					{
						if(count==2 || count>2 && (d3==firstyear || format==LOCAL_STATIC) && (format==LOCAL_STATIC || ifheader))	//added LOCAL_STATIC compatibility 091207 (not used, use SoilData instead for soilcode)
//						if(count==2 || count>2 && d3==firstyear && (format==LOCAL_STATIC || ifheader))	//added LOCAL_STATIC compatibility 091207 (not used, use SoilData instead for soilcode)
						{
							if(c.lon==d1 && c.lat==d2)
							{
								if(!ischeckingdata)
if(!SUPPRESSLARGEOUTPUT)
									dprintf("Coordinate <%.2f,%.2f> found in %s\n", d1, d2, fileName);
								found=1;
								break;
							}
							else if(ascendinglongitudes && c.lon<d1)		//Ny kod 091126: longitudes must be ascending in dataset for this to work !	
							{
								dprintf("c.lon<d1; rewinding...\n");
								break;
							}
						}
						else	
						{ 
							if(ifheader)
								dprintf("FORMAT ERROR in input file %s: FindRecord(), wrong firstyear, line %d\n",fileName,i);
							else
								dprintf("FORMAT ERROR in input file %s: FindRecord(), count!=2, line %d\n",fileName,i);
							error=1;
							break;
						}
					}
					else
					{
						dprintf("WARNING: blank line in file %s: FindRecord(), count==0\n",fileName);
						continue;
					}
				}		
				i++;
			}
		}

		if(!found)
		{
			lap++;
			rewind(ifp);
			if(error)
				break;
			else
				dprintf("Rewinding and searching from the beginning of the file...\n");
		}

	}while(!found && lap<2);

	if(found && !error)
	{
		if(ifheader)
			fseek(ifp, oldpos, 0);	//The found line needs to be read again in Load()
//			fseek(ifp, oldpos-newpos, SEEK_CUR);

		return 1;
	}
	else
		return FindRecord2(c);	//If not found, try FindRecord2()
}

//This version should handle blank or missing lines at all positions.
int TimeDataD::FindRecord2(Coord c) const	//No need for FindRecord2()
{

	int i=0, count=0, n=0, lap=0, lastyear;
	char line[MAXLINE], *p=NULL;
	double d1=0.0,d2=0.0,d3=0.0;
	bool found=0, error=0, start=true;
	long int oldpos;

//	printf("Inside FindRecord2(), looking for coordinate <%.3f, %.3f> in %s\n", c.lon, c.lat, fileName);	// Test

	lastyear=firstyear+nYears-1;

	do
	{	
		i=0;
		start=true;							//110531

		while(!feof(ifp))
		{
			if(ifheader)
			{
//				oldpos=ftell(ifp);				// Slow !

				if(d3==lastyear || start)  // Quicker !
					oldpos=ftell(ifp);				
			}

			if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
			{
				if(ifheader)
				{
					if(start==true)
					{
						start=false;

						if(oldpos==0)
							continue;
					}
				}	

				count=sscanf(line,"%lf%lf%lf", &d1, &d2, &d3);

				if(count>0)	// Avoid blank lines at the end of the file
				{
//					if(!(i%(nYears+1)) && !ifheader || !(i%nYears) && ifheader)		//fix 101111	; commented out to cope with missing lines
					{
						if(count==2 || count>2 && d3==firstyear && (format==LOCAL_STATIC || ifheader))	//added LOCAL_STATIC compatibility 091207 (not used, use SoilData instead for soilcode)
						{
							if(c.lon==d1 && c.lat==d2)
							{
								if(!ischeckingdata)
if(!SUPPRESSLARGEOUTPUT)
									dprintf("Coordinate <%.2f,%.2f> found in %s\n", d1, d2, fileName);
								found=1;
								break;
							}
							else if(ascendinglongitudes && c.lon<d1)		//Ny kod 091126: longitudes must be ascending in dataset for this to work !	
							{
								dprintf("c.lon<d1; rewinding...\n");
								break;
							}
						}
/*						else																		// commented out to cope with missing lines
						{ 
							dprintf("FORMAT ERROR in input file %s: FindRecord2(), count!=2, line %d\n",fileName,i);
							error=1;
							break;
						}
*/					}
				}
				else
				{
					dprintf("WARNING: blank line in file %s: FindRecord2(), count==0\n",fileName);
					continue;
				}	
				i++;
			}
		}

		if(!found)
		{
			lap++;
			rewind(ifp);
			if(lap<2)
				dprintf("Rewinding and searching from the beginning of the file...\n");
		}

	}while(!found && lap<2);

	if(found && !error)
	{
		if(ifheader)
			fseek(ifp, oldpos, 0);
		return 1;
	}
	else
		return 0;
}

void TimeDataD::Output(char *name)
{
	int i=0, j=0;
	FILE *ofp;

//	dprintf("Inside Output()\n");	// Test
	if(isfirstgrid)
		remove(name);

	if(format==GLOBAL_STATIC || format==GLOBAL_YEARLY)
		ofp=fopen(name, "w");
	else if(format==LOCAL_STATIC || format==LOCAL_YEARLY)
		ofp=fopen(name, "a");

	if(ifheader && header_arr && isfirstgrid)
	{
		switch (format)
		{
		case GLOBAL_STATIC:
			break;
		case GLOBAL_YEARLY:
			fprintf(ofp, "  year\t");
			break;
		case LOCAL_STATIC:
			fprintf(ofp, "   lon\t   lat\t");
			break;
		case LOCAL_YEARLY:
			fprintf(ofp, "   lon\t   lat\t  year\t");
			break;
		default:
			;
		}
			   
		for(int i=0; i<nRecords; i++)
			fprintf(ofp, "%8s\t", header_arr[i]);
		fprintf(ofp, "\n");
		isfirstgrid=false;
	}


	switch (format)
	{
	case GLOBAL_STATIC:
//		ofp=fopen(name, "w");
		fprintf(ofp, "%.3lf\n", *data);
		break;
	case GLOBAL_YEARLY:
//		ofp=fopen(name, "w");
		for(i=0;i<nYears;i++)
			fprintf(ofp, "%d\t%.3lf\n", year[i], data[i]);
		break;
	case LOCAL_STATIC:
		fprintf(ofp, "%6.2f\t%6.2f",currentStand.lon, currentStand.lat);
		for(j=0;j<nRecords;j++)
			fprintf(ofp, "\t%8.3f", data[nRecords*i+j]);
		fprintf(ofp, "\n");
		break;
	case LOCAL_YEARLY:
//		ofp=fopen(name, "a");
		
		if(!ifheader)
			fprintf(ofp, "%8.2f\t%8.2f\n",currentStand.lon, currentStand.lat);

		for(i=0;i<nYears;i++)
		{
			if(ifheader)
				fprintf(ofp, "%6.2f\t%6.2f\t",currentStand.lon, currentStand.lat);

			fprintf(ofp, "%6d ", year[i]);
			for(j=0;j<nRecords;j++)
				fprintf(ofp, "\t%8.3f", data[nRecords*i+j]);
			fprintf(ofp, "\n");
		}
		break;
	default:
		;
	}
	if(ofp)
		fclose(ofp);
}

int TimeDataD::OutputConvertedSpatial(char* outfilename)
{

	double *data_full=NULL, d[MAXRECORDS]={0.0};
	float d1=0,d2=0,d3=0, d1_prevLine=0, d2_prevLine=0, d3_prevLine=0;
	Coord *c_x=NULL;
	int count1=0, k=0, i=0, prevLine=0, cell_no=0, ncells1=0, ncells2=0, yearX, year_no=0;;
	char line[MAXLINE], *p=NULL;
	bool error=false;
	FILE *ofp_cft=NULL;

	if(!ifheader || format!=LOCAL_YEARLY)
	{
		printf("OutputConvertedSpatial() only works on files with LOCAL_YEARLY format with header.\n");

	}

	ofp_cft=fopen(outfilename, "w");

	data_full=new double[nRecords*nYears*nCells];
	if(data_full)
		memset(data_full, 0, nRecords*nYears*nCells*sizeof(double));
	else
	{
		printf("OutputConvertedSpatial(): not enough memory to create data structure, quitting.\n");
		return 0;
	}

	c_x=new Coord[nCells];
	
	rewind(ifp);

	while(!feof(ifp))
	{
		line[0]=0;

		if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
		{
			k=0;

			count1=sscanf(line,"%f%f%f", &d1,&d2,&d3);		//does not count header strings !

			if(count1>0)			// avoids blank lines
			{

				p=strtok(line," \t");	//lon
//				sscanf(p, "%f", &lonX);
				p=strtok(NULL, " \t");	//lat
//				sscanf(p, "%f", &latX);
				p=strtok(NULL, " \t");	//year
				sscanf(p, "%d", &yearX);
				if(!p)							//Fix for blank line 110531
					continue;

				count1=0;
				do
				{
					p=strtok(NULL, "\t\n ");
					if(p)
					{
						count1+=sscanf(p, "%lf", &d[k]);
					}
					k++;
				}
				while(p);

				if(d3!=d3_prevLine)	//New year.
				{
					ncells2=i-prevLine;
					if((ncells1!=ncells2) && year_no>1)
					{
						printf("FORMAT ERROR in input file %s !\n", fileName);
						error=true;
						break;
					}
					ncells1=ncells2;		//NB. not set if input file has data for only one coordinate !
					prevLine=i;
					cell_no=0;
					year_no++;
				}

				if(count1==nRecords)
				{
					for(int j=0;j<nRecords;j++)
					{
						data_full[cell_no*nYears*nRecords+nRecords*(yearX-firstyear)+j]=d[j];
					}
				}
				else
				{
					printf("FORMAT ERROR in input file %s: LoadNext(), count!=%d, year %d\n", fileName,nRecords+1,i+1);
					error=true;
					break;
				}

				if(year_no<2)
				{
					c_x[cell_no].lon=d1;
					c_x[cell_no].lat=d2;
				}
				else if(c_x[cell_no].lon!=d1 || c_x[cell_no].lat!=d2)
				{			
					printf("FORMAT ERROR in input file %s !\n", fileName);
					error=true;
					break;
				}

				i++;
				cell_no++;

				d1_prevLine=d1;
				d2_prevLine=d2;
				d3_prevLine=d3;
			}
		}
	}

	if(ofp_cft && !error)
	{
		fprintf(ofp_cft, "%8s%8s%8s", "Lon", "Lat", "year");
		for(int i=0;i<nRecords;i++)
			fprintf(ofp_cft, "%9s", GetHeader(i));
		fprintf(ofp_cft, "\n");

		for(int c=0;c<nCells;c++)
		{
			for(int y=0;y<nYears;y++)
			{
				fprintf(ofp_cft, "%8.1f%8.1f%8d", c_x[c].lon, c_x[c].lat, firstyear+y);
				for(int i=0;i<nRecords;i++)
					fprintf(ofp_cft, "%9f", data_full[c*nYears*nRecords+nRecords*y+i]);
				fprintf(ofp_cft, "\n");	
			}
		}
	}

	if(data_full)
		delete[] data_full;
	if(c_x)
		delete[] c_x;

	if(ofp_cft)
		fclose(ofp_cft);

	if(error)
		return 0;
	else
		return 1;
}

// Constructor
TimeDataD::TimeDataD(int formatX)
{
	ifp=NULL;
	fileName=NULL;
	ifheader=true;
	memset(header_arr,0,sizeof(char)*MAXRECORDS*MAXNAMESIZE);
//	for(int i=0;i<MAXRECORDS;i++)
//		printf("header_arr[i]=%s\n", header_arr[i]);
	currentStand.lon=0;
	currentStand.lat=0;
	data=NULL;
	checkdata=NULL;
	ischeckingdata=false;
	isfirstgrid=true;

	nRecords=0;
	nYears=0;
	nCells=0;
	year=NULL;
	format=formatX;
	fileopened=false;
}

//Deconstructor
TimeDataD::~TimeDataD()
{
//	printf("\nIn TimeDataD destructor\n");	
	if(ifp)
	{
		fclose(ifp);
		printf("Closing input file %s \n", fileName);	// Test
	}
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
//		printf("deleting fileName[] in TimeDataD::~TimeDataD()\n");	// Test
	}
	if(year)
	{
		delete[] year; 
		year=NULL;
//		printf("deleting year[] in TimeDataD::~TimeDataD()\n");	// Test
	}
	if(data)
	{
		delete[] data; 
		data=NULL;
//		printf("deleting data[] in TimeDataD::~TimeDataD()\n");	// Test
	}
	if(checkdata)
	{
		delete []checkdata;
		checkdata=NULL;
	}
}

void TimeDataD::Close()
{
	if(ifp)
	{
		fclose(ifp);
		printf("Closing input file %s \n", fileName);	// Test
	}
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
	}
	if(year)
	{
		delete[] year; 
		year=NULL;
	}
	if(data)
	{
		delete[] data; 
		data=NULL;
	}
	if(checkdata)
	{
		delete []checkdata;
		checkdata=NULL;
	}
}

int TimeDataDmem::CalenderYearToPosition(int calender_year) const {

	int yearX = calender_year - firstyear;

	// Use first or last year's data if calender year is not within data period.
	if(yearX < 0)
		yearX = 0;
	else if(yearX >= nYears)
		yearX = nYears -1;
	else
		yearX = calender_year - firstyear;

	return yearX;
}

double TimeDataDmem::Get(int calender_year, int column) const
{

	int yearX = CalenderYearToPosition(calender_year);

	if(currentCell >= 0 && column < nColumns)
		return data[currentCell][yearX * nColumns + column];
	else
		return 0.0;
}

double TimeDataDmem::Get(int calender_year, const char* name) const		//Returns a single value for column with header string name
{
	int column = -1;
	double dataX = -999;

	for(int i=0; i<nColumns; i++)
	{
		if(!strcmp(name, header_arr[i]))
		{
			column = i;
			break;
		}
	}

	if(column == -1)
	{
		if(calender_year == firstyear + 1)	//Set to 1 in crop branch, was 0.
		printf("WARNING: Value for %s not found in input file\n", name);
	}
	else
		dataX = Get(calender_year, column);

	return dataX;
}

int TimeDataDmem::Load(Coord c)
{
	bool error=true;

	if(gridlist[currentCell+1].lon==c.lon && gridlist[currentCell+1].lat==c.lat)	//In case gridlist cell order is same as in land use files.
	{
		currentCell++;
		error=false;
	}
	else
	{
		for(int i=0;i<nCells;i++)
		{
			if(gridlist[i].lon==c.lon && gridlist[i].lat==c.lat)
			{
				currentCell=i;
				error=false;
				break;
			}
		}
	}
	if(error) {
		loaded = false;
		return 0;
	}
	else {
		loaded = true;
		return 1;
	}
}
void TimeDataDmem::SetData(int index, double* dataX)
{
	if(data && data[index])
		memcpy(data[index], dataX, nColumns*nYears* sizeof(double));
}

void TimeDataDmem::SetCoord(int index, Coord c)
{
	gridlist[index].lon=c.lon;
	gridlist[index].lat=c.lat;
}

void TimeDataDmem::Open(int nCellsX, int nColumnsX, int nYearsX)
{
	nCells=nCellsX;
	nColumns=nColumnsX;
	nYears=nYearsX;
	gridlist=new Coord[nCells];
	data=new double*[nCells];
	for(int i=0;i<nCells;i++)
	{
		data[i]=new double[nColumns*nYears];
		if(data[i])
			memset(data[i], 0, nColumns*nYears*sizeof(double));
	}
}

void TimeDataDmem::Close()
{
	nCells = 0;
	nColumns = 0;
	nYears = 0;

	if(gridlist)
	{
		delete []gridlist;
		gridlist=NULL;
	}
	for(int i=0;i<nCells;i++)
	{
		if(data[i])
			delete[] data[i];
	}
	if(data)
	{
		delete[] data; 
		data=NULL;
	}
}

void TimeDataDmem::CopyFromTimeDataD(TimeDataD& Data, ListArray_id<Coord>& gridlistX)
{
	int cell_no=0;

	if(Data.GetHeader(header_arr))
		ifheader = true;

	firstyear = Data.GetFirstyear();

	double *celldata;
	celldata=new double[Data.nRecords*Data.nYears];

	while(Data.LoadNext()) {
		InData::Coord c;
		c=Data.GetCoord();

		gridlistX.firstobj();
		while(gridlistX.isobj) {
			Coord cc=gridlistX.getobj();
			if(c.lon==cc.lon && c.lat==cc.lat) {
				SetCoord(cell_no, c);
				Data.Get(celldata);
				SetData(cell_no, celldata);
				cell_no++;
				break;
			}
			gridlistX.nextobj();
		}
	}
	delete[] celldata;
}

int TimeDataDmem::GetFirstyear()
{
	return
		firstyear;
}

TimeDataDmem::TimeDataDmem()
{
	gridlist=NULL;
	data=NULL;
	nCells=0;
	ifheader=false;
	memset(header_arr,0,sizeof(char)*MAXRECORDS*MAXNAMESIZE);
	currentCell=-1;
}

TimeDataDmem::~TimeDataDmem()
{
	for(int i=0;i<nCells;i++)
	{
		if(data[i])
			delete[] data[i];
	}
	if(data)
		delete[] data;
	if(gridlist)
		delete[] gridlist;
}


///	class SoilData member function
SoilData::SoilData()
{
	ifp=NULL;
	fileName=NULL;
	soilcode=-1;
	currentStand.lon=0;
	currentStand.lat=0;
}

SoilData::~SoilData()
{
	if(ifp)
		fclose(ifp);
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
	}
}

void SoilData::Output(char *name)
{
	FILE *ofp;

//	dprintf("\nInside SoilData::Output()\n");	// Test
	ofp=fopen(name, "a");

	fprintf(ofp, "%7.2f\t%7.2f\t",currentStand.lon, currentStand.lat);
	fprintf(ofp, "%d\n", soilcode);

	if(ofp)
		fclose(ofp);
}

int SoilData::Open(char *name)
{
	if(ifp)
	{
		fclose(ifp);
		ifp=NULL;
	}
	if(fileName)
	{
		delete []fileName;
		fileName=NULL;
	}
	ifp=fopen(name, "r");

	if(ifp)
	{
//		dprintf("Opened file %s for input in SoilData::Open()\n\n", name);	// Test
		fileName=new char[strlen(name)+1];
		if(!fileName)
			exit(1);
		strcpy(fileName,name);
	}
	else
	{
		printf("File %s could not be opened for input ! Leaving program...\n\n", name);
		return 0;
	}

	return 1;
}

int SoilData::Load(Coord c)
{
	int i=0, count=0, lap=0, soil=-1;
	char line[MAXLINE];
	double d1=0.0,d2=0.0;
	float extra=0.0;
	bool found=0;

	currentStand.lon=c.lon;
	currentStand.lat=c.lat;

//	dprintf("Inside SoilData::Load(), looking for coordinate <%.3f, %.3f> in %s\n", c.lon, c.lat, fileName);	// Test
	if(ifp)
	{
		do
		{	
			i=0;
			while(!feof(ifp))
			{
				if(fgets(line,sizeof(line),ifp))	//OBS! behövs, annars läser scanf in gammal line igen efter sista raden
				{
					count=sscanf(line,"%lf%lf%d%f", &d1, &d2, &soil, &extra);
					if(count>0)	// skip empty lines
					{
						if(count==3)
						{
							if(c.lon==d1 && c.lat==d2)
							{
//								dprintf("Coordinate %6.2f,%6.2f found\n", d1, d2);
								found=1;
								soilcode=soil;
								break;
							}
						}
						else
							printf("Warning ! FORMAT ERROR in input file: wrong number of values on line %d\n",i);
					}
				i++;				
				}
			}

			if(!found)
			{
				lap++;
				rewind(ifp);
//				dprintf("Searching from the beginning of the file...\n");
			}

		}while(!found && lap<2);
	}
	else
	{
		printf("Cannot load from unopened file !\n");
		exit(1);
	}

	if(found)
	{
		printf("Data from %s loaded for coordinate <%.2f,%.2f>\n",fileName,c.lon,c.lat);
		return 1;
	}
	else
		return 0;
}

int SoilData::GetSoilcode(Coord c)
{
	int soilcodeX=0;

	if(currentStand.lon == c.lon && currentStand.lat == c.lat)
		soilcodeX=soilcode;
	else
	{
		if(Load(c))
		{
			soilcodeX=soilcode;
		}
		else
		{
			printf("No soilcode for this coordinate !\n");
			exit(1);
		}
	}

	return soilcodeX;
}
