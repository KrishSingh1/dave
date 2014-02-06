//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Classes for input data (used for landcover/crop fractions). File format can be either line 1:lon lat, line 2: year data OR line 1: header, line 2: lon lat year data   //
// Mats Lindeskog 101019                                                                                                                                                //  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INDATA_H
#define INDATA_H

namespace InData {

//#include <algorithm>
using std::min;
using std::max;

#define GUESS_VERSION
#define MAXLINE 20000	//Ändrat från 400 091227
#define MAXNAMESIZE 50
#define MAXRECORDS 500
#define MAXLINESPARSE 30000

enum {EMPTY, GLOBAL_STATIC, GLOBAL_YEARLY, LOCAL_STATIC, LOCAL_YEARLY};

struct Coord 
{
	// Type for storing grid cell longitude, latitude and description text
	int id;
	double lon;
	double lat;
};

#ifndef GUESS_VERSION

#define dprintf printf
#define SUPPRESSLARGEOUTPUT true

struct Coord 
{
	// Type for storing grid cell longitude, latitude and description text
	double lon;
	double lat;
};

struct Neighbour
{
	int layer;
	double dist;
	Coord coord;
};

class NeighbourList
{
	double resolution;
	Coord currentStand;
	Neighbour *neighbours;

public:

	int nNeighbours;
	void SetCurrentStand(Coord c) {currentStand=c;}
	void SetNeighbours(int layers);
	Neighbour *GetNeighbours() const;
	void NeighbourList::GetNeighbours(Neighbour*) const;
	NeighbourList(double resolutionX=0.5);
//	NeighbourList(double resolution);
	~NeighbourList();

};


#endif


#ifndef GUESS_VERSION
class Gridlist
{
	FILE *ifp;
	char *fileName;
	Coord currentStand;
	bool isobj;
public:
	int Open(char* name);						//Returns 0 if error; opens file, checks format, sets fileName, nRecords and nYears and allocates memory for data[] and year[].
	bool firstobj()
	{
	}
	bool nextobj()
	{
//		if(!eof)
//			isobj=true;

	}
	Coord getobj()
	{
		return currentStand;
	}
};
#endif

class TimeDataD									//Represents a set of double data over time (years).
{												//Data can be global or for a specific stand. Also static. set by format flag
	FILE *ifp;
	char *fileName;
	int nCells;									//Set in ParseNCells() or ParseNCellsSpatial()
	bool ifheader;
	char header_arr[MAXRECORDS][MAXNAMESIZE];
	Coord currentStand;
	double *data;								//allocated in Allocate()
	bool *checkdata;							//allocated in CheckIfPresent()
	bool ischeckingdata;
	int firstyear;								//110601; set in ParseNYears() or ParseNYearsSpatial() to be used in FindRecord()
	bool isfirstgrid;
	bool loaded;

	int ParseFormat();							//Called from Open(); Returns 0 if wrong format, sets nRecords, ifheader and header_arr[]
	int ParseNYears();							//Called from Open()
	int ParseNYearsGlobal();					//Called from ParseNYears()
	int ParseNYearsLocal();						//Called from ParseNYears()
	int ParseNYearsSpatial();					//Called from OpenSpatial()
	int Allocate();								//Called from Open() or OpenSpatial()
	int FindRecord(Coord c) const;				//Quick version
	int FindRecord2(Coord c) const;				//Slower version, can handle blank lines
	void ParseNCells();
	int ParseNCellsSpatial();

public:
	int nRecords;								//Set in ParseFormat()
	int nYears;									//Set in ParseNYears()
	int *year;									//allocated in Allocate(), set in Load(), Load(Coord) or LoadNext()
	int format;									//EMPTY, GLOBAL_STATIC, GLOBAL_YEARLY, LOCAL_STATIC, LOCAL_YEARLY
	bool *active;								//allocated in Allocate()
	bool fileopened;
	TimeDataD(int format=EMPTY);				//default format value can only be used with header version input files !
	~TimeDataD();
	int Open(char* name);						//Returns 0 if error; opens file, checks format, sets fileName, nRecords and nYears and allocates memory for data[] and year[].
	int OpenSpatial(char* name, bool replace_original_file=false);
	void Close();
	int OutputConvertedSpatial(char*);
	int Load();									//Loads global data
	int Load(Coord c);							//Loads local data for a certain coordinate. Returns 0 if coordinate not found.
	int LoadNext();								//For stepping through a data file, loading each coordinate data consecutively. Returns 0 if error.
	void Output(char*);	
	double Get(int year, int column) const;		// Returns a single value
	double Get(int year, const char* name) const;		// Returns a single value for column with header string name. Returns -999 if name not found.
	int Get(int year, double* dataX) const;		//Copies the values for one year data to the dataX array, returns 0 if wrong format.
	int Get(double* dataX) const;				//Copies all data to the dataX array, returns 0 if wrong format.
//	double* Get(int year) const;
	int GetnRecords() const {return nRecords;}
	int GetHeader(char cropnames[][MAXNAMESIZE]) const;
	int GetHeaderFull(char *header_line) const;			//120124
	int FindCoord(Coord c) const{return FindRecord(c);}	//120124
//	int GetActive(bool *activeX) const;
	char* GetHeader(int record) const;
	Coord& GetCoord() {return currentStand;} //added 100106, added to GUESS version 120123	; updated to compile on Linux
	void Rewind() {rewind(ifp);}
	int GetNCells();
	int GetFirstyear();
	bool isloaded() { return loaded;}

#if defined GUESS_VERSION
	void CheckIfPresent(ListArray_id<Coord>& gridlist);
#endif
	bool CFTPresent(int cft){return checkdata[cft];}

};

//Class for loading data from memory instead of file.

class TimeDataDmem
{
	Coord *gridlist;
	double **data;
	int nCells;
	int nColumns;
	int nYears;
	bool ifheader;
	char header_arr[MAXRECORDS][MAXNAMESIZE];
	int currentCell;
	bool loaded;
public:
	double Get(int year, int column) const;			// Returns a single value.
	double Get(int year, const char* name) const;	// Returns a single value for column with header string name. Returns -999 if name not found.
	int Load(Coord c);	// Returns 0 if coordinate not found.
	void SetCoord(int index, Coord c);
	void SetData(int index, double* data);
	void Open(int nCells, int nColumns, int nYears);
	void Close();
	bool isloaded() { return loaded;}
	void CopyFromTimeDataD(TimeDataD& Data, ListArray_id<Coord>& gridlistX);
	TimeDataDmem();
	~TimeDataDmem();
};

//

class SoilData
{
	FILE *ifp;
	char *fileName;
	int soilcode;
	Coord currentStand;
public:
	SoilData();
	~SoilData();
	int Open(char* name);						//changed to return 1 at success 091227
	int Load(Coord c);
	void Output(char *name);
	int GetSoilcode(Coord c);
	void Rewind() {rewind (ifp);}
};

}

#endif//INDATA_H
