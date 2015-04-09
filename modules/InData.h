////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file indata.h
/// \brief Classes for text input data (used mainly for landcover input).											  
/// File format can be either line 1:lon lat, line 2: year data OR line 1: header, line 2: lon lat year data  
/// \author Mats Lindeskog
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INDATA_H
#define INDATA_H

#define GUESS_VERSION

#ifdef GUESS_VERSION
#include "inputdefinitions.h"
#define MAPFILE	// Mapping of input file data when LUTOMEMORY not defined
#else
#include <stdio.h>
#include <math.h>
#include <algorithm>
#endif

using std::min;
using std::max;

namespace InData 
{

#define MAXLINE 20000
#define MAXNAMESIZE 50
#define MAXRECORDS 500
#define MAXLINESPARSE 30000
#define NOTFOUND -999
#define MAX_SEARCHRADIUS 1.0

/// Formats in text input file
enum {EMPTY, GLOBAL_STATIC, GLOBAL_YEARLY, LOCAL_STATIC, LOCAL_YEARLY};

/// Type for storing grid cell longitude, latitude and associated data position on disk
struct CoordPos 
{
	double lon;
	double lat;
	long int pos;
};

#ifndef GUESS_VERSION

#define dprintf printf

/// Type for storing grid cell longitude and latitude
struct Coord {

	double lon;
	double lat;
};

#endif

// Forward declaration of TimeDataDmem
class TimeDataDmem;

/// Class for reading a set of double data at coordinate positions over time (years), alternatively static or/and global.
class TimeDataD	{

	// PRIVATE VARIABLES	

	/// File pointer to input file
	FILE *ifp;
	/// File name
	char *fileName;
	/// Number of data columns
	int nColumns;								//Set in ParseFormat()
	/// Number of data years
	int nYears;									//Set in ParseNYears()
	/// Number of data gridcells
	int nCells;									//Set in ParseNCells() or ParseNCellsSpatial()
	// Spacial resolution of input data in degrees
	double spatial_resolution;
	/// Whether the input file structure includes a header line with column names and coordinates on each line of data
	bool ifheader;
	/// String array of data column names
	char header_arr[MAXRECORDS][MAXNAMESIZE];
	/// Coordinates for current gridcell
	Coord currentStand;
	/// Pointer to data for one (current) gridcell
	double *data;								//allocated in Allocate(), set in Load(), Load(Coord) or LoadNext()
	/// Pointer to array with data years
	int *year;									//allocated in Allocate(), set in Load(), Load(Coord) or LoadNext()
	/// Pointer to array with indication whether data column contains values > 0 or not
	bool *checkdata;							//allocated in CheckIfPresent()
	/// Format of input data
	int format;									//EMPTY, GLOBAL_STATIC, GLOBAL_YEARLY, LOCAL_STATIC, LOCAL_YEARLY
	/// First year of data
	int firstyear;								//set in ParseNYears() or ParseNYearsSpatial()
	/// Whether data is currently being checked by CheckIfPresent()
	bool ischeckingdata;
	/// Whether data file is opened
	bool fileopened;
	/// Whether data for the requested coordinates have been found and loaded
	bool loaded;

	/// Pointer to memory copy of all data for the gridlist
	TimeDataDmem *memory_copy;
	/// Pointer to map of file positions of data for all gridcells in the file
	CoordPos *filemap;

	// PRIVATE METHODS

	/// Methods for parsing input file data format and structure
	int ParseFormat();							//Called from Open(); Returns 0 if wrong format, sets nColumns, ifheader and header_arr[]
	int ParseNYears();							//Called from Open()
	int ParseNYearsGlobal();					//Called from ParseNYears()
	int ParseNYearsLocal();						//Called from ParseNYears()
	int ParseNYearsSpatial();					//Called from OpenSpatial()
	void ParseNCells();
	int ParseNCellsSpatial();
	double ParseSpatialResolution();			//Called from Open()

	/// Allocates memory for dynamic data structures
	int Allocate();								//Called from Open() or OpenSpatial()
	/// Finds data for a gridcell in input file. Quick version
	int FindRecord(Coord c) const;
	/// Finds data for a gridcell in input file. Slower version, can handle blank lines
	int FindRecord2(Coord c) const;
	/// Converts data column name to data column index
	int GetColumn(const char* name) const;		// Returns column number for header name.
	/// Converts calender year to valid year position in data array.
	int CalenderYearToPosition(int calender_year) const;	
	/// Creates map of the file positions of all gridcells' data
	void CreateFileMap();
	/// Sets the file pointer to required position (found in the file map)
	void SetPosition(long int pos) {fseek(ifp, pos, 0);}
	/// Rewinds the file pointer
	void Rewind() {rewind(ifp);}
	/// Loads local data for a certain coordinate from a file map. Returns 0 if coordinate not found.
	int LoadFromMap(Coord c);

#ifdef GUESS_VERSION
	/// Copies all data for the specified gridlist to memory
	void CopyToMemory(int ncells, ListArray_id<Coord>& lonlatlist);
#endif

public:

	// PUBLIC METHODS

	/// Constructor
	TimeDataD(int format=EMPTY);		// default format value can only be used with header version input files !
	/// Deconstructor
	~TimeDataD();

// Methods to open input files and access data
	/// Opens input file, checks format and allocates memory. Returns 0 if error
	int Open(char* name);
#ifdef GUESS_VERSION
	/// Opens input file, checks format and allocates memory. Copies all data for the gridlist into memory if LUTOMEMORY is defined. Returns 0 if error.
	int Open(char* name, ListArray_id<Coord>& gridlist);
#endif
	/// Releases dynamically allocated memory.
	void Close();
	/// Writes the data of the current coordinate to an output file
	void Output(char* outfile);
	/// Loads global data
	int Load();
	/// Loads data for a certain coordinate. Returns 0 if coordinate not found.
	int Load(Coord c, double offset = 0.0);
	/// Steps through a data file, loading each coordinate's data consecutively. Returns 0 if error.
	int LoadNext(long int *pos = NULL);
	/// Returns a single data value for a certain year and data column
	double Get(int calender_year, int column) const;
	/// Returns a single data value for column with header string name. Returns -999 if name not found.
	double Get(int calender_year, const char* name) const;
	/// Copies the data for the current gridcell for one year to an array, returns 0 if wrong format.
	int Get(int calender_year, double* dataX) const;
	/// Copies all data for the current gridcell to an array, returns 0 if wrong format.
	int Get(double* dataX) const;

// Methods used when reordering data from gridcell-timestep to timestep-gridcell structure
	int OpenSpatial(char* name, bool replace_original_file=false);
	int OutputConvertedSpatial(char*);

// Methods to access private data
	/// Returns 1 if data for a gridcell is found in the input file, 0 if not.
	int FindCoord(Coord c) const {return FindRecord(c);}
	/// Returns the number of data columns
	int GetnColumns() const {return nColumns;}
	/// Copies the data column names to a string array
	int GetHeader(char cropnames[][MAXNAMESIZE]) const;
	/// Copies the whole header to a string
	int GetHeaderFull(char *header_line) const;
	/// Returns a pointer to a data name string for a column by its index
	char* GetHeader(int record) const;
	/// Returns the current coordinates
	Coord& GetCoord() {return currentStand;}
	/// Returns the number of gridcells with data in the input file
	int GetNCells();	// Calls ParseNCells() if nCells not yet set
	/// Returns the number of years in the input data
	int GetnYears() const {return nYears;}
	/// Returns the first year in the input data
	int GetFirstyear() const {return firstyear;}
	/// Returns the data format (EMPTY, GLOBAL_STATIC, GLOBAL_YEARLY, LOCAL_STATIC, LOCAL_YEARLY)
	int GetFormat() const {return format;}
	/// Returns true if data for requested coordinates are found, false if not.
	bool isloaded();
	/// Sets spacial resolution
	void SetSpacialResolution(double resolution) {spatial_resolution = resolution;}
	/// Returns spacial resolution
	double GetSpacialResolution() const {return spatial_resolution;}

// Functions for finding out if data columns contain sensible data for a specified gridlist 
#ifdef GUESS_VERSION
	/// Checks if data column has any values > 0 in any of the gridcells in the gridlist
	void CheckIfPresent(ListArray_id<Coord>& gridlist, double offset = 0.0);
#endif
	// Returns true if data column has any values > 0 in any of the gridcells in the gridlist (after CheckIfPresent() call)
	bool item_has_data(char* name);
	// Returns true if data name is in header
	bool item_in_header(char* name);

	/// Used by TimeDataDmem class to set pointer to full data copy
	void register_memory_copy(TimeDataDmem* mem_copy) {memory_copy = mem_copy;}
};

/// Class for loading all data for a gridlist to memory.
class TimeDataDmem {

	// PRIVATE VARIABLES

	/// Pointer to gridlist
	Coord *gridlist;
	/// Pointer to data array
	double **data;
	/// Number of data columns
	int nColumns;
	/// Number of data years
	int nYears;
	/// Number of data gridcells
	int nCells;
	// Spacial resolution of input data in degrees
	double spatial_resolution;
	/// Whether the input file structure includes a header line with column names and coordinates on each line of data
	bool ifheader;
	/// String array of data column names
	char header_arr[MAXRECORDS][MAXNAMESIZE];
	/// Index of current gridcell in data array
	int currentCell;
	/// First year of data
	int firstyear;
	/// Whether data for the requested coordinates have been found and loaded
	bool loaded;

	// PRIVATE METHODS

	/// Converts calender year to valid year position in data array.
	int CalenderYearToPosition(int calender_year) const;
	/// Sets coord at index position
	void SetCoord(int index, Coord c);
	/// Sets data at index position
	void SetData(int index, double* data);

public:

	// PUBLIC METHODS

	/// Constructor
	TimeDataDmem();
	/// Deconstructor
	~TimeDataDmem();

	/// Allocates memory. Returns 0 if error
	void Open(int nCells, int nColumns, int nYears);
	/// Releases dynamically allocated memory.
	void Close();
#if defined GUESS_VERSION
	/// Copies all data for gridlist to memory
	void CopyFromTimeDataD(TimeDataD& Data, ListArray_id<Coord>& gridlistX);
#endif
	/// Loads data for a certain coordinate. Returns 0 if coordinate not found.
	int Load(Coord c);
	/// Returns a single data value for a certain year and data column
	double Get(int calender_year, int column) const;		// Returns a single value.
	/// Returns a single data value for column with header string name. Returns -999 if name not found.
	double Get(int calender_year, const char* name) const;

	/// Returns the first year in the input data
	int GetFirstyear() {return firstyear;}
	/// Returns true if data for requested coordinates are found, false if not.
	bool isloaded() const { return loaded;}
	/// Sets spacial resolution
	void SetSpacialResolution(double resolution) {spatial_resolution = resolution;}
};

} // namespace InData

#endif//INDATA_H
