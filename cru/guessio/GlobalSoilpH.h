//////////////////////////////////////////////////////////////////////////////////////
// GLOBALSOILPH.H
// Header file for input from a fast data archive
// Created automatically by FastArchive on Mon Dec 12 14:59:19 2011
//
// The following #includes should appear in your source code file:
//
//   #include <stdio.h>
//   #include <stdlib.h>
//   #include <string.h>
//   #include "GlobalSoilpH.h"
//
// Functionality to retrieve data from the archive is provided by class GlobalSoilpHArchive.
// The following public functions are provided:
//
// bool open(char* filename)
//   Attempts to open the specified file as a fast data archive. The format must be
//   exactly compatible with this version of GlobalSoilpH.h (normally the archive and
//   header file should have been produced together by the same program using class
//   CFastArchive). Returns false if the file could not be opened or had format
//   errors. open() with no argument is equivalent to open("GlobalSoilpH.bin").
//
// void close()
//   Closes the archive (if open).
//
// bool rewind()
//   Sets the file pointer to the first record in the archive file. Returns false if
//   no archive file is currently open.
//
// bool getnext(GlobalSoilpH& obj)
//   Retrieves the next record in the archive file and advances the file pointer to
//   the next record. Data are written to the member variables of obj. Returns false if
//   no archive file is currently open or if the file pointer is beyond the last
//   record. Use rewind() and getnext() to retrieve data sequentially from the archive.
//
// bool getindex(GlobalSoilpH& obj)
//   Searches the archive for a record matching the values specified for the index
//   items (longitude and latitude) in obj. If a matching record is found, the data are
//   written to the member variables of obj. Returns true if the archive was open and
//   a matching record was found, otherwise false. The search is iterative and fast.
//
// Sample program:
//
//   GlobalSoilpHArchive ark;
//   GlobalSoilpH data;
//   bool success,flag;
//
//   // Retrieve all records in sequence and print values of longitude and latitude:
//
//   success=ark.open("GlobalSoilpH.bin");
//   if (success) {
//      flag=ark.rewind();
//      while (flag) {
//         flag=ark.getnext(data);
//         if (flag)
//            printf("Loaded record: longitude=%g, latitude=%g\n",data.longitude,data.latitude);
//      }
//   }
//   
//   // Look for a record with longitude=-180, latitude=-90:
//
//   data.longitude=-180;
//   data.latitude=-90;
//   success=ark.getindex(data);
//   if (success) printf("Found it!\n");
//   else printf("Not found\n");
//
//   ark.close();


struct GlobalSoilpH {

	// Index part

	double longitude;
	double latitude;

	// Data part

	double pH[2];
};


const long GLOBALSOILPH_NRECORD=64165;
const int GLOBALSOILPH_DATA_LENGTH=5;
const int GLOBALSOILPH_INDEX_LENGTH=7;
const int GLOBALSOILPH_HEADERSIZE=339;
unsigned char GLOBALSOILPH_HEADER[GLOBALSOILPH_HEADERSIZE-4]={
	0x01,0x01,0x53,0x00,0x00,0x00,0x07,0x00,0x02,0x0A,0x6C,0x6F,0x6E,0x67,0x69,0x74,0x75,0x64,0x65,0x00,
	0x2D,0x31,0x38,0x30,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x31,0x38,0x30,0x00,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0x30,0x2E,0x35,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x0A,0x09,0x6C,0x61,
	0x74,0x69,0x74,0x75,0x64,0x65,0x00,0x2D,0x39,0x30,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x39,
	0x30,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x30,0x2E,0x35,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0x09,0x00,0x00,0x00,0x05,0x00,0x01,0x03,0x70,0x48,0x00,0x2D,0x39,0x39,0x39,0x39,0x00,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x31,0x30,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x30,0x2E,
	0x30,0x31,0x00,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,
	0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0xCC,0x14,0x00,0x00,0x00,0x02};


class GlobalSoilpHArchive {

private:

	FILE* pfile;
	long recno;
	long datano;
	unsigned char pindex[GLOBALSOILPH_INDEX_LENGTH];
	unsigned char pdata[GLOBALSOILPH_DATA_LENGTH];
	bool iseof;

	long readbin(int nbyte) {

		unsigned char buf[4];
		long mult[4]={0x1000000,0x10000,0x100,1},val;
		int i;

		fread(buf,nbyte,1,pfile);

		val=0;
		for (i=0;i<nbyte;i++) {
			val+=buf[i]*mult[4-nbyte+i];
		}

		return val;
	}

	void getindex(long n) {

		fseek(pfile,GLOBALSOILPH_INDEX_LENGTH*(n-recno),SEEK_CUR);
		fread(pindex,GLOBALSOILPH_INDEX_LENGTH,1,pfile);
		datano=pindex[GLOBALSOILPH_INDEX_LENGTH-4]*0x1000000+pindex[GLOBALSOILPH_INDEX_LENGTH-3]*0x10000+
			pindex[GLOBALSOILPH_INDEX_LENGTH-2]*0x100+pindex[GLOBALSOILPH_INDEX_LENGTH-1];
		recno=n+1;
		iseof=(recno==GLOBALSOILPH_NRECORD);
	}

	void getdata() {

		fseek(pfile,GLOBALSOILPH_INDEX_LENGTH*-recno+(datano-GLOBALSOILPH_NRECORD)*GLOBALSOILPH_DATA_LENGTH,SEEK_CUR);
		fread(pdata,GLOBALSOILPH_DATA_LENGTH,1,pfile);
		fseek(pfile,GLOBALSOILPH_DATA_LENGTH*(GLOBALSOILPH_NRECORD-datano-1)+GLOBALSOILPH_INDEX_LENGTH*recno,SEEK_CUR);
	}

	double popreal(unsigned char* bits,int nbyte,int nbit,double scalar,double offset) {

		unsigned char buf;
		int nb=nbit/8,i;
		double rval=0.0;
		long mult[4]={1,0x100,0x10000,0x1000000};

		for (i=0;i<4;i++) {
			if (i<nb) rval+=bits[nbyte-i-1]*mult[i];
			else if (i==nb) {
				buf=bits[nbyte-i-1]<<(8-nbit%8);
				buf>>=8-nbit%8;
				rval+=buf*mult[i];
			}
		}

		for (i=nbyte-1;i>=0;i--) {
			if (i>=nb)
				bits[i]=bits[i-nb];
			else
				bits[i]=0;
		}

		nb=nbit%8;

		for (i=nbyte-1;i>=0;i--) {
			bits[i]>>=nb;
			if (i>0) {
				buf=bits[i-1];
				buf<<=8-nb;
				bits[i]|=buf;
			}
		}

		rval=rval*scalar+offset;

		return rval;
	}

	void bitify(unsigned char buf[4],double fval,double offset,double scalar) {

		long ival=(fval-offset)/scalar+0.5;
		buf[0]=ival/0x1000000;
		ival-=buf[0]*0x1000000;
		buf[1]=ival/0x10000;
		ival-=buf[1]*0x10000;
		buf[2]=ival/0x100;
		ival-=buf[2]*0x100;
		buf[3]=ival;
	}

	void merge(unsigned char ptarget[3],unsigned char buf[4],int bits) {

		int nb=bits/8;
		int i,j;
		unsigned char nib;
		for (i=0;i<3;i++) {

			if (i<3-nb)
				ptarget[i]=ptarget[i+nb];
			else
				ptarget[i]=0;
		}
		nb=bits%8;
		for (i=0;i<3;i++) {
			ptarget[i]<<=nb;
			if (i<3-1) {
				nib=ptarget[i+1]>>(8-nb);
				ptarget[i]|=nib;
			}
		}

		nb=bits/8;
		if (bits%8) nb++;
		for (i=1;i<=nb;i++)
			ptarget[3-i]|=buf[4-i];
	}

	int compare_index(unsigned char* a,unsigned char* b) {

		int i;
		for (i=0;i<3;i++) {
			if (a[i]<b[i]) return -1;
			else if (a[i]>b[i]) return +1;
		}

		return 0;
	}

	bool initialise(char* filename) {

		int i;
		unsigned char* pheader;

		if (pfile) fclose(pfile);
		pfile=fopen(filename,"rb");
		if (!pfile) {
			printf("Could not open %s for input\n",filename);
			return false;
		}

		pheader=new unsigned char[GLOBALSOILPH_HEADERSIZE-4];
		if (!pheader) {
			printf("Out of memory\n");
			fclose(pfile);
			pfile=NULL;
			return false;
		}
		::rewind(pfile);
		fread(pheader,GLOBALSOILPH_HEADERSIZE-4,1,pfile);
		for (i=0;i<GLOBALSOILPH_HEADERSIZE-4;i++) {
			if (pheader[i]!=GLOBALSOILPH_HEADER[i]) {
				printf("Format of %s incompatible with this version of GlobalSoilpH.h\n",filename);
				fclose(pfile);
				pfile=NULL;
				delete pheader;
				return false;
			}
		}
		delete pheader;

		::rewind(pfile);
		fseek(pfile,GLOBALSOILPH_HEADERSIZE+GLOBALSOILPH_DATA_LENGTH*GLOBALSOILPH_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

public:

	GlobalSoilpHArchive() {
		pfile=NULL;
	}

	~GlobalSoilpHArchive() {
		if (pfile) fclose(pfile);
	}

	bool open(char* filename) {
		return initialise(filename);
	}

	bool open() {
		return open("GlobalSoilpH.bin");
	}

	void close() {
		if (pfile) {
			fclose(pfile);
			pfile=NULL;
		}
	}

	bool rewind() {

		if (!pfile) return false;

		::rewind(pfile);
		fseek(pfile,GLOBALSOILPH_HEADERSIZE+GLOBALSOILPH_DATA_LENGTH*GLOBALSOILPH_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

	bool getnext(GlobalSoilpH& obj) {

		if (!pfile || iseof) return false;

		int i;

		getindex(recno);
		getdata();

		obj.latitude=popreal(pindex,3,9,0.5,-90);
		obj.longitude=popreal(pindex,3,10,0.5,-180);

		for (i=1;i>=0;i--) obj.pH[i]=popreal(pdata,GLOBALSOILPH_DATA_LENGTH,20,0.01,-9999);

		return true;
	}

	bool getindex(GlobalSoilpH& obj) {

		if (!GLOBALSOILPH_NRECORD || !pfile) return false;

		// else

		unsigned char ptarget[3]={0,0,0};
		unsigned char buf[4];
		bitify(buf,obj.longitude,-180,0.5);
		merge(ptarget,buf,10);
		bitify(buf,obj.latitude,-90,0.5);
		merge(ptarget,buf,9);

		int i,c;
		bool not_done=true;
		long direction=1;
		long offset=GLOBALSOILPH_NRECORD/2+GLOBALSOILPH_NRECORD%2;
		long thisrecord=0-GLOBALSOILPH_NRECORD%2;

		while (not_done) {
			thisrecord+=offset*direction;
			if (thisrecord>=GLOBALSOILPH_NRECORD) thisrecord=GLOBALSOILPH_NRECORD-1;
			else if (thisrecord<0) thisrecord=0;
			if (offset==1) not_done=false;

			getindex(thisrecord);
			getdata();

			c=compare_index(pindex,ptarget);
			if (c<0) direction=1;
			else if (c>0) direction=-1;
			else { // found

				for (i=1;i>=0;i--) obj.pH[i]=popreal(pdata,GLOBALSOILPH_DATA_LENGTH,20,0.01,-9999);

				return true;
			}
			offset=offset/2+offset%2;
		}

		return false;
	}
};
