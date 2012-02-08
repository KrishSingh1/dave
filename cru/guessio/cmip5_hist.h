//////////////////////////////////////////////////////////////////////////////////////
// CMIP5_HIST.H
// Header file for input from a fast data archive
// Created automatically by FastArchive on Fri Dec 02 21:10:02 2011
//
// The following #includes should appear in your source code file:
//
//   #include <stdio.h>
//   #include <stdlib.h>
//   #include <string.h>
//   #include "C://Data//Climate_Data//ipsl_cm5a_lr_hist_r1i1p1//cmip5_hist.h"
//
// Functionality to retrieve data from the archive is provided by class Cmip5_histArchive.
// The following public functions are provided:
//
// bool open(char* filename)
//   Attempts to open the specified file as a fast data archive. The format must be
//   exactly compatible with this version of cmip5_hist.h (normally the archive and
//   header file should have been produced together by the same program using class
//   CFastArchive). Returns false if the file could not be opened or had format
//   errors. open() with no argument is equivalent to open("cmip5_hist.bin").
//
// void close()
//   Closes the archive (if open).
//
// bool rewind()
//   Sets the file pointer to the first record in the archive file. Returns false if
//   no archive file is currently open.
//
// bool getnext(Cmip5_hist& obj)
//   Retrieves the next record in the archive file and advances the file pointer to
//   the next record. Data are written to the member variables of obj. Returns false if
//   no archive file is currently open or if the file pointer is beyond the last
//   record. Use rewind() and getnext() to retrieve data sequentially from the archive.
//
// bool getindex(Cmip5_hist& obj)
//   Searches the archive for a record matching the values specified for the index
//   items (lon and lat) in obj. If a matching record is found, the data are
//   written to the member variables of obj. Returns true if the archive was open and
//   a matching record was found, otherwise false. The search is iterative and fast.
//
// Sample program:
//
//   Cmip5_histArchive ark;
//   Cmip5_hist data;
//   bool success,flag;
//
//   // Retrieve all records in sequence and print values of lon and lat:
//
//   success=ark.open("cmip5_hist.bin");
//   if (success) {
//      flag=ark.rewind();
//      while (flag) {
//         flag=ark.getnext(data);
//         if (flag)
//            printf("Loaded record: lon=%g, lat=%g\n",data.lon,data.lat);
//      }
//   }
//   
//   // Look for a record with lon=-1800, lat=-900:
//
//   data.lon=-1800;
//   data.lat=-900;
//   success=ark.getindex(data);
//   if (success) printf("Found it!\n");
//   else printf("Not found\n");
//
//   ark.close();


struct Cmip5_hist {

	// Index part

	double lon;
	double lat;

	// Data part

	double mtemp[1872];
	double mprec[1872];
	double mswrad[1872];
};


const long CMIP5_HIST_NRECORD=59191;
const int CMIP5_HIST_DATA_LENGTH=14040;
const int CMIP5_HIST_INDEX_LENGTH=7;
const int CMIP5_HIST_HEADERSIZE=548;
unsigned char CMIP5_HIST_HEADER[CMIP5_HIST_HEADERSIZE-4]={
	0x01,0x02,0x24,0x00,0x00,0x00,0x07,0x00,0x02,0x04,0x6C,0x6F,0x6E,0x00,0x2D,0x31,0x38,0x30,0x30,0x00,
	0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0xC8,0x09,0x36,0x00,0x24,0x26,
	0x40,0x00,0xBC,0x0D,0x42,0x00,0x31,0x38,0x30,0x30,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,
	0x00,0x00,0x04,0x00,0x00,0x00,0xC8,0x09,0x36,0x00,0x24,0x26,0x40,0x00,0xBC,0x0D,0x42,0x00,0x31,0x00,
	0x30,0x30,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0xC8,0x09,
	0x36,0x00,0x24,0x26,0x40,0x00,0xBC,0x0D,0x42,0x00,0x0C,0x04,0x6C,0x61,0x74,0x00,0x2D,0x39,0x30,0x30,
	0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0x58,0x0A,0x36,0x00,
	0x24,0x26,0x40,0x00,0xB8,0x0D,0x42,0x00,0x39,0x30,0x30,0x00,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,
	0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0x58,0x0A,0x36,0x00,0x24,0x26,0x40,0x00,0xB8,0x0D,0x42,0x00,
	0x31,0x00,0x30,0x00,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,
	0x58,0x0A,0x36,0x00,0x24,0x26,0x40,0x00,0xB8,0x0D,0x42,0x00,0x0B,0x00,0x00,0x36,0xD8,0x00,0x03,0x06,
	0x6D,0x74,0x65,0x6D,0x70,0x00,0x2D,0x31,0x30,0x30,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,
	0x00,0x00,0x06,0x00,0x00,0x00,0xA0,0x0A,0x36,0x00,0x3F,0x27,0x40,0x00,0xB2,0x0D,0x42,0x00,0x31,0x30,
	0x30,0x30,0x30,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xA0,0x0A,
	0x36,0x00,0x3F,0x27,0x40,0x00,0xB2,0x0D,0x42,0x00,0x30,0x2E,0x30,0x31,0x00,0x00,0x00,0x00,0xB0,0xF5,
	0x18,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xA0,0x0A,0x36,0x00,0x3F,0x27,0x40,0x00,0xB2,0x0D,
	0x42,0x00,0x14,0x00,0x00,0x07,0x50,0x06,0x6D,0x70,0x72,0x65,0x63,0x00,0x2D,0x31,0x30,0x30,0x00,0x00,
	0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0x48,0x0B,0x36,0x00,0x3F,0x27,
	0x40,0x00,0xAA,0x0D,0x42,0x00,0x31,0x30,0x30,0x30,0x30,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,
	0x00,0x00,0x06,0x00,0x00,0x00,0x48,0x0B,0x36,0x00,0x3F,0x27,0x40,0x00,0xAA,0x0D,0x42,0x00,0x30,0x2E,
	0x30,0x31,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0x48,0x0B,
	0x36,0x00,0x3F,0x27,0x40,0x00,0xAA,0x0D,0x42,0x00,0x14,0x00,0x00,0x07,0x50,0x07,0x6D,0x73,0x77,0x72,
	0x61,0x64,0x00,0x2D,0x31,0x30,0x30,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x07,
	0x00,0x00,0x00,0x90,0x0B,0x36,0x00,0x3F,0x27,0x40,0x00,0xA3,0x0D,0x42,0x00,0x31,0x30,0x30,0x30,0x30,
	0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,0x00,0x00,0x00,0x07,0x00,0x00,0x00,0x90,0x0B,0x36,0x00,0x3F,
	0x27,0x40,0x00,0xA3,0x0D,0x42,0x00,0x30,0x2E,0x30,0x31,0x00,0x00,0x00,0x00,0xB0,0xF5,0x18,0x00,0x00,
	0x00,0x00,0x00,0x07,0x00,0x00,0x00,0x90,0x0B,0x36,0x00,0x3F,0x27,0x40,0x00,0xA3,0x0D,0x42,0x00,0x14,
	0x00,0x00,0x07,0x50};


class Cmip5_histArchive {

private:

	FILE* pfile;
	long recno;
	long datano;
	unsigned char pindex[CMIP5_HIST_INDEX_LENGTH];
	unsigned char pdata[CMIP5_HIST_DATA_LENGTH];
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

		fseek(pfile,CMIP5_HIST_INDEX_LENGTH*(n-recno),SEEK_CUR);
		fread(pindex,CMIP5_HIST_INDEX_LENGTH,1,pfile);
		datano=pindex[CMIP5_HIST_INDEX_LENGTH-4]*0x1000000+pindex[CMIP5_HIST_INDEX_LENGTH-3]*0x10000+
			pindex[CMIP5_HIST_INDEX_LENGTH-2]*0x100+pindex[CMIP5_HIST_INDEX_LENGTH-1];
		recno=n+1;
		iseof=(recno==CMIP5_HIST_NRECORD);
	}

	void getdata() {

		fseek(pfile,CMIP5_HIST_INDEX_LENGTH*-recno+(datano-CMIP5_HIST_NRECORD)*CMIP5_HIST_DATA_LENGTH,SEEK_CUR);
		fread(pdata,CMIP5_HIST_DATA_LENGTH,1,pfile);
		fseek(pfile,CMIP5_HIST_DATA_LENGTH*(CMIP5_HIST_NRECORD-datano-1)+CMIP5_HIST_INDEX_LENGTH*recno,SEEK_CUR);
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

		pheader=new unsigned char[CMIP5_HIST_HEADERSIZE-4];
		if (!pheader) {
			printf("Out of memory\n");
			fclose(pfile);
			pfile=NULL;
			return false;
		}
		::rewind(pfile);
		fread(pheader,CMIP5_HIST_HEADERSIZE-4,1,pfile);
		for (i=0;i<CMIP5_HIST_HEADERSIZE-4;i++) {
			/*if (pheader[i]!=CMIP5_HIST_HEADER[i]) {
				printf("Format of %s incompatible with this version of cmip5_hist.h\n",filename);
				fclose(pfile);
				pfile=NULL;
				delete pheader;
				return false;
			}*/ //AA CMIP5 now the same header file can be used for all archives.
		}
		delete pheader;

		::rewind(pfile);
		fseek(pfile,CMIP5_HIST_HEADERSIZE+CMIP5_HIST_DATA_LENGTH*CMIP5_HIST_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

public:

	Cmip5_histArchive() {
		pfile=NULL;
	}

	~Cmip5_histArchive() {
		if (pfile) fclose(pfile);
	}

	bool open(char* filename) {
		return initialise(filename);
	}

	bool open() {
		return open("cmip5_hist.bin");
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
		fseek(pfile,CMIP5_HIST_HEADERSIZE+CMIP5_HIST_DATA_LENGTH*CMIP5_HIST_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

	bool getnext(Cmip5_hist& obj) {

		if (!pfile || iseof) return false;

		int i;

		getindex(recno);
		getdata();

		obj.lat=popreal(pindex,3,11,1,-900);
		obj.lon=popreal(pindex,3,12,1,-1800);

		for (i=1871;i>=0;i--) obj.mswrad[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);
		for (i=1871;i>=0;i--) obj.mprec[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);
		for (i=1871;i>=0;i--) obj.mtemp[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);

		return true;
	}

	bool getindex(Cmip5_hist& obj) {

		if (!CMIP5_HIST_NRECORD || !pfile) return false;

		// else

		unsigned char ptarget[3]={0,0,0};
		unsigned char buf[4];
		bitify(buf,obj.lon,-1800,1);
		merge(ptarget,buf,12);
		bitify(buf,obj.lat,-900,1);
		merge(ptarget,buf,11);

		int i,c;
		bool not_done=true;
		long direction=1;
		long offset=CMIP5_HIST_NRECORD/2+CMIP5_HIST_NRECORD%2;
		long thisrecord=0-CMIP5_HIST_NRECORD%2;

		while (not_done) {
			thisrecord+=offset*direction;
			if (thisrecord>=CMIP5_HIST_NRECORD) thisrecord=CMIP5_HIST_NRECORD-1;
			else if (thisrecord<0) thisrecord=0;
			if (offset==1) not_done=false;

			getindex(thisrecord);
			getdata();

			c=compare_index(pindex,ptarget);
			if (c<0) direction=1;
			else if (c>0) direction=-1;
			else { // found

				for (i=1871;i>=0;i--) obj.mswrad[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);
				for (i=1871;i>=0;i--) obj.mprec[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);
				for (i=1871;i>=0;i--) obj.mtemp[i]=popreal(pdata,CMIP5_HIST_DATA_LENGTH,20,0.01,-100);

				return true;
			}
			offset=offset/2+offset%2;
		}

		return false;
	}
};
