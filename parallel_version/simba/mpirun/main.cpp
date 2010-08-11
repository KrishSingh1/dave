///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Main module for command line version of LPJ-GUESS
//                        Version adapted for parallel runs on Simba
// Header file name:      main.h
// Source code file name: main.cpp
// Written by:            Ben Smith, Lund University
// Version dated:         2005-01-26

#include "config.h"
#include "guess.h"
#include <stdarg.h>
#include <unistd.h>
#include <mpi.h>


///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

FILE* logfile;


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL FUNCTIONS
// These functions are declared in the framework header file and are (therefore)
// accessible throughout the model code

void fail(xtring format,...) {

	// printf-style function accessible throughout the model code.
	// Sends text to stdio (screen) and log file, then terminates program

	va_list v;
	va_start(v,format);

	xtring output;
	formatf(output,format,v);

	// Produce output
	// No output to stdout - assume batch processing
	// (comment out one or both of these if output to the screen and/or log file
	// is not required)

	fprintf(stdout,"%s\n",(char*)output);
	fprintf(logfile,"%s\n",(char*)output);
	
	exit(99);
}


void dprintf(xtring format,...) {

	// printf-style function accessible throughout the model code.
	// Sends text to stdio (screen) and log file.

	va_list v;
	va_start(v,format);
	
	xtring output;
	formatf(output,format,v);

	// Produce output
	// No output to stdin - assume batch processing
	// (comment out one or both of these if output to the screen and/or log file
	// is not required)

	fprintf(stdout,"%s",(char*)output);
	fprintf(logfile,"%s",(char*)output);
	fflush(logfile);
}


void plot(xtring window_name,xtring series_name,double x,double y) {

	// Can't do anything here
}


void resetwindow(xtring window_name) {

	// Can't do anything here
}


void clear_all_graphs() {

	// Can't do anything here
}


bool abort_request_received() {

	// Can't do anything here
	
	return false;
}


void splitgrid(xtring gridfile,int nprocess) {

        // Opens the specified gridlist file and splits it into <nprocess> roughly equal
        // segments, placing one in each run directory

        int gridno;
        int ngrid,p,fac,count;
        xtring outfile,line,dummy;
        FILE* out;

        FILE* in=fopen(gridfile,"rt");
        if (!in) {
                printf("splitgrid: could not open %s for input\n",(char*)gridfile);
                exit(99);
        }

        // Work out number of lines in gridfile (should be one for each lon/lat pair)

        ngrid=0;
        while (!feof(in)) {
                readfor(in,"A#",&dummy);
                if (!feof(in)) ngrid++;
        }

        fac=ngrid/nprocess;
        if (ngrid%nprocess) fac++;

        rewind(in);

        // Get filename part of gridfile (which might be a path)

        p=gridfile.find('/');
        while (p!=-1) {
                gridfile=gridfile.mid(p+1);
                p=gridfile.find('/');
        }

        for (gridno=0;gridno<nprocess;gridno++) {

                outfile.printf("./run%d/%s",gridno+1,(char*)gridfile);
                out=fopen(outfile,"wt");
                if (!out) {
                        printf("splitgrid: could not open %s for output\n",(char*)outfile);
                        exit(99);
                }

                count=0;

                while (!feof(in) && ((count%fac) || count==0 ||  gridno==nprocess-1)) {
                        readfor(in,"A#",&line);
                        if (!feof(in)) fprintf(out,"%s\n",(char*)line);
                        count++;
                }

                fclose(out);

                printf("Created %s\n",(char*)outfile);
        }
        fclose(in);
}


void append(int nfile,char* filename[],int nprocess) {

        // Extracts and appends output from process run directories to main directory
        // This will only work for: text files with or without a single header line,
        // labels line, or leading blank line
        // Other files you will have to write your own software and extract from each
        // individual run directory

        xtring text,outfile,infile;
        int i,j,pos;
        FILE* in,*out;
        bool firstline;

        for (i=0;i<nfile;i++) {
                outfile=filename[i];
                out=fopen(outfile,"wt");
                if (out) {
                        for (j=0;j<nprocess;j++) {
                                infile.printf("run%d/%s",j+1,filename[i]);
                                in=fopen(infile,"rt");
                                if (!in) {
                                        printf("append: could not open %s for input, aborting for %s\n",
						(char*)infile,(char*)outfile);
                                        goto abort;
                                }
                                printf("Reading data from %s\n",(char*)infile);
                                firstline=true;
                                while (!feof(in)) {
                                        readfor(in,"A#",&text);
                                        if (firstline && text!="") {
                                                pos=text.findnotoneof(" .-\t");
                                                if (pos)
                                                        if (text[pos]>='0' && text[pos]<='9') firstline=false;
                                        }
                                        if (firstline && j) readfor(in,"A#",&text);
                                        firstline=false;
                                        if (text!="")
                                                fprintf(out,"%s\n",(char*)text);
                                }
                                fclose(in);
                        }
                }
                else {
                        printf("append: warning - could not open %s for output\n",(char*)outfile);
                }

                fclose(out);
                printf("Written %s\n",(char*)outfile);
abort:
        }
}


///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run
// Expected arguments <insfile> <gridlistfile> <outfile1> <outfile2> ... <outfile_n>

int instno;

int main(int argc,char* argv[]) {

	// Parallel
	
	const int MAXPATH=1024;
	int myrank,ntasks,nfile,p;
	xtring workdir,maindir,insfile;
	char* fargv[2];
	MPI_Status stat;
	
	MPI_Init(&argc,&argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&ntasks);

	instno=myrank+1;
		
	// Record current directory which will be the main output directory
	
	if (myrank==0) { // Master instance only

		splitgrid(argv[2],ntasks);
		char buf[MAXPATH];
		if (getcwd(buf,MAXPATH))
			maindir=buf;
		else
			fail("main: could not get path to working directory (pathname may be too long)");
	}
	
	// Change working directory to "./runx" (x=instance)
	// Each instance runs and writes output to a different "run" directory

	// Wait here to make sure splitgrid finishes
	MPI_Barrier(MPI_COMM_WORLD);	

	workdir.printf("./run%d",instno);
	printf("Going to %s\n",(char*)workdir);

	if (chdir(workdir))
		fail("main: instance %d could not go to directory %s",instno,(char*)workdir);

        // Get filename part of insfile (which might be a path)

	insfile=argv[1];
        p=insfile.find('/');
        while (p!=-1) {
                insfile=insfile.mid(p+1);
                p=insfile.find('/');
        }

	// Open log file if possible
	// (comment out if log file output not desired - you will have to comment
	// out the corresponding fprintf's in functions dprintf and fail above also)

	logfile=fopen(file_log,"wt");
	if (!logfile) {
		printf("main: could not open log file %s for output",(char*)file_log);
		exit(99);
	}	

	printf("%s : Started process %d of %d\n",argv[0],instno,ntasks);

	// Call the framework

	fargv[0]=argv[0];
	fargv[1]=(char*)insfile;

	printf("Calling framework with insfile=%s\n",fargv[1]);
	framework(2,fargv);
	
	// Wait for all processes to finish, then append output files
	// and write to main directory

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (myrank==0) { // Master instance
		chdir(maindir);
		nfile=argc-3;
		printf("\n");
		append(nfile,&argv[3],ntasks);
	}

	// Say goodbye
	dprintf("\nFinished\n");

	MPI_Finalize();

	return 0;
}
