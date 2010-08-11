///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Main module for command line version of LPJ-GUESS
//                        Version adapted for parallel runs on Simba
//                        Version adapted for parallel runs on Monlith
// Header file name:      main.h
// Source code file name: main.cpp
// Written by:            Ben Smith, Lund University
// Version dated:         2005-01-21

#include "config.h"
#include "guess.h"
#include <stdarg.h>
#include <unistd.h>
#include <mpi.h>

#define WORKTAG 1

///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

FILE* logfile;
xtring logfile_submitdir; // BEN 2007-08-16


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

void dprintf_mon(xtring format,...) {

        // printf-style function accessible throughout the model code.
        // Sends text to log file on submit directory (for monitoring progress)
	// BEN 2007-08-16.

        va_list v;
        va_start(v,format);

        xtring output;
        formatf(output,format,v);

        // Produce output
        // No output to stdin - assume batch processing
        // (comment out one or both of these if output to the screen and/or log file
        // is not required)

	FILE* out=fopen(logfile_submitdir,"a");
	if (out) {
        	fprintf(out,"%s",(char*)output);
	}

	fclose(out);
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


///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run
// Expected arguments on milleotto (NOTE CHANGE!! Ben 2007-08-16):
// <insfile> <submit-dir> 
// --> submit-dir = directory from which the job was submitted (= PBS_O_WORKDIR)

int instno;

int main(int argc,char* argv[]) {

	// Parallel
	
	const int MAXPATH=1024;
	int myrank,ntasks,nfile,p;
	xtring workdir,maindir,insfile,submitdir;
	char* fargv[2];
	MPI_Status stat;
	
	MPI_Init(&argc,&argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&ntasks);

	instno=myrank+1;
		
	// Change working directory to "runXXX" 
	
	workdir.printf("run%d",instno);
	if (chdir(workdir))
		fail("main: instance %d could not go to directory %s",instno,(char*)workdir);

        // Get filename part of insfile (which might be a path)

	insfile=argv[1];
        p=insfile.find('/');
        while (p!=-1) {
                insfile=insfile.mid(p+1);
                p=insfile.find('/');
        }

	// BEN 2007-08-16
	// Get submit directory
	submitdir=argv[2];

	// Open log file if possible
	// (comment out if log file output not desired - you will have to comment
	// out the corresponding fprintf's in functions dprintf and fail above also)

	logfile=fopen(file_log,"wt");
	if (!logfile) {
		printf("main: could not open log file %s for output",(char*)file_log);
		exit(99);
	}

	// BEN 2007-08-16
	// Special log file on submit directory (for checking progress)

	if (submitdir[submitdir.len()-1]=='/')
		logfile_submitdir.printf("%srun%d/guess.log",(char*)submitdir,instno);
	else
		logfile_submitdir.printf("%s/run%d/guess.log",(char*)submitdir,instno);	

	FILE* out=fopen(logfile_submitdir,"w");
	if (!out) {
		printf("Could not open %s for output",(char*)logfile_submitdir);
		exit(99);
	}

	if (myrank==0) { // Master instance

		// Seed the slaves
		int i;
		//for (i=1;i<ntasks;i++)
			//MPI_Send(&work,1, MPI_INT,i,WORKTAG,MPI_COMM_WORLD);
	}
	
	printf("%s : Started process %d of %d\n",argv[0],instno,ntasks);

	// Call the framework

	fargv[0]=argv[0];
	fargv[1]=(char*)insfile;
	framework(2,fargv);
	
	// Wait for all processes to finish, then append output files
	// and write to main directory

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Say goodbye
	dprintf("\nFinished\n");

	MPI_Finalize();

	return 0;
}
