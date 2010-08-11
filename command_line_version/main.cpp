///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Main module for command line version of LPJ-GUESS
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
// Header file name:      main.h
// Source code file name: main.cpp
// Written by:            Ben Smith
// Version dated:         2002-09-22

#include "config.h"
#include "guess.h"
#include <stdarg.h>

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
	// (comment out one or both of these if output to the screen and/or log file
	// is not required)

	fprintf(stdout,"%s",(char*)output);
	fprintf(logfile,"%s",(char*)output);
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

int main(int argc,char* argv[]) {

	// Open log file if possible
	// (comment out if log file output not desired - you will have to comment
	// out the corresponding fprintf's in functions dprintf and fail above also)

	logfile=fopen(file_log,"wt");
	if (!logfile) {
		printf("main: could not open log file %s for output",(char*)file_log);
		exit(99);
	}

	// Call the framework
	framework(argc,argv);

	// Say goodbye
	dprintf("\nFinished\n");

	return 0;
}
