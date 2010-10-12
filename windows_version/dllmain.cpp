///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Main module for interface to Windows shell
// Header file name:      main.h
// Source code file name: main.cpp
// Written by:            Ben Smith
// Version dated:         2001-09-05
//
// The shell should call function dll_main, passing a GuessParam object containing run
// time argument(s) for GUESS and pointers to the executable's own callback functions.

#include "config.h"
#include "dllmain.h"

#include <process.h>
#include <stdarg.h>

///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";


///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

FILE* logfile=false;
bool waiting;
bool ifabort;
xtring* poutput;
PlotArgs* pplotargs;
MessagePrintString* message_print_string;
MessagePlot* message_plot;
MessageFinished* message_finished;
MessageResetwindow* message_resetwindow;
MessageClearGraphs* message_clear_graphs;


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL FUNCTIONS
// These functions are declared in the framework header file and are (therefore)
// accessible throughout the model code

void fail(xtring format,...) {

	// printf-style function accessible throughout the model code.
	// Sends text to Windows shell and log file, then terminates program

	xtring output;
	va_list v;
	va_start(v,format);
	formatf(output,format,v);

	xtring* pbuf=new xtring;
	*pbuf=output;
	*pbuf+="\n";
	message_print_string(pbuf);

	if (logfile) fprintf(logfile,"%s\n",(char*)output);
	
	message_finished();

	_endthread();
}


void dprintf(xtring format,...) {

	// printf-style function accessible throughout the model code.
	// Sends text to Windows shell and log file.

	xtring output;
	va_list v;
	va_start(v,format);
	formatf(output,format,v);

	xtring* pbuf=new xtring;
	*pbuf=output;
	message_print_string(pbuf);

	if (logfile) fprintf(logfile,"%s",(char*)output);
}


void plot(xtring window_name,xtring series_name,double x,double y) {

	PlotArgs* pplotargs=new PlotArgs;

	pplotargs->window_name=window_name;
	pplotargs->series_name=series_name;
	pplotargs->x=x;
	pplotargs->y=y;
	pplotargs->rescale=true;
	
	message_plot(pplotargs);
}


void resetwindow(xtring window_name) {

	xtring* pxtring=new xtring;
	*pxtring=window_name;
	message_resetwindow(pxtring);
}


void clear_all_graphs() {

	waiting=true;
	message_clear_graphs();
}


bool abort_request_received() {

	// May be called by framework to respond to abort request from Windows shell
	// (returns true if shell has sent an abort request, otherwise false)

	return ifabort;
}


__declspec(dllexport) void cleanup_print_string(xtring* pxtring) {

	// To be called by shell to deallocate memory after GUESS sends a
	// message_print_string message

	delete pxtring;
}


__declspec(dllexport) void cleanup_plot(PlotArgs* pplotargs) {

	// To be called by shell to deallocate memory after GUESS sends a
	// message_plot message

	delete pplotargs;
}


///////////////////////////////////////////////////////////////////////////////////////
// DLL_MAIN
// This is the function called by the Windows shell to run the model

__declspec(dllexport) int dll_main(GuessParam param) {

	// Store parameters sent from shell as file scope global variables

	poutput=param.poutput;
	pplotargs=param.pplotargs;
	message_print_string=param.message_print_string;
	message_plot=param.message_plot;
	message_finished=param.message_finished;
	message_resetwindow=param.message_resetwindow;
	message_clear_graphs=param.message_clear_graphs;

	ifabort=false;

	// Open log file if possible
	if (!logfile) logfile=fopen(file_log,"wt");

	// Call the framework
	framework(param.argc,param.argv);

	// Say goodbye
	message_finished();

	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////
// ABORT_RUN
// Called by windows shell to request abort of current run

__declspec(dllexport) void abort_run() {


	dprintf("\nAbort request received\n");
	ifabort=true;
}

