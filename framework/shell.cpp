///////////////////////////////////////////////////////////////////////////////////////
/// \file shell.cpp
/// \brief The "shell" is the model's interface to the world
///
/// \author Joe Siltberg
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "shell.h"
#include "stacktrace.h"
#include "guess.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include <parallel.h>
#ifdef __unix__
#include <csignal>
#endif
#ifdef _MSC_VER
#include <windows.h>
#endif

namespace {

/// The global Shell object
std::unique_ptr<Shell> current_shell;

}

void dprintf(xtring format,...) {
	// printf-style function accessible throughout the model code.
	// Sends text to currently chosen Shell object, which in turn
	// sends it to the screen and/or log file.

	va_list v;
	va_start(v,format);

	xtring output;
	formatf(output,format,v);

	current_shell->log_message(output);

	va_end(v);
}

void lprintf(xtring format, ...) {
	va_list v;
	va_start(v, format);

	vlprintf(format, v);
	va_end(v);
}

void vlprintf(xtring format, va_list args) {
	// Prefix message with a timestamp.
	current_shell->vwrite_log((const char *)format, args);
}

[[noreturn]]
void fail(xtring format,...) {
	// printf-style function accessible throughout the model code.
	// Sends text to currently chosen Shell, which prints it and
	// terminates the program.

	va_list v;
	va_start(v,format);

	xtring output;
	formatf(output,format,v);
	va_end(v);

	current_shell->fail(output);

	// This is never hit, and Shell::fail() is [[noreturn]], but some compilers
	// complain nonetheless.
	std::abort();
}


// Do the last things before exiting, e.g. termiante MPI
void cleanupAtExit() {
	current_shell->cleanupAtExit();
}


// Handle console control events, e.g. a ctrl-C signal, sighup.
// These handlers are registred in the CommandLineShell constructor.
#ifdef __unix__
void sig_handler(int signum) {

	// All linux signals that stops the program execution.
	// SIGHUP and SIGINT will result in the desired outcome
	// (i.e. calling cleanupAtExit()) while SIGQUIT, SIGTERM, and
	// SIGKILL aborts the program too abruptly for it to happen,
	// at least on simba2.
	switch (signum) {
	case SIGHUP:	
	case SIGINT:	
	case SIGQUIT:	
	case SIGTERM:	
	case SIGKILL:
		cleanupAtExit();		// Do the last things before exiting.
		exit(99);
		break;
	}
}
#endif
#ifdef _MSC_VER
BOOL CtrlHandler(DWORD ctrlType) {

	// The windows signals are: CTRL_C_EVENT (ctrl-c), CTRL_BREAK_EVENT (similar), CTRL_CLOSE_EVENT,
	// CTRL_LOGOFF_EVENT, CTRL_SHUTDOWN_EVENT. I.e. all of the possible signals should unload MPI before 
	// termination of the program. So no need for a switch (ctrlType). Just do the same action for all of them.
	// N.B. This is not used by the Graphical shell: see instead abort_request_received() (dllmain.cpp).

	// Do the last things before exiting
	cleanupAtExit();	

	// Return FALSE to let the system go on with its action in response to the signal,
	// i.e. to terminate LPJ-GUESS as per the users' ctrl-C request or other ctrl-signal,
	// e.g. let LPJ-GUESS graphical shell do the termination rather than doing it here.
	// (TRUE would mean that our code would handle the signal instead of the system.)
	return FALSE;		
}
#endif


void plot(xtring window_name,xtring series_name,double x,double y) {
	current_shell->plot(window_name, series_name, x, y);
}


void resetwindow(xtring window_name) {
	current_shell->resetwindow(window_name);
}


void clear_all_graphs() {
	current_shell->clear_all_graphs();
}

void open3d() {
	current_shell->open3d();
}

void plot3d() {
	current_shell->plot3d();
}

bool abort_request_received() {
	return current_shell->abort_request_received();
}

void set_shell(Shell* s) {
	current_shell = std::unique_ptr<Shell>(s);
}

void plot3d_fileopen() {
	current_shell->plot3d_fileopen();
}

void plot3d_fileclose() {
	current_shell->plot3d_fileclose();
}

FILE* plot3d_getfilehandle() {
	return current_shell->plot3d_getfilehandle();
}

////////////////////////////////////////////////////////////////////////////////
//// Implementation of CommandLineShell
////////////////////////////////////////////////////////////////////////////////

CommandLineShell::CommandLineShell(const char* logfile_path) {
	logfile = fopen(logfile_path, "wt");

	if (!logfile) {
		printf("Could not open log file %s for output\n", logfile_path);
		exit(99);
	}

	// Register the SIGHUP, SIGINT or ctrl-C signal handler
#ifdef __unix__
	signal(SIGHUP, sig_handler);
	signal(SIGINT, sig_handler);
#endif
#ifdef _MSC_VER
	if (!SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE)) {
		printf("Could not register the ctrl-C handler.\n");
	}
#endif
}

CommandLineShell::~CommandLineShell() {
	if (logfile) {
		fclose(logfile);
	}
}

[[noreturn]]
void CommandLineShell::fail(const char* message) {
	// The leading newline finishes an in-place progress message ending in '\r'.
	log_message("\n\nFATAL ERROR:\n\n");
	log_message(message);
	log_message("\n\nSimulation terminated.\n");
	write_stack_trace(logfile);
	cleanupAtExit();
	exit(99);
}

// Do the last things before exiting, e.g. termiante MPI
void CommandLineShell::cleanupAtExit() {

	// Terminate parallel communication if available.
	GuessParallel::finalize();
}

void CommandLineShell::log_message(const char* message) {
	fprintf(stdout,"%s", message);
	fflush(stdout);
	write_log(message);
}

void CommandLineShell::write_log(const char *fmt, ...) {
	va_list v;
	va_start(v, fmt);

	vwrite_log(fmt, v);

	va_end(v);
}

void CommandLineShell::vwrite_log(const char *fmt, va_list args) {
	vfprintf(logfile, fmt, args);
	fflush(logfile);
}

void CommandLineShell::plot(const char* window_name, 
                            const char* series_name, 
                            double x, 
                            double y) {
	// Can't do anything here	 
}

void CommandLineShell::open3d() {
	// Can't do anything here	 
}

void CommandLineShell::plot3d() {
	// Can't do anything here	 
}

void CommandLineShell::resetwindow(const char* window_name) {
	// Can't do anything here
}

void CommandLineShell::clear_all_graphs() {
	// Can't do anything here
}

bool CommandLineShell::abort_request_received() {
	// Can't do anything here

	return false;
}

void CommandLineShell::plot3d_fileopen() {
	// Can't do anything here
}

void CommandLineShell::plot3d_fileclose() {
	// Can't do anything here
}

FILE* CommandLineShell::plot3d_getfilehandle() { 
	// Can't do anything here

	return 0; 
}
