///////////////////////////////////////////////////////////////////////////////////////
/// \file main.cpp
/// \brief Main module for command line version of LPJ-GUESS
///
/// \author Ben Smith
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "framework.h"
#include "parallel.h"

///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";



///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run

int main(int argc,char* argv[]) {

	// Set our shell for the model to communicate with the world
	set_shell(new CommandLineShell(file_log));

	// Initialize parallel communication if available
	GuessParallel::init(argc, argv);

	// Call the framework
	framework(argc,argv);

	// Say goodbye
	dprintf("\nFinished\n");

	return 0;
}
