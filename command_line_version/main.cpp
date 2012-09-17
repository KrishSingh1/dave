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
#include "guessio.h"
#include "commandlinearguments.h"
#include "parallel.h"
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////////////////
// LOG FILE
// The name of the log file to which output from all dprintf and fail calls is sent is
// set here

xtring file_log="guess.log";



///////////////////////////////////////////////////////////////////////////////////////
// MAIN
// This is the function called when the executable is run

int main(int argc,char* argv[]) {
	// Initialize parallel communication if available.
	// This needs to be done before command line parsing since some MPI
	// implementations put their own arguments in our command line
	// (which should then be removed after the MPI initialization).
	GuessParallel::init(argc, argv);

	// Parse command line arguments
	CommandLineArguments args(argc, argv);

	// Change working directory according to rank if requested
	if (args.get_goto_rundir()) {
		xtring path;
		path.printf("./run%d", GuessParallel::get_rank()+1);

		if (change_directory(path) != 0) {
			fprintf(stderr, "Failed to change to run directory\n");
			return EXIT_FAILURE;
		}
	}

	// Set our shell for the model to communicate with the world
	set_shell(new CommandLineShell(file_log));

	if (args.get_help()) {
		printhelp();
		return EXIT_SUCCESS;
	}

	// Call the framework
	framework(args);

	// Say goodbye
	dprintf("\nFinished\n");

	return EXIT_SUCCESS;
}
