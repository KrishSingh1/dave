///////////////////////////////////////////////////////////////////////////////////////
/// \file commandlinearguments.cpp
/// \brief Takes care of the command line arguments to LPJ-GUESS
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "commandlinearguments.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>

namespace {

std::string tolower(const char* str) {
	std::string result(str);
	std::transform(result.begin(), result.end(), result.begin(), ::tolower);
	return result;
}

}

CommandLineArguments::CommandLineArguments(int argc, char** argv) 
: help(false),
  goto_rundir(false) {

	if (!parse_arguments(argc, argv)) {
		print_usage(argv[0]);
	}
}

bool CommandLineArguments::parse_arguments(int argc, char** argv) {
	if (argc < 2) {
		return false;
	}

	// For now, just consider anything starting with '-' to be an option,
	// and anything else to be the ins file
	for (int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			std::string option = tolower(argv[i]);
			if (option == "-help") {
				help = true;
			}
			else if (option == "-goto-rundir") {
				goto_rundir = true;
			}
			else {
				fprintf(stderr, "Unknown option: \"%s\"\n", argv[i]);
				return false;
			}
		}
		else {
			if (insfile.empty()) {
				insfile = argv[i];
			}
			else {
				fprintf(stderr, "Two arguments parsed as ins file: %s and %s\n",
						  insfile.c_str(), argv[i]);
				return false;
			}
		}
	}

	return true;
}

void CommandLineArguments::print_usage(const char* command_name) const {
	fprintf(stderr, "\nUsage: %s [-goto-rundir] <instruction-script-filename> | -help", 
			  command_name);
	exit(EXIT_FAILURE);
}

bool CommandLineArguments::get_help() const {
	return help;
}

bool CommandLineArguments::get_goto_rundir() const {
	return goto_rundir;
}

const char* CommandLineArguments::get_instruction_file() const {
	return insfile.c_str();
}
