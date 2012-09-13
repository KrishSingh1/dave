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
: help(false) {

	if (!parse_arguments(argc, argv)) {
		print_usage(argv[0]);
	}
}

bool CommandLineArguments::parse_arguments(int argc, char** argv) {
	if (argc < 2) {
		return false;
	}

	// For now, just use argv[1], and require it to be either -help or the insfile
	if (tolower(argv[1]) == "-help") {
		help = true;
	}
	else if (argv[1][0] == '-') {
		fprintf(stderr, "Unknown option: %s\n", argv[1]);
		return false;
	}
	else {
		insfile = argv[1];
	}

	return true;
}

void CommandLineArguments::print_usage(const char* command_name) const {
	fprintf(stderr, "\nUsage: %s <instruction-script-filename> | -help", command_name);
	exit(EXIT_FAILURE);
}

bool CommandLineArguments::get_help() const {
	return help;
}

const char* CommandLineArguments::get_instruction_file() const {
	return insfile.c_str();
}
