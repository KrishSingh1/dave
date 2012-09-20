///////////////////////////////////////////////////////////////////////////////////////
/// \file commandlinearguments.h
/// \brief Takes care of the command line arguments to LPJ-GUESS
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H
#define LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H

#include <string>

/// Parses and stores the user's command line arguments
class CommandLineArguments {
public:
	/// Send in the arguments from main() to this constructor
	/** If the arguments are malformed, this constructor will
	 *  print out usage information and exit the program.
	 */
	CommandLineArguments(int argc, char** argv);

	/// Returns the instruction filename
	const char* get_instruction_file() const;

	/// Returns true if the user has specified the help option
	bool get_help() const;

	/// Returns true if the user has specified the goto-rundir option
	bool get_goto_rundir() const;

private:
	/// Does the actual parsing of the arguments
	bool parse_arguments(int argc, char** argv);

	/// Prints out usage information and exits the program
	void print_usage(const char* command_name) const;

	/// Instruction filename specified on the command line
	std::string insfile;

	/// Whether the user wants help on how to run the LPJ-GUESS command
	bool help;

	/// Whether we should step into a run directory before starting
	bool goto_rundir;
};

#endif // LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H
