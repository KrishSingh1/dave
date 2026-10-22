#pragma once

#include <string>
#include "temporary_directory.h"

/*
Given a fixture .ins file basename, copy it to the temporary directory, insert
absolute paths to other dependent files, and return the absolute path to the
modified .ins file.

@param temp_directory: The temporary directory to which the .ins file will be
                       copied.
@param basename: The base name of the fixture .ins file, e.g.
                 "cfinput_test.ins"
@return: The absolute path to the modified .ins file.
*/
std::string ins_file(const TemporaryDirectory &temp_directory,
                     const std::string &basename);

/*
Run a simulation using the specified fixture .ins file basename. The .ins file
will be copied to a temporary directory, and absolute paths to other dependent
files will be inserted.

@param temp_directory: The temporary directory to which the .ins file will be
                       copied.
@param ins_basename: The base name of the fixture .ins file, e.g.
                      "cfinput_test.ins"
*/
void run_simulation(const TemporaryDirectory &temp_directory,
                    const std::string &ins_basename);

/*
Get the path to an output file in the specified output directory within the
temporary directory.

@param temp_directory: The temporary directory under which the simulation was
                       run.
@param out_dir: The output directory in which outputs were written. Typically
                this is the value of the outputdirectory parameter.
@param filename: The basename of the output file (e.g. "cpool.out").

@return: The absolute path to the output file.
*/
std::string out_file_path(const TemporaryDirectory &temp_directory,
                                 const std::string &out_dir,
                                 const std::string &filename);
