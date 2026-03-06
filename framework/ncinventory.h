#ifndef LPJ_GUESS_NCINVENTORY_H
#define LPJ_GUESS_NCINVENTORY_H

#ifdef HAVE_NETCDF

#include <memory>
#include <string>
#include <vector>

/*
Encapsulates the name of a variable in a specific netcdf file.
*/
struct VariableInfo {
    // Path to the file.
    std::string file;
    // Name of the variable.
    std::string name;
};

/*
NCInventory is responsible for mapping standard names to variable names among
a collection of NetCDF files.
*/
class NCInventory {
public:
    /*
    Initialise the inventory with a list of NetCDF files. This will populate the
    inventory with the variables found in these files.
    */
    NCInventory(std::vector<std::string> files);

    /*
    Get the variable name corresponding to a standard name. If multiple
    variables (across all files) have the same standard name, an exception will
    be thrown. If no variable has the given standard name, an exception will be
    thrown.

    @param std_name: The standard name to look up.
    @return: The variable name corresponding to the given standard name, and the
    file in which it was found.
    */
    VariableInfo get_variable_name(const std::string &std_name) const;

    /*
    Get the variable name corresponding to a standard name from the specified
    file. If multiple variables in the file have the same standard name, an
    exception will be thrown. If no variable in the file has the given standard
    name, an exception will be thrown.

    @param std_name: The standard name to look up.
    @param file_name: The name of the file to look in.
    @return: The variable name corresponding to the given standard name in the
    specified file, and the file in which it was found.
    */
    VariableInfo get_variable_name(const std::string &std_name, const std::string &file_name) const;

    /*
    Try to get the variable name corresponding to a standard name. If multiple
    variables have the same standard name, an exception will be thrown. If no
    variable has the given standard name, false will be returned and
    variable_name will be left unchanged.

    @param std_name: The standard name to look up.
    @param binding: (out parameter) The variable name corresponding to the given
    standard name, if found, and the file in which it was found.
    @return: True if a variable with the given standard name was found, false
    otherwise.
    */
    bool try_get_variable_name(const std::string &std_name, VariableInfo &variable_name) const;
private:
    /*
    Encapsulates metadata about a variable in a NetCDF file.
    */
    class VariableMetadata {
public:
        VariableMetadata(std::string name, std::string std_name, std::string units) :
            name(name), std_name(std_name), units(units) {}

        const std::string &get_name() const { return name; }
        const std::string &get_std_name() const { return std_name; }
        const std::string &get_units() const { return units; }
private:
        std::string name;
        std::string std_name;
        std::string units;
    };

    /*
    Encapsulates metadata about a NetCDF file, including its variables.
    */
    class FileInfo {
public:
        FileInfo(std::string file_name, std::vector<std::unique_ptr<VariableMetadata>> variables) :
            file_name(file_name), variables(std::move(variables)) {}
        const std::string &get_file_name() const { return file_name; }
        const std::vector<std::unique_ptr<VariableMetadata>> &get_variables() const { return variables; }
        bool try_get_variable_name(const std::string &std_name, std::string &variable_name) const;
private:
        std::string file_name;
        std::vector<std::unique_ptr<VariableMetadata>> variables;
    };

    /*
    Parse a NetCDF file and return a FileInfo object containing metadata about
    the file and its variables.

    @param file: The path to the NetCDF file to parse.
    @return: A FileInfo object containing metadata about the file and its variables.
    */
    std::unique_ptr<FileInfo> parse(const std::string &file) const;

    /*
    Parse a variable in a NetCDF file and return a VariableInfo object
    containing metadata about the variable.

    @param ncid: The NetCDF file ID.
    @param varid: The variable ID.
    @return: A VariableInfo object containing metadata about the variable.
    */
    std::unique_ptr<VariableMetadata> get_variable_info(int ncid, int varid) const;

    /*
    The parsed inventory of files and variables.
    */
    std::vector<std::unique_ptr<FileInfo>> catalog;
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCINVENTORY_H
