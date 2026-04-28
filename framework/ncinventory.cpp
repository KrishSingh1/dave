#ifdef HAVE_NETCDF

#include "shell.h"
#include "ncinventory.h"
#include "guessnc.h"
#include "logging.h"

static const char *ATTR_STD_NAME = "standard_name";
static const char *ATTR_UNITS = "units";

NCInventory::NCInventory(std::vector<std::string> files) {
    for (const std::string &file : files) {
        catalog.push_back(parse(file));
    }
}

std::unique_ptr<NCInventory::FileInfo> NCInventory::parse(const std::string &file) const {
    int ncid = GuessNC::open_ncdf((const char *)file.c_str());
    std::vector<std::unique_ptr<VariableMetadata>> variables;
    try
    {
        int nvars = GuessNC::get_nvars(ncid);
        for (int i = 0; i < nvars; i++) {
            auto var_info = get_variable_info(ncid, i);
            if (var_info) {
                variables.push_back(std::move(var_info));
            }
        }
        GuessNC::close_ncdf(ncid);
    }
    catch (const std::exception &e)
    {
        // TODO: RAII via a wrapper around ncid? Probably overkill...
        GuessNC::close_ncdf(ncid);
        throw;
    }

    return std::make_unique<FileInfo>(file, std::move(variables));
}

std::unique_ptr<NCInventory::VariableMetadata> NCInventory::get_variable_info(int ncid, int varid) const {
    std::string std_name;
    std::string units;
    std::string name = GuessNC::get_variable_name(ncid, varid);
    bool has_std_name = GuessNC::get_attribute(ncid, varid, ATTR_STD_NAME, std_name);
    bool has_units = GuessNC::get_attribute(ncid, varid, ATTR_UNITS, units);
    if (has_std_name && has_units) {
        return std::make_unique<VariableMetadata>(name, std_name, units);
    }
    // TODO: true logging
    log_warning("Ignoring variable '%s' in file '%s' because it doesn't have both a standard_name and units attribute\n",
            name.c_str(), GuessNC::get_variable_name(ncid, varid).c_str());
    return nullptr;
}

bool NCInventory::try_get_variable_name(const std::string &std_name, VariableInfo &info) const {
    bool found = false;
    for (const auto &file : catalog) {
        std::string variable_name;
        if (file->try_get_variable_name(std_name, variable_name)) {
            if (found) {
                throw std::runtime_error("Multiple variables found with standard name '" + std_name + "' across multiple files");
            }
            info.file = file->get_file_name();
            info.name = variable_name;
            found = true;
        }
    }
    return found;
}

VariableInfo NCInventory::get_variable_name(const std::string &std_name) const {
    VariableInfo binding;
    if (!try_get_variable_name(std_name, binding)) {
        throw std::runtime_error("No variable found with standard name '" + std_name + "'");
    }
    return binding;
}

VariableInfo NCInventory::get_variable_name(const std::string &std_name, const std::string &file_name) const {
    for (const auto &file : catalog) {
        if (file->get_file_name() == file_name) {
            std::string name;
            if (!file->try_get_variable_name(std_name, name)) {
                throw std::runtime_error("No variable found with standard name '" + std_name + "' in file '" + file_name + "'");
            }
            VariableInfo info;
            info.file = file_name;
            info.name = name;
            return info;
        }
    }
    throw std::runtime_error("File not found in inventory: '" + file_name + "'");
}

bool NCInventory::FileInfo::try_get_variable_name(const std::string &std_name, std::string &variable_name) const {
    bool found = false;
    for (const auto &var : variables) {
        if (var->get_std_name() == std_name) {
            if (found) {
                throw std::runtime_error("Multiple variables found with standard name '" + std_name + "' in file '" + file_name + "'");
            }
            variable_name = var->get_name();
            found = true;
        }
    }
    return found;
}

#endif // HAVE_NETCDF
