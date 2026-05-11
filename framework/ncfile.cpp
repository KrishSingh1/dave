#ifdef HAVE_NETCDF

#include "guessnc.h"
#include "ncfile.h"
#include "logging.h"

NCFile::NCFile(const char *file_path) {
    id_ = GuessNC::open_ncdf(file_path);
    this->file_path = std::string(file_path);
}

NCFile::~NCFile() noexcept {
    try {
        GuessNC::close_ncdf(id_);
    } catch (const std::exception &e) {
        // Log the error but don't throw from the destructor.
        log_error("Failed to close NetCDF file %s with ID %d: %s",
                  file_path.c_str(), id_, e.what());
    } catch (...) {
        log_error("Failed to close NetCDF file %s with ID %d: unknown error",
                  file_path.c_str(), id_);
    }
}

#endif // HAVE_NETCDF
