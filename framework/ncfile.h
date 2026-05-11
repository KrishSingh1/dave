#ifdef HAVE_NETCDF

#ifndef LPJ_GUESS_NCFILE_H
#define LPJ_GUESS_NCFILE_H

#include <string>

/*
Wrapper around a NetCDF file, responsible for opening and closing the file.

TODO:

- Expand API so that we don't need to expose the file ID? This would require
  wrapping most of the nc_* functions in GuessNC.
*/
class NCFile {
    // Allow NCVariable and NCInventory to access id(). By convention they
    // should avoid accessing id_ directly.
    friend class NCVariable;
    friend class NCInventory;
public:
    explicit NCFile(std::string file_path) : NCFile(file_path.c_str()) {}
    explicit NCFile(const char *file_path);
    ~NCFile() noexcept;

	// Delete copy constructors to prevent accidental copying of this class.
	NCFile(const NCFile&) = delete;
	NCFile& operator=(const NCFile&) = delete;
	NCFile(NCFile&&) = delete;            // optional, but safest
	NCFile& operator=(NCFile&&) = delete; // optional, but safest

private:
    int id() const { return id_; }
    int id_;
    std::string file_path;
};

#endif // LPJ_GUESS_NCFILE_H

#endif // HAVE_NETCDF
