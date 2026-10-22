#pragma once

#include <string>

class TemporaryDirectory {
public:
    TemporaryDirectory();
    ~TemporaryDirectory() noexcept;

    // Delete copy constructors
	TemporaryDirectory(const TemporaryDirectory&) = delete;
	TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
	TemporaryDirectory(TemporaryDirectory&&) = delete;            // optional, but safest
	TemporaryDirectory& operator=(TemporaryDirectory&&) = delete; // optional, but safest

    /*
    Return the absolute path to an instruction file in the fixtures directory.
    The returned file path is guaranteed to exist, and will contain absolute
    paths internally, so may be used/parsed from anywhere.

    @param basename The base name of the instruction file, e.g.
                    "cfinput_test.ins"
    @return The absolute path to an equivalent instruction file.
    */
    std::string ins_file(const std::string &basename);

    /*
    Get the absolute path to the temporary directory. The directory is
    guaranteed to exist, and will be removed when the TemporaryDirectory object
    is destroyed.

    @return The absolute path to the temporary directory.
    */
    std::string get_path() const;
private:
    std::string tmp_path;
};
