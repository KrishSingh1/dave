#include <filesystem>

#include "temporary_directory.h"
#include "util.h"

std::string create_temp_directory(std::string_view prefix) {
    const size_t MAX_ATTEMPTS = 65536;

    std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
    temp_dir = std::filesystem::absolute(temp_dir);
    std::string pfx = std::string(prefix) + "_";

    std::string dir_name;
    std::filesystem::path result;

    for (size_t i = 0; i < MAX_ATTEMPTS; i++) {
        dir_name = pfx + std::to_string(std::rand());
        result = temp_dir / dir_name;
        if (std::filesystem::create_directory(result)) {
            return result.string();
        }
    }

    throw std::runtime_error("Failed to create temporary directory after " +
                             std::to_string(MAX_ATTEMPTS) + " attempts");
}

TemporaryDirectory::TemporaryDirectory() {
    tmp_path = create_temp_directory("lpj-guess-tests");
}

TemporaryDirectory::~TemporaryDirectory() noexcept {
    try {
        std::filesystem::remove_all(tmp_path);
    } catch (const std::exception& e) {
        // Suppress exceptions in destructor.
        // Log the error message to stderr.
        std::fprintf(stderr, "Error removing temporary directory '%s': %s\n",
                     tmp_path.c_str(), (const char *)get_nested_error_msg(e));
    }
}

std::string TemporaryDirectory::get_path() const {
    return tmp_path;
}
