#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <cerrno>
#include <cstring>

#include "commandlinearguments.h"
#include "config.hpp"
#include <filesystem>
#include "framework.h"
#include "guess.h" // terrible
#include "temporary_directory.h"

std::string read_file(const std::string &path) {
    errno = 0;
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open file '" + path + "': " +
                                 std::strerror(errno));
    }

    std::ostringstream contents;
    contents << in.rdbuf();

    if (!in.good() && !in.eof()) {
        throw std::runtime_error("Failed to read file '" + path + "': " +
                                 std::strerror(errno));
    }

    return contents.str();
}

void write_file(const std::string &path, const std::string &contents) {
    errno = 0;
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file '" + path + "' for writing: " +
                                 std::strerror(errno));
    }

    out << contents;

    if (!out.good()) {
        throw std::runtime_error("Failed to write file '" + path + "': " +
                                 std::strerror(errno));
    }
}

bool valid_placeholder_name(const std::string &name) {
    if (name.empty()) {
        return false;
    }

    for (char c : name) {
        if (!(std::isupper(c) || std::isdigit(c) || c == '_')) {
            return false;
        }
    }

    return true;
}

std::string render_template(const std::string &input,
                            const std::unordered_map<std::string, std::string> &replacements) {
    std::string output;
    output.reserve(input.size());

    std::size_t pos = 0;

    while (pos < input.size()) {
        // Find the next opening delimiter.
        const size_t open = input.find("{{", pos);

        if (open == std::string::npos) {
            // No more placeholders.
            output.append(input.substr(pos));
            break;
        }

        // Copy literal text before the placeholder.
        output.append(input.substr(pos, open - pos));

        const std::size_t name_begin = open + 2;
        const std::size_t close = input.find("}}", name_begin);

        if (close == std::string::npos) {
            throw std::runtime_error("Unterminated template placeholder at "
                                     "offset: " + std::to_string(open));
        }

        std::size_t name_len = close - name_begin;
        const std::string name = input.substr(name_begin, name_len);

        if (!valid_placeholder_name(name)) {
            throw std::runtime_error("Invalid template placeholder '{{" +
                                     name + "}}' at offset: " +
                                     std::to_string(open));
        }

        auto it = replacements.find(name);

        if (it == replacements.end()) {
            throw std::runtime_error("Invalid template placeholder '{{" +
                                     name + "}}' at offset: " +
                                     std::to_string(open) + " (no replacement" +
                                     " provided)");
        }

        output.append(it->second);

        // Continue after the closing delimiter.
        pos = close + 2;
    }

    return output;
}

std::string ins_file(const TemporaryDirectory &temp_directory,
                     const std::string &basename) {
    std::string input = (std::filesystem::path(config::source_root) /
                            "tests" /
                            "test_insfiles" /
                            basename).string();

    std::string contents = read_file(input);
    contents = render_template(contents, {
        {"SOURCE_ROOT", config::source_root},
        {"WORK_DIR", temp_directory.get_path()}
    });

    std::string temp_path = (std::filesystem::path(temp_directory.get_path()) / basename).string();
    write_file(temp_path, contents);
    return temp_path;
}

static char GUESS[] = "guess";
static char INPUT[] = "-input";
static char CF[] = "cf";

void run_simulation(const TemporaryDirectory &temp_directory,
                    const std::string &ins_basename) {
    std::string ins_path = ins_file(temp_directory, ins_basename);

    const int NARGS = 4;

    // Use ins_restart.data() after C++17.
    char *args[NARGS] = { GUESS, INPUT, CF, &ins_path[0] };

    // This is a global variable. We should remove all pfts because it could
    // interfere with other tests.
    pftlist.killall();

    framework(CommandLineArguments(NARGS, args));
}

std::string out_file_path(const TemporaryDirectory &temp_directory,
                                 const std::string &out_dir,
                                 const std::string &filename) {
    std::filesystem::path path = temp_directory.get_path();
    path = path / out_dir / filename;
    return path.string();
}
