#include <algorithm>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <glob.h>
#include <sys/stat.h>
#endif

#include "pathresolver.h"
#include "shell.h"

namespace {
std::string parent_path(const std::string& path) {
	size_t slash = path.find_last_of("/\\");
	if (slash == std::string::npos) {
		return "";
	}
	return path.substr(0, slash);
}

std::string join_path(const std::string& parent, const std::string& child) {
	if (parent.empty()) {
		return child;
	}
	char last = parent[parent.size() - 1];
	if (last == '/' || last == '\\') {
		return parent + child;
	}
	return parent + "/" + child;
}
}

std::vector<std::string> PathResolver::resolve(const std::string &path) const {
#ifdef _WIN32
    return resolve_windows(path);
#else
    return resolve_posix(path);
#endif
}

#ifdef _WIN32
std::vector<std::string> PathResolver::resolve_windows(const std::string &path) const {
    std::vector<std::string> resolved_paths;
	const std::string parent = parent_path(path);

    WIN32_FIND_DATAA find_file_data;
    HANDLE hFind = FindFirstFileA(path.c_str(), &find_file_data);
    if (hFind == INVALID_HANDLE_VALUE) {
		fail("Error: Unable to resolve path '%s'", path.c_str());
    }
    do {
        if (!(find_file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            resolved_paths.push_back(join_path(parent, find_file_data.cFileName));
        }
    } while (FindNextFile(hFind, &find_file_data) != 0);
    FindClose(hFind);

	if (resolved_paths.empty()) {
		fail("Error: Path does not match any files: '%s'", path.c_str());
	}

	std::sort(resolved_paths.begin(), resolved_paths.end());
    return resolved_paths;
}
#else
void PathResolver::handle_error_posix(const char *path, int error_code) const {
    // Handle errors based on the error code.
    switch (error_code) {
        case GLOB_NOSPACE:
            fail("Error: Out of memory while resolving path '%s'", path);
			break;
        case GLOB_ABORTED:
            fail("Error: Read error while resolving path '%s'", path);
			break;
        case GLOB_NOMATCH:
            fail("Error: Path does not match any files: '%s'", path);
			break;
        default:
            fail("Error: Unknown error code %d while resolving path '%s'", error_code, path);
			break;
    }
}

std::vector<std::string> PathResolver::resolve_posix(const std::string &path) const {
    std::vector<std::string> resolved_paths;
    int flags = GLOB_ERR;
    glob_t result{};
    int res = glob(path.c_str(), flags, nullptr, &result);
    if (res != 0) {
		globfree(&result);
        handle_error_posix(path.c_str(), res);
    }

    for (size_t i = 0; i < result.gl_pathc; i++) {
		struct stat st{};
		if (stat(result.gl_pathv[i], &st) == 0 && S_ISREG(st.st_mode)) {
			resolved_paths.push_back(result.gl_pathv[i]);
		}
    }
    globfree(&result);

	if (resolved_paths.empty()) {
		fail("Error: Path does not match any files: '%s'", path.c_str());
	}

	std::sort(resolved_paths.begin(), resolved_paths.end());
    return resolved_paths;
}
#endif // _WIN32
