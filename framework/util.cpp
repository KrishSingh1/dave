#include "util.h"

#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#ifdef _WIN32
#define mkdir(path, mode) _mkdir((path))
#endif

/* Make a directory; already existing dir okay */
void maybe_mkdir(const char* path, mode_t mode)
{
    errno = 0;

    /* Try to make the directory */
    if (mkdir(path, mode) == 0)
        return;

    /* If it fails for any reason but EEXIST, fail */
    if (errno != EEXIST) {
		char *reason = strerror(errno);
		errno = 0;
		throw_xtring("Unable to create directory %s: %s", path, reason);
	}

    /* Check if the existing path is a directory */
    struct stat st;
    if (stat(path, &st) != 0) {
		char *reason = strerror(errno);
		errno = 0;
		throw_xtring("Unable to stat %s: %s", path, reason);
	}

    /* If not, fail with ENOTDIR */
    if (!S_ISDIR(st.st_mode)) {
		std::string reason = "Path already exists but is not a directory";
        throw_xtring("Unable to create directory %s: %s", path, reason.c_str());
    }
}

void mkdir_p(const char *path)
{
    /*
	Adapted from
	https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
	*/
    mode_t mode = 0777;

    /* Copy string so it's mutable */
    char *_path = strdup(path);
	if (_path == NULL) {
		throw std::runtime_error("Out of memory");
	}

	try {
		// Iterate over the path.
		for (char *p = _path + 1; *p; p++) {
			if (*p == '/') {
				// Temporarily truncate the string.
				*p = '\0';

				// Attempt to create this directory (this can throw).
				maybe_mkdir(_path, mode);

				// Un-truncate the path.
				*p = '/';
			}
		}   

		// This is probably redundant.
		maybe_mkdir(_path, mode);

		free(_path);

	} catch (...) {
		free(_path);
		throw;
	}
}
