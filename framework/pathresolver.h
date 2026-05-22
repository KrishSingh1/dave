#ifndef LPJ_GUESS_PATH_RESOLVER_H
#define LPJ_GUESS_PATH_RESOLVER_H

#include <string>
#include <vector>

/*
 * PathResolver is responsible for resolving file paths from shell globs. It
 * takes a base path and resolves it to a list of file paths that match the
 * given pattern.
 *
 * Supported constructs in the input path include:
 * - Wildcards: '*' to match any sequence of characters, '?' to match any single character
 */
class PathResolver {
    public:
        /*
        Resolve the given path to a list of file paths.

        @param path: The input path to resolve. This can be a single file path
                     or a pattern with wildcards.
        @return: A vector of resolved file paths that match the input path.
        */
        std::vector<std::string> resolve(const std::string &path) const;
    private:
#ifdef _WIN32
        /*
        Resolve the given path on Windows systems.

        @param path: The input path to resolve. This can be a single file path
                     or a pattern with wildcards.
        @return: A vector of resolved file paths that match the input path.
        */
        std::vector<std::string> resolve_windows(const std::string &path) const;
#else
        /*
        Resolve the given path on POSIX-compliant systems.

        @param path: The input path to resolve. This can be a single file path
                     or a pattern with wildcards.
        @return: A vector of resolved file paths that match the input path.
        */
        std::vector<std::string> resolve_posix(const std::string &path) const;

        /*
        Handle errors that occur during path resolution on POSIX-compliant systems.
        @param path: The input path that was being resolved when the error occurred.
        @param error_code: The error code returned by the glob function.
        */
        void handle_error_posix(const char *path, int error_code) const;
#endif // _WIN32
};

#endif // LPJ_GUESS_PATH_RESOLVER_H
