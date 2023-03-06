#ifndef _UTIL_H_
#define _UTIL_H_

#include "gutil.h"
#include <stdexcept>

#ifndef PATH_MAX
// E.g. msvc
#define PATH_MAX 4096
#endif

/*
Throw an exception using the specified format string and arguments.

@param fmt: printf-style format string.
@param ...: Arguments to the format string.
*/
#define throw_xtring(fmt, ...) \
	xtring tmp_msg; \
	tmp_msg.printf((fmt), __VA_ARGS__); \
	throw std::runtime_error((const char *)tmp_msg);

/*
Throw a nested exception using the specified format string and arguments.

@param fmt: printf-style format string.
@param ...: Arguments to the format string.
*/
#define throw_xtring_nested(fmt, ...) \
	xtring tmp_msg; \
	tmp_msg.printf((fmt), __VA_ARGS__); \
	std::throw_with_nested(std::runtime_error((const char *)tmp_msg));

/*
Attempt to create the specified directory tree (equivalent to mkdir -p). Will
throw upon failure, but "directory already exists" is considered success.

@param path: The directory tree to be created.
*/
void mkdir_p(const char *path);

/*
Return current time (in μs since unix epoch).
*/
double get_wall_time();

#endif // _UTIL_H_
