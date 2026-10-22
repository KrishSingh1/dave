#ifndef _UTIL_H_
#define _UTIL_H_

#include "gutil.h"
#include <exception>
#include <stdexcept>

#ifndef PATH_MAX
#ifdef MAX_PATH
// E.g. msvc
#define PATH_MAX MAX_PATH
#else
// Fallback
#define PATH_MAX 4096
#endif
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
Format an exception and any nested exceptions as an indented message chain.
*/
xtring get_nested_error_msg(const std::exception& e);

#endif // _UTIL_H_
