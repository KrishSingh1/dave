#ifndef LPJ_GUESS_LOGGING_H
#define LPJ_GUESS_LOGGING_H

#include <string>

/*
Defines the level of a log message.
*/
enum class LogLevel {
    Fatal,
    Error,
    Warning,
    Information,
    Diagnostic,
    Debug,
    Trace
};

/*
Write a formatted log message.

@param level: The level of the log message.
@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log(LogLevel level, const char *fmt, ...);

/*
Log a fatal error. This will cause the simulation to terminate.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
[[ noreturn ]]
void log_fatal(const char *fmt, ...);

/*
Log a recoverable error. Simulation execution will not terminate.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_error(const char *fmt, ...);

/*
Log a warning message.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_warning(const char *fmt, ...);

/*
Log an informational message.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_info(const char *fmt, ...);

/*
Log a diagnostic message.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_diagnostic(const char *fmt, ...);

/*
Log a debug message.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_debug(const char *fmt, ...);

/*
Log a trace message.

@param fmt: A printf-style format string.
@param ...: Arguments corresponding to the format string.
*/
void log_trace(const char *fmt, ...);

#endif // LPJ_GUESS_LOGGING_H
