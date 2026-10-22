#ifndef LPJ_GUESS_LOGGING_H
#define LPJ_GUESS_LOGGING_H

#include <string>
#include <stdint.h>

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

/*
Resets the run statistics, including log message count, error count, and warning
count.
*/
void logging_reset_run_stats();

/*
Prints a summary of the run, including total log messages, errors, and warnings.
*/
void logging_print_run_summary();

#endif // LPJ_GUESS_LOGGING_H
