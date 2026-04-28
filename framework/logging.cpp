#include <stdarg.h>
#include <stdexcept>

#include "logging.h"
#include "shell.h"
#include "guess.h"

#define write_message(level, fmt) \
    do { \
        va_list args; \
        va_start(args, fmt); \
        vlog(level, fmt, args); \
        va_end(args); \
    } while (0)

const char *PFX_FATAL = "FATAL";
const char *PFX_ERROR = "ERROR";
const char *PFX_WARNING = " WARN";
const char *PFX_INFORMATION = " INFO";
const char *PFX_DIAGNOSTIC = " DIAG";
const char *PFX_DEBUG = " DEBUG";
const char *PFX_TRACE = " TRACE";

void get_pfx(LogLevel level, const char **pfx) {
    switch (level) {
        case LogLevel::Fatal:
            *pfx = PFX_FATAL;
            break;
        case LogLevel::Error:
            *pfx = PFX_ERROR;
            break;
        case LogLevel::Warning:
            *pfx = PFX_WARNING;
            break;
        case LogLevel::Information:
            *pfx = PFX_INFORMATION;
            break;
        case LogLevel::Diagnostic:
            *pfx = PFX_DIAGNOSTIC;
            break;
        case LogLevel::Debug:
            *pfx = PFX_DEBUG;
            break;
        case LogLevel::Trace:
            *pfx = PFX_TRACE;
            break;
        default:
            throw std::runtime_error("Invalid log level");
    }
}

bool should_log(LogLevel level) {
    // fixme
    return (int)level <= verbosity;
}

[[ noreturn ]]
void vlog_fatal(const char *fmt, va_list args) {
    fprintf(stderr, "\n\nFATAL ERROR:\n\n");
    vfprintf(stderr, fmt, args);
    exit(99);
}

void vlog(LogLevel level, const char *fmt, va_list args) {
    if (!should_log(level)) {
        return;
    }

    const char *pfx;
    get_pfx(level, &pfx);
    lprintf("[");
    if (date.is_initialised()) {
        lprintf("%04d-%02d-%02d ", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1);
    }
    lprintf("%s] ", pfx);

    vlprintf(fmt, args);
    lprintf("\n");
}

void log(LogLevel level, const char *fmt, ...) {
    write_message(level, fmt);

    if (level == LogLevel::Fatal) {
        va_list args;
        va_start(args, fmt);
        vlog_fatal(fmt, args);
        va_end(args); // Unreachable
    }
}

void log_fatal(const char *fmt, ...) {
    write_message(LogLevel::Fatal, fmt);

    va_list args;
    va_start(args, fmt);
    vlog_fatal(fmt, args);
    va_end(args); // Unreachable
}

void log_error(const char *fmt, ...) {
    write_message(LogLevel::Error, fmt);
}

void log_warning(const char *fmt, ...) {
    write_message(LogLevel::Warning, fmt);
}

void log_info(const char *fmt, ...) {
    write_message(LogLevel::Information, fmt);
}

void log_diagnostic(const char *fmt, ...) {
    write_message(LogLevel::Diagnostic, fmt);
}

void log_debug(const char *fmt, ...) {
    write_message(LogLevel::Debug, fmt);
}

void log_trace(const char *fmt, ...) {
    write_message(LogLevel::Trace, fmt);
}
