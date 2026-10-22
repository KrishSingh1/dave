#include <stdarg.h>
#include <stdexcept>
#include <stdint.h>

#include "logging.h"
#include "parallel.h"
#include "log_scope.h"
#include "shell.h"
#include "guess.h"
#include "util.h"
#include "stopwatch.h"

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

uint64_t log_message_count = 0;
uint64_t error_count = 0;
uint64_t warning_count = 0;
Stopwatch stopwatch;

void logging_reset_run_stats() {
    log_message_count = 0;
    error_count = 0;
    warning_count = 0;
    stopwatch.reset();
}

void logging_get_run_stats(uint64_t *out_log_message_count,
    uint64_t *out_error_count, uint64_t *out_warning_count,
    double *out_elapsed_time, double *out_cpu_time) {
    *out_log_message_count = log_message_count;
    *out_error_count = error_count;
    *out_warning_count = warning_count;
    *out_elapsed_time = stopwatch.get_elapsed_walltime();
    *out_cpu_time = stopwatch.get_elapsed_cputime();
}

void logging_print_run_summary() {
    double elapsed_time = stopwatch.get_elapsed_walltime();
    double cpu_time = stopwatch.get_elapsed_cputime();

    log_info("Wall time: %.2f seconds", elapsed_time);
    log_info("CPU time: %.2f seconds", cpu_time);

    printf("\n\nSimulation completed in %.2f seconds with %lu errors and %lu warnings.\n",
             elapsed_time, error_count, warning_count);
    printf("See guess.log for full details\n");
}

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

void write_pfx(LogLevel level);
void write_log_scope();

[[ noreturn ]]
void vlog_fatal(const char *fmt, va_list args) {
    write_pfx(LogLevel::Fatal);
    write_log_scope();

    xtring format = fmt;
    xtring message;
    va_list args_copy;
    va_copy(args_copy, args);
    formatf(message, format, args_copy);
    va_end(args_copy);
    fail("%s", (const char *)message);
}

void write_log_scope() {
    if (log_context.depth == ContextDepth::None ||
        (log_context.depth == ContextDepth::Gridcell && !GuessParallel::parallel)) {
        return;
    }

    lprintf("[");

    // We could always include gridcell coordinates in the log message, but in
    // practice it adds a lot of unnecessary noise. This is mainly useful in MPI
    // runs, in order to disambiguate output from different nodes.
    bool wrote_gridcell = false;
    if (GuessParallel::parallel && log_context.depth >= ContextDepth::Gridcell) {
        wrote_gridcell = true;
        lprintf("lon=%.4f, lat=%.4f", log_context.lon, log_context.lat);
    }

    if (log_context.depth >= ContextDepth::Stand) {
        lprintf("%ss%d", wrote_gridcell ? " " : "", log_context.stand_id);
    }

    if (log_context.depth >= ContextDepth::Patch) {
        lprintf(" p%d", log_context.patch_id);
    }

    if (log_context.depth >= ContextDepth::Individual) {
        lprintf(" i%d", log_context.individual_id);
    }

    lprintf("] ");
}

void write_pfx(LogLevel level) {
    const char *pfx;
    get_pfx(level, &pfx);
    lprintf("[");
    if (date.is_initialised()) {
        lprintf("%04d-%02d-%02d ", date.get_calendar_year(), date.month + 1, date.dayofmonth + 1);
    }
    lprintf("%s] ", pfx);
}

void vlog(LogLevel level, const char *fmt, va_list args) {
    if (!should_log(level)) {
        return;
    }

    if (level == LogLevel::Error) {
        error_count++;
    } else if (level == LogLevel::Warning) {
        warning_count++;
    }

    write_pfx(level);
    write_log_scope();

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
