///////////////////////////////////////////////////////////////////////////////////////
/// \file stacktrace.h
/// \brief Platform-specific stack trace support for fatal errors
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_STACKTRACE_H
#define LPJ_GUESS_STACKTRACE_H

#include <stdio.h>

/// Write a stack trace to an already-open log file. This never writes to stdout.
void write_stack_trace(FILE* output);

#endif // LPJ_GUESS_STACKTRACE_H
