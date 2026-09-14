///////////////////////////////////////////////////////////////////////////////////////
/// \file stacktrace.cpp
/// \brief Platform-specific stack trace support for fatal errors
///////////////////////////////////////////////////////////////////////////////////////

#include "stacktrace.h"

#if defined(HAVE_EXECINFO_BACKTRACE)

#include <execinfo.h>
#include <stdlib.h>
#include <string>

#if defined(__GNUG__) || defined(__clang__)
#include <cxxabi.h>
#endif

namespace {

std::string demangle_symbol(const char* symbol) {
	std::string result(symbol);
	std::string::size_type begin = result.find('(');
	std::string::size_type end = result.find('+', begin);
	if (begin == std::string::npos || end == std::string::npos || begin + 1 >= end) {
		return result;
	}

#if defined(__GNUG__) || defined(__clang__)
	std::string mangled = result.substr(begin + 1, end - begin - 1);
	int status = 0;
	char* demangled = abi::__cxa_demangle(mangled.c_str(), 0, 0, &status);
	if (status == 0 && demangled) {
		result.replace(begin + 1, end - begin - 1, demangled);
	}
	free(demangled);
#endif

	return result;
}

}

void write_stack_trace(FILE* output) {
	void* frames[64];
	int count = backtrace(frames, 64);

	fprintf(output, "\nStack trace:\n");
	if (count <= 2) {
		fprintf(output, "  No stack frames available.\n");
		fflush(output);
		return;
	}

	char** symbols = backtrace_symbols(frames, count);
	if (!symbols) {
		fprintf(output, "  Unable to resolve stack frame symbols.\n");
		fflush(output);
		return;
	}

	// Omit stack-trace helper frames.
	for (int i = 2; i < count; ++i) {
		fprintf(output, "  %s\n", demangle_symbol(symbols[i]).c_str());
	}
	free(symbols);
	fflush(output);
}

#elif defined(_WIN32)

#include <windows.h>

#ifdef HAVE_DBGHELP
#include <dbghelp.h>
#endif

void write_stack_trace(FILE* output) {
	void* frames[64];
	USHORT count = CaptureStackBackTrace(2, 64, frames, 0);

	fprintf(output, "\nStack trace:\n");
	if (count == 0) {
		fprintf(output, "  No stack frames available.\n");
		fflush(output);
		return;
	}

#ifdef HAVE_DBGHELP
	HANDLE process = GetCurrentProcess();
	bool symbols_available = SymInitialize(process, 0, TRUE) != FALSE;
	if (!symbols_available && GetLastError() == ERROR_INVALID_PARAMETER) {
		// Another component has already initialised DbgHelp for this process.
		symbols_available = true;
	}
#endif

	for (USHORT i = 0; i < count; ++i) {
#ifdef HAVE_DBGHELP
		if (symbols_available) {
			char symbol_buffer[sizeof(SYMBOL_INFO) + MAX_SYM_NAME];
			PSYMBOL_INFO symbol = reinterpret_cast<PSYMBOL_INFO>(symbol_buffer);
			symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
			symbol->MaxNameLen = MAX_SYM_NAME;
			DWORD64 displacement = 0;
			if (SymFromAddr(process, reinterpret_cast<DWORD64>(frames[i]),
				&displacement, symbol)) {
				fprintf(output, "  %s + 0x%llx [%p]\n", symbol->Name,
					static_cast<unsigned long long>(displacement), frames[i]);
				continue;
			}
		}
#endif
		fprintf(output, "  %p\n", frames[i]);
	}
	fflush(output);
}

#else

void write_stack_trace(FILE* output) {
	fprintf(output, "\nStack trace unavailable in this build.\n");
	fflush(output);
}

#endif
