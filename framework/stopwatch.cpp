#include "stopwatch.h"

#ifdef _WIN32

#include <processthreadsapi.h>
#include <windows.h>

std::string GetLastErrorAsString() {
	// Get the error ID.
	DWORD errorMessageID = ::GetLastError();
	if (errorMessageID == 0) {
		return std::string();
	}

	LPSTR messageBuffer = nullptr;

	size_t size = FormatMessageA(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		errorMessageID,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPSTR)&messageBuffer,
		0,
		NULL);

    if (size == 0) {
        // If FormatMessage fails, return a generic error message.
        return "Unknown error";
    }

	std::string message(messageBuffer, size);

	LocalFree(messageBuffer);
	return message;
}

#else

#include <sys/resource.h>

#endif

/*
Return current CPU time (in seconds) used by the process. This is the sum of
user and system time.
*/
double get_cpu_time() {
#ifdef _WIN32
	FILETIME creation_time, exit_time, kernel_time, user_time;
	if (GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time,
	                    &kernel_time, &user_time) == 0) {
		std::string message = GetLastErrorAsString();
		throw std::runtime_error("GetProcessTimes failed: " + message);
		return 0;
	}

	ULARGE_INTEGER user, kernel;
    user.LowPart = user_time.dwLowDateTime;
    user.HighPart = user_time.dwHighDateTime;

    kernel.LowPart = kernel_time.dwLowDateTime;
    kernel.HighPart = kernel_time.dwHighDateTime;

	return (static_cast<double>(user.QuadPart) +
            static_cast<double>(kernel.QuadPart)) * 1e-7;
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        // Handle error
        return 0;
    }
	// Get ru_stime

    return (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec * .000001 +
		   (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec * .000001;
#endif
}

StopwatchScope::StopwatchScope(Stopwatch &sw) : stopwatch(sw) {
}

StopwatchScope::~StopwatchScope() {
    stopwatch.stop();
}

Stopwatch::Stopwatch() {
    reset();
}

StopwatchScope Stopwatch::start() {
    unsafe_start();
    return StopwatchScope{*this};
}

void Stopwatch::unsafe_start() {
    if (!running) {
        running = true;
        start_time_wall = std::chrono::steady_clock::now();
        start_time_cpu = get_cpu_time();
    }
}

void Stopwatch::stop() {
    if (running) {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = now - start_time_wall;
        accumulated_walltime += elapsed.count();

        double now_cpu = get_cpu_time();
        accumulated_cputime += (now_cpu - start_time_cpu);

        running = false;
    }
}

void Stopwatch::reset() {
    accumulated_cputime = 0;
    accumulated_walltime = 0;
    running = true;
    start_time_wall = std::chrono::steady_clock::now();
    start_time_cpu = get_cpu_time();
}

double Stopwatch::get_elapsed_walltime() const {
    if (running) {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = now - start_time_wall;
        return accumulated_walltime + elapsed.count();
    } else {
        return accumulated_walltime;
    }
}

double Stopwatch::get_elapsed_cputime() const {
    if (running) {
        double now_cpu = get_cpu_time();
        return accumulated_cputime + (now_cpu - start_time_cpu);
    } else {
        return accumulated_cputime;
    }
}
