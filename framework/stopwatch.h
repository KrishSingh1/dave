#pragma once

#include <chrono>

class Stopwatch;

class StopwatchScope {
public:
    StopwatchScope(Stopwatch &sw);
    ~StopwatchScope();

    StopwatchScope(const StopwatchScope&) = delete;
    StopwatchScope& operator=(const StopwatchScope&) = delete;
private:
    Stopwatch &stopwatch;
};

class Stopwatch {
public:
    /*
    Create the stopwatch and start it immediately.
    */
    Stopwatch();

    /*
    Start the stopwatch. If it is already running, this will have no effect.
    This returns an object which will stop the stopwatch when it goes out of
    scope.
    */
    [[nodiscard]] StopwatchScope start();

    /*
    Start the stopwatch. If it is already running, this will have no effect.
    This is unsafe only in the sense that the stopwatch must later be manually
    stopped.
    */
    void unsafe_start();

    /*
    Pause the stopwatch. If it is already paused, this will have no effect.
    */
    void stop();

    /*
    Reset accumulated time, and start the stopwatch.
    */
    void reset();

    /*
    Get elapsed wall time in seconds.
    */
    double get_elapsed_walltime() const;

    /*
    Get elapsed CPU time in seconds.
    */
    double get_elapsed_cputime() const;
private:
    std::chrono::steady_clock::time_point start_time_wall;
    double start_time_cpu;

    bool running;

    double accumulated_walltime;
    double accumulated_cputime;
};
