#ifndef LPJ_GUESS_LOG_SCOPE_H
#define LPJ_GUESS_LOG_SCOPE_H

#include "guess.h"

// TBI: soil layer, canopy layer
enum class ContextDepth { None, Gridcell, Stand, Patch, Individual };

struct LogContext {
    ContextDepth depth = ContextDepth::None;
    double lon = 0.0;
    double lat = 0.0;
    int stand_id = 0;
    int patch_id = 0;
    int individual_id = 0;
};

extern thread_local LogContext log_context;

class LogScope {
public:
    explicit LogScope(const Gridcell &);
    explicit LogScope(const Stand &);
    explicit LogScope(const Patch &);
    explicit LogScope(const Individual &);

    ~LogScope();

	LogScope(const LogScope &) = delete;
	LogScope& operator=(const LogScope &) = delete;
	LogScope(LogScope &&) = delete;            // optional, but safest
	LogScope& operator=(LogScope &&) = delete; // optional, but safest
private:
    LogScope();
    LogContext previous;
};

#endif // LPJ_GUESS_LOG_SCOPE_H
