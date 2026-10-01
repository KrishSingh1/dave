#include "log_scope.h"

thread_local LogContext log_context;

LogScope::~LogScope() {
    log_context = previous;
}

LogScope::LogScope() : previous(log_context) {
    log_context = LogContext();
}

LogScope::LogScope(const Gridcell &gridcell) : LogScope() {
    log_context.lon = gridcell.get_lon();
    log_context.lat = gridcell.get_lat();
    log_context.depth = ContextDepth::Gridcell;
}

LogScope::LogScope(const Stand &stand) : LogScope(stand.get_gridcell()) {
    log_context.stand_id = stand.id;
    log_context.depth = ContextDepth::Stand;
}

LogScope::LogScope(const Patch &patch) : LogScope(patch.stand) {
    log_context.patch_id = patch.id;
    log_context.depth = ContextDepth::Patch;
}

LogScope::LogScope(const Individual &individual) : LogScope(individual.vegetation.patch) {
    log_context.individual_id = individual.id;
    log_context.depth = ContextDepth::Individual;
}
