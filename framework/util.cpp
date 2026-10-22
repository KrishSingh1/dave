#include "util.h"

#ifndef _WIN32
#include <time.h>
#include <sys/time.h>
#endif

double get_wall_time() {
#ifdef _WIN32
	return clock() / CLOCKS_PER_SEC;
#else
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
#endif
}

namespace {

xtring get_nested_error_msg(const std::exception& e, int level) {
       xtring msg;
       std::string spaces(level, ' ');
       msg.printf("%s%s\n", spaces.c_str(), e.what());
       try {
               std::rethrow_if_nested(e);
       } catch (const std::exception& nested) {
               msg += get_nested_error_msg(nested, level + 1);
       } catch (...) {
               // Preserve the original error even if a nested exception cannot be read.
       }
       return msg;
}

} // namespace

xtring get_nested_error_msg(const std::exception& e) {
       return get_nested_error_msg(e, 0);
}
