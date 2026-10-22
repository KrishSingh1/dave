#include "util.h"

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

xtring join_items(const std::set<xtring> &items) {
	xtring result;
	for (std::set<xtring>::const_iterator it = items.begin(); it != items.end(); it++) {
		if (it != items.begin()) {
			result += ", ";
		}
		result += *it;
	}
	return result;
}

xtring join_items(const std::set<xtring> &items) {
	xtring result;
	for (std::set<xtring>::const_iterator it = items.begin(); it != items.end(); it++) {
		if (it != items.begin()) {
			result += ", ";
		}
		result += *it;
	}
	return result;
}
