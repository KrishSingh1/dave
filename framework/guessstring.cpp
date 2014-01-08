////////////////////////////////////////////////////////////////////////////////
/// \file guessstring.h
/// \brief Utility functions for working with strings (std::string and char*)
///
/// \author Joe Siltberg
/// $Date$
///
////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guessstring.h"
#include <cctype>

std::string trim(const std::string& str) {
	size_t start_pos = 0;

	while (start_pos < str.size() && isspace(str[start_pos])) {
		++start_pos;
	}

	size_t end_pos = str.size();

	while (end_pos > 0 && isspace(str[end_pos-1])) {
		--end_pos;
	}

	if (start_pos < end_pos) {
		std::string result(str.begin()+start_pos, str.begin()+end_pos);
		return result;
	}
	else {
		return "";
	}
}

std::string to_upper(const std::string& str) {
	std::string result = str;

	for (size_t i = 0; i < result.size(); ++i) {
		result[i] = toupper(result[i]);
	}
	
	return result;
}

std::string to_lower(const std::string& str) {
	std::string result = str;

	for (size_t i = 0; i < result.size(); ++i) {
		result[i] = tolower(result[i]);
	}
	
	return result;
}
