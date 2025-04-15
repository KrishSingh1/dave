///////////////////////////////////////////////////////////////////////////////////////
/// \file guessstring.cpp
/// \brief Utility functions for working with strings (std::string and char*)
///
/// \author Joe Siltberg
/// $Date$
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guessstring.h"
#include <cctype>
#include <stdarg.h>
#include <stdio.h>

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

std::string format_string(const char* format, ...) {
	const size_t buffer_size = 4096;
	char buffer[buffer_size];
	va_list args;
	va_start(args, format);
	vsnprintf(buffer, buffer_size, format, args);
	return std::string(buffer);
}

int split_string(char* str) {

	char *p = strtok(str, "\t\n ");
	int count = 0;
	while(p) {
		count++;
		p = strtok(NULL, "\t\n ");
	}

	return count;
}

bool issubstring(const char* string, const char* substring) {

	bool found = false;

	char *p = NULL, string_copy[200] = {0};

	strcpy(string_copy, string);
	p = strtok(string_copy, "\t\n ");
	if(p) {
		if(!strcmp(substring, p)) {
			found = true;
		}
	}

	do {
		p = strtok(NULL, "\t\n ");
		if(p) {
			if(!strcmp(substring, p)) {
				found = true;
			}
		}
	}
	while(p && !found);

	return found;
}

/// @brief Expand environment variables inside string str
/// @param str string with potentially ${envar} in it
/// @return string with ${envvar} expanded
std::string expand_environment_variables( const std::string &str ) {
	
	// easy case, if no ${envvar} present, just return the input
	// this also terminates the recursive calls
    if( str.find( "${" ) == std::string::npos ) return str;

    std::string pre  = str.substr( 0, str.find( "${" ) );

	// we start with everything behind the starting '${'
    std::string post = str.substr( str.find( "${" ) + 2 );

	// strange case, if there is just a '${' in it, not sure what to do, just return the input
    if( post.find( '}' ) == std::string::npos ) return str;

	// extract the envar name
    std::string variable = post.substr( 0, post.find( '}' ) );
    // std::string value = ""; // this would delete the emvar name
    std::string value = variable; // keep the variable, in case it can not get expanded

	// advances post behind '}'
    post = post.substr( post.find( '}' ) + 1 );

    const char *v = getenv( variable.c_str() );
	// if it is defined expand to the value of the env variable
    if( v != NULL ) value = std::string( v );

	// recursive if more ${envar} might be present
    return expand_environment_variables( pre + value + post );
}

// wrapper for xtring around the std::string expand_environment_variables(str)
xtring expand_environment_variables( xtring &str ) {
	std::string tmpstr=std::string((char*) str);
	tmpstr=expand_environment_variables(tmpstr);
	return xtring(tmpstr.c_str());
}


