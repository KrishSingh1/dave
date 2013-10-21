////////////////////////////////////////////////////////////////////////////////
/// \file partitionedmapserializer.cpp
/// \brief Implementation file for PartitionedMapSerializer/Deserializer
///
/// Since these classes are templates, most of their implementation is in the
/// header.
///
/// \author Joe Siltberg
/// $Date$
///
////////////////////////////////////////////////////////////////////////////////

#include "partitionedmapserializer.h"
#include <sstream>

std::string create_path(const char* directory, int rank) {
	std::ostringstream os;
	os << directory << "/" << rank << ".state";
	return os.str();
}
