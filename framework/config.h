///////////////////////////////////////////////////////////////////////////////////////
/// \file config.h
/// \brief Configuration header file, accessible throughout the model code.
///
/// This header file contains preprocessor definitions accessible throughout the
/// model code. Each implementation file (.cpp file) in LPJ-GUESS should include
/// this header before anything else. Definitions should be placed here if they
/// are of interest when configuring how to build LPJ-GUESS, for instance
/// definitions for enabling/disabling a module, or changing a behaviour which
/// for some reason isn't configurable from the instruction file.
///
/// This file may also contain non-model related code for working around platform
/// specific issues, such as non-standard conforming compilers.
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CONFIG_H
#define LPJ_GUESS_CONFIG_H

// Defines for landcover version:
#define DYNAMIC_LANDCOVER_INPUT		// Reads landcover data from text files, using the TimeDataD class.
#define LUTOMEMORY					// Write land use fraction data to memory; enables efficient usage of randomized gridlists for parallell runs on Simba.
#define NEWSOWINGDATE				// Use sowing date method based on climate seasonality (modified version of Waha et al. 2012), as opposed to old method used in Bondeau et al. 2007.
#define IRRIGATION					// Crop irrigation on
#define NOPASTURESTOCH				// Undefine for fire and disturbance for pasture grass. Number of patches will be the same as for natural stands.
//#define GRASSFORCROP				// Transfer cropland to pasture landcover for simplified crop definition (harvested competing c3/c4 grass).
#define HARVEST_GRSC				// Harvest and/or turnover at the end of the growing season.
//#define NATURALPFTSINFOREST

const bool SUPPRESSLARGEOUTPUT=true;

#define CMASS_SEED 0.01	// 10g/m2;	// Initial carbon allocated to crop organs at sowing.

//#define PRINTFIRSTSTANDFROM1901	// Only on when printout of whole period of first stand wanted.
#define MAXNUMBER_STANDS 100		// Upper limit for multiple stand printout
#define MAXNUMBER_GRIDCELLS 50		// To make sure that muliple stand printout is not active when running large simulations

// Compiler specific checks, for instance for disabling specific warnings

// All versions of Microsoft's compiler
#ifdef _MSC_VER

// 'this' : used in base member initializer list
#pragma warning (disable: 4355)
// long name
#pragma warning (disable: 4786)

#endif

// min and max functions for MS Visual C++ 6.0
#if defined(_MSC_VER) && _MSC_VER == 1200
namespace std {
template <class T> inline T max(const T& a, const T& b) {
    return (a > b) ? a : b;
}

template <class T> inline T min(const T& a, const T& b) {
    return (a < b) ? a : b;
}
}
#else
// modern compilers define min and max in <algorithm>
#include <algorithm>
#endif

using std::min;
using std::max;

// platform independent function for changing working directory
// we'll call our new function change_directory
#ifdef _MSC_VER
// The Microsoft way
#include <direct.h>
#define change_directory _chdir
#else
// The POSIX way
#include <unistd.h>
#define change_directory chdir
#endif

#endif // LPJ_GUESS_CONFIG_H
