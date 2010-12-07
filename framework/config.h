///////////////////////////////////////////////////////////////////////////////////////
// CONFIG HEADER FILE
//
// WHAT SHOULD THIS FILE CONTAIN?
// This header file contains preprocessor definitions accessible throughout the
// model code. Each implementation file (.cpp file) in LPJ-GUESS should include
// this header before anything else. Definitions should be placed here if they
// are of interest when configuring how to build LPJ-GUESS, for instance
// definitions for enabling/disabling a module, or changing a behaviour which
// for some reason isn't configurable from the instruction file.
//
// This file may also contain non-model related code for working around platform
// specific issues, such as non-standard conforming compilers.

// Version dated:         2010-11-22


#ifndef LPJ_GUESS_CONFIG_H
#define LPJ_GUESS_CONFIG_H

// Make sure we use the CRU I/O module unless DEMO I/O has been chosen
#ifndef USE_DEMO_IO
#define USE_CRU
#endif

#endif // LPJ_GUESS_CONFIG_H
