///////////////////////////////////////////////////////////////////////////////////////
/// \file framework.h
/// \brief The 'mission control' of the model
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_FRAMEWORK_H
#define LPJ_GUESS_FRAMEWORK_H

/// The 'mission control' of the model
/** 
 *  Responsible for maintaining the primary model data structures and containing
 *  all explicit loops through space (grid cells/stands) and time (days and 
 *  years).
 *
 *  \param argc Number of arguments to the program, sent in from main
 *  \param argv The arguments to the program, sent in from main
 *
 */
int framework(int argc,char* argv[]);

#endif // LPJ_GUESS_FRAMEWORK_H
