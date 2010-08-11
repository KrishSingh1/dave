///////////////////////////////////////////////////////////////////////////////////////
// MODULE HEADER FILE
//
// Module:                LPJ-GUESS input/output module with input from instruction
//                        script
//                        Version adapted for analyses for joint conceptual paper for
//                        OECD conference with Wolfgang Knorr & Jean-Luc Widlowski
// Header file name:      guessio.h
// Source code file name: guessio.cpp
// Written by:            Ben Smith
// Version dated:         2002-08-11
//
// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GUESSIO_H
#define LPJ_GUESS_GUESSIO_H

#include "guess.h"

void initio(int argc,char* argv[],Pftlist& pftlist);
bool getstand(Stand& stand);
bool getclimate(Stand& stand);
void outannual(Stand& stand,Pftlist& pftlist);
void termio();

#endif // LPJ_GUESS_GUESSIO_H
