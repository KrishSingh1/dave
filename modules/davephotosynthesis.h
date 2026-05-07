#ifndef LPJ_GUESS_DAVEPHOTOSYNTHESIS_H
#define LPJ_GUESS_DAVEPHOTOSYNTHESIS_H

#include "guess.h"

/*
Top-level DAVE canopy exhange function, called once per day.
*/
void canopy_exchange_dave(Patch &patch, Climate &climate);

#endif // LPJ_GUESS_DAVEPHOTOSYNTHESIS_H
