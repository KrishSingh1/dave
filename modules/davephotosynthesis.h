#ifndef LPJ_GUESS_DAVEPHOTOSYNTHESIS_H
#define LPJ_GUESS_DAVEPHOTOSYNTHESIS_H

#include "guess.h"

// Leaf transmission coefficient (Sands, 1995).
const double M_TRANS = 0.0;

// Light extinction coefficient used for grasses in DAVE photosynthesis.
const double K_BEER_GRASS = 0.5;

// Light extinction coefficient used for trees in DAVE photosynthesis.
const double K_BEER_TREE = 0.5;

/*
Top-level DAVE canopy exhange function, called once per day.
*/
void canopy_exchange_dave(Patch &patch);

#endif // LPJ_GUESS_DAVEPHOTOSYNTHESIS_H
