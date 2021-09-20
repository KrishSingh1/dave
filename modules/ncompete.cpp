///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.cpp
/// \brief Distribution of nutrient among individuals according to supply, demand and
///        the individuals' uptake strength
///
/// \author David Wårlind
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "ncompete.h"
#include "guessmath.h"
#include <assert.h>
#include "driver.h"

void ncompete(std::vector<NCompetingIndividual>& individuals, double nutrient_mass) {
	double nsupply = nutrient_mass;		// nutrient available for uptake
	double total_ups = 0.0;				// Total uptake strength

	// calculate total nutrient uptake strength, and set all uptake to
	// zero initially
	for (size_t i = 0; i < individuals.size(); ++i) {
		total_ups += individuals[i].strength;
		individuals[i].fuptake = 0;
	}

	bool full_uptake = true;			// If an individual could get more than its demand, then
										// that indiv gets fnuptake = 1 and everything has to be 
										// redone for all indiv with fnuptake < 1 as more nutrient 
										// could be taken up per unit strength (starts with true to
										// get into while loop)

	while (full_uptake) {
		full_uptake = false;

		double ratio_uptake;				// nutrient per uptake strength

		// decide how much nutrient that will be taken up by each uptake strength 
		if (total_ups > 0.0) {
			ratio_uptake = nsupply / total_ups;
		}
		else {
			ratio_uptake = 0.0;
		}

		// Go through individuals
		for (size_t i = 0; i < individuals.size() && !full_uptake; ++i) {
			NCompetingIndividual& indiv = individuals[i];

			// if fnuptake doesn't meet indiv nutrient demand, then calculate a new value for fnuptake
			if (indiv.fuptake != 1.0) {

				// how much is the individual entitled to take up due to its strength?
				double entitlement = ratio_uptake * indiv.strength;

				// if indiv has the strength to take up more than its nutrient demand
				if (entitlement > indiv.demand && !negligible(indiv.demand)) {

					indiv.fuptake = 1.0;

					// nutrient demand is subtract from nutrient supply
					nsupply -= indiv.demand;

					// strength is subtracted from uptake strengths
					total_ups -= indiv.strength;

					// and redo indiv fnuptake calc for the rest of the indiv as this indiv probably could
					// take up more than its nutrient demand -> more available for the rest of the indiv
					full_uptake = true;
				}
				// normal nutrient limited uptake (0.0 < fnuptake < 1.0)
				else if (indiv.demand > 0.0) {
					indiv.fuptake = min(1.0, entitlement / indiv.demand);
				}
				else {
					indiv.fuptake = 0.0;
				}
			}
		}
	}
}
