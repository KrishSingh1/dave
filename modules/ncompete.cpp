///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.cpp
/// \brief Distribution of N among individuals according to supply, demand and the individuals' uptake strength
///
/// \author David Wårlind
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
		individuals[i].fnuptake = 0;
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
			if (indiv.fnuptake != 1.0) {

				// how much is the individual entitled to take up due to its strength?
				double entitlement = ratio_uptake * indiv.strength;

				// if indiv has the strength to take up more than its nutrient demand
				if (entitlement > indiv.ndemand && !negligible(indiv.ndemand)) {

					indiv.fnuptake = 1.0;

					// nutrient demand is subtract from nutrient supply
					nsupply -= indiv.ndemand;

					// strength is subtracted from uptake strengths
					total_ups -= indiv.strength;

					// and redo indiv fnuptake calc for the rest of the indiv as this indiv probably could
					// take up more than its nutrient demand -> more available for the rest of the indiv
					full_uptake = true;
				}
				// normal nutrient limited uptake (0.0 < fnuptake < 1.0)
				else if (indiv.ndemand > 0.0) {
					indiv.fnuptake = min(1.0, entitlement / indiv.ndemand);
				}
				else {
					indiv.fnuptake = 0.0;
				}
			}
		}
	}
}
