///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.cpp
/// \brief Distribution of N among individuals according to supply, demand and
///        the individuals' uptake strength
///
/// \author David Wårlind
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "ncompete.h"
#include "guessmath.h"

void ncompete(std::vector<NCompetingIndividual>& individuals, double nmass_avail) {
	double nsupply = nmass_avail;		// Nitrogen available for uptake
	double grassmin_nsupply = nmass_avail * 0.05; // Minimum grass nitrogen supply. 
										// Grass should at least get 5% of total available nitrogen
										// This will be changed in a future version where more soil
										// layers are introduced and then grass gets this advantage
										// by having higher fraction of roots in top layer
	double grass_ndemand = 0.0;			// Grass total nitrogen demand
	bool grass_satisfied = false;		// If grass gets what it demands from grassmin_nsupply 
	bool grass_competitive = true;		// Keeping track of if grass can compete with trees for 
										// more than what is especially assigned for grass
	double grass_ups = 0.0;				// Uptake strength of grasses
	double total_ups = 0.0;				// Total uptake strength
	double ratio_uptake;				// Nitrogen per uptake strength
	bool full_uptake = true;			// If an individual could get more than its demand, then
										// that indiv gets fnuptake = 1 and everything has to be 
										// redone for all indiv with fnuptake < 1 as more nitrogen 
										// could be taken up per unit strength (starts with true to
										// get into while loop)

	double total_ndemand = 0;
	for (size_t i = 0; i < individuals.size(); ++i) {
		total_ndemand += individuals[i].ndemand;
	}

	// calculate a starting value for the individuals' fnuptake
	double fnuptake = total_ndemand > 0.0 ? nmass_avail / total_ndemand : 0.0;

	// Determine strength and demand
	for (size_t i = 0; i < individuals.size(); ++i) {
		NCompetingIndividual& indiv = individuals[i];

		// All indiv
		// Sum uptake strengths
		if (!negligible(indiv.ndemand)) {
			total_ups += indiv.strength;
		}

		// Grass
		// Sum gras demand demand and uptake strengths
		if (indiv.isgrass && !negligible(indiv.ndemand)) {
			grass_ndemand += indiv.ndemand;
			grass_ups += indiv.strength;
		}
	}

	// Grass
	// Does grass get enough nitrogen from its part of the total
	if (grass_ndemand < grassmin_nsupply) 
		grass_satisfied = true;
	else 
		grass_satisfied = false;

	for (size_t i = 0; i < individuals.size(); ++i) {
		NCompetingIndividual& indiv = individuals[i];

		if (!negligible(total_ups)) {
		 
			indiv.fnuptake = fnuptake;

			// GRASS
			// if grasses gets enough from its part, then take it up
			if (indiv.isgrass && grass_satisfied && !negligible(indiv.ndemand)) {

				// when grass part of total N is enough, then subtract it from total
				nsupply -= indiv.ndemand;
				// set uptake to meet demand
				indiv.fnuptake = 1.0;
				// and subtract uptake strength as it will be added further down
				total_ups -= indiv.strength;
			}
		}
		else {
			indiv.fnuptake = 0.0;
			full_uptake = false;
		}
	}
	
	// Loop through indiv and decide their fnuptake
	while (full_uptake){

		full_uptake = false;	
		
		// restore nitrogen supply and uptake decider if grass_competitive == true
		// (which can happen if there is a full_uptake)
		// so that it can be calculated if they might be able to take up more 
		// than just the grass part
		if (!grass_competitive) {
			nsupply += grassmin_nsupply;
			total_ups += grass_ups;
			grass_competitive = true;
		}

		// decide how much nitrogen that will be taken up by each uptake strength 
		if (total_ups > 0.0 && nsupply > 0.0)
			ratio_uptake = nsupply / total_ups;
		else
			ratio_uptake = 0.0;

		// Grass
		// Grass minimum of avail nitrogen is not enough
		if (!grass_satisfied && grass_competitive) {

			// If grass can't get more than grassmin_nsupply
			if (ratio_uptake * grass_ups < grassmin_nsupply) {

				grass_competitive = false;
				// then grass takes grassmin_nsupply of total nitrogen supply
				nsupply -= grassmin_nsupply;
				// and grass strength is subtracted from total uptake strength
				total_ups -= grass_ups;
				// and a new ratio uptake is calculated for trees 
				ratio_uptake = nsupply / total_ups; 
			}
			else {
				// Grass can compite for more than grassmin_nsupply
				grass_competitive = true;
			}
		}

		for (size_t i = 0; i < individuals.size() && !full_uptake; ++i) {
			NCompetingIndividual& indiv = individuals[i];

			// if lifeform is grass and they can't compete with trees for more than grassmin_nsupply of the total nitrogen supply
			if (indiv.isgrass && !grass_competitive && indiv.fnuptake != 1.0) {
				if (!negligible(indiv.ndemand)) {

					// Calculate new fnuptake
					indiv.fnuptake = grassmin_nsupply * (indiv.strength
						/ grass_ups) / indiv.ndemand;

					// If indiv can get more nitrogen the its demand
					if (indiv.fnuptake > 1.0) {

						indiv.fnuptake = 1.0;

						// subtract nitrogen demand from grass nitrogen supply
						grassmin_nsupply -= indiv.ndemand;

						// and take away this indiv uptake strength from grass total
						grass_ups -= indiv.strength;

						// and redo indiv fnuptake calc for the rest of the indivs as this indiv probably could
						// take up more than its nitrogen demand -> more available for the rest of the indiv
						full_uptake = true;
					}
				}
				else
					indiv.fnuptake = 0.0;
			}

			// if lifeform is tree and grass if it can compete with trees
			else {

				// if fnuptake doesn't meet indiv nitrogen demand, then calculate a new value for fnuptake
				if (indiv.fnuptake != 1.0) {

					// if indiv has the strenght to take up more than nitrogen demand
					if (ratio_uptake * indiv.strength > indiv.ndemand && !negligible(indiv.ndemand)){
						
						indiv.fnuptake = 1.0;
						
						// subtract nitrogen demand from nitrogen supply
						nsupply -= indiv.ndemand;

						// and take away this indiv uptake strength from total
						total_ups -= indiv.strength;

						// and redo indiv fnuptake calc for the rest of the indiv as this indiv probably could
						// take up more than its nitrogen demand -> more available for the rest of the indiv
						full_uptake = true;
					}
					// normal nitrogen limited uptake (0.0 < fnuptake < 1.0)
					else if (indiv.ndemand > 0.0)
						indiv.fnuptake = min(1.0, (ratio_uptake * indiv.strength) / indiv.ndemand);
					else
						indiv.fnuptake = 0.0;
				}
			}
		}
	}
}
