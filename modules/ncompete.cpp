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
#include <assert.h>
#include "driver.h"

void ncompete(std::vector<NCompetingIndividual>& individuals, double nmass_avail) {

	double nsupply = nmass_avail;		// Nitrogen available for uptake

	//double grassmin_nsupply = nmass_avail * 0.05; // Minimum grass nitrogen supply. 
										// Grass should at least get 5% of total available nitrogen
										// This will be changed in a future version where more soil
										// layers are introduced and then grass gets this advantage
										// by having higher fraction of roots in top layer

	//double grass_ndemand = 0.0;			// Grass total nitrogen demand
	//bool grass_competitive = true;		// Keeping track of if grass can compete with trees for 
										// more than what is especially assigned for grass

	//double grass_ups = 0.0;				// Uptake strength of grasses
	double total_ups = 0.0;				// Total uptake strength

	double ratio_uptake;				// Nitrogen per uptake strength
	bool full_uptake = true;			// If an individual could get more than its demand, then
										// that indiv gets fnuptake = 1 and everything has to be 
										// redone for all indiv with fnuptake < 1 as more nitrogen 
										// could be taken up per unit strength (starts with true to
										// get into while loop)

	/*double temp_ndemand = 0.0;
	double temp_ups = 0.0;*/

	// calculate total nitrogen demand
	double total_ndemand = 0;
	for (size_t i = 0; i < individuals.size(); ++i) {
		total_ndemand += individuals[i].ndemand;
	}

	// calculate a starting value for the individuals' fnuptake
	double fnuptake = total_ndemand > 0.0 ? min(1.0, nmass_avail / total_ndemand) : 0.0;

	// Determine strength and set initial individual fnuptake
	for (size_t i = 0; i < individuals.size(); ++i) {
		NCompetingIndividual& indiv = individuals[i];

		indiv.fnuptake = fnuptake;

		// All indiv
		// Sum uptake strengths
		if (!negligible(indiv.ndemand)) {
			total_ups += indiv.strength;
		}

		//// Grass
		//// Sum gras demand demand and uptake strengths
		//if (indiv.isgrass && !negligible(indiv.ndemand)) {
		//	grass_ndemand += indiv.ndemand;
		//	grass_ups += indiv.strength;
		//}
	}

	//// GRASSES
	//// Does grass get enough nitrogen from its part of the total
	//if (grass_ndemand < grassmin_nsupply) {
	//	for (size_t i = 0; i < individuals.size(); ++i) {
	//		NCompetingIndividual& indiv = individuals[i];

	//		if (indiv.isgrass) {

	//			indiv.fnuptake = 1.0;

	//			// subtract nitrogen demand from nitrogen supply
	//			nsupply -= indiv.ndemand;
	//			// and grass strength is subtracted from total uptake strength
	//			total_ups -= indiv.strength;
	//		}
	//	}

	//	// no need of considering grasses any more
	//	grassmin_nsupply = 0.0;
	//	grass_ups = 0.0;
	//}
	
	// Loop through indiv and decide their fnuptake
	while (full_uptake){

		full_uptake = false;	

		//// reset value as grass can now be competitive
		//nsupply += temp_ndemand;
		//grassmin_nsupply+= temp_ndemand;
		//temp_ndemand = 0.0;
		//total_ups += temp_ups;
		//grass_ups += temp_ups;
		//temp_ups = 0.0;

		// decide how much nitrogen that will be taken up by each uptake strength 
		if (total_ups > 0.0) 
			ratio_uptake = nsupply / total_ups;
		else
			ratio_uptake = 0.0;

		//// Is grass competitive and able to take up more than grassmin_nsupply
		//if (ratio_uptake * grass_ups >= grassmin_nsupply && !negligible(grassmin_nsupply)) 
		//	grass_competitive = true;
		//else
		//	grass_competitive = false;

		//bool grass_full = true;

		//// Check grasses
		//while(grass_full && !grass_competitive) {

		//	grass_full = false;

		//	for (size_t i = 0; i < individuals.size() && !grass_full; ++i) {
		//		NCompetingIndividual& indiv = individuals[i];

		//		if (indiv.isgrass && indiv.fnuptake != 1.0) {

		//			if (!negligible(indiv.ndemand)) {

		//				// Calculate new fnuptake
		//				if (!negligible(grassmin_nsupply)) {
		//					indiv.fnuptake = grassmin_nsupply * 
		//						(indiv.strength	/ grass_ups) / indiv.ndemand;
		//				}
		//				else {	// no grassmin_nsupply mean other grass indiv have taken all, now have to compete with rest of indivs
		//					indiv.fnuptake = min(1.0, (ratio_uptake * indiv.strength) / indiv.ndemand);
		//				}

		//				if (indiv.fnuptake > 1.0) {

		//					grass_full = true;
		//					indiv.fnuptake = 1.0;

		//					// subtract nitrogen demand from nitrogen supplys
		//					nsupply          -= indiv.ndemand;
		//					grassmin_nsupply -= indiv.ndemand;

		//					// subtract strength from uptake strengths
		//					total_ups        -= indiv.strength;
		//					grass_ups        -= indiv.strength;
		//				}
		//				else {

		//					// subtract nitrogen demand from nitrogen supplys
		//					// and save if some indiv have fnuptake > 1
		//					double demand_temp = indiv.fnuptake * indiv.ndemand;
		//					nsupply           -= demand_temp;
		//					grassmin_nsupply  -= min(grassmin_nsupply, demand_temp);
		//					temp_ndemand      += demand_temp;

		//					// subtract strength from uptake strengths
		//					// and save if some indiv have fnuptake > 1
		//					total_ups         -= indiv.strength;
		//					grass_ups         -= indiv.strength;
		//					temp_ups          += indiv.strength;
		//				}
		//			}
		//			else {
		//				indiv.fnuptake = 1.0;
		//			}
		//		}
		//	}
		//}

		//// update ratio_uptake if grasses has taken up nitrogen from their minimum portion
		//if (!grass_competitive) {
		//	if (total_ups > 0.0)
		//		ratio_uptake = nsupply / total_ups;
		//	else
		//		ratio_uptake = 0.0;
		//}
		
		if (!full_uptake) {

			// Go through individuals
			for (size_t i = 0; i < individuals.size() && !full_uptake; ++i) {
				NCompetingIndividual& indiv = individuals[i];

				//// Only calculate for trees and grass that can compete
				//if (!indiv.isgrass || (indiv.isgrass && grass_competitive)) {
					
					// if fnuptake doesn't meet indiv nitrogen demand, then calculate a new value for fnuptake
					if (indiv.fnuptake != 1.0) {

						// if indiv has the strenght to take up more than nitrogen demand
						if (ratio_uptake * indiv.strength > indiv.ndemand && !negligible(indiv.ndemand)){

							indiv.fnuptake = 1.0;

							// nitrogen demand is subtract from nitrogen supply
							nsupply -= indiv.ndemand;

							// strength is subtracted from uptake strengths
							total_ups -= indiv.strength;

							//if (indiv.isgrass) {

							//	// nitrogen demand is subtract from grass nitrogen supply
							//	grassmin_nsupply -= min(grassmin_nsupply, indiv.ndemand);
							//	// strength is subtracted from grass uptake strengths
							//	grass_ups -= indiv.strength;
							//}

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
				//}
			}
		}
	}

	// check so nitrogen uptake doesn't exceed nitrogen supply
	double nuptake = 0.0;
	for (size_t i = 0; i < individuals.size(); ++i) {
		NCompetingIndividual& indiv = individuals[i];

		nuptake += indiv.fnuptake * indiv.ndemand;
	}

	// maximum nitrogen uptake mismatch
	double EPS = 1.0e-13;

	if (nmass_avail < total_ndemand) // tests might have higher nmass_avail then total_ndemand
		assert(abs(nuptake - nmass_avail) < EPS);
	else
		assert(abs(nuptake - total_ndemand) < EPS);
}
