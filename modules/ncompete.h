///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.h
/// \brief Distribution of nutrients among individuals according to supply, demand and
///        the individuals' uptake strength
///
/// \author David Wårlind
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_NCOMPETE_H
#define LPJ_GUESS_NCOMPETE_H

#include <vector>

/// Represents an individual competing for nutrient uptake
/** Contains what the ncompete function below needs to know about
 *  an individual in order to do the distribution of nutrient.
 */
struct NCompetingIndividual {
	/// This individual's nutrient demand
	double demand;

	/// A meassure of this individual's uptake strength
	double strength;

	/// Output from ncompete - fraction of the demand satisfied by the distribution
	double fuptake;
};

/// Distributes nutrients among individuals according to supply, demand and uptake strength
/** Grasses should get at least 5% and no
 *  individual should get more than 100% of its nutrient demand. 
 */
void ncompete(std::vector<NCompetingIndividual>& individuals, 
              double nutrient_mass);

#endif // LPJ_GUESS_NCOMPETE_H
