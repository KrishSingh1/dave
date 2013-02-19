///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete_test.cpp
/// \brief Unit tests for nitrogen uptake competition
///
/// \author Joe Lindström
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "ncompete.h"

TEST_CASE("ncompete/single", "Testing a single individual") {
	std::vector<NCompetingIndividual> indivs(1);

	indivs[0].ndemand = 1;
	indivs[0].strength = 1;
	indivs[0].isgrass = false;

	// more N available than needed
	ncompete(indivs, 2.0, 2.0);

	REQUIRE(indivs[0].fnuptake == Approx(1));


	// less N available than needed
	ncompete(indivs, 0.5, 0.5);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));

	// change to grass

	indivs[0].isgrass = true;

	// more N available than needed
	ncompete(indivs, 2.0, 2.0);

	REQUIRE(indivs[0].fnuptake == Approx(1));


	// less N available than needed
	ncompete(indivs, 0.5, 0.5);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));	
}


TEST_CASE("ncompete/double", "Testing two individuals") {
	std::vector<NCompetingIndividual> indivs(2);

	// two equal individuals

	indivs[0].ndemand = 1;
	indivs[0].strength = 1;
	indivs[0].isgrass = false;
	
	indivs[1] = indivs[0];

	// more N available than needed
	ncompete(indivs, 3.0, 3.0/2.0);
	
	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(1));

	// less N available than needed
	ncompete(indivs, 1.0, 1.0/2.0);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));
	REQUIRE(indivs[1].fnuptake == Approx(0.5));

	// change to grass

	indivs[0].isgrass = true;
	indivs[1].isgrass = true;

	// more N available than needed
	ncompete(indivs, 3.0, 3.0/2.0);
	
	REQUIRE(indivs[0].fnuptake == Approx(1));
	REQUIRE(indivs[1].fnuptake == Approx(1));

	// less N available than needed
	ncompete(indivs, 1.0, 1.0/2.0);

	REQUIRE(indivs[0].fnuptake == Approx(0.5));
	REQUIRE(indivs[1].fnuptake == Approx(0.5));

	// make the second indiv twice as strong
	indivs[1].strength *= 2;

	// less N available than needed
	ncompete(indivs, 1.0, 1.0/2.0);

	REQUIRE(indivs[0].fnuptake == Approx(1.0/3.0));
	REQUIRE(indivs[1].fnuptake == Approx(2.0/3.0));

	// reduce the demand for the stronger individual
	indivs[1].ndemand = 0.1;

	// less N available than needed
	ncompete(indivs, 1.0, 1.0/1.1);

	REQUIRE(indivs[0].fnuptake == Approx(0.9));
	REQUIRE(indivs[1].fnuptake == Approx(1));
}
