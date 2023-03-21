#include "daveoutputhelpers.h"
#include "guess.h"
#include "commonoutput.h"

namespace GuessOutput {

void ensure_channel_open() {
	if (!output_channel) {
		fail("Output channel is NULL");
	}
}

double get_lai(Individual &indiv) {
	return indiv.lai;
}

std::vector<Individual *> get_indivs(Vegetation &vegetation, const char *name) {
	vegetation.firstobj();
	std::vector<Individual *> result;
	while (vegetation.isobj) {
		Individual &indiv = vegetation.getobj();

		if (indiv.pft.name == name) {
			result.push_back(&indiv);
		}

		vegetation.nextobj();
	}
	return result;
}

std::vector<std::string> get_pfts() {
	std::vector<std::string> pfts;
	pftlist.firstobj();
	while (pftlist.isobj) {
		pfts.push_back((char *)pftlist.getobj().name);
		pftlist.nextobj();
	}
	return pfts;
}

} // namespace GuessOutput
