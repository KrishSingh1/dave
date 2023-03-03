#include "daveoutputhelpers.h"
#include "guess.h"
#include "commonoutput.h"

namespace GuessOutput {

/*
Ensure the output channel is open for writing. Abort the simualtion if not.
*/
void ensure_channel_open() {
	if (!output_channel) {
		fail("Output channel is NULL");
	}
}

/*
This is an example of a helper function which may be used instead of a lambda
expression when defining an OutputFile.
*/
double get_lai(Individual &indiv) {
	return indiv.lai;
}

Individual *get_indiv(Vegetation &vegetation, const char *name) {
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual &indiv = vegetation.getobj();

		if (indiv.pft.name == name) {
			return &indiv;
		}

		vegetation.nextobj();
	}
	return NULL;
}

Individual *get_first_indiv(Patch &patch) {
	Vegetation &vegetation = patch.vegetation;
	vegetation.firstobj();
	if (vegetation.isobj) {
		return &(vegetation.getobj());
	}
}

int get_nlayer_fpar(Patch &patch) {
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->fpar_layerwise.size();
	}
	// No individuals are established in this patch.
	return 0;
}


int get_nlayer_lai(Patch &patch) {
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->lai_layerwise.size();
	}
	return 0;
}

} // namespace GuessOutput
