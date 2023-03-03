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

} // namespace GuessOutput
