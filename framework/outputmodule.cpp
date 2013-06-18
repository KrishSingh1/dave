///////////////////////////////////////////////////////////////////////////////////////
/// \file outputmodule.cpp
/// \brief Implementation of OutputModule and its container class
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "outputmodule.h"
#include "parameters.h"
#include "guess.h"

namespace GuessOutput {

OutputChannel* output_channel;

///////////////////////////////////////////////////////////////////////////////////////
/// OutputModule
///

void OutputModule::create_output_table(Table& table, const char* file, const ColumnDescriptors& columns) {
	 table = output_channel->create_table(TableDescriptor(file, columns));
}


///////////////////////////////////////////////////////////////////////////////////////
/// OutputModuleContainer
///

OutputModuleContainer::OutputModuleContainer() {
	declare_parameter("outputdirectory", &outputdirectory, 300, "Directory for the output files");
}

OutputModuleContainer::~OutputModuleContainer() {
	for (size_t i = 0; i < modules.size(); ++i) {
		delete modules[i];
	}

	delete output_channel;
}

void OutputModuleContainer::add(OutputModule* output_module) {
	modules.push_back(output_module);
}

void OutputModuleContainer::init() {
	// We MUST have an output directory
	if (outputdirectory=="") {
		fail("No output directory given in the .ins file!");
	}

	// Create the output channel
	const int COORDINATES_PRECISION = 1; // decimal places for coords in output
	output_channel = new FileOutputChannel(outputdirectory.c_str(),
	                                       COORDINATES_PRECISION);

	for (size_t i = 0; i < modules.size(); ++i) {
		modules[i]->init();
	}
}

void OutputModuleContainer::outannual(Gridcell& gridcell) {
	for (size_t i = 0; i < modules.size(); ++i) {
		modules[i]->outannual(gridcell);
	}
}

void OutputModuleContainer::outdaily(Gridcell& gridcell) {
	for (size_t i = 0; i < modules.size(); ++i) {
		modules[i]->outdaily(gridcell);
	}
}

} // namespace
