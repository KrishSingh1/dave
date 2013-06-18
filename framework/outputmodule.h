///////////////////////////////////////////////////////////////////////////////////////
/// \file outputmodule.h
/// \brief Base class for output modules and a container class for output modules
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_OUTPUT_MODULE_H
#define LPJ_GUESS_OUTPUT_MODULE_H

#include <vector>
#include <string>
#include "outputchannel.h"

class Gridcell;

namespace GuessOutput {

/// Base class for output modules
/** 
 *  An output module should inherit from this class and implement
 *  the pure virtual functions.
 *
 *  In order for the output module to be used by the framework,
 *  it must be added to a container (\see OutputModuleContainer), 
 *  which is done in the function framework().
 */
class OutputModule {
public:
	/// Called after the instruction file has been read
	/**
	 *  If an output module needs to declare its own instruction
	 *  file parameters, it should declare them in its constructor
	 *  (which will be called before the instrction file is read).
	 *
	 *  Typically, some of the output module's initialization
	 *  can't be done until the instruction file parameters are 
	 *  available, and is therefore done in this init function.
	 */
	virtual void init() = 0;

	/// Called by the framework at the end of the last day of each simulation year
	/** The function should generate simulation results for the year, typically
	 *  by getting values from the gridcell, doing various calculations, and
	 *  sending the results along to the output channel (\see OutputChannel) */
	virtual void outannual(Gridcell& gridcell) = 0;

	/// Called by the framework at the end of each simulation day
	/** Similar to outannual but called every day */
	virtual void outdaily(Gridcell& gridcell) = 0;

protected:

	/// Help function to define_output_tables, creates one output table
	void create_output_table(Table& table, 
	                         const char* file, 
	                         const ColumnDescriptors& columns);
};


/// Manages a list of output modules
/**
 *  Apart from simply storing the output modules, 
 *  this class is also a place for things used by all output
 *  modules (for instance creating the OutputChannel or
 *  declaring the 'outputdirectory' parameter).
 */
class OutputModuleContainer {
public:
	OutputModuleContainer();

	~OutputModuleContainer();

	/// Adds an output module to the container
	/** The container will deallocate the module when
	 *  the container is destructed.
	 */
	void add(OutputModule* output_module);

	/// Calls init on all output modules
	/** Should be called after the instruction file has been read */
	void init();

	/// Calls outannual on all output modules
	void outannual(Gridcell& gridcell);
	
	/// Calls outdaily on all output modules
	void outdaily(Gridcell& gridcell);

private:
	
	/// The output modules
	std::vector<OutputModule*> modules;

	/// Instruction file parameter deciding where to create output files
	std::string outputdirectory;
};

/// The output channel through which all output is sent
/** Currently a global for legacy reasons (in order to not break
 *  existing outannual functions). Should be a member of 
 *  OutputModuleContainer, which is already responsible for
 *  creating and destroying the output_channel.
 */
extern OutputChannel* output_channel;

}

#endif // LPJ_GUESS_OUTPUT_MODULE_H
