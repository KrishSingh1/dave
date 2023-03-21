///////////////////////////////////////////////////////////////////////////////////////
/// \file dave.h
/// \brief Output module for DAVE outputs, including per-patch outputs.
///
/// \author Drew Holzworth
/// $Date: 2023-02-17 10:58:54 +1100 (Fri, 17 Feb 2023) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_DAVE_OUTPUT_H
#define LPJ_GUESS_DAVE_OUTPUT_H

#include <memory>

#include "outputmodule.h"
#include "guess.h"

namespace GuessOutput {

const int COLUMN_WIDTH = 14;
const int NUMERIC_PRECISION = 8;
const char *COLNAME_TOTAL = "total";
const char *COLNAME_PATCH = "patch";
const char *COLNAME_TIMESTEP = "timestep";

// Output module for DAVE outputs, including per-patch outputs.
class DaveOutput : public OutputModule {
public:
	/*
	Constructor.
	*/
	DaveOutput();

	/*
	Destructor.
	*/
	~DaveOutput();

	// OutputModule implementation. See base class for documentation.

	void init();
	void outannual(Gridcell&);
	void outdaily(Gridcell&);
	void openlocalfiles(Gridcell &, int) {}
	void closelocalfiles(Gridcell&) {}
private:
	class OutputFile {
		public:
			/*
			Virtual destructor.
			*/
			virtual ~OutputFile() = 0;

			/*
			Register parameters in the instruction file which will be required
			by this output file instance.
			*/
			void declare_parameters();

			/*
			Open the file connection and register this output table with the
			specified output channel.

			@param channel: The output channel.
			*/
			void open(OutputChannel &channel);

			/*
			Close the file connection.

			@param channel: The output channel.
			*/
			void close(OutputChannel &channel);

			/*
			Write all rows of output to the output file for the current
			timestep. No default implementation is provided, because the
			implementation depends on the output - e.g. annual, daily, daily
			per-patch, etc.

			@param OutputRows: The output channel used to write output.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel &channel, Gridcell &gridcell) = 0;

			/*
			Check if this output file is enabled or not. It is enabled iff the
			user has provided a file name in the instruction file.
			*/
			bool enabled();
		protected:
			/*
			Create a new OutputFile instance. This is protected because an
			OutputFile is abstract.

			@param param: Name of the paramter in the instruction file which determines the output file name.
			@param desc: Description of the output file (please include units).
			@param spinup_output: True to write outputs during spinup. False otherwise.
			*/
			OutputFile(std::string param, std::string desc);

			/*
			Write the next output value to the file.

			@param rows: Output rows manager object.
			@param value: The value to be written.
			*/
			void write_output(OutputRows &rows, double value);

			/*
			Write the next output value to the file.

			@param rows: Output rows manager object.
			@param value: The value to be written.
			*/
			void write_output(OutputRows &rows, int value);

			/*
			Finish the current row of output.

			@param rows: Output rows manager object.
			@param lon: Longitude of gridcell.
			@param lat: Latitude of gridcell.
			@param year: Current calendar year.
			@param day: Current day of year (0-364).
			*/
			void finish_row(OutputChannel &rows, double lon, double lat, int year, int day);

			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_columns() = 0;
		private:
			/*
			Name of the instruction file parameter controlling the file name.
			*/
			std::string param_name;

			/*
			Output file name.
			*/
			std::string file_name;

			/*
			Description of the output file (please put units in this!!).
			*/
			std::string description;

			/*
			Columns in the output file.
			*/
			ColumnDescriptors columns;

			/*
			The created output table. This will be NULL until open() is called.
			Note that this is allocated on the stack and shouldn't be deleted.
			*/
			Table *table;

			/*
			Check if an output should be written to the output rows object.

			@param out: The output rows object.
			*/
			bool should_write(OutputRows &out);

			/*
			Assert that the output table is open. Terminate the simulation if not.
			*/
			void assert_table_open();
	};

	/*
	Base class for patch-level outputs.
	*/
	class PatchLevelOutput : public OutputFile {
		public:
			/*
			Destructor.
			*/
			virtual ~PatchLevelOutput() = 0;

			/*
			Write all rows of output to the output file for the current
			timestep. No default implementation is provided, because the
			implementation depends on the output - e.g. annual, daily, daily
			per-patch, etc.

			@param OutputRows: The output channel used to write output.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel &channel, Gridcell &gridcell);
		protected:
			/*
			Create a new PatchLevelOutput instance. This is protected because
			PatchLevelOutput is abstract.

			@param param: Name of the paramter in the instruction file which determines the output file name.
			@param desc: Description of the output file (please include units).
			@param cols: Columns in the output file.
			@param spinup_output: True to write outputs during spinup. False otherwise.
			*/
			PatchLevelOutput(std::string param, std::string desc);

			/*
			Write all required row for the specified patch.

			@param patch: The patch for which outputs should be written.
			@param out: Output rows instance to which outputs will be written.
			*/
			virtual void write_patch(Patch &patch, OutputRows &out) = 0;
	};

	class IndividualOutput : public PatchLevelOutput {
		public:
			/*
			Create a new individual output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_output: A function which returns the output value from an individual.
			@param description: A description of the output (including units!!).
			@param spinup_output: True to write outputs during spinup. False otherwise.
			*/
			IndividualOutput(std::string param_name, double (*get_output)(Individual &), std::string description);

			/*
			Destructor.
			*/
			~IndividualOutput();

			/*
			Write a row of outputs containing one cell per individual.

			@param patch: The patch for which outputs should be written.
			@param out: Output rows instance to which outputs will be written.
			*/
			virtual void write_patch(Patch &patch, OutputRows &out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			virtual ColumnDescriptors get_columns();

		private:

			/*
			Pointer to a function which reads an output from a pointer to an
			object of the template type.

			@param indiv: An individual.
			@return: The value which appears in each column of the output table.
			*/
			double (*read_output)(Individual &indiv);
	};

	/*
	This class writes an array to an output file. It writes one row of output
	per array element. The array can (potentially) have a different number of
	elements for each patch. One column of output is written per individual.
	*/
	class IndivLevelArrayOutput : public OutputFile {
		public:
			/*
			Create a new individual output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param element_name: A name which describes the elements of the array (e.g. "layer"). This will be used as a column name.
			@param get_nelem: Pointer to a function which returns the number of array elements for the specified patch. The assumption is that every individual in the patch has an output array of the same length.
			@param get_output: A function which returns the output value from an individual.
			@param description: A description of the output (including units!!).
			@param spinup_output: True to write outputs during spinup. False otherwise.
			*/
			IndivLevelArrayOutput(std::string param_name, std::string element_name, int (*get_nelem)(Patch &), double (*get_output)(Individual &, int), std::string description);

			/*
			Destructor.
			*/
			~IndivLevelArrayOutput() { }

			/*
			Write all rows of output to the output file for the current
			timestep. One row will be written for each array element. One column
			will be written for each individual.

			@param OutputRows: The output channel used to write output.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel &channel, Gridcell &gridcell) override;
		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			ColumnDescriptors get_columns() override;

		private:
			/*
			Name of each element in the array (e.g. "layer").
			*/
			std::string elem_name;

			/*
			Get the number of elements in the output array.
			*/
			int (*get_n)(Patch &);

			/*
			Read the i-th element of the output array.
			*/
			double (*read_elem)(Individual &, int);
	};

	/*
	This class writes an array output file. One row is written per timestep
	(which can be yearly, daily, etc). One column is written per array element.
	*/
	class PatchLevelArrayOutput : public PatchLevelOutput {
		public:
			/*
			Create a new patch-level array output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param col_name: Name of the columns in the output file (e.g. sw => sw0, sw1, sw2, ...).
			@param n: Number of elements in the array.
			@param get_output: A function which returns the n-th element of the array from the patch.
			@param description: A description of the output (including units!!).
			@param spinup_output: True to write outputs during spinup. False otherwise.
			@param sum: If true, the sum of the data columns will be stored in the total column. If false, the mean will be used instead.
			*/
			PatchLevelArrayOutput(std::string param_name, std::string col_name, int n, double (*get_output)(Patch &, int), std::string description, bool sum = true);

			/*
			Destructor.
			*/
			~PatchLevelArrayOutput();

			/*
			Write a row of outputs containing one cell per array element.

			@param patch: The patch for which outputs should be written.
			@param out: Output rows instance to which outputs will be written.
			*/
			void write_patch(Patch &patch, OutputRows &out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			ColumnDescriptors get_columns();

		private:
			/*
			Number of elements in the output array.
			*/
			int nelem;

			/*
			Name of the columns in the output file.
			*/
			std::string column_name;

			/*
			Pointer to a function which reads an element of the output array
			from the specified patch.

			@param patch: The patch.
			@param i: Array index in range [0, n)
			@return: The i-th element of the output array.
			*/
			double (*read_output)(Patch &patch, int i);

			/*
			If true, store the sum of the data columns in the total column. If
			false, the mean will be stored instead.
			*/
			bool total_is_sum;
	};

	/*
	Write the specified outputs.
	*/
	void write_outputs(Gridcell&, std::vector<std::unique_ptr<OutputFile>>&);

	/*
	Daily outputs.
	*/
	std::vector<std::unique_ptr<OutputFile>> daily_outputs;

	/*
	Annual outputs.
	*/
	std::vector<std::unique_ptr<OutputFile>> annual_outputs;

	/*
	Add a daily individual-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param elem_name: Name of the index column (e.g. "layer", "timestep", etc).
	@param get_nelem: Get the number of array elements for a given patch.
	@param read_output: A function which returns the i-th element of the output property of an individual.
	@param description: Description of the output (with units!!).
	*/
	void add_daily_indiv_array_output(std::string name, std::string elem_name, int (*get_nelem)(Patch &), double (*read_output)(Individual&, int), std::string description);

	/*
	Add an annual individual-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param elem_name: Name of the index column (e.g. "layer", "timestep", etc).
	@param get_nelem: Get the number of array elements for a given patch.
	@param read_output: A function which returns the i-th element of the output property of an individual.
	@param description: Description of the output (with units!!).
	*/
	void add_annual_indiv_array_output(std::string name, std::string elem_name, int (*get_nelem)(Patch &), double (*read_output)(Individual&, int), std::string description);

	/*
	Add a daily individual-level array output.

	This is really just a wrapper around add_daily_indiv_array_output().

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns the i-th subdaily element of the output property of an individual.
	@param description: Description of the output (with units!!).
	*/
	void add_subdaily_indiv_output(std::string name, double (*read_output)(Individual&, int), std::string description);

	/*
	Add a daily individual-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	*/
	void add_daily_indiv_output(std::string name, double (*read_output)(Individual&), std::string description);

	/*
	Add an annual individual-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	*/
	void add_annual_indiv_output(std::string name, double (*read_output)(Individual&), std::string description);

	/*
	Add a daily patch-level array output.

	This is really just a wrapper around add_daily_patch_level_array_output().

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	*/
	void add_subdaily_patch_level_output(std::string name, double (*read_output)(Patch&, int), std::string description);

	/*
	Add a daily patch-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param n: Number of array elements.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param total_is_sum: If true, the sum of all data columns will be written to the total column. If false, the mean will be written instead.
	*/
	void add_daily_patch_level_array_output(std::string name, std::string colname, int nelem, double (*read_output)(Patch&, int), std::string description, bool total_is_sum = true);

	/*
	Add an annual patch-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param n: Number of array elements.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	*/
	void add_annual_patch_level_array_output(std::string name, std::string colname, int nelem, double (*read_output)(Patch&, int), std::string description);

	/*
	Register instruction file parameters for all output files defined by this
	output module.
	*/
	void declare_parameters();

	/*
	Open all enabled output files ready for writing.
	*/
	void open_outfiles();

	/*
	Open all enabled output files in the specified vector.

	@param outputs: the output files to be opened.
	*/
	void open_outfiles(std::vector<std::unique_ptr<OutputFile>> &outputs);

	/*
	Register all instruction file parameters for the output files defined in the
	given vector.

	@param outputs: The outputs to be registered.
	*/
	void declare_parameters(std::vector<std::unique_ptr<OutputFile>> &outputs);

	/*
	Create an individual output pointer using the specified name, function
	pointer and description.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	*/
	std::unique_ptr<OutputFile> create_individual_output(std::string name, double (*read_output)(Individual&), std::string description);

	/*
	Create a pointer to a patch-level array output. This is owned by the caller.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param n: Number of array elements.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param total_is_sum: If true, the sum of all data columns will be written to the total column. If false, the mean will be written instead.
	*/
	std::unique_ptr<OutputFile> create_array_output(std::string name, std::string colname, int n, double (*read_output)(Patch&, int), std::string description, bool total_is_sum = true);

	/*
	Create a pointer to an individual-level array output. This is owned by the caller.

	@param name: Parameter in the instruction file which determines output file name.
	@param elem_name: Name of the index column (e.g. "layer", "timestep", etc).
	@param get_nelem: Get the number of array elements for a given patch.
	@param read_output: A function which returns the i-th element of the output property of an individual.
	@param description: Description of the output (with units!!).
	*/
	std::unique_ptr<OutputFile> create_indiv_array_output(std::string name, std::string elem_name, int (*get_nelem)(Patch &), double (*read_output)(Individual&, int), std::string description);

	/*
	Total time spent in dave output code.
	*/
	double total_time;
};

/*
True iff outputs should be written during spinup.
*/
extern bool spinup_output;

}

#endif // LPJ_GUESS_DAVE_OUTPUT_H
