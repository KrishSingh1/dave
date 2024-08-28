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

extern const int COLUMN_WIDTH;
extern const int COLUMN_WIDTH_PFT;
extern const int NUMERIC_PRECISION;
extern const char *COLNAME_TOTAL;
extern const char *COLNAME_MEAN;
extern const char *COLNAME_PATCH;
extern const char *COLNAME_STAND;
extern const char *COLNAME_TIMESTEP;

template<class T>
class ColumnDefinition {
	public:
		/*
		Create a new column definition.

		@param name: Column name.
		@param get_output: Pointer to a function which returns the value to
							be written to the rows of this column.
		*/
		ColumnDefinition(std::string name, double (*get_output)(const T&)) {
			this->name = name;
			read_output = get_output;
		}

		/*
		Get the name of the column.
		*/
		std::string get_name() const {
			return name;
		}

		/*
		Read the current value of the column.
		*/
		double read_value(const T &t) const {
			return read_output(t);
		}

	private:
		std::string name;

		/*
		Pointer to a function which reads an output from a pointer to an
		object of the template type.

		@param patch: The patch.
		@return: The value which appears in each column of the output table.
		*/
		double (*read_output)(const T &t);
};

// Output module for DAVE outputs, including per-patch outputs.
template<class T>
class OutputReader {
	public:
		OutputReader(double (*get_value)(const T&)) {
			this->get_value = get_value;
		}
		double read(const T& t) const {
			return get_value(t);
		}
	private:
		double (*get_value)(const T&);
};

template<class T>
class ArrayOutputReader {
	public:
		ArrayOutputReader(double (*get_value)(const T&, int)) {
			this->get_value = get_value;
			this->get_value_simple = NULL;
		}
		ArrayOutputReader(double (*get_value)(const T&)) {
			this->get_value = NULL;
			this->get_value_simple = get_value;
		}
		virtual double read(const T& t, int i) const {
			return get_value ? get_value(t, i) : get_value_simple(t);
		}
	private:
		double (*get_value)(const T&, int);
		double (*get_value_simple)(const T&);
};

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
	typedef enum {SUBDAILY, DAILY, ANNUAL} timestep_t;

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

			@param out: Pointer to an output channel.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel *channel, Gridcell &gridcell) = 0;

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

			@param channel: Pointer to an output channel.
			@param value: The value to be written.
			*/
			void write_output(OutputChannel *channel, double value);

			/*
			Write the next output value to the file.

			@param channel: Pointer to an output channel.
			@param value: The value to be written.
			*/
			void write_output(OutputChannel *channel, int value);

			/*
			Write the text output value to the output file.

			@param channel: Pointer to an output channel.
			@param value: The value to be written.
			*/
			void write_output(OutputChannel *channel, const char *value);

			/*
			Finish the current row of output.

			@param channel: Pointer to an output channel.
			@param lon: Longitude of gridcell.
			@param lat: Latitude of gridcell.
			@param year: Current calendar year.
			@param day: Current day of year (0-364).
			*/
			void finish_row(OutputChannel *channel, double lon, double lat, int year, int day);

			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_cols() = 0;
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
			Check if an output should be written to the output channel.

			@param out: Pointer to an output channel.
			*/
			bool should_write(OutputChannel *out);

			/*
			Assert that the output table is open. Terminate the simulation if not.
			*/
			void assert_table_open();
	};

	/*
	Base class for patch-level outputs.
	*/
	class StandLevelOutput : public OutputFile {
		public:
			/*
			Destructor.
			*/
			~StandLevelOutput();

			/*
			Create a new StandLevelOutput instance.

			@param param: Name of the paramter in the instruction file which determines the output file name.
			@param get_columns: A function which defines column definitions.
			@param desc: Description of the output file (please include units).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			StandLevelOutput(std::string param, std::vector<ColumnDefinition<Stand>> (*get_columns)(), std::string desc, bool sum, bool mean);

			/*
			Create a new StandLevelOutput instance.

			@param param: Name of the paramter in the instruction file which determines the output file name.
			@param columns: The columns in this output file.
			@param desc: Description of the output file (please include units).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			StandLevelOutput(std::string param, std::vector<ColumnDefinition<Stand>> columns, std::string desc, bool sum, bool mean);

			/*
			Write all rows of output to the output file for the current
			timestep. No default implementation is provided, because the
			implementation depends on the output - e.g. annual, daily, daily
			per-patch, etc.

			@param out: Pointer to an output channel.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel *channel, Gridcell &gridcell);

			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_cols();

		private:

			/*
			Write all required row for the specified patch.

			@param stand: The stand for which outputs should be written.
			*/
			std::vector<double> write_stand(const Stand &stand);

			/*
			Pointer to a function which generates column definitions to be used
			for this output file.

			@param patch: The patch.
			@param pft: The PFT.
			@return: The value which appears in each column of the output table.
			*/
			std::vector<ColumnDefinition<Stand>> (*get_column_definitions)();

			/*
			The column definitions. This will be empty until get_columns() is
			called, when the output file is opened.
			*/
			std::vector<ColumnDefinition<Stand>> columns;

			/*
			Iff true, a column containing the sum of all other columns will be
			written to the output file.
			*/
			bool write_sum;

			/*
			Iff true, a column containing the mean of all other columns will be
			written to the output file.
			*/
			bool write_mean;
	};

	/*
	Base class for patch-level outputs.

	This class produces one row of output for each patch in each stand.
	*/
	class PatchLevelOutput : public OutputFile {
		public:
			/*
			Destructor.
			*/
			virtual ~PatchLevelOutput() = 0;

			/*
			Write all rows of output to the output file for the current
			timestep. This will write a single row of output for each patch in
			each stand.

			@param out: Pointer to an output channel.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel *channel, Gridcell &gridcell);

			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_cols();

		protected:
			/*
			Create a new PatchLevelOutput instance. This is protected because
			PatchLevelOutput is abstract.

			@param param: Name of the paramter in the instruction file which determines the output file name.
			@param desc: Description of the output file (please include units).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			PatchLevelOutput(std::string param, std::string desc, bool sum, bool mean);

			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_columns() = 0;

			/*
			Write all required row for the specified patch.

			@param patch: The patch for which outputs should be written.
			@param out: Pointer to an output channel.
			*/
			virtual std::vector<double> write_patch(Patch &patch, OutputChannel *out) = 0;
		private:
			/*
			Iff true, a column containing the sum of all other columns will be
			written to the output file.
			*/
			bool write_sum;

			/*
			Iff true, a column containing the mean of all other columns will be
			written to the output file.
			*/
			bool write_mean;
	};

	/*
	This class outputs non-aggregated individual-level data.
	
	This will produce one row of output for each individual in each patch on
	each timestep. Use with caution.
	*/
	class IndividualOutput : public OutputFile {
		public:
			/*
			Create a new individual output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_columns: Pointer to a function which returns the columns to be used in this output file.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			IndividualOutput(std::string param_name, std::vector<ColumnDefinition<Individual>> (*get_columns)(), std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~IndividualOutput();

			/*
			Write all rows of output to the output file for the current
			timestep. This will write one row per individual per patch.

			@param out: Pointer to an output channel.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel *channel, Gridcell &gridcell);

		protected:
			/*
			Function which gets columns used by this output file.
			*/
			virtual ColumnDescriptors get_cols();

		private:

			/*
			Pointer to a function which generates column definitions to be used
			for this output file.

			@param patch: The patch.
			@param pft: The PFT.
			@return: The value which appears in each column of the output table.
			*/
			std::vector<ColumnDefinition<Individual>> (*get_column_definitions)();

			/*
			The column definitions. This will be empty until get_columns() is
			called, when the output file is opened.
			*/
			std::vector<ColumnDefinition<Individual>> columns;

			/*
			Iff true, a column containing the sum of all other columns will be
			written to the output file.
			*/
			bool write_sum;

			/*
			Iff true, a column containing the mean of all other columns will be
			written to the output file.
			*/
			bool write_mean;

			/*
			Name of the indiv ID column in individual output files.
			*/
			const char *COLNAME_INDIV_ID = "indiv";

			/*
			Name of the PFT column in individual output files.
			*/
			const char *COLNAME_PFT = "pft";
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
			*/
			IndivLevelArrayOutput(std::string param_name, std::string element_name, int (*get_nelem)(const Patch &), ArrayOutputReader<Individual> *get_output, std::string description);

			/*
			Destructor.
			*/
			~IndivLevelArrayOutput() { }

			/*
			Write all rows of output to the output file for the current
			timestep. One row will be written for each array element. One column
			will be written for each individual.

			@param out: Pointer to an output channel.
			@param gridcell: The current gridcell.
			*/
			virtual void write_rows(OutputChannel *channel, Gridcell &gridcell) override;
		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			ColumnDescriptors get_cols() override;

		private:
			/*
			Get the output values for given row of a particular patch.
			*/
			std::vector<double> write_row(const Patch &patch, int i, const std::vector<std::string> &pfts);

			/*
			Name of each element in the array (e.g. "layer").
			*/
			std::string elem_name;

			/*
			Get the number of elements in the output array.
			*/
			int (*get_n)(const Patch &);

			/*
			Read the i-th element of the output array.
			*/
			std::unique_ptr<ArrayOutputReader<Individual>> read_elem;
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
			@param col_name: Column name prefix. Will be suffixed with column index (0-n).
			@param n: Number of elements in the array.
			@param get_output: A function which returns the n-th element of the array from the patch.
			@param description: A description of the output (including units!!).
			@param spinup_output: True to write outputs during spinup. False otherwise.
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			PatchLevelArrayOutput(std::string param_name, std::string col_name, int n, ArrayOutputReader<Patch> *get_output, std::string description, bool sum, bool mean);

			/*
			Create a new patch-level array output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param col_name: Column name prefix. Will be suffixed with column index (0-n).
			@param n: Pointer to a function which gets the number of columns for a given patch.
			@param get_output: A function which returns the n-th element of the array from the patch.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			PatchLevelArrayOutput(std::string param_name, std::string col_name, int (*get_nelem)(), ArrayOutputReader<Patch> *get_output, std::string description, bool sum, bool mean);

			/*
			Create a new patch-level array output using predefined column names.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param col_names: Predefined column names.
			@param n: Number of elements in the array.
			@param get_output: A function which returns the n-th element of the array from the patch.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			PatchLevelArrayOutput(std::string param_name, std::vector<std::string> col_names, int n, ArrayOutputReader<Patch> *get_output, std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~PatchLevelArrayOutput();

			/*
			Write a row of outputs containing one cell per array element.

			@param patch: The patch for which outputs should be written.
			@param out: Pointer to an output channel.
			*/
			std::vector<double> write_patch(Patch &patch, OutputChannel *out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			ColumnDescriptors get_columns() override;

		private:
			/*
			Number of elements in the output array.
			*/
			int nelem;

			/*
			Vector of column names.
			*/
			std::vector<std::string> column_names;

			/*
			Pointer to a function which reads an element of the output array
			from the specified patch.

			@param patch: The patch.
			@param i: Array index in range [0, n)
			@return: The i-th element of the output array.
			*/
			std::unique_ptr<ArrayOutputReader<Patch>> read_output;

			/*
			Optional pointer to a function which gets the number of columns for
			a given patch. This may be NULL at runtime.
			*/
			int (*get_nelem)();
	};

	class PftOutput : public PatchLevelOutput {
		public:
			/*
			Create a new pft output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_output: A function which returns the output value from an individual.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			PftOutput(std::string param_name, double (*get_output)(const Patch &, const Pft &), std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~PftOutput();

			/*
			Write a row of outputs containing one cell per individual.

			@param patch: The patch for which outputs should be written.
			@param out: Pointer to an output channel.
			*/
			virtual std::vector<double> write_patch(Patch &patch, OutputChannel *out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			virtual ColumnDescriptors get_columns() override;

			/*
			Get the output value for a particular pft in a particular patch.

			@param patch: The patch.
			@param pft: The pft.
			*/
			virtual double get_value(const Patch &patch, const Pft &pft) const;

		private:

			/*
			Pointer to a function which reads an output from a pointer to an
			object of the template type.

			@param patch: The patch.
			@param pft: The PFT.
			@return: The value which appears in each column of the output table.
			*/
			double (*read_output)(const Patch &patch, const Pft &pft);
	};

	/*
	This class outputs one row per timestep and one column per pft. The output
	value is per-indiv but for each pft is aggregated across all individuals of
	that pft.
	*/
	class AggregatedIndivOutput : public PftOutput {
		public:
			/*
			Create a new aggregated indiv output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_output: A function which returns the output value from an individual.
			@param agg: Pointer to a function which takes an array, aggregates the array in some way and returns a scalar.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			AggregatedIndivOutput(std::string param_name, OutputReader<Individual> *get_output, double (*agg)(double *, int), std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~AggregatedIndivOutput() { }

		protected:

			/*
			Get the output value for a particular pft in a particular patch.

			@param patch: The patch.
			@param pft: The pft.
			*/
			double get_value(const Patch &patch, const Pft &pft) const override;

		private:

			/*
			Pointer to a function which reads an output from an individual. This
			will then be aggregated across all individuals for each pft.

			@param indiv: The individual.
			@return: The value which appears in each column of the output table.
			*/
			std::unique_ptr<OutputReader<Individual>> read_output;

			/*
			Pointer to a function which aggregates the output value across individuals.
			*/
			double (*aggregator)(double *, int);
	};

	class CustomPatchLevelOutput : public PatchLevelOutput {
		public:
			/*
			Create a new patch-level output with custom columns.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_columns: A function which defines column definitions.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			CustomPatchLevelOutput(std::string param_name, std::vector<ColumnDefinition<Patch>> (*get_columns)(), std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~CustomPatchLevelOutput();

			/*
			Write a row of outputs containing one cell per individual.

			@param patch: The patch for which outputs should be written.
			@param out: Pointer to an output channel.
			*/
			virtual std::vector<double> write_patch(Patch &patch, OutputChannel *out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			virtual ColumnDescriptors get_columns() override;

		private:

			/*
			Pointer to a function which generates column definitions to be used
			for this output file.

			@param patch: The patch.
			@param pft: The PFT.
			@return: The value which appears in each column of the output table.
			*/
			std::vector<ColumnDefinition<Patch>> (*get_column_definitions)();

			/*
			The column definitions. This will be empty until get_columns() is
			called, when the output file is opened.
			*/
			std::vector<ColumnDefinition<Patch>> columns;
	};

	class StandPftOutput : public PatchLevelOutput {
		public:
			/*
			Create a new standpft output.

			@param param_name: Name of the parameter in the instruction file which contains the output file name.
			@param get_output: A function which returns the output value from an standpft.
			@param description: A description of the output (including units!!).
			@param sum: Iff true, a total column will be written after all other columns.
			@param mean: Iff true, a mean column will be written after all other columns.
			*/
			StandPftOutput(std::string param_name, OutputReader<Standpft> *get_output, std::string description, bool sum, bool mean);

			/*
			Destructor.
			*/
			~StandPftOutput();

			/*
			Write a row of outputs containing one cell per standpft.

			@param patch: The patch for which outputs should be written.
			@param out: Pointer to an output channel.
			*/
			std::vector<double> write_patch(Patch &patch, OutputChannel *out);

		protected:
			/*
			Get column descriptors to be used for the output table.
			*/
			virtual ColumnDescriptors get_columns() override;

		private:

			/*
			Pointer to a function which reads an output from a pointer to an
			object of the template type.

			@param indiv: An individual.
			@return: The value which appears in each column of the output table.
			*/
			std::unique_ptr<OutputReader<Standpft>> read_output;
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
	Add an individual-level output which will, for each pft, be aggregated over
	all cohorts/individuals of that pft.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	@param aggregator: Function which aggregates the output value across individuals.
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_aggregated_indiv_output(std::string name, timestep_t timestep, double (*read_output)(const Individual&), std::string description, double (*aggregator)(double*, int), bool sum = false, bool mean = false);

	/*
	Add an individual-level output. This will produce an output file with one
	row per individual in the patch on each timestep.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	@param aggregator: Function which aggregates the output value across individuals.
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_indiv_output(std::string name, timestep_t timestep, std::vector<ColumnDefinition<Individual>> (*get_columns)(), std::string description, double (*aggregator)(double*, int), bool sum = false, bool mean = false);

	/*
	Add a pft-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param read_output: A function which returns a property of an individual to be written to each column.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_pft_output(std::string name, timestep_t timestep, double (*read_output)(const Patch&, const Pft&), std::string description, bool sum = false, bool mean = false);

	/*
	Add an individual-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param elem_name: Name of the index column (e.g. "layer", "timestep", etc).
	@param get_nelem: Get the number of array elements for a given patch.
	@param read_output: A function which returns the i-th element of the output property of an individual.
	@param description: Description of the output (with units!!).
	*/
	void add_indiv_array_output(std::string name, timestep_t timestep, std::string elem_name, int (*get_nelem)(const Patch &), double (*read_output)(const Individual&, int), std::string description);

	/*
	Add a patch-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_patch_output(std::string name, timestep_t timestep, std::string colname, double (*read_output)(const Patch&), std::string description, bool sum = false, bool mean = false);

	/*
	Add a patch-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param timestep_t: Timestep at which the output file should be written.
	@param get_columns: Pointer to a function which returns the columns used in this output file.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_patch_output(std::string name, timestep_t timestep, std::vector<ColumnDefinition<Patch>> (*get_columns)(), std::string description, bool sum = false, bool mean = false);

	/*
	Add a stand-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param read_output: A function which returns an output value for a given stand.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_stand_output(std::string name, timestep_t timestep, std::string colname, double (*read_output)(const Stand &), std::string description, bool sum = false, bool mean = false);

	/*
	Add a stand-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param get_columns: Pointer to a function which returns the columns used in this output file.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_stand_output(std::string name, timestep_t timestep, std::vector<ColumnDefinition<Stand>> (*get_columns)(), std::string description, bool sum = false, bool mean = false);

	/*
	Add a patch-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name prefix. Will be suffixed with column index.
	@param n: Number of array elements.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_patch_array_output(std::string name, timestep_t timestep, std::string colname, int nelem, double (*read_output)(const Patch&, int), std::string description, bool sum = false, bool mean = false);

	/*
	Add a patch-level array output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name prefix. Will be suffixed with column index.
	@param get_nelem: Pointer to a function which gets the number of columns for the given patch.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_patch_array_output(std::string name, timestep_t timestep, std::string colname, int (*get_nelem)(), double (*read_output)(const Patch&, int), std::string description, bool sum = false, bool mean = false);

	/*
	Add a patch-level array output with a predetermined list of column names.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column names.
	@param read_output: A function which returns the i-th element of the output array.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_patch_array_output(std::string name, timestep_t timestep, std::vector<std::string> colnames, double (*read_output)(const Patch&, int), std::string description, bool sum = false, bool mean = false);

	/*
	Add a standpft-level output.

	@param name: Parameter in the instruction file which determines output file name.
	@param colname: Column name in the output file.
	@param read_output: A function which returns the output value for a particular standpft.
	@param description: Description of the output (with units!!).
	@param sum: Iff true, a total column will be written after all other columns.
	@param mean: Iff true, a mean column will be written after all other columns.
	*/
	void add_spft_output(std::string name, timestep_t timestep, double (*read_output)(const Standpft&), std::string description, bool sum = false, bool mean = false);

	/*
	Add an output file at the appropriate timestep.

	@param file: The output file to be registered.
	@param timestep: The timestep at which the file should be written.
	*/
	void register_output(OutputFile *file, timestep_t timestep);

	/*
	Get the vector of output files which corresponds to a particular timestep.

	@param timestep: The timestep.
	*/
	std::vector<std::unique_ptr<OutputFile>> &get_outputs_vector(timestep_t timestep);

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
	Total time spent in dave output code.
	*/
	double total_time;
};

/*
True iff outputs should be written during spinup.
*/
extern bool spinup_output;

/*
Iff true, mean/std dev across patches will be written instead of raw outputs.
*/
extern bool patch_output_stats;

}

#endif // LPJ_GUESS_DAVE_OUTPUT_H
