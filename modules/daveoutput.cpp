///////////////////////////////////////////////////////////////////////////////////////
/// \file daveoutput.cpp
/// \brief Implementation of the DAVE output module
///
/// \author Drew Holzworth
/// $Date: 2023-02-17 10:58:54 +1100 (Fri, 17 Feb 2023) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <limits.h>

#include "commonoutput.h"
#include "daveoutput.h"
#include "gutil.h"
#include "outputchannel.h"
#include "outputmodule.h"
#include "parameters.h"
#include "util.h"

namespace GuessOutput {

// Register the DAVE output module.
REGISTER_OUTPUT_MODULE("dave", DaveOutput)

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

// #############################################################################
// DaveOutput implementation.
// #############################################################################

DaveOutput::DaveOutput() {
	double t_start = get_wall_time();

	// Annual individual-level outputs.
	add_annual_indiv_output("file_dave_alai", get_lai, "Annual LAI(m2/m2)");

	// Daily individiaul-level outputs.
	add_daily_indiv_output("file_dave_lai", [](Individual &i) { return i.lai; }, "LAI by individual (m2/m2)");
	add_daily_indiv_output("file_dave_sla", [](Individual &i) { return i.sla(); }, "SLA by individual (m2/kgC)");
	add_daily_indiv_output("file_dave_cmass_leaf", [](Individual &i) { return i.get_cmass_leaf(); }, "Leaf C mass (kgC/m2)");
	add_daily_indiv_output("file_dave_nmass_leaf", [](Individual &i) { return i.get_nmass_leaf(); }, "N mass available (kgN/m2)");
	add_daily_indiv_output("file_dave_height", [](Individual &i) { return i.height; }, "Plant height (m)");
	add_daily_indiv_output("file_dave_diameter", [](Individual &i) { return i.diam; }, "Stem diameter (m)");
	add_daily_indiv_output("file_dave_density", [](Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)");
	add_daily_indiv_output("file_dave_fpar", [](Individual &i) { return i.fpar; }, "Daily fpar [0-1]");
	add_daily_indiv_output("file_dave_assim", [](Individual &i) { return i.dgpp; }, "Daily net C assimilation (gross photosynthesis minus leaf respiration (kgC/m2/day)");
	add_daily_indiv_output("file_dave_resp", [](Individual &i) { return i.dgpp - i.dnpp; }, "Daily growth respiration (kgC/m2/day)");
	add_daily_indiv_output("file_dave_npp", [](Individual &i) { return i.dnpp; }, "Daily NPP (kgC/m2/day)");

	// Annual patch-level outputs.
	// ...

	// Daily patch-level outputs.
	add_daily_patch_level_array_output("file_dave_sw", "sw", NSOILLAYER, [](Patch &p, int i) { return p.soil.get_layer_soil_water(i); }, "Soil water content as fraction of available water holding capacity");

	// Register instruction file parameters.
	declare_parameters();

	total_time = get_wall_time() - t_start;
}

DaveOutput::~DaveOutput() {
	printf("Total time spent in dave input: %.2fs\n", total_time);
}

void DaveOutput::init() {
	double t_start = get_wall_time();
	open_outfiles();
	total_time += get_wall_time() - t_start;
}

void DaveOutput::write_outputs(Gridcell &gridcell, std::vector<std::unique_ptr<OutputFile>> &outputs) {
	ensure_channel_open();

	for (std::unique_ptr<OutputFile> &output : outputs) {
		if (output->enabled()) {
			output->write_rows(*output_channel, gridcell);
		}
	}
}

void DaveOutput::outannual(Gridcell &gridcell) {
	double t0 = get_wall_time();
	write_outputs(gridcell, annual_outputs);
	total_time += get_wall_time() - t0;
}

void DaveOutput::outdaily(Gridcell &gridcell) {
	double t0 = get_wall_time();
	write_outputs(gridcell, daily_outputs);
	total_time += get_wall_time() - t0;
}

void DaveOutput::add_daily_patch_level_array_output(std::string name, std::string colname, int nelem, double (*read_output)(Patch&, int), std::string description) {
	daily_outputs.push_back(create_array_output(name, colname, nelem, read_output, description));
}

void DaveOutput::add_annual_patch_level_array_output(std::string name, std::string colname, int nelem, double (*read_output)(Patch&, int), std::string description) {
	annual_outputs.push_back(create_array_output(name, colname, nelem, read_output, description));
}

void DaveOutput::add_daily_indiv_output(std::string name, double (*read_output)(Individual&), std::string description) {
	daily_outputs.push_back(create_individual_output(name, read_output, description));
}

void DaveOutput::add_annual_indiv_output(std::string name, double (*read_output)(Individual&), std::string description) {
	annual_outputs.push_back(create_individual_output(name, read_output, description));
}

void DaveOutput::declare_parameters() {
	declare_parameters(daily_outputs);
	declare_parameters(annual_outputs);
}

void DaveOutput::declare_parameters(std::vector<std::unique_ptr<OutputFile>> &outputs) {
	for (auto &file : outputs) {
		file->declare_parameters();
	}
}

void DaveOutput::open_outfiles() {
	open_outfiles(daily_outputs);
	open_outfiles(annual_outputs);
}

void DaveOutput::open_outfiles(std::vector<std::unique_ptr<OutputFile>> &outputs) {
	ensure_channel_open();

	for (std::unique_ptr<OutputFile> &output : outputs) {
		if (output->enabled()) {
			output->open(*output_channel);
		}
	}
}

std::unique_ptr<DaveOutput::OutputFile> DaveOutput::create_individual_output(
	std::string name,
	double (*read_output)(Individual&),
	std::string description) {

	return std::make_unique<IndividualOutput>(name, read_output, description);
}

std::unique_ptr<DaveOutput::OutputFile> DaveOutput::create_array_output(
	std::string name,
	std::string colname,
	int n,
	double (*read_output)(Patch&, int),
	std::string description) {

	return std::make_unique<PatchLevelArrayOutput>(name, colname, n, read_output, description);
}

// #############################################################################
// DaveOutput::OutputFile implementation.
// #############################################################################

DaveOutput::OutputFile::~OutputFile() {
	// table will be NULL if this output file is disabled.
	if (table) {
		delete table;
	}
	table = NULL;
}

DaveOutput::OutputFile::OutputFile(std::string param, std::string desc) {
	param_name = param;
	description = desc;
	table = NULL;
}

void DaveOutput::OutputFile::declare_parameters() {
	declare_parameter(param_name.c_str(), &file_name, PATH_MAX, description.c_str());
}

void DaveOutput::OutputFile::open(OutputChannel &channel) {
	const char *name = file_name.c_str();
	Table tbl = channel.create_table(TableDescriptor(name, get_columns()));
	table = new Table(tbl.id());
}

void DaveOutput::OutputFile::close(OutputChannel &channel) {
	assert_table_open();
	channel.close_table(*table);
}

bool DaveOutput::OutputFile::should_write(OutputRows &rows) {
	// if (date.year < nyear_spinup) {
	// 	return false;
	// }
	return true;
}

void DaveOutput::OutputFile::write_output(OutputRows &rows, double value)
{
	assert_table_open();
	if (should_write(rows)) {
		rows.add_value(*table, value);
	}
}

void DaveOutput::OutputFile::write_output(OutputRows &rows, int value)
{
	assert_table_open();
	if (should_write(rows)) {
		rows.add_value(*table, value);
	}
}

void DaveOutput::OutputFile::finish_row(OutputChannel &rows, double lon, double lat, int year, int day) {
	assert_table_open();
	rows.finish_row(*table, lon, lat, year, day);
}

void DaveOutput::OutputFile::assert_table_open() {
	if (!table) {
		fail("Output table for '%s' is not open", param_name.c_str());
	}
}

bool DaveOutput::OutputFile::enabled() {
	// Output file is enabled iff output file name is not empty.
	return !file_name.empty();
}

// #############################################################################
// PatchLevelOutput implementation.
// #############################################################################

DaveOutput::PatchLevelOutput::~PatchLevelOutput() { }

DaveOutput::PatchLevelOutput::PatchLevelOutput(std::string param, std::string desc)
: OutputFile(param, desc) { }

void DaveOutput::PatchLevelOutput::write_rows(OutputChannel &channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;
		stand.firstobj();

		while (stand.isobj) {
			Patch &patch = stand.getobj();

			// Write patch ID.
			write_output(rows, patch.id);

			write_patch(patch, rows);

			finish_row(channel, lon, lat, year, day);

			stand.nextobj();
		}

		gc_itr++;
	}
}

// #############################################################################
// IndividualOutput implementation.
// #############################################################################

DaveOutput::IndividualOutput::IndividualOutput(
	std::string param_name,
	double (*get_output)(Individual &),
	std::string description) : PatchLevelOutput(param_name, description) {

	read_output = get_output;
}

DaveOutput::IndividualOutput::~IndividualOutput() { }

void DaveOutput::IndividualOutput::write_patch(Patch &patch, OutputRows &out) {
	Vegetation &vegetation = patch.vegetation;

	std::vector<std::string> pfts = get_pfts();

	// Running total.
	double total = 0;

	// Write one column per PFT - note that there will not necessarily be a 1:1
	// mapping of PFTs to individuals, as some PFTs may not be established.
	vegetation.firstobj();
	for (int i = 0; i < pfts.size(); i++) {
		Individual *indiv = get_indiv(vegetation, pfts[i].c_str());

		double value = indiv ? read_output(*indiv) : 0;
		write_output(out, value);
		total += value;

		vegetation.nextobj();
	}

	// Write total.
	write_output(out, total);
}

ColumnDescriptors DaveOutput::IndividualOutput::get_columns() {
	ColumnDescriptors columns;
	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	columns += ColumnDescriptors(get_pfts(), COLUMN_WIDTH, NUMERIC_PRECISION);
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);

	return columns;
}

std::vector<std::string> DaveOutput::IndividualOutput::get_pfts() {
	std::vector<std::string> pfts;
	pftlist.firstobj();
	while (pftlist.isobj) {
		pfts.push_back((char *)pftlist.getobj().name);
		pftlist.nextobj();
	}
	return pfts;
}

// #############################################################################
// PatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int n,
		double (*get_output)(Patch &, int),
		std::string description) : PatchLevelOutput(param_name, description) {
	column_name = col_name;
	nelem = n;
	read_output = get_output;
}

DaveOutput::PatchLevelArrayOutput::~PatchLevelArrayOutput() { }

void DaveOutput::PatchLevelArrayOutput::write_patch(
		Patch &patch, OutputRows &out) {
	double total = 0;
	for (int i = 0; i < nelem; i++) {
		double value = read_output(patch, i);
		total += value;
		write_output(out, value);
	}
	write_output(out, total);
}

ColumnDescriptors DaveOutput::PatchLevelArrayOutput::get_columns() {
	ColumnDescriptors columns;

	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	for (int i = 0; i < nelem; i++) {
		std::string name = column_name + std::to_string(i);
		columns += ColumnDescriptor(name.c_str(), COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

} // namespace GuessOutput
