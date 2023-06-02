///////////////////////////////////////////////////////////////////////////////////////
/// \file daveoutput.cpp
/// \brief Implementation of the DAVE output module
///
/// \author Drew Holzworth
/// $Date: 2023-02-17 10:58:54 +1100 (Fri, 17 Feb 2023) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <limits.h>

#include "config.h"
#include "commonoutput.h"
#include "daveoutput.h"
#include "daveoutputhelpers.h"
#include "gutil.h"
#include "outputchannel.h"
#include "outputmodule.h"
#include "parameters.h"
#include "util.h"

namespace GuessOutput {

// Register the DAVE output module.
REGISTER_OUTPUT_MODULE("dave", DaveOutput)

// #############################################################################
// DaveOutput implementation.
// #############################################################################

DaveOutput::DaveOutput() {
	double t_start = get_wall_time();

	// Annual per-patch individual-level outputs.
	add_annual_indiv_output("file_dave_alai", get_lai, "Annual LAI(m2/m2)");
	add_annual_indiv_output("file_dave_afpc", [](Individual &i) { return i.fpc; }, "Annual FPC ()");

	// Daily per-patch individual-level outputs.
	add_daily_indiv_output("file_dave_lai", [](Individual &i) { return i.lai; }, "LAI by individual (m2/m2)");
	add_daily_indiv_output("file_dave_sla", [](Individual &i) { return i.sla(); }, "SLA by individual (m2/kgC)");
	add_daily_indiv_output("file_dave_cmass_leaf", [](Individual &i) { return i.get_cmass_leaf(); }, "Leaf C mass (kgC/m2)");
	add_daily_indiv_output("file_dave_nmass_leaf", [](Individual &i) { return i.get_nmass_leaf(); }, "Leaf N mass (kgN/m2)");
	add_daily_indiv_output("file_dave_cmass_root", [](Individual &i) { return i.cmass_root; }, "Root C mass (kgC/m2)");
	add_daily_indiv_output("file_dave_cmass", [](Individual &i) { return i.ccont(); }, "Total C biomass per individual per patch (kgC/m2)");
	add_daily_indiv_output("file_dave_nmass", [](Individual &i) { return i.ncont(); }, "Total N biomass per individual per patch (kgN/m2)");
	add_daily_indiv_output("file_dave_cmass_storage", [](Individual &i) { return i.cmass_storage; }, "Storage C mass (kgC/m2)");
	add_daily_indiv_output("file_dave_height", [](Individual &i) { return i.height; }, "Plant height (m)");
	add_daily_indiv_output("file_dave_diameter", [](Individual &i) { return i.diam; }, "Stem diameter (m)");
	add_daily_indiv_output("file_dave_density", [](Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)");
	add_daily_indiv_output("file_dave_fpar", [](Individual &i) { return i.fpar; }, "Daily fpar [0-1]");
	add_daily_indiv_output("file_dave_indiv_gpp", [](Individual &i) { return i.dgpp * G_PER_KG; }, "Daily GPP (gC/m2/day)");
	add_daily_indiv_output("file_dave_indiv_resp", [](Individual &i) { return (i.dgpp - i.dnpp) * G_PER_KG; }, "Daily growth respiration (gC/m2/day)");
	add_daily_indiv_output("file_dave_indiv_npp", [](Individual &i) { return i.dnpp * G_PER_KG; }, "Daily NPP (gC/m2/day)");
	add_daily_indiv_output("file_dave_fpc", [](Individual &i) { return i.fpc; }, "Daily FPC ()");
	add_daily_indiv_output("file_dave_crownarea", [](Individual &i) { return i.crownarea; }, "Daily Crown Area (m2)");
	// fixme: wscal should really be a ppft-level output. Which is still TBI...
	add_daily_indiv_output("file_dave_wscal", [](Individual &i) { return i.patchpft().wscal; }, "Water scalar [0-1], 1 = min stress");
	add_daily_indiv_output("file_dave_agd_g", [](Individual &i) { return i.photosynthesis.agd_g; }, "Daily gross photosynthesis (gC/m2/day)");
	add_daily_indiv_output("file_dave_rd_g", [](Individual &i) { return i.photosynthesis.rd_g; }, "Daily leaf respiration (gC/m2/day)");
	add_daily_indiv_output("file_dave_je", [](Individual &i) { return i.photosynthesis.je; }, "Daily PAR-limited photosynthesis rate (gC/m2/h)");
	add_daily_indiv_output("file_dave_vm", [](Individual &i) { return i.photosynthesis.vm; }, "RuBisCO capacity (gC/m2/day)");
	add_daily_indiv_output("file_dave_growth_stage", [](Individual &i) { return i.growth_stage; }, "Growth Stage (0-1)");
	add_daily_indiv_output("file_dave_growth_demand", [](Individual &i) { return i.growth_demand; }, "Growth Demand (0-1)");

	// Daily individual-level array outputs.
	add_daily_indiv_array_output("file_dave_layerwise_fpar", "layer", get_nlayer_fpar, [](Individual &indiv, int i) { return indiv.fpar_layerwise[i]; }, "Daily layered fpar [0-1]");
	add_daily_indiv_array_output("file_dave_layerwise_lai", "layer", get_nlayer_lai, [](Individual &indiv, int i) { return indiv.lai_layerwise[i]; }, "Daily layered LAI (m2/m2)");

	// Subdaily individual-level outputs.
	add_subdaily_indiv_output("file_dave_subdaily_an", [](Individual &indiv, int i) { return indiv.ps_subdaily.An[i]; }, "Subdaily net photosynthesis (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_rd", [](Individual &indiv, int i) { return indiv.ps_subdaily.Rd[i]; }, "Subdaily leaf respiration (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_anc", [](Individual &indiv, int i) { return indiv.ps_subdaily.Anc[i]; }, "Subdaily RuBisCO-limited net photosynthesis (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_anj", [](Individual &indiv, int i) { return indiv.ps_subdaily.Anj[i]; }, "Subdaily RuBP regeneration-limited photoysnthesis (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_gsw", [](Individual &indiv, int i) { return indiv.ps_subdaily.gsw[i]; }, "Subdaily stomatal conductance to water vapour (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_ci", [](Individual &indiv, int i) { return indiv.ps_subdaily.ci[i]; }, "Subdaily intercellular CO2 concentration (mol/mol)");
	add_subdaily_indiv_output("file_dave_subdaily_vcmax", [](Individual &indiv, int i) { return indiv.ps_subdaily.vcmax[i]; }, "Subdaily max carboxylation rate (mol/m2/s)");
	add_subdaily_indiv_output("file_dave_subdaily_jmax", [](Individual &indiv, int i) { return indiv.ps_subdaily.jmax[i]; }, "Subdaily max electron transport rate (mol/m2/s)");

	// Annual patch-level outputs.
	// ...

	// Daily patch-level outputs.
	add_daily_patch_level_array_output("file_dave_soil_nmass_avail", "nmass", 1, [](Patch &p, int i) { return p.soil.nmass_avail() * min(1.0, p.fpc_total); }, "Available Soil N for plant uptake", false);
	add_daily_patch_level_array_output("file_dave_sw", "sw_", NSOILLAYER, [](Patch &p, int i) { return p.soil.get_layer_soil_water(i); }, "Soil water content as fraction of available water holding capacity", false);
	add_daily_patch_level_array_output("file_dave_swmm", "swmm_", NSOILLAYER, [](Patch &p, int i) { return p.soil.get_layer_soil_water(i) * p.soil.soiltype.wsats[i]; }, "Soil water content (mm)");
	add_daily_patch_level_array_output("file_dave_swvol", "swvol_", NSOILLAYER, [](Patch &p, int i) { return p.soil.get_layer_soil_water(i) * p.soil.soiltype.wsats[i] / Dz_soil; }, "Volumetric soil water content (m3/m3)", false);
	// add_daily_patch_level_array_output("file_dave_cfluxes_patch", get_patch_flux_names(), (int)Fluxes::NPERPATCHFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPatchFluxType)i, date.day) * G_PER_KG; }, "Daily patch-level fluxes (g/m2/day)");
	// add_daily_patch_level_array_output("file_dave_cfluxes_pft", get_pft_flux_names(), (int)Fluxes::NPERPFTFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPFTFluxType)i, date.day) * G_PER_KG; }, "Daily PFT-level fluxes (g/m2/day)");
	// -1 * NEE + ER = -(rh - npp) + ER = -(rh - (gpp - ra)) + ER = -rh + gpp - ra + ER
	// npp = gpp - ra
	// nee = rh - npp // todo: annual considers fire + est fluxes
	// gpp = patch.fluxes.get_daily_flux(Fluxes::GPP, date.day)
	// ra = patch.fluxes.get_daily_flux(Fluxes::RA, date.day)
	// rh = patch.fluxes.get_daily_flux(Fluxes::SOILC, date.day)
	add_daily_patch_level_array_output("file_dave_resp", "resp", 1, get_er, "Daily patch-level ecosystem respiration (gC/m2/day)");
	add_daily_patch_level_array_output("file_dave_gpp", "gpp", 1, get_gpp, "Daily patch-level GPP (gC/m2/day)");
	add_daily_patch_level_array_output("file_dave_npp", "npp", 1, get_npp, "Daily patch-level NPP (gC/m2/day)");
	add_daily_patch_level_array_output("file_dave_nee", "nee", 1, get_nee, "Daily NEE (gC/m2/day)");
	add_daily_patch_level_array_output("file_dave_et", "et", 1, get_et, "Daily patch-level evapotranspiration (mm/day)");

	// Subdaily patch-level outputs.
	add_subdaily_patch_level_output("file_dave_met_temp", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.temps[i]; }, "Subdaily temperature (°C)");
	add_subdaily_patch_level_output("file_dave_met_par", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.pars[i]; }, "Subdaily PAR (J/m2/timestep)");
	add_subdaily_patch_level_output("file_dave_met_vpd", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.vpds[i]; }, "Subdaily VPD (kPa)");
	add_subdaily_patch_level_output("file_dave_met_insol", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.insols[i]; }, "Subdaily insolation (units depend on instype)");
	add_subdaily_patch_level_output("file_dave_met_precip", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.precs[i]; }, "Subdaily precipitation (mm)");
	add_subdaily_patch_level_output("file_dave_met_pressure", [](Patch &patch, int i) { return patch.stand.get_gridcell().climate.pressure[i]; }, "Subdaily air pressure (Pa)");

	// Register instruction file parameters.
	declare_parameters();

	total_time = get_wall_time() - t_start;
}

DaveOutput::~DaveOutput() {
	printf("Total time spent writing dave outputs: %.2fs\n", total_time);
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
	if (spinup_output || date.year >= nyear_spinup) {
		write_outputs(gridcell, daily_outputs);
	}
	total_time += get_wall_time() - t0;
}

void DaveOutput::add_subdaily_indiv_output(std::string name, double (*read_output)(Individual&, int), std::string description) {
	add_daily_indiv_array_output(name, COLNAME_TIMESTEP, [](Patch &p) { return date.subdaily; }, read_output, description);
}

void DaveOutput::add_daily_indiv_array_output(std::string name, std::string elem_name, int (*get_nelem)(Patch &), double (*read_output)(Individual&, int), std::string description) {
	daily_outputs.push_back(create_indiv_array_output(name, elem_name, get_nelem, read_output, description));
}

void DaveOutput::add_annual_indiv_array_output(std::string name, std::string elem_name, int (*get_nelem)(Patch &), double (*read_output)(Individual&, int), std::string description) {
	annual_outputs.push_back(create_indiv_array_output(name, elem_name, get_nelem, read_output, description));
}

void DaveOutput::add_subdaily_patch_level_output(std::string name, double (*read_output)(Patch&, int), std::string description) {
	add_daily_patch_level_array_output(name, COLNAME_TIMESTEP, date.subdaily, read_output, description);
}

void DaveOutput::add_daily_patch_level_array_output(std::string name, std::string colname, int nelem, double (*read_output)(Patch&, int), std::string description, bool total_is_sum) {
	daily_outputs.push_back(create_array_output(name, colname, nelem, read_output, description, total_is_sum));
}

void DaveOutput::add_daily_patch_level_array_output(std::string name, std::vector<std::string> colnames, int nelem, double (*read_output)(Patch&, int), std::string description, bool total_is_sum) {
	daily_outputs.push_back(create_array_output(name, colnames, nelem, read_output, description, total_is_sum));
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
	declare_parameter("spinup_output", &spinup_output, "True to write DAVE outputs during the spinup, false otherwise.");

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
	std::string description,
	bool total_is_sum) {

	return std::make_unique<PatchLevelArrayOutput>(name, colname, n, read_output, description, total_is_sum);
}

std::unique_ptr<DaveOutput::OutputFile> DaveOutput::create_array_output(
	std::string name,
	std::vector<std::string> colnames,
	int n,
	double (*read_output)(Patch&, int),
	std::string description,
	bool total_is_sum) {

	return std::make_unique<PatchLevelArrayOutput>(name, colnames, n, read_output, description, total_is_sum);
}

std::unique_ptr<DaveOutput::OutputFile> DaveOutput::create_indiv_array_output(
		std::string name,
		std::string elem_name,
		int (*get_nelem)(Patch &),
		double (*read_output)(Individual&, int),
		std::string description) {
	return std::make_unique<IndivLevelArrayOutput>(name, elem_name, get_nelem, read_output, description);
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
	if (date.year < nyear_spinup && !spinup_output) {
		return false;
	}
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
	std::string description)
		: PatchLevelOutput(param_name, description) {

	read_output = get_output;
}

DaveOutput::IndividualOutput::~IndividualOutput() { }

void DaveOutput::IndividualOutput::write_patch(Patch &patch, OutputRows &out) {
	Vegetation &vegetation = patch.vegetation;

	std::vector<std::string> pfts = get_pfts(false);

	// Running total.
	double total = 0;

	// Write one column per PFT - note that there will not necessarily be a 1:1
	// mapping of PFTs to individuals, as some PFTs may not be established.
	vegetation.firstobj();
	for (int i = 0; i < pfts.size(); i++) {
		std::vector<Individual *> indivs = get_indivs(vegetation, pfts[i].c_str());

		double result = 0;
		if (indivs.size() > 0) {
			for (Individual *indiv : indivs) {
				result += indiv ? read_output(*indiv) : 0;
			}
			// result /= indivs.size();
		}
		write_output(out, result);
		total += result;

		vegetation.nextobj();
	}

	// Write total.
	write_output(out, total);
}

ColumnDescriptors DaveOutput::IndividualOutput::get_columns() {
	ColumnDescriptors columns;
	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	columns += ColumnDescriptors(get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);

	return columns;
}

// #############################################################################
// PatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int n,
		double (*get_output)(Patch &, int),
		std::string description,
		bool sum) : PatchLevelOutput(param_name, description) {

	nelem = n;
	read_output = get_output;
	total_is_sum = sum;

	for (int i = 0; i < n; i++) {
		if (n == 1) {
			column_names.push_back(col_name);
		} else {
			column_names.push_back(col_name + std::to_string(i));
		}
	}
}

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(std::string param_name, std::vector<std::string> col_names, int n, double (*get_output)(Patch &, int), std::string description, bool sum)
 : PatchLevelOutput(param_name, description) {
	column_names = col_names;
	nelem = n;
	read_output = get_output;
	total_is_sum = sum;
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
	// If total is not a sum, it is a mean. See documentation of total_is_sum.
	if (!total_is_sum) {
		total /= nelem;
	}
	write_output(out, total);
}

ColumnDescriptors DaveOutput::PatchLevelArrayOutput::get_columns() {
	ColumnDescriptors columns;

	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	for (int i = 0; i < nelem; i++) {
		std::string name = column_names[i];
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		columns += ColumnDescriptor(name.c_str(), width, NUMERIC_PRECISION);
	}
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// IndivLevelArrayOutput implementation.
// #############################################################################

DaveOutput::IndivLevelArrayOutput::IndivLevelArrayOutput(
		std::string param_name,
		std::string element_name,
		int (*get_nelem)(Patch &),
		double (*get_output)(Individual &, int),
		std::string description)
			: OutputFile(param_name, description)
{
	elem_name = element_name;
	get_n = get_nelem;
	read_elem = get_output;
}

void DaveOutput::IndivLevelArrayOutput::write_rows(OutputChannel &channel, Gridcell &gridcell) {
	double lon = gridcell.get_lon();
	double lat = gridcell.get_lat();
	int year = date.get_calendar_year();
	int day = date.day;
	OutputRows rows(&channel, lon, lat, year, day);

	std::vector<std::string> pfts = get_pfts(false);

	GuessContainer<Stand>::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {

		Stand &stand = *gc_itr;
		stand.firstobj();

		while (stand.isobj) {
			Patch &patch = stand.getobj();

			// Write patch ID.
			write_output(rows, patch.id);

			int nelem = get_n(patch);
			Vegetation &vegetation = patch.vegetation;
			for (int i = 0; i < nelem; i++) {

				// Write array index.
				write_output(rows, i);

				vegetation.firstobj();
				double total = 0;

				for (int i = 0; i < pfts.size(); i++) {
					std::vector<Individual *> indivs = get_indivs(vegetation, pfts[i].c_str());
						double indiv_total = 0;
					for (Individual *indiv : indivs) {
						// Read array element for this individual.
						indiv_total += read_elem(*indiv, i);
					}
					total += indiv_total;

					// Write current array element for this individual.
					write_output(rows, indiv_total);
				}

				// Write total value across all individuals for this array elem.
				write_output(rows, total);

				// Finish the row, so next array element is written to next row.
				finish_row(channel, lon, lat, year, day);
			}
		}
	}
}

ColumnDescriptors DaveOutput::IndivLevelArrayOutput::get_columns() {
	ColumnDescriptors columns;

	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	columns += ColumnDescriptor(elem_name.c_str(), COLUMN_WIDTH, 0);
	columns += ColumnDescriptors(get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);

	return columns;
}

bool spinup_output;

} // namespace GuessOutput
