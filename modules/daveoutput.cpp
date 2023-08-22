///////////////////////////////////////////////////////////////////////////////////////
/// \file daveoutput.cpp
/// \brief Implementation of the DAVE output module
///
/// \author Drew Holzworth
/// $Date: 2023-02-17 10:58:54 +1100 (Fri, 17 Feb 2023) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include <limits.h>
#include <sstream>

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
	add_indiv_output("file_dave_alai", ANNUAL, get_lai, "Annual LAI(m2/m2)", sum, true);
	add_indiv_output("file_dave_afpc", ANNUAL, [](const Individual &i) { return i.fpc; }, "Annual FPC ()", mean, false, true);
	add_indiv_output("file_dave_acmass", ANNUAL, [](const Individual &i) { return i.ccont(); }, "Annual cmass (kgC/m2)", sum, true);
	add_indiv_output("file_dave_anmass", ANNUAL, [](const Individual &i) { return i.ncont(); }, "Annual nmass (kgN/m2)", sum, true);
	add_indiv_output("file_dave_aheight", ANNUAL, [](const Individual &i) { return i.height; }, "Annual plant height (m)", mean, false, true);
	add_indiv_output("file_dave_aaet", ANNUAL, [](const Individual &i) { return i.aaet; }, "Annual actual evapotranspiration (mm)", sum, true);
	add_indiv_output("file_dave_adensity", ANNUAL, [](const Individual &i) { return i.densindiv; }, "Annual density of individuals over patch (/m2)", sum, true);

	// Annual per-patch pft-level outputs.
	add_pft_output("file_dave_anpp", ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id); }, "Total individual-level Annual NPP (kgC/m2/year)");
	add_pft_output("file_dave_agpp", ANNUAL, [](const Patch &patch, const Pft &pft) { return patch.fluxes.get_annual_flux(Fluxes::GPP, pft.id); }, "Total individual-level Annual GPP (kgC/m2/year)");
	add_pft_output("file_dave_aclitter", ANNUAL, get_pft_cmass_litter, "Total individual-level Annual GPP (kgC/m2/year)");

	// Daily per-patch individual-level outputs.
	add_indiv_output("file_dave_lai", DAILY, [](const Individual &i) { return i.lai; }, "LAI by individual (m2/m2)", sum, true);
	add_indiv_output("file_dave_sla", DAILY, [](const Individual &i) { return i.sla(); }, "SLA by individual (m2/kgC)", mean, false, true);
	add_indiv_output("file_dave_cmass_leaf", DAILY, [](const Individual &i) { return i.get_cmass_leaf(); }, "Leaf C mass (kgC/m2)", sum, true);
	add_indiv_output("file_dave_nmass_leaf", DAILY, [](const Individual &i) { return i.get_nmass_leaf(); }, "Leaf N mass (kgN/m2)", sum, true);
	add_indiv_output("file_dave_cmass_root", DAILY, [](const Individual &i) { return i.cmass_root; }, "Root C mass (kgC/m2)", sum, true);
	add_indiv_output("file_dave_cmass", DAILY, [](const Individual &i) { return i.ccont(); }, "Total C biomass per individual per patch (kgC/m2)", sum, true);
	add_indiv_output("file_dave_nmass", DAILY, [](const Individual &i) { return i.ncont(); }, "Total N biomass per individual per patch (kgN/m2)", sum, true);
	add_indiv_output("file_dave_cmass_storage", DAILY, [](const Individual &i) { return i.cmass_storage; }, "Storage C mass (kgC/m2)", sum, true);
	add_indiv_output("file_dave_height", DAILY, [](const Individual &i) { return i.height; }, "Plant height (m)", mean, false, true);
	add_indiv_output("file_dave_diameter", DAILY, [](const Individual &i) { return i.diam; }, "Stem diameter (m)", mean, false, true);
	add_indiv_output("file_dave_density", DAILY, [](const Individual &i) { return i.densindiv; }, "Density of individuals over patch (/m2)", sum, true);
	add_indiv_output("file_dave_fpar", DAILY, [](const Individual &i) { return i.fpar; }, "Daily fpar [0-1]", sum, true);
	add_indiv_output("file_dave_indiv_gpp", DAILY, [](const Individual &i) { return i.dgpp * G_PER_KG; }, "Daily GPP (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_indiv_resp", DAILY, [](const Individual &i) { return (i.dgpp - i.dnpp) * G_PER_KG; }, "Daily growth respiration (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_indiv_npp", DAILY, [](const Individual &i) { return i.dnpp * G_PER_KG; }, "Daily NPP (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_fpc", DAILY, [](const Individual &i) { return i.fpc; }, "Daily FPC ()", mean, false, true);
	add_indiv_output("file_dave_crownarea", DAILY, [](const Individual &i) { return i.crownarea; }, "Daily Crown Area (m2)", sum, true);
	// fixme: wscal should really be a ppft-level output. Which is still TBI...
	add_indiv_output("file_dave_wscal", DAILY, [](const Individual &i) { return i.patchpft().wscal; }, "Water scalar [0-1], 1 = min stress", mean, false, true);
	add_indiv_output("file_dave_agd_g", DAILY, [](const Individual &i) { return i.photosynthesis.agd_g; }, "Daily gross photosynthesis (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_rd_g", DAILY, [](const Individual &i) { return i.photosynthesis.rd_g; }, "Daily leaf respiration (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_je", DAILY, [](const Individual &i) { return i.photosynthesis.je; }, "Daily PAR-limited photosynthesis rate (gC/m2/h)", sum, true);
	add_indiv_output("file_dave_vm", DAILY, [](const Individual &i) { return i.photosynthesis.vm; }, "RuBisCO capacity (gC/m2/day)", sum, true);
	add_indiv_output("file_dave_fphen_activity", DAILY, [](const Individual &i) { return i.fphen_activity; }, "Dormancy downregulation (0-1)", mean, false, true);
	add_indiv_output("file_dave_growth_demand", DAILY, [](const Individual &i) { return i.growth_demand; }, "Growth Demand (0-1)", mean, false, true);

	// Daily individual-level array outputs.
	add_indiv_array_output("file_dave_layerwise_fpar", "layer", DAILY, get_nlayer_fpar, [](const Individual &indiv, int i) { return indiv.fpar_layerwise[i]; }, "Daily layered fpar [0-1]");
	add_indiv_array_output("file_dave_layerwise_lai", "layer", DAILY, get_nlayer_lai, [](const Individual &indiv, int i) { return indiv.lai_layerwise[i]; }, "Daily layered LAI (m2/m2)");

	// Subdaily individual-level outputs.
	//                                              column name              function which gets number of array elems     function which returns the i-th output value
	add_indiv_array_output("file_dave_subdaily_an", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.An[i]; }, "Subdaily net photosynthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_rd", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Rd[i]; }, "Subdaily leaf respiration (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_anc", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anc[i]; }, "Subdaily RuBisCO-limited net photosynthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_anj", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.Anj[i]; }, "Subdaily RuBP regeneration-limited photoysnthesis (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_gsw", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.gsw[i]; }, "Subdaily stomatal conductance to water vapour (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_ci", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.ci[i]; }, "Subdaily intercellular CO2 concentration (mol/mol)");
	add_indiv_array_output("file_dave_subdaily_vcmax", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.vcmax[i]; }, "Subdaily max carboxylation rate (mol/m2/s)");
	add_indiv_array_output("file_dave_subdaily_jmax", COLNAME_TIMESTEP, DAILY, [](const Patch &p) { return date.subdaily; }, [](const Individual &indiv, int i) { return indiv.ps_subdaily.jmax[i]; }, "Subdaily max electron transport rate (mol/m2/s)");

	// Annual patch-level outputs.
	add_patch_output("file_dave_patchage", "age", ANNUAL, [](const Patch &p) { return 1.0 * p.age; }, "Annual patch-level age since last disturbance (years)");
	// ...

	// Daily patch-level outputs.
	add_patch_output("file_dave_daylength", "daylength", DAILY, [](const Patch &p) { return p.get_climate().daylength; }, "Daily patch-level day-length (h)");
	add_patch_output("file_dave_soil_nmass_avail", "nmass", DAILY, [](const Patch &p) { return p.soil.get_nmass_avail() * min(1.0, p.fpc_total); }, "Available Soil N for plant uptake");

	// Daily patch-level array outputs.
	add_patch_array_output("file_dave_sw", "sw_", DAILY, NSOILLAYER, [](const Patch &p, int i) { return p.soil.get_layer_soil_water(i); }, "Soil water content as fraction of available water holding capacity", false);
	add_patch_array_output("file_dave_swmm", "swmm_", DAILY, NSOILLAYER, [](const Patch &p, int i) { return p.soil.get_layer_soil_water(i) * p.soil.soiltype.wsats[i]; }, "Soil water content (mm)");
	add_patch_array_output("file_dave_swvol", "swvol_", DAILY, NSOILLAYER, [](const Patch &p, int i) { return p.soil.get_layer_soil_water(i) * p.soil.soiltype.wsats[i] / Dz_soil; }, "Volumetric soil water content (m3/m3)", false);
	// add_daily_patch_level_array_output("file_dave_cfluxes_patch", get_patch_flux_names(), (int)Fluxes::NPERPATCHFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPatchFluxType)i, date.day) * G_PER_KG; }, "Daily patch-level fluxes (g/m2/day)");
	// add_daily_patch_level_array_output("file_dave_cfluxes_pft", get_pft_flux_names(), (int)Fluxes::NPERPFTFLUXTYPES, [](Patch &p, int i) { return p.fluxes.get_daily_flux((Fluxes::PerPFTFluxType)i, date.day) * G_PER_KG; }, "Daily PFT-level fluxes (g/m2/day)");
	add_patch_array_output("file_dave_resp", "resp", DAILY, 1, get_er, "Daily patch-level ecosystem respiration (gC/m2/day)");
	add_patch_array_output("file_dave_gpp", "gpp", DAILY, 1, get_gpp, "Daily patch-level GPP (gC/m2/day)");
	add_patch_array_output("file_dave_npp", "npp", DAILY, 1, get_npp, "Daily patch-level NPP (gC/m2/day)");
	add_patch_array_output("file_dave_nee", "nee", DAILY, 1, get_nee, "Daily patch-level NEE (gC/m2/day)");
	add_patch_array_output("file_dave_et", "et", DAILY, 1, get_et, "Daily patch-level evapotranspiration (mm/day)");

	// Subdaily patch-level outputs.
	//                                           column name              function which gets number of array elems     function which returns the i-th output value
	add_patch_array_output("file_dave_met_temp", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.temps[i]; }, "Subdaily temperature (°C)");
	add_patch_array_output("file_dave_met_par", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pars[i]; }, "Subdaily PAR (J/m2/timestep)");
	add_patch_array_output("file_dave_met_vpd", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.vpds[i]; }, "Subdaily VPD (kPa)");
	add_patch_array_output("file_dave_met_insol", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.insols[i]; }, "Subdaily insolation (units depend on instype)");
	add_patch_array_output("file_dave_met_precip", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.precs[i]; }, "Subdaily precipitation (mm)");
	add_patch_array_output("file_dave_met_pressure", COLNAME_TIMESTEP, DAILY, []() { return date.subdaily; }, [](const Patch &patch, int i) { return patch.stand.get_gridcell().climate.pressure[i]; }, "Subdaily air pressure (Pa)");

	// Standpft-level outputs.
	add_spft_output("file_dave_cmass_repr", DAILY, [](const Standpft &spft) { return spft.cmass_repr; }, "Daily reproductive cmass at standpft-level (gC/m2/day)");

	// Register instruction file parameters.
	declare_parameters();

	total_time = get_wall_time() - t_start;
}

DaveOutput::~DaveOutput() {
	printf("Total time spent in dave output: %.2fs\n", total_time);
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

bool should_write_today() {
	return spinup_output || date.year >= nyear_spinup;
}

void DaveOutput::outannual(Gridcell &gridcell) {
	double t0 = get_wall_time();
	if (should_write_today()) {
		write_outputs(gridcell, annual_outputs);
	}
	total_time += get_wall_time() - t0;
}

void DaveOutput::outdaily(Gridcell &gridcell) {
	double t0 = get_wall_time();
	if (should_write_today()) {
		write_outputs(gridcell, daily_outputs);
	}
	total_time += get_wall_time() - t0;
}

void DaveOutput::add_patch_output(std::string name, std::string colname, timestep_t timestep, double (*read_output)(const Patch&), std::string description, bool sum, bool mean) {
	std::vector<std::string> colnames;
	colnames.push_back(colname);
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colnames, 1, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_array_output(std::string name, std::vector<std::string> colnames, timestep_t timestep, double (*read_output)(const Patch&, int), std::string description, bool sum, bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colnames, colnames.size(), reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_array_output(std::string name, std::string colname, timestep_t timestep, int nelem, double (*read_output)(const Patch&, int), std::string description, bool sum, bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, nelem, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_patch_array_output(std::string name, std::string colname, timestep_t timestep, int (*get_nelem)(), double (*read_output)(const Patch&, int), std::string description, bool sum, bool mean) {
	ArrayOutputReader<Patch> *reader = new ArrayOutputReader<Patch>(read_output);
	PatchLevelArrayOutput *output = new PatchLevelArrayOutput(name, colname, get_nelem, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_output(std::string name, timestep_t timestep, double (*read_output)(const Individual&), std::string description, double (*aggregator)(double*, int), bool sum, bool mean) {
	OutputReader<Individual> *reader = new OutputReader<Individual>(read_output);
	auto output = new IndividualOutput(name, reader, description, aggregator, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_pft_output(std::string name, timestep_t timestep, double (*read_output)(const Patch &, const Pft&), std::string description, bool sum, bool mean) {
	PftOutput *output = new PftOutput(name, read_output, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::add_indiv_array_output(std::string name, std::string elem_name, timestep_t timestep, int (*get_nelem)(const Patch &), double (*read_output)(const Individual&, int), std::string description) {
	ArrayOutputReader<Individual> *reader = new ArrayOutputReader<Individual>(read_output);
	IndivLevelArrayOutput *output = new IndivLevelArrayOutput(name, elem_name, get_nelem, reader, description);
	register_output(output, timestep);
}

void DaveOutput::add_spft_output(
		std::string name,
		timestep_t timestep,
		double (*read_output)(const Standpft&),
		std::string description,
		bool sum,
		bool mean) {
	OutputReader<Standpft> *reader = new OutputReader<Standpft>(read_output);
	StandPftOutput *output = new StandPftOutput(name, reader, description, sum, mean);
	register_output(output, timestep);
}

void DaveOutput::register_output(OutputFile *file, timestep_t timestep) {
	auto &vec = get_outputs_vector(timestep);
	// vec.push_back(file);
	vec.emplace_back(file);
	// vec.push_back(std::move(file));
}

std::vector<std::unique_ptr<DaveOutput::OutputFile>> &DaveOutput::get_outputs_vector(timestep_t timestep) {
	if (timestep == ANNUAL) {
		return annual_outputs;
	}
	return daily_outputs;
}

void DaveOutput::declare_parameters() {
	declare_parameter("spinup_output", &spinup_output, "True to write DAVE outputs during the spinup, false otherwise.");
	declare_parameter("patch_output_stats", &patch_output_stats, "True to ONLY write stats (mean/std. deviation) for patch-level outputs. False to write the full patch-level outputs.");

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
	Table tbl = channel.create_table(TableDescriptor(name, get_cols()));
	table = new Table(tbl.id());
}

void DaveOutput::OutputFile::close(OutputChannel &channel) {
	assert_table_open();
	channel.close_table(*table);
}

bool DaveOutput::OutputFile::should_write(OutputRows &rows) {
	if (!should_write_today()) {
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

DaveOutput::PatchLevelOutput::PatchLevelOutput(std::string param, std::string desc, bool sum, bool mean)
: OutputFile(param, desc) {
	write_sum = sum;
	write_mean = mean;
}

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

		std::vector<std::vector<double>> data;

		while (stand.isobj) {
			Patch &patch = stand.getobj();

			std::vector<double> patch_data = write_patch(patch, rows);
			if (patch_output_stats) {
				data.push_back(patch_data);
			} else {
				// The code would be cleaner if we were to use the data 2d array
				// in both cases, but I suspect that performance will be better
				// with this shortcut when the number of patches is high.

				// Write patch ID.
				write_output(rows, patch.id);

				// Write patch data.
				double total = 0;
				for (double x : patch_data) {
					write_output(rows, x);
					total += x;
				}
				if (write_sum) {
					write_output(rows, total);
				}
				if (write_mean) {
					write_output(rows, total / patch_data.size());
				}
			}

			if (!patch_output_stats) {
				finish_row(channel, lon, lat, year, day);
			}

			stand.nextobj();
		}

		if (patch_output_stats && data.size()) {
			// Each element of data contains all columns for each patch. We want
			// to write one column at a time.

			// First, get the number of columns.
			int ncol = data[0].size();

			// Now iterate over the columns that we want to write.
			for (size_t j = 0; j < ncol; j++) {
				// Get this column's value for each patch.
				std::vector<double> column_data;
				double total = 0;
				for (size_t i = 0; i < data.size(); i++) {
					column_data.push_back(data[i][j]);
				}
				// Compute summary statistics for this column.
				double mean = 0, sd = 0, median = 0;
				compute_stats(column_data, &mean, &sd, &median);
				write_output(rows, mean);
				write_output(rows, median);
				write_output(rows, sd);

			}
			if (write_sum || write_mean) {
				std::vector<double> totals;
				// Iterate over patches.
				for (int i = 0; i < data.size(); i++) {
					// Get the total for this patch.
					double total = 0;
					for (size_t j = 0; j < data[i].size(); j++) {
						total += data[i][j];
					}
					totals.push_back(total);
				}
				double mean = 0, sd = 0, median = 0;
				compute_stats(totals, &mean, &sd, &median);
				if (write_sum) {
					write_output(rows, mean);
					write_output(rows, median);
					write_output(rows, sd);
				}
				if (write_mean) {
					write_output(rows, mean / ncol);
					write_output(rows, median / ncol);
					write_output(rows, sd / ncol);
				}
			}
			finish_row(channel, lon, lat, year, day);
		}

		gc_itr++;
	}
}

void add_column_with_suffix(ColumnDescriptors &cols, const char *colname, const char *sfx, int width, int precision) {
	std::stringstream ss;
	ss << colname << "_" << sfx;
	std::string sname = ss.str();
	const char *name = sname.c_str();
	width = max((size_t)width, sname.length() + 1);
	cols += ColumnDescriptor(name, width, precision);
}

void add_column(ColumnDescriptors &cols, const char *colname, int width, int precision) {
	if (patch_output_stats) {
		add_column_with_suffix(cols, colname, "mean", width, precision);
		add_column_with_suffix(cols, colname, "median", width, precision);
		add_column_with_suffix(cols, colname, "sd", width, precision);
	} else {
		cols += ColumnDescriptor(colname, width, precision);
	}
}

void add_columns(ColumnDescriptors &cols, std::vector<std::string> names, int width, int precision) {
	for (std::string name : names) {
		add_column(cols, name.c_str(), width, precision);
	}
}

ColumnDescriptors create_column_descriptors() {
	ColumnDescriptors columns;
	if (!patch_output_stats) {
		add_column(columns, COLNAME_PATCH, COLUMN_WIDTH, 0);
	}
	return columns;
}

ColumnDescriptors DaveOutput::PatchLevelOutput::get_cols() {
	ColumnDescriptors columns = get_columns();
	if (write_sum) {
		add_column(columns, COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	if (write_mean) {
		add_column(columns, COLNAME_MEAN, COLUMN_WIDTH, NUMERIC_PRECISION);
	}
	return columns;
}

// #############################################################################
// IndividualOutput implementation.
// #############################################################################

DaveOutput::IndividualOutput::IndividualOutput(
		std::string param_name,
		OutputReader<Individual> *get_output,
		std::string description,
		double (*agg)(double*, int),
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	read_output = std::unique_ptr<OutputReader<Individual>>(get_output);
	aggregator = agg;
}

DaveOutput::IndividualOutput::~IndividualOutput() { }

std::vector<double> DaveOutput::IndividualOutput::write_patch(Patch &patch, OutputRows &out) {
	Vegetation &vegetation = patch.vegetation;

	std::vector<std::string> pfts = get_pfts(false);
	std::vector<double> results;

	// Write one column per PFT - note that there will not necessarily be a 1:1
	// mapping of PFTs to individuals, as some PFTs may not be established.
	for (int i = 0; i < pfts.size(); i++) {
		std::vector<Individual *> indivs = get_indivs(vegetation, pfts[i].c_str());

		double result = 0;
		if (indivs.size() > 0) {
			for (Individual *indiv : indivs) {
				result += indiv ? read_output->read(*indiv) : 0;
			}
			// result /= indivs.size();
		}
		results.push_back(result);
	}

	return results;
}

ColumnDescriptors DaveOutput::IndividualOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	add_columns(columns, get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// PftOutput implementation.
// #############################################################################

DaveOutput::PftOutput::PftOutput(
	std::string param_name,
	double (*get_output)(const Patch &, const Pft &),
	std::string description,
	bool sum,
	bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	read_output = get_output;
}

ColumnDescriptors DaveOutput::PftOutput::get_columns() {
	ColumnDescriptors columns;

	std::vector<std::string> pfts = get_pfts(true);

	if (!patch_output_stats) {
		add_column(columns, COLNAME_PATCH, COLUMN_WIDTH, 0);
	}
	add_columns(columns, pfts, COLUMN_WIDTH, NUMERIC_PRECISION);

	return columns;
}

std::vector<double> DaveOutput::PftOutput::write_patch(Patch &patch, OutputRows &out) {
	std::vector<double> results;

	for (int i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];

		double result = read_output(patch, pft);
		results.push_back(result);
	}

	return results;
}

DaveOutput::PftOutput::~PftOutput() { }

// #############################################################################
// CustomPatchLevelOutput implementation.
// #############################################################################

DaveOutput::CustomPatchLevelOutput::CustomPatchLevelOutput(
	std::string param_name,
	std::vector<ColumnDefinition> (*get_columns)(),
	std::string description,
	bool sum,
	bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
		get_column_definitions = get_columns;
}

std::vector<double> DaveOutput::CustomPatchLevelOutput::write_patch(Patch &patch, OutputRows &out) {
	std::vector<double> result;

	for (const ColumnDefinition &col : columns) {
		result.push_back(col.read_value(patch));
	}

	return result;
}

ColumnDescriptors DaveOutput::CustomPatchLevelOutput::get_columns() {
	this->columns = this->get_column_definitions();

	ColumnDescriptors columns;

	if (!patch_output_stats) {
		columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	}

	for (const ColumnDefinition &column : this->columns) {
		std::string name = column.get_name();
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(columns, name.c_str(), width, NUMERIC_PRECISION);
	}

	return columns;
}

DaveOutput::CustomPatchLevelOutput::~CustomPatchLevelOutput() { }

// #############################################################################
// StandPftOutput implementation.
// #############################################################################

DaveOutput::StandPftOutput::StandPftOutput(
		std::string param_name,
		OutputReader<Standpft> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	this->read_output = std::unique_ptr<OutputReader<Standpft>>(get_output);
}

DaveOutput::StandPftOutput::~StandPftOutput() { }

std::vector<double> DaveOutput::StandPftOutput::write_patch(Patch &patch, OutputRows &out) {
	Vegetation &vegetation = patch.vegetation;
	Stand &stand = patch.stand;

	std::vector<double> results;

	// Note that patch id is written before this function is called, so we
	// shouldn't write it ourselves (that would result in a 'too many columns'
	// error).

	// Write one column per PFT.
	for (int i = 0; i < pftlist.nobj; i++) {
		Pft &pft = pftlist[i];
		const Standpft &spft = stand.pft[pft.id];
		double result = read_output->read(spft);
		results.push_back(result);
	}

	return results;
}

ColumnDescriptors DaveOutput::StandPftOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();
	add_columns(columns, get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// PatchLevelArrayOutput implementation.
// #############################################################################

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int n,
		ArrayOutputReader<Patch> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {

	nelem = n;
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
	get_nelem = NULL;

	for (int i = 0; i < n; i++) {
		if (n == 1) {
			column_names.push_back(col_name);
		} else {
			column_names.push_back(col_name + std::to_string(i));
		}
	}
}

DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
	std::string param_name,
	std::vector<std::string> col_names,
	int n,
	ArrayOutputReader<Patch> *get_output,
	std::string description,
	bool sum,
	bool mean)
 : PatchLevelOutput(param_name, description, sum, mean) {
	column_names = col_names;
	nelem = n;
	get_nelem = NULL;
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
}
DaveOutput::PatchLevelArrayOutput::PatchLevelArrayOutput(
		std::string param_name,
		std::string col_name,
		int (*get_nelem)(),
		ArrayOutputReader<Patch> *get_output,
		std::string description,
		bool sum,
		bool mean) : PatchLevelOutput(param_name, description, sum, mean) {
	get_nelem = get_nelem;
	column_names.push_back(col_name);
	read_output = std::unique_ptr<ArrayOutputReader<Patch>>(get_output);
}

DaveOutput::PatchLevelArrayOutput::~PatchLevelArrayOutput() { }

std::vector<double> DaveOutput::PatchLevelArrayOutput::write_patch(
		Patch &patch, OutputRows &out) {
	std::vector<double> results;
	for (int i = 0; i < nelem; i++) {
		double value = read_output->read(patch, i);
		results.push_back(value);
	}
	return results;
}

ColumnDescriptors DaveOutput::PatchLevelArrayOutput::get_columns() {
	ColumnDescriptors columns = create_column_descriptors();

	if (get_nelem) {
		nelem = get_nelem();
		std::string basename = column_names[0];
		column_names.clear();
		for (int i = 0; i < nelem; i++) {
			if (nelem == 1) {
				column_names.push_back(basename);
			} else {
				column_names.push_back(basename + std::to_string(i));
			}
		}
	}

	for (int i = 0; i < nelem; i++) {
		std::string name = column_names[i];
		int width = max((size_t)COLUMN_WIDTH, name.length() + 1);
		add_column(columns, name.c_str(), width, NUMERIC_PRECISION);
	}
	return columns;
}

// #############################################################################
// IndivLevelArrayOutput implementation.
// #############################################################################

DaveOutput::IndivLevelArrayOutput::IndivLevelArrayOutput(
		std::string param_name,
		std::string element_name,
		int (*get_nelem)(const Patch &),
		ArrayOutputReader<Individual> *get_output,
		std::string description) : OutputFile(param_name, description)
{
	elem_name = element_name;
	get_n = get_nelem;
	read_elem = std::unique_ptr<ArrayOutputReader<Individual>>(get_output);
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

		for (int j = 0; j < stand.nobj; j++) {
			Patch &patch = stand[j];

			int nelem = get_n(patch);
			Vegetation &vegetation = patch.vegetation;
			for (int i = 0; i < nelem; i++) {

				// Write patch ID.
				write_output(rows, patch.id);

				// Write array index.
				write_output(rows, i);

				vegetation.firstobj();
				double total = 0;

				for (int k = 0; k < pfts.size(); k++) {
					std::vector<Individual *> indivs = get_indivs(vegetation, pfts[k].c_str());
					double indiv_total = 0;
					for (Individual *indiv : indivs) {
						// Read array element for this individual.
						indiv_total += read_elem->read(*indiv, i);
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
		gc_itr++;
	}
}

ColumnDescriptors DaveOutput::IndivLevelArrayOutput::get_cols() {
	ColumnDescriptors columns;
	columns += ColumnDescriptor(COLNAME_PATCH, COLUMN_WIDTH, 0);
	columns += ColumnDescriptor(elem_name.c_str(), COLUMN_WIDTH, 0);
	columns += ColumnDescriptors(get_pfts(true), COLUMN_WIDTH, NUMERIC_PRECISION);
	columns += ColumnDescriptor(COLNAME_TOTAL, COLUMN_WIDTH, NUMERIC_PRECISION);
	return columns;
}

// #############################################################################
// ColumnDefinition implementation.
// #############################################################################

ColumnDefinition::ColumnDefinition(std::string name, double (*get_output)(const Patch &)) {
	this->name = name;
	read_output = get_output;
}

std::string ColumnDefinition::get_name() const {
	return name;
}

double ColumnDefinition::read_value(const Patch &patch) const {
	return read_output(patch);
}

bool spinup_output;
bool patch_output_stats;

} // namespace GuessOutput
