#include "daveoutputhelpers.h"
#include "guess.h"
#include "commonoutput.h"
#include "simfire.h"
#include "blaze.h"
#include "ncutil.h"

namespace GuessOutput {

void ensure_channel_open() {
	if (!output_channel) {
		fail("Output channel is NULL");
	}
}

double get_lai(const Individual &indiv) {
	return indiv.lai;
}

double get_er(const Patch &patch) {
	double ra_kg = patch.fluxes.get_daily_flux(Fluxes::RA, date.day);
	double rh_kg = patch.fluxes.get_daily_flux(Fluxes::SOILC, date.day);
	return (ra_kg + rh_kg) * G_PER_KG;
}

double get_gpp(const Patch &patch) {
	double gpp_kg = patch.fluxes.get_daily_flux(Fluxes::GPP, date.day);
	return gpp_kg * G_PER_KG;
}

double get_npp(const Patch &patch) {
	double npp_kg = patch.fluxes.get_daily_flux(Fluxes::NPP, date.day);
	return npp_kg * G_PER_KG;
}

double get_nee(const Patch &patch) {
	double rh_kg = patch.fluxes.get_daily_flux(Fluxes::SOILC, date.day);
	double npp_kg = patch.fluxes.get_daily_flux(Fluxes::NPP, date.day);
	return (rh_kg - npp_kg) * G_PER_KG;
}

double get_anee(const Patch &patch) {
	double rh_kg = patch.fluxes.get_annual_flux(Fluxes::SOILC);
	double npp_kg = patch.fluxes.get_annual_flux(Fluxes::NPP);
	return (rh_kg - npp_kg);
}

double get_et(const Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;
	double transpiration = 0;
	for (int i = 0; i < vegetation.nobj; i++)
		transpiration += vegetation[i].aet;
	return transpiration + patch.soil.get_layer_soil_water_evap();
}

double get_patch_ba(const Patch &patch) {
	double ba = 0;
	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		const Individual &indiv = patch.vegetation[i];
		if (indiv.pft.lifeform == TREE) {
			ba += indiv.basalarea;
		}
	}
	return ba * 10000.0 / patcharea;
}

std::vector<Individual *> get_indivs(const Vegetation &vegetation, const char *name) {
	std::vector<Individual *> result;
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		if (indiv.pft.name == name) {
			result.push_back(&indiv);
		}
	}
	return result;
}

void remove_spaces(char *x) {
	char *y = x;
	char nxt;
	do
	{
		bool space = false;
		while (isspace(*y)) {
			y++;
			space = true;
		}
		nxt = *y++;
		if (space) {
			nxt = toupper(nxt);
		}
	} while (*x++ = nxt);
}

std::vector<std::string> get_pfts(bool remove_whitespace) {
	std::vector<std::string> pfts;
	pftlist.firstobj();
	while (pftlist.isobj) {
		char *name = pftlist.getobj().name;
		if (remove_whitespace) {
			remove_spaces(name);
		}
		pfts.push_back(name);
		pftlist.nextobj();
	}
	return pfts;
}

std::vector<std::string> get_patch_flux_names() {
	std::vector<std::string> fluxes;
	fluxes.push_back("firec");
	fluxes.push_back("soilc");
	fluxes.push_back("estc");
	fluxes.push_back("harvestc");
	fluxes.push_back("seedc");
	fluxes.push_back("manurec");
	fluxes.push_back("manuren");
	fluxes.push_back("nfert");
	fluxes.push_back("estn");
	fluxes.push_back("harvestn");
	fluxes.push_back("seedn");
	fluxes.push_back("nh3_fire");
	fluxes.push_back("nox_fire");
	fluxes.push_back("n2o_fire");
	fluxes.push_back("n2_fire");
	fluxes.push_back("nh3_soil");
	fluxes.push_back("no_soil");
	fluxes.push_back("n2o_soil");
	fluxes.push_back("n2_soil");
	fluxes.push_back("doc_flux");
	fluxes.push_back("net_nitrif");
	fluxes.push_back("net_denitrif");
	fluxes.push_back("gross_nitrif");
	fluxes.push_back("gross_denitrif");
	fluxes.push_back("reprc");
	fluxes.push_back("ch4c");
	fluxes.push_back("ch4c_diff");
	fluxes.push_back("ch4c_plan");
	fluxes.push_back("ch4c_ebul");
	return fluxes;
}

std::vector<std::string> get_pft_flux_names() {
	std::vector<std::string> fluxes;
		/// NPP (kgC/m2)
		fluxes.push_back("npp");
		/// GPP (kgC/m2)
		fluxes.push_back("gpp");
		/// Autotrophic respiration (kgC/m2)
		fluxes.push_back("ra");
		/// Isoprene (mgC/m2)
		fluxes.push_back("iso");
		/// Monoterpene (mgC/m2)
		fluxes.push_back("mt_apin");
		fluxes.push_back("mt_bpin");
		fluxes.push_back("mt_limo");
		fluxes.push_back("mt_myrc");
		fluxes.push_back("mt_sabi");
		fluxes.push_back("mt_camp");
		fluxes.push_back("mt_tric");
		fluxes.push_back("mt_tboc");
		fluxes.push_back("mt_othr");
	return fluxes;
}

size_t get_nlayer_canopy(const Patch &patch) {
	size_t nlayer = 0;
	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		const Individual &indiv = patch.vegetation[i];
		nlayer = max(nlayer, indiv.canopy.size());
	}
	return nlayer;
}

std::vector<ColumnDefinition<Individual>> get_pft_litter_pools() {
	std::vector<ColumnDefinition<Individual>> columns;
	columns.push_back(ColumnDefinition<Individual>("Leaf", [](const Individual &i) { return i.patchpft().cmass_litter_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("Root", [](const Individual &i) { return i.patchpft().cmass_litter_root; }));
	columns.push_back(ColumnDefinition<Individual>("Sap", [](const Individual &i) { return i.patchpft().cmass_litter_sap; }));
	columns.push_back(ColumnDefinition<Individual>("Heart", [](const Individual &i) { return i.patchpft().cmass_litter_heart; }));
	columns.push_back(ColumnDefinition<Individual>("Repr", [](const Individual &i) { return i.patchpft().cmass_litter_repr; }));
	return columns;
}

double get_pft_standing_dead_cmass(const Patch &patch, const Pft &pft) {
	const Patchpft &ppft = patch.pft[pft.id];
	return ppft.cmass_litter_leaf + ppft.cmass_litter_root + ppft.cmass_litter_sap + ppft.cmass_litter_heart + ppft.cmass_litter_repr;
}

double get_pft_standing_dead_nmass(const Patch &patch, const Pft &pft) {
	const Patchpft &ppft = patch.pft[pft.id];
	return ppft.nmass_litter_leaf + ppft.nmass_litter_root + ppft.nmass_litter_sap + ppft.nmass_litter_heart + ppft.nmass_litter_repr;
}

double get_soil_moisture(const Soil &soil, int start_layer, int end_layer) {
	assert(start_layer < end_layer);
	double swmm = 0;
	for (size_t i = start_layer; i < end_layer; i++) {
		swmm += soil.get_layer_soil_water(i) * soil.soiltype.awc[i] + soil.soiltype.wp[i];
	}
	return swmm;
}

double get_available_soil_moisture(const Soil &soil, int start_layer, int end_layer) {
	assert(start_layer < end_layer);
	double swmm = 0;
	for (size_t i = start_layer; i < end_layer; i++) {
		swmm += soil.get_layer_soil_water(i) * soil.soiltype.awc[i];
	}
	return swmm;
}

double get_6m_soil_moisture_estimate(const Soil &soil) {
	// Soil goes down to 1.5m. Get fractional content of the bottom 0.5m (ie
	// 1m-1.5m), then multiply this 0.5m by 12 to get an esimate of water
	// content to 6m.
	// TODO: need to double check this. Does this even make physical sense?
	return get_soil_moisture(soil, 10, NSOILLAYER) * 12;
}

void compute_stats(std::vector<double> data, double *mean, double *sd, double *median) {
	double sum = 0;
	for (double x : data) {
		sum += x;
	}
	*mean = sum / data.size();;
	double var = 0;
	for (double x : data) {
		var += (x - *mean) * (x - *mean);
	}
	*sd = data.size() > 1 ? sqrt(var / (data.size() - 1)) : 0;
	std::sort(data.begin(), data.end());
	int n = data.size() / 2;
	if (data.size() % 2 == 0) {
		*median = (data[n] + data[n - 1]) / 2;
	} else {
		*median = data[n];
	}
}

double sum(const double *array, size_t n) {
	double total = 0;
	for (int i = 0; i < n; i++) {
		total += array[i];
	}
	return total;
}

double mean(const double *array, size_t n) {
	return n ? sum(array, n) / n : 0;
}

double weight_by_lai(const Individual &indiv, double (*selector)(const Individual&)) {
	double total_lai = 0;
	for (size_t i = 0; i < indiv.vegetation.nobj; i++) {
		if (indiv.pft.id == indiv.vegetation[i].id) {
			total_lai += indiv.vegetation[i].lai;
		}
	}
	return selector(indiv) * indiv.lai / total_lai;
}

double get_patch_nee(const Patch &patch) {
	// out_cflux, 
	double flux_veg = -patch.fluxes.get_annual_flux(Fluxes::NPP);
	double flux_repr = -patch.fluxes.get_annual_flux(Fluxes::REPRC);
	double flux_soil = patch.fluxes.get_annual_flux(Fluxes::SOILC);
	double flux_fire = patch.fluxes.get_annual_flux(Fluxes::FIREC);
	double flux_est = patch.fluxes.get_annual_flux(Fluxes::ESTC);
	double flux_seed = patch.fluxes.get_annual_flux(Fluxes::SEEDC);
	double flux_charvest = patch.fluxes.get_annual_flux(Fluxes::HARVESTC);
	double c_org_leach_gridcell = patch.soil.aorgCleach;
	Landcover &lc = patch.get_climate().gridcell.landcover;
	return flux_veg - flux_repr + flux_soil + flux_fire + flux_est + c_org_leach_gridcell +
	flux_seed + flux_charvest + lc.acflux_landuse_change + lc.acflux_harvest_slow;
}

std::vector<ColumnDefinition<Patch>> runoff_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("Surf", [](const Patch &p) { return p.asurfrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("Drain", [](const Patch &p) { return p.adrainrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("Base", [](const Patch &p) { return p.abaserunoff; }));
	columns.push_back(ColumnDefinition<Patch>("Precip", [](const Patch &p) { return p.get_climate().aprec; }));
	columns.push_back(ColumnDefinition<Patch>("NetInf", [](const Patch &p) { return p.stand.get_gridcell().climate.aprec - (p.asurfrunoff + p.adrainrunoff + p.abaserunoff); }));

	return columns;
}

double get_cmass_patch(const Patch &patch) {
	double ccont = 0;
	for (unsigned int i = 0; i < patch.vegetation.nobj; i++) {
		ccont += patch.vegetation[i].ccont();
	}
	return ccont;
}

double get_standing_c(const Patch &patch) {
	double standing = 0;

	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];
		standing += get_pft_standing_dead_cmass(patch, pft);
	}

	return standing;
}

double get_standing_n(const Patch &patch) {
	double standing = 0;
	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];
		standing += get_pft_standing_dead_nmass(patch, pft);
	}
	return standing;
}

double get_cmass_century(const Patch &patch) {
	double centuryc = 0;
	for (int i = 0; i < NSOMPOOL; i++) {
		if (i != SURFMETA && i != SURFSTRUCT && i != SURFFWD && i != SURFCWD)
			centuryc += patch.soil.sompool[i].cmass;
	}
	return centuryc;
}

double get_nmass_century(const Patch &patch) {
	double centuryn = 0;
	for (int i = 0; i < NSOMPOOL; i++) {
		if (i != SURFMETA && i != SURFSTRUCT && i != SOILMETA && i != SOILSTRUCT && i != SURFFWD && i != SURFCWD)
		centuryn += patch.soil.sompool[i].nmass;
	}
	return centuryn;
}

double get_cmass_harv_slowc(const Patch &patch) {
	double c_harv_slow = 0;
	for (int i = 0; i < npft; i++) {
		c_harv_slow += patch.pft[i].cmass_harvested_products_slow;
	}
	return c_harv_slow;
}

double get_cflux_veg(const Patch &patch) {
	return -patch.fluxes.get_annual_flux(Fluxes::NPP);
}

double get_litter_n(const Patch &patch) {
	return patch.soil.sompool[SURFMETA].nmass +
		   patch.soil.sompool[SURFSTRUCT].nmass +
           patch.soil.sompool[SURFFWD].nmass +
           patch.soil.sompool[SURFCWD].nmass;
}

double get_litter_c(const Patch &patch) {
	return patch.soil.sompool[SURFMETA].cmass +
	       patch.soil.sompool[SURFSTRUCT].cmass +
	       patch.soil.sompool[SURFFWD].cmass +
	       patch.soil.sompool[SURFCWD].cmass;
}

double get_nmass_litter(const Patch &patch) {
	double nmass_litter = 0;
	for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		const Patchpft &ppft = patch.pft[pftlist.getobj().id];
		nmass_litter += ppft.nmass_litter_leaf + ppft.nmass_litter_root + ppft.nmass_litter_sap + ppft.nmass_litter_heart + ppft.nmass_litter_repr;
	}
	return nmass_litter;
}

double get_nmass_vegetation(const Patch &patch) {
	double nmass = 0;
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		nmass += patch.vegetation[i].ncont();
	}
	return nmass;
}

double get_nmass_patch(const Patch &patch) {
	return get_nmass_vegetation(patch) + get_nmass_litter(patch);
}

double get_avail_n(const Patch &patch) {
	return patch.soil.NH4_mass + patch.soil.NO3_mass + patch.soil.snowpack_NH4_mass + patch.soil.snowpack_NO3_mass;
}

double get_centuryn(const Patch &patch) {
	return get_nmass_century(patch) + get_avail_n(patch);
}

double get_nmass_harv_slowc(const Patch &patch) {
	double n_harv_slow = 0;

	for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		Patchpft &ppft = patch.pft[pftlist.getobj().id];
		n_harv_slow += ppft.nmass_harvested_products_slow;
	}

	return n_harv_slow;
}

/*
Aggregate a value across all ppfts in the patch.
*/
double aggregate_ppft(const Patch &patch, double (*read_ppft)(const Patchpft &), double (*aggregate)(const double*, int)) {
	std::vector<double> elems(patch.pft.nobj);
	for (int i = 0; i < patch.pft.nobj; i++) {
		const Patchpft &ppft = patch.pft[i];
		elems[i] = read_ppft(ppft);
	}
	return elems.size() ? aggregate(&elems.front(), patch.pft.nobj) : 0;
}

/*
Aggregate a value across all individuals in the patch.
*/
double aggregate_indiv(const Patch &patch, double (*read_indiv)(const Individual &), double (*aggregate)(const double*, int)) {
	const Vegetation &vegetation = patch.vegetation;
	std::vector<double> elems(vegetation.nobj);
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		elems[i] = read_indiv(indiv);
	}
	return aggregate(&elems.front(), vegetation.nobj);
}

/*
Return the sum of all elements in the array.
*/
double sum_array(const double *array, int nelem) {
	double sum = 0;
	for (int i = 0; i < nelem; i++) {
		sum += array[i];
	}
	return sum;
}

/*
Get the sum of a value for each ppft.
*/
double sum_ppft(const Patch &patch, double (*read_ppft)(const Patchpft &)) {
	return aggregate_ppft(patch, read_ppft, sum_array);
}

/*
Get the sum of a value for each individual.
*/
double sum_indiv(const Patch &patch, double (*read_indiv)(const Individual &)) {
	return aggregate_indiv(patch, read_indiv, sum_array);
}

double get_fine_litter(const Patch &patch) {
	if (patch.avail_fuel == 0) {
		return 0;
	}
	double litter = sum_ppft(patch, [](const Patchpft &ppft) { return ppft.cmass_litter_leaf; });
	const Sompool *sompool = patch.soil.sompool;
	double som = sompool[SURFSTRUCT].cmass + sompool[SURFMETA].cmass;
	double fine_litter = som + litter;
	return TURNOVERFRACT[11][0] * fine_litter * G_PER_KG;
}

double get_fine_woody_debris(const Patch &patch) {
	if (patch.avail_fuel == 0) {
		return 0;
	}
	double trans_litter_sap = sum_ppft(patch, [](const Patchpft &ppft) { return ppft.cmass_litter_sap; });
	const Sompool *sompool = patch.soil.sompool;
	return TURNOVERFRACT[10][0] * (sompool[SURFFWD].cmass + trans_litter_sap) * G_PER_KG;
}

const Gridcell &get_gridcell(const Patch &p) {
	return p.get_climate().gridcell;
}

double get_coarse_woody_debris(const Patch &patch) {
	if (patch.avail_fuel == 0) {
		return 0;
	}
	double trans_litter_heart = sum_ppft(patch, [](const Patchpft &ppft) { return ppft.cmass_litter_heart; });
	const Sompool *sompool = patch.soil.sompool;
	double k = get_gridcell(patch).k_tun_litter;
	return TURNOVERFRACT[ 9][0] * k * (sompool[SURFCWD].cmass + trans_litter_heart) * k * G_PER_KG;
}

double get_half_cmass_grass(const Patch &patch) {
	if (patch.avail_fuel == 0) {
		return 0;
	}
	return sum_indiv(patch, [](const Individual &i) {
		return i.pft.is_grass() ? 0.5 * i.cmass_leaf_today() * G_PER_KG : 0;
	});
}

double get_heartwood_litter(const Patch &patch) {
	double k = get_gridcell(patch).k_tun_litter;
	double litter = sum_ppft(patch, [](const Patchpft &ppft) { return ppft.cmass_litter_heart; });
	// litter *= k * TURNOVERFRACT[ 9][0] * k
	return litter * G_PER_KG;
}

double get_coarse_woody_som(const Patch &patch) {
	double k = get_gridcell(patch).k_tun_litter;
	double som = patch.soil.sompool[SURFCWD].cmass;
	// som *= k * TURNOVERFRACT[ 9][0] * k;
	return som * G_PER_KG;
}

std::vector<ColumnDefinition<Patch>> cpool_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("VegC", get_cmass_patch));

	if (!ifcentury) {
		columns.push_back(ColumnDefinition<Patch>("StandingC", get_standing_c));
		columns.push_back(ColumnDefinition<Patch>("LittC", get_litter_c));
		columns.push_back(ColumnDefinition<Patch>("SoilfC", [](const Patch &p) { return p.soil.cpool_fast; }));
		columns.push_back(ColumnDefinition<Patch>("SoilsC", [](const Patch &p) { return p.soil.cpool_slow; }));
	} else {
		columns.push_back(ColumnDefinition<Patch>("StandingC", get_standing_c));
		columns.push_back(ColumnDefinition<Patch>("LitterC", get_litter_c));
		columns.push_back(ColumnDefinition<Patch>("SoilC", get_cmass_century));
	}
	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition<Patch>("HarvSlowC", get_cmass_harv_slowc));
	}

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_dbhcoloumns_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("diam", [](const Individual &i) { return i.diam; }));
	columns.push_back(ColumnDefinition<Individual>("deltadiam", [](const Individual &i) { return i.d_diam; }));

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_cpool_tree_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("cmass_leaf", [](const Individual &i) { return i.cmass_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_root", [](const Individual &i) { return i.cmass_root; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_sap", [](const Individual &i) { return i.cmass_sap; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_heart", [](const Individual &i) { return i.cmass_heart; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_repr", [](const Individual &i) { return i.cmass_repr; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_storage", [](const Individual &i) { return i.cmass_storage; }));

	columns.push_back(ColumnDefinition<Individual>("cmass_leaf_limit", [](const Individual &i) { return i.cmass_leaf_limit; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_root_limit", [](const Individual &i) { return i.cmass_root_limit; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_sap_limit", [](const Individual &i) { return i.cmass_sap_limit; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_storage_limit", [](const Individual &i) { return i.cmass_storage_limit; }));


	return columns;
}

std::vector<ColumnDefinition<Patch>> cflux_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("Veg", [](const Patch &p) { return -p.fluxes.get_annual_flux(Fluxes::NPP); }));
	columns.push_back(ColumnDefinition<Patch>("Repr", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::REPRC); }));
	columns.push_back(ColumnDefinition<Patch>("Soil", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::SOILC) + p.soil.aorgCleach; }));
	columns.push_back(ColumnDefinition<Patch>("Fire", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::FIREC); }));
	columns.push_back(ColumnDefinition<Patch>("Est", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::ESTC); }));
	if (run_landcover) {
		columns.push_back(ColumnDefinition<Patch>("Seed", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::SEEDC); }));
		columns.push_back(ColumnDefinition<Patch>("Harvest", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::HARVESTC); }));
		columns.push_back(ColumnDefinition<Patch>("LU_ch", [](const Patch &p) { return p.get_climate().gridcell.landcover.acflux_landuse_change; }));
		columns.push_back(ColumnDefinition<Patch>("Slow_h", [](const Patch &p) { return p.get_climate().gridcell.landcover.acflux_harvest_slow; }));
	}
	columns.push_back(ColumnDefinition<Patch>("NEE", get_patch_nee));

	return columns;
}

std::vector<ColumnDefinition<Patch>> nflux_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
		/// Carbon flux to atmosphere from burnt vegetation and li
		/// Flux to vegetation associated with manure addition (kgN/m2)
		columns.push_back(ColumnDefinition<Patch>("MANUREN", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::MANUREN); }));
		/// Flux to vegetation associated with N addition (kgN/m2) 
		columns.push_back(ColumnDefinition<Patch>("NFERT", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NFERT); }));
		/// Nitrogen flux from atmosphere to vegetation associated with establishment (kgN/m2)
		columns.push_back(ColumnDefinition<Patch>("ESTN", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::ESTN); }));
		/// Nitrogen flux to atmosphere from consumed harvested products (kgN/m2)
		columns.push_back(ColumnDefinition<Patch>("HARVESTN", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::HARVESTN); }));
		/// Nitrogen flux from atmosphere to vegetation associated with sowing (kgC/m2)
		columns.push_back(ColumnDefinition<Patch>("SEEDN", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::SEEDN); }));
		/// NH3 flux to atmosphere from fire
		columns.push_back(ColumnDefinition<Patch>("NH3_FIRE", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NH3_FIRE); }));
		/// NOx flux to atmosphere from fire
		columns.push_back(ColumnDefinition<Patch>("NOx_FIRE", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NOx_FIRE); }));
		/// N2O flux to atmosphere from fire
		columns.push_back(ColumnDefinition<Patch>("N2O_FIRE", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::N2O_FIRE); }));
		/// N2 flux to atmosphere from fire
		columns.push_back(ColumnDefinition<Patch>("N2_FIRE", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::N2_FIRE); }));
		//---- Soil N transformation -----
		/// NH3 flux from soil (ntransform)
		columns.push_back(ColumnDefinition<Patch>("NH3_SOIL", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NH3_SOIL); }));
		/// NO flux from soil (ntransform)
		columns.push_back(ColumnDefinition<Patch>("NO_SOIL", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NO_SOIL); }));
		/// N2O flux in soil (ntransform)
		columns.push_back(ColumnDefinition<Patch>("N2O_SOIL", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::N2O_SOIL); }));
		/// N2 flux from soil (ntransform)
		columns.push_back(ColumnDefinition<Patch>("N2_SOIL", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::N2_SOIL); }));
		/// Net nitrification (ntransform)
		columns.push_back(ColumnDefinition<Patch>("NET_NITRIF", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NET_NITRIF); }));
		/// Net denitrification (ntransform)
		columns.push_back(ColumnDefinition<Patch>("NET_DENITRIF", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NET_DENITRIF); }));
		/// Gross nitrification (ntransform)
		columns.push_back(ColumnDefinition<Patch>("GROSS_NITRIF", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::GROSS_NITRIF); }));
		/// Gross denitrification (ntransform)
		columns.push_back(ColumnDefinition<Patch>("GROSS_DENITRIF", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::GROSS_DENITRIF); }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> globfirm_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("fireprob", [](const Patch &p) { return p.fireprob; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> npool_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("VegN", get_nmass_patch));

	columns.push_back(ColumnDefinition<Patch>("StandingN", get_standing_n));
	columns.push_back(ColumnDefinition<Patch>("LitterN", get_litter_n));
	columns.push_back(ColumnDefinition<Patch>("SoilN", get_centuryn));

	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition<Patch>("HarvSlowN", get_nmass_harv_slowc));
	}
	columns.push_back(ColumnDefinition<Patch>("AvailN", [](const Patch &p) { return p.soil.get_nmass_avail(); }));

	// tbi: total

	return columns;
}

std::vector<ColumnDefinition<Patch>> npool_veg_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("leaf", get_nmass_patch));

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_lai_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("lai", [](const Individual &i) { return i.lai; }));
	columns.push_back(ColumnDefinition<Individual>("dnpp", [](const Individual &i) { return i.dnpp; }));
	columns.push_back(ColumnDefinition<Individual>("dgpp", [](const Individual &i) { return i.dgpp; }));
	columns.push_back(ColumnDefinition<Individual>("dresp", [](const Individual &i) { return i.dgpp - i.dnpp; }));
	columns.push_back(ColumnDefinition<Individual>("resp_growth", [](const Individual &i) { return i.resp_growth; }));
	columns.push_back(ColumnDefinition<Individual>("resp_maintenance", [](const Individual &i) { return i.resp_maintenance; }));
	columns.push_back(ColumnDefinition<Individual>("cstorage_fraction", [](const Individual &i) { return i.cmass_storage / i.cstorage_max; }));
	columns.push_back(ColumnDefinition<Individual>("rpc", [](const Individual &i) { return i.rpc; }));
	columns.push_back(ColumnDefinition<Individual>("height", [](const Individual &i) { return i.height; }));
	columns.push_back(ColumnDefinition<Individual>("diam", [](const Individual &i) { return i.diam; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_sap", [](const Individual &i) { return i.cmass_sap; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_heart", [](const Individual &i) { return i.cmass_heart; }));
	columns.push_back(ColumnDefinition<Individual>("cgrow", [](const Individual &i) { return i.cgrow; }));
	columns.push_back(ColumnDefinition<Individual>("cgrow_leaf", [](const Individual &i) { return i.cgrow_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_leaf", [](const Individual &i) { return i.cmass_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cturnover_leaf", [](const Individual &i) { return i.cturnover_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("leafCbalance", [](const Individual &i) { return i.cgrow_leaf - i.cturnover_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cgrow_sap", [](const Individual &i) { return i.cgrow_sap; }));
	columns.push_back(ColumnDefinition<Individual>("cturnover_sap", [](const Individual &i) { return i.cturnover_sap; }));
	columns.push_back(ColumnDefinition<Individual>("sapwoodCbalance", [](const Individual &i) { return i.cgrow_sap - i.cturnover_sap; }));
	columns.push_back(ColumnDefinition<Individual>("cgrow_root", [](const Individual &i) { return i.cgrow_root; }));
	columns.push_back(ColumnDefinition<Individual>("cturnover_root", [](const Individual &i) { return i.cturnover_root; }));
	columns.push_back(ColumnDefinition<Individual>("rootCbalance", [](const Individual &i) { return i.cgrow_root - i.cturnover_root; }));
	columns.push_back(ColumnDefinition<Individual>("cgrowrootminusleaf", [](const Individual &i) { return i.cgrow_root - i.cgrow_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cgrowminusleafturnover", [](const Individual &i) { return i.cgrow - i.cturnover_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("dnppminuscgrow", [](const Individual &i) { return i.dnpp - i.cgrow; }));
	columns.push_back(ColumnDefinition<Individual>("cstorage_dynam_14", [](const Individual &i) { return i.cstorage_dynam_14.sum(); }));
	columns.push_back(ColumnDefinition<Individual>("meanleafage", [](const Individual &i) { return i.meanleafage; }));
	columns.push_back(ColumnDefinition<Individual>("meansapage", [](const Individual &i) { return i.meansapage; }));
	columns.push_back(ColumnDefinition<Individual>("meanrootage", [](const Individual &i) { return i.meanrootage; }));
	columns.push_back(ColumnDefinition<Individual>("leaf_growth_demand", [](const Individual &i) { return i.absolute_alpha_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("root_growth_demand", [](const Individual &i) { return i.absolute_alpha_root; }));
	columns.push_back(ColumnDefinition<Individual>("sap_growth_demand", [](const Individual &i) { return i.absolute_alpha_sap; }));
	columns.push_back(ColumnDefinition<Individual>("density", [](const Individual &i) { return i.densindiv; }));
	// columns.push_back(ColumnDefinition<Individual>("fpar", [](const Individual &i) { return i.fpar; }));
	columns.push_back(ColumnDefinition<Individual>("latosa", [](const Individual &i) { return i.latosa; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> mwcont_upper_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("Jan", [](const Patch &p) { return p.soil.mwcont[0][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Feb", [](const Patch &p) { return p.soil.mwcont[1][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Mar", [](const Patch &p) { return p.soil.mwcont[2][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Apr", [](const Patch &p) { return p.soil.mwcont[3][0]; }));
	columns.push_back(ColumnDefinition<Patch>("May", [](const Patch &p) { return p.soil.mwcont[4][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Jun", [](const Patch &p) { return p.soil.mwcont[5][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Jul", [](const Patch &p) { return p.soil.mwcont[6][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Aug", [](const Patch &p) { return p.soil.mwcont[7][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Sep", [](const Patch &p) { return p.soil.mwcont[8][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Oct", [](const Patch &p) { return p.soil.mwcont[9][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Nov", [](const Patch &p) { return p.soil.mwcont[10][0]; }));
	columns.push_back(ColumnDefinition<Patch>("Dec", [](const Patch &p) { return p.soil.mwcont[11][0]; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> mwcont_lower_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("Jan", [](const Patch &p) { return p.soil.mwcont[0][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Feb", [](const Patch &p) { return p.soil.mwcont[1][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Mar", [](const Patch &p) { return p.soil.mwcont[2][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Apr", [](const Patch &p) { return p.soil.mwcont[3][1]; }));
	columns.push_back(ColumnDefinition<Patch>("May", [](const Patch &p) { return p.soil.mwcont[4][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Jun", [](const Patch &p) { return p.soil.mwcont[5][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Jul", [](const Patch &p) { return p.soil.mwcont[6][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Aug", [](const Patch &p) { return p.soil.mwcont[7][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Sep", [](const Patch &p) { return p.soil.mwcont[8][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Oct", [](const Patch &p) { return p.soil.mwcont[9][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Nov", [](const Patch &p) { return p.soil.mwcont[10][1]; }));
	columns.push_back(ColumnDefinition<Patch>("Dec", [](const Patch &p) { return p.soil.mwcont[11][1]; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> pet_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("pet", [](const Patch &p) { return p.apet; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> simfire_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("Biome", [](const Patch &p) { return (double)get_gridcell(p).simfire_biome; }));
	columns.push_back(ColumnDefinition<Patch>("MxNest", [](const Patch &p) { return get_gridcell(p).nesterov_max; }));
	columns.push_back(ColumnDefinition<Patch>("CurrentNesterov", [](const Patch &p) { return get_gridcell(p).nesterov_cur; }));
	columns.push_back(ColumnDefinition<Patch>("PopDens", [](const Patch &p) { return get_gridcell(p).pop_density; }));
	columns.push_back(ColumnDefinition<Patch>("AMxFApar", [](const Patch &p) { return get_gridcell(p).fapar_annual_max; }));
	columns.push_back(ColumnDefinition<Patch>("FireProb", [](const Patch &p) { return get_gridcell(p).simfire_annual_burned_area; }));
	columns.push_back(ColumnDefinition<Patch>("Region", [](const Patch &p) { return (double)get_gridcell(p).simfire_region; }));

	columns.push_back(ColumnDefinition<Patch>("PBiome", [](const Patch &p) { return get_p_biome(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition<Patch>("PFpar", [](const Patch &p) { return get_p_fpar(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition<Patch>("PMxNest", [](const Patch &p) { return get_p_nesterov(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition<Patch>("PPopDens", [](const Patch &p) { return get_p_popdens(get_gridcell(p)); }));

	columns.push_back(ColumnDefinition<Patch>("MFireChance", [](const Patch &p) { return get_gridcell(p).monthly_fire_risk[date.month]; }));
	columns.push_back(ColumnDefinition<Patch>("DFireChance", [](const Patch &p) { return get_gridcell(p).monthly_fire_risk[date.month] / (double)date.ndaymonth[date.month]; }));

	columns.push_back(ColumnDefinition<Patch>("avail_fuel", [](const Patch &p) { return p.avail_fuel; }));

	columns.push_back(ColumnDefinition<Patch>("SimfireBurntArea", [](const Patch &p) { return get_gridcell(p).simfire_annual_burned_area; }));
	columns.push_back(ColumnDefinition<Patch>("ABurntArea", [](const Patch &p) { return get_gridcell(p).annual_burned_area; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> fuel_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("FineLitter", get_fine_litter));
	columns.push_back(ColumnDefinition<Patch>("FineWood", get_fine_woody_debris));
	columns.push_back(ColumnDefinition<Patch>("CoarseWood", get_coarse_woody_debris));
	columns.push_back(ColumnDefinition<Patch>("GrassCMass", get_half_cmass_grass));
	columns.push_back(ColumnDefinition<Patch>("avail_fuel", [](const Patch &p) { return p.avail_fuel; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> coarse_woody_debris_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("LitterHeart", get_heartwood_litter));
	columns.push_back(ColumnDefinition<Patch>("som", get_coarse_woody_som));
	columns.push_back(ColumnDefinition<Patch>("k_tun_litter", [](const Patch &p) { return get_gridcell(p).k_tun_litter; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> met_year_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("met_year", [](const Patch &p) { return (double)p.get_climate().met_year; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> aet_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("evap", [](const Patch &p) { return p.evap; }));
	columns.push_back(ColumnDefinition<Patch>("transp", [](const Patch &p) { return p.transp; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> fevap_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("fevap", [](const Patch &p) { 
		
		double fpc_phen_total = 0.0;	// phenology-weighted FPC sum for patch

		const Vegetation& vegetation = p.vegetation;
		
		for (int i=0; i<vegetation.nobj; i++) {
			Individual& indiv = vegetation[i];
			fpc_phen_total += indiv.fpc_today();
		}
	
		double bare_ground = std::max(1.0 - fpc_phen_total, 0.0);

		return bare_ground; 
	}));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_cmass_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.sompool[SOILSTRUCT].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.sompool[SOILMICRO].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.sompool[SURFHUMUS].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.sompool[SURFMICRO].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.sompool[SOILMETA].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.sompool[SLOWSOM].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.sompool[PASSIVESOM].cmass; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_nmass_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.sompool[SOILSTRUCT].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.sompool[SOILMICRO].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.sompool[SURFHUMUS].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.sompool[SURFMICRO].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.sompool[SOILMETA].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.sompool[SLOWSOM].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.sompool[PASSIVESOM].nmass; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_nimmob_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.nimmob_pool[SURFSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.nimmob_pool[SOILSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.nimmob_pool[SOILMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.nimmob_pool[SURFHUMUS]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.nimmob_pool[SURFMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.nimmob_pool[SURFMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.nimmob_pool[SURFFWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.nimmob_pool[SURFCWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.nimmob_pool[SOILMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.nimmob_pool[SLOWSOM]; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.nimmob_pool[PASSIVESOM]; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_nmin_gross_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.nmin_gross_pool[SOILSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.nmin_gross_pool[SOILMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFHUMUS]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFFWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.nmin_gross_pool[SURFCWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.nmin_gross_pool[SOILMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.nmin_gross_pool[SLOWSOM]; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.nmin_gross_pool[PASSIVESOM]; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> ninput_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("dNO3dep", [](const Patch &p) { return get_gridcell(p).dNO3dep; }));
	columns.push_back(ColumnDefinition<Patch>("dNH4dep", [](const Patch &p) { return get_gridcell(p).dNH4dep; }));
	columns.push_back(ColumnDefinition<Patch>("nfert", [](const Patch &p) { return p.fluxes.get_daily_flux(Fluxes::NFERT, date.day); }));
	columns.push_back(ColumnDefinition<Patch>("total", [](const Patch &p) { return get_gridcell(p).dNO3dep + get_gridcell(p).dNH4dep; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> evap_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("evaporation", [](const Patch &p) { return p.evap; }));
	columns.push_back(ColumnDefinition<Patch>("interception", [](const Patch &p) { return p.intercep; }));
	columns.push_back(ColumnDefinition<Patch>("transpiration", [](const Patch &p) { return p.transp; }));
	columns.push_back(ColumnDefinition<Patch>("evapotranspiration", [](const Patch &p) { return p.evap + p.transp; }));
	return columns;
}

/*
for (int r = 0; r < NSOMPOOL; r++) {

	if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
		surfsoillitterc += patch.soil.sompool[r].cmass * to_gridcell_average;
		surfsoillittern += patch.soil.sompool[r].nmass * to_gridcell_average;
	}
	else if (r == SURFFWD || r == SURFCWD) {
		cwdc += patch.soil.sompool[r].cmass * to_gridcell_average;
		cwdn += patch.soil.sompool[r].nmass * to_gridcell_average;
	}
	else {
		centuryc += patch.soil.sompool[r].cmass * to_gridcell_average;
		centuryn += patch.soil.sompool[r].nmass * to_gridcell_average;
	}
}
*/

double get_soilc_surfsoillitterc(const Patch &patch) {
	double cmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
			cmass += patch.soil.sompool[r].cmass;
		}
	}
	return cmass;
}

double get_soilc_cwdc(const Patch &patch) {
	double cmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if (r == SURFFWD || r == SURFCWD) {
			cmass += patch.soil.sompool[r].cmass;
		}
	}
	return cmass;
}

double get_soilc_centuryc(const Patch &patch) {
	double cmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if(r != SURFMETA && r != SURFSTRUCT && r != SOILMETA && r != SOILSTRUCT && r != SURFFWD && r != SURFCWD){
			cmass += patch.soil.sompool[r].cmass;
		}
	}
	return cmass;
}

double get_soiln_surfsoillittern(const Patch &patch) {
	double nmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
			nmass += patch.soil.sompool[r].nmass;
		}
	}
	return nmass;
}

double get_soiln_cwdn(const Patch &patch) {
	double nmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if (r == SURFFWD || r == SURFCWD) {
			nmass += patch.soil.sompool[r].nmass;
		}
	}
	return nmass;
}

double get_soiln_centuryn(const Patch &patch) {
	double nmass = 0;
	for (int r = 0; r < NSOMPOOL; r++) {
		if(r != SURFMETA && r != SURFSTRUCT && r != SOILMETA && r != SOILSTRUCT && r != SURFFWD && r != SURFCWD){
			nmass += patch.soil.sompool[r].nmass;
		}
	}
	return nmass;
}

std::vector<ColumnDefinition<Patch>> soil_cmass_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("surfsoillitterc", get_soilc_surfsoillitterc));
	columns.push_back(ColumnDefinition<Patch>("cwdc", get_soilc_cwdc));
	columns.push_back(ColumnDefinition<Patch>("centuryc", get_soilc_centuryc));

	return columns;
}

std::vector<ColumnDefinition<Patch>> soil_nmass_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("surfsoillittern", get_soiln_surfsoillittern));
	columns.push_back(ColumnDefinition<Patch>("cwdn", get_soiln_cwdn));
	columns.push_back(ColumnDefinition<Patch>("centuryn", get_soiln_centuryn));
	columns.push_back(ColumnDefinition<Patch>("mineraln", [](const Patch &p) { return p.soil.NH4_mass + p.soil.NO3_mass; }));

	return columns;
}

double get_anuptake(const Patch &patch) {
	double nuptake = 0;
	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		const Individual &indiv = patch.vegetation[i];
		nuptake += indiv.anuptake - indiv.antoavail;
	}
	return nuptake;
}

double get_nuptake(const Patch &patch) {
	double nuptake = 0;
	for (size_t i = 0; i < patch.vegetation.nobj; i++) {
		nuptake += patch.vegetation[i].dnuptake - patch.vegetation[i].ntoavail;
	}
	return nuptake;
}

std::vector<ColumnDefinition<Patch>> aninput_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("aNO3dep", [](const Patch &p) { return get_gridcell(p).aNO3dep; }));
	columns.push_back(ColumnDefinition<Patch>("aNH4dep", [](const Patch &p) { return get_gridcell(p).aNH4dep; }));
	columns.push_back(ColumnDefinition<Patch>("nfert", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NFERT); }));
	columns.push_back(ColumnDefinition<Patch>("total", [](const Patch &p) { return get_gridcell(p).aNO3dep + get_gridcell(p).aNH4dep; }));
	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_cpool_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("cmass_leaf", [](const Individual &i) { return i.cmass_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_root", [](const Individual &i) { return i.cmass_root; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_crown", [](const Individual &i) { return i.cmass_crown; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_sap", [](const Individual &i) { return i.cmass_sap; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_heart", [](const Individual &i) { return i.cmass_heart; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_repr", [](const Individual &i) { return i.cmass_repr; }));
	columns.push_back(ColumnDefinition<Individual>("cmass_storage", [](const Individual &i) { return i.cmass_storage; }));

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_npool_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("nmass_leaf", [](const Individual &i) { return i.nmass_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_root", [](const Individual &i) { return i.nmass_root; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_crown", [](const Individual &i) { return i.nmass_crown; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_sap", [](const Individual &i) { return i.nmass_sap; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_heart", [](const Individual &i) { return i.nmass_heart; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_repr", [](const Individual &i) { return i.nmass_repr; }));
	columns.push_back(ColumnDefinition<Individual>("nmass_storage", [](const Individual &i) { return i.nmass_storage; }));

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_hydraulics_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("kmax", [](const Individual &i) { return i.kmax_drought_affected; }));
	columns.push_back(ColumnDefinition<Individual>("lai", [](const Individual &i) { return i.lai; }));
	columns.push_back(ColumnDefinition<Individual>("height", [](const Individual &i) { return i.height; }));
	columns.push_back(ColumnDefinition<Individual>("transp", [](const Individual &i) { return i.aet; }));
	columns.push_back(ColumnDefinition<Individual>("minimum_lwp", [](const Individual &i) { return min(i.lwp); }));
	columns.push_back(ColumnDefinition<Individual>("soil_water_potential", [](const Individual &i) { return i.weighted_psi; }));

	return columns;
}

std::vector<ColumnDefinition<Individual>> indiv_resp_columns() {
	std::vector<ColumnDefinition<Individual>> columns;

	columns.push_back(ColumnDefinition<Individual>("resp_leaf", [](const Individual &i) { return i.resp_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("resp_root", [](const Individual &i) { return i.resp_root; }));
	columns.push_back(ColumnDefinition<Individual>("resp_sap", [](const Individual &i) { return i.resp_sap; }));
	columns.push_back(ColumnDefinition<Individual>("resp_crown", [](const Individual &i) { return i.resp_crown; }));
	columns.push_back(ColumnDefinition<Individual>("resp_repr", [](const Individual &i) { return i.resp_repr; }));

	columns.push_back(ColumnDefinition<Individual>("resp_maintenance", [](const Individual &i) { return i.resp_maintenance; }));
	columns.push_back(ColumnDefinition<Individual>("resp_growth", [](const Individual &i) { return i.resp_growth; }));

	columns.push_back(ColumnDefinition<Individual>("total", [](const Individual &i) { return i.resp_maintenance + i.resp_growth; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> water_budget_cols_annual() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("precip", [](const Patch &p) { return p.get_climate().aprec; }));
	columns.push_back(ColumnDefinition<Patch>("evap", [](const Patch &p) { return p.aevap; }));
	columns.push_back(ColumnDefinition<Patch>("transp", [](const Patch &p) { return p.aaet; }));
	columns.push_back(ColumnDefinition<Patch>("interception", [](const Patch &p) { return p.aintercep; }));
	columns.push_back(ColumnDefinition<Patch>("surfrunoff", [](const Patch &p) { return p.asurfrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("drainrunoff", [](const Patch &p) { return p.adrainrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("baserunoff", [](const Patch &p) { return p.abaserunoff; }));
	columns.push_back(ColumnDefinition<Patch>("totrunoff", [](const Patch &p) { return p.arunoff; }));
	columns.push_back(ColumnDefinition<Patch>("flux", [](const Patch &p) { return p.get_climate().aprec - p.aevap - p.aaet - p.aintercep - p.arunoff; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> water_budget_cols_daily() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("precip", [](const Patch &p) { return p.get_climate().prec; }));
	columns.push_back(ColumnDefinition<Patch>("evap", [](const Patch &p) { return p.evap; }));
	columns.push_back(ColumnDefinition<Patch>("transp", [](const Patch &p) { return p.transp; }));
	columns.push_back(ColumnDefinition<Patch>("interception", [](const Patch &p) { return p.intercep; }));
	columns.push_back(ColumnDefinition<Patch>("surfrunoff", [](const Patch &p) { return p.soil.surfrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("drainrunoff", [](const Patch &p) { return p.soil.drainrunoff; }));
	columns.push_back(ColumnDefinition<Patch>("baserunoff", [](const Patch &p) { return p.soil.baserunoff; }));
	columns.push_back(ColumnDefinition<Patch>("totrunoff", [](const Patch &p) { return p.soil.runoff; }));
	columns.push_back(ColumnDefinition<Patch>("flux", [](const Patch &p) { return p.get_climate().prec - p.evap - p.transp - p.intercep - p.soil.runoff; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> nitrogen_budget_cols_annual() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("deposition", [](const Patch &p) { return p.stand.get_gridcell().aNO3dep + p.stand.get_gridcell().aNH4dep; }));
	columns.push_back(ColumnDefinition<Patch>("fixation", [](const Patch &p) { return p.soil.anfix; }));
	columns.push_back(ColumnDefinition<Patch>("uptake", get_anuptake));
	columns.push_back(ColumnDefinition<Patch>("net_mineralisation", [](const Patch &p) { return p.soil.anmin - p.soil.animmob; }));
	columns.push_back(ColumnDefinition<Patch>("gross_mineralisation", [](const Patch &p) { return p.soil.anmin; }));
	columns.push_back(ColumnDefinition<Patch>("immobilisation", [](const Patch &p) { return p.soil.animmob; }));
	columns.push_back(ColumnDefinition<Patch>("min_leach", [](const Patch &p) { return p.soil.aminleach; }));
	columns.push_back(ColumnDefinition<Patch>("org_leach", [](const Patch &p) { return p.soil.aorgNleach; }));
	columns.push_back(ColumnDefinition<Patch>("volatilisation", [](const Patch &p) { return p.soil.anvol; }));
	columns.push_back(ColumnDefinition<Patch>("gas_emissions", [](const Patch &p) { return p.soil.anloss_gas; }));
	// Annual flux in available N: deposition + fixation + mineralisation - uptake - leaching - volatilisation - gas emissions
	columns.push_back(ColumnDefinition<Patch>("flux_avail", [](const Patch &p) { return p.stand.get_gridcell().aNO3dep + p.stand.get_gridcell().aNH4dep + p.soil.anfix + p.soil.anmin - p.soil.animmob - get_anuptake(p) - p.soil.aminleach - p.soil.anvol - p.soil.anloss_gas; }));
	// Annual flux in total N: deposition + fixation - gas emissions
	columns.push_back(ColumnDefinition<Patch>("flux", [](const Patch &p) { return p.stand.get_gridcell().aNO3dep + p.stand.get_gridcell().aNH4dep + p.soil.anfix - p.soil.anloss_gas; }));
	columns.push_back(ColumnDefinition<Patch>("avail", [](const Patch &p) { return p.soil.get_nmass_avail(); }));
	columns.push_back(ColumnDefinition<Patch>("totaln", [](const Patch &p) { return p.ncont(); }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> nitrogen_budget_cols_daily() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("deposition", [](const Patch &p) { return p.stand.get_gridcell().dNO3dep + p.stand.get_gridcell().dNH4dep; }));
	columns.push_back(ColumnDefinition<Patch>("fixation", [](const Patch &p) { return p.soil.dnfix; }));
	columns.push_back(ColumnDefinition<Patch>("uptake", get_nuptake));
	columns.push_back(ColumnDefinition<Patch>("net_mineralisation", [](const Patch &p) { return p.soil.net_mineralisation; }));
	columns.push_back(ColumnDefinition<Patch>("gross_mineralisation", [](const Patch &p) { return p.soil.nmin_gross; }));
	columns.push_back(ColumnDefinition<Patch>("immobilisation", [](const Patch &p) { return p.soil.nimmob; }));
	columns.push_back(ColumnDefinition<Patch>("min_leach", [](const Patch &p) { return p.soil.dminleach; }));
	columns.push_back(ColumnDefinition<Patch>("org_leach", [](const Patch &p) { return p.soil.dorgNleach; }));
	columns.push_back(ColumnDefinition<Patch>("volatilisation", [](const Patch &p) { return p.soil.nvol; }));
	columns.push_back(ColumnDefinition<Patch>("gas_emissions", [](const Patch &p) { return p.soil.nloss_gas; }));
	// flux in available N
	columns.push_back(ColumnDefinition<Patch>("flux_avail", [](const Patch &p) { return p.stand.get_gridcell().dNO3dep + p.stand.get_gridcell().dNH4dep + p.soil.dnfix + p.soil.net_mineralisation - get_nuptake(p) - p.soil.dminleach - p.soil.nvol - p.soil.nloss_gas; }));
	// flux in total N
	columns.push_back(ColumnDefinition<Patch>("flux", [](const Patch &p) { return p.stand.get_gridcell().dNO3dep + p.stand.get_gridcell().dNH4dep + p.soil.dnfix - p.soil.nloss_gas; }));
	columns.push_back(ColumnDefinition<Patch>("avail", [](const Patch &p) { return p.soil.get_nmass_avail(); }));
	columns.push_back(ColumnDefinition<Patch>("totaln", [](const Patch &p) { return p.ncont(); }));
	return columns;
}

std::vector<ColumnDefinition<Individual>> mortality_cols() {
	std::vector<ColumnDefinition<Individual>> columns;
	columns.push_back(ColumnDefinition<Individual>("mort_age", [](const Individual &i) { return i.mort_age; }));
	columns.push_back(ColumnDefinition<Individual>("mort_greff", [](const Individual &i) { return i.mort_greff; }));
	columns.push_back(ColumnDefinition<Individual>("mort_cstarve", [](const Individual &i) { return i.mort_cstarve; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_decay_reduction_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.decay_reduction[SURFSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.decay_reduction[SOILSTRUCT]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.decay_reduction[SOILMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.decay_reduction[SURFHUMUS]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.decay_reduction[SURFMICRO]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.decay_reduction[SURFMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.decay_reduction[SURFFWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.decay_reduction[SURFCWD]; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.decay_reduction[SOILMETA]; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.decay_reduction[SLOWSOM]; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.decay_reduction[PASSIVESOM]; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_frac_remaining_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.sompool[SOILSTRUCT].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.sompool[SOILMICRO].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.sompool[SURFHUMUS].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.sompool[SURFMICRO].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.sompool[SOILMETA].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.sompool[SLOWSOM].fracremain; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.sompool[PASSIVESOM].fracremain; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> sompool_ntoc_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SOILSTRUCT", [](const Patch &p) { return p.soil.sompool[SOILSTRUCT].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMICRO", [](const Patch &p) { return p.soil.sompool[SOILMICRO].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SURFHUMUS", [](const Patch &p) { return p.soil.sompool[SURFHUMUS].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMICRO", [](const Patch &p) { return p.soil.sompool[SURFMICRO].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SOILMETA", [](const Patch &p) { return p.soil.sompool[SOILMETA].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("SLOWSOM", [](const Patch &p) { return p.soil.sompool[SLOWSOM].ntoc; }));
	columns.push_back(ColumnDefinition<Patch>("PASSIVESOM", [](const Patch &p) { return p.soil.sompool[PASSIVESOM].ntoc; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> clitter_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].cmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].cmass; }));
	return columns;
}


std::vector<ColumnDefinition<Patch>> nlitter_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("SURFSTRUCT", [](const Patch &p) { return p.soil.sompool[SURFSTRUCT].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFMETA", [](const Patch &p) { return p.soil.sompool[SURFMETA].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFFWD", [](const Patch &p) { return p.soil.sompool[SURFFWD].nmass; }));
	columns.push_back(ColumnDefinition<Patch>("SURFCWD", [](const Patch &p) { return p.soil.sompool[SURFCWD].nmass; }));
	return columns;
}

double get_nindiv(const Patch &patch, const Pft &pft) {
	double n = 0;
	const Vegetation &vegetation = patch.vegetation;
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		if (indiv.pft.id == pft.id) {
			n += 1;
		}
	}
	return n;
}

std::vector<ColumnDefinition<Gridcell>> get_gridcell_patchpft_columns(
	double (*read_value)(const Patch &, const Pft &)
) {
	std::vector<ColumnDefinition<Gridcell>> columns;
	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];
		std::string name = std::string((const char*)pft.name);
		columns.push_back(AggregatedPatchPftGridcellColumnDefinition(
			name,
			pft,
			read_value
		));
	}
	return columns;
}

std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_indiv_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t),
	bool (*pft_filter)(const Pft &)
) {
	std::vector<ColumnDefinition<Gridcell>> columns;
	for (size_t i = 0; i < pftlist.nobj; i++) {
		const Pft &pft = pftlist[i];
		if (!pft_filter(pft)) {
			continue;
		}
		const char *name = (const char *)pft.name;

		columns.push_back(AggregatedIndivGridcellColumnDefinition(
			name,
			pft,
			read_value,
			aggregate
		));
	}
	return columns;
}

std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_indiv_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t)
) {
	return get_gridcell_aggregated_indiv_columns(
		read_value,
		aggregate,
		[](const Pft &p) { return true; }
	);
}

std::vector<ColumnDefinition<Gridcell>> get_gridcell_aggregated_tree_columns(
    double (*read_value)(const Individual &),
    double (*aggregate)(const double *, size_t)
) {
	return get_gridcell_aggregated_indiv_columns(
		read_value,
		aggregate,
		[](const Pft &p) { return p.lifeform == TREE; }
	);
}

std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> ps_input_cols() {
	return {
		{"gdark",        [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].gdark; }},
		{"lai",          [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].lai; }},
		{"fphen",        [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].fphen_activity; }},
		{"ps_sink_lim",  [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ps_sink_lim; }},
		{"vcmax25",      [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].vcmax25; }},
		{"vpmax25",      [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].vpmax25; }},
		{"jmax25",       [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].jmax25; }},
		{"nscal",        [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].nscal; }},
		{"wscal",        [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].wscal; }},
		{"tleaf",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].tleaf; }},
		{"ha_v",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ha_v; }},
		{"ha_j",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ha_j; }},
		{"ds_v",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ds_v; }},
		{"ds_j",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ds_j; }},
		{"par",          [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].par; }},
		{"apar",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].apar(); }},
		{"ipar",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ipar(); }},
		{"fipar",         [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].fipar; }},
		{"vpd",          [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].vpd; }},
		{"ca",           [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].ca; }},
		{"X",            [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].X; }},
		{"stomata_open", [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].stomata_open ? 1.0 : 0.0; }},
		{"predawn_lwp",  [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].lwp; }},
		{"s_f",          [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].s_f; }},
		{"psi_f",        [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].psi_f; }},
		{"tgrowth",      [](const Individual &i, int t, int l) { return i.ps_inputs_layerwise[t][l].tgrowth; }},
	};
}

std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> ps_output_cols() {
	return {
		{"rd",    [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].rd; }},
		{"an",    [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].an; }},
		{"anc",   [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].anc; }},
		{"anj",   [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].anj; }},
		{"anp",   [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].anp; }},
		{"gsc",   [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].gsc; }},
		{"ci",    [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].ci; }},
		{"tleaf", [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].tleaf; }},
		{"vcmax", [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].vcmax; }},
		{"jmax",  [](const Individual &i, int t, int l) { return i.ps_outputs_layerwise[t][l].jmax; }},
	};
}

std::vector<std::pair<std::string, double (*)(const Individual&, int, int)>> energy_balance_cols() {
	return {
		{"converged",    [](const Individual &i, int t, int l) { return i.energy_balance_layerwise[t][l].get_converged() ? 1.0 : 0.0; }},
		{"delta_tleaf",  [](const Individual &i, int t, int l) { return i.energy_balance_layerwise[t][l].get_delta_tleaf(); }},
		{"niter",        [](const Individual &i, int t, int l) { return i.energy_balance_layerwise[t][l].get_niter() * 1.0; }},
		{"tleaf",        [](const Individual &i, int t, int l) { return i.energy_balance_layerwise[t][l].get_tleaf(); }},
		{"tair",         [](const Individual &i, int t, int l) { return i.vegetation.patch.get_climate().temps[t]; }},
	};
}

} // namespace GuessOutput
