#include "daveoutputhelpers.h"
#include "guess.h"
#include "commonoutput.h"
#include "simfire.h"
#include "blaze.h"

namespace GuessOutput {

void ensure_channel_open() {
	if (!output_channel) {
		fail("Output channel is NULL");
	}
}

double get_lai(const Individual &indiv) {
	return indiv.lai;
}

double get_er(const Patch &patch, int i) {
	double ra_kg = patch.fluxes.get_daily_flux(Fluxes::RA, date.day);
	double rh_kg = patch.fluxes.get_daily_flux(Fluxes::SOILC, date.day);
	return (ra_kg + rh_kg) * G_PER_KG;
}

double get_gpp(const Patch &patch, int i) {
	double gpp_kg = patch.fluxes.get_daily_flux(Fluxes::GPP, date.day);
	return gpp_kg * G_PER_KG;
}

double get_npp(const Patch &patch, int i) {
	double npp_kg = patch.fluxes.get_daily_flux(Fluxes::NPP, date.day);
	return npp_kg * G_PER_KG;
}

double get_nee(const Patch &patch, int i) {
	double rh_kg = patch.fluxes.get_daily_flux(Fluxes::SOILC, date.day);
	double npp_kg = patch.fluxes.get_daily_flux(Fluxes::NPP, date.day);
	return (rh_kg - npp_kg) * G_PER_KG;
}

double get_et(const Patch &patch, int i) {
	const Vegetation &vegetation = patch.vegetation;
	double transpiration = 0;
	for (int i = 0; i < vegetation.nobj; i++)
		transpiration += vegetation[i].aet;
	return transpiration + patch.soil.get_layer_soil_water_evap();
}

std::vector<Individual *> get_indivs(const Vegetation &vegetation, const char *name) {
	std::vector<Individual *> result;
	for (int i = 0; i < vegetation.nobj; i++) {
		Individual &indiv = vegetation[i];

		if (indiv.pft.name == name && indiv.alive) {
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

Individual *get_first_indiv(const Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;
	if (vegetation.nobj > 0) {
		return &(vegetation[0]);
	}
	return NULL;
}

int get_nlayer_fpar(const Patch &patch) {
	size_t nlayer = 0;
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		nlayer = std::max(nlayer, patch.vegetation[i].fpar_layerwise.size());
	}
	return nlayer;
}

double get_fpar_layer(const Individual &indiv, int i) {
	// i = 0 means bottom layer (so last element in the array), i = 1 means
	// second layer, etc.
	int index = indiv.fpar_layerwise.size() - i - 1;
	if (index < 0 || index > indiv.fpar_layerwise.size()) {
		return 0;
	}
	return indiv.fpar_layerwise[i];
}

int get_nlayer_lai(const Patch &patch) {
	size_t nlayer = 0;
	for (int i = 0; i < patch.vegetation.nobj; i++) {
		nlayer = std::max(nlayer, patch.vegetation[i].lai_layerwise.size());
	}
	return nlayer;
}

double get_lai_layer(const Individual &indiv, int i) {
	// i = 0 means bottom layer (so last element in the array), i = 1 means
	// second layer, etc.
	int index = indiv.lai_layerwise.size() - i - 1;
	if (index < 0 || index > indiv.lai_layerwise.size()) {
		return 0;
	}
	return indiv.lai_layerwise[i];
}

double get_pft_cmass_litter(const Patch &patch, const Pft &pft) {
	const Patchpft &patchpft = patch.pft[pft.id];
	return patchpft.cmass_litter_leaf + patchpft.cmass_litter_root + patchpft.cmass_litter_sap + patchpft.cmass_litter_heart + patchpft.cmass_litter_repr;
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

double sum(double *array, int n) {
	double total = 0;
	for (int i = 0; i < n; i++) {
		total += array[i];
	}
	return total;
}

double mean(double *array, int n) {
	return divide(sum(array, n), n);
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

std::vector<ColumnDefinition> runoff_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("Surf", [](const Patch &p) { return p.asurfrunoff; }));
	columns.push_back(ColumnDefinition("Drain", [](const Patch &p) { return p.adrainrunoff; }));
	columns.push_back(ColumnDefinition("Base", [](const Patch &p) { return p.abaserunoff; }));
	columns.push_back(ColumnDefinition("Total", [](const Patch &p) { return p.asurfrunoff + p.adrainrunoff + p.abaserunoff; }));

	return columns;
}

double get_cmass_patch(const Patch &patch) {
	double ccont = 0;
	for (unsigned int i = 0; i < patch.vegetation.nobj; i++) {
		ccont += patch.vegetation[i].ccont();
	}
	return ccont;
}

// clitter_gridcell <- standpft_clitter
double get_cmass_litter(const Patch &patch) {
	double cmass_litter = 0;

	for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		cmass_litter += get_pft_cmass_litter(patch, pftlist.getobj());
	}

	return cmass_litter;
}

double get_cmass_century(const Patch &patch) {
	double centuryc = 0;
	for (int i = 0; i < NSOMPOOL; i++) {
		if (i != SURFMETA && i != SURFSTRUCT && i != SOILMETA && i != SOILSTRUCT && i != SURFFWD && i != SURFCWD)
		centuryc += patch.soil.sompool[i].cmass;
		// centuryn += patch.soil.sompool[r].nmass * to_gridcell_average;
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

double get_surf_soil_litter_c(const Patch &patch) {
	double litter_c = 0;
	litter_c += patch.soil.sompool[SURFMETA].cmass;
	litter_c += patch.soil.sompool[SURFSTRUCT].cmass;
	litter_c += patch.soil.sompool[SOILMETA].cmass;
	litter_c += patch.soil.sompool[SOILSTRUCT].cmass;
	return litter_c;
}

double get_surf_soil_litter_n(const Patch &patch) {
	double litter_n = 0;
	litter_n += patch.soil.sompool[SURFMETA].nmass;
	litter_n += patch.soil.sompool[SURFSTRUCT].nmass;
	litter_n += patch.soil.sompool[SOILMETA].nmass;
	litter_n += patch.soil.sompool[SOILSTRUCT].nmass;
	return litter_n;
}

double get_cwdc(const Patch &patch) {
	double cwdc = 0;
	cwdc += patch.soil.sompool[SURFFWD].cmass;
	cwdc += patch.soil.sompool[SURFCWD].cmass;
	return cwdc;
}

double get_cwdn(const Patch &patch) {
	double cwdn = 0;
	cwdn += patch.soil.sompool[SURFFWD].nmass;
	cwdn += patch.soil.sompool[SURFCWD].nmass;
	return cwdn;
}

double get_littern(const Patch &patch) {
	return get_surf_soil_litter_n(patch) + get_cwdn(patch);
}

double get_litter_c(const Patch &patch) {
	double veg_c = get_cmass_litter(patch);
	double surfsoillitterc = get_surf_soil_litter_c(patch);
	double cwdc = get_cwdc(patch);
	return veg_c + surfsoillitterc + cwdc;
}

double get_nmass_litter(const Patch &patch) {
	double nmass_litter = 0;
	for (pftlist.firstobj(); pftlist.isobj; pftlist.nextobj()) {
		const Patchpft &ppft = patch.pft[pftlist.getobj().id];
		nmass_litter += ppft.nmass_litter_leaf + ppft.nmass_litter_root + ppft.nmass_litter_sap + ppft.nmass_litter_heart;
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
	double elems[patch.pft.nobj];
	for (int i = 0; i < patch.pft.nobj; i++) {
		const Patchpft &ppft = patch.pft[i];
		elems[i] = read_ppft(ppft);
	}
	return aggregate(elems, patch.pft.nobj);
}

/*
Aggregate a value across all individuals in the patch.
*/
double aggregate_indiv(const Patch &patch, double (*read_indiv)(const Individual &), double (*aggregate)(const double*, int)) {
	const Vegetation &vegetation = patch.vegetation;
	double elems[vegetation.nobj];
	for (int i = 0; i < vegetation.nobj; i++) {
		const Individual &indiv = vegetation[i];
		elems[i] = read_indiv(indiv);
	}
	return aggregate(elems, vegetation.nobj);
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

std::vector<ColumnDefinition> cpool_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("VegC", get_cmass_patch));

	if (!ifcentury) {
		columns.push_back(ColumnDefinition("LittC", get_litter_c));
		columns.push_back(ColumnDefinition("SoilfC", [](const Patch &p) { return p.soil.cpool_fast; }));
		columns.push_back(ColumnDefinition("SoilsC", [](const Patch &p) { return p.soil.cpool_slow; }));
	} else {
		columns.push_back(ColumnDefinition("LitterC", get_litter_c));
		columns.push_back(ColumnDefinition("SoilC", get_cmass_century));
	}
	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition("HarvSlowC", get_cmass_harv_slowc));
	}

	return columns;
}

std::vector<ColumnDefinition> cflux_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("Veg", [](const Patch &p) { return -p.fluxes.get_annual_flux(Fluxes::NPP); }));
	columns.push_back(ColumnDefinition("Repr", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::REPRC); }));
	columns.push_back(ColumnDefinition("Soil", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::SOILC) + p.soil.aorgCleach; }));
	columns.push_back(ColumnDefinition("Fire", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::FIREC); }));
	columns.push_back(ColumnDefinition("Est", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::ESTC); }));
	if (run_landcover) {
		columns.push_back(ColumnDefinition("Seed", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::SEEDC); }));
		columns.push_back(ColumnDefinition("Harvest", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::HARVESTC); }));
		columns.push_back(ColumnDefinition("LU_ch", [](const Patch &p) { return p.get_climate().gridcell.landcover.acflux_landuse_change; }));
		columns.push_back(ColumnDefinition("Slow_h", [](const Patch &p) { return p.get_climate().gridcell.landcover.acflux_harvest_slow; }));
	}
	columns.push_back(ColumnDefinition("NEE", get_patch_nee));

	return columns;
}

std::vector<ColumnDefinition> globfirm_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("fireprob", [](const Patch &p) { return p.fireprob; }));

	return columns;
}

std::vector<ColumnDefinition> npool_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("VegN", get_nmass_patch));

	columns.push_back(ColumnDefinition("LitterN", get_littern));
	columns.push_back(ColumnDefinition("SoilN", get_centuryn));

	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition("HarvSlowN", get_nmass_harv_slowc));
	}
	// tbi: total

	return columns;
}
std::vector<ColumnDefinition> mwcont_upper_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("Jan", [](const Patch &p) { return p.soil.mwcont[0][0]; }));
	columns.push_back(ColumnDefinition("Feb", [](const Patch &p) { return p.soil.mwcont[1][0]; }));
	columns.push_back(ColumnDefinition("Mar", [](const Patch &p) { return p.soil.mwcont[2][0]; }));
	columns.push_back(ColumnDefinition("Apr", [](const Patch &p) { return p.soil.mwcont[3][0]; }));
	columns.push_back(ColumnDefinition("May", [](const Patch &p) { return p.soil.mwcont[4][0]; }));
	columns.push_back(ColumnDefinition("Jun", [](const Patch &p) { return p.soil.mwcont[5][0]; }));
	columns.push_back(ColumnDefinition("Jul", [](const Patch &p) { return p.soil.mwcont[6][0]; }));
	columns.push_back(ColumnDefinition("Aug", [](const Patch &p) { return p.soil.mwcont[7][0]; }));
	columns.push_back(ColumnDefinition("Sep", [](const Patch &p) { return p.soil.mwcont[8][0]; }));
	columns.push_back(ColumnDefinition("Oct", [](const Patch &p) { return p.soil.mwcont[9][0]; }));
	columns.push_back(ColumnDefinition("Nov", [](const Patch &p) { return p.soil.mwcont[10][0]; }));
	columns.push_back(ColumnDefinition("Dec", [](const Patch &p) { return p.soil.mwcont[11][0]; }));

	return columns;
}

std::vector<ColumnDefinition> mwcont_lower_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("Jan", [](const Patch &p) { return p.soil.mwcont[0][1]; }));
	columns.push_back(ColumnDefinition("Feb", [](const Patch &p) { return p.soil.mwcont[1][1]; }));
	columns.push_back(ColumnDefinition("Mar", [](const Patch &p) { return p.soil.mwcont[2][1]; }));
	columns.push_back(ColumnDefinition("Apr", [](const Patch &p) { return p.soil.mwcont[3][1]; }));
	columns.push_back(ColumnDefinition("May", [](const Patch &p) { return p.soil.mwcont[4][1]; }));
	columns.push_back(ColumnDefinition("Jun", [](const Patch &p) { return p.soil.mwcont[5][1]; }));
	columns.push_back(ColumnDefinition("Jul", [](const Patch &p) { return p.soil.mwcont[6][1]; }));
	columns.push_back(ColumnDefinition("Aug", [](const Patch &p) { return p.soil.mwcont[7][1]; }));
	columns.push_back(ColumnDefinition("Sep", [](const Patch &p) { return p.soil.mwcont[8][1]; }));
	columns.push_back(ColumnDefinition("Oct", [](const Patch &p) { return p.soil.mwcont[9][1]; }));
	columns.push_back(ColumnDefinition("Nov", [](const Patch &p) { return p.soil.mwcont[10][1]; }));
	columns.push_back(ColumnDefinition("Dec", [](const Patch &p) { return p.soil.mwcont[11][1]; }));

	return columns;
}

std::vector<ColumnDefinition> pet_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("pet", [](const Patch &p) { return p.apet; }));

	return columns;
}

std::vector<ColumnDefinition> simfire_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("Biome", [](const Patch &p) { return (double)get_gridcell(p).simfire_biome; }));
	columns.push_back(ColumnDefinition("MxNest", [](const Patch &p) { return get_gridcell(p).nesterov_max; }));
	columns.push_back(ColumnDefinition("CurrentNesterov", [](const Patch &p) { return get_gridcell(p).nesterov_cur; }));
	columns.push_back(ColumnDefinition("PopDens", [](const Patch &p) { return get_gridcell(p).pop_density; }));
	columns.push_back(ColumnDefinition("AMxFApar", [](const Patch &p) { return get_gridcell(p).fapar_annual_max; }));
	columns.push_back(ColumnDefinition("FireProb", [](const Patch &p) { return get_gridcell(p).simfire_annual_burned_area; }));
	columns.push_back(ColumnDefinition("Region", [](const Patch &p) { return (double)get_gridcell(p).simfire_region; }));

	columns.push_back(ColumnDefinition("PBiome", [](const Patch &p) { return get_p_biome(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition("PFpar", [](const Patch &p) { return get_p_fpar(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition("PMxNest", [](const Patch &p) { return get_p_nesterov(get_gridcell(p)); }));
	columns.push_back(ColumnDefinition("PPopDens", [](const Patch &p) { return get_p_popdens(get_gridcell(p)); }));

	columns.push_back(ColumnDefinition("MFireChance", [](const Patch &p) { return get_gridcell(p).monthly_fire_risk[date.month]; }));
	columns.push_back(ColumnDefinition("DFireChance", [](const Patch &p) { return get_gridcell(p).monthly_fire_risk[date.month] / (double)date.ndaymonth[date.month]; }));

	columns.push_back(ColumnDefinition("avail_fuel", [](const Patch &p) { return p.avail_fuel; }));

	columns.push_back(ColumnDefinition("SimfireBurntArea", [](const Patch &p) { return get_gridcell(p).simfire_annual_burned_area; }));
	columns.push_back(ColumnDefinition("ABurntArea", [](const Patch &p) { return get_gridcell(p).annual_burned_area; }));

	return columns;
}

std::vector<ColumnDefinition> fuel_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("FineLitter", get_fine_litter));
	columns.push_back(ColumnDefinition("FineWood", get_fine_woody_debris));
	columns.push_back(ColumnDefinition("CoarseWood", get_coarse_woody_debris));
	columns.push_back(ColumnDefinition("GrassCMass", get_half_cmass_grass));
	columns.push_back(ColumnDefinition("avail_fuel", [](const Patch &p) { return p.avail_fuel; }));

	return columns;
}

std::vector<ColumnDefinition> coarse_woody_debris_columns() {
	std::vector<ColumnDefinition> columns;

	columns.push_back(ColumnDefinition("LitterHeart", get_heartwood_litter));
	columns.push_back(ColumnDefinition("som", get_coarse_woody_som));
	columns.push_back(ColumnDefinition("k_tun_litter", [](const Patch &p) { return get_gridcell(p).k_tun_litter; }));

	return columns;
}

std::vector<ColumnDefinition> met_year_columns() {
	std::vector<ColumnDefinition> columns;
	columns.push_back(ColumnDefinition("met_year", [](const Patch &p) { return (double)p.get_climate().met_year; }));
	return columns;
}

} // namespace GuessOutput
