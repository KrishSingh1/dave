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

double get_et(const Patch &patch) {
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

Individual *get_first_indiv(const Patch &patch) {
	const Vegetation &vegetation = patch.vegetation;
	if (vegetation.nobj > 0) {
		return &(vegetation[0]);
	}
	return NULL;
}

int get_nlayer_fpar(const Patch &patch) {
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->fipar_layerwise.size();
	}
	// No individuals are established in this patch.
	return 0;
}


int get_nlayer_lai(const Patch &patch) {
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->lai_layerwise.size();
	}
	return 0;
}

int get_nlayer_an(const Patch &patch) {
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->an_layerwise.size();
	}
	return 0;
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
	columns.push_back(ColumnDefinition<Patch>("Total", [](const Patch &p) { return p.asurfrunoff + p.adrainrunoff + p.abaserunoff; }));

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
	return aggregate(&elems.front(), patch.pft.nobj);
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
		columns.push_back(ColumnDefinition<Patch>("LittC", get_litter_c));
		columns.push_back(ColumnDefinition<Patch>("SoilfC", [](const Patch &p) { return p.soil.cpool_fast; }));
		columns.push_back(ColumnDefinition<Patch>("SoilsC", [](const Patch &p) { return p.soil.cpool_slow; }));
	} else {
		columns.push_back(ColumnDefinition<Patch>("LitterC", get_litter_c));
		columns.push_back(ColumnDefinition<Patch>("SoilC", get_cmass_century));
	}
	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition<Patch>("HarvSlowC", get_cmass_harv_slowc));
	}

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

std::vector<ColumnDefinition<Patch>> globfirm_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("fireprob", [](const Patch &p) { return p.fireprob; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> npool_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("VegN", get_nmass_patch));

	columns.push_back(ColumnDefinition<Patch>("LitterN", get_littern));
	columns.push_back(ColumnDefinition<Patch>("SoilN", get_centuryn));

	if (run_landcover && ifslowharvestpool) {
		columns.push_back(ColumnDefinition<Patch>("HarvSlowN", get_nmass_harv_slowc));
	}
	// tbi: total

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
	columns.push_back(ColumnDefinition<Individual>("leaf_growth_demand", [](const Individual &i) { return i.relative_alpha_leaf; }));
	columns.push_back(ColumnDefinition<Individual>("root_growth_demand", [](const Individual &i) { return i.relative_alpha_root; }));
	columns.push_back(ColumnDefinition<Individual>("sap_growth_demand", [](const Individual &i) { return i.relative_alpha_sap; }));
	columns.push_back(ColumnDefinition<Individual>("repr_growth_demand", [](const Individual &i) { return i.relative_alpha_repr; }));
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

/*
Sum a quantity across all SOM pools.
*/
double sum_sompool(const Patch &patch, double (*read_pool)(const Sompool &pool)) {
	double total = 0;
	for (size_t i = 0; i < NSOMPOOL; i++) {
		total += read_pool(patch.soil.sompool[i]);
	}
	return total;
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
	columns.push_back(ColumnDefinition<Patch>("total", [](const Patch &p) { return sum_sompool(p, [](const Sompool &pool) { return pool.cmass; }); }));
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
	columns.push_back(ColumnDefinition<Patch>("total", [](const Patch &p) { return sum_sompool(p, [](const Sompool &pool) { return pool.cmass; }); }));
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

std::vector<ColumnDefinition<Patch>> soil_nflux_columns() {
	std::vector<ColumnDefinition<Patch>> columns;

	columns.push_back(ColumnDefinition<Patch>("nmin_net", [](const Patch &p) { return p.soil.net_mineralisation; }));
	columns.push_back(ColumnDefinition<Patch>("nmin_leach", [](const Patch &p) { return p.soil.dminleach; }));
	columns.push_back(ColumnDefinition<Patch>("norg_leach", [](const Patch &p) { return p.soil.dorgNleach; }));
	columns.push_back(ColumnDefinition<Patch>("nvolatilisation", [](const Patch &p) { return p.soil.nvol; }));
	columns.push_back(ColumnDefinition<Patch>("nloss_gas", [](const Patch &p) { return p.soil.nloss_gas; }));

	return columns;
}

std::vector<ColumnDefinition<Patch>> aninput_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	columns.push_back(ColumnDefinition<Patch>("aNO3dep", [](const Patch &p) { return get_gridcell(p).aNO3dep; }));
	columns.push_back(ColumnDefinition<Patch>("aNH4dep", [](const Patch &p) { return get_gridcell(p).aNH4dep; }));
	columns.push_back(ColumnDefinition<Patch>("nfert", [](const Patch &p) { return p.fluxes.get_annual_flux(Fluxes::NFERT); }));
	columns.push_back(ColumnDefinition<Patch>("total", [](const Patch &p) { return get_gridcell(p).aNO3dep + get_gridcell(p).aNH4dep; }));
	return columns;
}

std::vector<ColumnDefinition<Patch>> establishment_columns() {
	std::vector<ColumnDefinition<Patch>> columns;
	// columns.push_back(ColumnDefinition<Patch>("anetps_ff_frac", [](const Patch &p) { return ; }));
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
	columns.push_back(ColumnDefinition<Individual>("minimum_lwp", [](const Individual &i) { return min(i.lwp); }));

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

} // namespace GuessOutput
