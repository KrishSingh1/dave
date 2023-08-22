#include "daveoutputhelpers.h"
#include "guess.h"
#include "commonoutput.h"

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
	Individual *indiv = get_first_indiv(patch);
	if (indiv) {
		return indiv->fpar_layerwise.size();
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

} // namespace GuessOutput
