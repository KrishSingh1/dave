#include "catch.hpp"
#include "guess.h"
#include <filesystem>
#include "simulation_helper.h"
#include "temporary_directory.h"

#include <fstream>
#include <map>
#include <sstream>
#include <tuple>

namespace {

template<typename T>
void roundtrip(T& source, T& restored) {
    std::stringstream bytes(std::ios::in | std::ios::out | std::ios::binary);
    ArchiveOutStream output(bytes);
    output & source;
    ArchiveInStream input(bytes);
    input & restored;
    REQUIRE(bytes.good());
}

using OutputKey = std::tuple<double, double, int>;
using OutputRows = std::map<OutputKey, std::vector<double>>;

OutputRows read_annual_output(const std::string& path) {
    std::ifstream file(path);
    REQUIRE(file.is_open());
    std::string line;
    REQUIRE(static_cast<bool>(std::getline(file, line))); // column headings
    OutputRows rows;
    while (std::getline(file, line)) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        std::istringstream record(line);
        double lon, lat;
        int year;
        REQUIRE(static_cast<bool>(record >> lon >> lat >> year));
        std::vector<double> values;
        double value;
        while (record >> value) {
            REQUIRE(std::isfinite(value));
            values.push_back(value);
        }
        REQUIRE(record.eof());
        REQUIRE_FALSE(values.empty());
        REQUIRE(rows.emplace(OutputKey(lon, lat, year), values).second);
    }
    REQUIRE_FALSE(rows.empty());
    return rows;
}

}

TEST_CASE("Dynamic histories retain their next-day behavior after restart", "[serialization][historic]") {
    HistoricDynamic<double> source;
    HistoricDynamic<double> restored;
    restored.set_capacity(2);
    restored.add(99.0);

    SECTION("empty history") {
        roundtrip(source, restored);
        REQUIRE(restored.get_capacity() == 0);
        REQUIRE(restored.size() == 0);
        restored.set_capacity(3);
        restored.add(0.5);
        REQUIRE(restored.mean() == 0.5);
    }

    SECTION("partially filled history") {
        source.set_capacity(5);
        source.add(0.1);
        source.add(0.8);
        roundtrip(source, restored);
        REQUIRE(restored.get_capacity() == 5);
        REQUIRE(restored.size() == 2);
        for (int day = 0; day < 8; ++day) {
            source.add(day * 0.125);
            restored.add(day * 0.125);
            REQUIRE(restored.mean() == source.mean());
            REQUIRE(restored.weighted_mean() == source.weighted_mean());
        }
    }

    SECTION("wrapped history") {
        source.set_capacity(3);
        for (int day = 0; day < 5; ++day) source.add(day * 0.125);
        roundtrip(source, restored);
        REQUIRE(restored.get_capacity() == 3);
        REQUIRE(restored.size() == 3);
        for (int day = 5; day < 12; ++day) {
            source.add(day * 0.125);
            restored.add(day * 0.125);
            REQUIRE(restored.mean() == source.mean());
            REQUIRE(restored.weighted_mean() == source.weighted_mean());
        }
    }
}

TEST_CASE("Restart preserves plant and patch carry-over state", "[serialization]") {
    date.init(1);
    pftlist.killall();
    npft = 1;
    Pft& pft = pftlist.createobj();
    pft.lifeform = TREE;
    pft.landcover = NATURAL;
    pft.kmax = 1.0;
    pft.g1 = 2.0;
    pft.root_beta = 0.95;
    pft.cton_leaf_min = 20.0;
    pft.cton_leaf_max = 40.0;
    std::fill_n(pft.input_rootdist, NSOILLAYERS_MAX, 0.0);
    pft.input_rootdist[0] = 0.25;
    pft.input_rootdist[1] = 0.75;

    Soiltype soiltype;
    Gridcell gridcell;
    Stand stand(0, &gridcell, soiltype, NATURAL, 1);
    Patch& patch = stand.getobj();
    Patchpft& source_pft = patch.pft[0];
    Patchpft restored_pft(0, pft);
    source_pft.rootdist[0] = 0.3;
    source_pft.rootdist[1] = 0.7;
    source_pft.wcont_germ_4.add(0.25);
    source_pft.wcont_germ_4.add(0.75);
    source_pft.tsoil_germ_4.add(12.0);
    source_pft.wscal_historical.add(0.25);
    roundtrip(source_pft, restored_pft);
    REQUIRE(restored_pft.rootdist[0] == 0.3);
    REQUIRE(restored_pft.rootdist[1] == 0.7);
    source_pft.wcont_germ_4.add(0.5);
    restored_pft.wcont_germ_4.add(0.5);
    REQUIRE(restored_pft.wcont_germ_4.mean() == source_pft.wcont_germ_4.mean());
    REQUIRE(restored_pft.tsoil_germ_4.mean() == 12.0);
    REQUIRE(restored_pft.wscal_historical.mean() == source_pft.wscal_historical.mean());

    Individual& source = patch.vegetation.createobj(pft, patch.vegetation);
    Individual restored(1, pft, patch.vegetation);
    source.growth_stage_days = 23;
    source.fdev_growth = 0.4;
    source.cmass_crown_flowering = 0.125;
    source.meanleafage = 1.5;
    source.meanrootage = 2.5;
    source.meansapage = 3.5;
    source.c_d_diam = 0.015;
    source.cmass_repr_10.add(0.25);
    source.root_beta = 0.97;
    source.rootdist[0] = 0.4;
    source.rootdist[1] = 0.6;
    source.set_wscal_dynamic_capacity(3);
    source.soiltemp_dynamic.set_capacity(3);
    source.set_wscal(0.25);
    source.set_wscal(0.75);
    source.soiltemp_dynamic.add(12.0);
    source.alpha_leaf_trend.add(0.3);
    source.kmax_drought_affected = 0.4;
    source.stomata_open = false;
    source.lwp.back() = -3.0;
    source.lwp_previous = -3.0;
    source.hydraulic_failure = true;
    source.sapwood_depth = 0.02;

    roundtrip(source, restored);
    REQUIRE(restored.growth_stage_days == 23);
    REQUIRE(restored.fdev_growth == 0.4);
    REQUIRE(restored.cmass_crown_flowering == 0.125);
    REQUIRE(restored.meanleafage == 1.5);
    REQUIRE(restored.meanrootage == 2.5);
    REQUIRE(restored.meansapage == 3.5);
    REQUIRE(restored.c_d_diam == 0.015);
    REQUIRE(restored.cmass_repr_10.mean() == 0.25);
    REQUIRE(restored.root_beta == 0.97);
    REQUIRE(restored.rootdist[0] == 0.4);
    REQUIRE(restored.rootdist[1] == 0.6);
    REQUIRE(restored.soiltemp_dynamic.mean() == 12.0);
    REQUIRE(restored.get_wscal_dynamic_capacity() == 3);
    source.set_wscal(0.5);
    restored.set_wscal(0.5);
    REQUIRE(restored.get_wscal_moving_average() == source.get_wscal_moving_average());
    REQUIRE(restored.alpha_leaf_trend.mean() == 0.3);
    REQUIRE(restored.kmax_drought_affected == 0.4);
    REQUIRE_FALSE(restored.stomata_open);
    REQUIRE(restored.lwp.back() == -3.0);
    REQUIRE(restored.lwp_previous == -3.0);
    REQUIRE(restored.hydraulic_failure);
    REQUIRE(restored.sapwood_depth == 0.02);
}

TEST_CASE("Restart preserves hydrology initialization and soil state", "[serialization]") {
    date.init(1);
    npft = 0;
    pftlist.killall();
    Gridcell gridcell;
    Soiltype soiltype;
    Stand stand(0, &gridcell, soiltype, NATURAL, 1);
    Soil& source = stand.getobj().soil;
    Soil restored(stand.getobj(), soiltype);
    source.nsoillayers = 1;
    source.soillayers.emplace_back(&source, 100.0, 50.0, false, false);
    source.soillayers[0].wmass = 25.0;
    source.soillayers[0].temp_old = 7.0;
    source.soillayers[0].sompool[0].cmass = 0.125;
    source.soillayers[0].sompool[0].cmass_init = 0.1;
    source.soillayers[0].sompool[0].multi = 8.0;
    source.soillayers[0].sompool[0].sign_trend = -1;
    source.soillayers[0].sompool[0].same_trend = 2;
    source.soillayers[0].sompool[0].changed_trend = 1;
    source.ssw.wmass_frac_thresh = 0.4;
    source.ssw.wmass_thresh = 12.5;
    roundtrip(source, restored);
    REQUIRE(restored.nsoillayers == 1);
    REQUIRE(restored.soillayers[0].soil == &restored);
    REQUIRE(restored.soillayers[0].wmass == 25.0);
    REQUIRE(restored.soillayers[0].temp_old == 7.0);
    REQUIRE(restored.soillayers[0].sompool[0].cmass == 0.125);
    REQUIRE(restored.soillayers[0].sompool[0].cmass_init == 0.1);
    REQUIRE(restored.soillayers[0].sompool[0].multi == 8.0);
    REQUIRE(restored.soillayers[0].sompool[0].sign_trend == -1);
    REQUIRE(restored.soillayers[0].sompool[0].same_trend == 2);
    REQUIRE(restored.soillayers[0].sompool[0].changed_trend == 1);
    REQUIRE(restored.ssw.wmass_frac_thresh == 0.4);
    REQUIRE(restored.ssw.wmass_thresh == 12.5);

    stand.days_grazed = 4;
    stand.patches_grazed = {0};
    Stand restored_stand(1, &gridcell, soiltype, NATURAL, 1);
    roundtrip(stand, restored_stand);
    REQUIRE(restored_stand.days_grazed == 4);
    REQUIRE(restored_stand.patches_grazed == stand.patches_grazed);
}

TEST_CASE("Restart reproduces the uninterrupted annual trajectory", "[serialization][integrationtest]") {
    TemporaryDirectory tmp;
    run_simulation(tmp, "cfinput_save_state_test.ins");
    run_simulation(tmp, "cfinput_restart_state_test.ins");

    // Compare all pools/PFTs and fluxes in the common years, rather than
    // allowing compensating errors in a single final total carbon value.
    for (const char* filename : {"cpool.out", "npool.out", "cmass.out", "nmass.out",
                                  "anpp.out", "agpp.out", "aaet.out", "dens.out"}) {
        INFO(filename);
        OutputRows uninterrupted = read_annual_output(out_file_path(tmp, "out_save", filename));
        OutputRows restarted = read_annual_output(out_file_path(tmp, "out_restart", filename));
        REQUIRE(restarted.size() < uninterrupted.size());
        REQUIRE(restarted.rbegin()->first == uninterrupted.rbegin()->first);
        for (const auto& row : restarted) {
            INFO("Year " << std::get<2>(row.first));
            auto expected = uninterrupted.find(row.first);
            REQUIRE(expected != uninterrupted.end());
            REQUIRE(row.second == expected->second);
        }
        size_t common_years = 0;
        for (const auto& row : uninterrupted) {
            if (row.first >= restarted.begin()->first) ++common_years;
        }
        REQUIRE(restarted.size() == common_years);
    }
}
