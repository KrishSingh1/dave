#include "catch.hpp"

#include "commandlinearguments.h"
#include "framework.h"
#include "guess.h"
#include "config.h"
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>

#include "temporary_directory.h"
#include "simulation_helper.h"

static const char *const SAVE_STATE_OUT_DIR = "out_save";
static const char *const RESTART_STATE_OUT_DIR = "out_restart";
static const double TOLERANCE_kgC = 0.001;

static const char *const SAVE_STATE_INSFILE = "cfinput_save_state_test.ins";
static const char *const RESTART_STATE_INSFILE = "cfinput_restart_state_test.ins";

static std::string getLastLineOfFile(const char *filepath);
static double parseCpoolFromLastLine(const std::string& filepath);

TEST_CASE("Simple state save and restart test works for one grid cell", "[integrationtest]"){

    TemporaryDirectory tmp;
    run_simulation(tmp, SAVE_STATE_INSFILE);

    std::string cpool = (std::filesystem::path(tmp.get_path()) / SAVE_STATE_OUT_DIR / "cpool.out").string();
    double cpoolTotalAtSaveStateEnd = parseCpoolFromLastLine(cpool);

    // std::string lastLine = getLastLineOfFile(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
    // std::istringstream iss(lastLine);
    // std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    // printf("Last line tokens size: %zu\n", tokens.size());
    // printf("Last line: %s\n", lastLine.c_str());
    // double cpoolTotalAtSaveStateEnd = std::stof(tokens.back());
    std::cout << "Total carbon pool at end of save state simulation: " << cpoolTotalAtSaveStateEnd << " kgC" << std::endl;

    std::cout << "LPJ-GUESS save state has finished." << std::endl;
    std::cout << std::endl;

	std::cout << "Starting LPJ-GUESS with cf input and restart from state" << std::endl;
    std::string ins_restart = ins_file(tmp, RESTART_STATE_INSFILE);

    run_simulation(tmp, RESTART_STATE_INSFILE);
    // check value of total carbon
    // std::string lastLine = getLastLineOfFile(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
    // std::istringstream iss(lastLine);
    // std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    cpool = (std::filesystem::path(tmp.get_path()) / RESTART_STATE_OUT_DIR / "cpool.out").string();
	double cpoolTotalAtRestartStateEnd = parseCpoolFromLastLine(cpool);
    std::cout << "Total carbon pool at end of restart from state simulation: " << cpoolTotalAtRestartStateEnd << " kgC" << std::endl;

    SECTION("A simple LPJ-GUESS simulation with cf input and 1 patch should produce the same data for save state and restart from state runs."){
        // If this check fails when you added code to the model, it means that
        // there was probably no plants growing. Maybe you need to adapt the
        // insfile with which this test is run, located at
        // ../tests/test_insfiles/cfinput_test.ins
        REQUIRE(fabs(cpoolTotalAtSaveStateEnd - cpoolTotalAtRestartStateEnd) < TOLERANCE_kgC);
    }
    std::cout << std::endl;
}

static std::string getLastLineOfFile(const char *filepath) {
    std::ifstream file(filepath);
    std::string lastLine;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + std::string(filepath));
    }

    std::string line;
    while (std::getline(file, line)) {
        lastLine = line;
    }
    file.close();

    return lastLine;
}

static double parseCpoolFromLastLine(const std::string& filepath) {
    std::string lastLine = getLastLineOfFile(const_cast<char *>(filepath.c_str()));
    std::istringstream iss(lastLine);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    if (tokens.empty()) {
        throw std::runtime_error("Failed to parse cpool from last line of file: " + filepath);
    }

    double cpoolTotalAtSimulationEnd = std::stof(tokens.back());
	return cpoolTotalAtSimulationEnd;
}
