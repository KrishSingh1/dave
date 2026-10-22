#include "catch.hpp"

#include "commandlinearguments.h"
#include "framework.h"
#include "guess.h"
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>

#include "temporary_directory.h"
#include "simulation_helper.h"
#include <filesystem>

static const char *const CPOOL_OUTPUT_FILEPATH = "cpool.out";
static const double TOLERANCE_kgC = 5.0;

static std::string getLastLineOfFile(const char *filepath);

TEST_CASE("Simple cf input run works for one grid cell", "[integrationtest]"){

    TemporaryDirectory tmp;
    run_simulation(tmp, "cfinput_test.ins");

    // check value of total carbon
    std::string out_dir = "test_outputs";
    std::string cpool = (std::filesystem::path(tmp.get_path()) / out_dir / CPOOL_OUTPUT_FILEPATH).string();
    std::string lastLine = getLastLineOfFile(cpool.c_str());
    std::istringstream iss(lastLine);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    double cpoolTotalAtSimulationEnd = std::stof(tokens.back());

    SECTION("A simple LPJ-GUESS simulation with nc input and 1 patch without disturbances of stochastics should simulate some plants growing, so the carbon pool should be some positive number."){
        // If this check fails when you added code to the model, it means that there was probably no plants growing.
        // Maybe you need to adapt the insfile with which this test is run, located at ../tests/test_insfiles/cfinput_test.ins
        REQUIRE(cpoolTotalAtSimulationEnd > 0);
    }
}

static std::string getLastLineOfFile(const char *filepath) {
    std::ifstream file(filepath);
    std::string lastLine;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            lastLine = line;
        }
        file.close();
    }
    return lastLine;
}
