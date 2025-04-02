#include "catch.hpp"

#include "commandlinearguments.h"
#include "framework.h"
#include "guess.h"
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>

static const char *const CPOOL_OUTPUT_FILEPATH = "../tests/test_outputs/cpool.out";
static const double TOLERANCE_kgC = 5.0;

static std::string getLastLineOfFile(const char *filepath);
static void cleanUpOutput();


TEST_CASE("Simple cf input run works for one grid cell", "[integrationtest]"){

    cleanUpOutput();

    std::cout << "Starting LPJ-GUESS" << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();
    char* castedArgs[4] = { const_cast<char*>("guess"), const_cast<char*>("-input"), const_cast<char*>("cf"), const_cast<char*>("../tests/test_insfiles/cfinput_test.ins") };
    framework(CommandLineArguments(4, castedArgs));

    std::cout << "LPJ-GUESS has finished. Checking outputs..." << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();

    // check value of total carbon
    std::string lastLine = getLastLineOfFile(CPOOL_OUTPUT_FILEPATH);
    std::istringstream iss(lastLine);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    double cpoolTotalAtSimulationEnd = std::stof(tokens.back());

    SECTION("A simple LPJ-GUESS simulation with nc input and 1 patch without disturbances of stochastics should simulate some plants growing, so the carbon pool should be some positive number."){
        // If this check fails when you added code to the model, it means that there was probably no plants growing.
        // Maybe you need to adapt the insfile with which this test is run, located at ../tests/test_insfiles/cfinput_test.ins
        REQUIRE(cpoolTotalAtSimulationEnd > 0);
    }
}

static void cleanUpOutput() {
    // remove cpool.out from test outputs, to make sure that we are really writing a cpool.out file and not checking a file that is lying around there from previous runs.
    std::ifstream file(CPOOL_OUTPUT_FILEPATH);
    if(file.good()){
        if(remove(CPOOL_OUTPUT_FILEPATH) == 0){
            std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " deleted. Ready to start the test." << std::endl;
        } else {
            std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " could not be deleted." << std::endl;
        }
    } else {
        std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " does not exist. This is not an error." << std::endl;
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