#include "catch.hpp"

#include "commandlinearguments.h"
#include "framework.h"
#include "guess.h"
#include "config.h"
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>

static const char *const SAVE_STATE_CPOOL_OUTPUT_FILEPATH = "../tests/test_outputs/save_state/cpool.out";
static const char *const RESTART_STATE_CPOOL_OUTPUT_FILEPATH = "../tests/test_outputs/restart_state/cpool.out";
static const double TOLERANCE_kgC = 0.001;


static const char *const STATE_FILEPATH = "../tests/test_outputs/state";
static const char *const SAVE_STATE_INSFILE_FILEPATH = "../tests/test_insfiles/cfinput_save_state_test.ins";
static const char *const RESTART_STATE_INSFILE_FILEPATH = "../tests/test_insfiles/cfinput_restart_state_test.ins";

static std::string getLastLineOfFile(const char *filepath);
static void cleanUpOutput();
static double parseCpoolFromLastLine(const std::string& filepath);
static void removeFileIfExists(const std::string& filepath);

TEST_CASE("Simple state save and restart test works for one grid cell", "[integrationtest]"){

    std::cout << "Making state test directories" << std::endl;
    make_directory(STATE_FILEPATH);
    make_directory("../tests/test_outputs/save_state");
    make_directory("../tests/test_outputs/restart_state");
    
    cleanUpOutput();

    std::cout << "Starting LPJ-GUESS with cf input to save a state" << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();
    char* castedArgs[4] = { const_cast<char*>("guess"), const_cast<char*>("-input"), const_cast<char*>("cf"), const_cast<char*>(SAVE_STATE_INSFILE_FILEPATH) };
    framework(CommandLineArguments(4, castedArgs));

    double cpoolTotalAtSaveStateEnd = parseCpoolFromLastLine(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
    // std::string lastLine = getLastLineOfFile(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
    // std::istringstream iss(lastLine);
    // std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
    // printf("Last line tokens size: %zu\n", tokens.size());
    // printf("Last line: %s\n", lastLine.c_str());
    // double cpoolTotalAtSaveStateEnd = std::stof(tokens.back());
    std::cout << "Total carbon pool at end of save state simulation: " << cpoolTotalAtSaveStateEnd << " kgC" << std::endl;

    std::cout << "LPJ-GUESS save state has finished." << std::endl;
    std::cout << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();
	
	std::cout << "Starting LPJ-GUESS with cf input and restart from state" << std::endl;
    char* castedArgsRestart[4] = { const_cast<char*>("guess"), const_cast<char*>("-input"), const_cast<char*>("cf"), const_cast<char*>(RESTART_STATE_INSFILE_FILEPATH) };
    framework(CommandLineArguments(4, castedArgsRestart));
	std::cout << "LPJ-GUESS restart from state has finished." << std::endl;
		
    // check value of total carbon
    // std::string lastLine = getLastLineOfFile(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
    // std::istringstream iss(lastLine);
    // std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

	double cpoolTotalAtRestartStateEnd = parseCpoolFromLastLine(RESTART_STATE_CPOOL_OUTPUT_FILEPATH);
    std::cout << "Total carbon pool at end of restart from state simulation: " << cpoolTotalAtRestartStateEnd << " kgC" << std::endl;

    SECTION("A simple LPJ-GUESS simulation with cf input and 1 patch should produce the same data for save state and restart from state runs."){
        // If this check fails when you added code to the model, it means that there was probably no plants growing.
        // Maybe you need to adapt the insfile with which this test is run, located at ../tests/test_insfiles/cfinput_test.ins
        REQUIRE(fabs(cpoolTotalAtSaveStateEnd - cpoolTotalAtRestartStateEnd) < TOLERANCE_kgC);
    }
    std::cout << std::endl;

}

static void removeFileIfExists(const std::string& filepath) {
	std::ifstream file(filepath);
	if(file.good()){
		if(remove(filepath.c_str()) == 0){
			std::cout << "Output file " << filepath << " deleted. Ready to start the test." << std::endl;
		} else {	
			std::cout << "Output file " << filepath << " could not be deleted." << std::endl;
		}
	} else {
		std::cout << "Output file " << filepath << " does not exist. This is not an error." << std::endl;
	}
}

static void cleanUpOutput() {
	removeFileIfExists(SAVE_STATE_CPOOL_OUTPUT_FILEPATH);
	removeFileIfExists(RESTART_STATE_CPOOL_OUTPUT_FILEPATH);
    removeFileIfExists(std::string(STATE_FILEPATH) + "/meta.bin");
    removeFileIfExists(std::string(STATE_FILEPATH) + "/0.state");		
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

static double parseCpoolFromLastLine(const std::string& filepath) {
    std::string lastLine = getLastLineOfFile(const_cast<char *>(filepath.c_str()));
    std::istringstream iss(lastLine);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    double cpoolTotalAtSimulationEnd = std::stof(tokens.back());
	return cpoolTotalAtSimulationEnd;
}
