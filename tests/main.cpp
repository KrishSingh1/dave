///////////////////////////////////////////////////////////////////////////////////////
/// \file main.cpp
/// \brief Main function for the unit tests
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "shell.h"

int main(int argc, char** argv) {

	// Set a shell so we have working dprintf, fail etc.
	set_shell(new CommandLineShell("tests.log"));

	// Let CATCH do the rest
	int result = Catch::Main(argc, argv);

	return result;
}
