///////////////////////////////////////////////////////////////////////////////////////
/// \file parallel.cpp
/// \brief Functionality for parallel computation
///
/// \author Joe Siltberg
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

// mpi.h needs to be the first include (or specifically, before stdio.h),
// due to a bug in the MPI-2 standard (both stdio.h and the C++ MPI API
// defines SEEK_SET, SEEK_CUR and SEEK_SET(!))
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "config.h"
#include "parallel.h"
#include "shell.h"
#include <memory>
#include <string>

namespace GuessParallel {

bool parallel = false;
static bool barrier_sync_hazard = false;

void init(int& argc, char**& argv) {
#ifdef HAVE_MPI

	// Only initialize MPI is the command line options explicitly
	// asks for a parallel run (with the -parallel option).
	// In most cases it wouldn't hurt to initialize MPI even if
	// it's not used, but apparently some implementations start
	// by going to the users home directory if the program isn't
	// launched by an mpiexec/mpirun launcher.

	// Unfortunately, since MPI initialization must be done before
	// we parse our options with CommandLineArguments, we need to
	// look for the -parallel option here by ourselves.
	// The file-global variable parallel is initiated = false;
	for (int i = 0; i < argc; ++i) {
		if (std::string(argv[i]) == "-parallel") {
			parallel = true;
			break;
		}
	}

	if (parallel) {
		MPI_Init(&argc, &argv);
	}
#endif
}

int get_rank() {
#ifdef HAVE_MPI
	if (parallel) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}
	else { 
		return 0;
	}
#else
	return 0;
#endif
}

int get_num_processes() {
#ifdef HAVE_MPI
	if (parallel) {
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		return size;
	}
	else {
		return 1;	
	}
#else
	return 1;
#endif
}

void finalize() {
#ifdef HAVE_MPI
	if (parallel) {
		MPI_Finalize();
	}
#endif
}

void begin_barrier_sync_section() {
	barrier_sync_hazard = true;
}

void end_barrier_sync_section() {
	barrier_sync_hazard = false;
}

void barrier() {
#ifdef HAVE_MPI
	if (parallel) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
}

bool barrier_sync_hazard_active() {
	return barrier_sync_hazard;
}

void abort_parallel(int errorcode) {
#ifdef HAVE_MPI
	if (parallel) {
		MPI_Abort(MPI_COMM_WORLD, errorcode);
	}
#endif
	// If we get here (no MPI or MPI_Abort didn't terminate), fall back.
	exit(errorcode);
}

}
