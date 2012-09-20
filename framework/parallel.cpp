///////////////////////////////////////////////////////////////////////////////////////
/// \file parallel.cpp
/// \brief Functionality for parallel computation
///
/// \author Joe Lindström
/// $Date$
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

namespace GuessParallel {

#ifdef HAVE_MPI

/// A class whose only purpose is to terminate the MPI library when deleted
class FinalizeCaller {
public:
	~FinalizeCaller() {
		MPI_Finalize();
	}
};

/// The auto pointer will delete the object some time after main() is finished
std::auto_ptr<FinalizeCaller> destructor;

#endif

void init(int& argc, char**& argv) {
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);

	// Make sure the MPI_Finalize function is called at program termination
	destructor = std::auto_ptr<FinalizeCaller>(new FinalizeCaller());
#endif
}

int get_rank() {
#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
#else
	return 0;
#endif
}

}
