///////////////////////////////////////////////////////////////////////////////////////
/// \file parallel.h
/// \brief Functionality for parallel computation
///
///        Currently a layer over MPI with functionality LPJ-GUESS might need,
///        but could be implemented as a layer over something else in future.
///
///        The fact that MPI might not be available on the system where LPJ-GUESS
///        was compiled is hidden by this module by simply pretending to be a
///        "parallel" run with just one process.
///
/// \author Joe Siltberg
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_PARALLEL_H
#define LPJ_GUESS_PARALLEL_H

namespace GuessParallel {

/// Initializes the module
/** For instance by initializing underlying communications library.
 *
 *  \param argc Passed along from main() to underlying communications library
 *  \param argv Passed along from main() to underlying communications library
 */
void init(int& argc, char**& argv);

/// Terminate the MPI library
void finalize();

/// The process id of this process in the parallel run
/** Returns zero when no MPI library is available/used. */
int get_rank();

/// The number of processes involved in the parallel run
/** Returns 1 when no MPI library is available/used. */
int get_num_processes();

/// Mark the beginning of a section where some ranks may be in an MPI_Barrier
/// while others are still doing work.
void begin_barrier_sync_section();

/// Mark the end of a section started with begin_barrier_sync_section().
void end_barrier_sync_section();

/// Wrapper for MPI_Barrier used together with begin_barrier_sync_section() /
/// end_barrier_sync_section().
void barrier();

/// Returns true if a barrier-sync section is currently active on this rank.
bool barrier_sync_hazard_active();

/// Abort the parallel run immediately (e.g. via MPI_Abort when available).
void abort_parallel(int errorcode);

/*
 * Example usage for the "rank 0 does extra work while other ranks wait in a
 * barrier" pattern that avoids deadlocks on fail():
 *
 *   GuessParallel::begin_barrier_sync_section();
 *   int rank = GuessParallel::get_rank();
 *   if (rank != 0) {
 *       GuessParallel::barrier();  // ranks 1..N wait here
 *   } else {
 *       rank0_only_work();
 *       GuessParallel::barrier();  // rank 0 joins the same barrier
 *   }
 *   GuessParallel::end_barrier_sync_section();
 *
 * Any fail() occurring between begin_barrier_sync_section() and
 * end_barrier_sync_section() will then abort the MPI job instead of trying to
 * finalize cleanly, preventing deadlock if some ranks are inside the barrier.
 */

}

#endif // LPJ_GUESS_PARALLEL_H
