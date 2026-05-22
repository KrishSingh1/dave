# Running in Parallel (MPI)

LPJ-GUESS can distribute a simulation across multiple CPU cores using MPI.
This page covers general parallel execution. For running on the Gadi HPC
system specifically, see [Running on Gadi](run-gadi.md).

## Prerequisites

The binary must have been compiled with MPI support. CMake detects MPI
automatically if an MPI implementation (e.g. OpenMPI, MPICH) is installed.
See [Dependencies and Toolchains](../development/toolchains.md) for details.

## Parallelisation Strategy

Parallelism operates at the **gridcell level**. Each MPI rank is assigned a
disjoint subset of gridcells from a separate gridlist file and processes them
independently. There is no communication between ranks during the simulation.

Consequences:

- Performance scales linearly with the number of processes (tested to tens of
  thousands of CPUs on NCI Gadi).
- The number of processes must not exceed the number of gridcells.
- Each rank writes its own output files, which must be merged after the run
  (see [Collating Outputs](#collating-outputs) below).

## Before Running: Directory Setup

Each rank expects a subdirectory named `run<N>` (1-indexed) relative to the
working directory where `mpirun` is invoked:

- Rank 0 → `./run1/`
- Rank 1 → `./run2/`
- ...
- Rank N−1 → `./runN/`

Each `runN/` directory must contain a gridlist file with the subset of
gridcells assigned to that rank. A typical layout for 3 ranks:

```
run_dir/
├── run1/
│   └── gridlist.txt
├── run2/
│   └── gridlist.txt
└── run3/
    └── gridlist.txt
```

Split the master gridlist into N roughly equal parts — one per directory. On
Gadi, the `submit_to_gadi.sh` script handles this automatically (see
[Running on Gadi](run-gadi.md)).

## Path Considerations

When `-parallel` is specified, each rank changes its working directory to its
`runN/` directory before the model starts. Paths in the instruction file are
then resolved relative to that directory.

**Use absolute paths for all parameters except the gridlist.** The gridlist
parameter must be relative (e.g. `./gridlist.txt`) because each rank resolves it
inside its own `runN/` directory. Additionally, the output directory parameter
must also be relative, so that each worker writes to its own output directory.
All other paths (forcing files, soil data, etc.) should be absolute to avoid
hard-to-diagnose path resolution failures.

## Running

From the prepared working directory:

```bash
mpirun -np <N> /absolute/path/to/build/guess -parallel -input nc /absolute/path/to/simulation.ins
```

The `-parallel` flag is required to activate the MPI code paths. Omitting it
while using `mpirun` will cause each rank to run the full simulation
independently rather than cooperating.

## Collating Outputs

Each rank writes output files (tab-separated values, `.out`) into its `runN/`
directory. After the run, these must be merged. Since all files share the same
header row, keep the header from the first rank and concatenate data rows from
all others:

```bash
# Example: merge aaet.out from all ranks into a single file
head -1 run1/aaet.out > aaet_merged.out
for dir in run*/; do
    tail -n +2 "${dir}aaet.out"
done >> aaet_merged.out
```

Repeat for each output variable. On Gadi, the `submit_to_gadi.sh` script
generates and submits a merge job automatically.

## Output During Parallel Runs

By default, progress messages from all ranks are interleaved on stdout. To
force each progress line onto its own line, add to the instruction file:

```
progress_newline 1
```
