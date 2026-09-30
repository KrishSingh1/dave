# Checking and Troubleshooting Gadi Runs

This page explains how to check a completed LPJ-GUESS DAVE job on Gadi and
diagnose common failures. For building and submitting a job, see
[Running on Gadi](run-gadi.md).

## Start with the Job and Worker Logs

The main job writes `<job name>.log` in the run directory. The bottom of this
file contains the PBS resource summary, including the resources requested and
used and the total job cost in service units (SU). Check this summary after
each run, especially if the job took substantially more or less time than
expected.

For a model error or an unexpected scientific result, open the affected
worker's `runs/runX/guess.log` first. The top-level job log is especially
useful for PBS startup, resource, and walltime problems. The instruction
file's `verbosity` controls each worker's model log: use at least `3` for
informational messages, or `4`/`5` for more detail, and rerun if needed. See
[Start with `guess.log`](common-failures.md#start-with-guesslog) for a search
strategy.

As a rough guide, LPJ-GUESS jobs should cost approximately **0.005–0.01 SU per
patch**, where:

```text
total patches = npatch × number of gridcells

lower cost estimate = total patches × 0.005 SU
upper cost estimate = total patches × 0.01 SU
```

A cost significantly outside this range should be investigated. An unusually
high cost may mean that work was duplicated between workers; an unusually low
cost may indicate that the run did not process everything expected. See
[Estimating LPJ-GUESS Run Cost](gadi-resources.md#estimating-lpj-guess-run-cost)
for more detail.

## Inspect Individual Worker Logs

Each worker writes its own log at `runs/runX/guess.log`, for example
`runs/run1/guess.log`. These logs show which gridcells that worker processed,
how long they took, and any warnings or errors encountered. When a run is
unexpectedly slow, incomplete, or unusually expensive, compare several worker
logs to determine how the work was divided and where time was spent.

In particular, if every `guess.log` lists every gridcell, check for the
[absolute gridlist path problem](#job-takes-much-longer-than-expected).

## Final Output Files and `append.cmd`

The submission script creates a second PBS job from `append.cmd`. PBS schedules
this job to run after the main job, but only if the main job succeeds. It
combines each worker's output files into the final `output/` directory. For
example:

```text
runs/run1/out/lai.out ┐
runs/run2/out/lai.out ├──> output/lai.out
runs/run3/out/lai.out ┘
```

The same process is applied to every output file. If the main job fails, the
append job does not run and the individual files remain under
`runs/runX/out/` for diagnosis or manual recovery.

## Common Errors

### Job Takes Much Longer Than Expected

One possible cause is that every worker has processed the complete gridlist
instead of only its assigned gridcells. The run may otherwise appear to work,
so this is easy to miss without checking the job cost and worker logs.

!!! warning "The gridlist path must be relative"
    When using `submit_to_gadi.sh`, never give the gridlist an absolute path in
    the instruction file. Use `./<gridlist-basename>` instead. This is an
    exception to the general recommendation to use absolute paths for input
    files.

For example, if the gridlist is named `gridlist.txt`:

```text
# Wrong: every worker reads the complete, original gridlist
file_gridlist "/scratch/pt17/my/input/gridlist.txt"

# Correct: every worker reads its assigned part of the gridlist
file_gridlist "./gridlist.txt"
```

The parameter name depends on the input module; for example, it may be
`file_gridlist` or `file_gridlist_cf`. The same rule applies in either case,
and the relative filename must match the basename of the original gridlist.

The submission script splits the gridlist between the
`runs/run1/`...`runs/runN/` directories. During a parallel run, each worker
changes into its own directory and should read the split gridlist stored there.
An absolute path bypasses that file, causing every worker to process the
complete original gridlist. This can multiply the job's cost and cause it to
exceed its requested walltime. See
[Path Considerations](run-parallel.md#path-considerations) for more detail.

### `libnetcdf.so.15: cannot open shared object file`

Example:

```text
/home/590/am8995/code/lpjg_aaron/build/guess: error while loading shared libraries: libnetcdf.so.15: cannot open shared object file: No such file or directory
```

This means the executable was compiled against a NetCDF library that is not
available in the environment where the model is being run. The most common
cause is forgetting to load the `netcdf` module before running. It can also
happen if the model was compiled with one NetCDF module version but run with a
different one.

Fix:

1. Use `./build_gadi` to compile the model. It loads the required Gadi modules,
   including the default `netcdf` module.
2. Submit runs with `submit_to_gadi.sh`, which generates a PBS script that
   loads the same kind of runtime environment before launching `guess`.
3. If running manually for a quick test, load NetCDF first:

   ```bash
   module purge
   module load openmpi gcc cmake netcdf
   ./build/guess -input nc /path/to/simulation.ins
   ```

If the error persists, check which modules are loaded in the job output:

```bash
module list
```

Then rebuild and rerun in a clean, consistent environment:

```bash
module purge
./build_gadi
```

### `submit_to_gadi.sh: command not found`

The DAVE helper scripts are probably not on your `PATH`. Add the shared DAVE
startup file to `~/.bashrc` as described in
[DAVE Shared Environment](gadi-setup.md#dave-shared-environment), then log out
and back in or run:

```bash
source ~/.bashrc
```

### Paths Containing `~` Do Not Work in the Configuration File

The submission script reads paths from the configuration file literally. Use
`${HOME}` or a full path instead of `~`:

```bash
BINARY=${HOME}/lpj-guess-dave/build/guess
```

### Job Starts but Cannot Read or Write Files

First check that input and output paths exist and are on a filesystem available
to the PBS job. Use `/scratch/<project>/$USER` for run outputs, not `$HOME`.

If you wrote your own PBS script, make sure it declares the filesystems it
uses, for example:

```bash
#PBS -l storage=scratch/pt17
```

The `submit_to_gadi.sh` script handles this for standard runs. See
[Storage](gadi-setup.md#storage) for more detail.

### Job Stops Near the Requested Walltime

PBS kills jobs that exceed their walltime request. Gadi also has a maximum
48-hour walltime for a single job. Increase `NPROCESS` so the run finishes
sooner, increase `WALLTIME` if it was set too low, or split the simulation into
smaller jobs. See
[Choosing CPU Count and Walltime](gadi-resources.md#choosing-cpu-count-and-walltime).
