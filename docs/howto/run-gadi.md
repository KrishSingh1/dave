# Running on Gadi

This page is the practical workflow for building and submitting an LPJ-GUESS
DAVE job on the [Gadi HPC system](https://nci.org.au/our-systems/hpc-systems).
If this is your first time using Gadi, complete
[Setting Up Gadi](gadi-setup.md) before continuing.

Related guides provide more detail when you need it:

- [Setting Up Gadi](gadi-setup.md)—accounts, SSH, PBS, modules, and storage
- [Checking and Troubleshooting Gadi Runs](gadi-troubleshooting.md)—logs,
  output collation, and common errors
- [Choosing Gadi Resources and Estimating Costs](gadi-resources.md)—CPU,
  memory, walltime, and SU estimates

## Before You Start

You should have:

- access to Gadi and an NCI project such as `pt17`
- the LPJ-GUESS-DAVE repository and your input files on Gadi
- the [DAVE shared environment](gadi-setup.md#dave-shared-environment) enabled
- a run output location under `/scratch`, not `$HOME`

Connect to Gadi:

```bash
ssh <nci_username>@gadi.nci.org.au
```

## Build LPJ-GUESS

The repository includes a `build_gadi` script that loads the correct modules
and compiles the model. From the repository root:

```bash
./build_gadi
```

The script clears previously loaded modules, loads OpenMPI, GCC, CMake, and
NetCDF, and creates a Release build in `build/`. Use this script to avoid
mismatches between the modules used to build and run LPJ-GUESS.

## Create a Submission Configuration

Copy [`example.conf`](https://github.com/hie-dave/lpj-scripts/blob/master/example.conf)
from the `hie-dave/lpj-scripts` repository as a starting point.

| Parameter | Description | Example |
|---|---|---|
| `OUT_DIR` | Root output directory. Use `/scratch`, not `$HOME` | `/scratch/pt17/some/directory` |
| `BINARY` | Path to the compiled `guess` binary | `${HOME}/lpj-guess-dave/build/guess` |
| `EXPERIMENT` | Subdirectory name within `OUT_DIR` for this run | `barra2-dave` |
| `NPROCESS` | Number of parallel workers (CPUs) | `48` |
| `WALLTIME` | Maximum job duration (`hh:mm:ss`); no more than 48 hours | `02:00:00` |
| `MEMORY` | Total memory allocation; use 4 GiB per CPU on `normal` | `192GB` |
| `QUEUE` | PBS queue | `normal` |
| `PROJECT` | NCI project to charge SUs to | `pt17` |
| `EMAIL` | Address for PBS job-status notifications | `you@uni.edu.au` |
| `EMAIL_NOTIFICATIONS` | `1` to receive email, `0` to suppress | `1` |
| `INSFILE` | Path to the instruction file | `/scratch/pt17/my/simulation.ins` |
| `INPUT_MODULE` | LPJ-GUESS input module; use `nc` for DAVE | `nc` |

Do not use `~` for home-directory expansion in this file. Use `${HOME}` or an
absolute path instead. See
[Choosing Gadi Resources and Estimating Costs](gadi-resources.md) before
selecting `NPROCESS`, `MEMORY`, and `WALLTIME` for a large run.

### Use a Relative Gridlist Path

!!! warning "The gridlist is the exception to the absolute-path rule"
    In the instruction file, the gridlist must use a relative path such as
    `./gridlist.txt`. An absolute gridlist path can cause every worker to run
    every gridcell, making the job take much longer and cost far more than
    expected.

The parameter may be named `file_gridlist`, `file_gridlist_cf`, or something
else depending on the input module. Its value must be
`./<gridlist-basename>`—for example:

```text
file_gridlist "./gridlist.txt"
```

All other input files should generally use absolute paths. See
[Job Takes Much Longer Than Expected](gadi-troubleshooting.md#job-takes-much-longer-than-expected)
for an explanation of this special case.

## Prepare a Dry Run

For parallel production runs, use `submit_to_gadi.sh` from the
[`hie-dave/lpj-scripts`](https://github.com/hie-dave/lpj-scripts/blob/master/submit_to_gadi.sh)
repository. Start with the `-d` option:

```bash
submit_to_gadi.sh -d -s my_run.conf
```

The required `-s` option identifies the configuration file. The `-d` option
performs all preparation but does not submit jobs to PBS. The script:

1. Creates the `runs/run1/`…`runs/runN/` working directories required by the
   parallel launcher (see [Running in Parallel](run-parallel.md)).
2. Splits the gridlist into approximately equal chunks, one per worker.
3. Generates the MPI execution command file, which loads the build-time
   modules and launches LPJ-GUESS with `mpirun`.
4. Generates `append.cmd`, which will combine the workers' output files.

Inspect the generated files and, if needed, manually test a single gridcell on
the login node. Keep login-node tests short and limited to one gridcell;
production runs must use PBS. A brief test will usually expose invalid paths,
missing instruction-file parameters, and similar startup errors immediately.

## Submit the Jobs

Once the dry run and any single-gridcell checks succeed, submit normally:

```bash
submit_to_gadi.sh -s my_run.conf
```

The script submits two PBS jobs:

1. The main job runs LPJ-GUESS in parallel.
2. The append job is scheduled to run only if the main job succeeds. It
   combines the workers' individual output files into `output/`.

You can monitor or cancel the main job with standard PBS commands:

```bash
qstat -u $USER
qstat -f <jobid>
qdel <jobid>
```

## Check the Result

After the jobs finish, check these files:

| Path | Purpose |
|---|---|
| `<job name>.log` | Top-level job log; its final lines report requested and used PBS resources and total SU cost |
| `runs/runX/guess.log` | Per-worker logs showing which gridcells ran, their timing, warnings, and errors |
| `runs/runX/out/` | Individual worker outputs |
| `output/` | Combined outputs created by `append.cmd` after a successful main job |

For example, the append job combines:

```text
runs/run1/out/lai.out ┐
runs/run2/out/lai.out ├──> output/lai.out
runs/run3/out/lai.out ┘
```

Check the worker `guess.log` files first when the model fails or produces
unexpected results. `verbosity` in the `.ins` file controls their detail;
use at least `3` for investigation, or `4`/`5` for more detail. Check the
top-level job log for PBS errors and resource usage after every run.
LPJ-GUESS jobs typically cost
around **0.005–0.01 SU per patch**, where the total number of patches is
`npatch × number of gridcells`. A substantially higher or lower cost should be
investigated.

See [Checking and Troubleshooting Gadi Runs](gadi-troubleshooting.md) for a
diagnostic workflow and [Estimating LPJ-GUESS Run Cost](gadi-resources.md#estimating-lpj-guess-run-cost)
for the full cost calculation.
