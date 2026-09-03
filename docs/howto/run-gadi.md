# Running on Gadi

[Gadi](https://nci.org.au/our-systems/hpc-systems) is the HPC system at the
National Computational Infrastructure (NCI) in Canberra. This page is a
practical guide to building and submitting LPJ-GUESS DAVE jobs on Gadi. It
covers the most important concepts but is not a substitute for the
[Gadi User Guide](https://opus.nci.org.au/spaces/Help/pages/90308778/0.+Welcome+to+Gadi),
which is comprehensive and well-maintained. Reading at least the sections on
PBS, storage, and job debiting before your first run is strongly recommended.

## Connecting to Gadi

```bash
ssh <nci_username>@gadi.nci.org.au
```

Password authentication is possible but requires typing a password on every
login, which is inconvenient for repeated sessions or scripted transfers.
Setting up an SSH key pair is recommended for day-to-day use.

To do this:

1. Generate an ed25519 key pair if you don't already have one (see
   [SSH Key Management](ssh-keys.md)).
2. Register the **public key** (`id_ed25519.pub`)

Log in to Gadi with your password, then append your public key to
`~/.ssh/authorized_keys` on the server. The easiest way is `ssh-copy-id`
(Linux/macOS):

```bash
ssh-copy-id -i ~/.ssh/id_ed25519.pub user@gadi.nci.org.au
```

You will be prompted for your NCI password once. After that, SSH and SFTP
connections will use the key instead.

If `ssh-copy-id` is not available (e.g. on Windows), log in manually and
run:

```bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
cat >> ~/.ssh/authorized_keys << 'EOF'
<paste the contents of id_ed25519.pub here>
EOF
chmod 600 ~/.ssh/authorized_keys
```


Once the key is registered, `ssh <nci_username>@gadi.nci.org.au` will connect
without a password prompt.

## Background: PBS

Gadi uses **PBS Pro** (Portable Batch System) as its job scheduler. Rather than
running programs directly, you describe the resources your job needs and submit
it to PBS, which then queues and executes it when resources become available.

Key concepts:

| Concept | Description |
|---|---|
| **Job** | A self-contained unit of work submitted via `qsub` |
| **Queue** | A named pool of compute nodes (e.g. `normal`, `hugemem`, `express`) |
| **Walltime** | The maximum real-world time a job may run before PBS kills it |
| **Service Unit (SU)** | The unit of compute cost; jobs are billed to a project |
| **Project** | A budget account (e.g. `pt17`) that SUs are charged to |

Useful PBS commands:

| Command | Description |
|---|---|
| `qstat -u $USER` | List your running and queued jobs |
| `qstat -f <jobid>` | Detailed status for a specific job |
| `qdel <jobid>` | Cancel a job |

**48-hour walltime limit.** No single PBS job on Gadi may run for longer than
48 hours. For LPJ-GUESS this means that running a large domain on very few
CPUs is not feasible — even setting aside patience, the job would be killed
before it finishes. Use enough CPUs to comfortably fit within the limit (see
[`choose_num_cpus.py`](#choose_num_cpuspy-choose-cpu-count) below).

## Environment Modules

Software on Gadi is managed via the **environment modules** system. Rather than
having one global version of each tool, multiple versions coexist and users
activate the ones they need:

```bash
module load netcdf          # load the default netcdf version
module load gcc/12.1.0      # load a specific version
module list                 # list currently loaded modules
module avail netcdf         # list available versions of netcdf
module show netcdf          # show what a module does (paths it sets, etc.)
module purge                # unload all modules
```

!!! warning "Module consistency"
    The exact module versions loaded when **compiling** must also be loaded when
    **running**. Mismatching MPI implementations between build and runtime is a
    particularly common and hard-to-diagnose error. The `build_gadi` script and
    `submit_to_gadi.sh` (described below) are designed to ensure consistency —
    use them rather than loading modules manually.

### DAVE Shared Environment

DAVE helper scripts and custom module files are available to all members of the
`pt17` group on Gadi. To make them available automatically at login, add this
line to your own `~/.bashrc` file:

```bash
source /scratch/pt17/.bin/.bashrc
```

That shared file adds `/scratch/pt17/.bin` to `PATH`, making the DAVE helper
scripts available as normal shell commands. It also adds the DAVE custom module
directory to `MODULEPATH`, so `module avail` and `module load` can see the
project-specific modules. One of these modules is `lpj-guess-utils`, which
provides Gadi-compiled versions of the standard LPJ-GUESS utilities such as
`aslice`, `tslice`, `gmap`, etc.

After editing `~/.bashrc`, either log out and back in, or source it in the
current shell:

```bash
source ~/.bashrc
```

## Storage

Gadi provides three types of storage, with very different properties.

### Home Directory (`$HOME`)

`$HOME` (`/home/xxx/<nci_username>`) is limited to **10 GiB**. It is intended
for configuration files, scripts, and source code — **never for job outputs**. A
run writing even modest output to `$HOME` will quickly exhaust the quota and
cause jobs to fail.

### Scratch (`/scratch/<project>/$USER`)

`/scratch` is the primary working space for job input and output. It is fast
and has a larger quota, but files that have not been accessed or modified for
**100 days** are automatically deleted without warning.

Files are first moved to a quarantine area before permanent deletion — consult
the [Gadi User Guide](https://opus.nci.org.au/spaces/Help/pages/90308778/0.+Welcome+to+Gadi)
for how to list and recover quarantined files if needed.

!!! warning "No permanent storage"
    The DAVE project does not have permanent storage on Gadi (only scratch).
    Important results should be backed up to [HIE-Storage](hie-storage.md) as
    soon as possible after a run completes.

### Gdata (`/g/data/<project>`)

`/g/data` is persistent project storage not subject to the 100-day deletion
policy. The DAVE project does not currently have any gdata allocation.

**PBS filesystem access.** PBS jobs must explicitly declare which filesystems
they use (e.g. `#PBS -l storage=scratch/pt17+gdata/pt17`); accessing an
undeclared path fails at runtime even if you can reach it interactively. The
`submit_to_gadi.sh` script handles this automatically for standard runs —
only manually-written PBS scripts need to consider it.

## Building on Gadi

The LPJ-GUESS-DAVE repository includes a `build_gadi` script that loads the
correct modules and compiles the model. From the repository root:

```bash
./build_gadi
```

This script:

1. Runs `module purge` to clear any previously loaded modules
2. Loads `openmpi`, `gcc`, `cmake`, and `netcdf`
3. Sets `CMAKE_PREFIX_PATH` so CMake can locate NetCDF and OpenMPI
4. Runs a Release build with C++14 into the `build/` directory

## Submitting Jobs with `submit_to_gadi.sh`

For parallel production runs, use
[`submit_to_gadi.sh`](https://github.com/hie-dave/lpj-scripts/blob/master/submit_to_gadi.sh)
from the `hie-dave/lpj-scripts` repository. It automates the full submission
workflow:

1. Creates the `runs/run1/`…`runs/runN/` working directories required by the
   parallel launcher (see [Running in Parallel](run-parallel.md))
2. Splits the gridlist into N approximately equal chunks, one per rank
3. Generates an MPI execution script that loads the same modules used at build
   time and launches the job with `mpirun`
4. Generates `append.cmd`, which combines the per-worker output files
5. Unless `-d` is used, submits both jobs to PBS, with the append job scheduled
   to run only after the main job completes successfully

```bash
submit_to_gadi.sh -s my_run.conf
```

The `-s` option is required and identifies the configuration file.

### Do a Dry-run First

It is strongly recommended to prepare a run with `-d` before submitting it:

```bash
submit_to_gadi.sh -d -s my_run.conf
```

This performs the same setup as a normal invocation—including creating the
worker directories, splitting the gridlist, and generating the command files—but
does not submit either job to PBS. This gives you an opportunity to inspect the
generated files and, if needed, manually test a single gridcell on the login
node.

Keep login-node tests short and limited to a single gridcell; production runs
must still be submitted through PBS. A brief test will usually expose invalid
paths, missing instruction file parameters, and similar startup errors
immediately, instead of making you wait in the PBS queue for a job that fails as
soon as it starts.

### Configuration File

Copy [`example.conf`](https://github.com/hie-dave/lpj-scripts/blob/master/example.conf)
as a starting point. **Do not use `~`** for home directory expansion — use
`${HOME}` or an absolute path.

| Parameter | Description | Example |
|---|---|---|
| `OUT_DIR` | Root output directory. Use `/scratch`, not `$HOME` | `/scratch/pt17/some/directory` |
| `BINARY` | Path to compiled `guess` binary | `${HOME}/lpj-guess-dave/build/guess` |
| `EXPERIMENT` | Subdirectory name within `OUT_DIR` for this run | `barra2-dave` |
| `NPROCESS` | Number of MPI ranks (CPUs) | `48` |
| `WALLTIME` | Maximum job duration (`hh:mm:ss`). Must not exceed 48 hours | `02:00:00` |
| `MEMORY` | Total memory allocation. Recommend 4 GiB per CPU on `normal` queue | `192GB` |
| `QUEUE` | PBS queue | `normal` |
| `PROJECT` | NCI project to charge SUs to | `pt17` |
| `EMAIL` | Address for PBS job status notifications | `you@uni.edu.au` |
| `EMAIL_NOTIFICATIONS` | `1` to receive email, `0` to suppress | `1` |
| `INSFILE` | Path to the instruction file | `/scratch/pt17/my/simulation.ins` |
| `INPUT_MODULE` | LPJ-GUESS input module — use `nc` for DAVE | `nc` |

## After a Job Finishes

### Check the Top-level Job Log

The main job writes `<job name>.log` in the run directory. The bottom of this
file contains the PBS resource summary, including the resources requested and
used and the total job cost in service units (SU). Check this summary after
each run, especially if the job took substantially more or less time than
expected.

As a rough guide, LPJ-GUESS jobs should cost approximately **0.005–0.01 SU per
patch**, where the total number of patches is:

```text
total patches = npatch × number of gridcells
```

Therefore:

```text
lower estimate = total patches × 0.005 SU
upper estimate = total patches × 0.01 SU
```

A cost significantly outside this range should be investigated. An unusually
high cost may mean that work was duplicated between workers; an unusually low
cost may indicate that the run did not process everything expected. See
[Estimating Run Cost](#estimating-run-cost) for more detail.

### Inspect Individual Worker Logs

Each worker writes its own log at `runs/runX/guess.log`, for example
`runs/run1/guess.log`. These logs show which gridcells that worker processed,
how long they took, and any warnings or errors encountered. When a run is
unexpectedly slow, incomplete, or unusually expensive, compare several worker
logs to determine how the work was divided and where time was spent.

In particular, if every `guess.log` lists every gridcell, check the gridlist
path in the instruction file as described under
[Job takes much longer than expected](#job-takes-much-longer-than-expected).

### Combined Output Files

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

### Job takes much longer than expected

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
[DAVE Shared Environment](#dave-shared-environment), then log out and back in
or run:

```bash
source ~/.bashrc
```

### Paths containing `~` do not work in the configuration file

The submission script reads paths from the configuration file literally. Use
`${HOME}` or a full path instead of `~`:

```bash
BINARY=${HOME}/lpj-guess-dave/build/guess
```

### Job starts but cannot read or write files

First check that input and output paths exist and are on a filesystem available
to the PBS job. Use `/scratch/<project>/$USER` for run outputs, not `$HOME`.

If you wrote your own PBS script, make sure it declares the filesystems it
uses, for example:

```bash
#PBS -l storage=scratch/pt17
```

The `submit_to_gadi.sh` script handles this for standard runs.

### Job stops near the requested walltime

PBS kills jobs that exceed their walltime request. Gadi also has a maximum
48-hour walltime for a single job. Increase `NPROCESS` so the run finishes
sooner, increase `WALLTIME` if it was set too low, or split the simulation into
smaller jobs.

## Helper Scripts

Several helper scripts are available to all DAVE team members on Gadi under
`/scratch/pt17/.bin`. If you have sourced the shared startup file described in
[DAVE Shared Environment](#dave-shared-environment), they will be available on
`PATH`. Many of these scripts are maintained in the
[`hol430/gadi-scripts`](https://github.com/hol430/gadi-scripts) repository.

### `nci_account` — Check SU Balance

Prints KSU (kilo service unit) usage and availability for your projects:

```bash
nci_account
```

Run this before large jobs to confirm sufficient allocation remains.

### `lquota` — Check Storage Quotas

Prints storage usage and limits for each project you are a member of:

```bash
lquota
```

### `qcost` — Estimate Job Cost

Estimates the SU cost of a PBS job before submitting it:

```bash
qcost -c <cpus> -m <memory_GiB> -w <walltime_hours> -q <queue>
```

Example — 48 CPUs, 192 GiB, 2 hours, `normal` queue:

```bash
qcost -c 48 -m 192 -w 2 -q normal
```

### `choose_num_cpus.py` — Choose CPU Count

Helps choose an efficient number of CPUs for a given gridcell count:

```bash
choose_num_cpus.py <total_gridcells> <queue>
```

Adding more CPUs reduces walltime for the same SU cost (billing is per
CPU-hour). However, if the gridcell count does not divide evenly across ranks,
some ranks finish early and sit idle while the slowest rank completes its last
gridcell. This script finds a CPU count that minimises that imbalance as well
as execution time (walltime).

## Job Charging

Understanding how Gadi bills jobs helps avoid inadvertent overspending.

### Billing Formula

From the Gadi User Guide:

> All compute jobs are charged based on the resources **reserved** for the
> job, and the amount of walltime **used** by the job. Resources reserved are
> calculated based on the maximum of either CPUs requested, or the proportion
> of node memory requested.
>
> **Job Cost (SU) = Queue Charge Rate × max(NCPUs, Memory Proportion) × Walltime Used (hours)**
>
> - **Queue Charge Rate** — listed in the queue limits table; express queues
>   cost more per SU but increase priority.
> - **NCPUs** — number of CPUs requested via `#PBS -l ncpus`.
> - **Memory Proportion** — memory requested ÷ memory-per-core for the queue
>   (memory-per-node ÷ CPUs-per-node).

The `memory proportion` term means that requesting an outsized amount of
memory — even if you don't use it — can increase your cost: a job requesting
enough memory to fill a node effectively reserves the whole node and is charged
accordingly.

### Practical Implications for LPJ-GUESS

LPJ-GUESS uses well under 1 GiB of memory per CPU. Our costs are therefore
**always dominated by CPU count**, not memory.

The `normal` queue nodes have 48 CPUs and 192 GiB RAM each, giving a
memory-per-core ratio of 4 GiB/CPU. Because our memory usage is below this
threshold, the `Memory Proportion` term never exceeds `NCPUs`, and the billing
formula simplifies to:

```
Job Cost (SU) = Queue Charge Rate × NCPUs × Walltime Used
```

Requesting less than 4 GiB per CPU on the `normal` queue costs the same as
requesting 4 GiB per CPU. **Request 4 GiB per CPU as a safe default** — it
will not cost more and avoids unexpected out-of-memory failures.

### Estimating Run Cost

LPJ-GUESS cost on Gadi is approximately **0.005 – 0.01 SU per patch**. For a
run over G gridcells each containing P patches:

```
estimated SUs ≈ G × P × 0.005   (lower bound)
estimated SUs ≈ G × P × 0.01    (upper bound)
```

The cost per SU is approximately **$0.04 AUD** (indicative — verify against
current NCI pricing). These are first-order estimates assuming uniform patch
runtime. More productive gridcells (higher biomass turnover) take slightly
longer, but the effect is small and `choose_num_cpus.py` can minimise
imbalance across ranks.
