# Choosing Gadi Resources and Estimating Costs

This page explains how to choose CPU, memory, and walltime requests for an
LPJ-GUESS DAVE job on Gadi and how to estimate its cost. For the submission
workflow, see [Running on Gadi](run-gadi.md).

## Choosing CPU Count and Walltime

No single PBS job on Gadi may run for longer than **48 hours**. Running a large
domain on very few CPUs is therefore not feasible: the job will be killed if it
reaches its walltime limit.

The shared DAVE environment provides `choose_num_cpus.py`, which helps choose
an efficient CPU count for a given gridcell count and queue:

```bash
choose_num_cpus.py <total_gridcells> <queue>
```

Adding more CPUs reduces elapsed time for the same approximate SU cost, since
billing is per CPU-hour. However, if the gridcell count does not divide evenly
across workers, some finish early and sit idle while the slowest worker
completes its last gridcell. The script finds a CPU count that minimises that
imbalance as well as execution time.

Choose a `WALLTIME` with enough margin for the expected run duration. PBS kills
jobs that exceed their requested walltime even when the queue's 48-hour maximum
has not been reached.

## Check the Available Allocation

Run `nci_account` before large jobs to confirm that the project has sufficient
allocation remaining. It prints KSU (kilo service unit) usage and availability
for your projects:

```bash
nci_account
```

These helper scripts are available under `/scratch/pt17/.bin`. If the commands
are not found, set up the
[DAVE shared environment](gadi-setup.md#dave-shared-environment).

## Estimate a Job Before Submission

The `qcost` helper estimates the SU cost of a PBS job:

```bash
qcost -c <cpus> -m <memory_GiB> -w <walltime_hours> -q <queue>
```

For example, for 48 CPUs, 192 GiB, two hours, and the `normal` queue:

```bash
qcost -c 48 -m 192 -w 2 -q normal
```

## Job Charging

Gadi charges jobs according to the resources **reserved** and the walltime
**used**. The general formula is:

```text
Job Cost (SU) = Queue Charge Rate
              × max(NCPUs, Memory Proportion)
              × Walltime Used (hours)
```

Where:

- **Queue Charge Rate** depends on the selected queue; express queues cost
  more per SU but increase priority.
- **NCPUs** is the number of requested CPUs.
- **Memory Proportion** is the requested memory divided by the queue's
  memory-per-core allowance (memory per node divided by CPUs per node).

Requesting an outsized amount of memory can therefore increase the cost even
if the job does not use it. A job requesting enough memory to fill a node
effectively reserves the whole node and is charged accordingly.

See the PBS job-cost documentation in the
[Gadi User Guide](https://opus.nci.org.au/spaces/Help/pages/90308778/0.+Welcome+to+Gadi)
for the current queue charge rates and limits.

### Practical Implications for LPJ-GUESS

LPJ-GUESS uses well under 1 GiB of memory per CPU. Its costs are therefore
normally dominated by CPU count rather than memory.

The `normal` queue nodes have 48 CPUs and 192 GiB RAM, giving a
memory-per-core ratio of 4 GiB/CPU. Because LPJ-GUESS memory use is below this
threshold, the billing formula simplifies to:

```text
Job Cost (SU) = Queue Charge Rate × NCPUs × Walltime Used
```

Requesting less than 4 GiB per CPU on the `normal` queue costs the same as
requesting 4 GiB per CPU. **Request 4 GiB per CPU as a safe default**; it does
not cost more and avoids unexpected out-of-memory failures.

## Estimating LPJ-GUESS Run Cost

LPJ-GUESS cost on Gadi is approximately **0.005–0.01 SU per patch**. For a run
over G gridcells, each containing P patches:

```text
estimated SUs ≈ G × P × 0.005   (lower bound)
estimated SUs ≈ G × P × 0.01    (upper bound)
```

The cost per SU is approximately **$0.04 AUD** (indicative—verify against
current NCI pricing). These are first-order estimates assuming uniform patch
runtime. More productive gridcells (higher biomass turnover) take slightly
longer, but the effect is small and `choose_num_cpus.py` can minimise imbalance
across workers.

After a job finishes, compare this estimate with the actual cost recorded in
the top-level job log. See
[Checking and Troubleshooting Gadi Runs](gadi-troubleshooting.md#start-with-the-job-and-worker-logs).
