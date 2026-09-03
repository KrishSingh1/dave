# Setting Up Gadi

[Gadi](https://nci.org.au/our-systems/hpc-systems) is the HPC system at the
National Computational Infrastructure (NCI) in Canberra. This page covers the
one-time setup and background needed to run LPJ-GUESS DAVE jobs on Gadi. Once
your account and environment are ready, continue with
[Running on Gadi](run-gadi.md).

This guide is not a substitute for the comprehensive
[Gadi User Guide](https://opus.nci.org.au/spaces/Help/pages/90308778/0.+Welcome+to+Gadi).
Reading at least its sections on PBS, storage, and job debiting before your
first run is strongly recommended.

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
2. Register the **public key** (`id_ed25519.pub`).

Log in to Gadi with your password, then append your public key to
`~/.ssh/authorized_keys` on the server. The easiest way is `ssh-copy-id`
(Linux/macOS):

```bash
ssh-copy-id -i ~/.ssh/id_ed25519.pub <nci_username>@gadi.nci.org.au
```

You will be prompted for your NCI password once. After that, SSH and SFTP
connections will use the key instead.

If `ssh-copy-id` is not available (e.g. on Windows), log in manually and run:

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
it to PBS, which queues and executes it when resources become available.

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
| `qstat -f <jobid>` | Show detailed status for a specific job |
| `qdel <jobid>` | Cancel a job |

No single PBS job on Gadi may run for longer than 48 hours. See
[Choosing Gadi Resources and Estimating Costs](gadi-resources.md) for guidance
on selecting CPUs, memory, and walltime.

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
    particularly common and hard-to-diagnose error. The `build_gadi` and
    `submit_to_gadi.sh` scripts ensure consistency; use them rather than loading
    modules manually.

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
provides Gadi-compiled versions of standard LPJ-GUESS utilities such as
`aslice`, `tslice`, and `gmap`.

Many of the shared helper scripts are maintained in the
[`hol430/gadi-scripts`](https://github.com/hol430/gadi-scripts) repository.

After editing `~/.bashrc`, either log out and back in, or source it in the
current shell:

```bash
source ~/.bashrc
```

## Storage

Gadi provides three types of storage, with very different properties.

### Home Directory (`$HOME`)

`$HOME` (`/home/xxx/<nci_username>`) is limited to **10 GiB**. It is intended
for configuration files, scripts, and source code—**never for job outputs**. A
run writing even modest output to `$HOME` will quickly exhaust the quota and
cause jobs to fail.

### Scratch (`/scratch/<project>/$USER`)

`/scratch` is the primary working space for job input and output. It is fast
and has a larger quota, but files that have not been accessed or modified for
**100 days** are automatically deleted without warning.

Files are first moved to a quarantine area before permanent deletion. Consult
the Gadi User Guide linked above for instructions on listing and recovering
quarantined files.

!!! warning "No permanent storage"
    The DAVE project does not have permanent storage on Gadi (only scratch).
    Important results should be backed up to [HIE-Storage](hie-storage.md) as
    soon as possible after a run completes.

### Gdata (`/g/data/<project>`)

`/g/data` is persistent project storage not subject to the 100-day deletion
policy. The DAVE project does not currently have any gdata allocation.

### PBS Filesystem Access

PBS jobs must explicitly declare which filesystems they use, for example:

```bash
#PBS -l storage=scratch/pt17+gdata/pt17
```

Accessing an undeclared path fails at runtime even if you can reach it
interactively. The `submit_to_gadi.sh` script handles this automatically for
standard runs; only manually written PBS scripts need to consider it.

### Check Storage Quotas with `lquota`

The shared DAVE environment provides `lquota`, which prints storage usage and
limits for each project you are a member of:

```bash
lquota
```
