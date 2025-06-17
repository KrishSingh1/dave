# The Global Tellme Script
# Matt Forrest 2024-04-07


## About

The global Tellme R markdown script (tellme_global.Rmd) is a benchmarking routine to evaluate the  `tellus` global benchmark simulation (which is intended to be the best possible global run with all features turned on).  As such, it benchmarks global GPP, LAI, burnt area, biomass, etc., and also compares the `tellus` simulation to another reference `tellus` simulation.  In addition to this, the script also optionally uses some other simulations for some other benchmarks: the `fluxnet` and `regrowth` simulations for their respective benchmarks (although at present neither of these are operational) and the `global` simulation for PNV and applying a biome specific correction to biomass (total to AGB) in the global biomass benchmark.

Within the LPJ-GUESS benchmarking system the script is run as a batch job on an HPC as a "summary job".  Typically this will be submitted at the same time as a suite of bencmarking runs (ie. `tellus` and others) and then is automatically run afterwards.  The script can also be run offline (for example in RStudio on your local PC).  This is straightforward but that usage is not covered here.

The script is controlled by a combination of command line arguments which specific so key inputs, and a YAML file (`tellme_global_config.yml`) which specifies a lot of additional customisability.  The  `tellme_global_config.yml` included here contains reasonable default values, but can be modified by the user.

## Usage 

#### Usage 1 - as a batch job as part of the full benchmarking suite on an HPC.

To run the `tellme` benchmark and the appropriate simulations.
1.  Download the code to a supported HPC (currently Lunarc and Levante)
1.  In the `benchmarks` folder run:
```
# to run the tellus simulation, the global simulation (also used by tellme script) and the tellme script.
./benchmarks -i "tellus global" -s "tellme <MyRunNameTag> <ReferenceRunNameTag> <Path/to/ReferenceSimulations>" <Path/to/Benchmaks/OutputDirectory>
# e.g.
./benchmarks -i "tellus global" -s "tellme new_feature trunk /data/mforrest/benchmarks/trunk_v4.1" /data/mforrest/benchmarks/new_feature
```
Note that you have free choice when it comes to `<MyRunNameTag>` and `<ReferenceRunNameTag>`.  <Path/to/ReferenceSimulations> can be ommited but otherwise is the path to the top level directory of a previous benchmark run.  <Path/to/Benchmarks/OutputDirectory>  is where the sumilations are performed and the output generated - *it is not recommended to do this on `/home` as likely you will use up your disck space quota.*

#### Usage 2 - as a single threaded job on any machine.

It is also possible to run the script through the benchmarking suite as a stand alone post-processing job *after* the `tellus` (and other simulations) have been performed.  To do this: 

```
# to run the tellme script stand alone after the tellus simulation as completed:
./benchmarks -1 -P -s "tellme <MyRunNameTag> <ReferenceRunNameTag> <Path/to/ReferenceSimulations>" <Path/to/Benchmaks/OutputDirectory>
# e.g.
./benchmarks -1 -P -s "tellme new_feature trunk /data/mforrest/benchmarks/trunk_v4.1" /data/mforrest/benchmarks/new_feature
```
The `-1` means run as a single local process (ie don't submit as a batch job) and `P` means do summary post-processing jobs only.  Please note that the tellme script is quite memory demanding, it currently needs something like 12 GB of RAM.

#### Usage 3 - as a batch job on an HPC, but after the main benchmark suite has been completed.

This is kind of a hybrid between the two above cases, where the script is run as a batch job, but after the benchmarks runs have already been done:
```
# to run the tellme script as a batch job, but after the simualtions have been completed:
./benchmarks -P  -s "tellme <MyRunNameTag> <ReferenceRunNameTag> <Path/to/ReferenceSimulations>" <Path/to/Benchmaks/OutputDirectory>
```
Here we just omit the `-1` flag and the script will be submitted as a batch job.


## Known issues/upcoming development

1. FLUXNET benchmark is not currently working and will be overhauled.
1. The `regrowth` benchmark simulation is not currently implemented. 
1. At time of writing, the current of version DGVMBenchmarks has a regression whereby the difference and percentage difference row are not being displayed in some comparison tables.  This is under investigation and should be fixed soon.
1. The actual `tellme_global.Rmd` script is currently under redevelopment to achieve two things:
   * Remove all the per-benchmark hard coded settings to the YAML configuration file (currently done for three benchmarks). 
   * Factorise the code out from the script into functions from the script into functions in the DGVMBenckmarks package.  This has only been done for one benchmark so far.
Completing these two goals is planned but with no fixed timetable.  As the things stand the benchmarks are operational and the planned changes are purely technical rather than scientific, so they don't prohoibit the use of the benchmarks now.


## Setup Instructions for a new HPC

In order to run the Tellme benchmarking within the benchmarking framework of LPJ-GUESS it is neccesary to install some packages R and make sure that some additional software is available, this is usually done by loading 'modules' with `module` commands, but this may depend on the HPC in question.

This guide will outline the general principles, but also the specific commands for certain machines, initially Simba, DKRZ's Levante and Lunarc's Cosmos.


### Load system modules/external software


First we need to make sure that certain software is available to us: 
 * R
 * the GCC compiler
 * Pandoc
 * OpenMPI (although I can't remember why at this stage)
 

These will likely be available as modules on the HPC.  But they might simply need to be installed on some machines.
 

**Cosmos (Lunarc)**:
```
> module purge
> module load GCC/11.3.0
> module load OpenMPI/4.1.4
> module load R/4.2.1
> module load Pandoc/3.1.2
```

**Simba**: 
```
fill_me_in_Johan!
```

**Levante (DKRZ)**:
```
> module purge
> module load gcc/11.2.0-gcc-11.2.0
> module load r/4.1.2-gcc-11.2.0
> module load openmpi/4.1.2-gcc-11.2.0
```


### Install or update R-packages, DGVMTools and DGVMBechmarks

Here we work within R to install first some all the R package dependencies of the tellme_global.Rmd script.
The same procedure is used when you update DGVMTools and DGVMBechmarks. 
The installation/update of DGVMTools takes a long time.

##### R and devtools

First, start up R and install devtools: 
 
``` 
> R
> install.packages("devtools") # just in case not already installed
```
On cosmos (lunarc) we have installed R in the common directory
/lunarc/nobackup/projects/snic2020-6-23/lpjguess/R/4.1/

##### DGVMTools and DGVMBenchmarks

Now install DGVMTools and DGVMBencharks.


```
# choose a mirror to avoid the popup window
> options(repos = 'http://ftp5.gwdg.de/pub/misc/cran/')
# this is for DGVMTools
> devtools::install_github("MagicForrest/DGVMTools", ref = "master", dependencies = c("Depends", "Imports"), build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
# This is for DGVMBenchmarks:
> devtools::install_github("MagicForrest/DGVMBenchmarks", ref = "main", dependencies = c("Depends", "Imports"), build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
```

The commands above will ask about which dependencies to update.   I recommend installing updated versions of all the packages when asked, but be warned this might take a while.  Also there might be issues, see **Troubleshooting** later.


##### Dependencies for the tellus.rmd markdown script

For making the benchmark script work, we also need to install the following packages (still inside R):

```
> install.packages("kableExtra")
> install.packages("tictoc")
> install.packages("pals")
> install.packages("dplyr")
> install.packages("showtext")

```

All going well we can exit R:

```
> q()
```



##### Troubleshooting R package installs

Unfortunately, there might be some complications with installing some R packages.  For example, some of the required R-packages require linking to external software (such as udunits2 and GDAL) at installation time in a rather specific manner.  This may be highly system specific.

Example: Levante (DKRZ) needs a special install command for the `units` package:
```
install.packages("units", configure.args = list(units=c("--with-udunits2-lib=/sw/spack-levante/udunits-2.2.28-da6pla/lib", "--with-udunits2-include=/sw/spack-levante/udunits-2.2.28-da6pla/include", "LIBS=-Wl,-rpath,/sw/spack-levante/udunits-2.2.28-da6pla/lib")))
```
In some cases it might even be necessary to install certain packages and libraries that the R packages depend upon.  But all this needs to be solved on case-by-case basis, likely with the help of your HPC admin.






