# Setup Instructions for Tellme/Tellus on an HPC
# Matt Forrest 2024-02-07.

In order to run the Tellme benchmarking within the benchmarking framework of LPJ-GUESS it is neccesary to install some packages R and make sure that some additional software is available, this is usually done by loading 'modules' with `module` commands, but this may depend on the HPC in question.

This guide will outline the general principles, but also the specific commands for certain machines, initially Simba and Aurora.


## Load system modules/extrenal software


First we need to load make sure that certain software is available to us: 
 * R
 * the GCC compiler
 * Pandoc
 * OpenMPI (although I can't remember why at this stage)
 

These will likely be available as modules on the HPC.  But they might simply need to be installed on some machines.
 

**Aurora**:
```
> module load GCC/10.3.0  OpenMPI/4.1.1  R/4.1.0  Pandoc
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


## Install R-packages, DGVMTools and DGVMBechmarks

Here we work within R to install 

##### R and devtools

First, start up R and install devtools: 
 
``` 
> R
> install.packages("devtools") # just in case not already installed
```

##### DGVMTools and DGVMBenchmarks

Now install DGVMTools and DGVMBencharks.


```
# choose a mirror to avoid the popup window
> options(repos = 'http://ftp5.gwdg.de/pub/misc/cran/')
# this is for DGVMTools
> devtools::install_github("MagicForrest/DGVMTools", ref = "master", dependencies = c("Depends", "Imports"), build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
# This is for DGVMBenchmarks:
> devtools::install_github("MagicForrest/DGVMBenchmarks", ref = "master", dependencies = c("Depends", "Imports"), build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
```

The commands above will ask about which dependencies to update.   I recommend installing updated versions of all the packages when asked, but be warned this might take a while.  Also there might be issues, see **Troubleshooting** later.


##### Dependencies for the tellus.rmd markdown script

For making the benckmake script we also need to install the following packages (still inside R):

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

Unfortunately, there might be some complications with installing the dependencies.  For example, some of the required R-packages might require specific linking to external software (such as udunits2 and GDAL).

Example: Levante (DKRZ) needs a special install command for the `units` package:
```
install.packages("units", configure.args = list(units=c("--with-udunits2-lib=/sw/spack-levante/udunits-2.2.28-da6pla/lib", "--with-udunits2-include=/sw/spack-levante/udunits-2.2.28-da6pla/include", "LIBS=-Wl,-rpath,/sw/spack-levante/udunits-2.2.28-da6pla/lib")))
```
In some cases it might even be necessary to install certain packages and libraries.  But all this needs to be solved on case-by-case basis.



### Run the Tellme/Tellus benchmark

Now we can go ahead and try to run the benchmark.  Checkout the relevent branch and change `benchmarks` directory.

```
./benchmarks -i "tellus" -e "tellme <MyRunNameTag>  <ReferenceRunNameTag> <Path/to/reference-run/output> <Path/to/where/I/want/MyRun/to/be/outputted>
```
Note that you have free choice when it comes to `<MyRunNameTag>` and `<ReferenceRunNameTag>`.






