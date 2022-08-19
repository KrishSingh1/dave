## Eucface Benchmark

The LPJ-Guess Eucface Benchmark.

### Configuration

All configuration is in the .ins files inside the config directory. The main
.ins file which should be passed as a command-line argument to LPJ-Guess is
called `eucface.ins`. This .ins file contains only 3 directives, importing the
other 3 .ins files in the config directory, which contain the actual
configuration required by LPJ-Guess.

**eucface_files.ins** contains the path to the input (data) files

**eucface_settings.ins** contains the global parameters and PFT definitions

**eucface_soiln.ins** contains some soil N config

### Data

The input files specified by `eucface_files.ins` should be located in a
directory called data, at the same level as the config directory
(ie `benchmarks/eucface/data`). The input files are all stored on our
hie-storage directory, so I'd recommend mounting that via sshfs, or, if
you've already mounted it somewhere else on your computer (even better!),
then just create a symlink to the mountpoint:

```bash
ln -s /path/to/mounted/hie-storage /path/to/lpjguess_dave_eucface/benchmarks/eucface/data
```

**Soil Data**

LPJ-Guess requires the following soil data:

**lon** - Provided in eucface MIP (-33.61°)
**lat** - Provided in eucface MIP (150.73°)
**clay** - Provided in eucface MIP (0.11)
**silt** - Provided in eucface MIP (0.09)
**sand** - Provided in eucface MIP (0.8)
**orgC** - Provided in eucface MIP (Soil carbon fraction for top 0-10, 10-20, 20-30cm: 1.55, 0.71, 0.42. Avr 0.89)
**CN** - Not actually used anywhere afaik
**pH** - (Provided in eucface MIP (top 0-10, 10-20, 20-30cm: 5.52, 5.39 , 5.28. Avr 5.40)
**cellfraction** - Not actually used anywhere afaik

In that sense,the soil input data file is not actually necessary when running with `ifeucface 1`.

I note that LPJ-Guess supports setting initial soil C content and bulk density
in the soil map file. We're currently not doing this.

### Running

You need to use the siteinput module, which was created by David Warlind
specifically for the Eucface run:

```bash
guess -input site /path/to/benchmarks/eucface/config/eucface.ins
```

This should fail quickly if your input data is not setup correctly.

### Debugging

I've created a debug profile for VS code and committed it to this branch. To
debug the eucface run, open the repository directory **as a folder** in vscode,
then hit f5. This will launch the eucface benchmark.

I've configured this so that lpj-guess will be rebuilt every time
you do a run. Thanks to the magic of makefiles, only the changed code files will
rebuilt. So if you haven't changed anything, this will be very fast (after the
first compilation).
