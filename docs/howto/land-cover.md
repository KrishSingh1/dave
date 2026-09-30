# Using land-cover change in DAVE

This guide configures a DAVE site simulation with natural vegetation, managed
forest, and optionally pasture. It assumes the usual non-crop DAVE workflow.
For the model concepts, see [Land cover and land-use change](../science/land-use-change.md).

## Choose the mode

The committed `benchmarks/ozflux/settings.ins` imports
`benchmarks/ozflux/landcover.ins`, which sets `run_landcover 0`. Use that setting
for a site with one ordinary stand and no prescribed land-cover fractions.
The `file_lu` parameter is ignored in this mode.

For a run with prescribed fractions, add a site-specific instruction file
*after* the shared OzFlux settings import and set `run_landcover 1`. Turn on
only the classes to be simulated, and define at least one included `st` for
each active class. A class can have several stand types; the `st` name is used
in stand-type fraction and gross-transition column headers. Keep crop disabled
for a DAVE run without crops.

```ins
! In the site-specific instruction file, imported after ../settings.ins
run_landcover 1
run_natural 1
run_forest 1
run_pasture 0
run_crop 0
run_urban 0
run_peatland 0
run_barren 0

st "Forest" (
    common_stand
    stinclude 1
    landcover "forest"
    naturalveg "all"
    harvest_system "continuous"
)
```

The shared OzFlux land-cover file already defines `st "Natural"`. Its exact
name matters for a `Natural_Forest` gross-transition column. If you enable
pasture, also define an included pasture stand type and its appropriate PFTs
and grazing settings. A managed forest stand type needs management choices
suited to the experiment; the snippet is an area-configuration example.

To initialise a forest stand from observed tree cohorts in a chosen year, see
[Prescribed tree establishment](prescribed-tree-establishment.md).

## Supply the net fractions

Set `file_lu` to a text file of grid-cell fractions. For an annual site input,
use a `Lon Lat Year` header followed by land-cover names in uppercase. Match
the longitude and latitude to the site's `file_gridlist`. Values are fractions
of the entire grid cell, from 0 to 1. Give one row per year and location; keep
the classes summing to 1. This minimal non-crop example moves 0.10 of the grid
cell from natural vegetation to forest in 1902:

```text
Lon Lat Year NATURAL FOREST
150.723611 -33.615278 1901 1.00 0.00
150.723611 -33.615278 1902 0.90 0.10
150.723611 -33.615278 1903 0.90 0.10
```

```ins
param "file_lu" (str "my-landcover.txt")
```

Omitted inactive classes contribute zero. Include a column for every active
class, even when its fraction is zero in some years. The reader supports other
text input layouts, including year-only data shared by grid cells; the
`Lon Lat Year` form makes a site file self-contained. File paths are resolved
from the run's working directory, so run from the site directory as described
in [Running Serially](run-serial.md).

When `run_landcover 1` has no `file_lu`, the model gives the active classes
equal fixed fractions. That is useful for a controlled test, but is usually not
useful behaviour. If only `NATURAL` is active, it occupies the whole grid cell.
When active fractions do not sum to one, the reader may rescale them or assign
unrepresented area to `NATURAL`; do not rely on that correction for intended
land use. Inspect log warnings and `frac_sts.out`.

## Decide whether to give gross transitions

For a straightforward conversion with one donor and one receiver, `file_lu`
alone can be sufficient. DAVE infers a transfer from the net change. For
several simultaneous changes, or for turnover with unchanged net fractions,
give an explicit gross-transition file.

For stand-type transitions, set `file_grossSTC` and use columns named
`FromStandType_ToStandType`, with the exact `st` names and case from the
instruction file. Its `Year` is the year **before** the target net-fraction
row:

```text
Lon Lat Year Natural_Forest Forest_Natural
150.723611 -33.615278 1901 0.10 0.00
150.723611 -33.615278 1902 0.00 0.00
```

```ins
param "file_grossSTC" (str "my-stand-transfers.txt")
param "file_grossLUC" (str "")
```

The 1901 gross row creates the 1902 net fractions in the example above. A
missing transfer column is treated as zero, so check spelling carefully.
Each transfer is a fraction of the entire grid cell. Before running, verify
for each year and stand type that **incoming minus outgoing transfers equals
the change in its net fraction**. For a class with several stand types, also
check that their target fractions add to the class fraction.

`file_grossLUC` is a different, LUH-style, **land-cover-level** transfer
format with short category columns such as `vp` and `pv`; it is not a table
of `Natural_Forest` columns. DAVE maps these class transfers to stand types
and may adjust inconsistent gross transfers to match net fractions. Use
`file_grossSTC` when the donor and receiver stand types are known. If both
gross paths are configured, `file_grossLUC` wins, so keep the unused path
empty.

If a land-cover class contains multiple stand types, the optional
`file_lunatural`, `file_luforest`, or `file_lupasture` supplies their fractions
over time. Without such a file, active stand types in that class have equal
shares. A stand-type fraction file uses stand-type names as columns. The
reader infers whether its values are shares *within* the class or absolute
grid-cell fractions by checking whether the numeric columns sum to one. Use
one convention consistently across cells and years, and confirm the reported
stand-type areas in `frac_sts.out`. For a single stand type in each class,
these extra files are unnecessary.

## Stand transfer settings

The shared `benchmarks/ozflux/landcover.ins` defines the following settings.
Override only those needed by the experiment:

| Setting | Effect |
| --- | --- |
| `npatch_secondarystand` | Number of replicate patches in a newly created secondary stand. |
| `reduce_all_stands` | If `0`, use age and destination rules to choose donor stands; if `1`, reduce stands more evenly. |
| `age_limit_reduce` | Age threshold used when selecting stands to shrink. |
| `transfer_level` | Pool transferred state across the grid cell (`0`), land-cover class (`1`), or stand type (`2`). |
| `iftransfer_to_new_stand` | Allow eligible converted areas to become separate stands, with vegetation carried over or cleared according to the receiving stand's settings. |
| `ifslowharvestpool` | Include the slow harvested-product pool in land-cover harvest accounting. |

`ifprimary_lc_transfer` and `ifprimary_to_secondary_transfer` are for detailed
primary/secondary transition data; leave them at `0` for the simple
stand-type example. `harvest_secondary_to_new_stand` controls a wood-harvest
path that can create secondary stands, even with no change in broad class
area. Harvest, planting, and disturbance settings can therefore change the
result of an area transition. See [Disturbance](../science/disturbance.md) when
prescribed fire, wind, bark beetle, or clearcut input is also active.

## Run and inspect

Run the site instruction file from its directory. With the standard build and
NetCDF climate input, the form is:

```bash
cd benchmarks/ozflux/<SiteName>
../../../build/guess -input nc <SiteName>.ins
```

Add output settings to the site-specific instruction file. These are useful
starting points:

```ins
file_frac_sts "frac_sts.out"
file_nstand_sts "nstand_sts.out"
file_cmass_sts "cmass_sts.out"
file_luc_cflux_sts "luc_cflux_sts.out"
file_cflux_natural "cflux_natural.out"
file_cflux_forest "cflux_forest.out"
```

`frac_sts.out` reports stand-type fractions, class fractions, and a grid-cell
total. `nstand_sts.out` reports how many runtime stands represent each type;
its count can increase after a transfer even if the type's area does not.
`cmass_sts.out` helps distinguish area changes from changes in vegetation
carbon. `luc_cflux_sts.out` and the class carbon-flux files help check the
carbon consequences of conversion. Set `printstandtypes 1` for stand-type
columns; `printseparatestands 1` requests individual stand output where
supported, but the model disables it for land-cover files covering more than
100 grid cells to avoid too many files. Output is written under the configured
`outputdirectory`.

Check the first input year, the first year in which fractions change, and the
year after the last supplied row. Compare `frac_sts.out` with the intended
fractions, and check that an explicit gross row dated $y-1$ produces the
expected transition into year $y$. During spin-up, land-cover changes can be
held until the soil spin-up threshold; read the log for the first active year.

## Common input problems

If the fractions, transfers, or carbon fluxes look wrong, check `guess.log`
in the site run directory **first**. Search for the affected year, stand type,
rescaling messages, warnings, and errors. Set `verbosity 3` in the `.ins` file
for informational land-cover messages; `4` or `5` can provide more detail.
See [Start with `guess.log`](common-failures.md#start-with-guesslog).

| Symptom | What to check |
| --- | --- |
| `All active land covers must have at least one stand type` | Define an included `st` with the matching `landcover` for every `run_* 1` class. |
| Missing input or a site excluded | Check `file_lu` path, `Lon`/`Lat` against `file_gridlist`, and that the file has data for the location. |
| Fraction-sum or rescaling warning | Check all active columns, their units, and that each year's fractions sum to one. |
| No change in the expected year | Check the initial row, spin-up timing, and the one-year offset of the gross file. |
| Wrong or missing donor stand | Check gross column names against `st` names and balance gross transfers against net stand-type fractions. |
| Unexpected carbon flux or number of stands | Inspect harvest, disturbance, `transfer_level`, and `iftransfer_to_new_stand` as well as the area files. |

The repository's `benchmarks/ozflux/CumberlandPlain/` directory contains a
larger `cup.luc.txt` and `cup.gross_stc.txt` example with both directions of
natural/forest conversion. Its site instruction imports `cup.luc.ins` after
the shared settings; use it as a worked example and adjust the years, location,
management, and input paths for your own site.
