# Prescribed tree establishment (treeinit)

Treeinit places specified tree cohorts into a stand type in a chosen calendar
year. Use it to initialise a stand from an observed tree size and density
structure. It does not prescribe land-cover area or replace ordinary recruitment
in later years. See [Establishment](../science/establishment.md) for annual
recruitment and [Land cover and land-use change](../science/land-use-change.md)
for stand areas.

## What happens

At the **end of `init_year`**, after annual vegetation dynamics and prescribed
disturbance, LPJ-GUESS removes the existing vegetation in each affected patch
and creates the trees described by the input rows. Existing trees are cut and
their biomass is divided among exported material and litter; existing grasses
become litter. [Harvesting and vegetation removal](../science/harvesting.md)
explains this accounting. The new trees receive size, density, age, carbon,
and nitrogen state. Their added carbon and nitrogen are recorded as prescribed
inputs through seed flux accounting, rather than as growth during that year.
Ordinary annual establishment is skipped for this stand type in the
initialisation year. The affected stands then move to `management2`, the
regime governing later cutting, planting, and establishment. See
[Management types](management.md) for the concepts and parameter reference.

With the default `forest_init_poolsoil 0`, existing stands keep their area,
patches, and soil state while their vegetation is replaced. With
`forest_init_poolsoil 1`, LPJ-GUESS pools the soil and other non-living state
across the stands being initialised, removes those stands, and creates one new
stand per affected stand type. The new stands have `npatch_treeinit` replicate
patches. Their area can be redistributed among affected stand types using
next-year fractions while their combined area is preserved. All affected stand
types must have the same initialisation year when pooling is enabled.

## Configure the stand type

Set `file_treeinit_st` inside the **stand type** to be initialised. In the
land-cover workflow, use `run_landcover 1`, ensure the type has positive area
in `init_year`, and provide valid `management1` and `management2` regimes.
Treeinit updates the second regime's PFT selection to include the PFTs in the
input file. For example, in a site-specific instruction file imported after
the shared OzFlux settings:

```ins
run_landcover 1
run_natural 1
run_forest 1
run_crop 0

mt "Forest_before_init" (
    harvest_system "continuous"
    cutfirstyear 0
)
mt "Forest_after_init" (
    harvest_system "continuous"
    cutfirstyear 0
)

st "Forest" (
    common_stand
    stinclude 1
    landcover "forest"
    naturalveg "all"
    management1 "Forest_before_init"
    management2 "Forest_after_init"
    file_treeinit_st "forest-treeinit.txt"
)

forest_init_from_biomass 1
forest_init_poolsoil 0
```

The shared OzFlux file already defines `st "Natural"`. Set forest area with
`file_lu` as described in [Using land-cover change](land-cover.md). Ensure the
stand permits each tree PFT named in the treeinit file. Choose the two
management regimes for the experiment; treeinit does not prescribe future
cutting or recruitment.

## Prepare the input file

Use a whitespace-delimited file with one row per prescribed tree **cohort**.
The same cohort description is applied to every affected patch of the stand
type. `Tree_no` is a record index, starting at `0` and increasing by one; it
need not identify an observed individual tree. Match `Lon` and `Lat` to the
site grid list and `str_pft` to an included tree PFT exactly, including case.

```text
       Lon        Lat Tree_no init_year str_pft biomass  dbh dens tree_age
150.723611 -33.615278       0      2015     MRS   0.180 18.0  450       20
150.723611 -33.615278       1      2015     MRS   0.025  8.0  900        8
```

| Column | Meaning |
| --- | --- |
| `Lon`, `Lat` | Grid-cell location. Each grid cell using this stand type needs a matching block of rows. |
| `Tree_no` | Consecutive cohort record index, starting at `0`; the limit is 100 rows per location and stand type. |
| `init_year` | Calendar year of the one-time event. Use the same year in all rows; the first row sets it. |
| `str_pft` | Exact name of a tree PFT included in the run. |
| `biomass` | Aboveground dry biomass per representative tree, in tonnes. Required with `forest_init_from_biomass 1`. |
| `dbh` | Diameter at breast height, in centimetres. Used as the size input with `forest_init_from_biomass 0`. |
| `dens` | Number of trees per hectare represented by the cohort. |
| `tree_age` | Age in years. If absent or nonpositive, DAVE estimates it from diameter. |

Rows after the first describe additional cohorts for the **same event**, not
later initialisation years. Use a separate stand type for areas that require
different event years or management histories. Give every cohort a positive
size and density. The file must cover every simulated grid cell for which this
stand type has `file_treeinit_st`.

`forest_init_from_biomass` defaults to `1`. In that mode, DAVE derives wood
carbon and diameter from `biomass`; the supplied `dbh` is not the controlling
size value. With `forest_init_from_biomass 0`, DAVE uses `dbh` and PFT allometry
to calculate size and wood carbon. Supplying both columns makes the file
easier to inspect and use with either setting. Supply `dens` explicitly: its
units are converted internally from trees per hectare to trees per square
metre, and the fallback for a missing value is not a site observation.

## Run and verify

[Run the model](../howto/ozflux-simulations.md) with the stand type and input
file. The log reports the initialisation year and each prescribed cohort
created. Enable stand-type outputs to check the area, structure, and carbon
around the event:

```ins
printstandtypes 1
file_frac_sts "frac_sts.out"
file_nstand_sts "nstand_sts.out"
file_dens_sts "dens_sts.out"
file_diam_sts "diam_sts.out"
file_cmass_sts "cmass_sts.out"
file_dave_stand_frac "dave_stand_frac.out"
file_dave_stand_type "dave_stand_type.out"
```

Check that the stand has positive area in `init_year`, every intended cohort
appears in the log, and density, diameter, and carbon reflect the input.
Inspect the following year too, when `management2` is active. If soil pooling
is enabled, also inspect stand count and soil pools. The repository's
`benchmarks/ozflux/CumberlandPlain/cup.treeinit.txt` and `cup.luc.ins` provide
a worked example.

## Common problems

Start with `guess.log` in the directory from which you ran the site (for
example, `benchmarks/ozflux/CumberlandPlain/guess.log`). Search for the stand
type, PFT, `init_year`, warnings, and errors. If the reason is still unclear,
set `verbosity 3` or higher in the `.ins` file and rerun; levels `4` and `5`
give progressively more detail. The shared OzFlux settings already use `3`.
See [Start with `guess.log`](common-failures.md#start-with-guesslog) for a
general log-reading workflow. Check the log before using the symptom-specific
checks below.

| Symptom | Check |
| --- | --- |
| No prescribed trees appear | Positive stand area in `init_year`, simulation reaches that year, and `str_pft` matches an included tree PFT. |
| Input file fails to load | File path relative to the run directory, site coordinates, and row block for every simulated grid cell. |
| Missing biomass error | Supply `biomass`, or set `forest_init_from_biomass 0` and supply `dbh`. |
| Missing `management2` error | Define and assign a valid second management regime. |
| Cohorts appear wrong | Use consecutive `Tree_no` indices from `0`, one row per cohort, and no more than 100 rows. |
| Unexpected stand count or soil change | Check `forest_init_poolsoil` and `npatch_treeinit`. |
