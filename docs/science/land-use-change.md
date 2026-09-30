# Land cover and land-use change

This page describes how DAVE represents land cover and changes in land use. For a configuration
guide and input examples, see [Using land-cover input](../howto/land-cover.md).

## The three levels of area

| Level | Meaning | Example |
| --- | --- | --- |
| Grid cell | The location and its total simulated area. Area fractions are measured against this total. | One OzFlux site. |
| Land-cover class | A broad use or cover category within the grid cell. | `NATURAL`, `FOREST`, or `PASTURE`. |
| Stand type and stand | A **stand type** (`st` in an instruction file) is a named template assigned to a land-cover class. A **stand** is a runtime area with that template, its own patches, vegetation, soil, and history. | A `Forest` stand type may have several stands of different ages. |

A stand contains one or more replicate **patches**. Patch count and stand area
are different things: doubling `npatch_secondarystand` adds replicates to a new
stand; it does not double the area transferred. Plant functional types (PFTs)
determine what can grow within each stand, while the stand type can specify
management and establishment rules.

For land-cover class $c$, let $f_c$ be its fraction of the grid cell. The
fractions of all included land-cover classes normally add to one. For a stand
type $s$ within class $c$, let $f_s$ be its grid-cell fraction. Then

\[
\sum_{s \in c} f_s = f_c, \qquad \sum_c f_c = 1.
\]

The grid-cell fractions of runtime stands belonging to a stand type add up to
that stand-type fraction. This is why
the number of stands can grow even when the area of a land-cover class stays
constant.

The code recognises `URBAN`, `CROPLAND`, `PASTURE`, `FOREST`, `NATURAL`,
`PEATLAND`, and `BARREN`. A DAVE run can enable only the classes it needs.
`FOREST` is the managed-forest class; `NATURAL` is natural vegetation, which can
also contain trees. A stand type's `landcover` setting assigns it to one of
these classes. Enabling a class requires at least one included stand type for
that class.

## Static sites and dynamic land cover

`run_landcover 0` disables landcover change over time. The committed ozflux
settings use this mode.

`run_landcover 1` activates the land-cover workflow. At initialisation,
LPJ-GUESS reads the class and stand-type fractions and creates stands with
positive area. At the start of each subsequent model year it updates the
fractions, determines transfers, reduces donor stands, expands or creates
receiving stands, and moves the relevant carbon, nitrogen, water, litter, and
soil state. If the fractions and gross transfers do not change, the land-cover
area transfer can be skipped. The switch alone does **not** make land cover
change: constant input fractions give a constant land-cover map, although other
management pathways may still act on stands.

The fraction file describes **how much area each class has**. A gross transfer
file, if supplied, describes **where area came from and where it went**. These
are separate inputs. A transition of 0.1 means 10% of the *whole grid cell*,
not 10% of the donor class.

## Net change and gross transfers

For year $y$, the net change of class $c$ is its new fraction minus its old
fraction. If $T_{a,b}$ is the gross grid-cell fraction transferred from class
or stand type $a$ to $b$, then a consistent set of inputs satisfies

\[
f_c(y)-f_c(y-1) = \sum_a T_{a,c} - \sum_b T_{c,b}.
\]

For example, moving 0.10 from `NATURAL` to `FOREST` changes their fractions
from 1.00/0.00 to 0.90/0.10. A simultaneous 0.05 transfer in *each* direction
has zero net class change but can still create or replace stands and affect
vegetation and carbon fluxes. Net fractions alone cannot reconstruct that
two-way turnover.

When no gross transfer file is provided, LPJ-GUESS constructs transfers from net
changes. A case with one shrinking and one growing class is direct. With several
donors and receivers, the model uses hardcoded preferences to assign origins and
destinations. This is a model assumption, so use explicit gross transfers when
the origin of converted land matters. Land-cover gross input (`file_grossLUC`)
uses the model's LUH-style category columns and is mapped to stand types.
Stand-type gross input (`file_grossSTC`) names the donor and receiver stand
types directly. If both paths are specified, the land-cover gross file takes
precedence.

Gross transitions can exist within one land-cover class: for example, between
two forest stand types, or between old and newly regrowing natural stands.
These can change stand history while total class area remains fixed.

## What happens to a transferring stand?

The yearly update selects area from donor stands and assigns it to receiving
stands. With multiple donors of the same stand type, the selection depends on
the destination and the configured age and primary/secondary rules. By default,
`reduce_all_stands 0` uses that selection logic; `age_limit_reduce` influences
which stands are eligible first. A new receiving stand normally uses
`npatch_secondarystand` patches. Depending on transfer settings and stand types,
LPJ-GUESS may create a fresh stand, expand an existing one, or clone donor
vegetation into a new stand. A stand whose area falls to zero is removed.

Land transfer also has a **state and flux** side. Donor vegetation may be
harvested or killed according to the transition and management settings.
Litter, soil carbon and nitrogen, and water are transferred or pooled into the
receiving area. `transfer_level` selects the pooling scale: `0` pools over the
grid cell, `1` by land-cover class, and `2` by stand type. The exact biomass
fate therefore depends on the chosen transition, stand type, and harvest
settings; an area fraction change is not equivalent to simply relabelling
vegetation.

`iftransfer_to_new_stand 1` allows qualifying land transfers to form distinct
new stands. This keeps the converted area and its subsequent history separate
from other stands of the receiving land-cover class. Depending on the source
and the receiving stand's management settings, existing vegetation may carry
over, or trees may be cleared and the new stand may regrow. The setting does
not force every transfer to create a separate stand.

`ifprimary_lc_transfer` and
`ifprimary_to_secondary_transfer` enable the more detailed primary/secondary
transition interpretation used with suitable gross land-cover input. A simple
site file with only `Natural_Forest` columns does not need these options.

Wood harvest and prescribed disturbance can also split or create stands
without a net change of land-cover class. They are separate drivers from the
`file_lu` class fractions; see [Disturbance](disturbance.md) for the prescribed
disturbance pathway.

Prescribed tree establishment can replace vegetation in a stand type in a
chosen year, and can optionally pool soil and recreate its stands. See
[Prescribed tree establishment](../howto/prescribed-tree-establishment.md)
for its input format and effect on stands.

## Timing and spin-up

Land-cover dynamics run on the first day of a model year, before the year's
vegetation and forest management steps. For annual net fractions, the row for
calendar year $y$ gives the target area used in year $y$. Explicit gross
transition input is read from row $y-1$ when the model enters year $y$:
the transition is considered to have happened at the end of the preceding
year. Thus a move from 1.00 natural in 1901 to 0.90 natural and 0.10 forest
in 1902 belongs in the **1901** gross-transition row.

During spin-up, the text input machinery holds an available boundary year's
fractions outside the input range. The land-cover code also prevents changes
before its soil spin-up reduction threshold. Check the run log to see when
changes first become active; this depends on the simulation's spin-up and
input years. A run that begins at the first historical input year should have
an initial net-fraction row, with gross transitions starting on the move to
the following year.

## Interpretation

- Fraction outputs measure area, while carbon and water outputs measure the
  states or fluxes of the resulting stands. A fraction change may alter a
  grid-cell average even if per-area vegetation behavior is unchanged.
- Class fractions, stand-type fractions, and runtime stand counts answer
  different questions. A single `FOREST` class can contain several stands of
  one `Forest` stand type.
- Land-cover conversion, clearcut, fire, and generic disturbance have distinct
  rules for allocating area and biomass. Compare their output streams before
  attributing a carbon flux to a particular cause.

The relevant implementation is in `framework/externalinput.cpp` (fraction and
gross-input reading), `modules/landcover.cpp` (stand transfers), and
`framework/framework.cpp` (yearly call order).
