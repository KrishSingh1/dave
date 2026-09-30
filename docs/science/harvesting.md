# Harvesting and vegetation removal

Harvesting in DAVE is an **accounting event** as well as a change to the trees
in a stand. Cutting removes some or all of the represented trees and divides
their carbon and nitrogen among material exported from the ecosystem, material
left as litter, and, when enabled, a longer lived harvested-product pool.
Exported carbon is reported as a harvest flux in the year of removal. The
product pool releases its carbon later. These bookkeeping categories do not
specify a real-world use for the wood.

## What determines the fate of a cut tree?

| Material | Typical treatment |
| --- | --- |
| Stem wood | A harvest-efficiency fraction is exported. Uncollected wood remains as woody litter. |
| Twigs and branches | A specified residue fraction can be removed; the rest becomes litter. Some leaves attached to removed twigs leave with them. |
| Coarse roots and stumps | A specified residue fraction can be removed; the rest becomes woody litter. |
| Leaves and fine roots | Usually enter litter, except for leaves removed with twigs. |
| Exported stem wood | If `ifslowharvestpool` is enabled, a PFT-specific share enters the harvested-product pool; the remainder is accounted for as a harvest carbon flux in that year. |

The fractions apply to the part of a tree that is cut. A thinning can cut a
fraction of the represented trees, whereas a clearcut normally removes all
eligible trees in the patch. Nitrogen follows corresponding exported or litter
paths. Stem, twig, and coarse-root shares depend on PFT traits, so equal
harvest settings can produce different fluxes for different species.

For ordinary forest management, `harv_eff_thin` and `harv_eff_cc` control the
stem harvest fraction during thinning and clearcutting. The matching
`res_outtake_twig_*` and `res_outtake_coarse_root_*` settings control residue
removal. If a management value is unset, stem efficiency and twig removal use
PFT settings, while coarse-root removal defaults to 0.1. The first scheduled
thinning in a clearcut rotation is precommercial: its cut biomass goes to
litter instead of exported products. See
[Management types](../howto/management.md) for the full
parameter reference. A cut from land-cover conversion is also accounted for
as harvest, but is reported with land-use-change fluxes rather than ordinary
wood-harvest fluxes; see [Land cover and land-use change](land-use-change.md).

## Treeinit removal

At the end of its prescribed establishment year, treeinit removes the
vegetation already present in affected patches before placing the prescribed
cohorts. For existing **trees**, it cuts the whole represented cohort, exports
all eligible stem wood and removable twig and coarse-root material, and removes
fine roots. Some leaves still become litter. This event does **not** put
removed material into the slow harvested-product pool, even when that pool is
enabled for ordinary harvests. Existing **grasses** are removed through
mortality accounting, which adds their biomass to litter rather than treating
them as harvested wood.

Treeinit's removal settings are fixed by that operation: changing the
`management2` harvest-efficiency or residue settings changes later cuts, not
this initial replacement. The added trees and their prescribed carbon and
nitrogen are explained in [Prescribed tree establishment](../howto/prescribed-tree-establishment.md).

## Reading outputs

Interpret a decrease in live biomass together with litter, harvested-product,
and harvest-flux outputs. A cut can transfer substantial carbon into litter
without exporting it; decomposition then affects later carbon fluxes. The
product pool can defer a portion of the export-associated carbon flux. If
treeinit is combined with `forest_init_poolsoil 1`, soil and litter pools from
the affected stands are pooled before the replacement, so compare the stand
and grid-cell totals as well as individual patch values.
