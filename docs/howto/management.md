# Management types for natural and forest stands

A **stand type** describes an area of a grid cell, its land cover, and the
management regimes it may use. A **management type** (`mt`) is a named set of
rules for that stand: which plant functional types (PFTs) may be planted or
establish, when and how trees are cut, and whether disturbance is suppressed.
Each stand contains one or more replicate patches. Patch ages and vegetation
can differ even when their stand type and management type are the same.

Management rules become relevant when land-cover mode (`run_landcover 1`) is
used. Normal forest cutting is evaluated for natural and forest stands at the
start of each year. Treeinit can use these same regimes to specify what
happens **after** its one-time replacement. [Harvesting and vegetation
removal](../science/harvesting.md) explains where biomass goes when a cut
occurs; [land-cover change](land-cover.md) explains how stand area changes.

## Assign a regime

Define an `mt` block and assign it to a stand type with `management1`.
`management2` and `management3` are optional later regimes. In a treeinit
stand, `management2` is required: treeinit switches affected stands to it
after installing the prescribed trees. The second regime can allow subsequent
growth without scheduled cutting, or prescribe later thinning and cutting.
Treeinit also adds the PFTs named in its input file to the second regime's
selection. Set `cutfirstyear 0` if the new regime should preserve the trees
through its first year.

```ins
mt "Before_initialisation" (
    harvest_system "continuous"
    cutfirstyear 0
)

mt "After_initialisation" (
    planting_system "selection"
    selection "MRS"
    reestab "restricted"
    harvest_system "continuous"
    cutfirstyear 0
)

st "Forest" (
    common_stand
    stinclude 1
    landcover "forest"
    naturalveg "all"
    management1 "Before_initialisation"
    management2 "After_initialisation"
    file_treeinit_st "forest-treeinit.txt"
)
```

An alternative is to place management settings directly in `st`; these apply
to its first regime. Prefer named `mt` blocks when using several regimes, and
avoid mixing inline first-regime settings with a named `management1` because
the named regime takes precedence. `mtinclude 0` excludes a named management
type from the run. Stand-type keys `mt1_startyear`, `mt2_startyear`, and
`mt3_startyear` set calendar years for normal regime rotation;
`rotation_wait_for_cc 1` waits for a clearcut before switching under a
clearcut system. Treeinit makes its own switch to regime 2 at the event.

## Core choices

`planting_system` determines which tree PFTs are planted after cutting. An
empty value leaves the stand's eligible tree PFTs available. `"selection"`
uses the names in `selection`; `"monoculture"` uses `pft`. Other supported
values select broad groups such as evergreen needleleaf or deciduous
broadleaf trees. `reestab "all"` permits eligible PFTs to establish later;
`"restricted"` limits tree reestablishment to planted PFTs, and `"none"`
prevents tree reestablishment. The stand type's `landcover` and `naturalveg`
settings still define the outer set of eligible vegetation. See
[Establishment](../science/establishment.md) for climatic establishment limits.

`harvest_system "clearcut"` uses rotations and optional intermediate
thinnings. `"continuous"` schedules repeated partial cutting and may use a
second cutting period. `"probability_input"` uses externally supplied
cutting probabilities and intensities. An empty harvest system has no normal
scheduled wood cutting; other events, including land-cover conversion and
treeinit, may still remove vegetation. Parameters below are activated by the
system or feature they describe; specifying a thinning strength alone does
not schedule a cut.

## Management-type parameter reference

Values in the same row are related but remain separate input parameters.
Years described as *calendar years* use the model's calendar, while intervals
and patch ages count elapsed years. Fractions generally range from 0 to 1.
Unless stated otherwise, switches use `0` for off and `1` for on.

### Selection, establishment, and grazing

| Parameter | Meaning |
| --- | --- |
| `planting_system` | Planting rule: empty, `selection`, `monoculture`, or a supported tree group. |
| `selection` | Space-separated PFT names used for selected planting and preferential cutting. |
| `pft` | Single PFT for `monoculture` planting; also used by crop regimes. |
| `plantdensity`, `plantdensity_pft` | Seedlings per hectare for selected PFTs, or the single `pft`, respectively. |
| `set_planting_density` | Use specified tree planting densities after a clearcut and bypass environmental establishment limits for that planting. |
| `reestab` | `all`, `restricted`, or `none`: controls later tree establishment after planting. |
| `relaxed_establishment` | Ignore climatic limits on establishment once management starts. |
| `stochmort`, `stochestab` | Enable stochastic mortality or establishment, respectively. Both default to on. |
| `grazeintens` | Daily grazing fraction for managed grass vegetation. |

For managed grassland, `mowdates` (day-of-year values) and `mowintens`
(fraction removed) can also be set **inline in `st`**. They are not accepted
inside a named `mt` block.

### Start and rotation timing

| Parameter | Meaning |
| --- | --- |
| `firstmanageyear` | First calendar year of management behavior, including relevant establishment and disturbance rules. |
| `firstcutyear` | Earliest calendar year of normal wood cutting. |
| `firstclearcutyear` | Earliest calendar year for a clearcut. |
| `clearcut_first_manageyear` | Clearcut when management first begins on an existing stand. |
| `firstcutyear_is_referenceyear` | Use years since `firstcutyear`, rather than patch age, to time thinning. |
| `cutfirstyear` | On entry to a new management regime: `0` preserve trees, `1` cut without cloning, `2` cut with cloning. The default is `1`; set this explicitly for treeinit stands. |
| `cutfirstyear_unsel` | Also cut PFTs outside the new selection on entry or cloning. |
| `cutinterval` | Main rotation or cutting interval in years. |
| `delayduecutting` | Spread overdue density-triggered clearcuts over this many years after `firstclearcutyear`. |
| `distribute_patch_ages` | Stagger patch ages in a newly managed forest stand. |
| `distribute_cuttings_among_patches` | Stagger cutting events among patches. |

### Cutting systems and thinnings

| Parameter | Meaning |
| --- | --- |
| `harvest_system` | `clearcut`, `continuous`, `probability_input`, or empty. External probability inputs require their own forcing files. |
| `thintime` | One or more thinning times as fractions of the first rotation or cutting period. |
| `thinstrength`, `thinstrength_unsel` | Fractions cut at corresponding `thintime` entries for selected and unselected PFTs. |
| `thinselectpft` | PFT preference at thinning: `0` none, `1` unselected, `2` selected, `3` separate selected/unselected strengths, `4` shrubs and shade-intolerant PFTs. |
| `thinselectage` | Age preference: `0` none, `1` younger, `2` older trees. |
| `thinselectdiam` | Diameter preference: `0` none, `1` smaller, `2` larger, `3` diameter-limit cutting. |
| `secondintervalstart` | Years after the first period begins at which continuous cutting enters its second period. |
| `secondcutinterval` | Cutting interval in that second period. |
| `secondthintime`, `secondthinstrength`, `secondthinstrength_unsel` | Second-period counterparts of the first-period thinning times and strengths. |
| `secondthinselectpft`, `secondthinselectage`, `secondthinselectdiam` | Second-period counterparts of the PFT, age, and diameter preferences. |
| `diam_cut_low`, `diam_cut_high` | Diameter thresholds in centimetres for diameter-limit cutting; between them partial cutting applies and above the high limit all eligible trees are cut. |
| `adapt_diam_limit` | Adjust the low diameter threshold for stands with small trees. |
| `ifthin_reineke` | Use automatic thinning based on Reineke's stand-density rule. |
| `alpha_st`, `rdi_target` | Self-thinning coefficient and target relative density for Reineke thinning; a lower density target means stronger thinning. |
| `ifclearcut_by_density`, `dens_target_cc` | Enable clearcut when tree density falls below the specified trees-per-hectare threshold. |
| `ifclearcut_optimal_age` | Use the model's optimal rotation age as a clearcut trigger. |
| `killgrass_at_cc` | Remove grass vegetation during clearcut. |

The `thin*` and `secondthin*` parameters accept corresponding arrays: the
first value in each array describes the first thinning, the second value the
second thinning, and so on. A relative thinning time of `0.5` means halfway
through the relevant period. Density and optimal-age triggers apply to
clearcut management; if several triggers are set, the model resolves their
priority internally, so avoid assuming all will cause independent cuts.

### Target PFT composition

Target cutting aims to move the **carbon biomass shares** of selected PFTs
towards prescribed fractions. It is separate from a rotation thinning.

| Parameter | Meaning |
| --- | --- |
| `targetfrac` | Target biomass fractions corresponding to PFTs in `selection`. |
| `file_targetfrac_pft_mt` | File supplying time-varying PFT target fractions. |
| `targetfrac_input_mode` | `0` use file fractions for cutting, retain instruction-file selection; `1` use file fractions and restrict selection to PFTs with nonzero fractions; `2` use file fractions only to define selection, without target cutting. |
| `firsttargetyear`, `lasttargetyear` | Calendar-year window for target cutting. |
| `targetstartage`, `targetcutinterval` | Minimum patch age and years between target cuts. |
| `targetcutmode` | `1` calculate deviations at patch level; `2` or `3` use stand-level deviations. |
| `targetthinselectage`, `targetthinselectdiam` | Prefer younger/older (`0`/`1`/`2`) or smaller/larger (`0`/`1`/`2`) individuals during target cutting. Diameter preference can override age preference. |
| `suppress_second_target` | Stop target cutting when continuous management enters its second period. |

### Disturbance and harvest accounting

| Parameter | Meaning |
| --- | --- |
| `suppress_fire` | Suppress prognostic fire while the regime is active and within its year window. |
| `suppress_fire_firstyear`, `suppress_fire_lastyear` | Calendar-year window for fire suppression. |
| `suppress_disturbance` | Suppress other prognostic disturbance while the regime is active and within its year window. |
| `suppress_disturbance_firstyear`, `suppress_disturbance_lastyear` | Calendar-year window for disturbance suppression. |
| `harv_eff_thin`, `harv_eff_cc` | Fractions of eligible stem wood exported at thinning and clearcut. |
| `res_outtake_twig_thin`, `res_outtake_twig_cc` | Fractions of twigs and attached leaves removed at thinning and clearcut. |
| `res_outtake_coarse_root_thin`, `res_outtake_coarse_root_cc` | Fractions of coarse roots and stumps removed at thinning and clearcut. |

The six harvest and residue values default to *unset*. Stem efficiency and
twig removal then use PFT-level settings; coarse-root removal defaults to
0.1. They affect normal management cuts, not the one-time
vegetation removal performed by treeinit. See [Harvesting and vegetation
removal](../science/harvesting.md) for the resulting biomass pathways.

### Crop-specific keys

The parser also accepts `hydrology`, `sdate`, `hdate`, `nfert`, `tillage_fact`,
and `fallow` in an `mt` block. They describe crop water supply, sowing and
harvest dates, nitrogen fertilisation, tillage, and fallow. They are listed
here for completeness but are outside the DAVE natural/forest workflow.
