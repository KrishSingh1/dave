# Disturbance

This page documents the disturbance pathways used by the DAVE fork of
LPJ-GUESS. The word "disturbance" is used for several related but distinct
mechanisms, so the important thing is to know which pathway is active in a
given run.

## Disturbance Concepts

### Generic Patch-Destroying Disturbance

The generic patch-destroying disturbance is implemented by
`disturbance(Patch&, double)` in `modules/vegdynam.cpp`.

It is controlled by `ifdisturb`, `disttype`, and either `distinterval` or
`file_disturb`, depending on the disturbance type. It applies in individual
and cohort modes.

When this pathway disturbs a patch, it:

- loops over every individual or cohort in the patch
- reports disturbance mortality fluxes with `MORT_DIST`
- optionally routes biomass through harvest-style accounting when
  `ifback_disturb_as_harv` is enabled
- kills and removes every vegetation object from the patch
- sets `patch.disturbed = true`
- resets `patch.age = 0`

This is the trunk-style disturbance mechanism where a disturbance destroys all
vegetation in the patch.

### Dave Deterministic Disturbance

Dave deterministic disturbance is a fixed-year wrapper around the generic
patch-destroying disturbance.

It is selected with the following combination of parameters:

```ins
ifdisturb 1
disttype "DETERMINISTIC"
file_disturb "..."
```

The `file_disturb` file contains one calendar year per line. In those years,
every patch is disturbed by calling `disturbance(patch, 1.0)` (ie disturbance is
triggered with probability 1).

This pathway is mainly useful for controlled debugging or reproducibility runs
where the desired behaviour is a known forest age at output time.

### Stochastic Generic Disturbance

Stochastic generic disturbance uses the same patch-destroying mechanism as the
Dave deterministic pathway. The difference is only the probability source:
stochastic disturbance calls `disturbance(patch, 1.0 / distinterval)`, using
the gridcell or stand disturbance interval.

### Prescribed Disturbance Inputs

Prescribed disturbance inputs are the disturbance machinery ported from
`european_applications`. They are activated by one or more of:

```ins
file_dist_clearcut
file_dist_fire
file_dist_wind
file_dist_barkbeetle
```

These inputs are handled by `prescribed_disturbance()` and
`force_prescribed_disturbance()` in `modules/vegdynam.cpp`.

This pathway is not the generic patch-destroying disturbance. The input is an
area fraction for the whole gridcell, not a list of particular patches to
disturb. For example, a prescribed wind input may say, in effect, "wind
disturbed 8% of the forest/natural area in this gridcell this year". The model
then has to decide which simulated stand area represents that disturbed 8%.

The prescribed-disturbance code does this as an allocation problem for each
disturbance type:

- start with the prescribed disturbed fraction
- convert it to an area over applicable stands
- rank the existing stands by a disturbance-specific susceptibility metric
- take stands from the high-susceptibility end of that ranking until the
  prescribed disturbed area has been allocated
- apply the disturbance treatment to the selected stand area

The stands are sorted because some stands may be more susceptible to a given
disturbance type than other stands, and if the disturbed area is less than the
stands' total area, only some stands are affected. The most susceptible stands
are disturbed first, so that, of all the stands, those that are affected by a
given disturbance type are those most susceptible to that type of disturbance.

Different disturbance types use different rankings because "most susceptible"
has different meanings for different processes:

- prescribed clearcut uses stand diameter structure, so larger/more harvestable
  stands are cut first
- prescribed fire uses fire sensitivity, so the stands which are burnt are those
  most susceptible to fire
- prescribed wind uses storm sensitivity
- prescribed bark-beetle disturbance uses a bark-beetle sensitivity metric based
  on tree size, density, and PFT weighting

These rankings do not make the total disturbed area larger or smaller. They
only decide which simulated stand area receives the prescribed disturbance.

If the remaining disturbance area is larger than a selected stand, the whole
stand is disturbed. If the remaining disturbance area is smaller than the next
selected stand, the code splits the stand area. It reduces the original stand's
area, clones a new stand to represent the disturbed fraction, and applies the
disturbance treatment to that cloned stand. This is why prescribed disturbance
uses land-cover/stand-transfer machinery even though the land-cover class may
not change: the model is carving out a fraction of a stand so that only the
prescribed amount of area is disturbed.

## Precedence

During active prescribed-disturbance input years, the generic patch-destroying
disturbance pathway is suppressed in `vegetation_dynamics()`.

This means that Dave deterministic disturbance and stochastic generic
disturbance do not run in years covered by active prescribed clearcut, fire,
wind, or bark-beetle disturbance inputs.

The intent is to avoid two independent disturbance systems acting on the same
patches in the same year.

## Practical Guidance

Use Dave deterministic disturbance when the desired behaviour is a controlled
whole-patch vegetation reset in known calendar years. This is really only useful
for single-patch debugging simulations.

Use stochastic generic disturbance when the desired behaviour is the original
LPJ-GUESS random patch-destroying disturbance.

Use prescribed disturbance inputs when the desired behaviour is to impose
external disturbed area fractions for clearcut, fire, wind, or bark beetle and
allocate those disturbances across stands.

Do not treat prescribed disturbance inputs as a drop-in replacement for Dave
deterministic disturbance. They differ in allocation, biomass routing,
treatment of non-tree vegetation, and patch-age handling.
