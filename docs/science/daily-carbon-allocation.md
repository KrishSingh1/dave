## Daily Carbon Allocation

C allocation is inspired by the method proposed by Grote (REF),. In this
approach, daily C sink, which is constrained by Net Primary Productivity
and N availability, is allocated to four tree C pools (i.e., leaves,
roots, stem, storage) based on their sink strength ($\phi_{p}$) , which
is determined by the difference between the actual C of each of the
pools and its potential target:

  --------------------------------------------------------------------------------------------------
$$\phi_{p} = \max\left( 0,\ \frac{C_{p.target(t)} - C_{p(t)}}{C_{p.target(t)}} \right) \tag{1} $$
  ------------------------------------------------------------------------------------------ -------

  --------------------------------------------------------------------------------------------------

where $C_{p}$ and $C_{p.target}$ are the daily actual and target carbon
mass (kgC m^-2^) of pool $p$ at time $t$. The proportion of C that is
allocated to each of the pools ($\alpha_{p}$) is calculated by the
proportion of the sink strength, relative to that of the other pools:

  ------------------------------------------------------------------------
$$\alpha_{p} = \frac{\phi_{p}}{\sum_{p = 0}^{n_{p}}\phi_{p}} \tag{2} $$
  ---------------------------------------------------------------- -------

  ------------------------------------------------------------------------

where $n_{p}$number of carbon pools. Sink strengths are therefore key to
the allocation process, as these are affected by the C targets (which
are calculated differently for each of the pools) and to the turnover
rates.

Daily change in C pools is therefore calculated using

  ---------------------------------------------------------------------------------
$$\frac{dC_{p}}{dt} = {\overline{\alpha}}_{p}C_{grow} - C_{p}\tau_{p} \tag{3} $$
  ------------------------------------------------------------------------- -------

  ---------------------------------------------------------------------------------

where , $\overline{\alpha_{p}}$ is a 30-day moving average of daily
$\alpha_{p}$ that is calculated, which intends to tune down further
erratic daily fluctuations in C allocation for each of the pools;
$C_{grow}$ is daily NPP; mass of pool p at day t-1 after x% was taken
for allocation to reproductive pool, which is a PFT parameter which is
*zero* below the age of reproduction maturity, assumed to be 10 years in
Eucalyptus; and $\tau_{p}$ is the daily turnover rates (1/day) of pool
p.

Potential C targets for the different pools are calculated to fulfil the
assumptions that trees would opt to: (i) grow at the minimum to pay back
the daily turnover of sapwood; (ii) maintain a constant ratio of leaf
area to sapwood area (constant k_latosa); (iii) alllocate to
below-ground pools (i.e., fine roots) depeinding on the availability (or
deficit) of stored water and Nitrogen in the soil.

Sapwood targe is calculated as:

  -------------------------------------------------------------------------------
$$C_{Sap.target} = C_{Sap}\left( 1 + \tau_{sap} \right)^{n.outlook} \tag{4} $$
  ----------------------------------------------------------------------- -------

  -------------------------------------------------------------------------------

where $n.outlook$ is the number of days in which turnover needs to be
replenished by the daily allocation to sapwood. $n.outlook$ is used to
prevent erratic daily fluctuations in $\alpha$, and result with a more
realistic and smoot allocation pattern. Here we set $n.outlook = 365$.

Potential C limit for leaves is calculated using the following equation,
fulfilling assumption *ii*, which states that the tree will opt to
maintain the ratio between tree leaf area to sapwood area, set by the
ratio between leaf and sapwood area, assuming the pipe model.

  ---------------------------------------------------------------------------------------------------------------------------
$$C_{Leaf.target} = \left( C_{Leaf},\ {\ \ \ C}_{Sap.target}\left( \frac{{\ k}_{LA:SA}}{SLA\ h} \right) \right) \tag{5} $$
  ------------------------------------------------------------------------------------------------------------------- -------

  ---------------------------------------------------------------------------------------------------------------------------

where SLA is the specific leaf area (m^2^ kgC^-1^); and ${\ k}_{LA:SA}$
is the ratio between leaf and sapwood area (Table X), both are PFT
parameters. Potential root target is calculated using:

  ------------------------------------------------------------------------------------------
$$C_{root.target} = max\left( C_{Root},\ \ \frac{C_{Leaf.target}}{l:r} \right) \tag{6} $$
  ---------------------------------------------------------------------------------- -------

  ------------------------------------------------------------------------------------------

where $k_{l:r}$ is the ratio between leaf to root, which decreases when
an increasing soil water and N stress, which is calculated using:

  -----------------------------------------------------------------------------------------------------------
$$k_{l:r} = k_{l:r.\ max}\ \ max\left( \min\left( s_{w},\ s_{N} \right)\ ,\ k_{l:r.min} \right) \tag{7} $$
  --------------------------------------------------------------------------------------------------- -------

  -----------------------------------------------------------------------------------------------------------

where $k_{l:r.\ max}$ and $k_{l:r.min}$ are the maximum and minimum leaf
to root ratio; and $s_{w}$ and $s_{N}$ are scalars that represent soil
moisture and nitrogen deficit, reactively.

Maximum capacity of a cohort to store none-structural carbohydrates (kgC
m^-2^) is assumed to depend on the capacity of the leaves, sapwood and
coarse roots (e.g., lignotuber) using:

  ----------------------------------------------------------------------------------------------------------------------------------------------------
$$C_{storage.max} = \ C_{Leaf}k_{Leaf.NSC.max} + \ \left( {p_{AG}C}_{Sap}k_{AG.Sap.NSC.max} + \ {P_{BG}C}_{Sap}k_{BG.Sap.NSC.max} \right)$$    (8a)
  --------------------------------------------------------------------------------------------------------------------------------------------- ------
$$C_{storage.target} = \ k_{NSC.max}\ C_{storage.max}$$                                                                                        (8b)

  ----------------------------------------------------------------------------------------------------------------------------------------------------

where $k_{Leaf.NSC.max}$ and $k_{Sap.NSC.max}$ (kgC m^-2^) are PFT
parameters that define the maximum possible none-structural carbohydrate
mass that the leaves and sapwood can store; and $p_{AG}$ and $p_{BG}$
are the proportions of the sapwood that is above and below ground,
respectively. These are calculated from the PFT-specific root to shoot
ratio, $k_{r:s}$ using

  -----------------------------------------------------------------------
$$p_{AG} = \frac{1}{\left( k_{r:s} + 1 \right)}$$                 (9a)
  ---------------------------------------------------------------- ------
$$p_{BG} = 1 - p_{AG}$$                                           (9b)

  -----------------------------------------------------------------------

To avoid that NSC storage exceeding its capacity, maximum storage
($C_{storage.max}$) is scaled down using the a fixed proportion
$p_{NSC.max}$.

Following fire, large proportion of Australian tree species resprout
epicormically or basally through buds located under the bark of the stem
and branches and on the lignotuber (respectively). Other tree species,
however, lack the ability to resprout, and rely on the canopy seedbank
to cause a major recruitment event following fire (REF). We therefore
assume that in order to maintain these two different strategies,
resprouting tree species would: i) prioritise carbon allocation to NSC,
when possible; ii) will have larger storage capacity per unit mass of
sapwood, leaves and coarse roots.

## Allometry

Change in tree height and diameter is calculated from the daily sapwood
carbon increment, $C_{grow.Sap}$, by taking a partial derivative of the
change

## Carbon starvation mortality
