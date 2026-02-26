# Photosynthesis

This page describes the photosynthesis implementation used in DAVE simulations,
enabled by the `ifdavephotosyn` instruction file flag. The implementation
replaces the standard LPJ-GUESS photosynthesis with a coupled
photosynthesis–stomatal conductance scheme operating at subdaily time steps.

## Overview

The DAVE photosynthesis scheme has three major components:

1. **Leaf-level photosynthesis** — separate models for C3 (Farquhar–von
   Caemmerer–Berry) and C4 (von Caemmerer 2000) plants
2. **Stomatal conductance** — the Medlyn et al. (2011) optimal conductance
   model, coupled to the photosynthesis model through the intercellular CO₂
   concentration
3. **Canopy integration** — multi-layer Beer–Lambert scaling from leaf to
   canopy, with a within-canopy nitrogen gradient

Inputs at each subdaily timestep include air temperature, VPD, photosynthetically
active radiation (PAR), and CO₂ concentration. Outputs include gross
assimilation, leaf respiration, transpiration, and net primary productivity.

---

## Leaf Nitrogen and Photosynthetic Capacity

Maximum carboxylation rate ($V_{cmax}$) and maximum electron transport rate
($J_{max}$) are derived from leaf nitrogen content, which varies with canopy
depth.

### Top-of-canopy nitrogen

The nitrogen mass per unit leaf area at the top of the canopy,
$N_{top}$ (g N m⁻² leaf), is calculated from the individual's total leaf
nitrogen mass divided by projected leaf area. The minimum value is bounded
by the PFT parameter `ncleaf_min`.

### Within-canopy nitrogen gradient

Nitrogen decreases exponentially through the canopy following Beer–Lambert
extinction:

$$
N(L) = N_{top} \cdot e^{-k_N \cdot L}
$$

where $L$ is the cumulative leaf area index from the top (m² m⁻²) and $k_N$
is the nitrogen extinction coefficient set by the `extkn` instruction file
parameter.

### Vcmax and Jmax from leaf nitrogen (C3)

For C3 plants, $V_{cmax,25}$ and $J_{max,25}$ (capacities at 25 °C) are
estimated from $N_{top}$ using the empirical relationships of Ellsworth et
al. (2022):

$$
V_{cmax,25} = a_V \cdot N_{top} + b_V
$$

$$
J_{max,25} = a_J \cdot N_{top} + b_J
$$

where the intercept and slope coefficients come from a meta-analysis of
global plant data. Alternatively, if the PFT parameter `ifnlim_photosyn` is
disabled, fixed $V_{cmax,25}$ and $J_{max,25}$ values from the PFT
definition are used.

### Vcmax and Jmax from leaf nitrogen (C4)

For C4 plants the relationship between leaf N and photosynthetic capacity
follows the Yuzhen Fan dataset. The ratio of Rubisco carboxylation rate to
maximum PEP carboxylation rate ($V_{cmax}/V_{pmax}$) is fixed at
$r_{VC4} = 0.56$. The maximum electron transport rate is scaled relative
to $V_{cmax,25}$ via $r_{JC4} = 6.2$.

---

## Temperature Responses

### Modified Arrhenius equation (C3 Vcmax, Jmax; C4 Vpmax, Vcmax, Jmax)

Temperature responses of enzyme-limited capacities follow a peaked modified
Arrhenius equation (Johnson et al. 1942):

$$
k(T_l) = k_{25} \cdot \exp\!\left(\frac{H_a (T_l - T_{ref})}{T_{ref} \, R \, T_l}\right)
          \cdot \frac{1 + \exp\!\left(\dfrac{T_{ref} \, \Delta S - H_d}{T_{ref} \, R}\right)}
                     {1 + \exp\!\left(\dfrac{T_l \, \Delta S - H_d}{T_l \, R}\right)}
$$

where $T_l$ is leaf temperature (K), $T_{ref} = 298.15$ K, $R$ is the
universal gas constant, $H_a$ is the activation energy, $H_d = 200$
kJ mol⁻¹ is the deactivation energy, and $\Delta S$ is the entropy term.
Parameters for each enzyme are listed in `canexch.h`.

C3 activation energies and Michaelis–Menten constants follow Bernacchi et
al. (2001); C4 Vpmax, Vcmax, and Jmax use the NAD-ME parameters from
Sonawane et al. (2017) and Yamori et al. (2010).

### Optional thermal acclimation (C3)

When `iftempacclimation` is enabled, the optimum temperature for $V_{cmax}$
and $J_{max}$ shifts with the climatological mean of the warmest month
($T_{home}$, read from `file_thome`) following Kumarathunge et al. (2019).
This moves the peaked response to track long-term thermal environment.

### Leaf respiration (C3 and C4)

Leaf mitochondrial respiration ($R_d$) uses the non-Arrhenius temperature
response of Heskel et al. (2016):

$$
R_d(T_l) = R_{d,25} \cdot \exp\!\left(a \, (T_l - 25) + b \, (T_l^2 - 25^2)\right)
$$

where $a$ and $b$ are empirically fitted coefficients from the global
dataset in Heskel et al. (2016).

The dark respiration rate $R_{d,25}$ is calculated from $V_{cmax,25}$ using
the fraction $b_{C3} = 0.015$ (C3) or $b_{C4} = 0.025$ (C4).

### Light inhibition of leaf respiration

Leaf respiration is reduced by 30 % during daylight hours when
$\text{PAR} > 2\ \text{W m}^{-2}$, following Huntingford et al. (2017).
This prevents overestimation of daytime carbon loss.

---

## C3 Photosynthesis — Farquhar–von Caemmerer–Berry Model

The C3 assimilation rate is the minimum of three co-limiting rates
(Farquhar et al. 1980):

$$
A_n = \min(A_{nc},\; A_{nj},\; A_{np})
$$

### Rubisco-limited rate

$$
A_{nc} = V_{cmax} \frac{c_i - \Gamma^*}{c_i + K_c \left(1 + \dfrac{O_i}{K_O}\right)} - R_d
$$

where $c_i$ is intercellular CO₂ (mol mol⁻¹), $\Gamma^*$ is the
photorespiratory compensation point, $K_c$ and $K_O$ are the Michaelis–Menten
constants for CO₂ and O₂, and $O_i = 0.2095$ mol mol⁻¹ is ambient O₂.
Constants at 25 °C are $K_{c,25} = 404.9 \times 10^{-6}$ mol mol⁻¹,
$K_{O,25} = 278.4 \times 10^{-3}$ mol mol⁻¹, and
$\Gamma^*_{25} = 42.75 \times 10^{-6}$ mol mol⁻¹ (Bernacchi et al. 2001).

### RuBP-regeneration-limited rate

$$
A_{nj} = \frac{J}{4} \frac{c_i - \Gamma^*}{c_i + 2\Gamma^*} - R_d
$$

where $J$ is the electron transport rate (see below).

### Triose phosphate utilisation (TPU) limitation

$$
A_{np} = 3\,V_{pmax} - R_d
$$

where $V_{pmax}$ is the TPU rate (a fixed fraction of $V_{cmax,25}$).

### Electron transport rate

The response of $J$ to absorbed PAR follows the non-rectangular hyperbola
of von Caemmerer (2013):

$$
J = \frac{\text{PAR} + J_{max} - \sqrt{(\text{PAR} + J_{max})^2 - 4\,\theta\,\text{PAR}\,J_{max}}}{2\theta}
$$

where $\theta = 0.7$ is the curvature parameter and PAR is converted from
W m⁻² to mol(quanta) m⁻² s⁻¹ using $\text{JTOMOL} = 4.6 \times 10^{-6}$.

---

## C4 Photosynthesis — von Caemmerer Two-Cell Model

The C4 model follows von Caemmerer (2000), partitioning photosynthesis
between mesophyll (M) and bundle-sheath (BS) cells. The three potential
limiting rates in the BS are:

- **Rubisco-limited:** $V_c$ controlled by the BS CO₂ concentration $c_s$
  and the C4-specific $K_c$, $K_O$ (Boyd et al. 2015)
- **PEP carboxylase-limited:** $V_p = \min(V_{pmax} \cdot c_m/(c_m + K_P),\; V_{pr})$
  where $c_m$ is mesophyll CO₂ and $K_P$ is the Michaelis constant for PEP
  carboxylase (Sonawane et al. 2017; DiMario & Cousins 2019)
- **Electron transport-limited:** constraints on $J$ are split between M and
  BS cells using the partition factor $x_j = 0.4$, with a fraction of
  cyclic electron transport $f_{cyc} = 0.3$ (von Caemmerer 2021)

Leakage of CO₂ from BS back to M cells is parameterised by the bundle-sheath
conductance $g_{bs} = 0.003$ mol m⁻² s⁻¹ (von Caemmerer 2000).

The net assimilation $A_n$, stomatal conductance $g_s$, and intercellular
CO₂ concentration $c_i$ are solved simultaneously via an iterative loop
(maximum 20 iterations) until convergence.

An approximate formulation for C4 (`ifapproximate_C4`) is available that
simplifies the BS CO₂ balance, closely following von Caemmerer (2000,
Chapter 4).

---

## Stomatal Conductance

Stomatal conductance is modelled using the optimal conductance model of
Medlyn et al. (2011):

$$
g_s = g_0 + \left(1 + \frac{g_1}{\sqrt{D}}\right)\frac{A_n}{c_a}
$$

where $D$ is leaf-to-air vapour pressure deficit (kPa), $c_a$ is ambient
CO₂ concentration (mol mol⁻¹), $g_0$ is the minimum (cuticular)
conductance, and $g_1$ is the slope parameter (PFT-specific).

The conductance to CO₂ is related to water vapour conductance via
$D_{w}/D_{CO_2} = 1.6$ (ratio of molecular diffusivities, `DWDC`).

### Hydraulic downregulation (optional)

When `ifdavehydraulics` is enabled, stomatal conductance is additionally
limited by pre-dawn leaf water potential following Zhou et al. (2014). The
effective $g_1$ is reduced as a function of soil water status, and minimum
conductance is clamped to zero under severe drought.

### Coupled photosynthesis–stomatal conductance solver

For C3 plants, the demand curve (photosynthesis as a function of $c_i$) and
the supply curve (stomatal conductance linking $c_i$ to $c_a$ and $A_n$) are
coupled analytically. Substituting the Medlyn form into the Fick's law
relationship

$$
A_n = g_s \cdot (c_a - c_i)
$$

yields a quadratic in $c_i$ that is solved analytically following
Wang & Leuning (1998) (`photosyn_solver()`).

---

## Canopy Integration

Photosynthesis is calculated for each canopy layer and subdaily timestep
and then summed to give daily canopy-level fluxes.

### Canopy layers

The canopy is divided into discrete layers. PAR at the top of each layer is
attenuated following Beer–Lambert:

$$
I(L) = I_0 \cdot e^{-k \cdot L}
$$

where $k$ is the light extinction coefficient. Leaf nitrogen also follows
the same Beer–Lambert profile (with extinction coefficient `extkn`), so that
photosynthetic capacity declines in parallel with light availability.

### Vegetation structure

Within a patch, each individual (cohort) competes for light. The scheme
accounts for the position of each cohort within the canopy based on its
projected leaf area and height-structured competition.

### Transpiration

After photosynthesis, transpiration is calculated from the summed canopy
stomatal conductance using the Monteith approach (`aet_monteith`). Actual
evapotranspiration is used to update the soil water balance.

---

## Summary of Model Parameters

Key parameters controlling DAVE photosynthesis behaviour:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ifdavephotosyn` | — | Enable DAVE photosynthesis |
| `ifdavehydraulics` | — | Enable hydraulic stomatal downregulation |
| `iftempacclimation` | — | Enable Vcmax/Jmax thermal acclimation |
| `ifnlim_photosyn` | — | Derive Vcmax/Jmax from leaf nitrogen |
| `ifapproximate_C4` | — | Use approximate C4 photosynthesis formulation |
| `extkn` | — | Within-canopy N extinction coefficient |
| `g1` | PFT-specific | Medlyn stomatal slope parameter |

---

## References

Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis, A.R. and Long, S.P.
(2001). Improved temperature response functions for models of
Rubisco-limited photosynthesis. *Plant, Cell and Environment*, **24**,
253–259.

Bernacchi, C.J., Pimentel, C. and Long, S.P. (2003). In vivo temperature
response functions of parameters required to model RuBP-limited
photosynthesis. *Plant, Cell and Environment*, **26**, 1419–1430.

Boyd, R.A., Gandin, A. and Cousins, A.B. (2015). Temperature responses of C4
photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate
carboxylase, and carbonic anhydrase in Setaria viridis. *Plant Physiology*,
**169**, 1850–1861.

DiMario, R.J. and Cousins, A.B. (2019). Temperature response of
phosphoenolpyruvate carboxylase in C4 grasses. *Plant Science*, **281**,
23–28.

Ellsworth, D.S., Crous, K.Y., De Kauwe, M.G., Verryckt, L.T., Goll, D.,
Zaehle, S., et al. (2022). Convergence in phosphorus constraints to
photosynthesis in forests around the world. *Nature Communications*, **13**,
5005.

Farquhar, G.D., von Caemmerer, S. and Berry, J.A. (1980). A biochemical
model of photosynthetic CO₂ assimilation in leaves of C3 species. *Planta*,
**149**, 78–90.

Heskel, M.A., O'Sullivan, O.S., Reich, P.B., Tjoelker, M.G., Weerasinghe,
L.K., Penillard, A., et al. (2016). Convergence in the temperature response
of leaf respiration across biomes and plant functional types. *Proceedings
of the National Academy of Sciences*, **113**, 3832–3837.

Huntingford, C., Atkin, O.K., Martinez-de la Torre, A., Mercado, L.M.,
Heskel, M.A., Harper, A.B., et al. (2017). Implications of improved
representations of plant respiration in a changing climate. *Nature
Communications*, **8**, 1602.

Johnson, F.H., Eyring, H. and Williams, R.W. (1942). The nature of enzyme
inhibitions in bacterial luminescence: sulfanilamide, urethane, temperature
and pressure. *Journal of Cellular and Comparative Physiology*, **20**,
247–268.

Kumarathunge, D.P., Medlyn, B.E., Drake, J.E., Tjoelker, M.G., Aspinwall,
M.J., Battaglia, M., et al. (2019). Acclimation and adaptation components
of the temperature dependence of plant photosynthesis at the global scale.
*New Phytologist*, **222**, 768–784.

Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Prentice, I.C.,
Barton, C.V.M., et al. (2011). Reconciling the optimal and empirical
approaches to modelling stomatal conductance. *Global Change Biology*,
**17**, 2134–2144.

Sonawane, B.V., Sharwood, R.E., von Caemmerer, S., Whitney, S.M. and
Cousins, A.B. (2017). Short-term thermal photosynthetic responses of C4
grasses are independent of the biochemical subtype. *Journal of
Experimental Botany*, **68**, 5583–5597.

von Caemmerer, S. (2000). *Biochemical Models of Leaf Photosynthesis*.
CSIRO Publishing, Collingwood, Australia.

von Caemmerer, S. (2013). Steady-state models of photosynthesis. *Plant,
Cell and Environment*, **36**, 1617–1630.

von Caemmerer, S. (2021). Updating the steady-state model of C4
photosynthesis. *Journal of Experimental Botany*, **72**, 6003–6017.

Wang, Y.-P. and Leuning, R. (1998). A two-leaf model for canopy conductance,
photosynthesis and partitioning of available energy I: Model description and
comparison with a multi-layered model. *Agricultural and Forest
Meteorology*, **91**, 89–111.

Yamori, W., Noguchi, K., Hanba, Y.T. and Terashima, I. (2010). Effects of
internal conductance on the temperature dependence of the photosynthetic
rate in spinach leaves from contrasting growth temperatures. *Plant and Cell
Physiology*, **47**, 1069–1080.

Yamori, W., Hikosaka, K. and Way, D.A. (2014). Temperature response of
photosynthesis in C3, C4, and CAM plants: temperature acclimation and
temperature adaptation. *Photosynthesis Research*, **119**, 101–117.

Zhou, S., Medlyn, B.E., Sabaté, S., Sperlich, D. and Prentice, I.C.
(2014). Short-term water stress impacts on stomatal, mesophyll and
biochemical limitations to photosynthesis differ consistently among tree
species from contrasting climates. *Tree Physiology*, **34**, 1035–1046.
