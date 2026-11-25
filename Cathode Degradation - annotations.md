# Cathode Degradation of PEMFC — Annotation Extract

Extracted from the highlighted copy `Cathode Degradation (noted).pdf` in `references/`.
Source paper: Fensterle et al., *Spatially resolved electrochemical characterization of a PEM fuel cell
stack during durability testing*, **Int. J. Hydrogen Energy 250 (2026) 156060**.

The highlights were flattened into the PDF (no annotation layer left), so they were recovered by matching
the colored highlight regions to the text underneath. Line-break hyphenation has been repaired for
readability; wording is otherwise as highlighted. The three colors used in the PDF are noted where it helps
(🟡 yellow = context/definitions, 🟢 green = method/mechanism, 🟣 magenta = quantitative findings), but the
grouping below is by **topic**, as requested.

---

## 1. Lab setup

- The study targets the **durability of proton exchange membrane fuel cells (PEMFC)** for heavy-duty
  applications, and the need to **identify and understand the aging mechanisms**. 🟡
- Hardware is from the **RealHyFC project**: the generic **HyFaB stack** with an **active area of 280 cm²**
  per cell and a **project-specific MEA**. Two builds were used — a **10-cell stack for the durability test**,
  and an **8-cell stack for the BoL EIS sensitivity analysis**. 🟡
- Investigating the **membrane electrode assembly (MEA)** is of particular importance, since the MEA drives
  fuel-cell performance degradation during aging. 🟡
- **Coolant flow was kept constant**, chosen so the **inlet-to-outlet temperature difference was 10 K at the
  highest load point (1.8 A/cm²)**. 🟡
- Spatially resolved EIS was done with the **CDDM (segmented current-density measurement board) and special
  expansion boards** that transmit the current from the measurement segments to the **Zahner device**. 🟡
- Focus of the work is on **spatially resolved aging effects across the active area**. 🟡

## 2. Testing procedure

- A **heavy-duty-relevant load cycle was applied for more than 2000 h** to the 280 cm² stack. 🟡
- At regular intervals the stack was characterized with **polarization curves, sensitivity tests, and
  electrochemical impedance spectroscopy (EIS) at cell level and spatially resolved**. 🟢
- Full characterization sequence (after break-in, which aims to reach full cell performance): **six
  polarization curves, the partial sensitivity test, EIS at cell level and spatially resolved, and cyclic
  voltammetry**. 🟢
- **Polarization curves:** 16 load points from **0 to 1.8 A/cm²**, measured in both **descending and ascending**
  order. The minimum current for stoichiometric operation was **0.3 A/cm²** for the durability test and
  **0.5 A/cm²** for the full sensitivity test. 🟡
- **Sensitivity tests** study the influence of operating conditions on performance — **temperature, humidity,
  pressure, and stoichiometry, at four different current densities**. 🟢
- **EIS** is a powerful method for separating the resistances of **transport, conduction, and reaction
  processes, which occur at different frequencies**. 🟡
- Both cell and spatially resolved impedance were measured **galvanostatically**. For the cell impedance the
  **AC was 5 % of the DC**; for the spatially resolved measurements a **slightly higher AC gave less noise**. 🟢/🟣

## 3. Results

**Overall performance loss**
- **Cell voltage decreases over time**, and the **largest jump is between the BoL curve and the curve after
  300 h (25 cycles)** — about half the total loss happens in the first 300 h. 🟣
- Degradation rate **15.7 µV/h** at 1.8 A/cm² overall; excluding the reversible loss (comparing the
  after-air/air-start curves) it drops to **13.0 µV/h** — i.e. ~2.7 µV/h of the loss is reversible. 🟣

**Effect of operating conditions (sensitivity)**
- **Degradation impact increases with decreasing pressure**, especially in the **resistance of the oxygen
  reduction reaction (ORR)**. 🟣
- Temperature, humidity and stoichiometry variations did **not** show a strong influence of aging; the aging
  effect on sensitivity **becomes evident under varying pressure**. 🟣
- **Increasing operating pressure enhances reactant availability** and cell operation (as at BoL), but with a
  **much more beneficial effect after aging**: higher partial pressure **improves oxygen transport to the
  catalyst layer**, mitigating the impact of aging and supporting more uniform, robust performance. 🟣
- Conclusions: **aging effects are strongest at low pressures and elevated temperatures**. Heterogeneity of
  the current-density distribution **decreases with increasing humidity, pressure and stoichiometry**, but
  **increases with temperature**. **Mass-transport resistance is most strongly affected by stoichiometry**,
  with pronounced increases in the **outlet region at low stoichiometric ratios**. 🟢

**Current-density distribution (spatially resolved)**
- The active area is divided into **27 columns in the flow direction** and averaged per column. 🟢
- Over time there is a **visible shift of the current-density maximum toward the cathode inlet**; the
  **inhomogeneity increases** and the **maximum reaches higher values after 2300 h**. 🟢
- At the lowest load, local current densities at the **air inlet become even lower**, compensated by higher
  values in the middle; at higher loads no air-inlet effect, but an **increased slope between the maximum
  (near the air inlet) and the middle–outlet area**. 🟢
- Degradation therefore has a **greater influence on the current-density distribution at lower temperatures**
  (and at low stoichiometry). 🟢

**EIS / resistance breakdown (at 1.8 A/cm²)**
- The paper focuses on EIS at **1.8 A/cm²**, where **mass-transport resistance** matters most for both aging
  and sensitivity. 🟢
- **HFR (high-frequency resistance) = the ohmic losses** — mainly the **membrane resistance**, plus electronic
  contact resistances and electrode-related effects (ionomer distribution, possible swelling). 🟡
- **LFR (low-frequency resistance) = the sum of all losses.** The arcs between the HFR, IFR and LFR intercepts
  correspond to two processes: **mass transport (low frequency)** and **ORR charge transfer on the cathode
  (intermediate frequency)**. 🟡
- The **membrane-resistance increase is attributed to membrane degradation**: **chemical aging by radical
  attack reduces proton conductivity and raises membrane resistance**; rising **interfacial contact
  resistance** (e.g. GDL–bipolar plate) may also contribute. 🟢
- Mass-transport resistance would **normally increase** over time (blocked gas-diffusion paths from liquid
  water accumulation, loss of GDL hydrophobicity, pore-structure changes, carbon-support corrosion) — but here
  it stayed roughly constant. 🟢
- Catalyst degradation primarily affects the **electrochemically active surface area (ECSA)**, and therefore
  the **charge-transfer and mass-transport resistances**, while the **ohmic (membrane) resistance remains
  largely unaffected**. 🟢

---

## 4. The membrane part — what it's actually saying

*(This section is an explanation, not from your highlights — it unpacks the membrane discussion you flagged.)*

**First, a naming note.** This is a **fuel-cell stack, not a battery**. A battery stores energy; this PEMFC
*converts* it, producing electricity continuously from H₂ + O₂ → water for as long as it's fed. The paper even
makes this contrast — PEMFCs are attractive for heavy-duty transport because, compared with battery-electric
systems, they give more range and faster refueling. So "heavy-duty battery stack" is really a heavy-duty
**fuel-cell** stack.

**What the membrane is.** The membrane is the **proton exchange membrane** (a thin polymer film, Nafion-type),
and it is the heart of each cell's MEA (membrane electrode assembly = membrane + the two catalyst layers). Its
job is threefold: conduct **protons (H⁺)** from the anode to the cathode, **block electrons** so they're forced
through the external circuit (that electron flow *is* your usable current), and **keep the H₂ and O₂ gases
physically separated**. Proton conduction only works when the membrane is **well hydrated**, which is why
humidity and temperature matter so much.

**Why the paper calls membrane degradation "the most critical."** The distinction the intro is drawing is
*failure vs. fading*:
- Most components (catalyst, gas-diffusion layer) degrade **gradually** — you slowly lose voltage.
- The membrane can fail **catastrophically**. If it develops **cracks or pinholes**, H₂ and O₂ leak straight
  across (**gas crossover**), which is a safety problem and can **end the cell's life outright**, not just
  lower its performance.

So "most critical" means *highest-stakes*, not *most likely* — it's the one failure mode that turns a slope
into a cliff.

**How the membrane degrades (the mechanisms the paper lists).**
- **Mechanical.** The membrane swells when wet and shrinks when dry, so repeated **hydration/dehydration
  cycles** (especially during start-up and shutdown) build up mechanical stress and can open cracks or
  pinholes. Defects can also originate during manufacturing or assembly.
- **Chemical.** **Chemical decomposition / radical attack**: reactive radicals (largely from hydrogen peroxide
  formed by gas crossover) attack the polymer backbone, causing **membrane thinning** and loss of the acidic
  sites that carry protons. This lowers **proton conductivity** and shows up as a higher ohmic resistance.
- **Delamination.** The **catalyst layer detaching from the membrane** opens cavities where water collects,
  raising charge-transfer resistance.

The first route mostly causes *failure*; the second and third mostly cause *performance loss*.

**How they actually measured membrane health — this is the part on page 8.** In EIS, the membrane shows up as
the **high-frequency resistance (HFR)**, i.e. the ohmic resistance **R₀** — the leftmost point of the Nyquist
plot, where all the curves cross the real axis. Ohmic loss is dominated by the *ionic* part, which is proton
transport through the membrane. So **as the membrane chemically ages and loses proton conductivity, R₀ / HFR
goes up** — that rise is the measurable fingerprint of membrane degradation.

**What this particular test found (and the likely source of your confusion).** The intro builds the membrane
up as the most critical risk, but the **results show the membrane held up well** over these 2300 h:
- The **membrane resistance increased only slightly**, with a modest jump after **1150 h** — no cracks,
  pinholes, or catastrophic failure.
- At cell level, the aging test caused **only minor effects** on the average impedance.
- The degradation that *did* dominate was on the **cathode side** — the **ORR charge-transfer resistance**
  rising from **catalyst/ECSA loss**, concentrated in the **cathode outlet region** — while the **ohmic
  (membrane) resistance stayed largely unaffected**.

In short: the membrane is framed as the component you most need to protect *because* its failure is the worst
case, but in this study it was the **cathode catalyst**, not the membrane, that degraded most. That's the whole
reason the project (and the dataset) is called **cathode** degradation.
