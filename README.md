# Cathode Degradation of PEMFC — Durability Dataset (RealHyFC stack)

## Overview

This dataset is a durability test of a 10-cell PEM fuel cell stack (metallic bipolar plates, 280 cm² active area per cell) run for
about 2300 h on a heavy-duty load cycle (RealHyFC cycle, based on PEMTASTIC; ~12 h per cycle, 200 cycles). At regular intervals
the stack was characterized with nominal polarization curves, partial sensitivity tests (temperature, humidity, pressure,
stoichiometry), cell-level EIS, and spatially resolved EIS from a segmented current-density measurement board (CDDM) under the
reference cell. Degradation is concentrated in the cathode catalyst layer and in the cathode outlet region of the active area.

Tests were performed at ZSW (Ulm, Germany) within the RealHyFC project. Methods and results are described in the two Fensterle
papers in [`../references`](../references).

**Merged dataset.** This folder combines the data published with both papers. They describe the same stack: the IJHE paper
(durability) and its companion Fuel Cells paper (beginning-of-life characterization, ref. [27] of the IJHE paper). Files
that appeared in both downloads were checked byte-for-byte and kept once. The companion study adds:
- longer recordings of the BoL nominal polarization curve and the BoL pressure sweep (all rows of the shorter IJHE versions
  are included, with identical values)
- cell-level and spatially resolved EIS at 25 A and 40 A external load (126 A and 252 A nominal, 0.45 and 0.9 A/cm²) at BoL

## Nominal operating conditions (RealHyFC)

| | Inlet temperature / °C | Humidity / % | Outlet pressure / kPaG | Stoichiometry / - | Gas |
|---|---|---|---|---|---|
| Anode | 75 | 50 | 160 | 1.5 | 80 % H2 / 20 % N2 |
| Cathode | 75 | 35 | 150 | 1.7 | Air |
| Coolant | 70 | - | - | - | - |

## Tests

1. [Polarization Curve Nominal](Polarization%20Curve%20Nominal) — 11 polarization curves from BoL to 2300 h, incl. 2 after air/air start
2. [Partial Sensitivity Test BoL](Partial%20Sensitivity%20Test%20BoL) — T, RH, p, stoich sweeps at beginning of life
3. [Partial Sensitivity Test 2300h](Partial%20Sensitivity%20Test%202300h) — the same sweeps after 2300 h
4. [EIS Cell](EIS%20Cell) — cell-level impedance spectra (BoL at 4 loads, 1150 h, 2300 h, and under sweeps)
5. [EIS Spatially Resolved](EIS%20Spatially%20Resolved) — 16 local spectra + reference cell (BoL at 3 loads, 1150 h, 2300 h, and under sweeps)
6. [Load Cycle](Load%20Cycle) — the raw durability-test log, cycles 0–200 (5.4 GB)

File naming of the EIS measurements is explained in [`ReadMe_EIS_measurements.txt`](ReadMe_EIS_measurements.txt) (from the data authors).

## Data format (all CSV files)

- Delimiter is **semicolon** (`;`), decimal point is `.`, file starts with a UTF-8 BOM.
- Files in *Partial Sensitivity Test 2300h* and some polarization curves write hour values with a **thousands comma**
  (`2,281.63`). Strip the comma before converting to a number.
- *Partial Sensitivity Test 2300h* has one extra column (`OpHrs2 [h]`), so always select columns **by header name**, not position.
- Sampling is about one row per second; each load point is held for several minutes. Extract steady-state points
  (e.g. last 60 s of each hold) before fitting.
- The hour labels in file names (300h, 1150h, …) are load-cycle hours; the `OpHrs [h]` column also counts break-in and
  characterization time, so it reads higher.

Main columns (Polarization Curve and Sensitivity files, 904 columns):

| Column | Meaning |
|---|---|
| `DateTime`, `CalHrs [h]`, `OpHrs [h]` | timestamp, calendar hours, operating hours |
| `I.S.set [A]`, `I.S.act [A]`, `j.S.act [A/cm2]` | stack current set/actual, current density |
| `Stoic.S.A`, `Stoic.S.C` (+`.set`) | anode / cathode stoichiometry |
| `U.S.AveCell [V]`, `MAX(U.S.MaxCell)`, `MIN(U.S.MinCell)` | average / max / min cell voltage |
| `U1 [V]` … `U11 [V]` | single-cell voltages. `U6` reads ~0 V throughout (unused channel), so leave it out |
| `F.Si.H.A`, `F.Si.N.A`, `F.Si.Air.C` [nl/min] | anode H2 / N2 and cathode air inlet flow |
| `F.Si.CL [l/min]` | coolant flow |
| `p.Si.*`, `p.So.*` [kPag] | inlet / outlet pressure (`A` anode, `C` cathode, `CL` coolant); `p.S.*.set` setpoints |
| `T.Si.*`, `T.So.*`, `DPT.Si.*` [°C] | inlet / outlet temperature and inlet dew point (sets humidity) |
| `Current_0` … `Current_666` | local currents of the 667 CDDM segments (reference cell) |
| `Temperature_0` … `Temperature_179` | local temperatures from the CDDM |

Abbreviations: `S` stack, `Si` stack inlet, `So` stack outlet, `A` anode, `C` cathode, `CL` coolant loop, `DPT` dew-point temperature.

## Notebooks

Analysis scripts and notebooks for this dataset go in [Notebooks](Notebooks).

## Cite

<pre>
E. Fensterle, J. Hunger, F. Rabsahl, F. Wilhelm, C. K. Christensen, A. Benson, R. Al-Khamissi, M. Tregaro, S. Escribano,
J. Scholta, M. Hölzle, Spatially resolved electrochemical characterization of a PEM fuel cell stack during durability testing,
International Journal of Hydrogen Energy 250 (2026) 156060. doi:10.1016/j.ijhydene.2026.156060
</pre>
<pre>
E. Fensterle, F. Rabsahl, J. Hunger, F. Wilhelm, J. Scholta, M. Hölzle, Spatially resolved current density and electrochemical
impedance spectroscopy in PEM fuel cells, Fuel Cells 26(2) (2026) e70057. doi:10.1002/fuce.70057
</pre>
