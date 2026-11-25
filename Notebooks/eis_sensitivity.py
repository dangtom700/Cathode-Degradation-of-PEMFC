"""Tier-1 EIS features for the sensitivity tests, and their join to steady_state.

Builds two tables in cathode_degradation.db (needs eis_tier1 and steady_state):

  eis_sensitivity   one row per (age_h, parameter, level, scope); scope = 'cell' (reference-cell
                    spectrum) or a plate region ('1 inlet' ... '5 outlet'). Tier-1 resistances in
                    mOhm cm2 (cell: x 280 cm2; regions: zone impedance x zone area, as in the paper)
  steady_state_eis  the steady_state sensitivity holds at ~1.8 A/cm2, each with the reference-cell
                    and region features of the matching (age_h, parameter, level)

Conventions (details in eis_sensitivity.ipynb):
  * Pressure labels in the EIS file names are absolute kPa. Table 2 of the papers only spans
    50-200 kPaG, so p160 / p1_160_150 -> 50 kPaG (cathode outlet), p210 -> 100, p310 -> 200.
  * Humidity level = cathode inlet RH in % referenced to 70 C (nominal 35): RH_30_30 -> 30,
    RH_70_60 -> 60, RH_90_90 -> 90. In steady_state the level is recovered from the dew point.
  * No variation set contains its nominal point. The nominal spatial EIS of the same age
    (504 A) is added as the nominal anchor (is_nominal = 1).
  * Stacks: the BoL *variation* EIS was measured on a fresh 8-cell short stack (Oct-Nov 2025).
    Everything else - BoL nominal EIS (Jan 2025), all aged EIS, and all polarization /
    sensitivity data behind steady_state - comes from the 10-cell durability stack.
    same_stack = 1 when an EIS row comes from the 10-cell durability stack. Absolute values
    differ between stacks (~13 % in HFR), so compare BoL 8-cell values only within that stack.
  * Region ORR/MT use the single-arc rule from eis_analysis.ipynb, applied within each
    (age, parameter) set: zones whose resolved mass transport stays < 5 mOhm cm2 are MT-free
    (MT = 0, ORR = LFR - HFR_420); other zones without a valley are left unresolved (NaN).

CLI:
    python eis_sensitivity.py
"""
import os
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")

REGION = {1: "1 inlet", 3: "1 inlet", 10: "1 inlet", 2: "2 upper", 4: "2 upper", 9: "2 upper",
          11: "2 upper", 5: "3 middle", 12: "3 middle", 6: "4 lower", 8: "4 lower", 13: "4 lower",
          15: "4 lower", 7: "5 outlet", 14: "5 outlet", 16: "5 outlet"}
NOMINAL = {"Pressure": 150.0, "Temperature": 70.0, "Stoichiometry": 1.7, "Humidity": 35.0}
LEVELS = {"Pressure": [50.0, 100.0, 150.0, 200.0], "Temperature": [55.0, 60.0, 65.0, 70.0, 75.0],
          "Stoichiometry": [1.6, 1.7, 2.1], "Humidity": [30.0, 35.0, 60.0, 90.0]}
MT_FREE_MAX = 5.0
LOAD_A = 505          # 1.8 A/cm2 (50 A external + 454 A test bench)


def eis_level(parameter, condition):
    """EIS condition label -> level in steady_state units."""
    c = condition.lower()
    if parameter == "Pressure":                       # p160 | p1_160_150 (absolute kPa, anode first)
        a = c.split("_")[1] if "_" in c else c[1:]
        return {"160": 50.0, "210": 100.0, "310": 200.0}[a]
    if parameter == "Temperature":                    # t55 | T1_55
        return float(c.split("_")[1] if "_" in c else c[1:])
    if parameter == "Stoichiometry":                  # stoic14 | Stoic1_14_16 (anode 1.4 -> cathode 1.6)
        a = c.split("_")[1] if "_" in c else c.replace("stoic", "")
        return {"14": 1.6, "19": 2.1}[a]
    if parameter == "Humidity":                       # RH_30_30 | RH1_30_30 -> cathode RH
        return float(c.split("_")[2])
    raise ValueError(parameter)


def _psat_kpa(t_c):                                   # Magnus formula, water over liquid
    return 0.61094 * np.exp(17.625 * t_c / (t_c + 243.04))


def steady_level(row):
    """steady_state hold -> (parameter, level) on the same grid as eis_level."""
    par = {"Stoic": "Stoichiometry"}.get(row["parameter"], row["parameter"])
    if par == "Pressure":
        x = row["p_cathode"]
    elif par == "Temperature":
        x = row["T_coolant_in"]
    elif par == "Stoichiometry":
        x = row["stoic_C"]
    else:                                             # humidity: RH of the cathode gas at 70 C
        x = 100 * _psat_kpa(row["dpt_C"]) / _psat_kpa(70.0)
    grid = np.array(LEVELS[par])
    return par, float(grid[np.argmin(np.abs(grid - x))])


def _regions(S):
    """single-arc rule within one (age, parameter) set, then aggregate zones to regions."""
    S = S.copy()
    res = S[S.single_arc == 0].groupby("zone").R_MT_t1_asr.max()
    free = {z for z in S.zone.unique() if z not in res.index or res[z] < MT_FREE_MAX}
    two, fr = S.single_arc == 0, S.zone.isin(free)
    S["ORR"] = np.where(two, S.IFR_asr - S.HFR_420_asr, np.where(fr, S.LFR_asr - S.HFR_420_asr, np.nan))
    S["MT"] = np.where(two, S.R_MT_t1_asr.clip(lower=0), np.where(fr, 0.0, np.nan))
    S["region"] = S.zone.map(REGION)
    g = S.groupby(["level", "is_nominal", "same_stack", "region"])
    out = g.agg(HFR=("HFR_420_asr", "mean"), ORR=("ORR", "mean"), MT=("MT", "mean"),
                LFR=("LFR_asr", "mean"), n_zones=("zone", "size"),
                n_unresolved=("ORR", lambda s: int(s.isna().sum()))).reset_index()
    cov = S.groupby(["level", "is_nominal"]).LFR_asr.agg(lambda s: s.std() / s.mean()).rename("zone_LFR_CoV")
    return out.rename(columns={"region": "scope"}), cov.reset_index()


def build_features(T1):
    var = T1[T1.measurement != "nominal"].copy()
    var["parameter"] = var.measurement
    var["level"] = [eis_level(p, c) for p, c in zip(var.parameter, var.condition)]
    var["is_nominal"] = 0
    var["same_stack"] = (var.age_h != 0).astype(int)          # BoL variations: 8-cell short stack
    nom = T1[(T1.measurement == "nominal") & (T1.kind == "spatial") & (T1.dc_current_A.round() == LOAD_A)].copy()
    nom["is_nominal"], nom["same_stack"] = 1, 1

    rows, covs = [], []
    for (age, par), g in var.groupby(["age_h", "parameter"]):
        n = nom[nom.age_h == age].copy()
        n["parameter"], n["level"] = par, NOMINAL[par]
        S = pd.concat([g, n], ignore_index=True)
        # reference cell: variation spectra live in 'EIS Cell', the anchor is spatial tst00
        cell = S[(S.kind == "cell") | (S.is_reference == 1)].copy()
        cell = cell.assign(scope="cell", HFR=cell.HFR_1k_asr, ORR=cell.IFR_asr - cell.HFR_1k_asr,
                           MT=cell.R_MT_t1_asr, LFR=cell.LFR_asr, n_zones=np.nan,
                           n_unresolved=cell.single_arc)
        cols = ["level", "is_nominal", "same_stack", "scope", "HFR", "ORR", "MT", "LFR", "n_zones", "n_unresolved"]
        reg, cov = _regions(S[(S.kind == "spatial") & (S.is_reference == 0)])
        block = pd.concat([cell[cols], reg[cols]], ignore_index=True)
        block.insert(0, "parameter", par); block.insert(0, "age_h", age)
        cov.insert(0, "parameter", par); cov.insert(0, "age_h", age)
        rows.append(block); covs.append(cov)
    return pd.concat(rows, ignore_index=True), pd.concat(covs, ignore_index=True)


def join_steady(ss, feats, covs, j_target=1.8, j_tol=0.05):
    holds = ss[(ss.dataset == "sensitivity") & ((ss.j_stack - j_target).abs() < j_tol)].copy()
    pl = holds.apply(steady_level, axis=1, result_type="expand")
    holds["eis_parameter"], holds["level"] = pl[0], pl[1]
    holds["is_nominal"] = [int(abs(l - NOMINAL[p]) < 1e-9) for p, l in zip(holds.eis_parameter, holds.level)]
    keys = ["age_h", "parameter", "level"]
    # union of keys: humidity has zone spectra only (no reference-cell file), so it gets region features
    base = feats.groupby(keys).same_stack.max().rename("eis_same_stack").reset_index()
    cell = feats[feats.scope == "cell"].rename(columns={"HFR": "cell_HFR", "ORR": "cell_ORR", "MT": "cell_MT",
                                                        "LFR": "cell_LFR"})
    out_ = feats[feats.scope == "5 outlet"][keys + ["LFR", "MT", "ORR"]].rename(
        columns={"LFR": "outlet_LFR", "MT": "outlet_MT", "ORR": "outlet_ORR"})
    in_ = feats[feats.scope == "1 inlet"][keys + ["LFR"]].rename(columns={"LFR": "inlet_LFR"})
    f = (base.merge(cell[keys + ["cell_HFR", "cell_ORR", "cell_MT", "cell_LFR"]], on=keys, how="left")
         .merge(out_, on=keys, how="left").merge(in_, on=keys, how="left")
         .merge(covs[keys + ["zone_LFR_CoV"]], on=keys, how="left"))
    f["outlet_inlet_LFR"] = f.outlet_LFR / f.inlet_LFR
    f = f.rename(columns={"parameter": "eis_parameter"})
    return holds.merge(f, on=["age_h", "eis_parameter", "level"], how="left")


def build(db_path=DB_PATH, write=True):
    with closing(sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)) as con:
        T1 = pd.read_sql_query("SELECT * FROM eis_tier1", con)
        ss = pd.read_sql_query("SELECT * FROM steady_state", con)
    feats, covs = build_features(T1)
    joined = join_steady(ss, feats, covs)
    if write:
        with closing(sqlite3.connect(db_path)) as con:
            for name, df in (("eis_sensitivity", feats), ("steady_state_eis", joined)):
                con.execute(f"DROP TABLE IF EXISTS {name}")
                df.to_sql(name, con, index=False)
        print(f"Done: eis_sensitivity {len(feats)} rows, steady_state_eis {len(joined)} rows "
              f"({joined.outlet_LFR.notna().sum()} with EIS features) -> {db_path}")
    return feats, covs, joined


if __name__ == "__main__":
    build()
