"""Export the tables the MATLAB optimization needs as CSV files.

Writes to  <project>/MATLAB/data/:
- holds.csv        every steady-state hold from steady_state (polarization curves at 9 ages, sensitivity sweeps at BoL and 2300 h)
- hfr_anchor.csv   reference values for the ohmic resistance from Tier-1 EIS (10-cell durability stack only)

Run after steady_state.py, eis_tier1.py and eis_sensitivity.py:   python export_matlab.py
"""
import os, sqlite3
from contextlib import closing
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")
OUT_DIR = os.path.join(PROJECT_DIR, "MATLAB")
CELL_AREA_CM2 = 280.0
NOMINAL = {"p_c": 150.0, "T_cin": 70.0, "lambda_c": 1.7}


def psat_kpa(t_c):
    """Saturation pressure of water (Magnus form), kPa."""
    return 0.61094 * np.exp(17.625 * t_c / (t_c + 243.04))


def read(sql, db_path=DB_PATH):
    with closing(sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)) as con:
        return pd.read_sql_query(sql, con)


def export_holds(db_path=DB_PATH):
    s = read("SELECT * FROM steady_state", db_path)
    out = pd.DataFrame({
        "dataset": s.dataset, "age_h": s.age_h, "air_air_start": s.air_air_start.fillna(0).astype(int),
        "parameter": s.parameter.fillna("none"), "j": s.j_stack, "U": s.U_avecell,
        "p_c": s.p_cathode, "p_a": s.p_anode, "lambda_c": s.stoic_C, "lambda_a": s.stoic_A,
        "dpt_c": s.dpt_C, "dpt_a": s.dpt_A, "T_cin": s.T_coolant_in, "T_mean": s.T_mean, "T_std": s.T_std,
    })
    # cathode inlet RH referenced to the coolant inlet temperature (the test held this ratio during the temperature sweep)
    out["RH"] = psat_kpa(out.dpt_c) / psat_kpa(out.T_cin)
    return out


def export_hfr(holds, db_path=DB_PATH):
    # dew point used at each age on the nominal polarization curve
    pol = holds[(holds.dataset == "polarization") & (holds.air_air_start == 0)]
    dpt_by_age = pol.groupby("age_h").dpt_c.median()

    # 1) nominal cell spectra: mean over the 10 cells at each age and current
    t = read("SELECT age_h, dc_current_A, HFR_1k_asr FROM eis_tier1 "
             "WHERE kind='cell' AND measurement='nominal' AND is_multi_cell=0", db_path)
    t["I"] = t.dc_current_A.round(-1)
    a = t.groupby(["age_h", "I"]).agg(HFR=("HFR_1k_asr", "mean"), n=("HFR_1k_asr", "size")).reset_index()
    a = a.assign(source="cell_nominal", j=a.I / CELL_AREA_CM2, **NOMINAL)
    a["RH"] = [psat_kpa(dpt_by_age.reindex([g], method="nearest").iloc[0]) / psat_kpa(70.0) for g in a.age_h]

    # 2) reference-cell sensitivity spectra on the durability stack (2300 h variations)
    e = read("SELECT * FROM eis_sensitivity WHERE scope='cell' AND same_stack=1 AND is_nominal=0", db_path)
    rows = []
    for _, r in e.iterrows():
        c = dict(NOMINAL, RH=psat_kpa(dpt_by_age.loc[2300]) / psat_kpa(70.0))
        key = {"Pressure": "p_c", "Temperature": "T_cin", "Stoichiometry": "lambda_c"}.get(r.parameter)
        if key is None:
            continue
        c[key] = r.level
        rows.append(dict(source="cell_sensitivity", age_h=r.age_h, j=1.8, HFR=r.HFR, n=1, **c))
    b = pd.DataFrame(rows)

    cols = ["source", "age_h", "j", "p_c", "T_cin", "lambda_c", "RH", "HFR", "n"]
    out = pd.concat([a[cols], b[cols]], ignore_index=True)
    out["HFR"] = out.HFR / 1e3            # mOhm cm2 -> Ohm cm2
    return out.dropna(subset=["HFR"])


def build(db_path=DB_PATH, out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    holds = export_holds(db_path)
    hfr = export_hfr(holds, db_path)
    holds.to_csv(os.path.join(out_dir, "holds.csv"), index=False)
    hfr.to_csv(os.path.join(out_dir, "hfr_anchor.csv"), index=False)
    print(f"holds.csv: {len(holds)} rows | hfr_anchor.csv: {len(hfr)} rows -> {out_dir}")
    return holds, hfr


if __name__ == "__main__":
    build()
