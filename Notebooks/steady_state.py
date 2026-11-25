"""Build the steady-state modelling table from polarization_raw + sensitivity_raw.

Each current hold (constant setpoint for a sustained run) is collapsed to one row, averaged
over the last `WINDOW` seconds of the hold. The 667 current segments and 180 temperature
segments are reduced to the locked compact features; raw segment columns are not carried over.

Output table `steady_state`, one row per hold:
  identity : dataset, age_h, air_air_start, parameter, hold_index, n_window
  op. point: j_stack, I_stack, p_cathode, p_anode, stoic_A, stoic_C, dpt_C, dpt_A,
             T_coolant_in, T_cath_in
  response : U_avecell, P_density
  thermal  : T_mean, T_std, T_p95_p05            (180 CDDM temperature sensors)
  current  : heterogeneity_CoV, flow_gradient, peak_pos, lateral_gradient, lateral_CoV
             (667 CDDM current segments, 29x23, boundary ring excluded -> 27x21)

CLI:
    python steady_state.py                 # process everything
    python steady_state.py --limit-groups 2   # quick test: first 2 source groups
"""
import argparse
import os
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")
TABLE = "steady_state"

NROW, NCOL = 29, 23
CUR = [f"Current_{i}" for i in range(NROW * NCOL)]     # 667
TMP = [f"Temperature_{i}" for i in range(180)]

SET_COL = "I.S.set [A]"
TOL      = 0.5      # A: a setpoint change beyond this starts a new segment
MIN_HOLD = 60      # rows: shorter runs are ramps/transients, skipped
WINDOW   = 60      # rows: average the last WINDOW seconds of each hold
MIN_J    = 0.05    # A/cm2: below this the current map is ~zero -> distribution features undefined

OP = {  # output name : source column
    "j_stack": "j.S.act [A/cm2]", "I_stack": "I.S.act [A]", "U_avecell": "U.S.AveCell [V]",
    "p_cathode": "p.So.C [kPag]", "p_anode": "p.So.A [kPag]",
    "stoic_A": "Stoic.S.A [-]", "stoic_C": "Stoic.S.C [-]",
    "dpt_C": "DPT.Si.C [°C]", "dpt_A": "DPT.Si.A [°C]",
    "T_coolant_in": "T.Si.CL [°C]", "T_cath_in": "T.Si.C [°C]",
}


# --- current-map features -----------------------------------------------------
def interior(flat):
    return np.asarray(flat, float).reshape(NROW, NCOL)[1:NROW-1, 1:NCOL-1]

def current_features(seg_mean, j_stack):
    if j_stack < MIN_J:        # near open-circuit: no meaningful current distribution
        return {k: float("nan") for k in
                ("heterogeneity_CoV", "flow_gradient", "peak_pos", "lateral_gradient", "lateral_CoV")}
    g = interior(seg_mean)
    corr = g.sum() / np.asarray(seg_mean, float).reshape(NROW, NCOL).sum()  # paper Eq.1-2
    g = g / corr
    flow, lat = g.mean(axis=1), g.mean(axis=0)
    return {
        "heterogeneity_CoV": g.std() / g.mean(),
        "flow_gradient":     np.polyfit(np.arange(len(flow)), flow / flow.mean(), 1)[0],
        "peak_pos":          int(np.argmax(flow)),
        "lateral_gradient":  np.polyfit(np.arange(len(lat)), lat / lat.mean(), 1)[0],
        "lateral_CoV":       lat.std() / lat.mean(),
    }

def thermal_features(temp_mean):
    v = np.asarray(temp_mean, float)
    v = v[(v > 30) & (v < 100)]        # drop any dead/implausible sensor
    return {
        "T_mean":     v.mean(),
        "T_std":      v.std(),
        "T_p95_p05":  np.percentile(v, 95) - np.percentile(v, 5),
    }


# --- hold detection -----------------------------------------------------------
def find_holds(setpoint):
    s = np.asarray(setpoint, float)
    bounds = [0] + [i for i in range(1, len(s)) if abs(s[i] - s[i-1]) > TOL] + [len(s)]
    holds = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a >= MIN_HOLD:
            holds.append((max(a, b - WINDOW), b))
    return holds


def process_group(df, ident):
    df = df.reset_index(drop=True)
    rows = []
    for k, (lo, hi) in enumerate(find_holds(df[SET_COL].values)):
        w = df.iloc[lo:hi]
        rec = dict(ident); rec["hold_index"] = k; rec["n_window"] = hi - lo
        for out, src in OP.items():
            rec[out] = float(w[src].mean())
        rec["P_density"] = rec["U_avecell"] * rec["j_stack"]
        rec.update(thermal_features(w[TMP].mean().values))
        rec.update(current_features(w[CUR].mean().values, rec["j_stack"]))
        rows.append(rec)
    return rows


def build(db_path=DB_PATH, limit_groups=None):
    base = ["j.S.act [A/cm2]"] + [SET_COL] + list(OP.values())
    cols = sorted(set(base)) + CUR + TMP
    sel = ",".join(f'"{c}"' for c in cols)
    all_rows = []
    with closing(sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)) as con:
        groups = []
        for age, air in con.execute(
            "SELECT DISTINCT age_h, air_air_start FROM polarization_raw ORDER BY age_h, air_air_start"):
            groups.append(("polarization_raw", f'age_h={age} AND air_air_start={air}',
                           {"dataset": "polarization", "age_h": age, "air_air_start": air, "parameter": None}))
        for age, par in con.execute(
            "SELECT DISTINCT age_h, parameter FROM sensitivity_raw ORDER BY age_h, parameter"):
            groups.append(("sensitivity_raw", f"age_h={age} AND parameter='{par}'",
                           {"dataset": "sensitivity", "age_h": age, "air_air_start": None, "parameter": par}))
        if limit_groups:
            groups = groups[:limit_groups]
        for tbl, where, ident in groups:
            df = pd.read_sql_query(f'SELECT {sel} FROM {tbl} WHERE {where}', con)
            r = process_group(df, ident)
            all_rows += r
            print(f"  {ident.get('parameter') or ('airair' if ident['air_air_start'] else 'nominal'):12} "
                  f"age {ident['age_h']:>4}: {len(df):>6} rows -> {len(r):>3} holds")
    out = pd.DataFrame(all_rows)
    with closing(sqlite3.connect(db_path)) as con:
        con.execute(f"DROP TABLE IF EXISTS {TABLE}")
        out.to_sql(TABLE, con, index=False)
    print(f"Done: {len(out)} steady-state rows -> {TABLE} ({db_path})")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit-groups", type=int, default=None)
    args = ap.parse_args()
    build(limit_groups=args.limit_groups)
