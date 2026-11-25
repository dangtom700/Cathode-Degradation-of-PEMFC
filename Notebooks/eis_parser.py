"""Parse the Zahner Zennium EIS `.txt` files into the SQLite database.

The EIS measurements (``EIS Cell/`` and ``EIS Spatially Resolved/``) are raw ASCII
files, not the semicolon CSVs the other datasets use. Each file is a Zahner export:
a ``Key...: value`` header block, then a whitespace-delimited table

    Number  Frequency/Hz  impedance/Ohm  Phase/deg  Significance  Time/s  Voltage/V  Current/A

Impedance is magnitude |Z| and phase phi (deg); this module adds the Cartesian form
Z' = |Z| cos(phi), Z'' = |Z| sin(phi) used for Nyquist plots and circuit fitting.

Two entry points:
    parse_eis_file(path)   -> (meta: dict, df: DataFrame)   # the text parser itself
    load_eis_data()        -> builds table ``eis_raw`` (one row per frequency point)

CLI:
    python eis_parser.py                 # load every EIS file into eis_raw
    python eis_parser.py --limit 20      # quick test: first 20 files only
"""
import argparse
import math
import os
import re
import sqlite3

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
EIS_DIRS = {
    "cell": os.path.join(PROJECT_DIR, "EIS Cell"),
    "spatial": os.path.join(PROJECT_DIR, "EIS Spatially Resolved"),
}
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")
TABLE = "eis_raw"

DATA_COLS = ["number", "frequency_Hz", "Zmag_Ohm", "phase_deg",
             "significance", "time_s", "voltage_V", "current_A"]
MEASUREMENTS = ["Temperature", "Pressure", "Stoichiometry", "Humidity"]

_FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _first_float(text):
    m = re.search(_FLOAT, text or "")
    return float(m.group()) if m else None


# --- the text parser ----------------------------------------------------------

def parse_eis_file(path):
    """Read one Zahner EIS .txt -> (meta dict, DataFrame with Z', Z'')."""
    with open(path, encoding="latin-1") as f:
        lines = f.read().splitlines()

    # split header block from the data table (data header starts with "Number")
    data_start = next(
        i for i, ln in enumerate(lines)
        if ln.strip().startswith("Number") and "Frequency" in ln
    )
    header = {}
    for ln in lines[:data_start]:
        if ":" in ln:
            key, _, val = ln.partition(":")
            header[key.strip(" .")] = val.strip()

    current_line = header.get("Current", "")
    ampl = re.search(r"Ampl:\s*(" + _FLOAT + ")", current_line)
    meta = {
        "potential_V": _first_float(header.get("Potential")),
        "dc_current_A": _first_float(current_line),
        "amplitude_A": float(ampl.group(1)) if ampl else None,
        "area_cm2": _first_float(header.get("(area/cm^2)")),
        "temperature_raw": header.get("Temperature"),
        "meas_date": (header.get("File", "").split("-", 1)[-1].strip() or None),
        "meas_time": header.get("Time"),
    }

    rows = []
    for ln in lines[data_start + 1:]:
        parts = ln.split()
        if len(parts) >= len(DATA_COLS) and re.match(r"^\d+$", parts[0]):
            rows.append([float(x) for x in parts[:len(DATA_COLS)]])
    df = pd.DataFrame(rows, columns=DATA_COLS)
    phi = df["phase_deg"] * math.pi / 180.0
    df["Zreal_Ohm"] = df["Zmag_Ohm"] * phi.apply(math.cos)
    df["Zimag_Ohm"] = df["Zmag_Ohm"] * phi.apply(math.sin)
    return meta, df


# --- metadata from the folder / file name -------------------------------------

def eis_metadata_from_path(kind, root, path):
    rel = os.path.relpath(path, root).replace("\\", "/")
    fname = os.path.basename(path)

    age = re.search(r"(?:^|[_/])(BoL|\d+h)(?:[_/]|$)", rel)
    age_h = 0 if (age and age.group(1) == "BoL") else (
        int(age.group(1).rstrip("h")) if age else None)

    measurement = next((m for m in MEASUREMENTS if f"{m}Variation" in rel), "nominal")

    # condition label: nested dir (spatial) or file-name prefix (cell); None if nominal
    condition = None
    if measurement != "nominal":
        parent = os.path.basename(os.path.dirname(path))
        if parent.endswith(("_BoL",)) or re.search(r"_\d+h$", parent) or "Variation" in parent:
            condition = fname.split("_")[0]        # cell: 't55', 'p210', 'stoic14'
        else:
            condition = parent                     # spatial: 'T1_55', 'RH_30_30', 'p2_210_200'

    tst = re.search(r"tst(\d+)", fname)
    tst_index = int(tst.group(1)) if tst else None

    return {
        "kind": kind,
        "age_h": age_h,
        "measurement": measurement,
        "condition": condition,
        "tst_index": tst_index,
        # spatial tst00 is the reference cell; tst01-16 are the segments
        "is_reference": int(kind == "spatial" and tst_index == 0),
        "source_path": rel,
    }


# --- loader -------------------------------------------------------------------

META_COLS = [
    ("kind", "TEXT"), ("age_h", "INTEGER"), ("measurement", "TEXT"),
    ("condition", "TEXT"), ("tst_index", "INTEGER"), ("is_reference", "INTEGER"),
    ("source_path", "TEXT"), ("potential_V", "REAL"), ("dc_current_A", "REAL"),
    ("amplitude_A", "REAL"), ("area_cm2", "REAL"), ("temperature_raw", "TEXT"),
    ("meas_date", "TEXT"), ("meas_time", "TEXT"),
]
DATA_SQL = [("number", "INTEGER"), ("frequency_Hz", "REAL"), ("Zmag_Ohm", "REAL"),
            ("phase_deg", "REAL"), ("significance", "REAL"), ("time_s", "REAL"),
            ("voltage_V", "REAL"), ("current_A", "REAL"),
            ("Zreal_Ohm", "REAL"), ("Zimag_Ohm", "REAL")]


def _iter_files():
    for kind, root in EIS_DIRS.items():
        for dirpath, _, filenames in os.walk(root):
            for fn in sorted(filenames):
                if fn.lower().endswith(".txt"):
                    yield kind, root, os.path.join(dirpath, fn)


def load_eis_data(limit=None):
    columns = [c for c, _ in META_COLS] + [c for c, _ in DATA_SQL]
    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute(f"DROP TABLE IF EXISTS {TABLE}")
    con.execute(f"CREATE TABLE {TABLE} ("
                + ", ".join(f'"{c}" {t}' for c, t in META_COLS + DATA_SQL) + ")")
    insert = f"INSERT INTO {TABLE} VALUES ({', '.join('?' * len(columns))})"

    n_files = n_rows = 0
    for kind, root, path in _iter_files():
        if limit is not None and n_files >= limit:
            break
        meta = eis_metadata_from_path(kind, root, path)
        hdr, df = parse_eis_file(path)
        meta.update(hdr)
        base = [meta[c] for c, _ in META_COLS]
        con.executemany(insert, [base + list(r) for r in df.itertuples(index=False)])
        n_files += 1
        n_rows += len(df)
    con.commit()
    con.execute(f"CREATE INDEX IF NOT EXISTS ix_{TABLE} ON {TABLE}(kind, age_h, measurement)")
    con.commit()
    con.close()
    print(f"Done: {n_rows:,} points from {n_files} files in {TABLE} ({DB_PATH})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None,
                    help="load only the first N files (quick test)")
    args = ap.parse_args()
    load_eis_data(limit=args.limit)
