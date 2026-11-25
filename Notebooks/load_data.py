"""Load the Cathode Degradation CSVs into one SQLite database, ``cathode_degradation.db``.

Three tables, one per dataset group (raw rows, column names kept exactly as in the
CSV headers with the BOM stripped, plus metadata columns parsed from the file names):

    load_cycle_raw     raw durability-test log, cycles 0-200 (~9.1 M rows, a few GB)
    polarization_raw   11 nominal polarization curves, BoL -> 2300 h
                       + age_h, air_air_start, source_file
    sensitivity_raw    partial sensitivity sweeps (T/RH/p/stoich) at BoL and 2300 h
                       + age_h, parameter, source_file

Usage (from any folder):
    python load_data.py                      # polarization + sensitivity (fast)
    python load_data.py --datasets load_cycle polarization sensitivity   # everything
    python load_data.py --main-only          # skip the CDDM Current_*/Temperature_* segment
                                             # columns (667+180 of them) -> much smaller tables
    python load_data.py --nrows 500          # quick test: first 500 rows per file

Each table is dropped and rebuilt on every run. Files with differing headers
(e.g. the 2300 h sensitivity files carry an extra ``OpHrs2 [h]`` column) are merged
by column *name*, so missing columns are stored as NULL rather than mis-aligned.
"""
import argparse
import os
import re
import sqlite3
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
LOAD_CYCLE_DIR = os.path.join(PROJECT_DIR, "Load Cycle")
POLARIZATION_DIR = os.path.join(PROJECT_DIR, "Polarization Curve Nominal")
SENSITIVITY_DIRS = {
    0: os.path.join(PROJECT_DIR, "Partial Sensitivity Test BoL"),
    2300: os.path.join(PROJECT_DIR, "Partial Sensitivity Test 2300h"),
}
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")
CHUNK_ROWS = 100_000

TEXT_COLS = {"DateTime", "source_file", "parameter"}
INT_COLS = {"Number of cells", "operation_point", "age_h", "air_air_start"}

# CDDM per-segment channels: Current_0..666 (667) and Temperature_0..179 (180).
# Only needed for spatially resolved analysis; --main-only drops them.
SEGMENT_RE = re.compile(r"^(Current|Temperature)_\d+$")


def first_cycle(fname):
    """'Loadcycle_51-70.csv' -> 51, so files load in chronological order."""
    return int(fname.split("_")[1].split("-")[0])


def read_header(path):
    with open(path, encoding="utf-8-sig") as f:
        return f.readline().rstrip("\r\n").split(";")


def union_columns(headers):
    """Merge headers in order, placing new columns after their predecessor."""
    cols = []
    for header in headers:
        for i, col in enumerate(header):
            if col in cols:
                continue
            if i == 0:
                cols.insert(0, col)
            else:
                cols.insert(cols.index(header[i - 1]) + 1, col)
    return cols


def quote(name):
    return '"' + name.replace('"', '""') + '"'


def sql_type(col):
    if col in TEXT_COLS:
        return "TEXT"
    if col in INT_COLS:
        return "INTEGER"
    return "REAL"


def load_cycle_data(nrows=None):
    files = sorted(
        (f for f in os.listdir(LOAD_CYCLE_DIR) if f.lower().endswith(".csv")),
        key=first_cycle,
    )
    paths = [os.path.join(LOAD_CYCLE_DIR, f) for f in files]
    columns = union_columns(read_header(p) for p in paths) + ["source_file"]

    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute("DROP TABLE IF EXISTS load_cycle_raw")
    con.execute(
        "CREATE TABLE load_cycle_raw ("
        + ", ".join(f"{quote(c)} {sql_type(c)}" for c in columns)
        + ")"
    )
    insert = f"INSERT INTO load_cycle_raw VALUES ({', '.join('?' * len(columns))})"

    total = 0
    t0 = time.time()
    for fname, path in zip(files, paths):
        n_file = 0
        reader = pd.read_csv(
            path,
            sep=";",
            encoding="utf-8-sig",
            thousands=",",          # '2,678.46' -> 2678.46 (safe: delimiter is ';')
            dtype={"DateTime": str},
            chunksize=CHUNK_ROWS,
            nrows=nrows,
            low_memory=False,
        )
        for chunk in reader:
            chunk["source_file"] = fname
            chunk = chunk.reindex(columns=columns)
            # astype(object) turns NaN into None-compatible values; sqlite3
            # stores float NaN as NULL anyway.
            con.executemany(insert, chunk.astype(object).itertuples(index=False))
            n_file += len(chunk)
        con.commit()
        total += n_file
        print(f"{fname:<24} {n_file:>9,} rows   total {total:>10,}   {time.time() - t0:6.0f} s")

    con.close()
    print(f"Done: {total:,} rows in load_cycle_raw ({DB_PATH})")


# --- characterization datasets (polarization + sensitivity) --------------------

def polarization_meta(fname):
    """'Polcurve_Nominal_RHFC_1500h_afterAirAirStart.csv' -> age_h=1500, air_air_start=1."""
    core = fname[:-4].replace("Polcurve_Nominal_RHFC_", "")
    air = int(core.endswith("_afterAirAirStart"))
    core = core.replace("_afterAirAirStart", "")
    age = 0 if core == "BoL" else int(core.rstrip("h"))
    return {"age_h": age, "air_air_start": air}


def sensitivity_meta(fname):
    """'PartialSensitivityTest_Humidity_2300h.csv' -> parameter='Humidity', age_h=2300."""
    _, parameter, age = fname[:-4].split("_")
    return {"parameter": parameter, "age_h": 0 if age == "BoL" else int(age.rstrip("h"))}


def _build_table(table, items, nrows, main_only, index_cols=None):
    """items: list of (fname, path, meta_dict). Metadata keys must match across items."""
    paths = [p for _, p, _ in items]
    data_cols = union_columns(read_header(p) for p in paths)
    if main_only:
        data_cols = [c for c in data_cols if not SEGMENT_RE.match(c)]
    meta_keys = list(items[0][2].keys())
    columns = data_cols + meta_keys + ["source_file"]

    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute(f"DROP TABLE IF EXISTS {table}")
    con.execute(
        f"CREATE TABLE {table} ("
        + ", ".join(f"{quote(c)} {sql_type(c)}" for c in columns)
        + ")"
    )
    insert = f"INSERT INTO {table} VALUES ({', '.join('?' * len(columns))})"

    total = 0
    t0 = time.time()
    for fname, path, meta in items:
        n_file = 0
        reader = pd.read_csv(
            path,
            sep=";",
            encoding="utf-8-sig",
            thousands=",",
            dtype={"DateTime": str},
            chunksize=CHUNK_ROWS,
            nrows=nrows,
            low_memory=False,
        )
        for chunk in reader:
            # Attach metadata in one concat: inserting columns one by one into a
            # ~900-column frame fragments it (pandas PerformanceWarning).
            extra = pd.DataFrame({**meta, "source_file": fname}, index=chunk.index)
            chunk = pd.concat(
                [chunk.reindex(columns=data_cols), extra], axis=1
            )   # reindex drops segment cols when main_only
            con.executemany(insert, chunk.astype(object).itertuples(index=False))
            n_file += len(chunk)
        con.commit()
        total += n_file
        print(f"{fname:<48} {n_file:>7,} rows   total {total:>9,}   {time.time() - t0:6.0f} s")

    if index_cols:
        idx = ", ".join(quote(c) for c in index_cols)
        con.execute(f"CREATE INDEX IF NOT EXISTS ix_{table} ON {table}({idx})")
        con.commit()
    con.close()
    print(f"Done: {total:,} rows in {table} ({DB_PATH})")


def load_polarization_data(nrows=None, main_only=False):
    files = [f for f in os.listdir(POLARIZATION_DIR) if f.lower().endswith(".csv")]
    items = [(f, os.path.join(POLARIZATION_DIR, f), polarization_meta(f)) for f in files]
    items.sort(key=lambda it: (it[2]["age_h"], it[2]["air_air_start"]))
    _build_table("polarization_raw", items, nrows, main_only,
                 index_cols=["age_h", "air_air_start"])


def load_sensitivity_data(nrows=None, main_only=False):
    items = []
    for age_h, directory in SENSITIVITY_DIRS.items():
        for f in os.listdir(directory):
            if f.lower().endswith(".csv"):
                items.append((f, os.path.join(directory, f), sensitivity_meta(f)))
    items.sort(key=lambda it: (it[2]["age_h"], it[2]["parameter"]))
    _build_table("sensitivity_raw", items, nrows, main_only,
                 index_cols=["age_h", "parameter"])


LOADERS = {
    "load_cycle": load_cycle_data,
    "polarization": load_polarization_data,
    "sensitivity": load_sensitivity_data,
}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--nrows", type=int, default=None,
                    help="rows per file (quick test; default: all)")
    ap.add_argument("--main-only", action="store_true",
                    help="skip CDDM Current_*/Temperature_* segment columns "
                         "(polarization/sensitivity only)")
    ap.add_argument("--datasets", nargs="+", choices=list(LOADERS),
                    default=["polarization", "sensitivity"],
                    help="which datasets to load (default: polarization sensitivity; "
                         "load_cycle is ~9 M rows)")
    args = ap.parse_args()

    for name in args.datasets:
        if name == "load_cycle":
            LOADERS[name](nrows=args.nrows)
        else:
            LOADERS[name](nrows=args.nrows, main_only=args.main_only)
