"""Tier-1 (fit-free) EIS resistances from the Nyquist intercepts, after Fensterle et al. (Fig. 9a).

For every spectrum in ``eis_raw`` this reads three real-axis values directly off the Nyquist plot
(no equivalent-circuit fit):

    HFR  high-frequency resistance  -> ohmic / membrane
         HFR_1k  : Z' interpolated at 1 kHz (fixed-frequency proxy; industry convention)
         HFR_420 : Z' at 420 Hz (common proxy for all zones, whose window stops at 421 Hz)
         HFR_ext : the ORR arc's high-frequency side extrapolated to Z''=0 with a circle fit
                   (the measured window stops before the arc reaches the real axis)
    IFR  intermediate-frequency resistance: Z' at the valley of -Z'' between the ORR and
         mass-transport arcs (NaN when the arcs are not separated, e.g. at low current)
    LFR  low-frequency resistance: Z' where -Z'' crosses zero on the low-frequency side

    R_ORR_t1 = IFR - HFR_ext      (ORR charge transfer)
    R_MT_t1  = LFR - IFR          (mass transport)

Corrections applied (see eis_analysis.ipynb):
  * spatial zone spectra are recorded ~180 deg rotated -> negated before reading intercepts
  * frequency windows from the paper: cell 0.1-1180 Hz; zones 1, 3, 10 (cathode inlet)
    0.1-1180 Hz; all other zones 0.1-421 Hz
  * cell-folder files spanning several cells (is_multi_cell=1, the 1150 h tst00) are skipped

Output table ``eis_tier1`` (one row per spectrum, resistances in Ohm; ``*_asr`` = x area, mOhm cm2).

CLI:
    python eis_tier1.py
"""
import os
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(HERE)
DB_PATH = os.path.join(PROJECT_DIR, "cathode_degradation.db")
TABLE = "eis_tier1"

CELL_AREA_CM2 = 280.0
FMIN = 0.1
FMAX_CELL = 1185.0          # include the 1180.4 Hz grid point (paper: "1180 Hz")
INLET_ZONES = {1, 3, 10}
FMAX_ZONE_INLET, FMAX_ZONE_OTHER = 1185.0, 425.0   # paper: 1180 Hz / 421 Hz

ZONES = {1: (24, 75), 2: (116, 213), 3: (31, 80), 4: (123, 218), 5: (305, 355), 6: (438, 534),
         7: (576, 626), 8: (444, 539), 9: (127, 222), 10: (40, 90), 11: (132, 228), 12: (316, 366),
         13: (448, 543), 14: (586, 635), 15: (453, 550), 16: (591, 642)}
SEG_AREA_CM2 = CELL_AREA_CM2 / (27 * 21)      # interior segment area (567 interior segments = active area)


def zone_area_cm2(z):
    s, e = ZONES[z]
    (r0, c0), (r1, c1) = divmod(s, 23), divmod(e, 23)
    return (r1 - r0 + 1) * (c1 - c0 + 1) * SEG_AREA_CM2


# --- intercept readers (inputs sorted high -> low frequency) ----------------------
def _prep(f, zr, zi, fmax, negate=False):
    f, zr, zi = map(np.asarray, (f, zr, zi))
    m = (f >= FMIN) & (f <= fmax)
    f, zr, zi = f[m], zr[m], zi[m]
    if negate:
        zr, zi = -zr, -zi
    o = np.argsort(f)[::-1]
    return f[o], zr[o], -zi[o]                 # return -Z'' (positive on capacitive arcs)


def hfr_fixed(f, zr, f0=1000.0):
    lf = np.log10(f); o = np.argsort(lf)
    if not (lf.min() <= np.log10(f0) <= lf.max()):
        return np.nan
    return float(np.interp(np.log10(f0), lf[o], zr[o]))


def _first_peak(nz):
    """index of the first local maximum of -Z'' coming from the high-frequency end."""
    for i in range(1, len(nz) - 1):
        if nz[i] >= nz[i - 1] and nz[i] > nz[i + 1]:
            return i
    return int(np.argmax(nz))


def hfr_extrap(f, zr, nz):
    """circle (Kasa) fit through the high-f side of the ORR arc, left real-axis root."""
    k = _first_peak(nz)
    x, y = zr[:k + 1], nz[:k + 1]
    if len(x) < 4:
        return np.nan
    A = np.column_stack([x, y, np.ones_like(x)])
    D, E, F = np.linalg.lstsq(A, -(x**2 + y**2), rcond=None)[0]
    disc = D**2 - 4 * F
    if disc < 0:
        return np.nan
    root = (-D - np.sqrt(disc)) / 2
    # sanity: left of the measured arc, and the extrapolated gap no longer than the measured
    # high-frequency half of the arc (rejects degenerate circle fits)
    if not (0 < root <= x.min() and x.min() - root <= (x.max() - x.min())):
        return np.nan
    return float(root)


def ifr_valley(f, zr, nz, band=(3.0, 200.0)):
    m = np.where((f >= band[0]) & (f <= band[1]))[0]
    if len(m) < 3:
        return np.nan, np.nan
    i = m[np.argmin(nz[m])]
    # require a genuine valley: higher -Z'' on both sides somewhere in the spectrum
    if i == 0 or i == len(nz) - 1 or nz[:i].max() <= nz[i] or nz[i + 1:].max() <= nz[i]:
        return np.nan, np.nan
    if 0 < i < len(f) - 1:                                  # parabolic refinement in log f
        x = np.log10(f[i - 1:i + 2]); y = nz[i - 1:i + 2]
        a, b, _ = np.polyfit(x, y, 2)
        if a > 0:
            xv = np.clip(-b / (2 * a), x.min(), x.max())
            o = np.argsort(x)
            return float(np.interp(xv, x[o], zr[i - 1:i + 2][o])), float(10**xv)
    return float(zr[i]), float(f[i])


def lfr_crossing(f, zr, nz):
    start = len(nz) // 2
    k0 = start + int(np.argmax(nz[start:]))                 # mass-transport arc apex
    for k in range(k0, len(nz) - 1):
        if nz[k] > 0 >= nz[k + 1]:
            t = nz[k] / (nz[k] - nz[k + 1])
            return float(zr[k] + t * (zr[k + 1] - zr[k])), True
    return float(zr[-1]), False                             # no crossing: lowest-frequency value


def tier1(f, zr, zi, fmax, negate=False):
    f, zr, nz = _prep(f, zr, zi, fmax, negate)
    if len(f) < 8:
        return None
    h1k = hfr_fixed(f, zr)
    h420 = hfr_fixed(f, zr, 420.0)       # common-frequency proxy usable for every zone (zone window: 421 Hz)
    hx = hfr_extrap(f, zr, nz)
    ifr, f_ifr = ifr_valley(f, zr, nz)
    lfr, crossed = lfr_crossing(f, zr, nz)
    return {"HFR_1k": h1k, "HFR_420": h420, "HFR_ext": hx, "IFR": ifr, "f_valley_Hz": f_ifr, "LFR": lfr,
            "lfr_crossed": int(crossed),
            "single_arc": int(not np.isfinite(ifr)),   # no valley: ORR and MT arcs not separated
            "R_ORR_t1": ifr - hx if np.isfinite(ifr) and np.isfinite(hx) else np.nan,
            "R_MT_t1": lfr - ifr if np.isfinite(ifr) else np.nan,
            "f_max_used": float(f.max())}


# --- table builder ----------------------------------------------------------------
def build(db_path=DB_PATH, write=True):
    with closing(sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)) as con:
        d = pd.read_sql_query(
            "SELECT kind, age_h, measurement, condition, tst_index, is_reference, is_multi_cell, channel, "
            "source_path, dc_current_A, potential_V, data_format, frequency_Hz, Zreal_Ohm, Zimag_Ohm "
            "FROM eis_raw", con)
    rows = []
    keys = ["kind", "age_h", "measurement", "condition", "tst_index", "is_reference",
            "is_multi_cell", "channel", "source_path", "dc_current_A", "potential_V", "data_format"]
    for k, g in d.groupby("source_path"):
        meta = g.iloc[0][keys].to_dict()
        if meta["is_multi_cell"]:
            continue
        spatial_zone = meta["kind"] == "spatial" and not meta["is_reference"]
        z = int(meta["tst_index"]) if spatial_zone else None
        fmax = (FMAX_ZONE_INLET if z in INLET_ZONES else FMAX_ZONE_OTHER) if spatial_zone else FMAX_CELL
        r = tier1(g.frequency_Hz.values, g.Zreal_Ohm.values, g.Zimag_Ohm.values, fmax, negate=spatial_zone)
        if r is None:
            continue
        area = zone_area_cm2(z) if spatial_zone else CELL_AREA_CM2
        rec = {**meta, "zone": z, "area_cm2": area, **r}
        for c in ("HFR_1k", "HFR_420", "HFR_ext", "IFR", "LFR", "R_ORR_t1", "R_MT_t1"):
            rec[c + "_asr"] = rec[c] * area * 1e3 if np.isfinite(rec[c]) else np.nan   # mOhm cm2
        rows.append(rec)
    out = pd.DataFrame(rows)
    if write:
        with closing(sqlite3.connect(db_path)) as con:
            con.execute(f"DROP TABLE IF EXISTS {TABLE}")
            out.to_sql(TABLE, con, index=False)
        print(f"Done: {len(out)} spectra -> {TABLE} ({db_path})")
    return out


if __name__ == "__main__":
    build()
