#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

baseline_dir = "/home/pchakrab/input/radiation/C12L181/after/new-style"
comp = "IRRAD"
states = ["import", "internal", "export"]

HUGE = np.finfo(np.float32).max
SKIP = {"time", "lons", "lats", "corner_lons", "corner_lats", "anchor", "contacts", "orientation"}


def undef(x, fills):
    # NaN, MAPL_UNDEF (+huge in MAPL2, -huge in MAPL3), or either file's _FillValue
    u = ~np.isfinite(x) | (np.abs(x) == HUGE)
    for f in fills:
        u |= x == f
    return u


for state in states:
    baseline = f"{baseline_dir}/{comp}_{state}_after_runPhase1.nc"
    current = f"checkpoints/last/{comp}_{state}.nc"
    print(f"\nSTATE: {state}\nBASELINE: {baseline}\nCURRENT: {current}")
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        bas.set_auto_mask(False)
        cur.set_auto_mask(False)
        for var in bas.variables:
            if var not in cur.variables or var in bas.dimensions or var in SKIP:
                continue
            if not np.issubdtype(bas[var].dtype, np.number):
                continue
            b, c = bas[var][:], cur[var][0]
            fills = [getattr(v, "_FillValue") for v in (bas[var], cur[var]) if hasattr(v, "_FillValue")]
            bu, cu = undef(b, fills), undef(c, fills)
            if (bu != cu).any():
                print(f" {var:<20} UNDEF MISMATCH: {(bu & ~cu).sum()} baseline-only, {(cu & ~bu).sum()} current-only")
            diff = np.where(bu | cu, 0.0, np.abs(b.astype(np.float64) - c))
            print(f" {var:<20} diff (min/max/norm): {diff.min()}, {diff.max()}, {np.linalg.norm(diff)}")
