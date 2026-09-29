#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

baseline_dir = "/home/pchakrab/input/solar/C12L181/after/new-style"
for state in ["import", "internal", "export"]:
    print("\nSTATE:", state)
    baseline = f"{baseline_dir}/SOLAR_{state}_after_runPhase1.nc"
    current = f"checkpoints/last/SOLAR_{state}.nc"
    print("BASELINE:", baseline)
    print("CURRENT:", current)
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if var in cur.variables:
                if not np.issubdtype(bas[var].dtype, np.number):
                    continue
                if var in bas.dimensions or var in ("time", "lons", "lats", "corner_lons", "corner_lats", "anchor", "contacts", "orientation"):
                    continue
                huge = np.finfo(np.float32).max
                bas_raw = np.ma.masked_invalid(np.ma.asarray(bas[var][:], dtype=np.float64))
                cur_raw = np.ma.masked_invalid(np.ma.asarray(cur[var][0], dtype=np.float64))
                # undef on either side: baseline masked (_FillValue) or |value| == huge
                # (MAPL2 MAPL_UNDEF = +huge, MAPL3 MAPL_UNDEF = -huge)
                bas_undef = np.ma.getmaskarray(bas_raw) | (np.abs(np.ma.filled(bas_raw, 0.0)) == huge)
                cur_undef = np.ma.getmaskarray(cur_raw) | (np.abs(np.ma.filled(cur_raw, 0.0)) == huge)
                bas_var = np.ma.filled(bas_raw, 0.0)
                cur_var = np.ma.filled(cur_raw, 0.0)
                # a baseline _FillValue-masked cell holding the same literal value on
                # the current side (e.g. RI/RL = 1e15) is not a real mismatch
                bas_fill = getattr(bas[var], "_FillValue", None)
                if bas_fill is not None:
                    same_fill = np.ma.getmaskarray(bas_raw) & (cur_var == bas_fill)
                    cur_undef = cur_undef | same_fill
                undef_mismatch = bas_undef != cur_undef
                if undef_mismatch.any():
                    print(f" {var:<20} UNDEF MISMATCH at {undef_mismatch.sum()} cells "
                          f"(bas undef only: {(bas_undef & ~cur_undef).sum()}, cur undef only: {(cur_undef & ~bas_undef).sum()})")
                diff = np.where(bas_undef | cur_undef, 0.0, np.abs(bas_var - cur_var))
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<20} diff (min/max/norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
