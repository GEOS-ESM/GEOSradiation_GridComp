#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

baseline_dir = "/home/pchakrab/input/irrad/C12L181/after/new-style"
for state in ["import", "internal", "export"]:
    print("\nSTATE:", state)
    baseline = f"{baseline_dir}/IRRAD_{state}_after_runPhase1.nc"
    current = f"checkpoints/last/IRRAD_{state}.nc"
    print("BASELINE:", baseline)
    print("CURRENT:", current)
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if var in cur.variables:
                bas_var = bas[var][:]
                if np.ma.is_masked(bas_var):
                    fill_value = bas[var].getncattr("_FillValue")
                    bas_var = bas_var.filled(fill_value)
                cur_var = cur[var][0]
                diff = np.abs(bas_var - cur_var)
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<20} diff (min/max/norm): {np.min(diff)}, {np.max(diff)}, {diffnorm}")
