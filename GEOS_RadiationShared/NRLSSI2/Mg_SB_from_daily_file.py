'''
Get Solar Mg and SB indices from NRLSSI2 daily auxiliary text file.

IMPORTANT: The filenames list below will be read in order, and any date that
already has final data will not be overwritten by data for that date in a later
file. For this reason, the filenames list should ONLY be APPENDED to, which
will ensure that any older final data is not overwritten, thus ensuring
historical reproducibility.

v03r00 support:
- v03 ancillary text columns:
  time (yyyy-MM-dd),Bolfac,Bolfac_unc,Bolfac_sm,Bolspot,Bolspot_unc,Ftcomp,Ftcomp_unc,source
- We DO NOT hard-code a conversion. Instead, when a v03 file is read we learn
  a linear mapping (least squares) from overlapping dates between:
    MgIndex_v02      = a_Mg * Bolfac + b_Mg
    SBindex_v02      = a_SB * Bolspot + b_SB
  using the Mg/SB already loaded from earlier v02 files in the filenames list.
- This keeps continuity across the version change and avoids manual tuning.
'''

import os
import numpy as np
from datetime import datetime

class Mg_SB_Daily(object):
    def __init__(self, DATADIR,
                 filenames=[
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20170630_c20170928.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20180331_c20180418.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20191231_c20200225.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20201231_c20210204.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20220630_c20220803.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20230331_c20230411.txt',
                     'tsi-ssi_v02r01_model-input-time-series_s18820101_e20231231_c20240221.txt',
                     # Append your v03 file(s) for 2024+ at the end of this list.
                     'tsi-ssi_v03r00_model-input-time-series_s18740509_e20250331_c20250723.txt',
                     'tsi-ssi_v03r00_model-input-time-series_s18740509_e20251231_c20260305.txt',
                 ],
                 verbose=True):

        self._data = {}         # yyyymmdd -> (rc, Mg, SB)
        self.nfinal = 0
        self.verbose = verbose

        # These will be learned the first time we see a v03 file
        self._v03_fitted = False
        self._a_Mg = None; self._b_Mg = None; self._r2_Mg = None
        self._a_SB = None; self._b_SB = None; self._r2_SB = None

        for filename in filenames:
            fullpath = os.path.join(DATADIR, filename)
            if self.verbose:
                print(f"[Mg_SB] Reading {fullpath}")
            if not os.path.exists(fullpath):
                if self.verbose:
                    print(f"[Mg_SB] File not found, skipping: {fullpath}")
                continue

            with open(fullpath, 'r') as f:
                lines = f.readlines()
            if not lines:
                continue

            header = lines[0].strip()
            is_v03 = ('Bolfac' in header and 'Bolspot' in header)

            if is_v03 and not self._v03_fitted:
                # Learn linear mapping using overlaps with whatever v02 data is already loaded
                self._fit_v03_to_v02(lines)

            # Now parse the file proper and insert/merge into _data
            if is_v03:
                self._ingest_v03(lines)
            else:
                self._ingest_v02(lines)

    # ---------- ingestion helpers ----------

    def _ingest_v02(self, lines):
        # header: time (yyyy-MM-dd), sunspot_darkening_function, facular_brightening_function, source
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.strip().split(',')]
            try:
                datestr, SB, Mg, status = parts[0], float(parts[1]), float(parts[2]), parts[3]
            except (IndexError, ValueError):
                continue
            yyyymmdd = datestr.replace('-', '')
            rc = self._status_to_rc(status)
            # final beats prelim
            if yyyymmdd in self._data and self._data[yyyymmdd][0] == 0:
                continue
            self._data[yyyymmdd] = (rc, Mg, SB)

    def _ingest_v03(self, lines):
        # header: time, Bolfac, Bolfac_unc, Bolfac_sm, Bolspot, Bolspot_unc, Ftcomp, Ftcomp_unc, source
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.strip().split(',')]
            try:
                datestr = parts[0]
                TF = float(parts[1])   # Bolfac (W/m^2-like units used by v03 table)
                TS = float(parts[4])   # Bolspot (W/m^2-like units used by v03 table)
                status = parts[-1] if parts else 'preliminary'
            except (IndexError, ValueError):
                continue

            yyyymmdd = datestr.replace('-', '')
            rc = self._status_to_rc(status)

            # Convert using learned linear maps; if not fitted, skip (or use identity as last resort)
            if self._v03_fitted:
                Mg = self._a_Mg * TF + self._b_Mg
                SB = self._a_SB * TS + self._b_SB
            else:
                # Should not happen if _fit_v03_to_v02 ran; fall back to passing through
                Mg = TF
                SB = TS
                if self.verbose:
                    print("[Mg_SB] WARNING: v03 coefficients not fitted; passing through raw values")

            # final beats prelim
            if yyyymmdd in self._data and self._data[yyyymmdd][0] == 0:
                continue
            self._data[yyyymmdd] = (rc, Mg, SB)

    def _status_to_rc(self, status):
        s = status.lower()
        if s == 'final':
            return 0
        if s.startswith('prelim'):
            return -1
        if s == 'preliminary':
            return -1
        raise RuntimeError(f'invalid status detected: {status}')

    # ---------- fitting: learn v03 -> legacy units ----------

    def _fit_v03_to_v02(self, v03_lines):
        # Build arrays for overlap where v02 data already exists
        x_TF = []
        x_TS = []
        y_Mg = []
        y_SB = []

        for line in v03_lines[1:]:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.strip().split(',')]
            try:
                datestr = parts[0]
                TF = float(parts[1])   # Bolfac
                TS = float(parts[4])   # Bolspot
            except (IndexError, ValueError):
                continue
            yyyymmdd = datestr.replace('-', '')
            if yyyymmdd in self._data and self._data[yyyymmdd][0] == 0:  # need finalized v02 values to anchor
                _, Mg_v02, SB_v02 = self._data[yyyymmdd]
                x_TF.append(TF); y_Mg.append(Mg_v02)
                x_TS.append(TS); y_SB.append(SB_v02)

        x_TF = np.asarray(x_TF); y_Mg = np.asarray(y_Mg)
        x_TS = np.asarray(x_TS); y_SB = np.asarray(y_SB)

        if len(x_TF) < 30 or len(x_TS) < 30:  # need a reasonable number of days
            # Not enough overlap; do a simple proportional guess based on medians to avoid blow-ups
            if self.verbose:
                print(f"[Mg_SB] Not enough overlap to fit (Mg n={len(x_TF)}, SB n={len(x_TS)}); using robust ratios")
            self._a_Mg, self._b_Mg, self._r2_Mg = self._robust_ratio(x_TF, y_Mg)
            self._a_SB, self._b_SB, self._r2_SB = self._robust_ratio(x_TS, y_SB)
        else:
            # Ordinary least squares with intercept
            self._a_Mg, self._b_Mg, self._r2_Mg = self._ols_with_r2(x_TF, y_Mg)
            self._a_SB, self._b_SB, self._r2_SB = self._ols_with_r2(x_TS, y_SB)

        self._v03_fitted = True
        if self.verbose:
            print(f"[Mg_SB] Fitted v03->legacy:")
            print(f"  Mg:  Mg_v02 ≈ {self._a_Mg:.6f} * Bolfac + {self._b_Mg:.6f}  (R^2={self._r2_Mg:.4f}, n={len(x_TF)})")
            print(f"  SB:  SB_v02 ≈ {self._a_SB:.6f} * Bolspot + {self._b_SB:.6f}  (R^2={self._r2_SB:.4f}, n={len(x_TS)})")

    @staticmethod
    def _ols_with_r2(x, y):
        A = np.vstack([x, np.ones_like(x)]).T
        coef, resid, _, _ = np.linalg.lstsq(A, y, rcond=None)
        a, b = coef
        # R^2
        yhat = a * x + b
        ss_res = np.sum((y - yhat)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1.0 - (ss_res / ss_tot if ss_tot > 0 else 0.0)
        return a, b, r2

    @staticmethod
    def _robust_ratio(x, y):
        # fallback: scale and small intercept from medians
        if len(x) == 0:
            return 1.0, 0.0, 0.0
        a = np.median(y) / max(np.median(x), 1e-9)
        b = np.median(y) - a * np.median(x)
        # no real R^2 here
        return a, b, 0.0

    # ---------- dict-like / legacy helpers ----------

    def keys(self):
        return self._data.keys()

    def __getitem__(self, yyyymmdd):
        return self._data[yyyymmdd]

    def getday(self, yyyymmdd):
        if hasattr(yyyymmdd, "strftime"):  # datetime/date
            yyyymmdd = yyyymmdd.strftime("%Y%m%d")
        return self._data[yyyymmdd]

