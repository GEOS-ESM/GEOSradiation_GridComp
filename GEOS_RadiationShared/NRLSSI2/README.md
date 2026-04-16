# Pre-processing NRLSSI2 Data for GEOS5

*Peter Norris, April 2024 — Matt Thompson, April 2025*

> **Note:** Starting with 2024 data, NRL released version v03r00 of the ancillary text files,
> which use a different column format than the legacy v02r01 files. See [Section I](#i-handling-nrlssi2-v03r00-input-data) for details.

---

## A. Introduction

In the GEOS5 RRTMG and RRTMGP solar codes, the Solar Spectral Irradiance (SSI) at TOA is
calculated from NRLSSI2 Total Solar Irradiance (TSI) and Mg (facular) and SB (sunspot) indices.
This data is currently updated yearly, when the finalized yearly products become available, which
is generally in late January to late February of the following year. For example, the finalized
2023 data came out on 2024-02-23. This README describes how to get the NRLSSI2 source data and
how to pre-process it into a text file read by the SolarGC.

---

## B. Get the Data for the New Year

The source data is obtained from <https://www.ncei.noaa.gov/data/total-solar-irradiance/access/>.
You can use `wget https://...` from Discover.

Two data files are required for each new year:

### B.1. Ancillary file (from `ancillary-data/`)

> **Before downloading, read the docstring in `Mg_SB_from_daily_file.py`.**

Example filenames:
- v02r01 (through 2023): `tsi-ssi_v02r01_model-input-time-series_s18820101_e20231231_c20240221.txt`
- v03r00 (2024 onward):  `tsi-ssi_v03r00_model-input-time-series_s18740509_e20251231_c20260305.txt`

Key rules:
- The `filenames` list in the code must **only ever be appended to** — never reordered or
  shortened — to ensure existing final data is never overwritten (historical reproducibility).
- v02r01 files are **no longer updated past 2023**. For 2024 onward, download v03r00 files.
- v02r01 files **must appear before** any v03r00 files in the list (see [Section I](#i-handling-nrlssi2-v03r00-input-data)).

### B.2. Daily TSI file (from `daily/`)

> **Before downloading, read the docstring in `TSI_from_daily_files.py`.**

Example filename: `tsi_v02r01_daily_s20170101_e20171231_c20180227.nc`

Key rules:
- There must be **no files with overlapping time periods**.
- **Never overwrite** any existing non-preliminary file, as this may cause historical
  non-reproducibility.

> **Note:** Only two files are needed for the new year. Existing data has already been obtained
> for previous years. Even though the ancillary file contains historical data, only the data for
> the new year will be used if you follow these instructions.

---

## C. Decide Where to Put the Source Data

The source data goes in `DATADIR`, which is currently:

```
/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/solar/NRLSSI2/data
```

This directory is managed by the GMAO Operations Team. If you are not they, download the data
somewhere else first and then work with "Ops" to add it to that directory. `DATADIR` is not an
environment variable — just a shorthand. Whatever path is used must also be explicitly set in
the `__main__` block of `TSI_Mg_SB_merged_from_daily.py`.

---

## D. Preparing the Preprocessor

1. **Append** the new ancillary filename (from B.1) to the end of the `filenames` list in the
   `__init__` of class `Mg_SB_Daily` in `Mg_SB_from_daily_file.py`.

   > **Important:** All v02r01 files must appear before any v03r00 files in the list. The v03
   > conversion relies on learning a mapping from the overlap period, which requires v02 data
   > to already be loaded.

2. Set `DATADIR` and `OUTDIR` in the `__main__` block of `TSI_Mg_SB_merged_from_daily.py` to
   where you want to read the source data and write the pre-processed output.

---

## E. Running the Preprocessor

```bash
python TSI_Mg_SB_merged_from_daily.py
```

After the script completes, rename `NRLSSI2.vYYYY.txt` in `OUTDIR` to reflect the correct year.

---

## F. Storing the Output File

Ask "Ops" to place the new `NRLSSI2.vYYYY.txt` from your `OUTDIR` into:

```
/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/solar
```

---

## G. Symlinks

There is a symlink in the "Ops solar" directory named `NRLSSI2.txt`. It will typically be
pointed to `NRLSSI2.vYYYY.txt` for the new year, but that decision must be approved by Ops
and by the Modelling Team and Data Assimilation Team heads.

---

## H. Using the Data

The GEOS5 run should have the following lines in its `ASGCM.rc`:

```
USE_NRLSSI2: .TRUE.
SOLAR_CYCLE_FILE_NAME: ExtData/g5gcm/solar/NRLSSI2.txt
```

> These settings are only used if the solar code is RRTMG or RRTMGP. Chou-Suarez runs do not
> use NRLSSI2 data.

You may replace the symlink `NRLSSI2.txt` with a specific `NRLSSI2.vYYYY.txt` if you need
reproducibility to a particular historical run.

---

## I. Handling NRLSSI2 v03r00 Input Data

### Background

The GEOS radiation code expects solar variability inputs in the legacy NRLSSI2 v02r01 format:

| Variable   | Description                                          |
|------------|------------------------------------------------------|
| `MgIndex`  | Facular brightening proxy (unitless)                 |
| `SBindex`  | Sunspot darkening proxy (millionths of a hemisphere) |

The newer NRLSSI2 v03r00 ancillary text file instead provides:

| Variable   | Description                          |
|------------|--------------------------------------|
| `Bolfac`   | Bolometric facular contribution      |
| `Bolspot`  | Bolometric sunspot contribution      |

These are physically related to the older indices but use different scaling. The v03r00 file
also extends the historical record back to 1874-05-09 (vs. 1882-01-01 for v02r01), though
only the 2024+ data is new.

### Approach

Rather than modifying the upstream NRLSSI2 files or the GEOS radiation code, the conversion
is performed in `Mg_SB_from_daily_file.py` at read time using an empirical linear mapping:

1. Read the historical v02r01 files first (they must come first in the `filenames` list).
2. When a v03r00 file is encountered, identify the overlapping dates between v02 and v03.
3. Fit linear least-squares relationships mapping v03 quantities to the legacy format:
   ```
   MgIndex_v02 = a_Mg * Bolfac  + b_Mg
   SBindex_v02 = a_SB * Bolspot + b_SB
   ```
4. Apply those fitted relationships to convert all v03-only dates into legacy Mg/SB values.

### Fitted Coefficients

Computed from the 1882–2023 overlap period:

```
Mg:  Mg_v02 ≈ 0.007506 * Bolfac + 0.150331        (R² = 0.9049, n=51864)
SB:  SB_v02 ≈ 1941.389784 * Bolspot + 42.906433    (R² = 0.9940, n=51864)
```

The coefficients are fitted automatically at runtime and printed to stdout. They will update
if NOAA/NCEI revises the underlying datasets.

### Interpretation

- **SB mapping** (R² = 0.994): Very tight fit. `SBindex` is well represented by a linear
  scaling of `Bolspot`.
- **Mg mapping** (R² = 0.905): Clear linear relationship, but with more scatter. Post-2023
  `MgIndex` values carry more uncertainty than pre-2024 values.

### Output Format

The final GEOS-facing output file remains in the same legacy format as before:

```
# yyyy doy TSI:W/m2 MgIndex   SBindex
```

### Notes and Caveats

- This conversion is **empirical**, not a first-principles physical derivation.
- Post-2023 Mg/SB values are approximations derived from v03 `Bolfac`/`Bolspot` quantities.
- The diagnostic plot produced by the preprocessor shades the post-2024 region to indicate
  where v03-derived mapped values begin.
- If NOAA/NCEI revises either the v02 or v03 datasets, the fitted coefficients should be
  recomputed and verified.
