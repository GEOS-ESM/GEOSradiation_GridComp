=========================================================
Instructions on how to pre-process NRLSSI2 data for GEOS5
=========================================================
Peter Norris, April 2024
Matt Thompson, April 2025

A. Introduction.
In the GEOS5 RRTMG and RRTMGP solar codes, the Solar Spectral Irradiance (SSI) at TOA is calculated from
NRLSSI2 Total Solar Irradiance (TSI) and Mg (facular) and SB (sunspot) indices. This data is currently
updated yearly, when the finalized yearly products become available, which is generally in late January
to late February of the following year. For example, the finalized 2023 data came out on 2024-02-23. This
README describes how to get the NRLSSI2 source data and how to pre-process it into a text file read by
the SolarGC.

NOTE: Starting with 2024 data, NRL released version v03r00 of the ancillary text files, which use a
different column format than the legacy v02r01 files. See section I for details on how this is handled.

B. Get the data for the NEW year.
The source data is obtained from https://www.ncei.noaa.gov/data/total-solar-irradiance/access/
Note: you can use wget https://... from discover.
There are two data files required for each new year:
1. from ancillary-data/, e.g.:
   - v02r01 (through 2023): tsi-ssi_v02r01_model-input-time-series_s18820101_e20231231_c20240221.txt
   - v03r00 (2024 onward):  tsi-ssi_v03r00_model-input-time-series_s18740509_e20251231_c20260305.txt
  !!! NOTE: Please read comments in docstring of Mg_SB_from_daily_file.py BEFORE downloading !!!
  In particular,
    -- the <filenames> list in that python code should ONLY be APPENDED to,
         which will ensure that any existing final data is not overwritten,
         thus ensuring historical reproducibility.
    -- v02r01 files are no longer updated past 2023. For 2024 onward, download v03r00 files.
    -- v02r01 files MUST appear before any v03r00 files in the filenames list (see section I).
2. from daily/, e.g., tsi_v02r01_daily_s20170101_e20171231_c20180227.nc
  !!! NOTE: Please read comments in docstring of TSI_from_daily_files.py BEFORE downloading !!!
  In particular,
    -- there should be no files with overlapping time periods, and
    -- you should NEVER overwrite any existing non-preliminary file,
         since this may cause historical non-reproducibility.
NOTE: Only two files are needed for the new year. Existing data has already been obtained
for previous years. There is no need to get historical data for past years. Even though the
file in part 1 above contains historical data, only the data for the new year inside it will
be used if you follow these instructions, in order to maintain historical reproducibility.

C. Decide where to put the source data.
The source data is put in the directory DATADIR, where DATADIR is currently
   /discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/solar/NRLSSI2/data
This directory is managed by the GMAO Operations Team, so if you is not they, you will need to download
the data somewhere else first and then work with "Ops" to add it to that directory. DATADIR here is NOT
an environmental variable, just a shorthand. But whatever DATADIR is used must also be explicitly set in
the "main" of the preprocessing driver: TSI_Mg_SB_merged_from_daily.py.

D. Preparing the preprocessor.
1. Add the filename from B.1. to the filenames list in Mg_SB_from_daily_file.py. Specifically APPEND
it to the END of the <filenames> list in the default argument of the __init__ of class Mg_SB_Daily
in Mg_SB_from_daily_file.py.
IMPORTANT: All v02r01 files must appear before any v03r00 files in the list. The v03 conversion
relies on learning a mapping from the overlap period, which requires v02 data to already be loaded.
2. Set the DATADIR and OUTDIR in the "main" of TSI_Mg_SB_merged_from_daily.py to where you want to
read the source data and store the pre-processed file for the new year.

E. Running the preprocessor.
python TSI_Mg_SB_merged_from_daily.py
After you exit the program, rename the NRLSSI2.VYYYY.txt in OUTDIR to the correct year.

F. Storing the Output file
Get "Ops" to put the new NRLSSI2.VYYYY.txt from your OUTDIR into their directory:
  /discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/solar

G. Symlinks
There is a symlink in later "Ops solar" directory named NRLSSI2.txt.
It will typically be pointed to the file NRLSSI2.vYYYY.txt for the new year, but that decision must
be approved by Ops and by the Modelling Team and Data Assimilation Team heads.

H. Using the data.
The GEOS5 run should have the following lines in its ASGCM.rc:
USE_NRLSSI2: .TRUE.
SOLAR_CYCLE_FILE_NAME: ExtData/g5gcm/solar/NRLSSI2.txt
(Note: These will only be used if the solar code being run is RRTMG or RRTMGP. Chou-Suarez runs do
not use NRLSSI2 data).
PS: You are free to change the use of the symlink NRLSSI2.txt to any specific NRLSSI2.vYYYY.txt 
if you want reproducibility to a specific historical run.

I. Handling NRLSSI2 v03r00 Input Data (conversion to legacy Mg/SB format).

Background:
The GEOS radiation code expects solar variability inputs in the legacy NRLSSI2 v02r01 format:
  - MgIndex: facular brightening proxy (unitless)
  - SBindex: sunspot darkening proxy (in millionths of a hemisphere)

The newer NRLSSI2 v03r00 ancillary text file instead provides:
  - Bolfac:  bolometric facular contribution
  - Bolspot: bolometric sunspot contribution

These are related to the older indices but are not provided in the same scaling as the legacy
GEOS input format. The v03r00 file also extends the historical record back further, to 1874-05-09
(compared to 1882-01-01 for v02r01), though only the 2024+ data is new.

Approach:
Rather than modifying the upstream NRLSSI2 files or the GEOS radiation code, the conversion is
done in Mg_SB_from_daily_file.py at read time using an empirical linear mapping:

  1. Read the historical v02r01 files first (these must come first in the filenames list).
  2. When a v03r00 file is encountered, use the overlapping dates between v02 and v03.
  3. Fit linear relationships (least squares) that map the v03 quantities onto the legacy format:
       MgIndex_v02 = a_Mg * Bolfac  + b_Mg
       SBindex_v02 = a_SB * Bolspot + b_SB
  4. Apply those fitted relationships to convert all v03-only dates into legacy Mg/SB values.

Current fitted coefficients (computed from the 1882-2023 overlap period):
  Mg:  Mg_v02 ≈ 0.007506 * Bolfac  + 0.150331     (R² ≈ 0.905, n=51864)
  SB:  SB_v02 ≈ 1941.389784 * Bolspot + 42.906433  (R² ≈ 0.994, n=51864)

Interpretation:
  - The sunspot (SB) mapping is very tight (R² ≈ 0.994): SBindex is well represented by a
    linear scaling of Bolspot.
  - The facular (Mg) mapping has more scatter (R² ≈ 0.905): MgIndex follows a clear linear
    relationship with Bolfac, but post-2023 Mg values carry more uncertainty than pre-2024 values.

The coefficients are fitted automatically at runtime and printed to stdout. If NOAA/NCEI revises
either the v02 or v03 datasets, the fitted coefficients will change accordingly.

The final GEOS-facing output file remains in the same legacy format as before:
  # yyyy doy TSI:W/m2 MgIndex   SBindex

Notes and caveats:
  - This conversion is empirical, not a first-principles physical derivation.
  - Post-2023 Mg/SB values are approximations derived from the v03 Bolfac/Bolspot quantities.
  - The diagnostic plot produced by the preprocessor shades the post-2024 region to indicate
    where v03-derived mapped values begin.
  - If NOAA/NCEI revises either the v02 or v03 datasets, the fitted coefficients should be
    recomputed and verified.
