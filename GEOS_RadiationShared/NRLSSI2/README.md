=========================================================
Instructions on how to pre-process NRLSSI2 data for GEOS5
=========================================================
Peter Norris, April 2024

A. Introduction.
In the GEOS5 RRTMG and RRTMGP solar codes, the Solar Spectral Irradiance (SSI) at TOA is calculated from
NRLSSI2 Total Solar Irradiance (TSI) and Mg (facular) and SB (sunspot) indices. This data is currently
updated yearly, when the finalized yearly products become available, which is generally in late January
to late February of the following year. For example, the finalized 2023 data came out on 2024-02-23. This
README describes how to get the NRLSSI2 source data and how to pre-process it into a text file read by
the SolarGC.

B. Get the data for the NEW year.
The source data is obtained from https://www.ncei.noaa.gov/data/total-solar-irradiance/access/
Note: you can use wget https://... from discover.
There are two data files required for each new year:
1. from ancilliary-data/, e.g., tsi-ssi_v02r01_model-input-time-series_s18820101_e20180331_c20180418.txt
  !!! NOTE: Please read comments in docstring of Mg_SB_from_daily_file.py BEFORE downloading !!!
  In particular,
    -- the <filenames> list in that python code should ONLY be APPENDED to,
         which will ensure that any existing final data is not overwritten,
         thus ensuring historical reproducibility.
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
an environmental variable, just a shorthand. But whatever DATADIR is used must also be explicity set in
the "main" of the preprocessing driver: TSI_Mg_SB_merged_from_daily.py.

D. Preparing the preprocessor.
1. Add the filename from B.1. to the filenames list in Mg_SB_from_daily_file.py. Specifically APPEND
it to the END of the <filenames> list in the default argument of the __init__ of class Mg_SB_Daily
in Mg_SB_from_daily_file.py.
2. Set the OUTDIR in the "main" of TSI_Mg_SB_merged_from_daily.py to where you want to store the 
pre-processed file for the neaw year.

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
