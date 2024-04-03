get data files from https://www.ncei.noaa.gov/data/total-solar-irradiance/access/
Note: can use wget https://... from discover
1. from ancilliary-data/, e.g., tsi-ssi_v02r01_model-input-time-series_s18820101_e20180331_c20180418.txt
  !!! NOTE: Please read comments in docstring of Mg_SB_from_daily_file.py BEFORE downloading !!!
  In particular,
    -- the filenames list in that file should ONLY be APPENDED to,
         which will ensure that any existing final data is not overwritten,
         thus ensuring historical reproduceability.
2. from daily/, e.g., tsi_v02r01_daily_s20170101_e20171231_c20180227.nc
  !!! NOTE: Please read comments in docstring of TSI_from_daily_files.py BEFORE downloading !!!
  In particular,
    -- there should be no files with overlapping time periods, and
    -- you should NEVER overwrite any existing non-preliminary file,
         since this may cause historical non-reproduceability.
