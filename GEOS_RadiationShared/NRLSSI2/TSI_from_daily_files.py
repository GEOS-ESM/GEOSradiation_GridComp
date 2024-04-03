'''
Get TSI from daily NRLSSI2 netCDF files
============================================================
IMPORTANT RULES THAT GOVERN data/tsi_*_daily_*.nc FILES:
============================================================
In order to maintain historical reproducabilty and not fail
the error checking below, please ensure the following for
data/tsi_*_daily_*.nc files:
1. files specified in the file_template argument must
obey 'data/tsi_*_daily_s????????_e????????_c????????.nc'
2. the first '*' is the version e.g., 'v02r01', with an
optional trailing '-preliminary'
3. The starting date s???????? should be of the form
     sYYYYMMDD, where MMDD will typically be 0101.
4. The ending date e???????? should be of the form
     eYYYYMMDD, where MMDD will typically be 1231,
     but may be before then for preliminary files.
3. The creation date c???????? should be of the form
     cYYYYMMDD. 
5. It is up to the downloader and maintainer of these files
to ENSURE that
  a. There should be no files with overlapping time periods
     (For example, multiple versions / creation dates for
     the same start and end date are not allowed. You must
     keep only one, and probably the ORIGINAL(!) one, as
     discussed in point (c.) below.)
  b. You should NEVER overwrite an existing non-preliminary
     file, since this may lose historical reproducability.
     (Preliminary files may be overwritten with non-prelim-
     inary ones, or with preliminary files of a later
     creation date.)
  c. Even if a new version comes along, that should only be
     used for new years (or else you should maintain data
     directories for different versions). Overwriting an 
     older version with a new one for an existing period in
     this data directory, will destroy historical reproduc-
     ability.
============================================================
'''

import os
import sys
import numpy as np
from glob import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from calendar import monthrange
from datetime import datetime, timedelta
from operator import itemgetter

def _yyyymmdd(t):
  ''' because datetime.strftime doesnt work before 1900 '''
  return('%04d%02d%02d' % (t.year, t.month, t.day))

class TSI_Daily:
  ''' Daily time series for NRLSSI2 TSI '''

  def __init__(self,
      file_template='data/tsi_*_daily_s????????_e????????_c????????.nc',
      verbose=False):

    self.file_template = file_template
    self.verbose  = verbose

    if verbose:
      print('reading TSI from', file_template)

    # ==============
    # Validate files
    # ==============

    NB = os.environ['NOBACKUP']

    # collect templated files and validate
    files = []
    for path in glob(os.sep.join((NB,'NRLSSI2',file_template))):
      filename = os.path.basename(path)
      # validate filename
      v = filename.split('_')
      if v[0] != 'tsi' or v[2] != 'daily':
        raise RuntimeError('%s not an NRLSSI2 tsi daily file' % file)
      # version preliminary?
      version = v[1]
      if '-preliminary' in version:
        rc = -1
      else:
        rc = 0
      # extract inclusive date range of file
      ds = datetime.strptime(v[3],'s%Y%m%d').date()
      de = datetime.strptime(v[4],'e%Y%m%d').date()
      if ds > de:
        raise RuntimeError('bad timestamp ordering in %s' % filename)
      # add to files to process
      files.append((ds, de, rc, version, path))

    if len(files) == 0:
      raise RuntimeError('no TSI files to process')

    
    # ========================
    # Fail if any time overlap
    # ========================

    # sort the files by start date
    files.sort(key=itemgetter(0))

    # detect any overlapping date ranges
    if len(files) > 1:
      ds0, de0, r0, v0, f0 = files[0]
      for ds1, de1, r1, v1, f1 in files[1:]:
        if ds0 == ds1:
          # same start so overlapping
          raise RuntimeError('date overlap: %s, %s' % (f0, f1))
        else:
          # now ds0 < ds1 because of sort
          if de0 >= ds1:
            raise RuntimeError('date overlap: %s, %s' % (f0, f1))
        ds0, de0, r0, v0, f0 = ds1, de1, r1, v1, f1

    # =================
    # Process the files
    # =================

    self._data = {}
    self.nfinal = 0
    first = True
    for ds, de, rc, version, path in files:

      # open file
      if self.verbose:
        print('processing path %s' % path)
      hf = Dataset(path,'r')

      # extract dates
      htime = hf.variables['time']
      tref = datetime.strptime(htime.units,'days since %Y-%m-%d %H:%M:%S')
      dates = [(tref+timedelta(days=float(days))).date() for days in htime[:]]

      # read TSI, TSI_UNC
      TSIs     = hf.variables['TSI'    ][:].astype(float)
      TSI_UNCs = hf.variables['TSI_UNC'][:].astype(float)

      # close file
      hf.close()
      if not self.verbose:
        sys.stdout.write('.')
        sys.stdout.flush()

      # process file
      for d, TSI, TSI_UNC in zip(dates, TSIs, TSI_UNCs):

        # load data dictionary
        yyyymmdd = _yyyymmdd(d)
        if yyyymmdd in self._data:
          raise RuntimeError('date duplicates detected for %s' % yyyymmdd)
        else:
          self._data[yyyymmdd] = (rc, TSI, TSI_UNC)

        # keep track of min and max final times
        if rc == 0:
          self.nfinal += 1
          if first:
            self.date_min_final = d
            self.date_max_final = d
            first = False
          else:
            if d < self.date_min_final: self.date_min_final = d
            if d > self.date_max_final: self.date_max_final = d

    if not self.verbose:
      sys.stdout.write('\n')
      sys.stdout.flush()

  def viewkeys(self):
    return self._data.viewkeys()

  def getday(self, yyyymmdd):
    '''
    get daily TSI values for daily string yyyymmdd
    returns (rc, TSI, TSI_UNC)
    rc = 0 (final), -1 (prelim), -2 (unavailable)
    '''
    if yyyymmdd not in self._data:
      return -2, np.nan, np.nan
    else:
      return self._data[yyyymmdd]

  def gettime(self, t):
    '''
    get time-interpolated TSI values for datetime t
    returns (rc, TSI, TSI_UNC)
    rc = 0 (final), -1 (prelim), -2 (unavailable)
    '''
    
    # assume daily average valid at noon
    tnoon = datetime(t.year,t.month,t.day,12)
    vnoon = self.getday(_yyyymmdd(tnoon))
    if t == tnoon: return vnoon
    
    # other noon bracketing t
    tother = tnoon + timedelta(days=(-1 if t < tnoon else +1))
    vother = self.getday(_yyyymmdd(tother))

    # fraction that the other daily average contributes
    fother = abs((t-tnoon).total_seconds()) / 86400.

    # only interpolate if both days exist
    if vnoon[0] == -2 or vother[0] == -2:
      return (-2, np.nan, np.nan)
    else:
      # now both days available
      TSI     = vnoon[1] * (1-fother) + vother[1] * fother
      TSI_UNC = vnoon[2] * (1-fother) + vother[2] * fother
      rc = -1 if vnoon[0] == -1 or vother[0] == -1 else 0
      return (rc, TSI, TSI_UNC)

  def _plot_all_final(self):
    spyear = 365.25 * 86400.
    d0 = self.date_min_final; d = d0
    dy = []; TSI = []; TSI_UNC = []
    while (d <= self.date_max_final):
      v = self.getday(_yyyymmdd(d))
      if v[0] == 0:
        # only plot finals
        dy.append((d-d0).total_seconds()/spyear)
        TSI.append(v[1])
        TSI_UNC.append(v[2])
      d += timedelta(days=1)
    plt.subplot(211); plt.plot(dy,TSI,'b-'); plt.ylabel('TSI [W/m2]')
    plt.title(r'%s'%self.file_template,fontsize=10)
    plt.subplot(212); plt.plot(dy,TSI_UNC,'b-'); plt.ylabel('TSI_UNC [W/m2]')
    plt.xlabel('years since %s' % self.date_min_final)

if __name__ == '__main__':

  TSI = TSI_Daily()
  print('16000101', TSI.getday('16000101'))
  print('20161130', TSI.getday('20161130'))
  print('20161201', TSI.getday('20161201'))
  print('20161202', TSI.getday('20161202'))
  t = datetime.strptime('2016-12-01 10:00:00','%Y-%m-%d %H:%M:%S')
  print(t, TSI.gettime(t))

  plt.figure()
  TSI._plot_all_final()
  plt.savefig(os.sep.join(('gx','TSI_plot_all_final.png')),
    pad_inches=0.33,bbox_inches='tight',dpi=100)
  plt.show()

