'''
Get Solar Mg and SB indices from NRLSSI2 daily auxilliary text file
IMPORTANT: The filenames list below will be read in order, and any date that
already has final data will not be overwriten by data for that date in a later
file. For this reason, the filenames list should ONLY be APPENDED to, which
will ensure that any older final data is not overwritten, thus ensuring
historical reproducability.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from calendar import monthrange
from datetime import datetime, timedelta

# because datetime.strftime doesnt work before 1900
def _yyyymmdd(t):
  return('%04d%02d%02d' % (t.year, t.month, t.day))

class Mg_SB_Daily:
  ''' Daily time series for Solar Mg and SB indices from NRLSSI2 '''

  def __init__(self,
      filenames=[
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20170630_c20170928.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20180331_c20180418.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20191231_c20200225.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20201231_c20210204.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20220630_c20220803.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20230331_c20230411.txt',
        'data/tsi-ssi_v02r01_model-input-time-series_s18820101_e20231231_c20240221.txt',
        ],
      verbose=False):

    self.filenames = filenames
    self.verbose   = verbose

    NB = os.environ['NOBACKUP'] 

    if verbose:
      print('reading Mg and SB indicies from', filenames)

    # read the files
    self._data = {}
    self.nfinal = 0
    first = True
    for filename in filenames:
      with open(os.sep.join((NB,'NRLSSI2',filename))) as f:
        lines = f.readlines()

      # process the data, ignoring the header (top line)
      for line in lines[1:]:
        datestr, SB, Mg, status = line.strip().split(',')
        # yyyy-mm-dd to yyyymmdd
        yyyymmdd = datestr.replace('-','')
        # make data real
        SB = float(SB)
        Mg = float(Mg)
        # validate status to rc
        if status == 'final':
          rc = 0
        elif status == 'prelim':
          rc = -1
        else:
          raise RuntimeError('invalid status detected: %s' % status)
        # if a date is already finalized, ignore further records
        #   for that date to ensure historical reproducability
        if yyyymmdd in self._data:
          rc_existing = self._data[yyyymmdd][0]
          if rc_existing == 0: continue
        # load data dictionary
        self._data[yyyymmdd] = (rc, Mg, SB)
        # keep track of min and max final times
        if rc == 0:
          self.nfinal += 1
          d = datetime.strptime(yyyymmdd,'%Y%m%d').date()
          if first:
            self.date_min_final = d
            self.date_max_final = d
            first = False
          else:
            if d < self.date_min_final: self.date_min_final = d
            if d > self.date_max_final: self.date_max_final = d

  def keys(self):
    return self._data.keys()

  def getday(self, yyyymmdd):
    '''
    get daily Mg and SB values for daily string yyyymmdd
    returns (rc, Mg, SB)
    rc = 0 (final), -1 (prelim), -2 (unavailable)
    '''
    if yyyymmdd not in self._data:
      return -2, np.nan, np.nan
    else:
      return self._data[yyyymmdd]

  def gettime(self, t):
    '''
    get time-interpolated Mg and SB values for datetime t
    returns (rc, Mg, SB)
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
      Mg = vnoon[1] * (1-fother) + vother[1] * fother
      SB = vnoon[2] * (1-fother) + vother[2] * fother
      rc = -1 if vnoon[0] == -1 or vother[0] == -1 else 0
      return (rc, Mg, SB)

  def _plot_all_final(self):
    spyear = 365.25 * 86400.
    d0 = self.date_min_final; d = d0
    dy = []; Mg = []; SB = []
    while (d <= self.date_max_final):
      v = self.getday(_yyyymmdd(d))
      if v[0] == 0:
        # only plot finals
        dy.append((d-d0).total_seconds()/spyear)
        Mg.append(v[1])
        SB.append(v[2])
      d += timedelta(days=1)
    plt.subplot(211); plt.plot(dy,Mg,'b-'); plt.ylabel('Mg')
    plt.subplot(212); plt.plot(dy,SB,'b-'); plt.ylabel('SB')
    plt.xlabel('years since %s' % self.date_min_final)

if __name__ == '__main__':

  MgSB = Mg_SB_Daily()
# print('16000101', MgSB.getday('16000101'))
# print('20161130', MgSB.getday('20161130'))
# print('20161201', MgSB.getday('20161201'))
# print('20161202', MgSB.getday('20161202'))
# t = datetime.strptime('2016-12-01 10:00:00','%Y-%m-%d %H:%M:%S')
# print(t, MgSB.gettime(t))

  plt.figure()
  MgSB._plot_all_final()
# plt.savefig(os.sep.join(('gx','MgSB_plot_all_final.png')),
#   pad_inches=0.33,bbox_inches='tight',dpi=100)
  plt.show()

