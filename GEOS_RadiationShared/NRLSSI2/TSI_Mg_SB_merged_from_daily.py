'''
Get TSI and Solar Mg and SB indices from NRLSSI2 daily data
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from calendar import monthrange
from datetime import datetime, timedelta
from TSI_from_daily_files import TSI_Daily
from Mg_SB_from_daily_file import Mg_SB_Daily

# because datetime.strftime doesnt work before 1900
def _yyyymmdd(t):
  return('%04d%02d%02d' % (t.year, t.month, t.day))

class TSI_Mg_SB_Merged_Daily:
  ''' Daily time series for TSI and Solar Mg and SB indices from NRLSSI2 '''

  def __init__(self, verbose=False):

    self.verbose  = verbose

    # get each data set
    zM = Mg_SB_Daily(verbose=verbose)
    zT =   TSI_Daily(verbose=verbose)

    # form the INTERSECTION:
    # (each set is indexed by yyyymmdd, unique in each set)
    self._data = {}
    self.nfinal = 0
    first = True
    for yyyymmdd in (zM.keys() & zT.keys()):

        # load data dictionary
        rc1, Mg, SB       = zM.getday(yyyymmdd)
        rc2, TSI, TSI_UNC = zT.getday(yyyymmdd)
        rc = -1 if rc1 == -1 or rc2 == -1 else 0
        self._data[yyyymmdd] = (rc, TSI, Mg, SB, TSI_UNC)

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

  def getday(self, yyyymmdd):
    '''
    get daily TSI, Mg and SB values for daily string yyyymmdd
    returns (rc, TSI, Mg, SB, TSI_UNC)
    rc = 0 (final), -1 (prelim), -2 (unavailable)
    '''
    if yyyymmdd not in self._data:
      return -2, np.nan, np.nan, np.nan
    else:
      return self._data[yyyymmdd]

  def gettime(self, t):
    '''
    get time-interpolated TSI, Mg, and SB values for datetime t
    returns (rc, TSI, Mg, SB, TSI_UNC)
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
      return (-2, np.nan, np.nan, np.nan)
    else:
      # now both days available
      TSI     = vnoon[1] * (1-fother) + vother[1] * fother
      Mg      = vnoon[2] * (1-fother) + vother[2] * fother
      SB      = vnoon[3] * (1-fother) + vother[3] * fother
      TSI_UNC = vnoon[4] * (1-fother) + vother[4] * fother
      rc = -1 if vnoon[0] == -1 or vother[0] == -1 else 0
      return (rc, TSI, Mg, SB, TSI_UNC)

  def final_date_range(self,
      yyyymmdd0, yyyymmdd1,
      no_gaps=True):
    ''' final only lists over date range '''

    d0 = datetime.strptime(yyyymmdd0,'%Y%m%d').date()
    d1 = datetime.strptime(yyyymmdd1,'%Y%m%d').date()
    d = d0; dates = []; TSI = []; Mg = []; SB = []
    while (d <= d1):
      v = self.getday(_yyyymmdd(d))
      if v[0] == 0:
        # only include finals
        dates.append(d)
        TSI.append(v[1])
        Mg.append( v[2])
        SB.append( v[3])
      elif no_gaps:
        raise RuntimeError('no_gaps: intervening non-final date found')
      d += timedelta(days=1)
    return dates, TSI, Mg, SB

  def _plot_all_final(self):
    dates, TSI, Mg, SB = self.final_date_range(
      _yyyymmdd(self.date_min_final),
      _yyyymmdd(self.date_max_final),
      no_gaps=False)
    spyear = 365.25 * 86400.
    dy = [(d-self.date_min_final).total_seconds()/spyear for d in dates]
    plt.subplot(311); plt.plot(dy,TSI,'b-'); plt.ylabel('TSI')
    plt.subplot(312); plt.plot(dy, Mg,'b-'); plt.ylabel('Mg')
    plt.subplot(313); plt.plot(dy, SB,'b-'); plt.ylabel('SB')
    plt.xlabel('years since %s' % self.date_min_final)

  def output_final_textfile(self, filename):
    f = open(filename,'w')
    f.write('# NRLSSI2 daily input\n')
    f.write('# treat daily values as valid at 12:00 GMT\n')
    f.write('# yyyy doy TSI:W/m2 MgIndex   SBindex\n')
    for d, TSI, Mg, SB in zip(*self.final_date_range(
        _yyyymmdd(self.date_min_final), _yyyymmdd(self.date_max_final))):
      f.write('  %04d %03d %8.3f %8.6f %9.4f\n' % (
        d.year, d.timetuple().tm_yday, TSI, Mg, SB))
    f.close()

if __name__ == '__main__':

  NB    = os.environ['NOBACKUP']

  z = TSI_Mg_SB_Merged_Daily()
# print('16000101', z.getday('16000101'))
# print('20161130', z.getday('20161130'))
# print('20161201', z.getday('20161201'))
# print('20161202', z.getday('20161202'))
# t = datetime.strptime('2016-12-01 10:00:00','%Y-%m-%d %H:%M:%S')
# print(t, z.gettime(t))

  z.output_final_textfile(os.sep.join((NB,'NRLSSI2','output','NRLSSI2.txt')))

  plt.figure(figsize=(6,12))
  z._plot_all_final()
# plt.savefig(os.sep.join(('gx','TSInMgSB_plot_all_final.png')),
#   pad_inches=0.33,bbox_inches='tight',dpi=100)
  plt.show()

