"""
Get TSI and Solar Mg and SB indices from NRLSSI2 daily data
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from calendar import monthrange
from datetime import datetime, timedelta
from TSI_from_daily_files import TSI_Daily
from Mg_SB_from_daily_file import Mg_SB_Daily


# because datetime.strftime doesnt work before 1900
def _yyyymmdd(t):
    return "%04d%02d%02d" % (t.year, t.month, t.day)


class TSI_Mg_SB_Merged_Daily:
    """Daily time series for TSI and Solar Mg and SB indices from NRLSSI2"""

    def __init__(self, DATADIR, verbose=True):

        self.verbose = verbose

        # get each data set
        zM = Mg_SB_Daily(DATADIR, verbose=verbose)
        zT = TSI_Daily(DATADIR, verbose=verbose)

        # form the INTERSECTION:
        # (each set is indexed by yyyymmdd, unique in each set)
        self._data = {}
        self.nfinal = 0
        first = True
        for yyyymmdd in zM.keys() & zT.keys():
            # load data dictionary
            rc1, Mg, SB = zM.getday(yyyymmdd)
            rc2, TSI, TSI_UNC = zT.getday(yyyymmdd)
            rc = -1 if rc1 == -1 or rc2 == -1 else 0
            self._data[yyyymmdd] = (rc, TSI, Mg, SB, TSI_UNC)

            # keep track of min and max final times
            if rc == 0:
                self.nfinal += 1
                d = datetime.strptime(yyyymmdd, "%Y%m%d").date()
                if first:
                    self.date_min_final = d
                    self.date_max_final = d
                    first = False
                else:
                    if d < self.date_min_final:
                        self.date_min_final = d
                    if d > self.date_max_final:
                        self.date_max_final = d

    def getday(self, yyyymmdd):
        """
        get daily TSI, Mg and SB values for daily string yyyymmdd
        returns (rc, TSI, Mg, SB, TSI_UNC)
        rc = 0 (final), -1 (prelim), -2 (unavailable)
        """
        if yyyymmdd not in self._data:
            return -2, np.nan, np.nan, np.nan
        else:
            return self._data[yyyymmdd]

    def gettime(self, t):
        """
        get time-interpolated TSI, Mg, and SB values for datetime t
        returns (rc, TSI, Mg, SB, TSI_UNC)
        rc = 0 (final), -1 (prelim), -2 (unavailable)
        """

        # assume daily average valid at noon
        tnoon = datetime(t.year, t.month, t.day, 12)
        vnoon = self.getday(_yyyymmdd(tnoon))
        if t == tnoon:
            return vnoon

        # other noon bracketing t
        tother = tnoon + timedelta(days=(-1 if t < tnoon else +1))
        vother = self.getday(_yyyymmdd(tother))

        # fraction that the other daily average contributes
        fother = abs((t - tnoon).total_seconds()) / 86400.0

        # only interpolate if both days exist
        if vnoon[0] == -2 or vother[0] == -2:
            return (-2, np.nan, np.nan, np.nan)
        else:
            # now both days available
            TSI = vnoon[1] * (1 - fother) + vother[1] * fother
            Mg = vnoon[2] * (1 - fother) + vother[2] * fother
            SB = vnoon[3] * (1 - fother) + vother[3] * fother
            TSI_UNC = vnoon[4] * (1 - fother) + vother[4] * fother
            rc = -1 if vnoon[0] == -1 or vother[0] == -1 else 0
            return (rc, TSI, Mg, SB, TSI_UNC)

    def final_date_range(self, yyyymmdd0, yyyymmdd1, no_gaps=True):
        """final only lists over date range"""

        d0 = datetime.strptime(yyyymmdd0, "%Y%m%d").date()
        d1 = datetime.strptime(yyyymmdd1, "%Y%m%d").date()
        d = d0
        dates = []
        TSI = []
        Mg = []
        SB = []
        while d <= d1:
            v = self.getday(_yyyymmdd(d))
            if v[0] == 0:
                # only include finals
                dates.append(d)
                TSI.append(v[1])
                Mg.append(v[2])
                SB.append(v[3])
            elif no_gaps:
                raise RuntimeError("no_gaps: intervening non-final date found")
            d += timedelta(days=1)
        return dates, TSI, Mg, SB

    def _plot_all_final(self):
        """
        Recreates the original 3-panel plot (TSI, Mg, SB) but computes
        the x-axis from the class's final date range. It then shades the
        post-2024 region using _shade_v03_block.
        Draws on the CURRENT figure (3 rows, 1 column).
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from datetime import datetime, date

        # Pull the final date/series from the class API you already use in the writer
        start = _yyyymmdd(self.date_min_final)
        end = _yyyymmdd(self.date_max_final)
        dates, TSI, Mg, SB = self.final_date_range(start, end)

        # Normalize date_min_final -> date
        d0 = self.date_min_final
        if isinstance(d0, str):
            d0 = datetime.strptime(d0, "%Y-%m-%d").date()
        elif hasattr(d0, "date"):  # datetime -> date
            d0 = d0.date()

        # x-axis: years since first final date
        dy = np.array([(dd - d0).days / 365.25 for dd in dates], dtype=float)
        TSI = np.asarray(TSI, dtype=float)
        Mg = np.asarray(Mg, dtype=float)
        SB = np.asarray(SB, dtype=float)

        # Plot on 3 subplots, like the original
        plt.subplot(311)
        plt.plot(dy, TSI, "-")
        plt.ylabel("TSI")
        plt.subplot(312)
        plt.plot(dy, Mg, "-")
        plt.ylabel("Mg")
        plt.subplot(313)
        plt.plot(dy, SB, "-")
        plt.ylabel("SB")
        plt.xlabel(f"years since {d0.isoformat()}")

        # Shade post-2024 (v03-mapped) region
        self._shade_v03_block(dy)

    def _shade_v03_block(self, dy):
        """
        Decorate the current 3-panel figure by shading the post-2024 region and
        adding a small label. Uses ONLY the provided dy (years since start) and
        self.date_min_final to compute the boundary position.
        Call this AFTER you have already plotted TSI, Mg, SB on 3 subplots.
        """
        from datetime import datetime, date
        from matplotlib.patches import Patch
        from matplotlib.transforms import blended_transform_factory
        import numpy as np
        import matplotlib.pyplot as plt

        # Normalize date_min_final -> date
        d0 = self.date_min_final
        if isinstance(d0, str):
            d0 = datetime.strptime(d0, "%Y-%m-%d").date()
        elif hasattr(d0, "date"):  # datetime -> date
            d0 = d0.date()
        elif not hasattr(d0, "toordinal"):
            # Very defensive fallback (won't affect data, only label position)
            d0 = date(1882, 1, 1)

        boundary = date(2024, 1, 1)
        bdy_years = (boundary - d0).days / 365.25

        fig = plt.gcf()
        axes = fig.axes
        if not axes:
            return

        x_end = float(np.asarray(dy)[-1])
        for ax in axes:
            # shade post-2024 region
            ax.axvspan(
                bdy_years, x_end, facecolor="lightgoldenrodyellow", alpha=0.6, zorder=0
            )
            # vertical boundary line
            ax.axvline(bdy_years, color="k", linestyle="--", lw=1.1, zorder=3)
            # thin top ribbon + centered label above data
            trans = blended_transform_factory(ax.transData, ax.transAxes)
            ax.axvspan(
                bdy_years,
                x_end,
                ymin=0.965,
                ymax=0.985,
                facecolor="#ffd54f",
                alpha=0.9,
                zorder=4,
            )
            ax.text(
                (bdy_years + x_end) / 2.0,
                0.982,
                "v03r00 mapped to legacy units",
                transform=trans,
                ha="center",
                va="top",
                fontsize=8,
                zorder=5,
            )

        # legend in bottom panel
        legend_patch = Patch(
            facecolor="lightgoldenrodyellow",
            edgecolor="none",
            label="v03r00 mapped (post-2024)",
        )
        axes[-1].legend(handles=[legend_patch], loc="upper left", frameon=False)

    def output_final_textfile(self, filename):
        f = open(filename, "w")
        f.write("# NRLSSI2 daily input\n")
        f.write("# treat daily values as valid at 12:00 GMT\n")
        f.write("# yyyy doy TSI:W/m2 MgIndex   SBindex\n")
        for d, TSI, Mg, SB in zip(
            *self.final_date_range(
                _yyyymmdd(self.date_min_final), _yyyymmdd(self.date_max_final)
            )
        ):
            f.write(
                "  %04d %03d %8.3f %8.6f %9.4f\n"
                % (d.year, d.timetuple().tm_yday, TSI, Mg, SB)
            )
        f.close()


if __name__ == "__main__":
    DATADIR = "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/solar/NRLSSI2/data"
    OUTDIR = os.sep.join((os.environ["NOBACKUP"], "NRLSSI2", "output"))

    # mathomp4 personal paths (for reference):
    # DATADIR = '/discover/nobackup/mathomp4/NRLSSI2-Construct/GEOSradiation_GridComp/GEOS_RadiationShared/NRLSSI2/NRLSSI2/data'
    # OUTDIR = '/discover/nobackup/mathomp4/NRLSSI2-Construct/GEOSradiation_GridComp/GEOS_RadiationShared/NRLSSI2/NRLSSI2/output'

    z = TSI_Mg_SB_Merged_Daily(DATADIR)
    # print('16000101', z.getday('16000101'))
    # print('20161130', z.getday('20161130'))
    # print('20161201', z.getday('20161201'))
    # print('20161202', z.getday('20161202'))
    # t = datetime.strptime('2016-12-01 10:00:00','%Y-%m-%d %H:%M:%S')
    # print(t, z.gettime(t))

    z.output_final_textfile(os.sep.join((OUTDIR, "NRLSSI2.vYYYY.txt")))

    plt.figure(figsize=(6, 12))
    z._plot_all_final()
    plt.savefig(
        os.sep.join(("gx", "TSInMgSB_plot_all_final.png")),
        pad_inches=0.33,
        bbox_inches="tight",
        dpi=100,
    )
    plt.show()
