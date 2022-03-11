"""
plots/nwis.py
=================

Residual histogram plot
--------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
import matplotlib
import datetime
import numpy as np
from numpy import ceil
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.dates as mdates
import matplotlib.ticker as tkr
import copy

consistent_date_axes = True
presentation_style = False

if consistent_date_axes:
    x_min = datetime.datetime(2016, 1, 1)
    x_max = datetime.datetime(2021, 10, 1)


def func(x, pos):
    s = "{}".format(x)
    return s


myFmt = mdates.DateFormatter("%Y")
y_format = tkr.FuncFormatter(func)


class PlotNwis(QtWidgets.QDialog):
    """
    Two-panel plot showing gravity and gwl time series (top) and
    specific yield plot (storage vs. gwl change, bottom)
    """

    def __init__(self, nwis_station, g_station, nwis_data, g_data,
                 t_offset=None, optimize=False, parent=None):
        super(PlotNwis, self).__init__(parent)
        self.nwis_data = nwis_data
        self.grav_data = g_data
        self.nwis_station = nwis_station
        self.g_station = g_station
        self.t_offset = t_offset
        self.optimize = optimize
        self.setWindowTitle("Gravity - Groundwater level comparison")
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

        self.convert_to_water = False
        self.meters = True
        self.a10sd = 5.0
        self.plot()
        self.resize(700, 800)

    def plot(self):
        self.figure.clear()
        ax = self.figure.add_subplot(211)
        grav_x = self.grav_data[0]
        grav_y = self.grav_data[1]
        y0 = grav_y[0]

        if self.convert_to_water:
            if self.meters:
                ytemp = [(p - y0) / 41.9 for p in grav_y]
                a10sd = self.a10_sd / 41.9
            else:
                ytemp = [(p - y0) / 12.77 for p in grav_y]
                a10sd = self.a10_sd / 12.77
            rng = max(ytemp) - min(ytemp)
            self.half_rng = ceil(rng)
        else:
            ytemp = [(p - y0) for p in grav_y]
            a10sd = self.a10sd
        grav_y = ytemp
        self.plot_wl_and_g(ax, grav_x, grav_y, a10sd)
        if self.optimize:
            plot_x, plot_y, self.t_offset = self.find_best_offset()
        else:
            plot_x, plot_y = self.get_sy_data(self.nwis_data, self.grav_data)
        ax = self.figure.add_subplot(212)
        self.plot_sy(ax, plot_x, plot_y)
        self.figure.subplots_adjust(
            bottom=0.15, top=0.85, hspace=0.4, left=0.20, right=0.85
        )
        self.canvas.draw()

    def find_best_offset(self):
        offsets = range(60)
        best_r2, best_offset = 0, 0
        best_x, best_y = self.get_sy_data(self.nwis_data, self.grav_data)
        for offset in offsets:
            # offset *= -1
            trial_nwis = copy.copy(self.nwis_data)
            trial_nwis['continuous_x'] = [a + datetime.timedelta(offset) for a in
                                          self.nwis_data['continuous_x']]
            trial_nwis['discrete_x'] = [a + datetime.timedelta(offset) for a in
                                          self.nwis_data['discrete_x']]
            plot_x, plot_y = self.get_sy_data(trial_nwis, self.grav_data)
            try:
                poly, cov = np.polyfit(plot_x, plot_y, 1, cov=True)
            except TypeError as e:
                continue
            r2 = np.corrcoef(plot_x, plot_y)[0, 1] ** 2
            if r2 > best_r2:
                best_r2 = r2
                best_offset = offset
                best_x, best_y = plot_x, plot_y

        return (best_x, best_y, best_offset * -1)


    def plot_wl_and_g(self, ax, grav_x, grav_y, a10sd):
        ax.errorbar(grav_x, grav_y, yerr=a10sd, fmt="kd")
        ax.set_title(f"NWIS: {self.nwis_station}  g: {self.g_station}")
        # ax = plt.gca()
        if self.convert_to_water:
            ax.set_ylim(0 - self.half_rng, 0 + self.half_rng)
        ax2 = ax.twinx()

        if self.nwis_data["continuous_x"]:
            nwis_x = self.nwis_data["continuous_x"]
            nwis_y = self.nwis_data["continuous_y"]
        else:
            nwis_x = self.nwis_data["discrete_x"]
            nwis_y = self.nwis_data["discrete_y"]

        if self.meters:  # NWIS default is feet
            nwis_y = [meas * 0.3048 for meas in nwis_y]
        try:
            ax2.plot(nwis_x, nwis_y)
        except:
            return
        ax2.invert_yaxis()

        # Remove scientific notation from axes labels
        # ax.yaxis.get_major_formatter().set_useOffset(False)

        # Add commas to y-axis tick mark labels
        ax.yaxis.set_major_formatter(y_format)
        # Set x-axis tick labels to just show year
        ax.xaxis.set_major_formatter(myFmt)

        if not consistent_date_axes:
            # Adjust ticks so they fall on Jan 1 and extend past the range of the data. If there
            # are data in January and December, add another year so that there is plenty of space.
            start_month = self.grav_data[0][0].month
            start_year = self.grav_data[0][0].year
            end_month = self.grav_data[0][-1].month
            end_year = self.grav_data[0][-1].year
            if start_month == 1:
                start_year = start_year - 1
            if end_month == 12:
                end_year = end_year + 1
            xticks = []
            for iii in range(start_year, end_year + 2):
                xticks.append(datetime.datetime(iii, 1, 1))
            ax.set_xticks(xticks)
        else:
            ax.set_xlim(x_min, x_max)
        if self.convert_to_water:
            if self.meters:
                ax.set_ylabel("Storage change, in m of water")
            else:
                ax.set_ylabel("Storage change, in ft of water")
        else:
            ax.set_ylabel("Gravity change, in microGal")
        if self.meters:
            ax2.set_ylabel("Depth to groundwater, in meters")
        else:
            ax2.set_ylabel("Depth to groundwater, in feet")

    @staticmethod
    def find_closest(g_date, data):
        repdate = np.repeat(g_date, len(data))
        delta = np.asarray(data) - repdate
        min_delta = min(np.absolute(delta))
        idx = np.argmin(np.absolute(delta))
        return delta, min_delta, idx

    @staticmethod
    def find_gap(delta):
        if any(i < datetime.timedelta(days=0) for i in delta) and any(
            i > datetime.timedelta(days=0) for i in delta
        ):
            closest_neg = max(
                [i for i in delta if i <= datetime.timedelta(days=0)]
            )  # time delta to closest negative diff
            closest_pos = min([i for i in delta if i >= datetime.timedelta(days=0)])
            (idx_closest_neg,) = np.nonzero(delta == closest_neg)[0]
            (idx_closest_pos,) = np.nonzero(delta == closest_pos)[0]
            gap = np.absolute(closest_neg) + closest_pos
            return idx_closest_neg, idx_closest_pos, gap
        else:
            return -999, -999, datetime.timedelta(days=1000000)

    @staticmethod
    def interpolate_gap(data, idx_neg, idx_pos, g_date):
        x1 = data["continuous_x"][idx_neg]
        x2 = data["continuous_x"][idx_pos]
        y1 = data["continuous_y"][idx_neg]
        y2 = data["continuous_y"][idx_pos]

        x1 = x1.toordinal()
        x2 = x2.toordinal()
        poly = np.polyfit([x1, x2], [y1, y2], 1)
        p = np.poly1d(poly)
        interpolated_dtw = p(g_date.toordinal())
        return interpolated_dtw

    def get_sy_data(self, nwis_data, grav_data):
        plot_x, plot_y = [], []
        min_delta_cont, min_delta_disc = (
            datetime.timedelta(days=1000000),
            datetime.timedelta(days=1000000),
        )
        # If within a threshold, just take the nearest data point
        noninterpolate_threshold = datetime.timedelta(days=5)
        # Otherwise, interpolate a value if there is data within interpolate_threshold
        interpolate_threshold = datetime.timedelta(days=50)

        if (
            nwis_data["continuous_x"] is None
                and nwis_data["discrete_x"] is None
        ):
            return

        # Iterate through the gravity values for a given station
        for g_idx, g_date in enumerate(grav_data[0]):
            # find closest continuous data
            if nwis_data["continuous_x"]:
                delta_cont, min_delta_cont, closest_cont_idx = self.find_closest(
                    g_date, nwis_data["continuous_x"]
                )
            if nwis_data["discrete_x"]:
                delta_disc, min_delta_disc, closest_disc_idx = self.find_closest(
                    g_date, nwis_data["discrete_x"]
                )

            # check threshold
            if (
                min_delta_cont < noninterpolate_threshold and min_delta_cont < min_delta_disc
            ):
                plot_x.append(nwis_data["continuous_y"][closest_cont_idx])
                plot_y.append(grav_data[1][g_idx])
                continue
            elif (
                min_delta_disc < noninterpolate_threshold and min_delta_disc < min_delta_cont
            ):
                plot_x.append(nwis_data["discrete_y"][closest_disc_idx])
                plot_y.append(grav_data[1][g_idx])
                continue

            # No water-level measurements are very close. Check if we can interpolate.
            cont_gap, disc_gap = (
                datetime.timedelta(days=1000000),
                datetime.timedelta(days=1000000),
            )
            interp_dtw = None
            if nwis_data["continuous_x"]:  # calculate continuous gap
                (
                    idx_closest_neg_cont,
                    idx_closest_pos_cont,
                    cont_gap,
                ) = self.find_gap(delta_cont)
            if nwis_data["discrete_x"]:  # calculate discrete gap
                (
                    idx_closest_neg_disc,
                    idx_closest_pos_disc,
                    disc_gap,
                ) = self.find_gap(delta_disc)
            if (
                cont_gap < disc_gap and cont_gap < interpolate_threshold
            ):  # interpolate the data type with the smaller gap
                interp_dtw = self.interpolate_gap(
                    nwis_data,
                    idx_closest_neg_cont,
                    idx_closest_pos_cont,
                    g_date,
                )
            elif disc_gap < cont_gap and disc_gap < interpolate_threshold:
                interp_dtw = self.interpolate_gap(
                    nwis_data,
                    idx_closest_neg_disc,
                    idx_closest_pos_disc,
                    g_date,
                )
            if interp_dtw:
                plot_x.append(interp_dtw)
                plot_y.append(grav_data[1][g_idx])

        if len(plot_y) > 1:
            plot_y = [(y - plot_y[0]) / 41.9 for y in plot_y]
            plot_x = [(x - plot_x[0]) * -0.3048 for x in plot_x] # Negative to convert depth-to-water to elevation
        return (plot_x, plot_y)

    def plot_sy(self, ax, plot_x, plot_y):
        if presentation_style:
            font = {"family": "normal", "weight": "bold", "size": 16}

            # plt.rc('font', **font)
            # ax.subplots_adjust(bottom=0.15, top=0.85, hspace=0.4,
            #                     left=0.25, right=0.85)

        try:  # Sometimes polyfit fails, even if there's 3 points?



            poly, cov = np.polyfit(plot_x, plot_y, 1, cov=True)
            r2 = np.corrcoef(plot_x, plot_y)[0, 1] ** 2
            line_x = np.linspace(min(plot_x) - 0.2, max(plot_x) + 0.2, 10)


            p = np.poly1d(poly)
            line_y = p(line_x)
            ax.plot(plot_x, plot_y, ".")
            ax.plot(line_x, line_y)

            ax.set_ylabel(
                "Change in water storage,\nmeters of free-standing water\n(from gravity data)"
            )
            ax.set_xlabel("Change in groundwater level (meters)")
            ax.text(
                0.70,
                0.17,
                "Sy = %0.2f ± %0.02f" % (poly[0], np.sqrt(cov[0, 0])),
                transform=ax.transAxes,
            )
            ax.text(0.75, 0.12, "r^2 = %0.2f" % r2, transform=ax.transAxes)
            if self.t_offset:
                ax.set_title("Specific Yield, GWL time offset: {} days".format(
                    abs(self.t_offset)
                ))
            else:
                ax.set_title("Specific Yield")

        except (ValueError, TypeError) as e:
            print(e)
