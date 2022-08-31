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
import numpy.linalg
from numpy import ceil
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.dates as mdates
import matplotlib.ticker as tkr
import copy

from ..gui.messages import MessageBox

consistent_date_axes = True
presentation_style = False

if consistent_date_axes:
    x_min = datetime.datetime(2016, 1, 1)
    x_max = datetime.datetime(2022, 10, 1)


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
                 t_offset=None, optimize=False,
                 t_threshold=0, rising='all',
                 parent=None):
        super(PlotNwis, self).__init__(parent)
        self.nwis_data = nwis_data
        self.grav_data = g_data
        self.nwis_station = nwis_station
        self.g_station = g_station
        self.t_offset = t_offset
        self.t_threshold = datetime.timedelta(days=t_threshold)
        self.optimize = optimize
        self.rising = rising
        self.setWindowTitle("Gravity - Groundwater level comparison")
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
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
        if plot_x is None:
            MessageBox.critical("Error",
                                f"Error retrieving plot data")
        ax = self.figure.add_subplot(212)
        self.plot_sy(ax, plot_x, plot_y)
        self.figure.subplots_adjust(
            bottom=0.15, top=0.85, hspace=0.4, left=0.20, right=0.85
        )
        self.canvas.draw()

    def find_best_offset(self):
        offsets = range(120)
        offsets = [o - 60 for o in offsets]
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
            except ValueError as e:
                return None, None, None
            r2 = np.corrcoef(plot_x, plot_y)[0, 1] ** 2
            if r2 > best_r2:
                best_r2 = r2
                best_offset = offset
                best_x, best_y = plot_x, plot_y

        return best_x, best_y, best_offset * -1


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
            if self.meters:  # NWIS default is feet
                nwis_y = [meas * 0.3048 for meas in nwis_y]
            ax2.plot(nwis_x, nwis_y)
        if self.nwis_data["discrete_x"]:
            nwis_x_discrete = self.nwis_data["discrete_x"]
            nwis_y_discrete = self.nwis_data["discrete_y"]
            if self.meters:  # NWIS default is feet
                nwis_y_discrete = [meas * 0.3048 for meas in nwis_y_discrete]
            ax2.plot(nwis_x_discrete,
                     nwis_y_discrete,
                     "s",
                     color=(0.12156, 0.46667, 0.705882))

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
    def interpolate_gap(data_x, data_y, idx_neg, idx_pos, g_date):
        x1 = data_x[idx_neg]
        x2 = data_x[idx_pos]
        y1 = data_y[idx_neg]
        y2 = data_y[idx_pos]

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
        noninterpolate_threshold = datetime.timedelta(days=2)
        # Otherwise, interpolate a value if there is data within interpolate_threshold
        # interpolate_threshold = datetime.timedelta(days=50)

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
                cont_gap < disc_gap and cont_gap < self.t_threshold
            ):  # interpolate the data type with the smaller gap
                interp_dtw = self.interpolate_gap(
                    nwis_data['continuous_x'],
                    nwis_data['continuous_y'],
                    idx_closest_neg_cont,
                    idx_closest_pos_cont,
                    g_date,
                )
            elif disc_gap < cont_gap and disc_gap < self.t_threshold:
                interp_dtw = self.interpolate_gap(
                    nwis_data['discrete_x'],
                    nwis_data['discrete_y'],
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

        fx, fy = self.parse_rise_or_fall(plot_x, plot_y)
        return fx, fy

    def parse_rise_or_fall(self, plot_x, plot_y):
        fx, fy = [], []
        if self.rising == 'all':
            return plot_x, plot_y

        for idx, val in enumerate(plot_x[1:]):
            if self.rising == 'rise':
                if val - plot_x[0] > 0:
                    fx.append(plot_x[idx])
                    fx.append(val)
                    fy.append(plot_y[idx])
                    fy.append(plot_y[idx+1])
            elif self.rising == 'fall':
                if val - plot_x[0] < 0:
                    fx.append(plot_x[idx])
                    fx.append(val)
                    fy.append(plot_y[idx])
                    fy.append(plot_y[idx+1])

        # This should maintain order, vs. coverting to a set
        fx, fy = self.remove_dups(fx, fy)
        return fx, fy

    def remove_dups(self, lst_x, lst_y):
        new_x, new_y, prev_item = [], [], None
        new_x += lst_x[:2]
        new_y += lst_y[:2]
        eps = 0.000001
        for idx, item in enumerate(lst_x[2:], start=2):
            if abs(lst_x[idx] - new_x[-1]) > eps and abs(lst_x[idx] - new_y[-1]) > eps:
                new_x.append(lst_x[idx])
                new_y.append(lst_y[idx])
        return new_x, new_y

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
                "SY = %0.2f ± %0.02f" % (poly[0], np.sqrt(cov[0, 0])),
                transform=ax.transAxes,
            )
            ax.text(0.75, 0.10, "r² = %0.2f" % r2, transform=ax.transAxes)
            if self.t_offset:
                ax.set_title("Specific Yield, GWL time offset: {} days".format(
                    abs(self.t_offset)
                ))
            else:
                ax.set_title("Specific Yield")

        except (ValueError, TypeError, numpy.linalg.LinAlgError) as e:
            print(e)
