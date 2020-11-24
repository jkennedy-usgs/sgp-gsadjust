"""
plots/datum.py
==============

Plots for Datum objects.
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
import datetime as dt

import matplotlib
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.dates import date2num


class PlotDatumComparisonTimeSeries(QtWidgets.QDialog):
    """
    Plots time series of observed vs. adjusted values at datum stations.
    """
    def __init__(self, obsTreeModel, parent=None):
        super(PlotDatumComparisonTimeSeries, self).__init__(parent)
        self.setWindowTitle("GSadjust results")
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        xdata_all, ydata_obs_all, ydata_adj_all, unique_datum_names = self.get_data(
            obsTreeModel
        )
        self.plot(xdata_all, ydata_obs_all, ydata_adj_all, unique_datum_names)

    def get_data(self, obsTreeModel):
        datum_names = []
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            for datum in survey.datums:
                datum_names.append(datum.station)

        unique_datum_names = list(set(datum_names))

        xdata_all, ydata_obs_all, ydata_adj_all = [], [], []
        for name in unique_datum_names:
            xdata, ydata_obs, ydata_adj = [], [], []
            for i in range(obsTreeModel.rowCount()):
                survey = obsTreeModel.invisibleRootItem().child(i)
                for datum in survey.datums:
                    if (
                        datum.station == name
                        and datum.residual > -998
                        and datum.checked == 2
                    ):
                        xdata.append(
                            date2num(dt.datetime.strptime(survey.name, "%Y-%m-%d"))
                        )
                        ydata_obs.append(datum.g)
                        ydata_adj.append(datum.g + datum.residual)
            ydata_obs = [i - ydata_obs[0] for i in ydata_obs]
            ydata_adj = [i - ydata_adj[0] for i in ydata_adj]
            xdata_all.append(xdata)
            ydata_obs_all.append(ydata_obs)
            ydata_adj_all.append(ydata_adj)
        return xdata_all, ydata_obs_all, ydata_adj_all, unique_datum_names

    def plot(self, xdata_all, ydata_obs_all, ydata_adj_all, names):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        i = 0
        for xdata, ydata_obs, ydata_adj in zip(xdata_all, ydata_obs_all, ydata_adj_all):
            a = ax.plot(xdata, ydata_obs, "-o", label=names[i] + "_obs")
            line_color = a[0].get_color()
            b = ax.plot(xdata, ydata_adj, "--o", c=line_color, label=names[i] + "_adj")
            i += 1
        ax.legend()
        ax.set_ylabel("Relative gravity, in \u00b5Gal")
        ax.xaxis_date()
        self.canvas.draw()


class PlotDatumCompare(QtWidgets.QDialog):
    """
    Bar plot of difference between specified datum (in datum table) and adjustment
    result.
    """

    def __init__(self, survey, parent=None):
        super(PlotDatumCompare, self).__init__(parent)
        self.setWindowTitle("GSadjust results")
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.survey = survey
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        diff, lbl = self.get_data()
        self.plot(diff, lbl)

    def get_data(self):
        survey = self.survey
        diff, lbl = [], []
        for datum in survey.datums:
            for result in survey.results:
                if result.station == datum.station:
                    diff.append(result.g - datum.g + datum.meas_height * datum.gradient)
                    lbl.append(datum.station)

        return diff, lbl

    def plot(self, diff, lbl):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ind = range(len(diff))
        ax.bar(ind, diff)
        ax.set_title("Adjusted g minus measured g (microGal)")
        ax.set_xticks(ind)
        ax.set_xticklabels(lbl)
        self.canvas.draw()
