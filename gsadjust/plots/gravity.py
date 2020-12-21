"""
plots/gravity.py
===============

Gravity change time-series plot.
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
import numpy as np
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class PlotGravityChange(QtWidgets.QDialog):
    def __init__(self, dates, table, parent=None):
        super(PlotGravityChange, self).__init__(parent)
        self.setWindowTitle("Gravity Change Time Series")
        self.setWhatsThis(
            "Click on a line in the legend to toggle visibility. Right-click anywhere"
            " to hide all lines. Middle-click to show all lines."
        )
        self.figure = matplotlib.figure.Figure(figsize=(10, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        ax = self.plot(dates, table)

        self.leg = self.interactive_legend(ax)
        self.canvas.draw()

    def plot(self, dates, table):
        ncols = len(dates)
        nstations = len(table[0])
        stations = table[0]
        ncol = int(np.ceil(nstations / 24))
        ax = self.figure.add_subplot(111)
        right_margin = 1 - ncol / 8
        self.figure.subplots_adjust(right=right_margin)
        for i in range(nstations):
            xdata, ydata = [], []
            for idx, col in enumerate(table[1:ncols]):
                if not col[i] == "-999":
                    if not ydata:
                        ydata.append(0)
                        ydata.append(float(col[i]))
                        xdata.append(dates[idx])
                        xdata.append(dates[idx + 1])
                    else:
                        ydata.append(float(col[i]) + ydata[-1])
                        xdata.append(dates[idx + 1])
            cmap = matplotlib.cm.get_cmap("gist_ncar")
            ax.plot(xdata, ydata, "-o", color=cmap(i / nstations), label=stations[i])
        ax.set_ylabel("Gravity change, in µGal")
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=ncol)
        self.figure.autofmt_xdate()
        return ax

    def interactive_legend(self, ax=None):
        if ax is None:
            return
        if ax.legend_ is None:
            ax.legend()
        return InteractiveLegend(ax.get_legend())


class InteractiveLegend(object):
    """
    Allows clicking on the legend to show/hide lines
    """
    def __init__(self, legend):
        self.legend = legend
        self.fig = legend.axes.figure

        self.lookup_artist, self.lookup_handle = self._build_lookups(legend)
        self._setup_connections()

        self.update()

    def _setup_connections(self):
        for artist in self.legend.texts + self.legend.legendHandles:
            artist.set_picker(10)  # 10 points tolerance

        self.fig.canvas.mpl_connect("pick_event", self.on_pick)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)

    def _build_lookups(self, legend):
        labels = [t.get_text() for t in legend.texts]
        handles = legend.legendHandles
        label2handle = dict(zip(labels, handles))
        handle2text = dict(zip(handles, legend.texts))

        lookup_artist = {}
        lookup_handle = {}
        for artist in legend.axes.get_children():
            if artist.get_label() in labels:
                handle = label2handle[artist.get_label()]
                lookup_handle[artist] = handle
                lookup_artist[handle] = artist
                lookup_artist[handle2text[handle]] = artist

        lookup_handle.update(zip(handles, handles))
        lookup_handle.update(zip(legend.texts, handles))

        return lookup_artist, lookup_handle

    def on_pick(self, event):
        handle = event.artist
        if handle in self.lookup_artist:
            artist = self.lookup_artist[handle]
            artist.set_visible(not artist.get_visible())
            self.update()

    def on_click(self, event):
        if event.button == 3:
            visible = False
        elif event.button == 2:
            visible = True
        else:
            return

        for artist in self.lookup_artist.values():
            artist.set_visible(visible)
        self.update()

    def update(self):
        for artist in self.lookup_artist.values():
            handle = self.lookup_handle[artist]
            if artist.get_visible():
                handle.set_visible(True)
            else:
                handle.set_visible(False)
        self.fig.canvas.draw()
