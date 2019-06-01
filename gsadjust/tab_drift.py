#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tab_drift.py
===============

PyQt graphical elements on the drift tab of GSadjust.
--------------------------------------------------------------------------------------------------------------------


This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
import datetime as dt
import logging

import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.dates import DateFormatter
from matplotlib.dates import date2num
from matplotlib.figure import Figure
from scipy.interpolate import UnivariateSpline

from data_objects import Delta
from gui_objects import IncrMinuteTimeEdit, show_message
from pyqt_models import DeltaTableModel, RomanTableModel, ObsTreeLoop


###########################################################################
# GSadjust drift tab
###########################################################################
class TabDrift(QtWidgets.QWidget):
    station_label = None

    def __init__(self, parent):
        super(TabDrift, self).__init__()
        self.parent = parent
        self.dpi = 100
        self.popup_menu = QtWidgets.QMenu(None)
        # Main window
        layout_main = QtWidgets.QVBoxLayout()
        main_vsplitter_window = QtWidgets.QSplitter(QtCore.Qt.Vertical, self)

        # Setup drift figures (Roman and Continuous). Only one will be shown at a time.
        # Drift figure: default and roman
        self.drift_window = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self)
        self.drift_fig = Figure((3.0, 5.0), dpi=self.dpi, facecolor='white')
        self.drift_single_canvas = FigureCanvas(self.drift_fig)
        self.drift_fig.subplots_adjust(wspace=0.3)
        self.axes_drift_single = self.drift_fig.add_subplot(111)

        # Drift figure: continuous
        # Plot pannel
        self.drift_cont_plotpanel = QtWidgets.QSplitter(QtCore.Qt.Vertical, self)
        # Top plot - lines
        self.drift_cont_figtop = Figure((3.0, 5.0), dpi=self.dpi, facecolor='white')
        self.drift_cont_canvastop = FigureCanvas(self.drift_cont_figtop)
        self.drift_cont_figtop.subplots_adjust(wspace=0.3)
        self.axes_drift_cont_upper = self.drift_cont_figtop.add_subplot(111)
        # Bottom plot - drift curves
        self.drift_cont_figbot = Figure((3.0, 5.0), dpi=self.dpi, facecolor='white')
        self.drift_cont_canvasbot = FigureCanvas(self.drift_cont_figbot)
        self.drift_cont_figbot.subplots_adjust(wspace=0.3)
        self.axes_drift_cont_lower = self.drift_cont_figbot.add_subplot(111)

        # Drift tab tables
        self.dg_samples_table = DeltaTableModel()
        self.delta_view = QtWidgets.QTableView()
        # Hide std_for_adj and residual columns
        self.cont_label_widget = QtWidgets.QWidget()
        self.dg_avg_model = RomanTableModel()
        self.dg_samples_view = QtWidgets.QTableView()

        #######################################################################
        # Widgets for right-hand display of drift controls/options
        #######################################################################
        # Drift method widget
        self.driftmethod_comboboxbox_key = {0: "None",
                                            1: "Network adjustment",
                                            2: "Roman (interpolate)",
                                            3: "Continuous model"}
        self.driftmethod_comboboxbox = QtWidgets.QComboBox()
        self.driftmethod_comboboxbox.activated.connect(self.set_drift_method)
        for item in self.driftmethod_comboboxbox_key.values():
            self.driftmethod_comboboxbox.addItem(item)

        # Widget to remove dg-observations with a long elapsed time in between
        self.drift_screen_elapsed_time = QtWidgets.QCheckBox('Max. time between repeats (hh:mm)')
        self.drift_screen_elapsed_time.setChecked(False)
        self.drift_screen_elapsed_time.stateChanged.connect(self.time_extent_changed)
        self.drift_time_spinner = IncrMinuteTimeEdit(QtCore.QTime(1, 0))
        self.drift_time_spinner.timeChanged.connect(self.time_extent_changed)
        self.drift_time_spinner.setDisplayFormat("hh:mm")

        # Widget to add horizontal-extent lines to drift-rate plot
        self.drift_plot_hz_extent = QtWidgets.QCheckBox('Show time-extent of drift observation')
        self.drift_plot_hz_extent.setChecked(False)
        self.drift_plot_hz_extent.stateChanged.connect(self.plot_drift)

        self.tension_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.tension_slider.setRange(10, 2500)
        self.tension_slider.setValue(1250)
        self.tension_slider.setEnabled(False)
        self.tension_slider.valueChanged.connect(self.update_tension)
        self.tension_label = QtWidgets.QLabel()
        self.tension_label.setText('{:2.2f}'.format(self.tension_slider.value()))
        self.tension_label.setAlignment(QtCore.Qt.AlignCenter)
        self.drift_cont_methodcombobox_key = {0: "None",
                                              1: "Spline",
                                              2: "1st order polynomial",
                                              3: "2nd order polynomial",
                                              4: "3rd order polynomial"}
        self.drift_polydegree_combobox = QtWidgets.QComboBox()
        self.drift_polydegree_combobox.activated.connect(self.drift_combobox_updated)
        for item in self.drift_cont_methodcombobox_key.values():
            self.drift_polydegree_combobox.addItem(item)
        self.drift_cont_behaviorcombobox_key = {0: "Extrapolate",
                                                1: "Constant"}
        self.drift_cont_startendcombobox = QtWidgets.QComboBox()
        self.drift_cont_startendcombobox.activated.connect(self.drift_combobox_updated)
        for item in self.drift_cont_behaviorcombobox_key.values():
            self.drift_cont_startendcombobox.addItem(item)

        self.offset_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.offset_slider.setRange(0, 10)
        self.offset_slider.setValue(0)
        self.offset_slider.valueChanged.connect(self.plot_drift)
        drift_controls = QtWidgets.QWidget()
        drift_cont_control_layout = QtWidgets.QHBoxLayout()
        drift_control_sublayout = QtWidgets.QVBoxLayout()
        grid_widget = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Drift correction method'), 0, 0)
        grid.addWidget(self.driftmethod_comboboxbox, 0, 1)
        grid.addWidget(QtWidgets.QLabel('Drift model type'), 1, 0)
        grid.addWidget(self.drift_polydegree_combobox, 1, 1)
        grid.addWidget(QtWidgets.QLabel('Behavior at start/end:'), 2, 0)
        grid.addWidget(self.drift_cont_startendcombobox, 2, 1)
        grid.addWidget(self.drift_screen_elapsed_time, 3, 0)
        grid.addWidget(self.drift_plot_hz_extent, 4, 0)
        grid.addWidget(self.drift_time_spinner, 3, 1)
        grid.addWidget(QtWidgets.QLabel('Spline tension:'), 6, 0)
        grid.addWidget(self.tension_slider, 6, 1)
        grid.addWidget(self.tension_label, 6, 2)
        grid.addWidget(QtWidgets.QLabel('Vertical line offset'), 5, 0)
        grid.addWidget(self.offset_slider, 5, 1)
        grid_widget.setLayout(grid)
        drift_control_sublayout.addWidget(grid_widget)

        self.tare_view = QtWidgets.QTableView()
        self.tare_view.clicked.connect(self.update_tares)
        self.tare_view.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tare_view.customContextMenuRequested.connect(self.tare_context_menu)
        self.resultsProxyModel = QtCore.QSortFilterProxyModel(self)

        self.tare_popup_menu = QtWidgets.QMenu("tare Popup Menu", self)
        self.mnDeleteTare = QtWidgets.QAction('Delete tare', self)
        self.mnDeleteTare.triggered.connect(self.parent.delete_tare)

        lbl = QtWidgets.QLabel("Tares")
        lbl.setFont(QtGui.QFont("Times", 11, QtGui.QFont.Bold))
        lbl.setFixedHeight(30)
        drift_control_sublayout.addItem(QtWidgets.QSpacerItem(40, 42))
        drift_control_sublayout.addWidget(lbl)
        drift_control_sublayout.addWidget(self.tare_view)

        control_subwidget = QtWidgets.QWidget()
        control_subwidget.setLayout(drift_control_sublayout)
        drift_cont_control_layout.addWidget(control_subwidget)
        drift_cont_control_layout.addStretch()
        drift_controls.setLayout(drift_cont_control_layout)

        self.drift_cont_plotpanel.addWidget(self.drift_cont_canvastop)
        self.drift_cont_plotpanel.addWidget(self.drift_cont_canvasbot)
        self.drift_window.addWidget(self.drift_cont_plotpanel)

        self.drift_window.addWidget(self.drift_single_canvas)
        self.drift_window.addWidget(self.drift_cont_plotpanel)
        self.drift_window.addWidget(drift_controls)
        main_vsplitter_window.addWidget(self.drift_window)
        self.drift_single_canvas.hide()

        lbls = QtWidgets.QHBoxLayout()
        lbl1 = QtWidgets.QLabel("Relative-gravity differences (delta-g's)", self)
        lbls.addWidget(lbl1)
        self.cont_label_widget.setLayout(lbls)
        self.cont_label_widget.setFixedHeight(30)
        main_vsplitter_window.addWidget(self.cont_label_widget)
        self.cont_label_widget.hide()

        self.roman_label_widget = QtWidgets.QWidget()
        lbls = QtWidgets.QHBoxLayout()
        lbl1 = QtWidgets.QLabel("Relative-gravity differences (delta-g's)", self)
        lbl2 = QtWidgets.QLabel("Average gravity differences", self)
        lbls.addWidget(lbl1)
        lbls.addWidget(lbl2)
        self.roman_label_widget.setLayout(lbls)
        self.roman_label_widget.setFixedHeight(30)
        main_vsplitter_window.addWidget(self.roman_label_widget)
        self.roman_label_widget.hide()

        # dg table (Roman method)
        self.dg_samples_view.setModel(self.dg_avg_model)
        self.delta_view.setModel(self.dg_samples_table)
        main_hsplitter_window = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self)
        main_hsplitter_window.addWidget(self.dg_samples_view)
        main_hsplitter_window.addWidget(self.delta_view)
        main_hsplitter_window.setMinimumHeight(300)
        main_vsplitter_window.addWidget(main_hsplitter_window)
        self.delta_view.show()
        self.dg_samples_view.hide()

        layout_main.addWidget(main_vsplitter_window)
        self.setLayout(layout_main)

    def time_extent_changed(self):
        self.parent.deltas_update_required()
        self.plot_drift()

    def drift_newpoint_picked(self, event):
        if event.button == 3:
            self.drift_rate_context_menu()

    def drift_point_picked(self, event):
        if event.mouseevent.button == 3:
            self.drift_rate_context_menu(from_pick=True)

    def drift_rate_context_menu(self, from_pick=False):
        """
        Not functional (other than showing the menu). Should allow points to be excluded, or artificial points added,
        to the continuous drift correction.
        :param from_pick: Boolean, True if a point was picked
        """
        if from_pick:
            add = QtWidgets.QAction(QtGui.QIcon(""), "Add point to drift model", self,
                                    triggered=self.drift_cont_addpoint,
                                    enabled=False)
            remove = QtWidgets.QAction(QtGui.QIcon(""), "Remove point from model", self,
                                       triggered=self.drift_cont_removepoint)
            self.popup_menu.addAction(remove)
        else:
            add = QtWidgets.QAction(QtGui.QIcon(""), "Add point to drift model", self,
                                    triggered=self.drift_cont_addpoint)
            remove = QtWidgets.QAction(QtGui.QIcon(""), "Remove point from model", self,
                                       triggered=self.drift_cont_removepoint,
                                       enabled=False)

        self.popup_menu.addAction(add)
        self.popup_menu.addAction(remove)
        cursor = QtGui.QCursor()
        self.popup_menu.popup(cursor.pos())

    def drift_cont_removepoint(self):
        pass

    def drift_cont_addpoint(self):
        pass

    def show_line_label(self, event, axes):
        """
        Shows the station name in the upper left of the drift plot when a line is clicked.
        :param event: Matplotlib event
        :param axes: Current axes (differs for none|netadj|roman vs continuous)
        """
        thisline = event.artist
        if self.station_label is not None:
            self.station_label.set_text('')
        self.station_label = axes.text(0.05, 0.95, thisline.name, horizontalalignment='center',
                                       verticalalignment='center',
                                       transform=axes.transAxes)
        axes.figure.canvas.draw()

    @staticmethod
    def screen_for_elapsed_time(plot_data, elapsed_time):
        """
        We may want to exclude repeat observations with a lot of elapsed time between occupations.
        :param plot_data: input data (list of lists)
        :param elapsed_time: maximum time to be considered a repeat, in minutes
        :return: list with same format as plot_data
        """
        new_data = []
        for line in plot_data:
            x = [x for x in line[0]]
            y = [y for y in line[1]]
            new_x, new_y = [], []
            # first = True
            # start_idx = 0
            i = 0
            for i in range(1, len(x)):
                # if first:
                #     first = False
                #     pass
                # else:
                x_diff = x[i] - x[i - 1]
                if x_diff * 1440 < elapsed_time:
                    # Check that there's at least two points in the new line segment
                    if len(new_x) == 0:
                        new_x += [x[i - 1], x[i]]
                        new_y += [y[i - 1], y[i]]
                    elif abs(new_x[-1] - x[i - 1]) < 0.0001:
                        new_x.append(x[i])
                        new_y.append(y[i])
                    else:
                        new_data.append([new_x, new_y, line[2]])
                        new_x = [x[i - 1], x[i]]
                        new_y = [y[i - 1], y[i]]
            #             new_data.append([x[start_idx:], y[start_idx:], line[2]])
            #             if i - start_idx > 1:
            #                 new_data.append([x[start_idx:i - 1], y[start_idx:i - 1], line[2]])
            #             start_idx = i
            #         elif len(x) == 2:
            #             new_data.append([x, y, line[2]])
            # if i - start_idx > 1:
            if len(new_x) > 0:
                new_data.append([new_x, new_y, line[2]])

        return new_data

    def update_tension(self):
        """
        Callback for spline tension slider
        """
        self.tension_label.setText(str(self.tension_slider.value()))
        model = self.plot_drift()
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
        self.update_delta_model(obstreeloop.drift_method, model)

    @staticmethod
    def calc_none_dg(data, loop_name):
        """
        Calculates delta-g's from successive gravity observations
        :param data: list of stations from which to calculate delta-g
        :return: PyQt DeltaTableModel
        """
        first = True
        delta_model = DeltaTableModel()
        for station in data:
            if first:
                first = False
                prev_station = station
                continue
            delta = Delta(prev_station,
                          station,
                          driftcorr=0.0,
                          loop=loop_name)
            delta_model.insertRows(delta, 0)
            prev_station = station
        return delta_model

    def calc_netadj_dg(self, data, loop_name):
        """
        Calculates delta-g's from successive gravity observations
        :param data: list of stations from which to calculate delta-g
        :param loop_name: stored with Delta object, used later in network adjustment
        :return: PyQt DeltaTableModel
        """
        first = True
        delta_model = DeltaTableModel()
        for station in data:
            if first:
                first = False
                prev_station = station
                continue
            delta = Delta(prev_station,
                          station,
                          driftcorr=0.0,
                          ls_drift=(loop_name, self.drift_polydegree_combobox.currentIndex() - 1),
                          loop=loop_name)
            delta_model.insertRows(delta, 0)
            prev_station = station
        return delta_model

    @staticmethod
    def calc_cont_dg(xp, yp, data, loop_name):
        """
        Calculates delta-g's while removing drift using the input drift model
        :param xp: times of continuous drift model
        :param yp: continuous drift model
        :param data: plot_data list
        :return: PyQt DeltaTableModel
        """
        delta_model = DeltaTableModel()
        first = True
        ypsum = [0]
        delta_list = []
        for x, drift_rate in zip(xp, yp):
            if first:
                first = False
                prev_x = x
            else:
                prev_sum = ypsum[-1]
                interval = (x - prev_x) * 24
                prev_x = x
                ypsum.append(prev_sum + drift_rate * interval)

        xp = xp.tolist()
        yp = ypsum  # yp = yp.tolist()
        first = True
        for station in data:
            if first:
                first = False
                prev_station = station
                continue
            drift1_idx = min(range(len(xp)), key=lambda i: abs(xp[i] - prev_station.tmean))
            drift1 = yp[drift1_idx]
            drift2_idx = min(range(len(xp)), key=lambda i: abs(xp[i] - station.tmean))
            drift2 = yp[drift2_idx]
            delta = Delta(prev_station,
                          station,
                          driftcorr=drift2 - drift1,
                          loop=loop_name)
            delta_list.append(delta)
            delta_model.insertRows(delta, 0)
            prev_station = station
        return delta_model

    @staticmethod
    def calc_roman_dg(data, loop_name, time_threshold=None):
        """
        Caculates delta-g between three station occupations (one station visited once, one station visited twice) by
        interpolating drift at the latter station.

        Accommodating the time threshold is tricky. for the plotting to be correct the initial g subtracted from
        each measurement has to vary.
    
        :param data: list of stations
        :return: tuple with 2 pyqt models (for dg samples and average dg) and plot data for vertical lines
        """
        # assumes stations in data are in chronological order

        roman_dg_model = RomanTableModel()
        vert_lines = []
        station_list = [i.station_name for i in data]
        unique_stations = list(set(station_list))

        # store initial value at each station
        initial_g = dict()

        # Easiest case: all g values are relative to the initial g at that station
        if time_threshold is None:
            for station_name in unique_stations:
                for station in data:
                    if station.station_name == station_name:
                        initial_g[station_name] = station.gmean
                        break
        else:
            # If time_threshold is specified (checked in the GUI) we need to build a list of possible initial g values
            # for each station. Possible values are those occurring after a gap >= time_threshold (i.e., if there is a
            # gap, reset the initial g that's subtracted from the measurements. This makes the lines start at y = 0 on
            # the plots.
            #
            # Builds the dictionary:
            # initial_g{key:station_name value:(time, g)}
            for station_name in unique_stations:
                stations = []
                for station in data:
                    if station.station_name == station_name:
                        stations.append(station)
                iter_stations = iter(stations)
                first_station = next(iter_stations)
                initial_xy = [(first_station.tmean, first_station.gmean)]
                for station in iter_stations:
                    if (station.tmean - first_station.tmean) * 1440 > time_threshold:
                        initial_xy.append((station.tmean, station.gmean))
                    first_station = station
                initial_g[station_name] = initial_xy

        # For each station in data, loop over all the other stations looking for two observations that bracket the
        # first station
        for station in data:
            for other_station in unique_stations:
                # Ignore it if its the same station
                if other_station == station.station_name:
                    continue
                else:
                    # get all occurrences of the other station
                    other_stations = [i for i in data if i.station_name == other_station]
                    if len(other_stations) > 1:
                        iter_stations = iter(other_stations)
                        other1 = next(iter_stations)
                        for other2 in iter_stations:
                            # Check for 3-point configuration (2 observations at other station bracket the initial obs)
                            if other1.tmean < station.tmean < other2.tmean:
                                # Check that time_threshold is met, or not set
                                if time_threshold is None or \
                                        (other2.tmean - other1.tmean) * 1440 < time_threshold:
                                    delta = Delta(station,
                                                  (other1, other2),
                                                  delta_type='three_point',
                                                  loop=loop_name)
                                    sta2_dg = other2.gmean - other1.gmean
                                    # this is the drift correction
                                    time_prorate = (station.tmean - other1.tmean) / (
                                            other2.tmean - other1.tmean)
                                    # Look for previous occupation at same station. If there is a break > time_threshold
                                    # between the previous and current occupation, we need to account for the shift in
                                    # initial g. Each station has a unique initial g (that might change,
                                    # depending on the time_threshold).
                                    if time_threshold is not None:
                                        initial_gees = initial_g[other1.station_name]
                                        other_initial_g = initial_gees[0][1]
                                        if len(initial_gees) > 1:
                                            for initial_xy in initial_gees[1:]:
                                                if other1.tmean >= initial_xy[0]:
                                                    other_initial_g = initial_xy[1]
                                        initial_gees = initial_g[station.station_name]
                                        station_initial_g = initial_gees[0][1]
                                        if len(initial_gees) > 1:
                                            for initial_xy in initial_gees[1:]:
                                                if station.tmean >= initial_xy[0]:
                                                    station_initial_g = initial_xy[1]
                                    # Easy case: everything relative to the initial observation.
                                    else:
                                        other_initial_g = initial_g[other1.station_name]
                                        station_initial_g = initial_g[station.station_name]

                                    vert_lines.append([(station.tmean, station.tmean),
                                                       ((other_initial_g -
                                                         other1.gmean -
                                                         (sta2_dg * time_prorate)) * -1,
                                                        station.gmean - station_initial_g)])

                                    roman_dg_model.insertRows(delta, 0)
                            other1 = other2

        # If there is more than one delta-g between a given station pair, average them
        # Setup dict to store averages '(sta1, sta2)':[g]
        avg_dg = dict()
        unique_pairs = set()
        for i in range(roman_dg_model.rowCount()):
            delta = roman_dg_model.data(roman_dg_model.index(i, 0), role=QtCore.Qt.UserRole)
            delta_key1 = (delta.station1.station_name, delta.station2[0].station_name)
            delta_key2 = (delta.station2[0].station_name, delta.station1.station_name)
            if delta_key1 not in unique_pairs and delta_key2 not in unique_pairs:
                unique_pairs.add(delta_key1)
                avg_dg[delta_key1] = [delta]
                for ii in range(i + 1, roman_dg_model.rowCount()):
                    testdelta = roman_dg_model.data(roman_dg_model.index(ii, 0), role=QtCore.Qt.UserRole)
                    testdelta_key1 = (testdelta.station1.station_name, testdelta.station2[0].station_name)
                    testdelta_key2 = (testdelta.station2[0].station_name, testdelta.station1.station_name)
                    if delta_key1 == testdelta_key1 or delta_key1 == testdelta_key2:
                        avg_dg[delta_key1].append(testdelta)

        roman_avg_dg_model = DeltaTableModel()
        for station_pair in avg_dg.items():
            # just send list of deltas, not key (station info is already in the deltas)
            avg_delta = Delta.from_list(station_pair[1])
            roman_avg_dg_model.insertRows(avg_delta, 0)

        return roman_dg_model, roman_avg_dg_model, vert_lines

    @staticmethod
    def plot_tares(axes, obstreeloop):
        """
        Plots a vertical line at the time of a tare
        """
        ylim = axes.get_ylim()
        if obstreeloop.tare_model:
            for i in range(obstreeloop.tare_model.rowCount()):
                idx = obstreeloop.tare_model.createIndex(i, 0)
                tare = obstreeloop.tare_model.data(idx, role=QtCore.Qt.UserRole)
                tm = tare.datetime.time()
                x_time = tare.datetime.toordinal() + tm.hour / 24 + tm.minute / 1440 + tm.second / 86400
                axes.plot([x_time, x_time], [ylim[0], ylim[1]], 'gray')
                axes.set_ylim(ylim)
                axes.figure.canvas.draw()

    def clear_axes(self):
        """
        Clears plot axes
        """
        self.axes_drift_single.cla()
        self.axes_drift_cont_lower.clear()
        self.axes_drift_cont_upper.clear()
        self.drift_single_canvas.draw()
        self.drift_cont_canvasbot.draw()
        self.drift_cont_canvastop.draw()

    def plot_drift(self, obstreeloop=None, update=True):
        """
        Catch-all function to plot drift
        :param obstreeloop: Can either specify a loop, or by default use the currentLoopIndex.
        """
        # TODO: plotting and calculating delta-gs is entertwined. To be efficient when loading many loops,
        # I use update to indicate lines that are run only if the plots are visible. If the plotting and
        # delta-g code were better separated, update wouldn't be needed.
        offset = 0
        if type(obstreeloop) is not ObsTreeLoop:
            obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
            obstreesurvey = obstreeloop.parent()
        drift_type = obstreeloop.drift_method
        plot_data = obstreeloop.get_data_for_plot()

        # Check that there's station repeats. If there isn't, skip the plotting but we still want to calculate
        # delta-g's (except for Roman correction).
        no_data = True
        if any([True for x in plot_data if len(x[0]) > 1]):
            no_data = False

        data = obstreeloop.checked_stations()

        # Only include drift observations that meet time criteria
        time_threshold = None
        if self.drift_screen_elapsed_time.isChecked():
            hour = self.drift_time_spinner.dateTime().time().hour()
            minute = self.drift_time_spinner.dateTime().time().minute()
            time_threshold = hour * 60 + minute
            plot_data = self.screen_for_elapsed_time(plot_data, elapsed_time=time_threshold)
        if drift_type == 'none' or drift_type == 'netadj':
            # none, netadj, and roman all use axes_drift_single
            delta_model = None
            if update:
                self.axes_drift_single.cla()
            logging.info('Plotting drift - no correction, Loop ' + obstreeloop.name)
            # Get data for plotting
            for line in plot_data:
                if len(line[0]) > 1:
                    # Make values relative to first station value
                    y = [f - line[1][0] + offset for f in line[1]]
                    x = [f for f in line[0]]
                    if update:
                        a = self.axes_drift_single.plot(x, y, '.-', picker=5)
                        a[0].name = line[2]
                        offset += self.offset_slider.value()

            # Plot
            if plot_data and not no_data and update:
                self.axes_drift_single.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                self.axes_drift_single.yaxis.set_label_text('Change in gravity since initial \nstation occupation, ' +
                                                            'in microGal')
                self.drift_fig.canvas.mpl_connect('pick_event',
                                                  lambda event: self.show_line_label(event, self.axes_drift_single))
                self.plot_tares(self.axes_drift_single, obstreeloop)
            elif update:
                self.axes_drift_single.cla()
                self.axes_drift_single.text(0.35, 0.5, 'NO STATION REPEATS')
            if update:
                self.axes_drift_single.set_title('Survey ' + obstreesurvey.name + ', Loop ' + obstreeloop.name)
                self.drift_single_canvas.draw()

            if drift_type == 'none':
                delta_model = self.calc_none_dg(data, obstreeloop.name)
            elif drift_type == 'netadj':
                delta_model = self.calc_netadj_dg(data, obstreeloop.name)
            return delta_model

        if drift_type == 'continuous':
            logging.info('Plotting continuous drift, Loop ' + obstreeloop.name)
            self.axes_drift_cont_lower.clear()
            self.axes_drift_cont_upper.clear()
            drift_rate, drift_time, drift_x = [], [], []

            # Get data for plotting
            min_time = 100000000
            max_time = 0
            for line in plot_data:
                # x and y are the time and g values for each station.
                # Make values relative to first station value
                y = [f - line[1][0] + offset for f in line[1]]
                x = [f for f in line[0]]
                if min(x) < min_time:
                    min_time = min(x)
                if max(x) > max_time:
                    max_time = max(x)
                # Only bother plotting if there's more than one station (don't do this otherwise, otherwise singleton
                # stations at the start or end of a survey won't be included when setting the min_time/max_time
                if len(line[0]) > 1:
                    # Loop over the line vertices
                    for idx, obs in enumerate(y):
                        y[idx] = obs + offset
                        # get drift rate for bottom plot
                        if idx >= 1:
                            dr = (y[idx] - y[idx - 1]) / ((x[idx] - x[idx - 1]) * 24)  # drift rate
                            drift_rate.append(dr)
                            xmean = np.mean([x[idx], x[idx - 1]])
                            drift_x.append(xmean)
                            drift_time.append(dt.datetime.utcfromtimestamp((xmean - 719163) * 86400.))
                            # Plot horizontal extent
                            if self.drift_plot_hz_extent.isChecked() and update:
                                self.axes_drift_cont_lower.plot([x[idx], x[idx - 1]], [dr, dr], '-', color='0.5')
                    if update:
                        a = self.axes_drift_cont_upper.plot(x, y, '.-', picker=5)
                        a[0].name = line[2]
                        offset += self.offset_slider.value()

            # Plot
            if plot_data:
                if update:
                    self.axes_drift_cont_upper.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                    self.axes_drift_cont_lower.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                    self.axes_drift_cont_lower.plot(drift_time, drift_rate, '.', picker=2)
                    xticks = self.axes_drift_cont_upper.get_xticks()
                    self.axes_drift_cont_lower.set_xticks(xticks)
                    xlims = self.axes_drift_cont_upper.get_xlim()
                    self.axes_drift_cont_lower.set_xlim(xlims)
                    self.axes_drift_cont_lower.yaxis.set_label_text('Drift rate,\nin microGal/hr')
                    self.axes_drift_cont_upper.yaxis.set_label_text('Drift, in microGal\n(arbitrary offset)')
                    self.drift_cont_figtop.canvas.mpl_connect('pick_event',
                                                              lambda event: self.show_line_label
                                                              (event, self.axes_drift_cont_upper))
                    self.drift_cont_figbot.canvas.mpl_connect('pick_event', self.drift_point_picked)
                    self.drift_cont_figbot.canvas.mpl_connect('button_release_event', self.drift_newpoint_picked)

                # Interpolate drift model: polynomial, spline, etc. at xp number of points. xp needs to remain
                # relatively low to maintain performance.
                xp = np.linspace(min(drift_x), max(drift_x), 50)  # constant
                method_key = self.drift_polydegree_combobox.currentIndex()
                if method_key == 0:  # no drift correction
                    yp = np.zeros(xp.size)
                else:
                    x0 = [f - np.min(drift_x) for f in drift_x]
                    xp0 = [f - np.min(xp) for f in xp]
                    idx = sorted(range(len(x0)), key=lambda xpt: x0[xpt])
                    x_sorted, drift_rate_sorted = [], []
                    for i in idx:
                        x_sorted.append(x0[i])
                        drift_rate_sorted.append(drift_rate[i])
                    x0 = x_sorted
                    drift_rate = drift_rate_sorted
                    if method_key == 9:
                        pass
                    if method_key == 1:  # spline
                        try:
                            s = UnivariateSpline(x0, drift_rate, k=3, s=self.tension_slider.value())
                            xs = np.linspace(x0[0], x0[-1], 50)
                            yp = s(xs)
                            logging.info('Spline drift correction, tension={}'.format(self.tension_slider.value()))
                        except:
                            self.msg = show_message('Insufficient drift observations for spline method', 'Error')
                    elif method_key == 2:  # 1st order poly
                        z = np.polyfit(x0, drift_rate, 1)
                        p = np.poly1d(z)
                        yp = p(xp0)
                        logging.info('First-order poly drift correction')
                    elif method_key == 3:  # 2nd order poly
                        z = np.polyfit(x0, drift_rate, 2)
                        p = np.poly1d(z)
                        yp = p(xp0)
                        logging.info('Second-order poly drift correction')
                    elif method_key == 4:  # 3rd order poly
                        z = np.polyfit(x0, drift_rate, 3)
                        p = np.poly1d(z)
                        yp = p(xp0)
                        logging.info('Third-order poly drift correction')

                # Method for extrapolating beyond fitted drift curve extene
                if self.drift_cont_startendcombobox.currentIndex() == 1:  # constant
                    new_xp = np.linspace(min_time, min(drift_x), 20)
                    new_xp = np.append(new_xp, xp)
                    new_xp = np.append(new_xp, np.linspace(max(drift_x), max_time, 20))
                    xp = new_xp
                    new_yp = np.ones(20) * yp[0]
                    new_yp = np.append(new_yp, yp)
                    new_yp = np.append(new_yp, np.ones(20) * yp[-1])
                    yp = new_yp
                else:  # linear extrapolation from first two (and last two) points
                    # get first two points
                    x = xp[:2]
                    y = yp[:2]
                    z = np.polyfit(x, y, 1)
                    p = np.poly1d(z)
                    new_xp1 = np.linspace(min_time, min(drift_x), 20)
                    yp1 = p(new_xp1)
                    x = xp[-2:]
                    y = yp[-2:]
                    z = np.polyfit(x, y, 1)
                    p = np.poly1d(z)
                    new_xp2 = np.linspace(max(drift_x), max_time, 20)
                    yp2 = p(new_xp2)
                    xp_temp = np.append(new_xp1, xp)
                    xp = np.append(xp_temp, new_xp2)
                    yp_temp = np.append(yp1, yp)
                    yp = np.append(yp_temp, yp2)
                delta_model = self.calc_cont_dg(xp, yp, data, obstreeloop.name)
                if update:
                    self.plot_tares(self.axes_drift_cont_lower, obstreeloop)
                    self.plot_tares(self.axes_drift_cont_upper, obstreeloop)
                    self.axes_drift_cont_lower.plot(xp, yp, 'k-')
                    self.axes_drift_cont_upper.set_title('Survey ' + obstreesurvey.name + ', Loop ' + obstreeloop.name)
                    self.drift_cont_canvasbot.draw()
                    self.drift_cont_canvastop.draw()
                return delta_model
            else:
                self.msg = show_message('No data available for plotting', 'Plot error')

        # Plots vertical dashed lines showing delta-g's
        elif drift_type == 'roman':
            logging.info('Plotting Roman drift, Loop ' + obstreeloop.name)
            if update:
                self.axes_drift_single.cla()
            models = self.calc_roman_dg(data, obstreeloop.name, time_threshold)

            for line in plot_data:
                if len(line[0]) > 1:
                    # Make values relative to first station value
                    y = [f - line[1][0] for f in line[1]]
                    # TODO: store dates in station object in datetime format, to avoid this conversion?
                    x = [dt.datetime.utcfromtimestamp((f - 719163) * 86400.) for f in line[0]]
                    a = self.axes_drift_single.plot(x, y, '.-', picker=5)
                    a[0].name = line[2]

            for line in models[2]:
                if update:
                    self.axes_drift_single.plot(line[0], line[1], '--')
            if plot_data and update:
                self.axes_drift_single.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            if update:
                self.axes_drift_single.yaxis.set_label_text('Change in gravity since initial \nstation occupation, ' +
                                                            'in microGal')
                self.drift_fig.canvas.mpl_connect('pick_event',
                                                  lambda event: self.show_line_label(event, self.axes_drift_single))
                self.axes_drift_single.set_title('Survey ' + obstreesurvey.name + ', Loop ' + obstreeloop.name)
                self.drift_single_canvas.draw()
            return models

    @staticmethod
    def show_all_columns(delta_view):
        """
        Helper function to reset columns in the delta_model
        :param delta_view: view to reset
        """
        model = delta_view.model()
        for i in range(model.columnCount()):
            delta_view.showColumn(i)

    def set_drift_method(self, update=True):
        """
        Called from update_drift_tables_and_plots + callback from GUI. Initiates plotting on drift tab.
        :param update: Boolean or int, controls if plots are updated. For performance, it's set to false when loading a file
        """
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
        method_key = self.driftmethod_comboboxbox.currentIndex()

        tare_model = obstreeloop.tare_model
        if tare_model:
            tare_model.dataChanged.connect(self.update_tares)
        self.tare_view.setModel(tare_model)

        inv_drift_lookup = {v: k for k, v in self.parent.drift_lookup.items()}
        method = inv_drift_lookup[method_key]
        logging.info('Drift method set to ' + method)
        orig_method = obstreeloop.drift_method
        obstreeloop.drift_method = method

        # These control the visibility of different tables
        # update is an int (index of menu item) when this function is called from the
        #   menu-item callback
        if type(update) is int or update is True:
            if method == 'none':
                self.drift_none()
            if method == 'netadj':
                self.drift_polydegree_combobox.setCurrentIndex(obstreeloop.drift_netadj_method)
                self.drift_adjust()
            if method == 'roman':
                self.drift_roman()
            if method == 'continuous':
                self.drift_polydegree_combobox.setCurrentIndex(obstreeloop.drift_cont_method)
                self.drift_cont_startendcombobox.setCurrentIndex(obstreeloop.drift_cont_startend)
                self.drift_continuous()
        model = self.plot_drift(update=update)

        if method != orig_method or self.parent.workspace_loaded:
            # Leave the status bar hollow if a workspace was loaded
            if not self.parent.workspace_loaded:
                self.parent.deltas_update_required()
            self.parent.workspace_loaded = False
        self.update_delta_model(method, model)

    def update_delta_model(self, method, model):
        """
        Show appropriate delta model for the selected loop.

        This method takes the model generated by plot_drift() and assigns it to the delta_view on the drift tab.
        :param method: If 'roman', show sample and average models. Otherwise, show a single model.
        :param model: a PyQt model or list of models (Roman method).
        """
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
        if method == 'roman':
            models = model
            self.dg_samples_view.setModel(models[0])
            obstreeloop.delta_model = models[1]
            self.delta_view.setModel(models[1])
            # Hide drift correction, std_for_adj, and residual columns
            self.show_all_columns(self.delta_view)
            self.delta_view.hideColumn(2)
            self.delta_view.hideColumn(5)
            self.delta_view.hideColumn(7)
            self.delta_view.hideColumn(8)
        else:
            obstreeloop.delta_model = model
            self.delta_view.setModel(model)
            # Hide std_for_adj and residual columns
            self.show_all_columns(self.delta_view)
            self.delta_view.hideColumn(2)
            self.delta_view.hideColumn(7)
            self.delta_view.hideColumn(8)

    def tare_context_menu(self, point):
        """
        Right-click context menu on tare table
        :param point: PyQt reference to click point, determines where to show popup.
        """
        selected = self.tare_view.selectedIndexes()
        if selected:
            self.tare_popup_menu.addAction(self.mnDeleteTare)
            self.tare_popup_menu.exec_(self.tare_view.mapToGlobal(point))

    #     self.tarepopMenu = QtWidgets.QMenu("Menu", self)
    #     self.tarepopMenu.addAction(self.mnDeleteTare)
    #     self.tarepopMenu.exec_(self.tare_view.mapToGlobal(point))
    #
    # def datum_context_menu(self, point):
    #     selected = self.datum_view.selectedIndexes()
    #     if selected:
    #         self.datum_popup_menu.addAction(self.mnDeleteDatum)
    #         self.datum_popup_menu.exec_(self.datum_view.mapToGlobal(point))

    # def init_tares_popup_menu(self):
    #     """
    #     Right-click menu on results table
    #     """
    #     self.mnDeleteTare = QtWidgets.QAction('Delete tare', self)
    #     self.mnDeleteTare.triggered.connect(self.delete_tare)
    #
    # def delete_tare(self, idxs):
    #     obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.currentLoopIndex)
    #     idxs = self.tare_view.selectedIndexes()
    #     if len(idxs) > 1:
    #         return
    #     else:
    #         obstreeloop.tare_model.removeRow(idxs[0].row(), idxs[0].parent())
    #         obstreeloop.tare_model.deleteTare(idxs[0])
    #         # tare = obstreeloop.tare_model.data(idxs[0], role=QtCore.Qt.UserRole)
    #         # obstreeloop.tare_mode.deleteTare(idxs[0])
    #     return

    @staticmethod
    def process_tares(obstreeloop):
        """
        Apply tares in the tare table.
        :param obstreeloop: Loop shown in drift tab
        """

        for i in range(obstreeloop.rowCount()):
            obstreestation = obstreeloop.child(i)
            for idx, t in enumerate(obstreestation.t):
                obstreestation.tare[idx] = 0
            for tare in obstreeloop.tare_model.tares:
                if tare.checked == 2:
                    qdt = QtCore.QDateTime(tare.date, tare.time)
                    tare_dt = date2num(qdt.toPyDateTime())
                    for idx, t in enumerate(obstreestation.t):
                        if t > tare_dt:
                            obstreestation.tare[idx] += tare.tare

    def update_tares(self, selected):
        """
        Respond to tare check/uncheck.
        """
        if hasattr(selected, 'model'):
            obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
            self.process_tares(obstreeloop)
            model = self.plot_drift()
            method = obstreeloop.drift_method
            self.update_delta_model(method, model)

    def drift_adjust(self):
        """
        Update which PyQt tables are shown
        """
        self.drift_single_canvas.show()
        self.drift_single_canvas.setMinimumWidth(700)
        self.drift_cont_plotpanel.hide()
        self.cont_label_widget.show()
        self.cont_label_widget.setMinimumHeight(50)
        self.roman_label_widget.hide()
        self.drift_window.setMinimumHeight(200)
        self.tension_slider.setEnabled(False)
        # Disable the 'none' and 'spline' options, they're not relevant
        self.drift_polydegree_combobox.model().item(0).setEnabled(False)
        self.drift_polydegree_combobox.model().item(1).setEnabled(False)
        self.drift_polydegree_combobox.setEnabled(True)
        self.drift_cont_startendcombobox.setEnabled(False)
        self.drift_plot_hz_extent.setEnabled(False)
        # self.delta_view.show()
        self.dg_samples_view.hide()

    def drift_roman(self):
        """
        Update which PyQt tables are shown
        """
        # self.plot_drift()
        self.roman_label_widget.show()
        self.drift_single_canvas.show()
        self.drift_single_canvas.setMinimumWidth(700)
        self.drift_cont_plotpanel.hide()
        self.dg_samples_view.show()
        # self.delta_view.show()
        self.cont_label_widget.hide()
        self.roman_label_widget.show()
        self.roman_label_widget.setMinimumHeight(50)
        self.drift_window.setMinimumHeight(200)
        self.tension_slider.setEnabled(False)
        self.offset_slider.setEnabled(False)
        self.drift_plot_hz_extent.setEnabled(False)
        self.drift_cont_startendcombobox.setEnabled(False)
        self.drift_polydegree_combobox.setEnabled(False)

    def drift_continuous(self):
        """
        Update which PyQt tables are shown
        """
        self.drift_single_canvas.hide()
        self.drift_cont_plotpanel.show()
        self.drift_cont_plotpanel.setMinimumWidth(700)
        self.dg_samples_view.hide()
        # Hide std_for_adj and residual columns
        self.show_all_columns(self.delta_view)
        self.delta_view.hideColumn(8)
        self.delta_view.hideColumn(9)
        self.cont_label_widget.show()
        self.cont_label_widget.setMinimumHeight(50)
        self.roman_label_widget.hide()
        self.drift_window.setMinimumHeight(200)
        # Re-enable these options (they're disabled if netadj drift was selected)
        self.drift_polydegree_combobox.model().item(0).setEnabled(True)
        self.drift_polydegree_combobox.model().item(1).setEnabled(True)
        self.tension_slider.setEnabled(True)
        self.offset_slider.setEnabled(True)
        self.drift_plot_hz_extent.setEnabled(True)
        self.drift_cont_startendcombobox.setEnabled(True)
        self.drift_polydegree_combobox.setEnabled(True)

    def drift_none(self):
        """
        Update which PyQt tables are shown
        """
        self.drift_single_canvas.show()
        self.drift_single_canvas.setMinimumWidth(700)
        self.drift_cont_plotpanel.hide()
        self.cont_label_widget.show()
        self.cont_label_widget.setMinimumHeight(50)
        self.roman_label_widget.hide()
        self.drift_window.setMinimumHeight(200)
        self.tension_slider.setEnabled(False)
        self.offset_slider.setEnabled(True)
        self.drift_polydegree_combobox.setEnabled(False)
        self.drift_cont_startendcombobox.setEnabled(False)
        self.drift_plot_hz_extent.setEnabled(False)
        self.dg_samples_view.hide()

    def drift_combobox_updated(self):
        """
        Called when either the continuous method or extrapolate/constant combobox is changed.
        """
        method_key = self.drift_polydegree_combobox.currentIndex()
        startend_key = self.drift_cont_startendcombobox.currentIndex()
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(self.parent.index_current_loop)
        drift_method = obstreeloop.drift_method
        if drift_method == 'continuous':
            obstreeloop.drift_cont_method = method_key
            obstreeloop.drift_cont_startend = startend_key
            if method_key == 1:
                self.tension_slider.setEnabled(True)
            else:
                self.tension_slider.setEnabled(False)
        elif drift_method == 'netadj':
            obstreeloop.drift_netadj_method = method_key

        model = self.plot_drift()
        self.update_delta_model(drift_method, model)
