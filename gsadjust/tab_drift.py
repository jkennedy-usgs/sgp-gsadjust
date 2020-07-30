#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tab_drift.py
===============

PyQt graphical elements on the drift tab of GSadjust.
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
import logging

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.dates import DateFormatter, date2num
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt

from data_objects import Delta
from drift_continuous import drift_continuous
from drift_roman import drift_roman
from gui_objects import IncrMinuteTimeEdit, show_message
from pyqt_models import (
    DeltaTableModel,
    NoCheckDeltaTableModel,
    ObsTreeLoop,
    RomanTableModel,
    TareTableModel,
)


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
        main_vsplitter_window = QtWidgets.QSplitter(Qt.Vertical, self)

        # Setup drift figures (Roman and Continuous). Only one will be shown at a time.
        # Drift figure: default and roman
        self.drift_window = QtWidgets.QSplitter(Qt.Horizontal, self)
        self.drift_fig = Figure((3.0, 5.0), dpi=self.dpi, facecolor='white')
        self.drift_single_canvas = FigureCanvas(self.drift_fig)
        self.drift_fig.subplots_adjust(wspace=0.3)
        self.axes_drift_single = self.drift_fig.add_subplot(111)

        # Drift figure: continuous
        # Plot panel
        self.drift_cont_plotpanel = QtWidgets.QSplitter(Qt.Vertical, self)
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
        self.driftmethod_combobox_key = {
            0: "None",
            1: "Network adjustment",
            2: "Roman (interpolate)",
            3: "Continuous model",
        }
        self.driftmethod_combobox = QtWidgets.QComboBox()
        self.driftmethod_combobox.activated.connect(self.set_drift_method)
        for item in self.driftmethod_combobox_key.values():
            self.driftmethod_combobox.addItem(item)

        # Widget to remove dg-observations with a long elapsed time in between
        self.drift_screen_elapsed_time = QtWidgets.QCheckBox(
            'Max. time between repeats (hh:mm)'
        )
        self.drift_screen_elapsed_time.setChecked(False)
        self.drift_screen_elapsed_time.stateChanged.connect(self.time_extent_changed)
        self.drift_time_spinner = IncrMinuteTimeEdit(QtCore.QTime(1, 0))
        self.drift_time_spinner.timeChanged.connect(self.time_extent_changed)
        self.drift_time_spinner.setDisplayFormat("hh:mm")

        # Widget to add horizontal-extent lines to drift-rate plot
        self.drift_plot_hz_extent = QtWidgets.QCheckBox(
            'Show time-extent of drift observation'
        )
        self.drift_plot_hz_extent.setChecked(False)
        self.drift_plot_hz_extent.stateChanged.connect(self.plot_drift)

        self.drift_plot_weighted = QtWidgets.QCheckBox('Weight drift observations')
        self.drift_plot_weighted.setChecked(False)
        self.drift_plot_weighted.stateChanged.connect(self.update_weighted)

        self.tension_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.tension_slider.setRange(10, 2500)
        self.tension_slider.setValue(1250)
        self.tension_slider.setEnabled(False)
        self.tension_slider.valueChanged.connect(self.update_tension)
        self.tension_label = QtWidgets.QLabel()
        self.tension_label.setText('{:2.2f}'.format(self.tension_slider.value()))
        self.tension_label.setAlignment(Qt.AlignCenter)
        self.drift_cont_methodcombobox_key = {
            0: "Constant",
            1: "Spline",
            2: "1st order polynomial",
            3: "2nd order polynomial",
            4: "3rd order polynomial",
        }
        self.drift_polydegree_combobox = QtWidgets.QComboBox()
        self.drift_polydegree_combobox.activated.connect(self.drift_combobox_updated)
        for item in self.drift_cont_methodcombobox_key.values():
            self.drift_polydegree_combobox.addItem(item)
        self.drift_cont_behaviorcombobox_key = {0: "Extrapolate", 1: "Constant"}
        self.drift_cont_startendcombobox = QtWidgets.QComboBox()
        self.drift_cont_startendcombobox.activated.connect(self.drift_combobox_updated)
        for item in self.drift_cont_behaviorcombobox_key.values():
            self.drift_cont_startendcombobox.addItem(item)

        self.offset_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.offset_slider.setRange(0, 10)
        self.offset_slider.setValue(0)
        self.offset_slider.valueChanged.connect(self.plot_drift)
        drift_controls = QtWidgets.QWidget()
        drift_cont_control_layout = QtWidgets.QHBoxLayout()
        drift_control_sublayout = QtWidgets.QVBoxLayout()
        grid_widget = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Drift correction method'), 0, 0)
        grid.addWidget(self.driftmethod_combobox, 0, 1)
        grid.addWidget(QtWidgets.QLabel('Drift model type'), 1, 0)
        grid.addWidget(self.drift_polydegree_combobox, 1, 1)
        grid.addWidget(QtWidgets.QLabel('Behavior at start/end:'), 2, 0)
        grid.addWidget(self.drift_cont_startendcombobox, 2, 1)
        grid.addWidget(self.drift_screen_elapsed_time, 3, 0)
        grid.addWidget(self.drift_time_spinner, 3, 1)
        grid.addWidget(self.drift_plot_hz_extent, 4, 0)
        grid.addWidget(self.drift_plot_weighted, 5, 0)
        grid.addWidget(QtWidgets.QLabel('Vertical line offset'), 6, 0)
        grid.addWidget(self.offset_slider, 6, 1)
        grid.addWidget(QtWidgets.QLabel('Spline tension:'), 7, 0)
        grid.addWidget(self.tension_slider, 7, 1)
        grid.addWidget(self.tension_label, 7, 2)

        grid_widget.setLayout(grid)
        drift_control_sublayout.addWidget(grid_widget)

        self.tare_view = QtWidgets.QTableView()
        self.tare_view.clicked.connect(self.update_tares)
        self.tare_view.setContextMenuPolicy(Qt.CustomContextMenu)
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
        main_hsplitter_window = QtWidgets.QSplitter(Qt.Horizontal, self)
        main_hsplitter_window.addWidget(self.dg_samples_view)
        main_hsplitter_window.addWidget(self.delta_view)
        main_hsplitter_window.setMinimumHeight(300)
        main_vsplitter_window.addWidget(main_hsplitter_window)
        self.delta_view.show()
        self.dg_samples_view.hide()

        layout_main.addWidget(main_vsplitter_window)
        self.setLayout(layout_main)

    def reset(self):
        self.driftmethod_combobox.setCurrentIndex(0)
        self.drift_polydegree_combobox.setCurrentIndex(0)
        self.axes_drift_single.cla()
        self.axes_drift_cont_lower.clear()
        self.axes_drift_cont_upper.clear()
        self.axes_drift_cont_upper.figure.canvas.draw()
        self.axes_drift_cont_lower.figure.canvas.draw()
        self.axes_drift_single.figure.canvas.draw()
        self.delta_view.setModel(DeltaTableModel())
        self.dg_samples_view.setModel(DeltaTableModel())
        self.tare_view.setModel(TareTableModel())

    def time_extent_changed(self):
        """
        Called when "Max time between repeats" is checked/unchecked
        """
        self.parent.deltas_update_required()
        self.plot_drift()

    # This section provides the right-click context menu in the continuous drift lower plot - not implemented
    # def drift_newpoint_picked(self, event):
    #     if event.button == 3:
    #         self.drift_rate_context_menu()
    #
    # def drift_point_picked(self, event):
    #     if event.mouseevent.button == 3:
    #         self.drift_rate_context_menu(from_pick=True)
    #
    # def drift_rate_context_menu(self, from_pick=False):
    #     """
    #     Not functional (other than showing the menu). Should allow points to be excluded, or artificial points added,
    #     to the continuous drift correction.
    #     :param from_pick: Boolean, True if a point was picked
    #     """
    #     if from_pick:
    #         add = QtWidgets.QAction(QtGui.QIcon(""), "Add point to drift model", self,
    #                                 triggered=self.drift_cont_addpoint,
    #                                 enabled=False)
    #         remove = QtWidgets.QAction(QtGui.QIcon(""), "Remove point from model", self,
    #                                    triggered=self.drift_cont_removepoint)
    #         self.popup_menu.addAction(remove)
    #     else:
    #         add = QtWidgets.QAction(QtGui.QIcon(""), "Add point to drift model", self,
    #                                 triggered=self.drift_cont_addpoint)
    #         remove = QtWidgets.QAction(QtGui.QIcon(""), "Remove point from model", self,
    #                                    triggered=self.drift_cont_removepoint,
    #                                    enabled=False)
    #
    #     self.popup_menu.addAction(add)
    #     self.popup_menu.addAction(remove)
    #     cursor = QtGui.QCursor()
    #     self.popup_menu.popup(cursor.pos())
    #
    # def drift_cont_removepoint(self):
    #     pass
    #
    # def drift_cont_addpoint(self):
    #     pass

    def show_line_label(self, event, axes):
        """
        Shows the station name in the upper left of the drift plot when a 
        line is clicked.
        :param event: Matplotlib event
        :param axes: Current axes (differs for none|netadj|roman vs continuous)
        """
        thisline = event.artist
        if self.station_label is not None:
            self.station_label.set_text('')
        self.station_label = axes.text(
            0.05,
            0.95,
            thisline.name,
            horizontalalignment='center',
            verticalalignment='center',
            transform=axes.transAxes,
        )
        axes.figure.canvas.draw()

    @staticmethod
    def screen_for_elapsed_time(plot_data, elapsed_time):
        """
        We may want to exclude repeat observations with a lot of elapsed time 
        between occupations.
        :param plot_data: input data (list of lists)
        :param elapsed_time: maximum time to be considered a repeat, in minutes
        :return: list with same format as plot_data
        """
        new_data = []
        for line in plot_data:
            x = [x for x in line[0]]
            y = [y for y in line[1]]
            new_x, new_y = [], []
            i = 0
            for i in range(1, len(x)):
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
            if len(new_x) > 0:
                new_data.append([new_x, new_y, line[2]])

        return new_data

    def update_weighted(self):
        """
        Callback for weight drift observations
        """
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_loop
        )
        if obstreeloop:
            obstreeloop.drift_cont_weighting = self.drift_plot_weighted.checkState()
            model = self.plot_drift()
            if model:
                self.update_delta_model(obstreeloop.drift_method, model)

    def update_tension(self):
        """
        Callback for spline tension slider
        """
        self.tension_label.setText(str(self.tension_slider.value()))
        model = self.plot_drift()
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_loop
        )
        self.update_delta_model(obstreeloop.drift_method, model)

    @staticmethod
    def calc_none_dg(data, loop_name):
        """
        Calculates delta-g's from successive gravity observations
        :param data: list of stations from which to calculate delta-g
        :return: PyQt DeltaTableModel
        """
        delta_model = DeltaTableModel()
        # Take the first station from the list, or None if there aren't any.
        prev_station = data.pop(0) if data else None

        for station in data:
            delta = Delta(prev_station, station, driftcorr=0.0, loop=loop_name)
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
        delta_model = DeltaTableModel()
        # Take the first station from the list, or None if there aren't any.
        prev_station = data.pop(0) if data else None

        for station in data:
            delta = Delta(
                prev_station,
                station,
                driftcorr=0.0,
                ls_drift=(loop_name, self.drift_polydegree_combobox.currentIndex() - 1),
                loop=loop_name,
            )
            delta_model.insertRows(delta, 0)
            prev_station = station
        return delta_model

    @staticmethod
    def populate_continuous_deltamodel(deltas):
        """
        Calculates delta-g's while removing drift using the input drift model
        :param xp: times of continuous drift model
        :param yp: continuous drift model
        :param data: plot_data list
        :return: PyQt DeltaTableModel
        """
        delta_model = DeltaTableModel()
        for delta in deltas:
            delta_model.insertRows(delta, 0)
        return delta_model

    @staticmethod
    def calc_roman_dg(data, loop_name, time_threshold=None):
        """
        Caculates delta-g between three station occupations (one station visited 
        once, one station visited twice) by interpolating drift at the latter station.

        Accommodating the time threshold is tricky. for the plotting to be 
        correct the initial g subtracted from each measurement has to vary.
    
        :param data: list of stations
        :return: tuple with 2 pyqt models (for dg samples and average dg) and 
            plot data for vertical lines
        """
        # assumes stations in data are in chronological order
        roman_dg_model = RomanTableModel()
        deltas, vert_lines = drift_roman(data, loop_name, time_threshold=None)

        for delta in deltas:
            roman_dg_model.insertRows(delta, 0)
        # If there is more than one delta-g between a given station pair, average them
        # Setup dict to store averages '(sta1, sta2)':[g]
        avg_dg = dict()
        unique_pairs = set()
        for i in range(roman_dg_model.rowCount()):
            delta = roman_dg_model.data(roman_dg_model.index(i, 0), role=Qt.UserRole)
            delta_key1 = (delta.station1.station_name, delta.station2[0].station_name)
            delta_key2 = (delta.station2[0].station_name, delta.station1.station_name)
            if delta_key1 not in unique_pairs and delta_key2 not in unique_pairs:
                unique_pairs.add(delta_key1)
                avg_dg[delta_key1] = [delta]
                for ii in range(i + 1, roman_dg_model.rowCount()):
                    testdelta = roman_dg_model.data(
                        roman_dg_model.index(ii, 0), role=Qt.UserRole
                    )
                    testdelta_key1 = (
                        testdelta.station1.station_name,
                        testdelta.station2[0].station_name,
                    )
                    testdelta_key2 = (
                        testdelta.station2[0].station_name,
                        testdelta.station1.station_name,
                    )
                    if delta_key1 == testdelta_key1 or delta_key1 == testdelta_key2:
                        avg_dg[delta_key1].append(testdelta)

        roman_avg_dg_model = NoCheckDeltaTableModel()
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
                tare = obstreeloop.tare_model.data(idx, role=Qt.UserRole)
                tm = tare.datetime.time()
                x_time = (
                    tare.datetime.toordinal()
                    + tm.hour / 24
                    + tm.minute / 1440
                    + tm.second / 86400
                )
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
        :param obstreeloop: Can either specify a loop, or by default use the 
            currentLoopIndex.
        """
        # I use update to indicate lines that are run only if the plots are visible. If the plotting and
        # delta-g code were better separated, update wouldn't be needed.
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        offset = 0
        if type(obstreeloop) is not ObsTreeLoop:
            obstreeloop = self.parent.obsTreeModel.itemFromIndex(
                self.parent.index_current_loop
            )
            obstreesurvey = obstreeloop.parent()
        drift_type = obstreeloop.drift_method
        plot_data = obstreeloop.get_data_for_plot()

        self.parent.obsTreeModel.resetStationAsd()

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
            plot_data = self.screen_for_elapsed_time(
                plot_data, elapsed_time=time_threshold
            )
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
                self.axes_drift_single.yaxis.set_label_text(
                    'Change in gravity since initial \nstation occupation, '
                    + 'in microGal'
                )
                self.drift_fig.canvas.mpl_connect(
                    'pick_event',
                    lambda event: self.show_line_label(event, self.axes_drift_single),
                )
                self.plot_tares(self.axes_drift_single, obstreeloop)
            elif update:
                self.axes_drift_single.cla()
                self.axes_drift_single.text(0.35, 0.5, 'NO STATION REPEATS')
            if update:
                self.axes_drift_single.set_title(
                    'Survey ' + obstreesurvey.name + ', Loop ' + obstreeloop.name
                )
                self.drift_single_canvas.draw()

            if drift_type == 'none':
                delta_model = self.calc_none_dg(data, obstreeloop.name)
            elif drift_type == 'netadj':
                delta_model = self.calc_netadj_dg(data, obstreeloop.name)
            QtWidgets.QApplication.restoreOverrideCursor()
            return delta_model

        if drift_type == 'continuous':
            logging.info('Plotting continuous drift, Loop ' + obstreeloop.name)
            self.axes_drift_cont_lower.clear()
            self.axes_drift_cont_upper.clear()
            # Get data for plotting
            min_time = 100000000
            max_time = 0
            drift_rate, drift_time, drift_x = [], [], []
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
                            dr = (y[idx] - y[idx - 1]) / (
                                (x[idx] - x[idx - 1]) * 24
                            )  # drift rate
                            drift_rate.append(dr)
                            xmean = np.mean([x[idx], x[idx - 1]])
                            drift_x.append(xmean)
                            drift_time.append(
                                dt.datetime.utcfromtimestamp((xmean - 719163) * 86400.0)
                            )
                            # Plot horizontal extent
                            if self.drift_plot_hz_extent.isChecked() and update:
                                self.axes_drift_cont_lower.plot(
                                    [x[idx], x[idx - 1]], [dr, dr], '-', color='0.5'
                                )
                    if update:
                        a = self.axes_drift_cont_upper.plot(x, y, '.-', picker=5)
                        a[0].name = line[2]
                        offset += self.offset_slider.value()

            # Plot
            if plot_data:
                if update:
                    self.axes_drift_cont_upper.xaxis.set_major_formatter(
                        DateFormatter('%H:%M')
                    )
                    self.axes_drift_cont_lower.xaxis.set_major_formatter(
                        DateFormatter('%H:%M')
                    )
                    self.axes_drift_cont_lower.plot(
                        drift_time, drift_rate, '.', picker=2
                    )
                    xticks = self.axes_drift_cont_upper.get_xticks()
                    self.axes_drift_cont_lower.set_xticks(xticks)
                    xlims = self.axes_drift_cont_upper.get_xlim()
                    self.axes_drift_cont_lower.set_xlim(xlims)
                    self.axes_drift_cont_lower.yaxis.set_label_text(
                        'Drift rate,\nin microGal/hr'
                    )
                    self.axes_drift_cont_upper.yaxis.set_label_text(
                        'Drift, in microGal\n(arbitrary offset)'
                    )
                    self.drift_cont_figtop.canvas.mpl_connect(
                        'pick_event',
                        lambda event: self.show_line_label(
                            event, self.axes_drift_cont_upper
                        ),
                    )
                    # drift_point_picked and drift_newpoint_picked are for adding/removing points to continuous drift
                    # curve - not yet implemented.
                    # self.drift_cont_figbot.canvas.mpl_connect('pick_event', self.drift_point_picked)
                    # self.drift_cont_figbot.canvas.mpl_connect('button_release_event', self.drift_newpoint_picked)

                try:
                    z = []
                    deltas, xp, yp, z = drift_continuous(
                        data,
                        plot_data,
                        drift_x,
                        drift_rate,
                        self.drift_polydegree_combobox.currentIndex(),
                        self.tension_slider.value(),
                        self.drift_cont_startendcombobox.currentIndex(),
                        self.drift_plot_weighted.checkState(),
                        min_time,
                        max_time,
                        obstreeloop.name,
                    )
                    delta_model = self.populate_continuous_deltamodel(deltas)

                    if update:
                        self.plot_tares(self.axes_drift_cont_lower, obstreeloop)
                        self.plot_tares(self.axes_drift_cont_upper, obstreeloop)
                        ln = self.axes_drift_cont_lower.plot(xp, yp, 'k-')
                        if any(z):
                            textcolor = 'k'
                            if len(z) == 1:
                                if type(z[0]) is tuple:
                                    mean_drift, sigma = z[0][0], z[0][1]
                                    tstat = mean_drift / sigma
                                    if (
                                        np.abs(tstat) > 4.303
                                    ):  # Critical value for 95% CI, 2 DOF, 2-tailed t-test
                                        textcolor = 'r'
                                    z = [mean_drift]
                                annot_text = "{:.2f} µGal/hr".format(*z)
                            elif len(z) == 2:
                                annot_text = "{:.2f} µGal/hr per day".format(*z)
                            elif len(z) == 3:
                                annot_text = "{:.2f}*t^2 {:+.2f}*t {:+.2f}".format(*z)
                            elif len(z) == 4:
                                annot_text = "{:.2f}*t^3 {:+.2f}*t^2 {:+.2f}*t {:+.2f}".format(
                                    *z
                                )
                            else:
                                annot_text = ""
                            annot = self.axes_drift_cont_lower.annotate(
                                annot_text,
                                xy=(737287, 45),
                                xytext=(-20, 20),
                                textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                color=textcolor,
                            )
                            # arrowprops=dict(arrowstyle="->"))
                            annot.set_visible(False)

                            def update_annot(ind):
                                x, y = ln[0].get_data()
                                annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])

                            def hover(event):
                                vis = annot.get_visible()
                                if event.inaxes == self.axes_drift_cont_lower:
                                    cont, ind = ln[0].contains(event)
                                    if cont:
                                        update_annot(ind)
                                        annot.set_visible(True)
                                        # fig.canvas.draw_idle()
                                    else:
                                        if vis:
                                            annot.set_visible(False)
                                    self.drift_cont_figbot.canvas.draw_idle()
                                    # fig.canvas.draw_idle()

                            self.drift_cont_figbot.canvas.mpl_connect(
                                'motion_notify_event', hover
                            )
                        self.axes_drift_cont_lower.set_ylim(
                            np.round(min(drift_rate), 0) - 5,
                            np.round(max(drift_rate), 0) + 5,
                        )
                        self.axes_drift_cont_upper.set_title(
                            'Survey '
                            + obstreesurvey.name
                            + ', Loop '
                            + obstreeloop.name
                        )
                        self.drift_cont_canvasbot.draw()
                        self.drift_cont_canvastop.draw()
                    QtWidgets.QApplication.restoreOverrideCursor()
                    return delta_model
                except IndexError as e:
                    if self.drift_polydegree_combobox.currentIndex() == 1:
                        self.msg = show_message(
                            'Insufficient drift observations for spline method', 'Error'
                        )
                    self.msg = show_message('Unknown error', 'Unknown error')
                    self.drift_polydegree_combobox.setCurrentIndex(0)
                except np.linalg.LinAlgError as e:
                    logging.error(e)
                    self.msg = show_message(
                        'Insufficient drift observations for ' 'polynomial method',
                        'Error',
                    )
                    self.drift_polydegree_combobox.setCurrentIndex(0)
                    obstreeloop.drift_cont_method = 0
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
                    x = [
                        dt.datetime.utcfromtimestamp((f - 719163) * 86400.0)
                        for f in line[0]
                    ]
                    a = self.axes_drift_single.plot(x, y, '.-', picker=5)
                    a[0].name = line[2]

            for line in models[2]:
                if update:
                    self.axes_drift_single.plot(line[0], line[1], '--')
            if plot_data and update:
                self.axes_drift_single.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            if update:
                self.axes_drift_single.yaxis.set_label_text(
                    'Change in gravity since initial \nstation occupation, '
                    + 'in microGal'
                )
                self.drift_fig.canvas.mpl_connect(
                    'pick_event',
                    lambda event: self.show_line_label(event, self.axes_drift_single),
                )
                self.axes_drift_single.set_title(
                    'Survey ' + obstreesurvey.name + ', Loop ' + obstreeloop.name
                )
                self.drift_single_canvas.draw()
            QtWidgets.QApplication.restoreOverrideCursor()
            return models
        QtWidgets.QApplication.restoreOverrideCursor()

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
        Called from update_drift_tables_and_plots + callback from GUI. 
        Initiates plotting on drift tab.
        :param update: Boolean or int, controls if plots are updated. For 
        performance, it's set to false when loading a file
        """

        if type(update) is int:
            update = True
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_loop
        )
        method_key = self.driftmethod_combobox.currentIndex()

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
        # menu-item callback
        if update:
            width = self.drift_window.sizes()
            if method == 'none':
                self.drift_none()
            if method == 'netadj':
                self.drift_polydegree_combobox.setCurrentIndex(
                    obstreeloop.drift_netadj_method
                )
                self.drift_adjust()
            if method == 'roman':
                self.drift_roman()
            if method == 'continuous':
                self.drift_polydegree_combobox.setCurrentIndex(
                    obstreeloop.drift_cont_method
                )
                self.drift_cont_startendcombobox.setCurrentIndex(
                    obstreeloop.drift_cont_startend
                )
                self.drift_plot_weighted.setCheckState(obstreeloop.drift_cont_weighting)
                self.drift_continuous()
            else:
                self.disable_weighted_checkbox()
            self.set_width(width, method)

        model = self.plot_drift(update=update)
        if method == 'roman':
            obstreeloop.delta_model = model[1]
        else:
            obstreeloop.delta_model = model

        if method != orig_method or self.parent.workspace_loaded:
            # Leave the status bar hollow if a workspace was loaded
            if not self.parent.workspace_loaded:
                self.parent.deltas_update_required()
            self.parent.workspace_loaded = False

        if update:
            self.update_delta_model(method, model)

    def set_width(self, width, method):
        """
        Maintains relative width of plot windows when switching between drift-correction methods.
        """
        if method == 'none' or method == 'netadj' or method == 'roman':
            if width[0] > width[1]:
                self.drift_window.setSizes([width[0], width[1], width[2]])
            else:
                self.drift_window.setSizes([width[1], width[0], width[2]])
        else:
            if width[0] > width[1]:
                self.drift_window.setSizes([width[1], width[0], width[2]])
            else:
                self.drift_window.setSizes([width[0], width[1], width[2]])

    def update_delta_model(self, method, model):
        """
        Show appropriate delta model for the selected loop.

        This method takes the model generated by plot_drift() and assigns it to 
        the delta_view on the drift tab.
        :param method: If 'roman', show sample and average models. Otherwise, 
            show a single model.
        :param model: a PyQt model or list of models (Roman method).
        """
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_loop
        )
        if method == 'roman':
            # Hide drift correction, std_for_adj, and residual columns
            self.dg_samples_view.setModel(model[0])
            self.delta_view.setModel(model[1])
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
            obstreeloop = self.parent.obsTreeModel.itemFromIndex(
                self.parent.index_current_loop
            )
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
        self.drift_plot_weighted.setEnabled(False)
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
        self.drift_plot_weighted.setEnabled(False)
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
        if self.drift_polydegree_combobox.currentIndex() == 0:
            self.enable_weighted_checkbox()
        else:
            self.disable_weighted_checkbox()

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
        self.drift_plot_weighted.setEnabled(False)
        self.dg_samples_view.hide()

    def disable_weighted_checkbox(self):
        self.drift_plot_weighted.setEnabled(False)
        self.drift_plot_weighted.setToolTip(
            'Weighted observations is only enabled when Continuous '
            'model drift correction method and Constant drift model type are selected.'
        )

    def enable_weighted_checkbox(self):
        self.drift_plot_weighted.setEnabled(True)
        self.drift_plot_weighted.setToolTip('')

    def drift_combobox_updated(self):
        """
        Called when either the drift poly degree or extrapolate/constant 
        combobox is changed.
        """

        method_key = self.drift_polydegree_combobox.currentIndex()
        startend_key = self.drift_cont_startendcombobox.currentIndex()
        obstreeloop = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_loop
        )
        drift_method = obstreeloop.drift_method
        self.parent.deltas_update_required()

        if drift_method == 'continuous':
            obstreeloop.drift_cont_method = method_key
            obstreeloop.drift_cont_startend = startend_key
            if method_key == 1:
                self.tension_slider.setEnabled(True)
            else:
                self.tension_slider.setEnabled(False)
            if method_key == 0:
                self.enable_weighted_checkbox()
            else:
                self.disable_weighted_checkbox()

        elif drift_method == 'netadj':
            obstreeloop.drift_netadj_method = method_key

        model = self.plot_drift()
        self.update_delta_model(drift_method, model)
