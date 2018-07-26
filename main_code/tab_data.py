#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tab_data.py
===============

PyQt graphical elements on the data tab of GSadjust.
--------------------------------------------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtWidgets, QtCore
import numpy as np
from matplotlib.dates import DateFormatter
import matplotlib.pyplot as plt
from gui_objects import show_message


###########################################################################
# GSadjust data selection tab
###########################################################################
# noinspection PyUnresolvedReferences
class TabData(QtWidgets.QWidget):
    """
    station data selection

    On the left panel is a treeview of the data.
    TreeWidgetItems are defined in pyqt_models.py.

    On the middle panel is a tableview of the selected station (within
    a survey/loop, clicked on the left panel). Table items are also defined
    in pyqt_models.py.

    On the right panel is a matplotlib-type figure of relevant data.

    Some options for automatic selection of data are also available, either
    to select data among the whole data set, or among the current displayed
    table.

    """

    dpi = 100


    def __init__(self, parent):
        super(TabData, self).__init__()
        self.parent = parent

        # Set up data tree view
        # right panel: plot and some options
        main_frame = QtWidgets.QWidget()
        fig = Figure((3.0, 2.0), dpi=self.dpi)
        # make some room for the axes labels (set more place between subplots)
        fig.subplots_adjust(right=0.95, wspace=0.3, hspace=0.35)

        self.data_canvas = FigureCanvas(fig)
        self.data_canvas.setParent(main_frame)
        self.axes_data_UL = fig.add_subplot(221)
        self.axes_data_UR = fig.add_subplot(222)
        self.axes_data_LR = fig.add_subplot(224)
        self.axes_data_LL = fig.add_subplot(223)

        self.data_mpl_toolbar = NavigationToolbar(self.data_canvas, main_frame)

        # Plot panel: define some buttons with actions (signal/slot events) and line edits
        checkselected_button = QtWidgets.QPushButton("&check selected", self)
        checkselected_button.clicked.connect(self.checkselected)
        uncheckselected_button = QtWidgets.QPushButton("&uncheck selected", self)
        uncheckselected_button.clicked.connect(self.uncheckselected)
        checkall_button = QtWidgets.QPushButton("&check all", self)
        checkall_button.clicked.connect(self.checkall)
        uncheckall_button = QtWidgets.QPushButton("&uncheck all", self)
        uncheckall_button.clicked.connect(self.uncheckall)
        autoselec_tilts = QtWidgets.QWidget()
        autoselec_sd = QtWidgets.QWidget()
        autoselec_grav = QtWidgets.QWidget()
        autoselec_dur = QtWidgets.QWidget()
        autoselec_all = QtWidgets.QWidget()
        autoselec_tilts.button = QtWidgets.QPushButton("&auto uncheck tilts >", self)
        autoselec_tilts.button.clicked.connect(lambda: self.autoselect_tilt(autoselec_tilts))
        autoselec_sd.button = QtWidgets.QPushButton("&auto uncheck SD >", self)
        autoselec_sd.button.clicked.connect(lambda: self.autoselect_sd(autoselec_sd))
        autoselec_grav.button = QtWidgets.QPushButton("&auto uncheck g >", self)
        autoselec_grav.button.clicked.connect(lambda: self.autoselect_grav(autoselec_grav))
        autoselec_dur.button = QtWidgets.QPushButton("&auto uncheck dur <>", self)
        autoselec_dur.button.clicked.connect(lambda: self.autoselect_dur(autoselec_dur))
        autoselec_all.button = QtWidgets.QPushButton("&apply to all data", self)
        autoselec_all.button.clicked.connect(lambda: self.autoselect_all(autoselec_tilts.val,
                                                                         autoselec_sd.val, autoselec_grav.val,
                                                                         autoselec_dur.val))
        autoselec_tilts.val = QtWidgets.QLineEdit()
        autoselec_sd.val = QtWidgets.QLineEdit()
        autoselec_grav.val = QtWidgets.QLineEdit()
        autoselec_dur.val = QtWidgets.QLineEdit()

        # for the right panel (options & display):
        # create sublayouts (allow the line edit to be of finite extent)
        layout_options = QtWidgets.QGridLayout()
        grid = QtWidgets.QGridLayout()

        # fill subplayouts:
        grid.addWidget(self.data_canvas, 3, 0, 19, 4)
        grid.addWidget(self.data_mpl_toolbar, 2, 0, 1, 4)
        layout_options.addWidget(checkselected_button, 0, 0)
        layout_options.addWidget(uncheckselected_button, 0, 1)
        layout_options.addWidget(autoselec_tilts.button, 0, 2)
        layout_options.addWidget(autoselec_tilts.val, 0, 3, 1, 1)
        layout_options.addWidget(autoselec_grav.button, 0, 4)
        layout_options.addWidget(autoselec_grav.val, 0, 5, 1, 1)
        layout_options.addWidget(autoselec_dur.button, 1, 4)
        layout_options.addWidget(autoselec_dur.val, 1, 5, 1, 1)
        layout_options.addWidget(checkall_button, 1, 0)
        layout_options.addWidget(uncheckall_button, 1, 1)
        layout_options.addWidget(autoselec_sd.button, 1, 2)
        layout_options.addWidget(autoselec_sd.val, 1, 3, 1, 1)
        # layout_options.addWidget(updateplot_button, 0, 6)
        layout_options.addWidget(autoselec_all.button, 1, 6)
        # add subplayouts to the main layout
        grid.addLayout(layout_options, 0, 0, 1, 2)

        self.data_view = QtWidgets.QTableView()
        self.data_view.setTabKeyNavigation(False)

        splitter_right = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter_right.addWidget(self.data_view)

        tmp = QtWidgets.QWidget()
        tmp.setLayout(grid)
        splitter_right.addWidget(tmp)
        window_right = QtWidgets.QHBoxLayout()
        window_right.addWidget(splitter_right)
        self.setLayout(window_right)

    def update_station_plot(self, station, meter_type):
        """
        Updates relative-gravity data plots shown on 'Data' tab by creating a new PyQt model array for the specified
        station.
        """
        # Store new station key
        self.parent.station_model.layoutAboutToBeChanged.emit()
        self.parent.station_model.createArrayData(station)
        self.parent.station_model.layoutChanged.emit()

        t = station.t
        keepdata = station.keepdata
        t_selec = [t[i] for i in range(len(t)) if keepdata[i] == 1]
        # gravity channel (convert to microgals for display)
        series = np.array(station.grav)
        series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
        if meter_type == 'Scintrex' or meter_type == 'CG6':
            self.set_plot(self.axes_data_UL, t, series, t_selec, series_selec, 'gravity', '$\mu$gal', '1')
            # tiltx channel
            series = np.array(station.tiltx)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_UR, t, series, t_selec, series_selec, 'tilt X', 'arcsec', '2')
            # tilty channel
            series = np.array(station.tilty)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_LR, t, series, t_selec, series_selec, 'tilt Y', 'arcsec', '3')
            # SD channel
            series = np.array(station.sd)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_LL, t, series, t_selec, series_selec, 'standard deviation', '$\mu$gal', '4')
        elif meter_type == 'Burris':
            self.set_plot(self.axes_data_UL, t, series, t_selec, series_selec, 'gravity', '$\mu$gal', '1')
            # tiltx channel
            series = np.array(station.feedback)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_UR, t, series, t_selec, series_selec, 'Feedback', 'mV', '2')
            # tilty channel
            series = np.array(station.etc)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_LR, t, series, t_selec, series_selec, 'Earth TIde Correction', 'microGal', '3')
            # SD channel
            series = np.array(station.tiltx)
            series_selec = [series[i] for i in range(len(series)) if keepdata[i] == 1]
            self.set_plot(self.axes_data_LL, t, series, t_selec, series_selec, 'Tilt Correction', 'microGal', '4')
            self.data_canvas.draw()
        return True

    def clear_axes(self):
        self.axes_data_LL.clear()
        self.axes_data_LR.clear()
        self.axes_data_UL.clear()
        self.axes_data_UR.clear()
        self.data_canvas.draw()

    def set_plot(self, axe, seriex, seriey, seriex_selec, seriey_selec, serie_type, serie_unit, plot_location):
        """
        plot data samples from relative-gravity meter at a single station (gravity, tilt, temp., etc.).
        :param axe: Axes on which to plot
        :param seriex: X data, all points
        :param seriey: Y data, all points
        :param seriex_selec: X data, highlighted (blue) points
        :param seriey_selec: Y data, highlighted (blue) points
        :param serie_type: Name of data to plot
        :param serie_unit: Data units for y-axis label
        :return:
        """
        axe.cla()
        axe.grid(True)
        xfmt = DateFormatter('%H:%M')
        # It should be possible to just set_data on the plot object returned from axe.plot, but I couldn't figure it out
        if serie_type == 'gravity' and seriey_selec:  # Plot horizontal line at mean g
            mean_g = np.mean(seriey_selec)
            axe.plot([seriex[0], seriex[-1]], [mean_g, mean_g], 'o-', color='b', label=serie_type)

        setattr(self, 'plot1_' + plot_location, axe.plot(seriex, seriey, 'o-', color='k', label=serie_type))
        setattr(self, 'plot2_' + plot_location, axe.plot(seriex_selec, seriey_selec, 'o-', color='b', label=serie_type))
        axe.set_ylabel(serie_unit, size='x-small')
        axe.set_title(serie_type, size='x-small')
        labels = axe.get_xticklabels() + axe.get_yticklabels()
        for label in labels:
            label.set_size('x-small')
        xfmt = DateFormatter('%H:%M')
        axe.xaxis.set_major_formatter(xfmt)
        plt.setp(axe.get_xticklabels(), rotation=30, horizontalalignment='right')
        # Scale y axis to mean +/- 10. If a larger range is needed, color the axis red as a warning.
        if serie_type == 'gravity':
            g_range_pos = max(seriey) - mean_g
            g_range_neg = mean_g - min(seriey)
            if g_range_neg <= 10 and g_range_pos <= 10:
                axe.set_ylim(mean_g-10, mean_g+10)
                axe.yaxis.label.set_color('black')
                axe.spines['left'].set_color('black')
                axe.tick_params(axis='y', colors='black')
            else:
                axe.yaxis.label.set_color('red')
                axe.spines['left'].set_color('red')
                axe.tick_params(axis='y', colors='red')


    # All of these check/uncheck routines are very slow. Why?
    def autoselect_tilt(self, autoselec):
        """
        function for automatic selection of data based on simple thresholds:
        Tilts: absolute value higher than threshold are set to keepdata=0
        """
        if obstreestation.meter_type == 'Burris':
            self.msg = show_message('Not implemented for Burris data', 'Data selection error')
            return
        tilt_column1 = 3
        tilt_column2 = 4
        tilt_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx1 = self.parent.station_model.index(i, tilt_column1)
            indx2 = self.parent.station_model.index(i, tilt_column2)
            datax = float(self.parent.station_model.data(indx1, role=QtCore.Qt.DisplayRole))
            datay = float(self.parent.station_model.data(indx2, role=QtCore.Qt.DisplayRole))
            if abs(datax) > tilt_threshold or abs(datay) > tilt_threshold:
                idx_chk = indx1.sibling(i, 0)
                self.parent.station_model.setData(idx_chk, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_sd(self, autoselec):
        """
        function for automatic selection of data based on simple thresholds
        sd: sd values higher than threshold are set to keepdata=0
        """
        if obstreestation.meter_type == 'Burris':
            self.msg = show_message('Not implemented for Burris data', 'Data selection error')
            return
        sd_column = 2
        sd_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx = self.parent.station_model.index(i, sd_column)
            data = float(self.parent.station_model.data(indx, role=QtCore.Qt.DisplayRole))
            if abs(data) > sd_threshold:
                idx_chk = indx.sibling(i, 0)
                self.parent.station_model.setData(idx_chk, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_dur(self, autoselec):
        """
        function for automatic selection of data based on simple thresholds
        dur: duration different than given value are set to keepdata=0
        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(self.parent.currentStationIndex)
        if obstreestation.meter_type == 'Burris':
            self.msg = show_message('Not implemented for Burris data', 'Data selection error')
            return
        dur_column = 6
        dur_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx = self.parent.station_model.index(i, dur_column)
            data = float(self.parent.station_model.data(indx, role=QtCore.Qt.DisplayRole))
            if abs(data) > dur_threshold:
                idx_chk = indx.sibling(i, 0)
                self.parent.station_model.setData(idx_chk, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_grav(self, autoselec):
        """
        function for automatic selection of data based on simple thresholds
        grav: absolute values higher than threshold offset with respect to the
        mean value from 3 last points are set to keepdata=0
        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(self.parent.currentStationIndex)
        if obstreestation.meter_type == 'Scintrex':
            g_column = 1
        elif obstreestation.meter_type == 'Burris':
            g_column = 4
        g_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        g = obstreestation.grav
        stabilized_value = np.mean(g[len(g) - 3:len(g)])
        for i in range(self.parent.station_model.rowCount()):
            indx = self.parent.station_model.index(i, g_column)
            data = float(self.parent.station_model.data(indx, role=QtCore.Qt.DisplayRole))
            if abs(data - stabilized_value) > g_threshold:
                idx_chk = indx.sibling(i, 0)
                self.parent.station_model.setData(idx_chk, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_all(self, tilts_thrshld, sd_thrshld, grav_thrshld, dur_thrshld):
        """
        function for automatic selection of data based on simple thresholds
        apply all selection critera which have been input by the user to all
        the data set.
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        selec_grav, selec_sd, selec_tilts, selec_dur = False, False, False, False
        g_threshold, sd_threshold, tilt_threshold, dur_threshold = 0, 0, 0, 0
        if grav_thrshld.text():
            g_threshold = float(grav_thrshld.text())
            selec_grav = True
        if sd_thrshld.text():
            sd_threshold = float(sd_thrshld.text())
            selec_sd = True
        if tilts_thrshld.text():
            tilt_threshold = float(tilts_thrshld.text())
            selec_tilts = True
        if dur_thrshld.text():
            dur_threshold = float(dur_thrshld.text())
            selec_dur = True

        for i in range(self.parent.obsTreeModel.invisibleRootItem().rowCount()):
            obstreesurvey = self.parent.obsTreeModel.invisibleRootItem().child(i)
            for ii in range(obstreesurvey.rowCount()):
                obstreeloop = obstreesurvey.child(ii)
                for iii in range(obstreeloop.rowCount()):
                    obstreestation = obstreeloop.child(iii)
                    g = obstreestation.grav
                    stabilized_value = np.mean(g[len(g) - 3:len(g)])
                    for iiii in range(len(obstreestation.keepdata)):
                        indx = self.parent.station_model.index(iiii, 0)
                        if selec_grav and abs(g[iiii]) > g_threshold + stabilized_value:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(indx, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
                            continue
                        if selec_sd and obstreestation.sd[iiii] > sd_threshold:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(indx, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
                            continue
                        if selec_dur and obstreestation.dur[iiii] > dur_threshold:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(indx, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
                            continue
                        if selec_tilts:
                            if obstreestation.meter_type == 'Scintrex':
                                if abs(obstreestation.tiltx[iiii]) > tilt_threshold and \
                                   abs(obstreestation.tilty[iiii]) > tilt_threshold:
                                    obstreestation.keepdata[iiii] = 0
                                    self.parent.station_model.setData(indx, QtCore.Qt.Unchecked,
                                                                      QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def uncheckall(self):
        """
        uncheck all items in the displayed table when button is clicked
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx = self.parent.station_model.index(i, 0)
            self.parent.station_model.setData(indx, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def checkall(self):
        """
        check all items in the displayed table when button is clicked
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx = self.parent.station_model.index(i, 0)
            self.parent.station_model.setData(indx, QtCore.Qt.Checked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def checkselected(self):
        """
        check all selected items (select the station column)
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        qmodelindex = self.data_view.selectedIndexes()
        for indx in list(qmodelindex):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx2 = indx.sibling(indx.row(), 0)
            self.parent.station_model.setData(indx2, QtCore.Qt.Checked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()

    def uncheckselected(self):
        """
        check all selected items (select the station column)
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        qmodelindex = self.data_view.selectedIndexes()
        for indx in list(qmodelindex):
            # because the setData function currently only works for column 0,
            # then pass the index of the element with same row as selected (if
            # the user does not select element from first column) and column=0
            indx2 = indx.sibling(indx.row(), 0)
            self.parent.station_model.setData(indx2, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()
