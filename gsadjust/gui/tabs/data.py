"""
gui/tabs/data.py
================

PyQt graphical elements on the data tab of GSadjust.
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

import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.dates import DateFormatter
from matplotlib.figure import Figure

from ..messages import MessageBox


###########################################################################
# GSadjust data selection tab
###########################################################################
# noinspection PyUnresolvedReferences
class TabData(QtWidgets.QWidget):
    """
    Tab for station data selection

    On the left panel is a tableview of the selected station (within
    a survey/loop, clicked on the left-hand tree view). Table items are also defined
    in models.

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

        # Plot panel: define some buttons with actions (signal/slot events) and line
        # edits
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
        autoselec_alldata = QtWidgets.QWidget()
        autoselec_tilts.button = QtWidgets.QPushButton("&auto uncheck |tilt| >", self)
        autoselec_tilts.button.clicked.connect(
            lambda: self.autoselect_tilt(autoselec_tilts)
        )
        autoselec_sd.button = QtWidgets.QPushButton("&auto uncheck SD >", self)
        autoselec_sd.button.clicked.connect(lambda: self.autoselect_sd(autoselec_sd))
        autoselec_grav.button = QtWidgets.QPushButton("&auto uncheck g >", self)
        autoselec_grav.button.clicked.connect(
            lambda: self.autoselect_grav(autoselec_grav)
        )
        autoselec_dur.button = QtWidgets.QPushButton("&auto uncheck dur >", self)
        autoselec_dur.button.clicked.connect(lambda: self.autoselect_dur(autoselec_dur))
        autoselec_all.button = QtWidgets.QPushButton("&apply all filters", self)
        autoselec_all.button.clicked.connect(
            lambda: self.autoselect_all(
                autoselec_tilts.val,
                autoselec_sd.val,
                autoselec_grav.val,
                autoselec_dur.val,
            )
        )
        autoselec_alldata.button = QtWidgets.QPushButton(
            "&apply all filters to all data", self
        )
        autoselec_alldata.button.clicked.connect(
            lambda: self.autoselect_alldata(
                autoselec_tilts.val,
                autoselec_sd.val,
                autoselec_grav.val,
                autoselec_dur.val,
            )
        )
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
        layout_options.addWidget(autoselec_all.button, 0, 6)
        layout_options.addWidget(autoselec_dur.button, 1, 4)
        layout_options.addWidget(autoselec_dur.val, 1, 5, 1, 1)
        layout_options.addWidget(checkall_button, 1, 0)
        layout_options.addWidget(uncheckall_button, 1, 1)
        layout_options.addWidget(autoselec_sd.button, 1, 2)
        layout_options.addWidget(autoselec_sd.val, 1, 3, 1, 1)
        # layout_options.addWidget(updateplot_button, 0, 6)
        layout_options.addWidget(autoselec_alldata.button, 1, 6)
        # add subplayouts to the main layout
        grid.addLayout(layout_options, 0, 0, 1, 2)

        self.data_view = QtWidgets.QTableView()
        self.data_view.setTabKeyNavigation(False)

        splitter_right = QtWidgets.QSplitter(Qt.Horizontal)
        splitter_right.addWidget(self.data_view)

        tmp = QtWidgets.QWidget()
        tmp.setLayout(grid)
        splitter_right.addWidget(tmp)
        window_right = QtWidgets.QHBoxLayout()
        window_right.addWidget(splitter_right)
        self.setLayout(window_right)

    def update_station_plot(self, station, meter_type):
        """
        Updates relative-gravity data plots shown on 'Data' tab by creating a
        new PyQt model array for the specified station.

        Parameters
        ----------
        station : ObsTreeStation
        meter_type : {'csv', 'CG5', 'CG6', 'CG6TSoft', 'Burris'}

        Returns
        -------

        """
        self.parent.station_model.layoutAboutToBeChanged.emit()
        self.parent.station_model.createArrayData(station)
        self.parent.station_model.layoutChanged.emit()

        t = station.t
        keepdata = station.keepdata
        t_selec = [t[i] for i in range(len(t)) if keepdata[i] == 1]
        # gravity channel (convert to microgals for display)
        series = np.array(station.grav())
        series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
        if meter_type == "CG5" or meter_type == "CG6" or meter_type == "csv":
            self.set_plot(
                self.axes_data_UL,
                t,
                series,
                t_selec,
                series_selec,
                "Gravity",
                r"µgal",
            )
            # tiltx channel
            series = np.array(station.tiltx)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_UR,
                t,
                series,
                t_selec,
                series_selec,
                "Tilt X",
                "arcsec",
            )
            # tilty channel
            series = np.array(station.tilty)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LR,
                t,
                series,
                t_selec,
                series_selec,
                "Tilt Y",
                "arcsec",
            )
            # SD channel
            series = np.array(station.sd)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LL,
                t,
                series,
                t_selec,
                series_selec,
                "Standard deviation",
                r"µgal",
            )

        elif meter_type == "CG6Tsoft":
            self.set_plot(
                self.axes_data_UL,
                t,
                series,
                t_selec,
                series_selec,
                "Gravity",
                r"µgal",
            )
            # tiltx channel
            series = np.array(station.tiltx)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_UR,
                t,
                series,
                t_selec,
                series_selec,
                "Tilt X",
                "arcsec",
            )
            # tilty channel
            series = np.array(station.tilty)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LR,
                t,
                series,
                t_selec,
                series_selec,
                "Tilt Y",
                "arcsec",
            )
            # SD channel
            series = np.array(station.etc)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LL,
                t,
                series,
                t_selec,
                series_selec,
                "Earth Tide Correction",
                r"µgal",
            )

        elif meter_type == "Burris":
            self.set_plot(
                self.axes_data_UL,
                t,
                series,
                t_selec,
                series_selec,
                "Gravity",
                r"µgal",
            )
            # tiltx channel
            series = np.array(station.feedback)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_UR,
                t,
                series,
                t_selec,
                series_selec,
                "Feedback",
                "mV",
            )
            # tilty channel
            series = np.array(station.etc)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LR,
                t,
                series,
                t_selec,
                series_selec,
                "Earth Tide Correction",
                r"µgal",
            )
            # SD channel
            series = np.array(station.tiltx)
            series_selec = [s for i, s in enumerate(series) if keepdata[i] == 1]
            self.set_plot(
                self.axes_data_LL,
                t,
                series,
                t_selec,
                series_selec,
                "Tilt Correction",
                r"µgal",
            )
        self.data_canvas.draw()
        return True

    def clear_axes(self):
        """
        Called when clearing a workspace.

        """
        self.axes_data_LL.clear()
        self.axes_data_LR.clear()
        self.axes_data_UL.clear()
        self.axes_data_UR.clear()
        self.data_canvas.draw()

    def set_plot(
        self, axe, seriex, seriey, seriex_selec, seriey_selec, serie_type, serie_unit
    ):
        """
        Plot data samples from relative-gravity meter at a single station
        (gravity, tilt, temp., etc.).

        Parameters
        ----------
        axe : AxesSubplot
            Axes on which to plot
        seriex : list
            List of data times (float). Plots in black
        seriey : ndarray
            Y data to plot. Plots in black
        seriex_selec : list
            List of selected data times. Plots in blue
        seriey_selec : list
            List of selected data. Plots in blue
        serie_type : str
            Data type, varies based on meter: 'Gravity', 'Tilt X', etc.
        serie_unit: str
            Data units, e.g. µGal

        """
        axe.cla()
        axe.grid(True)
        xfmt = DateFormatter("%H:%M")
        mean_g = np.mean(seriey_selec)
        if serie_type == "Gravity" and seriey_selec:  # Plot horizontal line at mean g
            axe.plot(
                [seriex[0], seriex[-1]],
                [mean_g, mean_g],
                "-",
                color="b",
                label=serie_type,
            )
        axe.plot(seriex, seriey, "o-", color="k", label=serie_type)
        axe.plot(seriex_selec, seriey_selec, "o-", color="b", label=serie_type)
        # axe.set_ylabel(serie_unit, size="x-small")
        # axe.set_title(serie_type, size="x-small")
        labels = axe.get_xticklabels() + axe.get_yticklabels()
        # for label in labels:
        #     label.set_size("x-small")
        axe.xaxis.set_major_formatter(xfmt)
        axe.set_xticklabels(
            axe.get_xticklabels(), rotation=30, horizontalalignment="right"
        )
        # Scale y axis to mean +/- 10. If a larger range is needed, color the axis
        # red as a warning.
        if serie_type == "Gravity":
            g_range_pos = max(seriey) - mean_g
            g_range_neg = mean_g - min(seriey)
            if g_range_neg <= 10 and g_range_pos <= 10:
                axe.set_ylim(mean_g - 10, mean_g + 10)
                axe.yaxis.label.set_color("black")
                axe.spines["left"].set_color("black")
                axe.tick_params(axis="y", colors="black")
            else:
                axe.yaxis.label.set_color("red")
                axe.spines["left"].set_color("red")
                axe.tick_params(axis="y", colors="red")

    # All of these check/uncheck routines are very slow. Why?
    def autoselect_tilt(self, autoselec):
        """
        function for automatic selection of data based on tilt deviation.
        Scintrex only.

        Parameters
        ----------
        autoselec : QWidget
            Tilt values higher than threshold are set to keepdata=0 (unchecked)
        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_station
        )
        if obstreestation.meter_type == "Burris":
            MessageBox.warning(
                "Data selection error",
                "Not implemented for Burris data",
            )
            return
        tilt_column1 = 4
        tilt_column2 = 5
        tilt_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx1 = self.parent.station_model.index(i, tilt_column1)
            indx2 = self.parent.station_model.index(i, tilt_column2)
            datax = float(self.parent.station_model.data(indx1, role=Qt.DisplayRole))
            datay = float(self.parent.station_model.data(indx2, role=Qt.DisplayRole))
            if abs(datax) > tilt_threshold or abs(datay) > tilt_threshold:
                idx_chk = indx1.sibling(i, 0)
                self.parent.station_model.setData(
                    idx_chk, Qt.Unchecked, Qt.CheckStateRole
                )
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_sd(self, autoselec):
        """
        Function for automatic selection of data based on standard deviation threshold
        Scintrex only.

        Parameters
        ----------
        autoselec : QWidget
            Standard deviations higher than this value are set to keepdata=0 (unchecked)
        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_station
        )
        if obstreestation.meter_type == "Burris":
            MessageBox.warning(
                "Data selection error",
                "Not implemented for Burris data",
            )
            return
        sd_column = 3
        sd_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx = self.parent.station_model.index(i, sd_column)
            data = float(self.parent.station_model.data(indx, role=Qt.DisplayRole))
            if abs(data) > sd_threshold:
                idx_chk = indx.sibling(i, 0)
                self.parent.station_model.setData(
                    idx_chk, Qt.Unchecked, Qt.CheckStateRole
                )
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_dur(self, autoselec):
        """
        function for automatic selection of data based on duration

        Parameters
        ----------
        autoselec : QWidget
            Duration values higher than threshold are set to keepdata=0 (unchecked)
        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_station
        )
        if obstreestation.meter_type == "Burris":
            MessageBox.warning(
                "Data selection error",
                "Not implemented for Burris data",
            )
            return
        dur_column = 7
        dur_threshold = float(autoselec.val.text())
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        for i in range(self.parent.station_model.rowCount()):
            indx = self.parent.station_model.index(i, dur_column)
            data = float(self.parent.station_model.data(indx, role=Qt.DisplayRole))
            if abs(data) > dur_threshold:
                idx_chk = indx.sibling(i, 0)
                self.parent.station_model.setData(
                    idx_chk, Qt.Unchecked, Qt.CheckStateRole
                )
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_grav(self, autoselec):
        """
        function for automatic selection of data based on deviation from the mean

        Parameters
        ----------
        autoselec : QWidget
            absolute values higher than threshold offset with respect to the mean
            value from 3 last points are unchecked

        """
        obstreestation = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_station
        )
        if obstreestation.meter_type == "Burris":
            g_column = 4
        else:
            g_column = 2
        try:
            g_threshold = float(autoselec.val.text())
            QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
            g = obstreestation.grav()
            stabilized_value = obstreestation.gmean()
            for i in range(self.parent.station_model.rowCount()):
                indx = self.parent.station_model.index(i, g_column)
                data = float(self.parent.station_model.data(indx, role=Qt.DisplayRole))
                if abs(data - stabilized_value) > g_threshold:
                    idx_chk = indx.sibling(i, 0)
                    self.parent.station_model.setData(
                        idx_chk, Qt.Unchecked, Qt.CheckStateRole
                    )
            QtWidgets.QApplication.restoreOverrideCursor()
        except ValueError:  # If box is empty
            return

    def autoselect_all(self, tilts_thrshld, sd_thrshld, grav_thrshld, dur_thrshld):
        """
        Function for automatic selection of data based on simple thresholds.

        Apply all selection criteria which have been input by the user to all
        of the samples of the current station.

        Parameters
        ----------
        tilts_thrshld : QWidget
        sd_thrshld : QWidget
        grav_thrshld : QWidget
        dur_thrshld : QWidget

        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
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

        obstreestation = self.parent.obsTreeModel.itemFromIndex(
            self.parent.index_current_station
        )
        g = obstreestation.grav()
        stabilized_value = obstreestation.gmean()
        for iiii in range(len(obstreestation.keepdata)):
            indx = self.parent.station_model.index(iiii, 0)
            if selec_grav and abs(g[iiii] - stabilized_value) > g_threshold:
                obstreestation.keepdata[iiii] = 0
                self.parent.station_model.setData(indx, Qt.Unchecked, Qt.CheckStateRole)
                continue
            if selec_sd and obstreestation.sd[iiii] > sd_threshold:
                obstreestation.keepdata[iiii] = 0
                self.parent.station_model.setData(indx, Qt.Unchecked, Qt.CheckStateRole)
                continue
            if selec_dur and obstreestation.dur[iiii] > dur_threshold:
                obstreestation.keepdata[iiii] = 0
                self.parent.station_model.setData(indx, Qt.Unchecked, Qt.CheckStateRole)
                continue
            if selec_tilts:
                if (
                    obstreestation.meter_type == "CG5"
                    or obstreestation.meter_type == "CG6"
                    or obstreestation.meter_type == "CG6Tsoft"
                ):
                    if (
                        abs(obstreestation.tiltx[iiii]) > tilt_threshold
                        and abs(obstreestation.tilty[iiii]) > tilt_threshold
                    ):
                        obstreestation.keepdata[iiii] = 0
                        self.parent.station_model.setData(
                            indx, Qt.Unchecked, Qt.CheckStateRole
                        )
        QtWidgets.QApplication.restoreOverrideCursor()

    def autoselect_alldata(self, tilts_thrshld, sd_thrshld, grav_thrshld, dur_thrshld):
        """
        Function for automatic selection of data based on simple thresholds.

        Apply all selection criteria which have been input by the user to all
        of the data in the entire campaign.

        Pretty slow!

        Parameters
        ----------
        tilts_thrshld : QWidget
        sd_thrshld : QWidget
        grav_thrshld : QWidget
        dur_thrshld : QWidget

        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
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
                    g = obstreestation.grav()
                    stabilized_value = obstreestation.gmean()
                    for iiii in range(len(obstreestation.keepdata)):
                        indx = self.parent.station_model.index(iiii, 0)
                        if selec_grav and abs(g[iiii] - stabilized_value) > g_threshold:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(
                                indx, Qt.Unchecked, Qt.CheckStateRole
                            )
                            continue
                        if selec_sd and obstreestation.sd[iiii] > sd_threshold:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(
                                indx, Qt.Unchecked, Qt.CheckStateRole
                            )
                            continue
                        if selec_dur and obstreestation.dur[iiii] > dur_threshold:
                            obstreestation.keepdata[iiii] = 0
                            self.parent.station_model.setData(
                                indx, Qt.Unchecked, Qt.CheckStateRole
                            )
                            continue
                        if selec_tilts:
                            if (
                                obstreestation.meter_type == "CG5"
                                or obstreestation.meter_type == "CG6"
                                or obstreestation.meter_type == "CG6Tsoft"
                            ):
                                if (
                                    abs(obstreestation.tiltx[iiii]) > tilt_threshold
                                    and abs(obstreestation.tilty[iiii]) > tilt_threshold
                                ):
                                    obstreestation.keepdata[iiii] = 0
                                    self.parent.station_model.setData(
                                        indx,
                                        Qt.Unchecked,
                                        Qt.CheckStateRole,
                                    )
        QtWidgets.QApplication.restoreOverrideCursor()

    def uncheckall(self):
        """
        uncheck all items in the displayed table when button is clicked
        """
        self.parent.station_model.uncheckAll()

    def checkall(self):
        """
        check all items in the displayed table when button is clicked
        """
        self.parent.station_model.checkAll()

    def checkselected(self):
        """
        check all selected items (select the station column)
        """
        self.check_or_uncheck_selected(Qt.Checked)

    def uncheckselected(self):
        """
        check all selected items (select the station column)
        """
        self.check_or_uncheck_selected(Qt.Unchecked)

    def check_or_uncheck_selected(self, check_type):
        """
        Check/Uncheck selected samples in the table

        Parameters
        ----------
        check_type: 0 or 2 (unchecked or checked)

        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        qmodelindex = self.data_view.selectedIndexes()
        col0_indexes = [x for x in list(qmodelindex) if x.column() == 0]
        for indx in col0_indexes[:-1]:
            # only update on the last one ("silent=False"). Saves a lot of time
            self.parent.station_model.setData(
                indx, check_type, Qt.CheckStateRole, silent=True
            )
        indx = col0_indexes[-1]
        self.parent.station_model.setData(indx, check_type, Qt.CheckStateRole)
        QtWidgets.QApplication.restoreOverrideCursor()
