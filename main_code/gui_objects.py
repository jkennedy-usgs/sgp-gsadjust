#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gui_objects.py
==============

Miscellaneous GUI objects for GSadjust
--------------------------------------------------------------------------------------------------------------------

Major GUI objects (tabs, table views) are in the tab_... files. This module has primarily pop-up dialogs used to set
network adjustment settings, show gravity change over time, etc. Major dialogs are written as classes and
instantiated in GSadjust.py. Minor dialogs are written as functions and called directly.

This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
import os
import datetime as dt
import numpy as np
from PyQt5 import QtGui, QtCore, QtWidgets
import logging
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from data_objects import Datum, AdjustmentOptions
from pyqt_models import GravityChangeModel, DatumTableModel, CustomSortingModel
import a10


class ApplyTimeCorrection(QtWidgets.QDialog):
    def __init__(self):
        super(ApplyTimeCorrection, self).__init__()
        self.time_correction_type = False
        self.msg = QtWidgets.QMessageBox()
        self.msg.setText("Apply time correction to:")
        self.msg.addButton(QtWidgets.QPushButton('Current station(s)'), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('Current loop'), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('Current survey'), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('All data'), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('Cancel'), QtWidgets.QMessageBox.RejectRole)
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        if btn.text() == 'Cancel':
            self.reject()
        else:
            options = {'Current station(s)': 'station',
                       'Current loop': 'loop',
                       'Current survey': 'survey',
                       'All data': 'all'}
            self.time_correction_type = options[btn.text()]
            self.accept()


class Overwrite(QtWidgets.QDialog):
    def __init__(self):
        super(Overwrite, self).__init__()
        self.overwrite = False
        self.msg = QtWidgets.QMessageBox()
        self.msg.setText("Overwrite or append to existing data?")
        self.msg.addButton(QtWidgets.QPushButton('Overwrite'), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('Append'), QtWidgets.QMessageBox.YesRole)
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        if btn.text() == 'Overwrite':
            self.overwrite = True
        else:
            self.overwrite = False
        self.accept()


class MeterType(QtWidgets.QDialog):
    def __init__(self):
        super(MeterType, self).__init__()
        self.msg = QtWidgets.QMessageBox()
        self.msg.setText("Choose meter file to import")
        self.msg.addButton(QtWidgets.QPushButton(' CG-3, CG-5 '), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton(' CG-6 '), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton(' Burris '), QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton(' Cancel '), QtWidgets.QMessageBox.RejectRole)
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        meter = {' CG-3, CG-5 ': 'Scintrex',
                 ' CG-6 ': 'CG6',
                 ' Burris ': 'Burris',
                 ' Cancel ': 'Cancel'}
        self.meter_type = meter[btn.text()]
        if self.meter_type == 'Cancel':
            self.reject()
        else:
            self.accept()


class AdjustOptions(QtWidgets.QDialog):
    """
    Dialog to set network adjustment options.
    """
    def __init__(self, survey_str, options):
        super(AdjustOptions, self).__init__()
        self.setGeometry(50, 50, 350, 350)
        self.ao = options
        self.surveys_to_update = ''
        self.update_options = []
        self.drift_temp_chk = QtWidgets.QCheckBox('Model temperature drift, polynomial degree:')
        self.sigma_factor_chk = QtWidgets.QCheckBox('Std. dev. multiplier')
        self.sigma_add_chk = QtWidgets.QCheckBox('Add to std. dev.')
        self.sigma_min_chk = QtWidgets.QCheckBox('Minimum std. dev.')
        self.cal_coeff_chk = QtWidgets.QCheckBox('Include relative meter calibration coefficient')
        self.alpha_text = QtWidgets.QLabel('Significance level for global model test')
        self.woutfiles_chk = QtWidgets.QCheckBox('Write output files')
        self.drift_temp_edit = QtWidgets.QLineEdit(str(self.ao.model_temp_degree))
        self.sigma_factor_edit = QtWidgets.QLineEdit(str(self.ao.sigma_factor))
        self.sigma_add_edit = QtWidgets.QLineEdit(str(self.ao.sigma_add))
        self.sigma_min_edit = QtWidgets.QLineEdit(str(self.ao.sigma_min))
        self.alpha_edit = QtWidgets.QLineEdit(str(self.ao.alpha))
        self.init_ui(survey_str)

    def init_ui(self, survey_name):
        if survey_name is not None:
            self.drift_temp_chk.setChecked(self.ao.use_model_temp)
            self.sigma_factor_chk.setChecked(self.ao.use_sigma_factor)
            self.sigma_add_chk.setChecked(self.ao.use_sigma_add)
            self.sigma_min_chk.setChecked(self.ao.use_sigma_min)
            self.cal_coeff_chk.setChecked(self.ao.cal_coeff)
            self.woutfiles_chk.setChecked(self.ao.woutfiles)

            # create buttons and actions
            cancel_button = QtWidgets.QPushButton('Cancel')
            cancel_button.clicked.connect(self.close)
            current_button = QtWidgets.QPushButton('Apply to current survey (' + survey_name + ')')
            current_button.clicked.connect(self.apply_current)
            all_button = QtWidgets.QPushButton('Apply to all surveys')
            all_button.clicked.connect(self.apply_all)

            buttonBox = QtWidgets.QDialogButtonBox(QtCore.Qt.Horizontal)
            buttonBox.addButton(current_button, QtWidgets.QDialogButtonBox.ActionRole)
            buttonBox.addButton(all_button, QtWidgets.QDialogButtonBox.ActionRole)
            buttonBox.addButton(cancel_button, QtWidgets.QDialogButtonBox.ActionRole)

            grid = QtWidgets.QGridLayout()
            grid.addWidget(self.drift_temp_chk, 1, 0)
            grid.addWidget(self.drift_temp_edit, 1, 1)
            grid.addWidget(self.sigma_factor_chk, 2, 0)
            grid.addWidget(self.sigma_factor_edit, 2, 1)
            grid.addWidget(self.sigma_add_chk, 3, 0)
            grid.addWidget(self.sigma_add_edit, 3, 1)
            grid.addWidget(self.sigma_min_chk, 4, 0)
            grid.addWidget(self.sigma_min_edit, 4, 1)
            grid.addWidget(self.cal_coeff_chk, 5, 0)
            grid.addWidget(self.alpha_text, 6, 0)
            grid.addWidget(self.alpha_edit, 6, 1)
            grid.addWidget(self.woutfiles_chk, 7, 0)
            grid.addWidget(buttonBox, 8, 0)

            self.setLayout(grid)
            self.setWindowTitle('Network adjustment options')
            self.setWindowModality(QtCore.Qt.ApplicationModal)
        else:
            self.msg = show_message('Please load a survey first', 'Network adjustment options')

    def set_adjust_options(self):
        if self.drift_temp_chk.isChecked():
            self.ao.use_model_temp = True
            self.ao.model_temp_degree = int(self.ao.model_temp_degree.text())
        else:
            self.ao.use_model_temp = False
        if self.sigma_factor_chk.isChecked():
            self.ao.use_sigma_factor = True
            self.ao.sigma_factor = float(self.sigma_factor_edit.text())
        else:
            self.ao.use_sigma_factor = False
        if self.sigma_add_chk.isChecked():
            self.ao.use_sigma_add = True
            self.ao.sigma_add = float(self.sigma_add_edit.text())
        else:
            self.ao.use_sigma_add = False
        if self.sigma_min_chk.isChecked():
            self.ao.use_sigma_min = True
            self.ao.sigma_min = float(self.sigma_min_edit.text())
        else:
            self.ao.use_sigma_min = False
        if self.cal_coeff_chk.isChecked():
            self.ao.cal_coeff = True
        else:
            self.ao.cal_coeff = False
        self.ao.alpha = float(self.alpha_edit.text())
        if self.woutfiles_chk.isChecked():
            self.ao.woutfiles = True
        else:
            self.ao.woutfiles = False

    def apply_current(self):
        self.set_adjust_options()
        self.surveys_to_update = 'single'
        self.accept()

    def apply_all(self):
        self.set_adjust_options()
        self.surveys_to_update = 'all'
        self.accept()


class IncrMinuteTimeEdit(QtWidgets.QTimeEdit):
    """
    Provides a QTimeEdit with that increments at some minute interval.

    The PyQt Time edit doesn't allow for rolling-over the minutes to hours, apparently. I.e., the hours and minutes
    boxes are limited to 0-60. This implements that somewhat, but one can't decrement minutes below 00 in either the
    MinuteSection or HourSection. This manifests as when at a time, e.g., 3:00, stepBy isn't called when in the
    Minute Section because time is already at 0.
    """
    def __init__(self, time):
        super(IncrMinuteTimeEdit, self).__init__(time)
        self.step = 10

    def stepBy(self, steps):
        if self.currentSection() == self.HourSection:
            self.incrTime(steps)
        if self.currentSection() == self.MinuteSection:
            self.incrTime(steps)

    def incrTime(self, steps):
        hours = self.dateTime().time().hour()
        minutes = self.dateTime().time().minute() + steps * self.step
        if minutes < 0:
            self.setTime(QtCore.QTime(hours-1, 60 + minutes))
        if minutes < 60:
            self.setTime(QtCore.QTime(hours, minutes))
        else:
            self.setTime(QtCore.QTime(hours+1, 60 % minutes))


class ProgressBar(QtWidgets.QWidget):
    """
    define progress bar
    """

    def __init__(self, parent=None, total=20, textmess='Progress'):
        super(ProgressBar, self).__init__(parent)
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)
        main_layout = QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 1)
        self.setLayout(main_layout)
        self.setWindowTitle(textmess)


def copy_cells_to_clipboard(table):
    if len(table.selectedIndexes()) > 0:
        # sort select indexes into rows and columns
        previous = table.selectedIndexes()[0]
        columns = []
        rows = []
        for index in table.selectedIndexes():  # columns first , then rows
            if previous.row() != index.row():
                columns.append(rows)
                rows = []
            rows.append(index.data())
            previous = index
        columns.append(rows)

        # add rows and columns to clipboard
        clipboard = ""
        ncols = len(columns[0])
        nrows = len(columns)

        # add header to clipboard
        for c in range(ncols):
            clipboard += table.model().headerData(c, QtCore.Qt.Horizontal)
            clipboard += '\t'
        clipboard += '\n'

        for r in range(nrows):
            for c in range(ncols):
                clipboard += columns[r][c]
                if c != (ncols - 1):
                    clipboard += '\t'
            clipboard += '\n'

        # copy to the system clipboard
        sys_clip = QtWidgets.QApplication.clipboard()
        sys_clip.setText(clipboard)
    else:
        msg = show_message('No rows selected (Ctrl-a to select all)', 'Copy warning')


def DateMethodDialog():
    """
    Dialog to select data import method - all data into one survey of with dates contained in a second file.
    :return: String indicating survey import type ('choose' or 'single')
    """
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Question)

    msg.setText("Start/end date file not found.")
    msg.setWindowTitle("Specify start/end dates")
    msg.addButton(QtWidgets.QPushButton('Choose File'), QtWidgets.QMessageBox.YesRole)
    msg.addButton(QtWidgets.QPushButton('Single Survey'), QtWidgets.QMessageBox.NoRole)
    msg.addButton(QtWidgets.QPushButton('Cancel'), QtWidgets.QMessageBox.RejectRole)
    method = msg.exec_()
    if method == 0:
        return 'choose'
    elif method == 1:
        return 'single'
    else:
        return False


def about_dialog():
    msg1 = '<html>GSadjust, a product of the USGS Southwest Gravity Program<br>' + \
           '<a href ="http://go.usa.gov/xqBnQ">http://go.usa.gov/xqBnQ</a>' + \
           '<br><br><a href ="https://github.com/jkennedy-usgs/sgp-gsadjust">' + \
           'https://github.com/jkennedy-usgs/sgp-gsadjust</a>' + \
           '<br><a href="mailto:jkennedy@usgs.gov">jkennedy@usgs.gov</a>'
    _, ok = QtWidgets.QMessageBox.about(None, "GSadust", msg1)


def VerticalGradientDialog(default_interval):
        text, ok = QtWidgets.QInputDialog.getDouble(None, "Vertical-gradient interval",
                                                    "Interval, in cm:",
                                                    default_interval,
                                                    0, 200, 1)
        if ok:
            return float(text)
        else:
            return default_interval


class DatumComparisonFigure(QtWidgets.QDialog):
    def __init__(self, obsTreeModel):
        super(DatumComparisonFigure, self).__init__()
        datum_names = []
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            for ii in range(survey.datum_model.rowCount()):
                datum = survey.datum_model.data(survey.datum_model.index(ii, 0), role=QtCore.Qt.UserRole)
                datum_names.append(datum.station)

        unique_datum_names = list(set(datum_names))

        xdata_all, ydata_obs_all, ydata_adj_all = [], [], []
        for name in unique_datum_names:
            xdata, ydata_obs, ydata_adj = [], [], []
            for i in range(obsTreeModel.rowCount()):
                survey = obsTreeModel.invisibleRootItem().child(i)
                for ii in range(survey.datum_model.rowCount()):
                    datum = survey.datum_model.data(survey.datum_model.index(ii, 0), role=QtCore.Qt.UserRole)
                    if datum.station == name:
                        xdata.append(date2num(dt.datetime.strptime(survey.name, '%Y-%m-%d')))
                        ydata_obs.append(datum.g)
                        ydata_adj.append(datum.g + datum.residual)
            ydata_obs = [i - ydata_obs[0] for i in ydata_obs]
            ydata_adj = [i - ydata_adj[0] for i in ydata_adj]
            xdata_all.append(xdata)
            ydata_obs_all.append(ydata_obs)
            ydata_adj_all.append(ydata_adj)
        self.plot_datum_window(xdata_all, ydata_obs_all, ydata_adj_all, unique_datum_names)

    def plot_datum_window(self, xdata_all, ydata_obs_all, ydata_adj_all, names):
        fig = plt.figure(figsize=(13, 8))
        i = 0
        for xdata, ydata_obs, ydata_adj in zip(xdata_all, ydata_obs_all, ydata_adj_all):
            a = plt.plot(xdata, ydata_obs, '-o', label=names[i] + '_obs')
            line_color = a[0].get_color()
            b = plt.plot(xdata, ydata_adj, '--o', c=line_color, label=names[i] + '_adj')
            i += 1
        fig.show()
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))


class GravityChangeTable(QtWidgets.QDialog):
    """
    Floating window to show gravity-change results
    :param MainProg:
    :param table:
    :param header:
    """
    def __init__(self, MainProg, table, header, dates=None, full_table=False):
        super(GravityChangeTable, self).__init__()
        self.table = table
        self.header = header
        self.dates = dates
        gravity_change_window = QtWidgets.QWidget()
        gravity_change_window.model = GravityChangeModel(table, header)
        gravity_change_window.table = QtWidgets.QTableView()
        gravity_change_window.table.setModel(gravity_change_window.model)

        # Create buttons and actions
        gravity_change_window.btn0 = QtWidgets.QPushButton('Plot')
        gravity_change_window.btn0.clicked.connect(lambda: self.plot_change_window())
        gravity_change_window.btn1 = QtWidgets.QPushButton('Copy to clipboard')
        gravity_change_window.btn1.clicked.connect(lambda: copy_cells_to_clipboard(MainProg.popup.table))
        if not full_table:
            gravity_change_window.btn2 = QtWidgets.QPushButton('Show full table')
            gravity_change_window.btn2.clicked.connect(lambda: show_full_table(MainProg))
        else:
            gravity_change_window.btn2 = QtWidgets.QPushButton('Show simple table')
            gravity_change_window.btn2.clicked.connect(lambda: MainProg.compute_gravity_change())

        # Locations
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(gravity_change_window.table)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(0)
        # Only show plot button on simple table
        if self.dates:
            hbox.addWidget(gravity_change_window.btn0)
        hbox.addWidget(gravity_change_window.btn1)
        hbox.addWidget(gravity_change_window.btn2)
        vbox.addLayout(hbox)
        gravity_change_window.setLayout(vbox)
        gravity_change_window.setWindowTitle('Gravity Change')
        gravity_change_window.setGeometry(200, 200, 600, 800)
        MainProg.popup = gravity_change_window
        MainProg.popup.show()

    def plot_change_window(self):
        plt.figure(figsize=(12,8))
        cmap = plt.cm.get_cmap('gist_ncar')
        nrows = len(self.dates)
        nstations = len(self.table[0])
        stations = self.table[0]
        # iterate through each station
        for i in range(nstations):
            xdata, ydata = [self.dates[0]], [0]
            for idx, row in enumerate(self.table[1:nrows]):
                if not row[i] == '-999':
                    xdata.append(self.dates[idx+1])
                    ydata.append(float(row[i])+ydata[-1])

            plt.plot(xdata,ydata,'-o', color=cmap(i/(nstations-1)), label=stations[i])
            plt.hold
        plt.ylabel('Gravity change, in µGal')
        plt.show()
        plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        return


def show_full_table(MainProg):
    MainProg.popup.close()
    table = MainProg.compute_gravity_change(full_table=True)
    header = table[0]
    tp_table = list(zip(*table[1:]))
    # GravityChangeTable(MainProg, tp_table, header, full_table=True)
    gravity_change_table = GravityChangeTable(MainProg, tp_table, header, full_table=True)
    return


def show_message(message, title, icon=QtWidgets.QMessageBox.Warning):
    """
    Generic dialog to show a message, with a single 'OK' button.
    :param message: string shown in dialog
    :param title:  string shown in dialog title bar.
    :param icon: Qt icon
    """
    msg = QtWidgets.QMessageBox()
    msg.setAttribute(QtCore.Qt.WA_DeleteOnClose)
    msg.setIcon(icon)
    msg.setText(message)
    msg.setWindowTitle(title)
    msg.show()
    return msg


def rename_dialog(old_name, new_name):
    """
    Dialog called after renaming station in treeview.

    Gives the option to rename stations in the current loop, survey, or throughout the campaign.
    :param old_name: string, old station name
    :param new_name: string, new station name
    :return: integer indicating extent of station rename.
    """
    msg = QtWidgets.QMessageBox()
    q_string = 'Rename all stations from {} to {} in...'.format(old_name, new_name)
    msg.setText(q_string)
    msg.addButton(QtWidgets.QPushButton('Campaign'), 0)
    msg.addButton(QtWidgets.QPushButton('Survey'), 0)
    msg.addButton(QtWidgets.QPushButton('Loop'), 0)
    msg.addButton(QtWidgets.QPushButton('Just this station'),0)
    msg.addButton(QtWidgets.QPushButton('Cancel'), 1)
    method = msg.exec_()
    methods = {0: 'Campaign',
               1: 'Survey',
               2: 'Loop',
               3: 'Station',
               4: 'Cancel'}

    return methods[method]


class TideCoordinatesDialog(QtWidgets.QDialog):

    def __init__(self, lat, lon, elev):
        super(TideCoordinatesDialog, self).__init__()
        self.init_ui(lat, lon, elev)

    def init_ui(self, lat, lon, elev):
        """
        Get coordinates for tide correction.
        """
        latLabel = QtWidgets.QLabel('Latitude')
        lonLabel = QtWidgets.QLabel('Longitude')
        elevLabel = QtWidgets.QLabel('Elevation')
        self.lat = QtWidgets.QLineEdit()
        self.lon = QtWidgets.QLineEdit()
        self.elev = QtWidgets.QLineEdit()

        self.lat.setText('{:.5f}'.format(lat))
        self.lon.setText('{:.5f}'.format(lon))
        self.elev.setText('{:.1f}'.format(elev))

        # create buttons and actions
        self.btn1 = QtWidgets.QPushButton('OK')
        self.btn1.clicked.connect(self.onClicked)
        # locations
        grid = QtWidgets.QGridLayout()
        grid.addWidget(latLabel, 1, 0)
        grid.addWidget(self.lat, 1, 1)
        grid.addWidget(lonLabel, 2, 0)
        grid.addWidget(self.lon, 2, 1)
        grid.addWidget(elevLabel, 3, 0)
        grid.addWidget(self.elev, 3, 1)
        grid.addWidget(self.btn1, 4, 0)
        self.setLayout(grid)
        self.setWindowTitle('Survey coordinates')
        # enterCoordinates.show()
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.accept()


class TideCorrectionDialog(QtWidgets.QDialog):
    def __init__(self):
        super(TideCorrectionDialog, self).__init__()
        self.init_ui()

    def init_ui(self):
        self.msg = QtWidgets.QMessageBox()
        self.msg.setText("Select tide-correction method")
        self.msg.addButton(QtWidgets.QPushButton('Meter-supplied'), QtWidgets.QMessageBox.YesRole)
        btn1 = QtWidgets.QPushButton('Predict')
        btn1.setEnabled(False)
        self.msg.addButton(btn1, QtWidgets.QMessageBox.YesRole)
        self.msg.addButton(QtWidgets.QPushButton('Agnew'), QtWidgets.QMessageBox.YesRole)
        btn2 = QtWidgets.QPushButton('Cancel')
        self.msg.addButton(btn2, QtWidgets.QMessageBox.RejectRole)
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.correction_type = btn.text()
        if self.correction_type == 'Cancel':
            self.reject()
        else:
            self.accept()


class AddDatumFromList(QtWidgets.QInputDialog):

    @classmethod
    def add_datum(cls, station_list):
        dialog = cls()
        text, ok = dialog.getItem(None,
                                  'Input Dialog',
                                  'Datum station:',
                                  station_list)
        if ok:
            return text
        else:
            return None


class LoopTimeThresholdDialog(QtWidgets.QDialog):
    """
    Dialog to specify elapsed time criteria for diving survey into loops.
    """

    def __init__(self):
        super(LoopTimeThresholdDialog, self).__init__()
        self.title = 'Specify elapsed time that indicates a new loop'
        self.left = 100
        self.top = 100
        self.width = 300
        self.height = 140
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Buttons
        ok_button = QtWidgets.QPushButton("Ok")
        ok_button.clicked.connect(self.return_time)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        # Button box
        button_box = QtWidgets.QHBoxLayout()
        button_box.addStretch(1)
        button_box.addWidget(cancel_button)
        button_box.addWidget(ok_button)
        button_widget = QtWidgets.QWidget()
        button_widget.setLayout(button_box)

        main_layout = QtWidgets.QGridLayout()
        main_layout.addWidget(QtWidgets.QLabel("Elapsed time"), 0, 0)
        self.dt_edit = IncrMinuteTimeEdit(QtCore.QTime(8, 0))
        self.dt_edit.setDisplayFormat("hh:mm")

        # self.dt_edit.setDateTime(default_time)
        main_layout.addWidget(self.dt_edit, 0, 1)
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(main_layout)
        # Layout
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addWidget(grid_widget)
        final_layout.addWidget(button_widget)
        self.setLayout(final_layout)

        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def return_time(self):
        self.accept()


class AddTareDialog(QtWidgets.QDialog):
    """
    Dialog to specify date, time and tare.
    """
    def __init__(self, default_time):
        super(AddTareDialog, self).__init__()
        self.title = 'Specify date, time, and magnitude of tare'
        self.left = 100
        self.top = 100
        self.width = 200
        self.height = 140
        self.init_ui(default_time)

    def init_ui(self, default_time):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Buttons
        ok_button = QtWidgets.QPushButton("Ok")
        ok_button.clicked.connect(self.return_date)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        # Button box
        button_box = QtWidgets.QHBoxLayout()
        button_box.addStretch(1)
        button_box.addWidget(cancel_button)
        button_box.addWidget(ok_button)
        button_widget = QtWidgets.QWidget()
        button_widget.setLayout(button_box)

        main_layout = QtWidgets.QGridLayout()
        main_layout.addWidget(QtWidgets.QLabel("Date"), 0, 0)
        self.dt_edit = QtWidgets.QDateTimeEdit()

        self.dt_edit.setDateTime(default_time)
        main_layout.addWidget(self.dt_edit, 0, 1)
        main_layout.addWidget(QtWidgets.QLabel("Tare (microGal)"), 1, 0)
        self.edit_box = QtWidgets.QLineEdit()
        main_layout.addWidget(self.edit_box, 1, 1)
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(main_layout)
        # Layout
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addWidget(grid_widget)
        final_layout.addWidget(button_widget)
        self.setLayout(final_layout)

        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def return_date(self):
        self.accept()


class SelectAbsg(QtWidgets.QDialog):
    """
    Dialog to show absolute-gravity values from *.project.txt files. The user can select the files to import as Datums.
    """
    def __init__(self, path):
        super(SelectAbsg, self).__init__()
        self.title = 'Select directory with Abs. g files (.project.txt)'
        self.left = 100
        self.top = 100
        self.width = 800
        self.height = 480
        self.new_datums = []
        self.path = path

        self.splitter_window = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self)
        self.tree_model = QtWidgets.QDirModel()
        self.tree = QtWidgets.QTreeView()
        self.table_model = DatumTableModel()
        self.tree_model.setFilter(QtCore.QDir.Dirs | QtCore.QDir.NoDotAndDotDot)
        self.table = QtWidgets.QTableView()

        self.ProxyModel = CustomSortingModel(self)
        self.table.setModel(self.ProxyModel)
        self.table.setSortingEnabled(True)
        self.init_ui()

    def export_and_close(self):
        for i in range(self.ProxyModel.rowCount()):
            ndi = self.ProxyModel.index(i, 0)
            nd = self.ProxyModel.data(ndi, role=QtCore.Qt.UserRole)
            chk = self.ProxyModel.data(ndi, role=QtCore.Qt.CheckStateRole)
            if chk == 2:
                self.new_datums.append(nd)
        self.accept()

    def closeEvent(self, QCloseEvent):
        return self.close()

    def init_ui(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # File-tree view and model
        self.tree.setModel(self.tree_model)
        self.tree.expand(self.tree_model.index(self.path))
        self.tree.scrollTo(self.tree_model.index(self.path))
        self.tree.setAnimated(False)
        self.tree.setIndentation(20)
        self.tree.setSortingEnabled(True)
        self.tree.setWindowTitle("Dir View")
        self.tree.resize(800, 480)

        # Buttons
        load_button = QtWidgets.QPushButton("Load")
        load_button.clicked.connect(self.load_a10_data)
        ok_button = QtWidgets.QPushButton("Import")
        ok_button.clicked.connect(self.export_and_close)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        # Button under tree view
        button_box_left = QtWidgets.QHBoxLayout()
        button_box_left.addStretch(1)
        button_box_left.addWidget(load_button)

        # Button under table
        button_box_right = QtWidgets.QHBoxLayout()
        button_box_right.addStretch(1)
        button_box_right.addWidget(cancel_button)
        button_box_right.addWidget(ok_button)

        self.ProxyModel.setSourceModel(self.table_model)
        self.tree.resizeColumnToContents(0)

        # Hide file size, date modified columns
        self.tree.setColumnHidden(1, True)
        self.tree.setColumnHidden(2, True)
        self.tree.setColumnHidden(3, True)

        # Hide column of residuals
        self.table.setColumnHidden(6, True)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        # Move date column to second from left
        self.table.horizontalHeader().moveSection(3, 1)

        # Layout
        final_layout = QtWidgets.QVBoxLayout()
        sublayout_left = QtWidgets.QVBoxLayout()
        sublayout_right = QtWidgets.QVBoxLayout()
        sublayout_left.addWidget(self.tree)
        sublayout_left.addLayout(button_box_left)
        sublayout_right.addWidget(self.table)
        sublayout_right.addLayout(button_box_right)
        left_widget = QtWidgets.QWidget()
        left_widget.setLayout(sublayout_left)
        right_widget = QtWidgets.QWidget()
        right_widget.setLayout(sublayout_right)
        self.splitter_window.addWidget(left_widget)
        self.splitter_window.addWidget(right_widget)
        self.splitter_window.setSizes([240, 560])
        final_layout.addWidget(self.splitter_window)
        self.setLayout(final_layout)

        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def load_a10_data(self):
        """
        Parses *.project.txt files in the selected paths. Populates the dialog table model directly.
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        files_found = False
        idxs = self.tree.selectedIndexes()
        self.table_model.clearDatums()
        for i in idxs:
            if i.model().isDir(i):
                path = str(i.model().filePath(i))
                for dirname, _, fileList in os.walk(path):
                    for name in fileList:
                        if '.project.txt' in name:
                            files_found = True
                            d = a10.A10(os.path.join(dirname, name))
                            datum = Datum(d.stationname,
                                          g=float(d.gravity),
                                          sd=float(d.setscatter),
                                          date=d.date,
                                          meas_height=float(d.transferht),
                                          gradient=float(d.gradient),
                                          checked=0)
                            self.table_model.insertRows(datum, 1)
                            self.path = path
        QtWidgets.QApplication.restoreOverrideCursor()
        if not files_found:
            self.msg = show_message('No *.project.txt files found in the selected directories.', 'Import error')
