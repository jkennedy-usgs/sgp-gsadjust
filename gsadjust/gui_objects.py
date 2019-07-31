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
import datetime as dt
import os

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from PyQt5 import QtGui, QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.dates import date2num
from matplotlib.figure import Figure

import a10
import data_analysis
from data_objects import Datum
from pyqt_models import GravityChangeModel, DatumTableModel, MeterCalibrationModel
from utils import *

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


class Overwrite(QtWidgets.QMessageBox):
    def __init__(self):
        super(Overwrite, self).__init__()
        self.overwrite = False
        # self.msg = QtWidgets.QMessageBox()
        self.setText("Overwrite existing data?")
        self.addButton(QtWidgets.QPushButton('Overwrite'), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton('Cancel'), QtWidgets.QMessageBox.RejectRole)
        self.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        if btn.text() == 'Overwrite':
            self.accept()
        else:
            self.reject()


class CoordinatesTable(QtWidgets.QDialog):
    """
    Shows table of coordinates, which can be edited or copy/pasted. Coordinates provided by self.coords(), called
    by the calling routine.
    """

    def __init__(self, coords):
        super(CoordinatesTable, self).__init__()
        vlayout = QtWidgets.QVBoxLayout()

        self.setWindowModality(QtCore.Qt.ApplicationModal)
        self.setWindowTitle('Station coordinates')

        ok_button = QtWidgets.QPushButton("Ok")
        ok_button.clicked.connect(self.accept)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        # Button under table
        button_box = QtWidgets.QWidget()
        bbox_layout = QtWidgets.QHBoxLayout()
        bbox_layout.addStretch(1)
        bbox_layout.addWidget(cancel_button)
        bbox_layout.addWidget(ok_button)
        button_box.setLayout(bbox_layout)

        self.table = QtWidgets.QTableWidget()
        self.table.setSizeAdjustPolicy(
            QtWidgets.QAbstractScrollArea.AdjustToContents)
        vlayout.addWidget(self.table)
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['Station', 'Latitude', 'Longitude', 'Elevation'])
        for k, v in coords.items():
            row = self.table.rowCount()
            self.table.insertRow(row)
            table_item_station = QtWidgets.QTableWidgetItem(k)
            table_item_lat = QtWidgets.QTableWidgetItem(str(v[1]))
            table_item_long = QtWidgets.QTableWidgetItem(str(v[0]))
            table_item_elev = QtWidgets.QTableWidgetItem(str(v[2]))
            self.table.setItem(row, 0, table_item_station)
            self.table.setItem(row, 1, table_item_lat)
            self.table.setItem(row, 2, table_item_long)
            self.table.setItem(row, 3, table_item_elev)

        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        # self.adjustSize()

        vlayout.addWidget(button_box)
        self.setLayout(vlayout)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.setSizePolicy(sizePolicy)

    def keyPressEvent(self, e):
        """
        Handle copy/paste. Without this, the default is to copy only one cell. No sanity checks.
        :param e: key pressed
        :return: None
        """
        if e.modifiers() & QtCore.Qt.ControlModifier:
            selected = self.table.selectedRanges()
            if e.key() == QtCore.Qt.Key_C:  # copy
                s = ''
                for r in range(selected[0].topRow(), selected[0].bottomRow() + 1):
                    row_text = ''
                    for c in range(self.table.columnCount()):
                        row_text += self.table.item(r, c).text() + '\t'
                    rt = row_text[:-1]  # Remove trailing \t
                    rt += '\n'
                    s += rt
                self.sys_clip = QtWidgets.QApplication.clipboard()
                self.sys_clip.setText(s)

            if e.key() == QtCore.Qt.Key_V:
                self.sys_clip = QtWidgets.QApplication.clipboard()
                s = self.sys_clip.text()
                rows = s.split('\n')
                for idx, r in enumerate(rows):
                    if r is not '':
                        elems = r.split('\t')
                        table_item_station = QtWidgets.QTableWidgetItem(elems[0])
                        table_item_lat = QtWidgets.QTableWidgetItem(elems[1])
                        table_item_long = QtWidgets.QTableWidgetItem(elems[2])
                        table_item_elev = QtWidgets.QTableWidgetItem(elems[3])
                        self.table.setItem(idx, 0, table_item_station)
                        self.table.setItem(idx, 1, table_item_lat)
                        self.table.setItem(idx, 2, table_item_long)
                        self.table.setItem(idx, 3, table_item_elev)

    def coords(self):
        """
        Puts table coordinates into a dict.
        :return: A dict with key = station name, value = tuple[lat, long, elev)
        """
        c = dict()
        for i in range(self.table.rowCount()):
            c[self.table.item(i, 0).text()] = (float(self.table.item(i, 2).text()),
                                               float(self.table.item(i, 1).text()),
                                               float(self.table.item(i, 3).text()))
        return c


class MeterType(QtWidgets.QMessageBox):
    def __init__(self):
        super(MeterType, self).__init__()
        self.setText("Choose meter file to import")
        self.addButton(QtWidgets.QPushButton(' CG-3, CG-5 '), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(' CG-6 '), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(' CG-6 Tsoft '), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(' Burris '), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(' CSV '), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(' Cancel '), QtWidgets.QMessageBox.RejectRole)
        self.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        meter = {' CG-3, CG-5 ': 'Scintrex',
                 ' CG-6 ': 'CG6',
                 ' CG-6 Tsoft ': 'CG6Tsoft',
                 ' Burris ': 'Burris',
                 ' CSV ': 'csv',
                 ' Cancel ': 'Cancel'}
        self.meter_type = meter[btn.text()]
        if self.meter_type == 'Cancel':
            self.reject()
        else:
            self.accept()


class LoopOptions(QtWidgets.QDialog):
    def __init__(self, loops, parent=None):
        super(LoopOptions, self).__init__(parent)
        self.setGeometry(50, 50, 350, 350)
        self.loops = loops
        self.init_ui()

    def init_ui(self):
        # create buttons and actions
        cancel_button = QtWidgets.QPushButton('Cancel')
        cancel_button.clicked.connect(self.close)
        ok_button = QtWidgets.QPushButton('OK')
        ok_button.clicked.connect(self.set_loop_options)

        buttonBox = QtWidgets.QDialogButtonBox(QtCore.Qt.Horizontal)
        buttonBox.addButton(cancel_button, QtWidgets.QDialogButtonBox.ActionRole)
        buttonBox.addButton(ok_button, QtWidgets.QDialogButtonBox.ActionRole)

        self.operator_edit = QtWidgets.QLineEdit(self.loops[0].oper)
        self.meter_edit = QtWidgets.QLineEdit(self.loops[0].meter)
        self.comment_edit = QtWidgets.QTextEdit()
        self.comment_edit.setPlainText(self.loops[0].comment)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Operator'), 1, 0)
        grid.addWidget(self.operator_edit, 1, 1)
        grid.addWidget(QtWidgets.QLabel('Meter ID'), 2, 0)
        grid.addWidget(self.meter_edit, 2, 1)
        grid.addWidget(QtWidgets.QLabel('Comment'), 3, 0)
        grid.addWidget(self.comment_edit, 4, 0, 1, 2)
        grid.addWidget(buttonBox, 5, 0, 1, 2)
        self.setLayout(grid)

        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def set_loop_options(self):
        self.accept()


class AdjustOptions(QtWidgets.QDialog):
    """
    Dialog to set network adjustment options.
    """

    def __init__(self, survey_str, options, parent=None):
        super(AdjustOptions, self).__init__(parent)
        self.setGeometry(50, 50, 350, 350)
        self.ao = options
        self.surveys_to_update = ''
        self.update_options = []
        self.drift_temp_chk = QtWidgets.QCheckBox('Model temperature drift, polynomial degree:')
        self.sigma_factor_chk = QtWidgets.QCheckBox('Std. dev. multiplier')
        self.sigma_add_chk = QtWidgets.QCheckBox('Add to std. dev.')
        self.sigma_min_chk = QtWidgets.QCheckBox('Minimum std. dev.')
        self.cal_coeff_chk = QtWidgets.QCheckBox('Calculate relative meter calibration coefficient')
        self.cal_coeff_chk.stateChanged.connect(self.calc_cal_coeff_checked_or_unchecked)
        self.alpha_text = QtWidgets.QLabel('Significance level for global model test')
        self.woutfiles_chk = QtWidgets.QCheckBox('Write output files')
        self.cal_coeff_specify_chk = QtWidgets.QCheckBox('Specify calibration coefficient')
        self.cal_coeff_specify_chk.stateChanged.connect(self.specify_cal_coeff_checked_or_unchecked)

        self.cal_coeff_table = QtWidgets.QTableView()

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
            try:
                self.cal_coeff_specify_chk.setChecked(self.ao.specify_cal_coeff)
            except:
                self.ao.specify_cal_coeff = False
            if not hasattr(self.ao, 'meter_cal_dict'):
                self.ao.meter_cal_dict = init_cal_coeff_dict(self.parent().obsTreeModel)
            elif self.ao.meter_cal_dict is None:
                self.ao.meter_cal_dict = init_cal_coeff_dict(self.parent().obsTreeModel)
            self.cal_coeff_model = MeterCalibrationModel()
            for k, v in self.ao.meter_cal_dict.items():
                self.cal_coeff_model.appendRow([QtGui.QStandardItem(k), QtGui.QStandardItem(str(v))])

            self.cal_coeff_table.setModel(self.cal_coeff_model)
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
            grid.addWidget(self.cal_coeff_specify_chk, 6, 0)
            grid.addWidget(self.cal_coeff_table, 7, 0)
            grid.addWidget(self.alpha_text, 8, 0)
            grid.addWidget(self.alpha_edit, 8, 1)
            grid.addWidget(self.woutfiles_chk, 9, 0)
            grid.addWidget(buttonBox, 10, 0)

            self.setLayout(grid)
            self.setWindowTitle('Network adjustment options')
            self.setWindowModality(QtCore.Qt.ApplicationModal)
        else:
            self.msg = show_message('Please load a survey first', 'Network adjustment options')

    def calc_cal_coeff_checked_or_unchecked(self, state):
        if state == 2:  # checked
            self.cal_coeff_specify_chk.setEnabled(False)
            self.cal_coeff_table.setEnabled(False)
        else:
            self.cal_coeff_specify_chk.setEnabled(True)
            self.cal_coeff_table.setEnabled(True)

    def specify_cal_coeff_checked_or_unchecked(self, state):
        if state == 2:  # checked
            self.cal_coeff_chk.setEnabled(False)
        else:
            self.cal_coeff_chk.setEnabled(True)

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
            self.ao.specify_cal_coeff = False
        else:
            self.ao.cal_coeff = False
        if self.cal_coeff_specify_chk.isChecked():
            self.ao.specify_cal_coeff = True
            self.ao.cal_coeff = False
            for i in range(self.cal_coeff_model.rowCount()):
                meter = self.cal_coeff_model.itemFromIndex(self.cal_coeff_model.index(i, 0))
                calval = self.cal_coeff_model.itemFromIndex(self.cal_coeff_model.index(i, 1))
                self.ao.meter_cal_dict[meter.text()] = float(calval.text())
        else:
            self.ao.specify_cal_coeff = False
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
            self.setTime(QtCore.QTime(hours - 1, 60 + minutes))
        if minutes < 60:
            self.setTime(QtCore.QTime(hours, minutes))
        else:
            self.setTime(QtCore.QTime(hours + 1, 60 % minutes))


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


class AboutDialog(QtWidgets.QDialog):
    def __init__(self, obsTreeModel):
        super(AboutDialog, self).__init__()

        msg1 = '<html>GSadjust, a product of the USGS Southwest Gravity Program<br>' + \
               '<a href ="http://go.usa.gov/xqBnQ">http://go.usa.gov/xqBnQ</a>' + \
               '<br><br><a href ="https://github.com/jkennedy-usgs/sgp-gsadjust">' + \
               'https://github.com/jkennedy-usgs/sgp-gsadjust</a>' + \
               '<br><a href="mailto:jkennedy@usgs.gov">jkennedy@usgs.gov</a>'
        _, ok = QtWidgets.QMessageBox.about(None, "GSadust", msg1)


class VerticalGradientDialog(QtWidgets.QInputDialog):
    def __init__(self, default_interval):
        super(VerticalGradientDialog, self).__init__()
        self.show_dialog(default_interval)

    def show_dialog(self, default_interval):
        text, ok = self.getDouble(None, "Vertical-gradient interval",
                                                    "Interval, in cm:",
                                                    default_interval,
                                                    0, 200, 1)



class FigureDatumComparisonTimeSeries(QtWidgets.QDialog):
    def __init__(self, obsTreeModel):
        super(FigureDatumComparisonTimeSeries, self).__init__()
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
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
        plt.ylabel('Relative gravity, in \u00b5Gal')
        plt.axes().xaxis_date()
        fig.show()


class GravityChangeTable(QtWidgets.QDialog):
    """
    Floating window to show gravity-change results
    :param MainProg:
    :param table:
    :param header:
    """

    def __init__(self, MainProg, table, header, dates=None, full_table=False):
        super(GravityChangeTable, self).__init__()
        self.header, self.table, self.dates = data_analysis.compute_gravity_change(MainProg.obsTreeModel)
        self.full_header, self.full_table, _ = data_analysis.compute_gravity_change(MainProg.obsTreeModel,
                                                                                    full_table=True)
        self.coords = MainProg.obsTreeModel.station_coords
        gravity_change_window = QtWidgets.QWidget()
        if full_table:
            gravity_change_window.model = GravityChangeModel(self.full_header, self.full_table, full_table=True)
        else:
            gravity_change_window.model = GravityChangeModel(self.header, self.table)
        gravity_change_window.table = QtWidgets.QTableView()
        gravity_change_window.table.setModel(gravity_change_window.model)

        # Create buttons and actions
        gravity_change_window.btn00 = QtWidgets.QPushButton('Map')
        gravity_change_window.btn00.clicked.connect(self.map_change_window)
        gravity_change_window.btn0 = QtWidgets.QPushButton('Plot')
        gravity_change_window.btn0.clicked.connect(lambda: self.plot_change_window())
        gravity_change_window.btn1 = QtWidgets.QPushButton('Copy to clipboard')
        gravity_change_window.btn1.clicked.connect(lambda: copy_cells_to_clipboard(MainProg.popup.table))
        if not full_table:
            gravity_change_window.btn2 = QtWidgets.QPushButton('Show full table')
            gravity_change_window.btn2.clicked.connect(lambda: show_full_table(MainProg))
        else:
            gravity_change_window.btn2 = QtWidgets.QPushButton('Show simple table')
            gravity_change_window.btn2.clicked.connect(lambda: show_simple_table(MainProg))

        # Locations
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(gravity_change_window.table)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(0)
        # Only show plot button on simple table
        if self.dates:
            hbox.addWidget(gravity_change_window.btn0)
            hbox.addWidget(gravity_change_window.btn00)
        hbox.addWidget(gravity_change_window.btn1)
        hbox.addWidget(gravity_change_window.btn2)
        vbox.addLayout(hbox)
        gravity_change_window.setLayout(vbox)
        gravity_change_window.setWindowTitle('Gravity Change')
        gravity_change_window.setGeometry(200, 200, 600, 800)
        MainProg.popup = gravity_change_window
        MainProg.popup.show()

    def map_change_window(self):

        self.win = GravityChangeMap(self.table, self.header, self.coords, self.full_table, self.full_header)
        self.win.show()

    def plot_change_window(self):
        plt.figure(figsize=(12, 8))
        cmap = plt.cm.get_cmap('gist_ncar')
        ncols = len(self.dates)
        nstations = len(self.table[0])
        stations = self.table[0]
        # iterate through each station
        for i in range(nstations):
            xdata, ydata = [], []
            for idx, col in enumerate(self.table[1:ncols]):
                if not col[i] == '-999':
                    if not ydata:
                        ydata.append(0)
                        ydata.append(float(col[i]))
                        xdata.append(self.dates[idx])
                        xdata.append(self.dates[idx + 1])
                    else:
                        ydata.append(float(col[i]) + ydata[-1])
                        xdata.append(self.dates[idx + 1])

            plt.plot(xdata, ydata, '-o', color=cmap(i / (nstations - 1)), label=stations[i])
            # plt.hold
        plt.ylabel('Gravity change, in µGal')

        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
        plt.show()


class GravityChangeMap(QtWidgets.QDialog):
    station_label = None

    def __init__(self, table, header, coords, full_table, full_header, parent=None):
        super(GravityChangeMap, self).__init__(parent)
        # a figure instance to plot on
        self.table = table
        self.header = header
        self.coords = coords
        self.full_table = full_table
        self.full_header = full_header
        self.figure = Figure(figsize=(10, 8))
        self.n_surveys = (len(self.table) - 1) / 2
        self.surveys = self.get_survey_dates(header)
        self.axlim = self.get_axis_lims(coords)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.time_widget = QtWidgets.QWidget()
        self.btnIncremental = QtWidgets.QRadioButton("Incremental", self)
        self.btnIncremental.setChecked(True)
        self.btnReference = QtWidgets.QRadioButton("Change relative to", self)
        self.drpReference = QtWidgets.QComboBox()
        self.btnIncremental.toggled.connect(self.plot)
        self.btnReference.toggled.connect(self.plot)
        bbox = QtWidgets.QHBoxLayout()
        bbox.addWidget(self.btnIncremental)
        bbox.addSpacing(10)
        bbox.addWidget(self.btnReference)
        bbox.addWidget(self.drpReference)
        bbox.addStretch(1)
        self.time_widget.setLayout(bbox)

        self.basemap_widget = QtWidgets.QWidget()
        self.cbBasemap = QtWidgets.QCheckBox("Show Basemap", self)
        self.cbBasemap.setChecked(False)
        self.cbBasemap.stateChanged.connect(self.plot)
        self.drpBasemap = QtWidgets.QComboBox()
        self.drpBasemap.currentIndexChanged.connect(self.plot)
        bbox = QtWidgets.QHBoxLayout()
        bbox.addWidget(self.cbBasemap)
        bbox.addSpacing(10)
        bbox.addWidget(self.drpBasemap)
        bbox.addStretch(1)
        self.basemap_widget.setLayout(bbox)

        self.units_widget = QtWidgets.QWidget()
        self.cbUnits = QtWidgets.QCheckBox("Show change in meters of water", self)
        self.cbUnits.setChecked(False)
        self.cbUnits.stateChanged.connect(self.plot)
        self.sliderColorRange = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sliderColorRange.valueChanged.connect(self.update_plot)
        self.sliderColorRange.setMinimum(1)
        self.sliderColorRange.setMaximum(15)
        self.sliderColorRange.setValue(100)
        self.sliderColorRange.resize(100,20)
        self.sliderColorRange.setTickInterval(10)
        bbox = QtWidgets.QHBoxLayout()
        bbox.addWidget(self.cbUnits)
        bbox.addSpacing(40)
        bbox.addWidget(QtWidgets.QLabel('Color range'))
        bbox.addWidget(self.sliderColorRange)
        bbox.addStretch(1)
        self.units_widget.setLayout(bbox)

        # Slider
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.valueChanged.connect(self.update_plot)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksBelow)

        for s in self.surveys:
            self.drpReference.addItem(s)

        self.maps = [
            'NatGeo_World_Map',
            'USA_Topo_Maps',
            'World_Imagery',
            'World_Shaded_Relief',
            'World_Street_Map',
            'World_Topo_Map']

        for map in self.maps:
            self.drpBasemap.addItem(map)

        self.slider_label = QtWidgets.QLabel(self.header[1])
        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.time_widget)
        layout.addWidget(self.basemap_widget)
        layout.addWidget(self.units_widget)
        layout.addStretch(1)
        layout.addWidget(self.canvas)
        layout.addWidget(self.slider_label)

        layout.addWidget(self.slider)
        layout.addStretch(1)
        self.setLayout(layout)

        self.plot()

    def get_survey_dates(self, header):
        dates = []
        first = True
        for title in header:
            if first:
                first = False
                continue
            d = title.split('_to_')
            dates.append(d[0])
            dates.append(d[1])
        return sorted(list(set(dates)))

    def update_plot(self):
        if hasattr(self, 'points'):
            self.cb.remove()
            self.points.remove()
            x, y, d, name = self.get_data()
            clim = self.get_color_lims(self.table)
            self.points = self.ax.scatter(x, y, c=d, s=200, vmin=clim[0], vmax=clim[1], cmap="RdBu",
                                          picker=5,
                                          zorder=10,
                                          transform=ccrs.Geodetic())
            self.cb = self.figure.colorbar(self.points)
            if not self.cbUnits.isChecked():
                self.cb.set_label('Gravity change, in µGal', fontsize=16)
            elif self.cbUnits.isChecked():
                self.cb.set_label('Aquifer-storage change,\n in meters of water', fontsize=16)
            # self.cb.set_clim(vmin=clim[0], vmax=clim[1])
            # self.cb.draw_all()
            self.points.name = name
            self.slider_label.setText(self.get_name())
            self.ax.set_title(self.get_name(), fontsize=16, fontweight='bold')
            self.canvas.draw()

    def get_data(self):
        x, y, d, name = [], [], [], []

        if self.btnIncremental.isChecked():
            data = self.table[self.slider.value()]
            stations = self.table[0]
            for sta_idx, sta in enumerate(stations):
                datum = float(data[sta_idx])
                if datum > -998:
                    x.append(self.coords[sta][0])
                    y.append(self.coords[sta][1])
                    if self.cbUnits.isChecked():
                        datum /= 41.9
                    d.append(datum)
                    name.append(sta)

        elif self.btnReference.isChecked():
            ref_survey = self.drpReference.currentData(role=QtCore.Qt.DisplayRole)
            ref_col_idx = self.full_header.index(ref_survey)
            current_survey = self.surveys[self.slider.value() - 1]
            current_col_idx = self.full_header.index(current_survey)
            for r in self.full_table:
                sta = r[0]
                ref_g = float(r[ref_col_idx])
                surv_g = float(r[current_col_idx])
                if ref_g > -998 and surv_g > -998:
                    x.append(self.coords[sta][0])
                    y.append(self.coords[sta][1])
                    if current_survey < ref_survey:
                        datum = (ref_g - surv_g)
                    else:
                        datum = (surv_g - ref_g)
                    if self.cbUnits.isChecked():
                        d.append(datum / 41.9)
                    else:
                        d.append(datum)
                    name.append(sta)

        return x, y, d, name

    def get_name(self):
        if self.btnIncremental.isChecked():
            name = self.header[self.slider.value()]
            return name.replace('_', ' ')
        elif not self.btnIncremental.isChecked():
            ref_survey = self.drpReference.currentData(role=QtCore.Qt.DisplayRole)
            current_survey = self.surveys[self.slider.value() - 1]
            if current_survey < ref_survey:
                return current_survey + ' to ' + ref_survey
            else:
                return ref_survey + ' to ' + current_survey

    def plot(self):
        clim = self.get_color_lims(self.table)

        if self.btnIncremental.isChecked():
            self.slider.setRange(1, self.n_surveys)
        elif not self.btnIncremental.isChecked():
            self.slider.setRange(1, self.n_surveys + 1)
        self.slider.setValue(1)

        self.figure.clf()

        map_center = (self.axlim[0] + self.axlim[1]) / 2
        self.ax = self.figure.add_subplot(1, 1, 1,
                                          projection=ccrs.AlbersEqualArea(map_center))  # self.stamen_terrain.crs)
        self.ax.clear()

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.points = []
        x, y, d, names = self.get_data()
        self.points = self.ax.scatter(x, y, c=d, s=200, vmin=clim[0], vmax=clim[1], cmap="RdBu",
                                      picker=5,
                                      zorder=10,
                                      transform=ccrs.Geodetic())
        self.points.name = names

        QtWidgets.QApplication.restoreOverrideCursor()
        # self.figure.colorbar(point)
        self.ax.set_title(self.get_name(), fontsize=16, fontweight='bold')
        self.cb = self.figure.colorbar(self.points)
        self.ax.set_xlabel('Distance, in meters', fontsize=16)
        if not self.cbUnits.isChecked():
            self.cb.set_label('Gravity change, in µGal', fontsize=16)
        elif self.cbUnits.isChecked():
            self.cb.set_label('Aquifer-storage change,\n in meters of water', fontsize=16)

        self.show_background()
        self.ax.set_extent(self.axlim, crs=ccrs.Geodetic())
        self.slider_label.setText(self.get_name())

        # refresh canvas
        self.figure.canvas.mpl_connect('pick_event', self.show_point_label)
        self.canvas.draw()

    def show_background(self):
        if self.cbBasemap.isChecked():
            self.stamen_terrain = cimgt.GoogleTiles(url='https://server.arcgisonline.com/ArcGIS/rest/services/' +
                                                        self.maps[self.drpBasemap.currentIndex()] +
                                                        '/MapServer/tile/{z}/{y}/{x}.jpg')  # cimgt.Stamen('terrain-background')
            self.ax.add_image(self.stamen_terrain, 12)

    def show_point_label(self, event):
        """
        Shows the station name in the upper left of the drift plot when a line is clicked.
        :param event: Matplotlib event
        :param axes: Current axes (differs for none|netadj|roman vs continuous)
        """
        thispoint = event.artist
        if self.station_label is not None:
            self.station_label.set_text('')
        self.station_label = self.ax.text(0.1, 0.95, thispoint.name[event.ind[0]], horizontalalignment='center',
                                          verticalalignment='center',
                                          transform=self.ax.transAxes)
        self.canvas.draw()

    def get_color_lims(self, table):
        # margin = 0.1
        # row_minmax = []
        # first = True
        # for row in table:
        #     if first:
        #         first = False
        #         continue
        #     valid_data = [float(d) for d in row if float(d) > -998]
        #     row_minmax.append(abs(min(valid_data)))
        #     row_minmax.append(abs(max(valid_data)))
        #
        # cmax = max(row_minmax)
        # clim = (cmax * -1 - cmax * margin, cmax + cmax * margin)
        slider_val = self.sliderColorRange.value() * 10
        clim = (slider_val * -1, slider_val)
        if self.cbUnits.isChecked():
            clim = (clim[0] / 41.9, clim[1] / 41.9)
        return clim

    def get_axis_lims(self, coords):
        margin = 0.25
        x, y = [], []
        for c in coords.values():
            x.append(c[0])
            y.append(c[1])

        xrange = abs(max(x) - min(x))
        yrange = abs(max(y) - min(y))
        xmin = min(x) - xrange * margin
        xmax = max(x) + xrange * margin
        ymin = min(y) - yrange * margin
        ymax = max(y) + yrange * margin
        return (xmin, xmax, ymin, ymax)


def show_full_table(MainProg):
    MainProg.popup.close()
    header, table, dates = data_analysis.compute_gravity_change(MainProg.obsTreeModel, full_table=True)
    # header = table[0]
    tp_table = list(zip(*table[1:]))
    # GravityChangeTable(MainProg, tp_table, header, full_table=True)
    GravityChangeTable(MainProg, tp_table, header, dates, full_table=True)
    return


def show_simple_table(MainProg):
    MainProg.popup.close()
    header, table, dates = data_analysis.compute_gravity_change(MainProg.obsTreeModel, full_table=False)
    GravityChangeTable(MainProg, table, header, dates, full_table=False)
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
    msg.addButton(QtWidgets.QPushButton('Just this station'), 0)
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
    Dialog to show absolute-gravity values from *.project.txt files. The user
    can select the files to import as Datums.
    """

    def __init__(self, path):
        super(SelectAbsg, self).__init__()
        self.title = 'Select directory with Abs. g files (.project.txt)'
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 800
        self.new_datums = []
        self.path = path

        self.splitter_window = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self)
        self.tree_model = QtWidgets.QDirModel()
        self.tree = QtWidgets.QTreeView()
        self.table_model = DatumTableModel()
        self.tree_model.setFilter(QtCore.QDir.Dirs | QtCore.QDir.NoDotAndDotDot)
        self.table = QtWidgets.QTableView()

        self.ProxyModel = QtCore.QSortFilterProxyModel()
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
        self.tree.setCurrentIndex(self.tree_model.index(self.path))
        self.tree.setAnimated(False)
        self.tree.setIndentation(20)
        self.tree.setSortingEnabled(True)
        self.tree.setWindowTitle("Dir View")
        self.tree.resize(800, 480)

        # Buttons and checkbox
        self.load_unpublished_cb = QtWidgets.QCheckBox('Ignore unpublished')
        self.load_unpublished_cb.setChecked(True)
        self.load_button = QtWidgets.QPushButton("Load")
        self.load_button.clicked.connect(self.load_a10_data)
        self.ok_button = QtWidgets.QPushButton("Import")
        self.ok_button.clicked.connect(self.export_and_close)
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.close)

        # Button under tree view
        button_box_left = QtWidgets.QHBoxLayout()
        button_box_left.addWidget(self.load_unpublished_cb)
        button_box_left.addStretch(1)
        button_box_left.addWidget(self.load_button)

        # Button under table
        button_box_right = QtWidgets.QHBoxLayout()
        button_box_right.addStretch(1)
        button_box_right.addWidget(self.cancel_button)
        button_box_right.addWidget(self.ok_button)

        self.ProxyModel.setSourceModel(self.table_model)
        self.tree.resizeColumnToContents(0)

        # Hide file size, date modified columns
        self.tree.setColumnHidden(1, True)
        self.tree.setColumnHidden(2, True)
        self.tree.setColumnHidden(3, True)

        # Hide column of residuals
        self.table.setColumnHidden(7, True)
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
        idxs = self.tree.selectedIndexes()
        self.table_model.clearDatums()
        for i in idxs:
            if i.model().isDir(i):
                path = str(i.model().filePath(i))
                files_found = self.append_datums(path)
        QtWidgets.QApplication.restoreOverrideCursor()
        if not files_found:
            self.msg = show_message('No *.project.txt files found in the selected directories.', 'Import error')

    def append_datums(self, path):
        files_found = False
        for dirname, _, fileList in os.walk(path):
            if self.load_unpublished_cb.isChecked() and dirname.find('unpublished') != -1:
                continue
            else:
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
                        datum.n_sets = d.processed
                        datum.time = d.time
                        self.table_model.insertRows(datum, 1)
                        self.path = path
        return files_found
