"""
gui/dialogs.py
==============

Dialog GUI objects.
--------------------------------------------------------------------------------

Major GUI objects (tabs, table views) are in the tab_... files. This module has
primarily pop-up dialogs used to set network adjustment settings, show gravity
change over time, etc. Major dialogs are written as classes and instantiated in
GSadjust.py. Minor dialogs are written as functions and called directly.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""


import logging
import os
import datetime

from PyQt5 import QtCore, QtGui, QtWidgets
from .messages import MessageBox
from .utils import copy_cells_to_clipboard
from .widgets import IncrMinuteTimeEdit, ProgressBar
from .map import GravityChangeMap
from ..data import Datum, analysis
from ..file import a10
from ..models import (
    DatumTableModel,
    GravityChangeModel,
    MeterCalibrationModel,
    CheckBoxDelegate,
)
from ..utils import init_cal_coeff_dict
from ..data.nwis import search_nwis, nwis_get_data, filter_ts
from ..data.analysis import compute_gravity_change
from ..plots import PlotNwis


class CoordinatesTable(QtWidgets.QDialog):
    """
    Shows table of coordinates, which can be edited or copy/pasted. Coordinates
    provided by self.coords(), called by the calling routine.

    Parameters
    ----------
    coords: dict (str: tuple)
        station: (lat, lon, elev)

    """

    def __init__(self, coords):
        super(CoordinatesTable, self).__init__()
        vlayout = QtWidgets.QVBoxLayout()
        self.sys_clip = QtWidgets.QApplication.clipboard()
        self.setWindowModality(QtCore.Qt.ApplicationModal)
        self.setWindowTitle("Station coordinates")

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
        self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        vlayout.addWidget(self.table)
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(
            ["Station", "Latitude", "Longitude", "Elevation"]
        )
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

        vlayout.addWidget(button_box)
        self.setLayout(vlayout)
        policy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )
        self.setSizePolicy(policy)

    def keyPressEvent(self, e):
        """
        Handle copy/paste. Without this, the default is to copy only one cell.
        No sanity checks.

        Parameters
        ----------
        e : key pressed

        """
        if e.modifiers() & QtCore.Qt.ControlModifier:
            selected = self.table.selectedRanges()
            if e.key() == QtCore.Qt.Key_C:  # copy
                s = ""
                for r in range(selected[0].topRow(), selected[0].bottomRow() + 1):
                    row_text = ""
                    for c in range(self.table.columnCount()):
                        row_text += self.table.item(r, c).text() + "\t"
                    rt = row_text[:-1]  # Remove trailing \t
                    rt += "\n"
                    s += rt
                clipboard = QtWidgets.QApplication.clipboard()
                clipboard.setText(s)

            if e.key() == QtCore.Qt.Key_V:
                clipboard = QtWidgets.QApplication.clipboard()
                s = clipboard.text()
                rows = s.split("\n")[:-1]  # Excel includes trailing "\n"
                # Check the clipboard is the right number of columns. The number of rows
                # can be different: rows will be added or deleted as needed.
                if len(rows[0].split('\t')) != 4:
                    return
                for idx, r in enumerate(rows):
                    if r == "":
                        continue
                    elems = r.split("\t")
                    table_item_station = QtWidgets.QTableWidgetItem(elems[0])
                    table_item_lat = QtWidgets.QTableWidgetItem(elems[1])
                    table_item_long = QtWidgets.QTableWidgetItem(elems[2])
                    table_item_elev = QtWidgets.QTableWidgetItem(elems[3])
                    if idx >= self.table.rowCount():
                        self.table.insertRow(idx)
                    self.table.setItem(idx, 0, table_item_station)
                    self.table.setItem(idx, 1, table_item_lat)
                    self.table.setItem(idx, 2, table_item_long)
                    self.table.setItem(idx, 3, table_item_elev)
                # Delete excesss rows
                while idx < self.table.rowCount() - 1:
                    self.table.removeRow(self.table.rowCount())

    def coords(self):
        """
        Puts table coordinates into a dict.

        Returns
        -------
        dict
            key = station name, value = tuple (lat, long, elev)

        """
        c = dict()
        for i in range(self.table.rowCount()):
            c[self.table.item(i, 0).text()] = (
                float(self.table.item(i, 2).text()),
                float(self.table.item(i, 1).text()),
                float(self.table.item(i, 3).text()),
            )
        return c


class DialogLoopProperties(QtWidgets.QDialog):
    """
    Dialog for editing basic loop properties (operator, meter, comments)

    Parameters
    ----------
    loops : list
        Loops to apply settings, can be more than one

    """

    def __init__(self, loops, parent=None):
        super(DialogLoopProperties, self).__init__(parent)
        self.setGeometry(50, 50, 350, 350)
        self.loops = loops
        self.setWindowTitle("Loop properties")
        self.init_ui()

    def init_ui(self):
        # create buttons and actions
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.close)
        self.ok_button = QtWidgets.QPushButton("OK")
        self.ok_button.clicked.connect(self.set_loop_options)

        button_box = QtWidgets.QDialogButtonBox(QtCore.Qt.Horizontal)
        button_box.addButton(self.cancel_button, QtWidgets.QDialogButtonBox.ActionRole)
        button_box.addButton(self.ok_button, QtWidgets.QDialogButtonBox.ActionRole)
        # Sometimes these are undefined:
        try:
            self.operator_edit = QtWidgets.QLineEdit(self.loops[0].oper)
        except AttributeError:
            self.operator_edit = QtWidgets.QLineEdit("")
        try:
            self.meter_edit = QtWidgets.QLineEdit(self.loops[0].meter)
        except AttributeError:
            self.meter_edit = QtWidgets.QLineEdit("")
        self.comment_edit = QtWidgets.QTextEdit()
        try:
            self.comment_edit.setPlainText(self.loops[0].comment)
        except (IndexError, AttributeError):
            self.comment_edit.setPlainText("")
        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel("Operator"), 1, 0)
        grid.addWidget(self.operator_edit, 1, 1)
        grid.addWidget(QtWidgets.QLabel("Meter ID"), 2, 0)
        grid.addWidget(self.meter_edit, 2, 1)
        grid.addWidget(QtWidgets.QLabel("Comment"), 3, 0)
        grid.addWidget(self.comment_edit, 4, 0, 1, 2)
        grid.addWidget(button_box, 5, 0, 1, 2)
        self.setLayout(grid)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def set_loop_options(self):
        self.accept()


class AdjustOptions(QtWidgets.QDialog):
    """
    Dialog to set network adjustment options.

    Parameters
    ----------
    survey_str : str
        Used to display the survey name in the apply button
    options : AdjustOptions

    """

    def __init__(self, survey_str, options, parent=None):
        super(AdjustOptions, self).__init__(parent)
        self.setGeometry(50, 50, 350, 350)
        self.ao = options
        self.surveys_to_update = ""
        self.update_options = []
        self.sigma_prefactor_chk = QtWidgets.QCheckBox(
            "Standard deviation multiplier: pre-minimum"
        )
        self.sigma_postfactor_chk = QtWidgets.QCheckBox(
            "Standard deviation multiplier: post-minimum"
        )
        self.sigma_prefactor_chk.setToolTip(
            "This multiplier is applied prior to enforcing the minimum standard "
            "deviation (if checked)"
        )
        self.sigma_postfactor_chk.setToolTip(
            "This multiplier is applied after enforcing the minimum standard "
            "deviation (if checked)"
        )
        self.sigma_add_chk = QtWidgets.QCheckBox(
            "Add to standard deviation: pre-minimum"
        )
        self.sigma_min_chk = QtWidgets.QCheckBox("Minimum standard deviation")
        self.sigma_min_chk.stateChanged.connect(self.min_sd_checked_or_unchecked)
        self.cal_coeff_chk = QtWidgets.QCheckBox(
            "Calculate relative meter calibration coefficient"
        )
        self.cal_coeff_chk.stateChanged.connect(
            self.calc_cal_coeff_checked_or_unchecked
        )
        self.alpha_text = QtWidgets.QLabel("Significance level for global model test")
        self.cal_coeff_specify_chk = QtWidgets.QCheckBox(
            "Specify calibration coefficient"
        )
        self.cal_coeff_specify_chk.stateChanged.connect(
            self.specify_cal_coeff_checked_or_unchecked
        )

        self.cal_coeff_table = QtWidgets.QTableView()

        # self.drift_temp_edit = QtWidgets.QLineEdit(str(self.ao.model_temp_degree))
        self.sigma_prefactor_edit = QtWidgets.QLineEdit(
            "{:.4f}".format(self.ao.sigma_prefactor)
        )
        self.sigma_postfactor_edit = QtWidgets.QLineEdit(
            "{:.4f}".format(self.ao.sigma_postfactor)
        )
        self.sigma_add_edit = QtWidgets.QLineEdit(str(self.ao.sigma_add))
        self.sigma_min_edit = QtWidgets.QLineEdit(str(self.ao.sigma_min))
        self.alpha_edit = QtWidgets.QLineEdit(str(self.ao.alpha))
        self.init_ui(survey_str)

    def init_ui(self, survey_name):
        if survey_name is not None:
            self.sigma_prefactor_chk.setChecked(self.ao.use_sigma_prefactor)
            self.sigma_postfactor_chk.setChecked(self.ao.use_sigma_postfactor)
            self.sigma_add_chk.setChecked(self.ao.use_sigma_add)
            self.sigma_min_chk.setChecked(self.ao.use_sigma_min)
            self.sigma_postfactor_chk.setEnabled(self.ao.use_sigma_min)
            self.cal_coeff_chk.setChecked(self.ao.cal_coeff)
            try:
                self.cal_coeff_specify_chk.setChecked(self.ao.specify_cal_coeff)
            except Exception:
                self.ao.specify_cal_coeff = False
            if not hasattr(self.ao, "meter_cal_dict"):
                self.ao.meter_cal_dict = init_cal_coeff_dict(self.parent().obsTreeModel)
            elif self.ao.meter_cal_dict is None:
                self.ao.meter_cal_dict = init_cal_coeff_dict(self.parent().obsTreeModel)
            self.cal_coeff_model = MeterCalibrationModel()
            for k, v in self.ao.meter_cal_dict.items():
                self.cal_coeff_model.appendRow(
                    [QtGui.QStandardItem(k), QtGui.QStandardItem(str(v))]
                )

            self.cal_coeff_table.setModel(self.cal_coeff_model)
            self.cal_coeff_table.setFixedWidth(260)

            # create buttons and actions
            btn_restore_default = QtWidgets.QPushButton("Restore defaults")
            btn_restore_default.clicked.connect(self.restore_default)
            btn_cancel = QtWidgets.QPushButton("Cancel")
            btn_cancel.clicked.connect(self.close)
            btn_current = QtWidgets.QPushButton(
                "Apply to current survey (" + survey_name + ")"
            )
            btn_current.clicked.connect(self.apply_current)
            btn_all = QtWidgets.QPushButton("Apply to all surveys")
            btn_all.clicked.connect(self.apply_all)

            buttonBox = QtWidgets.QDialogButtonBox(QtCore.Qt.Horizontal)
            buttonBox.addButton(
                btn_restore_default, QtWidgets.QDialogButtonBox.ActionRole
            )
            buttonBox.addButton(btn_current, QtWidgets.QDialogButtonBox.ActionRole)
            buttonBox.addButton(btn_all, QtWidgets.QDialogButtonBox.ActionRole)
            buttonBox.addButton(btn_cancel, QtWidgets.QDialogButtonBox.ActionRole)

            vlayout = QtWidgets.QVBoxLayout()
            gridwidget = QtWidgets.QWidget()
            grid = QtWidgets.QGridLayout()
            vlayout.addWidget(gridwidget)
            # grid.addWidget(self.drift_temp_chk, 1, 0)
            # grid.addWidget(self.drift_temp_edit, 1, 1)
            grid.addWidget(self.sigma_prefactor_chk, 0, 0)
            grid.addWidget(self.sigma_prefactor_edit, 0, 1)
            grid.addWidget(self.sigma_postfactor_chk, 1, 0)
            grid.addWidget(self.sigma_postfactor_edit, 1, 1)
            grid.addWidget(self.sigma_add_chk, 2, 0)
            grid.addWidget(self.sigma_add_edit, 2, 1)
            grid.addWidget(self.sigma_min_chk, 3, 0)
            grid.addWidget(self.sigma_min_edit, 3, 1)
            grid.addWidget(self.cal_coeff_chk, 4, 0)
            grid.addWidget(self.cal_coeff_specify_chk, 5, 0)
            grid.addWidget(self.cal_coeff_table, 6, 0)
            grid.addWidget(self.alpha_text, 7, 0)
            grid.addWidget(self.alpha_edit, 7, 1)

            # This isn't working? Column 1 always seems to be stretching
            grid.setColumnStretch(0, 1)
            grid.setColumnStretch(1, 0)
            gridwidget.setLayout(grid)

            vlayout.addWidget(buttonBox)

            self.setLayout(vlayout)
            self.setWindowTitle("Network adjustment options")
            self.setWindowModality(QtCore.Qt.ApplicationModal)
        else:
            MessageBox.warning(
                "Network adjustment options",
                "Please load a survey first",
            )

    def restore_default(self):
        cbs = [
            self.sigma_add_chk,
            self.sigma_postfactor_chk,
            self.sigma_prefactor_chk,
            self.sigma_min_chk,
        ]
        edits = [
            self.sigma_add_edit,
            self.sigma_postfactor_edit,
            self.sigma_prefactor_edit,
            self.sigma_min_edit,
        ]
        values = ["0.0", "1.0", "1.0", "3.0"]
        for cb, ed, v in zip(cbs, edits, values):
            self.set_value_and_uncheck(cb, ed, v)
        self.cal_coeff_specify_chk.setChecked(False)
        self.cal_coeff_chk.setChecked(False)
        self.alpha_edit.setText("0.05")

    def set_value_and_uncheck(self, cb, ed, value="1.0"):
        cb.setChecked(False)
        ed.setText(value)

    def min_sd_checked_or_unchecked(self, state):
        if state == 2:  # checked
            self.sigma_postfactor_chk.setEnabled(True)
        else:
            self.sigma_postfactor_chk.setEnabled(False)

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
        try:
            self.ao.use_sigma_prefactor = self.sigma_prefactor_chk.isChecked()
            self.ao.use_sigma_postfactor = self.sigma_postfactor_chk.isChecked()
            self.ao.sigma_prefactor = float(self.sigma_prefactor_edit.text())
            self.ao.sigma_postfactor = float(self.sigma_postfactor_edit.text())
            self.ao.use_sigma_add = self.sigma_add_chk.isChecked()
            self.ao.sigma_add = float(self.sigma_add_edit.text())
            self.ao.use_sigma_min = self.sigma_min_chk.isChecked()
            self.ao.sigma_min = float(self.sigma_min_edit.text())
            self.ao.cal_coeff = self.cal_coeff_chk.isChecked()
            self.ao.specify_cal_coeff = self.cal_coeff_specify_chk.isChecked()
            if self.cal_coeff_specify_chk.isChecked():
                for i in range(self.cal_coeff_model.rowCount()):
                    meter = self.cal_coeff_model.itemFromIndex(
                        self.cal_coeff_model.index(i, 0)
                    )
                    calval = self.cal_coeff_model.itemFromIndex(
                        self.cal_coeff_model.index(i, 1)
                    )
                    self.ao.meter_cal_dict[meter.text()] = float(calval.text())
            self.ao.alpha = float(self.alpha_edit.text())

        except ValueError as e:  # caught if invalid number is entered in any box
            logging.error("Error setting adjustment options: %s", e)

    def apply_current(self):
        self.set_adjust_options()
        self.surveys_to_update = "single"
        self.accept()

    def apply_all(self):
        self.set_adjust_options()
        self.surveys_to_update = "all"
        self.accept()


class DialogMeterType(QtWidgets.QMessageBox):
    """
    Dialog to specify meter type when appending a survey or loop.
    """

    def __init__(self):
        super(DialogMeterType, self).__init__()
        self.setText("Choose meter file to import")
        self.addButton(
            QtWidgets.QPushButton(" CG-3, CG-5 "), QtWidgets.QMessageBox.YesRole
        )
        self.addButton(QtWidgets.QPushButton(" CG-6 "), QtWidgets.QMessageBox.YesRole)
        self.addButton(
            QtWidgets.QPushButton(" CG-6 Tsoft "), QtWidgets.QMessageBox.YesRole
        )
        self.addButton(QtWidgets.QPushButton(" Burris "), QtWidgets.QMessageBox.YesRole)
        self.addButton(QtWidgets.QPushButton(" CSV "), QtWidgets.QMessageBox.YesRole)
        self.addButton(
            QtWidgets.QPushButton(" Cancel "), QtWidgets.QMessageBox.RejectRole
        )
        self.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        meter = {
            " CG-3, CG-5 ": "CG5",
            " CG-6 ": "CG6",
            " CG-6 Tsoft ": "CG6Tsoft",
            " Burris ": "Burris",
            " CSV ": "csv",
            " Cancel ": "Cancel",
        }
        self.meter_type = meter[btn.text()]
        if self.meter_type == "Cancel":
            self.reject()
        else:
            self.accept()


class DialogOverwrite(QtWidgets.QMessageBox):
    """
    Dialog to confirm workspace overwrite
    """

    def __init__(self):
        super(DialogOverwrite, self).__init__()
        self.overwrite = False
        self.setText("Overwrite existing data?")
        self.setWindowTitle("GSadjust")
        ow_button = QtWidgets.QPushButton("Overwrite")
        ow_button.clicked.connect(self.onClicked)

        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)
        self.addButton(cancel_button, 0)
        self.addButton(ow_button, 1)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.accept()


class PlotGravityAndWLs(QtWidgets.QDialog):
    def __init__(self, MainProg):
        super(PlotGravityAndWLs, self).__init__(MainProg)
        self.obstreemodel = MainProg.obsTreeModel
        self.stations = MainProg.obsTreeModel.results_stations()
        self.coords = MainProg.obsTreeModel.station_coords
        self.setWindowModality(QtCore.Qt.NonModal)
        self.init_gui()

    def init_gui(self):
        self.setWindowTitle("Groundwater level/storage-change comparison")
        self.station_comboBox = QtWidgets.QComboBox(self)
        for s in sorted(self.stations):
            self.station_comboBox.addItem(s)
        self.station_comboBox.currentIndexChanged.connect(self.populateGravTable)
        layout = QtWidgets.QVBoxLayout()

        criteria_box = QtWidgets.QVBoxLayout()
        criteria_box_parent = QtWidgets.QWidget()

        interpolate_threshold_box = QtWidgets.QHBoxLayout()
        interpolate_threshold_box.addWidget(
            QtWidgets.QLabel("Interpolate gaps less than (days):")
        )
        self.threshold_spinbox = QtWidgets.QSpinBox()
        self.threshold_spinbox.setValue(30)
        self.threshold_spinbox.setMaximum(1000)
        self.threshold_spinbox.setMinimum(0)
        interpolate_threshold_box.addWidget(self.threshold_spinbox)
        criteria_box.addLayout(interpolate_threshold_box)

        gwl_timelag_box = QtWidgets.QHBoxLayout()
        gwl_timelag_box.addWidget(
            QtWidgets.QLabel("Time lag for groundwater levels (days):")
        )
        self.timelag_spinbox = QtWidgets.QSpinBox()
        self.timelag_spinbox.setValue(0)
        self.timelag_spinbox.setMinimum(-365)
        gwl_timelag_box.addWidget(self.timelag_spinbox)
        criteria_box.addLayout(gwl_timelag_box)

        find_optimal_timelag_box = QtWidgets.QHBoxLayout()
        find_optimal_timelag_box.addWidget(
            QtWidgets.QLabel("Find optimal time lag (maximize r<sup>2</sup>)")
        )
        self.cb_find_optimal = QtWidgets.QCheckBox()
        find_optimal_timelag_box.addWidget(self.cb_find_optimal)
        criteria_box.addLayout(find_optimal_timelag_box)

        filter_data_box = QtWidgets.QHBoxLayout()
        filter_data_box.addWidget(QtWidgets.QLabel("Filter groundwater-level data"))
        self.cb_filter_data = QtWidgets.QCheckBox()
        filter_data_box.addWidget(self.cb_filter_data)
        criteria_box.addLayout(filter_data_box)

        filter_edit_box = QtWidgets.QHBoxLayout()
        filter_edit_box.addWidget(QtWidgets.QLabel("Window (days):"))
        self.filter_spinbox = QtWidgets.QSpinBox()
        self.filter_spinbox.setMaximum(2000)
        self.filter_spinbox.setValue(700)
        filter_edit_box.addWidget(self.filter_spinbox)
        criteria_box.addLayout(filter_edit_box)

        divider = QtWidgets.QFrame()
        divider.setFrameShape(QtWidgets.QFrame.HLine)
        criteria_box.addWidget(divider)

        # Consider all water-level changes, only rising, or only falling
        bg = QtWidgets.QButtonGroup()
        self.button_all = QtWidgets.QRadioButton("All groundwater levels")
        self.button_rising = QtWidgets.QRadioButton("Rising groundwater level only")
        self.button_falling = QtWidgets.QRadioButton("Falling groundwater level only")
        bg.addButton(self.button_all)
        bg.addButton(self.button_rising)
        bg.addButton(self.button_falling)
        self.button_all.setChecked(True)

        button_box = QtWidgets.QVBoxLayout()
        button_box.addWidget(self.button_all)
        button_box.addWidget(self.button_rising)
        button_box.addWidget(self.button_falling)
        criteria_box.addLayout(button_box)
        criteria_box.addWidget(divider)

        g_station_box = QtWidgets.QHBoxLayout()
        g_station_box.addWidget(QtWidgets.QLabel("Gravity Station"))
        g_station_box.addWidget(self.station_comboBox)
        criteria_box.addLayout(g_station_box)

        criteria_box_parent.setLayout(criteria_box)
        layout.addWidget(criteria_box_parent)

        self.createGravTable()
        layout.addWidget(self.tableWidgetGrav)

        layout.addStretch()
        nwis_station_box = QtWidgets.QHBoxLayout()
        nwis_station_box_widget = QtWidgets.QWidget()
        nwis_station_box.addWidget(QtWidgets.QLabel("NWIS Station ID"))
        self.nwis_station_line_edit = QtWidgets.QLineEdit()
        btn_search_nwis = QtWidgets.QToolButton()
        btn_search_nwis.setIcon(QtGui.QIcon(":/icons/mag.png"))
        btn_search_nwis.clicked.connect(self.search_for_nwis_stations)
        nwis_station_box.addWidget(self.nwis_station_line_edit)
        nwis_station_box.addWidget(btn_search_nwis)
        nwis_station_box_widget.setLayout(nwis_station_box)
        layout.addWidget(nwis_station_box_widget)

        self.createNWISTable()
        layout.addWidget(self.tableWidgetNWIS)

        # Button box
        ok_button = QtWidgets.QPushButton("Ok")
        ok_button.clicked.connect(self.process)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        button_box = QtWidgets.QHBoxLayout()
        button_box.addStretch(1)
        button_box.addWidget(cancel_button)
        button_box.addWidget(ok_button)
        button_widget = QtWidgets.QWidget()
        button_widget.setLayout(button_box)

        layout.addWidget(button_widget)
        self.setLayout(layout)
        self.resize(750, 700)

        self.populateGravTable()

    def process(self):
        if self.nwis_station_line_edit.text() == "":
            MessageBox.warning(
                "Invalid NWIS station", "Please enter an NWIS station ID"
            )
            return
        for row_idx in range(self.tableWidgetNWIS.rowCount()):
            item = self.tableWidgetNWIS.item(row_idx, 0)
            if item.checkState() == 2:
                nwis_station = item.text()
                break
        else:
            nwis_station = self.nwis_station_line_edit.text()[:16]
        g_station = self.station_comboBox.currentText()
        dates, g = self.get_checked_data()
        # sd = [data[1][3][idx] for (idx, sta) in enumerate(data[1][0]) if
        #       sta == g_station]

        nwis_data = nwis_get_data(nwis_station)
        if self.cb_filter_data.isChecked():
            nwis_data = filter_ts(nwis_data, int(self.filter_spinbox.text()))
        if self.timelag_spinbox.text() != 0 and not self.cb_find_optimal.isChecked():
            t_offset = int(self.timelag_spinbox.text()) * -1
            nwis_data["continuous_x"] = [
                a + datetime.timedelta(t_offset) for a in nwis_data["continuous_x"]
            ]
            nwis_data["discrete_x"] = [
                a + datetime.timedelta(t_offset) for a in nwis_data["discrete_x"]
            ]
        else:
            t_offset = None
        if dates:
            try:
                if self.button_all.isChecked():
                    rise_or_fall = "all"
                elif self.button_rising.isChecked():
                    rise_or_fall = "rise"
                elif self.button_falling.isChecked():
                    rise_or_fall = "fall"
                optimize = self.cb_find_optimal.isChecked()
                plt = PlotNwis(
                    nwis_station,
                    g_station,
                    nwis_data,
                    (dates, g),
                    t_offset=t_offset,
                    optimize=optimize,
                    t_threshold=int(self.threshold_spinbox.text()),
                    rising=rise_or_fall,
                    parent=self,
                )
                plt.show()
            except IndexError as e:
                pass

    def get_checked_data(self):
        dates, g = [], []
        for r in range(self.tableWidgetGrav.rowCount()):
            item = self.tableWidgetGrav.item(r, 0)
            if item.checkState() == 2:
                d = self.tableWidgetGrav.item(r, 1)
                dates.append(datetime.datetime.strptime(d.text(), "%Y-%m-%d"))
                g_val = self.tableWidgetGrav.item(r, 2)
                g.append(float(g_val.text()))
        return dates, g

    def get_g_data(self, g_station):
        header, data, _ = compute_gravity_change(self.obstreemodel, table_type="list")
        dates = [
            datetime.datetime.strptime(date, "%Y-%m-%d")
            for station, date, g, s in data
            if station == g_station
        ]
        g = [g for station, date, g, sd in data if station == g_station]
        return data, dates, g

    # def plot_nwis(self, nwis_station, g_station, nwis, data):
    #     plt = PlotNwis(nwis_station, g_station, nwis, data, parent=self)
    #     return plt
    # return

    def search_for_nwis_stations(self):
        coords = self.coords[self.station_comboBox.currentText()]
        ID_dict = search_nwis(coords)
        if len(ID_dict) > 0:
            try:
                self.populateTable(ID_dict)
                nwis_string = next(iter(ID_dict.items()))[0]
                self.nwis_station_line_edit.setText(nwis_string)
            except:
                self.nwis_station_line_edit.setText("Error")
        else:
            self.nwis_station_line_edit.setText("No stations found")

    def populateTable(self, ID_dict):
        self.tableWidgetNWIS.setRowCount(len(ID_dict))
        row = 0
        for k, v in ID_dict.items():
            item = QtWidgets.QTableWidgetItem(k)
            item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            if row == 0:
                item.setCheckState(QtCore.Qt.Checked)
            else:
                item.setCheckState(QtCore.Qt.Unchecked)
            self.tableWidgetNWIS.setItem(row, 0, item)

            for idx, content in enumerate(v):
                item = QtWidgets.QTableWidgetItem(str(content))
                item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                self.tableWidgetNWIS.setItem(row, idx + 1, item)
            row += 1
        self.tableWidgetNWIS.cellChanged.connect(self.onCellChanged)
        self.tableWidgetNWIS.setSortingEnabled(True)

    def populateGravTable(self):
        self.tableWidgetGrav.clearContents()
        row = 0
        g_station = self.station_comboBox.currentText()
        data, dates, g = self.get_g_data(g_station)
        dg = [(x - g[0]) / 41.9 for x in g]
        self.tableWidgetGrav.setRowCount(len(dates))
        for idx, v in enumerate(dates):
            item = QtWidgets.QTableWidgetItem(g_station)
            item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            item.setCheckState(QtCore.Qt.Checked)
            self.tableWidgetGrav.setItem(row, 0, item)

            item = QtWidgets.QTableWidgetItem(str(v.strftime("%Y-%m-%d")))
            self.tableWidgetGrav.setItem(row, 1, item)
            item = QtWidgets.QTableWidgetItem(f"{g[idx]:0.1f}")
            self.tableWidgetGrav.setItem(row, 2, item)
            item = QtWidgets.QTableWidgetItem(f"{dg[idx]:0.2f}")
            self.tableWidgetGrav.setItem(row, 3, item)
            # item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            # self.tableWidgetNWIS.setItem(row, idx + 1, item)
            row += 1
        self.tableWidgetGrav.setSortingEnabled(True)

    # Create table
    def createNWISTable(self):
        self.tableWidgetNWIS = QtWidgets.QTableWidget()

        # Row count
        self.tableWidgetNWIS.setRowCount(0)

        # Column count
        self.tableWidgetNWIS.setColumnCount(6)
        # Table will fit the screen horizontally
        self.tableWidgetNWIS.horizontalHeader().setStretchLastSection(True)
        # self.tableWidget.horizontalHeader().setSectionResizeMode(
        #     QtWidgets.QHeaderView.Stretch)
        self.tableWidgetNWIS.setColumnWidth(0, 120)
        self.tableWidgetNWIS.setColumnWidth(1, 180)
        self.tableWidgetNWIS.setColumnWidth(2, 80)
        self.tableWidgetNWIS.setColumnWidth(3, 90)
        self.tableWidgetNWIS.setColumnWidth(4, 90)
        self.tableWidgetNWIS.setHorizontalHeaderLabels(
            [
                "Station",
                "Name",
                "Distance (m)",
                "First data",
                "Last data",
                "Well depth (ft)",
            ]
        )

    def createGravTable(self):
        self.tableWidgetGrav = QtWidgets.QTableWidget()

        # Row count
        self.tableWidgetGrav.setRowCount(10)

        # Column count
        self.tableWidgetGrav.setColumnCount(5)
        # Table will fit the screen horizontally
        self.tableWidgetGrav.horizontalHeader().setStretchLastSection(True)
        # self.tableWidget.horizontalHeader().setSectionResizeMode(
        #     QtWidgets.QHeaderView.Stretch)
        self.tableWidgetGrav.setColumnWidth(0, 100)
        self.tableWidgetGrav.setColumnWidth(1, 100)
        self.tableWidgetGrav.setColumnWidth(2, 100)
        self.tableWidgetGrav.setColumnWidth(3, 100)
        self.tableWidgetGrav.setHorizontalHeaderLabels(
            ["Station", "Date", "g", "Storage change", ""]
        )

    def onCellChanged(self, row, column):
        item = self.tableWidgetNWIS.item(row, column)
        # lastState = item.data(LastStateRole)

        if column != 0:
            return

        currentState = item.checkState()
        try:
            if currentState == QtCore.Qt.Checked:
                self.nwis_station_line_edit.setText(item.text())
                for row_idx in range(self.tableWidgetNWIS.rowCount()):
                    if row_idx != row:
                        temp_item = self.tableWidgetNWIS.item(row_idx, 0)
                        temp_item.setCheckState(QtCore.Qt.Unchecked)
        except:
            return


class GravityChangeTable(QtWidgets.QDialog):
    """
    Floating window to show gravity-change results

    Parameters
    ----------
    MainProg : MainProg object
    data : tuple
        header, table, dates - return value from compute_gravity_change
    table_type : {"simple", "full", "list"}
        controls how data are displayed

    """

    def __init__(self, MainProg, data, table_type="simple"):
        super(GravityChangeTable, self).__init__(MainProg)
        self.header, self.table, self.dates = data[0], data[1], data[2]
        self.full_header, self.full_table, _ = analysis.compute_gravity_change(
            MainProg.obsTreeModel, table_type="full"
        )
        self.list_header, self.list_table, _ = analysis.compute_gravity_change(
            MainProg.obsTreeModel, table_type="list"
        )
        self.coords = MainProg.obsTreeModel.station_coords
        layout = QtWidgets.QVBoxLayout()
        self.simple_model = GravityChangeModel(
            self.header, self.table, table_type="simple"
        )
        self.full_model = GravityChangeModel(
            self.full_header, self.full_table, table_type="full"
        )
        self.list_model = GravityChangeModel(
            self.list_header, self.list_table, table_type="list"
        )
        self.table_view = QtWidgets.QTableView()
        self.table_view.setModel(self.simple_model)

        # Create buttons and actions
        self.type_comboBox = QtWidgets.QComboBox(self)
        self.type_comboBox.addItem("simple dg")
        self.type_comboBox.addItem("full dg")
        self.type_comboBox.addItem("list")
        self.type_comboBox.activated.connect(self.table_changed)
        btn_map = QtWidgets.QPushButton("Map")
        btn_map.clicked.connect(self.map_change_window)
        btn_plot = QtWidgets.QPushButton("Plot")
        btn_plot.clicked.connect(
            lambda: MainProg.plot_gravity_change(self.list_table, self)
        )
        btn1 = QtWidgets.QPushButton("Copy to clipboard")
        btn1.clicked.connect(lambda: copy_cells_to_clipboard(self.table_view))

        # Locations
        layout.addWidget(self.table_view)
        hbox, hbox_widget = QtWidgets.QHBoxLayout(), QtWidgets.QWidget()
        hbox.addStretch(0)
        # Only show plot button on simple table
        hbox.addWidget(QtWidgets.QLabel("Table type:"))
        hbox.addWidget(self.type_comboBox)
        if self.dates:
            hbox.addWidget(btn_plot)
            hbox.addWidget(btn_map)
        hbox.addWidget(btn1)
        hbox_widget.setLayout(hbox)
        layout.addWidget(hbox_widget)
        self.setLayout(layout)
        self.setWindowTitle("Gravity Change")
        self.setGeometry(200, 200, 600, 800)

    def table_changed(self, index):
        if self.type_comboBox.itemText(index) == "simple dg":
            self.table_view.setModel(self.simple_model)
        if self.type_comboBox.itemText(index) == "full dg":
            self.table_view.setModel(self.full_model)
        if self.type_comboBox.itemText(index) == "list":
            self.table_view.setModel(self.list_model)

    def map_change_window(self):
        try:
            logging.info("try importing cartopy")
            # This SSL section needed to avoid certificate errors on Mac (?)
            import ssl

            if hasattr(ssl, "_create_unverified_context"):
                ssl._create_default_https_context = ssl._create_unverified_context
            import cartopy.crs as ccrs
            import cartopy.io.img_tiles as cimgt

            self.win = GravityChangeMap(
                self.table, self.header, self.coords, self.full_table, self.full_header
            )
            self.win.show()
        except (OSError, ImportError) as e:
            logging.info("Cartopy import error")
            MessageBox.warning(
                "Import error",
                "Map view plots on Mac or Linux platforms requires installation of the"
                ' Geos and Proj libraries. Please install with homebrew ("brew install'
                ' geos proj").',
            )


class VerticalGradientDialog(QtWidgets.QInputDialog):
    def __init__(self, default_interval):
        super(VerticalGradientDialog, self).__init__()
        self.text, self.ok = self.getDouble(
            None,
            "Vertical-gradient interval",
            "Interval, in cm:",
            default_interval,
            0,
            200,
            3,
        )


class TideCoordinatesDialog(QtWidgets.QDialog):
    """
    Get coordinates for tide correction.
    """

    def __init__(self, lat, lon, elev):
        super(TideCoordinatesDialog, self).__init__()
        self.init_ui(lat, lon, elev)

    def init_ui(self, lat, lon, elev):
        latLabel = QtWidgets.QLabel("Latitude")
        lonLabel = QtWidgets.QLabel("Longitude")
        elevLabel = QtWidgets.QLabel("Elevation")
        self.lat = QtWidgets.QLineEdit()
        self.lon = QtWidgets.QLineEdit()
        self.elev = QtWidgets.QLineEdit()

        self.lat.setText("{:.5f}".format(lat))
        self.lon.setText("{:.5f}".format(lon))
        self.elev.setText("{:.1f}".format(elev))

        # create buttons and actions
        self.btn1 = QtWidgets.QPushButton("OK")
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
        self.setWindowTitle("Survey coordinates")
        # enterCoordinates.show()
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.accept()


class TideCorrectionDialog(QtWidgets.QDialog):
    """
    Get method for tide correction.
    """

    def __init__(self):
        super(TideCorrectionDialog, self).__init__()
        self.init_ui()

    def init_ui(self):
        self.msg = QtWidgets.QMessageBox()
        self.msg.setText("Select tide-correction method")
        self.msg.addButton(
            QtWidgets.QPushButton("Meter-supplied"), QtWidgets.QMessageBox.YesRole
        )
        self.msg.addButton(
            QtWidgets.QPushButton("Agnew"), QtWidgets.QMessageBox.YesRole
        )
        btn2 = QtWidgets.QPushButton("Cancel")
        self.msg.addButton(btn2, QtWidgets.QMessageBox.RejectRole)
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.correction_type = btn.text()
        if self.correction_type == "Cancel":
            self.reject()
        else:
            self.accept()


class AddDatumFromList(QtWidgets.QInputDialog):
    """
    Dialog to add existing station as a datum station.
    """

    @classmethod
    def add_datum(cls, station_list):
        dialog = cls()
        text, ok = dialog.getItem(None, "Input Dialog", "Datum station:", station_list)
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
        self.title = "Specify elapsed time that indicates a new loop"
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
        self.title = "Specify date, time, and magnitude of tare"
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

    Parameters
    ----------
    path : str
        Path to open. Should be same as previous, if dialog was previously opened.
    datum_table_model : list
        If not None, it's the datums the previous time this dialog was opened. Saves
        time for large directories.

    """

    def __init__(self, path, datum_table_model=None, parent=None):
        super(SelectAbsg, self).__init__(parent)
        self.title = "Select directory with Abs. g files (.project.txt)"
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 800
        self.new_datums = []
        self.path = path

        self.splitter_window = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self)
        self.tree_model = QtWidgets.QDirModel()
        self.tree = QtWidgets.QTreeView()
        self.tree_model.setFilter(QtCore.QDir.Dirs | QtCore.QDir.NoDotAndDotDot)
        self.table = QtWidgets.QTableView()

        self.ProxyModel = QtCore.QSortFilterProxyModel()

        self.table.setModel(self.ProxyModel)
        self.table.setSortingEnabled(True)

        if datum_table_model is not None:
            self.table_model = datum_table_model
            for i in range(self.table_model.rowCount()):
                idx = self.table_model.index(i, 0)
                self.table_model.setData(
                    idx, QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole
                )
        else:
            self.table_model = DatumTableModel()
        self.ProxyModel.setSourceModel(self.table_model)
        delegate = CheckBoxDelegate(None)
        self.table.setItemDelegateForColumn(0, delegate)
        QtWidgets.QApplication.processEvents()
        self.init_ui()

    def export_and_close(self):
        for i in range(self.ProxyModel.rowCount()):
            ndi = self.ProxyModel.index(i, 0)
            nd = self.ProxyModel.data(ndi, role=QtCore.Qt.UserRole)
            chk = self.ProxyModel.data(ndi, role=QtCore.Qt.CheckStateRole)
            if chk == 2 or chk == 1:
                self.new_datums.append(nd)
        self.accept()

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
        self.load_unpublished_cb = QtWidgets.QCheckBox("Ignore unpublished")
        self.load_unpublished_cb.setToolTip(
            'Files inside a directory named "unpublished" (anywhere on the path) will'
            " be ignored"
        )
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

        # self.ProxyModel.setSourceModel(self.table_model)
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
        Parses *.project.txt files in the selected paths. Populates the dialog
        table model directly.
        """
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        idxs = self.tree.selectedIndexes()
        self.table_model.clearDatums()
        files_found = None
        for i in idxs:
            if i.model().isDir(i):
                path = str(i.model().filePath(i))
                files_found = self.append_datums(path)
        QtWidgets.QApplication.restoreOverrideCursor()
        if not files_found:
            MessageBox.warning(
                "Import error",
                "No *.project.txt files found in the selected directories.",
            )

    def append_datums(self, path):
        files_found = False
        n_files = 100
        pbar = ProgressBar(total=n_files, textmess="Scanning directory...")
        ctr = 1
        pbar.show()
        for dirname, _, fileList in os.walk(path):
            if (
                self.load_unpublished_cb.isChecked()
                and dirname.find("unpublished") != -1
            ):
                continue
            for name in fileList:
                pbar.progressbar.setValue(ctr)
                QtWidgets.QApplication.processEvents()

                if ".project.txt" not in name:
                    continue
                ctr += 1
                if ctr > 100:
                    ctr = 1
                files_found = True
                d = a10.A10(os.path.join(dirname, name))
                datum = Datum(
                    d.stationname,
                    g=float(d.gravity),
                    sd=float(d.setscatter),
                    date=d.date,
                    meas_height=float(d.transferht),
                    gradient=float(d.gradient),
                    checked=0,
                )
                datum.coord = [d.long, d.lat, d.elev]
                datum.n_sets = d.processed
                datum.time = d.time
                self.table_model.insertRows(datum, 0)
                self.path = path
        return files_found


class DialogUpdateDatums(QtWidgets.QMessageBox):
    """
    Dialog to confirm workspace overwrite
    """

    def __init__(self, dlg_text):
        super(DialogUpdateDatums, self).__init__()
        self.overwrite = False
        self.setText("\n".join(dlg_text))
        self.setWindowTitle("GSadjust")

        ok_button = QtWidgets.QPushButton("Update")
        ok_button.clicked.connect(self.onClicked)

        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)
        self.addButton(cancel_button, 0)
        self.addButton(ok_button, 1)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        self.accept()


class AboutDialog(QtWidgets.QDialog):
    def __init__(self, version):
        super(AboutDialog, self).__init__()

        msg1 = (
            "<html>GSadjust, a product of the USGS Southwest Gravity Program<br>"
            + '<a href ="http://go.usa.gov/xqBnQ">http://go.usa.gov/xqBnQ</a>'
            + "<br><br>Commit "
            + version
            + '<br><br><a href ="https://github.com/jkennedy-usgs/sgp-gsadjust">'
            + "https://github.com/jkennedy-usgs/sgp-gsadjust</a>"
            + '<br><a href="mailto:jkennedy@usgs.gov">jkennedy@usgs.gov</a>'
        )
        QtWidgets.QMessageBox.about(None, "GSadust", msg1)


class ShowCalCoeffs(QtWidgets.QDialog):
    """
    Dialog to show table of calibration coefficients (specified or calculated.

    Parameters
    ----------
    cal_coeffs : dict
        key = meter name (str), value = list of tuples: (survey name, cal coeff, s.d.)

    """

    def __init__(self, cal_coeffs, parent=None):
        super(ShowCalCoeffs, self).__init__(parent)
        self.setWindowTitle("Calibration coefficients")
        # self.resize(20, 20)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )
        vlayout = QtWidgets.QVBoxLayout()
        for meter, data in cal_coeffs.items():
            h = 0
            vlayout.addWidget(QtWidgets.QLabel(meter))
            view = QtWidgets.QTableView()
            vlayout.addWidget(view)
            cal_model = QtGui.QStandardItemModel()
            cal_model.setColumnCount(3)
            cal_model.setHorizontalHeaderLabels(["Date", "Coeff.", "S.D."])
            for row in data:
                cal_model.appendRow(
                    [
                        QtGui.QStandardItem(row[0]),
                        QtGui.QStandardItem(f"{row[1]:.6f}"),
                        QtGui.QStandardItem(f"{row[2]:.6f}"),
                    ]
                )
                h += 30
            view.setModel(cal_model)
            view.setFixedSize(340, h + 40)
        self.setLayout(vlayout)
        self.resize(self.sizeHint())


class PreferencesDialog(QtWidgets.QDialog):
    """
    Dialog to show all preferences
    """

    def __init__(self, settings):
        super(PreferencesDialog, self).__init__()
        self.left = 100
        self.top = 100
        self.width = 400
        self.height = 200
        self.settings = settings
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("GSadjust Preferences")
        self.setGeometry(self.left, self.top, self.width, self.height)

        layout = QtWidgets.QVBoxLayout()

        include_datums = QtWidgets.QHBoxLayout()
        include_datums.addWidget(
            QtWidgets.QLabel("Include non-network datums in results")
        )
        include_datums.addStretch(1)
        self.cb_include_datums = QtWidgets.QCheckBox()
        if self.settings.value("include_all_datums_in_results"):
            self.cb_include_datums.setChecked(True)
        else:
            self.cb_include_datums.setChecked(False)
        include_datums.addWidget(self.cb_include_datums)
        layout.addLayout(include_datums)
        layout.addStretch(1)

        # Button box
        ok_button = QtWidgets.QPushButton("Ok")
        ok_button.clicked.connect(self.set_settings)
        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        button_box = QtWidgets.QHBoxLayout()
        button_box.addStretch(1)
        button_box.addWidget(cancel_button)
        button_box.addWidget(ok_button)
        button_widget = QtWidgets.QWidget()
        button_widget.setLayout(button_box)

        layout.addWidget(button_widget)
        self.setLayout(layout)

    def set_settings(self):
        if self.cb_include_datums.isChecked():
            self.settings.setValue("include_all_datums_in_results", True)
        else:
            self.settings.setValue("include_all_datums_in_results", False)
        self.close()


class DialogApplyTimeCorrection(QtWidgets.QDialog):
    """
    Dialog to apply time shift to observed times
    """

    def __init__(self):
        super(DialogApplyTimeCorrection, self).__init__()
        self.time_correction_type = False
        self.msg = QtWidgets.QMessageBox()
        self.msg.setWindowTitle("GSadjust")
        self.msg.setText("Apply time correction to:")
        self.msg.addButton(
            QtWidgets.QPushButton("Current station(s)"), QtWidgets.QMessageBox.YesRole
        )
        self.msg.addButton(
            QtWidgets.QPushButton("Current loop"), QtWidgets.QMessageBox.YesRole
        )
        self.msg.addButton(
            QtWidgets.QPushButton("Current survey"), QtWidgets.QMessageBox.YesRole
        )
        self.msg.addButton(
            QtWidgets.QPushButton("All data"), QtWidgets.QMessageBox.YesRole
        )
        self.msg.addButton(
            QtWidgets.QPushButton("Cancel"), QtWidgets.QMessageBox.RejectRole
        )
        self.msg.buttonClicked.connect(self.onClicked)
        self.setWindowModality(QtCore.Qt.ApplicationModal)

    def onClicked(self, btn):
        if btn.text() == "Cancel":
            self.reject()
        else:
            options = {
                "Current station(s)": "station",
                "Current loop": "loop",
                "Current survey": "survey",
                "All data": "all",
            }
            self.time_correction_type = options[btn.text()]
            self.accept()
