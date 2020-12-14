#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSadjust, for the interactive network adjustment of relative-gravity networks
=============================================================================
Jeff Kennedy, USGS
jkennedy@usgs.gov

GSadjust is a Python/PyQt5 program that provides a graphical user interface
for post-processing combined relative- and absolute- gravity surveys.

GSadjust was originally based on PyGrav:
Hector, B. and Hinderer, J.: pyGrav, a Python-based program for handling and
processing relative gravity data, Computers & Geosciences
doi:10.1016/j.cageo.2016.03.010, 2016.

GSadjust is distributed with the network adjustment software Gravnet
(Windows executable) (Hwang, C., C. Wang and L. Lee. 2002. Adjustment of
relative gravity measurements using weighted and datum-free constraints.
Computers & Geosciences 28 1005-1015)

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.

GSadjust includes the following submodules:

data
  adjustment
  analysis
  channel
  correction
  datum
  delta
  tare
drift
  continuous
  roman
file
  a10
  read
  write
gui
  dialogs
  logger
  menus
  messages
  utils
  widgets
  tabs
    data
    drift
    network
models
  burris
  calibration
  datum
  delta
  delta_proxy
  gravity
  results
  roman
  scintrex
  sort
  tare
  utils
obstree
  base
  loop
  model
  station
  survey
plots
  datum
  gravity
  loop
  network
  residual
resources
tides
  correction
  synthetic

PyQt models follow the PyQt CamelCase naming convention. The other
methods/functions in GSadjust use PEP-8 lowercase_underscore convention.

The general data structure of the program is stored as PyQt objects inherited
from QStandardItem:

+------------+
| Campaign   |
++-----------+
 |  +---------------------------+
 +--+ Surveys (ObsTreeSurvey)   |
    ++--------------------------+
     |  +-------------+
     +--+ Adjustment  |
     |  +-------------+
     |  +-------------+
     +--+ datums      |
     |  +-------------+
     |  +-------------+
     +--+ deltas      |
     |  +-------------+
     |  +---------------+
     +--+ results       |
     |  +---------------+
     |  +----------------------+
     |--+ Loops (ObsTreeLoop)  |
        ++---------------------+
         |  +-------------+
         +--+ deltas      |
            +-------------+
         |  +-------------+
         +--+ tares       |
            +-------------+
         |  +----------------------------+
         +--+ Stations (ObsTreeStation)  |
            ++---------------------------+
             |  +--------------------+
             +--+ g, lat, long, etc. |
                +--------------------+

The Adjustment object in each ObsTreeSurvey holds the network adjustment input,
output, and options. There is one Adjustment per Survey.

"""

import copy
import logging
import os

# Standard library modules
import sys
import time
import traceback
import webbrowser

# Modules that must be installed
import matplotlib
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QSettings, Qt
from matplotlib.dates import num2date, date2num

from . import resources
from .data import (
    ChannelList,
    Datum,
    Tare,
)
from .data.analysis import compute_gravity_change
from .data.correction import time_correction
from .file import (
    InvalidMeterException,
    export_data,
    export_metadata,
    export_summary,
    file_reader,
    import_abs_g_complete,
    import_abs_g_simple,
)
from .gui.dialogs import (
    AboutDialog,
    AddDatumFromList,
    AddTareDialog,
    AdjustOptions,
    CoordinatesTable,
    DialogApplyTimeCorrection,
    DialogLoopProperties,
    DialogMeterType,
    DialogOverwrite,
    GravityChangeTable,
    LoopTimeThresholdDialog,
    SelectAbsg,
    ShowCalCoeffs,
    TideCoordinatesDialog,
    TideCorrectionDialog,
    VerticalGradientDialog,
)
from .gui.logger import LoggerWidget
from .gui.menus import MENU_STATE, Menus
from .gui.messages import MessageBox
from .gui.tabs import TabAdjust, TabData, TabDrift
from .gui.widgets import ProgressBar
from .models import (
    BurrisTableModel,
    DatumTableModel,
    DeltaTableModel,
    ResultsTableModel,
    ScintrexTableModel,
)
from .obstree import ObsTreeLoop, ObsTreeModel, ObsTreeStation, ObsTreeSurvey
from .plots import (
    PlotDatumCompare,
    PlotDatumComparisonTimeSeries,
    PlotDgResidualHistogram,
    PlotGravityChange,
    PlotLoopAnimation,
    PlotNetworkGraph,
)
from .tides import tide_correction_agnew, tide_correction_meter
from .utils import (
    init_station_coords_dict,
)

matplotlib.use("qt5agg")

DOWN = 1
UP = -1
DIRECTIONS = (UP, DOWN)


class MainProg(QtWidgets.QMainWindow):
    """
    GSadjust main routine
    """

    path_output = None
    drift_lookup = {"none": 0, "netadj": 1, "roman": 2, "continuous": 3}
    obsTreeModel = ObsTreeModel()
    previous_loop = None
    previous_survey = None
    vertical_gradient_interval = 64.2
    workspace_loaded = False
    DPI = 60  # resolution for data plots
    all_survey_data = None

    # PyQt indexes
    index_current_survey = None
    index_current_loop = None
    index_current_station = None
    index_current_loop_survey = None
    index_current_station_loop = None
    index_current_station_survey = None

    label_adjust_update_required_set = False

    gui_obstreeview_popup_menu = None
    menus, selection_model = None, None
    tab_data, tab_drift, tab_adjust = None, None, None
    station_model = None
    workspace_savename = None

    def __init__(self, splash=None):
        super(MainProg, self).__init__()

        self.settings = QSettings("SGP", "GSADJUST")
        self.init_settings()
        self.settings.sync()

        self.menus = Menus(self)
        self.setWindowTitle("GSadjust")

        # Set up logging to use custom widget (no parent, so window).
        self.logview = LoggerWidget()
        logging.getLogger().setLevel(logging.INFO)
        logging.getLogger().addHandler(self.logview.handler)
        logging.info("Logger initialized.")

        # Create model objects for main views.
        self.delta_model = DeltaTableModel()
        self.datum_model = DatumTableModel()
        self.results_model = ResultsTableModel()

        # tab_....s are populated with GUI elements in the tab_...() functions
        self.tab_data = TabData(self)
        self.tab_drift = TabDrift(self)
        self.tab_adjust = TabAdjust(self)
        self.tab_widget = QtWidgets.QTabWidget()

        # Connect signals.
        self.delta_model.signal_adjust_update_required.connect(
            self.adjust_update_required
        )
        self.delta_model.tried_to_update_list_delta.connect(
            self.show_delta_update_message
        )
        self.datum_model.signal_adjust_update_required.connect(
            self.adjust_update_required
        )
        self.tab_drift.drift_plot_weighted.update_drift_plots.connect(
            self.update_drift_tables_and_plots
        )
        self.tab_drift.drift_screen_elapsed_time.update_drift_plots.connect(
            self.update_drift_tables_and_plots
        )
        self.tab_drift.tare_view.model().dataChanged.connect(self.tares_updated)

        self.tab_widget = QtWidgets.QTabWidget()
        self.tab_widget.addTab(self.tab_data, "Data")
        self.tab_widget.addTab(self.tab_drift, "Drift")
        self.tab_widget.addTab(self.tab_adjust, "Network Adjustment")
        self.tab_widget.setStyleSheet(
            "QTabBar::tab { font-size: 12pt; height: 40px; width: 250px }"
        )

        self.gui_treeview_widget = QtWidgets.QWidget()
        self.gui_treeview_box = QtWidgets.QVBoxLayout()
        self.gui_data_treeview = QtWidgets.QTreeView()

        self.qtaction_move_station_up = QtWidgets.QAction(
            QtGui.QIcon(":/icons/up.png"), "Move survey up", self
        )
        self.qtaction_move_station_up.triggered.connect(
            slot=lambda: self.move_survey(UP)
        )
        self.qtaction_move_station_down = QtWidgets.QAction(
            QtGui.QIcon(":/icons/down.png"), "Move survey down", self
        )
        self.qtaction_move_station_down.triggered.connect(
            slot=lambda: self.move_survey(DOWN)
        )
        self.gui_toolbar = QtWidgets.QToolBar()
        self.gui_toolbar.addAction(self.qtaction_move_station_up)
        self.gui_toolbar.addAction(self.qtaction_move_station_down)

        self.qtaction_collapse_all = QtWidgets.QAction(
            QtGui.QIcon(":/icons/ca.png"), "Collapse all", self
        )
        self.qtaction_collapse_all.triggered.connect(
            slot=self.gui_data_treeview.collapseAll
        )
        self.qtaction_expand_all = QtWidgets.QAction(
            QtGui.QIcon(":/icons/ea.png"), "Expand all", self
        )
        self.qtaction_expand_all.triggered.connect(
            slot=self.gui_data_treeview.expandAll
        )

        self.gui_toolbar.addAction(self.qtaction_collapse_all)
        self.gui_toolbar.addAction(self.qtaction_expand_all)
        self.gui_toolbar.addAction(self.menus.mnAdjAdjustCurrent)
        self.gui_toolbar.addAction(self.menus.mnAdjAdjust)
        self.gui_toolbar.addAction(self.menus.mnAdjUpdateSD)

        self.gui_treeview_box.addWidget(self.gui_toolbar)
        self.gui_treeview_box.addWidget(self.gui_data_treeview)
        self.gui_treeview_widget.setLayout(self.gui_treeview_box)

        self.gui_main_window_splitter = QtWidgets.QSplitter(Qt.Horizontal)
        self.gui_main_window_splitter.addWidget(self.gui_treeview_widget)
        self.gui_main_window_splitter.addWidget(self.tab_widget)
        self.setCentralWidget(self.gui_main_window_splitter)

        # Setup statusbar icons
        self.update_adjust_text = QtWidgets.QLabel("    Update adjustment:", self)
        self.update_not_needed_icon = QtGui.QPixmap(":/icons/ico3.png")
        self.update_adjust_icon = QtGui.QPixmap(":/icons/ico2.png")
        self.label_adjust_update_required = QtWidgets.QLabel()

        self.statusBar().addPermanentWidget(self.update_adjust_text)
        self.statusBar().addPermanentWidget(self.label_adjust_update_required)

        # Right-click tree view context menu
        self.menus.mnDeleteSurvey = self.menus.create_action(
            "Delete survey", slot=self.delete_survey
        )
        self.menus.mnDeleteLoop = self.menus.create_action(
            "Delete loop", slot=self.delete_loop
        )
        self.menus.mnLoopProperties = self.menus.create_action(
            "Loop properties...", slot=self.properties_loop
        )
        self.menus.mnRename = self.menus.create_action("Rename", slot=self.rename)
        self.menus.mnDeleteStation = self.menus.create_action(
            "Delete station(s)", slot=self.delete_station
        )
        self.menus.mnStationDuplicate = self.menus.create_action(
            "Duplicate station", slot=self.duplicate_station
        )
        self.menus.mnDataNewLoop = self.menus.create_action(
            "Move stations to new loop", slot=self.new_loop
        )
        self.menus.mnVerticalGradientWriteAction = self.menus.create_action(
            "Write vertical gradient file...", slot=self.vertical_gradient_write
        )
        self.menus.mnLoopAnimate = self.menus.create_action(
            "Animate loop", slot=self.animate_loop
        )

        self.update_menus()
        self.path_install = os.path.abspath(os.path.join(os.path.dirname(__file__),'..')) # os.getcwd()
        self.menus.set_state(MENU_STATE.UNINIT)

    def init_gui(self):
        """
        Called after loading a data file.
        """

        # Left panel: tree with data hierarchy (surveys, loops, stations)
        self.obsTreeModel.setHorizontalHeaderLabels(["Name", "Date", "g (\u00b5Gal)"])

        # Enable menus
        self.menus.set_state(MENU_STATE.INIT)

        # Resize, expand tree view
        self.gui_data_treeview.setModel(self.obsTreeModel)
        self.obsTreeModel.dataChanged.connect(self.on_obs_checked_change)

        self.selection_model = self.gui_data_treeview.selectionModel()
        self.selection_model.selectionChanged.connect(self.on_obs_tree_change)

        self.gui_data_treeview.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        self.gui_data_treeview.setContextMenuPolicy(Qt.CustomContextMenu)
        self.gui_data_treeview.customContextMenuRequested.connect(
            self.treeview_context_menu
        )
        self.gui_data_treeview.setItemDelegate(BoldDelegate(self))
        self.gui_data_treeview.doubleClicked.connect(self.activate_survey_or_loop)
        self.gui_data_treeview.setObjectName("data")
        self.gui_data_treeview.setEditTriggers(QtWidgets.QTreeView.EditKeyPressed)
        self.gui_data_treeview.setExpandsOnDoubleClick(False)
        self.gui_data_treeview.expandAll()
        self.gui_data_treeview.resizeColumnToContents(0)
        self.gui_data_treeview.resizeColumnToContents(1)
        self.gui_data_treeview.resizeColumnToContents(2)

        # Highlight first tree-view item
        self.select_first_treeview_item()

        # Set models for the tab views.
        self.tab_adjust.delta_proxy_model.setSourceModel(self.delta_model)
        self.tab_adjust.datum_proxy_model.setSourceModel(self.datum_model)
        self.tab_adjust.results_proxy_model.setSourceModel(self.results_model)

        # Activate first tree view item
        self.activate_survey_or_loop(self.index_current_loop)
        self.activate_survey_or_loop(self.index_current_survey)
        self.label_adjust_update_required_set = False
        # Set data plot
        self.update_data_tab()
        self.gui_data_treeview.setFocus()

    def toggle_logview(self):
        if self.logview.isVisible():
            self.logview.hide()
        else:
            self.logview.show()

    def adjust_update_required(self):
        """
        Updates status bar icon.
        """
        self.label_adjust_update_required_set = True
        self.label_adjust_update_required.setPixmap(self.update_adjust_icon)
        self.label_adjust_update_required.setToolTip("Update network adjustment")
        self.set_window_title_asterisk()
        self.update_menus()

    def adjust_update_not_required(self):
        """
        Clears status bar icon.
        """
        self.label_adjust_update_required_set = False
        self.label_adjust_update_required.setPixmap(self.update_not_needed_icon)
        self.label_adjust_update_required.setToolTip("Network adjustment is up to date")

    def select_first_treeview_item(self):
        """
        Selects the first item in the treeview
        """
        obstreesurvey = self.obsTreeModel.itemFromIndex(self.obsTreeModel.index(0, 0))
        obstreeloop = obstreesurvey.child(0)
        station = obstreeloop.child(0)

        self.index_current_survey = obstreesurvey.index()
        self.index_current_loop = obstreeloop.index()
        self.index_current_loop_survey = obstreesurvey.index()
        self.index_current_station_loop = obstreeloop.index()
        self.index_current_station_survey = obstreesurvey.index()
        self.index_current_station = station.index()

    def update_data_tab(self):
        """
        Get station to plot, update station table model if necessary.
        """
        obstreestation = self.obsTreeModel.itemFromIndex(self.index_current_station)
        obstreeloop = obstreestation.parent()
        station = obstreestation
        if obstreeloop.meter_type in ["CG5", "Scintrex", "CG6", "csv", "CG6Tsoft"]:
            self.station_model = ScintrexTableModel(station)
        elif obstreeloop.meter_type == "Burris":
            self.station_model = BurrisTableModel(station)
        self.station_model.dataChanged.connect(self.update_data_tab)
        self.station_model.signal_update_coordinates.connect(
            self.populate_station_coords
        )
        self.station_model.signal_adjust_update_required.connect(
            self.adjust_update_required
        )
        self.station_model.signal_uncheck_station.connect(self.uncheck_station)
        self.station_model.signal_check_station.connect(self.check_station)
        self.tab_data.data_view.setModel(self.station_model)
        self.index_current_station_loop = obstreeloop.index()
        self.index_current_station_survey = obstreeloop.parent().index()
        self.obsTreeModel.dataChanged.emit(QtCore.QModelIndex(), QtCore.QModelIndex())
        self.tab_data.update_station_plot(station, obstreeloop.meter_type)

    def uncheck_station(self):
        """
        Unchecks station in the tree view.

        Called when all samples at a station are unchecked; in that case we want to exclude
        that station from the deltas.

        Returns
        -------
        None
        """
        obstreestation = self.obsTreeModel.itemFromIndex(self.index_current_station)
        obstreestation.setCheckState(0)

    def check_station(self):
        """
        Checks station in the tree view.

        Called when one or more samples at a station are unchecked; in that case we want to exclude
        that station from the deltas.

        Returns
        -------
        None
        """
        obstreestation = self.obsTreeModel.itemFromIndex(self.index_current_station)
        obstreestation.setCheckState(2)

    ###########################################################################
    # Load/Open/Save routines
    ###########################################################################
    def open_file_dialog(self, open_type):
        """
        Called from menus.py functions
        Parameters
        ----------
        open_type

        Returns
        -------

        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            caption="Open file",
            directory=self.settings.value("current_dir"),
            filter="Data files (*.csv *.txt)",
        )

        if fname:
            self.settings.setValue("current_dir", os.path.dirname(fname))
            if self.open_raw_data(fname, open_type):
                self.init_gui()

    def open_raw_data(self, fname, open_type):
        """
        - Display a file opening window
        - Populate obsTreeModel

        Parameters
        ----------
        fname : str
            Filename to open
        open_type : {'choose', 'loop', meter type}
            'choose' - Choose meter-style data format
            'loop' - if appending loop to survey, otherwise assume appending
            survey to campaign (can be both choose and loop, e.g. 'chooseloop')
            'CG-6', 'Burris', or 'CG-5' - reading a raw data file, no appending
        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        # open file
        append_loop = True if open_type == "loop" else False

        # When "append survey' or 'append loop' is called: accommodate the rare
        # instance of combining meter types on a single survey
        if open_type in ["loop", "survey"]:
            QtWidgets.QApplication.restoreOverrideCursor()
            meter_type_dialog = DialogMeterType()
            test = meter_type_dialog.exec_()
            if test < 5:  # 5 = cancel  (accept/reject not working?)
                meter_type = meter_type_dialog.meter_type
            else:
                return False
        else:
            meter_type = open_type
            if self.obsTreeModel.invisibleRootItem().rowCount() > 0:
                overwrite_tree_dialog = DialogOverwrite()
                if overwrite_tree_dialog.exec_():
                    self.workspace_clear(confirm=False)
                else:
                    return False

        if fname:
            self.path_output = os.path.dirname(fname)
            logging.info("Loading data file: %s", fname)

            # populate a Campaign object
            e = None
            try:
                # TODO: this overwrites all_survey_data, is it used anywhere else?
                self.all_survey_data = self.read_raw_data_file(fname, meter_type)
                if append_loop:
                    obstreesurvey = self.obsTreeModel.itemFromIndex(
                        self.index_current_survey
                    )
                    # Loads all survey data into a single loop.
                    new_loop_idx = obstreesurvey.populate(
                        self.all_survey_data,
                        name=str(obstreesurvey.loop_count),
                        source=os.path.basename(fname),
                    )
                else:
                    obstreesurvey = ObsTreeSurvey(
                        str(num2date(self.all_survey_data.t[0]).date())
                    )
                    new_loop_idx = obstreesurvey.populate(
                        self.all_survey_data, source=os.path.basename(fname)
                    )
                    self.obsTreeModel.appendRow(
                        [
                            obstreesurvey,
                            QtGui.QStandardItem("a"),
                            QtGui.QStandardItem("a"),
                        ]
                    )
                self.activate_survey_or_loop(new_loop_idx)
            except IOError as err:
                e = err
                MessageBox.warning("File error", "No file : {}".format(fname))

            except (IndexError, ValueError) as err:
                stream = QtCore.QFile(":/text/err_{}.txt".format(meter_type))
                stream.open(QtCore.QIODevice.ReadOnly)
                text = QtCore.QTextStream(stream).readAll()
                stream.close()
                if hasattr(err, "i") and hasattr(err, "line"):
                    help_message = text.format(err.i, err.line)
                    logging.warning(help_message)
                    MessageBox.warning(
                        "File error",
                        f"Error reading file at line {err.i} ({help_message})",
                    )
                else:
                    help_message = text.format("NA", "NA")
                    MessageBox.warning(
                        "File error", f"Error reading file ({help_message})"
                    )
                e = err
            if e:
                logging.exception(e, exc_info=True)
                return False

            if self.obsTreeModel.rowCount() > 0:
                if open_type != "CG5":
                    self.populate_station_coords()
                self.workspace_loaded = True
                self.set_window_title_asterisk()
                self.update_menus()
                QtWidgets.QApplication.processEvents()
                QtWidgets.QApplication.restoreOverrideCursor()
            else:
                QtWidgets.QApplication.restoreOverrideCursor()
                MessageBox.warning("Unknown import error", "File error")
        else:
            QtWidgets.QApplication.restoreOverrideCursor()
            return False

        return True

    @staticmethod
    def read_raw_data_file(fname, meter_type):
        """
        Read raw relative-gravity text file in the format exported from meter
        (Scintrex or Burris). Data are returned to the calling function.

        Parameters
        ----------
        fname : str
            Full path to file
        meter_type : {'Burris', 'CG5', 'CG6', 'CG6TSoft, 'csv'}

        Returns
        -------
        all_survey_data : ChannelList
            Object with all survey data

        """
        try:
            with open(fname, "r") as fh:
                logging.info(
                    "number of lines: {:d}".format(
                        len([1 for line in open(fname, "r")])
                    )
                )
                return file_reader(meter_type, fh)

        # Returning e like this allows exceptions to be tested in pytest
        except (InvalidMeterException, IOError, ValueError, IndexError) as e:
            raise e

    def workspace_append(self):
        """
        Append previously-saved workspace to current workspace.
        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            None,
            "Open File",
            self.settings.value("current_dir"),
            filter="Workspace files (*.gsa)",
        )
        if fname:
            self.settings.setValue("current_dir", os.path.dirname(fname))
            if fname[-1] == "p":
                MessageBox.warning(
                    "Import error",
                    "If trying to append a .p file, please save it as a .gsa file"
                    " first.",
                )
                return
            obstreesurveys, coords = self.obsTreeModel.load_workspace(fname)
            # TODO: Do something with coords (append to table if not already there)
            for survey in obstreesurveys:
                self.obsTreeModel.appendRow(
                    [survey, QtGui.QStandardItem("0"), QtGui.QStandardItem("0")]
                )
            self.update_all_drift_plots()
            self.populate_station_coords()
            self.workspace_loaded = True
            QtWidgets.QApplication.restoreOverrideCursor()
            self.set_window_title_asterisk()
            self.update_menus()

    def workspace_clear(self, confirm=True):
        """
        Clears all models and refreshes view.
        """
        if confirm:
            if self.windowTitle()[-1] == "*":
                quit_msg = (
                    "The workspace isn't saved. Are you sure you want to clear all"
                    " data?"
                )
                reply = MessageBox.question(
                    "Message",
                    quit_msg,
                )
                if reply == QtWidgets.QMessageBox.No:
                    return

        logging.info("Workspace cleared")
        self.obsTreeModel = ObsTreeModel()
        self.gui_data_treeview.setModel(None)
        self.gui_data_treeview.update()
        self.tab_data.clear_axes()
        self.tab_data.data_view.setModel(None)
        self.tab_data.data_view.update()
        self.tab_drift.delta_view.model().init_data([])
        self.tab_drift.delta_view.update()
        self.tab_drift.dg_samples_view.model().sourceModel().init_data([])
        self.tab_drift.dg_samples_view.update()
        self.tab_drift.tare_view.model().init_data([])
        self.tab_drift.tare_view.update()
        self.tab_drift.clear_axes()
        self.tab_drift.drift_plot_weighted.setCheckState(0)
        self.tab_drift.drift_cont_startendcombobox.setCurrentIndex(0)
        self.tab_drift.drift_polydegree_combobox.setCurrentIndex(0)
        self.tab_drift.driftmethod_combobox.setCurrentIndex(0)
        self.tab_drift.tension_slider.setValue(1250)
        self.tab_adjust.delta_view.model().sourceModel().init_data([])
        self.tab_adjust.delta_view.update()
        self.tab_adjust.datum_view.model().sourceModel().init_data([])
        self.tab_adjust.datum_view.update()
        self.tab_adjust.results_view.model().sourceModel().init_data([])
        self.clear_adjustment_text()
        self.setWindowTitle("GSadjust")
        self.update_menus()

    def workspace_save(self):
        """
        Saves data if a workspace has already been saved
        """
        fname = self.obsTreeModel.save_workspace(self.workspace_savename)
        if not fname:
            logging.info("Workspace save error.")
            MessageBox.warning("Error", "Workspace save error")
            return

        logging.info(f"Workspace saved: {fname}.")
        MessageBox.information(
            "GSadjust",
            "Workspace saved",
        )
        self.set_window_title(fname)
        return True

    def workspace_save_as(self):
        """
        Saves data object using json.dump()
        """
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(
            None,
            "Save workspace as",
            self.settings.value("current_dir"),
            filter="Workspace files (*.gsa)",
        )
        if fname:
            save_name = None
            self.settings.setValue("current_dir", os.path.dirname(fname))
            try:
                save_name = self.obsTreeModel.save_workspace(fname)
            except Exception as e:
                logging.exception(e, exc_info=True)
            else:
                self.menus.mnFileSaveWorkspace.setEnabled(False)

            if not save_name:
                MessageBox.warning("Error", "Workspace save error")
                return

            self.set_window_title(fname)
            MessageBox.information(
                "GSadjust",
                "Workspace saved",
            )
            self.workspace_savename = fname
            self.update_menus()

    def workspace_open_getjson(self):
        """
        Gets filename to open and asks whether to append or overwrite, if applicable.
        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            None,
            "Open File",
            self.settings.value("current_dir"),
            filter="Workspace files (*.gsa)",
        )
        if not fname:
            return

        self.settings.setValue("current_dir", os.path.dirname(fname))
        if self.obsTreeModel.invisibleRootItem().rowCount() > 0:
            overwrite_tree_dialog = DialogOverwrite()
            if overwrite_tree_dialog.exec_():
                self.workspace_clear()
                self.workspace_open_json(fname)
            else:
                return
        else:
            self.workspace_open_json(fname)

    def workspace_open_json(self, fname):
        """
        Loads data from JSON file. Restores PyQt tables to Survey object
        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        QtWidgets.QApplication.processEvents()
        start = time.time()
        obstreesurveys, coords = self.obsTreeModel.load_workspace(fname)
        if obstreesurveys:
            self.workspace_savename = fname
            self.populate_obstreemodel(obstreesurveys)
            self.adjust_update_required()
            self.set_window_title(fname)
        end = time.time()
        logging.info("Load duration %s: %.2f secs", fname, end - start)
        if coords:
            self.obsTreeModel.station_coords = coords
        self.update_menus()

    def populate_obstreemodel(self, obstreesurveys):
        """
        Only called for loading workspace, not when loading raw data.

        Parameters
        ----------
        obstreesurveys : list
            List of ObsTreeSurvey objects

        """
        pbar = QtWidgets.QProgressDialog(
            labelText="Loading workspace", minimum=0, maximum=4
        )
        pbar.setWindowTitle("GSadjust")
        pbar.setLabelText("Building Observation Tree")
        pbar.show()
        saved_deltas = []
        for survey in obstreesurveys:
            self.obsTreeModel.appendRow(
                [survey, QtGui.QStandardItem("0"), QtGui.QStandardItem("0")]
            )
            saved_deltas += survey.deltas
        # For opening legacy workspaces. Exceptions can occur in various places:
        try:
            saved_deltas_dict = {
                d["key"]: (d["checked"], d["adj_sd"], d["assigned_dg"])
                for d in saved_deltas
            }
            saved_assigned_deltas_dict = {
                d["key"]: (d["checked"], d["adj_sd"], d["assigned_dg"])
                for d in saved_deltas
                if d["type"] == "assigned"
            }
            pbar.setLabelText("Building Observation Tree")
            QtWidgets.QApplication.processEvents()

            i = 0
            firststation = None
            # This avoids an error when the first loop (or subsequent loops) are empty
            firstsurvey = self.obsTreeModel.itemFromIndex(self.obsTreeModel.index(0, 0))
            while firststation is None:
                firstloop = firstsurvey.child(i)
                firststation = firstloop.child(0)
                i += 1
            pbar.setValue(1)
            pbar.setLabelText("Updating drift plots")
            QtWidgets.QApplication.processEvents()
            self.select_first_treeview_item()
            try:
                self.populate_station_coords()
            except Exception as e:
                logging.exception(str(e))
                # sometimes coordinates aren't valid
                pass

            # This is going to create new deltas on the drift and adjust tabs:
            self.update_all_drift_plots()

            # After creating the new deltas, we want to apply these attributes
            # that might have been user-specified:
            # - a delta was unchecked
            # - the delta was type='assigned' and assigned_dg = float
            # - the adj_sd was modified
            #
            # We do this by matching up the new delta with the
            # corresponding old delta and copying over those items.
            new_deltas = self.obsTreeModel.deltas()
            new_deltas_dict = {d.key: d for d in new_deltas}

            # Update adj_sd and checked
            for key, delta in saved_deltas_dict.items():
                new_deltas_dict[key].adj_sd = delta[1]
                new_deltas_dict[key].checked = delta[0]

            # Create type = 'assigned' deltas
            for key, delta in saved_assigned_deltas_dict.items():
                new_deltas_dict[key].type = "assigned"
                new_deltas_dict[key].assigned_dg = delta[2]
        except (KeyError, AttributeError) as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            MessageBox.warning(
                "Import warning",
                "This workspaced was saved by an "
                "earlier version of GSadjust. The station data and "
                "drift options will be loaded, but not the deltas "
                "on the adjust tab, which will need to be re-created"
                " with one of the 'Populate delta table' options. "
                "Please verify workspace for accuracy.",
            )
            # The time format for some old workspaces is different:
            self.temp_shifttime(obstreesurveys)
            QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)

        pbar.setValue(2)
        pbar.setLabelText("Populating delta tables")
        QtWidgets.QApplication.processEvents()

        pbar.setValue(3)
        pbar.setLabelText("Initializing GUI")
        QtWidgets.QApplication.processEvents()
        self.init_gui()
        self.update_adjust_tables()
        QtWidgets.QApplication.restoreOverrideCursor()
        pbar.setValue(4)
        self.workspace_loaded = True
        pbar.close()

    def temp_shifttime(self, obstreesurveys):
        """
        For some old workspaces we need to shift time to matplotlib format.
        """
        for survey in obstreesurveys:
            for loop_idx in range(len(survey.loops())):
                loop = survey.child(loop_idx)
                for i in range(loop.rowCount()):
                    obstreestation = loop.child(i)
                    if obstreestation.t[0] > 50000:
                        obstreestation.t = [s - 719163 for s in obstreestation.t]
                for tare in loop.tares:
                    if tare.t[0] > 50000:
                        tare.t = [s - 719163 for s in tare.t]
                self.process_tares(loop)
            survey.deltas = []
        return

    def populate_station_coords(self):
        """
        Stores a single set of coordinates for each station with the
        obsTreeModel object. The coordinates of the last Station in the
        Survey > Loop > Station hierarchy are used.

        Used as a slot for the self.station_model.signal_update_coordinates
        signal (thus init_station_coords can't be called directly)
        """
        self.obsTreeModel.station_coords = init_station_coords_dict(self.obsTreeModel)

    ###########################################################################
    # General routines
    ###########################################################################
    def populate_deltamodel(self, populate_type):
        """
        Called from menu item

        Parameters
        ----------
        populate_type : {'all', 'selectedLoop', 'selectedSurvey'}

        """
        table_updated = False

        if populate_type == "all":
            self.update_all_drift_plots()
            for survey in self.obsTreeModel.checked_surveys():
                table_updated = survey.populate_delta_model(clear=True)
                self.set_adj_sd(survey, survey.adjustment.adjustmentoptions)
        elif populate_type == "selectedsurvey":
            survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
            self.update_survey_drift_plots(survey)
            table_updated = survey.populate_delta_model(clear=True)
            self.set_adj_sd(survey, survey.adjustment.adjustmentoptions)
        elif populate_type == "selectedloop":
            selected_idx = self.gui_data_treeview.selectedIndexes()
            # There may be one, or multiple loops selected. If only one is
            # selected, we'll populate the delta table based on the
            # currentLoopIndex (which will be in bold but not necessarily highlighted).
            if len(selected_idx) >= 4:
                selected_items = []
                for i in selected_idx:
                    selected_items.append(self.obsTreeModel.itemFromIndex(i))
                # selected_items will contain 3 entries for every tree view
                # item (one for the name, plus 2 for
                # g and std. dev. First, decimate to just the name entries
                selected_items = selected_items[::3]
                loops = [item for item in selected_items if type(item) == ObsTreeLoop]
                first = True
                for loop in loops:
                    self.update_loop_drift_plots(loop)
                    survey = loop.parent()
                    if first:
                        table_updated = survey.populate_delta_model(loop, clear=True)
                        first = False
                    else:
                        table_updated = survey.populate_delta_model(loop, clear=False)
            else:
                loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
                survey = loop.parent()
                table_updated = survey.populate_delta_model(loop, clear=True)
            self.set_adj_sd(survey, survey.adjustment.adjustmentoptions)

        if table_updated:
            self.adjust_update_required()
            self.update_adjust_tables()
            self.delta_model.layoutChanged.emit()

        self.menus.set_state(MENU_STATE.DELTA_MODEL)

    def activate_survey_or_loop(self, index):
        """
        Highlights active survey or loop in tree view.

        Parameters
        ----------
        index : PyQt index
            Newly-highlighted tree item, sent by doubleClicked event
        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        item = self.obsTreeModel.itemFromIndex(index)
        try:
            # If a loop:
            if type(item) == ObsTreeLoop:
                if self.previous_loop is not None:
                    self.previous_loop.fontweight = QtGui.QFont.Normal
                    self.previous_loop.cellcolor = Qt.white
                self.previous_loop = item
                item.cellcolor = Qt.lightGray
                item.fontweight = QtGui.QFont.Bold
                self.index_current_loop = index
                self.update_drift_tables_and_plots(update_adjust_tables=False)

            # If a survey
            elif type(item) == ObsTreeSurvey:
                if self.previous_survey is not None:
                    self.previous_survey.fontweight = QtGui.QFont.Normal
                self.previous_survey = item
                self.index_current_survey = index
                item.fontweight = QtGui.QFont.Bold
                self.update_adjust_tables()
        except Exception as e:
            logging.exception(e, exc_info=True)

        self.obsTreeModel.layoutChanged.emit()
        self.adjust_update_not_required()
        QtWidgets.QApplication.restoreOverrideCursor()

    def set_window_title(self, fname):
        self.setWindowTitle("GSadjust - " + fname)

    def set_window_title_asterisk(self):
        """
        Indicates a workspace has not been saved.
        """
        title = self.windowTitle()
        last_char = title[-1]
        if last_char != "*":
            title += "*"
        self.setWindowTitle(title)

    def update_menus(self):
        """
        Set enabled/disabled states for menus.
        """
        if self.workspace_savename and self.windowTitle()[-1] == "*":
            self.menus.set_state(MENU_STATE.ACTIVE_WORKSPACE)
        else:
            self.menus.set_state(MENU_STATE.NO_ACTIVE_WORKSPACE)
        if self.obsTreeModel.rowCount() >= 1:
            self.menus.set_state(MENU_STATE.AT_LEAST_ONE_SURVEY)
            if self.obsTreeModel.invisibleRootItem().rowCount() > 1:
                self.menus.set_state(MENU_STATE.MORE_THAN_ONE_SURVEY)
                if not self.label_adjust_update_required_set:
                    self.menus.set_state(MENU_STATE.CALCULATE_CHANGE)
            try:
                current_survey = self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
                )
                if len(current_survey.deltas) > 0:
                    self.menus.set_state(MENU_STATE.SURVEY_HAS_DELTAS)
                else:
                    self.menus.set_state(MENU_STATE.SURVEY_HAS_NO_DELTAS)
                if (
                    len(current_survey.results) > 0
                    and not self.label_adjust_update_required_set
                ):
                    self.menus.set_state(MENU_STATE.SURVEY_HAS_RESULTS)
                else:
                    self.menus.set_state(MENU_STATE.SURVEY_HAS_NO_RESULTS)
            except TypeError:
                # catches during PyTest
                return
            except AttributeError:
                # catches if no delta_model
                self.menus.set_state(MENU_STATE.UNINIT)
        else:
            self.menus.set_state(MENU_STATE.UNINIT)

        QtWidgets.QApplication.restoreOverrideCursor()

    def update_all_drift_plots(self):
        """
        Updates drift_tab plots and delta_models, even if not in view.
        """
        orig_loop_index = self.index_current_loop
        for survey in self.obsTreeModel.surveys():
            self.update_survey_drift_plots(survey)
        self.index_current_loop = orig_loop_index

    def update_survey_drift_plots(self, survey):
        for loop in survey.loops():
            self.update_loop_drift_plots(loop)

    def update_loop_drift_plots(self, loop):
        self.index_current_loop = loop.index()
        self.update_drift_tables_and_plots(update=False)

    def show_delta_update_message(self):
        MessageBox.warning(
            "Delta error",
            "The Delta-g value can only be edited if the drift-correction method for"
            " the respective loop is a method other than the Roman method.",
        )

    def update_adjust_tables(self):
        """
        Update delta-g and datum tables after selecting a new survey in the tree
        view, after a network adjustment, or after changing the drift-correction method
        """
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)

        if survey:
            self.delta_model.init_data(survey.deltas)
            self.datum_model.init_data(survey.datums)
            self.results_model.init_data(survey.results)
            stats_model = QtGui.QStandardItemModel()
            if survey.adjustment.adjustmentresults.n_unknowns > 0:  # Numpy adjustment
                stats_model.setColumnCount(2)
                stats_model.setHorizontalHeaderLabels(["", ""])
                for line in survey.adjustment.results_string():
                    try:
                        line_elems = line.split(":")
                        if len(line_elems) == 2:  # normal line, stat: value
                            stats_model.appendRow(
                                [
                                    QtGui.QStandardItem(line_elems[0]),
                                    QtGui.QStandardItem(line_elems[1].strip()),
                                ]
                            )
                        else:  # Chi test accepted or rejected. No ":"
                            text = line_elems[0].strip()
                            qt_item = QtGui.QStandardItem(text)
                            if "rejected" in text:
                                qt_item.setForeground(Qt.red)
                            stats_model.appendRow([qt_item, QtGui.QStandardItem("")])
                    except:
                        pass
                self.tab_adjust.stats_view.setModel(stats_model)
                self.tab_adjust.stats_view.setColumnWidth(0, 350)
                self.tab_adjust.stats_view.setColumnWidth(1, 150)
            elif survey.adjustment.adjustmentresults.text:  # Gravnet adjustment
                for line in survey.adjustment.adjustmentresults.text:
                    stats_model.appendRow([QtGui.QStandardItem(line)])
                self.tab_adjust.stats_view.setModel(stats_model)
                self.tab_adjust.stats_view.setColumnWidth(0, 600)
                stats_model.setColumnCount(1)
                stats_model.setHorizontalHeaderLabels([""])
            self.tab_adjust.update_col_widths()

    def update_drift_tables_and_plots(self, update=True, update_adjust_tables=True):
        """
        First updates the drift_method combobox, then calls set_drift_method
        to update plots.

        Parameters
        ----------
        update: bool
            Plots are only updated if True. Saves time when loading a workspace.
        update_adjust_tables : bool
            When loading a workspace, we don't want to update adjust tables until
            after the deltas have been created on the Drift tab.

        """
        current_loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        current_survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)

        drift_method = current_loop.drift_method
        self.tab_drift.driftmethod_combobox.setCurrentIndex(
            self.drift_lookup[drift_method]
        )
        self.tab_drift.tare_view.model().init_data(current_loop.tares)
        self.tab_drift.set_drift_method(update, update_adjust_tables)
        # self.set_adj_sd(current_survey, current_survey.adjustment.ao, loop=current_loop)
        self.adjust_update_required()

    def on_obs_checked_change(self, selected):
        """
        Called when a checkbox state is changed, but not when a new item selected.
        Should update drift plots if on drift tab, but otherwise do nothing.

        Parameters
        ----------
        selected : list
            List of selected indexes
        """
        if selected.model() is not None:
            self.adjust_update_required()
            # if self.tab_widget.currentIndex() == 0:
            self.update_drift_tables_and_plots()
            # elif self.tab_widget.currentIndex() == 1:
            self.update_drift_tables_and_plots()
            # elif self.tab_widget.currentIndex() == 2:
            self.update_drift_tables_and_plots()
            # self.tab_adjust.delta_view.update()
        self.set_window_title_asterisk()
        self.update_menus()

    def on_obs_tree_change(self, selected):
        """
        Called when the selection model changes.

        Parameters
        ----------
        selected : list
            List of selected indexes
        """
        indexes = selected.indexes()
        if indexes:
            item = self.obsTreeModel.itemFromIndex(indexes[0])
            if item:
                if isinstance(item, ObsTreeStation):
                    self.index_current_station = indexes[0]
                    if self.tab_widget.currentIndex() == 0:
                        self.update_data_tab()

    def correct_recorded_time(self):
        """
        Correct all times from an offset, e.g. when local time is used when UTC time
        should have. Can be applied to a station, loop, etc.
        """
        text, ok = QtWidgets.QInputDialog.getText(
            self, "Input parameters", "time offset to apply (min)?"
        )
        if ok:
            try:
                text = float(text)
            except ValueError:
                MessageBox.warning(
                    "Time correction error", "Invalid value entered for time correction"
                )
                return
            time_correction_dialog = DialogApplyTimeCorrection()
            time_correction_dialog.msg.exec_()
            correction_type = time_correction_dialog.time_correction_type
            if correction_type:
                time_correction(
                    self.obsTreeModel,
                    correction_type,
                    text,
                    self.index_current_survey,
                    self.index_current_loop,
                    self.gui_data_treeview.selectedIndexes(),
                )
                self.update_data_tab()
                self.set_window_title_asterisk()

    def set_vertical_gradient_interval(self):
        """
        Dialog that queries user for the distance over which vertical gradient
        is measured.
        """
        interval = VerticalGradientDialog(self.vertical_gradient_interval)
        if interval.ok:
            self.vertical_gradient_interval = interval.text

    def vertical_gradient_write(self):
        """
        Writes a .grd file with two values: gradient and standard deviation.
        Only works when Roman drift method is used.
        """
        # current_loop_index = self.obsTreeModel.selectedIndexes()[0]
        current_loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        deltas = current_loop.deltas
        n_stations = current_loop.n_unique_stations()

        dg, sd = None, None
        if n_stations == 2:
            stationname = (
                self.obsTreeModel.itemFromIndex(self.index_current_loop)
                .child(0)
                .station_name
            )
            defaultfile = os.path.join(
                self.settings.value("current_dir"), stationname + ".grd"
            )
            try:
                if current_loop.drift_method == "roman":
                    filename, _ = QtWidgets.QFileDialog.getSaveFileName(
                        None, "Vertical gradient file to write", defaultfile
                    )
                    if filename:
                        delta = deltas[0]
                        dg = delta.dg
                        sd = delta.sd
                elif (
                    current_loop.drift_method == "none"
                    or current_loop.drift_method == "continuous"
                ):
                    filename, _ = QtWidgets.QFileDialog.getSaveFileName(
                        None, "Vertical gradient file to write", defaultfile
                    )
                    if filename:
                        dg_list, sd_list = [], []
                        for delta in current_loop.deltas:
                            dg_list.append(np.abs(delta.dg))
                            sd_list.append(delta.sd)
                        dg = np.mean(dg_list)
                        sd = np.mean(sd_list)
                else:
                    MessageBox.warning(
                        "Vertical gradient error",
                        "When writing a vertical gradient the drift correction method"
                        ' must be "None", "Roman", or "Continuous"',
                    )
                if dg:
                    with open(filename, "w") as fid:
                        fid.write(
                            "{:0.2f}".format(-1 * dg / self.vertical_gradient_interval)
                        )
                        fid.write(
                            " +/- {:0.2f}".format(sd / self.vertical_gradient_interval)
                        )
            except:
                MessageBox.warning(
                    "Vertical gradient error",
                    "Error writing gradient file (try visiting the Drift tab then"
                    " re-running this command).",
                )
        else:
            MessageBox.warning(
                "Vertical gradient error",
                "Writing a vertical gradient file requires that the respective loop has"
                ' only two stations. The gradient interval is set by the "Vertical'
                ' gradient interval..." menu command (instrument height in the meter'
                " file is ignored).",
            )

    def tares_updated(self):
        """
        Called when tab_drift.tare_view.model() is updated.
        """
        current_loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        self.process_tares(current_loop)
        self.update_drift_tables_and_plots()
        self.set_window_title_asterisk()

    def add_tare(self):
        """
        Opens a dialog to add tare to loop tare_model.
        """
        new_tare_date, new_tare_value = 0, 0
        current_loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        obstreestation = current_loop.child(0)
        default_time = num2date(obstreestation.t[0])
        taredialog = AddTareDialog(default_time)
        if taredialog.exec_():
            new_tare_date = taredialog.dt_edit.dateTime()
            new_tare_value = taredialog.edit_box.text()
        try:
            tare = Tare(date2num(new_tare_date.toPyDateTime()), new_tare_value)
        except:
            return

        current_loop.tares.append(tare)
        self.tab_drift.tare_view.model().init_data(current_loop.tares)
        self.process_tares(current_loop)
        self.update_drift_tables_and_plots()
        self.set_window_title_asterisk()

    @staticmethod
    def process_tares(obstreeloop):
        """
        Apply tares in the tare table.

        Parameters
        ----------
        obstreeloop : ObsTreeLoop
            Loop shown in drift tab

        """
        for obstreestation in obstreeloop.stations():
            for idx, t in enumerate(obstreestation.tare):
                obstreestation.tare[idx] = 0
            for tare in obstreeloop.tares:
                if tare.checked == 2:
                    for idx, t in enumerate(obstreestation.t):
                        if t > tare.datetime:
                            obstreestation.tare[idx] += tare.tare

    def clear_delta_model(self):
        """
        Remove all deltas from survey delta model shown on network adjustment tab.
        """
        self.tab_adjust.delta_view.model().sourceModel().init_data([])
        self.tab_adjust.delta_view.update()
        self.tab_adjust.results_view.model().sourceModel().init_data([])
        self.tab_adjust.results_view.update()
        obstreesurvey = self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
            )
        obstreesurvey.deltas = []
        self.clear_adjustment_text()
        self.set_window_title_asterisk()

    def clear_adjustment_text(self):
        """
        Called when delta or datum tables are cleared.
        """
        self.tab_adjust.stats_view.setModel(None)
        self.tab_adjust.stats_view.update()
        try:
            survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
            survey.adjustment.adjustmentresults.text = []
            survey.adjustment.adjustmentresults.n_unknowns = 0
        except (AttributeError, TypeError):
            return

    def clear_datum_model(self):
        """
        Remove all datums from datum model shown on network adjustment tab.
        """
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        survey.datums = []
        self.tab_adjust.datum_view.model().sourceModel().init_data([])
        self.tab_adjust.datum_view.update()
        self.tab_adjust.results_view.model().sourceModel().init_data([])
        obstreesurvey = self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
            )
        obstreesurvey.datums = []
        self.clear_adjustment_text()
        self.set_window_title_asterisk()

    def clear_results_model(self):
        """
        Remove all results from results model shown on network adjustment tab.
        """
        self.tab_adjust.results_view.model().sourceModel().init_data([])
        obstreesurvey = self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
            )
        obstreesurvey.results = []
        self.update_adjust_tables()

    def animate_loop(self):
        """
        Animated map of loop observations, called from context menu
        """
        loop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        coords = self.obsTreeModel.station_coords
        if not coords:
            MessageBox.warning("GSadjust error", "No station coordinates")
        else:
            lat, lon, dates = [], [], []
            try:
                for station in loop.checked_stations():
                    lat.append(coords[station.station_name][1])
                    lon.append(coords[station.station_name][0])
                    dates.append(station.tmean())
                self._plot_loop_animation = PlotLoopAnimation([lat, lon, dates])
                self._plot_loop_animation.show()
            except Exception:
                MessageBox.warning(
                    "GSadjust error", "Unknown error, station {}".format(station)
                )

    def properties_loop(self):
        """
        Show popup dialog to specify loop properties.
        """
        indexes = self.gui_data_treeview.selectedIndexes()
        loops = []
        for idx in indexes:
            if idx.column() == 0:
                loops.append(self.obsTreeModel.itemFromIndex(idx))
        loop_options = DialogLoopProperties(loops, parent=self)
        if loop_options.exec_():
            # Sync new meter numbers with station objects
            for loop in loops:
                loop.comment = loop_options.comment_edit.toPlainText()
                for i in range(loop.rowCount()):
                    obstreestation = loop.child(i)
                    obstreestation.meter = [loop_options.meter_edit.text()] * len(
                        obstreestation.meter
                    )
                    obstreestation.oper = [loop_options.operator_edit.text()] * len(
                        obstreestation.oper
                    )
            self.set_window_title_asterisk()
        self.update_data_tab()

    def delete_survey(self):
        """
        Remove loop or survey from tree view
        """
        index = self.gui_data_treeview.selectedIndexes()

        if index[0] == self.index_current_station_survey:
            update_selected_station = True
        else:
            update_selected_station = False

        if index[0] == self.index_current_survey:
            update_selected_survey = True
            old_row = index[0].row()
        else:
            update_selected_survey = False

        if len(self.obsTreeModel.surveys()) == 1:
            self.workspace_clear(confirm=False)
            self.update_menus()
            return

        self.obsTreeModel.beginRemoveRows(
            index[0].parent(), index[0].row(), index[0].row() + 1
        )
        self.obsTreeModel.removeRows(index[0].row(), 1, index[0].parent())
        self.obsTreeModel.endRemoveRows()

        if update_selected_survey:
            if old_row == 0:
                self.index_current_survey = self.obsTreeModel.index(0, 0)
            else:
                self.index_current_survey = self.obsTreeModel.index(old_row - 1, 0)
            if update_selected_station:
                obstreesurvey = self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
                )
                first_loop = obstreesurvey.child(0)
                first_station = first_loop.child(0)
                self.index_current_station = first_station.index()
                self.update_data_tab()
        self.activate_survey_or_loop(self.index_current_survey)
        self.set_window_title_asterisk()
        self.update_menus()

    def delete_loop(self):
        """
        Remove loop or survey from tree view
        """
        index = self.gui_data_treeview.selectedIndexes()
        obstreesurvey = self.obsTreeModel.itemFromIndex(index[0]).parent()
        loop_to_delete = self.obsTreeModel.itemFromIndex(index[0]).name
        # store bool True / False.
        update_selected_station = index[0] == self.index_current_station_loop

        if index[0] == self.index_current_loop:
            update_selected_loop = True
            old_row = index[0].row()
        else:
            update_selected_loop = False

        self.obsTreeModel.beginRemoveRows(
            index[0].parent(), index[0].row(), index[0].row() + 1
        )
        self.obsTreeModel.removeRows(index[0].row(), 1, index[0].parent())
        self.obsTreeModel.endRemoveRows()

        if obstreesurvey.rowCount() == 0:
            self.index_current_loop = None
            self.tab_drift.reset()
        else:
            if update_selected_loop:
                if old_row == 0:
                    loop = obstreesurvey.child(0)
                    self.index_current_loop = loop.index()
                else:
                    loop = obstreesurvey.child(old_row - 1)
                    self.index_current_loop = loop.index()
                if update_selected_station:
                    first_station = loop.child(0)
                    self.index_current_station = first_station.index()
                    self.update_data_tab()
                self.activate_survey_or_loop(self.index_current_loop)
        obstreesurvey.deltas = [
            d for d in obstreesurvey.deltas if d.loop != loop_to_delete
        ]
        self.update_adjust_tables()
        self.set_window_title_asterisk()
        self.update_menus()

    def delete_station(self):
        """
        Remove station from tree view
        """
        indexes = self.gui_data_treeview.selectedIndexes()

        # Because each tree item has three columns, len(indexes) equals the
        # number of items selected * 3. The next line takes every 3rd index.
        indexes = indexes[0::3]
        for index in reversed(indexes):
            self.obsTreeModel.removeRow(index.row(), index.parent())
            self.set_window_title_asterisk()
        if index.row() > 0:
            self.index_current_station = index.sibling(index.row() - 1, 0)
        else:
            self.index_current_station = index.sibling(0, 0)

        first_index = indexes[0]
        if first_index.row() > 0:
            row = first_index.row() - 1
        else:
            row = 0
        new_selection_index = self.obsTreeModel.index(row, 0, first_index.parent())
        self.selection_model.select(
            new_selection_index, QtCore.QItemSelectionModel.SelectCurrent
        )
        self.update_data_tab()
        self.update_drift_tables_and_plots()
        self.update_adjust_tables()
        self.update_menus()

    def rename(self):
        """
        Rename station; same as F2.
        """
        indexes = self.gui_data_treeview.selectedIndexes()
        index = indexes[0]
        trigger = self.gui_data_treeview.EditKeyPressed
        event = None
        self.gui_data_treeview.edit(index, trigger, event)
        self.set_window_title_asterisk()

    def delete_tare(self):
        """
        Called when user right-clicks a tare and selects delete from the context menu
        """
        index = self.tab_drift.tare_view.selectedIndexes()
        for idx in reversed(index):
            self.tab_drift.tare_view.model().removeRow(idx)
        self.tab_drift.tare_view.update()
        self.process_tares(
            self.obsTreeModel.itemFromIndex(self.index_current_loop)
        )
        self.update_drift_tables_and_plots()
        self.set_window_title_asterisk()

    def delete_datum(self):
        """
        Called when user right-clicks a datum and selects delete from the context menu
        """
        index = self.tab_adjust.datum_view.selectedIndexes()
        i = [self.tab_adjust.datum_proxy_model.mapToSource(idx) for idx in index]
        i.sort(key=lambda x: x.row(), reverse=True)
        for idx in i:
            self.tab_adjust.datum_proxy_model.sourceModel().removeRow(idx)
        self.tab_adjust.datum_view.update()
        self.set_window_title_asterisk()

    def divide_by_time(self):
        """
        Divide selected loop into shorter loops; a period less than the user-selected
        threshold indicates a break between loops.
        """
        # Prompt user to select time threshold
        loopdialog = LoopTimeThresholdDialog()
        if loopdialog.exec_():
            loop_thresh = loopdialog.dt_edit.dateTime()
            # Convert to days. Subtract one from the date because the default is 1 (
            # i.e., if the time set in the loop dialog is 8:00, loop thresh is a
            # Qdatetime equal to (2000,1,1,8,0).
            loop_thresh = (
                int(loop_thresh.toString("H")) / 24 + int(loop_thresh.toString("d")) - 1
            )
        else:
            return
        self.divide_survey(loop_thresh)
        self.update_survey_drift_plots(
            self.obsTreeModel.itemFromIndex(self.index_current_survey)
        )
        self.activate_survey_or_loop(self.index_current_loop)

    def divide_by_height(self):
        """
        Called from "Divide loop..." menu command. Shows a dialog to specify a
        time interval, then scans the current loop and divides station
        occupations separated by the time interval (or greater) into a loop.
        Useful primarily when several day's data is in a single file.
        """
        # Clear survey delta table, it causes problems otherwise
        self.clear_delta_model()

        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        obstreeloop = self.obsTreeModel.itemFromIndex(self.index_current_loop)

        if obstreeloop.rowCount() > 1:
            MessageBox.warning(
                "GSadjust error",
                "Loop must have a single station to divide by height.",
            )
            return

        station = obstreeloop.child(0)
        height = station.height[0]
        height_idx = 0
        n_samples = len(station.height)
        count_dict = dict()
        heights = set(station.height)
        for h in heights:
            count_dict[h] = 0
        for idx, h in enumerate(station.height):
            if h != height or idx == len(station.height) - 1:
                new_station = ChannelList()
                for k, v in station.__dict__.items():
                    try:  # some fields have no len()
                        if len(v) == n_samples:
                            setattr(new_station, k, v[height_idx:idx])
                        else:
                            setattr(new_station, k, v)
                    except Exception:
                        continue
                if height == 0:
                    name = station.station_name
                else:
                    name = station.station_name + "_{:.3f}".format(height)
                new_obstreestation = ObsTreeStation(
                    new_station, name, "{}".format(int(count_dict[height]))
                )
                count_dict[height] += 1
                obstreeloop.appendRow(
                    [
                        new_obstreestation,
                        QtGui.QStandardItem("a"),
                        QtGui.QStandardItem("a"),
                    ]
                )
                height_idx = idx
                height = h

        obstreeloop.removeRow(0)
        self.index_current_station = obstreeloop.child(0).index()
        self.update_drift_tables_and_plots()
        self.obsTreeModel.layoutChanged.emit()
        QtWidgets.QApplication.restoreOverrideCursor()
        self.set_window_title_asterisk()

    def divide_survey(self, loop_thresh):
        """
        Called from self.divide_by_time. Shows a dialog to specify a
        time interval, then scans the current loop and divides station occupations
        separated by the time interval (or greater) into a loop. Useful primarily
        when several day's data is in a single file.
        """
        # Clear survey delta table, it causes problems otherwise
        self.clear_delta_model()

        # Store the original current loop index so it can be restored.
        original_loop_index = self.index_current_loop

        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        obstreeloop = self.obsTreeModel.itemFromIndex(self.index_current_loop)
        indexes = []
        pbar = ProgressBar(total=obstreeloop.rowCount() - 1, textmess="Divide loop")
        pbar.show()
        pbar.progressbar.setValue(1)
        QtWidgets.QApplication.processEvents()
        # Step through loop backward
        pbar_idx = list(range(obstreeloop.rowCount()))
        pbar_idx.reverse()
        for i in range(obstreeloop.rowCount() - 1, 1, -1):
            station2 = obstreeloop.child(i)
            station1 = obstreeloop.child(i - 1)
            pbar.progressbar.setValue(pbar_idx[i])
            QtWidgets.QApplication.processEvents()
            # Check time difference between successive stations
            tdiff = station2.tmean() - station1.tmean()
            if tdiff > loop_thresh:
                for ii in range(i, obstreeloop.rowCount()):
                    indexes.append(obstreeloop.child(ii).index())
                self.new_loop_from_indexes(indexes)
                indexes = []
        self.index_current_loop = original_loop_index
        self.update_all_drift_plots()
        pbar.close()
        self.obsTreeModel.layoutChanged.emit()
        QtWidgets.QApplication.restoreOverrideCursor()
        self.set_window_title_asterisk()

    def duplicate_station(self):
        """
        Create a duplicate of a station in the tree view. Useful when the same
        station is observed at the end of one day and the start of the next day:
        when imported, it will appear as one station, but it should be two.
        """
        indexes = self.gui_data_treeview.selectedIndexes()
        if len(indexes) > 3:
            MessageBox.warning(
                "GSadjust error", "Please select a single station when duplicating."
            )
            return
        index = indexes[0]
        model = indexes[0].model()
        obstreeloop = model.itemFromIndex(indexes[0].parent())
        obstreesurvey = obstreeloop.parent()
        obstreestation = self.obsTreeModel.itemFromIndex(index)
        new_station_count = float(obstreestation.station_count) + 0.1
        new_obstreestation = ObsTreeStation(
            obstreestation,
            obstreestation.station_name,
            "{:.1f}".format(new_station_count),
        )

        obstreeloop.insertRow(
            index.row() + 1,
            [new_obstreestation, QtGui.QStandardItem("a"), QtGui.QStandardItem("a")],
        )
        logging.info(
            "Station duplicated: {}, Survey: {}, Loop: {} ".format(
                new_obstreestation.station_name, obstreesurvey.name, obstreeloop.name
            )
        )
        self.set_window_title_asterisk()

    def move_survey(self, direction=UP):
        """
        Used to move survey up or down in the tree view

        Parameters
        ----------
            direction: UP or DOWN (macros for 1 and -1)

        """
        if direction not in DIRECTIONS:
            return

        model = self.obsTreeModel
        index = self.index_current_survey
        row_number = index.row()
        new_row = row_number + direction
        if not (0 <= new_row < model.rowCount()):
            return False
        survey = model.takeRow(row_number)
        model.insertRow(new_row, survey)
        self.index_current_survey = model.indexFromItem(survey[0])
        self.set_window_title_asterisk()
        return True

    def new_loop(self):
        """
        Creates a new loop in tree view
        """
        indexes = self.gui_data_treeview.selectedIndexes()
        self.new_loop_from_indexes(indexes)
        self.activate_survey_or_loop(self.index_current_loop)
        self.update_drift_tables_and_plots(update=True)
        self.update_data_tab()
        self.set_window_title_asterisk()

    def new_loop_from_indexes(self, indexes):
        """
        Moves stations at the specified indexes to a new loop.

        Parameters
        ----------
        indexes : list
            List of PyQt indexes of selected stations
        """
        if len(indexes) > 0:
            model = indexes[0].model()
            obstreeloop = model.itemFromIndex(indexes[0].parent())
            obstreesurvey = obstreeloop.parent()
            new_loop_name = str(obstreesurvey.rowCount())
            # new loop, increment from loop parent
            new_obstreeloop = ObsTreeLoop(new_loop_name)

            logging.info("Loop {} added".format(new_loop_name))

            for idx in indexes:
                if idx.column() == 0:
                    obstreestation = model.itemFromIndex(idx)
                    new_obstreestation = ObsTreeStation(
                        obstreestation,
                        obstreestation.station_name,
                        obstreestation.station_count,
                    )
                    logging.info(
                        "Station added to new loop: {}".format(
                            obstreestation.station_name
                        )
                    )
                    new_obstreeloop.appendRow(
                        [
                            new_obstreestation,
                            QtGui.QStandardItem("a"),
                            QtGui.QStandardItem("a"),
                        ]
                    )
            copy_fields = ["source"]
            for field in copy_fields:
                try:
                    setattr(new_obstreeloop, field, getattr(obstreeloop, field))
                except:
                    setattr(new_obstreeloop, field, "")

            for idx in reversed(indexes):
                if idx.column() == 0:
                    self.obsTreeModel.beginRemoveRows(
                        idx.parent(), idx.row(), idx.row() + 1
                    )
                    self.obsTreeModel.removeRow(idx.row(), idx.parent())
                    self.obsTreeModel.endRemoveRows()
            obstreesurvey.appendRow(new_obstreeloop)
            self.gui_data_treeview.expand(new_obstreeloop.index())
            self.obsTreeModel.layoutChanged.emit()
            self.index_current_loop = new_obstreeloop.index()
            self.index_current_station = obstreeloop.child(0, 0).index()
            self.set_window_title_asterisk()

    def treeview_context_menu(self, point):
        """
        Right-click context menu on tree view
        """
        self.gui_data_treeview_popup_menu = QtWidgets.QMenu("Menu", self)
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnRename)
        self.gui_data_treeview_popup_menu.addSeparator()
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnDeleteSurvey)
        self.gui_data_treeview_popup_menu.addSeparator()
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnDeleteLoop)
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnLoopProperties)
        self.gui_data_treeview_popup_menu.addAction(
            self.menus.mnVerticalGradientWriteAction
        )
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnLoopAnimate)
        self.gui_data_treeview_popup_menu.addSeparator()
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnDeleteStation)
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnStationDuplicate)
        self.gui_data_treeview_popup_menu.addAction(self.menus.mnDataNewLoop)
        # enable as appropriate
        indexes = self.gui_data_treeview.selectedIndexes()

        if len(indexes) == 3:  # One item is selected
            index = indexes[0]
            item = index.model().itemFromIndex(index)
            if type(item) is ObsTreeStation:
                self.menus.set_state(MENU_STATE.OBS_TREE_STATION)
            elif type(item) is ObsTreeLoop:
                self.menus.set_state(MENU_STATE.OBS_TREE_LOOP)
            elif type(item) is ObsTreeSurvey:
                self.menus.set_state(MENU_STATE.OBS_TREE_SURVEY)
        elif len(indexes) > 3:
            index_types = [
                type(index.model().itemFromIndex(index)) for index in indexes[0::3]
            ]
            if all([it == ObsTreeStation for it in index_types]):
                self.menus.set_state(MENU_STATE.MULTIPLE_STATION)
            elif all([it == ObsTreeLoop for it in index_types]):
                self.menus.set_state(MENU_STATE.MULTIPLE_LOOP)
            else:
                self.menus.set_state(MENU_STATE.UNENABLE_ALL)
        else:
            self.menus.set_state(MENU_STATE.UNENABLE_ALL)
        self.gui_data_treeview_popup_menu.exec_(
            self.gui_data_treeview.mapToGlobal(point)
        )

    def adjust_network(self, how_many="all"):
        """
        Carries out network adjustment, updates output tables
        """
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self.menus.mnAdjPyLSQ.isChecked():
            adj_type = "PyLSQ"
        else:
            adj_type = "Gravnet"
        # Collect checked items into adjustment object
        if how_many == "all":
            for obstreesurvey in self.obsTreeModel.checked_surveys():
                obstreesurvey.run_inversion(adj_type)
        elif how_many == "current":
            obstreesurvey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
            obstreesurvey.run_inversion(adj_type)
        self.statusBar().showMessage("Network adjustment complete")
        self.update_adjust_tables()
        self.adjust_update_not_required()
        self.update_menus()
        QtWidgets.QApplication.restoreOverrideCursor()

    def update_SD_and_run_adjustment(self):
        """
        Only enabled after an initial adjustment is run.

        Adjust adj_sd so the posterior S.D equals one.
        """
        obstreesurvey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        ao = copy.deepcopy(obstreesurvey.adjustment.adjustmentoptions)
        if ao.use_sigma_min:
            ao.use_sigma_postfactor = True
            ao.sigma_postfactor = ao.sigma_postfactor * float(
                obstreesurvey.adjustment.SDaposteriori[0]
            )
        else:
            ao.use_sigma_prefactor = True
            ao.sigma_prefactor = ao.sigma_prefactor * float(
                obstreesurvey.adjustment.SDaposteriori[0]
            )
        obstreesurvey.adjustment.adjustmentoptions = ao
        self.set_adj_sd(obstreesurvey, ao)
        obstreesurvey.run_inversion()
        self.update_adjust_tables()
        self.adjust_update_not_required()
        self.update_menus()

    def show_gravity_change_table(self):
        """
        Shows table of gravity differences between surveys.
        """
        data = compute_gravity_change(self.obsTreeModel)
        if data:
            win = GravityChangeTable(self, data, table_type="simple")
            win.show()

    def show_cal_coeff(self):
        """
        Shows table of gravimeter calibration coefficients.
        """
        cal_coeffs = self.obsTreeModel.get_cal_coeffs()
        if cal_coeffs:
            win = ShowCalCoeffs(cal_coeffs, parent=self)
            win.show()

    def plot_network_graph_circular(self):
        """
        Show circular network graph
        """
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        coords = self.obsTreeModel.station_coords
        plt = PlotNetworkGraph(survey, coords, shape="circular", parent=self)
        plt.show()

    def plot_network_graph_map(self):
        """
        Show map-view network graph
        """
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        coords = self.obsTreeModel.station_coords
        plt = PlotNetworkGraph(survey, coords, shape="map", parent=self)
        plt.show()

    def plot_datum_vs_adjusted(self):
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        plt = PlotDatumCompare(survey, self)
        plt.show()

    def plot_datum_comparison_timeseries(self):
        plt = PlotDatumComparisonTimeSeries(self.obsTreeModel, self)
        plt.show()

    def plot_histogram(self):
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        plt = PlotDgResidualHistogram(survey, self)
        plt.show()

    def plot_gravity_change(self, dates, table, parent):
        plt = PlotGravityChange(dates, table, parent)
        plt.show()

    def set_adj_sd(self, survey, ao, loop=None):
        """
        Update delta table based on parameters in net. adj. options

        Parameters
        ----------
        survey : ObsTreeSurvey
        ao : AdjustmentOptions
        """
        if loop:
            for delta in loop.deltas:
                try:
                    delta.adj_sd = self.calc_adj_sd(ao, delta.sd)
                except AttributeError:
                    pass
        else:
            for delta in survey.deltas:
                try:
                    delta.adj_sd = self.calc_adj_sd(ao, delta.sd)
                except AttributeError:
                    pass

    def calc_adj_sd(self, ao, sd):
        """
        Calculates S.D. to show in "SD for adj." column, based on AdjustOptions
        values

        Parameters
        ----------
        ao : AdjustOptions
        sd : float
            Original delta s.d.

        Returns
        -------
        float
            Adjusted standard deviation

        """
        additive = 0
        factor = 1
        if ao.use_sigma_add:
            additive = ao.sigma_add
        if ao.use_sigma_prefactor:
            factor = ao.sigma_prefactor
        if ao.use_sigma_min:
            sigma = max(ao.sigma_min, sd * factor + additive)
            if ao.use_sigma_postfactor:
                sigma *= ao.sigma_postfactor
        else:
            sigma = sd * factor + additive
        return sigma

    def menu_import_abs_g_simple(self):
        """
        Slot for 'Import abs. g (simple)' menu command
        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            None,
            "Open file (3 columns, space delimited, station-g-std. dev.)",
            self.settings.value("abs_g_path"),
        )
        if fname:
            self.settings.setValue("abs_g_path", os.path.dirname(fname))
            logging.info("Importing absolute gravity data from {}".format(fname))
            try:
                datums = import_abs_g_simple(fname)
            except IndexError:
                MessageBox.warning(
                    "File read error",
                    "Error reading absolute gravity file. Is it three columns (station,"
                    " g, std. dev.), space delimited",
                )
            for datum in datums:
                self.obsTreeModel.itemFromIndex(
                    self.index_current_survey
                ).datum_model.insertRows(datum, 0)
                logging.info("Datum imported: {}".format(datum.__str__()))
            self.set_window_title_asterisk()

    def menu_import_abs_g_complete(self):
        """
        Slot for 'Import abs. g (complete)' menu command
        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            None, "Open A10_parse.py output file", self.settings.value("abs_g_path")
        )

        if fname:
            self.settings.setValue("abs_g_path", os.path.dirname(fname))
            logging.info("Importing absolute gravity data from {}".format(fname))
            datums = import_abs_g_complete(fname)
            survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
            survey.datums += datums
            for datum in datums:
                logging.info("Datum imported: {}".format(datum.__str__()))
            self.update_adjust_tables()
            self.set_window_title_asterisk()
        else:
            return

    def dialog_import_abs_g_direct(self):
        """
        Opens a PyQt dialog to select a directory with .project.txt files.
        """
        if self.obsTreeModel.rowCount() == 0:
            MessageBox(
                "Import error",
                "Please load a survey before loading absolute-gravity data.",
            )
            return
        if hasattr(self, "abs_import_table_model"):
            selectabsg = SelectAbsg(
                self.settings.value("abs_g_path"),
                datum_table_model=self.abs_import_table_model,
            )
        else:
            selectabsg = SelectAbsg(self.settings.value("abs_g_path"))
        if selectabsg.exec_():
            nds = selectabsg.new_datums
            self.abs_import_table_model = selectabsg.table_model
            self.settings.setValue("abs_g_path", selectabsg.path)
            if nds:
                survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
                survey.datums += nds
                self.set_window_title_asterisk()
        else:
            self.abs_import_table_model = selectabsg.table_model
        self.update_adjust_tables()
        self.update_menus()

    def write_metadata_text(self):
        """
        Exports processing summary for metadata file.
        """
        fn = export_metadata(self.obsTreeModel, self.settings.value("current_dir"))
        if fn:
            MessageBox.information("Data export", "Metadata written to {}".format(fn))
        else:
            MessageBox.warning("Write error", "No network adjustment results")

    def write_summary_text(self):
        """
        Write complete summary of data and adjustment, with the intent that the
        processing can be re-created later
        """
        fn = export_summary(self.obsTreeModel, self.settings.value("current_dir"))
        if fn:
            MessageBox.information(
                "Data export", "Text summary written to {}".format(fn)
            )

    def write_tabular_data(self):
        """
        Write data to file
        """
        fn = export_data(self.obsTreeModel, self.settings.value("current_dir"))
        if fn:
            MessageBox.information(
                "Data export", "Tabular data written to {}".format(fn)
            )

    def dialog_add_datum(self):
        """
        Opens PyQt dialog to select an existing station to assign a datum value
        """
        stations = []
        # for d in self.obsTreeModel.deltas():
        for i in range(self.obsTreeModel.invisibleRootItem().rowCount()):
            survey = self.obsTreeModel.invisibleRootItem().child(i)
            for d in survey.deltas:
                stations.append(d.sta1)
                stations.append(d.sta2)
        station_list = list(set(stations))
        station_list.sort()
        station = AddDatumFromList.add_datum(station_list)
        if station:
            d = Datum(str(station))
            survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
            survey.datums.append(d)
            self.tab_adjust.datum_view.model().sourceModel().init_data(survey.datums)
            logging.info("Datum station added: {}".format(station))
            self.set_window_title_asterisk()

    def dialog_adjustment_properties(self):
        """
        Opens PyQt dialog to set adjustment options
        """
        survey = self.obsTreeModel.itemFromIndex(self.index_current_survey)
        adjust_options = AdjustOptions(
            survey.__str__(), survey.adjustment.adjustmentoptions, parent=self
        )
        if adjust_options.exec_():
            if adjust_options.surveys_to_update == "single":
                # Not sure why the deepcopy is necessary. Without it, all of the
                # survey.adjustmentoptions reference the same object.
                ao = copy.deepcopy(adjust_options.ao)
                survey.adjustment.adjustmentoptions = ao
                self.set_adj_sd(survey, adjust_options.ao)
            elif adjust_options.surveys_to_update == "all":
                for i in range(self.obsTreeModel.invisibleRootItem().rowCount()):
                    ao = copy.deepcopy(adjust_options.ao)
                    survey = self.obsTreeModel.invisibleRootItem().child(i)
                    survey.adjustment.adjustmentoptions = ao
                    self.set_adj_sd(survey, adjust_options.ao)
            self.set_window_title_asterisk()

    def dialog_tide_correction(self):
        """
        Opens PyQt dialog to specify correction type
        """
        tide_correction_dialog = TideCorrectionDialog()
        tide_correction_dialog.msg.exec_()
        correction_type = tide_correction_dialog.correction_type
        if correction_type == "Cancel":
            return
        elif correction_type == "Meter-supplied":
            tide_correction_meter(self)
        elif correction_type == "Agnew":
            lat, lon, elev = [], [], []
            # Get mean coordinates as default position
            for i in range(self.obsTreeModel.invisibleRootItem().rowCount()):
                survey = self.obsTreeModel.invisibleRootItem().child(i)
                for ii in range(survey.rowCount()):
                    loop = survey.child(ii)
                    for iii in range(loop.rowCount()):
                        station = loop.child(iii)
                        lat += station.lat
                        lon += station.long
                        elev += station.elev
            try:
                lat.remove(0.0)
            except ValueError:
                pass
            try:
                lon.remove(0.0)
            except ValueError:
                pass
            try:
                elev.remove(0.0)
            except ValueError:
                pass

            lat = np.mean(list(lat))
            lon = np.mean(list(lon))
            elev = np.mean(list(elev))

            tc = TideCoordinatesDialog(lat, lon, elev)
            if tc.exec_():
                tide_correction_agnew(
                    self,
                    float(tc.lat.text()),
                    float(tc.lon.text()),
                    float(tc.elev.text()),
                )
            self.adjust_update_required()
        self.update_data_tab()

    def dialog_station_coordinates(self):
        """
        Shows station coordinates dialog.
        """
        if not self.obsTreeModel.station_coords:
            init_station_coords_dict(self.obsTreeModel)
        coordinates_dialog = CoordinatesTable(self.obsTreeModel.station_coords)
        accept = coordinates_dialog.exec_()
        if accept == 1:
            self.obsTreeModel.station_coords = coordinates_dialog.coords()

    def load_station_coordinates(self):
        """
        Load station coordinates from file
        """
        # TODO: load station coordinates dialog
        return

    def closeEvent(self, event):
        self.settings.sync()
        if self.windowTitle()[-1] == "*":
            quit_msg = (
                "The workspace isn't saved. Are you sure you want to exit the program?"
            )
            reply = MessageBox.question("Message", quit_msg)

            if reply == QtWidgets.QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

    def show_help(self):
        """
        Shows compiled help file created using Dr. Explain in default browser
        """
        help_path = os.path.join(self.path_install, "docs", "index.htm")
        webbrowser.open(help_path)

    def dialog_about(self):
        if hasattr(self, "commit"):
            AboutDialog(self.commit)
        else:
            AboutDialog("")

    def check_for_updates(self, show_uptodate_msg):
        """
        Check if the usgs root repo is at the same commit as this installation
        Parameters
        ----------
        show_uptodate_msg : bool
            Whether to display a msg if no updates found
        Returns
        -------
        bool
            whether or not to start GSadjust (True=yes)
        """
        try:
            gitpath = (
                os.path.dirname(self.path_install) + "\\gsadjust-env\\Lib\\mingw64\\bin"
            )
            os.environ["PATH"] += os.pathsep + gitpath
            from git import Repo

            logging.info("Checking for updates")
            repo = Repo(self.path_install)
            if not repo.active_branch.name == "master":
                return True
            fetch = [r for r in repo.remotes if r.name == "origin"][0].fetch()
            master = [f for f in fetch if f.name == "origin/master"][0]
            for f in fetch:
                logging.info("Git fetched: {}".format(f))
            self.commit = str(repo.head.commit)[:5]
            if repo.head.commit != master.commit:
                msg = (
                    "An update is available for GSadjust.\nWould you like to install"
                    " now?"
                )
                confirm = MessageBox.question(
                    "Update Available",
                    msg,
                )
                if confirm == QtWidgets.QMessageBox.Yes:
                    return self.update_from_github()
            elif show_uptodate_msg:
                logging.info("Git checked, GSadjust is up to date.")
                msg = "GSadjust is up to date."
                MessageBox.information("No Update Needed", msg)
                return True
            return True

        except BaseException as e:
            logging.info("Git update failed: {}".format(e))
            if show_uptodate_msg:
                msg = "Problem Encountered Updating from GitHub\n\nError Message:\n"
                msg += str(e)
                MessageBox.information("Update results", msg)
            return True  # Update didn't work, start GSadjust anyway

    def update_from_github(self):
        """
        Merge the latest version into the local repo
        """

        try:
            from git import Repo

            repo = Repo(
                self.path_install,
            )
            fetch = [r for r in repo.remotes if r.name == "origin"][0].fetch()
            master = [f for f in fetch if f.name == "origin/master"][0]
            repo.git.reset("--hard")
            repo.git.merge(master.name)
            logging.info("Git update successful")
            msg = (
                "Updated successfully downloaded from GitHub. Please\nrestart GSadjust."
            )
            MessageBox.information("Update results", msg)
            return False  # Don't launch GSadjust, need to restart to install updates

        except BaseException as e:
            logging.info("Git update failed: {}".format(e))
            msg = (
                "Problem Encountered Updating from GitHub\n\n"
                "Please upgrade to the latest release by reinstalling the "
                "application from GitHub "
                "\n(https://github.com/jkennedy-usgs/sgp-gsadjust/releases)\n\n"
                "Error Message:\n"
            )
            msg += str(e)
            MessageBox.information("Update results", msg)
            return True  # Update didn't work, launch anyway

    def init_settings(self):
        """
        If it's the first time a user has opened GSadjust, populate the
        directories/files
        """
        if self.settings.value("current_dir") is None:
            self.settings.setValue(
                "current_dir", os.path.join(os.getcwd(), "test_data")
            )
        if self.settings.value("abs_g_path") is None:
            self.settings.setValue("abs_g_path", os.getcwd())
        if self.settings.value("tide_path") is None:
            self.settings.setValue("tide_path", os.getcwd())
        if self.settings.value("delta_table_column_widths") is None:
            self.settings.setValue("delta_table_column_widths", None)
        if self.settings.value("datum_table_column_widths") is None:
            self.settings.setValue("datum_table_column_widths", None)
        if self.settings.value("results_table_column_widths") is None:
            self.settings.setValue("results_table_column_widths", None)
        if self.settings.value("data_table_column_widths") is None:
            self.settings.setValue("data_table_column_widths", None)


class BoldDelegate(QtWidgets.QStyledItemDelegate):
    """
    Makes selected item in tree view bold.
    See     http://www.qtcentre.org/threads/61716-Set-the-color-of-a-row-in-a-qtreeview
    """

    def paint(self, painter, option, index):
        m = index.model().itemFromIndex(index)
        # decide here if item should be bold and set font weight to bold if needed
        if not hasattr(m, "fontweight"):
            option.font.setWeight(QtGui.QFont.Normal)
            m.cellcolor = Qt.white
        else:
            option.font.setWeight(m.fontweight)
        painter.fillRect(option.rect, m.cellcolor)
        QtWidgets.QStyledItemDelegate.paint(self, painter, option, index)


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Sends exceptions to log file
    """
    # KeyboardInterrupt is a special case, don't raise the error dialog when it occurs.
    if issubclass(exc_type, KeyboardInterrupt):
        if QtWidgets.qApp:
            QtWidgets.qApp.quit()
        return

    filename, line, dummy, dummy = traceback.extract_tb(exc_traceback).pop()
    filename = os.path.basename(filename)
    error = "{}: {}".format(exc_type.__name__, exc_value)
    logging.error(error + " at line {:d} of file {}".format(line, filename))


DEBUG = False


def except_hook2(cls, exception, traceback):
    logging.exception("%s %s", cls, exception, exc_info=True)
    sys.__excepthook__(cls, exception, traceback)


def main():
    # QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = QtWidgets.QApplication(sys.argv)

    # Needed to show icon in Windows Taskbar:
    import ctypes
    myappid = "usgs.sgp.gsadjust.1.0"  # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    # start log file
    fn = os.path.join(
        os.path.dirname(__file__),'..','..',
        "GSadjustLog_{}.txt".format(time.strftime("%Y%m%d-%H%M")),
    )
    try:
        logging.basicConfig(
            filename=fn, format="%(levelname)s:%(message)s", level=logging.INFO
        )
    except PermissionError:
        MessageBox.warning(
            "GSadjust error",
            "Please install GSadjust somewhere where admin rights are not required.",
        )
    splash_pix = QtGui.QPixmap(":/icons/Splash.png")
    splash = QtWidgets.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    ex = MainProg(splash=splash)
    app.processEvents()
    splash.finish(ex)
    app.setWindowIcon(QtGui.QIcon(":/icons/app.ico"))
    if not DEBUG:
        if ex.check_for_updates(False):
            ex.showMaximized()
            app.processEvents()
            app.exec_()
        else:
            ex.close()
    else:
        ex.showMaximized()
        app.processEvents()
        sys.exit(app.exec_())
