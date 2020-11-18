#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
menus.py
===============

Menus for GSadjust
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
from PyQt5 import QtGui, QtWidgets


class MENU_STATE:
    UNINIT = -2
    INIT = -1
    CLEAR = 0
    ACTIVE_WORKSPACE = 1
    NO_ACTIVE_WORKSPACE = 12
    AT_LEAST_ONE_SURVEY = 2
    CALCULATE_CHANGE = 3
    MORE_THAN_ONE_SURVEY = 4
    SURVEY_HAS_DELTAS = 5
    SURVEY_HAS_NO_DELTAS = 6
    SURVEY_HAS_RESULTS = 7
    SURVEY_HAS_NO_RESULTS = 8
    # Context Menus
    OBS_TREE_STATION = 9
    OBS_TREE_LOOP = 10
    OBS_TREE_SURVEY = 13
    DELTA_MODEL = 11
    MULTIPLE_STATION = 14
    UNENABLE_ALL = 15


_ENABLED_MENUS = {
    MENU_STATE.UNINIT: [
        ('mnFileAppendLoop', False),
        ('mnFileAppendSurvey', False),
        ('mnFileAppendWorkspace', False),
        ('mnFileSaveWorkspaceAs', False),
        ('mnDivideLoopByTime', False),
        ('mnDivideLoopByHeight', False),
        ('mnEditAddTareDialog', False),
        ('mnEditVerticalGradientIntervalAction', False),
        ('mnEditVerticalGradientWriteAction', False),
        ('mnEditTideCorrection', False),
        ('mnEditCorrectRecordedTimeAction', False),
        ('mnEditShowCoordinates', False),
        ('mnEditLoadCoordinates', False),
        ('mnAdjUpdateDeltas', False),
        ('mnAdjUpdateDeltasCurrentSurvey', False),
        ('mnAdjUpdateDeltasCurrentLoop', False),
        ('mnAdjClearDeltaTable', False),
        ('mnAdjClearDatumTable', False),
        ('mnAdjUpdateDeltas', False),
        ('mnAdjUpdateDeltasCurrentSurvey', False),
        ('mnAdjOptions', False),
        ('mnAdjAdjustCurrent', False),
        ('mnAdjAdjust', False),
        ('mnAdjClearDeltaTable', False),
        ('mnAdjAddDatum', False),
        ('mnAdjImportAbsSimple', False),
        ('mnAdjImportAbsFull', False),
        ('mnAdjImportAbsDatabase', False),
        ('mnAdjClearDatumTable', False),
        ('mnAdjPlotHist', False),
        ('mnAdjPlotCompareDatum', False),
        ('mnAdjPlotObservedAdjustedAbs', False),
        ('mnToolsNGCircular', False),
        ('mnToolsNGMap', False),
        ('mnToolsComputeGravityChangeAction', False),
        ('mnToolsWriteMetadataText', False),
        ('mnToolsWriteTabularOutput', False),
        ('mnToolsWriteSummary', False),
    ],
    MENU_STATE.INIT: [
        ('mnFileAppendLoop', True),
        ('mnFileAppendSurvey', True),
        ('mnFileAppendWorkspace', True),
        ('mnFileSaveWorkspaceAs', True),
        ('mnEditTideCorrection', True),
        ('mnEditCorrectRecordedTimeAction', True),
        ('mnEditShowCoordinates', True),
        ('mnEditLoadCoordinates', True),
        ('mnAdjClearDatumTable', True),
        ('mnAdjUpdateDeltas', True),
        ('mnAdjUpdateDeltasCurrentSurvey', True),
        ('mnAdjOptions', True),
        ('mnAdjAdjustCurrent', True),
        ('mnAdjAdjust', True),
        ('mnAdjUpdateDeltas', True),
        ('mnAdjUpdateDeltasCurrentSurvey', True),
        ('mnAdjUpdateDeltasCurrentLoop', True),
        ('mnAdjClearDeltaTable', True),
        ('mnAdjAddDatum', True),
        ('mnAdjImportAbsSimple', True),
        ('mnAdjImportAbsFull', True),
        ('mnAdjImportAbsDatabase', True),
        ('mnAdjClearDatumTable', True),
    ],
    MENU_STATE.ACTIVE_WORKSPACE: [('mnFileSaveWorkspace', True),],
    MENU_STATE.NO_ACTIVE_WORKSPACE: [('mnFileSaveWorkspace', False),],
    MENU_STATE.SURVEY_HAS_DELTAS: [
        ('mnToolsNGCircular', True),
        ('mnToolsNGMap', True),
        ('mnToolsComputeGravityChangeAction', True),
    ],
    MENU_STATE.SURVEY_HAS_NO_DELTAS: [
        ('mnToolsNGCircular', False),
        ('mnToolsNGMap', False),
        ('mnToolsComputeGravityChangeAction', False),
    ],
    MENU_STATE.SURVEY_HAS_RESULTS: [
        ('mnAdjPlotHist', True),
        ('mnAdjPlotCompareDatum', True),
        ('mnAdjPlotObservedAdjustedAbs', True),
        ('mnToolsWriteSummary', True),
        ('mnToolsWriteTabularOutput', True),
        ('mnToolsWriteMetadataText', True),
        ('mnAdjPlotObservedAdjustedAbs', True),
        ('mnAdjUpdateSD', True),
        ('mnToolsShowCalCoeffTimeSeries', True),
    ],
    MENU_STATE.SURVEY_HAS_NO_RESULTS: [
        ('mnAdjPlotHist', False),
        ('mnAdjPlotCompareDatum', False),
        ('mnAdjPlotObservedAdjustedAbs', False),
        ('mnToolsWriteSummary', False),
        ('mnToolsWriteTabularOutput', False),
        ('mnToolsWriteMetadataText', False),
        ('mnAdjPlotObservedAdjustedAbs', False),
        ('mnAdjUpdateSD', False),
    ],
    MENU_STATE.CALCULATE_CHANGE: [('mnToolsComputeGravityChangeAction', True),],
    MENU_STATE.AT_LEAST_ONE_SURVEY: [
        ('mnAdjAdjust', True),
        ('mnAdjAdjustCurrent', True),
        ('mnDivideLoopByHeight', True),
        ('mnDivideLoopByTime', True),
        ('mnEditAddTareDialog', True),
        ('mnEditVerticalGradientIntervalAction', True),
        ('mnEditVerticalGradientWriteAction', True),
        ('mnToolsComputeGravityChangeAction', False),
    ],
    MENU_STATE.MORE_THAN_ONE_SURVEY: [
        ('mnAdjPlotObservedAdjustedAbs', True),
        ('mnAdjPlotCompareDatum', True),
        ('mnEditVerticalGradientWriteAction', False),
    ],
    # Right-click context menus below here
    MENU_STATE.OBS_TREE_STATION: [
        ('mnDeleteSurvey', False),
        ('mnDeleteLoop', False),
        ('mnLoopProperties', False),
        ('mnDeleteStation', True),
        ('mnStationDuplicate', True),
        ('mnRename', True),
        ('mnDataNewLoop', True),
        ('mnLoopAnimate', False),
        ('mnVerticalGradientWriteAction', False),
    ],
    MENU_STATE.OBS_TREE_LOOP: [
        ('mnDeleteSurvey', False),
        ('mnDeleteLoop', True),
        ('mnLoopProperties', True),
        ('mnDeleteStation', False),
        ('mnStationDuplicate', False),
        ('mnRename', True),
        ('mnDataNewLoop', False),
        ('mnLoopAnimate', True),
        ('mnVerticalGradientWriteAction', True),
    ],
    MENU_STATE.OBS_TREE_SURVEY: [
        ('mnDeleteSurvey', True),
        ('mnDeleteLoop', False),
        ('mnLoopProperties', False),
        ('mnDeleteStation', False),
        ('mnStationDuplicate', False),
        ('mnRename', True),
        ('mnDataNewLoop', False),
        ('mnLoopAnimate', False),
        ('mnVerticalGradientWriteAction', False),
    ],
    MENU_STATE.UNENABLE_ALL: [
        ('mnDeleteSurvey', False),
        ('mnDeleteLoop', False),
        ('mnLoopProperties', False),
        ('mnDeleteStation', False),
        ('mnStationDuplicate', False),
        ('mnRename', False),
        ('mnDataNewLoop', False),
        ('mnLoopAnimate', False),
        ('mnVerticalGradientWriteAction', False),
    ],
    MENU_STATE.MULTIPLE_STATION: [
        ('mnDeleteSurvey', False),
        ('mnDeleteLoop', False),
        ('mnLoopProperties', False),
        ('mnDeleteStation', True),
        ('mnStationDuplicate', False),
        ('mnRename', False),
        ('mnDataNewLoop', True),
        ('mnLoopAnimate', False),
        ('mnVerticalGradientWriteAction', False),
    ],
    MENU_STATE.DELTA_MODEL: [
        ('mnAdjClearDeltaTable', True),
        ('mnAdjAdjust', True),
        ('mnAdjAdjustCurrent', True),
    ],
}


class Menus:
    def __init__(self, mainProg):
        self.mainProg = mainProg
        """
        Menu creation
        """
        #######################################################################
        # File Menu:
        #######################################################################
        self.mnFile = self.mainProg.menuBar().addMenu("&File")
        # for file opening
        self.mnFileOpenScintrexFile = self.create_action(
            "&Load raw CG-3/CG-5 data...",
            shortcut="Ctrl+5",
            slot=lambda meter='CG5': self.mainProg.open_file_dialog('CG5'),
            tip="Open data file",
            enabled=True,
        )
        self.mnFileOpenBurrisFile = self.create_action(
            "Load raw &Burris data...",
            shortcut="Ctrl+b",
            slot=lambda meter='Burris': self.mainProg.open_file_dialog('Burris'),
            tip="Open data file",
            enabled=True,
        )
        self.mnFileOpenCG6File = self.create_action(
            "Load raw CG-6 data...",
            shortcut="Ctrl+6",
            slot=lambda meter='CG6': self.mainProg.open_file_dialog('CG6'),
            tip="Open data file",
            enabled=True,
        )
        self.mnFileOpenCG6FileTsoft = self.create_action(
            "Load raw CG-6 data (Tsoft format)...",
            slot=lambda meter='CG6Tsoft': self.mainProg.open_file_dialog('CG6Tsoft'),
            tip="Open data file",
            enabled=True,
        )
        self.mnFileOpenCsvFile = self.create_action(
            "Import csv data...",
            shortcut="Ctrl+7",
            slot=lambda meter='csv': self.mainProg.open_file_dialog('csv'),
            tip="Open comma-separated value data file",
            enabled=True,
        )
        self.mnFileAppendSurvey = self.create_action(
            "Append survey to campaign...",
            slot=lambda meter='choose': self.mainProg.open_file_dialog('survey'),
            enabled=False,
        )
        self.mnFileAppendLoop = self.create_action(
            "Append loop to current survey...",
            slot=lambda meter='choose': self.mainProg.open_file_dialog('loop'),
            enabled=False,
        )
        self.mnFileAppendWorkspace = self.create_action(
            "Append workspace...", slot=self.mainProg.workspace_append, enabled=False
        )
        self.mnFileClearWorkspace = self.create_action(
            "Clear workspace",
            slot=lambda clear=True: mainProg.workspace_clear(confirm=True),
            enabled=True,
        )
        self.mnFileSaveWorkspace = self.create_action(
            "Save workspace",
            shortcut="Ctrl+s",
            slot=self.mainProg.workspace_save,
            enabled=False,
        )
        self.mnFileSaveWorkspaceAs = self.create_action(
            "Save workspace as...",
            shortcut="F12",
            slot=self.mainProg.workspace_save_as,
            enabled=False,
        )
        self.mnFileLoadJSONAction = self.create_action(
            "&Open workspace...", slot=self.mainProg.workspace_open_getjson,
            shortcut="Ctrl+o"
        )
        self.mnFileExitAction = self.create_action(
            "&Exit", shortcut="Alt+F4", slot=self.mainProg.close, tip="Exit App"
        )

        # add actions to menu
        self.add_actions(
            self.mnFile,
            (
                self.mnFileOpenScintrexFile,
                self.mnFileOpenCG6File,
                self.mnFileOpenCG6FileTsoft,
                self.mnFileOpenBurrisFile,
                self.mnFileOpenCsvFile,
                None,
                self.mnFileAppendSurvey,
                self.mnFileAppendLoop,
                self.mnFileAppendWorkspace,
                None,
                self.mnFileSaveWorkspace,
                self.mnFileSaveWorkspaceAs,
                self.mnFileLoadJSONAction,
                self.mnFileClearWorkspace,
                None,
                self.mnFileExitAction,
            ),
        )

        #######################################################################
        # Edit Menu:
        #######################################################################
        self.mnEdit = self.mainProg.menuBar().addMenu("&Edit")
        self.mnEditTideCorrection = self.create_action(
            "&Tide correction...",
            shortcut="Ctrl+T",
            slot=self.mainProg.dialog_tide_correction,
            tip="Choose tide correction method",
            enabled=False,
        )
        self.mnEditOceanLoadingCorrectionAction = self.create_action(
            "&Ocean loading correction",
            shortcut="Ctrl+L",
            slot=self.mainProg.correction_ocean_loading,
            tip="Ocean loading corrections",
            enabled=False,
        )
        self.mnEditAtmosphericCorrectionAction = self.create_action(
            "&Atmospheric correction",
            shortcut="Ctrl+P",
            slot=self.mainProg.correction_atmospheric,
            tip="Apply atmospheric correction by loading a data file",
            enabled=False,
        )
        self.mnEditShowCoordinates = self.create_action(
            "Show station coordinates...",
            slot=self.mainProg.dialog_station_coordinates,
            enabled=False,
        )
        self.mnEditLoadCoordinates = self.create_action(
            "Load station coordinates",
            slot=self.mainProg.load_station_coordinates,
            enabled=False,
        )
        self.mnEditCorrectRecordedTimeAction = self.create_action(
            "&Correct recorded time...",
            slot=self.mainProg.correct_recorded_time,
            tip="Correct recorded time",
            enabled=False,
        )
        self.mnEditVerticalGradientIntervalAction = self.create_action(
            "&Vertical gradient interval...",
            slot=self.mainProg.set_vertical_gradient_interval,
            enabled=True,
        )
        self.mnEditVerticalGradientWriteAction = self.create_action(
            "Write vertical gradient file...",
            slot=self.mainProg.vertical_gradient_write,
            enabled=True,
        )
        self.mnEditAddTareDialog = self.create_action(
            "Add tare...", slot=self.mainProg.add_tare, enabled=True
        )
        self.mnDivideLoopByTime = self.create_action(
            "Divide selected loop by time...",
            slot=self.mainProg.divide_by_time,
            enabled=True,
        )
        self.mnDivideLoopByHeight = self.create_action(
            "Divide station by height...",
            slot=self.mainProg.divide_by_height,
            enabled=True,
        )
        # add actions to menu
        self.add_actions(
            self.mnEdit,
            (
                self.mnEditTideCorrection,
                self.mnEditCorrectRecordedTimeAction,
                None,
                self.mnEditShowCoordinates,
                # self.mnEditLoadCoordinates,
                None,
                self.mnDivideLoopByTime,
                self.mnDivideLoopByHeight,
                self.mnEditAddTareDialog,
                None,
                self.mnEditVerticalGradientIntervalAction,
            ),
        )

        #######################################################################
        # Adjust Menu
        #######################################################################
        self.mnAdj = self.mainProg.menuBar().addMenu("&Adjustment")
        self.mnAdjtype = self.mnAdj.addMenu('Adjustment type')
        actiongroup_adjustmenttype = QtWidgets.QActionGroup(self.mainProg)

        self.mnAdjPyLSQ = actiongroup_adjustmenttype.addAction(
            self.create_action('Numpy least squares', checkable=True)
        )
        self.mnAdjtype.addAction(self.mnAdjPyLSQ)
        self.mnAdjGravnet = actiongroup_adjustmenttype.addAction(
            self.create_action('Gravnet', checkable=True)
        )
        self.mnAdjtype.addAction(self.mnAdjGravnet)
        self.mnAdjPyLSQ.setChecked(True)

        self.mnAdjOptions = self.create_action(
            "Adjustment options...",
            slot=self.mainProg.dialog_adjustment_properties,
            enabled=True,
        )
        self.mnAdjAdjustCurrent = self.create_action(
            "&Adjust current survey",
            shortcut="Ctrl+2",
            slot=lambda: self.mainProg.adjust_network('current'),
            tip="Adjust current survey",
            enabled=False,
            icon=QtGui.QIcon(':/icons/ac.png'),
        )
        self.mnAdjAdjust = self.create_action(
            "&Adjust all surveys",
            shortcut="Ctrl+1",
            slot=lambda: self.mainProg.adjust_network('all'),
            tip="Adjust all surveys",
            enabled=False,
            icon=QtGui.QIcon(':/icons/aa.png'),
        )
        self.mnAdjUpdateDeltas = self.create_action(
            "&Populate delta table - all surveys",
            shortcut="Ctrl+A",
            slot=lambda: self.mainProg.populate_deltamodel('all'),
            enabled=False,
        )
        self.mnAdjUpdateDeltasCurrentSurvey = self.create_action(
            "Populate delta table - selected survey",
            slot=lambda: self.mainProg.populate_deltamodel('selectedsurvey'),
            enabled=False,
        )
        self.mnAdjUpdateDeltasCurrentLoop = self.create_action(
            "Populate delta table - selected loop(s)",
            slot=lambda: self.mainProg.populate_deltamodel('selectedloop'),
            enabled=False,
        )
        self.mnAdjClearDeltaTable = self.create_action(
            "Clear delta table", slot=self.mainProg.clear_delta_model, enabled=False
        )
        self.mnAdjAddDatum = self.create_action(
            "Add datum observation...",
            shortcut="Ctrl+D",
            slot=self.mainProg.dialog_add_datum,
            enabled=True,
        )
        self.mnAdjImportAbsSimple = self.create_action(
            "Import abs. g (simple)...",
            slot=self.mainProg.menu_import_abs_g_simple,
            enabled=True,
        )
        self.mnAdjImportAbsFull = self.create_action(
            "Import abs. g (complete)...",
            slot=self.mainProg.menu_import_abs_g_complete,
            enabled=True,
        )
        self.mnAdjImportAbsDatabase = self.create_action(
            "Import abs. g from project files...",
            slot=self.mainProg.dialog_import_abs_g_direct,
            enabled=True,
        )
        self.mnAdjClearDatumTable = self.create_action(
            "Clear datum table", slot=self.mainProg.clear_datum_model, enabled=False
        )

        self.mnAdjPlotHist = self.create_action(
            "Plot residual histogram",
            shortcut="Ctrl+H",
            slot=mainProg.plot_histogram,
            enabled=True,
        )
        self.mnAdjPlotCompareDatum = self.create_action(
            "Plot adjusted datum vs. measured",
            slot=self.mainProg.plot_datum_vs_adjusted,
            enabled=False,
        )
        self.mnAdjPlotObservedAdjustedAbs = self.create_action(
            "Plot adjusted datum vs. measured (time series)",
            slot=self.mainProg.plot_datum_comparison_timeseries,
            enabled=False,
        )
        self.mnAdjUpdateSD = self.create_action(
            "Scale std. dev. from results",
            slot=self.mainProg.update_SD_and_run_adjustment,
            enabled=False,
            icon=QtGui.QIcon(':/icons/ua.png'),
        )
        # self.mnAdjCompareAllDatum = self.create_action("Plot absolute vs. relative offset",
        #                  slot=self.mainProg.plot_compare_datum_all,
        #                  enabled=True)

        self.add_actions(
            self.mnAdj,
            (
                self.mnAdjOptions,
                self.mnAdjAdjustCurrent,
                self.mnAdjAdjust,
                self.mnAdjUpdateSD,
                None,
                self.mnAdjUpdateDeltas,
                self.mnAdjUpdateDeltasCurrentSurvey,
                self.mnAdjUpdateDeltasCurrentLoop,
                self.mnAdjClearDeltaTable,
                None,
                self.mnAdjAddDatum,
                self.mnAdjImportAbsSimple,
                self.mnAdjImportAbsFull,
                self.mnAdjImportAbsDatabase,
                self.mnAdjClearDatumTable,
                None,
                self.mnAdjPlotHist,
                self.mnAdjPlotCompareDatum,
                self.mnAdjPlotObservedAdjustedAbs,
            ),
        )

        #######################################################################
        # Tools Menu
        #######################################################################
        self.mnTools = self.mainProg.menuBar().addMenu('&Tools')
        self.mnToolsNGCircular = self.create_action(
            'Network graph - circular', slot=self.mainProg.plot_network_graph_circular
        )
        self.mnToolsNGMap = self.create_action(
            'Network graph - map view', slot=self.mainProg.plot_network_graph_map
        )
        self.mnToolsWriteMetadataText = self.create_action(
            'Write metadata text', slot=self.mainProg.write_metadata_text
        )
        self.mnToolsWriteTabularOutput = self.create_action(
            'Write tabular data', slot=self.mainProg.write_tabular_data
        )
        self.mnToolsWriteSummary = self.create_action(
            'Write adjustment summary', slot=self.mainProg.write_summary_text
        )
        self.mnToolsComputeGravityChangeAction = self.create_action(
            "&Compute gravity change",
            slot=self.mainProg.show_gravity_change_table,
            tip="Compute gravity change",
            enabled=False,
        )
        self.mnToolsShowCalCoeffTimeSeries = self.create_action(
            "Show calibration coefficients",
            slot=self.mainProg.show_cal_coeff,
            enabled=False,
        )

        self.add_actions(
            self.mnTools,
            (
                self.mnToolsNGCircular,
                self.mnToolsNGMap,
                None,
                self.mnToolsComputeGravityChangeAction,
                self.mnToolsShowCalCoeffTimeSeries,
                # self.mnToolsLOO,
                None,
                self.mnToolsWriteMetadataText,
                self.mnToolsWriteTabularOutput,
                self.mnToolsWriteSummary,
            ),
        )

        #######################################################################
        # Help Menu
        #######################################################################
        self.mnHelp = self.mainProg.menuBar().addMenu('&Help')
        self.mnHelpHelp = self.create_action('Help', slot=self.mainProg.show_help)
        self.mnHelpAbout = self.create_action('About', slot=self.mainProg.dialog_about)
        self.mnHelpLog = self.create_action('Log window', slot=self.mainProg.toggle_logview)
        self.mnHelpCheckForUpdates = self.create_action(
            'Check for updates...', slot=lambda: self.mainProg.check_for_updates(True)
        )

        self.add_actions(
            self.mnHelp,
            (self.mnHelpHelp, self.mnHelpAbout, self.mnHelpLog, None, self.mnHelpCheckForUpdates),
        )

    def create_action(
        self,
        text,
        slot=None,
        shortcut=None,
        icon=None,
        tip=None,
        checkable=False,
        signal="triggered()",
        enabled=True,
    ):
        """
        Simplify action creation
        """
        action = QtWidgets.QAction(text, self.mainProg)
        if icon is not None:
            action.setIcon(icon)
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        if enabled is not True:
            action.setEnabled(False)
        return action

    def add_actions(self, target, actions):
        """
        Add actions to a target
        """
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def set_state(self, *args):
        for state in args:
            if state not in _ENABLED_MENUS:
                raise Exception("Invalid menu state: " % state)

            for menu, enabled in _ENABLED_MENUS[state]:
                getattr(self, menu).setEnabled(enabled)
