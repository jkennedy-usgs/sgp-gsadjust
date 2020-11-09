#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
models/obstree
=======

PyQt models for GSadjust. Handles assembling input matrices for
network adjustment.
--------------------------------------------------------------------------------

NB: PyQt models follow the PyQt CamelCase naming convention. All other
methods/functions in GSadjust use PEP-8 lowercase_underscore convention.

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
import json
import logging

import jsons
from matplotlib.dates import date2num
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt

from .loop import ObsTreeLoop
from .station import ObsTreeStation
from .survey import ObsTreeSurvey

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0


class ObsTreeModel(QtGui.QStandardItemModel):
    """
    Tree model that shows station name, date, and average g value.

    The model is populated by appending stations to loops, and loops to surveys.
    The parent-child relationship is not explicitly stored.
    """

    # TODO: Implement drag and drop.
    def __init__(self):
        super(ObsTreeModel, self).__init__()
        self.setColumnCount(3)
        self.setHorizontalHeaderLabels(['Name', 'Date', 'g (\u00b5Gal)'])
        self.station_coords = None

    #
    # signal_refresh_view = QtCore.pyqtSignal()
    signal_name_changed = QtCore.pyqtSignal()

    def columnCount(self, QModelIndex_parent=None, *args, **kwargs):
        return 3

    def flags(self, QModelIndex):
        if not QModelIndex.isValid():
            return Qt.NoItemFlags
        return (
            Qt.ItemIsEnabled
            | Qt.ItemIsSelectable
            | Qt.ItemIsEditable
            | Qt.ItemIsUserCheckable
        )

    def data(self, index, role=Qt.DisplayRole):
        if index.model() is not None:
            column = index.column()
            if role == Qt.DisplayRole:

                if column > 0:
                    m = index.model().itemFromIndex(index.sibling(index.row(), 0))
                else:
                    m = index.model().itemFromIndex(index)
                try:  # Was getting "AttributeError: 'QStandardItem' object has no attribute 'display_column_map'"
                    # after deleting a survey
                    fn, *args = m.display_column_map.get(
                        column, (format_numeric_column, column)
                    )
                    return fn(*args)
                except AttributeError:
                    return ''

            elif role == Qt.CheckStateRole:
                if column == 0:
                    m = index.model().itemFromIndex(index)
                    return m.checkState()
            elif role == Qt.ToolTipRole:
                m = index.model().itemFromIndex(index)
                try:
                    return m.tooltip
                except AttributeError:
                    return ''

    def setData(self, index, value, role):
        if role == Qt.CheckStateRole and index.column() == 0:
            m = index.model().itemFromIndex(index)
            m.setCheckState(value)
            self.dataChanged.emit(index, index)
            return True

        if role == Qt.EditRole:
            if index.isValid() and index.column() == 0:
                item = index.model().itemFromIndex(index)
                return {
                    ObsTreeStation: self._handler_edit_ObsTreeStation,
                    ObsTreeLoop: self._handler_edit_ObsTreeLoop,
                    ObsTreeSurvey: self._handler_edit_ObsTreeSurvey,
                }[type(item)](item, value)

    def _handler_edit_ObsTreeStation(self, item, value):
        """
        Handle editing of ObsTreeStation objects (renaming).
        """
        old_name = item.station_name
        new_name = str(value)
        if new_name is not item.station_name:
            rename_type = rename_dialog(old_name, new_name)
            if rename_type == 'Loop':
                loop = item.parent()
                loop.rename(old_name, new_name)

            if rename_type == 'Survey':
                loop = item.parent()
                survey = loop.parent()
                survey.rename(old_name, new_name)

            if rename_type == 'Campaign':
                campaign = item.model().invisibleRootItem()
                for i in range(campaign.rowCount()):
                    survey = campaign.child(i, 0)
                    survey.rename(old_name, new_name)

            if rename_type == 'Station':
                item.station_name = new_name
                for i in range(len(item.station)):
                    item.station[i] = new_name

            logging.info(
                'Stations renamed from {} to {} in {}'.format(
                    old_name, new_name, rename_type
                )
            )
            self.signal_name_changed.emit()
        return True

    def _handler_edit_ObsTreeLoop(self, item, value):
        """
        Handle editing of ObsTreeLoop objects (renaming).
        """
        new_name = str(value)
        old_name = item.name
        logging.info('Loop renamed from {} to {}'.format(old_name, new_name))
        item.name = new_name
        return True

    def _handler_edit_ObsTreeSurvey(self, item, value):
        new_name = str(value)
        try:
            name_as_date = dt.datetime.strptime(new_name, '%Y-%m-%d')
            old_name = name_as_date
            logging.info('Loop renamed from {} to {}'.format(old_name, new_name))
            item.name = new_name
            return True
        except Exception as e:
            return False

    def checked_surveys(self):
        """
        Retrieves checked data from loop; used for calculating delta-gs in the
        "None" and "Continuous" drift options.
        :return: list of checked stations
        """
        data = []
        for i in range(self.invisibleRootItem().rowCount()):
            if self.invisibleRootItem().child(i).checkState() == 2:
                data.append(self.invisibleRootItem().child(i))
        return data

    def unique_meters(self):
        meters = []
        for survey in self.checked_surveys():
            meters += survey.unique_meters
        return list(set(meters))

    def get_cal_coeffs(self):
        cal_dict = {}
        for m in self.unique_meters():
            cal_dict[m] = []
        for survey in self.checked_surveys():
            if survey.adjustment.adjustmentoptions.cal_coeff:
                for k, v in survey.adjustment.adjustmentresults.cal_dic.items():
                    cal_dict[k].append((survey.name, v[0], v[1]))
            elif survey.adjustment.adjustmentoptions.specify_cal_coeff:
                for k, v in survey.adjustment.adjustmentoptions.meter_cal_dict.items():
                    cal_dict[k].append((survey.name, v, 0))
        return cal_dict

    def checkState(self, index):
        if index.column > 0:
            m = index.model().itemFromIndex(index.sibling(index.row(), 0))
        else:
            m = index.model().itemFromIndex(index)
        return m.checkState()

    def insertRows(self, station, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(index, position, position + rows - 1)
        self.appendRow(station)
        self.endInsertRows()

    def load_workspace(self, fname):
        """
        Load previously-save workspace. Need to recreate PyQt models.

        Importantly, there are two types of deltas: those on the Drift tab and
        those on the network adjustment tab. deltas are COPIED from one to
        another by discrete menu commands, it's not automatic.

        The workflow:
        1) Create obstreesurvey object so we have somewhere to store loop and
            station objects.
        2) deltas on the drift tab are created strictly from the station objects
            and specified options (drift correction method, etc.). We don't store
            a corresponding delta object in the saved workspace/json.
        3) deltas on the network adjustment tab have additional information that
            must be stored in the saved workspace/json (checked state, std. dev.
            for adj., etc.). When loading a workspace, this delta is INDEPENDENT of
            the delta stored in the delta table on the Drift tab. The Net Adj. tab
            delta table is re-created later, not here, because Roman-method deltas
            depend on the drift-tab deltas.
        :param fname:
        :return: (ObsTreeSurvey, delta_table, coords)
        """
        logging.info("Workspace loaded: " + fname)
        delta_tables, obstreesurveys = [], []
        coords, surveys = None, None
        with open(fname, "r") as f:
            data = jsons.load(json.load(f))
            if all(isinstance(x, ObsTreeSurvey) for x in data):
                surveys = data
            elif len(data) > 1:
                coords = data[1]
                surveys = data[0]
        # Populate PyQt objects. First, create obstreesurvey
        for survey in surveys:
            obstreesurvey = ObsTreeSurvey.from_json(survey)
            for loop in survey['loops']:
                obstreeloop = ObsTreeLoop.from_json(loop)
                for station in loop['stations']:
                    if (
                        'station_name' in station
                    ):  # Sometimes blank stations are generated, not sure why?
                        temp_station = tempStation(station)
                        obstreestation = ObsTreeStation(
                            temp_station,
                            temp_station.station_name,
                            temp_station.station_count,
                        )
                        if type(obstreestation.t[0]) == dt.datetime:
                            obstreestation.t = [date2num(i) for i in obstreestation.t]
                        obstreeloop.appendRow(
                            [
                                obstreestation,
                                QtGui.QStandardItem('0'),
                                QtGui.QStandardItem('0'),
                            ]
                        )
                # stations is provided by loop.stations(), it doesn't need to be stored separately in the ObsTreeLoop
                # object
                if hasattr(obstreeloop, 'stations'):
                    del obstreeloop.stations
                obstreesurvey.appendRow(
                    [obstreeloop, QtGui.QStandardItem('0'), QtGui.QStandardItem('0')]
                )
            obstreesurveys.append(obstreesurvey)
            delta_tables.append(obstreesurvey.deltas)

        return (obstreesurveys, delta_tables, coords)

    def datums(self):
        datum_list = []
        for i in range(self.rowCount()):
            obstreesurvey = self.itemFromIndex(self.index(i, 0))
            for ii in range(obstreesurvey.datum_model.rowCount()):
                idx = obstreesurvey.datum_model.index(ii, 0)
                datum = obstreesurvey.datum_model.data(idx, role=Qt.UserRole)
                datum_list.append(datum.station)
        return list(set(datum_list))

    def surveys(self):
        survey_list = []
        for i in range(self.rowCount()):
            obstreesurvey = self.itemFromIndex(self.index(i, 0))
            survey_list.append(obstreesurvey)
        return survey_list

    def deltas(self):
        delta_list = []
        for obstreesurvey in self.surveys():
            delta_list += obstreesurvey.delta_list()
        return delta_list

    def stations(self):
        station_list = []
        for obstreesurvey in self.surveys():
            for obstreeloop in obstreesurvey.loops():
                station_list += obstreeloop.stations()
        return list(set(station_list))

    def resetStationAsd(self):
        for survey in self.surveys():
            for loop in survey.loops():
                for station in loop.stations():
                    station.assigned_sd = None

    def save_workspace(self, fname):
        surveys = []
        for i in range(self.rowCount()):
            obstreesurvey = self.itemFromIndex(self.index(i, 0))
            surveys.append(obstreesurvey)
        workspace_data = [surveys, self.station_coords]
        if fname[-4:] != '.gsa':
            fname += '.gsa'
        with open(fname, "w") as f:
            json.dump(jsons.dump(workspace_data), f)
        logging.info('Saving JSON workspace to {}'.format(fname))
        return fname


class tempStation:
    def __init__(self, station):
        self.__dict__ = station


def rename_dialog(old_name, new_name):
    """
    Dialog called after renaming station in treeview.

    Gives the option to rename stations in the current loop, survey, or
    throughout the campaign.
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
    methods = {0: 'Campaign', 1: 'Survey', 2: 'Loop', 3: 'Station', 4: 'Cancel'}

    return methods[method]


def survey_serializer(obj, cls, **kwargs):
    """
    Handle serialization of ObsTreeSurvey via .to_json() method.
    """
    return obj.to_json()


def format_numeric_column(column):
    """
    Format fn for simple numeric columns.
    """
    return column + 1


jsons.set_serializer(survey_serializer, ObsTreeSurvey)
