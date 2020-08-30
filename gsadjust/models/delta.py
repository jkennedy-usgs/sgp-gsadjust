#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
pyqt_modules.py
===============

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

import jsons
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, QVariant

from ..data import AdjustmentOptions, Datum

# Constants for column headers
(
    DELTA_STATION1,
    DELTA_STATION2,
    LOOP,
    DELTA_TIME,
    DELTA_G,
    DELTA_DRIFTCORR,
    DELTA_SD,
    DELTA_ADJ_SD,
    DELTA_RESIDUAL,
) = range(9)


class DeltaTableModel(QtCore.QAbstractTableModel):
    """
    Model that stores delta-g's used in network adjustment. Used on the network
    adjustment tab and as the Roman average table (bottom right table when using
    Roman method drift correction).

    There is one DeltaTableModel per survey. The delta-g's that populate the
    model depend on the drift method used; if the Roman method is used, the
    delta-g's are the average of the individual delta-g's in the RomanTableModel.
    """

    _headers = {
        DELTA_STATION1: 'From',
        DELTA_STATION2: 'To',
        LOOP: 'Loop',
        DELTA_TIME: 'Time',
        DELTA_G: 'Delta-g',
        DELTA_DRIFTCORR: 'Drift corr.',
        DELTA_SD: 'Std. dev.',
        DELTA_ADJ_SD: 'SD for adj.',
        DELTA_RESIDUAL: 'Residual',
    }

    _is_checkable = True

    def __init__(self):
        super(DeltaTableModel, self).__init__(parent=None)
        self._data = []

    tried_to_update_list_delta = QtCore.pyqtSignal()
    signal_adjust_update_required = QtCore.pyqtSignal()

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.TextAlignmentRole:
            if orientation == Qt.Horizontal:
                return QVariant(int(Qt.AlignLeft | Qt.AlignVCenter))
            return QVariant(int(Qt.AlignRight | Qt.AlignVCenter))

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, delta, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self._data.append(delta)
        self.endInsertRows()

    def rowCount(self, parent=None):
        if parent is None:
            return len(self._data)
        elif parent.isValid():
            return 0
        return len(self._data)

    def columnCount(self, parent=None):
        if parent is None:
            return 9
        elif parent.isValid():
            return 0
        return 9

    def data(self, index, role):
        if index.isValid():
            delta = self._data[index.row()]
            column = index.column()
            if role == Qt.ForegroundRole:
                brush = QtGui.QBrush(Qt.black)
                if delta.type == 'normal':
                    try:
                        if (
                            delta.station1.data(role=Qt.CheckStateRole) == 2
                            and delta.station2.data(role=Qt.CheckStateRole) == 2
                        ):
                            brush = QtGui.QBrush(Qt.black)
                        else:
                            if delta.station1.data(role=Qt.CheckStateRole) == 2:
                                if column == 0:
                                    brush = QtGui.QBrush(Qt.darkGray)
                                else:
                                    brush = QtGui.QBrush(Qt.lightGray)
                            elif delta.station2.data(role=Qt.CheckStateRole) == 2:
                                if column == 1:
                                    brush = QtGui.QBrush(Qt.darkGray)
                                else:
                                    brush = QtGui.QBrush(Qt.lightGray)
                            else:
                                brush = QtGui.QBrush(Qt.lightGray)
                    except:
                        catch = 1
                elif delta.type == 'assigned':
                    if column == DELTA_G:
                        brush = QtGui.QBrush(Qt.red)
                return brush
            if role == Qt.DisplayRole:

                def delta_station_loop():
                    if delta.loop is None:
                        if type(delta.station2) == list:
                            return delta.station2[0].loop
                        else:
                            return "NA"
                    else:
                        return delta.loop

                def format_datetime(date):
                    return dt.datetime.utcfromtimestamp(
                        date * 86400.0
                    ).strftime('%Y-%m-%d %H:%M:%S')

                fn, *args = {
                    DELTA_STATION1: (str, delta.sta1),
                    DELTA_STATION2: (str, delta.sta2),
                    LOOP: (str, delta_station_loop()),
                    DELTA_TIME: (format_datetime, delta.time()),
                    DELTA_G: (format, delta.dg, "0.1f"),
                    DELTA_DRIFTCORR: (format, delta.driftcorr, "0.1f"),
                    DELTA_SD: (format, delta.sd, "0.1f"),
                    DELTA_ADJ_SD: (format, delta.sd_for_adjustment, "0.1f"),
                    DELTA_RESIDUAL: (format, delta.residual, "0.1f"),
                }.get(column, (str, "NA"))

                return fn(*args)

            elif role == Qt.CheckStateRole:
                if self._is_checkable and index.column() == 0:
                    return self.checkState(delta)
                # fall out, return None for all other cases.

            elif role == Qt.UserRole:
                return delta

    def setData(self, index, value, role):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        # self.layoutAboutToBeChanged.emit()
        if role == Qt.CheckStateRole and index.column() == 0:
            delta = self._data[index.row()]
            if value == Qt.Checked:
                delta.checked = 2
            elif value == Qt.Unchecked:
                delta.checked = 0
            self.signal_adjust_update_required.emit()
            return True
        if role == Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                delta = self._data[index.row()]
                column = index.column()
                if len(str(value)) > 0:
                    if column == DELTA_ADJ_SD:
                        delta.adj_sd = float(value)
                    if column == DELTA_G:
                        if delta.type == 'list':
                            self.tried_to_update_list_delta.emit()
                            return False
                        else:
                            delta.type = 'assigned'
                            delta.assigned_dg = float(value)
                    self.dataChanged.emit(index, index)
                    self.signal_adjust_update_required.emit()
                return True
        if role == Qt.UserRole:
            self._data[index.row()] = value

        # self.layoutChanged.emit()

    def clearDeltas(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        for delta in self._data:
            delta.residual = -999
        self._data = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        self.layoutChanged.emit()
        return QVariant()

    def updateTable(self):
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def flags(self, QModelIndex):
        if not QModelIndex.isValid():
            return Qt.NoItemFlags
        return (
            Qt.ItemIsUserCheckable
            | Qt.ItemIsEnabled
            | Qt.ItemIsEditable
            | Qt.ItemIsSelectable
        )

    def checkState(self, delta):
        """
        By default, everything is checked. If keepdata property from the
        ChannelList object is 0, it is unchecked
        """
        if delta.checked == 0:
            return Qt.Unchecked
        else:
            return Qt.Checked

    @classmethod
    def from_json(cls, data):
        """
        When loading a workspace, repopulate PyQt models
        """
        temp = cls(data['name'])

        deltas = []
        for delta in data['deltas']:
            from types import SimpleNamespace

            sd = SimpleNamespace()
            sd.adj_sd = delta['adj_sd']
            sd.checked = delta['checked']
            sd.driftcorr = delta['driftcorr']
            sd.loop = delta['loop']
            sd.ls_drift = delta['ls_drift']
            sd.type = delta['type']
            try:
                sd.assigned_dg = delta['assigned_dg']
            except:
                pass
            # d = data_objects.SimpleDelta(sd)
            try:
                sd.sta1 = delta['sta1']
                sd.sta2 = delta['sta2']
            except KeyError as e:
                # Raised if delta type is 'assigned'
                pass
            deltas.append(sd)

        temp.init_data(deltas)

        for datum in data['datums']:
            d = Datum(datum['station'])
            d.__dict__ = datum
            temp.datum_model.insertRows(d, 0)

        ao = AdjustmentOptions()
        ao.__dict__ = data['adjoptions']
        temp.adjustment.adjustmentoptions = ao

        if 'checked' in data:
            temp.setCheckState(data['checked'])
        return temp

    def to_json(self):
        # Normal delta
        json_data = []
        for delta in self._data:
            if delta.type == 'normal' or delta.type == 'assigned':
                sta1 = delta.station1.key
                sta2 = delta.station2.key
            elif delta.type == 'list':
                sta1 = None
                sta2 = []
                for threepoint in delta.station2:
                    stations = [
                        threepoint.station1.key,
                        threepoint.station2[0].key,
                        threepoint.station2[1].key,
                    ]
                    sta2.append(stations)
            adj_sd = delta.adj_sd
            delta_type = delta.type
            ls_drift = delta.ls_drift
            if delta.loop is None:
                if isinstance(delta.station2, list):
                    loop = delta.station2[0].loop
                else:
                    loop = "NA"
            else:
                loop = delta.loop
            driftcorr = delta.driftcorr
            checked = delta.checked
            temp_delta = {
                'sta1': sta1,
                'sta2': sta2,
                'adj_sd': adj_sd,
                'ls_drift': ls_drift,
                'driftcorr': driftcorr,
                'checked': checked,
                'loop': loop,
                'type': delta_type,
            }
            json_data.append(temp_delta)
        return json_data

    def init_data(self, data):
        self._data = data
        self.layoutChanged.emit()  # Refresh whole view.


class NoCheckDeltaTableModel(DeltaTableModel):
    _is_checkable = False


def delta_serializer(obj, cls, **kwargs):
    """
    Handle serialization of DeltaTableModel via .to_json() method.
    """
    return obj.to_json()


jsons.set_serializer(delta_serializer, DeltaTableModel)
