"""
models/delta.py
===============

Delta model for GSadjust.
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
import logging

from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, QVariant

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

    """
    _headers = {
        DELTA_STATION1: "From",
        DELTA_STATION2: "To",
        LOOP: "Loop",
        DELTA_TIME: "Time",
        DELTA_G: "Delta-g",
        DELTA_DRIFTCORR: "Drift corr.",
        DELTA_SD: "Std. dev.",
        DELTA_ADJ_SD: "SD for adj.",
        DELTA_RESIDUAL: "Residual",
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

                if delta.type == "normal":
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
                        return False
                elif delta.type == "assigned":
                    if column == DELTA_G:
                        brush = QtGui.QBrush(Qt.red)
                return brush
            if role == Qt.DisplayRole or role == Qt.EditRole:

                def delta_station_loop():
                    if delta.loop is None:
                        if type(delta.station2) == list:
                            return delta.station2[0].loop
                        else:
                            return "NA"
                    else:
                        return delta.loop

                def format_datetime(date):
                    return dt.datetime.utcfromtimestamp(date * 86400.0).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )

                def get_sd():
                    # Temporarily store previous value for logging, if the value
                    # is changed
                    delta._edited_sd = delta.adj_sd
                    return delta.adj_sd

                def get_g():
                    delta.edited_dg = delta.dg
                    if delta.type == "assigned":
                        return delta.assigned_dg
                    else:
                        return delta.dg

                def get_driftcorr():
                    if delta.driftcorr == "Roman" or delta.driftcorr == "Adj.":
                        return delta.driftcorr, ""
                    else:
                        return delta.driftcorr, "0.1f"

                fn, *args = {
                    DELTA_STATION1: (str, delta.sta1),
                    DELTA_STATION2: (str, delta.sta2),
                    LOOP: (str, delta_station_loop()),
                    DELTA_TIME: (format_datetime, delta.time()),
                    DELTA_G: (format, get_g(), "0.1f"),
                    DELTA_DRIFTCORR: (format, *get_driftcorr()),
                    DELTA_SD: (format, delta.sd, "0.1f"),
                    DELTA_ADJ_SD: (format, get_sd(), "0.1f"),
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
                try:
                    delta = self._data[index.row()]
                    column = index.column()
                    if len(str(value)) > 0:
                        if column == DELTA_ADJ_SD:
                            delta.adj_sd = float(value)
                            logging.info(
                                "delta {}, adj_sd changed from {} to {}".format(
                                    delta, delta._edited_sd, value
                                )
                            )
                        if column == DELTA_G:
                            if delta.type == "list":
                                self.tried_to_update_list_delta.emit()
                                return False
                            else:
                                delta.type = "assigned"
                                delta.assigned_dg = float(value)
                                logging.info(
                                    "delta {}, g changed from {} to {}".format(
                                        delta, delta.edited_dg, value
                                    )
                                )
                        self.dataChanged.emit(index, index)
                        self.signal_adjust_update_required.emit()
                    return True
                except ValueError:
                    return False
        if role == Qt.UserRole:
            self._data[index.row()] = value

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

    def init_data(self, data):
        self.beginResetModel()
        self._data = data
        self.endResetModel()
        self.layoutChanged.emit()  # Refresh whole view.


class NoCheckDeltaTableModel(DeltaTableModel):
    _is_checkable = False
