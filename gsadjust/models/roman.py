"""
models/roman.py
===============

PyQt model for showing Roman-method gravity differences.
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
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, QVariant
import datetime as dt
# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0

ROMAN_FROM, ROMAN_TO, ROMAN_DELTA, TIME = range(4)
# Constants for column headers


class SamplesTableModel(QtCore.QAbstractTableModel):
    """
    Model to store the individual delta-g's calculated using the Roman-method.
    Shown on the drift tab.
    """

    _headers = {
        ROMAN_FROM: "From",
        ROMAN_TO: "To",
        ROMAN_DELTA: "Delta g",
        TIME: "Time",
    }

    def __init__(self):
        super(SamplesTableModel, self).__init__()
        self._data = []

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
        return len(self._data)

    def columnCount(self, parent=None):
        return 4

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            delta = self._data[index.row()]

            def format_datetime(date):
                return dt.datetime.utcfromtimestamp(date * 86400.0).strftime(
                    "%Y-%m-%d %H:%M:%S"
                )

            if role == Qt.DisplayRole:
                column = index.column()
                fn, *args = {
                    ROMAN_FROM: (str, delta.sta1),
                    ROMAN_TO: (str, delta.sta2),
                    ROMAN_DELTA: (format, delta.dg, "0.1f"),
                    TIME: (format_datetime, delta.time),
                }.get(column)
                return fn(*args)

            elif role == Qt.UserRole:
                # check status definition
                return delta

    def setData(self, index, value, role):
        if role == Qt.CheckStateRole and index.column() == 0:
            datum = self.datums[index.row()]
            datum.checked = value
            return True

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def init_data(self, data):
        self.beginResetModel()
        self._data = data
        self.endResetModel()
        self.layoutChanged.emit()  # Refresh whole view.
