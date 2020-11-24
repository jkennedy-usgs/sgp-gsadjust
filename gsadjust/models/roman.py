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
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, QVariant

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0

ROMAN_FROM, ROMAN_TO, ROMAN_DELTA = range(3)
# Constants for column headers
(
    DATUM_STATION,
    DATUM_DATE,
    DATUM_G,
    DATUM_SD,
    N_SETS,
    MEAS_HEIGHT,
    GRADIENT,
    DATUM_RESIDUAL,
) = range(8)


class RomanTableModel(QtCore.QAbstractTableModel):
    """
    Model to store the individual delta-g's calculated using the Roman-method.
    Shown on the drift tab.
    """

    _headers = {ROMAN_FROM: 'From', ROMAN_TO: 'To', ROMAN_DELTA: 'Delta g'}

    def __init__(self):
        super(RomanTableModel, self).__init__()
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
        return 3

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            delta = self._data[index.row()]

            if role == Qt.DisplayRole:
                column = index.column()
                fn, *args = {
                    ROMAN_FROM: (str, delta.sta1),
                    ROMAN_TO: (str, delta.sta2),
                    ROMAN_DELTA: (format, delta.dg, "0.1f"),
                }.get(column)
                return fn(*args)

            elif role == Qt.UserRole:
                # check status definition
                return delta

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        If a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        if role == Qt.CheckStateRole and index.column() == 0:
            datum = self._data[index.row()]
            datum.checked = value
            return True

        if role == Qt.EditRole:
            if index.isValid() and index.row() >= 0 and value:
                datum = self._data[index.row()]
                column = index.column()

                if column == DATUM_STATION:
                    datum.station = value
                elif column == DATUM_G:
                    datum.g = float(value)
                elif column == DATUM_SD:
                    datum.sd = float(value)
                elif column == DATUM_DATE:
                    datum.date = value
                self.dataChanged.emit(index, index)

            return True

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def init_data(self, data):
        self._data = data
        self.layoutChanged.emit()  # Refresh whole view.
