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
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, QVariant

# Constants for column headers
TARE_DATE, TARE_TIME, TARE_TARE = range(3)


def format_numeric_column(column):
    """
    Format fn for simple numeric columns.
    """
    return column + 1


class tempStation:
    def __init__(self, station):
        self.__dict__ = station


class TareTableModel(QAbstractTableModel):
    """
    Model to store tares (offsets)
    """

    _headers = {TARE_DATE: 'Date', TARE_TIME: 'Time', TARE_TARE: 'Tare (\u00b5Gal)'}

    def __init__(self):
        super(TareTableModel, self).__init__()
        self._data = []

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.TextAlignmentRole:
            if orientation == Qt.Horizontal:
                return QVariant(int(Qt.AlignLeft | Qt.AlignVCenter))
            return QVariant(int(Qt.AlignRight | Qt.AlignVCenter))

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, tare, position=0, rows=1, index=QModelIndex()):
        self.beginInsertRows(QModelIndex(), position, position + rows - 1)
        self._data.append(tare)
        self.endInsertRows()

    def removeRow(self, index):
        tare = self.data(index, role=Qt.UserRole)
        self.beginRemoveRows(index, index.row(), 1)
        self._data.remove(tare)
        self.endRemoveRows()
        self.beginResetModel()
        self.endResetModel()

    def rowCount(self, parent=None):
        return len(self._data)

    def columnCount(self, parent=None):
        return 3

    def data(self, index, role):
        if index.isValid():
            tare = self._data[index.row()]

            if role == Qt.DisplayRole:
                column = index.column()
                fn, *args = {
                    TARE_DATE: (str, tare.date.toString('yyyy-MM-dd')),
                    TARE_TIME: (str, tare.time.toString()),
                    TARE_TARE: (format, tare.tare, "0.1f"),
                }.get(column)
                return fn(*args)

            elif role == Qt.CheckStateRole:
                if index.column() == 0:
                    return self.checkState(tare)

            elif role == Qt.UserRole:
                return tare

    def setData(self, index, value, role):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        if role == Qt.CheckStateRole and index.column() == 0:
            tare = self._data[index.row()]
            if value == Qt.Checked:
                tare.checked = 2
            elif value == Qt.Unchecked:
                tare.checked = 0
            return True
        if role == Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                tare = self._data[index.row()]
                column = index.column()
                if len(str(value)) > 0:
                    if column == 0:
                        tare.date = float(value)
                    elif column == 1:
                        tare.time = float(value)
                    elif column == 2:
                        tare.tare = float(value)
                    self.dataChanged.emit(index, index)
                return True
        if role == Qt.UserRole:
            self._data[index.row()] = value

    def deleteTare(self, idx):
        self.beginRemoveRows(idx, 0, 1)
        self.endRemoveRows

    def clearTares(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self._data = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def flags(self, index):
        return (
            Qt.ItemIsUserCheckable
            | Qt.ItemIsEnabled
            | Qt.ItemIsEditable
            | Qt.ItemIsSelectable
        )

    def checkState(self, tare):
        """
        By default, everything is checked. If keepdata property from the
        ChannelList object is 0, it is unchecked
        """
        if tare.checked == 0:
            return Qt.Unchecked
        else:
            return Qt.Checked

    def init_data(self, data):
        self.beginResetModel()
        self._data = data
        self.endResetModel()
        self.layoutChanged.emit()  # Refresh whole view.
