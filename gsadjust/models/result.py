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
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt, QVariant

from .utils import format_numeric_column

# Constants for column headers
ADJSTA_STATION, ADJSTA_G, ADJSTA_SD = range(3)


class ResultsTableModel(QtCore.QAbstractTableModel):
    """
    Model to store network-adjusted gravity values.

    There is one ResultsTableModel per survey.
    """

    _headers = {ADJSTA_STATION: 'Station', ADJSTA_G: 'g', ADJSTA_SD: 'Std. dev.'}

    def __init__(self):
        super(ResultsTableModel, self).__init__()
        self.adjusted_stations = []

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.TextAlignmentRole:
            if orientation == Qt.Horizontal:
                return QVariant(int(Qt.AlignLeft | Qt.AlignVCenter))
            return QVariant(int(Qt.AlignRight | Qt.AlignVCenter))

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(
        self, adjusted_station, position, rows=1, index=QtCore.QModelIndex()
    ):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self.adjusted_stations.append(adjusted_station)
        self.endInsertRows()

    def rowCount(self, parent=None):
        return len(self.adjusted_stations)

    def columnCount(self, parent=None):
        return 3

    def data(self, index, role):
        if index.isValid():
            sta = self.adjusted_stations[index.row()]

            if role == Qt.DisplayRole:
                column = index.column()
                fn, *args = {
                    ADJSTA_STATION: (str, sta.station),
                    ADJSTA_G: (format, sta.g, '8.1f'),
                    ADJSTA_SD: (format, sta.sd, '1.1f'),
                }.get(column, (format_numeric_column, column))
                return fn(*args)

            if role == Qt.UserRole:
                return sta

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        If a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        if role == Qt.UserRole:
            self.adjusted_stations = []
            return QVariant()
        # return QtCore.QAbstractTableModel.setData(self, index, value, role)

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def copyToClipboard(self):
        """
        Copies results table to clipboard. Table needs to be selected first by
        clicking in the upper left corner.
        """
        clipboard = ""

        for r in range(self.rowCount()):
            for c in range(self.columnCount()):
                idx = self.index(r, c)
                clipboard += str(self.data(idx, role=Qt.DisplayRole))
                if c != (self.columnCount() - 1):
                    clipboard += '\t'
            clipboard += '\n'

        # copy to the system clipboard
        sys_clip = QtWidgets.QApplication.clipboard()
        sys_clip.setText(clipboard)

    def clearResults(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self.adjusted_stations = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()
