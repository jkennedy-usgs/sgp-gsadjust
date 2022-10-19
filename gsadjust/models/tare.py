"""
models/tare.py
==============

PyQt models for showing Tares on the Drift tab.
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

from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, QVariant
from matplotlib.dates import num2date, date2num

# Constants for column headers
TARE_DATETIME, TARE_TARE = range(2)


class TareTableModel(QAbstractTableModel):
    """
    Model to store tares (offsets)
    """

    _headers = {TARE_DATETIME: "Date", TARE_TARE: "Tare (\u00b5Gal)"}

    def __init__(self):
        super(TareTableModel, self).__init__(parent=None)
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
        return 2

    def data(self, index, role):
        if index.isValid():
            tare = self._data[index.row()]

            if role == Qt.DisplayRole or role == Qt.EditRole:
                if role == Qt.EditRole:
                    jeff = 1
                column = index.column()
                fn, *args = {
                    TARE_DATETIME: (
                        str,
                        dt.datetime.strftime(
                            num2date(tare.datetime), "%Y-%m-%d %H:%M:%S"
                        ),
                    ),
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
            self.dataChanged.emit(index, index)
            return True
        if role == Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                tare = self._data[index.row()]
                column = index.column()
                if len(str(value)) > 0:
                    if column == 0:
                        try:
                            tare.datetime = date2num(value)
                        except ValueError:
                            return False
                    elif column == 1:
                        tare.tare = float(value)
                    self.dataChanged.emit(index, index)
                return True
        if role == Qt.UserRole:
            self._data[index.row()] = value

    def deleteTare(self, idx):
        self.beginRemoveRows(idx, 0, 1)
        self.endRemoveRows()

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
