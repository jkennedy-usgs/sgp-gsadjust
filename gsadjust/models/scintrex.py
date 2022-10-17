"""
models/scintrex.py
==================

PyQt model for Scintrex relative gravimeter data.
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
import numpy as np
from matplotlib.dates import num2date
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, QModelIndex

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0

ADJSTA_STATION, ADJSTA_G, ADJSTA_SD = range(3)
(
    SCINTREX_STATION,
    SCINTREX_DATE,
    SCINTREX_G,
    SCINTREX_SD,
    SCINTREX_X_TILT,
    SCINTREX_Y_TILT,
    SCINTREX_TEMP,
    SCINTREX_DUR,
    SCINTREX_REJ,
) = range(9)


class ScintrexTableModel(QtCore.QAbstractTableModel):
    """
    Model to store Scintrex data.

    There is one ScintrexTableModel (or BurrisTableModel) per station occupation.
    The model is created dynamically each time a station is selected in the tree
    view on the data tab, rather than stored in memory.

    The station position in the data hierarchy are stored, so that if a
    modification is triggered, the original data can be accessed and changed
    accordingly (keysurv,keyloop,keysta)

    by default, all table entries are checked (this can be modified to allow
    pre-check based on user criteria (tiltx, tilty,...)). Then, if one is
    unchecked, the keepdata property of the ChannelList object at the table
    row position is set to 0

    Attributes
    ----------
    _headers:             table header
    unchecked:            a dictionary of unchecked items. Keys are item
                          indexes, entries are states
    ChannelList_obj:      an object of ChannelList-type: used to store the
                          table data as the structured data
    arraydata:            an array representation of the data from the
                          ChannelList_obj
    """

    _headers = {
        SCINTREX_STATION: "Station",
        SCINTREX_DATE: "Date",
        SCINTREX_G: "g (\u00b5gal)",
        SCINTREX_SD: "sd (\u00b5gal)",
        SCINTREX_X_TILT: "X Tilt",
        SCINTREX_Y_TILT: "Y Tilt",
        SCINTREX_TEMP: "Temp (K)",
        SCINTREX_DUR: "dur (s)",
        SCINTREX_REJ: "rej",
    }

    signal_update_coordinates = QtCore.pyqtSignal()
    signal_adjust_update_required = QtCore.pyqtSignal()
    signal_uncheck_station = QtCore.pyqtSignal()
    signal_check_station = QtCore.pyqtSignal()

    def __init__(self, ChannelList_obj, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.unchecked = {}
        self.createArrayData(ChannelList_obj)

    def createArrayData(self, ChannelList_obj):
        """
        Create the np array data for table display, and update the
        ChannelList_obj.
        """
        self.ChannelList_obj = ChannelList_obj
        self.arraydata = np.concatenate(
            (
                ChannelList_obj.station,
                np.array(ChannelList_obj.t),
                np.array(ChannelList_obj.grav()),
                np.array(ChannelList_obj.sd),
                ChannelList_obj.tiltx,
                ChannelList_obj.tilty,
                ChannelList_obj.temp,
                ChannelList_obj.dur,
                ChannelList_obj.rej,
            )
        ).reshape(len(ChannelList_obj.t), 9, order="F")

    def rowCount(self, parent=None):
        return len(self.ChannelList_obj.t)

    def columnCount(self, parent):
        return len(self._headers)

    def flags(self, index):
        return Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            # view definition
            row = index.row()
            column = index.column()
            try:
                value = float(self.arraydata[row][column])
            except ValueError:
                value = self.arraydata[row][column]

            def format_datetime(dt):
                return num2date(float(dt)).strftime("%Y-%m-%d %H:%M:%S")

            fn, *args = {
                SCINTREX_DATE: (format_datetime, value),
                SCINTREX_REJ: (format, value, "2.0f"),
                SCINTREX_DUR: (format, value, "3.0f"),
                SCINTREX_G: (format, value, "8.1f"),
                SCINTREX_SD: (format, value, "8.1f"),
            }.get(column, (str, value))

            return fn(*args)

        if role == Qt.CheckStateRole:
            # check status definition
            if index.column() == 0:
                return self.checkState(index)

    def checkAll(self):
        self.ChannelList_obj.keepdata = [1] * len(self.ChannelList_obj.raw_grav)
        self.signal_adjust_update_required.emit()
        self.layoutChanged.emit()
        self.signal_check_station.emit()
        self.dataChanged.emit(QModelIndex(), QModelIndex())

    def uncheckAll(self):
        self.ChannelList_obj.keepdata = [0] * len(self.ChannelList_obj.raw_grav)
        self.signal_adjust_update_required.emit()
        self.layoutChanged.emit()
        self.signal_uncheck_station.emit()
        self.dataChanged.emit(QModelIndex(), QModelIndex())

    def checkState(self, index):
        """
        By default, everything is checked. If keepdata property from the
        ChannelList object is 0, it is unchecked
        """
        if self.ChannelList_obj.keepdata[index.row()] == 0:
            self.unchecked[index] = Qt.Unchecked
            return self.unchecked[index]
        else:
            return Qt.Checked

    def setData(self, index, value, role, silent=False):
        """
        if a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        if role == Qt.CheckStateRole and index.column() == 0:
            if value == Qt.Checked:
                self.ChannelList_obj.keepdata[index.row()] = 1
                self.signal_check_station.emit()
            elif value == Qt.Unchecked:
                self.unchecked[index] = value
                self.ChannelList_obj.keepdata[index.row()] = 0
                if not any(self.ChannelList_obj.keepdata):
                    self.signal_uncheck_station.emit()
            if not silent:
                self.signal_adjust_update_required.emit()
                self.dataChanged.emit(index, index)
            return True

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)
