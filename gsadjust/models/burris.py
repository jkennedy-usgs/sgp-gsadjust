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
import numpy as np
from matplotlib.dates import num2date
from PyQt5.QtCore import QAbstractTableModel, Qt, pyqtSignal, QModelIndex

# Constants for column headers
(
    BURRIS_STATION,
    BURRIS_OPER,
    BURRIS_METER,
    BURRIS_DATE,
    BURRIS_G,
    BURRIS_DIAL,
    BURRIS_FEEDBACK,
    BURRIS_TIDE,
    BURRIS_ELEVATION,
    BURRIS_LAT,
    BURRIS_LONG,
) = range(11)


# noinspection PyUnresolvedReferences
class BurrisTableModel(QAbstractTableModel):
    """
    Model to store Burris data.

    There is one BurrisTableModel (or ScintrexTableModel) per station occupation.
    The model is created dynamically each time a station is selected in the tree
    view on the data tab, rather than stored in memory.

    The station position in the data hierarchy are stored, so that if a
    modification is triggered, the original data can be accessed and changed
    accordingly (keysurv,keyloop,keysta)

    by default, all table entries are checked (this can be modified to allow
    pre-check based on user criteria (tiltx, tilty,...)). Then, if one is
    unchecked, the keepdata property of the ChannelList object at the table
    row position is set to 0

    properties:
    _headers:             table header
    unchecked:            a dictionary of unchecked items. Keys are item
                          indexes, entries are states
    ChannelList_obj:      an object of ChannelList-type: used to store the
                          table data as the structured data
    arraydata:            an array representation of the data from the
                          ChannelList_obj
    keysurv:              the key of a survey object within the grav_obj
                          hierarchy
    keyloop:              the key of a loop object within the grav_obj
                          hierarchy
    keysta:               the key of a station object within the grav_obj
                          hierarchy
    """

    _headers = {
        BURRIS_STATION: "Station",
        BURRIS_OPER: "Oper",
        BURRIS_METER: "Meter",
        BURRIS_DATE: "Date",
        BURRIS_G: "g (\u00b5gal)",
        BURRIS_DIAL: "Dial setting",
        BURRIS_FEEDBACK: "Feedback",
        BURRIS_TIDE: "Tide corr.",
        BURRIS_ELEVATION: "Elev.",
        BURRIS_LAT: "Lat",
        BURRIS_LONG: "Long",
    }

    signal_update_coordinates = pyqtSignal()
    signal_adjust_update_required = pyqtSignal()
    signal_uncheck_station = pyqtSignal()
    signal_check_station = pyqtSignal()

    def __init__(self, ChannelList_obj, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.unchecked = {}
        self.ChannelList_obj = None
        self.createArrayData(ChannelList_obj)
        self.arraydata = None

    def createArrayData(self, ChannelList_obj):
        """
        Create the np array data for table display, and update the
        ChannelList_obj. This function can be called from outside to update the
        table display
        """
        # channel list attributes aren't specified in advance, so don't worry whether it's Burris or CG5 data
        self.ChannelList_obj = ChannelList_obj
        try:
            self.arraydata = np.concatenate(
                (
                    ChannelList_obj.station,
                    ChannelList_obj.oper,
                    ChannelList_obj.meter,
                    ChannelList_obj.t,
                    np.array(ChannelList_obj.grav()),
                    np.array(ChannelList_obj.dial),
                    np.array(ChannelList_obj.feedback),
                    np.array(ChannelList_obj.etc),
                    np.array(ChannelList_obj.elev),
                    np.array(ChannelList_obj.lat),
                    np.array(ChannelList_obj.long),
                )
            ).reshape(len(ChannelList_obj.t), 11, order='F')
        except Exception:
            return

    def rowCount(self, parent=None):
        return len(self.ChannelList_obj.t)

    def columnCount(self, parent=None):
        return len(self._headers)

    def flags(self, index):
        return (
            Qt.ItemIsUserCheckable
            | Qt.ItemIsEnabled
            | Qt.ItemIsSelectable
            | Qt.ItemIsEditable
        )

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            row = index.row()
            column = index.column()
            try:
                value = float(self.arraydata[row][column])
            except ValueError:
                value = self.arraydata[row][column]

            def format_datetime(dt):
                if dt > 50000:
                    return num2date(dt - 719163).strftime('%Y-%m-%d %H:%M:%S')
                else:
                    return num2date(dt).strftime('%Y-%m-%d %H:%M:%S')

            fn, *args = {
                BURRIS_DATE: (format_datetime, value),
                BURRIS_TIDE: (format, value, "0.1f"),
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
        If a row is unchecked, update the keepdata value to 0 setData launched
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

        elif role == Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                if len(str(value)) > 0:
                    attr = None
                    if index.column() == 1:  # Oper
                        self.arraydata[index.row()][index.column()] = value
                        attr = 'oper'
                    if index.column() == 2:  # Meter
                        self.arraydata[index.row()][index.column()] = value
                        attr = 'meter'
                    if index.column() == 9:  # Lat
                        self.arraydata[index.row()][index.column()] = float(value)
                        attr = 'lat'
                        self.signal_update_coordinates.emit()
                    if index.column() == 10:  # Long
                        self.arraydata[index.row()][index.column()] = float(value)
                        attr = 'long'
                        self.signal_update_coordinates.emit()
                    if attr is not None:
                        self.dataChanged.emit(index, index)
                        return True

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)
