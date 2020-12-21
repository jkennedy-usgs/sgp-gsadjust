"""
models/gravity.py
=================

PyQt model for gravity change table.
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
from PyQt5 import QtCore
from PyQt5.QtCore import Qt


class GravityChangeModel(QtCore.QAbstractTableModel):
    """
    Model to store gravity change between surveys.

    There is only one such model per campaign. Gravity change is calculated when the
    respective menu item is chosen.
    """

    def __init__(self, header, table, table_type="simple", parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._headers = {n: col for n, col in enumerate(header)}
        self.createArrayData(table, table_type)
        self.table_type = table_type

    def createArrayData(self, table, table_type):
        if table_type == "simple" or table_type == "list":
            array = np.array(table).transpose()
        elif table_type == "full":
            array = np.array(table)
        self.arraydata = array

    def rowCount(self, parent=None):
        return self.arraydata.shape[0]

    def columnCount(self, parent=None):
        return len(self._headers)

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            row = index.row()
            column = index.column()
            try:
                value = float(self.arraydata[row][column])
                if self.table_type != 'full':
                    return format(value, "0.1f")
                elif column != 1 and column != 2:
                    return format(value, "0.1f")
                else:
                    return format(value, "0.5f")
            except ValueError:
                return str(self.arraydata[row][column])

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)
