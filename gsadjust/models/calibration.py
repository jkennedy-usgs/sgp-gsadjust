"""
models/calibration.py
===============

PyQt model for showing meter calibration coefficients in AdjustOptions dialog.
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
from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtCore import Qt


class MeterCalibrationModel(QStandardItemModel):
    def __init__(self):
        super(MeterCalibrationModel, self).__init__()
        self.setColumnCount(2)
        self.setHorizontalHeaderLabels(["Meter", "Calibration factor"])

    def data(self, index, role):
        if index.isValid():
            if role == Qt.DisplayRole:
                column = index.column()
                # if role == Qt.DisplayRole or role == Qt.EditRole:
                if column == 0:
                    return self.itemFromIndex(index).text()
                if column == 1:
                    return format(float(self.itemFromIndex(index).text()), "0.6f")

