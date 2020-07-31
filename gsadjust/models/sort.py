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
from PyQt5.QtCore import QDate, QSortFilterProxyModel


class CustomSortingModel(QSortFilterProxyModel):
    """
    Used to sort by date in importAbsG dialog
    """

    def lessThan(self, left, right):
        col = left.column()
        dataleft = left.data()
        dataright = right.data()

        if col == 3:
            dataleft = QDate.fromString(dataleft, "MM/dd/yy").addYears(100)
            dataright = QDate.fromString(dataright, "MM/dd/yy").addYears(100)
        elif not col == 0:
            dataleft = float(dataleft)
            dataright = float(dataright)

        return dataleft < dataright
