"""
obstree/base.py
===============

PyQt models for GSadjust data tree.
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
from PyQt5 import QtGui
from PyQt5.QtCore import Qt


class ObsTreeItemBase(QtGui.QStandardItem):
    """
    Base tree-view item used to populate data tree, used for Surveys, Loops,
    and Stations. Not used directly but inherited by ObsTreeStation, ...Loop,
    and ...Survey
    """

    def __init__(self):
        super(ObsTreeItemBase, self).__init__()
        self.setFlags(
            self.flags() | Qt.ItemIsEnabled | Qt.ItemIsEditable | Qt.ItemIsUserCheckable
        )

        self.setCheckState(Qt.Checked)
        self.fontweight = "normal" #QtGui.QFont.Normal
        self.cellcolor = Qt.white
        self.highlight = ''
    #
    def setHighlight(self, bool):
        self.highlight = bool