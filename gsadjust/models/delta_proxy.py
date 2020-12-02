import datetime as dt
import logging

from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt, QVariant


class CustomProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self):
        super(CustomProxyModel, self).__init__()

    def lessThan(self, left, right):

        s1 = self.sourceModel().data(left, role=Qt.UserRole + 1)
        s2 = self.sourceModel().data(right, role=Qt.UserRole + 1)

        # Compare, accounting for None values.

        if s1 is None:
            return s2 is not None  # False if both None, True if only s1 is None.
        if s2 is None:
            return True  # False if only s2 is None.

        return s1 < s2
