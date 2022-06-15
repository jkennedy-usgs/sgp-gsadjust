import datetime as dt
import logging

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt, QEvent


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


class CustomDatumProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self):
        super(CustomDatumProxyModel, self).__init__()

    def lessThan(self, left, right):
        if self.sortColumn() == 0:
            s1 = self.sourceModel()._data[left.row()].checked
            s2 = self.sourceModel()._data[right.row()].checked

        else:
            s1 = self.sourceModel().data(left)
            s2 = self.sourceModel().data(right)

        # Compare, accounting for None values.

        if s1 is None:
            return s2 is not None  # False if both None, True if only s1 is None.
        if s2 is None:
            return True  # False if only s2 is None.

        return s1 < s2


class CheckBoxDelegate(QtWidgets.QItemDelegate):
    """
    A delegate that places a fully functioning QCheckBox cell of the column to which it's applied.
    """
    def __init__(self, parent):
        QtWidgets.QItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        """
        Important, otherwise an editor is created if the user clicks in this cell.
        """
        return None

    def paint(self, painter, option, index):
        """
        Paint a checkbox without the label.
        """
        checkstate = index.model().data(index)
        self.drawCheck(painter, option, option.rect, checkstate)

    def editorEvent(self, event, model, option, index):
        '''
        Change the data in the model and the state of the checkbox
        if the user presses the left mousebutton and this cell is editable. Otherwise do nothing.
        '''
        if not int(index.flags() & Qt.ItemIsEditable) > 0:
            return False

        if event.type() == QEvent.MouseButtonRelease and event.button() == Qt.LeftButton:
            # Change the checkbox-state
            self.setModelData(None, model, index)
            return True

        return False


    def setModelData (self, editor, model, index):
        '''
        The user wanted to change the old state in the opposite.
        '''
        model.setData(index, 1 if int(index.data()) == 0 else 0, Qt.EditRole)