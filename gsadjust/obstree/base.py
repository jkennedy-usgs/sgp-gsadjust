from PyQt5 import QtGui
from PyQt5.QtCore import Qt


class ObsTreeItemBase(QtGui.QStandardItem):
    """
    Base tree-view item used to populate data tree, used for Surveys, Loops, and Stations.
    Not used directly but inherited by ObsTreeStation, ...Loop, and ...Survey
    """

    def __init__(self):
        super(ObsTreeItemBase, self).__init__()
        self.setFlags(
            self.flags() | Qt.ItemIsEnabled | Qt.ItemIsEditable | Qt.ItemIsUserCheckable
        )

        self.setCheckState(Qt.Checked)
        self.fontweight = QtGui.QFont.Normal
        self.cellcolor = Qt.white
