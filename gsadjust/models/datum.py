"""
models/datum.py
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

from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, QVariant, pyqtSignal

from .utils import format_numeric_column

# Constants for column headers
(
    DATUM_STATION,
    DATUM_DATE,
    DATUM_G,
    DATUM_SD,
    N_SETS,
    MEAS_HEIGHT,
    GRADIENT,
    DATUM_RESIDUAL,
) = range(8)


class tempStation:
    def __init__(self, station):
        self.__dict__ = station


# noinspection PyUnresolvedReferences
class DatumTableModel(QAbstractTableModel):
    """
    Model to store Datums, shown on the adjust tab.
    """

    _headers = {  # As map, so do not need to be kept in order with the above.
        DATUM_STATION: "Station",
        DATUM_G: "g",
        DATUM_SD: "Std. dev.",
        DATUM_DATE: "Date",
        MEAS_HEIGHT: "Meas. height",
        GRADIENT: "Gradient",
        DATUM_RESIDUAL: "Residual",
        N_SETS: "# sets",
    }

    _attrs = {  # From column constants to object attributes, for setting.
        DATUM_STATION: ("station", str),
        DATUM_G: ("g", float),
        DATUM_SD: ("sd", float),
        DATUM_DATE: ("date", lambda x: x),  # pass through
        MEAS_HEIGHT: ("meas_height", float),
        GRADIENT: ("gradient", float),
    }

    signal_adjust_update_required = pyqtSignal()
    signal_datum_table_updated = pyqtSignal()

    def __init__(self):
        super(DatumTableModel, self).__init__()
        self._data = []

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.TextAlignmentRole:
            if orientation == Qt.Horizontal:
                return QVariant(int(Qt.AlignLeft | Qt.AlignVCenter))
            return QVariant(int(Qt.AlignRight | Qt.AlignVCenter))

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, datum, position, rows=1, index=QModelIndex()):
        self.beginInsertRows(QModelIndex(), position, position + rows - 1)
        self._data.append(datum)
        self.endInsertRows()

    def removeRow(self, index):
        datum = self.data(index, role=Qt.UserRole)
        self.beginRemoveRows(index, index.row(), 1)
        self._data.remove(datum)
        self.endRemoveRows()
        self.beginResetModel()
        self.endResetModel()

    def rowCount(self, parent=None):
        return len(self._data)

    def columnCount(self, parent=None):
        return len(self._headers)

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            datum = self._data[index.row()]
            column = index.column()
            if role == Qt.DisplayRole or role == Qt.EditRole:
                # To accommodate old save files
                def get_nsets():
                    try:
                        return datum.n_sets
                    except:
                        return "NA"

                fn, *args = {
                    DATUM_SD: (format, datum.sd, "0.2f"),
                    DATUM_G: (format, datum.g, "8.1f"),
                    DATUM_STATION: (str, datum.station),
                    DATUM_DATE: (str, datum.date),
                    MEAS_HEIGHT: (format, datum.meas_height, "0.2f"),
                    GRADIENT: (format, datum.gradient, "0.2f"),
                    DATUM_RESIDUAL: (format, datum.residual, "0.1f"),
                    N_SETS: (get_nsets,),
                }.get(column, (format_numeric_column, column))

                return fn(*args)

            elif role == Qt.CheckStateRole:
                # check status definition
                if index.column() == 0:
                    return self.checkState(datum)

            elif role == Qt.UserRole:
                # check status definition
                return datum

    def checkState(self, datum):
        """
        By default, everything is checked. If keepdata property from the
        ChannelList object is 0, it is unchecked
        """
        if datum.checked == 0:
            return Qt.Unchecked
        else:
            return Qt.Checked

    def setData(self, index, value, role):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched
        when role is acting value is Qt.Checked or Qt.Unchecked
        """
        if role == Qt.CheckStateRole and index.column() == 0:
            datum = self._data[index.row()]
            if value == Qt.Checked:
                datum.checked = 2
            elif value == Qt.Unchecked:
                datum.checked = 0
            self.dataChanged.emit(index, index, [])
            self.signal_adjust_update_required.emit()
            return True

        if role == Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                if value:
                    try:
                        datum = self._data[index.row()]
                        column = index.column()
                        # Ideally they other columns wouldn't be editable at all, but
                        # the user can select them and enter new values. Here we
                        # discard them unless they're in an editable column.
                        #
                        # Should me able to make non-editable columns readonly using a
                        # proxy model, e.g.
                        # https://stackoverflow.com/questions/22886912
                        if column in [
                            DATUM_DATE,
                            DATUM_G,
                            DATUM_SD,
                            MEAS_HEIGHT,
                            GRADIENT,
                        ]:
                            attr, vartype = self._attrs.get(column, (None, None))
                            if attr:
                                setattr(datum, attr, vartype(value))
                            self.dataChanged.emit(index, index, [Qt.EditRole])
                    except ValueError:
                        pass
            return True

        if role == Qt.UserRole:
            self._data[index.row()] = value
            self.dataChanged.emit(index, index, [])

    def flags(self, index):
        return (
            Qt.ItemIsUserCheckable
            | Qt.ItemIsEnabled
            | Qt.ItemIsSelectable
            | Qt.ItemIsEditable
        )

    def clearDatums(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self._data = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def datum_names(self):
        dn = []
        for datum in self._data:
            dn.append(datum.station)
        return dn

    def init_data(self, data):
        self.beginResetModel()
        self._data = data
        self.endResetModel()
        self.layoutChanged.emit()  # Refresh whole view.)
