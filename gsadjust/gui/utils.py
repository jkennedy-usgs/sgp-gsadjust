"""
gui/utils.py
============

Utility helper functions
------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""

from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from .messages import MessageBox


def copy_cells_to_clipboard(table):
    """
    Copies selected cells in table to clipboard.

    Parameters
    ----------
    table : PyQt table

    """

    if len(table.selectedIndexes()) > 0:
        # sort select indexes into rows and columns
        previous = table.selectedIndexes()[0]
        columns = []
        rows = []
        for index in table.selectedIndexes():  # columns first , then rows
            if previous.row() != index.row():
                columns.append(rows)
                rows = []
            rows.append(index.data())
            previous = index
        columns.append(rows)

        # add rows and columns to clipboard
        clipboard = ""
        ncols = len(columns[0])
        nrows = len(columns)

        # add header to clipboard
        for c in range(ncols):
            clipboard += table.model().headerData(c, Qt.Horizontal)
            clipboard += "\t"
        clipboard += "\n"

        for r in range(nrows):
            for c in range(ncols):
                clipboard += columns[r][c]
                if c != (ncols - 1):
                    clipboard += "\t"
            clipboard += "\n"

        # copy to the system clipboard
        sys_clip = QtWidgets.QApplication.clipboard()
        sys_clip.setText(clipboard)
    else:
        MessageBox.warning("Copy warning", "No rows selected (Ctrl-a to select all)")
