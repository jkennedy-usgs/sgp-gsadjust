"""
A simple logging window to show messages logged using Python logging.
Martin Fitzgerald, LearnPyQt

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""

import logging

from PyQt5.QtCore import QSize
from PyQt5.QtGui import QColor, QFont
from PyQt5.QtWidgets import QTextEdit

LEVEL_COLORS = {
    logging.DEBUG: QColor("grey"),
    logging.INFO: QColor("black"),
    logging.WARNING: QColor("green"),
    logging.ERROR: QColor("red"),
}


class LoggerWidget(QTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)

        font = QFont("Courier New")
        self.setFont(font)
        self.setMinimumSize(QSize(600, 400))
        self.setReadOnly(True)
        self.handler = CustomLogger(widget=self)


class CustomLogger(logging.Handler):
    def __init__(self, *args, widget=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._widget = widget

    def emit(self, record):
        msg = self.format(record)
        color = LEVEL_COLORS.get(record.levelno, QColor("black"))

        self._widget.setTextColor(color)
        self._widget.append(msg)

    def write(self, m):
        pass
