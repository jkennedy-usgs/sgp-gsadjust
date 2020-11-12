"""
A simple logging window to show messages logged using Python logging.
"""

from PyQt5.QtWidgets import QTextEdit
import logging
from PyQt5.QtGui import QColor, QFont
from PyQt5.QtCore import QSize

LEVEL_COLORS = {
    logging.DEBUG: QColor('grey'),
    logging.INFO: QColor('black'),
    logging.WARNING: QColor('yellow'),
    logging.ERROR: QColor('red'),
}

class LoggerWidget(QTextEdit):

    def __init__(self, parent=None):
        super().__init__(parent)

        font = QFont('Courier New')
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
        color = LEVEL_COLORS.get(record.levelno, QColor('black'))

        self._widget.setTextColor(color)
        self._widget.append(msg)

    def write(self, m):
        pass

