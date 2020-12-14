"""
gui/widgets.py
==============

Custom PyQt widgets
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

from PyQt5 import QtCore, QtGui, QtWidgets


class IncrMinuteTimeEdit(QtWidgets.QTimeEdit):
    """
    Provides a QTimeEdit with that increments at some minute interval.

    The PyQt Time edit doesn't allow for rolling-over the minutes to hours,
    apparently. I.e., the hours and minutes boxes are limited to 0-60. This
    implements that somewhat, but one can't decrement minutes below 00 in either
    the MinuteSection or HourSection. This manifests as when at a time, e.g.,
    3:00, stepBy isn't called when in the Minute Section because time is
    already at 0.

    Parameters
    ----------
    time : QTime
        Default time to show in dialog.

    """

    def __init__(self, time):
        super(IncrMinuteTimeEdit, self).__init__(time)
        self.step = 10

    def stepBy(self, steps):
        if self.currentSection() == self.HourSection:
            self.incrTime(steps)
        if self.currentSection() == self.MinuteSection:
            self.incrTime(steps)

    def incrTime(self, steps):
        hours = self.dateTime().time().hour()
        minutes = self.dateTime().time().minute() + steps * self.step
        if minutes < 0:
            self.setTime(QtCore.QTime(hours - 1, 60 + minutes))
        if minutes < 60:
            self.setTime(QtCore.QTime(hours, minutes))
        else:
            self.setTime(QtCore.QTime(hours + 1, 60 % minutes))


class ProgressBar(QtWidgets.QWidget):
    """
    Define progress bar
    """

    def __init__(self, parent=None, total=20, textmess="Progress"):
        super(ProgressBar, self).__init__(parent)
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)
        main_layout = QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 1)
        self.setLayout(main_layout)
        self.setWindowTitle(textmess)


class helpButton(QtWidgets.QAbstractButton):
    """
    Not used yet.
    """

    def __init__(self, pixmap, parent=None):
        super(helpButton, self).__init__(parent)
        self.pixmap = pixmap

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.drawPixmap(event.rect(), self.pixmap)

    def sizeHint(self):
        return self.pixmap.size()

    def doNothing(self):
        return
