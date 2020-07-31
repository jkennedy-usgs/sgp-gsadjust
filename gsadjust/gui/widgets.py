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
    define progress bar
    """

    def __init__(self, parent=None, total=20, textmess='Progress'):
        super(ProgressBar, self).__init__(parent)
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)
        main_layout = QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 1)
        self.setLayout(main_layout)
        self.setWindowTitle(textmess)


class helpButton(QtWidgets.QAbstractButton):
    def __init__(self, pixmap, parent=None):
        super(helpButton, self).__init__(parent)
        self.pixmap = pixmap
        # self.clicked.connect(self.doNothing)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.drawPixmap(event.rect(), self.pixmap)

    def sizeHint(self):
        return self.pixmap.size()

    def doNothing(self):
        return
