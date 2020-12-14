"""
threads.py
=============

Used for loop animation.
Martin Fitzgerald, LearnPyQt
--------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""

from PyQt5.QtCore import (
    QObject,
    QRunnable,
    QThread,
    pyqtSignal,
    pyqtSlot,
)


class RunnerKilledException(Exception):
    pass


class GenericSignals(QObject):
    finished = pyqtSignal()
    error = pyqtSignal(object)
    data = pyqtSignal(object)
    result = pyqtSignal(object)


class RunnerBase(QRunnable):

    signals = GenericSignals()

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.args = args
        self.kwargs = args
        self.is_paused = False
        self.is_killed = False

    @pyqtSlot()
    def run(self):
        raise NotImplementedError

    def kill(self):
        self.is_killed = True

    def pause(self, pause=True):
        self.is_paused = pause


class ThreadBase(QThread):

    signals = GenericSignals()

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.args = args
        self.kwargs = args
        self.is_paused = False
        self.is_killed = False

    @pyqtSlot()
    def run(self):
        raise NotImplementedError

    def kill(self):
        self.is_killed = True

    def pause(self, pause=True):
        self.is_paused = pause
