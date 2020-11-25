import sys
import time

from PyQt5.QtCore import (
    QObject,
    QRunnable,
    Qt,
    QThread,
    QThreadPool,
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


# threadpool = QThreadPool()
