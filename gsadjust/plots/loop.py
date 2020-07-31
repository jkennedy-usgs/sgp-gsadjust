"""
gsa_plots.py
===============

GSadjust plotting module.
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
import threading
import time

import matplotlib
import numpy as np
from matplotlib.animation import TimedAnimation
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from matplotlib.dates import num2date
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from PyQt5 import QtCore, QtWidgets

############################################################################
# Loop animation
# - called from right-click context menu on tree view
# - animates stations in order they were observed
# - appears in a separate pop-up window
############################################################################

# You need to setup a signal slot mechanism, to
# send data to your GUI in a thread-safe way.
# Believe me, if you don't do this right, things
# go very very wrong..
class Communicate(QtCore.QObject):
    data_signal = QtCore.pyqtSignal(float, float, float)


def dataSendLoop(addData_callbackFunc, data):
    # Setup the signal-slot mechanism.
    mySrc = Communicate()
    mySrc.data_signal.connect(addData_callbackFunc)

    # Simulate some data
    n = data[0]
    y = data[1]
    dates = data[2]
    i = 0

    while True:
        if i > len(data[0]) - 1:
            time.sleep(2)
            mySrc.data_signal.emit(
                float(y[0]), n[0], dates[0]
            )  # <- Here you emit a signal!
            i = 0
        else:
            time.sleep(0.4)
            mySrc.data_signal.emit(
                float(y[i]), n[i], dates[i]
            )  # <- Here you emit a signal!
            i += 1


class PlotLoopAnimation(QtWidgets.QMainWindow):
    def __init__(self, data):
        super(PlotLoopAnimation, self).__init__()

        # Define the geometry of the main window
        self.setGeometry(200, 200, 800, 800)
        self.setWindowTitle("Loop animation")

        # Create FRAME_A
        self.FRAME_A = QtWidgets.QFrame(self)
        # self.FRAME_A.setStyleSheet("QWidget { background-color: %s }" % QtGui.QColor(210, 210, 235, 255).name())
        self.LAYOUT_A = QtWidgets.QGridLayout()
        self.FRAME_A.setLayout(self.LAYOUT_A)
        self.setCentralWidget(self.FRAME_A)

        lat_range = max(data[0]) - min(data[0])
        lon_range = max(data[1]) - min(data[1])
        ylim = [min(data[0]) - lat_range * 0.1, max(data[0]) + lat_range * 0.1]
        xlim = [min(data[1]) - lon_range * 0.1, max(data[1]) + lon_range * 0.1]
        # Place the matplotlib figure
        self.figure = CustomFigCanvas(xlim, ylim, len(data[0]))
        self.LAYOUT_A.addWidget(self.figure, *(0, 1))
        myDataLoop = threading.Thread(
            name='myDataLoop',
            target=dataSendLoop,
            daemon=True,
            args=(self.addData_callbackFunc, data),
        )
        myDataLoop.start()

    def addData_callbackFunc(self, x, value, date):
        # print("Add data: " + str(value))
        if x > -900:
            self.figure.addData(x, value, date)
        else:
            time.sleep(1)
            self.figure.__init__(self.figure.xlim, self.figure.ylim, len(self.figure.n))


class CustomFigCanvas(FigureCanvas, TimedAnimation):
    def __init__(self, xlim, ylim, n):
        n_head = 10
        self.addedY, self.addedX = [], []
        print('Matplotlib Version:', matplotlib.__version__)
        self.xlim, self.ylim = xlim, ylim
        # The data
        self.n = np.linspace(0, n - 1, n)
        # The window
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.ax1 = self.fig.add_subplot(111)

        # self.ax1 settings
        self.ax1.set_xlabel('Longitude')
        self.ax1.set_ylabel('Latitude')
        # self.points_blue = plt.Line2D([], [], marker='o', linewidth=0, color='0.5')
        self.points_blue = Line2D([], [], marker='o', linewidth=0, color='0.5')
        self.ax1.add_line(self.points_blue)
        self.lines_red = []
        a = list(np.logspace(1, 0, 5) / 10)
        a += [0, 0, 0, 0, 0]
        a.reverse()
        for i in range(n_head):
            self.lines_red.append(Line2D([], [], color='red', linewidth=4, alpha=a[i]))
        self.lines_gray = Line2D([], [], color='0.5', linewidth=1)
        self.points_red = Line2D(
            [], [], color='red', marker='o', markeredgecolor='r', linewidth=0
        )

        self.ax1.add_line(self.points_red)
        self.ax1.add_line(self.points_blue)
        self.ax1.add_line(self.lines_gray)
        for line in self.lines_red:
            self.ax1.add_line(line)
        self.ax1.set_xlim(xlim[0], xlim[1])
        self.ax1.set_ylim(ylim[0], ylim[1])
        self.title = self.ax1.text(
            0.15,
            0.95,
            "",
            bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
            transform=self.ax1.transAxes,
            ha="center",
        )
        ratio = 1.0
        xleft, xright = self.ax1.get_xlim()
        ybottom, ytop = self.ax1.get_ylim()
        self.ax1.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
        FigureCanvas.__init__(self, self.fig)
        TimedAnimation.__init__(self, self.fig, interval=400, blit=True)

    def new_frame_seq(self):
        return iter(range(self.n.size))

    def _init_draw(self):
        lines = [
            self.points_blue,
            self.points_red,
            self.lines_gray,
        ]  # , self.line1_tail]  #, self.line1_head]
        for l in lines:
            l.set_data([], [])
        for l in self.lines_red:
            l.set_data([], [])

    def addData(self, x, value, date):
        self.addedX.append(x)
        self.addedY.append(value)
        self.title.set_text(num2date(date).strftime('%Y-%m-%d %H:%M:%S'))

    def _step(self, *args):
        # Extends the _step() method for the TimedAnimation class.
        try:
            TimedAnimation._step(self, *args)
        except Exception as e:
            TimedAnimation._stop(self)
            pass

    def _draw_frame(self, framedata):
        MAP = 'winter'
        RESFACT = 10
        if len(self.addedX) > 0:
            self.points_blue.set_data(self.addedX, self.addedY)
            self.points_red.set_data(self.addedX[-1], self.addedY[-1])
        if len(self.addedX) > 1:
            self.lines_gray.set_data(self.addedX, self.addedY)
            # hrp = highResPoints(self.addedX[-2:], self.addedY[-2:])
            xHiRes, yHiRes = highResPoints(self.addedX[-2:], self.addedY[-2:], 1)
            npointsHiRes = len(xHiRes)
            # cm = plt.get_cmap(MAP)

            for idx, line in enumerate(self.lines_red):
                line.set_data(xHiRes[idx : idx + 2], yHiRes[idx : idx + 2])
            #
            # self.ax1.set_color_cycle([cm(1. * i / (npointsHiRes - 1))
            #                      for i in range(npointsHiRes - 1)])
            # for i in range(len(self.lines_red) - 1):
            #     self.ax1(xHiRes[i:i + 2], yHiRes[i:i + 2],
            #              alpha=float(i) / (npointsHiRes - 1),
            #              color=COLOR)
        # This is the drawing order
        self._drawn_artists = []
        self._drawn_artists.append(self.lines_gray)
        for line in self.lines_red:
            self._drawn_artists.append(line)
        for p in [self.points_blue, self.points_red, self.title]:
            self._drawn_artists.append(p)


def highResPoints(x, y, factor=10):
    '''
    Take points listed in two vectors and return them at a higher
    resolution. Create at least factor*len(x) new points that include the
    original points and those spaced in between.

    Returns new x and y arrays as a tuple (x,y).
    '''
    NPOINTS = 2
    RESFACT = 5
    # r is the distance spanned between pairs of points
    r = [0]
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        dy = y[i] - y[i - 1]
        r.append(np.sqrt(dx * dx + dy * dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())

    dr = rtot[-1] / (NPOINTS * RESFACT - 1)
    xmod = [x[0]]
    ymod = [y[0]]
    rPos = 0  # current point on walk along data
    rcount = 1
    while rPos < r.sum():
        x1, x2 = x[rcount - 1], x[rcount]
        y1, y2 = y[rcount - 1], y[rcount]
        dpos = rPos - rtot[rcount]
        theta = np.arctan2((x2 - x1), (y2 - y1))
        rx = np.sin(theta) * dpos + x1
        ry = np.cos(theta) * dpos + y1
        xmod.append(rx)
        ymod.append(ry)
        rPos += dr
        while rPos > rtot[rcount + 1]:
            rPos = rtot[rcount + 1]
            rcount += 1
            if rcount > rtot[-1]:
                break
    return xmod, ymod
