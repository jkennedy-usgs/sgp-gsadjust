import matplotlib
from matplotlib.dates import num2date

matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import networkx as nx
from PyQt5 import QtCore, QtWidgets, QtGui
from gui_objects import show_message, FigureDatumComparisonTimeSeries
from matplotlib.figure import Figure
from matplotlib.animation import TimedAnimation
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time
import threading
import matplotlib

def plot_compare_datum_to_adjusted(survey):
    """
    Bar plot of difference between specified datum (in datum table) and adjustment result.
    """
    results = survey.results_model
    diff, lbl = [], []
    for i in range(survey.datum_model.rowCount()):
        idx = survey.datum_model.index(i, 0)
        input_datum = survey.datum_model.data(idx, role=QtCore.Qt.UserRole)
        input_name = input_datum.station
        for ii in range(results.rowCount()):
            if results.data(results.index(ii, 0), role=QtCore.Qt.DisplayRole) == input_name:
                adj_g = results.data(results.index(ii, 1), role=QtCore.Qt.DisplayRole)
                diff.append(float(adj_g) - (input_datum.g - input_datum.meas_height * input_datum.gradient))
                lbl.append(input_name)
    fig, ax = plt.subplots()
    ind = range(len(diff))
    ax.bar(ind, diff)
    ax.set_title('Adj. datum - input datum (microGal)')
    ax.set_xticks(ind)
    ax.set_xticklabels(lbl)
    plt.show()


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

    while(True):
        if i > len(data[0]) - 1:
            time.sleep(2)
            # mySrc.data_signal.emit(-999., -999., -999.)
            mySrc.data_signal.emit(float(y[0]), n[0], dates[0])  # <- Here you emit a signal!
            i = 0
        else:
            time.sleep(0.4)
            mySrc.data_signal.emit(float(y[i]), n[i], dates[i])  # <- Here you emit a signal!
            i += 1


class CustomMainWindow(QtWidgets.QMainWindow):
    def __init__(self, data):
        super(CustomMainWindow, self).__init__()

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
        ylim=[min(data[0]) - lat_range * 0.1, max(data[0]) + lat_range * 0.1]
        xlim=[min(data[1]) - lon_range * 0.1, max(data[1]) + lon_range * 0.1]
        # Place the matplotlib figure
        self.myFig = CustomFigCanvas(xlim, ylim, len(data[0]))
        self.LAYOUT_A.addWidget(self.myFig, *(0, 1))
        myDataLoop = threading.Thread(name='myDataLoop', target=dataSendLoop, daemon=True,
                                      args=(self.addData_callbackFunc, data))
        myDataLoop.start()
        self.show()

    def addData_callbackFunc(self, x, value, date):
        # print("Add data: " + str(value))
        if x > -900:
            self.myFig.addData(x, value, date)
        else:
            time.sleep(1)
            self.myFig.__init__(self.myFig.xlim, self.myFig.ylim, len(self.myFig.n))

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
        self.points_blue = plt.Line2D([], [], marker='o', linewidth=0, color='0.5')
        self.lines_red = []
        a = list(np.logspace(1,0,5)/10)
        a += [0, 0, 0, 0, 0]
        a.reverse()
        for i in range(n_head):
            self.lines_red.append(Line2D([], [], color='red', linewidth=4, alpha = a[i]))
        self.lines_gray = Line2D([], [], color='0.5', linewidth=1)
        self.points_red = Line2D([], [], color='red', marker='o', markeredgecolor='r', linewidth=0)

        self.ax1.add_line(self.points_red)
        self.ax1.add_line(self.points_blue)
        self.ax1.add_line(self.lines_gray)
        for line in self.lines_red:
            self.ax1.add_line(line)
        self.ax1.set_xlim(xlim[0], xlim[1])
        self.ax1.set_ylim(ylim[0], ylim[1])
        self.title = self.ax1.text(0.15, 0.95, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax1.transAxes, ha="center")
        ratio = 1.0
        xleft, xright = self.ax1.get_xlim()
        ybottom, ytop = self.ax1.get_ylim()
        self.ax1.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
        FigureCanvas.__init__(self, self.fig)
        TimedAnimation.__init__(self, self.fig, interval=400, blit=True)

    def new_frame_seq(self):
        return iter(range(self.n.size))

    def _init_draw(self):
        lines = [self.points_blue, self.points_red, self.lines_gray]  #, self.line1_tail]  #, self.line1_head]
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
        # try:
        TimedAnimation._step(self, *args)
        # except Exception as e:
        #     self.abc += 1
        #     print(str(self.abc))
        #     TimedAnimation._stop(self)
        #     pass

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
            cm = plt.get_cmap(MAP)
            for idx, line in enumerate(self.lines_red):
                line.set_data(xHiRes[idx:idx+2], yHiRes[idx:idx+2])
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


def highResPoints(x,y,factor=10):
    '''
    Take points listed in two vectors and return them at a higher
    resultion. Create at least factor*len(x) new points that include the
    original points and those spaced in between.

    Returns new x and y arrays as a tuple (x,y).
    '''
    NPOINTS = 2
    RESFACT = 5
    # r is the distance spanned between pairs of points
    r = [0]
    for i in range(1,len(x)):
        dx = x[i]-x[i-1]
        dy = y[i]-y[i-1]
        r.append(np.sqrt(dx*dx+dy*dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())

    dr = rtot[-1]/(NPOINTS*RESFACT-1)
    xmod=[x[0]]
    ymod=[y[0]]
    rPos = 0 # current point on walk along data
    rcount = 1
    while rPos < r.sum():
        x1,x2 = x[rcount-1],x[rcount]
        y1,y2 = y[rcount-1],y[rcount]
        dpos = rPos-rtot[rcount]
        theta = np.arctan2((x2-x1),(y2-y1))
        rx = np.sin(theta)*dpos+x1
        ry = np.cos(theta)*dpos+y1
        xmod.append(rx)
        ymod.append(ry)
        rPos+=dr
        while rPos > rtot[rcount+1]:
            rPos = rtot[rcount+1]
            rcount+=1
            if rcount>rtot[-1]:
                break

    return xmod,ymod

def plot_loop_animation(obstreeloop):
    lat, lon, dates = [], [], []
    for i in range(obstreeloop.rowCount()):
        station = obstreeloop.child(i)
        lat.append(station.lat[0])
        lon.append(station.long[0])
        dates.append(station.tmean)

    myGUI = CustomMainWindow([lat, lon, dates])

def plot_adjust_residual_histogram(survey):
    """
    Matplotlib histogram of delta-g residuals
    """
    # the histogram of the data
    try:
        nrows = survey.delta_model.rowCount()
        rlist = []
        for i in range(nrows):
            idx = survey.delta_model.index(i, 8)
            results = survey.delta_model.data(idx, role=QtCore.Qt.DisplayRole)
            d = float(results.value())
            if d > -998:
                rlist.append(float(results.value()))

        plt.hist(rlist, 20, facecolor='green', alpha=0.75)
        plt.xlabel('Residual (microGal)')
        plt.ylabel('Count')
        plt.grid(True)
        plt.show()

        return plt.gcf()
    except Exception as e:
        return e


def plot_LOO_analysis(self, x_all, adj_g_all, obs_g_all, datums):
    """
    Plot results from leave-one-out analysis.
    """
    fig = plt.figure(figsize=(13, 8))
    fig.hold
    i = 0
    for x, adj_g, obs_g in zip(x_all, adj_g_all, obs_g_all):
        y1 = [y - obs_g[0] for y in obs_g]
        y2 = [y - adj_g[0] for y in adj_g]

        a = plt.plot(x, y1, '-o', label=datums[i] + '_obs')
        line_color = a[0].get_color()
        plt.plot(x, y2, '--o', c=line_color, label=datums[i] + '_adj')
        i += 1
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    fig.show()


def plot_network_graph(survey, coords, shape='circular'):
    """
    Networkx plot of network. If shape == 'map', accurate coordinates must be present in the input file.
    :param shape: 'circular' or 'map'
    """
    edges = nx.MultiGraph()
    disabled_edges = nx.MultiGraph()
    datum_nodelist, nondatum_nodelist = [], []

    delta_model = survey.delta_model
    if delta_model.rowCount() == 0:
        msg = show_message('Delta table is empty. Unable to plot network graph', 'Plot error')
    else:
        for i in range(delta_model.rowCount()):
            ind = delta_model.index(i, 0)
            delta = delta_model.data(ind, role=QtCore.Qt.UserRole)
            chk = delta_model.data(ind, QtCore.Qt.CheckStateRole)
            if chk == 2:
                edges.add_edge(delta.sta1, delta.sta2, key=ind)
            else:
                disabled_edges.add_edge(delta.sta1, delta.sta2, key=ind)
            datum_names = survey.datum_model.datum_names()
            for station_name in [delta.sta1, delta.sta2]:
                if station_name in datum_names:
                    if station_name not in datum_nodelist:
                        datum_nodelist.append(station_name)
                        continue
                elif station_name not in nondatum_nodelist:
                    nondatum_nodelist.append(station_name)
                    edges.add_node(station_name)
                    disabled_edges.add_node(station_name)

        H = nx.Graph(edges)

        if shape == 'circular':
            pos = nx.circular_layout(H)
        elif shape == 'map':
            new_dict = {}
            for k, v in coords.items():
                new_dict[k] = (v[0], v[1])
            pos = new_dict
            pos_label = {}
            y_off = 0.001  # offset on the y axis

            for k, v in pos.items():
                pos_label[k] = (v[0], v[1] + y_off)

        if not nx.is_connected(H):
            gs = [H.subgraph(c) for c in nx.connected_components(H)]
            for idx, g in enumerate(gs):
                plt.subplot(1, len(gs), idx + 1)
                nx.draw_networkx_edges(g, pos, width=1, alpha=0.4, node_size=0, edge_color='k')
                nx.draw_networkx_nodes(g, pos, node_color='w', alpha=0.4, with_labels=True)
                nx.draw_networkx_labels(g, pos, font_color='orange')
                plt.title('Networks are disconnected!')
        else:
            # edge width is proportional to number of delta-g's
            edgewidth = []
            for (u, v, d) in H.edges(data=True):
                edgewidth.append(len(edges.get_edge_data(u, v)) * 2 - 1)

            nx.draw_networkx_edges(H, pos, width=edgewidth, alpha=0.4, node_size=0, edge_color='k')
            nx.draw_networkx_edges(disabled_edges, pos, width=1, alpha=0.4, node_size=0, edge_color='r')
            nx.draw_networkx_nodes(H, pos, node_size=120, nodelist=datum_nodelist, node_color="k", node_shape='^',
                                   with_labels=True, alpha=0.8)
            nx.draw_networkx_nodes(H, pos, node_size=120, nodelist=nondatum_nodelist, node_color="k",
                                   node_shape='o',
                                   with_labels=True, alpha=0.3)
            nx.draw_networkx_labels(H, pos, font_color='r')

        plt.autoscale(enable=True, axis='both', tight=True)
        plt.axis('off')
        plt.show()
#
#
# def plot_datum_comparison_timeseries(self):
#     """
#     Creating a PyQt window heree in preparation for eventually having a tabbed plot window with various diagnostics
#     """
#     _ = FigureDatumComparisonTimeSeries(self.obsTreeModel)
