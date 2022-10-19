"""
plots/network.py
================

Plots network graphs using the Networkx library.
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
import matplotlib
import networkx as nx
import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as Toolbar

from ..gui.messages import MessageBox


class PlotNetworkGraph(QtWidgets.QDialog):
    """
    Networkx plot of network. If shape == 'map', accurate coordinates must be
    present in the input file.

    Parameters
    ----------
    survey
    coords
    shape : {'Circular', 'map'}
    """

    def __init__(self, survey, coords, shape="circular", parent=None):
        super(PlotNetworkGraph, self).__init__(parent)
        self.setWindowTitle("Network graph, Survey " + survey.name)
        self.survey = survey
        self.coords = coords
        self.shape = shape
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = Toolbar(self.canvas, self)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)
        self.plot_network()

    def plot_network(self):
        try:
            edges, disabled_edge, datum_nodelist, nondatum_nodelist = self.get_data()
            self.plot(edges, disabled_edge, datum_nodelist, nondatum_nodelist)
        except KeyError as e:
            MessageBox.warning(
                "Plot error", f"Error plotting network graph: {e.args[0]}",
            )

    def get_data(self):
        edges = nx.MultiGraph()
        disabled_edges = nx.MultiGraph()
        datum_nodelist, nondatum_nodelist = [], []

        deltas = self.survey.deltas
        if len(deltas) == 0:
            MessageBox.warning(
                "Plot error",
                "Delta table is empty. Unable to plot network graph",
            )
        else:
            for i, delta in enumerate(deltas):
                key = f"delta_{i}"
                if delta.checked:
                    edges.add_edge(delta.sta1, delta.sta2, key=key)
                else:
                    disabled_edges.add_edge(delta.sta1, delta.sta2, key=key)

                datum_names = [datum.station for datum in self.survey.datums]
                for station_name in [delta.sta1, delta.sta2]:
                    if station_name in datum_names:
                        if station_name not in datum_nodelist:
                            datum_nodelist.append(station_name)
                            continue
                    elif station_name not in nondatum_nodelist:
                        nondatum_nodelist.append(station_name)
                        edges.add_node(station_name)
                        disabled_edges.add_node(station_name)
        return edges, disabled_edges, datum_nodelist, nondatum_nodelist

    def format_map_axis(self, ax, shape):
        if shape == "circular":
            ax.set_xlim(-1.2, 1.2)
            ax.set_ylim(-1.2, 1.2)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        elif shape == "map":
            border = 0.1
            self.figure.tight_layout()
            xrange = np.abs(self.xmax - self.xmin)
            yrange = np.abs(self.ymax - self.ymin)
            ax.set_xlim(self.xmin - xrange * border, self.xmax + xrange * border)
            ax.set_ylim(self.ymin - yrange * border, self.ymax + yrange * border)
            ax.ticklabel_format(useOffset=False)
            ax.set_xlabel("(Coordinates are not projected)")

    def plot(self, edges, disabled_edges, datum_nodelist, nondatum_nodelist):
        self.figure.clear()
        H = nx.Graph(edges)

        if self.shape == "circular":
            pos = nx.circular_layout(H)
        elif self.shape == "map":
            pos = {}
            for k, v in self.coords.items():
                pos[k] = (v[0], v[1])

            self.xmin = min([x[0] for x in pos.values()])
            self.ymin = min([x[1] for x in pos.values()])
            self.xmax = max([x[0] for x in pos.values()])
            self.ymax = max([x[1] for x in pos.values()])

        if not nx.is_connected(H):
            gs = [H.subgraph(c) for c in nx.connected_components(H)]
            for idx, g in enumerate(gs):
                ax = self.figure.add_subplot(1, len(gs), idx + 1)
                nx.draw_networkx_edges(
                    g, pos, ax=ax, width=1, alpha=0.4, node_size=0, edge_color="k"
                )
                nx.draw_networkx_nodes(
                    g, pos, ax=ax, node_color="w", alpha=0.4, with_labels=True
                )
                nx.draw_networkx_labels(g, pos, ax=ax, font_color="orange")
                ax.set_title("Networks are disconnected!")
                self.format_map_axis(ax, self.shape)
        else:
            # edge width is proportional to number of delta-g's
            edgewidth = []
            ax = self.figure.add_subplot(111)
            for (u, v, d) in H.edges(data=True):
                edgewidth.append(len(edges.get_edge_data(u, v)) * 2 - 1)

            nx.draw_networkx_edges(
                H, pos, ax=ax, width=edgewidth, alpha=0.4, node_size=0, edge_color="k"
            )
            nx.draw_networkx_edges(
                disabled_edges,
                pos,
                ax=ax,
                width=1,
                alpha=0.4,
                node_size=0,
                edge_color="r",
            )
            nx.draw_networkx_nodes(
                H,
                pos,
                ax=ax,
                node_size=120,
                nodelist=datum_nodelist,
                node_color="k",
                node_shape="^",
                with_labels=True,
                alpha=0.8,
            )
            nx.draw_networkx_nodes(
                H,
                pos,
                ax=ax,
                node_size=120,
                nodelist=nondatum_nodelist,
                node_color="k",
                node_shape="o",
                with_labels=True,
                alpha=0.3,
            )
            nx.draw_networkx_labels(H, pos, ax=ax, font_color="r")
            self.format_map_axis(ax, self.shape)

        self.canvas.draw()
