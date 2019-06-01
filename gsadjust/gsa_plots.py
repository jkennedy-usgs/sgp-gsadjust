import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt

import networkx as nx
from PyQt5 import QtCore
from gui_objects import show_message, FigureDatumComparisonTimeSeries


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


def plot_datum_comparison_timeseries(self):
    """
    Creating a PyQt window heree in preparation for eventually having a tabbed plot window with various diagnostics
    """
    _ = FigureDatumComparisonTimeSeries(self.obsTreeModel)
