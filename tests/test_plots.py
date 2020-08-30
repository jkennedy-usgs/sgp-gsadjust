import pytest
import pytestqt
import gsadjust
import numpy as np
from PyQt5 import QtCore
from test_fixture_pyqt import mainprog, obstreesurvey_with_results, obstreesurvey
from gsadjust.plots import (
    PlotGravityChange,
    PlotDatumCompare,
    PlotDgResidualHistogram,
    PlotNetworkGraph,
    PlotLoopAnimation,
)


def test_PlotLoopAnimation(mainprog):
    survey = mainprog.obsTreeModel.invisibleRootItem().child(0)
    loop = survey.child(0)
    coords = mainprog.obsTreeModel.station_coords
    lat, lon, dates = [], [], []
    for i in range(loop.rowCount()):
        station = loop.child(i)
        lat.append(coords[station.station_name][1])
        lon.append(coords[station.station_name][0])
        dates.append(station.tmean())
    test_plot = PlotLoopAnimation([lat, lon, dates])
    assert test_plot.figure is not None
    assert len(test_plot.figure.n) == 86
    assert test_plot.figure.ax1 is not None


def test_PlotDatumCompare(obstreesurvey_with_results):
    test_plot = PlotDatumCompare(obstreesurvey_with_results)
    assert test_plot.figure is not None
    assert test_plot.survey.name == '2020-05-11'
    assert test_plot.get_data()[0][0] < 0.001
    assert test_plot.get_data()[1][0] == 'rg37'


def test_PlotDgResidualHistogram(obstreesurvey_with_results):
    test_plot = PlotDgResidualHistogram(obstreesurvey_with_results)
    jeff = 1
    assert test_plot.figure is not None
    assert test_plot.get_data()[0] - 4.8 < 0.001
    assert len(test_plot.get_data()) == 45


def test_PlotNetworkGraph(mainprog):
    test_plot = PlotNetworkGraph(
        mainprog.obsTreeModel.invisibleRootItem().child(0),
        mainprog.obsTreeModel.station_coords,
        shape='circular',
        parent=None,
    )
    assert test_plot.figure is not None
    assert len(test_plot.get_data()) == 4
    assert len(test_plot.get_data()[0]) == 31
    test_plot = PlotNetworkGraph(
        mainprog.obsTreeModel.invisibleRootItem().child(0),
        mainprog.obsTreeModel.station_coords,
        shape='map',
        parent=None,
    )
    assert test_plot.figure is not None
    assert np.abs(test_plot.figure.axes[0].get_xlim()[0] + 115.05) < 0.01
    assert len(test_plot.get_data()) == 4
    assert len(test_plot.get_data()[0]) == 31
