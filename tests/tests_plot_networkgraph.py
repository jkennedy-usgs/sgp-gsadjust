import numpy as np

from gsadjust.plots.network import PlotNetworkGraph


def test_PlotNetworkGraph(mainprog, obstreemodel):
    test_plot = PlotNetworkGraph(
        mainprog.obsTreeModel.invisibleRootItem().child(0),
        mainprog.obsTreeModel.station_coords,
        shape="circular",
        parent=mainprog,
    )
    assert test_plot.figure is not None
    assert len(test_plot.get_data()) == 4
    assert len(test_plot.get_data()[0]) == 29
    test_plot = PlotNetworkGraph(
        mainprog.obsTreeModel.invisibleRootItem().child(0),
        mainprog.obsTreeModel.station_coords,
        shape="map",
        parent=mainprog,
    )
    assert test_plot.figure is not None
    assert np.abs(test_plot.figure.axes[0].get_xlim()[0] + 106.678) < 0.01
    assert len(test_plot.get_data()) == 4
    assert len(test_plot.get_data()[0]) == 29
