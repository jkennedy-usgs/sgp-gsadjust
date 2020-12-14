from gsadjust.plots.loop import PlotLoopAnimation

def test_PlotLoopAnimation(obstreesurvey, mainprog):
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